"""
Phase 8 Script 2: boundary-gap decomposition.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from phase_8.phase8_common import (
    DEFAULT_RESULTS_DIR,
    DEFAULT_S,
    compute_matching_error_blocks,
    compute_phase8_matching_triplet,
    ensure_results_dirs,
    make_interior_sets,
    prepare_case,
    project_interior_boundary_errors as _project_interior_boundary_errors,
    save_csv,
    save_json,
    validate_problem_inputs,
)


def project_interior_boundary_errors(
    segment_errors: np.ndarray,
    interior_idx: list[int],
    boundary_idx: list[int],
) -> dict[str, float]:
    return _project_interior_boundary_errors(segment_errors, interior_idx, boundary_idx)


def _plot_boundary_energy_ratio(records: list[dict[str, Any]], save_path: str | Path) -> None:
    plt.figure(figsize=(8.2, 4.8))
    plt.plot(
        [record["k"] for record in records],
        [record["boundary_energy_ratio"] for record in records],
        marker="o",
        label="boundary energy ratio",
    )
    plt.xlabel("k")
    plt.ylabel("boundary energy ratio")
    plt.ylim(0.0, 1.05)
    plt.title("Phase 8 boundary energy ratio")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_boundary_profile(records: list[dict[str, Any]], save_path: str | Path) -> None:
    representative = [records[0], records[len(records) // 2], records[-1]] if len(records) >= 3 else records
    fig, axes = plt.subplots(len(representative), 1, figsize=(9.0, 3.0 * len(representative)), sharex=True)
    if len(representative) == 1:
        axes = [axes]
    for axis, record in zip(axes, representative):
        errors = np.asarray(record["segment_errors"], dtype=float)
        x = np.arange(1, errors.size + 1)
        axis.plot(x, errors, marker="o", label=f"k={record['k']}")
        for idx in record["boundary_idx"]:
            axis.axvspan(idx - 0.5, idx + 0.5, color="tab:red", alpha=0.08)
        axis.set_ylabel("block error")
        axis.grid(True, alpha=0.2)
        axis.legend()
    axes[-1].set_xlabel("segment index")
    fig.suptitle("Phase 8 boundary-gap profiles")
    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    plt.close(fig)


def _plot_stacked_error(records: list[dict[str, Any]], save_path: str | Path) -> None:
    k_values = [record["k"] for record in records]
    interior = np.asarray([record["interior_l2"] for record in records], dtype=float)
    boundary = np.asarray([record["boundary_l2"] for record in records], dtype=float)
    plt.figure(figsize=(8.2, 4.8))
    plt.bar(k_values, interior, label="interior")
    plt.bar(k_values, boundary, bottom=interior, label="boundary")
    plt.xlabel("k")
    plt.ylabel("stacked L2 matching error")
    plt.title("Phase 8 interior/boundary stacked error")
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_radius_sensitivity(records_by_mode: dict[str, list[dict[str, Any]]], save_path: str | Path) -> None:
    plt.figure(figsize=(8.2, 4.8))
    for mode, records in records_by_mode.items():
        plt.plot(
            [record["k"] for record in records],
            [record["boundary_energy_ratio"] for record in records],
            marker="o",
            label=f"r(k) mode = {mode}",
        )
    plt.xlabel("k")
    plt.ylabel("boundary energy ratio")
    plt.ylim(0.0, 1.05)
    plt.title("Phase 8 boundary-radius sensitivity")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _write_summary(records: list[dict[str, Any]], sensitivity: dict[str, list[dict[str, Any]]], save_path: str | Path) -> None:
    mean_boundary_ratio = float(np.mean([record["boundary_energy_ratio"] for record in records])) if records else 0.0
    ratio_drop = float(records[0]["boundary_energy_ratio"] - records[-1]["boundary_energy_ratio"]) if len(records) >= 2 else 0.0
    raw_vs_ref = float(np.mean([record["raw_vs_reference_l2"] for record in records])) if records else 0.0
    ref_vs_ideal = float(np.mean([record["reference_vs_ideal_l2"] for record in records])) if records else 0.0
    lines = [
        "# Phase 8 Boundary Gap Summary",
        "",
        f"- Mean boundary energy ratio across tested k values: `{mean_boundary_ratio:.6f}`.",
        f"- Boundary energy ratio change from smallest to largest k: `{ratio_drop:.6f}`.",
        f"- Mean raw-vs-reference window gap: `{raw_vs_ref:.6e}`.",
        f"- Mean reference-window-vs-ideal gap: `{ref_vs_ideal:.6e}`.",
        "",
        "Interpretation:",
        "- A boundary energy ratio above one half means the full matching gap is mostly carried by the boundary layer.",
        "- If the ratio decreases only slowly with k, a dedicated boundary-layer theorem is likely needed.",
        "- Comparing raw-vs-reference and reference-vs-ideal separates boundary-condition mismatch from window-model mismatch.",
        f"- Radius sensitivity modes evaluated: `{', '.join(sorted(sensitivity))}`.",
    ]
    Path(save_path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_boundary_gap_decomposition(
    q: np.ndarray,
    T: np.ndarray,
    k_values: list[int],
    s: int = DEFAULT_S,
    radius_mode: str = "default",
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    q, T = validate_problem_inputs(q, T, s=s)
    results_dirs = ensure_results_dirs(DEFAULT_RESULTS_DIR if save_dir is None else save_dir)
    base_dir = results_dirs["boundary_gap"]

    records: list[dict[str, Any]] = []
    for k in k_values:
        triplet = compute_phase8_matching_triplet(q, T, k, s=s)
        actual = triplet["actual"]
        ideal = triplet["ideal"]
        reference_window = triplet["reference_window"]
        segment_errors = compute_matching_error_blocks(actual["c_vec"], ideal["c_vec"], s=s)
        sets = make_interior_sets(T.size, k, radius_mode=radius_mode)
        projection = project_interior_boundary_errors(segment_errors, sets["interior_idx"], sets["boundary_idx"])
        raw_vs_reference = compute_matching_error_blocks(actual["c_vec"], reference_window["c_vec"], s=s)
        reference_vs_ideal = compute_matching_error_blocks(reference_window["c_vec"], ideal["c_vec"], s=s)
        records.append(
            {
                "k": int(k),
                "segment_errors": segment_errors,
                "interior_idx": sets["interior_idx"],
                "boundary_idx": sets["boundary_idx"],
                "r_k": sets["r_k"],
                **projection,
                "raw_vs_reference_l2": float(np.linalg.norm(raw_vs_reference)),
                "reference_vs_ideal_l2": float(np.linalg.norm(reference_vs_ideal)),
                "gamma_norm_mean": float(np.mean(reference_window["gamma_norms"])),
            }
        )

    sensitivity_modes = ["half", "default", "one_and_half"]
    sensitivity: dict[str, list[dict[str, Any]]] = {}
    for mode in sensitivity_modes:
        mode_records: list[dict[str, Any]] = []
        for k in k_values:
            triplet = compute_phase8_matching_triplet(q, T, k, s=s)
            segment_errors = compute_matching_error_blocks(triplet["actual"]["c_vec"], triplet["ideal"]["c_vec"], s=s)
            sets = make_interior_sets(T.size, k, radius_mode=mode)
            projection = project_interior_boundary_errors(segment_errors, sets["interior_idx"], sets["boundary_idx"])
            mode_records.append({"k": int(k), **projection, "r_k": sets["r_k"]})
        sensitivity[mode] = mode_records

    _plot_boundary_energy_ratio(records, base_dir / "boundary_energy_ratio_vs_k.png")
    _plot_boundary_profile(records, base_dir / "boundary_gap_profile.png")
    _plot_stacked_error(records, base_dir / "interior_boundary_stacked_error.png")
    _plot_radius_sensitivity(sensitivity, base_dir / "boundary_radius_sensitivity.png")

    save_csv(
        base_dir / "boundary_gap_decomposition.csv",
        [
            "k",
            "r_k",
            "interior_l2",
            "boundary_l2",
            "interior_l1",
            "boundary_l1",
            "interior_energy_ratio",
            "boundary_energy_ratio",
            "raw_vs_reference_l2",
            "reference_vs_ideal_l2",
            "gamma_norm_mean",
        ],
        [
            [
                record["k"],
                record["r_k"],
                record["interior_l2"],
                record["boundary_l2"],
                record["interior_l1"],
                record["boundary_l1"],
                record["interior_energy_ratio"],
                record["boundary_energy_ratio"],
                record["raw_vs_reference_l2"],
                record["reference_vs_ideal_l2"],
                record["gamma_norm_mean"],
            ]
            for record in records
        ],
    )

    sensitivity_rows: list[list[Any]] = []
    for mode, mode_records in sensitivity.items():
        for record in mode_records:
            sensitivity_rows.append([mode, record["k"], record["r_k"], record["boundary_energy_ratio"]])
    save_csv(
        base_dir / "boundary_radius_sensitivity.csv",
        ["radius_mode", "k", "r_k", "boundary_energy_ratio"],
        sensitivity_rows,
    )
    payload = {"records": records, "radius_sensitivity": sensitivity, "seed": int(seed), "radius_mode": radius_mode}
    save_json(base_dir / "summary_boundary_gap.json", payload)
    _write_summary(records, sensitivity, base_dir / "phase8_boundary_gap_summary.md")
    return payload


def main() -> None:
    q, T = prepare_case()
    run_boundary_gap_decomposition(q, T, [2, 4, 6, 8], save_dir=DEFAULT_RESULTS_DIR)


if __name__ == "__main__":
    main()

