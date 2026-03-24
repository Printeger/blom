"""
Phase 8 Script 1: interior-first matching check.
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
    compute_ideal_truncated_blom_k,
    compute_matching_error_blocks as _compute_matching_error_blocks,
    compute_minco_reference,
    compute_actual_blom_k,
    ensure_results_dirs,
    interior_fit_summary,
    make_interior_sets as _make_interior_sets,
    prepare_case,
    save_csv,
    save_json,
    summarize_interior_matching as _summarize_interior_matching,
    validate_problem_inputs,
)


def make_interior_sets(M: int, k: int, radius_mode: str = "default") -> dict[str, Any]:
    return _make_interior_sets(M, k, radius_mode=radius_mode)


def compute_matching_error_blocks(c_actual: np.ndarray, c_ideal: np.ndarray, s: int = DEFAULT_S) -> np.ndarray:
    return _compute_matching_error_blocks(c_actual, c_ideal, s=s)


def summarize_interior_matching(
    c_actual: np.ndarray,
    c_ideal: np.ndarray,
    M: int,
    k: int,
    s: int = DEFAULT_S,
    radius_mode: str = "default",
) -> dict[str, Any]:
    return _summarize_interior_matching(c_actual, c_ideal, M, k, s=s, radius_mode=radius_mode)


def _plot_full_vs_interior(records: list[dict[str, Any]], save_path: str | Path) -> None:
    k_values = [record["k"] for record in records]
    plt.figure(figsize=(8.2, 4.8))
    plt.plot(k_values, [record["full_l2"] for record in records], marker="o", label="full matching")
    plt.plot(k_values, [record["interior_l2"] for record in records], marker="s", label="interior matching")
    plt.plot(k_values, [record["boundary_l2"] for record in records], marker="^", label="boundary matching")
    plt.xlabel("k")
    plt.ylabel("L2 matching error")
    plt.title("Phase 8 full vs interior matching")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_logfit(records: list[dict[str, Any]], fits: dict[str, Any], save_path: str | Path) -> None:
    k_values = np.asarray([record["k"] for record in records], dtype=float)
    interior = np.asarray([record["interior_l2"] for record in records], dtype=float)
    full = np.asarray([record["full_l2"] for record in records], dtype=float)
    boundary = np.asarray([record["boundary_l2"] for record in records], dtype=float)
    plt.figure(figsize=(8.2, 4.8))
    plt.scatter(k_values, np.log(np.maximum(full, 1e-15)), label="log full", marker="o")
    plt.scatter(k_values, np.log(np.maximum(interior, 1e-15)), label="log interior", marker="s")
    plt.scatter(k_values, np.log(np.maximum(boundary, 1e-15)), label="log boundary", marker="^")
    fit = fits["interior"]
    if fit["valid_count"] >= 2:
        plt.plot(k_values, fit["fit_log_errors"], linestyle="--", label="interior fit")
    plt.xlabel("k")
    plt.ylabel("log error")
    plt.title("Phase 8 interior matching log-fit")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_segmentwise_heatmap(records: list[dict[str, Any]], save_path: str | Path) -> None:
    matrix = np.vstack([np.asarray(record["segment_errors"], dtype=float) for record in records])
    plt.figure(figsize=(9.0, 4.8))
    im = plt.imshow(matrix, aspect="auto", cmap="magma")
    plt.colorbar(im, label="segmentwise matching error")
    plt.yticks(np.arange(len(records)), [str(record["k"]) for record in records])
    plt.xlabel("segment index")
    plt.ylabel("k")
    plt.title("Phase 8 segmentwise matching heatmap")
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_boxplot(records: list[dict[str, Any]], save_path: str | Path) -> None:
    boundary_values: list[float] = []
    interior_values: list[float] = []
    for record in records:
        errors = np.asarray(record["segment_errors"], dtype=float)
        interior_idx = np.asarray(record["interior_idx"], dtype=int)
        boundary_idx = np.asarray(record["boundary_idx"], dtype=int)
        if interior_idx.size:
            interior_values.extend(errors[interior_idx - 1].tolist())
        if boundary_idx.size:
            boundary_values.extend(errors[boundary_idx - 1].tolist())
    plt.figure(figsize=(6.6, 4.6))
    plt.boxplot([interior_values, boundary_values], tick_labels=["interior", "boundary"])
    plt.ylabel("segmentwise matching error")
    plt.title("Phase 8 boundary vs interior matching")
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _write_summary(records: list[dict[str, Any]], fits: dict[str, Any], save_path: str | Path) -> None:
    latest = records[-1]
    interior_better = sum(record["interior_l2"] < record["full_l2"] for record in records)
    boundary_dominates = sum(record["boundary_l2"] > record["interior_l2"] for record in records)
    lines = [
        "# Phase 8 Interior Matching Summary",
        "",
        f"- Interior-only matching beats full matching for `{interior_better}/{len(records)}` tested k values.",
        f"- Boundary-only matching exceeds interior-only matching for `{boundary_dominates}/{len(records)}` tested k values.",
        f"- Interior log-fit slope: `{fits['interior']['slope']:.6f}` with `R^2 = {fits['interior']['r2']:.6f}`.",
        f"- Full log-fit slope: `{fits['full']['slope']:.6f}`; boundary log-fit slope: `{fits['boundary']['slope']:.6f}`.",
        f"- Largest tested k = `{latest['k']}` gives interior/full ratio `{latest['interior_l2'] / max(latest['full_l2'], 1e-15):.6f}`.",
        "",
        "Interpretation:",
        "- Interior-only matching is numerically cleaner whenever its curve sits well below the full-domain curve.",
        "- A more negative interior slope than the full-domain slope supports an interior-first matching theorem.",
        "- Persistent boundary dominance suggests the remaining gap should be treated as a dedicated boundary-layer remainder, not hidden inside a single full-domain statement.",
    ]
    Path(save_path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_interior_matching_check(
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
    base_dir = results_dirs["interior_matching"]
    reference = compute_minco_reference(q, T, s=s)
    records: list[dict[str, Any]] = []
    for k in k_values:
        actual = compute_actual_blom_k(q, T, k, s=s, scheme="C")
        ideal = compute_ideal_truncated_blom_k(q, T, k, s=s, reference=reference)
        summary = summarize_interior_matching(actual["c_vec"], ideal["c_vec"], T.size, k, s=s, radius_mode=radius_mode)
        records.append(summary)

    fits = interior_fit_summary(records)
    _plot_full_vs_interior(records, base_dir / "interior_matching_full_vs_interior.png")
    _plot_logfit(records, fits, base_dir / "interior_matching_logfit.png")
    _plot_segmentwise_heatmap(records, base_dir / "segmentwise_matching_heatmap.png")
    _plot_boxplot(records, base_dir / "boundary_vs_interior_matching_boxplot.png")

    save_csv(
        base_dir / "interior_matching_by_k.csv",
        [
            "k",
            "r_k",
            "full_l2",
            "interior_l2",
            "boundary_l2",
            "full_linf",
            "interior_linf",
            "boundary_linf",
            "boundary_energy_ratio",
        ],
        [
            [
                record["k"],
                record["r_k"],
                record["full_l2"],
                record["interior_l2"],
                record["boundary_l2"],
                record["full_linf"],
                record["interior_linf"],
                record["boundary_linf"],
                record["boundary_energy_ratio"],
            ]
            for record in records
        ],
    )
    save_csv(
        base_dir / "interior_matching_fit_summary.csv",
        ["curve", "slope", "intercept", "r2", "k_min", "k_max", "valid_count"],
        [
            [
                name,
                fit["slope"],
                fit["intercept"],
                fit["r2"],
                fit["k_min"],
                fit["k_max"],
                fit["valid_count"],
            ]
            for name, fit in fits.items()
        ],
    )

    payload = {"records": records, "fits": fits, "seed": int(seed), "radius_mode": radius_mode}
    save_json(base_dir / "summary_interior_matching.json", payload)
    _write_summary(records, fits, base_dir / "phase8_interior_matching_summary.md")
    return payload


def main() -> None:
    q, T = prepare_case()
    run_interior_matching_check(q, T, [2, 4, 6, 8], save_dir=DEFAULT_RESULTS_DIR)


if __name__ == "__main__":
    main()
