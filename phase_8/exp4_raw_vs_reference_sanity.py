"""
Phase 8 supplementary experiment 4: raw-vs-reference sanity checks.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from phase_8.phase8_common import extract_reference_gamma, solve_augmented_local_window
from phase_8.phase8_supplementary_common import (
    DEFAULT_K_VALUES,
    DEFAULT_LARGE_M_VALUES,
    DEFAULT_SUPPLEMENTARY_RESULTS_DIR,
    DEFAULT_RADIUS_MODES,
    build_case,
    compute_case_metrics,
    ensure_supplementary_results_dirs,
    lightweight_record,
    representative_segment_index,
    save_csv,
    save_json,
    save_markdown,
)


def _plot_raw_reference_error(rows: list[dict[str, Any]], save_path: str | Path) -> None:
    plt.figure(figsize=(8.8, 5.0))
    regimes = sorted({row["regime"] for row in rows})
    for regime in regimes:
        items = [row for row in rows if row["regime"] == regime]
        grouped: dict[int, list[float]] = {}
        for item in items:
            grouped.setdefault(int(item["k"]), []).append(float(item["raw_to_ref_l2"]))
        xs = sorted(grouped)
        ys = [float(np.mean(grouped[x])) for x in xs]
        plt.plot(xs, ys, marker="o", label=regime)
    plt.xlabel("k")
    plt.ylabel("raw vs reference-window L2")
    plt.title("Phase 8 supplementary: raw vs reference-window error")
    plt.yscale("log")
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_raw_reference_boxplot(rows: list[dict[str, Any]], save_path: str | Path) -> None:
    regimes = sorted({row["regime"] for row in rows})
    values = [[float(row["raw_to_ref_l2"]) for row in rows if row["regime"] == regime] for regime in regimes]
    plt.figure(figsize=(8.8, 4.8))
    plt.boxplot(values, tick_labels=regimes)
    plt.yscale("log")
    plt.xticks(rotation=20, ha="right")
    plt.ylabel("raw vs reference-window L2")
    plt.title("Phase 8 supplementary: raw/reference boxplot")
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_gamma_sensitivity(rows: list[dict[str, Any]], save_path: str | Path, loglog: bool = False) -> None:
    plt.figure(figsize=(8.8, 5.0))
    grouped: dict[str, dict[float, list[float]]] = {}
    for row in rows:
        grouped.setdefault(row["regime"], {}).setdefault(float(row["delta_scale"]), []).append(float(row["delta_response_l2"]))
    for regime, scale_map in sorted(grouped.items()):
        xs = sorted(scale_map)
        ys = [float(np.mean(scale_map[x])) for x in xs]
        if loglog:
            plt.loglog(xs, ys, marker="o", label=regime)
        else:
            plt.plot(xs, ys, marker="o", label=regime)
    plt.xlabel("||delta gamma|| scale")
    plt.ylabel("reference-window response L2")
    plt.title("Phase 8 supplementary: perturbed gamma sensitivity")
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def run_raw_vs_reference_sanity(
    *,
    M_values: list[int] | None = None,
    k_values: list[int] | None = None,
    n_trials: int = 12,
    h_values: list[float] | None = None,
    nonuniform_boxes: list[tuple[float, float]] | None = None,
    delta_scales: list[float] | None = None,
    radius_mode: str = "k",
    s: int = 4,
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    results_dirs = ensure_supplementary_results_dirs(DEFAULT_SUPPLEMENTARY_RESULTS_DIR if save_dir is None else save_dir)
    base_dir = results_dirs["exp4_raw_vs_reference_sanity"]
    M_values = list(DEFAULT_LARGE_M_VALUES[:3] if M_values is None else M_values)
    k_values = list(DEFAULT_K_VALUES if k_values is None else k_values)
    h_values = [1.0] if h_values is None else list(h_values)
    nonuniform_boxes = [(0.5, 2.0)] if nonuniform_boxes is None else list(nonuniform_boxes)
    delta_scales = [1e-6, 1e-5, 1e-4, 1e-3] if delta_scales is None else list(delta_scales)

    rng = np.random.default_rng(seed)
    records: list[dict[str, Any]] = []
    gamma_rows: list[dict[str, Any]] = []
    regimes = [("uniform", h) for h in h_values] + [("bounded_nonuniform", box) for box in nonuniform_boxes]

    for regime_kind, param in regimes:
        regime_name = f"uniform_h={float(param):.2f}" if regime_kind == "uniform" else f"bounded_[{float(param[0]):.2f},{float(param[1]):.2f}]"
        for M in M_values:
            for trial in range(n_trials):
                q, T = build_case(M, regime=regime_kind, rng=rng, h=float(param) if regime_kind == "uniform" else 1.0, nonuniform_box=param if regime_kind != "uniform" else (0.5, 2.0))
                rows = compute_case_metrics(q, T, k_values, radius_mode=radius_mode, include_reference=True, s=s)
                for row in rows:
                    row["experiment_name"] = "exp4_raw_vs_reference_sanity"
                    row["regime"] = regime_name
                    row["trial"] = int(trial)
                    records.append(row)

                    seg_i = representative_segment_index(M)
                    gamma_data = extract_reference_gamma(row["reference"], T, seg_i, int(row["k"]), s=s)
                    base = solve_augmented_local_window(q, T, seg_i, int(row["k"]), gamma=gamma_data["gamma"], s=s)
                    for delta_scale in delta_scales:
                        direction = rng.normal(size=gamma_data["gamma"].shape)
                        direction /= max(float(np.linalg.norm(direction)), 1e-12)
                        delta_gamma = float(delta_scale) * direction
                        perturbed = solve_augmented_local_window(q, T, seg_i, int(row["k"]), gamma=gamma_data["gamma"] + delta_gamma, s=s)
                        response = float(np.linalg.norm(perturbed["center_c_vec"] - base["center_c_vec"]))
                        gamma_rows.append(
                            {
                                "regime": regime_name,
                                "M": int(M),
                                "trial": int(trial),
                                "k": int(row["k"]),
                                "segment_i": int(seg_i),
                                "delta_scale": float(delta_scale),
                                "delta_response_l2": response,
                                "amplification_ratio": response / max(float(np.linalg.norm(delta_gamma)), 1e-15),
                            }
                        )

    _plot_raw_reference_error(records, base_dir / "phase8_raw_vs_reference_error.png")
    _plot_raw_reference_boxplot(records, base_dir / "phase8_raw_vs_reference_boxplot.png")
    _plot_gamma_sensitivity(gamma_rows, base_dir / "phase8_reference_gamma_sensitivity.png", loglog=False)
    _plot_gamma_sensitivity(gamma_rows, base_dir / "phase8_reference_gamma_loglog.png", loglog=True)

    save_csv(
        base_dir / "phase8_raw_vs_reference_stats.csv",
        ["experiment_name", "regime", "M", "trial", "k", "raw_to_ref_l2", "ref_to_ideal_l2", "gap_ratio"],
        [
            [
                row["experiment_name"],
                row["regime"],
                row["M"],
                row["trial"],
                row["k"],
                row["raw_to_ref_l2"],
                row["ref_to_ideal_l2"],
                row["gap_ratio"],
            ]
            for row in records
        ],
    )
    save_csv(
        base_dir / "phase8_reference_gamma_sensitivity.csv",
        ["regime", "M", "trial", "k", "segment_i", "delta_scale", "delta_response_l2", "amplification_ratio"],
        [
            [
                row["regime"],
                row["M"],
                row["trial"],
                row["k"],
                row["segment_i"],
                row["delta_scale"],
                row["delta_response_l2"],
                row["amplification_ratio"],
            ]
            for row in gamma_rows
        ],
    )
    payload = {
        "records": [lightweight_record(row) for row in records],
        "gamma_rows": gamma_rows,
        "config": {
            "M_values": M_values,
            "k_values": k_values,
            "n_trials": int(n_trials),
            "h_values": h_values,
            "nonuniform_boxes": [list(box) for box in nonuniform_boxes],
            "delta_scales": delta_scales,
            "radius_mode": radius_mode,
            "seed": int(seed),
        },
    }
    save_json(base_dir / "summary_raw_vs_reference_sanity.json", payload)
    raw_ref_mean = float(np.mean([row["raw_to_ref_l2"] for row in records])) if records else 0.0
    gap_ratio_mean = float(np.mean([row["gap_ratio"] for row in records])) if records else 0.0
    amp_mean = float(np.mean([row["amplification_ratio"] for row in gamma_rows])) if gamma_rows else 0.0
    save_markdown(
        base_dir / "phase8_raw_vs_reference_summary.md",
        [
            "# Phase 8 Supplementary Experiment 4",
            "",
            f"- Mean raw-to-reference error: `{raw_ref_mean:.6e}`.",
            f"- Mean bridge gap ratio ref/raw: `{gap_ratio_mean:.6e}`.",
            f"- Mean gamma-perturbation amplification ratio: `{amp_mean:.6e}`.",
            "",
            "Interpretation:",
            "- If raw-to-reference remains near machine precision across larger M and both time regimes, the previously observed agreement is credible rather than accidental.",
            "- A nontrivial gamma sensitivity confirms that the reference-window solver still reacts to inherited boundary traces instead of being numerically locked.",
        ],
    )
    return payload


def main() -> None:
    run_raw_vs_reference_sanity()


if __name__ == "__main__":
    main()
