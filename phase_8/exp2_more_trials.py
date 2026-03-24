"""
Phase 8 supplementary experiment 2: more random-trial robustness checks.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from phase_8.phase8_supplementary_common import (
    DEFAULT_K_VALUES,
    DEFAULT_SUPPLEMENTARY_RESULTS_DIR,
    aggregate_numeric_records,
    band_plot,
    boxplot_by_group,
    build_case,
    compute_case_metrics,
    compute_fit_row,
    ensure_supplementary_results_dirs,
    lightweight_record,
    save_csv,
    save_json,
    save_markdown,
)


def _plot_slope_scatter(fit_rows: list[dict[str, Any]], save_path: str | Path) -> None:
    plt.figure(figsize=(8.0, 5.0))
    for regime in sorted({row["regime"] for row in fit_rows}):
        rows = [row for row in fit_rows if row["regime"] == regime]
        plt.scatter(
            [row["full_slope"] for row in rows],
            [row["interior_slope"] for row in rows],
            label=regime,
            alpha=0.75,
        )
    plt.xlabel("full log-fit slope")
    plt.ylabel("interior log-fit slope")
    plt.title("Phase 8 supplementary: full vs interior slopes")
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def run_more_trials(
    *,
    M: int = 20,
    k_values: list[int] | None = None,
    h_values: list[float] | None = None,
    nonuniform_boxes: list[tuple[float, float]] | None = None,
    n_trials: int = 50,
    radius_mode: str = "k",
    s: int = 4,
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    results_dirs = ensure_supplementary_results_dirs(DEFAULT_SUPPLEMENTARY_RESULTS_DIR if save_dir is None else save_dir)
    base_dir = results_dirs["exp2_more_trials"]
    k_values = list(DEFAULT_K_VALUES if k_values is None else k_values)
    h_values = [1.0] if h_values is None else list(h_values)
    nonuniform_boxes = [(0.5, 2.0)] if nonuniform_boxes is None else list(nonuniform_boxes)

    rng = np.random.default_rng(seed)
    trial_records: list[dict[str, Any]] = []
    fit_rows: list[dict[str, Any]] = []
    regimes = [("uniform", h) for h in h_values] + [("bounded_nonuniform", box) for box in nonuniform_boxes]

    for regime_kind, param in regimes:
        regime_name = f"uniform_h={float(param):.2f}" if regime_kind == "uniform" else f"bounded_[{float(param[0]):.2f},{float(param[1]):.2f}]"
        for trial in range(n_trials):
            q, T = build_case(M, regime=regime_kind, rng=rng, h=float(param) if regime_kind == "uniform" else 1.0, nonuniform_box=param if regime_kind != "uniform" else (0.5, 2.0))
            records = compute_case_metrics(q, T, k_values, radius_mode=radius_mode, include_reference=False, s=s)
            fit = compute_fit_row(records)
            fit_rows.append({"regime": regime_name, "trial": int(trial), **fit})
            for row in records:
                row["experiment_name"] = "exp2_more_trials"
                row["regime"] = regime_name
                row["trial"] = int(trial)
                trial_records.append(row)

    aggregated = aggregate_numeric_records(
        trial_records,
        group_keys=["regime", "k"],
        metric_keys=["full_matching_l2", "interior_matching_l2", "boundary_matching_l2", "boundary_energy_ratio"],
    )

    slope_groups = {regime: [row["interior_slope"] for row in fit_rows if row["regime"] == regime] for regime in sorted({row["regime"] for row in fit_rows})}
    boxplot_by_group(
        slope_groups,
        ylabel="interior log-fit slope",
        title="Phase 8 supplementary: interior slope distribution",
        save_path=base_dir / "phase8_interior_slope_boxplot.png",
    )
    ratio_groups = {regime: [float(np.mean([item["boundary_energy_ratio"] for item in trial_records if item["regime"] == regime and item["trial"] == trial])) for trial in range(n_trials)] for regime in sorted({row["regime"] for row in fit_rows})}
    boxplot_by_group(
        ratio_groups,
        ylabel="mean boundary energy ratio",
        title="Phase 8 supplementary: boundary ratio distribution",
        save_path=base_dir / "phase8_boundary_ratio_boxplot.png",
    )
    _plot_slope_scatter(fit_rows, base_dir / "phase8_full_vs_interior_slope_scatter.png")

    band_series: dict[str, tuple[list[float], list[float], list[float]]] = {}
    for regime in sorted({row["regime"] for row in aggregated}):
        rows = [row for row in aggregated if row["regime"] == regime]
        band_series[regime] = (
            [int(row["k"]) for row in rows],
            [float(row["interior_matching_l2_mean"]) for row in rows],
            [float(row["interior_matching_l2_std"]) for row in rows],
        )
    band_plot(
        band_series,
        xlabel="k",
        ylabel="interior matching L2",
        title="Phase 8 supplementary: mean matching curve with std band",
        save_path=base_dir / "phase8_trial_mean_curves.png",
    )

    save_csv(
        base_dir / "phase8_more_trials_aggregated.csv",
        [
            "regime",
            "k",
            "full_matching_l2_mean",
            "full_matching_l2_median",
            "full_matching_l2_std",
            "full_matching_l2_q25",
            "full_matching_l2_q75",
            "interior_matching_l2_mean",
            "interior_matching_l2_median",
            "interior_matching_l2_std",
            "interior_matching_l2_q25",
            "interior_matching_l2_q75",
            "boundary_matching_l2_mean",
            "boundary_matching_l2_median",
            "boundary_matching_l2_std",
            "boundary_matching_l2_q25",
            "boundary_matching_l2_q75",
            "boundary_energy_ratio_mean",
            "boundary_energy_ratio_median",
            "boundary_energy_ratio_std",
            "boundary_energy_ratio_q25",
            "boundary_energy_ratio_q75",
        ],
        [
            [
                row["regime"],
                row["k"],
                row["full_matching_l2_mean"],
                row["full_matching_l2_median"],
                row["full_matching_l2_std"],
                row["full_matching_l2_q25"],
                row["full_matching_l2_q75"],
                row["interior_matching_l2_mean"],
                row["interior_matching_l2_median"],
                row["interior_matching_l2_std"],
                row["interior_matching_l2_q25"],
                row["interior_matching_l2_q75"],
                row["boundary_matching_l2_mean"],
                row["boundary_matching_l2_median"],
                row["boundary_matching_l2_std"],
                row["boundary_matching_l2_q25"],
                row["boundary_matching_l2_q75"],
                row["boundary_energy_ratio_mean"],
                row["boundary_energy_ratio_median"],
                row["boundary_energy_ratio_std"],
                row["boundary_energy_ratio_q25"],
                row["boundary_energy_ratio_q75"],
            ]
            for row in aggregated
        ],
    )
    save_csv(
        base_dir / "phase8_more_trials_fit_rows.csv",
        ["regime", "trial", "full_slope", "full_r2", "interior_slope", "interior_r2", "boundary_slope", "boundary_r2"],
        [
            [
                row["regime"],
                row["trial"],
                row["full_slope"],
                row["full_r2"],
                row["interior_slope"],
                row["interior_r2"],
                row["boundary_slope"],
                row["boundary_r2"],
            ]
            for row in fit_rows
        ],
    )
    payload = {
        "aggregated": aggregated,
        "fit_rows": fit_rows,
        "records": [lightweight_record(row) for row in trial_records],
        "config": {
            "M": int(M),
            "k_values": k_values,
            "h_values": h_values,
            "nonuniform_boxes": [list(box) for box in nonuniform_boxes],
            "n_trials": int(n_trials),
            "radius_mode": radius_mode,
            "seed": int(seed),
        },
    }
    save_json(base_dir / "summary_more_trials.json", payload)

    interior_slope_mean = float(np.mean([row["interior_slope"] for row in fit_rows])) if fit_rows else 0.0
    full_slope_mean = float(np.mean([row["full_slope"] for row in fit_rows])) if fit_rows else 0.0
    boundary_ratio_mean = float(np.mean([row["boundary_energy_ratio"] for row in trial_records])) if trial_records else 0.0
    save_markdown(
        base_dir / "phase8_more_trials_summary.md",
        [
            "# Phase 8 Supplementary Experiment 2",
            "",
            f"- Number of random trials per regime: `{n_trials}`.",
            f"- Mean interior slope across all trials: `{interior_slope_mean:.6f}`.",
            f"- Mean full slope across all trials: `{full_slope_mean:.6f}`.",
            f"- Mean boundary energy ratio across all per-k samples: `{boundary_ratio_mean:.6f}`.",
            "",
            "Interpretation:",
            "- If the interior slope distribution stays consistently below the full slope distribution, interior-first matching is a stable empirical law rather than a single-case artifact.",
            "- A large boundary-ratio distribution across regimes indicates that boundary-layer dominance persists beyond the representative Phase 8 case.",
        ],
    )
    return payload


def main() -> None:
    run_more_trials(n_trials=50)


if __name__ == "__main__":
    main()
