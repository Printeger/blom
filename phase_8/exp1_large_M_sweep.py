"""
Phase 8 supplementary experiment 1: larger-M sweep with nonempty interior.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from phase_8.phase8_supplementary_common import (
    DEFAULT_K_VALUES,
    DEFAULT_LARGE_M_VALUES,
    DEFAULT_SUPPLEMENTARY_RESULTS_DIR,
    build_case,
    ensure_supplementary_results_dirs,
    heatmap,
    lightweight_record,
    line_plot,
    save_csv,
    save_json,
    save_markdown,
    compute_case_metrics,
)


def _plot_matching_curves(records: list[dict[str, Any]], save_path: str | Path) -> None:
    full_series: dict[str, tuple[list[float], list[float]]] = {}
    interior_series: dict[str, tuple[list[float], list[float]]] = {}
    boundary_series: dict[str, tuple[list[float], list[float]]] = {}
    for M in sorted({int(record["M"]) for record in records}):
        rows = [record for record in records if int(record["M"]) == M]
        k_values = [int(row["k"]) for row in rows]
        full_series[f"M={M} full"] = (k_values, [float(row["full_matching_l2"]) for row in rows])
        interior_series[f"M={M} interior"] = (k_values, [float(row["interior_matching_l2"]) for row in rows])
        boundary_series[f"M={M} boundary"] = (k_values, [float(row["boundary_matching_l2"]) for row in rows])

    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(9.4, 10.5), sharex=True)
    for axis, title, series in zip(
        axes,
        ["Full matching", "Interior-only matching", "Boundary-layer matching"],
        [full_series, interior_series, boundary_series],
    ):
        for label, (xs, ys) in series.items():
            axis.plot(xs, ys, marker="o", label=label)
        axis.set_ylabel("L2 error")
        axis.set_title(title)
        axis.grid(True, alpha=0.25)
        axis.legend(fontsize=7, ncols=2)
    axes[-1].set_xlabel("k")
    fig.suptitle("Phase 8 supplementary: larger-M matching curves")
    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    plt.close(fig)


def _plot_interior_count(records: list[dict[str, Any]], save_path: str | Path) -> None:
    series = {}
    for M in sorted({int(record["M"]) for record in records}):
        rows = [record for record in records if int(record["M"]) == M]
        series[f"M={M}"] = (
            [int(row["k"]) for row in rows],
            [int(row["interior_count"]) for row in rows],
        )
    line_plot(
        series,
        xlabel="k",
        ylabel="interior count",
        title="Phase 8 supplementary: interior count vs k",
        save_path=save_path,
    )


def _plot_boundary_ratio(records: list[dict[str, Any]], save_path: str | Path) -> None:
    series = {}
    for M in sorted({int(record["M"]) for record in records}):
        rows = [record for record in records if int(record["M"]) == M]
        series[f"M={M}"] = (
            [int(row["k"]) for row in rows],
            [float(row["boundary_energy_ratio"]) for row in rows],
        )
    line_plot(
        series,
        xlabel="k",
        ylabel="boundary energy ratio",
        title="Phase 8 supplementary: boundary ratio under larger M",
        save_path=save_path,
        ylim=(0.0, 1.05),
    )


def _plot_heatmaps(records: list[dict[str, Any]], save_path: str | Path) -> None:
    import matplotlib.pyplot as plt

    M_values = sorted({int(record["M"]) for record in records})
    fig, axes = plt.subplots(len(M_values), 1, figsize=(10.0, 2.8 * len(M_values)), sharex=False)
    if len(M_values) == 1:
        axes = [axes]
    for axis, M in zip(axes, M_values):
        rows = [record for record in records if int(record["M"]) == M]
        matrix = np.vstack([np.asarray(row["segment_errors"], dtype=float) for row in rows])
        im = axis.imshow(matrix, aspect="auto", cmap="magma")
        axis.set_title(f"M={M}")
        axis.set_ylabel("k")
        axis.set_yticks(np.arange(len(rows)), [str(row["k"]) for row in rows])
        fig.colorbar(im, ax=axis, fraction=0.02, pad=0.02)
    axes[-1].set_xlabel("segment index")
    fig.suptitle("Phase 8 supplementary: segmentwise matching heatmaps")
    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    plt.close(fig)


def run_large_M_sweep(
    M_values: list[int] | None = None,
    k_values: list[int] | None = None,
    *,
    radius_mode: str = "k",
    h: float = 1.0,
    s: int = 4,
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    results_dirs = ensure_supplementary_results_dirs(DEFAULT_SUPPLEMENTARY_RESULTS_DIR if save_dir is None else save_dir)
    base_dir = results_dirs["exp1_large_M"]
    M_values = list(DEFAULT_LARGE_M_VALUES if M_values is None else M_values)
    k_values = list(DEFAULT_K_VALUES if k_values is None else k_values)

    rng = np.random.default_rng(seed)
    records: list[dict[str, Any]] = []
    for M in M_values:
        q, T = build_case(M, regime="uniform", rng=rng, h=h)
        rows = compute_case_metrics(q, T, k_values, radius_mode=radius_mode, include_reference=False, s=s)
        for row in rows:
            row["experiment_name"] = "exp1_large_M"
            row["regime"] = f"uniform_h={h:.2f}"
            records.append(row)

    _plot_matching_curves(records, base_dir / "phase8_largeM_matching_curves.png")
    _plot_interior_count(records, base_dir / "phase8_interior_count_vs_k.png")
    _plot_boundary_ratio(records, base_dir / "phase8_boundary_ratio_largeM.png")
    _plot_heatmaps(records, base_dir / "phase8_largeM_heatmaps.png")

    save_csv(
        base_dir / "summary_large_M.csv",
        [
            "experiment_name",
            "regime",
            "M",
            "k",
            "radius_mode",
            "interior_count",
            "boundary_count",
            "full_matching_l2",
            "interior_matching_l2",
            "boundary_matching_l2",
            "boundary_energy_ratio",
            "interior_fraction",
            "is_empty_interior",
        ],
        [
            [
                row["experiment_name"],
                row["regime"],
                row["M"],
                row["k"],
                row["radius_mode"],
                row["interior_count"],
                row["boundary_count"],
                row["full_matching_l2"],
                row["interior_matching_l2"],
                row["boundary_matching_l2"],
                row["boundary_energy_ratio"],
                row["interior_fraction"],
                row["is_empty_interior"],
            ]
            for row in records
        ],
    )
    payload = {
        "records": [lightweight_record(row) for row in records],
        "config": {
            "M_values": M_values,
            "k_values": k_values,
            "radius_mode": radius_mode,
            "h": float(h),
            "seed": int(seed),
        },
    }
    save_json(base_dir / "summary_large_M.json", payload)
    mean_boundary_ratio = float(np.mean([row["boundary_energy_ratio"] for row in records])) if records else 0.0
    nonempty_rows = [row for row in records if not row["is_empty_interior"]]
    save_markdown(
        base_dir / "phase8_large_M_summary.md",
        [
            "# Phase 8 Supplementary Experiment 1",
            "",
            f"- Tested larger M values: `{M_values}`.",
            f"- Tested k values: `{k_values}`.",
            f"- Mean boundary energy ratio: `{mean_boundary_ratio:.6f}`.",
            f"- Nonempty-interior cases: `{len(nonempty_rows)}/{len(records)}`.",
            "",
            "Interpretation:",
            "- If interior-only matching stays well below full matching when interior_count > 0, the Phase 8 interior-first observation is not an empty-set artifact.",
            "- Stable boundary-energy ratios across larger M support the claim that the residual gap is still boundary-layer dominated.",
        ],
    )
    return payload


def main() -> None:
    run_large_M_sweep()


if __name__ == "__main__":
    main()
