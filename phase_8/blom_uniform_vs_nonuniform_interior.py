"""
Phase 8 Script 3: uniform-time vs bounded-nonuniform interior matching.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from phase_8.phase8_common import (
    DEFAULT_RESULTS_DIR,
    DEFAULT_S,
    default_q_sampler,
    ensure_results_dirs,
    interior_fit_summary,
    make_k_grid,
    prepare_case,
    sample_bounded_nonuniform_time as _sample_bounded_nonuniform_time,
    sample_uniform_time as _sample_uniform_time,
    save_csv,
    save_json,
    summarize_interior_matching,
)
from phase_7.blom_convergence_vs_k import compute_actual_blom_k, compute_ideal_truncated_blom_k, compute_minco_reference


def sample_uniform_time(M: int, h: float) -> np.ndarray:
    return _sample_uniform_time(M, h)


def sample_bounded_nonuniform_time(M: int, T_min: float, T_max: float, rng: np.random.Generator) -> np.ndarray:
    return _sample_bounded_nonuniform_time(M, T_min, T_max, rng)


def _regime_label(kind: str, param: Any) -> str:
    if kind == "uniform":
        return f"uniform_h={float(param):.2f}"
    low, high = param
    return f"bounded_[{float(low):.2f},{float(high):.2f}]"


def _aggregate_by_regime(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, int], list[dict[str, Any]]] = {}
    for record in records:
        grouped.setdefault((record["regime"], record["k"]), []).append(record)
    aggregated: list[dict[str, Any]] = []
    for (regime, k), items in sorted(grouped.items(), key=lambda item: (item[0][0], item[0][1])):
        aggregated.append(
            {
                "regime": regime,
                "k": k,
                "full_matching_l2": float(np.mean([item["full_l2"] for item in items])),
                "interior_matching_l2": float(np.mean([item["interior_l2"] for item in items])),
                "boundary_matching_l2": float(np.mean([item["boundary_l2"] for item in items])),
                "boundary_energy_ratio": float(np.mean([item["boundary_energy_ratio"] for item in items])),
            }
        )
    return aggregated


def _plot_matching_compare(aggregated: list[dict[str, Any]], save_path: str | Path) -> None:
    plt.figure(figsize=(9.2, 5.0))
    regimes = sorted({row["regime"] for row in aggregated})
    for regime in regimes:
        rows = [row for row in aggregated if row["regime"] == regime]
        plt.plot([row["k"] for row in rows], [row["interior_matching_l2"] for row in rows], marker="o", label=regime)
    plt.xlabel("k")
    plt.ylabel("mean interior matching error")
    plt.title("Phase 8 uniform vs bounded-nonuniform interior matching")
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=8, ncols=2)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_slope_boxplot(fit_rows: list[dict[str, Any]], save_path: str | Path) -> None:
    regimes = sorted({row["regime"] for row in fit_rows})
    values = [[row["interior_slope"] for row in fit_rows if row["regime"] == regime] for regime in regimes]
    plt.figure(figsize=(9.2, 4.8))
    plt.boxplot(values, tick_labels=regimes)
    plt.ylabel("interior log-fit slope")
    plt.xticks(rotation=20, ha="right")
    plt.title("Phase 8 interior slopes by time regime")
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_boundary_ratio(aggregated: list[dict[str, Any]], save_path: str | Path) -> None:
    plt.figure(figsize=(9.2, 5.0))
    regimes = sorted({row["regime"] for row in aggregated})
    for regime in regimes:
        rows = [row for row in aggregated if row["regime"] == regime]
        plt.plot([row["k"] for row in rows], [row["boundary_energy_ratio"] for row in rows], marker="o", label=regime)
    plt.xlabel("k")
    plt.ylabel("mean boundary energy ratio")
    plt.ylim(0.0, 1.05)
    plt.title("Phase 8 boundary ratios by time regime")
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=8, ncols=2)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_regime_heatmaps(records: list[dict[str, Any]], save_path: str | Path) -> None:
    regimes = []
    for record in records:
        if record["trial"] == 0 and record["regime"] not in regimes:
            regimes.append(record["regime"])
    chosen = regimes[:4]
    fig, axes = plt.subplots(len(chosen), 1, figsize=(9.0, 2.8 * len(chosen)), sharex=True)
    if len(chosen) == 1:
        axes = [axes]
    for axis, regime in zip(axes, chosen):
        regime_rows = [row for row in records if row["regime"] == regime and row["trial"] == 0]
        matrix = np.vstack([np.asarray(row["segment_errors"], dtype=float) for row in regime_rows])
        im = axis.imshow(matrix, aspect="auto", cmap="magma")
        axis.set_ylabel(regime)
        axis.set_yticks(np.arange(len(regime_rows)), [str(row["k"]) for row in regime_rows])
        fig.colorbar(im, ax=axis, fraction=0.02, pad=0.02)
    axes[-1].set_xlabel("segment index")
    fig.suptitle("Phase 8 representative regime heatmaps")
    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    plt.close(fig)


def _write_summary(aggregated: list[dict[str, Any]], fit_rows: list[dict[str, Any]], save_path: str | Path) -> None:
    uniform_slopes = [row["interior_slope"] for row in fit_rows if row["regime"].startswith("uniform")]
    bounded_slopes = [row["interior_slope"] for row in fit_rows if row["regime"].startswith("bounded")]
    uniform_mean = float(np.mean(uniform_slopes)) if uniform_slopes else 0.0
    bounded_mean = float(np.mean(bounded_slopes)) if bounded_slopes else 0.0
    lines = [
        "# Phase 8 Uniform vs Nonuniform Summary",
        "",
        f"- Mean uniform-time interior slope: `{uniform_mean:.6f}`.",
        f"- Mean bounded-nonuniform interior slope: `{bounded_mean:.6f}`.",
        f"- Uniform-time is more negative than bounded-nonuniform: `{uniform_mean < bounded_mean}`.",
        "",
        "Interpretation:",
        "- More negative uniform-time slopes indicate that the cleanest first theorem should still start from the uniform regime.",
        "- If bounded-nonuniform slopes remain negative with decent R^2, the theorem bridge is still numerically supported beyond the uniform case.",
        "- Boundary-energy statistics provide the practical limit on how far the theorem can be pushed without an explicit boundary remainder term.",
    ]
    Path(save_path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_uniform_vs_nonuniform_interior(
    q_sampler: Callable[[np.random.Generator, int], np.ndarray],
    M: int,
    k_values: list[int],
    h_values: list[float],
    nonuniform_boxes: list[tuple[float, float]],
    n_trials: int = 50,
    s: int = DEFAULT_S,
    radius_mode: str = "default",
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    results_dirs = ensure_results_dirs(DEFAULT_RESULTS_DIR if save_dir is None else save_dir)
    base_dir = results_dirs["uniform_vs_nonuniform"]

    records: list[dict[str, Any]] = []
    fit_rows: list[dict[str, Any]] = []
    regime_specs = [("uniform", h) for h in h_values] + [("bounded", box) for box in nonuniform_boxes]
    for regime_kind, param in regime_specs:
        regime = _regime_label(regime_kind, param)
        for trial in range(n_trials):
            q = q_sampler(rng, M)
            if regime_kind == "uniform":
                T = sample_uniform_time(M, float(param))
            else:
                low, high = param
                T = sample_bounded_nonuniform_time(M, float(low), float(high), rng)
            reference = compute_minco_reference(q, T, s=s)
            trial_records: list[dict[str, Any]] = []
            for k in k_values:
                actual = compute_actual_blom_k(q, T, k, s=s, scheme="C")
                ideal = compute_ideal_truncated_blom_k(q, T, k, s=s, reference=reference)
                summary = summarize_interior_matching(actual["c_vec"], ideal["c_vec"], M, k, s=s, radius_mode=radius_mode)
                record = {**summary, "regime": regime, "trial": int(trial)}
                records.append(record)
                trial_records.append(record)
            fits = interior_fit_summary(trial_records)
            fit_rows.append(
                {
                    "regime": regime,
                    "trial": int(trial),
                    "interior_slope": float(fits["interior"]["slope"]),
                    "interior_r2": float(fits["interior"]["r2"]),
                }
            )

    aggregated = _aggregate_by_regime(records)
    _plot_matching_compare(aggregated, base_dir / "uniform_vs_nonuniform_interior_matching.png")
    _plot_slope_boxplot(fit_rows, base_dir / "uniform_vs_nonuniform_slope_boxplot.png")
    _plot_boundary_ratio(aggregated, base_dir / "uniform_vs_nonuniform_boundary_ratio.png")
    _plot_regime_heatmaps(records, base_dir / "time_regime_matching_heatmaps.png")

    save_csv(
        base_dir / "uniform_vs_nonuniform_summary.csv",
        [
            "regime",
            "k",
            "full_matching_l2",
            "interior_matching_l2",
            "boundary_matching_l2",
            "boundary_energy_ratio",
        ],
        [
            [
                row["regime"],
                row["k"],
                row["full_matching_l2"],
                row["interior_matching_l2"],
                row["boundary_matching_l2"],
                row["boundary_energy_ratio"],
            ]
            for row in aggregated
        ],
    )
    save_csv(
        base_dir / "uniform_vs_nonuniform_fit_summary.csv",
        ["regime", "trial", "interior_slope", "interior_r2"],
        [[row["regime"], row["trial"], row["interior_slope"], row["interior_r2"]] for row in fit_rows],
    )
    payload = {
        "aggregated": aggregated,
        "fit_rows": fit_rows,
        "seed": int(seed),
        "n_trials": int(n_trials),
        "h_values": list(h_values),
        "nonuniform_boxes": [list(box) for box in nonuniform_boxes],
    }
    save_json(base_dir / "summary_uniform_vs_nonuniform.json", payload)
    _write_summary(aggregated, fit_rows, base_dir / "phase8_uniform_vs_nonuniform_summary.md")
    return payload


def main() -> None:
    run_uniform_vs_nonuniform_interior(
        q_sampler=default_q_sampler,
        M=10,
        k_values=[2, 4, 6, 8, 10],
        h_values=[0.5, 1.0, 2.0],
        nonuniform_boxes=[(0.9, 1.1), (0.5, 2.0), (0.2, 3.0)],
        n_trials=6,
        save_dir=DEFAULT_RESULTS_DIR,
    )


if __name__ == "__main__":
    main()

