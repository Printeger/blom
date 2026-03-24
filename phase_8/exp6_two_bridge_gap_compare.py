"""
Phase 8 supplementary experiment 6: two-bridge gap comparison.
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
    build_case,
    compute_case_metrics,
    ensure_supplementary_results_dirs,
    lightweight_record,
    save_csv,
    save_json,
    save_markdown,
)


def _plot_two_bridge_gaps(records: list[dict[str, Any]], save_path: str | Path) -> None:
    plt.figure(figsize=(8.8, 5.0))
    regimes = sorted({row["regime"] for row in records})
    for regime in regimes:
        rows = [row for row in records if row["regime"] == regime]
        grouped_raw: dict[int, list[float]] = {}
        grouped_ref: dict[int, list[float]] = {}
        for row in rows:
            grouped_raw.setdefault(int(row["k"]), []).append(float(row["raw_to_ref_l2"]))
            grouped_ref.setdefault(int(row["k"]), []).append(float(row["ref_to_ideal_l2"]))
        xs = sorted(grouped_raw)
        plt.plot(xs, [float(np.mean(grouped_raw[x])) for x in xs], marker="o", label=f"{regime} raw->ref")
        plt.plot(xs, [float(np.mean(grouped_ref[x])) for x in xs], marker="s", linestyle="--", label=f"{regime} ref->ideal")
    plt.xlabel("k")
    plt.ylabel("L2 gap")
    plt.yscale("log")
    plt.title("Phase 8 supplementary: two bridge gaps")
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=7, ncols=2)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_gap_ratio(records: list[dict[str, Any]], save_path: str | Path) -> None:
    plt.figure(figsize=(8.8, 5.0))
    regimes = sorted({row["regime"] for row in records})
    for regime in regimes:
        rows = [row for row in records if row["regime"] == regime]
        grouped: dict[int, list[float]] = {}
        for row in rows:
            grouped.setdefault(int(row["k"]), []).append(float(row["gap_ratio"]))
        xs = sorted(grouped)
        plt.plot(xs, [float(np.mean(grouped[x])) for x in xs], marker="o", label=regime)
    plt.xlabel("k")
    plt.ylabel("gap ratio")
    plt.yscale("log")
    plt.title("Phase 8 supplementary: gap ratio vs k")
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_gap_ratio_boxplot(records: list[dict[str, Any]], save_path: str | Path) -> None:
    regimes = sorted({row["regime"] for row in records})
    values = [[float(row["gap_ratio"]) for row in records if row["regime"] == regime] for regime in regimes]
    plt.figure(figsize=(8.8, 4.8))
    plt.boxplot(values, tick_labels=regimes)
    plt.yscale("log")
    plt.xticks(rotation=20, ha="right")
    plt.ylabel("gap ratio")
    plt.title("Phase 8 supplementary: gap ratio by regime")
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def run_two_bridge_gap_compare(
    *,
    M: int = 20,
    k_values: list[int] | None = None,
    h_values: list[float] | None = None,
    nonuniform_boxes: list[tuple[float, float]] | None = None,
    n_trials: int = 20,
    radius_mode: str = "k",
    s: int = 4,
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    results_dirs = ensure_supplementary_results_dirs(DEFAULT_SUPPLEMENTARY_RESULTS_DIR if save_dir is None else save_dir)
    base_dir = results_dirs["exp6_two_bridge_gaps"]
    k_values = list(DEFAULT_K_VALUES if k_values is None else k_values)
    h_values = [1.0] if h_values is None else list(h_values)
    nonuniform_boxes = [(0.5, 2.0)] if nonuniform_boxes is None else list(nonuniform_boxes)

    rng = np.random.default_rng(seed)
    records: list[dict[str, Any]] = []
    regimes = [("uniform", h) for h in h_values] + [("bounded_nonuniform", box) for box in nonuniform_boxes]
    for regime_kind, param in regimes:
        regime_name = f"uniform_h={float(param):.2f}" if regime_kind == "uniform" else f"bounded_[{float(param[0]):.2f},{float(param[1]):.2f}]"
        for trial in range(n_trials):
            q, T = build_case(M, regime=regime_kind, rng=rng, h=float(param) if regime_kind == "uniform" else 1.0, nonuniform_box=param if regime_kind != "uniform" else (0.5, 2.0))
            rows = compute_case_metrics(q, T, k_values, radius_mode=radius_mode, include_reference=True, s=s)
            for row in rows:
                row["experiment_name"] = "exp6_two_bridge_gap_compare"
                row["regime"] = regime_name
                row["trial"] = int(trial)
                records.append(row)

    _plot_two_bridge_gaps(records, base_dir / "phase8_two_bridge_gaps.png")
    _plot_gap_ratio(records, base_dir / "phase8_gap_ratio_vs_k.png")
    _plot_gap_ratio_boxplot(records, base_dir / "phase8_gap_ratio_boxplot.png")

    save_csv(
        base_dir / "phase8_two_bridge_gap_stats.csv",
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
    payload = {
        "records": [lightweight_record(row) for row in records],
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
    save_json(base_dir / "summary_two_bridge_gaps.json", payload)
    gap_ratio_mean = float(np.mean([row["gap_ratio"] for row in records])) if records else 0.0
    save_markdown(
        base_dir / "phase8_two_bridge_gaps_summary.md",
        [
            "# Phase 8 Supplementary Experiment 6",
            "",
            f"- Mean bridge gap ratio `ref->ideal / raw->ref`: `{gap_ratio_mean:.6e}`.",
            "",
            "Interpretation:",
            "- A consistently large gap ratio means the dominant unresolved bridge is reference-window versus ideal truncation, not raw Scheme C versus its own inherited-boundary reference window.",
            "- If the ratio drops toward one in some regime, that regime deserves extra caution before promoting the observation to a strong paper statement.",
        ],
    )
    return payload


def main() -> None:
    run_two_bridge_gap_compare()


if __name__ == "__main__":
    main()
