"""
Phase 8 supplementary experiment 5: empty-interior false-positive risk.
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
    DEFAULT_LARGE_M_VALUES,
    DEFAULT_RADIUS_MODES,
    DEFAULT_SUPPLEMENTARY_RESULTS_DIR,
    build_case,
    compute_case_metrics,
    ensure_supplementary_results_dirs,
    heatmap,
    lightweight_record,
    save_csv,
    save_json,
    save_markdown,
)


def _plot_interior_fraction(records: list[dict[str, Any]], save_path: str | Path) -> None:
    plt.figure(figsize=(9.2, 5.0))
    labels = sorted({f"M={row['M']} / {row['radius_mode']}" for row in records})
    for label in labels:
        rows = [row for row in records if f"M={row['M']} / {row['radius_mode']}" == label]
        plt.plot([row["k"] for row in rows], [row["interior_fraction"] for row in rows], marker="o", label=label)
    plt.xlabel("k")
    plt.ylabel("interior fraction")
    plt.title("Phase 8 supplementary: interior fraction vs k")
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=7, ncols=2)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_interior_error_flags(records: list[dict[str, Any]], save_path: str | Path) -> None:
    plt.figure(figsize=(9.2, 5.0))
    labels = sorted({f"M={row['M']} / {row['radius_mode']}" for row in records})
    for label in labels:
        rows = [row for row in records if f"M={row['M']} / {row['radius_mode']}" == label]
        plt.plot([row["k"] for row in rows], [row["interior_matching_l2"] for row in rows], marker="o", label=label)
        empty_rows = [row for row in rows if row["is_empty_interior"]]
        if empty_rows:
            plt.scatter(
                [row["k"] for row in empty_rows],
                [row["interior_matching_l2"] for row in empty_rows],
                marker="x",
                s=80,
                linewidths=2.0,
                color="black",
            )
    plt.xlabel("k")
    plt.ylabel("interior matching L2")
    plt.title("Phase 8 supplementary: interior error with empty-set flags")
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=7, ncols=2)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def run_empty_interior_risk(
    *,
    M_values: list[int] | None = None,
    k_values: list[int] | None = None,
    radius_modes: list[str] | None = None,
    h: float = 1.0,
    s: int = 4,
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    results_dirs = ensure_supplementary_results_dirs(DEFAULT_SUPPLEMENTARY_RESULTS_DIR if save_dir is None else save_dir)
    base_dir = results_dirs["exp5_empty_interior_risk"]
    M_values = list([10] + DEFAULT_LARGE_M_VALUES if M_values is None else M_values)
    k_values = list(DEFAULT_K_VALUES if k_values is None else k_values)
    radius_modes = list(DEFAULT_RADIUS_MODES if radius_modes is None else radius_modes)

    rng = np.random.default_rng(seed)
    records: list[dict[str, Any]] = []
    for M in M_values:
        q, T = build_case(M, regime="uniform", rng=rng, h=h)
        valid_k_values = [k for k in k_values if int(k) <= M]
        for mode in radius_modes:
            rows = compute_case_metrics(q, T, valid_k_values, radius_mode=mode, include_reference=False, s=s)
            for row in rows:
                row["experiment_name"] = "exp5_empty_interior_risk"
                row["regime"] = f"uniform_h={h:.2f}"
                records.append(row)

    _plot_interior_fraction(records, base_dir / "phase8_interior_fraction_vs_k.png")
    _plot_interior_error_flags(records, base_dir / "phase8_interior_error_with_empty_flags.png")
    heatmap_rows = [f"M={M}/{mode}" for M in M_values for mode in radius_modes]
    heatmap_matrix = np.asarray(
        [
            [
                (
                    1.0
                    if any(
                        row["M"] == M and row["radius_mode"] == mode and row["k"] == k and row["is_empty_interior"]
                        for row in records
                    )
                    else 0.0
                )
                for k in k_values
            ]
            for M in M_values
            for mode in radius_modes
        ],
        dtype=float,
    )
    heatmap(
        heatmap_matrix,
        x_labels=[str(k) for k in k_values],
        y_labels=heatmap_rows,
        title="Phase 8 supplementary: empty-interior risk heatmap",
        colorbar_label="1 = empty interior",
        save_path=base_dir / "phase8_empty_interior_risk.png",
        cmap="viridis",
    )

    save_csv(
        base_dir / "phase8_empty_interior_risk.csv",
        [
            "experiment_name",
            "regime",
            "M",
            "k",
            "radius_mode",
            "interior_count",
            "interior_fraction",
            "is_empty_interior",
            "interior_matching_l2",
            "full_matching_l2",
        ],
        [
            [
                row["experiment_name"],
                row["regime"],
                row["M"],
                row["k"],
                row["radius_mode"],
                row["interior_count"],
                row["interior_fraction"],
                row["is_empty_interior"],
                row["interior_matching_l2"],
                row["full_matching_l2"],
            ]
            for row in records
        ],
    )
    payload = {
        "records": [lightweight_record(row) for row in records],
        "config": {"M_values": M_values, "k_values": k_values, "radius_modes": radius_modes, "h": float(h), "seed": int(seed)},
    }
    save_json(base_dir / "summary_empty_interior_risk.json", payload)
    empty_count = sum(1 for row in records if row["is_empty_interior"])
    save_markdown(
        base_dir / "phase8_empty_interior_summary.md",
        [
            "# Phase 8 Supplementary Experiment 5",
            "",
            f"- Empty-interior cases detected: `{empty_count}/{len(records)}`.",
            "",
            "Interpretation:",
            "- Any interior-matching plot that includes empty-interior points must be annotated explicitly and must not be used as proof of a genuine interior advantage.",
            "- If strong interior-only gains survive after filtering out empty-interior cases, the Phase 8 observation is substantially more credible.",
        ],
    )
    return payload


def main() -> None:
    run_empty_interior_risk()


if __name__ == "__main__":
    main()
