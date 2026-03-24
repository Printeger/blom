"""
Phase 8 supplementary experiment 3: boundary-radius sensitivity.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from phase_8.phase8_supplementary_common import (
    DEFAULT_K_VALUES,
    DEFAULT_SUPPLEMENTARY_RESULTS_DIR,
    DEFAULT_RADIUS_MODES,
    build_case,
    compute_case_metrics,
    compute_fit_row,
    ensure_supplementary_results_dirs,
    heatmap,
    line_plot,
    radius_mode_label,
    save_csv,
    save_json,
    save_markdown,
    lightweight_record,
)


def run_radius_sensitivity(
    *,
    M: int = 40,
    k_values: list[int] | None = None,
    radius_modes: list[str] | None = None,
    h: float = 1.0,
    s: int = 4,
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    results_dirs = ensure_supplementary_results_dirs(DEFAULT_SUPPLEMENTARY_RESULTS_DIR if save_dir is None else save_dir)
    base_dir = results_dirs["exp3_radius_sensitivity"]
    k_values = list(DEFAULT_K_VALUES if k_values is None else k_values)
    radius_modes = list(DEFAULT_RADIUS_MODES if radius_modes is None else radius_modes)

    rng = np.random.default_rng(seed)
    q, T = build_case(M, regime="uniform", rng=rng, h=h)
    records: list[dict[str, Any]] = []
    fit_rows: list[dict[str, Any]] = []
    for mode in radius_modes:
        rows = compute_case_metrics(q, T, k_values, radius_mode=mode, include_reference=False, s=s)
        fit = compute_fit_row(rows)
        fit_rows.append({"radius_mode": mode, **fit})
        for row in rows:
            row["experiment_name"] = "exp3_radius_sensitivity"
            row["regime"] = f"uniform_h={h:.2f}"
            records.append(row)

    line_plot(
        {
            radius_mode_label(mode): (
                [int(row["k"]) for row in records if row["radius_mode"] == mode],
                [float(row["interior_matching_l2"]) for row in records if row["radius_mode"] == mode],
            )
            for mode in radius_modes
        },
        xlabel="k",
        ylabel="interior matching L2",
        title="Phase 8 supplementary: radius-mode interior matching",
        save_path=base_dir / "phase8_radius_mode_interior_matching.png",
    )
    line_plot(
        {
            radius_mode_label(mode): (
                [int(row["k"]) for row in records if row["radius_mode"] == mode],
                [float(row["boundary_energy_ratio"]) for row in records if row["radius_mode"] == mode],
            )
            for mode in radius_modes
        },
        xlabel="k",
        ylabel="boundary energy ratio",
        title="Phase 8 supplementary: radius-mode boundary ratio",
        save_path=base_dir / "phase8_radius_mode_boundary_ratio.png",
        ylim=(0.0, 1.05),
    )
    line_plot(
        {
            radius_mode_label(mode): (
                [int(row["k"]) for row in records if row["radius_mode"] == mode],
                [int(row["interior_count"]) for row in records if row["radius_mode"] == mode],
            )
            for mode in radius_modes
        },
        xlabel="k",
        ylabel="interior count",
        title="Phase 8 supplementary: radius-mode interior count",
        save_path=base_dir / "phase8_radius_mode_interior_count.png",
    )

    heatmap_matrix = np.asarray(
        [[float(row["interior_matching_l2"]) for row in records if row["radius_mode"] == mode] for mode in radius_modes],
        dtype=float,
    )
    heatmap(
        heatmap_matrix,
        x_labels=[str(k) for k in k_values],
        y_labels=[radius_mode_label(mode) for mode in radius_modes],
        title="Phase 8 supplementary: radius sensitivity heatmap",
        colorbar_label="interior matching L2",
        save_path=base_dir / "phase8_radius_sensitivity_heatmap.png",
    )

    save_csv(
        base_dir / "phase8_radius_sensitivity.csv",
        [
            "experiment_name",
            "regime",
            "M",
            "k",
            "radius_mode",
            "r_k",
            "interior_count",
            "boundary_count",
            "interior_matching_l2",
            "boundary_matching_l2",
            "boundary_energy_ratio",
        ],
        [
            [
                row["experiment_name"],
                row["regime"],
                row["M"],
                row["k"],
                row["radius_mode"],
                row["r_k"],
                row["interior_count"],
                row["boundary_count"],
                row["interior_matching_l2"],
                row["boundary_matching_l2"],
                row["boundary_energy_ratio"],
            ]
            for row in records
        ],
    )
    save_csv(
        base_dir / "phase8_radius_sensitivity_fit_rows.csv",
        ["radius_mode", "full_slope", "full_r2", "interior_slope", "interior_r2", "boundary_slope", "boundary_r2"],
        [
            [
                row["radius_mode"],
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
        "records": [lightweight_record(row) for row in records],
        "fit_rows": fit_rows,
        "config": {"M": int(M), "k_values": k_values, "radius_modes": radius_modes, "h": float(h), "seed": int(seed)},
    }
    save_json(base_dir / "summary_radius_sensitivity.json", payload)
    best_mode = min(fit_rows, key=lambda row: row["interior_slope"])["radius_mode"] if fit_rows else "n/a"
    save_markdown(
        base_dir / "phase8_radius_sensitivity_summary.md",
        [
            "# Phase 8 Supplementary Experiment 3",
            "",
            f"- Tested radius modes: `{radius_modes}`.",
            f"- Most negative interior slope observed under: `{best_mode}`.",
            "",
            "Interpretation:",
            "- If all radius modes preserve a clear interior/full separation, the interior-first observation is structurally robust rather than a special choice of r(k).",
            "- Large variation across radius modes would mean the future theorem should define the interior set more conservatively.",
        ],
    )
    return payload


def main() -> None:
    run_radius_sensitivity()


if __name__ == "__main__":
    main()
