"""
Phase 8 unified validation suite.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from phase_8.blom_boundary_gap_decomposition import run_boundary_gap_decomposition
from phase_8.blom_interior_matching_check import run_interior_matching_check
from phase_8.blom_uniform_vs_nonuniform_interior import run_uniform_vs_nonuniform_interior
from phase_8.phase8_common import (
    DEFAULT_RESULTS_DIR,
    DEFAULT_S,
    default_q_sampler,
    ensure_results_dirs,
    make_k_grid,
    prepare_case,
    save_csv,
    save_json,
)


def _write_interpretation_summary(
    interior_payload: dict[str, Any],
    boundary_payload: dict[str, Any],
    regime_payload: dict[str, Any],
    save_path: str | Path,
) -> None:
    interior_fit = interior_payload["fits"]["interior"]
    boundary_ratios = [record["boundary_energy_ratio"] for record in boundary_payload["records"]]
    regime_slopes = [row["interior_slope"] for row in regime_payload["fit_rows"]]
    uniform_slopes = [row["interior_slope"] for row in regime_payload["fit_rows"] if row["regime"].startswith("uniform")]
    bounded_slopes = [row["interior_slope"] for row in regime_payload["fit_rows"] if row["regime"].startswith("bounded")]
    lines = [
        "# Phase 8 Interpretation Summary",
        "",
        f"- Interior-first matching significant: `{interior_fit['slope'] < interior_payload['fits']['full']['slope']}`.",
        f"- Mean boundary energy ratio: `{sum(boundary_ratios) / max(len(boundary_ratios), 1):.6f}`.",
        f"- Uniform mean interior slope: `{sum(uniform_slopes) / max(len(uniform_slopes), 1):.6f}`.",
        f"- Bounded mean interior slope: `{sum(bounded_slopes) / max(len(bounded_slopes), 1):.6f}`.",
        "",
        "Recommended next theorem target:",
        "- If the interior slope remains clearly negative while the boundary ratio stays large, prioritize a boundary-layer remainder theorem.",
        "- If raw-vs-reference gaps fall faster than reference-vs-ideal gaps, prioritize a reference-window / ideal-truncation consistency proposition.",
        "- If both uniform and bounded-nonuniform regimes remain stable, use the uniform-time case as the clean first theorem and then generalize.",
    ]
    Path(save_path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_phase8_validation_suite(
    q: Any = None,
    T: Any = None,
    M: int = 20,
    k_values: list[int] | None = None,
    h_values: list[float] | None = None,
    nonuniform_boxes: list[tuple[float, float]] | None = None,
    n_trials: int = 50,
    s: int = DEFAULT_S,
    radius_mode: str = "default",
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    results_dirs = ensure_results_dirs(DEFAULT_RESULTS_DIR if save_dir is None else save_dir)
    q_case, T_case = prepare_case(q, T, M=M, seed=seed)
    k_values = [2, 4, 6, 8] if k_values is None and T_case.size == 8 else (make_k_grid(T_case.size, start=2, step=2, include_M=True) if k_values is None else list(k_values))
    h_values = [0.5, 1.0, 2.0] if h_values is None else list(h_values)
    nonuniform_boxes = [(0.9, 1.1), (0.5, 2.0), (0.2, 3.0)] if nonuniform_boxes is None else list(nonuniform_boxes)

    interior_payload = run_interior_matching_check(
        q_case,
        T_case,
        k_values,
        s=s,
        radius_mode=radius_mode,
        save_dir=results_dirs["base"],
        seed=seed,
    )
    boundary_payload = run_boundary_gap_decomposition(
        q_case,
        T_case,
        k_values,
        s=s,
        radius_mode=radius_mode,
        save_dir=results_dirs["base"],
        seed=seed,
    )
    regime_payload = run_uniform_vs_nonuniform_interior(
        q_sampler=default_q_sampler,
        M=max(M, 10),
        k_values=make_k_grid(max(M, 10), start=2, step=2, include_M=False)[: min(len(k_values), 5)],
        h_values=h_values,
        nonuniform_boxes=nonuniform_boxes,
        n_trials=n_trials,
        s=s,
        radius_mode=radius_mode,
        save_dir=results_dirs["base"],
        seed=seed,
    )

    overview_rows: list[list[Any]] = []
    interior_r2 = interior_payload["fits"]["interior"]["r2"]
    interior_slope = interior_payload["fits"]["interior"]["slope"]
    boundary_map = {record["k"]: record["boundary_energy_ratio"] for record in boundary_payload["records"]}
    for record in interior_payload["records"]:
        overview_rows.append(
            [
                "single_case",
                "deterministic",
                record["k"],
                record["full_l2"],
                record["interior_l2"],
                record["boundary_l2"],
                interior_slope,
                interior_r2,
                boundary_map[record["k"]],
            ]
        )
    for row in regime_payload["aggregated"]:
        overview_rows.append(
            [
                "regime_split",
                row["regime"],
                row["k"],
                row["full_matching_l2"],
                row["interior_matching_l2"],
                row["boundary_matching_l2"],
                float("nan"),
                float("nan"),
                row["boundary_energy_ratio"],
            ]
        )

    save_csv(
        results_dirs["compare"] / "phase8_overview.csv",
        [
            "experiment_name",
            "regime",
            "k",
            "full_matching_l2",
            "interior_matching_l2",
            "boundary_matching_l2",
            "interior_slope",
            "interior_r2",
            "boundary_energy_ratio",
        ],
        overview_rows,
    )
    _write_interpretation_summary(
        interior_payload,
        boundary_payload,
        regime_payload,
        results_dirs["compare"] / "phase8_interpretation_summary.md",
    )
    payload = {
        "single_case": interior_payload,
        "boundary_gap": boundary_payload,
        "regime_split": regime_payload,
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
    save_json(results_dirs["compare"] / "summary_phase8_suite.json", payload)
    return payload


def main() -> None:
    run_phase8_validation_suite(M=10, n_trials=6, save_dir=DEFAULT_RESULTS_DIR)


if __name__ == "__main__":
    main()

