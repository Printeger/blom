"""
Unified supplementary experiment suite for Phase 8.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from phase_8.exp1_large_M_sweep import run_large_M_sweep
from phase_8.exp2_more_trials import run_more_trials
from phase_8.exp3_radius_sensitivity import run_radius_sensitivity
from phase_8.exp4_raw_vs_reference_sanity import run_raw_vs_reference_sanity
from phase_8.exp5_empty_interior_risk import run_empty_interior_risk
from phase_8.exp6_two_bridge_gap_compare import run_two_bridge_gap_compare
from phase_8.phase8_supplementary_common import (
    DEFAULT_K_VALUES,
    DEFAULT_LARGE_M_VALUES,
    DEFAULT_RADIUS_MODES,
    DEFAULT_SUPPLEMENTARY_RESULTS_DIR,
    ensure_supplementary_results_dirs,
    save_csv,
    save_json,
    save_markdown,
)


def _collect_overview_rows(payload: dict[str, Any], experiment_name: str) -> list[list[Any]]:
    rows: list[list[Any]] = []
    for record in payload.get("records", []):
        rows.append(
            [
                experiment_name,
                record.get("regime", "n/a"),
                record.get("M", "n/a"),
                record.get("k", "n/a"),
                record.get("radius_mode", "n/a"),
                record.get("interior_count", "n/a"),
                record.get("full_matching_l2", "n/a"),
                record.get("interior_matching_l2", "n/a"),
                record.get("boundary_matching_l2", "n/a"),
                record.get("boundary_energy_ratio", "n/a"),
                record.get("raw_to_ref_l2", "n/a"),
                record.get("ref_to_ideal_l2", "n/a"),
                record.get("gap_ratio", "n/a"),
            ]
        )
    return rows


def _write_suite_summary(payload: dict[str, Any], save_path: str | Path) -> None:
    exp1_rows = payload["exp1_large_M"]["records"]
    exp2_fits = payload["exp2_more_trials"]["fit_rows"]
    exp4_rows = payload["exp4_raw_vs_reference"]["records"]
    exp5_rows = payload["exp5_empty_interior"]["records"]
    exp6_rows = payload["exp6_two_bridge"]["records"]

    nonempty_rows = [row for row in exp1_rows if not row["is_empty_interior"]]
    exp1_supports = False
    if nonempty_rows:
        exp1_supports = float(np.mean([row["interior_matching_l2"] < row["full_matching_l2"] for row in nonempty_rows])) > 0.7
    exp2_supports = float(np.mean([row["interior_slope"] < row["full_slope"] for row in exp2_fits])) > 0.7 if exp2_fits else False
    exp4_raw_ref = float(np.mean([row["raw_to_ref_l2"] for row in exp4_rows])) if exp4_rows else 0.0
    exp5_empty_fraction = float(np.mean([row["is_empty_interior"] for row in exp5_rows])) if exp5_rows else 0.0
    exp6_gap_ratio = float(np.mean([row["gap_ratio"] for row in exp6_rows])) if exp6_rows else 0.0
    bridge_choice = "reference-window / ideal-truncation consistency" if exp6_gap_ratio > 10.0 else "boundary attenuation"
    strong_observation = exp1_supports and exp2_supports and exp4_raw_ref < 1e-8 and exp5_empty_fraction < 0.25

    save_markdown(
        save_path,
        [
            "# Phase 8 Supplementary Summary",
            "",
            f"1. Larger M and more trials remain supportive: `{exp1_supports and exp2_supports}`.",
            f"2. Interior-first still holds without collapsing into empty-interior false positives: `{exp5_empty_fraction < 0.25}`.",
            f"3. `raw vs reference-window ≈ 0` remains credible: `{exp4_raw_ref < 1e-8}`.",
            f"4. Recommended next bridge theorem: `{bridge_choice}`.",
            f"5. Strong paper-level observation already supported: `{strong_observation}`.",
            "",
            "Interpretation:",
            "- Experiment 1 checks whether the interior advantage survives when the interior set stays nonempty at larger M.",
            "- Experiment 2 turns the single-case observation into a trial-level robustness statement.",
            "- Experiment 4 is the main sanity check for the surprisingly small raw-vs-reference gap.",
            "- Experiment 5 protects against the empty-interior false-positive explanation.",
            "- Experiment 6 directly identifies which bridge remains dominant numerically.",
        ],
    )


def run_phase8_supplementary_suite(
    M_values: list[int] | None = None,
    k_values: list[int] | None = None,
    radius_modes: list[str] | None = None,
    h_values: list[float] | None = None,
    nonuniform_boxes: list[tuple[float, float]] | None = None,
    n_trials: int = 50,
    s: int = 4,
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    results_dirs = ensure_supplementary_results_dirs(DEFAULT_SUPPLEMENTARY_RESULTS_DIR if save_dir is None else save_dir)
    M_values = list(DEFAULT_LARGE_M_VALUES if M_values is None else M_values)
    k_values = list(DEFAULT_K_VALUES if k_values is None else k_values)
    radius_modes = list(DEFAULT_RADIUS_MODES if radius_modes is None else radius_modes)
    h_values = [1.0] if h_values is None else list(h_values)
    nonuniform_boxes = [(0.5, 2.0)] if nonuniform_boxes is None else list(nonuniform_boxes)

    exp1 = run_large_M_sweep(M_values=M_values, k_values=k_values, radius_mode="k", h=h_values[0], s=s, save_dir=results_dirs["base"], seed=seed)
    exp2 = run_more_trials(M=max(M_values[0], 20), k_values=k_values, h_values=h_values, nonuniform_boxes=nonuniform_boxes, n_trials=n_trials, radius_mode="k", s=s, save_dir=results_dirs["base"], seed=seed)
    exp3 = run_radius_sensitivity(M=max(M_values[1] if len(M_values) > 1 else M_values[0], 40), k_values=k_values, radius_modes=radius_modes, h=h_values[0], s=s, save_dir=results_dirs["base"], seed=seed)
    exp4 = run_raw_vs_reference_sanity(M_values=M_values[: min(len(M_values), 3)], k_values=k_values, n_trials=max(8, min(n_trials, 20)), h_values=h_values, nonuniform_boxes=nonuniform_boxes, radius_mode="k", s=s, save_dir=results_dirs["base"], seed=seed)
    exp5 = run_empty_interior_risk(M_values=[10] + [M for M in M_values if M != 10], k_values=k_values, radius_modes=radius_modes, h=h_values[0], s=s, save_dir=results_dirs["base"], seed=seed)
    exp6 = run_two_bridge_gap_compare(M=max(M_values[0], 20), k_values=k_values, h_values=h_values, nonuniform_boxes=nonuniform_boxes, n_trials=max(8, min(n_trials, 20)), radius_mode="k", s=s, save_dir=results_dirs["base"], seed=seed)

    payload = {
        "exp1_large_M": exp1,
        "exp2_more_trials": exp2,
        "exp3_radius_sensitivity": exp3,
        "exp4_raw_vs_reference": exp4,
        "exp5_empty_interior": exp5,
        "exp6_two_bridge": exp6,
        "config": {
            "M_values": M_values,
            "k_values": k_values,
            "radius_modes": radius_modes,
            "h_values": h_values,
            "nonuniform_boxes": [list(box) for box in nonuniform_boxes],
            "n_trials": int(n_trials),
            "s": int(s),
            "seed": int(seed),
        },
    }

    overview_rows: list[list[Any]] = []
    overview_rows.extend(_collect_overview_rows(exp1, "exp1_large_M"))
    overview_rows.extend(_collect_overview_rows(exp2, "exp2_more_trials"))
    overview_rows.extend(_collect_overview_rows(exp3, "exp3_radius_sensitivity"))
    overview_rows.extend(_collect_overview_rows(exp4, "exp4_raw_vs_reference_sanity"))
    overview_rows.extend(_collect_overview_rows(exp5, "exp5_empty_interior_risk"))
    overview_rows.extend(_collect_overview_rows(exp6, "exp6_two_bridge_gap_compare"))
    save_csv(
        results_dirs["summary"] / "phase8_supplementary_overview.csv",
        [
            "experiment_name",
            "regime",
            "M",
            "k",
            "radius_mode",
            "interior_count",
            "full_matching_l2",
            "interior_matching_l2",
            "boundary_matching_l2",
            "boundary_energy_ratio",
            "raw_to_ref_l2",
            "ref_to_ideal_l2",
            "gap_ratio",
        ],
        overview_rows,
    )
    save_json(results_dirs["summary"] / "summary_phase8_supplementary_suite.json", payload)
    _write_suite_summary(payload, results_dirs["summary"] / "phase8_supplementary_summary.md")
    return payload


def main() -> None:
    run_phase8_supplementary_suite()


if __name__ == "__main__":
    main()
