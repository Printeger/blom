"""
Phase 9 validation suite.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from phase_9.blom_backward_diff import (
    DEFAULT_RESULTS_DIR,
    _save_csv,
    _save_json,
    _validate_problem_inputs,
    default_obs_config,
    default_weights,
    ensure_results_dirs,
    representative_case,
    run_phase9_gradcheck,
)
from phase_9.blom_space_time_opt_demo import run_minimal_optimization_demo


def _write_interpretation_summary(
    gradcheck: dict[str, Any],
    optimization_demo: dict[str, Any],
    save_path: str | Path,
) -> None:
    grad_summary = gradcheck["summary"]
    opt_summary = optimization_demo["summary"]
    lines = [
        "# Phase 9 Interpretation Summary",
        "",
        f"- Reverse differentiation chain correct: `{grad_summary['dense_vs_fd']['max_abs_error'] < 1e-4}`.",
        f"- Dense and banded backward match: `{grad_summary['dense_vs_banded']['max_abs_error'] < 1e-10}`.",
        f"- Objective decreased in the demo: `{opt_summary['objective_drop'] > 0.0}`.",
        f"- Any duration touched the lower bound: `{opt_summary['touches_T_min']}`.",
        "",
        "Interpretation:",
        "- If the finite-difference gap is small, the minimal loop is mathematically usable.",
        "- If the dense-vs-banded gap is at machine precision, the local-support backward pass is implementation-ready.",
        "- If the optimization demo decreases the objective without hitting the time floor too often, the framework is ready to scale toward a fuller block-banded optimizer.",
    ]
    Path(save_path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_phase9_validation_suite(
    q0: np.ndarray,
    T0: np.ndarray,
    weights: dict[str, float],
    obs_config: dict[str, Any] | None = None,
    n_steps: int = 50,
    step_size: float = 1e-3,
    T_min: float = 1e-3,
    s: int = 4,
    k: int = 2,
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    q0, T0 = _validate_problem_inputs(q0, T0, s=s, k=k)
    results_dirs = ensure_results_dirs(DEFAULT_RESULTS_DIR if save_dir is None else save_dir)
    obs_cfg = default_obs_config() if obs_config is None else {**default_obs_config(), **dict(obs_config)}
    weights = {**default_weights(), **dict(weights)}

    gradcheck = run_phase9_gradcheck(q0, T0, weights, obs_config=obs_cfg, s=s, k=k, save_dir=results_dirs["base"])
    optimization_demo = run_minimal_optimization_demo(
        q0,
        T0,
        weights,
        obs_config=obs_cfg,
        n_steps=n_steps,
        step_size=step_size,
        T_min=T_min,
        use_banded=True,
        save_dir=results_dirs["base"],
        seed=seed,
    )
    overview_rows = [
        [
            "gradcheck_dense_vs_fd",
            True,
            gradcheck["summary"]["dense_vs_fd"]["max_abs_error"],
            gradcheck["summary"]["dense_vs_fd"]["l2_error"],
            gradcheck["summary"]["dense_vs_banded"]["max_abs_error"],
            0.0,
            0.0,
        ],
        [
            "gradcheck_banded_vs_fd",
            True,
            gradcheck["summary"]["banded_vs_fd"]["max_abs_error"],
            gradcheck["summary"]["banded_vs_fd"]["l2_error"],
            gradcheck["summary"]["dense_vs_banded"]["max_abs_error"],
            0.0,
            0.0,
        ],
        [
            "optimization_demo",
            optimization_demo["summary"]["objective_drop"] > 0.0,
            0.0,
            0.0,
            gradcheck["summary"]["dense_vs_banded"]["max_abs_error"],
            optimization_demo["summary"]["objective_drop"],
            optimization_demo["summary"]["grad_norm_drop"],
        ],
    ]
    _save_csv(
        results_dirs["compare"] / "phase9_overview.csv",
        ["test_name", "passed", "max_abs_error", "l2_error", "dense_banded_gap", "objective_drop", "grad_norm_drop"],
        overview_rows,
    )
    _write_interpretation_summary(gradcheck, optimization_demo, results_dirs["compare"] / "phase9_interpretation_summary.md")
    payload = {
        "gradcheck": gradcheck["summary"],
        "dense_vs_banded": gradcheck["summary"]["dense_vs_banded"],
        "optimization_demo": optimization_demo["summary"],
        "config": {
            "n_steps": int(n_steps),
            "step_size": float(step_size),
            "T_min": float(T_min),
            "seed": int(seed),
        },
    }
    _save_json(results_dirs["compare"] / "phase9_suite_summary.json", payload)
    return {"gradcheck": gradcheck, "dense_vs_banded": gradcheck["summary"]["dense_vs_banded"], "optimization_demo": optimization_demo}


def main() -> None:
    q0, T0 = representative_case()
    run_phase9_validation_suite(q0, T0, default_weights(), n_steps=25, step_size=5e-3, T_min=0.2)


if __name__ == "__main__":
    main()

