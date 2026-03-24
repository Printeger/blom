"""
Phase 9 minimal space-time optimization demo.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from phase_1.minco_scalar_baseline import sample_trajectory
from phase_9.blom_backward_diff import (
    DEFAULT_K,
    DEFAULT_RESULTS_DIR,
    DEFAULT_S,
    _save_csv,
    _save_json,
    _validate_problem_inputs,
    backward_diff_banded,
    backward_diff_dense,
    default_obs_config,
    default_weights,
    ensure_results_dirs,
    evaluate_minimal_objective,
)


def gradient_descent_step(
    q: np.ndarray,
    T: np.ndarray,
    weights: dict[str, float],
    obs_config: dict[str, Any] | None = None,
    step_size: float = 1e-3,
    T_min: float = 1e-3,
    use_banded: bool = True,
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    q, T = _validate_problem_inputs(q, T, s=DEFAULT_S, k=DEFAULT_K)
    if step_size <= 0.0:
        raise ValueError("step_size must be positive.")
    if T_min <= 0.0:
        raise ValueError("T_min must be positive.")
    obs_cfg = default_obs_config() if obs_config is None else {**default_obs_config(), **dict(obs_config)}
    backward = backward_diff_banded if use_banded else backward_diff_dense
    current = backward(q, T, weights, obs_config=obs_cfg, config=config)
    grad_q = np.asarray(current["grad_q"], dtype=float)
    grad_T = np.asarray(current["grad_T"], dtype=float)
    obj_old = float(current["objective"]["value"])

    q_interior = np.asarray(q[1:-1], dtype=float)
    base_step = float(step_size)
    backtrack = float((config or {}).get("backtrack_factor", 0.5))
    max_backtracks = int((config or {}).get("max_backtracks", 12))
    max_rel_T_change = float((config or {}).get("max_rel_T_change", 0.15))

    q_new = q.copy()
    T_new = T.copy()
    obj_new = obj_old
    accepted = False
    used_step = base_step
    for _ in range(max_backtracks + 1):
        trial_q = np.concatenate(([q[0]], q_interior - used_step * grad_q, [q[-1]]))
        delta_T = -used_step * grad_T
        delta_T = np.clip(delta_T, -max_rel_T_change * T, max_rel_T_change * T)
        trial_T = np.maximum(T + delta_T, T_min)
        trial_obj = evaluate_minimal_objective(trial_q, trial_T, weights, obs_config=obs_cfg, config=config)
        if float(trial_obj["value"]) <= obj_old + 1e-12:
            q_new = trial_q
            T_new = trial_T
            obj_new = float(trial_obj["value"])
            accepted = True
            break
        used_step *= backtrack

    grad_norm = float(np.linalg.norm(np.concatenate((grad_q, grad_T))))
    return {
        "q_new": q_new,
        "T_new": T_new,
        "obj_old": obj_old,
        "obj_new": obj_new,
        "grad_norm": grad_norm,
        "step_used": float(used_step),
        "accepted": bool(accepted),
        "grad_q": grad_q,
        "grad_T": grad_T,
        "objective_old": current["objective"],
    }


def _plot_objective_history(history: dict[str, list[float]], save_path: str | Path) -> None:
    plt.figure(figsize=(8.4, 5.0))
    plt.plot(history["objective"], label="total objective")
    plt.plot(history["ctrl"], label="control")
    plt.plot(history["time"], label="time")
    plt.plot(history["obs"], label="obs")
    plt.xlabel("step")
    plt.ylabel("value")
    plt.title("Phase 9 objective curve")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_gradnorm(history: dict[str, list[float]], save_path: str | Path) -> None:
    plt.figure(figsize=(8.4, 4.8))
    plt.plot(history["grad_norm"], marker="o")
    plt.xlabel("step")
    plt.ylabel("gradient norm")
    plt.title("Phase 9 gradient norm curve")
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_traj_before_after(
    coeffs_before: np.ndarray,
    coeffs_after: np.ndarray,
    T_before: np.ndarray,
    T_after: np.ndarray,
    obs_config: dict[str, Any],
    save_path: str | Path,
) -> None:
    before = sample_trajectory(coeffs_before, T_before, num_per_segment=40, orders=(0,))
    after = sample_trajectory(coeffs_after, T_after, num_per_segment=40, orders=(0,))
    plt.figure(figsize=(8.8, 4.8))
    plt.plot(before["t_global"], before["order_0"], label="before")
    plt.plot(after["t_global"], after["order_0"], label="after")
    center = float(obs_config.get("center", 0.0))
    radius = float(obs_config.get("radius", 0.25))
    plt.axhspan(center - radius, center + radius, color="tab:red", alpha=0.12, label="obstacle band")
    plt.xlabel("global time")
    plt.ylabel("trajectory value")
    plt.title("Phase 9 trajectory before vs after")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_time_before_after(T_before: np.ndarray, T_after: np.ndarray, save_path: str | Path) -> None:
    indices = np.arange(1, T_before.size + 1)
    width = 0.38
    plt.figure(figsize=(8.4, 4.8))
    plt.bar(indices - width / 2.0, T_before, width=width, label="before")
    plt.bar(indices + width / 2.0, T_after, width=width, label="after")
    plt.xlabel("segment index")
    plt.ylabel("duration")
    plt.title("Phase 9 time allocation before vs after")
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_obj_delta(history: dict[str, list[float]], save_path: str | Path) -> None:
    objective = np.asarray(history["objective"], dtype=float)
    delta = objective[:-1] - objective[1:]
    plt.figure(figsize=(8.4, 4.8))
    plt.bar(np.arange(delta.size), delta)
    plt.xlabel("step")
    plt.ylabel("objective decrease")
    plt.title("Phase 9 objective decrease per step")
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def run_minimal_optimization_demo(
    q0: np.ndarray,
    T0: np.ndarray,
    weights: dict[str, float],
    obs_config: dict[str, Any] | None = None,
    n_steps: int = 50,
    step_size: float = 1e-3,
    T_min: float = 1e-3,
    use_banded: bool = True,
    save_dir: str | Path | None = None,
    seed: int = 42,
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    q, T = _validate_problem_inputs(q0, T0, s=DEFAULT_S, k=DEFAULT_K)
    obs_cfg = default_obs_config() if obs_config is None else {**default_obs_config(), **dict(obs_config)}
    results_dirs = ensure_results_dirs(DEFAULT_RESULTS_DIR if save_dir is None else save_dir)
    opt_dir = results_dirs["opt"]
    weights = {**default_weights(), **dict(weights)}

    history: dict[str, list[float]] = {
        "objective": [],
        "ctrl": [],
        "time": [],
        "obs": [],
        "grad_norm": [],
        "step_used": [],
        "accepted": [],
        "T_min": [],
    }
    q_history = [q.copy()]
    T_history = [T.copy()]

    current = evaluate_minimal_objective(q, T, weights, obs_config=obs_cfg, config=config)
    history["objective"].append(float(current["value"]))
    history["ctrl"].append(float(current["parts"]["ctrl"]))
    history["time"].append(float(current["parts"]["time"]))
    history["obs"].append(float(current["parts"]["obs"]))
    history["grad_norm"].append(float("nan"))
    history["step_used"].append(0.0)
    history["accepted"].append(True)
    history["T_min"].append(float(np.min(T)))

    for _ in range(int(n_steps)):
        step = gradient_descent_step(
            q,
            T,
            weights,
            obs_config=obs_cfg,
            step_size=step_size,
            T_min=T_min,
            use_banded=use_banded,
            config=config,
        )
        q = step["q_new"]
        T = step["T_new"]
        current = evaluate_minimal_objective(q, T, weights, obs_config=obs_cfg, config=config)
        history["objective"].append(float(current["value"]))
        history["ctrl"].append(float(current["parts"]["ctrl"]))
        history["time"].append(float(current["parts"]["time"]))
        history["obs"].append(float(current["parts"]["obs"]))
        history["grad_norm"].append(float(step["grad_norm"]))
        history["step_used"].append(float(step["step_used"]))
        history["accepted"].append(bool(step["accepted"]))
        history["T_min"].append(float(np.min(T)))
        q_history.append(q.copy())
        T_history.append(T.copy())

    coeffs_before = evaluate_minimal_objective(q0, T0, weights, obs_config=obs_cfg, config=config)["coeff_blocks"]
    coeffs_after = current["coeff_blocks"]
    _plot_objective_history(history, opt_dir / "phase9_objective_curve.png")
    _plot_gradnorm(history, opt_dir / "phase9_gradnorm_curve.png")
    _plot_traj_before_after(coeffs_before, coeffs_after, T0, T, obs_cfg, opt_dir / "phase9_traj_before_after.png")
    _plot_time_before_after(T0, T, opt_dir / "phase9_time_before_after.png")
    _plot_obj_delta(history, opt_dir / "phase9_obj_delta_per_step.png")

    rows = []
    for idx in range(len(history["objective"])):
        rows.append(
            [
                idx,
                history["objective"][idx],
                history["ctrl"][idx],
                history["time"][idx],
                history["obs"][idx],
                history["grad_norm"][idx],
                history["step_used"][idx],
                history["accepted"][idx],
                history["T_min"][idx],
            ]
        )
    _save_csv(
        opt_dir / "phase9_opt_history.csv",
        ["step", "objective", "ctrl", "time", "obs", "grad_norm", "step_used", "accepted", "T_min"],
        rows,
    )

    summary = {
        "objective_drop": float(history["objective"][0] - history["objective"][-1]),
        "grad_norm_drop": float(
            (history["grad_norm"][1] if len(history["grad_norm"]) > 1 and np.isfinite(history["grad_norm"][1]) else 0.0)
            - (history["grad_norm"][-1] if np.isfinite(history["grad_norm"][-1]) else 0.0)
        ),
        "ctrl_drop": float(history["ctrl"][0] - history["ctrl"][-1]),
        "time_drop": float(history["time"][0] - history["time"][-1]),
        "obs_drop": float(history["obs"][0] - history["obs"][-1]),
        "touches_T_min": bool(np.any(np.asarray(history["T_min"], dtype=float) <= T_min + 1e-12)),
        "final_objective": float(history["objective"][-1]),
        "initial_objective": float(history["objective"][0]),
    }
    _save_json(opt_dir / "phase9_opt_summary.json", summary)
    lines = [
        "# Phase 9 Optimization Summary",
        "",
        f"- Objective drop: `{summary['objective_drop']:.6e}`.",
        f"- Control drop: `{summary['ctrl_drop']:.6e}`.",
        f"- Time drop: `{summary['time_drop']:.6e}`.",
        f"- Obstacle drop: `{summary['obs_drop']:.6e}`.",
        f"- Gradient-norm drop proxy: `{summary['grad_norm_drop']:.6e}`.",
        f"- Any duration touched the lower bound: `{summary['touches_T_min']}`.",
        "",
        "Interpretation:",
        "- A negative objective drop would indicate the loop is unstable or the step size is too large.",
        "- When the obstacle term decreases while durations stay feasible, the minimal loop is already useful as an optimization primitive.",
        "- If durations repeatedly hit the lower bound, the next phase should add better time-parameterization handling rather than a larger objective.",
    ]
    (opt_dir / "phase9_opt_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {
        "history": history,
        "q_history": q_history,
        "T_history": T_history,
        "summary": summary,
        "coeffs_before": coeffs_before,
        "coeffs_after": coeffs_after,
    }


def main() -> None:
    q0 = np.asarray([0.0, 0.6, -0.2, 0.85, -0.1, 0.55, 0.0], dtype=float)
    T0 = np.asarray([0.9, 1.1, 0.85, 1.2, 0.95, 1.05], dtype=float)
    run_minimal_optimization_demo(q0, T0, default_weights(), n_steps=25, step_size=5e-3, T_min=0.2)


if __name__ == "__main__":
    main()
