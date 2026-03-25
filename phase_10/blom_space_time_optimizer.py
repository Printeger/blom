"""
Phase 10 space-time optimizer on xi = (q_bar, tau).
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from phase_1.minco_scalar_baseline import evaluate_trajectory
from phase_10.blom_full_backward_diff import (
    DEFAULT_RESULTS_DIR,
    DEFAULT_S,
    default_bc_config,
    default_dyn_config,
    default_obs_config,
    default_reg_config,
    default_weights,
    ensure_results_dirs,
    evaluate_full_objective,
    full_backward_diff_reparam,
    q_tau_to_xi,
    qbar_to_full,
    save_csv,
    save_json,
    tau_to_T,
    xi_to_qbar_tau,
)


def _infer_M_from_xi(xi: np.ndarray) -> int:
    xi = np.asarray(xi, dtype=float).reshape(-1)
    if xi.size < 3 or xi.size % 2 == 0:
        raise ValueError("xi must have odd length 2*M-1 with M >= 2.")
    return (xi.size + 1) // 2


def _boundaries_from_config(config: dict[str, Any] | None) -> tuple[float, float]:
    cfg = {} if config is None else dict(config)
    if "q_start" not in cfg or "q_end" not in cfg:
        raise ValueError("optimizer config must provide q_start and q_end for xi-based optimization.")
    return float(cfg["q_start"]), float(cfg["q_end"])


def _objective_from_xi(
    xi: np.ndarray,
    weights: dict[str, float],
    *,
    T_min: float,
    obs_config: dict[str, Any] | None,
    dyn_config: dict[str, Any] | None,
    bc_config: dict[str, Any] | None,
    reg_config: dict[str, Any] | None,
    s: int,
    k: int,
    config: dict[str, Any] | None,
) -> dict[str, Any]:
    q_start, q_end = _boundaries_from_config(config)
    q_bar, tau = xi_to_qbar_tau(xi)
    q = qbar_to_full(q_bar, q_start, q_end)
    T = tau_to_T(tau, T_min=T_min)
    return evaluate_full_objective(q, T, weights, obs_config=obs_config, dyn_config=dyn_config, bc_config=bc_config, reg_config=reg_config, s=s, k=k)


def optimizer_step(
    xi: np.ndarray,
    weights: dict[str, float],
    T_min: float = 1e-3,
    obs_config: dict[str, Any] | None = None,
    dyn_config: dict[str, Any] | None = None,
    bc_config: dict[str, Any] | None = None,
    reg_config: dict[str, Any] | None = None,
    method: str = "gd",
    step_rule: str = "armijo",
    step_size: float = 1e-3,
    s: int = 4,
    k: int = 2,
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    if method != "gd":
        raise NotImplementedError("Phase 10 currently implements gradient descent only.")
    if step_size <= 0.0:
        raise ValueError("step_size must be positive.")
    cfg = {} if config is None else dict(config)
    q_start, q_end = _boundaries_from_config(cfg)

    q_bar, tau = xi_to_qbar_tau(xi)
    q = qbar_to_full(q_bar, q_start, q_end)
    backward = full_backward_diff_reparam(
        q,
        tau,
        weights,
        T_min=T_min,
        obs_config=obs_config,
        dyn_config=dyn_config,
        bc_config=bc_config,
        reg_config=reg_config,
        s=s,
        k=k,
        config=cfg,
    )
    grad = np.asarray(backward["grad_xi"], dtype=float)
    obj_old = float(backward["objective"]["value"])
    grad_norm = float(np.linalg.norm(grad))

    base_step = float(step_size)
    backtrack = float(cfg.get("armijo_backtrack", 0.5))
    armijo_c = float(cfg.get("armijo_c", 1e-4))
    max_backtracks = int(cfg.get("max_backtracks", 16))

    accepted = False
    xi_new = np.asarray(xi, dtype=float).copy()
    obj_new = obj_old
    used_step = base_step

    if step_rule == "fixed":
        candidate = np.asarray(xi, dtype=float) - used_step * grad
        trial = _objective_from_xi(candidate, weights, T_min=T_min, obs_config=obs_config, dyn_config=dyn_config, bc_config=bc_config, reg_config=reg_config, s=s, k=k, config=cfg)
        xi_new = candidate
        obj_new = float(trial["value"])
        accepted = bool(obj_new <= obj_old + 1e-12)
    elif step_rule == "armijo":
        for _ in range(max_backtracks + 1):
            candidate = np.asarray(xi, dtype=float) - used_step * grad
            trial = _objective_from_xi(candidate, weights, T_min=T_min, obs_config=obs_config, dyn_config=dyn_config, bc_config=bc_config, reg_config=reg_config, s=s, k=k, config=cfg)
            if float(trial["value"]) <= obj_old - armijo_c * used_step * float(grad @ grad):
                xi_new = candidate
                obj_new = float(trial["value"])
                accepted = True
                break
            used_step *= backtrack
        if not accepted:
            xi_new = np.asarray(xi, dtype=float).copy()
            obj_new = obj_old
    else:
        raise ValueError(f"Unsupported step_rule {step_rule!r}.")

    q_bar_new, tau_new = xi_to_qbar_tau(xi_new)
    q_new = qbar_to_full(q_bar_new, q_start, q_end)
    T_new = tau_to_T(tau_new, T_min=T_min)
    return {
        "xi_new": xi_new,
        "q_new": q_new,
        "tau_new": tau_new,
        "T_new": T_new,
        "obj_old": obj_old,
        "obj_new": obj_new,
        "grad_norm": grad_norm,
        "step_size_used": used_step,
        "accepted": bool(accepted),
    }


def _plot_objective_curve(history: dict[str, list[float]], save_path: str | Path) -> None:
    plt.figure(figsize=(8.8, 5.2))
    plt.plot(history["objective"], label="total")
    for key in ["ctrl", "time", "obs", "dyn", "bc", "reg"]:
        plt.plot(history[key], label=key)
    plt.xlabel("step")
    plt.ylabel("objective")
    plt.title("Phase 10 objective curve")
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=8, ncols=2)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_gradnorm(history: dict[str, list[float]], save_path: str | Path) -> None:
    plt.figure(figsize=(8.4, 4.8))
    plt.plot(history["grad_norm"], marker="o")
    plt.xlabel("step")
    plt.ylabel("gradient norm")
    plt.title("Phase 10 gradient norm curve")
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_step_size(history: dict[str, list[float]], save_path: str | Path) -> None:
    plt.figure(figsize=(8.4, 4.8))
    plt.plot(history["step_size"], marker="o")
    plt.xlabel("step")
    plt.ylabel("step size")
    plt.title("Phase 10 step size curve")
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_traj_progress(coeffs_hist: list[np.ndarray], T_hist: list[np.ndarray], obs_config: dict[str, Any], save_path: str | Path) -> None:
    indices = sorted(set([0, len(coeffs_hist) // 2, len(coeffs_hist) - 1]))
    plt.figure(figsize=(9.0, 5.0))
    for idx in indices:
        T_local = np.asarray(T_hist[idx], dtype=float)
        coeff_local = np.asarray(coeffs_hist[idx], dtype=float)
        total_time = float(np.sum(T_local))
        t_global = np.linspace(0.0, total_time, max(2, 40 * T_local.size), endpoint=True)
        t_global[-1] = np.nextafter(total_time, 0.0)
        values = np.asarray([evaluate_trajectory(coeff_local, T_local, float(t), order=0) for t in t_global], dtype=float)
        label = "start" if idx == 0 else ("end" if idx == len(coeffs_hist) - 1 else f"step {idx}")
        plt.plot(t_global, values, label=label)
    center = float(obs_config.get("center", 0.0))
    radius = float(obs_config.get("radius", 0.25))
    plt.axhspan(center - radius, center + radius, color="tab:red", alpha=0.12, label="obstacle band")
    plt.xlabel("global time")
    plt.ylabel("trajectory value")
    plt.title("Phase 10 trajectory progress")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_time_progress(T_before: np.ndarray, T_after: np.ndarray, save_path: str | Path) -> None:
    idx = np.arange(1, T_before.size + 1)
    width = 0.38
    plt.figure(figsize=(8.8, 4.8))
    plt.bar(idx - width / 2.0, T_before, width=width, label="before")
    plt.bar(idx + width / 2.0, T_after, width=width, label="after")
    plt.xlabel("segment index")
    plt.ylabel("duration")
    plt.title("Phase 10 time allocation before vs after")
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
    plt.title("Phase 10 objective decrease per step")
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def run_space_time_optimization(
    q0: np.ndarray,
    tau0: np.ndarray,
    weights: dict[str, float],
    obs_config: dict[str, Any] | None = None,
    dyn_config: dict[str, Any] | None = None,
    bc_config: dict[str, Any] | None = None,
    reg_config: dict[str, Any] | None = None,
    method: str = "gd",
    step_rule: str = "armijo",
    n_steps: int = 100,
    step_size: float = 1e-3,
    T_min: float = 1e-3,
    s: int = 4,
    k: int = 2,
    config: dict[str, Any] | None = None,
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    del seed
    q0 = np.asarray(q0, dtype=float).reshape(-1)
    tau0 = np.asarray(tau0, dtype=float).reshape(-1)
    if tau0.size != q0.size - 1:
        raise ValueError("tau0 must have one entry per segment.")

    weights = {**default_weights(), **dict(weights)}
    obs_cfg = default_obs_config() if obs_config is None else {**default_obs_config(), **dict(obs_config)}
    dyn_cfg = default_dyn_config() if dyn_config is None else {**default_dyn_config(), **dict(dyn_config)}
    bc_cfg = default_bc_config(q0) if bc_config is None else {**default_bc_config(q0), **dict(bc_config)}
    reg_cfg = default_reg_config() if reg_config is None else {**default_reg_config(), **dict(reg_config)}
    cfg = {} if config is None else dict(config)
    cfg["q_start"] = float(q0[0])
    cfg["q_end"] = float(q0[-1])
    results_dirs = ensure_results_dirs(DEFAULT_RESULTS_DIR if save_dir is None else save_dir)
    opt_dir = results_dirs["optimizer"]

    xi = q_tau_to_xi(q0, tau0)
    q = q0.copy()
    tau = tau0.copy()
    T = tau_to_T(tau, T_min=T_min)
    current = evaluate_full_objective(q, T, weights, obs_config=obs_cfg, dyn_config=dyn_cfg, bc_config=bc_cfg, reg_config=reg_cfg, s=s, k=k)

    history = {key: [] for key in ["objective", "ctrl", "time", "obs", "dyn", "bc", "reg", "grad_norm", "step_size", "runtime_sec"]}
    accepted_hist = []
    coeffs_hist = []
    T_hist = []

    for _step in range(int(n_steps) + 1):
        current = evaluate_full_objective(q, T, weights, obs_config=obs_cfg, dyn_config=dyn_cfg, bc_config=bc_cfg, reg_config=reg_cfg, s=s, k=k)
        history["objective"].append(float(current["value"]))
        for key in ["ctrl", "time", "obs", "dyn", "bc", "reg"]:
            history[key].append(float(current["parts"][key]))
        coeffs_hist.append(np.asarray(current["coeff_blocks"], dtype=float))
        T_hist.append(T.copy())
        if _step == int(n_steps):
            history["grad_norm"].append(float("nan"))
            history["step_size"].append(0.0)
            history["runtime_sec"].append(0.0)
            accepted_hist.append(True)
            break
        start = time.perf_counter()
        step = optimizer_step(
            xi,
            weights,
            T_min=T_min,
            obs_config=obs_cfg,
            dyn_config=dyn_cfg,
            bc_config=bc_cfg,
            reg_config=reg_cfg,
            method=method,
            step_rule=step_rule,
            step_size=step_size,
            s=s,
            k=k,
            config=cfg,
        )
        elapsed = time.perf_counter() - start
        xi = np.asarray(step["xi_new"], dtype=float)
        q = np.asarray(step["q_new"], dtype=float)
        tau = np.asarray(step["tau_new"], dtype=float)
        T = np.asarray(step["T_new"], dtype=float)
        history["grad_norm"].append(float(step["grad_norm"]))
        history["step_size"].append(float(step["step_size_used"]))
        history["runtime_sec"].append(float(elapsed))
        accepted_hist.append(bool(step["accepted"]))

    _plot_objective_curve(history, opt_dir / "phase10_objective_curve.png")
    _plot_gradnorm(history, opt_dir / "phase10_gradnorm_curve.png")
    _plot_step_size(history, opt_dir / "phase10_stepsize_curve.png")
    _plot_traj_progress(coeffs_hist, T_hist, obs_cfg, opt_dir / "phase10_traj_progress.png")
    _plot_time_progress(T_hist[0], T_hist[-1], opt_dir / "phase10_time_progress.png")
    _plot_obj_delta(history, opt_dir / "phase10_obj_delta_per_step.png")

    rows = []
    for idx in range(len(history["objective"])):
        rows.append(
            [
                idx,
                history["objective"][idx],
                history["ctrl"][idx],
                history["time"][idx],
                history["obs"][idx],
                history["dyn"][idx],
                history["bc"][idx],
                history["reg"][idx],
                history["grad_norm"][idx],
                history["step_size"][idx],
                history["runtime_sec"][idx],
                accepted_hist[idx],
            ]
        )
    save_csv(
        opt_dir / "phase10_opt_history.csv",
        ["step", "objective", "ctrl", "time", "obs", "dyn", "bc", "reg", "grad_norm", "step_size", "runtime_sec", "accepted"],
        rows,
    )
    summary = {
        "objective_initial": float(history["objective"][0]),
        "objective_final": float(history["objective"][-1]),
        "objective_drop": float(history["objective"][0] - history["objective"][-1]),
        "grad_norm_peak": float(np.nanmax(history["grad_norm"])),
        "grad_norm_final": float(history["grad_norm"][-2]) if len(history["grad_norm"]) >= 2 else float("nan"),
        "accepted_ratio": float(np.mean(accepted_hist)),
        "touches_T_min": bool(np.any(np.asarray(T_hist[-1], dtype=float) <= T_min + 1e-12)),
        "runtime_total": float(np.sum(history["runtime_sec"])),
        "runtime_per_iter": float(np.mean(history["runtime_sec"][1:])) if len(history["runtime_sec"]) > 1 else 0.0,
    }
    save_json(opt_dir / "phase10_opt_summary.json", {"summary": summary, "config": {"n_steps": int(n_steps), "step_rule": step_rule, "method": method, "T_min": float(T_min)}})
    Path(opt_dir / "phase10_opt_summary.md").write_text(
        "\n".join(
            [
                "# Phase 10 Optimization Summary",
                "",
                f"- Total objective decreased: `{summary['objective_drop'] > 0.0}`.",
                f"- Objective drop: `{summary['objective_drop']:.6e}`.",
                f"- Accepted-step ratio: `{summary['accepted_ratio']:.6f}`.",
                f"- Any duration touched the floor: `{summary['touches_T_min']}`.",
                f"- Runtime per iteration: `{summary['runtime_per_iter']:.6e}` sec.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return {
        "summary": summary,
        "history": history,
        "accepted": accepted_hist,
        "q_final": q,
        "tau_final": tau,
        "T_final": T,
        "coeffs_final": coeffs_hist[-1],
        "weights": weights,
    }


def main() -> None:
    from phase_10.blom_full_backward_diff import representative_case, T_to_tau

    q0, T0 = representative_case()
    tau0 = T_to_tau(T0, T_min=0.2)
    run_space_time_optimization(q0, tau0, default_weights(), n_steps=25, step_size=5e-3, T_min=0.2)


if __name__ == "__main__":
    main()
