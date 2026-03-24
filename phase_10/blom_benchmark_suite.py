"""
Phase 10 benchmark and ablation framework.
"""

from __future__ import annotations

import math
import time
import tracemalloc
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from phase_1.minco_scalar_baseline import solve_minco_coefficients
from phase_4.utils.hermite_utils import hermite_reconstruct_center_segment
from phase_5.blom_boundary_jump_check import assemble_scheme_A
from phase_10.blom_full_backward_diff import (
    DEFAULT_K,
    DEFAULT_RESULTS_DIR,
    DEFAULT_S,
    T_to_tau,
    compute_raw_schemeC_jacobians_general_k,
    compute_support_width_stats,
    default_bc_config,
    default_dyn_config,
    default_obs_config,
    default_reg_config,
    default_weights,
    ensure_results_dirs,
    evaluate_full_objective,
    evaluate_full_objective_from_coeffs,
    representative_case,
    save_csv,
    save_json,
)
from phase_10.blom_space_time_optimizer import run_space_time_optimization


def _sample_case(M: int, regime_box: tuple[float, float], seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    increments = rng.normal(loc=0.0, scale=0.9, size=M)
    q = np.zeros((M + 1,), dtype=float)
    q[1:] = np.cumsum(increments)
    T = rng.uniform(float(regime_box[0]), float(regime_box[1]), size=M).astype(float)
    return q, T


def _measure_memory_and_time(fn, *args, **kwargs) -> tuple[Any, float, float]:
    tracemalloc.start()
    start = time.perf_counter()
    result = fn(*args, **kwargs)
    elapsed = time.perf_counter() - start
    _current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return result, float(elapsed), float(peak)


def _support_stats_or_default(q: np.ndarray, T: np.ndarray, k: int) -> dict[str, float]:
    try:
        jac = compute_raw_schemeC_jacobians_general_k(q, T, k=k)
        return compute_support_width_stats(jac, T.size)
    except NotImplementedError:
        return {"support_width_mean": float("nan"), "support_width_max": float("nan")}


def _constraint_violation(
    coeff_blocks: np.ndarray,
    T: np.ndarray,
    *,
    q_boundary: tuple[float, float],
    obs_config: dict[str, Any],
    dyn_config: dict[str, Any],
    bc_config: dict[str, Any],
) -> float:
    max_violation = 0.0
    center = float(obs_config.get("center", 0.0))
    radius = float(obs_config.get("radius", 0.25))
    alphas_obs = np.asarray(obs_config.get("sample_alphas", [0.5]), dtype=float)
    for seg_idx, duration in enumerate(T):
        coeff = coeff_blocks[seg_idx]
        for alpha in alphas_obs:
            value = evaluate_full_objective_from_coeffs(
                coeff_blocks[seg_idx : seg_idx + 1],
                np.asarray([duration], dtype=float),
                {"lambda_T": 0.0, "lambda_obs": 0.0, "lambda_dyn": 0.0, "lambda_bc": 0.0, "lambda_reg": 0.0},
                q_boundary=q_boundary,
                obs_config=obs_config,
                dyn_config=dyn_config,
                bc_config=bc_config,
                reg_config=default_reg_config(),
            )
            del value
            x = float(np.asarray([sum(coeff[p] * ((alpha * duration) ** p) for p in range(coeff.size))], dtype=float)[0])
            max_violation = max(max_violation, max(radius - abs(x - center), 0.0))
    orders = [int(order) for order in dyn_config.get("derivative_orders", [1, 2, 3])]
    limits = np.asarray(dyn_config.get("limits", np.ones(len(orders))), dtype=float).reshape(-1)
    if limits.size == 1:
        limits = np.full((len(orders),), float(limits[0]), dtype=float)
    alphas_dyn = np.asarray(dyn_config.get("sample_alphas", [0.5]), dtype=float)
    for seg_idx, duration in enumerate(T):
        coeff = coeff_blocks[seg_idx]
        for alpha in alphas_dyn:
            tau = float(alpha * duration)
            for order, limit in zip(orders, limits):
                value = abs(sum(coeff[p] * math.factorial(p) / math.factorial(p - order) * (tau ** (p - order)) for p in range(order, coeff.size)))
                max_violation = max(max_violation, max(value - limit, 0.0))
    desired_start = np.asarray(bc_config.get("desired_start"), dtype=float)
    desired_end = np.asarray(bc_config.get("desired_end"), dtype=float)
    for order in range(desired_start.size):
        start_val = sum(coeff_blocks[0, p] * math.factorial(p) / math.factorial(p - order) * (0.0 ** max(p - order, 0)) for p in range(order, coeff_blocks.shape[1]))
        end_val = sum(
            coeff_blocks[-1, p] * math.factorial(p) / math.factorial(p - order) * (float(T[-1]) ** (p - order))
            for p in range(order, coeff_blocks.shape[1])
        )
        max_violation = max(max_violation, abs(start_val - desired_start[order]), abs(end_val - desired_end[order]))
    return float(max_violation)


def _heuristic_coeffs(q: np.ndarray, T: np.ndarray) -> np.ndarray:
    coeffs = np.zeros((T.size, 2 * DEFAULT_S), dtype=float)
    zeros = np.zeros((3,), dtype=float)
    for idx, duration in enumerate(T):
        coeffs[idx] = hermite_reconstruct_center_segment(float(q[idx]), float(q[idx + 1]), zeros, zeros, float(duration))["coeff"]
    return coeffs


def run_baseline_raw_schemeC(
    q0: np.ndarray,
    T0: np.ndarray,
    weights: dict[str, float],
    *,
    k: int = DEFAULT_K,
    n_steps: int = 40,
    T_min: float = 0.2,
    obs_config: dict[str, Any] | None = None,
    dyn_config: dict[str, Any] | None = None,
    bc_config: dict[str, Any] | None = None,
    reg_config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    tau0 = T_to_tau(T0, T_min=T_min)
    result, elapsed, peak = _measure_memory_and_time(
        run_space_time_optimization,
        q0,
        tau0,
        weights,
        obs_config=obs_config,
        dyn_config=dyn_config,
        bc_config=bc_config,
        reg_config=reg_config,
        n_steps=n_steps,
        step_size=5e-3,
        T_min=T_min,
        k=k,
    )
    return {
        "method": "raw_schemeC",
        "final_objective": float(result["summary"]["objective_final"]),
        "runtime_total": elapsed,
        "runtime_per_iter": float(result["summary"]["runtime_per_iter"]),
        "memory_peak": peak,
        "constraint_violation": float("nan"),
        "history": result["history"],
        "support_stats": _support_stats_or_default(q0, T0, k),
    }


def run_baseline_minco(
    q0: np.ndarray,
    T0: np.ndarray,
    weights: dict[str, float],
    *,
    obs_config: dict[str, Any] | None = None,
    dyn_config: dict[str, Any] | None = None,
    bc_config: dict[str, Any] | None = None,
    reg_config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    zeta_start = np.zeros((DEFAULT_S,), dtype=float)
    zeta_end = np.zeros((DEFAULT_S,), dtype=float)
    zeta_start[0] = float(q0[0])
    zeta_end[0] = float(q0[-1])
    solved, elapsed, peak = _measure_memory_and_time(solve_minco_coefficients, q0, T0, zeta_start, zeta_end)
    objective = evaluate_full_objective_from_coeffs(
        np.asarray(solved["coeffs"], dtype=float),
        T0,
        weights,
        q_boundary=(float(q0[0]), float(q0[-1])),
        obs_config=default_obs_config() if obs_config is None else obs_config,
        dyn_config=default_dyn_config() if dyn_config is None else dyn_config,
        bc_config=default_bc_config(q0) if bc_config is None else bc_config,
        reg_config=default_reg_config() if reg_config is None else reg_config,
    )
    return {
        "method": "minco",
        "final_objective": float(objective["value"]),
        "runtime_total": elapsed,
        "runtime_per_iter": elapsed,
        "memory_peak": peak,
        "constraint_violation": float("nan"),
        "history": {"objective": [float(objective["value"])]},
        "support_stats": {"support_width_mean": float(T0.size), "support_width_max": float(T0.size)},
    }


def run_baseline_schemeA(
    q0: np.ndarray,
    T0: np.ndarray,
    weights: dict[str, float],
    *,
    obs_config: dict[str, Any] | None = None,
    dyn_config: dict[str, Any] | None = None,
    bc_config: dict[str, Any] | None = None,
    reg_config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    result, elapsed, peak = _measure_memory_and_time(assemble_scheme_A, q0, T0, DEFAULT_S, DEFAULT_K)
    objective = evaluate_full_objective_from_coeffs(
        np.asarray(result["coeffs"], dtype=float),
        T0,
        weights,
        q_boundary=(float(q0[0]), float(q0[-1])),
        obs_config=default_obs_config() if obs_config is None else obs_config,
        dyn_config=default_dyn_config() if dyn_config is None else dyn_config,
        bc_config=default_bc_config(q0) if bc_config is None else bc_config,
        reg_config=default_reg_config() if reg_config is None else reg_config,
    )
    return {
        "method": "schemeA",
        "final_objective": float(objective["value"]),
        "runtime_total": elapsed,
        "runtime_per_iter": elapsed,
        "memory_peak": peak,
        "constraint_violation": float("nan"),
        "history": {"objective": [float(objective["value"])]},
        "support_stats": {"support_width_mean": float(T0.size), "support_width_max": float(T0.size)},
    }


def run_baseline_heuristic(
    q0: np.ndarray,
    T0: np.ndarray,
    weights: dict[str, float],
    *,
    obs_config: dict[str, Any] | None = None,
    dyn_config: dict[str, Any] | None = None,
    bc_config: dict[str, Any] | None = None,
    reg_config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    coeffs, elapsed, peak = _measure_memory_and_time(_heuristic_coeffs, q0, T0)
    objective = evaluate_full_objective_from_coeffs(
        coeffs,
        T0,
        weights,
        q_boundary=(float(q0[0]), float(q0[-1])),
        obs_config=default_obs_config() if obs_config is None else obs_config,
        dyn_config=default_dyn_config() if dyn_config is None else dyn_config,
        bc_config=default_bc_config(q0) if bc_config is None else bc_config,
        reg_config=default_reg_config() if reg_config is None else reg_config,
    )
    return {
        "method": "heuristic",
        "final_objective": float(objective["value"]),
        "runtime_total": elapsed,
        "runtime_per_iter": elapsed,
        "memory_peak": peak,
        "constraint_violation": float("nan"),
        "history": {"objective": [float(objective["value"])]},
        "support_stats": {"support_width_mean": 1.0, "support_width_max": 1.0},
    }


def _compute_objective_gaps(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[int, int], list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault((int(row["M"]), int(row["k"])), []).append(row)
    for key, items in grouped.items():
        ref_val = min(float(item["final_objective"]) for item in items if np.isfinite(item["final_objective"]))
        for item in items:
            item["objective_gap"] = (float(item["final_objective"]) - ref_val) / max(abs(ref_val), 1e-12)
    return rows


def _plot_runtime_vs_M(rows: list[dict[str, Any]], save_path: str | Path) -> None:
    plt.figure(figsize=(8.8, 5.0))
    methods = sorted({row["method"] for row in rows})
    for method in methods:
        items = [row for row in rows if row["method"] == method]
        plt.plot([row["M"] for row in items], [row["runtime_total"] for row in items], marker="o", label=method)
    plt.xlabel("M")
    plt.ylabel("runtime total (sec)")
    plt.title("Phase 10 runtime vs M")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_obj_gap_vs_k(rows: list[dict[str, Any]], save_path: str | Path) -> None:
    plt.figure(figsize=(8.8, 5.0))
    methods = sorted({row["method"] for row in rows})
    for method in methods:
        items = [row for row in rows if row["method"] == method]
        plt.plot([row["k"] for row in items], [row["objective_gap"] for row in items], marker="o", label=method)
    plt.xlabel("k")
    plt.ylabel("objective gap")
    plt.title("Phase 10 objective gap vs k")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_memory_vs_k(rows: list[dict[str, Any]], save_path: str | Path) -> None:
    plt.figure(figsize=(8.8, 5.0))
    methods = sorted({row["method"] for row in rows})
    for method in methods:
        items = [row for row in rows if row["method"] == method]
        plt.plot([row["k"] for row in items], [row["memory_peak"] / (1024.0**2) for row in items], marker="o", label=method)
    plt.xlabel("k")
    plt.ylabel("peak memory (MB)")
    plt.title("Phase 10 memory vs k")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_baseline_compare(rows: list[dict[str, Any]], save_path: str | Path) -> None:
    methods = [row["method"] for row in rows]
    values = [row["final_objective"] for row in rows]
    plt.figure(figsize=(8.8, 4.8))
    plt.bar(methods, values)
    plt.ylabel("final objective")
    plt.title("Phase 10 baseline comparison")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_ablation_heatmap(rows: list[dict[str, Any]], save_path: str | Path) -> None:
    labels = [row["ablation_name"] for row in rows]
    matrix = np.asarray([[row["final_objective"], row["runtime_total"], row["objective_gap"]] for row in rows], dtype=float)
    plt.figure(figsize=(8.8, 4.8))
    im = plt.imshow(matrix, aspect="auto", cmap="viridis")
    plt.colorbar(im, label="value")
    plt.xticks(np.arange(3), ["objective", "runtime", "obj gap"])
    plt.yticks(np.arange(len(labels)), labels)
    plt.title("Phase 10 ablation heatmap")
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_quality_speed_pareto(rows: list[dict[str, Any]], save_path: str | Path) -> None:
    plt.figure(figsize=(8.8, 5.0))
    for row in rows:
        plt.scatter(row["runtime_total"], row["final_objective"], label=row["method"])
        plt.text(row["runtime_total"], row["final_objective"], row["method"], fontsize=8)
    plt.xlabel("runtime total (sec)")
    plt.ylabel("final objective")
    plt.title("Phase 10 quality-speed Pareto")
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_sparsity_vs_k(rows: list[dict[str, Any]], save_path: str | Path) -> None:
    items = [row for row in rows if row["method"] == "raw_schemeC"]
    plt.figure(figsize=(8.8, 5.0))
    plt.plot([row["k"] for row in items], [row["support_width_mean"] for row in items], marker="o", label="mean support width")
    plt.plot([row["k"] for row in items], [row["support_width_max"] for row in items], marker="s", label="max support width")
    plt.xlabel("k")
    plt.ylabel("support width")
    plt.title("Phase 10 sparsity vs k")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def run_benchmark_suite(
    q0: np.ndarray | None = None,
    T0: np.ndarray | None = None,
    weights: dict[str, float] | None = None,
    obs_config: dict[str, Any] | None = None,
    dyn_config: dict[str, Any] | None = None,
    bc_config: dict[str, Any] | None = None,
    reg_config: dict[str, Any] | None = None,
    M_values: list[int] | None = None,
    k_values: list[int] | None = None,
    obs_weight_values: list[float] | None = None,
    time_boxes: list[tuple[float, float]] | None = None,
    benchmark_methods: list[str] | None = None,
    n_steps: int = 30,
    T_min: float = 0.2,
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    q_rep, T_rep = representative_case() if q0 is None or T0 is None else (np.asarray(q0, dtype=float), np.asarray(T0, dtype=float))
    weights = {**default_weights(), **({} if weights is None else dict(weights))}
    obs_cfg = default_obs_config() if obs_config is None else {**default_obs_config(), **dict(obs_config)}
    dyn_cfg = default_dyn_config() if dyn_config is None else {**default_dyn_config(), **dict(dyn_config)}
    bc_cfg = default_bc_config(q_rep) if bc_config is None else {**default_bc_config(q_rep), **dict(bc_config)}
    reg_cfg = default_reg_config() if reg_config is None else {**default_reg_config(), **dict(reg_config)}
    M_values = [10, 20, 40, 80] if M_values is None else list(M_values)
    k_values = [2, 4, 6, 8] if k_values is None else list(k_values)
    obs_weight_values = [0.0, 0.1, 1.0, 10.0] if obs_weight_values is None else list(obs_weight_values)
    time_boxes = [(0.9, 1.1), (0.5, 2.0), (0.2, 3.0)] if time_boxes is None else list(time_boxes)
    benchmark_methods = ["raw_schemeC", "minco", "schemeA", "heuristic"] if benchmark_methods is None else list(benchmark_methods)
    results_dirs = ensure_results_dirs(DEFAULT_RESULTS_DIR if save_dir is None else save_dir)

    baseline_runners = {
        "raw_schemeC": lambda q, T, w, k: run_baseline_raw_schemeC(q, T, w, k=k, n_steps=n_steps, T_min=T_min, obs_config=obs_cfg, dyn_config=dyn_cfg, bc_config=bc_cfg, reg_config=reg_cfg),
        "minco": lambda q, T, w, k: run_baseline_minco(q, T, w, obs_config=obs_cfg, dyn_config=dyn_cfg, bc_config=bc_cfg, reg_config=reg_cfg),
        "schemeA": lambda q, T, w, k: run_baseline_schemeA(q, T, w, obs_config=obs_cfg, dyn_config=dyn_cfg, bc_config=bc_cfg, reg_config=reg_cfg),
        "heuristic": lambda q, T, w, k: run_baseline_heuristic(q, T, w, obs_config=obs_cfg, dyn_config=dyn_cfg, bc_config=bc_cfg, reg_config=reg_cfg),
    }

    benchmark_rows: list[dict[str, Any]] = []
    for idx, M in enumerate(M_values):
        q_case, T_case = _sample_case(M, time_boxes[min(idx, len(time_boxes) - 1)], seed + idx)
        for method in benchmark_methods:
            result = baseline_runners[method](q_case, T_case, weights, 2)
            benchmark_rows.append(
                {
                    "method": method,
                    "M": int(M),
                    "k": 2,
                    "final_objective": float(result["final_objective"]),
                    "runtime_total": float(result["runtime_total"]),
                    "runtime_per_iter": float(result["runtime_per_iter"]),
                    "memory_peak": float(result["memory_peak"]) if result["memory_peak"] is not None else float("nan"),
                    "constraint_violation": float(result["constraint_violation"]) if result["constraint_violation"] is not None else float("nan"),
                    "support_width_mean": float(result["support_stats"]["support_width_mean"]),
                    "support_width_max": float(result["support_stats"]["support_width_max"]),
                    "history": result["history"],
                }
            )

    k_rows: list[dict[str, Any]] = []
    for k in k_values:
        try:
            result = baseline_runners["raw_schemeC"](q_rep, T_rep, weights, k)
            k_rows.append(
                {
                    "method": "raw_schemeC",
                    "M": int(T_rep.size),
                    "k": int(k),
                    "final_objective": float(result["final_objective"]),
                    "runtime_total": float(result["runtime_total"]),
                    "runtime_per_iter": float(result["runtime_per_iter"]),
                    "memory_peak": float(result["memory_peak"]) if result["memory_peak"] is not None else float("nan"),
                    "constraint_violation": float(result["constraint_violation"]) if result["constraint_violation"] is not None else float("nan"),
                    "support_width_mean": float(result["support_stats"]["support_width_mean"]),
                    "support_width_max": float(result["support_stats"]["support_width_max"]),
                }
            )
        except NotImplementedError:
            k_rows.append(
                {
                    "method": "raw_schemeC",
                    "M": int(T_rep.size),
                    "k": int(k),
                    "final_objective": float("nan"),
                    "runtime_total": float("nan"),
                    "runtime_per_iter": float("nan"),
                    "memory_peak": float("nan"),
                    "constraint_violation": float("nan"),
                    "support_width_mean": float("nan"),
                    "support_width_max": float("nan"),
                }
            )

    _compute_objective_gaps(benchmark_rows)
    _compute_objective_gaps(k_rows)

    baseline_compare_rows = [row for row in benchmark_rows if row["M"] == benchmark_rows[0]["M"]]
    for row in baseline_compare_rows:
        row.setdefault("objective_gap", 0.0)

    ablation_rows: list[dict[str, Any]] = []
    for weight_name in ["lambda_obs", "lambda_dyn", "lambda_T", "lambda_reg"]:
        ablated = dict(weights)
        ablated[weight_name] = 0.0
        result = run_baseline_raw_schemeC(q_rep, T_rep, ablated, k=2, n_steps=n_steps, T_min=T_min, obs_config=obs_cfg, dyn_config=dyn_cfg, bc_config=bc_cfg, reg_config=reg_cfg)
        ablation_rows.append(
            {
                "ablation_name": f"drop_{weight_name}",
                "method": "raw_schemeC",
                "M": int(T_rep.size),
                "k": 2,
                "final_objective": float(result["final_objective"]),
                "runtime_total": float(result["runtime_total"]),
                "objective_gap": 0.0,
            }
        )
    ref_val = min(row["final_objective"] for row in ablation_rows)
    for row in ablation_rows:
        row["objective_gap"] = (row["final_objective"] - ref_val) / max(abs(ref_val), 1e-12)

    _plot_runtime_vs_M(benchmark_rows, results_dirs["benchmark_M"] / "phase10_runtime_vs_M.png")
    _plot_obj_gap_vs_k(k_rows, results_dirs["benchmark_k"] / "phase10_obj_gap_vs_k.png")
    _plot_memory_vs_k(k_rows, results_dirs["benchmark_k"] / "phase10_memory_vs_k.png")
    _plot_baseline_compare(baseline_compare_rows, results_dirs["baseline"] / "phase10_baseline_compare.png")
    _plot_ablation_heatmap(ablation_rows, results_dirs["ablation"] / "phase10_ablation_heatmap.png")
    _plot_quality_speed_pareto(baseline_compare_rows, results_dirs["baseline"] / "phase10_quality_speed_pareto.png")
    _plot_sparsity_vs_k(k_rows, results_dirs["benchmark_k"] / "phase10_sparsity_vs_k.png")

    overview_rows = benchmark_rows + k_rows
    save_csv(
        results_dirs["summary"] / "phase10_overview.csv",
        ["method", "M", "k", "final_objective", "objective_gap", "runtime_total", "runtime_per_iter", "memory_peak", "constraint_violation", "support_width_mean", "support_width_max"],
        [
            [
                row["method"],
                row["M"],
                row["k"],
                row["final_objective"],
                row.get("objective_gap", float("nan")),
                row["runtime_total"],
                row["runtime_per_iter"],
                row["memory_peak"],
                row["constraint_violation"],
                row["support_width_mean"],
                row["support_width_max"],
            ]
            for row in overview_rows
        ],
    )
    save_csv(
        results_dirs["summary"] / "phase10_benchmark_summary.csv",
        ["method", "M", "k", "final_objective", "objective_gap", "runtime_total", "runtime_per_iter", "memory_peak", "constraint_violation", "support_width_mean", "support_width_max"],
        [
            [
                row["method"],
                row["M"],
                row["k"],
                row["final_objective"],
                row.get("objective_gap", float("nan")),
                row["runtime_total"],
                row["runtime_per_iter"],
                row["memory_peak"],
                row["constraint_violation"],
                row["support_width_mean"],
                row["support_width_max"],
            ]
            for row in benchmark_rows
        ],
    )
    save_csv(
        results_dirs["summary"] / "phase10_ablation_summary.csv",
        ["ablation_name", "method", "M", "k", "final_objective", "runtime_total", "objective_gap"],
        [[row["ablation_name"], row["method"], row["M"], row["k"], row["final_objective"], row["runtime_total"], row["objective_gap"]] for row in ablation_rows],
    )
    summary = {
        "benchmark_rows": overview_rows,
        "ablation_rows": ablation_rows,
        "config": {
            "M_values": M_values,
            "k_values": k_values,
            "obs_weight_values": obs_weight_values,
            "time_boxes": [list(box) for box in time_boxes],
            "benchmark_methods": benchmark_methods,
        },
    }
    save_json(results_dirs["summary"] / "phase10_benchmark_summary.json", summary)
    Path(results_dirs["summary"] / "phase10_benchmark_summary.md").write_text(
        "\n".join(
            [
                "# Phase 10 Benchmark Summary",
                "",
                f"- Raw Scheme C objective at k=2 available: `{np.isfinite(k_rows[0]['final_objective'])}`.",
                f"- Any general-k rows not yet implemented: `{any(not np.isfinite(row['final_objective']) for row in k_rows if row['k'] != 2)}`.",
                f"- Fastest baseline on representative case: `{min(baseline_compare_rows, key=lambda row: row['runtime_total'])['method']}`.",
                f"- Best objective on representative case: `{min(baseline_compare_rows, key=lambda row: row['final_objective'])['method']}`.",
                "",
                "Interpretation:",
                "- This benchmark layer is designed to compare quality, runtime, and sparsity under one common interface.",
                "- When k > 2 is not yet fully implemented, the framework still records this explicitly instead of silently fabricating comparison rows.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return {"overview_rows": overview_rows, "benchmark_rows": benchmark_rows, "ablation_rows": ablation_rows}


def main() -> None:
    q0, T0 = representative_case()
    run_benchmark_suite(q0=q0, T0=T0, weights=default_weights())


if __name__ == "__main__":
    main()
