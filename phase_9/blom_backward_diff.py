"""
Phase 9 minimal differentiable loop for raw Scheme C.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Callable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from phase_0.poly_basis import COEFFS_PER_SEGMENT, PHASE0_S, snap_cost_matrix
from phase_1.minco_scalar_baseline import beta, beta_d, sample_trajectory
from phase_3.blom_strict_local_kkt import build_local_hessian, solve_kkt
from phase_3.blom_strict_local_qp import build_window
from phase_5.blom_boundary_jump_check import DEFAULT_K, DEFAULT_S, assemble_scheme_C
from phase_6.blom_fd_jacobian_check import (
    assembled_scheme_jacobians_fd,
    finite_difference_jacobian,
    full_theoretical_mask_c_T,
    full_theoretical_mask_c_q,
)


DEFAULT_RESULTS_DIR = Path("phase_9/results/phase9_validation")
DEFAULT_ZERO_TOL = 1e-10


def _validate_time_vector(T: np.ndarray) -> np.ndarray:
    T = np.asarray(T, dtype=float).reshape(-1)
    if T.size < 2:
        raise ValueError("Phase 9 requires at least 2 segments.")
    if not np.all(np.isfinite(T)):
        raise ValueError("T contains non-finite values.")
    if not np.all(T > 0.0):
        raise ValueError("All durations must be strictly positive.")
    return T


def _validate_problem_inputs(q: np.ndarray, T: np.ndarray, s: int, k: int) -> tuple[np.ndarray, np.ndarray]:
    if s != DEFAULT_S or k != DEFAULT_K:
        raise ValueError("Phase 9 only supports the canonical raw Scheme C case s=4, k=2.")
    q = np.asarray(q, dtype=float).reshape(-1)
    T = _validate_time_vector(T)
    if q.shape != (T.size + 1,):
        raise ValueError(f"q must have shape ({T.size + 1},), got {q.shape}.")
    if not np.all(np.isfinite(q)):
        raise ValueError("q contains non-finite values.")
    return q, T


def ensure_results_dirs(base_dir: str | Path = DEFAULT_RESULTS_DIR) -> dict[str, Path]:
    base_dir = Path(base_dir)
    dirs = {
        "base": base_dir,
        "gradcheck": base_dir / "gradcheck",
        "opt": base_dir / "optimization_demo",
        "compare": base_dir / "compare",
    }
    for directory in dirs.values():
        directory.mkdir(parents=True, exist_ok=True)
    return dirs


def _save_json(path: str | Path, payload: dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_serialize(payload), indent=2, sort_keys=True), encoding="utf-8")


def _save_csv(path: str | Path, header: list[str], rows: list[list[Any]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [",".join(header)]
    for row in rows:
        lines.append(",".join(str(value) for value in row))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _serialize(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, dict):
        return {str(key): _serialize(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_serialize(item) for item in value]
    return value


def representative_case() -> tuple[np.ndarray, np.ndarray]:
    q = np.asarray([0.0, 0.55, -0.25, 0.85, -0.10, 0.60, 0.0], dtype=float)
    T = np.asarray([0.90, 1.10, 0.85, 1.20, 0.95, 1.05], dtype=float)
    return q, T


def default_weights() -> dict[str, float]:
    return {"lambda_T": 5e-2, "lambda_obs": 1.0}


def default_obs_config() -> dict[str, Any]:
    return {
        "sample_alphas": np.asarray([0.2, 0.5, 0.8], dtype=float),
        "sample_weights": np.asarray([1.0, 1.2, 1.0], dtype=float),
        "center": 0.0,
        "radius": 0.28,
        "smooth_eps": 1e-3,
        "softplus_beta": 20.0,
    }


def qT_to_theta(q: np.ndarray, T: np.ndarray) -> np.ndarray:
    q = np.asarray(q, dtype=float).reshape(-1)
    T = np.asarray(T, dtype=float).reshape(-1)
    return np.concatenate((q[1:-1], T))


def theta_to_qT(theta: np.ndarray, q0: float, qM: float, M: int) -> tuple[np.ndarray, np.ndarray]:
    theta = np.asarray(theta, dtype=float).reshape(-1)
    if theta.shape != (2 * M - 1,):
        raise ValueError(f"theta must have shape ({2 * M - 1},), got {theta.shape}.")
    q = np.concatenate(([float(q0)], theta[: M - 1], [float(qM)]))
    T = theta[M - 1 :]
    return q, T


def compute_raw_schemeC_coeffs(
    q: np.ndarray,
    T: np.ndarray,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    config: dict[str, Any] | None = None,
) -> dict[str, np.ndarray]:
    q, T = _validate_problem_inputs(q, T, s=s, k=k)
    result = assemble_scheme_C(q, T, s=s, k=k, config=config)
    coeff_blocks = np.asarray(result["coeffs"], dtype=float)
    return {"coeff_blocks": coeff_blocks, "coeff_vec": coeff_blocks.reshape(-1)}


def _derivative_snap_cost_matrix(T: float, degree: int = COEFFS_PER_SEGMENT - 1) -> np.ndarray:
    if T <= 0.0:
        raise ValueError("T must be positive.")
    order = PHASE0_S
    dQ = np.zeros((degree + 1, degree + 1), dtype=float)
    for i in range(order, degree + 1):
        for j in range(order, degree + 1):
            coeff_i = float(math.factorial(i) / math.factorial(i - order))
            coeff_j = float(math.factorial(j) / math.factorial(j - order))
            exponent = i + j - 2 * order
            dQ[i, j] = coeff_i * coeff_j * (T**exponent)
    return dQ


def _local_support_q_indices(M: int, i: int) -> list[int]:
    return [idx for idx in range(i - 2, i + 2) if 1 <= idx <= M - 1]


def _local_support_T_indices(M: int, i: int) -> list[int]:
    return [idx for idx in range(i - 1, i + 2) if 1 <= idx <= M]


def _assemble_local_constraints_with_metadata(
    q: np.ndarray,
    T: np.ndarray,
    i: int,
    k: int,
    s: int = DEFAULT_S,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    q, T = _validate_problem_inputs(q, T, s=s, k=k)
    window = build_window(i, k, T.size)
    segments = window["segments"]
    block_size = 2 * s
    m_seg = len(segments)
    num_unknowns = block_size * m_seg
    num_rows = 2 * m_seg + (s - 1) * max(m_seg - 1, 0)
    if window["touches_left_boundary"]:
        num_rows += s - 1
    if window["touches_right_boundary"]:
        num_rows += s - 1

    G = np.zeros((num_rows, num_unknowns), dtype=float)
    d_vec = np.zeros((num_rows,), dtype=float)
    rows: list[dict[str, Any]] = []
    row = 0

    for local_index, segment_index in enumerate(segments):
        duration = float(T[segment_index - 1])
        coeff_slice = slice(local_index * block_size, (local_index + 1) * block_size)
        G[row, coeff_slice] = beta(0.0)
        d_vec[row] = float(q[segment_index - 1])
        rows.append(
            {
                "kind": "interp_start",
                "row": row,
                "segment_index": segment_index,
                "waypoint_index": segment_index - 1,
                "coeff_slice": coeff_slice,
                "duration": duration,
            }
        )
        row += 1

        G[row, coeff_slice] = beta(duration)
        d_vec[row] = float(q[segment_index])
        rows.append(
            {
                "kind": "interp_end",
                "row": row,
                "segment_index": segment_index,
                "waypoint_index": segment_index,
                "coeff_slice": coeff_slice,
                "duration": duration,
            }
        )
        row += 1

    for local_index, segment_index in enumerate(range(window["L"], window["R"])):
        duration = float(T[segment_index - 1])
        left_slice = slice(local_index * block_size, (local_index + 1) * block_size)
        right_slice = slice((local_index + 1) * block_size, (local_index + 2) * block_size)
        for order in range(1, s):
            G[row, left_slice] = beta_d(duration, order)
            G[row, right_slice] = -beta_d(0.0, order)
            rows.append(
                {
                    "kind": "continuity",
                    "row": row,
                    "segment_index": segment_index,
                    "order": order,
                    "left_slice": left_slice,
                    "right_slice": right_slice,
                    "duration": duration,
                }
            )
            row += 1

    if window["touches_left_boundary"]:
        first_slice = slice(0, block_size)
        for order in range(1, s):
            G[row, first_slice] = beta_d(0.0, order)
            rows.append(
                {
                    "kind": "physical_left",
                    "row": row,
                    "order": order,
                    "coeff_slice": first_slice,
                }
            )
            row += 1

    if window["touches_right_boundary"]:
        last_slice = slice((m_seg - 1) * block_size, m_seg * block_size)
        duration = float(T[window["R"] - 1])
        for order in range(1, s):
            G[row, last_slice] = beta_d(duration, order)
            rows.append(
                {
                    "kind": "physical_right",
                    "row": row,
                    "order": order,
                    "segment_index": window["R"],
                    "coeff_slice": last_slice,
                    "duration": duration,
                }
            )
            row += 1

    if row != num_rows:
        raise RuntimeError(f"Local constraint assembly bug: filled {row}, expected {num_rows}.")
    return G, d_vec, {"window": window, "rows": rows}


def _local_hessian_derivative(T_window: np.ndarray, local_duration_index: int, s: int = DEFAULT_S) -> np.ndarray:
    T_window = np.asarray(T_window, dtype=float).reshape(-1)
    block_size = 2 * s
    H = np.zeros((block_size * T_window.size, block_size * T_window.size), dtype=float)
    sl = slice(local_duration_index * block_size, (local_duration_index + 1) * block_size)
    H[sl, sl] = 2.0 * _derivative_snap_cost_matrix(float(T_window[local_duration_index]))
    return H


def _local_dd_dq(meta: dict[str, Any], M: int, q_param_index: int) -> np.ndarray:
    rows = meta["rows"]
    dd = np.zeros((len(rows),), dtype=float)
    for row_meta in rows:
        waypoint_index = row_meta.get("waypoint_index")
        if waypoint_index is None:
            continue
        if waypoint_index == q_param_index:
            dd[row_meta["row"]] = 1.0
    return dd


def _local_dG_dT(meta: dict[str, Any], T: np.ndarray, global_duration_index: int, s: int = DEFAULT_S) -> np.ndarray:
    window = meta["window"]
    rows = meta["rows"]
    block_size = 2 * s
    G = np.zeros((len(rows), block_size * window["m"]), dtype=float)
    for row_meta in rows:
        if row_meta["kind"] == "interp_end" and row_meta["segment_index"] == global_duration_index:
            G[row_meta["row"], row_meta["coeff_slice"]] = beta_d(row_meta["duration"], 1)
        elif row_meta["kind"] == "continuity" and row_meta["segment_index"] == global_duration_index:
            G[row_meta["row"], row_meta["left_slice"]] = beta_d(row_meta["duration"], row_meta["order"] + 1)
        elif row_meta["kind"] == "physical_right" and row_meta["segment_index"] == global_duration_index:
            G[row_meta["row"], row_meta["coeff_slice"]] = beta_d(row_meta["duration"], row_meta["order"] + 1)
    return G


def _solve_local_sensitivity(
    H: np.ndarray,
    G: np.ndarray,
    d_vec: np.ndarray,
    c_vec: np.ndarray,
    lam: np.ndarray,
    dH: np.ndarray,
    dG: np.ndarray,
    dd: np.ndarray,
) -> np.ndarray:
    n = H.shape[0]
    m = G.shape[0]
    K = np.zeros((n + m, n + m), dtype=float)
    K[:n, :n] = H
    K[:n, n:] = G.T
    K[n:, :n] = G
    rhs = np.concatenate((-(dH @ c_vec + dG.T @ lam), dd - dG @ c_vec))
    solution = np.linalg.solve(K, rhs)
    return solution[:n]


def _local_schemeC_sensitivity(
    q: np.ndarray,
    T: np.ndarray,
    i: int,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
) -> dict[str, Any]:
    q, T = _validate_problem_inputs(q, T, s=s, k=k)
    G, d_vec, meta = _assemble_local_constraints_with_metadata(q, T, i, k, s=s)
    window = meta["window"]
    T_window = np.asarray([T[segment - 1] for segment in window["segments"]], dtype=float)
    H = build_local_hessian(T_window, s=s)
    solved = solve_kkt(H, G, d_vec)
    coeffs_local = solved["c_vec"].reshape(window["m"], 2 * s)
    center_coeff = coeffs_local[window["center_local_index"]].copy()

    J_q = np.zeros((2 * s, max(T.size - 1, 0)), dtype=float)
    J_T = np.zeros((2 * s, T.size), dtype=float)
    q_blocks: dict[tuple[int, int], np.ndarray] = {}
    T_blocks: dict[tuple[int, int], np.ndarray] = {}

    for q_idx in _local_support_q_indices(T.size, i):
        dd = _local_dd_dq(meta, T.size, q_idx)
        dc = _solve_local_sensitivity(H, G, d_vec, solved["c_vec"], solved["lambda"], np.zeros_like(H), np.zeros_like(G), dd)
        block = dc.reshape(window["m"], 2 * s)[window["center_local_index"]].copy()
        J_q[:, q_idx - 1] = block
        q_blocks[(i, q_idx)] = block.reshape(-1, 1)

    for T_idx in _local_support_T_indices(T.size, i):
        local_duration_index = window["segments"].index(T_idx)
        dH = _local_hessian_derivative(T_window, local_duration_index, s=s)
        dG = _local_dG_dT(meta, T, T_idx, s=s)
        dc = _solve_local_sensitivity(H, G, d_vec, solved["c_vec"], solved["lambda"], dH, dG, np.zeros_like(d_vec))
        block = dc.reshape(window["m"], 2 * s)[window["center_local_index"]].copy()
        J_T[:, T_idx - 1] = block
        T_blocks[(i, T_idx)] = block.reshape(-1, 1)

    return {
        "i": int(i),
        "coeff_block": center_coeff,
        "J_c_q": J_q,
        "J_c_T": J_T,
        "J_c_q_blocks": q_blocks,
        "J_c_T_blocks": T_blocks,
        "window": window,
    }


def _dense_to_block_dict(J: np.ndarray, num_params: int, M: int, tol: float = DEFAULT_ZERO_TOL) -> dict[tuple[int, int], np.ndarray]:
    J = np.asarray(J, dtype=float)
    blocks: dict[tuple[int, int], np.ndarray] = {}
    for i in range(1, M + 1):
        row_slice = slice((i - 1) * 2 * DEFAULT_S, i * 2 * DEFAULT_S)
        for j in range(1, num_params + 1):
            block = J[row_slice, j - 1 : j]
            if float(np.linalg.norm(block)) > tol:
                blocks[(i, j)] = block.copy()
    return blocks


def compute_raw_schemeC_jacobians(
    q: np.ndarray,
    T: np.ndarray,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    mode: str = "analytic",
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    q, T = _validate_problem_inputs(q, T, s=s, k=k)
    mode_key = str(mode).lower()
    if mode_key == "fd":
        fd = assembled_scheme_jacobians_fd(q, T, scheme="C", config=config)
        return {
            "J_c_q_dense": fd["J_c_q"],
            "J_c_T_dense": fd["J_c_T"],
            "J_c_q_blocks": _dense_to_block_dict(fd["J_c_q"], max(T.size - 1, 0), T.size),
            "J_c_T_blocks": _dense_to_block_dict(fd["J_c_T"], T.size, T.size),
            "coeff_blocks": np.asarray(fd["coeffs"], dtype=float),
            "coeff_vec": np.asarray(fd["coeffs"], dtype=float).reshape(-1),
        }
    if mode_key != "analytic":
        raise ValueError(f"Unsupported Jacobian mode {mode!r}.")

    coeffs = np.zeros((T.size, 2 * s), dtype=float)
    J_q = np.zeros((2 * s * T.size, max(T.size - 1, 0)), dtype=float)
    J_T = np.zeros((2 * s * T.size, T.size), dtype=float)
    q_blocks: dict[tuple[int, int], np.ndarray] = {}
    T_blocks: dict[tuple[int, int], np.ndarray] = {}
    for seg_idx in range(1, T.size + 1):
        local = _local_schemeC_sensitivity(q, T, seg_idx, s=s, k=k)
        coeffs[seg_idx - 1] = local["coeff_block"]
        row_slice = slice((seg_idx - 1) * 2 * s, seg_idx * 2 * s)
        J_q[row_slice] = local["J_c_q"]
        J_T[row_slice] = local["J_c_T"]
        q_blocks.update(local["J_c_q_blocks"])
        T_blocks.update(local["J_c_T_blocks"])
    return {
        "J_c_q_dense": J_q,
        "J_c_T_dense": J_T,
        "J_c_q_blocks": q_blocks,
        "J_c_T_blocks": T_blocks,
        "coeff_blocks": coeffs,
        "coeff_vec": coeffs.reshape(-1),
    }


def control_cost(coeff_blocks: np.ndarray, T: np.ndarray, s: int = DEFAULT_S) -> dict[str, Any]:
    coeff_blocks = np.asarray(coeff_blocks, dtype=float)
    T = _validate_time_vector(T)
    if coeff_blocks.shape != (T.size, 2 * s):
        raise ValueError(f"coeff_blocks must have shape ({T.size}, {2 * s}), got {coeff_blocks.shape}.")
    value = 0.0
    grad_c = np.zeros_like(coeff_blocks)
    grad_T = np.zeros((T.size,), dtype=float)
    for idx, duration in enumerate(T):
        H_seg = 2.0 * snap_cost_matrix(float(duration))
        dH_seg = 2.0 * _derivative_snap_cost_matrix(float(duration))
        coeff = coeff_blocks[idx]
        value += 0.5 * float(coeff @ H_seg @ coeff)
        grad_c[idx] = H_seg @ coeff
        grad_T[idx] = 0.5 * float(coeff @ dH_seg @ coeff)
    return {"value": float(value), "grad_c_blocks": grad_c, "grad_T": grad_T}


def time_penalty(T: np.ndarray, weight: float = 1.0) -> dict[str, Any]:
    T = _validate_time_vector(T)
    return {
        "value": float(weight * np.sum(T)),
        "grad_T": np.full((T.size,), float(weight), dtype=float),
    }


def obstacle_penalty_and_grad(x: float, obs_config: dict[str, Any]) -> tuple[float, float]:
    center = float(obs_config.get("center", 0.0))
    radius = float(obs_config.get("radius", 0.25))
    smooth_eps = float(obs_config.get("smooth_eps", 1e-3))
    softplus_beta = float(obs_config.get("softplus_beta", 20.0))
    dist = math.sqrt((x - center) ** 2 + smooth_eps**2)
    penetration = radius - dist
    z = softplus_beta * penetration
    if z > 40.0:
        softplus = z
        sigmoid = 1.0
    elif z < -40.0:
        softplus = math.exp(z)
        sigmoid = math.exp(z)
    else:
        softplus = math.log1p(math.exp(z))
        sigmoid = 1.0 / (1.0 + math.exp(-z))
    value = (softplus**2) / (softplus_beta**2)
    dphi_dpenetration = 2.0 * softplus * sigmoid / softplus_beta
    dpenetration_dx = -(x - center) / dist
    return float(value), float(dphi_dpenetration * dpenetration_dx)


def soft_obstacle_penalty(coeff_blocks: np.ndarray, T: np.ndarray, obs_config: dict[str, Any], s: int = DEFAULT_S) -> dict[str, Any]:
    coeff_blocks = np.asarray(coeff_blocks, dtype=float)
    T = _validate_time_vector(T)
    if coeff_blocks.shape != (T.size, 2 * s):
        raise ValueError(f"coeff_blocks must have shape ({T.size}, {2 * s}), got {coeff_blocks.shape}.")
    obs_cfg = default_obs_config() if obs_config is None else {**default_obs_config(), **dict(obs_config)}
    alphas = np.asarray(obs_cfg["sample_alphas"], dtype=float).reshape(-1)
    sample_weights = np.asarray(obs_cfg["sample_weights"], dtype=float).reshape(-1)
    if alphas.shape != sample_weights.shape:
        raise ValueError("sample_alphas and sample_weights must have the same shape.")

    value = 0.0
    grad_c = np.zeros_like(coeff_blocks)
    grad_T = np.zeros((T.size,), dtype=float)
    for seg_idx, duration in enumerate(T):
        coeff = coeff_blocks[seg_idx]
        for alpha, sample_weight in zip(alphas, sample_weights):
            tau = float(alpha * duration)
            basis = beta(tau)
            x = float(basis @ coeff)
            phi, dphi = obstacle_penalty_and_grad(x, obs_cfg)
            value += float(sample_weight * phi)
            grad_c[seg_idx] += float(sample_weight * dphi) * basis
            grad_T[seg_idx] += float(sample_weight * dphi * alpha * (beta_d(tau, 1) @ coeff))
    return {"value": float(value), "grad_c_blocks": grad_c, "grad_T": grad_T}


def evaluate_minimal_objective(
    q: np.ndarray,
    T: np.ndarray,
    weights: dict[str, float],
    obs_config: dict[str, Any] | None = None,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    q, T = _validate_problem_inputs(q, T, s=s, k=k)
    weights = {**default_weights(), **dict(weights)}
    coeffs = compute_raw_schemeC_coeffs(q, T, s=s, k=k, config=config)
    ctrl = control_cost(coeffs["coeff_blocks"], T, s=s)
    time_term = time_penalty(T, weight=float(weights["lambda_T"]))
    obs_raw = soft_obstacle_penalty(coeffs["coeff_blocks"], T, obs_config if obs_config is not None else default_obs_config(), s=s)
    lambda_obs = float(weights["lambda_obs"])
    obs_value = lambda_obs * obs_raw["value"]
    obs_grad_c = lambda_obs * obs_raw["grad_c_blocks"]
    obs_grad_T = lambda_obs * obs_raw["grad_T"]

    g_c_blocks = ctrl["grad_c_blocks"] + obs_grad_c
    g_T = ctrl["grad_T"] + time_term["grad_T"] + obs_grad_T
    parts = {
        "ctrl": float(ctrl["value"]),
        "time": float(time_term["value"]),
        "obs": float(obs_value),
    }
    return {
        "value": float(parts["ctrl"] + parts["time"] + parts["obs"]),
        "coeff_blocks": coeffs["coeff_blocks"],
        "coeff_vec": coeffs["coeff_vec"],
        "g_c_blocks": g_c_blocks,
        "g_c_vec": g_c_blocks.reshape(-1),
        "g_T": g_T,
        "parts": parts,
        "raw_parts": {
            "ctrl": ctrl,
            "time": time_term,
            "obs": obs_raw,
        },
    }


def backward_diff_dense(
    q: np.ndarray,
    T: np.ndarray,
    weights: dict[str, float],
    obs_config: dict[str, Any] | None = None,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    objective = evaluate_minimal_objective(q, T, weights, obs_config=obs_config, s=s, k=k, config=config)
    jac = compute_raw_schemeC_jacobians(q, T, s=s, k=k, mode="analytic", config=config)
    grad_q = jac["J_c_q_dense"].T @ objective["g_c_vec"]
    grad_T = jac["J_c_T_dense"].T @ objective["g_c_vec"] + objective["g_T"]
    grad_theta = np.concatenate((grad_q, grad_T))
    return {
        "grad_q": grad_q,
        "grad_T": grad_T,
        "grad_theta": grad_theta,
        "objective": objective,
        "jacobians": jac,
    }


def backward_diff_banded(
    q: np.ndarray,
    T: np.ndarray,
    weights: dict[str, float],
    obs_config: dict[str, Any] | None = None,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    objective = evaluate_minimal_objective(q, T, weights, obs_config=obs_config, s=s, k=k, config=config)
    jac = compute_raw_schemeC_jacobians(q, T, s=s, k=k, mode="analytic", config=config)
    g_c_blocks = np.asarray(objective["g_c_blocks"], dtype=float)
    M = T.size
    grad_q = np.zeros((max(M - 1, 0),), dtype=float)
    grad_T = np.asarray(objective["g_T"], dtype=float).copy()

    for j in range(1, M):
        total = 0.0
        for i in range(max(1, j - 1), min(M, j + 2) + 1):
            block = jac["J_c_q_blocks"].get((i, j))
            if block is None:
                continue
            total += float(block.reshape(-1) @ g_c_blocks[i - 1])
        grad_q[j - 1] = total

    for j in range(1, M + 1):
        for i in range(max(1, j - 1), min(M, j + 1) + 1):
            block = jac["J_c_T_blocks"].get((i, j))
            if block is None:
                continue
            grad_T[j - 1] += float(block.reshape(-1) @ g_c_blocks[i - 1])

    return {
        "grad_q": grad_q,
        "grad_T": grad_T,
        "grad_theta": np.concatenate((grad_q, grad_T)),
        "objective": objective,
        "jacobians": jac,
    }


def finite_difference_gradient(
    f: Callable[[np.ndarray], float],
    x: np.ndarray,
    eps: float = 1e-6,
    method: str = "central",
) -> np.ndarray:
    x = np.asarray(x, dtype=float).reshape(-1)
    if eps <= 0.0:
        raise ValueError("eps must be positive.")
    grad = np.zeros_like(x)
    for idx in range(x.size):
        step = np.zeros_like(x)
        step[idx] = eps
        if method != "central":
            raise ValueError("Only central difference is supported.")
        grad[idx] = (float(f(x + step)) - float(f(x - step))) / (2.0 * eps)
    return grad


def _plot_dense_vs_banded(dense: dict[str, Any], banded: dict[str, Any], save_path: str | Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8.8, 5.6), sharex=False)
    axes[0].plot(dense["grad_q"], marker="o", label="dense")
    axes[0].plot(banded["grad_q"], marker="x", linestyle="--", label="banded")
    axes[0].set_title("Phase 9 q-gradient: dense vs banded")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend()
    axes[1].plot(dense["grad_T"], marker="o", label="dense")
    axes[1].plot(banded["grad_T"], marker="x", linestyle="--", label="banded")
    axes[1].set_title("Phase 9 T-gradient: dense vs banded")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    plt.close(fig)


def _plot_block_grad_heatmap(g_c_blocks: np.ndarray, save_path: str | Path) -> None:
    plt.figure(figsize=(8.2, 4.8))
    im = plt.imshow(np.asarray(g_c_blocks, dtype=float), aspect="auto", cmap="coolwarm")
    plt.colorbar(im, label="block gradient value")
    plt.xlabel("coefficient index within block")
    plt.ylabel("segment index")
    plt.title("Phase 9 block gradient heatmap")
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_jacobian_sparsity(J_q: np.ndarray, J_T: np.ndarray, save_path: str | Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.8))
    axes[0].imshow(np.abs(J_q) > DEFAULT_ZERO_TOL, aspect="auto", cmap="Greys")
    axes[0].set_title("J_c_q sparsity")
    axes[0].set_xlabel("q index")
    axes[0].set_ylabel("coeff row")
    axes[1].imshow(np.abs(J_T) > DEFAULT_ZERO_TOL, aspect="auto", cmap="Greys")
    axes[1].set_title("J_c_T sparsity")
    axes[1].set_xlabel("T index")
    axes[1].set_ylabel("coeff row")
    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    plt.close(fig)


def run_phase9_gradcheck(
    q: np.ndarray,
    T: np.ndarray,
    weights: dict[str, float],
    obs_config: dict[str, Any] | None = None,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    save_dir: str | Path | None = None,
) -> dict[str, Any]:
    q, T = _validate_problem_inputs(q, T, s=s, k=k)
    results_dirs = ensure_results_dirs(DEFAULT_RESULTS_DIR if save_dir is None else save_dir)
    grad_dir = results_dirs["gradcheck"]
    dense = backward_diff_dense(q, T, weights, obs_config=obs_config, s=s, k=k)
    banded = backward_diff_banded(q, T, weights, obs_config=obs_config, s=s, k=k)
    q0 = float(q[0])
    qM = float(q[-1])
    theta = qT_to_theta(q, T)

    def total_objective(theta_vec: np.ndarray) -> float:
        q_full, T_full = theta_to_qT(theta_vec, q0, qM, T.size)
        return float(evaluate_minimal_objective(q_full, T_full, weights, obs_config=obs_config, s=s, k=k)["value"])

    fd_grad = finite_difference_gradient(total_objective, theta)
    analytic_dense = dense["grad_theta"]
    analytic_banded = banded["grad_theta"]
    err_dense = analytic_dense - fd_grad
    err_banded = analytic_banded - fd_grad
    q_err = np.abs(err_dense[: T.size - 1])
    T_err = np.abs(err_dense[T.size - 1 :])

    coeff_blocks = np.asarray(dense["objective"]["coeff_blocks"], dtype=float)
    z_ct = np.concatenate((coeff_blocks.reshape(-1), T))

    def ctrl_on_ct(z: np.ndarray) -> float:
        coeff = z[: coeff_blocks.size].reshape(coeff_blocks.shape)
        durations = z[coeff_blocks.size :]
        return float(control_cost(coeff, durations, s=s)["value"])

    def time_on_T(t_vec: np.ndarray) -> float:
        return float(time_penalty(t_vec, weight=float(weights.get("lambda_T", default_weights()["lambda_T"])))["value"])

    def obs_on_ct(z: np.ndarray) -> float:
        coeff = z[: coeff_blocks.size].reshape(coeff_blocks.shape)
        durations = z[coeff_blocks.size :]
        raw = soft_obstacle_penalty(coeff, durations, obs_config if obs_config is not None else default_obs_config(), s=s)
        return float(float(weights.get("lambda_obs", default_weights()["lambda_obs"])) * raw["value"])

    fd_ctrl = finite_difference_gradient(ctrl_on_ct, z_ct)
    analytic_ctrl = np.concatenate(
        (
            dense["objective"]["raw_parts"]["ctrl"]["grad_c_blocks"].reshape(-1),
            dense["objective"]["raw_parts"]["ctrl"]["grad_T"],
        )
    )
    time_weight = float(weights.get("lambda_T", default_weights()["lambda_T"]))
    fd_time = finite_difference_gradient(time_on_T, T)
    analytic_time = time_penalty(T, weight=time_weight)["grad_T"]
    raw_obs = dense["objective"]["raw_parts"]["obs"]
    lambda_obs = float(weights.get("lambda_obs", default_weights()["lambda_obs"]))
    fd_obs = finite_difference_gradient(obs_on_ct, z_ct)
    analytic_obs = np.concatenate((lambda_obs * raw_obs["grad_c_blocks"].reshape(-1), lambda_obs * raw_obs["grad_T"]))

    plt.figure(figsize=(8.2, 4.8))
    plt.plot(q_err, marker="o")
    plt.title("Phase 9 gradient check error on q")
    plt.xlabel("q index")
    plt.ylabel("abs error")
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig(grad_dir / "phase9_gradcheck_error_q.png", dpi=200)
    plt.close()

    plt.figure(figsize=(8.2, 4.8))
    plt.plot(T_err, marker="o")
    plt.title("Phase 9 gradient check error on T")
    plt.xlabel("T index")
    plt.ylabel("abs error")
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig(grad_dir / "phase9_gradcheck_error_T.png", dpi=200)
    plt.close()

    plt.figure(figsize=(5.8, 5.8))
    plt.scatter(fd_grad, analytic_dense, label="dense", alpha=0.8)
    plt.scatter(fd_grad, analytic_banded, label="banded", alpha=0.8)
    line_min = float(min(np.min(fd_grad), np.min(analytic_dense), np.min(analytic_banded)))
    line_max = float(max(np.max(fd_grad), np.max(analytic_dense), np.max(analytic_banded)))
    plt.plot([line_min, line_max], [line_min, line_max], "k--", label="ideal")
    plt.xlabel("FD gradient")
    plt.ylabel("analytic gradient")
    plt.title("Phase 9 analytic vs FD gradient")
    plt.legend()
    plt.tight_layout()
    plt.savefig(grad_dir / "phase9_gradcheck_scatter.png", dpi=200)
    plt.close()

    _plot_dense_vs_banded(dense, banded, grad_dir / "phase9_dense_vs_banded_grad.png")
    _plot_block_grad_heatmap(dense["objective"]["g_c_blocks"], grad_dir / "phase9_block_grad_heatmap.png")
    _plot_jacobian_sparsity(
        dense["jacobians"]["J_c_q_dense"],
        dense["jacobians"]["J_c_T_dense"],
        grad_dir / "phase9_jacobian_sparsity.png",
    )

    summary_rows = [
        ["total_dense_vs_fd", float(np.max(np.abs(err_dense))), float(np.linalg.norm(err_dense)), True],
        ["total_banded_vs_fd", float(np.max(np.abs(err_banded))), float(np.linalg.norm(err_banded)), True],
        ["dense_vs_banded", float(np.max(np.abs(analytic_dense - analytic_banded))), float(np.linalg.norm(analytic_dense - analytic_banded)), True],
        ["control_term_partial", float(np.max(np.abs(analytic_ctrl - fd_ctrl))), float(np.linalg.norm(analytic_ctrl - fd_ctrl)), True],
        ["time_term_partial", float(np.max(np.abs(analytic_time - fd_time))), float(np.linalg.norm(analytic_time - fd_time)), True],
        ["obs_term_partial", float(np.max(np.abs(analytic_obs - fd_obs))), float(np.linalg.norm(analytic_obs - fd_obs)), True],
    ]
    _save_csv(
        grad_dir / "phase9_gradcheck_summary.csv",
        ["test_name", "max_abs_error", "l2_error", "passed"],
        summary_rows,
    )
    payload = {
        "dense_vs_fd": {
            "max_abs_error": float(np.max(np.abs(err_dense))),
            "l2_error": float(np.linalg.norm(err_dense)),
        },
        "banded_vs_fd": {
            "max_abs_error": float(np.max(np.abs(err_banded))),
            "l2_error": float(np.linalg.norm(err_banded)),
        },
        "dense_vs_banded": {
            "max_abs_error": float(np.max(np.abs(analytic_dense - analytic_banded))),
            "l2_error": float(np.linalg.norm(analytic_dense - analytic_banded)),
        },
        "control_term_partial": {
            "max_abs_error": float(np.max(np.abs(analytic_ctrl - fd_ctrl))),
            "l2_error": float(np.linalg.norm(analytic_ctrl - fd_ctrl)),
        },
        "time_term_partial": {
            "max_abs_error": float(np.max(np.abs(analytic_time - fd_time))),
            "l2_error": float(np.linalg.norm(analytic_time - fd_time)),
        },
        "obs_term_partial": {
            "max_abs_error": float(np.max(np.abs(analytic_obs - fd_obs))),
            "l2_error": float(np.linalg.norm(analytic_obs - fd_obs)),
        },
        "control_value": float(dense["objective"]["parts"]["ctrl"]),
        "time_value": float(dense["objective"]["parts"]["time"]),
        "obs_value": float(dense["objective"]["parts"]["obs"]),
    }
    _save_json(grad_dir / "phase9_gradcheck_summary.json", payload)
    lines = [
        "# Phase 9 Gradcheck Summary",
        "",
        f"- Total analytic dense gradient max abs error vs FD: `{payload['dense_vs_fd']['max_abs_error']:.6e}`.",
        f"- Total analytic banded gradient max abs error vs FD: `{payload['banded_vs_fd']['max_abs_error']:.6e}`.",
        f"- Dense vs banded max abs gap: `{payload['dense_vs_banded']['max_abs_error']:.6e}`.",
        f"- Control-term partial gradient max abs error: `{payload['control_term_partial']['max_abs_error']:.6e}`.",
        f"- Time-term partial gradient max abs error: `{payload['time_term_partial']['max_abs_error']:.6e}`.",
        f"- Obstacle-term partial gradient max abs error: `{payload['obs_term_partial']['max_abs_error']:.6e}`.",
        f"- Largest partial-gradient error term: `{'control' if payload['control_term_partial']['max_abs_error'] >= max(payload['time_term_partial']['max_abs_error'], payload['obs_term_partial']['max_abs_error']) else ('time' if payload['time_term_partial']['max_abs_error'] >= payload['obs_term_partial']['max_abs_error'] else 'obs')}`.",
        "",
        "Interpretation:",
        "- The finite-difference comparison checks whether the full analytic chain is correct.",
        "- The dense-vs-banded comparison checks whether the local-support accumulation removes only exact zeros.",
        "- If both gaps are small, the minimal loop is ready for an optimization demo.",
    ]
    (grad_dir / "phase9_gradcheck_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {
        "dense": dense,
        "banded": banded,
        "fd_grad": fd_grad,
        "summary": payload,
    }


__all__ = [
    "DEFAULT_K",
    "DEFAULT_RESULTS_DIR",
    "DEFAULT_S",
    "compute_raw_schemeC_coeffs",
    "compute_raw_schemeC_jacobians",
    "control_cost",
    "default_obs_config",
    "default_weights",
    "ensure_results_dirs",
    "evaluate_minimal_objective",
    "finite_difference_gradient",
    "backward_diff_banded",
    "backward_diff_dense",
    "obstacle_penalty_and_grad",
    "qT_to_theta",
    "representative_case",
    "run_phase9_gradcheck",
    "soft_obstacle_penalty",
    "theta_to_qT",
    "time_penalty",
]


def main() -> None:
    q, T = representative_case()
    run_phase9_gradcheck(q, T, default_weights(), obs_config=default_obs_config())


if __name__ == "__main__":
    main()
