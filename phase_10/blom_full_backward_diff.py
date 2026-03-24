"""
Phase 10 full block-banded backward differentiation framework.
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Any, Callable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from phase_0.poly_basis import COEFFS_PER_SEGMENT, PHASE0_S, snap_cost_matrix
from phase_1.minco_scalar_baseline import beta, beta_d, evaluate_segment, sample_trajectory
from phase_3.blom_strict_local_kkt import build_local_hessian, solve_kkt
from phase_3.blom_strict_local_qp import build_window
from phase_5.blom_boundary_jump_check import DEFAULT_K, DEFAULT_S, assemble_scheme_C


DEFAULT_RESULTS_DIR = Path("phase_10/results/phase10_framework")
DEFAULT_ZERO_TOL = 1e-10


def _validate_time_vector(T: np.ndarray) -> np.ndarray:
    T = np.asarray(T, dtype=float).reshape(-1)
    if T.size < 2:
        raise ValueError("Phase 10 requires at least 2 segments.")
    if not np.all(np.isfinite(T)):
        raise ValueError("T contains non-finite values.")
    if not np.all(T > 0.0):
        raise ValueError("All segment durations must be strictly positive.")
    return T


def _validate_problem_inputs(q: np.ndarray, T: np.ndarray, s: int = DEFAULT_S) -> tuple[np.ndarray, np.ndarray]:
    if s != DEFAULT_S:
        raise ValueError("Phase 10 currently supports only the canonical setting s=4.")
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
        "backward": base_dir / "backward_diff",
        "optimizer": base_dir / "optimizer",
        "benchmark_M": base_dir / "benchmark_M_sweep",
        "benchmark_k": base_dir / "benchmark_k_sweep",
        "baseline": base_dir / "baseline_compare",
        "ablation": base_dir / "ablation",
        "summary": base_dir / "summary",
    }
    for directory in dirs.values():
        directory.mkdir(parents=True, exist_ok=True)
    return dirs


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


def save_json(path: str | Path, payload: dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_serialize(payload), indent=2, sort_keys=True), encoding="utf-8")


def save_csv(path: str | Path, header: list[str], rows: list[list[Any]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [",".join(header)]
    for row in rows:
        lines.append(",".join(str(value) for value in row))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def representative_case() -> tuple[np.ndarray, np.ndarray]:
    q = np.asarray([0.0, 0.65, -0.30, 0.90, -0.15, 0.72, 0.0], dtype=float)
    T = np.asarray([0.95, 1.15, 0.85, 1.25, 1.00, 1.10], dtype=float)
    return q, T


def default_weights() -> dict[str, float]:
    return {
        "lambda_T": 5e-2,
        "lambda_obs": 1.0,
        "lambda_dyn": 1e-1,
        "lambda_bc": 5e-1,
        "lambda_reg": 1e-1,
    }


def default_obs_config() -> dict[str, Any]:
    return {
        "sample_alphas": np.asarray([0.2, 0.5, 0.8], dtype=float),
        "sample_weights": np.asarray([1.0, 1.2, 1.0], dtype=float),
        "derivative_orders": [0],
        "component_index": 0,
        "center": 0.0,
        "radius": 0.28,
        "smooth_eps": 1e-3,
        "softplus_beta": 20.0,
    }


def default_dyn_config() -> dict[str, Any]:
    return {
        "sample_alphas": np.asarray([0.25, 0.5, 0.75], dtype=float),
        "sample_weights": np.asarray([1.0, 1.0, 1.0], dtype=float),
        "derivative_orders": [1, 2, 3],
        "limits": np.asarray([1.6, 2.8, 7.5], dtype=float),
    }


def default_bc_config(q: np.ndarray | None = None) -> dict[str, Any]:
    q0 = 0.0 if q is None else float(np.asarray(q, dtype=float).reshape(-1)[0])
    qM = 0.0 if q is None else float(np.asarray(q, dtype=float).reshape(-1)[-1])
    return {
        "orders": [0, 1, 2, 3],
        "desired_start": np.asarray([q0, 0.0, 0.0, 0.0], dtype=float),
        "desired_end": np.asarray([qM, 0.0, 0.0, 0.0], dtype=float),
        "weight_start": 1.0,
        "weight_end": 1.0,
    }


def default_reg_config() -> dict[str, Any]:
    return {"mode": "time_smoothing"}


def full_q_to_qbar(q: np.ndarray) -> np.ndarray:
    q = np.asarray(q, dtype=float).reshape(-1)
    return q[1:-1].copy()


def qbar_to_full(q_bar: np.ndarray, q_start: float, q_end: float) -> np.ndarray:
    q_bar = np.asarray(q_bar, dtype=float).reshape(-1)
    return np.concatenate(([float(q_start)], q_bar, [float(q_end)]))


def q_tau_to_xi(q: np.ndarray, tau: np.ndarray) -> np.ndarray:
    q = np.asarray(q, dtype=float).reshape(-1)
    tau = np.asarray(tau, dtype=float).reshape(-1)
    if tau.size != q.size - 1:
        raise ValueError("tau must have one entry per segment.")
    return np.concatenate((q[1:-1], tau))


def xi_to_qbar_tau(xi: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    xi = np.asarray(xi, dtype=float).reshape(-1)
    if xi.size < 3 or xi.size % 2 == 0:
        raise ValueError("xi must have odd length 2*M-1 with M >= 2.")
    M = (xi.size + 1) // 2
    return xi[: M - 1].copy(), xi[M - 1 :].copy()


def tau_to_T(tau: np.ndarray, T_min: float = 1e-3) -> np.ndarray:
    tau = np.asarray(tau, dtype=float).reshape(-1)
    if T_min <= 0.0:
        raise ValueError("T_min must be positive.")
    softplus = np.log1p(np.exp(-np.abs(tau))) + np.maximum(tau, 0.0)
    return float(T_min) + softplus


def dT_dtau(tau: np.ndarray) -> np.ndarray:
    tau = np.asarray(tau, dtype=float).reshape(-1)
    sigma = np.empty_like(tau)
    positive = tau >= 0.0
    sigma[positive] = 1.0 / (1.0 + np.exp(-tau[positive]))
    exp_tau = np.exp(tau[~positive])
    sigma[~positive] = exp_tau / (1.0 + exp_tau)
    return sigma


def T_to_tau(T: np.ndarray, T_min: float = 1e-3) -> np.ndarray:
    T = np.asarray(T, dtype=float).reshape(-1)
    if np.any(T <= T_min):
        raise ValueError("T must satisfy T > T_min to invert softplus.")
    x = T - float(T_min)
    return x + np.log(-np.expm1(-x))


def _require_supported_k(k: int) -> None:
    if k < 0 or k % 2 != 0:
        raise ValueError("Phase 10 expects nonnegative even k.")


def compute_raw_schemeC_coeffs_general_k(
    q: np.ndarray,
    T: np.ndarray,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    config: dict[str, Any] | None = None,
) -> dict[str, np.ndarray]:
    _require_supported_k(k)
    q, T = _validate_problem_inputs(q, T, s=s)
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


def _assemble_local_constraints_with_metadata(
    q: np.ndarray,
    T: np.ndarray,
    i: int,
    k: int,
    s: int = DEFAULT_S,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    q, T = _validate_problem_inputs(q, T, s=s)
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
            rows.append({"kind": "physical_left", "row": row, "order": order, "coeff_slice": first_slice})
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
        raise RuntimeError(f"Constraint assembly bug: filled {row}, expected {num_rows}.")
    return G, d_vec, {"window": window, "rows": rows}


def _local_hessian_derivative(T_window: np.ndarray, local_duration_index: int, s: int = DEFAULT_S) -> np.ndarray:
    T_window = np.asarray(T_window, dtype=float).reshape(-1)
    block_size = 2 * s
    H = np.zeros((block_size * T_window.size, block_size * T_window.size), dtype=float)
    sl = slice(local_duration_index * block_size, (local_duration_index + 1) * block_size)
    H[sl, sl] = 2.0 * _derivative_snap_cost_matrix(float(T_window[local_duration_index]))
    return H


def _local_dd_dq(meta: dict[str, Any], q_param_index: int) -> np.ndarray:
    rows = meta["rows"]
    dd = np.zeros((len(rows),), dtype=float)
    for row_meta in rows:
        waypoint_index = row_meta.get("waypoint_index")
        if waypoint_index is None:
            continue
        if waypoint_index == q_param_index:
            dd[row_meta["row"]] = 1.0
    return dd


def _local_dG_dT(meta: dict[str, T], global_duration_index: int, s: int = DEFAULT_S) -> np.ndarray:
    rows = meta["rows"]
    window = meta["window"]
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


def _local_support_q_indices(meta: dict[str, Any], M: int) -> list[int]:
    indices = sorted(
        {
            int(row_meta["waypoint_index"])
            for row_meta in meta["rows"]
            if row_meta.get("waypoint_index") is not None and 1 <= int(row_meta["waypoint_index"]) <= M - 1
        }
    )
    return indices


def _local_support_T_indices(meta: dict[str, Any]) -> list[int]:
    return sorted(int(segment) for segment in meta["window"]["segments"])


def _local_schemeC_sensitivity(q: np.ndarray, T: np.ndarray, i: int, s: int = DEFAULT_S, k: int = DEFAULT_K) -> dict[str, Any]:
    q, T = _validate_problem_inputs(q, T, s=s)
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

    for q_idx in _local_support_q_indices(meta, T.size):
        dd = _local_dd_dq(meta, q_idx)
        dc = _solve_local_sensitivity(H, G, d_vec, solved["c_vec"], solved["lambda"], np.zeros_like(H), np.zeros_like(G), dd)
        block = dc.reshape(window["m"], 2 * s)[window["center_local_index"]].copy()
        J_q[:, q_idx - 1] = block
        q_blocks[(i, q_idx)] = block.reshape(-1, 1)

    for T_idx in _local_support_T_indices(meta):
        local_duration_index = window["segments"].index(T_idx)
        dH = _local_hessian_derivative(T_window, local_duration_index, s=s)
        dG = _local_dG_dT(meta, T_idx, s=s)
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


def compute_raw_schemeC_jacobians_general_k(
    q: np.ndarray,
    T: np.ndarray,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    mode: str = "analytic",
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    del config
    _require_supported_k(k)
    q, T = _validate_problem_inputs(q, T, s=s)
    mode_key = str(mode).lower()
    if mode_key != "analytic":
        raise NotImplementedError("Phase 10 currently provides the analytic checker path only.")

    coeffs = np.zeros((T.size, 2 * s), dtype=float)
    J_q = np.zeros((2 * s * T.size, max(T.size - 1, 0)), dtype=float)
    J_T = np.zeros((2 * s * T.size, T.size), dtype=float)
    q_blocks: dict[tuple[int, int], np.ndarray] = {}
    T_blocks: dict[tuple[int, int], np.ndarray] = {}
    support_q: dict[int, list[int]] = {}
    support_T: dict[int, list[int]] = {}
    for seg_idx in range(1, T.size + 1):
        local = _local_schemeC_sensitivity(q, T, seg_idx, s=s, k=k)
        coeffs[seg_idx - 1] = local["coeff_block"]
        row_slice = slice((seg_idx - 1) * 2 * s, seg_idx * 2 * s)
        J_q[row_slice] = local["J_c_q"]
        J_T[row_slice] = local["J_c_T"]
        q_blocks.update(local["J_c_q_blocks"])
        T_blocks.update(local["J_c_T_blocks"])
        support_q[seg_idx] = sorted(j for (i, j) in local["J_c_q_blocks"] if i == seg_idx)
        support_T[seg_idx] = sorted(j for (i, j) in local["J_c_T_blocks"] if i == seg_idx)
    return {
        "J_c_q_dense": J_q,
        "J_c_T_dense": J_T,
        "J_c_q_blocks": q_blocks,
        "J_c_T_blocks": T_blocks,
        "support_q": support_q,
        "support_T": support_T,
        "coeff_blocks": coeffs,
        "coeff_vec": coeffs.reshape(-1),
    }


def _stacked_operator(t: float, derivative_orders: list[int]) -> np.ndarray:
    return np.vstack([beta_d(float(t), int(order)) for order in derivative_orders])


def _soft_band_penalty(z: np.ndarray, T_val: float, config: dict[str, Any]) -> dict[str, Any]:
    del T_val
    center = float(config.get("center", 0.0))
    radius = float(config.get("radius", 0.25))
    smooth_eps = float(config.get("smooth_eps", 1e-3))
    beta_sp = float(config.get("softplus_beta", 20.0))
    component_index = int(config.get("component_index", 0))
    x = float(np.asarray(z, dtype=float).reshape(-1)[component_index])
    dist = math.sqrt((x - center) ** 2 + smooth_eps**2)
    penetration = radius - dist
    scaled = beta_sp * penetration
    if scaled > 40.0:
        softplus = scaled
        sigmoid = 1.0
    elif scaled < -40.0:
        softplus = math.exp(scaled)
        sigmoid = math.exp(scaled)
    else:
        softplus = math.log1p(math.exp(scaled))
        sigmoid = 1.0 / (1.0 + math.exp(-scaled))
    value = (softplus**2) / (beta_sp**2)
    dvalue_dpenetration = 2.0 * softplus * sigmoid / beta_sp
    dpenetration_dx = -(x - center) / dist
    grad_z = np.zeros_like(np.asarray(z, dtype=float).reshape(-1))
    grad_z[component_index] = dvalue_dpenetration * dpenetration_dx
    return {"value": float(value), "grad_z": grad_z, "grad_T": 0.0}


def _quadratic_excess_penalty(z: np.ndarray, T_val: float, config: dict[str, Any]) -> dict[str, Any]:
    del T_val
    values = np.asarray(z, dtype=float).reshape(-1)
    limits = np.asarray(config.get("limits", np.ones_like(values)), dtype=float).reshape(-1)
    if limits.size == 1:
        limits = np.full(values.shape, float(limits[0]), dtype=float)
    if limits.shape != values.shape:
        raise ValueError("limits must broadcast to the sampled jet shape.")
    excess = np.maximum(np.abs(values) - limits, 0.0)
    value = 0.5 * float(np.sum(excess**2))
    grad_z = np.sign(values) * excess
    return {"value": value, "grad_z": grad_z, "grad_T": 0.0}


def control_cost_full(coeff_blocks: np.ndarray, T: np.ndarray, s: int = DEFAULT_S) -> dict[str, Any]:
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


def time_penalty_full(T: np.ndarray, psi_type: str = "linear", psi_config: dict[str, Any] | None = None) -> dict[str, Any]:
    T = _validate_time_vector(T)
    psi_config = {} if psi_config is None else dict(psi_config)
    mode = str(psi_type).lower()
    if mode == "linear":
        return {"value": float(np.sum(T)), "grad_T": np.ones_like(T)}
    if mode == "quadratic":
        return {"value": float(np.sum(T**2)), "grad_T": 2.0 * T}
    if mode == "custom":
        fn = psi_config.get("callable")
        if not callable(fn):
            raise ValueError("custom psi_type requires psi_config['callable'].")
        value = 0.0
        grad = np.zeros_like(T)
        for idx, duration in enumerate(T):
            out = fn(float(duration), dict(psi_config))
            value += float(out["value"])
            grad[idx] = float(out["grad"])
        return {"value": float(value), "grad_T": grad}
    raise ValueError(f"Unsupported psi_type {psi_type!r}.")


def _sampled_penalty_common(
    coeff_blocks: np.ndarray,
    T: np.ndarray,
    config: dict[str, Any],
    *,
    default_config: dict[str, Any],
    penalty_builder: Callable[[np.ndarray, float, dict[str, Any]], dict[str, Any]],
    s: int = DEFAULT_S,
) -> dict[str, Any]:
    coeff_blocks = np.asarray(coeff_blocks, dtype=float)
    T = _validate_time_vector(T)
    cfg = {**default_config, **({} if config is None else dict(config))}
    alphas = np.asarray(cfg.get("sample_alphas", default_config.get("sample_alphas", [0.5])), dtype=float).reshape(-1)
    weights = np.asarray(cfg.get("sample_weights", np.ones_like(alphas)), dtype=float).reshape(-1)
    derivative_orders = [int(order) for order in cfg.get("derivative_orders", [0])]
    if alphas.shape != weights.shape:
        raise ValueError("sample_alphas and sample_weights must have the same shape.")

    value = 0.0
    grad_c = np.zeros_like(coeff_blocks)
    grad_T = np.zeros((T.size,), dtype=float)
    for seg_idx, duration in enumerate(T):
        coeff = coeff_blocks[seg_idx]
        for alpha, weight in zip(alphas, weights):
            tau = float(alpha * duration)
            op = _stacked_operator(tau, derivative_orders)
            z = op @ coeff
            penalty = penalty_builder(z, float(duration), cfg)
            grad_z = np.asarray(penalty["grad_z"], dtype=float).reshape(-1)
            dz_dT = np.asarray(
                [float(alpha * (beta_d(tau, order + 1) @ coeff)) for order in derivative_orders],
                dtype=float,
            )
            value += float(weight) * float(penalty["value"])
            grad_c[seg_idx] += float(weight) * (op.T @ grad_z)
            grad_T[seg_idx] += float(weight) * (grad_z @ dz_dT + float(penalty.get("grad_T", 0.0)))
    return {"value": float(value), "grad_c_blocks": grad_c, "grad_T": grad_T}


def obstacle_penalty_full(coeff_blocks: np.ndarray, T: np.ndarray, obs_config: dict[str, Any], s: int = DEFAULT_S) -> dict[str, Any]:
    return _sampled_penalty_common(coeff_blocks, T, obs_config, default_config=default_obs_config(), penalty_builder=_soft_band_penalty, s=s)


def dynamics_penalty_full(coeff_blocks: np.ndarray, T: np.ndarray, dyn_config: dict[str, Any], s: int = DEFAULT_S) -> dict[str, Any]:
    return _sampled_penalty_common(coeff_blocks, T, dyn_config, default_config=default_dyn_config(), penalty_builder=_quadratic_excess_penalty, s=s)


def boundary_penalty_full(coeff_blocks: np.ndarray, T: np.ndarray, bc_config: dict[str, Any], s: int = DEFAULT_S) -> dict[str, Any]:
    coeff_blocks = np.asarray(coeff_blocks, dtype=float)
    T = _validate_time_vector(T)
    cfg = {**default_bc_config(), **({} if bc_config is None else dict(bc_config))}
    orders = [int(order) for order in cfg.get("orders", [0, 1, 2, 3])]
    desired_start = np.asarray(cfg.get("desired_start"), dtype=float).reshape(-1)
    desired_end = np.asarray(cfg.get("desired_end"), dtype=float).reshape(-1)
    if desired_start.size != len(orders) or desired_end.size != len(orders):
        raise ValueError("desired_start and desired_end must match the number of bc orders.")
    weight_start = float(cfg.get("weight_start", 1.0))
    weight_end = float(cfg.get("weight_end", 1.0))

    value = 0.0
    grad_c = np.zeros_like(coeff_blocks)
    grad_T = np.zeros((T.size,), dtype=float)

    start_op = _stacked_operator(0.0, orders)
    start_res = start_op @ coeff_blocks[0] - desired_start
    value += 0.5 * weight_start * float(start_res @ start_res)
    grad_c[0] += weight_start * (start_op.T @ start_res)

    end_tau = float(T[-1])
    end_op = _stacked_operator(end_tau, orders)
    end_res = end_op @ coeff_blocks[-1] - desired_end
    value += 0.5 * weight_end * float(end_res @ end_res)
    grad_c[-1] += weight_end * (end_op.T @ end_res)
    d_end_dT = np.asarray([float(beta_d(end_tau, order + 1) @ coeff_blocks[-1]) for order in orders], dtype=float)
    grad_T[-1] += weight_end * float(end_res @ d_end_dT)

    return {"value": float(value), "grad_c_blocks": grad_c, "grad_T": grad_T}


def regularization_penalty_full(T: np.ndarray, reg_config: dict[str, Any] | None = None) -> dict[str, Any]:
    T = _validate_time_vector(T)
    if T.size < 2:
        return {"value": 0.0, "grad_T": np.zeros_like(T)}
    diff = T[1:] - T[:-1]
    value = 0.5 * float(np.sum(diff**2))
    grad = np.zeros_like(T)
    grad[:-1] -= diff
    grad[1:] += diff
    return {"value": value, "grad_T": grad}


def evaluate_full_objective_from_coeffs(
    coeff_blocks: np.ndarray,
    T: np.ndarray,
    weights: dict[str, float],
    *,
    q_boundary: tuple[float, float] | None = None,
    obs_config: dict[str, Any] | None = None,
    dyn_config: dict[str, Any] | None = None,
    bc_config: dict[str, Any] | None = None,
    reg_config: dict[str, Any] | None = None,
    psi_type: str = "linear",
    psi_config: dict[str, Any] | None = None,
    s: int = DEFAULT_S,
) -> dict[str, Any]:
    coeff_blocks = np.asarray(coeff_blocks, dtype=float)
    T = _validate_time_vector(T)
    weights = {**default_weights(), **dict(weights)}

    ctrl = control_cost_full(coeff_blocks, T, s=s)
    time_term = time_penalty_full(T, psi_type=psi_type, psi_config=psi_config)
    obs_raw = obstacle_penalty_full(coeff_blocks, T, default_obs_config() if obs_config is None else obs_config, s=s)
    dyn_raw = dynamics_penalty_full(coeff_blocks, T, default_dyn_config() if dyn_config is None else dyn_config, s=s)
    if bc_config is None:
        if q_boundary is not None:
            q0, qM = q_boundary
            bc_cfg = default_bc_config(np.asarray([q0, qM], dtype=float))
            bc_cfg["desired_start"][0] = float(q0)
            bc_cfg["desired_end"][0] = float(qM)
        else:
            bc_cfg = default_bc_config()
    else:
        bc_cfg = dict(bc_config)
    bc_raw = boundary_penalty_full(coeff_blocks, T, bc_cfg, s=s)
    reg_raw = regularization_penalty_full(T, reg_config=reg_config)

    g_c_blocks = (
        ctrl["grad_c_blocks"]
        + float(weights["lambda_obs"]) * obs_raw["grad_c_blocks"]
        + float(weights["lambda_dyn"]) * dyn_raw["grad_c_blocks"]
        + float(weights["lambda_bc"]) * bc_raw["grad_c_blocks"]
    )
    g_T = (
        ctrl["grad_T"]
        + float(weights["lambda_T"]) * time_term["grad_T"]
        + float(weights["lambda_obs"]) * obs_raw["grad_T"]
        + float(weights["lambda_dyn"]) * dyn_raw["grad_T"]
        + float(weights["lambda_bc"]) * bc_raw["grad_T"]
        + float(weights["lambda_reg"]) * reg_raw["grad_T"]
    )
    parts = {
        "ctrl": float(ctrl["value"]),
        "time": float(weights["lambda_T"]) * float(time_term["value"]),
        "obs": float(weights["lambda_obs"]) * float(obs_raw["value"]),
        "dyn": float(weights["lambda_dyn"]) * float(dyn_raw["value"]),
        "bc": float(weights["lambda_bc"]) * float(bc_raw["value"]),
        "reg": float(weights["lambda_reg"]) * float(reg_raw["value"]),
    }
    return {
        "value": float(sum(parts.values())),
        "coeff_blocks": coeff_blocks,
        "coeff_vec": coeff_blocks.reshape(-1),
        "g_c_blocks": g_c_blocks,
        "g_c_vec": g_c_blocks.reshape(-1),
        "g_T": g_T,
        "parts": parts,
        "raw_parts": {
            "ctrl": ctrl,
            "time": time_term,
            "obs": obs_raw,
            "dyn": dyn_raw,
            "bc": bc_raw,
            "reg": reg_raw,
        },
    }


def evaluate_full_objective(
    q: np.ndarray,
    T: np.ndarray,
    weights: dict[str, float],
    obs_config: dict[str, Any] | None = None,
    dyn_config: dict[str, Any] | None = None,
    bc_config: dict[str, Any] | None = None,
    reg_config: dict[str, Any] | None = None,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    del config
    q, T = _validate_problem_inputs(q, T, s=s)
    coeffs = compute_raw_schemeC_coeffs_general_k(q, T, s=s, k=k)
    return evaluate_full_objective_from_coeffs(
        coeffs["coeff_blocks"],
        T,
        weights,
        q_boundary=(float(q[0]), float(q[-1])),
        obs_config=obs_config,
        dyn_config=dyn_config,
        bc_config=bc_config,
        reg_config=reg_config,
        s=s,
    )


def full_backward_diff_dense(
    q: np.ndarray,
    T: np.ndarray,
    weights: dict[str, float],
    obs_config: dict[str, Any] | None = None,
    dyn_config: dict[str, Any] | None = None,
    bc_config: dict[str, Any] | None = None,
    reg_config: dict[str, Any] | None = None,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    objective = evaluate_full_objective(q, T, weights, obs_config=obs_config, dyn_config=dyn_config, bc_config=bc_config, reg_config=reg_config, s=s, k=k, config=config)
    jac = compute_raw_schemeC_jacobians_general_k(q, T, s=s, k=k, mode="analytic", config=config)
    grad_q = jac["J_c_q_dense"].T @ objective["g_c_vec"]
    grad_T = jac["J_c_T_dense"].T @ objective["g_c_vec"] + objective["g_T"]
    return {
        "grad_q": grad_q,
        "grad_T": grad_T,
        "grad_theta": np.concatenate((grad_q, grad_T)),
        "objective": objective,
        "jacobians": jac,
    }


def full_backward_diff_sparse(
    q: np.ndarray,
    T: np.ndarray,
    weights: dict[str, float],
    obs_config: dict[str, Any] | None = None,
    dyn_config: dict[str, Any] | None = None,
    bc_config: dict[str, Any] | None = None,
    reg_config: dict[str, Any] | None = None,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    objective = evaluate_full_objective(q, T, weights, obs_config=obs_config, dyn_config=dyn_config, bc_config=bc_config, reg_config=reg_config, s=s, k=k, config=config)
    jac = compute_raw_schemeC_jacobians_general_k(q, T, s=s, k=k, mode="analytic", config=config)
    g_c_blocks = np.asarray(objective["g_c_blocks"], dtype=float)
    M = T.size
    grad_q = np.zeros((max(M - 1, 0),), dtype=float)
    grad_T = np.asarray(objective["g_T"], dtype=float).copy()

    q_to_segments: dict[int, list[int]] = {j: [] for j in range(1, M)}
    T_to_segments: dict[int, list[int]] = {j: [] for j in range(1, M + 1)}
    for seg_idx, support in jac["support_q"].items():
        for q_idx in support:
            q_to_segments[q_idx].append(seg_idx)
    for seg_idx, support in jac["support_T"].items():
        for T_idx in support:
            T_to_segments[T_idx].append(seg_idx)

    for q_idx, segments in q_to_segments.items():
        total = 0.0
        for seg_idx in segments:
            block = jac["J_c_q_blocks"].get((seg_idx, q_idx))
            if block is not None:
                total += float(block.reshape(-1) @ g_c_blocks[seg_idx - 1])
        grad_q[q_idx - 1] = total

    for T_idx, segments in T_to_segments.items():
        for seg_idx in segments:
            block = jac["J_c_T_blocks"].get((seg_idx, T_idx))
            if block is not None:
                grad_T[T_idx - 1] += float(block.reshape(-1) @ g_c_blocks[seg_idx - 1])

    return {
        "grad_q": grad_q,
        "grad_T": grad_T,
        "grad_theta": np.concatenate((grad_q, grad_T)),
        "objective": objective,
        "jacobians": jac,
    }


def full_backward_diff_reparam(
    q: np.ndarray,
    tau: np.ndarray,
    weights: dict[str, float],
    T_min: float = 1e-3,
    obs_config: dict[str, Any] | None = None,
    dyn_config: dict[str, Any] | None = None,
    bc_config: dict[str, Any] | None = None,
    reg_config: dict[str, Any] | None = None,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    q = np.asarray(q, dtype=float).reshape(-1)
    tau = np.asarray(tau, dtype=float).reshape(-1)
    if tau.size != q.size - 1:
        raise ValueError("tau must have one entry per segment.")
    T = tau_to_T(tau, T_min=T_min)
    backward = full_backward_diff_sparse(q, T, weights, obs_config=obs_config, dyn_config=dyn_config, bc_config=bc_config, reg_config=reg_config, s=s, k=k, config=config)
    grad_tau = dT_dtau(tau) * backward["grad_T"]
    grad_xi = np.concatenate((backward["grad_q"], grad_tau))
    return {
        "T": T,
        "grad_q": backward["grad_q"],
        "grad_tau": grad_tau,
        "grad_xi": grad_xi,
        "objective": backward["objective"],
        "dense_sparse_payload": backward,
    }


def compute_support_width_stats(jac: dict[str, Any], M: int) -> dict[str, float]:
    widths = []
    for i in range(1, M + 1):
        q_support = jac["support_q"].get(i, [])
        T_support = jac["support_T"].get(i, [])
        widths.append(float(max(len(q_support), len(T_support))))
    return {"support_width_mean": float(np.mean(widths)), "support_width_max": float(np.max(widths))}


def finite_difference_gradient(f: Callable[[np.ndarray], float], x: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    x = np.asarray(x, dtype=float).reshape(-1)
    grad = np.zeros_like(x)
    for idx in range(x.size):
        step = np.zeros_like(x)
        step[idx] = eps
        grad[idx] = (float(f(x + step)) - float(f(x - step))) / (2.0 * eps)
    return grad


def _plot_dense_vs_sparse(dense: dict[str, Any], sparse: dict[str, Any], save_path: str | Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8.8, 5.8))
    axes[0].plot(dense["grad_q"], marker="o", label="dense")
    axes[0].plot(sparse["grad_q"], marker="x", linestyle="--", label="sparse")
    axes[0].set_title("Phase 10 q-gradient")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend()
    axes[1].plot(dense["grad_T"], marker="o", label="dense")
    axes[1].plot(sparse["grad_T"], marker="x", linestyle="--", label="sparse")
    axes[1].set_title("Phase 10 T-gradient")
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
    plt.title("Phase 10 objective block gradient heatmap")
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_jacobian_sparsity_k_sweep(sparsity_rows: list[dict[str, Any]], save_path: str | Path) -> None:
    valid = [row for row in sparsity_rows if row["implemented"]]
    plt.figure(figsize=(8.4, 4.8))
    if valid:
        plt.plot([row["k"] for row in valid], [row["nnz_q"] for row in valid], marker="o", label="J_c_q nnz")
        plt.plot([row["k"] for row in valid], [row["nnz_T"] for row in valid], marker="s", label="J_c_T nnz")
    for row in sparsity_rows:
        if not row["implemented"]:
            plt.axvline(row["k"], color="tab:red", linestyle=":", alpha=0.5)
    plt.xlabel("k")
    plt.ylabel("nonzero entries")
    plt.title("Phase 10 Jacobian sparsity across k")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def _plot_support_width_vs_k(sparsity_rows: list[dict[str, Any]], save_path: str | Path) -> None:
    valid = [row for row in sparsity_rows if row["implemented"]]
    plt.figure(figsize=(8.4, 4.8))
    if valid:
        plt.plot([row["k"] for row in valid], [row["support_width_mean"] for row in valid], marker="o", label="mean support width")
        plt.plot([row["k"] for row in valid], [row["support_width_max"] for row in valid], marker="s", label="max support width")
    plt.xlabel("k")
    plt.ylabel("support width")
    plt.title("Phase 10 support width vs k")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def run_full_backward_diff_demo(
    q: np.ndarray | None = None,
    T: np.ndarray | None = None,
    weights: dict[str, float] | None = None,
    obs_config: dict[str, Any] | None = None,
    dyn_config: dict[str, Any] | None = None,
    bc_config: dict[str, Any] | None = None,
    reg_config: dict[str, Any] | None = None,
    k_values: list[int] | None = None,
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    del seed
    q, T = representative_case() if q is None or T is None else _validate_problem_inputs(q, T)
    weights = {**default_weights(), **({} if weights is None else dict(weights))}
    obs_cfg = default_obs_config() if obs_config is None else {**default_obs_config(), **dict(obs_config)}
    dyn_cfg = default_dyn_config() if dyn_config is None else {**default_dyn_config(), **dict(dyn_config)}
    bc_cfg = default_bc_config(q) if bc_config is None else {**default_bc_config(q), **dict(bc_config)}
    reg_cfg = default_reg_config() if reg_config is None else {**default_reg_config(), **dict(reg_config)}
    results_dirs = ensure_results_dirs(DEFAULT_RESULTS_DIR if save_dir is None else save_dir)
    back_dir = results_dirs["backward"]
    k_values = [2, 4, 6, 8] if k_values is None else list(k_values)

    dense = full_backward_diff_dense(q, T, weights, obs_config=obs_cfg, dyn_config=dyn_cfg, bc_config=bc_cfg, reg_config=reg_cfg, k=2)
    sparse = full_backward_diff_sparse(q, T, weights, obs_config=obs_cfg, dyn_config=dyn_cfg, bc_config=bc_cfg, reg_config=reg_cfg, k=2)
    tau = T_to_tau(T, T_min=1e-3)
    xi = q_tau_to_xi(q, tau)
    q0 = float(q[0])
    qM = float(q[-1])

    def xi_objective(xi_vec: np.ndarray) -> float:
        q_bar, tau_vec = xi_to_qbar_tau(xi_vec)
        q_full = qbar_to_full(q_bar, q0, qM)
        T_full = tau_to_T(tau_vec, T_min=1e-3)
        return float(evaluate_full_objective(q_full, T_full, weights, obs_config=obs_cfg, dyn_config=dyn_cfg, bc_config=bc_cfg, reg_config=reg_cfg, k=2)["value"])

    fd_grad = finite_difference_gradient(xi_objective, xi)
    reparam = full_backward_diff_reparam(q, tau, weights, T_min=1e-3, obs_config=obs_cfg, dyn_config=dyn_cfg, bc_config=bc_cfg, reg_config=reg_cfg, k=2)
    fd_error = reparam["grad_xi"] - fd_grad

    _plot_dense_vs_sparse(dense, sparse, back_dir / "phase10_dense_vs_sparse_grad.png")
    _plot_block_grad_heatmap(dense["objective"]["g_c_blocks"], back_dir / "phase10_block_grad_heatmap.png")

    sparsity_rows: list[dict[str, Any]] = []
    for k in k_values:
        try:
            jac = compute_raw_schemeC_jacobians_general_k(q, T, k=k)
            width_stats = compute_support_width_stats(jac, T.size)
            sparsity_rows.append(
                {
                    "k": int(k),
                    "implemented": True,
                    "nnz_q": int(np.count_nonzero(np.abs(jac["J_c_q_dense"]) > DEFAULT_ZERO_TOL)),
                    "nnz_T": int(np.count_nonzero(np.abs(jac["J_c_T_dense"]) > DEFAULT_ZERO_TOL)),
                    **width_stats,
                }
            )
        except NotImplementedError:
            sparsity_rows.append(
                {
                    "k": int(k),
                    "implemented": False,
                    "nnz_q": float("nan"),
                    "nnz_T": float("nan"),
                    "support_width_mean": float("nan"),
                    "support_width_max": float("nan"),
                }
            )
    _plot_jacobian_sparsity_k_sweep(sparsity_rows, back_dir / "phase10_jacobian_sparsity_k_sweep.png")
    _plot_support_width_vs_k(sparsity_rows, back_dir / "phase10_support_width_vs_k.png")

    summary = {
        "dense_vs_sparse_max_abs": float(np.max(np.abs(dense["grad_theta"] - sparse["grad_theta"]))),
        "dense_vs_sparse_l2": float(np.linalg.norm(dense["grad_theta"] - sparse["grad_theta"])),
        "reparam_vs_fd_max_abs": float(np.max(np.abs(fd_error))),
        "reparam_vs_fd_l2": float(np.linalg.norm(fd_error)),
        "objective_value": float(dense["objective"]["value"]),
    }
    save_csv(
        back_dir / "phase10_backward_summary.csv",
        ["k", "implemented", "nnz_q", "nnz_T", "support_width_mean", "support_width_max"],
        [[row["k"], row["implemented"], row["nnz_q"], row["nnz_T"], row["support_width_mean"], row["support_width_max"]] for row in sparsity_rows],
    )
    save_json(back_dir / "phase10_backward_summary.json", {"summary": summary, "sparsity_rows": sparsity_rows})
    Path(back_dir / "phase10_backward_summary.md").write_text(
        "\n".join(
            [
                "# Phase 10 Backward Summary",
                "",
                f"- Dense vs sparse max abs gap: `{summary['dense_vs_sparse_max_abs']:.6e}`.",
                f"- Reparameterized gradient vs FD max abs gap: `{summary['reparam_vs_fd_max_abs']:.6e}`.",
                f"- Objective value on representative case: `{summary['objective_value']:.6e}`.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return {
        "dense": dense,
        "sparse": sparse,
        "reparam": reparam,
        "fd_grad": fd_grad,
        "summary": summary,
        "sparsity_rows": sparsity_rows,
    }


def main() -> None:
    run_full_backward_diff_demo()


if __name__ == "__main__":
    main()
