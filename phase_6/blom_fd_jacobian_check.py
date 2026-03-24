"""
Phase 6 Jacobian validation for BLOM basic properties.

This module validates:

- raw BLOM-Analytic local-support sparsity for the canonical (s=4, k=2) map
- analytic-vs-finite-difference Jacobian consistency
- weak time sensitivity through d c / d T statistics
- locality-vs-continuity trade-offs across assembly schemes A/B/C
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

from phase_2.phase2_validation import estimate_effective_bandwidth
from phase_4.blom_k2_s4_numeric import solve_local_system_s4_k2
from phase_4.utils.hermite_utils import (
    D_matrix,
    Lambda8,
    build_C4,
    build_G4,
    build_Pq,
    build_Px,
    build_R4,
    build_m_minus,
    build_m_plus,
    hermite_reconstruct_center_segment,
)
from phase_5.blom_boundary_jump_check import (
    DEFAULT_K,
    DEFAULT_S,
    assemble_scheme_A,
    assemble_scheme_B,
    assemble_scheme_C,
)


DEFAULT_RESULTS_DIR = Path("phase_6/results/phase6_fd_jacobian_check")
DEFAULT_ZERO_TOL = 1e-10
DEFAULT_EPS_LIST = (1e-4, 1e-5, 1e-6, 1e-7)


def _validate_time_vector(T: np.ndarray) -> np.ndarray:
    T = np.asarray(T, dtype=float).reshape(-1)
    if T.size < 3:
        raise ValueError("Phase 6 raw local-support checks require at least 3 segments.")
    if not np.all(np.isfinite(T)):
        raise ValueError("T contains non-finite values.")
    if not np.all(T > 0.0):
        raise ValueError("All durations must be strictly positive.")
    return T


def _validate_problem_inputs(q: np.ndarray, T: np.ndarray, s: int, k: int) -> tuple[np.ndarray, np.ndarray]:
    q = np.asarray(q, dtype=float).reshape(-1)
    T = _validate_time_vector(T)
    if s != 4 or k != 2:
        raise ValueError("Phase 6 currently validates only the canonical case s=4, k=2.")
    if q.shape != (T.size + 1,):
        raise ValueError(f"q must have shape ({T.size + 1},), got {q.shape}.")
    if not np.all(np.isfinite(q)):
        raise ValueError("q contains non-finite values.")
    return q, T


def representative_case() -> tuple[np.ndarray, np.ndarray]:
    q = np.asarray([0.0, 1.05, -0.55, 1.35, -0.25, 0.82, 0.15], dtype=float)
    T = np.asarray([0.85, 1.35, 0.95, 1.25, 0.75, 1.10], dtype=float)
    return q, T


def ensure_results_dirs(base_dir: str | Path = DEFAULT_RESULTS_DIR) -> dict[str, Path]:
    base_dir = Path(base_dir)
    raw_dir = base_dir / "raw_scheme_C"
    scheme_a = base_dir / "scheme_A"
    scheme_b = base_dir / "scheme_B"
    compare_dir = base_dir / "compare"
    for directory in (base_dir, raw_dir, scheme_a, scheme_b, compare_dir):
        directory.mkdir(parents=True, exist_ok=True)
    return {
        "base": base_dir,
        "raw": raw_dir,
        "scheme_A": scheme_a,
        "scheme_B": scheme_b,
        "compare": compare_dir,
    }


def _save_json(path: str | Path, payload: dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _save_csv(path: str | Path, header: list[str], rows: list[list[Any]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [",".join(header)]
    for row in rows:
        rendered = []
        for value in row:
            if isinstance(value, bool):
                rendered.append("true" if value else "false")
            else:
                rendered.append(str(value))
        lines.append(",".join(rendered))
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
    if isinstance(value, list):
        return [_serialize(item) for item in value]
    if isinstance(value, tuple):
        return [_serialize(item) for item in value]
    return value


def _central_interior_segment(M_seg: int) -> int:
    return max(2, min(M_seg - 1, (M_seg + 1) // 2))


def theoretical_mask_c_q(M: int, i: int, s: int = DEFAULT_S, k: int = DEFAULT_K) -> np.ndarray:
    del k
    if s != 4:
        raise ValueError("Phase 6 canonical mask expects s=4.")
    if not 1 <= i <= M:
        raise ValueError(f"i must be in [1, {M}], got {i}.")
    mask = np.zeros((2 * s, max(M - 1, 0)), dtype=bool)
    allowed = {idx for idx in (i - 2, i - 1, i, i + 1) if 1 <= idx <= M - 1}
    for idx in allowed:
        mask[:, idx - 1] = True
    return mask


def theoretical_mask_c_T(M: int, i: int, s: int = DEFAULT_S, k: int = DEFAULT_K) -> np.ndarray:
    del k
    if s != 4:
        raise ValueError("Phase 6 canonical mask expects s=4.")
    if not 1 <= i <= M:
        raise ValueError(f"i must be in [1, {M}], got {i}.")
    mask = np.zeros((2 * s, M), dtype=bool)
    allowed = {idx for idx in (i - 1, i, i + 1) if 1 <= idx <= M}
    for idx in allowed:
        mask[:, idx - 1] = True
    return mask


def theoretical_mask_x_q(M: int, i: int) -> np.ndarray:
    mask = np.zeros((6, max(M - 1, 0)), dtype=bool)
    allowed = {idx for idx in (i - 2, i - 1, i, i + 1) if 1 <= idx <= M - 1}
    for idx in allowed:
        mask[:, idx - 1] = True
    return mask


def theoretical_mask_x_T(M: int, i: int) -> np.ndarray:
    mask = np.zeros((6, M), dtype=bool)
    allowed = {idx for idx in (i - 1, i, i + 1) if 1 <= idx <= M}
    for idx in allowed:
        mask[:, idx - 1] = True
    return mask


def full_theoretical_mask_c_q(M: int, s: int = DEFAULT_S, k: int = DEFAULT_K) -> np.ndarray:
    return np.vstack([theoretical_mask_c_q(M, i, s=s, k=k) for i in range(1, M + 1)])


def full_theoretical_mask_c_T(M: int, s: int = DEFAULT_S, k: int = DEFAULT_K) -> np.ndarray:
    return np.vstack([theoretical_mask_c_T(M, i, s=s, k=k) for i in range(1, M + 1)])


def finite_difference_jacobian(
    f: Callable[[np.ndarray], np.ndarray],
    x: np.ndarray,
    eps: float = 1e-6,
    method: str = "central",
    mask: np.ndarray | None = None,
) -> np.ndarray:
    """Finite-difference Jacobian for vector-valued functions."""
    x = np.asarray(x, dtype=float).reshape(-1)
    base = np.asarray(f(x), dtype=float).reshape(-1)
    J = np.zeros((base.size, x.size), dtype=float)
    if mask is None:
        active = np.ones((x.size,), dtype=bool)
    else:
        active = np.asarray(mask, dtype=bool).reshape(-1)
        if active.shape != (x.size,):
            raise ValueError(f"mask must have shape ({x.size},), got {active.shape}.")
    if eps <= 0.0:
        raise ValueError("eps must be positive.")

    for idx in range(x.size):
        if not active[idx]:
            continue
        step = np.zeros_like(x)
        step[idx] = eps
        if method == "central":
            forward = np.asarray(f(x + step), dtype=float).reshape(-1)
            backward = np.asarray(f(x - step), dtype=float).reshape(-1)
            J[:, idx] = (forward - backward) / (2.0 * eps)
        elif method == "forward":
            forward = np.asarray(f(x + step), dtype=float).reshape(-1)
            J[:, idx] = (forward - base) / eps
        else:
            raise ValueError(f"Unsupported finite-difference method {method!r}.")
    return J


def compute_jacobian_errors(J_ref: np.ndarray, J_test: np.ndarray) -> dict[str, Any]:
    """Return standard norm comparisons between two Jacobians."""
    J_ref = np.asarray(J_ref, dtype=float)
    J_test = np.asarray(J_test, dtype=float)
    if J_ref.shape != J_test.shape:
        raise ValueError(f"Jacobians must have the same shape, got {J_ref.shape} vs {J_test.shape}.")
    diff = J_ref - J_test
    ref_norm = float(np.linalg.norm(J_ref))
    test_norm = float(np.linalg.norm(J_test))
    diff_norm = float(np.linalg.norm(diff))
    return {
        "fro_error": diff_norm,
        "max_abs_error": float(np.max(np.abs(diff))) if diff.size else 0.0,
        "relative_error": diff_norm / max(ref_norm, 1e-15),
        "ref_fro_norm": ref_norm,
        "test_fro_norm": test_norm,
    }


def _d_D_matrix(h: float) -> np.ndarray:
    h = float(h)
    return np.diag([1.0, 2.0 * h, 3.0 * h * h]).astype(float)


def _d_Lambda8(h: float) -> np.ndarray:
    h = float(h)
    diagonal = [0.0]
    diagonal.extend((-power) * h ** (-power - 1) for power in range(1, 8))
    return np.diag(diagonal).astype(float)


def _d_build_m_plus(h: float) -> np.ndarray:
    h = float(h)
    return np.asarray([1.0, h, 0.5 * h * h], dtype=float)


def _d_build_m_minus(h: float) -> np.ndarray:
    h = float(h)
    return np.asarray([1.0, -h, 0.5 * h * h], dtype=float)


def _d_build_Px(h: float) -> np.ndarray:
    mat = np.zeros((8, 6), dtype=float)
    dD = _d_D_matrix(h)
    mat[1:4, 0:3] = dD
    mat[5:8, 3:6] = dD
    return mat


def _d_build_H_mid(h: float) -> np.ndarray:
    h = float(h)
    R4 = build_R4()
    Px = build_Px(h)
    dPx = _d_build_Px(h)
    return 2.0 * ((-7.0) * h ** -8 * (Px.T @ R4 @ Px) + h ** -7 * (dPx.T @ R4 @ Px + Px.T @ R4 @ dPx))


def _d_build_g_mid(h: float) -> np.ndarray:
    h = float(h)
    R4 = build_R4()
    Px = build_Px(h)
    dPx = _d_build_Px(h)
    Pq = build_Pq()
    return -2.0 * ((-7.0) * h ** -8 * (Px.T @ R4 @ Pq) + h ** -7 * (dPx.T @ R4 @ Pq))


def _d_outer_rank_one_hessian(side: str, h: float) -> np.ndarray:
    h = float(h)
    if side == "left":
        vec = build_m_minus(h)
        dvec = _d_build_m_minus(h)
    elif side == "right":
        vec = build_m_plus(h)
        dvec = _d_build_m_plus(h)
    else:
        raise ValueError(f"Unsupported side {side!r}.")
    return 504.0 * ((-7.0) * h ** -8 * np.outer(vec, vec) + h ** -7 * (np.outer(dvec, vec) + np.outer(vec, dvec)))


def _d_outer_linear_term(side: str, q_left: float, q_right: float, h: float) -> np.ndarray:
    delta = float(q_right - q_left)
    h = float(h)
    if side == "left":
        vec = build_m_minus(h)
        dvec = _d_build_m_minus(h)
    elif side == "right":
        vec = build_m_plus(h)
        dvec = _d_build_m_plus(h)
    else:
        raise ValueError(f"Unsupported side {side!r}.")
    return 504.0 * delta * ((-7.0) * h ** -8 * vec + h ** -7 * dvec)


def _qT_to_theta(q: np.ndarray, T: np.ndarray) -> np.ndarray:
    q = np.asarray(q, dtype=float).reshape(-1)
    T = np.asarray(T, dtype=float).reshape(-1)
    return np.concatenate((q[1:-1], T))


def _theta_to_qT(theta: np.ndarray, q0: float, qM: float, M_seg: int) -> tuple[np.ndarray, np.ndarray]:
    theta = np.asarray(theta, dtype=float).reshape(-1)
    if theta.shape != (2 * M_seg - 1,):
        raise ValueError(f"theta must have shape ({2 * M_seg - 1},), got {theta.shape}.")
    q_interior = theta[: M_seg - 1]
    T = theta[M_seg - 1 :]
    q = np.concatenate(([float(q0)], q_interior, [float(qM)]))
    return q, T


def _raw_local_data(q: np.ndarray, T: np.ndarray, i: int) -> tuple[np.ndarray, np.ndarray]:
    if not 2 <= i <= T.size - 1:
        raise ValueError(f"raw local analytic map is implemented only for interior segments 2..{T.size - 1}, got {i}.")
    q_local = np.asarray([q[i - 2], q[i - 1], q[i], q[i + 1]], dtype=float)
    T_local = np.asarray([T[i - 2], T[i - 1], T[i]], dtype=float)
    return q_local, T_local


def raw_local_coefficient_map(q: np.ndarray, T: np.ndarray, i: int, s: int = DEFAULT_S, k: int = DEFAULT_K) -> np.ndarray:
    q, T = _validate_problem_inputs(q, T, s=s, k=k)
    q_local, T_local = _raw_local_data(q, T, i)
    solved = solve_local_system_s4_k2(*q_local, *T_local)
    coeff = hermite_reconstruct_center_segment(
        q_local[1],
        q_local[2],
        solved["x_opt"][:3],
        solved["x_opt"][3:],
        T_local[1],
    )["coeff"]
    return np.asarray(coeff, dtype=float)


def raw_local_jet_state_map(q: np.ndarray, T: np.ndarray, i: int, s: int = DEFAULT_S, k: int = DEFAULT_K) -> np.ndarray:
    q, T = _validate_problem_inputs(q, T, s=s, k=k)
    q_local, T_local = _raw_local_data(q, T, i)
    solved = solve_local_system_s4_k2(*q_local, *T_local)
    return np.asarray(solved["x_opt"], dtype=float)


def _raw_local_analytic_jacobians_local(
    q_local: np.ndarray,
    T_local: np.ndarray,
) -> dict[str, np.ndarray]:
    q_local = np.asarray(q_local, dtype=float).reshape(4)
    T_local = np.asarray(T_local, dtype=float).reshape(3)
    q_im2, q_im1, q_i, q_ip1 = q_local
    T_im1, T_i, T_ip1 = T_local

    solved = solve_local_system_s4_k2(q_im2, q_im1, q_i, q_ip1, T_im1, T_i, T_ip1)
    A2 = solved["A2"]
    x_opt = solved["x_opt"]
    g_mid = solved["g_mid"]
    H4 = np.linalg.inv(build_C4())
    D_i = D_matrix(T_i)
    Lambda_i = Lambda8(T_i)
    y = np.concatenate(([q_im1], D_i @ x_opt[:3], [q_i], D_i @ x_opt[3:]))

    left_scale = 504.0 * T_im1 ** -7
    right_scale = 504.0 * T_ip1 ** -7
    J_B_q_local = np.zeros((6, 4), dtype=float)
    J_B_q_local[:3, 0] = -left_scale * build_m_minus(T_im1)
    J_B_q_local[:, 1] += g_mid[:, 0]
    J_B_q_local[:3, 1] += left_scale * build_m_minus(T_im1)
    J_B_q_local[:, 2] += g_mid[:, 1]
    J_B_q_local[3:, 2] += -right_scale * build_m_plus(T_ip1)
    J_B_q_local[3:, 3] = right_scale * build_m_plus(T_ip1)

    J_x_q_local = np.linalg.solve(A2, J_B_q_local)

    J_y_q_local = np.zeros((8, 4), dtype=float)
    J_y_q_local[0, 1] = 1.0
    J_y_q_local[4, 2] = 1.0
    J_y_q_local[1:4, :] = D_i @ J_x_q_local[:3, :]
    J_y_q_local[5:8, :] = D_i @ J_x_q_local[3:, :]
    J_c_q_local = Lambda_i @ H4 @ J_y_q_local

    dA_left = np.zeros((6, 6), dtype=float)
    dA_left[:3, :3] = _d_outer_rank_one_hessian("left", T_im1)
    dA_mid = _d_build_H_mid(T_i)
    dA_right = np.zeros((6, 6), dtype=float)
    dA_right[3:, 3:] = _d_outer_rank_one_hessian("right", T_ip1)
    dA_local = [dA_left, dA_mid, dA_right]

    dB_left = np.zeros((6,), dtype=float)
    dB_left[:3] = _d_outer_linear_term("left", q_im2, q_im1, T_im1)
    dB_mid = _d_build_g_mid(T_i) @ np.asarray([q_im1, q_i], dtype=float)
    dB_right = np.zeros((6,), dtype=float)
    dB_right[3:] = _d_outer_linear_term("right", q_i, q_ip1, T_ip1)
    dB_local = [dB_left, dB_mid, dB_right]

    J_x_T_local = np.column_stack([np.linalg.solve(A2, dB_local[col] - dA_local[col] @ x_opt) for col in range(3)])

    J_y_T_local = np.zeros((8, 3), dtype=float)
    J_y_T_local[1:4, :] = D_i @ J_x_T_local[:3, :]
    J_y_T_local[5:8, :] = D_i @ J_x_T_local[3:, :]
    J_y_T_local[1:4, 1] += _d_D_matrix(T_i) @ x_opt[:3]
    J_y_T_local[5:8, 1] += _d_D_matrix(T_i) @ x_opt[3:]
    J_c_T_local = Lambda_i @ H4 @ J_y_T_local
    J_c_T_local[:, 1] += _d_Lambda8(T_i) @ H4 @ y

    return {
        "x": x_opt,
        "c": Lambda_i @ H4 @ y,
        "J_x_q_local": J_x_q_local,
        "J_x_T_local": J_x_T_local,
        "J_c_q_local": J_c_q_local,
        "J_c_T_local": J_c_T_local,
    }


def raw_local_jacobians(
    q: np.ndarray,
    T: np.ndarray,
    i: int,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    mode: str = "analytic",
    eps: float = 1e-6,
) -> dict[str, np.ndarray]:
    q, T = _validate_problem_inputs(q, T, s=s, k=k)
    M = T.size
    q_local, T_local = _raw_local_data(q, T, i)

    if mode == "analytic":
        local = _raw_local_analytic_jacobians_local(q_local, T_local)
        J_c_q = np.zeros((2 * s, max(M - 1, 0)), dtype=float)
        J_c_T = np.zeros((2 * s, M), dtype=float)
        J_x_q = np.zeros((6, max(M - 1, 0)), dtype=float)
        J_x_T = np.zeros((6, M), dtype=float)

        q_indices = [i - 2, i - 1, i, i + 1]
        for local_col, q_idx in enumerate(q_indices):
            if 1 <= q_idx <= M - 1:
                J_c_q[:, q_idx - 1] = local["J_c_q_local"][:, local_col]
                J_x_q[:, q_idx - 1] = local["J_x_q_local"][:, local_col]
        t_indices = [i - 1, i, i + 1]
        for local_col, t_idx in enumerate(t_indices):
            if 1 <= t_idx <= M:
                J_c_T[:, t_idx - 1] = local["J_c_T_local"][:, local_col]
                J_x_T[:, t_idx - 1] = local["J_x_T_local"][:, local_col]
        return {
            "c": local["c"],
            "x": local["x"],
            "J_c_q": J_c_q,
            "J_c_T": J_c_T,
            "J_x_q": J_x_q,
            "J_x_T": J_x_T,
        }

    if mode == "fd":
        q0 = float(q[0])
        qM = float(q[-1])
        theta = _qT_to_theta(q, T)
        q_mask = np.zeros_like(theta, dtype=bool)
        q_mask[: M - 1] = True
        T_mask = np.zeros_like(theta, dtype=bool)
        T_mask[M - 1 :] = True

        def coeff_map(theta_vec: np.ndarray) -> np.ndarray:
            q_full, T_full = _theta_to_qT(theta_vec, q0, qM, M)
            return raw_local_coefficient_map(q_full, T_full, i, s=s, k=k)

        def jet_map(theta_vec: np.ndarray) -> np.ndarray:
            q_full, T_full = _theta_to_qT(theta_vec, q0, qM, M)
            return raw_local_jet_state_map(q_full, T_full, i, s=s, k=k)

        J_c = finite_difference_jacobian(coeff_map, theta, eps=eps, method="central")
        J_x = finite_difference_jacobian(jet_map, theta, eps=eps, method="central")
        return {
            "c": coeff_map(theta),
            "x": jet_map(theta),
            "J_c_q": J_c[:, : M - 1],
            "J_c_T": J_c[:, M - 1 :],
            "J_x_q": J_x[:, : M - 1],
            "J_x_T": J_x[:, M - 1 :],
        }

    raise ValueError(f"Unsupported raw Jacobian mode {mode!r}.")


def _stack_raw_jacobians(q: np.ndarray, T: np.ndarray, mode: str, eps: float = 1e-6) -> dict[str, np.ndarray]:
    M = T.size
    rows_c_q: list[np.ndarray] = []
    rows_c_T: list[np.ndarray] = []
    rows_x_q: list[np.ndarray] = []
    rows_x_T: list[np.ndarray] = []
    for i in range(2, M):
        item = raw_local_jacobians(q, T, i, mode=mode, eps=eps)
        rows_c_q.append(item["J_c_q"])
        rows_c_T.append(item["J_c_T"])
        rows_x_q.append(item["J_x_q"])
        rows_x_T.append(item["J_x_T"])
    if not rows_c_q:
        raise ValueError("Need at least one interior segment to stack raw local Jacobians.")
    return {
        "J_c_q": np.vstack(rows_c_q),
        "J_c_T": np.vstack(rows_c_T),
        "J_x_q": np.vstack(rows_x_q),
        "J_x_T": np.vstack(rows_x_T),
    }


def _autodiff_status(use_autodiff: bool) -> dict[str, Any]:
    status = {"requested": bool(use_autodiff), "available": False, "backend": None, "reason": ""}
    if not use_autodiff:
        status["reason"] = "autodiff disabled by caller"
        return status
    try:
        __import__("jax")
        status["available"] = True
        status["backend"] = "jax"
        return status
    except Exception:
        pass
    try:
        __import__("autograd")
        status["available"] = True
        status["backend"] = "autograd"
        return status
    except Exception:
        pass
    status["reason"] = "neither jax nor autograd is installed in the execution environment"
    return status


def assembled_scheme_A_map(q: np.ndarray, T: np.ndarray, i: int, config: dict[str, Any] | None = None) -> np.ndarray:
    result = assemble_scheme_A(q, T, config=config)
    return np.asarray(result["coeffs"][i - 1], dtype=float)


def assembled_scheme_B_map(q: np.ndarray, T: np.ndarray, i: int, config: dict[str, Any] | None = None) -> np.ndarray:
    result = assemble_scheme_B(q, T, config=config)
    return np.asarray(result["coeffs"][i - 1], dtype=float)


def assembled_scheme_C_map(q: np.ndarray, T: np.ndarray, i: int, config: dict[str, Any] | None = None) -> np.ndarray:
    result = assemble_scheme_C(q, T, config=config)
    return np.asarray(result["coeffs"][i - 1], dtype=float)


def _assemble_full_coeffs(q: np.ndarray, T: np.ndarray, scheme: str, config: dict[str, Any] | None = None) -> np.ndarray:
    scheme = scheme.upper()
    if scheme == "A":
        result = assemble_scheme_A(q, T, config=config)
    elif scheme == "B":
        result = assemble_scheme_B(q, T, config=config)
    elif scheme == "C":
        result = assemble_scheme_C(q, T, config=config)
    else:
        raise ValueError(f"Unsupported assembly scheme {scheme!r}.")
    return np.asarray(result["coeffs"], dtype=float)


def assembled_scheme_jacobians_fd(
    q: np.ndarray,
    T: np.ndarray,
    scheme: str,
    eps_q: float = 1e-6,
    eps_T: float = 1e-6,
    config: dict[str, Any] | None = None,
) -> dict[str, np.ndarray]:
    q, T = _validate_problem_inputs(q, T, s=DEFAULT_S, k=DEFAULT_K)
    M = T.size
    q0 = float(q[0])
    qM = float(q[-1])
    theta = _qT_to_theta(q, T)

    def coeff_map(theta_vec: np.ndarray) -> np.ndarray:
        q_full, T_full = _theta_to_qT(theta_vec, q0, qM, M)
        return _assemble_full_coeffs(q_full, T_full, scheme=scheme, config=config).reshape(-1)

    mask_q = np.zeros_like(theta, dtype=bool)
    mask_q[: M - 1] = True
    mask_T = np.zeros_like(theta, dtype=bool)
    mask_T[M - 1 :] = True

    J_q = finite_difference_jacobian(coeff_map, theta, eps=eps_q, method="central", mask=mask_q)[:, : M - 1]
    J_T = finite_difference_jacobian(coeff_map, theta, eps=eps_T, method="central", mask=mask_T)[:, M - 1 :]
    coeffs = coeff_map(theta).reshape(M, 2 * DEFAULT_S)
    return {"coeffs": coeffs, "J_c_q": J_q, "J_c_T": J_T}


def _segment_parameter_bandwidth(J: np.ndarray, M_seg: int, num_params: int, tol: float) -> dict[str, Any]:
    J = np.asarray(J, dtype=float)
    block_size = 2 * DEFAULT_S
    if J.shape != (block_size * M_seg, num_params):
        raise ValueError(f"Unexpected Jacobian shape {J.shape}; expected ({block_size * M_seg}, {num_params}).")
    per_param = []
    for param_index in range(num_params):
        distances = np.abs(np.arange(M_seg) - param_index)
        significant = distances[np.linalg.norm(J.reshape(M_seg, block_size, num_params)[:, :, param_index], axis=1) > tol]
        per_param.append(int(np.max(significant)) if significant.size else 0)
    return {
        "max_effective_bandwidth": int(max(per_param, default=0)),
        "mean_effective_bandwidth": float(np.mean(per_param)) if per_param else 0.0,
        "per_parameter_bandwidth": per_param,
    }


def _sparsity_rows(
    scheme: str,
    J: np.ndarray,
    theory_mask: np.ndarray,
    jacobian_type: str,
    tol: float,
) -> list[list[Any]]:
    J = np.asarray(J, dtype=float)
    theory_mask = np.asarray(theory_mask, dtype=bool)
    block_size = 2 * DEFAULT_S
    M_seg = J.shape[0] // block_size
    rows = []
    for segment_idx in range(1, M_seg + 1):
        block = J[(segment_idx - 1) * block_size : segment_idx * block_size]
        block_mask = theory_mask[(segment_idx - 1) * block_size : segment_idx * block_size]
        abs_block = np.abs(block)
        in_band = abs_block[block_mask]
        out_band = abs_block[~block_mask]
        rows.append(
            [
                scheme,
                segment_idx,
                jacobian_type,
                int(np.sum(abs_block > tol)),
                int(np.sum(in_band > tol)),
                int(np.sum(out_band > tol)),
                float(np.max(out_band)) if out_band.size else 0.0,
                float(np.mean(out_band)) if out_band.size else 0.0,
            ]
        )
    return rows


def _plot_heatmap(mat: np.ndarray, save_path: str | Path, title: str, xlabel: str, ylabel: str) -> None:
    mat = np.asarray(mat, dtype=float)
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    display = np.log10(np.abs(mat) + 1e-16) if mat.size else np.zeros((1, 1), dtype=float)
    image = ax.imshow(display, origin="lower", aspect="auto", cmap="viridis")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("log10(|value| + 1e-16)")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_binary_vs_heatmap(mask: np.ndarray, mat: np.ndarray, save_path: str | Path, title: str) -> None:
    mask = np.asarray(mask, dtype=float)
    mat = np.asarray(mat, dtype=float)
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].imshow(mask, origin="lower", aspect="auto", cmap="gray_r")
    axes[0].set_title(f"{title}: Theory Mask")
    axes[0].set_xlabel("Parameter Index")
    axes[0].set_ylabel("Coefficient Row")
    image = axes[1].imshow(np.log10(np.abs(mat) + 1e-16), origin="lower", aspect="auto", cmap="magma")
    axes[1].set_title(f"{title}: FD |Jacobian|")
    axes[1].set_xlabel("Parameter Index")
    axes[1].set_ylabel("Coefficient Row")
    colorbar = fig.colorbar(image, ax=axes[1])
    colorbar.set_label("log10(|value| + 1e-16)")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_error_heatmap(mat: np.ndarray, save_path: str | Path, title: str) -> None:
    mat = np.asarray(mat, dtype=float)
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    image = ax.imshow(np.log10(np.abs(mat) + 1e-16), origin="lower", aspect="auto", cmap="magma")
    ax.set_title(title)
    ax.set_xlabel("Parameter Index")
    ax.set_ylabel("Coefficient Row")
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("log10(|analytic - FD| + 1e-16)")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_stepsize_sweep(records: list[dict[str, Any]], save_path: str | Path, title: str) -> None:
    eps = np.asarray([record["eps"] for record in records], dtype=float)
    fro_err = np.asarray([record["fro_error"] for record in records], dtype=float)
    max_err = np.asarray([record["max_abs_error"] for record in records], dtype=float)
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(eps, fro_err, marker="o", linewidth=1.8, label="Frobenius error")
    ax.plot(eps, max_err, marker="s", linewidth=1.8, label="Max abs error")
    ax.set_title(title)
    ax.set_xlabel("Finite-difference step eps")
    ax.set_ylabel("Error")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_bandwidth_compare(records: list[dict[str, Any]], save_path: str | Path) -> None:
    schemes = [record["scheme"] for record in records]
    q_bw = np.asarray([record["q_bandwidth"] for record in records], dtype=float)
    T_bw = np.asarray([record["T_bandwidth"] for record in records], dtype=float)
    x = np.arange(len(schemes))
    width = 0.36
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    ax.bar(x - width / 2, q_bw, width, label="q bandwidth")
    ax.bar(x + width / 2, T_bw, width, label="T bandwidth")
    ax.set_xticks(x)
    ax.set_xticklabels(schemes)
    ax.set_ylabel("Effective Bandwidth")
    ax.set_title("Jacobian Locality Width Comparison")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_boxplot(groups: dict[str, list[float]], save_path: str | Path, title: str, ylabel: str) -> None:
    labels = list(groups)
    values = [groups[label] for label in labels]
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.boxplot(values, tick_labels=labels)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_yscale("log")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _zero_pattern_pass(J: np.ndarray, mask: np.ndarray, tol: float) -> bool:
    outside = np.abs(np.asarray(J, dtype=float))[~np.asarray(mask, dtype=bool)]
    return bool(np.max(outside) <= tol) if outside.size else True


def _stepsize_sweep_for_raw(
    q: np.ndarray,
    T: np.ndarray,
    segment_idx: int,
    analytic: dict[str, np.ndarray],
    eps_list: tuple[float, ...] = DEFAULT_EPS_LIST,
) -> dict[str, list[dict[str, Any]]]:
    q_records = []
    T_records = []
    for eps in eps_list:
        fd = raw_local_jacobians(q, T, segment_idx, mode="fd", eps=eps)
        q_err = compute_jacobian_errors(analytic["J_c_q"], fd["J_c_q"])
        T_err = compute_jacobian_errors(analytic["J_c_T"], fd["J_c_T"])
        q_records.append({"eps": float(eps), **q_err})
        T_records.append({"eps": float(eps), **T_err})
    return {"q": q_records, "T": T_records}


def _raw_support_and_fd_check(
    q: np.ndarray,
    T: np.ndarray,
    use_autodiff: bool,
    save_root: dict[str, Path],
    tol: float = DEFAULT_ZERO_TOL,
) -> dict[str, Any]:
    M = T.size
    segment_idx = _central_interior_segment(M)
    analytic = raw_local_jacobians(q, T, segment_idx, mode="analytic")
    fd = raw_local_jacobians(q, T, segment_idx, mode="fd", eps=1e-6)
    raw_stack_analytic = _stack_raw_jacobians(q, T, mode="analytic")
    raw_stack_fd = _stack_raw_jacobians(q, T, mode="fd", eps=1e-6)

    mask_q_single = theoretical_mask_c_q(M, segment_idx)
    mask_T_single = theoretical_mask_c_T(M, segment_idx)
    mask_x_q_single = theoretical_mask_x_q(M, segment_idx)
    mask_x_T_single = theoretical_mask_x_T(M, segment_idx)
    mask_q_stack = np.vstack([theoretical_mask_c_q(M, i) for i in range(2, M)])
    mask_T_stack = np.vstack([theoretical_mask_c_T(M, i) for i in range(2, M)])

    q_error = compute_jacobian_errors(analytic["J_c_q"], fd["J_c_q"])
    T_error = compute_jacobian_errors(analytic["J_c_T"], fd["J_c_T"])
    x_q_error = compute_jacobian_errors(analytic["J_x_q"], fd["J_x_q"])
    x_T_error = compute_jacobian_errors(analytic["J_x_T"], fd["J_x_T"])
    sweep = _stepsize_sweep_for_raw(q, T, segment_idx, analytic)

    _plot_heatmap(raw_stack_analytic["J_c_q"], save_root["raw"] / "jacobian_mask_q_raw.png", "Raw Scheme C J_c_q", "Interior Waypoint Index", "Stacked Coefficient Row")
    _plot_heatmap(raw_stack_analytic["J_c_T"], save_root["raw"] / "jacobian_mask_T_raw.png", "Raw Scheme C J_c_T", "Duration Index", "Stacked Coefficient Row")
    _plot_binary_vs_heatmap(mask_q_single, fd["J_c_q"], save_root["raw"] / "jacobian_theory_vs_fd_q_raw.png", "Raw J_c_q")
    _plot_binary_vs_heatmap(mask_T_single, fd["J_c_T"], save_root["raw"] / "jacobian_theory_vs_fd_T_raw.png", "Raw J_c_T")
    _plot_error_heatmap(analytic["J_c_q"] - fd["J_c_q"], save_root["compare"] / "jacobian_error_analytic_vs_fd_q.png", "Raw J_c_q Analytic vs FD Error")
    _plot_error_heatmap(analytic["J_c_T"] - fd["J_c_T"], save_root["compare"] / "jacobian_error_analytic_vs_fd_T.png", "Raw J_c_T Analytic vs FD Error")
    _plot_stepsize_sweep(sweep["q"], save_root["compare"] / "fd_stepsize_sweep_q.png", "FD Step Sweep for Raw J_c_q")
    _plot_stepsize_sweep(sweep["T"], save_root["compare"] / "fd_stepsize_sweep_T.png", "FD Step Sweep for Raw J_c_T")

    sparsity_rows = []
    sparsity_rows.extend(_sparsity_rows("raw", raw_stack_fd["J_c_q"], mask_q_stack, "c_q", tol))
    sparsity_rows.extend(_sparsity_rows("raw", raw_stack_fd["J_c_T"], mask_T_stack, "c_T", tol))
    _save_csv(
        save_root["raw"] / "jacobian_sparsity_stats_raw.csv",
        [
            "scheme",
            "segment_idx",
            "jacobian_type",
            "nnz_total",
            "nnz_theory_band",
            "nnz_outside_band",
            "max_abs_outside_band",
            "mean_abs_outside_band",
        ],
        sparsity_rows,
    )

    error_rows = [
        ["raw", segment_idx, "c_q", "fro_error", q_error["fro_error"], math.nan, math.nan, 1e-6],
        ["raw", segment_idx, "c_q", "max_abs_error", q_error["max_abs_error"], math.nan, math.nan, 1e-6],
        ["raw", segment_idx, "c_T", "fro_error", T_error["fro_error"], math.nan, math.nan, 1e-6],
        ["raw", segment_idx, "c_T", "max_abs_error", T_error["max_abs_error"], math.nan, math.nan, 1e-6],
        ["raw", segment_idx, "x_q", "fro_error", x_q_error["fro_error"], math.nan, math.nan, 1e-6],
        ["raw", segment_idx, "x_T", "fro_error", x_T_error["fro_error"], math.nan, math.nan, 1e-6],
    ]

    autodiff = _autodiff_status(use_autodiff)
    result = {
        "segment_idx": segment_idx,
        "analytic": analytic,
        "fd": fd,
        "mask_pass": {
            "J_c_q": _zero_pattern_pass(fd["J_c_q"], mask_q_single, tol=1e-5),
            "J_c_T": _zero_pattern_pass(fd["J_c_T"], mask_T_single, tol=1e-5),
            "J_x_q": _zero_pattern_pass(fd["J_x_q"], mask_x_q_single, tol=1e-5),
            "J_x_T": _zero_pattern_pass(fd["J_x_T"], mask_x_T_single, tol=1e-5),
        },
        "errors": {
            "J_c_q": q_error,
            "J_c_T": T_error,
            "J_x_q": x_q_error,
            "J_x_T": x_T_error,
        },
        "stepsize_sweep": sweep,
        "sparsity_rows": sparsity_rows,
        "error_rows": error_rows,
        "autodiff": autodiff,
        "figure_paths": {
            "jacobian_mask_q_raw": str(save_root["raw"] / "jacobian_mask_q_raw.png"),
            "jacobian_mask_T_raw": str(save_root["raw"] / "jacobian_mask_T_raw.png"),
            "jacobian_theory_vs_fd_q_raw": str(save_root["raw"] / "jacobian_theory_vs_fd_q_raw.png"),
            "jacobian_theory_vs_fd_T_raw": str(save_root["raw"] / "jacobian_theory_vs_fd_T_raw.png"),
            "jacobian_error_analytic_vs_fd_q": str(save_root["compare"] / "jacobian_error_analytic_vs_fd_q.png"),
            "jacobian_error_analytic_vs_fd_T": str(save_root["compare"] / "jacobian_error_analytic_vs_fd_T.png"),
            "fd_stepsize_sweep_q": str(save_root["compare"] / "fd_stepsize_sweep_q.png"),
            "fd_stepsize_sweep_T": str(save_root["compare"] / "fd_stepsize_sweep_T.png"),
        },
        "table_paths": {
            "jacobian_sparsity_stats_raw": str(save_root["raw"] / "jacobian_sparsity_stats_raw.csv"),
        },
    }
    return result


def _scheme_nonlocality_check(
    q: np.ndarray,
    T: np.ndarray,
    scheme: str,
    save_dir: Path,
    tol: float = DEFAULT_ZERO_TOL,
) -> dict[str, Any]:
    scheme = scheme.upper()
    assembled = assembled_scheme_jacobians_fd(q, T, scheme)
    M = T.size
    mask_q = full_theoretical_mask_c_q(M)
    mask_T = full_theoretical_mask_c_T(M)
    _plot_heatmap(
        assembled["J_c_q"],
        save_dir / f"jacobian_mask_q_scheme_{scheme}.png",
        f"Scheme {scheme} Final J_c_q",
        "Interior Waypoint Index",
        "Stacked Coefficient Row",
    )
    _plot_heatmap(
        assembled["J_c_T"],
        save_dir / f"jacobian_mask_T_scheme_{scheme}.png",
        f"Scheme {scheme} Final J_c_T",
        "Duration Index",
        "Stacked Coefficient Row",
    )

    q_bandwidth = estimate_effective_bandwidth(assembled["J_c_q"], M, tol=tol, block_size=2 * DEFAULT_S)
    T_bandwidth = _segment_parameter_bandwidth(assembled["J_c_T"], M, M, tol=tol)
    sparsity_rows = []
    sparsity_rows.extend(_sparsity_rows(f"scheme_{scheme}", assembled["J_c_q"], mask_q, "c_q", tol))
    sparsity_rows.extend(_sparsity_rows(f"scheme_{scheme}", assembled["J_c_T"], mask_T, "c_T", tol))
    csv_name = f"jacobian_sparsity_stats_scheme_{scheme}.csv"
    _save_csv(
        save_dir / csv_name,
        [
            "scheme",
            "segment_idx",
            "jacobian_type",
            "nnz_total",
            "nnz_theory_band",
            "nnz_outside_band",
            "max_abs_outside_band",
            "mean_abs_outside_band",
        ],
        sparsity_rows,
    )
    outside_q = np.abs(assembled["J_c_q"])[~mask_q]
    outside_T = np.abs(assembled["J_c_T"])[~mask_T]
    return {
        "coeffs": assembled["coeffs"],
        "J_c_q": assembled["J_c_q"],
        "J_c_T": assembled["J_c_T"],
        "q_bandwidth": q_bandwidth,
        "T_bandwidth": T_bandwidth,
        "outside_band": {
            "q_nonzero_ratio": float(np.mean(outside_q > tol)) if outside_q.size else 0.0,
            "q_max_abs": float(np.max(outside_q)) if outside_q.size else 0.0,
            "T_nonzero_ratio": float(np.mean(outside_T > tol)) if outside_T.size else 0.0,
            "T_max_abs": float(np.max(outside_T)) if outside_T.size else 0.0,
        },
        "sparsity_rows": sparsity_rows,
        "figure_paths": {
            "jacobian_mask_q": str(save_dir / f"jacobian_mask_q_scheme_{scheme}.png"),
            "jacobian_mask_T": str(save_dir / f"jacobian_mask_T_scheme_{scheme}.png"),
        },
        "table_paths": {
            "jacobian_sparsity_stats": str(save_dir / csv_name),
        },
    }


def run_fd_jacobian_check(
    q: np.ndarray,
    T: np.ndarray,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    scheme: str = "C",
    use_autodiff: bool = True,
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    del seed
    q, T = _validate_problem_inputs(q, T, s=s, k=k)
    scheme = scheme.upper()
    if save_dir is None:
        save_dir = DEFAULT_RESULTS_DIR
    dirs = ensure_results_dirs(save_dir)
    start_time = time.perf_counter()

    raw = _raw_support_and_fd_check(q, T, use_autodiff=use_autodiff, save_root=dirs)
    if scheme == "C":
        scheme_result = _scheme_nonlocality_check(q, T, "C", dirs["raw"])
    elif scheme == "A":
        scheme_result = _scheme_nonlocality_check(q, T, "A", dirs["scheme_A"])
    elif scheme == "B":
        scheme_result = _scheme_nonlocality_check(q, T, "B", dirs["scheme_B"])
    else:
        raise ValueError(f"Unsupported scheme {scheme!r}.")

    elapsed = time.perf_counter() - start_time
    summary = {
        "inputs": {"q": q, "T": T, "s": s, "k": k},
        "scheme": scheme,
        "raw": raw,
        "scheme_result": scheme_result,
        "autodiff": raw["autodiff"],
        "elapsed_sec": elapsed,
    }
    _save_json(Path(save_dir) / "summary_phase6_fd_check.json", _serialize(summary))
    return summary


def run_compare_all_schemes(
    q: np.ndarray,
    T: np.ndarray,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    use_autodiff: bool = True,
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    del seed
    q, T = _validate_problem_inputs(q, T, s=s, k=k)
    if save_dir is None:
        save_dir = DEFAULT_RESULTS_DIR
    dirs = ensure_results_dirs(save_dir)
    start_time = time.perf_counter()

    raw = _raw_support_and_fd_check(q, T, use_autodiff=use_autodiff, save_root=dirs)
    scheme_a = _scheme_nonlocality_check(q, T, "A", dirs["scheme_A"])
    scheme_b = _scheme_nonlocality_check(q, T, "B", dirs["scheme_B"])
    scheme_c = _scheme_nonlocality_check(q, T, "C", dirs["raw"])

    bandwidth_records = [
        {
            "scheme": "C",
            "q_bandwidth": scheme_c["q_bandwidth"]["max_effective_bandwidth"],
            "T_bandwidth": scheme_c["T_bandwidth"]["max_effective_bandwidth"],
        },
        {
            "scheme": "A",
            "q_bandwidth": scheme_a["q_bandwidth"]["max_effective_bandwidth"],
            "T_bandwidth": scheme_a["T_bandwidth"]["max_effective_bandwidth"],
        },
        {
            "scheme": "B",
            "q_bandwidth": scheme_b["q_bandwidth"]["max_effective_bandwidth"],
            "T_bandwidth": scheme_b["T_bandwidth"]["max_effective_bandwidth"],
        },
    ]
    _plot_bandwidth_compare(bandwidth_records, dirs["compare"] / "jacobian_bandwidth_compare.png")

    bound_groups = {
        "T_{i-1}": [float(np.linalg.norm(raw["analytic"]["J_c_T"][:, 0]))],
        "T_i": [float(np.linalg.norm(raw["analytic"]["J_c_T"][:, 1]))],
        "T_{i+1}": [float(np.linalg.norm(raw["analytic"]["J_c_T"][:, 2]))],
    }
    _plot_boxplot(
        bound_groups,
        dirs["compare"] / "dc_dT_bound_statistics.png",
        "Representative ||dc/dT|| Statistics",
        "Jacobian Column Norm",
    )

    error_rows = list(raw["error_rows"])
    _save_csv(
        dirs["compare"] / "jacobian_error_stats.csv",
        [
            "scheme",
            "segment_idx",
            "jacobian_type",
            "error_metric",
            "analytic_vs_fd",
            "ad_vs_fd",
            "analytic_vs_ad",
            "eps",
        ],
        error_rows,
    )

    elapsed = time.perf_counter() - start_time
    result = {
        "inputs": {"q": q, "T": T, "s": s, "k": k},
        "raw": raw,
        "scheme_A": scheme_a,
        "scheme_B": scheme_b,
        "scheme_C": scheme_c,
        "bandwidth_records": bandwidth_records,
        "autodiff": raw["autodiff"],
        "figure_paths": {
            "jacobian_bandwidth_compare": str(dirs["compare"] / "jacobian_bandwidth_compare.png"),
            "dc_dT_bound_statistics": str(dirs["compare"] / "dc_dT_bound_statistics.png"),
            "jacobian_error_boxplot_random_trials": str(dirs["compare"] / "jacobian_error_boxplot_random_trials.png"),
        },
        "table_paths": {
            "jacobian_error_stats": str(dirs["compare"] / "jacobian_error_stats.csv"),
        },
        "elapsed_sec": elapsed,
    }
    _save_json(dirs["compare"] / "summary_phase6_fd_check.json", _serialize(result))
    return result


def run_random_trials(
    n_trials: int = 100,
    M: int = 20,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    use_autodiff: bool = True,
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    if n_trials < 1:
        raise ValueError("n_trials must be >= 1.")
    if M < 4:
        raise ValueError("Phase 6 random trials require M >= 4.")
    if save_dir is None:
        save_dir = DEFAULT_RESULTS_DIR
    dirs = ensure_results_dirs(save_dir)
    rng = np.random.default_rng(seed)

    raw_rows = []
    scheme_rows = []
    bound_norms = []
    error_groups = {"analytic_vs_fd_q": [], "analytic_vs_fd_T": []}
    for trial in range(n_trials):
        q = np.zeros((M + 1,), dtype=float)
        q[1:] = np.cumsum(rng.normal(scale=0.8, size=M))
        T = rng.uniform(0.5, 2.0, size=M).astype(float)
        segment_idx = _central_interior_segment(M)
        analytic = raw_local_jacobians(q, T, segment_idx, mode="analytic")
        fd = raw_local_jacobians(q, T, segment_idx, mode="fd", eps=1e-6)
        q_err = compute_jacobian_errors(analytic["J_c_q"], fd["J_c_q"])
        T_err = compute_jacobian_errors(analytic["J_c_T"], fd["J_c_T"])
        q_mask = theoretical_mask_c_q(M, segment_idx)
        T_mask = theoretical_mask_c_T(M, segment_idx)
        raw_rows.append(
            {
                "trial": trial,
                "mask_pass_q": _zero_pattern_pass(fd["J_c_q"], q_mask, tol=1e-5),
                "mask_pass_T": _zero_pattern_pass(fd["J_c_T"], T_mask, tol=1e-5),
                "analytic_fd_q_max": q_err["max_abs_error"],
                "analytic_fd_T_max": T_err["max_abs_error"],
            }
        )
        error_groups["analytic_vs_fd_q"].append(q_err["max_abs_error"])
        error_groups["analytic_vs_fd_T"].append(T_err["max_abs_error"])
        bound_norms.extend(float(np.linalg.norm(analytic["J_c_T"][:, col])) for col in range(3))

        for scheme in ("A", "B", "C"):
            scheme_jac = _scheme_nonlocality_check(q, T, scheme, dirs["scheme_A"] if scheme == "A" else dirs["scheme_B"] if scheme == "B" else dirs["raw"])
            scheme_rows.append(
                {
                    "trial": trial,
                    "scheme": scheme,
                    "outside_band_q_ratio": scheme_jac["outside_band"]["q_nonzero_ratio"],
                    "outside_band_T_ratio": scheme_jac["outside_band"]["T_nonzero_ratio"],
                    "q_bandwidth": scheme_jac["q_bandwidth"]["max_effective_bandwidth"],
                    "T_bandwidth": scheme_jac["T_bandwidth"]["max_effective_bandwidth"],
                }
            )

    _plot_boxplot(
        error_groups,
        dirs["compare"] / "jacobian_error_boxplot_random_trials.png",
        "Random-Trial Jacobian Error Distribution",
        "Max Abs Error",
    )
    _plot_boxplot(
        {"||dc/dT||": bound_norms},
        dirs["compare"] / "dc_dT_bound_statistics.png",
        "Random-Trial ||dc/dT|| Bound Statistics",
        "Jacobian Column Norm",
    )

    _save_csv(
        dirs["compare"] / "random_trials_raw_stats.csv",
        ["trial", "mask_pass_q", "mask_pass_T", "analytic_fd_q_max", "analytic_fd_T_max"],
        [
            [
                row["trial"],
                row["mask_pass_q"],
                row["mask_pass_T"],
                row["analytic_fd_q_max"],
                row["analytic_fd_T_max"],
            ]
            for row in raw_rows
        ],
    )
    _save_csv(
        dirs["compare"] / "random_trials_scheme_locality.csv",
        ["trial", "scheme", "outside_band_q_ratio", "outside_band_T_ratio", "q_bandwidth", "T_bandwidth"],
        [
            [
                row["trial"],
                row["scheme"],
                row["outside_band_q_ratio"],
                row["outside_band_T_ratio"],
                row["q_bandwidth"],
                row["T_bandwidth"],
            ]
            for row in scheme_rows
        ],
    )

    scheme_summary: dict[str, dict[str, Any]] = {}
    for scheme in ("A", "B", "C"):
        rows = [row for row in scheme_rows if row["scheme"] == scheme]
        scheme_summary[scheme] = {
            "mean_outside_band_q_ratio": float(np.mean([row["outside_band_q_ratio"] for row in rows])),
            "mean_outside_band_T_ratio": float(np.mean([row["outside_band_T_ratio"] for row in rows])),
            "mean_q_bandwidth": float(np.mean([row["q_bandwidth"] for row in rows])),
            "mean_T_bandwidth": float(np.mean([row["T_bandwidth"] for row in rows])),
        }
    raw_summary = {
        "mask_pass_rate_q": float(np.mean([row["mask_pass_q"] for row in raw_rows])),
        "mask_pass_rate_T": float(np.mean([row["mask_pass_T"] for row in raw_rows])),
        "mean_analytic_fd_q_max": float(np.mean([row["analytic_fd_q_max"] for row in raw_rows])),
        "mean_analytic_fd_T_max": float(np.mean([row["analytic_fd_T_max"] for row in raw_rows])),
        "max_dc_dT_norm": float(np.max(bound_norms)) if bound_norms else 0.0,
    }
    summary = {
        "raw": raw_summary,
        "schemes": scheme_summary,
        "autodiff": _autodiff_status(use_autodiff),
        "figure_paths": {
            "jacobian_error_boxplot_random_trials": str(dirs["compare"] / "jacobian_error_boxplot_random_trials.png"),
            "dc_dT_bound_statistics": str(dirs["compare"] / "dc_dT_bound_statistics.png"),
        },
    }
    _save_json(dirs["compare"] / "random_trials_summary.json", _serialize(summary))
    return {"raw_rows": raw_rows, "scheme_rows": scheme_rows, "summary": summary}


def main(results_dir: str | Path = DEFAULT_RESULTS_DIR, seed: int = 42) -> dict[str, Any]:
    q, T = representative_case()
    compare = run_compare_all_schemes(q, T, save_dir=results_dir, seed=seed)
    random_trials = run_random_trials(n_trials=12, M=8, save_dir=results_dir, seed=seed)
    print("Phase 6 FD Jacobian check")
    raw = compare["raw"]
    print(f"raw q max abs error: {raw['errors']['J_c_q']['max_abs_error']:.3e}")
    print(f"raw T max abs error: {raw['errors']['J_c_T']['max_abs_error']:.3e}")
    print(f"scheme C q bandwidth: {compare['scheme_C']['q_bandwidth']['max_effective_bandwidth']}")
    print(f"scheme A q bandwidth: {compare['scheme_A']['q_bandwidth']['max_effective_bandwidth']}")
    print(f"scheme B q bandwidth: {compare['scheme_B']['q_bandwidth']['max_effective_bandwidth']}")
    print(f"random raw mask pass rate q: {random_trials['summary']['raw']['mask_pass_rate_q']:.2f}")
    print(f"results dir: {results_dir}")
    return {"compare": compare, "random_trials": random_trials}


if __name__ == "__main__":
    main()
