"""
Linear-algebra back-end for the Phase 3 BLOM-Strict local QP.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from phase_1.minco_scalar_baseline import beta, beta_d


DEFAULT_S = 4


def _validate_time_vector(T_window: np.ndarray) -> np.ndarray:
    T_window = np.asarray(T_window, dtype=float).reshape(-1)
    if T_window.size < 1:
        raise ValueError("Local window must contain at least one segment.")
    if not np.all(np.isfinite(T_window)):
        raise ValueError("T_window contains non-finite values.")
    if not np.all(T_window > 0.0):
        raise ValueError("All local segment durations must be strictly positive.")
    return T_window


def _normalize_problem_inputs(
    q: np.ndarray,
    T: np.ndarray,
    zeta_start: np.ndarray | None,
    zeta_end: np.ndarray | None,
    s: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    q = np.asarray(q, dtype=float).reshape(-1)
    T = np.asarray(T, dtype=float).reshape(-1)
    if q.shape != (T.size + 1,):
        raise ValueError(f"q must have shape ({T.size + 1},), got {q.shape}.")
    if not np.all(np.isfinite(q)):
        raise ValueError("q contains non-finite values.")
    T = _validate_time_vector(T)

    if zeta_start is None:
        zeta_start = np.zeros((s,), dtype=float)
        zeta_start[0] = q[0]
    else:
        zeta_start = np.asarray(zeta_start, dtype=float).reshape(-1)
    if zeta_end is None:
        zeta_end = np.zeros((s,), dtype=float)
        zeta_end[0] = q[-1]
    else:
        zeta_end = np.asarray(zeta_end, dtype=float).reshape(-1)

    if zeta_start.shape != (s,):
        raise ValueError(f"zeta_start must have shape ({s},), got {zeta_start.shape}.")
    if zeta_end.shape != (s,):
        raise ValueError(f"zeta_end must have shape ({s},), got {zeta_end.shape}.")
    if not np.isclose(zeta_start[0], q[0]):
        raise ValueError("zeta_start[0] must equal q[0].")
    if not np.isclose(zeta_end[0], q[-1]):
        raise ValueError("zeta_end[0] must equal q[-1].")
    return q, T, zeta_start, zeta_end


def _window_bounds(i: int, k: int, M_seg: int) -> tuple[int, int]:
    if M_seg < 1:
        raise ValueError("M_seg must be >= 1.")
    if not 1 <= i <= M_seg:
        raise ValueError(f"window center i must be in [1, {M_seg}], got {i}.")
    if k < 0:
        raise ValueError("k must be >= 0.")
    half = k // 2
    left = max(1, i - half)
    right = min(M_seg, i + half)
    return left, right


def _nullspace_basis(G: np.ndarray, tol: float = 1e-12) -> tuple[np.ndarray, int]:
    G = np.asarray(G, dtype=float)
    if G.ndim != 2:
        raise ValueError("G must be a 2D array.")
    _, singular_values, vh = np.linalg.svd(G, full_matrices=True)
    rank = int(np.sum(singular_values > tol))
    basis = vh[rank:].T.copy()
    return basis, rank


def _stationarity_multiplier(H: np.ndarray, G: np.ndarray, c_vec: np.ndarray) -> np.ndarray:
    rhs = -(H @ c_vec)
    if G.size == 0:
        return np.zeros((0,), dtype=float)
    multiplier, *_ = np.linalg.lstsq(G.T, rhs, rcond=None)
    return multiplier


def _kkt_metrics(H: np.ndarray, G: np.ndarray, d_vec: np.ndarray, c_vec: np.ndarray) -> dict[str, float]:
    lam = _stationarity_multiplier(H, G, c_vec)
    stationarity = H @ c_vec + G.T @ lam
    feasibility = G @ c_vec - d_vec
    return {
        "stationarity_residual_norm": float(np.linalg.norm(stationarity)),
        "constraint_residual_norm": float(np.linalg.norm(feasibility)),
        "kkt_residual": float(np.linalg.norm(np.concatenate((stationarity, feasibility)))),
    }


def _snap_cost_block(duration: float, s: int = DEFAULT_S) -> np.ndarray:
    degree = 2 * s - 1
    block = np.zeros((degree + 1, degree + 1), dtype=float)
    for alpha in range(degree + 1):
        if alpha < s:
            continue
        coeff_alpha = float(math.factorial(alpha) / math.factorial(alpha - s))
        for beta_index in range(degree + 1):
            if beta_index < s:
                continue
            coeff_beta = float(math.factorial(beta_index) / math.factorial(beta_index - s))
            power = alpha + beta_index - 2 * s + 1
            block[alpha, beta_index] = coeff_alpha * coeff_beta * (float(duration) ** power) / power
    return block


def build_local_hessian(T_window: np.ndarray, s: int = DEFAULT_S) -> np.ndarray:
    """
    Build the local quadratic objective matrix H for

        min 0.5 * c^T H c

    corresponding to the BLOM-Strict energy integral of |p^(s)|^2.
    """
    if s < 2:
        raise ValueError("s must be >= 2.")
    T_window = _validate_time_vector(T_window)
    block_size = 2 * s
    H = np.zeros((block_size * T_window.size, block_size * T_window.size), dtype=float)
    for local_index, duration in enumerate(T_window):
        block = 2.0 * _snap_cost_block(float(duration), s=s)
        row_slice = slice(local_index * block_size, (local_index + 1) * block_size)
        H[row_slice, row_slice] = block
    return H


def build_local_constraints(
    q: np.ndarray,
    T: np.ndarray,
    i: int,
    k: int,
    zeta_start: np.ndarray | None = None,
    zeta_end: np.ndarray | None = None,
    s: int = DEFAULT_S,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build the affine local constraint system Gc = d.

    The row order is:
    1. segmentwise start interpolation
    2. segmentwise end interpolation
    3. interior derivative continuity for orders 1..s-1
    4. inherited physical start jet if the window touches the left boundary
    5. inherited physical end jet if the window touches the right boundary
    """
    q, T, zeta_start, zeta_end = _normalize_problem_inputs(q, T, zeta_start, zeta_end, s=s)
    left, right = _window_bounds(i, k, T.size)
    m_seg = right - left + 1
    block_size = 2 * s
    num_unknowns = block_size * m_seg

    num_rows = 2 * m_seg + (s - 1) * max(m_seg - 1, 0)
    if left == 1:
        num_rows += s - 1
    if right == T.size:
        num_rows += s - 1

    G = np.zeros((num_rows, num_unknowns), dtype=float)
    d_vec = np.zeros((num_rows,), dtype=float)
    row = 0

    for local_index, segment_index in enumerate(range(left, right + 1)):
        duration = float(T[segment_index - 1])
        coeff_slice = slice(local_index * block_size, (local_index + 1) * block_size)

        G[row, coeff_slice] = beta(0.0)
        d_vec[row] = float(q[segment_index - 1])
        row += 1

        G[row, coeff_slice] = beta(duration)
        d_vec[row] = float(q[segment_index])
        row += 1

    for local_index, segment_index in enumerate(range(left, right)):
        duration = float(T[segment_index - 1])
        left_slice = slice(local_index * block_size, (local_index + 1) * block_size)
        right_slice = slice((local_index + 1) * block_size, (local_index + 2) * block_size)
        for order in range(1, s):
            G[row, left_slice] = beta_d(duration, order)
            G[row, right_slice] = -beta_d(0.0, order)
            row += 1

    if left == 1:
        first_slice = slice(0, block_size)
        for order in range(1, s):
            G[row, first_slice] = beta_d(0.0, order)
            d_vec[row] = float(zeta_start[order])
            row += 1

    if right == T.size:
        last_slice = slice((m_seg - 1) * block_size, m_seg * block_size)
        duration = float(T[right - 1])
        for order in range(1, s):
            G[row, last_slice] = beta_d(duration, order)
            d_vec[row] = float(zeta_end[order])
            row += 1

    if row != num_rows:
        raise RuntimeError(f"Constraint assembly bug: filled {row} rows, expected {num_rows}.")
    return G, d_vec


def solve_kkt(H: np.ndarray, G: np.ndarray, d_vec: np.ndarray) -> dict[str, Any]:
    """Solve the equality-constrained QP through the saddle-point KKT system."""
    H = np.asarray(H, dtype=float)
    G = np.asarray(G, dtype=float)
    d_vec = np.asarray(d_vec, dtype=float).reshape(-1)
    if H.ndim != 2 or H.shape[0] != H.shape[1]:
        raise ValueError("H must be square.")
    if G.ndim != 2 or G.shape[1] != H.shape[0]:
        raise ValueError("G must have the same number of columns as H.")
    if d_vec.shape != (G.shape[0],):
        raise ValueError(f"d_vec must have shape ({G.shape[0]},), got {d_vec.shape}.")

    num_unknowns = H.shape[0]
    num_constraints = G.shape[0]
    K = np.zeros((num_unknowns + num_constraints, num_unknowns + num_constraints), dtype=float)
    K[:num_unknowns, :num_unknowns] = H
    K[:num_unknowns, num_unknowns:] = G.T
    K[num_unknowns:, :num_unknowns] = G

    rhs = np.concatenate((np.zeros((num_unknowns,), dtype=float), d_vec))
    solution = np.linalg.solve(K, rhs)
    c_vec = solution[:num_unknowns]
    multipliers = solution[num_unknowns:]

    metrics = _kkt_metrics(H, G, d_vec, c_vec)
    return {
        "c_vec": c_vec,
        "lambda": multipliers,
        "objective": float(0.5 * c_vec @ H @ c_vec),
        **metrics,
    }


def solve_reduced_qp(
    H: np.ndarray,
    G: np.ndarray,
    d_vec: np.ndarray,
    *,
    z0: np.ndarray | None = None,
    iterative: bool = False,
    max_iters: int = 2000,
    tol: float = 1e-12,
) -> dict[str, Any]:
    """
    Solve the local QP in null-space coordinates.

    When `iterative=False`, the reduced linear system is solved directly.
    When `iterative=True`, steepest descent is run in reduced coordinates to
    support multistart uniqueness experiments.
    """
    H = np.asarray(H, dtype=float)
    G = np.asarray(G, dtype=float)
    d_vec = np.asarray(d_vec, dtype=float).reshape(-1)

    c_particular, *_ = np.linalg.lstsq(G, d_vec, rcond=None)
    particular_residual = float(np.linalg.norm(G @ c_particular - d_vec))
    null_basis, rank = _nullspace_basis(G)
    null_dim = null_basis.shape[1]

    if null_dim == 0:
        c_vec = c_particular
        reduced_hessian = np.zeros((0, 0), dtype=float)
        z_star = np.zeros((0,), dtype=float)
        iterations = 0
        converged = True
        min_eig = float("inf")
        cond_number = 1.0
    else:
        reduced_hessian = null_basis.T @ H @ null_basis
        reduced_gradient = null_basis.T @ H @ c_particular
        eigvals = np.linalg.eigvalsh(reduced_hessian)
        min_eig = float(np.min(eigvals))
        cond_number = float(np.linalg.cond(reduced_hessian))

        if z0 is None:
            z = np.zeros((null_dim,), dtype=float)
        else:
            z = np.asarray(z0, dtype=float).reshape(-1)
            if z.shape != (null_dim,):
                raise ValueError(f"z0 must have shape ({null_dim},), got {z.shape}.")

        if iterative:
            rhs = -reduced_gradient
            residual = rhs - reduced_hessian @ z
            direction = residual.copy()
            converged = False
            iterations = 0
            residual_norm_sq = float(residual @ residual)
            if residual_norm_sq <= tol * tol:
                converged = True
            for iterations in range(1, max_iters + 1):
                if converged:
                    break
                matvec = reduced_hessian @ direction
                curvature = float(direction @ matvec)
                if curvature <= 0.0:
                    raise np.linalg.LinAlgError("Reduced Hessian is not positive definite.")
                step = residual_norm_sq / curvature
                z = z + step * direction
                new_residual = residual - step * matvec
                new_norm_sq = float(new_residual @ new_residual)
                if new_norm_sq <= tol * tol:
                    residual = new_residual
                    converged = True
                    break
                beta_coeff = new_norm_sq / residual_norm_sq
                direction = new_residual + beta_coeff * direction
                residual = new_residual
                residual_norm_sq = new_norm_sq
            z_star = z
            if not converged:
                gradient = reduced_hessian @ z_star + reduced_gradient
                converged = float(np.linalg.norm(gradient)) <= 10.0 * tol
        else:
            z_star = -np.linalg.solve(reduced_hessian, reduced_gradient)
            iterations = 1
            converged = True

        c_vec = c_particular + null_basis @ z_star

    metrics = _kkt_metrics(H, G, d_vec, c_vec)
    return {
        "c_vec": c_vec,
        "objective": float(0.5 * c_vec @ H @ c_vec),
        "c_particular": c_particular,
        "particular_residual_norm": particular_residual,
        "nullspace_basis": null_basis,
        "constraint_rank": rank,
        "reduced_hessian": reduced_hessian,
        "reduced_hessian_min_eig": min_eig,
        "reduced_hessian_cond": cond_number,
        "z_star": z_star,
        "iterations": iterations,
        "converged": converged,
        **metrics,
    }
