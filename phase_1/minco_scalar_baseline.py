"""
Phase 1 scalar MINCO baseline with fixed boundary jets.

This module implements the strengthened global mother problem described in
`REQ_phase_1_minco_scalar_baseline.md`: the trajectory is scalar, each segment
is degree 7, interior waypoints constrain only position, and the boundary jets
up to order 3 are fixed at the start and end.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from phase_0.poly_basis import monomial_row, derivative_row


PHASE1_DEGREE = 7
PHASE1_JET_ORDER = 3
PHASE1_CONTINUITY_ORDER = 6
COEFFS_PER_SEGMENT = PHASE1_DEGREE + 1


def beta(t: float) -> np.ndarray:
    """Return the monomial basis beta(t) = [1, t, ..., t^7]."""
    return monomial_row(float(t), degree=PHASE1_DEGREE)


def beta_d(t: float, order: int) -> np.ndarray:
    """Return the derivative basis beta^(order)(t)."""
    if not 0 <= order <= PHASE1_DEGREE:
        raise ValueError(f"order must be in [0, {PHASE1_DEGREE}], got {order}.")
    return derivative_row(float(t), degree=PHASE1_DEGREE, order=order)


def _validate_time_vector(T: np.ndarray) -> np.ndarray:
    T = np.asarray(T, dtype=float).reshape(-1)
    if T.size < 1:
        raise ValueError("Phase 1 requires at least one segment, so len(T) must be >= 1.")
    if not np.all(np.isfinite(T)):
        raise ValueError("T contains non-finite values.")
    if not np.all(T > 0.0):
        raise ValueError("All segment durations must be strictly positive.")
    return T


def _validate_inputs(
    q: np.ndarray,
    T: np.ndarray,
    zeta_start: np.ndarray,
    zeta_end: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Validate the Phase 1 MINCO baseline inputs and return normalized arrays."""
    q = np.asarray(q, dtype=float).reshape(-1)
    T = _validate_time_vector(T)
    zeta_start = np.asarray(zeta_start, dtype=float).reshape(-1)
    zeta_end = np.asarray(zeta_end, dtype=float).reshape(-1)

    if q.shape != (T.size + 1,):
        raise ValueError(f"q must have shape ({T.size + 1},), got {q.shape}.")
    if zeta_start.shape != (PHASE1_JET_ORDER + 1,):
        raise ValueError(
            f"zeta_start must have shape ({PHASE1_JET_ORDER + 1},), got {zeta_start.shape}."
        )
    if zeta_end.shape != (PHASE1_JET_ORDER + 1,):
        raise ValueError(
            f"zeta_end must have shape ({PHASE1_JET_ORDER + 1},), got {zeta_end.shape}."
        )
    if not np.all(np.isfinite(q)):
        raise ValueError("q contains non-finite values.")
    if not np.all(np.isfinite(zeta_start)):
        raise ValueError("zeta_start contains non-finite values.")
    if not np.all(np.isfinite(zeta_end)):
        raise ValueError("zeta_end contains non-finite values.")
    if not np.isclose(zeta_start[0], q[0]):
        raise ValueError("zeta_start[0] must equal q[0].")
    if not np.isclose(zeta_end[0], q[-1]):
        raise ValueError("zeta_end[0] must equal q[-1].")
    return q, T, zeta_start, zeta_end


def _coefficient_slice(segment_index: int) -> slice:
    start = segment_index * COEFFS_PER_SEGMENT
    stop = start + COEFFS_PER_SEGMENT
    return slice(start, stop)


def build_system_matrix(T: np.ndarray) -> np.ndarray:
    """
    Assemble the square linear system M(T) for the strengthened MINCO baseline.

    The row order is:
    1. segment start interpolation p_i(0) = q_{i-1}
    2. segment end interpolation p_i(T_i) = q_i
    3. interior continuity for derivative orders 1..6
    4. initial boundary jet constraints for derivative orders 1..3
    5. terminal boundary jet constraints for derivative orders 1..3
    """
    T = _validate_time_vector(T)
    M_seg = T.size
    num_unknowns = COEFFS_PER_SEGMENT * M_seg
    M_mat = np.zeros((num_unknowns, num_unknowns), dtype=float)

    row = 0

    for segment_index, duration in enumerate(T):
        segment_slice = _coefficient_slice(segment_index)

        # Row block for endpoint interpolation p_i(0) = q_{i-1}.
        M_mat[row, segment_slice] = beta(0.0)
        row += 1

        # Row block for endpoint interpolation p_i(T_i) = q_i.
        M_mat[row, segment_slice] = beta(float(duration))
        row += 1

    for segment_index, duration in enumerate(T[:-1]):
        left_slice = _coefficient_slice(segment_index)
        right_slice = _coefficient_slice(segment_index + 1)
        for order in range(1, PHASE1_CONTINUITY_ORDER + 1):
            # Continuity of derivative order `order` at the interior knot.
            M_mat[row, left_slice] = beta_d(float(duration), order)
            M_mat[row, right_slice] = -beta_d(0.0, order)
            row += 1

    for order in range(1, PHASE1_JET_ORDER + 1):
        # Row block for initial jet constraints p^(order)(0) = zeta^-_order.
        M_mat[row, _coefficient_slice(0)] = beta_d(0.0, order)
        row += 1

    for order in range(1, PHASE1_JET_ORDER + 1):
        # Row block for terminal jet constraints p^(order)(T_total) = zeta^+_order.
        M_mat[row, _coefficient_slice(M_seg - 1)] = beta_d(float(T[-1]), order)
        row += 1

    if row != num_unknowns:
        raise RuntimeError(
            f"System assembly bug: filled {row} rows for a {num_unknowns}x{num_unknowns} matrix."
        )
    return M_mat


def build_rhs(
    q: np.ndarray,
    zeta_start: np.ndarray,
    zeta_end: np.ndarray,
) -> np.ndarray:
    """
    Assemble b(q, zeta_start, zeta_end) using the same row order as `build_system_matrix`.
    """
    q = np.asarray(q, dtype=float).reshape(-1)
    zeta_start = np.asarray(zeta_start, dtype=float).reshape(-1)
    zeta_end = np.asarray(zeta_end, dtype=float).reshape(-1)

    M_seg = q.size - 1
    if M_seg < 1:
        raise ValueError("q must describe at least one segment.")

    rhs: list[float] = []

    for segment_index in range(M_seg):
        rhs.append(float(q[segment_index]))
        rhs.append(float(q[segment_index + 1]))

    rhs.extend([0.0] * (PHASE1_CONTINUITY_ORDER * max(M_seg - 1, 0)))
    rhs.extend(float(zeta_start[order]) for order in range(1, PHASE1_JET_ORDER + 1))
    rhs.extend(float(zeta_end[order]) for order in range(1, PHASE1_JET_ORDER + 1))
    return np.asarray(rhs, dtype=float)


def solve_minco_coefficients(
    q: np.ndarray,
    T: np.ndarray,
    zeta_start: np.ndarray,
    zeta_end: np.ndarray,
    *,
    return_system: bool = False,
) -> dict[str, Any]:
    """
    Solve M(T) c = b(q, zeta_start, zeta_end) and return the unique coefficients.

    The system is solved with `numpy.linalg.solve` in Phase 1 for transparency.
    The matrix itself is assembled separately so a future banded solver can
    replace this dense solve without changing the public API.
    """
    q, T, zeta_start, zeta_end = _validate_inputs(q, T, zeta_start, zeta_end)
    M_mat = build_system_matrix(T)
    b_vec = build_rhs(q, zeta_start, zeta_end)

    try:
        c_vec = np.linalg.solve(M_mat, b_vec)
    except np.linalg.LinAlgError as exc:
        raise np.linalg.LinAlgError(
            "Failed to solve the Phase 1 MINCO baseline system; "
            "the assembled matrix may be singular for the provided inputs."
        ) from exc

    coeffs = c_vec.reshape(T.size, COEFFS_PER_SEGMENT)
    result: dict[str, Any] = {
        "coeffs": coeffs,
        "c_vec": c_vec,
        "residual_norm": float(np.linalg.norm(M_mat @ c_vec - b_vec)),
    }
    if return_system:
        result["M"] = M_mat
        result["b"] = b_vec
    return result


def evaluate_segment(coeff: np.ndarray, t: float, order: int = 0) -> float:
    """Evaluate one segment derivative p_i^(order)(t)."""
    coeff = np.asarray(coeff, dtype=float).reshape(-1)
    if coeff.shape != (COEFFS_PER_SEGMENT,):
        raise ValueError(f"coeff must have shape ({COEFFS_PER_SEGMENT},), got {coeff.shape}.")
    if not np.isfinite(t):
        raise ValueError("t must be finite.")
    return float(coeff @ beta_d(float(t), order))


def _locate_segment(T: np.ndarray, t_global: float) -> tuple[int, float]:
    T = _validate_time_vector(T)
    T_total = float(np.sum(T))
    if not (0.0 <= t_global <= T_total):
        raise ValueError(f"t_global must lie in [0, {T_total}], got {t_global}.")
    cumulative = np.cumsum(T)
    if np.isclose(t_global, T_total):
        return T.size - 1, float(T[-1])
    segment_index = int(np.searchsorted(cumulative, t_global, side="right"))
    start_time = 0.0 if segment_index == 0 else float(cumulative[segment_index - 1])
    return segment_index, float(t_global - start_time)


def evaluate_trajectory(
    coeffs: np.ndarray,
    T: np.ndarray,
    t_global: float,
    order: int = 0,
) -> float:
    """Evaluate the global trajectory at one global time."""
    coeffs = np.asarray(coeffs, dtype=float)
    T = _validate_time_vector(T)
    if coeffs.shape != (T.size, COEFFS_PER_SEGMENT):
        raise ValueError(f"coeffs must have shape ({T.size}, {COEFFS_PER_SEGMENT}), got {coeffs.shape}.")
    segment_index, local_time = _locate_segment(T, float(t_global))
    return evaluate_segment(coeffs[segment_index], local_time, order=order)


def sample_trajectory(
    coeffs: np.ndarray,
    T: np.ndarray,
    num_per_segment: int = 50,
    orders: tuple[int, ...] = (0, 1, 2, 3, 4),
) -> dict[str, np.ndarray]:
    """
    Sample the global trajectory on a per-segment grid for the requested derivative orders.
    """
    coeffs = np.asarray(coeffs, dtype=float)
    T = _validate_time_vector(T)
    if coeffs.shape != (T.size, COEFFS_PER_SEGMENT):
        raise ValueError(f"coeffs must have shape ({T.size}, {COEFFS_PER_SEGMENT}), got {coeffs.shape}.")
    if num_per_segment < 2:
        raise ValueError("num_per_segment must be >= 2.")

    times: list[float] = []
    cursor = 0.0
    for segment_index, duration in enumerate(T):
        local_times = np.linspace(
            0.0,
            float(duration),
            num=num_per_segment,
            endpoint=(segment_index == T.size - 1),
        )
        times.extend((cursor + local_times).tolist())
        cursor += float(duration)

    t_global = np.asarray(times, dtype=float)
    samples: dict[str, np.ndarray] = {"t_global": t_global}
    for order in orders:
        samples[f"order_{order}"] = np.asarray(
            [evaluate_trajectory(coeffs, T, time_value, order=order) for time_value in t_global],
            dtype=float,
        )
    return samples


def interpolation_errors(coeffs: np.ndarray, T: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Return start/end interpolation errors for every segment."""
    coeffs = np.asarray(coeffs, dtype=float)
    T = _validate_time_vector(T)
    q = np.asarray(q, dtype=float).reshape(-1)
    errors = []
    for segment_index, duration in enumerate(T):
        errors.append(evaluate_segment(coeffs[segment_index], 0.0, order=0) - q[segment_index])
        errors.append(
            evaluate_segment(coeffs[segment_index], float(duration), order=0) - q[segment_index + 1]
        )
    return np.asarray(errors, dtype=float)


def continuity_jumps(coeffs: np.ndarray, T: np.ndarray, max_order: int = PHASE1_CONTINUITY_ORDER) -> np.ndarray:
    """Return interior continuity jumps for derivative orders 0..max_order."""
    coeffs = np.asarray(coeffs, dtype=float)
    T = _validate_time_vector(T)
    jumps = np.zeros((max_order + 1, max(T.size - 1, 0)), dtype=float)
    for segment_index, duration in enumerate(T[:-1]):
        for order in range(max_order + 1):
            left_value = evaluate_segment(coeffs[segment_index], float(duration), order=order)
            right_value = evaluate_segment(coeffs[segment_index + 1], 0.0, order=order)
            jumps[order, segment_index] = left_value - right_value
    return jumps


def boundary_jet_errors(
    coeffs: np.ndarray,
    T: np.ndarray,
    zeta_start: np.ndarray,
    zeta_end: np.ndarray,
) -> np.ndarray:
    """Return start/end jet errors for derivative orders 0..3."""
    coeffs = np.asarray(coeffs, dtype=float)
    T = _validate_time_vector(T)
    zeta_start = np.asarray(zeta_start, dtype=float).reshape(-1)
    zeta_end = np.asarray(zeta_end, dtype=float).reshape(-1)

    errors = []
    for order in range(PHASE1_JET_ORDER + 1):
        errors.append(evaluate_segment(coeffs[0], 0.0, order=order) - zeta_start[order])
        errors.append(evaluate_segment(coeffs[-1], float(T[-1]), order=order) - zeta_end[order])
    return np.asarray(errors, dtype=float)


def system_residual(M_mat: np.ndarray, c_vec: np.ndarray, b_vec: np.ndarray) -> float:
    """Return ||M c - b||_2 for the assembled Phase 1 system."""
    M_mat = np.asarray(M_mat, dtype=float)
    c_vec = np.asarray(c_vec, dtype=float).reshape(-1)
    b_vec = np.asarray(b_vec, dtype=float).reshape(-1)
    return float(np.linalg.norm(M_mat @ c_vec - b_vec))


def solve_adjoint(M_mat: np.ndarray, dK_dc: np.ndarray) -> np.ndarray:
    """
    Solve M(T)^T lambda = dK/dc for future gradient propagation.

    Phase 1 only needs the interface and linear solve entry point; the actual
    high-level gradient formulas are intentionally deferred to later phases.
    """
    M_mat = np.asarray(M_mat, dtype=float)
    dK_dc = np.asarray(dK_dc, dtype=float).reshape(-1)
    return np.linalg.solve(M_mat.T, dK_dc)


def grad_wrt_q(*_args: Any, **_kwargs: Any) -> np.ndarray:
    """
    Placeholder for dW/dq.

    TODO: derive and implement the explicit Phase 1 gradient propagation rules
    that reuse the MINCO system solve instead of forming any dense inverse.
    """
    raise NotImplementedError("Phase 1 reserves this interface for later gradient propagation.")


def grad_wrt_T(*_args: Any, **_kwargs: Any) -> np.ndarray:
    """
    Placeholder for dW/dT.

    TODO: derive and implement the explicit time-sensitivity formulas for the
    strengthened MINCO mother problem.
    """
    raise NotImplementedError("Phase 1 reserves this interface for later gradient propagation.")

