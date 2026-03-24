"""
Piecewise trajectory evaluation for Phase 0 degree-7 trajectories.
"""

from __future__ import annotations

import numpy as np

from phase_0.poly_basis import eval_poly


def locate_segment(T: np.ndarray, t: float) -> tuple[int, float]:
    """Locate the segment index and local time for one global time value."""
    T = np.asarray(T, dtype=float).reshape(-1)
    total_time = float(np.sum(T))
    if not (0.0 <= t <= total_time):
        raise ValueError(f"t must lie in [0, {total_time}], got {t}.")
    cumulative = np.cumsum(T)
    if np.isclose(t, total_time):
        return T.size - 1, float(T[-1])
    segment = int(np.searchsorted(cumulative, t, side="right"))
    start = 0.0 if segment == 0 else float(cumulative[segment - 1])
    return segment, float(t - start)


def eval_piecewise(
    coeffs: np.ndarray,
    T: np.ndarray,
    t: float | np.ndarray,
    order: int = 0,
) -> float | np.ndarray:
    """Evaluate a piecewise polynomial trajectory or one derivative order."""
    coeffs = np.asarray(coeffs, dtype=float)
    T = np.asarray(T, dtype=float).reshape(-1)
    t_array = np.asarray(t, dtype=float)

    if t_array.ndim == 0:
        segment, local_t = locate_segment(T, float(t_array))
        return eval_poly(coeffs[segment], local_t, order=order)

    values = np.empty_like(t_array, dtype=float)
    for index, time_value in np.ndenumerate(t_array):
        segment, local_t = locate_segment(T, float(time_value))
        values[index] = eval_poly(coeffs[segment], local_t, order=order)
    return values


def sample_trajectory(
    coeffs: np.ndarray,
    T: np.ndarray,
    num: int = 200,
    order: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample a trajectory on an evenly spaced global time grid."""
    T = np.asarray(T, dtype=float).reshape(-1)
    total_time = float(np.sum(T))
    times = np.linspace(0.0, total_time, num=num)
    values = eval_piecewise(coeffs, T, times, order=order)
    return times, np.asarray(values, dtype=float)


def junction_jumps(
    coeffs: np.ndarray,
    T: np.ndarray,
    max_order: int = 6,
) -> np.ndarray:
    """
    Return jumps across internal junctions.

    Output shape is (max_order + 1, M - 1), where entry [l, i] equals
    p_i^(l)(T_i) - p_{i+1}^(l)(0).
    """
    coeffs = np.asarray(coeffs, dtype=float)
    T = np.asarray(T, dtype=float).reshape(-1)
    num_segments = T.size
    jumps = np.zeros((max_order + 1, max(num_segments - 1, 0)), dtype=float)
    for segment in range(num_segments - 1):
        for order in range(max_order + 1):
            left = eval_poly(coeffs[segment], float(T[segment]), order=order)
            right = eval_poly(coeffs[segment + 1], 0.0, order=order)
            jumps[order, segment] = left - right
    return jumps
