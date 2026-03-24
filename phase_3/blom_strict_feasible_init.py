"""
Explicit feasible initialization for the Phase 3 BLOM-Strict local problem.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from phase_1.minco_scalar_baseline import beta_d
from phase_3.blom_strict_local_kkt import DEFAULT_S, build_local_constraints


def _window_bounds(i: int, k: int, M_seg: int) -> tuple[int, int]:
    half = k // 2
    left = max(1, i - half)
    right = min(M_seg, i + half)
    return left, right


def _normalize_boundary_jets(
    q: np.ndarray,
    zeta_start: np.ndarray | None,
    zeta_end: np.ndarray | None,
    s: int,
) -> tuple[np.ndarray, np.ndarray]:
    if zeta_start is None:
        zeta_start = np.zeros((s,), dtype=float)
        zeta_start[0] = float(q[0])
    else:
        zeta_start = np.asarray(zeta_start, dtype=float).reshape(-1)
    if zeta_end is None:
        zeta_end = np.zeros((s,), dtype=float)
        zeta_end[0] = float(q[-1])
    else:
        zeta_end = np.asarray(zeta_end, dtype=float).reshape(-1)
    return zeta_start, zeta_end


def _segment_from_endpoint_jets(
    duration: float,
    left_jet: np.ndarray,
    right_jet: np.ndarray,
    s: int,
) -> np.ndarray:
    block_size = 2 * s
    hermite = np.zeros((block_size, block_size), dtype=float)
    rhs = np.zeros((block_size,), dtype=float)

    for order in range(s):
        hermite[order] = beta_d(0.0, order)
        rhs[order] = float(left_jet[order])

    for order in range(s):
        hermite[s + order] = beta_d(float(duration), order)
        rhs[s + order] = float(right_jet[order])

    return np.linalg.solve(hermite, rhs)


def build_feasible_local_spline(
    q: np.ndarray,
    T: np.ndarray,
    i: int,
    k: int,
    zeta_start: np.ndarray | None = None,
    zeta_end: np.ndarray | None = None,
    s: int = DEFAULT_S,
) -> dict[str, Any]:
    """
    Construct a feasible local spline by assigning jets at all local knots and
    solving one Hermite interpolation problem per segment.
    """
    q = np.asarray(q, dtype=float).reshape(-1)
    T = np.asarray(T, dtype=float).reshape(-1)
    if q.shape != (T.size + 1,):
        raise ValueError(f"q must have shape ({T.size + 1},), got {q.shape}.")
    zeta_start, zeta_end = _normalize_boundary_jets(q, zeta_start, zeta_end, s=s)
    left, right = _window_bounds(i, k, T.size)

    knot_jets: dict[int, np.ndarray] = {}
    for knot_index in range(left - 1, right + 1):
        jet = np.zeros((s,), dtype=float)
        jet[0] = float(q[knot_index])
        if knot_index == 0 and left == 1:
            jet[:] = zeta_start
        elif knot_index == T.size and right == T.size:
            jet[:] = zeta_end
        knot_jets[knot_index] = jet

    coeffs = []
    for segment_index in range(left, right + 1):
        duration = float(T[segment_index - 1])
        left_jet = knot_jets[segment_index - 1]
        right_jet = knot_jets[segment_index]
        coeff = _segment_from_endpoint_jets(duration, left_jet, right_jet, s=s)
        coeffs.append(coeff)

    coeffs_arr = np.asarray(coeffs, dtype=float)
    c_loc = coeffs_arr.reshape(-1)
    G, d_vec = build_local_constraints(q, T, i, k, zeta_start=zeta_start, zeta_end=zeta_end, s=s)
    residual = float(np.linalg.norm(G @ c_loc - d_vec))
    return {
        "coeffs": coeffs_arr,
        "c_loc": c_loc,
        "knot_jets": {key: value.copy() for key, value in knot_jets.items()},
        "constraint_residual_norm": residual,
    }

