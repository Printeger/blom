"""
Constraint assembly for the canonical Phase 0 minimum-snap problems.

This module keeps the linear algebra visible: each row corresponds to a
mathematical condition on degree-7 segment coefficients.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from phase_0.blom_problem_setup import BLOMProblemSetup
from phase_0.poly_basis import COEFFS_PER_SEGMENT, PHASE0_DEGREE, derivative_row, snap_cost_matrix


CONTINUITY_ORDERS = tuple(range(1, 2 * 4 - 1))
NATURAL_BOUNDARY_ORDERS = tuple(range(4, 2 * 4 - 1))


@dataclass
class ConstraintSystem:
    """Linear equality system A x = b for either the global or local problem."""

    A: np.ndarray
    b: np.ndarray
    segment_indices: tuple[int, ...]
    durations: np.ndarray

    def residual(self, x: np.ndarray) -> np.ndarray:
        return self.A @ x - self.b


def coefficient_slice(local_segment_index: int) -> slice:
    """Return the slice for one segment inside the flattened coefficient vector."""
    start = local_segment_index * COEFFS_PER_SEGMENT
    stop = start + COEFFS_PER_SEGMENT
    return slice(start, stop)


def build_snap_hessian(durations: np.ndarray) -> np.ndarray:
    """
    Build H so that 0.5 * x.T @ H @ x is the summed minimum-snap cost.
    """
    durations = np.asarray(durations, dtype=float).reshape(-1)
    size = durations.size * COEFFS_PER_SEGMENT
    H = np.zeros((size, size), dtype=float)
    for local_index, duration in enumerate(durations):
        sl = coefficient_slice(local_index)
        H[sl, sl] = 2.0 * snap_cost_matrix(float(duration), degree=PHASE0_DEGREE)
    return H


def build_global_minco_constraints(setup: BLOMProblemSetup) -> ConstraintSystem:
    """Build the canonical global Phase 0 spline system."""
    setup.validate()
    return _build_constraint_system(setup, tuple(range(setup.M)))


def build_local_blom_constraints(
    setup: BLOMProblemSetup,
    center_segment: int,
) -> ConstraintSystem:
    """Build the canonical local BLOM window system for W(i, 2)."""
    setup.validate()
    return _build_constraint_system(setup, setup.window(center_segment))


def _build_constraint_system(
    setup: BLOMProblemSetup,
    segment_indices: tuple[int, ...],
) -> ConstraintSystem:
    rows: list[np.ndarray] = []
    rhs: list[float] = []
    num_segments = len(segment_indices)
    total_coeffs = num_segments * COEFFS_PER_SEGMENT

    def add_row(local_segment_index: int, time_value: float, order: int, value: float) -> None:
        row = np.zeros(total_coeffs, dtype=float)
        row[coefficient_slice(local_segment_index)] = derivative_row(
            time_value,
            degree=PHASE0_DEGREE,
            order=order,
        )
        rows.append(row)
        rhs.append(float(value))

    for local_index, global_index in enumerate(segment_indices):
        duration = float(setup.T[global_index])
        add_row(local_index, 0.0, order=0, value=float(setup.q[global_index]))
        add_row(local_index, duration, order=0, value=float(setup.q[global_index + 1]))

    for local_index in range(num_segments - 1):
        left_global_index = segment_indices[local_index]
        duration = float(setup.T[left_global_index])
        for order in CONTINUITY_ORDERS:
            row = np.zeros(total_coeffs, dtype=float)
            row[coefficient_slice(local_index)] = derivative_row(
                duration,
                degree=PHASE0_DEGREE,
                order=order,
            )
            row[coefficient_slice(local_index + 1)] = -derivative_row(
                0.0,
                degree=PHASE0_DEGREE,
                order=order,
            )
            rows.append(row)
            rhs.append(0.0)

    first_local = 0
    last_local = num_segments - 1
    last_duration = float(setup.T[segment_indices[-1]])
    for order in NATURAL_BOUNDARY_ORDERS:
        add_row(first_local, 0.0, order=order, value=0.0)
        add_row(last_local, last_duration, order=order, value=0.0)

    return ConstraintSystem(
        A=np.vstack(rows),
        b=np.asarray(rhs, dtype=float),
        segment_indices=segment_indices,
        durations=setup.T[list(segment_indices)].copy(),
    )
