"""
Local BLOM window solver for the canonical Phase 0 setting (s=4, k=2).

The result is intentionally local: it solves only the window problem around
one center segment and does not claim that independent windows glue into a
globally correct BLOM trajectory.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from phase_0.blom_problem_setup import BLOMProblemSetup
from phase_0.constraint_builders import (
    ConstraintSystem,
    build_local_blom_constraints,
    build_snap_hessian,
)
from phase_0.minco_scalar_baseline import solve_equality_constrained_qp


@dataclass
class BlomLocalResult:
    center_segment: int
    active_segments: tuple[int, ...]
    window_coeffs: np.ndarray
    center_coeffs: np.ndarray
    cost: float
    setup: BLOMProblemSetup
    constraint_system: ConstraintSystem


def solve_blom_local_qp(
    setup: BLOMProblemSetup,
    center_segment: int,
) -> BlomLocalResult:
    """Solve the canonical local window QP around one center segment."""
    setup.validate()
    constraints = build_local_blom_constraints(setup, center_segment=center_segment)
    H = build_snap_hessian(constraints.durations)
    x = solve_equality_constrained_qp(H, constraints.A, constraints.b)
    window_coeffs = x.reshape(len(constraints.segment_indices), -1)
    center_local_index = constraints.segment_indices.index(center_segment)
    center_coeffs = window_coeffs[center_local_index].copy()
    cost = 0.5 * float(x @ H @ x)
    return BlomLocalResult(
        center_segment=center_segment,
        active_segments=constraints.segment_indices,
        window_coeffs=window_coeffs,
        center_coeffs=center_coeffs,
        cost=cost,
        setup=setup,
        constraint_system=constraints,
    )
