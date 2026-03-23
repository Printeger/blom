"""
Global scalar Phase 0 minimum-snap baseline.

This is a numerical reference solver for the canonical 1D problem. It uses
the shared constraint builder and a dense equality-constrained QP solve.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from phase_0.blom_problem_setup import BLOMProblemSetup
from phase_0.constraint_builders import (
    ConstraintSystem,
    build_global_minco_constraints,
    build_snap_hessian,
)


@dataclass
class MincoScalarResult:
    coeffs: np.ndarray
    cost: float
    setup: BLOMProblemSetup
    constraint_system: ConstraintSystem


def solve_equality_constrained_qp(
    H: np.ndarray,
    A: np.ndarray,
    b: np.ndarray,
) -> np.ndarray:
    """
    Solve min 0.5 x^T H x subject to A x = b.

    H is positive semidefinite; the Phase 0 constraints make the solution unique.
    """
    H = np.asarray(H, dtype=float)
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float).reshape(-1)
    num_variables = H.shape[0]
    num_constraints = A.shape[0]
    kkt = np.zeros((num_variables + num_constraints, num_variables + num_constraints))
    kkt[:num_variables, :num_variables] = H
    kkt[:num_variables, num_variables:] = A.T
    kkt[num_variables:, :num_variables] = A
    rhs = np.concatenate((np.zeros(num_variables, dtype=float), b))

    try:
        solution = np.linalg.solve(kkt, rhs)
    except np.linalg.LinAlgError:
        solution, *_ = np.linalg.lstsq(kkt, rhs, rcond=None)
    return solution[:num_variables]


def solve_minco_scalar(setup: BLOMProblemSetup) -> MincoScalarResult:
    """Solve the global scalar Phase 0 minimum-snap reference problem."""
    setup.validate()
    constraints = build_global_minco_constraints(setup)
    H = build_snap_hessian(constraints.durations)
    x = solve_equality_constrained_qp(H, constraints.A, constraints.b)
    cost = 0.5 * float(x @ H @ x)
    coeffs = x.reshape(setup.M, -1)
    return MincoScalarResult(
        coeffs=coeffs,
        cost=cost,
        setup=setup,
        constraint_system=constraints,
    )
