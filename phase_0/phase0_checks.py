"""
Centralized validation utilities for the Phase 0 BLOM code.
"""

from __future__ import annotations

import numpy as np

from phase_0.blom_local_qp import BlomLocalResult
from phase_0.blom_problem_setup import BLOMProblemSetup
from phase_0.minco_scalar_baseline import MincoScalarResult
from phase_0.poly_basis import eval_poly
from phase_0.trajectory_eval import junction_jumps


def _max_abs(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return 0.0 if values.size == 0 else float(np.max(np.abs(values)))


def check_setup(setup: BLOMProblemSetup) -> dict[str, float | bool]:
    """Validate canonical setup data and return a small report."""
    setup.validate()
    return {
        "passed": True,
        "M": setup.M,
        "theta_size": int(setup.theta().size),
        "total_time": setup.total_time,
        "min_duration": float(np.min(setup.T)),
    }


def check_minco_result(
    result: MincoScalarResult,
    tol: float = 1e-8,
) -> dict[str, float | bool]:
    """Check interpolation, continuity, boundary residuals, and cost."""
    setup = result.setup
    setup.validate()
    coeffs = np.asarray(result.coeffs, dtype=float)

    interpolation_errors = []
    for segment in range(setup.M):
        interpolation_errors.append(eval_poly(coeffs[segment], 0.0, order=0) - setup.q[segment])
        interpolation_errors.append(
            eval_poly(coeffs[segment], float(setup.T[segment]), order=0) - setup.q[segment + 1]
        )

    jumps = junction_jumps(coeffs, setup.T, max_order=6)
    boundary_residuals = []
    for order in range(4, 7):
        boundary_residuals.append(eval_poly(coeffs[0], 0.0, order=order))
        boundary_residuals.append(
            eval_poly(coeffs[-1], float(setup.T[-1]), order=order)
        )

    max_interp = _max_abs(np.asarray(interpolation_errors))
    max_jump = _max_abs(jumps)
    max_bc = _max_abs(np.asarray(boundary_residuals))
    return {
        "passed": (
            max_interp < tol
            and max_jump < tol
            and max_bc < tol
            and np.isfinite(result.cost)
            and result.cost >= -tol
        ),
        "max_interpolation_error": max_interp,
        "max_junction_jump": max_jump,
        "max_boundary_residual": max_bc,
        "cost": float(result.cost),
    }


def check_blom_local_result(
    result: BlomLocalResult,
    tol: float = 1e-8,
) -> dict[str, float | bool]:
    """Check local interpolation, local continuity, and natural BC residuals."""
    setup = result.setup
    setup.validate()
    coeffs = np.asarray(result.window_coeffs, dtype=float)
    active_segments = result.active_segments

    interpolation_errors = []
    for local_index, global_index in enumerate(active_segments):
        interpolation_errors.append(
            eval_poly(coeffs[local_index], 0.0, order=0) - setup.q[global_index]
        )
        interpolation_errors.append(
            eval_poly(coeffs[local_index], float(setup.T[global_index]), order=0)
            - setup.q[global_index + 1]
        )

    continuity_errors = []
    for local_index in range(len(active_segments) - 1):
        duration = float(setup.T[active_segments[local_index]])
        for order in range(7):
            continuity_errors.append(
                eval_poly(coeffs[local_index], duration, order=order)
                - eval_poly(coeffs[local_index + 1], 0.0, order=order)
            )

    boundary_residuals = []
    right_duration = float(setup.T[active_segments[-1]])
    for order in range(4, 7):
        boundary_residuals.append(eval_poly(coeffs[0], 0.0, order=order))
        boundary_residuals.append(eval_poly(coeffs[-1], right_duration, order=order))

    max_interp = _max_abs(np.asarray(interpolation_errors))
    max_continuity = _max_abs(np.asarray(continuity_errors))
    max_bc = _max_abs(np.asarray(boundary_residuals))
    return {
        "passed": (
            max_interp < tol
            and max_continuity < tol
            and max_bc < tol
            and np.isfinite(result.cost)
            and result.cost >= -tol
        ),
        "max_interpolation_error": max_interp,
        "max_continuity_error": max_continuity,
        "max_boundary_residual": max_bc,
        "cost": float(result.cost),
    }
