"""
Shared 1D polynomial utilities for Phase 0.

Each segment uses a degree-7 monomial basis because Phase 0 fixes s = 4,
so 2s - 1 = 7.
"""

from __future__ import annotations

import numpy as np


PHASE0_DEGREE = 7
PHASE0_S = 4
COEFFS_PER_SEGMENT = PHASE0_DEGREE + 1


def _falling_factorial(n: int, order: int) -> float:
    if order < 0:
        raise ValueError("order must be non-negative.")
    if order > n:
        return 0.0
    value = 1.0
    for k in range(order):
        value *= n - k
    return value


def monomial_row(t: float, degree: int = PHASE0_DEGREE) -> np.ndarray:
    """Return [1, t, ..., t^degree]."""
    return np.array([t**power for power in range(degree + 1)], dtype=float)


def derivative_row(
    t: float,
    degree: int = PHASE0_DEGREE,
    order: int = 0,
) -> np.ndarray:
    """
    Return a row r such that p^(order)(t) = r @ coeffs.

    coeffs has shape (degree + 1,) in ascending monomial order.
    """
    if order < 0:
        raise ValueError("order must be non-negative.")
    row = np.zeros(degree + 1, dtype=float)
    for power in range(order, degree + 1):
        row[power] = _falling_factorial(power, order) * (t ** (power - order))
    return row


def eval_poly(coeffs: np.ndarray, t: float, order: int = 0) -> float:
    """Evaluate a degree-7 polynomial or one of its derivatives."""
    coeffs = np.asarray(coeffs, dtype=float).reshape(-1)
    degree = coeffs.size - 1
    return float(derivative_row(t, degree=degree, order=order) @ coeffs)


def snap_cost_matrix(T: float, degree: int = PHASE0_DEGREE) -> np.ndarray:
    """
    Return Q such that coeffs.T @ Q @ coeffs = integral |p''''(t)|^2 dt.
    """
    if T <= 0.0:
        raise ValueError(f"T must be positive, got {T}.")
    order = PHASE0_S
    Q = np.zeros((degree + 1, degree + 1), dtype=float)
    for i in range(order, degree + 1):
        for j in range(order, degree + 1):
            numerator = _falling_factorial(i, order) * _falling_factorial(j, order)
            exponent = i + j - 2 * order + 1
            Q[i, j] = numerator * (T**exponent) / exponent
    return Q


def hermite_endpoint_matrix(
    T: float,
    degree: int = PHASE0_DEGREE,
    max_order: int = PHASE0_S - 1,
) -> np.ndarray:
    """
    Map coefficients to endpoint derivatives [p(0), ..., p^(r)(0), p(T), ..., p^(r)(T)].
    """
    rows = [derivative_row(0.0, degree=degree, order=order) for order in range(max_order + 1)]
    rows.extend(
        derivative_row(T, degree=degree, order=order) for order in range(max_order + 1)
    )
    return np.vstack(rows)


def coeffs_from_endpoint_derivatives(
    endpoint_values: np.ndarray,
    T: float,
    degree: int = PHASE0_DEGREE,
    max_order: int = PHASE0_S - 1,
) -> np.ndarray:
    """
    Recover degree-7 coefficients from Hermite endpoint data.

    For Phase 0, max_order = 3 gives 8 scalar conditions for 8 coefficients.
    """
    endpoint_values = np.asarray(endpoint_values, dtype=float).reshape(-1)
    expected_size = 2 * (max_order + 1)
    if endpoint_values.size != expected_size:
        raise ValueError(
            f"endpoint_values must have size {expected_size}, got {endpoint_values.size}."
        )
    matrix = hermite_endpoint_matrix(T, degree=degree, max_order=max_order)
    return np.linalg.solve(matrix, endpoint_values)
