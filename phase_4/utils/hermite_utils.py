"""
Shared Hermite and quadratic-form utilities for Phase 4.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np


DEFAULT_S2_DEGREE = 3
DEFAULT_S4_DEGREE = 7


def normalized_monomial_row(tau: float, degree: int) -> np.ndarray:
    tau = float(tau)
    return np.asarray([tau**power for power in range(degree + 1)], dtype=float)


def normalized_derivative_row(tau: float, order: int, degree: int) -> np.ndarray:
    tau = float(tau)
    if order < 0:
        raise ValueError("order must be >= 0.")
    row = np.zeros((degree + 1,), dtype=float)
    for power in range(order, degree + 1):
        coeff = math.factorial(power) / math.factorial(power - order)
        row[power] = coeff * (tau ** (power - order))
    return row


def build_C4() -> np.ndarray:
    rows = [normalized_derivative_row(0.0, order, DEFAULT_S4_DEGREE) for order in range(4)]
    rows.extend(normalized_derivative_row(1.0, order, DEFAULT_S4_DEGREE) for order in range(4))
    return np.vstack(rows)


def build_G4() -> np.ndarray:
    degree = DEFAULT_S4_DEGREE
    G4 = np.zeros((degree + 1, degree + 1), dtype=float)
    for alpha in range(4, degree + 1):
        coeff_alpha = math.factorial(alpha) / math.factorial(alpha - 4)
        for beta in range(4, degree + 1):
            coeff_beta = math.factorial(beta) / math.factorial(beta - 4)
            power = alpha + beta - 8
            G4[alpha, beta] = coeff_alpha * coeff_beta / (power + 1)
    return G4


def build_R4() -> np.ndarray:
    C4 = build_C4()
    H4 = np.linalg.inv(C4)
    G4 = build_G4()
    return H4.T @ G4 @ H4


def D_matrix(h: float) -> np.ndarray:
    h = float(h)
    return np.diag([h, h**2, h**3]).astype(float)


def Lambda8(h: float) -> np.ndarray:
    h = float(h)
    return np.diag([h ** (-power) for power in range(8)]).astype(float)


def build_m_plus(h: float) -> np.ndarray:
    h = float(h)
    return np.asarray([h, 0.5 * h**2, (h**3) / 6.0], dtype=float)


def build_m_minus(h: float) -> np.ndarray:
    h = float(h)
    return np.asarray([h, -0.5 * h**2, (h**3) / 6.0], dtype=float)


def build_Pq() -> np.ndarray:
    Pq = np.zeros((8, 2), dtype=float)
    Pq[0, 0] = 1.0
    Pq[4, 1] = 1.0
    return Pq


def build_Px(h: float) -> np.ndarray:
    D_h = D_matrix(h)
    Px = np.zeros((8, 6), dtype=float)
    Px[1:4, 0:3] = D_h
    Px[5:8, 3:6] = D_h
    return Px


def build_H_mid(h: float) -> np.ndarray:
    h = float(h)
    R4 = build_R4()
    Px = build_Px(h)
    return 2.0 * (h ** -7) * Px.T @ R4 @ Px


def build_g_mid(h: float) -> np.ndarray:
    h = float(h)
    R4 = build_R4()
    Px = build_Px(h)
    Pq = build_Pq()
    return -2.0 * (h ** -7) * Px.T @ R4 @ Pq


def build_gamma_mid(q_left: float, q_right: float, h: float) -> float:
    h = float(h)
    q_vec = np.asarray([q_left, q_right], dtype=float)
    R4 = build_R4()
    Pq = build_Pq()
    return float((h ** -7) * q_vec @ (Pq.T @ R4 @ Pq) @ q_vec)


def outer_rank_one_hessian(side: str, h: float) -> np.ndarray:
    if side not in {"left", "right"}:
        raise ValueError(f"Unsupported side {side!r}.")
    vec = build_m_minus(h) if side == "left" else build_m_plus(h)
    return 504.0 * (float(h) ** -7) * np.outer(vec, vec)


def outer_linear_term(side: str, q_left: float, q_right: float, h: float) -> np.ndarray:
    if side not in {"left", "right"}:
        raise ValueError(f"Unsupported side {side!r}.")
    delta = float(q_right - q_left)
    vec = build_m_minus(h) if side == "left" else build_m_plus(h)
    return 504.0 * (float(h) ** -7) * delta * vec


def outer_constant_term(q_left: float, q_right: float, h: float) -> float:
    delta = float(q_right - q_left)
    return 252.0 * (float(h) ** -7) * delta**2


def y_vector(q_left: float, q_right: float, x_left: np.ndarray, x_right: np.ndarray, h: float) -> np.ndarray:
    x_left = np.asarray(x_left, dtype=float).reshape(3)
    x_right = np.asarray(x_right, dtype=float).reshape(3)
    return np.concatenate(([float(q_left)], D_matrix(h) @ x_left, [float(q_right)], D_matrix(h) @ x_right))


def hermite_reconstruct_center_segment(
    q_left: float,
    q_right: float,
    x_left: np.ndarray,
    x_right: np.ndarray,
    h: float,
) -> dict[str, np.ndarray]:
    C4 = build_C4()
    H4 = np.linalg.inv(C4)
    y = y_vector(q_left, q_right, x_left, x_right, h)
    alpha = H4 @ y
    coeff = Lambda8(h) @ alpha
    return {
        "y": y,
        "alpha": alpha,
        "coeff": coeff,
        "C4": C4,
        "H4": H4,
        "Lambda8": Lambda8(h),
    }


def evaluate_physical_polynomial(coeff: np.ndarray, t: float, order: int = 0) -> float:
    coeff = np.asarray(coeff, dtype=float).reshape(-1)
    t = float(t)
    value = 0.0
    for power, coefficient in enumerate(coeff):
        if power < order:
            continue
        scale = math.factorial(power) / math.factorial(power - order)
        value += coefficient * scale * (t ** (power - order))
    return float(value)


def quadratic_objective_value(A_mat: np.ndarray, B_vec: np.ndarray, C_scalar: float, x_vec: np.ndarray) -> float:
    A_mat = np.asarray(A_mat, dtype=float)
    B_vec = np.asarray(B_vec, dtype=float).reshape(-1)
    x_vec = np.asarray(x_vec, dtype=float).reshape(-1)
    return float(0.5 * x_vec @ A_mat @ x_vec - x_vec @ B_vec + C_scalar)


def serialize_array(arr: np.ndarray) -> list[Any]:
    return np.asarray(arr, dtype=float).tolist()

