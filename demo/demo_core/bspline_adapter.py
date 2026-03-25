from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


def open_uniform_knots(n_ctrl: int, degree: int) -> np.ndarray:
    n_knots = n_ctrl + degree + 1
    knots = np.zeros((n_knots,), dtype=float)
    knots[degree : n_ctrl + 1] = np.linspace(0.0, 1.0, n_ctrl - degree + 1)
    knots[n_ctrl + 1 :] = 1.0
    return knots


def basis_value(i: int, degree: int, u: float, knots: np.ndarray) -> float:
    if degree == 0:
        if knots[i] <= u < knots[i + 1] or (u == 1.0 and knots[i + 1] == 1.0):
            return 1.0
        return 0.0
    left_denom = knots[i + degree] - knots[i]
    right_denom = knots[i + degree + 1] - knots[i + 1]
    left = 0.0 if left_denom == 0.0 else (u - knots[i]) / left_denom * basis_value(i, degree - 1, u, knots)
    right = 0.0 if right_denom == 0.0 else (knots[i + degree + 1] - u) / right_denom * basis_value(i + 1, degree - 1, u, knots)
    return left + right


def basis_matrix(control_points: np.ndarray, degree: int = 3, num_samples: int = 300) -> tuple[np.ndarray, np.ndarray]:
    ctrl = np.asarray(control_points, dtype=float).reshape(-1)
    knots = open_uniform_knots(ctrl.size, degree)
    u = np.linspace(0.0, 1.0, num_samples)
    B = np.zeros((num_samples, ctrl.size), dtype=float)
    for row, uu in enumerate(u):
        for col in range(ctrl.size):
            B[row, col] = basis_value(col, degree, float(uu), knots)
    return u, B


def evaluate_bspline(control_points: np.ndarray, degree: int = 3, num_samples: int = 300) -> dict[str, np.ndarray]:
    ctrl = np.asarray(control_points, dtype=float).reshape(-1)
    u, B = basis_matrix(ctrl, degree=degree, num_samples=num_samples)
    y = B @ ctrl
    return {"u": u, "y": y, "basis": B, "control_points": ctrl}


def perturb_control_point(control_points: np.ndarray, index: int, delta: float) -> np.ndarray:
    ctrl = np.asarray(control_points, dtype=float).reshape(-1).copy()
    idx = int(np.clip(index, 0, ctrl.size - 1))
    ctrl[idx] += float(delta)
    return ctrl


def default_control_points_from_waypoints(q: np.ndarray) -> np.ndarray:
    q = np.asarray(q, dtype=float).reshape(-1)
    return q.copy()


def sampled_response(control_points: np.ndarray, degree: int = 3, num_samples: int = 300) -> np.ndarray:
    _, B = basis_matrix(control_points, degree=degree, num_samples=num_samples)
    return B
