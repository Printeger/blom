from __future__ import annotations

import time
from typing import Any

import numpy as np

from phase_1.minco_scalar_baseline import evaluate_trajectory, solve_minco_coefficients
from phase_7.blom_convergence_vs_k import compute_minco_reference


def compute_minco_data(q: np.ndarray, T: np.ndarray) -> dict[str, Any]:
    reference = compute_minco_reference(q, T)
    return {
        "mode": "strict",
        "coeffs": np.asarray(reference["coeffs"], dtype=float),
        "c_vec": np.asarray(reference["c_vec"], dtype=float),
        "J_q": np.asarray(reference["J_q"], dtype=float),
        "cost": float(reference["cost"]),
        "meta": reference,
    }


def sample_minco_trajectory(coeffs: np.ndarray, T: np.ndarray, num_samples: int = 300) -> dict[str, np.ndarray]:
    T = np.asarray(T, dtype=float).reshape(-1)
    total_time = float(np.sum(T))
    t = np.linspace(0.0, total_time, num_samples)
    if t.size:
        t[-1] = np.nextafter(total_time, 0.0)
    y = np.asarray([evaluate_trajectory(coeffs, T, float(tt), order=0) for tt in t], dtype=float)
    return {"t": t, "y": y}


def compute_minco_fd_jacobian_T(q: np.ndarray, T: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    q = np.asarray(q, dtype=float).reshape(-1)
    T = np.asarray(T, dtype=float).reshape(-1)
    zeta_start = np.zeros((4,), dtype=float)
    zeta_end = np.zeros((4,), dtype=float)
    zeta_start[0] = float(q[0])
    zeta_end[0] = float(q[-1])
    base = solve_minco_coefficients(q, T, zeta_start, zeta_end, return_system=False)
    c0 = np.asarray(base["c_vec"], dtype=float)
    J = np.zeros((c0.size, T.size), dtype=float)
    for idx in range(T.size):
        step = np.zeros_like(T)
        step[idx] = eps
        plus = solve_minco_coefficients(q, T + step, zeta_start, zeta_end, return_system=False)
        minus = solve_minco_coefficients(q, T - step, zeta_start, zeta_end, return_system=False)
        J[:, idx] = (np.asarray(plus["c_vec"], dtype=float) - np.asarray(minus["c_vec"], dtype=float)) / (2.0 * eps)
    return J


def runtime_minco_reference(q: np.ndarray, T: np.ndarray) -> tuple[dict[str, Any], float]:
    start = time.perf_counter()
    data = compute_minco_data(q, T)
    return data, time.perf_counter() - start
