"""
Main Phase 3 BLOM-Strict local QP construction and solution API.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from phase_1.minco_scalar_baseline import evaluate_segment
from phase_3.blom_strict_feasible_init import build_feasible_local_spline
from phase_3.blom_strict_local_kkt import (
    DEFAULT_S,
    build_local_constraints,
    build_local_hessian,
    solve_kkt,
    solve_reduced_qp,
)


DEFAULT_K = 2
DEFAULT_CONTINUITY_ORDER = 2 * DEFAULT_S - 2


def build_window(i: int, k: int, M: int) -> dict[str, Any]:
    """Build the Phase 3 local window metadata."""
    if M < 1:
        raise ValueError("M must be >= 1.")
    if not 1 <= i <= M:
        raise ValueError(f"i must be in [1, {M}], got {i}.")
    if k < 0:
        raise ValueError("k must be >= 0.")

    half = k // 2
    left = max(1, i - half)
    right = min(M, i + half)
    segments = list(range(left, right + 1))
    if left == 1 and right == M:
        window_type = "full"
    elif left == 1:
        window_type = "left"
    elif right == M:
        window_type = "right"
    else:
        window_type = "interior"
    return {
        "center": int(i),
        "k": int(k),
        "L": int(left),
        "R": int(right),
        "m": int(right - left + 1),
        "segments": segments,
        "knots": list(range(left - 1, right + 1)),
        "center_local_index": int(i - left),
        "touches_left_boundary": left == 1,
        "touches_right_boundary": right == M,
        "window_type": window_type,
    }


def _normalize_problem_inputs(
    q: np.ndarray,
    T: np.ndarray,
    zeta_start: np.ndarray | None,
    zeta_end: np.ndarray | None,
    s: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    q = np.asarray(q, dtype=float).reshape(-1)
    T = np.asarray(T, dtype=float).reshape(-1)
    if q.shape != (T.size + 1,):
        raise ValueError(f"q must have shape ({T.size + 1},), got {q.shape}.")
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
    return q, T, zeta_start, zeta_end


def build_local_problem(
    q: np.ndarray,
    T: np.ndarray,
    i: int,
    k: int,
    zeta_start: np.ndarray | None = None,
    zeta_end: np.ndarray | None = None,
    s: int = DEFAULT_S,
) -> dict[str, Any]:
    """Construct the full local BLOM-Strict equality-constrained QP data."""
    q, T, zeta_start, zeta_end = _normalize_problem_inputs(q, T, zeta_start, zeta_end, s=s)
    window = build_window(i, k, T.size)
    T_window = np.asarray([T[segment - 1] for segment in window["segments"]], dtype=float)
    sigma = np.concatenate(([0.0], np.cumsum(T_window)))
    feasible = build_feasible_local_spline(q, T, i, k, zeta_start=zeta_start, zeta_end=zeta_end, s=s)
    H = build_local_hessian(T_window, s=s)
    G, d_vec = build_local_constraints(q, T, i, k, zeta_start=zeta_start, zeta_end=zeta_end, s=s)

    window["sigma"] = sigma
    window["T_window"] = T_window
    window["T_local"] = float(np.sum(T_window))
    return {
        "q": q,
        "T": T,
        "zeta_start": zeta_start,
        "zeta_end": zeta_end,
        "s": s,
        "i": int(i),
        "k": int(k),
        "window": window,
        "H": H,
        "G": G,
        "d": d_vec,
        "feasible_init": feasible,
    }


def compute_interpolation_errors(
    coeffs: np.ndarray,
    q: np.ndarray,
    window: dict[str, Any],
) -> np.ndarray:
    """Return start/end interpolation errors for every segment in the window."""
    q = np.asarray(q, dtype=float).reshape(-1)
    coeffs = np.asarray(coeffs, dtype=float)
    errors = []
    for local_index, segment_index in enumerate(window["segments"]):
        duration = float(window["T_window"][local_index])
        errors.append(abs(evaluate_segment(coeffs[local_index], 0.0, order=0) - q[segment_index - 1]))
        errors.append(abs(evaluate_segment(coeffs[local_index], duration, order=0) - q[segment_index]))
    return np.asarray(errors, dtype=float)


def compute_continuity_jumps(
    coeffs: np.ndarray,
    window: dict[str, Any],
    max_order: int = DEFAULT_CONTINUITY_ORDER,
) -> np.ndarray:
    """Return derivative jumps across interior window knots for orders 0..max_order."""
    coeffs = np.asarray(coeffs, dtype=float)
    num_junctions = max(window["m"] - 1, 0)
    jumps = np.zeros((num_junctions, max_order + 1), dtype=float)
    for local_index in range(num_junctions):
        duration = float(window["T_window"][local_index])
        left_coeff = coeffs[local_index]
        right_coeff = coeffs[local_index + 1]
        for order in range(max_order + 1):
            left_val = evaluate_segment(left_coeff, duration, order=order)
            right_val = evaluate_segment(right_coeff, 0.0, order=order)
            jumps[local_index, order] = left_val - right_val
    return jumps


def compute_boundary_jet_errors(
    coeffs: np.ndarray,
    window: dict[str, Any],
    zeta_start: np.ndarray,
    zeta_end: np.ndarray,
    s: int = DEFAULT_S,
) -> dict[str, np.ndarray]:
    """Return physical boundary jet residuals on the sides where they apply."""
    coeffs = np.asarray(coeffs, dtype=float)
    result: dict[str, np.ndarray] = {
        "left": np.zeros((0,), dtype=float),
        "right": np.zeros((0,), dtype=float),
    }
    if window["touches_left_boundary"]:
        result["left"] = np.asarray(
            [abs(evaluate_segment(coeffs[0], 0.0, order=order) - zeta_start[order]) for order in range(s)],
            dtype=float,
        )
    if window["touches_right_boundary"]:
        last_coeff = coeffs[-1]
        duration = float(window["T_window"][-1])
        result["right"] = np.asarray(
            [abs(evaluate_segment(last_coeff, duration, order=order) - zeta_end[order]) for order in range(s)],
            dtype=float,
        )
    return result


def compute_natural_boundary_residuals(
    coeffs: np.ndarray,
    window: dict[str, Any],
    s: int = DEFAULT_S,
) -> dict[str, np.ndarray]:
    """Return natural-boundary residuals p^(s)..p^(2s-2) at artificial boundaries."""
    coeffs = np.asarray(coeffs, dtype=float)
    result: dict[str, np.ndarray] = {
        "left": np.zeros((0,), dtype=float),
        "right": np.zeros((0,), dtype=float),
    }
    if not window["touches_left_boundary"]:
        result["left"] = np.asarray(
            [abs(evaluate_segment(coeffs[0], 0.0, order=order)) for order in range(s, 2 * s - 1)],
            dtype=float,
        )
    if not window["touches_right_boundary"]:
        last_coeff = coeffs[-1]
        duration = float(window["T_window"][-1])
        result["right"] = np.asarray(
            [abs(evaluate_segment(last_coeff, duration, order=order)) for order in range(s, 2 * s - 1)],
            dtype=float,
        )
    return result


def sample_local_trajectory(
    coeffs: np.ndarray,
    window: dict[str, Any],
    *,
    num_per_segment: int = 80,
    orders: tuple[int, ...] = (0, 1, 2, 3, 4),
) -> dict[str, np.ndarray]:
    """Sample the local trajectory on a piecewise-uniform grid."""
    coeffs = np.asarray(coeffs, dtype=float)
    times: list[float] = []
    values: dict[str, list[float]] = {f"order_{order}": [] for order in orders}
    cursor = 0.0
    for local_index, duration in enumerate(window["T_window"]):
        local_times = np.linspace(
            0.0,
            float(duration),
            num=num_per_segment,
            endpoint=(local_index == window["m"] - 1),
        )
        coeff = coeffs[local_index]
        for local_time in local_times:
            times.append(cursor + float(local_time))
            for order in orders:
                values[f"order_{order}"].append(evaluate_segment(coeff, float(local_time), order=order))
        cursor += float(duration)
    sampled = {"sigma": np.asarray(times, dtype=float)}
    for key, value in values.items():
        sampled[key] = np.asarray(value, dtype=float)
    return sampled


def summarize_solution(solution: dict[str, Any]) -> dict[str, Any]:
    """Collect the main verification metrics for one local solution."""
    problem = solution["problem"]
    coeffs = solution["coeffs"]
    continuity = compute_continuity_jumps(coeffs, problem["window"])
    interpolation = compute_interpolation_errors(coeffs, problem["q"], problem["window"])
    boundary = compute_boundary_jet_errors(
        coeffs,
        problem["window"],
        problem["zeta_start"],
        problem["zeta_end"],
        s=problem["s"],
    )
    natural = compute_natural_boundary_residuals(coeffs, problem["window"], s=problem["s"])

    max_cs_minus_1 = 0.0
    if continuity.size:
        max_cs_minus_1 = float(np.max(np.abs(continuity[:, : problem["s"]])))
    max_c2s_minus_2 = float(np.max(np.abs(continuity))) if continuity.size else 0.0
    max_boundary = max(
        float(np.max(boundary["left"])) if boundary["left"].size else 0.0,
        float(np.max(boundary["right"])) if boundary["right"].size else 0.0,
    )
    max_natural = max(
        float(np.max(natural["left"])) if natural["left"].size else 0.0,
        float(np.max(natural["right"])) if natural["right"].size else 0.0,
    )
    return {
        "max_interp_error": float(np.max(interpolation)) if interpolation.size else 0.0,
        "max_Cs_minus_1_jump": max_cs_minus_1,
        "max_C2s_minus_2_jump": max_c2s_minus_2,
        "max_boundary_jet_error": max_boundary,
        "max_natural_bc_residual": max_natural,
        "continuity_jumps": continuity,
        "interpolation_errors": interpolation,
        "boundary_jet_errors": boundary,
        "natural_bc_residuals": natural,
    }


def solve_blom_strict_local_qp(problem: dict[str, Any], method: str = "kkt", **kwargs: Any) -> dict[str, Any]:
    """Solve one BLOM-Strict local problem using KKT or reduced-space linear algebra."""
    H = problem["H"]
    G = problem["G"]
    d_vec = problem["d"]

    if method == "kkt":
        raw = solve_kkt(H, G, d_vec)
    elif method == "reduced":
        raw = solve_reduced_qp(H, G, d_vec, iterative=False, **kwargs)
    elif method == "reduced_iterative":
        raw = solve_reduced_qp(H, G, d_vec, iterative=True, **kwargs)
    else:
        raise ValueError(f"Unsupported method {method!r}.")

    coeffs = raw["c_vec"].reshape(problem["window"]["m"], 2 * problem["s"])
    solution = {
        **raw,
        "problem": problem,
        "coeffs": coeffs,
        "method": method,
    }
    solution["summary"] = summarize_solution(solution)
    return solution


def extract_segment_coeff(solution: dict[str, Any], seg_idx: int) -> np.ndarray:
    """Extract one segment's coefficient block by global segment index."""
    window = solution["problem"]["window"]
    if seg_idx not in window["segments"]:
        raise KeyError(f"seg_idx {seg_idx} is not in window segments {window['segments']}.")
    local_index = window["segments"].index(seg_idx)
    return solution["coeffs"][local_index].copy()

