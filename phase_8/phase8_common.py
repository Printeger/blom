"""
Shared utilities for Phase 8 interior-first matching validation.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

from phase_0.poly_basis import COEFFS_PER_SEGMENT
from phase_1.minco_scalar_baseline import evaluate_segment
from phase_3.blom_strict_local_kkt import build_local_hessian, solve_kkt
from phase_3.blom_strict_local_qp import DEFAULT_S, build_window
from phase_5.blom_boundary_jump_check import default_boundary_jets
from phase_7.blom_convergence_vs_k import (
    compute_actual_blom_k,
    compute_ideal_truncated_blom_k,
    compute_minco_reference,
    fit_log_error_vs_k,
    make_k_grid,
    representative_case,
)


DEFAULT_RESULTS_DIR = Path("phase_8/results/phase8_validation")


def _validate_time_vector(T: np.ndarray) -> np.ndarray:
    T = np.asarray(T, dtype=float).reshape(-1)
    if T.size < 1:
        raise ValueError("T must describe at least one segment.")
    if not np.all(np.isfinite(T)):
        raise ValueError("T contains non-finite values.")
    if not np.all(T > 0.0):
        raise ValueError("All durations must be strictly positive.")
    return T


def validate_problem_inputs(q: np.ndarray, T: np.ndarray, s: int = DEFAULT_S) -> tuple[np.ndarray, np.ndarray]:
    if s != DEFAULT_S:
        raise ValueError("Phase 8 currently supports only the canonical setting s=4.")
    q = np.asarray(q, dtype=float).reshape(-1)
    T = _validate_time_vector(T)
    if q.shape != (T.size + 1,):
        raise ValueError(f"q must have shape ({T.size + 1},), got {q.shape}.")
    if not np.all(np.isfinite(q)):
        raise ValueError("q contains non-finite values.")
    return q, T


def ensure_results_dirs(base_dir: str | Path = DEFAULT_RESULTS_DIR) -> dict[str, Path]:
    base_dir = Path(base_dir)
    directories = {
        "base": base_dir,
        "interior_matching": base_dir / "interior_matching",
        "boundary_gap": base_dir / "boundary_gap",
        "uniform_vs_nonuniform": base_dir / "uniform_vs_nonuniform",
        "compare": base_dir / "compare",
    }
    for directory in directories.values():
        directory.mkdir(parents=True, exist_ok=True)
    return directories


def save_json(path: str | Path, payload: dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(serialize(payload), indent=2, sort_keys=True), encoding="utf-8")


def save_csv(path: str | Path, header: list[str], rows: list[list[Any]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [",".join(header)]
    for row in rows:
        lines.append(",".join(str(value) for value in row))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def serialize(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, dict):
        return {str(key): serialize(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(item) for item in value]
    return value


def boundary_radius(k: int, mode: str = "default") -> int:
    if k < 1:
        raise ValueError("k must be >= 1.")
    mode_key = str(mode).lower()
    if mode_key in {"minimal", "half", "k_over_2"}:
        return max(int(math.ceil(0.5 * k)), 1)
    if mode_key in {"default", "k"}:
        return int(k)
    if mode_key in {"one_and_half", "three_halves", "1.5k"}:
        return max(int(math.ceil(1.5 * k)), 1)
    raise ValueError(f"Unsupported radius mode {mode!r}.")


def make_interior_sets(M: int, k: int, radius_mode: str = "default") -> dict[str, Any]:
    if M < 1:
        raise ValueError("M must be >= 1.")
    r_k = boundary_radius(k, mode=radius_mode)
    interior_idx = [idx for idx in range(1, M + 1) if idx - r_k >= 1 and idx + r_k <= M]
    boundary_idx = [idx for idx in range(1, M + 1) if idx not in set(interior_idx)]
    return {
        "interior_idx": interior_idx,
        "boundary_idx": boundary_idx,
        "r_k": r_k,
    }


def compute_matching_error_blocks(c_actual: np.ndarray, c_ideal: np.ndarray, s: int = DEFAULT_S) -> np.ndarray:
    c_actual = np.asarray(c_actual, dtype=float).reshape(-1)
    c_ideal = np.asarray(c_ideal, dtype=float).reshape(-1)
    block_size = 2 * s
    if c_actual.shape != c_ideal.shape:
        raise ValueError(f"Coefficient vectors must match, got {c_actual.shape} and {c_ideal.shape}.")
    if c_actual.size % block_size != 0:
        raise ValueError(f"Coefficient vector length must be divisible by {block_size}.")
    diff_blocks = (c_actual - c_ideal).reshape(-1, block_size)
    return np.linalg.norm(diff_blocks, axis=1)


def project_interior_boundary_errors(
    segment_errors: np.ndarray,
    interior_idx: list[int],
    boundary_idx: list[int],
) -> dict[str, float]:
    segment_errors = np.asarray(segment_errors, dtype=float).reshape(-1)
    interior_mask = np.zeros((segment_errors.size,), dtype=bool)
    boundary_mask = np.zeros((segment_errors.size,), dtype=bool)
    for idx in interior_idx:
        if 1 <= idx <= segment_errors.size:
            interior_mask[idx - 1] = True
    for idx in boundary_idx:
        if 1 <= idx <= segment_errors.size:
            boundary_mask[idx - 1] = True

    interior_vals = segment_errors[interior_mask]
    boundary_vals = segment_errors[boundary_mask]
    total_energy = float(np.sum(segment_errors**2))
    interior_energy = float(np.sum(interior_vals**2))
    boundary_energy = float(np.sum(boundary_vals**2))
    return {
        "interior_l2": float(np.linalg.norm(interior_vals)),
        "boundary_l2": float(np.linalg.norm(boundary_vals)),
        "interior_l1": float(np.sum(np.abs(interior_vals))),
        "boundary_l1": float(np.sum(np.abs(boundary_vals))),
        "interior_energy_ratio": interior_energy / max(total_energy, 1e-15),
        "boundary_energy_ratio": boundary_energy / max(total_energy, 1e-15),
    }


def summarize_interior_matching(
    c_actual: np.ndarray,
    c_ideal: np.ndarray,
    M: int,
    k: int,
    s: int = DEFAULT_S,
    radius_mode: str = "default",
) -> dict[str, Any]:
    segment_errors = compute_matching_error_blocks(c_actual, c_ideal, s=s)
    if segment_errors.shape != (M,):
        raise ValueError(f"Expected {M} segment errors, got {segment_errors.shape}.")
    sets = make_interior_sets(M, k, radius_mode=radius_mode)
    decomposition = project_interior_boundary_errors(segment_errors, sets["interior_idx"], sets["boundary_idx"])

    def _linf(indices: list[int]) -> float:
        if not indices:
            return 0.0
        return float(np.max(segment_errors[np.asarray(indices, dtype=int) - 1]))

    return {
        "k": int(k),
        "full_l2": float(np.linalg.norm(segment_errors)),
        "interior_l2": decomposition["interior_l2"],
        "boundary_l2": decomposition["boundary_l2"],
        "full_linf": float(np.max(segment_errors)) if segment_errors.size else 0.0,
        "interior_linf": _linf(sets["interior_idx"]),
        "boundary_linf": _linf(sets["boundary_idx"]),
        "segment_errors": segment_errors,
        "interior_idx": sets["interior_idx"],
        "boundary_idx": sets["boundary_idx"],
        "r_k": sets["r_k"],
        "interior_energy_ratio": decomposition["interior_energy_ratio"],
        "boundary_energy_ratio": decomposition["boundary_energy_ratio"],
    }


def sample_uniform_time(M: int, h: float) -> np.ndarray:
    if M < 1:
        raise ValueError("M must be >= 1.")
    if h <= 0.0:
        raise ValueError("h must be strictly positive.")
    return np.full((M,), float(h), dtype=float)


def sample_bounded_nonuniform_time(M: int, T_min: float, T_max: float, rng: np.random.Generator) -> np.ndarray:
    if M < 1:
        raise ValueError("M must be >= 1.")
    if T_min <= 0.0 or T_max <= 0.0:
        raise ValueError("T_min and T_max must be strictly positive.")
    if T_min > T_max:
        raise ValueError("T_min must be <= T_max.")
    return rng.uniform(T_min, T_max, size=M).astype(float)


def default_q_sampler(rng: np.random.Generator, M: int) -> np.ndarray:
    increments = rng.normal(loc=0.0, scale=0.9, size=M)
    q = np.zeros((M + 1,), dtype=float)
    q[1:] = np.cumsum(increments)
    return q


def representative_phase8_case() -> tuple[np.ndarray, np.ndarray]:
    return representative_case()


def _gamma_split(gamma: np.ndarray, s: int = DEFAULT_S) -> tuple[np.ndarray, np.ndarray]:
    gamma = np.asarray(gamma, dtype=float).reshape(2 * (s - 1))
    return gamma[: s - 1], gamma[s - 1 :]


def extract_reference_gamma(reference: dict[str, Any], T: np.ndarray, i: int, k: int, s: int = DEFAULT_S) -> dict[str, Any]:
    coeffs = np.asarray(reference["coeffs"], dtype=float)
    T = _validate_time_vector(T)
    window = build_window(i, k, T.size)
    gamma_left = np.zeros((s - 1,), dtype=float)
    gamma_right = np.zeros((s - 1,), dtype=float)
    if not window["touches_left_boundary"]:
        left_seg = window["L"] - 1
        duration = float(T[left_seg - 1])
        gamma_left = np.asarray(
            [evaluate_segment(coeffs[left_seg - 1], duration, order=order) for order in range(s, 2 * s - 1)],
            dtype=float,
        )
    if not window["touches_right_boundary"]:
        right_seg = window["R"]
        duration = float(T[right_seg - 1])
        gamma_right = np.asarray(
            [evaluate_segment(coeffs[right_seg - 1], duration, order=order) for order in range(s, 2 * s - 1)],
            dtype=float,
        )
    gamma = np.concatenate((gamma_left, gamma_right))
    return {
        "gamma": gamma,
        "gamma_left": gamma_left,
        "gamma_right": gamma_right,
        "window": window,
    }


def build_augmented_local_constraints(
    q: np.ndarray,
    T: np.ndarray,
    i: int,
    k: int,
    *,
    gamma: np.ndarray | None = None,
    zeta_start: np.ndarray | None = None,
    zeta_end: np.ndarray | None = None,
    s: int = DEFAULT_S,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    from phase_1.minco_scalar_baseline import beta, beta_d

    q, T = validate_problem_inputs(q, T, s=s)
    window = build_window(i, k, T.size)
    gamma_vec = np.zeros((2 * (s - 1),), dtype=float) if gamma is None else np.asarray(gamma, dtype=float).reshape(2 * (s - 1))
    gamma_left, gamma_right = _gamma_split(gamma_vec, s=s)
    zeta_start, zeta_end = default_boundary_jets(q, s=s) if zeta_start is None or zeta_end is None else (
        np.asarray(zeta_start, dtype=float).reshape(s),
        np.asarray(zeta_end, dtype=float).reshape(s),
    )

    segments = window["segments"]
    m_seg = len(segments)
    block_size = 2 * s
    num_unknowns = block_size * m_seg

    num_rows = 2 * m_seg + (s - 1) * max(m_seg - 1, 0)
    num_rows += s - 1
    num_rows += s - 1
    G = np.zeros((num_rows, num_unknowns), dtype=float)
    d_vec = np.zeros((num_rows,), dtype=float)
    row = 0

    for local_index, segment_index in enumerate(segments):
        duration = float(T[segment_index - 1])
        coeff_slice = slice(local_index * block_size, (local_index + 1) * block_size)
        G[row, coeff_slice] = beta(0.0)
        d_vec[row] = float(q[segment_index - 1])
        row += 1
        G[row, coeff_slice] = beta(duration)
        d_vec[row] = float(q[segment_index])
        row += 1

    for local_index, segment_index in enumerate(range(window["L"], window["R"])):
        duration = float(T[segment_index - 1])
        left_slice = slice(local_index * block_size, (local_index + 1) * block_size)
        right_slice = slice((local_index + 1) * block_size, (local_index + 2) * block_size)
        for order in range(1, s):
            G[row, left_slice] = beta_d(duration, order)
            G[row, right_slice] = -beta_d(0.0, order)
            row += 1

    first_slice = slice(0, block_size)
    if window["touches_left_boundary"]:
        for order in range(1, s):
            G[row, first_slice] = beta_d(0.0, order)
            d_vec[row] = float(zeta_start[order])
            row += 1
    else:
        for order in range(s, 2 * s - 1):
            G[row, first_slice] = beta_d(0.0, order)
            d_vec[row] = float(gamma_left[order - s])
            row += 1

    last_slice = slice((m_seg - 1) * block_size, m_seg * block_size)
    duration = float(T[window["R"] - 1])
    if window["touches_right_boundary"]:
        for order in range(1, s):
            G[row, last_slice] = beta_d(duration, order)
            d_vec[row] = float(zeta_end[order])
            row += 1
    else:
        for order in range(s, 2 * s - 1):
            G[row, last_slice] = beta_d(duration, order)
            d_vec[row] = float(gamma_right[order - s])
            row += 1

    if row != num_rows:
        raise RuntimeError(f"Augmented local constraint assembly bug: filled {row}, expected {num_rows}.")

    meta = {
        "window": window,
        "gamma": gamma_vec,
        "gamma_left": gamma_left,
        "gamma_right": gamma_right,
        "touches_left_boundary": window["touches_left_boundary"],
        "touches_right_boundary": window["touches_right_boundary"],
    }
    return G, d_vec, meta


def solve_augmented_local_window(
    q: np.ndarray,
    T: np.ndarray,
    i: int,
    k: int,
    *,
    gamma: np.ndarray | None = None,
    zeta_start: np.ndarray | None = None,
    zeta_end: np.ndarray | None = None,
    s: int = DEFAULT_S,
) -> dict[str, Any]:
    q, T = validate_problem_inputs(q, T, s=s)
    window = build_window(i, k, T.size)
    T_window = np.asarray([T[segment - 1] for segment in window["segments"]], dtype=float)
    H = build_local_hessian(T_window, s=s)
    G, d_vec, meta = build_augmented_local_constraints(
        q,
        T,
        i,
        k,
        gamma=gamma,
        zeta_start=zeta_start,
        zeta_end=zeta_end,
        s=s,
    )
    solved = solve_kkt(H, G, d_vec)
    coeffs_local = solved["c_vec"].reshape(window["m"], 2 * s)
    center_local_index = window["center_local_index"]
    center_coeff = coeffs_local[center_local_index].copy()

    left_trace_error = np.zeros((0,), dtype=float)
    right_trace_error = np.zeros((0,), dtype=float)
    if not window["touches_left_boundary"]:
        left_trace_error = np.asarray(
            [
                abs(evaluate_segment(coeffs_local[0], 0.0, order=order) - meta["gamma_left"][order - s])
                for order in range(s, 2 * s - 1)
            ],
            dtype=float,
        )
    if not window["touches_right_boundary"]:
        duration = float(T_window[-1])
        right_trace_error = np.asarray(
            [
                abs(evaluate_segment(coeffs_local[-1], duration, order=order) - meta["gamma_right"][order - s])
                for order in range(s, 2 * s - 1)
            ],
            dtype=float,
        )

    return {
        "i": int(i),
        "k": int(k),
        "s": s,
        "window": window,
        "T_window": T_window,
        "H": H,
        "G": G,
        "d": d_vec,
        "gamma": meta["gamma"],
        "coeffs_local": coeffs_local,
        "center_coeff": center_coeff,
        "center_c_vec": center_coeff.reshape(-1),
        "objective": float(solved["objective"]),
        "residual_norm": float(solved["kkt_residual"]),
        "trace_error_left": left_trace_error,
        "trace_error_right": right_trace_error,
        "raw": solved,
    }


def compute_boundary_response_matrix(
    q: np.ndarray,
    T: np.ndarray,
    i: int,
    k: int,
    s: int = DEFAULT_S,
) -> dict[str, Any]:
    base = solve_augmented_local_window(q, T, i, k, gamma=np.zeros((2 * (s - 1),), dtype=float), s=s)
    matrix = np.zeros((2 * s, 2 * (s - 1)), dtype=float)
    for col in range(2 * (s - 1)):
        gamma = np.zeros((2 * (s - 1),), dtype=float)
        gamma[col] = 1.0
        perturbed = solve_augmented_local_window(q, T, i, k, gamma=gamma, s=s)
        matrix[:, col] = perturbed["center_c_vec"] - base["center_c_vec"]
    return {
        "base": base,
        "M_matrix": matrix,
        "operator_norm": float(np.linalg.norm(matrix, ord=2)),
    }


def compute_reference_window_coefficients(
    q: np.ndarray,
    T: np.ndarray,
    k: int,
    *,
    reference: dict[str, Any] | None = None,
    s: int = DEFAULT_S,
) -> dict[str, Any]:
    q, T = validate_problem_inputs(q, T, s=s)
    reference = compute_minco_reference(q, T, s=s) if reference is None else reference
    coeffs = np.zeros((T.size, 2 * s), dtype=float)
    gamma_norms = np.zeros((T.size,), dtype=float)
    operator_norms = np.zeros((T.size,), dtype=float)
    raw_gap = np.zeros((T.size,), dtype=float)
    for seg_idx in range(1, T.size + 1):
        gamma_data = extract_reference_gamma(reference, T, seg_idx, k, s=s)
        solved = solve_augmented_local_window(q, T, seg_idx, k, gamma=gamma_data["gamma"], s=s)
        coeffs[seg_idx - 1] = solved["center_coeff"]
        gamma_norms[seg_idx - 1] = float(np.linalg.norm(gamma_data["gamma"]))
        response = compute_boundary_response_matrix(q, T, seg_idx, k, s=s)
        operator_norms[seg_idx - 1] = response["operator_norm"]
        raw_gap[seg_idx - 1] = float(np.linalg.norm(coeffs[seg_idx - 1] - response["base"]["center_coeff"]))
    return {
        "kind": "reference_window",
        "k": int(k),
        "coeffs": coeffs,
        "c_vec": coeffs.reshape(-1),
        "gamma_norms": gamma_norms,
        "operator_norms": operator_norms,
        "raw_vs_reference_segment_l2": raw_gap,
    }


def compute_phase8_matching_triplet(
    q: np.ndarray,
    T: np.ndarray,
    k: int,
    *,
    s: int = DEFAULT_S,
    reference: dict[str, Any] | None = None,
) -> dict[str, Any]:
    q, T = validate_problem_inputs(q, T, s=s)
    reference = compute_minco_reference(q, T, s=s) if reference is None else reference
    actual = compute_actual_blom_k(q, T, k, s=s, scheme="C")
    ideal = compute_ideal_truncated_blom_k(q, T, k, s=s, reference=reference)
    reference_window = compute_reference_window_coefficients(q, T, k, reference=reference, s=s)
    return {
        "reference": reference,
        "actual": actual,
        "ideal": ideal,
        "reference_window": reference_window,
    }


def prepare_case(
    q: np.ndarray | None = None,
    T: np.ndarray | None = None,
    *,
    M: int = 8,
    seed: int = 42,
    q_sampler: Callable[[np.random.Generator, int], np.ndarray] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    if q is not None or T is not None:
        if q is None or T is None:
            raise ValueError("q and T must either both be provided or both be omitted.")
        return validate_problem_inputs(q, T)
    if M == 8:
        return representative_phase8_case()
    rng = np.random.default_rng(seed)
    q_sampler = default_q_sampler if q_sampler is None else q_sampler
    q_val = q_sampler(rng, M)
    T_val = rng.uniform(0.6, 1.6, size=M).astype(float)
    return validate_problem_inputs(q_val, T_val)


def interior_fit_summary(records: list[dict[str, Any]]) -> dict[str, Any]:
    k_values = [int(record["k"]) for record in records]
    full_fit = fit_log_error_vs_k(k_values, [float(record["full_l2"]) for record in records])
    interior_fit = fit_log_error_vs_k(k_values, [float(record["interior_l2"]) for record in records])
    boundary_fit = fit_log_error_vs_k(k_values, [float(record["boundary_l2"]) for record in records])
    return {
        "full": full_fit,
        "interior": interior_fit,
        "boundary": boundary_fit,
    }


__all__ = [
    "DEFAULT_RESULTS_DIR",
    "DEFAULT_S",
    "COEFFS_PER_SEGMENT",
    "boundary_radius",
    "compute_boundary_response_matrix",
    "compute_matching_error_blocks",
    "compute_phase8_matching_triplet",
    "compute_reference_window_coefficients",
    "default_q_sampler",
    "ensure_results_dirs",
    "extract_reference_gamma",
    "fit_log_error_vs_k",
    "interior_fit_summary",
    "make_interior_sets",
    "make_k_grid",
    "prepare_case",
    "project_interior_boundary_errors",
    "representative_phase8_case",
    "sample_bounded_nonuniform_time",
    "sample_uniform_time",
    "save_csv",
    "save_json",
    "serialize",
    "solve_augmented_local_window",
    "summarize_interior_matching",
    "validate_problem_inputs",
]
