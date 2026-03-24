"""
Phase 2 validation utilities for BLOM necessity theory.

This module treats the Phase 1 MINCO solver as a fixed backend and builds
Jacobian extraction, finite-difference checks, influence analysis, bandwidth
statistics, and artifact export on top of it.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import json
from typing import Any

import numpy as np

from phase_1.minco_scalar_baseline import (
    COEFFS_PER_SEGMENT,
    build_system_matrix,
    build_rhs,
    solve_minco_coefficients,
)


DEFAULT_SEED = 42
DEFAULT_TOL = 1e-10


@dataclass
class Phase2Case:
    case_name: str
    q: np.ndarray
    T: np.ndarray
    zeta_start: np.ndarray
    zeta_end: np.ndarray


def ensure_results_dirs(base_dir: str | Path = Path("phase_2/results")) -> dict[str, Path]:
    """Create the standard Phase 2 results directories and return them."""
    base_dir = Path(base_dir)
    figures = base_dir / "figures"
    tables = base_dir / "tables"
    logs = base_dir / "logs"
    for directory in (base_dir, figures, tables, logs):
        directory.mkdir(parents=True, exist_ok=True)
    return {"base": base_dir, "figures": figures, "tables": tables, "logs": logs}


def default_boundary_jets(q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Use zero velocity/acceleration/jerk boundary jets around fixed endpoint positions."""
    q = np.asarray(q, dtype=float).reshape(-1)
    zeta_start = np.array([q[0], 0.0, 0.0, 0.0], dtype=float)
    zeta_end = np.array([q[-1], 0.0, 0.0, 0.0], dtype=float)
    return zeta_start, zeta_end


def make_uniform_time_case(M_seg: int, case_name: str = "uniform_time") -> Phase2Case:
    """Create a smooth deterministic case with uniform segment durations."""
    waypoints = np.linspace(0.0, 2.0 * np.pi, M_seg + 1)
    q = 0.8 * np.sin(waypoints) + 0.15 * np.cos(2.0 * waypoints)
    T = np.ones(M_seg, dtype=float)
    zeta_start, zeta_end = default_boundary_jets(q)
    return Phase2Case(case_name=case_name, q=q, T=T, zeta_start=zeta_start, zeta_end=zeta_end)


def make_nonuniform_time_case(
    M_seg: int,
    *,
    seed: int = DEFAULT_SEED,
    case_name: str = "nonuniform_time",
) -> Phase2Case:
    """Create a deterministic random-looking case with nonuniform durations."""
    rng = np.random.default_rng(seed)
    q = np.concatenate(([0.0], rng.normal(scale=0.8, size=M_seg)))
    T = rng.uniform(0.5, 1.5, size=M_seg)
    zeta_start, zeta_end = default_boundary_jets(q)
    return Phase2Case(case_name=case_name, q=q, T=T, zeta_start=zeta_start, zeta_end=zeta_end)


def make_symmetric_case(M_seg: int, case_name: str = "symmetric") -> Phase2Case:
    """Create a time-reversal-symmetric case for sanity checks."""
    if M_seg < 2:
        raise ValueError("Symmetric case requires M_seg >= 2.")
    midpoint = np.linspace(0.0, 1.0, M_seg + 1)
    q = 1.0 - 4.0 * (midpoint - 0.5) ** 2
    q[0] = 0.0
    q[-1] = 0.0
    T = np.ones(M_seg, dtype=float)
    zeta_start, zeta_end = default_boundary_jets(q)
    return Phase2Case(case_name=case_name, q=q, T=T, zeta_start=zeta_start, zeta_end=zeta_end)


def build_waypoint_selector(M_seg: int, s: int = 4) -> np.ndarray:
    """
    Build the RHS selector S_q for interior waypoint variables.

    For the Phase 1 row ordering, each interior waypoint q_j appears in two
    interpolation rows: p_j(T_j) and p_{j+1}(0).
    """
    if M_seg < 1:
        raise ValueError("M_seg must be >= 1.")
    total_rows = 2 * s * M_seg
    num_interior = max(M_seg - 1, 0)
    S_q = np.zeros((total_rows, num_interior), dtype=float)
    for interior_index in range(num_interior):
        S_q[2 * interior_index + 1, interior_index] = 1.0
        S_q[2 * interior_index + 2, interior_index] = 1.0
    return S_q


def compute_exact_jacobian_q(M_mat: np.ndarray, S_q: np.ndarray) -> np.ndarray:
    """Compute J_q = A(T)^{-1} S_q without forming an explicit inverse."""
    M_mat = np.asarray(M_mat, dtype=float)
    S_q = np.asarray(S_q, dtype=float)
    if S_q.size == 0:
        return np.zeros((M_mat.shape[0], 0), dtype=float)
    return np.linalg.solve(M_mat, S_q)


def finite_difference_jacobian_q(
    q: np.ndarray,
    T: np.ndarray,
    zeta_start: np.ndarray,
    zeta_end: np.ndarray,
    eps: float = 1e-6,
) -> np.ndarray:
    """
    Estimate D_q c by perturbing interior waypoints one column at a time.
    """
    q = np.asarray(q, dtype=float).reshape(-1)
    T = np.asarray(T, dtype=float).reshape(-1)
    zeta_start = np.asarray(zeta_start, dtype=float).reshape(-1)
    zeta_end = np.asarray(zeta_end, dtype=float).reshape(-1)

    num_interior = max(q.size - 2, 0)
    base = solve_minco_coefficients(q, T, zeta_start, zeta_end)["c_vec"]
    if num_interior == 0:
        return np.zeros((base.size, 0), dtype=float)

    J_fd = np.zeros((base.size, num_interior), dtype=float)
    for column in range(num_interior):
        q_perturbed = q.copy()
        q_perturbed[column + 1] += eps
        perturbed = solve_minco_coefficients(q_perturbed, T, zeta_start, zeta_end)["c_vec"]
        J_fd[:, column] = (perturbed - base) / eps
    return J_fd


def compare_exact_vs_fd_jacobian(
    J_exact: np.ndarray,
    J_fd: np.ndarray,
) -> dict[str, Any]:
    """Compare exact and finite-difference Jacobians with scalar summary metrics."""
    J_exact = np.asarray(J_exact, dtype=float)
    J_fd = np.asarray(J_fd, dtype=float)
    diff = J_exact - J_fd
    exact_norm = float(np.linalg.norm(J_exact))
    fd_norm = float(np.linalg.norm(J_fd))
    diff_norm = float(np.linalg.norm(diff))
    denominator = max(exact_norm, 1e-15)
    per_column_error = (
        np.linalg.norm(diff, axis=0) if diff.ndim == 2 and diff.shape[1] > 0 else np.zeros((0,), dtype=float)
    )
    return {
        "fro_error": diff_norm,
        "max_abs_error": float(np.max(np.abs(diff))) if diff.size else 0.0,
        "relative_error": diff_norm / denominator,
        "exact_fro_norm": exact_norm,
        "fd_fro_norm": fd_norm,
        "per_column_error": per_column_error.tolist(),
    }


def compute_segmentwise_influence_norms(
    J_q: np.ndarray,
    M_seg: int,
    block_size: int = COEFFS_PER_SEGMENT,
) -> np.ndarray:
    """
    Compute block-level influence strengths I_ij = ||dc_i / dq_j||_2.
    """
    J_q = np.asarray(J_q, dtype=float)
    num_interior = max(M_seg - 1, 0)
    if J_q.shape != (block_size * M_seg, num_interior):
        raise ValueError(
            f"J_q must have shape ({block_size * M_seg}, {num_interior}), got {J_q.shape}."
        )
    influence = np.zeros((M_seg, num_interior), dtype=float)
    for segment_index in range(M_seg):
        block = J_q[segment_index * block_size : (segment_index + 1) * block_size, :]
        influence[segment_index] = np.linalg.norm(block, axis=0)
    return influence


def compute_waypoint_influence_profile(
    J_q: np.ndarray,
    M_seg: int,
    target_waypoint_idx: int,
    block_size: int = COEFFS_PER_SEGMENT,
) -> dict[str, Any]:
    """Return the segmentwise influence profile for one interior waypoint column."""
    influence = compute_segmentwise_influence_norms(J_q, M_seg, block_size=block_size)
    if not 0 <= target_waypoint_idx < influence.shape[1]:
        raise IndexError(
            f"target_waypoint_idx must be in [0, {max(influence.shape[1] - 1, 0)}], got {target_waypoint_idx}."
        )
    distances = np.abs(np.arange(M_seg) - target_waypoint_idx)
    norms = influence[:, target_waypoint_idx]
    return {
        "target_waypoint_idx": int(target_waypoint_idx),
        "segment_indices": np.arange(M_seg, dtype=int),
        "distances": distances,
        "influence_norms": norms,
        "max_influence": float(np.max(norms)) if norms.size else 0.0,
        "farthest_nonzero_distance": int(np.max(distances[norms > 0.0])) if np.any(norms > 0.0) else 0,
    }


def estimate_effective_bandwidth(
    J_q: np.ndarray,
    M_seg: int,
    tol: float,
    block_size: int = COEFFS_PER_SEGMENT,
) -> dict[str, Any]:
    """
    Estimate the apparent block bandwidth of J_q under a numerical threshold.
    """
    influence = compute_segmentwise_influence_norms(J_q, M_seg, block_size=block_size)
    num_interior = influence.shape[1]
    if num_interior == 0:
        return {
            "tol": float(tol),
            "max_effective_bandwidth": 0,
            "mean_effective_bandwidth": 0.0,
            "far_nonzero_ratio": 0.0,
            "per_waypoint_bandwidth": [],
            "influence_matrix": influence,
        }

    per_waypoint_bandwidth = []
    far_entries = 0
    far_entries_above_tol = 0
    max_far_influence = 0.0
    mean_far_accum = 0.0

    for waypoint_index in range(num_interior):
        distances = np.abs(np.arange(M_seg) - waypoint_index)
        significant = distances[influence[:, waypoint_index] > tol]
        bandwidth = int(np.max(significant)) if significant.size else 0
        per_waypoint_bandwidth.append(bandwidth)

        far_mask = distances > 0
        far_values = influence[far_mask, waypoint_index]
        far_entries += far_values.size
        far_entries_above_tol += int(np.sum(far_values > tol))
        if far_values.size:
            max_far_influence = max(max_far_influence, float(np.max(far_values)))
            mean_far_accum += float(np.sum(far_values))

    mean_far_influence = mean_far_accum / max(far_entries, 1)
    return {
        "tol": float(tol),
        "max_effective_bandwidth": int(max(per_waypoint_bandwidth, default=0)),
        "mean_effective_bandwidth": float(np.mean(per_waypoint_bandwidth)) if per_waypoint_bandwidth else 0.0,
        "far_nonzero_ratio": far_entries_above_tol / max(far_entries, 1),
        "max_far_influence": max_far_influence,
        "mean_far_influence": mean_far_influence,
        "num_far_entries_above_tol": int(far_entries_above_tol),
        "per_waypoint_bandwidth": per_waypoint_bandwidth,
        "influence_matrix": influence,
    }


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def run_phase2_validation_suite(
    q: np.ndarray,
    T: np.ndarray,
    zeta_start: np.ndarray,
    zeta_end: np.ndarray,
    *,
    case_name: str,
    eps: float = 1e-6,
    tol: float = DEFAULT_TOL,
    target_waypoint_idx: int | None = None,
    results_dir: str | Path = Path("phase_2/results"),
) -> dict[str, Any]:
    """
    Run the full Phase 2 validation flow for one case and save figures/tables/logs.
    """
    from phase_2.phase2_plotting import (
        visualize_block_sparsity_Jq,
        visualize_effective_bandwidth_vs_M,
        visualize_heatmap,
        visualize_jacobian_fd_compare,
        visualize_matrix_sparsity,
        visualize_waypoint_influence_profile,
    )

    directories = ensure_results_dirs(results_dir)
    q = np.asarray(q, dtype=float).reshape(-1)
    T = np.asarray(T, dtype=float).reshape(-1)
    zeta_start = np.asarray(zeta_start, dtype=float).reshape(-1)
    zeta_end = np.asarray(zeta_end, dtype=float).reshape(-1)
    M_seg = T.size

    system = solve_minco_coefficients(q, T, zeta_start, zeta_end, return_system=True)
    A_mat = system["M"]
    b_vec = system["b"]
    c_vec = system["c_vec"]
    S_q = build_waypoint_selector(M_seg)
    J_exact = compute_exact_jacobian_q(A_mat, S_q)
    J_fd = finite_difference_jacobian_q(q, T, zeta_start, zeta_end, eps=eps)
    jacobian_comparison = compare_exact_vs_fd_jacobian(J_exact, J_fd)
    influence_matrix = compute_segmentwise_influence_norms(J_exact, M_seg)
    if target_waypoint_idx is None and M_seg > 1:
        target_waypoint_idx = max((M_seg - 1) // 2 - 1, 0)
    elif target_waypoint_idx is None:
        target_waypoint_idx = 0
    profile = compute_waypoint_influence_profile(J_exact, M_seg, target_waypoint_idx)
    bandwidth = estimate_effective_bandwidth(J_exact, M_seg, tol=tol)

    figures_dir = directories["figures"]
    tables_dir = directories["tables"]
    logs_dir = directories["logs"]

    visualize_matrix_sparsity(
        A_mat,
        figures_dir / f"A_sparsity_M{M_seg}_case{case_name}.png",
        title=f"A(T) Sparsity, M={M_seg}, {case_name}",
    )
    visualize_heatmap(
        J_exact,
        figures_dir / f"Jq_heatmap_M{M_seg}_case{case_name}.png",
        title=f"J_q Heatmap, M={M_seg}, {case_name}",
        log_scale=True,
    )
    visualize_block_sparsity_Jq(
        J_exact,
        M_seg,
        figures_dir / f"block_influence_M{M_seg}_case{case_name}.png",
    )
    visualize_waypoint_influence_profile(
        profile,
        figures_dir / f"influence_profile_q{target_waypoint_idx + 1}_M{M_seg}_case{case_name}.png",
        title=f"Influence Profile for q_{target_waypoint_idx + 1}, M={M_seg}, {case_name}",
    )
    visualize_jacobian_fd_compare(
        J_exact,
        J_fd,
        figures_dir / f"jacobian_fd_compare_M{M_seg}_case{case_name}.png",
        title=f"Exact vs FD Jacobian, M={M_seg}, {case_name}",
    )

    jacobian_error_rows = [
        {
            "case_name": case_name,
            "M": M_seg,
            "eps": eps,
            "fro_error": jacobian_comparison["fro_error"],
            "max_abs_error": jacobian_comparison["max_abs_error"],
            "relative_error": jacobian_comparison["relative_error"],
        }
    ]
    _write_csv(
        tables_dir / "jacobian_error_summary.csv",
        ["case_name", "M", "eps", "fro_error", "max_abs_error", "relative_error"],
        jacobian_error_rows,
    )

    effective_bandwidth_rows = [
        {
            "case_name": case_name,
            "M": M_seg,
            "tol": tol,
            "max_effective_bandwidth": bandwidth["max_effective_bandwidth"],
            "mean_effective_bandwidth": bandwidth["mean_effective_bandwidth"],
            "far_nonzero_ratio": bandwidth["far_nonzero_ratio"],
        }
    ]
    _write_csv(
        tables_dir / "effective_bandwidth_summary.csv",
        ["case_name", "M", "tol", "max_effective_bandwidth", "mean_effective_bandwidth", "far_nonzero_ratio"],
        effective_bandwidth_rows,
    )

    far_field_rows = []
    for waypoint_index in range(max(M_seg - 1, 0)):
        profile_row = compute_waypoint_influence_profile(J_exact, M_seg, waypoint_index)
        distances = profile_row["distances"]
        norms = profile_row["influence_norms"]
        far_mask = distances > 1
        far_norms = norms[far_mask]
        far_field_rows.append(
            {
                "case_name": case_name,
                "M": M_seg,
                "target_waypoint_idx": waypoint_index,
                "max_far_influence": float(np.max(far_norms)) if far_norms.size else 0.0,
                "mean_far_influence": float(np.mean(far_norms)) if far_norms.size else 0.0,
                "num_far_entries_above_tol": int(np.sum(far_norms > tol)),
            }
        )
    _write_csv(
        tables_dir / "far_field_influence_summary.csv",
        [
            "case_name",
            "M",
            "target_waypoint_idx",
            "max_far_influence",
            "mean_far_influence",
            "num_far_entries_above_tol",
        ],
        far_field_rows,
    )

    log_payload = {
        "case_name": case_name,
        "M": M_seg,
        "eps": eps,
        "tol": tol,
        "system_residual_norm": float(np.linalg.norm(A_mat @ c_vec - b_vec)),
        "jacobian_comparison": jacobian_comparison,
        "bandwidth": {
            key: value
            for key, value in bandwidth.items()
            if key != "influence_matrix"
        },
        "target_waypoint_idx": target_waypoint_idx,
    }
    _write_json(logs_dir / f"phase2_validation_M{M_seg}_case{case_name}.json", log_payload)

    return {
        "case_name": case_name,
        "M": M_seg,
        "A_mat": A_mat,
        "b_vec": b_vec,
        "c_vec": c_vec,
        "coeffs": system["coeffs"],
        "S_q": S_q,
        "J_q": J_exact,
        "J_fd": J_fd,
        "jacobian_comparison": jacobian_comparison,
        "influence_matrix": influence_matrix,
        "profile": profile,
        "bandwidth": bandwidth,
        "results_dir": directories,
    }


def run_scaling_experiment(
    M_values: tuple[int, ...] = (4, 8, 16, 32),
    *,
    case_name: str = "uniform_time",
    tol: float = DEFAULT_TOL,
    eps: float = 1e-6,
    results_dir: str | Path = Path("phase_2/results"),
) -> dict[str, Any]:
    """Run the effective-bandwidth experiment as M grows and save summary artifacts."""
    from phase_2.phase2_plotting import visualize_effective_bandwidth_vs_M

    directories = ensure_results_dirs(results_dir)
    records = []
    for M_seg in M_values:
        if case_name == "uniform_time":
            case = make_uniform_time_case(M_seg, case_name=case_name)
        elif case_name == "nonuniform_time":
            case = make_nonuniform_time_case(M_seg, seed=DEFAULT_SEED, case_name=case_name)
        elif case_name == "symmetric":
            case = make_symmetric_case(M_seg, case_name=case_name)
        else:
            raise ValueError(f"Unsupported scaling case_name {case_name!r}.")

        system = solve_minco_coefficients(
            case.q,
            case.T,
            case.zeta_start,
            case.zeta_end,
            return_system=True,
        )
        S_q = build_waypoint_selector(M_seg)
        J_q = compute_exact_jacobian_q(system["M"], S_q)
        bandwidth = estimate_effective_bandwidth(J_q, M_seg, tol=tol)
        J_fd = finite_difference_jacobian_q(case.q, case.T, case.zeta_start, case.zeta_end, eps=eps)
        comparison = compare_exact_vs_fd_jacobian(J_q, J_fd)

        records.append(
            {
                "case_name": case_name,
                "M": M_seg,
                "tol": tol,
                "max_effective_bandwidth": bandwidth["max_effective_bandwidth"],
                "mean_effective_bandwidth": bandwidth["mean_effective_bandwidth"],
                "far_nonzero_ratio": bandwidth["far_nonzero_ratio"],
                "fro_error": comparison["fro_error"],
                "max_abs_error": comparison["max_abs_error"],
                "relative_error": comparison["relative_error"],
            }
        )

    _write_csv(
        directories["tables"] / "effective_bandwidth_summary.csv",
        [
            "case_name",
            "M",
            "tol",
            "max_effective_bandwidth",
            "mean_effective_bandwidth",
            "far_nonzero_ratio",
            "fro_error",
            "max_abs_error",
            "relative_error",
        ],
        records,
    )
    visualize_effective_bandwidth_vs_M(
        records,
        directories["figures"] / f"effective_bandwidth_vs_M_case{case_name}.png",
        title=f"Effective Bandwidth vs M ({case_name})",
    )
    _write_json(
        directories["logs"] / f"effective_bandwidth_scaling_case{case_name}.json",
        {"records": records, "tol": tol, "eps": eps},
    )
    return {"records": records, "results_dir": directories}
