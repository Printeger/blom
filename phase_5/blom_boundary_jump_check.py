"""
Phase 5 global assembly validation for BLOM.

This module implements three assembly schemes on top of the Phase 1/3/4
building blocks and evaluates their boundary derivative jumps:

- Scheme A: shared junction states
- Scheme B: overlapping consensus
- Scheme C: raw central-segment extraction
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Any, Callable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from phase_3.blom_strict_local_qp import build_local_problem, extract_segment_coeff, solve_blom_strict_local_qp
from phase_4.blom_k2_s4_numeric import middle_segment_cost_matrix
from phase_4.utils.hermite_utils import normalized_derivative_row


DEFAULT_S = 4
DEFAULT_K = 2
DEFAULT_RESULTS_DIR = Path("phase_5/results/phase5_boundary_jump_check")
DEFAULT_ZERO_TOL = 1e-10


def _validate_time_vector(T: np.ndarray) -> np.ndarray:
    T = np.asarray(T, dtype=float).reshape(-1)
    if T.size < 1:
        raise ValueError("T must describe at least one segment.")
    if not np.all(np.isfinite(T)):
        raise ValueError("T contains non-finite values.")
    if not np.all(T > 0.0):
        raise ValueError("All segment durations must be strictly positive.")
    return T


def _validate_problem_inputs(q: np.ndarray, T: np.ndarray, s: int, k: int) -> tuple[np.ndarray, np.ndarray]:
    q = np.asarray(q, dtype=float).reshape(-1)
    T = _validate_time_vector(T)
    if s < 2:
        raise ValueError("s must be >= 2.")
    if k < 0:
        raise ValueError("k must be >= 0.")
    if q.shape != (T.size + 1,):
        raise ValueError(f"q must have shape ({T.size + 1},), got {q.shape}.")
    if not np.all(np.isfinite(q)):
        raise ValueError("q contains non-finite values.")
    return q, T


def default_boundary_jets(q: np.ndarray, s: int = DEFAULT_S) -> tuple[np.ndarray, np.ndarray]:
    q = np.asarray(q, dtype=float).reshape(-1)
    zeta_start = np.zeros((s,), dtype=float)
    zeta_end = np.zeros((s,), dtype=float)
    zeta_start[0] = float(q[0])
    zeta_end[0] = float(q[-1])
    return zeta_start, zeta_end


def _ensure_dir(path: str | Path | None) -> Path | None:
    if path is None:
        return None
    target = Path(path)
    target.mkdir(parents=True, exist_ok=True)
    return target


def _save_json(path: str | Path, payload: dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _save_csv(path: str | Path, header: list[str], rows: list[list[Any]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [",".join(header)]
    for row in rows:
        rendered = []
        for value in row:
            if isinstance(value, bool):
                rendered.append("true" if value else "false")
            else:
                rendered.append(str(value))
        lines.append(",".join(rendered))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _serialize(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, dict):
        return {str(key): _serialize(val) for key, val in value.items()}
    if isinstance(value, list):
        return [_serialize(item) for item in value]
    return value


def _poly_derivative_value(coeff: np.ndarray, t: float, order: int) -> float:
    coeff = np.asarray(coeff, dtype=float).reshape(-1)
    degree = coeff.size - 1
    row = normalized_derivative_row(float(t), order, degree)
    return float(coeff @ row)


def compute_jumps(
    coeffs: np.ndarray,
    T: np.ndarray,
    s: int,
    max_order: int | None = None,
) -> dict[int, np.ndarray]:
    """Compute derivative jumps across all interior knots up to `max_order`."""
    coeffs = np.asarray(coeffs, dtype=float)
    T = _validate_time_vector(T)
    if coeffs.shape != (T.size, 2 * s):
        raise ValueError(f"coeffs must have shape ({T.size}, {2 * s}), got {coeffs.shape}.")
    if max_order is None:
        max_order = 2 * s - 2
    if max_order < 0:
        raise ValueError("max_order must be >= 0.")

    jumps: dict[int, np.ndarray] = {}
    for order in range(max_order + 1):
        values = np.zeros((max(T.size - 1, 0),), dtype=float)
        for knot_index in range(T.size - 1):
            left = _poly_derivative_value(coeffs[knot_index], float(T[knot_index]), order)
            right = _poly_derivative_value(coeffs[knot_index + 1], 0.0, order)
            values[knot_index] = left - right
        jumps[order] = values
    return jumps


def summarize_jumps(jumps: dict[int, np.ndarray], tol: float = DEFAULT_ZERO_TOL) -> dict[int, dict[str, Any]]:
    """Summarize jump statistics for every derivative order."""
    summary: dict[int, dict[str, Any]] = {}
    for order, values in sorted(jumps.items()):
        values = np.asarray(values, dtype=float).reshape(-1)
        abs_values = np.abs(values)
        if values.size == 0:
            summary[order] = {
                "max_abs": 0.0,
                "mean_abs": 0.0,
                "median_abs": 0.0,
                "rms": 0.0,
                "q95_abs": 0.0,
                "is_zero_tol": True,
            }
            continue
        summary[order] = {
            "max_abs": float(np.max(abs_values)),
            "mean_abs": float(np.mean(abs_values)),
            "median_abs": float(np.median(abs_values)),
            "rms": float(np.sqrt(np.mean(values**2))),
            "q95_abs": float(np.quantile(abs_values, 0.95)),
            "is_zero_tol": bool(np.max(abs_values) <= tol),
        }
    return summary


def _jump_matrix(jumps: dict[int, np.ndarray]) -> np.ndarray:
    if not jumps:
        return np.zeros((0, 0), dtype=float)
    orders = sorted(jumps)
    if jumps[orders[0]].size == 0:
        return np.zeros((len(orders), 0), dtype=float)
    return np.vstack([np.asarray(jumps[order], dtype=float) for order in orders])


def _plot_jump_heatmap(jumps: dict[int, np.ndarray], save_path: str | Path, title: str) -> None:
    mat = _jump_matrix(jumps)
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    display = np.log10(np.abs(mat) + 1e-16) if mat.size else np.zeros((1, 1), dtype=float)
    image = ax.imshow(display, aspect="auto", origin="lower", cmap="magma")
    ax.set_title(title)
    ax.set_xlabel("Interior Knot Index")
    ax.set_ylabel("Derivative Order")
    if mat.shape[1] > 0:
        ax.set_xticks(np.arange(mat.shape[1]))
        ax.set_xticklabels([str(index + 1) for index in range(mat.shape[1])])
    if mat.shape[0] > 0:
        ax.set_yticks(np.arange(mat.shape[0]))
        ax.set_yticklabels([str(order) for order in sorted(jumps)])
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("log10(|jump| + 1e-16)")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_jump_maxbar(stats: dict[int, dict[str, Any]], save_path: str | Path, title: str) -> None:
    orders = np.asarray(sorted(stats), dtype=int)
    max_vals = np.asarray([stats[int(order)]["max_abs"] for order in orders], dtype=float)
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.bar(orders.astype(str), max_vals, color="#1f77b4")
    ax.set_title(title)
    ax.set_xlabel("Derivative Order")
    ax.set_ylabel("max_i |jump|")
    ax.set_yscale("log")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_lower_order_compare(
    results: dict[str, dict[str, Any]],
    s: int,
    save_path: str | Path,
) -> None:
    orders = np.arange(s, dtype=int)
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for scheme, result in sorted(results.items()):
        values = [result["stats"][int(order)]["max_abs"] for order in orders]
        ax.plot(orders, values, marker="o", linewidth=1.8, label=f"Scheme {scheme}")
    ax.set_title("Lower-Order Jump Comparison")
    ax.set_xlabel("Derivative Order")
    ax.set_ylabel("max_i |jump|")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_eta_mismatch(norms: np.ndarray, save_path: str | Path, title: str) -> None:
    norms = np.asarray(norms, dtype=float).reshape(-1)
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4))
    if norms.size:
        ax.plot(np.arange(1, norms.size + 1), norms, marker="o", linewidth=1.8)
    ax.set_title(title)
    ax.set_xlabel("Interior Knot Index")
    ax.set_ylabel("||eta^- - eta^+||_2")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_consensus_improvement(
    dispersion: np.ndarray,
    post_jump: np.ndarray,
    save_path: str | Path,
) -> None:
    dispersion = np.asarray(dispersion, dtype=float).reshape(-1)
    post_jump = np.asarray(post_jump, dtype=float).reshape(-1)
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4))
    if dispersion.size:
        knots = np.arange(1, dispersion.size + 1)
        ax.plot(knots, dispersion, marker="o", linewidth=1.8, label="Pre-consensus dispersion")
        ax.plot(knots, post_jump, marker="s", linewidth=1.8, label="Post-consensus ||jump_1:3||_2")
    ax.set_title("Scheme B Consensus Improvement")
    ax.set_xlabel("Interior Knot Index")
    ax.set_ylabel("Mismatch Magnitude")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_random_boxplot(records: list[dict[str, Any]], s: int, save_path: str | Path) -> None:
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    schemes = sorted({record["scheme"] for record in records})
    lower_data = [
        [record["lower_order_max"] for record in records if record["scheme"] == scheme] for scheme in schemes
    ]
    higher_data = [
        [record["higher_order_max"] for record in records if record["scheme"] == scheme] for scheme in schemes
    ]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    axes[0].boxplot(lower_data, tick_labels=schemes)
    axes[0].set_title(f"Lower-Order Jump Max, orders 0..{s - 1}")
    axes[0].set_ylabel("max jump")
    axes[0].set_yscale("log")
    axes[0].grid(True, axis="y", alpha=0.3)

    axes[1].boxplot(higher_data, tick_labels=schemes)
    axes[1].set_title(f"Higher-Order Jump Max, orders {s}..{2 * s - 2}")
    axes[1].set_ylabel("max jump")
    axes[1].set_yscale("log")
    axes[1].grid(True, axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _eta_slice(knot_index: int, s: int) -> slice:
    width = s - 1
    start = (knot_index - 1) * width
    return slice(start, start + width)


def _segment_state_maps(
    segment_index: int,
    M: int,
    s: int,
    zeta_start: np.ndarray,
    zeta_end: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    n_eta = max(M - 1, 0) * (s - 1)
    local_A = np.zeros((2 * (s - 1), n_eta), dtype=float)
    local_b = np.zeros((2 * (s - 1),), dtype=float)
    if segment_index == 1:
        local_b[: s - 1] = zeta_start[1:s]
    else:
        local_A[: s - 1, _eta_slice(segment_index - 1, s)] = np.eye(s - 1)
    if segment_index == M:
        local_b[s - 1 :] = zeta_end[1:s]
    else:
        local_A[s - 1 :, _eta_slice(segment_index, s)] = np.eye(s - 1)
    return local_A, local_b


def _segment_quadratic(
    q_left: float,
    q_right: float,
    duration: float,
) -> tuple[np.ndarray, np.ndarray, float]:
    if duration <= 0.0:
        raise ValueError("Segment duration must be positive.")
    if DEFAULT_S != 4:
        raise RuntimeError("Phase 5 canonical segment quadratic assumes s=4.")
    middle = middle_segment_cost_matrix(q_left, q_right, duration)
    H = middle["H_mid"]
    B = middle["g_mid"] @ np.asarray([q_left, q_right], dtype=float)
    C = float(middle["constant"])
    return H, B, C


def _reconstruct_segment_coeff(
    q_left: float,
    q_right: float,
    duration: float,
    eta_left: np.ndarray,
    eta_right: np.ndarray,
    s: int,
) -> np.ndarray:
    degree = 2 * s - 1
    rows = [normalized_derivative_row(0.0, order, degree) for order in range(s)]
    rows.extend(normalized_derivative_row(1.0, order, degree) for order in range(s))
    C_s = np.vstack(rows)
    H_s = np.linalg.inv(C_s)
    duration_powers = np.asarray([duration**order for order in range(1, s)], dtype=float)
    y = np.concatenate(
        (
            [float(q_left)],
            duration_powers * np.asarray(eta_left, dtype=float).reshape(s - 1),
            [float(q_right)],
            duration_powers * np.asarray(eta_right, dtype=float).reshape(s - 1),
        )
    )
    alpha = H_s @ y
    scaling = np.diag([duration ** (-power) for power in range(2 * s)])
    return scaling @ alpha


def _extract_eta_from_coeff(
    coeff: np.ndarray,
    duration: float,
    s: int,
    side: str,
) -> np.ndarray:
    t = 0.0 if side == "left" else float(duration)
    return np.asarray([_poly_derivative_value(coeff, t, order) for order in range(1, s)], dtype=float)


def _build_scheme_a_system(
    q: np.ndarray,
    T: np.ndarray,
    s: int,
    zeta_start: np.ndarray,
    zeta_end: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    M = T.size
    n_eta = max(M - 1, 0) * (s - 1)
    H_global = np.zeros((n_eta, n_eta), dtype=float)
    g_global = np.zeros((n_eta,), dtype=float)
    for segment_index in range(1, M + 1):
        local_A, local_b = _segment_state_maps(segment_index, M, s, zeta_start, zeta_end)
        H_seg, B_seg, _ = _segment_quadratic(q[segment_index - 1], q[segment_index], float(T[segment_index - 1]))
        H_global += local_A.T @ H_seg @ local_A
        g_global += local_A.T @ (B_seg - H_seg @ local_b)
    return H_global, g_global


def _reconstruct_from_eta(
    q: np.ndarray,
    T: np.ndarray,
    eta_shared: np.ndarray,
    s: int,
    zeta_start: np.ndarray,
    zeta_end: np.ndarray,
) -> np.ndarray:
    M = T.size
    coeffs = np.zeros((M, 2 * s), dtype=float)
    for segment_index in range(1, M + 1):
        if segment_index == 1:
            eta_left = zeta_start[1:s]
        else:
            eta_left = eta_shared[segment_index - 2]
        if segment_index == M:
            eta_right = zeta_end[1:s]
        else:
            eta_right = eta_shared[segment_index - 1]
        coeffs[segment_index - 1] = _reconstruct_segment_coeff(
            q[segment_index - 1],
            q[segment_index],
            float(T[segment_index - 1]),
            eta_left,
            eta_right,
            s,
        )
    return coeffs


def _format_stats_for_csv(stats: dict[int, dict[str, Any]]) -> list[list[Any]]:
    rows = []
    for order in sorted(stats):
        record = stats[order]
        rows.append(
            [
                order,
                record["max_abs"],
                record["mean_abs"],
                record["median_abs"],
                record["rms"],
                record["q95_abs"],
                record["is_zero_tol"],
            ]
        )
    return rows


def assemble_scheme_A(
    q: np.ndarray,
    T: np.ndarray,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Assemble the trajectory via globally shared interior junction states."""
    del k
    q, T = _validate_problem_inputs(q, T, s=s, k=DEFAULT_K)
    config = {} if config is None else dict(config)
    zeta_start = np.asarray(config.get("zeta_start", default_boundary_jets(q, s)[0]), dtype=float).reshape(s)
    zeta_end = np.asarray(config.get("zeta_end", default_boundary_jets(q, s)[1]), dtype=float).reshape(s)
    start_time = time.perf_counter()

    n_eta = max(T.size - 1, 0) * (s - 1)
    if n_eta == 0:
        eta_vec = np.zeros((0,), dtype=float)
        H_global = np.zeros((0, 0), dtype=float)
        g_global = np.zeros((0,), dtype=float)
    else:
        H_global, g_global = _build_scheme_a_system(q, T, s, zeta_start, zeta_end)
        eta_vec = np.linalg.solve(H_global, g_global)
    eta_shared = eta_vec.reshape(max(T.size - 1, 0), s - 1)
    coeffs = _reconstruct_from_eta(q, T, eta_shared, s, zeta_start, zeta_end)
    jumps = compute_jumps(coeffs, T, s)
    stats = summarize_jumps(jumps, tol=float(config.get("zero_tol", DEFAULT_ZERO_TOL)))
    elapsed = time.perf_counter() - start_time

    lower_max = max((stats[order]["max_abs"] for order in range(s)), default=0.0)
    higher_orders = range(s, 2 * s - 1)
    higher_max = max((stats[order]["max_abs"] for order in higher_orders), default=0.0)
    result = {
        "scheme": "A",
        "eta_shared": eta_shared,
        "coeffs": coeffs,
        "jumps": jumps,
        "stats": stats,
        "meta": {
            "q": q,
            "T": T,
            "s": s,
            "k": DEFAULT_K,
            "elapsed_sec": elapsed,
            "success": True,
        },
        "system": {
            "H_global": H_global,
            "g_global": g_global,
            "symmetry_error": float(np.max(np.abs(H_global - H_global.T))) if H_global.size else 0.0,
            "condition_number": float(np.linalg.cond(H_global)) if H_global.size else 1.0,
        },
        "stats_overview": {
            "lower_order_max": lower_max,
            "higher_order_max": higher_max,
            "reaches_Cs_minus_1": bool(all(stats[order]["is_zero_tol"] for order in range(s))),
            "needs_global_variables": True,
            "needs_consensus": False,
            "locality_score": "low",
        },
    }
    return result


def _collect_window_solution(
    q: np.ndarray,
    T: np.ndarray,
    center_segment: int,
    s: int,
    k: int,
    config: dict[str, Any],
    zeta_start: np.ndarray,
    zeta_end: np.ndarray,
) -> dict[str, Any]:
    problem = build_local_problem(q, T, i=center_segment, k=k, zeta_start=zeta_start, zeta_end=zeta_end, s=s)
    method = str(config.get("local_solver", "kkt"))
    solution = solve_blom_strict_local_qp(problem, method=method)
    return solution


def assemble_scheme_C(
    q: np.ndarray,
    T: np.ndarray,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Assemble the trajectory by raw central-segment extraction only."""
    q, T = _validate_problem_inputs(q, T, s=s, k=k)
    config = {} if config is None else dict(config)
    zeta_start, zeta_end = default_boundary_jets(q, s)
    start_time = time.perf_counter()

    coeffs = np.zeros((T.size, 2 * s), dtype=float)
    local_summaries: list[dict[str, Any]] = []
    for seg_idx in range(1, T.size + 1):
        solution = _collect_window_solution(q, T, seg_idx, s, k, config, zeta_start, zeta_end)
        coeffs[seg_idx - 1] = extract_segment_coeff(solution, seg_idx)
        local_summaries.append(
            {
                "center_segment": seg_idx,
                "window_type": solution["problem"]["window"]["window_type"],
                "objective": float(solution["objective"]),
            }
        )

    eta_left = np.vstack([_extract_eta_from_coeff(coeffs[idx], T[idx], s, "right") for idx in range(T.size - 1)])
    eta_right = np.vstack([_extract_eta_from_coeff(coeffs[idx + 1], T[idx + 1], s, "left") for idx in range(T.size - 1)])
    eta_delta = eta_left - eta_right if T.size > 1 else np.zeros((0, s - 1), dtype=float)

    jumps = compute_jumps(coeffs, T, s)
    stats = summarize_jumps(jumps, tol=float(config.get("zero_tol", DEFAULT_ZERO_TOL)))
    elapsed = time.perf_counter() - start_time
    lower_orders = [stats[order]["max_abs"] for order in range(1, s)] if s > 1 else []
    higher_orders = [stats[order]["max_abs"] for order in range(s, 2 * s - 1)]

    return {
        "scheme": "C",
        "coeffs": coeffs,
        "eta_left": eta_left,
        "eta_right": eta_right,
        "eta_delta": eta_delta,
        "jumps": jumps,
        "stats": stats,
        "local_windows": local_summaries,
        "meta": {
            "q": q,
            "T": T,
            "s": s,
            "k": k,
            "elapsed_sec": elapsed,
            "success": True,
        },
        "stats_overview": {
            "lower_order_max": max([stats[0]["max_abs"], *lower_orders], default=0.0),
            "higher_order_max": max(higher_orders, default=0.0),
            "reaches_Cs_minus_1": bool(all(stats[order]["is_zero_tol"] for order in range(s))),
            "needs_global_variables": False,
            "needs_consensus": False,
            "locality_score": "high",
            "eta_mismatch_norms": np.linalg.norm(eta_delta, axis=1) if eta_delta.size else np.zeros((0,), dtype=float),
        },
    }


def assemble_scheme_B(
    q: np.ndarray,
    T: np.ndarray,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Assemble the trajectory via local predictions followed by knot consensus."""
    q, T = _validate_problem_inputs(q, T, s=s, k=k)
    config = {} if config is None else dict(config)
    zeta_start, zeta_end = default_boundary_jets(q, s)
    start_time = time.perf_counter()

    prediction_lists: dict[int, list[dict[str, Any]]] = {knot_index: [] for knot_index in range(1, T.size)}
    local_windows: list[dict[str, Any]] = []
    for center_segment in range(1, T.size + 1):
        solution = _collect_window_solution(q, T, center_segment, s, k, config, zeta_start, zeta_end)
        window = solution["problem"]["window"]
        local_windows.append(
            {
                "center_segment": center_segment,
                "segments": list(window["segments"]),
                "window_type": window["window_type"],
                "objective": float(solution["objective"]),
            }
        )
        for knot_index in range(window["L"], window["R"]):
            if not (1 <= knot_index <= T.size - 1):
                continue
            coeff_left = extract_segment_coeff(solution, knot_index)
            eta = _extract_eta_from_coeff(coeff_left, T[knot_index - 1], s, "right")
            prediction_lists[knot_index].append(
                {
                    "window_center": center_segment,
                    "eta": eta,
                    "weight": float(config.get("consensus_weight", 1.0)),
                }
            )

    eta_consensus = np.zeros((max(T.size - 1, 0), s - 1), dtype=float)
    counts = np.zeros((max(T.size - 1, 0),), dtype=int)
    dispersion = np.zeros((max(T.size - 1, 0),), dtype=float)
    weights_dict: dict[int, list[float]] = {}
    for knot_index in range(1, T.size):
        predictions = prediction_lists[knot_index]
        if not predictions:
            raise RuntimeError(f"No local predictions collected for interior knot {knot_index}.")
        weights = np.asarray([entry["weight"] for entry in predictions], dtype=float)
        values = np.vstack([entry["eta"] for entry in predictions])
        eta_mean = np.average(values, axis=0, weights=weights)
        eta_consensus[knot_index - 1] = eta_mean
        counts[knot_index - 1] = len(predictions)
        dispersion[knot_index - 1] = float(np.max(np.linalg.norm(values - eta_mean[None, :], axis=1)))
        weights_dict[knot_index] = weights.tolist()

    coeffs = _reconstruct_from_eta(q, T, eta_consensus, s, zeta_start, zeta_end)
    jumps = compute_jumps(coeffs, T, s)
    stats = summarize_jumps(jumps, tol=float(config.get("zero_tol", DEFAULT_ZERO_TOL)))
    post_lower_norm = np.zeros((max(T.size - 1, 0),), dtype=float)
    for knot_index in range(1, T.size):
        components = np.asarray([jumps[order][knot_index - 1] for order in range(1, s)], dtype=float)
        post_lower_norm[knot_index - 1] = float(np.linalg.norm(components))
    elapsed = time.perf_counter() - start_time

    return {
        "scheme": "B",
        "eta_predictions": {
            knot_index: [
                {
                    "window_center": entry["window_center"],
                    "eta": entry["eta"],
                    "weight": entry["weight"],
                }
                for entry in predictions
            ]
            for knot_index, predictions in prediction_lists.items()
        },
        "eta_consensus": eta_consensus,
        "weights": weights_dict,
        "coeffs": coeffs,
        "jumps": jumps,
        "stats": stats,
        "local_windows": local_windows,
        "consensus": {
            "prediction_counts": counts,
            "pre_dispersion": dispersion,
            "post_lower_order_jump_norm": post_lower_norm,
        },
        "meta": {
            "q": q,
            "T": T,
            "s": s,
            "k": k,
            "elapsed_sec": elapsed,
            "success": True,
        },
        "stats_overview": {
            "lower_order_max": max((stats[order]["max_abs"] for order in range(s)), default=0.0),
            "higher_order_max": max((stats[order]["max_abs"] for order in range(s, 2 * s - 1)), default=0.0),
            "reaches_Cs_minus_1": bool(all(stats[order]["is_zero_tol"] for order in range(s))),
            "needs_global_variables": False,
            "needs_consensus": True,
            "locality_score": "medium-high",
        },
    }


def _write_interpretation_summary(
    compare_result: dict[str, Any],
    save_path: str | Path,
    s: int,
) -> None:
    scheme_a = compare_result["scheme_results"]["A"]
    scheme_b = compare_result["scheme_results"]["B"]
    scheme_c = compare_result["scheme_results"]["C"]

    def lower_zero(result: dict[str, Any]) -> bool:
        return all(result["stats"][order]["is_zero_tol"] for order in range(s))

    lines = [
        "# Phase 5 Interpretation Summary",
        "",
        "## Scheme A",
        f"- lower-order jumps reach machine tolerance: `{lower_zero(scheme_a)}`",
        "- requires extra global shared knot-state variables",
        f"- max lower-order jump: `{scheme_a['stats_overview']['lower_order_max']:.3e}`",
        f"- max higher-order jump: `{scheme_a['stats_overview']['higher_order_max']:.3e}`",
        "",
        "## Scheme B",
        f"- lower-order jumps reach machine tolerance: `{lower_zero(scheme_b)}`",
        f"- consensus max pre-dispersion: `{np.max(scheme_b['consensus']['pre_dispersion']) if scheme_b['consensus']['pre_dispersion'].size else 0.0:.3e}`",
        f"- consensus max post-jump norm: `{np.max(scheme_b['consensus']['post_lower_order_jump_norm']) if scheme_b['consensus']['post_lower_order_jump_norm'].size else 0.0:.3e}`",
        "- preserves independent sliding-window solves and then projects to a shared knot state",
        "",
        "## Scheme C",
        f"- exact C0 only: `order 0 is_zero_tol = {scheme_c['stats'][0]['is_zero_tol']}`",
        f"- any lower-order derivative beyond position guaranteed zero: `{all(scheme_c['stats'][order]['is_zero_tol'] for order in range(1, s))}`",
        f"- max ||eta^- - eta^+||_2: `{np.max(scheme_c['stats_overview']['eta_mismatch_norms']) if scheme_c['stats_overview']['eta_mismatch_norms'].size else 0.0:.3e}`",
        "- keeps the cleanest local-map interpretation but does not automatically guarantee C^{s-1}",
        "",
        "## Recommendation",
        "- Scheme B is the safest assembly choice when exact global C^{s-1} continuity and local-window independence are both desired.",
        "- Scheme A is also smooth up to order s-1, but it introduces a true global shared-state optimization layer.",
        "- Scheme C remains useful as a diagnostic local-map baseline, not as the default globally smooth assembly mechanism.",
    ]
    Path(save_path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_boundary_jump_check(
    q: np.ndarray,
    T: np.ndarray,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    scheme: str = "A",
    config: dict[str, Any] | None = None,
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    """Run one assembly scheme and optionally save its artifacts."""
    del seed
    scheme = str(scheme).upper()
    if scheme not in {"A", "B", "C"}:
        raise ValueError(f"Unsupported scheme {scheme!r}.")
    assemblers: dict[str, Callable[..., dict[str, Any]]] = {
        "A": assemble_scheme_A,
        "B": assemble_scheme_B,
        "C": assemble_scheme_C,
    }
    result = assemblers[scheme](q, T, s=s, k=k, config=config)
    result["figure_paths"] = {}
    result["table_paths"] = {}

    if save_dir is None:
        return result

    scheme_dir = _ensure_dir(Path(save_dir) / f"scheme_{scheme}")
    assert scheme_dir is not None
    heatmap_path = scheme_dir / f"jump_heatmap_scheme_{scheme}.png"
    maxbar_path = scheme_dir / f"jump_maxbar_scheme_{scheme}.png"
    summary_json_path = scheme_dir / f"summary_scheme_{scheme}.json"
    stats_csv_path = scheme_dir / f"jump_stats_scheme_{scheme}.csv"

    _plot_jump_heatmap(result["jumps"], heatmap_path, f"Scheme {scheme} Jump Heatmap")
    _plot_jump_maxbar(result["stats"], maxbar_path, f"Scheme {scheme} Jump Maxima")
    _save_json(summary_json_path, _serialize(result))
    _save_csv(
        stats_csv_path,
        ["order", "max_abs", "mean_abs", "median_abs", "rms", "q95_abs", "is_zero_tol"],
        _format_stats_for_csv(result["stats"]),
    )
    result["figure_paths"]["heatmap"] = str(heatmap_path)
    result["figure_paths"]["maxbar"] = str(maxbar_path)
    result["table_paths"]["summary_json"] = str(summary_json_path)
    result["table_paths"]["jump_stats_csv"] = str(stats_csv_path)

    if scheme == "C":
        mismatch_path = scheme_dir / "scheme_C_eta_mismatch.png"
        mismatch_norms = result["stats_overview"]["eta_mismatch_norms"]
        _plot_eta_mismatch(mismatch_norms, mismatch_path, "Scheme C Eta Mismatch")
        result["figure_paths"]["eta_mismatch"] = str(mismatch_path)
    if scheme == "B":
        consensus_path = scheme_dir / "scheme_B_consensus_improvement.png"
        _plot_consensus_improvement(
            result["consensus"]["pre_dispersion"],
            result["consensus"]["post_lower_order_jump_norm"],
            consensus_path,
        )
        result["figure_paths"]["consensus_improvement"] = str(consensus_path)
    return result


def run_compare_all_schemes(
    q: np.ndarray,
    T: np.ndarray,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    config: dict[str, Any] | None = None,
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    """Run schemes A/B/C on one case and save all comparison artifacts."""
    q, T = _validate_problem_inputs(q, T, s=s, k=k)
    if save_dir is None:
        save_dir = DEFAULT_RESULTS_DIR
    save_root = _ensure_dir(save_dir)
    assert save_root is not None

    results = {
        scheme: run_boundary_jump_check(q, T, s=s, k=k, scheme=scheme, config=config, save_dir=save_root, seed=seed)
        for scheme in ("A", "B", "C")
    }

    lower_compare_path = save_root / "jump_lower_orders_compare.png"
    _plot_lower_order_compare(results, s=s, save_path=lower_compare_path)

    scheme_c_mismatch_path = save_root / "scheme_C_eta_mismatch.png"
    _plot_eta_mismatch(
        results["C"]["stats_overview"]["eta_mismatch_norms"],
        scheme_c_mismatch_path,
        "Scheme C Eta Mismatch",
    )

    scheme_b_consensus_path = save_root / "scheme_B_consensus_improvement.png"
    _plot_consensus_improvement(
        results["B"]["consensus"]["pre_dispersion"],
        results["B"]["consensus"]["post_lower_order_jump_norm"],
        scheme_b_consensus_path,
    )

    comparison_rows = []
    for scheme in ("A", "B", "C"):
        overview = results[scheme]["stats_overview"]
        comparison_rows.append(
            [
                scheme,
                overview["lower_order_max"],
                max((results[scheme]["stats"][order]["rms"] for order in range(s)), default=0.0),
                overview["higher_order_max"],
                overview["reaches_Cs_minus_1"],
                overview["needs_global_variables"],
                overview["needs_consensus"],
                overview["locality_score"],
                results[scheme]["meta"]["elapsed_sec"],
            ]
        )
    comparison_csv_path = save_root / "scheme_comparison_summary.csv"
    _save_csv(
        comparison_csv_path,
        [
            "scheme",
            "max_jump_lower_orders",
            "rms_jump_lower_orders",
            "max_jump_higher_orders",
            "reaches_Cs_minus_1",
            "needs_global_variables",
            "needs_consensus",
            "locality_score",
            "elapsed_sec",
        ],
        comparison_rows,
    )

    interpretation_path = save_root / "phase5_interpretation_summary.md"
    compare_result = {
        "inputs": {"q": q, "T": T, "s": s, "k": k},
        "scheme_results": results,
        "figure_paths": {
            "jump_lower_orders_compare": str(lower_compare_path),
            "scheme_C_eta_mismatch": str(scheme_c_mismatch_path),
            "scheme_B_consensus_improvement": str(scheme_b_consensus_path),
        },
        "table_paths": {
            "scheme_comparison_summary": str(comparison_csv_path),
        },
    }
    _write_interpretation_summary(compare_result, interpretation_path, s=s)
    compare_result["table_paths"]["interpretation_summary"] = str(interpretation_path)
    return compare_result


def _default_q_sampler(rng: np.random.Generator, M: int) -> np.ndarray:
    increments = rng.normal(loc=0.0, scale=0.9, size=M)
    q = np.zeros((M + 1,), dtype=float)
    q[1:] = np.cumsum(increments)
    return q


def _default_T_sampler(rng: np.random.Generator, M: int) -> np.ndarray:
    return rng.uniform(0.5, 2.0, size=M).astype(float)


def run_random_trials(
    n_trials: int = 100,
    M: int = 20,
    s: int = DEFAULT_S,
    k: int = DEFAULT_K,
    q_sampler: Callable[[np.random.Generator, int], np.ndarray] | None = None,
    T_sampler: Callable[[np.random.Generator, int], np.ndarray] | None = None,
    schemes: tuple[str, ...] = ("A", "B", "C"),
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    """Run randomized stability trials for the requested assembly schemes."""
    if n_trials < 1:
        raise ValueError("n_trials must be >= 1.")
    rng = np.random.default_rng(seed)
    q_sampler = _default_q_sampler if q_sampler is None else q_sampler
    T_sampler = _default_T_sampler if T_sampler is None else T_sampler

    records: list[dict[str, Any]] = []
    for trial_idx in range(n_trials):
        q = np.asarray(q_sampler(rng, M), dtype=float).reshape(M + 1)
        T = np.asarray(T_sampler(rng, M), dtype=float).reshape(M)
        for scheme in schemes:
            start_time = time.perf_counter()
            try:
                result = run_boundary_jump_check(q, T, s=s, k=k, scheme=scheme, config=None, save_dir=None, seed=seed)
                elapsed = time.perf_counter() - start_time
                lower = max((result["stats"][order]["max_abs"] for order in range(s)), default=0.0)
                higher = max((result["stats"][order]["max_abs"] for order in range(s, 2 * s - 1)), default=0.0)
                records.append(
                    {
                        "trial": trial_idx,
                        "scheme": scheme,
                        "success": True,
                        "lower_order_max": lower,
                        "higher_order_max": higher,
                        "elapsed_sec": elapsed,
                    }
                )
            except Exception as exc:  # pragma: no cover - defensive path for batch logging
                elapsed = time.perf_counter() - start_time
                records.append(
                    {
                        "trial": trial_idx,
                        "scheme": scheme,
                        "success": False,
                        "lower_order_max": math.nan,
                        "higher_order_max": math.nan,
                        "elapsed_sec": elapsed,
                        "error": str(exc),
                    }
                )

    summary: dict[str, dict[str, Any]] = {}
    for scheme in schemes:
        scheme_records = [record for record in records if record["scheme"] == scheme]
        successes = [record for record in scheme_records if record["success"]]
        lower_vals = np.asarray([record["lower_order_max"] for record in successes], dtype=float)
        higher_vals = np.asarray([record["higher_order_max"] for record in successes], dtype=float)
        time_vals = np.asarray([record["elapsed_sec"] for record in successes], dtype=float)
        summary[scheme] = {
            "success_rate": float(len(successes) / max(len(scheme_records), 1)),
            "mean_lower_order_max": float(np.mean(lower_vals)) if lower_vals.size else math.nan,
            "mean_higher_order_max": float(np.mean(higher_vals)) if higher_vals.size else math.nan,
            "mean_elapsed_sec": float(np.mean(time_vals)) if time_vals.size else math.nan,
        }

    if save_dir is not None:
        save_root = _ensure_dir(save_dir)
        assert save_root is not None
        boxplot_path = save_root / "jump_boxplot_random_trials.png"
        _plot_random_boxplot([record for record in records if record["success"]], s=s, save_path=boxplot_path)
        _save_csv(
            save_root / "random_trials_summary.csv",
            ["trial", "scheme", "success", "lower_order_max", "higher_order_max", "elapsed_sec", "error"],
            [
                [
                    record["trial"],
                    record["scheme"],
                    record["success"],
                    record.get("lower_order_max", ""),
                    record.get("higher_order_max", ""),
                    record.get("elapsed_sec", ""),
                    record.get("error", ""),
                ]
                for record in records
            ],
        )
        _save_json(save_root / "random_trials_aggregate.json", _serialize(summary))

    return {
        "records": records,
        "summary": summary,
    }


def main() -> None:
    q = np.asarray([0.0, 1.1, -0.6, 1.4, -0.2, 0.9, 0.1], dtype=float)
    T = np.asarray([0.8, 1.4, 0.9, 1.3, 0.7, 1.1], dtype=float)
    compare = run_compare_all_schemes(q, T, save_dir=DEFAULT_RESULTS_DIR)
    random_result = run_random_trials(n_trials=24, M=12, save_dir=DEFAULT_RESULTS_DIR, seed=42)

    print("Phase 5 boundary jump check")
    for scheme in ("A", "B", "C"):
        stats = compare["scheme_results"][scheme]["stats_overview"]
        print(
            f"scheme {scheme}: lower max={stats['lower_order_max']:.3e}, "
            f"higher max={stats['higher_order_max']:.3e}"
        )
    print(
        "random success rates: "
        + ", ".join(
            f"{scheme}={random_result['summary'][scheme]['success_rate']:.2f}" for scheme in ("A", "B", "C")
        )
    )
    print(f"results dir: {DEFAULT_RESULTS_DIR}")


if __name__ == "__main__":
    main()

