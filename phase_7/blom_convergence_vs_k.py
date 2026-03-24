"""
Phase 7 convergence validation for idealized and actual BLOM-k families.

This module compares three related objects against the exact global Phase 1
MINCO reference:

- the exact MINCO coefficient map
- the idealized truncated kernel surrogate
- the actual assembled BLOM-k trajectories from Schemes A/B/C
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

from phase_0.poly_basis import COEFFS_PER_SEGMENT, snap_cost_matrix
from phase_1.minco_scalar_baseline import solve_minco_coefficients
from phase_2.phase2_validation import build_waypoint_selector, compute_exact_jacobian_q
from phase_5.blom_boundary_jump_check import (
    assemble_scheme_A,
    assemble_scheme_B,
    assemble_scheme_C,
)


DEFAULT_S = 4
DEFAULT_RESULTS_DIR = Path("phase_7/results/phase7_convergence_vs_k")
DEFAULT_ZERO_TOL = 1e-12


def _validate_time_vector(T: np.ndarray) -> np.ndarray:
    T = np.asarray(T, dtype=float).reshape(-1)
    if T.size < 1:
        raise ValueError("T must describe at least one segment.")
    if not np.all(np.isfinite(T)):
        raise ValueError("T contains non-finite values.")
    if not np.all(T > 0.0):
        raise ValueError("All durations must be strictly positive.")
    return T


def _validate_problem_inputs(q: np.ndarray, T: np.ndarray, s: int) -> tuple[np.ndarray, np.ndarray]:
    if s != DEFAULT_S:
        raise ValueError("Phase 7 currently supports only the canonical case s=4.")
    q = np.asarray(q, dtype=float).reshape(-1)
    T = _validate_time_vector(T)
    if q.shape != (T.size + 1,):
        raise ValueError(f"q must have shape ({T.size + 1},), got {q.shape}.")
    if not np.all(np.isfinite(q)):
        raise ValueError("q contains non-finite values.")
    return q, T


def representative_case() -> tuple[np.ndarray, np.ndarray]:
    """Deterministic nonuniform case used by demos and RESULT_PHASE7.md."""
    q = np.asarray([0.0, 1.10, -0.70, 1.55, -0.40, 0.95, -0.20, 1.05, 0.10], dtype=float)
    T = np.asarray([0.80, 1.35, 0.95, 1.25, 0.75, 1.10, 0.90, 1.40], dtype=float)
    return q, T


def default_boundary_jets(q: np.ndarray, s: int = DEFAULT_S) -> tuple[np.ndarray, np.ndarray]:
    q = np.asarray(q, dtype=float).reshape(-1)
    zeta_start = np.zeros((s,), dtype=float)
    zeta_end = np.zeros((s,), dtype=float)
    zeta_start[0] = float(q[0])
    zeta_end[0] = float(q[-1])
    return zeta_start, zeta_end


def ensure_results_dirs(base_dir: str | Path = DEFAULT_RESULTS_DIR) -> dict[str, Path]:
    base_dir = Path(base_dir)
    ideal_dir = base_dir / "ideal"
    scheme_a = base_dir / "scheme_A"
    scheme_b = base_dir / "scheme_B"
    scheme_c = base_dir / "scheme_C"
    compare_dir = base_dir / "compare"
    for directory in (base_dir, ideal_dir, scheme_a, scheme_b, scheme_c, compare_dir):
        directory.mkdir(parents=True, exist_ok=True)
    return {
        "base": base_dir,
        "ideal": ideal_dir,
        "scheme_A": scheme_a,
        "scheme_B": scheme_b,
        "scheme_C": scheme_c,
        "compare": compare_dir,
    }


def _save_json(path: str | Path, payload: dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _save_csv(path: str | Path, header: list[str], rows: list[list[Any]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rendered = [",".join(header)]
    for row in rows:
        rendered.append(",".join(str(value) for value in row))
    path.write_text("\n".join(rendered) + "\n", encoding="utf-8")


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
    if isinstance(value, tuple):
        return [_serialize(item) for item in value]
    return value


def make_k_grid(M: int, start: int = 2, step: int = 2, include_M: bool = True) -> list[int]:
    """Construct the BLOM-k sweep grid."""
    if M < 1:
        raise ValueError("M must be >= 1.")
    if start < 1:
        raise ValueError("start must be >= 1.")
    if step < 1:
        raise ValueError("step must be >= 1.")
    values: list[int] = []
    current = start
    while current <= M:
        values.append(int(current))
        current += step
    if include_M and M not in values:
        values.append(int(M))
    if not values:
        values.append(int(M))
    return sorted(set(values))


def compute_total_snap_cost(coeffs: np.ndarray, T: np.ndarray) -> float:
    coeffs = np.asarray(coeffs, dtype=float)
    T = _validate_time_vector(T)
    if coeffs.shape != (T.size, COEFFS_PER_SEGMENT):
        raise ValueError(f"coeffs must have shape ({T.size}, {COEFFS_PER_SEGMENT}), got {coeffs.shape}.")
    total = 0.0
    for segment_index, duration in enumerate(T):
        total += float(coeffs[segment_index] @ snap_cost_matrix(float(duration)) @ coeffs[segment_index])
    return total


def _zero_centered_waypoints(q: np.ndarray) -> np.ndarray:
    q = np.asarray(q, dtype=float).reshape(-1)
    centered = q.copy()
    if centered.size > 2:
        centered[1:-1] = 0.0
    return centered


def _row_block(segment_index: int) -> slice:
    start = segment_index * COEFFS_PER_SEGMENT
    stop = start + COEFFS_PER_SEGMENT
    return slice(start, stop)


def compute_minco_reference(q: np.ndarray, T: np.ndarray, s: int = DEFAULT_S) -> dict[str, Any]:
    """Compute the exact Phase 1 global MINCO reference and its exact waypoint kernel."""
    q, T = _validate_problem_inputs(q, T, s=s)
    zeta_start, zeta_end = default_boundary_jets(q, s=s)

    solved = solve_minco_coefficients(q, T, zeta_start, zeta_end, return_system=True)
    q_zero = _zero_centered_waypoints(q)
    zero_solution = solve_minco_coefficients(q_zero, T, zeta_start, zeta_end, return_system=False)
    S_q = build_waypoint_selector(T.size, s=s)
    J_q = compute_exact_jacobian_q(solved["M"], S_q)

    c_reconstructed = zero_solution["c_vec"] + J_q @ q[1:-1]
    return {
        "q": q,
        "T": T,
        "s": s,
        "zeta_start": zeta_start,
        "zeta_end": zeta_end,
        "coeffs": solved["coeffs"],
        "c_vec": solved["c_vec"],
        "cost": compute_total_snap_cost(solved["coeffs"], T),
        "system": {"M": solved["M"], "b": solved["b"], "residual_norm": float(solved["residual_norm"])},
        "S_q": S_q,
        "J_q": J_q,
        "q_interior": q[1:-1].copy(),
        "c_zero": zero_solution["c_vec"],
        "centered_reconstruction_error": float(np.linalg.norm(c_reconstructed - solved["c_vec"])),
    }


def _ideal_mask(M: int, k: int) -> np.ndarray:
    num_interior = max(M - 1, 0)
    mask = np.zeros((COEFFS_PER_SEGMENT * M, num_interior), dtype=bool)
    for segment_index in range(1, M + 1):
        allowed = {waypoint_index for waypoint_index in range(1, M) if abs(waypoint_index - segment_index) <= k}
        row_slice = _row_block(segment_index - 1)
        for waypoint_index in allowed:
            mask[row_slice, waypoint_index - 1] = True
    return mask


def compute_ideal_truncated_blom_k(
    q: np.ndarray,
    T: np.ndarray,
    k: int,
    s: int = DEFAULT_S,
    *,
    reference: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Compute the idealized truncated-kernel surrogate defined in Phase 7."""
    q, T = _validate_problem_inputs(q, T, s=s)
    if k < 1:
        raise ValueError("k must be >= 1.")
    reference = compute_minco_reference(q, T, s=s) if reference is None else reference
    mask = _ideal_mask(T.size, k)
    J_trunc = np.where(mask, reference["J_q"], 0.0)
    c_vec = reference["c_zero"] + J_trunc @ reference["q_interior"]
    coeffs = c_vec.reshape(T.size, COEFFS_PER_SEGMENT)
    centered_tail = reference["c_vec"] - c_vec
    return {
        "kind": "ideal",
        "k": int(k),
        "mask": mask,
        "J_trunc": J_trunc,
        "coeffs": coeffs,
        "c_vec": c_vec,
        "cost": compute_total_snap_cost(coeffs, T),
        "tail_norm": float(np.linalg.norm(centered_tail)),
        "tail_max_abs": float(np.max(np.abs(centered_tail))) if centered_tail.size else 0.0,
    }


def compute_actual_blom_k(
    q: np.ndarray,
    T: np.ndarray,
    k: int,
    s: int = DEFAULT_S,
    scheme: str = "C",
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Compute the actual BLOM-k object from one of the assembly schemes."""
    q, T = _validate_problem_inputs(q, T, s=s)
    if k < 1:
        raise ValueError("k must be >= 1.")
    scheme = str(scheme).upper()
    assemblers: dict[str, Callable[..., dict[str, Any]]] = {
        "A": assemble_scheme_A,
        "B": assemble_scheme_B,
        "C": assemble_scheme_C,
    }
    if scheme not in assemblers:
        raise ValueError(f"Unsupported scheme {scheme!r}.")

    result = assemblers[scheme](q, T, s=s, k=k, config=config)
    coeffs = np.asarray(result["coeffs"], dtype=float)
    meta = dict(result.get("meta", {}))
    meta["requested_k"] = int(k)
    meta["is_k_dependent"] = bool(scheme in {"B", "C"})
    meta["comparison_role"] = "shared_state_baseline" if scheme == "A" else "actual_blom_k"
    return {
        **result,
        "scheme": scheme,
        "coeffs": coeffs,
        "c_vec": coeffs.reshape(-1),
        "cost": compute_total_snap_cost(coeffs, T),
        "meta": meta,
    }


def compute_convergence_errors(
    c_ref: np.ndarray,
    c_test: np.ndarray,
    c_ideal: np.ndarray | None = None,
    cost_ref: float | None = None,
    cost_test: float | None = None,
) -> dict[str, Any]:
    """Compute Phase 7 coefficient, matching, and cost-gap metrics."""
    c_ref = np.asarray(c_ref, dtype=float).reshape(-1)
    c_test = np.asarray(c_test, dtype=float).reshape(-1)
    if c_ref.shape != c_test.shape:
        raise ValueError(f"c_ref and c_test must match, got {c_ref.shape} and {c_test.shape}.")
    diff = c_test - c_ref
    ref_norm = float(np.linalg.norm(c_ref))
    global_l2 = float(np.linalg.norm(diff))
    segmentwise = np.linalg.norm(diff.reshape(-1, COEFFS_PER_SEGMENT), axis=1)

    metrics = {
        "global_l2": global_l2,
        "relative_l2": global_l2 / max(ref_norm, 1e-15),
        "max_abs": float(np.max(np.abs(diff))) if diff.size else 0.0,
        "segmentwise_l2": segmentwise.tolist(),
        "max_segmentwise_l2": float(np.max(segmentwise)) if segmentwise.size else 0.0,
        "mean_segmentwise_l2": float(np.mean(segmentwise)) if segmentwise.size else 0.0,
    }
    if c_ideal is not None:
        c_ideal = np.asarray(c_ideal, dtype=float).reshape(-1)
        if c_ideal.shape != c_ref.shape:
            raise ValueError(f"c_ideal must have shape {c_ref.shape}, got {c_ideal.shape}.")
        match = c_test - c_ideal
        metrics["matching_l2"] = float(np.linalg.norm(match))
        metrics["matching_max_abs"] = float(np.max(np.abs(match))) if match.size else 0.0
    if cost_ref is not None and cost_test is not None:
        metrics["abs_cost_gap"] = float(cost_test - cost_ref)
        metrics["rel_cost_gap"] = float(cost_test / max(cost_ref, 1e-15) - 1.0)
    return metrics


def fit_log_error_vs_k(k_values: list[int] | np.ndarray, errors: list[float] | np.ndarray, eps: float = 1e-15) -> dict[str, Any]:
    """Fit log(error) ~ alpha + beta k by least squares."""
    x = np.asarray(k_values, dtype=float).reshape(-1)
    y_raw = np.asarray(errors, dtype=float).reshape(-1)
    if x.shape != y_raw.shape:
        raise ValueError(f"k_values and errors must have the same shape, got {x.shape} and {y_raw.shape}.")
    valid = np.isfinite(x) & np.isfinite(y_raw)
    if not np.any(valid):
        return {
            "slope": math.nan,
            "intercept": math.nan,
            "r2": math.nan,
            "k_min": math.nan,
            "k_max": math.nan,
            "valid_count": 0,
            "log_errors": [],
            "fit_log_errors": [],
        }
    x_valid = x[valid]
    y_log = np.log(np.maximum(y_raw[valid], eps))
    design = np.column_stack([np.ones_like(x_valid), x_valid])
    coeffs, *_ = np.linalg.lstsq(design, y_log, rcond=None)
    intercept = float(coeffs[0])
    slope = float(coeffs[1])
    fitted = design @ coeffs
    ss_res = float(np.sum((y_log - fitted) ** 2))
    ss_tot = float(np.sum((y_log - np.mean(y_log)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else 1.0
    return {
        "slope": slope,
        "intercept": intercept,
        "r2": r2,
        "k_min": int(np.min(x_valid)),
        "k_max": int(np.max(x_valid)),
        "valid_count": int(x_valid.size),
        "log_errors": y_log.tolist(),
        "fit_log_errors": fitted.tolist(),
    }


def _plot_error_curves(
    k_values: list[int],
    curves: list[tuple[str, list[float]]],
    save_path: str | Path,
    *,
    title: str,
    ylabel: str,
    logy: bool = True,
) -> None:
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    for label, values in curves:
        ax.plot(k_values, values, marker="o", linewidth=1.9, label=label)
    ax.set_title(title)
    ax.set_xlabel("k")
    ax.set_ylabel(ylabel)
    if logy:
        ax.set_yscale("log")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_logfit(
    k_values: list[int],
    curves: list[tuple[str, list[float], dict[str, Any]]],
    save_path: str | Path,
) -> None:
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    x = np.asarray(k_values, dtype=float)
    for label, errors, fit in curves:
        y = np.log(np.maximum(np.asarray(errors, dtype=float), 1e-15))
        ax.scatter(x, y, s=28, alpha=0.75, label=f"{label} data")
        if fit["valid_count"] > 0:
            fitted = fit["intercept"] + fit["slope"] * x
            ax.plot(x, fitted, linewidth=1.8, linestyle="--", label=f"{label} fit")
    ax.set_title("log Error vs k Linear Fits")
    ax.set_xlabel("k")
    ax.set_ylabel("log(error)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_segmentwise_heatmap(
    mat: np.ndarray,
    k_values: list[int],
    save_path: str | Path,
    *,
    title: str,
) -> None:
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    display = np.log10(np.asarray(mat, dtype=float) + 1e-16)
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    image = ax.imshow(display, aspect="auto", origin="lower", cmap="magma")
    ax.set_title(title)
    ax.set_xlabel("Segment Index")
    ax.set_ylabel("k sweep index")
    ax.set_xticks(np.arange(mat.shape[1]))
    ax.set_xticklabels([str(index + 1) for index in range(mat.shape[1])])
    ax.set_yticks(np.arange(len(k_values)))
    ax.set_yticklabels([str(k) for k in k_values])
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("log10(segmentwise error + 1e-16)")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_slope_boxplot(records: list[dict[str, Any]], save_path: str | Path) -> None:
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    curves = sorted({record["curve"] for record in records if record.get("success", False)})
    values = [
        [record["slope"] for record in records if record.get("success", False) and record["curve"] == curve]
        for curve in curves
    ]
    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    ax.boxplot(values, tick_labels=curves)
    ax.set_title("Random-Trial log-fit slope distribution")
    ax.set_xlabel("Curve")
    ax.set_ylabel("Fitted slope beta")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _interpret_curve(fit: dict[str, Any], label: str) -> str:
    slope = fit["slope"]
    r2 = fit["r2"]
    if not np.isfinite(slope):
        return f"- {label}: fit unavailable."
    if abs(slope) < 1e-3:
        return f"- {label}: essentially flat across the tested k range (slope={slope:.3e}, R^2={r2:.3f})."
    if slope < -0.05 and r2 >= 0.8:
        return f"- {label}: numerically consistent with approximate exponential decay (slope={slope:.3e}, R^2={r2:.3f})."
    if slope < 0.0:
        return f"- {label}: decreases with k, but the log-linear fit is only moderate (slope={slope:.3e}, R^2={r2:.3f})."
    return f"- {label}: no convincing decay trend in this run (slope={slope:.3e}, R^2={r2:.3f})."


def _write_interpretation_summary(result: dict[str, Any], save_path: str | Path) -> None:
    fits = result["fits"]
    final_by_curve = result["final_by_curve"]
    best_scheme = min(
        (scheme for scheme in result["schemes"]),
        key=lambda scheme: final_by_curve[scheme]["global_l2"],
    )
    strict_k_schemes = [scheme for scheme in result["schemes"] if scheme in {"B", "C"}]
    best_strict_scheme = (
        min(strict_k_schemes, key=lambda scheme: final_by_curve[scheme]["global_l2"]) if strict_k_schemes else None
    )
    scheme_a_note = ""
    if "A" in result["schemes"]:
        scheme_a_note = (
            "\n- Scheme A is currently treated as a shared-state baseline. "
            "In the present Phase 5 implementation its assembly ignores k, so its Phase 7 trend should be read as a baseline, not as a strict BLOM-k family."
        )

    lines = [
        "# Phase 7 Interpretation Summary",
        "",
        "## Exponential-Decay Reading",
        _interpret_curve(fits["ideal"], "Ideal truncation"),
    ]
    for scheme in result["schemes"]:
        lines.append(_interpret_curve(fits[scheme], f"Scheme {scheme}"))
    lines.extend(
        [
            "",
            "## Best MINCO Approximation in This Run",
            f"- Best final actual scheme by coefficient error: `Scheme {best_scheme}`",
            (
                f"- Best final strict k-family by coefficient error: `Scheme {best_strict_scheme}`"
                if best_strict_scheme is not None
                else "- No strict k-family scheme was included in this run."
            ),
            f"- Final ideal error: `{final_by_curve['ideal']['global_l2']:.3e}`",
        ]
    )
    for scheme in result["schemes"]:
        final_metrics = final_by_curve[scheme]
        lines.append(
            f"- Scheme {scheme}: final error `{final_metrics['global_l2']:.3e}`, "
            f"matching `{final_metrics.get('matching_l2', math.nan):.3e}`, "
            f"relative cost gap `{final_metrics.get('rel_cost_gap', math.nan):.3e}`"
        )
    lines.extend(
        [
            "",
            "## Matching Assessment",
        ]
    )
    for scheme in result["schemes"]:
        fit = result["matching_fits"][scheme]
        lines.append(_interpret_curve(fit, f"Scheme {scheme} vs ideal"))
    lines.extend(
        [
            "",
            "## Theory Recommendation",
            "- If the ideal truncation decays cleanly while actual Scheme B/C track it closely, the next theory target should remain the missing matching theorem.",
            "- If actual errors stagnate while ideal truncation still decays, the priority should shift from theorem completion toward revising the assembly or the comparison target.",
            "- Cost-gap trends should be read as diagnostic evidence only; they do not by themselves prove near-optimality because the actual assembled candidates are not the full Phase 1 affine-feasible family.",
            scheme_a_note,
        ]
    )
    Path(save_path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_convergence_vs_k(
    q: np.ndarray,
    T: np.ndarray,
    s: int = DEFAULT_S,
    k_values: list[int] | None = None,
    schemes: tuple[str, ...] = ("A", "B", "C"),
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    """Run one Phase 7 convergence sweep over k."""
    del seed
    q, T = _validate_problem_inputs(q, T, s=s)
    schemes = tuple(str(scheme).upper() for scheme in schemes)
    if k_values is None:
        k_values = make_k_grid(T.size)
    k_values = [int(k) for k in k_values]
    start_time = time.perf_counter()

    reference = compute_minco_reference(q, T, s=s)
    ideal_records: list[dict[str, Any]] = []
    actual_records: dict[str, list[dict[str, Any]]] = {scheme: [] for scheme in schemes}
    segmentwise_errors: dict[str, list[list[float]]] = {"ideal": []}
    for scheme in schemes:
        segmentwise_errors[scheme] = []

    for k in k_values:
        ideal = compute_ideal_truncated_blom_k(q, T, k, s=s, reference=reference)
        ideal_metrics = compute_convergence_errors(
            reference["c_vec"],
            ideal["c_vec"],
            cost_ref=reference["cost"],
            cost_test=ideal["cost"],
        )
        ideal_record = {"k": int(k), "kind": "ideal", **ideal_metrics}
        ideal_records.append(ideal_record)
        segmentwise_errors["ideal"].append(ideal_metrics["segmentwise_l2"])

        for scheme in schemes:
            actual = compute_actual_blom_k(q, T, k, s=s, scheme=scheme, config=None)
            actual_metrics = compute_convergence_errors(
                reference["c_vec"],
                actual["c_vec"],
                c_ideal=ideal["c_vec"],
                cost_ref=reference["cost"],
                cost_test=actual["cost"],
            )
            record = {
                "k": int(k),
                "kind": "actual",
                "scheme": scheme,
                "is_k_dependent": bool(actual["meta"]["is_k_dependent"]),
                **actual_metrics,
            }
            actual_records[scheme].append(record)
            segmentwise_errors[scheme].append(actual_metrics["segmentwise_l2"])

    fits = {
        "ideal": fit_log_error_vs_k(k_values, [record["global_l2"] for record in ideal_records]),
    }
    matching_fits: dict[str, dict[str, Any]] = {}
    for scheme in schemes:
        fits[scheme] = fit_log_error_vs_k(k_values, [record["global_l2"] for record in actual_records[scheme]])
        matching_fits[scheme] = fit_log_error_vs_k(
            k_values,
            [record["matching_l2"] for record in actual_records[scheme]],
        )

    final_by_curve = {"ideal": ideal_records[-1]}
    for scheme in schemes:
        final_by_curve[scheme] = actual_records[scheme][-1]

    elapsed = time.perf_counter() - start_time
    result = {
        "inputs": {"q": q, "T": T, "s": s, "k_values": k_values},
        "reference": {
            "cost": reference["cost"],
            "residual_norm": reference["system"]["residual_norm"],
            "centered_reconstruction_error": reference["centered_reconstruction_error"],
        },
        "ideal_records": ideal_records,
        "actual_records": actual_records,
        "fits": fits,
        "matching_fits": matching_fits,
        "schemes": list(schemes),
        "final_by_curve": final_by_curve,
        "segmentwise_errors": segmentwise_errors,
        "elapsed_sec": elapsed,
    }

    if save_dir is None:
        return result

    dirs = ensure_results_dirs(save_dir)
    compare_dir = dirs["compare"]
    rows: list[list[Any]] = []
    for record in ideal_records:
        rows.append(
            [
                "ideal",
                "",
                record["k"],
                record["global_l2"],
                record["relative_l2"],
                "",
                record["abs_cost_gap"],
                record["rel_cost_gap"],
                record["max_segmentwise_l2"],
            ]
        )
    for scheme in schemes:
        for record in actual_records[scheme]:
            rows.append(
                [
                    "actual",
                    scheme,
                    record["k"],
                    record["global_l2"],
                    record["relative_l2"],
                    record["matching_l2"],
                    record["abs_cost_gap"],
                    record["rel_cost_gap"],
                    record["max_segmentwise_l2"],
                ]
            )

    _save_csv(
        compare_dir / "convergence_errors_by_k.csv",
        [
            "kind",
            "scheme",
            "k",
            "global_l2",
            "relative_l2",
            "matching_l2",
            "abs_cost_gap",
            "rel_cost_gap",
            "max_segmentwise_l2",
        ],
        rows,
    )

    fit_rows = [
        ["ideal", "", fits["ideal"]["slope"], fits["ideal"]["intercept"], fits["ideal"]["r2"], fits["ideal"]["k_min"], fits["ideal"]["k_max"], True]
    ]
    for scheme in schemes:
        fit_rows.append(
            [
                "actual",
                scheme,
                fits[scheme]["slope"],
                fits[scheme]["intercept"],
                fits[scheme]["r2"],
                fits[scheme]["k_min"],
                fits[scheme]["k_max"],
                bool(scheme != "A"),
            ]
        )
    _save_csv(
        compare_dir / "logfit_summary.csv",
        ["kind", "scheme", "slope", "intercept", "r2", "k_min", "k_max", "is_k_dependent"],
        fit_rows,
    )

    _plot_error_curves(
        k_values,
        [
            ("ideal", [record["global_l2"] for record in ideal_records]),
            *[(f"Scheme {scheme}", [record["global_l2"] for record in actual_records[scheme]]) for scheme in schemes],
        ],
        compare_dir / "coef_error_vs_k_all.png",
        title="Coefficient Error vs k",
        ylabel="||c^(k) - c*||_2",
        logy=True,
    )

    _plot_logfit(
        k_values,
        [
            ("ideal", [record["global_l2"] for record in ideal_records], fits["ideal"]),
            *[
                (f"Scheme {scheme}", [record["global_l2"] for record in actual_records[scheme]], fits[scheme])
                for scheme in schemes
            ],
        ],
        compare_dir / "log_coef_error_fit_all.png",
    )

    _plot_error_curves(
        k_values,
        [(f"Scheme {scheme}", [record["matching_l2"] for record in actual_records[scheme]]) for scheme in schemes],
        compare_dir / "matching_error_vs_k.png",
        title="Matching Error vs k",
        ylabel="||c_actual^(k) - c_ideal^(k)||_2",
        logy=True,
    )

    _plot_error_curves(
        k_values,
        [
            ("ideal", [record["rel_cost_gap"] for record in ideal_records]),
            *[(f"Scheme {scheme}", [record["rel_cost_gap"] for record in actual_records[scheme]]) for scheme in schemes],
        ],
        compare_dir / "relative_cost_gap_vs_k.png",
        title="Relative Cost Gap vs k",
        ylabel="J(c^(k), T) / J(c*, T) - 1",
        logy=False,
    )

    _plot_error_curves(
        k_values,
        [
            ("ideal", [record["global_l2"] for record in ideal_records]),
            *[(f"Scheme {scheme}", [record["global_l2"] for record in actual_records[scheme]]) for scheme in schemes],
        ],
        compare_dir / "ideal_vs_actual_compare.png",
        title="Ideal vs Actual BLOM-k",
        ylabel="Coefficient error to MINCO",
        logy=True,
    )

    _plot_segmentwise_heatmap(
        np.asarray(segmentwise_errors["ideal"], dtype=float),
        k_values,
        dirs["ideal"] / "segmentwise_error_heatmap_ideal.png",
        title="Ideal truncation segmentwise error",
    )
    for scheme in schemes:
        _plot_segmentwise_heatmap(
            np.asarray(segmentwise_errors[scheme], dtype=float),
            k_values,
            dirs[f"scheme_{scheme}"] / f"segmentwise_error_heatmap_scheme_{scheme}.png",
            title=f"Scheme {scheme} segmentwise error",
        )
        _save_json(
            dirs[f"scheme_{scheme}"] / f"summary_scheme_{scheme}.json",
            _serialize(
                {
                    "scheme": scheme,
                    "records": actual_records[scheme],
                    "fit": fits[scheme],
                    "matching_fit": matching_fits[scheme],
                }
            ),
        )
    _save_json(
        dirs["ideal"] / "summary_ideal.json",
        _serialize({"records": ideal_records, "fit": fits["ideal"]}),
    )

    summary_json = {
        "reference": result["reference"],
        "fits": fits,
        "matching_fits": matching_fits,
        "final_by_curve": final_by_curve,
        "elapsed_sec": elapsed,
        "notes": {
            "scheme_A_role": "shared_state_baseline_not_strict_k_family" if "A" in schemes else "not_run",
        },
    }
    _save_json(compare_dir / "summary_phase7_convergence.json", _serialize(summary_json))
    _write_interpretation_summary(result, compare_dir / "phase7_interpretation_summary.md")
    result["paths"] = {
        "compare": compare_dir,
        "summary_json": compare_dir / "summary_phase7_convergence.json",
        "interpretation": compare_dir / "phase7_interpretation_summary.md",
    }
    return result


def _default_q_sampler(rng: np.random.Generator, M: int) -> np.ndarray:
    increments = rng.normal(loc=0.0, scale=0.85, size=M)
    q = np.zeros((M + 1,), dtype=float)
    q[1:] = np.cumsum(increments)
    return q


def _default_T_sampler(rng: np.random.Generator, M: int) -> np.ndarray:
    return rng.uniform(0.5, 2.0, size=M).astype(float)


def run_random_trials(
    n_trials: int = 100,
    M: int = 20,
    s: int = DEFAULT_S,
    k_values: list[int] | None = None,
    schemes: tuple[str, ...] = ("A", "B", "C"),
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    """Run randomized convergence sweeps and summarize slope/matching/cost statistics."""
    if n_trials < 1:
        raise ValueError("n_trials must be >= 1.")
    if M < 3:
        raise ValueError("M must be >= 3 for a meaningful Phase 7 sweep.")
    if k_values is None:
        k_values = make_k_grid(M)
    rng = np.random.default_rng(seed)
    records: list[dict[str, Any]] = []

    for trial_idx in range(n_trials):
        q = _default_q_sampler(rng, M)
        T = _default_T_sampler(rng, M)
        start = time.perf_counter()
        try:
            result = run_convergence_vs_k(q, T, s=s, k_values=k_values, schemes=schemes, save_dir=None, seed=seed)
            elapsed = time.perf_counter() - start
            records.append(
                {
                    "trial": trial_idx,
                    "curve": "ideal",
                    "scheme": "",
                    "success": True,
                    "slope": result["fits"]["ideal"]["slope"],
                    "r2": result["fits"]["ideal"]["r2"],
                    "final_global_l2": result["final_by_curve"]["ideal"]["global_l2"],
                    "final_rel_cost_gap": result["final_by_curve"]["ideal"]["rel_cost_gap"],
                    "final_matching_l2": math.nan,
                    "elapsed_sec": elapsed,
                }
            )
            for scheme in schemes:
                records.append(
                    {
                        "trial": trial_idx,
                        "curve": scheme,
                        "scheme": scheme,
                        "success": True,
                        "slope": result["fits"][scheme]["slope"],
                        "r2": result["fits"][scheme]["r2"],
                        "final_global_l2": result["final_by_curve"][scheme]["global_l2"],
                        "final_rel_cost_gap": result["final_by_curve"][scheme]["rel_cost_gap"],
                        "final_matching_l2": result["final_by_curve"][scheme]["matching_l2"],
                        "elapsed_sec": elapsed,
                    }
                )
        except Exception as exc:  # pragma: no cover - defensive batch logging
            elapsed = time.perf_counter() - start
            records.append(
                {
                    "trial": trial_idx,
                    "curve": "ideal",
                    "scheme": "",
                    "success": False,
                    "slope": math.nan,
                    "r2": math.nan,
                    "final_global_l2": math.nan,
                    "final_rel_cost_gap": math.nan,
                    "final_matching_l2": math.nan,
                    "elapsed_sec": elapsed,
                    "error": str(exc),
                }
            )
            for scheme in schemes:
                records.append(
                    {
                        "trial": trial_idx,
                        "curve": scheme,
                        "scheme": scheme,
                        "success": False,
                        "slope": math.nan,
                        "r2": math.nan,
                        "final_global_l2": math.nan,
                        "final_rel_cost_gap": math.nan,
                        "final_matching_l2": math.nan,
                        "elapsed_sec": elapsed,
                        "error": str(exc),
                    }
                )

    curves = ["ideal", *schemes]
    summary: dict[str, Any] = {"curves": {}}
    for curve in curves:
        curve_records = [record for record in records if record["curve"] == curve]
        successes = [record for record in curve_records if record["success"]]
        slopes = np.asarray([record["slope"] for record in successes], dtype=float)
        r2_vals = np.asarray([record["r2"] for record in successes], dtype=float)
        match_vals = np.asarray([record["final_matching_l2"] for record in successes], dtype=float)
        cost_vals = np.asarray([record["final_rel_cost_gap"] for record in successes], dtype=float)
        summary["curves"][curve] = {
            "success_rate": float(len(successes) / max(len(curve_records), 1)),
            "mean_slope": float(np.nanmean(slopes)) if slopes.size else math.nan,
            "mean_r2": float(np.nanmean(r2_vals)) if r2_vals.size else math.nan,
            "median_final_matching_l2": float(np.nanmedian(match_vals)) if np.any(np.isfinite(match_vals)) else math.nan,
            "median_final_rel_cost_gap": float(np.nanmedian(cost_vals)) if np.any(np.isfinite(cost_vals)) else math.nan,
        }

    if save_dir is not None:
        dirs = ensure_results_dirs(save_dir)
        compare_dir = dirs["compare"]
        _save_csv(
            compare_dir / "random_trials_summary.csv",
            [
                "trial",
                "curve",
                "scheme",
                "success",
                "slope",
                "r2",
                "final_global_l2",
                "final_rel_cost_gap",
                "final_matching_l2",
                "elapsed_sec",
            ],
            [
                [
                    record["trial"],
                    record["curve"],
                    record["scheme"],
                    record["success"],
                    record["slope"],
                    record["r2"],
                    record["final_global_l2"],
                    record["final_rel_cost_gap"],
                    record["final_matching_l2"],
                    record["elapsed_sec"],
                ]
                for record in records
            ],
        )
        _plot_slope_boxplot([record for record in records if record["success"]], compare_dir / "logfit_slope_boxplot_random_trials.png")
        _save_json(compare_dir / "random_trials_aggregate.json", _serialize(summary))

    return {"records": records, "summary": summary}


def main() -> None:
    q, T = representative_case()
    convergence = run_convergence_vs_k(q, T, save_dir=DEFAULT_RESULTS_DIR)
    random_result = run_random_trials(n_trials=12, M=12, save_dir=DEFAULT_RESULTS_DIR, seed=42)

    print("Phase 7 convergence validation")
    print(f"k grid: {convergence['inputs']['k_values']}")
    print(f"reference cost: {convergence['reference']['cost']:.6e}")
    print(f"ideal final error: {convergence['final_by_curve']['ideal']['global_l2']:.6e}")
    for scheme in convergence["schemes"]:
        final_metrics = convergence["final_by_curve"][scheme]
        print(
            f"scheme {scheme}: final error={final_metrics['global_l2']:.6e}, "
            f"matching={final_metrics['matching_l2']:.6e}, "
            f"rel cost gap={final_metrics['rel_cost_gap']:.6e}"
        )
    print(
        "random success rates: "
        + ", ".join(
            f"{curve}={random_result['summary']['curves'][curve]['success_rate']:.2f}"
            for curve in ["ideal", *convergence["schemes"]]
        )
    )
    print(f"results dir: {DEFAULT_RESULTS_DIR}")


if __name__ == "__main__":
    main()
