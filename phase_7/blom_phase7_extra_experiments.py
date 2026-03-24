"""
Phase 7 extra experiments for post-convergence validation.

This module adds four follow-up experiment families on top of the main
Phase 7 convergence toolbox:

- larger-M sweep
- uniform-time vs bounded-nonuniform time-regime split
- interior-only vs full-error separation
- light corrected variants of Scheme C
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Callable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from phase_2.phase2_validation import estimate_effective_bandwidth
from phase_4.utils.hermite_utils import normalized_derivative_row
from phase_5.blom_boundary_jump_check import compute_jumps, summarize_jumps
from phase_6.blom_fd_jacobian_check import finite_difference_jacobian
from phase_7.blom_convergence_vs_k import (
    COEFFS_PER_SEGMENT,
    DEFAULT_S,
    compute_actual_blom_k,
    compute_convergence_errors,
    compute_ideal_truncated_blom_k,
    compute_minco_reference,
    compute_total_snap_cost,
    fit_log_error_vs_k,
    make_k_grid,
    representative_case,
    run_convergence_vs_k,
)


DEFAULT_RESULTS_DIR = Path("phase_7/results/phase7_extra_experiments")
DEFAULT_BOUNDARY_TRIM = 2


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
        raise ValueError("Phase 7 extra experiments currently support only s=4.")
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
        "exp1": base_dir / "exp1_large_M_sweep",
        "exp2": base_dir / "exp2_time_regime_split",
        "exp3": base_dir / "exp3_interior_vs_full",
        "exp4": base_dir / "exp4_schemeC_light_assembly",
        "compare": base_dir / "compare_summary",
    }
    for directory in directories.values():
        directory.mkdir(parents=True, exist_ok=True)
    return directories


def _save_json(path: str | Path, payload: dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _save_csv(path: str | Path, header: list[str], rows: list[list[Any]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [",".join(header)]
    for row in rows:
        lines.append(",".join(str(value) for value in row))
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
    if isinstance(value, tuple):
        return [_serialize(item) for item in value]
    return value


def _default_q_sampler(rng: np.random.Generator, M: int) -> np.ndarray:
    increments = rng.normal(loc=0.0, scale=0.85, size=M)
    q = np.zeros((M + 1,), dtype=float)
    q[1:] = np.cumsum(increments)
    return q


def _default_T_sampler(rng: np.random.Generator, M: int, low: float = 0.5, high: float = 2.0) -> np.ndarray:
    return rng.uniform(low, high, size=M).astype(float)


def _k_values_for_M(M: int, k_max: int = 20) -> list[int]:
    max_k = min(M, k_max)
    return make_k_grid(max_k, start=2, step=2, include_M=True)


def _curve_label_order() -> list[str]:
    return ["ideal", "A", "B", "C", "C1", "C2"]


def _segment_slice(M: int, boundary_trim: int) -> slice:
    start = max(boundary_trim, 0)
    stop = max(M - max(boundary_trim, 0), start)
    if start >= stop:
        return slice(0, M)
    return slice(start, stop)


def compute_full_and_interior_errors(
    c_ref: np.ndarray,
    c_test: np.ndarray,
    *,
    c_ideal: np.ndarray | None = None,
    cost_ref: float | None = None,
    cost_test: float | None = None,
    boundary_trim: int = DEFAULT_BOUNDARY_TRIM,
) -> dict[str, Any]:
    """Extend Phase 7 metrics with interior-only statistics."""
    base = compute_convergence_errors(c_ref, c_test, c_ideal=c_ideal, cost_ref=cost_ref, cost_test=cost_test)
    ref_blocks = np.asarray(c_ref, dtype=float).reshape(-1, COEFFS_PER_SEGMENT)
    test_blocks = np.asarray(c_test, dtype=float).reshape(-1, COEFFS_PER_SEGMENT)
    diff_blocks = test_blocks - ref_blocks
    M = diff_blocks.shape[0]
    interior = _segment_slice(M, boundary_trim)
    interior_diff = diff_blocks[interior]
    boundary_mask = np.ones((M,), dtype=bool)
    boundary_mask[interior] = False
    boundary_diff = diff_blocks[boundary_mask]

    base["full_coef_error"] = base["global_l2"]
    base["interior_coef_error"] = float(np.linalg.norm(interior_diff.reshape(-1))) if interior_diff.size else 0.0
    base["boundary_coef_error"] = float(np.linalg.norm(boundary_diff.reshape(-1))) if boundary_diff.size else 0.0
    base["interior_segmentwise_l2"] = np.linalg.norm(interior_diff, axis=1).tolist() if interior_diff.size else []
    base["boundary_segmentwise_l2"] = np.linalg.norm(boundary_diff, axis=1).tolist() if boundary_diff.size else []

    if c_ideal is not None:
        ideal_blocks = np.asarray(c_ideal, dtype=float).reshape(-1, COEFFS_PER_SEGMENT)
        match_blocks = test_blocks - ideal_blocks
        base["full_matching_error"] = base["matching_l2"]
        base["interior_matching_error"] = float(np.linalg.norm(match_blocks[interior].reshape(-1)))
        base["boundary_matching_error"] = float(np.linalg.norm(match_blocks[boundary_mask].reshape(-1)))

    return base


def _reconstruct_segment_coeff(
    q_left: float,
    q_right: float,
    duration: float,
    eta_left: np.ndarray,
    eta_right: np.ndarray,
    s: int = DEFAULT_S,
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


def _poly_derivative_value(coeff: np.ndarray, t: float, order: int) -> float:
    coeff = np.asarray(coeff, dtype=float).reshape(-1)
    degree = coeff.size - 1
    row = normalized_derivative_row(float(t), order, degree)
    return float(coeff @ row)


def _extract_eta_from_coeff(coeff: np.ndarray, duration: float, s: int, side: str) -> np.ndarray:
    t = 0.0 if side == "left" else float(duration)
    return np.asarray([_poly_derivative_value(coeff, t, order) for order in range(1, s)], dtype=float)


def _default_boundary_jets(q: np.ndarray, s: int = DEFAULT_S) -> tuple[np.ndarray, np.ndarray]:
    q = np.asarray(q, dtype=float).reshape(-1)
    zeta_start = np.zeros((s,), dtype=float)
    zeta_end = np.zeros((s,), dtype=float)
    zeta_start[0] = float(q[0])
    zeta_end[0] = float(q[-1])
    return zeta_start, zeta_end


def assemble_scheme_C_light(
    q: np.ndarray,
    T: np.ndarray,
    k: int,
    *,
    mode: str = "C1",
    s: int = DEFAULT_S,
) -> dict[str, Any]:
    """Build corrected Scheme C coefficients from raw left/right endpoint jets."""
    q, T = _validate_problem_inputs(q, T, s=s)
    mode = str(mode).upper()
    raw = compute_actual_blom_k(q, T, k, s=s, scheme="C", config=None)
    eta_left = np.asarray(raw["eta_left"], dtype=float)
    eta_right = np.asarray(raw["eta_right"], dtype=float)
    zeta_start, zeta_end = _default_boundary_jets(q, s=s)

    if eta_left.size == 0:
        coeffs = np.asarray(raw["coeffs"], dtype=float)
        jumps = compute_jumps(coeffs, T, s)
        return {
            "scheme": mode,
            "coeffs": coeffs,
            "eta_corrected": np.zeros((0, s - 1), dtype=float),
            "jumps": jumps,
            "stats": summarize_jumps(jumps),
            "cost": compute_total_snap_cost(coeffs, T),
            "meta": {"requested_k": int(k), "is_k_dependent": True, "correction_mode": mode},
        }

    if mode == "C1":
        eta_corrected = 0.5 * (eta_left + eta_right)
        weights_left = np.ones((T.size - 1,), dtype=float)
        weights_right = np.ones((T.size - 1,), dtype=float)
    elif mode == "C2":
        weights_left = 1.0 / np.asarray(T[:-1], dtype=float)
        weights_right = 1.0 / np.asarray(T[1:], dtype=float)
        denom = (weights_left + weights_right)[:, None]
        eta_corrected = (weights_left[:, None] * eta_left + weights_right[:, None] * eta_right) / denom
    else:
        raise ValueError(f"Unsupported corrected Scheme C mode {mode!r}.")

    coeffs = np.zeros((T.size, 2 * s), dtype=float)
    for segment_index in range(T.size):
        if segment_index == 0:
            eta_l = zeta_start[1:s]
        else:
            eta_l = eta_corrected[segment_index - 1]
        if segment_index == T.size - 1:
            eta_r = zeta_end[1:s]
        else:
            eta_r = eta_corrected[segment_index]
        coeffs[segment_index] = _reconstruct_segment_coeff(
            q[segment_index],
            q[segment_index + 1],
            float(T[segment_index]),
            eta_l,
            eta_r,
            s=s,
        )

    jumps = compute_jumps(coeffs, T, s)
    stats = summarize_jumps(jumps)
    eta_delta_raw = eta_left - eta_right
    return {
        "scheme": mode,
        "coeffs": coeffs,
        "c_vec": coeffs.reshape(-1),
        "eta_corrected": eta_corrected,
        "eta_left_raw": eta_left,
        "eta_right_raw": eta_right,
        "eta_delta_raw": eta_delta_raw,
        "weights_left": weights_left,
        "weights_right": weights_right,
        "jumps": jumps,
        "stats": stats,
        "cost": compute_total_snap_cost(coeffs, T),
        "meta": {"requested_k": int(k), "is_k_dependent": True, "correction_mode": mode},
    }


def _curve_result_to_dict(result: dict[str, Any], curve: str) -> dict[str, dict[int, dict[str, Any]]]:
    out: dict[str, dict[int, dict[str, Any]]] = {curve: {}}
    if curve == "ideal":
        for record in result["ideal_records"]:
            out[curve][int(record["k"])] = record
    else:
        for record in result["actual_records"][curve]:
            out[curve][int(record["k"])] = record
    return out


def _run_case_curves(
    q: np.ndarray,
    T: np.ndarray,
    *,
    k_values: list[int] | None = None,
    boundary_trim: int = DEFAULT_BOUNDARY_TRIM,
) -> dict[str, Any]:
    q, T = _validate_problem_inputs(q, T, DEFAULT_S)
    if k_values is None:
        k_values = _k_values_for_M(T.size)
    reference = compute_minco_reference(q, T, s=DEFAULT_S)

    ideal_by_k: dict[int, dict[str, Any]] = {}
    actual_by_curve: dict[str, dict[int, dict[str, Any]]] = {curve: {} for curve in ("A", "B", "C", "C1", "C2")}

    for k in k_values:
        ideal = compute_ideal_truncated_blom_k(q, T, k, s=DEFAULT_S, reference=reference)
        ideal_metrics = compute_full_and_interior_errors(
            reference["c_vec"],
            ideal["c_vec"],
            cost_ref=reference["cost"],
            cost_test=ideal["cost"],
            boundary_trim=boundary_trim,
        )
        ideal_by_k[int(k)] = {**ideal_metrics, "k": int(k), "coeffs": ideal["coeffs"], "c_vec": ideal["c_vec"], "cost": ideal["cost"]}

        for scheme in ("A", "B", "C"):
            actual = compute_actual_blom_k(q, T, k, s=DEFAULT_S, scheme=scheme, config=None)
            metrics = compute_full_and_interior_errors(
                reference["c_vec"],
                actual["c_vec"],
                c_ideal=ideal["c_vec"],
                cost_ref=reference["cost"],
                cost_test=actual["cost"],
                boundary_trim=boundary_trim,
            )
            actual_by_curve[scheme][int(k)] = {**metrics, "k": int(k), **_serialize(actual)}

        for mode in ("C1", "C2"):
            corrected = assemble_scheme_C_light(q, T, k, mode=mode, s=DEFAULT_S)
            metrics = compute_full_and_interior_errors(
                reference["c_vec"],
                corrected["c_vec"],
                c_ideal=ideal["c_vec"],
                cost_ref=reference["cost"],
                cost_test=corrected["cost"],
                boundary_trim=boundary_trim,
            )
            actual_by_curve[mode][int(k)] = {
                **metrics,
                "k": int(k),
                "lower_order_jump_max": max((corrected["stats"][order]["max_abs"] for order in range(DEFAULT_S)), default=0.0),
                **_serialize(corrected),
            }

    fits = {
        "ideal_full": fit_log_error_vs_k(k_values, [ideal_by_k[k]["full_coef_error"] for k in k_values]),
        "ideal_interior": fit_log_error_vs_k(k_values, [ideal_by_k[k]["interior_coef_error"] for k in k_values]),
    }
    for curve in ("A", "B", "C", "C1", "C2"):
        fits[f"{curve}_full"] = fit_log_error_vs_k(k_values, [actual_by_curve[curve][k]["full_coef_error"] for k in k_values])
        fits[f"{curve}_interior"] = fit_log_error_vs_k(k_values, [actual_by_curve[curve][k]["interior_coef_error"] for k in k_values])

    return {
        "q": q,
        "T": T,
        "k_values": list(k_values),
        "reference": reference,
        "ideal": ideal_by_k,
        "actual": actual_by_curve,
        "fits": fits,
        "boundary_trim": boundary_trim,
    }


def _plot_lines(
    x_values: list[int] | np.ndarray,
    curves: list[tuple[str, list[float]]],
    save_path: str | Path,
    *,
    title: str,
    ylabel: str,
    logy: bool = False,
) -> None:
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    x = np.asarray(x_values, dtype=float)
    for label, values in curves:
        ax.plot(x, np.asarray(values, dtype=float), marker="o", linewidth=1.9, label=label)
    ax.set_title(title)
    ax.set_xlabel("k" if len(x.shape) == 1 else "index")
    ax.set_ylabel(ylabel)
    if logy:
        ax.set_yscale("log")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_multi_panel_error_by_M(
    records: dict[int, dict[str, Any]],
    save_path: str | Path,
    *,
    key: str,
    title: str,
    ylabel: str,
) -> None:
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 8.0), sharex=False)
    for M, record in sorted(records.items()):
        k_values = record["k_values"]
        axes[0].plot(k_values, [record["ideal"][k][key] for k in k_values], marker="o", linewidth=1.7, label=f"M={M}")
        axes[1].plot(k_values, [record["actual"]["C"][k][key] for k in k_values], marker="o", linewidth=1.7, label=f"M={M}")
    axes[0].set_title(f"{title}: ideal truncation")
    axes[1].set_title(f"{title}: Scheme C")
    for ax in axes:
        ax.set_xlabel("k")
        ax.set_ylabel(ylabel)
        ax.set_yscale("log")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_heatmap(
    mat: np.ndarray,
    save_path: str | Path,
    *,
    title: str,
    xlabel: str,
    ylabel: str,
) -> None:
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    display = np.log10(np.asarray(mat, dtype=float) + 1e-16)
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    image = ax.imshow(display, aspect="auto", origin="lower", cmap="magma")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("log10(value + 1e-16)")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_boxplot(groups: dict[str, list[float]], save_path: str | Path, *, title: str, ylabel: str, logy: bool = False) -> None:
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    labels = list(groups)
    values = [groups[label] for label in labels]
    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    ax.boxplot(values, tick_labels=labels)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    if logy:
        ax.set_yscale("log")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_bar(
    labels: list[str],
    values: list[float],
    save_path: str | Path,
    *,
    title: str,
    ylabel: str,
) -> None:
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.bar(x, values, color="#4c78a8", edgecolor="black")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_scatter(
    points: list[tuple[str, float, float]],
    save_path: str | Path,
    *,
    title: str,
    xlabel: str,
    ylabel: str,
    logy: bool = False,
) -> None:
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    for label, x_value, y_value in points:
        ax.scatter([x_value], [y_value], s=60, label=label)
        ax.annotate(label, (x_value, y_value), xytext=(4, 4), textcoords="offset points")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if logy:
        ax.set_yscale("log")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _supports_matching(best_slope: float, final_match_error: float) -> bool:
    return bool(np.isfinite(best_slope) and best_slope < -0.1 and np.isfinite(final_match_error) and final_match_error < 1.0)


def _supports_feasibility(final_cost_rel: float, jump_metric: float | None = None) -> bool:
    if not np.isfinite(final_cost_rel):
        return False
    if jump_metric is None:
        return abs(final_cost_rel) < 5e-2
    return abs(final_cost_rel) < 5e-2 and jump_metric < 1e-8


def _suggests_model_revision(best_slope: float, final_coef_error: float, final_match_error: float | None = None) -> bool:
    if not np.isfinite(best_slope):
        return True
    if best_slope >= -0.05:
        return True
    if final_coef_error > 10.0:
        return True
    if final_match_error is not None and np.isfinite(final_match_error) and final_match_error > 10.0:
        return True
    return False


def _safe_nanmedian(values: list[float]) -> float:
    if not values:
        return math.nan
    arr = np.asarray(values, dtype=float)
    if not np.any(np.isfinite(arr)):
        return math.nan
    return float(np.nanmedian(arr))


def run_large_M_sweep(
    *,
    M_values: list[int] | None = None,
    k_max: int = 20,
    seed: int = 42,
    save_dir: str | Path | None = None,
) -> dict[str, Any]:
    if M_values is None:
        M_values = [10, 20, 40, 80, 120]
    records: dict[int, dict[str, Any]] = {}
    rows: list[list[Any]] = []
    for M in M_values:
        rng = np.random.default_rng(seed + int(M))
        q = _default_q_sampler(rng, int(M))
        T = _default_T_sampler(rng, int(M))
        k_values = _k_values_for_M(int(M), k_max=k_max)
        result = _run_case_curves(q, T, k_values=k_values)
        records[int(M)] = result
        for curve in ("ideal", "A", "B", "C"):
            fit = result["fits"][f"{curve}_full"]
            final = result["ideal"][k_values[-1]] if curve == "ideal" else result["actual"][curve][k_values[-1]]
            rows.append(
                [
                    M,
                    curve,
                    fit["slope"],
                    fit["r2"],
                    final["full_coef_error"],
                    final.get("full_matching_error", math.nan),
                    final.get("rel_cost_gap", math.nan),
                ]
            )

    output = {"records": records, "summary_rows": rows}
    if save_dir is None:
        return output

    save_dir = Path(save_dir)
    _save_csv(
        save_dir / "summary_large_M.csv",
        ["M", "curve", "slope", "r2", "final_coef_error", "final_match_error", "final_cost_rel"],
        rows,
    )
    _save_json(
        save_dir / "summary_large_M.json",
        _serialize({"summary_rows": rows, "M_values": M_values}),
    )
    _plot_multi_panel_error_by_M(records, save_dir / "coef_error_vs_k_by_M.png", key="full_coef_error", title="Coefficient error vs k by M", ylabel="full coefficient error")
    _plot_multi_panel_error_by_M(records, save_dir / "logfit_vs_k_by_M.png", key="interior_coef_error", title="Interior coefficient error vs k by M", ylabel="interior coefficient error")

    labels = [str(M) for M in M_values]
    for metric_name, filename, title in (
        ("slope", "slope_vs_M.png", "Log-fit slope vs M"),
        ("final_match_error", "matching_vs_M.png", "Final matching error vs M"),
        ("final_cost_rel", "cost_gap_vs_M.png", "Final relative cost gap vs M"),
    ):
        curves = []
        for curve in ("ideal", "A", "B", "C"):
            values = []
            for M in M_values:
                result = records[M]
                fit = result["fits"][f"{curve}_full"]
                final = result["ideal"][result["k_values"][-1]] if curve == "ideal" else result["actual"][curve][result["k_values"][-1]]
                if metric_name == "slope":
                    values.append(fit["slope"])
                elif metric_name == "final_match_error":
                    values.append(final.get("full_matching_error", math.nan))
                else:
                    values.append(final.get("rel_cost_gap", math.nan))
            curves.append((curve, values))
        _plot_lines(M_values, curves, save_dir / filename, title=title, ylabel=metric_name.replace("_", " "), logy=(metric_name != "slope"))
    return output


def run_time_regime_split(
    *,
    M: int = 20,
    n_trials: int = 12,
    k_max: int = 20,
    seed: int = 42,
    save_dir: str | Path | None = None,
) -> dict[str, Any]:
    regimes: list[tuple[str, Callable[[np.random.Generator], np.ndarray]]] = []
    for h in (0.5, 1.0, 2.0):
        regimes.append((f"uniform_h={h:.1f}", lambda rng, h=h: np.full((M,), h, dtype=float)))
    for low, high in ((0.9, 1.1), (0.5, 2.0), (0.2, 3.0)):
        regimes.append((f"bounded_{low:.1f}_{high:.1f}", lambda rng, low=low, high=high: _default_T_sampler(rng, M, low=low, high=high)))

    k_values = _k_values_for_M(M, k_max=k_max)
    records: list[dict[str, Any]] = []
    error_curves_by_regime: dict[str, list[np.ndarray]] = {}
    for regime_index, (label, sampler) in enumerate(regimes):
        error_curves_by_regime[label] = []
        for trial in range(n_trials):
            rng = np.random.default_rng(seed + 100 * regime_index + trial)
            q = _default_q_sampler(rng, M)
            T = np.asarray(sampler(rng), dtype=float).reshape(M)
            result = _run_case_curves(q, T, k_values=k_values)
            error_curves_by_regime[label].append(np.asarray([result["actual"]["C"][k]["full_coef_error"] for k in k_values], dtype=float))
            for curve in ("ideal", "B", "C"):
                fit = result["fits"][f"{curve}_full"]
                final = result["ideal"][k_values[-1]] if curve == "ideal" else result["actual"][curve][k_values[-1]]
                records.append(
                    {
                        "regime": label,
                        "trial": trial,
                        "curve": curve,
                        "slope": fit["slope"],
                        "r2": fit["r2"],
                        "final_coef_error": final["full_coef_error"],
                        "final_match_error": final.get("full_matching_error", math.nan),
                        "final_cost_rel": final.get("rel_cost_gap", math.nan),
                    }
                )

    if save_dir is None:
        return {"records": records, "k_values": k_values}

    save_dir = Path(save_dir)
    _save_csv(
        save_dir / "summary_time_regime.csv",
        ["regime", "trial", "curve", "slope", "r2", "final_coef_error", "final_match_error", "final_cost_rel"],
        [
            [
                record["regime"],
                record["trial"],
                record["curve"],
                record["slope"],
                record["r2"],
                record["final_coef_error"],
                record["final_match_error"],
                record["final_cost_rel"],
            ]
            for record in records
        ],
    )
    _save_json(save_dir / "summary_time_regime.json", _serialize({"records": records, "k_values": k_values}))

    median_curves = []
    for label, mats in error_curves_by_regime.items():
        median_curves.append((label, np.median(np.vstack(mats), axis=0).tolist()))
    _plot_lines(k_values, median_curves, save_dir / "uniform_vs_bounded_coef_error.png", title="Uniform vs bounded coefficient error", ylabel="median Scheme C full coefficient error", logy=True)

    slope_groups = {label: [record["slope"] for record in records if record["regime"] == label and record["curve"] == "C"] for label, _ in regimes}
    _plot_boxplot(slope_groups, save_dir / "uniform_vs_bounded_slope_boxplot.png", title="Scheme C slope distribution by time regime", ylabel="fitted slope", logy=False)

    match_groups = {label: [record["final_match_error"] for record in records if record["regime"] == label and record["curve"] == "C"] for label, _ in regimes}
    _plot_boxplot(match_groups, save_dir / "time_regime_matching_compare.png", title="Scheme C final matching error by time regime", ylabel="final matching error", logy=True)

    cost_groups = {label: [max(record["final_cost_rel"], 1e-15) for record in records if record["regime"] == label and record["curve"] == "C"] for label, _ in regimes}
    _plot_boxplot(cost_groups, save_dir / "time_regime_cost_gap_compare.png", title="Scheme C final relative cost gap by time regime", ylabel="final relative cost gap", logy=True)
    return {"records": records, "k_values": k_values}


def run_interior_vs_full(
    *,
    q: np.ndarray | None = None,
    T: np.ndarray | None = None,
    M: int = 24,
    boundary_trim: int = DEFAULT_BOUNDARY_TRIM,
    n_trials: int = 16,
    seed: int = 42,
    save_dir: str | Path | None = None,
) -> dict[str, Any]:
    if q is None or T is None:
        rng = np.random.default_rng(seed)
        q = _default_q_sampler(rng, M)
        T = _default_T_sampler(rng, M)
    q, T = _validate_problem_inputs(q, T, DEFAULT_S)
    k_values = _k_values_for_M(T.size)
    result = _run_case_curves(q, T, k_values=k_values, boundary_trim=boundary_trim)

    rows = []
    for curve in ("ideal", "A", "B", "C"):
        for k in k_values:
            record = result["ideal"][k] if curve == "ideal" else result["actual"][curve][k]
            rows.append(
                [
                    curve,
                    k,
                    record["full_coef_error"],
                    record["interior_coef_error"],
                    record.get("full_matching_error", math.nan),
                    record.get("interior_matching_error", math.nan),
                ]
            )

    random_records: list[dict[str, Any]] = []
    for trial in range(n_trials):
        rng = np.random.default_rng(seed + 1000 + trial)
        q_trial = _default_q_sampler(rng, T.size)
        T_trial = _default_T_sampler(rng, T.size)
        case = _run_case_curves(q_trial, T_trial, k_values=k_values, boundary_trim=boundary_trim)
        final_k = k_values[-1]
        for curve in ("ideal", "B", "C"):
            record = case["ideal"][final_k] if curve == "ideal" else case["actual"][curve][final_k]
            random_records.append(
                {
                    "curve": curve,
                    "boundary_coef_error": record["boundary_coef_error"],
                    "interior_coef_error": record["interior_coef_error"],
                }
            )

    if save_dir is None:
        return {"case": result, "summary_rows": rows, "random_records": random_records}

    save_dir = Path(save_dir)
    _save_csv(
        save_dir / "summary_interior_vs_full.csv",
        ["curve", "k", "full_coef_error", "interior_coef_error", "full_matching_error", "interior_matching_error"],
        rows,
    )
    _save_json(save_dir / "summary_interior_vs_full.json", _serialize({"summary_rows": rows, "boundary_trim": boundary_trim}))

    curves = []
    for curve in ("ideal", "A", "B", "C"):
        curves.append((f"{curve} full", [(result["ideal"][k] if curve == "ideal" else result["actual"][curve][k])["full_coef_error"] for k in k_values]))
        curves.append((f"{curve} interior", [(result["ideal"][k] if curve == "ideal" else result["actual"][curve][k])["interior_coef_error"] for k in k_values]))
    _plot_lines(k_values, curves, save_dir / "full_vs_interior_coef_error.png", title="Full vs interior coefficient error", ylabel="coefficient error", logy=True)

    labels = []
    values = []
    for curve in ("ideal", "A", "B", "C"):
        labels.extend([f"{curve} full", f"{curve} int"])
        values.extend([result["fits"][f"{curve}_full"]["slope"], result["fits"][f"{curve}_interior"]["slope"]])
    _plot_bar(labels, values, save_dir / "full_vs_interior_logfit.png", title="Full vs interior fitted slope", ylabel="slope")

    heatmap = np.asarray([result["actual"]["C"][k]["segmentwise_l2"] for k in k_values], dtype=float)
    _plot_heatmap(heatmap, save_dir / "segmentwise_error_heatmap.png", title="Scheme C segmentwise error heatmap", xlabel="segment index", ylabel="k sweep index")

    groups = {}
    for curve in ("ideal", "B", "C"):
        groups[f"{curve} boundary"] = [record["boundary_coef_error"] for record in random_records if record["curve"] == curve]
        groups[f"{curve} interior"] = [record["interior_coef_error"] for record in random_records if record["curve"] == curve]
    _plot_boxplot(groups, save_dir / "boundary_vs_interior_boxplot.png", title="Boundary vs interior error distribution", ylabel="error", logy=True)
    return {"case": result, "summary_rows": rows, "random_records": random_records}


def _scheme_c_q_bandwidth(
    q: np.ndarray,
    T: np.ndarray,
    assembler: Callable[[np.ndarray], np.ndarray],
    *,
    eps: float = 1e-6,
) -> dict[str, Any]:
    q = np.asarray(q, dtype=float).reshape(-1)
    T = _validate_time_vector(T)
    q0 = float(q[0])
    qM = float(q[-1])
    q_interior = q[1:-1].copy()

    def coeff_map(q_inner: np.ndarray) -> np.ndarray:
        q_full = np.concatenate(([q0], np.asarray(q_inner, dtype=float).reshape(-1), [qM]))
        return np.asarray(assembler(q_full), dtype=float).reshape(-1)

    J_q = finite_difference_jacobian(coeff_map, q_interior, eps=eps, method="central")
    return estimate_effective_bandwidth(J_q, T.size, tol=1e-8, block_size=2 * DEFAULT_S)


def run_schemeC_light_assembly(
    *,
    q: np.ndarray | None = None,
    T: np.ndarray | None = None,
    M: int = 20,
    seed: int = 42,
    save_dir: str | Path | None = None,
) -> dict[str, Any]:
    if q is None or T is None:
        rng = np.random.default_rng(seed)
        q = _default_q_sampler(rng, M)
        T = _default_T_sampler(rng, M)
    q, T = _validate_problem_inputs(q, T, DEFAULT_S)
    k_values = _k_values_for_M(T.size)
    reference = compute_minco_reference(q, T, s=DEFAULT_S)

    ideal_records: dict[int, dict[str, Any]] = {}
    curve_records: dict[str, dict[int, dict[str, Any]]] = {curve: {} for curve in ("C", "C1", "C2")}
    for k in k_values:
        ideal = compute_ideal_truncated_blom_k(q, T, k, s=DEFAULT_S, reference=reference)
        ideal_records[k] = compute_full_and_interior_errors(reference["c_vec"], ideal["c_vec"], cost_ref=reference["cost"], cost_test=ideal["cost"])
        ideal_records[k]["k"] = k

        raw = compute_actual_blom_k(q, T, k, s=DEFAULT_S, scheme="C", config=None)
        raw_metrics = compute_full_and_interior_errors(reference["c_vec"], raw["c_vec"], c_ideal=ideal["c_vec"], cost_ref=reference["cost"], cost_test=raw["cost"])
        raw_metrics["lower_order_jump_max"] = max((raw["stats"][order]["max_abs"] for order in range(DEFAULT_S)), default=0.0)
        curve_records["C"][k] = {**raw_metrics, **_serialize(raw), "k": k}

        for mode in ("C1", "C2"):
            corrected = assemble_scheme_C_light(q, T, k, mode=mode, s=DEFAULT_S)
            corrected_metrics = compute_full_and_interior_errors(reference["c_vec"], corrected["c_vec"], c_ideal=ideal["c_vec"], cost_ref=reference["cost"], cost_test=corrected["cost"])
            corrected_metrics["lower_order_jump_max"] = max((corrected["stats"][order]["max_abs"] for order in range(DEFAULT_S)), default=0.0)
            curve_records[mode][k] = {**corrected_metrics, **_serialize(corrected), "k": k}

    locality_points: list[tuple[str, float, float]] = []
    for curve in ("C", "C1", "C2"):
        final_k = k_values[-1]
        if curve == "C":
            assembler = lambda q_full, final_k=final_k: compute_actual_blom_k(q_full, T, final_k, s=DEFAULT_S, scheme="C", config=None)["coeffs"]
        else:
            assembler = lambda q_full, curve=curve, final_k=final_k: assemble_scheme_C_light(q_full, T, final_k, mode=curve, s=DEFAULT_S)["coeffs"]
        bandwidth = _scheme_c_q_bandwidth(q, T, lambda q_full, assembler=assembler: assembler(q_full))
        curve_records[curve][final_k]["q_bandwidth"] = bandwidth
        locality_points.append((curve, bandwidth["max_effective_bandwidth"], curve_records[curve][final_k]["full_coef_error"]))

    fits = {curve: fit_log_error_vs_k(k_values, [curve_records[curve][k]["full_coef_error"] for k in k_values]) for curve in ("C", "C1", "C2")}
    rows = []
    for curve in ("C", "C1", "C2"):
        for k in k_values:
            record = curve_records[curve][k]
            rows.append(
                [
                    curve,
                    k,
                    record["full_coef_error"],
                    record["full_matching_error"],
                    record["rel_cost_gap"],
                    record["lower_order_jump_max"],
                ]
            )

    if save_dir is None:
        return {"ideal": ideal_records, "curves": curve_records, "fits": fits, "rows": rows}

    save_dir = Path(save_dir)
    _save_csv(
        save_dir / "summary_schemeC_correction.csv",
        ["curve", "k", "full_coef_error", "full_matching_error", "rel_cost_gap", "lower_order_jump_max"],
        rows,
    )
    _save_json(save_dir / "summary_schemeC_correction.json", _serialize({"rows": rows, "fits": fits}))

    _plot_lines(
        k_values,
        [
            ("ideal", [ideal_records[k]["full_coef_error"] for k in k_values]),
            ("raw C", [curve_records["C"][k]["full_coef_error"] for k in k_values]),
            ("C1", [curve_records["C1"][k]["full_coef_error"] for k in k_values]),
            ("C2", [curve_records["C2"][k]["full_coef_error"] for k in k_values]),
        ],
        save_dir / "schemeC_raw_vs_corrected_error.png",
        title="Scheme C raw vs corrected coefficient error",
        ylabel="coefficient error",
        logy=True,
    )
    _plot_lines(
        k_values,
        [
            ("raw C", [curve_records["C"][k]["full_matching_error"] for k in k_values]),
            ("C1", [curve_records["C1"][k]["full_matching_error"] for k in k_values]),
            ("C2", [curve_records["C2"][k]["full_matching_error"] for k in k_values]),
        ],
        save_dir / "schemeC_raw_vs_corrected_matching.png",
        title="Scheme C raw vs corrected matching error",
        ylabel="matching error",
        logy=True,
    )
    _plot_lines(
        k_values,
        [
            ("raw C", [curve_records["C"][k]["lower_order_jump_max"] for k in k_values]),
            ("C1", [curve_records["C1"][k]["lower_order_jump_max"] for k in k_values]),
            ("C2", [curve_records["C2"][k]["lower_order_jump_max"] for k in k_values]),
        ],
        save_dir / "schemeC_raw_vs_corrected_jumps.png",
        title="Scheme C raw vs corrected lower-order jumps",
        ylabel="max lower-order jump",
        logy=True,
    )
    _plot_scatter(
        locality_points,
        save_dir / "schemeC_raw_vs_corrected_locality_tradeoff.png",
        title="Scheme C locality tradeoff",
        xlabel="q effective bandwidth",
        ylabel="final coefficient error",
        logy=True,
    )
    return {"ideal": ideal_records, "curves": curve_records, "fits": fits, "rows": rows}


def _build_overview_rows(exp1: dict[str, Any], exp2: dict[str, Any], exp3: dict[str, Any], exp4: dict[str, Any]) -> list[list[Any]]:
    rows: list[list[Any]] = []

    for curve in ("ideal", "A", "B", "C"):
        matching = [row[5] for row in exp1["summary_rows"] if row[1] == curve]
        slopes = [row[2] for row in exp1["summary_rows"] if row[1] == curve]
        r2_vals = [row[3] for row in exp1["summary_rows"] if row[1] == curve]
        coef_vals = [row[4] for row in exp1["summary_rows"] if row[1] == curve]
        cost_vals = [row[6] for row in exp1["summary_rows"] if row[1] == curve]
        best_slope = _safe_nanmedian(slopes)
        best_r2 = _safe_nanmedian(r2_vals)
        final_coef = _safe_nanmedian(coef_vals)
        final_match = _safe_nanmedian(matching)
        final_cost = _safe_nanmedian(cost_vals)
        rows.append(
            [
                "exp1_large_M_sweep",
                curve,
                best_slope,
                best_r2,
                final_coef,
                final_match,
                final_cost,
                _supports_matching(best_slope, final_match),
                _supports_feasibility(final_cost),
                _suggests_model_revision(best_slope, final_coef, final_match),
            ]
        )

    records2 = exp2["records"]
    for curve in ("ideal", "B", "C"):
        subset = [record for record in records2 if record["curve"] == curve]
        best_slope = _safe_nanmedian([record["slope"] for record in subset])
        best_r2 = _safe_nanmedian([record["r2"] for record in subset])
        final_coef = _safe_nanmedian([record["final_coef_error"] for record in subset])
        final_match = _safe_nanmedian([record["final_match_error"] for record in subset])
        final_cost = _safe_nanmedian([record["final_cost_rel"] for record in subset])
        rows.append(
            [
                "exp2_time_regime_split",
                curve,
                best_slope,
                best_r2,
                final_coef,
                final_match,
                final_cost,
                _supports_matching(best_slope, final_match),
                _supports_feasibility(final_cost),
                _suggests_model_revision(best_slope, final_coef, final_match),
            ]
        )

    case3 = exp3["case"]
    for curve in ("ideal", "A", "B", "C"):
        fit = case3["fits"][f"{curve}_interior"]
        final_k = case3["k_values"][-1]
        final = case3["ideal"][final_k] if curve == "ideal" else case3["actual"][curve][final_k]
        rows.append(
            [
                "exp3_interior_vs_full",
                curve,
                fit["slope"],
                fit["r2"],
                final["interior_coef_error"],
                final.get("interior_matching_error", math.nan),
                final.get("rel_cost_gap", math.nan),
                _supports_matching(fit["slope"], final.get("interior_matching_error", math.nan)),
                _supports_feasibility(final.get("rel_cost_gap", math.nan)),
                _suggests_model_revision(fit["slope"], final["interior_coef_error"], final.get("interior_matching_error", math.nan)),
            ]
        )

    final_k = max(exp4["curves"]["C"])
    for curve in ("C", "C1", "C2"):
        fit = exp4["fits"][curve]
        final = exp4["curves"][curve][final_k]
        rows.append(
            [
                "exp4_schemeC_light_assembly",
                curve,
                fit["slope"],
                fit["r2"],
                final["full_coef_error"],
                final["full_matching_error"],
                final["rel_cost_gap"],
                _supports_matching(fit["slope"], final["full_matching_error"]),
                _supports_feasibility(final["rel_cost_gap"], final["lower_order_jump_max"]),
                _suggests_model_revision(fit["slope"], final["full_coef_error"], final["full_matching_error"]),
            ]
        )
    return rows


def _write_interpretation_summary(rows: list[list[Any]], save_path: str | Path) -> None:
    def _subset(experiment: str) -> list[list[Any]]:
        return [row for row in rows if row[0] == experiment]

    exp1 = _subset("exp1_large_M_sweep")
    exp2 = _subset("exp2_time_regime_split")
    exp3 = _subset("exp3_interior_vs_full")
    exp4 = _subset("exp4_schemeC_light_assembly")

    def _best(subset: list[list[Any]]) -> str:
        if not subset:
            return "n/a"
        return min(subset, key=lambda row: float(row[4]))[1]

    lines = [
        "# Phase 7 Extra Interpretation Summary",
        "",
        "## Large-M Sweep",
        f"- Best curve by median final coefficient error: `{_best(exp1)}`.",
        "- If Scheme C keeps a stable negative slope while M grows, the matching-theorem path remains numerically plausible.",
        "",
        "## Time-Regime Split",
        f"- Best curve by median final coefficient error: `{_best(exp2)}`.",
        "- Strong degradation under bounded-nonuniform regimes would suggest shrinking the next theorem to the uniform-time regime first.",
        "",
        "## Interior vs Full Error",
        f"- Best interior-focused curve by median interior coefficient error: `{_best(exp3)}`.",
        "- If interior errors are much smaller than full errors, boundary effects likely dominate and an interior theorem becomes a sensible next milestone.",
        "",
        "## Scheme C Light Assembly",
        f"- Best corrected Scheme C variant by final coefficient error: `{_best(exp4)}`.",
        "- If a very light correction sharply improves matching while preserving low locality width, it is a strong candidate for promotion to the formal BLOM object.",
        "",
        "## Research Direction",
        "- Continue toward the matching theorem if ideal truncation remains strong, Scheme C or corrected Scheme C keeps a stable negative slope, and final matching errors become small.",
        "- Prioritize feasibility / cost theorems only when cost gaps are already consistently near zero under the same candidates.",
        "- Revisit the BLOM definition or assembly if large-M scaling degrades sharply or if corrected Scheme C still cannot close the ideal-vs-actual gap.",
    ]
    Path(save_path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_phase7_extra_experiments(
    q: np.ndarray | None = None,
    T: np.ndarray | None = None,
    M_values: list[int] | None = None,
    n_trials: int = 100,
    s: int = DEFAULT_S,
    k_values: list[int] | None = None,
    schemes: tuple[str, ...] = ("A", "B", "C"),
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    del k_values, schemes
    if s != DEFAULT_S:
        raise ValueError("Phase 7 extra experiments currently support only s=4.")
    if save_dir is None:
        save_dir = DEFAULT_RESULTS_DIR
    dirs = ensure_results_dirs(save_dir)

    exp1 = run_large_M_sweep(M_values=M_values, seed=seed, save_dir=dirs["exp1"])
    exp2 = run_time_regime_split(n_trials=min(n_trials, 16), seed=seed, save_dir=dirs["exp2"])
    exp3 = run_interior_vs_full(q=q, T=T, n_trials=min(n_trials, 24), seed=seed, save_dir=dirs["exp3"])
    exp4 = run_schemeC_light_assembly(q=q, T=T, seed=seed, save_dir=dirs["exp4"])

    overview_rows = _build_overview_rows(exp1, exp2, exp3, exp4)
    _save_csv(
        dirs["compare"] / "final_experiment_overview.csv",
        [
            "experiment_name",
            "scheme",
            "best_slope",
            "best_r2",
            "final_coef_error",
            "final_match_error",
            "final_cost_rel",
            "supports_matching_theorem",
            "supports_feasibility_theorem",
            "suggests_model_revision",
        ],
        overview_rows,
    )
    _save_json(
        dirs["compare"] / "phase7_extra_overview.json",
        _serialize({"overview_rows": overview_rows}),
    )
    _write_interpretation_summary(overview_rows, dirs["compare"] / "phase7_extra_interpretation_summary.md")
    return {
        "exp1": exp1,
        "exp2": exp2,
        "exp3": exp3,
        "exp4": exp4,
        "overview_rows": overview_rows,
        "paths": {
            "overview_csv": dirs["compare"] / "final_experiment_overview.csv",
            "overview_json": dirs["compare"] / "phase7_extra_overview.json",
            "summary_md": dirs["compare"] / "phase7_extra_interpretation_summary.md",
        },
    }


def main() -> None:
    q, T = representative_case()
    result = run_phase7_extra_experiments(q=q, T=T, M_values=[10, 20, 40], n_trials=8, save_dir=DEFAULT_RESULTS_DIR, seed=42)
    print("Phase 7 extra experiments")
    print(f"overview rows: {len(result['overview_rows'])}")
    print(f"overview csv: {result['paths']['overview_csv']}")


if __name__ == "__main__":
    main()
