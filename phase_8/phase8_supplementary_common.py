"""
Shared helpers for Phase 8 supplementary experiments.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from phase_7.blom_convergence_vs_k import compute_actual_blom_k, compute_ideal_truncated_blom_k, compute_minco_reference
from phase_8.phase8_common import (
    DEFAULT_S,
    compute_matching_error_blocks,
    default_q_sampler,
    extract_reference_gamma,
    fit_log_error_vs_k,
    make_interior_sets,
    sample_bounded_nonuniform_time,
    sample_uniform_time,
    save_csv,
    save_json,
    serialize,
    solve_augmented_local_window,
    summarize_interior_matching,
    validate_problem_inputs,
)


DEFAULT_SUPPLEMENTARY_RESULTS_DIR = Path("phase_8/results/phase8_supplementary")
DEFAULT_LARGE_M_VALUES = [20, 40, 80, 120]
DEFAULT_K_VALUES = [2, 4, 6, 8, 10, 12]
DEFAULT_RADIUS_MODES = ["half_k", "k", "three_half_k"]
EPSILON = 1e-12


def ensure_supplementary_results_dirs(base_dir: str | Path = DEFAULT_SUPPLEMENTARY_RESULTS_DIR) -> dict[str, Path]:
    base = Path(base_dir)
    directories = {
        "base": base,
        "exp1_large_M": base / "exp1_large_M",
        "exp2_more_trials": base / "exp2_more_trials",
        "exp3_radius_sensitivity": base / "exp3_radius_sensitivity",
        "exp4_raw_vs_reference_sanity": base / "exp4_raw_vs_reference_sanity",
        "exp5_empty_interior_risk": base / "exp5_empty_interior_risk",
        "exp6_two_bridge_gaps": base / "exp6_two_bridge_gaps",
        "summary": base / "summary",
    }
    for directory in directories.values():
        directory.mkdir(parents=True, exist_ok=True)
    return directories


def normalize_radius_mode(mode: str) -> str:
    key = str(mode).strip().lower()
    aliases = {
        "half_k": "half",
        "half": "half",
        "k_over_2": "half",
        "minimal": "half",
        "k": "default",
        "default": "default",
        "three_half_k": "one_and_half",
        "three_halves": "one_and_half",
        "1.5k": "one_and_half",
        "one_and_half": "one_and_half",
    }
    if key not in aliases:
        raise ValueError(f"Unsupported supplementary radius mode {mode!r}.")
    return aliases[key]


def radius_mode_label(mode: str) -> str:
    key = normalize_radius_mode(mode)
    labels = {
        "half": "r(k)=ceil(k/2)",
        "default": "r(k)=k",
        "one_and_half": "r(k)=ceil(3k/2)",
    }
    return labels[key]


def build_case(
    M: int,
    *,
    regime: str,
    rng: np.random.Generator,
    h: float = 1.0,
    nonuniform_box: tuple[float, float] = (0.5, 2.0),
) -> tuple[np.ndarray, np.ndarray]:
    q = default_q_sampler(rng, M)
    if regime == "uniform":
        T = sample_uniform_time(M, h)
    elif regime == "bounded_nonuniform":
        T = sample_bounded_nonuniform_time(M, float(nonuniform_box[0]), float(nonuniform_box[1]), rng)
    else:
        raise ValueError(f"Unsupported regime {regime!r}.")
    return validate_problem_inputs(q, T)


def compute_reference_window_coefficients_fast(
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
    for seg_idx in range(1, T.size + 1):
        gamma_data = extract_reference_gamma(reference, T, seg_idx, k, s=s)
        solved = solve_augmented_local_window(q, T, seg_idx, k, gamma=gamma_data["gamma"], s=s)
        coeffs[seg_idx - 1] = solved["center_coeff"]
        gamma_norms[seg_idx - 1] = float(np.linalg.norm(gamma_data["gamma"]))
    return {
        "kind": "reference_window_fast",
        "k": int(k),
        "coeffs": coeffs,
        "c_vec": coeffs.reshape(-1),
        "gamma_norms": gamma_norms,
    }


def compute_case_metrics(
    q: np.ndarray,
    T: np.ndarray,
    k_values: list[int],
    *,
    radius_mode: str = "k",
    include_reference: bool = False,
    s: int = DEFAULT_S,
) -> list[dict[str, Any]]:
    q, T = validate_problem_inputs(q, T, s=s)
    reference = compute_minco_reference(q, T, s=s)
    M = T.size
    normalized_mode = normalize_radius_mode(radius_mode)
    records: list[dict[str, Any]] = []
    for k in k_values:
        actual = compute_actual_blom_k(q, T, k, s=s, scheme="C")
        ideal = compute_ideal_truncated_blom_k(q, T, k, s=s, reference=reference)
        summary = summarize_interior_matching(actual["c_vec"], ideal["c_vec"], M, k, s=s, radius_mode=normalized_mode)
        sets = make_interior_sets(M, k, radius_mode=normalized_mode)
        record: dict[str, Any] = {
            "k": int(k),
            "r_k": int(summary["r_k"]),
            "M": int(M),
            "radius_mode": radius_mode,
            "full_matching_l2": float(summary["full_l2"]),
            "interior_matching_l2": float(summary["interior_l2"]),
            "boundary_matching_l2": float(summary["boundary_l2"]),
            "boundary_energy_ratio": float(summary["boundary_energy_ratio"]),
            "interior_energy_ratio": float(summary["interior_energy_ratio"]),
            "segment_errors": np.asarray(summary["segment_errors"], dtype=float),
            "interior_idx": list(summary["interior_idx"]),
            "boundary_idx": list(summary["boundary_idx"]),
            "interior_count": len(summary["interior_idx"]),
            "boundary_count": len(summary["boundary_idx"]),
            "interior_fraction": len(summary["interior_idx"]) / max(M, 1),
            "is_empty_interior": len(summary["interior_idx"]) == 0,
            "reference": reference,
            "actual": actual,
            "ideal": ideal,
            "sets": sets,
        }
        if include_reference:
            reference_window = compute_reference_window_coefficients_fast(q, T, k, reference=reference, s=s)
            raw_to_ref_blocks = compute_matching_error_blocks(actual["c_vec"], reference_window["c_vec"], s=s)
            ref_to_ideal_blocks = compute_matching_error_blocks(reference_window["c_vec"], ideal["c_vec"], s=s)
            raw_to_ref = float(np.linalg.norm(raw_to_ref_blocks))
            ref_to_ideal = float(np.linalg.norm(ref_to_ideal_blocks))
            record.update(
                {
                    "reference_window": reference_window,
                    "raw_to_ref_l2": raw_to_ref,
                    "ref_to_ideal_l2": ref_to_ideal,
                    "gap_ratio": ref_to_ideal / max(raw_to_ref, EPSILON),
                    "raw_to_ref_segment": raw_to_ref_blocks,
                    "ref_to_ideal_segment": ref_to_ideal_blocks,
                }
            )
        records.append(record)
    return records


def compute_fit_row(records: list[dict[str, Any]]) -> dict[str, Any]:
    k_values = [int(record["k"]) for record in records]
    full_fit = fit_log_error_vs_k(k_values, [float(record["full_matching_l2"]) for record in records])
    interior_fit = fit_log_error_vs_k(k_values, [float(record["interior_matching_l2"]) for record in records])
    boundary_fit = fit_log_error_vs_k(k_values, [float(record["boundary_matching_l2"]) for record in records])
    return {
        "full_slope": float(full_fit["slope"]),
        "full_r2": float(full_fit["r2"]),
        "interior_slope": float(interior_fit["slope"]),
        "interior_r2": float(interior_fit["r2"]),
        "boundary_slope": float(boundary_fit["slope"]),
        "boundary_r2": float(boundary_fit["r2"]),
    }


def aggregate_numeric_records(
    records: list[dict[str, Any]],
    *,
    group_keys: list[str],
    metric_keys: list[str],
) -> list[dict[str, Any]]:
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for record in records:
        key = tuple(record[name] for name in group_keys)
        grouped.setdefault(key, []).append(record)
    rows: list[dict[str, Any]] = []
    for key, items in sorted(grouped.items(), key=lambda item: item[0]):
        row = {name: value for name, value in zip(group_keys, key)}
        for metric in metric_keys:
            values = np.asarray([float(item[metric]) for item in items], dtype=float)
            row[f"{metric}_mean"] = float(np.mean(values))
            row[f"{metric}_median"] = float(np.median(values))
            row[f"{metric}_std"] = float(np.std(values))
            row[f"{metric}_q25"] = float(np.quantile(values, 0.25))
            row[f"{metric}_q75"] = float(np.quantile(values, 0.75))
        rows.append(row)
    return rows


def line_plot(
    series: dict[str, tuple[list[float], list[float]]],
    *,
    xlabel: str,
    ylabel: str,
    title: str,
    save_path: str | Path,
    ylim: tuple[float, float] | None = None,
) -> None:
    plt.figure(figsize=(8.8, 5.0))
    for label, (xs, ys) in series.items():
        plt.plot(xs, ys, marker="o", label=label)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if ylim is not None:
        plt.ylim(*ylim)
    plt.title(title)
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def band_plot(
    series: dict[str, tuple[list[float], list[float], list[float]]],
    *,
    xlabel: str,
    ylabel: str,
    title: str,
    save_path: str | Path,
) -> None:
    plt.figure(figsize=(8.8, 5.0))
    for label, (xs, mean_vals, std_vals) in series.items():
        xs_arr = np.asarray(xs, dtype=float)
        mean_arr = np.asarray(mean_vals, dtype=float)
        std_arr = np.asarray(std_vals, dtype=float)
        plt.plot(xs_arr, mean_arr, marker="o", label=label)
        plt.fill_between(xs_arr, mean_arr - std_arr, mean_arr + std_arr, alpha=0.18)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def boxplot_by_group(
    grouped_values: dict[str, list[float]],
    *,
    ylabel: str,
    title: str,
    save_path: str | Path,
    rotate_xticks: bool = True,
) -> None:
    labels = list(grouped_values)
    values = [grouped_values[label] for label in labels]
    plt.figure(figsize=(9.0, 4.8))
    plt.boxplot(values, tick_labels=labels)
    if rotate_xticks:
        plt.xticks(rotation=20, ha="right")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def heatmap(
    matrix: np.ndarray,
    *,
    x_labels: list[str],
    y_labels: list[str],
    title: str,
    colorbar_label: str,
    save_path: str | Path,
    cmap: str = "magma",
) -> None:
    plt.figure(figsize=(9.2, 4.8))
    im = plt.imshow(np.asarray(matrix, dtype=float), aspect="auto", cmap=cmap)
    plt.colorbar(im, label=colorbar_label)
    plt.xticks(np.arange(len(x_labels)), x_labels, rotation=20, ha="right")
    plt.yticks(np.arange(len(y_labels)), y_labels)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()


def save_markdown(path: str | Path, lines: list[str]) -> None:
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def lightweight_record(record: dict[str, Any]) -> dict[str, Any]:
    keep = {}
    for key, value in record.items():
        if key in {"reference", "actual", "ideal", "reference_window", "sets"}:
            continue
        keep[key] = serialize(value)
    return keep


def representative_segment_index(M: int) -> int:
    return max(1, min(M, int(math.ceil(0.5 * M))))


__all__ = [
    "DEFAULT_K_VALUES",
    "DEFAULT_LARGE_M_VALUES",
    "DEFAULT_RADIUS_MODES",
    "DEFAULT_SUPPLEMENTARY_RESULTS_DIR",
    "EPSILON",
    "aggregate_numeric_records",
    "band_plot",
    "boxplot_by_group",
    "build_case",
    "compute_case_metrics",
    "compute_fit_row",
    "compute_reference_window_coefficients_fast",
    "ensure_supplementary_results_dirs",
    "heatmap",
    "lightweight_record",
    "line_plot",
    "normalize_radius_mode",
    "radius_mode_label",
    "representative_segment_index",
    "save_csv",
    "save_json",
    "save_markdown",
]
