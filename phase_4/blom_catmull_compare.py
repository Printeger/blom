"""
Comparison between exact BLOM analytic systems and Catmull-style heuristics.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from phase_4.blom_k2_s2_sympy import _catmull_velocity_numeric, _exact_velocity_numeric
from phase_4.blom_k2_s4_numeric import (
    compare_with_blom_strict_qp,
    hermite_reconstruct_center_segment,
    local_objective_s4_k2,
    solve_local_system_s4_k2,
)
from phase_4.utils.io_utils import ensure_dir, save_csv, save_json
from phase_4.utils.plotting_utils import save_histogram, save_overlay_curves, save_scatter_with_diagonal


DEFAULT_RESULTS_DIR = Path("phase_4/results/catmull_compare")


def catmull_velocity_time_weighted(q_im1: float, q_i: float, q_ip1: float, T_i: float, T_ip1: float) -> float:
    return _catmull_velocity_numeric(q_im1, q_i, q_ip1, T_i, T_ip1)


def heuristic_local_state_s4(
    q_im2: float,
    q_im1: float,
    q_i: float,
    q_ip1: float,
    T_im1: float,
    T_i: float,
    T_ip1: float,
) -> np.ndarray:
    v_left = catmull_velocity_time_weighted(q_im2, q_im1, q_i, T_im1, T_i)
    v_right = catmull_velocity_time_weighted(q_im1, q_i, q_ip1, T_i, T_ip1)
    return np.asarray([v_left, 0.0, 0.0, v_right, 0.0, 0.0], dtype=float)


def compare_s2_exact_vs_catmull(q_im1: float, q_i: float, q_ip1: float, T_i: float, T_ip1: float) -> dict[str, Any]:
    v_exact = _exact_velocity_numeric(q_im1, q_i, q_ip1, T_i, T_ip1)
    v_catmull = catmull_velocity_time_weighted(q_im1, q_i, q_ip1, T_i, T_ip1)
    return {
        "v_exact": v_exact,
        "v_catmull": v_catmull,
        "abs_error": abs(v_exact - v_catmull),
    }


def compare_s4_exact_vs_catmull(
    q_im2: float,
    q_im1: float,
    q_i: float,
    q_ip1: float,
    T_im1: float,
    T_i: float,
    T_ip1: float,
) -> dict[str, Any]:
    exact = solve_local_system_s4_k2(q_im2, q_im1, q_i, q_ip1, T_im1, T_i, T_ip1)
    x_exact = exact["x_opt"]
    x_heuristic = heuristic_local_state_s4(q_im2, q_im1, q_i, q_ip1, T_im1, T_i, T_ip1)
    c_exact = hermite_reconstruct_center_segment(x_exact, q_im1, q_i, T_i)["coeff"]
    c_heuristic = hermite_reconstruct_center_segment(x_heuristic, q_im1, q_i, T_i)["coeff"]
    J_exact = local_objective_s4_k2(x_exact, exact)
    J_heuristic = local_objective_s4_k2(x_heuristic, exact)
    return {
        "x_exact": x_exact,
        "x_heuristic": x_heuristic,
        "c_exact": c_exact,
        "c_heuristic": c_heuristic,
        "state_error_norm": float(np.linalg.norm(x_exact - x_heuristic)),
        "coeff_error_norm": float(np.linalg.norm(c_exact - c_heuristic)),
        "objective_gap": float(J_heuristic - J_exact),
        "v_exact": np.asarray([x_exact[0], x_exact[3]], dtype=float),
        "v_catmull": np.asarray([x_heuristic[0], x_heuristic[3]], dtype=float),
    }


def _time_imbalance_metric(T_im1: float, T_i: float, T_ip1: float) -> float:
    values = np.asarray([T_im1, T_i, T_ip1], dtype=float)
    return float(np.max(values) / np.min(values))


def run_random_benchmark(num_samples: int = 180, seed: int = 42) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    rows = []
    s2_exact_values = []
    s2_catmull_values = []
    s4_exact_vel = []
    s4_catmull_vel = []
    coeff_errors = []
    objective_gaps = []
    imbalance = []

    geometry_modes = ["near_linear", "monotone", "sharp_turn", "asymmetric"]
    for sample_index in range(num_samples):
        mode = geometry_modes[sample_index % len(geometry_modes)]
        if sample_index % 3 == 0:
            T_vals = np.ones((3,), dtype=float)
            time_case = "uniform"
        else:
            T_vals = rng.uniform(0.4, 1.8, size=3)
            time_case = "nonuniform"

        if mode == "near_linear":
            base = rng.normal()
            step = rng.uniform(0.3, 1.1)
            q_vals = np.asarray([base, base + step, base + 2 * step + 0.02 * rng.normal(), base + 3 * step], dtype=float)
        elif mode == "monotone":
            increments = np.abs(rng.normal(loc=0.8, scale=0.3, size=3))
            q_vals = np.concatenate(([0.0], np.cumsum(increments)))
        elif mode == "sharp_turn":
            q_vals = np.asarray([0.0, 1.2, -1.0, 0.6], dtype=float) + 0.1 * rng.normal(size=4)
        else:
            q_vals = np.asarray([-0.5, 0.7, 1.6, -0.2], dtype=float) + 0.2 * rng.normal(size=4)

        s2_cmp = compare_s2_exact_vs_catmull(q_vals[1], q_vals[2], q_vals[3], T_vals[1], T_vals[2])
        s4_cmp = compare_s4_exact_vs_catmull(*q_vals, *T_vals)
        time_imbalance = _time_imbalance_metric(*T_vals)

        rows.append(
            {
                "sample_index": sample_index,
                "time_case": time_case,
                "geometry_mode": mode,
                "time_imbalance": time_imbalance,
                "s2_abs_error": s2_cmp["abs_error"],
                "s4_state_error_norm": s4_cmp["state_error_norm"],
                "s4_coeff_error_norm": s4_cmp["coeff_error_norm"],
                "s4_objective_gap": s4_cmp["objective_gap"],
            }
        )

        s2_exact_values.append(s2_cmp["v_exact"])
        s2_catmull_values.append(s2_cmp["v_catmull"])
        s4_exact_vel.extend(s4_cmp["v_exact"].tolist())
        s4_catmull_vel.extend(s4_cmp["v_catmull"].tolist())
        coeff_errors.append(s4_cmp["coeff_error_norm"])
        objective_gaps.append(s4_cmp["objective_gap"])
        imbalance.append(time_imbalance)

    return {
        "rows": rows,
        "s2_exact": np.asarray(s2_exact_values, dtype=float),
        "s2_catmull": np.asarray(s2_catmull_values, dtype=float),
        "s4_exact_vel": np.asarray(s4_exact_vel, dtype=float),
        "s4_catmull_vel": np.asarray(s4_catmull_vel, dtype=float),
        "coeff_errors": np.asarray(coeff_errors, dtype=float),
        "objective_gaps": np.asarray(objective_gaps, dtype=float),
        "time_imbalance": np.asarray(imbalance, dtype=float),
    }


def plot_comparison_figures(benchmark: dict[str, Any], results_dir: str | Path) -> dict[str, str]:
    results_dir = ensure_dir(results_dir)
    fig_paths = {}
    save_scatter_with_diagonal(
        benchmark["s2_exact"],
        benchmark["s2_catmull"],
        results_dir / "s2_exact_vs_catmull_scatter.png",
        title="s=2 Exact vs Catmull Velocity",
        xlabel="v_exact",
        ylabel="v_catmull",
    )
    fig_paths["s2_exact_vs_catmull_scatter"] = str(results_dir / "s2_exact_vs_catmull_scatter.png")

    save_scatter_with_diagonal(
        benchmark["s4_exact_vel"],
        benchmark["s4_catmull_vel"],
        results_dir / "s4_velocity_exact_vs_catmull_scatter.png",
        title="s=4 Exact vs Catmull-like Velocity Components",
        xlabel="v_exact component",
        ylabel="v_catmull component",
    )
    fig_paths["s4_velocity_exact_vs_catmull_scatter"] = str(results_dir / "s4_velocity_exact_vs_catmull_scatter.png")

    save_histogram(
        benchmark["coeff_errors"],
        results_dir / "s4_coeff_error_hist.png",
        title="s=4 Center-Coefficient Error Histogram",
        xlabel="||c_exact - c_heuristic||",
    )
    fig_paths["s4_coeff_error_hist"] = str(results_dir / "s4_coeff_error_hist.png")

    save_histogram(
        benchmark["objective_gaps"],
        results_dir / "s4_objective_gap_hist.png",
        title="s=4 Objective Gap Histogram",
        xlabel="J_heuristic - J_exact",
    )
    fig_paths["s4_objective_gap_hist"] = str(results_dir / "s4_objective_gap_hist.png")

    representative = compare_s4_exact_vs_catmull(0.2, 1.1, -0.6, 0.8, 0.6, 1.4, 0.5)
    grid = np.linspace(0.0, 1.4, num=250)
    exact_curve = np.asarray([np.polyval(list(reversed(representative["c_exact"])), t) if False else 0.0 for t in grid], dtype=float)
    heuristic_curve = np.asarray([np.polyval(list(reversed(representative["c_heuristic"])), t) if False else 0.0 for t in grid], dtype=float)
    exact_curve = np.asarray([sum(coef * (t ** power) for power, coef in enumerate(representative["c_exact"])) for t in grid], dtype=float)
    heuristic_curve = np.asarray([sum(coef * (t ** power) for power, coef in enumerate(representative["c_heuristic"])) for t in grid], dtype=float)
    save_overlay_curves(
        grid,
        [("exact analytic", exact_curve), ("catmull heuristic", heuristic_curve)],
        results_dir / "representative_curve_overlay.png",
        title="Representative Center Segment Overlay",
        xlabel="Local Time",
        ylabel="p(t)",
    )
    fig_paths["representative_curve_overlay"] = str(results_dir / "representative_curve_overlay.png")

    order = np.argsort(benchmark["time_imbalance"])
    save_overlay_curves(
        benchmark["time_imbalance"][order],
        [("coeff error", benchmark["coeff_errors"][order]), ("objective gap", benchmark["objective_gaps"][order])],
        results_dir / "error_vs_time_imbalance.png",
        title="s=4 Error vs Time Imbalance",
        xlabel="Time Imbalance max(T)/min(T)",
        ylabel="Error / Gap",
    )
    fig_paths["error_vs_time_imbalance"] = str(results_dir / "error_vs_time_imbalance.png")
    return fig_paths


def main(results_dir: str | Path = DEFAULT_RESULTS_DIR, seed: int = 42) -> dict[str, Any]:
    results_dir = ensure_dir(results_dir)
    benchmark = run_random_benchmark(seed=seed)
    fig_paths = plot_comparison_figures(benchmark, results_dir)

    s2_metrics = {
        "max_abs_error": float(np.max(np.abs(benchmark["s2_exact"] - benchmark["s2_catmull"]))),
        "mean_abs_error": float(np.mean(np.abs(benchmark["s2_exact"] - benchmark["s2_catmull"]))),
    }
    s4_metrics = {
        "mean_coeff_error": float(np.mean(benchmark["coeff_errors"])),
        "max_coeff_error": float(np.max(benchmark["coeff_errors"])),
        "mean_objective_gap": float(np.mean(benchmark["objective_gaps"])),
        "min_objective_gap": float(np.min(benchmark["objective_gaps"])),
        "max_objective_gap": float(np.max(benchmark["objective_gaps"])),
    }
    save_json(results_dir / "s2_compare_metrics.json", s2_metrics)
    save_json(results_dir / "s4_compare_metrics.json", s4_metrics)
    save_csv(
        results_dir / "random_benchmark_table.csv",
        [
            "sample_index",
            "time_case",
            "geometry_mode",
            "time_imbalance",
            "s2_abs_error",
            "s4_state_error_norm",
            "s4_coeff_error_norm",
            "s4_objective_gap",
        ],
        benchmark["rows"],
    )
    return {
        "inputs": {"seed": seed, "num_samples": len(benchmark["rows"])},
        "solution": {},
        "metrics": {"s2": s2_metrics, "s4": s4_metrics},
        "fig_paths": fig_paths,
    }


if __name__ == "__main__":
    summary = main()
    print("Phase 4 Catmull comparison")
    print(f"s2 max abs error: {summary['metrics']['s2']['max_abs_error']:.3e}")
    print(f"s4 mean coeff error: {summary['metrics']['s4']['mean_coeff_error']:.3e}")
    print(f"s4 mean objective gap: {summary['metrics']['s4']['mean_objective_gap']:.3e}")
    print(f"results dir: {DEFAULT_RESULTS_DIR}")

