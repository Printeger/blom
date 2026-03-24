"""
Numeric analytic-core validation for the canonical BLOM case (s=4, k=2).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from phase_1.minco_scalar_baseline import evaluate_segment
from phase_3.blom_strict_local_qp import build_local_problem, extract_segment_coeff, solve_blom_strict_local_qp
from phase_4.utils.hermite_utils import (
    build_C4 as _build_C4,
    build_G4 as _build_G4,
    build_H_mid,
    build_R4 as _build_R4,
    build_gamma_mid,
    build_g_mid,
    build_m_minus as _build_m_minus,
    build_m_plus as _build_m_plus,
    evaluate_physical_polynomial,
    hermite_reconstruct_center_segment as _hermite_reconstruct_center_segment,
    outer_constant_term,
    outer_linear_term,
    outer_rank_one_hessian,
    quadratic_objective_value,
)
from phase_4.utils.io_utils import ensure_dir, save_json, save_npy
from phase_4.utils.plotting_utils import save_bar_chart, save_contour, save_heatmap, save_overlay_curves


DEFAULT_RESULTS_DIR = Path("phase_4/results/s4_numeric")


def build_m_plus(h: float) -> np.ndarray:
    return _build_m_plus(h)


def build_m_minus(h: float) -> np.ndarray:
    return _build_m_minus(h)


def left_outer_cost_matrix(qL: float, qR: float, h: float) -> dict[str, Any]:
    return {
        "hessian": outer_rank_one_hessian("left", h),
        "linear": outer_linear_term("left", qL, qR, h),
        "constant": outer_constant_term(qL, qR, h),
    }


def right_outer_cost_matrix(qL: float, qR: float, h: float) -> dict[str, Any]:
    return {
        "hessian": outer_rank_one_hessian("right", h),
        "linear": outer_linear_term("right", qL, qR, h),
        "constant": outer_constant_term(qL, qR, h),
    }


def build_C4() -> np.ndarray:
    return _build_C4()


def build_G4() -> np.ndarray:
    return _build_G4()


def build_R4() -> np.ndarray:
    return _build_R4()


def middle_segment_cost_matrix(qL: float, qR: float, h: float) -> dict[str, Any]:
    H_mid = build_H_mid(h)
    g_mid = build_g_mid(h)
    return {
        "H_mid": H_mid,
        "g_mid": g_mid,
        "constant": build_gamma_mid(qL, qR, h),
    }


def build_local_quadratic_s4_k2(
    q_im2: float,
    q_im1: float,
    q_i: float,
    q_ip1: float,
    T_im1: float,
    T_i: float,
    T_ip1: float,
) -> dict[str, Any]:
    h_minus = float(T_im1)
    h_mid = float(T_i)
    h_plus = float(T_ip1)

    left = left_outer_cost_matrix(q_im2, q_im1, h_minus)
    right = right_outer_cost_matrix(q_i, q_ip1, h_plus)
    middle = middle_segment_cost_matrix(q_im1, q_i, h_mid)
    H_mid = middle["H_mid"]
    g_mid = middle["g_mid"]
    q_mid = np.asarray([q_im1, q_i], dtype=float)
    R4 = build_R4()
    C_mid = float(middle["constant"])

    A2 = H_mid.copy()
    A2[:3, :3] += left["hessian"]
    A2[3:, 3:] += right["hessian"]
    A2 = 0.5 * (A2 + A2.T)

    B2 = g_mid @ q_mid
    B2[:3] += left["linear"]
    B2[3:] += right["linear"]

    C2 = left["constant"] + C_mid + right["constant"]
    return {
        "A2": A2,
        "B2": B2,
        "C2": float(C2),
        "H_mid": H_mid,
        "g_mid": g_mid,
        "left": left,
        "right": right,
        "R4": R4,
        "C4": build_C4(),
        "G4": build_G4(),
        "inputs": {
            "q_local": [q_im2, q_im1, q_i, q_ip1],
            "T_local": [T_im1, T_i, T_ip1],
        },
    }


def local_objective_s4_k2(x_vec: np.ndarray, quadratic: dict[str, Any]) -> float:
    return quadratic_objective_value(quadratic["A2"], quadratic["B2"], quadratic["C2"], x_vec)


def solve_local_system_s4_k2(
    q_im2: float,
    q_im1: float,
    q_i: float,
    q_ip1: float,
    T_im1: float,
    T_i: float,
    T_ip1: float,
) -> dict[str, Any]:
    quadratic = build_local_quadratic_s4_k2(q_im2, q_im1, q_i, q_ip1, T_im1, T_i, T_ip1)
    A2 = quadratic["A2"]
    B2 = quadratic["B2"]
    x_opt = np.linalg.solve(A2, B2)
    residual = A2 @ x_opt - B2
    objective = local_objective_s4_k2(x_opt, quadratic)
    return {
        **quadratic,
        "x_opt": x_opt,
        "residual_norm": float(np.linalg.norm(residual)),
        "symmetry_error": float(np.max(np.abs(A2 - A2.T))),
        "eigvals": np.linalg.eigvalsh(A2),
        "objective": objective,
    }


def hermite_reconstruct_center_segment(x_local: np.ndarray, q_im1: float, q_i: float, T_i: float) -> dict[str, np.ndarray]:
    x_local = np.asarray(x_local, dtype=float).reshape(6)
    return _hermite_reconstruct_center_segment(q_im1, q_i, x_local[:3], x_local[3:], T_i)


def _synthetic_global_for_phase3(
    q_im2: float,
    q_im1: float,
    q_i: float,
    q_ip1: float,
    T_im1: float,
    T_i: float,
    T_ip1: float,
) -> tuple[np.ndarray, np.ndarray]:
    q_global = np.asarray([q_im2, q_im2, q_im1, q_i, q_ip1, q_ip1], dtype=float)
    T_global = np.asarray([1.0, T_im1, T_i, T_ip1, 1.0], dtype=float)
    return q_global, T_global


def _extract_local_state_from_center_coeff(coeff: np.ndarray, T_i: float) -> np.ndarray:
    return np.asarray(
        [
            evaluate_segment(coeff, 0.0, order=1),
            evaluate_segment(coeff, 0.0, order=2),
            evaluate_segment(coeff, 0.0, order=3),
            evaluate_segment(coeff, float(T_i), order=1),
            evaluate_segment(coeff, float(T_i), order=2),
            evaluate_segment(coeff, float(T_i), order=3),
        ],
        dtype=float,
    )


def compare_with_blom_strict_qp(
    q_im2: float,
    q_im1: float,
    q_i: float,
    q_ip1: float,
    T_im1: float,
    T_i: float,
    T_ip1: float,
) -> dict[str, Any]:
    analytic = solve_local_system_s4_k2(q_im2, q_im1, q_i, q_ip1, T_im1, T_i, T_ip1)
    reconstructed = hermite_reconstruct_center_segment(analytic["x_opt"], q_im1, q_i, T_i)

    q_global, T_global = _synthetic_global_for_phase3(q_im2, q_im1, q_i, q_ip1, T_im1, T_i, T_ip1)
    problem = build_local_problem(q_global, T_global, i=3, k=2)
    qp_solution = solve_blom_strict_local_qp(problem, method="kkt")
    center_coeff_qp = extract_segment_coeff(qp_solution, seg_idx=3)
    x_qp = _extract_local_state_from_center_coeff(center_coeff_qp, T_i)

    objective_qp_from_x = local_objective_s4_k2(x_qp, analytic)
    center_coeff_exact = reconstructed["coeff"]
    coeff_error = float(np.linalg.norm(center_coeff_exact - center_coeff_qp))
    state_error = float(np.linalg.norm(analytic["x_opt"] - x_qp))

    endpoint_residuals = {
        "p0": abs(evaluate_physical_polynomial(center_coeff_exact, 0.0, order=0) - q_im1),
        "p1": abs(evaluate_physical_polynomial(center_coeff_exact, T_i, order=0) - q_i),
        "v0": abs(evaluate_physical_polynomial(center_coeff_exact, 0.0, order=1) - analytic["x_opt"][0]),
        "a0": abs(evaluate_physical_polynomial(center_coeff_exact, 0.0, order=2) - analytic["x_opt"][1]),
        "j0": abs(evaluate_physical_polynomial(center_coeff_exact, 0.0, order=3) - analytic["x_opt"][2]),
        "v1": abs(evaluate_physical_polynomial(center_coeff_exact, T_i, order=1) - analytic["x_opt"][3]),
        "a1": abs(evaluate_physical_polynomial(center_coeff_exact, T_i, order=2) - analytic["x_opt"][4]),
        "j1": abs(evaluate_physical_polynomial(center_coeff_exact, T_i, order=3) - analytic["x_opt"][5]),
    }

    return {
        "analytic": analytic,
        "x_exact": analytic["x_opt"],
        "x_qp": x_qp,
        "c_exact": center_coeff_exact,
        "c_qp": center_coeff_qp,
        "state_error_norm": state_error,
        "coeff_error_norm": coeff_error,
        "analytic_objective": analytic["objective"],
        "qp_objective": float(qp_solution["objective"]),
        "qp_objective_from_x": objective_qp_from_x,
        "endpoint_residuals": endpoint_residuals,
        "qp_summary": qp_solution["summary"],
    }


def _sample_center_curve(coeff: np.ndarray, duration: float, num: int = 200) -> tuple[np.ndarray, np.ndarray]:
    grid = np.linspace(0.0, float(duration), num=num)
    values = np.asarray([evaluate_physical_polynomial(coeff, t, order=0) for t in grid], dtype=float)
    return grid, values


def _energy_landscape(quadratic: dict[str, Any], x_opt: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    idx_a, idx_b = 0, 3
    span_a = np.linspace(x_opt[idx_a] - 2.0, x_opt[idx_a] + 2.0, 120)
    span_b = np.linspace(x_opt[idx_b] - 2.0, x_opt[idx_b] + 2.0, 120)
    X, Y = np.meshgrid(span_a, span_b)
    Z = np.zeros_like(X)
    for row in range(X.shape[0]):
        for col in range(X.shape[1]):
            trial = x_opt.copy()
            trial[idx_a] = X[row, col]
            trial[idx_b] = Y[row, col]
            Z[row, col] = local_objective_s4_k2(trial, quadratic)
    return X, Y, Z


def main(results_dir: str | Path = DEFAULT_RESULTS_DIR, seed: int = 42) -> dict[str, Any]:
    results_dir = ensure_dir(results_dir)
    q_local = np.asarray([0.1, 0.85, -0.2, 0.65], dtype=float)
    T_local = np.asarray([0.9, 1.25, 0.7], dtype=float)
    comparison = compare_with_blom_strict_qp(*q_local, *T_local)
    analytic = comparison["analytic"]

    save_npy(results_dir / "A2_matrix.npy", analytic["A2"])
    save_npy(results_dir / "B2_vector.npy", analytic["B2"])
    save_npy(results_dir / "x_local_opt.npy", comparison["x_exact"])
    save_npy(results_dir / "center_coeffs.npy", comparison["c_exact"])
    save_npy(results_dir / "R4_matrix.npy", analytic["R4"])
    save_npy(results_dir / "C4_matrix.npy", analytic["C4"])
    save_npy(results_dir / "G4_matrix.npy", analytic["G4"])

    metrics = {
        "symmetry_error": analytic["symmetry_error"],
        "min_eigenvalue": float(np.min(analytic["eigvals"])),
        "max_eigenvalue": float(np.max(analytic["eigvals"])),
        "residual_norm": analytic["residual_norm"],
        "state_error_norm": comparison["state_error_norm"],
        "coeff_error_norm": comparison["coeff_error_norm"],
        "analytic_objective": comparison["analytic_objective"],
        "qp_objective": comparison["qp_objective"],
        "qp_objective_from_x": comparison["qp_objective_from_x"],
        "endpoint_residuals": comparison["endpoint_residuals"],
    }
    save_json(results_dir / "analytic_vs_qp_metrics.json", metrics)

    save_heatmap(
        analytic["A2"],
        results_dir / "A2_heatmap.png",
        title="A2(T) Heatmap for Representative Sample",
        xlabel="Column",
        ylabel="Row",
    )

    rng = np.random.default_rng(seed)
    eig_samples = []
    for _ in range(64):
        q_rand = np.concatenate(([0.0], rng.normal(scale=0.8, size=3)))
        T_rand = rng.uniform(0.5, 1.6, size=3)
        solved = solve_local_system_s4_k2(*q_rand, *T_rand)
        eig_samples.append(np.sort(solved["eigvals"]))
    eig_array = np.asarray(eig_samples, dtype=float)
    save_overlay_curves(
        np.arange(1, eig_array.shape[0] + 1),
        [(f"lambda_{j + 1}", eig_array[:, j]) for j in range(eig_array.shape[1])],
        results_dir / "A2_eigenvalues.png",
        title="A2(T) Eigenvalues Across Random Samples",
        xlabel="Sample Index",
        ylabel="Eigenvalue",
    )

    save_bar_chart(
        ["||x_exact-x_qp||", "||c_exact-c_qp||"],
        np.asarray([comparison["state_error_norm"], comparison["coeff_error_norm"]], dtype=float),
        results_dir / "analytic_vs_qp_error_bar.png",
        title="Representative Analytic-vs-QP Errors",
        xlabel="Error Type",
        ylabel="Norm",
    )

    sigma, exact_curve = _sample_center_curve(comparison["c_exact"], T_local[1])
    _, qp_curve = _sample_center_curve(comparison["c_qp"], T_local[1])
    save_overlay_curves(
        sigma,
        [("analytic", exact_curve), ("qp", qp_curve)],
        results_dir / "center_segment_overlay.png",
        title="Center Segment: Analytic Reconstruction vs BLOM-Strict QP",
        xlabel="Local Time",
        ylabel="p(t)",
    )

    X, Y, Z = _energy_landscape(analytic, comparison["x_exact"])
    save_contour(
        X,
        Y,
        Z,
        optimum=(comparison["x_exact"][0], comparison["x_exact"][3]),
        save_path=results_dir / "local_energy_landscape_2d.png",
        title="2D Slice of Local Analytic Energy Landscape",
        xlabel="v_{i-1}",
        ylabel="v_i",
    )

    return {
        "inputs": {
            "q_local": q_local.tolist(),
            "T_local": T_local.tolist(),
        },
        "matrices": {
            "A2_shape": list(analytic["A2"].shape),
            "B2_shape": list(analytic["B2"].shape),
            "C4_shape": list(analytic["C4"].shape),
            "G4_shape": list(analytic["G4"].shape),
            "R4_shape": list(analytic["R4"].shape),
        },
        "solution": {
            "x_exact": comparison["x_exact"].tolist(),
            "x_qp": comparison["x_qp"].tolist(),
            "c_exact": comparison["c_exact"].tolist(),
            "c_qp": comparison["c_qp"].tolist(),
        },
        "metrics": metrics,
        "fig_paths": {
            "A2_heatmap": str(results_dir / "A2_heatmap.png"),
            "A2_eigenvalues": str(results_dir / "A2_eigenvalues.png"),
            "analytic_vs_qp_error_bar": str(results_dir / "analytic_vs_qp_error_bar.png"),
            "center_segment_overlay": str(results_dir / "center_segment_overlay.png"),
            "local_energy_landscape_2d": str(results_dir / "local_energy_landscape_2d.png"),
        },
    }


if __name__ == "__main__":
    summary = main()
    print("Phase 4 s=4 numeric validation")
    print(f"symmetry error: {summary['metrics']['symmetry_error']:.3e}")
    print(f"min eigenvalue: {summary['metrics']['min_eigenvalue']:.3e}")
    print(f"state error vs QP: {summary['metrics']['state_error_norm']:.3e}")
    print(f"coeff error vs QP: {summary['metrics']['coeff_error_norm']:.3e}")
    print(f"results dir: {DEFAULT_RESULTS_DIR}")
