"""
Batch validation runner for Phase 3 BLOM-Strict.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import json
from typing import Any

import numpy as np

from phase_1.minco_scalar_baseline import solve_minco_coefficients
from phase_3.blom_strict_local_qp import (
    build_local_problem,
    sample_local_trajectory,
    solve_blom_strict_local_qp,
)
from phase_3.phase3_plotting import (
    plot_continuity_jumps,
    plot_local_trajectory,
    plot_natural_bc_residuals,
    plot_perturbation_response,
    plot_uniqueness_multistart,
    plot_window_layout,
)


DEFAULT_SEED = 42
DEFAULT_S = 4
DEFAULT_K = 2
DEFAULT_RESULTS_DIR = Path("phase_3/results")


@dataclass
class Phase3Case:
    M: int
    q: np.ndarray
    T: np.ndarray
    zeta_start: np.ndarray
    zeta_end: np.ndarray


def ensure_results_dirs(base_dir: str | Path = DEFAULT_RESULTS_DIR) -> dict[str, Path]:
    base = Path(base_dir)
    figures = base / "figures"
    tables = base / "tables"
    logs = base / "logs"
    for directory in (base, figures, tables, logs):
        directory.mkdir(parents=True, exist_ok=True)
    return {"base": base, "figures": figures, "tables": tables, "logs": logs}


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)


def make_phase3_case(M: int, *, seed: int = DEFAULT_SEED) -> Phase3Case:
    rng = np.random.default_rng(seed + M)
    q = np.concatenate(([0.0], rng.normal(scale=0.7, size=M)))
    T = rng.uniform(0.65, 1.35, size=M)
    zeta_start = np.array([q[0], 0.0, 0.0, 0.0], dtype=float)
    zeta_end = np.array([q[-1], 0.0, 0.0, 0.0], dtype=float)
    return Phase3Case(M=M, q=q, T=T, zeta_start=zeta_start, zeta_end=zeta_end)


def _window_positions(M: int) -> list[tuple[str, int]]:
    if M < 3:
        return [("left", 1), ("right", M)]
    interior = max(2, min(M - 1, (M + 1) // 2))
    return [("left", 1), ("interior", interior), ("right", M)]


def _pairwise_distances(vectors: list[np.ndarray]) -> np.ndarray:
    count = len(vectors)
    matrix = np.zeros((count, count), dtype=float)
    for row in range(count):
        for col in range(count):
            matrix[row, col] = float(np.linalg.norm(vectors[row] - vectors[col]))
    return matrix


def _choose_perturb_target(problem: dict[str, Any]) -> int:
    window = problem["window"]
    candidate_knots = [knot for knot in window["knots"] if 0 < knot < problem["T"].size]
    if not candidate_knots:
        return 1
    center_knot = window["center"]
    return min(candidate_knots, key=lambda knot: abs(knot - center_knot))


def _trajectory_sup_norm(solution_a: dict[str, Any], solution_b: dict[str, Any]) -> float:
    samples_a = sample_local_trajectory(solution_a["coeffs"], solution_a["problem"]["window"], num_per_segment=120)
    samples_b = sample_local_trajectory(solution_b["coeffs"], solution_b["problem"]["window"], num_per_segment=120)
    return float(np.max(np.abs(samples_a["order_0"] - samples_b["order_0"])))


def _run_multistart(problem: dict[str, Any], num_restarts: int, seed: int) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    baseline = solve_blom_strict_local_qp(problem, method="reduced")
    null_dim = baseline["nullspace_basis"].shape[1]

    solutions = []
    for _ in range(num_restarts):
        z0 = rng.normal(size=null_dim) if null_dim > 0 else np.zeros((0,), dtype=float)
        solution = solve_blom_strict_local_qp(
            problem,
            method="reduced_iterative",
            z0=z0,
            max_iters=4000,
            tol=1e-13,
        )
        solutions.append(solution)

    coeff_vectors = [solution["c_vec"] for solution in solutions]
    objectives = np.asarray([solution["objective"] for solution in solutions], dtype=float)
    distances = _pairwise_distances(coeff_vectors)
    return {
        "solutions": solutions,
        "pairwise_coeff_diff": distances,
        "max_pairwise_coeff_diff": float(np.max(distances)),
        "objective_std": float(np.std(objectives)),
        "baseline": baseline,
    }


def _run_perturbation_continuity(
    case: Phase3Case,
    i: int,
    k: int,
    epsilons: tuple[float, ...],
) -> dict[str, Any]:
    base_problem = build_local_problem(case.q, case.T, i, k, zeta_start=case.zeta_start, zeta_end=case.zeta_end)
    base_solution = solve_blom_strict_local_qp(base_problem, method="kkt")
    target = _choose_perturb_target(base_problem)
    rows = []
    coeff_diffs = []
    cond_numbers = []

    for eps in epsilons:
        q_perturbed = case.q.copy()
        q_perturbed[target] += float(eps)
        perturbed_problem = build_local_problem(
            q_perturbed,
            case.T,
            i,
            k,
            zeta_start=case.zeta_start,
            zeta_end=case.zeta_end,
        )
        perturbed_solution = solve_blom_strict_local_qp(perturbed_problem, method="reduced")
        coeff_diff = float(np.linalg.norm(perturbed_solution["c_vec"] - base_solution["c_vec"]))
        traj_diff = _trajectory_sup_norm(base_solution, perturbed_solution)
        cond_number = float(perturbed_solution["reduced_hessian_cond"])
        coeff_diffs.append(coeff_diff)
        cond_numbers.append(cond_number)
        rows.append(
            {
                "M": case.M,
                "i": i,
                "k": k,
                "perturb_target": target,
                "eps": eps,
                "coeff_diff_norm": coeff_diff,
                "traj_diff_sup_norm": traj_diff,
                "cond_number": cond_number,
            }
        )

    return {
        "rows": rows,
        "target": target,
        "epsilons": np.asarray(epsilons, dtype=float),
        "coeff_diff_norms": np.asarray(coeff_diffs, dtype=float),
        "cond_numbers": np.asarray(cond_numbers, dtype=float),
    }


def _full_window_phase1_comparison(case: Phase3Case) -> dict[str, Any]:
    center = max(1, min(case.M, (case.M + 1) // 2))
    problem = build_local_problem(
        case.q,
        case.T,
        center,
        2 * case.M,
        zeta_start=case.zeta_start,
        zeta_end=case.zeta_end,
    )
    local_solution = solve_blom_strict_local_qp(problem, method="reduced")
    global_solution = solve_minco_coefficients(
        case.q,
        case.T,
        case.zeta_start,
        case.zeta_end,
        return_system=False,
    )
    diff = float(np.linalg.norm(local_solution["c_vec"] - global_solution["c_vec"]))
    return {
        "window_segments": problem["window"]["segments"],
        "coeff_diff_norm": diff,
        "local_objective": float(local_solution["objective"]),
    }


def run_phase3_validation(
    *,
    M_values: tuple[int, ...] = (5, 6),
    k: int = DEFAULT_K,
    num_restarts: int = 6,
    epsilons: tuple[float, ...] = (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6),
    seed: int = DEFAULT_SEED,
    results_dir: str | Path = DEFAULT_RESULTS_DIR,
) -> dict[str, Any]:
    directories = ensure_results_dirs(results_dir)
    figures_dir = directories["figures"]
    tables_dir = directories["tables"]
    logs_dir = directories["logs"]

    feasibility_rows: list[dict[str, Any]] = []
    uniqueness_rows: list[dict[str, Any]] = []
    perturbation_rows: list[dict[str, Any]] = []
    case_logs: list[dict[str, Any]] = []

    first_interior_perturbation: dict[str, Any] | None = None
    first_interior_problem: dict[str, Any] | None = None

    for M in M_values:
        case = make_phase3_case(M, seed=seed)
        for window_type, i in _window_positions(M):
            problem = build_local_problem(
                case.q,
                case.T,
                i,
                k,
                zeta_start=case.zeta_start,
                zeta_end=case.zeta_end,
            )
            solution_kkt = solve_blom_strict_local_qp(problem, method="kkt")
            solution_red = solve_blom_strict_local_qp(problem, method="reduced")
            summary = solution_kkt["summary"]
            solver_status = "ok" if np.linalg.norm(solution_kkt["c_vec"] - solution_red["c_vec"]) <= 1e-10 else "mismatch"

            multistart = _run_multistart(problem, num_restarts=num_restarts, seed=seed + 10 * M + i)
            if window_type == "interior" and first_interior_perturbation is None:
                first_interior_perturbation = _run_perturbation_continuity(case, i, k, epsilons)
                first_interior_problem = problem
                perturbation_rows.extend(first_interior_perturbation["rows"])

            sigma_knots = problem["window"]["sigma"]
            knot_values = np.asarray([case.q[index] for index in problem["window"]["knots"]], dtype=float)
            samples = sample_local_trajectory(solution_kkt["coeffs"], problem["window"], num_per_segment=120)

            suffix = f"M{M}_i{i}_k{k}"
            plot_window_layout(M, problem["window"], figures_dir / f"fig_window_layout_{suffix}.png")
            plot_local_trajectory(
                samples,
                sigma_knots,
                knot_values,
                figures_dir / f"fig_local_traj_{suffix}.png",
                title=f"Local Trajectory, M={M}, i={i}, k={k}",
            )
            plot_continuity_jumps(
                summary["continuity_jumps"],
                figures_dir / f"fig_continuity_jump_{suffix}.png",
                title=f"Continuity Jumps, M={M}, i={i}, k={k}",
            )
            plot_natural_bc_residuals(
                summary["natural_bc_residuals"],
                problem["s"],
                figures_dir / f"fig_natural_bc_residual_{suffix}.png",
                title=f"Natural BC Residuals, M={M}, i={i}, k={k}",
            )
            plot_uniqueness_multistart(
                multistart["pairwise_coeff_diff"],
                figures_dir / f"fig_uniqueness_multistart_{suffix}.png",
                title=f"Multistart Coefficient Distance, M={M}, i={i}, k={k}",
            )

            if window_type == "interior" and first_interior_perturbation is not None and first_interior_problem is problem:
                plot_perturbation_response(
                    first_interior_perturbation["epsilons"],
                    first_interior_perturbation["coeff_diff_norms"],
                    first_interior_perturbation["cond_numbers"],
                    figures_dir / f"fig_perturbation_response_{suffix}.png",
                    title=f"Perturbation Response, M={M}, i={i}, k={k}",
                )
            else:
                local_perturb = _run_perturbation_continuity(case, i, k, epsilons[:4])
                plot_perturbation_response(
                    local_perturb["epsilons"],
                    local_perturb["coeff_diff_norms"],
                    local_perturb["cond_numbers"],
                    figures_dir / f"fig_perturbation_response_{suffix}.png",
                    title=f"Perturbation Response, M={M}, i={i}, k={k}",
                )
                perturbation_rows.extend(local_perturb["rows"])

            feasibility_rows.append(
                {
                    "M": M,
                    "i": i,
                    "k": k,
                    "window_type": window_type,
                    "max_interp_error": summary["max_interp_error"],
                    "max_Cs_minus_1_jump": summary["max_Cs_minus_1_jump"],
                    "max_C2s_minus_2_jump": summary["max_C2s_minus_2_jump"],
                    "max_boundary_jet_error": summary["max_boundary_jet_error"],
                    "max_natural_bc_residual": summary["max_natural_bc_residual"],
                    "solver_status": solver_status,
                }
            )
            uniqueness_rows.append(
                {
                    "M": M,
                    "i": i,
                    "k": k,
                    "num_restarts": num_restarts,
                    "max_pairwise_coeff_diff": multistart["max_pairwise_coeff_diff"],
                    "objective_std": multistart["objective_std"],
                    "reduced_hessian_min_eig": solution_red["reduced_hessian_min_eig"],
                    "reduced_hessian_cond": solution_red["reduced_hessian_cond"],
                    "kkt_residual": solution_kkt["kkt_residual"],
                }
            )
            case_log = {
                "M": M,
                "i": i,
                "k": k,
                "window_type": window_type,
                "segments": problem["window"]["segments"],
                "kkt_objective": solution_kkt["objective"],
                "reduced_objective": solution_red["objective"],
                "coeff_diff_norm": float(np.linalg.norm(solution_kkt["c_vec"] - solution_red["c_vec"])),
                "summary": {
                    key: value
                    for key, value in summary.items()
                    if key not in {"continuity_jumps", "interpolation_errors", "boundary_jet_errors", "natural_bc_residuals"}
                },
                "reduced_hessian_min_eig": solution_red["reduced_hessian_min_eig"],
                "reduced_hessian_cond": solution_red["reduced_hessian_cond"],
                "multistart_max_pairwise_coeff_diff": multistart["max_pairwise_coeff_diff"],
                "multistart_objective_std": multistart["objective_std"],
            }
            case_logs.append(case_log)
            _write_json(logs_dir / f"phase3_case_{suffix}.json", case_log)

    full_window_checks = [
        _full_window_phase1_comparison(make_phase3_case(M, seed=seed)) for M in M_values
    ]

    _write_csv(
        tables_dir / "table_feasibility_summary.csv",
        [
            "M",
            "i",
            "k",
            "window_type",
            "max_interp_error",
            "max_Cs_minus_1_jump",
            "max_C2s_minus_2_jump",
            "max_boundary_jet_error",
            "max_natural_bc_residual",
            "solver_status",
        ],
        feasibility_rows,
    )
    _write_csv(
        tables_dir / "table_uniqueness_summary.csv",
        [
            "M",
            "i",
            "k",
            "num_restarts",
            "max_pairwise_coeff_diff",
            "objective_std",
            "reduced_hessian_min_eig",
            "reduced_hessian_cond",
            "kkt_residual",
        ],
        uniqueness_rows,
    )
    _write_csv(
        tables_dir / "table_perturbation_continuity.csv",
        [
            "M",
            "i",
            "k",
            "perturb_target",
            "eps",
            "coeff_diff_norm",
            "traj_diff_sup_norm",
            "cond_number",
        ],
        perturbation_rows,
    )

    summary = {
        "cases": case_logs,
        "full_window_phase1_checks": full_window_checks,
        "results_dir": {key: str(value) for key, value in directories.items()},
    }
    _write_json(logs_dir / "phase3_validation_summary.json", summary)
    return summary


def main() -> None:
    summary = run_phase3_validation()
    max_interp = max(case["summary"]["max_interp_error"] for case in summary["cases"])
    max_natural = max(case["summary"]["max_natural_bc_residual"] for case in summary["cases"])
    max_pairwise = max(case["multistart_max_pairwise_coeff_diff"] for case in summary["cases"])
    full_window_diff = max(check["coeff_diff_norm"] for check in summary["full_window_phase1_checks"])
    print("Phase 3 BLOM-Strict validation")
    print(f"cases: {len(summary['cases'])}")
    print(f"max interpolation error: {max_interp:.3e}")
    print(f"max natural BC residual: {max_natural:.3e}")
    print(f"max multistart coeff diff: {max_pairwise:.3e}")
    print(f"full-window vs Phase 1 coeff diff: {full_window_diff:.3e}")
    print(f"results dir: {summary['results_dir']['base']}")


if __name__ == "__main__":
    main()

