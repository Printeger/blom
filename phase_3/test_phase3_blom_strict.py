"""
Tests for Phase 3 BLOM-Strict local validation tools.
"""

from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np

from phase_1.minco_scalar_baseline import solve_minco_coefficients
from phase_3.blom_strict_feasible_init import build_feasible_local_spline
from phase_3.blom_strict_local_kkt import build_local_constraints
from phase_3.blom_strict_local_qp import (
    build_local_problem,
    build_window,
    solve_blom_strict_local_qp,
)
from phase_3.validate_phase3_blom_strict import run_phase3_validation


class Phase3ValidationTest(unittest.TestCase):
    def setUp(self) -> None:
        self.q = np.array([0.0, 0.25, -0.4, 0.7, -0.3, 0.1], dtype=float)
        self.T = np.array([0.9, 1.1, 0.8, 1.0, 1.2], dtype=float)
        self.zeta_start = np.array([self.q[0], 0.0, 0.0, 0.0], dtype=float)
        self.zeta_end = np.array([self.q[-1], 0.0, 0.0, 0.0], dtype=float)

    def test_build_window_canonical_interior(self) -> None:
        window = build_window(i=3, k=2, M=5)
        self.assertEqual(window["segments"], [2, 3, 4])
        self.assertEqual(window["L"], 2)
        self.assertEqual(window["R"], 4)
        self.assertEqual(window["m"], 3)
        self.assertEqual(window["window_type"], "interior")

    def test_feasible_init_satisfies_constraints(self) -> None:
        feasible = build_feasible_local_spline(
            self.q,
            self.T,
            i=3,
            k=2,
            zeta_start=self.zeta_start,
            zeta_end=self.zeta_end,
        )
        G, d_vec = build_local_constraints(
            self.q,
            self.T,
            i=3,
            k=2,
            zeta_start=self.zeta_start,
            zeta_end=self.zeta_end,
        )
        residual = np.linalg.norm(G @ feasible["c_loc"] - d_vec)
        self.assertLess(residual, 1e-10)

    def test_kkt_and_reduced_agree_and_natural_bc_hold(self) -> None:
        problem = build_local_problem(
            self.q,
            self.T,
            i=3,
            k=2,
            zeta_start=self.zeta_start,
            zeta_end=self.zeta_end,
        )
        solution_kkt = solve_blom_strict_local_qp(problem, method="kkt")
        solution_red = solve_blom_strict_local_qp(problem, method="reduced")
        self.assertLess(np.linalg.norm(solution_kkt["c_vec"] - solution_red["c_vec"]), 1e-10)
        self.assertLess(solution_kkt["summary"]["max_natural_bc_residual"], 1e-9)
        self.assertGreater(solution_red["reduced_hessian_min_eig"], 0.0)

    def test_full_window_matches_phase1_global(self) -> None:
        local_problem = build_local_problem(
            self.q,
            self.T,
            i=3,
            k=20,
            zeta_start=self.zeta_start,
            zeta_end=self.zeta_end,
        )
        local_solution = solve_blom_strict_local_qp(local_problem, method="reduced")
        global_solution = solve_minco_coefficients(
            self.q,
            self.T,
            self.zeta_start,
            self.zeta_end,
        )
        self.assertLess(np.linalg.norm(local_solution["c_vec"] - global_solution["c_vec"]), 1e-9)

    def test_validation_runner_generates_required_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            summary = run_phase3_validation(
                M_values=(5,),
                num_restarts=3,
                epsilons=(1e-2, 1e-4),
                results_dir=Path(tempdir),
            )
            self.assertEqual(len(summary["cases"]), 3)
            self.assertTrue((Path(tempdir) / "tables" / "table_feasibility_summary.csv").exists())
            self.assertTrue((Path(tempdir) / "tables" / "table_uniqueness_summary.csv").exists())
            self.assertTrue((Path(tempdir) / "tables" / "table_perturbation_continuity.csv").exists())
            self.assertTrue((Path(tempdir) / "logs" / "phase3_validation_summary.json").exists())


if __name__ == "__main__":
    unittest.main()

