from __future__ import annotations

import unittest

from phase_0.blom_local_qp import solve_blom_local_qp
from phase_0.blom_problem_setup import BLOMProblemSetup
from phase_0.minco_scalar_baseline import solve_minco_scalar
from phase_0.phase0_checks import check_blom_local_result, check_minco_result, check_setup


class Phase0ChecksTests(unittest.TestCase):
    def test_check_setup_reports_basic_metrics(self) -> None:
        setup = BLOMProblemSetup.make_random(M=5, seed=42)
        report = check_setup(setup)
        self.assertTrue(report["passed"])
        self.assertEqual(report["theta_size"], 2 * setup.M - 1)

    def test_minco_check_reports_small_residuals(self) -> None:
        setup = BLOMProblemSetup.make_random(M=5, seed=42)
        result = solve_minco_scalar(setup)
        report = check_minco_result(result)
        self.assertLess(report["max_interpolation_error"], 1e-8)
        self.assertLess(report["max_junction_jump"], 1e-8)
        self.assertLess(report["max_boundary_residual"], 1e-8)

    def test_local_check_reports_small_residuals(self) -> None:
        setup = BLOMProblemSetup.make_random(M=5, seed=42)
        result = solve_blom_local_qp(setup, center_segment=2)
        report = check_blom_local_result(result)
        self.assertLess(report["max_interpolation_error"], 1e-8)
        self.assertLess(report["max_continuity_error"], 1e-8)
        self.assertLess(report["max_boundary_residual"], 1e-8)


if __name__ == "__main__":
    unittest.main()
