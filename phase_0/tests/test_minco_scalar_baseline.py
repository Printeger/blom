from __future__ import annotations

import unittest

from phase_0.blom_problem_setup import BLOMProblemSetup
from phase_0.minco_scalar_baseline import solve_minco_scalar
from phase_0.phase0_checks import check_minco_result
from phase_0.trajectory_eval import eval_piecewise, sample_trajectory


class MincoScalarBaselineTests(unittest.TestCase):
    def test_minco_solution_passes_phase0_checks(self) -> None:
        setup = BLOMProblemSetup.make_random(M=6, seed=42)
        result = solve_minco_scalar(setup)
        checks = check_minco_result(result)
        self.assertTrue(checks["passed"])
        self.assertGreaterEqual(checks["cost"], 0.0)

    def test_trajectory_eval_matches_boundary_waypoints(self) -> None:
        setup = BLOMProblemSetup.make_random(M=4, seed=42)
        result = solve_minco_scalar(setup)
        self.assertAlmostEqual(eval_piecewise(result.coeffs, setup.T, 0.0), setup.q[0])
        self.assertAlmostEqual(
            eval_piecewise(result.coeffs, setup.T, setup.total_time),
            setup.q[-1],
        )

    def test_sampling_returns_requested_grid_size(self) -> None:
        setup = BLOMProblemSetup.make_random(M=3, seed=42)
        result = solve_minco_scalar(setup)
        times, values = sample_trajectory(result.coeffs, setup.T, num=25, order=2)
        self.assertEqual(times.shape, (25,))
        self.assertEqual(values.shape, (25,))


if __name__ == "__main__":
    unittest.main()
