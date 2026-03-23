from __future__ import annotations

import unittest

import numpy as np

from phase_0.blom_local_qp import solve_blom_local_qp
from phase_0.blom_problem_setup import BLOMProblemSetup
from phase_0.phase0_checks import check_blom_local_result


class BlomLocalQPTests(unittest.TestCase):
    def test_local_qp_solution_passes_phase0_checks(self) -> None:
        setup = BLOMProblemSetup.make_random(M=6, seed=42)
        result = solve_blom_local_qp(setup, center_segment=3)
        checks = check_blom_local_result(result)
        self.assertTrue(checks["passed"])
        self.assertEqual(result.active_segments, (2, 3, 4))

    def test_local_qp_is_deterministic(self) -> None:
        setup = BLOMProblemSetup.make_random(M=6, seed=42)
        first = solve_blom_local_qp(setup, center_segment=2)
        second = solve_blom_local_qp(setup, center_segment=2)
        np.testing.assert_allclose(first.window_coeffs, second.window_coeffs)
        self.assertAlmostEqual(first.cost, second.cost)

    def test_edge_window_uses_two_segments(self) -> None:
        setup = BLOMProblemSetup.make_random(M=6, seed=42)
        result = solve_blom_local_qp(setup, center_segment=0)
        self.assertEqual(result.active_segments, (0, 1))
        self.assertEqual(result.window_coeffs.shape[0], 2)


if __name__ == "__main__":
    unittest.main()
