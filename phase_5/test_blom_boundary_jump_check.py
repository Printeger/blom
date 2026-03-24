from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from phase_5.blom_boundary_jump_check import (
    compute_jumps,
    run_boundary_jump_check,
    run_compare_all_schemes,
    run_random_trials,
    summarize_jumps,
)


class BoundaryJumpCheckTests(unittest.TestCase):
    def setUp(self) -> None:
        self.q = np.asarray([0.0, 0.9, -0.7, 1.3, 0.1, 0.8], dtype=float)
        self.T = np.asarray([0.7, 1.2, 0.8, 1.1, 0.9], dtype=float)

    def test_compute_jumps_reports_zero_for_identical_segments(self) -> None:
        coeff_left = np.asarray([0.0, 1.0, 0.0, 0.0], dtype=float)
        coeff_right = np.asarray([1.0, 1.0, 0.0, 0.0], dtype=float)
        coeffs = np.vstack((coeff_left, coeff_right))
        T = np.asarray([1.0, 1.0], dtype=float)
        jumps = compute_jumps(coeffs, T, s=2)
        summary = summarize_jumps(jumps)
        self.assertLessEqual(summary[0]["max_abs"], 1e-14)
        self.assertLessEqual(summary[1]["max_abs"], 1e-14)

    def test_scheme_a_reaches_lower_order_continuity(self) -> None:
        result = run_boundary_jump_check(self.q, self.T, scheme="A")
        for order in range(4):
            self.assertLessEqual(result["stats"][order]["max_abs"], 1e-9)
        self.assertGreater(result["system"]["condition_number"], 0.0)

    def test_scheme_b_reaches_lower_order_continuity_and_reduces_dispersion(self) -> None:
        result = run_boundary_jump_check(self.q, self.T, scheme="B")
        for order in range(4):
            self.assertLessEqual(result["stats"][order]["max_abs"], 1e-9)
        if result["consensus"]["pre_dispersion"].size:
            self.assertGreater(np.max(result["consensus"]["pre_dispersion"]), 0.0)
            self.assertLessEqual(np.max(result["consensus"]["post_lower_order_jump_norm"]), 1e-9)

    def test_scheme_c_is_exactly_c0_but_not_generally_c3(self) -> None:
        result = run_boundary_jump_check(self.q, self.T, scheme="C")
        self.assertLessEqual(result["stats"][0]["max_abs"], 1e-9)
        lower_nonzero = max(result["stats"][order]["max_abs"] for order in range(1, 4))
        self.assertGreater(lower_nonzero, 1e-6)

    def test_compare_and_random_trials_generate_outputs(self) -> None:
        with TemporaryDirectory() as tmp_dir:
            compare = run_compare_all_schemes(self.q, self.T, save_dir=tmp_dir)
            random_result = run_random_trials(n_trials=6, M=8, save_dir=tmp_dir, seed=7)
            self.assertTrue(Path(tmp_dir, "jump_lower_orders_compare.png").exists())
            self.assertTrue(Path(tmp_dir, "scheme_comparison_summary.csv").exists())
            self.assertTrue(Path(tmp_dir, "phase5_interpretation_summary.md").exists())
            self.assertTrue(Path(tmp_dir, "jump_boxplot_random_trials.png").exists())
            self.assertIn("A", compare["scheme_results"])
            self.assertEqual(set(random_result["summary"]), {"A", "B", "C"})


if __name__ == "__main__":
    unittest.main()
