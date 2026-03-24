from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from phase_7.blom_convergence_vs_k import (
    compute_actual_blom_k,
    compute_ideal_truncated_blom_k,
    compute_minco_reference,
    make_k_grid,
    representative_case,
    run_convergence_vs_k,
    run_random_trials,
)


class Phase7ConvergenceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.q, self.T = representative_case()

    def test_make_k_grid_includes_terminal_value(self) -> None:
        self.assertEqual(make_k_grid(8), [2, 4, 6, 8])
        self.assertEqual(make_k_grid(7), [2, 4, 6, 7])

    def test_reference_and_full_ideal_truncation_match(self) -> None:
        reference = compute_minco_reference(self.q, self.T)
        full = compute_ideal_truncated_blom_k(self.q, self.T, k=self.T.size, reference=reference)
        self.assertLess(np.linalg.norm(full["c_vec"] - reference["c_vec"]), 1e-10)

    def test_actual_schemes_run_and_return_expected_shapes(self) -> None:
        for scheme in ("A", "B", "C"):
            result = compute_actual_blom_k(self.q, self.T, k=2, scheme=scheme)
            self.assertEqual(result["coeffs"].shape, (self.T.size, 8))
            self.assertTrue(np.isfinite(result["cost"]))

    def test_run_convergence_vs_k_generates_required_artifacts(self) -> None:
        with TemporaryDirectory() as tmp_dir:
            result = run_convergence_vs_k(self.q, self.T, k_values=[2, 4, self.T.size], save_dir=tmp_dir)
            self.assertIn("ideal", result["fits"])
            self.assertIn("A", result["fits"])
            self.assertTrue(Path(tmp_dir, "compare", "coef_error_vs_k_all.png").exists())
            self.assertTrue(Path(tmp_dir, "compare", "matching_error_vs_k.png").exists())
            self.assertTrue(Path(tmp_dir, "compare", "logfit_summary.csv").exists())
            self.assertTrue(Path(tmp_dir, "compare", "summary_phase7_convergence.json").exists())

    def test_random_trials_produce_summary_and_boxplot(self) -> None:
        with TemporaryDirectory() as tmp_dir:
            summary = run_random_trials(n_trials=3, M=8, k_values=[2, 4, 8], save_dir=tmp_dir, seed=7)
            self.assertIn("ideal", summary["summary"]["curves"])
            self.assertIn("B", summary["summary"]["curves"])
            self.assertTrue(Path(tmp_dir, "compare", "logfit_slope_boxplot_random_trials.png").exists())
            self.assertTrue(Path(tmp_dir, "compare", "random_trials_summary.csv").exists())


if __name__ == "__main__":
    unittest.main()

