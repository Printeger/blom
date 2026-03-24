"""
Tests for Phase 4 analytic-validation scripts.
"""

from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np

from phase_4.blom_k2_s2_sympy import verify_catmull_equivalence_symbolic
from phase_4.blom_k2_s4_numeric import compare_with_blom_strict_qp, solve_local_system_s4_k2
from phase_4.blom_catmull_compare import compare_s4_exact_vs_catmull, main as compare_main


class Phase4AnalyticTest(unittest.TestCase):
    def test_s2_symbolic_catmull_equivalence_is_exact(self) -> None:
        verification = verify_catmull_equivalence_symbolic()
        self.assertEqual(str(verification["difference"]), "0")

    def test_s4_A2_is_symmetric_positive_definite(self) -> None:
        solved = solve_local_system_s4_k2(0.0, 0.8, -0.25, 0.5, 0.9, 1.2, 0.7)
        self.assertLess(np.max(np.abs(solved["A2"] - solved["A2"].T)), 1e-12)
        self.assertGreater(np.min(solved["eigvals"]), 0.0)
        self.assertLess(solved["residual_norm"], 1e-10)

    def test_s4_analytic_matches_phase3_qp(self) -> None:
        comparison = compare_with_blom_strict_qp(0.0, 0.8, -0.25, 0.5, 0.9, 1.2, 0.7)
        self.assertLess(comparison["state_error_norm"], 1e-8)
        self.assertLess(comparison["coeff_error_norm"], 1e-8)
        self.assertLess(max(comparison["endpoint_residuals"].values()), 1e-10)

    def test_s4_catmull_heuristic_is_not_exact(self) -> None:
        comparison = compare_s4_exact_vs_catmull(0.0, 0.8, -0.25, 0.5, 0.5, 1.4, 0.7)
        self.assertGreater(comparison["coeff_error_norm"], 1e-6)
        self.assertGreaterEqual(comparison["objective_gap"], -1e-10)

    def test_compare_script_generates_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            summary = compare_main(results_dir=Path(tempdir), seed=7)
            self.assertIn("s2", summary["metrics"])
            self.assertTrue((Path(tempdir) / "s2_compare_metrics.json").exists())
            self.assertTrue((Path(tempdir) / "s4_compare_metrics.json").exists())
            self.assertTrue((Path(tempdir) / "random_benchmark_table.csv").exists())


if __name__ == "__main__":
    unittest.main()

