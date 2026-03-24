from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from phase_7.blom_phase7_extra_experiments import (
    assemble_scheme_C_light,
    representative_case,
    run_interior_vs_full,
    run_large_M_sweep,
    run_phase7_extra_experiments,
    run_schemeC_light_assembly,
    run_time_regime_split,
)


class Phase7ExtraExperimentTests(unittest.TestCase):
    def setUp(self) -> None:
        self.q, self.T = representative_case()

    def test_large_M_sweep_runs_and_saves_artifacts(self) -> None:
        with TemporaryDirectory() as tmp_dir:
            result = run_large_M_sweep(M_values=[8, 12], k_max=8, seed=7, save_dir=tmp_dir)
            self.assertEqual(sorted(result["records"]), [8, 12])
            self.assertTrue(Path(tmp_dir, "summary_large_M.csv").exists())
            self.assertTrue(Path(tmp_dir, "coef_error_vs_k_by_M.png").exists())

    def test_time_regime_split_runs_and_saves_summary(self) -> None:
        with TemporaryDirectory() as tmp_dir:
            result = run_time_regime_split(M=8, n_trials=2, k_max=8, seed=7, save_dir=tmp_dir)
            self.assertGreater(len(result["records"]), 0)
            self.assertTrue(Path(tmp_dir, "summary_time_regime.csv").exists())
            self.assertTrue(Path(tmp_dir, "uniform_vs_bounded_slope_boxplot.png").exists())

    def test_interior_errors_are_available(self) -> None:
        with TemporaryDirectory() as tmp_dir:
            result = run_interior_vs_full(q=self.q, T=self.T, n_trials=2, save_dir=tmp_dir)
            final_k = result["case"]["k_values"][-1]
            final_c = result["case"]["actual"]["C"][final_k]
            self.assertIn("interior_coef_error", final_c)
            self.assertTrue(Path(tmp_dir, "summary_interior_vs_full.csv").exists())

    def test_scheme_c_light_preserves_c0_and_saves_artifacts(self) -> None:
        corrected = assemble_scheme_C_light(self.q, self.T, k=2, mode="C1")
        self.assertLess(corrected["stats"][0]["max_abs"], 1e-8)
        with TemporaryDirectory() as tmp_dir:
            result = run_schemeC_light_assembly(q=self.q, T=self.T, save_dir=tmp_dir)
            self.assertIn("C1", result["curves"])
            self.assertTrue(Path(tmp_dir, "schemeC_raw_vs_corrected_error.png").exists())

    def test_run_all_extra_experiments_generates_overview(self) -> None:
        with TemporaryDirectory() as tmp_dir:
            result = run_phase7_extra_experiments(q=self.q, T=self.T, M_values=[8, 12], n_trials=2, save_dir=tmp_dir, seed=7)
            self.assertGreater(len(result["overview_rows"]), 0)
            self.assertTrue(Path(tmp_dir, "compare_summary", "final_experiment_overview.csv").exists())


if __name__ == "__main__":
    unittest.main()

