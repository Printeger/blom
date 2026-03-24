from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from phase_8.blom_phase8_supplementary_suite import run_phase8_supplementary_suite
from phase_8.exp4_raw_vs_reference_sanity import run_raw_vs_reference_sanity
from phase_8.exp5_empty_interior_risk import run_empty_interior_risk
from phase_8.phase8_supplementary_common import normalize_radius_mode


class Phase8SupplementarySuiteTests(unittest.TestCase):
    def test_normalize_radius_mode_aliases(self) -> None:
        self.assertEqual(normalize_radius_mode("half_k"), "half")
        self.assertEqual(normalize_radius_mode("k"), "default")
        self.assertEqual(normalize_radius_mode("three_half_k"), "one_and_half")

    def test_raw_vs_reference_sanity_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            payload = run_raw_vs_reference_sanity(
                M_values=[20],
                k_values=[2, 4],
                n_trials=1,
                delta_scales=[1e-6, 1e-4],
                save_dir=tmp_dir,
                seed=7,
            )
            self.assertGreater(len(payload["records"]), 0)
            self.assertGreater(len(payload["gamma_rows"]), 0)
            self.assertTrue(Path(tmp_dir, "exp4_raw_vs_reference_sanity", "phase8_raw_vs_reference_error.png").exists())

    def test_empty_interior_risk_marks_empty_cases(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            payload = run_empty_interior_risk(
                M_values=[10],
                k_values=[6, 8, 10, 12],
                radius_modes=["three_half_k"],
                save_dir=tmp_dir,
                seed=3,
            )
            empty_flags = [row["is_empty_interior"] for row in payload["records"]]
            self.assertTrue(any(empty_flags))
            self.assertTrue(Path(tmp_dir, "exp5_empty_interior_risk", "phase8_empty_interior_risk.png").exists())

    def test_full_suite_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            payload = run_phase8_supplementary_suite(
                M_values=[20, 40],
                k_values=[2, 4],
                radius_modes=["half_k", "k"],
                h_values=[1.0],
                nonuniform_boxes=[(0.7, 1.3)],
                n_trials=2,
                save_dir=tmp_dir,
                seed=5,
            )
            self.assertIn("exp6_two_bridge", payload)
            self.assertTrue(Path(tmp_dir, "summary", "phase8_supplementary_overview.csv").exists())
            self.assertTrue(Path(tmp_dir, "summary", "phase8_supplementary_summary.md").exists())


if __name__ == "__main__":
    unittest.main()
