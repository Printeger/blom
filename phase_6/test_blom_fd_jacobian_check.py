from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from phase_6.blom_fd_jacobian_check import (
    compute_jacobian_errors,
    raw_local_jacobians,
    representative_case,
    run_compare_all_schemes,
    run_random_trials,
    theoretical_mask_c_T,
    theoretical_mask_c_q,
)


class Phase6JacobianCheckTests(unittest.TestCase):
    def setUp(self) -> None:
        self.q, self.T = representative_case()
        self.segment_idx = 3

    def test_theoretical_masks_have_expected_support(self) -> None:
        mask_q = theoretical_mask_c_q(self.T.size, self.segment_idx)
        mask_T = theoretical_mask_c_T(self.T.size, self.segment_idx)
        self.assertEqual(mask_q.shape, (8, self.T.size - 1))
        self.assertEqual(mask_T.shape, (8, self.T.size))
        self.assertTrue(np.all(mask_q[:, 0:4]))
        self.assertFalse(np.any(mask_q[:, 4:]))
        self.assertTrue(np.all(mask_T[:, 1:4]))
        self.assertFalse(np.any(mask_T[:, :1]))
        self.assertFalse(np.any(mask_T[:, 4:]))

    def test_raw_analytic_and_fd_jacobians_match_for_q(self) -> None:
        analytic = raw_local_jacobians(self.q, self.T, self.segment_idx, mode="analytic")
        fd = raw_local_jacobians(self.q, self.T, self.segment_idx, mode="fd", eps=1e-6)
        errors = compute_jacobian_errors(analytic["J_c_q"], fd["J_c_q"])
        self.assertLess(errors["max_abs_error"], 1e-5)
        self.assertLess(errors["relative_error"], 1e-5)

    def test_raw_analytic_and_fd_jacobians_match_for_T(self) -> None:
        analytic = raw_local_jacobians(self.q, self.T, self.segment_idx, mode="analytic")
        fd = raw_local_jacobians(self.q, self.T, self.segment_idx, mode="fd", eps=1e-6)
        errors = compute_jacobian_errors(analytic["J_c_T"], fd["J_c_T"])
        self.assertLess(errors["max_abs_error"], 1e-4)
        self.assertLess(errors["relative_error"], 1e-4)

    def test_scheme_c_stays_narrower_than_A_and_B(self) -> None:
        with TemporaryDirectory() as tmp_dir:
            result = run_compare_all_schemes(self.q, self.T, save_dir=tmp_dir)
            bw_c = result["scheme_C"]["q_bandwidth"]["max_effective_bandwidth"]
            bw_a = result["scheme_A"]["q_bandwidth"]["max_effective_bandwidth"]
            bw_b = result["scheme_B"]["q_bandwidth"]["max_effective_bandwidth"]
            self.assertLessEqual(bw_c, 2)
            self.assertGreaterEqual(bw_a, bw_c)
            self.assertGreaterEqual(bw_b, bw_c)

    def test_artifacts_and_random_trials_are_generated(self) -> None:
        with TemporaryDirectory() as tmp_dir:
            compare = run_compare_all_schemes(self.q, self.T, save_dir=tmp_dir)
            random_result = run_random_trials(n_trials=4, M=6, save_dir=tmp_dir, seed=7)
            self.assertTrue(Path(tmp_dir, "compare", "jacobian_bandwidth_compare.png").exists())
            self.assertTrue(Path(tmp_dir, "compare", "jacobian_error_stats.csv").exists())
            self.assertTrue(Path(tmp_dir, "raw_scheme_C", "jacobian_theory_vs_fd_q_raw.png").exists())
            self.assertIn("A", random_result["summary"]["schemes"])
            self.assertIn("raw", random_result["summary"])
            self.assertIn("scheme_A", compare)


if __name__ == "__main__":
    unittest.main()

