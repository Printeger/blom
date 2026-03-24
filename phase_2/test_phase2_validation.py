from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from phase_2.phase2_validation import (
    build_waypoint_selector,
    compare_exact_vs_fd_jacobian,
    compute_exact_jacobian_q,
    compute_segmentwise_influence_norms,
    compute_waypoint_influence_profile,
    ensure_results_dirs,
    estimate_effective_bandwidth,
    finite_difference_jacobian_q,
    make_nonuniform_time_case,
    run_phase2_validation_suite,
    run_scaling_experiment,
)
from phase_1.minco_scalar_baseline import build_system_matrix


class Phase2ValidationTests(unittest.TestCase):
    def test_waypoint_selector_has_expected_shape_and_pattern(self) -> None:
        S_q = build_waypoint_selector(4)
        self.assertEqual(S_q.shape, (32, 3))
        expected_rows = [(1, 2), (3, 4), (5, 6)]
        for column, (r0, r1) in enumerate(expected_rows):
            self.assertEqual(S_q[r0, column], 1.0)
            self.assertEqual(S_q[r1, column], 1.0)
            self.assertEqual(int(np.count_nonzero(S_q[:, column])), 2)

    def test_exact_and_fd_jacobian_match(self) -> None:
        case = make_nonuniform_time_case(5, seed=42, case_name="fd_match")
        A_mat = build_system_matrix(case.T)
        S_q = build_waypoint_selector(case.T.size)
        J_exact = compute_exact_jacobian_q(A_mat, S_q)
        J_fd = finite_difference_jacobian_q(
            case.q,
            case.T,
            case.zeta_start,
            case.zeta_end,
            eps=1e-7,
        )
        comparison = compare_exact_vs_fd_jacobian(J_exact, J_fd)
        self.assertLess(comparison["relative_error"], 1e-5)
        self.assertLess(comparison["max_abs_error"], 1e-5)

    def test_influence_and_bandwidth_helpers_behave_consistently(self) -> None:
        J_q = np.zeros((24, 2), dtype=float)
        J_q[0:8, 0] = 2.0
        J_q[8:16, 0] = 0.5
        J_q[16:24, 1] = 3.0

        influence = compute_segmentwise_influence_norms(J_q, 3)
        self.assertEqual(influence.shape, (3, 2))

        profile = compute_waypoint_influence_profile(np.zeros((24, 2)), 3, 1)
        self.assertEqual(profile["segment_indices"].shape, (3,))

        synthetic = np.zeros((24, 2))
        synthetic[0:8, 0] = 1.0
        synthetic[8:16, 0] = 1e-3
        synthetic[16:24, 1] = 2.0
        bandwidth = estimate_effective_bandwidth(synthetic, 3, tol=1e-4)
        self.assertGreaterEqual(bandwidth["max_effective_bandwidth"], 1)

    def test_validation_suite_generates_required_artifacts(self) -> None:
        case = make_nonuniform_time_case(4, seed=123, case_name="smoke")
        with tempfile.TemporaryDirectory() as tmpdir:
            result = run_phase2_validation_suite(
                case.q,
                case.T,
                case.zeta_start,
                case.zeta_end,
                case_name=case.case_name,
                results_dir=tmpdir,
            )
            figures = Path(tmpdir) / "figures"
            tables = Path(tmpdir) / "tables"
            logs = Path(tmpdir) / "logs"

            self.assertTrue((figures / "A_sparsity_M4_casesmoke.png").exists())
            self.assertTrue((figures / "Jq_heatmap_M4_casesmoke.png").exists())
            self.assertTrue((figures / "block_influence_M4_casesmoke.png").exists())
            self.assertTrue((figures / "jacobian_fd_compare_M4_casesmoke.png").exists())
            self.assertTrue((tables / "jacobian_error_summary.csv").exists())
            self.assertTrue((tables / "effective_bandwidth_summary.csv").exists())
            self.assertTrue((tables / "far_field_influence_summary.csv").exists())
            self.assertTrue(any(logs.iterdir()))
            self.assertIn("J_q", result)

    def test_scaling_experiment_generates_summary_plot(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            summary = run_scaling_experiment(
                M_values=(4, 8),
                case_name="uniform_time",
                results_dir=tmpdir,
            )
            self.assertEqual(len(summary["records"]), 2)
            self.assertTrue(
                (Path(tmpdir) / "figures" / "effective_bandwidth_vs_M_caseuniform_time.png").exists()
            )


if __name__ == "__main__":
    unittest.main()
