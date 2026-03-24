from __future__ import annotations

import unittest

import numpy as np

from phase_0.poly_basis import coeffs_from_endpoint_derivatives
from phase_1.minco_scalar_baseline import (
    PHASE1_CONTINUITY_ORDER,
    boundary_jet_errors,
    build_rhs,
    build_system_matrix,
    continuity_jumps,
    evaluate_trajectory,
    interpolation_errors,
    sample_trajectory,
    solve_minco_coefficients,
    system_residual,
)


class Phase1MincoScalarBaselineTests(unittest.TestCase):
    def test_single_segment_matches_closed_form_hermite_solution(self) -> None:
        q = np.array([0.5, -0.2])
        T = np.array([1.7])
        zeta_start = np.array([0.5, 0.3, -0.4, 0.2])
        zeta_end = np.array([-0.2, 0.1, 0.2, -0.1])

        result = solve_minco_coefficients(q, T, zeta_start, zeta_end, return_system=True)
        expected = coeffs_from_endpoint_derivatives(
            np.concatenate((zeta_start, zeta_end)),
            float(T[0]),
        )

        np.testing.assert_allclose(result["coeffs"][0], expected, atol=1e-11, rtol=1e-11)
        self.assertLess(system_residual(result["M"], result["c_vec"], result["b"]), 1e-11)

    def test_random_case_satisfies_interpolation_continuity_and_boundary_jets(self) -> None:
        rng = np.random.default_rng(42)
        M_seg = 5
        q = np.concatenate(([0.0], rng.normal(size=M_seg)))
        T = rng.uniform(0.5, 1.5, size=M_seg)
        zeta_start = np.array([q[0], 0.0, 0.0, 0.0])
        zeta_end = np.array([q[-1], 0.0, 0.0, 0.0])

        result = solve_minco_coefficients(q, T, zeta_start, zeta_end, return_system=True)
        interp = interpolation_errors(result["coeffs"], T, q)
        jumps = continuity_jumps(result["coeffs"], T, max_order=PHASE1_CONTINUITY_ORDER)
        bc = boundary_jet_errors(result["coeffs"], T, zeta_start, zeta_end)

        self.assertLess(np.max(np.abs(interp)), 1e-10)
        self.assertLess(np.max(np.abs(jumps)), 1e-10)
        self.assertLess(np.max(np.abs(bc)), 1e-10)
        self.assertLess(system_residual(result["M"], result["c_vec"], result["b"]), 1e-10)

    def test_solver_is_deterministic_for_identical_inputs(self) -> None:
        q = np.array([0.0, 0.7, -0.3, 1.2])
        T = np.array([0.8, 1.1, 0.9])
        zeta_start = np.array([q[0], 0.0, 0.0, 0.0])
        zeta_end = np.array([q[-1], 0.0, 0.0, 0.0])

        first = solve_minco_coefficients(q, T, zeta_start, zeta_end)
        second = solve_minco_coefficients(q, T, zeta_start, zeta_end)

        np.testing.assert_allclose(first["coeffs"], second["coeffs"])
        np.testing.assert_allclose(first["c_vec"], second["c_vec"])
        self.assertAlmostEqual(first["residual_norm"], second["residual_norm"])

    def test_two_segment_symmetric_case_is_time_reversible(self) -> None:
        q = np.array([0.0, 1.0, 0.0])
        T = np.array([1.0, 1.0])
        zeta_start = np.array([0.0, 0.0, 0.0, 0.0])
        zeta_end = np.array([0.0, 0.0, 0.0, 0.0])

        result = solve_minco_coefficients(q, T, zeta_start, zeta_end)
        total_time = float(np.sum(T))
        times = np.linspace(0.0, total_time, 41)
        forward = np.asarray(
            [evaluate_trajectory(result["coeffs"], T, t, order=0) for t in times],
            dtype=float,
        )
        backward = np.asarray(
            [evaluate_trajectory(result["coeffs"], T, total_time - t, order=0) for t in times],
            dtype=float,
        )

        np.testing.assert_allclose(forward, backward, atol=1e-10, rtol=1e-10)

    def test_single_segment_time_scaling_preserves_normalized_position_profile(self) -> None:
        q = np.array([0.0, 1.0])
        zeta_start = np.array([0.0, 0.0, 0.0, 0.0])
        zeta_end = np.array([1.0, 0.0, 0.0, 0.0])

        short = solve_minco_coefficients(q, np.array([1.0]), zeta_start, zeta_end)
        long = solve_minco_coefficients(q, np.array([2.0]), zeta_start, zeta_end)

        taus = np.linspace(0.0, 1.0, 9)
        short_samples = sample_trajectory(short["coeffs"], np.array([1.0]), num_per_segment=9, orders=(0, 4))
        long_samples = sample_trajectory(long["coeffs"], np.array([2.0]), num_per_segment=9, orders=(0, 4))

        np.testing.assert_allclose(short_samples["t_global"], taus, atol=1e-12)
        np.testing.assert_allclose(long_samples["t_global"], 2.0 * taus, atol=1e-12)
        np.testing.assert_allclose(short_samples["order_0"], long_samples["order_0"], atol=1e-10)
        np.testing.assert_allclose(short_samples["order_4"], 16.0 * long_samples["order_4"], atol=1e-8)

    def test_random_stability_over_multiple_cases(self) -> None:
        rng = np.random.default_rng(123)
        for M_seg in (1, 2, 5):
            for _ in range(4):
                q = np.concatenate(([0.0], rng.normal(size=M_seg)))
                T = rng.uniform(0.2, 1.7, size=M_seg)
                zeta_start = np.array([q[0], 0.0, 0.0, 0.0])
                zeta_end = np.array([q[-1], 0.0, 0.0, 0.0])
                result = solve_minco_coefficients(q, T, zeta_start, zeta_end, return_system=True)

                self.assertTrue(np.all(np.isfinite(result["coeffs"])))
                self.assertTrue(np.all(np.isfinite(result["c_vec"])))
                self.assertLess(system_residual(result["M"], result["c_vec"], result["b"]), 1e-9)

    def test_build_system_matrix_and_rhs_have_expected_shapes(self) -> None:
        q = np.array([0.0, 0.5, -0.1, 0.3])
        T = np.array([0.7, 1.2, 0.9])
        zeta_start = np.array([0.0, 0.0, 0.0, 0.0])
        zeta_end = np.array([0.3, 0.0, 0.0, 0.0])

        M_mat = build_system_matrix(T)
        b_vec = build_rhs(q, zeta_start, zeta_end)

        self.assertEqual(M_mat.shape, (24, 24))
        self.assertEqual(b_vec.shape, (24,))


if __name__ == "__main__":
    unittest.main()
