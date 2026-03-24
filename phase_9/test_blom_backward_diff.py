from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from phase_9.blom_backward_diff import (
    backward_diff_banded,
    backward_diff_dense,
    compute_raw_schemeC_coeffs,
    compute_raw_schemeC_jacobians,
    control_cost,
    default_obs_config,
    default_weights,
    evaluate_minimal_objective,
    finite_difference_gradient,
    qT_to_theta,
    representative_case,
    soft_obstacle_penalty,
    theta_to_qT,
    time_penalty,
)
from phase_9.blom_phase9_validation_suite import run_phase9_validation_suite
from phase_9.blom_space_time_opt_demo import run_minimal_optimization_demo


class Phase9BackwardDiffTests(unittest.TestCase):
    def setUp(self) -> None:
        self.q, self.T = representative_case()
        self.weights = default_weights()
        self.obs = default_obs_config()

    def test_raw_schemeC_interface_shapes(self) -> None:
        coeffs = compute_raw_schemeC_coeffs(self.q, self.T)
        jac = compute_raw_schemeC_jacobians(self.q, self.T)
        M = self.T.size
        self.assertEqual(coeffs["coeff_blocks"].shape, (M, 8))
        self.assertEqual(coeffs["coeff_vec"].shape, (8 * M,))
        self.assertEqual(jac["J_c_q_dense"].shape, (8 * M, M - 1))
        self.assertEqual(jac["J_c_T_dense"].shape, (8 * M, M))

    def test_dense_equals_banded(self) -> None:
        dense = backward_diff_dense(self.q, self.T, self.weights, obs_config=self.obs)
        banded = backward_diff_banded(self.q, self.T, self.weights, obs_config=self.obs)
        np.testing.assert_allclose(dense["grad_theta"], banded["grad_theta"], atol=1e-10, rtol=1e-8)

    def test_control_term_gradient(self) -> None:
        coeffs = compute_raw_schemeC_coeffs(self.q, self.T)["coeff_blocks"]
        result = control_cost(coeffs, self.T)
        z = np.concatenate((coeffs.reshape(-1), self.T))

        def f(z_vec: np.ndarray) -> float:
            coeff = z_vec[: coeffs.size].reshape(coeffs.shape)
            durations = z_vec[coeffs.size :]
            return float(control_cost(coeff, durations)["value"])

        fd = finite_difference_gradient(f, z)
        analytic = np.concatenate((result["grad_c_blocks"].reshape(-1), result["grad_T"]))
        np.testing.assert_allclose(analytic, fd, atol=1e-5, rtol=1e-4)

    def test_time_term_gradient(self) -> None:
        result = time_penalty(self.T, weight=0.7)
        self.assertAlmostEqual(result["value"], float(0.7 * np.sum(self.T)))
        np.testing.assert_allclose(result["grad_T"], 0.7)

    def test_soft_obstacle_gradient(self) -> None:
        coeffs = compute_raw_schemeC_coeffs(self.q, self.T)["coeff_blocks"]
        result = soft_obstacle_penalty(coeffs, self.T, self.obs)
        z = np.concatenate((coeffs.reshape(-1), self.T))

        def f(z_vec: np.ndarray) -> float:
            coeff = z_vec[: coeffs.size].reshape(coeffs.shape)
            durations = z_vec[coeffs.size :]
            return float(soft_obstacle_penalty(coeff, durations, self.obs)["value"])

        fd = finite_difference_gradient(f, z)
        analytic = np.concatenate((result["grad_c_blocks"].reshape(-1), result["grad_T"]))
        np.testing.assert_allclose(analytic, fd, atol=1e-5, rtol=1e-4)

    def test_total_gradient_matches_fd(self) -> None:
        dense = backward_diff_dense(self.q, self.T, self.weights, obs_config=self.obs)
        q0 = float(self.q[0])
        qM = float(self.q[-1])
        theta = qT_to_theta(self.q, self.T)

        def objective(theta_vec: np.ndarray) -> float:
            q_full, T_full = theta_to_qT(theta_vec, q0, qM, self.T.size)
            return float(evaluate_minimal_objective(q_full, T_full, self.weights, obs_config=self.obs)["value"])

        fd = finite_difference_gradient(objective, theta, eps=1e-6)
        np.testing.assert_allclose(dense["grad_theta"], fd, atol=5e-4, rtol=5e-3)
        banded = backward_diff_banded(self.q, self.T, self.weights, obs_config=self.obs)
        np.testing.assert_allclose(banded["grad_theta"], fd, atol=5e-4, rtol=5e-3)

    def test_minimal_optimization_demo_decreases_objective(self) -> None:
        demo = run_minimal_optimization_demo(self.q, self.T, self.weights, obs_config=self.obs, n_steps=8, step_size=5e-3, T_min=0.2)
        self.assertGreaterEqual(demo["summary"]["objective_drop"], 0.0)

    def test_phase9_suite_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            payload = run_phase9_validation_suite(
                self.q,
                self.T,
                self.weights,
                obs_config=self.obs,
                n_steps=6,
                step_size=5e-3,
                T_min=0.2,
                save_dir=tmp_dir,
            )
            self.assertTrue(Path(tmp_dir, "compare", "phase9_overview.csv").exists())
            self.assertIn("gradcheck", payload)


if __name__ == "__main__":
    unittest.main()
