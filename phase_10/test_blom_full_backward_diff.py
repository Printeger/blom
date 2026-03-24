from __future__ import annotations

import unittest

import numpy as np

from phase_10.blom_full_backward_diff import (
    default_weights,
    evaluate_full_objective,
    finite_difference_gradient,
    full_backward_diff_dense,
    full_backward_diff_reparam,
    full_backward_diff_sparse,
    q_tau_to_xi,
    qbar_to_full,
    representative_case,
    tau_to_T,
    xi_to_qbar_tau,
)


class Phase10FullBackwardDiffTests(unittest.TestCase):
    def test_dense_and_sparse_match(self) -> None:
        q, T = representative_case()
        weights = default_weights()
        dense = full_backward_diff_dense(q, T, weights, k=2)
        sparse = full_backward_diff_sparse(q, T, weights, k=2)
        self.assertLess(np.max(np.abs(dense["grad_theta"] - sparse["grad_theta"])), 1e-8)

    def test_reparam_matches_fd(self) -> None:
        q, T = representative_case()
        tau = np.log(np.expm1(T - 0.2))
        weights = default_weights()
        payload = full_backward_diff_reparam(q, tau, weights, T_min=0.2, k=2)
        xi = q_tau_to_xi(q, tau)

        def objective_from_xi(xi_vec: np.ndarray) -> float:
            q_bar, tau_vec = xi_to_qbar_tau(xi_vec)
            q_full = qbar_to_full(q_bar, float(q[0]), float(q[-1]))
            T_full = tau_to_T(tau_vec, T_min=0.2)
            return float(evaluate_full_objective(q_full, T_full, weights, k=2)["value"])

        fd = finite_difference_gradient(objective_from_xi, xi, eps=1e-6)
        self.assertLess(np.max(np.abs(payload["grad_xi"] - fd)), 1e-4)


if __name__ == "__main__":
    unittest.main()
