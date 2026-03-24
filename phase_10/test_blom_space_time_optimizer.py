from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from phase_10.blom_full_backward_diff import T_to_tau, default_weights, q_tau_to_xi, representative_case
from phase_10.blom_space_time_optimizer import optimizer_step, run_space_time_optimization


class Phase10SpaceTimeOptimizerTests(unittest.TestCase):
    def test_optimizer_step_preserves_positive_durations(self) -> None:
        q, T = representative_case()
        tau = T_to_tau(T, T_min=0.2)
        step = optimizer_step(
            xi=q_tau_to_xi(q, tau),
            weights=default_weights(),
            T_min=0.2,
            step_rule="armijo",
            step_size=5e-3,
            config={"q_start": float(q[0]), "q_end": float(q[-1])},
        )
        self.assertTrue((step["T_new"] > 0.2).all())

    def test_run_space_time_optimization_smoke(self) -> None:
        q, T = representative_case()
        tau = T_to_tau(T, T_min=0.2)
        with tempfile.TemporaryDirectory() as tmpdir:
            result = run_space_time_optimization(
                q,
                tau,
                default_weights(),
                n_steps=3,
                step_rule="armijo",
                step_size=5e-3,
                T_min=0.2,
                save_dir=tmpdir,
            )
            self.assertLessEqual(result["summary"]["objective_final"], result["summary"]["objective_initial"] + 1e-8)
            self.assertTrue((Path(tmpdir) / "optimizer" / "phase10_opt_summary.json").exists())


if __name__ == "__main__":
    unittest.main()
