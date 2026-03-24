from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from phase_8.blom_uniform_vs_nonuniform_interior import (
    run_uniform_vs_nonuniform_interior,
    sample_bounded_nonuniform_time,
    sample_uniform_time,
)
from phase_8.phase8_common import default_q_sampler


class UniformVsNonuniformInteriorTests(unittest.TestCase):
    def test_sample_uniform_time(self) -> None:
        T = sample_uniform_time(5, 1.25)
        self.assertTrue(np.allclose(T, 1.25))

    def test_sample_bounded_nonuniform_time(self) -> None:
        rng = np.random.default_rng(0)
        T = sample_bounded_nonuniform_time(6, 0.7, 1.3, rng)
        self.assertEqual(T.shape, (6,))
        self.assertTrue(np.all(T >= 0.7))
        self.assertTrue(np.all(T <= 1.3))

    def test_run_uniform_vs_nonuniform_smoke(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            payload = run_uniform_vs_nonuniform_interior(
                q_sampler=default_q_sampler,
                M=8,
                k_values=[2, 4, 6],
                h_values=[1.0],
                nonuniform_boxes=[(0.8, 1.2)],
                n_trials=2,
                save_dir=tmp_dir,
            )
            self.assertGreater(len(payload["aggregated"]), 0)
            self.assertTrue(Path(tmp_dir, "uniform_vs_nonuniform", "uniform_vs_nonuniform_interior_matching.png").exists())


if __name__ == "__main__":
    unittest.main()

