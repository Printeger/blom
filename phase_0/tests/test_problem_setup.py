from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from phase_0.blom_problem_setup import BLOMProblemSetup


class BLOMProblemSetupTests(unittest.TestCase):
    def test_make_random_is_deterministic(self) -> None:
        first = BLOMProblemSetup.make_random(M=4, seed=42)
        second = BLOMProblemSetup.make_random(M=4, seed=42)
        np.testing.assert_allclose(first.q, second.q)
        np.testing.assert_allclose(first.T, second.T)
        np.testing.assert_allclose(first.theta(), second.theta())

    def test_window_and_theta_shapes_follow_phase0_spec(self) -> None:
        setup = BLOMProblemSetup.make_random(M=5, seed=42)
        self.assertEqual(setup.window(0), (0, 1))
        self.assertEqual(setup.window(2), (1, 2, 3))
        self.assertEqual(setup.window(4), (3, 4))
        self.assertEqual(setup.theta().shape, (2 * setup.M - 1,))

    def test_json_round_trip(self) -> None:
        setup = BLOMProblemSetup.make_random(M=6, seed=42)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "setup.json"
            setup.save_json(path)
            loaded = BLOMProblemSetup.load_json(path)
        np.testing.assert_allclose(setup.q, loaded.q)
        np.testing.assert_allclose(setup.T, loaded.T)
        self.assertEqual(setup.window(3), loaded.window(3))

    def test_validate_rejects_invalid_duration(self) -> None:
        setup = BLOMProblemSetup(
            M=2,
            q=np.array([0.0, 1.0, 2.0]),
            T=np.array([1.0, 0.0]),
        )
        with self.assertRaises(ValueError):
            setup.validate()


if __name__ == "__main__":
    unittest.main()
