from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from phase_8.blom_boundary_gap_decomposition import project_interior_boundary_errors, run_boundary_gap_decomposition
from phase_8.phase8_common import compute_boundary_response_matrix, prepare_case


class BoundaryGapDecompositionTests(unittest.TestCase):
    def test_projection_exact_energy_split(self) -> None:
        errors = np.asarray([1.0, 2.0, 3.0, 4.0], dtype=float)
        projection = project_interior_boundary_errors(errors, [2, 3], [1, 4])
        total = float(np.sum(errors**2))
        self.assertAlmostEqual(projection["interior_energy_ratio"] + projection["boundary_energy_ratio"], 1.0)
        self.assertAlmostEqual(projection["interior_l2"] ** 2 + projection["boundary_l2"] ** 2, total)

    def test_boundary_response_matrix_shape(self) -> None:
        q, T = prepare_case()
        response = compute_boundary_response_matrix(q, T, i=4, k=4)
        self.assertEqual(response["M_matrix"].shape, (8, 6))
        self.assertGreaterEqual(response["operator_norm"], 0.0)

    def test_run_boundary_gap_smoke(self) -> None:
        q, T = prepare_case()
        with tempfile.TemporaryDirectory() as tmp_dir:
            payload = run_boundary_gap_decomposition(q, T, [2, 4], save_dir=tmp_dir)
            self.assertEqual(len(payload["records"]), 2)
            self.assertTrue(Path(tmp_dir, "boundary_gap", "boundary_energy_ratio_vs_k.png").exists())


if __name__ == "__main__":
    unittest.main()

