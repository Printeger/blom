from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from phase_8.blom_interior_matching_check import run_interior_matching_check
from phase_8.phase8_common import make_interior_sets, prepare_case


class InteriorMatchingCheckTests(unittest.TestCase):
    def test_make_interior_sets_default(self) -> None:
        sets = make_interior_sets(M=10, k=4, radius_mode="default")
        self.assertEqual(sets["r_k"], 4)
        self.assertEqual(sets["interior_idx"], [5, 6])
        self.assertEqual(sets["boundary_idx"], [1, 2, 3, 4, 7, 8, 9, 10])

    def test_run_interior_matching_smoke(self) -> None:
        q, T = prepare_case()
        with tempfile.TemporaryDirectory() as tmp_dir:
            payload = run_interior_matching_check(q, T, [2, 4, 6, 8], save_dir=tmp_dir)
            self.assertEqual(len(payload["records"]), 4)
            for record in payload["records"]:
                self.assertLessEqual(record["interior_l2"], record["full_l2"] + 1e-12)
            self.assertTrue(Path(tmp_dir, "interior_matching", "interior_matching_full_vs_interior.png").exists())


if __name__ == "__main__":
    unittest.main()

