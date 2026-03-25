from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from phase_10.blom_benchmark_suite import run_benchmark_suite
from phase_10.blom_full_backward_diff import default_weights, representative_case


class Phase10BenchmarkSuiteTests(unittest.TestCase):
    def test_benchmark_suite_smoke(self) -> None:
        q, T = representative_case()
        with tempfile.TemporaryDirectory() as tmpdir:
            result = run_benchmark_suite(
                q0=q,
                T0=T,
                weights=default_weights(),
                M_values=[6],
                k_values=[2],
                benchmark_methods=["raw_schemeC", "heuristic"],
                n_steps=2,
                T_min=0.2,
                save_dir=tmpdir,
                seed=7,
            )
            self.assertGreaterEqual(len(result["benchmark_rows"]), 2)
            self.assertTrue((Path(tmpdir) / "summary" / "phase10_benchmark_summary.json").exists())


if __name__ == "__main__":
    unittest.main()
