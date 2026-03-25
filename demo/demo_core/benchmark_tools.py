from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from demo.demo_core.data_interface import load_json, phase_result_path


def load_phase10_benchmark_artifacts() -> dict[str, Any]:
    summary_dir = phase_result_path("phase_10", "results", "phase10_framework", "summary")
    benchmark_json = load_json(summary_dir / "phase10_benchmark_summary.json")
    benchmark_rows = None
    if benchmark_json is not None and "benchmark_rows" in benchmark_json:
        benchmark_rows = pd.DataFrame(benchmark_json["benchmark_rows"])
    return {
        "suite_overview": pd.read_csv(summary_dir / "phase10_overview.csv") if (summary_dir / "phase10_overview.csv").exists() else None,
        "benchmark_rows": benchmark_rows,
        "ablation": pd.read_csv(summary_dir / "phase10_ablation_summary.csv") if (summary_dir / "phase10_ablation_summary.csv").exists() else None,
        "suite": load_json(summary_dir / "phase10_suite_summary.json"),
        "benchmark": benchmark_json,
    }


def load_phase8_artifacts() -> dict[str, Any]:
    compare_dir = phase_result_path("phase_8", "results", "phase8_validation", "compare")
    return {
        "suite": load_json(compare_dir / "phase8_suite_summary.json"),
        "summary": (compare_dir / "phase8_interpretation_summary.md").read_text(encoding="utf-8") if (compare_dir / "phase8_interpretation_summary.md").exists() else None,
    }
