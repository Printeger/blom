from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def demo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def results_cache_dir() -> Path:
    path = demo_root() / "results_cache"
    path.mkdir(parents=True, exist_ok=True)
    return path


def default_case(M: int = 6, regime: str = "bounded-nonuniform", seed: int = 42) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    q = np.zeros((M + 1,), dtype=float)
    q[1:] = np.cumsum(rng.normal(loc=0.0, scale=0.85, size=M))
    if regime == "uniform":
        T = np.full((M,), 1.0, dtype=float)
    else:
        T = rng.uniform(0.6, 1.6, size=M).astype(float)
    return q, T


def representative_case() -> tuple[np.ndarray, np.ndarray]:
    from phase_10.blom_full_backward_diff import representative_case as phase10_case

    return phase10_case()


def case_for_demo(M: int, regime: str, seed: int = 42) -> tuple[np.ndarray, np.ndarray]:
    if M == 6 and regime == "bounded-nonuniform":
        return representative_case()
    return default_case(M=M, regime=regime, seed=seed)


def case_signature(q: np.ndarray, T: np.ndarray, **kwargs: Any) -> str:
    q = np.asarray(q, dtype=float).reshape(-1)
    T = np.asarray(T, dtype=float).reshape(-1)
    head = {
        "q": np.round(q, 6).tolist(),
        "T": np.round(T, 6).tolist(),
        **kwargs,
    }
    return json.dumps(head, sort_keys=True)


def save_cache_json(name: str, payload: dict[str, Any]) -> Path:
    path = results_cache_dir() / f"{name}.json"
    path.write_text(json.dumps(_serialize(payload), indent=2, sort_keys=True), encoding="utf-8")
    return path


def load_json(path: str | Path) -> dict[str, Any] | None:
    path = Path(path)
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def load_csv_frame(path: str | Path) -> pd.DataFrame | None:
    path = Path(path)
    if not path.exists():
        return None
    return pd.read_csv(path)


def phase_result_path(*parts: str) -> Path:
    return repo_root().joinpath(*parts)


def load_phase10_summary() -> dict[str, Any] | None:
    return load_json(phase_result_path("phase_10", "results", "phase10_framework", "summary", "phase10_suite_summary.json"))


def load_phase10_benchmark_frame() -> pd.DataFrame | None:
    return load_csv_frame(phase_result_path("phase_10", "results", "phase10_framework", "summary", "phase10_overview.csv"))


def load_phase10_ablation_frame() -> pd.DataFrame | None:
    return load_csv_frame(phase_result_path("phase_10", "results", "phase10_framework", "summary", "phase10_ablation_summary.csv"))


def load_phase8_suite() -> dict[str, Any] | None:
    return load_json(phase_result_path("phase_8", "results", "phase8_validation", "compare", "phase8_suite_summary.json"))


def _serialize(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, dict):
        return {str(k): _serialize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_serialize(v) for v in value]
    return value


def comparison_table_frame() -> pd.DataFrame:
    rows = [
        {
            "Object": "B-spline",
            "Primary variables": "control points",
            "Interpolation": "not required",
            "Local support": "exact local basis support",
            "Jacobian pattern": "banded in basis space",
            "Variational meaning": "not intrinsic here",
            "Main role": "local geometric-control baseline",
        },
        {
            "Object": "MINCO",
            "Primary variables": "interpolation waypoints",
            "Interpolation": "yes",
            "Local support": "global coupling",
            "Jacobian pattern": "typically dense / wide",
            "Variational meaning": "global minimum-snap reference",
            "Main role": "global optimal reference",
        },
        {
            "Object": "raw Scheme C",
            "Primary variables": "interpolation waypoints",
            "Interpolation": "segment-wise local extraction",
            "Local support": "exact local support in canonical analysis",
            "Jacobian pattern": "block-banded",
            "Variational meaning": "local BLOM object",
            "Main role": "primary BLOM object",
        },
        {
            "Object": "Scheme A / B",
            "Primary variables": "shared or consensus junction states",
            "Interpolation": "assembled",
            "Local support": "weaker than raw Scheme C",
            "Jacobian pattern": "broader than raw Scheme C",
            "Variational meaning": "continuity-vs-locality trade-off",
            "Main role": "comparison baselines",
        },
    ]
    return pd.DataFrame(rows)
