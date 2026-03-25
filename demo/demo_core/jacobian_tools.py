from __future__ import annotations

from typing import Any

import numpy as np


def threshold_mask(J: np.ndarray, threshold: float) -> np.ndarray:
    J = np.asarray(J, dtype=float)
    return np.abs(J) >= float(threshold)


def block_row_norms(J: np.ndarray, block_rows: int) -> np.ndarray:
    J = np.asarray(J, dtype=float)
    if J.shape[0] % block_rows != 0:
        raise ValueError("Row count must be divisible by block_rows.")
    blocks = J.reshape(J.shape[0] // block_rows, block_rows, J.shape[1])
    return np.linalg.norm(blocks, axis=1)


def support_width_from_blocks(block_norms: np.ndarray, tol: float = 1e-10) -> dict[str, Any]:
    block_norms = np.asarray(block_norms, dtype=float)
    widths = []
    active_spans = []
    for col in range(block_norms.shape[1]):
        active = np.where(block_norms[:, col] > tol)[0]
        if active.size == 0:
            widths.append(0)
            active_spans.append((None, None))
        else:
            widths.append(int(active[-1] - active[0] + 1))
            active_spans.append((int(active[0]), int(active[-1])))
    return {
        "widths": widths,
        "mean_width": float(np.mean(widths)) if widths else 0.0,
        "max_width": float(np.max(widths)) if widths else 0.0,
        "active_spans": active_spans,
    }


def outside_band_energy(J: np.ndarray, expected_band: int, block_rows: int = 8) -> float:
    J = np.asarray(J, dtype=float)
    norms = block_row_norms(J, block_rows=block_rows)
    energy_total = float(np.sum(norms**2))
    if energy_total <= 1e-15:
        return 0.0
    energy_out = 0.0
    for seg_idx in range(norms.shape[0]):
        for col_idx in range(norms.shape[1]):
            if abs(seg_idx - col_idx) > expected_band:
                energy_out += float(norms[seg_idx, col_idx] ** 2)
    return energy_out / energy_total


def summarize_dense_jacobian(J: np.ndarray, block_rows: int = 8, tol: float = 1e-10, expected_band: int | None = None) -> dict[str, Any]:
    norms = block_row_norms(J, block_rows=block_rows)
    support = support_width_from_blocks(norms, tol=tol)
    summary = {
        "mean_width": support["mean_width"],
        "max_width": support["max_width"],
        "nnz": int(np.count_nonzero(np.abs(J) > tol)),
    }
    if expected_band is not None:
        summary["outside_band_energy"] = outside_band_energy(J, expected_band=expected_band, block_rows=block_rows)
    return summary
