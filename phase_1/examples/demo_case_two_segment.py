from __future__ import annotations

import numpy as np

from phase_1.minco_scalar_baseline import evaluate_trajectory, solve_minco_coefficients


def main() -> None:
    q = np.array([0.0, 1.0, 0.0])
    T = np.array([1.0, 1.0])
    zeta_start = np.array([0.0, 0.0, 0.0, 0.0])
    zeta_end = np.array([0.0, 0.0, 0.0, 0.0])

    result = solve_minco_coefficients(q, T, zeta_start, zeta_end)
    total_time = float(np.sum(T))
    times = np.linspace(0.0, total_time, 41)
    forward = np.asarray(
        [evaluate_trajectory(result["coeffs"], T, t, order=0) for t in times],
        dtype=float,
    )
    backward = np.asarray(
        [evaluate_trajectory(result["coeffs"], T, total_time - t, order=0) for t in times],
        dtype=float,
    )
    symmetry_error = float(np.max(np.abs(forward - backward)))

    print("Two-segment symmetric Phase 1 demo")
    print(f"coeffs shape: {result['coeffs'].shape}")
    print(f"max symmetry error: {symmetry_error:.3e}")
    print(f"system residual: {result['residual_norm']:.3e}")


if __name__ == "__main__":
    main()
