from __future__ import annotations

import numpy as np

from phase_1.minco_scalar_baseline import solve_minco_coefficients


def main() -> None:
    q = np.array([0.0, 1.0])
    T = np.array([1.5])
    zeta_start = np.array([0.0, 0.0, 0.0, 0.0])
    zeta_end = np.array([1.0, 0.0, 0.0, 0.0])

    result = solve_minco_coefficients(q, T, zeta_start, zeta_end, return_system=True)
    print("Single-segment Phase 1 demo")
    print(f"coeffs shape: {result['coeffs'].shape}")
    print(f"coefficients: {result['coeffs'][0]}")
    print(f"system residual: {result['residual_norm']:.3e}")


if __name__ == "__main__":
    main()

