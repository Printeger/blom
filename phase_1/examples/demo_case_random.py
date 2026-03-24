from __future__ import annotations

import numpy as np

from phase_1.minco_scalar_baseline import (
    boundary_jet_errors,
    continuity_jumps,
    interpolation_errors,
    solve_minco_coefficients,
)


def main() -> None:
    rng = np.random.default_rng(42)
    M_seg = 5
    q = np.concatenate(([0.0], rng.normal(scale=0.8, size=M_seg)))
    T = rng.uniform(0.5, 1.5, size=M_seg)
    zeta_start = np.array([q[0], 0.0, 0.0, 0.0])
    zeta_end = np.array([q[-1], 0.0, 0.0, 0.0])

    result = solve_minco_coefficients(q, T, zeta_start, zeta_end, return_system=True)
    interp = interpolation_errors(result["coeffs"], T, q)
    jumps = continuity_jumps(result["coeffs"], T)
    bc = boundary_jet_errors(result["coeffs"], T, zeta_start, zeta_end)

    print("Random Phase 1 demo")
    print(f"q: {q}")
    print(f"T: {T}")
    print(f"system residual: {result['residual_norm']:.3e}")
    print(f"max interpolation error: {np.max(np.abs(interp)):.3e}")
    print(f"max continuity jump: {np.max(np.abs(jumps)):.3e}")
    print(f"max boundary jet error: {np.max(np.abs(bc)):.3e}")


if __name__ == "__main__":
    main()

