from __future__ import annotations

import numpy as np

from phase_5.blom_boundary_jump_check import DEFAULT_RESULTS_DIR, run_boundary_jump_check


def main() -> None:
    q = np.asarray([0.0, 1.0, -0.8, 1.5, 0.2, 0.7], dtype=float)
    T = np.asarray([0.8, 1.4, 0.9, 1.1, 0.75], dtype=float)
    result = run_boundary_jump_check(q, T, scheme="C", save_dir=DEFAULT_RESULTS_DIR)
    print("Phase 5 demo scheme C")
    print(f"position jump max: {result['stats'][0]['max_abs']:.3e}")
    print(f"velocity/accel/jerk jump max: {max(result['stats'][order]['max_abs'] for order in range(1, 4)):.3e}")


if __name__ == "__main__":
    main()

