from __future__ import annotations

import numpy as np

from phase_5.blom_boundary_jump_check import DEFAULT_RESULTS_DIR, run_boundary_jump_check


def main() -> None:
    q = np.asarray([0.0, 1.0, -0.8, 1.5, 0.2, 0.7], dtype=float)
    T = np.asarray([0.8, 1.4, 0.9, 1.1, 0.75], dtype=float)
    result = run_boundary_jump_check(q, T, scheme="B", save_dir=DEFAULT_RESULTS_DIR)
    print("Phase 5 demo scheme B")
    print(f"max lower-order jump: {result['stats_overview']['lower_order_max']:.3e}")
    print(f"max pre-consensus dispersion: {np.max(result['consensus']['pre_dispersion']):.3e}")


if __name__ == "__main__":
    main()

