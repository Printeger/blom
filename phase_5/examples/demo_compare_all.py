from __future__ import annotations

import numpy as np

from phase_5.blom_boundary_jump_check import DEFAULT_RESULTS_DIR, run_compare_all_schemes


def main() -> None:
    q = np.asarray([0.0, 1.1, -0.6, 1.4, -0.2, 0.9, 0.1], dtype=float)
    T = np.asarray([0.8, 1.4, 0.9, 1.3, 0.7, 1.1], dtype=float)
    result = run_compare_all_schemes(q, T, save_dir=DEFAULT_RESULTS_DIR)
    print("Phase 5 demo compare all")
    for scheme in ("A", "B", "C"):
        stats = result["scheme_results"][scheme]["stats_overview"]
        print(
            f"scheme {scheme}: lower={stats['lower_order_max']:.3e}, "
            f"higher={stats['higher_order_max']:.3e}"
        )


if __name__ == "__main__":
    main()

