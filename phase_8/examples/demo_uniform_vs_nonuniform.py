from __future__ import annotations

from phase_8.blom_uniform_vs_nonuniform_interior import run_uniform_vs_nonuniform_interior
from phase_8.phase8_common import DEFAULT_RESULTS_DIR, default_q_sampler


def main() -> None:
    run_uniform_vs_nonuniform_interior(
        q_sampler=default_q_sampler,
        M=10,
        k_values=[2, 4, 6, 8, 10],
        h_values=[0.5, 1.0, 2.0],
        nonuniform_boxes=[(0.9, 1.1), (0.5, 2.0), (0.2, 3.0)],
        n_trials=6,
        save_dir=DEFAULT_RESULTS_DIR,
    )


if __name__ == "__main__":
    main()

