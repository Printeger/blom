from __future__ import annotations

from phase_9.blom_backward_diff import default_obs_config, default_weights, representative_case
from phase_9.blom_space_time_opt_demo import run_minimal_optimization_demo


def main() -> None:
    q, T = representative_case()
    run_minimal_optimization_demo(q, T, default_weights(), obs_config=default_obs_config(), n_steps=25, step_size=5e-3, T_min=0.2)


if __name__ == "__main__":
    main()

