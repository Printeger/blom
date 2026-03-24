from __future__ import annotations

from phase_10.blom_full_backward_diff import T_to_tau, default_weights, representative_case
from phase_10.blom_space_time_optimizer import run_space_time_optimization


def main() -> None:
    q, T = representative_case()
    tau = T_to_tau(T, T_min=0.2)
    run_space_time_optimization(q, tau, default_weights(), n_steps=12, step_size=5e-3, T_min=0.2)


if __name__ == "__main__":
    main()
