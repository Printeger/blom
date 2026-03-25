from __future__ import annotations

from phase_10.blom_full_backward_diff import T_to_tau, default_weights, representative_case
from phase_10.blom_phase10_framework_suite import run_phase10_framework_suite


def main() -> None:
    q, T = representative_case()
    tau = T_to_tau(T, T_min=0.2)
    run_phase10_framework_suite(q0=q, tau0=tau, weights=default_weights(), M_values=[6, 10], k_values=[2, 4], benchmark_methods=["raw_schemeC", "heuristic"], n_steps=8, T_min=0.2)


if __name__ == "__main__":
    main()
