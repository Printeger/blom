from __future__ import annotations

from phase_10.blom_benchmark_suite import run_benchmark_suite
from phase_10.blom_full_backward_diff import default_weights, representative_case


def main() -> None:
    q, T = representative_case()
    run_benchmark_suite(q0=q, T0=T, weights=default_weights(), M_values=[6], k_values=[2], benchmark_methods=["raw_schemeC", "minco", "schemeA", "heuristic"], n_steps=4, T_min=0.2)


if __name__ == "__main__":
    main()
