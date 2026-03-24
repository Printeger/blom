from __future__ import annotations

from phase_6.blom_fd_jacobian_check import DEFAULT_RESULTS_DIR, representative_case, run_compare_all_schemes


def main() -> None:
    q, T = representative_case()
    result = run_compare_all_schemes(q, T, save_dir=DEFAULT_RESULTS_DIR)
    print("Phase 6 demo compare all")
    print(f"scheme C q bandwidth: {result['scheme_C']['q_bandwidth']['max_effective_bandwidth']}")
    print(f"scheme A q bandwidth: {result['scheme_A']['q_bandwidth']['max_effective_bandwidth']}")
    print(f"scheme B q bandwidth: {result['scheme_B']['q_bandwidth']['max_effective_bandwidth']}")


if __name__ == "__main__":
    main()

