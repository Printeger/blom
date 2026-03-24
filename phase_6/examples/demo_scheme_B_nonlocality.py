from __future__ import annotations

from phase_6.blom_fd_jacobian_check import DEFAULT_RESULTS_DIR, representative_case, run_compare_all_schemes


def main() -> None:
    q, T = representative_case()
    result = run_compare_all_schemes(q, T, save_dir=DEFAULT_RESULTS_DIR)
    scheme = result["scheme_B"]
    print("Phase 6 demo scheme B nonlocality")
    print(f"q bandwidth: {scheme['q_bandwidth']['max_effective_bandwidth']}")
    print(f"outside-band q max abs: {scheme['outside_band']['q_max_abs']:.3e}")


if __name__ == "__main__":
    main()

