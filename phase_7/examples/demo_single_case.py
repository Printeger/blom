from __future__ import annotations

from phase_7.blom_convergence_vs_k import DEFAULT_RESULTS_DIR, representative_case, run_convergence_vs_k


def main() -> None:
    q, T = representative_case()
    result = run_convergence_vs_k(q, T, save_dir=DEFAULT_RESULTS_DIR)
    print("Phase 7 single-case convergence demo")
    print(f"k grid: {result['inputs']['k_values']}")
    print(f"ideal final error: {result['final_by_curve']['ideal']['global_l2']:.6e}")
    for scheme in result["schemes"]:
        print(f"scheme {scheme} final error: {result['final_by_curve'][scheme]['global_l2']:.6e}")


if __name__ == "__main__":
    main()

