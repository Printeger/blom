from __future__ import annotations

from phase_7.blom_convergence_vs_k import DEFAULT_RESULTS_DIR, representative_case, run_convergence_vs_k


def main() -> None:
    q, T = representative_case()
    result = run_convergence_vs_k(q, T, schemes=("B", "C"), save_dir=DEFAULT_RESULTS_DIR)
    print("Phase 7 ideal-vs-actual demo")
    print(f"ideal final error: {result['final_by_curve']['ideal']['global_l2']:.6e}")
    for scheme in result["schemes"]:
        final_metrics = result["final_by_curve"][scheme]
        print(
            f"scheme {scheme}: final error={final_metrics['global_l2']:.6e}, "
            f"matching={final_metrics['matching_l2']:.6e}"
        )


if __name__ == "__main__":
    main()

