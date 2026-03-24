from __future__ import annotations

from phase_7.blom_phase7_extra_experiments import DEFAULT_RESULTS_DIR, representative_case, run_interior_vs_full


def main() -> None:
    q, T = representative_case()
    result = run_interior_vs_full(q=q, T=T, n_trials=4, save_dir=DEFAULT_RESULTS_DIR / "exp3_interior_vs_full")
    print("Phase 7 extra exp3: interior vs full")
    print(f"k grid: {result['case']['k_values']}")


if __name__ == "__main__":
    main()

