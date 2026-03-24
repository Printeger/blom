from __future__ import annotations

from phase_7.blom_phase7_extra_experiments import DEFAULT_RESULTS_DIR, representative_case, run_phase7_extra_experiments


def main() -> None:
    q, T = representative_case()
    result = run_phase7_extra_experiments(q=q, T=T, M_values=[10, 20, 40], n_trials=4, save_dir=DEFAULT_RESULTS_DIR, seed=42)
    print("Phase 7 extra experiments")
    print(f"overview rows: {len(result['overview_rows'])}")


if __name__ == "__main__":
    main()
