from __future__ import annotations

from phase_7.blom_convergence_vs_k import DEFAULT_RESULTS_DIR, run_random_trials


def main() -> None:
    result = run_random_trials(n_trials=8, M=10, save_dir=DEFAULT_RESULTS_DIR, seed=42)
    print("Phase 7 random-trials demo")
    for curve, summary in result["summary"]["curves"].items():
        print(
            f"{curve}: success={summary['success_rate']:.2f}, "
            f"mean slope={summary['mean_slope']:.6e}, "
            f"mean R^2={summary['mean_r2']:.6f}"
        )


if __name__ == "__main__":
    main()

