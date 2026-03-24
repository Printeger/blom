from __future__ import annotations

from phase_7.blom_phase7_extra_experiments import DEFAULT_RESULTS_DIR, run_time_regime_split


def main() -> None:
    result = run_time_regime_split(M=12, n_trials=4, seed=42, save_dir=DEFAULT_RESULTS_DIR / "exp2_time_regime_split")
    print("Phase 7 extra exp2: time-regime split")
    print(f"records: {len(result['records'])}")


if __name__ == "__main__":
    main()

