from __future__ import annotations

from phase_7.blom_phase7_extra_experiments import DEFAULT_RESULTS_DIR, run_large_M_sweep


def main() -> None:
    result = run_large_M_sweep(M_values=[10, 20, 40], seed=42, save_dir=DEFAULT_RESULTS_DIR / "exp1_large_M_sweep")
    print("Phase 7 extra exp1: large-M sweep")
    print(f"cases: {sorted(result['records'])}")


if __name__ == "__main__":
    main()
