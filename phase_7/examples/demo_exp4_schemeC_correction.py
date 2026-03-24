from __future__ import annotations

from phase_7.blom_phase7_extra_experiments import DEFAULT_RESULTS_DIR, representative_case, run_schemeC_light_assembly


def main() -> None:
    q, T = representative_case()
    result = run_schemeC_light_assembly(q=q, T=T, save_dir=DEFAULT_RESULTS_DIR / "exp4_schemeC_light_assembly")
    print("Phase 7 extra exp4: Scheme C correction")
    print(f"rows: {len(result['rows'])}")


if __name__ == "__main__":
    main()

