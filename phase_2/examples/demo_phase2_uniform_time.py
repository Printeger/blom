from __future__ import annotations

from phase_2.phase2_validation import make_uniform_time_case, run_phase2_validation_suite


def main() -> None:
    case = make_uniform_time_case(8, case_name="uniform_time")
    summary = run_phase2_validation_suite(
        case.q,
        case.T,
        case.zeta_start,
        case.zeta_end,
        case_name=case.case_name,
    )
    print("Phase 2 uniform-time demo")
    print(f"M: {summary['M']}")
    print(f"fro error (exact vs fd): {summary['jacobian_comparison']['fro_error']:.3e}")
    print(f"max effective bandwidth: {summary['bandwidth']['max_effective_bandwidth']}")
    print(f"results dir: {summary['results_dir']['base']}")


if __name__ == "__main__":
    main()

