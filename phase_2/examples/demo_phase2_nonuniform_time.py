from __future__ import annotations

from phase_2.phase2_validation import make_nonuniform_time_case, run_phase2_validation_suite


def main() -> None:
    case = make_nonuniform_time_case(8, seed=42, case_name="nonuniform_time")
    summary = run_phase2_validation_suite(
        case.q,
        case.T,
        case.zeta_start,
        case.zeta_end,
        case_name=case.case_name,
    )
    print("Phase 2 nonuniform-time demo")
    print(f"M: {summary['M']}")
    print(f"relative Jacobian error: {summary['jacobian_comparison']['relative_error']:.3e}")
    print(f"far-field ratio: {summary['bandwidth']['far_nonzero_ratio']:.3e}")
    print(f"results dir: {summary['results_dir']['base']}")


if __name__ == "__main__":
    main()

