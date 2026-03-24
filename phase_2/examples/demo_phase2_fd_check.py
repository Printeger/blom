from __future__ import annotations

from phase_2.phase2_validation import make_nonuniform_time_case, run_phase2_validation_suite


def main() -> None:
    case = make_nonuniform_time_case(6, seed=7, case_name="fd_check")
    summary = run_phase2_validation_suite(
        case.q,
        case.T,
        case.zeta_start,
        case.zeta_end,
        case_name=case.case_name,
        eps=1e-7,
    )
    comparison = summary["jacobian_comparison"]
    print("Phase 2 finite-difference check demo")
    print(f"fro error: {comparison['fro_error']:.3e}")
    print(f"max abs error: {comparison['max_abs_error']:.3e}")
    print(f"relative error: {comparison['relative_error']:.3e}")
    print(f"results dir: {summary['results_dir']['base']}")


if __name__ == "__main__":
    main()

