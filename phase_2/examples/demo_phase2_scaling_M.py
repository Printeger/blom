from __future__ import annotations

from phase_2.phase2_validation import run_scaling_experiment


def main() -> None:
    summary = run_scaling_experiment(M_values=(4, 8, 16, 32), case_name="uniform_time")
    print("Phase 2 scaling-M demo")
    for record in summary["records"]:
        print(
            f"M={record['M']}: max_bw={record['max_effective_bandwidth']}, "
            f"mean_bw={record['mean_effective_bandwidth']:.2f}, "
            f"far_ratio={record['far_nonzero_ratio']:.3e}"
        )
    print(f"results dir: {summary['results_dir']['base']}")


if __name__ == "__main__":
    main()

