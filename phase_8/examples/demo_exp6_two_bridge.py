from phase_8.exp6_two_bridge_gap_compare import run_two_bridge_gap_compare


if __name__ == "__main__":
    run_two_bridge_gap_compare(M=20, k_values=[2, 4, 6, 8], n_trials=8, seed=42)
