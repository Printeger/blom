from phase_8.exp4_raw_vs_reference_sanity import run_raw_vs_reference_sanity


if __name__ == "__main__":
    run_raw_vs_reference_sanity(M_values=[20, 40], k_values=[2, 4, 6], n_trials=4, seed=42)
