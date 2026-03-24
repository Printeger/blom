from phase_8.exp1_large_M_sweep import run_large_M_sweep


if __name__ == "__main__":
    run_large_M_sweep(M_values=[20, 40], k_values=[2, 4, 6, 8], seed=42)
