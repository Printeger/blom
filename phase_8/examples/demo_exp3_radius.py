from phase_8.exp3_radius_sensitivity import run_radius_sensitivity


if __name__ == "__main__":
    run_radius_sensitivity(M=40, k_values=[2, 4, 6, 8], radius_modes=["half_k", "k", "three_half_k"], seed=42)
