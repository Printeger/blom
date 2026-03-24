from __future__ import annotations

from phase_6.blom_fd_jacobian_check import DEFAULT_RESULTS_DIR, representative_case, run_fd_jacobian_check


def main() -> None:
    q, T = representative_case()
    result = run_fd_jacobian_check(q, T, scheme="C", save_dir=DEFAULT_RESULTS_DIR)
    raw = result["raw"]
    print("Phase 6 demo raw local support")
    print(f"raw J_c_q max abs error: {raw['errors']['J_c_q']['max_abs_error']:.3e}")
    print(f"raw J_c_T max abs error: {raw['errors']['J_c_T']['max_abs_error']:.3e}")
    print(f"raw q mask pass: {raw['mask_pass']['J_c_q']}")


if __name__ == "__main__":
    main()

