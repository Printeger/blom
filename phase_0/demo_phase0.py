"""
Minimal end-to-end Phase 0 demo.
"""

from __future__ import annotations

from pprint import pprint

from phase_0.blom_local_qp import solve_blom_local_qp
from phase_0.blom_problem_setup import BLOMProblemSetup, DEFAULT_SEED
from phase_0.minco_scalar_baseline import solve_minco_scalar
from phase_0.phase0_checks import check_blom_local_result, check_minco_result, check_setup


def main() -> None:
    setup = BLOMProblemSetup.make_random(M=6, seed=DEFAULT_SEED, q_scale=0.8)
    minco_result = solve_minco_scalar(setup)
    center_segment = setup.M // 2
    local_result = solve_blom_local_qp(setup, center_segment=center_segment)

    print("Phase 0 setup")
    print(f"q = {setup.q}")
    print(f"T = {setup.T}")
    print()
    print(f"MINCO total cost: {minco_result.cost:.12f}")
    print(f"BLOM local cost (segment {center_segment}): {local_result.cost:.12f}")
    print()
    print("Checks")
    pprint(check_setup(setup))
    pprint(check_minco_result(minco_result))
    pprint(check_blom_local_result(local_result))


if __name__ == "__main__":
    main()
