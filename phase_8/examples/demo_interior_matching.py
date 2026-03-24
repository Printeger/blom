from __future__ import annotations

from phase_8.blom_interior_matching_check import run_interior_matching_check
from phase_8.phase8_common import DEFAULT_RESULTS_DIR, prepare_case


def main() -> None:
    q, T = prepare_case()
    run_interior_matching_check(q, T, [2, 4, 6, 8], save_dir=DEFAULT_RESULTS_DIR)


if __name__ == "__main__":
    main()

