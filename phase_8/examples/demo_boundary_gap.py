from __future__ import annotations

from phase_8.blom_boundary_gap_decomposition import run_boundary_gap_decomposition
from phase_8.phase8_common import DEFAULT_RESULTS_DIR, prepare_case


def main() -> None:
    q, T = prepare_case()
    run_boundary_gap_decomposition(q, T, [2, 4, 6, 8], save_dir=DEFAULT_RESULTS_DIR)


if __name__ == "__main__":
    main()

