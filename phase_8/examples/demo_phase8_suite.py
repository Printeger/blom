from __future__ import annotations

from phase_8.blom_phase8_validation_suite import run_phase8_validation_suite
from phase_8.phase8_common import DEFAULT_RESULTS_DIR


def main() -> None:
    run_phase8_validation_suite(M=10, n_trials=6, save_dir=DEFAULT_RESULTS_DIR)


if __name__ == "__main__":
    main()

