from __future__ import annotations

from phase_9.blom_backward_diff import default_obs_config, default_weights, representative_case, run_phase9_gradcheck


def main() -> None:
    q, T = representative_case()
    run_phase9_gradcheck(q, T, default_weights(), obs_config=default_obs_config())


if __name__ == "__main__":
    main()

