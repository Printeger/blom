# Phase 5: Boundary Jump Check

Phase 5 closes the gap between local BLOM windows and a full assembled global trajectory.
The implementation in [blom_boundary_jump_check.py](/home/dev/code/blom/phase_5/blom_boundary_jump_check.py) compares three assembly mechanisms:

- Scheme A: global shared junction states
- Scheme B: independent local windows plus consensus
- Scheme C: raw central-segment extraction

## Main API

- `compute_jumps(coeffs, T, s, max_order=None)`
- `summarize_jumps(jumps)`
- `assemble_scheme_A(...)`
- `assemble_scheme_B(...)`
- `assemble_scheme_C(...)`
- `run_boundary_jump_check(...)`
- `run_compare_all_schemes(...)`
- `run_random_trials(...)`

## Reused building blocks

- Phase 1: polynomial derivative evaluation and trajectory conventions
- Phase 3: BLOM-Strict local window solve for independent window predictions
- Phase 4: exact Hermite local-state reconstruction and segment-energy quadratic forms

## Expected mathematical behavior

- Scheme A and Scheme B should achieve exact global `C^{s-1}` continuity after assembly.
- Scheme C only guarantees exact waypoint continuity `C^0`; higher-order jumps are diagnostic outputs.

## How To Run

```bash
python3 -m phase_5.blom_boundary_jump_check
python3 -m phase_5.examples.demo_scheme_A
python3 -m phase_5.examples.demo_scheme_B
python3 -m phase_5.examples.demo_scheme_C
python3 -m phase_5.examples.demo_compare_all
python3 -m unittest phase_5.test_blom_boundary_jump_check
```

## Artifacts

Artifacts are written under [phase_5/results/phase5_boundary_jump_check](/home/dev/code/blom/phase_5/results/phase5_boundary_jump_check):

- per-scheme heatmaps and max-bar figures
- lower-order comparison plot
- Scheme C eta mismatch plot
- Scheme B consensus-improvement plot
- random-trial boxplot
- JSON/CSV summaries and interpretation markdown

## Interpretation

Scheme B is the preferred practical assembly mechanism in this repository: it keeps the local-window solve phase independent, then restores exact global `C^{s-1}` continuity with a closed-form consensus step.
