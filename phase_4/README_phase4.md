# Phase 4 Analytic Validation

This directory implements the Phase 4 validation layer that turns the Phase 3
local variational definition into a concrete analytic local solver prototype.

## Scope

Phase 4 has three roles:

- verify the exact warm-up model `(s, k) = (2, 2)` symbolically,
- build the exact canonical local analytic system for `(s, k) = (4, 2)`,
- compare the exact analytic local model against Catmull-style heuristics.

The goal is not to build a full BLOM planner yet. The goal is to validate the
analytic core that later BLOM-Analytic implementations can reuse.

## Main Scripts

- `blom_k2_s2_sympy.py`: symbolic derivation and verification that time-weighted Catmull--Rom is exact for `s = 2`
- `blom_k2_s4_numeric.py`: exact `6x6` local analytic system for the canonical `s = 4, k = 2` case
- `blom_catmull_compare.py`: exact-vs-heuristic comparison across random and representative samples

## Shared Utilities

- `utils/hermite_utils.py`: Hermite matrices, Gram matrices, scaling operators, reconstruction helpers
- `utils/io_utils.py`: text / json / npy / csv saving
- `utils/plotting_utils.py`: white-background plotting helpers for exported figures

## Running the Scripts

```bash
python3 -m phase_4.blom_k2_s2_sympy
python3 -m phase_4.blom_k2_s4_numeric
python3 -m phase_4.blom_catmull_compare
```

## Running the Tests

```bash
python3 -m unittest phase_4.test_phase4_analytic
```

## Result Directories

```text
phase_4/results/
├── s2_sympy/
├── s4_numeric/
└── catmull_compare/
```

## How to Read the Outputs

- `s2_sympy`: if `v_exact - v_catmull` simplifies to zero, then Catmull is exactly the BLOM local optimum in the `s = 2` warm-up model.
- `s4_numeric`: if `A2(T)` is symmetric positive definite and the analytic solution matches the Phase 3 local QP, then the analytic `6x6` core is correct.
- `catmull_compare`: if `s = 2` errors stay at machine precision while `s = 4` errors and objective gaps stay nonzero, then Catmull is only exact in the warm-up model and must be replaced in the canonical minimum-snap case.

## Important Scope Boundary

The heuristic completion used in `blom_catmull_compare.py` for `s = 4` is only
a baseline:

- velocities come from time-weighted Catmull-style formulas,
- accelerations and jerks are set to zero,
- this baseline is **not** exact theory.

It exists only to make the exact analytic local system visibly comparable to a
familiar heuristic.

