# Phase 0 Traceability

This document maps the Phase 0 requirements in `phase_0/REQ-Phase0.md` to the implemented code, demo, and verification artifacts.

## Scope

Canonical setting implemented in this repository:

- `s = 4`
- `k = 2`
- `d_i = 1`
- `dim = 1`
- natural boundary conditions only

Design note:

- Order-0 continuity is enforced through shared endpoint interpolation on adjacent segments.
- `constraint_builders.py` only adds continuity rows for orders `1..6` to avoid redundant linear constraints in the numerical system.

## Requirement Mapping

| Requirement | Implementation | Verification |
| --- | --- | --- |
| Unified canonical setup object with `(q, T)` interface, validation, `theta`, window access, and JSON round-trip | `phase_0/blom_problem_setup.py` -> `BLOMProblemSetup`, `make_random`, `theta`, `window`, `to_dict`, `from_dict`, `save_json`, `load_json` | `phase_0/tests/test_problem_setup.py`, `phase_0/tests/test_phase0_checks.py` |
| Degree-7 polynomial basis, derivative evaluation, and snap cost matrix | `phase_0/poly_basis.py` -> `monomial_row`, `derivative_row`, `eval_poly`, `snap_cost_matrix`, `hermite_endpoint_matrix`, `coeffs_from_endpoint_derivatives` | Transitively exercised by `test_minco_scalar_baseline.py`, `test_blom_local_qp.py`, `test_phase0_checks.py` |
| Shared global/local constraint assembly for Phase 0 | `phase_0/constraint_builders.py` -> `build_global_minco_constraints`, `build_local_blom_constraints`, `build_snap_hessian` | Transitively exercised by `test_minco_scalar_baseline.py`, `test_blom_local_qp.py`, `test_phase0_checks.py` |
| Global scalar MINCO baseline solver returning coefficients and cost | `phase_0/minco_scalar_baseline.py` -> `MincoScalarResult`, `solve_minco_scalar` | `phase_0/tests/test_minco_scalar_baseline.py`, `phase_0/tests/test_phase0_checks.py` |
| Local Phase 0 BLOM window solver for `s=4, k=2` | `phase_0/blom_local_qp.py` -> `BlomLocalResult`, `solve_blom_local_qp` | `phase_0/tests/test_blom_local_qp.py`, `phase_0/tests/test_phase0_checks.py` |
| Piecewise trajectory evaluation, sampling, and junction jumps | `phase_0/trajectory_eval.py` -> `locate_segment`, `eval_piecewise`, `sample_trajectory`, `junction_jumps` | `phase_0/tests/test_minco_scalar_baseline.py`, `phase_0/tests/test_phase0_checks.py` |
| Centralized setup, MINCO, and local BLOM checks | `phase_0/phase0_checks.py` -> `check_setup`, `check_minco_result`, `check_blom_local_result` | `phase_0/tests/test_phase0_checks.py` |
| End-to-end demo from random setup through MINCO, local BLOM, and checks | `phase_0/demo_phase0.py` | `python3 -m phase_0.demo_phase0` |
| Canonical sample JSON matching implemented setup interface | `phase_0/blom_phase0_problem.json` regenerated from `python3 -m phase_0.blom_problem_setup` | Manual regeneration command executed during verification |

## Verification Summary

Executed during implementation:

- `python3 -m compileall phase_0`
- `python3 -m unittest discover -s phase_0/tests -t .`
- `python3 -m phase_0.blom_problem_setup`
- `python3 -m phase_0.demo_phase0`

Observed acceptance-level results:

- Global MINCO interpolation error: `2.89e-15`
- Global MINCO max junction jump: `1.15e-14`
- Global MINCO max boundary residual: `2.66e-15`
- Local BLOM interpolation error: `4.16e-16`
- Local BLOM max continuity error: `4.16e-16`
- Local BLOM max boundary residual: `0.0`
- Unit tests: `13 passed`

## Deliverables

Implemented deliverables under `phase_0/`:

- `blom_problem_setup.py`
- `poly_basis.py`
- `constraint_builders.py`
- `minco_scalar_baseline.py`
- `blom_local_qp.py`
- `trajectory_eval.py`
- `phase0_checks.py`
- `demo_phase0.py`
- `tests/`
- `blom_phase0_problem.json`
