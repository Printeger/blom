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

## Phase 1 Traceability

This section maps the Phase 1 requirements in `phase_1/REQ_phase_1_minco_scalar_baseline.md`
to the implemented Phase 1 MINCO baseline, demos, and verification assets.

### Scope

Implemented Phase 1 object:

- scalar output only
- degree-7 polynomial per segment
- fixed boundary jets up to order 3 at the start and end
- interior waypoint positions from `q`
- interior `C^6` continuity
- system form `M(T) c = b(q, zeta_start, zeta_end)`

### Requirement Mapping

| Requirement | Implementation | Verification |
| --- | --- | --- |
| Phase 1 scalar MINCO mother problem with fixed boundary jets | `phase_1/minco_scalar_baseline.py` -> `solve_minco_coefficients`, `build_system_matrix`, `build_rhs` | `phase_1/test_minco_scalar_baseline.py` |
| Monomial basis and derivative basis APIs `beta`, `beta_d` | `phase_1/minco_scalar_baseline.py` -> `beta`, `beta_d` | `phase_1/test_minco_scalar_baseline.py` |
| Deterministic coefficient solve returning `coeffs`, `c_vec`, optional system data | `phase_1/minco_scalar_baseline.py` -> `solve_minco_coefficients(..., return_system=True)` | `phase_1/test_minco_scalar_baseline.py` |
| Trajectory evaluation API for segment/global evaluation and sampling | `phase_1/minco_scalar_baseline.py` -> `evaluate_segment`, `evaluate_trajectory`, `sample_trajectory` | `phase_1/test_minco_scalar_baseline.py`, `phase_1/examples/demo_case_two_segment.py` |
| Verification helpers for interpolation, continuity, boundary jets, and system residual | `phase_1/minco_scalar_baseline.py` -> `interpolation_errors`, `continuity_jumps`, `boundary_jet_errors`, `system_residual` | `phase_1/test_minco_scalar_baseline.py`, `phase_1/examples/demo_case_random.py` |
| Sensitivity-ready design with future gradient interfaces | `phase_1/minco_scalar_baseline.py` -> `solve_adjoint`, `grad_wrt_q`, `grad_wrt_T` | Code inspection and API availability |
| Single-segment, symmetric two-segment, and random demos | `phase_1/examples/demo_case_single_segment.py`, `phase_1/examples/demo_case_two_segment.py`, `phase_1/examples/demo_case_random.py` | Direct execution of each demo script |
| Minimal Phase 1 README with math object, shapes, example, and limitations | `phase_1/README_phase1.md` | Manual review |

### Verification Summary

Executed for Phase 1 implementation:

- `python3 -m compileall phase_1`
- `python3 -m unittest phase_1.test_minco_scalar_baseline`
- `python3 -m unittest discover -s . -p 'test*.py'`
- `python3 -m phase_1.examples.demo_case_single_segment`
- `python3 -m phase_1.examples.demo_case_two_segment`
- `python3 -m phase_1.examples.demo_case_random`

Observed Phase 1 results:

- Phase 1 unit tests: `7 passed`
- Full repository test suite: `20 passed`
- Single-segment demo residual: `9.875e-14`
- Two-segment symmetry demo max mirrored-time error: `3.275e-15`
- Random demo max interpolation error: `1.232e-14`
- Random demo max continuity jump: `6.484e-14`
- Random demo max boundary jet error: `1.625e-14`

### Deliverables

Implemented deliverables under `phase_1/`:

- `REQ_phase_1_minco_scalar_baseline.md`
- `__init__.py`
- `minco_scalar_baseline.py`
- `test_minco_scalar_baseline.py`
- `examples/__init__.py`
- `examples/demo_case_single_segment.py`
- `examples/demo_case_two_segment.py`
- `examples/demo_case_random.py`
- `README_phase1.md`
