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

## Phase 2 Traceability

This section maps the Phase 2 requirements in `phase_2/REQ-Phase2-validation.md`
to the implemented validation toolbox, result artifacts, and verification assets.

### Scope

Implemented Phase 2 object:

- exact waypoint Jacobian extraction on top of the fixed Phase 1 global system
- finite-difference Jacobian validation
- block-level coefficient influence analysis
- effective-bandwidth estimation as `M` scales
- automatic saving of figures, CSV tables, and JSON logs

### Requirement Mapping

| Requirement | Implementation | Verification |
| --- | --- | --- |
| Waypoint selector `S_q` compatible with the Phase 1 system row ordering | `phase_2/phase2_validation.py` -> `build_waypoint_selector` | `phase_2/test_phase2_validation.py` |
| Exact Jacobian extraction `J_q = A(T)^{-1} S_q` without explicit inverse | `phase_2/phase2_validation.py` -> `compute_exact_jacobian_q` | `phase_2/test_phase2_validation.py`, demo scripts |
| Finite-difference Jacobian and exact-vs-FD comparison | `phase_2/phase2_validation.py` -> `finite_difference_jacobian_q`, `compare_exact_vs_fd_jacobian` | `phase_2/test_phase2_validation.py`, `phase_2/examples/demo_phase2_fd_check.py` |
| Segment-wise block influence norms and single-waypoint influence profile | `phase_2/phase2_validation.py` -> `compute_segmentwise_influence_norms`, `compute_waypoint_influence_profile` | `phase_2/test_phase2_validation.py`, uniform/nonuniform demos |
| Effective-bandwidth estimation and scaling-in-M analysis | `phase_2/phase2_validation.py` -> `estimate_effective_bandwidth`, `run_scaling_experiment` | `phase_2/test_phase2_validation.py`, `phase_2/examples/demo_phase2_scaling_M.py` |
| Automatic validation suite with figure/table/log export | `phase_2/phase2_validation.py` -> `run_phase2_validation_suite`, `ensure_results_dirs` | `phase_2/test_phase2_validation.py`, demo scripts |
| Matrix sparsity plot, Jacobian heatmap, block heatmap, influence profile, FD comparison, scaling plot | `phase_2/phase2_plotting.py` -> `visualize_matrix_sparsity`, `visualize_heatmap`, `visualize_block_sparsity_Jq`, `visualize_waypoint_influence_profile`, `visualize_jacobian_fd_compare`, `visualize_effective_bandwidth_vs_M` | saved artifacts under `phase_2/results/figures/` |
| Example scripts for uniform time, nonuniform time, FD check, and scaling in `M` | `phase_2/examples/demo_phase2_uniform_time.py`, `phase_2/examples/demo_phase2_nonuniform_time.py`, `phase_2/examples/demo_phase2_fd_check.py`, `phase_2/examples/demo_phase2_scaling_M.py` | direct execution of each script |
| Phase 2 README with purpose, interfaces, outputs, and run instructions | `phase_2/README_phase2.md` | manual review |

### Verification Summary

Executed for Phase 2 implementation:

- `python3 -m compileall phase_2`
- `python3 -m unittest phase_2.test_phase2_validation`
- `python3 -m unittest discover -s . -p 'test*.py'`
- `python3 -m phase_2.examples.demo_phase2_uniform_time`
- `python3 -m phase_2.examples.demo_phase2_nonuniform_time`
- `python3 -m phase_2.examples.demo_phase2_fd_check`
- `python3 -m phase_2.examples.demo_phase2_scaling_M`

Observed Phase 2 results:

- Phase 2 unit tests: `5 passed`
- Full repository test suite: `25 passed`
- Uniform-time demo Frobenius Jacobian error: `1.179e-07`
- Uniform-time demo max effective bandwidth at `M=8`: `7`
- Nonuniform-time demo relative Jacobian error: `2.397e-08`
- Nonuniform-time demo far-field ratio: `1.000e+00`
- FD-check demo Frobenius Jacobian error: `2.168e-05`
- FD-check demo max absolute Jacobian error: `1.069e-05`
- Scaling demo max effective bandwidth sequence: `3, 7, 15, 31` for `M = 4, 8, 16, 32`

### Deliverables

Implemented deliverables under `phase_2/`:

- `REQ-Phase2-validation.md`
- `__init__.py`
- `phase2_validation.py`
- `phase2_plotting.py`
- `test_phase2_validation.py`
- `examples/__init__.py`
- `examples/demo_phase2_uniform_time.py`
- `examples/demo_phase2_nonuniform_time.py`
- `examples/demo_phase2_scaling_M.py`
- `examples/demo_phase2_fd_check.py`
- `README_phase2.md`
- `results/figures/`
- `results/tables/`
- `results/logs/`
