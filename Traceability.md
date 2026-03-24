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

## Phase 3 Traceability

This section maps the Phase 3 requirements in `phase_3/REQ-Phase3-validation.md`
to the implemented BLOM-Strict local QP toolbox, validation artifacts, and verification assets.

### Scope

Implemented Phase 3 object:

- BLOM-Strict local variational problem under the canonical setting `s = 4`, `k = 2`
- explicit local feasible-set construction
- equality-constrained local QP solve by both KKT and reduced-space methods
- natural-boundary, higher-continuity, perturbation-continuity, and multistart-uniqueness validation
- automatic export of Phase 3 figures, CSV tables, and JSON logs

### Requirement Mapping

| Requirement | Implementation | Verification |
| --- | --- | --- |
| Local window definition `W(i,k)` and main local problem assembly | `phase_3/blom_strict_local_qp.py` -> `build_window`, `build_local_problem`, `extract_segment_coeff` | `phase_3/test_phase3_blom_strict.py`, `python3 -m phase_3.validate_phase3_blom_strict` |
| Explicit feasible local spline showing admissible-set nonemptiness | `phase_3/blom_strict_feasible_init.py` -> `build_feasible_local_spline` | `phase_3/test_phase3_blom_strict.py`, feasibility rows in `phase_3/results/tables/table_feasibility_summary.csv` |
| Local Hessian, affine constraints, KKT solve, and reduced-space solve | `phase_3/blom_strict_local_kkt.py` -> `build_local_hessian`, `build_local_constraints`, `solve_kkt`, `solve_reduced_qp` | `phase_3/test_phase3_blom_strict.py`, case logs under `phase_3/results/logs/` |
| Natural-boundary, interpolation, boundary-jet, and continuity diagnostics | `phase_3/blom_strict_local_qp.py` -> `compute_interpolation_errors`, `compute_continuity_jumps`, `compute_boundary_jet_errors`, `compute_natural_boundary_residuals`, `summarize_solution` | feasibility table, case logs, `RESULT_PHASE3.md` |
| Batch validation over left / interior / right windows, perturbation continuity, multistart uniqueness, and full-window Phase 1 comparison | `phase_3/validate_phase3_blom_strict.py` -> `run_phase3_validation`, `_run_perturbation_continuity`, `_run_multistart`, `_full_window_phase1_comparison` | `python3 -m phase_3.validate_phase3_blom_strict`, `phase_3/results/tables/*.csv`, `phase_3/results/logs/phase3_validation_summary.json` |
| Required Phase 3 figures for window layout, local trajectory, continuity jumps, natural BC residuals, perturbation response, and multistart uniqueness | `phase_3/phase3_plotting.py` -> `plot_window_layout`, `plot_local_trajectory`, `plot_continuity_jumps`, `plot_natural_bc_residuals`, `plot_perturbation_response`, `plot_uniqueness_multistart` | saved artifacts under `phase_3/results/figures/` |
| Phase 3 README describing BLOM-Strict scope, interfaces, run commands, artifact locations, and interpretation of the four diagnostic families | `phase_3/README_phase3.md` | manual review |

### Verification Summary

Executed for Phase 3 implementation:

- `python3 -m compileall phase_3`
- `python3 -m unittest phase_3.test_phase3_blom_strict`
- `python3 -m unittest discover -s . -p 'test*.py'`
- `python3 -m phase_3.validate_phase3_blom_strict`

Observed Phase 3 results:

- Phase 3 unit tests: `5 passed`
- Full repository test suite: `30 passed`
- Batch validation cases: `6`
- Max interpolation error across Phase 3 cases: `1.701e-14`
- Max natural-boundary residual across Phase 3 cases: `8.983e-11`
- Max multistart coefficient difference across Phase 3 cases: `4.609e-12`
- Full-window local-vs-Phase-1 coefficient difference: `3.781e-11`

### Deliverables

Implemented deliverables under `phase_3/`:

- `REQ-Phase3-validation.md`
- `__init__.py`
- `blom_strict_local_qp.py`
- `blom_strict_local_kkt.py`
- `blom_strict_feasible_init.py`
- `phase3_plotting.py`
- `validate_phase3_blom_strict.py`
- `test_phase3_blom_strict.py`
- `README_phase3.md`
- `results/figures/`
- `results/tables/`
- `results/logs/`

## Phase 4 Traceability

This section maps the Phase 4 requirements in `phase_4/REQ-Phase4-analytic-validation.md`
to the implemented analytic-validation scripts, utilities, and verification assets.

### Scope

Implemented Phase 4 object:

- symbolic exact verification of the warm-up model `(s, k) = (2, 2)`
- exact numeric analytic core for the canonical local system `(s, k) = (4, 2)`
- explicit comparison between the exact analytic model and Catmull-style heuristics
- artifact export for text expressions, matrices, arrays, CSV tables, JSON summaries, and figures

### Requirement Mapping

| Requirement | Implementation | Verification |
| --- | --- | --- |
| Symbolic one-sided natural cubic derivation, local energy derivation, exact velocity solve, and Catmull equivalence proof for `s=2` | `phase_4/blom_k2_s2_sympy.py` -> `build_left_natural_cubic_symbolic`, `build_right_natural_cubic_symbolic`, `derive_s2_local_energy_symbolic`, `derive_s2_exact_velocity_symbolic`, `verify_catmull_equivalence_symbolic` | `phase_4/test_phase4_analytic.py`, `python3 -m phase_4.blom_k2_s2_sympy`, files under `phase_4/results/s2_sympy/` |
| Shared Hermite matrices, Gram matrices, scaling operators, one-sided rank-one costs, and reconstruction utilities | `phase_4/utils/hermite_utils.py` | Transitively exercised by all Phase 4 scripts and tests |
| Exact `6x6` local analytic system for canonical `s=4,k=2` and Hermite center-segment reconstruction | `phase_4/blom_k2_s4_numeric.py` -> `left_outer_cost_matrix`, `right_outer_cost_matrix`, `build_C4`, `build_G4`, `build_R4`, `middle_segment_cost_matrix`, `build_local_quadratic_s4_k2`, `solve_local_system_s4_k2`, `hermite_reconstruct_center_segment` | `phase_4/test_phase4_analytic.py`, `python3 -m phase_4.blom_k2_s4_numeric`, arrays and figures under `phase_4/results/s4_numeric/` |
| Direct numeric comparison between Phase 4 analytic core and Phase 3 BLOM-Strict local QP | `phase_4/blom_k2_s4_numeric.py` -> `compare_with_blom_strict_qp` | `phase_4/test_phase4_analytic.py`, `phase_4/results/s4_numeric/analytic_vs_qp_metrics.json` |
| Time-weighted Catmull baseline and exact-vs-heuristic comparison for `s=2` and `s=4` | `phase_4/blom_catmull_compare.py` -> `catmull_velocity_time_weighted`, `heuristic_local_state_s4`, `compare_s2_exact_vs_catmull`, `compare_s4_exact_vs_catmull`, `run_random_benchmark`, `plot_comparison_figures` | `phase_4/test_phase4_analytic.py`, `python3 -m phase_4.blom_catmull_compare`, files under `phase_4/results/catmull_compare/` |
| Shared artifact saving and white-background plotting utilities | `phase_4/utils/io_utils.py`, `phase_4/utils/plotting_utils.py` | Artifact generation during all Phase 4 script runs |
| Phase 4 README documenting scope, scripts, outputs, and interpretation | `phase_4/README_phase4.md` | manual review |

### Verification Summary

Executed for Phase 4 implementation:

- `python3 -m compileall phase_4`
- `python3 -m unittest phase_4.test_phase4_analytic`
- `python3 -m unittest discover -s . -p 'test*.py'`
- `python3 -m phase_4.blom_k2_s2_sympy`
- `python3 -m phase_4.blom_k2_s4_numeric`
- `python3 -m phase_4.blom_catmull_compare`

Observed Phase 4 results:

- Phase 4 unit tests: `5 passed`
- Full repository test suite: `35 passed`
- `s=2` symbolic difference `v_exact - v_catmull`: `0`
- `s=2` random max exact-vs-Catmull velocity error: `1.332e-15`
- `s=4` analytic matrix symmetry error: `0.0`
- `s=4` minimum eigenvalue of `A2(T)` on representative sample: `3.139`
- `s=4` analytic-vs-Phase-3 local-state error: `8.565e-11`
- `s=4` analytic-vs-Phase-3 center-coefficient error: `1.292e-11`
- Catmull comparison `s=4` mean center-coefficient error: `4.766e+02`
- Catmull comparison `s=4` minimum objective gap `J_heuristic - J_exact`: `0.611`

### Deliverables

Implemented deliverables under `phase_4/`:

- `REQ-Phase4-analytic-validation.md`
- `__init__.py`
- `blom_k2_s2_sympy.py`
- `blom_k2_s4_numeric.py`
- `blom_catmull_compare.py`
- `test_phase4_analytic.py`
- `README_phase4.md`
- `utils/`
- `results/s2_sympy/`
- `results/s4_numeric/`
- `results/catmull_compare/`

## Phase 5 Traceability

This section maps the Phase 5 requirements in `phase_5/REQ-Phase5-boundary-jump-check.md`
to the implemented global-assembly validation code, tests, and generated artifacts.

### Scope

Implemented Phase 5 object:

- a unified jump-computation and jump-summary layer for assembled piecewise polynomials
- three global assembly schemes: shared junction states, overlapping consensus, and raw central extraction
- per-scheme plotting and CSV/JSON summary export
- comparison and random-trial tooling for smoothness, locality, and stability diagnostics

### Requirement Mapping

| Requirement | Implementation | Verification |
| --- | --- | --- |
| Unified derivative jump computation `[[p^(ell)]](t_i)` up to order `2s-2` | `phase_5/blom_boundary_jump_check.py` -> `compute_jumps` | `phase_5/test_blom_boundary_jump_check.py`, scheme summaries under `phase_5/results/phase5_boundary_jump_check/` |
| Unified jump statistics interface with max/mean/median/RMS/q95/zero-tolerance flags | `phase_5/blom_boundary_jump_check.py` -> `summarize_jumps` | `phase_5/test_blom_boundary_jump_check.py`, `jump_stats_scheme_*.csv` |
| Scheme A shared-state assembly with global quadratic solve and Hermite reconstruction | `phase_5/blom_boundary_jump_check.py` -> `assemble_scheme_A`, `_build_scheme_a_system`, `_reconstruct_from_eta` | `phase_5/test_blom_boundary_jump_check.py`, `scheme_A/summary_scheme_A.json`, `scheme_A/jump_heatmap_scheme_A.png` |
| Scheme B local-window predictions plus closed-form consensus projection | `phase_5/blom_boundary_jump_check.py` -> `assemble_scheme_B` | `phase_5/test_blom_boundary_jump_check.py`, `scheme_B/summary_scheme_B.json`, `scheme_B/scheme_B_consensus_improvement.png` |
| Scheme C raw central-segment extraction and eta mismatch diagnostics | `phase_5/blom_boundary_jump_check.py` -> `assemble_scheme_C` | `phase_5/test_blom_boundary_jump_check.py`, `scheme_C/summary_scheme_C.json`, `scheme_C/scheme_C_eta_mismatch.png` |
| High-level single-scheme runner and all-schemes comparison runner | `phase_5/blom_boundary_jump_check.py` -> `run_boundary_jump_check`, `run_compare_all_schemes` | `phase_5/examples/demo_compare_all.py`, `scheme_comparison_summary.csv`, `phase5_interpretation_summary.md` |
| Randomized multi-trial evaluation | `phase_5/blom_boundary_jump_check.py` -> `run_random_trials` | `random_trials_summary.csv`, `random_trials_aggregate.json`, `jump_boxplot_random_trials.png` |
| Example entry points for schemes A/B/C and the all-schemes comparison | `phase_5/examples/demo_scheme_A.py`, `phase_5/examples/demo_scheme_B.py`, `phase_5/examples/demo_scheme_C.py`, `phase_5/examples/demo_compare_all.py` | direct script execution |
| Phase 5 README documenting scope, interfaces, and artifacts | `phase_5/README_phase5.md` | manual review |

### Verification Summary

Executed for Phase 5 implementation:

- `python3 -m compileall phase_5`
- `python3 -m phase_5.blom_boundary_jump_check`
- `python3 -m phase_5.examples.demo_scheme_A`
- `python3 -m phase_5.examples.demo_scheme_B`
- `python3 -m phase_5.examples.demo_scheme_C`
- `python3 -m phase_5.examples.demo_compare_all`
- `python3 -m unittest phase_5.test_blom_boundary_jump_check`
- `python3 -m unittest discover -s . -p 'test*.py'`

Observed Phase 5 results:

- Phase 5 unit tests: `5 passed`
- Full repository test suite: `40 passed`
- Scheme A representative max lower-order jump: `6.580e-12`
- Scheme B representative max lower-order jump: `1.014e-11`
- Scheme C representative position-jump max: `3.997e-15`
- Scheme C representative velocity/acceleration/jerk jump max: `3.952e+01`
- Scheme B representative max pre-consensus dispersion: `2.008e+01`
- Random-trial success rate: `A=1.00, B=1.00, C=1.00`
- Random-trial mean lower-order jump max: `A=1.307e-11, B=1.331e-11, C=1.976e+01`

### Deliverables

Implemented deliverables under `phase_5/`:

- `REQ-Phase5-boundary-jump-check.md`
- `__init__.py`
- `blom_boundary_jump_check.py`
- `test_blom_boundary_jump_check.py`
- `examples/`
- `README_phase5.md`
- `results/phase5_boundary_jump_check/`

## Phase 6 Traceability

This section maps the Phase 6 requirements in `phase_6/REQ-Phase6-fd-jacobian-check.md`
to the implemented Jacobian-validation toolbox, tests, and generated artifacts.

### Scope

Implemented Phase 6 object:

- theoretical sparsity masks for the canonical raw BLOM-Analytic coefficient map
- analytic-vs-finite-difference validation for raw `J_c_q`, `J_c_T`, `J_x_q`, and `J_x_T`
- scheme-level locality comparison for final assembled Scheme A / B / C coefficients
- random-trial verification of mask stability, Jacobian error levels, and `||dc/dT||` boundedness statistics

### Requirement Mapping

| Requirement | Implementation | Verification |
| --- | --- | --- |
| Theoretical raw waypoint and duration support masks | `phase_6/blom_fd_jacobian_check.py` -> `theoretical_mask_c_q`, `theoretical_mask_c_T`, `theoretical_mask_x_q`, `theoretical_mask_x_T` | `phase_6/test_blom_fd_jacobian_check.py`, raw sparsity CSV and theory-vs-FD plots under `phase_6/results/phase6_fd_jacobian_check/raw_scheme_C/` |
| Unified vector-output finite-difference Jacobian helper | `phase_6/blom_fd_jacobian_check.py` -> `finite_difference_jacobian` | `phase_6/test_blom_fd_jacobian_check.py`, all Phase 6 FD comparisons |
| Raw canonical BLOM-Analytic coefficient and local-jet maps | `phase_6/blom_fd_jacobian_check.py` -> `raw_local_coefficient_map`, `raw_local_jet_state_map`, `raw_local_jacobians` | `phase_6/test_blom_fd_jacobian_check.py`, `python3 -m phase_6.examples.demo_raw_local_support` |
| Analytic derivative chain for `d x / dT` and `d c / dT` in the canonical local system | `phase_6/blom_fd_jacobian_check.py` -> `_raw_local_analytic_jacobians_local`, `_d_build_H_mid`, `_d_build_g_mid`, `_d_outer_rank_one_hessian`, `_d_outer_linear_term`, `_d_D_matrix`, `_d_Lambda8` | raw analytic-vs-FD error heatmaps and step-size sweeps in `phase_6/results/phase6_fd_jacobian_check/compare/` |
| Scheme A / B / C assembled coefficient maps and nonlocality diagnostics | `phase_6/blom_fd_jacobian_check.py` -> `assembled_scheme_A_map`, `assembled_scheme_B_map`, `assembled_scheme_C_map`, `assembled_scheme_jacobians_fd`, `_scheme_nonlocality_check` | `phase_6/test_blom_fd_jacobian_check.py`, scheme heatmaps and sparsity CSVs under `phase_6/results/phase6_fd_jacobian_check/` |
| High-level single-run, compare-all, and random-trial runners | `phase_6/blom_fd_jacobian_check.py` -> `run_fd_jacobian_check`, `run_compare_all_schemes`, `run_random_trials` | `phase_6/examples/demo_compare_all.py`, random-trial summaries and boxplots |
| Phase 6 example entry points | `phase_6/examples/demo_raw_local_support.py`, `phase_6/examples/demo_scheme_A_nonlocality.py`, `phase_6/examples/demo_scheme_B_nonlocality.py`, `phase_6/examples/demo_compare_all.py` | direct script execution |
| Phase 6 README documenting interfaces, outputs, and interpretation | `phase_6/README_phase6.md` | manual review |

### Verification Summary

Executed for Phase 6 implementation:

- `python3 -m compileall phase_6`
- `python3 -m phase_6.blom_fd_jacobian_check`
- `python3 -m phase_6.examples.demo_raw_local_support`
- `python3 -m phase_6.examples.demo_scheme_A_nonlocality`
- `python3 -m phase_6.examples.demo_scheme_B_nonlocality`
- `python3 -m phase_6.examples.demo_compare_all`
- `python3 -m unittest phase_6.test_blom_fd_jacobian_check`
- `python3 -m unittest discover -s . -p 'test*.py'`

Observed Phase 6 results:

- Phase 6 unit tests: `5 passed`
- Full repository test suite: `45 passed`
- Raw `J_c_q` analytic-vs-FD max abs error: `4.063e-08`
- Raw `J_c_T` analytic-vs-FD max abs error: `1.209e-07`
- Raw `J_c_q` mask pass: `True`
- Raw `J_c_T` mask pass: `True`
- Scheme C final `q` bandwidth: `2`
- Scheme A final `q` bandwidth: `5`
- Scheme B final `q` bandwidth: `3`
- Scheme C final outside-band `q` max abs: `0.0`
- Scheme A final outside-band `q` max abs: `1.317`
- Scheme B final outside-band `q` max abs: `8.843`
- Random-trial raw `q` mask pass rate: `1.00`
- Random-trial raw `T` mask pass rate: `1.00`
- Random-trial mean `q` bandwidth: `A=7.0, B=3.0, C=2.0`
- Random-trial max `||dc/dT||`: `1.765`

### Deliverables

Implemented deliverables under `phase_6/`:

- `REQ-Phase6-fd-jacobian-check.md`
- `__init__.py`
- `blom_fd_jacobian_check.py`
- `test_blom_fd_jacobian_check.py`
- `examples/`
- `README_phase6.md`
- `results/phase6_fd_jacobian_check/`

## Phase 7 Traceability

This section maps the Phase 7 requirements in `phase_7/REQ-Phase7-convergence-validation.md`
to the implemented convergence-validation toolbox, tests, and generated artifacts.

### Scope

Implemented Phase 7 object:

- exact Phase 1 MINCO reference construction and centered kernel extraction
- idealized truncated BLOM-k built by truncating `G(T)=A(T)^{-1}S_q`
- actual BLOM-k comparison using assembly Schemes A / B / C
- coefficient-error, matching-error, and cost-gap sweeps over `k`
- log-linear fitting of `log ||c^(k) - c*||`
- randomized batch validation and automatic interpretation summaries

### Requirement Mapping

| Requirement | Implementation | Verification |
| --- | --- | --- |
| Exact MINCO reference plus centered kernel `G(T)=A(T)^{-1}S_q` | `phase_7/blom_convergence_vs_k.py` -> `compute_minco_reference`, `compute_total_snap_cost` | `phase_7/test_blom_convergence_vs_k.py`, `python3 -m phase_7.blom_convergence_vs_k` |
| k-grid generation for `k = 2, 4, ...` with terminal inclusion | `phase_7/blom_convergence_vs_k.py` -> `make_k_grid` | `phase_7/test_blom_convergence_vs_k.py` |
| Idealized truncated BLOM-k from blockwise kernel truncation | `phase_7/blom_convergence_vs_k.py` -> `compute_ideal_truncated_blom_k`, `_ideal_mask` | `phase_7/test_blom_convergence_vs_k.py`, `phase_7/results/phase7_convergence_vs_k/ideal/summary_ideal.json` |
| Actual BLOM-k for Schemes A / B / C | `phase_7/blom_convergence_vs_k.py` -> `compute_actual_blom_k` | `phase_7/test_blom_convergence_vs_k.py`, `phase_7/results/phase7_convergence_vs_k/scheme_A/summary_scheme_A.json`, `phase_7/results/phase7_convergence_vs_k/scheme_B/summary_scheme_B.json`, `phase_7/results/phase7_convergence_vs_k/scheme_C/summary_scheme_C.json` |
| Unified coefficient-error, matching-error, and cost-gap metrics | `phase_7/blom_convergence_vs_k.py` -> `compute_convergence_errors` | `phase_7/results/phase7_convergence_vs_k/compare/convergence_errors_by_k.csv`, `RESULT_PHASE7.md` |
| Least-squares fit of `log(error)` against `k` | `phase_7/blom_convergence_vs_k.py` -> `fit_log_error_vs_k` | `phase_7/results/phase7_convergence_vs_k/compare/logfit_summary.csv`, `phase_7/results/phase7_convergence_vs_k/compare/log_coef_error_fit_all.png` |
| Single-case Phase 7 sweep runner with mandatory plots and summary files | `phase_7/blom_convergence_vs_k.py` -> `run_convergence_vs_k`, `_write_interpretation_summary` | `python3 -m phase_7.blom_convergence_vs_k`, compare artifacts under `phase_7/results/phase7_convergence_vs_k/compare/` |
| Randomized multi-trial convergence statistics | `phase_7/blom_convergence_vs_k.py` -> `run_random_trials` | `python3 -m phase_7.examples.demo_random_trials`, `phase_7/results/phase7_convergence_vs_k/compare/random_trials_summary.csv`, `phase_7/results/phase7_convergence_vs_k/compare/random_trials_aggregate.json` |
| Example entry points for single-case, scheme comparison, ideal-vs-actual, and random-trial workflows | `phase_7/examples/demo_single_case.py`, `phase_7/examples/demo_compare_schemes.py`, `phase_7/examples/demo_ideal_vs_actual.py`, `phase_7/examples/demo_random_trials.py` | direct script execution |
| Phase 7 README documenting scope, outputs, and the Scheme A caveat | `phase_7/README_phase7.md` | manual review |

### Verification Summary

Executed for Phase 7 implementation:

- `python3 -m compileall phase_7`
- `python3 -m phase_7.blom_convergence_vs_k`
- `python3 -m phase_7.examples.demo_single_case`
- `python3 -m phase_7.examples.demo_compare_schemes`
- `python3 -m phase_7.examples.demo_ideal_vs_actual`
- `python3 -m phase_7.examples.demo_random_trials`
- `python3 -m unittest phase_7.test_blom_convergence_vs_k`
- `python3 -m unittest discover -s . -p 'test*.py'`

Observed Phase 7 results:

- Phase 7 unit tests: `5 passed`
- Full repository test suite: `50 passed`
- Ideal-truncation final coefficient error at `k = 8`: `2.169e-13`
- Scheme A final coefficient error at `k = 8`: `7.400e-12`
- Scheme B final coefficient error at `k = 8`: `1.417e+02`
- Scheme C final coefficient error at `k = 8`: `3.469e+00`
- Ideal log-fit slope: `-4.686`, `R^2 = 0.686`
- Scheme A log-fit slope: `-3.553e-15`, `R^2 = 1.000`
- Scheme B log-fit slope: `-2.070e-01`, `R^2 = 0.922`
- Scheme C log-fit slope: `-3.160e-01`, `R^2 = 0.996`
- Random-trial success rates: `ideal=1.00, A=1.00, B=1.00, C=1.00`

### Deliverables

Implemented deliverables under `phase_7/`:

- `REQ-Phase7-convergence-validation.md`
- `__init__.py`
- `blom_convergence_vs_k.py`
- `test_blom_convergence_vs_k.py`
- `examples/`
- `README_phase7.md`
- `results/phase7_convergence_vs_k/`

## Phase 7 Extra Experiments Traceability

This section maps the extra Phase 7 experiment requirements in
`phase_7/REQ-Phase7-extra-experiments-v2.md` to the implemented follow-up
experiment toolbox, tests, and generated artifacts.

### Scope

Implemented Phase 7 extra object:

- large-`M` convergence sweep over ideal truncation and actual BLOM families
- time-regime split between uniform and bounded-nonuniform duration regimes
- full-vs-interior coefficient and matching error separation
- light corrected Scheme C variants with jump and locality tradeoff checks
- cross-experiment overview table and automated interpretation summary

### Requirement Mapping

| Requirement | Implementation | Verification |
| --- | --- | --- |
| Standalone extra experiment module with unified runners | `phase_7/blom_phase7_extra_experiments.py` -> `run_large_M_sweep`, `run_time_regime_split`, `run_interior_vs_full`, `run_schemeC_light_assembly`, `run_phase7_extra_experiments` | `phase_7/test_blom_phase7_extra_experiments.py`, `python3 -m phase_7.blom_phase7_extra_experiments` |
| Large-`M` sweep with saved figures, CSV, and JSON | `phase_7/blom_phase7_extra_experiments.py` -> `run_large_M_sweep` | `phase_7/examples/demo_exp1_large_M.py`, artifacts under `phase_7/results/phase7_extra_experiments/exp1_large_M_sweep/` |
| Uniform-time vs bounded-nonuniform split with slope, matching, and cost-gap summaries | `phase_7/blom_phase7_extra_experiments.py` -> `run_time_regime_split` | `phase_7/examples/demo_exp2_time_regime.py`, artifacts under `phase_7/results/phase7_extra_experiments/exp2_time_regime_split/` |
| Interior-only vs full-error separation with boundary-trim support | `phase_7/blom_phase7_extra_experiments.py` -> `compute_full_and_interior_errors`, `run_interior_vs_full` | `phase_7/examples/demo_exp3_interior_vs_full.py`, `summary_interior_vs_full.csv`, heatmap and boxplot artifacts |
| Light corrected Scheme C variants `C1` and `C2` | `phase_7/blom_phase7_extra_experiments.py` -> `assemble_scheme_C_light`, `run_schemeC_light_assembly` | `phase_7/examples/demo_exp4_schemeC_correction.py`, `summary_schemeC_correction.csv`, correction comparison figures |
| Locality-width comparison for raw and corrected Scheme C | `phase_7/blom_phase7_extra_experiments.py` -> `_scheme_c_q_bandwidth` | `schemeC_raw_vs_corrected_locality_tradeoff.png`, `RESULT_PHASE7_extra.md` |
| Global cross-experiment overview and automated research-direction summary | `phase_7/blom_phase7_extra_experiments.py` -> `_build_overview_rows`, `_write_interpretation_summary`, `run_phase7_extra_experiments` | `compare_summary/final_experiment_overview.csv`, `compare_summary/phase7_extra_interpretation_summary.md`, `demo_all_phase7_extra.py` |
| Extra experiment demos and regression tests | `phase_7/examples/demo_exp1_large_M.py`, `phase_7/examples/demo_exp2_time_regime.py`, `phase_7/examples/demo_exp3_interior_vs_full.py`, `phase_7/examples/demo_exp4_schemeC_correction.py`, `phase_7/examples/demo_all_phase7_extra.py`, `phase_7/test_blom_phase7_extra_experiments.py` | direct script execution and unit-test run |

### Verification Summary

Executed for Phase 7 extra experiments:

- `python3 -m compileall phase_7`
- `python3 -m phase_7.blom_phase7_extra_experiments`
- `python3 -m phase_7.examples.demo_exp1_large_M`
- `python3 -m phase_7.examples.demo_exp2_time_regime`
- `python3 -m phase_7.examples.demo_exp3_interior_vs_full`
- `python3 -m phase_7.examples.demo_exp4_schemeC_correction`
- `python3 -m phase_7.examples.demo_all_phase7_extra`
- `python3 -m unittest phase_7.test_blom_phase7_extra_experiments`
- `python3 -m unittest phase_7.test_blom_convergence_vs_k`
- `python3 -m unittest discover -s . -p 'test*.py'`

Observed extra-experiment results:

- Phase 7 extra unit tests: `5 passed`
- Phase 7 original unit tests re-run: `5 passed`
- Full repository test suite: `55 passed`
- Large-`M` sweep median final Scheme C coefficient error: `6.263e-02`
- Time-regime split median final Scheme C coefficient error: `1.713e-02`
- Interior-only final Scheme C coefficient error at `k = 8`: `1.564`
- Scheme C final coefficient error at `k = 8`: `3.469`
- Corrected `C1` final coefficient error at `k = 8`: `5.984e+01`
- Corrected `C2` final coefficient error at `k = 8`: `4.990e+01`
- Raw Scheme C final `q` bandwidth: `5`
- Corrected `C1` final `q` bandwidth: `6`
- Corrected `C2` final `q` bandwidth: `6`

### Deliverables

Implemented extra-experiment deliverables under `phase_7/`:

- `REQ-Phase7-extra-experiments-v2.md`
- `blom_phase7_extra_experiments.py`
- `test_blom_phase7_extra_experiments.py`
- `examples/demo_exp1_large_M.py`
- `examples/demo_exp2_time_regime.py`
- `examples/demo_exp3_interior_vs_full.py`
- `examples/demo_exp4_schemeC_correction.py`
- `examples/demo_all_phase7_extra.py`
- `results/phase7_extra_experiments/`

## Phase 8 Traceability

This section maps the Phase 8 requirements in
`phase_8/REQ-Phase8-interior-matching-validation.md` to the implemented
interior-first matching toolbox, tests, and generated artifacts.

### Scope

Implemented Phase 8 object:

- interior-first matching validation for raw Scheme C versus ideal truncated BLOM-k
- explicit interior / boundary block projectors with configurable radius `r(k)`
- boundary-gap decomposition with projector statistics and radius sensitivity
- augmented local window solve with artificial high-order boundary traces for reference-window experiments
- uniform-time versus bounded-nonuniform regime split
- unified suite runner, overview table, and automated interpretation summary

### Requirement Mapping

| Requirement | Implementation | Verification |
| --- | --- | --- |
| Shared Phase 8 utilities for `r(k)`, interior/boundary sets, block errors, and result serialization | `phase_8/phase8_common.py` -> `boundary_radius`, `make_interior_sets`, `compute_matching_error_blocks`, `project_interior_boundary_errors`, `ensure_results_dirs`, `save_csv`, `save_json` | `phase_8/test_blom_interior_matching_check.py`, `phase_8/test_blom_boundary_gap_decomposition.py` |
| Augmented local window solve with artificial boundary traces `gamma` and reference-window extraction | `phase_8/phase8_common.py` -> `build_augmented_local_constraints`, `solve_augmented_local_window`, `extract_reference_gamma`, `compute_reference_window_coefficients`, `compute_boundary_response_matrix` | `phase_8/test_blom_boundary_gap_decomposition.py`, `boundary_gap_decomposition.csv` |
| Script 1 interior-first matching runner and mandatory artifacts | `phase_8/blom_interior_matching_check.py` -> `make_interior_sets`, `compute_matching_error_blocks`, `summarize_interior_matching`, `run_interior_matching_check` | `python3 -m phase_8.blom_interior_matching_check`, `phase_8/results/phase8_validation/interior_matching/` |
| Script 2 boundary-gap decomposition runner and radius sensitivity analysis | `phase_8/blom_boundary_gap_decomposition.py` -> `project_interior_boundary_errors`, `run_boundary_gap_decomposition` | `python3 -m phase_8.blom_boundary_gap_decomposition`, `phase_8/results/phase8_validation/boundary_gap/` |
| Script 3 uniform-time versus bounded-nonuniform regime split | `phase_8/blom_uniform_vs_nonuniform_interior.py` -> `sample_uniform_time`, `sample_bounded_nonuniform_time`, `run_uniform_vs_nonuniform_interior` | `python3 -m phase_8.blom_uniform_vs_nonuniform_interior`, `phase_8/results/phase8_validation/uniform_vs_nonuniform/` |
| Unified Phase 8 suite runner with overview CSV and interpretation summary | `phase_8/blom_phase8_validation_suite.py` -> `run_phase8_validation_suite` | `python3 -m phase_8.blom_phase8_validation_suite`, `phase_8/results/phase8_validation/compare/` |
| Phase 8 examples for each experiment family and the full suite | `phase_8/examples/demo_interior_matching.py`, `phase_8/examples/demo_boundary_gap.py`, `phase_8/examples/demo_uniform_vs_nonuniform.py`, `phase_8/examples/demo_phase8_suite.py` | direct script execution |
| Phase 8 README documenting interfaces, outputs, and interpretation target | `phase_8/README_phase8.md` | manual review |

### Verification Summary

Executed for Phase 8 implementation:

- `python3 -m compileall phase_8`
- `python3 -m phase_8.blom_interior_matching_check`
- `python3 -m phase_8.blom_boundary_gap_decomposition`
- `python3 -m phase_8.blom_uniform_vs_nonuniform_interior`
- `python3 -m phase_8.blom_phase8_validation_suite`
- `python3 -m phase_8.examples.demo_interior_matching`
- `python3 -m phase_8.examples.demo_boundary_gap`
- `python3 -m phase_8.examples.demo_uniform_vs_nonuniform`
- `python3 -m phase_8.examples.demo_phase8_suite`
- `python3 -m unittest phase_8.test_blom_interior_matching_check phase_8.test_blom_boundary_gap_decomposition phase_8.test_blom_uniform_vs_nonuniform_interior`
- `python3 -m unittest discover -s . -p 'test*.py'`

Observed Phase 8 results:

- Phase 8 unit tests: `8 passed`
- Full repository test suite: `63 passed`
- Single-case interior matching full slope: `-0.446569`, `R^2 = 0.990902`
- Single-case interior matching slope: `-5.265687`, `R^2 = 0.761420`
- Mean boundary energy ratio: `0.898402`
- Mean raw-vs-reference window gap: `4.042120e-14`
- Mean reference-window-vs-ideal gap: `1.528579e+00`
- Uniform-time mean interior slope: `-5.385500`
- Bounded-nonuniform mean interior slope: `-5.309442`

### Deliverables

Implemented deliverables under `phase_8/`:

- `REQ-Phase8-interior-matching-validation.md`
- `__init__.py`
- `phase8_common.py`
- `blom_interior_matching_check.py`
- `blom_boundary_gap_decomposition.py`
- `blom_uniform_vs_nonuniform_interior.py`
- `blom_phase8_validation_suite.py`
- `test_blom_interior_matching_check.py`
- `test_blom_boundary_gap_decomposition.py`
- `test_blom_uniform_vs_nonuniform_interior.py`
- `examples/`
- `README_phase8.md`
- `results/phase8_validation/`

## Phase 9 Traceability

## Phase 8 Supplementary Traceability

This section maps the supplementary Phase 8 requirements in
`phase_8/REQ-Phase8-supplementary-experiments.md` to the implemented
large-scale robustness checks, tests, and generated artifacts.

### Scope

Implemented supplementary Phase 8 object:

- larger-`M` sweep with explicit nonempty-interior accounting
- higher-trial robustness checks under uniform-time and bounded-nonuniform regimes
- radius-mode sensitivity analysis for the interior/boundary split
- raw-vs-reference-window sanity checks with inherited-`gamma` perturbation response
- explicit empty-interior false-positive detection
- raw-to-reference vs reference-to-ideal two-bridge comparison
- unified supplementary suite runner with overview CSV and interpretation summary

### Requirement Mapping

| Requirement | Implementation | Verification |
| --- | --- | --- |
| Shared supplementary utilities for result directories, radius-mode aliases, case generation, fast reference-window extraction, and cross-experiment aggregation | `phase_8/phase8_supplementary_common.py` -> `ensure_supplementary_results_dirs`, `normalize_radius_mode`, `build_case`, `compute_reference_window_coefficients_fast`, `compute_case_metrics`, `aggregate_numeric_records` | `phase_8/test_phase8_supplementary_suite.py` |
| Experiment 1 larger-`M` sweep with interior count, boundary ratio, and segmentwise heatmaps | `phase_8/exp1_large_M_sweep.py` -> `run_large_M_sweep` | `python3 -m phase_8.blom_phase8_supplementary_suite`, `phase_8/results/phase8_supplementary/exp1_large_M/` |
| Experiment 2 more-trials robustness with slope boxplots, boundary-ratio boxplots, and mean-curve bands | `phase_8/exp2_more_trials.py` -> `run_more_trials` | `python3 -m phase_8.blom_phase8_supplementary_suite`, `phase_8/results/phase8_supplementary/exp2_more_trials/` |
| Experiment 3 radius sensitivity across `r(k)=ceil(k/2)`, `r(k)=k`, and `r(k)=ceil(3k/2)` | `phase_8/exp3_radius_sensitivity.py` -> `run_radius_sensitivity` | `python3 -m phase_8.blom_phase8_supplementary_suite`, `phase_8/results/phase8_supplementary/exp3_radius_sensitivity/` |
| Experiment 4 raw-vs-reference sanity checks with large-`M`, regime split, and perturbed inherited-`gamma` response | `phase_8/exp4_raw_vs_reference_sanity.py` -> `run_raw_vs_reference_sanity` | `python3 -m phase_8.blom_phase8_supplementary_suite`, `phase_8/results/phase8_supplementary/exp4_raw_vs_reference_sanity/` |
| Experiment 5 explicit empty-interior risk detection and annotation outputs | `phase_8/exp5_empty_interior_risk.py` -> `run_empty_interior_risk` | `python3 -m phase_8.blom_phase8_supplementary_suite`, `phase_8/results/phase8_supplementary/exp5_empty_interior_risk/` |
| Experiment 6 raw-to-reference vs reference-to-ideal bridge-gap comparison | `phase_8/exp6_two_bridge_gap_compare.py` -> `run_two_bridge_gap_compare` | `python3 -m phase_8.blom_phase8_supplementary_suite`, `phase_8/results/phase8_supplementary/exp6_two_bridge_gaps/` |
| Unified supplementary suite with total overview CSV and mandatory high-level summary answers | `phase_8/blom_phase8_supplementary_suite.py` -> `run_phase8_supplementary_suite`, `_write_suite_summary` | `phase_8/results/phase8_supplementary/summary/phase8_supplementary_overview.csv`, `phase_8/results/phase8_supplementary/summary/phase8_supplementary_summary.md` |
| Supplementary tests covering radius aliases, raw/reference sanity smoke run, empty-interior marking, and suite smoke run | `phase_8/test_phase8_supplementary_suite.py` | `python3 -m unittest phase_8.test_phase8_supplementary_suite`, `python3 -m unittest discover -s phase_8 -p 'test*.py'` |
| Supplementary examples and README updates | `phase_8/examples/demo_exp1_large_M.py`, `phase_8/examples/demo_exp2_more_trials.py`, `phase_8/examples/demo_exp3_radius.py`, `phase_8/examples/demo_exp4_raw_ref.py`, `phase_8/examples/demo_exp5_empty_interior.py`, `phase_8/examples/demo_exp6_two_bridge.py`, `phase_8/examples/demo_phase8_supplementary_suite.py`, `phase_8/README_phase8.md` | direct script execution and manual review |

### Verification Summary

Executed for supplementary Phase 8 implementation:

- `python3 -m compileall phase_8`
- `python3 -m phase_8.blom_phase8_supplementary_suite`
- `python3 -m unittest discover -s phase_8 -p 'test*.py'`

Observed supplementary Phase 8 results:

- Phase 8 full test subset: `12 passed`
- Larger-`M` sweep nonempty-interior cases: `22 / 24`
- Larger-`M` sweep mean boundary energy ratio: `0.424942`
- More-trials mean interior slope: `-4.147938`
- More-trials mean full slope: `-0.369023`
- Share of trials with interior slope more negative than full slope: `1.0`
- Raw-vs-reference mean gap: `2.187265e-13`
- Raw-vs-reference max gap: `6.702577e-12`
- Empty-interior cases flagged: `13 / 87`
- Mean two-bridge gap ratio: `3.541086e+12`
- Recommended next bridge from the suite summary: `reference-window / ideal-truncation consistency`

### Deliverables

Implemented supplementary deliverables under `phase_8/`:

- `REQ-Phase8-supplementary-experiments.md`
- `phase8_supplementary_common.py`
- `exp1_large_M_sweep.py`
- `exp2_more_trials.py`
- `exp3_radius_sensitivity.py`
- `exp4_raw_vs_reference_sanity.py`
- `exp5_empty_interior_risk.py`
- `exp6_two_bridge_gap_compare.py`
- `blom_phase8_supplementary_suite.py`
- `test_phase8_supplementary_suite.py`
- `examples/demo_exp1_large_M.py`
- `examples/demo_exp2_more_trials.py`
- `examples/demo_exp3_radius.py`
- `examples/demo_exp4_raw_ref.py`
- `examples/demo_exp5_empty_interior.py`
- `examples/demo_exp6_two_bridge.py`
- `examples/demo_phase8_supplementary_suite.py`
- `results/phase8_supplementary/`

## Phase 9 Traceability

This section maps the Phase 9 requirements in
`phase_9/REQ-Phase9-minimal-differentiable-loop.md` to the implemented
minimal differentiable-loop toolbox, tests, and generated artifacts.

### Scope

Implemented Phase 9 object:

- raw Scheme C coefficient and Jacobian interfaces in the canonical case `s=4, k=2`
- exact control / time / smooth obstacle objective terms with block gradients
- dense checker and exact block-banded backward accumulation
- finite-difference gradient validation for the total objective and each objective term
- minimal projected space-time gradient-descent demo
- unified Phase 9 suite with summary tables, plots, and interpretation output

### Requirement Mapping

| Requirement | Implementation | Verification |
| --- | --- | --- |
| Raw Scheme C coefficient interface returning block and flattened coefficients | `phase_9/blom_backward_diff.py` -> `compute_raw_schemeC_coeffs` | `phase_9/test_blom_backward_diff.py` |
| Raw Scheme C dense Jacobian and block-Jacobian interfaces | `phase_9/blom_backward_diff.py` -> `compute_raw_schemeC_jacobians`, `_local_schemeC_sensitivity` | `phase_9/test_blom_backward_diff.py`, `phase9_jacobian_sparsity.png` |
| Control cost, time penalty, and smooth obstacle penalty with exact block gradients | `phase_9/blom_backward_diff.py` -> `control_cost`, `time_penalty`, `soft_obstacle_penalty`, `obstacle_penalty_and_grad` | `phase_9/test_blom_backward_diff.py`, `phase9_gradcheck_summary.csv` |
| Total minimal objective assembly and exact block gradients `g_c`, `g_T` | `phase_9/blom_backward_diff.py` -> `evaluate_minimal_objective` | `phase_9/test_blom_backward_diff.py`, `phase9_gradcheck_summary.md` |
| Dense reverse differentiation and exact banded accumulation | `phase_9/blom_backward_diff.py` -> `backward_diff_dense`, `backward_diff_banded` | `phase_9/test_blom_backward_diff.py`, `phase9_dense_vs_banded_grad.png` |
| Gradient-check runner with FD comparison, plots, CSV, JSON, and summary | `phase_9/blom_backward_diff.py` -> `run_phase9_gradcheck` | `python3 -m phase_9.blom_backward_diff`, `phase_9/results/phase9_validation/gradcheck/` |
| Minimal projected gradient-descent demo for raw Scheme C | `phase_9/blom_space_time_opt_demo.py` -> `gradient_descent_step`, `run_minimal_optimization_demo` | `python3 -m phase_9.blom_space_time_opt_demo`, `phase_9/results/phase9_validation/optimization_demo/` |
| Unified Phase 9 validation suite with overview CSV and interpretation summary | `phase_9/blom_phase9_validation_suite.py` -> `run_phase9_validation_suite` | `python3 -m phase_9.blom_phase9_validation_suite`, `phase_9/results/phase9_validation/compare/` |
| Phase 9 tests for FD consistency, dense-vs-banded equality, per-term gradients, and demo smoke runs | `phase_9/test_blom_backward_diff.py` | `python3 -m unittest phase_9.test_blom_backward_diff` |
| Phase 9 examples and README | `phase_9/examples/`, `phase_9/README_phase9.md` | direct script execution and manual review |

### Verification Summary

Executed for Phase 9 implementation:

- `python3 -m compileall phase_9`
- `python3 -m phase_9.blom_backward_diff`
- `python3 -m phase_9.blom_space_time_opt_demo`
- `python3 -m phase_9.blom_phase9_validation_suite`
- `python3 -m phase_9.examples.demo_backward_diff`
- `python3 -m phase_9.examples.demo_gradcheck`
- `python3 -m phase_9.examples.demo_minimal_opt`
- `python3 -m phase_9.examples.demo_phase9_suite`
- `python3 -m unittest phase_9.test_blom_backward_diff`
- `python3 -m unittest discover -s . -p 'test*.py'`

Observed Phase 9 results:

- Phase 9 unit tests: `8 passed`
- Full repository test suite: `71 passed`
- Total dense-vs-FD gradient max abs error: `3.755835e-05`
- Total banded-vs-FD gradient max abs error: `3.755835e-05`
- Dense-vs-banded max abs gap: `1.818989e-12`
- Control-term partial gradient max abs error: `1.485118e-05`
- Time-term partial gradient max abs error: `2.631780e-11`
- Obstacle-term partial gradient max abs error: `4.200990e-09`
- Minimal optimization objective drop: `2.688224e+03`
- Minimal optimization gradient-norm drop proxy: `1.243217e+04`
- Any duration touching the lower bound: `False`

### Deliverables

Implemented deliverables under `phase_9/`:

- `REQ-Phase9-minimal-differentiable-loop.md`
- `__init__.py`
- `blom_backward_diff.py`
- `blom_space_time_opt_demo.py`
- `blom_phase9_validation_suite.py`
- `test_blom_backward_diff.py`
- `examples/`
- `README_phase9.md`
- `results/phase9_validation/`

## Phase 10 Traceability

This section maps the Phase 10 requirements in `phase_10/REQ-Phase10-full-space-time-framework.md`
to the implemented full objective / optimizer / benchmark framework, result artifacts, and verification assets.

### Scope

Implemented Phase 10 object:

- raw Scheme C full objective layer on top of the canonical `s = 4` local family
- general-even-`k` public interfaces for coefficient and Jacobian access
- dense checker and sparse backward dual implementation
- `xi = (q_bar, tau)` reparameterized optimizer
- benchmark, baseline comparison, and ablation framework
- unified suite aggregating backward, optimizer, and benchmark evidence

### Requirement Mapping

| Requirement | Implementation | Verification |
| --- | --- | --- |
| General-`k` raw Scheme C coefficient API returning block and flattened coefficients | `phase_10/blom_full_backward_diff.py` -> `compute_raw_schemeC_coeffs_general_k` | `phase_10/test_blom_full_backward_diff.py`, `python3 -m phase_10.blom_full_backward_diff` |
| General-`k` Jacobian API with dense matrices, block dictionaries, and support-set metadata | `phase_10/blom_full_backward_diff.py` -> `compute_raw_schemeC_jacobians_general_k`, `_local_schemeC_sensitivity` | `phase_10/test_blom_full_backward_diff.py`, `python3 -m phase_10.blom_full_backward_diff`, sparsity/support outputs under `phase_10/results/phase10_framework/backward_diff/` |
| Time reparameterization `tau -> T` and chain-rule support for `dT/dtau` | `phase_10/blom_full_backward_diff.py` -> `tau_to_T`, `dT_dtau`, `T_to_tau`, `q_tau_to_xi`, `xi_to_qbar_tau` | `phase_10/test_blom_full_backward_diff.py`, `phase_10/test_blom_space_time_optimizer.py` |
| Full objective layer `K_ctrl + lambda_T K_time + lambda_obs K_obs + lambda_dyn K_dyn + lambda_bc K_bc + lambda_reg K_reg` with value/gradient outputs | `phase_10/blom_full_backward_diff.py` -> `control_cost_full`, `time_penalty_full`, `obstacle_penalty_full`, `dynamics_penalty_full`, `boundary_penalty_full`, `regularization_penalty_full`, `evaluate_full_objective`, `evaluate_full_objective_from_coeffs` | `phase_10/test_blom_full_backward_diff.py`, `python3 -m phase_10.blom_full_backward_diff`, `python3 -m phase_10.blom_space_time_optimizer` |
| Dense backward, sparse backward, and reparameterized backward under `xi=(q_bar,tau)` | `phase_10/blom_full_backward_diff.py` -> `full_backward_diff_dense`, `full_backward_diff_sparse`, `full_backward_diff_reparam`, `finite_difference_gradient` | `phase_10/test_blom_full_backward_diff.py`, `python3 -m phase_10.blom_full_backward_diff`, `phase_10/results/phase10_framework/backward_diff/phase10_backward_summary.json` |
| Optimizer-level framework in `xi` space with step rules and optimization history export | `phase_10/blom_space_time_optimizer.py` -> `optimizer_step`, `run_space_time_optimization` | `phase_10/test_blom_space_time_optimizer.py`, `python3 -m phase_10.blom_space_time_optimizer`, optimizer artifacts under `phase_10/results/phase10_framework/optimizer/` |
| Framework-level benchmark, baseline comparison, `M`-sweep, `k`-sweep, and ablation study | `phase_10/blom_benchmark_suite.py` -> `run_baseline_raw_schemeC`, `run_baseline_minco`, `run_baseline_schemeA`, `run_baseline_heuristic`, `run_benchmark_suite` | `phase_10/test_blom_benchmark_suite.py`, `python3 -m phase_10.blom_benchmark_suite`, benchmark artifacts under `phase_10/results/phase10_framework/summary/` |
| Unified Phase 10 suite integrating backward, optimizer, benchmark, and top-level interpretation summary | `phase_10/blom_phase10_framework_suite.py` -> `run_phase10_framework_suite` | `python3 -m phase_10.blom_phase10_framework_suite`, `phase_10/results/phase10_framework/summary/phase10_suite_summary.json` |
| Example entry points for backward, optimizer, `k`-sweep, baseline compare, ablation, and full suite replay | `phase_10/examples/demo_full_backward_diff.py`, `phase_10/examples/demo_optimizer.py`, `phase_10/examples/demo_k_sweep.py`, `phase_10/examples/demo_baseline_compare.py`, `phase_10/examples/demo_ablation.py`, `phase_10/examples/demo_phase10_suite.py` | direct execution of each example entry point |
| Phase 10 README describing scope, main entry points, and run commands | `phase_10/README_phase10.md` | manual review |

### Verification Summary

Executed for Phase 10 implementation:

- `python3 -m compileall phase_10`
- `python3 -m phase_10.blom_full_backward_diff`
- `python3 -m phase_10.blom_space_time_optimizer`
- `python3 -m phase_10.blom_benchmark_suite`
- `python3 -m phase_10.blom_phase10_framework_suite`
- `python3 -m phase_10.examples.demo_full_backward_diff`
- `python3 -m phase_10.examples.demo_optimizer`
- `python3 -m phase_10.examples.demo_k_sweep`
- `python3 -m phase_10.examples.demo_baseline_compare`
- `python3 -m phase_10.examples.demo_ablation`
- `python3 -m phase_10.examples.demo_phase10_suite`
- `python3 -m unittest discover -s phase_10 -p 'test*.py'`
- `python3 -m unittest discover -s . -p 'test*.py'`

Observed Phase 10 results:

- Phase 10 unit tests: `5 passed`
- Full repository test suite: `80 passed`
- Dense-vs-sparse backward max abs gap: `9.094947e-13`
- Reparameterized gradient vs FD max abs gap: `1.326091e-05` in standalone backward run
- Standalone optimizer objective drop: `2.702673e+03`
- Standalone optimizer final objective: `1.851143e+01`
- Framework-suite optimizer objective drop: `2.710814e+03`
- Representative-case `k`-sweep support-width means: `3.0, 4.0, 5.0, 5.6667` for `k = 2, 4, 6, 8`
- Representative-case benchmark best objective baseline: `raw_schemeC`
- Fastest representative-case benchmark baseline: `heuristic`

### Deliverables

Implemented deliverables under `phase_10/`:

- `REQ-Phase10-full-space-time-framework.md`
- `__init__.py`
- `blom_full_backward_diff.py`
- `blom_space_time_optimizer.py`
- `blom_benchmark_suite.py`
- `blom_phase10_framework_suite.py`
- `test_blom_full_backward_diff.py`
- `test_blom_space_time_optimizer.py`
- `test_blom_benchmark_suite.py`
- `examples/`
- `README_phase10.md`
- `results/phase10_framework/`
