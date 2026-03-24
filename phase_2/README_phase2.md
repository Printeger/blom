# Phase 2 Validation Toolbox

This directory implements the Phase 2 numerical evidence pipeline for BLOM necessity theory.

## Purpose

Phase 2 does not solve BLOM itself. Instead, it validates the structural claim that:

- the global Phase 1 MINCO system matrix `A(T)` is banded,
- but the exact waypoint Jacobian `J_q = A(T)^{-1} S_q` is generally not strictly banded,
- so the exact coefficient map does not admit uniform finite local support under the sparse waypoint parameterization.

## Main Files

- `phase2_validation.py`: selector construction, exact Jacobian, finite difference, influence analysis, batch runners
- `phase2_plotting.py`: sparsity plots, heatmaps, profile plots, scaling plots
- `test_phase2_validation.py`: numerical and plotting smoke tests
- `examples/`: demo scripts for the required Phase 2 case families
- `results/`: saved figures, tables, and logs

## Core Interfaces

```python
build_waypoint_selector(M_seg: int, s: int = 4) -> np.ndarray
compute_exact_jacobian_q(M_mat: np.ndarray, S_q: np.ndarray) -> np.ndarray
finite_difference_jacobian_q(q, T, zeta_start, zeta_end, eps=1e-6) -> np.ndarray
compare_exact_vs_fd_jacobian(J_exact, J_fd) -> dict
compute_segmentwise_influence_norms(J_q, M_seg, block_size=8) -> np.ndarray
compute_waypoint_influence_profile(J_q, M_seg, target_waypoint_idx, block_size=8) -> dict
estimate_effective_bandwidth(J_q, M_seg, tol, block_size=8) -> dict
run_phase2_validation_suite(...) -> dict
run_scaling_experiment(...) -> dict
```

## Generated Artifacts

The validation suite saves the required artifact types:

- `A_sparsity_M{M}_case{case_name}.png`
- `Jq_heatmap_M{M}_case{case_name}.png`
- `block_influence_M{M}_case{case_name}.png`
- `influence_profile_q{j}_M{M}_case{case_name}.png`
- `jacobian_fd_compare_M{M}_case{case_name}.png`
- `effective_bandwidth_vs_M_case{case_name}.png`

It also writes:

- `jacobian_error_summary.csv`
- `effective_bandwidth_summary.csv`
- `far_field_influence_summary.csv`
- JSON logs under `results/logs/`

## Minimal Example

```python
from phase_2.phase2_validation import make_uniform_time_case, run_phase2_validation_suite

case = make_uniform_time_case(8, case_name="uniform_time")
summary = run_phase2_validation_suite(
    case.q,
    case.T,
    case.zeta_start,
    case.zeta_end,
    case_name=case.case_name,
)
print(summary["jacobian_comparison"])
```

## Running Demos

```bash
python3 -m phase_2.examples.demo_phase2_uniform_time
python3 -m phase_2.examples.demo_phase2_nonuniform_time
python3 -m phase_2.examples.demo_phase2_scaling_M
python3 -m phase_2.examples.demo_phase2_fd_check
```

## Running Tests

```bash
python3 -m unittest phase_2.test_phase2_validation
```

## Notes

- Phase 2 reuses `phase_1/minco_scalar_baseline.py` directly and does not modify the Phase 1 mother problem.
- Random cases use deterministic seeds.
- Plotting uses the non-interactive `Agg` backend so results can be generated in headless environments.

