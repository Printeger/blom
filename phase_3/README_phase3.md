# Phase 3 BLOM-Strict Validation

This directory implements the Phase 3 numerical validation suite for the
**BLOM-Strict local variational problem**.

## What Phase 3 Proves

Phase 3 does not assemble a full global BLOM trajectory. Its goal is narrower
and more foundational:

- define the local BLOM-Strict problem as a legitimate equality-constrained QP,
- show that the admissible set is nonempty,
- show that the local minimizer is unique,
- verify natural boundary conditions at artificial window boundaries,
- verify automatic higher continuity inside the local window,
- verify continuity of the local solution under small perturbations.

In other words, Phase 3 proves **local well-posedness**, not global assembly.

## Main Files

- `blom_strict_local_qp.py`: window construction, local problem assembly, solver wrapper, local diagnostics
- `blom_strict_local_kkt.py`: local Hessian, affine constraints, direct KKT solve, reduced-space solve
- `blom_strict_feasible_init.py`: explicit feasible local spline construction
- `phase3_plotting.py`: required Phase 3 figure generation
- `validate_phase3_blom_strict.py`: batch validation runner and artifact export
- `test_phase3_blom_strict.py`: regression and smoke tests
- `results/`: figures, tables, and logs

## Inputs and Outputs

Inputs follow the Phase 1 notation:

- `q = (q_0, ..., q_M)`
- `T = (T_1, ..., T_M)`
- `zeta_start`, `zeta_end` for the physical boundary jets
- canonical setting `s = 4`, `k = 2`

The main entry points are:

```python
build_window(i: int, k: int, M: int) -> dict
build_local_problem(q, T, i, k, zeta_start=None, zeta_end=None) -> dict
solve_blom_strict_local_qp(problem: dict, method: str = "kkt") -> dict
extract_segment_coeff(solution: dict, seg_idx: int) -> np.ndarray
```

## Running the Validation Script

```bash
python3 -m phase_3.validate_phase3_blom_strict
```

## Running the Tests

```bash
python3 -m unittest phase_3.test_phase3_blom_strict
```

## Result Directories

```text
phase_3/results/
├── figures/
├── tables/
└── logs/
```

Generated figures include:

- window layout
- local trajectory and derivatives
- continuity jumps
- natural boundary residuals
- perturbation response
- multistart uniqueness

Generated CSV tables include:

- `table_feasibility_summary.csv`
- `table_uniqueness_summary.csv`
- `table_perturbation_continuity.csv`

## How to Read the Main Diagnostics

- `uniqueness`: small pairwise coefficient differences and tiny objective spread across multistart runs mean the local optimum is not an initialization artifact.
- `natural BC`: small artificial-boundary residuals for orders `4,5,6` support the variational natural boundary conditions.
- `continuity`: tiny jumps up to order `6` support the predicted automatic `C^6` regularity for `s = 4`.
- `perturbation`: if `||delta c||` shrinks with `eps`, the local map is continuous under small waypoint perturbations.

## Scope Boundary

Phase 3 intentionally does **not**:

- derive the explicit analytic local system of Phase 4,
- prove global assembly consistency across neighboring BLOM windows,
- recover the exact global MINCO coefficient map from local windows.

That last point matters: Phase 2 showed that exact finite local support is not
available for the global MINCO map. Phase 3 instead shows that a **local
variational replacement** is still a legitimate, unique, and numerically stable
mathematical object.

