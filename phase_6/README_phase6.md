# Phase 6: FD Jacobian Check

Phase 6 turns the basic BLOM Jacobian propositions into a reproducible validation toolbox.

The main implementation is [blom_fd_jacobian_check.py](/home/dev/code/blom/phase_6/blom_fd_jacobian_check.py). It verifies:

- raw BLOM-Analytic local-support masks for `J_c_q` and `J_c_T`
- analytic-vs-finite-difference Jacobian agreement
- weak time sensitivity through `||dc/dT||` statistics
- locality changes after Scheme A / B / C assembly

## Main API

- `theoretical_mask_c_q(M, i, s=4, k=2)`
- `theoretical_mask_c_T(M, i, s=4, k=2)`
- `finite_difference_jacobian(f, x, eps=1e-6, method="central", mask=None)`
- `raw_local_coefficient_map(q, T, i, s=4, k=2)`
- `raw_local_jacobians(q, T, i, s=4, k=2, mode="analytic")`
- `assembled_scheme_A_map(q, T, i, config=None)`
- `assembled_scheme_B_map(q, T, i, config=None)`
- `assembled_scheme_C_map(q, T, i, config=None)`
- `run_fd_jacobian_check(...)`
- `run_compare_all_schemes(...)`
- `run_random_trials(...)`

## Reused Building Blocks

- Phase 2: Jacobian error metrics and bandwidth ideas
- Phase 4: canonical raw BLOM-Analytic local map
- Phase 5: assembled Scheme A / B / C coefficient maps

## How To Run

```bash
python3 -m phase_6.blom_fd_jacobian_check
python3 -m phase_6.examples.demo_raw_local_support
python3 -m phase_6.examples.demo_scheme_A_nonlocality
python3 -m phase_6.examples.demo_scheme_B_nonlocality
python3 -m phase_6.examples.demo_compare_all
python3 -m unittest phase_6.test_blom_fd_jacobian_check
```

## Artifacts

Artifacts are written under [phase_6/results/phase6_fd_jacobian_check](/home/dev/code/blom/phase_6/results/phase6_fd_jacobian_check):

- `raw_scheme_C/`: raw local-support masks and analytic-vs-FD comparisons
- `scheme_A/`: final assembled Scheme A Jacobian heatmaps and sparsity CSV
- `scheme_B/`: final assembled Scheme B Jacobian heatmaps and sparsity CSV
- `compare/`: bandwidth plots, step-size sweeps, bound statistics, random-trial summaries

## Interpretation

Phase 6 is expected to show a clean split:

- raw Scheme C / BLOM-Analytic stays exactly local in the canonical stencil
- Scheme A and Scheme B gain stronger global continuity but broaden the effective Jacobian footprint
- `dc/dT` remains well behaved on compact admissible samples, supporting the weak time-sensitivity proposition
