# Phase 7 Convergence Validation

Phase 7 validates how BLOM-k behaves as the window parameter grows, using the
exact Phase 1 MINCO solution as the reference object.

The implementation separates three objects:

- `MINCO` exact global coefficients from `phase_1/minco_scalar_baseline.py`
- `ideal truncated BLOM-k`, built by truncating the exact kernel `G(T)=A(T)^{-1}S_q`
- `actual BLOM-k`, built from the Phase 5 assembly schemes

## Main module

`phase_7/blom_convergence_vs_k.py`

Key APIs:

- `make_k_grid(M, start=2, step=2, include_M=True)`
- `compute_minco_reference(q, T, s=4)`
- `compute_ideal_truncated_blom_k(q, T, k, s=4)`
- `compute_actual_blom_k(q, T, k, s=4, scheme="C", config=None)`
- `compute_convergence_errors(...)`
- `fit_log_error_vs_k(...)`
- `run_convergence_vs_k(...)`
- `run_random_trials(...)`

## Artifacts

Saved under `phase_7/results/phase7_convergence_vs_k/`:

- `ideal/`
- `scheme_A/`
- `scheme_B/`
- `scheme_C/`
- `compare/`

The `compare/` directory contains the required cross-scheme figures and summary
tables, including:

- `coef_error_vs_k_all.png`
- `log_coef_error_fit_all.png`
- `matching_error_vs_k.png`
- `relative_cost_gap_vs_k.png`
- `convergence_errors_by_k.csv`
- `logfit_summary.csv`
- `summary_phase7_convergence.json`
- `phase7_interpretation_summary.md`

## Important note on Scheme A

The current Phase 5 implementation of Scheme A ignores the requested `k` and
acts as a shared-state global baseline. Phase 7 therefore reports Scheme A as a
useful baseline, but not as a strict BLOM-k family in the same sense as
Schemes B and C.

## Run commands

```bash
python3 -m phase_7.blom_convergence_vs_k
python3 -m phase_7.examples.demo_single_case
python3 -m phase_7.examples.demo_compare_schemes
python3 -m phase_7.examples.demo_ideal_vs_actual
python3 -m phase_7.examples.demo_random_trials
python3 -m unittest phase_7.test_blom_convergence_vs_k
```
