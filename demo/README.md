# BLOM Research Demo

This directory contains a Streamlit + Plotly research demo aligned with the
existing BLOM phases in this repository.

## Run

From the repository root:

```bash
streamlit run demo/app.py
```

## Pages

1. `Overview`
   - Structural comparison between B-spline, MINCO, raw Scheme C, and Scheme A/B.
   - This page is a comparison page, not a theorem page.

2. `Single Perturbation Response`
   - Geometry-level perturbation view.
   - B-spline uses control-point perturbation.
   - MINCO / BLOM use interpolation-waypoint perturbation.

3. `Jacobian Sparsity / Support Pattern`
   - Main locality evidence page.
   - Shows `d c / d q_bar` and `d c / d T` heatmaps, thresholded masks, support width, and outside-band energy.

4. `k-sweep`
   - Locality / approximation / runtime trade-off across `k`.
   - Compares raw Scheme C and ideal truncation.

5. `Scheme A / B / C Compare`
   - Continuity-vs-locality trade-off.
   - Shows jumps, support width, outside-band sensitivity, and matching-to-ideal.

6. `Phase 8 Matching Bridge`
   - Real Phase 8 triplet comparison:
     - actual raw Scheme C
     - reference-window
     - ideal truncation

7. `Optimizer Demo`
   - Phase 9 minimal objective and Phase 10 full objective.
   - Shows objective history, before/after trajectory, and dense-vs-sparse consistency.

8. `Benchmark / Ablation Summary`
   - Loads real Phase 10 benchmark outputs.
   - Shows runtime vs `M`, objective gap vs `k`, support-width/memory trends, baselines, and ablations.

## Backend Interface Assumptions

The demo reuses existing repository interfaces via adapters rather than
rewriting BLOM mathematics. In particular it expects:

- `phase_7.blom_convergence_vs_k.compute_minco_reference`
- `phase_7.blom_convergence_vs_k.compute_ideal_truncated_blom_k`
- `phase_7.blom_convergence_vs_k.compute_actual_blom_k`
- `phase_5.blom_boundary_jump_check.assemble_scheme_A/B/C`
- `phase_9.blom_backward_diff.compute_raw_schemeC_jacobians`
- `phase_10.blom_full_backward_diff.compute_raw_schemeC_jacobians_general_k`
- `phase_10.blom_full_backward_diff.evaluate_full_objective`
- `phase_10.blom_full_backward_diff.full_backward_diff_dense`
- `phase_10.blom_full_backward_diff.full_backward_diff_sparse`
- `phase_10.blom_space_time_optimizer.run_space_time_optimization`
- `phase_9.blom_space_time_opt_demo.run_minimal_optimization_demo`

## Strict vs Illustrative Modes

- `Overview`, `Jacobian`, `k-sweep`, `Scheme compare`, `Matching bridge`,
  `Optimizer`, and `Benchmark` pages all prefer strict numerical backends.
- `B-spline` views are geometric baselines and are explicitly labeled as such.
- If a required result artifact is missing, the app reports it explicitly instead
  of silently fabricating a result.

## Notes

- Static Plotly image export is not bundled here because it would require
  additional image-export dependencies such as `kaleido`.
- CSV download is supported on the benchmark page.
