# Phase 8 Interior-First Matching Validation

Phase 8 adds an interior-first matching toolbox on top of the Phase 7
convergence experiments.

Implemented objects:

- raw Scheme C coefficients `c^{C,k}`
- ideal truncated coefficients `\tilde c^{(k)}`
- interior and boundary block partitions with configurable radius `r(k)`
- augmented local window solves with artificial high-order boundary traces
- boundary-gap decomposition and regime-split experiments

Main entry points:

- `phase_8/blom_interior_matching_check.py`
- `phase_8/blom_boundary_gap_decomposition.py`
- `phase_8/blom_uniform_vs_nonuniform_interior.py`
- `phase_8/blom_phase8_validation_suite.py`

Default result directory:

- `phase_8/results/phase8_validation/`

Recommended commands:

```bash
python3 -m phase_8.blom_interior_matching_check
python3 -m phase_8.blom_boundary_gap_decomposition
python3 -m phase_8.blom_uniform_vs_nonuniform_interior
python3 -m phase_8.blom_phase8_validation_suite
```

Scope note:

- Phase 8 does not claim a new theorem.
- It builds the numerical bridge needed to decide whether the next proof target
  should be an interior-first theorem, a boundary-layer remainder theorem, or a
  reference-window consistency proposition.

## Supplementary Experiments

Phase 8 also now includes a supplementary experiment layer for stronger
robustness checks on top of the main validation suite.

Supplementary entry points:

- `phase_8/exp1_large_M_sweep.py`
- `phase_8/exp2_more_trials.py`
- `phase_8/exp3_radius_sensitivity.py`
- `phase_8/exp4_raw_vs_reference_sanity.py`
- `phase_8/exp5_empty_interior_risk.py`
- `phase_8/exp6_two_bridge_gap_compare.py`
- `phase_8/blom_phase8_supplementary_suite.py`

Supplementary result directory:

- `phase_8/results/phase8_supplementary/`

Recommended supplementary commands:

```bash
python3 -m phase_8.exp1_large_M_sweep
python3 -m phase_8.exp2_more_trials
python3 -m phase_8.exp3_radius_sensitivity
python3 -m phase_8.exp4_raw_vs_reference_sanity
python3 -m phase_8.exp5_empty_interior_risk
python3 -m phase_8.exp6_two_bridge_gap_compare
python3 -m phase_8.blom_phase8_supplementary_suite
```

Supplementary interpretation target:

- test whether the interior-first effect survives larger `M` and many trials,
- rule out empty-interior false positives,
- stress-test the very small `raw vs reference-window` gap,
- identify whether the dominant unresolved bridge is still
  `reference-window -> ideal truncation`.
