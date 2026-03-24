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
