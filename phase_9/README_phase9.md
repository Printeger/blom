# Phase 9 Minimal Differentiable Loop

Phase 9 verifies that raw Scheme C already works as a differentiable
trajectory representation in the canonical case `s=4, k=2`.

Implemented objects:

- raw Scheme C coefficients and dense Jacobians
- control / time / smooth obstacle objective terms
- dense and block-banded reverse differentiation
- finite-difference gradient checks
- a minimal projected gradient-descent demo
- a unified validation suite with figures, CSV tables, JSON logs, and summaries

Main entry points:

- `phase_9/blom_backward_diff.py`
- `phase_9/blom_space_time_opt_demo.py`
- `phase_9/blom_phase9_validation_suite.py`

Recommended commands:

```bash
python3 -m phase_9.blom_backward_diff
python3 -m phase_9.blom_space_time_opt_demo
python3 -m phase_9.blom_phase9_validation_suite
python3 -m unittest phase_9.test_blom_backward_diff
```

Scope note:

- This is a minimal differentiable loop, not a full planner.
- The main deliverable is a trustworthy backward pass and a stable descent demo.
