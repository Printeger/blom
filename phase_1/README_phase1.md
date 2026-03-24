# Phase 1 MINCO Scalar Baseline

This directory implements the strengthened scalar MINCO mother problem for BLOM Phase 1.

## Mathematical Object

The implemented object is a scalar, degree-7, multi-segment trajectory with:

- waypoint positions `q_0, ..., q_M`
- positive segment durations `T_1, ..., T_M`
- fixed boundary jets up to order 3 at the start and end
- interior `C^6` continuity

The core linear system is

```text
M(T) c = b(q, zeta_start, zeta_end)
```

where:

- `c` is the stacked coefficient vector of shape `(8*M,)`
- `M(T)` depends only on segment durations
- `b(...)` contains waypoint and boundary-jet data

## Shapes

- `q`: `(M+1,)`
- `T`: `(M,)`
- `zeta_start`: `(4,)`
- `zeta_end`: `(4,)`
- `coeffs`: `(M, 8)`
- `c_vec`: `(8*M,)`
- `M`: `(8*M, 8*M)`
- `b`: `(8*M,)`

## Minimal Example

```python
import numpy as np
from phase_1.minco_scalar_baseline import solve_minco_coefficients

q = np.array([0.0, 1.0, 0.0])
T = np.array([1.0, 1.0])
zeta_start = np.array([0.0, 0.0, 0.0, 0.0])
zeta_end = np.array([0.0, 0.0, 0.0, 0.0])

result = solve_minco_coefficients(q, T, zeta_start, zeta_end, return_system=True)
print(result["coeffs"])
print(result["residual_norm"])
```

## Verification Checklist

The module exposes helpers for the key validation tasks:

- `interpolation_errors(coeffs, T, q)`
- `continuity_jumps(coeffs, T, max_order=6)`
- `boundary_jet_errors(coeffs, T, zeta_start, zeta_end)`
- `system_residual(M, c_vec, b)`

Run the Phase 1 tests with:

```bash
python3 -m unittest phase_1.test_minco_scalar_baseline
```

Run the example scripts with:

```bash
python3 -m phase_1.examples.demo_case_single_segment
python3 -m phase_1.examples.demo_case_two_segment
python3 -m phase_1.examples.demo_case_random
```

## Sensitivity-Ready Interfaces

The Phase 1 module already separates:

- system assembly
- right-hand-side assembly
- coefficient solve
- trajectory evaluation

It also reserves the following future interfaces for gradient propagation:

- `solve_adjoint(...)`
- `grad_wrt_q(...)`
- `grad_wrt_T(...)`

At this stage, `grad_wrt_q` and `grad_wrt_T` intentionally remain TODO placeholders.

## Known Limitations

- scalar output only
- fixed `s = 4`
- dense `numpy.linalg.solve` instead of a banded solver
- no BLOM local-window logic in this directory
