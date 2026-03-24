# Phase 5 Interpretation Summary

## Scheme A
- lower-order jumps reach machine tolerance: `True`
- requires extra global shared knot-state variables
- max lower-order jump: `6.580e-12`
- max higher-order jump: `1.123e-09`

## Scheme B
- lower-order jumps reach machine tolerance: `True`
- consensus max pre-dispersion: `2.008e+01`
- consensus max post-jump norm: `1.034e-11`
- preserves independent sliding-window solves and then projects to a shared knot state

## Scheme C
- exact C0 only: `order 0 is_zero_tol = True`
- any lower-order derivative beyond position guaranteed zero: `False`
- max ||eta^- - eta^+||_2: `4.016e+01`
- keeps the cleanest local-map interpretation but does not automatically guarantee C^{s-1}

## Recommendation
- Scheme B is the safest assembly choice when exact global C^{s-1} continuity and local-window independence are both desired.
- Scheme A is also smooth up to order s-1, but it introduces a true global shared-state optimization layer.
- Scheme C remains useful as a diagnostic local-map baseline, not as the default globally smooth assembly mechanism.
