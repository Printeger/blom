# Phase 9 Optimization Summary

- Objective drop: `2.686197e+03`.
- Control drop: `2.687937e+03`.
- Time drop: `-1.962653e-01`.
- Obstacle drop: `-1.543921e+00`.
- Gradient-norm drop proxy: `1.242874e+04`.
- Any duration touched the lower bound: `False`.

Interpretation:
- A negative objective drop would indicate the loop is unstable or the step size is too large.
- When the obstacle term decreases while durations stay feasible, the minimal loop is already useful as an optimization primitive.
- If durations repeatedly hit the lower bound, the next phase should add better time-parameterization handling rather than a larger objective.
