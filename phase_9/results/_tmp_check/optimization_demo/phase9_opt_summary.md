# Phase 9 Optimization Summary

- Objective drop: `2.688224e+03`.
- Control drop: `2.688807e+03`.
- Time drop: `-1.802551e-01`.
- Obstacle drop: `-4.023269e-01`.
- Gradient-norm drop proxy: `1.243217e+04`.
- Any duration touched the lower bound: `False`.

Interpretation:
- A negative objective drop would indicate the loop is unstable or the step size is too large.
- When the obstacle term decreases while durations stay feasible, the minimal loop is already useful as an optimization primitive.
- If durations repeatedly hit the lower bound, the next phase should add better time-parameterization handling rather than a larger objective.
