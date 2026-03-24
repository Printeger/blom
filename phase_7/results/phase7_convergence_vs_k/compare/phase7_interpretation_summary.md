# Phase 7 Interpretation Summary

## Exponential-Decay Reading
- Ideal truncation: decreases with k, but the log-linear fit is only moderate (slope=-4.686e+00, R^2=0.686).
- Scheme A: essentially flat across the tested k range (slope=-3.553e-15, R^2=1.000).
- Scheme B: numerically consistent with approximate exponential decay (slope=-2.070e-01, R^2=0.922).
- Scheme C: numerically consistent with approximate exponential decay (slope=-3.160e-01, R^2=0.996).

## Best MINCO Approximation in This Run
- Best final actual scheme by coefficient error: `Scheme A`
- Best final strict k-family by coefficient error: `Scheme C`
- Final ideal error: `2.169e-13`
- Scheme A: final error `7.400e-12`, matching `7.455e-12`, relative cost gap `5.351e-14`
- Scheme B: final error `1.417e+02`, matching `1.417e+02`, relative cost gap `2.088e+00`
- Scheme C: final error `3.469e+00`, matching `3.469e+00`, relative cost gap `-1.419e-01`

## Matching Assessment
- Scheme A vs ideal: decreases with k, but the log-linear fit is only moderate (slope=-4.156e+00, R^2=0.697).
- Scheme B vs ideal: numerically consistent with approximate exponential decay (slope=-2.071e-01, R^2=0.922).
- Scheme C vs ideal: numerically consistent with approximate exponential decay (slope=-2.939e-01, R^2=0.990).

## Theory Recommendation
- If the ideal truncation decays cleanly while actual Scheme B/C track it closely, the next theory target should remain the missing matching theorem.
- If actual errors stagnate while ideal truncation still decays, the priority should shift from theorem completion toward revising the assembly or the comparison target.
- Cost-gap trends should be read as diagnostic evidence only; they do not by themselves prove near-optimality because the actual assembled candidates are not the full Phase 1 affine-feasible family.

- Scheme A is currently treated as a shared-state baseline. In the present Phase 5 implementation its assembly ignores k, so its Phase 7 trend should be read as a baseline, not as a strict BLOM-k family.
