# RESULT_PHASE7

## Scope

This document records the actual Phase 7 execution results for the convergence-validation toolbox in `phase_7/`.
The goal of Phase 7 is to compare:

- the exact global MINCO reference,
- the idealized truncated BLOM-k surrogate,
- the actual assembled BLOM-k trajectories from Schemes A / B / C.

All commands below were executed in this repository workspace.

## Step 1: Compile The Phase 7 Package

### Command

```bash
python3 -m compileall phase_7
```

### Function

This checks that the new Phase 7 package, examples, and tests are syntactically valid and importable.

### Result

- `phase_7/__init__.py`
- `phase_7/blom_convergence_vs_k.py`
- `phase_7/examples/*.py`
- `phase_7/test_blom_convergence_vs_k.py`

all compiled successfully.

### Analysis

This is the first sanity gate. Passing `compileall` means the package structure and module imports are coherent enough to proceed to numerical validation.

## Step 2: Run The Main Phase 7 Convergence Sweep

### Command

```bash
python3 -m phase_7.blom_convergence_vs_k
```

### Function

This is the main Phase 7 entry point. It performs:

1. exact MINCO reference construction,
2. ideal truncated BLOM-k computation over the `k` grid,
3. actual Scheme A / B / C computation over the same `k` grid,
4. coefficient-error, matching-error, and cost-gap evaluation,
5. random-trial summary generation,
6. figure, CSV, JSON, and interpretation-summary export.

### Result

Console output:

```text
Phase 7 convergence validation
k grid: [2, 4, 6, 8]
reference cost: 3.994198e+04
ideal final error: 2.168534e-13
scheme A: final error=7.399695e-12, matching=7.454979e-12, rel cost gap=5.351275e-14
scheme B: final error=1.416847e+02, matching=1.416847e+02, rel cost gap=2.088343e+00
scheme C: final error=3.469311e+00, matching=3.469311e+00, rel cost gap=-1.419301e-01
random success rates: ideal=1.00, A=1.00, B=1.00, C=1.00
results dir: phase_7/results/phase7_convergence_vs_k
```

Saved artifacts:

- `phase_7/results/phase7_convergence_vs_k/compare/coef_error_vs_k_all.png`
- `phase_7/results/phase7_convergence_vs_k/compare/log_coef_error_fit_all.png`
- `phase_7/results/phase7_convergence_vs_k/compare/matching_error_vs_k.png`
- `phase_7/results/phase7_convergence_vs_k/compare/relative_cost_gap_vs_k.png`
- `phase_7/results/phase7_convergence_vs_k/compare/convergence_errors_by_k.csv`
- `phase_7/results/phase7_convergence_vs_k/compare/logfit_summary.csv`
- `phase_7/results/phase7_convergence_vs_k/compare/summary_phase7_convergence.json`
- `phase_7/results/phase7_convergence_vs_k/compare/phase7_interpretation_summary.md`
- per-object summaries in `ideal/`, `scheme_A/`, `scheme_B/`, `scheme_C/`

Core `k`-sweep results from `convergence_errors_by_k.csv`:

| Curve | k=2 | k=4 | k=6 | k=8 |
| --- | ---: | ---: | ---: | ---: |
| Ideal error `||c^(k)-c*||_2` | `4.451e+00` | `1.052e+00` | `1.795e-01` | `2.169e-13` |
| Scheme A error | `7.400e-12` | `7.400e-12` | `7.400e-12` | `7.400e-12` |
| Scheme B error | `4.836e+02` | `4.300e+02` | `2.722e+02` | `1.417e+02` |
| Scheme C error | `2.316e+01` | `1.345e+01` | `7.193e+00` | `3.469e+00` |

Log-fit results from `logfit_summary.csv`:

| Curve | Slope | R^2 |
| --- | ---: | ---: |
| Ideal | `-4.686e+00` | `0.686` |
| Scheme A | `-3.553e-15` | `1.000` |
| Scheme B | `-2.070e-01` | `0.922` |
| Scheme C | `-3.160e-01` | `0.996` |

### Analysis

- The ideal truncated object behaves exactly as Phase 7 theory hopes: its coefficient error drops rapidly with `k`, and by `k = 8` it is numerically identical to the exact MINCO reference.
- Scheme C shows a clean empirical decay trend and is the best strict `k`-dependent family in this implementation.
- Scheme B also improves with `k`, but the gap to MINCO stays very large, which is evidence that its current assembly does not numerically match the ideal truncation well.
- Scheme A is nearly exact relative to MINCO here, but it is effectively flat across `k` because the current Phase 5 implementation ignores `k`. In Phase 7 it should therefore be interpreted as a shared-state baseline, not as a genuine BLOM-k family.
- Scheme C has a negative relative cost gap at `k = 8`. This does not mean it "beats" MINCO. It means the candidate lies outside the full Phase 1 affine-feasible family, so its snap cost can be lower while still violating the exact global continuity structure that defines the MINCO mother problem.

## Step 3: Run The Single-Case Demo

### Command

```bash
python3 -m phase_7.examples.demo_single_case
```

### Function

This is the smallest user-facing Phase 7 demo. It runs the deterministic representative case and reports the final coefficient error for the ideal truncation and each assembly scheme.

### Result

```text
Phase 7 single-case convergence demo
k grid: [2, 4, 6, 8]
ideal final error: 2.168534e-13
scheme A final error: 7.399695e-12
scheme B final error: 1.416847e+02
scheme C final error: 3.469311e+00
```

### Analysis

This step confirms that the single-case summary matches the saved CSV tables. It is the quickest way to see the ranking:

- ideal truncation is effectively exact at the terminal `k`,
- Scheme A baseline is extremely close,
- Scheme C is the best actual `k`-dependent implementation,
- Scheme B is still far from MINCO.

## Step 4: Run The Scheme-Comparison Demo

### Command

```bash
python3 -m phase_7.examples.demo_compare_schemes
```

### Function

This demo emphasizes the actual Scheme A / B / C comparison at the largest tested `k`, including matching error and relative cost gap.

### Result

```text
Phase 7 scheme comparison demo
scheme A: final error=7.399695e-12, matching=7.454979e-12, rel cost gap=5.351275e-14
scheme B: final error=1.416847e+02, matching=1.416847e+02, rel cost gap=2.088343e+00
scheme C: final error=3.469311e+00, matching=3.469311e+00, rel cost gap=-1.419301e-01
```

### Analysis

- Scheme A aligns almost exactly with MINCO on this case, but again this is a baseline effect, not evidence of `k`-driven convergence.
- Scheme B has both large approximation error and large matching error, so the current consensus assembly is not numerically close to the ideal truncated model in this experiment.
- Scheme C has much smaller error than B and a much smaller matching gap, which makes it the most promising strict local family among the actual assemblies.

## Step 5: Run The Ideal-Vs-Actual Demo

### Command

```bash
python3 -m phase_7.examples.demo_ideal_vs_actual
```

### Function

This demo isolates the comparison between the ideal truncated object and the actual `k`-dependent assemblies B and C.

### Result

```text
Phase 7 ideal-vs-actual demo
ideal final error: 2.168534e-13
scheme B: final error=1.416847e+02, matching=1.416847e+02
scheme C: final error=3.469311e+00, matching=3.469311e+00
```

### Analysis

At the terminal `k`, the ideal truncation is already exact, so the matching error and the MINCO error coincide for Schemes B and C. This is useful because it cleanly exposes the missing bridge:

- the theory for the ideal truncation is working,
- the actual assembly still needs a stronger matching theorem or a revised construction,
- Scheme C is much closer to that target than Scheme B.

## Step 6: Run The Random-Trials Demo

### Command

```bash
python3 -m phase_7.examples.demo_random_trials
```

### Function

This demo runs a smaller randomized batch and reports success rate, mean fitted slope, and mean `R^2` for each curve.

### Result

```text
Phase 7 random-trials demo
ideal: success=1.00, mean slope=-3.191321e+00, mean R^2=0.661783
A: success=1.00, mean slope=-1.588822e-15, mean R^2=0.325000
B: success=1.00, mean slope=-2.746292e-01, mean R^2=0.949765
C: success=1.00, mean slope=-3.577788e-01, mean R^2=0.990983
```

Additional aggregate values from the saved `random_trials_aggregate.json` generated by the main run:

- success rates: `ideal=1.00`, `A=1.00`, `B=1.00`, `C=1.00`
- median final matching error: `A=3.114e-12`, `B=1.986e+01`, `C=1.171e-01`
- median final relative cost gap: `A=-4.724e-14`, `B=1.516e-01`, `C=1.136e-03`

### Analysis

- Random trials reinforce the single-case picture: Scheme C is the most stable strict `k`-dependent family.
- Scheme B usually decays, but not nearly as strongly as C.
- Scheme A remains flat because it is not really changing with `k`.
- The ideal truncation continues to show strong negative slopes, but the `R^2` is more variable because the error often collapses toward machine precision near the largest `k`.

## Step 7: Run Phase 7 Unit Tests

### Command

```bash
python3 -m unittest phase_7.test_blom_convergence_vs_k
```

### Function

This validates the implementation-level contract:

- `k` grid generation,
- MINCO reference construction,
- full-window ideal truncation recovering the exact solution,
- actual Scheme A / B / C execution,
- artifact generation,
- random-trial summary export.

### Result

```text
.....
----------------------------------------------------------------------
Ran 5 tests in 1.076s

OK
```

### Analysis

The targeted Phase 7 tests passed without failure, which means the new interfaces and required artifact outputs are stable enough for repository integration.

## Step 8: Run Full Repository Regression Tests

### Command

```bash
python3 -m unittest discover -s . -p 'test*.py'
```

### Function

This checks that Phase 7 did not break earlier phases.

### Result

```text
.....................
----------------------------------------------------------------------
Ran 50 tests in 12.385s

OK
```

### Analysis

Phase 7 integrates cleanly with the earlier phases. This is important because Phase 7 reuses:

- Phase 1 global MINCO construction,
- Phase 2 kernel extraction logic,
- Phase 5 assembly schemes.

Passing the full test suite means the new convergence layer did not regress the older mathematical infrastructure.

## Final Takeaway

The Phase 7 data supports three clear conclusions:

1. The ideal truncated BLOM-k surrogate converges strongly toward the exact MINCO reference as `k` grows.
2. Among the actual `k`-dependent assemblies, Scheme C is the most convincing approximation family in the current codebase.
3. The main missing bridge is still the matching theorem, especially for explaining why the ideal truncation behaves so well while the actual assemblies, particularly Scheme B, remain much farther away.
