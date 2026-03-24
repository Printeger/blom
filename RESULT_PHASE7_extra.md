# RESULT_PHASE7_extra

## Scope

This document records the actual execution results for the extra Phase 7 experiments defined in `phase_7/REQ-Phase7-extra-experiments-v2.md`.

The extra experiment framework adds four follow-up studies on top of the main
Phase 7 convergence validation:

1. large-`M` sweep,
2. uniform-time vs bounded-nonuniform split,
3. interior-only vs full-error separation,
4. light corrected Scheme C variants.

All commands below were executed in this repository workspace.

## Step 1: Compile The Updated Phase 7 Package

### Command

```bash
python3 -m compileall phase_7
```

### Function

This checks the new extra-experiment module, its demos, and its tests for syntax and import correctness.

### Result

`phase_7/blom_phase7_extra_experiments.py`, all new demo scripts, and the new
test module compiled successfully.

### Analysis

This confirms that the extra experiments are integrated into the existing Phase
7 package cleanly enough to proceed to numerical execution.

## Step 2: Run The Main Extra-Experiment Entry Point

### Command

```bash
python3 -m phase_7.blom_phase7_extra_experiments
```

### Function

This is the top-level driver for all four new experiment families. It runs:

1. the large-`M` sweep,
2. the time-regime split,
3. the interior-vs-full study,
4. the Scheme C correction study,
5. the final cross-experiment overview table and interpretation summary.

### Result

Console output:

```text
Phase 7 extra experiments
overview rows: 14
overview csv: phase_7/results/phase7_extra_experiments/compare_summary/final_experiment_overview.csv
```

Saved directories:

- `phase_7/results/phase7_extra_experiments/exp1_large_M_sweep/`
- `phase_7/results/phase7_extra_experiments/exp2_time_regime_split/`
- `phase_7/results/phase7_extra_experiments/exp3_interior_vs_full/`
- `phase_7/results/phase7_extra_experiments/exp4_schemeC_light_assembly/`
- `phase_7/results/phase7_extra_experiments/compare_summary/`

Key rows from `compare_summary/final_experiment_overview.csv`:

| Experiment | Curve | Best slope | Best R^2 | Final coef error | Final match error | Final cost rel |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| exp1_large_M_sweep | ideal | `-1.212e+00` | `0.690` | `4.582e-13` | `nan` | `1.976e-14` |
| exp1_large_M_sweep | C | `-3.244e-01` | `0.998` | `6.263e-02` | `6.263e-02` | `-1.380e-03` |
| exp2_time_regime_split | C | `-3.253e-01` | `0.997` | `1.713e-02` | `1.713e-02` | `-1.156e-04` |
| exp3_interior_vs_full | C | `-3.271e-01` | `0.978` | `1.564e+00` | `1.564e+00` | `-1.419e-01` |
| exp4_schemeC_light_assembly | C | `-3.160e-01` | `0.996` | `3.469e+00` | `3.469e+00` | `-1.419e-01` |
| exp4_schemeC_light_assembly | C1 | `-3.411e-01` | `0.969` | `5.984e+01` | `5.984e+01` | `3.231e-01` |
| exp4_schemeC_light_assembly | C2 | `-3.457e-01` | `0.979` | `4.990e+01` | `4.990e+01` | `2.544e-01` |

### Analysis

- The ideal truncation remains the strongest object numerically, exactly as the
  Phase 7 theory suggests.
- Across the actual `k`-dependent families, raw Scheme C remains the only
  consistently strong candidate in the extra experiments.
- The large-`M` and time-regime overview rows both still support the
  "continue toward the matching theorem" direction for Scheme C.
- The corrected Scheme C variants `C1/C2` are the most important negative
  result in this batch: they improve continuity, but badly damage approximation
  quality and cost behavior.

## Step 3: Run Experiment 1 Demo, Large-M Sweep

### Command

```bash
python3 -m phase_7.examples.demo_exp1_large_M
```

### Function

This demo runs the large-`M` sweep on a representative scale set and saves the
per-`M` figures and summary table.

### Result

```text
Phase 7 extra exp1: large-M sweep
cases: [10, 20, 40]
```

Selected values from `exp1_large_M_sweep/summary_large_M.csv`:

| M | Curve | Slope | Final coef error |
| --- | --- | ---: | ---: |
| 10 | ideal | `-3.160e+00` | `1.274e-13` |
| 10 | C | `-4.158e-01` | `2.532e-01` |
| 20 | ideal | `-1.212e+00` | `4.582e-13` |
| 20 | C | `-3.244e-01` | `1.325e-02` |
| 40 | ideal | `-6.730e-01` | `7.088e-05` |
| 40 | C | `-3.015e-01` | `6.263e-02` |

### Analysis

- Scheme C keeps a stable negative slope across `M = 10, 20, 40`, which is a
  strong positive signal for the current theory direction.
- The ideal truncation remains excellent, but at `M = 40` its final error is
  no longer machine-zero because this experiment capped the sweep at `k_max = 20`.
  That is expected and actually useful: it shows the truncated-kernel picture
  still improving rather than saturating artificially.
- Scheme B remains far away from MINCO as `M` grows, so increasing problem
  size does not rescue the current consensus assembly.

## Step 4: Run Experiment 2 Demo, Time-Regime Split

### Command

```bash
python3 -m phase_7.examples.demo_exp2_time_regime
```

### Function

This demo compares uniform-time and bounded-nonuniform time allocations, with
multiple small random trials per regime.

### Result

```text
Phase 7 extra exp2: time-regime split
records: 72
```

Representative rows from `exp2_time_regime_split/summary_time_regime.csv`:

| Regime | Curve | Slope | Final coef error |
| --- | --- | ---: | ---: |
| `uniform_h=0.5` | C | `-3.387e-01` | `4.124e-01` |
| `uniform_h=1.0` | C | `-3.379e-01` | `7.675e-03` |
| `uniform_h=2.0` | C | `-3.427e-01` | `1.802e-03` |
| `bounded_0.9_1.1` | C | `-3.308e-01` | `2.507e-02` |
| `bounded_0.5_2.0` | C | `-3.151e-01` | `1.063e-02` |
| `bounded_0.2_3.0` | C | `-3.356e-01` | `1.343e-02` |

### Analysis

- Scheme C stays remarkably stable across both uniform and bounded-nonuniform
  regimes. Its slopes remain close to `-0.32` to `-0.35` in all six groups.
- Uniform-time is still the cleanest regime numerically, especially at larger
  step size (`h = 2.0`), but bounded-nonuniform behavior does not collapse.
- This means the next theorem can still reasonably aim beyond the strictly
  uniform regime, even if uniform-time remains the easiest place to prove it first.

## Step 5: Run Experiment 3 Demo, Interior Vs Full Error

### Command

```bash
python3 -m phase_7.examples.demo_exp3_interior_vs_full
```

### Function

This demo separates boundary-driven and interior-driven error by trimming two
 segments at each end when computing the interior-only metrics.

### Result

```text
Phase 7 extra exp3: interior vs full
k grid: [2, 4, 6, 8]
```

Selected rows from `exp3_interior_vs_full/summary_interior_vs_full.csv`:

| Curve | k | Full coef error | Interior coef error |
| --- | ---: | ---: | ---: |
| ideal | 2 | `4.451e+00` | `3.613e+00` |
| ideal | 6 | `1.795e-01` | `1.749e-13` |
| ideal | 8 | `2.169e-13` | `1.749e-13` |
| C | 2 | `2.316e+01` | `1.110e+01` |
| C | 6 | `7.193e+00` | `3.723e+00` |
| C | 8 | `3.469e+00` | `1.564e+00` |

### Analysis

- Boundary effects are clearly important. The strongest evidence is the ideal
  truncation at `k = 6`: the full error is still `1.795e-01`, but the
  interior-only error has already collapsed to `1.749e-13`.
- For raw Scheme C, trimming the boundary also cuts the error substantially,
  but not completely. At `k = 8`, the full error is `3.469`, while the
  interior-only error is still `1.564`.
- This means the next theorem path should likely include an interior-first
  version, but boundary effects are not the whole story for actual BLOM.

## Step 6: Run Experiment 4 Demo, Light Corrected Scheme C

### Command

```bash
python3 -m phase_7.examples.demo_exp4_schemeC_correction
```

### Function

This demo compares raw Scheme C with two lightweight corrected variants:

- `C1`: endpoint-jet averaging
- `C2`: weighted endpoint-jet reconciliation

### Result

```text
Phase 7 extra exp4: Scheme C correction
rows: 12
```

Key rows from `exp4_schemeC_light_assembly/summary_schemeC_correction.csv`:

| Curve | k | Full coef error | Full matching error | Rel cost gap | Lower-order jump max |
| --- | ---: | ---: | ---: | ---: | ---: |
| C | 8 | `3.469e+00` | `3.469e+00` | `-1.419e-01` | `6.592e+00` |
| C1 | 8 | `5.984e+01` | `5.984e+01` | `3.231e-01` | `6.015e-12` |
| C2 | 8 | `4.990e+01` | `4.990e+01` | `2.544e-01` | `4.731e-12` |

Locality-width check for the final `k`:

| Curve | q bandwidth |
| --- | ---: |
| C | `5` |
| C1 | `6` |
| C2 | `6` |

### Analysis

- This is the most important negative result in the whole extra-experiment batch.
- `C1` and `C2` do exactly what they were designed to do for continuity:
  lower-order jumps drop from `6.592e+00` in raw C to about `1e-12`.
- But that continuity gain comes at a steep price: both corrected variants are
  much worse than raw C in coefficient error, matching error, and relative cost gap.
- They also widen the `q`-Jacobian bandwidth from `5` to `6`, so the correction
  is not only less accurate, but also slightly less local.
- The current conclusion is therefore clear: this very light correction should
  not be promoted to the formal BLOM object in its current form.

## Step 7: Run The All-In-One Demo

### Command

```bash
python3 -m phase_7.examples.demo_all_phase7_extra
```

### Function

This demo runs the whole extra-experiment pipeline on a moderate default setup
and confirms that the overview outputs are generated end-to-end.

### Result

```text
Phase 7 extra experiments
overview rows: 14
```

### Analysis

This confirms the full integrated workflow is operational, not just the
individual experiment components.

## Step 8: Run Phase 7 Extra Unit Tests

### Command

```bash
python3 -m unittest phase_7.test_blom_phase7_extra_experiments
```

### Function

This validates:

- large-`M` sweep execution,
- time-regime split execution,
- interior/full metric generation,
- corrected Scheme C reconstruction,
- overview artifact generation.

### Result

```text
----------------------------------------------------------------------
Ran 5 tests in 21.288s

OK
```

### Analysis

The new experiment layer is stable enough to run end-to-end and generate its
required artifacts without manual intervention.

## Step 9: Re-Run The Original Phase 7 Unit Tests

### Command

```bash
python3 -m unittest phase_7.test_blom_convergence_vs_k
```

### Function

This checks that the new extra-experiment layer did not break the original
Phase 7 convergence-validation module.

### Result

```text
----------------------------------------------------------------------
Ran 5 tests in 1.052s

OK
```

### Analysis

The original convergence-validation interface still works after adding the
follow-up experiment framework.

## Step 10: Run Full Repository Regression Tests

### Command

```bash
python3 -m unittest discover -s . -p 'test*.py'
```

### Function

This confirms that the extra Phase 7 experiments did not regress earlier phases.

### Result

```text
----------------------------------------------------------------------
Ran 55 tests in 32.958s

OK
```

### Analysis

The repository remains consistent after adding the new experiment family.

## Final Takeaway

The extra Phase 7 experiments sharpen the research picture substantially:

1. large-`M` scaling still supports Scheme C as the most credible actual
   `k`-dependent family;
2. bounded-nonuniform time allocations do not destroy the main trend, so the
   theory path does not need to retreat immediately to a uniform-only claim;
3. boundary effects are real and strong, so an interior theorem would be a
   sensible intermediate milestone;
4. the current `C1/C2` light corrections are not a good replacement for raw
   Scheme C, because they improve continuity but significantly worsen matching,
   cost, and locality.
