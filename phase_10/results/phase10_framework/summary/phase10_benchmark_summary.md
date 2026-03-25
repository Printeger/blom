# Phase 10 Benchmark Summary

- Raw Scheme C objective at k=2 available: `True`.
- Any general-k rows not yet implemented: `False`.
- Fastest baseline on representative case: `heuristic`.
- Best objective on representative case: `raw_schemeC`.

Interpretation:
- This benchmark layer is designed to compare quality, runtime, and sparsity under one common interface.
- When k > 2 is not yet fully implemented, the framework still records this explicitly instead of silently fabricating comparison rows.
