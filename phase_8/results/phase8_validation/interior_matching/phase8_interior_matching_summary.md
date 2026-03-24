# Phase 8 Interior Matching Summary

- Interior-only matching beats full matching for `5/5` tested k values.
- Boundary-only matching exceeds interior-only matching for `5/5` tested k values.
- Interior log-fit slope: `-5.265687` with `R^2 = 0.761420`.
- Full log-fit slope: `-0.446569`; boundary log-fit slope: `-0.418230`.
- Largest tested k = `10` gives interior/full ratio `0.000000`.

Interpretation:
- Interior-only matching is numerically cleaner whenever its curve sits well below the full-domain curve.
- A more negative interior slope than the full-domain slope supports an interior-first matching theorem.
- Persistent boundary dominance suggests the remaining gap should be treated as a dedicated boundary-layer remainder, not hidden inside a single full-domain statement.
