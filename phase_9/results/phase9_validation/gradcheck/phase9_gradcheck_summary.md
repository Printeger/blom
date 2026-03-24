# Phase 9 Gradcheck Summary

- Total analytic dense gradient max abs error vs FD: `3.755835e-05`.
- Total analytic banded gradient max abs error vs FD: `3.755835e-05`.
- Dense vs banded max abs gap: `1.818989e-12`.
- Control-term partial gradient max abs error: `1.485118e-05`.
- Time-term partial gradient max abs error: `2.631780e-11`.
- Obstacle-term partial gradient max abs error: `4.200990e-09`.
- Largest partial-gradient error term: `control`.

Interpretation:
- The finite-difference comparison checks whether the full analytic chain is correct.
- The dense-vs-banded comparison checks whether the local-support accumulation removes only exact zeros.
- If both gaps are small, the minimal loop is ready for an optimization demo.
