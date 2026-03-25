from __future__ import annotations

import sys
from pathlib import Path

import streamlit as st

DEMO_ROOT = Path(__file__).resolve().parent
REPO_ROOT = DEMO_ROOT.parent
for candidate in (str(REPO_ROOT), str(DEMO_ROOT)):
    if candidate not in sys.path:
        sys.path.insert(0, candidate)

from demo_core.data_interface import repo_root


st.set_page_config(
    page_title="BLOM Research Demo",
    page_icon="📈",
    layout="wide",
    initial_sidebar_state="expanded",
)


def main() -> None:
    st.title("BLOM Research Demo")
    st.caption("Scientifically aligned with the existing BLOM phases, focused on real numerical objects rather than schematic curves.")

    st.markdown(
        """
This app is a research-oriented visualization layer over the existing BLOM codebase.
It is designed to compare:

- `raw Scheme C` as the primary BLOM object
- `MINCO` as the global interpolating reference
- `B-spline` as a local geometric-control baseline
- optional `Scheme A / B / C` continuity-vs-locality trade-offs
"""
    )

    col1, col2, col3 = st.columns(3)
    col1.metric("Backend Root", str(repo_root().name))
    col2.metric("Primary Locality Evidence", "Jacobian sparsity")
    col3.metric("Primary Optimization Evidence", "Phase 9/10 results")

    st.info(
        "Use the page list in the sidebar to explore perturbation response, Jacobian sparsity, k-sweeps, matching bridges, optimizer behavior, and benchmark summaries."
    )

    st.subheader("Scientific Guardrails")
    st.markdown(
        """
- BLOM locality claims are supported mainly by Jacobian sparsity, support width, outside-band sensitivity, and scaling summaries.
- B-spline control points are not treated as BLOM/MINCO waypoints.
- If a page falls back to cached or illustrative data, it is labeled explicitly.
- Phase 8 and Phase 10 pages prefer loading real result artifacts already produced in this repository.
"""
    )

    st.subheader("Run Command")
    st.code("streamlit run demo/app.py", language="bash")

    st.subheader("Project Paths")
    st.markdown(
        f"""
- Demo root: [{Path(__file__).resolve().parent}]({Path(__file__).resolve().parent})
- Results cache: [{Path(__file__).resolve().parent / "results_cache"}]({Path(__file__).resolve().parent / "results_cache"})
- Existing BLOM traceability: [{repo_root() / "Traceability.md"}]({repo_root() / "Traceability.md"})
"""
    )


if __name__ == "__main__":
    main()
