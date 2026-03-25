from __future__ import annotations

import sys
from pathlib import Path

import streamlit as st

DEMO_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = DEMO_ROOT.parent
for candidate in (str(REPO_ROOT), str(DEMO_ROOT)):
    if candidate not in sys.path:
        sys.path.insert(0, candidate)

from demo_core.blom_adapter import compute_matching_bridge
from demo_core.data_interface import case_for_demo
from demo_core.plotting import PAPER_COLORS, line_metric_figure


st.title("Phase 8 Matching Bridge")
st.caption("复用 Phase 8 的真实对象：actual raw Scheme C, reference-window, ideal truncation。")

col1, col2, col3 = st.columns(3)
M = col1.slider("M", 4, 10, 6, 1, help="Number of segments.")
regime = col2.selectbox("regime", ["bounded-nonuniform", "uniform"], index=0, help="Duration regime.")
k = col3.slider("k", 2, 8, 2, 2, help="BLOM window size.")
radius_mode = st.selectbox("radius mode", ["default", "minimal", "one_and_half"], index=0, help="How aggressively the interior set is trimmed away from the boundary layer.")

with st.expander("参数说明 / 名词解释 / 图表说明"):
    st.markdown(
        """
- `actual -> ideal`: raw Scheme C 到 ideal truncation 的差。
- `reference -> ideal`: reference-window 到 ideal truncation 的差。
- `actual -> reference`: raw Scheme C 到 reference-window 的差。
- `full l2`: 全部 segment blocks 的整体误差。
- `interior l2`: 只在 interior blocks 上统计的误差。
- `boundary l2`: 只在 boundary-layer blocks 上统计的误差。
- `boundary energy ratio`: 误差能量中有多少比例来自 boundary layer。
"""
    )

q, T = case_for_demo(M=M, regime=regime)
bridge = compute_matching_bridge(q, T, k=k, radius_mode=radius_mode)

records = [
    ("actual -> ideal", bridge["actual_vs_ideal"]),
    ("reference -> ideal", bridge["reference_vs_ideal"]),
    ("actual -> reference", bridge["actual_vs_reference"]),
]

metrics = {
    "full_l2": [rec["full_l2"] for _, rec in records],
    "interior_l2": [rec["interior_l2"] for _, rec in records],
    "boundary_l2": [rec["boundary_l2"] for _, rec in records],
    "boundary_energy_ratio": [rec["boundary_energy_ratio"] for _, rec in records],
}

labels = [name for name, _ in records]
st.plotly_chart(
    line_metric_figure(
        labels,
        [
            {"name": "full l2", "y": metrics["full_l2"], "color": PAPER_COLORS["minco"]},
            {"name": "interior l2", "y": metrics["interior_l2"], "color": PAPER_COLORS["blom"]},
            {"name": "boundary l2", "y": metrics["boundary_l2"], "color": PAPER_COLORS["schemeA"]},
        ],
        "Bridge error decomposition",
        "bridge",
        "error",
    ),
    use_container_width=True,
)
st.plotly_chart(
    line_metric_figure(
        labels,
        [{"name": "boundary energy ratio", "y": metrics["boundary_energy_ratio"], "color": PAPER_COLORS["schemeB"]}],
        "Boundary-layer dominance",
        "bridge",
        "ratio",
    ),
    use_container_width=True,
)

sets = bridge["interior_sets"]
st.info(
    f"Interior indices: {sets['interior_idx'] or 'empty'} | Boundary indices: {sets['boundary_idx']} | r(k) = {sets['r_k']}"
)

st.subheader("结论")
actual_full = bridge["actual_vs_ideal"]["full_l2"]
actual_interior = bridge["actual_vs_ideal"]["interior_l2"]
ref_gap = bridge["reference_vs_ideal"]["full_l2"]
raw_ref_gap = bridge["actual_vs_reference"]["full_l2"]
st.markdown(
    f"""
- 当前 case 下，`actual -> ideal` 的 full error 是 `{actual_full:.6e}`，其中 interior error 是 `{actual_interior:.6e}`。
- `actual -> reference` 的 gap 是 `{raw_ref_gap:.6e}`，`reference -> ideal` 的 gap 是 `{ref_gap:.6e}`。
- 因此这页主要帮助判断主桥梁问题更偏向 `raw -> reference-window`，还是更偏向 `reference-window -> ideal truncation`。
"""
)
