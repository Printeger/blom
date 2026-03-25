from __future__ import annotations

import numpy as np
import sys
from pathlib import Path
import streamlit as st

DEMO_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = DEMO_ROOT.parent
for candidate in (str(REPO_ROOT), str(DEMO_ROOT)):
    if candidate not in sys.path:
        sys.path.insert(0, candidate)

from demo_core.blom_adapter import assemble_scheme_compare, compute_scheme_fd_jacobian
from demo_core.data_interface import case_for_demo
from demo_core.jacobian_tools import summarize_dense_jacobian
from demo_core.plotting import PAPER_COLORS, bar_metric_figure
from phase_7.blom_convergence_vs_k import compute_ideal_truncated_blom_k


st.title("Scheme A / B / C Compare")
st.caption("展示 continuity-vs-locality trade-off。raw Scheme C 仍是主对象，A/B 作为对照。")

col1, col2 = st.columns(2)
M = col1.slider("M", 4, 8, 6, 1, help="Number of segments.")
regime = col2.selectbox("regime", ["bounded-nonuniform", "uniform"], index=0, help="Duration regime.")

with st.expander("参数说明 / 名词解释 / 图表说明"):
    st.markdown(
        """
- `Scheme A`: 共享全局 junction states，连续性强，但局部性较弱。
- `Scheme B`: 局部预测后再做 consensus，介于 A 和 C 之间。
- `Scheme C`: 直接取局部窗口中心段，是当前 demo 的主 BLOM 对象。
- `lower-order jumps`: 低阶导数在拼接点的跳变，越小表示连续性越强。
- `higher-order jumps`: 更高阶导数的跳变。
- `outside-band sensitivity`: 超出预期局部带宽之外的灵敏度能量。
- `matching vs ideal truncation`: 与 ideal truncated BLOM-k 的系数差。
"""
    )

q, T = case_for_demo(M=M, regime=regime)
schemes = assemble_scheme_compare(q, T, k=2)
ideal = compute_ideal_truncated_blom_k(q, T, k=2)

lower_jump = []
higher_jump = []
support_width = []
outside_band = []
matching = []
labels = ["A", "B", "C"]

for scheme in labels:
    result = schemes[scheme]
    stats = result["stats"]
    lower_jump.append(max(float(stats[order]["max_abs"]) for order in range(4)))
    higher_jump.append(max(float(stats[order]["max_abs"]) for order in range(4, 7)))
    fd = compute_scheme_fd_jacobian(tuple(q.tolist()), tuple(T.tolist()), scheme)
    jac_summary = summarize_dense_jacobian(fd["J_c_q"], block_rows=8, tol=1e-8, expected_band=2)
    support_width.append(jac_summary["mean_width"])
    outside_band.append(jac_summary["outside_band_energy"])
    matching.append(float(np.linalg.norm(np.asarray(result["coeffs"], dtype=float).reshape(-1) - np.asarray(ideal["c_vec"], dtype=float))))

cols = st.columns(2)
cols[0].plotly_chart(bar_metric_figure(labels, lower_jump, "Lower-order jumps", "max abs jump", color=PAPER_COLORS["schemeB"]), use_container_width=True)
cols[1].plotly_chart(bar_metric_figure(labels, higher_jump, "Higher-order jumps", "max abs jump", color=PAPER_COLORS["schemeA"]), use_container_width=True)
cols = st.columns(3)
cols[0].plotly_chart(bar_metric_figure(labels, support_width, "Support width", "mean support width", color=PAPER_COLORS["schemeC"]), use_container_width=True)
cols[1].plotly_chart(bar_metric_figure(labels, outside_band, "Outside-band sensitivity", "energy ratio", color=PAPER_COLORS["diff"]), use_container_width=True)
cols[2].plotly_chart(bar_metric_figure(labels, matching, "Matching vs ideal truncation", "l2 error", color=PAPER_COLORS["minco"]), use_container_width=True)

st.warning("The locality comparison on this page is based on finite-difference assembly Jacobians and is therefore an empirical result page, not a theorem page.")

best_locality = labels[int(np.argmin(outside_band))]
best_matching = labels[int(np.argmin(matching))]
best_continuity = labels[int(np.argmin(lower_jump))]
st.subheader("结论")
st.markdown(
    f"""
- 当前 case 下，连续性最好的是 `{best_continuity}`，因为它的 lower-order jumps 最小。
- 局部性最好的是 `{best_locality}`，因为它的 outside-band sensitivity 最小。
- 与 ideal truncation 最接近的是 `{best_matching}`。
- 这页的核心结论是：A/B/C 不是简单的优劣排序，而是在连续性、局部性和 matching 之间做不同取舍。
"""
)
