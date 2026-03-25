from __future__ import annotations

import sys
from pathlib import Path

import streamlit as st

DEMO_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = DEMO_ROOT.parent
for candidate in (str(REPO_ROOT), str(DEMO_ROOT)):
    if candidate not in sys.path:
        sys.path.insert(0, candidate)

from demo_core.blom_adapter import compute_actual, sample_actual
from demo_core.bspline_adapter import default_control_points_from_waypoints, evaluate_bspline
from demo_core.data_interface import case_for_demo, comparison_table_frame
from demo_core.minco_adapter import compute_minco_data, sample_minco_trajectory
from demo_core.plotting import PAPER_COLORS, trajectory_figure


st.title("Overview / 方法总览")
st.caption("结构对比页，不把视觉现象写成 theorem。")

col1, col2, col3 = st.columns(3)
M = col1.slider("M", min_value=4, max_value=12, value=6, step=1, help="Number of polynomial segments / waypoint intervals.")
regime = col2.selectbox("regime", ["bounded-nonuniform", "uniform"], index=0, help="`uniform` uses equal durations; `bounded-nonuniform` samples durations from a bounded interval.")
k = col3.slider("BLOM k", min_value=2, max_value=8, value=2, step=2, help="BLOM locality window size. Larger `k` means a wider local neighborhood.")

with st.expander("参数说明 / 图表说明"):
    st.markdown(
        """
- `M`: 轨迹段数，越大表示问题规模越大。
- `regime`: 时间配置；`uniform` 更规则，`bounded-nonuniform` 更接近一般情况。
- `k`: BLOM 的局部窗口尺度。

图中三条曲线都被放在同一个横轴上：
- MINCO / BLOM 使用真实全局时间轴；
- B-spline 的参数 `u in [0,1]` 会被线性映射到同样的总时长，只是为了可视化对齐；
- 这不意味着 B-spline 和 BLOM/MINCO 共享同一类变量。
"""
    )

with st.expander("名词解释"):
    st.markdown(
        """
- `raw Scheme C`: 当前 demo 里的主 BLOM 对象，来自局部窗口求解后提取中心段。
- `MINCO`: 全局 minimum-snap 插值参考解。
- `B-spline`: 这里作为局部几何控制 baseline，不代表与 BLOM 完全同类的优化对象。
- `local support`: 某个参数只影响局部系数块或局部曲线区域。
- `variational meaning`: 该对象是否直接对应一个明确的优化问题或极小化原理。
"""
    )

q, T = case_for_demo(M=M, regime=regime)
blom = compute_actual(tuple(q.tolist()), tuple(T.tolist()), k=k, scheme="C")
minco = compute_minco_data(q, T)
bspline = evaluate_bspline(default_control_points_from_waypoints(q))
total_time = float(T.sum())

st.subheader("Structure Comparison")
st.dataframe(comparison_table_frame(), use_container_width=True, hide_index=True)

st.subheader("Representative Trajectories")
blom_curve = sample_actual(blom["coeffs"], T)
minco_curve = sample_minco_trajectory(minco["coeffs"], T)
fig = trajectory_figure(
    [
        {"x": blom_curve["t"], "y": blom_curve["y"], "name": f"raw Scheme C (k={k})", "color": PAPER_COLORS["blom"]},
        {"x": minco_curve["t"], "y": minco_curve["y"], "name": "MINCO", "color": PAPER_COLORS["minco"]},
        {"x": bspline["u"] * total_time, "y": bspline["y"], "name": "B-spline control baseline", "color": PAPER_COLORS["bspline"]},
    ],
    "Representative geometry across the three objects",
    x_label="aligned time axis",
    y_label="trajectory value",
)
st.plotly_chart(fig, use_container_width=True)

st.info(
    "B-spline is shown as a control-polygon baseline. Its control points are not treated as BLOM/MINCO interpolation waypoints, and its horizontal axis is only time-aligned for visual comparison."
)

st.subheader("结论")
st.markdown(
    """
- 这一页的结论是“对象定位差异”而不是数值优劣。
- B-spline 更像局部几何控制基线；MINCO 是全局参考；raw Scheme C 介于两者之间，重点在局部性和可优化性。
- 如果后续要论证 BLOM 的优势，应主要看 Jacobian、support width、matching bridge 和 optimizer/benchmark 页面，而不是这页的几何外观。
"""
)
