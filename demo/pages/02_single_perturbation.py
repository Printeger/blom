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

from demo_core.blom_adapter import compute_actual, sample_actual
from demo_core.bspline_adapter import default_control_points_from_waypoints, evaluate_bspline, perturb_control_point
from demo_core.data_interface import case_for_demo
from demo_core.minco_adapter import compute_minco_data, sample_minco_trajectory
from demo_core.plotting import PAPER_COLORS, difference_heatmap, trajectory_figure


st.title("Single Perturbation Response")
st.caption("几何直观页。主要用于帮助理解，不是 BLOM 局部性的主证据页。")

row1 = st.columns(5)
M = row1[0].slider("M", 4, 12, 6, 1, help="Number of segments / intervals.")
regime = row1[1].selectbox("regime", ["bounded-nonuniform", "uniform"], index=0, help="Time-allocation regime for BLOM/MINCO.")
method = row1[2].selectbox("method", ["BLOM", "MINCO", "B-spline"], index=0, help="Choose which object to perturb.")
k = row1[3].slider("BLOM k", 2, 8, 2, 2, help="Only used for BLOM.")
delta = row1[4].slider("delta", -0.8, 0.8, 0.25, 0.05, help="Signed perturbation magnitude.")

with st.expander("参数说明 / 图表说明"):
    st.markdown(
        """
- 这页主要是几何直观页，不是 locality 的主证据页。
- `B-spline` 扰动的是 control point。
- `MINCO / BLOM` 扰动的是 interpolation waypoint。
- 热图显示的是 `perturbed - original` 的函数值差，而不是 Jacobian 证据。
"""
    )

with st.expander("名词解释"):
    st.markdown(
        """
- `control point`: B-spline 的控制点，通常不要求曲线经过它。
- `interpolation waypoint`: MINCO / BLOM 的插值点，轨迹需要穿过这些点或以这些点作为核心数据对象。
- `delta`: 单个控制变量的扰动幅度。
- `difference heatmap`: 展示扰动前后轨迹函数值的差异大小与符号。
"""
    )

q, T = case_for_demo(M=M, regime=regime)
if method == "B-spline":
    max_index = q.size - 1
    index = st.slider("control-point index", 0, max_index, min(max_index, max_index // 2))
    ctrl = default_control_points_from_waypoints(q)
    ctrl_p = perturb_control_point(ctrl, index, delta)
    base = evaluate_bspline(ctrl)
    perturbed = evaluate_bspline(ctrl_p)
    total_time = float(T.sum())
    x_base, y_base = base["u"] * total_time, base["y"]
    x_new, y_new = perturbed["u"] * total_time, perturbed["y"]
    note = "B-spline uses control-point perturbation, not interpolation waypoint perturbation."
else:
    max_index = q.size - 2
    index = st.slider("interior waypoint index", 1, max_index, min(max_index, max_index // 2))
    q_p = q.copy()
    q_p[index] += delta
    if method == "BLOM":
        base = compute_actual(tuple(q.tolist()), tuple(T.tolist()), k=k, scheme="C")
        perturbed = compute_actual(tuple(q_p.tolist()), tuple(T.tolist()), k=k, scheme="C")
        base_curve = sample_actual(base["coeffs"], T)
        pert_curve = sample_actual(perturbed["coeffs"], T)
        x_base, y_base = base_curve["t"], base_curve["y"]
        x_new, y_new = pert_curve["t"], pert_curve["y"]
    else:
        base = compute_minco_data(q, T)
        perturbed = compute_minco_data(q_p, T)
        base_curve = sample_minco_trajectory(base["coeffs"], T)
        pert_curve = sample_minco_trajectory(perturbed["coeffs"], T)
        x_base, y_base = base_curve["t"], base_curve["y"]
        x_new, y_new = pert_curve["t"], pert_curve["y"]
    note = "MINCO/BLOM use interpolation-waypoint perturbation."

curve_fig = trajectory_figure(
    [
        {"x": x_base, "y": y_base, "name": "original", "color": PAPER_COLORS["diff"]},
        {"x": x_new, "y": y_new, "name": "perturbed", "color": PAPER_COLORS["blom"] if method == "BLOM" else PAPER_COLORS["minco"] if method == "MINCO" else PAPER_COLORS["bspline"]},
    ],
    f"{method} perturbation response",
    x_label="aligned time axis",
    y_label="trajectory value",
)
st.plotly_chart(curve_fig, use_container_width=True)
st.plotly_chart(difference_heatmap(np.asarray(x_new), np.asarray(y_base), np.asarray(y_new), f"{method} response heatmap"), use_container_width=True)
st.info(note)

st.subheader("结论")
if method == "B-spline":
    st.markdown(
        """
- 这页说明 B-spline 的控制点扰动会产生局部几何变化，但它和 BLOM/MINCO 的 waypoint 扰动不是同类变量。
- 因此这页只能帮助建立几何直观，不能拿来单独证明 BLOM 的 Jacobian 局部性。
"""
    )
else:
    st.markdown(
        """
- 这页可以帮助观察单个 waypoint 扰动后的几何响应范围。
- 但 BLOM 的局部性结论仍然应以 Jacobian 稀疏性和 support-width 指标为主，而不是只靠这张曲线图。
"""
    )
