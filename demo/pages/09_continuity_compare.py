from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

DEMO_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = DEMO_ROOT.parent
for candidate in (str(REPO_ROOT), str(DEMO_ROOT)):
    if candidate not in sys.path:
        sys.path.insert(0, candidate)

from demo_core.data_interface import case_for_demo
from demo_core.plotting import PAPER_COLORS
from phase_1.minco_scalar_baseline import continuity_jumps
from phase_5.blom_boundary_jump_check import assemble_scheme_A, assemble_scheme_B, assemble_scheme_C, compute_jumps, summarize_jumps
from phase_7.blom_convergence_vs_k import compute_minco_reference


def _max_abs_per_order_from_minco(coeffs: np.ndarray, T: np.ndarray, max_order: int = 4) -> dict[int, float]:
    jumps = continuity_jumps(coeffs, T, max_order=max_order)
    return {order: float(np.max(np.abs(jumps[:, order]))) if jumps.size else 0.0 for order in range(max_order + 1)}


def _max_abs_per_order_from_scheme(coeffs: np.ndarray, T: np.ndarray, max_order: int = 4) -> dict[int, float]:
    jumps = compute_jumps(coeffs, T, s=4, max_order=max_order)
    summary = summarize_jumps(jumps)
    return {order: float(summary[order]["max_abs"]) for order in range(max_order + 1)}


def _guarantee_map() -> dict[str, dict[str, bool]]:
    return {
        "MINCO": {"C0": True, "C1": True, "C2": True, "C3": True, "C4": True},
        "Scheme A": {"C0": True, "C1": True, "C2": True, "C3": True, "C4": False},
        "Scheme B": {"C0": True, "C1": True, "C2": True, "C3": True, "C4": False},
        "raw Scheme C": {"C0": True, "C1": False, "C2": False, "C3": False, "C4": False},
        "B-spline (cubic)": {"C0": True, "C1": True, "C2": True, "C3": False, "C4": False},
    }


st.title("Continuity Compare")
st.caption("比较不同对象在 `C0` 到 `C4` 上的连续性：既看理论保证，也看当前数值 case 的实际 jump。")

row = st.columns(4)
M = row[0].slider("M", 4, 10, 6, 1, help="Number of segments.")
regime = row[1].selectbox("regime", ["bounded-nonuniform", "uniform"], index=0, help="Duration regime.")
k = row[2].slider("BLOM k", 2, 8, 2, 2, help="Used for Scheme A/B/C assembly.")
tol = row[3].slider("zero tolerance", 1e-12, 1e-4, 1e-8, format="%.1e", help="Threshold for labeling a jump as numerically zero.")

with st.expander("参数说明 / 名词解释 / 图表说明"):
    st.markdown(
        """
- `C0` 连续表示函数值连续，也就是轨迹没有位置跳变。
- `C1` 连续表示一阶导数连续，也就是速度连续。
- `C2` 连续表示二阶导数连续，也就是加速度连续。
- `C3` 连续表示三阶导数连续，也就是 jerk 连续。
- `C4` 连续表示四阶导数连续。
- 上面的 `zero tolerance` 用来判断“数值上是否可视为 0”，不是理论定义。
- 本页分两部分：
  - `理论保证表`：回答“能否保证”
  - `数值 jump 图`：回答“在当前 case 上观察到什么”
"""
    )

q, T = case_for_demo(M=M, regime=regime)
minco = compute_minco_reference(q, T)
scheme_a = assemble_scheme_A(q, T, k=k)
scheme_b = assemble_scheme_B(q, T, k=k)
scheme_c = assemble_scheme_C(q, T, k=k)

observed = {
    "MINCO": _max_abs_per_order_from_minco(np.asarray(minco["coeffs"], dtype=float), T, max_order=4),
    "Scheme A": _max_abs_per_order_from_scheme(np.asarray(scheme_a["coeffs"], dtype=float), T, max_order=4),
    "Scheme B": _max_abs_per_order_from_scheme(np.asarray(scheme_b["coeffs"], dtype=float), T, max_order=4),
    "raw Scheme C": _max_abs_per_order_from_scheme(np.asarray(scheme_c["coeffs"], dtype=float), T, max_order=4),
}

guarantee = _guarantee_map()
guarantee_rows = []
for method, states in guarantee.items():
    row_dict = {"method": method}
    row_dict.update({key: ("guaranteed" if val else "not guaranteed") for key, val in states.items()})
    guarantee_rows.append(row_dict)

st.subheader("Theoretical Guarantee Table")
st.dataframe(pd.DataFrame(guarantee_rows), use_container_width=True, hide_index=True)

orders = [f"C{order}" for order in range(5)]
fig = go.Figure()
for method, color in [
    ("MINCO", PAPER_COLORS["minco"]),
    ("Scheme A", PAPER_COLORS["schemeA"]),
    ("Scheme B", PAPER_COLORS["schemeB"]),
    ("raw Scheme C", PAPER_COLORS["schemeC"]),
]:
    fig.add_trace(
        go.Bar(
            name=method,
            x=orders,
            y=[observed[method][order] for order in range(5)],
            marker_color=color,
        )
    )
fig.update_layout(
    title="Observed max abs jump by continuity order",
    template="plotly_white",
    barmode="group",
    xaxis_title="continuity order",
    yaxis_title="max abs jump",
    yaxis_type="log",
)
st.plotly_chart(fig, use_container_width=True)

heatmap_methods = ["MINCO", "Scheme A", "Scheme B", "raw Scheme C"]
heatmap = np.asarray([[observed[m][order] for order in range(5)] for m in heatmap_methods], dtype=float)
heatmap_fig = go.Figure(
    data=go.Heatmap(
        z=np.log10(heatmap + 1e-16),
        x=orders,
        y=heatmap_methods,
        colorscale="Viridis",
        colorbar=dict(title="log10(max jump + 1e-16)"),
    )
)
heatmap_fig.update_layout(
    title="Observed continuity heatmap on the current case",
    template="plotly_white",
    xaxis_title="continuity order",
    yaxis_title="method",
)
st.plotly_chart(heatmap_fig, use_container_width=True)

numerical_rows = []
for method in heatmap_methods:
    record = {"method": method}
    for order in range(5):
        value = observed[method][order]
        record[f"C{order}_jump"] = value
        record[f"C{order}_numerically_zero"] = bool(value <= tol)
    numerical_rows.append(record)

st.subheader("Numerical Zero Check")
st.dataframe(pd.DataFrame(numerical_rows), use_container_width=True, hide_index=True)

st.subheader("结论")
st.markdown(
    f"""
- 在理论层面，MINCO 在当前比较的 `C0` 到 `C4` 上都可以保证连续；Scheme A / B 主要保证到 `C3`；raw Scheme C 通常只保证 `C0`；cubic B-spline 通常保证到 `C2`。
- 在当前数值 case 中，MINCO 的 `C0` 到 `C4` jump 都非常小，说明它确实表现为高阶连续对象。
- Scheme A / B 在 `C0` 到 `C3` 上通常会接近数值零，但 `C4` 不应被解释为理论保证。
- raw Scheme C 如果 `C1` 到 `C4` jump 明显大于容差 `{tol:.1e}`，这正说明它的优势主要是局部性，而不是高阶全局连续性。
"""
    )
