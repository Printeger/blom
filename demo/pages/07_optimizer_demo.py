from __future__ import annotations

import numpy as np
import sys
from pathlib import Path
import pandas as pd
import streamlit as st

DEMO_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = DEMO_ROOT.parent
for candidate in (str(REPO_ROOT), str(DEMO_ROOT)):
    if candidate not in sys.path:
        sys.path.insert(0, candidate)

from demo_core.blom_adapter import run_optimizer_cached
from demo_core.data_interface import case_for_demo
from demo_core.plotting import PAPER_COLORS, trajectory_figure


st.title("Phase 9 / 10 Optimizer Demo")
st.caption("展示 BLOM 已经是一个可优化对象，而不是静态轨迹表示。")

row = st.columns(5)
M = row[0].slider("M", 4, 8, 6, 1, help="Number of segments in the optimization case.")
regime = row[1].selectbox("regime", ["bounded-nonuniform", "uniform"], index=0, help="Initial duration regime.")
objective_mode = row[2].selectbox("objective", ["full", "minimal"], index=0, help="`minimal` uses the Phase 9 objective; `full` uses the Phase 10 objective.")
n_steps = row[3].slider("n_steps", 5, 30, 12, 1, help="Number of optimization iterations to run.")
step_size = row[4].slider("step_size", 1e-3, 1e-2, 5e-3, format="%.3e", help="Gradient-descent base step size.")
k = st.slider("k (full objective only)", 2, 8, 2, 2, help="Local BLOM window size for the full objective path.")

with st.expander("参数说明 / 图表说明"):
    st.markdown(
        """
- `objective = minimal` 对应 Phase 9 的最小可微闭环。
- `objective = full` 对应 Phase 10 的完整目标层。
- 第一张轨迹 before / after 图使用归一化进度轴 `t / total_time`，这样优化前后总时长不同也不会把 `before` 曲线压得很短。
- 第二张轨迹 before / after 图使用真实时间轴，用来展示绝对时间分配变化后曲线是如何被拉伸或压缩的。
- 时间分配变化单独用柱状图展示，避免和轨迹形状混在一起。
- `dense vs sparse gap` 只在 full objective 下显示，用来说明 backward consistency。
"""
    )

with st.expander("名词解释"):
    st.markdown(
        """
- `trajectory value` / `p(t)`: 标量轨迹在时刻 `t` 的函数值，也就是这条 1D 轨迹本身的纵坐标。
- `ctrl`: 控制代价，主要对应高阶导数或 snap 相关代价。
- `time`: 时间项，对应总时长或时间分配代价。
- `obs`: 障碍相关软惩罚项。
- `dyn`: 动力学约束相关软惩罚项。
- `bc`: 边界条件相关惩罚项。
- `reg`: 规则化项，通常用于时间平滑等。
- `dense vs sparse gap`: dense backward checker 与 sparse backward 实现的梯度差。
"""
    )

q, T = case_for_demo(M=M, regime=regime)
payload = run_optimizer_cached(tuple(q.tolist()), tuple(T.tolist()), objective_mode, n_steps, float(step_size), k)

before = payload["before"]
after = payload["after"]
before_x = before["t"] / max(float(before["t"][-1]), 1e-12) if len(before["t"]) else before["t"]
after_x = after["t"] / max(float(after["t"][-1]), 1e-12) if len(after["t"]) else after["t"]
fig = trajectory_figure(
    [
        {"x": before_x, "y": before["y"], "name": "before", "color": PAPER_COLORS["diff"]},
        {"x": after_x, "y": after["y"], "name": "after", "color": PAPER_COLORS["blom"]},
    ],
    f"{objective_mode} objective: trajectory before / after (normalized progress)",
    x_label="normalized progress",
    y_label="trajectory value p(t)",
)
st.plotly_chart(fig, use_container_width=True)

fig_abs = trajectory_figure(
    [
        {"x": before["t"], "y": before["y"], "name": "before", "color": PAPER_COLORS["diff"]},
        {"x": after["t"], "y": after["y"], "name": "after", "color": PAPER_COLORS["blom"]},
    ],
    f"{objective_mode} objective: trajectory before / after (real time axis)",
    x_label="real time",
    y_label="trajectory value p(t)",
)
st.plotly_chart(fig_abs, use_container_width=True)

st.info(
    "两张 before / after 图展示的是同一个任务实例在优化前后的同一条标量轨迹。第一张强调形状变化，第二张强调真实时间分配变化。"
)

if objective_mode == "minimal":
    history = payload["result"]["history"]
    obj = history["objective"]
    parts = {"ctrl": history["ctrl"], "time": history["time"], "obs": history["obs"]}
    st.line_chart({"objective": obj, **parts})
    st.metric("objective drop", f"{payload['result']['summary']['objective_drop']:.6e}")
    T_before = np.asarray(T, dtype=float)
    T_after = np.asarray(payload["result"]["T_history"][-1], dtype=float)
else:
    history = payload["result"]["history"]
    st.line_chart(
        {
            "objective": history["objective"],
            "ctrl": history["ctrl"],
            "time": history["time"],
            "obs": history["obs"],
            "dyn": history["dyn"],
            "bc": history["bc"],
            "reg": history["reg"],
        }
    )
    st.metric("objective drop", f"{payload['result']['summary']['objective_drop']:.6e}")
    dense = payload["dense"]["grad_theta"]
    sparse = payload["sparse"]["grad_theta"]
    st.metric("dense vs sparse gap", f"{np.max(np.abs(dense - sparse)):.3e}")
    T_before = np.asarray(T, dtype=float)
    T_after = np.asarray(payload["result"]["T_final"], dtype=float)

time_df = pd.DataFrame({"before": T_before, "after": T_after})
st.subheader("Time Allocation Before / After")
st.bar_chart(time_df)
col_a, col_b = st.columns(2)
col_a.metric("total time before", f"{T_before.sum():.4f}")
col_b.metric("total time after", f"{T_after.sum():.4f}")

st.info("This page uses true Phase 9 / Phase 10 optimizers with Streamlit caching to keep interaction responsive.")

st.subheader("结论")
if objective_mode == "minimal":
    summary = payload["result"]["summary"]
    st.markdown(
        f"""
- 当前 minimal objective 下，objective drop 为 `{summary['objective_drop']:.6e}`。
- 这说明 Phase 9 的最小可微闭环在当前 case 上可以形成稳定下降，而不是静态表示。
"""
    )
else:
    summary = payload["result"]["summary"]
    gap = float(np.max(np.abs(payload["dense"]["grad_theta"] - payload["sparse"]["grad_theta"])))
    st.markdown(
        f"""
- 当前 full objective 下，objective drop 为 `{summary['objective_drop']:.6e}`。
- dense vs sparse backward gap 为 `{gap:.3e}`，说明前向优化和反向传播链路在当前 case 上是一致的。
- 这页的核心结论是：BLOM 在本仓库里已经不仅能表示轨迹，也能作为真实优化对象参与时空联合优化。
"""
    )
