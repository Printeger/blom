from __future__ import annotations

import time
import sys
from pathlib import Path

import streamlit as st

DEMO_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = DEMO_ROOT.parent
for candidate in (str(REPO_ROOT), str(DEMO_ROOT)):
    if candidate not in sys.path:
        sys.path.insert(0, candidate)

from demo_core.blom_adapter import compute_actual, compute_ideal, compute_raw_jacobians
from demo_core.data_interface import case_for_demo
from demo_core.jacobian_tools import summarize_dense_jacobian
from demo_core.plotting import PAPER_COLORS, line_metric_figure
from phase_7.blom_convergence_vs_k import compute_convergence_errors


st.title("k-sweep / Locality-Approximation Trade-off")
st.caption("比较 raw Scheme C 和 ideal truncation 随 k 变化时的 matching、support width 和 runtime。")

col1, col2, col3 = st.columns(3)
M = col1.slider("M", 4, 10, 6, 1, help="Number of segments.")
regime = col2.selectbox("regime", ["bounded-nonuniform", "uniform"], index=0, help="Duration regime.")
k_max = col3.slider("max k", 2, 8, 8, 2, help="Largest locality window size to compare.")

with st.expander("参数说明 / 名词解释 / 图表说明"):
    st.markdown(
        """
- `matching l2`: raw Scheme C 与 ideal truncation 之间的系数差。
- `ideal truncation tail norm`: ideal truncation 相对全局参考被截断掉的尾部量级。
- `support width`: Jacobian 在段级 block 上的平均影响宽度。
- `runtime`: 当前 demo 前端调用该对象时的实际运行时间。
- 第一张图更偏近似误差，第二张图更偏 locality / cost-of-locality。
"""
    )

q, T = case_for_demo(M=M, regime=regime)
k_values = list(range(2, k_max + 1, 2))
matching = []
support = []
runtime = []
ideal_tail = []

for k in k_values:
    start = time.perf_counter()
    actual = compute_actual(tuple(q.tolist()), tuple(T.tolist()), k=k, scheme="C")
    runtime.append(time.perf_counter() - start)
    ideal = compute_ideal(tuple(q.tolist()), tuple(T.tolist()), k=k)
    metrics = compute_convergence_errors(ideal["c_vec"], actual["c_vec"], c_ideal=ideal["c_vec"], cost_ref=ideal["cost"], cost_test=actual["cost"])
    matching.append(float(metrics["matching_l2"]))
    ideal_tail.append(float(ideal["tail_norm"]))
    jac = compute_raw_jacobians(tuple(q.tolist()), tuple(T.tolist()), k=k)
    support.append(summarize_dense_jacobian(jac["J_c_q_dense"], block_rows=8, tol=1e-8, expected_band=k)["mean_width"])

st.plotly_chart(
    line_metric_figure(
        k_values,
        [
            {"name": "raw-vs-ideal matching l2", "y": matching, "color": PAPER_COLORS["blom"]},
            {"name": "ideal truncation tail norm", "y": ideal_tail, "color": PAPER_COLORS["minco"]},
        ],
        "Approximation metrics vs k",
        "k",
        "metric",
    ),
    use_container_width=True,
)

st.subheader("结论")
if len(k_values) >= 2:
    trend_matching = "下降" if matching[-1] <= matching[0] else "上升"
    trend_support = "上升" if support[-1] >= support[0] else "下降"
    st.markdown(
        f"""
- 在当前 case 下，raw-vs-ideal matching 从 `k={k_values[0]}` 到 `k={k_values[-1]}` 整体呈 `{trend_matching}` 趋势。
- support width 随 `k` 整体呈 `{trend_support}` 趋势，这体现了 locality 与近似质量之间的基本 trade-off。
- 因此这页的核心结论不是“`k` 越大一定越好”，而是：更大的局部窗口通常意味着更宽的依赖范围和更高的计算开销。
"""
    )
st.plotly_chart(
    line_metric_figure(
        k_values,
        [
            {"name": "support width", "y": support, "color": PAPER_COLORS["schemeC"]},
            {"name": "runtime (sec)", "y": runtime, "color": PAPER_COLORS["diff"]},
        ],
        "Locality / runtime vs k",
        "k",
        "value",
    ),
    use_container_width=True,
)
