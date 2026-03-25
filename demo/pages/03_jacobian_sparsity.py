from __future__ import annotations

import sys
from pathlib import Path

import streamlit as st

DEMO_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = DEMO_ROOT.parent
for candidate in (str(REPO_ROOT), str(DEMO_ROOT)):
    if candidate not in sys.path:
        sys.path.insert(0, candidate)

from demo_core.blom_adapter import compute_raw_jacobians
from demo_core.data_interface import case_for_demo
from demo_core.jacobian_tools import summarize_dense_jacobian, threshold_mask
from demo_core.minco_adapter import compute_minco_data, compute_minco_fd_jacobian_T
from demo_core.plotting import binary_mask_heatmap, matrix_heatmap


st.title("Jacobian Sparsity / Support Pattern")
st.caption("这是 BLOM 局部性的主证据页之一，优先展示真实 Jacobian、support width 和 outside-band 指标。")

col1, col2, col3, col4 = st.columns(4)
M = col1.slider("M", 4, 10, 6, 1, help="Number of segments.")
regime = col2.selectbox("regime", ["bounded-nonuniform", "uniform"], index=0, help="Time-allocation regime.")
k = col3.slider("k", 2, 8, 2, 2, help="BLOM locality window size.")
threshold = col4.slider("threshold", 1e-12, 1e-2, 1e-8, format="%.1e", help="Values below this threshold are treated as visually zero in the sparsity mask.")

obj = st.radio("object", ["raw Scheme C", "MINCO"], horizontal=True)
jac_kind = st.radio("Jacobian", ["d c / d q_bar", "d c / d T"], horizontal=True)

with st.expander("参数说明 / 名词解释 / 图表说明"):
    st.markdown(
        """
- `d c / d q_bar`: 系数向量对 interior waypoints 的 Jacobian。
- `d c / d T`: 系数向量对时间分配的 Jacobian。
- `support width`: 单个参数能影响多少个 segment block。
- `outside-band energy`: 超出预期局部带宽之外的灵敏度能量占比。
- 上面的热图显示 `log10(|Jacobian entry| + 1e-16)`。
- 下面的 mask 图只显示“大于 threshold 的项是否为非零”。
"""
    )

q, T = case_for_demo(M=M, regime=regime)

if obj == "raw Scheme C":
    data = compute_raw_jacobians(tuple(q.tolist()), tuple(T.tolist()), k=k)
    J = data["J_c_q_dense"] if jac_kind == "d c / d q_bar" else data["J_c_T_dense"]
    expected_band = k if jac_kind == "d c / d q_bar" else max(1, k // 2 + 1)
else:
    minco = compute_minco_data(q, T)
    J = minco["J_q"] if jac_kind == "d c / d q_bar" else compute_minco_fd_jacobian_T(q, T)
    expected_band = T.size

stats = summarize_dense_jacobian(J, block_rows=8, tol=threshold, expected_band=expected_band)

colA, colB, colC = st.columns(3)
colA.metric("mean support width", f"{stats['mean_width']:.3f}")
colB.metric("max support width", f"{stats['max_width']:.3f}")
colC.metric("outside-band energy", f"{stats.get('outside_band_energy', 0.0):.3e}")

st.plotly_chart(matrix_heatmap(J, f"{obj}: {jac_kind} absolute heatmap", "parameter index", "coefficient row"), use_container_width=True)
st.plotly_chart(binary_mask_heatmap(threshold_mask(J, threshold), f"{obj}: thresholded sparsity mask", "parameter index", "coefficient row"), use_container_width=True)

if obj == "MINCO":
    st.warning("For MINCO, `d c / d q_bar` uses the exact Phase 7 kernel; `d c / d T` is shown via finite differences because the repository does not expose an exact dense T-kernel API.")

st.subheader("结论")
if obj == "raw Scheme C":
    st.markdown(
        f"""
- 当前设置下，`{jac_kind}` 的 mean support width 是 `{stats['mean_width']:.3f}`，max support width 是 `{stats['max_width']:.3f}`。
- outside-band energy 为 `{stats.get('outside_band_energy', 0.0):.3e}`，这页正是 BLOM 局部性的主要数值证据之一。
- 这里如果看到明显的 block-banded 结构，说明局部性并不是只靠曲线外观“看起来像局部”，而是体现在真实 Jacobian 上。
"""
    )
else:
    st.markdown(
        f"""
- 当前设置下，MINCO 的 `{jac_kind}` support width 明显更宽，mean support width 为 `{stats['mean_width']:.3f}`。
- 这说明 MINCO 更适合作为全局 reference，而不是局部支持的主对象。
"""
    )
