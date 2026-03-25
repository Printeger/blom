from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import streamlit as st

DEMO_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = DEMO_ROOT.parent
for candidate in (str(REPO_ROOT), str(DEMO_ROOT)):
    if candidate not in sys.path:
        sys.path.insert(0, candidate)

from demo_core.benchmark_tools import load_phase10_benchmark_artifacts
from demo_core.plotting import PAPER_COLORS, dataframe_download, line_metric_figure


st.title("Benchmark / Ablation Summary")
st.caption("适合组会和论文补充材料的结果汇总页，优先复用 Phase 10 已生成的真实结果文件。")

artifacts = load_phase10_benchmark_artifacts()
benchmark_rows = artifacts["benchmark_rows"]
ablation = artifacts["ablation"]
suite = artifacts["suite"]
benchmark_json = artifacts["benchmark"]

with st.expander("参数说明 / 图表说明"):
    st.markdown(
        """
- `Runtime vs M`: 比较不同方法在规模 `M` 增大时的总运行时间。
- `Objective gap vs k`: 比较 raw Scheme C 随 `k` 变化时的目标差异。
- `Support width / memory vs k`: 展示 locality 与资源占用的关系。
- `Term ablation`: 依次去掉 `obs / dyn / time / reg` 项后的目标变化。
- 本页读取的是 Phase 10 已生成的真实 benchmark 结果，不是前端临时重算。
"""
    )

with st.expander("名词解释"):
    st.markdown(
        """
- `objective gap`: 相对当前比较组里最优目标值的差距。
- `support width mean`: 参数对 segment blocks 的平均影响宽度。
- `memory peak`: 运行时峰值内存。
- `baseline compare`: raw Scheme C、MINCO、Scheme A 和 heuristic 的横向对比。
- `term ablation`: 去掉某个目标项后，整体目标与行为的变化。
"""
    )

if benchmark_rows is None or ablation is None or suite is None or benchmark_json is None:
    st.error("Phase 10 benchmark artifacts were not found. Please run `python3 -m phase_10.blom_phase10_framework_suite` first.")
else:
    m_values = set(benchmark_json.get("config", {}).get("M_values", []))
    m_rows = benchmark_rows[benchmark_rows["M"].isin(m_values)].copy()
    k_rows = benchmark_rows[(benchmark_rows["method"] == "raw_schemeC") & (~benchmark_rows["M"].isin(m_values))].copy()
    if k_rows.empty:
        k_rows = benchmark_rows[benchmark_rows["method"] == "raw_schemeC"].copy()

    st.subheader("Runtime vs M")
    runtime_fig = line_metric_figure(
        sorted(m_rows["M"].unique().tolist()),
        [
            {"name": method, "y": m_rows[m_rows["method"] == method].sort_values("M")["runtime_total"].tolist(), "color": color}
            for method, color in [
                ("raw_schemeC", PAPER_COLORS["blom"]),
                ("minco", PAPER_COLORS["minco"]),
                ("schemeA", PAPER_COLORS["schemeA"]),
                ("heuristic", PAPER_COLORS["bspline"]),
            ]
        ],
        "Runtime vs M",
        "M",
        "runtime total (sec)",
    )
    st.plotly_chart(runtime_fig, use_container_width=True)

    st.subheader("Objective gap / support width vs k")
    cols = st.columns(2)
    cols[0].plotly_chart(
        line_metric_figure(
            k_rows["k"].tolist(),
            [{"name": "objective gap", "y": k_rows["objective_gap"].tolist(), "color": PAPER_COLORS["minco"]}],
            "Objective gap vs k",
            "k",
            "objective gap",
        ),
        use_container_width=True,
    )
    cols[1].plotly_chart(
        line_metric_figure(
            k_rows["k"].tolist(),
            [
                {"name": "support width mean", "y": k_rows["support_width_mean"].tolist(), "color": PAPER_COLORS["schemeC"]},
                {"name": "memory peak", "y": k_rows["memory_peak"].tolist(), "color": PAPER_COLORS["diff"]},
            ],
            "Support width / memory vs k",
            "k",
            "value",
        ),
        use_container_width=True,
    )

    st.subheader("Baseline compare")
    st.dataframe(m_rows[["method", "M", "final_objective", "runtime_total", "support_width_mean"]], use_container_width=True, hide_index=True)

    st.subheader("Term ablation")
    st.dataframe(ablation, use_container_width=True, hide_index=True)

    st.subheader("Downloads")
    st.download_button("Download benchmark rows CSV", dataframe_download(benchmark_rows), file_name="phase10_benchmark_rows.csv", mime="text/csv")
    st.download_button("Download ablation CSV", dataframe_download(ablation), file_name="phase10_ablation_summary.csv", mime="text/csv")

    st.info(
        f"Suite summary: optimizer drop = {suite['optimizer']['objective_drop']:.6e}, backward dense-sparse gap = {suite['backward']['dense_vs_sparse_max_abs']:.3e}."
    )

    best_method_row = m_rows.loc[m_rows["final_objective"].idxmin()]
    fastest_row = m_rows.loc[m_rows["runtime_total"].idxmin()]
    st.subheader("结论")
    st.markdown(
        f"""
- 在当前 benchmark 明细中，目标值最优的方法是 `{best_method_row['method']}`，其代表行 objective 为 `{best_method_row['final_objective']:.6e}`。
- 运行最快的方法是 `{fastest_row['method']}`，其 runtime 为 `{fastest_row['runtime_total']:.6e} sec`。
- 因此这页的主要结论是：Phase 10 已经可以把质量、运行时间和局部性指标放到同一框架里讨论，而不是只盯单一目标值。
"""
    )
