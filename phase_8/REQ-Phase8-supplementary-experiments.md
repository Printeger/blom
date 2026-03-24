# REQ-Phase8-supplementary-experiments.md

## 1. 文档目的

本需求文档用于指导实现 **Phase 8 还需要补充的一整套验证实验**，用于判断当前 Phase 8 的结论是否足以从“方向性成立”提升到“更强的数值封口”。

当前已有 Phase 8 结果已经支持以下方向性结论：

1. interior-first matching 现象存在；
2. full matching gap 的主要能量集中在 boundary layer；
3. uniform-time 与 bounded-nonuniform 下都能观察到 interior-first 现象；
4. 下一条最值得优先补的理论桥更可能是：
   \[
   \text{reference-window / ideal-truncation consistency}
   \]
   而不是首先补 boundary attenuation。

但是，当前结果仍存在以下不足：

- \(M\) 偏小；
- `n_trials` 偏少；
- \(r(k)\) 选择固定为 `r(k)=k`，结论可能受定义影响；
- `raw vs reference-window \approx 0` 这一现象过于理想，需要更强 sanity check；
- interior error 在较大 \(k\) 下可能因为 interior set 为空而“看起来很好”。

因此，本需求文档的目标是：**一次性给出 Phase 8 还需要补做的补充实验集**，帮助判断：

- Phase 8 是否已经足以支撑论文中的强 observation；
- 是否还需要继续补更多数值支撑，才能写得更强；
- 后续 theorem 应优先补哪一个。

本需求文档默认兼容以下脚本：

- `blom_interior_matching_check.py`
- `blom_boundary_gap_decomposition.py`
- `blom_uniform_vs_nonuniform_interior.py`
- `blom_phase8_validation_suite.py`

---

## 2. 与前面 Phases 完全统一的符号要求

实现必须与 Phase 0--9 的符号保持一致。

### 2.1 基本对象

- 标量轨迹：
  \[
  p:[0,\mathcal T]\to\mathbb R
  \]
- 段数：
  \[
  M
  \]
- 航点：
  \[
  q_0,q_1,\dots,q_M\in\mathbb R
  \]
- interior waypoint vector：
  \[
  \bar q=(q_1,\dots,q_{M-1})^\top
  \]
- 时间向量：
  \[
  T=(T_1,\dots,T_M)^\top,\qquad T_i>0
  \]
- 总时长：
  \[
  \mathcal T=\sum_{i=1}^{M}T_i
  \]

### 2.2 canonical setting

本阶段仍只考虑 canonical matching 验证环境：

\[
s=4,\qquad k\in 2\mathbb N
\]

### 2.3 三个核心系数对象

- 全局 MINCO 参考系数：
  \[
  c^\star(q,T)
  \]
- idealized truncated BLOM-\(k\)：
  \[
  \tilde c^{(k)}(q,T)
  \]
- actual raw Scheme C：
  \[
  c^{\mathrm{C},k}(q,T)
  \]

### 2.4 matching error 记号

定义：

\[
e^{(k)}:=c^{\mathrm{C},k}-\tilde c^{(k)}
\]

以及其 block 误差：

\[
e_i^{(k)}:=c_i^{\mathrm{C},k}-\tilde c_i^{(k)}.
\]

full matching error：
\[
E_{\mathrm{full}}(k):=\|e^{(k)}\|_2
\]

interior-only matching error：
\[
E_{\mathrm{int}}(k):=\|P_{\mathrm{int}}^{(k)}e^{(k)}\|_2
\]

boundary-layer matching error：
\[
E_{\mathrm{bnd}}(k):=\|P_{\mathrm{bnd}}^{(k)}e^{(k)}\|_2
\]

---

## 3. 总体补充实验思路

Phase 8 当前还需要补的实验，建议一次性分成以下六类：

1. **更大 \(M\) sweep with nonempty interior**
2. **更多 random trials 的统计稳健性实验**
3. **boundary radius \(r(k)\) 敏感性实验**
4. **`raw vs reference-window` sanity check 实验**
5. **interior cardinality / empty-interior 假阳性排查实验**
6. **reference-window vs ideal-truncation 主导性验证实验**

这六类实验组合起来，能回答：

- 现有 Phase 8 结论是否稳健；
- 哪部分是定义效应，哪部分是真实结构；
- 下一步 theorem 到底该先补哪一条桥。

---

## 4. 实验 1：更大 \(M\) sweep with nonempty interior

## 4.1 目的

当前 `M=10` 且 `r(k)=k` 时，较大 \(k\) 下 interior set 可能为空。  
这会让 interior error 看起来“直接为 0”，从而夸大 interior-first 现象。

因此必须在更大 \(M\) 下重复 Phase 8 实验，使得：

- 在较大 \(k\) 下 interior set 仍然非空；
- 观察到的 interior improvement 不是由集合为空导致。

---

## 4.2 自变量设计

建议：

```python
M_values = [20, 40, 80, 120]
k_values = [2, 4, 6, 8, 10, 12]
```

要求：

- 对每个 \(M\) 和 \(k\)，必须记录 interior set 大小
- 若 interior set 为空，必须显式标注，不允许静默记为 0

---

## 4.3 必须计算的量

对每个 \((M,k)\)：

1. `interior_count`
2. `boundary_count`
3. `full_matching_l2`
4. `interior_matching_l2`
5. `boundary_matching_l2`
6. `boundary_energy_ratio`

---

## 4.4 必须保存的图片

### 图 1：不同 \(M\) 下 full / interior / boundary matching 曲线（必须）

文件名建议：

```text
phase8_largeM_matching_curves.png
```

### 图 2：interior cardinality vs \(k\)（必须）

文件名建议：

```text
phase8_interior_count_vs_k.png
```

### 图 3：boundary energy ratio vs \(M\) 与 \(k\)（推荐）

文件名建议：

```text
phase8_boundary_ratio_largeM.png
```

### 图 4：不同 \(M\) 下的 segmentwise matching heatmap（推荐）

文件名建议：

```text
phase8_largeM_heatmaps.png
```

---

## 4.5 实验结论目标

这一实验应回答：

- interior-first 现象是否在 interior 非空时仍然明显成立；
- 当前 strong interior result 是否不是由“小 \(M\) + 空 interior”造成的；
- boundary dominance 是否在更大 \(M\) 下仍稳定。

---

## 5. 实验 2：更多 random trials 的统计稳健性实验

## 5.1 目的

当前 Phase 8 的方向判断已成立，但 trial 数偏少。  
必须提高样本量，避免把偶然样本特性误当结构规律。

---

## 5.2 设计要求

建议：

```python
n_trials = 50 or 100
```

至少对两类 regime 都做：

- uniform-time
- bounded-nonuniform

每组 trial 要随机生成：

- `q`
- `T`（若 nonuniform）
- optional fixed seed list

---

## 5.3 必须统计的量

对每个 trial / each \(k\)：

1. full matching error
2. interior-only matching error
3. boundary-only matching error
4. interior log-fit slope
5. full log-fit slope
6. boundary energy ratio

再跨 trial 统计：

- mean
- median
- std
- q25 / q75

---

## 5.4 必须保存的图片

### 图 1：interior slope 分布箱线图（必须）

文件名建议：

```text
phase8_interior_slope_boxplot.png
```

### 图 2：boundary energy ratio 分布箱线图（必须）

文件名建议：

```text
phase8_boundary_ratio_boxplot.png
```

### 图 3：full vs interior slope 散点图（推荐）

文件名建议：

```text
phase8_full_vs_interior_slope_scatter.png
```

### 图 4：随机 trial 下 matching curve 的均值 ± 方差图（推荐）

文件名建议：

```text
phase8_trial_mean_curves.png
```

---

## 5.5 实验结论目标

这一实验应回答：

- interior-first 结论是否是稳定统计规律；
- boundary dominance 是否是普遍现象；
- current Phase 8 的结论是否足以从“单例 observation”提升到“稳健 empirical law”。

---

## 6. 实验 3：boundary radius \(r(k)\) 敏感性实验

## 6.1 目的

当前结论依赖 interior/boundary 的划分，而划分由 \(r(k)\) 决定。  
必须测试结论是否对 \(r(k)\) 的定义稳健。

---

## 6.2 建议的 radius family

至少测试以下三种：

### 模式 A
\[
r(k)=\frac{k}{2}
\]

### 模式 B
\[
r(k)=k
\]

### 模式 C
\[
r(k)=\left\lceil \frac{3k}{2}\right\rceil
\]

也可允许额外自定义：

```python
radius_modes = ["half_k", "k", "three_half_k"]
```

---

## 6.3 必须计算的量

对每个 radius mode：

1. interior cardinality
2. interior matching error
3. boundary matching error
4. boundary energy ratio
5. interior log-fit slope

---

## 6.4 必须保存的图片

### 图 1：不同 \(r(k)\) 下的 interior matching 曲线（必须）

文件名建议：

```text
phase8_radius_mode_interior_matching.png
```

### 图 2：不同 \(r(k)\) 下的 boundary ratio 曲线（必须）

文件名建议：

```text
phase8_radius_mode_boundary_ratio.png
```

### 图 3：不同 \(r(k)\) 下的 interior cardinality 曲线（必须）

文件名建议：

```text
phase8_radius_mode_interior_count.png
```

### 图 4：radius sensitivity summary heatmap（推荐）

文件名建议：

```text
phase8_radius_sensitivity_heatmap.png
```

---

## 6.5 实验结论目标

这一实验应回答：

- Phase 8 的 interior-first 结论是否对 \(r(k)\) 稳健；
- 现有 strongest observation 是否只是某个特定 radius 选择造成的；
- 后续 theorem 中 interior set 的定义是否必须保守选取。

---

## 7. 实验 4：`raw vs reference-window` sanity check

## 7.1 目的

当前结果中：

\[
\|c^{\mathrm{C},k}-\hat c^{(k)}\|
\approx 0
\]

这是非常强的信号，但也因为“太漂亮”而必须做额外 sanity check，避免：

- 实现路径共享过多导致的人为趋同；
- reference-window 构造本身存在隐式退化；
- 结论只在某个狭窄配置下成立。

---

## 7.2 必须设计的 sanity checks

### Check A：更大 \(M\) 与更多 trial 下重复
检查 `raw vs reference-window` 是否仍接近机器精度。

### Check B：uniform-time 与 bounded-nonuniform 分开看
检查该结论是否只在某一 regime 下成立。

### Check C：不同 \(k\) 下单独画曲线
检查是否所有 \(k\) 都接近 0，还是仅大 \(k\) 才接近 0。

### Check D：人工扰动 \(\gamma^\star\)
对 reference-window 人工加入一个小 perturbation：
\[
\gamma^\star \mapsto \gamma^\star+\delta\gamma
\]
看 central coefficient 是否有合理响应。

这个 check 的目的不是改理论，而是确认实现真的对 boundary trace 有感应，而不是路径上被“锁死”。

---

## 7.3 必须计算的量

1. raw vs reference-window error
   \[
   \|c^{\mathrm{C},k}-\hat c^{(k)}\|_2
   \]
2. perturbed reference-window sensitivity
   \[
   \|c^{\mathrm{ref},k}(\gamma^\star+\delta\gamma)-c^{\mathrm{ref},k}(\gamma^\star)\|_2
   \]
3. perturbation amplification ratio

---

## 7.4 必须保存的图片

### 图 1：raw vs reference-window error vs \(k\)（必须）

文件名建议：

```text
phase8_raw_vs_reference_error.png
```

### 图 2：different regimes raw-vs-reference boxplot（必须）

文件名建议：

```text
phase8_raw_vs_reference_boxplot.png
```

### 图 3：perturbed gamma sensitivity curve（必须）

文件名建议：

```text
phase8_reference_gamma_sensitivity.png
```

### 图 4：perturbed gamma scatter / log-log 图（推荐）

文件名建议：

```text
phase8_reference_gamma_loglog.png
```

---

## 7.5 实验结论目标

这一实验应回答：

- `raw vs reference-window ≈ 0` 是否是稳健现象；
- 该现象是否只是实现偶然；
- 这是否足以把理论重点从 boundary attenuation 转到 reference-window / ideal-truncation consistency。

---

## 8. 实验 5：empty-interior 假阳性排查实验

## 8.1 目的

专门排查以下风险：

> interior error 很小，是因为 interior 真的更好，还是因为 interior set 太小甚至为空？

这是 Phase 8 当前最需要防的解释风险之一。

---

## 8.2 实验设计

对每个 \((M,k,r(k))\)，必须记录：

- `interior_count`
- `interior_fraction`
- `is_empty_interior`

并对所有绘图和 summary 强制加上以下规则：

- 如果 interior set 为空，不允许把 interior error 直接当成“成功结果”
- 必须在 summary 中标注 `"empty interior"` 状态

---

## 8.3 必须保存的图片

### 图 1：interior fraction vs \(k\)（必须）

文件名建议：

```text
phase8_interior_fraction_vs_k.png
```

### 图 2：interior error with empty-set annotation（必须）

文件名建议：

```text
phase8_interior_error_with_empty_flags.png
```

### 图 3：empty-interior risk table visualization（推荐）

文件名建议：

```text
phase8_empty_interior_risk.png
```

---

## 8.4 实验结论目标

这一实验应回答：

- 现有 strongest interior result 是否部分来自 empty-interior 假阳性；
- 后续写论文时哪些图必须加注释；
- 后续 theorem 是否必须限定 interior cardinality 非空的 regime。

---

## 9. 实验 6：reference-window vs ideal-truncation 主导性验证实验

## 9.1 目的

Phase 8 当前最关键的数值判断是：

\[
\|c^{\mathrm{C},k}-\hat c^{(k)}\| \ll \|\hat c^{(k)}-\tilde c^{(k)}\|
\]

即真正大的 gap 主要来自：
\[
\hat c^{(k)}-\tilde c^{(k)}
\]

必须系统验证这一结论在多条件下是否仍成立。

---

## 9.2 必须计算的量

对每个 \((q,T,k)\)：

1. raw-to-reference gap
   \[
   E_{\mathrm{raw\to ref}}(k)
   :=
   \|c^{\mathrm{C},k}-\hat c^{(k)}\|_2
   \]
2. reference-to-ideal gap
   \[
   E_{\mathrm{ref\to ideal}}(k)
   :=
   \|\hat c^{(k)}-\tilde c^{(k)}\|_2
   \]
3. gap ratio
   \[
   R_{\mathrm{gap}}(k)
   :=
   \frac{E_{\mathrm{ref\to ideal}}(k)}{E_{\mathrm{raw\to ref}}(k)+\varepsilon}
   \]

---

## 9.3 必须保存的图片

### 图 1：两类桥梁 gap 对比图（必须）

文件名建议：

```text
phase8_two_bridge_gaps.png
```

### 图 2：gap ratio vs \(k\)（必须）

文件名建议：

```text
phase8_gap_ratio_vs_k.png
```

### 图 3：不同 regime 下的 gap ratio boxplot（推荐）

文件名建议：

```text
phase8_gap_ratio_boxplot.png
```

---

## 9.4 实验结论目标

这一实验应回答：

- 当前下一条应优先补的理论桥是否确实是：
  \[
  \hat c^{(k)} \leftrightarrow \tilde c^{(k)}
  \]
- 还是 boundary attenuation 其实仍不可忽视。

---

## 10. 建议的统一主运行器

建议新增统一脚本：

```python
run_phase8_supplementary_suite(
    M_values=None,
    k_values=None,
    radius_modes=None,
    h_values=None,
    nonuniform_boxes=None,
    n_trials=50,
    s=4,
    save_dir=None,
    seed=42
)
```

必须支持：

1. large-\(M\) sweep
2. more-trials robustness
3. radius sensitivity
4. raw-vs-reference sanity check
5. empty-interior risk analysis
6. two-bridge gap comparison

---

## 11. 统一结果目录要求

所有结果必须保存到：

```text
results/phase8_supplementary/
```

建议目录结构：

```text
results/phase8_supplementary/
├── exp1_large_M/
├── exp2_more_trials/
├── exp3_radius_sensitivity/
├── exp4_raw_vs_reference_sanity/
├── exp5_empty_interior_risk/
├── exp6_two_bridge_gaps/
└── summary/
```

---

## 12. 必须保存的总表和总 summary

### 12.1 总对比表（必须）

文件名建议：

```text
phase8_supplementary_overview.csv
```

至少包含：

- experiment_name
- regime
- M
- k
- radius_mode
- interior_count
- full_matching_l2
- interior_matching_l2
- boundary_matching_l2
- boundary_energy_ratio
- raw_to_ref_l2
- ref_to_ideal_l2
- gap_ratio

### 12.2 总结性 summary（必须）

文件名建议：

```text
phase8_supplementary_summary.md
```

必须自动回答：

1. Phase 8 现有结论是否在更大 \(M\) 与更多 trial 下稳健；
2. interior-first 是否仍成立且不是 empty-interior 假阳性；
3. `raw vs reference-window ≈ 0` 是否可信；
4. 下一条 bridge theorem 最应优先补哪一条；
5. Phase 8 是否已经足以支撑论文中的强 observation，还是还需继续补实验。

---

## 13. 推荐代码结构

```text
phase_8/
├── blom_phase8_supplementary_suite.py
├── exp1_large_M_sweep.py
├── exp2_more_trials.py
├── exp3_radius_sensitivity.py
├── exp4_raw_vs_reference_sanity.py
├── exp5_empty_interior_risk.py
├── exp6_two_bridge_gap_compare.py
├── test_phase8_supplementary_suite.py
├── examples/
│   ├── demo_exp1_large_M.py
│   ├── demo_exp2_more_trials.py
│   ├── demo_exp3_radius.py
│   ├── demo_exp4_raw_ref.py
│   ├── demo_exp5_empty_interior.py
│   ├── demo_exp6_two_bridge.py
│   └── demo_phase8_supplementary_suite.py
└── results/
    └── phase8_supplementary/
```

---

## 14. 验收标准

只有满足以下条件，本需求对应的实现才算完成：

1. 与前面 phases 的符号系统完全统一；
2. 六类补充实验都能独立运行；
3. 能保存所有要求的图片、表格、JSON、summary；
4. 能明确区分“真实 interior 优势”和“empty-interior 假阳性”；
5. 能对 `raw vs reference-window` 做额外 sanity check；
6. 能明确回答下一条 theorem 该优先补哪一条；
7. 结果足以作为 Phase 8 是否进一步补实验的直接决策依据。

---

## 15. 给实现 AI 的最终指令

你要实现的是一套 **Phase 8 补充实验框架**，目的是判断当前 Phase 8 的 strongest observation 是否已经足够稳健。

请优先保证：

1. 结果必须能区分“真实结构现象”与“定义/样本造成的假象”；
2. \(M\)、trial 数、radius mode 三个维度都必须系统 sweep；
3. `raw vs reference-window` 必须单独做 sanity check；
4. empty-interior 风险必须显式检测，不允许被忽略；
5. 所有图片与 summary 必须足够直观，能直接用于论文与组会汇报；
6. 不要把数值现象直接写成 theorem。

本阶段的核心任务是：**判断当前 Phase 8 是否已经足够强，还是还需要更多实验才能安全地支撑下一条理论桥梁的选择。**
