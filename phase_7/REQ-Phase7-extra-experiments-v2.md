# REQ-Phase7-extra-experiments.md

## 1. 文档目的

本需求文档用于指导实现 **Phase 7 之后必须补做的四类关键实验**，作为 BLOM 理论分叉决策的统一数值验证补充。

该文档的目标是把目前“还需要更多证据”的四类实验整理为可直接执行的实验任务，供实现 AI 参考并用于 vibe coding。

这四类实验分别是：

1. **更大 \(M\) 的 sweep 实验**
2. **uniform-time vs bounded-nonuniform 分层实验**
3. **interior-only error vs full error 分离实验**
4. **Scheme C 的轻装配修正版实验**

这些实验的共同目标是回答：

- actual BLOM-\(k\) 的指数衰减趋势是否随着问题规模增加仍然稳定；
- 当前数值现象是否依赖某类时间分配；
- 误差是否主要由 boundary effects 驱动；
- Scheme C 的 gap 是否可以通过轻量一致化显著缩小；
- 后续理论应优先继续补：
  - `matching theorem`
  - `feasibility / cost theorem`
  
  还是应当回头修改：
  - BLOM 定义
  - 装配方案
  - 比较对象
  - 主 claim

该文档默认与以下脚本兼容：

- `minco_scalar_baseline.py`
- `blom_k2_s4_numeric.py`
- `blom_boundary_jump_check.py`
- `blom_fd_jacobian_check.py`
- `blom_convergence_vs_k.py`

---

## 2. 与前面 Phases 完全统一的符号要求

实现必须保持以下符号体系与 Phase 0--7 一致：

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
- 时间向量：
  \[
  T=(T_1,\dots,T_M)^\top,\qquad T_i>0
  \]
- 总时长：
  \[
  \mathcal T=\sum_{i=1}^{M}T_i
  \]
- canonical setting：
  \[
  s=4,\qquad k=2
  \]
- 全局 MINCO 参考系数：
  \[
  c^\star(q,T)
  \]
- 理想截断 BLOM-\(k\)：
  \[
  \tilde c^{(k)}(q,T)
  \]
- actual BLOM-\(k\)：
  \[
  c^{\mathrm{BLOM},k}(q,T)
  \]
- 控制代价：
  \[
  J(c,T)
  \]
- 外层参数：
  \[
  \theta=(q_1,\dots,q_{M-1},T_1,\dots,T_M)^\top
  \]

误差默认统一为：

- 系数误差：
  \[
  E_{\mathrm{coef}}(k)=\|c^{(k)}-c^\star\|_2
  \]
- matching 误差：
  \[
  E_{\mathrm{match}}(k)=\|c^{\mathrm{BLOM},k}-\tilde c^{(k)}\|_2
  \]
- 相对代价 gap：
  \[
  \Delta J_{\mathrm{rel}}(k)=\frac{J(c^{(k)},T)}{J(c^\star,T)}-1
  \]

---

## 3. 统一结果目录要求

所有结果必须保存到：

```text
results/phase7_extra_experiments/
```

目录建议为：

```text
results/phase7_extra_experiments/
├── exp1_large_M_sweep/
├── exp2_time_regime_split/
├── exp3_interior_vs_full/
├── exp4_schemeC_light_assembly/
└── compare_summary/
```

每类实验都必须同时保存：

- 图片（PNG）
- 原始统计表（CSV）
- 结构化摘要（JSON）
- 自动解释性结论（Markdown）

---

## 4. 实验 1：更大 \(M\) 的 sweep 实验

### 4.1 实验目的

检验以下问题：

1. actual BLOM-\(k\) 的误差随 \(k\) 下降的趋势在更大 \(M\) 下是否仍稳定；
2. ideal truncation 的“漂亮收敛”是否仍然成立；
3. Scheme C 的 slope 是否随 \(M\) 增大仍保持稳定负值；
4. boundary effects 是否随着 \(M\) 增大被稀释，还是反而更严重。

### 4.2 自变量设计

必须 sweep 多个问题规模，例如：

```python
M_values = [10, 20, 40, 80, 120]
```

每个 \(M\) 下都要生成对应的：

- `q`
- `T`
- `k_values`

建议：

```python
k_values = [2, 4, 6, ..., min(M, k_max)]
k_max = min(M, 20)
```

### 4.3 必须计算的量

对每个 \(M\)、每个 \(k\)、每个方案 A/B/C，至少计算：

1. 全局系数误差
2. matching 误差
3. 相对代价 gap
4. \(\log\)-线性拟合斜率
5. \(R^2\)

### 4.4 必须保存的图片

- `coef_error_vs_k_by_M.png`
- `logfit_vs_k_by_M.png`
- `slope_vs_M.png`
- `matching_vs_M.png`
- `cost_gap_vs_M.png`

### 4.5 实验结论目标

这一实验应帮助回答：

- 若 \(M\) 增大时 slope 仍稳定为负，说明继续补 spectral / matching theorem 是有意义的；
- 若 \(M\) 一大，趋势就崩，说明理论主线可能不稳，应重新审查 BLOM 定义或装配方案。

---

## 5. 实验 2：uniform-time vs bounded-nonuniform 分层实验

### 5.1 实验目的

检验当前数值现象是否依赖时间分配类型。

### 5.2 实验组设计

#### Group U：uniform-time

设
\[
T_i\equiv h
\]
例如：

```python
h_values = [0.5, 1.0, 2.0]
```

#### Group B：bounded-nonuniform

例如：
\[
T_i\sim \mathrm{Uniform}(T_{\min}, T_{\max})
\]

建议多组：

```python
[(0.9, 1.1), (0.5, 2.0), (0.2, 3.0)]
```

### 5.3 必须计算的量

对 uniform 和 bounded-nonuniform 两类分别计算：

1. 误差曲线
2. matching 误差曲线
3. slope 分布
4. \(R^2\) 分布
5. relative cost gap

### 5.4 必须保存的图片

- `uniform_vs_bounded_coef_error.png`
- `uniform_vs_bounded_slope_boxplot.png`
- `time_regime_matching_compare.png`
- `time_regime_cost_gap_compare.png`

### 5.5 实验结论目标

这一实验应帮助回答：

- uniform-time 下是否明显更支持当前理论；
- bounded-nonuniform 是否也稳定；
- 是否需要把后续 theorem 先收缩到 uniform regime。

---

## 6. 实验 3：interior-only error vs full error 分离实验

### 6.1 实验目的

检验当前误差是否主要由边界段驱动。

### 6.2 指标定义

#### full error
\[
E_{\mathrm{full}}(k)=\|c^{(k)}-c^\star\|_2
\]

#### interior-only error
去掉前后 \(b\) 个边界段，例如默认 \(b=2\)：
\[
E_{\mathrm{int}}(k)
=
\left\|
\{c_i^{(k)}-c_i^\star\}_{i=b+1}^{M-b}
\right\|_2
\]

### 6.3 必须计算的量

对 ideal / A / B / C 各自计算：

1. full coefficient error
2. interior-only coefficient error
3. full matching error
4. interior-only matching error
5. full vs interior slope

### 6.4 必须保存的图片

- `full_vs_interior_coef_error.png`
- `full_vs_interior_logfit.png`
- `segmentwise_error_heatmap.png`
- `boundary_vs_interior_boxplot.png`

### 6.5 实验结论目标

这一实验应帮助回答：

- boundary effect 是否主导误差；
- 后续 theorem 是否更应先写成 interior theorem；
- 是否需要优先修边界，而不是重做整个 BLOM 核心。

---

## 7. 实验 4：Scheme C 的轻装配修正版实验

### 7.1 实验目的

检验 Scheme C 当前与 ideal truncation / MINCO 的 gap，是否主要来自“完全不做轻量一致化”。

### 7.2 允许的轻装配修正类型

至少实现一种，建议两种：

#### 方案 C1：endpoint jet averaging
\[
\eta_i^{\mathrm{corr}}=\frac{\eta_i^-+\eta_i^+}{2}
\]

#### 方案 C2：local weighted correction
\[
\eta_i^{\mathrm{corr}}
=
\frac{w_i^- \eta_i^- + w_i^+ \eta_i^+}{w_i^-+w_i^+}
\]

#### 方案 C3（可选）：small regularized reconciliation
\[
\min_{\eta}
\sum_i \|\eta_i-\eta_i^-\|^2+\|\eta_i-\eta_i^+\|^2+\lambda\|\nabla \eta_i\|^2
\]

### 7.3 必须比较的对象

至少比较：

- raw Scheme C
- corrected Scheme C1
- corrected Scheme C2
- ideal truncation
- MINCO reference

### 7.4 必须计算的量

1. 系数误差
2. matching 误差
3. relative cost gap
4. lower-order jump
5. Jacobian locality 宽度（推荐）

### 7.5 必须保存的图片

- `schemeC_raw_vs_corrected_error.png`
- `schemeC_raw_vs_corrected_matching.png`
- `schemeC_raw_vs_corrected_jumps.png`
- `schemeC_raw_vs_corrected_locality_tradeoff.png`

### 7.6 实验结论目标

这一实验应帮助回答：

- very light correction 是否足以显著改善 matching / continuity / cost；
- 是否值得把“修正版 Scheme C”提升为正式理论对象；
- correction 是否破坏了 raw C 最重要的局部性优势。

---

## 8. 统一主运行器要求

建议总控脚本：

```python
run_phase7_extra_experiments(
    q=None,
    T=None,
    M_values=None,
    n_trials=100,
    s=4,
    k_values=None,
    schemes=("A", "B", "C"),
    save_dir=None,
    seed=42
)
```

必须支持：

- 单案例实验
- 多 \(M\) sweep
- random trials
- 各实验可单独运行
- 自动保存所有结果

---

## 9. 必须保存的表和摘要文件

### 9.1 每个实验的 CSV

建议：

```text
exp1_large_M_sweep/summary_large_M.csv
exp2_time_regime_split/summary_time_regime.csv
exp3_interior_vs_full/summary_interior_vs_full.csv
exp4_schemeC_light_assembly/summary_schemeC_correction.csv
```

### 9.2 总对比表（必须）

文件名建议：

```text
compare_summary/final_experiment_overview.csv
```

至少包含：

- experiment_name
- scheme
- best_slope
- best_r2
- final_coef_error
- final_match_error
- final_cost_rel
- supports_matching_theorem (bool-like)
- supports_feasibility_theorem (bool-like)
- suggests_model_revision (bool-like)

### 9.3 总解释性 summary（必须）

文件名建议：

```text
compare_summary/phase7_extra_interpretation_summary.md
```

必须自动生成面向研究决策的总结，至少包含：

1. 哪些现象在更大 \(M\) 下稳定；
2. uniform-time 与 bounded-nonuniform 是否有本质差异；
3. boundary effect 是否是主要误差来源；
4. Scheme C 的轻装配修正是否值得纳入正式定义；
5. 后续理论应该往哪个方向推进。

---

## 10. 推荐代码结构

```text
phase_7/
├── blom_phase7_extra_experiments.py
├── test_blom_phase7_extra_experiments.py
├── examples/
│   ├── demo_exp1_large_M.py
│   ├── demo_exp2_time_regime.py
│   ├── demo_exp3_interior_vs_full.py
│   ├── demo_exp4_schemeC_correction.py
│   └── demo_all_phase7_extra.py
└── results/
    └── phase7_extra_experiments/
```

核心函数建议至少包括：

```python
def run_large_M_sweep(...): ...
def run_time_regime_split(...): ...
def run_interior_vs_full(...): ...
def run_schemeC_light_assembly(...): ...
def run_phase7_extra_experiments(...): ...
```

---

## 11. 验收标准

只有满足以下条件，本需求对应的实现才算完成：

1. 与 Phase 0--7 的符号体系完全统一；
2. 四类实验都能独立运行；
3. 每类实验都能保存必须的图、表、JSON、summary；
4. 能比较 ideal truncation / actual BLOM-\(k\) / A/B/C；
5. 能支持多组随机测试；
6. 结果足以直接支撑下一步理论分叉决策；
7. 结果能直观展示并可用于论文图表。

---

## 12. 给实现 AI 的最终指令

你要实现的是一套 **Phase 7 之后的补充实验框架**，不是单个小 demo。

请优先保证：

1. 四类实验逻辑清晰分离；
2. 全部结果自动保存；
3. 图表足够直观、适合论文展示；
4. 结论性 summary 自动生成；
5. 不要偷换 ideal truncation 与 actual BLOM-\(k\)；
6. 不要只做单案例，必须支持随机试验和规模 sweep。

这套实验的核心任务是：**为后续是继续补 theorem，还是回头改 BLOM 本体，提供足够清楚的数值证据。**
