# REQ-Phase8-interior-matching-validation.md

## 1. 文档目的

本需求文档用于指导实现 **Phase 8：interior-first matching 与 boundary-gap 分解** 的统一数值验证代码。

目标不是直接证明完整的 full-domain matching theorem，而是把 Phase 8 的理论对象落到统一可复现的数值验证框架中，重点回答：

1. raw Scheme C 与 idealized truncated BLOM-\(k\) 的差异，是否在 **interior segments** 上明显更小；
2. full matching error 是否可以数值上稳定分解为：
   - interior matching part
   - boundary-layer matching part
3. 在 uniform-time 与 bounded-nonuniform 两种 regime 下，interior-first 现象是否都稳定存在；
4. boundary gap 是否足够显著，值得后续单独补 theorem；
5. 后续理论应优先继续：
   - interior-first matching theorem
   - boundary-layer remainder theorem

   还是需要回头修改：
   - raw Scheme C 定义
   - window 设计
   - boundary handling

本需求文档对应以下三个脚本：

- `blom_interior_matching_check.py`
- `blom_boundary_gap_decomposition.py`
- `blom_uniform_vs_nonuniform_interior.py`

---

## 2. 与前面 Phases 完全统一的符号要求

实现必须与 Phase 0--8 的符号保持一致。

### 2.1 轨迹与段

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

### 2.2 canonical setting

默认 canonical case：

\[
s=4,\qquad d_i=1,\qquad k\in 2\mathbb N
\]

其中 \(k\) 是 Phase 7 / Phase 8 讨论中的 window-size parameter。

### 2.3 三个核心系数对象

- 全局 MINCO 参考系数：
  \[
  c^\star(q,T)\in\mathbb R^{2sM}
  \]
- idealized truncated BLOM-\(k\)：
  \[
  \tilde c^{(k)}(q,T)\in\mathbb R^{2sM}
  \]
- actual raw Scheme C：
  \[
  c^{\mathrm{C},k}(q,T)\in\mathbb R^{2sM}
  \]

### 2.4 block 表示

必须支持 block form：

\[
c^\star =
\begin{bmatrix}
c_1^\star\\ \vdots\\ c_M^\star
\end{bmatrix},
\qquad
\tilde c^{(k)} =
\begin{bmatrix}
\tilde c_1^{(k)}\\ \vdots\\ \tilde c_M^{(k)}
\end{bmatrix},
\qquad
c^{\mathrm{C},k} =
\begin{bmatrix}
c_1^{\mathrm{C},k}\\ \vdots\\ c_M^{\mathrm{C},k}
\end{bmatrix},
\]
其中每个 block 满足
\[
c_i \in \mathbb R^{2s}.
\]

canonical \(s=4\) 下：
\[
c_i\in \mathbb R^8.
\]

---

## 3. 本 Phase 8 数值验证要对应的理论问题

Phase 8 的核心不是 full-domain convergence，而是 **matching gap 的结构拆分**。

### 3.1 interior-first matching 问题

要研究：

\[
\|c_i^{\mathrm{C},k}-\tilde c_i^{(k)}\|
\]

对 interior segment \(i\) 是否比 full-domain 全体 segment 的 matching error 更小、更干净、更接近指数衰减。

### 3.2 full matching error 的 interior / boundary 分解

要研究：

\[
e^{(k)}:=c^{\mathrm{C},k}-\tilde c^{(k)}.
\]

并把它数值上拆成：

\[
e^{(k)} = P_{\mathrm{int}}^{(k)} e^{(k)} + P_{\mathrm{bnd}}^{(k)} e^{(k)},
\]

其中：

- \(P_{\mathrm{int}}^{(k)}\)：interior block projector
- \(P_{\mathrm{bnd}}^{(k)}\)：boundary-layer block projector

然后分别统计：

\[
\|P_{\mathrm{int}}^{(k)} e^{(k)}\|_2,
\qquad
\|P_{\mathrm{bnd}}^{(k)} e^{(k)}\|_2.
\]

### 3.3 uniform-time vs bounded-nonuniform

要研究在两种 regime 下：

- \(T_i\equiv h\)
- \(0<T_{\min}\le T_i\le T_{\max}\)

interior-first 现象是否都稳定。

---

## 4. interior set 与 boundary-layer set 的定义要求

必须支持一个明确的 boundary-layer radius：

\[
r(k)\in \mathbb N,\qquad r(k)\ge \frac{k}{2}.
\]

建议默认定义：

```python
def boundary_radius(k: int, mode: str = "default") -> int:
    # default: r(k)=k
```

默认建议先取：

\[
r(k)=k
\]

因为这比最小窗口半宽 \(k/2\) 更保守，更适合先做实验。

### interior index set

\[
\mathcal I_{\mathrm{int}}(k)
=
\{\, i\in\{1,\dots,M\}: i-r(k)\ge 1,\ i+r(k)\le M \,\}
\]

### boundary-layer index set

\[
\mathcal I_{\mathrm{bnd}}(k)
=
\{1,\dots,M\}\setminus \mathcal I_{\mathrm{int}}(k)
\]

代码中必须把这两个集合显式输出并可复用。

---

## 5. 脚本 1：`blom_interior_matching_check.py`

### 5.1 目的

对 fixed \((q,T)\) 与多组 \(k\)：

- 比较 full matching error
- interior-only matching error
- per-segment matching error
- \(\log\) error vs \(k\) 的拟合

核心问题：

> interior-only matching 是否显著优于 full matching？

### 5.2 必须实现的函数

```python
def make_interior_sets(M, k, radius_mode="default"):
    ...
```

输出：

```python
{
    "interior_idx": list[int],
    "boundary_idx": list[int],
    "r_k": int,
}
```

```python
def compute_matching_error_blocks(c_actual, c_ideal, s=4):
    ...
```

输出每段 block error：
\[
\|c_i^{\mathrm{C},k}-\tilde c_i^{(k)}\|_2
\]

```python
def summarize_interior_matching(c_actual, c_ideal, M, k, s=4, radius_mode="default"):
    ...
```

至少输出：

```python
{
    "k": int,
    "full_l2": float,
    "interior_l2": float,
    "boundary_l2": float,
    "full_linf": float,
    "interior_linf": float,
    "boundary_linf": float,
    "segment_errors": np.ndarray,  # shape (M,)
    "interior_idx": list[int],
    "boundary_idx": list[int],
}
```

```python
def run_interior_matching_check(q, T, k_values, s=4, radius_mode="default",
                                save_dir=None, seed=42):
    ...
```

---

### 5.3 必须保存的图片

#### 图 1：full vs interior matching error 曲线（必须）

文件名建议：

```text
interior_matching_full_vs_interior.png
```

横轴：\(k\)  
纵轴：
- full matching error
- interior-only matching error
- boundary-only matching error

#### 图 2：\(\log\) matching error vs \(k\) 拟合图（必须）

文件名建议：

```text
interior_matching_logfit.png
```

至少拟合：

\[
\log \|P_{\mathrm{int}}^{(k)}(c^{\mathrm{C},k}-\tilde c^{(k)})\|_2
\]
与 \(k\) 的关系。

#### 图 3：per-segment matching heatmap（必须）

文件名建议：

```text
segmentwise_matching_heatmap.png
```

横轴：segment index \(i\)  
纵轴：\(k\)  
颜色：
\[
\|c_i^{\mathrm{C},k}-\tilde c_i^{(k)}\|_2
\]

#### 图 4：boundary vs interior 箱线图（推荐）

文件名建议：

```text
boundary_vs_interior_matching_boxplot.png
```

---

### 5.4 必须保存的表和摘要

- `interior_matching_by_k.csv`
- `interior_matching_fit_summary.csv`
- `summary_interior_matching.json`
- `phase8_interior_matching_summary.md`

#### `phase8_interior_matching_summary.md` 至少应自动回答：

1. interior-only matching 是否明显优于 full；
2. interior slope 是否更负；
3. boundary term 是否主导了 full matching；
4. 该结果是否支持继续补 interior-first matching theorem。

---

## 6. 脚本 2：`blom_boundary_gap_decomposition.py`

### 6.1 目的

把 full matching error 更明确地拆成：

- interior part
- boundary-layer part

并展示 boundary gap 的空间分布和尺度。

目标不是“直接消灭边界误差”，而是把它从一个模糊现象变成清晰的可视对象。

### 6.2 必须实现的函数

```python
def project_interior_boundary_errors(segment_errors, interior_idx, boundary_idx):
    ...
```

输出：

```python
{
    "interior_l2": float,
    "boundary_l2": float,
    "interior_l1": float,
    "boundary_l1": float,
    "interior_energy_ratio": float,
    "boundary_energy_ratio": float,
}
```

其中建议定义：

\[
\text{boundary energy ratio}
=
\frac{\|P_{\mathrm{bnd}}^{(k)}e^{(k)}\|_2^2}{\|e^{(k)}\|_2^2}
\]

```python
def run_boundary_gap_decomposition(q, T, k_values, s=4, radius_mode="default",
                                   save_dir=None, seed=42):
    ...
```

---

### 6.3 必须保存的图片

#### 图 1：boundary energy ratio vs \(k\)（必须）

文件名建议：

```text
boundary_energy_ratio_vs_k.png
```

#### 图 2：segmentwise boundary-gap profile（必须）

文件名建议：

```text
boundary_gap_profile.png
```

对固定几个代表性的 \(k\)，画：
\[
i \mapsto \|c_i^{\mathrm{C},k}-\tilde c_i^{(k)}\|_2
\]

并标记 interior region / boundary region。

#### 图 3：interior / boundary stacked bar chart（必须）

文件名建议：

```text
interior_boundary_stacked_error.png
```

#### 图 4：不同 \(r(k)\) 选择下的敏感性图（推荐）

文件名建议：

```text
boundary_radius_sensitivity.png
```

即比较 \(r(k)=k/2, k, 3k/2\) 等选择下的结论是否稳健。

---

### 6.4 必须保存的表和摘要

- `boundary_gap_decomposition.csv`
- `boundary_radius_sensitivity.csv`
- `summary_boundary_gap.json`
- `phase8_boundary_gap_summary.md`

#### `phase8_boundary_gap_summary.md` 至少应自动回答：

1. boundary-layer term 是否占 full matching 的主要比例；
2. 这种 dominance 是否随 \(k\) 减弱；
3. 结论是否对 boundary radius 的选取敏感；
4. 这是否支持后续补 boundary-layer theorem / remainder theorem。

---

## 7. 脚本 3：`blom_uniform_vs_nonuniform_interior.py`

### 7.1 目的

对应 Phase 8 理论中两层版本：

- uniform-time
- bounded-nonuniform time

检查 interior-first matching 是否在两种 regime 下都稳定。

### 7.2 实验组设计

#### Group U：uniform-time

\[
T_i \equiv h
\]

建议多个 \(h\)：

```python
h_values = [0.5, 1.0, 2.0]
```

#### Group B：bounded-nonuniform

\[
T_i\sim \mathrm{Uniform}(T_{\min},T_{\max})
\]

建议多组：

```python
[(0.9, 1.1), (0.5, 2.0), (0.2, 3.0)]
```

---

### 7.3 必须实现的函数

```python
def sample_uniform_time(M, h):
    ...
```

```python
def sample_bounded_nonuniform_time(M, T_min, T_max, rng):
    ...
```

```python
def run_uniform_vs_nonuniform_interior(
    q_sampler,
    M,
    k_values,
    h_values,
    nonuniform_boxes,
    n_trials=50,
    s=4,
    radius_mode="default",
    save_dir=None,
    seed=42
):
    ...
```

---

### 7.4 必须计算的量

对每个 regime / 每组参数 / 每个 \(k\)，至少统计：

1. full matching error
2. interior-only matching error
3. boundary-only matching error
4. interior log-fit slope
5. interior fit \(R^2\)

---

### 7.5 必须保存的图片

#### 图 1：uniform vs bounded-nonuniform 的 interior matching 对比图（必须）

文件名建议：

```text
uniform_vs_nonuniform_interior_matching.png
```

#### 图 2：uniform vs bounded-nonuniform 的 interior slope 箱线图（必须）

文件名建议：

```text
uniform_vs_nonuniform_slope_boxplot.png
```

#### 图 3：uniform vs bounded-nonuniform 的 boundary energy ratio 图（推荐）

文件名建议：

```text
uniform_vs_nonuniform_boundary_ratio.png
```

#### 图 4：各 time regime 的 representative heatmap（推荐）

文件名建议：

```text
time_regime_matching_heatmaps.png
```

---

### 7.6 必须保存的表和摘要

- `uniform_vs_nonuniform_summary.csv`
- `uniform_vs_nonuniform_fit_summary.csv`
- `summary_uniform_vs_nonuniform.json`
- `phase8_uniform_vs_nonuniform_summary.md`

#### `phase8_uniform_vs_nonuniform_summary.md` 至少应自动回答：

1. interior-first 现象是否在 uniform-time 下更明显；
2. bounded-nonuniform 是否仍支持继续补 theorem；
3. 后续 theorem 是否应先以 uniform-time 为主，再推广。

---

## 8. 三个脚本的统一主运行器要求

建议再提供一个总控脚本：

```python
run_phase8_validation_suite(
    q=None,
    T=None,
    M=20,
    k_values=None,
    h_values=None,
    nonuniform_boxes=None,
    n_trials=50,
    s=4,
    radius_mode="default",
    save_dir=None,
    seed=42
)
```

它至少要完成：

1. single-case interior matching check
2. single-case boundary-gap decomposition
3. regime split random trials
4. 自动汇总总结

---

## 9. 必须保存的总表和总 summary

### 9.1 总对比表（必须）

文件名建议：

```text
phase8_overview.csv
```

至少包含：

- experiment_name
- regime
- k
- full_matching_l2
- interior_matching_l2
- boundary_matching_l2
- interior_slope
- interior_r2
- boundary_energy_ratio

### 9.2 总 summary（必须）

文件名建议：

```text
phase8_interpretation_summary.md
```

必须自动生成结论，至少包括：

1. interior-first matching 是否显著成立；
2. boundary gap 是否必须单独理论处理；
3. uniform-time 与 bounded-nonuniform 是否都支持继续理论推进；
4. 下一步理论最值得补哪一个：
   - boundary attenuation proposition
   - reference-window / ideal-truncation consistency proposition
   - boundary-layer remainder theorem

---

## 10. 推荐代码结构

```text
phase_8/
├── blom_interior_matching_check.py
├── blom_boundary_gap_decomposition.py
├── blom_uniform_vs_nonuniform_interior.py
├── blom_phase8_validation_suite.py
├── test_blom_interior_matching_check.py
├── test_blom_boundary_gap_decomposition.py
├── test_blom_uniform_vs_nonuniform_interior.py
├── examples/
│   ├── demo_interior_matching.py
│   ├── demo_boundary_gap.py
│   ├── demo_uniform_vs_nonuniform.py
│   └── demo_phase8_suite.py
└── results/
    └── phase8_validation/
```

---

## 11. 验收标准

只有满足以下条件，本需求对应的实现才算完成：

1. 与 Phase 0--8 的符号系统完全统一；
2. 三个脚本都能独立运行；
3. 能明确构造 interior / boundary index sets；
4. 能分别统计 full / interior / boundary matching error；
5. 能保存必须的图、表、JSON、summary；
6. 能对 uniform-time 与 bounded-nonuniform 做分层比较；
7. 结果足以直接支撑 Phase 8 论文中的结论与后续 theorem 选择。

---

## 12. 给实现 AI 的最终指令

你要实现的是一套 **Phase 8 interior-first matching 数值验证框架**，不是单个示意脚本。

请优先保证：

1. interior set 与 boundary set 的定义清楚且可重复；
2. raw Scheme C 与 idealized truncation 的比较严格一致；
3. full / interior / boundary 三类误差必须分开报告；
4. 图片必须足够直观，适合直接进入论文或组会汇报；
5. 解释性 summary 必须自动生成；
6. 不要把数值现象直接写成 theorem；
7. 结果必须能告诉我们：下一步该优先补哪个 bridge proposition。

本阶段的核心任务是：**把 matching 问题从“模糊地还差一点”变成“明确知道 interior 和 boundary 各差多少、下一步该补哪条理论桥梁”。**
