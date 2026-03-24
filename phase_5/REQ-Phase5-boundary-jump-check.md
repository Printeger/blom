# REQ-Phase5-boundary-jump-check.md

## 1. 文档目的

本需求文档用于指导实现 `blom_boundary_jump_check.py`，作为 BLOM 理论 **Phase 5：全局装配机制** 的统一验证脚本。

该脚本的目标不是重新求解 MINCO 或 BLOM-Strict/BLOM-Analytic，而是：

1. 对 **方案 A / 方案 B / 方案 C** 的全局装配结果进行统一验证；
2. 在与 Phase 0--5 完全一致的符号体系下，计算并统计各 interior knot 处的导数 jump：
   \[
   \llbracket p^{(\ell)} \rrbracket(t_i)
   \]
3. 用数值证据回答以下问题：
   - 哪种装配方案可以严格达到 \(C^{s-1}\)？
   - 哪种方案只保证 \(C^0\)？
   - 哪种方案在更高阶导数上 jump 更小？
   - 在随机 \(q,T\) 下，各方案的稳定性、局部性和光滑性表现如何？
4. 自动保存必要的结果图和统计表，用于论文展示与方案选择。

本脚本应与 Phase 1、Phase 2、Phase 3、Phase 4 的实现保持接口兼容，尤其要兼容：

- `minco_scalar_baseline.py`
- `blom_strict_local_qp.py`
- `blom_k2_s4_numeric.py`

---

## 2. 与前面 Phases 的符号统一要求

实现必须保持以下符号与前文一致：

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
- 控制阶次：
  \[
  s
  \]
  canonical case 默认
  \[
  s=4
  \]
- canonical bandwidth：
  \[
  k=2
  \]
- 内点 knot：
  \[
  t_i=\sum_{r=1}^{i}T_r,\qquad i=1,\dots,M-1
  \]
- 导数 jump 定义：
  \[
  \llbracket p^{(\ell)}\rrbracket(t_i)
  :=
  p_i^{(\ell)}(T_i)-p_{i+1}^{(\ell)}(0)
  \]
- lower-order knot-state：
  \[
  \eta_i=
  \begin{bmatrix}
  p^{(1)}(t_i)\\
  \vdots\\
  p^{(s-1)}(t_i)
  \end{bmatrix}
  \]
- 在 canonical \(s=4\) 下：
  \[
  \eta_i=\chi_i=
  \begin{bmatrix}
  v_i\\ a_i\\ j_i
  \end{bmatrix}
  \in\mathbb R^3
  \]

不允许在代码中改名为与前文冲突的另一套核心符号体系。

---

## 3. 本脚本要验证的三种装配方案

## 3.1 方案 A：Shared junction states

### 数学对象

引入全局共享 knot-state
\[
\eta=(\eta_1^\top,\dots,\eta_{M-1}^\top)^\top
\]
并通过全局二次目标
\[
\mathcal J_A(\eta;q,T)
\]
求得唯一装配状态 \(\eta^\star\)，然后用统一 Hermite 重构每一段：
\[
c_i^A=\mathcal H_i(q_{i-1},q_i,T_i,\eta_{i-1}^\star,\eta_i^\star).
\]

### 理论预期

- 必须严格满足
  \[
  \llbracket p^{(\ell)}\rrbracket(t_i)=0,\qquad \ell=0,1,\dots,s-1
  \]
- 不要求自动满足
  \[
  \ell=s,\dots,2s-2
  \]

### 要输出的对象

至少输出：

```python
{
    "scheme": "A",
    "eta_shared": np.ndarray,      # shape ((M-1), s-1)
    "coeffs": np.ndarray,          # shape (M, 2*s)
    "jumps": dict,                 # order -> ndarray shape (M-1,)
    "stats": dict,
}
```

---

## 3.2 方案 B：Overlapping consensus

### 数学对象

每个窗口独立求解 BLOM-Strict / BLOM-Analytic，得到对每个 interior knot 的局部预测
\[
\widehat\eta_i^{(r)}.
\]
然后做 weighted consensus：
\[
\eta_i^\star
=
\frac{\sum_{r\in\mathcal I(i)}\omega_i^{(r)}\widehat\eta_i^{(r)}}
{\sum_{r\in\mathcal I(i)}\omega_i^{(r)}}.
\]
再统一 Hermite 重构：
\[
c_i^B=\mathcal H_i(q_{i-1},q_i,T_i,\eta_{i-1}^\star,\eta_i^\star).
\]

### 理论预期

- 必须严格满足
  \[
  \llbracket p^{(\ell)}\rrbracket(t_i)=0,\qquad \ell=0,1,\dots,s-1
  \]
- 不要求自动满足更高阶连续性；
- 与方案 A 相比，更接近 sliding-window 思想。

### 要输出的对象

```python
{
    "scheme": "B",
    "eta_predictions": list,       # 每个 knot 的窗口预测列表
    "eta_consensus": np.ndarray,   # shape ((M-1), s-1)
    "weights": dict,
    "coeffs": np.ndarray,          # shape (M, 2*s)
    "jumps": dict,
    "stats": dict,
}
```

---

## 3.3 方案 C：Central-segment extraction

### 数学对象

直接采用每个窗口的中心段：
\[
c_i^C:=c_i^{\mathrm{BLOM},k}.
\]
不做共享变量、不做 consensus、不做额外装配。

### 理论预期

- 必须严格满足
  \[
  \llbracket p\rrbracket(t_i)=0
  \]
  即 \(C^0\)
- 对于
  \[
  \ell=1,\dots,s-1
  \]
  一般**不保证** jump 为零
- 对于高阶导数 jump，通常也不为零

### 要输出的对象

```python
{
    "scheme": "C",
    "coeffs": np.ndarray,          # shape (M, 2*s)
    "eta_left": np.ndarray,        # shape ((M-1), s-1)
    "eta_right": np.ndarray,       # shape ((M-1), s-1)
    "jumps": dict,
    "stats": dict,
}
```

---

## 4. 脚本核心功能要求

## 4.1 统一 jump 计算接口

必须实现统一函数：

```python
compute_jumps(coeffs, T, s, max_order=None) -> dict
```

### 功能

对每个 interior knot \(t_i\) 和每个导数阶 \(\ell\)，计算
\[
\llbracket p^{(\ell)}\rrbracket(t_i)
=
p_i^{(\ell)}(T_i)-p_{i+1}^{(\ell)}(0)
\]

### 输入

- `coeffs`: shape `(M, 2*s)`
- `T`: shape `(M,)`
- `s`: int
- `max_order`: 默认 `2*s-2`

### 输出

```python
{
    0: np.ndarray,   # shape (M-1,)
    1: np.ndarray,
    ...
    2*s-2: np.ndarray
}
```

---

## 4.2 统一统计接口

必须实现：

```python
summarize_jumps(jumps: dict) -> dict
```

### 对每一阶 \(\ell\) 必须统计

1. 最大绝对 jump
   \[
   \max_i |\llbracket p^{(\ell)}\rrbracket(t_i)|
   \]
2. 平均绝对 jump
   \[
   \frac{1}{M-1}\sum_i |\llbracket p^{(\ell)}\rrbracket(t_i)|
   \]
3. 中位数绝对 jump
4. RMS jump
   \[
   \sqrt{\frac{1}{M-1}\sum_i \llbracket p^{(\ell)}\rrbracket(t_i)^2}
   \]
5. 95% 分位数绝对 jump
6. 是否接近机器精度为零（容忍度默认 \(10^{-10}\)）

### 输出格式

```python
{
    order: {
        "max_abs": ...,
        "mean_abs": ...,
        "median_abs": ...,
        "rms": ...,
        "q95_abs": ...,
        "is_zero_tol": bool,
    }
}
```

---

## 4.3 三种方案的统一运行器

必须实现高层接口：

```python
run_boundary_jump_check(
    q, T, s=4, k=2,
    scheme="A",
    config=None,
    save_dir=None,
    seed=42
) -> dict
```

### 行为要求

- 根据 `scheme` 调用 A / B / C 对应装配逻辑
- 输出统一结构：
  - 系数
  - jump
  - 统计量
  - 结果图路径
  - 元信息
- 若 `save_dir` 不为空，必须保存图与表

---

## 5. 依赖前面 phases 的接口要求

本脚本不能重复发明前面 phase 已经实现的对象，而必须尽量复用。

## 5.1 来自 Phase 1

应复用：

- 单段多项式评估接口
- 导数评估接口
- Hermite / monomial 基函数接口
- 轨迹采样接口（如果已有）

## 5.2 来自 Phase 3

应复用或兼容：

- BLOM-Strict 局部问题求解器
- 局部窗口定义
- 自然边界设置
- 局部唯一解输出格式

## 5.3 来自 Phase 4

应复用或兼容：

- `blom_k2_s4_numeric.py`
- 局部状态
  \[
  x_i^{\mathrm{loc}}
  \]
- Hermite 重构
  \[
  c_i=\mathcal H_i(\cdot)
  \]
- central segment coefficient extraction

---

## 6. 方案 A / B / C 的具体验证要求

## 6.1 方案 A：必须验证的命题

### 理论命题

应验证：

\[
\llbracket p^{(\ell)}\rrbracket(t_i)=0,\qquad \ell=0,1,\dots,s-1
\]

### 数值检查

对所有 interior knot 和所有
\[
\ell=0,\dots,s-1
\]
统计：

- `max_abs`
- `mean_abs`
- `rms`

并判断是否小于给定阈值，例如：

```python
tol_exact = 1e-10
```

### 额外检查

- Hessian / assembled system 是否对称
- 若实现了求解矩阵，检查条件数
- 若求解失败，记录失败案例

---

## 6.2 方案 B：必须验证的命题

### 理论命题

也应验证：

\[
\llbracket p^{(\ell)}\rrbracket(t_i)=0,\qquad \ell=0,\dots,s-1
\]

### 数值检查

同方案 A。

### 额外比较项

方案 B 还必须额外输出：

- 每个 knot 的局部预测个数
- 每个 knot 的 consensus 前分散度，例如：
  \[
  \max_{r\in\mathcal I(i)} \|\widehat\eta_i^{(r)}-\bar\eta_i\|_2
  \]
- consensus 前后 jump 的改善程度

这部分是方案 B 的核心价值展示。

---

## 6.3 方案 C：必须验证的命题

### 理论命题

必须验证：

\[
\llbracket p\rrbracket(t_i)=0
\]

但对
\[
\ell=1,\dots,s-1
\]
不能默认设为 0，而是要当作诊断量。

### 数值检查

至少计算：

- \( \ell=0 \) 时 jump 是否机器精度为零
- \( \ell=1,\dots,s-1 \) 时 jump 分布
- \( \ell=s,\dots,2s-2 \) 时高阶 jump 分布

### 额外输出

必须显式输出：

\[
\eta_i^-,
\qquad
\eta_i^+
\]
以及它们的差：
\[
\eta_i^- - \eta_i^+
\]
因为这正是方案 C 在 lower-order continuity 上失败或接近成功的直接证据。

---

## 7. 必须保存的结果图片

所有图片必须保存到：

```text
results/phase5_boundary_jump_check/
```

并按方案分目录：

```text
results/phase5_boundary_jump_check/
├── scheme_A/
├── scheme_B/
└── scheme_C/
```

## 7.1 图 1：各阶导数 jump 的热图（必须）

文件名建议：

```text
jump_heatmap_scheme_A.png
jump_heatmap_scheme_B.png
jump_heatmap_scheme_C.png
```

### 内容

横轴：interior knot index \(i\)  
纵轴：导数阶 \(\ell\)  
颜色：\(\log_{10}(|\llbracket p^{(\ell)}\rrbracket|+\varepsilon)\)

### 目的

直观看：
- 哪些阶 jump 最严重
- 哪些方案能把 \(\ell\le s-1\) 压到数值零
- 哪些 knot 最不稳定

---

## 7.2 图 2：各阶导数 jump 最大值柱状图（必须）

文件名建议：

```text
jump_maxbar_scheme_A.png
jump_maxbar_scheme_B.png
jump_maxbar_scheme_C.png
```

### 内容

横轴：导数阶 \(\ell\)  
纵轴：
\[
\max_i |\llbracket p^{(\ell)}\rrbracket(t_i)|
\]

### 目的

用于论文中做方案对比，最直观。

---

## 7.3 图 3：lower-order jump 对比图（必须）

文件名建议：

```text
jump_lower_orders_compare.png
```

### 内容

在一张图里同时对比方案 A / B / C 的
\[
\ell=0,\dots,s-1
\]
的最大 jump 或 RMS jump。

### 目的

直接展示：

- A、B 是否达到理论上的 \(C^{s-1}\)
- C 是否仅在 \(C^0\) 成立

---

## 7.4 图 4：方案 C 的左右 knot-state 差异图（必须）

文件名建议：

```text
scheme_C_eta_mismatch.png
```

### 内容

横轴：interior knot index \(i\)  
纵轴：\(\|\eta_i^- - \eta_i^+\|_2\)

### 目的

直接可视化说明为什么方案 C 一般无法自动达到 \(C^{s-1}\)。

---

## 7.5 图 5：方案 B 的 consensus 改善图（推荐）

文件名建议：

```text
scheme_B_consensus_improvement.png
```

### 内容

对每个 knot 比较：

- consensus 前局部预测分散度
- consensus 后 jump / mismatch

### 目的

展示方案 B 相对于纯 raw extraction 的好处。

---

## 7.6 图 6：随机样本统计箱线图（推荐）

文件名建议：

```text
jump_boxplot_random_trials.png
```

### 内容

在多组随机 \(q,T\) 下，对比方案 A / B / C 的：

- lower-order jump
- higher-order jump

### 目的

展示方案稳健性，不只看单个案例。

---

## 8. 必须保存的表格/数据文件

## 8.1 每次运行的 JSON 摘要

文件名建议：

```text
summary_scheme_A.json
summary_scheme_B.json
summary_scheme_C.json
```

必须保存：

- 输入 `q,T,s,k`
- scheme 名称
- 每阶 jump 统计
- 求解是否成功
- 耗时
- 若失败，失败原因

## 8.2 CSV 统计表（必须）

文件名建议：

```text
jump_stats_scheme_A.csv
jump_stats_scheme_B.csv
jump_stats_scheme_C.csv
```

每行一个导数阶 \(\ell\)，列包含：

- `order`
- `max_abs`
- `mean_abs`
- `median_abs`
- `rms`
- `q95_abs`
- `is_zero_tol`

## 8.3 总对比表（必须）

文件名建议：

```text
scheme_comparison_summary.csv
```

至少包含：

- scheme
- max jump over \(\ell=0,\dots,s-1\)
- rms jump over \(\ell=0,\dots,s-1\)
- max jump over \(\ell=s,\dots,2s-2\)
- 是否严格达到 \(C^{s-1}\)
- 是否需要全局变量
- 是否需要 consensus
- 局部性评分（人工标签即可）
- 运行时间

---

## 9. 随机测试要求

必须支持多组随机 \(q,T\) 自动测试。

建议函数：

```python
run_random_trials(
    n_trials=100,
    M=20,
    s=4,
    k=2,
    q_sampler=None,
    T_sampler=None,
    schemes=("A", "B", "C"),
    save_dir=None,
    seed=42
)
```

### 默认采样建议

- `q`：高斯随机游走或均匀随机 waypoint
- `T`：正值随机，例如
  \[
  T_i \sim \mathrm{Uniform}(0.5, 2.0)
  \]

### 必须统计

- 各方案成功率
- 各方案 lower-order jump 最大值分布
- 各方案 higher-order jump 最大值分布
- 各方案平均运行时间

---

## 10. 方案 A / B / C 的对比结论输出要求

脚本最后必须自动生成一个简短文字摘要（txt 或 markdown），内容包括：

### 对方案 A
- 是否在 lower-order jump 上达到机器精度
- 是否需要额外全局共享变量
- 相比 B/C 的优缺点

### 对方案 B
- 是否在 lower-order jump 上达到机器精度
- consensus 是否明显降低 mismatch
- 是否比 A 更保留 sliding-window 风格

### 对方案 C
- 是否仅在 \(C^0\) 成立
- lower-order jump 是否显著非零
- 是否仍然具有最干净的局部映射解释

文件名建议：

```text
phase5_interpretation_summary.md
```

---

## 11. 推荐代码结构

```text
phase_5/
├── blom_boundary_jump_check.py
├── test_blom_boundary_jump_check.py
├── examples/
│   ├── demo_scheme_A.py
│   ├── demo_scheme_B.py
│   ├── demo_scheme_C.py
│   └── demo_compare_all.py
└── results/
    └── phase5_boundary_jump_check/
```

`blom_boundary_jump_check.py` 至少应包含：

```python
def compute_jumps(coeffs, T, s, max_order=None): ...
def summarize_jumps(jumps): ...

def assemble_scheme_A(q, T, s=4, k=2, config=None): ...
def assemble_scheme_B(q, T, s=4, k=2, config=None): ...
def assemble_scheme_C(q, T, s=4, k=2, config=None): ...

def run_boundary_jump_check(q, T, s=4, k=2, scheme="A", config=None,
                            save_dir=None, seed=42): ...

def run_compare_all_schemes(q, T, s=4, k=2, config=None,
                            save_dir=None, seed=42): ...

def run_random_trials(n_trials=100, M=20, s=4, k=2,
                      q_sampler=None, T_sampler=None,
                      schemes=("A", "B", "C"),
                      save_dir=None, seed=42): ...
```

---

## 12. 验收标准

只有当以下条件同时满足时，本脚本才算完成：

1. 与 Phase 0--5 的符号和接口一致；
2. 能独立运行方案 A / B / C；
3. 能统一计算
   \[
   \llbracket p^{(\ell)}\rrbracket(t_i)
   \]
   到 \(2s-2\) 阶；
4. 能保存所有必须图片；
5. 能保存 CSV / JSON 统计结果；
6. 能自动给出三种方案的数值对比；
7. 能在 canonical case \(s=4,k=2\) 下运行稳定；
8. 能在多组随机 \(q,T\) 下批量测试；
9. 结果足以直接用于论文图表和方案选择讨论。

---

## 13. 给实现 AI 的最终指令

你要实现的是一个**统一的 Phase 5 全局装配验证脚本**，不是单一方案 demo。

请优先保证：

1. jump 定义和计算严格正确；
2. 方案 A / B / C 都能运行；
3. 统计量和图像输出完整；
4. lower-order / higher-order jump 分开报告；
5. 与前面 phase 的 BLOM 局部求解和 Hermite 重构接口兼容；
6. 结果目录清晰、命名统一、可复现；
7. 所有关键图都自动保存。

不要为了“写得快”而把方案 A / B / C 混在一起，或只做单次案例不做统计。
本脚本的核心任务是：**给出装配机制的数学验证证据，并形成可用于论文展示的统一结果。**
