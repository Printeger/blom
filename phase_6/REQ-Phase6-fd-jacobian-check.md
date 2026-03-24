# REQ-Phase6-fd-jacobian-check.md

## 1. 文档目的

本需求文档用于指导实现 `blom_fd_jacobian_check.py`，作为 BLOM 理论 **Phase 6：基础性质验证** 的统一 Jacobian 验证脚本。

该脚本的目标不是重新推导 BLOM 理论，而是把 Phase 6 已经给出的基础性质，特别是：

- Proposition 1：局部支撑（local support）
- Proposition 4：弱版时间感知（局部系数关于 \(T\) 连续可微，局部 Jacobian 在紧致可容许集合上有界）

转化为**统一可复现的数值验证程序**。

该脚本必须完成以下任务：

1. 检查 raw BLOM-Analytic 局部系数映射
   \[
   c_i^{\mathrm{BLOM},2}
   \]
   对 waypoint 和 duration 的 Jacobian 稀疏模式是否与理论一致；
2. 用 finite difference 对解析 Jacobian 进行逐项核对；
3. 在可行时，用 autodiff 或符号微分结果对 finite difference 再做交叉核对；
4. 对方案 A / B / C 的装配后对象分别验证“哪些 Jacobian 还应保持局部、哪些不应再保持严格局部”；
5. 自动保存必要的图、表和摘要文件，用于论文展示与性质证明。

该脚本应与前面 phases 的实现兼容，尤其要兼容：

- `minco_scalar_baseline.py`
- `blom_strict_local_qp.py`
- `blom_k2_s4_numeric.py`
- `blom_boundary_jump_check.py`

---

## 2. 与前面 Phases 完全统一的符号要求

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
- canonical control order：
  \[
  s=4
  \]
- canonical bandwidth：
  \[
  k=2
  \]
- 外层稀疏参数：
  \[
  \theta=(q_1,\dots,q_{M-1},T_1,\dots,T_M)^\top
  \]
- raw BLOM-Analytic 局部映射：
  \[
  c_i^{\mathrm{BLOM},2}
  =
  \mathcal M_{2,i}(q_{i-2},q_{i-1},q_i,q_{i+1},T_{i-1},T_i,T_{i+1})
  \]
- 内点 knot-state：
  \[
  \eta_i=
  \begin{bmatrix}
  p^{(1)}(t_i)\\
  \vdots\\
  p^{(s-1)}(t_i)
  \end{bmatrix}
  \]
- canonical \(s=4\) 下：
  \[
  \eta_i=\chi_i=
  \begin{bmatrix}
  v_i\\ a_i\\ j_i
  \end{bmatrix}\in\mathbb R^3
  \]

不允许在代码中另起一套与论文冲突的命名体系。

---

## 3. 本脚本需要验证的数学命题

## 3.1 命题 A：raw BLOM-Analytic 的 waypoint Jacobian 稀疏模式

对 canonical interior segment \(i\)，Phase 6 的理论命题是：

\[
\frac{\partial c_i^{\mathrm{BLOM},2}}{\partial q_j}=0,
\qquad
j\notin\{i-2,i-1,i,i+1\}.
\]

即：

- 只允许对 \(q_{i-2},q_{i-1},q_i,q_{i+1}\) 非零敏感；
- 对其他 waypoint 必须严格为零（理论上）。

代码中必须把这一理论模式显式编码成 mask。

---

## 3.2 命题 B：raw BLOM-Analytic 的 duration Jacobian 稀疏模式

对 canonical interior segment \(i\)，Phase 6 的理论命题是：

\[
\frac{\partial c_i^{\mathrm{BLOM},2}}{\partial T_j}=0,
\qquad
j\notin\{i-1,i,i+1\}.
\]

即：

- 只允许对 \(T_{i-1},T_i,T_{i+1}\) 非零敏感；
- 对其余 duration 必须严格为零（理论上）。

---

## 3.3 命题 C：raw Jacobian 的 finite difference / 解析导数一致性

Phase 6 还给出了 canonical analytic map 的导数链条：

\[
x_i^{\mathrm{loc},\star}
=
A_2^{(i)}(T)^{-1}B_2^{(i)}(q,T),
\]
\[
\frac{\partial x_i^{\mathrm{loc},\star}}{\partial T_j}
=
A_2^{(i)}(T)^{-1}
\left(
\frac{\partial B_2^{(i)}}{\partial T_j}
-
\frac{\partial A_2^{(i)}}{\partial T_j}x_i^{\mathrm{loc},\star}
\right),
\]
\[
c_i^{\mathrm{BLOM},2}
=
\Lambda_8(T_i)C_4^{-1}y_i.
\]

因此，本脚本必须验证：

1. 有解析 Jacobian 的位置，用 finite difference 检查是否一致；
2. 有 autodiff 的位置，用 autodiff 与 finite difference 和解析结果三方对照；
3. 如果某一部分只有 finite difference，也要在文档中明确标注。

---

## 3.4 命题 D：装配后 Jacobian 是否仍保持局部

Phase 6 明确指出：

- raw BLOM / 方案 C：局部支撑最强；
- 方案 A / 方案 B：由于共享变量或 consensus 引入全局/准全局耦合，装配后的系数或最终轨迹对远处参数一般不再严格局部。

因此本脚本必须同时验证：

### 对方案 C

- raw coefficient map 的 Jacobian 稀疏模式是否正确；
- 若直接以 raw central extraction 为最终对象，则 Jacobian 是否仍保持局部带状。

### 对方案 A

- 最终装配后的 segment coefficient / knot-state 对远处 waypoint / duration 是否出现非零敏感度；
- 也就是说，要验证“Proposition 1 不再对装配后 Scheme A 成立”这一事实。

### 对方案 B

- consensus 后的对象是否比 raw 更非局部；
- 远处变量是否可通过重叠传播到某个 segment coefficient。

---

## 4. 核心验证对象

本脚本必须至少支持以下四类 Jacobian：

## 4.1 raw coefficient Jacobian（必须）

对每个 segment \(i\)，至少支持：

\[
J_{c,q}^{(i)}
:=
\frac{\partial c_i}{\partial \bar q},
\qquad
\bar q=(q_1,\dots,q_{M-1})^\top
\]
以及
\[
J_{c,T}^{(i)}
:=
\frac{\partial c_i}{\partial T}.
\]

### 输出形状

- \(c_i\in\mathbb R^{2s}\)
- `J_c_q_i`: shape `(2*s, M-1)`
- `J_c_T_i`: shape `(2*s, M)`

canonical \(s=4\) 时：
- `J_c_q_i`: shape `(8, M-1)`
- `J_c_T_i`: shape `(8, M)`

---

## 4.2 local jet-state Jacobian（推荐）

若实现方便，还应支持：

\[
J_{x,q}^{(i)}
:=
\frac{\partial x_i^{\mathrm{loc},\star}}{\partial \bar q},
\qquad
J_{x,T}^{(i)}
:=
\frac{\partial x_i^{\mathrm{loc},\star}}{\partial T}.
\]

canonical \(s=4,k=2\) 下
\[
x_i^{\mathrm{loc},\star}\in\mathbb R^6.
\]

输出形状：

- `J_x_q_i`: shape `(6, M-1)`
- `J_x_T_i`: shape `(6, M)`

这有助于定位误差到底来自 local linear solve 还是 Hermite reconstruction。

---

## 4.3 assembled knot-state Jacobian（方案 A/B，推荐）

### 方案 A

支持：
\[
\frac{\partial \eta^\star}{\partial \bar q},
\qquad
\frac{\partial \eta^\star}{\partial T}.
\]

### 方案 B

支持：
\[
\frac{\partial \eta^{\mathrm{consensus}}}{\partial \bar q},
\qquad
\frac{\partial \eta^{\mathrm{consensus}}}{\partial T}.
\]

这些 Jacobian 不要求带状，但必须能被统计和可视化。

---

## 4.4 final trajectory coefficient Jacobian（方案 A/B/C，推荐）

若实现可承受，建议还验证最终装配后：

\[
\frac{\partial c_i^{A}}{\partial \bar q},\quad
\frac{\partial c_i^{A}}{\partial T},
\qquad
\frac{\partial c_i^{B}}{\partial \bar q},\quad
\frac{\partial c_i^{B}}{\partial T},
\qquad
\frac{\partial c_i^{C}}{\partial \bar q},\quad
\frac{\partial c_i^{C}}{\partial T}.
\]

这是最直接展示“局部性 vs 连续性 trade-off”的数值证据。

---

## 5. finite difference 验证要求

## 5.1 必须实现统一 finite difference 接口

```python
finite_difference_jacobian(
    f, x, eps=1e-6, method="central", mask=None
) -> np.ndarray
```

### 要求

- 默认使用 central difference：
  \[
  \frac{\partial f}{\partial x_j}
  \approx
  \frac{f(x+\varepsilon e_j)-f(x-\varepsilon e_j)}{2\varepsilon}
  \]
- 支持 `mask`，只在需要检查的变量上做 perturbation
- 支持向量输出 \(f:\mathbb R^n\to\mathbb R^m\)

### 输出

shape `(m, n)` 的 Jacobian 近似

---

## 5.2 finite difference 步长要求

默认：

```python
eps_q = 1e-6
eps_T = 1e-6
```

但要允许 sweep 多个步长，例如：

- \(10^{-4}\)
- \(10^{-5}\)
- \(10^{-6}\)
- \(10^{-7}\)

并比较误差随步长变化的趋势。

这样可以避免单个 \(\varepsilon\) 选得不好导致误判。

---

## 5.3 finite difference 检查内容

### 对 raw BLOM-Analytic

必须检查：

1. `J_c_q_i` 的零模式是否与理论 mask 一致；
2. `J_c_T_i` 的零模式是否与理论 mask 一致；
3. 非零条目的 numerical value 是否与解析 Jacobian 接近。

### 对方案 A/B/C

至少检查：

1. 最终 `c_i` 对远处参数是否出现非零敏感度；
2. A/B 是否比 C 更“不局部”；
3. C 是否保持最窄的 Jacobian 带宽。

---

## 6. autodiff / symbolic 对照要求

## 6.1 首选方案

若实现环境允许，优先支持以下之一：

- `jax.jacrev`
- `autograd`
- `sympy.diff`（仅适合小尺寸 symbolic 检查）

## 6.2 对照要求

对 canonical \(s=4,k=2\) raw map：

- 解析 Jacobian vs finite difference
- autodiff Jacobian vs finite difference
- 解析 Jacobian vs autodiff Jacobian

至少输出三类误差：

\[
\|J_{\mathrm{analytic}}-J_{\mathrm{FD}}\|_\infty,
\quad
\|J_{\mathrm{AD}}-J_{\mathrm{FD}}\|_\infty,
\quad
\|J_{\mathrm{analytic}}-J_{\mathrm{AD}}\|_\infty.
\]

如果某种 autodiff 框架不适用，必须在 summary 中说明原因。

---

## 7. 方案 A / B / C 分别需要验证什么

## 7.1 方案 C（必须最完整）

因为 Phase 6 的局部支撑 Proposition 1 就是针对 raw BLOM-Analytic / Scheme C 的。

### 必须验证

1. `J_c_q_i` 的理论带状模式：
   \[
   j\notin\{i-2,i-1,i,i+1\}\Rightarrow \frac{\partial c_i}{\partial q_j}=0
   \]
2. `J_c_T_i` 的理论带状模式：
   \[
   j\notin\{i-1,i,i+1\}\Rightarrow \frac{\partial c_i}{\partial T_j}=0
   \]
3. 有解析 Jacobian 时，数值值是否正确；
4. 紧致 admissible 区间内，\(\partial c_i/\partial T_j\) 是否保持有界；
5. 多组随机测试下，带状模式是否稳定成立。

---

## 7.2 方案 A（必须做“反证式”验证）

理论上，装配后 Scheme A 一般不再保持严格局部支撑。

### 必须验证

1. 装配后 `c_i^A` 对远处 waypoint / duration 是否出现非零敏感度；
2. 这些非零敏感度是否显著大于数值噪声；
3. lower-order continuity 的获得是否以牺牲局部性为代价。

### 结论目标

用图和统计表明确展示：

- Scheme A 的 Jacobian 比 Scheme C 更“宽”
- 这是 Shared junction states 带来的合理后果

---

## 7.3 方案 B（必须做“传播范围”验证）

理论上，Scheme B 通过 consensus 会把局部窗口预测耦合起来，因此最终对象也会比 raw Scheme C 更不局部。

### 必须验证

1. consensus 前后的 Jacobian 带宽变化；
2. consensus 后远处参数对某个 segment coefficient 的影响是否出现；
3. 与 Scheme A 相比，它是否仍相对更“局部”。

### 结论目标

展示：

- Scheme B 比 C 更宽
- Scheme B 通常比 A 更保留 sliding-window 局部性（若数值上成立）

---

## 8. 必须保存的结果图片

所有图片必须保存到：

```text
results/phase6_fd_jacobian_check/
```

并按方案组织目录：

```text
results/phase6_fd_jacobian_check/
├── raw_scheme_C/
├── scheme_A/
├── scheme_B/
└── compare/
```

## 8.1 图 1：Jacobian 稀疏模式热图（必须）

文件名建议：

```text
jacobian_mask_q_raw.png
jacobian_mask_T_raw.png
jacobian_mask_q_scheme_A.png
jacobian_mask_T_scheme_A.png
jacobian_mask_q_scheme_B.png
jacobian_mask_T_scheme_B.png
```

### 内容

横轴：参数索引  
纵轴：系数索引  
颜色：\(\log_{10}(|J_{ab}|+\varepsilon)\)

### 目的

直观看理论带状模式是否成立，以及 A/B 的带宽是否更宽。

---

## 8.2 图 2：理论 mask vs 数值 Jacobian 对照图（必须）

文件名建议：

```text
jacobian_theory_vs_fd_q_raw.png
jacobian_theory_vs_fd_T_raw.png
```

### 内容

左右并排：

- 左：理论零/非零 mask
- 右：finite difference 的实际绝对值热图

### 目的

直接证明 Proposition 1 的稀疏模式。

---

## 8.3 图 3：解析 Jacobian 与 FD 误差热图（必须）

文件名建议：

```text
jacobian_error_analytic_vs_fd_q.png
jacobian_error_analytic_vs_fd_T.png
```

### 内容

\[
|J_{\mathrm{analytic}}-J_{\mathrm{FD}}|
\]
的热图。

### 目的

直接展示解析公式正确性。

---

## 8.4 图 4：误差随步长变化曲线（必须）

文件名建议：

```text
fd_stepsize_sweep_q.png
fd_stepsize_sweep_T.png
```

### 内容

横轴：\(\varepsilon\)  
纵轴：Jacobian 误差范数

### 目的

说明 finite difference 验证不是偶然结果。

---

## 8.5 图 5：方案 A / B / C 局部性宽度对比图（必须）

文件名建议：

```text
jacobian_bandwidth_compare.png
```

### 内容

对每个方案统计“超出理论 local support 邻域的非零条目比例”或“有效带宽”。

### 目的

直观证明：

- C 最局部
- A / B 更宽
- 这是 continuity 与 locality 的 trade-off

---

## 8.6 图 6：\(\partial c / \partial T\) 有界性统计图（推荐）

文件名建议：

```text
dc_dT_bound_statistics.png
```

### 内容

在随机样本下，统计
\[
\left\|\frac{\partial c_i}{\partial T_j}\right\|
\]
的最大值、分位数、箱线图。

### 目的

支撑 Phase 6 Proposition 4 的“紧致 admissible 集合上有界”。

---

## 8.7 图 7：random trials 误差分布箱线图（推荐）

文件名建议：

```text
jacobian_error_boxplot_random_trials.png
```

### 内容

多组随机试验中：
- 解析 vs FD
- AD vs FD
- 解析 vs AD

误差分布。

### 目的

证明验证结果不是单个案例巧合。

---

## 9. 必须保存的表和数据文件

## 9.1 Jacobian 稀疏统计 CSV（必须）

文件名建议：

```text
jacobian_sparsity_stats_raw.csv
jacobian_sparsity_stats_scheme_A.csv
jacobian_sparsity_stats_scheme_B.csv
```

列至少包含：

- `scheme`
- `segment_idx`
- `jacobian_type` (`c_q`, `c_T`, `x_q`, `x_T`)
- `nnz_total`
- `nnz_theory_band`
- `nnz_outside_band`
- `max_abs_outside_band`
- `mean_abs_outside_band`

---

## 9.2 Jacobian 误差统计 CSV（必须）

文件名建议：

```text
jacobian_error_stats.csv
```

列至少包含：

- `scheme`
- `segment_idx`
- `jacobian_type`
- `error_metric`
- `analytic_vs_fd`
- `ad_vs_fd`
- `analytic_vs_ad`
- `eps`

---

## 9.3 总结 JSON（必须）

文件名建议：

```text
summary_phase6_fd_check.json
```

必须包含：

- 输入 `q,T,s,k`
- 方案名称
- 检查了哪些 Jacobian
- 理论 mask 是否通过
- outside-band 最大误差
- 解析 vs FD 是否通过
- autodiff 是否可用
- 运行时间
- 图路径
- 文件路径

---

## 10. 推荐代码结构

```text
phase_6/
├── blom_fd_jacobian_check.py
├── test_blom_fd_jacobian_check.py
├── examples/
│   ├── demo_raw_local_support.py
│   ├── demo_scheme_A_nonlocality.py
│   ├── demo_scheme_B_nonlocality.py
│   └── demo_compare_all.py
└── results/
    └── phase6_fd_jacobian_check/
```

`blom_fd_jacobian_check.py` 至少包含：

```python
def theoretical_mask_c_q(M, i, s=4, k=2): ...
def theoretical_mask_c_T(M, i, s=4, k=2): ...

def finite_difference_jacobian(f, x, eps=1e-6, method="central", mask=None): ...
def compute_jacobian_errors(J_ref, J_test): ...

def raw_local_coefficient_map(q, T, i, s=4, k=2): ...
def raw_local_jacobians(q, T, i, s=4, k=2, mode="analytic"): ...

def assembled_scheme_A_map(q, T, i, config=None): ...
def assembled_scheme_B_map(q, T, i, config=None): ...
def assembled_scheme_C_map(q, T, i, config=None): ...

def run_fd_jacobian_check(q, T, s=4, k=2, scheme="C",
                          use_autodiff=True,
                          save_dir=None,
                          seed=42): ...

def run_compare_all_schemes(q, T, s=4, k=2,
                            use_autodiff=True,
                            save_dir=None,
                            seed=42): ...

def run_random_trials(n_trials=100, M=20, s=4, k=2,
                      use_autodiff=True,
                      save_dir=None,
                      seed=42): ...
```

---

## 11. 随机测试要求

必须支持多组随机 \(q,T\) 自动测试。

建议：

```python
run_random_trials(
    n_trials=100,
    M=20,
    s=4,
    k=2,
    use_autodiff=True,
    save_dir=None,
    seed=42
)
```

### 默认采样建议

- `q`：高斯随机游走
- `T`：独立正值样本，例如
  \[
  T_i\sim \mathrm{Uniform}(0.5,2.0)
  \]

### 必须统计

- raw Jacobian 的理论 mask 通过率
- 解析 Jacobian vs FD 的通过率
- 方案 A / B / C 的 outside-band 非零率
- \(\partial c/\partial T\) 的最大范数分布
- 平均运行时间

---

## 12. 验收标准

只有满足以下条件，本脚本才算完成：

1. 与 Phase 0--6 的符号和数学对象一致；
2. 能对 raw BLOM-Analytic 的 `J_c_q`, `J_c_T` 做有限差分验证；
3. 能显式检查理论局部支撑 mask；
4. 能保存必须的热图、误差图、对比图；
5. 能输出 CSV / JSON 统计结果；
6. 能区分 Scheme A / B / C 的 Jacobian 局部性差异；
7. 在 canonical \(s=4,k=2\) 下运行稳定；
8. 能批量随机测试；
9. 结果足以直接用于论文中支撑 Proposition 1 和 Proposition 4。

---

## 13. 给实现 AI 的最终指令

你要实现的是一个**Phase 6 Jacobian 基础性质验证脚本**。

请优先保证：

1. 理论稀疏 mask 构造正确；
2. finite difference 计算稳定；
3. 解析 Jacobian 与 FD 对照清楚；
4. 方案 A / B / C 的对比逻辑清楚；
5. 所有关键图和表自动保存；
6. 结果结构统一，便于直接写论文；
7. 必须把“C 最局部、A/B 更非局部”这一核心结论数值化。

不要把 raw local map 与 assembled scheme 混为一谈。
本脚本的核心任务是：**用数值证据把 Phase 6 的局部支撑与弱时间感知命题落地，并同时说明不同装配方案对 Jacobian 局部性的影响。**
