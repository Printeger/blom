# REQ-Phase7-convergence-validation.md

## 1. 文档目的

本需求文档用于指导实现 `blom_convergence_vs_k.py`，作为 BLOM 理论 **Phase 7：收敛理论数值验证版** 的统一实验脚本。

该脚本的目标不是直接证明所有收敛定理，而是用系统性的数值证据回答以下问题：

1. actual BLOM-\(k\) 的误差是否随 \(k\) 近似指数衰减；
2. 哪种装配方案（A / B / C）更接近全局 MINCO；
3. actual BLOM-\(k\) 与 idealized truncated BLOM-\(k\) 是否数值上接近，即 matching 看起来是否成立；
4. 代价误差是否也呈现随 \(k\) 衰减的趋势；
5. 根据数值结果，后续理论应如何分叉：
   - 若数值非常漂亮，则继续补：
     - matching theorem
     - feasibility / cost theorem
   - 若数值一般，甚至不支持当前想法，则不要硬补 theorem，而应回头修改：
     - BLOM 定义
     - 装配方案
     - 比较对象
     - 主 claim

本脚本应与前面各 phase 保持接口兼容，尤其应兼容：

- `minco_scalar_baseline.py`
- `blom_strict_local_qp.py`
- `blom_k2_s4_numeric.py`
- `blom_boundary_jump_check.py`
- `blom_fd_jacobian_check.py`

---

## 2. 与前面 Phases 统一的符号要求

必须与 Phase 0--7 完全统一：

- 轨迹：
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
- Phase 1 全局 MINCO 最优系数：
  \[
  c^\star(q,T)\in\mathbb R^{2sM}
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
- 误差范数默认使用：
  \[
  \| \cdot \|_2
  \]

禁止在代码中换成与前文不兼容的新符号系统。

---

## 3. 本脚本需要区分的三个核心对象

### 3.1 全局 MINCO 参考解（reference baseline）

这是 Phase 1 的全局最优系数：
\[
c^\star(q,T).
\]

### 3.2 理想截断 BLOM-\(k\)（idealized truncated BLOM-\(k\)）

这是 Phase 7 理论里先证明收敛的对象：
\[
\tilde c_i^{(k)}
=
c_i^\star(0,T)
+
\sum_{j\in\mathcal N_k(i)} G_{ij}(T)\,q_j.
\]

### 3.3 actual BLOM-\(k\)

这是你真正定义并实现的 BLOM-\(k\)，来自：

- Phase 3：BLOM-Strict 局部问题
- Phase 4：BLOM-Analytic 局部解析系统
- Phase 5：装配方案 A / B / C

---

## 4. 本脚本必须回答的核心数值问题

### 4.1 误差是否随 \(k\) 近似指数衰减

对固定一组 \(q,T\)，必须计算：
\[
k=2,4,6,\dots,k_{\max}
\]
对应的误差，例如：
\[
E_{\mathrm{actual}}(k):=\|c^{\mathrm{BLOM},k}-c^\star\|_2,
\]
\[
E_{\mathrm{ideal}}(k):=\|\tilde c^{(k)}-c^\star\|_2.
\]

然后拟合：
\[
\log E(k)
\]
与 \(k\) 的线性关系。

### 4.2 actual BLOM-\(k\) 与 idealized truncated BLOM-\(k\) 是否接近

必须计算 matching 误差：
\[
E_{\mathrm{match}}(k):=
\|c^{\mathrm{BLOM},k}-\tilde c^{(k)}\|_2.
\]

### 4.3 哪种装配方案更接近 MINCO

对方案 A / B / C 分别计算：

\[
E_A(k):=\|c^{A,k}-c^\star\|_2,
\qquad
E_B(k):=\|c^{B,k}-c^\star\|_2,
\qquad
E_C(k):=\|c^{C,k}-c^\star\|_2.
\]

### 4.4 代价误差是否也呈现衰减趋势

必须计算：

\[
\Delta J_{\mathrm{abs}}(k)
:=
J(c^{(k)},T)-J(c^\star,T),
\]
\[
\Delta J_{\mathrm{rel}}(k)
:=
\frac{J(c^{(k)},T)}{J(c^\star,T)}-1.
\]

---

## 5. k 序列的选取要求

必须支持：

```python
k_values = [2, 4, 6, ..., k_max]
```

建议接口：

```python
make_k_grid(M, start=2, step=2, include_M=True)
```

---

## 6. 必须计算的误差指标

### 6.1 系数误差（必须）

\[
E_{\mathrm{global}}(k)
=
\|c^{(k)}-c^\star\|_2
\]

### 6.2 matching 误差（必须）

\[
E_{\mathrm{match}}(k)
=
\|c^{\mathrm{BLOM},k}-\tilde c^{(k)}\|_2
\]

### 6.3 代价误差（必须）

\[
\Delta J_{\mathrm{abs}}(k)
=
J(c^{(k)},T)-J(c^\star,T),
\]
\[
\Delta J_{\mathrm{rel}}(k)
=
\frac{J(c^{(k)},T)}{J(c^\star,T)}-1.
\]

---

## 7. 必须实现的核心函数

### 7.1 reference MINCO 解

```python
compute_minco_reference(q, T, s=4) -> dict
```

### 7.2 理想截断 BLOM-\(k\)

```python
compute_ideal_truncated_blom_k(q, T, k, s=4) -> dict
```

### 7.3 actual BLOM-\(k\)

```python
compute_actual_blom_k(q, T, k, s=4, scheme="C", config=None) -> dict
```

### 7.4 误差计算

```python
compute_convergence_errors(c_ref, c_test, c_ideal=None, cost_ref=None, cost_test=None) -> dict
```

### 7.5 主运行器

```python
run_convergence_vs_k(
    q, T, s=4,
    k_values=None,
    schemes=("A", "B", "C"),
    save_dir=None,
    seed=42
) -> dict
```

---

## 8. 拟合要求：如何画 \(\log\|c^{(k)}-c^\star\|\) 对 \(k\) 的图

### 8.1 必须拟合的对象

至少包括：

1. ideal truncation:
   \[
   \log E_{\mathrm{ideal}}(k)
   \]
2. actual BLOM-\(k\), A:
   \[
   \log E_A(k)
   \]
3. actual BLOM-\(k\), B:
   \[
   \log E_B(k)
   \]
4. actual BLOM-\(k\), C:
   \[
   \log E_C(k)
   \]

### 8.2 拟合方式

默认用最小二乘直线拟合：
\[
\log E(k)\approx \alpha + \beta k.
\]

必须输出：

- slope \(\beta\)
- intercept \(\alpha\)
- \(R^2\)
- 拟合区间

---

## 9. 必须保存的结果图片

所有图片必须保存到：

```text
results/phase7_convergence_vs_k/
```

并分为：

```text
results/phase7_convergence_vs_k/
├── ideal/
├── scheme_A/
├── scheme_B/
├── scheme_C/
└── compare/
```

### 图 1：全局误差 vs \(k\) 曲线（必须）

文件名建议：

```text
coef_error_vs_k_all.png
```

### 图 2：\(\log\) 误差 vs \(k\) 拟合图（必须）

文件名建议：

```text
log_coef_error_fit_all.png
```

### 图 3：matching 误差 vs \(k\)（必须）

文件名建议：

```text
matching_error_vs_k.png
```

### 图 4：relative cost gap vs \(k\)（必须）

文件名建议：

```text
relative_cost_gap_vs_k.png
```

### 图 5：每段误差热图（推荐）

文件名建议：

```text
segmentwise_error_heatmap_scheme_A.png
segmentwise_error_heatmap_scheme_B.png
segmentwise_error_heatmap_scheme_C.png
```

### 图 6：random trials 下拟合 slope 分布图（推荐）

文件名建议：

```text
logfit_slope_boxplot_random_trials.png
```

### 图 7：ideal vs actual 对比图（推荐）

文件名建议：

```text
ideal_vs_actual_compare.png
```

---

## 10. 必须保存的表和数据文件

### 10.1 每个 \(k\) 的误差明细 CSV（必须）

文件名建议：

```text
convergence_errors_by_k.csv
```

### 10.2 拟合结果 CSV（必须）

文件名建议：

```text
logfit_summary.csv
```

### 10.3 总结 JSON（必须）

文件名建议：

```text
summary_phase7_convergence.json
```

### 10.4 自动解释性 summary（必须）

文件名建议：

```text
phase7_interpretation_summary.md
```

内容至少包括：

1. actual BLOM-\(k\) 是否近似指数衰减；
2. 哪种装配方案更接近 MINCO；
3. matching 看起来是否成立；
4. 后续理论建议：继续补 theorem 还是改模型。

---

## 11. 随机测试要求

必须支持多组随机 \(q,T\) 自动测试。

建议函数：

```python
run_random_trials(
    n_trials=100,
    M=20,
    s=4,
    k_values=None,
    schemes=("A", "B", "C"),
    save_dir=None,
    seed=42
)
```

### 默认采样建议

- `q`：高斯随机游走
- `T`：正值随机，例如
  \[
  T_i\sim \mathrm{Uniform}(0.5, 2.0)
  \]

### 必须统计

- 各方案成功率
- 各方案 slope 分布
- 各方案 \(R^2\) 分布
- matching 误差分布
- relative cost gap 分布

---

## 12. 推荐代码结构

```text
phase_7/
├── blom_convergence_vs_k.py
├── test_blom_convergence_vs_k.py
├── examples/
│   ├── demo_single_case.py
│   ├── demo_compare_schemes.py
│   ├── demo_ideal_vs_actual.py
│   └── demo_random_trials.py
└── results/
    └── phase7_convergence_vs_k/
```

`blom_convergence_vs_k.py` 至少应包含：

```python
def make_k_grid(M, start=2, step=2, include_M=True): ...

def compute_minco_reference(q, T, s=4): ...
def compute_ideal_truncated_blom_k(q, T, k, s=4): ...
def compute_actual_blom_k(q, T, k, s=4, scheme="C", config=None): ...

def compute_convergence_errors(c_ref, c_test, c_ideal=None,
                               cost_ref=None, cost_test=None): ...

def fit_log_error_vs_k(k_values, errors, eps=1e-15): ...

def run_convergence_vs_k(q, T, s=4, k_values=None,
                         schemes=("A", "B", "C"),
                         save_dir=None, seed=42): ...

def run_random_trials(n_trials=100, M=20, s=4,
                      k_values=None,
                      schemes=("A", "B", "C"),
                      save_dir=None, seed=42): ...
```

---

## 13. 验收标准

只有满足以下条件，本脚本才算完成：

1. 与 Phase 0--7 的数学符号统一；
2. 能构造全局 MINCO 参考解；
3. 能区分 idealized truncated BLOM-\(k\) 与 actual BLOM-\(k\)；
4. 能对 \(k=2,4,6,\dots,M\) 自动计算误差；
5. 能画 \(\log\|c^{(k)}-c^\star\|\) vs \(k\) 并给线性拟合；
6. 能分别比较方案 A/B/C；
7. 能保存必须的图、表、JSON、summary；
8. 能批量随机测试；
9. 能自动输出“后续理论该如何分叉”的建议。

---

## 14. 给实现 AI 的最终指令

你要实现的是一个**Phase 7 收敛理论数值验证脚本**，不是单个示意图脚本。

请优先保证：

1. global MINCO reference 正确；
2. ideal truncation 与 actual BLOM-\(k\) 明确区分；
3. 误差指标定义统一；
4. A/B/C 三种方案都能跑；
5. \(\log\) 拟合图和 matching 图必须自动保存；
6. 结论摘要必须自动生成；
7. 随机试验必须支持。

不要把“数值上看起来像指数衰减”直接写成 theorem。
本脚本的核心任务是：**用数据告诉我们 actual BLOM-\(k\) 是否真的像理论希望的那样工作，并据此决定下一步该补 theorem 还是改模型。**
