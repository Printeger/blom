# REQ-Phase4-analytic-validation.md

## 1. 文档目的

本需求文档用于指导实现并验证 Phase 4 的三个核心脚本：

- `blom_k2_s2_sympy.py`
- `blom_k2_s4_numeric.py`
- `blom_catmull_compare.py`

其目标不是直接构建完整 BLOM 规划器，而是把 **Phase 4 的数学结论** 变成可检查、可复现实验代码，服务于：

1. 验证 Phase 4 推导是否正确；
2. 为后续 BLOM-Analytic 的实现提供统一接口；
3. 为论文和汇报生成直观结果图；
4. 与 Phase 1 的全局 `minco_scalar_baseline.py`、Phase 2 的必要性验证代码保持输入输出兼容。

本阶段的代码任务分为两层：

- **符号层验证**：确认最小模型 `s=2, k=2` 的解析推导严格正确；
- **数值层验证**：在 canonical case `s=4, k=2` 下搭建局部解析系统，并验证它与局部 BLOM-Strict 问题一致；
- **对比层验证**：定量比较时间加权 Catmull--Rom 与精确局部解析系统在 `s=4` 下的差异。

---

## 2. 与 Phase 0 / Phase 1 / Phase 2 / Phase 3 的统一要求

实现必须保持以下统一约定：

### 2.1 外层参数与输入记号

- 全局航点：
  \[
  q_0, q_1, \dots, q_M \in \mathbb R
  \]
- 时间向量：
  \[
  T=(T_1,\dots,T_M)^\top, \qquad T_i > 0
  \]
- 总时长：
  \[
  \mathcal T = \sum_{i=1}^{M} T_i
  \]
- 外层稀疏参数：
  \[
  \theta = (q_1,\dots,q_{M-1}, T_1,\dots,T_M)^\top
  \]
- 当前只做 **标量轨迹**。

### 2.2 窗口记号

对中央段 `i` 的三段窗口：
\[
W(i,2)=\{i-1,i,i+1\}
\]
对应 knot 数据：
\[
q_{i-2}, q_{i-1}, q_i, q_{i+1}
\]
对应时间数据：
\[
T_{i-1}, T_i, T_{i+1}
\]
并采用简写：
\[
h_-:=T_{i-1}, \qquad h_0:=T_i, \qquad h_+:=T_{i+1}.
\]

### 2.3 Phase 4 中的局部未知量

- 在 `s=2, k=2` 中，局部未知量为单个 knot 速度：
  \[
  v_i
  \]
- 在 `s=4, k=2` 中，局部未知量为中央段两端 knot jet：
  \[
  x_i^{\mathrm{loc}} = [v_{i-1}, a_{i-1}, j_{i-1}, v_i, a_i, j_i]^\top \in \mathbb R^6.
  \]

### 2.4 与 Phase 1 / 2 / 3 的接口兼容要求

- 脚本必须能读取或复用 `phase_1/minco_scalar_baseline.py` 中的 waypoint / duration 数据结构；
- 若 Phase 3 的 `blom_strict_local_qp.py` 已实现，则 `blom_k2_s4_numeric.py` 必须提供与其可对比的局部目标值和最优解接口；
- 输出格式必须允许后续拼接进 BLOM-Analytic 的局部映射：
  \[
  x \to c_i.
  \]

---

## 3. 本阶段要验证的数学结论

## 3.1 Script A：`blom_k2_s2_sympy.py`

### 目标

对最小可解析模型 `s=2, k=2` 做**严格符号推导验证**，确认以下公式正确：

1. 单段 one-sided natural cubic 的解析形式；
2. 左段和右段代价写成关于公共速度 `v_i` 的一维二次函数；
3. 一阶最优性条件可化为
   \[
   A_2^{(i)}(T) v_i = B_2^{(i)}(q,T)
   \]
4. 精确最优速度为
   \[
   v_i^* = \frac{h_+ s_- + h_- s_+}{h_- + h_+},
   \]
   其中
   \[
   s_- = \frac{q_i-q_{i-1}}{T_i}, \qquad s_+ = \frac{q_{i+1}-q_i}{T_{i+1}}.
   \]
5. 该公式与时间加权 Catmull--Rom 完全一致。

### 必须完成的数学任务

1. 用 SymPy 建立三次多项式：
   \[
   p(t)=a_0+a_1 \tau + a_2 \tau^2 + a_3 \tau^3, \qquad \tau=t/h.
   \]
2. 显式施加约束：
   - 左 natural: \(p''(0)=0\)
   - 右端插值: \(p(h)=q_R\)
   - 右端速度: \(p'(h)=v_R\)
3. 求解系数 `a0,a1,a2,a3`；
4. 计算
   \[
   J = \int_0^h |p''(t)|^2 dt
   \]
   并化简；
5. 拼成两段局部总代价；
6. 对 `v_i` 求导并解最优点；
7. 证明结果等于时间加权 Catmull--Rom 公式。

### 代码要求

脚本应包含以下函数：

```python
build_left_natural_cubic_symbolic()
build_right_natural_cubic_symbolic()
derive_s2_local_energy_symbolic()
derive_s2_exact_velocity_symbolic()
verify_catmull_equivalence_symbolic()
main()
```

### 输出要求

必须输出：

1. 符号化的多项式系数；
2. 单段代价表达式；
3. 两段总代价表达式；
4. 最优速度闭式表达；
5. 与 Catmull--Rom 公式差值化简结果（应为 0）。

### 必须保存的结果文件

保存在：

```text
phase_4/results/s2_sympy/
```

必须包含：

- `s2_left_segment_coeffs.txt`
- `s2_right_segment_coeffs.txt`
- `s2_energy_expression.txt`
- `s2_optimal_velocity_expression.txt`
- `s2_catmull_equivalence.txt`

### 必须保存的结果图片

脚本应生成并保存以下图片：

1. `s2_energy_vs_v.png`
   - 固定一个代表性例子，绘制局部代价关于 `v_i` 的曲线；
   - 在图中标出最优 `v_i^*`；
2. `s2_catmull_match_demo.png`
   - 对若干随机样本，比较解析最优速度与 Catmull 速度，画散点图 `v_exact` vs `v_catmull`；
   - 理想结果应落在对角线上。

### 验收标准

1. SymPy 化简后，`v_exact - v_catmull == 0`；
2. 代价函数关于 `v_i` 的二次项系数严格为正；
3. 数值随机样本中，解析值与 Catmull 值误差在机器精度范围内。

---

## 3.2 Script B：`blom_k2_s4_numeric.py`

### 目标

实现 canonical case `s=4, k=2` 的**精确局部解析系统数值版**，并验证：

1. Phase 4 中的局部目标可以写成
   \[
   J_i^{(4,2)}(x) = \frac12 x^\top A_2^{(i)}(T)x - x^\top B_2^{(i)}(q,T) + C_2^{(i)}(q,T)
   \]
2. 最优解满足
   \[
   A_2^{(i)}(T)x_i^{\mathrm{loc},*} = B_2^{(i)}(q,T)
   \]
3. `A_2^{(i)}(T)` 对 admissible 输入是对称正定；
4. 得到的最优 `x_i^{loc,*}` 通过 Hermite 重构后，中央段系数与局部 BLOM-Strict QP 的最优中央段一致。

### 必须完成的数学对象

#### (a) one-sided outer segment cost

实现左 outer segment 代价：
\[
J_{4,\mathrm{L}}(x_R;q_L,q_R,h)=252h^{-7}(\Delta - m_-(h)^\top x_R)^2
\]
和右 outer segment 代价：
\[
J_{4,\mathrm{R}}(x_L;q_L,q_R,h)=252h^{-7}(\Delta - m_+(h)^\top x_L)^2.
\]
其中
\[
m_-(h)=\begin{bmatrix}h\\ -h^2/2\\ h^3/6\end{bmatrix},
\qquad
m_+(h)=\begin{bmatrix}h\\ h^2/2\\ h^3/6\end{bmatrix}.
\]

#### (b) middle segment Hermite energy

实现中央段 exact Hermite energy：
\[
J_{4,\mathrm{mid}}(x_L,x_R;q_L,q_R,h)=h^{-7} y^\top R_4 y
\]
其中：

- `C4`：归一化 Hermite 约束矩阵；
- `H4 = C4^{-1}`；
- `G4`：归一化四阶导 Gram 矩阵；
- `R4 = H4^T G4 H4`；
- `y = [q_L, D(h)x_L, q_R, D(h)x_R]^T`；
- `D(h)=diag(h,h^2,h^3)`。

#### (c) total local quadratic

将三段窗口总代价拼成：
\[
J_i^{(4,2)}(x_i^{loc}) = J_{4,L}+J_{4,mid}+J_{4,R}
\]
并从中自动提取：

- `A2_i(T)`：6x6 Hessian
- `B2_i(q,T)`：6x1 线性项
- `C2_i(q,T)`：常数项

#### (d) exact local solve

实现
\[
x_i^{loc,*} = A_2^{(i)}(T)^{-1} B_2^{(i)}(q,T)
\]
并检查：

- `A2_i` 对称；
- `A2_i` 正定；
- 梯度在最优点为 0。

#### (e) Hermite reconstruction

实现中央段系数重构：
\[
c_i^{\mathrm{BLOM},2} = \Lambda_8(T_i) C_4^{-1} y_i
\]
其中
\[
\Lambda_8(T_i)=\operatorname{diag}(1,T_i^{-1},\dots,T_i^{-7}).
\]

### 代码要求

脚本应包含以下函数：

```python
build_m_plus(h)
build_m_minus(h)
left_outer_cost_matrix(qL, qR, h)
right_outer_cost_matrix(qL, qR, h)
build_C4()
build_G4()
build_R4()
middle_segment_cost_matrix(qL, qR, h)
build_local_quadratic_s4_k2(q_im2, q_im1, q_i, q_ip1, T_im1, T_i, T_ip1)
solve_local_system_s4_k2(...)
hermite_reconstruct_center_segment(...)
compare_with_blom_strict_qp(...)
main()
```

### 与 Phase 3 的接口要求

如果 `phase_3/blom_strict_local_qp.py` 存在，则必须提供一个对比函数：

```python
compare_with_blom_strict_qp(...)
```

比较以下量：

1. analytic local solution `x_exact`
2. local QP solution `x_qp`
3. central segment coefficients `c_exact`
4. QP 中央段系数 `c_qp`
5. local objective values

### 输出要求

必须输出：

- `A2_i(T)` 数值矩阵；
- `B2_i(q,T)` 数值向量；
- `x_i^{loc,*}`；
- `c_i` 中央段 8 个系数；
- 与 BLOM-Strict QP 的误差对比。

### 必须保存的结果文件

保存在：

```text
phase_4/results/s4_numeric/
```

必须包含：

- `A2_matrix.npy`
- `B2_vector.npy`
- `x_local_opt.npy`
- `center_coeffs.npy`
- `R4_matrix.npy`
- `C4_matrix.npy`
- `G4_matrix.npy`
- `analytic_vs_qp_metrics.json`

### 必须保存的结果图片

必须生成：

1. `A2_heatmap.png`
   - 画出 `A2_i(T)` 的热力图；
2. `A2_eigenvalues.png`
   - 画出若干随机样本下最小特征值 / 全部特征值，展示正定性；
3. `analytic_vs_qp_error_bar.png`
   - 展示 `||x_exact-x_qp||` 和 `||c_exact-c_qp||`；
4. `center_segment_overlay.png`
   - 对某个代表性样本，将 analytic 重构中央段和 QP 中央段在同一图中叠加；
5. `local_energy_landscape_2d.png`
   - 固定其余变量，选取两个代表分量，对局部代价画二维切片等高线图，并标记最优点。

### 验收标准

1. `A2_i` 必须数值上对称；
2. `A2_i` 必须正定（随机样本下最小特征值严格正）；
3. `A2_i x_exact - B2_i` 残差足够小；
4. analytic 与 local QP 的 `x` 和 `c` 误差足够小；
5. Hermite 重构中央段满足两端位置、速度、加速度、jerk 的精确插值。

---

## 3.3 Script C：`blom_catmull_compare.py`

### 目标

验证并直观展示：

1. 在 `s=2, k=2` 中，时间加权 Catmull--Rom 与精确局部解完全一致；
2. 在 `s=4, k=2` 中，时间加权 Catmull--Rom 不再是精确局部最优解；
3. 这种偏差如何随 waypoint 几何关系和时间非均匀性变化。

### 要比较的对象

#### (a) `s=2`
- `v_catmull`
- `v_exact_s2`

#### (b) `s=4`
构造一个“Catmull-like heuristic jet”基线。至少要求：
- 用 Catmull 时间加权公式给出速度；
- 加速度、jerk 可设为 0，或采用一个清楚说明的启发式补全方案；
- 将该 heuristic jet 送入 central Hermite reconstruction，形成 heuristic 中央段；
- 与 exact analytic local system 得到的 `x_exact`、`c_exact` 比较。

### 代码要求

脚本应包含以下函数：

```python
catmull_velocity_time_weighted(q_im1, q_i, q_ip1, T_i, T_ip1)
heuristic_local_state_s4(...)
compare_s2_exact_vs_catmull(...)
compare_s4_exact_vs_catmull(...)
run_random_benchmark(...)
plot_comparison_figures(...)
main()
```

### 数值实验要求

至少包含三类实验：

#### Case A：均匀时间
- `T_{i-1}=T_i=T_{i+1}`
- 检查 `s=4` 下 heuristic 是否仍与 exact 有偏差；

#### Case B：非均匀时间
- 随机采样显著不均匀时长；
- 这是最关键的对比情形；

#### Case C：局部几何变化
- 改变 waypoint 的曲率模式，例如：
  - 近线性点列；
  - 单调上升；
  - 尖锐折线；
  - 非对称转折；
- 比较 heuristic 与 exact 的差异大小。

### 输出要求

必须输出：

- `s2` 中 exact 与 Catmull 的误差统计（应接近 0）；
- `s4` 中 heuristic 与 exact 的误差统计；
- 中央段系数误差、轨迹误差、局部目标值差异。

### 必须保存的结果文件

保存在：

```text
phase_4/results/catmull_compare/
```

必须包含：

- `s2_compare_metrics.json`
- `s4_compare_metrics.json`
- `random_benchmark_table.csv`

### 必须保存的结果图片

必须生成：

1. `s2_exact_vs_catmull_scatter.png`
   - `v_exact` vs `v_catmull`；
2. `s4_velocity_exact_vs_catmull_scatter.png`
   - `v_exact` vs `v_catmull`；
3. `s4_coeff_error_hist.png`
   - `||c_exact-c_heuristic||` 的直方图；
4. `s4_objective_gap_hist.png`
   - `J_heuristic - J_exact` 的直方图；
5. `representative_curve_overlay.png`
   - 选一个非均匀时间样本，叠加 heuristic 和 exact 的中央段轨迹；
6. `error_vs_time_imbalance.png`
   - 横轴为时间不均匀程度指标，纵轴为 heuristic-exact 偏差。

### 验收标准

1. `s=2` 时，Catmull 与 exact 完全一致；
2. `s=4` 时，存在非零差异，且至少在部分非均匀时间样本下明显可见；
3. exact analytic 的目标值不高于 heuristic 目标值；
4. 图像和统计结果足以支持论文中的一句核心结论：
   > 时间加权 Catmull--Rom 只在 `s=2` 情形下是 exact local optimum，在 `s=4` 中必须被真正的局部解析系统替代。

---

## 4. 统一目录结构要求

建议目录：

```text
phase_4/
├── blom_k2_s2_sympy.py
├── blom_k2_s4_numeric.py
├── blom_catmull_compare.py
├── utils/
│   ├── hermite_utils.py
│   ├── plotting_utils.py
│   └── io_utils.py
├── results/
│   ├── s2_sympy/
│   ├── s4_numeric/
│   └── catmull_compare/
└── README_phase4.md
```

---

## 5. 统一输入输出格式要求

### 输入

所有脚本统一使用：

- waypoint 标量数组 `q`
- duration 数组 `T`
- 中央段索引 `i`
- 随机种子 `seed`

其中 `s=4, k=2` 脚本中，窗口输入统一裁剪为：

```python
q_local = [q[i-2], q[i-1], q[i], q[i+1]]
T_local = [T[i-1], T[i], T[i+1]]
```

### 输出

统一输出一个字典：

```python
{
    "inputs": ...,
    "matrices": ...,
    "solution": ...,
    "metrics": ...,
    "fig_paths": ...,
}
```

这样便于后续 BLOM 总实验或论文图表脚本统一读取。

---

## 6. 统一绘图风格要求

所有图应满足：

1. 白底；
2. 坐标轴、标题、图例完整；
3. 文件格式优先 `png`，必要时同时保存 `pdf`；
4. 文件名固定、可复现；
5. 图中必须标出样本参数或关键结论；
6. 不要使用过度花哨配色，保证论文插图可读性。

---

## 7. 实现限制与注意事项

1. 先保证**数学对象正确**，再考虑工程优化；
2. `s=2` 脚本以符号验证为主，不追求性能；
3. `s=4` 脚本以数值稳定和可检查为主，不要求第一版就写最优实现；
4. 大矩阵条目不要手工硬编码，应通过基函数 / Gram 矩阵 / Hermite 矩阵自动生成；
5. 所有关键步骤都必须有清晰注释，对应到 Phase 4 中的具体公式；
6. 如果有与 Phase 3 local QP 的偏差，必须保存误差信息，不允许静默忽略。

---

## 8. 最终验收总标准

只有当以下条件全部满足时，Phase 4 验证代码才算完成：

1. `blom_k2_s2_sympy.py` 严格推出并验证 Catmull 等价公式；
2. `blom_k2_s4_numeric.py` 成功构造并求解 6x6 精确局部系统；
3. `blom_k2_s4_numeric.py` 能与 Phase 3 local QP 数值吻合；
4. `blom_catmull_compare.py` 明确展示 `s=2` 与 `s=4` 的本质差异；
5. 三个脚本都保存了必要的结果文件和图片；
6. 输出格式统一，可直接服务后续论文表格、图像与 BLOM-Analytic 实现。

---

## 9. 给实现 AI 的最终指令

你要实现的不是“一个局部规划器”，而是 **Phase 4 理论的验证器与解析内核原型**。

请优先保证：

1. 数学推导对象与 Phase 4 严格一致；
2. 所有矩阵均由程序自动生成，可检查、可重现；
3. 与 Phase 1 / 2 / 3 的接口统一；
4. 必须保存足够的中间结果与图片，以便论文中直接展示；
5. 不要跳过 symbolic / numeric / compare 三层验证中的任何一层。

如果某些启发式补全（例如 `s=4` 的 heuristic jerk/acceleration）存在多种可能方案，代码必须：

- 明确标注当前方案只是 baseline；
- 不得把 heuristic baseline 写成 exact theory。
