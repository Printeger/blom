# REQ-Phase9-minimal-differentiable-loop.md

## 1. 文档目的

本需求文档用于指导实现 **Phase 9：raw Scheme C 的最小可微优化闭环** 的统一数值验证与代码框架。

本阶段的目标不是直接实现完整 planner，而是严格围绕以下最小闭环：

\[
\theta=(\bar q,T)
\ \longmapsto\
c(\theta)=c^{\mathrm{C},2}(\bar q,T)
\ \longmapsto\
K(c,T)
\ \longmapsto\
\nabla_\theta J(\theta),
\qquad
J(\theta)=K(c(\theta),T),
\]

验证以下问题：

1. raw Scheme C 的系数映射 \(c(\theta)\) 是否能稳定支持反向微分；
2. 解析梯度
   \[
   \nabla_\theta J
   \]
   是否与 finite difference 一致；
3. banded backward accumulation 是否与 dense checker 完全一致；
4. 最小优化 demo 中，目标函数是否能稳定下降；
5. 这一最小优化闭环是否足以作为进入下一阶段完整优化框架的数学与工程基础。

本需求文档统一覆盖以下三个文件：

- `blom_backward_diff.py`
- `blom_space_time_opt_demo.py`
- `test_blom_backward_diff.py`

---

## 2. 与前面 Phases 完全统一的符号要求

实现必须与 Phase 0--9 的符号保持完全一致。

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
  \bar q=(q_1,\dots,q_{M-1})^\top\in\mathbb R^{M-1}
  \]
- 时间向量：
  \[
  T=(T_1,\dots,T_M)^\top,\qquad T_i>0
  \]
- 总时长：
  \[
  \mathcal T=\sum_{i=1}^{M}T_i
  \]
- 外层参数：
  \[
  \theta=(\bar q,T)
  \]

### 2.2 canonical setting

本阶段只允许 canonical minimal loop：

\[
s=4,\qquad k=2
\]

且只允许使用：

- raw Scheme C
- 不允许使用 Scheme A / B
- 不允许使用 C1 / C2 修正版

### 2.3 系数对象

raw Scheme C 的系数向量记为：

\[
c(\theta)=c^{\mathrm{C},2}(\bar q,T)\in\mathbb R^{2sM}=\mathbb R^{8M}
\]

并按段 block 存储为：

\[
c(\theta)=
\begin{bmatrix}
c_1(\theta)\\
\vdots\\
c_M(\theta)
\end{bmatrix},
\qquad
c_i(\theta)\in\mathbb R^{2s}=\mathbb R^8.
\]

### 2.4 Jacobian 记号

必须统一使用：

\[
J_{c\bar q}(\theta):=\frac{\partial c}{\partial \bar q}(\theta)\in\mathbb R^{2sM\times(M-1)},
\]
\[
J_{cT}(\theta):=\frac{\partial c}{\partial T}(\theta)\in\mathbb R^{2sM\times M}.
\]

block 形式记为：

\[
\frac{\partial c_i}{\partial q_j}\in\mathbb R^{2s\times 1},
\qquad
\frac{\partial c_i}{\partial T_j}\in\mathbb R^{2s\times 1}.
\]

---

## 3. 本阶段的理论对象与工程范围

## 3.1 只做最小目标函数

本阶段只考虑：

\[
J(\theta):=
K(c(\theta),T)
=
K_{\mathrm{ctrl}}(c,T)
+
\lambda_T K_{\mathrm{time}}(T)
+
\lambda_{\mathrm{obs}} K_{\mathrm{soft\mbox{-}obs}}(c,T).
\]

### 3.1.1 控制代价

\[
K_{\mathrm{ctrl}}(c,T)
:=
\frac12\sum_{i=1}^{M} c_i^\top H_i(T_i)c_i.
\]

### 3.1.2 时间正则项

建议默认：

\[
K_{\mathrm{time}}(T)=\sum_{i=1}^{M}T_i.
\]

### 3.1.3 软障碍项

必须采用最小、光滑、可导的 soft penalty，不引入复杂 ESDF 系统。

建议采用：

\[
K_{\mathrm{soft\mbox{-}obs}}(c,T)
=
\sum_{i=1}^{M}\sum_{\nu=1}^{N_s}
w_\nu\,\varphi_{\mathrm{obs}}(x_{i,\nu}(c_i,T_i)),
\]
其中
\[
x_{i,\nu}(c_i,T_i)=B^{(0)}(\alpha_\nu T_i)c_i.
\]

默认要求：

- \(\varphi_{\mathrm{obs}}\in C^1\)
- \(N_s\) 固定
- \(\alpha_\nu\in[0,1]\)

---

## 3.2 本阶段不做的内容

以下内容必须明确排除：

1. 不做完整多约束 planner
2. 不做硬动力学约束
3. 不做复杂障碍系统 / ESDF-heavy system
4. 不做完整 general-\(k\) 优化框架
5. 不做 Scheme A/B/C 三套并行优化
6. 不证明 full \(O(k^2M)\) complexity theorem
7. 不证明全局 nonconvex optimization convergence theorem

---

## 4. `blom_backward_diff.py` 的需求

## 4.1 文件目标

实现 raw Scheme C 下最小可微闭环的核心微分接口：

\[
(\bar q,T)\mapsto c(\theta),\qquad
(c,T)\mapsto K(c,T),\qquad
(\bar q,T)\mapsto \nabla_\theta J(\theta).
\]

这个文件是本阶段最核心的底层文件。

---

## 4.2 必须实现的函数

### 4.2.1 raw Scheme C 系数接口

```python
def compute_raw_schemeC_coeffs(q, T, s=4, k=2, config=None):
    ...
```

要求：

- 返回 block coeffs 与 flattened coeff vector
- 只支持 canonical \(s=4,k=2\)
- 必须与前面 Phase 4 / Phase 6 的 raw Scheme C 接口兼容

输出至少包括：

```python
{
    "coeff_blocks": np.ndarray,   # shape (M, 2*s)
    "coeff_vec": np.ndarray,      # shape (2*s*M,)
}
```

---

### 4.2.2 raw Scheme C Jacobian 接口

```python
def compute_raw_schemeC_jacobians(q, T, s=4, k=2, mode="analytic", config=None):
    ...
```

必须输出：

```python
{
    "J_c_q_dense": np.ndarray,      # shape (2*s*M, M-1)
    "J_c_T_dense": np.ndarray,      # shape (2*s*M, M)
    "J_c_q_blocks": dict,           # (i,j) -> block
    "J_c_T_blocks": dict,           # (i,j) -> block
}
```

要求：

- dense 版本用于 checker
- block 版本用于真正 banded backward
- block 稀疏模式必须和 Phase 6 一致

---

### 4.2.3 control cost

```python
def control_cost(coeff_blocks, T, s=4):
    ...
```

必须同时返回：

```python
{
    "value": float,
    "grad_c_blocks": np.ndarray,   # g_{c,i}^{ctrl}, shape (M, 2*s)
    "grad_T": np.ndarray,          # shape (M,)
}
```

要求与公式一致：

\[
\frac{\partial K_{\mathrm{ctrl}}}{\partial c_i}=H_i(T_i)c_i,
\qquad
\frac{\partial K_{\mathrm{ctrl}}}{\partial T_i}
=
\frac12 c_i^\top H_i'(T_i)c_i.
\]

---

### 4.2.4 time penalty

```python
def time_penalty(T, weight=1.0):
    ...
```

默认实现：

\[
K_{\mathrm{time}}(T)=\sum_i T_i
\]

输出：

```python
{
    "value": float,
    "grad_T": np.ndarray,
}
```

---

### 4.2.5 soft obstacle penalty

```python
def soft_obstacle_penalty(coeff_blocks, T, obs_config, s=4):
    ...
```

必须输出：

```python
{
    "value": float,
    "grad_c_blocks": np.ndarray,   # shape (M, 2*s)
    "grad_T": np.ndarray,          # shape (M,)
}
```

要求：

- 用固定 sample points \(\alpha_\nu\)
- \(\varphi_{\mathrm{obs}}\) 和 \(\varphi'_{\mathrm{obs}}\) 必须显式可调用
- 输出必须与链式法则一致

---

### 4.2.6 总目标函数与 block 梯度

```python
def evaluate_minimal_objective(q, T, weights, obs_config=None, s=4, k=2, config=None):
    ...
```

必须输出：

```python
{
    "value": float,
    "coeff_blocks": np.ndarray,
    "coeff_vec": np.ndarray,
    "g_c_blocks": np.ndarray,      # total ∂K/∂c_i
    "g_c_vec": np.ndarray,         # flattened
    "g_T": np.ndarray,             # total ∂K/∂T
    "parts": {
        "ctrl": float,
        "time": float,
        "obs": float,
    }
}
```

---

### 4.2.7 dense backward gradient

```python
def backward_diff_dense(q, T, weights, obs_config=None, s=4, k=2, config=None):
    ...
```

必须使用：

\[
\nabla_{\bar q}J = J_{c\bar q}^\top g_c,
\qquad
\nabla_T J = J_{cT}^\top g_c + g_T.
\]

输出：

```python
{
    "grad_q": np.ndarray,   # shape (M-1,)
    "grad_T": np.ndarray,   # shape (M,)
    "grad_theta": np.ndarray,
}
```

---

### 4.2.8 banded backward gradient

```python
def backward_diff_banded(q, T, weights, obs_config=None, s=4, k=2, config=None):
    ...
```

必须使用 raw Scheme C 的 banded accumulation：

\[
\frac{\partial J}{\partial q_j}
=
\sum_{i=\max(1,j-1)}^{\min(M,j+2)}
\left(\frac{\partial c_i}{\partial q_j}\right)^\top g_{c,i},
\]
\[
\frac{\partial J}{\partial T_j}
=
g_{T,j}
+
\sum_{i=\max(1,j-1)}^{\min(M,j+1)}
\left(\frac{\partial c_i}{\partial T_j}\right)^\top g_{c,i}.
\]

输出格式与 dense 版本保持一致。

---

## 4.3 必须保存的图片

若 `save_dir` 不为空，必须保存：

### 图 1：dense vs banded gradient 对比图（必须）

文件名建议：

```text
phase9_dense_vs_banded_grad.png
```

内容：
- \(\nabla_{\bar q}J\) dense vs banded
- \(\nabla_TJ\) dense vs banded

### 图 2：block gradient heatmap（推荐）

文件名建议：

```text
phase9_block_grad_heatmap.png
```

内容：
- \(g_{c,i}\) 的 block 热图
- 可显示哪些段贡献最大

### 图 3：Jacobian sparsity visualization（推荐）

文件名建议：

```text
phase9_jacobian_sparsity.png
```

内容：
- `J_c_q_dense`
- `J_c_T_dense`

---

## 5. `test_blom_backward_diff.py` 的需求

## 5.1 文件目标

这个文件负责数值验证，不做业务逻辑实现。必须验证：

1. analytic gradient vs finite difference
2. dense backward vs banded backward
3. 各目标项单独梯度正确
4. 对 \(q\) 与 \(T\) 的梯度都正确

---

## 5.2 必须实现的测试函数

### 5.2.1 finite difference gradient checker

```python
def finite_difference_gradient(f, x, eps=1e-6, method="central"):
    ...
```

要求：

- central difference
- 支持向量输入
- 支持输出标量目标函数

---

### 5.2.2 对总目标的梯度检查

```python
def test_total_gradient_matches_fd():
    ...
```

必须比较：

- analytic dense gradient vs FD
- analytic banded gradient vs FD

至少报告误差：

\[
\|\nabla_{\theta}^{\mathrm{analytic}} - \nabla_{\theta}^{\mathrm{FD}}\|_\infty,
\qquad
\|\cdot\|_2
\]

---

### 5.2.3 dense vs banded 一致性测试

```python
def test_dense_equals_banded():
    ...
```

必须验证：

\[
\nabla_\theta^{\mathrm{dense}}J
=
\nabla_\theta^{\mathrm{banded}}J
\]

数值容忍度建议：

```python
atol = 1e-10
rtol = 1e-8
```

---

### 5.2.4 各项梯度单独测试

```python
def test_control_term_gradient():
    ...
def test_time_term_gradient():
    ...
def test_soft_obstacle_gradient():
    ...
```

目的：

- 一旦总梯度错了，方便定位具体是哪一部分出了问题。

---

## 5.3 必须保存的图片

### 图 1：gradient check 误差图（必须）

文件名建议：

```text
phase9_gradcheck_error_q.png
phase9_gradcheck_error_T.png
```

### 图 2：analytic vs FD 散点图（推荐）

文件名建议：

```text
phase9_gradcheck_scatter.png
```

内容：
- analytic gradient
- FD gradient
- 理想对角线

---

## 5.4 必须保存的表和摘要

- `phase9_gradcheck_summary.csv`
- `phase9_gradcheck_summary.json`
- `phase9_gradcheck_summary.md`

`phase9_gradcheck_summary.md` 至少自动回答：

1. 总梯度是否通过 finite difference；
2. dense vs banded 是否完全一致；
3. 哪一项（control / time / obs）误差最大；
4. 现有解析链是否足以进入最小优化 demo。

---

## 6. `blom_space_time_opt_demo.py` 的需求

## 6.1 文件目标

实现一个 **最小优化 demo**，证明 raw Scheme C 不只是“能算梯度”，而是真的能进入一个稳定下降的优化闭环。

这个 demo 不是 full planner，只是 minimal example。

---

## 6.2 优化问题要求

建议默认优化：

\[
\min_{\bar q,T}\ 
K_{\mathrm{ctrl}}(c(\theta),T)
+
\lambda_T \sum_i T_i
+
\lambda_{\mathrm{obs}}K_{\mathrm{soft\mbox{-}obs}}(c(\theta),T)
\]

其中：

- raw Scheme C
- \(s=4,k=2\)
- obstacle 使用最小光滑 penalty
- \(T_i>T_{\min}\) 必须保持可行

---

## 6.3 必须实现的函数

### 6.3.1 单步更新函数

```python
def gradient_descent_step(q, T, weights, obs_config=None,
                          step_size=1e-3, T_min=1e-3,
                          use_banded=True, config=None):
    ...
```

要求：

- 默认沿 \(-\nabla J\) 更新
- 对 \(T\) 的更新后必须投影或裁剪到 \(T_i>T_{\min}\)
- 返回更新前后目标值

输出至少包括：

```python
{
    "q_new": np.ndarray,
    "T_new": np.ndarray,
    "obj_old": float,
    "obj_new": float,
    "grad_norm": float,
}
```

---

### 6.3.2 多步优化 demo

```python
def run_minimal_optimization_demo(q0, T0, weights, obs_config=None,
                                  n_steps=50, step_size=1e-3,
                                  T_min=1e-3, use_banded=True,
                                  save_dir=None, seed=42):
    ...
```

必须记录：

- objective history
- control/time/obs 各子项 history
- gradient norm history
- parameter history（至少可选）

---

## 6.4 必须保存的图片

### 图 1：objective curve（必须）

文件名建议：

```text
phase9_objective_curve.png
```

内容：
- total objective
- control term
- time term
- obs term

### 图 2：gradient norm curve（必须）

文件名建议：

```text
phase9_gradnorm_curve.png
```

### 图 3：trajectory before vs after（必须）

文件名建议：

```text
phase9_traj_before_after.png
```

内容：
- 初始轨迹
- 优化后轨迹
- obstacle region（若配置中有）

### 图 4：time allocation before vs after（推荐）

文件名建议：

```text
phase9_time_before_after.png
```

### 图 5：objective decrease per step（推荐）

文件名建议：

```text
phase9_obj_delta_per_step.png
```

---

## 6.5 必须保存的表和摘要

- `phase9_opt_history.csv`
- `phase9_opt_summary.json`
- `phase9_opt_summary.md`

`phase9_opt_summary.md` 至少自动回答：

1. 目标函数是否稳定下降；
2. 哪个子项下降最明显；
3. 梯度范数是否下降；
4. 是否出现 \(T_i\) 触碰下界；
5. 该最小优化闭环是否支持进入下一阶段完整优化框架。

---

## 7. 推荐的总控脚本接口

建议额外提供一个总控函数，统一跑：

- gradient check
- dense vs banded 对比
- 最小优化 demo

```python
def run_phase9_validation_suite(
    q0,
    T0,
    weights,
    obs_config=None,
    n_steps=50,
    step_size=1e-3,
    T_min=1e-3,
    s=4,
    k=2,
    save_dir=None,
    seed=42
):
    ...
```

输出至少包括：

```python
{
    "gradcheck": {...},
    "dense_vs_banded": {...},
    "optimization_demo": {...},
}
```

---

## 8. 必须保存的总表和总 summary

### 8.1 总对比表（必须）

文件名建议：

```text
phase9_overview.csv
```

至少包含：

- test_name
- passed
- max_abs_error
- l2_error
- dense_banded_gap
- objective_drop
- grad_norm_drop

### 8.2 总结性 summary（必须）

文件名建议：

```text
phase9_interpretation_summary.md
```

必须自动回答：

1. 反向微分链条是否正确；
2. block-banded backward 是否与 dense checker 一致；
3. 最小优化闭环是否稳定；
4. 是否已经足以支撑进入下一阶段：
   - full \(O(k^2M)\) block-banded optimization
   - 更复杂约束与完整优化框架

---

## 9. 推荐代码结构

```text
phase_9/
├── blom_backward_diff.py
├── blom_space_time_opt_demo.py
├── test_blom_backward_diff.py
├── blom_phase9_validation_suite.py
├── examples/
│   ├── demo_backward_diff.py
│   ├── demo_gradcheck.py
│   ├── demo_minimal_opt.py
│   └── demo_phase9_suite.py
└── results/
    └── phase9_validation/
```

---

## 10. 验收标准

只有满足以下条件，本需求对应的实现才算完成：

1. 与前面 phases 的符号系统完全统一；
2. raw Scheme C 的系数、Jacobian、目标项梯度接口完整可用；
3. analytic gradient 与 finite difference 一致；
4. dense backward 与 banded backward 完全一致；
5. 最小优化 demo 中目标函数稳定下降；
6. 必须保存所有要求的图、表、JSON、summary；
7. 结果足以直接支撑 Phase 9 论文中的核心结论。

---

## 11. 给实现 AI 的最终指令

你要实现的是一套 **Phase 9 最小可微优化闭环验证框架**，不是完整 planner。

请优先保证：

1. 数学链条正确；
2. 符号与前面 phases 完全一致；
3. dense checker 与 banded implementation 双轨并存；
4. gradient check 严格；
5. 图片与 summary 足够直观，可直接用于论文或汇报；
6. 不要提前扩展到 Scheme A/B/C 全框架；
7. 不要提前写 full \(O(k^2M)\) theorem 对应的实现。

本阶段的核心任务是：**证明 raw Scheme C 已经是一个真正可微、可反传、可下降的表示对象，从而为下一阶段完整优化框架提供最小而坚实的基础。**
