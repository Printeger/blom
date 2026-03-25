# REQ-Phase10-full-space-time-framework.md

## 1. 文档目的

本需求文档用于指导实现 **Phase 10：基于 raw Scheme C 的完整 block-banded 时空优化框架** 的统一代码与数值验证系统。

本阶段不再回答：

- BLOM 是否接近 global MINCO（这属于 Phase 7 / Phase 8）
- raw Scheme C 是否能反传（这属于 Phase 9）

本阶段要完成的是：

> 把 raw Scheme C 从“最小可微优化闭环”扩展成一个完整、模块化、可 benchmark、复杂度清晰的时空联合优化框架。

本需求文档统一覆盖以下三个核心文件：

- `blom_full_backward_diff.py`
- `blom_space_time_optimizer.py`
- `blom_benchmark_suite.py`

同时要求结果能直接服务于：

1. Phase 10 理论中的 full objective layer；
2. general-\(k\) block-banded backward differentiation；
3. framework-level complexity evaluation；
4. optimizer-level benchmark；
5. ablation study 与 baseline comparison。

---

## 2. 与前面 Phases 完全统一的符号要求

实现必须与 Phase 0--10 的符号保持完全一致。

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
- 外层参数（原变量）：
  \[
  \theta=(\bar q,T)
  \]
- 外层参数（重参数化变量）：
  \[
  \xi=(\bar q,\tau)
  \]

### 2.2 时间重参数化

必须统一使用：

\[
T_i(\tau_i)=T_{\min}+\operatorname{softplus}(\tau_i),
\qquad
\operatorname{softplus}(x)=\log(1+e^x)
\]

并且：

\[
\frac{dT_i}{d\tau_i}=\sigma(\tau_i)=\frac{1}{1+e^{-\tau_i}}
\]

### 2.3 系数对象

raw Scheme C 的 general-\(k\) 系数记为：

\[
c^{\mathrm{C},k}(\bar q,T)\in\mathbb R^{2sM}
\]

并按 block 写成：

\[
c^{\mathrm{C},k}=
\begin{bmatrix}
c_1\\
\vdots\\
c_M
\end{bmatrix},
\qquad
c_i\in\mathbb R^{2s}
\]

在 canonical \(s=4\) 下：
\[
c_i\in\mathbb R^8.
\]

### 2.4 Jacobian 记号

必须统一使用：

\[
J_{c\bar q}^{(k)}(\theta)
:=
\frac{\partial c^{\mathrm{C},k}}{\partial \bar q}(\theta)
\in\mathbb R^{2sM\times(M-1)},
\]
\[
J_{cT}^{(k)}(\theta)
:=
\frac{\partial c^{\mathrm{C},k}}{\partial T}(\theta)
\in\mathbb R^{2sM\times M}.
\]

block 形式记为：

\[
\frac{\partial c_i}{\partial q_j}\in\mathbb R^{2s\times 1},
\qquad
\frac{\partial c_i}{\partial T_j}\in\mathbb R^{2s\times 1}.
\]

---

## 3. 本阶段的理论对象与工程范围

## 3.1 本阶段要实现的总目标函数

必须统一围绕 Phase 10 理论章节中的 full objective layer：

\[
K(c,T)
=
K_{\mathrm{ctrl}}(c,T)
+\lambda_T K_{\mathrm{time}}(T)
+\lambda_{\mathrm{obs}} K_{\mathrm{obs}}(c,T)
+\lambda_{\mathrm{dyn}} K_{\mathrm{dyn}}(c,T)
+\lambda_{\mathrm{bc}} K_{\mathrm{bc}}(c,T)
+\lambda_{\mathrm{reg}} K_{\mathrm{reg}}(T)
\]

定义总目标：

\[
J(\bar q,T)=K(c^{\mathrm{C},k}(\bar q,T),T),
\qquad
\widetilde J(\bar q,\tau)=J(\bar q,T(\tau))
\]

## 3.2 本阶段允许的范围

### 必须支持
1. raw Scheme C
2. general even \(k\)
3. dense checker 与 sparse backward 双轨并存
4. full objective layer
5. optimizer
6. benchmark / ablation

### 不要求本阶段完成的
1. matching theorem
2. feasibility theorem
3. cost theorem 回到 Phase 1 exact feasible family
4. global nonconvex convergence theorem
5. full multidimensional planner deployment

---

## 4. `blom_full_backward_diff.py` 的需求

## 4.1 文件目标

这个文件是 Phase 10 的底层核心，负责实现：

1. general-\(k\) raw Scheme C 系数接口
2. full objective layer 的各项值与梯度
3. dense backward differentiation
4. block-banded backward differentiation
5. 在重参数化变量 \(\xi=(\bar q,\tau)\) 下的总梯度

---

## 4.2 必须实现的核心函数

### 4.2.1 general-\(k\) raw Scheme C 系数接口

```python
def compute_raw_schemeC_coeffs_general_k(q, T, s=4, k=2, config=None):
    ...
```

要求：

- 支持 general even \(k\)
- 返回 block coeffs 与 flattened coeff vector
- 与前面 Phase 4 / Phase 6 / Phase 9 的 raw Scheme C 接口兼容
- 若 general \(k\) 尚未实现，应在代码中显式抛出 `NotImplementedError`，但接口形式必须预留好

输出至少包括：

```python
{
    "coeff_blocks": np.ndarray,   # shape (M, 2*s)
    "coeff_vec": np.ndarray,      # shape (2*s*M,)
}
```

---

### 4.2.2 general-\(k\) Jacobian 接口

```python
def compute_raw_schemeC_jacobians_general_k(q, T, s=4, k=2, mode="analytic", config=None):
    ...
```

必须输出：

```python
{
    "J_c_q_dense": np.ndarray,      # shape (2*s*M, M-1)
    "J_c_T_dense": np.ndarray,      # shape (2*s*M, M)
    "J_c_q_blocks": dict,           # (i,j) -> np.ndarray
    "J_c_T_blocks": dict,           # (i,j) -> np.ndarray
    "support_q": dict,              # i -> list of j
    "support_T": dict,              # i -> list of j
}
```

要求：

- dense 版本用于 checker
- block 版本用于 sparse backward
- 必须保留 support sets 信息，供 benchmark 与 complexity 验证
- 若目前只实现 \(k=2\)，也必须把 general-\(k\) 接口预留完整

---

### 4.2.3 时间重参数化接口

```python
def tau_to_T(tau, T_min=1e-3):
    ...
```

```python
def dT_dtau(tau):
    ...
```

要求：

- 使用 softplus
- 显式返回 \(\sigma(\tau)\)
- 必须数值稳定，避免 \(\exp(\tau)\) 溢出

---

### 4.2.4 control cost

```python
def control_cost_full(coeff_blocks, T, s=4):
    ...
```

输出：

```python
{
    "value": float,
    "grad_c_blocks": np.ndarray,
    "grad_T": np.ndarray,
}
```

必须与理论公式一致：

\[
\frac{\partial K_{\mathrm{ctrl}}}{\partial c_i}=H_i(T_i)c_i,
\qquad
\frac{\partial K_{\mathrm{ctrl}}}{\partial T_i}
=
\frac12 c_i^\top H_i'(T_i)c_i.
\]

---

### 4.2.5 time penalty

```python
def time_penalty_full(T, psi_type="linear", psi_config=None):
    ...
```

至少支持：

- linear：
  \[
  \psi_T(T)=T
  \]
- quadratic：
  \[
  \psi_T(T)=T^2
  \]
- custom callable

输出：

```python
{
    "value": float,
    "grad_T": np.ndarray,
}
```

---

### 4.2.6 obstacle penalty

```python
def obstacle_penalty_full(coeff_blocks, T, obs_config, s=4):
    ...
```

必须支持 sampled penalty：

\[
K_{\mathrm{obs}}(c,T)
=
\sum_{i,\nu}
w_\nu^{\mathrm{obs}}
\Phi_{\mathrm{obs},\nu}(z_{i,\nu}^{\mathcal R_{\mathrm{obs}}},T_i)
\]

输出：

```python
{
    "value": float,
    "grad_c_blocks": np.ndarray,
    "grad_T": np.ndarray,
}
```

要求：

- sample points \(\alpha_\nu^{\mathrm{obs}}\)
- derivative index set \(\mathcal R_{\mathrm{obs}}\)
- \(\Phi_{\mathrm{obs},\nu}\) 可调用且返回 value/grad

---

### 4.2.7 dynamics penalty

```python
def dynamics_penalty_full(coeff_blocks, T, dyn_config, s=4):
    ...
```

形式与 obstacle penalty 类似。

输出：

```python
{
    "value": float,
    "grad_c_blocks": np.ndarray,
    "grad_T": np.ndarray,
}
```

要求：

- 支持 sampled velocity / acceleration / jerk penalty
- \(\mathcal R_{\mathrm{dyn}}\) 必须可配置

---

### 4.2.8 boundary penalty

```python
def boundary_penalty_full(coeff_blocks, T, bc_config, s=4):
    ...
```

形式应对应：

\[
K_{\mathrm{bc}}(c,T)
=
\Phi_{\mathrm{bc}}^-(\cdot)
+
\Phi_{\mathrm{bc}}^+(\cdot,T_M)
\]

输出：

```python
{
    "value": float,
    "grad_c_blocks": np.ndarray,
    "grad_T": np.ndarray,
}
```

---

### 4.2.9 regularization penalty

```python
def regularization_penalty_full(T, reg_config=None):
    ...
```

至少支持 time smoothing：

\[
K_{\mathrm{reg}}(T)=\frac12\sum_{i=1}^{M-1}(T_{i+1}-T_i)^2
\]

输出：

```python
{
    "value": float,
    "grad_T": np.ndarray,
}
```

---

### 4.2.10 full objective layer

```python
def evaluate_full_objective(q, T, weights, obs_config=None, dyn_config=None,
                            bc_config=None, reg_config=None, s=4, k=2, config=None):
    ...
```

必须输出：

```python
{
    "value": float,
    "coeff_blocks": np.ndarray,
    "coeff_vec": np.ndarray,
    "g_c_blocks": np.ndarray,
    "g_c_vec": np.ndarray,
    "g_T": np.ndarray,
    "parts": {
        "ctrl": float,
        "time": float,
        "obs": float,
        "dyn": float,
        "bc": float,
        "reg": float,
    }
}
```

---

### 4.2.11 dense backward differentiation

```python
def full_backward_diff_dense(q, T, weights, obs_config=None, dyn_config=None,
                             bc_config=None, reg_config=None, s=4, k=2, config=None):
    ...
```

必须使用：

\[
\nabla_{\bar q}J
=
\bigl(J_{c\bar q}^{(k)}\bigr)^\top g_c,
\qquad
\nabla_T J
=
\bigl(J_{cT}^{(k)}\bigr)^\top g_c+g_T.
\]

输出：

```python
{
    "grad_q": np.ndarray,
    "grad_T": np.ndarray,
    "grad_theta": np.ndarray,
}
```

---

### 4.2.12 sparse backward differentiation

```python
def full_backward_diff_sparse(q, T, weights, obs_config=None, dyn_config=None,
                              bc_config=None, reg_config=None, s=4, k=2, config=None):
    ...
```

必须使用 support sets 做 block-banded accumulation：

\[
\frac{\partial J}{\partial q_j}
=
\sum_{i:\,j\in\mathcal S_q(i,k)}
\left(\frac{\partial c_i}{\partial q_j}\right)^\top g_{c,i},
\]
\[
\frac{\partial J}{\partial T_j}
=
g_{T,j}
+
\sum_{i:\,j\in\mathcal S_T(i,k)}
\left(\frac{\partial c_i}{\partial T_j}\right)^\top g_{c,i}.
\]

---

### 4.2.13 重参数化变量 \(\xi=(\bar q,\tau)\) 下的梯度

```python
def full_backward_diff_reparam(q, tau, weights, T_min=1e-3,
                               obs_config=None, dyn_config=None,
                               bc_config=None, reg_config=None,
                               s=4, k=2, config=None):
    ...
```

必须实现：

\[
\nabla_{\bar q}\widetilde J=\nabla_{\bar q}J,
\qquad
\nabla_\tau \widetilde J=D_\sigma(\tau)\nabla_T J.
\]

输出：

```python
{
    "T": np.ndarray,
    "grad_q": np.ndarray,
    "grad_tau": np.ndarray,
    "grad_xi": np.ndarray,
}
```

---

## 4.3 必须保存的图片

若 `save_dir` 不为空，必须保存：

### 图 1：dense vs sparse gradient comparison（必须）

文件名建议：

```text
phase10_dense_vs_sparse_grad.png
```

### 图 2：general-\(k\) Jacobian sparsity（必须）

文件名建议：

```text
phase10_jacobian_sparsity_k_sweep.png
```

### 图 3：objective block gradient heatmap（推荐）

文件名建议：

```text
phase10_block_grad_heatmap.png
```

### 图 4：support width vs \(k\)（推荐）

文件名建议：

```text
phase10_support_width_vs_k.png
```

---

## 5. `blom_space_time_optimizer.py` 的需求

## 5.1 文件目标

实现基于重参数化变量
\[
\xi=(\bar q,\tau)
\]
的完整时空优化器。

本阶段不要求 full planner 级部署，但必须支持：

1. projected-free unconstrained optimization（通过 \(\tau\) 重参数化）
2. line search / Armijo
3. 多目标项联合优化
4. objective / gradient / timing logging
5. ablation-ready 接口

---

## 5.2 必须实现的优化器接口

### 5.2.1 单步更新接口

```python
def optimizer_step(xi, weights, T_min=1e-3,
                   obs_config=None, dyn_config=None, bc_config=None, reg_config=None,
                   method="gd", step_rule="armijo", step_size=1e-3,
                   s=4, k=2, config=None):
    ...
```

必须支持：

- `method="gd"`：梯度下降
- `step_rule="fixed"`：固定步长
- `step_rule="armijo"`：回溯线搜索

输出：

```python
{
    "xi_new": np.ndarray,
    "q_new": np.ndarray,
    "tau_new": np.ndarray,
    "T_new": np.ndarray,
    "obj_old": float,
    "obj_new": float,
    "grad_norm": float,
    "step_size_used": float,
    "accepted": bool,
}
```

---

### 5.2.2 多步优化主循环

```python
def run_space_time_optimization(
    q0, tau0, weights,
    obs_config=None, dyn_config=None, bc_config=None, reg_config=None,
    method="gd", step_rule="armijo",
    n_steps=100, step_size=1e-3, T_min=1e-3,
    s=4, k=2, config=None, save_dir=None, seed=42
):
    ...
```

必须记录：

- total objective history
- 各子项 history
- gradient norm history
- step size history
- support width / sparsity stats（可选但推荐）
- wall-clock time per iteration

---

### 5.2.3 支持多种 \(k\)

必须允许：

```python
k_values = [2, 4, 6, ...]
```

用于后续 benchmark 和 ablation。

如果目前 general-\(k\) 尚未完整实现，也必须允许 benchmark framework 按 `try/except + explicit not-implemented` 的方式组织。

---

## 5.3 必须保存的图片

### 图 1：objective curve（必须）

文件名建议：

```text
phase10_objective_curve.png
```

内容：
- total objective
- ctrl / time / obs / dyn / bc / reg 各子项

### 图 2：gradient norm curve（必须）

文件名建议：

```text
phase10_gradnorm_curve.png
```

### 图 3：step size curve（推荐）

文件名建议：

```text
phase10_stepsize_curve.png
```

### 图 4：trajectory before / after / during（必须）

文件名建议：

```text
phase10_traj_progress.png
```

### 图 5：time allocation before / after（推荐）

文件名建议：

```text
phase10_time_progress.png
```

### 图 6：objective decrease per step（推荐）

文件名建议：

```text
phase10_obj_delta_per_step.png
```

---

## 5.4 必须保存的表和摘要

- `phase10_opt_history.csv`
- `phase10_opt_summary.json`
- `phase10_opt_summary.md`

`phase10_opt_summary.md` 至少自动回答：

1. 总目标是否稳定下降；
2. 哪些项是主要下降来源；
3. Armijo / fixed step 是否稳定；
4. \(\tau\)-parameterization 是否避免了 \(T_i\le T_{\min}\) 问题；
5. 当前 optimizer 是否足以进入 benchmark 阶段。

---

## 6. `blom_benchmark_suite.py` 的需求

## 6.1 文件目标

实现完整 benchmark / ablation framework，比较：

1. raw Scheme C optimizer
2. global MINCO baseline
3. Scheme A baseline
4. heuristic baseline

---

## 6.2 benchmark 维度要求

必须至少支持以下 sweep：

### 6.2.1 \(M\) sweep
例如：
```python
M_values = [10, 20, 40, 80]
```

### 6.2.2 \(k\) sweep
例如：
```python
k_values = [2, 4, 6, 8]
```

### 6.2.3 obstacle density / penalty strength sweep
例如：
```python
obs_weight_values = [0.0, 0.1, 1.0, 10.0]
```

### 6.2.4 time-range spread sweep
例如 bounded-nonuniform time boxes：
```python
[(0.9, 1.1), (0.5, 2.0), (0.2, 3.0)]
```

### 6.2.5 initialization quality sweep
例如：
- warm start
- random start
- perturbed reference start

---

## 6.3 baseline 接口要求

至少预留以下统一接口：

```python
def run_baseline_minco(...): ...
def run_baseline_schemeA(...): ...
def run_baseline_heuristic(...): ...
def run_baseline_raw_schemeC(...): ...
```

输出必须统一，便于比较：

```python
{
    "method": str,
    "final_objective": float,
    "runtime_total": float,
    "runtime_per_iter": float,
    "memory_peak": float | None,
    "constraint_violation": float | None,
    "history": dict,
}
```

---

## 6.4 必须计算的 benchmark 指标

至少包括：

### 6.4.1 目标差异

\[
\Delta_{\mathrm{obj}}
=
\frac{K_{\mathrm{final}}-K_{\mathrm{ref}}}{|K_{\mathrm{ref}}|+\varepsilon}
\]

### 6.4.2 单次迭代时间

\[
\Delta_{\mathrm{time}}
=
\text{wall-clock time per iteration}
\]

### 6.4.3 峰值内存

\[
\Delta_{\mathrm{mem}}
=
\text{peak memory}
\]

### 6.4.4 约束 violation（若定义）

\[
\Delta_{\mathrm{viol}}
=
\text{sampled constraint violation}
\]

---

## 6.5 ablation study 要求

必须至少支持：

### 6.5.1 \(k\)-sweep ablation
看随着 \(k\) 增大：

- quality
- runtime
- memory
- gradient cost

如何变化。

### 6.5.2 term ablation
分别去掉：

- \(K_{\mathrm{obs}}\)
- \(K_{\mathrm{dyn}}\)
- \(K_{\mathrm{time}}\)
- \(K_{\mathrm{reg}}\)

看目标、轨迹与 runtime 如何变化。

### 6.5.3 baseline ablation
raw Scheme C vs MINCO vs Scheme A vs heuristic。

---

## 6.6 必须保存的图片

### 图 1：runtime vs \(M\)（必须）

文件名建议：

```text
phase10_runtime_vs_M.png
```

### 图 2：objective gap vs \(k\)（必须）

文件名建议：

```text
phase10_obj_gap_vs_k.png
```

### 图 3：memory vs \(k\)（推荐）

文件名建议：

```text
phase10_memory_vs_k.png
```

### 图 4：baseline comparison bar chart（必须）

文件名建议：

```text
phase10_baseline_compare.png
```

### 图 5：ablation comparison heatmap（必须）

文件名建议：

```text
phase10_ablation_heatmap.png
```

### 图 6：quality-vs-speed Pareto 图（推荐）

文件名建议：

```text
phase10_quality_speed_pareto.png
```

### 图 7：support width / actual sparsity vs \(k\)（推荐）

文件名建议：

```text
phase10_sparsity_vs_k.png
```

---

## 6.7 必须保存的表和摘要

- `phase10_benchmark_summary.csv`
- `phase10_ablation_summary.csv`
- `phase10_benchmark_summary.json`
- `phase10_benchmark_summary.md`

`phase10_benchmark_summary.md` 至少自动回答：

1. raw Scheme C 相比 baseline 的速度 / 质量表现如何；
2. 是否观察到随 \(k\) 增大带来的质量-速度 trade-off；
3. 哪些 penalty term 对最终表现影响最大；
4. Phase 10 是否已经足以支撑“完整 block-banded framework”的经验结论。

---

## 7. 推荐的统一总控脚本接口

建议额外提供一个 Phase 10 总控函数：

```python
def run_phase10_framework_suite(
    q0=None,
    tau0=None,
    weights=None,
    obs_config=None,
    dyn_config=None,
    bc_config=None,
    reg_config=None,
    M_values=None,
    k_values=None,
    benchmark_methods=None,
    n_steps=100,
    T_min=1e-3,
    save_dir=None,
    seed=42
):
    ...
```

必须完成：

1. dense vs sparse backward consistency check
2. optimizer demo
3. \(M\)-sweep / \(k\)-sweep
4. baseline compare
5. ablation suite
6. 自动生成总 summary

---

## 8. 统一结果目录要求

所有结果必须保存到：

```text
results/phase10_framework/
```

建议目录结构：

```text
results/phase10_framework/
├── backward_diff/
├── optimizer/
├── benchmark_M_sweep/
├── benchmark_k_sweep/
├── baseline_compare/
├── ablation/
└── summary/
```

---

## 9. 必须保存的总表和总 summary

### 9.1 总对比表（必须）

文件名建议：

```text
phase10_overview.csv
```

至少包含：

- method
- M
- k
- final_objective
- objective_gap
- runtime_total
- runtime_per_iter
- memory_peak
- constraint_violation
- support_width_mean
- support_width_max

### 9.2 总结性 summary（必须）

文件名建议：

```text
phase10_interpretation_summary.md
```

必须自动回答：

1. full objective layer 是否稳定工作；
2. general-\(k\) block-banded backward 是否数值一致；
3. runtime / memory 是否符合 \(O(k^2M)\) 与 \(O(kM)\) 的经验趋势；
4. raw Scheme C optimizer 是否已经显著优于 heuristic；
5. 与 MINCO / Scheme A 相比，是否观察到明确的 speed-quality trade-off；
6. 现有结果是否足以支撑论文中“完整 block-banded space-time optimization framework”的表述。

---

## 10. 推荐代码结构

```text
phase_10/
├── blom_full_backward_diff.py
├── blom_space_time_optimizer.py
├── blom_benchmark_suite.py
├── blom_phase10_framework_suite.py
├── test_blom_full_backward_diff.py
├── test_blom_space_time_optimizer.py
├── test_blom_benchmark_suite.py
├── examples/
│   ├── demo_full_backward_diff.py
│   ├── demo_optimizer.py
│   ├── demo_k_sweep.py
│   ├── demo_baseline_compare.py
│   ├── demo_ablation.py
│   └── demo_phase10_suite.py
└── results/
    └── phase10_framework/
```

---

## 11. 验收标准

只有满足以下条件，本需求对应的实现才算完成：

1. 与前面 phases 的符号系统完全统一；
2. full objective layer 各项值与梯度都可计算；
3. dense backward 与 sparse backward 数值一致；
4. optimizer 可在 \(\xi=(\bar q,\tau)\) 空间稳定运行；
5. benchmark / baseline / ablation 三大模块都能独立运行；
6. 必须保存所有要求的图、表、JSON、summary；
7. 结果足以直接支撑 Phase 10 论文中的核心结论。

---

## 12. 给实现 AI 的最终指令

你要实现的是一套 **Phase 10 完整 block-banded 时空优化框架**，不是最小 demo。

请优先保证：

1. full objective layer 完整；
2. dense checker 与 sparse implementation 双轨并存；
3. \(\tau\)-重参数化稳定、统一；
4. optimizer 与 benchmark 解耦，接口清晰；
5. baseline compare 与 ablation 必须可复用；
6. 图片和 summary 适合直接进入论文、补充材料或组会汇报；
7. 不要把还没证明的 theorem 强行写进代码注释或 summary。

本阶段的核心任务是：**把 raw Scheme C 从“最小可微闭环”推进成“完整的 block-banded 时空联合优化框架”，并用系统 benchmark 证明它在工程与算法层面是一个严肃可比较的对象。**
