# REQ-Phase0.md

## 1. 文档目的

本文档定义 BLOM 理论构建 **Phase 0** 的最小代码需求，用于指导 AI 或开发者进行实现。

本阶段目标不是完成完整轨迹优化器，而是建立一个**统一、稳定、可扩展**的理论研究代码骨架，为后续工作提供一致的数据接口与基础求解能力。

本阶段重点支持以下后续任务：

- 全局 MINCO 标量参考问题求解
- BLOM 局部窗口问题数值求解
- BLOM-Analytic 局部线性系统实现
- 局部支撑、插值性、连续性等基础性质验证
- 后续 BLOM 与 MINCO 的对比实验、收敛实验与 Jacobian 验证

---

## 2. Canonical Setting

Phase 0 固定为最小可研究版本，不追求 general case。

### 2.1 固定参数

- 控制代价阶次：`s = 4`
- 局部窗口参数：`k = 2`
- 中间点约束阶数：`d_i = 1`
- 轨迹维数：`dim = 1`
- 边界条件：仅考虑**自然边界条件**

### 2.2 含义说明

- `s = 4`：minimum snap，即代价为 4 阶导数平方积分
- `k = 2`：3-point stencil 的最小局部窗口
- `d_i = 1`：中间点只约束位置，不额外约束速度/加速度
- `dim = 1`：先做标量 1D，后续再扩展到 `R^m`

### 2.3 当前不做的内容

以下内容明确不属于 Phase 0：

- 任意 `s`
- 任意 `k`
- 任意 `d_i`
- 周期边界条件
- 多维 `R^m` 一般情况
- 障碍物约束、碰撞代价、时间正则项
- planner loop
- GTSAM / CasADi / cvxpy / 图优化器

---

## 3. Phase 0 总体目标

Phase 0 需要完成以下四个基础层：

1. **统一输入层**：统一的 `(q, T)` 数据结构和窗口接口
2. **全局参考层**：标量 1D 的全局 MINCO 参考求解器
3. **局部窗口层**：`s=4, k=2` 的 BLOM 局部窗口数值求解器
4. **验证层**：插值、连续性、边界条件残差等基础检查模块

完成后，后续 Phase 1–2 应能直接复用这套输入/输出接口，而不需要推倒重来。

---

## 4. 代码范围

Phase 0 应包含三类文件：

### 4.1 核心模块（必须实现）

- `blom_problem_setup.py`
- `minco_scalar_baseline.py`
- `blom_local_qp.py`

### 4.2 支撑模块（强烈建议同时实现）

- `poly_basis.py`
- `constraint_builders.py`
- `trajectory_eval.py`
- `phase0_checks.py`

### 4.3 Demo / 入口模块

- `demo_phase0.py`

---

## 5. 目录结构要求

建议目录结构如下：

```text
blom_phase0/
├── blom_problem_setup.py
├── minco_scalar_baseline.py
├── blom_local_qp.py
├── poly_basis.py
├── constraint_builders.py
├── trajectory_eval.py
├── phase0_checks.py
├── demo_phase0.py
└── tests/
    ├── test_problem_setup.py
    ├── test_minco_scalar_baseline.py
    ├── test_blom_local_qp.py
    └── test_phase0_checks.py
```

---

## 6. 文件级需求

## 6.1 `blom_problem_setup.py`

### 目标

统一生成 Phase 0 的 canonical 输入数据，作为后续所有模块的标准输入源。

### 必须支持

1. 生成 waypoint 序列 `q_0, ..., q_M`
2. 生成时间向量 `T_1, ..., T_M`
3. 固定 canonical 参数：
   - `s = 4`
   - `k = 2`
   - `d_i = 1`
   - `dim = 1`
4. 提供统一参数向量：
   
   ```math
   \theta = (q_1, \dots, q_{M-1}, T_1, \dots, T_M)
   ```

5. 提供局部窗口函数：

   ```math
   W(i, k) = \{i - \lfloor k/2 \rfloor, \dots, i + \lfloor k/2 \rfloor\} \cap [1, M]
   ```

6. 提供 JSON 序列化 / 反序列化
7. 提供输入合法性检查：
   - `M >= 2`
   - 所有 `T_i > 0`
   - waypoint 维度一致

### 建议数据结构

```python
class BLOMProblemSetup:
    M: int
    q: np.ndarray      # shape (M+1,)
    T: np.ndarray      # shape (M,)
    s: int             # fixed 4
    k: int             # fixed 2
    d_i: int           # fixed 1
    dim: int           # fixed 1
```

### 必须提供的方法

- `make_random(...)`
- `theta()`
- `window(i)`
- `to_dict()`
- `from_dict()`
- `validate()`

### 输出要求

- 返回对象必须可被后续 MINCO 与 BLOM solver 直接消费
- 不允许每个 solver 再自行定义一套 `q, T`

---

## 6.2 `poly_basis.py`

### 目标

统一提供 1D 多项式基函数及其导数评估，避免 MINCO 和 BLOM 两边重复写基础数学代码。

### 必须支持

1. 单段 `2s - 1` 次多项式的系数表示
   - 对 Phase 0：`s = 4`，每段是 **7 次多项式**
2. 在任意 `t ∈ [0, T]` 处评估：
   - `p(t)`
   - `p'(t)`
   - `p''(t)`
   - `p'''(t)`
   - `p''''(t)`
3. 计算最小 snap 代价：

   ```math
   \int_0^T \left(p^{(4)}(t)\right)^2 dt
   ```

4. 提供 Hermite 风格端点导数到多项式系数的转换接口
   - 当前阶段允许先保留函数壳

### 建议函数

- `monomial_row(t, degree)`
- `derivative_row(t, degree, order)`
- `eval_poly(coeffs, t, order=0)`
- `snap_cost_matrix(T, degree=7)`

---

## 6.3 `constraint_builders.py`

### 目标

统一构造约束矩阵，避免在不同 solver 中复制粘贴约束逻辑。

### 必须支持

#### A. 全局 MINCO 约束装配

用于构造标量 1D 全局参考问题的线性约束系统。

应支持：

1. 每段端点插值：

   ```math
   p_i(0) = q_{i-1}, \quad p_i(T_i) = q_i
   ```

2. 相邻段在连接点的 `C^(2s-2)` 连续性
   - 对 Phase 0：`C^6`

3. 起点终点自然边界条件

#### B. BLOM 局部窗口约束装配

用于构造 `s=4, k=2` 局部窗口问题的约束矩阵。

应支持：

1. 窗口内各段端点插值
2. 窗口内 `C^6` 连续性
3. 自然边界条件

自然边界条件可先采用：

```math
p_{\text{left}}^{(l)}(0) = 0, \quad p_{\text{right}}^{(l)}(T) = 0,
\qquad l = s, \dots, 2s-2
```

### 建议函数

- `build_global_minco_constraints(setup)`
- `build_local_blom_constraints(setup, center_segment)`

---

## 6.4 `minco_scalar_baseline.py`

### 目标

实现一个**标量 1D 的全局 MINCO 参考求解器**，用于作为 BLOM 的全局比较基准。

### 作用

1. 统一全局参考母问题的数据接口
2. 为后续不可能性分析提供 ground truth
3. 为 BLOM-vs-MINCO 的数值对比提供参考解

### 必须支持

1. 输入：`BLOMProblemSetup`
2. 输出：
   - 每段多项式系数
   - 总控制代价
   - 全局轨迹评估接口
3. 使用统一的 constraint builder
4. 至少支持一种数值求解方式：
   - 线性方程组求解
   - KKT 求解

### 建议结果对象

```python
class MincoScalarResult:
    coeffs: np.ndarray
    cost: float
    setup: BLOMProblemSetup
```

### 验收要求

- 所有中间点插值误差 `< 1e-8`
- 相邻段在连接点处 0–6 阶导连续误差 `< 1e-8`
- 结果能被 `trajectory_eval.py` 正确评估
- 代价必须为非负实数

### 非目标

- 不要求完整复刻 GCOPTER 的线性复杂度实现
- 不要求多维版本
- 不要求梯度

---

## 6.5 `blom_local_qp.py`

### 目标

只针对 **`s=4, k=2`** 的局部窗口问题实现数值求解。

### 作用

这是 BLOM-Strict 的第一个可运行版本。

它不要求给出解析解，但必须能根据局部 QP 定义稳定求出中心段对应的局部系数。

### 必须支持

1. 输入：`BLOMProblemSetup` 与 `center_segment`
2. 自动构造对应窗口 `W(i, 2)`
3. 构造局部目标：

   ```math
   \sum_{j \in W(i,2)} \int_0^{T_j} \left(p_j^{(4)}(t)\right)^2 dt
   ```

4. 构造局部约束：
   - 端点插值
   - 窗口内 `C^6` 连续性
   - 自然边界条件
5. 返回：
   - 窗口内所有段系数
   - 中心段系数
   - 局部最优代价
   - 活动窗口索引

### 重要说明

本文件只负责**窗口内问题求解**，不负责整条轨迹装配的全局正确性。

当前 BLOM 理论中，窗口边界高阶连续性不是自动成立的，因此不能在本文件里默认“局部解拼起来就全局正确”。

### 验收要求

- 局部约束残差 `< 1e-8`
- 同一输入重复运行结果一致
- 随机输入下可稳定求解
- 可提取中心段系数作为后续 BLOM 装配原子对象

---

## 6.6 `trajectory_eval.py`

### 目标

统一评估全局 MINCO 和局部 BLOM 输出的 piecewise polynomial trajectory。

### 必须支持

1. 在任意 global time 上定位对应 segment
2. 评估：
   - 位置
   - 1–4 阶导数
3. 支持批量采样
4. 计算连接点各阶导 jump：

   ```math
   \llbracket p^{(l)} \rrbracket
   ```

### 建议函数

- `eval_piecewise(coeffs, T, t, order=0)`
- `sample_trajectory(coeffs, T, num=200, order=0)`
- `junction_jumps(coeffs, T, max_order=6)`

---

## 6.7 `phase0_checks.py`

### 目标

集中放置 Phase 0 的基础正确性检查逻辑。

### 必须支持

#### A. 输入检查

- 时间是否全正
- 段数是否合理
- 数据 shape 是否一致

#### B. MINCO 结果检查

- 插值误差
- 全局连续性误差
- 代价非负

#### C. BLOM 局部结果检查

- 局部插值误差
- 窗口内连续性误差
- 自然边界条件残差

#### D. 后续预留

- 局部支撑数值检查
- BLOM-vs-MINCO 差异统计

### 建议函数

- `check_setup(setup)`
- `check_minco_result(result)`
- `check_blom_local_result(local_result)`

---

## 6.8 `demo_phase0.py`

### 目标

提供一个最小端到端 demo，用于快速确认 Phase 0 的骨架已经连通。

### 必须执行的流程

1. 生成随机 `BLOMProblemSetup`
2. 调用 `minco_scalar_baseline.py`
3. 选一个中心段，调用 `blom_local_qp.py`
4. 打印：
   - `q`
   - `T`
   - MINCO 总代价
   - BLOM 局部代价
   - 基础残差检查结果
5. 可选：输出采样值，但当前不要求画图

---

## 7. 统一设计要求

## 7.1 数据结构统一

所有 solver 必须共用同一个 `BLOMProblemSetup`。

不允许每个文件自行定义独立的 `q, T` 接口。

## 7.2 先 1D，不提前做泛化设计

Phase 0 先只做 `dim = 1`。

代码可以预留未来扩展到 `dim > 1` 的空间，但不要一开始就写成复杂泛型框架。

## 7.3 依赖约束

只允许使用：

- `numpy`
- `scipy`
- Python 标准库

不允许引入：

- `casadi`
- `cvxpy`
- `gtsam`
- `jax`
- `torch`

## 7.4 随机种子

所有 demo / 测试统一使用：

```python
rng = np.random.default_rng(42)
```

## 7.5 注释要求

注释必须说明：

- 当前数学对象是什么
- 向量 / 矩阵 shape 是什么
- 当前实现对应哪条理论假设
- 哪些地方只是数值版，不是最终解析版

## 7.6 风格要求

- 使用简单、清晰、可读的实现
- 避免炫技式抽象
- 优先保证接口稳定和数学语义清楚
- 允许后续最小修改替换为 BLOM-Analytic 版

---

## 8. Phase 0 之后建议继续写的代码

以下内容不要求在 Phase 0 立即完成，但建议在 Phase 0 完成后优先接上。

### 8.1 `blom_k2_exact_solver.py`

目标：将 `blom_local_qp.py` 的数值局部求解，进一步替换成 BLOM-Analytic 需要的局部线性系统：

```math
A_k(T)x = B_k(q, T)
```

### 8.2 `blom_boundary_jump_check.py`

目标：显式检查窗口边界高阶导 jump。

### 8.3 `blom_fd_jacobian_check.py`

目标：对后续局部支撑 Proposition 做数值准备，检查：

```math
\frac{\partial c_i}{\partial q_j}, \qquad \frac{\partial c_i}{\partial T_j}
```

是否呈现预期带状稀疏模式。

---

## 9. 验收标准

Phase 0 完成的最低标准如下：

### A. 统一输入层成立

`blom_problem_setup.py` 能稳定生成后续所有模块共享的数据结构。

### B. 全局参考层成立

`minco_scalar_baseline.py` 能给出稳定的 1D 全局参考解。

### C. 局部窗口层成立

`blom_local_qp.py` 能对 `s=4, k=2` 窗口问题稳定求解。

### D. 验证层成立

`phase0_checks.py` 能自动报告：

- 插值误差
- 连续性误差
- 边界条件残差

若以上四条未同时成立，则 Phase 0 不应视为完成。

---

## 10. 给 AI 的实现原则

请严格按以下原则实现：

1. 严格遵守 canonical setting：只做 `s=4, k=2, d_i=1, dim=1`
2. 优先保证接口统一，不优先追求高性能
3. 优先保证数学对象清晰，不隐藏 shape 与约束定义
4. 数值求解先行，解析推导后置
5. 所有模块必须可独立测试
6. 不要提前实现 obstacle cost / planner loop / gradient descent
7. 尽量保持模块边界清楚，便于后续用解析实现替换数值实现

---

## 11. 最小任务清单

Phase 0 的最小任务清单如下：

1. 完成 `blom_problem_setup.py`
2. 完成 `poly_basis.py`
3. 完成 `constraint_builders.py`
4. 完成 `minco_scalar_baseline.py`
5. 完成 `blom_local_qp.py`
6. 完成 `trajectory_eval.py`
7. 完成 `phase0_checks.py`
8. 完成 `demo_phase0.py`

---

## 12. 交付要求

最终交付应至少包含：

- 一份完整可运行的 Phase 0 代码目录
- 一个最小 demo
- 一组基础测试
- 清晰的 README 或注释说明每个模块的责任边界

如无特殊说明，默认 Python 版本为 `3.10+`。
