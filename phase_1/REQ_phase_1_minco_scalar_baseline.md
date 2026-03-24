# REQ-Phase1-minco_scalar_baseline.md

## 1. 文档目的

本需求文档用于指导实现 `minco_scalar_baseline.py`，作为 BLOM 理论 Phase 1 的**全局参考母问题求解器**。

该模块**不是** BLOM 本身，也**不是** GCOPTER 工程复现。  
它的作用是：

1. 在与 Phase 0 完全一致的外层符号和输入格式下，建立一个严格可解的全局 MINCO 型参考问题；
2. 提供从参数 `theta=(q,T)` 到全局最优多项式系数 `c(q,T)` 的确定性映射；
3. 提供后续 BLOM 所需的“参照链条”：
   - 母问题
   - 最优性与唯一性
   - 带状线性系统
   - 参数化
   - 梯度传播
4. 为后续 BLOM 的局部化改写、误差比较、收敛验证提供统一 baseline。

该设计遵循 MINCO 的系统构建纪律：先固定统一的多阶段最小控制问题，再由最优性条件得到唯一解、带状系统与线性复杂度操作，随后才引入稀疏参数化与时空形变。

---

## 2. 与 Phase 0 的统一要求

实现必须严格继承 Phase 0 的 canonical setting：

- 标量输出：
  \[
  p:[0,\mathcal T]\to\mathbb R
  \]
- 段数：\(M\)
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
- 控制阶次固定为：
  \[
  s=4
  \]
  即 minimum snap
- 中间点只固定 0 阶导：
  \[
  d_i=1,\qquad i=1,\dots,M-1
  \]
- 外层稀疏参数向量固定为：
  \[
  \theta=(q_1,\dots,q_{M-1},T_1,\dots,T_M)^\top
  \]
- \(q_0,q_M\) 视为固定数据，而不是优化变量。

---

## 3. 数学对象：本模块实际求解的“母问题”

### 3.1 关键修正

虽然 Phase 0 的全局参考问题只写了 waypoint + 全局 \(C^6\) 连续，但如果只固定位置而不固定初末端高阶 jet，则 minimum snap 问题一般**不唯一**。  
因此，本模块必须实现的是一个**加强版 MINCO 母问题**：

- 中间点仍然只固定位置；
- 但初末端必须额外给定到 \(s-1=3\) 阶导数的 boundary jet。

也就是说，本模块的输入除了 `q` 和 `T` 外，还必须包含固定的初末端边界数据：
\[
\zeta^-=(p(0),p'(0),p''(0),p^{(3)}(0)),
\qquad
\zeta^+=(p(\mathcal T),p'(\mathcal T),p''(\mathcal T),p^{(3)}(\mathcal T)).
\]

这一步不是实现细节，而是保证唯一性和后续梯度光滑性的必要条件。MINCO 的理论结构中，“存在唯一性”正是带状系统与参数灵敏度平滑性的基础。

### 3.2 母问题定义

实现对象为如下全局 minimum snap 问题：

\[
\min_{\{p_i\}_{i=1}^{M}}
\sum_{i=1}^{M}\int_0^{T_i}\left|p_i^{(4)}(t)\right|^2dt
\]

subject to

1. 每段 \(p_i\) 为次数不超过 7 的多项式；
2. 段端点插值：
   \[
   p_i(0)=q_{i-1},\qquad p_i(T_i)=q_i,\qquad i=1,\dots,M;
   \]
3. 中间结点全局 \(C^6\) 连续：
   \[
   p_i^{(\ell)}(T_i)=p_{i+1}^{(\ell)}(0),\qquad \ell=0,1,\dots,6,\quad i=1,\dots,M-1;
   \]
4. 初末端 boundary jet 固定：
   \[
   p^{(\ell)}(0)=\zeta^-_\ell,\qquad p^{(\ell)}(\mathcal T)=\zeta^+_\ell,\qquad \ell=0,1,2,3.
   \]

其中 \(\zeta^-_0=q_0,\ \zeta^+_0=q_M\)。  
该问题在当前设置下是 `minco_scalar_baseline.py` 的唯一数学求解对象。

---

## 4. 模块输出的数学意义

本模块输出的不是离散采样轨迹，而是**全局最优多项式系数矩阵**。

对每一段，采用标准单项式基：
\[
\beta(t)=(1,t,t^2,\dots,t^7)^\top.
\]
写成
\[
p_i(t)=c_i^\top \beta(t),\qquad c_i\in\mathbb R^8.
\]

将所有段系数堆叠为
\[
c=(c_1^\top,\dots,c_M^\top)^\top\in\mathbb R^{8M}.
\]

模块的核心输出应为：

- `coeffs`: 形状 `(M, 8)` 的 ndarray，对应每段 8 个系数；
- `c_vec`: 展平后的 `(8*M,)` 向量；
- 若启用梯度相关接口，则还应输出系统矩阵和分解缓存，以供后续 `W(q,T)=K(c(q,T),T)` 的梯度传播复用。

---

## 5. 代码必须实现的四个核心部分

## 5.1 Part A：问题构造（Mother Problem Builder）

### 目标

将输入的 `q, T, zeta_start, zeta_end` 转化为全局带状线性系统
\[
M(T)c(q,T)=b(q,\zeta^-,\zeta^+).
\]

### 要求

实现以下基础构件：

1. 单项式基函数
   ```python
   beta(t) -> shape (8,)
   ```
2. 各阶导数基函数
   ```python
   beta_d(t, order) -> shape (8,)
   ```
   其中 `order = 0,1,...,7`
3. 每段的系数到边界 jet 的线性映射
4. 全局矩阵 `M(T)` 的组装
5. 右端项 `b(...)` 的组装

### 数学对应

矩阵必须体现“相邻两段耦合”的带状结构。  
这是 MINCO 线性复杂度成立的根本结构基础。

### 最低验收要求

- 对任意 `M >= 1`，矩阵尺寸正确；
- 对任意 `T_i > 0`，矩阵条目有限；
- 所有约束行的物理意义在代码注释中必须清楚标明。

---

## 5.2 Part B：全局系数求解（Global Coefficient Solver）

### 目标

求解
\[
M(T)c=b
\]
并返回唯一的全局最优系数。

### 要求

实现函数：

```python
solve_minco_coefficients(q, T, zeta_start, zeta_end, *, return_system=False)
```

### 输入规范

- `q`: shape `(M+1,)`
- `T`: shape `(M,)`
- `zeta_start`: shape `(4,)`
- `zeta_end`: shape `(4,)`

其中必须检查：

- `len(q) == len(T) + 1`
- `np.all(T > 0)`
- `zeta_start[0] == q[0]`
- `zeta_end[0] == q[-1]`

### 输出规范

至少返回：

```python
{
    "coeffs": np.ndarray,   # shape (M, 8)
    "c_vec": np.ndarray,    # shape (8*M,)
    "M": np.ndarray,        # optional
    "b": np.ndarray,        # optional
}
```

### 求解要求

首版允许用 `numpy.linalg.solve`，但代码结构必须预留将来替换成 banded solver 的接口。  
理由：本模块首先服务于 Phase 1 理论验证，不必一开始就过度工程化；但系统矩阵本质上是带状的，后续 BLOM 和大规模 MINCO baseline 必须利用这一结构。

### 最低验收要求

- 小规模测试中求解成功；
- 残差
  \[
  \|Mc-b\|
  \]
  在数值精度范围内接近 0；
- 对相同输入多次运行输出一致。

---

## 5.3 Part C：轨迹评估接口（Trajectory Evaluation API）

### 目标

把“系数”转成“可查询轨迹对象”，便于后续验证和对比 BLOM。

### 要求

实现以下接口：

```python
evaluate_segment(coeff, t, order=0)
evaluate_trajectory(coeffs, T, t_global, order=0)
sample_trajectory(coeffs, T, num_per_segment=50, orders=(0,1,2,3,4))
```

### 数学对应

- `evaluate_segment` 对应
  \[
  p_i^{(\ell)}(t)=c_i^\top \beta^{(\ell)}(t)
  \]
- `evaluate_trajectory` 负责把全局时刻映射到对应段内局部时间

### 最低验收要求

- 可正确输出位置、速度、加速度、jerk、snap；
- 在结点左右极限上，可用于验证连续性。

---

## 5.4 Part D：参数敏感度与梯度传播准备（Sensitivity-Ready Design）

### 目标

虽然本阶段不一定完整实现高层优化，但代码必须**为后续梯度传播预留接口**。

MINCO 原始框架中，任何用户定义目标都写成
\[
W(q,T)=K(c(q,T),T),
\]
再由
\[
M(T)c(q,T)=b(q)
\]
推导
\[
\frac{\partial W}{\partial q},\qquad \frac{\partial W}{\partial T}.
\]

其关键是复用 `M` 的带状分解，避免显式求逆，从而保持线性复杂度。

### 要求

即使首版不完整实现梯度，也必须在代码结构上预留：

```python
build_system_matrix(T)
build_rhs(q, zeta_start, zeta_end)
solve_adjoint(M, dK_dc)
grad_wrt_q(...)
grad_wrt_T(...)
```

### 最低验收要求

- 文档中必须明确说明：`minco_scalar_baseline.py` 的系统矩阵接口未来会被 BLOM 用来做比较；
- 代码中必须有占位函数和 TODO 注释；
- 不允许把系统构造、求解、采样全部揉成一个不可拆的大函数。

---

## 6. 数学—代码对照表（实现时必须遵守）

| 数学对象 | 符号 | 代码对象 | 说明 |
|---|---|---|---|
| 段数 | \(M\) | `M_seg` | `len(T)` |
| 航点 | \(q_0,\dots,q_M\) | `q` | shape `(M+1,)` |
| 时间向量 | \(T_1,\dots,T_M\) | `T` | shape `(M,)` |
| 总时长 | \(\mathcal T\) | `T_total` | `T.sum()` |
| 段多项式 | \(p_i\) | `coeffs[i]` | 第 `i` 段 8 个系数 |
| 单项式基 | \(\beta(t)\) | `beta(t)` | shape `(8,)` |
| 导数基 | \(\beta^{(\ell)}(t)\) | `beta_d(t, ell)` | shape `(8,)` |
| 全局系数向量 | \(c\) | `c_vec` | shape `(8*M,)` |
| 系统矩阵 | \(M(T)\) | `M_mat` | shape `(8*M, 8*M)` |
| 右端项 | \(b(q)\) | `b_vec` | shape `(8*M,)` |
| 外层稀疏参数 | \(\theta=(q_1,\dots,q_{M-1},T)\) | 不直接展开 | 保持 Phase 0 记号一致 |
| 用户目标 | \(K(c,T)\) | `objective_fn(c, T)` | 后续阶段使用 |
| 外层目标 | \(W(q,T)\) | `W(theta)` | 后续阶段使用 |

---

## 7. 非目标（明确禁止）

本模块首版 **不做** 以下内容：

1. 不实现 BLOM 局部窗口求解；
2. 不实现自然边界的局部窗口 QP；
3. 不实现 Catmull-Rom 型局部速度公式；
4. 不实现障碍物约束、动力学约束、ESDF、平坦映射；
5. 不实现多维 \(\mathbb R^m\) 版本；
6. 不实现时间均匀 MINCO 的专门优化版；
7. 不直接写高层优化器（如 LBFGS）；
8. 不做 GPU、C++、并行化。

原因：本模块唯一任务是提供一个**严格、干净、可验证的全局参考母问题**。  
其他内容属于 Phase 2 以后。

---

## 8. 数值验证要求（必须实现）

实现完成后，至少提供如下测试函数或脚本：

### 8.1 插值验证

检查：
\[
p_i(0)=q_{i-1},\qquad p_i(T_i)=q_i
\]
最大误差应接近机器精度。

### 8.2 连续性验证

对每个内点 \(t_i\)，检查
\[
p_i^{(\ell)}(T_i)=p_{i+1}^{(\ell)}(0),\qquad \ell=0,\dots,6
\]
输出每一阶的最大 jump。

### 8.3 边界 jet 验证

检查
\[
p^{(\ell)}(0)=\zeta^-_\ell,\qquad
p^{(\ell)}(\mathcal T)=\zeta^+_\ell,\qquad \ell=0,1,2,3.
\]

### 8.4 系统残差验证

检查
\[
\|M(T)c-b\|
\]
是否足够小。

### 8.5 稳定性测试

随机生成多组 `q,T`，反复运行：

- 无 NaN / Inf
- 无奇异矩阵错误（在合理 `T_i>0` 条件下）
- 输出连续

### 8.6 与简单基准对比

至少包含以下 sanity cases：

#### Case A：单段

`M=1`  
给定初末端全 jet，返回唯一 7 次多项式。

#### Case B：两段对称

对称 waypoint 与对称时间配置下，轨迹应表现出可解释对称性。

#### Case C：时间缩放

将 `T` 按比例缩放时，检查系数与高阶导数尺度变化是否合理。

---

## 9. 建议的代码结构

```text
phase_1/
├── minco_scalar_baseline.py
├── test_minco_scalar_baseline.py
├── examples/
│   ├── demo_case_single_segment.py
│   ├── demo_case_two_segment.py
│   └── demo_case_random.py
└── README_phase1.md
```

`minco_scalar_baseline.py` 内建议至少包含：

```python
def beta(t: float) -> np.ndarray: ...
def beta_d(t: float, order: int) -> np.ndarray: ...

def build_system_matrix(T: np.ndarray) -> np.ndarray: ...
def build_rhs(q: np.ndarray,
              zeta_start: np.ndarray,
              zeta_end: np.ndarray) -> np.ndarray: ...

def solve_minco_coefficients(q: np.ndarray,
                             T: np.ndarray,
                             zeta_start: np.ndarray,
                             zeta_end: np.ndarray,
                             *,
                             return_system: bool = False) -> dict: ...

def evaluate_segment(coeff: np.ndarray,
                     t: float,
                     order: int = 0) -> float: ...

def evaluate_trajectory(coeffs: np.ndarray,
                        T: np.ndarray,
                        t_global: float,
                        order: int = 0) -> float: ...

def sample_trajectory(coeffs: np.ndarray,
                      T: np.ndarray,
                      num_per_segment: int = 50,
                      orders: tuple = (0,1,2,3,4)) -> dict: ...
```

---

## 10. 文档与注释要求

### 10.1 注释要求

注释必须说明“这一行对应哪条数学约束”，而不是只写“build matrix”。

例如：

- `# row block for initial jet constraints`
- `# continuity of derivative order ell at interior knot`
- `# endpoint interpolation p_i(T_i)=q_i`

### 10.2 README 必须包含

1. 数学对象说明  
2. 输入/输出形状说明  
3. 一个最小可复现示例  
4. 如何验证插值/连续性/边界 jet  
5. 已知限制：当前仅支持 scalar + \(s=4\)

---

## 11. Acceptance Criteria（验收标准）

只有同时满足以下条件，`minco_scalar_baseline.py` 才算完成：

1. **符号统一**：与 Phase 0 的 \(q,T,\theta,s,d_i\) 记号完全一致；
2. **对象正确**：实现的是“加强版 MINCO 母问题”，即含初末端 jet；
3. **求解正确**：对随机测试用例，插值误差、连续性误差、边界误差都接近机器精度；
4. **结构清晰**：系统构造、求解、轨迹评估、验证函数彼此分离；
5. **可扩展**：为后续 \(W(q,T)=K(c(q,T),T)\) 梯度传播预留接口；
6. **可对比**：输出格式足以作为 BLOM 后续局部窗口解的 reference baseline；
7. **可复现**：附带最小 demo 和测试脚本。

---

## 12. 给实现 AI 的最终指令

你要实现的是：

- 一个**标量、全局、minimum snap、带 boundary jet 的 MINCO baseline**
- 其外层输入格式仍然是 Phase 0 风格的 `q` 与 `T`
- 其输出是全局最优多项式系数 `c(q,T)`
- 它将作为后续 BLOM 所有局部化理论与数值实验的比较基准

请优先保证：

1. 数学对象定义正确；
2. 线性系统组装正确；
3. 插值与连续性验证正确；
4. 代码结构清晰；
5. 不要提前引入 BLOM 的局部窗口逻辑。

不要为了“看起来高级”而过早写复杂优化器、自动微分框架或工程封装。  
这一阶段最重要的是**正确、透明、可验证**。

