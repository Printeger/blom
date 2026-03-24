# REQ-Phase2-validation.md

## 1. 文档目的

本需求文档用于指导实现 **Phase 2：BLOM 必要性理论验证代码**，作为 Phase 2 理论章节的数值与结构性支撑工具。

本模块的目标不是完成 BLOM 本体求解，而是围绕下列两个结论提供可复现、可视化、可量化的验证：

1. **弱版结论**：对于 Phase 1 的全局 MINCO 型精确系数映射
   \[
   c^\star(q,T)=A(T)^{-1}b(q,T),
   \]
   若
   \[
   D_q c^\star = A(T)^{-1}S_q
   \]
   不具有严格有限带宽，则在当前稀疏 waypoint 参数化下不存在 exact finite local support。
2. **强版结论**：在 Phase 0 / Phase 1 的统一参数化下，若 exact global optimum 由全局带状系统唯一确定，则“精确全局最优 + 稀疏 waypoint 参数化 + 统一有限局部支撑”不能同时成立。

本模块的主要作用是：

- 从 Phase 1 的 `minco_scalar_baseline.py` 出发，计算精确全局系数映射对 waypoint 的灵敏度；
- 通过解析 Jacobian、有限差分、影响衰减曲线、热力图与带宽统计，验证 Phase 2 的必要性理论；
- 保存所有关键结果图片，供论文图、汇报图、debug 图使用；
- 为后续 BLOM 与 MINCO 的对比实验提供结构性证据。

---

## 2. 与 Phase 0 / Phase 1 的统一要求

实现必须严格继承前两阶段的记号与对象：

- 标量输出：
  \[
  p:[0,\mathcal T]\to\mathbb R
  \]
- 控制阶次固定为：
  \[
  s=4
  \]
- 中间点仅固定 0 阶导：
  \[
  d_i=1
  \]
- 外层稀疏参数向量：
  \[
  \theta=(q_1,\dots,q_{M-1},T_1,\dots,T_M)^\top
  \]
- 全局参考母问题由 Phase 1 的 `minco_scalar_baseline.py` 提供：
  \[
  A(T)c^\star(q,T)=b(q,T)
  \]
- 其中 `q_0, q_M` 以及初末端 boundary jets 为固定数据，不视为优化变量。

**禁止** 在 Phase 2 验证代码中擅自改动 Phase 1 的数学对象；Phase 2 仅允许在其上做灵敏度分析、带宽分析与可视化，不允许把问题改成局部窗口版本。

---

## 3. Phase 2 需要验证的核心问题

## 3.1 弱版定理对应的数值验证目标

需要验证：

1. 精确系数映射对内点 waypoint 的 Jacobian
   \[
   D_{\bar q}c^\star = A(T)^{-1}S_q
   \]
   是否严格有限带宽；
2. 对于固定一个 waypoint 扰动 \(\Delta q_j\)，其影响是否传播到远处段系数 \(c_i\)；
3. 这种传播是否只是“衰减”而不是“严格为零”；
4. 结论是否随不同 `M`、不同 `T`、不同 waypoint 分布保持稳定。

## 3.2 强版定理对应的数值验证目标

需要验证：

1. 不存在一个与 `M` 无关的统一带宽 `w`，使得所有测试样例中
   \[
   \frac{\partial c_i}{\partial q_j}=0 \quad \text{for } |i-j|>w
   \]
   都严格成立；
2. 当 `M` 增大时，远距非零影响依然存在；
3. 这种现象不是数值误差造成的，而是精确全局系统的结构性结果；
4. `A(T)` 自身是带状的，但 `A(T)^{-1}S_q` 一般不是严格带状。

---

## 4. 实现范围

本模块应实现以下四类功能：

1. **系统级灵敏度提取**
2. **有限差分对照验证**
3. **局部支撑失效的可视化**
4. **结果保存与批量实验**

本模块不实现：

- BLOM 局部窗口求解；
- 局部 QP；
- Catmull-Rom 或 BLOM 解析局部公式；
- 障碍物、动力学约束、ESDF；
- 多维轨迹版本；
- 高层优化器。

---

## 5. 代码模块要求

## 5.1 Part A：解析 Jacobian 提取

### 目标

基于 Phase 1 的全局系统
\[
A(T)c^\star(q,T)=b(q,T)
\]
提取 waypoint 方向的精确 Jacobian：
\[
D_{\bar q}c^\star = A(T)^{-1}S_q.
\]

### 要求

实现函数：

```python
build_waypoint_selector(M: int, s: int = 4) -> np.ndarray
compute_exact_jacobian_q(M_mat: np.ndarray, S_q: np.ndarray) -> np.ndarray
```

### 说明

- `S_q` 是 Phase 1 右端项对 interior waypoints 的选择矩阵；
- `J_q = A^{-1}S_q` 的形状应为 `(2*s*M, M-1)`；
- 对当前 canonical setting，形状应为 `(8*M, M-1)`。

### 最低验收要求

- 维度正确；
- 与直接逐列求解 `A x = e_j` 的结果一致；
- 不显式求逆，优先通过线性求解实现；
- 保持与 Phase 1 中 `build_system_matrix` 的接口兼容。

---

## 5.2 Part B：有限差分 Jacobian 对照验证

### 目标

通过有限差分验证解析 Jacobian 的正确性，避免将后续现象误判为实现错误。

### 要求

实现函数：

```python
finite_difference_jacobian_q(
    q: np.ndarray,
    T: np.ndarray,
    zeta_start: np.ndarray,
    zeta_end: np.ndarray,
    eps: float = 1e-6,
) -> np.ndarray
compare_exact_vs_fd_jacobian(...)
```

### 说明

对每个 interior waypoint `q_j` 做微小扰动：
\[
q_j \mapsto q_j + \varepsilon,
\]
用 Phase 1 baseline 重算 `c`，估计
\[
\frac{\partial c}{\partial q_j}
\approx
\frac{c(q+\varepsilon e_j,T)-c(q,T)}{\varepsilon}.
\]

### 最低验收要求

- 解析 Jacobian 与有限差分在合理精度范围内一致；
- 可输出逐列误差、Frobenius 范数误差、相对误差。

---

## 5.3 Part C：影响传播与局部支撑失效分析

### 目标

直观展示：“改变一个 waypoint 会影响多远的段系数”。

### 要求

实现函数：

```python
compute_waypoint_influence_profile(J_q: np.ndarray, M_seg: int, target_waypoint_idx: int) -> dict
compute_segmentwise_influence_norms(J_q: np.ndarray, M_seg: int) -> np.ndarray
estimate_effective_bandwidth(J_q: np.ndarray, tol: float) -> dict
```

### 建议定义

对给定 interior waypoint `q_j`，定义第 `i` 段的影响强度为：
\[
I_{ij} := \left\|\frac{\partial c_i}{\partial q_j}\right\|_2.
\]
其中 `c_i` 表示第 `i` 段 8 维系数块。

### 需要输出的量

1. 单个 waypoint 的 influence profile：
   \[
   i \mapsto I_{ij}
   \]
2. 全部 waypoint–segment 影响矩阵：
   \[
   \mathcal I = (I_{ij})
   \]
3. 在给定阈值 `tol` 下的“表观带宽”统计；
4. `M` 增大时最大影响距离是否增长。

### 最低验收要求

- 至少支持 `M=4, 8, 16, 32`；
- 支持均匀时间与非均匀时间样例；
- 支持不同随机 waypoint 样例；
- 能够证明“远距影响往往衰减但不严格为 0”。

---

## 5.4 Part D：带状系统 vs 稠密逆效应展示

### 目标

把 Phase 2 的核心逻辑可视化：

- `A(T)` 是带状的；
- 但 `A(T)^{-1}S_q` 一般不是严格带状；
- 因此 exact coefficient map 不具有统一 finite local support。

### 要求

实现函数：

```python
visualize_matrix_sparsity(mat: np.ndarray, save_path: str, title: str) -> None
visualize_block_sparsity_Jq(J_q: np.ndarray, M_seg: int, save_path: str) -> None
```

### 需要生成的图

1. `A(T)` 的稀疏模式图；
2. `J_q = A^{-1}S_q` 的稀疏/幅值热力图；
3. block-level 影响热力图（segment × waypoint）。

### 最低验收要求

- 图像分辨率足够用于论文插图与汇报；
- 文件名规范；
- 图上必须有标题、坐标轴或 colorbar；
- 图像可自动保存。

---

## 6. 建议实现的 Python 文件

```text
phase_2/
├── phase2_validation.py
├── phase2_plotting.py
├── test_phase2_validation.py
├── examples/
│   ├── demo_phase2_uniform_time.py
│   ├── demo_phase2_nonuniform_time.py
│   ├── demo_phase2_scaling_M.py
│   └── demo_phase2_fd_check.py
├── results/
│   ├── figures/
│   ├── tables/
│   └── logs/
└── README_phase2.md
```

### 文件职责建议

#### `phase2_validation.py`
负责：
- `S_q` 构造
- 解析 Jacobian 提取
- 有限差分 Jacobian
- influence profile 计算
- 有效带宽估计

#### `phase2_plotting.py`
负责：
- 矩阵稀疏图
- 热力图
- 影响曲线图
- `M` 扩展实验结果图

#### `test_phase2_validation.py`
负责：
- 解析 Jacobian 与有限差分一致性测试
- 维度与索引正确性测试
- 图像生成 smoke test

---

## 7. 必须保存的结果图片

本模块必须自动保存至少以下图片，供论文、汇报和 debug 使用。

## 7.1 图 1：全局系统矩阵稀疏模式图

### 文件名

```text
A_sparsity_M{M}_case{case_name}.png
```

### 内容

展示 `A(T)` 的稀疏模式，突出其带状结构。

### 目的

支持“全局最优系统本身是局部耦合”的论述。

---

## 7.2 图 2：解析 Jacobian 热力图

### 文件名

```text
Jq_heatmap_M{M}_case{case_name}.png
```

### 内容

展示
\[
J_q=A^{-1}S_q
\]
的元素幅值热力图，建议使用 `log10(abs(J_q)+eps)`。

### 目的

直观展示“逆诱导的全局耦合”。

---

## 7.3 图 3：block-level influence 热力图

### 文件名

```text
block_influence_M{M}_case{case_name}.png
```

### 内容

展示
\[
I_{ij}=\left\|\frac{\partial c_i}{\partial q_j}\right\|_2
\]
的 segment × waypoint 热力图。

### 目的

比原始 Jacobian 更适合解释 local support 是否失效。

---

## 7.4 图 4：单个 waypoint 的影响衰减曲线

### 文件名

```text
influence_profile_q{j}_M{M}_case{case_name}.png
```

### 内容

固定一个 interior waypoint `q_j`，画出其对各段系数块的影响范数：
\[
i \mapsto I_{ij}.
\]

### 目的

展示“远处影响通常衰减但不严格为 0”。

---

## 7.5 图 5：有限差分 vs 解析 Jacobian 对照图

### 文件名

```text
jacobian_fd_compare_M{M}_case{case_name}.png
```

### 内容

展示解析 Jacobian 与有限差分 Jacobian 的元素散点图或误差图。

### 目的

证明数值现象不是代码错误。

---

## 7.6 图 6：随 M 扩展的“表观带宽”统计图

### 文件名

```text
effective_bandwidth_vs_M_case{case_name}.png
```

### 内容

在给定阈值 `tol` 下，统计 `J_q` 的有效带宽随 `M` 的变化趋势。

### 目的

支持强版定理的数值证据：不存在与 `M` 无关的统一有限带宽。

---

## 8. 建议保存的表格结果

除了图片，还应保存 CSV 或 JSON 结果表。

## 8.1 Jacobian 对照误差表

```text
jacobian_error_summary.csv
```

列建议包括：

- `case_name`
- `M`
- `eps`
- `fro_error`
- `max_abs_error`
- `relative_error`

## 8.2 有效带宽统计表

```text
effective_bandwidth_summary.csv
```

列建议包括：

- `case_name`
- `M`
- `tol`
- `max_effective_bandwidth`
- `mean_effective_bandwidth`
- `far_nonzero_ratio`

## 8.3 远距影响统计表

```text
far_field_influence_summary.csv
```

列建议包括：

- `case_name`
- `M`
- `target_waypoint_idx`
- `max_far_influence`
- `mean_far_influence`
- `num_far_entries_above_tol`

---

## 9. 建议的最小实验设计

## 9.1 均匀时间样例

- `T_i = 1`；
- waypoint 取平滑但非线性分布；
- 用于观察最基础结构。

## 9.2 非均匀时间样例

- `T_i` 在合理范围内随机采样；
- 验证现象不是均匀时间特例。

## 9.3 对称样例

- 对称 waypoint + 对称时间；
- 用于 sanity check 与图像解释。

## 9.4 随 M 扩展样例

- `M = 4, 8, 16, 32`；
- 核心用于强版结论的数值支撑。

---

## 10. 实现时必须提供的核心函数接口

```python
def build_waypoint_selector(M_seg: int, s: int = 4) -> np.ndarray: ...

def compute_exact_jacobian_q(M_mat: np.ndarray,
                             S_q: np.ndarray) -> np.ndarray: ...

def finite_difference_jacobian_q(q: np.ndarray,
                                 T: np.ndarray,
                                 zeta_start: np.ndarray,
                                 zeta_end: np.ndarray,
                                 eps: float = 1e-6) -> np.ndarray: ...

def compare_exact_vs_fd_jacobian(J_exact: np.ndarray,
                                 J_fd: np.ndarray) -> dict: ...

def compute_segmentwise_influence_norms(J_q: np.ndarray,
                                        M_seg: int,
                                        block_size: int = 8) -> np.ndarray: ...

def compute_waypoint_influence_profile(J_q: np.ndarray,
                                       M_seg: int,
                                       target_waypoint_idx: int,
                                       block_size: int = 8) -> dict: ...

def estimate_effective_bandwidth(J_q: np.ndarray,
                                 M_seg: int,
                                 tol: float,
                                 block_size: int = 8) -> dict: ...

def visualize_matrix_sparsity(mat: np.ndarray,
                              save_path: str,
                              title: str) -> None: ...

def visualize_heatmap(mat: np.ndarray,
                      save_path: str,
                      title: str,
                      log_scale: bool = True) -> None: ...

def run_phase2_validation_suite(...) -> dict: ...
```

---

## 11. 数学—代码对照表

| 数学对象 | 符号 | 代码对象 | 说明 |
|---|---|---|---|
| 全局系统矩阵 | \(A(T)\) | `M_mat` 或 `A_mat` | 直接复用 Phase 1 系统矩阵 |
| 右端项 | \(b(q,T)\) | `b_vec` | 与 Phase 1 保持一致 |
| waypoint 选择矩阵 | \(S_q\) | `S_q` | interior waypoint 对 RHS 的线性映射 |
| 精确系数向量 | \(c^\star(q,T)\) | `c_vec` | 由 Phase 1 baseline 得到 |
| waypoint Jacobian | \(D_{\bar q} c^\star\) | `J_q` | 理论上等于 `A^{-1}S_q` |
| 段块系数 | \(c_i\) | `coeffs[i]` | 每段 8 维 |
| 影响强度 | \(I_{ij}\) | `I[i, j]` | block-level 2 范数 |
| 表观带宽 | - | `effective_bandwidth` | 给定阈值下估计 |

---

## 12. 验收标准（Acceptance Criteria）

只有同时满足以下条件，Phase 2 验证代码才算完成：

1. **统一性**：完全兼容 Phase 1 的 `minco_scalar_baseline.py`；
2. **正确性**：解析 Jacobian 与有限差分 Jacobian 一致；
3. **可视化充分**：至少自动保存 6 类结果图；
4. **可复现**：固定随机种子，重复运行结果一致；
5. **结构清晰**：验证逻辑、绘图逻辑、批量实验逻辑彼此分离；
6. **可支撑理论**：能清楚展示 `A(T)` 带状，但 `A^{-1}S_q` 一般不严格带状；
7. **可汇报**：生成的图片和表格可以直接用于展示 Phase 2 结论。

---

## 13. 给实现 AI 的最终指令

你要实现的是：

- 一个 **Phase 2 必要性理论验证工具箱**；
- 它建立在 Phase 1 的全局 MINCO baseline 上；
- 它不求 BLOM，只验证“为什么 exact MINCO coefficient map 不能直接局部化”；
- 它必须同时提供：
  - 解析 Jacobian
  - 有限差分验证
  - 影响传播分析
  - 结构图与热力图
  - 自动保存图片与结果表

请优先保证：

1. 数学定义与 Phase 1 完全一致；
2. 代码接口干净、易复用；
3. 结果图自动保存；
4. 便于后续在论文和组会中直接使用。

不要提前引入 BLOM 局部窗口逻辑，也不要把本模块写成高层优化器。  
这一阶段最重要的是：**把必要性理论的数值证据做扎实、做透明、做可复现。**
