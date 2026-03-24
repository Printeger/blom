# REQ-Phase3-validation.md

## 1. 目标

本需求文档用于指导实现 **Phase 3 验证代码**，对应论文中的 **BLOM-Strict** 数学对象与 Theorem 2（局部问题存在唯一性）的数值验证。

本阶段代码不是为了完成 BLOM 全系统，而是为了验证以下命题是否与理论一致：

1. 局部问题的可行集构造正确；
2. 自然边界条件的数值表现与理论一致；
3. 在给定窗口、给定 `q,T`、给定边界 jet 下，局部问题解是唯一的；
4. 对输入的微小扰动，局部最优解连续变化；
5. 局部解满足插值与窗口内部连续性；
6. 局部问题可写成有限维 equality-constrained QP，并与 Phase 1 全局 MINCO baseline、Phase 2 必要性验证代码兼容。

---

## 2. 与 Phase 1 / Phase 2 的接口统一要求

### 2.1 与 Phase 1 的统一

必须复用或兼容 Phase 1 中的以下对象：

- `q = (q_0, ..., q_M)`
- `T = (T_1, ..., T_M)`
- `s = 4`（当前 canonical setting）
- 单项式基 `beta(t)` 与导数基 `beta_d(t, order)`
- 每段多项式系数表示：
  \[
  p_i(t)=c_i^\top \beta(t),\qquad c_i\in\mathbb R^{2s}
  \]
- 初末端 physical boundary jet：
  \[
  \zeta^-=(p(0),p'(0),\dots,p^{(s-1)}(0)),\qquad
  \zeta^+=(p(\mathcal T),p'(\mathcal T),\dots,p^{(s-1)}(\mathcal T))
  \]

### 2.2 与 Phase 2 的统一

必须复用或兼容 Phase 2 中的以下对象：

- 窗口定义
  \[
  W(i,k)=\{i-\lfloor k/2\rfloor,\dots,i+\lfloor k/2\rfloor\}\cap[1,M]
  \]
- canonical setting：
  \[
  s=4,\qquad k=2
  \]
- 内部窗口大小：
  \[
  |W(i,2)|=3=s-1
  \]
- Phase 2 中“不要试图直接局部化 exact MINCO map”的结论；Phase 3 验证代码只针对 **BLOM-Strict 的局部 variational problem**，不声称恢复全局 exact MINCO。

---

## 3. 本阶段需要验证的数学命题

### 3.1 命题 A：局部问题可行

对给定窗口 `W(i,k)`、给定 `q,T`、给定必要的物理边界 jet，局部 admissible set 非空。

需要数值上构造至少一个可行解，并验证其满足：

1. 窗口内所有 knot 的插值条件；
2. 窗口内部 `C^{s-1}` 连续；
3. 若窗口触碰真实全局边界，则满足 physical boundary jet；
4. 若窗口不触碰真实边界，则仅固定位置，不固定导数。

### 3.2 命题 B：局部问题唯一

对固定输入，局部 BLOM-Strict 问题有唯一最优解。

需要通过两条路径验证：

1. **优化角度**：同一输入、多次随机初始化，求得相同系数；
2. **线性代数角度**：reduced Hessian 正定，或 KKT 系统在约束子空间上非奇异。

### 3.3 命题 C：自然边界条件成立

对人工窗口边界，若该边界不是真实物理边界，则局部最优解应满足理论给出的自然边界条件。当前 canonical case 为 `s=4`，因此人工边界应数值接近：

\[
 p^{(4)}=p^{(5)}=p^{(6)}=0
\]

在窗口左右人工边界分别检查。

### 3.4 命题 D：窗口内部高阶连续性成立

局部问题定义时只强制 `C^{s-1}`，但理论上最优解应在窗口内部 knot 自动提升到
\[
C^{2s-2}
\]
。

当前 `s=4` 时，需检查内部 knot 处：

\[
\llbracket p^{(r)} \rrbracket = 0,
\qquad r=0,1,2,3,4,5,6.
\]

### 3.5 命题 E：局部映射对小扰动连续

对固定窗口和固定 `T`，对 waypoint 或 boundary jet 做小扰动，局部解系数应连续变化。

至少需要验证：

1. 对单个 waypoint 加入微扰 `eps` 时，系数变化量随 `eps` 缩小而缩小；
2. 不出现明显的跳变、分支切换、数值不稳定；
3. KKT / reduced system 的条件数在合理范围内。

---

## 4. 必须实现的代码模块

## 4.1 主模块

### 文件名

`blom_strict_local_qp.py`

### 作用

给定窗口 `W(i,k)`、全局 `q,T`、边界 jet，构造并求解 BLOM-Strict 局部 equality-constrained QP。

### 建议函数接口

```python
build_window(i: int, k: int, M: int) -> dict
build_local_problem(q, T, i, k, zeta_start=None, zeta_end=None) -> dict
solve_blom_strict_local_qp(problem: dict, method: str = "kkt") -> dict
extract_segment_coeff(solution: dict, seg_idx: int) -> np.ndarray
```

---

## 4.2 可行解构造模块

### 文件名

`blom_strict_feasible_init.py`

### 作用

显式构造一个满足局部约束的可行解，用于：

- 检查 admissible set 非空；
- 给局部 QP 提供初始点；
- 单独测试约束构造是否正确。

### 建议函数接口

```python
build_feasible_local_spline(q, T, i, k, zeta_start=None, zeta_end=None) -> dict
```

实现建议：

- 先给每个窗口 knot 指定 0 到 `s-1` 阶 jet；
- 内部 knot 的高阶导可先置零；
- 每段用 Hermite 多项式拼出 degree `2s-1` 的可行段。

---

## 4.3 局部 QP / KKT 模块

### 文件名

`blom_strict_local_kkt.py`

### 作用

将局部问题写成
\[
\min \frac12 c^\top H c \quad \text{s.t.}\quad Gc=d
\]
并支持以下两种解法：

1. 直接解 KKT 系统；
2. Null-space reduction 解 reduced Hessian。

### 建议函数接口

```python
build_local_hessian(T_window, s=4) -> np.ndarray
build_local_constraints(q, T, i, k, zeta_start=None, zeta_end=None, s=4) -> tuple[np.ndarray, np.ndarray]
solve_kkt(H, G, d) -> dict
solve_reduced_qp(H, G, d) -> dict
```

### 说明

必须保留两套解法，原因：

- KKT 解法便于调试与检查约束；
- reduced Hessian 解法便于验证“唯一性来自约束子空间上的正定性”。

---

## 4.4 理论验证脚本

### 文件名

`validate_phase3_blom_strict.py`

### 作用

批量运行所有 Phase 3 验证项，并保存数值结果与结果图片。

### 建议功能

1. 遍历多个 `M`、多个窗口位置 `i`；
2. 区分三类窗口：
   - 左边界窗口
   - 内部窗口
   - 右边界窗口
3. 输出定量指标表；
4. 自动生成图片；
5. 将结果保存到固定目录。

---

## 5. 数学—代码对照表

| 数学对象 | 符号 | 代码对象 | 说明 |
|---|---|---|---|
| 控制阶次 | \(s\) | `s` | 当前固定为 4 |
| 窗口索引 | \(W(i,k)\) | `window["segments"]` | 局部段索引集合 |
| 窗口左/右边界 | \(L_i,R_i\) | `window["L"]`, `window["R"]` | |
| 窗口大小 | \(m_i\) | `window["m"]` | `R-L+1` |
| 局部时间轴 | \(\sigma_r\) | `window["sigma"]` | 局部 knot 时间 |
| 局部 horizon | \(\mathcal T_{i,k}\) | `window["T_local"]` | |
| 局部 spline 空间 | \(\mathcal S_{i,k}^{(s-1)}\) | 隐式由 `Gc=d` 表示 | |
| 局部 admissible set | \(\mathcal A_{i,k}\) | `problem["G"]`, `problem["d"]` | affine set |
| 局部目标泛函 | \(J_{i,k}[u]\) | `H` | quadratic form |
| 局部系数 | \(c^{loc}_{i,k}\) | `c_loc` | shape `(2*s*m_i,)` |
| reduced Hessian | \(R_{i,k}\) | `R_red` | null-space reduction 后的 Hessian |
| 中心段系数提取 | \(\Pi_i\) | `extract_segment_coeff` | 提取窗口中心段 |

---

## 6. 必须保存的结果图片

所有图片必须自动保存，不能只在 notebook 中显示。

### 6.1 图片保存目录

```text
phase_3/results/
├── figures/
├── tables/
└── logs/
```

### 6.2 必须生成的图片

#### Figure 1：窗口示意图

文件名：

```text
fig_window_layout_M{M}_i{i}_k{k}.png
```

内容：

- 全局段编号；
- 当前窗口 `W(i,k)` 高亮；
- 左右是否触碰物理边界的标记；
- 中心段位置标记。

用途：直观展示当前局部问题对应哪一段。

---

#### Figure 2：局部轨迹及导数图

文件名：

```text
fig_local_traj_M{M}_i{i}_k{k}.png
```

内容：

- 局部窗口内的 `p, p', p'', p''', p''''` 至少 5 条曲线；
- knot 位置竖线标出；
- 插值点高亮；
- 左右人工边界或物理边界标记。

用途：直观看局部最优曲线是否合理。

---

#### Figure 3：内部连续性 jump 图

文件名：

```text
fig_continuity_jump_M{M}_i{i}_k{k}.png
```

内容：

- 横轴为导数阶数 `r = 0,...,2s-2`；
- 纵轴为每阶 jump 的绝对值；
- 对每个内部 knot 可单独一组柱状图或折线图。

用途：验证最优解是否确实自动提升为 `C^{2s-2}`。

---

#### Figure 4：自然边界条件残差图

文件名：

```text
fig_natural_bc_residual_M{M}_i{i}_k{k}.png
```

内容：

- 左右人工边界上 `p^{(s)},...,p^{(2s-2)}` 的残差大小；
- 若某侧是真实物理边界，则标记为“不适用”。

用途：验证 Proposition 级别的自然边界条件。

---

#### Figure 5：扰动连续性响应图

文件名：

```text
fig_perturbation_response_M{M}_i{i}_k{k}.png
```

内容：

- 横轴：扰动幅度 `eps`；
- 纵轴：局部系数变化范数 `||delta c||`；
- 至少用 log-log 图展示；
- 可附加条件数变化曲线。

用途：验证局部映射对输入小扰动连续。

---

#### Figure 6：唯一性/初始化鲁棒性对比图

文件名：

```text
fig_uniqueness_multistart_M{M}_i{i}_k{k}.png
```

内容：

- 多次不同初始值求解得到的系数差异；
- 或最终 objective 差异；
- 或 pairwise coefficient distance heatmap。

用途：从数值角度直观看“唯一解”不是初始化伪象。

---

## 7. 必须保存的表格结果

### 表格目录

```text
phase_3/results/tables/
```

### 必须输出 CSV

#### Table A：局部可行性与约束残差

文件名：

```text
table_feasibility_summary.csv
```

字段至少包括：

- `M`
- `i`
- `k`
- `window_type` (`left` / `interior` / `right`)
- `max_interp_error`
- `max_Cs_minus_1_jump`
- `max_C2s_minus_2_jump`
- `max_boundary_jet_error`
- `max_natural_bc_residual`
- `solver_status`

#### Table B：唯一性与数值稳定性

文件名：

```text
table_uniqueness_summary.csv
```

字段至少包括：

- `M`
- `i`
- `k`
- `num_restarts`
- `max_pairwise_coeff_diff`
- `objective_std`
- `reduced_hessian_min_eig`
- `reduced_hessian_cond`
- `kkt_residual`

#### Table C：扰动连续性

文件名：

```text
table_perturbation_continuity.csv
```

字段至少包括：

- `M`
- `i`
- `k`
- `perturb_target`
- `eps`
- `coeff_diff_norm`
- `traj_diff_sup_norm`
- `cond_number`

---

## 8. 建议目录结构

```text
phase_3/
├── blom_strict_local_qp.py
├── blom_strict_local_kkt.py
├── blom_strict_feasible_init.py
├── validate_phase3_blom_strict.py
├── test_phase3_blom_strict.py
├── results/
│   ├── figures/
│   ├── tables/
│   └── logs/
└── README_phase3.md
```

---

## 9. 最低数值验证要求

### 9.1 Case 1：canonical interior window

- `s=4, k=2`
- 取一个内部窗口，例如 `i=3`
- 验证：
  - 可行解存在；
  - 局部最优解唯一；
  - 内部 knot 达到 `C^6`；
  - 人工边界满足自然边界条件。

### 9.2 Case 2：left boundary window

- `i=1`
- 验证：
  - 左侧继承 physical boundary jet；
  - 右侧人工边界满足自然边界条件；
  - 唯一性依然成立。

### 9.3 Case 3：right boundary window

- `i=M`
- 验证右边界与左边界对称结果。

### 9.4 Case 4：扰动试验

- 对某一个 interior waypoint `q_j` 加入 `eps = 1e-1,1e-2,...,1e-6`
- 检查 `||delta c||` 与 `eps` 是否随 `eps -> 0` 一起趋于 0。

### 9.5 Case 5：多初始化唯一性

- 采用多个随机可行初值或多个随机初始化；
- 检查最终系数是否一致到数值精度。

---

## 10. 与 Phase 1 / Phase 2 的联动验证

### 10.1 与 Phase 1 对照

对于足够大的窗口（例如窗口覆盖全轨迹时），局部问题应退化到对应的全局问题。允许只做小规模数值检查，不要求在本阶段给出一般性定理。

### 10.2 与 Phase 2 对照

本阶段不验证“不可能性定理”，但需要体现其后果：

- Phase 2 说明不能直接从 exact MINCO map 获得 finite local support；
- Phase 3 则验证：只要改成局部 variational problem，局部问题本身是良定且唯一的。

因此 README 中必须明确说明：

> Phase 3 的数值验证目标是“局部问题成为合法数学对象”，而不是“恢复 exact global MINCO”。

---

## 11. 非目标

本阶段**不做**以下内容：

1. 不做 Phase 4 的解析显式公式推导；
2. 不做收敛到全局 MINCO 的误差界；
3. 不做全局装配后的整条 BLOM 轨迹一致性；
4. 不做 obstacle cost、dynamics feasibility、flatness backward pass；
5. 不做多维版本；
6. 不做自动微分框架优化；
7. 不做 BLOM-Analytic，只做 BLOM-Strict 的数值验证。

---

## 12. README 必须写清楚的内容

`README_phase3.md` 必须包含：

1. 当前数学对象：BLOM-Strict local variational problem；
2. 为什么只要求 local well-posedness，不讨论全局装配；
3. 输入输出说明；
4. 运行验证脚本的方法；
5. 所有结果图片与 CSV 的目录说明；
6. 如何解读 uniqueness / natural BC / continuity / perturbation 四类图。

---

## 13. Acceptance Criteria

只有同时满足以下条件，Phase 3 验证代码才算完成：

1. 能正确构造局部 admissible problem；
2. 至少两套求解方法（KKT / reduced）结果一致；
3. 能稳定验证 canonical `s=4,k=2` 情形；
4. 局部插值误差、约束残差、KKT 残差接近机器精度；
5. reduced Hessian 数值上正定；
6. 多初始化得到相同解；
7. 扰动响应连续；
8. 自动保存规定的 PNG / CSV 文件；
9. 与 Phase 1/2 的符号、接口、目录命名保持统一。

---

## 14. 给实现 AI 的最终指令

你要实现的是：

- 一个 **BLOM-Strict 局部 equality-constrained QP 求解与验证套件**；
- 它服务于 Phase 3：证明“局部问题是合法数学对象”；
- 它必须兼容 Phase 1 的全局 MINCO baseline 和 Phase 2 的必要性理论；
- 它必须自动保存所有关键数值结果和图片，便于后续论文展示；
- 当前优先级是：**正确、清晰、可验证**，不是性能极致或大规模工程化。
