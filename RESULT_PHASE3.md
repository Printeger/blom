# RESULT_PHASE3

本文档记录 2026-03-24 对 `phase_3/` 代码的实际运行结果。  
结构统一为“步骤 -> 功能说明 -> 实际结果 -> 分析”，所有数值均来自本次实际执行。

Phase 3 的目标不是恢复全局 exact MINCO，而是验证 **BLOM-Strict 局部 variational problem** 的数学良定性：

- 局部可行集非空；
- 局部 equality-constrained QP 可正确构造；
- KKT 与 reduced 两种解法一致；
- 人工边界满足自然边界条件；
- 窗口内部自动提升到高阶连续；
- 对小扰动连续；
- 多初始化不会导向不同局部最优解。

本次代表性 `M = 5` 案例使用：

- `q = [0.0, -0.4665058580777545, 0.06821567603647731, -1.329752780152007, -1.3580488311632326, -0.9873633963907998]`
- `T = [0.8745252246154968, 0.7902320858528984, 1.2504314324351662, 1.009202263409843, 0.7754680383521799]`
- `zeta_start = [q_0, 0, 0, 0]`
- `zeta_end = [q_M, 0, 0, 0]`

批量验证脚本额外覆盖了 `M = 6` 的 left / interior / right 窗口。

## Step 1. 静态可执行性检查

### 运行命令

```bash
python3 -m compileall phase_3
```

### 功能说明

确认 `phase_3/` 下的主模块、验证脚本、绘图模块和测试模块都能被 Python 正常解析和导入。

### 实际结果

```text
Listing 'phase_3'...
Compiling 'phase_3/blom_strict_local_kkt.py'...
Compiling 'phase_3/phase3_plotting.py'...
Listing 'phase_3/results'...
Listing 'phase_3/results/figures'...
Listing 'phase_3/results/logs'...
Listing 'phase_3/results/tables'...
```

### 分析

- `compileall` 正常结束，说明 Phase 3 代码的语法和模块组织是可执行的。
- `results/` 目录结构已就位，后续验证脚本可以直接产出图片、CSV 和日志。

## Step 2. canonical interior window 的局部问题构造与可行解

### 运行方式

通过 Python 调用：

```python
problem = build_local_problem(q, T, i=3, k=2, zeta_start=zeta_start, zeta_end=zeta_end)
```

### 功能说明

构造 canonical `s = 4, k = 2` interior window 的局部 QP 数据，包括：

- 窗口 `W(3,2) = {2,3,4}`
- 局部 Hessian `H`
- 局部约束 `Gc = d`
- 显式可行初值 `feasible_init`

### 实际结果

- `window = [2, 3, 4]`
- `H.shape = (24, 24)`
- `G.shape = (12, 24)`
- `feasible_init_residual = 9.77032209594446e-13`

### 分析

- interior window 含 3 段，每段 8 个系数，因此未知量是 `24`，尺寸正确。
- `G.shape = (12, 24)` 对应 6 条位置约束和 6 条 `C^3` 连续约束，符合 Phase 3 的局部 admissible set 定义。
- 显式构造的 feasible init 残差在 `1e-12` 量级，说明 admissible set 数值上非空，命题 A 得到了直接支撑。

## Step 3. interior window 的 KKT / reduced 一致性与唯一性

### 运行方式

```python
sol_kkt = solve_blom_strict_local_qp(problem, method="kkt")
sol_red = solve_blom_strict_local_qp(problem, method="reduced")
```

### 功能说明

验证同一个 interior 局部问题在两条路径下得到相同解：

- 直接解 KKT saddle-point system
- 通过 null-space reduction 解 reduced Hessian 系统

同时检查 reduced Hessian 是否正定，从线性代数角度支撑“唯一最优解”。

### 实际结果

- `objective_kkt = 0.0`
- `objective_reduced = 1.1396929361370126e-24`
- `coeff_diff_norm = 6.00667770482123e-13`
- `max_interp_error = 4.85722573273506e-16`
- `max_Cs_minus_1_jump = 4.85722573273506e-16`
- `max_C2s_minus_2_jump = 4.85722573273506e-16`
- `max_natural_bc_residual = 0.0`
- `reduced_hessian_min_eig = 5.426112236089494`
- `reduced_hessian_cond = 32471.877144141665`

### 分析

- 两套解法的目标值和系数差异几乎都在机器精度范围内，说明局部 QP 构造和求解是一致的。
- `reduced_hessian_min_eig > 0`，数值上直接支持“约束子空间上严格凸”，也就是命题 B。
- 这个 interior case 的目标值接近 `0`，说明该窗口数据恰好允许一个零 snap 的低阶多项式解；这并不破坏唯一性，反而说明 reduced system 成功锁定了唯一的零代价可行解。

## Step 4. left boundary window 的 physical boundary inheritance 与自然边界

### 运行方式

对 `M = 5, i = 1, k = 2` 运行 local solve。

### 功能说明

验证左边界窗口同时满足两类条件：

- 左侧继承 Phase 1 的 physical boundary jet
- 右侧人工边界自动满足自然边界条件

### 实际结果

- `window = [1, 2]`
- `H.shape = (16, 16)`
- `G.shape = (10, 16)`
- `feasible_init_residual = 3.6545954574843133e-13`
- `objective_kkt = 1119.350407832252`
- `objective_reduced = 1119.3504078321841`
- `coeff_diff_norm = 2.0548403847372576e-13`
- `max_interp_error = 6.6058269965196814e-15`
- `max_Cs_minus_1_jump = 7.105427357601002e-15`
- `max_C2s_minus_2_jump = 1.7465140444983263e-11`
- `max_boundary_jet_error = 2.842170943040401e-14`
- `max_natural_bc_residual = 5.0732751333271153e-11`
- `reduced_hessian_min_eig = 5.4341103532628745`

### 分析

- 位置插值、`C^3` 连续性和 physical boundary jet 都达到了 `1e-14` 级精度。
- 人工右边界上的自然边界残差是 `5.07e-11`，对 `p^{(4)}, p^{(5)}, p^{(6)}` 这类高阶量来说已经很接近机器精度。
- 这很好地对应了文档里的命题 C 和 Case 2。

## Step 5. right boundary window 的对称验证

### 运行方式

对 `M = 5, i = 5, k = 2` 运行 local solve。

### 功能说明

验证右边界窗口与左边界窗口的对称情形：

- 右侧继承 physical boundary jet
- 左侧人工边界自动满足自然边界条件

### 实际结果

- `window = [4, 5]`
- `H.shape = (16, 16)`
- `G.shape = (10, 16)`
- `feasible_init_residual = 1.159587956055221e-13`
- `objective_kkt = 931.2997365158136`
- `objective_reduced = 931.2997365156916`
- `coeff_diff_norm = 4.1915328768134175e-12`
- `max_interp_error = 3.6637359812630166e-15`
- `max_Cs_minus_1_jump = 3.552713678800501e-15`
- `max_C2s_minus_2_jump = 2.4570567802584264e-11`
- `max_boundary_jet_error = 6.217248937900877e-15`
- `max_natural_bc_residual = 2.6545247599153835e-12`
- `reduced_hessian_min_eig = 5.036610016405353`

### 分析

- 右边界结果与左边界窗口保持同等级精度，说明 Phase 3 的左右边界处理没有非对称 bug。
- `max_natural_bc_residual` 甚至下降到 `1e-12` 量级，说明人工左边界的自然边界条件数值表现非常好。

## Step 6. 扰动连续性试验

### 运行方式

通过验证脚本读取 `table_perturbation_continuity.csv`。

### 功能说明

对固定窗口中的某个 waypoint 加入微扰 `eps`，检查局部系数和局部轨迹是否随 `eps -> 0` 连续收缩。

### 实际结果

以 canonical interior case `M = 5, i = 3, k = 2, perturb_target = 3` 为例：

- `eps = 1e-1`: `coeff_diff_norm = 2.4965610352612957e-01`, `traj_diff_sup_norm = 1.0235945446582351e-01`
- `eps = 1e-2`: `coeff_diff_norm = 2.4965610352629430e-02`, `traj_diff_sup_norm = 1.0235945446573913e-02`
- `eps = 1e-3`: `coeff_diff_norm = 2.4965610352746847e-03`, `traj_diff_sup_norm = 1.0235945446503970e-03`
- `eps = 1e-4`: `coeff_diff_norm = 2.4965610354977090e-04`, `traj_diff_sup_norm = 1.0235945446068762e-04`
- `eps = 1e-5`: `coeff_diff_norm = 2.4965610422818650e-05`, `traj_diff_sup_norm = 1.0235945439696081e-05`
- `eps = 1e-6`: `coeff_diff_norm = 2.4965610302812780e-06`, `traj_diff_sup_norm = 1.0235945404613034e-06`

同一 case 的 `cond_number = 32471.877144141665`。

### 分析

- `coeff_diff_norm` 和 `traj_diff_sup_norm` 基本与 `eps` 成线性缩放，这正是“局部映射对小扰动连续”的数值表现。
- 条件数在整个扰动序列中保持稳定，说明没有出现明显的分支切换或奇异化迹象。
- 这一步直接支撑了命题 E。

## Step 7. 多初始化唯一性验证

### 运行方式

验证脚本对每个 case 运行 `num_restarts = 6` 个 reduced-space 随机初始化。

### 功能说明

从“优化角度”验证唯一性：若局部问题确实只有一个最优解，则不同随机初始化都应收敛到相同系数。

### 实际结果

`table_uniqueness_summary.csv` 的代表性结果：

- `M=5, i=1`: `max_pairwise_coeff_diff = 2.076328291207854e-12`, `objective_std = 8.12788010264417e-12`
- `M=5, i=3`: `max_pairwise_coeff_diff = 4.609052754863275e-12`, `objective_std = 9.988979692872001e-24`
- `M=5, i=5`: `max_pairwise_coeff_diff = 1.722541524051364e-14`, `objective_std = 6.261382790496671e-13`
- `M=6, i=1`: `max_pairwise_coeff_diff = 2.1959557935276977e-12`
- `M=6, i=3`: `max_pairwise_coeff_diff = 2.459240219180731e-12`
- `M=6, i=6`: `max_pairwise_coeff_diff = 1.051010673900849e-13`

### 分析

- 所有 case 的多初始化系数差异都落在 `1e-12` 到 `1e-14` 级，已经足够说明最终解不是初始化伪象。
- 结合 reduced Hessian 的正最小特征值，这一步形成了“优化角度 + 线性代数角度”的双重唯一性证据。

## Step 8. full-window 与 Phase 1 的对照

### 运行方式

将窗口扩到覆盖全轨迹，再把 local BLOM-Strict 解与 Phase 1 global MINCO baseline 做比较。

### 功能说明

验证需求文档 10.1：当窗口覆盖全轨迹时，局部 variational problem 应退化到对应的全局问题。

### 实际结果

对 `M = 5`：

- `window_segments = [1, 2, 3, 4, 5]`
- `coeff_diff_norm = 4.4602351968766654e-12`

对 `M = 6`：

- `window_segments = [1, 2, 3, 4, 5, 6]`
- `coeff_diff_norm = 3.781329409427601e-11`

### 分析

- 两个 full-window check 的系数差异都在 `1e-11` 左右，说明 Phase 3 的局部 variational定义与 Phase 1 全局 minimum-snap 母问题在全窗口情形下是数值一致的。
- 这很好地完成了 Phase 1 / Phase 3 的联动 sanity check。

## Step 9. 批量验证脚本与结果产物

### 运行命令

```bash
python3 -m phase_3.validate_phase3_blom_strict
```

### 功能说明

批量遍历 left / interior / right 三类窗口，自动生成：

- 图片
- CSV 表格
- JSON 日志

### 实际结果

控制台输出：

```text
Phase 3 BLOM-Strict validation
cases: 6
max interpolation error: 1.701e-14
max natural BC residual: 8.983e-11
max multistart coeff diff: 4.609e-12
full-window vs Phase 1 coeff diff: 3.781e-11
results dir: phase_3/results
```

已生成的主要文件：

- Figures: `fig_window_layout_*`, `fig_local_traj_*`, `fig_continuity_jump_*`, `fig_natural_bc_residual_*`, `fig_perturbation_response_*`, `fig_uniqueness_multistart_*`
- Tables: `table_feasibility_summary.csv`, `table_uniqueness_summary.csv`, `table_perturbation_continuity.csv`
- Logs: `phase3_case_*.json`, `phase3_validation_summary.json`

### 分析

- `6` 个 case 覆盖了 `M = 5, 6` 的 left / interior / right 窗口，已经满足“多个 M、多个窗口位置”的批量验证要求。
- 批量结果中最差的插值误差仍只有 `1.70e-14`，最差自然边界残差也在 `1e-10` 量级，整体数值精度是非常好的。

## Step 10. 自动化测试

### 运行命令

```bash
python3 -m unittest phase_3.test_phase3_blom_strict
python3 -m unittest discover -s . -p 'test*.py'
```

### 功能说明

分别验证：

- Phase 3 专项回归
- 全仓回归，确认没有破坏 Phase 0 / 1 / 2

### 实际结果

Phase 3 专项测试：

```text
.....
----------------------------------------------------------------------
Ran 5 tests in 2.015s

OK
```

全仓测试：

```text
..............................
----------------------------------------------------------------------
Ran 30 tests in 2.654s

OK
```

### 分析

- Phase 3 自身的 5 个测试全部通过，覆盖了窗口构造、可行解、KKT/reduced 一致性、full-window 对照和结果导出。
- 全仓 `30` 个测试通过，说明这次实现没有破坏之前的 Phase 0 / 1 / 2 基线。

## 总结

本次实际运行表明：

- BLOM-Strict 局部 admissible problem 可以稳定构造；
- KKT 与 reduced 两条求解路径数值一致；
- reduced Hessian 在所有验证 case 上都正定；
- artificial boundary 自然边界条件成立；
- interior window 自动提升到 `C^6`；
- 局部解对小扰动连续；
- full-window 情形与 Phase 1 global MINCO 一致；
- Phase 3 所要求的 PNG / CSV / JSON 产物都已自动保存。

因此，Phase 3 所要求的“局部问题是合法数学对象”的数值验证目标已经完成。
