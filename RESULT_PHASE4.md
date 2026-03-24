# RESULT_PHASE4

本文档记录 2026-03-24 对 `phase_4/` 代码的实际运行结果。  
结构统一为“步骤 -> 功能说明 -> 实际结果 -> 分析”，所有数值均来自本次实际执行。

Phase 4 的目标不是再证明局部问题存在，而是把 Phase 3 的 BLOM-Strict 局部 variational 对象转成可检查的 analytic local solver，并验证：

- `s=2,k=2` 时，时间加权 Catmull--Rom 公式就是 exact local optimum；
- `s=4,k=2` 时，存在一个真正的 `6x6` 精确局部系统；
- 该解析系统与 Phase 3 local QP 数值吻合；
- Catmull-style heuristic 在 `s=4` 中不再是精确局部最优。

## Step 1. 静态可执行性检查

### 运行命令

```bash
python3 -m compileall phase_4
```

### 功能说明

确认 `phase_4/` 的三份主脚本、工具模块和测试文件都能被 Python 正常解析与导入。

### 实际结果

```text
Listing 'phase_4'...
Compiling 'phase_4/blom_k2_s4_numeric.py'...
Listing 'phase_4/results'...
Listing 'phase_4/results/catmull_compare'...
Listing 'phase_4/results/s2_sympy'...
Listing 'phase_4/results/s4_numeric'...
Listing 'phase_4/utils'...
Compiling 'phase_4/utils/plotting_utils.py'...
```

### 分析

- `compileall` 正常结束，说明 Phase 4 的模块结构已经完整。
- 结果目录 `s2_sympy / s4_numeric / catmull_compare` 也已经就位，后续可以直接落地中间结果和图片。

## Step 2. Script A: `s=2, k=2` 符号验证

### 运行命令

```bash
python3 -m phase_4.blom_k2_s2_sympy
```

### 功能说明

用 SymPy 显式完成 warm-up 模型的全部推导：

- one-sided natural cubic 系数
- 单段和两段局部代价
- 对 `v_i` 的一阶最优性条件
- 闭式最优速度
- 与时间加权 Catmull--Rom 公式的严格等价

### 实际结果

控制台输出：

```text
Phase 4 s=2 symbolic validation
symbolic difference: 0
random max abs error: 1.332e-15
results dir: phase_4/results/s2_sympy
```

落地文件包括：

- `s2_left_segment_coeffs.txt`
- `s2_right_segment_coeffs.txt`
- `s2_energy_expression.txt`
- `s2_optimal_velocity_expression.txt`
- `s2_catmull_equivalence.txt`
- `s2_energy_vs_v.png`
- `s2_catmull_match_demo.png`

关键指标：

- `v_exact - v_catmull` 的符号化简结果：`0`
- 随机样本最大绝对误差：`1.332e-15`

### 分析

- `symbolic difference = 0` 说明 Catmull 等价不是“数值近似成立”，而是符号上严格成立。
- 随机数值误差也在机器精度范围内，进一步验证了符号推导和数值评估是一致的。
- 这一步直接完成了 Phase 4 最重要的 warm-up 结论：`s=2` 中 Catmull 就是 exact local optimum。

## Step 3. Script B: `s=4, k=2` 精确局部解析系统

### 运行命令

```bash
python3 -m phase_4.blom_k2_s4_numeric
```

### 功能说明

构造 canonical minimum-snap 情形下的局部解析系统：

\[
A_2^{(i)}(T)x_i^{loc} = B_2^{(i)}(q,T)
\]

同时实现：

- 左/右 outer segment rank-one cost
- 中央段精确 Hermite energy
- `C4 / G4 / R4` 自动生成
- `x_i^{loc,*}` 求解
- 中央段 8 个系数的 Hermite 重构
- 与 Phase 3 local QP 的数值对照

### 实际结果

控制台输出：

```text
Phase 4 s=4 numeric validation
symmetry error: 0.000e+00
min eigenvalue: 3.139e+00
state error vs QP: 8.565e-11
coeff error vs QP: 1.292e-11
results dir: phase_4/results/s4_numeric
```

代表性样本的结构化结果：

- `A2.shape = (6, 6)`
- `B2.shape = (6,)`
- `residual_norm = 3.6450773090749405e-12`
- `min_eigenvalue = 3.13943665709458`
- `max_eigenvalue = 36677.348907628606`
- `state_error_norm = 8.564831735037706e-11`
- `coeff_error_norm = 1.2915049525142093e-11`

代表性 `x_exact`：

```text
[-0.5902006357655678, -2.0064986557018525, 3.8563672150536545,
 -0.08553706866359802, 2.813960363038388, 3.856367214933717]
```

代表性中央段系数 `c_exact`：

```text
[0.85, -0.5902006357655678, -1.0032493278509262, 0.6427278691756091,
 -3.2402567740064114e-12, -6.599293556064368e-13,
 6.28642737865448e-15, 7.35744833946228e-14]
```

Hermite 端点重构残差最大值：

- `1.851852005074761e-13`

保存的关键文件包括：

- `A2_matrix.npy`
- `B2_vector.npy`
- `x_local_opt.npy`
- `center_coeffs.npy`
- `R4_matrix.npy`
- `C4_matrix.npy`
- `G4_matrix.npy`
- `analytic_vs_qp_metrics.json`
- `A2_heatmap.png`
- `A2_eigenvalues.png`
- `analytic_vs_qp_error_bar.png`
- `center_segment_overlay.png`
- `local_energy_landscape_2d.png`

### 分析

- `A2` 已经达到数值上严格对称，且最小特征值显著大于 0，这正是“解析局部系统对称正定”的核心验证。
- `A2 x - B2` 残差在 `1e-12` 级，说明 `6x6` 解析系统求解稳定。
- analytic 与 Phase 3 local QP 的状态误差约 `8.56e-11`、中央段系数误差约 `1.29e-11`，说明 Phase 4 的 analytic core 已经成功复现了 Phase 3 的 variational ground truth。
- 中央段重构的高阶系数后四项几乎为零，说明该代表样本下的中央段结构非常规整，但仍然由一般 Hermite / Gram 矩阵流程自动得到，并不是手工特化结果。

## Step 4. Script C: Catmull heuristic 与 exact analytic 的对比

### 运行命令

```bash
python3 -m phase_4.blom_catmull_compare
```

### 功能说明

分两层比较：

- `s=2`: exact local optimum vs Catmull
- `s=4`: exact analytic local state vs 一个明确标注为 heuristic baseline 的 Catmull-like local state

并在均匀时间、非均匀时间、不同几何模式下做随机 benchmark。

### 实际结果

控制台输出：

```text
Phase 4 Catmull comparison
s2 max abs error: 6.661e-16
s4 mean coeff error: 4.766e+02
s4 mean objective gap: 8.852e+05
results dir: phase_4/results/catmull_compare
```

`s2_compare_metrics.json`：

- `max_abs_error = 6.661338147750939e-16`
- `mean_abs_error = 1.1075004725031374e-16`

`s4_compare_metrics.json`：

- `mean_coeff_error = 476.61691491027193`
- `max_coeff_error = 11206.423784583396`
- `mean_objective_gap = 885187.2100506986`
- `min_objective_gap = 0.6111533079019864`
- `max_objective_gap = 28097520.840984255`

从随机 benchmark 表中统计：

- uniform-time 平均 `s4` 系数误差：`71.18068534014795`
- nonuniform-time 平均 `s4` 系数误差：`679.3350296953339`
- uniform-time 平均 objective gap：`90300.66831407842`
- nonuniform-time 平均 objective gap：`1282630.4809190084`
- 最大时间不均匀度：`3.996446719009407`

保存的关键文件包括：

- `s2_compare_metrics.json`
- `s4_compare_metrics.json`
- `random_benchmark_table.csv`
- `s2_exact_vs_catmull_scatter.png`
- `s4_velocity_exact_vs_catmull_scatter.png`
- `s4_coeff_error_hist.png`
- `s4_objective_gap_hist.png`
- `representative_curve_overlay.png`
- `error_vs_time_imbalance.png`

### 分析

- `s=2` 的 Catmull 误差仍然处于机器精度范围，和 Script A 的符号结论完全一致。
- `s=4` 中 heuristic 与 exact 的偏差不只是“轻微不同”，而是系统性且显著的。
- 更关键的是 `objective_gap` 的最小值仍然为正 `0.611...`，说明在全部随机样本里，exact analytic 解的局部目标值都不高于 heuristic。
- 非均匀时间下的平均系数误差和平均 objective gap 都显著大于均匀时间，这正符合 Phase 4 需求中“非均匀时间是最关键对比情形”的预期。

## Step 5. Phase 4 专项测试

### 运行命令

```bash
python3 -m unittest phase_4.test_phase4_analytic
```

### 功能说明

验证四类关键性质：

- `s=2` 符号 Catmull 等价
- `s=4` 的 `A2` 对称正定
- analytic vs Phase 3 local QP 一致
- `s=4` heuristic 的确不是 exact

### 实际结果

```text
.....
----------------------------------------------------------------------
Ran 5 tests in 0.869s

OK
```

### 分析

- 5 个专项测试全部通过，说明 Phase 4 的三层验证逻辑已经形成闭环。
- 这里尤其重要的是第三和第四类测试，它们分别把 Phase 4 锚定回 Phase 3，又把 Catmull baseline 与 exact theory 明确区分开了。

## Step 6. 全仓回归测试

### 运行命令

```bash
python3 -m unittest discover -s . -p 'test*.py'
```

### 功能说明

确认 Phase 4 的新增代码没有破坏 Phase 0 / 1 / 2 / 3 的既有实现。

### 实际结果

```text
...................................
----------------------------------------------------------------------
Ran 35 tests in 3.428s

OK
```

### 分析

- 全仓 `35` 个测试通过，说明当前仓库从 Phase 0 到 Phase 4 的代码仍然一致可用。
- 这也说明 Phase 4 与之前阶段的接口耦合方式是健康的，没有偷偷破坏 Phase 3 的 local QP ground truth。

## 总结

本次实际运行表明：

- `s=2,k=2` 时，时间加权 Catmull--Rom 在符号和数值两层都与 exact local optimum 完全一致；
- `s=4,k=2` 时，可以稳定构造一个对称正定的 `6x6` 精确局部解析系统；
- 该解析系统和 Phase 3 local QP 在状态与中央段系数上达到 `1e-10` 到 `1e-11` 级一致；
- Catmull-style heuristic 在 `s=4` 中不再是 exact，且在非均匀时间下偏差更明显；
- Phase 4 要求的 symbolic / numeric / compare 三层验证都已完成，并且结果文件和图片已自动保存。

因此，Phase 4 所要求的“BLOM-Strict 到 BLOM-Analytic 的验证器与解析内核原型”已经完成。
