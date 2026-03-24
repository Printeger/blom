# RESULT_PHASE5

本文档记录 2026-03-24 对 `phase_5/` 代码的实际运行结果。  
结构统一为“步骤 -> 功能说明 -> 实际结果 -> 分析”，所有数值均来自本次真实执行。

Phase 5 的目标不是重新推导局部 BLOM，而是验证三种全局装配机制在 interior knot 上的导数 jump 表现：

- Scheme A: shared junction states
- Scheme B: overlapping consensus
- Scheme C: central-segment extraction

核心诊断量是：

\[
\llbracket p^{(\ell)} \rrbracket(t_i)
=
p_i^{(\ell)}(T_i)-p_{i+1}^{(\ell)}(0),
\qquad \ell=0,\dots,2s-2.
\]

在 canonical case `s=4, k=2` 下，理论预期是：

- Scheme A / B 应该把 `order = 0,1,2,3` 的 jump 压到数值零；
- Scheme C 只保证 `order = 0` 为零，`order = 1,2,3` 需要显式诊断。

## Step 1. 静态可执行性检查

### 运行命令

```bash
python3 -m compileall phase_5
```

### 功能说明

确认 Phase 5 主模块、examples、tests 和结果目录结构都能被 Python 正常解析。

### 实际结果

```text
Listing 'phase_5'...
Listing 'phase_5/examples'...
Compiling 'phase_5/examples/demo_scheme_A.py'...
Compiling 'phase_5/examples/demo_scheme_B.py'...
Compiling 'phase_5/examples/demo_scheme_C.py'...
Listing 'phase_5/results'...
Listing 'phase_5/results/phase5_boundary_jump_check'...
Listing 'phase_5/results/phase5_boundary_jump_check/scheme_A'...
Listing 'phase_5/results/phase5_boundary_jump_check/scheme_B'...
Listing 'phase_5/results/phase5_boundary_jump_check/scheme_C'...
Compiling 'phase_5/test_blom_boundary_jump_check.py'...
```

### 分析

- `compileall` 正常结束，说明 Phase 5 的模块骨架、example 入口和测试文件都完整可导入。
- 结果目录 `phase_5/results/phase5_boundary_jump_check/` 以及 `scheme_A/B/C` 子目录都已就位，后续图表和统计文件可以稳定落地。

## Step 2. Scheme A: shared junction states

### 运行命令

```bash
python3 -m phase_5.examples.demo_scheme_A
```

### 功能说明

运行 Scheme A 的最小 demo：

- 用全局共享的 interior knot state `eta`
- 组装全局二次目标并求解共享状态
- 用统一 Hermite 重构所有段
- 统计 lower-order jump 与系统条件数

### 实际结果

```text
Phase 5 demo scheme A
max lower-order jump: 4.440e-12
hessian condition number: 1.907e+04
```

### 分析

- `4.440e-12` 的 lower-order jump 已经显著低于 `1e-10` 容忍度，满足 Phase 5 对 Scheme A 的核心验证目标。
- 全局 Hessian 条件数约 `1.907e+04`，说明共享状态系统在代表样本下可稳定求解，没有明显病态。
- 这一步直接验证了：只要左右段共享同一个 `eta_i`，`order = 0,1,2,3` 的 jump 会数值归零。

## Step 3. Scheme B: overlapping consensus

### 运行命令

```bash
python3 -m phase_5.examples.demo_scheme_B
```

### 功能说明

运行 Scheme B 的最小 demo：

- 先独立求解所有 local windows
- 对每个 interior knot 收集局部预测 `\hat eta_i^(r)`
- 做 closed-form weighted consensus
- 用 consensus 后的共享 knot-state 重构全局轨迹

### 实际结果

```text
Phase 5 demo scheme B
max lower-order jump: 2.551e-12
max pre-consensus dispersion: 2.005e+01
```

### 分析

- consensus 后的 lower-order jump 约 `2.551e-12`，同样处于机器精度层面，说明 Scheme B 成功恢复了全局 `C^{s-1}`。
- 更关键的是 consensus 前的局部预测分散度约 `2.005e+01`，说明 raw local windows 对同一个 knot 的预测差异可以非常明显。
- 这正体现了 Scheme B 的价值：它保留了 sliding-window 的独立求解，但通过一个显式 closed-form projection 把 lower-order continuity 恢复回来。

## Step 4. Scheme C: central-segment extraction

### 运行命令

```bash
python3 -m phase_5.examples.demo_scheme_C
```

### 功能说明

运行 Scheme C 的最小 demo：

- 对每个 segment 独立求解其中心窗口
- 直接抽取中心段系数作为全局装配结果
- 不引入共享变量，也不做 consensus
- 直接统计 position 和 lower-order derivative jump

### 实际结果

```text
Phase 5 demo scheme C
position jump max: 3.997e-15
velocity/accel/jerk jump max: 3.943e+01
```

### 分析

- `position jump max = 3.997e-15` 说明 Scheme C 的 `C^0` 性质确实严格成立。
- 但 `velocity/accel/jerk jump max = 3.943e+01` 也非常直接地说明：raw central extraction 一般并不会自动达到 `C^1` 到 `C^3`。
- 这一步和理论文档完全一致：Scheme C 是最“局部映射化”的装配方式，但更高阶连续性必须被视为诊断量，而不是自动性质。

## Step 5. 三方案统一比较与随机试验

### 运行命令

```bash
python3 -m phase_5.blom_boundary_jump_check
python3 -m phase_5.examples.demo_compare_all
```

### 功能说明

这一步会完整运行 Phase 5 的主流程：

- 对同一个代表样本同时运行 Scheme A / B / C
- 保存每个方案的 jump heatmap、max-bar、JSON、CSV
- 生成 lower-order jump 对比图、Scheme C eta mismatch 图、Scheme B consensus 改善图
- 再做一组随机 trial，统计各方案稳定性与 jump 分布

### 实际结果

主脚本输出：

```text
Phase 5 boundary jump check
scheme A: lower max=6.580e-12, higher max=1.123e-09
scheme B: lower max=1.014e-11, higher max=2.529e+05
scheme C: lower max=3.952e+01, higher max=5.248e+02
random success rates: A=1.00, B=1.00, C=1.00
results dir: phase_5/results/phase5_boundary_jump_check
```

compare-all demo 输出：

```text
Phase 5 demo compare all
scheme A: lower=6.580e-12, higher=1.123e-09
scheme B: lower=1.014e-11, higher=2.529e+05
scheme C: lower=3.952e+01, higher=5.248e+02
```

从保存的摘要文件中读取的关键统计量：

- Scheme A
  - `lower_order_max = 6.579625733138528e-12`
  - `higher_order_max = 1.1229985830141231e-09`
- Scheme B
  - `lower_order_max = 1.014299755297543e-11`
  - `higher_order_max = 252912.8821815843`
  - `max pre-consensus dispersion = 20.081789241754045`
- Scheme C
  - `lower_order_max = 39.51770066270133`
  - `higher_order_max = 524.7551686569486`
  - `max ||eta^- - eta^+||_2 = 40.16357848350809`

随机试验聚合统计：

- Scheme A
  - `success_rate = 1.0`
  - `mean_lower_order_max = 1.3070283281605791e-11`
  - `mean_higher_order_max = 1.19944772632626e-09`
- Scheme B
  - `success_rate = 1.0`
  - `mean_lower_order_max = 1.3307801619936535e-11`
  - `mean_higher_order_max = 225725.36478655436`
- Scheme C
  - `success_rate = 1.0`
  - `mean_lower_order_max = 19.764570017017107`
  - `mean_higher_order_max = 1197.9427582243081`

生成的关键文件包括：

- `scheme_A/jump_heatmap_scheme_A.png`
- `scheme_A/jump_maxbar_scheme_A.png`
- `scheme_A/jump_stats_scheme_A.csv`
- `scheme_B/jump_heatmap_scheme_B.png`
- `scheme_B/jump_maxbar_scheme_B.png`
- `scheme_B/scheme_B_consensus_improvement.png`
- `scheme_C/jump_heatmap_scheme_C.png`
- `scheme_C/jump_maxbar_scheme_C.png`
- `scheme_C/scheme_C_eta_mismatch.png`
- `jump_lower_orders_compare.png`
- `jump_boxplot_random_trials.png`
- `scheme_comparison_summary.csv`
- `phase5_interpretation_summary.md`

### 分析

- Scheme A 与 Scheme B 在 `order = 0..3` 上都稳定达到 `1e-11` 到 `1e-12` 级 jump，完全符合理论预期。
- Scheme C 的 `order 0` jump 为机器精度零，但 `order = 1..3` 的 jump 仍显著非零，这正是 Phase 5 最核心的诊断结论。
- Scheme B 的高阶 jump 很大并不是 bug，而是理论上就“不保证 `C^6`”的直接数值体现。Phase 5 文档只要求它恢复 `C^{s-1}`，也就是 `C^3`。
- 随机试验中三种方案成功率都为 `1.0`，说明当前实现是稳定的；但 lower-order smoothness 上只有 A/B 始终维持在数值零附近，C 则系统性偏离。

## Step 6. Phase 5 专项测试

### 运行命令

```bash
python3 -m unittest phase_5.test_blom_boundary_jump_check
```

### 功能说明

验证以下关键性质：

- `compute_jumps` / `summarize_jumps` 基础正确性
- Scheme A 达到 lower-order continuity
- Scheme B 达到 lower-order continuity 且 consensus 真的降低 mismatch
- Scheme C 严格 `C^0` 但一般不 `C^3`
- 结果文件能正确生成

### 实际结果

```text
.....
----------------------------------------------------------------------
Ran 5 tests in 1.161s

OK
```

### 分析

- 5 个专项测试全部通过，说明 Phase 5 的公共 jump 接口、三种装配方案和结果导出链路已经闭环。
- 这里最关键的是 Scheme B 和 Scheme C 的对照测试：前者验证“projection 可以恢复 smoothness”，后者验证“raw extraction 只能保证 `C^0`”。

## Step 7. 全仓回归测试

### 运行命令

```bash
python3 -m unittest discover -s . -p 'test*.py'
```

### 功能说明

确认 Phase 5 的新增代码没有破坏 Phase 0 到 Phase 4 的既有实现。

### 实际结果

```text
........................................
----------------------------------------------------------------------
Ran 40 tests in 4.452s

OK
```

### 分析

- 全仓 `40` 个测试全部通过，说明当前仓库从 Phase 0 到 Phase 5 仍然是一套一致可用的实现。
- 这也说明 Phase 5 与前面 phases 的接口复用方式是健康的：Phase 3 的 local QP 和 Phase 4 的 Hermite/core utilities 被复用后没有引入回归。

## 总结

本次实际运行表明：

- Scheme A 和 Scheme B 都可以稳定实现全局 `C^{s-1}` 装配；
- Scheme C 只严格保证 `C^0`，更高阶连续性必须显式诊断；
- Scheme B 在“保持 local-window independence”与“恢复全局 smoothness”之间给出了最平衡的方案；
- 当前实现已经生成了完整的 heatmap、bar chart、comparison plot、random-trial boxplot、CSV、JSON 和 interpretation summary；
- Phase 5 所要求的统一 boundary-jump checker 已经完成，并且结果足以直接用于论文展示与方案比较。

因此，Phase 5 要求的“全局装配机制统一验证器”已经完成。
