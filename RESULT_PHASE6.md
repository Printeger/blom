# RESULT_PHASE6

本文档记录 2026-03-24 对 `phase_6/` 代码的实际运行结果。  
结构统一为“步骤 -> 功能说明 -> 实际结果 -> 分析”，所有数值均来自本次真实执行。

Phase 6 的目标不是重新推导 BLOM，而是把论文里的 Jacobian 基础性质转成统一可复现的数值验证程序，重点覆盖：

- raw BLOM-Analytic 的 exact local support
- 解析 Jacobian 与 finite difference 的一致性
- `dc/dT` 的弱时间感知与有界性
- Scheme A / B / C 在装配后 locality 上的差异

## Step 1. 静态可执行性检查

### 运行命令

```bash
python3 -m compileall phase_6
```

### 功能说明

确认 Phase 6 主模块、examples、tests 和结果目录结构都能被 Python 正常解析与导入。

### 实际结果

```text
Listing 'phase_6'...
Listing 'phase_6/examples'...
Compiling 'phase_6/examples/__init__.py'...
Compiling 'phase_6/examples/demo_compare_all.py'...
Compiling 'phase_6/examples/demo_raw_local_support.py'...
Listing 'phase_6/results'...
Listing 'phase_6/results/phase6_fd_jacobian_check'...
Listing 'phase_6/results/phase6_fd_jacobian_check/compare'...
Listing 'phase_6/results/phase6_fd_jacobian_check/raw_scheme_C'...
Listing 'phase_6/results/phase6_fd_jacobian_check/scheme_A'...
Listing 'phase_6/results/phase6_fd_jacobian_check/scheme_B'...
```

### 分析

- `compileall` 正常结束，说明 Phase 6 的代码组织、结果目录和 examples/tests 都已成型。
- 结果目录结构与需求文档一致，已经区分 `raw_scheme_C`、`scheme_A`、`scheme_B` 和 `compare`。

## Step 2. raw local support 与解析 Jacobian 校验

### 运行命令

```bash
python3 -m phase_6.examples.demo_raw_local_support
```

### 功能说明

对 canonical interior raw map 做最直接的本地验证：

- 检查 `J_c_q` / `J_c_T` 的理论零模式
- 计算解析 Jacobian
- 用 central finite difference 对照解析 Jacobian

### 实际结果

```text
Phase 6 demo raw local support
raw J_c_q max abs error: 4.063e-08
raw J_c_T max abs error: 1.209e-07
raw q mask pass: True
```

从汇总文件读取的关键结果：

- `raw J_c_q max abs error = 4.062682501526638e-08`
- `raw J_c_T max abs error = 1.209002866708886e-07`
- `raw J_c_q mask pass = True`
- `raw J_c_T mask pass = True`

### 分析

- `10^-8` 到 `10^-7` 量级的误差对于一条含逆矩阵和时间缩放链条的 analytic-vs-FD 校验来说是非常健康的。
- `mask pass = True` 说明理论上的“只依赖局部 stencil”并不是近似现象，而是数值上可以稳定识别的 exact local support。
- 这一步已经直接支撑了 Phase 6 Proposition 1 的核心结论。

## Step 3. Scheme A 的非局部化验证

### 运行命令

```bash
python3 -m phase_6.examples.demo_scheme_A_nonlocality
```

### 功能说明

检查 shared junction states 装配后，最终 segment coefficient 的 Jacobian 是否比 raw Scheme C 更宽。

### 实际结果

```text
Phase 6 demo scheme A nonlocality
q bandwidth: 5
outside-band q max abs: 1.317e+00
```

### 分析

- Scheme A 的 `q bandwidth = 5`，明显大于 raw Scheme C 的理论局部带宽 `2`。
- `outside-band q max abs = 1.317` 远大于数值噪声，说明远场参数敏感度是真实存在的，不是 FD 误差造成的假象。
- 这正是 Phase 6 要验证的结论：Scheme A 用 shared state 换取了全局平滑性，也因此牺牲了 exact finite local support。

## Step 4. Scheme B 的传播范围验证

### 运行命令

```bash
python3 -m phase_6.examples.demo_scheme_B_nonlocality
```

### 功能说明

检查 overlapping consensus 装配后，最终 coefficient 的 Jacobian 是否比 raw 更宽、但仍比 Scheme A 更保留局部性。

### 实际结果

```text
Phase 6 demo scheme B nonlocality
q bandwidth: 3
outside-band q max abs: 8.843e+00
```

### 分析

- Scheme B 的 `q bandwidth = 3`，确实比 raw Scheme C 更宽，但又明显小于 Scheme A 的 `5`。
- `outside-band q max abs = 8.843` 说明 consensus 以后，远场影响已经真实传播到了本地段系数。
- 这与论文里的 Phase 6 叙述一致：Scheme B 比 C 更非局部，但通常仍比 A 更保留 sliding-window 风格。

## Step 5. 全部方案统一对比与随机试验

### 运行命令

```bash
python3 -m phase_6.blom_fd_jacobian_check
python3 -m phase_6.examples.demo_compare_all
```

### 功能说明

这一步运行 Phase 6 的完整主流程：

- raw analytic-vs-FD 校验
- Scheme A / B / C final coefficient Jacobian 对比
- bandwidth 比较
- `dc/dT` 统计
- 多组随机 trial 的 mask/误差/局部性验证

### 实际结果

主脚本输出：

```text
Phase 6 FD Jacobian check
raw q max abs error: 4.063e-08
raw T max abs error: 1.209e-07
scheme C q bandwidth: 2
scheme A q bandwidth: 5
scheme B q bandwidth: 3
random raw mask pass rate q: 1.00
results dir: phase_6/results/phase6_fd_jacobian_check
```

compare-all demo 输出：

```text
Phase 6 demo compare all
scheme C q bandwidth: 2
scheme A q bandwidth: 5
scheme B q bandwidth: 3
```

从 `summary_phase6_fd_check.json` 和 `random_trials_summary.json` 中读取的关键统计量：

- raw representative case
  - `RAW_Q_MAX = 4.062682501526638e-08`
  - `RAW_T_MAX = 1.209002866708886e-07`
  - `RAW_MASK_PASS_Q = True`
  - `RAW_MASK_PASS_T = True`
- representative bandwidth / outside-band values
  - `SCHEME_C_Q_BW = 2`
  - `SCHEME_A_Q_BW = 5`
  - `SCHEME_B_Q_BW = 3`
  - `SCHEME_C_OUT_Q_MAX = 0.0`
  - `SCHEME_A_OUT_Q_MAX = 1.3170582951471665`
  - `SCHEME_B_OUT_Q_MAX = 8.842513764761861`
- random trials
  - `RANDOM_RAW_Q_PASS = 1.0`
  - `RANDOM_RAW_T_PASS = 1.0`
  - `RANDOM_A_Q_BW = 7.0`
  - `RANDOM_B_Q_BW = 3.0`
  - `RANDOM_C_Q_BW = 2.0`
  - `RANDOM_MAX_DCDT = 1.7653138539427928`

生成的关键文件包括：

- `raw_scheme_C/jacobian_mask_q_raw.png`
- `raw_scheme_C/jacobian_mask_T_raw.png`
- `raw_scheme_C/jacobian_theory_vs_fd_q_raw.png`
- `raw_scheme_C/jacobian_theory_vs_fd_T_raw.png`
- `compare/jacobian_error_analytic_vs_fd_q.png`
- `compare/jacobian_error_analytic_vs_fd_T.png`
- `compare/fd_stepsize_sweep_q.png`
- `compare/fd_stepsize_sweep_T.png`
- `compare/jacobian_bandwidth_compare.png`
- `compare/dc_dT_bound_statistics.png`
- `compare/jacobian_error_boxplot_random_trials.png`
- `raw_scheme_C/jacobian_sparsity_stats_raw.csv`
- `scheme_A/jacobian_sparsity_stats_scheme_A.csv`
- `scheme_B/jacobian_sparsity_stats_scheme_B.csv`
- `compare/jacobian_error_stats.csv`
- `compare/summary_phase6_fd_check.json`

### 分析

- raw Scheme C 的 `q` 与 `T` mask 在随机试验中通过率都是 `1.00`，说明理论零模式在多组样本下稳定成立。
- representative case 中 C/A/B 的 `q bandwidth = 2 / 5 / 3`，很清晰地给出了 Phase 6 最关键的数值结论：C 最局部，B 次之，A 最宽。
- 在随机试验中，Scheme A 的 mean `q bandwidth = 7.0`，已经接近“全局敏感”的状态；B 维持在 `3.0`，而 C 稳定在 `2.0`。
- `RANDOM_MAX_DCDT = 1.765...` 也说明在当前 compact admissible 采样范围内，`dc/dT` 没有出现爆炸性增长，这为 Proposition 4 的“有界性”提供了直接数值支撑。

## Step 6. Phase 6 专项测试

### 运行命令

```bash
python3 -m unittest phase_6.test_blom_fd_jacobian_check
```

### 功能说明

验证以下关键性质：

- 理论 mask 构造正确
- raw analytic-vs-FD 的 `J_c_q` / `J_c_T` 一致
- Scheme C 的带宽最窄
- 结果文件能正确生成

### 实际结果

```text
.....
----------------------------------------------------------------------
Ran 5 tests in 7.278s

OK
```

### 分析

- 5 个专项测试全部通过，说明 Phase 6 的 Jacobian 验证底座、scheme locality 对比和 artifact 导出已经闭环。
- 这里最关键的是 locality ordering 测试，它把“C 最局部、A/B 更宽”的核心结论固定成了可回归的自动检查。

## Step 7. 全仓回归测试

### 运行命令

```bash
python3 -m unittest discover -s . -p 'test*.py'
```

### 功能说明

确认 Phase 6 的新增实现没有破坏 Phase 0 到 Phase 5 的既有功能。

### 实际结果

```text
.............................................
----------------------------------------------------------------------
Ran 45 tests in 11.475s

OK
```

### 分析

- 全仓 `45` 个测试全部通过，说明当前仓库从 Phase 0 到 Phase 6 仍是一套一致可用的实现。
- 这也说明 Phase 6 对 Phase 4 raw map、Phase 5 assembly 和 Phase 2 bandwidth 思路的复用是稳定的，没有引入回归。

## 总结

本次实际运行表明：

- raw BLOM-Analytic 的 Jacobian 稀疏模式与论文理论一致；
- raw `J_c_q` 与 `J_c_T` 的解析公式都能被 finite difference 稳定复现；
- raw Scheme C 保持最强的 exact finite local support；
- Scheme A / B 在获得更强装配连续性的同时，的确引入了更宽的 Jacobian；
- Scheme B 的 locality 介于 A 和 C 之间，符合“比 C 更宽、比 A 更局部”的预期；
- `dc/dT` 在当前紧致 admissible 采样上保持有界，支撑了弱时间感知命题。

因此，Phase 6 所要求的“FD Jacobian 基础性质验证器”已经完成。
