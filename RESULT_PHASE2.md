# RESULT_PHASE2

本文档记录 2026-03-24 对 `phase_2/` 代码的实际运行结果。  
结构统一为“步骤 -> 功能说明 -> 实际结果 -> 分析”，所有数值均来自本次实际执行。

Phase 2 的重点不再是重新发明求解器，而是把 Phase 1 的 `minco_scalar_baseline.py` 作为固定后端，验证需求文档里的必要性分析链路：

- 解析 Jacobian 提取
- 有限差分一致性检查
- 影响传播统计
- effective bandwidth 分析
- 图、表、日志产物导出

## Step 1. 静态可执行性检查

### 运行命令

```bash
python3 -m compileall phase_2
```

### 功能说明

确认 `phase_2/` 下的验证模块、绘图模块和示例脚本都能被 Python 正常解析和导入。

### 实际结果

```text
Listing 'phase_2'...
Compiling 'phase_2/__init__.py'...
Listing 'phase_2/examples'...
Compiling 'phase_2/examples/__init__.py'...
Compiling 'phase_2/examples/demo_phase2_fd_check.py'...
Compiling 'phase_2/examples/demo_phase2_nonuniform_time.py'...
Compiling 'phase_2/examples/demo_phase2_scaling_M.py'...
Compiling 'phase_2/examples/demo_phase2_uniform_time.py'...
Compiling 'phase_2/phase2_plotting.py'...
Compiling 'phase_2/phase2_validation.py'...
Listing 'phase_2/results'...
Listing 'phase_2/results/figures'...
Listing 'phase_2/results/logs'...
Listing 'phase_2/results/tables'...
```

### 分析

- `compileall` 成功，说明 Phase 2 代码在语法和导入层面是完整的。
- 这一步也证明 `results/` 目录结构已经就位，可供后续实验直接写入产物。

## Step 2. Uniform-time 验证套件

### 运行命令

```bash
python3 -m phase_2.examples.demo_phase2_uniform_time
```

### 功能说明

在一个平滑、均匀时间分配的 `M = 8` 案例上，执行完整验证流程：

- 从 Phase 1 后端求解系数
- 构造 `S_q`
- 计算解析 Jacobian `J_q = A(T)^{-1} S_q`
- 用有限差分构造 `J_fd`
- 比较 `J_q` 与 `J_fd`
- 统计 block influence 和 effective bandwidth
- 输出图、表、日志

### 实际结果

命令行输出：

```text
Phase 2 uniform-time demo
M: 8
fro error (exact vs fd): 1.179e-07
max effective bandwidth: 7
results dir: phase_2/results
```

结构化指标：

- `q = [0.15, 0.565685424949238, 0.65, 0.5656854249492381, 0.1500000000000001, -0.565685424949238, -0.9500000000000001, -0.5656854249492382, 0.1499999999999998]`
- `T = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]`
- `system_residual_norm = 1.0954794401904095e-13`
- `fro_error = 1.1794703841772456e-07`
- `max_abs_error = 3.09126098052559e-08`
- `relative_error = 6.04075768968253e-09`
- `max_effective_bandwidth = 7`
- `mean_effective_bandwidth = 5.285714285714286`
- `far_nonzero_ratio = 1.0`
- `max_far_influence = 4.818996211622124`
- `target_profile_max_influence = 2.2909009070209168`
- `target_profile_farthest_nonzero_distance = 5`

### 分析

- 解析 Jacobian 和有限差分 Jacobian 的相对误差只有 `6.04e-09`，说明解析推导和实现是一致的。
- `far_nonzero_ratio = 1.0` 表明远离目标 waypoint 的段上仍然存在明显影响，正好支持“非局部支撑”这一 Phase 2 核心结论。
- `max_effective_bandwidth = 7` 与 `M = 8` 的设置对应，说明在当前阈值 `1e-10` 下，影响并没有表现出真正的紧带状截断。

## Step 3. Nonuniform-time 验证套件

### 运行命令

```bash
python3 -m phase_2.examples.demo_phase2_nonuniform_time
```

### 功能说明

把同样的验证流程放到不均匀时间分配上，检查结论是否只在“均匀时间”这一理想条件下成立。

### 实际结果

命令行输出：

```text
Phase 2 nonuniform-time demo
M: 8
relative Jacobian error: 2.397e-08
far-field ratio: 1.000e+00
results dir: phase_2/results
```

结构化指标：

- `q = [0.0, 0.2437736638035451, -0.8319872849923965, 0.6003609566451659, 0.7524517731129712, -1.5608281509230693, -1.0417436054898546, 0.1022723225338283, -0.2529940738748658]`
- `T = [0.6281136326755459, 0.9503859378955671, 0.8707980242325812, 1.4267649888486018, 1.1438651200806644, 1.32276161327083, 0.9434141988273311, 0.7272387217847769]`
- `system_residual_norm = 1.8631960027791706e-12`
- `fro_error = 2.1779704695144896e-06`
- `max_abs_error = 1.1174495649912686e-06`
- `relative_error = 2.3967005153260997e-08`
- `max_effective_bandwidth = 7`
- `mean_effective_bandwidth = 5.285714285714286`
- `far_nonzero_ratio = 1.0`
- `max_far_influence = 16.27751785310837`
- `target_profile_max_influence = 3.034016929857525`
- `target_profile_farthest_nonzero_distance = 5`

### 分析

- 即使时间分配变成非均匀，解析 Jacobian 和有限差分仍然高度一致，相对误差仍在 `1e-8` 级。
- `far_nonzero_ratio = 1.0` 依旧成立，说明“影响不是局部支撑”的结论并不依赖均匀时间这一特殊情形。
- 非均匀时间下的远场 influence 明显更大，这和系统条件数变化、段尺度不一致带来的传播放大是吻合的。

## Step 4. 有限差分敏感性检查

### 运行命令

```bash
python3 -m phase_2.examples.demo_phase2_fd_check
```

### 功能说明

使用更敏感的 `eps = 1e-7` 和另一组 `M = 6` 非均匀时间样例，专门观察解析 Jacobian 与有限差分之间的数值吻合程度。

### 实际结果

命令行输出：

```text
Phase 2 finite-difference check demo
fro error: 2.168e-05
max abs error: 1.069e-05
relative error: 1.377e-07
results dir: phase_2/results
```

结构化指标：

- `q = [0.0, 0.0009841226859860594, 0.23899643000677592, -0.21931028428977406, -0.7124734710058194, -0.36373662813737806, -0.79331724399717]`
- `T = [0.5052653045655747, 1.3212284183827663, 1.2970694287520463, 0.9679349528437208, 0.8030324268193135, 0.7784256121007733]`
- `system_residual_norm = 1.077606224488413e-13`
- `fro_error = 2.1680506966869264e-05`
- `max_abs_error = 1.0692300695058776e-05`
- `relative_error = 1.3773347440725413e-07`
- `max_effective_bandwidth = 5`
- `mean_effective_bandwidth = 3.8`
- `far_nonzero_ratio = 1.0`
- `max_far_influence = 21.896499703478703`
- `target_profile_max_influence = 5.483298898481389`
- `target_profile_farthest_nonzero_distance = 4`

### 分析

- 绝对误差比前两个案例大，但相对误差依然只有 `1.377e-07`，仍在可接受范围内。
- 这说明有限差分本身在更小 `eps`、更敏感案例下会受数值舍入影响，但解析 Jacobian 主体仍然可信。
- 该案例进一步说明“远场影响非零”不是偶发现象，而是系统结构本身决定的结果。

## Step 5. Scaling-M 实验

### 运行命令

```bash
python3 -m phase_2.examples.demo_phase2_scaling_M
```

### 功能说明

研究随着分段数 `M` 增大，effective bandwidth 是否会保持局部支撑，还是会继续扩展。

### 实际结果

命令行输出：

```text
Phase 2 scaling-M demo
M=4: max_bw=3, mean_bw=2.33, far_ratio=1.000e+00
M=8: max_bw=7, mean_bw=5.29, far_ratio=1.000e+00
M=16: max_bw=15, mean_bw=11.27, far_ratio=1.000e+00
M=32: max_bw=31, mean_bw=23.26, far_ratio=1.000e+00
results dir: phase_2/results
```

结构化记录：

- `M = 4`: `max_bw = 3`, `mean_bw = 2.3333333333333335`, `relative_error = 1.9254070613621197e-08`
- `M = 8`: `max_bw = 7`, `mean_bw = 5.285714285714286`, `relative_error = 6.04075768968253e-09`
- `M = 16`: `max_bw = 15`, `mean_bw = 11.266666666666667`, `relative_error = 3.347853888907687e-09`
- `M = 32`: `max_bw = 31`, `mean_bw = 23.258064516129032`, `relative_error = 2.810124381954938e-09`

### 分析

- 最大 effective bandwidth 基本呈现 `M - 1` 的增长关系：`3, 7, 15, 31`。
- 这说明随着分段数增大，影响范围并没有收缩成固定宽度的局部带，而是继续扩展。
- `far_ratio` 始终为 `1.0`，进一步支持需求文档中“必要性意义下不满足严格局部支撑”的论点。

## Step 6. 结果产物检查

### 功能说明

确认 Phase 2 运行后，需求里要求的图、表和日志确实被生成。

### 实际结果

生成的图像文件：

- `phase_2/results/figures/A_sparsity_M6_casefd_check.png`
- `phase_2/results/figures/A_sparsity_M8_casenonuniform_time.png`
- `phase_2/results/figures/A_sparsity_M8_caseuniform_time.png`
- `phase_2/results/figures/Jq_heatmap_M6_casefd_check.png`
- `phase_2/results/figures/Jq_heatmap_M8_casenonuniform_time.png`
- `phase_2/results/figures/Jq_heatmap_M8_caseuniform_time.png`
- `phase_2/results/figures/block_influence_M6_casefd_check.png`
- `phase_2/results/figures/block_influence_M8_casenonuniform_time.png`
- `phase_2/results/figures/block_influence_M8_caseuniform_time.png`
- `phase_2/results/figures/effective_bandwidth_vs_M_caseuniform_time.png`
- `phase_2/results/figures/influence_profile_q2_M6_casefd_check.png`
- `phase_2/results/figures/influence_profile_q3_M8_casenonuniform_time.png`
- `phase_2/results/figures/influence_profile_q3_M8_caseuniform_time.png`
- `phase_2/results/figures/jacobian_fd_compare_M6_casefd_check.png`
- `phase_2/results/figures/jacobian_fd_compare_M8_casenonuniform_time.png`
- `phase_2/results/figures/jacobian_fd_compare_M8_caseuniform_time.png`

生成的表格文件：

- `phase_2/results/tables/effective_bandwidth_summary.csv`
- `phase_2/results/tables/far_field_influence_summary.csv`
- `phase_2/results/tables/jacobian_error_summary.csv`

生成的日志文件：

- `phase_2/results/logs/effective_bandwidth_scaling_caseuniform_time.json`
- `phase_2/results/logs/phase2_validation_M6_casefd_check.json`
- `phase_2/results/logs/phase2_validation_M8_casenonuniform_time.json`
- `phase_2/results/logs/phase2_validation_M8_caseuniform_time.json`

### 分析

- 图、表、日志均已落地，说明 Phase 2 不只是“能算”，而是已经具备结果沉淀和可视化能力。
- 当前实现中，部分共享 CSV 会按最近一次运行覆盖；但每个案例对应的 JSON 日志会保留下来，因此复盘和追溯仍然是完整的。

## Step 7. 单元测试

### 运行命令

```bash
python3 -m unittest phase_2.test_phase2_validation
```

### 功能说明

验证 Jacobian、selector、bandwidth 统计和结果导出逻辑是否能稳定回归。

### 实际结果

```text
.....
----------------------------------------------------------------------
Ran 5 tests in 0.710s

OK
```

### 分析

- `5` 个测试全部通过，说明 Phase 2 的主流程已经被自动化回归覆盖。
- 结合 demo、图表和日志产物，可以认为这一版 Phase 2 已经具备可复现实验条件。

## 总结

本次实际运行表明：

- Phase 2 代码可执行、可测试、可生成图表与日志；
- 解析 Jacobian 与有限差分 Jacobian 在多个案例下保持高一致性；
- 远场 influence 在均匀和非均匀时间下都明显存在；
- effective bandwidth 会随 `M` 增长扩展，而不是收缩到固定局部宽度。

因此，这一版 `phase_2/` 已经完成了需求文档要求的核心验证任务，并给出了支持“非局部支撑”结论的可复现实验结果。
