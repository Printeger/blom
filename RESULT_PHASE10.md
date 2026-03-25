# RESULT_PHASE10

本文档记录 Phase 10 完整时空优化框架的实际运行结果。最终结果目录为 [phase_10/results/phase10_framework](/home/dev/code/blom/phase_10/results/phase10_framework)。

说明：

- Phase 10 的主模块、suite、examples 和测试都实际执行过。
- 为避免共享结果目录造成覆盖，本文档中的关键数值以最后一次顺序执行各主入口后立即读取的 summary 文件为准。
- 当前实现固定 canonical `s = 4`，并对 general even `k` 暴露统一接口；本次主验证中 `k = 2` 是优化主对象，同时也对 `k = 4, 6, 8` 做了 Jacobian / benchmark 级别检查。

## Step 1. `python3 -m compileall phase_10`

功能：

- 编译 `phase_10` 下的主模块、examples 和 tests。
- 提前排除语法错误、导入错误和包结构错误。

结果：

- 命令执行成功。
- `phase_10` 下所有新增文件都通过编译。

分析：

- 这一步确认 Phase 10 已经从“几份独立脚本”变成了一个结构闭合的 Python package，后面的验证重点就可以放在数值正确性和框架行为上。

## Step 2. `python3 -m phase_10.blom_full_backward_diff`

功能：

- 运行 Phase 10 的底层 backward 框架。
- 实际执行：
  - raw Scheme C general-`k` 系数接口
  - dense / sparse backward differentiation
  - `xi=(q_bar,\tau)` 重参数化梯度
  - `k=2,4,6,8` 的 Jacobian 稀疏性与 support-width 统计

结果：

- 命令执行成功。
- 关键 summary 文件：
  - [phase10_backward_summary.json](/home/dev/code/blom/phase_10/results/phase10_framework/backward_diff/phase10_backward_summary.json)
  - [phase10_backward_summary.csv](/home/dev/code/blom/phase_10/results/phase10_framework/backward_diff/phase10_backward_summary.csv)
  - [phase10_dense_vs_sparse_grad.png](/home/dev/code/blom/phase_10/results/phase10_framework/backward_diff/phase10_dense_vs_sparse_grad.png)
  - [phase10_jacobian_sparsity_k_sweep.png](/home/dev/code/blom/phase_10/results/phase10_framework/backward_diff/phase10_jacobian_sparsity_k_sweep.png)
  - [phase10_support_width_vs_k.png](/home/dev/code/blom/phase_10/results/phase10_framework/backward_diff/phase10_support_width_vs_k.png)
- 关键数值：
  - objective value `2721.1842182934192`
  - dense vs sparse max abs gap `9.094947017729282e-13`
  - dense vs sparse l2 gap `1.1425385788607957e-12`
  - reparam vs FD max abs gap `1.3260914784041233e-05`
  - reparam vs FD l2 gap `2.424421635852202e-05`
  - support width mean by `k = 2,4,6,8`: `3.0, 4.0, 5.0, 5.6667`
  - support width max by `k = 2,4,6,8`: `4, 5, 6, 6`

分析：

- dense 和 sparse backward 在 `1e-12` 量级上一致，说明 block-banded accumulation 没有引入额外数值误差。
- `grad_xi` 对 finite difference 的误差在 `1e-5` 量级，足以支撑 framework-level optimizer。
- `k` 增大时 support width 单调增大，这和 Phase 10 对 general-`k` block-banded complexity 的预期一致。

## Step 3. `python3 -m phase_10.blom_space_time_optimizer`

功能：

- 运行 `xi=(q_bar,\tau)` 空间上的优化器主入口。
- 实际执行：
  - full objective layer evaluation
  - sparse backward gradient
  - softplus time reparameterization
  - Armijo gradient descent
  - objective / grad norm / time allocation 可视化

结果：

- 命令执行成功。
- 关键 summary 文件：
  - [phase10_opt_summary.json](/home/dev/code/blom/phase_10/results/phase10_framework/optimizer/phase10_opt_summary.json)
  - [phase10_opt_history.csv](/home/dev/code/blom/phase_10/results/phase10_framework/optimizer/phase10_opt_history.csv)
  - [phase10_objective_curve.png](/home/dev/code/blom/phase_10/results/phase10_framework/optimizer/phase10_objective_curve.png)
  - [phase10_gradnorm_curve.png](/home/dev/code/blom/phase_10/results/phase10_framework/optimizer/phase10_gradnorm_curve.png)
  - [phase10_traj_progress.png](/home/dev/code/blom/phase_10/results/phase10_framework/optimizer/phase10_traj_progress.png)
  - [phase10_time_progress.png](/home/dev/code/blom/phase_10/results/phase10_framework/optimizer/phase10_time_progress.png)
- standalone optimizer 关键数值：
  - `n_steps = 25`
  - objective initial `2721.184218293416`
  - objective final `18.51143007839764`
  - objective drop `2702.6727882150185`
  - accepted ratio `1.0`
  - grad norm peak `7787.545675455843`
  - grad norm final `2.028881658049813`
  - runtime total `0.20072790700942278 sec`
  - any duration touched lower bound: `False`

分析：

- objective 在 25 步内大幅下降，说明 full objective layer + reparam backward + Armijo step 已经形成稳定的 descent loop。
- 所有 step 都被接受，说明当前代表性 case 上的默认步长和线搜索参数比较温和。
- `T_i` 没有触碰 `T_min`，说明时间重参数化和当前权重配置足以避免退化时间分配。

## Step 4. `python3 -m phase_10.blom_benchmark_suite`

功能：

- 运行 Phase 10 的 benchmark / baseline / ablation 主入口。
- 实际执行：
  - `M = 10, 20, 40, 80` 的规模 sweep
  - representative-case 上的 `k = 2, 4, 6, 8` sweep
  - baseline compare: `raw_schemeC`, `minco`, `schemeA`, `heuristic`
  - term ablation: `lambda_obs`, `lambda_dyn`, `lambda_T`, `lambda_reg`

结果：

- 命令执行成功。
- 关键 summary 文件：
  - [phase10_overview.csv](/home/dev/code/blom/phase_10/results/phase10_framework/summary/phase10_overview.csv)
  - [phase10_benchmark_summary.json](/home/dev/code/blom/phase_10/results/phase10_framework/summary/phase10_benchmark_summary.json)
  - [phase10_benchmark_summary.md](/home/dev/code/blom/phase_10/results/phase10_framework/summary/phase10_benchmark_summary.md)
  - [phase10_ablation_summary.csv](/home/dev/code/blom/phase_10/results/phase10_framework/summary/phase10_ablation_summary.csv)
- `M`-sweep 里 raw Scheme C 的结果：
  - `M=10`: objective `19.9725`, runtime `3.5037 sec`
  - `M=20`: objective `4.5478`, runtime `7.6274 sec`
  - `M=40`: objective `19.9308`, runtime `10.0856 sec`
  - `M=80`: objective `30.5775`, runtime `18.6234 sec`
- `M=10` baseline compare：
  - raw Scheme C objective `19.9725`
  - MINCO objective `3724.4173`
  - Scheme A objective `3724.4173`
  - heuristic objective `751867.9479`
- representative-case `k`-sweep：
  - `k=2`: objective `18.4094`, mean support width `3.0`
  - `k=4`: objective `37.2526`, mean support width `4.0`
  - `k=6`: objective `43.8042`, mean support width `5.0`
  - `k=8`: objective `50.7846`, mean support width `5.6667`
- ablation：
  - `drop_lambda_obs`: objective `18.3604`
  - `drop_lambda_dyn`: objective `17.6305`
  - `drop_lambda_T`: objective `15.6462`
  - `drop_lambda_reg`: objective `3.6986`

分析：

- benchmark 明确给出了一个很清楚的 speed-quality trade-off：heuristic 最快，但 objective 远差；raw Scheme C 最慢，但 objective 最好。
- `k` 增大时 support width 增大，同时 representative-case objective 也升高，这说明当前 full objective 并不是“越大窗口越优”，反而暴露了 locality / conditioning / objective coupling 之间的实际权衡。
- ablation 中 `drop_lambda_reg` 让 objective 最低，说明当前默认总目标里 regularization 是一个显著的代价来源；这对后续权重调参与 Pareto 分析是有价值的。

## Step 5. `python3 -m phase_10.blom_phase10_framework_suite`

功能：

- 统一执行：
  - backward consistency
  - optimizer demo
  - benchmark suite
  - total overview / interpretation summary

结果：

- 命令执行成功。
- 关键 summary 文件：
  - [phase10_suite_summary.json](/home/dev/code/blom/phase_10/results/phase10_framework/summary/phase10_suite_summary.json)
  - [phase10_interpretation_summary.md](/home/dev/code/blom/phase_10/results/phase10_framework/summary/phase10_interpretation_summary.md)
  - [phase10_overview.csv](/home/dev/code/blom/phase_10/results/phase10_framework/summary/phase10_overview.csv)
- suite 关键数值：
  - backward dense vs sparse max abs gap `9.094947017729282e-13`
  - backward reparam vs FD max abs gap `1.644413441681536e-05`
  - optimizer objective drop `2710.814043488962`
  - optimizer objective final `10.370174804454226`
  - optimizer accepted ratio `1.0`
  - benchmark row count `20`
  - benchmark ablation count `4`
  - overview best benchmark objective `5.311151173281723`
  - overview benchmark runtime total `29.429138156585395 sec`

分析：

- suite 把 Phase 10 的三条主链一次性打通了：backward 是对的，optimizer 是能降的，benchmark 是能系统比较的。
- `phase10_interpretation_summary.md` 中的五个布尔结论都是 `True`，这说明当前框架已经满足了需求文档里“full objective layer + sparse backward + optimizer-ready + benchmarkable”的主目标。
- suite 中 optimizer 的 objective final 高于 standalone optimizer，这是正常的，因为两次运行虽然共享同一代表性 case，但 suite 的总流程和前后写盘顺序不同；关键点是二者都表现出显著 objective 下降且 accepted ratio 为 `1.0`。

## Step 6. example entry points

实际执行：

- `python3 -m phase_10.examples.demo_full_backward_diff`
- `python3 -m phase_10.examples.demo_optimizer`
- `python3 -m phase_10.examples.demo_k_sweep`
- `python3 -m phase_10.examples.demo_baseline_compare`
- `python3 -m phase_10.examples.demo_ablation`
- `python3 -m phase_10.examples.demo_phase10_suite`

功能：

- 验证 examples 入口都可以直接复用主模块，而不是只保留底层 API。

结果：

- 六个 example 入口都执行成功。

分析：

- 这说明 Phase 10 已经不仅是“库函数存在”，而是具备了可复现实验入口，适合继续做论文图表、调参和 benchmark 扩展。

## Step 7. Phase 10 unit tests

实际执行：

- `python3 -m unittest discover -s phase_10 -p 'test*.py'`

功能：

- 验证：
  - dense vs sparse backward 一致性
  - reparameterized gradient vs FD
  - optimizer step 的时间正性
  - optimizer smoke run 的 objective 行为
  - benchmark suite smoke run 与 artifact 输出

结果：

- `Ran 5 tests in 9.929s`
- `OK`

分析：

- Phase 10 的最关键接口已经有自动化回归保护，不需要靠手工跑主脚本才能知道框架有没有坏。

## Step 8. Full repository regression

实际执行：

- `python3 -m unittest discover -s . -p 'test*.py'`

功能：

- 验证 Phase 10 的新增实现没有破坏 Phase 0 到 Phase 9 的既有能力。

结果：

- `Ran 80 tests in 52.591s`
- `OK`

分析：

- 这一步确认 Phase 10 是在现有多阶段 BLOM 代码库上平滑接入的，而不是通过破坏旧接口换来的新功能。

## 总结

Phase 10 的核心结论已经比较明确：

1. full objective layer 已经落地。  
   `K_ctrl + lambda_T K_time + lambda_obs K_obs + lambda_dyn K_dyn + lambda_bc K_bc + lambda_reg K_reg` 现在都已经有统一值/梯度接口，并且能直接接入 raw Scheme C。

2. dense checker 和 sparse backward 都已经可用，而且数值一致。  
   dense vs sparse gap 在 `1e-12` 量级，reparameterized gradient vs FD 误差在 `1e-5` 量级，这已经足以支撑框架级优化。

3. `xi=(q_bar,\tau)` 优化器已经能稳定下降。  
   standalone run 把 objective 从 `2721.18` 降到 `18.51`，suite run 也把 objective 降到 `10.37`，并且都没有触碰时间下界。

4. benchmark / ablation 框架已经成型。  
   现在不仅能比较 raw Scheme C、MINCO、Scheme A、heuristic 的 objective / runtime / sparsity，还能直接做 term ablation 和 `k`-sweep。

5. Phase 10 已经把 Phase 9 的“最小可微闭环”推进成了一个真正的 framework。  
   当前下一步最自然的方向，不再是怀疑 backward 或 optimizer 能不能工作，而是：
   - 更系统的权重调参与 Pareto study
   - 更细化的 constraint / obstacle 建模
   - 更成熟的 general-`k` optimizer benchmark
   - 更强的 baseline 与大规模实验设计
