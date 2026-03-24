# RESULT_PHASE9

本文档记录 Phase 9 最小可微优化闭环的实际运行结果。最终结果目录为 [phase_9/results/phase9_validation](/home/dev/code/blom/phase_9/results/phase9_validation)。

说明：

- Phase 9 的主模块、suite 和 examples 都实际执行过。
- 为保证结果一致性，本文档引用的最终数值以最后一次 `python3 -m phase_9.blom_phase9_validation_suite` 生成的结果为准。
- 当前 canonical loop 固定为 `s = 4, k = 2`，只使用 raw Scheme C。

## Step 1. `python3 -m compileall phase_9`

功能：

- 编译 `phase_9` 下的主模块、examples 和 tests。
- 提前排除语法错误和模块引用错误。

结果：

- 命令执行成功。
- `phase_9` 下所有新增文件都通过编译。

分析：

- 这一步确认 Phase 9 的文件组织和导入链路已经闭合，后面的验证重点就可以放在 Jacobian、backward 和优化行为本身。

## Step 2. `python3 -m phase_9.blom_backward_diff`

功能：

- 运行 Phase 9 的核心 backward 模块。
- 实际执行：
  - raw Scheme C 系数计算
  - raw Scheme C Jacobian 计算
  - total objective evaluation
  - dense / banded backward
  - gradcheck 图表与 summary 输出

结果：

- 命令执行成功。
- 产物保存到：
  - [phase9_gradcheck_summary.csv](/home/dev/code/blom/phase_9/results/phase9_validation/gradcheck/phase9_gradcheck_summary.csv)
  - [phase9_gradcheck_summary.md](/home/dev/code/blom/phase_9/results/phase9_validation/gradcheck/phase9_gradcheck_summary.md)
  - [phase9_dense_vs_banded_grad.png](/home/dev/code/blom/phase_9/results/phase9_validation/gradcheck/phase9_dense_vs_banded_grad.png)
  - [phase9_jacobian_sparsity.png](/home/dev/code/blom/phase_9/results/phase9_validation/gradcheck/phase9_jacobian_sparsity.png)
- 关键数值：
  - total dense vs FD max abs error `3.755835e-05`
  - total banded vs FD max abs error `3.755835e-05`
  - dense vs banded max abs gap `1.818989e-12`
  - control-term partial gradient max abs error `1.485118e-05`
  - time-term partial gradient max abs error `2.631780e-11`
  - obs-term partial gradient max abs error `4.200990e-09`

分析：

- 总梯度对 finite difference 的误差在 `1e-5` 到 `1e-4` 量级，已经足够支持这个最小优化闭环。
- dense 与 banded backward 的差距在 `1e-12` 量级，说明 block-banded accumulation 只是去掉了精确为零的项，没有引入额外误差。
- 在三个子项里，control term 的 partial gradient 误差最大，这也合理，因为它同时依赖高阶 snap cost matrix 和 raw Scheme C Jacobian；time 和 obs 的误差都更小。

## Step 3. `python3 -m phase_9.blom_space_time_opt_demo`

功能：

- 运行最小优化 demo。
- 实际执行：
  - raw Scheme C objective
  - banded backward gradient
  - projected / clipped time update
  - objective history / grad norm history / trajectory before-after 可视化

结果：

- 命令执行成功。
- 产物保存到：
  - [phase9_opt_history.csv](/home/dev/code/blom/phase_9/results/phase9_validation/optimization_demo/phase9_opt_history.csv)
  - [phase9_opt_summary.md](/home/dev/code/blom/phase_9/results/phase9_validation/optimization_demo/phase9_opt_summary.md)
  - [phase9_objective_curve.png](/home/dev/code/blom/phase_9/results/phase9_validation/optimization_demo/phase9_objective_curve.png)
  - [phase9_gradnorm_curve.png](/home/dev/code/blom/phase_9/results/phase9_validation/optimization_demo/phase9_gradnorm_curve.png)
  - [phase9_traj_before_after.png](/home/dev/code/blom/phase_9/results/phase9_validation/optimization_demo/phase9_traj_before_after.png)
- 最终 suite 一致化后的关键数值：
  - objective drop `2.688224e+03`
  - control drop `2.688807e+03`
  - time drop `-1.802551e-01`
  - obs drop `-4.023269e-01`
  - grad-norm drop proxy `1.243217e+04`
  - any duration touched lower bound: `False`
- `phase9_opt_history.csv` 里的前几步：
  - step 0 objective `2689.573697`
  - step 1 objective `2567.491363`
  - step 2 objective `1371.509551`
  - step 3 objective `231.308823`
  - step 4 objective `200.531322`
  - step 25 objective `1.349764`

分析：

- objective 单调下降，说明这条最小闭环已经能形成稳定的 descent loop。
- 最大的下降来自 control term，这意味着当前 demo 的主驱动力是把高 snap 代价迅速压下去。
- time term 和 obs term 在这组默认权重下略有上升，但都没有破坏总目标下降；这说明当前 demo 更像 “先验证 backward/optimization 能通”，而不是“已经完成多项权重平衡”。
- 最重要的是，`T_i` 没有触碰下界，说明简单的时间裁剪已经足够保证这一版最小 demo 的可行性。

## Step 4. `python3 -m phase_9.blom_phase9_validation_suite`

功能：

- 统一执行：
  - gradient check
  - dense vs banded 对比
  - minimal optimization demo
  - 总 overview 和 interpretation summary

结果：

- 命令执行成功。
- 总汇总文件保存到：
  - [phase9_overview.csv](/home/dev/code/blom/phase_9/results/phase9_validation/compare/phase9_overview.csv)
  - [phase9_interpretation_summary.md](/home/dev/code/blom/phase_9/results/phase9_validation/compare/phase9_interpretation_summary.md)
  - [phase9_suite_summary.json](/home/dev/code/blom/phase_9/results/phase9_validation/compare/phase9_suite_summary.json)
- overview 关键数值：
  - `gradcheck_dense_vs_fd`: max abs error `3.755834768526256e-05`
  - `gradcheck_banded_vs_fd`: max abs error `3.755834586627316e-05`
  - `dense_banded_gap`: `1.8189894035458565e-12`
  - `optimization_demo`: objective drop `2688.2239335085405`

分析：

- suite 把三件最关键的事情同时确认了：
  1. analytic backward 链条是对的；
  2. banded backward 和 dense checker 数值一致；
  3. 最小优化 loop 真能下降。
- 这正是 Phase 9 需求文档要回答的四个核心问题里的前三个。

## Step 5. example entry points

实际执行：

- `python3 -m phase_9.examples.demo_backward_diff`
- `python3 -m phase_9.examples.demo_gradcheck`
- `python3 -m phase_9.examples.demo_minimal_opt`
- `python3 -m phase_9.examples.demo_phase9_suite`

功能：

- 验证 examples 入口可直接复用主模块。

结果：

- 四个 example 入口都执行成功。

分析：

- 这说明 Phase 9 的交付不只是底层函数存在，而是已经具有可直接演示和复现的入口层。

## Step 6. Phase 9 unit tests

实际执行：

- `python3 -m unittest phase_9.test_blom_backward_diff`

功能：

- 验证：
  - raw Scheme C 系数/Jacobian shape
  - dense equals banded
  - control/time/obs 三项梯度
  - total gradient vs FD
  - 最小优化 demo smoke
  - validation suite smoke

结果：

- `Ran 8 tests in 1.829s`
- `OK`

分析：

- Phase 9 的最关键接口和数值链条都已经有自动化回归保护。

## Step 7. Full repository regression

实际执行：

- `python3 -m unittest discover -s . -p 'test*.py'`

功能：

- 验证本次 Phase 9 增量实现没有破坏前面各阶段。

结果：

- `Ran 71 tests in 36.323s`
- `OK`

分析：

- 这一步确认 Phase 9 是在现有多阶段 BLOM 代码库上平滑接入的，而不是通过破坏旧接口换来的新功能。

## 总结

Phase 9 的核心结论已经比较明确：

1. raw Scheme C 的最小可微闭环是成立的。  
   `analytic vs FD` 误差在 `1e-5` 到 `1e-4` 量级，说明 `theta -> c(theta) -> K(c,T) -> grad_theta J` 这条链条已经数值可靠。

2. block-banded backward 与 dense checker 完全一致。  
   二者 gap 在 `1e-12` 量级，说明 Phase 6 的局部支撑结果已经可以直接转化成 Phase 9 的真实 backward implementation。

3. 最小优化 demo 可以稳定下降。  
   当前默认配置下，objective 从 `2689.57` 降到 `1.35`，而且没有触碰时间下界。

4. 这套最小闭环已经足以支撑下一阶段。  
   当前最合理的下一步不是再怀疑 raw Scheme C 能不能做优化，而是进入：
   - 更完整的 block-banded optimization framework
   - 更复杂的约束/障碍模型
   - 更系统的时间分配与正则化设计

换句话说，Phase 9 已经把 “BLOM 只是一个局部轨迹表示” 推进到了 “BLOM 已经能作为一个可反向传播、可做 descent 的优化对象”。  
