# RESULT_PHASE8

本文档记录 Phase 8 代码的实际运行结果。最终结果目录为 [phase_8/results/phase8_validation](/home/dev/code/blom/phase_8/results/phase8_validation)。

说明：

- `phase_8` 的单案例脚本和 examples 会刷新同一批结果文件。
- 本文档引用的最终 CSV / JSON 数值以最后一次 `python3 -m phase_8.blom_phase8_validation_suite` 和 `python3 -m phase_8.examples.demo_phase8_suite` 生成的结果为准。
- 最终 suite 配置是 `M = 10`、`k = [2, 4, 6, 8, 10]`、`n_trials = 6`。

## Step 1. `python3 -m compileall phase_8`

功能：

- 编译 `phase_8` 下的主模块、examples 和 tests。
- 用最便宜的方式先检查语法层面是否完整。

结果：

- 命令执行成功。
- `phase_8` 下所有新增文件都通过编译，没有语法错误。

分析：

- 这一步确认 Phase 8 的模块组织已经闭合，后续错误若出现，基本只会来自数值逻辑或接口对接，而不是文件结构问题。

## Step 2. `python3 -m phase_8.blom_interior_matching_check`

功能：

- 对固定 `(q, T)` 和多组 `k` 计算：
  - full matching error
  - interior-only matching error
  - boundary-only matching error
  - segmentwise matching heatmap
  - `log(error)` vs `k` 拟合

结果：

- 脚本执行成功。
- 产物保存到：
  - [interior_matching_by_k.csv](/home/dev/code/blom/phase_8/results/phase8_validation/interior_matching/interior_matching_by_k.csv)
  - [interior_matching_fit_summary.csv](/home/dev/code/blom/phase_8/results/phase8_validation/interior_matching/interior_matching_fit_summary.csv)
  - [phase8_interior_matching_summary.md](/home/dev/code/blom/phase_8/results/phase8_validation/interior_matching/phase8_interior_matching_summary.md)
- 最终 suite 刷新的关键数值：
  - `k=2`: full `4.8933`, interior `3.0907`, boundary `3.7937`
  - `k=4`: full `1.7308`, interior `0.5716`, boundary `1.6337`
  - `k=6`: full `0.5934`, interior `0.0`, boundary `0.5934`
  - `k=8`: full `0.2873`, interior `0.0`, boundary `0.2873`
  - `k=10`: full `0.1381`, interior `0.0`, boundary `0.1381`
- 拟合结果：
  - full slope `-0.446569`, `R^2 = 0.990902`
  - interior slope `-5.265687`, `R^2 = 0.761420`
  - boundary slope `-0.418230`, `R^2 = 0.995858`

分析：

- `interior-only matching` 在全部 `5/5` 个测试点上都明显优于 full matching。
- interior slope 比 full slope 更负很多，说明误差衰减在 interior 部分明显更快，这正是 “interior-first matching” 想要抓住的现象。
- 从 `k >= 6` 开始，当前 `r(k)=k` 定义下 interior 集合已经为空，所以最终误差完全落在 boundary layer。这个结果不是坏事，反而说明 boundary-layer projector 的定义在实验中确实起到了“隔离边界区”的作用。

## Step 3. `python3 -m phase_8.blom_boundary_gap_decomposition`

功能：

- 把 full matching gap 分解为 interior part 和 boundary-layer part。
- 额外计算 Phase 8 理论里更细的两阶段对象：
  - raw vs reference-window
  - reference-window vs ideal-truncated

结果：

- 脚本执行成功。
- 产物保存到：
  - [boundary_gap_decomposition.csv](/home/dev/code/blom/phase_8/results/phase8_validation/boundary_gap/boundary_gap_decomposition.csv)
  - [boundary_radius_sensitivity.csv](/home/dev/code/blom/phase_8/results/phase8_validation/boundary_gap/boundary_radius_sensitivity.csv)
  - [phase8_boundary_gap_summary.md](/home/dev/code/blom/phase_8/results/phase8_validation/boundary_gap/phase8_boundary_gap_summary.md)
- 关键数值：
  - mean boundary energy ratio `0.898402`
  - `k=2` boundary energy ratio `0.601060`
  - `k=4` boundary energy ratio `0.890952`
  - `k=6,8,10` boundary energy ratio `1.0`
  - mean raw-vs-reference gap `4.042120e-14`
  - mean reference-window-vs-ideal gap `1.528579e+00`

分析：

- 从 projector 分解看，full matching gap 的大部分能量确实在 boundary layer 中，这强烈支持后续把 boundary remainder 单独理论化。
- 一个很有意思的现象是：`raw vs reference-window` 几乎是机器精度，而 `reference-window vs ideal` 反而承担了主要差异。这说明在当前实现和样例上，真正大的 gap 更像是 “reference-window local model” 与 “ideal truncation operator” 的差，而不是“natural artificial boundary vs exact inherited gamma” 的差。
- 这给下一步理论一个很明确的提示：除了 boundary-layer theorem 之外，`reference-window / ideal-truncation consistency proposition` 也很值得优先推进。

## Step 4. `python3 -m phase_8.blom_uniform_vs_nonuniform_interior`

功能：

- 在 uniform-time 与 bounded-nonuniform 两类时间配置下，重复 interior-first matching 实验。
- 统计每个 regime 的：
  - full / interior / boundary matching error
  - interior slope
  - interior fit `R^2`
  - boundary energy ratio

结果：

- 脚本执行成功。
- 产物保存到：
  - [uniform_vs_nonuniform_summary.csv](/home/dev/code/blom/phase_8/results/phase8_validation/uniform_vs_nonuniform/uniform_vs_nonuniform_summary.csv)
  - [uniform_vs_nonuniform_fit_summary.csv](/home/dev/code/blom/phase_8/results/phase8_validation/uniform_vs_nonuniform/uniform_vs_nonuniform_fit_summary.csv)
  - [phase8_uniform_vs_nonuniform_summary.md](/home/dev/code/blom/phase_8/results/phase8_validation/uniform_vs_nonuniform/phase8_uniform_vs_nonuniform_summary.md)
- 最终 suite 刷新的代表性数值：
  - uniform mean interior slope `-5.385500`
  - bounded-nonuniform mean interior slope `-5.309442`
  - `uniform_h=0.50` 的 interior slopes 大约在 `-5.68` 到 `-5.82`
  - `uniform_h=1.00` 的 interior slopes 大约在 `-5.26` 到 `-5.36`
  - bounded cases 也都稳定为负，主要落在 `-5.06` 到 `-5.61`

分析：

- interior-first 现象不只存在于 uniform-time，也在 bounded-nonuniform 下稳定存在。
- uniform-time 的平均 slope 更负，说明如果要先写一个最干净的 theorem，确实应该先从 uniform regime 下手。
- bounded-nonuniform 的 slope 仍然稳健为负，没有出现“趋势消失”的情况，所以后续推广并不是拍脑袋，而是已经有比较扎实的数值支撑。

## Step 5. `python3 -m phase_8.blom_phase8_validation_suite`

功能：

- 统一执行：
  - single-case interior matching
  - single-case boundary-gap decomposition
  - uniform vs bounded-nonuniform regime split
  - 总 overview 表与自动解释 summary

结果：

- 脚本执行成功。
- 总汇总文件保存到：
  - [phase8_overview.csv](/home/dev/code/blom/phase_8/results/phase8_validation/compare/phase8_overview.csv)
  - [phase8_interpretation_summary.md](/home/dev/code/blom/phase_8/results/phase8_validation/compare/phase8_interpretation_summary.md)
  - [summary_phase8_suite.json](/home/dev/code/blom/phase_8/results/phase8_validation/compare/summary_phase8_suite.json)
- overview 总结出的关键结论：
  - `Interior-first matching significant: True`
  - mean boundary energy ratio `0.898402`
  - uniform mean interior slope `-5.385500`
  - bounded mean interior slope `-5.309442`

分析：

- 这一步把三个分散实验收束成了一致的判断：Phase 8 的主现象已经足够清楚，不再是“full-domain 还差一点”的模糊状态。
- 目前最值得继续补的桥梁不是重新定义 raw Scheme C，而是：
  - boundary-layer remainder theorem
  - reference-window / ideal-truncation consistency proposition
- 也就是说，Phase 8 的数值工作已经把“下一步理论应该补哪一条桥”这件事具体化了。

## Step 6. examples smoke run

实际执行：

- `python3 -m phase_8.examples.demo_interior_matching`
- `python3 -m phase_8.examples.demo_boundary_gap`
- `python3 -m phase_8.examples.demo_uniform_vs_nonuniform`
- `python3 -m phase_8.examples.demo_phase8_suite`

功能：

- 验证 examples 入口可以直接复用主模块运行。

结果：

- 四个 example 入口都执行成功，没有额外报错。

分析：

- 这说明 Phase 8 交付不只是“主脚本能跑”，而是已经具备了对外演示和复现实验的入口层。

## Step 7. Phase 8 unit tests

实际执行：

- `python3 -m unittest phase_8.test_blom_interior_matching_check phase_8.test_blom_boundary_gap_decomposition phase_8.test_blom_uniform_vs_nonuniform_interior`

功能：

- 验证 Phase 8 的核心函数：
  - interior / boundary index set
  - projector 分解
  - augmented local boundary response
  - 三个主脚本的 smoke run

结果：

- `Ran 8 tests in 1.315s`
- `OK`

分析：

- Phase 8 的新增公共层和三个实验脚本都至少通过了基本行为验证。

## Step 8. Full repository regression

实际执行：

- `python3 -m unittest discover -s . -p 'test*.py'`

功能：

- 验证本次 Phase 8 改动没有破坏前面各阶段。

结果：

- `Ran 63 tests in 35.014s`
- `OK`

分析：

- 这一步确认 Phase 8 是在现有多阶段代码库上“增量接入”成功，而不是靠破坏旧接口换来新实验。

## 总结

Phase 8 已经把 matching gap 从“总体上还存在差异”拆成了两个更可操作的层次：

1. `interior-first matching` 数值上显著成立。
2. `boundary-layer` 在 full matching gap 中占主导地位。

同时，两阶段分解还暴露出一个更具体的下一步问题：

- `raw vs reference-window` 在当前实验里几乎为零；
- `reference-window vs ideal-truncation` 才是更大的剩余项。

因此，当前最值得继续补的理论桥梁是：

1. boundary-layer remainder theorem
2. reference-window / ideal-truncation consistency proposition

而不是先回头重写 raw Scheme C。  
