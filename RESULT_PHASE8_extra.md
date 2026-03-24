# RESULT_PHASE8_extra

本文档记录 Phase 8 supplementary 实验代码的实际运行结果。最终结果目录为 [phase_8/results/phase8_supplementary](/home/dev/code/blom/phase_8/results/phase8_supplementary)。

说明：

- 这次实际执行的主命令是：
  - `python3 -m compileall phase_8`
  - `python3 -m phase_8.blom_phase8_supplementary_suite`
  - `python3 -m unittest discover -s phase_8 -p 'test*.py'`
- `phase_8.blom_phase8_supplementary_suite` 会顺序运行 6 组 supplementary experiments，因此下面的 Step 2.1 到 Step 2.6 都来自这次统一主运行器的真实输出文件，而不是手工推断。
- 最终 suite 配置是：
  - `M_values = [20, 40, 80, 120]`
  - `k_values = [2, 4, 6, 8, 10, 12]`
  - `radius_modes = ['half_k', 'k', 'three_half_k']`
  - `h_values = [1.0]`
  - `nonuniform_boxes = [(0.5, 2.0)]`
  - `n_trials = 50`

## Step 1. `python3 -m compileall phase_8`

功能：

- 编译 `phase_8` 下的 supplementary 模块、suite、examples 和 tests。
- 先确认新增文件的语法和导入闭合。

结果：

- 命令执行成功。
- 新增的 supplementary 文件全部通过编译，包括：
  - [phase8_supplementary_common.py](/home/dev/code/blom/phase_8/phase8_supplementary_common.py)
  - [exp1_large_M_sweep.py](/home/dev/code/blom/phase_8/exp1_large_M_sweep.py)
  - [exp2_more_trials.py](/home/dev/code/blom/phase_8/exp2_more_trials.py)
  - [exp3_radius_sensitivity.py](/home/dev/code/blom/phase_8/exp3_radius_sensitivity.py)
  - [exp4_raw_vs_reference_sanity.py](/home/dev/code/blom/phase_8/exp4_raw_vs_reference_sanity.py)
  - [exp5_empty_interior_risk.py](/home/dev/code/blom/phase_8/exp5_empty_interior_risk.py)
  - [exp6_two_bridge_gap_compare.py](/home/dev/code/blom/phase_8/exp6_two_bridge_gap_compare.py)
  - [blom_phase8_supplementary_suite.py](/home/dev/code/blom/phase_8/blom_phase8_supplementary_suite.py)
  - [test_phase8_supplementary_suite.py](/home/dev/code/blom/phase_8/test_phase8_supplementary_suite.py)

分析：

- 这一轮没有语法或导入错误，说明 supplementary 层已经和现有 Phase 8 主实现顺利接上。

## Step 2. `python3 -m phase_8.blom_phase8_supplementary_suite`

功能：

- 统一执行 6 组 supplementary experiments：
  1. larger-`M` sweep
  2. more-trials robustness
  3. radius sensitivity
  4. raw-vs-reference sanity
  5. empty-interior risk
  6. two-bridge gap comparison
- 自动生成总表 [phase8_supplementary_overview.csv](/home/dev/code/blom/phase_8/results/phase8_supplementary/summary/phase8_supplementary_overview.csv) 和总总结 [phase8_supplementary_summary.md](/home/dev/code/blom/phase_8/results/phase8_supplementary/summary/phase8_supplementary_summary.md)。

结果：

- 脚本执行成功。
- 产物目录完整生成：
  - [exp1_large_M](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp1_large_M)
  - [exp2_more_trials](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp2_more_trials)
  - [exp3_radius_sensitivity](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp3_radius_sensitivity)
  - [exp4_raw_vs_reference_sanity](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp4_raw_vs_reference_sanity)
  - [exp5_empty_interior_risk](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp5_empty_interior_risk)
  - [exp6_two_bridge_gaps](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp6_two_bridge_gaps)
  - [summary](/home/dev/code/blom/phase_8/results/phase8_supplementary/summary)

分析：

- supplementary suite 已经把需求文档要求的 6 类实验全部跑完，并且自动落盘了图、CSV、JSON 和 summary。

## Step 2.1. Exp 1: larger-`M` sweep with nonempty interior

功能：

- 在更大的 `M = [20, 40, 80, 120]` 下重复 Phase 8 的 matching 实验。
- 显式记录 `interior_count`、`boundary_count`、`full/interior/boundary matching` 和 `boundary_energy_ratio`。
- 核心是排查“小 `M` + 空 interior”带来的假象。

结果：

- 对应文件：
  - [summary_large_M.csv](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp1_large_M/summary_large_M.csv)
  - [summary_large_M.json](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp1_large_M/summary_large_M.json)
  - [phase8_largeM_matching_curves.png](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp1_large_M/phase8_largeM_matching_curves.png)
  - [phase8_interior_count_vs_k.png](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp1_large_M/phase8_interior_count_vs_k.png)
- 关键数值：
  - 总记录数 `24`
  - nonempty-interior cases `22 / 24`
  - mean boundary energy ratio `0.424942`
  - nonempty cases 的 mean `interior/full` ratio `0.773649`
  - 分 `M` 看：
    - `M=20`: nonempty `4/6`，mean boundary ratio `0.787945`
    - `M=40`: nonempty `6/6`，mean boundary ratio `0.530674`
    - `M=80`: nonempty `6/6`，mean boundary ratio `0.263034`
    - `M=120`: nonempty `6/6`，mean boundary ratio `0.118116`

分析：

- 最重要的结论是：即使在 interior 不为空的更大规模下，interior-only matching 仍然明显优于 full matching，所以 Phase 8 的 interior-first 现象不是单纯由空 interior 造成的。
- 另一个有意思的现象是 boundary ratio 会随着 `M` 增大而下降。这说明 boundary layer 仍然重要，但在更长轨迹上，它对总误差的相对占比不再像小规模案例里那么极端。
- 换句话说，Phase 8 原先的 boundary-dominance 结论在小规模代表例里更强，而 supplementary 的 larger-`M` 实验告诉我们：这个结论是成立的，但强度会随规模变化。

## Step 2.2. Exp 2: more random-trials robustness

功能：

- 在 `uniform_h=1.0` 和 `bounded_[0.5, 2.0]` 两种 regime 下，各做 `50` 个随机 trial。
- 对每个 trial 的 `full/interior/boundary` matching 曲线做 log-fit，再统计 slope 分布与 boundary ratio 分布。

结果：

- 对应文件：
  - [phase8_more_trials_aggregated.csv](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp2_more_trials/phase8_more_trials_aggregated.csv)
  - [phase8_more_trials_fit_rows.csv](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp2_more_trials/phase8_more_trials_fit_rows.csv)
  - [phase8_interior_slope_boxplot.png](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp2_more_trials/phase8_interior_slope_boxplot.png)
  - [phase8_boundary_ratio_boxplot.png](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp2_more_trials/phase8_boundary_ratio_boxplot.png)
- 关键数值：
  - fit rows `100`
  - mean interior slope `-4.147938`
  - mean full slope `-0.369023`
  - share of trials with `interior_slope < full_slope`: `1.0`
  - mean boundary energy ratio `0.736434`
  - 分 regime：
    - `uniform_h=1.00`: mean interior slope `-4.140161`, mean full slope `-0.373042`
    - `bounded_[0.50,2.00]`: mean interior slope `-4.155714`, mean full slope `-0.365005`

分析：

- 这是本轮 supplementary 里最强的稳健性证据之一。`100` 个 trial 里 interior slope 全部都比 full slope 更负，说明 interior-first 不是偶然样本，而是稳定统计规律。
- uniform 和 bounded-nonuniform 的结果非常接近，这也比原 Phase 8 更强地支持了“该现象不依赖极端理想时间配置”的判断。
- 同时，mean boundary ratio 仍然接近 `0.74`，说明 boundary layer 仍是主要误差来源之一。

## Step 2.3. Exp 3: boundary radius sensitivity

功能：

- 比较三种 interior/boundary 划分：
  - `r(k)=ceil(k/2)`
  - `r(k)=k`
  - `r(k)=ceil(3k/2)`
- 检查 interior-first 是否只是某个 radius 定义的人造结果。

结果：

- 对应文件：
  - [phase8_radius_sensitivity.csv](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp3_radius_sensitivity/phase8_radius_sensitivity.csv)
  - [phase8_radius_sensitivity_fit_rows.csv](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp3_radius_sensitivity/phase8_radius_sensitivity_fit_rows.csv)
  - [phase8_radius_mode_interior_matching.png](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp3_radius_sensitivity/phase8_radius_mode_interior_matching.png)
  - [phase8_radius_sensitivity_heatmap.png](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp3_radius_sensitivity/phase8_radius_sensitivity_heatmap.png)
- 关键数值：
  - `half_k`: interior slope `-0.326877`, full slope `-0.335389`
  - `k`: interior slope `-0.344816`, full slope `-0.335389`
  - `three_half_k`: interior slope `-0.412172`, full slope `-0.335389`

分析：

- `r(k)=k` 和 `r(k)=ceil(3k/2)` 下，interior slope 都比 full slope 更负，说明 interior-first 不是只在默认半径上成立。
- `r(k)=ceil(k/2)` 下 interior slope 反而略弱于 full slope，说明 radius 的确会影响现象强度，不能把 “interior 的定义” 当成无关选择。
- 这组实验给出的最稳妥结论是：Phase 8 的现象是存在的，但未来 theorem 里对 interior set 的定义需要更保守、更清楚。

## Step 2.4. Exp 4: `raw vs reference-window` sanity check

功能：

- 在更大的 `M`、更多 trial、不同 regime 下重复验证：
  - `||c^{C,k} - \hat c^{(k)}||`
- 同时人工扰动 inherited `gamma*`，验证 reference-window 解对边界 trace 不是“锁死不动”。

结果：

- 对应文件：
  - [phase8_raw_vs_reference_stats.csv](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp4_raw_vs_reference_sanity/phase8_raw_vs_reference_stats.csv)
  - [phase8_reference_gamma_sensitivity.csv](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp4_raw_vs_reference_sanity/phase8_reference_gamma_sensitivity.csv)
  - [phase8_raw_vs_reference_error.png](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp4_raw_vs_reference_sanity/phase8_raw_vs_reference_error.png)
  - [phase8_reference_gamma_loglog.png](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp4_raw_vs_reference_sanity/phase8_reference_gamma_loglog.png)
- 关键数值：
  - mean raw-to-reference error `2.187265e-13`
  - max raw-to-reference error `6.702577e-12`
  - mean bridge gap ratio `5.851218e+12`
  - mean gamma amplification ratio `9.371873e-10`
  - mean gamma response `3.400386e-15`

分析：

- 这组 sanity check 基本把 “raw vs reference-window ≈ 0” 坐实了。即使放大到更大 `M` 和更多 trial，这个 gap 仍然维持在机器精度附近。
- 同时，gamma 扰动不是完全零响应，所以 reference-window 不是因为实现退化才和 raw Scheme C 一样，而是两者在当前 canonical 设置下确实高度一致。
- 这也是为什么 `bridge gap ratio` 会巨大到 `1e12` 量级：真正大的 gap 根本不在 `raw -> ref` 这条桥上，而是在 `ref -> ideal` 上。

## Step 2.5. Exp 5: empty-interior false-positive risk

功能：

- 显式记录：
  - `interior_count`
  - `interior_fraction`
  - `is_empty_interior`
- 检查哪些 `(M, k, r(k))` 组合会让 interior 为空，从而导致 “interior 误差很好看” 的假阳性。

结果：

- 对应文件：
  - [phase8_empty_interior_risk.csv](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp5_empty_interior_risk/phase8_empty_interior_risk.csv)
  - [phase8_interior_fraction_vs_k.png](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp5_empty_interior_risk/phase8_interior_fraction_vs_k.png)
  - [phase8_interior_error_with_empty_flags.png](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp5_empty_interior_risk/phase8_interior_error_with_empty_flags.png)
  - [phase8_empty_interior_risk.png](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp5_empty_interior_risk/phase8_empty_interior_risk.png)
- 关键数值：
  - empty-interior cases `13 / 87`
  - min interior fraction `0.0`
  - max interior fraction `0.983333`
  - 按 radius mode 分：
    - `half_k`: `1 / 29`
    - `k`: `5 / 29`
    - `three_half_k`: `7 / 29`

分析：

- 这一步很重要，因为它明确告诉我们：empty-interior 风险是真实存在的，而且越保守的半径定义越容易把 interior 吃空。
- 但它同时也表明 supplementary 的主要正结论并不是靠 empty-interior 撑出来的，因为我们在 Exp 1 和 Exp 2 里已经看到了大量 nonempty 情况下的稳定 interior-first。
- 后续写论文时，凡是引用 interior-only 误差曲线的图，都应该加上 interior cardinality 或 empty-set 注释。

## Step 2.6. Exp 6: two-bridge gap comparison

功能：

- 系统比较两条桥上的 gap：
  - `E_raw_to_ref = ||c^{C,k} - \hat c^{(k)}||`
  - `E_ref_to_ideal = ||\hat c^{(k)} - \tilde c^{(k)}||`
- 直接回答下一条 theorem 应该优先补哪一条桥。

结果：

- 对应文件：
  - [phase8_two_bridge_gap_stats.csv](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp6_two_bridge_gaps/phase8_two_bridge_gap_stats.csv)
  - [phase8_two_bridge_gaps.png](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp6_two_bridge_gaps/phase8_two_bridge_gaps.png)
  - [phase8_gap_ratio_vs_k.png](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp6_two_bridge_gaps/phase8_gap_ratio_vs_k.png)
  - [phase8_gap_ratio_boxplot.png](/home/dev/code/blom/phase_8/results/phase8_supplementary/exp6_two_bridge_gaps/phase8_gap_ratio_boxplot.png)
- 关键数值：
  - mean gap ratio `3.541086e+12`
  - median gap ratio `1.408477e+12`
  - min gap ratio `9.264042e+10`
  - max gap ratio `5.398599e+13`
  - mean raw-to-reference gap `1.308464e-13`
  - mean reference-to-ideal gap `3.548109`
  - 分 regime：
    - `uniform_h=1.00`: mean raw->ref `6.968578e-14`, mean ref->ideal `2.942965`
    - `bounded_[0.50,2.00]`: mean raw->ref `1.920070e-13`, mean ref->ideal `4.153254`

分析：

- 这一组结果几乎是“直接判决”。`raw -> ref` 始终在机器精度附近，而 `ref -> ideal` 仍然是 `O(1)` 量级。
- 因此当前最应该优先补的理论桥，明确就是：
  - `reference-window / ideal-truncation consistency`
- 这也意味着如果接下来还要补 theorem，优先级已经不再模糊。

## Step 2.7. Supplementary suite 总结

功能：

- 把 6 组实验压缩成一个可以直接用于决策的总 summary。

结果：

- 总 summary 文件：
  - [phase8_supplementary_overview.csv](/home/dev/code/blom/phase_8/results/phase8_supplementary/summary/phase8_supplementary_overview.csv)
  - [phase8_supplementary_summary.md](/home/dev/code/blom/phase_8/results/phase8_supplementary/summary/phase8_supplementary_summary.md)
- suite 自动回答：
  - `Larger M and more trials remain supportive: True`
  - `Interior-first still holds without collapsing into empty-interior false positives: True`
  - ``raw vs reference-window ≈ 0` remains credible: True`
  - `Recommended next bridge theorem: reference-window / ideal-truncation consistency`
  - `Strong paper-level observation already supported: True`

分析：

- 这一步说明 supplementary 实验并没有推翻原 Phase 8，反而把原来的方向性判断进一步加固了。
- 最重要的新增信息不是 “现象还在”，而是：
  - interior-first 现象经受住了 larger-`M`、more-trials 和 empty-interior 风险排查；
  - `raw vs reference-window ≈ 0` 经受住了 sanity check；
  - 下一条 bridge theorem 的优先级已经非常明确。

## Step 3. `python3 -m unittest discover -s phase_8 -p 'test*.py'`

功能：

- 回归测试 Phase 8 主实验和 supplementary 扩展层。
- 确认新增实验没有破坏原有 Phase 8 验证链路。

结果：

- `Ran 12 tests in 7.078s`
- `OK`

分析：

- 说明当前 Phase 8 目录下的主脚本、supplementary 脚本和新增 smoke tests 都通过了。

## 总结

这次 Phase 8 supplementary 的核心结论可以压缩成三点：

1. `interior-first` 现象在 larger-`M`、more-trials 和多种 radius 定义下仍然成立，不是空 interior 的假阳性。
2. `raw Scheme C -> reference-window` 的 gap 在数值上几乎可以忽略；真正大的桥是 `reference-window -> ideal truncation`。
3. 就当前数值证据看，Phase 8 已经足以支撑一个比较强的 paper-level observation，而下一条最值得优先补的 theorem 是 `reference-window / ideal-truncation consistency`，而不是回头重新解释 raw Scheme C 本身。
