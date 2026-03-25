你现在是一个非常严格的“科学可视化 + 交互式研究 demo + Python 工程实现”开发者。请基于我前面已经完成的 BLOM 各 phases 的理论与代码框架，生成一个 **Python 版本的、科学准确的、可交互的 demo**。你的目标不是做一个好看的概念玩具，而是做一个**与现有 BLOM 数学对象、数值结果、phase 结论严格一致**的研究型演示系统。

【总目标】

请实现一个 Python 可视化 demo，用来展示并比较：

1. BLOM（优先 raw Scheme C）
2. MINCO（global reference）
3. B-spline（作为局部几何控制 baseline）
4. 如有必要，可选显示 Scheme A / Scheme B 作为 continuity-vs-locality 对照

这个 demo 必须能科学地展示 BLOM 的核心特性，而不是只做概念示意。必须尽量复用前面 phases 已经建立的 BLOM 理论对象、接口和结果，不允许随意发明新的数学定义或用手工权重伪装成求解结果。

【最重要的原则】

1. 所有可视化都必须尽量建立在“真实求解结果”之上，而不是拍脑袋的示意曲线。
2. 如果某个对象只能做近似或占位实现，必须在代码和界面中明确标注“示意模式”。
3. BLOM 的局部性展示，优先使用：
   - 系数块 \(c_i\)
   - Jacobian 稀疏模式
   - support width
   - perturbation response map
   而不是只靠“肉眼看曲线影响范围”。
4. 必须与前面 phases 的符号和对象保持一致：
   - \(\bar q\), \(T\), \(c^\star\), \(\tilde c^{(k)}\), \(c^{\mathrm{C},k}\)
   - raw Scheme C
   - reference-window
   - ideal truncation
   - support sets
   - block-banded Jacobian
5. 不允许把 B-spline 的 control points 和 MINCO/BLOM 的 waypoints 混为一谈；界面上必须明确区分变量层级。
6. demo 的结论表达必须克制，不能把“视觉现象”写成 theorem。

【实现形式要求】

请使用 Python 实现一个本地可运行的交互 demo，优先选择以下技术路线之一：

- 首选：Streamlit + Plotly
- 或：Dash + Plotly
- 如确有必要，也可用 Panel / Bokeh

要求：
- 交互顺畅
- 支持 slider / dropdown / checkbox / tabs
- 图可缩放、悬停查看数值
- 结构清晰，适合答辩演示和补充材料展示

【代码结构要求】

请生成一个完整但清晰的项目结构，建议类似：

demo_bloom/
├── app.py
├── demo_core/
│   ├── data_interface.py
│   ├── blom_adapter.py
│   ├── minco_adapter.py
│   ├── bspline_adapter.py
│   ├── jacobian_tools.py
│   ├── benchmark_tools.py
│   └── plotting.py
├── pages/
│   ├── 01_overview.py
│   ├── 02_single_perturbation.py
│   ├── 03_jacobian_sparsity.py
│   ├── 04_k_sweep.py
│   ├── 05_scheme_compare.py
│   ├── 06_matching_bridge.py
│   ├── 07_optimizer_demo.py
│   └── 08_benchmark_summary.py
├── assets/
├── results_cache/
├── requirements.txt
└── README.md

要求你在输出中给出：
1. 完整代码
2. 运行说明
3. 每个页面展示什么、为什么这样展示
4. 如果某部分必须调用我前面 phase 的已有接口，请明确写出接口假设

【必须复用/对齐的前置 BLOM 对象】

请假定我已经有前面 phases 中的以下函数或等价接口；你要优先适配它们，而不是重写 BLOM 数学：

1. MINCO 全局 reference：
   - compute_minco_reference(q, T, ...)
   - 返回 \(c^\star\)

2. ideal truncation：
   - compute_ideal_truncated_blom_k(q, T, k, ...)
   - 返回 \(\tilde c^{(k)}\)

3. actual raw Scheme C：
   - compute_actual_blom_k(q, T, k, scheme=\"C\", ...)
   - 返回 \(c^{\mathrm{C},k}\)

4. raw Scheme C Jacobians：
   - compute_raw_schemeC_jacobians(...)
   - 返回 \(\partial c/\partial \bar q\), \(\partial c/\partial T\)

5. optimizer / backward diff：
   - backward_diff_dense(...)
   - backward_diff_banded(...)
   - evaluate_full_objective(...)
   - run_space_time_optimization(...)

如果接口名略有不同，请你设计适配层（adapter），不要强行改动底层实现。

【demo 必须包含的核心页面】

请至少包含以下 8 类页面或 tab：

### 页面 1：Overview / 方法总览
展示：
- B-spline、MINCO、BLOM 的定位差别
- 哪些是插值，哪些不是
- 哪些有 exact local support
- 哪些 Jacobian 稠密 / block-banded
- 哪些有 global variational meaning

要求：
- 用表格 + 简图展示
- 明确注明这是“结构对比”，不是 theorem 图

### 页面 2：Single perturbation response（单点扰动响应）
目标：
展示当只扰动一个 waypoint（对 MINCO/BLOM）或一个 control point（对 B-spline）时，曲线如何变化。

要求：
- 支持调节总点数 \(M\)
- 支持选择被扰动索引
- 支持调节扰动大小
- 支持调节 BLOM 的 \(k\)
- 必须同时显示：
  - 原始轨迹
  - 扰动后轨迹
  - 差值热图 / 局部响应图
- 必须明确区分：
  - B-spline 用的是 control points
  - MINCO/BLOM 用的是 interpolation waypoints

注意：
- 这个页面只能作为“几何直观”
- 不要把它当 BLOM 局部性的主要证据页

### 页面 3：Jacobian sparsity / support pattern
这是最重要的页面之一。

必须展示：
- \(\partial c / \partial \bar q\)
- \(\partial c / \partial T\)

对比对象：
- raw Scheme C
- MINCO
- 如方便可加 B-spline

要求：
- 用 heatmap 展示绝对值
- 允许切换 full matrix / thresholded sparsity mask
- 允许查看 block-level sparsity
- 显示：
  - support width
  - bandwidth
  - outside-band energy
- 这页必须作为 BLOM 局部性优势的主要展示页

### 页面 4：k-sweep / locality-approximation trade-off
目标：
展示随 \(k\) 增大，BLOM 在 locality / approximation / support width / runtime 上的变化。

必须展示：
- matching error vs \(k\)
- support width vs \(k\)
- runtime vs \(k\)
- Jacobian width vs \(k\)
- 可选：trajectory deformation vs \(k\)

要求：
- 支持 fixed \(q,T\)
- 支持多个 \(M\)
- 支持同时显示 raw Scheme C 与 ideal truncation

### 页面 5：Scheme A / B / C 对比
目标：
展示 continuity vs locality trade-off。

必须展示：
- lower-order jumps
- higher-order jumps
- support width
- outside-band sensitivity
- matching error（若已有）

要求：
- 把 Scheme A / B / C 清楚区分
- 明确指出 raw Scheme C 是主对象
- A/B 更多作为 trade-off 对照

### 页面 6：Phase 8 matching bridge 可视化
目标：
把 actual BLOM、reference-window、ideal truncation 的关系讲清楚。

必须展示：
- full matching error
- interior-only matching error
- boundary-layer matching error
- raw -> reference-window gap
- reference-window -> ideal-truncation gap
- empty-interior 标记（若 applicable）

要求：
- 这一页要尽量复用 Phase 8 的真实结果
- 不允许只画概念箭头
- 必须让用户看出：
  “主桥梁问题在哪里”

### 页面 7：Phase 9/10 optimizer demo
目标：
展示 BLOM 已经不是静态表示，而是可优化对象。

必须展示：
- objective curve
- 各子项（ctrl/time/obs/dyn/bc/reg）曲线
- gradient norm curve
- 轨迹 before / after
- 时间分配 before / after
- dense vs sparse backward consistency（可做小面板）

要求：
- 支持切换最小 objective / full objective
- 若 runtime 较慢，可使用缓存

### 页面 8：Benchmark / ablation summary
目标：
展示 framework-level 结果。

必须展示：
- runtime vs \(M\)
- objective gap vs \(k\)
- memory / support width vs \(k\)
- baseline compare：
  - raw Scheme C
  - MINCO
  - Scheme A
  - heuristic
- term ablation：
  - drop obs
  - drop dyn
  - drop time
  - drop reg

要求：
- 图要适合论文 / 组会直接使用
- 指标要写清楚定义
- 不要夸大结论

【必须支持的交互控件】

至少包含：
1. \(M\) slider
2. \(k\) slider
3. 扰动点索引选择
4. 扰动大小 slider
5. method selector（B-spline / MINCO / BLOM）
6. scheme selector（A/B/C）
7. regime selector（uniform / bounded-nonuniform）
8. matrix threshold slider（查看稀疏 mask）
9. objective term toggles
10. baseline compare toggles

【可视化要求】

1. 图风格要学术、简洁、白底、适合论文补充材料
2. 颜色统一且可区分，但不要花哨
3. 所有图必须有清楚标题、坐标轴、图例
4. 如果某图来自“示意模式”而不是真实求解结果，必须在标题或角落显式标注：
   - “schematic”
   - “illustrative only”
5. 优先输出 Plotly figure，保证可互动
6. 同时支持导出静态 PNG / PDF（如果容易实现）

【科学准确性要求】

请务必遵守：

1. 不要手工捏造 BLOM / MINCO 的影响范围来假装是真实求解结果。
2. 若某页必须用占位示意，请明确：
   - 这是示意，不是数值证据。
3. BLOM 局部性优势的主要证据必须来自：
   - Jacobian sparsity
   - support width
   - outside-band sensitivity
   - benchmark scaling
   而不是只靠曲线肉眼观察。
4. 不要把 B-spline 的 control point 扰动结果和 BLOM/MINCO 的 waypoint 扰动结果写成“完全同类量”。
5. 若某个结论目前前面 phases 尚未严格证明，只能在页面说明中写成：
   - observation
   - empirical result
   - numerically supported
   不能写成 theorem

【你输出时必须包含】

1. 完整 Python 代码
2. 所有文件内容
3. 一个 README，说明：
   - 如何运行
   - 每页展示什么
   - 哪些接口依赖前面 phases
   - 哪些页面是“严格结果模式”，哪些是“示意模式”
4. 如有必要，给出一份 mock data fallback 机制：
   - 当底层 BLOM 接口不可用时，demo 仍可打开
   - 但必须显式标注“mock mode”

【额外要求】

请尽量让这个 demo 同时适合：
1. 自己理解 BLOM
2. 给导师汇报
3. 论文补充材料
4. 组会讲解
5. 答辩时对比 B-spline / MINCO / BLOM

不要只生成一个漂亮外壳。请优先保证：
- 科学准确
- 结构清晰
- 能解释 BLOM 的真正优势和局限
- 和前面 phases 的数学对象严格对齐

现在请开始生成完整项目。