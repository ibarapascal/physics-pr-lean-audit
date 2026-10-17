## idaholab/moose#33307 Use log10 in relevant SCM closures and update transient regime calculation

- What changed (1–3 sentences, naming files and formulas):
  `modules/subchannel/src/scmclosures/SCMFrictionUpgradedChengTodreas.C`（由 Updated 重命名）将绕丝扫流系数改为 \(W_s=19-11\log_{10}(H/D)\)，`SCMHTCBorishanskii.C` 将基础换热项改为 \(24.15\log_{10}(-8.12+12.76P/D-3.65(P/D)^2)\)；摩擦和 `SCMMixingChengTodreas.C` 的流态选择/插值改用入口整体 Reynolds 数，摩擦幂律的局部 Reynolds 数则截断为 \(\max(Re_i,1)\)。`SCMHTCGnielinski.C` 将修改版 Gnielinski 分母截断至至少 0.1；`SubChannel1PhaseProblem.C` 增加整体 Reynolds 数计算并修正隐式横流估计的 gap 索引，`SCMTriAssemblyMeshGenerator.C` 增加棒间距/丝径检查，另有测试基准及验证绘图更新。下面选择三个闭合层不变量，不把横流求解器的完整守恒性或经验相关式的实验准确度一并宣称已覆盖。

- Invariants (mathematical statements; write them as checkable equalities/inequalities, give variables, sign conventions, domain of validity):

  **I1 — primary：Cheng–Todreas 闭合应使用同一个整体流态参数，并产生非负摩擦耗散。**
  设棒径 \(D>0\)、节距 \(P>D\)、绕丝导程 \(H>0\)、丝径 \(0<d_w\le P-D\)，局部流通面积 \(S_i>0\)、湿周 \(p_{w,i}>0\)、\(D_{h,i}=4S_i/p_{w,i}\)。以正轴向质量流为正，入口 \(m_i\ge0\)、\(M=\sum_i m_i>0\)、动力黏度 \(\mu_i>0\)，要求
  \[
  \bar\mu={\sum_i m_i\mu_i\over M},\quad
  Re_b={M D_{h,b}\over A_b\bar\mu},\quad
  R_L=320\,10^{P/D-1},\quad R_T=10^4\,10^{0.7(P/D-1)}.
  \]
  对固定几何和局部 \(r_i=\max(Re_i,1)\)，\(f_L=C_L/r_i\)、\(f_T=C_T/r_i^{0.18}\) 应满足
  \[
  f_i=\begin{cases}
  f_L&Re_b\le R_L,\\
  f_L(1-\psi)^{1/3}(1-\psi^7)+f_T\psi^{1/3}&R_L<Re_b<R_T,\\
  f_T&Re_b\ge R_T,
  \end{cases}\qquad
  \psi={\ln(Re_b/R_L)\over\ln(R_T/R_L)}.
  \]
  因而 \(f_i>0\)、端点连续，且 Darcy 摩擦项 \(-dp/dz|_f=f_i\rho u|u|/(2D_{h,i})\) 对应的耗散 \((-dp/dz|_f)u=f_i\rho|u|^3/(2D_{h,i})\ge0\)。不要求过渡式是凸组合或整体单调；其权重和并不恒为 1。

  为直接识别本 PR 的 log10 修正，选三角晶格 EDGE 子通道、壁面到棒表面的距离 \(g\)，取 \(x=g/D\in[0.1,0.3]\)。定义
  \[
  c_L=44.4+256.7x-267.6x^2,\quad c_T=0.1430+0.04199x-0.04428x^2,
  \]
  \[
  \tan\theta={\pi(D+d_w)\over H},\quad
  A_r={\pi(D+d_w)d_w\over4},\quad
  A_p=S_i+{\pi d_w^2\over8\cos\theta},\quad
  B=1+[19-11\log_{10}(H/D)]{A_r\over A_p}\tan^2\theta.
  \]
  应有 \(C_L=c_L B\)、\(C_T=c_T B^{1.41}\)。在下面选择的 \(8\le H/D\le52\) 中，\(W_s>0\)，故 \(B\ge1\)。这给出不依赖整机 gold CSV 的明确系数 oracle；其他子通道不套用 EDGE 系数。

  同一不变量的 **secondary 混合一致性检查**：对于有 CENTER 邻接的 gap，令 \(\beta_L,\beta_T\) 为固定几何在纯层流/纯湍流整体入口条件下返回的混合系数，则
  \[
  \beta=\beta_L+(\beta_T-\beta_L)\psi^{2/3},\quad
  0\le\beta_L\le\beta\le\beta_T
  \]
  在过渡区成立；对纯外围 EDGE/CORNER gap 的 sweep 系数同样成立。保持几何、入口 \(M,\bar\mu\) 不变而重新分配下游局部流量，混合系数不变；摩擦的流态权重不变，但摩擦值允许因局部 \(Re_i\) 改变。纯外围 gap 的普通 mixing 返回 0，涉及 CENTER 的 gap 的 sweep 返回 0，避免把两种机制相加。

  **I2 — primary：Borishanskii 换热的无量纲公式、低 Pe 平台及连接点一致性。**
  限定燃料棒、三角排列、\(p=P/D\in[1.1,1.5]\)，首先取局部 \(Re\ge R_T(p)\)，\(Pe=Re\,Pr\in[0,2200]\)。令
  \[
  F(p)=-8.12+12.76p-3.65p^2,\quad
  N_0(p)=24.15\log_{10}F(p),\quad
  a(p)=0.0174[1-e^{6-6p}].
  \]
  应有
  \[
  Nu=N_0(p)+a(p)[\max(Pe-200,0)]^{0.9}>0,\qquad h={kNu\over D_h}>0\quad(k,D_h>0).
  \]
  即 \(Nu(Pe_1)=Nu(Pe_2)=N_0\) 对 \(0\le Pe_1,Pe_2\le200\) 成立；\(Nu\) 在 \(Pe=200\) 连续，且对 \(Pe\) 非减。这里不要求该点可微，右导数有奇异性。将所有长度乘 \(s>0\) 并保持 \(Re,Pr\) 不变，应有 \(Nu'=Nu\)；若 \(k\) 不变，则 \(h'=h/s\)。

  对低流量的附属库一致性检查必须使用公共基类的过渡约定：\(Re\le R_L\) 时 \(Nu=N_L\)，CENTER/EDGE/CORNER 的 \(N_L\) 分别为 3.73/3.59/3.52；\(R_L<Re<R_T\) 时令 \(w=(Re-R_L)/(R_T-R_L)\)，应有 \(Nu=(1-w)N_L+wN_B(p,R_TPr)\)。不能直接拿低 Reynolds 数返回值和未经混合的 Borishanskii 湍流式比较。

  **I3 — secondary：修改版 Gnielinski 的被动换热、分母截断连续性及未截断区域一致性。**
  对局部 \(Re\ge R_T\)、\(Pr>0\)、有限正 Darcy 摩擦因子 \(f\)，定义
  \[
  a=f/8,\quad Q=Pr+0.01,\quad d=1+12.7\sqrt a(Q^{2/3}-1),\quad n=a(Re-1000)Q.
  \]
  应有 \(Nu=n/\max(d,0.1)\)，故 \(0<Nu\le10n\)、\(h=kNu/D_h>0\)；在 \(d=0.1\) 处连续，\(d\ge0.1\) 时与未截断的修改版公式相同。令热流正方向为壁面到冷却剂，\(q''=h(T_w-T_c)\)，对 \(T_w,T_c>0\) 应满足
  \[
  q''(1/T_c-1/T_w)={h(T_w-T_c)^2\over T_wT_c}\ge0.
  \]
  截断分母只给出固定输入下的有限性和上述依赖 \(n\) 的上界，不给出与 \(Re,f\) 无关的上界，也不证明 Nu 对 Reynolds 数呈平台。在库覆盖的过渡区另检查最终混合返回值 \(Nu\ge0\)；此处不能用 \(n\ge0\) 推理，因为 \(Re<1000\) 时分子可能为负。

- Sources of the invariants (first principles / literature formula / consistency with other implementations in the library / dimensions):
  摩擦耗散、被动传热熵产生和 \(hD_h/k=Nu\) 分别来自机械耗散、热力学第二定律和量纲定义。具体系数属于经验闭合约定，并非第一性原理定理；本次只通过 GitHub 合并版本核对其公开定义，未独立获取/核验原论文全文。
  I1 的流态参数以合并版 [摩擦文档](https://github.com/idaholab/moose/blob/b5a1d9e3df14e98850e514f3e57d02c25c44fca9/modules/subchannel/doc/content/source/scmclosures/SCMFrictionUpgradedChengTodreas.md) 为依据；该文档引用 Chen/Todreas 2018 与 Todreas/Kazimi，并明确使用简化插值。EDGE 系数及扫流项来自 [公开摩擦实现](https://github.com/idaholab/moose/blob/b5a1d9e3df14e98850e514f3e57d02c25c44fca9/modules/subchannel/src/scmclosures/SCMFrictionUpgradedChengTodreas.C)，混合约定见 [混合实现](https://github.com/idaholab/moose/blob/b5a1d9e3df14e98850e514f3e57d02c25c44fca9/modules/subchannel/src/scmclosures/SCMMixingChengTodreas.C)。因此系数等式是实现契约审计，不能单凭和源码一致就断言经验模型正确。
  I2 的独立文档公式见 [Borishanskii 页面](https://github.com/idaholab/moose/blob/b5a1d9e3df14e98850e514f3e57d02c25c44fca9/modules/subchannel/doc/content/source/scmclosures/SCMHTCBorishanskii.md)，引用 Borishanskii 1969；层流值和过渡插值来自 [HTC 基类文档](https://github.com/idaholab/moose/blob/b5a1d9e3df14e98850e514f3e57d02c25c44fca9/modules/subchannel/doc/content/source/scmclosures/SCMHTCClosureBase.md)。
  I3 的移位 Pr 和 0.1 下限是 [本库修改版 Gnielinski 约定](https://github.com/idaholab/moose/blob/b5a1d9e3df14e98850e514f3e57d02c25c44fca9/modules/subchannel/doc/content/source/scmclosures/SCMHTCGnielinski.md)，不是把原始 Gnielinski 1975 公式或低 Pr 实验验证直接移植到当前修正；物理必要条件是非负传热，不是特定的 0.1 数字。

- Common operating domain (range of each input, with the basis):
  共同测试环境为单相、固定几何、有限正密度/黏度/导热率、非反向轴向流、无入口发展效应；不模拟沸腾、反向流或变形。闭合 API 的实际自变量是无量纲 \(Re,Pr\)、几何以及已初始化的通道状态，故无需虚构普适压力/温度范围；若使用流体属性对象，应另限制在该对象有效区间。
  I1 推荐三角绕丝束 \(P/D\in[1.1,1.35]\)、\(H/D\in[8,52]\)、\(d_w/D\in[0.02,\min(0.2,P/D-1)]\)、\(g/D\in[0.1,0.3]\)，7、19、37、61、91、127、169 或 217 根棒；必须由 mesh generator 生成所有正面积和湿周，不能任意独立指定互不相容的几何。前两区间及棒数采用摩擦/混合文档与警告范围的保守交集（摩擦源码注释与警告的导程/棒数界限存在差异，因此不采用争议边界）；丝径和壁隙下限是审计采样选择，不是新宣称的实验适用界限。取 \(D\in[0.003,0.02]\) m 作为数值尺度，\(Re_b,Re_i\in[50,10^6]\)，混合的常规采样要求局部 gap Reynolds 数也在 \([400,10^6]\)。
  I2 使用裸棒三角束（基类适用表列出的几何），\(P/D\in[1.1,1.5]\)、\(Re\in[R_T,10^6]\)、\(0<Pe\le2200\)，由 \(Pr=Pe/Re\) 构造输入；\(Pe=0\) 仅为数学边界/极限。层流及过渡样本单独按库混合规则检查。I3 先用 \(P/D\in[1.1,1.35]\)、\(Re\in[R_T,10^6]\)、\(Pr\in[10^{-5},2000]\) 检查库声称的数值覆盖范围；该范围不能等同于原始相关式的实验验证范围，尤其低 Pr 与绕丝场景只作库扩展的健全性检查。附加过渡样本为 \(Re\in[R_L,R_T]\)。各通道/棒索引必须有效，混合 API 的轴向索引 \(iz\ge1\)。

- Numerical check method (which public API of the merged commit, what to compute, what to compare against, suggested tolerance and why; if finite differences are needed, the step-size strategy):
  以下均为后续执行方案，本次没有编译或运行。使用 SubChannelApp 正常的 `[TriSubChannelMesh]`、`[SubChannel]`、`[SCMClosures]` 输入接口创建对象，通过正常初始化和 `externalSolve()` 形成入口状态；在后续 C++ 审计驱动/测试对象中调用公开方法，不修改 protected 状态、不开私有访问漏洞。对象需要真实 MOOSE problem/mesh 上下文，并非可以直接默认构造的纯函数。公开类型与 getter 见 [SubChannel1PhaseProblem.h](https://github.com/idaholab/moose/blob/b5a1d9e3df14e98850e514f3e57d02c25c44fca9/modules/subchannel/include/problems/SubChannel1PhaseProblem.h) 和 [SCMHTCClosureBase.h](https://github.com/idaholab/moose/blob/b5a1d9e3df14e98850e514f3e57d02c25c44fca9/modules/subchannel/include/scmclosures/SCMHTCClosureBase.h)。先用单 MPI rank，避免把 rank-local 状态与闭合公式问题混在一起。

  I1：调用 `getBulkReynoldsNumber()`，与入口每通道 `mdot`、`mu` 及 mesh 的 `getAssemblyHydraulicDiameter()`、`getAssemblyFlowArea()` 所得参考值比较。调用 `getFrictionClosure()->computeFrictionFactor(FrictionStruct(i_ch, Re_i, S_i, w_perim_i))`，对 EDGE 样本以以上显式公式独立计算全部参考值。固定同一局部参数，改变输入入口质量流以使 \(Re_b=R_L(R_T/R_L)^t\)，取 \(t=0,0.01,0.1,0.5,0.9,0.99,1\) 及区间外点；另构造局部流态与整体流态不同的样本，防止旧的局部 Re 选支逻辑也通过检查。通过命名 UserObject 的公共获取接口取得 `SCMMixingChengTodreas`，调用 `computeMixingParameter(gap, iz)` 和 `computeSweepFlowMixingParameter(gap, iz)`，检查同一 \(t\) 下的 \(t^{2/3}\) 插值和 gap 类型限制。流量重分配测试应维持总入口流量与加权黏度；使用等温、同一黏度入口可简化构造。

  I2/I3：从 `getPinHTCClosure()`（适用时也可用 `getDuctHTCClosure()`）取得 HTC 对象，调用 `computeNusseltNumber(FrictionStruct(...), NusseltStruct(Re, Pr, i_pin, iz, i_ch))` 及 `computeHTC(..., k)`。两个结构中的 Reynolds 数采用一致的物理状态；Borishanskii 使用有效棒索引，不用 duct sentinel；若有 `Dpin` 辅助变量，应将其正确初始化，使有效棒径与参考几何一致。I2 以对数均匀采样的 Re、若干 p，以及 \(Pe=1,50,199,200,201,1000,2200\) 检查公式；以 \(s=0.5,2,10\) 检查几何相似性，保持无量纲输入不变。

  I3 先调用实际摩擦 closure 获得 f，再计算参考分母和分子；不把任意虚构的 f 当作已暴露的 HTC 输入参数。扫描入口流量、通道种类、几何和 Pr 寻找 \(d<0.1\)、\(d=0.1\) 两侧样本；若某固定 f 能达到截断边界，可用
  \(Pr_*=[1-0.9/(12.7\sqrt{f/8})]^{3/2}-0.01\)
  定位边界，但仅在括号为正且 \(Pr_*\) 落在采样域内时使用。对 \(Re\ge R_T\) 检查精确截断公式、有限性和 \(Nu\le10n\)；对过渡区检查实际 \(Nu\ge0\)，并按源码使用局部 Re 的 Gnielinski 分子再作基类混合，不能误套 Borishanskii 的 \(R_T\) 端点评价策略。若合法输入导致 `computeHTC` 报负值错误，记为域内失败，不能将“成功拒绝负 HTC”视为物理通过。

  双精度直接闭合比较建议 \(|y-y_*|\le10^{-12}+10^{-10}|y_*|\)（均为无量纲量）；HTC 用 \(hD_h/k\) 比较以避免有量纲绝对容差。比较整体 Re 时按已读取的实际入口状态重算，以避免求解误差污染公式测试。若只能通过整机 CSV 检查，建议相对 \(10^{-5}\) 起步并收紧求解容差验证收敛，不用该宽松阈值替代闭合单点检查。
  无须有限差分求导。端点用相对扰动 \(10^{-3},10^{-5},10^{-7},10^{-9}\) 序列检查一侧极限并逐点与解析参考比较；不能要求端点误差按扰动线性缩小（摩擦含立方根、混合含 \(2/3\) 次幂，Borishanskii 在 Pe=200 右侧含 \(0.9\) 次幂）。随机采样同时记录种子与输入，单独报告每个约束。
  另作 I1 的域外数值边界诊断：逐渐将所有入口质量流趋于零，再试恰好为零。物理静止极限应有 \(Re_b\to0\) 且摩擦压降趋零；局部 \(\max(Re_i,1)\) 是数值正则化，不要求 \(Re_i<1\) 时还严格服从层流线性阻力定律。零流量不计入上述经验相关式的常规适用域，失败需单独标注。

- Do the PR's own tests already check this invariant (yes/partial/no, naming test files and assertions):
  **I1：partial。** `modules/subchannel/test/tests/problems/SFR/sodium-19pin/tests` 中 `test_explicit`、`test_implicit`、`test_monolithic` 是温度等输出的 CSVDiff（implicit 的 `rel_err=3e-3`，monolithic 的 `rel_err=7e-6`）；`test_implicit_tiny_dt` 等为 RunApp 收敛检查。它们间接使用修改后的闭合，但未直接断言整体/局部 Re 分离、扫流 log10 系数或端点恒等式。输入 diff 显示 `test/tests/scmclosures/friction_updated_cheng_todreas/tests` 及其非法 corner 输入被删除，不能把该测试算作合并版覆盖。
  **I2：partial。** [heat_transfer_correlations/tests](https://github.com/idaholab/moose/blob/b5a1d9e3df14e98850e514f3e57d02c25c44fca9/modules/subchannel/test/tests/problems/heat_transfer_correlations/tests) 中 `EBR-II_SHRT-17_SS/borishanskii` 用 `XX09_SCM_SS17.i` 比对 `borishanskii.csv`、`abs_zero=1e-5`；输入的 postprocessor 为三个棒表面温度位置。PR 更新了 gold，不存在由该 CSV 本身独立证明 log10 公式、Pe 平台或尺度不变性的断言。
  **I3：partial，截断边界/低 Pr 数学性质的专门检查为 no。** 同目录 `gnielinski` 是 CSV 回归；`negative_htc_error` 使用 `SCMHTCNegativeTest` 并匹配 `The heat transfer coefficient must be non-negative and finite.`，它验证通用防护而非 Gnielinski 域内始终非负。`invalid_wire_diameter_error` 匹配 `The wire diameter must not exceed the pin-to-pin gap`，只验证几何拒绝。sodium-19pin 的 duct closure 已改为 Dittus–Boelter，不能将这些运行成功作为低 Pr Gnielinski 截断的验证。以上是对本次变更及所查看测试文件的覆盖判断，不是整个仓库测试全集的穷尽声明。

- Do you expect the change to satisfy it (satisfies / may not satisfy / uncertain), with one or two sentences of reasoning (speculation, not a conclusion):
  **I1：satisfies（正流常规域，推测）；零流边界 may not satisfy。** 正流分支、十进制扫流项及混合公式与上述契约相符，但 `computeBulkReynoldsNumber()` 仍计算 \(\sum m_i\mu_i/M\)，全零流量会出现 0/0，局部 Re 截断并不能消除这条路径；是否在具体公共输入流程触发需后续数值验证。
  **I2：satisfies（推测）。** log10 修改与文档公式一致，正参数下 Pe 项的符号及分段连接正确；这不证明整个棒束经验模型或更新的 gold 对实验正确。
  **I3：satisfies（限定湍流代数条件，推测），过渡区及低 Pr 物理准确度 uncertain。** 下限为正且湍流分子为正可支持非负性，但仅分母截断不足以证明文档所说的低 Pr 平台，也不能单独证明过渡区的最终非负性。

- Environment: merge commit SHA, language and dependencies needed, estimated setup difficulty (low/medium/high):
  Merge SHA：`b5a1d9e3df14e98850e514f3e57d02c25c44fca9`；输入 PR URL：[idaholab/moose#33307](https://github.com/idaholab/moose/pull/33307)。实现为 C++，需要匹配该版本的 MOOSE framework、SubChannel 模块、流体属性依赖、libMesh、PETSc、MPI 及其构建工具；Python 用于 MOOSE 测试 harness，CSV-only 审计无需生成验证图。`modules/subchannel/Makefile` 通过 framework/build.mk、framework/moose.mk 和 modules/modules.mk 构建应用；新增薄 C++ 审计驱动应链接相同模块而非复制闭合源码。预计首次配置难度 **high**（完整 MOOSE 工具链及对象初始化），已有兼容环境时为 medium。本次只读 `./input` 与 GitHub 合并版本源文件，未查看其他本地文件、另一审计轨道、PR 评论/issue/fork，也未克隆、安装、编译或执行物理代码。

- Formalization (Lean) feasibility: can the invariant be written as a proposition over finite-dimensional reals, and what are the main obstacles:
  **可以。** 固定有限个通道和 gap，将几何、入口流量和物性表示为有限维实向量；I1 是加权平均、正系数分段幂律、共同插值参数与连续性命题；I2 是含 log/exp/实数幂的分段函数相等、单调性及尺度不变性；I3 是 max 截断函数的正性、连续性、上界和熵产生恒等式。主要障碍是 log/sqrt/非整数幂的域条件、网格几何保证正面积、分段边界、将 C++ 浮点与实数模型连接，以及证明正常初始化后的整体 Re 状态与各 closure 使用值一致。经验系数的实验有效性和低 Pr 改式的可信度不能通过实数代数证明；NaN、错误处理、MPI 状态与迭代收敛也不属于单纯有限维实数恒等式证明，需另设运行语义和数值误差层。
