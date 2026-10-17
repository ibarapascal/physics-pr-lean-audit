# MOOSE#33307 ΔQ1 抽查实施日志

> 目的：记录三路检查计划逐句到 `spotcheck.py` / `scm_port.py` 的映射与全部实施判断，供协议 §6/§7 审核。
> 结构：§0 环境（非原生）→ §1 第一路 → §2 第二路 → §3 第三路 → §4 耗时与内存 → §5 偏离候选 → §6 冒烟后修正 → §7 待裁定（冒烟拒绝的归类）。
> 更新规则：冒烟后只许修脚本错误；改计划映射须记入协议 §7，并在 §5 追加。
> 引用：上游 协议、track1.md/track2.md 的 #33307 节、`specs_codex/idaholab_moose_33307/spec.md` §2/4/5、合并提交 b5a1d9e 源码；下游 `dq1/out/mo33307_t{1,2,3}.json`。

## §0 环境：非原生重算（偏离候选）

- compute node A 与 compute node B 只读检查：`env/` 下无 MOOSE/libMesh/PETSc，`subchannel-opt`、`mpiexec` 均不在 PATH，家目录 4 层内无 moose 目录。按指示不构建 MOOSE。
- 改为按合并提交 b5a1d9e 的 C++ 源码在 `scm_port.py` 中逐行重算被改函数（文件头列出覆盖清单）：
  - 被改：`SCMFrictionUpgradedChengTodreas`（构造检查/警告、三角与方栅格摩擦因子）、`SCMMixingChengTodreas`（构造警告、`computeMixingParameter`、`computeSweepFlowMixingParameter`）、`SCMHTCBorishanskii`、`SCMHTCGnielinski`、`SubChannel1PhaseProblem::computeBulkReynoldsNumber`（单进程语义）、`SCMTriAssemblyMeshGenerator` 的 nrings/pitch/dwire 检查（`absoluteFuzzyGreaterThan` 默认容差 libMesh::TOLERANCE² = 1e-12，已读 framework `MooseUtils.h`）。
  - 未改但被调用：`SCMHTCClosureBase`（preinfo、blend、computeHTC）、`TriSubChannelMesh`/`QuadSubChannelMesh` 的子通道面积、湿周、组件水力参数。
  - 不能重算：`implicitPetscSolve` 横流估计索引（#33340，需 PETSc 求解器上下文）、`externalSolve` 全流程、MPI（`processor_id()!=0` 分支）。
- `mooseError`→`MooseError` 异常（计失败）；`flagSolutionWarning`→警告列表。
- 运行环境：任意 Python 3.10+＋numpy；JobSpec 用 `env/py2/bin/python`。
- 种子 20261017＋9×10＋路号。几何默认值取 `test/tests/problems/SFR/sodium-19pin/test19_implicit.i`（nrings=3、pitch=0.012、D=0.01、dwire=0.002、hwire=0.0833、flat_to_flat=0.056 → 管壁间隙 0.0022154 m）；方栅格侧隙比取 `verification/friction_model_verification/two_channel2.i`（side_gap/D=0.1）。
- 子通道局部 Re（计划未给区间）：按均匀质量通量取 Re_i = Re_b·D_h,i/D_h,b（第一、二路）；第三路计划要求局部与组件 Re 各自独立抽样。

## §1 第一路（track1.md）

分支：子通道类型 CENTER/EDGE/CORNER（k=3，每支 ⌊N/3⌋）。

| 计划原文 | 代码 |
|---|---|
| 常用工况域 P/D∈[1.05,1.42]、H/D∈[8,52]、D_w/D∈[0.05,0.25]、N_pin∈{7,19,37,61,217}、Re_b∈[50,1e6]、Pr∈[0.003,0.03]、Pe∈[10,2200] | `t1_sample`：连续量线性均匀；N_pin 在 5 个值中均匀取 |
| A：H/D 取 8、20、52，Re_b≈5×10⁴，带绕丝；`wire_diameter`=1e-6 m 跑近光棒；ΔP_wire ≥ ΔP_bare 且随 H/D 减小单调增大；严格不等式 | 闭合层替代：`t1_A_f_wire>=f_nearbare(Dw=1e-6m)_same_mdot`（同入口质量流量换算 Re_b、Re_i）、`t1_A_f_monotone_increasing_as_HoD_decreases(8,20,52)`；必测点 3 类型×H/D{8,20,52}×Re_b=5e4 |
| A（不变量）：修正因子 [1+Ws_T(A_r/A'_p)tan²θ]^1.41 ≥ 1、层流因子 ≥ 1 ⇔ Ws_T>0（边/角） | `t1_A_WsT>0`、`t1_A_turbulent/laminar_wire_factor>=1`（只在 EDGE/CORNER 支） |
| B：Re_b=Re_L(1±1e-3)、(1±1e-4)，ΔP 相对差 <10δ；Re_T 同理 | 闭合层替代：固定 Re_i，`t1_B_f_jump_at_ReL/ReT_delta=…<10delta` |
| B（不变量）：ψ=0 时 f=f_L(Re_i)，ψ=1 时 f=f_T(Re_i)；Cm、Cs 同理；所有子通道共享同一个 ψ | `t1_B_f(Reb=ReL)==fL`、`…==fT`；Cm（CENTER 支，center–center 间隙）/Cs（EDGE 支 edge–edge、CORNER 支 corner–edge）端点与跳变；`t1_B_same_psi_all_channel_types` |
| C：`pin_HTC_closure=borishanskii`，低流量 Pe<200，反推 Nu 与 24.15·log₁₀(poly) 比，相对 1e-6 | Pe≤200 的抽样点：Re = Pe/Pr，直接读 `computeNusseltNumber`；`t1_C_…`；P/D∈[1.1,1.5] 时另检 Nu∈[4.25,10.83] |
| D：explicit/implicit/monolithic 三版比较 | 不能重算，列入 `not_covered` |
| E（不变量）：Re>10³、Pr>0 时 Nu_T 有限且 >0，关于 Pr 连续 | Re_i>1e3 的点：`t1_E_Gnielinski_NuT_finite>0`、`…_Nu_finite>0`、`t1_E_Nu_continuous_in_Pr_delta=…<10delta`（摩擦因子取本闭合） |
| 建议：1 进程与 4 进程一致 | 不能重算，`not_covered` |

实施判断：
1. D_w/D 上限 0.25 大于 P/D−1 下限 0.05：D_w > P−D 的抽样点（约 27%）在网格生成器检查处报错，按原生行为计「失败」。
2. 近光棒与 H/D 换档时保持同一入口质量流量与粘度：Re_b' = Re_b·(A/D_h)·(D_h'/A')，Re_i' = Re_b'·D_h,i'/D_h,b'。
3. B 的「ΔP 相对差」换成摩擦因子相对差；E 的「连续」没有判据，借用同路 B 的 <10δ（δ=1e-3、1e-4）。
4. 混合闭合的局部场（只影响局部 Re 警告）取均匀质量通量、μ=2.5e-4 Pa·s。
5. 必测点里 Pe 取 1000（不触发 C），Pr=0.005，P/D=1.2、D_w/D=0.2、N_pin=19（模板值）。

## §2 第二路（track2.md）

分支：CENTER/EDGE/CORNER（k=3）。

| 计划原文 | 代码 |
|---|---|
| 工况域 P/D 1.0–1.42、H/D 8–52、N_pin 7–217、Re_b 50–1e6；液态金属 Pr 0.004–0.03 与水；Borishanskii 1.1≤P/D≤1.5、Pe 200–2000 | `t2_sample`：N_pin 取 7–217 内全部合法六角数 {7,19,37,61,91,127,169,217}；Borishanskii 的 P/D 单独抽 |
| I1：扫 H/D ∈ {8,15,25,52}，有/无绕丝（d_wire=0），压降有绕丝 ≥ 无绕丝；不出现 turbulent wire correction 报错 | 闭合层：`t2_I1_f_wire>=f_bare_same_Re`（边/角）、`t2_I1_no_turbulent_wire_correction_error`；必测点 3 类型×H/D{8,15,25,52} |
| I1（不变量）：Ws_T ∈ [0.12, 9.07] > 0；边/角修正因子 ≥1；中心 f_wire = f_bare·(P_w,p/P_w) + 正的拖曳项 | `t2_I1_WsT_*`、`t2_I1_*_wire_factor>=1`、`t2_I1_center_drag_terms>0`、`t2_I1_center_CfT==…` |
| I2：Re_b 在 Re_L、Re_T 两侧 ±0.1%，跳变 ≤1e-3 相对 | `t2_I2_f_jump_at_*_pm0.1pct<=1e-3`；Cm/Cs 同样 |
| I2（不变量）：ψ=0→f_L、ψ=1→f_T；f≥0；Re_L、Re_T 公式；同层同 ψ；Cm、Cs 端点连续 | `t2_I2_*` |
| I3：三种求解比较；单进程与 mpiexec -n 2 | 不能重算，`not_covered` |
| I3（不变量）：Pe≤200 时 Nu=24.15·log₁₀(…)，P/D=1.1 ≈4.25、1.5 ≈10.8 | 抽样点在 Pe=200（本路 Pe 区间下端）评估：`t2_I3_Borishanskii_Pe=200_…`；必测点 P/D∈{1.1,1.5}×Pr∈{0.004,0.03}×类型，容差为末位半个单位（0.005、0.05） |

实施判断：
1. 计划未给 D_w：取模板 test19 的关系 D_w = P−D（模板 0.002 = 0.012−0.01）；这样网格检查恒通过。
2. 「无绕丝（d_wire=0）」：只设 d_wire=0 会在摩擦闭合构造函数报错（要求两者同为 0），取 d_wire=hwire=0 的光棒；与有绕丝在同一 (Re_i, Re_b) 比较（计划式写的是 f(Re) 对 f(Re)）。
3. 计划内部不一致：Borishanskii 工况域 Pe 200–2000，而 I3 的式子与核对方法针对 Pe≤200（「设低 Pe」）。处理：Pe 照工况域抽样（只记录），Pe≤200 的式子在区间端点 Pe=200 上检查，Re = 200/Pr。
4. 水：工况域只写「与水」，没给 Pr 区间；Borishanskii 只用于液态金属，未加水分支。
5. 局部 Re 同第一路（均匀质量通量）。

## §3 第三路（spec.md §2/§4/§5）

分支（k=20，每支 500 点）：摩擦 9 支（三角光棒/三角绕丝/方栅格 × CENTER/EDGE/CORNER）、混合 2 支（Cm：center 参与的间隙；Cs：边/角间隙）、Borishanskii 3 支（类型）、Gnielinski 6 支（三角绕丝/方栅格 × 类型）。

| 计划原文 | 代码 |
|---|---|
| 「编一个小测试应用，调用公开 C++ 接口；不要用 Python 转写替代闭合」 | 不可行（无构建）：闭合按源码逐行重算，独立 oracle 另按文档公式＋系数表写（`o_*` 函数）；容差照计划 \|y−y_ref\| ≤ 1e-12+1e-10\|y_ref\| |
| I1：摩擦 CENTER/EDGE/CORNER；三角光棒与绕丝；表切换 1.1 两侧与相等（偏移 1e-6）；H/D=8,20,52；d_w/(P−D)=0.1,0.5,1；R_b ∈ {0.5R_L, R_L, √(R_LR_T), R_T, 2R_T} × R_i ∈ {0,0.1,1,50,10³,10⁵} | 抽样支 `F_*`；必测点 `grid_*`（810）、`switch_*`（27） |
| I1：端点值直接核；同几何同 R_b 改局部流量不改 ψ；过渡区 ψ∈[0,1]；不断言凸组合 | `t3_I1_f(Reb=ReL)==fL`、`…ReT…`、`t3_I1_psi_independent_of_local_Re`、`t3_I1_psi_in_[0,1]` |
| I1：独立累加入口 μ_i·m_i 比 R_b；均匀/非均匀粘度与流量分配在同 (M, μ_b) | 三角支：随机分配＋随机粘度算 R_b 与独立式比；另一分配＋均匀粘度应得同一 R_b |
| I1 混合：两个函数都调；中心/边/角组合；零分支；端点界；改局部场不改；环数 2,3,9，单环另作探针；Cm、Cs 常数 | `M_*` 支：`t3_I1_mixing_beta_vs_oracle`、`CL<=C<=CT`、端点、`nondecreasing`（R_b×1.001）、`other_mixing_function_exact_zero`、`independent_of_local_fields`；必测点 `mix_*`（60）、`mixing_one_ring_probe` |
| I2 Borishanskii：r=1.1,1.2,1.35,1.5；R_i>R_T 固定；Pe=1,199.999,200,200.001,500,2200；比绝对值；Pe>2200 只看警告；作 duct 闭合应报错；k 的量纲检查；长度缩放 0.1、10；名义与变形 Dpin；局部 R_i 扫过 R_L、R_T 端点 | `B_*` 支（r、Pe 均匀，R_i∈(R_T,1e6] 均匀，Pr=Pe/R_i）＋必测点 `boris_*`（84，含 Pe=2500 警告探针）、`boris_as_duct_closure` |
| I3 Gnielinski：pin 与 duct 哨兵；三角与方栅格；Pr = 1e-5,0.001,0.005,0.01,0.5,1,10,2000；层流/过渡/湍流局部 Re 与独立组件区段；下限门槛 f_c(1±1e-6) 与 f_c；R_L<R_i<1000 的被动性 | `G_*` 支（Pr∈[1e-5,2000] 线性均匀）：oracle、有限、Nu≥0（含混合区）、R_i≥max(R_T,1000) 时上界、未截断时等于原式、截断时有警告、h·D_h/k=Nu 或负值被 computeHTC 拒绝、duct 哨兵；必测点 `gn_*`（144）、`gn_floor_s=…`（反解 R_i 使 f=f_c(1+s)） |
| I4：横流估计关联矩阵/质量守恒（需跑 subchannel-opt） | 不能重算，`not_covered` |
| I5：P/D=1, 1−1e-6, 1+1e-6, 1.1 × d_w=0, P−D, P−D+0.5ε_g, P−D+2ε_g；先读 ε_g；尺寸缩放一例 | 必测点 `geom_*`（D=0.01 与 0.1 各 16 例），ε_g=1e-12 |
| I6：局部 R_i=0,1e-12,0.1,1,1+1e-6 的有限平台；入口整体缩放 t=1,1e-3,1e-6,0 | 每个摩擦抽样点检平台；必测点 `bulk_scale_t=…` |

实施判断：
1. 抽样区间：摩擦 P/D、W/D∈[1.01,1.4]（计划「表切换两侧」用的数学检查区间）、H/D∈[8,52]（两个来源的保守交集）、N_pin 7–217（交集），D∈[0.005,0.015] m，d_w/(P−D)∈(0,1]，R_b、R_i 各自 [50,1e6]；混合 P/D∈[1.067,1.35]、H/D∈[4,52]、R_b∈[400,1e6]、环数∈{2,3,9}；Borishanskii r∈[1.1,1.5]、Pe∈[0,2200]；Gnielinski Pr 线性均匀（计划未写对数均匀）。
2. 方栅格侧隙未给：取 two_channel2.i 的 side_gap/D=0.1；Borishanskii 所需面积/湿周用三角光棒几何、管壁间隙比取模板值。
3. 变形 Dpin 计划未给值：取 0.98·D。
4. Gnielinski 的摩擦闭合用三角绕丝（模板 test19 配置）或方栅格光棒。
5. R_b 单调性用 R_b 与 1.001·R_b 两点比较（计划没给步长）。
6. t=0 时 computeBulkReynoldsNumber 的 0/0 在 C++ 中得 NaN（不报错）；重算同样得 NaN → 按协议 §4 计「失败」，在报告中单列。
7. 计划中「transient inlet changes」「serial vs 2 ranks」「zero-flow channel in implicit solver」属解算器行为，未覆盖。

## §4 预计耗时与峰值内存（推测，未运行）

| 路 | 单点 | 10⁴ 点 | 依据 |
|---|---|---|---|
| t1 | ~1–2 ms（约 30 次闭合调用） | ~20 s | 纯 Python 标量 |
| t2 | ~1 ms | ~10 s | |
| t3 | ~1–3 ms | ~30 s＋必测约 1.2k 点 | |

峰值内存 <200 MB（只有 numpy）；JobSpec 申报 512 MiB。

## §5 偏离候选

- **非原生**：节点无 MOOSE 构建，被改函数按合并提交 C++ 源码在 Python 中重算（覆盖清单见 §0）；依赖完整解算器的检查（组件压降、三种求解一致、MPI、横流索引）未测。
- 第二路计划内部不一致（Pe 区间 vs Pe≤200 式子）的处理见 §2-3。

## §6 冒烟后修正（2026-10-17，N=100 冒烟之后；只修脚本错误，计划映射不变）

- `spotcheck.py` 第一路 B（约第 141–148 行）、第二路 I2（约第 239–245 行）：跳变检查里「固定的 Re_i」原先取 `re_local(g, t, Reb)`（抽样 Re_b 下的局部 Re，可达 1e6），却在 Re_b=Re_L(1±δ)≈300–800 处求 f，造成 Re_b≈Re_L 而 Re_i≫Re_L 的不自洽状态，与本日志 §0「Re_i = Re_b·D_h,i/D_h,b（均匀质量通量）」及两路计划「调入口质量流量使 Re_b = Re_L(1±δ)」（第一路 B、第二路 I2 数值核对方法）所描述的状态不符（比较对象错）。改为在每个端点 R0 上取 `Rei0 = re_local(g, t, R0)`，在 ±δ 两侧保持固定（「固定 Re_i」的映射不变）。端点等式检查 `f(Reb=ReL/ReT)==fL/fT` 同步用 `Rei0`（恒等式，结论不受影响）。
- 预期：首轮跳变残差（相对 0.3–28）主要来自该不自洽；修正后量级降到约 1e-2（δ=1e-3）/1e-3 量级（δ=1e-4），但 Re_L 处仍会拒绝，原因见 §7-1（计划本身）。验证交 `job_smoke2.json`（只跑 t1、t2，输出 `dq1/out/smoke2/`）。
- 「t2 拒绝 118 > n=100」不是计数错误：`total_rejected` = 抽样支 100 ＋ 必测点 18（必测 24 点：12 个摩擦必测点全在 Re_L 跳变处拒绝，6 个 Borishanskii Pr=0.03 必测点见 §7-2）。第一路 82 = 73 ＋ 必测 9 同理。未改。

## §7 待裁定（冒烟拒绝归类，脚本不改）

1. **(P) 第一路 B「ΔP 相对差 <10δ、按 δ 线性缩小」、第二路 I2「±0.1% 跳变 ≤1e-3」在 Re_L 处与库的插值式矛盾**。证据：首轮所有点 δ=1e-3 与 δ=1e-4 的跳变比恒为 2.15≈10^{1/3}，即跳变 ∝ δ^{1/3}；合并提交 `SCMFrictionUpgradedChengTodreas` 过渡区 f = f_L(1−ψ)^{1/3}(1−ψ⁷) + f_T ψ^{1/3}（第二路计划 I2 自己也写了此式），ψ = ln(Re_b/Re_L)/ln(Re_T/Re_L) ≈ δ/ln(Re_T/Re_L)，故 Re_L 右侧多出 f_T·ψ^{1/3} 项：函数连续（端点等式全过）但不是 Lipschitz，δ→0 时跳变按 δ^{1/3} 而非 δ 收敛。手算（CENTER 光棒系数、P/D=1.2、Re_i=Re_L，推测量级）：f_L≈0.195、f_T≈0.049、ψ(δ=1e-3)≈3.0e-4 → 跳变≈0.015＞0.01；δ=1e-4 → ≈0.008＞0.001。Re_T 侧 (1−ψ)^{1/3}(1−ψ⁷)≈7(1−ψ)^{4/3}，可微，首轮 Re_T 跳变 0 拒绝，与此一致。计划写错之处：第一路「δ→0 时线性收敛」「< 10·δ」，第二路「跳变 ≤1e-3 相对」（在 ±0.1% 处）。
2. **(P) 第二路 I2 的 Cm、Cs「±0.1% 跳变 ≤1e-3」**：混合参数 C = C_L + (C_T−C_L)ψ^{2/3}（第二路计划原文），δ=1e-3 时 ψ^{2/3}≈4.5e-3，首轮观测相对跳变 0.0036–0.0037（全部 Re_L 处），超 1e-3。与 1 同类：连续但按 δ^{2/3} 收敛。第一路容差 10δ 下 Cm/Cs 不拒绝，也与此一致。
3. **(P) 第一路 C、第二路 I3「Pe ≤ 200 时 Nu = 24.15·log₁₀(poly)」**：首轮第一路 6/8、第二路抽样 68/112、必测 6/12 拒绝（观测 Nu 低于式值，如 P/D=1.418、Pr=0.0247 时 7.22 对 10.16）。原因：合并提交未改动的 `SCMHTCClosureBase` 对局部 Re 分段混合——Re ≤ Re_L 取层流 Nu（3.52–3.73），Re_L<Re<Re_T 按 (Re−Re_L)/(Re_T−Re_L) 线性混合，并在过渡区用 Re_T 计湍流 Pe（Re_T·Pr 可 >200，再加 Pe 修正项）。计划给的 Pe（≤200）与 Pr（0.003–0.03）使 Re = Pe/Pr ∈ [333, 6.7e4]，大量落在 Re < Re_T（P/D 1.1–1.5 时 Re_T≈1.2e4–2.2e4），这时闭合返回的 Nu（也是全解算器里反推得到的 Nu）不是 24.15·log₁₀(poly)。只有 Re ≥ Re_T 的点（Pe ≤ 200 时不加修正项、混合权为 1）等式精确成立，首轮通过的点即此类（必测 Pr=0.004 → Re=5e4 全过，Pr=0.03 → Re≈6.7e3 全拒，与此一致）。计划写错之处：第一路 C、第二路 I3 把 NuT 的式子当作闭合整体输出，漏了库的层流/过渡混合（第二路连带「Nu∈[4.25,10.83]」的区间检查 3 次拒绝，同因）。
4. **(P) 第一路 27 次失败（非拒绝）**：全部为 `MooseError: The wire diameter must not exceed the pin-to-pin gap`。计划常用工况域把 P/D∈[1.05,1.42] 与 D_w/D∈[0.05,0.25] 独立给出，而本 PR 新加的网格生成器检查要求 D_w ≤ P−D；工况域的矩形约 27% 落在 PR 明确禁止的几何上。按协议 §4 计失败、不计违反；是否在 10⁴ 全跑时把这部分从分母中剔除，交主会话裁定（脚本未改）。
5. 本 PR 首轮未见 (V) 类。

