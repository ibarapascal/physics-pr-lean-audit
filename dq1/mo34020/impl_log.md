# MOOSE#34020 ΔQ1 抽查实施日志

> 目的：记录三路检查计划逐句到 `spotcheck.py` / `mortar_port.py` 的映射与全部实施判断，供协议 §6/§7 审核。
> 结构：§0 环境（非原生）→ §1 第一路 → §2 第二路 → §3 第三路 → §4 耗时与内存 → §5 偏离候选。
> 更新规则：冒烟后只许修脚本错误；改计划映射须记入协议 §7，并在 §5 追加。
> 引用：上游 协议、track1.md/track2.md 的 #34020 节、`specs_codex/idaholab_moose_34020/spec.md` §2/4/5、合并提交 f90f508 源码；下游 `dq1/out/mo34020_t{1,2,3}.json`。

## §0 环境：非原生重算（偏离候选）

- compute node A、compute node B 无 MOOSE/libMesh/PETSc 构建（只读检查，同 #33307）。不构建。
- `mortar_port.py` 按合并提交源码逐行重算 PR 改动的全部数值逻辑：
  - `ComputeDynamicWeightedGapLMMechanicalContact::computeQpProperties/computeQpIProperties`：`_qp_factor = JxW_msm·coord`，`weighted_gap[dof].second += test·_qp_factor`（合并后无条件累加）；
  - `ComputeDynamicFrictionalForceLMMechanicalContact::computeQpProperties/computeQpIProperties`：`qp_tv = −relative_velocity·JxW·coord`，`relative_velocity = 主面 − 从面`（两侧各按本侧单元形函数插值节点速度），`weighted_tv[dof][α] += test·qp_tv·t_α(node)`；
  - `timestepSetup`（本进程拥有的 dof：s = weighted_tv / weighted_gap.second）、`oldNodalSlipRate`（无历史 → 0）、`computeFrictionValue`（F(0, p_old, √(s0²+s1²+1e-24), 0) 或常数 mu）；
  - 构造函数参数校验（`mu > 0` 范围检查、两者都给/都不给报错）。
- 几何限定：平直界面；2D 为线段界面（1 个切向），3D 为轴对齐矩形面（2 个切向）；主从两侧张量积网格，mortar 段＝两侧断点并集划出的区间/矩形，每段 Gauss 积分。dual 基按 libMesh `compute_dual_shape_coeffs`（dual_coeff = A⁻¹D，A 单元质量阵、D=diag∫φ，ψ_i = Σ_j dual_coeff(j,i)φ_j；矩形单元为一维系数的张量积）。
- 不能重算：非线性接触求解（法向/切向 LM、|λ_t|/λ_n、位移）、Newmark 速度更新、MPI 通信本身（第三路以积分点随机分区求部分和代替）、edge dropping 的节点取舍、restep、elemental LM。
- 凡计划用 μ_eff = |λ_t|/λ_n 反推滑移速率处，改读传给摩擦函数的实参与返回值（滑移状态下 Coulomb 律取等号时 μ_eff 等于该返回值）。
- 速度场直接给定（不经 Newmark）；滑移场记为「从面速度 − 主面速度」＝代码中的 −v_rel，主面静止。
- 种子 20261017＋10×10＋路号。

## §1 第一路（track1.md，2D 平面问题）

分支：`use_dual` ∈ {false, true}（k=2）。

| 计划原文 | 代码 |
|---|---|
| 工况域：\|V\|∈[1e-4,10] m/s；p∈[1e3,1e9] Pa；界面 4–64 个单元；Newmark β=0.25、γ=0.5；μ=0.2+0.1·tanh(y/0.1) | `t1_sample`：\|V\| 线性均匀、符号随机；p 线性均匀；单元数整数均匀 4–64；Δt、β、γ 不进入重算 |
| A：常向量相对切向速度 V，每个活动节点 s_j = −V·t_j，\|s_j\|=\|V\|，与网格尺寸、单元阶次、积分阶次、dual 与否无关 | `t1_A_s_j==-V_*`（相对速度＝主面−从面＝V） |
| A 核对：μ_eff 在第 ≥2 步与 g(\|V\|) 比；g(y)=y 时直接读出 \|s_j\|；相对 1e-6；网格加密 2 次 μ_eff 不变 | `t1_A_mu_eff==g(|V|)_*`、`t1_A_identity_law_|s_j|==|V|_*`；每点对 加密 ×1/×2/×4、单元阶次 1/2、积分点数 p+1/p+2 共 12 种组合都做；`t1_A_mu_eff_mesh/order/quadrature_invariance` |
| B：同网格 `function_friction=g` 与 `mu=g(V)` 比，相对 1e-6 | `t1_B_function_vs_constant_mu_*`（只比摩擦系数，不比位移/LM） |
| C：上块匀速转动，切向速度沿界面线性分布，use_dual=true，μ_eff 推出的 \|s_j\| 与解析线性场比 | dual 支：`t1_C_dual_linear_field_*` |
| 第 1 步与刚进入界面的节点取 s=0，核对时排除 | 只看由上一步运动学建立的历史 |

实施判断：
1. 界面长度未给：取 1（从面 [0,1]）；主面取 [−0.5,1.5]、单元数与从面相同（单元尺寸 2 倍，非匹配），仿夹具上下块宽度不同。
2. 线性场：v(x) = V_lo + (V−V_lo)·x，V_lo = 1e-4·sign(V)，使全界面速度落在工况域 [1e-4, \|V\|]（若取 V·x，x=0 处速度为 0，√(·+1e-24) 的正则项会让 abs 1e-12 的默认容差恰在边界）。C 无容差，用默认 rel 1e-6、abs 1e-12。
3. 不变量 A 写了「单元阶次、积分阶次」无关，核对方法只写网格加密；按协议 §2 纳入：阶次 {1,2}、每方向积分点 {p+1, p+2}。
4. 2D 中「绕法向转动」不在平面内，按计划括注改为沿界面线性分布的切向速度。

## §2 第二路（track2.md）

分支（k=8）：{2D, 3D} × {standard, dual} × {1 阶, 2 阶}。

| 计划原文 | 代码 |
|---|---|
| 工况域：滑移速度 1e-3–10 m/s；μ₀ 0.1–0.8；h 三档加密；一/二阶单元；2D/3D；standard/dual | `t2_sample`：速度线性均匀，2D 符号随机，3D 方向角 U[0,2π)；μ₀ 线性均匀 |
| 模板：`frictional-mortar-3d-dynamics-light-function`（上块 0.5 宽、下块 1 宽，各 2×2 单元） | 从面 [−0.25,0.25]^d、主面 [−0.5,0.5]^d，基础 2 单元/方向，三档加密 = 2、4、8 单元/方向 |
| I1：恒定切向速度 v，`function_friction`=0.5−0.1·v，\|λ_t\|/λ_n 与 0.5−0.1\|v\| 比 rtol 1e-6；3 档网格×2 种积分阶不随之变化 | `t2_I1_s_i==v.t_i_*`、`t2_I1_mu==0.5-0.1|v|_*`（每点 3 档 × 积分点 {p+1,p+2}） |
| I2：use_dual=true，上块绕法向匀速转动，节点 μ 反推的 v 与 ω r_i 比，相对 1e-8 | dual 支：`t2_I2_dual_rotation_s_i==omega_x_r_*`（分量）与 `t2_I2_dual_rotation_mu_inferred_speed_*`（由 μ=0.5−0.1y 反推 y） |
| I3：μ(v)=μ₀−b·v，μ₀−b·v_max>0 时摩擦系数在 [μ₀−b·v_max, μ₀]、不为负；首步/新进入节点 μ=μ(p_old,0) | `t2_I3_mu_in_[mu0-b*vmax,mu0]_and_>0_*`（条件不满足的点不检）、`t2_I3_first_step_mu==mu(p_old,0)` |

实施判断：
1. I3 的 b 未给：取 I1 律的 0.1；v_max 取该点抽到的速度。
2. 旧压力未在工况域给出、且两条摩擦律都不依赖压力：固定 p_old=1。
3. 转动：3D 为 ω ẑ×r（中心在界面中点），ω = v/(L/√2) 使角点速度＝抽到的 v；2D 用线性场 ω·x，ω = v/(L/2)。容差相对 1e-8，计划未给绝对容差，按协议用 1e-12。中心节点理论滑移为 0，由 μ 反推的速度至少是 √1e-24 = 1e-12，与 abs 1e-12 恰在边界（可能出拒绝，属正则约定，交裁定）。
4. 滑移方向约定：计划写的 v 取为 PR 描述里的「−v_rel」（从面减主面），故期望 s_i = v·t_i。
5. 二阶 dual 基：EDGE3/QUAD9 的 ∫φ 都为正，libMesh 不走 QUAD8/TRI6 修正分支；按 A⁻¹D 计算。

## §3 第三路（spec.md §2/§4/§5）

分支（k=16）：{2D,3D} × normalize_c{false,true} × friction_projection_degree{ONE,TWO} × correct_edge_dropping{false,true}（计划 §5「核心矩阵」）。

| 计划原文 | 代码 |
|---|---|
| 平面 patch 长 0.1–10，切向每方向 1–8 单元，匹配与非匹配；速度分量有限带符号（列表 0…1 两种符号）；压力 0…1000；dual 乘子；积分阶「足够」 | `t3_sample`：L 线性均匀；从/主面单元数各自整数均匀 1–8（相等即匹配）；分量 U[−1,1]；p_old U[0,1000]；积分点 {2,3}；一阶单元 |
| I1：W_i s_i = V_i | `t3_I1_W*s==V` |
| I1：均匀相对速度 → s=b，与支撑长度×0.1/1/10、加密 1/2/4/8、非匹配主面细分、积分阶无关；容差 abs 1e-12 v*＋rel 1e-9 | `t3_I1_constant_reproduction_s==b`、`t3_I1_support_scale_*`、`t3_I1_refine_x*` |
| I1：仿射非均匀滑移场，用 dual 检验函数与 mortar 段测度独立积分 V_i/W_i，比投影后的范数；含支撑内反向滑移 | `t3_I1_affine_…`、`t3_I1_opposing_…`：独立参考 `oracle_affine`（只在从面单元上高阶 Gauss 直接积分解析场，不经主面插值与分段），比 √(\|s\|²+1e-24) |
| I2：≥3 个已接受步、加速滑移与变化法向载荷；第 n+1 步回调实参 = 第 n 步投影；首步速度 1e-12、旧压力、t 与 z 为 0；新进入的 dof 取零滑移 | `t3_I2_first_step_args_…`、`t3_I2_step{2,3}_speed_from_step…`、`t3_I2_step*_old_pressure`（容差 1e-9 p*＋1e-7\|p\|）、`t3_I2_entering_dof_zero_slip_convention`（上一步主面只覆盖前一半单元） |
| I2：`--test-restep` | 不能重算，`not_covered` |
| I3：分区求和与串行一致 1e-10 v*＋1e-8\|s\|；normalize_c、projection degree、edge dropping 不改变本构实参 | `t3_I3_partition_{2,4}_vs_serial`（积分点随机分给 2/4 个「进程」）；三个开关在重算中不进入 s，记一条恒真检查并在 not_covered 说明 |
| I4：2D ±b；3D (b,0),(0,b),(−b,0),(0,−b),(3b/5,4b/5) 及反向；几何与场一起刚性转动；速度相等 1e-10 相对＋I1 绝对 | `t3_I4_speed_isotropic_over_directions`、`t3_I4_rigid_rotation_speed_unchanged`（3D，切向基随转动） |
| I5：常数 mu=0.3 与恒为 0.3 的 ParsedFunction 结果一致；mu=0、−0.1、两者都不给、两者都给须报错，只给正常数成功 | `t3_I5_constant_mu_vs_constant_function`；必测点 `valid_*` |
| I5：投影恒等式、圆盘界、耗散、epsilon 门限、merge/base 回归 | 需 LM 解，`not_covered` |
| 速度分量列表 0,1e-12,1e-8,1e-4,0.01,0.02,0.1,1（两种符号），压力 0, ε/2, ε, 2ε, 1, 10, 100, 1000（ε=1e-7） | 必测点 `2D_b=…`、`3D_b=(…)`（96 例，每例 8 个压力）：常数再现、实参 (0,p,速度,0) 精确、所给律在速度 ≤0.02 时为正 |
| wear_depth 零与非零 | 重算中 W 只含 `_qp_factor`，与 wear 无关；未单列（恒真） |

实施判断：
1. 速度分量取值只给离散列表（0…1、两种符号）：抽样区间取 [−1, 1]，列表值作必测点。
2. 主面与从面同一外框（全支撑），只改主面细分；部分覆盖只在 I2「新进入」检查中用，覆盖边界对齐从面单元边界（保证每个节点 W_i>0）。
3. 仿射场梯度 G 的各分量 U[−1,1]/L；反向滑移场 b·(x−L/2)/(L/2)。
4. 「Use prescribed kinematics」：重算本就是给定运动学，不涉及非线性轨迹。
5. 一阶单元（夹具 HEX8）；dual 乘子（计划称已建立的配置）。

## §4 预计耗时与峰值内存（推测，未运行）

| 路 | 单点 | 10⁴ 点 | 依据 |
|---|---|---|---|
| t1 | ~20–50 ms（12 种组合，最多 256 单元×2 阶） | ~5–8 min | numpy 向量化累加＋逐节点 Python 摩擦调用 |
| t2 | ~10–30 ms（3D 二阶 8×8 单元最大） | ~3–5 min | |
| t3 | ~50–150 ms（约 23 次累加，3D 加密 ×8 时约 4.7 万积分点） | ~10–25 min | `np.add.at` |

峰值内存 <300 MB；JobSpec 申报 1 GiB。

## §5 偏离候选

- **非原生**：节点无 MOOSE 构建；被改逻辑在 Python 中重算（覆盖清单见 §0）。μ_eff=|λ_t|/λ_n 一律以摩擦函数实参/返回值代替；需要非线性求解、restep、真 MPI、edge dropping 的检查未测。
