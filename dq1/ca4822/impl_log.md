# ca4822 实施日志（ClimaAtmos#4822，ΔQ1 统一抽查）

> 目的：记录三路检查计划逐句对应到 `spotcheck.jl` 哪段代码，以及全部实施判断。
> 结构：§0 共用 → §1–§3 三路 → §4 未实现项（偏离候选）→ §5 耗时与内存 → §6 冒烟后修改。
> 更新规则：冒烟后只修脚本错误，记入 §6。引用：协议 `notes/20261017-prereg-dQ1-uniform-spotcheck.md`；共用库 `../jl_common/spotcheck_common.jl`。

## §0 共用

- 环境：compute node B `env/ca`，ClimaAtmos ba2bf122（即 #4822 合并提交；与原生重放同一环境）。有原生环境。单线程。
- 入口：`--tracks 1,2,3 --n N --outdir DIR [--budgets …]`，一进程三路（由 `--track/--out` 改）。种子 `20261017 + 3×10 + 路号`。
- 模型级配置（实施判断）：微物理 `1M`、`cloud_model=quadrature`、`use_sgs_quadrature=true`、`tq_correlation_model=diagnosed`（第三路主路径；第一、二路未指定微物理，统一用它，1M 下 T 由内能直接反解、无饱和调整，便于制造状态）。
  - plane：x_max=2000 km、8 谱元（Δx_h = 250/3 ≈ 83 km）、z_max=2 km、10 层均匀（Δz=200 m）。
  - box（第三路）：200 km × 200 km、4×4 谱元（Δx_h ≈ 16.7 km）、z_max=2 km、20 层（Δz=100 m）。
  - 精度 Float64、Float32 各一套（共 4 个配置，plane 另有 4/16 谱元两套供第二路分辨率收敛）。
- 制造状态（`manufacture!`）：给定 θ_li(x,y,z) 与 q_tot，ρ 沿用初值，凝结物与风为 0；由 θ_li 反解 T（不动点迭代，调用 `TD.liquid_ice_pottemp`），ρe_tot = ρ(e_int(T,q)+Φ)，再 `CA.set_precomputed_quantities!`。
- 水平场用「锯齿」：每个谱元内斜率 ±1、相邻谱元反号、在谱元边界连续（偶数谱元周期连续），故每个谱元内梯度为常数，谱元导数、`lumpedₕ` 限制与加权 DSS 对 |∇ψ|²、∇ψ₁·∇ψ₂ 都精确——解析梯度可作精确参照，也避开周期接缝的仿射跳变（第三路计划明写）。
- 参数开关：`with_sgs` 复制 `p.params.sgs_quadrature_params` 并替换字段（c_Δx、c_Δz、C_cov=`diagnostic_covariance_coeff`、ρ_turb=`Tq_correlation_coefficient`），`with_corr` 换 `AtmosWater.tq_correlation_model`；随后只调 `CA.set_covariance_cache!(Y, p2, thermo_params)` 重算协方差（precomputed 共用）。这样不必为每个参数值重建模拟。
- 层：排除首末层（竖直梯度在边界层为单侧/边界值），第三路排除首末两层（计划写 interior levels）。

## §1 第一路（track1.md）

| 计划原文 | 代码 |
|---|---|
| 逐点：直接调用 `CA.isentropic_slope_ratio`、`CA.isentropic_cell_cap`、`CA.isentropic_gradient_invariant` | `r_of`、`f_of`、`G_of` |
| A：随机抽 ∇_hθ、θ_z ∈ [θz_min,1e-2]、b，令 ∇_hq=b∇_hθ、q_z=bθ_z，检 I ≤ 10ε(inv_q + r²inv_θ) | 分支 `Ai_affine_{F64,F32}` |
| 再对任意 ∇q 检 I(∇q+b∇θ, q_z+bθ_z) 与 I(∇q,q_z) 之差 ≤ 10ε·(各项量级) | 分支 `Aii_gauge_*` |
| B：对 s 网格算 f、ℓ²=fΔx²，检 (sℓ)²+fΔz²=Δz²（1e-6/F32、1e-14/F64）与投影上界 | 分支 `B_cellcap_*` |
| 整场：plane 配置、diagnosed，状态 q_tot = a + b·θ_li，调 `CA.set_covariance_cache!`，比较 c_Δx 取 0 与非 0 的 q′q′，差 ≤ 100ε·q′q′ | 分支 `field_affine_cdx_*` |
| C：一般状态读 ᶜT′q′、T′T′、q′q′（clamp 前），检 (T′q′)² ≤ T′T′·q′q′·(1+100ε) | 分支 `field_CS_*` |
| 工况域 Δx 1–250 km、Δz 20–1000 m、θ_z −5e-3–1e-2、|∇_hθ| 0–1e-4、q_z −1e-5–0、|∇_hq| 0–1e-8；F32 与 F64 | 抽样 |

实施判断：
1. A 的 b 未给区间：由 q_z ∈ [−1e-5,0] 抽出再令 b = q_z/θ_z；A(ii) 的「任意 b」同法另抽一个 q_z。∇ 为二维向量，方向在 [0,2π) 均匀。
2. A(ii) 的「各项量级」取两次求值中 inv_q + 2|r·inv_θq| + r²inv_θ 的较大者。
3. B 的「投影上界」容差与恒等式同（相对 1e-6/1e-14）。
4. 整场：θ_z 取 [θz_min,1e-2]（A 的适用域；θ_z < θz_min 时计划自己写明 A 不成立）。为使状态自洽（q 非负、|∇_hq| 在域内），拒绝抽样 |b|·|∇_hθ| ≤ 1e-8；q 在地面参考 θ=300 K 处为 0.022。c_Δx「非 0」取库默认 1。比较量取所有内部节点 |Δq′q′|/q′q′ 的最大值。其余参数为库默认（湍流块开，但风为零时 Smagorinsky 混合长为 0）。
5. C：一般状态 θ、q 各自独立（q 的锯齿相位错开半个谱元、符号随机），取 max (T′q′)²/(T′T′q′q′) 与 1+100ε 比。读的是 precomputed 的 ᶜT′q′（相关系数 clamp 之前的量；`σ_q` 上限对 T′q′ 与 q′q′ 同比缩放，不改变这一比值）。

## §2 第二路（track2.md）

| 计划原文 | 代码 |
|---|---|
| 一层：用导出函数对域随机抽样（F64/F32），拼水平块 c_g c² Δx² cap·inv 与竖直块 c_g c² Δz² cap·q_z²，检和 = c_g c² Δz² q_z²，rtol 1e-12/1e-5 | 分支 `pt_I1_*`；另测 θ_z ≤ θz_min：`pt_I1_floor_*` |
| I2：随机 F′，inv_q=F′²|∇θ|²、inv_θq=F′|∇θ|²、q_z=F′θ_z，inv ≤ 4ε·inv_q | 分支 `pt_I2_*`（附 0<f≤1、(sℓ)²+Δz²f=Δz²） |
| 二层：plane、c_Δx=c_Δz、湍流块关；q_tot 只随 z 变、θ 面倾斜；调 `set_covariance_cache!` 读 ᶜq′q′ 与 c_g c² Δz² (∂_z q)² 比，相对 1e-3 并检随分辨率收敛；避开 σ_q 截断 | 分支 `model_I1_*`、`model_I1_floor_*` |
| I3：同上构造 q = a + bθ_li，读 T′T′/q′q′（θ 基对比），验证 q′q′/θ′θ′ = f·b²，并报告与 b² 的偏差 | 分支 `model_I3_*`（`report_I3_ratio_over_b2` 只记录不判定） |
| 工况域 Δx_h 25–250 km、Δz 50–1000 m、θ_z 1e-3–1e-2（另测 ≤θz_min）、|∇_hθ| 0–1e-4、q_z −5e-6–0、|∇_hq| 0–1e-7 | 抽样 |

实施判断：
1. c 未给值：取库默认 c_Δx=c_Δz=1；c_g 取默认 1/12。
2. floor 分支 θ_z 区间计划未给：取 [−1e-3, 1e-3]（量级与 θz_min 相同）。
3. I2 的 F′ 由抽到的 q_z/θ_z 给出（与工况域一致）。f 与恒等式的附加检查来自不变量 I2 原文，方法段未给容差，用相对 1e-6。
4. 二层：湍流块关用 `diagnostic_covariance_coeff = 0`；分辨率用 plane 4/8/16 个谱元（Δx_h ≈ 167/83/42 km，都在域内）；「收敛」判据为 16 谱元的误差不大于 4 谱元的（或已 ≤1e-12）。F32 只做 8 谱元一档。竖直梯度参照用解析 q_z（线性剖面在内部层精确）。σ_q 截断：预期值 ≥ 0.9·(0.5q)² 的节点排除。
5. I3：θ′θ′ = T′T′/(∂T/∂θ)²，雅可比用库内 `CA.compute_∂T_∂θ!`；f 用解析 |∇θ|、库内 Δx_h 与 Δz 算；容差沿用二层的 1e-3；取全体内部节点中相对偏差最大者。为保持 q 自洽，拒绝抽样 |b||∇θ| ≤ 1e-7。

## §3 第三路（spec.md §5，box 网格）

| 计划原文（摘） | 代码 |
|---|---|
| 公开构造（AtmosConfig/get_simulation），读注册诊断 `env_temperature_variance`、`env_q_tot_variance`、`env_q_tot_temperature_covariance`、`env_q_tot_temperature_correlation` 的 `.compute(Y,p,t)` | `diag_vals` |
| 内部谱元、内部层用局部仿射场，不在周期接缝处造仿射跳变；参照用解析梯度；τ=1e-10/2e-5；G 的尺度 S_G = C+2|rB|+r²A；PSD 用 C² ≤ V_TV_q + τ max(V_TV_q, C²) | 锯齿场、`interior_mask`、`τ_of` |
| I1：C_cov=0、c_Δz=0、c_Δx=1；a>d；b 两号、零梯度、q=αθ+β、近相消、0/π/2/π 与随机角；V_q = h f G、V_T = 0、诊断交叉协方差 = 0；旋转 90° 与轴反转；floor 以下用正则化式 | 分支 `I1_isentropic_contrast`；角度 0、π/2、π 为必测点 |
| I2：c_Δz=0、∇q≠0、b=0，由 V_q/(hC) 反推 f；x = 0,1e-3,0.1,1,10,100，a = −d,0,d(1±1e-3),2d,10d；检 f 与偏移不等式 | 分支 `I2_slope_cap`（工况域抽样）＋ 72 个必测点（x×a×精度） |
| I3：两种相关模型、ρ_turb=−1,−0.5,0,0.5,1、块开关、湍流有无；无截断/截断/q=0；与完整变换＋截断参照比较，再查 PSD 与相关界；跨 V_TV_q = eps² | 分支 `I3_realizable`；ρ_turb 五值×两模型、eps² 两侧为必测点 |
| I4：c_Δx、c_Δz 的四种 0/非 0 组合，平与陡 θ 面，c = 0.5,1,2；C_cov=0、∇_hq=0、两尺度=1 时 V_q = c_g b²Z²；c_Δx=0 全竖直块；两者为 0 时方差为零与回退相关；V_T 与 J²v a² 比 | 分支 `I4_partition` |
| I5：A>0、b≠0、c_Δz=1、C_cov=0；c_Δx = 1,1e-1,…,1e-6,0；比较极限与 0 处的 q 方差与协方差 | 分支 `I5_cdx_continuity`（`report_I5_relative_jump_qvar` 只记录） |

实施判断：
1. 参照式按 §2 的不变量与 §5 的描述独立组装（`reference_cov`）：V_θ = c_g(c_ΔzZ)²a²；V_q = h f G + c_g(c_ΔzZ)² f b²；C_θq = c_g(c_ΔzZ)² f a b；V_T = J²V_θ、C_Tq = J·C_θq；σ_q 上限 r_max=0.5 与 T′q′ 同比缩放。c_Δx=0 时 f 不作用。J 取库内 `compute_∂T_∂θ!`（计划写「finite J from the library's thermodynamics」）。
2. 「零梯度」等零输入：谱元导数对水平常值场只到舍入，按计划「otherwise the same scaled absolute tolerance」给一个舍入级绝对下限（100·h·(ε·θ/Δx_h)² 量级）。
3. 状态自洽：拒绝抽样使 θ ∈ (250, 320) K、q ∈ [0, 0.025]；|∇_hq| 上界按此收紧到 1e-7（计划域到 1e-6 时 q 会为负）。I2 的必测点若不自洽记为空点（x=100 配 a=10d 等）。
4. 旋转 90° 与轴反转：以 (α+π/2, β+π/2) 与 (−α, −β) 重造状态，比较 V_q（锯齿基下这正是对应的旋转/反射）。
5. I3 非零湍流：风为零时 Smagorinsky 混合长为 0，故加线性切变风 u = S·z（S ∈ U[0,1e-2] s⁻¹），能量中补动能后重算 precomputed；非零湍流情形无独立参照（需混合长），只检 PSD、界与回退。C_cov 非零取 U[0, 4.2]（默认 2.1 的两倍内）。q=0 情形约 10%。「刚好截断／强截断」未专门构造，由抽样自然覆盖并在观测里可见。
6. eps² 两侧：只留竖直块，取 a=b 使 V_TV_q 约为 eps(FT)² 的 0.5 与 2 倍（按 J≈1 粗估）。
7. I5 的判据「lim_{c_Δx→0⁺} Σ = Σ(0)」照写为检查（容差 τ），计划自己说预期跳变为 (1−f)vb²，属于单列报告的条件性判据。

可能的系统性拒绝（推测，供解读）：第三路 F32 的参照比较（`ref_VT` 等）中，θ_li 由 T 在 F32 下重算，单点舍入约 3e-5 K，在 Δz=100 m 上给竖直梯度约 3e-7 K/m 的绝对误差，相对于 |a| 只有在 |a| > 0.04 K/m 时才低于 τ₃₂=2e-5。若 F32 参照检查大面积拒绝，可能属于「计划容差在 F32 下不可达」，交裁定。

## §4 未实现项（偏离候选）

- 第三路：单柱与球面网格、平坦球面与地形跟随网格（只做了 box）；0M、2M、2MP3 微物理；`turbconv` 的 `edonly_edmfx`、`prognostic_edmfx`；`sgs_distribution` = gaussian/lognormal/mean；`uses_covariances` 的各种替代配置（含 MLCloud 网络）；d=0、d<0 的配置构造探测；一步 `CA.solve_atmos!` 冒烟；不同 c_Δx/c_Δz 下 I2 的「非零 b 推断下限后比值」；凝结物（液/冰/混合相/降水）状态——制造状态只用了无凝结物的 1M 状态。
- 第一路：sphere 配置（只做 plane）。
- 理由：每种配置都要重新编译 ClimaAtmos（约 1–2 分钟、数 GiB 内存），全部覆盖会超出 compute node B 24 GiB 与单 job 4 小时；是否补做请主会话裁定。

## §5 耗时与内存估计（推测）

- 依据：同节点原生重放 job 9cd55d6e 中 ClimaAtmos 部分（plane 与 EDMF 单柱两三个配置、各推进若干步）约 2.5 分钟、RSS 峰值约 11 GiB。
- 第一路：逐点检查约微秒级；整场每点约 0.1–0.3 s（制造状态＋两次协方差）；10 个分支中 4 个整场 → 10⁴ 点约 15–25 分钟（含编译），预算 2400 s。
- 第二路：model_I1 F64 每点三档分辨率约 0.5 s，其余约 0.2 s → 10⁴ 点约 40–60 分钟，预算 3000 s。
- 第三路：box 每点 1–3 次制造状态＋最多 8 次协方差 → 约 0.5–1.5 s → 10⁴ 点约 1.5–4 小时，预算 3000 s。
- 编译：4 个配置（plane/box × F64/F32）约 5–10 分钟。
- 峰值内存估 12–16 GiB（4 个较轻的非 EDMF 配置），申报 20 GiB。全量 job 估 9000 s，硬上限 14400 s。

## §6 冒烟后修改记录

（空）
