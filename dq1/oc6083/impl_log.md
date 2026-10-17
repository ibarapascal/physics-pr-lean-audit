# oc6083 实施日志（Oceananigans#6083，ΔQ1 统一抽查）

> 目的：记录三路检查计划逐句对应到 `spotcheck.jl` 哪段代码，以及全部实施判断。
> 结构：§0 共用 → §1–§3 三路 → §4 未实现项（偏离候选）→ §5 耗时与内存 → §6 冒烟后修改。
> 更新规则：冒烟后只修脚本错误，记入 §6。引用：协议 `notes/20261017-prereg-dQ1-uniform-spotcheck.md`；共用库 `../jl_common/spotcheck_common.jl`（参数、种子、轮转、预算、判定、重放、JSON 同 oc6024 §0）。

## §0 共用

- 环境：compute node B `env/ocean`，Oceananigans 24de7026（**晚于本 PR 合并提交 1734ff9e**；与原生重放同一环境，按协议 §1 用它）。有原生环境。单线程。
- 入口：`--tracks 1,2,3 --n N --outdir DIR [--budgets …]`，一个进程三路（由 `--track/--out` 改，省编译）。
- 种子 `20261017 + 5×10 + 路号`。
- 漂移速度扩展：照 PR 测试文件 `test/unit/adaptive_implicit_vertical_advection.jl` 的 `SinkingParticles`（`required_biogeochemical_tracers`、`biogeochemical_drift_velocity`、`Adapt.adapt_structure`）；另写 `TwoTracerBGC`（D 有漂移、E 无）供第三路查示踪物查找错位。
- 单柱：`RectilinearGrid(size=40, z=(-400,0), topology=(Flat,Flat,Bounded))`，`momentum_advection=nothing, free_surface=nothing`（同 PR 测试）；w 场内部面为 W、顶底面为 0（`wfield`）。

## §1 第一路（track1.md）

| 计划原文 | 代码 |
|---|---|
| A：WENO(order=5, AIVA(cfl=0.5))；(i) SinkingParticles 漂移 (ii) `PrescribedVelocityFields(w=…)` (iii) `forcing=(D=AdvectiveForcing(w=…),)`；QAB2 与 SplitRK3 各推进若干步，比较 `interior(D)`，相对 1e-12；(ii) 若走不同分支至少 1e-10 | 分支 `A_three_paths_{QAB2,SRK3}`：`A_forcing_vs_bgc_rel1e-12`、`A_prescribed_vs_bgc_rel1e-10` |
| B：同上设置，Gaussian 中心 −100 m、宽 20 m，C ∈ {0.4,4,10}；Σc 相对 1e-12；质心位移 vs w_s nΔt（建议 <1e-2）；cfl 设极小几乎全隐式，误差应降到 1e-6 量级 | 分支 `B_cfl0.5_*`、`B_nearly_implicit_*`；C 三值为必测点（12 个） |
| C：16×16×16 三维、倾斜等密面，`IsopycnalSkewSymmetricDiffusivity(κ_skew, slope_limiter=FluxTapering(S_max), skew_flux_formulation=AdvectiveFormulation())`；读 uₑ,vₑ,wₑ，离散散度机器精度为零（相对 1e-12）；局部陡坡检 |ψ| ≤ κS_max（直接调用 `ϵSxᶠᶜᶠ`）；b_z ≤ 0 处 ψ=0 | 分支 `C_GM_3d`（`divergence_stats`、`taper_stats`） |
| 工况域 w_s 1–200 m/day、C 0.1–10、GM κ_skew 100–2000、S_max 1e-3–1e-2 | 各分支均匀抽样 |

实施判断：
1. 「推进若干步」与 B 的步数：取 n = max(1, ⌊120 m/(CΔz)⌋)，即位移约 120 m（同 PR 测试），保持示踪物远离边界（B 的不变量前提）；A 用同一步数。
2. Gaussian「宽 20 m」取 exp(−(z+100)²/(2·20²))，与 PR 测试同式。A 的初值计划未写，用同一 Gaussian。
3. B 用漂移路径（PR 改动的路径）。「极小 cfl」取 1e-6；「1e-6 量级」取容差 1e-5（同一数量级）。
4. (ii) 的两档容差：超过 1e-10 才算拒绝，实测相对差记在观测里（若介于 1e-12 与 1e-10 之间可在 JSON 的 obs 范围里看到）。
5. 相对差定义 max|D_a−D_b| / max|D_b|。
6. C 的网格与场：Bounded³，L_x=L_y=100 km，H=1000 m（计划只给格点数）；b = N²z − 2N²·max(z−z_u,0) + Δb·tanh((x−L_x/2)/L_f) + y 向同类锋；N²=1e-5、L_f=10 km、z_u=−100 m（上方 b_z<0 为不稳定层）；Δb 使最大坡度 = steep·S_max，steep ∈ U[10,100]（计划只说 |S| ≫ S_max）。散度相对 max_cell Σ|A U|。ψ=κϵS，检 κ·max|ϵSx| ≤ κS_max（相对 1e-12 余量）；不稳定层判据：面高度 > z_u+Δz（两侧格心都在不稳定层）处 ϵS 须恰为 0。`minimum_bz` 用默认 0。

## §2 第二路（track2.md）

| 计划原文 | 代码 |
|---|---|
| I1：模型 A `PrescribedVelocityFields(w)`、模型 B `SinkingParticles`（输运速度 0），同一 Gaussian，推进 3–20 步，比较 interior(D)，rtol 1e-10；扫 C ∈ {0.25,0.5,1,4,10} | 分支 `I1I2_{QAB2,SRK3}_{WENO5,Centered2,Upwind3}`：`I1_prescribed_vs_drift_rel1e-10`；C 五值 × 6 配置为必测点（30 个） |
| I2：sum(D) 相对漂移 ≤1e-12；质心位移相对误差 ≤5% | `I2_inventory_*`、`I2_centroid_drift_rel5pct` |
| I3：16×16×16 倾斜 b 场＋随机扰动＋局部不稳定层，`compute_closure_fields!` 后读 uₑ,vₑ,wₑ，用 Operators 的 δ 与面积算子算散度，相对 max|A uₑ| ≤1e-12；核内调用 `ϵSxᶠᶜᶠ` 抽查坡度上限；等密面水平或 b_z ≤ 0 时 uₑ=vₑ=wₑ=0 | 分支 `I3_GM_tilted`（散度、|ϵSx|,|ϵSy| ≤ S_max、不稳定层为 0）、`I3_GM_flat`（全零） |
| 工况域：漂移 1–200 m/day，Δz 1–20，Δt 10 min–1 day，Courant 0.1–10；WENO5 与 Centered/UpwindBiased；κ_skew、S_max；b_z 含负值与零 | 抽样见 `t2_sampler`、`gmsamp` |

实施判断：
1. 网格按计划方法固定为 40 层、400 m（Δz=10 m，在 1–20 内），不另抽 Δz。抽 w_s 与 C，Δt = CΔz/w_s，拒绝抽样到 Δt ∈ [600 s, 86 400 s]。步数在 3–20 均匀整数。
2. Centered 取 order=2、UpwindBiased 取 order=3（库默认阶）；三者都套 AIVA(cfl=0.5)。
3. 质心检查只在预计位移 ≤ 200 m（终点质心不低于 −300 m）时做，否则该子检查不出现（「远离边界」是不变量前提）。
4. GM 扰动：幅度 1e-3·N²H 的均匀噪声（每点独立种子，种子记入输入）；倾斜 M²=slope·N²，slope ∈ U[1e-4,1e-1]，y 向分量比例 U[−1,1]；不稳定层顶 z_u 在柱深 50%–95% 处。N²=1e-5 固定（计划未给）。散度按计划用 `δxᶜᵃᵃ/δyᵃᶜᵃ/δzᵃᵃᶜ` 与 `Ax_qᶠᶜᶜ/Ay_qᶜᶠᶜ/Az_qᶜᶜᶠ`。
5. 「等密面水平」分支要求逐位为 0。

## §3 第三路（spec.md §5，CPU 部分）

| 计划原文（摘） | 代码 |
|---|---|
| PR 原设置：三步 0.4 天、−100 m/day，库存误差 < T，下沉 120 m 在 5% 内；镜像上升；两种时间推进、两种精度 | 必测点 `pr_setup`（8 个） |
| 两种时间推进、两种精度、零速、级 CFL 比 C/C* ∈ {0,0.1,0.99,1,1.01,2,8}；C* ∈ {0.25,0.5,1}；WENO5、UpwindBiased(1)、Centered(2) | 分支 `I1I2_source_matrix` 抽样维度 |
| 漂移单独、AdvectiveForcing 单独、两个平流强迫加普通零 Forcing、相消的反向漂移、阈值两侧、两示踪物不同分配；整数组比较，容差 T·max(1,‖q0‖∞) | 子模式 `bgc_vs_forcing`、`two_forcings`、`cancel`、`two_tracers` |
| 阈值以下 adaptive 与 `ExplicitTimeDiscretization()` 一致 | 子模式 `explicit_equiv`（只在 C/C* < 1 时） |
| 每次运行检 I2 库存 | `I2_inventory_*` |
| GM 闭合源 vs 用 AdvectiveForcing 提供观测到的涡致速度，Δt 减半，差/Δt 收敛到 0 | 分支 `I1_GM_source_vs_forcing` |
| 闭合 nothing、κ_skew=nothing、κ=0、单个 GM、长度 1–5 的元组中 GM 置于各位置（其余零系数） | 分支 `I1_closure_tuples` |
| I3：16×16×32 有界网格，线性 b=N²(z−S_xx−S_yy)，`SmallSlopeIsopycnalTensor(minimum_bz=m)`、`FluxTapering(S_m)`、AdvectiveFormulation、κ 为 z 的正线性函数；内部 uₑ=−κ′ϵSx 等与独立分段式比，相对 T；扫坡度范数比、象限、下限情形，S_x=S_y=0.8S_m；x/y 互换与符号反射；反推坡度范数界；Float64 参数用于 Float32 网格 | 分支 `I3_taper`（`taper_case`） |
| I4：周期水平变化、稳定层结、κ 可变；逐格带号面通量和 ≤ T·Σ|A U|；顶底无法向通量；x–z、y–z、三维；均匀与拉伸 z；常数被动示踪物在定常预设涡致速度下保持常数 | 分支 `I4_divergence` |
| I5（固定网格部分）：物理水平 b=N²z+b₀ 不产生 GM 输运 | 分支 `I5_fixed_grid_horizontal_b` |
| 容差 T₆₄=1e-11、T₃₂=5e-5，短程 ≤20 步 | `T3(FT)` |

实施判断：
1. 级 CFL 比→Δt：C = ratio·C*，Δt = CΔz/|W|；ratio=0 表示零速（W=0），Δt 取 C*Δz/|W_sampled|。|W| ∈ U[1e-6,1e-2] m/s，符号 ±（向上时初值改为中心 −300 m 的镜像 Gaussian）。步数 U{1..20}。SplitRK3 每级 Δt 的最大值等于 Δt（β 末级为 1），故 C 按整步 Δt 定。
2. 两个强迫的分配 λ ∈ U[−1,2]（W₁=λW、W₂=(1−λ)W）；零 Forcing 用离散形式返回 0。
3. GM 源比较：二维 x–z 网格 32×32（100 km × 1 km，Bounded, Flat, Bounded）、PrescribedVelocityFields() 冻结输运为零，κ_skew=1000、κ_symmetric=nothing；对照模型对 b 与 D 都加同一 AdvectiveForcing(u=uₑ, w=wₑ)（uₑ,wₑ 为 GM 模型初态涡致速度拷贝），使两模型「其余相同」。Δt₀ = 0.1·min(Δx/|uₑ|, Δz/|wₑ|)，再 /2、/4；「收敛到 0」判据：Δt₀/4 的 差/Δt 不大于 Δt₀ 的一半，或已 ≤ T/Δt₀。
4. 元组：零系数闭合用 `ScalarDiffusivity(ν=0, κ=0)`；推进 3 步、Δt=3600 s；比较 D 数组，容差 T·max(1,‖D0‖∞)。
5. taper：L_x=L_y=100 km、H=1000 m；κ = 1000·(1+0.5(z+H)/H)（κ′=0.5 m/s）；b 的竖直梯度按 b_z 情形取 {−1e-6, 0, m/2, m, U[1e-8,1e-3]}，水平梯度按目标「下限后坡度」S 设 b_x=−max(b_z,m)·S_x（下限后 ≤0 时 b_x 取 −1e-6·S_x，预期速度全零）；角度取八个象限方向或随机；参数精度 Float64 或 Float32 各半。内部比较区离边界 2 格。特征速度 κ′S_m 作零分支的绝对尺度。
6. I4：拉伸网格 z_k = −H(1−k/Nz)^{1.5}；b = 1e-5 z + amp·1e-3·sin(2πx/L)cos(2πy/L)·e^{z/300}；常数示踪物用 Centered 默认平流，5 步，Δt 由 CFL 0.1 定。
7. 每次单柱运行都检库存（I2）；质心只在 PR 原设置必测点检（计划另写的大域/40-80-160 加密质心检查未做，见 §4）。

## §4 未实现项（偏离候选）

- **CUDA 全部（I6 的 GPU 部分、CPU/GPU 一致性）**：compute node B 无 CUDA 设备。Float32/Float64 的 CPU 部分已覆盖。
- I1/I2：大域与 40/80/160 层加密的质心检查；UpwindBiased(1) 与独立组装的守恒三对角更新对照（AIVA 显式比例取决于 `clock.last_Δt` 初始化细节，独立组装难以保证与库内同一定义，未做）；「输运单独」的二维胞状环流分支；背景输运＋漂移＋强迫同时开、与细步参考比较；拉伸网格、`ImmersedBoundaryGrid(GridFittedBottom)`、可变体积（`MutableVerticalDiscretization`/`ZStarCoordinate`）上的库存检查。
- I3：constant / face-field κ 表示（只做了函数 κ）；DiffusiveFormulation 与 VerticallyImplicit 选项、κ_symmetric 正值、按示踪物的对称系数；三维变梯度情形对两种交错的独立插值；拒绝按示踪物 NamedTuple skew 系数的负检查。
- I4：在浸没边界外围的边界通量；可变网格；advective 与 diffusive 形式在 16³/32³/64³ 的收敛（< 2%）。
- I5：可变网格（z-star，`set!(model; η=…)`）与浸没底的部分（本质部分）；只做了固定网格「零回退」。
- 这些多数是另建网格/闭合组合即可实现的工作量问题，不是环境不可得；是否补做请主会话裁定。

## §5 耗时与内存估计（推测）

- 依据：同环境原生重放 job 9cd55d6e 中 Oceananigans 两份约 2.5 分钟、峰值约 3.5 GiB。
- 第一路：A 三次单柱、每次 ≤120 步（40 层）→ 约 0.3–1 s；B 约 0.2–0.5 s；C 建 16³ 模型并逐面调内核 → 约 0.5–1 s。10⁴ 点约 1–2.5 小时，预算 3000 s（预计 N≈4000–10⁴）。
- 第二路：两次单柱 ≤20 步 → 约 0.2 s；GM 约 0.5 s → 10⁴ 点约 40–60 分钟，预算 2400 s。
- 第三路：单柱矩阵约 0.3 s；GM 源/元组约 0.5–1 s；taper 三个 16×16×32 模型约 1 s；I4 约 1 s → 10⁴ 点约 2 小时，预算 3000 s。
- 编译：6 个格式×时间推进组合、两种精度、多种闭合与网格拓扑，首次合计估 10–15 分钟。
- 峰值内存估 4–6 GiB，申报 10 GiB。全量 job 估 9000 s，硬上限 14400 s。

## §6 冒烟后修改记录

（空）
