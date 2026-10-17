# ca4637 实施日志（ClimaAtmos#4637，ΔQ1 统一抽查）

> 目的：记录三路检查计划逐句对应到 `spotcheck.jl` 哪段代码，以及全部实施判断。
> 结构：§0 共用 → §1–§3 三路 → §4 未实现项（偏离候选）→ §5 耗时与内存 → §6 冒烟后修改。
> 更新规则：冒烟后只修脚本错误，记入 §6。引用：协议 `notes/20261017-prereg-dQ1-uniform-spotcheck.md`；共用库 `../jl_common/spotcheck_common.jl`。

## §0 共用

- 环境：compute node B `env/ca`，ClimaAtmos ba2bf122（**晚于本 PR 合并提交 a8efce5a**；与原生重放同一环境，按协议 §1 用它）。有原生环境。单线程。
- 入口：`--tracks 1,2,3 --n N --outdir DIR [--budgets …] [--precision Float64|Float32]`（由 `--track/--out` 改，一进程三路省编译）。种子 `20261017 + 2×10 + 路号`。
- **分两个进程**（同一 job 命令里串行）：`--precision Float64` 只跑 F64 分支，写 `ca4637_t<路>.json`；`--precision Float32` 只跑 F32 分支，写 `ca4637_f32_t<路>.json`。理由：F64、F32 各需约 4 个 ClimaAtmos 配置，一个进程装 8 个配置的编译产物估计超过 compute node B 约 24 GiB 的任务上限（此前 3 个配置实测峰值约 11 GiB）。两进程种子相同、分支集合不同，故抽样序列不同；两份 JSON 合起来是该路的全部点。**这是实施判断，需主会话认可。**
- 配置（照 PR 测试夹具）：
  - plane：`ConstantBuoyancyFrequencyProfile`、Schar、x_max=100 km、8 谱元、z_elem=20、dz_bottom=500 m、`disable_surface_flux_tendency=true`；干模式（`microphysics_model=dry`）用于质量守恒，1M（`implicit_microphysics=false`）用于 χ 检验。
  - EDMF：`DYCOMS_RF01` 单柱、prognostic EDMFX、1M、z_max=1.5 km、30 层、dt=120 s、ARS222、`toml/prognostic_edmfx.toml`，其余开关照 `tracer_mass_consistency_tests.jl`。
- 「已推进状态」：每个基态先推进若干步（plane 5 步、EDMF 10 步，同 PR 测试），再逐步存 20 个快照；每点抽一个快照编号，装回状态后调 `CA.set_precomputed_quantities!` 使缓存与该快照一致，再按计划改 Y、调倾向函数（与 PR 测试一样，改 Y 后不重算缓存）。这样拒绝点可凭快照编号单独重放。

## §1 第一路（track1.md）

| 计划原文 | 代码 |
|---|---|
| 仿 advection_tests 与 tracer_mass_consistency_tests 用 `CA.AtmosConfig` 建模拟（plane、Schar、Linear、z_max 21 km） | `plane_cfg`、`edmf_cfg` |
| A：读 `p.precomputed.ᶠu³` 顶层与 `CA.compute_ᶠuₕ³`，检比值 | `A_u3_{bottom,top}_le_100eps_uh3` |
| 再用干模型、关地表通量、无微物理推进 1000 步，每步 `sum(Y.c.ρ)` 记 M，相对漂移 ≤1e-12（F64）且不随时间线性增长；对照 SLEVE | 分支 `A_mass_{Linear,SLEVE}_*`：`A_mass_rel_drift_max`、`A_mass_drift_linear_trend_over_1000_steps` |
| B：在已推进状态上 ρχ ← cρ、χʲ ← c（c = 1 与 1e-3），分别调 `implicit_tendency!`、`explicit_vertical_advection_tendency!`、`edmfx_sgs_mass_flux_tendency!`，比较 Yₜ.ρχ 与 c·(相应 ρ 路径)，容差 100ε·max|ρ 倾向| | 分支 `B_grid_paths_*`（隐式、显式）、`B_sgs_*`（SGS） |
| 再做推进版：均匀 χ 初值、关微物理源汇，推进 N 步后 max|χ−c|/c ≤ N·100ε | 分支 `B_timestep_*`（**实施偏离，见判断 5**） |
| C：同一调用比较 `parent(Yₜ.c.ρ) == parent(Yₜ.c.ρq_tot)` | 分支 `C_sgs_rho_eq_rhoqtot_*` |
| 工况域：地形 0–3 km、z_top 20–45 km、风 0–30 m/s、Δt 1–120 s；EDMF 单柱 Δt 60–120 s；F64 与 F32 | 抽样 |

实施判断：
1. A 每点按抽到的 z_top ∈ U[20,45] km、Δt ∈ U{1…120} s（整数秒，配置只收「N secs」）新建一个干模拟，风速把初值的均匀 10 m/s 按 U/10 缩放（U ∈ U[0,30]）。Δt 大时水平声波显式 CFL 超限会发散 → 记失败（不算违反）。
2. **Schar 山高不可经配置改**（`ScharTopography` 默认 h_max=25 m，未暴露为配置项），工况域的 0–3 km 无法抽样，固定默认；Agnesi 山未做（换地形类型需重编译）。偏离候选。
3. F32 漂移容差计划未给：用协议默认相对 1e-6。「不随时间线性增长」判为：相对漂移对步数的最小二乘斜率 × 1000 步 ≤ 同一容差。
4. B：照 PR 测试把全部水物种（q_tot、q_lcl、q_icl、q_rai、q_sno）同时设成 cρ（EDMF 中各子域 χʲ 也设 c）。「相应 ρ 路径」：隐式用 `Yₜ.ρ`（隐式倾向的 ρ 行）；显式用 PR 测试的 −ᶜadvdivᵥ(ᶠinterp(ρJ)/ᶠJ·ᶠu³)；SGS 用同一调用的 `Yₜ.ρ`。SGS 的「max|ρ 倾向|」取**未改动状态**下 SGS 调用给出的 max|Yₜ.ρ|（改成均匀后它理应为 0，不能作尺度）。容差不随 c 缩放（照原文）。
   - 推测：c=1e-3 时隐式路径可能出现差 (1−c)×沉降项（库把水的沉降通量同等加到 ρ 与 ρq_tot），这会是「计划写法 vs 库约定」问题，交裁定。
5. **推进版无法照原文**：1M 没有关闭微物理源汇的开关，χ≡1 或 c=1e-3 均匀的凝结物会被微物理改变。改为：q_tot 设成极小常数 c=1e-7（各处次饱和，凝结/沉降源为零）、凝结物为 0，重算缓存后推进 N=20 步（N 计划未给），检 max|q_tot−c|/c ≤ N·100ε。这是偏离，记 §7 候选。
6. 海绵层：默认关（计划提醒须关闭或限制在海绵层以下）。

## §2 第二路（track2.md）

| 计划原文 | 代码 |
|---|---|
| I1：`CTS.step!` 推进，每步 ClimaCore `sum(Y.c.ρ)` 算总质量，100–1000 步后比较，容差 1e3·N·eps·M(0)；顶面每步 ≤100·eps·max|uₕ³| | 分支 `I1_mass_{Linear,SLEVE}_*` |
| 建议另测 Sleve、更高山、dt 1–10 s、Float32 | Sleve 与 F32 做了；dt ∈ U{1..10}；山高见 §1 判断 2 |
| I2：初始化后 ρq_* ← ρ、子域 q_* ← 1，关微物理与海绵，推进 N=50 步，检 max|ρq_x/ρ − 1|，容差 1e-10（F64） | 分支 `I2_timestep_{plane,edmf}_*`（**同 §1 判断 5 的偏离**） |
| I3：调 `edmfx_sgs_mass_flux_tendency!` 得 Yₜ，按柱竖直积分，相对 max|Yₜ.ρq_tot|·Σ JΔz ≤ 1e-12；且 Yₜ.ρ − Yₜ.ρq_tot = 0 | 分支 `I3_sgs_column_*` |

实施判断：
1. I1 用干模式（不变量要求关闭降水沉降和所有源汇）；PR 夹具的其余设置不变（z_max 21 km、风 10 m/s）；步数 U{100…1000}。「每步顶面检查」：每步后调 `set_implicit_precomputed_quantities!`（状态滤波幂等）再读 ᶠu³，取全程最坏比值与 1 比。
2. I2 推进版：同 §1 判断 5（只 q_tot、c=1e-7）；EDMF 用一套关地表通量的同型配置（`disable_surface_flux_tendency=true`，否则地表水汽通量本身破坏均匀），其余同夹具。F32 容差 1e-4（计划只给 F64；按 eps 比例放大，实施判断）。
3. I3：柱积分用 ClimaCore `sum`（WJ 加权）；Σ JΔz 用全 1 场的 `sum`。对 ρ、ρq_tot 及四种凝结物都检；第二式按「= 0」要求逐位相等。

## §3 第三路（spec.md §5）

| 计划原文（摘） | 代码 |
|---|---|
| 从 PR 两个夹具起步；关物理源测输运 | `plane1M`、`edmf` 基态 |
| I1：给两种符号的风后调 `CA.set_implicit_precomputed_quantities!`，在 `CA.half` 与 `nlevels+half` 读 ᶠu³ 与每个 ᶠu³ʲs，残差 ≤ 100ε·max(|uₕ³|,|g³³u₃|)，零尺度时恰为 0；扭曲回归须有非零水平投影 | 分支 `I1_impermeability_*`（Linear、SLEVE、EDMF 单柱）；U=0, ±10 为必测点 |
| I2：零缓冲；密度参照 −advdiv(ρu³)；对各被动场调 `explicit_vertical_advection_tendency!`，总水用 `CA.vertical_transport`；C=0,1,1e-6,0.01；‖A_χ−CA_ρ‖∞ ≤ 200εS·|C|；人为非零边界速度下 D₀ 仍结构上零通量 | 分支 `I2_constant_tracer_*` |
| I3：各子域与平均都设常数后调 `edmfx_sgs_mass_flux_tendency!`，500εS_SGS 内为零；不等速变面积；等速时非均匀标量也为零；关闭时缓冲不变 | 分支 `I3_null_sgs_*`（三种模式） |
| I4：非均匀湿度、非零 SGS 水通量，零缓冲出发 Yₜ.ρ 与 Yₜ.ρq_tot 逐位相等；相同非零初缓冲重复 | 分支 `I4_water_air_*` |
| I5：用 ClimaCore 体积分，|ΣVT| ≤ 500NεΣV|T|；20/40/80 层与拉伸网格 | 分支 `I5_closed_column_*` |

实施判断：
1. I1：plane 用干模式 Linear/SLEVE 夹具（Δt=1 s、z_max 21 km），每点从初值拷贝重来并把水平风按 U/U₀ 缩放（U ∈ U[−100,100]）；EDMF 用快照状态同样缩放。|g³³u₃| 取 |u³−uₕ³|。
2. I2：S 取 max|A_ρ|（计划另写「与面通量绝对值和/体积取大」，未单独算这一项，实施简化）；「configured reconstruction」取 `p.atmos.numerics.energy_q_tot_upwinding` 并另测中心格式；人为边界速度：把缓存 ᶠu³ 顶层设为 1e-2 后，检显式凝结物倾向的体积分 ≈ 0（容差同 I5）。数密度常数（2M）不适用于 1M，未测。
3. I3：S_SGS 取未改状态下 SGS 调用的 max(|Yₜ.ρq_tot|,|Yₜ.ρ|)，能量行用 max|Yₜ.ρe_tot|；条件数因子 max(1,ρ/m₀) 未计算与报告（简化）。「等速」把缓存的 ᶠu³ʲs 与 ᶠu³⁰ 设成 ᶠu³。「关闭」把 `edmfx_model.sgs_mass_flux` 换成 false（复制结构替换字段）、缓冲预置为 1 后检逐位不变。
4. I4：非零初缓冲为每点独立种子的 U[0,1) 随机场，ρ 与 ρq_tot 行相同。另记 `report_I4_flux_nonzero`（只记录，保证非空）。
5. I5：plane 检 ρ 平流、q_tot 中心输运、四种凝结物显式平流；EDMF 检 SGS 的 ρ 与各水物种行；EDMF 的层数取 {30（夹具）, 40, 80}、拉伸开/关（计划写 20/40/80，单柱夹具本身是 30 层，取 30 代替 20）。

## §4 未实现项（偏离候选）

- 第一路：地形高度抽样（配置不暴露）、Agnesi；推进版 χ 检验按原文（见 §1 判断 5，已用可行变体代替）。
- 第二路：「更高山」；I2 原文的 χ≡1、全部水物种推进 50 步（同上）。
- 第三路：平坦 plane、3-D 扭曲 box；`turbconv` 无/EDOnlyEDMFX、`updraft_number=2`；0M、2M、2MP3 微物理与数密度示踪物；自定义被动示踪物；`tracer_upwinding`/`energy_q_tot_upwinding`/SGS 重构/`edmfx_tracer_upwinding` 各选项扫描（每个都要新配置重编译）；20 次带随机密度与面积的种子实现（用 20 个快照代替）；I3 的环境小面积回退阈值附近取样；I6 Jacobian 全部（`CA.Jacobian`、`update_jacobian!`、`first_column_block_arrays`、方向差分）。
- 理由：配置数与编译内存/时间；I6 工作量大且依赖内部块名。是否补做请主会话裁定。

## §5 耗时与内存估计（推测）

- 依据：同节点原生重放 job 9cd55d6e：ClimaAtmos#4637 部分（plane 干/1M 与 EDMF 单柱、EDMF 推进 10 步）约 2.5 分钟，RSS 峰值约 11 GiB。
- 第一路：A 每点建模拟（约 10–30 s）＋1000 步（plane 8×4 列×20 层、ARS343 隐式，估每步 20–60 ms）→ 单点约 0.5–1.5 分钟；B/C 每点约 0.2–2 s。A 占 1/3 的点 → 10⁴ 点需数十小时，按预算截断：F64 2400 s、F32 1200 s（预计各约 60–150 点）。
- 第二路：I1 每点建模拟＋100–1000 步 → 约 0.3–1 分钟；I2 推进 50 步约 5–15 s；I3 约 0.2 s。预算 F64 2000 s、F32 1000 s。
- 第三路：多为单次倾向调用（约 0.1–0.5 s），I5 的 EDMF 层数变体首次各需建模拟＋30 步。预算 F64 1500 s、F32 800 s。
- 编译：每个进程约 4 个配置，各 1–2 分钟。
- 峰值内存：每个进程估 13–16 GiB（3 个配置 11 GiB 外推），申报 20 GiB（两进程串行，峰值取较大者）。全量 job 估 10 800 s，硬上限 14 400 s；冒烟（F64 各路 400 s、F32 各路 300 s）估 3600 s。

## §6 冒烟后修改记录

（空）
