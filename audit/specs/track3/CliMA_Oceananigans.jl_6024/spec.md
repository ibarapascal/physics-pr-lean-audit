## CliMA/Oceananigans.jl#6024 Fix CATKE's bottom mixing length coefficient Cᵇ

- What changed (1–3 sentences, naming files and formulas):
  `src/TurbulenceClosures/turbulence_closure_implementations/TKEBasedVerticalDiffusivities/catke_mixing_length.jl` 将默认底边界系数 `Cᵇ` 从 `0.28` 改为 `0.967`，并添加中性相似层关系 κ = C[(Cˡᵒu)³/CˡᵒD]^(1/4) 的说明；`docs/src/models/turbulence_closures.md` 同步展示值。此次没有改离散算子或增加测试。以下是条件性闭合模型不变量，不是任意海流的普适守恒定律。

- Invariants (mathematical statements; write them as checkable equalities/inequalities, give variables, sign conventions, domain of validity):
  约定 z 向上，平底 z = −H、静止表面 z = 0，d = z + H > 0；e 是单位质量 TKE（m²/s²），b 是浮力，N² = ∂b/∂z，S = ∂U/∂z > 0。令 a = Cˡᵒu = 0.361、q = Cˡᵒc = 0.369、r = Cˡᵒe = 7.863、D = CˡᵒD = 1.604、B = Cᵇ = 0.967、A = Cˢ = 1.131。ν、Kc、Ke 分别是动量、普通示踪剂及 TKE 的垂向扩散系数；ε ≥ 0 是耗散率，TKE 方程中的贡献为 −ε。

  1. **Primary — 底部中性平衡层的 κ 标定。** 在 N² = 0、无对流、底距离限制主导、e 高于数值下限、无扩散率截断的局地相似层，ℓ = Bd，ν = aℓ√e，ε = D e^(3/2)/ℓ。若剪切生产 P = νS² 与 ε 平衡，且正的运动学应力大小 u★² = νS，则必须有
     e = u★²/√(aD)，κ_eff = B(a³/D)^(1/4)，ν = κ_eff u★d，S = u★/(κ_eff d)。
     因而 κ_eff⁴D = B⁴a³，且默认参数应满足 |κ_eff − 0.4| ≤ 5×10⁻⁴。最后一个容差允许 B 仅保留三位小数；不能要求 κ_eff 与十进制 0.4 浮点严格相等。对应连续速度差 U(d₂) − U(d₁) = (u★/κ_eff) log(d₂/d₁)，但这里不要求整个有限水柱或最底层网格严格服从对数律。

  2. **Secondary — 中性底部混合的距离律、比例与非负性。** 在均匀正 e、N² = 0、零表面浮力通量、未触发 H 长度上限和扩散率上限时，底部控制区有
     (ν, Kc, Ke) = Bd√e (a, q, r)，ν ≥ 0，Kc/ν = q/a，Ke/ν = r/a。
     对任意两个内部面 d₁,d₂，若均处于该分支，ν(d₂)/ν(d₁) = d₂/d₁。在固定 e 的两组配置 B₁、B₂ 中，若两者都受底部控制，ν₂/ν₁ = Kc₂/Kc₁ = Ke₂/Ke₁ = B₂/B₁；特别是新旧默认系数的比为 0.967/0.28。一般中性水柱的参考式是 Kα = √e min[H, cα min(Bd,A(H−d))]，cα ∈ {a,q,r}，并再施加相应最大扩散率上限。该一般式解释为什么 Ke 在较远离底部处可能先饱和，不能在那里强制要求三个系数同比变化。边界面被实现显式置零，不用于内部面比例测试。

  3. **Secondary — 同一底长度必须一致地控制 TKE 耗散。** 在 N² = 0、无剪切、无表面强迫、初始均匀 e₀ > minimum_tke 的水柱内，初始 TKE 扩散为零。对于底部控制且 ℓ/D < H 的内部单元中心：
     ε(d) = D e₀^(3/2)/(Bd) > 0，lim(Δt→0⁺)[e(d,Δt)−e₀]/Δt = −ε(d)。
     固定 e₀、d 时 ε(B₂)/ε(B₁) = B₁/B₂。因此同一个 B 增大既增加给定 TKE 下的涡黏性，又减小给定 TKE 下的局地耗散；不能把“底拖曳增强”理解为所有反馈后的流态均有同一响应。

- Sources of the invariants (first principles / literature formula / consistency with other implementations in the library / dimensions):
  第一项来自上述涡黏闭合、正定耗散、局地 P = ε 及恒定应力假设的代数消元；κ = 0.4 是此 PR 采用的经验刚壁标定目标，不是由守恒律独立推导的常数。PR 输入将关系归于 Wagner et al. (2025), *Formulation and calibration of CATKE, a one-equation parameterization for microscale ocean mixing*, §3.1.3, Eq. (35), DOI 10.1029/2024MS004522；本审计只在允许的 GitHub 资料中核对了这一归属，未独立读取期刊全文，也不以 PR 自述替代证明。第二、三项是同一模型的构成关系及量纲一致性：ℓ 单位 m，ν 单位 m²/s，ε 单位 m²/s³。
  合并版源码依据：[混合长度与稳定函数](https://github.com/CliMA/Oceananigans.jl/blob/24de7026893af38df065920fbbf80f05024fd1b9/src/TurbulenceClosures/turbulence_closure_implementations/TKEBasedVerticalDiffusivities/catke_mixing_length.jl)、[扩散系数与默认值](https://github.com/CliMA/Oceananigans.jl/blob/24de7026893af38df065920fbbf80f05024fd1b9/src/TurbulenceClosures/turbulence_closure_implementations/TKEBasedVerticalDiffusivities/catke_vertical_diffusivity.jl)、[耗散长度 ℓ/D 与耗散率](https://github.com/CliMA/Oceananigans.jl/blob/24de7026893af38df065920fbbf80f05024fd1b9/src/TurbulenceClosures/turbulence_closure_implementations/TKEBasedVerticalDiffusivities/catke_equation.jl)、[TKE 下限、Ri=0 分支与边界掩膜](https://github.com/CliMA/Oceananigans.jl/blob/24de7026893af38df065920fbbf80f05024fd1b9/src/TurbulenceClosures/turbulence_closure_implementations/TKEBasedVerticalDiffusivities/TKEBasedVerticalDiffusivities.jl)。

- Common operating domain (range of each input, with the basis):
  数学域为 H > 0、0 < d < H、B,a,D > 0、有限 e > minimum_tke，且 Bd < A(H−d)。三种扩散系数一起检查时还需 rBd < H；排除浸没边界、湿干单元、波浪强迫、旋转、平流和其他闭合叠加。N² 精确设为零，表面浮力及动量通量设为零；不能用略为负的 N² 代替零，因为负 Ri 使用不同稳定系数。
  建议随机样本：H ∈ [10,1000] m；内部采样 d/H ∈ [0.01,0.05]；e₀ ∈ [10⁻⁶,10⁻²] m²/s²；B 测试取旧值 0.28、新值 0.967，或 [0.2,1.2] 内随机值。这些是便于隔离分支和避开默认 e 下限 10⁻⁹ 的审计范围，不是文献宣称的全适用范围。第一项取 u★ ∈ [10⁻³,10⁻¹] m/s，并由平衡式生成 e；U 的参考常数可取零。相似律只用于远离黏性/粗糙子层的湍流区域；本闭合不解析分子黏性子层。均匀网格建议 Nz ≥ 512，排除接触物理边界的差分模板；所有比较使用真实网格坐标。

- Numerical check method (which public API of the merged commit, what to compute, what to compare against, suggested tolerance and why; if finite differences are needed, the step-size strategy):
  **公共入口与观测量。** 使用 `using Oceananigans`，从 `Oceananigans.TurbulenceClosures` 导入 `CATKEVerticalDiffusivity`, `viscosity`, `diffusivity`。建立 `RectilinearGrid(CPU(), Float64; size=Nz, z=(-H,0), topology=(Flat,Flat,Bounded))`，再建立 `HydrostaticFreeSurfaceModel(grid; closure=CATKEVerticalDiffusivity(Float64), buoyancy=BuoyancyTracer(), tracers=(:b,:c), coriolis=nothing, momentum_advection=nothing, tracer_advection=nothing)`。**不要在 tracers 列表中提供 :e**：合并版自动添加它，显式提供反而报错。用 `set!(model, b=0, c=0, e=e₀, u=..., v=0)` 设置状态；此公共入口会更新闭合场，无需调用私有 kernel。用 `Array(interior(viscosity(model)))`、`Array(interior(diffusivity(model, Val(:c))))`、`Array(interior(diffusivity(model, Val(:e))))` 读取三个面场，用 `znodes(grid, Face())`/`znodes(grid, Center())` 取得相应位置。所有数组索引需与面/中心坐标匹配。
  已核对入口源码：[set!](https://github.com/CliMA/Oceananigans.jl/blob/24de7026893af38df065920fbbf80f05024fd1b9/src/Models/HydrostaticFreeSurfaceModels/set_hydrostatic_free_surface_model.jl)、[状态更新调用闭合计算](https://github.com/CliMA/Oceananigans.jl/blob/24de7026893af38df065920fbbf80f05024fd1b9/src/Models/HydrostaticFreeSurfaceModels/update_hydrostatic_free_surface_model_state.jl)。上述方案是待执行规范，本次未运行 Julia。

  **第一项。** 先从公开构造出的 `closure.mixing_length` 与 `closure.turbulent_kinetic_energy_equation` 读取系数，验证 κ_eff⁴D = B⁴a³ 及 |κ_eff−0.4| ≤ 5×10⁻⁴。为避免仅重复参数公式，再用独立目标 κ₀ = 0.4 制造状态：e₀ = u★²/√(aD)，在各中心设置 U(d) = (u★/κ₀)log(d/d_ref)，d_ref 可取 H。在选定内部面比较实际 ν/(u★d) 与 0.4，绝对误差 ≤ 5×10⁻⁴；比较 νS_h/u★² 与 1，容许系数舍入相对误差 1.25×10⁻³ 加空间误差。S_h 用相邻中心速度差除以中心间距。以 Nz、2Nz、4Nz 加密并固定物理采样带，要求对数曲线导数误差约按 O((Δz/d)²) 下降，最终空间误差 < 10⁻⁴；不要把固定前三个单元当作固定物理点。还可在相同内部位置检查 P/ε → 1，其中 P = νS_h²、ε 采用第三项独立检验过的耗散式。此快照检查检验构成关系，不宣称整柱已达到稳态。旧 B = 0.28 应明显不能通过 κ≈0.4 的检查。

  **第二项。** 以均匀 e₀、u=v=0 做快照，比较三个面场与明确的分支参考式及比值。Float64 使用 `rtol=1e-11`，绝对容差取 `1e-13 × 本次参考场最大值`；该测试没有时间积分误差，均匀 e 消除了 √e 插值误差。Float32 可用 `rtol=1e-5`。如需自定义旧系数，使用已公开文档化但未导出的 `Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities.CATKEMixingLength(Cᵇ=...)` 作为 `CATKEVerticalDiffusivity(mixing_length=...)` 的构造参数，其余默认值保持相同。另取表面限制主导的点，对比新旧 B 应相同；取分支转换两侧用一般 min 公式检查连续性，不要求导数连续。

  **第三项。** 使用均匀 e₀、u=v=0 的新模型，从 `Oceananigans.TimeSteppers` 导入 `time_step!`，执行一个很小的 Δt 后读取 `interior(model.tracers.e)`。对每个被测中心计算 R_h = −(e_new−e₀)/Δt，比较 D e₀^(3/2)/(Bd)。默认垂向隐式路径仍应在 Δt→0 时得到同一个极限；不能把私有 explicit `dissipation` 对隐式方案返回零误认为物理耗散为零。起始 Δt 取所有网格中心最短耗散时间 Bd/(D√e₀) 的 10⁻³ 倍，并限制 Δt·max(Ke)/min(Δz)² ≤ 10⁻³；随后逐次减半，每次从相同初态的新模型开始。要求相对误差降到 10⁻³ 以下且呈收敛趋势；在浮点相减噪声开始增大前停止，必要时用相邻两级线性外推。比较新旧 B 的 R_h 比例及 R_h ≥ 0。此法包含有限 Δt 后产生的 TKE 梯度扩散，但其初始贡献为零，误差在 Δt→0 消失。

- Do the PR's own tests already check this invariant (yes/partial/no, naming test files and assertions):
  **No，三个目标均没有 PR 新增断言。** 输入 diff 只含源码默认值/说明和文档。合并版 [test/turbulence_closures/closures.jl](https://github.com/CliMA/Oceananigans.jl/blob/24de7026893af38df065920fbbf80f05024fd1b9/test/turbulence_closures/closures.jl) 有相邻覆盖：`run_catke_tke_substepping_tests` 检查 Euler/AB 步进恒等式；`run_time_step_with_catke_tests` 检查构造、自动 e 示踪剂规则、表面边界条件与能否步进；`run_catke_closure_tuple_surface_tke_flux_tests` 检查表面通量与 closure tuple 等价性；另有面场位置及 CFL 数值断言。这些不检验底部 κ=0.4、距离比例或底部耗散的物理参考值。[validation/vertical_mixing_closures/test_single_column_model.jl](https://github.com/CliMA/Oceananigans.jl/blob/24de7026893af38df065920fbbf80f05024fd1b9/validation/vertical_mixing_closures/test_single_column_model.jl) 是带稳定分层和表面强迫的步进/计时例子，没有这些不变量的断言。

- Do you expect the change to satisfy it (satisfies / may not satisfy / uncertain), with one or two sentences of reasoning (speculation, not a conclusion):
  **satisfies（预测，未数值验证）**：在限定分支内，合并源码的 ν = aBd√e 与 ε = D e^(3/2)/(Bd) 代数上给出含 a³ 的四次根关系，新 B 应在舍入容差内恢复 κ=0.4，且相同 B 同时用于各混合长度。对一般三维、分层、带底拖曳反馈的模拟响应仍属 uncertain；本规范不把局地标定外推成全局拖曳或耗散一定增加的结论。

- Environment: merge commit SHA, language and dependencies needed, estimated setup difficulty (low/medium/high):
  合并 SHA `24de7026893af38df065920fbbf80f05024fd1b9`；Julia，库版本 `0.113.5`，`Project.toml` 的 Julia compatibility 为 `1.10`。未来执行需该提交的 Oceananigans 及项目依赖（包括 KernelAbstractions、Adapt、OffsetArrays、StaticArrays、FFTW、MPI、JLD2 等）；CPU 一维列测试不需要 CUDA、GPU 或绘图库。环境准备难度 **medium**：Julia 包解析及二进制依赖/预编译可能耗时，但计算本身很小。[合并版 Project.toml](https://github.com/CliMA/Oceananigans.jl/blob/24de7026893af38df065920fbbf80f05024fd1b9/Project.toml)。本次仅读取输入与 GitHub 官方仓库公开源码，未读取评论、issues、fork 或其他审计轨道，未安装、编译、运行项目代码。

- Formalization (Lean) feasibility: can the invariant be written as a proposition over finite-dimensional reals, and what are the main obstacles:
  **可以。** 固定一个或有限多个网格点，以 H,d,e,a,D,B,u★,ν,S,ε 为正实数，写出 min 分支假设及构成关系；由 P=ε、νS=u★² 证明 κ⁴D=B⁴a³ 与各比例关系。可用正根变量及多项式等式避免先引入一般实幂，三位小数视为有理数并用区间界证明 κ 接近 0.4。第二项的非负性、分支比例与 min 连续性也适合直接证明。第三项的初始连续模型关系可有限维化为半离散 ODE；若证明实际 Julia 更新，主要障碍是面/中心插值、边界掩膜、隐式线性求解、时间分裂及浮点误差的语义对应。Lean 代数证明不能证明经验目标 0.4 对所有真实海底都成立，也不能从局地平衡假设推出整个水柱必然收敛到该平衡。
