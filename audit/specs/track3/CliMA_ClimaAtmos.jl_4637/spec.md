## CliMA/ClimaAtmos.jl#4637 Fixes problems near upper boundary and tracer-mass continuity

- What changed (1–3 sentences, naming files and formulas):
  `src/cache/precomputed_quantities.jl` 将顶边界条件改为协变分量 `u₃ = -uₕ³/g³³`，使总逆变分量 `u³ = uₕ³ + g³³u₃` 为零；`src/prognostic_equations/implicit/implicit_tendency.jl` 的连续方程改用与示踪物相同的零边界通量散度 `ᶜadvdivᵥ`。`advection.jl` 恢复所有被发现示踪物的网格平均垂直平流，`edmfx_sgs_flux.jl` 将 SGS 修正改为 `ρᵏaᵏ(u³ᵏ−u³)(χᵏ−χ)` 的通量散度，并把总水 SGS 质量倾向同样加到湿空气密度；`implicit/manual_sparse_jacobian.jl` 改为匹配中央重构的上升气流近似线性化，省略环境项。`remaining_tendency.jl` 仅给现有凝结物海绵层添加 TODO，并未修复该项。

- Invariants (mathematical statements; write them as checkable equalities/inequalities, give variables, sign conventions, domain of validity):
  1. **Primary — 不可穿透边界与封闭垂直质量收支。** 令竖直坐标向上，`u₃` 为协变速度、`uₕ³` 为水平风对逆变竖直方向的投影，正定网格度量满足 `g³³ > 0`。非 prescribed-flow 模式下，完成边界过滤后，顶面及底面应满足 `u³_b = uₕ³_b + g³³_b u₃_b = 0`；共享水平速度的每个预报 EDMFX 子域也应满足 `u³ʲ_b = 0`。对单独的网格平均垂直质量平流 `Aρ = −D₀ Fρ`，其中 `Fρ_f = I(ρJ_c)_f u³_f/J_f`、`D₀` 在两端指定零通量、`I` 为面插值，应有 `Σ_i V_i Aρ_i = 0`，`V_i > 0` 是包含几何及水平求积权重的物理单元体积。零通量必须在散度算子上结构性成立，即使人为给缓存的端面速度一个非零值。这里不包括降水自由流出、地表源汇或其他倾向。
  2. **Primary — 常量示踪物与质量平流一致，SGS 不从零对比产生通量。** 对任一具体示踪物 `χ`，设 `ρχ_i = Cρ_i`，`C` 为常数，则网格平均垂直平流应逐点满足 `A(ρχ)_i = C Aρ_i`，包括首末单元。若所有上升气流及诊断环境也满足 `χᵏ_i = χ_i = C`，则单独 SGS 修正 `Sχ = −D₀ Σ_k Fχᵏ` 应为零；连续表达为 `Fχᵏ = ρᵏaᵏ(u³ᵏ−u³)(χᵏ−χ)`，离散表达应按实际重构 `aᵏ(χᵏ−χ)`，不能假定插值乘积等于乘积插值。该性质允许密度、面积比例和速度随高度变化，也允许网格平均与子域密度不同。它是分算子的恒等式：当总水 SGS 通量非零时，不能据此声称任意均匀被动示踪物在完整湿空气方程中始终保持均匀，因为分母 ρ 还接受水质量倾向。
  3. **Secondary — SGS 水质量与湿空气质量成对输运。** `ρ` 为湿空气质量密度，`q_t` 为总水质量分数，`ρ_d = ρ−ρq_t` 为干空气密度。仅打开并提取 `edmfx_sgs_mass_flux_tendency!` 的增量时，应逐点有 `Sρ = S(ρq_t)`，因而 `Sρ_d = 0`；封闭垂直边界还应满足 `Σ_i V_i Sρ_i = Σ_i V_i S(ρq_t)_i = 0`。这表达本库将此 SGS 水项视为水质量重新分布的约定，不主张所有湍流过程的干空气局地倾向都为零。质量分数无量纲，质量倾向为 kg m⁻³ s⁻¹；数浓度示踪物若采用每单位质量粒子数，其倾向单位另行对应，不能与水质量相加。

- Sources of the invariants (first principles / literature formula / consistency with other implementations in the library / dimensions):
  不变量 1 来自固定壁面的 `u·n = 0`、坐标变换和离散散度定理；不变量 2 来自连续性方程与常量标量输运的乘法恒等式，以及协方差通量在标量无子域对比时为零。这两者不依赖某个经验湍流参数。不变量 3 来自湿空气组成恒等式和本库已有水扩散处理的一致性；合并版 [SGS 通量实现及扩散对应项](https://raw.githubusercontent.com/CliMA/ClimaAtmos.jl/a8efce5a0040924094a58739db91a0f3162f779c/src/prognostic_equations/edmfx_sgs_flux.jl) 同时提供两条路径。几何边界、零通量散度及重构依据输入 diff 中的 `precomputed_quantities.jl`、`utils/abbreviations.jl`，以及合并版 [垂直输运实现](https://raw.githubusercontent.com/CliMA/ClimaAtmos.jl/a8efce5a0040924094a58739db91a0f3162f779c/src/prognostic_equations/implicit/implicit_tendency.jl)。这些是第一原理/实现约定推导，不把 PR 中的数值改善描述当作独立验证证据。

- Common operating domain (range of each input, with the basis):
  数学域为有限、正的 `ρ, ρᵏ, J_c, J_f, g³³`，有限速度，非退化网格；`0 ≤ aᵏ ≤ 1`，且诊断环境质量 `ρa⁰ = ρ−Σ_j ρʲaʲ > 0`。随机常规抽查建议采用正定且条件良好的网格，物理空气密度 0.1–1.5 kg m⁻³、温度 230–310 K、水平风 −30–30 m s⁻¹、内部竖直速度 −5–5 m s⁻¹，一条上升气流面积比例 0.01–0.2，环境质量占比至少 0.5；必须由一致状态构造并拒绝不满足热力学/面积约束的样本，不能独立乱配密度、温度、压力。这些是建议的温和大气抽样范围，不是代码声明的有效域。
  水分取 `0 ≤ q_t ≤ 0.03`、各凝结物非负且其和不超过 `q_t`；均匀标量分别用 `C = 0, 10⁻⁶, 10⁻³`，`C = 1` 只作为独立输运算子的代数探针，不能把所有水物种同时设成 1 后当作物理热力学状态重新初始化。数浓度可在有效 2M 状态上取非负常数作补充，不是本规格必需基线。
  有出处的基线来自新增测试：Schar 地形、Linear 网格扭曲、100 km 水平域/8 单元、21 km 高/20 层、底层 500 m、10 m s⁻¹ 初始风、dt=1 s；SGS 基线为 DYCOMS_RF01、1.5 km/30 等厚层、prognostic_edmfx、1M、一条上升气流、ARS222、dt=120 s。优先在这些配置上作小扰动及 20/40/80 层抽查；额外显式推进建议满足实际局部垂直 Courant 数 ≤0.2，此为审计控制条件而非 PR 的稳定性保证。Float64 为主，Float32 复核；SGS 默认中央重构，并补查源码支持的 first_order/third_order/vanleer_limiter（配置允许时）。排除 prescribed flow、源汇、海绵层、相变、沉降和扩散对所检查增量的污染。

- Numerical check method (which public API of the merged commit, what to compute, what to compare against, suggested tolerance and why; if finite differences are needed, the step-size strategy):
  使用 `import ClimaAtmos as CA`，以 `CA.AtmosConfig(Dict(...); job_id=...)` 和 `CA.get_simulation(config)` 创建配置，取 `simulation.integrator.u/p/t`。合并版 [test_helpers.jl](https://raw.githubusercontent.com/CliMA/ClimaAtmos.jl/a8efce5a0040924094a58739db91a0f3162f779c/test/test_helpers.jl) 证实测试 helper 只包装这一入口；无需复制测试 helper。物理配置照输入 diff 的两个测试完整字典。下列 `CA.` 倾向函数和缓存是合并提交可直接调用且测试使用的模块接口，但不宣称它们是跨版本稳定、导出的 API；审计必须固定 SHA。
  对 1：在正常初始化/预计算后，以 `ClimaCore.Fields.level` 读取 `p.precomputed.ᶠu³` 的 `CA.half` 与 `Spaces.nlevels(axes(Y.c))+CA.half` 两端面；用 `CA.compute_ᶠuₕ³(Y.c.uₕ,Y.c.ρ)`、几何度量和 `Y.f.u₃` 独立重建 `uₕ³+g³³u₃`。要求斜坡案例 `max|uₕ³|>0`，再检查网格平均及各 `p.precomputed.ᶠu³ʲs` 端面。将输出缓冲 `T=similar(Y)` 清零，调用 `CA.implicit_vertical_advection_tendency!(T,Y,p,t)`；为单独提取质量平流，使用无沉降的配置，或将此独立算子探针中的沉降速度缓存置零。比较质量行与独立按 `−[J_f F_f]_{下}^{上}/(J_c Δξ)` 组装的有限体积通量差；横向求积权重在体积和中保留。随后仅在复制的缓存中修改端面 `ᶠu³`，直接调用该倾向函数而不重新过滤，要求质量倾向及其体积和不变；此项区分速度过滤与结构性零通量边界。
  对 2：在物理有效的基线缓存上仅修改所测示踪物，令 `Y.c.ρχ=C*Y.c.ρ`，并在需要时把每个 `Y.c.sgsʲs[j].χ` 设为 C。分开清零缓冲并调用 `CA.explicit_vertical_advection_tendency!` 检查凝结物平流、`CA.implicit_vertical_advection_tendency!` 检查总水平流；显式函数不输出质量平流，所以参考值使用上述独立 `Aρ`，而非显式输出的 ρ 行。对总水比较必须先排除沉降质量贡献。SGS 实验先由有效 DYCOMS 状态经 `ClimaTimeSteppers.step!(integrator)` 推进 10 步得到非零循环，记录非均匀基线通量，再在固定密度/速度缓存上设置均匀标量并调用 `CA.edmfx_sgs_mass_flux_tendency!(T,Y,p,t,p.atmos.turbconv_model)`。检查诊断环境 `CA.ᶜspecific_env_value` 的相应值也等于 C，并要求输出 Sχ 为零。固定缓存的代数探针不得随后用作物理时间推进状态。
  对 3：使用未被均匀标量探针改写的有效 DYCOMS 状态，清零 T 后只调用 SGS mass-flux 函数，要求 `max|T.c.ρq_tot|>0`，检查 `T.c.ρ == T.c.ρq_tot`，以及两者带体积权重的和为零；再独立计算干空气倾向差。可从同一状态分别构造不同面积/速度/湿度对比，避免仅验证零流场。
  容差：设 ε 为实际浮点精度。端面检查 `|u³| ≤ 100ε max(|uₕ³|,|g³³u₃|)`；平流恒等式使用 `100ε` 乘以相关局部通量绝对值之和除以单元体积所给出的倾向尺度，而不是接近零的净倾向；SGS 零值使用未消去的子域输运尺度，容差可取 `500ε` 乘以该尺度与环境恢复的条件数 `κ=max(1,ρ/(ρa⁰))`。体积和检查 `|Σ V_i T_i| ≤ 100Nε Σ V_i B_i`，B_i 为局部绝对通量/体积尺度，N 为求和项数；两质量行从零缓冲按相同顺序累加应可直接要求位级相等。记录原始残差、尺度、κ、精度和配置，不能用固定宽松绝对阈值掩盖小流量错误。以上为代数消去与求和误差的建议阈值，实际尚未运行确认；不需要有限差分或长时间积分。
  不要求完整手写 Jacobian 等于完整倾向的精确导数：本 PR 明确采用中央、仅上升气流的近似并省略环境项，完整有限差分比较会混入有意省略的项；这一近似的非线性求解误差不属于上述瞬时输运恒等式。

- Do the PR's own tests already check this invariant (yes/partial/no, naming test files and assertions):
  **1：partial。** `test/prognostic_equations/advection_tests.jl` 验证 Schar/Linear 两端 `max|uₕ³|>0` 且 `max|u³|≤100eps(FT)max|uₕ³|`；没有显式检查 EDMFX 各子域、体积质量和，或人为污染端面缓存后的结构性守恒。
  **2：partial。** `test/prognostic_equations/tracer_mass_consistency_tests.jl` 检查 χ=1 时隐式 `ρq_tot` 与 ρ 倾向相等、显式四种凝结物倾向与构造的质量平流相等，容差 100eps×最大参考倾向，并要求参考非零；SGS 均匀测试检查 ρ、总水及四凝结物倾向 ≤`1e-10*flux_scale`。未覆盖多个 C、Float32、2M 数浓度、非默认重构或完整随机域；使用保持预计算不变的代数状态，不是物理水分组成测试。
  **3：partial。** 同一文件在推进 10 步后要求 SGS 总水倾向非零，且 `parent(T.c.ρ)==parent(T.c.ρq_tot)` 位级相等。未检查体积积分及更广状态域。`test/runtests.jl` 将该文件接入 dynamics/all。以上仅为源码中的断言，未运行测试。

- Do you expect the change to satisfy it (satisfies / may not satisfy / uncertain), with one or two sentences of reasoning (speculation, not a conclusion):
  **1：satisfies（推测）；2：satisfies（在所列分算子及非退化域内的推测）；3：satisfies（推测）。** 顶面直接抵消度量投影，统一零边界散度支持望远镜求和，差值通量在标量对比为零时消失，而总水和湿空气行累加同一个增量；环境恢复/面积正则化附近的精度、不同重构及全时间推进行为仍未验证，不能由这些局部恒等式推出整个湿物理闭合正确。

- Environment: merge commit SHA, language and dependencies needed, estimated setup difficulty (low/medium/high):
  SHA `a8efce5a0040924094a58739db91a0f3162f779c`，Julia，合并版包版本 0.41.1；[Project.toml](https://raw.githubusercontent.com/CliMA/ClimaAtmos.jl/a8efce5a0040924094a58739db91a0f3162f779c/Project.toml) 的兼容项包含 Julia `1.9`、ClimaCore `0.14.54`、ClimaComms `0.6.9`、ClimaTimeSteppers `0.9.0`、CloudMicrophysics `0.36`、Thermodynamics `1.2`、ClimaParams `1.1`，另需项目其余依赖；这些是 compat 约束，不是已解析版本锁。CPU 单进程足够；难度 **medium**，主要在 Julia 依赖/初次编译、配置和缓存一致性，完整测试环境依赖更多。未安装、编译或执行任何库代码，也未读取另一审计轨道、其他本地目录、GitHub 评论、issues 或 forks。

- Formalization (Lean) feasibility: can the invariant be written as a proposition over finite-dimensional reals, and what are the main obstacles:
  **可以。** 固定有限网格，将速度、密度、标量、面积、度量及求积权重写成有限维实向量，D₀ 和 I 为矩阵：由 `g³³>0` 证明投影抵消，由 `VᵀD₀` 只取边界通量及零端面证明质量和为零，由常量保持重构/零输入重构证明 `A(Cρ)=CAρ` 与 SGS 零值，由同一增量证明 `S(ρ−ρq_t)=0`。主要障碍是把 ClimaCore 的几何求积和具体重构对应到矩阵定理、形式化环境诊断的质量加权关系及正则化分支、区分有限差分边界与沉降边界，并另行证明浮点误差界；实数恒等式本身不证明 Julia 缓存调度、近似 Jacobian 或时间积分正确。
