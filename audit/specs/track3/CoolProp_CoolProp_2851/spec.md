## CoolProp/CoolProp#2851 fix(HEOS): correct first_saturation_deriv for mixtures (#2091)

- What changed (1–3 sentences, naming files and formulas):
  `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp` 将混合物饱和切线从纯流体公式 `dT/dp = T(vV−vL)/(hV−hL)` 改为 `−A/B`，其中 A、B 是两相 ln(逸度系数) 对压力、温度偏导之差，泡点用气相组成加权、露点用液相组成加权。混合物仅接受 Q=0 或 Q=1，无参数饱和相重载委托给显式 SatL/SatV 重载；纯流体及 pseudo-pure 分支保留原式。`src/Tests/CoolProp-Tests.cpp` 新增一个泡点差分测试和一个中间 Q 抛异常测试。

- Invariants (mathematical statements; write them as checkable equalities/inequalities, give variables, sign conventions, domain of validity):

  约定所有组成为摩尔分数：总组成 z，液相 x，气相 y，均归一化且分量为正；Q 是气相摩尔品质。T 用 K，p 用 Pa，h 用 J/mol，ρ 用 mol/m³，v=1/ρ；相差统一为气相减液相。固定 z，令 q∈{0,1}，`p_q(T;z)` 为同一连续汽液平衡端点分支：q=0 时 x=z，q=1 时 y=z。限于两相可区分、EOS 光滑、平衡解局部唯一、相关导数分母非零的状态；不含临界点、相分支跳跃、固相或三相平衡。

  1. **Primary — 混合物相平衡切线。** 令 Dq=`first_saturation_deriv(iT,iP)`，Eq=`first_saturation_deriv(iP,iT)`。应有
     \[
     E_q=\frac{d p_q}{dT}\bigg|_z,\qquad D_q E_q=1.
     \]
     更具体地，令 φᵅᵢ(T,p,cᵅ) 为相 α 的逸度系数，偏导时保持该相组成不变，定义
     \[
     A_q=\sum_i w_i^{(q)}[(\partial_p\ln\phi_i^L)_{T,x}-(\partial_p\ln\phi_i^V)_{T,y}],\quad
     B_q=\sum_i w_i^{(q)}[(\partial_T\ln\phi_i^L)_{p,x}-(\partial_T\ln\phi_i^V)_{p,y}],
     \]
     其中 w⁽⁰⁾=y，w⁽¹⁾=x。则 **Aq+Bq Dq=0**，等价于 Dq=−Aq/Bq。A 的单位为 Pa⁻¹，B 为 K⁻¹，D 为 K/Pa。不能要求泡点与露点斜率相等，也不能对所有混合物强加正斜率；正斜率只作为本 PR 的指定普通状态的局部检查。

  2. **Secondary — 固定组成端点物性的总导数。** 对 F∈{hmolar, smolar, rhomolar}，令 Fq(T;z) 为 `update(QT_INPUTS,q,T)` 后主状态的 F（泡点的主体液相或露点的主体气相）。在这些端点，主体组成 z 不随 T 改变，因此
     \[
     \operatorname{satderiv}(F,T)=\frac{dF_q}{dT}
       =(\partial_T F)_{p,z}+(\partial_p F)_{T,z}E_q,
     \qquad
     \operatorname{satderiv}(F,p)=D_q\frac{dF_q}{dT}.
     \]
     这里不声称泡点的微量气相物性或露点的微量液相物性满足这个无组成项公式：这些相的组成会随曲线变化。该不变量覆盖本 PR 新斜率进入既有物性链式求导分支后的实际公共 API 行为。

  3. **Secondary — 真正纯流体的 Clapeyron 一致性及单组分约化。** 对真正纯流体、同温同压的平衡液气相，Δh=hV−hL，Δv=vV−vL，应有
     \[
     D_q\Delta h=T\Delta v,\qquad
     E_q=\frac{\Delta h}{T\Delta v},\qquad D_0=D_1.
     \]
     普通汽液平衡区 Δh>0、Δv>0，故 Dq>0。混合物切线在光滑单组分极限应约化为该式；不强制含零摩尔分数的混合物 API 可用。**不把 pseudo-pure 的保留行为直接解释为严格物理定理**：有温度滑移的拟纯模型未必代表同温同压的真实单组分两相平衡。

  API 范围补充（不计作物理不变量）：混合物成功更新到 0<Q<1 后，本方法应抛异常，不能冒充泡点/露点切线。这不表示固定中间 Q 的曲线在数学上不存在，只表示该 API 未实现那种导数。

- Sources of the invariants (first principles / literature formula / consistency with other implementations in the library / dimensions):

  主不变量可直接从相平衡 `x_i φ_i^L p = y_i φ_i^V p` 推导，不依赖将新实现本身当作正确性标准。对其对数微分：泡点 dx=0，乘以 yᵢ 求和；利用 Σdyᵢ=0 及气相在固定 T,p 下的 Gibbs–Duhem 关系 Σyᵢ d_composition lnφᵢV=0，得到 A₀dp+B₀dT=0。露点 dy=0，改乘 xᵢ 后利用液相的同一关系，得到 A₁dp+B₁dT=0。该推导要求逸度来自热力学一致的混合物模型。

  第二项来自固定组成链式法则；第三项来自 dμ=−s dT+v dp、两相 μ 相等及 Δh=TΔs。单组分时 ∂p lnφ=v/(RT)−1/p，∂T lnφ=−hᵣ/(RT²)，因此 A=−Δv/(RT)、B=Δh/(RT²)，给出同一 Clapeyron 式。

  输入 PR 将混合物式归于 Gernert 学位论文 3.96/3.97；本审计未独立取得论文原文，不把该编号当作已核实的外部证据。已读取的固定版本依据为 [HEOS 实现](https://raw.githubusercontent.com/CoolProp/CoolProp/54a2e32764bfa0ae81877e835c33ccd13ba8f853/src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp)、[公共 AbstractState 接口](https://raw.githubusercontent.com/CoolProp/CoolProp/54a2e32764bfa0ae81877e835c33ccd13ba8f853/include/AbstractState.h) 和 [Python 接口](https://raw.githubusercontent.com/CoolProp/CoolProp/54a2e32764bfa0ae81877e835c33ccd13ba8f853/wrappers/Python/CoolProp/AbstractState.pyx)。未读取其他独立审计轨道、评论、issue 或 fork。

- Common operating domain (range of each input, with the basis):

  首轮取 `HEOS`, `Methane&Ethane`，z=(a,1−a)，a∈[0.6,0.8]，T∈[170,190] K，q 分别取 0、1。基准是 PR 明确测试的 a=0.7、T=180 K、q=0；区间是围绕该点的保守抽样建议，**尚未运行验证，也不是文献保证的整个可用范围**。p 是 QT flash 输出而非独立输入，初筛 10³≤p≤10⁷ Pa；只接纳中心和差分邻点都收敛、保持同一平衡分支的样本。

  额外要求所有相分数为正且归一化，ρL>ρV>0，(ρL−ρV)/ρL>0.05 以避开相密度合并；检查两相逸度一致性，而不单凭 flash 返回就认定平衡正确。此密度门只是数值筛选，不是临界区的严格充分判据。记录失败/未收敛样本，不将其计为不变量通过。中间 Q 异常检查使用原测试点 T=180 K、a=0.7、Q=0.5，可再抽 Q∈[0.1,0.9]。

  纯流体回归用 `HEOS`, `n-Propane`，T∈[280,320] K、q∈{0,1}，以已有回归测试的 300 K 为中心；同样要求平衡与非退化条件。纯组分不是把二元 z 数值置零的调用。

- Numerical check method (which public API of the merged commit, what to compute, what to compare against, suggested tolerance and why; if finite differences are needed, the step-size strategy):

  使用该合并提交的 C++ `AbstractState::factory("HEOS","Methane&Ethane")`、`set_mole_fractions({a,1-a})`、`update(QT_INPUTS,q,T)`、`p()`、`first_saturation_deriv`、`keyed_output`。Python 同名 AbstractState 包装也可实现检查。对中心和邻点使用独立状态或每次重设状态，保存中心解析值后再移动 flash 状态。所有差分保持 **z 和 q** 不变，不锁定微量相的组成。

  **主检查：** 计算 C(h)=[p_q(T+h)−p_q(T−h)]/(2h)，比较 Eq 与 C，以及 Dq·C 与 1；另外检查 Dq·Eq≈1，但这一代数检查不能替代独立 flash 差分。使用 h=0.1、0.05、0.025、0.0125 K 的阶梯，在平滑区寻找二阶收敛；可用 [4C(h/2)−C(h)]/3 外推。建议差分稳定到相对 10⁻⁴ 后要求解析值相对误差≤10⁻³（原 PR 仅要求 5×10⁻³）。若继续减小 h 后误差反增，采用稳定平台；若无平台，标记数值未决，不自动放宽到通过。互逆的代数检查可用 |DqEq−1|≤10⁻¹⁰。接近零斜率、转折或分母退化的样本移出本次非奇异域，另行记录。

  **可选的主检查交叉证据：** 从 `mole_fractions_liquid()`、`mole_fractions_vapor()` 取得 x、y，为两相分别创建单相状态并固定各自组成，用 `specify_phase(iphase_liquid/iphase_gas)` 和 `update(PT_INPUTS,p,T)` 保持相根；以 `fugacity_coefficient(i)` 的对数中心差分构造 A、B。初始温度步长 0.05 K、压力步长 max(1 Pa,10⁻⁴p)，分别反复减半寻找平台；这些偏导不能沿 QT 曲线差分。先要求重建相密度与饱和相密度相对一致至约 10⁻⁶，且 maxᵢ|ln(xᵢφᵢL/(yᵢφᵢV))|≤10⁻⁷；否则判为参考状态构造失败。检查 |A+BDq|/(|A|+|BDq|)≤10⁻³，分母须非零。这条路径交叉验证固定组成的逸度导数与权重，主判据仍为独立 flash 曲线切线。

  **端点物性检查：** 同一组 QT 邻点读取 hmolar、smolar、rhomolar，构造 CF(h)=[Fq(T+h)−Fq(T−h)]/(2h)，比较 `first_saturation_deriv(iF,iT)` 与 CF、`first_saturation_deriv(iF,iP)` 与 CF/C。采用同样的步长平台和 10⁻³ 相对容差。若导数近零，用混合误差界 10⁻³·max(|参考导数|,10⁻⁶SF)，其中温度导数的 SF=max(|F中心|,Funit)/T，Funit 分别为 1 J/mol、1 J/(mol K)、1 mol/m³；压力导数的尺度改为 SF/|C|。这些是待验证的数值容差，不是实验物性精度声明。

  **纯流体检查：** 通过 `saturated_liquid_keyed_output(iHmolar/iDmolar)` 与 `saturated_vapor_keyed_output(...)` 构造 Δh、Δv，比较 |DΔh−TΔv|/max(|DΔh|,|TΔv|)≤10⁻¹⁰；同时独立 QT 差分，建议 10⁻⁴ 相对容差，以防只检查了实现的原式。相同 T 下的 Q=0、1 切线以相对 10⁻⁸ 比较。中间 Q 另检查调用确实抛出异常。所有建议均为后续执行规范；本次没有运行代码或生成数值结果。

- Do the PR's own tests already check this invariant (yes/partial/no, naming test files and assertions):

  **Primary：partial。** 合并版本 `src/Tests/CoolProp-Tests.cpp:4815–4840` 的 `[first_saturation_deriv][2091]` 只检查 Methane/Ethane=(0.7,0.3)、T=180 K、Q=0，步长 0.05 K；断言 D>0，且与反演压力差分的相对误差<5e−3。没有露点、组成扫描、步长收敛、逸度残差或直接 Eq 检查。

  **端点物性：partial。** 同文件 `:2411–2449` 的 `[first_saturation_partial_deriv]` 对纯 n-Propane、300 K、Q=1 检查 P、T、Dmolar、Hmolar、Smolar、Umolar 的十组温度/压力导数，要求 |numerical/analytical−1|<1e−4；未覆盖新混合物分支的物性导数。

  **纯流体 Clapeyron：partial。** 上述旧测试间接检查纯流体的饱和切线，但未直接断言潜热关系，也未比较 Q=0 与 Q=1。`:4843–4849` 对 Q=0.5 的 `first_saturation_deriv(iT,iP)` 使用 `CHECK_THROWS`，完整覆盖该单一状态的新增异常行为。PR 描述声称相关回归通过；本审计只核对测试源码，未复跑。

- Do you expect the change to satisfy it (satisfies / may not satisfy / uncertain), with one or two sentences of reasoning (speculation, not a conclusion):

  **satisfies（推测，限上述非退化域）。** 新的对侧组成权重和负号与相平衡微分及 Gibbs–Duhem 消元一致，端点主体组成固定也支持保留的链式法则，纯流体公式未变。露点与混合物物性导数缺乏本 PR 的直接测试，因此仍需数值验证；对奇异区及 pseudo-pure 的严格物理一致性保持 uncertain。

- Environment: merge commit SHA, language and dependencies needed, estimated setup difficulty (low/medium/high):

  SHA：`54a2e32764bfa0ae81877e835c33ccd13ba8f853`。建议以 C++ 公共 API 为基准；[CMakeLists.txt](https://raw.githubusercontent.com/CoolProp/CoolProp/54a2e32764bfa0ae81877e835c33ccd13ba8f853/CMakeLists.txt) 指定 C++17、CMake≥3.14，依赖通过 CPM 管理。[cmake/dependencies.cmake](https://raw.githubusercontent.com/CoolProp/CoolProp/54a2e32764bfa0ae81877e835c33ccd13ba8f853/cmake/dependencies.cmake) 列出 Eigen 5.0.1、fmt 12.0.0、msgpack、rapidjson、IF97、REFPROP headers、boost headers、multicomplex；运行原生测试还需 Catch2 3.8.0。这里只用 HEOS，不需要 REFPROP 的专有运行库。Python 可作为可选驱动，但须绑定到此 SHA；Python 包构建依赖未完整核实，不能以任意已安装版本代替。预计搭建难度 **medium**：需要精确版本本地构建与依赖解析。本次未安装、编译或运行任何项目代码，未克隆仓库。

- Formalization (Lean) feasibility: can the invariant be written as a proposition over finite-dimensional reals, and what are the main obstacles:

  **可以。** 固定有限组分数 N，以正实数 T,p 及单纯形内部 x,y 为变量，以光滑实函数 lnφᵢᵅ 和局部平衡路径表示 EOS；假设逸度相等、归一化、Gibbs–Duhem、固定主体组成及非零分母，即可形式化推导 A dp+B dT=0、反函数斜率、物性链式法则和单组分约化。主要障碍是从实际多参数 Helmholtz EOS 证明这些热力学一致性前提、用隐函数定理证明局部平衡解存在且唯一、处理临界点及相选择，以及将实数证明连接到浮点计算/flash 求解器误差；证明代数公式并不能证明求解器找到了物理稳定相。中间 Q 抛异常属于程序契约，需独立的程序语义证明，不是实数热力学命题。
