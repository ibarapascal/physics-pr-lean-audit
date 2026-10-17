# 不变量清单（w1 路，独立撰写）

撰写日期：2026-10-17。来源仅限 pilot_list.csv 与 GitHub 公开信息（PR 页面、diff、合并提交上的源码）。
所有「预期是否满足」均为推测，未运行任何代码。

---

## CliMA/Oceananigans.jl#6024 Fix CATKE's bottom mixing length coefficient Cᵇ
- 改了什么：`src/TurbulenceClosures/.../TKEBasedVerticalDiffusivities/catke_mixing_length.jl` 中 `CATKEMixingLength` 的默认值 `Cᵇ` 由 0.28 改为 0.967（docstring 与文档同步）。`Cᵇ` 进入 `stable_length_scale`：`ℓ = min(Cˢ·depth, Cᵇ·height_above_bottom, ℓᴺ)`。
- 不变量：
  - 主（参数代数）：κ_b ≡ Cᵇ·(Cˡᵒu³/CˡᵒD)^{1/4} = 0.40，|κ_b − 0.4| ≤ 1e-3（默认 Cˡᵒu=0.361、CˡᵒD=1.604）。对照：Cˢ 给 κ_s ≈ 0.468。
  - 主（模型行为）：中性层结（N²=0 ⇒ Ri=0，取 Ri≥0 分支，稳定性函数取 Cˡᵒ 值）、无对流、等应力层 Kᵤ∂_z u = u★² 的稳态里，近底面（Cᵇ d < Cˢ·depth，且 ℓ < H）满足 Kᵤ(z) = κ_b·u★·d，d = 底面以上高度；等价地 ∂_z u = u★/(κ_b d)（对数律，κ_b=0.4）。推导：代码里 Kᵤ = Cˡᵒu ℓ √e，ε = e^{3/2}/ℓᴰ，ℓᴰ = ℓ/CˡᵒD（`catke_equation.jl` 的 `ℓ★ = ℓ★ / σᴰ`），局地平衡 Kᵤ S² = ε 与 Kᵤ S = u★² 联立得 e = u★²(CˡᵒD/Cˡᵒu)^{1/2}，Kᵤ = ℓ u★ (Cˡᵒu³/CˡᵒD)^{1/4}。
  - 次（近顶面同式）：Kᵤ = 0.468·u★·depth（本 PR 不改，作对照组，验证推导与代码的一致性）。
- 不变量的来源：文献式（Wagner et al. 2025, JAMES, doi:10.1029/2024MS004522, §3.1.3 式 35）＋第一性原理（壁面相似律、TKE 局地平衡）＋库内一致性（上面推导直接取自合并提交的 `catke_mixing_length.jl`、`catke_equation.jl`、`catke_vertical_diffusivity.jl` 中 `κuᶜᶜᶠ = ℓu·w★`）。
- 常用工况域：u★ ∈ [1e-3, 3e-2] m/s（海洋底边界层典型摩擦速度）；水深 H ∈ [50, 5000] m；底边界层厚 O(10–100 m)，检查区 d ∈ [2Δz, 0.2·H]；Δz ∈ [1, 10] m。
- 数值核对方法：
  1. 参数代数：`using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity`；`c = CATKEVerticalDiffusivity()`；算 `c.mixing_length.Cᵇ * (c.mixing_length.Cˡᵒu^3 / c.turbulent_kinetic_energy_equation.CˡᵒD)^(1/4)`，与 0.4 比，容差 1e-3（默认值只给 3 位有效数字，实际 0.40018）。
  2. 单柱稳态：`RectilinearGrid(size=Nz, z=(-H,0), topology=(Flat,Flat,Bounded))`，`HydrostaticFreeSurfaceModel(; grid, closure=CATKEVerticalDiffusivity(), buoyancy=BuoyancyTracer(), tracers=(:b,:e))`，b≡0；u 的顶、底动量通量均设为 −u★²（同号，使柱内应力处处等于 u★²），f=0。积分至 Kᵤ 稳定（相对变化 < 1e-4 / 惯性时间），取 `model.closure_fields.κu` 在面上的值，除以 u★·d_face，在 Cᵇ d < Cˢ(H−d) 的底层区间内应 ≈ 0.40，顶层区间应 ≈ 0.468。容差 5%（离散：√e 由中心向面插值，边界 TKE 通量 Cᵂϵ、Cᵂu★ 只影响紧邻边界 1–2 格，检查时排除这些格）。用 Δz 减半两次看是否收敛到 0.4，比单点容差更稳。
- PR 自带测试是否已检查这条不变量：否。diff 只改了默认值、docstring 与 `docs/src/models/turbulence_closures.md` 的打印样例，无新增断言。
- 你预期改动是否满足：满足（推测）。按合并提交的代码推导出的 κ 式与 PR 引用的式 35 一致，0.967 代入得 0.4002。附带观察：旧值 0.28 恰等于 0.4·(Cˡᵒu/CˡᵒD)^{1/4}（即 Cˡᵒu 与 CˡᵒD 位置颠倒且未立方），和 PR 描述的「Cˡᵒu 未立方」不完全一致；旧值对应的实际 κ_b 只有约 0.116。
- 环境：合并提交 24de7026893af38df065920fbbf80f05024fd1b9；Julia ≥1.10 ＋ Oceananigans.jl（CPU 即可）；难度 低（代数检查）／中（单柱稳态，需 Julia 预编译与积分到稳态）。
- 形式化（Lean）可行性：高。(a) 数值命题可写成 rpow 自由的形式：0.3995⁴·CˡᵒD ≤ Cᵇ⁴·Cˡᵒu³ ≤ 0.4005⁴·CˡᵒD，`norm_num` 直接可证。(b) 推导命题：对 ℓ,u★,e>0，若 K = Cu·ℓ·√e、K·S = u★²、K·S² = CD·e^{3/2}/ℓ，则 K⁴ = ℓ⁴u★⁴·Cu³/CD——有限维实数代数，主要障碍只是 √ 与 e^{3/2} 的改写。模型行为（离散稳态）不宜形式化。

---

## pybamm-team/PyBaMM#5888 fix: convert Ai2020 electrolyte diffusivity from cm2/s to m2/s
- 改了什么：`packages/pybamm/src/pybamm/input/parameters/lithium_ion/Ai2020.py` 的 `electrolyte_diffusivity_Ai2020(c_e, T)` 由 `10**(-4.43 - 54/(T - 229 - 5e-3 c_e) - 0.22e-3 c_e)` 改为前面乘 `1e-4`（Valøen–Reimers 2005 式 14 给 cm²/s，参数名要求 m²/s）；docstring 改正。
- 不变量：
  - 主（单位换算一致）：对域内任意 (c_e, T)，D_Ai2020(c_e,T) [m²/s] = 1e-4 · 10^{−4.43 − 54/(T − 229 − 5×10⁻³c_e) − 2.2×10⁻⁴c_e}，c_e 以 mol/m³、T 以 K 计；等价地 log10(D) + 4 = −4.43 − 54/(T−229−0.005c_e) − 0.00022c_e。
  - 主（库内一致性）：D_Ai2020(c_e,T) = D_Xu2019(c_e,T)（`Xu2019.electrolyte_diffusivity_Valoen2005`，同一拟合式），相对误差 ≤ 1e-12。
  - 次（物理量级与单调性）：在域内 1e-11 ≤ D ≤ 2e-9 m²/s（LiPF6/碳酸酯电解液的公认量级）；∂D/∂T > 0，∂D/∂c_e < 0（只要 T − 229 − 0.005c_e > 0，两项导数符号可直接由式子读出）。旧实现在 1000 mol/m³、298.15 K 给 3.2e-6，违反此界。
- 不变量的来源：文献式（Valøen & Reimers, J. Electrochem. Soc. 152(5) A882 (2005), eq. 14，cm²/s）＋量纲（参数名 `[m2.s-1]`，1 cm² = 1e-4 m²）＋库内其它实现一致性（Xu2019）。
- 常用工况域：c_e ∈ [100, 3000] mol/m³（锂离子电池放电中电解液浓度常见范围，PR 表中 1C DFN 末态 668–1683）；T ∈ [253, 333] K（−20…60 °C，电池工作温区；Valøen 拟合的测量温区约 −10…50 °C〔推测，未核原文〕）。须避开奇点 T = 229 + 0.005c_e（c_e=3000 时 244 K）。
- 数值核对方法：`pv = pybamm.ParameterValues("Ai2020")`；`f = pv["Electrolyte diffusivity [m2.s-1]"]`；对 (c_e,T) 网格（如 c_e 取 100…3000 步 100，T 取 253…333 步 5）算 `pv.evaluate(f(pybamm.Scalar(c), pybamm.Scalar(T)))`，与 numpy 直接按文献式（含 1e-4）计算比，rtol 1e-12（同一浮点式，只差运算顺序）；同网格与 `pybamm.ParameterValues("Xu2019")` 的同名参数比，rtol 1e-12；再查量级界与有限差分符号（步长 ΔT=0.01 K、Δc=0.1 mol/m³，中心差分，只看符号）。可选：1C DFN 放电（Ai2020）末态 c_e 应出现明显梯度（max−min ≫ 1 mol/m³），作为宏观旁证。
- PR 自带测试是否已检查这条不变量：是。`packages/pybamm/tests/unit/test_parameters/test_parameter_sets/test_Ai2020.py` 新增 `test_electrolyte_diffusivity`：Ai2020 与 Xu2019 在 (1000,298.15)、(500,318.15)、(1500,278.15) 三点 `assert_allclose(rtol=1e-12)`；`test_functions` 把容差由 `abs=1e-4` 改为 `rel=1e-3` 并检查 1000 mol/m³、298.15 K 处 3.2227e-10。量级界与单调性未测。
- 你预期改动是否满足：满足（推测）。手算 10^{−4.43−54/64.15−0.22}·1e-4 = 3.22e-10 m²/s，落在典型量级；与 Xu2019 是同式。
- 环境：合并提交 343ece2062dbb39789dbfa977e0f2daf96e9e370；Python ≥3.10 ＋ pybamm（从该提交 `packages/pybamm` 安装）＋ numpy；难度 低。
- 形式化（Lean）可行性：高。主不变量是实数上的恒等式（两边同一表达式乘常数），量级界可写成对 c_e∈[100,3000]、T∈[253,333] 的区间命题：指数项 −4.43−54/(T−229−0.005c)−0.00022c 的上下界由单调性得出，再用 10^x 单调；障碍只是 Real.rpow 底 10 的数值包络，需要 `interval_cases`／手写有理上下界。

---

## pybamm-team/PyBaMM#5889 fix: look up Chayambuka2022 NVPF diffusivity at the solid concentration
- 改了什么：`packages/pybamm/src/pybamm/input/parameters/sodium_ion/Chayambuka2022.py` 的 `NVPF_diffusivity_Chayambuka2022(sto, T)` 中插值自变量 `sto * c_max` 的 `c_max` 由 `"Initial concentration in electrolyte [mol.m-3]"`（1000）改为 `"Maximum concentration in positive electrode [mol.m-3]"`（15320）。`D_p.csv` 的横轴是正极颗粒浓度 131.6–15197.4 mol/m³。
- 不变量：
  - 主（自变量量纲与覆盖）：D_p(sto,T) = I_Dp(sto · c_{p,max})，其中 I_Dp 是 `D_p.csv` 的插值，c_{p,max} = 15320 mol/m³；当 sto 取遍 [0,1] 时，查表自变量 sto·c_{p,max} 覆盖 CSV 横轴 [131.6, 15197.4] 的全部（即 sto ∈ [0.00859, 0.99199] 落在表内）。
  - 主（与负极同构）：正、负极颗粒扩散率都按「固相浓度 = sto × 本电极 c_max」查表，即 D_n(sto) = I_Dn(sto·c_{n,max})、D_p(sto) = I_Dp(sto·c_{p,max})；对任意固相浓度 c_s ∈ 表内，D_p(c_s/c_{p,max}) = I_Dp(c_s)（与 c_{p,max} 的数值无关，换 c_{p,max} 只改变 sto 到 c_s 的映射，不改变 c_s 到 D 的映射）。
  - 次（量级）：sto ∈ [0.01,0.99] 内 D_p ∈ [1.9e-17, 2.6e-15] m²/s（CSV 列的值域），且在 sto ≥ 0.3 时 D_p ≤ 1e-16 m²/s。
- 不变量的来源：库内一致性（`HC_diffusivity_Chayambuka2022`、`*_exchange_current_density` 都以固相浓度为自变量）＋数据表的轴定义（`D_p.csv` 表头 "Positive particle concentration [mol.m-3]"）＋量纲（sto 无量纲，乘以电极固相最大浓度才得到与表轴同量纲同物理含义的量）。文献：Chayambuka et al., Electrochim. Acta 404 (2022) 139764。
- 常用工况域：sto ∈ [0.01, 0.99]（钠离子 NVPF 正极充放电区间）；T ∈ [273, 323] K（该函数不依赖 T，T 只是签名参数）。
- 数值核对方法：`pv = pybamm.ParameterValues("Chayambuka2022")`；`f = pv["Positive particle diffusivity [m2.s-1]"]`；对 sto 网格 0.01…0.99（步 0.01）算 `pv.evaluate(f(pybamm.Scalar(sto), pybamm.Scalar(298.15)))`，与 `numpy.interp(sto*15320, x_csv, y_csv)` 比（`x_csv,y_csv` 读合并提交的 `data/D_p.csv`）；PyBaMM 1D 插值默认 linear，rtol 1e-10。若默认插值器不是线性则改为与 `scipy.interpolate` 同种插值比，或只在 CSV 节点处比（节点处任何插值都应精确），rtol 1e-10。再把 `pv.update({"Maximum concentration in positive electrode [mol.m-3]": c})` 改成别的 c，检查 f(c_s/c) 在固定 c_s 下不变（rtol 1e-10）。宏观旁证：`pybamm.sodium_ion.BasicDFN()` 默认参数 2C 放电至 2.0 V 容量约 1.18 mAh（PR 自报，非物理不变量，只作回归参考）。
- PR 自带测试是否已检查这条不变量：部分。`test_Chayambuka2022.py` 新增 `test_positive_particle_diffusivity` 在 sto=0.3、0.9 检查定值（rel=1e-4），`test_functions` 在 sto=0.5 检查 8.3317e-17。只是三个点的定值回归，没有检查与 CSV 的整段一致、覆盖范围或与 c_max 的同构关系。另：同一测试把电解液扩散率期望值由 2.5061e-10 改为 2.5061e-11（代码未改，旧值靠 `abs=1e-4` 蒙混过去），说明该表的电解液扩散率实际是 2.5e-11 m²/s，偏低一个量级是否符合原文未核。
- 你预期改动是否满足：满足（推测）。修改后自变量同量纲、与负极同构，sto∈[0,1] 正好覆盖表轴。
- 环境：合并提交 eda03ed6d1003b186cc2d70ad8b3240a7b6350f1；Python ≥3.10 ＋ pybamm（`packages/pybamm`）＋ numpy；难度 低。
- 形式化（Lean）可行性：中高。覆盖命题是两个有理不等式（131.58/15320 ≥ 0 且 15197.37/15320 ≤ 1），平凡；「与 c_max 无关」是复合函数恒等式 I((c_s/c)·c) = I(c_s)（c≠0），平凡。真正有物理内容的部分（表的物理正确性）无法在 Lean 中表达，只能对数据点做有限检查。
---

## CoolProp/CoolProp#2851 fix(HEOS): correct first_saturation_deriv for mixtures (#2091)
- 改了什么：`src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp` 的 `calc_first_saturation_deriv`。纯/准纯流体仍用 Clapeyron 式 dT/dp = T(1/ρ_V − 1/ρ_L)/(h_V − h_L)；混合物改为 dT/dp = −Σ_i w_i Δ(∂lnφ_i/∂p)_{T,n} / Σ_i w_i Δ(∂lnφ_i/∂T)_{p,n}，Δ = 液相 − 气相，Q=0（泡点线）时 w = 气相组成 y，Q=1（露点线）时 w = 液相组成 x，其它 Q 抛异常（Gernert 学位论文式 3.96/3.97）。无参重载改为委托有参版本。
- 不变量：
  - 主 A（沿泡点线/露点线的切线）：总组成 z 固定，记 p_b(T) 为 Q=0 的饱和压力、p_d(T) 为 Q=1 的饱和压力。则 `first_saturation_deriv(iT,iP)` 在 Q=0 处 = (dp_b/dT)⁻¹，在 Q=1 处 = (dp_d/dT)⁻¹。理论依据：泡点线上 Σ_i z_i K_i = 1，K_i = φ_i^L/φ_i^V；全微分后气相组成变化项按 Gibbs–Duhem（Σ_i y_i d lnφ_i^V|_{T,p} = 0）消去，得到上式且权重恰为 y_i = z_i K_i；露点线同理权重为 x_i。
  - 主 B（纯流体极限）：二元混合物 z = (1−ε, ε)，ε→0 时 Q=0 与 Q=1 两个值都收敛到纯组分 1 的 Clapeyron 值 T Δv/Δh，|差| = O(ε)。（由 (∂lnφ_i/∂p)_T = v̄_i/(RT) − 1/p、(∂lnφ_i/∂T)_p = −h̄_i^res/(RT²) 代入，权重退化为 δ_{i1}。）
  - 次 C（链式导数）：`first_saturation_deriv(iP,iT)` = 1/`first_saturation_deriv(iT,iP)`；且对 Q=0（体相为组成固定 = z 的液相）`first_saturation_deriv(iHmolar,iT)` = d h_L/dT 沿泡点线（Q=1 时为气相 h_V 沿露点线）。
- 不变量的来源：第一性原理（相平衡条件 f_i^L = f_i^V 沿饱和线的全微分＋Gibbs–Duhem）；文献式（Gernert 学位论文 3.96/3.97，PR 引用）；库内一致性（纯流体 Clapeyron 分支）。
- 常用工况域：Methane&Ethane（z_CH4 ∈ [0.1, 0.9]，T ∈ [120, 240] K 且距临界/cricondentherm ≥ 5 K）；另取一组制冷剂二元 R32&R125 或 R32&R1234yf（z ∈ [0.2,0.8]，T ∈ [230, 330] K）。近临界点 Σ w_i Δ(∂lnφ/∂T) → 0、两相线合一，排除。依据：PR 测试工况（0.7/0.3 甲烷乙烷、180 K）与 CoolProp 常用制冷剂混合物。
- 数值核对方法：CoolProp Python 低层接口 `AS = CoolProp.AbstractState("HEOS","Methane&Ethane")`，`AS.set_mole_fractions([z,1-z])`，`AS.update(CoolProp.QT_INPUTS, Q, T)`，`AS.first_saturation_deriv(CoolProp.iT, CoolProp.iP)`。A：另建状态在 T±h 做 QT 闪蒸取 p，中心差分 2h/(p₊−p₋)；h 取 0.2、0.1、0.05 K 做 Richardson 外推（闪蒸收敛误差约 1e-8 相对，h 太小会被它淹没），对 Q=0 和 Q=1 都做，容差相对 1e-4（PR 自测用 5e-3 是偏松的）。B：ε 取 1e-2、1e-3、1e-4，与 `AbstractState("HEOS","Methane")` 在同 T、Q=0 的 `first_saturation_deriv(iT,iP)` 比，检查差值随 ε 线性下降。C：同一 QT 状态求 `first_saturation_deriv(iHmolar,iT)`，与沿 Q=0 线 h(T±h) 的中心差分比，相对 1e-4。另核 Q=0.5 时抛异常。
- PR 自带测试是否已检查这条不变量：部分。`src/Tests/CoolProp-Tests.cpp` 新增 `[2091]` 两例：Methane&Ethane 0.7/0.3、180 K、Q=0 与 dT=0.05 K 中心差分比，相对 5e-3，且 dT/dP>0；Q=0.5 抛异常。Q=1（露点分支）、纯流体极限、链式导数都没测。
- 你预期改动是否满足：满足（推测）。代码式子与上面由 Σ z_i K_i = 1 推出的结果在权重与符号上逐项一致。残余风险（不确定）：判断用 `this->_Q == 0` 精确比较，若别的输入对（如 PQ）闪蒸后 _Q 存为 1e-16 级非零值会误抛；有参重载用的是 `this` 的 Q 而非传入 SatL/SatV 的状态。
- 环境：合并提交 54a2e32764bfa0ae81877e835c33ccd13ba8f853；C++17 ＋ CMake 构建 CoolProp Python 扩展（或该提交之后的发布版 wheel，未核是否已发布）；难度 中（源码构建约十几分钟，内存需求不大）。
- 形式化（Lean）可行性：部分可行。最后一步代数（Σ w_i a_i dp + Σ w_i b_i dT = 0 ⇒ dT/dp = −Σwa/Σwb）与纯流体极限（权重 δ_{i1} 时化为 T Δv/Δh，给定 v̄、h̄^res 的偏导恒等式作为假设）是有限维实数命题，容易。核心的 Gibbs–Duhem 消去需要把 lnφ_i 作为 (T,p,n) 的齐次函数并用欧拉定理，可以在 Mathlib 里以「给定满足 Σ_i n_i ∂_j lnφ_i = 0 的导数矩阵」为假设来写，障碍是不对具体 Helmholtz 能 EOS 建模。

---

## CoolProp/CoolProp#3339 fix(transport): rhosr-CS viscosity ignored x_crossover from the fluid file
- 改了什么：`src/Backends/Helmholtz/TransportRoutines.cpp` 的 `viscosity_rhosr`：切换函数 ψ_liq = 1/(1+exp(−100(x − 2))) 中写死的 2 改为读流体文件的 `data.x_crossover`；陡度 100 改成具名常量。x = ρ R (τ α^r_τ − α^r)/(ρ s_r)_crit = −ρ s^r/(ρ s^r)_crit。
- 不变量：
  - 主（对已发布数据为恒等变换）：8 个 rhosr-CS 流体（R1233zd(E)、R1234yf、R1234ze(E)、R124、R152A、R22、R245fa、R32）的 `x_crossover` 均为 2，因此对任意 (T,ρ)，η_merge(T,ρ) 与父提交 afcc79091533dfbd258cd59c065c0765747fd230 的 η_parent(T,ρ) 按位相等。
  - 主（参数语义）：η(T,ρ; x_c) = η_dilute(T)·[1 + C(exp(ψ f_L(x) + (1−ψ) f_V(x)) − 1)]，ψ = 1/(1+exp(−100(x − x_c)))，f_L、f_V 为流体文件 `c_liq`、`c_vap` 的三次多项式；η_dilute 与 x_c 无关。故两克隆（只改 x_c=a,b）之比 η_a/η_b = [1 + C(e^{g_a}−1)]/[1 + C(e^{g_b}−1)]，g_k = ψ(x;k) f_L + (1−ψ(x;k)) f_V，完全由流体文件参数与可算的 x 决定。
  - 次（零密度极限，与本 PR 无直接关系但同一函数）：ρ→0 时 x→0、ψ→0（x_c=2 时 ψ≈e^{−200}），η → η_dilute·[1 + C(exp(c_vap[0]) − 1)]；若关联式设计为稀薄气体极限一致，应有 c_vap[0] = 0、η → η_dilute。
- 不变量的来源：库内一致性（解析器读入的字段应被使用；式子取自合并提交源码）＋ IEEE-754（JSON 整数 2 转 double 精确等于 2.0）。物理背景为残余熵标度（Bell 等人的 ρ·s_r 对应态黏度模型〔文献出处推测，未核〕）。
- 常用工况域：T ∈ [230, 420] K，ρ 从 1 mol/m³ 到饱和液密度（约 1.0–1.3e4 mol/m³），单相（这 8 种制冷剂的常用工作区）；x 由此覆盖约 [0, 10]。
- 数值核对方法：Python CoolProp（合并提交构建）`PropsSI("V","T",T,"Dmolar",ρ,fluid)`。主1：同一 (T,ρ) 网格在父提交与合并提交两版各算一次，比较 `float.hex` 是否相同（需要两次构建）。主2：仿 PR 测试用 `CoolProp.add_fluids_as_JSON("HEOS", json)` 注册 x_c=a、b 的克隆；x 由 `AbstractState.smolar_residual()` 与 `rhomolar()` 及流体文件 `rhosr_critical` 算出（x = −ρ s^r/rhosr_critical），用流体文件的 C、c_liq、c_vap 按上式算理论比值，与 `PropsSI` 比值比，rtol 1e-12。网格取 x ∈ {x_c ± 0.01, ±0.05, ±0.5} 覆盖过渡区内外。次：ρ = 1e-3 mol/m³ 时读 c_vap[0]，看 η/η_dilute 是否 → 1（η_dilute 未单独暴露，用 ρ 递减外推）。
- PR 自带测试是否已检查这条不变量：部分。`src/Tests/CoolProp-Tests.cpp` 新增「rhosr-CS viscosity honours x_crossover」：克隆 R1233zd(E)，x_c=2 与原版相对差 <1e-14，x_c=12 在 T=300 K、ρ=1e4 mol/m³（x≈9）处差 >1%。没有检查比值的定量公式，也没对 8 种流体做按位回归（PR 文字称全测试集结果不变）。
- 你预期改动是否满足：满足（推测）。改动对已发布数据是逐位恒等，参数语义也直接按式子实现。物理上可检验的内容很少，这份更接近接线修复。
- 环境：合并提交 5418cec58c52183e3a010b2894e08b818ed4cc81（父 afcc79091533dfbd258cd59c065c0765747fd230）；C++17＋CMake 构建 CoolProp Python 扩展，按位回归需构建两版；难度 中。
- 形式化（Lean）可行性：实数层面平凡（x_c = 2 时 x − x_c = x − 2）；按位相等是浮点命题，超出有限维实数，需要 Flocq 一类的浮点形式化，不建议。比值公式是实数恒等式，可证但没有物理内容。
---

## e2nIEE/pandapower#3016 Neutral earthing impedance (rn_ohm) and zigzag earthing transformers in the zero-sequence short circuit
- 改了什么：`pandapower/create/trafo_create.py` 增加 `rn_ohm` 参数；`pandapower/pd2ppc_zero.py::_add_trafo_sc_impedance_zero` 中 Z_N = rn_ohm + j·xn_ohm，(1) 对 dyn/ynd/yyn/ynyn/yny/yzn 非 PSU 变压器把 3Z_N（以接地绕组侧母线电压换算成标幺）串入零序支路，且在 IEC 修正系数 K_T（`z0_k *= kt`）之后加，不被 K_T 缩放；PSU 分支把原 3j·xn 改为 3(rn + j·xn)；(2) 新增 ZNyn/ZNd/ZNy/ZN（代码里 `zd`）零序模型：高压侧对地并联导纳 1/(si0_hv_partial·Z0_leak + 3Z_N)，高低压零序解耦，ZNyn 低压侧对地 1/(Z0_mag + (1−si0_hv_partial)·Z0_leak)，ZNd/ZNy 低压侧无零序通路。
- 不变量：
  - 主 A（中性点阻抗按 3Z_N 串联且不被 K_T 缩放）：设母线 b 的零序戴维南阻抗 Z0_th(b; Z_N)（Ω，`res_bus_sc.rk0_ohm + j·xk0_ohm`）。对 Dyn 变压器的低压母线（低压侧唯一零序通路就是该变压器）：Z0_th(b; Z_N) − Z0_th(b; 0) = 3Z_N（Ω，在 b 的电压等级上），对 case="max" 与 case="min" 都精确成立（两者 K_T 不同，但差值不应变）；正序 Z1_th(b) 与 Z_N 无关。一般网络（如 YNd 高压侧另有电网零序并联）：Z0_th(b; Z_N) = (Y_rest + 1/(K_T·Z0_T + 3Z_N))⁻¹，即 1/(1/Z0_th(b;Z_N) − Y_rest) − 3Z_N 对不同 Z_N 取同一常数 K_T·Z0_T（Y_rest 由断开该变压器的同网络算出）。
  - 主 B（IEC 60909 单相短路电流）：I″_k1 = √3·c·U_n / |2Z1_th + Z0_th|（Z2 = Z1）；故 |2Z1+Z0| 随 R_N≥0 单调不减，且在 Z_N 为纯电阻、Z0_th 实部≥0 时 I″_k1 随 R_N 严格单调递减。
  - 次 C（曲折变压器）：接地变压器 ZN* 接在无其它零序通路的母线 b_hv（三角形侧馈电）上时，Z0_th(b_hv) = si0_hv_partial·K_T·Z0_leak + 3Z_N（Ω 换算后），与次级绕组是 yn 还是 d 无关（ZNyn 与 ZNd 的 b_hv 结果按位相同）；ZNyn 的低压母线 Z0_th 与高压侧网络（电网零序参数、Z_N）无关；ZNd 低压母线 I″_k1 = 0。
- 不变量的来源：IEC 60909-0（单相短路电流式；中性点接地阻抗在零序网络中以 3Z_N 出现；K_T 只作用于变压器自身阻抗）＋第一性原理（对称分量法，中性线电流 = 3I0）＋库内一致性（PSU 分支早已用 3·xn、非缩放的写法）。曲折变压器零序等值电路依据 PR 所引 PowerFactory 技术参考与 J&P Transformer Book，si0 分摊方式是建模约定而非第一性原理。
- 常用工况域：中压 10–35 kV、低压 0.4 kV、高压 110–132 kV；变压器 0.1–100 MVA，vk 4–14 %，vkr 0.3–1.5 %；R_N ∈ [0, 1000] Ω（中压接地电阻常把接地故障电流限到 0.1–2 kA，对应 20 kV 时约 6–120 Ω）；X_N ∈ [0, 100] Ω；si0_hv_partial ∈ [0.5, 0.95]。依据：PR 测试与验证所用 132/33/0.8 kV 光伏汇集网及 IEC 60909 典型值。
- 数值核对方法：`pp.create_empty_network`、`create_ext_grid`（给 s_sc、rx、x0x、r0x0）、`create_transformer_from_parameters(..., vector_group=..., vk0_percent, vkr0_percent, mag0_percent, mag0_rx, si0_hv_partial, rn_ohm, xn_ohm)`，`pandapower.shortcircuit.calc_sc(net, fault="1ph", case=..., inverse_y=True/False)`，读 `net.res_bus_sc[["rk0_ohm","xk0_ohm","rk_ohm","xk_ohm","ikss_ka"]]`。A：Dyn 网络对 Z_N ∈ {0, 0.01, 0.1+0.05j, 10} Ω × case ∈ {max,min} 检查 ΔZ0 = 3Z_N，容差绝对 1e-9 Ω 或相对 1e-9（纯线性代数，误差只来自矩阵求逆）；YNd 网络按上面并联式反推常数，对不同 Z_N 一致到相对 1e-9。B：用同一结果与 √3·c·U_n/|2Z1+Z0| 比（c 取 `net.bus.c_max`/IEC 表，或从 `ikss_ka` 与 Z 反推 c 后检查各 Z_N 下 c 不变），相对 1e-9；R_N 扫描检查单调。C：三角形馈电中压母线＋ZNyn/ZNd，比较两者 b_hv 的 Z0 与 si0·K_T·Z0_leak+3Z_N（K_T 需按 IEC 式自己算，或用 case 切换验证 3Z_N 部分不随 K_T 变），再改电网零序参数检查 ZNyn 低压 Z0 不变。附加（推测风险）：`mode="pf_3ph"` 下 Dyn/YNyn/Yzn 也会走到新加的 3Z_N 代码（该段在 `if mode=="sc"` 外），而 pf_3ph 的变压器标幺用 3·sn_mva 基准、zn_pu 用 sn_mva 基准，建议用 `runpp_3ph` 在单相负荷下检查「Dyn＋Z_N」与「Dyn＋在低压中性点外接等值 Z_N」的中性点位移是否一致。
- PR 自带测试是否已检查这条不变量：部分。`pandapower/test/shortcircuit/test_1ph.py`：`test_trafo_neutral_earthing_impedance` 只在 Dyn、case="max" 检查 ΔR0 = 3·0.01、ΔX0 = 3·0.02（atol 1e-6）、正序不变、ikss 下降——没有测 case="min"（所以「不被 K_T 缩放」只被单点覆盖），没有测并联情形，也没有对 IEC 电流式做定量核对。`test_zigzag_earthing_transformer` 只检查 ZNyn 与 YNd 电流比在 0.5–2、ZNyn 与 ZNd 高压电流相同、ZNd 低压电流 0、R_N 使电流下降，未核 Z0 的定量式。
- 你预期改动是否满足：sc 模式下 A、B 满足（推测，代码在 kt 之后加 zn_pu，基准用接地侧母线电压）；C 按代码本身满足，但其物理正确性依赖 si0 分摊约定，不确定。pf_3ph 模式可能不满足（推测）：新代码对 pf_3ph 也生效且基准因子可能差 3 倍，且 YNyn 两侧都接地时只按高压侧换算一次 3Z_N，语义不明确。
- 环境：合并提交 0c6a1f0a65273af91dd1d5b854498cd89d457939；Python ≥3.10 ＋ pandapower（该提交源码 pip install）＋ numpy/scipy/pandas；难度 低。
- 形式化（Lean）可行性：A 的串联与并联式、B 的单调性都是复数（可拆成实部虚部的有限维实数）代数命题，可写可证；例如「R_N ≥ 0、Re(2Z1+Z0₀) ≥ 0 ⇒ |2Z1+Z0₀+3R_N| 关于 R_N 单调」只需展开模平方。障碍是要把 pandapower 的导纳矩阵求解抽象成「戴维南阻抗 = 网络其余部分 ∥ 变压器支路」这个假设，不能形式化库代码本身。
---

## idaholab/moose#34020 Pass the nodal slip rate to dynamic mortar function_friction
- 改了什么：`modules/contact/src/constraints/ComputeDynamicFrictionalForceLMMechanicalContact.C`：传给速度相关摩擦函数的滑移速率由「Σ_qp ψ_j(qp)·(−v_rel)·t_j」（未乘 JxW、未归一化）改为 s_j = ∫_Γ ψ_j (−v_rel)·t_j dA / ∫_Γ ψ_j dA（`timestepSetup` 中 weighted tangential velocity 除以 `_dof_to_weighted_gap[dof].second`，Wohlmuth 2011 式 3.10a 的 mortar 投影节点系数）；为此 `ComputeDynamicWeightedGapLMMechanicalContact.C` 无条件累加 ∫ψ_j dA。无历史的节点滑移速率取 0；常数 `mu` 加 `mu > 0` 范围检查。摩擦函数以 f(t=0, x=p_old, y=|s_j|, z=0) 调用。
- 不变量：
  - 主 A（常值场精确再现）：若接触面上相对切向速度处处等于常向量 V（刚体均匀滑动、平直界面、节点切向 t_j 常值），则每个活动接触节点 s_j = −V·t_j（2D 一分量，3D 两分量），|s_j| = |V|，且与网格尺寸、单元阶次、积分阶次、是否用 dual 基无关。原因：∫ψ_j c dA / ∫ψ_j dA = c 恒成立。
  - 主 B（等效性）：在 A 的工况下，从第 2 个时间步起（摩擦系数用上一步的 s_j），`function_friction = g(y)`（只依赖 y=|s|）的结果与常数 `mu = g(|V|)` 的结果相同（位移、法向/切向 LM 一致到非线性求解容差）。
  - 次 C（线性场，dual 基）：若用 dual Lagrange 乘子（`use_dual = true`），相对切向速度为界面坐标的线性函数时，s_j 等于该线性场在节点 j 的值（双正交性 ∫ψ_j φ_k = δ_jk ∫φ_k）；标准基时为集中质量式近似，误差 O(h²)。
- 不变量的来源：数值方法文献（Wohlmuth 2011 Acta Numerica「Variationally consistent discretization schemes and numerical algorithms for contact problems」式 3.10a 的 mortar 投影）＋库内一致性（PR 说准静态约束 #34016 已用同一量）＋量纲（s_j 必须是速度 m/s，旧量的量纲是速度乘以积分点个数，且随积分阶次改变）。
- 常用工况域：滑移速度 |V| ∈ [1e-4, 10] m/s（MOOSE 动力学接触测试与工程摩擦试验范围）；接触压力 p ∈ [1e3, 1e9] Pa；网格 h 使界面上 4–64 个单元；时间步 Δt 使每步滑移 < 0.1 h；Newmark β=0.25、γ=0.5。摩擦律取单调函数，如 μ = 0.2 + 0.1·tanh(y/0.1)，保证 μ>0 且对速度敏感。
- 数值核对方法：MOOSE contact 模块（合并提交构建的 `contact-opt` 或 `combined-opt`）。仿 `modules/contact/test/tests/mortar_dynamics/frictional-mortar-3d-dynamics-light-function.i` 改写一个 2D 平面问题：上块被位移 Dirichlet BC 驱动做 u_x = V·t 的匀速切向运动并施加法向压紧，下块固定；加 `AuxVariable` 输出摩擦 LM 与法向 LM（或用 `NodalValueSampler`/`VectorPostprocessor` 输出接触节点上的值）。A：在约束类里没有直接输出 s_j 的公开量，故用 μ 反推：μ_eff = |λ_t|/λ_n（滑移状态下 Coulomb 律取等号），在第 ≥2 步与 g(|V|) 比；同时把 `function_friction` 设为 g(y) = y（纯线性）可直接读出 |s_j| = μ_eff。容差：相对 1e-6（受非线性容差 nl_rel_tol 和半光滑牛顿收敛控制，需把 nl_rel_tol 设到 1e-10）。网格加密 2 次，μ_eff 不应变化（旧实现会随积分点数变化）。B：同网格分别跑 `function_friction=g` 与 `mu=g(V)`，exodiff 相对 1e-6。C：上块改为绕法向的匀速转动（2D 中让切向速度沿界面线性分布），`use_dual=true`，比较 μ_eff 推出的 |s_j| 与解析线性场。
- PR 自带测试是否已检查这条不变量：否（只有回归）。唯一相关测试 `frictional-mortar-3d-dynamics-light-function`（单元素弹跳块）被重新 gold，PR 文字称摩擦系数由约 −0.48 变为 0.48；exodiff 只比对新旧输出，不检查 s_j 与真实滑移速度的关系或网格无关性。
- 你预期改动是否满足：A、C 满足（推测，按定义恒成立）；B 大体满足，但「第 1 步与刚进入界面的节点取 s=0」意味着那些节点上 μ = g(0) 而非 g(|V|)，核对时要排除。附带风险（不确定）：归一化量 ∫ψ_j dA 来自另一个约束对象的 `_dof_to_weighted_gap` 在 `timestepSetup` 时刻的值，若该映射在 timestep 之间被清空或只含活动集，可能出现查不到 dof 的情况（代码用 `libmesh_map_find`，查不到会报错而不是给出错误值）。
- 环境：合并提交 f90f508580f7611f7d75ee458cd177affc7288d1；C++、PETSc、libMesh、MOOSE 框架＋contact 模块（conda 的 moose-dev 环境）；难度 高（完整构建需数 GB 内存与较长编译时间，不宜在 8GB 本机做）。
- 形式化（Lean）可行性：A 可写成有限维命题：给定权 w_q ≥ 0、测试函数值 ψ_j(q)、Σ_q w_q ψ_j(q) ≠ 0，则 (Σ_q w_q ψ_j(q) c)/(Σ_q w_q ψ_j(q)) = c——平凡。C 需要双正交性作假设，也是有限维线性代数。B 涉及非线性动力学求解，不可形式化。
---

## idaholab/moose#33307 Use log10 in relevant SCM closures and update transient regime calculation
- 改了什么（subchannel 模块，修 #33306、#33340）：(1) `src/scmclosures/SCMFrictionUpgradedChengTodreas.C`（由 `...UpdatedChengTodreas` 改名）：三角形栅格边/角子通道的绕丝扫流系数 Ws_T = −11·log₁₀(H/D) + 19（原为 ln）；过渡区插值因子 ψ = ln(Re_b/Re_L)/ln(Re_T/Re_L) 与区段判定改用组件整体 Re_b（新函数 `SubChannel1PhaseProblem::computeBulkReynoldsNumber`，Re_b = ṁ_in D_h,b/(μ̄_in A_b)），各子通道自身 Re 只进 f_L = C_fL·Re⁻¹、f_T = C_fT·Re^{−0.18}，并以 max(Re,1) 截断；Re_L = 320·10^{P/D−1}，Re_T = 10⁴·10^{0.7(P/D−1)}。(2) `SCMMixingChengTodreas.C` 的 Cm、Cs 区段判定与 ψ 同样改用 Re_b。(3) `SCMHTCBorishanskii.C`：Nu_T = 24.15·log₁₀(−8.12 + 12.76x − 3.65x²)，x = P/D（原为 ln）。(4) `SCMHTCGnielinski.C`：分母 1 + 12.7√f(…)(Pr_s^{2/3} − 1) 下限截到 0.1。(5) `SubChannel1PhaseProblem.C` 单体隐式求解里按 gap 索引读横流估计（`row_vec = i_gap + n_gaps·iz`）。(6) 网格生成器检查 P > D 且 D_w ≤ P − D。
- 不变量：
  - 主 A（绕丝只增不减阻力）：三角形栅格边/角子通道的绕丝修正因子 [1 + Ws_T·(A_r/A'_p)·tan²θ]^{1.41} ≥ 1 与层流因子 1 + Ws_L(A_r/A'_p)tan²θ ≥ 1，在关联式数据域 H/D ∈ [8, 52] 上成立，等价于 Ws_T = 19 − 11·log₁₀(H/D) > 0；该式在 H/D < 10^{19/11} ≈ 53.4 时为正（H/D=8 时 9.07，H/D=52 时 0.124）。旧的 ln 版本在 H/D > e^{19/11} ≈ 5.6 时为负（H/D=8 时 −3.87、52 时 −24.5），会让绕丝反而降低阻力或触发修正项为负的报错。可检形式：同一几何与流量下，带绕丝的组件压降 ≥ 光棒（D_w→0 但保持同 P、D）的组件压降。
  - 主 B（区段切换连续且全通道一致）：对固定子通道 Re_i，f(Re_b) 在 Re_b = Re_L 与 Re_b = Re_T 处连续：ψ=0 时 f = f_L(Re_i)，ψ=1 时 f = f_T(Re_i)；对 Cm、Cs 同理（ψ=0 → C_L，ψ=1 → C_T）。所有子通道共享同一个 ψ（依赖 Re_b 而非 Re_i）。数值形式：入口流量使 Re_b = Re_L(1 ± δ) 两次计算的组件压降相对差 → 0（δ→0 时线性收敛），Re_T 处同理。
  - 主 C（Borishanskii 低 Pe 极限）：Pe ≤ 200 时 Nu = 24.15·log₁₀(−8.12 + 12.76x − 3.65x²)；在适用域 x ∈ [1.1, 1.5] 上 Nu ∈ [4.25, 10.83]（x=1.2 时 6.93）。旧 ln 版本大 ln10 ≈ 2.303 倍。
  - 次 D（求解策略一致性）：#33340 修的是稳定化估计量的索引；若该稳定化在收敛时不改变离散方程，则 explicit、implicit 分离式、implicit 单体式三种求解得到同一收敛解（压降、出口温度一致到求解容差）。
  - 次 E（Gnielinski 有界性）：Re > 10³、Pr > 0 时 Nu_T 有限且 > 0，且 Nu_T 关于 Pr 连续（截断是 max(·,0.1)，连续）。
- 不变量的来源：A——第一性原理（绕丝增加形阻与扫流，不会降低轴向阻力）＋文献式（Cheng–Todreas 1986 与 Chen–Chen–Todreas 2018 的 Upgraded CTD；issue 引用 Pacio–Chen–Todreas 的 PCTD 论文，指出应为常用对数〔原文未核〕）＋数据域端点（52 < 53.4 这一吻合本身支持 log₁₀）。B——关联式构造（Cheng–Todreas 过渡区插值以组件整体 Re 定义）。C——文献式（Borishanskii 等，Todreas & Kazimi《Nuclear Systems》中液态金属棒束式，log₁₀）。D——数值方法一致性。E——数值稳健性约定。
- 常用工况域：钠冷快堆绕丝组件 P/D ∈ [1.05, 1.42]、H/D ∈ [8, 52]、D_w/D ∈ [0.05, 0.25]、N_pin ∈ {7, 19, 37, 61, 217}；Re_b ∈ [50, 10⁶]（过渡区约 400–2×10⁴）；液态金属 Pr ∈ [0.003, 0.03]、Pe ∈ [10, 2200]。依据：关联式自报数据域（源码注释与警告阈值）与 EBR-II、THORS、ORNL-19 等验证算例。
- 数值核对方法：MOOSE subchannel 模块（`subchannel-opt`，合并提交构建）。A：以 `modules/subchannel/test/tests/problems/SFR/sodium-19pin/test19_implicit.i` 为基础，`friction_closure = SCMFrictionUpgradedChengTodreas`，Re_b ≈ 5×10⁴（湍流），H/D 取 8、20、52 三档跑带绕丝；再把 `wire_diameter` 设为 1e-6 m 跑「近光棒」。读出口/入口压力 postprocessor 求 ΔP，检查 ΔP_wire ≥ ΔP_bare 且随 H/D 减小单调增大。容差：严格不等式，无需容差。B：调入口质量流量使 Re_b = Re_L·(1 ± 1e-3)、(1 ± 1e-4)，ΔP 相对差应约按 δ 线性缩小（< 10·δ），Re_T 处同理；用 `verbose_subchannel=true` 在控制台核 Re_b 值。C：`pin_HTC_closure = borishanskii`，低流量使 Pe < 200，读壁面温度与主流温度反推 Nu = h·D_h/k，与 24.15·log₁₀(poly) 比，相对 1e-6（单点代数，误差只来自温度反推的 FP）。D：同一输入的 explicit/implicit/monolithic 三版（仓库中 `test19_explicit/implicit/monolithic.i`）比较压降与出口温度，相对 1e-5（由 P_tol、T_tol 决定）。
- PR 自带测试是否已检查这条不变量：否（全部是回归）。大量 gold CSV 重生成；新增 `invalid_wire_diameter_error`（D_w > P−D 报错）；删除了 `scmclosures/friction_updated_cheng_todreas/invalid_corner_wire_correction.i` 测试（修正项为负的报错路径在 log₁₀ 下大概已无法触发，这与不变量 A 相符）。没有测 Ws>0、区段连续、Borishanskii 定值或三种求解一致。
- 你预期改动是否满足：A、C 满足（推测；有端点吻合和量级证据）；B 按式子满足。风险（不确定）：`computeBulkReynoldsNumber` 首行 `if (processor_id() != 0) return;`，其它 rank 上 `_bulk_Re` 保持初值 1.0（低于任何 Re_L，即全按层流）——若摩擦/混合闭合在非 0 rank 上也参与计算（并行运行），结果会依赖进程数；另外摩擦闭合的构造函数里读 `getBulkReynoldsNumber()` 做范围警告时，`externalSolve` 尚未调用该函数，读到的是 1.0，警告恒被触发。建议在 D 之外加一条「1 进程与 4 进程结果一致」的核对。
- 环境：合并提交 b5a1d9e3df14e98850e514f3e57d02c25c44fca9；C++、PETSc、libMesh、MOOSE＋subchannel（及 fluid_properties）模块；难度 高（完整构建，本机不宜）。
- 形式化（Lean）可行性：A 的核心是「h ∈ [8,52] ⇒ 19 − 11·log₁₀ h > 0」，等价于 52¹¹ < 10¹⁹（7.35×10¹⁸ < 10¹⁹），用 `Real.logb` 单调性＋`norm_num` 可证；进一步「Ws>0、A_r/A'_p>0 ⇒ 修正因子 ≥ 1」也容易。B 是 ψ ∈ {0,1} 端点处的代数恒等式，平凡。C 的区间界需要 log₁₀ 在有理点的包络，可行但繁琐。D、E 涉及完整求解器，不可形式化。
---

## CliMA/Oceananigans.jl#6083 Include drift velocities in the implicit part of adaptive implicit tracer advection
- 改了什么：(1) `src/Models/HydrostaticFreeSurfaceModels/hydrostatic_free_surface_{ab2,rk}_step.jl`：`AdaptiveVerticallyImplicitDiscretization`（AIVA）的 `implicit_step!` 原只收到 `model.transport_velocities`，现在收到与显式通量相同的总速度 `tracer_advecting_velocities(...)` = 输运速度 + 生化沉降速度 + 闭合辅助速度（如 GM 涡致速度）+ 平流强迫（新函数在 `hydrostatic_free_surface_tendency_kernel_functions.jl`）。(2) 夹带的修改：`advective_skew_diffusion.jl` 中 GM 平流形式的坡度改为与扩散形式同一套处理——ϵ 用 `calc_tapering`（全水平坡度 |S|² = Sx²+Sy²、b_z 下限 `minimum_bz`、b_z ≤ 0 时关闭），S = −b_h/max(b_z, b_z,min) + ∂_h z；新增 ∂x_zᶜᶠᶠ、∂y_zᶠᶜᶠ 算子。
- 不变量：
  - 主 A（速度分解无关性）：示踪物的时间推进只依赖总平流速度，不依赖它是作为输运速度、生化沉降速度还是 `AdvectiveForcing` 给出的。即对同一初值、同一 w(z) 与 Δt（包括垂向 Courant 数 > cfl 的情形），三种给法的 c(t_n) 相同（浮点舍入内）。
  - 主 B（匀速沉降的质心与质量）：单柱、w 在内部面上为常数 −w_s、上下边界面 w=0、示踪物远离边界时：Σ_k c_k Δz 守恒（相对 1e-12），质心位移 = w_s·nΔt。对 AIVA 隐式部分（一阶迎风）这一位移在离散层面严格成立（每个单元的迎风面值恰好把 Σ c_face = Σ c），显式 WENO 部分只近似成立。
  - 次 C（GM 平流形式的坡度限幅）：IsopycnalSkewSymmetricDiffusivity 平流形式下，涡流函数分量满足 |ψ_x| = |κ_skew·ϵ·S_x| ≤ κ_skew·S_max，|ψ_y| 同理（ϵ = min(1, S_max²/(S_x²+S_y²)) ⇒ ϵ·|S| ≤ S_max）；b_z ≤ 0（且 `minimum_bz` = 0）处 ψ = 0；涡致速度离散无散：δx(Δy Δz u_e) + δy(Δx Δz v_e) + δz(A_z w_e) = 0（流函数构造）。
- 不变量的来源：A——库内一致性（同一个物理平流过程，显式通量本来就用总速度，隐式部分应一致）＋第一性原理（示踪方程只含总速度 ∇·(u_total c)）。B——第一性原理（通量形式守恒；匀速平移的一阶矩方程 d/dt∫z c = ∫w c）。C——Gerdes 等（1991）式二次型限幅的构造性质＋流函数形式的无散性。
- 常用工况域：单柱 H ∈ [100, 5000] m，Δz ∈ [1, 50] m；沉降速度 w_s ∈ [1, 200] m/day（海洋颗粒有机物典型沉速）；Δt 使垂向 Courant 数 C = w_sΔt/Δz ∈ [0.1, 10]（AIVA 的设计用途就是允许 C > cfl，cfl 取 0.5）。GM：κ_skew ∈ [100, 2000] m²/s，S_max ∈ [1e-3, 1e-2]。
- 数值核对方法：Julia，Oceananigans 合并提交。A：`RectilinearGrid(size=40, z=(-400,0), topology=(Flat,Flat,Bounded))`；tracer_advection = `WENO(order=5, time_discretization=AdaptiveVerticallyImplicitDiscretization(cfl=0.5))`；(i) 仿 PR 测试的 `SinkingParticles` 生化模型给沉降速度；(ii) 同一 w 用 `HydrostaticFreeSurfaceModel(...; velocities=PrescribedVelocityFields(w=...))`；(iii) 用 `forcing=(D=AdvectiveForcing(w=...),)`。三者对 QuasiAdamsBashforth2 与 SplitRungeKutta3 各推进若干步，比较 `interior(model.tracers.D)`，容差相对 1e-12（同一算子的不同组装路径，只差浮点求和顺序；若(ii)因预设速度触发不同的时间步分支而有差异，至少应达 1e-10）。B：同上设置，Gaussian 初值中心 −100 m、宽 20 m，C ∈ {0.4, 4, 10}，检查 Σc 相对 1e-12、质心位移与 w_s nΔt 的相对误差（建议 < 1e-2；PR 用 5% 偏松），并把 cfl 设得极小使几乎全隐式，此时误差应降到 1e-6 量级（只剩边界影响）。C：三维小网格（如 16×16×16）加倾斜等密面，用 `IsopycnalSkewSymmetricDiffusivity(κ_skew=..., slope_limiter=FluxTapering(S_max), skew_flux_formulation=AdvectiveFormulation())`（参数名按合并提交核对），计算闭合场后读涡致速度 `model.closure_fields`（uₑ, vₑ, wₑ），检查离散散度在机器精度为零（相对 1e-12）；构造局部陡坡（|S| ≫ S_max）检查 |ψ| ≤ κS_max（直接调用模块内 `ϵSxᶠᶜᶠ` 等内部函数，或由 uₑ 的垂向积分反推 ψ）。
- PR 自带测试是否已检查这条不变量：部分。`test/unit/adaptive_implicit_vertical_advection.jl` 新增 testset：Courant 4、3 步，QAB2 与 SplitRK3，检查 `sum(D) ≈ sum(D₀) rtol=1e-12` 与质心下降 ≈ w_sΔt·n（rtol 0.05）——覆盖 B（容差较松），未覆盖 A 的三路等价，GM 平流形式坡度修改（C）没有新测试。
- 你预期改动是否满足：A、B 满足（推测：隐式与显式现在拿到同一个总速度对象）。C 按构造满足（推测）；但夹带的 GM 修改改变了平流形式的数值结果（限幅从只看分量改为看全水平坡度，并加了 b_z 下限），PR 描述没有提到，属于未声明的行为变化。
- 环境：合并提交 1734ff9e4509198f945c5e7740a0bccacc69a720；Julia ≥1.10＋Oceananigans.jl（CPU）；难度 中（Julia 预编译耗时、内存约 2–4 GB，本机勉强）。
- 形式化（Lean）可行性：B 的离散命题可形式化：对一维有限体积迎风格式（N 个单元、边界零通量、常速度），Σ c 守恒与一阶矩每步精确增加 w_sΔt·Σc（隐式迎风时由线性方程组 (I + ΔtA)c^{n+1} = c^n 的列和性质推出）是有限维线性代数命题，可行。C 的 ϵ|S| ≤ S_max 是简单实数不等式，平凡。A 是代码路径等价，不适合形式化。
---

## Cantera/cantera#2107 Fix thermo consistency errors
- 改了什么（一揽子修复，修 #1300、#1308、#1311、#1317、#1324、#1500、#1502）：`src/thermo/PDSS_HKFT.cpp` 新增 dVdT、dVdP（对 molarVolumeTP 中心差分，相对步长 1e-4）与 cv；`PDSS_SSVol.cpp` 中 cp 的压力修正项 −T·Δp·d²V/dT² 补除以 R，cv 改为 = cp；`IdealSolnGasVPSS.cpp` 的 cv 由 cp − R 改为 cp − T(Σx_k ∂V_k/∂T)²/(−Σx_k ∂V_k/∂p)（Σx_k∂V_k/∂p = 0 时 cv = cp），并新增 β、κ_T；`MargulesVPSSTP` 同类 cv 修复与超额体积项；`GibbsExcessVPSSTP` 密度计算去掉忽略压力项的覆写；`MolalityVPSSTP::getActivityCoefficients` 溶质项由 γ^□/X_o 改为 γ^□/(X_o M̃_o)（M̃_o 为溶剂摩尔质量 kg/mol）；`WaterProps::ADebye` 中 d(热膨胀系数)/dT 的差分由「除以 T+ΔT」的错式改为中心差分；`DebyeHuckel` 稀溶液极限下溶剂项的 ∂lnγ/∂P、∂²lnγ/∂T² 单独给式；`IdealMolalSoln::getPartialMolarEntropies` 标准态项补乘 R。测试侧：`test/data/consistency-cases.yaml` 删除大批 known-failures，`test/thermo_consistency/consistency.cpp` 新增 Mayer 关系测试，同时把 cv_eq_dudT 与 cv_eq_dsdT_const_v_times_T 对不可压缩相改为跳过。
- 不变量：
  - 主 A（定容热容，Mayer 关系）：对 κ_T > 0 的相，c_v = c_p − T·v·β²/κ_T，且 c_v = T(∂s/∂T)_v，其中 β = (1/v)(∂v/∂T)_p、κ_T = −(1/v)(∂v/∂p)_T 必须分别等于摩尔体积的有限差分；κ_T = 0（严格不可压缩）时按库约定 c_v = c_p。符号：v 为摩尔体积 [m³/kmol]，c 为摩尔热容 [J/kmol/K]。
  - 主 B（化学势—活度一致）：对所有物种 k，μ_k = μ_k^∘(T,p) + RT·ln(γ_k X_k)，μ_k^∘ 为 `standard_gibbs_RT`·RT（对 MolalityVPSSTP 族是质量摩尔浓度标准态 μ_k^□），γ_k 为 `activity_coefficients`（摩尔分数基）。等价地，溶质 γ_k X_k = γ_k^□·m_k/m°，m_k = X_k/(X_o M̃_o)，m° = 1 mol/kg。注意：教科书式 γ_x = γ_m/X_o 对应 Henry 定律摩尔分数标准态，两者差常数 M̃_o，本不变量只在「μ_k^∘ 取 μ_k^□」这一约定下成立。
  - 主 C（偏摩尔量与 Maxwell 关系）：s̄_k = −(∂μ_k/∂T)_{p,n}、Σ_k X_k s̄_k = s；c_p = (∂h/∂T)_p = T(∂s/∂T)_p；(∂s/∂p)_T = −(∂v/∂T)_p。分别对应 IdealMolalSoln（缺 R）、PDSS_SSVol 与 DebyeHuckel 可变 A_Debye 的 c_p、DebyeHuckel 稀溶液极限的压力导数这几处修复。
- 不变量的来源：第一性原理（经典热力学恒等式：Mayer 关系、Gibbs 函数的 Maxwell 关系、偏摩尔量的欧拉齐次性、化学势的活度定义）＋库内约定（ThermoPhase 的 μ_k = μ_k^∘ + RT ln(γ_k X_k) 定义）。文献：Denbigh（1981）关于质量摩尔浓度与摩尔分数活度系数换算；HKFT 模型（Helgeson–Kirkham–Flowers / Tanger–Helgeson）的标准态体积式。
- 常用工况域：水溶液相 T ∈ [273.15, 373.15] K（DebyeHuckel/HMW/IdealMolal 的常用温区；HKFT 可到 573 K）、p ∈ [1, 100] atm、离子强度 I ∈ [1e-4, 6] mol/kg（HMW-NaCl 可到饱和 ~6.1 mol/kg；测试用例取到 9.4）；熔盐 Margules（LiCl–KCl）T ∈ [700, 1400] K；固溶体 VPSS T ∈ [300, 1000] K。依据：`test/data/consistency-cases.yaml` 的状态点与各模型适用域。
- 数值核对方法：Cantera（合并提交构建）Python 接口，相文件取合并提交的 `test/data/thermo-models.yaml`、`pdss_hkft.yaml`、`debye-huckel-all.yaml`、`HMW_NaCl.yaml`、`LiKCl_liquid.yaml`、`Li_Liquid.yaml`、`IdealSolidSolnPhaseExample.yaml`。A：读 `cv_mole`、`cp_mole`、`thermal_expansion_coeff`、`isothermal_compressibility`、`volume_mole`；β、κ_T 各自与 `TPX` 设置下 v(T±h,p)、v(T,p±δ) 的中心差分比（h = 1e-4·T，δ = 1e-4·p；含 IAPWS 的相用 rtol 4e-6，其它 1e-7）；再做独立于该类自身 β、κ 的 c_v 核对：沿等体积路径 (T, p) → (T+h, p + (β/κ_T)h) 用 setState_TP 迭代把 v 拉回原值（牛顿 2–3 步），c_v ≈ T·Δs/ΔT，rtol 1e-5。B：`chemical_potentials` 与 `standard_gibbs_RT*RT + RT*log(activity_coefficients*X)` 比，rtol 1e-10；溶质项再与手算 m_k = X_k/(X_o M̃_o) 的 γ^□m_k 比（γ^□ 在 Python 层若无接口，用 C++ `MolalityVPSSTP::getMolalityActivityCoefficients`）。C：μ_k(T±h) 中心差分对比 −`partial_molar_entropies`（rtol 1e-6）；Σ X_k s̄_k 对比 `entropy_mole`（rtol 1e-12）；h(T±h) 中心差分对比 `cp_mole`；s(p±δ) 与 −v(T±h) 的差分互比。
- PR 自带测试是否已检查这条不变量：部分。`test/thermo_consistency/consistency.cpp` 参数化测试覆盖 B（chem_potentials_to_activities、activity_coeffs）、C（gk_eq_hk_minus_T_sk、s_eq_sum_sk_Xk、cp_eq_dhdT、dsdP_const_T_eq_minus_dV_dT_const_P 等）在 yaml 列出的状态点上通过，并新增 `cv_eq_cp_minus_T_V_beta2_over_kappa_T`。但 A 的覆盖有缺口：(1) 同一 PR 把 `cv_eq_dudT`、`cv_eq_dsdT_const_v_times_T` 对 `isCompressible()==false` 的相改为 GTEST_SKIP，这些相的 c_v 不再与能量/熵的有限差分比；(2) 新 Mayer 测试用的 β、κ_T 来自同一类中与 c_v 同源的 Σx_k∂V_k/∂T、Σx_k∂V_k/∂p，因此对 IdealSolnGasVPSS 近乎同义反复；β 有 `alphaV_eq_dmv_dT_const_P_div_mv` 的有限差分把关，κ_T 是否有同类测试未核。
- 你预期改动是否满足：B、C 满足（推测；测试在这些状态点已通过，修复都是明显的漏乘 R、错分母一类）。A 对可压缩相（含 IAPWS 水、HKFT）满足（推测）；对 κ_T=0 的相只按约定成立，没有物理可检的内容；对 IdealSolnGasVPSS 等「κ_T 很小但非零」的相，c_v 是否与等体积路径的 T·Δs/ΔT 一致是不确定的——v_mix = Σx_k v_k 的理想溶液假设下应一致，但该路径未被 PR 测试覆盖。
- 环境：合并提交 d58a23265414e1789536cf5332d9fb29de095a43；C++17＋SCons 构建 Cantera（依赖 Eigen、fmt、yaml-cpp、SUNDIALS、Boost 头文件）＋Python 接口（Cython）；难度 中-高（源码构建约 20–40 分钟，内存 ≥4 GB 较稳）。
- 形式化（Lean）可行性：A、C 的恒等式可在「给定满足 Maxwell 关系的偏导数实数」的抽象层写成有限维代数命题（例如由 ds = a dT + b dp、dv = c dT + e dp、b = −c、e < 0 推出 (∂s/∂T)_v = a + c²/e，再结合 a = c_p/T 得 c_v = c_p − T c²/(−e)），容易证明。要对具体模型（HKFT、IAPWS、Debye–Hückel）证明则需对其解析式做多元微分，工作量大且含有迭代求解（IAPWS 密度反演），不现实。B 在给定 γ 的定义下是对数恒等式，平凡。
---

## CliMA/ClimaAtmos.jl#4637 Fixes problems near upper boundary and tracer-mass continuity
- 改了什么：(1) `src/cache/precomputed_quantities.jl`：`set_velocity_at_top!` 由令顶层 u₃ = 0 改为 u₃ = −u_h³/g³³（新函数 `top_velocity`），使总逆变速度 u³ = u_h³ + u₃g³³ 在模式顶为 0，即使地形跟随坐标面在顶部仍倾斜（g³ʰ ≠ 0）；`implicit_tendency.jl` 质量方程改用顶底零通量的 `ᶜadvdivᵥ`。(2) `edmfx_sgs_flux.jl`：自动发现的 SGS 示踪物（q_lcl、q_icl、q_rai、q_sno、n_*）改为与 mse、q_tot 相同的差分形式 SGS 通量 ρᵏaᵏ(u³ᵏ − u³)(χᵏ − χ)，并统一用 `edmfx_sgsflux_upwinding`；`advection.jl` 中这些示踪物的格点平均垂直平流在 sgs_mass_flux 打开时也照常施加。(3) q_tot 的 SGS 质量通量同时加到 ρ（湿空气质量）。Jacobian 相应修改。
- 不变量：
  - 主 A（顶边界不可穿透）：每次状态滤波后，所有柱在顶层面上 u³ = 0（|u³_top| ≤ 100·ε·max|u_h³_top|），底层面同理；由此无地表通量、无降水的干绝热运行中全域总质量 M = Σ ρ J 守恒：|M(t_n) − M(0)|/M(0) ≤ C·n·ε（地形坐标面在顶部倾斜时也成立）。
  - 主 B（示踪物—质量一致，即「χ≡1 检验」）：若某比湿型示踪物在格点平均与所有子域中都等于常数 c（ρχ = cρ，χʲ = c），则每条输运路径给出 ∂(ρχ)/∂t = c·∂ρ/∂t（逐点，到舍入），其中差分形式 SGS 通量逐点为零；推论：无源汇时均匀的 χ 在时间推进后保持均匀。
  - 次 C（库约定：SGS 水通量只改变湿空气质量）：`edmfx_sgs_mass_flux_tendency!` 给出的 ∂ρ/∂t 与 ∂(ρq_tot)/∂t 逐位相等，即干空气质量 ρ_d = ρ − ρq_tot 不受 SGS 水通量影响。
- 不变量的来源：A——第一性原理（刚盖顶边界不可穿透：法向质量通量 ρu³J = 0）＋离散守恒（有限体积通量形式的望远镜求和、水平 DSS 守恒）。B——离散一致性（示踪物方程与连续方程在 χ 为常数时须退化为同一方程，Lin–Rood/Gross 等「free-stream preservation」要求）。C——库内一致性（PR 称与扩散通量的处理对称）；其物理正确性见下文。
- 常用工况域：平面/球面配置；地形高度 0–3 km（Schär 山、Agnesi 山），模式顶 z_top ∈ [20, 45] km，线性网格扭曲（顶部仍倾斜，最能触发 A）；水平风 0–30 m/s；Δt ∈ [1, 120] s；EDMF 单柱 DYCOMS/BOMEX 类（z ≤ 3 km，Δt 60–120 s）；Float64 与 Float32 都应检查（B 的 PR 自报量级在 Float32 下）。
- 数值核对方法：Julia，ClimaAtmos 合并提交，仿 `test/prognostic_equations/advection_tests.jl` 与 `tracer_mass_consistency_tests.jl` 用 `CA.AtmosConfig(Dict(...))` 建模拟（plane、Schar、mesh_warp_type="Linear"、z_max 21 km）。A：读 `p.precomputed.ᶠu³` 顶层与 `CA.compute_ᶠuₕ³`，检查比值；再用干模型、关地表通量（`disable_surface_flux_tendency=true`、无微物理）推进 1000 步，每步用 ClimaCore 的 `sum(Y.c.ρ)`（带 J 的体积分）记录 M，相对漂移应 ≤ 1e-12（Float64）且不随时间线性增长；对照：把 `mesh_warp_type` 改为顶部变平的 SLEVE，两者应同样守恒。B：在任意已推进的状态上把 ρχ ← cρ、χʲ ← c（c 取 1 与 1e-3 两档，c≠1 可排除「恰好等于 ρ」的巧合），分别调用 `CA.implicit_tendency!`、`CA.explicit_vertical_advection_tendency!`、`CA.edmfx_sgs_mass_flux_tendency!`，比较 Yₜ.c.ρχ 与 c·(相应 ρ 路径)，容差 100·ε·max|ρ 倾向|；再做时间推进版：均匀 χ 初值、关微物理源汇，推进 N 步后 max|χ − c|/c ≤ N·100ε。C：同一调用比较 `parent(Yₜ.c.ρ) == parent(Yₜ.c.ρq_tot)`。
- PR 自带测试是否已检查这条不变量：部分。A：`advection_tests.jl` 检查顶/底 u³ ≤ 100ε·max|u_h³|（单状态），未检查长时间总质量守恒。B：`tracer_mass_consistency_tests.jl` 对 χ≡1 检查隐式倾向 ρq_tot 与 ρ、显式平流各凝结物与 ρ 平流、EDMF SGS 通量为零（tol 1e-10·flux_scale），只用 c = 1，且只检查单次倾向、不检查推进后的均匀性保持。C：同文件 `@test parent(Yₜ.c.ρ) == parent(Yₜ.c.ρq_tot)`。
- 你预期改动是否满足：A、B 满足（推测；修改正是为此，测试已覆盖单步）。C 按代码满足，但物理上不确定：若格点平均速度取质量加权（重心）速度，则 Σ_k ρᵏaᵏ(u³ᵏ − u³) = 0，SGS 输运的总湿空气质量通量应为零、干空气的 SGS 通量等于水通量的相反数；把 q_tot 的 SGS 通量加到 ρ 相当于忽略干空气的 SGS 输运（与扩散通量的同类近似一致，是建模选择而非恒等式）。另：PR 自己在 `remaining_tendency.jl` 留 TODO 指出海绵层对凝结物的阻尼没有对应的 ρq_tot/ρ/ρe_tot 项，会破坏 B 的推论（均匀 χ 在海绵层内不保持），核对时需关闭海绵层或把检查区限制在海绵层以下。
- 环境：合并提交 a8efce5a0040924094a58739db91a0f3162f779c；Julia ≥1.10＋ClimaAtmos 及 ClimaCore、ClimaTimeSteppers 等；难度 中-高（依赖多、首次预编译很久，内存 4–8 GB，本机吃紧，建议放计算节点）。
- 形式化（Lean）可行性：A 的离散守恒可写成有限维命题：一维柱 N 个单元，质量倾向 −(F_{k+1/2} − F_{k−1/2})/ΔV_k，F_{1/2} = F_{N+1/2} = 0 ⇒ Σ ΔV_k ∂ρ_k/∂t = 0（望远镜求和）；顶层 u³ = u_h³ + u₃g³³ = 0 由 u₃ = −u_h³/g³³（g³³ > 0）直接得出。B：当 χ_k = c 时差分通量 ρᵏaᵏ(u³ᵏ−u³)(χᵏ−χ) 中 (χᵏ−χ)=0 的插值/迎风重构仍为 0，需要假设迎风重构对常数精确（线性重构的常数再现性），可在有限维上证明。整套 ClimaCore 算子不宜形式化。
---

## CliMA/ClimaAtmos.jl#4822 Compute horizontal gradient along isentropes in SGS variance
- 改了什么：`src/parameterized_tendencies/microphysics/cloud_fraction.jl::set_covariance_cache!`：水平已分辨梯度（几何）SGS 方差项由「θ′θ′ += geo_h|∇_hθ_li|²、q′q′ += geo_h|∇_hq|²、T′q′ += geo_h ∇_hθ·∇_hq」改为只加到 q′q′ 的等熵形式 q′q′ += geo_h·f·I，I = |∇_hq − r∇_hθ_li|²（展开式 `isentropic_gradient_invariant`，max(·,0) 防舍入），r = q_z/max(θ_z, θz_min)（`isentropic_slope_ratio`），f = 1/(1 + (sΔx/Δz)²)，s = |∇_hθ_li|/max(θ_z, θz_min)（`isentropic_cell_cap`）；垂直均匀盒块的 q′q′ 与 θ′q′ 项乘同一 f，θ′θ′ 不变。新参数 `sgs_variance_isentropic_min_dtheta_dz` 默认 1e-3 K/m。
- 不变量：
  - 主 A（沿等熵面导数，规范不变性）：在 θ_z ≥ θz_min 的稳定层内，I 等于 q 沿 θ_li 等值面的水平梯度平方：若等值面为 z = h(x,y)，则 ∇_s q = ∇_hq + q_z∇h = ∇_hq − (q_z/θ_z)∇_hθ。推论（可直接检）：(i) 若 q = F(θ)（q 只是 θ 的函数），I = 0；(ii) 对任意光滑 G，I(q + G(θ)) = I(q)；(iii) 取 F、G 为线性（q = a + bθ），离散梯度算子是线性的，故离散层面 I 也应精确为 0（到舍入）。
  - 主 B（不重复计数与几何上界）：对任意 s ≥ 0，f ∈ (0,1]、随 s 单调递减；记 ℓ² = fΔx²，则恒等式 (sℓ)² + fΔz² = Δz²（沿面段已采样的垂直位移与垂直块剩余份额之和恰为 Δz²）；且水平块中「投影部分」r²|∇_hθ|²ℓ² ≤ (q_zΔz)²。
  - 主 C（相关系数的 Cauchy–Schwarz 界）：`tq_correlation_model = diagnosed` 时，各块协方差矩阵半正定：水平块 [[0,0],[0,geo_h f I]]；垂直块 c[[θ_z², fθ_zq_z],[fθ_zq_z, fq_z²]]，其行列式 c²θ_z²q_z²(f − f²) ≥ 0；湍流块在 |ρ_turb| ≤ 1 时半正定；θ→T 的线性变换保持半正定。故钳位前 |T′q′| ≤ √(T′T′·q′q′)（只差舍入），下游 clamp 只应在 O(ε) 上起作用。
- 不变量的来源：A——几何恒等式（隐函数求导，等值面上的方向导数）＋物理论证（大尺度搅拌近似绝热，气块沿等熵面移动，PR 描述与 docstring）。B——构造恒等式（docstring 给出 Δz² − (sℓ)² = Δz²/(1 + (sΔx/Δz)²)）。C——线性代数（半正定矩阵之和半正定）＋ docstring 中的 Cauchy–Schwarz 论证。
- 常用工况域：Δx ∈ [1, 250] km（GCM 到区域模式；PR 测试用 2.5e5 m），Δz ∈ [20, 1000] m；θ_z ∈ [−5e-3, 1e-2] K/m（含中性/逆温层，θz_min = 1e-3 K/m 起作用），|∇_hθ| ∈ [0, 1e-4] K/m，q_z ∈ [−1e-5, 0] 1/m，|∇_hq| ∈ [0, 1e-8] 1/m；由此 s Δx/Δz 覆盖约 [0, 50]，f 从 1 到 ~4e-4。Float32 与 Float64。
- 数值核对方法：Julia，ClimaAtmos 合并提交。逐点（低难度）：直接调用 `CA.isentropic_slope_ratio`、`CA.isentropic_cell_cap`、`CA.isentropic_gradient_invariant`（PR 测试已通过 `CA.` 前缀访问）。A：随机抽 ∇_hθ、θ_z ∈ [θz_min, 1e-2]、b，令 ∇_hq = b∇_hθ、q_z = bθ_z，算 inv_q = |∇_hq|²、inv_θq = ∇_hθ·∇_hq、inv_θ = |∇_hθ|²，检查 I ≤ 10ε·(inv_q + r²inv_θ)；再对任意 ∇q 检查 I(∇q + b∇θ, q_z + bθ_z) 与 I(∇q, q_z) 的差 ≤ 10ε·(各项量级)（展开式有抵消误差，Float32 下 I ≪ |∇q|² 时相对误差会很大，这是需要报告的数值性质，不是不变量失败）。B：对 s 网格算 f、ℓ² = fΔx²，检查 (sℓ)² + fΔz² = Δz²（相对 1e-6/Float32、1e-14/Float64）与投影上界。整场（中难度）：建 plane/sphere 配置，`tq_correlation_model = "diagnosed"`，把状态设成 q_tot = a + b·θ_li（需按热力学参数由 θ_li 与 q 反解 T、ρ，或先设 θ_li 场再设 q），调用 `CA.set_covariance_cache!(Y, p, thermo_params)`，比较 `sgs_variance_horizontal_scale_factor` 取 0 与非 0 时的 q′q′：差应 ≤ 100ε·q′q′。C：在一般状态上读 `p.precomputed.ᶜT′q′`、T′T′、q′q′（clamp 之前的量），检查 (T′q′)² ≤ T′T′·q′q′·(1 + 100ε)。
- PR 自带测试是否已检查这条不变量：部分。`test/parameterized_tendencies/microphysics/cloud_fraction.jl` 新增两组：`isentropic_slope_ratio` 的取值与下限行为（稳定、临界、中性、逆温、符号、有界、类型），`isentropic_cell_cap` 的平坦/陡峭极限、单调性、下限、(0,1] 有界与类型。没有测 A（q = F(θ) ⇒ I = 0、规范不变性）、B 的恒等式、C 的相关系数界，也没有对 `isentropic_gradient_invariant` 的单测。
- 你预期改动是否满足：A、B、C 满足（推测；都是代数构造性质，代码与式子逐项对应）。需要注意（不确定）：(1) I 用展开式计算，有效数字在近乎对齐时损失，Float32 下 q = a + bθ 的整场检验可能只达到 1e-3 量级的相对精度；(2) θ_z < θz_min 时 r、s 都被下限替代，A 不再成立（设计如此）；(3) 水平块尺度用 c_Δx·Δx，而 f 里用 Δx_h，c_Δx ≠ 1 时两者不一致。物理上「等熵搅拌只产生 q 方差、不产生 θ 方差」是建模假设，没有可数值检验的第一性原理对照。
- 环境：合并提交 ba2bf1225fc48058d709853c065b1db21a412204；Julia ≥1.10＋ClimaAtmos（及 Thermodynamics、ClimaCore）；逐点检查难度 低（只需加载包），整场检查 中-高。
- 形式化（Lean）可行性：高。A(i)(ii) 是 ℝ² 向量的代数恒等式：|(g_q + b g_θ) − ((q_z + bθ_z)/θ_z) g_θ|² = |g_q − (q_z/θ_z)g_θ|²（θ_z ≠ 0），`field_simp; ring` 可证；B 的 (sℓ)² + fΔz² = Δz² 与投影上界是单变量有理函数恒等式/不等式；C 的 2×2 行列式非负（f ∈ [0,1] ⇒ f − f² ≥ 0）与半正定矩阵和的 Cauchy–Schwarz 也都是有限维实数命题。浮点抵消误差的界不在实数形式化范围内。
