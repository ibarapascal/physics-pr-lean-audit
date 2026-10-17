## Cantera/cantera#2107 Fix thermo consistency errors

- What changed (1–3 sentences, naming files and formulas):

  `src/thermo/IdealSolnGasVPSS.cpp`、`MargulesVPSSTP.cpp` 及 `PDSS_{ConstVol,SSVol,HKFT,Water}.cpp` 增加体积温压导数并修正相级热容为 \(c_v=c_p-TV\alpha^2/\kappa_T\)（零压缩率分支另行处理），`PDSS_SSVol.cpp` 将无量纲热容中的压力修正改为 \(-T(P-P_r)V''(T)/R\)，删除 `GibbsExcessVPSSTP::calcDensity` 覆盖以使用偏摩尔体积计算密度。`MargulesVPSSTP.cpp` 将 \((P-P_r)V^E\) 纳入活度及其导数，`MolalityVPSSTP.cpp` 修正两种活度系数基准的换算，`IdealMolalSoln.cpp` 补上标准态熵的 \(R\)，`DebyeHuckel.cpp` 与 `WaterProps.cpp` 修正溶剂活度的温压导数和水性质的中心差分。`BinarySolutionTabulatedThermo.cpp/.h` 统一端点对数截断并解释参考态约定，另有测试及构建路径修改；以下选择前三类关键物理不变量，不将表格模型端点约定或构建路径当成新的物理定律。

- Invariants (mathematical statements; write them as checkable equalities/inequalities, give variables, sign conventions, domain of validity):

  **共同记号。** \(T>0\) 为绝对温度，\(P>0\) 为绝对压力，\(x_k\geq0,\sum_kx_k=1\)；所有温压偏导固定组成，禁止在扰动过程中调用化学平衡。\(V\) 为相的摩尔体积（m³/kmol），\(\rho\) 为质量密度，\(\overline W=\sum_kx_kW_k\)（kg/kmol），\(g,h,u,\mu_k\) 的单位为 J/kmol，\(s,c_p,c_v\) 为 J/kmol/K，\(R\) 取 `Cantera::GasConstant`。上标 \(\circ\) 是当前 \(T,P\) 的标准态，横线是偏摩尔量。化学势采用 `getChemPotentials` 的化学部分，电势固定为零。导数要求局部单相、光滑且不跨截断边界。

  **I1 — primary：可压缩液相的热容差必须由实际体积响应决定。** 定义
  \[
  V=\overline W/\rho,\qquad V_T=(\partial_TV)_{P,x},\quad V_P=(\partial_PV)_{T,x},\qquad
  \alpha=V_T/V,\quad\kappa_T=-V_P/V.
  \]
  对 `ideal-solution-VPSS` 或 `Margules` 的机械稳定分支 \(V>0,V_P<0\)，应有
  \[
  c_p-c_v=-T\frac{V_T^2}{V_P}=TV\frac{\alpha^2}{\kappa_T}\geq0,
  \qquad c_v=\left(\frac{\partial u}{\partial T}\right)_{V,x}.
  \]
  不能把各物种的 \(c_p-c_v\) 简单按组成相加：先求混合物 \(V_T,V_P\)，再作平方和相除。对体积在邻域内同时与 \(T,P\) 无关的真正恒体积模型，退化关系是 \(c_v=c_p\)。如果 \(V_P=0,V_T\ne0\)（SSVol，或具有温度依赖过剩体积的恒标准体积 Margules），固定体积下不存在通常的局部变温路径；此时返回 `cp` 是库的约定，不能将其当作上述物理等式的证据，也不能把公式中的分母替换成小正数继续检验。单个 HKFT 溶质的标准态偏摩尔体积可以为负，不对它强加纯物质的稳定性不等式。

  **I2 — primary：压力体积项和活度项必须来自同一个 Gibbs 势。** 针对 Margules 过剩体积、SSVol 热容、Debye–Hückel 温压导数及 IdealMolal 熵修正，检查以下相互关联的 Gibbs 导数闭合：
  \[
  \bar V_k=(\partial_P\mu_k)_{T,x},\quad
  \bar s_k=-(\partial_T\mu_k)_{P,x},\quad
  \mu_k=\bar h_k-T\bar s_k,\quad
  \bar c_{p,k}=(\partial_T\bar h_k)_{P,x}=T(\partial_T\bar s_k)_{P,x},
  \]
  \[
  V=\sum_kx_k\bar V_k=\overline W/\rho,\qquad
  (\partial_Ps)_{T,x}=-V_T,\qquad c_p=(\partial_Th)_{P,x}=T(\partial_Ts)_{P,x}.
  \]
  为避免仅用互相共享错误的库函数作比较，使用两个明确的独立解析特例：

  对一个二元 Margules 相，取 \(x=x_A,y=x_B=1-x\)，\(P_r=101325\) Pa，并定义
  \[
  C_j=H_j-TS_j+(P-P_r)(V_{Hj}-TV_{Sj}),\quad
  g^E=xy(C_0+C_1y),\quad
  Q_j=V_{Hj}-TV_{Sj}.
  \]
  输入参数 \(H_j,S_j,V_{Hj},V_{Sj}\) 分别使用 J/kmol、J/kmol/K、m³/kmol、m³/kmol/K。则
  \[
  V^E=xy(Q_0+Q_1y),\quad
  \bar V_A-V_A^\circ=y^2[Q_0+Q_1(1-2x)],\quad
  \bar V_B-V_B^\circ=x^2[Q_0+2Q_1y],
  \]
  \[
  g-\sum_kx_kg_k^\circ-RT\sum_kx_k\ln x_k=g^E,\quad
  h-\sum_kx_kh_k^\circ=xy[H_0+H_1y+(P-P_r)(V_{H0}+V_{H1}y)],
  \]
  \[
  s-\sum_kx_ks_k^\circ+R\sum_kx_k\ln x_k
   =xy[S_0+S_1y+(P-P_r)(V_{S0}+V_{S1}y)].
  \]
  在恒标准体积特例中，\(V_T=-xy(V_{S0}+V_{S1}y)\)。体积收缩允许 \(V^E<0\)，只要求总 \(V>0\)；不能错误地要求每个过剩量为正。

  对 SSVol 的压力无关标准体积 \(V_k^\circ(T)\)，应有
  \[
  c_{p,k}^\circ(T,P)-c_{p,k}^\circ(T,P_r)
  =-T(P-P_r)\frac{d^2 V_k^\circ}{dT^2}.
  \]
  这里 \(P_r\) 使用物种实际参考压力，不能一概假定它等于 Margules 的 1 atm。对于 volume-temperature-polynomial，可直接对输入多项式求二阶导；对于 density-temperature-polynomial，应先写 \(V=W/\rho(T)\) 再求导。该等式直接检验遗漏 \(R\) 的量纲错误。

  对 Debye–Hückel 的淡溶液，记 \(L_k=\ln\gamma_k^{(m)}\)，固定组成时还应有
  \[
  \bar V_k-V_k^\circ=RT L_{k,P},\qquad
  \bar c_{p,k}-c_{p,k}^\circ=-R(2T L_{k,T}+T^2 L_{k,TT}).
  \]
  这两式尤其要检查溶剂 \(k=0\)：它的 \(L_0\) 包含与 Debye 参数 \(A\) 无关的组成项，不能将整个 \(L_0\) 乘 \(A_P/A\) 或 \(A_{TT}/A\) 作为导数。此检查通过相级性质覆盖 `WaterProps` 修正，无须访问私有导数数组。

  **I3 — secondary：molality 与 mole-fraction 活度基准转换不能改变化学势。** 溶剂为索引 0，\(M_0=W_0/1000\) 为 kg/mol，\(m^\circ=1\) mol/kg；在无溶剂下限截断的区域，对每个溶质 \(k>0\)，
  \[
  m_k=\frac{x_k}{M_0x_0},\quad
  \gamma_k^{(x)}=\frac{\gamma_k^{(m)}}{M_0m^\circ x_0},\quad
  a_k=x_k\gamma_k^{(x)}=\frac{m_k}{m^\circ}\gamma_k^{(m)},\quad
  \frac{\mu_k-\mu_k^\circ}{RT}=\ln a_k.
  \]
  溶剂则为 \(\gamma_0^{(x)}=\gamma_0^{(m)}\)、\(a_0=x_0\gamma_0^{(m)}\)，不能给溶剂套用溶质换算因子。`getMolalities` 的数值按 mol/kg 解释，尽管其头文件局部注释有不一致；以实现采用的 \(W_0/1000\) 及上述定义为准。理想 molal 溶质 \(\gamma_k^{(m)}=1\) 不意味着 \(\gamma_k^{(x)}=1\)；水溶剂无限稀释时后者约为 55.5，这不是异常。

- Sources of the invariants (first principles / literature formula / consistency with other implementations in the library / dimensions):

  I1 由固定组成的 \(dg=-s\,dT+V\,dP\)、Maxwell 关系及沿 \(dV=0\) 的链式法则导出，稳定性只用于 \(V_P<0\) 的符号判定。I2 来自同一 Gibbs 势的导数和总量对物质量的一次齐次性；Margules 解析式由输入定义的 \(g^E\) 对压力、温度和物质量求导，而非重复调用生产导数函数。I3 来自 molality 定义和无量纲活度约定；显式 \(m^\circ\) 避免对有量纲的 molality 取对数。

  对照合并提交的公开源文件：[Margules 实现](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/src/thermo/MargulesVPSSTP.cpp)、[继承后的密度实现](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/src/thermo/VPStandardStateTP.cpp)、[Debye–Hückel 实现](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/src/thermo/DebyeHuckel.cpp)、[molality 公共接口与约定](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/include/cantera/thermo/MolalityVPSSTP.h)。这些来源确定模型和 API 语义，不是实验真实性的独立验证。

- Common operating domain (range of each input, with the basis):

  所有测试只改变 \(T,P,x\) 或等价的 molality，物种与模型参数固定为合并提交的 YAML。以下范围是保守的抽查设计，不声称覆盖各状态方程的全部有效域；状态和差分邻域都必须落在物种温区且保持液态。

  | 用例与输入范围 | 依据与限制 |
  |---|---|
  | I1：`test/data/pdss_hkft.yaml`，相 `vpss_soln_pdss_hkft_phase`；\(T=300\)–350 K，\(P=1\)–10 atm；\(x_{Na+}=x_{Cl-}=q\)，\(q=10^{-4}\)–0.02，\(x_{H+}=x_{OH-}=10^{-7}\)，水补足至 1 | 水标准态声明 273.15–600 K；该窗口避开沸腾、临界点与 HKFT 的 228 K 奇点，选择稀水溶液以避免用高浓度理想溶液近似作实验判断。检查实际 \(V>0,V_P<0\) 后才用热容差不等式。 |
  | I1 恒体积对照：`IdealSolidSolnPhaseExample.yaml` / `VpssSolidSolutionExample`；\(T=300\)–400 K、\(P=1\)–5 atm、正组成且每种 \(x_k\geq0.01\) | 与既有一致性案例同一区域；仅在确认所有标准体积恒定且无温变过剩体积后采用 \(c_v=c_p\)。 |
  | I2 Margules：`thermo-models.yaml` / `ethanol-water-rough`；\(T=298.15\)–308.15 K、\(P=1\)–20 atm、\(x_A=0.05\)–0.95 | 温区来自新增用例；扩大压力范围专门激活 \(P-P_r\) 项。YAML 明确是虚构的粗略测试模型，不是乙醇水实验拟合。参数固定为文件值；例如 \(V_{H0}=-12.707\)、\(V_{H1}=0.476\) cm³/mol，\(V_{S0}=-0.030\)、\(V_{S1}=0\) cm³/mol/K，换成 m³/kmol 时乘 \(10^{-3}\)。 |
  | I2 SSVol：`Li_Liquid.yaml` / `Li(L)`；\(T=460\)–500 K、\(P=1\)–10 atm、\(x_{Li(L)}=0.1\)–0.9，另一物种补足 | 位于新增 300–500 K 测试区间的上部，避开低温液态锂外推；这是模型一致性测试，不保证模型参数的实测准确度。专门取 \(P\ne P_r\) 且 \(V''\ne0\) 的点。 |
  | I2/I3 电解质：`debye-huckel-all.yaml` / `debye-huckel-dilute-IAPWS`，并以 `debye-huckel-dilute` 作常数 \(A\) 对照；\(T=300\)–350 K，\(P=1\)–10 atm；\(m_{Na+}=m_{Cl-}=10^{-4}\)–\(10^{-2}\) mol/kg，\(m_{H+}=m_{OH-}=10^{-7}\)，其他溶质为零 | 满足电中性，离子强度 \(I\lesssim0.011\) mol/kg，符合 dilute-limit 的低浓度用途且远离文件设置的 \(I=3\) 截断。原测试中约 9 mol/kg 的盐浓度只能作回归点，不能作为淡溶液物理有效范围的依据。只对正浓度物种取对数。 |
  | I2/I3 理想 molal：`thermo-models.yaml` / `ideal-molal-aqueous`；\(T=300\)–330 K，\(P=1\)–5 atm；CO2(aq)、H2S(aq)、CH4(aq) 各 \(10^{-3}\)–0.1 mol/kg | 围绕原案例的 300/330 K；保持水为绝大多数、避开低溶剂分数截断。这里的假物种热力学数据用于单位与基准一致性验证。 |

  用例及原始状态依据：[consistency-cases.yaml](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/test/data/consistency-cases.yaml)、[新增 Margules 参数及假物种说明](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/test/data/thermo-models.yaml)、[HKFT 标准态定义](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/test/data/pdss_hkft.yaml)、[Debye–Hückel 配置](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/test/data/debye-huckel-all.yaml)。

- Numerical check method (which public API of the merged commit, what to compute, what to compare against, suggested tolerance and why; if finite differences are needed, the step-size strategy):

  **入口。** 后续验证使用合并提交的 C++ API：`Cantera::newThermo(file, phaseName)`；通过 `setMoleFractionsByName` 或转为 `MolalityVPSSTP&` 后的 `setMolalitiesByName` 设置组成，再用 `setState_TP(T,P)` 设置状态。记录实际返回的组成，避免归一化与电中性调整被误认为输入。使用文件的完整路径或 `addDataDirectory` 解析引用；每次扰动明确重设 \(T,P\)，每组计算后恢复中心状态，并交错属性查询顺序检查水状态缓存是否影响结果。不得仅因 `isCompressible()` 返回 false 就跳过实际有压力响应的 VPSS 模型。

  **I1。** 读取 `cp_mole()`、`cv_mole()`、`thermalExpansionCoeff()`、`isothermalCompressibility()`、`meanMolecularWeight()`、`density()`、`molarVolume()`。独立用 \(\overline W/\rho(T\pm\delta T,P)\)、\(\overline W/\rho(T,P\pm\delta P)\) 得到 \(V_T,V_P\)，先核对 \(\alpha,\kappa_T\)，再与 \(c_p-c_v\) 比较；仅将公开 \(\alpha,\kappa_T\) 代回生产 `cv` 同一公式是不够的。再独立在 \(T\pm\delta T\) 用一维求根调节 \(P_\pm\)，使 `molarVolume()` 等于中心 \(V\)，从 `intEnergy_mole()` 的中心差分核对 \(c_v\)。初始压力预测为 \(P_\pm=P\mp(V_T/V_P)\delta T\)，采用有界求根，必要时减小 \(\delta T\) 保持正压力及液态；不要依赖 VPSS 未必支持的 `setState_TD`。体积根残差目标相对 \(10^{-11}\)，若水求解器无法达到，记录误差界而非声称该导数已通过。真恒体积对照只检查 \(c_p=c_v\) 和 \(c_p=\partial_Tu\)。

  **I2。** 用 `getChemPotentials`、`getPartialMolarVolumes`、`getPartialMolarEntropies`、`getPartialMolarEnthalpies`、`getPartialMolarCp`、`getStandardVolumes`、`getStandardChemPotentials`、`getEnthalpy_RT`、`getEntropy_R`、`getCp_R` 取得中心态及扰动态数组；最后三者须分别乘 \(RT,R,R\) 还原单位。逐物种比较 \(\mu_P\) 与 \(\bar V\)、\(-\mu_T\) 与 \(\bar s\)、\(\partial_T\bar h\) 与 \(\bar c_p\)，同时检查相级 Maxwell 关系和体积加和。Margules 的 \(g^E,h^E,s^E,V^E,\bar V^E\) 用上述输入参数多项式独立计算，在 1 atm 和至少两个非参考压力比较；不能只检查总 \(h\) 而让小的过剩项被大标准态焓掩盖。SSVol 直接用输入多项式求 \(V''\)，比较两个压力的 `getCp_R()*R` 差。Debye 的两个额外等式用 `getMolalityActivityCoefficients` 的对数中心差分构造 \(L_P,L_T,L_{TT}\)，专门报告溶剂残差；相比只检查混合物 \(c_p\)，这样更能定位溶剂项错误。

  **I3。** 同一状态读取 `getActivities`、`getActivityCoefficients`、`MolalityVPSSTP::getMolalityActivityCoefficients`、`getMolalities`、`getMoleFractions`、`molecularWeight(0)` 以及两种化学势。直接验证 I3 的全部换算，化学势比较使用 \((\mu-\mu^\circ)/(RT)-\ln a\)，避免指数溢出。在至少两种稀释度重复；不能只测试纯溶剂，因为错误的溶质换算因子在那里不可见。

  **差分与容差。** 温度一阶导数从 \(\delta T/T=10^{-4}\) 开始，压力一阶导数从 \(\delta P/P=10^{-3}\) 开始，各试 \(h/2,h,2h\)；光滑情况下用 \([f(q+h)-f(q-h)]/(2h)\)，二阶导数用 \([f(q+h)-2f(q)+f(q-h)]/h^2\)，后者从 \(\delta T=0.1\) K 开始。只在邻域有效且出现收敛平台时接受结果；利用步长差估计截断误差，过小时出现消减误差则回退，不能机械缩步。HKFT 内部采用相对 \(10^{-4}\) 差分，WaterProps 内部温度步长为 0.1 K，故不要求机器精度导数。

  无差分的 I3 乘法关系建议 \(10^{-10}\) 相对误差加 \(10^{-14}\) 无量纲绝对误差；化学势对数残差建议 \(10^{-9}\)。能量代数闭合建议 \(|r|\le10^{-5}\) J/kmol 加 \(10^{-12}\) 乘参与运算各项绝对值之和，以容纳标准态大数消减。解析过剩体积与密度加和建议绝对 \(10^{-11}\) m³/kmol 加 \(10^{-9}\) 相对误差。一阶导数默认 \(10^{-6}\) 相对误差，含水/HKFT 放宽为 \(10^{-5}\)；热容绝对底限 \(10^{-3}\) J/kmol/K，体积导数所得 \(\mu_P\) 绝对底限 \(10^{-8}\) m³/kmol，\(\kappa_T\) 绝对底限 \(5\times10^{-14}\) Pa⁻¹，\(\alpha\) 绝对底限 \(10^{-9}\) K⁻¹。二阶活度导数所得热容修正可从 \(10^{-4}\) 相对、\(10^{-2}\) J/kmol/K 绝对容差起步，但必须同时通过一阶焓导数检查；近零的修正应选较高但仍稀的浓度，使信号高于数值误差。容差是后续试验建议，尚未运行校准；步长不收敛或信号低于误差底限时标记“数值未决”，不能计为通过。

- Do the PR's own tests already check this invariant (yes/partial/no, naming test files and assertions):

  **I1：partial。** `test/thermo_consistency/consistency.cpp` 新增 `cv_eq_cp_minus_T_V_beta2_over_kappa_T`，并有 `betaT_eq_minus_dmv_dP_const_T_div_mv`、`alphaV_eq_dmv_dT_const_P_div_mv` 的有限差分检查；新增 `MargulesSSVol`、`MargulesWithExcessVolume` 实例，HKFT 设置 `rtol_fd=4e-6`、`atol_c=5e-14`。但 `cv_eq_dudT` 与 `cv_eq_dsdT_const_v_times_T` 对 `!isCompressible()` 直接跳过，新增热容差断言与生产代码共享公式，而且零 \(\kappa_T\) 分支没有要求 \(V_T=0\)。因此不是本文限定域内独立恒体积路径的完整验证。

  **I2：partial。** 同文件有 `gk_eq_hk_minus_T_sk`、`v_eq_sum_vk_Xk`、`cp_eq_dhdT`、`cp_eq_dsdT_const_p_times_T`、`cpk0_eq_dhk0dT` 和 `dsdP_const_T_eq_minus_dV_dT_const_P`。`consistency-cases.yaml` 删除 IdealMolal、Debye 可变 \(A\) 与 HMW 的相关已知失败标记；`test_problems/cathermo/ims/output_blessed.txt` 更新熵回归值，`test_problems/cathermo/HMW_graph_CpvT/output_blessed.txt` 更新热容结果。新增 `ethanol-water-rough` 的三个中心态都在 1 atm（Maxwell 检查仍会微扰压力），未明确覆盖宽压力范围的解析过剩量，也未逐物种用 \(\mu_P\) 核验偏摩尔体积；更新参考输出本身不构成独立物理证明。

  **I3：yes（核心恒等式），partial（独立基准换算）。** `activity_coeffs` 断言 \(a_k=\gamma_k^{(x)}x_k\)，容差为 \(10^{-10}\) 相对加 \(10^{-20}\)；`chem_potentials_to_activities` 断言 \(\exp[(\mu_k-\mu_k^\circ)/(RT)]=a_k\)，容差 \(10^{-9}\) 相对加 \(10^{-14}\)。修订的 molal 相不再将此标为已知失败，但没有在这些断言中独立通过 \(W_0/1000\) 和 \(\gamma^{(m)}\) 重建完整换算。以上是源码覆盖判断，不是本次执行结果。测试依据：[合并版 consistency.cpp](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/test/thermo_consistency/consistency.cpp)。

- Do you expect the change to satisfy it (satisfies / may not satisfy / uncertain), with one or two sentences of reasoning (speculation, not a conclusion):

  **I1：satisfies（推测，仅上述稳定可压缩域及真恒体积域）**，新公式正确使用混合物体积导数；零压缩率但非零热膨胀的模型中，`cv=cp` 的物理解释 **uncertain**，不能由这段实现证明。**I2：satisfies（推测）**，加入 \((P-P_r)V^E\)、补 \(R\) 和拆出溶剂专属导数均符合 Gibbs 导数与量纲，剩余风险是内部水性质数值误差及未充分覆盖的非参考压力；**I3：satisfies（推测）**，除以 \(x_0W_0/1000\) 正是活度基准换算所需因子。

- Environment: merge commit SHA, language and dependencies needed, estimated setup difficulty (low/medium/high):

  合并 SHA：`d58a23265414e1789536cf5332d9fb29de095a43`。首选 C++20 的 Cantera 库和上述同提交 YAML（含被引用的物种文件），用 C++ 公共接口即可，不需要 Python 绑定；以后若从该 SHA 构建，一般需要 Python/SCons 构建驱动及 Cantera 的 Eigen、fmt、yaml-cpp、SUNDIALS 等依赖，重跑原 GTest 套件还需 GoogleTest/Boost 头文件，精确版本依该提交的构建检查解析。[SConstruct](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/SConstruct) 指定 C++20 并提供依赖配置；预计 **medium**，主要成本是取得与 SHA 匹配的开发版库和完整测试数据，而非数值计算量。

  以合并源码为准：PR 正文声称实现 HKFT 的 `cv_mole`，但该 SHA 的 `PDSS_HKFT.h/.cpp` 无此覆盖；不要把直接调用 `PDSS_HKFT::cv_mole()` 作为已实现的检查入口，本文使用含 HKFT 物种的 `ThermoPhase::cv_mole()`。本次仅阅读 `input/pr.json`、`input/pr.diff` 和 GitHub 上指定 SHA 的公开文件，并写出本文；没有克隆仓库、安装、编译或执行数值代码，没有读取其他审计轨道或工作目录外文件。

- Formalization (Lean) feasibility: can the invariant be written as a proposition over finite-dimensional reals, and what are the main obstacles:

  **可以。** 用 \((T,P,n_0,\ldots,n_{N-1})\in\mathbb R^{N+2}\)、\(n_k>0\)、\(x_k=n_k/\sum_jn_j\) 表示状态，定义 \(G=(\sum n_k)g\)，则偏摩尔量和 I2 的 Maxwell 关系是有限维实函数的偏导命题。Margules 的特例基本是多项式/有理函数求导及一次齐次性，SSVol 压力热容项是二阶导数与单位换算，I3 是正实数上的代数和对数等式，适合首先形式化。

  I1 需明确 \(V_P\ne0\) 并用隐函数/链式法则构造恒体积路径，再在 \(V_P<0\) 下证明热容差非负；\(V_P=0,V_T\ne0\) 的分支应作为不满足前提的模型退化处理，不能沿用除法恒等式。主要障碍是完整 HKFT/IAPWS 参数化与迭代求密度的数学对应、截断导致的分段光滑性，以及从实数模型证明连接到 C++ 浮点执行、缓存状态恢复和内部有限差分的误差界；抽象 Gibbs 恒等式的证明不能自动证明生产实现满足它。
