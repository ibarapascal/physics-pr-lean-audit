# ct2107 实施日志（ΔQ1 统一抽查，Cantera#2107）

> 目的：三路检查计划逐句→代码对应与全部实施判断，供协议 §6 核查与 §7 偏离登记。时间切片：2026-10-17 冒烟前写成；冒烟后只追加「冒烟后修正」节。
> 文件：`spotcheck.py`、`harness.py`（四份 PR 目录相同副本）、`job_smoke.json`、`job_full.json`。
> 读过的材料：协议；`track1.md`/`track2.md` 本 PR 节；`specs_codex/Cantera_cantera_2107/spec.md` §2/§4/§5；合并提交 d58a2326 源码树 `test/data/*.yaml`（含 consistency-cases.yaml）、`MargulesVPSSTP.h/.cpp`、`BinarySolutionTabulatedThermo.h/.cpp`；`env/build_cantera.sh`。

## 0. 环境与通用约定

- 原生环境：compute node A `env/ctb`，由 `env/build_cantera.sh` 从合并提交 d58a2326 源码 scons 构建（`ct.__version__`=4.0.0a1，`__git_commit__`=d58a232）。**有原生环境**（合并提交源码构建，非发布版）。相文件读同一源码树的 `env/cantera-src/test/data`。
- 只用 Python 接口（`ct.ThermoPhase`）；第三路计划要求的 C++ harness 未写（见 §5 D1）。
- 种子＝20261017＋1×10＋路号；分支＝相（配置），N 等分；判定、重放、失败计数规则同 pp3016 日志 §0。
- 质量摩尔浓度相的组成：由 m_k 换算 X（X_o=1/(1+M̃_oΣm)，X_k=m_kM̃_oX_o，M̃_o=分子量/1000 kg/mol）后 `TPX` 设定；固定 X 即固定 m。
- 计划只给相对容差时绝对项取 0；完全没给时用相对 1e-6、绝对 1e-12。
- 「含 IAPWS 的相」判定：相内任一物种 equation-of-state 为 liquid-water-IAPWS95 或 HKFT（运行时读 input_data）。
- YAML 改写（合成夹具、截断模型、标准浓度基准、物种顺序）用 ruamel.yaml 读源码树文件、改字段、`ct.ThermoPhase(yaml=…)`。

## 1. 第一路（track1.md）

相选择（判断 C1）：计划列了 7 个文件、没列相名 → 取这些文件里热力学模型属于本 PR 改动范围（Margules、ideal-molal-solution、Debye-Huckel、HMW-electrolyte、ideal-solution-VPSS，含 HKFT/SSVol 物种）的相，排除 `HMW-bogus`；共 23 个分支：debye-huckel-all 的 10 个相、thermo-models 的 debye-huckel-B-dot-ak/beta_ij、HMW-NaCl-electrolyte、HMW-NaCl-HKFT、ideal-molal-aqueous、ethanol-water-rough、molten-salt-Margules、IdealSolnGas-liquid，以及 HMW_NaCl、pdss_hkft、LiKCl_liquid、Li_Liquid、VpssSolidSolutionExample。

| 计划原句 | 代码 |
|---|---|
| 水溶液 T∈[273.15,373.15] K（HKFT 可到 573 K）、p∈[1,100] atm、I∈[1e-4,6] mol/kg（HMW 可到 ~6.1）；熔盐 Margules T∈[700,1400]；固溶体 VPSS T∈[300,1000] | `t1_domain`、`t1_sampler` |
| A：读 cv、cp、β、κ_T、v；β、κ_T 与 v(T±h,p)、v(T,p±δ) 中心差分比（h=1e-4·T，δ=1e-4·p；含 IAPWS 的相 rtol 4e-6，其它 1e-7） | `t1_eval` 的 `A.beta=…`、`A.kappaT=…` |
| A：c_v = c_p − T v β²/κ_T（κ_T>0），κ_T=0 时 c_v=c_p（库约定） | `A.Mayer…`（同 A 的 rtol）、`A.kappa=0 => cv=cp`（默认容差） |
| A：等体积路径 (T,p)→(T+h, p+(β/κ_T)h)，牛顿 2–3 步把 v 拉回，c_v≈T·Δs/ΔT，rtol 1e-5 | `isochoric_newton` 3 步；正负两侧各做一次取中心差分（判断 C3）；残差进 records |
| B：chemical_potentials 与 standard_gibbs_RT·RT＋RT·log(γX) 比，rtol 1e-10 | `B.mu=mu0+RTln(gX)` |
| B：溶质项再与手算 m_k=X_k/(X_oM̃_o) 的 γ^□m_k 比（γ^□ 若 Python 无接口用 C++ getMolalityActivityCoefficients） | Python 无 γ^□ 接口；取 `activities`（MolalityVPSSTP 的质量摩尔浓度基活度 a_k=γ^□m_k/m°）作 γ^□m_k，与 γ_kX_k 比，默认容差（判断 C4） |
| C：μ_k(T±h) 中心差分对 −partial_molar_entropies（rtol 1e-6）；ΣX_k s̄_k 对 entropy_mole（rtol 1e-12）；h(T±h) 中心差分对 cp_mole；s(p±δ) 与 −v(T±h) 互比 | `C.*`；后两项计划无容差→默认 |

实施判断（第一路）：
- C1 相的选取见上。
- C2 组成：计划只给 I 的区间 → 1:1 的 Na⁺=Cl⁻ 取 m=I（均匀抽），其余物种固定在 consistency-cases.yaml 状态点 0 的值（DH：H⁺ 1.05e-8、OH⁻ 1.3765e-6、NaCl(aq) 0.98492…；HMW：H⁺ 1.05e-5、OH⁻ 1e-5）；ideal-molal、ethanol-water、熔盐、固溶体、HKFT（摩尔分数相）组成固定在状态点 0。熔盐与固溶体计划没给压力 → 固定 1 atm。ethanol-water 归入「水溶液」类区间（计划没单列）。HKFT 类（pdss_hkft、HMW-NaCl-HKFT）用 T 到 573.15 K。
- C3 等体积 c_v：单侧 T·Δs/ΔT 的截断误差约 h/(2T)·(T s_TT/s_T) ~5e-5，超过计划的 1e-5；按「中心差分」两侧各走一次等容路径。
- C4 γ^□m_k 用 `activities`（见上表）；这比 C++ 接口少一层独立性，记为偏离候选。
- C5 κ_T、β 的比较绝对项取 0（计划只给相对）。p 低于饱和压力时（如 573 K、1 atm）水的 PDSS 可能给出气相分支或报错——照算，报错计失败。

## 2. 第二路（track2.md）

分支（17 个，取 consistency-cases.yaml 中本 PR 改动模型的用例相）：LiKCl_liquid、ethanol-water-rough、Li_Liquid、debye-huckel-all 10 个相、HMW_NaCl、ideal-molal-aqueous、VpssSolidSolutionExample、pdss_hkft。

| 计划原句 | 代码 |
|---|---|
| 常用工况域：LiKCl T 900–1400、P 1–20 atm；乙醇–水 T 298–308、X_eth 0.3–0.5；Li(L) T 300–500、P 1–10 atm；DH/HMW T 300–360、P 1–5 atm、m 0.01–9.4；HKFT T 300–360；「建议把压力扩到 1–100 atm」 | `t2_sampler`：T、X_eth、m 按表；P 一律 1–100 atm（判断 D1） |
| I1：TPX 设 P±dP，chemical_potentials 中心差分与 partial_molar_volumes 比；dP=max(1e-3P,1e3 Pa)；rtol 1e-6、atol 1e-12·max\|v̄\|，IAPWS 相 1e-5 | `I1.dmu/dP=vbar`（适用于计划列的相，即全部 17 个） |
| I2：二阶中心差分 (v(T+h)−2v+v(T−h))/h²，h=0.5–1 K；cp 对 P 一阶中心差分；rtol 1e-4 | h 取 0.5 K（判断 D2）；`I2.dcp/dP=-T d2v/dT2` |
| I2：SSVol 物种 c_p,k°(T,P)−c_p,k°(T,P_ref)=−2c·T·(P−P_ref) | 只对 molar-volume-temperature-polynomial 物种（Li(L)2）；数据是三次多项式，用 v°''=2c+6dT（判断 D3）；P_ref=1 atm |
| I3：cv、cp、v、β、κ_T 代数比较 rtol 1e-10；κ_T=0 时 c_v=c_p；c_v≤c_p | `I3.*` |
| I3：activities 与 activity_coefficients·X，rtol 1e-12（MolalityVPSSTP 子类） | DH、HMW、ideal-molal 分支 |
| I3：Gibbs–Duhem：组成微扰 δX（ΣX=1，\|δX\|~1e-6）算 ΣX_kΔμ_k，相对 max\|Δμ_k\| <1e-6 | δX_k=1e-6·X_k·(r_k−Σ_jX_jr_j)（r 每点随机）；Δμ 用 X±δX 的中心差（判断 D4） |

实施判断（第二路）：
- D1 压力：计划域里的「建议扩到 1–100 atm」按协议「建议也算」→ 所有相 P 均匀 1–100 atm（含 ethanol-water、HKFT、ideal-molal）。表里没列 T 的相：ideal-molal 用其状态点范围 300–340 K，VpssSolid 用 300–500 K。
- D2 二阶差分步长取区间下端 0.5 K。
- D3 计划把 SSVol 写成二次多项式（v''=2c），实际数据是 4 项三次多项式 → 用该多项式的精确二阶导；density-temperature-polynomial 物种（Li(L)）不做显式式，只走一般 I2。
- D4 Gibbs–Duhem：按 X_k 比例微扰保证正值且 ΣδX=0；单侧差分的二阶项相对量 ~1e-6·X，会与 1e-6 的判据同量级 → 改用中心差分消去二阶项（这是对计划「Δμ」的实施判断）。

## 3. 第三路（spec.md §5）

分支（26 个）：debye-huckel-all 10 相、HMW_NaCl、ideal-molal 截断模型 none/poly/polyexp（3，polyexp 即原文件）、ideal-molal 标准浓度 unity/species-molar-volume（2，solvent-molar-volume 即原文件）、pdss_hkft、VpssSolidSolutionExample、Li_Liquid、LiKCl_liquid、ethanol-water-rough、margules_synthetic（乙醇–水交互参数随机化）、Redlich-Kister-LiC6、Redlich-Kister-complex、anode（表格二元）、tab_simplified（零背景、unity 标准浓度的简化表格夹具）。必测点：DH/HMW 的 stress 质量摩尔浓度 {0.1,1,6,9.4}（298.15 K、1 atm）；两种表格夹具的全部节点、相邻中点、x=0、1、{1e-2,1e-6,1e-10,1e-14,1e-300}，以及显式给参考物种痕量 {1e-2,…,1e-300} 的对称痕量态，两种物种排列顺序各一遍。

| 计划原句 | 代码 |
|---|---|
| 水/HKFT/DH/HMW：T 280–360 K、P 0.101325–10 MPa；NaCl 1e-6–1e-2 mol/kg（0.1,1,6,9.4 只作 stress）；组成内部 1e-6–1−1e-6；Margules/SSVol：乙醇–水 280–330、Li 300–500、熔盐 900–1400、P 1–20 atm；表格 280–350 K、1–10 atm、x∈[0,1] | `t3_sampler` |
| 中心差分 h_T/T∈{1e-4,5e-5,2.5e-5}、h_P/P∈{1e-3,5e-4,2.5e-4}，先比相邻估计；不收敛＝不确定 | `deriv3`：后两档差超过容差记失败，否则用最细档比 |
| 代数恒等式 \|a−b\| ≤ 1e-5＋1e-12·max；导数相对 1e-6（解析）/1e-5（水/HKFT），绝对下限 cp 1e-3、体积 1e-10、Maxwell 1e-9、κ 5e-14、β 1e-9 | `alg`、`alg_v`、`deriv3`、`FLOOR` |
| I1：g=h−Ts、h=u+Pv、μ_k=h̄_k−Ts̄_k、f=ΣX_kf̄_k（v,h,s,g,c_p）、c_p=h_T=Ts_T、s_P=−v_T、v̄_k=∂_Pμ_k | `t3_I1` |
| I1：SSVol 两种多项式、非零曲率、在 P_r 与远离 P_r，与解析 SSVol 式比 | `t3_ssvol`（Li_Liquid 两个物种各一种多项式；P 取抽样值，P_r=1 atm） |
| I1：setTemperature 后 setPressure、反序、setState_TP 结果一致 | `t3_setorder`：Python 用 `TP=(T,None)`/`TP=(None,P)` 近似（判断 E2） |
| I1：ideal molal 的 setCutoffModel none/poly/polyexp；溶剂富集态与过渡两侧 | 改 YAML `cutoff.model` 三个分支；X_o 在 [1e-6,1−1e-6] 均匀抽 |
| I1：GibbsExcess 子类含 Margules（有/无超额体积）与 Redlich–Kister | 对应分支 |
| I2：五种 PDSS（定容、两种 SSVol、IAPWS 水、HKFT）；dVdT/dVdP 对外部差分的 molarVolume；导数调用前后父相/水的性质不变 | 相级：β、κ_T 对差分体积（`deriv3`）；反复调用导数后 h、s、v、μ 不变（`==`）；PDSS 级 dVdT/dVdP 不可达（D1） |
| I2：比较 ρ、偏摩尔体积和、β、κ_T、c_p−c_v 与独立差分体积 | `I2.v=Mbar/rho`、I1 的 v=ΣXv̄、`I2.cp-cv=-T vT^2/vP` |
| I2：v_P≠0 时用括区间标量根在 P 上解 v(T±h,P±)=v(T,P)，c_v 与 u、s 的等容导数比；根残差 <1e-11 相对体积；括不住＝不确定 | `_root_P`＋`t3_I2`（两档步长收敛再比） |
| I2：v_P=0、v_T≠0 的回退只作约定检查；HKFT 直接 cv_mole 记为不支持 | 约定只记录；pdss_hkft 的相级 cv_mole 可调用（PDSS 级不可达），出错才记录 |
| I3：a=Xγ（1e-10 相对＋1e-20 绝对）；化学势差（RT 单位）1e-9 绝对；ideal molal、HMW、10 个 DH；标准浓度 unity/species-molar-volume/solvent-molar-volume | `t3_I3`；标准浓度只对 ideal-molal 改 YAML（DH/HMW 无此选项） |
| I4：g^E,h^E,s^E,v^E 由独立多项式算，与相性质减标准态/理想混合比；P=P*、2P*、20P*、10 MPa；各系数含零值与正负号 | `margules_excess`（式取 MargulesVPSSTP.h 文档，参数从 input_data 读 SI 值）、`t3_I4`；系数随机化在 margules_synthetic 分支（每个系数 1/4 概率为 0，否则在 ±(数据量级) 均匀） |
| I4：归一化组成路径 Σ_kX_kD_k=0；J_kj 由 n_j↦n_je^{±ε} 归一后差分 ln γ（ε=1e-4,1e-5,1e-6） | `I4.Gibbs-Duhem…`、`I4.sum_k X_k dlngamma_k/dln n_j=0`（容差 1e-6 相对＋1e-8 绝对） |
| I5：DH 溶剂与各溶质偏摩尔 c_p、体积与化学势/焓导数比 | 体积在 I1（v̄=∂_Pμ），c_p 在 `t3_I5`（c̄_p=∂_T h̄） |
| I6：anode 与简化夹具；节点、中点、端点、痕量；两种排列；有限性与 I1 代数和；简化夹具上 μ_tab−μ_ref=h_tab−Ts_tab（独立插值）与 s=x s_tab−R ln(1−x)；换温度后回到同组成检查缓存 | `t3_I6`；插值用自写线性插值（表外取端值，同源码约定） |

## 4. 未实现项（第三路，Python 接口不可达；均为偏离候选）

- `MargulesVPSSTP::getdlnActCoeffdT / getd2lnActCoeffdT2 / getdlnActCoeffds / getdlnActCoeffdlnN / _diag / dlnX_diag` 与有限差分的对照（I4）。
- `WaterProps::ADebye(T,P,ifunc)` 四种模式的导数对照、饱和压力钳位约定（I5）。
- `providePDSS(k)` 的 dVdT/dVdP 直接对照、PDSS_HKFT 的直接 cv（I2）。
- `setTemperature` / `setPressure` 单独调用（用 TP 的 None 形式近似）、`setMolalities`（用 X 换算）、`getMolalityActivityCoefficients`/`getMolalities`（手算 m，γ^□ 只能由 a/m 得到，等于同义反复，不作检查）。
- I2 的「water+rigid」「water+HKFT 进 Margules」等计划要求而测试数据里没有的合成混合物；I4 的三元合成相（零/一/两个交互）未构造。
- I6 的非零背景/不等标准浓度/升压下「显式修正」的单独核对未做（计划说只标记量纲问题）。

## 5. 耗时与内存估计（推测）

单次 TPX＋性质：非水相 ~0.02–0.1 ms，IAPWS 相 ~0.2 ms（经验值）。第一路每点 ~40 个状态 → 5–20 ms；第二路 ~30 个状态 → 同量级；第三路每点 ~200–1000 个状态（等容求根最多 ~800 次）→ 0.05–0.3 s。

| 路 | N=100（含必测） | N=10⁴ |
|---|---|---|
| 1 | ~10 s | ~5 min |
| 2 | ~10 s | ~5 min |
| 3 | ~2 min（必测 ~350 个表格点很快） | ~30 min（上限估 1 h） |

峰值内存 ~300 MB（Cantera＋numpy；margules_synthetic 每点新建相，旧相随引用释放）。JobSpec 申报 1 GiB，全量串行。

## 6. 偏离候选（交主会话裁定）

- D1 第三路要求 C++ harness，本脚本只用 Python，§4 所列项未测。若要补，需在 compute node A 用 `env/cantera-src/build/lib/libcantera_shared.dylib`＋源码 include 编一个小 harness（需另行授权构建）。
- D2 第一路 B 的 γ^□m_k 用 `activities` 代替 C++ 接口（C4）。
- D3 第二路 Gibbs–Duhem 用中心差分（D4）。
