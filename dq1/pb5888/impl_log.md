# PyBaMM#5888 ΔQ1 抽查实施日志

> 目的：记录三路检查计划逐句到 `spotcheck.py` 的映射与全部实施判断，供协议 §6/§7 审核。
> 结构：§0 环境 → §1 第一路 → §2 第二路 → §3 第三路 → §4 耗时与内存 → §5 偏离候选 → §6 冒烟后修正 → §7 待裁定（冒烟拒绝的归类）。
> 更新规则：冒烟后只许修脚本错误；改计划映射须记入协议 §7，并在本文件 §5 追加一行。
> 引用：上游 协议 `notes/20261017-prereg-dQ1-uniform-spotcheck.md`、track1.md/track2.md 的 #5888 节、`specs_codex/pybamm-team_PyBaMM_5888/spec.md` §2/4/5；下游 `dq1/out/pb5888_t{1,2,3}.json`。

## §0 环境与通用约定

- 节点 compute node A，`env/py2`（`setup_pca_py.sh` 建，pybamm 26.10.1.0 来自 PyPI，不是从合并提交构建）。核对：`Ai2020.py`、`Xu2019.py`、`lithium_ion_parameters.py`、`settings.py` 与合并提交 343ece2 的 GitHub raw 逐字节相同（`diff -q` 无输出）。实施判断：视为合并提交原生环境。
- 种子 20261017＋11×10＋路号（t1=20261127，t2=20261128，t3=20261129）。每路只有 1 个分支（被改函数无选项/相/流态分支），N 点全部归这一支。
- 计划给出的网格、离散值、名义值另加为必测点（mandatory，不占 N），结果单独计数。
- 判定、重放与 JSON 格式见 `spotlib.py` 文件头：任一检查有限且超容差 → 拒绝候选 → 同输入重算一次逐位相同才计拒绝；异常或 NaN 计失败。
- 取值 API：`ParameterValues(...).evaluate(pv[key](pybamm.Scalar(c), pybamm.Scalar(T)))`（三路一致）。为防 10⁴ 点下内部表达式缓存无界增长，每 2000 次求值调用一次 `_processor.clear_cache()`（私有属性，hasattr 保护；不影响数值）。
- 实施者阅读范围说明：第三路 spec.md 曾一次性整文件输出（含 §1、§6–§9）；实现只依据 §2、§4、§5。

## §1 第一路（track1.md「数值核对方法」）

| 计划原文 | 代码 |
|---|---|
| `pv = ParameterValues("Ai2020")`；`f = pv["Electrolyte diffusivity [m2.s-1]"]`；`pv.evaluate(f(Scalar(c), Scalar(T)))` | `D_ai()` |
| 常用工况域 c_e∈[100,3000]、T∈[253,333] | `t1_sample`：两变量线性均匀 |
| 对 (c_e,T) 网格（c 100…3000 步 100，T 253…333 步 5） | 必测点 `grid_c=…_T=…`，30×17=510 点，跑与抽样点相同的全部检查 |
| 与 numpy 直接按文献式（含 1e-4）比，rtol 1e-12 | `t1_api_vs_literature_rtol1e-12`（`t1_lit`，atol=0） |
| 与 `ParameterValues("Xu2019")` 同名参数比，rtol 1e-12 | `t1_ai_vs_xu2019_rtol1e-12` |
| 量级界 1e-11 ≤ D ≤ 2e-9（不变量「次」） | `t1_magnitude_1e-11<=D<=2e-9` |
| 有限差分符号 ΔT=0.01 K、Δc=0.1 mol/m³，中心差分，只看符号；∂D/∂T>0，∂D/∂c<0 | `t1_fd_sign_dD_dT>0`、`t1_fd_sign_dD_dc<0`（用 API 值做差分） |
| 可选：1C DFN 放电（Ai2020）末态 c_e 明显梯度 max−min ≫ 1 mol/m³ | 必测点 `optional_dfn_1C`：`lithium_ion.DFN()`＋Ai2020 默认参数（Current function 2.28 A＝1C），`sim.solve([0,3600])`，取全电池 `Electrolyte concentration [mol.m-3]` 末时刻 |

实施判断：
1. 「可选」检查按协议 §2「计划写到的都算」纳入；它没有可抽样的输入，只跑一次。
2. 「≫ 1 mol/m³」取阈值 ≥10 mol/m³（一个数量级）。
3. 网格点与抽样点不去重。

## §2 第二路（track2.md「数值核对方法」）

| 计划原文 | 代码 |
|---|---|
| 常用工况域 c_e 0–3000、T 263–333、DFN 放电 C/10–2C | `t2_sample`：c_e、T、C_rate 三者线性均匀 |
| I1：与按原式独立用 numpy 算的值比 rtol 1e-12（原式按 mol/L：c=c_e/1000） | `t2_I1_api_vs_formula_rtol1e-12`（`t2_formula`） |
| I1：量级界 1e-11 ≤ D_e ≤ 1e-9 | `t2_I1_magnitude_1e-11<=D<=1e-9` |
| I1：1000 mol/m³、298.15 K 处与 PR 的 3.2227e-10 一致 | 必测点 `nominal_1000_298.15`：abs 容差 0.00005e-10（5 位有效数字的半个末位） |
| I2：与 Xu2019 同名参数比 rtol 1e-12 | `t2_I2_ai_vs_xu2019_rtol1e-12` |
| I2：单调性中心差分 ΔT=0.01、Δc=0.1「在 100 个随机点检符号」 | `t2_I2_fd_sign_*`，在全部 N 个抽样点上做 |
| I3：`DFN()`（Ai2020），1C 放电 3600 s；取隔膜网格上 c_e 与温度；二阶差分算 ∂c/∂x；与 (1−t⁺)i/(F·D_e,eff) 比，容差 5% | `DFNCase`＋`t2_dfn_checks`：`np.gradient(c, x_s, edge_order=2)`；D_e,eff = D_e(c,T)·ε_sep^b（`Separator porosity`=0.5、`Separator Bruggeman coefficient (electrolyte)`=1.5）；t⁺、F、i 取参数集/`pybamm.constants.F`/`Total current density [A.m-2]`；逐隔膜节点 rtol 0.05 |
| I3：准稳态（t ≳ 1000 s 后） | 若终止时刻 <1000 s，该条记「不适用」（记通过） |
| I3：∫εc_e dx 相对初值漂移 ≤1e-6 | `t2_I3_electrolyte_Li_drift<=1e-6`：用库变量 `Total lithium in electrolyte [mol]` 末/初 |
| 1C 这一离散值 | 必测点 `C_rate=1_c=1000_T=298.15` |

实施判断：
1. 「100 个随机点」被协议 §3 的 N 点均匀抽样取代（同一检查在全部 N 点上做）。
2. C 倍率从工况域 [0.1, 2] 均匀抽，1C 作必测点；放电时长取 3600/C（把「1C 放电 3600 s」按名义全放电时长推广），到 Ai2020 下截止 3.0 V 自动终止；在解的最后时刻取值。
3. 为复用已离散的模型，`Current function [A]` 设为 `[input]`，`sim.solve([0, 3600/C], inputs=...)`；电流 = C×`Nominal cell capacity [A.h]`(2.28)。
4. D_e,eff 中的 D_e 用模型侧 `LithiumIonParameters.D_e`（带 10 mol/m³ 下限；隔膜内 c_e 远高于下限，不起作用）；温度取 `Separator temperature [K]`，取不到用环境温度。
5. 符号约定：放电电流为正，−D_eff ∂c/∂x 与 (1−t⁺)i/F 按带符号比较。
6. 重放：被拒检查只涉及本构项时，重放重算本构、复用同点刚算的 DFN 结果；涉及 `t2_I3_*` 时整点（含 DFN）重算（`t2_replay`）。
7. c_e 接近 0 时 Δc=0.1 的中心差分会用到 −0.1 mol/m³，原式仍有定义，照算。

## §3 第三路（spec.md §2/§4/§5）

| 计划原文 | 代码 |
|---|---|
| 工作窗 500≤c≤2000、273.15≤T≤333.15 | `t3_sample` 线性均匀 |
| Cartesian 7×5 网格 | 必测点 `grid_c=…_T=…`（35 点，跑 `t3_eval` 全部检查） |
| 200 个 `default_rng(5888)` 均匀对 | 被协议 N 点、协议种子取代 |
| I1：与 C=c/1000、T_g=229+5C、`np.exp(np.log(10)*(…))*1e-4` 比，rtol 1e-12、atol 0，显式有限 | `t3_I1_api_vs_oracle_rtol1e-12`（`t3_oracle_D`；`SL.close` 对非有限值记失败） |
| I1：名义值 3.222722528605257e-10 同容差 | 必测点 `nominal` |
| I1：若另装基线则比 new/base=1e-4 | 未装基线（条件不满足），不做 |
| I2：Ai2020 与 Xu2019 每个样本 rtol 1e-12 | `t3_I2_ai_vs_xu2019_rtol1e-12` |
| I2：模型侧 `lp.D_e` 与施加浓度下限后的原函数比 | `t3_I2_model_De_vs_raw_floor_rtol1e-12`（`LithiumIonParameters().D_e`） |
| I2：c=0,5,9.999,10,10.001,20 × 每个列出温度；读设置并要求 =10 | 必测点 `boundary_c=…_T=…`（30 点），`t3_I2_setting_D_e__c_e==10` |
| I3：raw 与模型侧有限且 >0（所有正常点与边界点） | `t3_I3_raw_finite_positive`、`t3_I3_model_finite_positive` |
| I4：网格每行随 T 严格增、每列随 c 严格减 | 必测点 `grid_monotone` |
| I4：随机内点 c∈[501,1999]、T∈[273.16,333.14] 上对 ln D 中心差分 h_c=0.1、h_T=0.001，与解析对数导数比 rtol 1e-6 | 抽样点落在该内区时做 `t3_I4_dlnD_dT/dc_vs_analytic_rtol1e-6` |
| I4：包装函数下限以下为常数（rtol 1e-12）、跨下限连续（由公式对照覆盖）、下限以上严格减；不在拐点求导 | 边界点 `t3_I4_wrapper_constant_below_floor`、`t3_I4_wrapper_strict_decrease_above_floor`（与前一个边界浓度比） |
| I5：电导率与 §1 多项式比 rtol 1e-12；有限非负；raw c=0 恰为 0；名义 1.1943 rel 1e-3；`kappa_e` 带下限对照 | `t3_I5_*`（`t3_oracle_kappa`） |
| 网格上成对 Vector、标量/Vector 各组合逐元素与标量比 rtol 1e-12 | 必测点 `vector_grid`（扩散率与电导率两个函数都做） |

实施判断：
1. 电导率多项式的「§1 转写」：实施者未依据 §1，直接按合并提交 `electrolyte_conductivity_Ai2020` 源码与 §2 I5 的 κ=1e-4·c·P² 写参考式。
2. 抽样点都在下限以上，下限分支只由必测边界点覆盖（计划 §5 说这是唯一相关分支）。
3. 「跨下限连续」没有单独容差，由各边界点「模型侧＝原函数(max(c,10))」的 rtol 1e-12 对照代表。

## §4 预计耗时与峰值内存（推测，未运行；冒烟后更新）

| 路 | 单点 | 10⁴ 点 | 依据 |
|---|---|---|---|
| t1 | ~15 ms（7 次 evaluate） | ~3 min＋510 网格＋1 次 DFN（~10 s） | evaluate 单次约 1–3 ms（表达式处理＋求值） |
| t2 | ~0.5–1 s（1 次 DFN 求解＋本构） | ~1.5–3 h（约 3 成点因量级界被拒而重放本构，DFN 复用） | Ai2020 DFN 默认网格 IDAKLU 求解约 0.3–0.8 s，模型只离散一次 |
| t3 | ~25 ms | ~5 min | 约 10 次 evaluate |

峰值内存：约 1–1.5 GB（pybamm/casadi 导入约 0.4 GB，DFN 离散模型与解约 0.5 GB）；JobSpec 申报 3 GiB。
若冒烟显示 t2 单点 >1.4 s（10⁴ 点超 4 h），按协议 §3 把 t2 的 N 降到 4 h 能跑完的最大值并记 §7。

## §5 偏离候选

- （无）环境为 PyPI 发行版而非从合并提交构建：相关源文件逐字节相同，记作说明，不作偏离。

## §6 冒烟后修正（2026-10-17）

- 无。首轮拒绝均为计划本身的量级界问题（§7），脚本忠实于计划，未改。

## §7 待裁定（冒烟拒绝归类，脚本不改）

1. **(P) 第一路「次」不变量量级界 1e-11 ≤ D ≤ 2e-9（「在域内」）在计划自己的工况域 c_e∈[100,3000]、T∈[253,333] 上不成立**。首轮 77 次拒绝（抽样 14/100、网格必测 63/510）全部是 `t1_magnitude` 下界；同点 API＝文献式（rtol 1e-12，0 拒绝）、Ai2020＝Xu2019（0 拒绝）、导数符号全过——库的值就是计划给的式子。拒绝点集中在低温高浓度角：T ≲ 272 K 且 c_e ≳ 440 mol/m³，越近奇点 T = 229 + 0.005c_e 越小（如 c_e=2855、T=255.2 → 2.5e-14）。按计划原式手算：c_e=1000、T=253 → 3.2e-12；c_e=3000、T=253 → 8e-16，均低于 1e-11。对矩形工况域均匀抽样，下界被违反的面积占比约 10.4%（10⁶ 点蒙特卡罗，按原式），全跑预计约 1000 次同类拒绝；上界 2e-9 在本域内不触发（最大约 1.06e-9）。计划写错之处：「在域内 1e-11 ≤ D ≤ 2e-9」与同节给的 Valøen–Reimers 式及「须避开奇点」的工况域说明自相矛盾（第一性原理：式子在奇点附近 D→0）。
2. **(P) 第二路 I1 量级界 1e-11 ≤ D_e ≤ 1e-9（「在上述域内」，c_e∈[0,3000]、T∈[263,333]）同样不成立**。首轮 3 次拒绝（抽样）全为下界：c_e≈2820–2955、T≈266–271 → 3.8e-12～8.3e-12；I1 公式对照、Xu2019 对照、导数符号、DFN 隔膜通量与锂守恒全过。按原式，下界违反区 T ≲ 272 K 且 c_e ≳ 1830（占比约 2.5%），**上界 1e-9 在 c_e ≲ 200、T ≳ 324 K 也被违反**（如 c_e=0、T=333 → 1.12e-9，占比约 0.5%，首轮 100 点未抽到，全跑会出现）。计划写错之处：I1「在上述域内 1e-11 ≤ D_e ≤ 1e-9」与同句的原式在域角不符。
3. 本 PR 首轮未见 (S)、(V) 类。
