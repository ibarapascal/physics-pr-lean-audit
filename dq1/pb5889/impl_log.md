# PyBaMM#5889 ΔQ1 抽查实施日志

> 目的：记录三路检查计划逐句到 `spotcheck.py` 的映射与全部实施判断，供协议 §6/§7 审核。
> 结构：§0 环境 → §1 第一路 → §2 第二路 → §3 第三路 → §4 耗时与内存 → §5 偏离候选 → §6 冒烟后修正 → §7 待裁定（冒烟拒绝的归类）。
> 更新规则：冒烟后只许修脚本错误；改计划映射须记入协议 §7，并在 §5 追加。
> 引用：上游 协议、track1.md/track2.md 的 #5889 节、`specs_codex/pybamm-team_PyBaMM_5889/spec.md` §2/4/5；下游 `dq1/out/pb5889_t{1,2,3}.json`。

## §0 环境与通用约定

- 节点 compute node A，`env/py2`（pybamm 26.10.1.0，PyPI）。核对：`Chayambuka2022.py`、`data/D_p.csv`、`data/D_n.csv`、`lithium_ion_parameters.py`、`settings.py`、`expression_tree/interpolant.py` 与合并提交 eda03ed 逐字节相同。视为合并提交原生环境。
- 父提交 c235b90 的 `Chayambuka2022.py` 与合并提交只差 162 行一处（c_max 取 `Initial concentration in electrolyte`）；第二路「与修复前版本比」用这一处差异在同一环境里重建父提交函数（`base_NVPF_diffusivity`，参数覆盖），不另装父提交。
- 种子 20261017＋12×10＋路号。每路 1 个分支（单一 NVPF 相、单一线性插值、无选项）。必测点不占 N。
- 数据表由脚本独立读取：`np.loadtxt(D_p.csv, delimiter=",", skiprows=1)`（安装包内文件，已核与合并提交相同）。
- 每个 (C_max, E) 配置各建一个全新 `ParameterValues`，建好后只 `update` 一次、不再改动；每 2000 次求值清一次内部缓存。

## §1 第一路

| 计划原文 | 代码 |
|---|---|
| 常用工况域 sto∈[0.01,0.99]、T∈[273,323]（函数不依赖 T） | `t1_sample` 线性均匀 |
| sto 网格 0.01…0.99 步 0.01、T=298.15 | 必测点 `grid_sto=…`（99 点，跑全部检查） |
| `pv.evaluate(f(Scalar(sto), Scalar(T)))` 与 `numpy.interp(sto*15320, x_csv, y_csv)` 比，rtol 1e-10（默认线性） | `t1_api_vs_numpy_interp_rtol1e-10` |
| `pv.update({c_max: c})` 改成别的 c，固定 c_s 下 f(c_s/c) 不变，rtol 1e-10 | `t1_cmax_invariance_c=…`：c ∈ {7660, 30640} |
| 次（量级）：D_p ∈ [1.9e-17, 2.6e-15]；sto≥0.3 时 D_p ≤ 1e-16 | `t1_magnitude_*` |
| 宏观旁证：BasicDFN 默认参数 2C 放电至 2.0 V 容量约 1.18 mAh | 必测点 `macro_BasicDFN_2C`：`Experiment(["Discharge at 2C until 2.0 V"])`，读 `Discharge capacity [A.h]` 末值 |

实施判断：
1. 「别的 c」计划未给值：取 0.5×与 2×默认 c_max（7660、30640），固定 c_s = sto×15320；查表自变量 c_s/c×c 还原为 c_s，c_s 落在表内。
2. 「约 1.18 mAh」按 3 位有效数字的半个末位取 abs 容差 0.005 mAh。计划自注「非物理不变量，只作回归参考」，仍按协议 §2 纳入。
3. 默认插值器已核为 linear（`Interpolant(..., interpolator="linear", extrapolate=True)`），不需要计划里的替代比较。

## §2 第二路

| 计划原文 | 代码 |
|---|---|
| 常用工况域：sto 0.01–0.99，T 298.15 固定，倍率 C/50–3C，截止 2.0 V | `t2_sample`：sto、C_rate 线性均匀，T 固定 |
| 工况域「BasicDFN 与 DFN/SPMe 模型」 | 只用 BasicDFN：`pybamm.sodium_ion` 只提供 `BasicDFN`（`__init__.py` 只导出它），没有钠离子 DFN/SPMe 可用 |
| I1：表的每个节点 c_k，sto_k=c_k/15320，与表值 y_k 比 rtol 1e-10 | 必测点 `node_k=0…33` |
| I1：节点间随机点与 `numpy.interp` 比 rtol 1e-10 | 抽样点 `t2_I1_api_vs_numpy_interp_rtol1e-10` |
| I1：sto∈{0,1} 检正值 | 必测点 `endpoint_sto=0/1` |
| I1（不变量陈述）：sto∈[0,1] 上 1.7e-17 ≤ D_p ≤ 2.51e-15，且 D_p>0 | 抽样点与两个端点都检 `t2_I1_bounds_*`、`t2_I1_positive` |
| I1（不变量陈述）：负极同理 D_n(sto)=L_n(sto·14540) | 抽样点 `t2_I1_negative_Dn_vs_numpy_interp_rtol1e-10`（`D_n.csv`） |
| I2：BasicDFN＋`Experiment([f"Discharge at {C} until 2.0 V"])`，读 `Discharge capacity [A.h]` 末值；C/50 与修复前版本比 ≤2%；2C 与 1.18 mAh 比 rtol 2%；单调性 | 必测点 `anchors_C/50_1C_2C_parent`：四次 Experiment 求解（C/50 合并、C/50 父提交、1C、2C） |
| I2：容量随倍率单调不增 Q(C/50) ≥ Q(1C) ≥ Q(2C) | 锚点上检三者顺序；每个抽样倍率 C 再与锚点比：Q(C/50) ≥ Q(C)；C≤1 时 Q(C) ≥ Q(1C)，否则 ≤；C≤2 时 Q(C) ≥ Q(2C)，否则 ≤ |
| I3：取 `Positive particle radius [m]` 自行计算 τ_p = R_p²/D_p(sto=0.5)，应与 1800 s 同量级或更大 | 必测点 `tau_p_sto=0.5` |

实施判断：
1. 不变量 I1 里写的量级界、负极同构，虽不在「数值核对方法」逐字列出，但该节以「I1」引用不变量，按协议 §2 纳入。
2. 抽样倍率的容量用快速路线：同一 BasicDFN 只离散一次，`Current function [A]` 设为 `[input]`，`solve([0, 5·3600/C])`，由模型自带 `Minimum voltage [V]` 事件在 `Lower voltage cut-off [V]`=2.0 V（参数集值，等于计划的 2.0 V）终止。锚点比较用的 Q(C/50)、Q(1C)、Q(2C) 也用同一快速路线另算一份（`anchor_fast`），避免两条路线混比；计划点名的锚点检查用 Experiment 路线。
3. 单调性无容差，按默认 rel 1e-6、abs 1e-12 放宽比较。
4. 「同量级或更大」取 τ_p ≥ 180 s（不低于 1800 s 的 1/10）。
5. 「修复前版本（父提交）」见 §0：只覆盖该函数。

## §3 第三路

| 计划原文 | 代码 |
|---|---|
| 公开路线：每次覆盖后新建参数对象；D 用 rtol 1e-10、atol 0，要求有限 | `pv(C,E)` 按配置缓存全新对象；`Dp()` |
| 独立参考：选夹住的两行按两点公式求值，不用另一个 Interpolant；端段显式外延 | `L(c)` |
| I1：默认 C 下 34 个节点、33 个区间各取 0.01/0.25/0.5/0.75/0.99 分点 | 必测点 `knot_k`（34）、`interval_k_frac_f`（165），跑与抽样点相同的全部检查 |
| I1：在表上均匀抽 256 个浓度（`default_rng(5889)`） | 被协议 N 点、协议种子取代：`t3_sample` c∈[x₀, x_last] 线性均匀 |
| I1：θ=0.3、0.5、0.9 的独立值 5.6673280002e-17、8.3317028568e-17、3.5626719998e-17 | 必测点 `theta_*` |
| I1：样本作为 `pybamm.Vector(theta_array)` 传入，与标量结果逐元素比 | 必测点 `vector_samples`（全部抽样点）、`vector_knots_intervals` |
| I1 模型接线：`BasicDFN().param.p.prim.D(Scalar(c), Scalar(298.15))` 与 L(c) 比（内部点） | 抽样点 `t3_I1_model_adapter_vs_L(c)_rtol1e-10` |
| I1 边界：c/C = 0, τ/2, τ, 2τ, 1−2τ, 1−τ, 1−τ/2, 1 与 L(C·clip(c/C,τ,1−τ)) 比 | 必测点 `adapter_c/C=…`（τ = settings `D__c_s` = 1e-10） |
| 直接调用 θ=0、x₀/(2C)、(x_last/C+1)/2、1 与端段线性外延比；θ=0、1 约 2.7385714278e-15、1.7387466637e-17 | 必测点 `direct_theta=…` |
| I2：E=500/1000/1500 重复默认 C 抽样，固定 θ 比较 | 抽样点 `t3_I2_E=…_equals_default` |
| I2：C=12256/15320/18384 × 同组 E，θ=c/C，在公共可容许区间 [x₀, 12256(1−τ)] 的节点与分点上等于 L(c) | 抽样点（及节点/分点必测点）落在该区间时 `t3_I2_C=…_E=…_vs_L(c)` |
| I3：D>0，落在夹它的两个表值之间（余量 1e-12×较大者） | `t3_I3_positive`、`t3_I3_within_bracketing_ordinates` |
| I3：每个内部节点 c=x_i±δ_i，δ_i=1e-6·min(相邻间距)，与单侧线性式比、\|D−d_i\| ≤ \|slope\|δ_i+1e-12 d_i | 必测点 `knot_i_delta_±1`（32×2） |
| I3：通量推论 g=0、±1e9，−Dg² ≤ 0 | `t3_I3_flux_sign_g=…` |
| I4：每个列出温度与 298.15 K 结果比 rtol 1e-12 | `t3_I4_T=273.15/323.15_equals_298.15` |

实施判断：
1. 分支：计划写明「单一相、单一线性插值、无充放电/滞回/温度分支」，33 段与外延、截断情形由必测点覆盖；N 点不再按段等分（k=1）。
2. 模型接线用 BasicDFN 的 `param`（`LithiumIonParameters`），用默认配置的参数对象求值。
3. 向量检查对全部抽样点一次 Vector 调用（与抽样同种子重新生成同一序列），逐元素比，容差沿用 D 的 rtol 1e-10。
4. 负极 HC 函数（计划称可作对照）未做：计划写的是「can serve as a control」，没给检查。

## §4 预计耗时与峰值内存（推测，未运行）

| 路 | 单点 | 10⁴ 点 | 依据 |
|---|---|---|---|
| t1 | ~10 ms（3 次 evaluate） | ~2 min＋网格＋1 次 Experiment（~10–30 s） | |
| t2 | ~0.5–1 s（1 次 BasicDFN 快速求解＋2 次 evaluate） | ~1.5–3 h；另加 4 次 Experiment 锚点（C/50 约 50 h 模拟时长，求解约十几秒） | BasicDFN 只离散一次 |
| t3 | ~40 ms（约 13 次 evaluate） | ~7 min＋必测约 450 点＋向量对照（再 10⁴ 次标量） | |

峰值内存：约 1–1.5 GB；JobSpec 申报 3 GiB。若 t2 单点 >1.4 s，按协议 §3 降 N 并记 §7。

## §5 偏离候选

- 第二路工况域写了 DFN/SPMe，钠离子模块只有 BasicDFN，只做 BasicDFN（候选：记 §7「计划写到的配置在库内不存在」）。
- 第二路「修复前版本」用同环境覆盖单个函数重建，而非安装父提交（父提交仅此一处源码差异，已核）。

## §6 冒烟后修正（2026-10-17）

- 无。首轮唯一拒绝为计划的界写法问题（§7），脚本忠实于计划，未改。

## §7 待裁定（冒烟拒绝归类，脚本不改）

1. **(P) 第二路 I1「在 sto ∈ [0,1] 上 min(表) ≈ 1.7e-17 ≤ D_p ≤ max(表) = 2.51e-15（端部外推量很小）」在 sto=0 端点不成立**。首轮唯一拒绝：必测点 `endpoint_sto=0.0`，观测 D_p = 2.7386e-15 ＞ 2.51e-15（超 9.1%）；可重放（同输入重算逐位相同）。原因：合并提交 `D_p.csv` 首两节点为 (131.578947, 2.51e-15)、(592.105263, 1.71e-15)，PyBaMM 线性插值在表外线性外推，c=0 处 2.51e-15 + 0.80e-15×131.58/460.53 = 2.7386e-15，与观测一致到末位。另一端 sto=1（c=15320）外推得约 1.74e-17，计划下界 1.7e-17 正是按外推值写的；上界却取表内最大值，两端处理不一致。计划写错之处：「≤ max(表) = 2.51e-15」与同句「端部外推」及库的默认线性外推矛盾（应为 ≈2.74e-15）。抽样点（sto∈[0.01,0.99]，c ≥ 153 > 131.6）不受影响，全跑只会有这一个必测点拒绝。
2. 本 PR 首轮未见 (S)、(V) 类。
