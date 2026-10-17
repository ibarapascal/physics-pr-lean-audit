# cp2851 实施日志（ΔQ1 统一抽查，CoolProp#2851）

> 目的：三路检查计划逐句→代码对应与全部实施判断，供协议 §6 核查与 §7 偏离登记。时间切片：2026-10-17 冒烟前写成；冒烟后只追加「冒烟后修正」节。
> 文件：`spotcheck.py`、`harness.py`（四份相同副本）、`job_smoke.json`、`job_full.json`。
> 读过的材料：协议；`track1.md`/`track2.md` 本 PR 节；`specs_codex/CoolProp_CoolProp_2851/spec.md` §2/§4/§5；合并提交 54a2e327 与 v8.0.0（ae81610e）的 `HelmholtzEOSMixtureBackend.cpp`；`env/setup_pca_py.sh`。

## 0. 环境

- compute node A `env/py2`：CoolProp 8.0.0 wheel，`__gitrevision__`=ae81610e。GitHub compare：ae81610e 比合并提交 54a2e327 **新 291 个提交、旧 0 个**（含本 PR）；两处 `calc_first_saturation_deriv`（带 SatL/SatV 的与无参委托版）在 v8.0.0 与合并提交逐字相同（md5 一致）。**视为有原生环境**，但不是合并提交本身的构建，其余代码有 291 个提交的差异——记为偏离候选 D0。
- 种子＝20261017＋6×10＋路号；判定/重放/失败规则同 pp3016 日志 §0；计划只给相对容差时绝对项取 0。
- 每次闪蒸新建 AbstractState（第三路明写「separate state instances」；第一、二路同样做，避免上一状态作初值影响重放）。

## 1. 第一路（track1.md）

分支（k=5）：A_ME_Q0、A_ME_Q1、A_R32R125_Q0、A_R32R125_Q1、B_pure_limit_ME。

| 计划原句 | 代码 |
|---|---|
| Methane&Ethane（z_CH4∈[0.1,0.9]，T∈[120,240] K，距临界/cricondentherm ≥5 K）；另取一组 R32&R125 或 R32&R1234yf（z∈[0.2,0.8]，T∈[230,330]） | `T1_MIX`、`t1_sampler`；二选一取 R32&R125（列在前） |
| 近临界点排除 | `t1_accept`：T−5、T、T+5 三处 Q=0 与 Q=1 的 QT 闪蒸都成功才接受，否则重抽（≤60 次）——判断 E1 |
| A：AS.update(QT_INPUTS,Q,T)，first_saturation_deriv(iT,iP)；另建状态 T±h 闪蒸取 p，中心差分 2h/(p₊−p₋)；h=0.2,0.1,0.05 K 做 Richardson 外推；Q=0、Q=1 都做；相对 1e-4 | `fd_slope`、`richardson`（S(h)=S₀+ah²+bh⁴ 过三点取 S₀）、`A.dT/dp_sat=…` |
| B：ε=1e-2,1e-3,1e-4 与 Methane 同 T、Q=0 的值比，检查差值随 ε 线性下降 | `t1_eval_B`、`_linear_rate`：相邻两档 d/ε 相等（协议默认容差）＋d 递减；Q=0、Q=1 两个混合物值都对纯甲烷 Q=0 |
| C：first_saturation_deriv(iHmolar,iT) 与沿 Q=0 线 h(T±h) 的中心差分比，相对 1e-4 | 只在 Q=0 分支；h 用 A 的三档＋Richardson（判断 E2） |
| 另核 Q=0.5 时抛异常 | 每点 `expect_raise` |

实施判断（第一路）：
- E1 「距临界/cricondentherm ≥5 K」没有给判定办法 → 用 ±5 K 处两条饱和线闪蒸都成功作操作化定义；被拒的抽样重抽（等价于在允许集合上均匀）。
- E2 C 的差分步长计划未写 → 沿用 A 的 0.2/0.1/0.05 K＋Richardson；体相 h 用 QT 状态的 `keyed_output(iHmolar)`（Q=0 时即组成为 z 的饱和液）。
- E3 B 的「线性下降」没有容差 → 按协议默认（相对 1e-6）判 d/ε 相邻两档相等。O(ε²) 项会使 d/ε 相差 ~ε 量级，这条检查**很可能因容差而普遍拒绝**——拿不准，交主会话裁定是否归「计划写法」类（见报告）。
- E4 B 的 T：在 [120,240] 均匀抽，要求纯甲烷与 ε=1e-2 混合物在 T±5 K 都能闪蒸（同 E1）；甲烷 Tc≈190.6 K，故有效区间约 120–185 K。

## 2. 第二路（track2.md）

分支（k=11）：{ME, NM（Nitrogen&Methane）, R32R125, PB（Propane&n-Butane）, CM（CO2&Methane）}×{Q0,Q1}＋I2_pure_limit_ME。必测点：z₁∈{0.05,0.3,0.5,0.7,0.95}（计划离散值）×每个分支，T=0.75·T_c,mix。

| 计划原句 | 代码 |
|---|---|
| 五种混合物，z₁∈{0.05,…,0.95}，T 取 0.6–0.9 倍混合物临界温度 | `t2_sampler`：z₁∈U[0.05,0.95]、T=U[0.6,0.9]·T_c(z)；T_c 用 `T_critical()`，失败（多个临界点）时退回 `T_reducing()` 并记录（判断F1） |
| I1：中心差分 2ΔT/(p(T+ΔT)−p(T−ΔT))，ΔT=0.01,0.02,0.05 K 做 Richardson；容差 1e-3；Q=0、Q=1 | `I1.dT/dp_sat=Richardson FD` |
| 附带：first_saturation_deriv(iP,iT)·first_saturation_deriv(iT,iP)=1 | `I1.dP/dT*dT/dP=1`（无容差→默认） |
| I2：z₂→0 序列（1e-2,1e-3,1e-4）上与 Methane 比，偏差随 z₂ 线性下降 | `t2_eval_I2`，T=U[0.6,0.9]·T_c(甲烷)，Q=0、Q=1；线性判法同第一路 E3 |
| I3：first_saturation_deriv(iHmolar,iT) 与 (h(T+ΔT)−h(T−ΔT))/(2ΔT)（沿 Q=0）比，rtol 1e-3；抛异常或返回两相值记为发现 | Q=0 分支；步长与 Richardson 同 I1；抛异常→记录「finding」并按协议计失败（判断 F2） |

实施判断（第二路）：F1 见上；F2 计划把「抛异常」当发现，但协议规定异常不算违反 → 计失败、另行记录；F3 I2 只用 Methane&Ethane（计划只点名 Methane）；F4 Richardson 三档步长比不等（2、2.5），用 h²、h⁴ 两项多项式过三点精确外推。

## 3. 第三路（spec.md §5）

分支（k=10）：{ME, MEP（Methane&Ethane&Propane）, Water, nPropane, R410A}×{Q0,Q1}。必测点：§4 表的全部离散值（ME 140–220 K 每 10 K × 甲烷分数 0.1/0.3/0.5/0.7/0.9；MEP 160–220 K × 两组组成；Water 300/400/500/550；n-Propane 230/270/300/330；R410A 250/280/310；各两端点），另加 180 K、(0.7,0.3) 的内部干度拒绝测试点。随机点：ME T∈U[140,220]、z∈U[0.1,0.9]；MEP T∈U[160,220]、组成在两组中等概率取；纯流体按表的 T 区间均匀。

| 计划原句 | 代码 |
|---|---|
| \|a−b\| ≤ r·max(\|a\|,\|b\|)+ε_abs；斜率 ε_abs=1e-8·T/p；性质导数 ε_abs=1e-8·F*/X* | `rel_ok`、`_close3` |
| I1：δT=0.05,0.025,0.0125 K 的 QT 差分，相邻参考值相对 5e-4 内收敛，再以 r=1e-3 比；单独复现 PR 的 0.05 K、0.5% | `t3_I1`；纯流体也用 r=1e-3？——计划 I1 只写混合物 r=1e-3；纯流体在 I2 用 1e-4（代码对纯流体 I1 用 1e-4，判断 G1） |
| I1：PQ_INPUTS，δp/p=1e-4,5e-5,2.5e-5 的 [T(p+δp)−T(p−δp)]/(2δp) | 同函数 |
| I1 加权残差：液/气组成与饱和密度；以 specify_phase＋PT_INPUTS 建均相态；中心密度与平衡分支一致；有限差分 ln φ_i 得 a_{k,i}；\|A_T S+A_p\|/(\|A_T S\|+\|A_p\|)<2e-3；根切换使该档无效；组成归一 1e-10、max\|ln(f_L/f_V)\|<1e-6 | `t3_I1_weighted`；权重：泡点线用气相组成、露点线用液相组成（判断 G2）；根切换判据：差分态密度偏离平衡分支 >5% |
| I2：纯流体端点 S 与 TΔv/Δh（饱和 keyed 输出）r=1e-10；差分 r=1e-4；同温泡/露点斜率 r=1e-8；R410A 只做表达式检查，base vs merge 隔离构建比 | `t3_I2`；R410A 的 base vs merge 节点上没有隔离父提交构建 → 只记录未跑 |
| I3：iDmolar、iHmolar、iSmolar、iUmolar 及质量版，对 iT 与 iP；与固定 q,z 的 keyed_output 差分比（混合物 1e-3、纯 1e-4）；链式法则用 first_partial_deriv（1e-10）；倒数与恒等 1e-10 绝对 | `t3_I3` |
| I4：内部干度 0.1/0.5/0.9（180 K、0.7/0.3）四种形式都须抛；端点上 (iT,iDmolar) 须抛 | 拒绝测试在必测点；(iT,iDmolar) 每点核 |
| I5：二元反序、三元 6 种排列；p_sat、斜率、重排后相组成 r=1e-7，组成绝对 1e-9 | `t3_I5`（只对混合物） |
| I6：second_saturation_deriv(F,iP,iP)（F=T,Dmolar,Hmolar,Smolar,Umolar）与 first_saturation_deriv(F,iP) 的压力中心差分比；参考收敛 0.1%；r=5e-3，ε_abs=1e-8F*/p²；不支持的调用单独报告 | `t3_I6`；不支持→records |

实施判断（第三路）：
- G1 纯流体端点 I1 差分比较用 1e-4（计划 I2 写「independently use the flash finite differences above at r=1e-4」）。
- G2 A_T、A_p 的权重：spec 第 1 节（不在可读范围）里的定义我看不到 → 按相平衡全微分＋Gibbs–Duhem 推得的标准形式（泡点 y、露点 x）。拿不准。
- G3 均相态密度与平衡分支一致的容差计划未给 → 默认相对 1e-6。
- G4 F* 取 max(1, \|F 中心值\|)。

## 4. 耗时与内存估计（推测；compute node A 探针：混合物 QT 闪蒸 ~0.7 ms，AbstractState 新建 ~1–3 ms）

| 路 | 单点 | N=100（含必测） | N=10⁴ |
|---|---|---|---|
| 1 | ~50–80 ms（含接受性闪蒸） | ~10 s | ~15 min |
| 2 | ~30–50 ms | ~10 s | ~10 min |
| 3 | ~0.5 s（~250 次建态＋闪蒸） | ~2 min（必测 ~150 点） | ~1.5 h |

峰值内存 ~200–300 MB/进程。全量三路并行（申报 2 GiB、3 核）。

## 5. 偏离候选

- D0 环境是 v8.0.0 wheel（含本 PR、改动函数逐字相同），不是合并提交构建。
- D1 第一路 B、第二路 I2 的「线性下降」按协议默认容差判，预计普遍拒绝（E3）。
- D2 第三路 I2 的 R410A base-vs-merge 未跑（没有父提交隔离构建）。
