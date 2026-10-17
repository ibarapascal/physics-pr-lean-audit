# pp3016 实施日志（ΔQ1 统一抽查，pandapower#3016）

> 目的：记录三路检查计划逐句→代码的对应与全部实施判断，供协议 §6 核查与 §7 偏离登记。时间切片：2026-10-17 冒烟前写成；冒烟后只追加「冒烟后修正」节。
> 文件：`spotcheck.py`（入口）、`harness.py`（抽样/判定/重放/输出框架，四份 PR 目录各一份相同副本）、`job_smoke.json`、`job_full.json`。
> 读过的材料：协议、`track1.md`/`track2.md` 本 PR 节、`specs_codex/e2nIEE_pandapower_3016/spec.md` §2/§4/§5；合并提交 0c6a1f0a 的 `pd2ppc_zero.py`、`trafo_create.py`、`build_branch.py`、`build_bus.py`、`shortcircuit/*`、`test_1ph.py`、`test_iec60909_4.py`（只取夹具数据）；`env/setup_pca_py.sh`。

## 0. 环境与通用约定

- 原生环境：compute node A `env/pp3016`（`uv pip install pandapower@git+…@0c6a1f0a`，pandapower 3.4.0，numpy 2.3.5，scipy 1.16.3，pandas 2.3.3）。**有原生环境**。numba 未装（只影响速度）。
- 种子＝20261017＋8×10＋路号；每点随机数 `default_rng([seed, 分支号, 点号])`，点参数先抽好再求值，拒绝点用同一参数原样重算一遍，坏检查与观测值相同才计拒绝（不同则计「不可重放」，不计拒绝）。
- 分支等分：N 按分支数 k 等分（余数给前几个分支），每分支 ≥⌊N/k⌋。计划写明的离散值作为必测点另加（不占 N），JSON 里 `n_must` 单列。
- 判定：点内任一检查超容差＝拒绝；没有超容差但出现异常/NaN/不收敛＝失败（单列，不算违反）。
- 容差：用计划写的；计划只给相对容差时绝对项取 0；计划完全没给时用相对 1e-6、绝对 1e-12（协议默认）。复数：第二、三路写「逐分量」→实部虚部分别判；第一路未写→按复数模判。
- 网络构建一次、改表再算（`net.trafo.at[...]=`、`net.ext_grid.loc[...]=`）：实测单次建网 ~58 ms、单次 calc_sc ~7 ms（compute node A 秒级探针），改表比重建快一个数量级。第三路计划本身允许「public table assignments」；第一、二路记为实施判断。

## 1. 第一路（track1.md「数值核对方法」＋附加）

分支（k=4，各 2500 点）：`A_Dyn`、`A_YNd`、`C_zigzag`、`ADD_pf3ph_Dyn`。必测点：A 写明的 Z_N ∈ {0, 0.01, 0.1+0.05j, 10} Ω，在中点网络（MV 22.5 kV、50.05 MVA、vk 9%、vkr 0.9%、si0 0.725）上各跑一次（A_Dyn、A_YNd 各 4 点）。

| 计划原句 | 代码 |
|---|---|
| 常用工况域：中压 10–35 kV、低压 0.4 kV、高压 110–132 kV；0.1–100 MVA，vk 4–14 %，vkr 0.3–1.5 %；R_N ∈ [0,1000] Ω；X_N ∈ [0,100] Ω；si0 ∈ [0.5,0.95] | `t1_tr`、`t1_s_*`：各量均匀抽样 |
| A：Dyn 网络对 Z_N × case ∈ {max,min} 检查 ΔZ0 = 3Z_N，容差绝对 1e-9 Ω 或相对 1e-9 | `t1_e_dyn`：`A.dZ0=3ZN[case]`，`mode="max"`（取两者较宽） |
| 正序 Z1_th 与 Z_N 无关（不变量 A） | `A.Z1_indep[case]` 同容差 |
| YNd 网络按并联式反推常数（Y_rest 由断开该变压器的同网络算出），对不同 Z_N 一致到相对 1e-9 | `t1_e_ynd`：变压器 `in_service=False` 求 Z0_rest；`const=1/(1/Z0−1/Z0_rest)−3Z_N`，对 Z_N∈{(0,X),(R/2,X),(R,X)} 与 Z_N=0 比 |
| B：与 √3·c·U_n/\|2Z1+Z0\| 比，相对 1e-9 | 每次求解都核（`B.ikss=…`），c 用 IEC 表（下 §判断 J3） |
| B：R_N 扫描检查单调 | 每个 case 在 R_N∈{0, R/2, R}（X_N 固定为抽样值）核 ikss 不增 |
| C：三角形馈电中压母线＋ZNyn/ZNd，比较两者 b_hv 的 Z0（按位相同） | `t1_e_zz`：同一网只改 `vector_group`，`C.ZNyn_vs_ZNd_hv_bitwise` 用 `==` |
| C：与 si0·K_T·Z0_leak+3Z_N 比（K_T 按 IEC 式自己算，或用 case 切换验证 3Z_N 部分不随 K_T 变） | 两种都做：case=max 用 IEC 式 K_T 核绝对值；max、min 都核 ΔZ0=3Z_N |
| C：改电网零序参数检查 ZNyn 低压 Z0 不变 | 外网 x0x/r0x0 由 1.0/0.1 改为 2.0/0.2 再算，`C.ZNyn_lv_Z0_indep_of_grid0` |
| 不变量 C：ZNd 低压母线 I″_k1 = 0 | `C.ZNd_lv_ikss=0`（绝对 1e-12 kA） |
| 附加：runpp_3ph 单相负荷下比较「Dyn＋Z_N」与「Dyn＋在低压中性点外接等值 Z_N」的中性点位移 | `t1_e_pf`：网 A 变压器带 rn/xn；网 B 变压器 Z_N=0，低压母线后接一条正序≈0（1e-8 Ω/km）、零序 r0=3R_N、x0=3X_N 的 1 km 线路再挂负荷；比较负荷母线 V0=(Va+Vb+Vc)/3（标幺） |

实施判断（第一路）：
- J1 A 的 Dyn 网：中压侧 U∈[10,35] kV、低压 0.4 kV（计划写「Dyn 变压器的低压母线」，PR 测试为 20/0.4 kV）；YNd 网：高压 110–132 kV/中压 10–35 kV，故障在高压母线；曲折网：110–132 kV 外网→YNd 主变→中压 10–35 kV→ZN* 接地变→0.4 kV 辅助母线。
- J2 计划没给区间的输入按 PR 测试固定：外网 Dyn 用 1000 MVA、YNd/曲折用 5000 MVA，rx 0.1、x0x 1.0、r0x0 0.1；vk0=vk、vkr0=vkr；mag0_percent 100、mag0_rx 0；pfe_kw=i0_percent=0；shift_degree 按 PR 测试（Dyn/YNd/ZN* 330°）。主变与接地变的 sn、vk、vkr、si0 各自独立在同一区间抽样。
- J3 c 取 IEC 60909-0:2016 表 1（计划给的两种做法之一）：<1 kV 且 10% 容差 c_max 1.10/c_min 0.90，其余 1.10/1.00；lv_tol_percent 用默认 10。
- J4 IEC K_T＝0.95·c_max/(1+0.6·x_T)，x_T 用正序 vk/vkr，c_max 取接地变低压侧母线（0.4 kV→1.10）；只在 case=max 核绝对值（IEC §6.3.3 只对最大短路电流用 K_T，min 情形计划给了「case 切换」替代做法）。Z0_leak 用接地变高压侧额定电压归算（额定电压＝母线电压）。
- J5 C 节计划没写容差，用协议默认（相对 1e-6、绝对 1e-12）；A、B 用计划的 1e-9。
- J6 inverse_y：计划只在 API 里列 True/False，没区分分支 → 每点均匀抽一个。
- J7 附加检查的负荷大小计划未给、库文档 create_asymmetric_load 默认 0（无意义）→ 取 a 相 p=0.01·S_T、q=0.005·S_T；只加载 a 相（计划写「单相负荷」）。容差未给→默认。runpp_3ph 收敛容差取 1e-10 MVA（两网同设，降低数值噪声，记为判断）。R_N 抽到数百 Ω 时 0.4 kV 侧中性点位移很大、潮流可能不收敛 → 计失败。
- J8 「外接等值 Z_N」的建模：pandapower 没有中性点阻抗元件，用「零序 3Z_N、正序近 0 的串联线路」表示；该线路正序 1e-8 Ω 带来的电压差 ≪ 1e-6 相对（推测）。

## 2. 第二路（track2.md「数值核对方法」）

分支（k=20，各 500 点）：{Dyn 辐射、ZNd 辐射、ZNyn 辐射、YNd 并联、Zd}×case{max,min}×inverse_y{True,False}。计划把 case、inverse_y 列为工况域选项 → 按协议「选项值」作为分支等分。

| 计划原句 | 代码 |
|---|---|
| U_n 0.4–132 kV；配变 0.1–2.5 MVA、主变 10–100 MVA；vk0 4–14%，vkr0 0.3–1.5%，si0 0.5–0.95；R_N 0–200 Ω，X_N 0–50 Ω；S″_k 100–10000 MVA，x0x 0.5–3 | `t2_base`、`t2_s_*` 均匀抽样 |
| I1：随机 Z_N（R、X 同时非零），辐射型母线检差值 = 3Z_N，atol 1e-9 Ω | `_t2_I1_radial`：ΔR0、ΔX0 分别绝对 1e-9；正序 R1、X1 不变同容差 |
| I1：并联型母线用三组 Z_N 拟合 (Z0_p, Z0_branch) 后用第四组验证，rtol 1e-9 | `t2_make_e_par`：Z_N∈{0,Z_b,Z_c,Z_d}；由前两点解二次方程得 u=1/Z0_p 两根，用第三点挑根，求 Z0_branch，预测第四点 Z0，复数模相对 1e-9；第三点拟合残差只记录 |
| I2：按式子重算 ikss 与 ikss_ka 比，rtol 1e-9 | 每次求解核 |
| I2：R_N ∈ {10,100,1000} Ω 检单调性与渐近式，R_N=1000 Ω 时相对偏差 <1% | `_t2_I2`：X_N=0 扫描，严格递减；\|I(1000)/(cU/(√3·1000))−1\|<0.01 |
| I3(a)：ZNyn 改外网 x0x_max 两个值，LV rk0/xk0 差 ≤1e-9；改 HV 侧 R_N 不变 | ZNyn 分支：x0x 换成第二个抽样值 x0x2（max/min 两列同改）；Z_N 与 0 比 |
| I3(b)：接地变矢量组设为 "Zd"，MV 母线 ikss 应为 0（atol 1e-9 kA） | `Zd` 分支：Z_N=0 与抽样 Z_N 两次 |
| I3(c)：ZNd 的 LV 母线 I″_k1 = 0 | ZNd 分支，计划未给容差 → 默认绝对 1e-12 kA |

实施判断（第二路）：
- K1 电压只给一个区间 0.4–132 kV → 每个网需要的 2 或 3 个电压等级各在 [0.4,132] 均匀抽、降序排为 HV>MV>LV（所以会出现很接近的电压比，属于计划域内）。
- K2 未给的输入按 PR 测试固定：vk=vk0、vkr=vkr0；外网 rx 0.1、r0x0 0.1；mag0 100%、mag0_rx 0；Zd/ZNd/ZNyn 网的主变用 YNd（计划写「三角形供电」，PR 测试同）。
- K3 I2 的 R_N 扫描令 X_N=0（渐近式写的是纯电阻）。
- K4 I2 渐近与单调只在「唯一接地经 NER」的辐射型分支做（Dyn、ZNd、ZNyn 的中压母线），并联型只核电流式。

## 3. 第三路（spec.md §5，含 Additional changed branches、T 模型路径、pf_3ph）

分支（k=20，各 500 点）：`Dyn`、`YNd_oneport`、`Yyn`、`YNy`、`YNyn`、`Yzn`、`Yy`、`Yd`、`Dy`、`Dd`、`ZNyn`、`ZNd`、`ZNy`、`Zd`、`I5_parallel_service_api`、`PSU`、`Tmodel_path`、`PF3ph_Dyn`、`PF3ph_YNyn`、`PF3ph_Yzn`；另有 `API_contract`（只有 1 个确定性必测点，不占 N，结果并入同一 JSON 的 `api_contract`）。每个抽样分支另加 5 个确定性必测点（fixture 1/2 原值，s 取 {0.01,0.1,0.5,0.9,0.99}，3|Z_N|/|Z0T| 取 {0.1,1,10,100,1}，45° 角）。

| 计划原句 | 代码 |
|---|---|
| deterministic values first; 30 seeded log-uniform samples per connection over positive impedance/rating ranges, uniformly sampling s and dimensionless R/X ratios | 必测点＝确定性值；随机点数按协议改为 N/k；正的阻抗/额定量（sn、vk、vkr、vk0、vkr0、mag0、R_N、X_N、S″_k、电压）对数均匀，s、各 R/X 比、x0x 均匀（`t3_sample`）；vkr≥vk 的组合重抽 |
| For every new/affected connection, cross both solver choices and both SC cases | 每次求解 inverse_y True/False 都跑并直接比（`run3`）；case 每点均匀抽 max/min（见判断 L2） |
| I1 fixture 1：pure R、pure X、combined，exact complex 3Z_N increment and unchanged orthogonal component | `t3_e_dyn` 四次求解、三个增量逐分量（rtol 1e-7、atol 1e-8 Ω） |
| fixture 2 replace grounding unit by YNd → isolated HV one-port | `t3_e_ynd` |
| Yyn / YNy：default pi one-port Z_m+Z_ℓ+3Z_N；positive sequence from the opposite ungrounded winding | `t3_make_e_pi`：Yyn 外网在 HV、故障在 LV；YNy 外网在 LV（中压级）、故障在 HV |
| YNyn：两端驱动点阻抗对独立组装的 T 矩阵（3Z_N 全在 HV 臂）；改对侧电网零序终端 | `t3_e_ynyn`＋`_ynyn_oracle`：A=sK_TZ_ℓ+3Z_N、B=(1−s)K_TZ_ℓ、M=Z_m；HV 驱动点 Z0g∥(A+M)，LV 驱动点 B+M∥(A+Z0g)；外网 x0x/r0x0 换第二组再核 |
| Yzn：reverse zigzag one-port（接地低压漏抗＋全 3Z_N）；baseline calibration 单独处理 | 增量 3Z_N 计入判定；基线 (1−s)K_TZ_ℓ,L 只记录（`records`，不计拒绝） |
| four ungrounded groups：zero/nonzero neutral equal | `Yy/Yd/Dy/Dd` 分支，全母线全列比 |
| I2：ZNyn/ZNd/ZNy/Zd 的 Z0_H=sZ_ℓ,H+3Z_N；ZNyn 的 Z0_L=Z_m,L+(1−s)Z_ℓ,L | `t3_make_e_zz` |
| change only grid zero-sequence strength at the remote side → no cross-side change | 外网 x0x/r0x0 换组重建网，MV Z0（及 ZNyn 的 LV Z0）不变 |
| vary HV neutral → LV Z0 unchanged；vary m → HV Z0 unchanged | ZNyn：Z_N 与 0 比；改 mag0_percent 为 `mag0_alt` |
| no auxiliary source：LV ikss < 1e-6 kA for ZNd/ZNy/Zd | `I2.lv_ikss<1e-6` |
| Attempt ZN separately and record；mixed-case spellings；ZNyn11 → validation error | `API_contract`：znyn/ZNYN/Znyn/zNyN 与 ZNyn 结果相等；ZNyn11 须抛 ValueError；ZN 只记录结果 |
| I3：重算电流（c 表 1.1/1.0；<1 kV：10% 1.1/0.9，6% 1.05/0.95） | `i3_3` 的 `I3.ikss_recon`，lv_tol 每点抽 {6,10} |
| I3：3ph、2ph 在中性点改动前后电流与 Z1 不变 | `I3.3ph/2ph_unchanged_by_ZN` |
| I3：isolated return paths 上 R 或 X 单调扫描；大 R 对比渐近，相对误差界＝剩余串联阻抗/(3R_N) | R 扫描 3R/\|Z0T\|∈{0.1,1,10,100}、X 扫描 {0.1,1,10}，均不增；比值 100 处 \|I/渐近−1\| ≤ \|2Z1+Z0(0)\|/(3R)+1e-7。用于 Dyn、YNd、Yyn、YNy、Yzn、ZN*（YNyn 不是孤立回路，不扫） |
| I4：每个接地情形对网络基准 {1,10,100,0.1,123} MVA 重复，欧姆/kA 不变 | `bases3` |
| I4：HV/LV 额定/母线电压比 0.9–1.1 独立扫；HV shunt 用物理 HV 漏抗 K_T U_rH²(…)/(100 S_T p)；tap −2/0/+2、2.5%/档、各 tap 侧 | 曲折分支每点抽 ratio_h、ratio_l、tap_side、tap_pos；oracle 用 U_rH=vn_hv·(1+0.025·tap)（tap 在 HV 时） |
| I4：只改 ZNd 的辅助 LV 母线额定电压，HV 物理 Z0 不变 | ZNd 分支 `I4.ZNd_hv_Z0_indep_dummy_lv_vn`（母线电压×U[0.9,1.1]） |
| I5：parallel=p 与 p 个显式行；先同 Z_N 记录，再把成组行的中性点设为 Z_N/p 做等价比较；零中性点；普通/曲折路径 | `t3_e_i5` (a)(b)(c)，路径 Dyn 或 ZNd 每点抽，p∈{1,2,3,4} |
| switch one duplicate out of service vs delete；out-of-service near-boundary | (d) 仅 p≥2（「duplicate」）；退出那台的 s 取 0.01 或 0.99 |
| shuffle row/index order；multiple HV/LV shunts accumulate；mix groups | (e) 用置换顺序与非连续 index 建网；(f) 同母线两台不同参数（Dyn 路径两台 Dyn、ZNd 路径 ZNd+ZNyn）与并联式比 |
| four creator APIs（create_std_type 注册标准型）、scalar bulk、vectors、mixed zero/NaN/nonzero、omitted columns、combined R/X | (g) 三孤岛网，四种接口结果逐岛相等；NaN 岛＝零岛；标量批量＝逐个；完全不传 rn/xn＝零 |
| I6：finite nonnegative currents、Re Z0 ≥ −ε、两求解器直接比；two-port oracle Y 的 Hermitian 部分半正定 | `run3`；`I6.oracle_Y0_hermitian_psd`（YNyn） |
| Additional：iec_60909_4() 的 YNd PSU，rn 省略/0/非零 × xn 0/非零，增量 3Z_N（不是 k_ST·3Z_N） | `_psu_net` 照抄 t1＋G1 数据（只含这一台 PSU 与其发电机）；rn 省略＝rn 0；三种非零组合核增量 |
| PSU=True with each zigzag group：不支持应显式拒绝；接受则中性点只算一次、LV 独立 | 抛异常→记录为拒绝（通过）；不抛→核 ΔZ0_HV=3Z_N，ZNyn 再核 LV Z0 不随 HV 中性点变 |
| mixed PSU flags in one vector group | 同网两孤岛：PSU YNd＋普通 YNd（外网在其三角侧），两处增量都核 |
| T-model path：runpp(trafo_model="t") 后 calc_sc(bus=target, use_pre_fault_voltage=True)，再用 "pi" 重复 | `t3_e_tmodel`，每点抽一个联结组 {Dyn,YNd,Yyn,YNy,YNyn}；T/pi 各对各自的 oracle，并核增量 |
| pf_3ph：Dyn、YNyn、Yzn 不对称星形负荷，负荷相轮换 a/b/c；由电压与给定复功率算相电流，V0、I0 与物理 T 电路（全 3Z_N 在接地臂，无 IEC 修正）比 | `t3_make_e_pf`：V0_pred=−Z0_phys·I0（负荷电流方向），残差 ≤1e-6 pu＋1e-5 相对 |
| MVA base changes；balanced loads 与中性点无关；HV zigzag groups 在 PF 模式被拒绝 | 同函数；拒绝检查放在 `API_contract`（要求抛 NotImplementedError/ValueError） |
| SC 容差 rtol 1e-7、atol 1e-8 Ω/1e-8 kA；开路电流 1e-6 kA；PF 1e-5 相对、1e-6 pu | `RT3/AZ3/AI3` 与各处 |

实施判断（第三路）：
- L1 区间：LV 0.23–1 kV、MV 3.3–33 kV、HV 66–220 kV、sn 0.1–100 MVA、vk/vk0 2–20%、vkr/vkr0 0.1–2%、R_N/X_N 0.01–100 Ω、mag0 100–10000%、mag0_rx 0–1、s 0.01–0.99、S″_k 100–10000 MVA（min≤max）、rx/r0x0 0.05–0.3、x0x 0.5–3，均取 §4 表；R_N=0、X_N=0 作为基线每点都算。
- L2 「cross both SC cases」：每点 case 均匀抽一个（不是每点两种都跑），inverse_y 每次两种都跑（I6 要直接比）。如需每点两种 case 都跑，单点代价约翻倍。
- L3 IEC K_T 只在 case=max 乘到漏抗上（IEC §6.3.3），不乘励磁阻抗 Z_m；K_T 的 c_max 取该变压器低压侧母线。spec §2 没写 K_T 是否作用于 Z_m，这是判断（拿不准，见报告）。
- L4 Z_m＝mag0/100·\|Z_ℓ\|(rx+j)/√(1+rx²)，\|Z_ℓ\| 用未修正的零序短路阻抗（按 create 文档「z_mag0/z0」与「R/X」定义）。
- L5 外网零序阻抗（oracle 用）：Z_Q=c·U²/S″（max 用 c_max 与 *_max 参数，min 用 c_min 与 *_min）；pf_3ph 的 YNyn oracle 用 c=1.1 与 *_max（库在 pf 模式的外网约定）。
- L6 Yzn「grounded LV leakage」理解为 (1−s)K_TZ_ℓ,L＋3Z_N；按计划「separately」只记录基线，不计拒绝。
- L7 PSU 夹具只取 iec_60909_4() 里的 t1 与其发电机（计划要「valid YNd PSU fixture」且要能读出 HV 增量；整张网是网状、增量不是 3Z_N 的一端口量）。R_N、X_N 在 1–100 Ω 对数均匀（§4「MV examples 1–100」）。
- L8 pf 分支的中性点大小：§4 只说「sufficiently small」→ 取 3\|Z_N\|/\|Z_ℓ\| 对数均匀于 [0.1,1]、角度均匀于 [0,π/2]（YNyn 用 HV 侧 \|Z_ℓ\|）；负荷 1–20% 额定、功率因数 0.9–1、50/60 Hz 均匀抽。
- L9 I5 (d) 只在 p≥2 做；(e) 打乱时各行 R_N 乘 (1+0.1k) 使行可区分。
- L10 「Use these as arm-placement checks」：T 路径对 T 电路 oracle（Dyn：B+3Z_N+A∥M；YNd：A+3Z_N+B∥M；Yyn：B+3Z_N+M；YNy：A+3Z_N+M；YNyn：同 I1），pi 路径对 pi oracle。
- L11 每点的 Z_N 取「组合」值做绝对比较；纯 R、纯 X 增量只在 fixture 1（计划只在 fixture 1 写了三种）。

## 4. 耗时与内存估计（推测，依据 compute node A 秒级探针：建 2 母线网 ~58 ms、calc_sc ~7 ms、runpp_3ph ~27 ms）

| 路 | 单点（均值） | N=100（含必测） | N=10⁴ |
|---|---|---|---|
| 1 | ~0.13 s | ~15 s | ~0.4 h |
| 2 | ~0.08 s | ~10 s | ~0.25 h |
| 3 | ~0.5 s（I5 ~2.5 s、PSU ~1.5 s、其余 0.1–0.6 s） | ~1.5 min | ~1.5 h（上限估 2.5 h） |

峰值内存：单进程 ~300–400 MB（pandas/scipy＋小网，推测）；全量三路并行 ~1.2 GB → JobSpec 申报 2 GiB。全量 JobSpec 三路并行（各自日志），墙钟按第三路计。

## 5. 偏离候选（交主会话裁定）

- D1 第三路「30 samples per connection」按协议改为 N/k 均匀＋对数均匀抽样（协议 §3 优先）。
- D2 第三路 case 每点抽一个而非每点 max/min 都跑（L2）。
- D3 第一、二路用改表重算代替逐次重建网（性能；结果应等价，推测）。

## 6. 冒烟后修正（2026-10-17，首次冒烟 N=100 之后；协议 §6 只许修脚本错误）

首次冒烟输出：compute node A `dq1/out/pp3016_t{1,2,3}.json`。第一路只有 `ADD_pf3ph_Dyn` 异常（25/25 失败、n_checks=0），第二路 52/100 拒绝。诊断只用冒烟 JSON、库源码（compute node A `env/pp3016` 只读查看）与按种子重算抽样参数（纯随机数，不跑潮流）；修正效果由冒烟二验证，此前标「推测」。

- **S1 第一路附加检查 `t1_e_pf`（`spotcheck.py` 该函数全段）**
  - 症状：25 点全部 `LoadflowNotConverged: Power Flow nr did not converge after 30 iterations!`，抛自 `runpp_3ph.py:646`——外层序分量定点迭代 30 次未达 `outer_tolerance_mva=3e-8`（写死，标幺值，与 `tolerance_mva` 无关）。
  - 原因：J7 取负荷 0.01·S_T（1 kW–1 MW）单相接 0.4 kV，计划域 R_N∈[0,1000] Ω 下单相负荷回路要经中性点，回路可传最大功率约 V_ph²/(4|Z_N|)（1000 Ω 时约 13 W）。按种子重算这 25 点：负荷是可传最大功率的 116–70 800 倍，物理上无解，潮流必然不收敛。另外 J8 外接线路正序 1e-8 Ω 使导纳比负荷大十几个数量级，`tolerance_mva=1e-10` 低于双精度可达的失配下限（推测），小负荷时也收敛不了。
  - 改法（只改 J7/J8 的实施取值，计划映射不变：仍是同一 Dyn 网对比「变压器 rn/xn」与「低压中性点外接等值 Z_N」的负荷母线 V0，判据与容差不变）：
    1. 负荷 S=0.01·V_ph²/(|Z_N|+|Z_T|)，V_ph²=0.4²/3 kV²，|Z_T|=vk/100·0.4²/S_T，p:q 仍为 2:1、只接 a 相——中性点位移约 1%，单相负荷可供、定点迭代收缩因子约 0.01。
    2. 网络基准 `sn_mva=100·S`：外层判据 3e-8 是标幺值，基准随负荷缩放后失配相对负荷约 3e-6，同时高于浮点失配下限。查过 `pd2ppc_zero.py:253–270、324`：pf_3ph 变压器用 3·sn_mva、中性点用 sn_mva，两者比值与 sn_mva 取值无关，改基准不会掩盖或制造该处差异。
    3. 外接线路正序取 1e-6·(|Z_N|+|Z_T|) Ω（r=x 各取其 1/√2），零序仍为 3R_N、3X_N；它对 V0 的影响约 2e-8（相对），低于默认容差 1e-6。
    4. `runpp_3ph` 用库默认 `tolerance_mva`（1e-8），因为 1e-10 达不到，而且最终精度由写死的外层 3e-8 决定。
  - 对应计划句：track1「附加……建议用 `runpp_3ph` 在单相负荷下检查「Dyn＋Z_N」与「Dyn＋在低压中性点外接等值 Z_N」的中性点位移是否一致」。计划没给负荷大小、基准和收敛设置，这几项原属实施判断 J7/J8，现由本节替换。
- 第二路：无 S 类错误，代码未改。逐条归类见下节。

## 7. 待裁定（冒烟拒绝归类；只记录，不改脚本）

- **P1 第二路 I2「R_N=1000 Ω 时相对偏差应 <1%」（`_t2_I2` 的 `I2.asymptote_RN1000<1%`）**——Dyn_radial 13/20、ZNd_radial 10/20、ZNyn_radial 9/20 点拒绝；这三个分支没有其他检查拒绝。
  - 证据：Dyn 的 8 个拒绝样例，观测偏差与第一性原理式 1−3R_N/|2Z_1+Z_0(0)+3R_N| 吻合到 4 位（例：U_lv=54.37 kV、S_T=0.671 MVA、vk0 6.96%，|2Z_1+Z_0,T|=924 Ω，预测与观测偏差都是 0.0635）；同一点 `I2.ikss_formula`（rtol 1e-9）全部通过，库的电流与 IEC 式一致。
  - 归类 P（计划写错）：计划同节自己写了「R_N ≫ |Z_1|,|Z_0,T| 时……相对偏差 O(|2Z_1+Z_0,T|/(3R_N))」，<1% 需要 |2Z_1+Z_0,T| ≲ 30 Ω。但在计划自己的工况域里（U_n 到 132 kV、配变 0.1–2.5 MVA、vk0 到 14%），|Z_T| 可达数千 Ω，即使 20 kV、0.1 MVA、14% 也有约 560 Ω。固定 1000 Ω 与 <1% 这一句，和计划自己的误差式在其域内矛盾。K1（各级电压在 [0.4,132] kV 独立抽）让高压配变组合更常出现，但不是成因。
  - 可选裁定：照原样报拒绝并归「计划写错」；或者只在满足前提 |2Z_1+Z_0,T|/(3·1000)≤0.01 的点上判。后者改了计划映射，须记入协议 §7。
- **V1 第二路 I3(b)「接地变矢量组设为 Zd，MV 母线 ikss 应为 0（atol 1e-9 kA）」（`t2_make_e_zz("Zd")`）**——Zd 全部 20/20 点、两态（Z_N=0 与抽样 Z_N）都拒绝，ikss 在 0.06–2.9 kA，且随 Z_N 增大而下降，说明中性点阻抗确实计入了。
  - 不是脚本错误：MV 母线由 YNd 主变的三角侧供电（K2），唯一可能的零序通路就是被测接地变。库源码 `pd2ppc_zero.py:317–318` 把 `"zd"` 同时列入 `_earthed_hv` 与 `_zigzag`，即按 HV 侧接地的 zigzag 处理；计划依命名约定（Z 后无 N＝中性点不引出）认为它不应提供接地通路，并在同节预期中写了「I3(b) 可能不满足」。
  - 归类 V（可能真实违反，或依赖命名约定），交裁定。
  - 最小可重放输入（harness 已同输入重放两次，结果相同）：分支 `Zd[max,inv=True]` 点 0，种子 [20261099, 4, 0]；`_t2_zz_net(p,"Zd")`：hv=77.87916664712574、mv=76.6399406245256、lv=59.60668975629773 kV；外网 s_sc=7315.297108441367 MVA、rx 0.1、x0x=2.809690248697072、r0x0 0.1；主变 YNd：m_sn=56.957303183261956、m_vk0=vk=5.449961667644061、m_vkr0=vkr=1.4176092705741878、m_si0=0.5085547788556974；接地变 Zd：sn=1.5118253020672532、vk0=vk=8.437929867394306、vkr0=vkr=0.32555153669828374、si0=0.8755107551018427；其余按 K2。`calc_sc(fault="1ph", case="max", inverse_y=True, lv_tol_percent=10)`，MV 母线（bus 1）ikss_ka=0.48939681475729857（rn=xn=0）；rn=8.290917827516807、xn=0.16478871303208797 时为 0.48501016070336156；期望 0。
- **注 S1 的边界**：负荷大小、基准、收敛容差原是计划未给、由实施者取的值（J7/J8）。本次按「实施取值使计划在全域不可执行」归 S 修正，没有改判据和容差。若主会话认为这属于改计划映射，须补记协议 §7。

## 8. v2 第三路（`spotcheck_v2.py`；计划＝`specs_codex_v2/e2nIEE_pandapower_3016/spec.md`「Numerical check method」＋同节 Invariants、Common operating domain）

读过：v2 spec.md 全文、`input/pr.diff`（只取 PR 测试夹具数值）、合并提交 `calc_sc.py` 签名与 `shortcircuit/results.py` 的 infinity 替换规则（compute node A 只读）。没读第三路 a 规格以外的审计材料，也没读 `compare3_*`。

分支（k=6，N 等分，各 ⌊N/6⌋ 或多 1）：`I1_Dyn`、`I1_YNd`、`I1_YNyn`、`I2_ZNyn`、`I2_ZNd`、`I2_ZNy`。必测点 12 个（每分支 2 个，不占 N）。种子 20261017+8×10+3=20261100；JSON 同 harness 输出，另加 `track:"3"`、`plan:"codex_v2"`。

| 计划原句 | 代码 |
|---|---|
| 用 create_empty_network/create_bus/create_ext_grid/create_transformer_from_parameters 建夹具，calc_sc，仅读 res_bus_sc 的 rk0/xk0/rk/xk/ikss；保存每次结果再运行下一种故障 | `solve()`→`spotcheck.sc()` 返回字典副本；只用这五列 |
| 求解选项 fault 1ph、case max、lv_tol 10、use_pre_fault_voltage False、r/x_fault 0、branch_results False；两种 inverse_y | `solve()` 显式传全部选项 |
| I1：两母线 Dyn（HV 放 ext_grid、观察 LV）、反向 YNd（LV 放 ext_grid、观察 HV）；Z_N=0,R,jX,R+jX；复数增量等式、正序不变、3ph/2ph 电流不变 | `fx_i1`、`make_ev_i1`：每个 inverse_y 跑四态；`I1.dZ0=3dZN`（R、jX、RX 三个增量）、`I1.Z1_unchanged`、`I1.{3ph,2ph}_ikss_unchanged`（Z_N=0 对 R+jX） |
| 同时改变正序 vk/vkr 使 K_T 变化，在每一固定参数组内比较中性点增量；不跨 K_T 比总 Z0 | 每点第二组 (vk_b, vkr_b)，零序参数不变，组内重复上一行全部检查，组间不比 |
| YNyn 补充夹具：两母线，仅 LV 有 ext_grid，HV 无其他元件；比 HV 的 Z0(R+jX)−Z0(0) 与 3(R+jX)；用 T 电路式作独立绝对预言（Z_{k,r}、K_T、Z_m 式；x0x=1、r0x0=rx=ρ；Z_g0,H=1.1U_H²(ρ+j)/(S_g√(1+ρ²))） | `fx_ynyn`、`ynyn_pred`、`ev_ynyn`：四态绝对值 `I1.YNyn_T_pred_Z0H`＋增量＋正序/3ph/2ph 不变 |
| 扫描有限磁化阻抗及 α | mag0、mag0_rx、si0 每点在域内均匀抽 |
| I2：110 kV ext_grid→110/20 kV YNd→20 kV bus→20/0.4 kV zigzag→aux（aux 无负载/电源）；ZNyn/ZNd/ZNy；z_k=K_TZ_{k,r}/Z_B,H、z_m=Z_m/Z_B,H 与 HV/LV Ω 预言比较；LV 乘 Z_B,L | `fx_zz`、`zz_pred`：HV 预言 αK_TZ_{k,r}+3Z_N；ZNyn 的 LV 预言 (Z_m+(1−α)K_TZ_{k,r})·(U_L/U_H)²；四态×两 inverse_y |
| 分别扰动 rn/xn、mag0_percent/mag0_rx、次级接线，检查独立性 | rn/xn：四态都对 HV 预言，ZNyn LV 对 Z_N=0 与 R+jX 两态比；磁化：换 (mag0_alt, mag0_rx_alt)，HV Z0 不变、ZNyn LV 对新预言；次级：同网改成另外两种 ZN*，HV Z0 不变 |
| aux 加 ext_grid，正序 s_sc/rx 固定，只变 x0x/r0x0，HV Z0 不变 | 两组 (x0x, r0x0) 各建一网，`I2.Z0H_indep_LV_zero_seq_terminal` |
| ZNd/ZNy 无源 aux 只查 I_f≈0，允许 Z0 为 infinity，不做 inf−inf | `lv_check`：ikss≤1e-6 kA；`cmp_bus` 两边同为无穷就跳过 Z0，一边有限一边无穷记失败 |
| I3：每个有限结果用 I1/I2 独立预言的 Z0 及公开 Z1 算 √3·1.1·U_b/\|2Z1+Z0\| 比 ikss；另用公开 Z0 算一次（内部一致性） | `i3_formula`：Dyn/YNd 的独立预言＝Z0(0)+3Z_N（Z_N=0 态没有独立预言，只做内部一致性）；YNyn 用 T 式；zigzag 用 HV 预言，ZNyn 的 LV 也用 LV 预言 |
| 固定 X_N 递增 R_N 扫描查单调性和上界 | `i3_sweep`：R_N∈{0, R/2, R, 2R, R_lim}（排序），单调不增、R_N>0 时 I_f≤cU/(√3R_N)；用于 I1 三个夹具（故障母线只有该变压器一个零序入口）与 zigzag HV |
| 可借助同网 3ph 结果，用 I1ph/I3ph=3\|Z1\|/\|2Z1+Z0\| 消去电压系数 | `i3_ratio`（Z_N=0、R+jX 两态） |
| 大 NER：t=\|D\|/(3R_N)<1 时 1/(1+t)≤I_f/(C/(3R_N))≤1；用 t≤0.01；保证未进开路截断 | R_lim=\|D\|/0.03（t=0.01）；I_f<1e-6 kA 记失败 |
| 每种夹具分别用 inverse_y True/False，并更改 S_B=1,10,100 后检查相同 Ω/kA | `aux_checks`：Z_N=0 与 R+jX 两态 × 两 inverse_y × S_B∈{10,100}，对 S_B=1 比；两 inverse_y 互比；比 Z0、Z1（复数模）与 ikss |
| 容差：阻抗 1e-8 Ω+1e-7\|Z*\|（复数模），电流 1e-8 kA+1e-7\|I*\|，开路 1e-6 kA | `ZR/ZA/IR/IA/IOPEN` |
| 增量选到 3\|ΔZ_N\|≥0.1\|Z_{k,earth}\|；域：3R_N/\|Z_{k,earth}\|、3X_N/\|Z_{k,earth}\|∈[0.1,100]，并含 PR 的 0.01、0.02、5 Ω | 比值各在 [0.1,100] 均匀抽；必测点 (R_N,X_N)=(0.01,0.02) 与 (5,0) |
| 可选导数：中心差分 h~ε^{1/3}max(R_N,\|D\|/3)，R_N>h，比较 h、h/2、2h 稳定后对解析导数 | `i3_sweep` 末段：三个步长差商之差 >1e-6 相对记失败（病态，不算违反），稳定就对 −3C(A+3R)/((A+3R)²+B²)^{3/2} 按默认容差判 |
| 被动、无旋转电机/电流源；额定绕组电压＝母线额定；无 tap、parallel=1、非 PSU、pfe=i0=0、50 Hz；shift 330°，YNyn 0° | `spotcheck.trafo`（pfe/i0=0、SHIFT 表：dyn/ynd/zn* 330°、ynyn 不在表中即 0°） |

实施判断（v2）：
- M1 分支只取计划写到的六个夹具/联结组。inverse_y 与 S_B 计划要求同一夹具两种都跑、几种基准都跑再互比，所以作为点内交叉检查，不拆成分支。case 只有 max（计划求解选项）。
- M2 均匀抽样：U_H∈[10,110]、U_L∈[0.4,33] kV，要求 U_H>U_L（不满足重抽）；S_T∈[0.4,40] MVA（所有变压器）；vk、vk0∈[4,15]，vkr/vk、vkr0/vk0∈[0.02,0.4]；α∈[0.1,0.9]；mag0∈[10,1000]%；mag0_rx∈[0,1]；外网 s_sc∈[100,1e4]、rx∈[0.05,0.5]、x0x∈[0.5,3]、r0x0∈[0.05,0.5]，min 参数取与 max 相同的正值。计划没写对数均匀，故不用对数。
- M3 \|Z_{k,earth}\| 取接地侧电压下的原始零序短路阻抗模 vk0/100·U²/S_T：Dyn 用 U_L，YNd、YNyn、zigzag 用 U_H。
- M4 I2 夹具的电压：计划写 110/20/0.4 kV，域给 U_H、U_L 区间且「主变源端可用 110 kV」。按协议 §2，源端固定 110 kV、主变 110/U_H，U_H 与 U_L 在域内抽；110/20/0.4 作为必测点。
- M5 K_T 一律按计划式取 c=1.1（case max、lv_tol 10%，计划写明各电压级 c=1.1），Z_m 不乘 K_T（计划明写沿用库定义）。
- M6 Dyn/YNd 没有绝对 Z0 式，I3 的「独立预言」取 I1 预言 Z0(0)+3Z_N。
- M7 R_N 扫描点 {0, R/2, R, 2R, R_lim}：计划只写「递增 R_N 扫描」，点数和取值由实施者定。A 或 B<0 时计划前提不满足，只记录不判。
- M8 导数检查计划写为「可选」；协议 §2 规定「建议」类也算，故实现。稳定性与解析比较都没给容差，用协议默认（相对 1e-6、绝对 1e-12）；不稳定记失败、不算违反（计划写「遇病态先单独记录」）。
- M9 必测点：Dyn 与 YNyn 用 PR 的 20/0.4 kV、0.63 MVA、6%/1%、mag0 100、mag0_rx 0、si0 0.9、外网 1000 MVA；YNyn 外网在 LV，ρ=0.1。YNd 用 PR zigzag 测试的 110/20 kV、40 MVA、12%/0.4% 主变参数、外网 5000 MVA 接 LV。zigzag 用 PR 测试原夹具。下列是自取值（计划没给）：K_T 第二组 (4%,1%)，YNd 为 (6%,1%)；磁化扰动值 (1000%, 0.5)；aux 外网 1000 MVA、rx 0.1、(x0x,r0x0)=(1,0.1)/(3,0.5)。
- M10 aux 外网检查只比两组 x0x/r0x0 之间的 HV Z0，不和无 aux 外网时比：加外网会改变正序，计划只要求「保持正序不变、只改零序终端」时 HV Z0 不变。
- M11 3ph/2ph 不变只比 Z_N=0 与 R+jX 两态的 ikss（计划写「电流不变」）。

耗时与内存（推测）：首次冒烟里 2–3 母线网单次 calc_sc 约 10 ms（第三路 a 的 Dyn 分支约 23 次求解 0.235 s/点）。v2 每点求解次数：I1_Dyn/YNd 约 51 次＋2 次建网，I1_YNyn 约 35 次，I2 约 34 次＋3 次建网（3 母线）。单点约 0.35–0.55 s，均值约 0.45 s；N=100（另加 12 个必测点）约 50 s，N=10⁴ 约 1.3 h（上限估 2.5 h），以冒烟二实测为准。峰值内存同第一至三路，单进程约 300–400 MB（pandas/scipy＋小网）。全量四进程并行（第一、二路、v2、第三路 a）合计约 1.6 GB，JobSpec 申报 3 GiB。

JobSpec：`job_smoke2.json`（修过的第一、二路与 v2 各 N=100，`;` 串行，输出 `dq1/out/smoke2/`）；`job_full.json`（重写：四路各 N=10⁴ 并行后 wait，输出 `dq1/out/full/`，第三路 a 输出名 `pp3016_t3a.json`）。原 `job_smoke.json` 保留未动。
