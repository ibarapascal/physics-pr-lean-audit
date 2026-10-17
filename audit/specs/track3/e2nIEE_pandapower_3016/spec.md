## e2nIEE/pandapower#3016 Neutral earthing impedance (rn_ohm) and zigzag earthing transformers in the zero-sequence short circuit

- What changed (1–3 sentences, naming files and formulas):

  `pandapower/create/trafo_create.py` 为四个单个/批量变压器创建 API 增加 `rn_ohm`；`pandapower/pd2ppc_zero.py` 将中性点阻抗记为 (Z_N=R_N+jX_N)，普通接地绕组加入 (3Z_N)，且置于 IEC 变压器修正之后。新增 zigzag 分支以 HV 对地支路 (α Z_{0k}+3Z_N) 和（仅 `ZNyn`）LV 对地支路 (Z_{0m}+(1-\alpha)Z_{0k}) 表示两侧零序解耦；测试位于 `pandapower/test/shortcircuit/test_1ph.py`，说明位于 `CHANGELOG.rst`。注意输入 PR 描述中的 `ZN` 与合并源码不一致：实际分支接受大小写不敏感的 `znyn/znd/zny/zd`，没有 `zn`；下列物理检查使用明确的 `ZNyn/ZNd/ZNy`，不把名称支持问题混入物理结论。

- Invariants (mathematical statements; write them as checkable equalities/inequalities, give variables, sign conventions, domain of validity):

  统一约定：线性、工频、对称分量短路等值电路，RMS 相量，(j^2=-1)，被动阻抗 (Z=R+jX)，感性 (X\ge0)，端口电流以流入元件为正。(U_b) 是母线线间额定电压（kV），(S_B) 是任意计算基准（MVA），(Z_{B,b}=U_b^2/S_B)（Ω）；`rn_ohm/xn_ohm` 是接地绕组中性点本侧的实际 Ω，不是标幺值。(Z_{q,b}) 表示母线 b 的 q 序驱动点戴维南阻抗，(I_f) 是单相接地初始对称短路电流幅值（kA）。所有比较均保持拓扑、正序设备参数和电压系数不变，故障阻抗为零，无电流源型逆变器、发电机或电容对地支路。

  **I1 — primary：中性点三倍阻抗必须处于正确的绕组端口，且不改变正/负序。**

  由 (I_N=I_a+I_b+I_c=3I_0)，中性点压降为 (3Z_NI_0)，故该绕组零序端口增加的串联阻抗严格为 (3Z_N)，正、负序增加量为零。对于故障点在接地绕组母线、该绕组是该母线唯一零序返回入口的隔离夹具：

  \[
  Z_{0,b}(Z_N+\Delta Z_N)-Z_{0,b}(Z_N)=3\Delta Z_N,
  \quad Z_{1,b}(Z_N+\Delta Z_N)=Z_{1,b}(Z_N),
  \quad Z_{2,b}(Z_N+\Delta Z_N)=Z_{2,b}(Z_N).
  \]

  这意味着电阻扰动只改变零序实部，电抗扰动只改变虚部；系数不是 (3K_T)，也不是 (3\alpha)。先以 `Dyn` 的 LV 接地端、`YNd` 的 HV 接地端检验；不要把此驱动点等式推广到有多个并联接地路径的任意母线。换侧折算时必须遵守 (Z'=a^2Z)，其中 (a=U_H/U_L)；标幺中性点增量是 (3\Delta Z_N/Z_{B,\mathrm{earth}})。

  `YNyn` 的补充端口检验也属于 I1：按合并源码的 HV 侧选择语义，令 `rn_ohm/xn_ohm` 指 HV 中性点阻抗，LV 中性点直接接地。全部阻抗折算到 HV，设原始总零序漏阻抗 (Z_{k,r})、修正后 (Z_k=K_TZ_{k,r})、磁化阻抗 (Z_m)，以及 (α=\texttt{si0_hv_partial})。在 LV 外部电网供电、HV 为无其他接地路径的故障端时，应该有

  \[
  Z_{0,H}=\alpha Z_k+3Z_N+
    \frac{Z_m\big[(1-\alpha)Z_k+Z_{g0,H}\big]}
         {Z_m+(1-\alpha)Z_k+Z_{g0,H}}.
  \]

  因而对有限 (Z_m,Z_{g0,H}) 仍必须满足 (\Delta Z_{0,H}=3\Delta Z_N)。该式表示中性点阻抗加在 HV 支路之外，不能先加进 (Z_k) 再拆分到 HV/LV；这是在声明 HV 中性点语义后的物理要求，不是声称本 API 能独立指定两个中性点。

  **I2 — primary：zigzag 的零序端口解耦及其各自的接地阻抗。**

  对本 PR 声明的理想 zigzag 等值模型，采用共同功率基准、各自绕组电压基准，令 (z_k,z_m) 为相应标幺阻抗，(α\in(0,1))。应有

  \[
  \begin{bmatrix}i_{0,H}\\i_{0,L}\end{bmatrix}=
  \begin{bmatrix}1/z_H&0\\0&y_L\end{bmatrix}
  \begin{bmatrix}v_{0,H}\\v_{0,L}\end{bmatrix},\qquad
  z_H=\alpha z_k+3Z_N/Z_{B,H},
  \]
  \[
  y_L=\begin{cases}
    1/[z_m+(1-\alpha)z_k],&\texttt{ZNyn},\\
    0,&\texttt{ZNd},\texttt{ZNy}.
  \end{cases}
  \]

  零序交叉导纳 (Y_{0,HL}=Y_{0,LH}=0) 是核心结构约束。在主变 LV 为 delta、该系统仅由 zigzag 提供接地、aux 母线无其他设备的三母线夹具中，(Z_{0,H}=Z_{B,H}z_H)，`ZNyn` 时 (Z_{0,L}=Z_{B,L}[z_m+(1-\alpha)z_k])；`ZNd/ZNy` 的 aux 母线理想接地故障电流为零。HV 的 (Z_0) 不随 secondary 的 `yn/d/y` 类型或磁化阻抗变化；`ZNyn` 的 LV (Z_0) 不随 HV 中性点阻抗变化。若以外部电网改变 LV 的零序接地强度且保持正序强度不变，HV (Z_0) 仍应不变。此处不要求任意接线变更后的总故障电流相同，因为其正序网络也可能改变。

  其中“交叉导纳为零”是理想 zigzag 安匝抵消模型的结构要求；把输入 `vk0_percent` 的哪一部分解释为 zigzag 自身漏阻抗，则依赖该 PR 的参数化约定。不能仅凭代码采用了 (α z_k) 就宣称实际设备也符合这一参数映射，需单独保留其文献/实测验证边界。

  **I3 — secondary：三序串联故障公式与 NER 限流。**

  对上述无旋转电机、电流源的被动夹具，(Z_2=Z_1)，单相金属性接地故障满足

  \[
  I_0=I_1=I_2=I_f^{\mathrm{phasor}}/3,\qquad
  I_f=\frac{\sqrt3\,cU_b}{|2Z_{1,b}+Z_{0,b}|}.
  \]

  在 I1 的唯一返回路径夹具或 I2 的 zigzag HV 夹具中，固定 (X_N\ge0)，定义 (D=2Z_{1,b}+Z_{0,b}(R_N=0)=A+jB)，假定 (A,B\ge0)，则

  \[
  I_f(R_N)=\frac{\sqrt3\,cU_b}{\sqrt{(A+3R_N)^2+B^2}},\quad
  \frac{dI_f}{dR_N}=-\frac{3\sqrt3\,cU_b(A+3R_N)}{[(A+3R_N)^2+B^2]^{3/2}}\le0.
  \]

  对 (R_N>0)，(I_f\le cU_b/(\sqrt3R_N))，且 (R_NI_f\to cU_b/\sqrt3)。该单调性不推广到含谐振电容、并联接地回路或电流源的任意网络。

- Sources of the invariants (first principles / literature formula / consistency with other implementations in the library / dimensions):

  I1 来自 Kirchhoff 电流定律、对称分量定义与中性点 Ohm 定律；正确支路位置来自上述端口电路。I2 的解耦来自理想 zigzag 绕组的零序安匝抵消，其数值漏阻抗拆分是本 PR 声明的等值模型。I3 来自故障条件 (I_b=I_c=0,V_a=0) 下三序串联以及 (E_1=cU_b/\sqrt3)；限流不等式和导数直接由这一公式推导。单位关系 (\mathrm{kV}^2/\mathrm{MVA}=\Omega)、(\mathrm{kV}/\Omega=\mathrm{kA}) 为独立量纲检查。

  合并源码证据（均固定在本次 merge SHA）：[零序建模及磁化/漏阻抗定义](https://github.com/e2nIEE/pandapower/blob/0c6a1f0a65273af91dd1d5b854498cd89d457939/pandapower/pd2ppc_zero.py)、[IEC 修正因子实现](https://github.com/e2nIEE/pandapower/blob/0c6a1f0a65273af91dd1d5b854498cd89d457939/pandapower/build_branch.py)、[单相短路公式](https://github.com/e2nIEE/pandapower/blob/0c6a1f0a65273af91dd1d5b854498cd89d457939/pandapower/shortcircuit/currents.py)、[公开结果列提取](https://github.com/e2nIEE/pandapower/blob/0c6a1f0a65273af91dd1d5b854498cd89d457939/pandapower/shortcircuit/results.py)。这些是 API/参数语义及库内一致性证据，不充当独立物理正确性证明。

  输入 [PR #3016](https://github.com/e2nIEE/pandapower/pull/3016) 的正文及源码注释提及 IEC 60909、PowerFactory transformer technical reference Fig. 8.3、J&P Transformer Book Appendix 5；本次没有独立取得这些文献的相应版本和图式，因此其精确 zigzag 参数映射、IEC 条款符合性及 PR 所称 PowerFactory 实网验证均标为未验证。没有读取评论、issue、fork 或另一审计轨道。

- Common operating domain (range of each input, with the basis):

  采用有限、线性、无饱和的短路等值模型；不覆盖电弧、暂态峰值、热限值、三相不平衡潮流或实际铁芯结构推断。下述范围是用于稳定随机抽查的工程尺度与输入合法性限制，不是声称 IEC 对每个参数规定了这些上下限。

  | 输入 | 抽查范围/固定值 | 依据 |
  |---|---|---|
  | 电压 | (U_H=10\ldots110) kV，(U_L=0.4\ldots33) kV，(U_H>U_L)；主变源端可用 110 kV | 覆盖 PR 的 20/0.4、110/20 kV 夹具 |
  | 设备额定功率 | (S_T=0.4\ldots40) MVA | PR 小型接地变至主变尺度 |
  | 任意网络功率基准 | (S_B\in\{1,10,100\}) MVA | 更换计算基准不应改变 Ω/kA 物理结果 |
  | 正/零序短路参数 | `vk_percent`, `vk0_percent` 各 4…15；`vkr/vk`、`vkr0/vk0` 各 0.02…0.4 | 非零、被动且平方根严格为正；避免零参数回退逻辑 |
  | 零序拆分/磁化 | (α=0.1\ldots0.9)，`mag0_percent`=10…1000，`mag0_rx`=0…1 | 内部拆分非退化；测试覆盖值为 0.9、100、0 |
  | 中性点阻抗 | (R_N,X_N\ge0)；基线 0，扰动令 (3R_N/|Z_{k,\mathrm{earth}}|)、(3X_N/|Z_{k,\mathrm{earth}}|\in[0.1,100])，并含 PR 的 0.01 Ω、0.02 Ω、5 Ω | 覆盖轻微变化到 NER 主导，避免极端病态尺度 |
  | 外部电网 | `s_sc_max_mva`=100…10000 MVA，`rx_max`=0.05…0.5，`x0x_max`=0.5…3，`r0x0_max`=0.05…0.5；需要 min 时给对应正值 | 有限被动电源，覆盖 PR 的 1000/5000 MVA |
  | 固定建模条件 | 额定绕组电压等于接入母线额定值；无 tap、`parallel=1`、`in_service=True`、非 power-station unit；`pfe_kw=i0_percent=0`；50 Hz | 排除变比、并联中性点语义和电站修正的额外因素 |
  | 求解选项 | `fault="1ph"`, `case="max"`, `lv_tol_percent=10`, `use_pre_fault_voltage=False`, `r_fault_ohm=x_fault_ohm=0`, `branch_results=False`；两种 `inverse_y` | 已读公开 API 的标准短路路径；此组合各电压级 (c=1.1) |

  `shift_degree` 分别单独指定，vector group 不带时钟数字；可沿用新增测试的 330°，`YNyn` 夹具用 0°。本规范不推定 `parallel>1` 表示共用或独立 NER，也不覆盖 power-station-unit 复合修正。理想开路检查单独处理，不要求其结果为有限复阻抗。

- Numerical check method (which public API of the merged commit, what to compute, what to compare against, suggested tolerance and why; if finite differences are needed, the step-size strategy):

  **入口与观测量。** 用 `pandapower.create_empty_network`、`create_bus`、`create_ext_grid`、`create_transformer_from_parameters` 建立下述夹具，调用 `pandapower.shortcircuit.calc_sc`，仅读取公开 `net.res_bus_sc` 的 `rk0_ohm/xk0_ohm/rk_ohm/xk_ohm/ikss_ka`。令 (Z_0=\texttt{rk0_ohm}+j\texttt{xk0_ohm})、(Z_1=\texttt{rk_ohm}+j\texttt{xk_ohm})；不调用私有转换函数，不以 `_ppc0` 作为验收入口。保存每次结果后再运行下一种故障，防止结果表被覆盖。[API 签名及默认模型](https://github.com/e2nIEE/pandapower/blob/0c6a1f0a65273af91dd1d5b854498cd89d457939/pandapower/shortcircuit/calc_sc.py)。

  **I1。** 两母线 `Dyn` 夹具在 HV 放 ext_grid、观察 LV；反向夹具 `YNd` 在 LV 放 ext_grid、观察 HV。分别运行 (Z_N=0,R,jX,R+jX)，直接检查复数增量等式、正序不变，以及 `fault="3ph"`、`fault="2ph"` 的电流不变（负序影响的间接公开观测）。同时改变正序 `vk/vkr` 使 (K_T) 变化，再在每一固定参数组内比较中性点增量，验证 NER 不受 (K_T) 缩放。不要跨不同 (K_T) 直接比较总 (Z_0)。

  `YNyn` 补充夹具：两母线、仅 LV 有 ext_grid、HV 无其他元件，比较 HV 的 (Z_0(R+jX)-Z_0(0)) 与 (3(R+jX))，并用 I1 给出的 T 电路公式作独立绝对值预言。设 (p_0=\texttt{vk0_percent},r_0=\texttt{vkr0_percent},p=\texttt{vk_percent},r=\texttt{vkr_percent},q=\texttt{mag0_rx},m=\texttt{mag0_percent}/100)，则在 HV Ω 下

  \[
  Z_{k,r}=\frac{U_H^2}{S_T}\frac{r_0+j\sqrt{p_0^2-r_0^2}}{100},\quad
  K_T=\frac{0.95(1.1)}{1+0.6\sqrt{p^2-r^2}/100},\quad
  Z_m=m|Z_{k,r}|\frac{q+j}{\sqrt{1+q^2}}.
  \]

  为简化源阻抗预言，取 `x0x_max=1`、`r0x0_max=rx_max=ρ`，则 LV 电网折至 HV 的 (Z_{g0,H}=1.1U_H^2(\rho+j)/(S_g\sqrt{1+\rho^2}))。扫描有限磁化阻抗及 (α)，避免用无限磁化近似掩盖中性点放错支路。这里 (Z_m) 沿用库现有未乘 (K_T) 的定义，审计目标是中性点位置而不是另行裁定磁化修正规则。

  **I2。** 使用新增测试的三母线结构：110 kV ext_grid → 110/20 kV `YNd` 主变 → 20 kV bus → 20/0.4 kV 待测 zigzag → aux bus；aux 不接负载/电源。给所有变压器完整的正/零序参数；测试 `ZNyn/ZNd/ZNy`。用以上 (Z_{k,r},K_T,Z_m) 公式但换成接地变参数计算 (z_k=K_TZ_{k,r}/Z_{B,H})、(z_m=Z_m/Z_{B,H})，与 I2 的 HV/LV Ω 预言比较；主变 delta 阻断源侧零序，所以 HV 预言不含主电网零序阻抗。注意 LV Ω 结果须乘 (Z_{B,L})，不能直接用 HV Ω。

  保持漏阻抗不变，分别扰动 `rn_ohm/xn_ohm`、`mag0_percent/mag0_rx` 和次级接线，检查 I2 的独立性。再在 aux 增加 ext_grid，保持其正序 `s_sc_max_mva/rx_max` 固定而只变化 `x0x_max/r0x0_max`，检查 HV (Z_0) 不变，以实际改变 LV 零序终端检验解耦。对 `ZNd/ZNy` 的无源 aux 开路只检查 (I_f\approx0)，并允许 `rk0_ohm/xk0_ohm` 为 infinity；合并实现用大阻抗模拟开路，结果提取也会把很大阻抗转换成 infinity，不能对此做 `inf-inf`。

  **I3。** 每个有限结果用 I1/I2 独立预言的 (Z_0) 及公开结果 (Z_1) 计算 (\sqrt3(1.1)U_b/|2Z_1+Z_0|)，对比 `ikss_ka`；另用公开 (Z_0) 计算一次，只作为结果内部一致性检查，不能替代独立预言。对固定 (X_N) 做递增 (R_N) 扫描检查单调性和上界；可以借助同网 `fault="3ph"` 的结果，用 (I_{1ph}/I_{3ph}=3|Z_1|/|2Z_1+Z_0|) 消去电压系数。

  大 NER 极限不能以任意固定“足够大”电阻判断：令 (C=\sqrt3cU_b)、(t=|D|/(3R_N)<1)，则应有 
  \(1/(1+t)\le I_f/(C/(3R_N))\le1\)。用 (t\le0.01) 给出不超过约 1% 的理论偏差，并保证结果未进入开路截断范围。

  每种夹具分别使用 `inverse_y=True/False`，并更改 (S_B=1,10,100) 后检查相同物理 Ω/kA 结果；这些是辅助数值一致性检查，不是独立物理预言。核心等式使用复数模残差 ( |\hat Z-Z_*|\le10^{-8}\,\Omega+10^{-7}|Z_*|\)，电流使用 ( |\hat I-I_*|\le10^{-8}\,\mathrm{kA}+10^{-7}|I_*|\)；近零开路电流阈值可用新增测试的 (10^{-6}\) kA。依据是小规模、非退化双精度线性电路应远优于工程百分比误差；遇病态/开路结果先单独记录求解失败与尺度，不盲目放宽阈值。增量选到 (3|\Delta Z_N|\ge0.1|Z_{k,\mathrm{earth}}|)，避免大数相减淹没信号。不需要有限差分：增量关系和限流公式均可直接检验；可选导数检查用中心差分 (h\sim\epsilon^{1/3}\max(R_N,|D|/3))，要求 (R_N>h)，比较 (h,h/2,2h) 稳定性后再与解析导数比较。

- Do the PR's own tests already check this invariant (yes/partial/no, naming test files and assertions):

  **I1：partial。** [test_1ph.py](https://github.com/e2nIEE/pandapower/blob/0c6a1f0a65273af91dd1d5b854498cd89d457939/pandapower/test/shortcircuit/test_1ph.py) 的 `test_trafo_neutral_earthing_impedance` 在一个 `Dyn`、`parallel=1`、网络基准 1 MVA 的夹具中断言 (\Delta R_0=0.03\) Ω、(\Delta X_0=0.06\) Ω，另一分量不变；电阻案例正序 R/X 不变，`atol=1e-6, rtol=0`。两种求解路径均参数化；未覆盖 HV 接地组、`YNyn` 支路位置、混合 R+jX、不同修正因子与计算基准。

  **I2：partial。** 同文件 `test_zigzag_earthing_transformer` 检查 `ZNyn` HV 电流为正、与 `YNd` 的比在 0.5…2、`ZNyn/ZNd` HV 电流相等、`ZNyn` LV 电流为正及 `ZNd` LV 电流约零。没有核验 (α z_k) 和 LV 磁化支路的精确阻抗，没有测试 `ZNy`、磁化扰动、LV 终端扰动或 HV NER 对 LV 的独立性；“同数量级”不能确认漏阻抗公式正确。

  **I3：partial。** 两个新增测试只断言加入指定电阻后故障电流下降；未直接检验三序绝对值公式、连续单调性、解析上界或 NER 主导极限。以上是读取断言所得，不表示本次执行了测试，也未将 PR 正文“既有测试全通过”视作当前运行证据。

- Do you expect the change to satisfy it (satisfies / may not satisfy / uncertain), with one or two sentences of reasoning (speculation, not a conclusion):

  **I1：may not satisfy。** 推测 `Dyn/YNd` 的标准短路路径能满足三倍增量，但 `YNyn` 把 (3Z_N) 加进总漏阻抗后再按 (α,1-\alpha) 拆分，有限磁化支路下可能违背 HV 中性点应完整串入 HV 支路的端口公式；尚未数值验证。

  **I2：uncertain。** 推测额定变比、单台、非电站夹具可满足声明的对地支路及解耦结构，因为实现保留极大 HV–LV 阻抗并分别加入母线导纳；实际 zigzag 漏阻抗能否由这些输入唯一正确表征，仍缺独立文献/设备数据确认。

  **I3：satisfies（限已正确实现 I1/I2 的夹具，推测）。** 合并源码明确用 (2Z_1+Z_0) 计算电流，唯一被动接地路径应具有上述限流行为；若上游零序网络有错，电流公式内部一致并不能证明物理结果正确。

- Environment: merge commit SHA, language and dependencies needed, estimated setup difficulty (low/medium/high):

  Merge SHA：`0c6a1f0a65273af91dd1d5b854498cd89d457939`。Python；该提交 [pyproject.toml](https://github.com/e2nIEE/pandapower/blob/0c6a1f0a65273af91dd1d5b854498cd89d457939/pyproject.toml) 要求 Python >=3.10，运行依赖为 `numpy>=1.26,<2.4`、`scipy<1.17`、`pandas~=2.3`、`networkx~=3.4`、`packaging~=25.0`、`tqdm~=4.67`、`deepdiff~=8.6`、`geojson~=3.2`、`typing_extensions~=4.9`、`pandera~=0.26.1`。重复原测试还需 pytest 与其配置使用的 pytest-timeout；不需要 PowerFactory、GPU 或实际设备。预计 **low**（常见 CPU 环境、数值库 wheel），具体平台兼容性未验证。本次只读取输入和 GitHub 固定提交源码、撰写规范；未安装、编译或运行任何项目代码。

- Formalization (Lean) feasibility: can the invariant be written as a proposition over finite-dimensional reals, and what are the main obstacles:

  可以。用复数等价的二维实向量表示相量，用有限维实矩阵表示节点导纳；I1 是 (I_N=3I_0) 和串联阻抗的代数恒等式，`YNyn` 补充式是有限电路消元；I2 是块对角矩阵与两端口导纳恒等式；I3 的平方电流公式、单调不等式和有限误差界可归结为有理式、平方根及正分母假设。

  必须显式声明被动性、非零分母、可逆/接地条件、理想变压器电压折算和 zigzag 安匝假设。主要障碍是证明真实绕组到所声明等值电路的映射、明确 IEC 修正与输入参数含义，以及将该数学模型对应到 pandas/NumPy/SciPy 的浮点与稀疏求解实现；代码注释不足以证明这些前提。开路应以零导纳或极限表示，不将 IEEE infinity 当作实数；绝对数学命题与软件容差/大阻抗近似的误差证明需分开。
