## idaholab/moose#34020 Pass the nodal slip rate to dynamic mortar function_friction

- What changed (1–3 sentences, naming files and formulas):

  `modules/contact/src/constraints/ComputeDynamicFrictionalForceLMMechanicalContact.C` 将传给速度相关摩擦函数的旧时间步速度由不含积分权重的求和改为节点投影速度 \(\mathbf s_i=\mathbf W_i/D_i\)，只在自由度所属进程保存历史，并以零向量处理没有历史的自由度；常数参数 `mu` 新增 `mu > 0` 校验。`ComputeDynamicWeightedGapLMMechanicalContact.C` 无论 `normalize_c` 是否开启，都累加并通信分母 \(D_i\)；对应头文件、对象文档及 `mortar_dynamics/gold/frictional-mortar-3d-dynamics-light-function_out.e` 同步更新。以下以合并提交源码为准，所有数值检查均为待执行方案，本次未运行计算。

- Invariants (mathematical statements; write them as checkable equalities/inequalities, give variables, sign conventions, domain of validity):

  **公共符号及边界。** 在二维或三维接触中，取 \(k=1\) 或 \(2\) 个正交单位节点切向量 \(\mathbf t_{ia}\)。定义相对速度 \(\mathbf v_r=\mathbf v_{secondary}-\mathbf v_{primary}\)，压应力乘子 \(p_i>0\)，切向乘子向量 \(\boldsymbol\lambda_i\in\mathbb R^k\)；本库的切向乘子与滑移同向，作用在 secondary 上的物理摩擦牵引取负号。令 \(\psi_i\) 为实际使用的 dual mortar 测试函数，\(w_q=J_q\omega_q\,coord_q\) 为同一 mortar 接触积分域上的几何、求积及坐标因子，定义
  \[
  D_i^n=\sum_q\psi_i(q)w_q^n>0,\qquad
  W_{ia}^n=\sum_q\psi_i(q)w_q^n\,\mathbf v_r^n(q)\cdot\mathbf t_{ia}^n,
  \qquad s_{ia}^n=W_{ia}^n/D_i^n.
  \]
  仅对有非退化接触支撑、未被 edge dropping 去除的自由度提出这些断言。dual 测试函数可能有负值，不能据此要求投影速度总在所有求积点速度的最小值与最大值之间。

  **I1 — Primary：投影速度的常量再现、积分尺度及分区不变性。** 实际用于摩擦函数的节点历史速度必须满足上述 \(D_i\mathbf s_i=\mathbf W_i\)。在平面接触、切向基固定且 \(\mathbf v_r(q)\cdot\mathbf t_{ia}=U_a\) 时，必须有 \(s_{ia}=U_a\)，与节点相邻单元数、接触面面积、求积点数无关（前提是相应积分准确且物理速度场相同）。对共同的正因子 \(\alpha\)，\(\mathbf s_i(\alpha\mathbf W_i,\alpha D_i)=\mathbf s_i(\mathbf W_i,D_i)\)；任意不重不漏的进程分区满足 \(\mathbf s_i=(\sum_r\mathbf W_{ir})/(\sum_r D_{ir})\)，不能平均各进程的比值。该物理输入不应取决于 `normalize_c=true/false`。这些是常量场的精确再现关系，不声称任意非均匀解在不同网格上逐节点相等。

  **I2 — Secondary：各向同性速度范数和明确的时间滞后。** 定义 \(\delta=10^{-12}\)（当前速度单位下的源码正则化常数），在第 \(n\) 个时间步内：
  \[
  \mu_i^n=F\left(0,\ p_i^{n-1},\ \sqrt{\|\mathbf s_i^{n-1}\|^2+\delta^2},\ 0\right).
  \]
  这里函数的四个位置依次为时间参数、空间参数 `x`、`y`、`z`，即 `x` 被重用为旧接触压力，`y` 被重用为旧滑移速率，时间和 `z` 都为零。在没有速度历史时令 \(\mathbf s_i^{n-1}=0\)，所以函数收到的是 \(\delta\)，并非精确的零；压力仍来自 `solutionOld()`，不强制为零。同一接受状态后的重试不得改用被拒绝试算的速度历史。对任何 \(Q\in O(k)\)，\(\mu(p,Q\mathbf s)=\mu(p,\mathbf s)\)；二维速度反号或三维切向基旋转均不改变系数。无历史时置零是离散初始化约定，不能解释为一个刚进入界面的物体在物理上必然静止。

  **I3 — Secondary：使用正确历史系数的库仑界及非正摩擦功率。** 限定当前持续闭合、真实加权法向间隙为零、\(p_i^n>\epsilon\)、\(\mu_i^n\ge0\)、\(c_t>0\)、\(\Delta t>0\)。在收敛解上应满足
  \[
  \|\boldsymbol\lambda_i^n\|\le\mu_i^n p_i^n;
  \quad \|\mathbf s_i^n\|>0\ \Longrightarrow\
  \boldsymbol\lambda_i^n=\mu_i^n p_i^n\frac{\mathbf s_i^n}{\|\mathbf s_i^n\|};
  \quad \mathcal P_i=-\boldsymbol\lambda_i^n\cdot\mathbf W_i^n
  =-D_i^n\mu_i^n p_i^n\|\mathbf s_i^n\|\le0.
  \]
  在 sticking 时 \(\mathbf s_i=0\)，只要求圆盘界和零滑移功率。该功率是 dual mortar 离散配对意义下的节点贡献，不是任意插值点的点态断言。对正常分离、\(p_i<\epsilon\) 的分支，源码残差要求 \(\boldsymbol\lambda_i=0\)。动态法向 persistency 允许真实间隙有残差，故一般动态状态不能无条件使用半径 \(\mu p\)：源码实际半径为 \(\mu\max(0,p-c_{eff}\bar g)\)；上述物理界测试特意选 \(\bar g=0\) 或量化其误差。

  常数 `mu <= 0` 应被参数校验拒绝，是此 PR 的附加 API 合同；零摩擦本身并非违反物理。任意 `function_friction` 的输出没有在本次改动中被钳制或强制为正，故不能提出“所有允许函数都耗散”的无条件不变量。

- Sources of the invariants (first principles / literature formula / consistency with other implementations in the library / dimensions):

  I1 来源于投影对常量的再现以及积分的线性性：\(\sum\psi wU=U\sum\psi w\)。\(D\) 是接触测度，\(W\) 是接触测度乘速度，因此 \(W/D\) 才具有物理速度的尺度；仅靠量纲不能排除原来的无权求和错误，必须检查常量再现。合并版 [动态约束实现](https://github.com/idaholab/moose/blob/f90f508580f7611f7d75ee458cd177affc7288d1/modules/contact/src/constraints/ComputeDynamicFrictionalForceLMMechanicalContact.C) 的 `timestepSetup()` 及 [加权法向约束](https://github.com/idaholab/moose/blob/f90f508580f7611f7d75ee458cd177affc7288d1/modules/contact/src/constraints/ComputeDynamicWeightedGapLMMechanicalContact.C) 给出积分、符号和历史读取。[准静态约束](https://github.com/idaholab/moose/blob/f90f508580f7611f7d75ee458cd177affc7288d1/modules/contact/src/constraints/ComputeFrictionalForceLMMechanicalContact.C) 也使用该比值，但用的是当前状态，不能直接要求两种完整瞬态解一致。源码引用 Wohlmuth (2011), Eq. 3.10a；本次未独立核对论文原式，证明依据采用上述直接代数关系，不把源码的文献归属当作独立验证。

  I2 的范数不变性来自各向同性标量摩擦定律；使用上一接受时间步是该公开对象的离散合同，而不是连续物理定律。I3 来自库仑摩擦的最大耗散条件；[MortarContactUtils.h](https://github.com/idaholab/moose/blob/f90f508580f7611f7d75ee458cd177affc7288d1/modules/contact/include/utils/MortarContactUtils.h) 中公开的球投影残差给出独立可检查的代数表述，[TangentialMortarMechanicalContact.C](https://github.com/idaholab/moose/blob/f90f508580f7611f7d75ee458cd177affc7288d1/modules/contact/src/constraints/TangentialMortarMechanicalContact.C) 明确给出物理牵引与乘子相反的符号约定。

- Common operating domain (range of each input, with the basis):

  数学有效域为：二维/三维、有限速度与压力、dual LM 插值、正且非微小的 \(D_i\)、非退化网格、正时间步，以及 I3 所列的接触与系数限制。使用 `use_displaced_mesh=true`，平面固定法向和固定切向基；先不加入 wear、热依赖、旋转表面、重网格或恢复运行。源码要求 `newmark_beta > 0`、`newmark_gamma >= 0`，建议取测试已有的 \(\beta=0.25,\gamma=0.5\)，且动力学和约束的参数一致。

  建议随机抽样域（工程测试尺度，不是作者声明的全部有效范围）：以选定单位 \(L_0,T_0,P_0\) 无量纲化，面内长度 \(0.1\!\sim\!10L_0\)，每方向 1、2、4 个单元并采用非匹配的对面剖分，\(p/P_0\in[0.1,10]\)，\(\|\mathbf s\|/(L_0/T_0)\in[10^{-3},10^{-1}]\)，另加精确零速度，所有切向角度与正负号；\(\Delta t/T_0\in[10^{-3},2\times10^{-2}]\)，保持每步位移远小于接触支撑宽度。取 \(\epsilon/P_0\le10^{-7}\)，正 \(c,c_t\) 按残差量纲缩放，`normalize_c` 两种取值都测；常数 \(\mu\in[0.05,1]\)。默认正则化是硬编码数值，换速度单位时须同步解释 \(\delta\)，不宣称极小速度下严格单位不变。

  诊断函数可取 `0.2 + 0.02*x/P0 + 0.5*y/V0`（用实际数值替换 `P0,V0`），在以上域保持正值且对速度单调。已有测试使用 \(F(p,s)=0.3+0.5\,2.17^{-p/100}-10s\)，只在 \(0\le s<(0.3+0.5\,2.17^{-p/100})/10\) 时是正摩擦模型；该例的速度弱化表达式不能外推为任意速度下的物理定律。

- Numerical check method (which public API of the merged commit, what to compute, what to compare against, suggested tolerance and why; if finite differences are needed, the step-size strategy):

  **公开入口。** 后续在合并版本的 ContactApp 上通过 `.i` 输入配置 `[Constraints] type=ComputeDynamicFrictionalForceLMMechanicalContact`、`function_friction`、`friction_lm`、三维 `friction_lm_dir`，搭配公开的 `NormalMortarMechanicalContact` / `TangentialMortarMechanicalContact`、dual LM 变量、`ParsedFunction`、NewmarkBeta 时间积分。以现有 `frictional-mortar-3d-dynamics-light-function.i` 为结构参考，确保实际使用的 UserObject 的 boundary/subdomain 名称与所造网格匹配。通过 `NodalValueSampler` 导出法向与两个切向 LM，通过 `MortarFrictionalPressureVectorAux` 输出笛卡尔切向乘子向量，使用 Exodus/CSV 保存每个接受步。类内 `oldNodalSlipRate()` 和 `computeFrictionValue()` 都是 protected，不能把直接调用它们作为公开 API 检验方法。

  **I1/I2 的黑箱可观察量。** 构造有足够重叠的平面接触小块，法向加载建立正压力，法向位移保留求解自由度；用边界条件规定接触面的均匀切向运动，保证持续滑移，避免用一个顶部速度代替未测量的真实界面速度。由积分器实际速度或其 Newmark 递推重建每个接受步的界面速度；对均匀场，投影的精确参考即该速度，无需读内部 map。使用上面的正单调诊断函数，计算
  \[
  \mu_{obs,i}^n=\|\boldsymbol\lambda_i^n\|/p_i^n,
  \qquad \mu_{ref,i}^n=0.2+0.02p_i^{n-1}/P_0
       +0.5\sqrt{\|\mathbf U^{n-1}\|^2+\delta^2}/V_0.
  \]
  只对当前滑移且间隙误差足够小的节点比较两者。若要检查非均匀场，独立用同一接受步的位移、实际速度和 displaced mortar 几何重建 dual 测试函数积分 \(W,D\)，再将其比值代入参考函数；不能把另一个未核对的“real tangential velocity”输出当作真值。还可从单调函数反解 \(s_{reg,obs}=V_0(\mu_{obs}-0.2-0.02p^{n-1}/P_0)/0.5\)，直接比较参考范数。

  改变网格面积、相邻单元数、精确度足够的求积阶数、`normalize_c` 和 1/2/4 个 MPI 进程，比较相同物理场下的 \(\mu_{obs}\)；压力相关项须使用每个算例实际的旧压力扣除。改变两相邻接受步的速度幅值及压力，保证新旧参考系数有明显差异，验证采用旧值；先从静止初始化，再做非零且变化的滑移。新进入界面的自由度单独测试“无历史→零速度范数正则化”的约定，并区分它与本来已有历史的节点。利用合并源码支持的 `--test-restep` 比较相同接受时间点的历史/LM 输出，检查重试不污染历史；不使用恢复功能，源码明确不支持 recover。

  **I3 及方向检查。** 使用上述输出核验圆盘界、滑动等式、\(-\lambda\cdot W\le0\)；在 3D 使用两个切向分量的欧氏范数，不能分别要求每一分量饱和。反转速度、在固定平面旋转速度方向以及旋转完整坐标设置，检查系数保持不变、笛卡尔乘子按旋转变换。另以 `function_friction` 常量函数与相同正 `mu` 做配对计算。公开 C++ 模板 `Moose::Mortar::Contact::frictionalContactResidual<Real,k>` 可在将来小型测试中输入已知 \(\lambda,W,c_t,\Delta t,p,0,\mu,\epsilon\)，对 `ONE` 和 `TWO` 两种残差分别验证同一库仑根；该纯代数测试不能替代实际对象的归一化和历史检查。以 `mu=0`、负值测试输入拒绝，以小正值测试接受。

  **容差。** 恒定场参考积分的代数检查建议相对 \(10^{-11}\)、绝对 \(10^{-12}V_0\)，用于双精度小网格的求和误差；求解后的 \(\mu\) 比较先取 \(10^{-6}\) 相对加 \(10^{-8}\) 绝对，牵引误差按 \(P_0\) 缩放，功率按 \(D_iP_0V_0\) 缩放到 \(10^{-6}\)。要求 \(|c_{eff}\bar g|/p\ll10^{-6}\) 且 \(p\) 远离门限，否则将已量化的间隙/求解误差计入容差，不能把它误判为速度输入错误。收紧非线性容差一个数量级后核对误差是否稳定；不同进程结果允许浮点归约误差。无需有限差分；不从相邻位移输出的一阶差分推断 Newmark 速度，避免把时间离散误差混入速度归一化检验。跨面积、网格的比较只要求规定的恒定物理场与同一函数输入一致，不要求整个弹性动力学响应完全一样。

- Do the PR's own tests already check this invariant (yes/partial/no, naming test files and assertions):

  **I1：partial；I2：partial；I3：partial。** [mortar_dynamics/tests](https://github.com/idaholab/moose/blob/f90f508580f7611f7d75ee458cd177affc7288d1/modules/contact/test/tests/mortar_dynamics/tests) 中 `frictional-mortar-3d-dynamics-light-function` 用 `Exodiff` 比较此次更新的 `.e` 金标准，`rel_err=5e-4`、`abs_zero=1e-5`、要求 `ad_size>=64`。对应 [输入文件](https://github.com/idaholab/moose/blob/f90f508580f7611f7d75ee458cd177affc7288d1/modules/contact/test/tests/mortar_dynamics/frictional-mortar-3d-dynamics-light-function.i) 采用压力和速度相关函数，`dt=0.02`、`end_time=0.04`，间接覆盖三维瞬态结果，但没有显式断言 \(W/D\)、解析摩擦系数、圆盘界或耗散。

  本次 diff 只修改该二进制金标准，没有新增恒定场、面积/求积分区不变性、历史切换、重试或非正 `mu` 的测试断言；这些针对性性质在 PR 新增测试中为 **no**。现有同目录常数摩擦、edge-dropping 及其它动力学金标准提供背景回归覆盖，不能视为这些物理关系的独立证明。PR 描述中的摩擦系数约由 -0.48 变为 +0.48 是作者报告，本次未读取或数值复验二进制金标准。

- Do you expect the change to satisfy it (satisfies / may not satisfy / uncertain), with one or two sentences of reasoning (speculation, not a conclusion):

  **I1：satisfies（推测）**，分子分母使用同一测试函数和几何权重，且分母现在总在 owner 上聚合，直接支持常量再现。**I2：satisfies（推测）**，源码明确使用旧压力、归一化旧速度、正则化范数和重试保护；新界面节点及并行生命周期仍需实测。**I3：satisfies（在给定正系数、零间隙和收敛条件下的推测）**，共享残差的根满足库仑投影关系；任意函数返回负值的更强断言则 **may not satisfy**，本次改动没有提供这种保护。

- Environment: merge commit SHA, language and dependencies needed, estimated setup difficulty (low/medium/high):

  合并提交 `f90f508580f7611f7d75ee458cd177affc7288d1`。实现为 C++，公开运行入口是 MOOSE ContactApp 的输入文件；完整检查需要该版本 MOOSE framework/contact 及所用 solid mechanics 模块、匹配的 libMesh/PETSc、MetaPhysicL AD、TIMPI/MPI、C++ 工具链和 Make，输出检查需要 Exodus/NetCDF 读取器或 CSV 工具，官方测试框架另需 Python 与 Exodiff。构建关系见 [contact/Makefile](https://github.com/idaholab/moose/blob/f90f508580f7611f7d75ee458cd177affc7288d1/modules/contact/Makefile)；从零准备估计 **high**，已有匹配版本环境则较低。本次只读 `input` 与 GitHub 原始源码，未克隆、安装、编译或执行项目代码，也未访问另一审计轨道。

- Formalization (Lean) feasibility: can the invariant be written as a proposition over finite-dimensional reals, and what are the main obstacles:

  **可以。** 固定有限求积点、节点和进程分区后，I1 是实数有限和、非零除法及线性映射上的命题，常量再现和分区聚合证明很直接；I2 的正交范数不变性也属有限维实内积空间，历史约定可用有限状态转移关系表达。I3 可在 \(D,p,c_t,\Delta t>0\)、\(\mu\ge0\) 和零间隙假设下，将圆盘投影或两种残差为零分情况化简，推出滑动方向、大小和非正功率。

  主要障碍是把证明假设连接到真实实现：dual mortar 几何及 \(D_i>0\)、法向/切向基一致性、活动集与 edge dropping、owner 通信不重不漏、接受/拒绝步状态推进，以及浮点误差和非线性求解停止条件。无需先形式化整个 PDE，但仅证明有限和等式不能证明这些 C++ 生命周期条件；也不能从 lagged 摩擦耗散推出整个 Newmark 系统的总能量单调性。
