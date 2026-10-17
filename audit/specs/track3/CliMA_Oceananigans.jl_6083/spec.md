## CliMA/Oceananigans.jl#6083 Include drift velocities in the implicit part of adaptive implicit tracer advection

- What changed (1–3 sentences, naming files and formulas):
  `hydrostatic_free_surface_{ab2,rk}_step.jl` 与 `hydrostatic_free_surface_tendency_kernel_functions.jl` 现在共同组装示踪物总速度 $\boldsymbol U_c=\boldsymbol U_{transport}+\boldsymbol U_{bgc,c}+\boldsymbol U_{closure,c}+\boldsymbol U_{forcing,c}$，使自适应显式通量和隐式求解使用同一个速度。`advective_skew_diffusion.jl` 改为在各流函数交错位置用完整二维坡度计算 GM 截断，并应用 `minimum_bz`；`isopycnal_skew_symmetric_diffusivity.jl` 将截断参数转换到网格浮点类型。`derivative_operators.jl` 和 `mutable_immersed_grid.jl` 补充两个交错位置的坐标坡度算子。

- Invariants (mathematical statements; write them as checkable equalities/inequalities, give variables, sign conventions, domain of validity):

  1. **Primary — 总漂移速度不能在显隐分裂中丢失，且结果不依赖同一速度的来源分解。** 取向上为正的 $z$，浓度 $c$，无反应和物理扩散。在每个垂直面及实际子步长度 $\Delta\tau>0$ 上，设 $w=w_{transport}+w_{bgc}+w_{closure}+w_{forcing}$，$\alpha=|w|\Delta\tau/\Delta z$，$s=1$（$\alpha=0$）或 $s=\min(1,C_*/\alpha)$。必须满足
     \[
     w^e=sw,\qquad w^i=(1-s)w,\qquad w^e+w^i=w,\qquad |w^e|\Delta\tau/\Delta z\le C_*.
     \]
     更直接的公开 API 不变量是：在同一网格、初值、时间步、边界和各阶段总速度相同的两个模型中，$c_A^n=c_B^n$，即把固定漂移从生物地球化学接口移到 `AdvectiveForcing` 或分摊到两者，不应改变浓度数组。特别是两个相反速度抵消时，总速度为零，$c^{n+1}=c^n$。连续方程 $\partial_t c+\partial_z(wc)=0$ 对常数 $w$、远离边界的脉冲给出 $\bar z(t)-\bar z(0)=wt$，$\bar z=\int zc\,dz/\int c\,dz$；这是数值解应趋近的物理基准，**不是有限网格上要求逐位成立的恒等式**。不同 $C_*$ 的浓度形状可因隐式迎风扩散而不同。

  2. **Secondary — 封闭域示踪物总量守恒。** 对静止网格、固定体积 $V_i>0$、无反应/源汇/扩散边界通量且总法向速度在边界为零，
     \[
     M^n=\sum_i V_i c_i^n,\qquad M^{n+1}=M^n.
     \]
     周期面通量成对抵消也适用；一维非均匀柱用 $V_k=A\Delta z_k$，不能直接用未加权 `sum(c)`。此条件不要求漂移速度无散。它能检查新增隐式漂移是否错误地增减质量，但单独通过不能证明输运速度正确；不把 WENO 全离散正性或示踪物方差守恒强加为本 PR 的保证。

  3. **Secondary — GM 坡度截断具有完整水平范数依赖，并产生无散涡旋速度。** 限定静止笛卡尔网格、单个 GM closure、`κ_symmetric=nothing`。取浮力 $b$ 单位 m s⁻²，稳定层结 $b_z>0$，$m=\texttt{minimum_bz}\ge0$，$q=\max(b_z,m)$。在 $q>0$ 时，$\boldsymbol S=-(b_x,b_y)/q$，$S_m>0$，$\epsilon=\min(1,S_m^2/|\boldsymbol S|^2)$，零坡度处定义 $\epsilon=1$，$\boldsymbol P=\kappa\epsilon\boldsymbol S$。应满足
     \[
     0\le\epsilon\le1,\quad |\epsilon\boldsymbol S|\le S_m,\quad
     \boldsymbol P(R\boldsymbol S)=R\boldsymbol P(\boldsymbol S),\quad
     (u_e,v_e,w_e)=(-\partial_zP_x,-\partial_zP_y,\partial_xP_x+\partial_yP_y),\quad \nabla\cdot\boldsymbol U_e=0.
     \]
     $R$ 是水平正交旋转；在离散网格上仅要求对应的解析制造解/网格对称变换满足旋转关系，不要求任意旋转的采样场逐位相同。$\kappa$ 单位 m² s⁻¹，$P$ 为 m² s⁻¹，涡旋速度为 m s⁻¹。对 $m=0,b_z\le0$，预期 $\boldsymbol P=0$、所有速度有限；正 $m$ 会人为抬高分母，因此不能再要求不稳定层结自动关闭。笛卡尔交错差分应满足 $D_xu_e+D_yv_e+D_zw_e=0$ 至舍入误差。这里不把该离散恒等式直接推广到活动地形网格；该部分需要单独推导体积/面积度量与边界条件。

- Sources of the invariants (first principles / literature formula / consistency with other implementations in the library / dimensions):
  不变量 1 来自保守输运方程及速度可加性；$\bar z$ 关系来自对输运方程乘 $z$ 后积分，边界项须消失。分裂参数的确切定义见合并版本 [adaptive_implicit_vertical_advection.jl](https://github.com/CliMA/Oceananigans.jl/blob/1734ff9e4509198f945c5e7740a0bccacc69a720/src/Advection/adaptive_implicit_vertical_advection.jl)。不变量 2 来自有限体积共享面通量的望远镜求和，和物理粒子数/示踪物质量守恒。
  不变量 3 中坡度来自等浮力面 $db=0$，无散关系来自混合偏导交换；截断公式是库的模型约定，并非普适物理定律。[ISSD 构造器与截断公式]$https://github.com/CliMA/Oceananigans.jl/blob/1734ff9e4509198f945c5e7740a0bccacc69a720/src/TurbulenceClosures/turbulence_closure_implementations/isopycnal_skew_symmetric_diffusivity.jl) 引用 Gent–McWilliams (1990)、Redi (1982) 及 Gerdes, Koberle & Willebrand (1991)；此规格依据已读取的库公式，不声称已独立核查论文全文。涡旋速度的符号与交错位置见 [advective_skew_diffusion.jl](https://github.com/CliMA/Oceananigans.jl/blob/1734ff9e4509198f945c5e7740a0bccacc69a720/src/TurbulenceClosures/turbulence_closure_implementations/advective_skew_diffusion.jl)。[旋转张量源文件](https://github.com/CliMA/Oceananigans.jl/blob/1734ff9e4509198f945c5e7740a0bccacc69a720/src/TurbulenceClosures/isopycnal_rotation_tensor_components.jl) 的文档字符串把稳定层结写成 (b_z<0$，与实现的正分母分支不一致；本规格明确采用向上坐标下稳定浮力梯度 $b_z>0$，不照抄该文档符号。

- Common operating domain (range of each input, with the basis):
  **输运：** CPU Float64 为基线，另用 Float32 和可用 GPU 复核。均匀柱 $z\in[-H,0]$，$H=400$–$1600$ m，$N_z=40,80,160,320$，浓度峰值 0.1–10（任意一致浓度单位），平滑正脉冲宽度至少 2–4 个初始网格单元；加密时保持物理宽度不变。速度绝对值 10–200 m/day，正负两个方向及精确抵消；总 Courant 数取 0、0.1、0.49、0.51、1、4、8，$C_*=0.3,0.5$，时间步由 $\Delta t=C\Delta z/|w|$ 决定，零速单独选正时间步。以 3–10 步和固定物理终止时间加密两种方式测试。上述是围绕 PR 的 $H=400,N_z=40,w=-100$ m/day、$\Delta t=0.4$ day、三步算例设计的审计范围，不是声称所有组合均达到给定精度；质心检查只采用脉冲及数值尾部尚未明显触及边界的组合，必要时扩域验证。
  **GM：** 三维均匀有界 `RectilinearGrid`，每方向 16–64 单元，水平长度 10–100 km、深度 0.4–2 km；$b_z=B=10^{-6}$–$10^{-4}$ s⁻²，$m\in\{0,B/2,2B\}$，$S_m=10^{-3}$–$10^{-2}$，
  $\boldsymbol S=rS_m$\cos\theta,\sin\theta),\ r\in\{0,0.5,0.99,1,1.01,\sqrt2,2\},\ \theta\in[0,2\pi)$。这些是小坡度模型及默认 $S_m=10^{-2}$ 附近的分支覆盖范围。$\kappa\ge0$ 取 0–2000 m²/s，包括随深度线性变化的非负系数。另测 $m=0,b_z=0$ 和负值作为关闭分支，不将其解释为该闭合模型可准确模拟对流。所有输入有限，避免下溢、溢出；水平输运仍显式，演化测试须满足其 CFL。

- Numerical check method (which public API of the merged commit, what to compute, what to compare against, suggested tolerance and why; if finite differences are needed, the step-size strategy):
  **1. 公开构造与总速度检查。** 用 `RectilinearGrid`、`ZFaceField`、`set!`、`HydrostaticFreeSurfaceModel`，设置 `momentum_advection=nothing, free_surface=nothing`，并分别用 `timestepper=:QuasiAdamsBashforth2` 和 `:SplitRungeKutta3`。示踪物方案为 `WENO$order=5, time_discretization=AdaptiveVerticallyImplicitDiscretization(cfl=0.5))`。通过 `Oceananigans.Biogeochemistry.AbstractBiogeochemistry` 的 `required_biogeochemical_tracers` 和 `biogeochemical_drift_velocity(::..., ::Val{:D})` 扩展实现无反应漂移，照 PR 的 `SinkingParticles`；GPU 版本还需 `Adapt.adapt_structure`。`w` 在两个端面为零，内部为给定常数。另一模型通过公开 `forcing=(D=AdvectiveForcing(w=w),)` 施加完全相同的面速度；第三个把速度按 $\lambda$ 与 $1-\lambda$ 分配给两接口，取 $\lambda=0.25,0.5,1.5$，最后额外测两接口精确相反的速度。应给各模型独立的可变 advection 对象和初态，不能共享时步缓存。
  通过 `time_step!` 演化，`Array$interior(model.tracers.D))` 读取结果，比较 $\|c_A-c_B\|_\infty/c_{scale}\le10^{-11}$（Float64，十步内），Float32 初始阈值 $2\times10^{-5}$。精确抵消同时与初值比较。读取 `znodes$grid, Center())` 计算质心：PR 基线应下降约 120 m，先用 5% 误差限；扩域和同时减半 $\Delta z,\Delta t$（固定 Courant 和物理终止时间）后误差应下降，再以小时间步显式结果作独立收敛参考。不能要求高 Courant 的完整浓度形状等于显式解，也不能把边界滞留误判为速度损失。这个可直接构建的比较覆盖生物漂移与 advective forcing；closure 通道还应以单个 GM closure 生成的速度做冻结场重放比较，确保所有阶段速度相同才要求数组相等，不把随浮力演化的两个不同流场作严格等价比较。
  **2. 质量。** 在上述每一步计算 $M$，包括低/高 Courant、两个时间积分器、变号空间速度以及非均匀柱。均匀柱可复用 `sum$D)`，非均匀柱从定义网格时给定的面坐标计算 $\Delta z_k$。要求 $|M^n-M^0|/\sum V_i|c_i^0|\le10^{-12}$（三步 Float64 基线），十步扩展可用 $10^{-11}$；Float32 用 $2\times10^{-5}$。这些是约百至千次舍入量级的起始容差，需记录归约方式、规模与精度；不以百分比容差掩盖守恒错误。
  **3. GM 制造解与无散检查。** 用 `HydrostaticFreeSurfaceModel$grid; buoyancy=BuoyancyTracer(), tracers=(:b,:c), closure=IsopycnalSkewSymmetricDiffusivity(κ_skew=κfun, κ_symmetric=nothing, skew_flux_formulation=AdvectiveFormulation(), isopycnal_tensor=SmallSlopeIsopycnalTensor(minimum_bz=m), slope_limiter=FluxTapering(S_m)))`（closure 类型可从 `Oceananigans.TurbulenceClosures` 显式导入）。令 $b=G_xx+G_yy+Bz$，$\kappa$z)=\kappa_0+a(z+H/2)\ge0$，其中 $\kappa_0>|a|H/2$。`set!(model, b=...)` 会更新 closure fields，随后用 `interior(model.closure_fields.u/v/w)` 读取诊断速度，无需调用私有 slope/kernel 函数；这一更新路径由合并版本 [set! 实现](https://github.com/CliMA/Oceananigans.jl/blob/1734ff9e4509198f945c5e7740a0bccacc69a720/src/Models/HydrostaticFreeSurfaceModels/set_hydrostatic_free_surface_model.jl) 确认。
  在距所有边界至少三个单元的内部，独立计算 $q=\max$B,m),S_x=-G_x/q,S_y=-G_y/q$，预期 $u_e=-a\epsilon S_x,v_e=-a\epsilon S_y,w_e=0$。必须取非零 $a$：常系数、常坡度在内部给零速度，无法分辨错误截断。尤其 $S_x=S_y=S_m$ 时 $\epsilon=1/2$，能检出旧的逐分量截断；交换 x/y 以及旋转梯度复核范数依赖。零坡度、$\kappa=0$、$m=0,B\le0$ 分支须输出有限的零速度。将线性浮力改为平滑空间变化场后，用实际交错速度数组作相邻面差分，检查 $D_xu+D_yv+D_zw$，至少覆盖内部共享完整 stencil 的单元。采用误差限 $10^{-10}U_{scale}$（Float64）或 $10^{-4}U_{scale}$（Float32）检查制造速度；散度容差取相同相对系数乘 $U_{scale}$1/\Delta x+1/\Delta y+1/\Delta z)$。$U_{scale}$ 取该批非零制造解的速度尺度，零解也用该尺度作为绝对容差。Float32 网格保留默认 Float64 closure 参数再测一次，检查所有输出有限以及相对于 Float64 基准的精度；CUDA 是可选兼容性复核。
  上述线性制造解无需数值微分连续函数；无散检查使用有限体积定义本身的网格差分。若扩展到非线性制造解，用解析导数作基准，网格尺寸减半检查截断误差；若用 $$c(\delta t)-c(0))/\delta t$ 诊断倾向，则连续减半 $\delta t$，寻找误差下降后舍入主导前的平台，不使用单一极小步长。

- Do the PR's own tests already check this invariant (yes/partial/no, naming test files and assertions):
  **1: partial。** [test/unit/adaptive_implicit_vertical_advection.jl](https://github.com/CliMA/Oceananigans.jl/blob/1734ff9e4509198f945c5e7740a0bccacc69a720/test/unit/adaptive_implicit_vertical_advection.jl) 新增 `AIVA carries the biogeochemical drift velocity through the implicit solve`，对 `archs` 和两个时间积分器断言 `descent ≈ speed * Δt * nsteps rtol=0.05`；只覆盖一个向下、生物漂移、Courant=4 算例，未新增速度来源分解/抵消、advective forcing、closure 通道或收敛断言。同文件既有 `AIVA flux is CFL-limited when α > cfl` 断言显式通量缩放，并有低 CFL 显式一致性测试，但不能替代总速度端到端检查。
  **2: yes，限新增均匀柱算例；整体 partial。** 新测试明确断言 `sum(D) ≈ sum(D₀) rtol=1e-12`；未新增非均匀体积或更广速度来源覆盖。
  **3: no（本 PR 新增测试范围内）。** 输入 diff 没有新增 GM 截断、二维梯度旋转、散度、Float32 GM 或活动网格物理断言。未穷尽整库所有既有测试，不将此表述为整库完全没有相关测试。

- Do you expect the change to satisfy it (satisfies / may not satisfy / uncertain), with one or two sentences of reasoning (speculation, not a conclusion):
  **1: satisfies（推测）**，显式与两个隐式时间步入口现在共享总速度组装；仍须实测阶段时间步、多个来源及动态 closure 的一致性。**2: satisfies（推测）**，共享面通量的守恒结构未被替换，新增质量断言针对正确量，但没有运行证据。**3: satisfies（推测，限上述静止笛卡尔域）**，完整二维截断与流函数差分形式支持这些等式；活动地形网格的度量相容性及 GPU 混合类型执行仍为 uncertain，不能从这份源码阅读判定通过。

- Environment: merge commit SHA, language and dependencies needed, estimated setup difficulty (low/medium/high):
  合并 SHA `1734ff9e4509198f945c5e7740a0bccacc69a720`，Julia；[Project.toml](https://github.com/CliMA/Oceananigans.jl/blob/1734ff9e4509198f945c5e7740a0bccacc69a720/Project.toml) 标注 Oceananigans `0.113.4`、Julia compat `1.10`。以后复核应固定该 SHA，而不是仅装同名最新版本；核心依赖包括 Adapt 4.1.1 compat、KernelAbstractions 0.9.39 compat、OffsetArrays、StaticArrays、FFTW、MPI、SeawaterPolynomials，具体解析版本需保存。CPU 小网格设置难度 medium（Julia 包及二进制依赖），GPU 复核 high/取决于现成硬件及驱动，CUDA 是可选 weak dependency。此次没有克隆、安装、编译或执行任何模拟/测试；只读 `./input` 与 GitHub 合并提交的公开源码和文件树，没有读取另一轨结果、讨论评论或 issue。

- Formalization (Lean) feasibility: can the invariant be written as a proposition over finite-dimensional reals, and what are the main obstacles:
  **可以。** 速度分裂是 $\Delta z,\Delta\tau,C_*>0$ 下含绝对值、min/max 和零速分支的实数命题；质量守恒可表达为体积向量与显式通量差分、隐式矩阵作用之间的左不变量，证明边界面抵消及隐式矩阵可逆后可组合至 AB2/RK 子步。GM 截断界限与正交协变性是二维代数，交错网格散度为零是有限差分矩阵可交换所导致的有限维恒等式。
  主要障碍是将这些数学对象连接到 Julia 实际 dispatch、交错位置、halo/边界与阶段状态；要证明完整数组等价还须证明各来源在每阶段组装出相同总速度。浮点类型转换、除零中间表达式和 NaN 选择语义需另建 IEEE 浮点层；连续质心精确平移到离散误差界的连接需要收敛/边界尾部分析，不能仅凭有限维速度恒等式推出。活动网格还需独立的度量恒等式证明，本规格的笛卡尔证明不覆盖它。
