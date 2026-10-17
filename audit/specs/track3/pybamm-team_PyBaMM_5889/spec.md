## pybamm-team/PyBaMM#5889 fix: look up Chayambuka2022 NVPF diffusivity at the solid concentration

- What changed (1–3 sentences, naming files and formulas):
  `packages/pybamm/src/pybamm/input/parameters/sodium_ion/Chayambuka2022.py` 的 `NVPF_diffusivity_Chayambuka2022(sto, T)` 将查表坐标从 `sto * Initial concentration in electrolyte` 改为 `sto * Maximum concentration in positive electrode`，即由 θ cₑ,₀ 改为固相浓度 cₛ = θ C。默认 C = 15320 mol m⁻³，而 cₑ,₀ = 1000 mol m⁻³；`test_Chayambuka2022.py` 同时增加正极扩散系数回归点、改用相对容差并修正其他期望值，另更新 CHANGELOG 和 sodium-ion 示例 notebook 输出。电解液扩散系数期望值的修正不对应新的电解液物理实现。

- Invariants (mathematical statements; write them as checkable equalities/inequalities, give variables, sign conventions, domain of validity):
  令 θ 为无量纲正极粒子钠占位率，C > 0 为正极最大固相浓度，cₛ = θ C，e = cₑ,₀ 为初始电解液浓度，T₀ = 298.15 K。令合并提交 `data/D_p.csv` 的有序数据为 (xᵢ, dᵢ)，x 的单位为 mol m⁻³，d 的单位为 m² s⁻¹；定义 L(c) 为这张表的分段线性插值。用 D(θ,T;C,e) 表示经过参数处理的公开正极扩散系数函数输出。

  1. **Primary — 固相浓度坐标与表格一致。** 对 xᵢ ≤ θ C ≤ xᵢ₊₁，令 λ = (θ C − xᵢ)/(xᵢ₊₁ − xᵢ)，应有
     D(θ,T₀;C,e) = L(θ C) = (1−λ)dᵢ + λdᵢ₊₁。
     特别地 D(xᵢ/C,T₀;C,e) = dᵢ。归一化变换还应满足 D(c/C₁,T₀;C₁,e) = D(c/C₂,T₀;C₂,e)，前提是保持同一张浓度表且两个占位率均在 [0,1] 内。这是对参数坐标语义的检查，并非宣称不同材料具有同一扩散曲线。
  2. **Secondary — 固相材料扩散系数不依赖无关的电解液初始条件。** 固定 θ、C、T₀ 和数据表，对任意测试域内 e₁,e₂，应有 D(θ,T₀;C,e₁) = D(θ,T₀;C,e₂)。该条件仅针对这个局部本构函数；改变 e 可以通过电解液传输和反应动力学改变整个电池的轨迹，故不要求电压或容量保持不变。
  3. **Secondary — 表内正扩散与局部包络。** 对上述每个区间，0 < min(dᵢ,dᵢ₊₁) ≤ D ≤ max(dᵢ,dᵢ₊₁)。整张表内进一步有 1.95×10⁻¹⁷ ≤ D ≤ 2.51×10⁻¹⁵ m² s⁻¹。采用朝浓度增加方向为正的浓度梯度 g = ∇cₛ，以及 Fick 通量 J = −Dg，则 J·g = −D‖g‖² ≤ 0，g ≠ 0 时严格小于零。这个局部耗散符号不等式不等同于对整套 DFN 离散系统的稳定性或能量守恒证明。数据本身非单调，不应要求 D 随 θ 全局单调。

- Sources of the invariants (first principles / literature formula / consistency with other implementations in the library / dimensions):
  主要证据是合并提交的 [Chayambuka2022 参数源码](https://github.com/pybamm-team/PyBaMM/blob/eda03ed6d1003b186cc2d70ad8b3240a7b6350f1/packages/pybamm/src/pybamm/input/parameters/sodium_ion/Chayambuka2022.py) 和 [D_p.csv 原始表及单位表头](https://github.com/pybamm-team/PyBaMM/blob/eda03ed6d1003b186cc2d70ad8b3240a7b6350f1/packages/pybamm/src/pybamm/input/parameters/sodium_ion/data/D_p.csv)。占位率定义 θ = cₛ/C 给出第一条；同文件的负极 `HC_diffusivity_Chayambuka2022` 也以 θ 乘该电极最大浓度。仅凭量纲不能识别原错误，因为 e 与 C 的量纲相同，必须区分固相与电解液的物理身份。
  [Interpolant 源码](https://github.com/pybamm-team/PyBaMM/blob/eda03ed6d1003b186cc2d70ad8b3240a7b6350f1/packages/pybamm/src/pybamm/expression_tree/interpolant.py) 确认默认一维插值为 linear；因此凸组合公式和包络属于这个库的本构表示约定。第二条来自该参数模型将固相扩散表示为固相浓度函数的假设；第三条的通量符号来自 Fick 定律和正扩散要求。源码引用 K. Chayambuka et al., “Physics-based modeling of sodium-ion batteries part II. Model and validation”, Electrochimica Acta 404 (2022), 139764，DOI 10.1016/j.electacta.2021.139764，并注明参数来自 COMSOL 示例；本审计未独立读取论文或 COMSOL，故不将表格与原始实验的准确性视为已验证，也不推断额外温度依赖。

- Common operating domain (range of each input, with the basis):
  核心物理域固定 C = 15320 mol m⁻³，T = 298.15 K，e = 1000 mol m⁻³，依据该参数集默认值。只在已提供数据的 cₛ ∈ [131.578947, 15197.36842] mol m⁻³ 上主张表格精度，即 θ ∈ [131.578947/15320, 15197.36842/15320]，约 [0.0085887, 0.991995]。不把整个 [0,1] 的端部外推视为已验证物理域：该 Interpolant 默认允许外推，但缺少表外数据支撑。
  对第二条，用 e ∈ [500,1500] mol m⁻³ 作围绕默认值的有限扰动域；这是用于识别错误耦合的测试选择，不是文献验证的电池运行浓度范围。对第一条的坐标变换，用 C₁,C₂ ∈ [15320,18384] mol m⁻³（默认值的 1–1.2 倍），保持 cₛ 在原表内，再取 θⱼ = cₛ/Cⱼ；这一域保证 θⱼ ∈ (0,1)，仅用于代数一致性测试。若做通量符号检查，可选 g = (10⁶ mol m⁻⁴)u，u ∈ [−1,1]³，包含零与混合符号方向；幅度只是避免浮点极端值的测试尺度，不是额外的材料适用范围。T 参数虽被函数接受，但当前函数未使用它；不将真实材料的温度不变性列为物理不变量。

- Numerical check method (which public API of the merged commit, what to compute, what to compare against, suggested tolerance and why; if finite differences are needed, the step-size strategy):
  下面是后续执行方案，本次没有安装、导入或运行 PyBaMM，也没有数值验证结果。
  使用 `p = pybamm.ParameterValues("Chayambuka2022")`，按需要用 `p.update({...})` 设置完整名称 `"Maximum concentration in positive electrode [mol.m-3]"` 与 `"Initial concentration in electrolyte [mol.m-3]"`。构造 `expr = p["Positive particle diffusivity [m2.s-1]"](pybamm.Scalar(theta), pybamm.Scalar(298.15))`，通过 `p.evaluate(expr)` 取得数值（若返回单元素数组，提取该元素）。每组 C、e 使用新的 ParameterValues 和表达式，以免缓存或前次处理污染检查。
  独立读取固定合并提交的 CSV 两列，不从被测表达式的 children 或内部 x/y 反推出期望值，也不再调用 `pybamm.Interpolant` 作为 oracle。先确认 x 严格递增且 d > 0，再用上面的两点凸组合直接计算 L。
  第一条检查全部 34 个节点及 33 个区间中点，并对每个区间采样若干 λ ∈ (0,1)，令 θ = ((1−λ)xᵢ+λxᵢ₊₁)/C。覆盖端点和局部极值附近，并进行上述不同 C、相同 cₛ 的成对比较。建议对未舍入 CSV oracle 使用 `abs(D−L)/abs(L) ≤ 1e−10`，绝对容差为零；这是同一十进制表的双精度算术一致性检查，容差给浮点坐标换算和区间定位留余量，不代表实验精度。
  第二条对相同 θ、C 取 e = 500、1000、1500 及区间内随机值，两两比较相对误差 ≤ 1e−12，绝对容差为零。原实现会依赖 e，故该检验能区别浓度身份错误；正确表达式不引用 e，预期差异为零。
  第三条复用所有采样，先要求有限且 D > 0，再检查局部上下界，边界允许 `1e−10 * max(dᵢ,dᵢ₊₁)` 的浮点余量。用 API 返回的 D 和选定 g 构造 J，检查 J·g ≤ 0；该步骤只核查本构后果，不声称执行了库内完整通量离散。D 的严格正值检查不能被大绝对容差替代。
  额外的已公布回归点可取 θ = 0.3、0.5、0.9，对应约 5.6673e−17、8.3317e−17、3.5627e−17 m² s⁻¹；这些舍入值采用 1e−4 相对容差即可。无需有限差分：直接代入查表公式和参数成对测试足够，且分段线性插值在节点一般不可微。PR 描述的 1C/2C 截止容量仅是多因素模拟结果，不作为这里的硬性 oracle。

- Do the PR's own tests already check this invariant (yes/partial/no, naming test files and assertions):
  **第一条 partial。** [packages/pybamm/tests/unit/test_parameters/test_parameter_sets/test_Chayambuka2022.py](https://github.com/pybamm-team/PyBaMM/blob/eda03ed6d1003b186cc2d70ad8b3240a7b6350f1/packages/pybamm/tests/unit/test_parameters/test_parameter_sets/test_Chayambuka2022.py) 的 `test_positive_particle_diffusivity` 在 θ = 0.3、0.9 断言上述值，`pytest.approx(..., rel=1e-4, abs=0)`；`test_functions` 在 θ = 0.5 断言 8.3317e−17，`rel=1e-3, abs=0`。这些能检测默认参数下的原错误，但没有遍历节点、独立插值 oracle 或改变 C 的坐标一致性检查。
  **第二条 no。** PR 修改的测试均固定默认电解液初始浓度，没有参数独立性断言。
  **第三条 partial。** 三个正值期望隐含检验这三个点的正性；没有覆盖整表局部包络或 Fick 通量符号的断言。notebook 输出更新不是这些不等式的自动化断言。上述结论限于本 PR 自身改动和已读取测试文件，不声称扫描了整个仓库的测试。

- Do you expect the change to satisfy it (satisfies / may not satisfy / uncertain), with one or two sentences of reasoning (speculation, not a conclusion):
  **satisfies（推测，限上述域及这三个不变量）。** 合并实现明确在 θ C 上进行线性插值，不引用 e，且表内所有 dᵢ 均为正，因此三个条件从表达式结构上应成立；这不是数值运行结论，也不能证明经验曲线的实验准确性或表外、非等温行为。

- Environment: merge commit SHA, language and dependencies needed, estimated setup difficulty (low/medium/high):
  合并 SHA：`eda03ed6d1003b186cc2d70ad8b3240a7b6350f1`。Python；[该提交的 packages/pybamm/pyproject.toml](https://github.com/pybamm-team/PyBaMM/blob/eda03ed6d1003b186cc2d70ad8b3240a7b6350f1/packages/pybamm/pyproject.toml) 要求 Python >=3.10,<3.16，NumPy >=2.0.0，SciPy >=1.11.4 且 !=1.18.0，pybammsolvers >=0.11.0,<0.12.0，以及列出的其余运行依赖；读取的 Interpolant 模块也直接 import CasADi，后续环境需包含兼容版本。可使用 Python 3.12 与符合该提交约束的完整依赖环境；这里的参数求值不需要网格、DFN 求解、GPU 或 notebook。预计设置难度 **medium**：检查本身轻量，但固定历史合并版本及其 solver/CasADi 二进制依赖需匹配；本次未验证安装可用性。

- Formalization (Lean) feasibility: can the invariant be written as a proposition over finite-dimensional reals, and what are the main obstacles:
  可以。以 34 个有理数节点和正系数嵌入 ℝ，定义分段仿射 L；假设 xᵢ < xᵢ₊₁、C > 0、λ ∈ [0,1]，第一条及浓度重参数化化为乘除代数恒等式，第二条化为 e 不出现在定义中，第三条化为正数凸组合的界及有限维内积范数平方非负。无需形式化 PDE 或时间积分即可表达这三个命题。
  主要障碍是证明实际 Python 参数替换与 SciPy 插值执行精化该实数定义，以及处理 IEEE-754 舍入和节点分支；实数层证明本身不会自动证明机器精度误差界。数据表的实验来源和 θ 的材料物理解释仍是外部建模前提，无法由这些代数证明推出。
