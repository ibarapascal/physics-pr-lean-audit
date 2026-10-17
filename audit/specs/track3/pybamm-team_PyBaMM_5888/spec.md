## pybamm-team/PyBaMM#5888 fix: convert Ai2020 electrolyte diffusivity from cm2/s to m2/s

- What changed (1–3 sentences, naming files and formulas):
  `packages/pybamm/src/pybamm/input/parameters/lithium_ion/Ai2020.py` 的 `electrolyte_diffusivity_Ai2020(c_e,T)` 将原来的 `10**(-4.43-54/(T-229-0.005*c_e)-0.00022*c_e)` 乘以 `1e-4`，把 cm²/s 转成 m²/s；同时修正扩散率和电导率的文档单位，电导率公式未改变。`packages/pybamm/tests/unit/test_parameters/test_parameter_sets/test_Ai2020.py` 新增扩散率、电导率数值及与 Xu2019 的比较，并收紧原有参数测试的相对容差；`CHANGELOG.md` 记录修复。

- Invariants (mathematical statements; write them as checkable equalities/inequalities, give variables, sign conventions, domain of validity):
  这里的 \(c\) 是以 mol/m³ 表示的电解液盐浓度数值，\(T\) 是以 K 表示的绝对温度数值；\(D_A\)、\(D_X\) 分别是 Ai2020、Xu2019 公共参数接口返回的扩散系数，单位 m²/s。令 \(C=c/1000\) 为以 mol/L 表示的浓度数值，\(q=T-229-0.005c\)。以下取 \(c\geq0,q>0\)，数值检验进一步限制在下述矩形区域；这是经验本构式的不变量，不是对所有电解液的普遍定律。

  1. **Primary — SI 单位及本构式的绝对标度。**
     \[
     D_A(c,T)=10^{-4}\exp\!\left[\ln(10)\left(-4.43-\frac{54}{T-229-5C}-0.22C\right)\right]\quad\mathrm{m^2/s}.
     \]
     等价地，\(\log_{10}[D_A/(1\,\mathrm{m^2/s})]=-8.43-54/q-0.00022c\)。物理单位恒等式是 \(1\,\mathrm{cm^2/s}=10^{-4}\,\mathrm{m^2/s}\)，浓度单位转换是 \(1000\,\mathrm{mol/m^3}=1\,\mathrm{mol/L}\)。特别地，\(D_A(1000,298.15)\approx3.2227\times10^{-10}\,\mathrm{m^2/s}\)；这个近似值是舍入的检查点，完整等式才是基准。

  2. **Secondary — 正扩散率及此拟合的温度、浓度响应。** 对上述区域，\(0<D_A<\infty\)，且
     \[
     \partial_T\ln[D_A/(1\,\mathrm{m^2/s})]=\ln(10)\frac{54}{q^2}>0,
     \qquad
     \partial_c\ln[D_A/(1\,\mathrm{m^2/s})]=-\ln(10)\left(\frac{0.27}{q^2}+0.00022\right)<0.
     \]
     因此固定浓度升温时扩散率严格增加，固定温度增加盐浓度时严格减小。正号与 Fick 扩散通量 \(\mathbf N_{\rm diff}=-D_A\nabla c\) 一致，\(\mathbf N_{\rm diff}\cdot\nabla c=-D_A\lVert\nabla c\rVert^2\leq0\)；此处仅指扩散通量分量，不对包含迁移、对流的总通量作此断言。单独检验这些符号不能识别本 PR 修正的常数倍错误，必须同时检验 primary。

  3. **Secondary — 相同输运拟合在参数集之间一致。** 对相同的 SI 输入，\(D_A(c,T)=D_X(c,T)\)。仅比较两者使用的同一个 Valøen–Reimers 扩散率拟合，不要求两个完整电池参数集或放电轨迹相同。

- Sources of the invariants (first principles / literature formula / consistency with other implementations in the library / dimensions):
  Primary 的单位因子来自长度单位平方；经验公式来源是合并源码引用的 Valøen and Reimers (2005), *Transport properties of LiPF6-based Li-ion battery electrolytes*, JES 152, A882–A891, Eq. 14，DOI `10.1149/1.1872737`。本次没有读取论文全文，原始实验标定范围和论文符号未独立核实；公式归属以 [Ai2020 合并源码](https://raw.githubusercontent.com/pybamm-team/PyBaMM/343ece2062dbb39789dbfa977e0f2daf96e9e370/packages/pybamm/src/pybamm/input/parameters/lithium_ion/Ai2020.py) 为依据。温度与浓度导数是上述公式的直接数学推论；正扩散率与 Fick 定律的降梯度方向一致。[Xu2019 合并源码](https://raw.githubusercontent.com/pybamm-team/PyBaMM/343ece2062dbb39789dbfa977e0f2daf96e9e370/packages/pybamm/src/pybamm/input/parameters/lithium_ion/Xu2019.py) 先把浓度除以 1000，再计算同一拟合并乘以 `1e-4`，还明确说明采用以 10 为底的指数。它提供库内一致性证据，但两份实现相同不能独立证明经验拟合符合真实电解液。

- Common operating domain (range of each input, with the basis):
  建议共同抽检区域：\(c\in[500,1500]\,\mathrm{mol/m^3}\)（0.5–1.5 mol/L），\(T\in[278.15,318.15]\,\mathrm K\)（5–45 °C）。依据是本 PR 的三个测试点覆盖的浓度、温度端点，以及 Ai2020 的默认初始浓度 1000 mol/m³、初始温度 298.15 K；矩形内的其他组合属于本审计建议补充的抽检，不是 PR 已测组合，也不是已经核实的文献适用区间。此区间 \(q\in[41.65,86.65]\,\mathrm K\)，远离分母奇点。排除负浓度、摄氏度直接输入、\(q\leq0\) 和趋近奇点的外推；不要求接口对这些非法物理输入一定抛异常。由单调性还可检验 \(D_A(1500,278.15)\leq D_A(c,T)\leq D_A(500,318.15)\)，两个端点应由 primary 的独立公式计算。

- Numerical check method (which public API of the merged commit, what to compute, what to compare against, suggested tolerance and why; if finite differences are needed, the step-size strategy):
  以下是供后续执行的方案，本次未运行任何代码。使用合并版本的 `p = pybamm.ParameterValues("Ai2020")`，通过 `p.evaluate(p["Electrolyte diffusivity [m2.s-1]"](pybamm.Scalar(c), pybamm.Scalar(T)))` 得到标量值；以 `ParameterValues("Xu2019")` 同样求值。该调用方式由合并后的测试直接确认，不需要建立 DFN 或运行放电求解器。

  在矩形四角、默认点、PR 三个点及固定种子的至少 100 个均匀随机点求值。对 primary，用独立的 `exp(log(10)*...)` 形式和 mol/L 输入转换计算基准，不调用被测函数作为 oracle；先断言有限、严格为正，再比较相对误差 \(\lvert D_A/D_{\rm ref}-1\rvert\leq10^{-12}\)。用相同相对容差、显式 `atol=0` 比较 Xu2019；小且平滑的双精度表达式在该区间没有病态分母，`1e-12` 足以容纳运算顺序舍入，同时能严格检查 SI 标度。舍入检查点 `3.2227e-10` 只用 `rtol=1e-3, atol=0`，不要赋予它完整公式的精度。

  对第二条，令 \(L(c,T)=\ln[D_A(c,T)/(1\,\mathrm{m^2/s})]\)，以中心差分分别计算 \([L(c,T+h_T)-L(c,T-h_T)]/(2h_T)\) 和浓度对应式。初始取 \(h_T=0.01\,\mathrm K\)、\(h_c=0.1\,\mathrm{mol/m^3}\)，再减半；两个步长结果均与解析导数比较，建议 `rtol=1e-5, atol=0`，并核对正负号。该区间最小分母为 41.65 K，温度中心差分的主截断相对误差约为 \((h_T/q)^2<6\times10^{-8}\)，故该容差有充足余量。仅在扰动点仍在矩形内时做中心差分；边界改用向内点或只做函数值比较，若步长缩小时误差不稳定，应记录并检查消减误差，不能放宽标度测试。符号断言也可用相隔至少 1 K 或 10 mol/m³ 的点对直接检查单调性。

- Do the PR's own tests already check this invariant (yes/partial/no, naming test files and assertions):
  **Primary: partial。** [合并后的 `test_Ai2020.py`](https://raw.githubusercontent.com/pybamm-team/PyBaMM/343ece2062dbb39789dbfa977e0f2daf96e9e370/packages/pybamm/tests/unit/test_parameters/test_parameter_sets/test_Ai2020.py) 中 `TestAi2020.test_functions` 对 `(1000,298.15)` 的扩散率断言 `pytest.approx(3.2227e-10, rel=1e-3, abs=0)`，确实约束物理量级，能识别旧的 `1e4` 倍错误，但未覆盖区域内函数形状。它也新增电导率 `1.1943` 的单点断言；电导率实现没有行为改动，故不额外列为本审计的不变量。

  **正性、响应导数与区域界限: no。** 新增测试没有显式检查这些性质。**跨参数集一致性: partial。** `TestAi2020.test_electrolyte_diffusivity` 在 `(1000,298.15)`、`(500,318.15)`、`(1500,278.15)` 调用 `np.testing.assert_allclose(..., rtol=1e-12)`，但没有设置 `atol`。NumPy 的默认 `atol=1e-8` 大于被测扩散率，因而不能把这项断言解读成真实的 `1e-12` 相对一致性检验；例如零值也可通过该项比较。旧的约 `1e-6` 数值仍会失败，且上述独立单点测试确实防止默认点退化为零。这里只评价所读 PR 测试，未声称审查了整个测试套件。

- Do you expect the change to satisfy it (satisfies / may not satisfy / uncertain), with one or two sentences of reasoning (speculation, not a conclusion):
  **satisfies（推测，尚未经数值执行验证）。** 修改后的公式包含正确的 cm²/s 到 m²/s 因子，与 Xu2019 的 mol/L 写法可代数化为同一表达式，且在指定区域内直接推出上述正性和响应方向；对真实电解液及区域外的经验准确性仍未验证，PR 描述中的最终电压或浓度范围不能当成普遍物理不变量。

- Environment: merge commit SHA, language and dependencies needed, estimated setup difficulty (low/medium/high):
  合并提交：`343ece2062dbb39789dbfa977e0f2daf96e9e370`。Python；[合并提交的包配置](https://raw.githubusercontent.com/pybamm-team/PyBaMM/343ece2062dbb39789dbfa977e0f2daf96e9e370/packages/pybamm/pyproject.toml) 要求 Python `>=3.10,<3.16`、NumPy `>=2.0.0`、SciPy `>=1.11.4,!=1.18.0`、`pybammsolvers>=0.11.0,<0.12.0`，并声明 xarray、anytree、sympy、typing-extensions、pandas、pooch、posthog、pyyaml、platformdirs、black 等运行依赖；构建元数据使用 hatchling/hatch-vcs。只需完整安装的 PyBaMM 参数入口及随包数据即可作标量检查，pytest 仅在采用其断言或原测试时需要。预计设置难度 **medium**：应固定该提交且保留参数集注册元数据、CSV 数据及兼容的原生 solver 依赖；实际检查不调用求解器，计算量很低。本次只读取输入与 GitHub 源码并编写规范，没有克隆、安装、编译或执行测试。

- Formalization (Lean) feasibility: can the invariant be written as a proposition over finite-dimensional reals, and what are the main obstacles:
  **可以。** 把 \((c,T)\in\mathbb R^2\)、精确十进制常数作为有理数、\(10^x\) 写成 `Real.exp (Real.log 10 * x)`，并假设 \(q>0\)，即可形式化 SI 因子、两种浓度表达的相等、正性、偏导符号及矩形端点界限。主要证明工作是非零分母、指数/对数与微分规则、严格不等式；数值舍入检查点需要区间界而非等式。物理单位需要显式量纲封装或单位转换定义；经验拟合的实验正确性无法由实数代数证明，且实数公式证明不自动覆盖 PyBaMM 符号求值路径或 IEEE-754 舍入，后者需要另做实现关联与误差证明。
