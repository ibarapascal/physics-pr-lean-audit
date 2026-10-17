## CoolProp/CoolProp#3339 fix(transport): rhosr-CS viscosity ignored x_crossover from the fluid file

- What changed (1–3 sentences, naming files and formulas):
  `src/Backends/Helmholtz/TransportRoutines.cpp` 中 `viscosity_rhosr` 把液相分支权重从 `1/(1+exp(-100*(x-2)))` 改成 `1/(1+exp(-100*(x-x_crossover)))`，读取已经由 JSON 解析的切换中心，并将固定陡峭度 100 命名为局部常量。`src/Tests/CoolProp-Tests.cpp` 新增通过公共 JSON 注册接口改变切换中心的回归测试；分支多项式、残余熵定义和最终黏度缩放均未改变。

- Invariants (mathematical statements; write them as checkable equalities/inequalities, give variables, sign conventions, domain of validity):
  统一记号：纯流体 HEOS 状态的温度为 \(T>0\) [K]，摩尔密度为 \(\rho>0\) [mol/m³]，无量纲残余 Helmholtz 能为 \(\alpha^r\)，\(\tau=T_r/T\)。采用有符号的摩尔残余熵
  \[
  s^r=R(\tau\alpha^r_\tau-\alpha^r)\quad[\mathrm{J/(mol\,K)}],\qquad
  x=\rho s^r/q_c,
  \]
  其中 \(q_c\) 是 JSON 的 `rhosr_critical` [J/(m³ K)]，必须非零。不要额外对残余熵取负号；在选定流体常见致密态，\(s^r\) 与 \(q_c\) 都是负数，\(x\) 为正。令 \(a=\texttt{x_crossover}\)、\(k=100\)，二者对应无量纲 \(x\)；\(L(x)=\sum_{j=0}^3 c_{L,j}x^j\)、\(V(x)=\sum_{j=0}^3 c_{V,j}x^j\)，\(C\) 是无量纲缩放系数，\(\eta_0(T)>0\) 是同一流体的稀薄气体动黏性系数以外的**动力黏度** [Pa·s]。以下 \(\eta\) 均指动力黏度。

  1. **Primary — 残余熵坐标下的对数混合及中心定位。** 对有限且有效的模型参数与单相状态，定义 \(w=1/(1+e^{-k(x-a)})\)，应有
     \[
     \eta(T,\rho;a)=\eta_0(T)\{1+C[\exp(wL+(1-w)V)-1]\}.
     \]
     特别地，\(a=x\) 时，归一化参考黏度必须为两分支的几何平均：
     \[
     Z:=1+\frac{\eta/\eta_0-1}{C}=e^{(L+V)/2},\qquad C\ne0.
     \]
     这不是最终动力黏度两端值的几何平均，因为最后还有仿射缩放。远离中心时，\(x-a\to+\infty\) 选择液相分支，\(x-a\to-\infty\) 选择气相分支；这里的极限在固定状态下通过移动参数 \(a\) 定义，不要求把物理状态推到无限密度。

  2. **Secondary — 无超调与有条件的正耗散。** 令
     \(\eta_L=\eta_0[1+C(e^L-1)]\)、\(\eta_V=\eta_0[1+C(e^V-1)]\)。对于 \(C>0\)，应有
     \[
     \min(\eta_L,\eta_V)\le\eta(a)\le\max(\eta_L,\eta_V).
     \]
     若两端均为正，则 \(\eta(a)>0\)，符合剪切黏性耗散非负的要求；尤其 \(0<C\le1\) 时正性自动成立。不能对任意用户 JSON 系数断言正性：\(C>1\) 时还需 \(\min(e^L,e^V)>1-1/C\)。本检查选用 R1233zd(E) 的 \(C=0.8089\)，因此无需额外假设。

  3. **Secondary — 切换中心的定向、平滑响应。** 固定 \(T,\rho\) 与其他系数，仅改变 \(a\)，有
     \[
     \frac{\partial\log Z}{\partial a}=-k\,w(1-w)(L-V),\qquad
     \left.\frac{\partial\log Z}{\partial a}\right|_{a=x}=-25(L-V).
     \]
     因此若 \(L>V\)，增加中心必须使黏度减小；若 \(L<V\)，方向相反；若 \(L=V\)，改变中心不能改变黏度。这是对参数的偏导，不是沿状态路径对温度或密度的单调性断言。不得假设标为“液相”的多项式在所有 \(x\) 上都较大。

- Sources of the invariants (first principles / literature formula / consistency with other implementations in the library / dimensions):
  残余熵关系来自 Helmholtz 热力学恒等式 \(s=-(\partial a/\partial T)_\rho\)；\(\rho s^r\) 与 \(q_c\) 同量纲，所以指数中的 \(x-a\) 无量纲。正性与剪切耗散的联系来自牛顿流体的 \(2\eta D':D'\ge0\)。具体混合模型属于经验关联式，不由热力学唯一决定：上述中心、边界和导数是该关联式必须满足的数学后果，不能单凭它们证明实验准确性。

  合并版本的 [TransportRoutines.cpp](https://raw.githubusercontent.com/CoolProp/CoolProp/5418cec58c52183e3a010b2894e08b818ed4cc81/src/Backends/Helmholtz/TransportRoutines.cpp) 给出混合、最终缩放和稀薄气体公式；[FluidLibrary.h](https://raw.githubusercontent.com/CoolProp/CoolProp/5418cec58c52183e3a010b2894e08b818ed4cc81/src/Backends/Helmholtz/Fluids/FluidLibrary.h) 的 `parse_rhosr_viscosity` 确认全部参数来自 JSON。[R1233zd(E) 数据](https://raw.githubusercontent.com/CoolProp/CoolProp/5418cec58c52183e3a010b2894e08b818ed4cc81/dev/fluids/R1233zd%28E%29.json) 的文献键为 `Bell-PURDUE-2016-ETA`，[库内书目](https://raw.githubusercontent.com/CoolProp/CoolProp/5418cec58c52183e3a010b2894e08b818ed4cc81/CoolPropBibTeXLibrary.bib) 将其标为 Bell 与 Laesecke，2016，*Viscosity of refrigerants and other working fluids from residual entropy scaling*。本审阅未读取论文正文，故不将固定常数 100 或可任意更改中心的物理有效性宣称为已由原论文验证。

- Common operating domain (range of each input, with the basis):
  首选纯 R1233zd(E)，以及只修改中心和身份字段的克隆，保留 EOS 与四项分支系数。数据为 \(q_c=-48503.78752999625\) J/(m³ K)、\(C=0.8089\)、默认 \(a=2\)，\(c_L=(0.6100913843,0.4508958312,-0.017063,0.000564)\)，\(c_V=(0,1.2,-0.308598,0.035317)\)。

  建议抽样域为 \(250\le T\le433\) K、\(1\le\rho\le10000\) mol/m³（密度按对数抽样并增加致密点），只接受公共 API 判为稳定单相、\(0<p<100\) MPa、\(0\le x\le10\) 且分支指数有限的点，不强制相态来绕过两相区。该矩形是审计采样范围，不是每一点都有效的模型适用域；EOS 数据给出三相点 165.75 K、温度上限 450 K、压力上限 100 MPa。四个现有回归点中可直接选取 \((350,100),(300,10000),(400,8000)\)，单位依上；其单相性在实际检查时仍应确认。

  更接近实验约束的液体子域是 313.67–433.32 K、1.005–4.074 MPa，数据注释称 C 在这一区间的 61 个 Miyara 等人液体测量点上拟合，实验标准不确定度 3.0%；这只是源码记录，未独立核验实验表格。300 K、10000 mol/m³ 是约 18.2 MPa 的回归外推点，不应称为实验验证点。

  中心参数采样 \(a=x+d\)，\(d\in[-0.4,0.4]\)，重点取 \(0,\pm0.001,\pm0.01,\pm0.05,\pm0.2,\pm0.4\)，另保留 PR 的 \(a=2,12\) 比较。前者控制指数在 ±40 内并覆盖过渡；任意改动中心的克隆只用于关联式结构测试，不当作真实流体的新标定。默认中心为 2 的单一流体点不能检验数据驱动的中心是否生效。

- Numerical check method (which public API of the merged commit, what to compute, what to compare against, suggested tolerance and why; if finite differences are needed, the step-size strategy):
  建议后续使用 C++ 公共 API，声明可核对 [include/CoolProp/AbstractState.h](https://raw.githubusercontent.com/CoolProp/CoolProp/5418cec58c52183e3a010b2894e08b818ed4cc81/include/CoolProp/AbstractState.h)。本次仅给出方案，未运行。

  1. 用 `get_fluid_param_string("R1233zd(E)", "JSON")` 取得单元素数组。按 PR 测试修改 `INFO.NAME`、`INFO.CAS` 为唯一标识，清空 `ALIASES`、设 `REFPROP_NAME="N/A"` 并移除测试列出的身份字段；使用 `add_fluids_as_JSON("HEOS", json_text)` 注册。每种中心用新克隆，不能覆盖现有流体或复用旧 `AbstractState` 对象。
  2. 通过 `AbstractState::factory("HEOS", name)`、`update(DmolarT_INPUTS,rho,T)` 建状态，以 `phase()` 和 `p()` 筛选。由 `gas_constant()`、`tau()`、`dalphar_dTau()`、`alphar()` 重建 \(s^r\) 和 \(x\)，也可与 `smolar_residual()` 核对。以 `viscosity()` 得到 \(\eta\)；高层 `PropsSI("V","T",T,"Dmolar",rho,name)` 应给出相同量和单位。
  3. 为避免调用内部稀薄气体函数，注册一个额外的 **C=0 控制克隆**，其余物性参数保持一致，在同一状态调用 `viscosity()` 得到 \(\eta_0\)。这是公共 API 可实施的隔离控制，不是物理标定；仅在分支指数有限的域使用，避免 0 乘无穷。它隔离本 PR 的混合行为，不独立验证原有稀薄气体模型。
  4. 独立计算多项式与 logistic（可用稳定的分段实现），算出不变量 1 的完整 \(\eta_{pred}\)，与每个克隆的 \(\eta_{API}\) 比较：\(|\eta_{API}-\eta_{pred}|\le10^{-10}\max(|\eta_{pred}|,\eta_0)\)。另用 \(y=\log Z\) 检查中心等式与分支界，绝对误差界取 \(10^{-10}\max(1,|L|,|V|)\)。这些是双精度公式一致性容差，不是实验误差；远小于百分比但允许 EOS 求值、JSON 往返及指数舍入差异。取 \(d=\pm0.4\) 时，更严谨的尾端界是 \(|y-L|\le|L-V|e^{-40}\) 或 \(|y-V|\le|L-V|e^{-40}\)，再加上述数值误差，不要求浮点下严格到达极限。
  5. 检查不变量 2 时将同一尺度的绝对容差加入端点界，正性单独要求 \(\eta>0\) 且有限。不要用日志归一化后的正性替代最终动力黏度的正性。
  6. 对不变量 3，在固定状态为 \(a\pm h\) 各注册克隆，计算 \([y(a+h)-y(a-h)]/(2h)\)。中心附近从 \(h=10^{-4}\) 起，依次减半到约 \(10^{-6}\)，寻找二阶误差下降及稳定平台；以 \(10^{-5}\max(1,|\partial_a y|)\) 为建议导数误差界。步长依据切换宽度 \(1/k=0.01\)，过大将平均过渡，过小放大舍入；在饱和尾端不用相对导数误差。\(|L-V|\) 太小时仅做绝对误差检查，勿由噪声判定单调方向。

  不额外强加任意中心下精确的 \(\eta\to\eta_0\) 稀薄极限：这套系数 \(L(0)\ne0\)，有限 \(a\) 时 logistic 液相权重严格非零；默认 \(a=2\) 的偏差极小，但移到零附近会使该偏差明显。这是经验模型参数物理适用性的限制，不等同于本次中心读取修复出错。

- Do the PR's own tests already check this invariant (yes/partial/no, naming test files and assertions):
  **整体 partial。** 新增 `src/Tests/CoolProp-Tests.cpp` 的 `rhosr-CS viscosity honours x_crossover from the fluid file` 要求类型为 `rhosr-CS`、默认中心为 2、克隆注册不抛异常，检查三个黏度有限，要求 `abs(clone_same/shipped-1)<1e-14` 与 `abs(clone_moved/shipped-1)>0.01`。这部分覆盖了不变量 1 的“中心参数确实生效”，但不检查预测数值、等权中心、混合方向或陡峭度；不变量 2 与 3 均 **no**，有限性不等于正性或无超调。

  同文件已有 `R1233zd(E) has a viscosity model (#3330)` 检查一个饱和液态值在 \(10^{-4}\) 与 \(10^{-3}\) Pa·s 之间，以及四个固定点相对误差低于 \(10^{-6}\)，能保护默认参数的回归行为，仍不能验证可变中心的上述性质。[合并版本测试源码](https://raw.githubusercontent.com/CoolProp/CoolProp/5418cec58c52183e3a010b2894e08b818ed4cc81/src/Tests/CoolProp-Tests.cpp)。输入 PR 所称全套测试通过不是本次执行结果。

- Do you expect the change to satisfy it (satisfies / may not satisfy / uncertain), with one or two sentences of reasoning (speculation, not a conclusion):
  **satisfies（预期，未执行验证）**：对上述有限、单相且正端点的域，合并实现直接采用数据中的中心，余下结构与三项不变量的代数推导一致。任意中心的实验可靠性仍 **uncertain**；两克隆黏度相差很大只证明参数有效，不能证明新参数是正确物理标定。

- Environment: merge commit SHA, language and dependencies needed, estimated setup difficulty (low/medium/high):
  合并 SHA：`5418cec58c52183e3a010b2894e08b818ed4cc81`。建议 C++ 检查驱动，依赖此提交构建的 CoolProp HEOS 库及公共头文件、可用 C++ 编译器与构建系统，JSON 克隆可使用项目已有的 nlohmann/json；若运行原测试还需项目 Catch2 测试配置及构建依赖。预计 **medium**：公开安装包不保证对应这次合并，必须固定源码版本及构建所需依赖；本次没有安装、编译或执行任何项目代码，也未实测构建环境。

- Formalization (Lean) feasibility: can the invariant be written as a proposition over finite-dimensional reals, and what are the main obstacles:
  **可以。** 将 \(x,a,C,\eta_0\) 及两组四维系数作为实数参数，三个命题都只需多项式、指数、对数、序与单变量微分；先证明 \(0<w<1\)，即可推出几何中心、无超调与带端点假设的正性。导数命题需先证明 \(Z>0\)，再用链式法则。

  主要障碍是把 EOS 的残余熵输出与该实数模型连接、证明指定状态处于有效单相域，以及验证 C++ 浮点、指数溢出、JSON 解析和 API 分派与实数语义之间的对应。可将 EOS 输出作为假设输入先证明混合层；这不能替代对实验准确性、EOS 全域正确性或任意自定义 JSON 安全性的证明。
