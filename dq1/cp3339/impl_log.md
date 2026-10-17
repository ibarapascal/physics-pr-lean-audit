# cp3339 实施日志（ΔQ1 统一抽查，CoolProp#3339）

> 目的：三路检查计划逐句→代码对应与全部实施判断，供协议 §6 核查与 §7 偏离登记。时间切片：2026-10-17 冒烟前写成；冒烟后只追加「冒烟后修正」节。
> 文件：`spotcheck.py`、`harness.py`（四份相同副本）、`job_smoke.json`、`job_full.json`。
> 读过的材料：协议；`track1.md`/`track2.md` 本 PR 节；`specs_codex/CoolProp_CoolProp_3339/spec.md` §2/§4/§5；合并提交 5418cec5、父提交 afcc7909、v8.0.0（ae81610e）的 `TransportRoutines.cpp::viscosity_rhosr`；合并提交 `src/Tests/CoolProp-Tests.cpp` 的 PR 测试；合并提交 `dev/fluids/` 下 8 个流体 JSON 的 TRANSPORT.viscosity 块；`env/setup_pca_py.sh`。

## 0. 环境（重要：没有原生合并提交环境）

- compute node A/compute node B 的 `env/py2` 是 CoolProp 8.0.0（ae81610e）。GitHub compare：ae81610e 是合并提交 5418cec5 的**祖先，落后 63 个提交，不含本 PR**；其 `viscosity_rhosr` 与父提交 afcc7909 逐字相同（中心写死为 2）。另外 8.0.0 的 `get_fluid_param_string(name,"JSON")` 输出没有 TRANSPORT 段（探针：只有 ANCILLARIES/EOS/INFO/STATES）。
- 结论：**节点上没有本 PR 的原生环境**，只能按合并提交源码构建。按共同说明没有自己去构建。env/py2 可作为「父提交等价」的基线（`--base-python env/py2/bin/python`），但跨构建比较不是同一构建的按位比较。
- 脚本在当前解释器上探测 `honours_x_crossover`（克隆 x_c=2 与 12，在 PR 测试状态比较），写进 JSON 的 `versions.env_probe`；若当前 CoolProp 是 8.0.0，第一路主2、第三路 I1 会大面积「拒绝」，那是环境不对，不是实测结果。
- JobSpec 假定合并提交构建装在 `env/cp3339`（不存在，需主会话决定是否构建：`uv pip install` 从 git 5418cec5 源码编译 CoolProp Python 扩展，C++17＋CMake，约十几分钟——推测）。
- 种子＝20261017＋7×10＋路号；判定/重放/失败规则同 pp3016 日志 §0。
- 克隆注册：仿 PR 测试（改 NAME/CAS、清空 ALIASES、REFPROP_NAME="N/A"、删 INCHI/SMILES/CHEMSPIDER/2DPNG 字段），`add_fluids_as_JSON("HEOS", …)`；开启 `OVERWRITE_FLUIDS`，每个用途用固定名字反复覆盖，避免进程内注册数随点数增长。
- x 由公开量算：x = ρ·R·(τ·∂α^r/∂τ − α^r)/rhosr_critical（`rhomolar/gas_constant/tau/dalphar_dTau/alphar`），rhosr_critical、C、c_liq、c_vap 读流体 JSON。

## 1. 第一路（track1.md）

分支：8 个流体（R1233zd(E)、R1234yf、R1234ze(E)、R124、R152A、R22、R245fa、R32），各 1250 点。必测点：T=300 K 时 x∈{x_c±0.01, ±0.05, ±0.5}（x_c∈{2,12}）对应的单相密度（二分求出；够不着的跳过）。

| 计划原句 | 代码 |
|---|---|
| T∈[230,420] K，ρ 从 1 mol/m³ 到饱和液密度，单相 | `t1_sampler`、`sample_rho`：T<Tc 时 ρ∈U[1, ρ_L,sat(T)]、落在两相区重抽；T≥Tc 时上限取 ρ_L,sat(230 K)（判断 H1） |
| 主1：PropsSI("V","T",T,"Dmolar",ρ,fluid) 父提交与合并提交各算一次，比较 float.hex | 基线进程（`--base-python`，`--serve-base` 模式）返回 PropsSI 的 hex；`main1.parent==merge bitwise` |
| 主2：add_fluids_as_JSON 注册 x_c=a、b 的克隆；x 由 smolar_residual() 与 rhomolar() 及 rhosr_critical 算；按 C、c_liq、c_vap 算理论比值，与 PropsSI 比值比，rtol 1e-12 | a=2（随库值）、b=12（PR 测试值）；`main2.eta_a/eta_b=…` |
| 网格 x∈{x_c±0.01, ±0.05, ±0.5} | 作为必测点 |
| 次：ρ=1e-3 mol/m³ 读 c_vap[0]，看 η/η_dilute→1（η_dilute 用 ρ 递减外推） | η₀=2η(1e-3)−η(2e-3)；`minor.c_vap[0]==0`、`minor.eta/eta_dilute->1`（默认容差） |

实施判断（第一路）：
- H1 T>Tc 时「到饱和液密度」无定义 → 上限取区间内最大的饱和液密度 ρ_L(230 K)。
- H2 计划写 x=−ρs^r/(ρs^r)_crit，同句又写 x=ρR(τα_τ−α)/(ρs_r)_crit；两式符号是否一致取决于 smolar_residual 的符号约定 → 以后一式（与源码一致）为准；每点也用 smolar_residual 算一遍，不一致时只记录。
- H3 主2 选 a=2、b=12 两个固定克隆（而不是每点按 x 生成克隆），网格点用必测点覆盖过渡区。

## 2. 第二路（track2.md）

分支：8 个流体，各 1250 点。

| 计划原句 | 代码 |
|---|---|
| T 从三相点到 1.5 T_c，ρ 从 1e-6 ρ_c 到饱和液密度，单相 | `t2_sampler`（T≥Tc 时上限取三相点饱和液密度，判断 H4）；超出 EOS Tmax 的点照算，出错计失败 |
| I1：合并提交与父提交在 (T,ρ) 上 PropsSI("V",…) 逐位相等 | 同第一路主1 |
| I2：克隆 a=2、b∈{1,3,6,12}；x≤min(a,b)−δ 或 ≥max(a,b)+δ（δ=0.2）时 \|η_a/η_b−1\| ≤ \|f_L−f_V\|·e^{−100δ} 量级 | 「量级」按 10 倍以内判：界＝10·\|f_L−f_V\|·e^{−20}（判断 H5） |
| I2：对任意 x，L≡ln(1+(η/η₀−1)/C) ∈ [min(f_L,f_V), max(f_L,f_V)]；η₀ 用 ρ→0 外推或 JSON 稀薄参数自算 | η₀ 取 ρ→0 线性外推（1e-6ρ_c、2e-6ρ_c 两点，判断 H6）；a=2 与各 b 都核；边界容差默认 |
| I3：ρ=ρ_c·{1e-3,1e-4,1e-5,1e-6} 取 η，线性外推到 0 与最小 ρ 值比 rtol 1e-6；读 c_V[0] | 外推用最小两点；斜率序列只记录（「随 ρ 线性趋近」为定性，判断 H7） |

## 3. 第三路（spec.md §5）

分支：8 个流体，各 1250 点。必测点：每流体一个 I4 稀薄点组（T∈{250,300,min(400,Tmax)} × ρ/ρ_c∈{1e-6,1e-8,1e-10,1e-12} × a∈{2,12,0}）；QT 的 Q=0/1 于 T/Tc∈{0.70,0.85,0.95}；超临界 T/Tc∈{1.01,1.05} × ρ/ρ_c∈{0.1,0.5,1,1.5,2,3}（在 Tmax 内）；PR 状态 R1233zd(E) 300 K、1e4 mol/m³。这些必测点同时做 PropsSI 交叉核对。

| 计划原句 | 代码 |
|---|---|
| 12 个温度 × 24 个从 1 kPa 到压力上限的对数压力（PT_INPUTS），实用区 250 K–min(450 K,Tmax)、0<p≤min(20 MPa,pmax)，ρ≲3ρ_c | 随机点：T 均匀、p 对数均匀（计划明写对数）；ρ>3ρ_c 重抽（判断 H8）；计划的规则网格改由均匀抽样替代（协议 §3） |
| I1：C=0 克隆给 η₀；a=x+d，d∈{−0.2,−0.05,−0.01,0,0.01,0.05,0.2}，另加 a=2、12；独立写的溢出安全 logistic；相对误差 ≤5e-12 | `t3_state` |
| I1：恢复 g=log[1+(η/η₀−1)/C]，g(x+d)+g(x−d)=L+V，绝对 1e-10·max(1,\|L\|,\|V\|) | d∈{0.01,0.05,0.2} |
| I1：a 导数中心差分 h=1e-5、5e-6，相对 1e-5；退化时绝对容差按 η₀max(1,CE\|L−V\|) 缩放 | 在 a=x 处；退化（\|L−V\|<1e-6 或 C=0）时绝对容差 1e-10·η₀·max(1,CE\|L−V\|)（系数 1e-10 为判断 H9） |
| I1：a=x±8 尾部探针，内层指数溢出应给权重 0 而非非有限黏度 | `I1.tail` |
| I2：基线与合并两构建比较 8 种流体；同构建要按位，否则相对 1e-13；随库流体对 a=2 克隆 1e-14 | 跨构建 1e-13 判定、按位相等只记录；克隆比较 1e-14 |
| I3：η 有限且为正；两端点界（slack 5e-12）；x≥0 时 η≥η₀；负 x 单独报告；等多项式克隆与 C=0 对照下改 crossover 无效 | `t3_state` I3 段；等多项式克隆令 c_liq=c_vap |
| I4：250/300/min(400,Tmax) K、ρ/ρ_c=1e-6…1e-12，a=2,12,0：与有限密度式比 5e-12；最后两点对极限式 1e-8；a=2,12 再核 η/η₀≈1；a=0 核非 1 极限 | `t3_dilute`；a=0 的「核非 1 极限」＝对极限式比，另记录比值 |
| 每个流体用新进程隔离注册 | 未照做：单进程、固定名覆盖注册（判断 H10，偏离候选） |

## 4. 耗时与内存估计（推测）

PropsSI/AbstractState 单次 ~1 ms；克隆注册（深拷贝 ~70 KB JSON＋解析建流体）~5–20 ms。

| 路 | 单点 | N=100（含必测） | N=10⁴ |
|---|---|---|---|
| 1 | ~15 ms＋基线往返 | ~5 s | ~3–5 min |
| 2 | ~30 ms | ~5 s | ~5–10 min |
| 3 | ~0.2–0.3 s（每点 13 次覆盖注册） | ~1–2 min | ~45–60 min |

峰值内存：主进程 ~200–400 MB＋基线进程 ~150 MB；覆盖注册若有泄漏会线性增长（未验证）→ 申报 2 GiB，冒烟时看峰值。全量串行。

## 5. 偏离候选（交主会话裁定）

- D0 **没有合并提交原生环境**；需构建 `env/cp3339`（合并提交）才能跑；基线用 env/py2（父提交等价的 viscosity_rhosr，但不同构建）。第一路主1/第二路 I1 要求「按位」，跨构建按位比较可能因编译差异产生假拒绝。
- D1 第三路规则网格→均匀/对数均匀抽样；单进程覆盖注册代替每流体新进程（H10）。
- D2 第二路 I2「量级」按 10 倍判（H5）。
