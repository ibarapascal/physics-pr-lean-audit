# 盘点第一轮（pass 1）检索与初判日志：AI agent 署名的已合并 PR 盘点与「是否改动物理方程」初判

> 执行：LLM 代理（Claude Sonnet），第一轮。执行日 2026-10-16。只读 GitHub 公开信息（`gh api`），无任何写操作。
> 脚本与检索状态见 `search_scripts/`（`search.py`/`search2.py`/`search3.py`/`build.py`/`triage.py`/`finalize.py`，状态 `state2.json`/`state3.json`）；第一批检索状态 `state.json`、API 缓存与补丁缓存 `patches/` 未随包。

## 1. 汇总

- 候选 PR（去重、merged_at 在 2025-01-01 至 2026-10-16 内、署名已核验）：**1517**，涉及 19 个库。
- 初判 physics_change=yes：**68**，分布在 **7 个库**（CoolProp、Oceananigans、moose、PyBaMM、ClimaAtmos、Cantera、pandapower）。
- 判定依据分三类：读 patch 的 209 份；仅文档/测试/CI/构建/示例/配置文件的 403 份（自动排除）；有源文件但按标题与文件路径排除、未逐行读 patch 的 905 份（对其中小改动 PR 另做了一轮 diff 速览抽查，补出若干漏项；仍是本路径最大的不确定来源）。
- 其中 CoolProp 的 31 个 yes 里有 16 个带「[数据级]」标记（修改既有流体 JSON 的 EOS、粘度/热导率关联式、熔化线等，不是 C++ 公式改动）；若只要代码层方程，应先剔除这 16 个。
- 另有 15 个 yes 在 reason 中标了「边界案例」，是求解器分支、边界 mask、几何度量等难以判定的情形，第二路盲判最可能与此处分歧。

### 每库候选与判定计数

| 库 | 候选 | yes | 读patch | 自动排除(无源码) | 路径+标题排除 |
|---|---|---|---|---|---|
| CoolProp/CoolProp | 491 | 31 | 107 | 157 | 227 |
| CliMA/Oceananigans.jl | 251 | 10 | 36 | 47 | 168 |
| SciML/ModelingToolkit.jl | 242 | 0 | 1 | 91 | 150 |
| FEniCS/dolfinx | 158 | 0 | 0 | 30 | 128 |
| idaholab/moose | 150 | 5 | 14 | 20 | 116 |
| pybamm-team/PyBaMM | 92 | 10 | 23 | 17 | 52 |
| firedrakeproject/firedrake | 36 | 0 | 0 | 10 | 26 |
| CliMA/ClimaAtmos.jl | 28 | 9 | 10 | 5 | 13 |
| SciML/ModelingToolkitStandardLibrary.jl | 23 | 0 | 2 | 18 | 3 |
| Cantera/cantera | 15 | 2 | 7 | 1 | 7 |
| e2nIEE/pandapower | 8 | 1 | 6 | 0 | 2 |
| PyPSA/PyPSA | 7 | 0 | 0 | 2 | 5 |
| pvlib/pvlib-python | 6 | 0 | 2 | 1 | 3 |
| CliMA/Thermodynamics.jl | 3 | 0 | 0 | 2 | 1 |
| NatLabRockies/EnergyPlus | 2 | 0 | 0 | 2 | 0 |
| devitocodes/devito | 2 | 0 | 0 | 0 | 2 |
| ClapeyronThermo/Clapeyron.jl | 1 | 0 | 0 | 0 | 1 |
| OpenMDAO/OpenMDAO | 1 | 0 | 0 | 0 | 1 |
| modelica/ModelicaStandardLibrary | 1 | 0 | 1 | 0 | 0 |

## 2. 检索方法

- 每库每条检索：串行调用 `search/issues`（PR）或 `search/commits`，间隔约 2.2 s，403/429/secondary 限流时休眠 65 s 重试（实际未触发限流重试，retries 记录为 0）。结果逐条落盘，可续跑。
- 作者账号检索（`is:pr is:merged merged:2025-01-01..2026-10-16 author:X`）：`app/devin-ai-integration`、`app/copilot-swe-agent`、`app/chatgpt-codex-connector`、`app/google-labs-jules`、`openhands-agent`、`app/claude`、`sweep-ai`、`app/sweep-ai-deprecated`。
- 正文署名检索：`"Generated with Claude Code"`、`"Co-Authored-By: Claude"`、`"claude.ai/code"`、`"Devin run"`、`"chatgpt.com/codex"`、`"Codex"`。命中后在 PR 正文（含标题）里按正则复核署名；正文不含署名的，再取 PR 的提交信息复核（`pulls/N/commits`）；仍无署名的丢弃（103 条，见第 5 节）。
- 提交 trailer 检索（`search/commits`）：`"Co-Authored-By: Claude"`、`"copilot-swe-agent"`、`"devin-ai-integration"`（committer-date 2025-01-01..2026-10-16）。复核提交信息含署名后，用 `repos/O/R/commits/sha/pulls` 映射到 PR，只保留 merged_at 落在时间窗内的 PR。
- 署名复核正则：Claude（`Generated with [Claude Code`、`Co-Authored-By: Claude`、`claude.ai/code`、`noreply@anthropic.com`）、Devin（`Devin run`、`app.devin.ai`、`devin.ai/sessions`）、Codex（`chatgpt.com/codex`、`Codex task`、`chatgpt-codex-connector`）、Copilot（`copilot-swe-agent`、`Co-authored-by: Copilot`）、Jules。只看署名与账号，不按文风推断。
- 每个候选取 `pulls/N/files`（带 patch）；「源文件」指排除文档、测试、CI、示例、脚本、配置、构建文件后剩下的文件（CoolProp 的 `dev/fluids|mixtures|cubics|pcsaft|incompressible` 下 JSON 视为源数据保留）。

## 3. 每库每条检索的命中数（search 的 total_count，未经署名复核）

列含义：前 8 列为作者账号检索；`T:` 为 PR 正文检索；`C:` 为提交检索；`—` 表示未执行；`ERR` 表示该库名查不到（见第 4 节）。

| 库 | devin | copilot | codex-bot | jules | openhands | claude-app | sweep | sweep-dep | T:GenCC | T:CoAuth | T:claude.ai/code | T:DevinRun | T:codexlink | T:Codex | C:CoAuth | C:copilot | C:devin |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| modelica/ModelicaStandardLibrary | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 |
| lbl-srg/modelica-buildings | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ibpsa/modelica-ibpsa | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| RWTH-EBC/AixLib | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| UdK-VPT/BuildingSystems | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| TransiEnt-Library/TransiEnt | ERR | ERR | ERR | ERR | ERR | ERR | — | — | ERR | ERR | ERR | ERR | ERR | ERR | ERR | ERR | ERR |
| casella/ThermoPower | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| xrg-simulation/ClaRa-official | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| SciML/ModelingToolkitStandardLibrary.jl | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 15 | 2 | 1 | 0 | 0 | 2 | 18 | 0 | 1 |
| SciML/ModelingToolkit.jl | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 159 | 6 | 90 | 3 | 9 | 29 | 310 | 0 | 12 |
| ClapeyronThermo/Clapeyron.jl | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| CliMA/Thermodynamics.jl | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 3 | 0 | 0 |
| CliMA/Oceananigans.jl | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 103 | 1 | 1 | 0 | 0 | 4 | 238 | 9 | 0 |
| CliMA/ClimaAtmos.jl | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 31 | 1 | 0 |
| lanl-ansi/PowerModels.jl | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| pybamm-team/PyBaMM | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 70 | 0 | 0 | 0 | 0 | 3 | 89 | 0 | 0 |
| e2nIEE/pandapower | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 3 | 5 | 0 | 0 |
| PyPSA/PyPSA | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 1 | 6 | 1 | 0 |
| IDAES/idaes-pse | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Cantera/cantera | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 20 | 128 | 0 | 0 |
| CoolProp/CoolProp | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 449 | 13 | 13 | 0 | 0 | 2 | 478 | 0 | 0 |
| CalebBell/thermo | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| CalebBell/fluids | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| CalebBell/ht | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| CalebBell/chemicals | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| pvlib/pvlib-python | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 1 | 0 | 0 | 0 | 6 | 0 | 0 |
| wind-python/windpowerlib | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| OpenMDAO/OpenMDAO | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| usnistgov/fipy | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| sfepy/sfepy | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| devitocodes/devito | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 0 | 2 | 2 | 0 | 0 |
| FEniCS/dolfinx | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 56 | 20 | 1 | 0 | 0 | 19 | 154 | 1 | 0 |
| firedrakeproject/firedrake | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 15 | 1 | 5 | 0 | 0 | 5 | 32 | 0 | 0 |
| NREL/BLAST-Lite | ERR | ERR | ERR | ERR | ERR | ERR | — | — | ERR | ERR | ERR | ERR | ERR | ERR | ERR | ERR | ERR |
| NREL/ssc | ERR | ERR | ERR | ERR | ERR | ERR | — | — | ERR | ERR | ERR | ERR | ERR | ERR | ERR | ERR | ERR |
| NatLabRockies/EnergyPlus | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 2 | 0 | 0 |
| SINTEF/thermopack | ERR | ERR | ERR | ERR | ERR | ERR | — | — | ERR | ERR | ERR | ERR | ERR | ERR | ERR | ERR | ERR |
| SanPen/GridCal | ERR | ERR | ERR | ERR | ERR | ERR | ERR | ERR | ERR | ERR | ERR | ERR | ERR | ERR | ERR | ERR | ERR |
| idaholab/moose | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 27 | 1 | 0 | 0 | 0 | 67 | 439 | 0 | 0 |
| OpenFOAM/OpenFOAM-dev | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| NatLabRockies/BLAST-Lite | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| NatLabRockies/ssc | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| thermotools/thermopack | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| SanPen/VeraGrid | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

## 4. 失败、重试与未覆盖项

- 总调用：search 约 841 次（其中 157 次为 HTTP 422：无效作者名或无效库名，见下）；core API（`pulls`、`commits/.../pulls`、`pulls/N/files`、`pulls/N/commits`）约 3800 次。重试：限流重试 0 次（用时约 1 小时 45 分，主要是串行等待）。
- 库名失效（检索清单里的名字查不到）：
  - `TransiEnt-Library/TransiEnt` 返回 404（该组织不存在）。同名项目在 `ClaRaLibrary/TransiEnt`，已归档，最后推送 2022-06，时间窗内无活动，未检索。**该库未覆盖**。
  - `NREL/BLAST-Lite`、`NREL/ssc`：组织已改名为 `NatLabRockies`，搜索 API 不跟随重定向。已用 `NatLabRockies/BLAST-Lite`、`NatLabRockies/ssc` 补查，命中 0。`NREL/EnergyPlus` 清单里原已写作 `NatLabRockies/EnergyPlus`。
  - `SINTEF/thermopack` 已迁至 `thermotools/thermopack`，补查命中 0。
  - `SanPen/GridCal` 已更名 `SanPen/VeraGrid`，补查命中 0。
- 无效作者限定：`author:Copilot`、`author:app/sweep-ai` 在搜索 API 上直接 422（用户名不可搜索）。Copilot 改用 `app/copilot-swe-agent`；sweep-ai 改用 `sweep-ai` 与 `app/sweep-ai-deprecated`。这两类在 4 个库里没有命中（Copilot 账号在候选里出现 4 次，均不是物理改动）。
- 命中为 0 的库（所有检索都为 0，已核实库名可查）：ModelicaStandardLibrary 之外的全部 Modelica 库（modelica-buildings、modelica-ibpsa、AixLib、BuildingSystems、ThermoPower、ClaRa-official）、PowerModels.jl、idaes-pse、CalebBell 四库（thermo/fluids/ht/chemicals）、windpowerlib、sfepy、OpenFOAM-dev。这些库没有署名留痕的 AI PR，或者署名不留痕（见第 6 节）。
- 搜索结果上限：单条检索最多取前 1000 条；本次最大单条为 CoolProp 提交检索 478，未触顶。
- 搜索的 `"Codex"` 命中噪声大（moose 67、Cantera 20 等），逐个核对正文，没有 `chatgpt.com/codex` 链接或 Codex 署名的丢弃。

## 5. 被丢弃的命中

- 文本检索命中但正文与提交信息都无署名：101 条；提交命中但提交信息无署名或无 PR：2 条（`unverified.json`，未随包）。
- 提交映射到的 PR 若 merged_at 不在时间窗内，丢弃。

## 6. 方法性局限（会影响 N_phys 的解读）

1. **候选几乎全是「人写、带 Claude Code 署名」的 PR**：作者账号命中仅 4 份（Copilot，均为文档/测试/UI 类，无物理改动）；Devin、Codex 机器人、Jules、OpenHands、sweep-ai、Claude app 账号在 19 个有命中的库里一份都没有。纳入判据允许这类 PR 计入，但「AI agent 自主写」与「人用 Claude Code 写」不同，解读时要注意。
2. 530 份候选只靠合并提交的 `Co-Authored-By` trailer 映射到 PR，证据是 PR 里某个提交带署名，不一定覆盖整个 PR；pandapower 的 3041/3054/3076/3107 是 develop→release 的分支合并 PR，署名来自被合并的个别提交，已按「归因存疑」判 no。其他库里同类的分支同步 PR 若存在，也可能混在 yes 之外的 no 里或在 yes 里被误计，未逐个排查。
3. 905 份「路径+标题排除」未逐行读 patch。CoolProp、Oceananigans、moose、PyBaMM、ClimaAtmos 的全部标题我都过了一遍，并对改动 ≤10–14 行的未读 PR 做了 diff 速览（因此补出 CoolProp 3029/2771/2784/3326/3335/3399），但标题中性、改动较大的 PR 里仍可能藏着物理项改动，属漏检风险。
4. 署名不留痕的 AI PR（用户删署名、squash 后 trailer 被去掉、Cursor/Aider/Copilot 聊天等不在清单的工具）无法检出；本盘点只纳入留署名的。
5. 判据上的取舍（供第二路盲判对照）：求解器鲁棒性/根括区间/相态判定/输入校验/新增功能（原先抛错现在能算）判 no；公式、常数、符号、单位换算、守恒记账、闭合关联式、修改既有流体 JSON 数据判 yes；新增纯数据（新流体、新混合对）判 no；时间步进与时钟判 no。
