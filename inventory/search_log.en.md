# Inventory Round One (pass 1) Search and Preliminary Assessment Log: Inventory of Merged PRs Attributed to AI Agents and Preliminary Assessment of “Whether Physical Equations Were Changed”

> Executed by: LLM agent (Claude Sonnet), round one. Execution date 2026-10-16. Read-only access to public GitHub information (`gh api`), with no write operations.
> For scripts and search state, see `search_scripts/` (`search.py`/`search2.py`/`search3.py`/`build.py`/`triage.py`/`finalize.py`, state `state2.json`/`state3.json`); the first batch search state `state.json`, API cache, and patch cache `patches/` are not included in the package.

## 1. Summary

- Candidate PRs (deduplicated, merged_at within 2025-01-01 to 2026-10-16, attribution verified): **1517**, across 19 repositories.
- Preliminary assessment physics_change=yes: **68**, across **7 repositories** (CoolProp, Oceananigans, moose, PyBaMM, ClimaAtmos, Cantera, pandapower).
- The assessment basis falls into three categories: 209 with patches read; 403 containing only documentation/tests/CI/build/examples/configuration files (automatically excluded); and 905 containing source files but excluded based on titles and file paths, without reading patches line by line (an additional round of quick diff spot-checks was performed on PRs with small changes among these, recovering several omissions; this remains the largest source of uncertainty in this assessment track).
- Of the 31 yes cases in CoolProp, 16 carry the “[Data-level]” label (changes to EOS, viscosity/thermal conductivity correlations, melting lines, etc. in existing fluid JSON, rather than changes to C++ formulas); if only code-level equations are wanted, these 16 should be excluded first.
- Another 15 yes cases are marked “Borderline case” in reason, involving difficult-to-assess situations such as solver branches, boundary masks, and geometric metrics; the second independent blind assessment is most likely to disagree here.

### Candidate and Assessment Counts by Repository

| Repository | Candidates | yes | Patch read | Automatically excluded (no source code) | Path+title exclusion |
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

## 2. Search Method

- Each query for each repository: serial calls to `search/issues` (PRs) or `search/commits`, approximately 2.2 s apart, in the event of 403/429/secondary rate limiting, sleeping for 65 s before retrying (no rate-limit retries were actually triggered; retries was recorded as 0). Results were saved to disk individually, allowing resumption.
- Author account searches (`is:pr is:merged merged:2025-01-01..2026-10-16 author:X`): `app/devin-ai-integration`, `app/copilot-swe-agent`, `app/chatgpt-codex-connector`, `app/google-labs-jules`, `openhands-agent`, `app/claude`, `sweep-ai`, `app/sweep-ai-deprecated`.
- Body attribution searches: `"Generated with Claude Code"`, `"Co-Authored-By: Claude"`, `"claude.ai/code"`, `"Devin run"`, `"chatgpt.com/codex"`, `"Codex"`. After a hit, attribution was rechecked using regular expressions in the PR body (including the title); if the body contained no attribution, the PR commit messages were then fetched for verification (`pulls/N/commits`); those still lacking attribution were discarded (103 entries, see Section 5).
- Commit trailer searches (`search/commits`): `"Co-Authored-By: Claude"`, `"copilot-swe-agent"`, `"devin-ai-integration"` (committer-date 2025-01-01..2026-10-16). After verifying that the commit message contained attribution, `repos/O/R/commits/sha/pulls` was used to map it to PRs, retaining only PRs whose merged_at fell within the time window.
- Attribution verification regular expressions: Claude (`Generated with [Claude Code`, `Co-Authored-By: Claude`, `claude.ai/code`, `noreply@anthropic.com`), Devin (`Devin run`, `app.devin.ai`, `devin.ai/sessions`), Codex (`chatgpt.com/codex`, `Codex task`, `chatgpt-codex-connector`), Copilot (`copilot-swe-agent`, `Co-authored-by: Copilot`), Jules. Only attribution and accounts were considered; no inference was made from writing style.
- `pulls/N/files` (with patches) was fetched for each candidate; “source files” refers to the files remaining after excluding documentation, tests, CI, examples, scripts, configuration, and build files (JSON under CoolProp's `dev/fluids|mixtures|cubics|pcsaft|incompressible` was retained as source data).

## 3. Hit Counts for Each Query in Each Repository (search total_count, Before Attribution Verification)

Column meanings: the first 8 columns are author account searches; `T:` indicates PR body searches; `C:` indicates commit searches; `—` indicates not executed; `ERR` indicates that the repository name could not be found (see Section 4).

| Repository | devin | copilot | codex-bot | jules | openhands | claude-app | sweep | sweep-dep | T:GenCC | T:CoAuth | T:claude.ai/code | T:DevinRun | T:codexlink | T:Codex | C:CoAuth | C:copilot | C:devin |
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

## 4. Failures, Retries, and Uncovered Items

- Total calls: approximately 841 search calls (including 157 HTTP 422 responses: invalid author names or invalid repository names, see below); approximately 3800 core API calls (`pulls`, `commits/.../pulls`, `pulls/N/files`, `pulls/N/commits`). Retries: 0 rate-limit retries (elapsed time approximately 1 hour 45 minutes, mostly serial waiting).
- Invalid repository names (names in the search list could not be found):
  - `TransiEnt-Library/TransiEnt` returned 404 (the organization does not exist). The project of the same name is at `ClaRaLibrary/TransiEnt`, is archived, was last pushed to in 2022-06, and had no activity within the time window; it was not searched. **This repository was not covered**.
  - `NREL/BLAST-Lite`, `NREL/ssc`: the organization has been renamed to `NatLabRockies`; the search API does not follow redirects. Supplementary searches using `NatLabRockies/BLAST-Lite` and `NatLabRockies/ssc` yielded 0 hits. `NREL/EnergyPlus` was already listed as `NatLabRockies/EnergyPlus` in the original list.
  - `SINTEF/thermopack` has moved to `thermotools/thermopack`; the supplementary search yielded 0 hits.
  - `SanPen/GridCal` has been renamed to `SanPen/VeraGrid`; the supplementary search yielded 0 hits.
- Invalid author qualifiers: `author:Copilot` and `author:app/sweep-ai` directly returned 422 from the search API (usernames not searchable). Copilot was changed to `app/copilot-swe-agent`; sweep-ai was changed to `sweep-ai` and `app/sweep-ai-deprecated`. These two categories had no hits in 4 repositories (the Copilot account appeared 4 times among the candidates, none involving physical changes).
- Repositories with 0 hits (all searches returned 0, and the repository names were verified as searchable): all Modelica repositories other than ModelicaStandardLibrary (modelica-buildings, modelica-ibpsa, AixLib, BuildingSystems, ThermoPower, ClaRa-official), PowerModels.jl, idaes-pse, the four CalebBell repositories (thermo/fluids/ht/chemicals), windpowerlib, sfepy, OpenFOAM-dev. These repositories have no AI PRs with attribution traces, or attribution leaves no trace (see Section 6).
- Search result limit: at most the first 1000 results were retrieved per query; the largest single query this time was the CoolProp commit search with 478, which did not reach the limit.
- The `"Codex"` search produced many noisy hits (moose 67, Cantera 20, etc.); bodies were checked individually, and those without a `chatgpt.com/codex` link or Codex attribution were discarded.

## 5. Discarded Hits

- Text search hits with no attribution in either the body or commit messages: 101 entries; commit hits with no attribution in the commit message or no PR: 2 entries (`unverified.json`, not included in the package).
- PRs mapped from commits were discarded if merged_at fell outside the time window.

## 6. Methodological Limitations (Affecting the Interpretation of N_phys)

1. **Almost all candidates are PRs “written by humans, with Claude Code attribution”**: there were only 4 author account hits (Copilot, all documentation/tests/UI-related, with no physical changes); Devin, Codex bot, Jules, OpenHands, sweep-ai, and Claude app accounts had none in the 19 repositories with hits. The inclusion criteria allow such PRs to be counted, but “written autonomously by an AI agent” differs from “written by a human using Claude Code”; this distinction requires attention when interpreting the results.
2. 530 candidates were mapped to PRs solely through the `Co-Authored-By` trailer of merge commits; the evidence is that a commit within the PR carries attribution, which does not necessarily cover the entire PR; pandapower 3041/3054/3076/3107 are develop→release branch merge PRs, with attribution originating from individual merged commits, and were assessed as no on the basis of “uncertain attribution.” If similar branch synchronization PRs exist in other repositories, they may also be mixed into the no cases outside yes, or erroneously counted among yes cases; they were not investigated individually.
3. Patches were not read line by line for the 905 “path+title exclusions.” I went through all titles in CoolProp, Oceananigans, moose, PyBaMM, and ClimaAtmos, and performed quick diff scans of unread PRs with changes of ≤10–14 lines (thus recovering CoolProp 3029/2771/2784/3326/3335/3399), but PRs with neutral titles and larger changes may still contain changes to physical terms, posing a risk of missed detections.
4. AI PRs without attribution traces (users deleting attribution, trailers removed after squashing, tools not on the list such as Cursor/Aider/Copilot chat) cannot be detected; this inventory includes only those retaining attribution.
5. Trade-offs in the criteria (for comparison with the second independent blind assessment): solver robustness/root bracketing/phase determination/input validation/new functionality (previously raised an error, now can compute) were assessed as no; formulas, constants, signs, unit conversions, conservation bookkeeping, closure correlations, and modifications to existing fluid JSON data were assessed as yes; new pure data (new fluids, new mixture pairs) were assessed as no; time stepping and clocks were assessed as no.
