# Third-party material

This archive contains fragments of, transcriptions of, and unmodified copies of source code from the open-source
libraries listed below. That material remains under its upstream license; the MIT license (`LICENSE`) and the
CC BY 4.0 license (`LICENSE-DATA`) of this archive do not apply to it. Full upstream license texts are in
`third_party_licenses/` (taken from each project's default branch) and, for the controlled-comparison snapshots,
in `controlled/tasks/LICENSES/` (taken at the snapshot commits).

## Libraries whose code is transcribed in the Lean modules and replay scripts

The Lean modules in `lean/PhysAudit/` transcribe the changed functions of the audited pull requests into Lean
definitions and quote the corresponding upstream source lines (with file and line references) in comments. The
replay scripts in `replay/scripts/` call the native libraries and, where noted, re-implement or embed upstream code.

| Library | Upstream | Copyright holder (as stated upstream) | License | Files in this archive |
|---|---|---|---|---|
| Cantera | https://github.com/Cantera/cantera | California Institute of Technology; Sandia Corporation; Cantera Developers | BSD-3-Clause (`third_party_licenses/Cantera_License.txt`) | `lean/PhysAudit/Cantera2107.lean`; `replay/scripts/Cantera2107_replay.py` (embeds a phase definition from Cantera's `test/data/thermo-models.yaml`) |
| pandapower | https://github.com/e2nIEE/pandapower | University of Kassel; Fraunhofer Institute for Energy Economics and Energy System Technology (IEE), Kassel; individual contributors | BSD-3-Clause (`third_party_licenses/pandapower_LICENSE.txt`) | `lean/PhysAudit/Pandapower3016.lean`; `replay/scripts/pandapower3016_replay.py`, `replay/scripts/yzn_base_invariance.py` |
| PyBaMM | https://github.com/pybamm-team/PyBaMM | the PyBaMM team | BSD-3-Clause (`third_party_licenses/PyBaMM_LICENSE.txt`) | `lean/PhysAudit/PyBaMM5888.lean`, `lean/PhysAudit/PyBaMM5889.lean`; `replay/scripts/PyBaMM5888_replay.py`, `replay/scripts/PyBaMM5889_replay.py` |
| ClimaAtmos.jl | https://github.com/CliMA/ClimaAtmos.jl | California Institute of Technology (Climate Modeling Alliance) | Apache-2.0 (`third_party_licenses/ClimaAtmos_LICENSE.txt`) | `lean/PhysAudit/ClimaAtmos4637.lean`, `lean/PhysAudit/ClimaAtmos4822.lean`; `replay/scripts/ClimaAtmos4637_replay.jl`, `replay/scripts/climaatmos4822_diff.jl` |
| Oceananigans.jl | https://github.com/CliMA/Oceananigans.jl | Climate Modeling Alliance | MIT (`third_party_licenses/Oceananigans_LICENSE.txt`) | `lean/PhysAudit/Oceananigans6024.lean`, `lean/PhysAudit/Oceananigans6083.lean`; `replay/scripts/Oceananigans6024_replay.jl`, `replay/scripts/Oceananigans6083_replay.jl` |
| CoolProp | https://github.com/CoolProp/CoolProp | Ian H. Bell and other CoolProp developers | MIT (`third_party_licenses/CoolProp_LICENSE.txt`) | `lean/PhysAudit/CoolProp2851.lean`, `lean/PhysAudit/CoolProp3339.lean`; `replay/scripts/CoolProp2851_replay.py`, `replay/scripts/CoolProp3339_replay.py` |
| MOOSE | https://github.com/idaholab/moose | Battelle Energy Alliance, LLC and other contributors (`third_party_licenses/MOOSE_COPYRIGHT.txt`) | LGPL-2.1 (`third_party_licenses/MOOSE_LICENSE.txt`) | `lean/PhysAudit/Moose33307.lean`, `lean/PhysAudit/Moose34020.lean` (quoted fragments) |

**LGPL-2.1 files.** `replay/scripts/Moose33307_replay.py` and `replay/scripts/Moose34020_replay.py` contain
line-by-line Python re-implementations of MOOSE C++ source (`SCMFrictionUpgradedChengTodreas.C`,
`SCMMixingChengTodreas.C`, `SCMHTCBorishanskii.C`, `SubChannel1PhaseProblem`,
`ComputeDynamicFrictionalForceLMMechanicalContact.C`, `ComputeDynamicWeightedGapLMMechanicalContact.C`).
These two files as a whole are distributed under LGPL-2.1, not under the MIT license of this archive. The same
applies to `dq1/mo33307/scm_port.py` and `dq1/mo34020/mortar_port.py`, line-by-line Python re-implementations of
the same MOOSE classes used by the ΔQ1 spot checks.

## Pull-request diffs given to specification track 3

`audit/specs/track3/<PR>/input/pr.diff` are unmodified diffs of the audited pull requests, and `input/pr.json`
holds their public GitHub metadata (title, body, commit messages and file list; e-mail addresses other than
no-reply addresses replaced by `<email>`). The diffs remain under the license of the library in the table above.

## Library snapshots used in the controlled comparison

`controlled/tasks/<task>/snapshot/` holds unmodified subsets of the library trees at the pre-change commit of each
task. `controlled/judge/controls/*.diff` are patches derived from the historical upstream commits, and
`controlled/judge/extra/buildings_3027/` contains one unmodified upstream model file. The generated samples in
`controlled/samples/` quote fragments of these libraries.

| Library | Upstream | Snapshot commits | License at the snapshot commits |
|---|---|---|---|
| Modelica Buildings Library | https://github.com/lbl-srg/modelica-buildings | 309c814f (buildings_2725), dff251c2 (buildings_3027), b83b2278 (buildings_4775) | Revised BSD-3-Clause license of the Modelica Association, IBPSA, The Regents of the University of California (Lawrence Berkeley National Laboratory) and contributors, with an added paragraph on enhancements; text of `Buildings/legal.html` at each commit in `controlled/tasks/LICENSES/modelica-buildings_<commit>_legal.html` (the three texts differ only in the copyright years) |
| Modelica Standard Library | https://github.com/modelica/ModelicaStandardLibrary | a3041258 (msl_1335), 7c8ccafa (msl_1530) | Modelica License 2 (Modelica Association); text of `Modelica/Resources/Documentation/ModelicaLicense2.html` at these commits (identical for both) in `controlled/tasks/LICENSES/ModelicaStandardLibrary_ModelicaLicense2.html`. The library has been licensed under BSD-3-Clause since version 4.0.0; the current upstream `LICENSE` is included for reference as `controlled/tasks/LICENSES/ModelicaStandardLibrary_LICENSE_BSD-3-Clause_current.txt` |

## Software used but not redistributed

- Lean 4 (toolchain `leanprover/lean4:v4.34.1`) and Mathlib (`v4.34.1`, commit pinned in `lean/lake-manifest.json`)
  with its dependencies (batteries, aesop, Qq, ProofWidgets, plausible, LeanSearchClient, importGraph, Cli) are
  fetched by `lake` and are not included. Lean and Mathlib are licensed under Apache-2.0.
- The native libraries above (at the commits and versions listed in `README.md`), OpenModelica, Julia and the
  Python packages installed by `replay/env/` are downloaded by the setup scripts and are not included.
- The pull-request metadata in `inventory/` and `audit/pr_provenance.json` (repository names, PR numbers, titles,
  account names, dates, commit identifiers, file paths and short signature excerpts) was retrieved from the public
  GitHub API.
