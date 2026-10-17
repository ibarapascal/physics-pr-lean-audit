# AI-attributed physics-equation changes in open-source physics libraries: inventory, invariant specifications, Lean transcriptions, and native replays

This dataset contains:

1. **Inventory** (`inventory/`): merged pull requests (PRs), merged between 2025-01-01 and 2026-10-16, in a fixed list of 40
   open-source physics and engineering libraries (44 repository names including renamed successors; one library
   could not be searched, see `inventory/search_log.md`) that carry an AI-agent attribution (agent account, PR-body
   signature or commit trailer), each classified in two independent passes as changing or not changing a physics
   equation.
2. **Audit of 12 changes** (`audit/`, `lean/`, `replay/`): a stratified sample of 12 changes that both passes
   classified as code-level physics-equation changes. For each one: three specification tracks of library-level
   invariants (Tracks 1 and 2 from two independent Claude agents, Track 3 from OpenAI Codex), a Lean 4
   transcription of the changed code with proofs or counterexamples, a native replay that compares the installed
   library with the transcription, and per-PR judgments.
3. **ΔQ1 spot checks** (`dq1/`): each track's initial check plan implemented as written and run natively with
   10^4 points on pandapower#3016 (the only PR with violations caught at the Lean stage), plus the scripts written
   for the other 11 PRs.
4. **Controlled comparison** (`controlled/`): five historical physics-equation changes in the Modelica Buildings
   Library and the Modelica Standard Library, re-posed as change requests on a snapshot of the pre-change library;
   10 generated patches per task from a coding agent; and the results of library regression runs and invariant
   spot checks (10^4 points) for each patch, for the pre-change library and for known wrong/correct reference patches.

Citation: see `CITATION.cff`.

## Directory structure

```
inventory/
  candidates.csv            1517 candidate PRs (19 libraries with hits): repo, number, url, title, author,
                            merged_at, signature_evidence, files_changed, query
  pass1_classify.csv        pass 1 (all 1517): physics_change yes/no, invariant_class, reason
  pass2_input.csv           pass-2 input (316 PRs): repo, number, url, title, files_changed
  pass2_classify.csv        pass 2, blind to pass 1 (316 PRs): physics_change, data_level, invariant_class, reason
  pass2_strata.csv          stratum of each pass-2 PR: p1_read (209, patch read in pass 1) or
                            p1_title_excluded_sample (107, sample of PRs excluded in pass 1 by path and title)
  merged.csv                both passes side by side (316 PRs); n_phys = 1 if both passes say yes (55),
                            data_level = 1 if the change only edits property data files
  search_log.md             pass-1 search and classification log (queries, counts, failures, limitations)
  search_scripts/           pass-1 scripts (search*.py -> build.py -> triage.py -> finalize.py, plus inspection
                            helpers), their logs, search state for the supplementary searches (state2/3.json),
                            per-PR decisions read from the patches (decisions.txt) and hit counts (hits_table.md)
audit/
  sample.csv                the 12 audited PRs
  pr_provenance.json        author, association, merger, reviewers/approvers and per-commit AI trailers
  specs/track1.md           invariant specification, track 1
  specs/track2.md           invariant specification, track 2 (written independently, blind to track 1)
  specs/compare.md          comparison of tracks 1 and 2
  specs/track12_prompt.txt  prompt given to tracks 1 and 2 (original, Chinese; English: track12_prompt.en.txt)
  specs/track3/<PR>/        track 3 (OpenAI Codex CLI, gpt-6-astra, reasoning effort medium, one run per PR, run after
                            the Lean stage with no audit conclusions in its inputs): spec.md (original output; English:
                            spec.en.md), codex.log (raw CLI log, not translated), input/pr.json and input/pr.diff (the
                            PR metadata and diff it received); prompt.txt (English translation of the track 1-2 prompt)
                            and run.sh
  specs/track3_guided/<PR>/ an earlier track 3 run whose prompt, written with knowledge of the audit results, also asked
                            for every vector group, phase, regime and option value (spec.md, prompt.txt); kept as a
                            control and not used for ΔQ1
  specs/compare3.md         comparison of the primary invariants and check plans of tracks 1-3 (unguided track 3)
  specs/compare3_guided.md  the same comparison with the guided track 3 run
  per_pr_judgments.csv      per-PR judgments (record_type = pr, 12 rows) and the four inconsistency classes
                            found in pandapower#3016 (record_type = inconsistency, 4 rows); columns below
lean/                       Lean 4 project "PhysAudit" (one module per PR) and build/axiom logs
replay/
  scripts/                  native replay scripts, one per PR (Python or Julia)
  env/                      environment setup and run scripts
  outputs/                  replay outputs (JSON; ClimaAtmos4822 as text)
controlled/
  tasks/<task>/             request.md (request text), snapshot/ (pre-change library subset), meta.json
                            (commits, invariant, sampling domain with sources, regression examples)
  tasks/LICENSES/           upstream license texts for the snapshots
  generate.py               sample generator (Codex CLI)
  samples/<task>/codex/     sNN.md (agent reply) and sNN.json (record: model, commands, events, usage, diff)
  judge/                    judge_b.py, test harnesses, reference patches (controls/), library download script
  results/                  full.jsonl (one record per task x variant) and full.stdout (one line per variant)
dq1/
  <code>/                   ΔQ1 spot-check scripts per PR (pp3016 = pandapower#3016; ca4637/ca4822 = ClimaAtmos.jl;
                            cp2851/cp3339 = CoolProp; ct2107 = Cantera; mo33307/mo34020 = MOOSE; oc6024/oc6083 =
                            Oceananigans.jl; pb5888/pb5889 = PyBaMM), each with impl_log.md (implementer's plan-to-code
                            log, original; English: impl_log.en.md); jl_common/ holds shared Julia code
  outputs/                  pandapower#3016 full runs (10^4 points each): pp3016_t1.json, pp3016_t2.json, pp3016_t3.json
                            (tracks 1, 2, 3; t3 from spotcheck_v2.py, the unguided track 3 plan) and pp3016_t3a.json
                            (guided track 3 plan, control only); pp3016_t1_zz_probe.txt (output of zz_bitwise_probe.py)
controlled/patch_similarity.py  similarity of generated patches to the historical reference patches (difflib)
third_party_licenses/       upstream license texts for code transcribed or re-implemented in lean/ and replay/
```

The five controlled tasks are `buildings_2725`, `buildings_3027`, `buildings_4775` (Modelica Buildings Library
issues #2725, #3027, #4775) and `msl_1335`, `msl_1530` (Modelica Standard Library issues #1335, #1530).
`buildings_4775` uses a reconstructed request text (see its `meta.json`).

### Columns of `audit/per_pr_judgments.csv`

| Column | Meaning |
|---|---|
| `spec_primary_invariant_agreement` | agree / partial: whether the two tracks' primary invariants agree substantively |
| `pr_tests_check_primary_invariant` | yes / partial / no, judged from reading the PR's tests under one uniform standard (tests not run); `pr_test_scope` gives the scope |
| `verdict_class` | discriminating_certificate (universal theorems that hold for the transcription and distinguish pre- and post-change code), partial_certificate_plus_violation, partial_certificate_plus_unresolved_convention, partial_certificate_uncovered_branch (main theorem premises exclude a branch the code executes), non_discriminating, low_content |
| `strict_criterion`, `lenient_criterion` | 1 if the PR counts under the strict / lenient certificate criterion (5 and 10 of 12) |
| `strict_boundary` | 1 for the two strict-criterion PRs that do not count under a stricter reading (change distinguished only by a constant, or no paired violation theorem and a failed transcription check) |
| `main_theorems`, `paired_violation_theorems` | Lean theorem names (namespace = module name) |
| `old_error_reproduced_on` | where the pre-change error stated by the violation theorem was reproduced: native, python_reimplementation, lean_column_model_float64 |
| `transcription_reference`, `transcription_check`, `transcription_metrics` | what the Lean transcription was compared with, pass/fail/partial, and the measured maximum relative differences (`n` = number of compared points) |
| `inconsistency_*`, `parent_commit_behavior`, `origin_class`, `adjudication` | for the pandapower#3016 rows: path, Lean theorems, native representative value, behaviour on the parent commit, origin (incomplete_new_feature, error_to_silent_wrong_value) and adjudication |

In `replay/outputs/pandapower3016_replay*.json` and `replay/scripts/pandapower3016_replay.py`, keys prefixed
`w1_` / `w2_` refer to issues raised by specification track 1 / track 2.

## Language of the files

Several files were produced as-is by LLM agents in Chinese and are kept unchanged as the original record:
the classification reasons and class labels in the inventory tables, `inventory/search_log.md`,
`inventory/search_scripts/decisions.txt`, the specification tracks (including the Chinese parts of the track 3
output) and the track 1-2 prompt, the ΔQ1 implementation logs, and the Chinese descriptive fields of
`controlled/tasks/*/meta.json`. English machine translations, produced with OpenAI Codex (model gpt-6-astra,
reasoning effort medium) and checked for row and heading alignment, are provided alongside; the Chinese
originals remain authoritative:

| Original | Translation |
|---|---|
| `inventory/pass1_classify.csv`, `inventory/pass2_classify.csv` | columns `invariant_class_en`, `reason_en` |
| `inventory/merged.csv` | columns `p1_class_en`, `p2_class_en`, `p1_reason_en`, `p2_reason_en` |
| `inventory/search_log.md` | `inventory/search_log.en.md` |
| `inventory/search_scripts/decisions.txt` | `inventory/search_scripts/decisions.en.txt` |
| `audit/specs/track1.md`, `audit/specs/track2.md` | `audit/specs/track1.en.md`, `audit/specs/track2.en.md` |
| `controlled/tasks/<task>/meta.json` | `controlled/tasks/<task>/meta.en.json` |
| `audit/specs/track12_prompt.txt` | `audit/specs/track12_prompt.en.txt` |
| `audit/specs/track3/<PR>/spec.md` | `audit/specs/track3/<PR>/spec.en.md` |
| `dq1/<code>/impl_log.md` | `dq1/<code>/impl_log.en.md` |

Comments and docstrings in the Lean modules and scripts were translated into English in place; only comments
and docstrings changed, which was checked by comparing each file with its original after removing comments (Lean,
Julia, shell) or docstrings (Python abstract syntax tree), so the Lean modules were not recompiled.
For readability, Chinese error messages, check labels and `help=` texts in string literals of the ΔQ1 scripts
for the 11 PRs whose outputs are not included, and two `help=` texts of `replay/scripts/pandapower3016_replay.py`,
were also translated; nothing else in those files changed. Strings that are written into archived outputs or
matched against the original records were left as they are: the category patterns in
`inventory/search_scripts/finalize.py` (they match the Chinese classification records), and the `install` and
`spec` strings of `dq1/pp3016/spotcheck*.py`, which appear in `dq1/outputs/pp3016_t*.json` (`spec` names the
track 3 output by its working path: `specs_codex_v2/` is `audit/specs/track3/` and `specs_codex/` is
`audit/specs/track3_guided/` in this archive; `env/pp3016` and `env/setup_pca_py.sh` correspond to
`replay/env/pp3016` and `replay/env/setup_python.sh`).
The raw Codex logs (`audit/specs/track3/*/codex.log`) are not translated. In them and in `input/pr.json`, the
temporary working directory appears as `<workdir>` and e-mail addresses other than no-reply trailer addresses are
replaced by `<email>`.

## Reproducing

Measured figures below come from the original runs. Nothing in this archive needs to run to read the data.

### Inventory

Requires an authenticated GitHub CLI (`gh`). From `inventory/search_scripts/`: `python3 search.py` (then
`search2.py`, `search3.py` for renamed repositories), `python3 build.py` (signature verification, mapping commit
hits to PRs, file lists; writes `../candidates.csv`), `python3 triage.py` (source-file marking), then the
per-PR judgments (`decisions.txt`) and `python3 finalize.py` (writes `../pass1_classify.csv`). The search ran
serially with a 2.2 s pause between calls and took about 1 h 45 min (about 841 search and 3800 core API calls).
Pass-1 patch reading (`decisions.txt`) and pass 2 were done by LLM agents (Claude Sonnet for pass 1) and cannot be
re-run deterministically. The search depends on the live GitHub search API, so a re-run will not return exactly
the same candidates (PRs merged, edited or deleted after 2026-10-16; repository renames).

The audit sample was drawn from the 42 PRs with `n_phys = 1` and `data_level = 0`, stratified by library
(2 PRs from libraries with at least 3 such PRs, 1 from each other library, random seed 20261017).

### Lean

Requires elan. The toolchain is pinned in `lean/lean-toolchain` (`leanprover/lean4:v4.34.1`) and Mathlib in
`lean/lake-manifest.json` (tag `v4.34.1`, commit `d13f23b723b8a846827a245b89c10fc7d3f11612`). From `lean/`:

```
lake exe cache get        # prebuilt Mathlib
lake build                # root module PhysAudit (imports all 12 modules)
bash scripts/build.sh     # per-module build logs and logs/summary.txt
bash scripts/axioms.sh    # #print axioms for every theorem; logs/axioms_<module>.txt
python3 scripts/summarize_logs.py   # rebuilds logs/summary.txt and logs/axioms_summary.txt from the logs
```

With the Mathlib cache, each module compiled in 5–17 s (see `logs/summary.txt`). The 12 modules have 7592 lines
and 415 theorems and lemmas; no `sorry`, `admit` or `native_decide`. `logs/summary.txt` and
`logs/axioms_summary.txt` were generated by `scripts/summarize_logs.py` from the per-module logs: build status
from `Build completed successfully` and `error:` lines, the module compile time reported by lake, and axiom sets
parsed from the `#print axioms` output (which Lean may wrap over several lines). propext, Classical.choice and
Quot.sound are Lean's standard axioms; every theorem uses only these. `logs/root_build.log` was recorded with one
additional smoke-test module that is not part of this archive and is imported by no other module; the build job count of a rebuild therefore differs from the one in that log.

### ΔQ1 spot checks

From the archive root, with the pandapower environment of the native replays (merge commit 0c6a1f0a):

```
replay/env/pp3016/bin/python dq1/pp3016/spotcheck.py --track 1 --n 10000 --out dq1/outputs/pp3016_t1.json
replay/env/pp3016/bin/python dq1/pp3016/spotcheck.py --track 2 --n 10000 --out dq1/outputs/pp3016_t2.json
replay/env/pp3016/bin/python dq1/pp3016/spotcheck_v2.py --n 10000 --out dq1/outputs/pp3016_t3.json
replay/env/pp3016/bin/python dq1/pp3016/spotcheck.py --track 3 --n 10000 --out dq1/outputs/pp3016_t3a.json
```

Points are split equally over the branches of each plan and drawn with seeds derived from
`20261017 + 10 * PR index + track`, point and branch number (`harness.py`), so a re-run reproduces the archived
rejection counts up to platform floating-point differences. In `pp3016_t3a.json`, a working-directory prefix inside one recorded exception message was
replaced by `<workdir>/`; no data values were changed. The four runs took 25-82 min each on one node.
`zz_bitwise_probe.py` re-evaluates the 2 track 1 points (3 checks) where the bitwise ZNyn/ZNd identity failed;
its output, `dq1/outputs/pp3016_t1_zz_probe.txt`, shows residuals of 2-4 x 10^-16 Ohm. Scripts for the other 11 PRs were written and smoke-tested (three also completed full runs for
tracks 1 and 2); their outputs are not included and were not analyzed.

### Native replays

| PR | Script | Environment (`replay/env/`) | Library version | Command (from `replay/`) |
|---|---|---|---|---|
| Cantera#2107 | `Cantera2107_replay.py` | `build_cantera.sh` -> `env/ctb` | Cantera built from merge commit d58a2326 (4.0.0a1), Python 3.12 | `env/ctb/bin/python scripts/Cantera2107_replay.py > outputs/Cantera2107_replay.json` |
| CoolProp#2851 | `CoolProp2851_replay.py` | `setup_python.sh` -> `env/py2` | CoolProp 8.0.0 (PyPI) | `env/py2/bin/python -I scripts/CoolProp2851_replay.py --n 2000 --seed 0 > outputs/CoolProp2851_replay.json` |
| CoolProp#3339 | `CoolProp3339_replay.py` | `env/py2` | CoolProp 8.0.0 (does not contain the PR; see script) | `env/py2/bin/python -I scripts/CoolProp3339_replay.py 2000 3339 > outputs/CoolProp3339_replay.json` |
| PyBaMM#5888, #5889 | `PyBaMM5888_replay.py`, `PyBaMM5889_replay.py` | `env/py2` | PyBaMM 26.10.1.0 | `env/py2/bin/python scripts/PyBaMM5888_replay.py > outputs/PyBaMM5888_replay.json` (same for 5889) |
| MOOSE#33307, #34020 | `Moose33307_replay.py`, `Moose34020_replay.py` | any Python >= 3.8 (+ mpmath optional) | no native MOOSE; Python re-implementation of the C++ | `python3 -I scripts/Moose33307_replay.py > outputs/Moose33307_replay.json` (same for 34020) |
| pandapower#3016 | `pandapower3016_replay.py`, `yzn_base_invariance.py` | `env/pp3016` (merge 0c6a1f0a), `env/pp_parent` (parent be6b67a8) | pandapower 3.4.0, Python 3.12 | `bash env/run_pandapower_replays.sh`; library tests: `bash env/run_pp_tests.sh` |
| ClimaAtmos#4822, #4637 | `climaatmos4822_diff.jl`, `ClimaAtmos4637_replay.jl` | `setup_julia.sh` -> `env/ca` | Julia 1.11.7; ClimaAtmos.jl at ba2bf122 (0.42.14; contains the #4637 merge a8efce5a) | `bash env/run_julia_replays.sh` |
| Oceananigans#6024, #6083 | `Oceananigans6024_replay.jl`, `Oceananigans6083_replay.jl` | `env/ocean` | Julia 1.11.7; Oceananigans.jl at 24de7026 (contains the #6083 merge 1734ff9e) | `bash env/run_julia_replays.sh` |

`setup_python.sh` needs `uv`, `build_cantera.sh` needs `micromamba` (set `UV` / `MICROMAMBA` if they are not on
PATH), `setup_julia.sh` downloads the official Julia tarball (set `JULIA_URL` for platforms other than macOS
aarch64). Each script states its inputs, sampling domain and seed in its header; defaults reproduce the archived
runs (2000 points per check unless stated). The Julia replays run single-threaded; ClimaAtmos4637 needed more than
8 GiB of memory. The pandapower replay with n = 40 parameter sets runs in seconds; run times of the other replays
were not recorded. Replay logs (`*.log`, `*.err`) are not included.

Two edits were made to the archived outputs: absolute paths of the original working directory in Python
tracebacks were replaced by `<env>/` (site-packages of the virtual environment) and `<workdir>/` (directory of the
replay scripts); and `pandapower3016_replay_parent.json` has an added `meta.env_commit` field. In that file
`meta.merge_sha` is a constant written by the script (the PR merge commit), while the run used the parent-commit
environment `env/pp_parent` (be6b67a8562a448d538f385a43c46718956cbc31), as also recorded in
`yzn_base_invariance_pp_parent.json` (`direct_url.vcs_info.commit_id`). No data values were changed.

### Controlled comparison

- **Libraries**: `controlled/judge/prep_libs.sh` downloads the full library trees at the pre-change commits and
  the MSL versions they depend on into `controlled/libs/` and prints the SHA-256 of each tarball.
- **Generation**: `GEN_CODEX_HOME=<isolated home> GEN_WORK_ROOT=<tmp dir> python3 controlled/generate.py --k 10`
  (Codex CLI, model `gpt-6-astra`, reasoning effort medium, read-only sandbox, no web search). Session
  identifiers were removed from the archived event records, and the temporary working directory appears as
  `<workdir>` in recorded commands. Generation is stochastic and cannot be re-run identically.
- **Judging**: OpenModelica 1.26.3 (container image not recorded) with MSL 4.1.0 provided by the container, and
  MSL 3.2.3 from GitHub for `buildings_2725`. Mount `controlled/` at `/work` and run `sh judge/run_full.sh`; it
  writes `results/full.jsonl` and `results/full.stdout`. Stage (ii) runs 10^4 points per variant (about 20 s per
  variant for the boiler harness, plus regression simulations with a 1800 s timeout each). Results are cached per
  deduplicated patch content (`key`), so identical patches share one result.

## Comparing a re-run with the archived data

- Lean: compare the regenerated `lean/logs/summary.txt` and `lean/logs/axioms_summary.txt` with the archived ones
  (build status, theorem counts and axiom sets must match; compile times will differ).
- Native replays: the scripts are seeded, so a re-run with the same library versions should reproduce the
  archived JSON values up to floating-point differences between platforms; compare the `max_rel*` fields and the
  pass/fail flags of each check. Differences in the Julia replays can also come from package versions resolved
  at install time.
- Controlled comparison: compare `valid`, `reg_pass`, `samp.n_failed` and the violation counts per task and
  variant in `results/full.jsonl`; the deduplication `key` identifies identical patched trees.

## Known gaps

- Only the top-level Python packages are pinned (CoolProp, PyBaMM, pandapower commits, Python 3.12); transitive
  dependency versions were not recorded.
- The Julia `Manifest.toml` files of the original environments were not archived.
- The OpenModelica container image was not recorded.
- The inventory depends on the live GitHub API; the first-batch search state (`state.json`), the API cache and the
  cached patches are not included.
- The code that drew the 12-PR sample is not included; the sampling rule and seed are stated above.
- LLM steps (pass-1 patch reading, pass 2, the three specification tracks, the ΔQ1 implementations, sample
  generation) are not deterministic.
- No native MOOSE build was used; the MOOSE replays compare a Python re-implementation of the C++ with the Lean
  formulas. CoolProp 8.0.0 does not contain PR #3339, so that replay is partial.

## License

- Code (scripts, Lean sources, judge, harnesses): MIT, see `LICENSE`. Exception: `replay/scripts/Moose33307_replay.py`,
  `replay/scripts/Moose34020_replay.py`, `dq1/mo33307/scm_port.py` and `dq1/mo34020/mortar_port.py` are LGPL-2.1.
- Data, specifications, judgments and documentation: CC BY 4.0, see `LICENSE-DATA`.
- Third-party material (library snapshots, transcribed or quoted upstream code, reference patches): upstream
  licenses, see `THIRD_PARTY.md`.
