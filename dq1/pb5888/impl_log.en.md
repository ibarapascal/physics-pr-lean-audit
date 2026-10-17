# PyBaMM#5888 ΔQ1 spot-check implementation log

> Purpose: record the sentence-by-sentence mapping of the three tracks' check plans to `spotcheck.py`, and all implementation judgments, for review under §6/§7 of the ΔQ1 protocol.
> Structure: §0 environment → §1 Track 1 → §2 Track 2 → §3 Track 3 → §4 time and memory → §5 deviation candidates → §6 post-smoke corrections → §7 to be ruled on (classification of smoke rejections).
> Update rule: after the smoke run only script errors may be fixed; changing the plan mapping must be recorded in protocol §7 and a line appended to §5 of this file.
> References: upstream: the ΔQ1 protocol (preregistration), the #5888 sections of `audit/specs/track1.md`/`track2.md`, `audit/specs/track3_guided/pybamm-team_PyBaMM_5888/spec.md` §2/4/5; downstream: `dq1/outputs/pb5888_t{1,2,3}.json`.

## §0 Environment and general conventions

- Node: compute node A, `env/py2` (built by `setup_pca_py.sh`, pybamm 26.10.1.0 from PyPI, not built from the merge commit). Verified: `Ai2020.py`, `Xu2019.py`, `lithium_ion_parameters.py`, `settings.py` are byte-for-byte identical to the GitHub raw files of merge commit 343ece2 (`diff -q` prints nothing). Implementation judgment: treated as the native environment of the merge commit.
- Seed 20261017 + 11×10 + track number (t1=20261127, t2=20261128, t3=20261129). Each track has only 1 branch (the changed function has no option/phase/regime branches), and all N points belong to it.
- The grids, discrete values and nominal values given by the plan are added as mandatory points (not counted in N), with results counted separately.
- Verdict, replay and JSON format: see the header of `spotlib.py`: any check that is finite and exceeds its tolerance → rejection candidate → counted as a rejection only if recomputation with the same input is bitwise identical; an exception or NaN counts as a failure.
- Value API: `ParameterValues(...).evaluate(pv[key](pybamm.Scalar(c), pybamm.Scalar(T)))` (the same for all three tracks). To prevent unbounded growth of the internal expression cache over 10⁴ points, `_processor.clear_cache()` is called every 2000 evaluations (a private attribute, guarded by hasattr; does not affect the numerics).
- Note on the implementer's reading scope: the Track 3 spec.md was once output as a whole file (including §1, §6–§9); the implementation relies only on §2, §4, §5.

## §1 Track 1 (`audit/specs/track1.md`, "numerical check method")

| Plan text | Code |
|---|---|
| `pv = ParameterValues("Ai2020")`; `f = pv["Electrolyte diffusivity [m2.s-1]"]`; `pv.evaluate(f(Scalar(c), Scalar(T)))` | `D_ai()` |
| Usual operating domain c_e∈[100,3000], T∈[253,333] | `t1_sample`: both variables linearly uniform |
| Grid over (c_e,T) (c 100…3000 in steps of 100, T 253…333 in steps of 5) | mandatory points `grid_c=…_T=…`, 30×17=510 points, running all the same checks as the sampled points |
| Compare with numpy directly from the literature formula (including 1e-4), rtol 1e-12 | `t1_api_vs_literature_rtol1e-12` (`t1_lit`, atol=0) |
| Compare with the same-named parameter of `ParameterValues("Xu2019")`, rtol 1e-12 | `t1_ai_vs_xu2019_rtol1e-12` |
| Magnitude bound 1e-11 ≤ D ≤ 2e-9 (invariant "secondary") | `t1_magnitude_1e-11<=D<=2e-9` |
| Finite-difference sign ΔT=0.01 K, Δc=0.1 mol/m³, central difference, sign only; ∂D/∂T>0, ∂D/∂c<0 | `t1_fd_sign_dD_dT>0`, `t1_fd_sign_dD_dc<0` (differences of API values) |
| Optional: 1C DFN discharge (Ai2020), final-state c_e gradient max−min ≫ 1 mol/m³ | mandatory point `optional_dfn_1C`: `lithium_ion.DFN()` + Ai2020 default parameters (Current function 2.28 A = 1C), `sim.solve([0,3600])`, take the full-cell `Electrolyte concentration [mol.m-3]` at the last time |

Implementation judgments:
1. The "optional" check is included per protocol §2 ("everything the plan writes counts"); it has no sampleable input and is run once.
2. "≫ 1 mol/m³" takes the threshold ≥10 mol/m³ (one order of magnitude).
3. Grid points and sampled points are not deduplicated.

## §2 Track 2 (`audit/specs/track2.md`, "numerical check method")

| Plan text | Code |
|---|---|
| Usual operating domain c_e 0–3000, T 263–333, DFN discharge C/10–2C | `t2_sample`: c_e, T, C_rate each linearly uniform |
| I1: compare with the value computed independently by numpy from the original formula, rtol 1e-12 (the original formula in mol/L: c=c_e/1000) | `t2_I1_api_vs_formula_rtol1e-12` (`t2_formula`) |
| I1: magnitude bound 1e-11 ≤ D_e ≤ 1e-9 | `t2_I1_magnitude_1e-11<=D<=1e-9` |
| I1: agreement with the PR's 3.2227e-10 at 1000 mol/m³, 298.15 K | mandatory point `nominal_1000_298.15`: abs tolerance 0.00005e-10 (half of the last place of 5 significant digits) |
| I2: compare with the same-named parameter of Xu2019, rtol 1e-12 | `t2_I2_ai_vs_xu2019_rtol1e-12` |
| I2: monotonicity central difference ΔT=0.01, Δc=0.1 "check the sign at 100 random points" | `t2_I2_fd_sign_*`, done at all N sampled points |
| I3: `DFN()` (Ai2020), 1C discharge for 3600 s; take c_e and temperature on the separator grid; second-order difference for ∂c/∂x; compare with (1−t⁺)i/(F·D_e,eff), tolerance 5% | `DFNCase` + `t2_dfn_checks`: `np.gradient(c, x_s, edge_order=2)`; D_e,eff = D_e(c,T)·ε_sep^b (`Separator porosity`=0.5, `Separator Bruggeman coefficient (electrolyte)`=1.5); t⁺, F, i taken from the parameter set / `pybamm.constants.F` / `Total current density [A.m-2]`; rtol 0.05 at each separator node |
| I3: quasi-steady (after t ≳ 1000 s) | if the end time is <1000 s, this item is recorded as "not applicable" (recorded as a pass) |
| I3: relative drift of ∫εc_e dx from the initial value ≤1e-6 | `t2_I3_electrolyte_Li_drift<=1e-6`: library variable `Total lithium in electrolyte [mol]` final/initial |
| The discrete value 1C | mandatory point `C_rate=1_c=1000_T=298.15` |

Implementation judgments:
1. The "100 random points" are replaced by the uniform sampling of N points of protocol §3 (the same check is done at all N points).
2. The C rate is drawn uniformly from the operating domain [0.1, 2], with 1C as a mandatory point; the discharge duration is 3600/C (extending "1C discharge for 3600 s" to the nominal full-discharge time), terminating automatically at the 3.0 V cutoff under Ai2020; values are taken at the last time of the solution.
3. To reuse the already discretised model, `Current function [A]` is set to `[input]` and `sim.solve([0, 3600/C], inputs=...)` is called; current = C × `Nominal cell capacity [A.h]` (2.28).
4. The D_e in D_e,eff uses the model-side `LithiumIonParameters.D_e` (with a 10 mol/m³ floor; c_e in the separator is far above the floor, so it has no effect); the temperature is `Separator temperature [K]`, and the ambient temperature if that is not available.
5. Sign convention: discharge current is positive, and −D_eff ∂c/∂x and (1−t⁺)i/F are compared with sign.
6. Replay: when the rejected check involves only the constitutive term, the replay recomputes the constitutive term and reuses the DFN result just computed for the same point; when `t2_I3_*` is involved, the whole point (including DFN) is recomputed (`t2_replay`).
7. When c_e is close to 0 the central difference with Δc=0.1 uses −0.1 mol/m³; the original formula is still defined there and is computed anyway.

## §3 Track 3 (spec.md §2/§4/§5)

| Plan text | Code |
|---|---|
| Working window 500≤c≤2000, 273.15≤T≤333.15 | `t3_sample`, linearly uniform |
| Cartesian 7×5 grid | mandatory points `grid_c=…_T=…` (35 points, running all checks of `t3_eval`) |
| 200 uniform pairs from `default_rng(5888)` | replaced by the protocol's N points and protocol seed |
| I1: compare with C=c/1000, T_g=229+5C, `np.exp(np.log(10)*(…))*1e-4`, rtol 1e-12, atol 0, explicitly finite | `t3_I1_api_vs_oracle_rtol1e-12` (`t3_oracle_D`; `SL.close` records non-finite values as failures) |
| I1: nominal value 3.222722528605257e-10 with the same tolerance | mandatory point `nominal` |
| I1: if a baseline is installed separately, compare new/base=1e-4 | no baseline installed (condition not met), not done |
| I2: Ai2020 and Xu2019, rtol 1e-12 per sample | `t3_I2_ai_vs_xu2019_rtol1e-12` |
| I2: model-side `lp.D_e` against the original function after applying the concentration floor | `t3_I2_model_De_vs_raw_floor_rtol1e-12` (`LithiumIonParameters().D_e`) |
| I2: c=0,5,9.999,10,10.001,20 × each listed temperature; read the setting and require =10 | mandatory points `boundary_c=…_T=…` (30 points), `t3_I2_setting_D_e__c_e==10` |
| I3: raw and model-side finite and >0 (all normal points and boundary points) | `t3_I3_raw_finite_positive`, `t3_I3_model_finite_positive` |
| I4: each grid row strictly increasing in T, each column strictly decreasing in c | mandatory point `grid_monotone` |
| I4: at random interior points c∈[501,1999], T∈[273.16,333.14], central differences of ln D with h_c=0.1, h_T=0.001, compared with the analytic log derivative at rtol 1e-6 | done when the sampled point falls in this interior zone: `t3_I4_dlnD_dT/dc_vs_analytic_rtol1e-6` |
| I4: the wrapper function is constant below the floor (rtol 1e-12), continuous across the floor (covered by the formula comparison), strictly decreasing above the floor; no derivative at the kink | boundary points `t3_I4_wrapper_constant_below_floor`, `t3_I4_wrapper_strict_decrease_above_floor` (compared with the previous boundary concentration) |
| I5: conductivity compared with the §1 polynomial at rtol 1e-12; finite and non-negative; raw c=0 exactly 0; nominal 1.1943 rel 1e-3; `kappa_e` compared with its floor | `t3_I5_*` (`t3_oracle_kappa`) |
| Paired Vectors on the grid, scalar/Vector combinations compared elementwise with scalars at rtol 1e-12 | mandatory point `vector_grid` (both the diffusivity and conductivity functions) |

Implementation judgments:
1. "Transcription of §1" for the conductivity polynomial: the implementer did not rely on §1, and wrote the reference formula directly from the source of `electrolyte_conductivity_Ai2020` in the merge commit and κ=1e-4·c·P² of §2 I5.
2. All sampled points are above the floor, and the floor branch is covered only by the mandatory boundary points (plan §5 says this is the only relevant branch).
3. "Continuous across the floor" has no tolerance of its own and is represented by the rtol 1e-12 comparison "model side = original function(max(c,10))" at each boundary point.

## §4 Expected time and peak memory (speculation, not run; to be updated after the smoke run)

| Track | Per point | 10⁴ points | Basis |
|---|---|---|---|
| t1 | ~15 ms (7 evaluate calls) | ~3 min + 510 grid points + 1 DFN (~10 s) | one evaluate takes about 1–3 ms (expression processing + evaluation) |
| t2 | ~0.5–1 s (1 DFN solve + constitutive) | ~1.5–3 h (about 30% of the points are rejected by the magnitude bound and replay the constitutive term, reusing the DFN) | the Ai2020 DFN default-mesh IDAKLU solve takes about 0.3–0.8 s, and the model is discretised only once |
| t3 | ~25 ms | ~5 min | about 10 evaluate calls |

Peak memory: about 1–1.5 GB (pybamm/casadi import about 0.4 GB, the discretised DFN model and solution about 0.5 GB); the JobSpec declares 3 GiB.
If the smoke run shows t2 above 1.4 s per point (10⁴ points exceeding 4 h), N of t2 is reduced per protocol §3 to the largest value that can be completed in 4 h and recorded in §7.

## §5 Deviation candidates

- (None) The environment is the PyPI release rather than a build from the merge commit: the relevant source files are byte-for-byte identical, recorded as a note and not as a deviation.

## §6 Post-smoke corrections (2026-10-17)

- None. The rejections of the first round all come from the magnitude bounds of the plan itself (§7); the script is faithful to the plan and was not changed.

## §7 To be ruled on (classification of smoke rejections, script not changed)

1. **(P) The Track 1 "secondary" invariant magnitude bound 1e-11 ≤ D ≤ 2e-9 ("within the domain") does not hold on the plan's own operating domain c_e∈[100,3000], T∈[253,333]**. In the first round 77 rejections (sampled 14/100, grid mandatory 63/510) are all the lower bound of `t1_magnitude`; at the same points API = literature formula (rtol 1e-12, 0 rejections), Ai2020 = Xu2019 (0 rejections), derivative signs all pass: the library value is exactly the formula given in the plan. The rejected points concentrate in the low-temperature high-concentration corner: T ≲ 272 K and c_e ≳ 440 mol/m³, and the closer to the singularity T = 229 + 0.005c_e, the smaller (e.g. c_e=2855, T=255.2 → 2.5e-14). By hand calculation from the plan's original formula: c_e=1000, T=253 → 3.2e-12; c_e=3000, T=253 → 8e-16, both below 1e-11. With uniform sampling over the rectangular operating domain, the area fraction violating the lower bound is about 10.4% (Monte Carlo with 10⁶ points, using the original formula), so about 1000 rejections of the same kind are expected in the full run; the upper bound 2e-9 is not triggered in this domain (maximum about 1.06e-9). Where the plan is wrong: "within the domain 1e-11 ≤ D ≤ 2e-9" contradicts the Valøen–Reimers formula given in the same section and the operating-domain note "must avoid the singularity" (from first principles: the formula gives D→0 near the singularity).
2. **(P) The Track 2 I1 magnitude bound 1e-11 ≤ D_e ≤ 1e-9 ("within the domain above", c_e∈[0,3000], T∈[263,333]) likewise does not hold**. In the first round 3 rejections (sampled) are all the lower bound: c_e≈2820–2955, T≈266–271 → 3.8e-12 to 8.3e-12; the I1 formula comparison, the Xu2019 comparison, derivative signs, DFN separator flux and lithium conservation all pass. By the original formula, the region violating the lower bound is T ≲ 272 K and c_e ≳ 1830 (fraction about 2.5%), and **the upper bound 1e-9 is also violated at c_e ≲ 200, T ≳ 324 K** (e.g. c_e=0, T=333 → 1.12e-9, fraction about 0.5%, not drawn in the first 100 points and will appear in the full run). Where the plan is wrong: I1 "within the domain above 1e-11 ≤ D_e ≤ 1e-9" is inconsistent with the original formula in the same sentence at the corners of the domain.
3. No (S) or (V) class was seen for this PR in the first round.
