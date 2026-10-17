# Comparison of the three specification tracks (unguided Track 3 prompt)

Comparison of the primary invariants and check plans of Track 3 (Codex, gpt-6-astra, run with an English translation of the original Track 1/2 prompt, without per-branch guidance) against Track 1 and Track 2, and against the guided Track 3 run ("v1a", see `audit/specs/compare3_guided.md`). The purpose is to judge how much of v1a's branch-coverage advantage came from the prompt. This run is called "v2" below.

Inputs: `audit/specs/track3/<PR>/spec.md`, `prompt.txt` and `codex.log`; `audit/specs/track1.md`, `track2.md` and `compare.md`; the criteria and v1a conclusions are taken from `audit/specs/compare3_guided.md`; PR-test judgments are those in `audit/per_pr_judgments.csv`. Everything is a comparison made by reading documents; no code was run.

## Criteria

- Same as in `compare3_guided.md`. **Primary-invariant agreement** uses the lenient criterion: if the core checkable statements are substantively the same, the entry is "same"; if some part is shared but the core of one track is missing from the other, "partial"; if nothing is shared, "different". A strict-criterion note is also given (counted as the same only if every primary invariant of the other track is also a primary invariant in Codex).
- **PR's own tests**: for each track, the first primary invariant is taken (T2 takes I1) and compared with the judgment harmonised in Table 3.
- **Expected verdict**: Section 7 of the Codex spec is a prediction made by reading the source, not a measurement.
- **Contamination**: for each PR, all commands and URLs in the `exec` blocks of `codex.log`, the `web search` entries, and the contents of the working directory were checked.
- **Difference from v1a**: only the amount of branch coverage is compared, recorded as "fewer / same / more".

## Prompt differences and two circumstances to disclose

1. The v2 prompt is an English translation of the original prompt of the first two tracks. It lacks the v1a Section 5 requirement to "write out which vector groups, phases, regimes and option values are covered", and it does not permit reading the base commit. The remaining fields of the two versions are identical, so the v1a-to-v2 difference can mainly be attributed to these two prompt differences (speculation; the run-to-run variance of the same model at the same effort was not measured).
2. The prompt is internally inconsistent: the Independence rule lists "PR/issue discussion" as readable, while the Rules say "no comments, issues". All 12 logs follow the stricter of the two and did not access discussion pages.
3. For 4822, 2851 and 5888, `input/` contains a `source/` directory whose file list is identical to the v1a run directory: the merge/base commit source and git tree (`tree.json`) that Codex downloaded during v1a. It contains no specification or analysis text. Only 4822 actually read it (the spec states "inspected through the supplied input/source copies"; fetching raw files directly did not succeed); the logs of 2851 and 5888 show no read of it. This is file residue across runs: the content is public source code and does not contaminate T1/T2, but the files that 4822 read were selected during v1a.

---

## Cantera#2107

- **Primary invariants**: v2 primary = I1 Mayer relation c_p−c_v=TVα²/κ_T (first compute the mixture V_T and V_P; when V_P=0 and V_T≠0 it returns c_p, marked as a library convention), I2 closure of the derivatives of the same Gibbs potential (v̄_k=∂_Pμ_k, Maxwell relation, the c_p° pressure-correction formula of SSVol, Margules analytic excess quantities, DH solvent derivatives). I3 activity-basis conversion is demoted to secondary.
  - vs T1: **same** (T1's A and C correspond to I1 and I2). Strict criterion: partial (T1's primary B, chemical potential–activity, is the secondary I3 in v2).
  - vs T2: **same** (T2's I1 and I2 are both in v2's I2, and are compared at 1 atm and at several non-reference pressures). Strict criterion: same.
- **Branch coverage**: compressible VPSS (HKFT dilute solution), true constant-volume control, the branch V_P=0 and V_T≠0 (it states explicitly that returning c_p is a convention, cannot be used as evidence for a physical equality, and that the denominator must not be replaced by a small positive number), Margules at elevated pressure, SSVol (P≠P_r), both DH variants (IAPWS and constant A), two dilutions of IdealMolal. The cutoff region is explicitly excluded. T1 tests only compressible phases and T2 does only algebraic comparisons: **Codex more**.
- **PR's own tests**: C = Partial (I1; I2 Partial, the core identity of I3 yes); T1 = Partial; T2 = No; Table 3 = Partial. C agrees with Table 3. It also notes that the newly added heat-capacity difference assertion "shares the formula with production code", the same point as T1's "nearly tautological".
- **Expected violations**: none. The physical interpretation of the V_P=0 branch is marked uncertain. It still points out that `PDSS_HKFT::cv_mole`, mentioned in the PR text, does not exist in the merged source (same as v1a). The non-primary I4 of v1a (pre-existing defects such as `getdlnActCoeffds` apparently being an N-fold derivative) does not appear in v2.
- **Contamination**: none seen. Only files of the merge commit on raw.githubusercontent were fetched.
- **Difference from v1a**: branch coverage **fewer** (the 5 PDSS types, the 10 DH fixtures and the cutoff model are not listed one by one), but the branch κ_T=0 with β≠0 is still singled out.

## ClimaAtmos#4637

- **Primary invariants**: v2 primary = I1 u³=0 at top/bottom, plus the volume sum of the separate mass advection being zero (requires the zero flux to hold structurally in the divergence operator), I2 constant tracer A(ρχ)=CAρ, and zero SGS correction for a uniform scalar. I3 S_ρ=S(ρq_t) is secondary, marked as a library convention.
  - vs T1: **same**. Strict criterion: same (T1's A and B correspond to I1 and I2; in v1a mass conservation was non-primary, now it is merged into I1).
  - vs T2: **same** (lenient criterion). Strict criterion: partial (T2's I2 requires χ≡1 across time stepping; v2 states explicitly that it goes only as far as instantaneous tendencies and does not claim multi-step preservation for the full moist system).
- **Branch coverage**: Schar/Linear baseline, 20/40/80 levels, Float64/Float32, central/first_order/third_order/vanleer reconstructions of SGS, DYCOMS 1M single updraft, a check with artificially contaminated end-face caches for structural conservation, several values of C. No SLEVE, no 3-D warped box, no combinations of several microphysics and turbulence options. T1/T2 are mainly Schar Linear/SLEVE and a single EDMF example: **Codex slightly more** (more reconstructions and precisions, no SLEVE).
- **PR's own tests**: C = Partial; T1 = Partial; T2 = Partial; Table 3 = Partial. All three agree.
- **Expected violations**: none. It states explicitly that the hand-written Jacobian is not required to equal the exact derivative (the PR omits environment terms on purpose). For the sponge layer it only records that the PR added a TODO and did not fix it, and excludes it from the checks.
- **Contamination**: none seen.
- **Difference from v1a**: **fewer** (v1a had combinations of 3 geometry classes × 4 turbulence classes × 5 microphysics classes); v1a's non-primary I6 (Jacobian SGS velocity block) no longer appears.

## ClimaAtmos#4822

- **Primary invariants**: v2 primary = I1 moisture contrast along isentropes (V_q,h=0 when q=Aθ+B, invariant under horizontal rotation), **I3 continuity back to the purely vertical closure as c_x→0⁺**. I2 (f∈(0,1], s²ℓ²+fΔz²=Δz², covariance PSD, total-moisture budget H+Gfv²=c_gΔz²v² when b=0) is secondary.
  - vs T1: **same** (lenient: T1's A corresponds to I1). Strict criterion: partial (T1's B and C are both the secondary I2 in v2).
  - vs T2: **same** (lenient: T2's I2 corresponds to I1). Strict criterion: partial (T2's I1, "no double counting", is written into the budget formula of v2's secondary I2).
  - It still does **not** write T2's variance consistency q′q′=F′²θ′θ′ (the Lean theorem `isentropic_variance_inconsistency`); it only says that the closure's partitioning of θ and q is a modelling choice.
- **Branch coverage**: θ_z above and below the lower bound, c_x, c_z ∈ {0,1} and a c_x→0 scan, both signs of q=Aθ+B, b=0, v=0, rotated gradients, humidity limiter tested separately, Float32/64; the turbulence block is switched off for isolation. T1/T2 are mainly the lower-bound branch and function-level sampling: **Codex more**.
- **PR's own tests**: C = Partial (I1; I3 No); T1 = Partial; T2 = No; Table 3 = Partial. **C agrees with Table 3** (v1a recorded No, which disagrees with Table 3).
- **Expected violations**: **I3 (primary) may not satisfy**: f does not depend on a nonzero c_x, while at c_x=0 the code sets f to 1, so even an arbitrarily weak horizontal sampling removes a finite amount of vertical variance; it gives the analytic prediction V_q(ε)=V_q(0)(1+ε²)/2, C_θq(ε)=C_θq(0)/2 (worked example with f=1/2). In v1a this was the non-primary I5, and neither of the other two tracks wrote it.
- **Contamination**: no network overreach seen; it read the cross-run residue `input/source` (see item 3 at the beginning).
- **Difference from v1a**: **fewer** (v1a also covered two correlation modes, a ρ_turb scan, several microphysics and turbconv options, and sgs_distribution), but it promotes the c_x→0 discontinuity to primary.

## Oceananigans#6024

- **Primary invariants**: v2 primary = I1 bottom-control neutral equilibrium layer κ_eff=B(a³/D)^{1/4}≈0.4 and ν=κu★d (log law). The distance law and proportions, and TKE dissipation consistency, are secondary.
  - vs T1: **same**. Strict criterion: same.
  - vs T2: **same**. Strict criterion: same.
- **Branch coverage**: both sides of bottom control and surface control, both sides of the min-branch transition, old and new B and random B, three diffusivities, refined Nz, Float32/64, the Δt→0 dissipation limit; Ri branches other than N²=0, immersed boundaries and closure tuples are explicitly excluded. T1/T2 use only a single-column steady state: **Codex slightly more**.
- **PR's own tests**: C = No; T1 = No; T2 = No; Table 3 = No. All three agree.
- **Expected violations**: none. The old B=0.28 should clearly fail κ≈0.4. It did not read the full text of Wagner et al. 2025 (states explicitly that it was not obtained), so the typographical inconsistency in Eqs. (35)/(53) pointed out by v1a is absent from v2.
- **Contamination**: none seen. api.github.com was used only for the git tree of the merge SHA.
- **Difference from v1a**: **fewer** (v1a had four length paths, both sides of each Ri threshold, immersed boundaries, closure tuples, two time discretisations).

## Oceananigans#6083

- **Primary invariants**: v2 primary = I1 the total drift velocity is conserved in the explicit/implicit split (w^e+w^i=w) and the result is independent of the source decomposition (including exact cancellation, with the centroid displacement as the convergence benchmark). I2 conservation of the tracer total in a closed domain and I3 GM combined slope truncation and non-divergent eddy velocity are secondary (v1a had promoted GM to primary).
  - vs T1: **same**. Strict criterion: partial (T1's primary B, centroid and mass, is the secondary I2 in v2 for mass conservation).
  - vs T2: **same**. Strict criterion: partial (T2's primary I2, conservation, as above).
- **Branch coverage**: QAB2/SplitRK3 × Courant {0, 0.1, 0.49, 0.51, 1, 4, 8} × C_*∈{0.3, 0.5} × sources (biogeochemical drift, AdvectiveForcing, split by λ=0.25/0.5/1.5, exact cancellation, GM frozen-field replay) × uniform/non-uniform columns × Float32/64 × optional GPU; on the GM side m∈{0, B/2, 2B}, slope ratio r∈{0,…,2}, rotation angles. T1: three supply modes × 2 steppers; T2: two modes plus a Courant scan. **Codex more**: both sides of the threshold and GM as a velocity source are present.
- **PR's own tests**: C = Partial; T1 = Partial; T2 = No; Table 3 = Partial. C agrees with Table 3.
- **Expected violations**: none. A new documentation issue is found: the docstring of `isopycnal_rotation_tensor_components.jl` writes stable stratification as b_z<0, which is inconsistent with the positive-denominator branch of the implementation. The GM change is listed in Section 1 and its test judged No, but, as in v1a, it does not state explicitly that the PR description does not mention it.
- **Contamination**: none seen.
- **Difference from v1a**: **fewer**, by a small margin (only the WENO5 scheme is used, no immersed grid; both sides of the threshold are still present).

## CoolProp#2851

- **Primary invariants**: v2 primary = I1 bubble-point/dew-point (Q=0, Q=1) tangent E_q=dp_q/dT and A_q+B_qD_q=0. The chain-rule derivatives of the endpoint properties (I2) and pure-fluid Clapeyron (I3) are secondary.
  - vs T1: **same**. Strict criterion: partial (T1's primary B is the binary ε→0 convergence to the pure fluid; v2 mentions the single-component reduction in one sentence under I3 and does no ε→0 check).
  - vs T2: **same**. Strict criterion: partial (T2's primary I2, the pure-component limit, as above).
- **Branch coverage**: Q=0 and Q=1, composition a∈[0.6,0.8], temperature scan, residuals of A and B independently reconstructed from fugacity coefficients, derivatives of h/s/ρ endpoint properties, both ends of pure n-propane, exception for intermediate Q. No PQ parametrisation, no pseudo-pure fluid, no reversed component order, no ternary, no second derivatives. T1/T2: Q=0, Q=1 plus the pure-fluid limit. **Codex slightly more** (lacks ε→0).
- **PR's own tests**: C = Partial; T1 = Partial; T2 = Partial; Table 3 = Partial. All three agree.
- **Expected violations**: none. v1a's non-primary I6 (the mixture `second_saturation_deriv` still uses the pure-fluid formula) does not appear in v2.
- **Contamination**: none seen. curl results were written to `/tmp/cp2851-backend.cpp`, a temporary write allowed by the sandbox, not a read of local files.
- **Difference from v1a**: **fewer**.

## CoolProp#3339

- **Primary invariants**: v2 primary = I1 the logarithmic blend prescribed by the configured crossover centre (the geometric mean of the two branches when a=x). No overshoot and conditional positivity (I2), and the directional derivative of moving the centre (I3), are secondary. It does **not** write "the 8 shipped fluids are bitwise unchanged before and after the merge"; it only keeps the PR test's a=2 and 12 as comparison points.
  - vs T1: **partial**: T1's primary A is exactly the bitwise invariance for 8 fluids, which is missing in v2; T1's primary B (parameter semantics) corresponds to I1.
  - vs T2: **partial**: T2's primary I1 (bitwise invariance for 8 fluids) is missing; T2's primary I2 (interpolation bounds) corresponds to v2's secondary I2.
- **Branch coverage**: only one fluid, R1233zd(E); a=x+d scanned at several levels over d∈[−0.4,0.4], a=2/12, single-phase screening, a C=0 control clone to isolate η_0. T1/T2 both require 8 fluids. **Codex fewer** (1 fluid versus 8, with a centre scan added).
- **PR's own tests**: C = Partial; T1 = Partial; T2 = Partial; Table 3 = Partial. All three agree.
- **Expected violations**: none. It only points out that the dilute-limit deviation grows when a moves near zero, a model-applicability limit (same as v1a, and in the same direction as the Lean theorem `dilute_limit_broken_at_xc_zero`).
- **Contamination**: none seen. api.github.com was used only for the git tree of the merge SHA.
- **Difference from v1a**: **fewer** (v1a had 8 fluids, liquid/vapour/supercritical grids, a=0); v1a's primary invariant I2, "a=2 backward compatibility", was also dropped.

## pandapower#3016

- **Primary invariants**: v2 primary = I1 the neutral 3Z_N sits at the correct winding port, is not scaled by K_T, and does not change the positive/negative sequence (including, for YNyn, the T-equivalent port formula Z_{0,H}=αZ_k+3Z_N+Z_m‖[(1−α)Z_k+Z_{g0,H}], which requires 3Z_N to be added outside the HV branch and not first added to Z_k and then split); I2 zig-zag transformer zero-sequence decoupling (zero cross admittance, z_H=αz_k+3Z_N/Z_B). The IEC single-phase current formula and NER current limiting (I3) are secondary.
  - vs T1: **same** (core ΔZ₀=3Z_N, not multiplied by K_T). Strict criterion: partial (T1's primary B, the IEC current formula, is the secondary I3 in v2).
  - vs T2: **same**. Strict criterion: partial (T2's primary I2, IEC current and NER asymptotics, is the secondary I3 in v2).
- **Vector groups and modes covered by the check plan**:
  - Mode: only short circuit `calc_sc` (fault=1ph; 3ph and 2ph are also used to check that the positive sequence is unchanged), case=max, two values of `inverse_y`, calculation base S_B∈{1,10,100}. `pf_3ph`/`runpp_3ph` are **not covered**; the operating domain states explicitly that it does not cover "... unbalanced three-phase power flow".
  - Vector groups: Dyn (observed on LV), YNd (observed on HV), **YNyn** (grid only on LV, fault on the HV side, scanning finite Z_m and α, using the T-circuit formula as an independent absolute-value oracle), ZNyn, ZNd, ZNy (three-bus fixture; the YNd main transformer blocks source-side zero sequence, and a grid is added on the aux bus to change the zero-sequence strength to test decoupling).
  - Explicitly excluded: tap, `parallel>1` (it does not presume whether the NER is shared or per unit), PSU, pfe/i0. **Yzn is not in the plan**; **Zd is only listed in Section 1** as a name the code accepts (znyn/znd/zny/zd), with no checks; ungrounded groups, pi/T model comparison and case=min are all absent.
  - The four key points: YNyn **covered**; Yzn **not covered**; pf_3ph **not covered**; Zd **not covered**.
- **Branch coverage compared with the two tracks**: T1 plans Dyn, YNd, ZNyn/ZNd, plus a suggestion for Dyn under pf_3ph (and mentions that the semantics of YNyn with both sides grounded is unclear); T2 plans Dyn, YNd, ZNd/ZNyn, Zd/Zy. v2 adds quantitative checks for YNyn and ZNy, and lacks T1's pf_3ph and T2's Zd. **Roughly equal, each misses some**.
- **PR's own tests**: C = Partial (I1; I2 and I3 Partial); T1 = Partial; T2 = Partial; Table 3 = Partial. All three agree. The description of the tests (only Dyn tested, parallel=1, 1 MVA base, ΔR₀=0.03, ΔX₀=0.06; the zig-zag transformer only gets qualitative comparisons such as current sign and a ratio of 0.5–2) is the same as in Table 3.
- **Problems raised by the expected verdicts** (predictions from reading the source, not measured):
  - **YNyn: I1 may not satisfy**: the code adds 3Z_N to the total leakage impedance and then splits it by (α, 1−α); with a finite magnetising branch the HV neutral is not fully in series with the HV branch. This agrees with the mechanism of Lean theorem A7. It also states that the API cannot specify the two neutral points independently, and the specification is written "following the HV-side selection semantics of the merged source", so the convention dependence is identified as well.
  - Zig-zag I2: uncertain, because the parameter mapping of αz_k lacks confirmation from literature or equipment data; this is not a violation verdict.
  - Name `ZN`: the PR description says `ZN`, but the code tuple has only znyn/znd/zny/zd, no zn (T2 and v1a also point this out).
  - Yzn, pf_3ph and Zd: no prediction at all. Of the other 4 items of v1a (PSU with zig-zag counting the neutral twice, tap_corr using the wrong side, the semantics of parallel, ZN being rejected), only ZN remains.
  - Of the 4 inconsistency classes of pandapower#3016 (2 violations, 2 convention-dependent), v2 identifies 1 class (YNyn).
- **Contamination**: none seen. Only raw files of the merge commit were fetched; there was only one web search, which opened `pd2ppc_zero.py` at the merge SHA. Unlike v1a it did not search for the PowerFactory technical reference; the spec states that IEC 60909, PowerFactory Fig. 8.3 and J&P Appendix 5 were not obtained, and marks them unverified.
- **Difference from v1a**: **fewer**, and the most pronounced reduction of the 12: v1a covered the 10 vector groups of the pandapower#3016 native replay and pf_3ph, while v2 retains only 6 vector groups and a single mode, sc.

## MOOSE#34020

- **Primary invariants**: v2 primary = I1 constant reproduction of the projected velocity s=W/D, invariance under scaling by a positive factor, invariance under process partitioning (the per-process ratios must not be averaged), and independence of normalize_c. Time lag and isotropy (I2), and the Coulomb bound with the correct history coefficients and non-positive friction power (I3), are secondary.
  - vs T1: **same** (T1's primary A corresponds to I1). Strict criterion: partial (T1's primary B, equivalence of function friction and an equal constant μ, appears only as a paired computation in I3's check plan).
  - vs T2: **same** (T2's primary I1 corresponds to I1). Strict criterion: partial (T2's primary I2, exact projection of the dual basis, is not written as an invariant; v2 only says that dual test functions may be negative and that nodewise equality is not required for non-uniform fields).
- **Branch coverage**: 2D/3D, both values of normalize_c, MPI 1/2/4, mesh area/cell count/quadrature order, algebraic checks of the ONE/TWO residuals, paired constant and function friction, `--test-restep`, newly contacting nodes, velocity sign reversal and tangential basis rotation, rejection of mu≤0. Wear, edge dropping and recover are excluded. T1: 2D plane plus dual basis; T2: 3 meshes × 2 quadratures plus dual basis. **Codex more** (lacks the dual-basis linear field).
- **PR's own tests**: C = Partial (only an Exodiff golden standard, `rel_err=5e-4`); T1 = No; T2 = Partial; Table 3 = No. The original wordings disagree; rereading the v2 text under Table 3's rule that "only a golden regression counts as No", it should be No (this is a reclassification, not Codex's own wording).
- **Expected violations**: the primary is expected to be satisfied. I3's stronger assertion that "any function dissipates" **may not satisfy**: this PR does not clamp the output of `function_friction`; the F=0.3+0.5·2.17^{−p/100}−10s used in the PR test is positive only when s is below a threshold (same as v1a). New interface nodes and the parallel lifecycle are marked as to be measured.
- **Contamination**: none seen.
- **Difference from v1a**: **fewer**, by a small margin (lacks correct_edge_dropping, wear and the element-LM branch).

## MOOSE#33307

- **Primary invariants**: v2 primary = I1 the Cheng–Todreas friction uses the same bulk Re_b to select the flow regime, R_L/R_T endpoint continuity, positive friction and non-negative dissipation, with an explicit coefficient oracle for the EDGE subchannel, which states "W_s>0 for 8≤H/D≤52, hence B≥1"; I2 the Borishanskii log10 formula, the Pe≤200 plateau and continuity at Pe=200. The denominator truncation of the modified Gnielinski (I3) is secondary.
  - vs T1: **same** (T1's primary A, Ws_T>0, primary B, continuity and shared R_b, and primary C, Borishanskii, are all in v2's primary). Strict criterion: same, but W_s>0 is a one-sentence corollary within I1, and the check plan only does formula-fidelity comparison without a separate sign check.
  - vs T2: **same** (T2's primary I1, Ws_T>0, and primary I2, transition-region continuity, are both in I1). Strict criterion: same (same note).
  - This is a reversal relative to v1a: v1a happened to lack Ws_T>0 (the main Lean theorem `WsT_pos`), whereas v2 writes it into the primary.
- **Branch coverage**: EDGE subchannel coefficient oracle, triangular lattice, wire-wrapped and bare rods (I2 uses bare rods), Re_b at t∈{0, 0.01, …, 1} and outside the interval, local and bulk regimes deliberately different, mixing/sweep for CENTER-adjacent and peripheral gaps, Pe=1/50/199/200/201/1000/2200, geometric scaling s, both sides of Gnielinski d=0.1, 7–217 rods, zero-flow diagnostics. Only a single rank; no square lattice, no CENTER/CORNER coefficients, no multiple solvers. T1/T2: H/D scan, both sides of R_L/R_T, three solvers, multiple processes. **Each has strengths and weaknesses**: v2 has more closed branches, while the two tracks have the solver and parallel dimensions.
- **PR's own tests**: C = Partial (all are CSVDiff golden standards and RunApp); T1 = No; T2 = Partial; Table 3 = No. The original wordings disagree; rereading v2's text under the harmonised rule it should be No (reclassification).
- **Expected violations**: I1 is expected to be satisfied in the normal positive-flow domain; the **zero-flow boundary may not satisfy**: Σm_iμ_i/M in `computeBulkReynoldsNumber()` is 0/0 at all-zero flow (same as v1a; neither of the other two tracks wrote it). It does **not** point out what both other tracks wrote, that the bulk Re is computed only on rank 0: v2 only says to start with a single rank and places MPI state outside the real-number statement; T1's "the constructor reads the default value 1" is also absent.
- **Contamination**: none seen. One api.github.com contents query (the test directory under the merge SHA) returned 404.
- **Difference from v1a**: **fewer** (v1a covered CENTER/EDGE/CORNER × triangular/square lattices, ring counts 2/3/9, single-block and multi-block implicit, MPI 2 ranks).

## PyBaMM#5888

- **Primary invariants**: v2 primary = I1 SI conversion from cm²/s to m²/s and the absolute scale (D(1000, 298.15)≈3.2227e-10). Positivity and monotonicity (I2) and agreement with Xu2019 (I3) are secondary.
  - vs T1: **same**. Strict criterion: partial (T1's primary, agreement with Xu2019, is the secondary I3 in v2).
  - vs T2: **same**. Strict criterion: partial (T2's primary I2, Xu2019 plus monotonicity, is secondary in v2).
- **Branch coverage**: c∈[500,1500], the four corners of T∈[278.15,318.15], the default point, the three PR points and 100 random points, central differences in the temperature and concentration directions. T1/T2 are a grid plus monotonicity (T2's domain is wider: c up to 3000, T 263–333 K). **Roughly equal**.
- **PR's own tests**: C = Partial; T1 = Yes; T2 = Yes; Table 3 = Partial. C agrees with Table 3 and differs from the original wording of the two tracks. v2 additionally notes that `assert_allclose(..., rtol=1e-12)` sets no atol, and NumPy's default atol=1e-8 is larger than the diffusivity under test, so this assertion is not in fact a 1e-12 relative test (a zero value would also pass). Neither of the other two tracks nor v1a wrote this; it can serve as a supplementary reason for the Partial judgment in Table 3.
- **Expected violations**: none.
- **Contamination**: none seen. Unlike v1a it did not access the Valøen–Reimers original or the NumPy documentation; the spec states that the paper was not read.
- **Difference from v1a**: **fewer** (v1a also had the c_min=10 lower bound on three sides, scalar/vector inputs, conductivity).

## PyBaMM#5889

- **Primary invariants**: v2 primary = I1 lookup by the solid-phase concentration θC, nodes reproduce table values, and invariance under reparametrisation of C at fixed c. Independence of the initial electrolyte concentration e (I2) and in-table positivity and envelope (I3) are secondary.
  - vs T1: **same**. Strict criterion: same (both of T1's primary invariants are in I1).
  - vs T2: **same** (T2's primary I1 corresponds to I1). Strict criterion: partial (T2's primary I2, DFN low-rate capacity invariance, is not written; v2 states explicitly that voltage or capacity invariance is not required).
- **Branch coverage**: 34 nodes, the midpoints of 33 intervals with several λ, C varied between 1 and 1.2 times, e=500/1000/1500 and random values, flux sign; extrapolation outside the table, temperature, vector input and model-side clamping are explicitly not done. T1: grid plus c_max variation; T2: nodes plus DFN. **Codex slightly more**.
- **PR's own tests**: C = Partial; T1 = Partial; T2 = Partial; Table 3 = Partial. All three agree.
- **Expected violations**: none. T2 pointed out that D_e.csv is low by an order of magnitude; v2 only says that the correction of the electrolyte expected value does not correspond to a new physical implementation, and does not comment on the data itself.
- **Contamination**: none seen.
- **Difference from v1a**: **fewer** (v1a also had extrapolation at both ends, upper and lower clamping by the model adapter, temperature, vector input).

---

## Summary table

| PR | Primary vs T1 | vs T2 | Branch coverage (v2 vs the two tracks) | PR's own tests C / T1 / T2 (Table 3) | Expected violations (v2, prediction) | Contamination | Difference from v1a (branch coverage) |
|---|---|---|---|---|---|---|---|
| Cantera#2107 | same (strict: partial) | same (strict: same) | Codex more; singles out the V_P=0, V_T≠0 branch | P / P / N (P) | none | none seen | fewer |
| ClimaAtmos#4637 | same (strict: same) | same (strict: partial) | Codex slightly more (lacks SLEVE) | P / P / P (P) | none | none seen | fewer |
| ClimaAtmos#4822 | same (strict: partial) | same (strict: partial) | Codex more; still does not write variance consistency | P / P / N (P) | **primary I3: c_x→0⁺ discontinuity** | no network overreach seen; read the source copy left over from v1a | fewer |
| Oceananigans#6024 | same | same | Codex slightly more | N / N / N (N) | none | none seen | fewer |
| Oceananigans#6083 | same (strict: partial) | same (strict: partial) | Codex more | P / P / N (P) | none (documentation sign inconsistency) | none seen | fewer (slightly) |
| CoolProp#2851 | same (strict: partial) | same (strict: partial) | Codex slightly more (lacks ε→0) | P / P / P (P) | none | none seen | fewer |
| CoolProp#3339 | **partial** (lacks 8-fluid bitwise invariance) | **partial** (same) | Codex fewer (1 fluid versus 8) | P / P / P (P) | none | none seen | fewer |
| pandapower#3016 | same (strict: partial) | same (strict: partial) | roughly equal, each misses some: has YNyn, ZNy; lacks Yzn, pf_3ph, Zd | P / P / P (P) | **primary I1: YNyn split**; ZN name mismatch | none seen | fewer (most pronounced) |
| MOOSE#34020 | same (strict: partial) | same (strict: partial) | Codex more (lacks dual-basis linear field) | P / N / P (N) | primary none; an arbitrary friction function can be negative (conditional) | none seen | fewer (slightly) |
| MOOSE#33307 | same (strict: same) | same (strict: same) | each has strengths and weaknesses (no multiple solvers, no MPI) | P / N / P (N) | primary none in the normal domain; zero-flow 0/0 (outside the domain); rank-0 issue not pointed out | none seen | fewer |
| PyBaMM#5888 | same (strict: partial) | same (strict: partial) | roughly equal | P / Y / Y (P) | none | none seen | fewer |
| PyBaMM#5889 | same (strict: same) | same (strict: partial) | Codex slightly more | P / P / P (P) | none | none seen | fewer |

## Counts

- **Primary invariants (lenient criterion)**: vs T1 same 11, partial 1 (CoolProp#3339), different 0; vs T2 same 11, partial 1 (the same PR), different 0. The counts are the same as in v1a, but the "partial" PR changes from MOOSE#33307 (v1a lacked Ws_T>0) to CoolProp#3339 (v2 lacks 8-fluid bitwise invariance).
- **Primary invariants (strict criterion)**: vs T1 same 4 (4637, 6024, 33307, 5889), partial 8; vs T2 same 3 (Cantera, 6024, 33307), partial 9; different 0 in both.
- **Branch coverage, v2 vs the two tracks**: Codex more or slightly more in 8 (Cantera, 4637, 4822, 6024, 6083, 2851, 34020, 5889), roughly equal or each with strengths and weaknesses in 3 (pandapower, 33307, 5888), fewer in 1 (3339). v1a was "most" in 11 and "roughly equal" in 1.
- **Branch coverage, v2 vs v1a**: fewer 12, same 0, more 0 (6083 and 34020 are only slightly fewer).
- **PR's own tests**: v2 gives Yes 0 / Partial 11 / No 1. It agrees with the judgment harmonised in Table 3 in 10/12; the disagreements are the two MOOSE PRs (Codex records a golden regression as partial); rereading the two descriptions under the Table 3 rule gives 12/12 (reclassification, not Codex's own wording). v1a's original wording agreed in 9/12; the disagreeing ClimaAtmos#4822 now agrees with Table 3 in v2. The PRs for which the three tracks' original wordings on the first primary invariant are fully identical remain 6 (4637, 6024, 2851, 3339, pandapower, 5889).
- **Expected violations**: at the primary level, 2 PRs are expected to fail: pandapower#3016 (YNyn) and ClimaAtmos#4822 (c_x→0 discontinuity, non-primary in v1a). Two more PRs are expected to fail only on conditional or out-of-domain assertions (34020 arbitrary friction function, 33307 zero flow). Three non-primary expected failures of v1a (Cantera I4, 4637 I6, 2851 I6) do not appear in v2.
- **pandapower**: the check plan covers 6 vector groups (Dyn, YNd, YNyn, ZNyn, ZNd, ZNy) and only one mode, sc. Of the four key points, YNyn is covered and predicted (the mechanism agrees with Lean theorem A7, and the HV-side semantic dependence is flagged); Yzn, pf_3ph and Zd are neither covered nor predicted. v1a covered 4/4 and identified 4/4.
- **Contamination**: none seen in the 12. All execs were inside the working directory; network access was only files of the merge commit on raw.githubusercontent and the git tree/contents of the merge SHA on api.github.com (6024, 6083, 3339, 33307); the web search of each log has only one group of entries, and all opened raw files of the merge commit (33307 has none). No PR comments, issues, pulls/issues endpoints or any literature page were accessed. v1a had several web-search entries with empty content that cannot be fully excluded; v2 has no such gap. The only item to disclose is the v1a-downloaded source left in 3 directories, of which 4822 read it.

## Notes on interpretation

1. With the per-branch guidance removed, the branch coverage of all 12 is lower than in v1a, and pandapower drops from 4/4 to 1/4. This supports the judgment at the start of v1a: a considerable part of v1a's advantage in branch coverage, especially on pandapower, came from the prompt (the evidence is an A/B comparison with the same model at the same effort; a single run, run-to-run variance not measured, so still speculation).
2. The agreement of the primary invariants is insensitive to the prompt: both v1a and v2 give 11/12 same under the lenient criterion, and the only difference is which PR is "partial".
3. Even without guidance, v2 still has more branch coverage than the two tracks on 8/12 PRs, and pandapower still independently wrote the quantitative YNyn check and the violation mechanism; so changing the model family itself also brings part of the coverage difference (speculation; no same-prompt repeated runs were made on the Claude side).
4. All expected verdicts come from reading the source, with no compilation or execution; the agreement with the native replay measurements is an after-the-fact comparison, not an independent reproduction.
