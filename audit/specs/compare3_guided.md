# Comparison of the three specification tracks (guided Track 3 prompt, not counted in ΔQ1)

Comparison of the primary invariants and check plans of Track 3 (Codex, gpt-6-astra, run with the guided Track 3 prompt) against Track 1 and Track 2. This guided run is called "v1a" in `audit/specs/compare3.md`, which repeats the comparison with the unguided prompt.

Inputs: `audit/specs/track3_guided/<PR>/spec.md` (Track 3, Codex, with the common prompt `prompt.txt`); `audit/specs/track1.en.md`, `track2.en.md` and `compare.md`; `audit/per_pr_judgments.csv`; `audit/per_pr_judgments.csv` (PR-test scope and the pandapower#3016 inconsistency rows). Everything is a comparison made by reading documents; no code was run.

## Criteria

- **Primary-invariant agreement**: the same lenient criterion as the two-track comparison in `compare.md`. If the core checkable statements of both sides are substantively the same, the entry is "same"; if some part is shared but the core of one side is missing from the other or only secondary, "partial overlap"; if nothing is shared, "different". A **strict-criterion** note is also given: counted as the same only if every primary invariant of the other track is also a primary invariant in Codex (under the strict criterion T1 vs T2 itself would also disagree on several PRs, so the figures cannot be compared directly with the 11/12 Track 1–2 agreement in `compare.md`).
- **PR's own tests**: for each track the first primary invariant is taken (T2 takes I1, consistent with the raw track counts T1 Y1/P8/N3, T2 Y1/P7/N4); for Codex I1 is taken, with the other primary invariants noted in parentheses.
- **Expected verdict**: Section 7 of the Codex spec is a prediction made by reading the source, not a measurement.
- **Contamination**: for each PR, all `exec` commands in `codex.log`, the URLs in python snippets, and the `web search` entries were checked.

## A confounding factor that must be disclosed

The Track 3 prompt differs from that of the first two tracks: its Section 5 explicitly requires "Be explicit about which branches, cases, or configurations (for example every vector group, phase, regime, or option value the changed code distinguishes) the plan covers", its Section 2 requires "If the correct behaviour depends on a modelling convention ... say so", and it permits reading the base commit source. The field descriptions of T1/T2 contain no per-branch enumeration requirement. The example "vector group" points directly at pandapower. Therefore the advantage of Codex in branch coverage and convention annotation below cannot be attributed to the change of model family alone (speculation: the prompt plays a comparable role).

---

## Cantera#2107

- **Primary invariants**: Codex primary = I1 potential-function/heat consistency (including v̄_k=∂_Pμ_k and the c_p° pressure-correction formula of SSVol), I2 Mayer relation (including the κ_T=0 convention), I3 activity-basis consistency.
  - vs T1: **same**: I1/I2/I3 correspond one-to-one to T1's C/A/B. Strict criterion: same.
  - vs T2: **same**: both of T2's primary invariants (∂μ_k/∂P=v̄_k, ∂c_p/∂P=−T∂²v/∂T²) are written in Codex's I1, and the check plan also tests the pressure term at P_*, 2P_* and 20P_*, the same point as T2's "excess-volume term is zero at 1 atm". Codex's primary set is equivalent to the union of the two tracks, which were originally only partially consistent with each other. Strict criterion: same.
- **Branch coverage**: Codex covers three cases, κ_T≠0, rigid (v_T=v_P=0), and **v_P=0 but v_T≠0**, the last explicitly marked "a declared library convention, not a consequence of the identity or a continuous limit". This is the same point as the boundary of the Lean certificate (when κ_T=0 and β≠0 no differentiable isochoric path exists). It also covers all 5 PDSS types, Margules with and without excess volume at elevated pressure, 10 DH fixtures, and the cutoff model. T1 only says that for κ_T=0 it takes c_v=c_p by convention, and its plan tests only compressible phases; T2 does only algebraic comparisons. **Codex most**.
- **PR's own tests**: C = I1 "core identity yes / full domain partial", recorded as Partial (I2 Partial, I3 yes at a given state point); T1 = Partial; T2 = No; Table 3 = Partial. Codex also points out that the new Mayer test "can share the implementation's errors", consistent with T1's "nearly tautological".
- **Expected violations**: all primary invariants are expected to be satisfied. The non-primary I4 is expected to fail: `getdlnActCoeffds` is suspected to be an N-fold derivative when there is more than one interaction term, and `dlnX_diag` has a suspicious product. Both are defects that already existed before the PR, and Codex distinguishes them itself. It also points out that the PR does not actually implement `PDSS_HKFT::cv_mole`.
- **Contamination**: none seen. Only source of the merge commit on raw.githubusercontent was fetched.

## ClimaAtmos#4637

- **Primary invariants**: Codex primary = I1 u³=0 at top/bottom, I2 uniform tracer A_χ=CA_ρ, I3 zero SGS correction for a uniform scalar, I4 S_ρ=S_qtot (marked as a library convention). Column mass conservation is the non-primary I5.
  - vs T1: **same**: u³=0 and χ≡1 agree; T1 merges mass conservation into its primary A, while Codex demotes it to I5. Strict criterion: partial.
  - vs T2: **same** (lenient criterion): both share impermeability and χ≡1. But T2's main axis is column mass conservation after time stepping and multi-step χ≡1, while Codex explicitly does not claim multi-step preservation for the whole moist system (I2/I3 go only as far as instantaneous tendencies). Strict criterion: partial.
- **Branch coverage**: geometry (flat column/plane, Schar Linear and SLEVE, 3-D warped box) × turbulence (none, EDOnly, prognostic EDMFX SGS on/off, 1/2 updrafts) × microphysics (dry/0M/1M/2M/2MP3) × all upwinding options. T1/T2 are mainly Schar Linear/SLEVE and a single EDMF example. **Codex most**.
- **PR's own tests**: C = Partial (I2, I3 Partial, I4 Yes); T1 = Partial; T2 = Partial; Table 3 = Partial. All three agree.
- **Expected violations**: all primary invariants are expected to be satisfied. Non-primary I6: the SGS velocity block of the hand-written Jacobian is an "interpolation of a product" rather than a "product of interpolations", and is expected to be unequal for general profiles. It is the same block as T2's "Jacobian drops O(aʲ²) terms", but the mechanism differs. It also marks "SGS water flux counted into ρ" as a convention and raises the passive-tracer problem (the same doubt as T1's), and also mentions the sponge.
- **Contamination**: none seen. It read the merge commit source and listed files with the api.github.com git tree.

## ClimaAtmos#4822

- **Primary invariants**: Codex primary = I1 gradient contrast along isentropes (G=0 when q=αθ+β), I3 covariance realizability (PSD, |C|≤√(V_θV_q)). The "no double counting" partition identity and (sℓ)²+fΔz²=Δz² are demoted to the secondary I4/I2.
  - vs T1: **same**: T1's A and C correspond to Codex's I1 and I3; T1's B is secondary in Codex. Strict criterion: partial.
  - vs T2: **same** (lenient criterion): T2's I2 corresponds to Codex's I1; T2's main axis I1, "no double counting", is the secondary I4 in Codex. Strict criterion: partial.
  - Codex does **not** write T2's variance consistency q′q′=F′²θ′θ′ (the Lean theorem `isentropic_variance_inconsistency`, |ρ|=√f); it only says "the asymmetric treatment in which the θ variance is not multiplied by f is a modelling choice".
- **Branch coverage**: θ_z above/below the lower bound, the four zero/nonzero combinations of c_Δx and c_Δz, two correlation modes, a ρ_turb scan, humidity cap, several microphysics and turbconv options, sgs_distribution. T1/T2 are mainly the lower-bound branch and function-level sampling. **Codex most**.
- **PR's own tests**: C = No (I1 and I3 both No); T1 = Partial (but its text says A/B/C all have no tests); T2 = No; Table 3 = Partial.
- **Expected violations**: the primary invariants are expected to be satisfied. The non-primary I5: the limit c_Δx→0⁺ is discontinuous with c_Δx=0 (a finite jump (1−f)vb²), which is expected to fail if understood with continuous-amplitude semantics; Codex marks it as awaiting a modelling ruling. This is a new finding that neither of the other two tracks wrote.
- **Contamination**: none seen.

## Oceananigans#6024

- **Primary invariants**: Codex I1 = κ=C(a³/b)^{1/4}≈0.4 plus the neutral wall law K_u=κu★d.
  - vs T1: **same** (both primary invariants are in I1). Strict criterion: same.
  - vs T2: **same**. Strict criterion: same.
- **Branch coverage**: four length paths for centre/face, both sides of each Ri threshold, bottom/surface/stratification limits, immersed-boundary grids, closure tuples, two time discretisations. T1/T2 use only a single-column steady state. **Codex most** (the PR itself only changes a constant, so the practical significance of the extra coverage is limited).
- **PR's own tests**: C = No; T1 = No; T2 = No; Table 3 = No. All three agree.
- **Expected violations**: none. The old value gives κ_b of 0.1159, the same as the 0.116 of T1 and T2. It also points out that Eqs. (35)/(53) of Wagner et al. 2025 are printed as a product, inconsistent with the division in Eq. (32), judged a typographical problem in the paper, editing status unverified. T1 pointed out something different: the old value 0.28 corresponds to swapped coefficient positions.
- **Contamination**: none seen. It read the Wagner paper cited by the PR, which is within what is permitted.

## Oceananigans#6083

- **Primary invariants**: Codex primary = I1 consistency of the total velocity, with the result unchanged when an equal total velocity is redistributed among sources (including the centroid displacement), I2 conservation of the tracer total, **I3 GM combined slope reduction**.
  - vs T1: **same** (T1's A and B correspond to Codex's I1 and I2). Strict criterion: same.
  - vs T2: **same** (T2's I1 and I2 correspond to Codex's I1 and I2). Strict criterion: same.
  - Codex promotes the GM slope reduction to primary. Both other tracks treat it as secondary and point out "not mentioned in the PR description, not tested". Codex lists this change in full in Section 1 and judges the test of I3 = No, but does not state explicitly that the PR description does not mention it.
- **Branch coverage**: AB2/SplitRK3 × CFL ratio {0,0.1,0.99,1,1.01,2,8} × velocity sources (biogeochemical drift, AdvectiveForcing, transport velocity of a two-dimensional circulation, GM closure, combinations, mutual cancellation) × three schemes × stretched/immersed/variable grids × Float32/64 × CPU/GPU. T1: three supply modes × 2 steppers; T2: two modes plus a Courant scan. **Codex most**, especially covering both sides of the threshold and the GM closure as a velocity source.
- **PR's own tests**: C = Partial (I2 Partial, I3 No); T1 = Partial; T2 = No (T2's I2 is Yes); Table 3 = Partial.
- **Expected violations**: none. It explains that a positive `minimum_bz` means taking the lower bound first, which is a convention.
- **Contamination**: none seen.

## CoolProp#2851

- **Primary invariants**: Codex primary = I1 bubble-line/dew-line tangent (Q=0, Q=1), I2 pure-fluid Clapeyron consistency plus pseudo-pure-fluid compatibility.
  - vs T1: **same** (core I1 agrees). Strict criterion: partial: T1's B is binary-mixture ε→0 convergence to the pure fluid, whereas Codex's I2 checks only the pure fluid itself and does not write this limit.
  - vs T2: **same**. Strict criterion: partial (same reason).
- **Branch coverage**: Q=0 and Q=1, both QT and PQ flash parametrisations, independent weighted residuals based on fugacity coefficients, both ends of the pure fluid, the R410A pseudo-pure fluid, reversed binary order and all ternary permutations, rejection of intermediate Q, second derivatives. T1/T2: Q=0, Q=1 plus the pure-fluid limit. **Codex most**.
- **PR's own tests**: C = Partial; T1 = Partial; T2 = Partial; Table 3 = Partial. All three agree.
- **Expected violations**: the primary invariants are expected to be satisfied. Non-primary I6: `second_saturation_deriv` for mixtures is still built as TΔv/Δh, and is expected not to be the derivative of the new slope (failure or unsupported). Neither of the other two tracks wrote this.
- **Contamination**: none seen.

## CoolProp#3339

- **Primary invariants**: Codex primary = I1 the smooth logarithmic blend prescribed by the configured crossover control, I2 backward compatibility at a=2 (bitwise identical).
  - vs T1: **same** (one-to-one with T1's two primary invariants). Strict criterion: same.
  - vs T2: **same** (T2's I1 corresponds to Codex's I2; T2's I2, "locality and interpolation bounds", is covered by Codex's I1 formula plus the I3 bound). Strict criterion: same.
- **Branch coverage**: 8 fluids, liquid/vapour/supercritical grids, a seven-level a=x+d scan for each state, a=2/12/0, dilute limit. T1: bitwise comparison of 8 fluids plus both sides of x_c; T2 is similar. **Roughly equal**, with Codex adding the a=0 dilute limit.
- **PR's own tests**: C = Partial; T1 = Partial; T2 = Partial; Table 3 = Partial. All three agree.
- **Expected violations**: none. It only points out that at a=0 the dilute limit is not normalised, a pre-existing modelling limitation exposed once the change makes the parameter effective.
- **Contamination**: none seen.

## pandapower#3016

- **Primary invariants**: Codex primary = I1 neutral-point law and correct placement (exactly one 3Z_N in the zero-sequence arm of each grounded winding, added after the K_T correction; for YNyn a T-equivalent matrix is given), I2 zig-zag transformer decoupling and the leakage-reactance split convention, I4 dimensional and base covariance (including runpp_3ph). The IEC single-phase short-circuit current formula is demoted to the secondary I3.
  - vs T1: **same**: the core ΔZ₀=3Z_N (not multiplied by K_T) agrees. T1's B (IEC current formula) is secondary in Codex. Strict criterion: partial.
  - vs T2: **same**: T2's I1 (ΔZ₀=3Z_N, Möbius form) corresponds to Codex's I1; T2's I2 (IEC current and NER asymptotics) is the secondary I3 in Codex. Strict criterion: partial.
- **Branch coverage**: Codex's plan covers Dyn, YNd, Yyn, YNy, **YNyn (driving-point impedance at both ends compared with an independently assembled T matrix, scanning s to expose the wrong split s·3Z_N/(1−s)·3Z_N)**, **Yzn**, ZNyn/ZNd/ZNy/Zd, ZN (expected to be rejected), four ungrounded groups (Yy/Yd/Dy/Dd, where Z_N should have no effect), PSU (YNd PSU, PSU with zig-zag), pi/T models, **Dyn/YNyn/Yzn under pf_3ph**, case max/min, two solvers, several MVA bases, taps, number of parallel units and in-service status. All 10 vector groups tested in the pandapower#3016 native replay are included. T1's plan has only Dyn, YNd, ZNyn/ZNd (plus a suggestion of Dyn under pf_3ph); T2 has Dyn, YNd, ZNd/ZNyn, Zd/Zy. Neither track has YNyn or Yzn. **Codex most, and it fills in the two vector groups that both Claude tracks missed.**
- **PR's own tests**: C = Partial (I2 Partial, I4 No); T1 = Partial; T2 = Partial; Table 3 = Partial. All three agree. Codex's description of the tests (only Dyn tested, case=max, R 0.01/X 0.02; the zig-zag transformer only gets qualitative comparisons) is the same as in Table 3.
- **Expected violations** (predictions from reading the source, not measured):
  - **YNyn**: expected to fail under the HV-side neutral convention: 3Z_N is first added to the total leakage reactance and then split, so the HV arm gets only s·3Z_N and the LV arm gets an extra (1−s)·3Z_N. It also points out that the docstring does not specify the neutral-point assignment of YNyn and that the code treats it as HV (convention dependence). This agrees with the mechanism of Lean theorem A7 and with its "convention-dependent" ruling.
  - **Yzn**: placement "suspect", and "Yzn retains an extra baseMVA multiplier ... suggesting a base-dependence failure". This corresponds to the two sources identified in the audit (the wrong (1−si) split and the 1/sn base dependence); the base dependence is stated clearly, while the wrong split goes only as far as "suspect".
  - **pf_3ph (T model)**: "Likely failure": the leakage reactance uses 3S_b/U² while the neutral uses S_b/U², and the placement in the T model is a separate problem. This agrees with the mechanism of Lean theorem B1.
  - **Zd**: "Uncertain physical interpretation": the absence of N in the name does not show that the neutral is accessible, and an ungrounded zig-zag transformer should not acquire a path to ground from the winding connection alone. It is only marked as a convention to be settled, and no violation is judged. This agrees with the audit's ruling "convention-dependent, pending upstream confirmation".
  - Four further predictions that neither track wrote: with PSU plus zig-zag the neutral may be counted twice (HV side (s+1)·3Z_N); when the HV/LV rated-voltage ratio is inconsistent, the HV zig-zag shunt uses the LV-side tap_corr; under `parallel` the neutral is treated as "shared", conflicting with the "one per unit" interpretation; the `ZN` that the PR claims to support is actually rejected (T2 also noticed that there is no "zn" in the code tuple).
- **Contamination**: none seen. It read only the merge and base commit source. The web search was a single "DIgSILENT technical reference ... zig zag zero sequence", the PowerFactory technical reference cited by the PR, which is within what is permitted.

## MOOSE#34020

- **Primary invariants**: Codex primary = I1 consistent mortar projection and constant reproduction (including exact projection of the dual basis), I2 constitutive parameters taken from the last accepted step.
  - vs T1: **same** (T1's A corresponds to Codex's I1). Strict criterion: partial (T1's B, equivalence of function friction and an equal constant μ, is only a comparison within Codex's I5 plan).
  - vs T2: **same** (T2's I1 and I2 are both in Codex's I1). Strict criterion: same.
- **Branch coverage**: 2D/3D × normalize_c × projection order ONE/TWO × correct_edge_dropping × MPI 1/2/4 × constant/function friction law, retry/restep, newly contacting nodes, wear, element-LM branch. T1: 2D plane plus dual basis; T2: 3 meshes × 2 quadratures plus dual basis. **Codex most**.
- **PR's own tests**: C = Partial; T1 = No; T2 = Partial; Table 3 = No. Codex records the golden-regression-only case as partial (same as T2); under the harmonised rule of Table 3 (only a golden regression counts as No), Codex's textual description should be classified as No.
- **Expected violations**: the primary invariants are expected to be satisfied. Non-primary I5: the friction law F=0.3+0.5·2.17^{−p/100}−10s bundled with the PR becomes negative once the velocity exceeds a threshold, and a user function is not checked either, so the unconditional claim "the new code always gives a positive friction coefficient/dissipation" is expected to fail (conditional). Parallel ownership migration is marked uncertain, consistent with T2.
- **Contamination**: none seen.

## MOOSE#33307

- **Primary invariants**: Codex primary = I1 correlation and flow-regime fidelity (all channels/gaps share R_b, continuous at R_L and R_T), I2 Borishanskii normalisation (log10), I3 the regularised Gnielinski is finite and non-negative, I4 cross-flow correlation matrix and mass conservation (index fix #33340).
  - vs T1: **partial overlap**: shares B (flow-regime continuity, shared R_b) and C (Borishanskii). T1's main axis A (Ws_T=19−11·log10(H/D)>0, i.e. "wire wrap only increases resistance") is **not written by Codex as an invariant**; I1 only requires "equal to the corrected wire-wrap factor in Section 1", with formula-fidelity comparisons at H/D=8/20/52. Strict criterion: partial.
  - vs T2: **partial overlap**: shares T2's I2 (transition-region continuity); T2's I1 (Ws_T>0) is missing. Strict criterion: partial.
  - The main Lean theorem for this PR is exactly `WsT_pos`.
- **Branch coverage**: CENTER/EDGE/CORNER × triangular/square × bare rod/wire-wrapped × both sides of coefficient table 1.1 × R_b and R_i deliberately in different flow regimes × two blending functions × all gap types, ring counts 2/3/9, both sides of Pe=200, both sides of the Gnielinski lower bound 0.1, single-block implicit plus multi-block, MPI 2 ranks, zero flow. T1/T2: H/D scan, both sides of R_L/R_T, three solvers, multiple processes. **Codex most**, but it lacks the sign check of Ws_T.
- **PR's own tests**: C = Partial (I1–I4 all Partial); T1 = No; T2 = Partial; Table 3 = No. Codex likewise records the golden CSV regression as partial; under the harmonised rule it should be No.
- **Expected violations**: the primary invariants are expected to be satisfied, with the note "MPI/lifecycle coverage uncertain": the plan requires comparing serial and 2 ranks because the bulk Re is computed only on rank 0 (an issue both other tracks flagged; Codex lists it only as to be investigated). It points out that the constructor of the triangular friction closure reads the default value 1 before the solve computes R_b (same as T1). Non-primary I6: when the total inlet flow is 0, `computeBulkReynoldsNumber` gives NaN from 0/0, expected to fail (neither track wrote this). The non-negativity of the blended Gnielinski in the transition region is uncertain.
- **Contamination**: none seen.

## PyBaMM#5888

- **Primary invariants**: Codex primary is only I1 (SI conversion from cm²/s to m²/s, giving D(1000,298.15)=3.2227e-10). Agreement with Xu2019 and monotonicity are the secondary I2/I4.
  - vs T1: **same** (the core unit conversion agrees). Strict criterion: partial (T1's second primary invariant, agreement with Xu2019, is secondary in Codex).
  - vs T2: **same**. Strict criterion: partial (T2's I2 is secondary in Codex).
- **Branch coverage**: full Cartesian grid of c×T plus random points, three sides of the model-side c_min=10 lower bound, scalar/vector input, conductivity. T1/T2: grid plus monotonicity. **Codex slightly more** (adds the lower-bound branch).
- **PR's own tests**: C = **Partial**; T1 = Yes; T2 = Yes; Table 3 = Partial. Codex gave partial directly, without a harmonised-rule prompt, the same as the harmonised judgment of Table 3 and different from the original judgment of the two tracks.
- **Expected violations**: none.
- **Contamination**: none seen. The web search accessed Valøen–Reimers 2005 (the literature cited by the PR, via ResearchGate) and the NumPy documentation (general API documentation, not an analysis of this PR).

## PyBaMM#5889

- **Primary invariants**: Codex primary = I1 lookup using the solid-phase concentration coordinate, nodes reproduce table values; I2 independence of the initial electrolyte concentration E, and invariance under reparametrisation of c_max at fixed c.
  - vs T1: **same** (one-to-one with T1's two primary invariants). Strict criterion: same.
  - vs T2: **same** (T2's I1 corresponds to Codex's I1). Strict criterion: partial (T2's I2, DFN low-rate capacity invariance, is not written by Codex, which states explicitly that DFN is not needed).
- **Branch coverage**: 34 nodes, 5 quantile points in each of the 33 intervals, extrapolation at both ends, upper and lower clamping by the model adapter, C and E varied, temperature, vector input. T1: grid plus c_max variation; T2: nodes plus DFN. **Codex most**.
- **PR's own tests**: C = Partial (I2 No); T1 = Partial; T2 = Partial; Table 3 = Partial. All three agree.
- **Expected violations**: none. T2 pointed out that D_e.csv is low by an order of magnitude; Codex only recorded the correction of the test expected values and did not comment.
- **Contamination**: none seen. The COMSOL example page accessed is the source cited in the source docstring, which is within what is permitted.

---

## Summary table

| PR | Primary vs T1 | vs T2 | Branch coverage (who has more / what Codex missed) | PR's own tests C / T1 / T2 (Table 3) | Expected violations (Codex, prediction) | Contamination |
|---|---|---|---|---|---|---|
| Cantera#2107 | same (strict: same) | same (strict: same) | Codex most; singles out the κ_T=0, β≠0 branch | P / P / N (P) | primary none; I4 pre-existing derivative API defect | none seen |
| ClimaAtmos#4637 | same (strict: partial) | same (strict: partial) | Codex most | P / P / P (P) | primary none; Jacobian SGS velocity block (non-primary) | none seen |
| ClimaAtmos#4822 | same (strict: partial) | same (strict: partial) | Codex most; does not write variance consistency \|ρ\|=√f | N / P / N (P) | primary none; c_Δx→0 discontinuity (convention) | none seen |
| Oceananigans#6024 | same | same | Codex most | N / N / N (N) | none | none seen |
| Oceananigans#6083 | same | same | Codex most; GM reduction promoted to primary | P / P / N (P) | none | none seen |
| CoolProp#2851 | same (strict: partial) | same (strict: partial) | Codex most | P / P / P (P) | primary none; mixture second saturation derivative (non-primary) | none seen |
| CoolProp#3339 | same | same | roughly equal | P / P / P (P) | none | none seen |
| pandapower#3016 | same (strict: partial) | same (strict: partial) | Codex most; covers YNyn, Yzn, pf_3ph, Zd | P / P / P (P) | YNyn split, Yzn base (placement suspect), pf_3ph base, all in primary; Zd convention pending; 4 further items | none seen |
| MOOSE#34020 | same (strict: partial) | same | Codex most | P / N / P (N) | primary none; friction law can be negative (conditional) | none seen |
| MOOSE#33307 | **partial** (lacks Ws_T>0) | **partial** (lacks Ws_T>0) | Codex most but no Ws sign check | P / N / P (N) | primary none (MPI to be investigated); zero-flow NaN (non-primary) | none seen |
| PyBaMM#5888 | same (strict: partial) | same (strict: partial) | Codex slightly more | P / Y / Y (P) | none | none seen |
| PyBaMM#5889 | same | same (strict: partial) | Codex most | P / P / P (P) | none | none seen |

## Counts

- **Primary invariants (lenient criterion, same as compare.md)**: vs T1 same 11, partial 1 (MOOSE#33307), different 0; vs T2 same 11, partial 1 (the same PR), different 0.
- **Primary invariants (strict criterion)**: vs T1 same 5 (Cantera, 6024, 6083, 3339, 5889), partial 7, different 0; vs T2 same 5 (Cantera, 6024, 6083, 3339, 34020), partial 7, different 0.
- **PR's own tests**: Codex gives Yes 0 / Partial 10 / No 2. The PRs for which the three tracks' judgments on the first primary invariant are fully identical number 6 (6024, 5889, 2851, 3339, pandapower, 4637). Codex agrees with the harmonised judgment of Table 3 in 9/12; the disagreements are MOOSE#34020 and #33307 (Codex records a golden regression as partial) and ClimaAtmos#4822 (Codex records No). Rereading Codex's textual descriptions under the Table 3 rule "only a golden regression counts as No", the two MOOSE PRs turn into No, giving 11/12 agreement, with only #4822 remaining. This step is a reclassification, not Codex's own wording.
- **Expected violations**: only one PR, pandapower#3016, has a primary invariant that Codex expects to fail. Six more PRs (Cantera, 4637, 4822, 2851, 34020, 33307) have expected failures only in non-primary invariants, mostly pre-existing defects or convention issues.
- **The four pandapower classes**: YNyn (a clearly expected failure, with the mechanism, and flagged as convention-dependent), pf_3ph (likely failure, with the mechanism), Yzn (expected failure from base dependence, placement "suspect"), Zd (only marked as a convention to be settled, no violation judged). All four are independently identified, consistent in direction with the audit's ruling of "2 violation classes, 2 convention-dependent classes". The (1−si) split of Yzn is not written as an explicit formula the way YNyn is.
- **Contamination**: none seen in the 12. All execs stayed in the working directory; network access was files of the merge/base commit on raw.githubusercontent, git-tree listings on api.github.com, and literature/pages cited by the PR or the source (Wagner 2025, Valøen–Reimers 2005, PowerFactory TechRef, the COMSOL example), plus one NumPy documentation access. No PR comments, issue discussions, pulls/issues endpoints of the api, or searches for others' analyses of the PR were accessed. Note: the logs contain several `web search:` entries whose query is empty (4 in 6024, 5 in 5888, etc.) and the specific pages opened were not recorded, so these cannot be fully excluded; the entries with records are all within what is permitted.

## Notes on interpretation

1. Prompt confounding: see the beginning. "Codex has the most branch coverage in 11 of 12" should be read mainly as "what can be done under the same prompt", and not as "what a change of family alone achieves".
2. Timing and training data: these PRs were merged between 2026-05 and 2026-10, and it cannot be judged from the logs whether Codex's training data contains later discussion (speculation: low risk, but unverified).
3. All of Codex's expected verdicts come from reading the source, with no compilation or execution; the agreement of the pandapower predictions with the native replay measurements is an after-the-fact comparison and does not constitute an independent measured reproduction.
