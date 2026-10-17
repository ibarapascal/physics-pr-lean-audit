# MOOSE#33307 ΔQ1 spot-check implementation log

> Purpose: record the sentence-by-sentence mapping of the three tracks' check plans to `spotcheck.py` / `scm_port.py`, and all implementation judgments, for review under §6/§7 of the ΔQ1 protocol.
> Structure: §0 environment (non-native) → §1 Track 1 → §2 Track 2 → §3 Track 3 → §4 time and memory → §5 deviation candidates → §6 post-smoke corrections → §7 to be ruled on (classification of smoke rejections).
> Update rule: after the smoke run only script errors may be fixed; changing the plan mapping must be recorded in protocol §7 and appended to §5.
> References: upstream: the protocol, the #33307 sections of `audit/specs/track1.md`/`track2.md`, `audit/specs/track3_guided/idaholab_moose_33307/spec.md` §2/4/5, the source of merge commit b5a1d9e; downstream: `dq1/outputs/mo33307_t{1,2,3}.json`.

## §0 Environment: non-native recomputation (deviation candidate)

- Read-only check of compute node A and compute node B: there is no MOOSE/libMesh/PETSc under `env/`, neither `subchannel-opt` nor `mpiexec` is on PATH, and there is no moose directory within 4 levels of the home directory. Per instruction, MOOSE was not built.
- Instead, the changed functions are recomputed line by line in `scm_port.py` from the C++ source of merge commit b5a1d9e (the coverage list is in the file header):
  - Changed: `SCMFrictionUpgradedChengTodreas` (constructor checks/warnings, triangular and square lattice friction factors), `SCMMixingChengTodreas` (constructor warnings, `computeMixingParameter`, `computeSweepFlowMixingParameter`), `SCMHTCBorishanskii`, `SCMHTCGnielinski`, `SubChannel1PhaseProblem::computeBulkReynoldsNumber` (single-process semantics), the nrings/pitch/dwire checks of `SCMTriAssemblyMeshGenerator` (the default tolerance of `absoluteFuzzyGreaterThan` is libMesh::TOLERANCE² = 1e-12, `MooseUtils.h` of the framework was read).
  - Unchanged but called: `SCMHTCClosureBase` (preinfo, blend, computeHTC), subchannel areas, wetted perimeters and component hydraulic parameters of `TriSubChannelMesh`/`QuadSubChannelMesh`.
  - Cannot be recomputed: the cross-flow estimate indices of `implicitPetscSolve` (#33340, needs a PETSc solver context), the whole `externalSolve` flow, MPI (the `processor_id()!=0` branch).
- `mooseError` → `MooseError` exception (counted as a failure); `flagSolutionWarning` → warning list.
- Runtime environment: any Python 3.10+ with numpy; the JobSpec uses `env/py2/bin/python`.
- Seed 20261017 + 9×10 + track number. Geometry defaults are taken from `test/tests/problems/SFR/sodium-19pin/test19_implicit.i` (nrings=3, pitch=0.012, D=0.01, dwire=0.002, hwire=0.0833, flat_to_flat=0.056 → wall gap 0.0022154 m); the side-gap ratio of the square lattice is taken from `verification/friction_model_verification/two_channel2.i` (side_gap/D=0.1).
- Subchannel local Re (the plan gives no range): for uniform mass flux Re_i = Re_b·D_h,i/D_h,b (Tracks 1 and 2); the Track 3 plan requires local and component Re to be sampled independently of each other.

## §1 Track 1 (`audit/specs/track1.md`)

Branches: subchannel types CENTER/EDGE/CORNER (k=3, ⌊N/3⌋ each).

| Plan text | Code |
|---|---|
| Usual operating domain P/D∈[1.05,1.42], H/D∈[8,52], D_w/D∈[0.05,0.25], N_pin∈{7,19,37,61,217}, Re_b∈[50,1e6], Pr∈[0.003,0.03], Pe∈[10,2200] | `t1_sample`: continuous quantities linearly uniform; N_pin uniform over the 5 values |
| A: H/D = 8, 20, 52, Re_b≈5×10⁴, with wire wrap; `wire_diameter`=1e-6 m to run a nearly bare rod; ΔP_wire ≥ ΔP_bare and increasing monotonically as H/D decreases; strict inequalities | closure-level substitute: `t1_A_f_wire>=f_nearbare(Dw=1e-6m)_same_mdot` (Re_b, Re_i converted at the same inlet mass flow), `t1_A_f_monotone_increasing_as_HoD_decreases(8,20,52)`; mandatory points 3 types × H/D{8,20,52} × Re_b=5e4 |
| A (invariant): correction factor [1+Ws_T(A_r/A'_p)tan²θ]^1.41 ≥ 1, laminar factor ≥ 1 ⇔ Ws_T>0 (edge/corner) | `t1_A_WsT>0`, `t1_A_turbulent/laminar_wire_factor>=1` (only in the EDGE/CORNER branches) |
| B: Re_b=Re_L(1±1e-3), (1±1e-4), relative difference of ΔP <10δ; the same for Re_T | closure-level substitute: Re_i fixed, `t1_B_f_jump_at_ReL/ReT_delta=…<10delta` |
| B (invariant): f=f_L(Re_i) when ψ=0, f=f_T(Re_i) when ψ=1; the same for Cm, Cs; all subchannels share the same ψ | `t1_B_f(Reb=ReL)==fL`, `…==fT`; endpoints and jump of Cm (CENTER branch, center–center gap)/Cs (EDGE branch edge–edge, CORNER branch corner–edge); `t1_B_same_psi_all_channel_types` |
| C: `pin_HTC_closure=borishanskii`, low flow Pe<200, back out Nu and compare with 24.15·log₁₀(poly), relative 1e-6 | sampled points with Pe≤200: Re = Pe/Pr, `computeNusseltNumber` read directly; `t1_C_…`; for P/D∈[1.1,1.5] also check Nu∈[4.25,10.83] |
| D: compare the explicit/implicit/monolithic versions | cannot be recomputed, listed in `not_covered` |
| E (invariant): for Re>10³, Pr>0, Nu_T is finite and >0 and continuous in Pr | points with Re_i>1e3: `t1_E_Gnielinski_NuT_finite>0`, `…_Nu_finite>0`, `t1_E_Nu_continuous_in_Pr_delta=…<10delta` (friction factor taken from this closure) |
| Suggestion: 1 process and 4 processes agree | cannot be recomputed, `not_covered` |

Implementation judgments:
1. The upper limit 0.25 of D_w/D is greater than the lower limit 0.05 of P/D−1: sampled points with D_w > P−D (about 27%) raise an error at the mesh-generator check, and are counted as "failure" per native behaviour.
2. When switching between a nearly bare rod and H/D levels the same inlet mass flow and viscosity are kept: Re_b' = Re_b·(A/D_h)·(D_h'/A'), Re_i' = Re_b'·D_h,i'/D_h,b'.
3. The "relative difference of ΔP" of B is replaced by the relative difference of the friction factor; "continuous" in E has no criterion, so the <10δ of B of the same track (δ=1e-3, 1e-4) is borrowed.
4. The local field of the mixing closure (affects only local-Re warnings) is taken as uniform mass flux, μ=2.5e-4 Pa·s.
5. In the mandatory points, Pe is 1000 (does not trigger C), Pr=0.005, P/D=1.2, D_w/D=0.2, N_pin=19 (template values).

## §2 Track 2 (`audit/specs/track2.md`)

Branches: CENTER/EDGE/CORNER (k=3).

| Plan text | Code |
|---|---|
| Operating domain P/D 1.0–1.42, H/D 8–52, N_pin 7–217, Re_b 50–1e6; liquid metal Pr 0.004–0.03 and water; Borishanskii 1.1≤P/D≤1.5, Pe 200–2000 | `t2_sample`: N_pin taken from all valid hexagonal numbers within 7–217 {7,19,37,61,91,127,169,217}; P/D for Borishanskii is drawn separately |
| I1: scan H/D ∈ {8,15,25,52}, with/without wire wrap (d_wire=0), pressure drop with wire ≥ without; no turbulent wire correction error | closure level: `t2_I1_f_wire>=f_bare_same_Re` (edge/corner), `t2_I1_no_turbulent_wire_correction_error`; mandatory points 3 types × H/D{8,15,25,52} |
| I1 (invariant): Ws_T ∈ [0.12, 9.07] > 0; edge/corner correction factor ≥1; center f_wire = f_bare·(P_w,p/P_w) + a positive drag term | `t2_I1_WsT_*`, `t2_I1_*_wire_factor>=1`, `t2_I1_center_drag_terms>0`, `t2_I1_center_CfT==…` |
| I2: Re_b at ±0.1% on both sides of Re_L and Re_T, jump ≤1e-3 relative | `t2_I2_f_jump_at_*_pm0.1pct<=1e-3`; the same for Cm/Cs |
| I2 (invariant): ψ=0→f_L, ψ=1→f_T; f≥0; formulas for Re_L, Re_T; same ψ in the same layer; Cm, Cs continuous at the endpoints | `t2_I2_*` |
| I3: compare the three solves; single process and mpiexec -n 2 | cannot be recomputed, `not_covered` |
| I3 (invariant): Nu=24.15·log₁₀(…) for Pe≤200, ≈4.25 at P/D=1.1, ≈10.8 at 1.5 | the sampled points are evaluated at Pe=200 (the lower end of this track's Pe interval): `t2_I3_Borishanskii_Pe=200_…`; mandatory points P/D∈{1.1,1.5} × Pr∈{0.004,0.03} × type, tolerance half a unit of the last digit (0.005, 0.05) |

Implementation judgments:
1. The plan gives no D_w: the template test19 relation D_w = P−D is used (the template 0.002 = 0.012−0.01); this way the mesh check always passes.
2. "Without wire (d_wire=0)": setting only d_wire=0 raises an error in the friction-closure constructor (it requires both to be 0), so a bare rod with d_wire=hwire=0 is used; compared with the wire-wrapped case at the same (Re_i, Re_b) (the plan's formula is f(Re) against f(Re)).
3. The plan is internally inconsistent: the Borishanskii operating domain has Pe 200–2000, while the formula of I3 and the check method target Pe≤200 ("set a low Pe"). Handling: Pe is sampled per the operating domain (recorded only), and the Pe≤200 formula is checked at the interval end Pe=200, with Re = 200/Pr.
4. Water: the operating domain says only "and water" and gives no Pr interval; Borishanskii is used only for liquid metals, and no water branch was added.
5. Local Re as in Track 1 (uniform mass flux).

## §3 Track 3 (spec.md §2/§4/§5)

Branches (k=20, 500 points each): friction 9 branches (triangular bare rod/triangular wire-wrapped/square lattice × CENTER/EDGE/CORNER), mixing 2 branches (Cm: gaps with a center channel; Cs: edge/corner gaps), Borishanskii 3 branches (type), Gnielinski 6 branches (triangular wire-wrapped/square lattice × type).

| Plan text | Code |
|---|---|
| "Write a small test application that calls the public C++ interface; do not substitute a Python transcription for the closure" | Not feasible (no build): the closures are recomputed line by line from the source, and the independent oracle is written separately from the documented formulas + coefficient tables (`o_*` functions); tolerance per the plan \|y−y_ref\| ≤ 1e-12+1e-10\|y_ref\| |
| I1: friction CENTER/EDGE/CORNER; triangular bare rod and wire-wrapped; both sides of and equality with table switching 1.1 (offset 1e-6); H/D=8,20,52; d_w/(P−D)=0.1,0.5,1; R_b ∈ {0.5R_L, R_L, √(R_LR_T), R_T, 2R_T} × R_i ∈ {0,0.1,1,50,10³,10⁵} | sampled branches `F_*`; mandatory points `grid_*` (810), `switch_*` (27) |
| I1: endpoint values checked directly; with the same geometry and same R_b, changing the local flow does not change ψ; ψ∈[0,1] in the transition region; no assertion of a convex combination | `t3_I1_f(Reb=ReL)==fL`, `…ReT…`, `t3_I1_psi_independent_of_local_Re`, `t3_I1_psi_in_[0,1]` |
| I1: independently accumulate the inlet μ_i·m_i and compare with R_b; uniform/non-uniform viscosity and flow distribution at the same (M, μ_b) | triangular branch: a random distribution + random viscosity gives R_b, compared with the independent formula; another distribution + uniform viscosity should give the same R_b |
| I1 mixing: call both functions; center/edge/corner combinations; zero branches; endpoint bounds; changing the local field does not change it; ring counts 2,3,9, a single ring as a separate probe; Cm, Cs constants | `M_*` branches: `t3_I1_mixing_beta_vs_oracle`, `CL<=C<=CT`, endpoints, `nondecreasing` (R_b×1.001), `other_mixing_function_exact_zero`, `independent_of_local_fields`; mandatory points `mix_*` (60), `mixing_one_ring_probe` |
| I2 Borishanskii: r=1.1,1.2,1.35,1.5; R_i>R_T fixed; Pe=1,199.999,200,200.001,500,2200; compare absolute values; for Pe>2200 look only at the warning; used as a duct closure should raise an error; dimensional check of k; length scaling 0.1, 10; nominal and deformed Dpin; local R_i swept across the R_L, R_T endpoints | `B_*` branches (r, Pe uniform, R_i∈(R_T,1e6] uniform, Pr=Pe/R_i) + mandatory points `boris_*` (84, including the Pe=2500 warning probe), `boris_as_duct_closure` |
| I3 Gnielinski: pin and duct sentinels; triangular and square lattices; Pr = 1e-5,0.001,0.005,0.01,0.5,1,10,2000; laminar/transition/turbulent local Re and independent component segments; lower-bound threshold f_c(1±1e-6) and f_c; passivity for R_L<R_i<1000 | `G_*` branches (Pr∈[1e-5,2000] linearly uniform): oracle, finiteness, Nu≥0 (including the mixed region), upper bound for R_i≥max(R_T,1000), equal to the original formula when not truncated, a warning when truncated, rejection by computeHTC of h·D_h/k=Nu or negative values, duct sentinel; mandatory points `gn_*` (144), `gn_floor_s=…` (solve R_i backwards so that f=f_c(1+s)) |
| I4: cross-flow estimate correlation matrix/mass conservation (needs running subchannel-opt) | cannot be recomputed, `not_covered` |
| I5: P/D=1, 1−1e-6, 1+1e-6, 1.1 × d_w=0, P−D, P−D+0.5ε_g, P−D+2ε_g; read ε_g first; one size-scaling example | mandatory points `geom_*` (16 cases each for D=0.01 and 0.1), ε_g=1e-12 |
| I6: finite plateau for local R_i=0,1e-12,0.1,1,1+1e-6; overall inlet scaling t=1,1e-3,1e-6,0 | the plateau is checked at every friction sample point; mandatory points `bulk_scale_t=…` |

Implementation judgments:
1. Sampling intervals: friction P/D, W/D∈[1.01,1.4] (the mathematical check interval used by the plan for "both sides of table switching"), H/D∈[8,52] (the conservative intersection of the two sources), N_pin 7–217 (intersection), D∈[0.005,0.015] m, d_w/(P−D)∈(0,1], R_b and R_i each [50,1e6]; mixing P/D∈[1.067,1.35], H/D∈[4,52], R_b∈[400,1e6], ring counts∈{2,3,9}; Borishanskii r∈[1.1,1.5], Pe∈[0,2200]; Gnielinski Pr linearly uniform (the plan does not say log-uniform).
2. The side gap of the square lattice is not given: side_gap/D=0.1 of two_channel2.i is used; the area/wetted perimeter required by Borishanskii use the triangular bare-rod geometry, with the wall-gap ratio taken from the template.
3. The value of the deformed Dpin is not given in the plan: 0.98·D is used.
4. The friction closure for Gnielinski uses triangular wire-wrapped (template test19 configuration) or square-lattice bare rod.
5. R_b monotonicity compares the two points R_b and 1.001·R_b (the plan gives no step).
6. At t=0 the 0/0 of computeBulkReynoldsNumber gives NaN in C++ (no error); the recomputation likewise gives NaN → counted as "failure" per protocol §4, listed separately in the report.
7. "Transient inlet changes", "serial vs 2 ranks" and "zero-flow channel in implicit solver" in the plan are solver behaviours and are not covered.

## §4 Expected time and peak memory (speculation, not run)

| Track | Per point | 10⁴ points | Basis |
|---|---|---|---|
| t1 | ~1–2 ms (about 30 closure calls) | ~20 s | pure Python scalars |
| t2 | ~1 ms | ~10 s | |
| t3 | ~1–3 ms | ~30 s + about 1.2k mandatory points | |

Peak memory <200 MB (numpy only); the JobSpec declares 512 MiB.

## §5 Deviation candidates

- **Non-native**: the nodes have no MOOSE build, and the changed functions are recomputed in Python from the C++ source of the merge commit (coverage list in §0); checks that depend on the full solver (component pressure drop, agreement of the three solves, MPI, cross-flow indices) were not tested.
- For the handling of the internal inconsistency of the Track 2 plan (Pe interval vs the Pe≤200 formula) see §2-3.

## §6 Post-smoke corrections (2026-10-17, after the N=100 smoke run; only script errors are fixed, the plan mapping is unchanged)

- `spotcheck.py`, Track 1 B (about lines 141–148) and Track 2 I2 (about lines 239–245): in the jump check the "fixed Re_i" was originally taken as `re_local(g, t, Reb)` (the local Re at the sampled Re_b, which can reach 1e6), but f was evaluated at Re_b=Re_L(1±δ)≈300–800, producing an inconsistent state with Re_b≈Re_L and Re_i≫Re_L, which does not match the state described by "Re_i = Re_b·D_h,i/D_h,b (uniform mass flux)" in §0 of this log and by "adjust the inlet mass flow so that Re_b = Re_L(1±δ)" in the plans of both tracks (numerical check method of Track 1 B and Track 2 I2) (the wrong objects were compared). Changed to take `Rei0 = re_local(g, t, R0)` at each endpoint R0 and keep it fixed on both sides ±δ (the mapping "fixed Re_i" is unchanged). The endpoint equality check `f(Reb=ReL/ReT)==fL/fT` uses `Rei0` as well (an identity, the conclusion is unaffected).
- Expectation: the jump residuals of the first round (relative 0.3–28) mainly come from this inconsistency; after the correction the magnitude drops to about 1e-2 (δ=1e-3)/1e-3 (δ=1e-4), but rejections at Re_L will remain, for the reason in §7-1 (the plan itself). Verification is left to `job_smoke2.json` (runs only t1 and t2, output `dq1/outputs/smoke2/`).
- "t2 rejections 118 > n=100" is not a counting error: `total_rejected` = 100 sampled-branch points + 18 mandatory points (of the 24 mandatory points: all 12 friction mandatory points are rejected at the Re_L jump, and the 6 Borishanskii Pr=0.03 mandatory points see §7-2). Track 1's 82 = 73 + 9 mandatory is the same. Not changed.

## §7 To be ruled on (classification of smoke rejections, scripts not changed)

1. **(P) Track 1 B "relative difference of ΔP <10δ, shrinking linearly with δ" and Track 2 I2 "±0.1% jump ≤1e-3" contradict the library's interpolation formula at Re_L**. Evidence: at all points in the first round the ratio of the jumps at δ=1e-3 and δ=1e-4 is constantly 2.15≈10^{1/3}, i.e. the jump ∝ δ^{1/3}; in the transition region of `SCMFrictionUpgradedChengTodreas` of the merge commit f = f_L(1−ψ)^{1/3}(1−ψ⁷) + f_T ψ^{1/3} (the Track 2 plan I2 itself also writes this formula), ψ = ln(Re_b/Re_L)/ln(Re_T/Re_L) ≈ δ/ln(Re_T/Re_L), so there is an extra f_T·ψ^{1/3} term to the right of Re_L: the function is continuous (all endpoint equalities pass) but not Lipschitz, and as δ→0 the jump converges as δ^{1/3} rather than δ. Hand calculation (CENTER bare-rod coefficients, P/D=1.2, Re_i=Re_L, speculated magnitude): f_L≈0.195, f_T≈0.049, ψ(δ=1e-3)≈3.0e-4 → jump≈0.015>0.01; δ=1e-4 → ≈0.008>0.001. On the Re_T side (1−ψ)^{1/3}(1−ψ⁷)≈7(1−ψ)^{4/3}, which is differentiable, and in the first round the Re_T jumps have 0 rejections, consistent with this. Where the plan is wrong: Track 1 "converges linearly as δ→0", "< 10·δ"; Track 2 "jump ≤1e-3 relative" (at ±0.1%).
2. **(P) Track 2 I2 Cm, Cs "±0.1% jump ≤1e-3"**: the mixing parameter C = C_L + (C_T−C_L)ψ^{2/3} (original text of the Track 2 plan); at δ=1e-3, ψ^{2/3}≈4.5e-3, and in the first round the observed relative jumps are 0.0036–0.0037 (all at Re_L), exceeding 1e-3. The same kind as 1: continuous but converging as δ^{2/3}. Under Track 1's tolerance 10δ, Cm/Cs are not rejected, also consistent with this.
3. **(P) Track 1 C and Track 2 I3 "Nu = 24.15·log₁₀(poly) for Pe ≤ 200"**: in the first round Track 1 had 6/8 rejected, Track 2 68/112 sampled and 6/12 mandatory (observed Nu is below the formula value, e.g. 7.22 against 10.16 at P/D=1.418, Pr=0.0247). Cause: `SCMHTCClosureBase`, not changed by the merge commit, blends piecewise in the local Re: for Re ≤ Re_L the laminar Nu (3.52–3.73) is taken, for Re_L<Re<Re_T it blends linearly by (Re−Re_L)/(Re_T−Re_L), and in the transition region the turbulent Pe is computed with Re_T (Re_T·Pr can be >200, with the Pe correction term added on top). The plan's Pe (≤200) and Pr (0.003–0.03) make Re = Pe/Pr ∈ [333, 6.7e4], and a large part falls at Re < Re_T (Re_T≈1.2e4–2.2e4 for P/D 1.1–1.5), where the Nu returned by the closure (also the Nu backed out in the full solver) is not 24.15·log₁₀(poly). The equality holds exactly only for points with Re ≥ Re_T (no correction term at Pe ≤ 200, blending weight 1), and the points that passed in the first round are of this kind (mandatory Pr=0.004 → Re=5e4 all pass, Pr=0.03 → Re≈6.7e3 all rejected, consistent with this). Where the plan is wrong: Track 1 C and Track 2 I3 treat the NuT formula as the whole output of the closure and miss the library's laminar/transition blending (the interval check "Nu∈[4.25,10.83]" of Track 2 has 3 rejections for the same cause).
4. **(P) Track 1: 27 failures (not rejections)**: all are `MooseError: The wire diameter must not exceed the pin-to-pin gap`. The plan's usual operating domain gives P/D∈[1.05,1.42] and D_w/D∈[0.05,0.25] independently, whereas the mesh-generator check newly added by this PR requires D_w ≤ P−D; about 27% of the rectangle of the operating domain falls on geometry explicitly forbidden by the PR. Counted as failures per protocol §4, not as violations; whether to remove this part from the denominator in the full 10⁴ run is for the orchestrating agent to rule on (script not changed).
5. No (V) class was seen for this PR in the first round.
