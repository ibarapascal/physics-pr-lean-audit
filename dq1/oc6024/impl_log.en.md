# oc6024 implementation log (Oceananigans#6024, ΔQ1 unified spot check)

> Purpose: record which part of `spotcheck.jl` each sentence of the three tracks' check plans maps to, and all implementation judgments.
> Structure: §0 common → §1–§3 the three tracks (plan sentence → code; implementation judgments) → §4 unimplemented items (deviation candidates) → §5 time and memory estimate.
> Update rule: after the smoke run only script errors are fixed, with changes appended to §6; the plan mapping is not changed. References: the ΔQ1 protocol (preregistration); shared library `../jl_common/spotcheck_common.jl`.

## §0 Common

- Environment: compute node B `env/ocean` (Oceananigans installed with `Pkg.add(rev="24de702689")`, i.e. merge commit 24de7026 of #6024; the same environment as the native replay). Julia 1.11.7, `JULIA_NUM_THREADS=1`. A native environment is available.
- Entry point: `julia --project=env/ocean dq1/oc6024/spotcheck.jl --tracks 1,2,3 --n N --outdir DIR [--budgets s1,s2,s3]`; one process runs the three tracks in turn and each track writes `DIR/oc6024_t<track>.json` (`--track k --out PATH` is also accepted). **The entry arguments were changed from `--track/--out` to `--tracks/--outdir`, so that the three tracks do not each pay a compilation cost.**
- Seed: track t uses `20261017 + 4×10 + t` (Xoshiro); the inputs of all points of the same track are drawn in order from this single generator, so the input of point i is determined only by the seed and i.
- Equal split over branches: the N points rotate over the k branches by `mod1(i, k)` (⌊N/k⌋ or one more per branch). Deterministic checks with no input to draw (parameter algebra) do not enter the rotation and are run once as mandatory points (implementation judgment: repeating 5000 times gives exactly the same result and carries no information).
- Time budget: `--budgets` gives seconds per track, and the run stops when exceeded (what has been computed is a prefix of the same sampling sequence, and the branches are still roughly equally split); the JSON records `n` (actually run), `n_planned`, `stopped_by_budget`. This is how "reduce to the largest N that can be completed" of protocol §3 is implemented.
- Verdict: each point has several checks; any exceeding its tolerance = rejection; exception/NaN/non-convergence = failure (not a violation); nothing to check = empty (n_empty). The first 20 rejected points are recomputed once, and `replay_identical` records whether the two observations are bitwise identical. Each check also records the observed minimum/maximum and the maximum of residual/tolerance.
- When the plan gives no tolerance, relative 1e-6 and absolute 1e-12 are used.

## §1 Track 1 (`audit/specs/track1.md`, "numerical check method")

| Plan text | Code |
|---|---|
| 1. Parameter algebra: `c = CATKEVerticalDiffusivity()`, compute `Cᵇ (Cˡᵒu³/CˡᵒD)^{1/4}` and compare with 0.4, tolerance 1e-3 | `t1_param_eval`; mandatory point `param_algebra` |
| 2. Single column: `RectilinearGrid(size=Nz, z=(-H,0), topology=(Flat,Flat,Bounded))`, `HydrostaticFreeSurfaceModel(grid; closure=CATKEVerticalDiffusivity(), buoyancy=BuoyancyTracer(), tracers=…)`, b≡0; set the top and bottom momentum flux of u to −u★², f=0 | `stress_column(Float64, H, Nz, ustar)` |
| Integrate until Kᵤ is steady (relative change < 1e-4 / inertial time) | `t1_steady`: every `T1_CHUNK` seconds check the maximum relative change of κu in the check zone |
| Take the face values of `model.closure_fields.κu`, divide by u★·d_face; the bottom zone should be ≈0.40, the top ≈0.468, tolerance 5% | `t1_column_eval`: `bottom_ratio_0.40`, `top_ratio_0.468` (the face with the largest deviation in the zone) |
| Check zone d ∈ [2Δz, 0.2H], Cᵇd < Cˢ(H−d); exclude 1–2 cells next to the boundary | `t1_zones` (d ≥ 2Δz excludes the bottom two cells; the top layer is the mirror image) |
| Halve Δz twice and see whether it converges to 0.4 | `bottom_refine_err_nonincreasing` |
| Usual operating domain u★ ∈ [1e-3,3e-2], H ∈ [50,5000], Δz ∈ [1,10] | branch `column_steady`, uniform sampling |

Implementation judgments:
1. The plan writes `tracers=(:b,:e)`; the constructor of the merge commit adds `:e` automatically and **rejects it when given explicitly** (`ArgumentError`), so `tracers=:b` is used. The keyword `grid` written in the plan is a positional argument in this commit, and is passed positionally.
2. "Inertial time": with f=0 there is no inertial period; 2π/f≈62 832 s for f=1e-4 s⁻¹ is used as the check interval. If 1e-4 is still not reached after at most 100 segments (≈73 days), it is recorded as a failure.
3. Time step Δt=600 s (a value commonly used in Oceananigans CATKE single-column examples; not given in the plan); the time discretisation uses the CATKE default (vertically implicit) and the default QAB2.
4. The initial state is the library default (u=v=e=0); no equilibrium TKE is preset, so as not to bias the initial value towards the answer.
5. The "top zone" is not defined in the plan; depth ∈ [2Δz, 0.2H] symmetric to the bottom layer with Cˢ·depth < Cᵇ(H−depth) is used.
6. 5% is interpreted as the relative deviation of each face |K/(u★d)/0.40−1| ≤ 0.05 (over all faces in the zone).
7. "Halve Δz twice and see whether it converges": on the same physical check zone defined by the coarse grid, the maximum relative deviation of the Δz/4 grid is not larger than that of the original grid (or already ≤1e-3); all three grids are integrated to steady state on their own. Done only for the bottom layer (the plan only says converge to 0.4).
8. Nz = round(H/Δz), actual Δz = H/Nz; lower limit Nz 4. If both the bottom and top check zones are empty (e.g. H small and Δz large) the point is recorded as empty.

## §2 Track 2 (`audit/specs/track2.md`, "numerical check method")

| Plan text | Code |
|---|---|
| I1: `c = CATKEVerticalDiffusivity()`, read the parameters and compute directly, tolerance 5e-4 | `t2_param_eval` (mandatory point): κ_b against 0.4, κ_s against 0.468 |
| I2: the same single column as Track 1, `FluxBoundaryCondition(-u*²)` for u at top and bottom, b constant, advance to steady state (> 10 eddy turnover times H/u*) | `t2_run`: integrate 10·H/u★ |
| Take κu and S of ∂z(u), compute the median of K_u/(u*h) and u*/(hS) over 5–15 grid points near the bottom, tolerance ±0.02 | `bottom_K_over_ustar_h`, `bottom_ustar_over_hS` |
| Check that the deviation decreases as Δz is halved | `bottom_refine_dev_nonincreasing` |
| Also check the near-surface 0.468 as a control | `top_K_over_ustar_d_0.468`, `top_ustar_over_dS_0.468` |
| Usual operating domain u* 1e-3–3e-2, H 20–200, Δz 0.25–4, f=0 or 1e-4, Float32/Float64 | 4 branches `column_{Float64,Float32}_f{0,1e-4}` |

Implementation judgments:
1. `tracers=:b` (same as §1 judgment 1). With rotation `FPlane(f=1e-4)` is used.
2. Δt = clamp(H/(50u★), 1, 600) s, guaranteeing at least about 500 steps within 10H/u★; not given in the plan.
3. The initial state is all zero (same as §1 judgment 4).
4. "5–15 grid points near the bottom": the face height h = 5Δz…15Δz (coarse grid), which must be bottom-limited (Cᵇh < Cˢ(H−h)); refined grids take values at the same physical height. The top control is the mirror image.
5. With rotation the plan notes "the log layer holds only within h ≪ u*/f": faces with h ≤ 0.1·u★/f are taken ("≪" taken as 0.1 times); if no face can be checked, the point is recorded as empty. **Whether this sentence should serve as a screening condition is for the orchestrating agent to rule on** (without screening, the rotation branches would be rejected in large numbers, and the cause is already stated in the plan).
6. Invariant I2 also writes e = u*²/√(Cˡᵒu CˡᵒD), but the "numerical check method" does not list this item and gives no tolerance; not checked.
7. S is taken as √((∂z u)²+(∂z v)²) (v≠0 with rotation).

## §3 Track 3 (spec.md §5)

| Plan text (excerpt) | Code |
|---|---|
| I1 parameters: F64/F32 × explicit/vertically implicit, κ_b against 0.4 absolute 2e-4; against the rounded predicted value relative 1e-12/2e-6; surface separately; overriding Cᵇ=0.28 gives 0.1158743582 | `t3_param_eval`, 4 mandatory points |
| I1 constitutive: three log-equally-spaced u★, uniform e=u★²/√(ab), bottom-limited faces with d_b/H ≤ 0.01, API viscosity vs κ_b u★ d_b; near the top mirror with Cˢ | branch `I1_viscosity` (`viscosity(model)`) |
| I1 manufactured log law u=(u★/κ)log(d/d0), away from the first and last cells, K_u\|D_zu\|/u★² within 1% where Δz/d ≤ 0.05, error decreases with refinement | branch `I1_logprofile` (Nz=32/64/128) |
| Central dissipation: stationary neutral uniform TKE, zero flux, one very small step; Δt such that Δt√e0/ℓ_D and the diffusive Courant number are < 1e-3, then halve twice; F64 relative error < 0.5% and converging; exclude the bottom cell; two time discretisations × tke_time_step=nothing/two substeps | branch `I1_dissipation` (`t3_decay`, `decay_dt`) |
| I2: compare new/old Cᵇ in the same state and a scan of positive coefficients in [0.1,2]; three diagnostic fields and dissipation; strict bottom/surface/stratification limits (10% margin); no cap/depth cap/finite cap; ratio assertions only for no cap; relative 1e-10/1e-4, absolute 1e-13/1e-8; dissipation ratio 0.5% | branch `I2_compare` |
| Ri taken as −0.1,0,0.1,0.254,0.764,1.274,2 and 1e-5 on both sides of the three thresholds; stable zero shear, neutral zero shear; shear directions (S,0),(0,S),(S/√2,S/√2) | `RI_LIST`, `ri_state`, `shear_dir` |
| I3: all wet faces finite/non-negative/bounded, boundary mask zero; TKE floor cases and each cap; K_c(D_zc)² and K_uS² non-negative; closed column variance non-increasing in one step (explicit diffusive Courant ≤ 0.1, then implicit and half step), tolerance 100ε; frozen-coefficient variance identity relative 1e-10/1e-4; two time steppers | branch `I3_bounds_variance` |

Implementation judgments:
1. `set!(...; u=z->S*z, c=z->z/H, ...)` and `tracers=(:b,:c)` from plan §5 are adopted as written; `momentum_advection=nothing, tracer_advection=nothing`.
2. I1 constitutive: with Nz ∈ {32,64,128} only Nz=128 has interior faces with d_b/H ≤ 0.01, and the others are recorded as empty points (per the plan, the grid is not changed). No tolerance is written; relative 1e-6 is used.
3. Log law d0=1 m (the derivative is independent of d0); "error decreases with refinement" compares Nz=64 and 128 on the same interval d ≥ 20H/64 (with Nz=32 there is no bottom-limited face satisfying Δz/d ≤ 0.05).
4. The predicted central dissipation uses the plan's formula CˡᵒD e0^{3/2}/L_c, L_c = min(Cˢ d_s, Cᵇ d_b) (in this interval ℓ_D is not truncated by H); Δt is half of the smaller of the two constraints. "Two substeps" uses `tke_time_step = Δt/2`. Only F64 is done (the plan gives 0.5% only for F64). Convergence criterion: the error at Δt/4 is not larger than that at Δt (rounding allowance 1e-12).
5. The "§1 formula" of I2 is outside the readable range (only §2, 4, 5 may be read); the ratio formulas of §2, K₂/K₁ = ℓ_D₂/ℓ_D₁ = C₂/C₁, ε₂/ε₁ = C₁/C₂, and the monotonicity formula are used instead. Strict bottom limit: C₂d_b < 0.9·min(Cˢd_s, ℓᴺ), ℓᴺ = √e0/√N²; strict other limit: min(Cˢd_s, ℓᴺ) < 0.9·C₁d_b (in which case invariance is required). The "depth cap" is judged by whether the library output ℓ = K/w★ is close to H (> 0.9H), and only monotonicity is checked in this case. The finite cap is q×Cˡᵒu·C₂·0.5H·√e0, q ∈ U[0.2,5] (the plan says only "below/above the predicted value").
6. Discrete values of Ri: the plan writes discrete values and the operating domain gives only the intervals of N² and S, so for each point Ri is drawn first (uniform from the discrete table), then S is drawn in S ∈ [1e-5,1e-1] and N²=Ri·S² is set; if N² falls outside [−1e-5,1e-3], N² is drawn instead and S is solved backwards.
7. The dissipation ratio is done on a stationary neutral state with the same e0 (with shear the production term is nonzero, and the dissipation cannot be inferred from the decay).
8. I3: TKE floor cases e0 ∈ {0, emin(1−1e-3), emin, emin(1+1e-3)} and [1e-8,1e-2]; caps ∈ {none, finite, 0}. The bounded upper limit is min(K_max, H·w_f), with relative allowance 1e-6. The frozen-coefficient identity is done only on the explicit time discretisation: after one QAB2 step `model.timestepper.Gⁿ.c` is the initial-state tendency (verified: in this commit QAB2 computes the tracer tendency at the start of `ab2_step!`, and after the step `update_state!` only computes the momentum tendency). Closed column: explicit Δt = 0.1Δz²/max κ; implicit uses the same Δt and half of it; QAB2 and SplitRK3 are each run. The tracer initial value is c = z/H of the plan's `set!` example.

## §4 Unimplemented items (deviation candidates, for the orchestrating agent to rule on)

- Track 3 I2: `ImmersedBoundaryGrid(GridFittedBottom)` (including a raised flat bottom, two column depths, active-cell map on/off); horizontally arrayed closures (different Cᵇ per column); tuples of CATKE with a zero-diffusivity closure.
- Track 3 I3: the immersed-boundary fixture within "for each test fixture" (not built, as above).
- Track 3 "uniform and smoothly stretched" grids: only uniform grids were done.
- Track 3 I1 constitutive check: "two time discretisations" and F32 were both done; the dissipation check was done only for F64 (as in the plan).
- No GPU (compute node B has only CPUs); the plan for this PR does not require a GPU.

## §5 Time and memory estimate (speculation, to be replaced by measurements after the smoke run)

- Basis: in the earlier native-replay job (9cd55d6e) on this node in the same environment, the two Oceananigans PRs took about 2.5 minutes in total, RSS peak about 3.5 GiB (including compilation).
- Track 1: each point integrates three grids to steady state, mean Nz≈650, ×1/2/4, 1–3 ms per step, typically 10 segments × 105 steps → about 5–60 s per point; 10⁴ points would take about 50–150 hours, far over 4 hours, so truncated by budget at 3600 s (N≈100–400 expected).
- Track 2: each point has three integrations of about 500 steps, Nz ≤ 3200 → about 2–5 s per point; 10⁴ points about 6–14 hours, budget 3000 s (N≈600–1500 expected).
- Track 3: each point builds 3–16 small models (about 50–200 ms each) → about 0.3–2 s per point; 10⁴ points about 1–5 hours, budget 2400 s.
- Compilation: closure-type combinations of F64/F32, explicit/implicit and various caps; the first compilation is estimated at 5–10 minutes in total.
- Peak memory: estimated 3–5 GiB, declared 8 GiB.
- The full job is estimated at 9600 s (sum of the three tracks' budgets + compilation), hard limit 14400 s; the smoke run gives 600 s to each track. If the orchestrating agent wants each track to get 4 hours, it can be split into one job per track (the same command with `--tracks k` changed).

## §6 Post-smoke change record

(empty)
