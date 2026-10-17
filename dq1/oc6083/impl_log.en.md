# oc6083 implementation log (Oceananigans#6083, ΔQ1 unified spot check)

> Purpose: record which part of `spotcheck.jl` each sentence of the three tracks' check plans maps to, and all implementation judgments.
> Structure: §0 common → §1–§3 the three tracks → §4 unimplemented items (deviation candidates) → §5 time and memory → §6 post-smoke changes.
> Update rule: after the smoke run only script errors are fixed, recorded in §6. References: the ΔQ1 protocol (preregistration); shared library `../jl_common/spotcheck_common.jl` (parameters, seed, rotation, budget, verdict, replay and JSON as in §0 of oc6024).

## §0 Common

- Environment: compute node B `env/ocean`, Oceananigans 24de7026 (**later than the merge commit 1734ff9e of this PR**; the same environment as the native replay, used per protocol §1). A native environment is available. Single thread.
- Entry point: `--tracks 1,2,3 --n N --outdir DIR [--budgets …]`, one process for the three tracks (changed from `--track/--out`, to save compilation).
- Seed `20261017 + 5×10 + track number`.
- Drift-velocity extension: following the `SinkingParticles` of the PR test file `test/unit/adaptive_implicit_vertical_advection.jl` (`required_biogeochemical_tracers`, `biogeochemical_drift_velocity`, `Adapt.adapt_structure`); in addition `TwoTracerBGC` (D has drift, E has none) is written for Track 3 to check tracer-lookup misalignment.
- Single column: `RectilinearGrid(size=40, z=(-400,0), topology=(Flat,Flat,Bounded))`, `momentum_advection=nothing, free_surface=nothing` (as in the PR test); the w field is W on interior faces and 0 on the top and bottom faces (`wfield`).

## §1 Track 1 (`audit/specs/track1.md`)

| Plan text | Code |
|---|---|
| A: WENO(order=5, AIVA(cfl=0.5)); (i) SinkingParticles drift (ii) `PrescribedVelocityFields(w=…)` (iii) `forcing=(D=AdvectiveForcing(w=…),)`; advance a few steps with QAB2 and SplitRK3 each, compare `interior(D)`, relative 1e-12; (ii) at least 1e-10 if it takes a different branch | branch `A_three_paths_{QAB2,SRK3}`: `A_forcing_vs_bgc_rel1e-12`, `A_prescribed_vs_bgc_rel1e-10` |
| B: same setup, Gaussian centred at −100 m with width 20 m, C ∈ {0.4,4,10}; Σc relative 1e-12; centroid displacement vs w_s nΔt (suggested <1e-2); with a very small cfl (almost fully implicit) the error should drop to the order of 1e-6 | branches `B_cfl0.5_*`, `B_nearly_implicit_*`; the three values of C are mandatory points (12) |
| C: 16×16×16 three-dimensional, tilted isopycnals, `IsopycnalSkewSymmetricDiffusivity(κ_skew, slope_limiter=FluxTapering(S_max), skew_flux_formulation=AdvectiveFormulation())`; read uₑ,vₑ,wₑ, the discrete divergence is zero to machine precision (relative 1e-12); check \|ψ\| ≤ κS_max on steep local slopes (call `ϵSxᶠᶜᶠ` directly); ψ=0 where b_z ≤ 0 | branch `C_GM_3d` (`divergence_stats`, `taper_stats`) |
| Operating domain w_s 1–200 m/day, C 0.1–10, GM κ_skew 100–2000, S_max 1e-3–1e-2 | each branch sampled uniformly |

Implementation judgments:
1. "Advance a few steps" and the number of steps of B: n = max(1, ⌊120 m/(CΔz)⌋), i.e. a displacement of about 120 m (as in the PR test), keeping the tracer away from the boundary (a premise of the invariant of B); A uses the same number of steps.
2. The Gaussian "width 20 m" is exp(−(z+100)²/(2·20²)), the same expression as in the PR test. The plan does not give the initial value of A; the same Gaussian is used.
3. B uses the drift path (the path changed by the PR). The "very small cfl" is 1e-6; the "order of 1e-6" takes tolerance 1e-5 (the same order of magnitude).
4. Two tolerance levels for (ii): only exceeding 1e-10 counts as a rejection, and the measured relative difference is recorded in the observations (values between 1e-12 and 1e-10 can be seen in the obs range of the JSON).
5. The relative difference is defined as max|D_a−D_b| / max|D_b|.
6. Grid and field for C: Bounded³, L_x=L_y=100 km, H=1000 m (the plan only gives the number of grid points); b = N²z − 2N²·max(z−z_u,0) + Δb·tanh((x−L_x/2)/L_f) + a similar front in y; N²=1e-5, L_f=10 km, z_u=−100 m (above it b_z<0 is an unstable layer); Δb is such that the maximum slope = steep·S_max, steep ∈ U[10,100] (the plan only says |S| ≫ S_max). The divergence is relative to max_cell Σ|A U|. ψ=κϵS, and κ·max|ϵSx| ≤ κS_max is checked (relative margin 1e-12); unstable-layer criterion: ϵS must be exactly 0 at faces higher than z_u+Δz (both cell centres in the unstable layer). `minimum_bz` takes the default 0.

## §2 Track 2 (`audit/specs/track2.md`)

| Plan text | Code |
|---|---|
| I1: model A `PrescribedVelocityFields(w)`, model B `SinkingParticles` (transport velocity 0), the same Gaussian, advance 3–20 steps, compare interior(D), rtol 1e-10; scan C ∈ {0.25,0.5,1,4,10} | branch `I1I2_{QAB2,SRK3}_{WENO5,Centered2,Upwind3}`: `I1_prescribed_vs_drift_rel1e-10`; five values of C × 6 configurations are mandatory points (30) |
| I2: relative drift of sum(D) ≤1e-12; relative error of the centroid displacement ≤5% | `I2_inventory_*`, `I2_centroid_drift_rel5pct` |
| I3: 16×16×16 tilted b field + random perturbation + a local unstable layer, after `compute_closure_fields!` read uₑ,vₑ,wₑ, compute the divergence with the δ and area operators of Operators, relative to max\|A uₑ\| ≤1e-12; call `ϵSxᶠᶜᶠ` inside a kernel to spot-check the slope upper bound; uₑ=vₑ=wₑ=0 when the isopycnals are horizontal or b_z ≤ 0 | branch `I3_GM_tilted` (divergence, \|ϵSx\|,\|ϵSy\| ≤ S_max, zero in the unstable layer), `I3_GM_flat` (all zero) |
| Domain: drift 1–200 m/day, Δz 1–20, Δt 10 min–1 day, Courant 0.1–10; WENO5 and Centered/UpwindBiased; κ_skew, S_max; b_z including negative values and zero | sampling in `t2_sampler`, `gmsamp` |

Implementation judgments:
1. The grid is fixed by the plan's method to 40 levels, 400 m (Δz=10 m, within 1–20), and Δz is not drawn separately. w_s and C are drawn, Δt = CΔz/w_s, and samples are rejected unless Δt ∈ [600 s, 86 400 s]. The number of steps is a uniform integer in 3–20.
2. Centered takes order=2, UpwindBiased order=3 (library default orders); all three are wrapped with AIVA(cfl=0.5).
3. The centroid check is made only when the predicted displacement is ≤ 200 m (final centroid not lower than −300 m); otherwise this sub-check does not appear ("away from the boundary" is a premise of the invariant).
4. GM perturbation: uniform noise of amplitude 1e-3·N²H (independent seed per point, the seed is recorded in the input); tilt M²=slope·N², slope ∈ U[1e-4,1e-1], y-component ratio U[−1,1]; the top z_u of the unstable layer lies at 50%–95% of the column depth. N²=1e-5 fixed (not given in the plan). The divergence uses, as the plan says, `δxᶜᵃᵃ/δyᵃᶜᵃ/δzᵃᵃᶜ` and `Ax_qᶠᶜᶜ/Ay_qᶜᶠᶜ/Az_qᶜᶜᶠ`.
5. The "horizontal isopycnal" branch requires bitwise 0.

## §3 Track 3 (spec.md §5, CPU part)

| Plan text (excerpt) | Code |
|---|---|
| Original PR setup: three steps of 0.4 day, −100 m/day, inventory error < T, 120 m of sinking within 5%; mirrored rising; two time steppers, two precisions | mandatory points `pr_setup` (8) |
| Two time steppers, two precisions, zero velocity, stage CFL ratios C/C* ∈ {0,0.1,0.99,1,1.01,2,8}; C* ∈ {0.25,0.5,1}; WENO5, UpwindBiased(1), Centered(2) | sampling dimension of branch `I1I2_source_matrix` |
| Drift alone, AdvectiveForcing alone, two advective forcings plus an ordinary zero Forcing, cancelling opposite drifts, both sides of the threshold, two tracers with different allocations; whole-array comparison, tolerance T·max(1,‖q0‖∞) | submodes `bgc_vs_forcing`, `two_forcings`, `cancel`, `two_tracers` |
| Below the threshold adaptive agrees with `ExplicitTimeDiscretization()` | submode `explicit_equiv` (only when C/C* < 1) |
| Check the I2 inventory in every run | `I2_inventory_*` |
| GM closure source vs providing the observed eddy-induced velocity with AdvectiveForcing, halve Δt, difference/Δt converges to 0 | branch `I1_GM_source_vs_forcing` |
| Closure nothing, κ_skew=nothing, κ=0, a single GM, tuples of length 1–5 with GM at each position (the others with zero coefficients) | branch `I1_closure_tuples` |
| I3: 16×16×32 bounded grid, linear b=N²(z−S_xx−S_yy), `SmallSlopeIsopycnalTensor(minimum_bz=m)`, `FluxTapering(S_m)`, AdvectiveFormulation, κ a positive linear function of z; interior uₑ=−κ′ϵSx etc. compared with an independent piecewise formula, relative to T; scan slope-norm ratios, quadrants, floor cases, S_x=S_y=0.8S_m; x/y swap and sign reflection; back out the slope-norm bound; Float64 parameters used on a Float32 grid | branch `I3_taper` (`taper_case`) |
| I4: horizontally periodic variation, stable stratification, variable κ; signed face-flux sum cell by cell ≤ T·Σ\|A U\|; no normal flux at top and bottom; x–z, y–z, three-dimensional; uniform and stretched z; a constant passive tracer stays constant under a steady prescribed eddy-induced velocity | branch `I4_divergence` |
| I5 (fixed-grid part): physical horizontal b=N²z+b₀ produces no GM transport | branch `I5_fixed_grid_horizontal_b` |
| Tolerances T₆₄=1e-11, T₃₂=5e-5, short runs ≤20 steps | `T3(FT)` |

Implementation judgments:
1. Stage CFL ratio → Δt: C = ratio·C*, Δt = CΔz/|W|; ratio=0 means zero velocity (W=0), with Δt taken as C*Δz/|W_sampled|. |W| ∈ U[1e-6,1e-2] m/s, sign ± (for upward motion the initial value is changed to a mirrored Gaussian centred at −300 m). Number of steps U{1..20}. The maximum of the per-stage Δt of SplitRK3 equals Δt (β is 1 in the last stage), so C is set by the whole-step Δt.
2. The allocation λ of the two forcings is ∈ U[−1,2] (W₁=λW, W₂=(1−λ)W); the zero Forcing returns 0 in discrete form.
3. GM source comparison: two-dimensional x–z grid 32×32 (100 km × 1 km, Bounded, Flat, Bounded), PrescribedVelocityFields() freezes the transport to zero, κ_skew=1000, κ_symmetric=nothing; the control model adds the same AdvectiveForcing(u=uₑ, w=wₑ) to both b and D (uₑ,wₑ being a copy of the initial eddy-induced velocity of the GM model), so that the two models are "otherwise identical". Δt₀ = 0.1·min(Δx/|uₑ|, Δz/|wₑ|), then /2, /4; "converges to 0" criterion: the difference/Δt at Δt₀/4 is not greater than half of that at Δt₀, or already ≤ T/Δt₀.
4. Tuples: the zero-coefficient closure uses `ScalarDiffusivity(ν=0, κ=0)`; 3 steps, Δt=3600 s; the D arrays are compared, tolerance T·max(1,‖D0‖∞).
5. taper: L_x=L_y=100 km, H=1000 m; κ = 1000·(1+0.5(z+H)/H) (κ′=0.5 m/s); the vertical gradient of b takes {−1e-6, 0, m/2, m, U[1e-8,1e-3]} by b_z case, and the horizontal gradient is set by the target "slope after the floor" S as b_x=−max(b_z,m)·S_x (when ≤0 after the floor, b_x takes −1e-6·S_x, and the velocities are expected to be all zero); the angle is one of eight quadrant directions or random; the parameter precision is Float64 or Float32, half each. The interior comparison zone is 2 cells from the boundary. The characteristic velocity κ′S_m is the absolute scale for the zero branch.
6. I4: stretched grid z_k = −H(1−k/Nz)^{1.5}; b = 1e-5 z + amp·1e-3·sin(2πx/L)cos(2πy/L)·e^{z/300}; the constant tracer uses the default Centered advection, 5 steps, Δt set by CFL 0.1.
7. The inventory (I2) is checked in every single-column run; the centroid is checked only at the PR-setup mandatory points (the large-domain/40-80-160 refined centroid check also written in the plan was not done, see §4).

## §4 Unimplemented items (deviation candidates)

- **All CUDA (the GPU part of I6, CPU/GPU agreement)**: compute node B has no CUDA device. The CPU part for Float32/Float64 is covered.
- I1/I2: centroid checks on a large domain and with 40/80/160 levels; comparison of UpwindBiased(1) with an independently assembled conservative tridiagonal update (the explicit fraction of AIVA depends on the details of the `clock.last_Δt` initialisation, and an independent assembly can hardly guarantee the same definition as in the library, so it was not done); the "transport alone" two-dimensional cellular circulation branch; background transport + drift + forcing all on, compared with a fine-step reference; inventory checks on stretched grids, `ImmersedBoundaryGrid(GridFittedBottom)`, variable volumes (`MutableVerticalDiscretization`/`ZStarCoordinate`).
- I3: constant / face-field representations of κ (only function κ was done); DiffusiveFormulation and VerticallyImplicit options, positive κ_symmetric, per-tracer symmetric coefficients; independent interpolation for the two staggerings in three-dimensional variable-gradient cases; negative checks rejecting per-tracer NamedTuple skew coefficients.
- I4: boundary fluxes around immersed boundaries; variable grids; convergence of advective and diffusive formulations at 16³/32³/64³ (< 2%).
- I5: variable grids (z-star, `set!(model; η=…)`) and the immersed-bottom part (the essential part); only the fixed-grid "zero fallback" was done.
- Most of these are workload issues that can be implemented by building other grids/closure combinations, not an unavailable environment; whether to make up for them is for the orchestrating agent to rule on.

## §5 Time and memory estimate (speculation)

- Basis: in native-replay job 9cd55d6e in the same environment, the two Oceananigans PRs took about 2.5 minutes, peak about 3.5 GiB.
- Track 1: A has three single columns, each ≤120 steps (40 levels) → about 0.3–1 s; B about 0.2–0.5 s; C builds a 16³ model and calls kernels face by face → about 0.5–1 s. 10⁴ points take about 1–2.5 hours, budget 3000 s (N≈4000–10⁴ expected).
- Track 2: two single columns ≤20 steps → about 0.2 s; GM about 0.5 s → 10⁴ points about 40–60 minutes, budget 2400 s.
- Track 3: the single-column matrix about 0.3 s; GM source/tuples about 0.5–1 s; taper with three 16×16×32 models about 1 s; I4 about 1 s → 10⁴ points about 2 hours, budget 3000 s.
- Compilation: 6 schemes × time-stepper combinations, two precisions, several closures and grid topologies; the first compilation is estimated at 10–15 minutes in total.
- Peak memory estimated 4–6 GiB, declared 10 GiB. The full job is estimated at 9000 s, hard limit 14400 s.

## §6 Post-smoke change record

(empty)
