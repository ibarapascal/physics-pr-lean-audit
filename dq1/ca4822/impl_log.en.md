# ca4822 implementation log (ClimaAtmos#4822, ΔQ1 unified spot check)

> Purpose: record which part of `spotcheck.jl` each sentence of the three tracks' check plans maps to, and all implementation judgments.
> Structure: §0 common → §1–§3 the three tracks → §4 unimplemented items (deviation candidates) → §5 time and memory → §6 post-smoke changes.
> Update rule: after the smoke run only script errors are fixed, recorded in §6. References: the ΔQ1 protocol (preregistration); shared library `../jl_common/spotcheck_common.jl`.

## §0 Common

- Environment: compute node B `env/ca`, ClimaAtmos ba2bf122 (i.e. the merge commit of #4822; the same environment as the native replay). A native environment is available. Single thread.
- Entry point: `--tracks 1,2,3 --n N --outdir DIR [--budgets …]`, one process for the three tracks (changed from `--track/--out`). Seed `20261017 + 3×10 + track number`.
- Model-level configuration (implementation judgment): microphysics `1M`, `cloud_model=quadrature`, `use_sgs_quadrature=true`, `tq_correlation_model=diagnosed` (the main path of Track 3; Tracks 1 and 2 do not specify microphysics, so this is used uniformly; under 1M, T is solved directly from the internal energy with no saturation adjustment, which makes it easy to manufacture states).
  - plane: x_max=2000 km, 8 spectral elements (Δx_h = 250/3 ≈ 83 km), z_max=2 km, 10 uniform levels (Δz=200 m).
  - box (Track 3): 200 km × 200 km, 4×4 spectral elements (Δx_h ≈ 16.7 km), z_max=2 km, 20 levels (Δz=100 m).
  - One set each for Float64 and Float32 (4 configurations in all; plane has two more sets with 4/16 spectral elements for Track 2's resolution convergence).
- Manufacturing a state (`manufacture!`): given θ_li(x,y,z) and q_tot, ρ keeps its initial value, condensates and wind are 0; T is solved back from θ_li (fixed-point iteration calling `TD.liquid_ice_pottemp`), ρe_tot = ρ(e_int(T,q)+Φ), then `CA.set_precomputed_quantities!`.
- The horizontal fields use a "sawtooth": within each spectral element the slope is ±1, adjacent elements have opposite signs, and the field is continuous at element boundaries (periodically continuous for an even number of elements), so the gradient is constant within each element and the element derivative, the `lumpedₕ` restriction and the weighted DSS are exact for |∇ψ|² and ∇ψ₁·∇ψ₂; the analytic gradient can thus serve as an exact reference, and the affine jump at the periodic seam is also avoided (Track 3's plan states this explicitly).
- Parameter switches: `with_sgs` copies `p.params.sgs_quadrature_params` and replaces fields (c_Δx, c_Δz, C_cov=`diagnostic_covariance_coeff`, ρ_turb=`Tq_correlation_coefficient`), `with_corr` replaces `AtmosWater.tq_correlation_model`; afterwards only `CA.set_covariance_cache!(Y, p2, thermo_params)` is called to recompute the covariance (precomputed is shared). In this way the simulation need not be rebuilt for each parameter value.
- Levels: the first and last level are excluded (the vertical gradient at boundary layers is one-sided/a boundary value), and Track 3 excludes the first and last two levels (the plan says interior levels).

## §1 Track 1 (`audit/specs/track1.md`)

| Plan text | Code |
|---|---|
| Pointwise: call `CA.isentropic_slope_ratio`, `CA.isentropic_cell_cap`, `CA.isentropic_gradient_invariant` directly | `r_of`, `f_of`, `G_of` |
| A: draw ∇_hθ, θ_z ∈ [θz_min,1e-2], b at random, set ∇_hq=b∇_hθ, q_z=bθ_z, check I ≤ 10ε(inv_q + r²inv_θ) | branch `Ai_affine_{F64,F32}` |
| Then for arbitrary ∇q check that the difference between I(∇q+b∇θ, q_z+bθ_z) and I(∇q,q_z) is ≤ 10ε·(magnitude of each term) | branch `Aii_gauge_*` |
| B: compute f and ℓ²=fΔx² over a grid of s, check (sℓ)²+fΔz²=Δz² (1e-6/F32, 1e-14/F64) and the projection upper bound | branch `B_cellcap_*` |
| Whole field: plane configuration, diagnosed, state q_tot = a + b·θ_li, call `CA.set_covariance_cache!`, compare q′q′ for c_Δx equal to 0 and nonzero, difference ≤ 100ε·q′q′ | branch `field_affine_cdx_*` |
| C: for a general state read ᶜT′q′, T′T′, q′q′ (before clamp), check (T′q′)² ≤ T′T′·q′q′·(1+100ε) | branch `field_CS_*` |
| Operating domain Δx 1–250 km, Δz 20–1000 m, θ_z −5e-3–1e-2, \|∇_hθ\| 0–1e-4, q_z −1e-5–0, \|∇_hq\| 0–1e-8; F32 and F64 | sampling |

Implementation judgments:
1. The range of b in A is not given: q_z ∈ [−1e-5,0] is drawn and then b = q_z/θ_z is set; the "arbitrary b" of A(ii) draws another q_z in the same way. ∇ is a two-dimensional vector with direction uniform in [0,2π).
2. The "magnitude of each term" in A(ii) is the larger of inv_q + 2|r·inv_θq| + r²inv_θ over the two evaluations.
3. The tolerance of the "projection upper bound" in B is the same as for the identity (relative 1e-6/1e-14).
4. Whole field: θ_z is taken from [θz_min,1e-2] (the domain of validity of A; the plan itself states that A does not hold for θ_z < θz_min). To keep the state self-consistent (q non-negative, \|∇_hq\| within the domain), points are drawn by rejection with \|b\|·\|∇_hθ\| ≤ 1e-8; q is 0.022 at the ground reference θ=300 K. The "nonzero" value of c_Δx is the library default 1. The compared quantity is the maximum over all interior nodes of \|Δq′q′\|/q′q′. The other parameters are library defaults (turbulence block on, but with zero wind the Smagorinsky mixing length is 0).
5. C: for a general state θ and q are independent (the sawtooth phase of q is shifted by half an element, with random sign), max (T′q′)²/(T′T′q′q′) is compared with 1+100ε. What is read is the precomputed ᶜT′q′ (the quantity before clamping of the correlation coefficient; the upper limit on `σ_q` scales T′q′ and q′q′ in proportion and does not change this ratio).

## §2 Track 2 (`audit/specs/track2.md`)

| Plan text | Code |
|---|---|
| One layer: sample the domain at random with the exported functions (F64/F32), assemble the horizontal block c_g c² Δx² cap·inv and the vertical block c_g c² Δz² cap·q_z², check that the sum = c_g c² Δz² q_z², rtol 1e-12/1e-5 | branch `pt_I1_*`; also tested for θ_z ≤ θz_min: `pt_I1_floor_*` |
| I2: random F′, inv_q=F′²\|∇θ\|², inv_θq=F′\|∇θ\|², q_z=F′θ_z, inv ≤ 4ε·inv_q | branch `pt_I2_*` (with 0<f≤1, (sℓ)²+Δz²f=Δz²) |
| Two layers: plane, c_Δx=c_Δz, turbulence block off; q_tot varies only with z, θ surfaces tilted; call `set_covariance_cache!` and compare the read ᶜq′q′ with c_g c² Δz² (∂_z q)², relative 1e-3 and check convergence with resolution; avoid the σ_q cutoff | branch `model_I1_*`, `model_I1_floor_*` |
| I3: construct q = a + bθ_li as above, read T′T′/q′q′ (θ-based contrast), verify q′q′/θ′θ′ = f·b², and report the deviation from b² | branch `model_I3_*` (`report_I3_ratio_over_b2` only records, does not judge) |
| Operating domain Δx_h 25–250 km, Δz 50–1000 m, θ_z 1e-3–1e-2 (also tested ≤θz_min), \|∇_hθ\| 0–1e-4, q_z −5e-6–0, \|∇_hq\| 0–1e-7 | sampling |

Implementation judgments:
1. The value of c is not given: the library defaults c_Δx=c_Δz=1 are used; c_g takes the default 1/12.
2. The θ_z range of the floor branch is not given in the plan: [−1e-3, 1e-3] is used (the same order of magnitude as θz_min).
3. The F′ of I2 is given by the drawn q_z/θ_z (consistent with the operating domain). The additional checks of f and the identity come from the original text of invariant I2; the method section gives no tolerance, so relative 1e-6 is used.
4. Two layers: the turbulence block is turned off with `diagnostic_covariance_coeff = 0`; the resolutions are plane with 4/8/16 spectral elements (Δx_h ≈ 167/83/42 km, all within the domain); the "convergence" criterion is that the error at 16 elements is not larger than at 4 elements (or already ≤1e-12). F32 uses only the 8-element level. The vertical-gradient reference uses the analytic q_z (a linear profile is exact at interior levels). σ_q cutoff: nodes whose expected value is ≥ 0.9·(0.5q)² are excluded.
5. I3: θ′θ′ = T′T′/(∂T/∂θ)², with the Jacobian from the in-library `CA.compute_∂T_∂θ!`; f is computed from the analytic \|∇θ\|, the in-library Δx_h and Δz; the tolerance of the two-layer case, 1e-3, is reused; the worst relative deviation over all interior nodes is taken. To keep q self-consistent, points are drawn by rejection with \|b\|\|∇θ\| ≤ 1e-7.

## §3 Track 3 (spec.md §5, box grid)

| Plan text (excerpt) | Code |
|---|---|
| Public construction (AtmosConfig/get_simulation), read the `.compute(Y,p,t)` of the registered diagnostics `env_temperature_variance`, `env_q_tot_variance`, `env_q_tot_temperature_covariance`, `env_q_tot_temperature_correlation` | `diag_vals` |
| Interior elements and interior levels use locally affine fields, no affine jump at the periodic seam; the reference uses the analytic gradient; τ=1e-10/2e-5; the scale of G is S_G = C+2\|rB\|+r²A; PSD with C² ≤ V_TV_q + τ max(V_TV_q, C²) | sawtooth fields, `interior_mask`, `τ_of` |
| I1: C_cov=0, c_Δz=0, c_Δx=1; a>d; two signs of b, zero gradient, q=αθ+β, near cancellation, angles 0/π/2/π and random; V_q = h f G, V_T = 0, diagnosed cross-covariance = 0; rotation by 90° and axis reversal; below the floor use the regularised formula | branch `I1_isentropic_contrast`; angles 0, π/2, π are mandatory points |
| I2: c_Δz=0, ∇q≠0, b=0, back out f from V_q/(hC); x = 0,1e-3,0.1,1,10,100, a = −d,0,d(1±1e-3),2d,10d; check f and the offset inequality | branch `I2_slope_cap` (operating-domain sampling) + 72 mandatory points (x×a×precision) |
| I3: two correlation models, ρ_turb=−1,−0.5,0,0.5,1, block switches, turbulence on/off; no cutoff/cutoff/q=0; compare with the full transform + cutoff reference, then check PSD and correlation bounds; across V_TV_q = eps² | branch `I3_realizable`; five values of ρ_turb × two models and both sides of eps² are mandatory points |
| I4: the four zero/nonzero combinations of c_Δx and c_Δz, flat and steep θ surfaces, c = 0.5,1,2; V_q = c_g b²Z² when C_cov=0, ∇_hq=0, both scales =1; full vertical block when c_Δx=0; zero variance and fallback correlation when both are 0; compare V_T with J²v a² | branch `I4_partition` |
| I5: A>0, b≠0, c_Δz=1, C_cov=0; c_Δx = 1,1e-1,…,1e-6,0; compare the q variance and covariance at the limit and at 0 | branch `I5_cdx_continuity` (`report_I5_relative_jump_qvar` only records) |

Implementation judgments:
1. The reference formulas are assembled independently from the invariants of §2 and the description of §5 (`reference_cov`): V_θ = c_g(c_ΔzZ)²a²; V_q = h f G + c_g(c_ΔzZ)² f b²; C_θq = c_g(c_ΔzZ)² f a b; V_T = J²V_θ, C_Tq = J·C_θq; the σ_q upper limit r_max=0.5 scales in proportion with T′q′. When c_Δx=0, f does not act. J is taken from the in-library `compute_∂T_∂θ!` (the plan writes "finite J from the library's thermodynamics").
2. For zero inputs such as "zero gradient": the element derivative of a horizontally constant field is only accurate to rounding, so, following the plan's "otherwise the same scaled absolute tolerance", a rounding-level absolute floor is given (of the order 100·h·(ε·θ/Δx_h)²).
3. State self-consistency: points are drawn by rejection so that θ ∈ (250, 320) K and q ∈ [0, 0.025]; the upper bound of \|∇_hq\| is tightened to 1e-7 accordingly (with the plan domain up to 1e-6, q would be negative). Mandatory points of I2 that are not self-consistent are recorded as empty points (x=100 with a=10d, etc.).
4. Rotation by 90° and axis reversal: the state is rebuilt with (α+π/2, β+π/2) and (−α, −β), and V_q is compared (in the sawtooth basis this is exactly the corresponding rotation/reflection).
5. I3 nonzero turbulence: with zero wind the Smagorinsky mixing length is 0, so a linear shear wind u = S·z (S ∈ U[0,1e-2] s⁻¹) is added, the kinetic energy is added to the energy, and precomputed is recomputed; the nonzero-turbulence case has no independent reference (it needs the mixing length), so only PSD, bounds and fallback are checked. Nonzero C_cov is taken from U[0, 4.2] (within twice the default 2.1). The q=0 case is about 10%. "Just truncated / strongly truncated" is not constructed specifically; it is covered naturally by sampling and visible in the observations.
6. Both sides of eps²: only the vertical block is kept, with a=b so that V_TV_q is about 0.5 and 2 times eps(FT)² (rough estimate with J≈1).
7. The criterion of I5, "lim_{c_Δx→0⁺} Σ = Σ(0)", is written as a check (tolerance τ); the plan itself says the expected jump is (1−f)vb², which is a conditional criterion reported separately.

Possible systematic rejections (speculation, for interpretation): in the F32 reference comparison of Track 3 (`ref_VT` etc.), θ_li is recomputed from T in F32, with a single-point rounding of about 3e-5 K, which gives an absolute error of about 3e-7 K/m in the vertical gradient over Δz=100 m; relative to \|a\| it falls below τ₃₂=2e-5 only when \|a\| > 0.04 K/m. If the F32 reference checks are rejected on a large scale, this may be a case of "the plan tolerance is unattainable in F32", for ruling.

## §4 Unimplemented items (deviation candidates)

- Track 3: single-column and spherical grids, flat spheres and terrain-following grids (only box was done); 0M, 2M, 2MP3 microphysics; `edonly_edmfx` and `prognostic_edmfx` of `turbconv`; `sgs_distribution` = gaussian/lognormal/mean; the various alternative configurations of `uses_covariances` (including the MLCloud network); configuration-construction probes for d=0, d<0; a one-step `CA.solve_atmos!` smoke test; the "ratio after inferring the lower bound for nonzero b" of I2 under different c_Δx/c_Δz; condensate (liquid/ice/mixed phase/precipitation) states: the manufactured states use only 1M states without condensate.
- Track 1: sphere configuration (only plane was done).
- Reason: every configuration requires recompiling ClimaAtmos (about 1–2 minutes, several GiB of memory), and full coverage would exceed the 24 GiB of compute node B and the 4 hours of a single job; whether to make up for them is for the orchestrating agent to rule on.

## §5 Time and memory estimate (speculation)

- Basis: the ClimaAtmos part of native replay job 9cd55d6e on the same node (two or three configurations of plane and EDMF single column, each advanced a few steps) takes about 2.5 minutes, RSS peak about 11 GiB.
- Track 1: pointwise checks take about microseconds; each whole-field point about 0.1–0.3 s (manufacturing the state + two covariances); 4 of the 10 branches are whole-field → 10⁴ points take about 15–25 minutes (including compilation), budget 2400 s.
- Track 2: model_I1 F64 takes about 0.5 s per point for the three resolution levels, the rest about 0.2 s → 10⁴ points about 40–60 minutes, budget 3000 s.
- Track 3: box takes 1–3 state manufacturings + up to 8 covariances per point → about 0.5–1.5 s → 10⁴ points about 1.5–4 hours, budget 3000 s.
- Compilation: 4 configurations (plane/box × F64/F32) about 5–10 minutes.
- Peak memory is estimated at 12–16 GiB (4 lighter non-EDMF configurations), declared as 20 GiB. The full job is estimated at 9000 s, hard limit 14400 s.

## §6 Post-smoke change record

(empty)
