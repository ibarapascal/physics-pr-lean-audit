# ca4637 implementation log (ClimaAtmos#4637, ΔQ1 unified spot check)

> Purpose: record which part of `spotcheck.jl` each sentence of the three tracks' check plans maps to, and all implementation judgments.
> Structure: §0 common → §1–§3 the three tracks → §4 unimplemented items (deviation candidates) → §5 time and memory → §6 post-smoke changes.
> Update rule: after the smoke run only script errors are fixed, recorded in §6. References: the ΔQ1 protocol (preregistration); shared library `../jl_common/spotcheck_common.jl`.

## §0 Common

- Environment: compute node B `env/ca`, ClimaAtmos ba2bf122 (**later than the merge commit a8efce5a of this PR**; the same environment as the native replay, used per protocol §1). A native environment is available. Single thread.
- Entry point: `--tracks 1,2,3 --n N --outdir DIR [--budgets …] [--precision Float64|Float32]` (changed from `--track/--out`, one process for the three tracks to save compilation). Seed `20261017 + 2×10 + track number`.
- **Two processes** (run serially in the same job command): `--precision Float64` runs only the F64 branches and writes `ca4637_t<track>.json`; `--precision Float32` runs only the F32 branches and writes `ca4637_f32_t<track>.json`. Reason: F64 and F32 each need about 4 ClimaAtmos configurations, and one process holding the compiled artifacts of 8 configurations is estimated to exceed the roughly 24 GiB job limit of compute node B (the earlier measured peak with 3 configurations was about 11 GiB). The two processes have the same seed but different branch sets, so the sampling sequences differ; the two JSONs together are all the points of that track. **This is an implementation judgment that needs the orchestrating agent's approval.**
- Configurations (following the PR test fixtures):
  - plane: `ConstantBuoyancyFrequencyProfile`, Schar, x_max=100 km, 8 spectral elements, z_elem=20, dz_bottom=500 m, `disable_surface_flux_tendency=true`; the dry mode (`microphysics_model=dry`) is used for mass conservation, 1M (`implicit_microphysics=false`) for the χ test.
  - EDMF: `DYCOMS_RF01` single column, prognostic EDMFX, 1M, z_max=1.5 km, 30 levels, dt=120 s, ARS222, `toml/prognostic_edmfx.toml`, the other switches as in `tracer_mass_consistency_tests.jl`.
- "Advanced state": each base state is first advanced a few steps (5 for plane, 10 for EDMF, as in the PR test), then 20 snapshots are stored step by step; each point draws a snapshot number, restores the state, and calls `CA.set_precomputed_quantities!` so that the cache is consistent with that snapshot, then modifies Y as the plan says and calls the tendency functions (as in the PR test, the cache is not recomputed after Y is modified). In this way a rejected point can be replayed on its own from the snapshot number.

## §1 Track 1 (`audit/specs/track1.md`)

| Plan text | Code |
|---|---|
| Following advection_tests and tracer_mass_consistency_tests, build the simulation with `CA.AtmosConfig` (plane, Schar, Linear, z_max 21 km) | `plane_cfg`, `edmf_cfg` |
| A: read the top level of `p.precomputed.ᶠu³` and `CA.compute_ᶠuₕ³`, check the ratio | `A_u3_{bottom,top}_le_100eps_uh3` |
| Then advance 1000 steps with a dry model, surface flux off, no microphysics, record `sum(Y.c.ρ)` as M at each step, relative drift ≤1e-12 (F64) and not growing linearly with time; compare with SLEVE | branch `A_mass_{Linear,SLEVE}_*`: `A_mass_rel_drift_max`, `A_mass_drift_linear_trend_over_1000_steps` |
| B: on an advanced state set ρχ ← cρ, χʲ ← c (c = 1 and 1e-3), call `implicit_tendency!`, `explicit_vertical_advection_tendency!`, `edmfx_sgs_mass_flux_tendency!` respectively, compare Yₜ.ρχ with c·(corresponding ρ path), tolerance 100ε·max\|ρ tendency\| | branches `B_grid_paths_*` (implicit, explicit), `B_sgs_*` (SGS) |
| Then the stepping version: uniform χ initial value, microphysics sources and sinks off, after N steps max\|χ−c\|/c ≤ N·100ε | branch `B_timestep_*` (**implementation deviation, see judgment 5**) |
| C: in the same call compare `parent(Yₜ.c.ρ) == parent(Yₜ.c.ρq_tot)` | branch `C_sgs_rho_eq_rhoqtot_*` |
| Operating domain: terrain 0–3 km, z_top 20–45 km, wind 0–30 m/s, Δt 1–120 s; EDMF single column Δt 60–120 s; F64 and F32 | sampling |

Implementation judgments:
1. For A, at each point a new dry simulation is built with the drawn z_top ∈ U[20,45] km and Δt ∈ U{1…120} s (integer seconds, since the configuration only accepts "N secs"), and the wind speed rescales the uniform initial 10 m/s by U/10 (U ∈ U[0,30]). When Δt is large the explicit horizontal acoustic CFL is exceeded and the run diverges → recorded as a failure (not a violation).
2. **The Schar mountain height cannot be changed through the configuration** (`ScharTopography` defaults to h_max=25 m and is not exposed as a configuration item), so the 0–3 km of the operating domain cannot be sampled and the default is fixed; the Agnesi mountain was not done (changing the terrain type requires recompilation). Deviation candidate.
3. The plan gives no F32 drift tolerance: the protocol default relative 1e-6 is used. "Not growing linearly with time" is judged as: the least-squares slope of the relative drift against the step number × 1000 steps ≤ the same tolerance.
4. B: following the PR test, all water species (q_tot, q_lcl, q_icl, q_rai, q_sno) are set to cρ at the same time (in EDMF the subdomain χʲ is also set to c). "Corresponding ρ path": implicit uses `Yₜ.ρ` (the ρ row of the implicit tendency); explicit uses the PR test's −ᶜadvdivᵥ(ᶠinterp(ρJ)/ᶠJ·ᶠu³); SGS uses `Yₜ.ρ` of the same call. The "max\|ρ tendency\|" of SGS is the max\|Yₜ.ρ\| given by the SGS call on the **unmodified state** (after being made uniform it should be 0 and cannot serve as a scale). The tolerance is not scaled with c (as in the original text).
   - Speculation: at c=1e-3 the implicit path may show a difference of (1−c)× the settling term (the library adds the settling flux of water equally to ρ and ρq_tot); this would be a "plan wording vs library convention" issue, for ruling.
5. **The stepping version cannot follow the original text**: 1M has no switch to turn off microphysics sources and sinks, and condensate with χ≡1 or uniform c=1e-3 would be changed by microphysics. Changed to: q_tot is set to a tiny constant c=1e-7 (subsaturated everywhere, condensation/settling sources are zero), condensates are 0, the cache is recomputed and N=20 steps are taken (N is not given in the plan), checking max\|q_tot−c\|/c ≤ N·100ε. This is a deviation, a §7 candidate.
6. Sponge layer: off by default (the plan reminds that it must be off or restricted below the sponge layer).

## §2 Track 2 (`audit/specs/track2.md`)

| Plan text | Code |
|---|---|
| I1: advance with `CTS.step!`, compute the total mass with ClimaCore `sum(Y.c.ρ)` at each step, compare after 100–1000 steps, tolerance 1e3·N·eps·M(0); top face at each step ≤100·eps·max\|uₕ³\| | branch `I1_mass_{Linear,SLEVE}_*` |
| It is suggested also to test Sleve, higher mountains, dt 1–10 s, Float32 | Sleve and F32 done; dt ∈ U{1..10}; mountain height see §1 judgment 2 |
| I2: after initialisation ρq_* ← ρ, subdomain q_* ← 1, microphysics and sponge off, advance N=50 steps, check max\|ρq_x/ρ − 1\|, tolerance 1e-10 (F64) | branch `I2_timestep_{plane,edmf}_*` (**the same deviation as §1 judgment 5**) |
| I3: call `edmfx_sgs_mass_flux_tendency!` to get Yₜ, integrate vertically by column, relative to max\|Yₜ.ρq_tot\|·Σ JΔz ≤ 1e-12; and Yₜ.ρ − Yₜ.ρq_tot = 0 | branch `I3_sgs_column_*` |

Implementation judgments:
1. I1 uses the dry mode (the invariant requires precipitation settling and all sources and sinks to be off); the other settings of the PR fixture are unchanged (z_max 21 km, wind 10 m/s); the number of steps is U{100…1000}. "Top-face check at each step": after each step `set_implicit_precomputed_quantities!` is called (the state filter is idempotent) and then ᶠu³ is read, and the worst ratio over the whole run is compared with 1.
2. I2 stepping version: same as §1 judgment 5 (only q_tot, c=1e-7); EDMF uses a configuration of the same type with surface flux off (`disable_surface_flux_tendency=true`, otherwise the surface vapour flux itself breaks the uniformity), the rest as in the fixture. The F32 tolerance is 1e-4 (the plan gives only F64; scaled up in proportion to eps, an implementation judgment).
3. I3: the column integral uses ClimaCore `sum` (WJ-weighted); Σ JΔz uses `sum` of an all-ones field. ρ, ρq_tot and the four condensates are all checked; the second equation is required to be bitwise equal as "= 0".

## §3 Track 3 (spec.md §5)

| Plan text (excerpt) | Code |
|---|---|
| Start from the two PR fixtures; switch off physical sources to test transport | `plane1M`, `edmf` base states |
| I1: after giving winds of both signs call `CA.set_implicit_precomputed_quantities!`, read ᶠu³ and each ᶠu³ʲ at `CA.half` and `nlevels+half`, residual ≤ 100ε·max(\|uₕ³\|,\|g³³u₃\|), exactly 0 for zero scale; a warped-terrain regression must have a nonzero horizontal projection | branch `I1_impermeability_*` (Linear, SLEVE, EDMF single column); U=0, ±10 are mandatory points |
| I2: zero buffer; density reference −advdiv(ρu³); call `explicit_vertical_advection_tendency!` for each passive field, `CA.vertical_transport` for total water; C=0,1,1e-6,0.01; ‖A_χ−CA_ρ‖∞ ≤ 200εS·\|C\|; with an artificially nonzero boundary velocity D₀ is still structurally zero flux | branch `I2_constant_tracer_*` |
| I3: set all subdomains and the mean to constants and call `edmfx_sgs_mass_flux_tendency!`, zero within 500εS_SGS; unequal velocities with varying area; non-uniform scalars also zero at equal velocities; buffer unchanged when switched off | branch `I3_null_sgs_*` (three modes) |
| I4: non-uniform humidity, nonzero SGS water flux, starting from a zero buffer Yₜ.ρ and Yₜ.ρq_tot are bitwise equal; repeat with the same nonzero initial buffer | branch `I4_water_air_*` |
| I5: use the ClimaCore volume integral, \|ΣVT\| ≤ 500NεΣV\|T\|; 20/40/80 levels and a stretched grid | branch `I5_closed_column_*` |

Implementation judgments:
1. I1: plane uses the dry-mode Linear/SLEVE fixtures (Δt=1 s, z_max 21 km), each point restarts from a copy of the initial value and rescales the horizontal wind by U/U₀ (U ∈ U[−100,100]); EDMF uses the snapshot state with the same rescaling. \|g³³u₃\| is taken as \|u³−uₕ³\|.
2. I2: S is taken as max\|A_ρ\| (the plan also writes "the larger of the sum of absolute face fluxes/volume", which was not computed separately, an implementation simplification); "configured reconstruction" takes `p.atmos.numerics.energy_q_tot_upwinding` and a central scheme is also tested; artificial boundary velocity: after setting the top level of the cached ᶠu³ to 1e-2, the volume integral of the explicit condensate tendency is checked to be ≈ 0 (tolerance as in I5). The constant number density (2M) is not applicable to 1M and was not tested.
3. I3: S_SGS is the max(\|Yₜ.ρq_tot\|,\|Yₜ.ρ\|) of the SGS call on the unmodified state, and max\|Yₜ.ρe_tot\| for the energy row; the condition-number factor max(1,ρ/m₀) is neither computed nor reported (simplification). "Equal velocity" sets the cached ᶠu³ʲs and ᶠu³⁰ to ᶠu³. "Switched off" replaces `edmfx_model.sgs_mass_flux` with false (copying the structure and replacing the field), presets the buffer to 1, and checks that it is bitwise unchanged.
4. I4: the nonzero initial buffer is a U[0,1) random field with an independent seed per point, identical in the ρ and ρq_tot rows. `report_I4_flux_nonzero` is also recorded (record only, to guarantee non-emptiness).
5. I5: plane checks ρ advection, q_tot central transport and the explicit advection of the four condensates; EDMF checks the ρ and each water-species row of the SGS; the number of EDMF levels is taken from {30 (fixture), 40, 80}, stretching on/off (the plan writes 20/40/80, but the single-column fixture itself has 30 levels, so 30 is used in place of 20).

## §4 Unimplemented items (deviation candidates)

- Track 1: terrain-height sampling (not exposed in the configuration), Agnesi; the stepping-version χ test follows the original text (see §1 judgment 5; a feasible variant was used instead).
- Track 2: "higher mountains"; the original I2 with χ≡1 and all water species stepped for 50 steps (as above).
- Track 3: flat plane, 3-D warped box; `turbconv` none/EDOnlyEDMFX, `updraft_number=2`; 0M, 2M, 2MP3 microphysics and number-density tracers; custom passive tracers; scans of the options `tracer_upwinding`/`energy_q_tot_upwinding`/SGS reconstruction/`edmfx_tracer_upwinding` (each needs a new configuration and recompilation); 20 seeded realisations with random density and area (replaced by 20 snapshots); sampling near the environment small-area fallback threshold in I3; all of the I6 Jacobian (`CA.Jacobian`, `update_jacobian!`, `first_column_block_arrays`, directional differences).
- Reason: number of configurations and compilation memory/time; I6 is a large amount of work and depends on internal block names. Whether to make up for them is for the orchestrating agent to rule on.

## §5 Time and memory estimate (speculation)

- Basis: native replay job 9cd55d6e on the same node: the ClimaAtmos#4637 part (plane dry/1M and EDMF single column, EDMF stepping 10 steps) takes about 2.5 minutes, RSS peak about 11 GiB.
- Track 1: A builds a simulation per point (about 10–30 s) + 1000 steps (plane 8×4 columns×20 levels, ARS343 implicit, estimated 20–60 ms per step) → about 0.5–1.5 minutes per point; B/C about 0.2–2 s per point. A accounts for 1/3 of the points → 10⁴ points would need tens of hours, so truncated by budget: F64 2400 s, F32 1200 s (about 60–150 points each expected).
- Track 2: I1 builds a simulation per point + 100–1000 steps → about 0.3–1 minute; I2 stepping 50 steps about 5–15 s; I3 about 0.2 s. Budgets F64 2000 s, F32 1000 s.
- Track 3: mostly single tendency calls (about 0.1–0.5 s); the EDMF level-count variants of I5 each need a simulation build + 30 steps the first time. Budgets F64 1500 s, F32 800 s.
- Compilation: about 4 configurations per process, 1–2 minutes each.
- Peak memory: estimated 13–16 GiB per process (extrapolated from 11 GiB for 3 configurations), declared as 20 GiB (two processes run serially, the peak is the larger one). The full job is estimated at 10 800 s, hard limit 14 400 s; the smoke run (F64 400 s per track, F32 300 s per track) is estimated at 3600 s.

## §6 Post-smoke change record

(empty)
