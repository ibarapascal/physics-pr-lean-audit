## 1. What changed

PR #6024, “Fix CATKE's bottom mixing length coefficient Cᵇ”, changes the default `CATKEMixingLength{FT}.Cᵇ` from **0.28 to 0.967** in `src/TurbulenceClosures/turbulence_closure_implementations/TKEBasedVerticalDiffusivities/catke_mixing_length.jl`. It adds a docstring explaining the intended bottom similarity constant and updates the displayed default in `docs/src/models/turbulence_closures.md`. No numerical algorithm, signature, or test file is changed. Explicit user overrides of `Cᵇ` remain possible.

The affected, unchanged consumers in `catke_mixing_length.jl` are `stable_length_scaleᶜᶜᶠ` and `stable_length_scaleᶜᶜᶜ`. At faces and centers respectively they compute

\[
 L=\min\{C^s d_s,C^b d_b,L_N\},\qquad
 L_N=\begin{cases}w/\sqrt{N^2},&N^2>0,\\+\infty,&N^2\le0,\end{cases}
\]

where distances are to the surface and local bottom, and the turbulent speed is constructed from \(w=\sqrt{\max(e,e_{\min})}\), with face interpolation where required. A NaN length falls back to the wall-distance minimum.

`momentum_mixing_lengthᶜᶜᶠ`, `tracer_mixing_lengthᶜᶜᶠ`, and `TKE_mixing_lengthᶜᶜᶠ` then use

\[
 \ell_\psi=\min\{H,\max[\sigma_\psi(Ri)L,\ell_{\mathrm{conv},\psi}]\},
 \quad K_\psi=\min\{K_{\psi,\max},w_f\ell_\psi\},\quad\psi\in\{u,c,e\}.
\]

Here \(u\) denotes viscosity for both horizontal velocity components, \(c\) ordinary tracer diffusivity, and \(e\) TKE diffusivity. The diffusivity wrappers are in `catke_vertical_diffusivity.jl` in the same directory. In `catke_equation.jl`, `dissipation_length_scaleᶜᶜᶜ` uses

\[
 \ell_D=\min\{H,\max[L/\sigma_D(Ri),\ell_{\mathrm{conv},D}]\},
 \qquad \epsilon=e^{3/2}/\ell_D\quad(e>e_{\min}).
\]

Thus the change increases bottom-limited mixing and decreases local dissipation at a fixed state; it does not directly change the bottom momentum boundary condition. Clipping, stratification, surface limitation, and convection can suppress its effect.

## 2. Invariants

**I1 — Primary: neutral wall similarity and the intended bottom calibration.** Let \(a=C_u^{lo}>0\), \(b=C_D^{lo}>0\), and \(C>0\) denote the active wall coefficient. In a neutral, nonrotating, locally equilibrated, constant-stress layer, assume no buoyancy production, negligible TKE transport, \(e>e_{\min}\), no active length/diffusivity caps, and \(L=Cd\). With kinematic stress magnitude \(u_*^2\), the closure must satisfy

\[
 K_u=aCd\sqrt e,\quad \epsilon={b e^{3/2}\over Cd},\quad
 K_u S=u_*^2,\quad K_uS^2=\epsilon.
\]

These imply

\[
 e={u_*^2\over\sqrt{ab}},\quad K_u=\kappa u_*d,\quad
 S={u_*\over\kappa d},\quad
 \boxed{\kappa=C(a^3/b)^{1/4}}.
\]

For defaults \(a=0.361,b=1.604\), require \(|\kappa_b-0.4|\le2\times10^{-4}\), accommodating the rounded input coefficient. The actual values are \(\kappa_b=0.40018037285465263\) and \(\kappa_s=0.4680496398124221\). Exact equality to 0.4 is not the contract. The inverse coefficient is 0.9665641451648294. The base default gives 0.11587435822058197.

The value 0.4 is an empirical rigid-wall modelling target, not a conservation law or universal constant for every ocean boundary layer. The surface uses its separate calibration; it should not also be forced to 0.4. User-modified \(a,b,C\) obey the formula, not the default numerical target. No assertion is made at \(d=0\), in the viscous sublayer, in the bottommost cell with an extra TKE boundary sink, or where the stated equilibrium assumptions fail.

**I2 — Primary: consistent wall-distance response across all four length scales.** For a fixed valid state with convection/entrainment disabled, define \(L_j=\min(C^sd_s,C^b_jd_b,L_N)\). For \(C^b_2\ge C^b_1>0\), all other parameters fixed and positive,

\[
 L_2\ge L_1,\quad K_{\psi,2}\ge K_{\psi,1}\ (\psi=u,c,e),
 \quad \ell_{D,2}\ge\ell_{D,1},\quad\epsilon_2\le\epsilon_1.
\]

If both states are strictly bottom-limited and none of the outer caps is active, then

\[
 {K_{\psi,2}\over K_{\psi,1}}={\ell_{D,2}\over\ell_{D,1}}
 ={C^b_2\over C^b_1},\qquad
 {\epsilon_2\over\epsilon_1}={C^b_1\over C^b_2}.
\]

The new/old ratio is \(0.967/0.28\). If the surface or stratification limit is already smaller than both bottom limits, the state is unchanged. These are comparisons at identical fields, not monotonicity claims for evolved solutions. Domain: finite wet-cell data, \(e>e_{\min}\), finite Richardson number or its well-defined zero-shear stable limit, nonnegative diffusivity caps, and positive stability coefficients. Distances and caps must use their actual center/face locations. The column-depth cap is particularly important for TKE because \(C_e^{lo}=7.863\).

**I3 — Supporting: nonnegative, bounded, downgradient transport.** For finite valid wet-cell inputs, nonnegative \(e_{\min}\), nonnegative diffusivity caps and the default positive coefficients, with surface convection disabled,

\[
 0\le K_\psi\le\min(K_{\psi,\max},Hw_f),\qquad
 -J_c\,\partial_z c=K_c(\partial_zc)^2\ge0,
 \qquad P=K_u[(\partial_zu)^2+(\partial_zv)^2]\ge0.
\]

For positive TKE at interior cell centers, \(\epsilon\ge0\). In an impermeable column with no imposed tracer flux, the corresponding semidiscrete diffusion operator satisfies
\(\sum_k V_kc_k\dot c_k=-\sum_f A_f\Delta z_f K_{c,f}(D_zc)_f^2\le0\).
Boundary-face diffusivities are zero by masking. This does not prohibit explicitly imposed boundary fluxes, nor require inactive immersed cells to be finite: they are deliberately marked NaN. It is not a claim that arbitrary time steps, advection, or external forcing preserve tracer variance.

## 3. Sources of the invariants

I1 is independently obtained by eliminating \(S\) and \(e\) from the four balance/constitutive equations above. The implementation supplies \(\ell_D=Cd/b\), fixing the denominator convention unambiguously.

[Wagner et al. (2025), *Formulation and Calibration of CATKE*, §3.1.3 and §5.1.2](https://doi.org/10.1029/2024MS004522) supplies the similarity assumptions and surface interpretation. There is a source discrepancy: the accessible published HTML writes a **product** with \(C_D^0\) in equations (35) and (53), while equation (32) specifies \(\ell_D=C^sd/C_D^0\). The product gives about 0.593 with these defaults, not the reported 0.47. The division used by this PR follows from equation (32) and production–dissipation balance. This appears to be an algebraic/typesetting inconsistency in the paper, not evidence for an alternative implemented convention. Its editorial status is unverified. The paper's calibration does not establish the bottom target independently.

I2 and the bounds in I3 follow from monotonicity of positive multiplication, `min` and `max` in the merged implementation. The flux signs and variance identity in I3 follow from downgradient transport and discrete summation by parts on a closed finite-volume column.

Pinned implementation sources:

- [Mixing lengths and defaults](https://github.com/CliMA/Oceananigans.jl/blob/24de7026893af38df065920fbbf80f05024fd1b9/src/TurbulenceClosures/turbulence_closure_implementations/TKEBasedVerticalDiffusivities/catke_mixing_length.jl).
- [TKE equation and dissipation convention](https://github.com/CliMA/Oceananigans.jl/blob/24de7026893af38df065920fbbf80f05024fd1b9/src/TurbulenceClosures/turbulence_closure_implementations/TKEBasedVerticalDiffusivities/catke_equation.jl).
- [Diffusivity construction](https://github.com/CliMA/Oceananigans.jl/blob/24de7026893af38df065920fbbf80f05024fd1b9/src/TurbulenceClosures/turbulence_closure_implementations/TKEBasedVerticalDiffusivities/catke_vertical_diffusivity.jl) and [TKE stepping, including the bottom sink](https://github.com/CliMA/Oceananigans.jl/blob/24de7026893af38df065920fbbf80f05024fd1b9/src/TurbulenceClosures/turbulence_closure_implementations/TKEBasedVerticalDiffusivities/time_step_catke_equation.jl).

## 4. Operating domain

The following are representative ocean-column check ranges, not enforced API limits or a claim of empirical validation throughout their Cartesian product.

| Input | Range/configuration relevant to these checks |
|---|---|
| Column depth \(H\), distances \(d_b,d_s\) | \(H=10\)–1000 m; \(0\le d_b,d_s\le H\), summing to \(H\) in the flat fixed column. Similarity samples strictly inside the water and away from either boundary stencil. |
| Vertical resolution | 32, 64, 128 cells; uniform and smoothly stretched, with additional refinement for the log profile. |
| TKE \(e\) | \(10^{-8}\)–\(10^{-2}\) m² s⁻²; default floor \(10^{-9}\). Separately test zero and values just below/at/above the floor for I3, not I1/I2 dissipation ratios. |
| Horizontal velocity and shear | Velocity offsets of order 0–1 m s⁻¹; \(S=10^{-5}\)–\(10^{-1}\) s⁻¹, plus zero shear. Rotate the shear between \(u\) and \(v\). |
| Buoyancy | Use \(b(z)=N^2z\); \(N^2=-10^{-5}\)–\(10^{-3}\) s⁻², including exact zero. Only gradients matter. Select compatible \(N^2,S\) to exercise specified \(Ri=N^2/S^2\). |
| Wall stress/friction velocity | \(u_*=10^{-3}\)–\(3\times10^{-2}\) m s⁻¹; kinematic stress \(u_*^2\). Density is unnecessary when specifying kinematic fluxes. |
| Surface buoyancy flux and rotation | Set \(J_b=0\) and no Coriolis force for I1; keep \(J_b=0\) for I2/I3 to exclude convective and entrainment lengths. Default flux regularizer is \(10^{-11}\) m² s⁻³. |
| Dimensionless coefficients | Defaults; compare \(C^b=0.28,0.967\), and the exact inverse value above. For algebraic scaling sample positive \(C^b\) in 0.1–2. Stability transitions are \(Ri=0,0.254,1.274\). |
| Caps, tracer gradients, time step | Default caps are infinity; additionally use finite positive caps below/above the predicted diffusivities and zero caps. Passive tracer amplitudes may be normalized to unity. Time steps must resolve both diffusion and TKE decay as detailed below. |

## 5. Numerical check plan

Use the merged library, not a stand-alone reimplementation of its kernels. Import `CATKEVerticalDiffusivity` and the time-discretization types from `Oceananigans.TurbulenceClosures`, and the configurable `CATKEMixingLength{FT}` from `Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities`. Construct

```julia
grid = RectilinearGrid(CPU(), FT; size=Nz, z=(-H, 0),
                       topology=(Flat, Flat, Bounded))
closure = CATKEVerticalDiffusivity(td, FT;
    mixing_length=CATKEMixingLength{FT}(; Cᵇ=Cb))
model = HydrostaticFreeSurfaceModel(grid; closure,
    buoyancy=BuoyancyTracer(), tracers=(:b, :c),
    momentum_advection=nothing, tracer_advection=nothing)
set!(model; b=z -> N2*z, e=e0, u=z -> S*z, v=0, c=z -> z/H)
```

Do **not** include `:e` in `tracers`: the merged constructor adds it and rejects explicit duplication. `set!` calls the model state update. Obtain diagnostics with `viscosity(model)` and `diffusivity(model, Val(:c))`, `diffusivity(model, Val(:e))` from `Oceananigans.TurbulenceClosures`, then `Array(interior(...))`. Use `nodes` to identify their actual face positions. Model evolution uses `Simulation(model; Δt, stop_iteration=...)` and `run!`. No direct call to private mixing-length or dissipation kernels is necessary.

**I1 checks.** First read \(a,b,C\) from constructed closures for both `Float64` and `Float32`, and both explicit and vertically implicit discretizations. Evaluate the boxed relation, compare the bottom result to 0.4 with absolute tolerance \(2\times10^{-4}\), and compare to the rounded-default oracle with relative tolerances \(10^{-12}\) and \(2\times10^{-6}\), respectively. Check the surface separately and check that an explicit old-coefficient override yields 0.1158743582. This tests parameter delivery but is not sufficient alone.

Next choose three logarithmically spaced \(u_*\) values in the stated range, initialize constant \(e=u_*^2/\sqrt{ab}\), and sample bottom-limited faces with \(d_b/H\le0.01\). Compare the API viscosity to \(\kappa_bu_*d_b\). Mirror the experiment near the surface with \(C^s\). No boundary stress need be applied for this instantaneous constitutive check. For the equilibrium interpretation, initialize \(u=(u_*/\kappa)\log(d/d_0)\) at interior centers and evaluate the discrete shear away from the first/last two cells. Require \(K_u|D_zu|/u_*^2\) within 1% where \(\Delta z/d\le0.05\), and decreasing error on refinement. This is a local manufactured-profile check, not a claim that zero-flux boundary conditions sustain a global steady log layer.

Test the **center dissipation path** independently: initialize neutral, motionless, uniform positive TKE and zero boundary fluxes, then advance one very small step. At centers away from the bottommost cell, initial shear production and TKE diffusion vanish, so
\[
 -[e(\Delta t)-e(0)]/\Delta t\longrightarrow b e_0^{3/2}/L_c.
\]
Choose \(\Delta t\) so that both \(\Delta t\sqrt{e_0}/\ell_D\) and the diffusive Courant number are below \(10^{-3}\); halve it twice. For Float64 require relative error below 0.5% and convergence toward the oracle. Use this diagnostic together with the viscosity check to verify production–dissipation compatibility without accessing private functions. Exclude the bottommost cell's separate sink. Repeat with both time discretizations and `tke_time_step=nothing` versus two TKE substeps.

**I2 checks.** On identical states compare closures with the old/new \(C^b\) and a positive coefficient sweep. Test all three diagnostic fields, and the center dissipation diagnostic above. Construct strict bottom-limited, surface-limited, and stratification-limited samples with at least a 10% margin between competing lengths, as well as points at either side of each limiting-length crossover. Test uncapped, depth-capped, and finite-diffusivity-capped cases; only apply exact ratio assertions to uncapped samples. Compare fields to the formulas in §1 with relative tolerance \(10^{-10}\) (Float64) or \(10^{-4}\) (Float32), plus absolute \(10^{-13}\) or \(10^{-8}\) m² s⁻¹ near zero. The decay-derived dissipation ratios use the 0.5% tolerance above.

For the face branches sample \(Ri=-0.1,0,0.1,0.254,0.764,1.274,2\), and offsets of \(10^{-5}\) on each side of all three thresholds. Include stable zero shear and neutral zero shear. Use uniform TKE and linear buoyancy/velocity profiles so interpolation is controlled; vary \(e_0\) to choose the active length limit. Repeat shear directions \((S,0),(0,S),(S/\sqrt2,S/\sqrt2)\). Repeat on a small horizontally periodic grid with `ImmersedBoundaryGrid(..., GridFittedBottom(...))`, including a raised grid-aligned flat bottom and two different column depths, using active-cell maps both on and off. Compute distances from the represented wet bottom, and ignore dry NaNs. Repeat selected cases with a horizontal array of closures carrying different \(C^b\) values and with CATKE plus a zero-diffusivity closure tuple, to check coefficient selection.

**I3 checks.** Apply finite/nonnegative/bounded assertions to every wet face of those fixtures, and check boundary masks are exactly zero. Include TKE floor cases and each cap value. Compute \(K_c(D_zc)^2\) and \(K_uS^2\) from API fields and gradients, and require nonnegativity to roundoff. For a closed-column dynamics check, initialize a nonconstant smooth passive tracer, zero normal tracer fluxes and no advection; evolve a small step and compare its volume-weighted variance. Choose explicit diffusion Courant number at most 0.1 and repeat with implicit stepping and half the time step. Allow a variance increase no larger than \(100\,\epsilon_{\rm mach}\) times initial variance; a large increase fails. Verify the frozen-coefficient finite-volume variance identity separately with relative residual below \(10^{-10}\) or \(10^{-4}\). Test both `:QuasiAdamsBashforth2` and `:SplitRungeKutta3` by actually passing `timestepper` to the model.

Coverage is deliberately limited to the nonconvective domain of these invariants. Positive surface-flux convection, entrainment, negative-TKE damping, moving vertical coordinates, and arbitrary roughness/bottom-drag laws are not covered by the similarity claim. The PR adds no branch in those regimes. These checks have been specified, not executed against Julia here.

## 6. Do the PR's own tests check these invariants

**No tests are added or modified by this PR.** For I1, I2, and I3 individually, PR-specific coverage is **no**. Existing relevant tests in [test/turbulence_closures/closures.jl at the merge](https://github.com/CliMA/Oceananigans.jl/blob/24de7026893af38df065920fbbf80f05024fd1b9/test/turbulence_closures/closures.jl) provide the following narrower evidence:

| Invariant | Existing coverage | Actual assertions and limitations |
|---|---|---|
| I1 | No | `Time-stepping with CATKE closure` / `CATKE by itself [...]` checks construction, boundary-condition replacement, and successful stepping; `run_catke_tke_substepping_tests` checks Euler/Adams–Bashforth update identities. Neither computes κ, fits a log layer, nor checks the dissipation balance. |
| I2 | No | `CATKE in a 2-tuple with HorizontalScalarDiffusivity [...]` and `CATKE in a 3-tuple [...]` check tuple ordering and stepping identities. `CATKE with ImmersedBoundaryGrid and active_cells_map [...]` runs two steps and asserts model type. None compares old/new coefficients, four length paths, or limiting branches. |
| I3 | Partial, indirect only | `Vertical diffusive CFL diagnostics` asserts finite positive hard-coded explicit-CATKE diffusion diagnostics and infinite diffusion time/zero CFL for implicit CATKE. This exercises diffusivity aggregation, not pointwise signs, bounds, or tracer variance. No relevant variance assertion was found in these CATKE tests. |

`Surface TKE flux with CATKE in closure tuples of length 1 to 6 [...]` checks the imposed surface flux formula and equality to a single-closure TKE result; it does not check wall similarity. Also, `run_time_step_with_catke_tests(arch, closure, timestepper)` accepts a timestepper argument but does not forward it to its model constructor, so its caller's two-symbol loop is not evidence of both integrators being exercised there.

## 7. Expected verdict

- **I1: expected pass under the stated assumptions and rounding tolerance.** Independent scalar arithmetic confirms the new target. The denominator and cubic power follow from the actual dissipation law. Full finite-grid equilibrium checks remain unverified; no claim is made about a complete bottom boundary layer or stronger drag in an evolved OMIP simulation.
- **I2: expected pass for the fixed-state, nonconvective formulas.** The increased coefficient enters both center and face wall minima and therefore all four paths consistently. Immersed-grid, closure-array, and time-derived dissipation checks remain numerically unverified.
- **I3: expected pass in the stated domain with stable time steps.** Positive coefficients and clipping preserve instantaneous bounds. The time-integrated variance check is less certain because it also exercises time stepping and evolving TKE, not just this constant change.

The changed docstring's attribution to the printed paper equation needs the qualification in §3; the numerical correction itself is supported independently.

## 8. Environment

Pin Oceananigans to merge commit `24de7026893af38df065920fbbf80f05024fd1b9`; the checked `Project.toml` reports version **0.113.5** and Julia compatibility **1.10**. The base commit is `aca5970fa8ed41d5ef119b98b33b67bbcbb6c5da`. A base installation is optional because an explicit `Cᵇ=0.28` on the merge isolates the changed parameter.

Use a fresh Julia project, add Oceananigans from the repository at that revision, instantiate its declared dependencies, and retain the resulting manifest. `Test`, `Random` with a fixed seed, and ordinary array/reduction facilities suffice for these CPU column checks. The package resolves dependencies including KernelAbstractions, OffsetArrays, and its grid/solver dependencies; no LES dataset, plotting package, GPU, or remote compute node is required. Run Float64 first and Float32 as a separate precision check. Julia/runtime availability and package installation have not been tested in this review.

## 9. Formalization feasibility

**I1:** readily expressible over positive real numbers. Introduce positive fourth roots or prove \(\kappa^4=C^4a^3/b\), derive it algebraically from stress and production–dissipation balance, and establish uniqueness of the positive root. Rational constants plus interval bounds can certify the 0.4 tolerance. This proves conditional model consistency, not the empirical universality of 0.4, the existence of a turbulent similarity layer, or its applicability to waves/rough walls.

**I2:** suitable for a short Lean proof over transcribed scalar `min`/`max` functions, with positive-coefficient hypotheses and explicit cases for the selected limits. Exact ratio identities require strict active-branch and nonzero-denominator hypotheses. Extending the proof to grids requires formalizing staggered distances, interpolation, wet masks, and parameter dispatch.

**I3:** positivity and bounds are elementary real inequalities; the finite-volume variance identity can be proved by discrete summation by parts with zero boundary fluxes. A fully discrete theorem additionally requires explicit stability conditions or properties of the implicit solve. IEEE floating-point rounding, NaN/Inf sentinels, eager `ifelse` evaluation, TKE floors, GPU kernels, and Julia dispatch are substantial obstacles to proving the compiled implementation. A proof of a real-valued transcription alone would not certify those implementation details.
