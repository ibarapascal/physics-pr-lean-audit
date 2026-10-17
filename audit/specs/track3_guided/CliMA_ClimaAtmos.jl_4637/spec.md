## 1. What changed

PR #4637, “Fixes problems near upper boundary and tracer-mass continuity,” changes ClimaAtmos.jl from base `2d43bf4c1e40e4fff6085dcb0de6ace627375094` to merge `a8efce5a0040924094a58739db91a0f3162f779c`. This specification uses the supplied diff and source fetched at that merge; commit-message experiments are claims, not independently reproduced results.

- `src/cache/precomputed_quantities.jl`: new `top_velocity`, revised `set_velocity_at_top!(Y, ᶠuₕ³, turbconv_model)`, and its caller `set_implicit_precomputed_quantities!` replace the top condition `u₃ = 0` with `u₃ = -uₕ³/g³³`. Since `u³ = uₕ³ + g³³u₃`, this cancels flow normal to a coordinate surface even on a warped grid. Each prognostic EDMFX draft receives the same boundary covariant velocity. Prescribed flow bypasses this boundary-state adjustment.
- `src/prognostic_equations/implicit/implicit_tendency.jl`: `implicit_vertical_advection_tendency!` now uses `ᶜadvdivᵥ`, whose bottom and top fluxes are zero, for continuity, instead of unrestricted `ᶜdivᵥ`. `src/utils/abbreviations.jl` updates the operator's comment; its definition was already zero-flux.
- `src/prognostic_equations/advection.jl`: `explicit_vertical_advection_tendency!` always applies grid-mean vertical transport to discovered grid-scale tracers other than energy and total water. Previously it skipped this for prognostic EDMFX with SGS mass flux enabled.
- `src/prognostic_equations/edmfx_sgs_flux.jl`: `edmfx_sgs_mass_flux_tendency!` replaces absolute draft/environment transport of discovered tracers with deviation transport, using `edmfx_sgsflux_upwinding` instead of `edmfx_tracer_upwinding`. It also adds each draft and environment total-water tendency to moist-air density.

For precision, let \(I\) denote center-to-face interpolation, \(D_0\) the divergence with zero boundary flux, \(J_c,J_f\) geometry Jacobians, and \(R_s(v,b)\) the chosen velocity-weighted reconstruction, with \(R_{\mathrm{none}}(v,b)=vI(b)\). The actual discretization is

\[
 T_s(r,v,b)=-D_0\!\left[\frac{I(rJ_c)}{J_f}R_s(v,b)\right],\qquad
 A_\rho=-D_0\!\left[\frac{I(\rho J_c)}{J_f}u^3\right].
\]

With \(Q_\chi=\rho\chi\), \(a_k=(\rho a)_k/\rho_k\), and \(k=0\) for environment, the revised correction is

\[
 S_\chi=\sum_{k=0}^n T_s(\rho_k,u_k^3-u^3,a_k(\chi_k-Q_\chi/\rho)),
 \qquad S_\rho=S_{q_{tot}}.
\]

This is a discrete reconstruction of a covariance flux; it is not generally interpolation of the single product \(\rho_k a_k(u_k^3-u^3)(\chi_k-\chi)\).

- `src/prognostic_equations/implicit/manual_sparse_jacobian.jl`: `update_jacobian!` replaces absolute upwind tracer blocks, including environment contributions, with a central, first-draft-only approximation for `q_lcl`, `q_icl`, `q_rai`, `q_sno`, `n_lcl`, `n_rai`, conditional on field existence and the 1M/2M model types. Define
  \[
  B=dt\gamma D_0\operatorname{diag}\!\left(\frac{I(\rho_1J_c)}{J_f}(u_1^3-u^3)\right)I\operatorname{diag}(a_1).
  \]
  The new SGS residual blocks are \(W_{Q_\chi,\chi_1}=-B\), an addition \(B\operatorname{diag}(1/\rho)\) to \(W_{Q_\chi,Q_\chi}\), and
  \[
  W_{Q_\chi,u_3}=dt\gamma D_0\operatorname{diag}\!\left(\frac{I((\chi_1-\chi)\rho_1J_ca_1)}{J_f}g^{33}\right).
  \]
  The convention is a tendency derivative scaled by \(dt\gamma\), minus identity on diagonal residual blocks. Environment derivatives are deliberately omitted; comments motivate this by small-draft-area scaling. Two obsolete scratch matrices are removed in `src/cache/temporary_quantities.jl`.
- `src/prognostic_equations/remaining_tendency.jl` only adds a TODO about the existing condensate Rayleigh sponge; its behavior is not fixed. `NEWS.md`, the reproducibility counter (371 to 372), and regression tests are updated.

## 2. Invariants

**I1 — Primary: impermeability in the coordinate-normal direction.** At bottom and top, after boundary-state preparation,
\[
 u^3=u_h^3+g^{33}u_3=0,\qquad u_j^3=0\quad(j=1,\ldots,n).
\]
Domain: finite velocities, nonsingular static geometry with \(g^{33}>0\), no prescribed-flow override, and drafts sharing the grid-mean horizontal velocity/projection. The physical condition is zero normal velocity, not zero Cartesian vertical velocity or zero covariant component. Flat grids are the special case \(u_h^3=0\).

**I2 — Primary: grid-mean constant-tracer consistency.** For each transported mass-specific scalar and constant \(C\),
\[
 Q_\chi=C\rho\implies A_\chi=C A_\rho.
\]
Domain: isolated grid-mean vertical advection, positive density, matching velocity and geometry, and constant-preserving reconstruction. Include total water in its implicit pathway and discovered species/number tracers in their explicit pathway. Energy is excluded: the energy equation transports total enthalpy, not a passive specific energy. Sedimentation, conversion, forcing, surface fluxes and SGS corrections are separate operators.

**I3 — Primary: null SGS correction for a uniform scalar.**
\[
 Q_\chi=C\rho,\quad \chi_j=C\ (\forall j),\quad \chi_0=C
 \implies S_\chi=0.
\]
Domain: the difference-form SGS pathways, finite positive subdomain thermodynamic densities and admissible area masses. Environment reconstruction must preserve a constant, including its small-area fallback. The property also holds if every subdomain velocity equals the mean, irrespective of scalar differences. Disabled SGS flux and the generic non-prognostic fallback must contribute zero.

**I4 — Primary: water/air mass bookkeeping.** For SGS mass flux alone,
\[
 S_\rho=S_{q_{tot}},\qquad \partial_t(\rho-\rho q_{tot})\big|_{SGS}=0.
\]
Domain: moist prognostic EDMFX with SGS mass flux active; both increments must be compared from identical initial tendency accumulators. This encodes the library convention that turbulent water redistribution contributes moist-air mass while leaving dry-air mass unchanged under this operator. It is a modeling convention, not a universal statement about every barycentric turbulence formulation.

I2 and I3 do **not** jointly assert arbitrary passive-tracer mixing-ratio preservation under the entire new moist system. If \(S_\rho\ne0\) from nonuniform water but \(S_\chi=0\) for a uniform passive tracer, then \(\partial_t\chi|_{SGS}=-\chi S_\rho/\rho\). Whether passive species should accompany this water-associated mass flux requires an additional modeling choice. Likewise, defining vapor by \(q_v=q_{tot}-\sum q_{cond}\) gives a bookkeeping identity, not positivity or phase-equilibrium preservation.

**I5 — closed-column transport conservation.** For cell volume/quadrature weights \(V_i>0\), each isolated conservative vertical tendency satisfies
\[
 \sum_i V_i A_{\rho,i}=0,\quad \sum_i V_i A_{\chi,i}=0,
 \quad \sum_i V_i S_{\chi,i}=0,\quad \sum_i V_i S_{\rho,i}=0.
\]
Domain: matching conservative discrete divergence/weights and zero top/bottom flux. This applies to each column; global extensions require periodic/closed horizontal boundaries. Precipitation outflow and imposed surface fluxes are excluded and would require explicit boundary budgets.

**I6 — consistency of retained linearization.** With geometry, thermodynamic densities and areas frozen, the manual blocks intended to represent central draft transport should equal \(dt\gamma\) times the derivative of that retained transport, with the residual identity removed. In particular,
\[
 W_{Q,\chi_1}\delta c+W^{SGS}_{Q,Q}(\rho\delta c)=0
\]
for a common perturbation of mean and draft scalar at fixed density. The isolated advective density/velocity block must differentiate \(A_\rho\), including boundary rows. Domain: central reconstruction, a single retained draft, smooth states away from clipping/regularization transitions. Exact equality with the full nonlinear Jacobian is **not** required: environment terms, noncentral reconstruction and other couplings are omitted. A specific concern is that the implemented velocity block interpolates a product, whereas differentiating the actual central flux gives the product \(I(\rho_1J_c)I(a_1(\chi_1-\chi))\); these need not agree for vertically varying fields.

## 3. Sources of the invariants

I1 follows from the tensor identity \(u^3=g^{3i}u_i\) and no penetration through a surface of constant vertical coordinate. The implementation convention is documented directly in [boundary-state preparation](https://raw.githubusercontent.com/CliMA/ClimaAtmos.jl/a8efce5a0040924094a58739db91a0f3162f779c/src/cache/precomputed_quantities.jl).

I2 follows by substituting a constant concentration into the conservative continuity and scalar equations. I5 follows from the discrete divergence theorem: interior face fluxes cancel when weighted by cell volumes. Their particular discrete form is established by [vertical transport and continuity](https://raw.githubusercontent.com/CliMA/ClimaAtmos.jl/a8efce5a0040924094a58739db91a0f3162f779c/src/prognostic_equations/implicit/implicit_tendency.jl) and [operator boundary definitions](https://raw.githubusercontent.com/CliMA/ClimaAtmos.jl/a8efce5a0040924094a58739db91a0f3162f779c/src/utils/abbreviations.jl).

I3 and I4 follow algebraically from the [SGS flux implementation](https://raw.githubusercontent.com/CliMA/ClimaAtmos.jl/a8efce5a0040924094a58739db91a0f3162f779c/src/prognostic_equations/edmfx_sgs_flux.jl). The environment's regularized mean and its deliberate departure from exact domain decomposition at small area are documented in [variable manipulations](https://raw.githubusercontent.com/CliMA/ClimaAtmos.jl/a8efce5a0040924094a58739db91a0f3162f779c/src/utils/variable_manipulations.jl).

I6 is a differentiation/linear-algebra consistency requirement, qualified by the explicit approximation in [the manual Jacobian](https://raw.githubusercontent.com/CliMA/ClimaAtmos.jl/a8efce5a0040924094a58739db91a0f3162f779c/src/prognostic_equations/implicit/manual_sparse_jacobian.jl). No external standard or cited paper is necessary to derive these identities. The source does not establish universal physical correctness of the moist-air-mass convention or the untreated sponge.

## 4. Operating domain

The following are representative atmospheric verification ranges, not library-enforced bounds:

| Variable | Common range and required constraints |
|---|---|
| Height and grid spacing | Columns 0–1.5 km; mountain configurations to 21–30 km; vertical spacing 10–1000 m. Positive cell volumes; test uniform and stretched spacing. |
| Pressure, temperature, density | Approximately 1–105 kPa, 180–320 K, 0.01–1.5 kg m⁻³. Derive consistent states with the library thermodynamics rather than sample these independently. Density and subdomain densities must be strictly positive. |
| Velocity | Horizontal components −100 to 100 m s⁻¹; vertical/draft departures −20 to 20 m s⁻¹, including both signs and exact zero. Coordinate components depend on geometry and must be constructed from physical vectors. |
| Moisture | \(0\le q_{tot}\lesssim0.03\) kg kg⁻¹; each condensate typically 0–0.01, with their sum at most total water. Use zero, trace and nonzero values. \(C=1\) is an algebraic stress test, not a physical state with all phases equal to one. |
| Number tracers | Nonnegative counts per unit air mass; sample 0 and \(10^3\)–\(10^9\) kg⁻¹ where the selected 2M state supports them. No unit-interval constraint applies. |
| Area masses | \(m_j=(\rho a)_j\ge0\), \(m_0=\rho-\sum m_j\ge0\); typical draft area fractions 0.001–0.2. Include zero drafts and tiny environment mass as numerical edge cases, not ordinary conditions. |
| Geometry | Flat, Schar/Linear with nonzero boundary projection, and SLEVE with a flattened lid; positive \(J_c,J_f,g^{33}\), resolved slopes. Do not independently randomize metric entries. |
| Time step | 0.1–120 s. For explicit transport choose a vertical Courant number at most 0.5; 120 s is the PR's implicit column spin-up case, not a universal stable step. |

Use Float64 and Float32. Reject NaN/infinite inputs. Sample regularization near both \(m_0=\epsilon(FT)\) and the configured `a_half` transition separately; finite-difference derivative tests must not straddle either branch.

## 5. Numerical check plan

Construct configurations through `CA.AtmosConfig(Dict(...); job_id=...)` and `CA.get_simulation(config)`, then obtain `simulation.integrator.u`, `.p`, `.t`. Advance realistic states with `ClimaTimeSteppers.step!`. This is the same public construction route used by `generate_test_simulation`; no reimplementation of the dynamical core is needed. Operator isolation and matrix inspection below use module-qualified functions available in this revision and used by its tests; they are not promised stable external APIs. A public time-stepping-only test cannot isolate all these invariants.

Start with the PR's two fixtures: (a) plane, `ConstantBuoyancyFrequencyProfile`, Schar, Linear, 100 km, 8 horizontal elements, 21 km, 20 vertical elements, 500 m bottom spacing, 1 s step; (b) `DYCOMS_RF01` column, prognostic EDMFX, 1M, 1500 m, 30 levels, ARS222, 120 s step, and `toml/prognostic_edmfx.toml`, spun up for ten steps. Disable physical sources for isolated transport measurements; do not confuse disabling implicit microphysics with disabling sedimentation.

Coverage matrix:

- Grid geometry: flat column; flat plane; Schar with Linear and SLEVE; a small 3-D warped box to exercise both horizontal vector components. Check bottom and top separately.
- Turbulence: no turbulence, `EDOnlyEDMFX` fallback where constructible, prognostic EDMFX with SGS mass flux off/on; `updraft_number=1,2`, inspecting every draft. Multi-draft tendency tests are separate from the first-draft manual approximation.
- Microphysics: `dry`, `0M`, `1M`, `2M`, and `2MP3` when compatible with the chosen turbulence configuration. Enumerate actual fields using `CA.sgs_tracer_names(Y)` and the grid tracer traversal; do not silently equate 2MP3 with the changed 1M/2M Jacobian branch. Check all present liquid, ice, rain, snow mass fields and both number fields. For a supported custom passive-tracer extension, test one additional discovered field; absence of a public registration route must be reported as a coverage gap.
- Sweep `tracer_upwinding` and `energy_q_tot_upwinding` through `none`, `first_order`, `third_order`, `vanleer_limiter`; SGS reconstruction through its documented `none`, `first_order`, `third_order`. Independently vary draft-internal `edmfx_tracer_upwinding` over documented `none`, `first_order`, and exercise the explicit third-order matrix branch if accepted. An unsupported configuration is recorded, not counted as a pass.

Use deterministic constant, sinusoidal and sharp-but-admissible profiles, plus 20 seeded realizations with spatially varying density and draft area. Ensure at least one nonzero interior tendency or face flux per positive case. Rebuild precomputed fields for physical states. For unphysical constant-species probes, intentionally freeze a consistent physical cache, as the PR tests do; never claim such a probe tests thermodynamic admissibility.

**I1:** Call `CA.set_implicit_precomputed_quantities!(Y,p,t)` after assigning both signs of wind. Inspect `p.precomputed.ᶠu³` and every `ᶠu³ʲs` with `Fields.level` at `CA.half` and `Spaces.nlevels(axes(Y.c))+CA.half`. Require residual at most \(100\epsilon\max(|u_h^3|,|g^{33}u_3|)\) per boundary, with exact zero for an identically zero scale. Require nonzero horizontal projection in the warped regression. A prescribed-flow configuration is a negative control: state clamping is intentionally bypassed.

**I2:** Zero separate tendency buffers. Obtain the advective density reference from the dry/no-sedimentation `CA.implicit_vertical_advection_tendency!`, or the geometry-aware \(A_\rho\) above. Call `CA.explicit_vertical_advection_tendency!` for each discovered passive field and isolate total-water advection using `CA.vertical_transport` with the configured reconstruction. Test \(C=0,1,10^{-6},0.01\), and number-scale constants separately. Require 
\(\|A_\chi-CA_\rho\|_\infty\le200\epsilon S\), where \(S\) is the larger of the reference tendency and the sum of absolute face-flux contributions divided by cell volume, scaled by 
\(|C|\). The latter avoids cancellation-dependent tolerance. Test even with artificial nonzero cached boundary velocity: \(D_0\) must still enforce zero flux structurally.

**I3:** Call `CA.edmfx_sgs_mass_flux_tendency!` on zero buffers with every draft and the mean set to each constant. Verify the reconstructed environment constant separately. Test unequal, nonzero draft velocities and varying areas; then test equal velocities with nonuniform scalars. Require zero within \(500\epsilon S_{SGS}\), using the corresponding undeviated subdomain flux-divergence scale and an environment subtraction condition factor \(\max(1,\rho/m_0)\) when above the fallback threshold. Report this factor; very ill-conditioned cases cannot establish tight accuracy. In disabled/fallback cases require unchanged buffers exactly.

**I4:** Use realistic nonuniform humidity and nonzero SGS water flux. From zero buffers require bitwise equality of `Yₜ.c.ρ` and `Yₜ.c.ρq_tot`; repeat with identical nonzero starting buffers. Test dry and disabled cases for absence of this addition. Do not demand equality of complete tendencies containing unrelated physics.

**I5:** Integrate each isolated tendency using ClimaCore field integration/volume weights, never an unweighted array sum. Require 
\(|\sum V_i T_i|\le 500N\epsilon\sum V_i|T_i|\) for \(N\) vertical levels (exact zero when both sides vanish). Repeat with 20, 40 and 80 levels and stretched meshes, for all reconstructions. A short closed dry run provides a secondary time-step-level mass budget; phase-changing moist runs require accounting for precipitation and other sources and are not the primary test.

**I6:** Build `CA.Jacobian(CA.ManualSparseJacobian(),Y,p.atmos)`, call `CA.update_jacobian!(jac,Y,p,dtγ,t)`, and inspect named blocks with `CA.first_column_block_arrays`. Subtract identity and unrelated diffusion contributions consistently. Compare directional differences of the retained central draft operator with steps \(h/s=10^{-4},10^{-5},10^{-6}\) in Float64 (use nondimensional perturbation scales \(s\)); accept relative error \(10^{-6}\) when the difference converges. Check scalar cancellation to \(500\epsilon\) times the block-action scale. Perturb interior and boundary velocities separately for the density row, with SGS off. For the SGS velocity block compare constant and varying \(\rho_1J_c\) and \(a_1(\chi_1-\chi)\): this explicitly tests the interpolation discrepancy. Test diffusion explicit/implicit and all six present 1M/2M fields. Full-Jacobian comparisons at noncentral options or finite environment area are diagnostic only; an exact-match assertion there would contradict the documented approximation.

## 6. Do the PR's own tests check these invariants

| Invariant | Assessment | Added test and actual assertion |
|---|---|---|
| I1 | Partial | `Impenetrable boundary condition over topography` in `test/prognostic_equations/advection_tests.jl`: Float64 Schar/Linear, nonzero projection at both boundaries, grid-mean normal velocity within \(100\epsilon\) times projection. No draft, 3-D or Float32 check. |
| I2 | Partial | `Grid-mean transport consistency over topography`: assigns all five water mass fractions one with the existing cache; full implicit total-water tendency equals density tendency, including matched sedimentation; explicit four condensate advections equal a separately assembled continuity flux divergence. Requires nonzero references and \(100\epsilon\)-scaled tolerances. Only 1M/default reconstruction, without active EDMFX. |
| I3 | Partial | `EDMFX SGS mass-flux consistency`: after ten DYCOMS steps sets total water and four condensate fractions to one in mean and draft 1, retaining physical cached densities/velocities. All five tracer tendencies and density must be within \(10^{-10}\) times the original total-water flux scale. No number species, multi-draft, option sweep or small-environment case. |
| I4 | Yes for the tested branch | The same EDMFX test first requires nonzero SGS water tendency and bitwise equality of density and total-water tendency arrays. This directly checks the identity for the default moist one-draft 1M path; broader model coverage is absent. |
| I5 | No | Neither new test integrates column/global mass or tracer budgets. Commit-message mass-drift measurements are not test assertions. |
| I6 | No | No added finite-difference, matrix-action or Jacobian null-direction test. Successful spin-up is only an indirect solver exercise. |

The new consistency file is registered under `Tracer/mass transport consistency` in `test/runtests.jl`; both regression files run in the `dynamics`/`all` groups.

## 7. Expected verdict

- **I1: expected pass** for the stated domain, including shared-horizontal-velocity drafts. Division/multiplication cancellation permits roundoff rather than guaranteed bitwise zero.
- **I2: expected pass** for isolated advection with constant-preserving reconstructions. Branch-wide verification remains unperformed.
- **I3: expected pass** in ordinary admissible states. Tiny environment mass is **uncertain** numerically because environment subtraction can magnify roundoff, even though the exact-arithmetic constant identity and fallback are sound.
- **I4: expected pass**: the same computed increment is accumulated into both variables in the same order. This does not resolve the broader passive-tracer/modeling question in Section 2.
- **I5: expected pass** for isolated flux-divergence operators with the correct geometry weights. This is not a claim of closed mass budgets for the full atmosphere with precipitation or surface exchanges.
- **I6: mixed**. Scalar cancellation and the isolated advective density/velocity block should pass. The claimed exact retained central SGS velocity derivative is **expected to fail for general varying profiles**, because interpolation of a product differs from the product of interpolants. This is a source-derived expectation, not a measured failure. Full SGS density-row exactness is also not established by this PR: adding water flux to density introduces couplings beyond the advective density row. The manual matrix is expressly approximate, so these limitations do not alone demonstrate an incorrect time-integrated physical solution.

No numerical checks were executed for this specification. The unchanged condensate sponge precludes extending these conclusions to full mass/energy/phase consistency in sponge-enabled runs.

## 8. Environment

Use the exact merge commit, not merely package version `0.41.1` (the `Project.toml` version did not uniquely identify this change). The project declares Julia compatibility `1.9`; use a compatible Julia release and instantiate the merge's project/test environment, retaining the resolved manifest for repeatability.

Relevant declared compatibility ranges include ClimaCore `0.14.54`, ClimaComms `0.6.9`, ClimaTimeSteppers `0.9.0`, Thermodynamics `1.2`, CloudMicrophysics `0.36`, ClimaParams `1.1`, and AtmosphericProfilesLibrary `0.1.8`; these are compatibility entries, not a verified exact dependency lock. Test support includes Test, SafeTestsets and IntervalSets. Use CPU singleton communication first; GPU/distributed verification is optional additional coverage. Full source/configuration/TOML files and any initialization artifacts must be available. The small column and mountain checks do not require a production-scale run. Base-commit runs are useful regression controls but require adapting the changed boundary-helper signature.

## 9. Formalization feasibility

I1–I5 can be stated and proved over finite real-valued arrays and transcribed operators in Lean. I1 is cancellation with \(g^{33}\ne0\). I2 needs constant-preserving reconstruction and shared mass flux. I3 needs a lemma for constant preservation of environment reconstruction, including its fallback. I4 is equality of two updates. I5 is telescoping with geometry-weighted divergence and zero boundary terms.

I6's scalar blocks admit a short matrix proof. A full equality proof for the velocity block would require the generally false interpolation-product identity; instead one can formalize a counterexample or prove equality under additional constant-coefficient hypotheses. Small-area environment error estimates require explicit assumptions about regularization and which variables are differentiated; the source's order-of-area argument alone is not a proof for every derivative block.

Main obstacles are faithful transcription of ClimaCore staggered geometry, boundary operators and reconstruction branches; Julia generated field traversal and mutation; floating-point cancellation and division bounds; and nonlinear limiter/regularization nondifferentiability. Real-arithmetic proofs would establish structural invariants, not the numerical tolerances or correctness of the complete Julia solver. Proving end-to-end time integration additionally requires the IMEX method, inexact nonlinear/linear solves, state filtering and all active physical source terms.
