## 1. What changed

PR #4822, “Compute horizontal gradient along isentropes in SGS variance,” changes ClimaAtmos from base `e853363fa93882900b87e9572d8971f677246b7e` to merge `ba2bf1225fc48058d709853c065b1db21a412204`. This specification is based on the supplied diff and source files retrieved at the merge commit; no numerical execution of ClimaAtmos was performed.

In `src/parameterized_tendencies/microphysics/cloud_fraction.jl`, three new scalar functions, `isentropic_slope_ratio`, `isentropic_gradient_invariant`, and `isentropic_cell_cap`, feed the modified `set_covariance_cache!`.

Let θ = θ_li, q = q_tot, a = ∂zθ, b = ∂zq, X = Δx_h, Z = Δz, d = θz_min, D = max(a,d), and r = b/D. Define the *discrete* horizontal invariants

\[
 A=F(|\nabla_h\theta|^2),\quad B=F(\nabla_h\theta\cdot\nabla_h q),\quad C=F(|\nabla_hq|^2),
\]

where F is the existing lumped GLL{2} restriction followed by weighted DSS. `hgrad_invariant!` and `hgrad_cross_invariant!` compute these; their algorithms are unchanged. The new code forms

\[
 G=\max(C-2rB+r^2A,0),\qquad
 s=\sqrt{\max(A,0)}/D,\qquad f=(1+(sX/Z)^2)^{-1}.
\]

With h = c_g(c_Δx X)² and v = c_g(c_Δz Z)², the former horizontal covariance block in (θ,q) was h[[A,B],[B,C]]. It is now

\[
 \Sigma_h=\begin{pmatrix}0&0\\0&hfG\end{pmatrix}.
\]

The former vertical block v[[a²,ab],[ab,b²]] becomes

\[
 \Sigma_v=v\begin{pmatrix}a^2&fab\\fab&fb^2\end{pmatrix}.
\]

These vertical scalar expressions apply directly on flat grids; on general geometries the code uses `Geometry.WVector` inner products for that block, and `projected_vector_data(C3, ..., local_geometry)` for r and f. Their geometric consistency must be tested rather than assumed from flat-grid algebra.

If `c_Δx == 0`, horizontal computation is skipped and f is explicitly set to 1. If `c_Δz == 0`, the vertical block is skipped. The turbulent block is unchanged: with k = 2 C_cov L_mix², its diagonal entries are ka² and kb² on a flat grid, and its cross term is ρ_turb√(ka² kb²). Only the `diagnosed` correlation mode accumulates the cross terms. The θ-to-T transformation remains diag(J,1), J = ∂T/∂θ_li. Finally, q variance is capped at (r_max q)²; the diagnosed cross covariance is scaled with the moisture standard deviation before correlation is diagnosed. `constant` mode prescribes the sampled correlation instead.

`src/parameters/Parameters.jl` and `src/parameters/create_parameters.jl` add `sgs_variance_isentropic_min_dtheta_dz`, provisionally defaulting to 10⁻³ K m⁻¹ and overridable through parameter TOML. Although its docstring says positivity is checked in `get_sgs_quadrature`, the merged implementation in `src/config/model_getters.jl` contains no such check. Positive d is therefore a precondition, not an established configuration guarantee.

`docs/src/prophet_closures.md`, `NEWS.md`, and config help explain the new closure, but the documentation's displayed vertical covariance formula omits f. The cache docstring also retains an obsolete second horizontal-addition step; executable code adds it only once. The parameter test expects c_Δx = 1, whereas the pre-PR prose described an off-by-default term. Use actual parameter values, not that prose, in checks. Reproducibility reference counter changes from 426 to 427. These changes are visible in the [merged closure source](https://github.com/CliMA/ClimaAtmos.jl/blob/ba2bf1225fc48058d709853c065b1db21a412204/src/parameterized_tendencies/microphysics/cloud_fraction.jl) and [parameter construction](https://github.com/CliMA/ClimaAtmos.jl/blob/ba2bf1225fc48058d709853c065b1db21a412204/src/parameters/create_parameters.jl).

## 2. Invariants

**I1 — Primary: isentropic contrast and horizontal covariance allocation.** For differentiable fields with a ≥ d > 0, the horizontal-coordinate derivative of q restricted to a local θ surface is

\[
 \nabla_{h|\theta}q=\nabla_hq-(b/a)\nabla_h\theta.
\]

Consequently G must equal its squared norm for constant gradients, remain nonnegative, and vanish when q = αθ + β with constant α,β. It is invariant under orthogonal changes of the physical horizontal basis, including reflections. The horizontal block contributes zero θ variance and zero θ–q covariance. These are block-level statements, not a claim that switching on the block leaves the total covariance unchanged: f also changes the vertical block. Likewise, constant-correlation sampling can change the *sampled* cross covariance when q variance changes.

For filtered spatially varying fields the precise discrete requirement is G = C − 2rB + r²A (up to its roundoff floor), with local r applied **after** F. This is not generally F(|∇q − r∇θ|²) for spatially varying r. Exact cancellation for globally affine q(θ) requires that the vertical derivatives use consistent stencils and that a ≥ d wherever the check is made. When a < d, the code deliberately substitutes a fictitious minimally stable surface; exact physical isentropic cancellation is not required. For example, q = αθ + β then gives G = α²A(1 − a/d)² for constant gradients.

**I2 — Regularized slope and bounded cell excursion.** For finite a,b,A with A ≥ 0, unrestricted signs of a,b, d > 0 and X,Z > 0,

\[
 |r|\le |b|/d,\quad \operatorname{sign}(r)=\operatorname{sign}(b),\quad
 0<f\le1,\quad \ell^2=fX^2\le X^2,\quad s^2\ell^2\le Z^2,
\]

\[
 Z^2-s^2\ell^2=fZ^2.
\]

At A = 0, f = 1. At a ≤ d, r and f are independent of a. For x = sX/Z → ∞, f x² → 1. These are real-arithmetic properties; finite floating-point inputs must also avoid intermediate overflow/underflow to require strict f > 0.

**I3 — Primary: realizable covariance.** For c_g,C_cov ≥ 0, finite L_mix, |ρ_turb| ≤ 1, d,X,Z > 0 and physically admissible thermodynamics, the summed covariance must satisfy

\[
 V_\theta\ge0,\quad V_q\ge0,\quad C_{\theta q}^2\le V_\theta V_q.
\]

In particular det(Σ_v) = v²a²b²f(1−f) ≥ 0. The horizontal block is positive semidefinite even in the floored regime. A finite Jacobian J preserves positive semidefiniteness under diag(J,1). The moisture cap must give 0 ≤ V_q,out ≤ (r_max q)² for q ≥ 0, r_max ≥ 0. The sampled covariance must obey |C_Tq,sampled| ≤ √(V_T V_q,out) in both correlation modes with a prescribed coefficient in [−1,1]. Diagnosed correlation is invariant under a nonzero moisture-cap rescaling only when the code's epsilon floors and fallback threshold are inactive. Where V_T V_q,out ≤ eps(FT)², the specified result is the prescribed fallback, not the limiting ratio.

This is covariance realizability, not a guarantee of nonnegative support at every quadrature node or conservation of moisture by every downstream tendency.

**I4 — Vertical partition and limiting-case consistency.** For flat-grid gradients, the exact f-scaled vertical block above must be assembled. Disabling the horizontal block restores the full vertical block. With flat θ surfaces (A = 0), the original vertical uniform-box variance must be recovered; for c_g = 1/12 and c_Δz = 1 this is the covariance of linear profiles sampled uniformly over a layer of depth Z. Setting both scale factors and C_cov to zero gives zero variances and diagnosed fallback. As a useful partition check, when ∇h q = 0, c_Δx = c_Δz = 1, and C_cov = 0,

\[
 V_{q,h}+V_{q,v}=c_g b^2 Z^2
\]

before the moisture cap, for all s. This follows from I2. It is not a universal identity for unequal scale factors or nonzero horizontal q gradients. Scaling vertical q variance and covariance by f while leaving θ variance unchanged is a modelling choice; the code explicitly notes its asymmetric treatment of θ and q.

**I5 — Conditional modelling-consistency criterion: continuous removal of an amplitude-scaled term.** If c_Δx is intended as a continuous horizontal amplitude parameter, then at fixed inputs and c_Δz > 0,

\[
 \lim_{c_{\Delta x}\to0^+}\Sigma(c_{\Delta x})=\Sigma(0).
\]

This is not required if any nonzero coefficient intentionally selects a different vertical-partition model. The documentation describes a length scale factor and an off switch but does not settle this discontinuous-model convention. Report this criterion separately from I1–I4 rather than presenting it as a conservation law.

## 3. Sources of the invariants

I1 follows from the implicit-function chain rule: dθ = ∇hθ·δx + aδz = 0 gives δz = −∇hθ·δx/a, and substitution into dq gives the stated gradient. Orthogonal invariance follows from preservation of inner products. Conservation of θ_li under the represented displacement is an assumed quasi-adiabatic parcel model; it is not exact under arbitrary precipitation, diabatic heating, or irreversible mixing. The PR cites no external standard fixing the floor, smoothing, or asymmetric partition.

I2 follows algebraically from f = 1/(1+x²). The geometric length is a smooth bound, not literally min(X,Z/s), and the “along-surface length” here is a horizontal projection, not full arclength.

I3 follows from the defining covariance identity uᵀΣu = Var(uᵀZ) ≥ 0 and closure under sums and linear transformations. I4's 1/12 derives from Var(U) = Z²/12 for U uniform on [−Z/2,Z/2]; the f partition is a library convention derived from I2, not the full covariance of arbitrary linear fields over a rectangular cell. I5 is continuity of an amplitude-controlled closure, conditional on the intended parameter semantics.

The [closure documentation](https://github.com/CliMA/ClimaAtmos.jl/blob/ba2bf1225fc48058d709853c065b1db21a412204/docs/src/prophet_closures.md) provides the modelling rationale. The equations in section 1 resolve its omissions using executable source. The [diagnostic implementation](https://github.com/CliMA/ClimaAtmos.jl/blob/ba2bf1225fc48058d709853c065b1db21a412204/src/diagnostics/core_diagnostics.jl) establishes that reported covariance is reconstructed from the sampled correlation and final variances, rather than exposing the raw accumulated cross term.

## 4. Operating domain

The following are proposed atmospheric test ranges, not documented hard validity limits. Sample coherent profiles so that gradients, cell sizes, and mean values do not imply negative water or impossible temperatures.

| Quantity | Common range / required edge cases |
|---|---|
| T, θ_li | T = 200–310 K; θ_li approximately 250–400 K; include warm liquid, mixed-phase, and cold ice states |
| Density / pressure | ρ = 0.1–1.3 kg m⁻³; pressure approximately 10–105 kPa, related by thermodynamics rather than independently assigned |
| q_tot; liquid and ice fractions | q = 10⁻⁶–0.025 kg kg⁻¹, plus q = 0; q_liq,q_ice ≥ 0 and q_liq + q_ice ≤ q; include condensate-free and precipitating cases |
| a = ∂zθ_li | −10⁻² to 3×10⁻² K m⁻¹; explicitly a = 0,d and just above/below d |
| b = ∂zq_tot | −10⁻⁵ to 10⁻⁵ (kg kg⁻¹) m⁻¹; zero and both signs |
| Horizontal gradients | |∇hθ| = 0 or 10⁻⁷–10⁻³ K m⁻¹; |∇hq| = 0 or 10⁻¹¹–10⁻⁶ (kg kg⁻¹) m⁻¹; parallel, antiparallel, perpendicular and arbitrary angles |
| X,Z | X = 10–250,000 m; Z = 10–1,000 m; use actual node horizontal length scale and per-cell depth, not nominal element width |
| d | Default 10⁻³ K m⁻¹; positive sensitivity values 10⁻⁴,10⁻²; zero/negative are invalid-domain configuration probes |
| c_g,c_Δx,c_Δz | Default test expectations 1/12,1,1; c_g ≥ 0; scales 0,0.5,1,2, and a sequence approaching zero |
| C_cov,L_mix | Nonnegative closure coefficient; L_mix = 0–1,000 m as a representative diagnostic range; use C_cov = 0 to isolate resolved blocks |
| ρ_turb,r_max,J | ρ_turb ∈ [−1,1], r_max default 0.5; finite J from the library's thermodynamics, not an independently prescribed constant |
| Precision | Float32 and Float64, with identical physical states rounded to each type |

Representable extremes outside these ranges are robustness tests, not evidence of atmospheric validity. Finite valid inputs exclude NaNs, infinities, zero cell thickness and nonpositive floor values.

## 5. Numerical check plan

Use public simulation construction and registered diagnostics as the acceptance route. With `import ClimaAtmos as CA`, construct configurations using `CA.AtmosConfig(Dict(...); job_id=...)` and `CA.get_simulation(config)`, or construct a manufactured setup with `CA.AtmosModel(grid; setup=..., params=..., ...)` and `CA.AtmosSimulation(model; ...)`. Use a small `CA.BoxGrid`, then a single-column and a spherical grid. Supply compatible temperature, density, moisture and condensate initial conditions through the setup interface, using thermodynamics to realize θ_li. Parameter overrides belong in a TOML passed through `"toml" => [path]`.

At the initialized state, read Y = sim.integrator.u, p = sim.integrator.p and t = sim.integrator.t. Obtain registered variables through `CA.Diagnostics.get_diagnostic_variable(name)` and invoke their `.compute(Y,p,t)` for `env_temperature_variance`, `env_q_tot_variance`, `env_q_tot_temperature_covariance`, and `env_q_tot_temperature_correlation`. Materialize returned lazy fields before another cache evaluation. These are the public diagnostic observations; the scalar helpers and `set_covariance_cache!` are internal functions, not substitutes for an API-level acceptance test. A supplemental internal test may call them to localize a failure. Follow the [simulation API](https://github.com/CliMA/ClimaAtmos.jl/blob/ba2bf1225fc48058d709853c065b1db21a412204/docs/src/api.md) and [diagnostics registration interface](https://github.com/CliMA/ClimaAtmos.jl/blob/ba2bf1225fc48058d709853c065b1db21a412204/src/diagnostics/Diagnostics.jl).

Use locally affine fields on interior elements and interior vertical levels; do not impose an affine jump across a periodic seam. Add smooth periodic fields for full-domain/filter tests. Form independent reference invariants with the same metric, restriction and DSS weights, evaluated in high precision, and apply local r only after filtering. An analytic gradient is an appropriate oracle for affine interior tests but not an exact oracle for underresolved or spatially filtered fields. Include nodes adjoining elements, where DSS is active. Varying profiles must test spatially varying r as well as constant-r cancellation.

For an output F with characteristic uncancelled block scale S of the same units, require |F−F_ref| ≤ τ S, where τ = 10⁻¹⁰ for Float64 and 2×10⁻⁵ for Float32 in affine interior tests. For cancellation in G, use S_G = C + 2|rB| + r²A; do not demand relative accuracy against a zero residual. For PSD use C² ≤ V_T V_q + τ max(V_T V_q,C²) and separately nonnegative diagonals; for bounds use the corresponding variance scale. Require exact zero for deliberately zero algebraic inputs when thermodynamic/gradient reconstruction is exact, otherwise the same scaled absolute tolerance. For non-affine geometry tests compare to a discrete oracle at τ = 10⁻⁸ / 10⁻⁴, and separately demonstrate convergence toward the continuous reference over three mesh refinements. Do not loosen tolerances to hide a finite partition discrepancy.

| Invariant | API-level configurations, sampling, and assertions |
|---|---|
| I1 | Set C_cov = 0, c_Δz = 0 and c_Δx = 1; choose q means large enough to avoid the moisture cap. Sample a above d, both b signs, zero gradients, q = αθ+β, nearly cancelling gradients, and independent vectors at angles 0, π/2, π and seeded random angles (at least 100 valid profiles/type). Check V_q = hfG and V_T = 0. Check diagnosed sampled cross covariance = 0; a nonzero fallback correlation is allowed when V_T = 0. Rotate horizontal fields by 90° on a square grid and reverse axes. Repeat below the floor against the regularized, not physical-surface, formula. |
| I2 | With c_Δz = 0 and nonzero horizontal q gradient but b = 0, infer f from V_q/(hC). Choose x = 0,10⁻³,0.1,1,10,100 and a = −d,0,d(1±10⁻³),2d,10d. Compare f to the formula and check the excursion inequalities; with nonzero b and selected cross-gradients also infer the floored ratio through the predicted G. Supplemental direct helper calls check the ratio bound/type stability without inversion ambiguity. |
| I3 | Exercise both `tq_correlation_model` values, ρ_turb = −1,−0.5,0,0.5,1, and horizontal/vertical scales independently on/off. Test zero turbulence and nonzero turbulence. Construct uncapped, just-capped, strongly capped and q=0 cases; compare the full transformed-and-capped reference, then PSD and correlation bounds. Tune amplitudes across V_T V_q = eps(FT)² and V_q = eps(FT)²; expect fallback rather than correlation invariance there. Registered diagnostics alone cannot reveal an invalid raw covariance hidden by the clamp: compare to the independently assembled pre-clamp reference, and optionally inspect `p.precomputed.ᶜT′q′` in a supplemental internal test. |
| I4 | Test all four exact (c_Δx=0/nonzero, c_Δz=0/nonzero) combinations, flat and steep θ surfaces, and c_Δx,c_Δz = 0.5,1,2. With C_cov=0, ∇hq=0 and both scales=1, verify V_q=c_g b²Z² below the moisture cap. For c_Δx=0 verify the full uniform-box result; with both scales=0 verify zero variance and fallback. Check transformed V_T, not θ variance, against J²v a². |
| I5 | Hold A>0,b≠0,c_Δz=1,C_cov=0 and the same state fixed. Sweep c_Δx=1,10⁻¹,…,10⁻⁶,0 with inactive moisture cap. Compare the limiting q variance and diagnosed covariance to the value at exactly zero. The expected jump is (1−f)vb² in q variance and (1−f)vab before the J transform in cross covariance. |

For all changed geometric and partition branches use `microphysics_model="1M"`, `cloud_model="quadrature"`, `use_sgs_quadrature=true` as the main path. Repeat representative cancellation, cap, and PSD cases for `0M`, `2M`, and `2MP3` where supported by their own compatible setups; the changed closure has no species-specific formula. Cover condensate-free, liquid-only, ice-only, mixed-phase and precipitation-containing states because θ reconstruction and J read different moisture quantities. Repeat with `turbconv=nothing`, `edonly_edmfx`, and `prognostic_edmfx`: the latter two use the cached EDMF mixing length, and prognostic EDMFX still uses grid-mean gradients, not environmental gradients.

Exercise `sgs_distribution="gaussian"`, `"lognormal"`, and `"mean"` at a representative state; these should not change this covariance assembly for identical thermodynamic fields and zero turbulent coefficient. With turbulent feedback enabled, compare each run to its own realized mixing length, not blindly to another run. Also cover the `uses_covariances` alternatives: quadrature present; quadrature disabled but 1M/2M requiring covariances; quadrature disabled with `QuadratureCloud` or `MLCloud`; and a 0M/GridScaleCloud configuration with no consumer. For the last case assert successful no-op, not existence of unallocated covariance diagnostics. MLCloud coverage needs its network artifact. Use a single column with c_Δx=0 and separately attempt nonzero c_Δx only if its horizontal-scale operation is supported; an unsupported combination is not an isentropic-formula failure. Add a flat sphere and terrain-following grid to check horizontal metrics versus vertical projections.

Probe d=0 and d<0 through configuration construction separately. Record whether rejected; do not count subsequent NaNs as a violation under the valid d>0 domain. No claim is made that the existing code rejects them. A one-step `CA.solve_atmos!(sim)` smoke run can check finite diagnostics after cache refresh; algebraic acceptance should use initialization so dynamical evolution cannot mask errors.

## 6. Do the PR's own tests check these invariants

The PR adds two helper testsets to `test/parameterized_tendencies/microphysics/cloud_fraction.jl` and one default-value assertion to `test/parameter_tests.jl`. It adds no simulation/cache integration test. Existing cloud-fraction inversion and cover tests in that file do not directly test the changed covariance assembly. The following assessment concerns assertions attributable to this PR, not an assertion that no other repository test could ever exercise the code. See the [merged helper tests](https://github.com/CliMA/ClimaAtmos.jl/blob/ba2bf1225fc48058d709853c065b1db21a412204/test/parameterized_tendencies/microphysics/cloud_fraction.jl).

| Invariant | Coverage | Exact evidence |
|---|---|---|
| I1 | No | No added call to `isentropic_gradient_invariant`, no cancellation/rotation test, and no check that the assembled horizontal θ variance/cross term vanishes. |
| I2 | Partial | `isentropic_slope_ratio: plain ratio with a stratification floor` tests stable, at-floor, neutral, overturned and near-zero a, the sign and bound for one negative b, and Float32/64 type preservation. `isentropic_cell_cap: along-surface segment bounded by the cell` tests flat f=1, x=1 giving f≈1/2, a steep formula/asymptote, monotonicity with stratification, identical floored results, 0<f≤1 at five A values, and types. It does not check the length/share identity through the cache, positive b, or parameter rejection. |
| I3 | No | No added PSD, transformed covariance, cap-rescaling, correlation-mode or fallback test. |
| I4 | Partial | The cap helper's f=1 and steep-limit assertions check ingredients only. No test checks f's placement in vertical q variance/cross covariance, the switch combinations, or the total-q partition identity. |
| I5 | No | No test varies c_Δx toward zero. |

`Geometric SGS variance parameters` verifies the new d default is FT(1e−3) in both precisions. It does not establish positive-value validation or a new-floor override test. The reproducibility counter bump is not a numerical invariant assertion.

## 7. Expected verdict

- **I1: expected pass for the declared discrete model; integration unverified.** The scalar polynomial has the correct minus cross term and zero floor, and the horizontal θ/cross increments are removed. Near cancellation requires absolute, scale-aware tolerance. General spatially varying-r physical isentrope accuracy and terrain/phase consistency remain uncertain; the filtered discrete formula is not an exact pointwise derivative in all such cases.
- **I2: expected pass in the positive-floor, nonoverflowing domain.** The formulas establish the bounds directly. The stated constructor validation is missing; rejection of invalid d is not expected from `get_sgs_quadrature`.
- **I3: expected pass under valid coefficients and thermodynamics.** Positive semidefinite blocks, scalar Jacobian transformation and cap rescaling support this conclusion; diagnosed clamping further bounds the sampled correlation. Full integration has not been run. An out-of-range prescribed correlation is outside the invariant's domain, especially because the fallback/constant branches return it verbatim.
- **I4: expected pass for the executable f-partition convention on flat grids; general-geometry verification remains open.** The advertised absence of double counting is justified by the stated special partition identity, not by a universal proof for arbitrary unequal scale factors or arbitrary fields.
- **I5: expected fail if continuous amplitude semantics are intended.** For every nonzero c_Δx, f is computed from X,Z,A,D without that coefficient; at exactly zero f is replaced by 1. Hence the vertical reduction survives as the horizontal variance tends to zero. This is a finite jump, not roundoff. Whether it is a defect or an intentional regime switch requires a modelling decision.

The asymmetric vertical partition, fictitious stable surface below the floor, and omission of moisture sensitivity in a one-scalar θ-to-T transform are modelling conventions. The PR alone cannot establish that they faithfully represent all moist parcel thermodynamics.

## 8. Environment

Use ClimaAtmos at merge `ba2bf1225fc48058d709853c065b1db21a412204`, not simply the tag `v0.42.14`: the merged `Project.toml` still reports that version. The base commit is needed only for optional off-branch regression comparisons.

Use Julia compatible with the project's `julia = "1.10"` constraint and instantiate its project/test dependencies. Relevant compatibility entries include ClimaCore `1.0.3`, ClimaParams `1.3`, Thermodynamics `1.2`, CloudMicrophysics `0.43, 0.44`, ClimaComms `0.6.12`, ClimaDiagnostics `0.3.11`, and ClimaTimeSteppers `0.10.4, 1`; these are compatibility constraints, not a frozen resolved manifest. Record Julia and the resolved package manifest. The source of these constraints is the [merged Project.toml](https://github.com/CliMA/ClimaAtmos.jl/blob/ba2bf1225fc48058d709853c065b1db21a412204/Project.toml).

A CPU, single process, small grid and both floating types suffice for the primary checks. GPU/MPI runs are optional portability checks, not prerequisites for the mathematical specification. Standard Julia Test and high-precision reference arithmetic are sufficient for scalar oracles; no external experimental dataset is required. Full simulation initialization may retrieve standard package artifacts, and MLCloud needs its configured network artifact. Actual runtime and test success are unmeasured in this review.

## 9. Formalization feasibility

I1's chain rule can be stated in Lean for differentiable real fields under a≠0; its constant-gradient/affine-cancellation and orthogonal-invariance lemmas reduce to inner-product algebra. The floored version is a piecewise real function, and its nonnegativity follows from squared norms or from the Gram inequality B²≤AC. Lifting that inequality through F requires proving positivity and common weights of the actual restriction and DSS operators; it cannot be assumed for an arbitrary linear filter.

I2 is directly provable over reals with d,X,Z>0 by splitting max(a,d), proving nonnegative x² and clearing positive denominators. The share identity and steep-limit theorem are elementary algebra and limits. I3 follows from positive-semidefinite matrix closure under sums/congruences and f²≤f; the cap/fallback branches require explicit domain and threshold splits. I4 follows by substitution. I5 admits a short formal counterexample with A>0,b≠0,v>0 and c_Δx→0⁺.

These proofs would establish properties of a transcribed real-valued closure, not automatically verify the Julia implementation. Obstacles to a whole-library proof include floating-point cancellation/overflow, `ifelse` evaluation, mutable scratch aliases and broadcasts, discrete metric/projection identities, thermodynamic Jacobians across phases, and nonlinear cloud-fraction/mixing-length feedback. Proving the empirical closure physically correct is outside the scope of a proof assistant; its modelling assumptions must remain explicit axioms or hypotheses.
