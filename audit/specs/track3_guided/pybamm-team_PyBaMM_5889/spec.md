## 1. What changed

PR #5889 changes one constitutive lookup in `packages/pybamm/src/pybamm/input/parameters/sodium_ion/Chayambuka2022.py`: `NVPF_diffusivity_Chayambuka2022(sto, T)`. Write the dimensionless stoichiometry as θ, positive-electrode maximum solid concentration as C, initial electrolyte concentration as E, and the piecewise-linear interpolation of `data/D_p.csv` as L. The change is

\[
D_{p,\mathrm{base}}(\theta,T)=L(\theta E)
\quad\longrightarrow\quad
D_{p,\mathrm{merge}}(\theta,T)=L(\theta C).
\]

The defaults are C = 15,320 mol m⁻³ and E = 1,000 mol m⁻³. The CSV explicitly labels its abscissa as positive particle concentration, not stoichiometry or electrolyte concentration. Its 34 rows span 131.578947–15,197.36842 mol m⁻³. Neither the table, the interpolation algorithm, nor the unused temperature argument changes. The analogous negative-electrode function already uses its own maximum solid concentration.

In `packages/pybamm/tests/unit/test_parameters/test_parameter_sets/test_Chayambuka2022.py`, `test_functions` switches all eight constitutive-function checks to `rel=1e-3, abs=0`, updates the positive diffusivity expectation, corrects the electrolyte diffusivity expectation from 2.5061e-10 to 2.5061e-11 m² s⁻¹, and increases precision of several other reference values. No electrolyte constitutive formula changes. The new `test_positive_particle_diffusivity` checks θ = 0.3 and 0.9 with `rel=1e-4, abs=0`.

`CHANGELOG.md` records the fix. `docs/source/examples/notebooks/models/sodium-ion.ipynb` has refreshed outputs, including its voltage–capacity plot; its simulation code is unchanged. It runs `pybamm.sodium_ion.BasicDFN()` at C-rates 1/12, 5/12, 10/12, and 1. The PR body's reported 1C/2C capacities are simulation observations, not invariants or asserted tests.

## 2. Invariants

Let the committed table be ordered pairs (xᵢ, dᵢ), with strictly increasing xᵢ and positive dᵢ. For xᵢ ≤ c ≤ xᵢ₊₁, define

\[
L(c)=(1-\lambda)d_i+\lambda d_{i+1},\qquad
\lambda=\frac{c-x_i}{x_{i+1}-x_i}.
\]

**I1 — Correct solid-concentration coordinate (primary).** For C > 0, 0 ≤ θ ≤ 1, and c = θC in the tabulated interval,

\[
D_p(\theta,T;C,E)=L(c),\qquad D_p(x_i/C,T;C,E)=d_i.
\]

This is an identity of the parameterization, not a claim that the experimental table is exact physics. It includes consistency with the model-facing concentration argument. `ParticleLithiumIonParameters.D` in `packages/pybamm/src/pybamm/parameters/lithium_ion_parameters.py` divides c by C before dispatching to this function. It also clips stoichiometry: with τ = `pybamm.settings.tolerances["D__c_s"]` (default 10⁻¹⁰), the actual model-facing identity is

\[
D_{\mathrm{model}}(c,T)=L\!\left(C\,\operatorname{clip}(c/C,\tau,1-\tau)\right).
\]

Thus cancellation to L(c) applies only when clipping is inactive. The direct parameter callable does not clip.

**I2 — Independence from unrelated electrolyte initialization (primary).** At fixed θ, C, T and fixed solid-diffusivity table,

\[
D_p(\theta,T;C,E_1)=D_p(\theta,T;C,E_2)
\quad\text{for all }E_1,E_2>0.
\]

Equivalently, its partial derivative with respect to E is zero. As a complementary coordinate-consistency test, if c is fixed and θⱼ = c/Cⱼ is admissible, then Dₚ(θ₁,T;C₁,E₁) = Dₚ(θ₂,T;C₂,E₂) = L(c). Changing C here is a mathematical reparameterization with the table frozen, not a claim that different physical materials share a diffusivity. This invariance concerns the constitutive function at a fixed state; complete cell trajectories may depend strongly on E through other equations.

**I3 — Positive, continuous, non-overshooting interpolation.** On every closed tabulated interval,

\[
0<\min(d_i,d_{i+1})\le D_p(\theta,T)\le\max(d_i,d_{i+1}),
\qquad D_p\in C^0.
\]

Over the full table, 1.95×10⁻¹⁷ ≤ Dₚ ≤ 2.51×10⁻¹⁵ m² s⁻¹. Consequently the model's Fickian flux N = −Dₚ∇c obeys N·∇c = −Dₚ|∇c|² ≤ 0 for any finite concentration gradient, evaluated within this domain. This is local down-gradient diffusion; it is not a separate proof of full-cell thermodynamic consistency. The table is nonmonotonic, so neither global monotonicity nor differentiability at its knots is required.

**I4 — Temperature-independent closure.** For finite positive T₁,T₂ and a fixed admissible θ,C,

\[
D_p(\theta,T_1)=D_p(\theta,T_2).
\]

This is the implemented modeling convention: T is accepted but unused. It is not a physical assertion that NVPF diffusion is temperature-independent. The supplied parameterization gives no solid-diffusion activation energy; inventing an Arrhenius factor would change the model. Quantitative physical interpretation is restricted to the reference-temperature parameterization, 298.15 K.

The exact linear interpolation and endpoint extrapolation are library conventions, not uniquely prescribed by diffusion physics. `pybamm.Interpolant` defaults to linear interpolation and `extrapolate=True`; the changed function inherits both. Outside the table, its direct numerical value extends the first or last segment. I1–I3's table-based physical acceptance domain does not extend to arbitrary extrapolated concentrations. No endpoint clipping, zero diffusivity at saturation, or globally positive extrapolation should be presumed.

## 3. Sources of the invariants

- **I1 and I2:** the concentration definition θ = c/C and internal agreement between the [merged parameter functions](https://github.com/pybamm-team/PyBaMM/blob/eda03ed6d1003b186cc2d70ad8b3240a7b6350f1/packages/pybamm/src/pybamm/input/parameters/sodium_ion/Chayambuka2022.py), the [CSV axis and values](https://github.com/pybamm-team/PyBaMM/blob/eda03ed6d1003b186cc2d70ad8b3240a7b6350f1/packages/pybamm/src/pybamm/input/parameters/sodium_ion/data/D_p.csv), and the [model parameter adapter](https://github.com/pybamm-team/PyBaMM/blob/eda03ed6d1003b186cc2d70ad8b3240a7b6350f1/packages/pybamm/src/pybamm/parameters/lithium_ion_parameters.py). Dimensional analysis alone cannot detect the original bug: C and E have the same units but represent different physical quantities.
- **I3:** convexity of linear interpolation and positive committed ordinates; the down-gradient consequence follows from Fick's law. The [sodium-ion BasicDFN source](https://github.com/pybamm-team/PyBaMM/blob/eda03ed6d1003b186cc2d70ad8b3240a7b6350f1/packages/pybamm/src/pybamm/models/full_battery_models/sodium_ion/basic_dfn.py) explicitly uses `-param.p.prim.D(c_s_p, T) * pybamm.grad(c_s_p)` and its negative divergence. The [Interpolant implementation](https://github.com/pybamm-team/PyBaMM/blob/eda03ed6d1003b186cc2d70ad8b3240a7b6350f1/packages/pybamm/src/pybamm/expression_tree/interpolant.py) establishes the linear/extrapolation convention.
- **I4 and empirical scope:** the parameter source and example cite K. Chayambuka, G. Mulder, D. L. Danilov, and P. H. L. Notten, *Physics-based modeling of sodium-ion batteries part II. Model and validation*, Electrochimica Acta 404 (2022), 139764, [DOI](https://doi.org/10.1016/j.electacta.2021.139764). They identify the [COMSOL isothermal sodium-ion example](https://www.comsol.com/model/1d-isothermal-sodium-ion-battery-117341) as the source of the specific numerical parameter values and interpolant data. The COMSOL page confirms the isothermal application and sodium intercalation model context. The paper's full text was not accessible in this review; no independent extraction or experimental validation of its diffusivity curve is claimed. The committed CSV is the numerical oracle for this specification.

No external standard is needed to establish the coordinate error or interpolation identities.

## 4. Operating domain

| Variable | Relevant range and interpretation |
|---|---|
| θ, positive particle stoichiometry | Physical normalization 0–1; empirically tabulated range at default C is approximately 0.0085887041–0.9919953277. Initial θ = 3320/15320 ≈ 0.21671. |
| c = θC, solid sodium concentration | Table: 131.578947–15,197.36842 mol m⁻³. Nominal physical range 0–15,320 mol m⁻³; its short endpoint portions require extrapolation. |
| C, maximum positive concentration | Fixed material default 15,320 mol m⁻³. Use 0.8C, C, and 1.2C only as controlled parameter-binding tests, restricting comparisons to admissible c and θ. |
| E, initial electrolyte concentration | Default 1,000 mol m⁻³. Use 500, 1,000, 1,500 mol m⁻³ to challenge independence. These are diagnostic perturbations, not a validated electrolyte operating envelope. |
| T | Reference and default initial/ambient temperature 298.15 K. Use 273.15, 298.15, and 323.15 K to test the unused argument; this does not validate the table across that temperature range. |
| ∇c, only for the flux consequence | Any finite gradient, units mol m⁻⁴, with local c inside the table. No material-specific maximum gradient is established here. |

Current, voltage, particle radius, and time do not enter the changed constitutive function. For contextual cell simulations, the supplied example uses 1/12–1C discharge, nominal capacity 3 mAh, and voltage limits 2.0–4.2 V. The PR additionally discusses 2C. None of those observations establishes an invariant capacity or a universal relation between diffusivity and terminal capacity.

## 5. Numerical check plan

Use the pinned merged package and its committed CSV, parsed independently with `numpy.loadtxt(..., delimiter=",", skiprows=1)`. Construct an independent reference by selecting the bracketing rows and evaluating the displayed two-point formula; do not derive the oracle by evaluating another `pybamm.Interpolant`. Extend the endpoint segments explicitly for the separate extrapolation checks.

The public parameter-evaluation route is:

```python
import numpy as np
import pybamm

C_key = "Maximum concentration in positive electrode [mol.m-3]"
E_key = "Initial concentration in electrolyte [mol.m-3]"
D_key = "Positive particle diffusivity [m2.s-1]"

def evaluate(theta, T=298.15, C=15320.0, E=1000.0):
    p = pybamm.ParameterValues("Chayambuka2022")
    p.update({C_key: C, E_key: E})
    symbol = p[D_key](pybamm.Scalar(float(theta)), pybamm.Scalar(T))
    return p.evaluate(symbol)
```

Create fresh parameter objects after each override to avoid reusing processed expressions with stale parameter bindings. Unless specified otherwise, compare D using `rtol=1e-10, atol=0` and require finite outputs. This tolerance accommodates double-precision coordinate/interpolation roundoff and is much smaller than the original error; an absolute tolerance such as 10⁻⁴ m² s⁻¹ would render these checks meaningless.

- **I1:** at default C, evaluate every one of the 34 knots θᵢ = xᵢ/C, and fractions 0.01, 0.25, 0.5, 0.75, and 0.99 of each of the 33 intervals. Also sample 256 concentrations uniformly over the table with `numpy.random.default_rng(5889)`, and θ = 0.3, 0.5, 0.9. Independently computed values at those last three points are 5.6673280002e-17, 8.3317028568e-17, and 3.5626719998e-17 m² s⁻¹. Repeat the sample as a `pybamm.Vector(theta_array)` passed to `p[D_key]`, evaluated through `p.evaluate`, and compare flattened results with scalar evaluations. This covers scalar and field-like vector usage, not new physical phases.
- **I1, model wiring and boundaries:** instantiate `m = pybamm.sodium_ion.BasicDFN()` and evaluate `p.evaluate(m.param.p.prim.D(pybamm.Scalar(c), pybamm.Scalar(298.15)))`. For each interior concentration above, compare with L(c). Separately use c/C = 0, τ/2, τ, 2τ, 1−2τ, 1−τ, 1−τ/2, and 1 and compare with L(C clip(c/C,τ,1−τ)). Compare the direct callable at θ = 0, x₀/(2C), (x_last/C+1)/2, and 1 against endpoint linear extensions. These are convention checks, not additional empirically validated points. At θ = 0 and 1 the direct expected values are approximately 2.7385714278e-15 and 1.7387466637e-17 m² s⁻¹. Do not expect the clipped adapter and unclipped direct callable to be identical at the endpoints.
- **I2:** repeat the entire default-C sampling with E = 500, 1000, 1500 and compare at fixed θ. Then use C = 12,256, 15,320, 18,384 and the same E values with θ = c/C at all knots and interval samples in the common admissible concentration range [x₀, 12,256(1−τ)]. Every result must equal the independently computed L(c). This catches both retention of the electrolyte scale and hard-coding of the default C.
- **I3:** on the knot/interval/dense samples, require D > 0 and each result to lie between its bracketing ordinates, allowing a rounding margin of 10⁻¹² times their maximum. For each interior knot, evaluate c = xᵢ ± δᵢ with δᵢ = 10⁻⁶ min(xᵢ−xᵢ₋₁, xᵢ₊₁−xᵢ). Check agreement with the appropriate one-sided linear formula and the bound |D(c)−dᵢ| ≤ |slope|δᵢ + 10⁻¹²dᵢ. For the flux consequence, combine the API-returned D with gradients 0 and ±10⁹ mol m⁻⁴ and require −Dg² ≤ 0; this verifies the constitutive sign consequence, not a spatial discretization theorem.
- **I4:** repeat all default-C samples at each listed temperature, requiring agreement with the 298.15 K result at `rtol=1e-12, atol=0`.

The changed function has one positive NVPF phase, one 1D linear interpolant, and no charge/discharge, hysteresis, or selectable temperature branch. All 33 segments and both endpoint-extension cases are covered above. The model adapter adds lower-clipped, interior, and upper-clipped cases, also covered. Cubic/PCHIP interpolation, multiple active phases, and alternative diffusivity options are not selected by this parameter set. The unchanged negative HC function can serve as a control but must use its own table and C = 14,540; its values must not be equated to the positive values. No full DFN solve is necessary to decide I1–I4.

## 6. Do the PR's own tests check these invariants

| Invariant | Coverage | Test and actual assertion |
|---|---|---|
| I1 | Partial | `TestChayambuka2022.test_functions` checks Dₚ(0.5, 298.15) ≈ 8.3317e-17 with `rel=1e-3, abs=0`. `TestChayambuka2022.test_positive_particle_diffusivity` checks 0.3 → 5.6673e-17 and 0.9 → 3.5627e-17 with `rel=1e-4, abs=0`. These detect the old scale but do not cover every interval, parameter overrides, vectors, or the model adapter. |
| I2 | No | Both tests use the default E and C; neither varies an unrelated electrolyte parameter nor holds c fixed while changing C. |
| I3 | Partial | The three positive reference-value assertions indirectly ensure positivity at those three points. They do not test interval bounds, knot continuity, the complete table, or flux direction. |
| I4 | No | Both tests use only 298.15 K. |

`test_functions` also checks the negative diffusivity, two OCPs, two exchange-current densities, electrolyte diffusivity, and electrolyte conductivity at one state. These stronger numerical regressions do not add tests of I2–I4. The refreshed notebook contains plotting and solver calls, not assertions of these invariants. The author's report that seven relevant tests passed is not an independently reproduced result in this review.

## 7. Expected verdict

- **I1: expected pass.** The exact parameter-name substitution restores the c/C → θ → θC coordinate round trip, subject to the adapter's documented clipping. Independent arithmetic on the CSV reproduces the three new reference values.
- **I2: expected pass.** The merged expression contains C and no E; the same-coordinate identity follows algebraically. It would fail in the base implementation at generic states.
- **I3: expected pass within the stated table domain.** All 34 ordinates are positive and interpolation is linear. This property alone would not distinguish the bug: the old function could also return positive diffusivities at the wrong concentration.
- **I4: expected pass as an implementation convention.** Neither the old nor new solid-diffusion function uses T. Physical validity away from 298.15 K remains unvalidated.

These are source-based expectations, supported by independently calculated interpolation values, not a claim that the merged PyBaMM API or its test suite was executed. The quoted full-cell capacities and the experimental fidelity of the supplied table remain unverified here.

## 8. Environment

Pin PyBaMM to merge commit `eda03ed6d1003b186cc2d70ad8b3240a7b6350f1`; the comparison base is `c235b902faea805128b87f7581f882d965538ef9`. The changelog places the fix after v26.10.0.0, so that release alone is not an adequate pin. Keep the CSV from the same commit.

The merged `packages/pybamm/pyproject.toml` specifies Python ≥3.10,<3.16, NumPy ≥2.0.0, SciPy ≥1.11.4 excluding 1.18.0, and pybammsolvers ≥0.11.0,<0.12.0, among the normal package dependencies. Use a supported interpreter such as Python 3.12 with a compatible solver wheel; install the package's declared dependency set, including transitive CasADi support needed by its imports. Add pytest to run the PR tests. Install from the `packages/pybamm` project, since the repository-root pyproject is workspace orchestration rather than an installable package. If building the workspace solver from source, its separate native build requirements apply.

All proposed invariant checks are small CPU parameter evaluations and need no GPU, external datasets, COMSOL license, or nonlinear time integration. Matplotlib and notebook/widget packages are only needed to reproduce the example plots. Record the installed Python, PyBaMM, NumPy, SciPy, CasADi, and solver versions with any future execution results.

## 9. Formalization feasibility

I1–I4 are feasible to state and prove over a transcribed real-valued function in Lean. Encode the decimal table entries as exact rationals, prove strict ordering of the abscissae and positivity of the ordinates, and define interval selection and linear interpolation explicitly. I1 reduces to c/C · C = c under C > 0 and a separate clipping case split. I2 and I4 follow from absence of E and T in the expression. I3 follows from 0 ≤ λ ≤ 1, convex-combination bounds, equality of neighboring segment values at each knot, and nonnegativity of |∇c|².

The main obstacles are connecting this real-arithmetic transcription to Python symbolic parameter substitution, SciPy's floating-point interpolation, vector shapes, and parameter-processing caches. A proof of the transcription does not certify the executing package without a refinement/error argument. The CSV's provenance and empirical validity cannot be established by such a proof. Proving full DFN conservation, numerical stability, or capacity predictions would additionally require the PDE, boundary conditions, discretization, event handling, and solver analysis; none follows merely from this one-line constitutive correction.
