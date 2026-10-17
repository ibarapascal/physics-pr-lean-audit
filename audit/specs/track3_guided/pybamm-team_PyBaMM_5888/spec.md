## 1. What changed

PR #5888, “fix: convert Ai2020 electrolyte diffusivity from cm2/s to m2/s,” changes PyBaMM from base `9d8df9d039e6b8ec78f94db6a916d7f5c003fbd7` to merge `343ece2062dbb39789dbfa977e0f2daf96e9e370`.

In `packages/pybamm/src/pybamm/input/parameters/lithium_ion/Ai2020.py`, `electrolyte_diffusivity_Ai2020(c_e, T)` acquires a factor of `1e-4`. With numerical inputs (c=c_e\,[\mathrm{mol\,m^{-3}}]), (T\,[\mathrm K]), and (q=T-229-0.005c), the old and new numerical outputs are

\[
D_{\rm old}=10^{-4.43-54/q-0.00022c},\qquad
D_{\rm new}=10^{-4}D_{\rm old}
            =10^{-8.43-54/q-0.00022c}\;\mathrm{m^2\,s^{-1}}.
\]

The old numerical value was in cm²/s but was supplied under the SI parameter key. The corrected docstring identifies electrolyte diffusivity in m²/s and cites Valøen and Reimers (2005), equation 14.

In the same file, `electrolyte_conductivity_Ai2020(c_e, T)` receives documentation corrections only: concentration in mol/m³, temperature in K, and conductivity in S/m. Its unchanged formula is

\[
\begin{aligned}
P(c,T)={}&-10.5+0.668\,10^{-3}c+0.494\,10^{-6}c^2\\
 &+(0.074-1.78\,10^{-5}c-8.86\,10^{-10}c^2)T\\
 &+(-6.96\,10^{-5}+2.8\,10^{-8}c)T^2,\\
\kappa(c,T)={}&10^{-4}cP(c,T)^2\;\mathrm{S\,m^{-1}}.
\end{aligned}
\]

`get_parameter_values()` still registers these functions as `"Electrolyte diffusivity [m2.s-1]"` and `"Electrolyte conductivity [S.m-1]"`. No equation, option, electrode phase, or solver is added. `CHANGELOG.md` records the correction. `packages/pybamm/tests/unit/test_parameters/test_parameter_sets/test_Ai2020.py` adds two reference values, replaces the function-table absolute tolerance with `rel=1e-3, abs=0`, and adds an Ai2020/Xu2019 comparison.

## 2. Invariants

Define the physically intended liquid-fit domain (\mathcal U=\{c\geq0,\ T>229+0.005c\}), restricted in practice to the operating window in section 4. Mathematical statements about the fit do not establish experimental validity throughout this unbounded domain.

**I1 — Primary: dimensional calibration of salt diffusivity.** Let (C=c/1000) be concentration in mol/L and

\[
F(C,T)=10^{-4.43-54/(T-229-5C)-0.22C}.
\]

For every ((c,T)\in\mathcal U), the public Ai2020 parameter must return (10^{-4}F(c/1000,T)) in m²/s, because (1\,\mathrm{cm^2/s}=10^{-4}\,\mathrm{m^2/s}). In particular,

\[
D(1000,298.15)=3.222722528605257\times10^{-10}\;\mathrm{m^2/s}.
\]

The displayed number is an independent evaluation of the formula, not a measured value. Equivalently, the new/base ratio is (10^{-4}) wherever both evaluations are finite and nonzero. Empirical coefficients carry the implied units needed to make the exponent dimensionless.

**I2 — Shared-correlation and model-interface consistency.** On (\mathcal U), (D_{\rm Ai2020}(c,T)=D_{\rm Xu2019}(c,T)): both parameter sets specify the same electrolyte correlation. This does not imply equality of their complete cell models. The model-facing `LithiumIonParameters.D_e` must evaluate to

\[
D_{\rm model}(c,T)=D_{\rm Ai2020}(\max(c,c_{\min}),T),
\]

where (c_{\min}=\texttt{pybamm.settings.tolerances["D_e__c_e"]}=10\) mol/m³ at this commit. For this second identity require (c\geq0) and (T>229+0.005\max(c,c_{\min})). Below the threshold, raw-function and model-facing values need not coincide. This floor is a numerical convention already present in the library, not a physical assertion about dilute electrolytes.

**I3 — Positive, finite diffusion coefficient.** For finite inputs in (\mathcal U), exact-real (D(c,T)>0). On the bounded test window it must also be finite and strictly positive in floating point. Consequently, the library's diffusive flux law (\mathbf N_d=-\mathcal B D\nabla c\) obeys (\mathbf N_d\cdot\nabla c=-\mathcal B D\lVert\nabla c\rVert^2\leq0) for nonnegative transport efficiency (\mathcal B). This concerns the diffusion component; migration and convection need not be down the concentration gradient. The numerical check targets the constitutive coefficient, not a new PDE conservation test.

**I4 — Temperature and concentration dependence of this fit.** Throughout (\mathcal U),

\[
\frac{\partial\ln D}{\partial T}=\ln(10)\frac{54}{q^2}>0,\qquad
\frac{\partial\ln D}{\partial c}=-\ln(10)\left(\frac{0.27}{q^2}+0.00022\right)<0.
\]

Thus diffusivity increases with temperature at fixed concentration and decreases with concentration at fixed temperature, on intervals contained in the domain. These are properties of this chosen correlation, not universal laws for all electrolytes. For the model-facing floor, concentration dependence is constant below 10, decreasing above 10, and not classically differentiable at 10.

**I5 — Conductivity contract preserved.** The documented public conductivity remains (\kappa=10^{-4}cP^2\) in S/m. For finite (c\geq0,T), its exact-real value is nonnegative and the raw function has (\kappa(0,T)=0). The nominal value is approximately (1.1943\) S/m at ((1000,298.15)). The model wrapper similarly replaces (c) by `max(c, settings.tolerances["kappa_e__c_e"])`, also 10 at this commit; the raw zero-concentration limit is therefore not a requirement on the wrapper. No Nernst–Einstein proportionality between this salt diffusivity and conductivity is imposed.

**Conventions and limits.** The source uses base-10 exponentiation. Xu2019 explicitly explains that the paper's “ln” in this correlation is interpreted as log10. A literal natural-exponential reading would define a different calibration; the library convention and reported diffusivity scale resolve the implementation target. Neither the function documentation nor this PR specifies clipping, exceptions, or an extrapolation rule at (q\leq0). The pole (q=0), its low-temperature side, negative concentrations, NaN, and infinity are outside the physical acceptance domain. The library's EC:DMC description also does not prove that the original PC/EC/DMC fit is quantitatively accurate for every solvent blend.

## 3. Sources of the invariants

- **I1:** Exact SI length conversion and the empirical fit in Valøen and Reimers, *Transport Properties of LiPF₆-Based Li-Ion Battery Electrolytes*, JES 152, A882–A891 (2005), [DOI](https://doi.org/10.1149/1.1872737). Equations 14–16 and Table II give the correlation coefficients; the [author-uploaded paper](https://www.researchgate.net/publication/260511055_Transport_Properties_of_LiPF_6-Based_Li-Ion_Battery_Electrolytes) discusses the measured temperature/concentration trends. Its logarithm notation needs the convention explained above.
- **I1–I2:** The registered SI keys and independent molarity-form implementation in [merged Ai2020.py](https://raw.githubusercontent.com/pybamm-team/PyBaMM/343ece2062dbb39789dbfa977e0f2daf96e9e370/packages/pybamm/src/pybamm/input/parameters/lithium_ion/Ai2020.py) and [merged Xu2019.py](https://raw.githubusercontent.com/pybamm-team/PyBaMM/343ece2062dbb39789dbfa977e0f2daf96e9e370/packages/pybamm/src/pybamm/input/parameters/lithium_ion/Xu2019.py). Ai2020 also cites Ai et al. (2020), [DOI 10.1149/2.0122001JES](https://doi.org/10.1149/2.0122001JES); that citation establishes intended provenance, not an independently verified calibration in this review.
- **I2 and I5 wrapper behavior:** [lithium_ion_parameters.py](https://raw.githubusercontent.com/pybamm-team/PyBaMM/343ece2062dbb39789dbfa977e0f2daf96e9e370/packages/pybamm/src/pybamm/parameters/lithium_ion_parameters.py), methods `D_e` and `kappa_e`, and [settings.py](https://raw.githubusercontent.com/pybamm-team/PyBaMM/343ece2062dbb39789dbfa977e0f2daf96e9e370/packages/pybamm/src/pybamm/settings.py).
- **I3:** Positivity of a real exponential and the diffusive constitutive law used by `Full.get_coupled_variables` in [full_diffusion.py](https://raw.githubusercontent.com/pybamm-team/PyBaMM/343ece2062dbb39789dbfa977e0f2daf96e9e370/packages/pybamm/src/pybamm/models/submodels/electrolyte_diffusion/full_diffusion.py).
- **I4:** Direct differentiation of the specified fit; empirical monotonicity is not assumed outside this model.
- **I5:** The unchanged polynomial and SI contract, checked against the base file as well as the merge file. Nonnegativity follows algebraically from its square. The source cites Torchio et al. (2016), *LIONSIMBA*, JES 163, A1192, but no new conductivity calibration is introduced by this PR.

## 4. Operating domain

Only concentration and absolute temperature enter the changed constitutive functions. Use (500\leq c\leq2000\) mol/m³ (0.5–2 mol/L) and (273.15\leq T\leq333.15\) K (0–60 °C) as a representative battery-model checking window. This is a chosen practical window, not a claim that documentation guarantees experimental accuracy over the entire rectangle. The original paper discusses diffusion measurements from −10 to +60 °C; the precise validated concentration envelope is not established here.

The Ai2020 initial concentration is 1000 mol/m³, and reference, ambient, and initial temperatures are 298.15 K. In the proposed window (q\geq34.15\) K, safely away from the pole, and the fit spans approximately (3.54\times10^{-11}\) to (8.49\times10^{-10}\) m²/s. Additional (c=0\) and concentrations near 10 are mathematical/interface boundary checks, not calibrated dilute-solution predictions. Derivative probes remain within the window. No current, time, state of charge, particle stoichiometry, or geometry is an argument of these invariants; the full-cell voltage examples in the PR body are not acceptance thresholds.

## 5. Numerical check plan

Use the merged library's parameter API; do not substitute a standalone transcription for the function under test:

```python
import numpy as np
import pybamm

ai = pybamm.ParameterValues("Ai2020")
xu = pybamm.ParameterValues("Xu2019")
D_key = "Electrolyte diffusivity [m2.s-1]"
k_key = "Electrolyte conductivity [S.m-1]"

def value(params, key, c, T):
    return float(params.evaluate(
        params[key](pybamm.Scalar(c), pybamm.Scalar(T))
    ))

lp = pybamm.LithiumIonParameters()
# Model-facing calls:
# ai.evaluate(lp.D_e(pybamm.Scalar(c), pybamm.Scalar(T)))
# ai.evaluate(lp.kappa_e(pybamm.Scalar(c), pybamm.Scalar(T)))
```

Sample the full Cartesian product of concentrations `[500, 750, 1000, 1250, 1500, 1750, 2000]` and temperatures `[273.15, 278.15, 298.15, 318.15, 333.15]`, plus 200 uniformly sampled pairs with `np.random.default_rng(5888)`. Cartesian coverage is important: the PR's three pairs change both variables together.

- **I1:** Compare public outputs with an independently organized oracle using (C=c/1000), (T_g=229+5C), and `np.exp(np.log(10) * (-4.43 - 54/(T-T_g) - .22*C)) * 1e-4`. Require `rtol=1e-12, atol=0, equal_nan=False` and explicit finiteness. Check the precise nominal value at the same tolerance. If the base is installed separately, compare new/base values to (10^{-4}) with relative tolerance (10^{-12}); do not mix imports from two checkouts.
- **I2:** Compare Ai2020 and Xu2019 over every sample at `rtol=1e-12, atol=0, equal_nan=False`. Compare the model-facing `D_e` call to the raw Ai2020 function with the concentration floor applied. Cover below/equal/above threshold using (c=0,5,9.999,10,10.001,20) at every listed temperature, reading the setting before running and requiring its commit-default value 10. Include all normal samples, where the floor is inactive.
- **I3:** Require `np.isfinite(D)` and `D > 0` at every normal and boundary sample for both raw and model-facing evaluations, without a tolerance allowing zero or negatives. This checks the coefficient condition underlying down-gradient diffusion; it does not claim to test total electrolyte flux or the discretization.
- **I4:** Require strictly increasing values along every temperature row and strictly decreasing values along every raw concentration column. At the random interior samples with (c\in[501,1999]), (T\in[273.16,333.14]), differentiate `np.log(value(...))` by central differences with (h_c=0.1\) mol/m³ and (h_T=0.001\) K. Compare to the analytic logarithmic derivatives in section 2 at `rtol=1e-6, atol=0`. For the wrapper require constancy below the floor (`rtol=1e-12, atol=0`), continuity across it within the formula-based comparison, and strict decrease above it. Do not request a derivative at the kink.
- **I5:** Evaluate public conductivity at all normal and boundary samples. Compare to the independently transcribed polynomial in section 1 at `rtol=1e-12, atol=0, equal_nan=False`; require finite nonnegative results and exactly zero for the raw (c=0) case. Check the rounded nominal value with `rel=1e-3, abs=0`. Verify `kappa_e` against the same formula with its own concentration floor, covering below/equal/above 10. This also ensures that the diffusivity conversion was not accidentally applied to conductivity.

The changed functions have no internal branches, phase selectors, hysteresis, or option values. Scalar calls cover their full semantic case distinction; additionally evaluate paired `pybamm.Vector` inputs and each scalar/vector combination on the grid and compare each element to its scalar evaluation (`rtol=1e-12, atol=0`) to check field-wise use. The unchanged wrapper's concentration floor is the only relevant branch added to the coverage plan. There is no need to sweep unrelated battery-model options or run a DFN discharge to isolate this unit correction.

## 6. Do the PR's own tests check these invariants

The test names below are in `packages/pybamm/tests/unit/test_parameters/test_parameter_sets/test_Ai2020.py`.

| Invariant | Coverage | Actual assertions and limitations |
|---|---|---|
| I1: SI calibration | **Partial** | `TestAi2020.test_functions` checks (D(1000,298.15)\approx3.2227\times10^{-10}) with `rel=1e-3, abs=0`; this directly rejects the factor-of-10000 bug. `test_electrolyte_diffusivity` adds cross-set agreement at three points, but no independent formula sweep. |
| I2: correlation and interface consistency | **Partial** | `TestAi2020.test_electrolyte_diffusivity` compares Ai2020 to Xu2019 at `(1000,298.15)`, `(500,318.15)`, `(1500,278.15)`, using `np.testing.assert_allclose(..., rtol=1e-12)`. It does not check `D_e`, the floor, or vector inputs. |
| I3: finite positive diffusion | **Partial** | The nominal positive reference indirectly checks positivity and finiteness at one point. Neither test explicitly asserts them over a domain. The cross-comparison alone can accept equal negative values or matching NaNs. |
| I4: monotonicity and derivatives | **No** | There are no fixed-temperature/fixed-concentration sweeps or derivative assertions. The three cross-set pairs do not isolate either dependence. |
| I5: conductivity contract | **Partial** | `TestAi2020.test_functions` checks (\kappa(1000,298.15)\approx1.1943) with `rel=1e-3, abs=0`; no domain sweep, zero-concentration limit, or wrapper-floor check. |

`np.testing.assert_allclose` defaults to `atol=0`, so the new comparison really is a relative check, unlike `np.allclose`'s defaults. Its default `equal_nan=True` explains the finiteness limitation above. See the [NumPy API documentation](https://numpy.org/doc/stable/reference/generated/numpy.testing.assert_allclose.html).

The stricter table tolerance also improves regression checks of unchanged solid diffusivities and cracking rates. The PR body describes earlier testing and reported integration-test counts, but those reports are not additional invariant assertions or independently rerun evidence.

## 7. Expected verdict

- **I1: Expected pass.** The exact missing cm²-to-m² conversion is now present. The base implementation fails the calibrated nominal value by a factor of 10000.
- **I2: Expected pass.** Substituting (C=c/1000) makes the two fit expressions algebraically equal; the inspected model wrapper retains the concentration floor.
- **I3: Expected pass on the stated window.** A positive prefactor multiplies a positive exponential, with exponents comfortably representable there. Positivity alone would also have passed before the fix.
- **I4: Expected pass.** Multiplication by a positive constant preserves the existing derivative signs and logarithmic sensitivities. This invariant alone cannot detect the original unit error.
- **I5: Expected pass.** The computational body is unchanged and is a nonnegative concentration times a square.

These are source-based expectations, not results of executing the merged PyBaMM API. Experimental accuracy for a particular EC:DMC blend and behavior near the singularity remain **uncertain/outside scope**. No exact final voltage or concentration range from the PR narrative is promoted to an invariant.

## 8. Environment

Use the complete PyBaMM monorepo at merge `343ece2062dbb39789dbfa977e0f2daf96e9e370`, including packaged parameter CSV files and parameter-set entry points. Install `packages/pybamm` from that checkout in an isolated environment. The [package metadata at the merge](https://raw.githubusercontent.com/pybamm-team/PyBaMM/343ece2062dbb39789dbfa977e0f2daf96e9e370/packages/pybamm/pyproject.toml) requires Python `>=3.10,<3.16`, NumPy `>=2.0.0`, SciPy `>=1.11.4,!=1.18.0`, and `pybammsolvers>=0.11.0,<0.12.0`, plus its declared runtime dependencies. Use a supported 64-bit Python such as 3.12 and float64 inputs. Record resolved dependency versions and concentration-floor settings.

The constitutive checks need no time integration, GPU, or optional plotting package. Pytest is useful for assertions; running the repository test file under its own configuration additionally needs its declared pytest plugins, including xdist, mock, and benchmark. Run that test file with `-n 0` for a small local check. An optional base comparison needs a second isolated install at the base commit. This review downloaded and inspected source and evaluated the reference arithmetic, but did not install or execute PyBaMM.

## 9. Formalization feasibility

All five invariants can be stated over transcribed real-valued functions in Lean. Represent decimal coefficients as exact rationals and (10^x) as (\exp(x\log10)).

- **I1:** Prove the unit factor from the square of (1/100), and use exponential identities for the equivalent exponent −8.43. Quantities or explicit unit tags are needed to distinguish physical units from bare real numbers.
- **I2:** Rational algebra proves the molarity substitution; case splits on (c\leq c_{\min}) prove the wrapper identity.
- **I3:** Exponential positivity and positivity of (10^{-4}) are straightforward; the flux consequence follows from a nonnegative squared norm.
- **I4:** Chain/quotient rules with (q>0) prove the logarithmic derivatives and monotonicity. The floor requires separate constant and differentiable regions and a continuity proof at the threshold.
- **I5:** Polynomial identities and square nonnegativity prove the contract and dilute limit; the wrapper again needs a maximum-function case split.

The main obstacles are outside the elementary real mathematics: linking the transcription to Python/PyBaMM expression evaluation and parameter registration, formalizing IEEE-754 rounding and exponentiation error to justify numerical tolerances, and proving scientific validity of empirical coefficients. A proof assistant can verify consistency with the selected correlation and units; it cannot establish solvent-specific experimental accuracy from this code alone. The domain restriction is essential because proof-assistant division at zero may have totalized semantics unlike Python's evaluation.
