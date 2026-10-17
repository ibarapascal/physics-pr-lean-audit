## 1. What changed

CoolProp PR #2851, “fix(HEOS): correct first_saturation_deriv for mixtures (#2091),” merges as `54a2e32764bfa0ae81877e835c33ccd13ba8f853`, against base `b3360b8c6e540e88f8425a47086a93554870ed56`.

In `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp`, the overload `HelmholtzEOSMixtureBackend::calc_first_saturation_deriv(Of1, Wrt1, SatL, SatV)` now distinguishes pure/pseudo-pure fluids from mixtures. Previously both used

\[
S\equiv(dT/dp)_\sigma=T\frac{v_V-v_L}{h_V-h_L},\qquad v_\alpha=1/\rho_{m,\alpha}.
\]

Here density and enthalpy are molar. This expression is retained for `is_pure_or_pseudopure`. For mixtures, define

\[
a_{p,i}=\left(\frac{\partial\ln\phi_i^L}{\partial p}\right)_{T,x}-\left(\frac{\partial\ln\phi_i^V}{\partial p}\right)_{T,y},\quad
a_{T,i}=\left(\frac{\partial\ln\phi_i^L}{\partial T}\right)_{p,x}-\left(\frac{\partial\ln\phi_i^V}{\partial T}\right)_{p,y}.
\]

The new expression is

\[
S=-\frac{A_p}{A_T},\qquad A_k=\sum_i w_i a_{k,i},\qquad
w_i=\begin{cases}y_i&Q=0\text{ (bubble)},\\x_i&Q=1\text{ (dew)}.\end{cases}
\]

The partials use `MixtureDerivatives::dln_fugacity_coefficient_dp__constT_n` and `dln_fugacity_coefficient_dT__constp_n`, with `XN_DEPENDENT`; the last mole fraction is constrained by normalization. Despite the accumulator names `dQ_dPsat` and `dQ_dTsat`, these sums should not be interpreted as general quality derivatives inside the two-phase region. Every other mixture quality throws `ValueError`.

The two-parameter overload retains its saturation-state availability check and delegates to the four-parameter overload. The dispatch following calculation of \(S\) remains: return \(S\) for `(iT,iP)`, \(1/S\) for `(iP,iT)`, apply the homogeneous-property chain rule for other outputs with respect to `iT` or `iP`, and throw for other independent parameters.

`src/Tests/CoolProp-Tests.cpp` adds two tests, detailed in section 6. No second-derivative formula was changed, although `calc_second_saturation_deriv` calls the changed function. Its differentiation of the Clapeyron expression remains a relevant downstream limitation. [Merged implementation](https://github.com/CoolProp/CoolProp/blob/54a2e32764bfa0ae81877e835c33ccd13ba8f853/src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp#L3646).

## 2. Invariants

**I1 — Primary: the slope is tangent to the appropriate fixed-composition equilibrium boundary.** Let \(z\) be the overall mole-fraction vector. On the bubble curve \(x=z\); on the dew curve \(y=z\). At equilibrium,

\[
x_i\phi_i^L(T,p,x)p=y_i\phi_i^V(T,p,y)p,\quad \sum_i x_i=\sum_i y_i=1.
\]

For a smooth local boundary \(p=p_q(T;z)\), the returned slope must satisfy

\[
S_q=\left(\frac{dp_q}{dT}\right)^{-1},\qquad A_T S_q+A_p=0.
\]

The phase whose composition is not fixed must re-equilibrate along this boundary. Domain: positive temperature and pressure, positive component fractions, stable liquid–vapor coexistence, differentiable EOS, nonsingular local equilibrium problem, and \(A_T\ne0\). Reciprocal comparisons additionally require finite nonzero \(dp_q/dT\). Critical coalescence, spinodals, and turning points requiring a different curve parameter are excluded from regular finite-slope assertions. No universal positive-slope assertion is imposed on all mixture branches.

**I2 — Primary: pure-fluid thermodynamic consistency and pseudo-pure compatibility.** For a true pure fluid away from criticality,

\[
\frac{dp_\sigma}{dT}=\frac{s_V-s_L}{v_V-v_L}=\frac{h_V-h_L}{T(v_V-v_L)},\qquad S=\frac{T\Delta v}{\Delta h}.
\]

This is independent of which saturation endpoint is selected. Domain: a single-component equilibrium liquid–vapor curve with nonzero latent heat and volume jump. For a pseudo-pure fluid with distinct fitted bubble/dew curves, equality to the derivative of each fitted curve is **not** implied by pure-fluid thermodynamics. The compatibility invariant here is preservation of the base implementation's Clapeyron-style result using the supplied saturation states, separately at each endpoint. The documentation does not resolve an exact thermodynamic interpretation of that approximation for temperature-gliding pseudo-pure fluids.

**I3 — endpoint property chain rule and coordinate consistency.** For a supported homogeneous property \(F(T,p,z)\) of the endpoint phase,

\[
D_pF=F_p+F_TS,\qquad D_TF=F_T+F_p/S,\qquad D_TF=D_pF/S.
\]

Also \(D_pT\,D_Tp=1\), \(D_TT=1\), and \(D_pp=1\). Domain: I1's regular endpoints or I2's true pure-fluid curve, with nonzero finite \(S\) and differentiable \(F\). At a bubble endpoint this describes the liquid at fixed \(z\); at a dew endpoint it describes the vapor at fixed \(z\). A property of the *other*, composition-changing phase would require additional composition derivatives and is not covered by this formula. For pseudo-pure fluids the algebraic chain rule with the retained \(S\) is the compatibility target, not necessarily a derivative of an independently fitted boundary.

**I4 — endpoint-only mixture API contract.** For a successfully established mixture state with \(0<Q<1\), `first_saturation_deriv` must reject the request; for valid regular \(Q=0,1\) it must accept supported output/independent-variable combinations. Unsupported `Wrt1` must throw. This is a library contract, not a statement that all interior-quality derivatives are mathematically undefined: a specified fixed-\(Q\), fixed-\(z\) path can have a derivative, but this API does not implement that path. No new rejection requirement is imposed on pure-fluid interior states by this PR.

**I5 — component relabeling invariance.** For any permutation \(\pi\) of species and their fractions, \(S_q(T;\pi z,\pi\text{species})=S_q(T;z,\text{species})\), with correspondingly permuted phase compositions. Domain: the same physical EOS and interaction parameters, within I1's regular domain. Choosing the last fraction as dependent must not select a physically special component.

**I6 — downstream derivative consistency (secondary).** Wherever the public second saturation derivative is supported on a smooth endpoint curve,

\[
\texttt{second\_saturation\_deriv}(F,p,p)=\frac{d}{dp}\big[\texttt{first\_saturation\_deriv}(F,p)\big]_{q,z}.
\]

Domain: twice-differentiable regular true-pure or mixture boundaries, away from singularities. Pseudo-pure fitted-curve exactness is excluded. This is a consistency requirement on an existing caller, not a claim that the PR explicitly promises to repair mixture second derivatives.

## 3. Sources of the invariants

I1 follows from equality of component chemical potentials, \(f_i=x_i\phi_i p\), and Gibbs–Duhem. At the bubble point, differentiate \(\ln x_i+\ln\phi_i^L-\ln y_i-\ln\phi_i^V=0\) with \(dx=0\), multiply by \(y_i\), and sum. Normalization eliminates \(\sum_i dy_i\); Gibbs–Duhem eliminates the composition contribution \(\sum_i y_i d_y\ln\phi_i^V\). The remaining equation is \(A_TdT+A_pdp=0\). At the dew point the same argument uses \(dy=0\) and weights \(x_i\). Thus the opposite-phase weights follow from equilibrium thermodynamics, rather than just from the patch.

The code attributes the expressions to equations 3.96–3.97 of Georg Johannes Gernert, *A New Helmholtz Energy Model for Humid Gases and CCS Mixtures*, Ruhr-Universität Bochum, 2013. The [repository bibliography](https://github.com/CoolProp/CoolProp/blob/54a2e32764bfa0ae81877e835c33ccd13ba8f853/CoolPropBibTeXLibrary.bib#L870) confirms that citation; the thesis equations were not independently inspected here. The derivation above is the independent basis of I1.

I2 follows from \(d\mu=-s\,dT+v\,dp\), equality of pure-phase chemical potentials, and \(\Delta h=T\Delta s\). I3 and I6 follow from the chain rule and the definition of a derivative. The [public API header](https://github.com/CoolProp/CoolProp/blob/54a2e32764bfa0ae81877e835c33ccd13ba8f853/include/AbstractState.h#L1315) describes the saturation chain rule and cites Thorade and Saadat, “Partial derivatives of thermodynamic state properties for dynamic simulation,” *Environmental Earth Sciences* 70, 3497–3503 (2013). I4 is supported by the [endpoint-only documentation](https://github.com/CoolProp/CoolProp/blob/54a2e32764bfa0ae81877e835c33ccd13ba8f853/Web/coolprop/HighLevelAPI.rst#L325) and the new explicit exception. I5 expresses invariance to species indexing. No external standard mandates a particular numerical tolerance or pseudo-pure interpretation here.

## 4. Operating domain

Use SI units: \(T\) in K, \(p\) in Pa, \(\rho_m\) in mol/m³, \(h,u\) in J/mol, \(s\) in J/(mol K), \(S\) in K/Pa. \(Q\) is molar vapor quality; \(z,x,y\) are dimensionless normalized mole fractions. Pressure and phase compositions are equilibrium outputs when sampling `(Q,T)`, not independently selectable inputs.

A practical test domain is:

| Configuration | Temperature sampling | Composition / quality |
|---|---|---|
| HEOS Methane&Ethane | 140–220 K in 10 K increments, explicitly 180 K | Methane fraction 0.1, 0.3, 0.5, 0.7, 0.9; both endpoints |
| HEOS Methane&Ethane&Propane | 160–220 K in 10 K increments | (0.7,0.2,0.1), (0.3,0.4,0.3); both endpoints |
| HEOS Water | 300, 400, 500, 550 K | Both endpoints |
| HEOS n-Propane | 230, 270, 300, 330 K | Both endpoints |
| HEOS R410A as one pseudo-pure fluid | 250, 280, 310 K | Both endpoints; compatibility only |

These are proposed ordinary cryogenic/refrigeration/steam ranges, not verified universal flash-success ranges. Target positive saturation pressures roughly \(10^4\)–\(10^7\) Pa and remain within each EOS's validity bounds. Record unavailable or failed flashes rather than counting them as passing checks. Exclude coalescing phases and unresolved turning points from regular error statistics, while reporting the excluded configurations. Require converged stencils and a resolvable phase-density gap. For rejection tests use \(Q=0.1,0.5,0.9\) at the known regression configuration, 180 K and (0.7,0.3). Exact-zero fractions and pure-component limits represented as mixtures are outside the main matrix because absent-component handling introduces separate degeneracies.

## 5. Numerical check plan

Use the merged Python public API, `import CoolProp as CP` and `from CoolProp import AbstractState`. Create `A=AbstractState('HEOS', fluid_string)`, set mixture mole fractions, and call `A.update(CP.QT_INPUTS,q,T)`. Use separate state instances for finite differences; restore the center before each analytic measurement. Do not force a homogeneous phase on the saturation-flash state.

For every numerical comparison use \(|a-b|\le r\max(|a|,|b|)+\epsilon_{abs}\). For property derivatives choose \(\epsilon_{abs}=10^{-8}F_*/X_*\), with \(F_*\) the maximum of a characteristic SI scale (1 in the output's units) and the magnitudes of center/stencil outputs, and \(X_*=T\) or \(p\). For slope use \(\epsilon_{abs}=10^{-8}T/p\). Convergence of the numerical reference must be demonstrated before interpreting a disagreement.

- **I1:** For every mixture, composition, temperature, and both qualities in section 4, compare `A.first_saturation_deriv(CP.iT,CP.iP)` with \(2\delta T/[p(T+\delta T,q,z)-p(T-\delta T,q,z)]\). Use \(\delta T=0.05,0.025,0.0125\) K; require successive reference slopes to agree within \(5\times10^{-4}\) relative, then compare with \(r=10^{-3}\). Reproduce the PR's exact 0.05 K, 0.5% test separately. Also use `PQ_INPUTS,p±δp,q` at fixed composition, with \(\delta p/p=10^{-4},5\times10^{-5},2.5\times10^{-5}\), to check \([T(p+\delta p)-T(p-\delta p)]/(2\delta p)\). This checks both flash parameterizations and avoids treating reciprocal agreement as independent physical evidence.
- **I1, independent weighted residual:** Obtain `mole_fractions_liquid()`, `mole_fractions_vapor()`, and `saturated_liquid_keyed_output(CP.iDmolar)` / vapor equivalent. Construct separate homogeneous liquid/vapor states with those compositions; use `specify_phase(CP.iphase_liquid/gas)` and `PT_INPUTS` near the center. Verify center densities match the equilibrium branches. Finite-difference `log(B.fugacity_coefficient(i))` at fixed phase composition in \(T\) and \(p\) to estimate every \(a_{k,i}\). Check \(|A_TS+A_p|/(|A_TS|+|A_p|)<2\times10^{-3}\) with nonvanishing denominator and step refinement. This is entirely public API access; a root switch invalidates that stencil. Also verify normalized compositions to \(10^{-10}\) and \(\max_i|\ln(f_i^L/f_i^V)|<10^{-6}\) as equilibrium preconditions.
- **I2:** At every pure-fluid endpoint, read saturated molar densities and enthalpies with the saturated keyed-output methods and compare \(S\) with \(T\Delta v/\Delta h\) at \(r=10^{-10}\); independently use the flash finite differences above at \(r=10^{-4}\). Compare bubble/dew slopes at the same pure-fluid temperature at \(r=10^{-8}\). For R410A, repeat the expression check and compare base versus merge using isolated builds at \(r=10^{-10}\); do not demand coincident bubble/dew curves or exact finite-difference agreement.
- **I3:** At each regular endpoint test outputs `iDmolar`, `iHmolar`, `iSmolar`, `iUmolar` and mass counterparts, with both `iT` and `iP`. Compare `first_saturation_deriv` against finite differences of `keyed_output(F)` along fixed \(q,z\), using the preceding stencils, \(r=10^{-3}\) for mixtures and \(10^{-4}\) for true pure fluids. Test the algebraic chain rule using `first_partial_deriv(F,iT,iP)` and `first_partial_deriv(F,iP,iT)` at \(r=10^{-10}\). Check reciprocal and identity relations to \(10^{-10}\) absolute. This covers both direct return branches and both general-property branches.
- **I4:** At each specified interior quality call all four forms `(iT,iP)`, `(iP,iT)`, `(iHmolar,iT)`, `(iHmolar,iP)` and require exceptions. At endpoints require successful finite results for the regular cases. Use `(iT,iDmolar)` at endpoints to exercise the unsupported-independent-variable exception. Check exception type/message where exposed, without requiring an exact formatted floating-point string. Never classify a failed state update as successful derivative rejection.
- **I5:** Reverse each binary order and test all six ternary permutations, permuting fractions identically. Compare saturation pressure, slopes, and reindexed phase fractions at \(r=10^{-7}\), absolute fraction tolerance \(10^{-9}\). Any differing interaction-parameter selection must be diagnosed separately from summation roundoff.
- **I6:** On the same true-pure and mixture endpoints, compare `second_saturation_deriv(F,iP,iP)` for `F=iT,iDmolar,iHmolar,iSmolar,iUmolar` against central pressure differences of `first_saturation_deriv(F,iP)` at fixed \(q,z\). Use the pressure steps above, require reference convergence within 0.1%, and compare at \(r=5\times10^{-3}\), \(\epsilon_{abs}=10^{-8}F_*/p^2\). Report unsupported calls separately from wrong finite results. Test pure-fluid cases as controls.

All exceptions, rejected stencils, and singular cases must be reported with the configuration. Tolerance relaxation without a convergence study is not evidence of correctness. The internal explicit-SatL/SatV overload cannot be selected directly through `AbstractState`; public-API tests cover the delegating overload and reachable callers, not arbitrary internal state combinations.

## 6. Do the PR's own tests check these invariants

The two new tests and relevant existing tests were inspected in the [merged test source](https://github.com/CoolProp/CoolProp/blob/54a2e32764bfa0ae81877e835c33ccd13ba8f853/src/Tests/CoolProp-Tests.cpp#L4815).

| Invariant | Coverage | Actual assertions |
|---|---|---|
| I1 | **Partial** | New `first_saturation_deriv mixture: dT/dP via Gernert vs finite-difference along bubble curve` (`[first_saturation_deriv][2091]`) uses Methane&Ethane, (0.7,0.3), 180 K, Q=0. It checks positive slope and relative discrepancy below `5e-3` against the reciprocal central pressure difference with 0.05 K steps. No dew, other composition, ternary, pressure stencil, or independent fugacity residual. |
| I2 | **Partial** | Existing `Check the first partial derivatives` (`[first_saturation_partial_deriv]`) compares ten derivatives against finite differences for n-Propane, Q=1, 300 K at `1e-4` relative tolerance. It includes T/P and P/T, but no explicit latent-heat identity, liquid endpoint, or pseudo-pure compatibility check. No new pure-fluid test is added. |
| I3 | **Partial** | That existing test checks Dmolar, Hmolar, Smolar, Umolar with respect to T and P at the same pure-vapor state, using `1e-5` K or `1e-2` Pa perturbations. No mixture property derivatives, mass-output cases, or explicit reciprocal/identity assertions. |
| I4 | **Partial** | New `first_saturation_deriv mixture: throws at intermediate Q` (`[first_saturation_deriv][2091]`) checks only that `(iT,iP)` throws at Q=0.5, 180 K, (0.7,0.3). It does not inspect the exception reason or test other rejected qualities/dispatch branches. |
| I5 | **No** | Neither added test nor the cited saturation regression tests permutes species. |
| I6 | **Partial** | Existing `Check the second saturation derivatives` (`[second_saturation_partial_deriv]`) checks T, Dmolar, Hmolar, Smolar, Umolar differentiated twice with respect to pressure against differences of first saturation derivatives, for n-Propane, Q=1, 300 K, `1e-2` Pa steps and `1e-4` relative tolerance. No mixture case. |

The PR body reports successful regression runs and red/green testing; these claims were not reproduced for this specification. Its cited `[derivatives]` tests concern ordinary partial derivatives and do not establish the new bubble/dew tangent formula.

## 7. Expected verdict

- **I1: expected pass in the regular domain, not numerically verified here.** The opposite-phase weighting follows from Gibbs–Duhem and the fixed endpoint composition. Both quality branches implement the correct weighting. Conditioning near critical/turning points remains uncertain; there is no new denominator guard.
- **I2: expected pass for true pure fluids and for the stated pseudo-pure compatibility invariant.** The expression is unchanged. Exact pseudo-pure fitted-curve tangency is uncertain and is not claimed.
- **I3: expected pass for regular endpoint bulk properties.** The endpoint phase has fixed overall composition, so there is no missing composition term for these properties. This expectation does not extend to properties of the incipient opposite phase along the same boundary.
- **I4: expected pass.** The mixture branch explicitly rejects every quality other than exact 0 or 1 before dispatch. Invalid/uninitialized states are not evidence about valid endpoint behavior.
- **I5: expected pass, uncertain until permutation checks.** The formula is a symmetric component sum; helper derivatives and interaction-parameter lookup still need empirical verification.
- **I6: expected pass for true pure fluids; expected generic failure or unsupported behavior for mixtures.** The second-derivative implementation uses the corrected first slope but still constructs its derivatives from \(T\Delta v/\Delta h\). That is generally not the derivative of the new fugacity-based slope, including the varying incipient-phase composition. Accidental agreement at a special state would not resolve the inconsistency. This conclusion is a source-based expectation, not an observed test failure.

## 8. Environment

Build CoolProp at merge commit `54a2e32764bfa0ae81877e835c33ccd13ba8f853`, including its pinned submodules/dependency configuration, fluid data, mixture interaction data, and generated headers. An arbitrary released wheel is insufficient evidence that this merged code is running. Record `CoolProp.CoolProp.get_global_param_string('gitrevision')`, version, compiler, architecture, and configuration flags.

For the proposed public checks use the Python wrapper built from that checkout, Python 3 and its build dependencies (including Cython as required by the checkout); NumPy is convenient for sampling but not necessary. No REFPROP installation or proprietary data is needed. For the existing C++ tests use CMake, a supported C++ compiler, and the checkout's Catch2 dependencies, configuring `COOLPROP_CATCH_MODULE=ON` and running `CatchTestRunner` with the tags listed above. Use a separate base-commit build for pseudo-pure compatibility. The review retrieved source and documentation at pinned commits but did not build or execute CoolProp.

## 9. Formalization feasibility

The algebraic part is feasible in Lean over real-valued transcribed functions. With finite component index sets, normalized positive compositions, differentiable chemical potentials/fugacity coefficients, Gibbs–Duhem, equilibrium equalities, a differentiable local boundary, and nonzero denominators as hypotheses, one can prove the opposite-phase weighted tangent relation, the Clapeyron identity, the endpoint chain rule, reciprocal identities, and permutation invariance. The quality dispatch and delegation equivalence are straightforward finite case splits. I6 could be stated precisely and used to expose the difference between differentiating the new slope and the retained Clapeyron-based expression.

The main obstacles are proving those hypotheses for the actual Helmholtz EOS, existence and local uniqueness of the selected flash root, phase stability, derivative correctness of the mixture helpers, and handling critical degeneracies. Proving the real-number formula would not establish the C++ implementation's floating-point error bounds, cache/state validity, root selection, exceptions, or convergence. Pseudo-pure compatibility is a program-equivalence property; it cannot supply a missing physical modelling convention or prove exact tangency of independently fitted saturation curves.
