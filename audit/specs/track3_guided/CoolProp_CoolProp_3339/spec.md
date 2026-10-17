## 1. What changed

CoolProp PR #3339 changes `TransportRoutines::viscosity_rhosr` in `src/Backends/Helmholtz/TransportRoutines.cpp`. Base: `afcc79091533dfbd258cd59c065c0765747fd230`; merge: `5418cec58c52183e3a010b2894e08b818ed4cc81`.

For molar density \(\rho\), temperature \(T\), dimensionless residual Helmholtz energy \(\alpha^r\), and inverse reduced temperature \(\tau\), define

\[
s^r=R(\tau\alpha^r_\tau-\alpha^r),\quad
x=\rho s^r/q_c,\quad q_c=\mathtt{rhosr\_critical},
\]
\[
L(x)=\sum_{j=0}^3c_{L,j}x^j,\quad V(x)=\sum_{j=0}^3c_{V,j}x^j,
\quad w(x,a)=\frac1{1+\exp[-100(x-a)]}.
\]

The old implementation used \(a=2\); the merged implementation uses \(a=\mathtt{x\_crossover}\). The sharpness remains 100, now a named local constant. The unchanged surrounding formula is

\[
g=wL+(1-w)V,\quad E=\exp(g),\qquad
\eta=\eta_0(T)\{1+C[E-1]\}.
\]

Here \(\eta_0\) is the kinetic-theory dilute viscosity in Pa s, computed by `viscosity_dilute_kinetic_theory` in the same file. It uses positive Lennard–Jones parameters and the Neufeld collision-integral approximation. The logistic interpolation acts on **log reference reduced viscosity**, before the fluid-specific affine transformation by \(C\); it is not an arithmetic blend of viscosities.

`parse_rhosr_viscosity` in `src/Backends/Helmholtz/Fluids/FluidLibrary.h` already loaded the crossover field. `parse_viscosity` selects the first entry when viscosity data are an array. `HelmholtzEOSMixtureBackend::calc_viscosity_contributions` dispatches to this correlation and places its entire result in the contribution named `critical`; that name does not mean this formula is solely a critical enhancement. These routines are context, not changed code.

All eight relevant `dev/fluids/<name>.json` files use integer `x_crossover: 2`: R1233zd(E), R1234yf, R1234ze(E), R124, R152A, R22, R245fa and R32. The PR changes no coefficients or fluid files. It adds one test in `src/Tests/CoolProp-Tests.cpp`.

## 2. Invariants

**I1 — Primary: the configured crossover controls the prescribed smooth log blend.** For any valid homogeneous EOS state and finite model parameters, with \(q_c\ne0\), four coefficients in each branch, and representable intermediates, the returned viscosity must obey the equations in section 1 with the supplied \(a\). In real arithmetic,

\[
w(a,a)=\tfrac12,\quad w(a+d,a)+w(a-d,a)=1,\quad
\partial_xw=100w(1-w)>0,\quad \partial_aw=-100w(1-w).
\]

In particular, at a fixed state with \(a=x\),

\[
\eta=\eta_0[1+C(\exp((L+V)/2)-1)],
\quad
\partial_a\eta=-100\eta_0CEw(1-w)(L-V).
\]

Thus increasing the crossover reduces viscosity when \(C>0\) and \(L>V\), increases it when \(L<V\), and has no effect if \(L=V\) or \(C=0\). A universal claim that moving the crossover must change every result would be incorrect. Strict inequalities for the weight become non-strict after floating-point saturation. Smoothness in physical state variables is asserted only where the EOS is smooth, not across a two-phase flash boundary.

**I2 — Primary: backward compatibility at the shipped crossover.** For each of the eight unchanged fluid models, at every state where both revisions return a finite viscosity,

\[
\eta_{\rm merge}(T,\rho;a=2)=\eta_{\rm base}(T,\rho),
\qquad \eta_{\rm clone,a=2}=\eta_{\rm shipped}.
\]

This is a numerical compatibility invariant, not an experimental accuracy claim. Bitwise identity is expected when compiler, dependencies, floating-point mode, data and state are identical: the JSON integer 2 converts exactly to the double formerly supplied by the literal. Cross-build comparisons should instead allow rounding tolerance.

**I3 — Bounded blending and conditional physical positivity.** For finite real branch values,

\[
\min(L,V)\le g\le\max(L,V),\quad
\min(e^L,e^V)\le E\le\max(e^L,e^V).
\]

For \(C\ge0\) and \(\eta_0>0\), this also bounds \(\eta\) between the two branch endpoint values \(\eta_0[1+C(e^L-1)]\) and \(\eta_0[1+C(e^V-1)]\). Positive viscosity is required for passive Newtonian-fluid dissipation. Algebraically it is guaranteed for \(0\le C\le1\), and for \(C>1\) exactly when \(E>1-1/C\). R245fa has \(C=1.21398\), so arbitrary coefficient sets are not unconditionally positive. With the shipped coefficients and \(x\ge0\), both polynomials are nonnegative, giving \(E\ge1\) and \(\eta\ge\eta_0>0\) for all eight fluids. The polynomial derivatives have positive minima on this domain. No unconditional finiteness assertion is made for arbitrarily large \(x\): cubic exponentiation can overflow.

**I4 — Dilute-limit consistency, including its convention-dependent qualification.** Along a fixed positive temperature, as \(\rho\to0^+\), a regular residual EOS gives \(s^r\to0\) and \(x\to0\). Hence the exact implemented limit is

\[
\lim_{\rho\to0^+}\frac{\eta}{\eta_0}
=1+C\left[\exp\left(w(0,a)c_{L,0}+[1-w(0,a)]c_{V,0}\right)-1\right].
\]

For shipped coefficients, \(c_{V,0}=0\), \(c_{L,0}=0.6100913843\). With \(a=2\), the difference from unity is of order \(10^{-87}\), below double precision. For a custom \(a=0\), the limit is instead \(1+C[\exp(0.6100913843/2)-1]\). Therefore exact recovery of \(\eta_0\) for every user-supplied crossover is **not** an invariant of this model. The inspected code supplies no admissible range of custom crossover values or normalization adjustment. Whether such edits are physically calibrated models remains open; the parameter-wiring fix alone does not settle that modelling choice.

## 3. Sources of the invariants

- I1 and I2 follow from the library's data contract and internal consistency: [merged transport implementation](https://raw.githubusercontent.com/CoolProp/CoolProp/5418cec58c52183e3a010b2894e08b818ed4cc81/src/Backends/Helmholtz/TransportRoutines.cpp), [base implementation](https://raw.githubusercontent.com/CoolProp/CoolProp/afcc79091533dfbd258cd59c065c0765747fd230/src/Backends/Helmholtz/TransportRoutines.cpp), and [fluid parser](https://raw.githubusercontent.com/CoolProp/CoolProp/5418cec58c52183e3a010b2894e08b818ed4cc81/src/Backends/Helmholtz/Fluids/FluidLibrary.h). Logistic symmetry, derivatives and convexity are elementary mathematical consequences, not additional fitted assumptions.
- The active fluid blocks cite Bell and Laesecke, *Viscosity of refrigerants and other working fluids from residual entropy scaling*, 16th International Refrigeration and Air Conditioning Conference at Purdue (2016), key `Bell-PURDUE-2016-ETA` in the [repository bibliography](https://raw.githubusercontent.com/CoolProp/CoolProp/5418cec58c52183e3a010b2894e08b818ed4cc81/CoolPropBibTeXLibrary.bib). This citation identifies the model; the full paper was not used to assert independent parameter-validity bounds.
- I3's blend bounds follow from convexity and monotonicity of the exponential. Physical positivity follows from nonnegative viscous dissipation, \(2\eta D':D'\ge0\), for arbitrary deviatoric strain rate \(D'\). The stronger coefficient-specific result follows from the actual shipped polynomials.
- I4 follows from \(s^r=-(\partial a^r/\partial T)_\rho\), \(a^r=RT\alpha^r\), the regular ideal-gas limit of the residual EOS, and substitution into the implemented blend. Residual entropy is signed; replacing it by its magnitude without also changing the signed normalizer is wrong. No thermodynamic first principle fixes the empirical crossover at 2 or the sharpness at 100.

## 4. Operating domain

Target the pure-fluid HEOS backend. Mixtures, alternative ECS models in later array entries, other transport correlations, and bulk two-phase viscosity mixing rules are outside this change. Both saturated endpoints are meaningful checks; intermediate quality is excluded.

The shipped EOS limits below are thermodynamic envelopes, **not validated viscosity ranges**. All pressures are absolute.

| Fluid | Triple point / K | EOS maximum T / K | EOS maximum p / MPa |
|---|---:|---:|---:|
| R1233zd(E) | 165.75 | 450 | 100 |
| R1234yf | 121.60 | 410 | 100 |
| R1234ze(E) | 168.62 | 420 | 15 |
| R124 | 75.00 | 470 | 40 |
| R152A | 154.56 | 500 | 58 |
| R22 | 115.73 | 550 | 60 |
| R245fa | 171.05 | 440 | 200 |
| R32 | 136.34 | 435 | 70 |

These values were read from `EOS[0]` in the eight merge-commit JSON files. A practical refrigeration/working-fluid check envelope is 250 K to \(\min(450\ {
m K},T_{\max})\), with positive pressure up to \(\min(20\ {
m MPa},p_{\max})\). Sample density through the EOS, approximately dilute gas to \(3\rho_c\) (critical molar densities are about 3704–8150 mol/m³), rejecting solid, two-phase-interior, unstable, or out-of-envelope states. Density is not independently bounded by a universal validated maximum. Extend toward each triple point only as an explicitly labelled boundary test, with melting information where available.

Relevant model inputs: \(a=2\) in production, \(C=0.6682\)–1.21398, \(q_c=-107473.04\) to \(-48503.79\) J/(m³ K). All fluids use \(c_L=(0.6100913843,0.4508958312,-0.017063,0.000564)\) and \(c_V=(0,1.2,-0.308598,0.035317)\). Thus \(x,a,L,V,C\) are dimensionless. Kinetic inputs require \(T>0\), molar mass \(M>0\), collision diameter \(\sigma>0\), and \(\epsilon/k>0\); these are fixed by the selected fluid and parser defaults, not extra user operating controls. Custom crossovers below are numerical probes, not experimentally validated parameter ranges.

## 5. Numerical check plan

Use the merged Python wrapper's public APIs: `CoolProp.CoolProp.get_fluid_param_string(name, "JSON")`, `add_fluids_as_JSON("HEOS", json_text)`, `AbstractState("HEOS", name)`, `update(DmolarT_INPUTS, rho, T)`, `viscosity()`, `gas_constant()`, `tau()`, `alphar()`, and `dalphar_dTau()`. Compute signed \(x\) from these public EOS outputs and the JSON normalizer. Also cross-check selected states with `PropsSI("V", "T", T, "Dmolar", rho, name)`.

For arrays, edit the first viscosity element; for R1233zd(E), edit the object. Register uniquely named clones with unique CAS identifiers, empty aliases and stripped identifying metadata, following the PR test. Use a fresh process per fluid to contain registrations. Do not submit malformed JSON, short coefficient arrays or a zero normalizer: this PR adds no validation for them.

For each of the eight fluids, build a state grid with 12 temperatures across the practical envelope and 24 logarithmic pressures from 1 kPa to the pressure ceiling using `PT_INPUTS`; record the resulting density and phase. Supplement with `QT_INPUTS` at \(Q=0,1\) and \(T/T_c=0.70,0.85,0.95\) where allowed, and supercritical states at \(T/T_c=1.01,1.05\) where below \(T_{\max}\), scanning \(\rho/\rho_c=0.1,0.5,1,1.5,2,3\). Filter by the stated physical envelope and report rejected points and achieved phase coverage. Do not silently count rejected states as passes. Include the PR state, R1233zd(E), 300 K and 10000 mol/m³.

1. **I1:** Register a \(C=0\) clone to obtain \(\eta_0(T)\) through the same public API; this is a control for testing the changed blend, not an independent validation of kinetic theory. For each fixed state, clone with \(a=x+d\), \(d\in\{-0.2,-0.05,-0.01,0,0.01,0.05,0.2\}\), plus fixed \(a=2,12\). Evaluate the expected formula using an independently written, overflow-safe logistic and direct polynomial sums. Require relative viscosity error at most \(5\times10^{-12}\). The sweep covers vapor-dominated, liquid-dominated, midpoint and both transition shoulders independently of whether a physical density scan encounters the transition. At \(a=x\), explicitly test the geometric-reference midpoint. Recover \(g=\log[1+(\eta/\eta_0-1)/C]\) and test the complement identity \(g(x+d)+g(x-d)=L+V\), with absolute error \(10^{-10}\max(1,|L|,|V|)\). Here the arguments denote crossover values at a fixed state. Test the analytic \(a\)-derivative using centered differences with \(h=10^{-5}\) and \(5\times10^{-6}\); allow \(10^{-5}\) relative error away from cancellation/saturation. At degeneracy, use an absolute tolerance scaled by \(\eta_0\max(1,CE|L-V|)\). Include \(a=x\pm8\) as floating-point tail probes: under ordinary nontrapping IEEE arithmetic, overflow of the inner exponential should yield the correct zero weight, not a nonfinite viscosity.
2. **I2:** Run the same state list in isolated base and merge builds and compare all eight shipped fluids. Require identical bits for matched reproducible builds, otherwise relative error \(10^{-13}\). Compare each merge shipped fluid to its unchanged \(a=2\) JSON clone at \(10^{-14}\), matching the PR's assertion. Record build differences before interpreting a bitwise failure.
3. **I3:** On the state and crossover sweeps, require finite positive \(\eta\), check both transformed branch bounds with slack \(5\times10^{-12}\) times the largest absolute endpoint or returned viscosity, and check \(\eta\ge\eta_0\) wherever \(x\ge0\). Report negative-\(x\) states separately using the general positivity condition, especially for R245fa. Add equal-polynomial clones and \(C=0\) controls to check that crossover changes become inert in these degenerate cases. No monotonicity with density or temperature is assumed.
4. **I4:** At 250 K, 300 K and \(\min(400\ {
m K},T_{\max})\), use dilute vapor states with \(\rho/\rho_c=10^{-6},10^{-8},10^{-10},10^{-12}\). For each fluid and \(a=2,12,0\), compare to the finite-density formula at \(5\times10^{-12}\) relative tolerance, then require the last two points to agree with the exact limiting expression within \(10^{-8}\) relative error. For \(a=2,12\), additionally check \(\eta/\eta_0\approx1\); for \(a=0\), explicitly check the nonunit limit rather than imposing an incorrect normalization.

These are proposed checks; no compiled-library numerical run was performed for this specification.

## 6. Do the PR's own tests check these invariants

The sole added case is `rhosr-CS viscosity honours x_crossover from the fluid file`, tagged `[viscosity],[transport],[rhosr]`, in `src/Tests/CoolProp-Tests.cpp`.

| Invariant | Coverage | Actual assertions and gaps |
|---|---|---|
| I1 | Partial | Requires the source model to be `rhosr-CS` with crossover 2; registers crossovers 2 and 12 without throwing; at 300 K, 10000 mol/m³ requires a relative change greater than 0.01. Does not assert the formula, midpoint, sign, sharpness, symmetry, derivative, actual x, or other fluids. |
| I2 | Partial | Requires the crossover-2 clone to agree with the shipped fluid within \(10^{-14}\) at one state. Does not execute the base revision or cover all fluids/states. |
| I3 | Partial | Requires `ValidNumber` for the three viscosities at that state. Does not explicitly assert positivity or either blend bound. |
| I4 | No | No dilute state or dilute-limit assertion in the added test. |

The pre-existing `R1233zd(E) has a viscosity model (#3330)` checks saturated-liquid viscosity at 273.15 K is finite and between \(10^{-4}\) and \(10^{-3}\) Pa s, checks the bibliography key, and pins four viscosities to \(10^{-6}\) relative tolerance: dilute gas (250 K, \(10^{-10}\) mol/m³), vapor (350 K, 100 mol/m³), and two compressed-liquid states. It provides limited existing regression/positivity/dilute evidence for I2–I4, but does not test configurable crossover semantics. The PR body's reported suite success is not an independently reproduced result here.

## 7. Expected verdict

- **I1: expected pass** for valid finite-range inputs. The merged expression uses the parsed field with the correct sign and unchanged sharpness. Tail behavior with floating-point exceptions enabled is uncertain; the implementation does not use an overflow-safe logistic.
- **I2: expected pass**, including bitwise preservation under matching builds, because the shipped value is exactly 2. Cross-build last-bit behavior is uncertain until run.
- **I3: expected pass** on the stated shipped-coefficient, nonnegative-x domain. Arbitrary custom coefficients, extreme extrapolation and \(C>1\) without the stated positivity condition have no guarantee.
- **I4: expected pass** for the qualified limiting formula and effectively unit shipped limit. A stronger demand for exact unit normalization for every configurable crossover would fail, notably at \(a=0\); that is a pre-existing modelling limitation exposed by making the parameter operative.

## 8. Environment

Build CoolProp at merge `5418cec58c52183e3a010b2894e08b818ed4cc81` with its matching submodules and generated fluid data. Use a C++17 compiler (the checked CMake file sets this standard), CMake, and the repository's Python-wrapper build dependencies, including Python and Cython as required by that build. The Python check harness needs only standard-library JSON/math utilities beyond the compiled wrapper; a root finder is optional because crossover sweeps avoid solving for x. No REFPROP installation or external fluid database is needed.

Build the base revision separately for I2 and verify each wrapper's revision before importing it in its own process. Use ordinary double-precision, non-fast-math, nontrapping arithmetic for the main comparison. Build the repository test executable with `COOLPROP_CATCH_MODULE=ON` to run the added `[rhosr]` case and existing `[3330]` case. Capture compiler, dependency versions and floating-point settings. Obtaining source files and reviewing formulas was sufficient for this specification; the builds and numerical checks remain to be run.

## 9. Formalization feasibility

I1, I3 and I4 can be stated and proved in Lean over a real-valued transcription of the polynomial/logistic/exponential function. Logistic range, midpoint, complement symmetry and derivative are elementary; monotonicity of exp gives the blend bounds. Positivity requires the explicit hypotheses on \(C\), \(E\) and \(\eta_0\); shipped-polynomial nonnegativity can be proved using exact decimal rationals and quadratic derivative bounds. I4 needs an assumed or separately proved regular residual-EOS limit.

I2 is immediate by substituting \(a=2\) in the real formula, but bitwise equivalence additionally needs a floating-point semantics, exact JSON conversion, identical operation ordering and compiler assumptions. Proving the actual C++/Python path would require modelling JSON registration, backend dispatch, EOS state updates, Horner rounding, transcendental-library error, overflow and exception behavior. Short arrays and zero normalizers are outside the valid-input theorem and remain unchecked implementation hazards. Neither such a proof nor passing numerical checks would establish empirical accuracy of the correlation or physical validity of arbitrary crossover edits.
