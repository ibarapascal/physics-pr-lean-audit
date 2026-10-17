## CoolProp/CoolProp#2851 fix(HEOS): correct first_saturation_deriv for mixtures (#2091)

- What changed (1–3 sentences, naming files and formulas):
  `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp` changes the mixture saturation tangent from the pure-fluid formula `dT/dp = T(vV−vL)/(hV−hL)` to `−A/B`, where A and B are the differences between the two phases of the partial derivatives of ln(fugacity coefficient) with respect to pressure and temperature, weighted by the vapor-phase composition at the bubble point and by the liquid-phase composition at the dew point. For mixtures only Q=0 or Q=1 is accepted, and the parameterless saturated-phase overload is delegated to the explicit SatL/SatV overloads; the pure-fluid and pseudo-pure branches keep the original formula. `src/Tests/CoolProp-Tests.cpp` adds a bubble-point finite-difference test and an intermediate-Q exception-throwing test.

- Invariants (mathematical statements; write them as checkable equalities/inequalities, give variables, sign conventions, domain of validity):

  Convention: all compositions are mole fractions: overall composition z, liquid x, vapor y, all normalized with positive components; Q is the molar vapor quality. T in K, p in Pa, h in J/mol, ρ in mol/m³, v=1/ρ; phase differences are uniformly vapor minus liquid. For fixed z, let q∈{0,1}, and `p_q(T;z)` be the same continuous vapor–liquid-equilibrium endpoint branch: x=z for q=0, y=z for q=1. Restricted to states where the two phases are distinguishable, the EOS is smooth, the equilibrium solution is locally unique, and the denominators of the relevant derivatives are nonzero; this excludes critical points, phase-branch jumps, and solid or three-phase equilibria.

  1. **Primary — the phase-equilibrium tangent of a mixture.** Let Dq=`first_saturation_deriv(iT,iP)` and Eq=`first_saturation_deriv(iP,iT)`. We should have
     \[
     E_q=\frac{d p_q}{dT}\bigg|_z,\qquad D_q E_q=1.
     \]
     More specifically, let φᵅᵢ(T,p,cᵅ) be the fugacity coefficient of phase α, with the composition of that phase held fixed in the partial derivatives, and define
     \[
     A_q=\sum_i w_i^{(q)}[(\partial_p\ln\phi_i^L)_{T,x}-(\partial_p\ln\phi_i^V)_{T,y}],\quad
     B_q=\sum_i w_i^{(q)}[(\partial_T\ln\phi_i^L)_{p,x}-(\partial_T\ln\phi_i^V)_{p,y}],
     \]
     where w⁽⁰⁾=y and w⁽¹⁾=x. Then **Aq+Bq Dq=0**, equivalently Dq=−Aq/Bq. A has units Pa⁻¹, B has K⁻¹, and D has K/Pa. The bubble-point and dew-point slopes must not be required to be equal, nor may a positive slope be imposed on all mixtures; a positive slope serves only as a local check of the ordinary states designated by this PR.

  2. **Secondary — total derivative of the endpoint properties at fixed composition.** For F∈{hmolar, smolar, rhomolar}, let Fq(T;z) be the F of the primary state after `update(QT_INPUTS,q,T)` (the bulk liquid at the bubble point or the bulk vapor at the dew point). At these endpoints the bulk composition z does not change with T, therefore
     \[
     \operatorname{satderiv}(F,T)=\frac{dF_q}{dT}
       =(\partial_T F)_{p,z}+(\partial_p F)_{T,z}E_q,
     \qquad
     \operatorname{satderiv}(F,p)=D_q\frac{dF_q}{dT}.
     \]
     No claim is made here that the properties of the trace vapor phase at the bubble point or of the trace liquid phase at the dew point satisfy this composition-free formula: the compositions of those phases change along the curve. This invariant covers the actual public-API behavior after the PR's new slope enters the existing chain-rule branch for property derivatives.

  3. **Secondary — Clapeyron consistency of genuinely pure fluids and the single-component reduction.** For a genuinely pure fluid, with the liquid and vapor phases in equilibrium at the same temperature and pressure, Δh=hV−hL, Δv=vV−vL, we should have
     \[
     D_q\Delta h=T\Delta v,\qquad
     E_q=\frac{\Delta h}{T\Delta v},\qquad D_0=D_1.
     \]
     In the ordinary vapor–liquid-equilibrium region Δh>0, Δv>0, so Dq>0. The mixture tangent should reduce to this formula in the smooth single-component limit; usability of the mixture API with zero mole fractions is not enforced. **The retained behavior of pseudo-pure fluids is not directly interpreted as a strict physical theorem**: a pseudo-pure model with temperature glide does not necessarily represent genuine single-component two-phase equilibrium at the same temperature and pressure.

  API scope addendum (not counted as a physical invariant): after a mixture is successfully updated to 0<Q<1, this method should throw an exception, and must not pass itself off as a bubble-point/dew-point tangent. This does not mean that the curve at fixed intermediate Q does not exist mathematically, only that the API does not implement that derivative.

- Sources of the invariants (first principles / literature formula / consistency with other implementations in the library / dimensions):

  The primary invariant can be derived directly from the phase equilibrium `x_i φ_i^L p = y_i φ_i^V p`, without using the new implementation itself as the criterion of correctness. Differentiate its logarithm: at the bubble point dx=0, multiply by yᵢ and sum; using Σdyᵢ=0 and the Gibbs–Duhem relation of the vapor phase at fixed T,p, Σyᵢ d_composition lnφᵢV=0, one obtains A₀dp+B₀dT=0. At the dew point dy=0; multiplying instead by xᵢ and using the same relation for the liquid phase gives A₁dp+B₁dT=0. This derivation requires the fugacities to come from a thermodynamically consistent mixture model.

  The second item comes from the chain rule at fixed composition; the third comes from dμ=−s dT+v dp, equality of μ in the two phases, and Δh=TΔs. For a single component, ∂p lnφ=v/(RT)−1/p and ∂T lnφ=−hᵣ/(RT²), hence A=−Δv/(RT) and B=Δh/(RT²), which gives the same Clapeyron formula.

  The input PR attributes the mixture formula to Eqs. 3.96/3.97 of Gernert's dissertation; this audit did not independently obtain the original thesis and does not treat these equation numbers as verified external evidence. The fixed-version sources actually read are the [HEOS implementation](https://raw.githubusercontent.com/CoolProp/CoolProp/54a2e32764bfa0ae81877e835c33ccd13ba8f853/src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp), the [public AbstractState interface](https://raw.githubusercontent.com/CoolProp/CoolProp/54a2e32764bfa0ae81877e835c33ccd13ba8f853/include/AbstractState.h) and the [Python interface](https://raw.githubusercontent.com/CoolProp/CoolProp/54a2e32764bfa0ae81877e835c33ccd13ba8f853/wrappers/Python/CoolProp/AbstractState.pyx). Other independent audit tracks, comments, issues or forks were not read.

- Common operating domain (range of each input, with the basis):

  First round: `HEOS`, `Methane&Ethane`, z=(a,1−a), a∈[0.6,0.8], T∈[170,190] K, with q taking 0 and 1 respectively. The baseline is the a=0.7, T=180 K, q=0 explicitly tested by the PR; the interval is a conservative sampling suggestion around that point, **not yet verified by running, and not the whole usable range guaranteed by the literature**. p is a QT-flash output rather than an independent input, with an initial screen of 10³≤p≤10⁷ Pa; only samples are admitted for which the center and the finite-difference neighbors all converge and stay on the same equilibrium branch.

  Additionally require all phase fractions to be positive and normalized, ρL>ρV>0, and (ρL−ρV)/ρL>0.05 to avoid phase-density merging; check the consistency of the fugacities of the two phases, rather than deciding that the equilibrium is correct merely because the flash returned. This density gate is only a numerical screen, not a strict sufficient criterion for the critical region. Record failed/unconverged samples, and do not count them as passing the invariant. The intermediate-Q exception check uses the original test point T=180 K, a=0.7, Q=0.5, and may additionally sample Q∈[0.1,0.9].

  The pure-fluid regression uses `HEOS`, `n-Propane`, T∈[280,320] K, q∈{0,1}, centered on the 300 K of the existing regression test; the same equilibrium and non-degeneracy conditions are required. A pure component is not a call that numerically sets a binary z to zero.

- Numerical check method (which public API of the merged commit, what to compute, what to compare against, suggested tolerance and why; if finite differences are needed, the step-size strategy):

  Use the C++ `AbstractState::factory("HEOS","Methane&Ethane")`, `set_mole_fractions({a,1-a})`, `update(QT_INPUTS,q,T)`, `p()`, `first_saturation_deriv`, `keyed_output` of the merge commit. The Python AbstractState wrapper of the same name can also implement the checks. Use independent states for the center and the neighboring points, or reset the state each time, and save the analytic center values before moving the flash state. All finite differences keep **z and q** fixed, and do not lock the composition of the trace phase.

  **Primary check:** compute C(h)=[p_q(T+h)−p_q(T−h)]/(2h), and compare Eq with C, and Dq·C with 1; also check Dq·Eq≈1, but this algebraic check cannot replace the independent flash differences. Use a ladder of h=0.1, 0.05, 0.025, 0.0125 K and look for second-order convergence in the smooth region; [4C(h/2)−C(h)]/3 can be used for extrapolation. It is suggested that, once the differences are stable to a relative 10⁻⁴, the analytic value be required to have relative error ≤10⁻³ (the original PR requires only 5×10⁻³). If the error grows when h is reduced further, adopt the stable plateau; if there is no plateau, mark as numerically undecided, and do not automatically relax to a pass. The reciprocal algebraic check can use |DqEq−1|≤10⁻¹⁰. Samples near zero slope, turning points or denominator degeneracy are moved out of this non-singular domain and recorded separately.

  **Optional cross-evidence for the primary check:** obtain x, y from `mole_fractions_liquid()` and `mole_fractions_vapor()`, create a single-phase state for each phase with its own composition fixed, keep the phase root with `specify_phase(iphase_liquid/iphase_gas)` and `update(PT_INPUTS,p,T)`; construct A and B from central differences of the logarithm of `fugacity_coefficient(i)`. Start with a temperature step of 0.05 K and a pressure step of max(1 Pa,10⁻⁴p), and repeatedly halve each to look for a plateau; these partial derivatives must not be differenced along the QT curve. First require that the reconstructed phase densities agree with the saturated-phase densities to about 10⁻⁶ relative and that maxᵢ|ln(xᵢφᵢL/(yᵢφᵢV))|≤10⁻⁷; otherwise it is judged a failure of reference-state construction. Check |A+BDq|/(|A|+|BDq|)≤10⁻³, with a nonzero denominator. This path cross-validates the fixed-composition fugacity derivatives and the weights, while the primary criterion remains the tangent of the independent flash curve.

  **Endpoint property check:** from the same group of QT neighbors read hmolar, smolar, rhomolar, construct CF(h)=[Fq(T+h)−Fq(T−h)]/(2h), and compare `first_saturation_deriv(iF,iT)` with CF and `first_saturation_deriv(iF,iP)` with CF/C. Use the same step plateau and a 10⁻³ relative tolerance. If a derivative is near zero, use the mixed error bound 10⁻³·max(|reference derivative|,10⁻⁶SF), where for the temperature derivative SF=max(|F center|,Funit)/T, with Funit being 1 J/mol, 1 J/(mol K), 1 mol/m³ respectively; for the pressure derivative the scale is changed to SF/|C|. These are numerical tolerances to be verified, not claims of experimental property accuracy.

  **Pure-fluid check:** construct Δh and Δv through `saturated_liquid_keyed_output(iHmolar/iDmolar)` and `saturated_vapor_keyed_output(...)`, and compare |DΔh−TΔv|/max(|DΔh|,|TΔv|)≤10⁻¹⁰; also make independent QT differences, with a suggested relative tolerance of 10⁻⁴, to guard against checking only the implementation's own original formula. Compare the Q=0 and Q=1 tangents at the same T with a relative 10⁻⁸. For intermediate Q, additionally check that the call really throws an exception. All suggestions are specifications for later execution; no code was run and no numerical results were produced this time.

- Do the PR's own tests already check this invariant (yes/partial/no, naming test files and assertions):

  **Primary: partial.** `[first_saturation_deriv][2091]` at `src/Tests/CoolProp-Tests.cpp:4815–4840` in the merged version checks only Methane/Ethane=(0.7,0.3), T=180 K, Q=0, step 0.05 K; it asserts D>0 and a relative error <5e−3 against the finite difference of the inverted pressure. There is no dew point, composition sweep, step-size convergence, fugacity residual, or direct Eq check.

  **Endpoint properties: partial.** `[first_saturation_partial_deriv]` at `:2411–2449` of the same file checks ten temperature/pressure derivatives of P, T, Dmolar, Hmolar, Smolar, Umolar for pure n-Propane, 300 K, Q=1, requiring |numerical/analytical−1|<1e−4; it does not cover the property derivatives of the new mixture branch.

  **Pure-fluid Clapeyron: partial.** The old test above indirectly checks the saturation tangent of a pure fluid, but does not assert the latent-heat relation directly, nor compare Q=0 with Q=1. `:4843–4849` uses `CHECK_THROWS` for `first_saturation_deriv(iT,iP)` at Q=0.5, fully covering the newly added exception behavior for that single state. The PR description claims the related regressions pass; this audit only checked the test source and did not re-run it.

- Do you expect the change to satisfy it (satisfies / may not satisfy / uncertain), with one or two sentences of reasoning (speculation, not a conclusion):

  **satisfies (speculation, limited to the non-degenerate domain above).** The new opposite-side composition weights and negative sign are consistent with the phase-equilibrium differentiation and the Gibbs–Duhem elimination, the fixed bulk composition at the endpoints also supports the retained chain rule, and the pure-fluid formula is unchanged. The dew point and the mixture property derivatives lack direct tests in this PR and therefore still need numerical verification; strict physical consistency in singular regions and for pseudo-pure fluids remains uncertain.

- Environment: merge commit SHA, language and dependencies needed, estimated setup difficulty (low/medium/high):

  SHA: `54a2e32764bfa0ae81877e835c33ccd13ba8f853`. The C++ public API is suggested as the baseline; [CMakeLists.txt](https://raw.githubusercontent.com/CoolProp/CoolProp/54a2e32764bfa0ae81877e835c33ccd13ba8f853/CMakeLists.txt) specifies C++17, CMake≥3.14, with dependencies managed through CPM. [cmake/dependencies.cmake](https://raw.githubusercontent.com/CoolProp/CoolProp/54a2e32764bfa0ae81877e835c33ccd13ba8f853/cmake/dependencies.cmake) lists Eigen 5.0.1, fmt 12.0.0, msgpack, rapidjson, IF97, REFPROP headers, boost headers, multicomplex; running the native tests also needs Catch2 3.8.0. Only HEOS is used here, and the proprietary REFPROP runtime library is not needed. Python can serve as an optional driver, but must be bound to this SHA; the Python package build dependencies were not completely verified, and an arbitrary installed version must not be substituted. Estimated setup difficulty **medium**: a local build of the exact version and dependency resolution are needed. This time no project code was installed, compiled or run, and the repository was not cloned.

- Formalization (Lean) feasibility: can the invariant be written as a proposition over finite-dimensional reals, and what are the main obstacles:

  **Yes.** Fix a finite number of components N, take the positive reals T,p and the simplex-interior x,y as variables, and represent the EOS by smooth real functions lnφᵢᵅ and a local equilibrium path; assuming equal fugacities, normalization, Gibbs–Duhem, fixed bulk composition and nonzero denominators, one can formally derive A dp+B dT=0, the inverse-function slope, the property chain rule and the single-component reduction. The main obstacles are proving these thermodynamic-consistency premises from the actual multi-parameter Helmholtz EOS, proving by the implicit function theorem that the local equilibrium solution exists and is unique, handling critical points and phase selection, and connecting the real-number proof to floating-point computation/flash-solver errors; proving the algebraic formulas does not prove that the solver finds the physically stable phase. Throwing an exception at intermediate Q is a program contract, requiring a separate proof of program semantics, and is not a real-number thermodynamic proposition.
