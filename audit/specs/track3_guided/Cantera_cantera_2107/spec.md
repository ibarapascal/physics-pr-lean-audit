## 1. What changed

This specification concerns Cantera PR #2107, **Fix thermo consistency errors**, merged as `d58a23265414e1789536cf5332d9fb29de095a43`, relative to `a6feb655151a3740f83703be87522322855f4a75`. It is based on the supplied diff and the merged source, not on execution of Cantera. Quantities below use Cantera's molar units: J/kmol, J/(kmol K), m³/kmol, Pa and K; (R) has corresponding units. Composition is held fixed unless stated otherwise.

- `src/thermo/PDSS.cpp` and `include/cantera/thermo/PDSS*.h`: introduce public virtual `dVdT()` and `dVdP()`. The base implementations throw `NotImplementedError`. `PDSS_ConstVol.cpp` returns zero for both and changes `cv_mole()` to `cp_mole()`. `PDSS_SSVol.cpp` returns its cached temperature derivative and zero pressure derivative, also sets (c_v=c_p), and fixes both `setPressure` and `setTemperature` to use (c_p^\circ/R=c_{p,r}^\circ/R-T(P-P_r)V^{\circ\prime\prime}(T)/R).
- `src/thermo/PDSS_Water.cpp`: adds (V_T^\circ=\beta V^\circ\) and (V_P^\circ=-\kappa_T V^\circ\). `PDSS_HKFT.cpp` factors volume evaluation into `molarVolumeTP(T,P)` and adds centered differences with steps (10^{-4}T) and (10^{-4}P), restoring the shared water state afterward. **Despite the PR body's wording, the merged `PDSS_HKFT` does not override `cv_mole()`; direct calls still reach the throwing base implementation.**
- `src/thermo/IdealSolnGasVPSS.cpp`: replaces (c_p-R) by (c_p+T V_T^2/V_P), using composition-weighted PDSS derivatives, and adds `thermalExpansionCoeff()` and `isothermalCompressibility()`. The exactly zero (V_P) branch returns (c_p).
- `src/thermo/MargulesVPSSTP.cpp` and its header: add the analogous mixture mechanical properties and heat-capacity formula, including the temperature derivative of excess volume. For interaction (i=(A_i,B_i)), define
  \[
  q_{ji}=h_{ji}-Ts_{ji}+(P-P_*)[v^H_{ji}-Tv^S_{ji}],\qquad j=0,1,\quad P_*=101325\ {\rm Pa}.
  \]
  The molar excess Gibbs energy is (g^E=\sum_i X_{A_i}X_{B_i}(q_{0i}+q_{1i}X_{B_i})). These pressure terms are now included in `s_update_lnActCoeff`, `s_update_dlnActCoeff_dT`, `getdlnActCoeffds`, and the three updates for `dlnN`, `dlnN_diag`, and `dlnX_diag`. The volume contribution is (v^E=\sum_i X_{A_i}X_{B_i}[(v^H_{0i}-Tv^S_{0i})+(v^H_{1i}-Tv^S_{1i})X_{B_i}]).
- `src/thermo/GibbsExcessVPSSTP.cpp` and its header: remove the density override based only on standard volumes. The inherited `VPStandardStateTP::calcDensity()` uses partial molar volumes: (\rho=\overline M/\sum_kX_k\bar v_k). This includes excess volume as well as current standard-state pressure dependence.
- `src/thermo/MolalityVPSSTP.cpp` and its header: `getActivityCoefficients()` now divides solute molality coefficients by (\max(X_0,X_{0,\min})\widetilde M_0), rather than just the solvent mole fraction. Here (\widetilde M_0=M_0/1000) is the solvent molecular weight in kg/mol. The solvent coefficient is unchanged.
- `src/thermo/IdealMolalSoln.cpp`: `getPartialMolarEntropies()` multiplies the standard (s^\circ/R) values by (R) before applying the activity correction, in both cutoff and no-cutoff branches.
- `src/thermo/WaterProps.cpp`: `ADebye(...,2)` now approximates (\partial_T\beta\) by ([\beta(T+0.1,P)-\beta(T-0.1,P)]/0.2), replacing an erroneous division by (T+\Delta T). `src/thermo/DebyeHuckel.cpp`: in the dilute-limit branch, solvent derivatives become (\partial_P\ln\gamma_0^\triangle=(2/3)\widetilde M_0 I^{3/2}A_P) and (\partial_T^2\ln\gamma_0^\triangle=(2/3)\widetilde M_0 I^{3/2}A_{TT}); solutes retain the proportional-to-(A) formula.
- `src/thermo/BinarySolutionTabulatedThermo.cpp`: `_updateThermo()` replaces endpoint entropy sentinels by (R\ln[\max(10^{-300},x)/\max(10^{-300},1-x)]), retaining interpolation and the existing (R/F\) standard-concentration correction. The header explains the vacancy-reference convention and entropy divergence.
- Tests add Margules excess-volume and polynomial-volume cases, remove known-failure entries, add a (c_p-c_v) identity test, skip constant-volume differentiation for non-compressible API phases, adjust tolerances and reference outputs. Build scripts add test RPATH/data paths; these are execution changes, not physical laws.

## 2. Invariants

**I1 — Primary: thermodynamic potential and caloric consistency.** For a smooth, single-phase model at fixed composition,
\[
 g=h-Ts,\quad h=u+Pv,\quad \mu_k=\bar h_k-T\bar s_k,\quad
 c_p=h_T=T s_T,\quad s_P=-v_T,\quad \bar v_k=(\partial_P\mu_k)_{T,X}.
\]
Bulk extensive molar properties obey (f=\sum_kX_k\bar f_k), including (v), (h), (s), (g), and (c_p). Standard-state analogues apply species by species. For SSVol specifically, integration of (g_P^\circ=V^\circ(T)) gives
\[
 g^\circ=g_r^\circ+\Delta P V^\circ,\quad
 s^\circ=s_r^\circ-\Delta P V_T^\circ,\quad
 h^\circ=h_r^\circ+\Delta P(V^\circ-TV_T^\circ),\quad
 c_p^\circ=c_{p,r}^\circ-T\Delta P V_{TT}^\circ.
\]
These statements exclude phase transitions, polynomial/interpolation derivative discontinuities and clipping boundaries. For the binary tabulated convention, apply the algebraic and fixed-composition thermal identities; do not assume that its composition-dependent “standard states” define a conventional globally integrable mixture Gibbs energy.

**I2 — Primary: mechanical and heat-capacity consistency.** Let (v=\overline M/\rho>0). Then
\[
 \beta=v_T/v,\quad\kappa_T=-v_P/v,\quad
 c_p-c_v=-T v_T^2/v_P=T v\beta^2/\kappa_T.
\]
The last identity requires (v_P\ne0) and a differentiable local change of coordinates between ((T,P)) and ((T,v)). On stable branches with (\kappa_T>0), it implies (c_p\ge c_v). Do not require stability or positive individual HKFT solute volumes for arbitrary fitted species. For truly rigid (v_T=v_P=0), (c_v=c_p). For (v_P=0) but (v_T\ne0), a varying-temperature isochore is unavailable: the implementation's (c_v=c_p) is a **declared library convention**, not a consequence or continuous limit of the identity. This matters for SSVol and temperature-dependent Margules excess volume.

Each supported PDSS must return the actual derivatives of its own volume. Evaluating these derivatives must leave observable parent/water thermodynamic state unchanged. The direct HKFT (c_v) method is not part of the implemented guarantee.

**I3 — Primary: activity-basis consistency.** With the library's molal standard state, dimensionless activity uses the numerical molality in mol/kg relative to (m^\circ=1\) mol/kg:
\[
 a_k=X_k\gamma_k=(m_k/m^\circ)\gamma_k^\triangle,\quad
 m_k=X_k/(X_0\widetilde M_0),\quad
 \mu_k=\mu_k^\circ+RT\ln a_k\quad(k>0).
\]
For solvent, (a_0=X_0\gamma_0^\triangle\) and (\gamma_0=\gamma_0^\triangle). Thus the solute coefficient includes the numerical (1/(X_0\widetilde M_0)) conversion. These coefficients need not approach one under mole-fraction dilution: the chemical-potential reference has not been changed to a Raoult-law standard state. Require these equalities at positive composition above solvent clipping thresholds; below those thresholds test the documented regularization separately. At zero activity use the limiting/clamped chemical-potential rule rather than an unqualified logarithm of zero.

**I4 — Margules excess potential and derivative consistency.** With (g^E) in section 1, (RT\ln\gamma_k=\partial(ng^E)/\partial n_k), (v^E=g_P^E), (s^E=-g_T^E), (h^E=g^E+Ts^E), and (c_p^E=0) for temperature-independent interaction coefficients. The public derivative functions must differentiate this same activity function. In particular, for a normalized path,
\[
 D_k=\frac{d\ln\gamma_k(T(s),X(s),P)}{ds},\qquad\sum_jX'_j=0,
\]
and (J_{kj}=\partial\ln\gamma_k/\partial\ln n_j\) with other mole numbers fixed. Its diagonal must agree with `getdlnActCoeffdlnN_diag`; Gibbs–Duhem requires (\sum_kX_kD_k=0) for an isothermal, isobaric composition path. Domain: any interior composition, positive mixture volume, any finite collection of binary interactions. The `dlnX_diag` documentation does not precisely fix which other fractions are dependent. Specify a coordinate convention before treating a discrepancy in that API as a physical violation; test the unambiguous normalized-path and mole-number APIs independently.

**I5 — Debye-water derivative consistency.** On the liquid branch with imposed (P>p_{sat}(T)), write (A=C\rho_w^{1/2}(\epsilon_rT)^{-3/2}). Then
\[
 L=\partial_T\ln A=-\tfrac12\beta_w-\tfrac32/T-\tfrac32\epsilon_{r,T}/\epsilon_r,
 \quad A_{TT}=A(L^2+L_T),\quad
 A_P=A[\tfrac12\kappa_w-\tfrac32\epsilon_{r,P}/\epsilon_r].
\]
At fixed molality and below any ionic-strength cap, dilute-limit solutes have (\ln\gamma_k^\triangle=-z_k^2 A\sqrt I); solvent has (\ln\gamma_0^\triangle=C_0(X)+(2/3)\widetilde M_0 A I^{3/2}). Only the (A)-dependent term contributes to the two derivatives changed here. Consequently (\bar v_k-v_k^\circ=RT\partial_P\ln\gamma_k^\triangle) and
\(\bar c_{p,k}-c_{p,k}^\circ=-R[2T\partial_T\ln\gamma_k^\triangle+T^2\partial_T^2\ln\gamma_k^\triangle]\).
`ADebye` evaluates water at (max(P,p_{sat})). Its derivative options below saturation are local fixed-effective-pressure quantities, not necessarily derivatives of the composed, pressure-clamped value function; the distinction must not be hidden in a finite-difference test.

**I6 — Binary endpoint consistency under the vacancy convention.** With zero reference/background (h,s,c_p), equal standard concentrations, and reference pressure, the implemented reaction difference must satisfy
\[
 \mu_{tab}-\mu_{ref}=h_{tab}(x)-Ts_{tab}(x).
\]
It must remain finite and preserve I1's algebraic identities at (x=0,1) when both standard-state and activity logarithms use the same positive floor. For (0<x<1), under these simplifying conditions,
(s=x s_{tab}-R\ln(1-x)); thus divergence as (x\to1\) is expected before the numerical floor cuts it off. Do not demand a finite physical entropy limit. With nonzero species backgrounds, unequal standard concentrations or nonreference pressure, include their explicit chemical-potential contributions; the header's simple reaction equality is not unconditional. In particular the retained (R/F\) entropy correction is dimensionally questionable in SI, and its general interpretation is not settled by the endpoint fix.

## 3. Sources of the invariants

I1 and I2 follow from (dg=-s\,dT+v\,dP), extensivity, equality of mixed derivatives, and the chain rule, with a nonsingular equation of state for (c_v). The zero-compressibility fallback comes from the changed class contracts, not a universal thermodynamic theorem.

I3 follows directly from the definition of molality and the library-wide activity interface. The [merged MolalityVPSSTP header](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/include/cantera/thermo/MolalityVPSSTP.h) cites Denbigh (1981), p. 278; that external text was not independently consulted. I4 follows by differentiation of the potential in the [Margules header](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/include/cantera/thermo/MargulesVPSSTP.h), rather than empirical accuracy of a parameter fit. The new ethanol/water data explicitly disclaim scientific accuracy.

I5 is a calculus consistency requirement for the dielectric/water expression actually evaluated by [WaterProps.cpp](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/src/thermo/WaterProps.cpp) and the dilute Debye–Hückel expression in [DebyeHuckel.cpp](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/src/thermo/DebyeHuckel.cpp). It does not independently certify the underlying IAPWS or dielectric fits. I6 is specifically the [BinarySolutionTabulatedThermo contract](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/include/cantera/thermo/BinarySolutionTabulatedThermo.h) together with the actual `_updateThermo` implementation. No external standard resolves its unusual reference convention in the material reviewed.

## 4. Operating domain

Use ordinary single-phase states, and separate empirical validity from mathematical stress testing:

| Inputs/model | Main check range and restrictions |
|---|---|
| Water, HKFT, aqueous Debye–Hückel/HMW | (T=280\)–(360\) K, (P=0.101325\)–(10\) MPa, with pressure safely above saturation for the whole stencil. Extend to 500 K only at suitable liquid pressure and within species validity. Avoid critical/saturation boundaries, the HKFT pole at 228 K, and (P=0) where the HKFT difference step vanishes. |
| Molality and charge | Charge-balanced NaCl from (10^{-6}\) to (10^{-2}\) mol/kg for physically dilute DH checks; (0.1,1,6,9.4\) mol/kg only as formal model/regression stress cases. Include neutral solutes and monovalent ions, plus a valid multivalent fixture. (I=\tfrac12\sum m_kz_k^2\); stay below the configured cap for derivative invariants. |
| Composition | Nonnegative normalized fractions. Use interior fractions (10^{-6}\)–(1-10^{-6}\), pure vertices for defined algebraic quantities, solvent-rich aqueous states, and cutoff-adjacent cases separately. Solvent molar mass is positive; test water and an artificial different-mass solvent for the generic conversion. |
| Margules/SSVol | Ethanol/water fixture: 280–330 K; Li polynomial fixture: 300–500 K; molten LiCl/KCl: 900–1400 K. Pressure 1–20 atm plus elevated-pressure synthetic cases. Require (v>0\) and polynomial density nonzero. Test every coefficient type (h_0,h_1,s_0,s_1,v^H_0,v^H_1,v^S_0,v^S_1\), including zero values and both signs where admissible. Typical liquid molar volumes are (10^{-2}\)–(10^{-1}\) m³/kmol; artificial fixtures need not be realistic materials. |
| Binary tables | 280–350 K, 1–10 atm; (x\in[0,1]\). Finite enthalpy/entropy entries with ordered distinct interpolation nodes and positive volume. Include interior nodes, between-node states, exact endpoints, and trace fractions down to the logarithmic floor. Avoid crossing interpolation knots in derivative checks. |

All coefficient sets, standard-concentration choices, species reference thermochemistry and cutoff parameters are inputs to record. No universal admissible range for arbitrary fitted coefficients is claimed.

## 5. Numerical check plan

Use a C++ harness linked to the merged library, so all relevant public APIs are reachable. Construct phases with `newThermo(file, phase)`, add `data` and `test/data` via `addDataDirectory`, and set normalized composition followed by `setState_TP`. Use `MolalityVPSSTP::setMolalities` for aqueous cases. Never equilibrate during differentiation. Read `enthalpy_mole`, `intEnergy_mole`, `entropy_mole`, `gibbs_mole`, `cp_mole`, `cv_mole`, `molarVolume`, `density`, `getChemPotentials`, the partial-molar getters, and the standard-state getters as needed. Cast to `VPStandardStateTP` for `providePDSS(k)` and to `MargulesVPSSTP` for its derivative methods.

Use centered first differences at (h_T/T=10^{-4},5\times10^{-5},2.5\times10^{-5}\) and (h_P/P=10^{-3},5\times10^{-4},2.5\times10^{-4}\), retaining only stencils inside one regime. Compare successive estimates before assessing a residual. For ordinary algebraic identities use (|a-b|\le10^{-5}+10^{-12}\max(|a|,|b|)) in energy units. Derivative tolerances: relative (10^{-6}\) for analytic models, (10^{-5}\) for water/HKFT; absolute floors (10^{-3}\) J/(kmol K) for heat capacity, (10^{-10}\) m³/kmol for volume, (10^{-9}\) m³/(kmol K) for Maxwell/volume-temperature derivatives, (5\times10^{-14}\) Pa⁻¹ for compressibility and (10^{-9}\) K⁻¹ for expansion. Near-zero responses need an absolute floor; a nonconvergent stencil is inconclusive, not a pass. These are proposed acceptance tolerances, not measured results.

1. **I1:** Compare all algebraic identities and bulk/partial sums. Differentiate (h,s,\mu_k\) independently. For SSVol, use both `molar-volume-temperature-polynomial` and `density-temperature-polynomial`, with nonzero quadratic/cubic curvature, at (P_r\) and well away from it. Compare to the analytic SSVol equations as well as finite differences. Exercise `setTemperature` then `setPressure`, the reverse order, and `setState_TP`, comparing final properties. For ideal molal solution run `setCutoffModel("none")`, `"poly"`, `"polyexp"`; test solvent-rich states and states on either side of each configured transition, without straddling it. For GibbsExcess descendants include Margules with/without excess volume and the existing Redlich–Kister fixtures.
2. **I2:** Include all five registered PDSS models: constant volume, both SSVol polynomials, liquid-water-IAPWS95, and HKFT. For HKFT include charged and neutral solutes with their required water parent. Check `dVdT`/`dVdP` against externally differenced `molarVolume`; snapshot parent and water properties before/after repeated derivative calls. In both ideal VPSS and Margules, use rigid-only, SSVol-only, water plus rigid solute, water plus HKFT, and compressible mixtures with excess-volume coefficients. Compare (\rho\), partial-volume sums, (\beta\), (\kappa_T\) and (c_p-c_v\) against independently differenced volume. Exercise exact-zero and nonzero compressibility branches. Since VPSS reports `isCompressible()==false`, do not rely on `setState_TD`: when (v_P\ne0\), solve (v(T\pm h,P_\pm)=v(T,P)\) by a bracketed scalar pressure root using `setState_TP`, then compare `cv_mole` to the isochoric derivatives of (u\) and (s\). Require the root's residual below (10^{-11}\) relative volume; inability to bracket within one liquid branch is inconclusive. Test the (v_P=0,v_T\ne0\) fallback only as a convention. A direct HKFT `cv_mole` call must be recorded as unsupported.
3. **I3:** Use `getActivities`, `getActivityCoefficients`, `getMolalityActivityCoefficients`, `getMolalities`, and `getStandardChemPotentials`. Test (a=X\gamma\) at (10^{-10}\) relative plus (10^{-20}\) absolute, and chemical-potential differences in units of (RT\) at (10^{-9}\) absolute. Cover ideal molal solution, HMW, and all ten DH fixtures: dilute-limit, B-dot-with-variable-a, B-dot-with-common-a, Pitzer-with-beta_ij, beta_ij, each with constant and water-dependent (A\). Repeat standard-concentration choices `unity`, `species-molar-volume`, and `solvent-molar-volume` where supported. Include trace solutes and pure solvent; treat solvent-floor cases separately from the unclipped invariant.
4. **I4:** Compute (g^E,h^E,s^E,v^E\) from independent polynomial expressions and compare to phase properties minus standard-state/ideal-mixing terms. Test (P=P_*\), (2P_*\), (20P_*\), and 10 MPa with admissible coefficients, so pressure terms cannot disappear. Use zero, one, and at least two interactions in a ternary phase, including a species outside a particular interacting pair, with both (q_0\) and (q_1\) nonzero. Test `getdlnActCoeffdT`, `getd2lnActCoeffdT2`, `getdlnActCoeffds`, `getdlnActCoeffdlnN`, and its diagonal. For (J\), perturb (n_j\mapsto n_je^{\pm\epsilon}\), renormalize, and difference `getActivityCoefficients` in log space; use (epsilon=10^{-4},10^{-5},10^{-6}\). Check normalized composition directions, pure-temperature directions, and combined directions, including multiple interactions. Target (10^{-6}\) relative plus (10^{-8}\) absolute for dimensionless derivatives. For `dlnX_diag`, explicitly report both an independent-coordinate polynomial extension and a binary normalized path; do not silently substitute a mole-number derivative for this ambiguous interface.
5. **I5:** Directly use public `WaterProps::ADebye(T,P,ifunc)` for `ifunc=0,1,2,3`. Compare mode 2 with the temperature derivative of mode 1 and the second difference of mode 0, trying 0.05, 0.1 and 0.2 K steps; require (10^{-4}\) relative plus (10^{-9}\) in numerical (A\)-units/K². Modes 1 and 3 use first-difference tolerances (10^{-5}\) relative with scaled absolute floors. Also compare DH solvent and each solute partial (c_p\) and volume with derivatives of their chemical potentials/enthalpies, not just bulk averages. Constant-(A\) controls must have zero (A\)-derivative contributions. Repeat the caloric/Maxwell checks for `HMW_NaCl.yaml`. Probe below-saturation and saturation-clamp cases as convention tests, explicitly excluding them from the constant-imposed-pressure derivative claim.
6. **I6:** Use `BinarySolutionTabulatedThermo.yaml`, `anode`, and a simplified zero-background, equal-standard-concentration fixture. Test both species orderings/tabulated indices. Sample every table interval midpoint and node, (x=0,1\), (10^{-2},10^{-6},10^{-10},10^{-14},10^{-300}\), and symmetric trace-reference states. Set fractions explicitly to avoid losing tiny (1-x\) by subtraction. Check finite outputs and I1's algebraic sums, then the reaction difference against independently interpolated data, with energy tolerance above. For unequal standard concentrations, nonzero background thermo, and elevated pressure, verify the explicit implemented corrections separately and flag the dimensional issue rather than assume the simplified equality. Revisit the same composition after temperature changes to check cache consistency.

## 6. Do the PR's own tests check these invariants

The following names are `TestConsistency` tests in `test/thermo_consistency/consistency.cpp`; configurations are in `test/data/consistency-cases.yaml`. Coverage below describes assertions, not a locally observed passing run.

| Invariant | Coverage | Actual assertions and gaps |
|---|---|---|
| I1 | **Yes for core identities; partial for full domain** | `g_eq_h_minus_Ts`, `gk_eq_hk_minus_T_sk`, `gk0_eq_hk0_minus_T_sk0`, `h_eq_u_plus_Pv`, `cp_eq_dhdT`, `cp_eq_dsdT_const_p_times_T`, `cpk0_eq_dhk0dT`, `dsdP_const_T_eq_minus_dV_dT_const_P`, and bulk/partial sum tests compare the stated quantities. Ideal-molal known failures are removed; `MargulesSSVol` exercises both polynomial models. No broad cutoff/state-order matrix or individual μ-pressure derivative matrix is added. The updated IMS blessed output is regression evidence only. |
| I2 | **Partial** | New `cv_eq_cp_minus_T_V_beta2_over_kappa_T` checks the algebraic formula, or (c_v=c_p\) at exactly zero compressibility. Existing `betaT_eq_minus_dmv_dP_const_T_div_mv`, `alphaV_eq_dmv_dT_const_P_div_mv`, and `v_eq_sum_vk_Xk` compare mechanical derivatives/volume sums. `cv_eq_dudT` and `cv_eq_dsdT_const_v_times_T` now skip non-compressible API phases, including VPSS even when their physical volumes depend on pressure. Thus the new test can share the implementation's errors; no independent isochoric root check or direct PDSS restoration check is added. |
| I3 | **Yes at supplied states; partial at cutoffs/options** | `activity_coeffs` asserts (a_k=X_k\gamma_k\) for every species, tolerance (10^{-10}\) relative plus (10^{-20}\) absolute; `chem_potentials_to_activities` compares exponentiated potentials to activities. Enabled across the molality cases. `ThermoFromYaml.HMWSoln_HKFT` updates coefficient regression values (solute factor about 55.5), with absolute coefficient tolerance (10^{-6}\). No different-solvent-mass or full cutoff matrix. |
| I4 | **Partial** | New `MargulesWithExcessVolume` checks the general thermodynamic identities at two compositions and two temperatures, all centered at 1 atm; finite pressure perturbations can detect missing pressure derivatives, but not broadly validate finite-pressure values. The added fixture has one binary interaction and (v^S_1=0\). No new assertions for the public activity derivative matrix/path functions or multiple interactions. |
| I5 | **Partial** | Enabled DH variable-(A\) and HMW `cp_eq_dhdT`, `cp_eq_dsdT_const_p_times_T`, and Maxwell tests indirectly exercise the correction. DH IAPWS cases use (10^{-5}\) relative finite-difference tolerance; HMW uses (4\times10^{-6}\). Updated `HMW_graph_CpvT/output_blessed.txt` provides regressions, not a direct independent test of all `ADebye` derivative modes or the saturation clamp. |
| I6 | **Partial** | Binary cases include both exact endpoints and trace fractions (10^{-10}\); removed failures enable `g_eq_h_minus_Ts`, `g_eq_sum_gk_Xk`, `s_eq_sum_sk_Xk`, `gk_eq_hk_minus_T_sk`, `gk0_eq_hk0_minus_T_sk0`, `chem_potentials_to_activities`, and `gRef_eq_hRef_minus_T_sRef`. No dedicated assertion of the interpolated reaction driving force, species-order symmetry or general standard-concentration correction. |

Default finite-difference relative tolerance is (10^{-6}\). HKFT cases use (4\times10^{-6}\) and absolute compressibility tolerance (5\times10^{-14}\). Tests catching `NotImplementedError` report skips; removing a known-failure entry does not turn an unsupported property into a tested one.

## 7. Expected verdict

- **I1: expected to pass in the stated smooth domain**, especially the missing-(R\) entropy and SSVol heat-capacity corrections. Broad cutoff regimes and extreme numerical states remain unverified.
- **I2: expected to pass for the implemented mixture methods**, subject to finite-difference/solver tolerance and excluding singular states. The zero-compressibility fallback passes its stated convention but is not a physically derived (c_v\) when thermal expansion is nonzero. Any claim that this PR implements direct `PDSS_HKFT::cv_mole()` is false for the merged code.
- **I3: expected to pass above the solvent floor**. The added molecular-weight factor follows exactly from the selected standard state. Below-floor behavior is uncertain until each derived activity implementation is checked.
- **I4: expected to pass for the principal excess potential/volume changes, but not all derivative APIs.** Source inspection shows `getdlnActCoeffds` adding the *already summed* temperature-derivative vector inside the loop over interactions. With (N>1\) interactions, a pure-temperature direction appears to return (N\) times the required derivative; this is an expected failure, not a measured runtime result. `s_update_dlnActCoeff_dlnX_diag` also retains the suspicious product `2*g1*-2*g0` rather than a linear combination of interaction energies. Its interpretation needs the coordinate convention resolved. These expressions predate this PR; the pressure-term changes nevertheless make their consistency relevant. Distinguish these retained defects from regressions introduced by the PR.
- **I5: expected to pass away from saturation/critical/clipping boundaries**, with finite-step errors. The new solvent derivatives correctly discard composition-only terms. A naive derivative of the pressure-clamped `ADebye(...,0)` below saturation is not expected to match all derivative options.
- **I6: expected to pass endpoint algebraic consistency and the restricted reaction identity**. General unequal-concentration interpretation remains uncertain; the retained (R/F\) factor is not validated by the endpoint fix. Divergent unregularized entropy is intentional under the chosen convention.

These are source-based expectations; no build, numerical check or test suite was run for this specification.

## 8. Environment

Use the exact merge commit above, including its `data`, `test/data`, and referenced mechanism files; an older packaged release is insufficient because the PDSS derivative APIs are new (annotated “New in Cantera 4.0”). Build the C++ library and harness with the compiler/dependency versions accepted by that checkout's `SConstruct`, Python/SCons for building, Eigen, fmt, yaml-cpp and Sundials (bundled versions or supported system versions). GoogleTest is needed to execute the existing C++ tests. Build with `scons build`; use `scons test` for the repository suite, or the generated thermo-consistency binary from the repository root. The PR adds build-directory RPATH and relative test data paths to support direct binary execution.

The proposed harness needs only public C++ interfaces and a scalar root-finding routine; no GPU, remote machine or heavy simulation is required. Python bindings with NumPy can cover bulk checks, but are not a replacement for direct PDSS, WaterProps and Margules derivative API checks. Record compiler, floating-point mode, actual dependency versions and all skips. Avoid fast-math for the initial reference run.

## 9. Formalization feasibility

I1's algebraic relations, the SSVol corrections, I3's activity conversion, and the polynomial Margules part of I4 are practical Lean targets over real numbers. One can define (G=n g(T,P,n/\sum n)\), prove its partial derivatives and Gibbs–Duhem consistency with nonzero-denominator/positive-composition hypotheses, and prove the two activity representations equal after making the mol/kg versus kmol units explicit. A units-aware formalization would also expose the questionable binary (R/F\) term.

I2 follows from the chain rule and implicit-function theorem assuming (v_P\ne0\); the singular branch needs a separate specification, not division by zero followed by a limit. A transcribed multi-interaction directional-derivative function can be disproved by a simple two-interaction, pure-temperature counterexample if the inspected accumulation behavior is preserved. Ambiguous `dlnX_diag` semantics must be fixed before a meaningful proof can be attempted.

I5 can be proved symbolically assuming smooth water density and dielectric functions. Proving the actual implementation's accuracy is much harder: it requires finite-difference truncation bounds, IAPWS iterative-solver error bounds, saturation-branch conditions and floating-point rounding analysis. The HKFT derivatives similarly approximate rather than exactly equal mathematical derivatives. State restoration adds a mutable-state/caching proof obligation beyond scalar calculus.

I6's floor-based algebra can be proved piecewise, including exact endpoints, but logarithm laws require explicit positivity and consistent floors. Interpolation knots, floating-point composition normalization, reference-state choices and cancellation of large chemical potentials need separate treatment. A proof of the transcribed formulas alone would establish neither empirical validity of fitted thermodynamics nor equivalence to the compiled C++ implementation.
