## Cantera/cantera#2107 Fix thermo consistency errors

- What changed (1–3 sentences, naming files and formulas):

  `src/thermo/IdealSolnGasVPSS.cpp`, `MargulesVPSSTP.cpp` and `PDSS_{ConstVol,SSVol,HKFT,Water}.cpp` add the temperature and pressure derivatives of the volume and correct the phase-level heat capacity to \(c_v=c_p-TV\alpha^2/\kappa_T\) (the zero-compressibility branch is handled separately); `PDSS_SSVol.cpp` changes the pressure correction in the dimensionless heat capacity to \(-T(P-P_r)V''(T)/R\), and the `GibbsExcessVPSSTP::calcDensity` override is deleted so that the density is computed from the partial molar volumes. `MargulesVPSSTP.cpp` includes \((P-P_r)V^E\) in the activities and their derivatives, `MolalityVPSSTP.cpp` corrects the conversion between the two activity-coefficient bases, `IdealMolalSoln.cpp` adds the missing \(R\) to the standard-state entropy, and `DebyeHuckel.cpp` and `WaterProps.cpp` correct the temperature and pressure derivatives of the solvent activity and the central differences for water properties. `BinarySolutionTabulatedThermo.cpp/.h` unifies the endpoint logarithm truncation and explains the reference-state convention; there are also test and build-path changes. The following selects the three most important classes of physical invariants, and does not treat the endpoint convention of the tabulated model or the build path as new physical laws.

- Invariants (mathematical statements; write them as checkable equalities/inequalities, give variables, sign conventions, domain of validity):

  **Common notation.** \(T>0\) is the absolute temperature, \(P>0\) the absolute pressure, \(x_k\geq0,\sum_kx_k=1\); all temperature and pressure partial derivatives are taken at fixed composition, and chemical equilibrium must not be invoked during the perturbations. \(V\) is the molar volume of the phase (m³/kmol), \(\rho\) the mass density, \(\overline W=\sum_kx_kW_k\) (kg/kmol); \(g,h,u,\mu_k\) are in J/kmol, \(s,c_p,c_v\) in J/kmol/K, and \(R\) is `Cantera::GasConstant`. The superscript \(\circ\) denotes the standard state at the current \(T,P\), and an overbar denotes a partial molar quantity. The chemical potential is the chemical part returned by `getChemPotentials`, with the electric potential fixed at zero. Derivatives require a locally single-phase, smooth state that does not cross a truncation boundary.

  **I1 — primary: the heat-capacity difference of a compressible liquid phase must be determined by the actual volume response.** Define
  \[
  V=\overline W/\rho,\qquad V_T=(\partial_TV)_{P,x},\quad V_P=(\partial_PV)_{T,x},\qquad
  \alpha=V_T/V,\quad\kappa_T=-V_P/V.
  \]
  On the mechanically stable branch \(V>0,V_P<0\) of `ideal-solution-VPSS` or `Margules`, we must have
  \[
  c_p-c_v=-T\frac{V_T^2}{V_P}=TV\frac{\alpha^2}{\kappa_T}\geq0,
  \qquad c_v=\left(\frac{\partial u}{\partial T}\right)_{V,x}.
  \]
  The \(c_p-c_v\) of the individual species must not simply be added according to composition: first obtain the mixture \(V_T,V_P\), then take the ratio of the square to the other. For a genuinely constant-volume model whose volume is independent of both \(T\) and \(P\) in a neighborhood, the degenerate relation is \(c_v=c_p\). If \(V_P=0,V_T\ne0\) (SSVol, or a constant-standard-volume Margules model with a temperature-dependent excess volume), the usual local temperature-varying path at fixed volume does not exist; returning `cp` in this case is a library convention, and must not be taken as evidence for the above physical equality, nor may the denominator in the formula be replaced by a small positive number to continue the check. The standard-state partial molar volume of a single HKFT solute can be negative; the stability inequality for pure substances is not imposed on it.

  **I2 — primary: the pressure-volume terms and the activity terms must come from the same Gibbs potential.** For the Margules excess volume, the SSVol heat capacity, the Debye–Hückel temperature and pressure derivatives, and the IdealMolal entropy correction, check the following interrelated Gibbs-derivative closures:
  \[
  \bar V_k=(\partial_P\mu_k)_{T,x},\quad
  \bar s_k=-(\partial_T\mu_k)_{P,x},\quad
  \mu_k=\bar h_k-T\bar s_k,\quad
  \bar c_{p,k}=(\partial_T\bar h_k)_{P,x}=T(\partial_T\bar s_k)_{P,x},
  \]
  \[
  V=\sum_kx_k\bar V_k=\overline W/\rho,\qquad
  (\partial_Ps)_{T,x}=-V_T,\qquad c_p=(\partial_Th)_{P,x}=T(\partial_Ts)_{P,x}.
  \]
  To avoid comparing only library functions that share the same errors, use two explicit independent analytic special cases:

  For a binary Margules phase, take \(x=x_A,y=x_B=1-x\), \(P_r=101325\) Pa, and define
  \[
  C_j=H_j-TS_j+(P-P_r)(V_{Hj}-TV_{Sj}),\quad
  g^E=xy(C_0+C_1y),\quad
  Q_j=V_{Hj}-TV_{Sj}.
  \]
  The input parameters \(H_j,S_j,V_{Hj},V_{Sj}\) use J/kmol, J/kmol/K, m³/kmol, m³/kmol/K respectively. Then
  \[
  V^E=xy(Q_0+Q_1y),\quad
  \bar V_A-V_A^\circ=y^2[Q_0+Q_1(1-2x)],\quad
  \bar V_B-V_B^\circ=x^2[Q_0+2Q_1y],
  \]
  \[
  g-\sum_kx_kg_k^\circ-RT\sum_kx_k\ln x_k=g^E,\quad
  h-\sum_kx_kh_k^\circ=xy[H_0+H_1y+(P-P_r)(V_{H0}+V_{H1}y)],
  \]
  \[
  s-\sum_kx_ks_k^\circ+R\sum_kx_k\ln x_k
   =xy[S_0+S_1y+(P-P_r)(V_{S0}+V_{S1}y)].
  \]
  In the constant-standard-volume special case, \(V_T=-xy(V_{S0}+V_{S1}y)\). Volume contraction allows \(V^E<0\); only the total \(V>0\) is required, and each excess quantity must not wrongly be required to be positive.

  For the pressure-independent standard volume \(V_k^\circ(T)\) of SSVol, we must have
  \[
  c_{p,k}^\circ(T,P)-c_{p,k}^\circ(T,P_r)
  =-T(P-P_r)\frac{d^2 V_k^\circ}{dT^2}.
  \]
  Here \(P_r\) uses the actual reference pressure of the species, and must not be assumed in general to equal the 1 atm of Margules. For a volume-temperature-polynomial, the second derivative of the input polynomial can be taken directly; for a density-temperature-polynomial, first write \(V=W/\rho(T)\) and then differentiate. This equality directly tests for the dimensional error of the omitted \(R\).

  For a dilute Debye–Hückel solution, write \(L_k=\ln\gamma_k^{(m)}\); at fixed composition we must also have
  \[
  \bar V_k-V_k^\circ=RT L_{k,P},\qquad
  \bar c_{p,k}-c_{p,k}^\circ=-R(2T L_{k,T}+T^2 L_{k,TT}).
  \]
  These two relations should be checked especially for the solvent \(k=0\): its \(L_0\) contains a composition term that is independent of the Debye parameter \(A\), so the whole of \(L_0\) must not be multiplied by \(A_P/A\) or \(A_{TT}/A\) to serve as the derivative. This check covers the `WaterProps` correction through phase-level properties, without needing access to private derivative arrays.

  **I3 — secondary: the conversion between molality and mole-fraction activity bases must not change the chemical potential.** The solvent has index 0, \(M_0=W_0/1000\) is in kg/mol, and \(m^\circ=1\) mol/kg; in the region where the solvent-free lower-limit truncation is not active, for each solute \(k>0\),
  \[
  m_k=\frac{x_k}{M_0x_0},\quad
  \gamma_k^{(x)}=\frac{\gamma_k^{(m)}}{M_0m^\circ x_0},\quad
  a_k=x_k\gamma_k^{(x)}=\frac{m_k}{m^\circ}\gamma_k^{(m)},\quad
  \frac{\mu_k-\mu_k^\circ}{RT}=\ln a_k.
  \]
  For the solvent, \(\gamma_0^{(x)}=\gamma_0^{(m)}\) and \(a_0=x_0\gamma_0^{(m)}\); the solute conversion factor must not be applied to the solvent. The values of `getMolalities` are to be interpreted in mol/kg, despite an inconsistency in a local comment in its header file; the \(W_0/1000\) used by the implementation and the definitions above take precedence. An ideal molal solute with \(\gamma_k^{(m)}=1\) does not imply \(\gamma_k^{(x)}=1\); at infinite dilution in a water solvent the latter is about 55.5, which is not an anomaly.

- Sources of the invariants (first principles / literature formula / consistency with other implementations in the library / dimensions):

  I1 follows from \(dg=-s\,dT+V\,dP\) at fixed composition, the Maxwell relations, and the chain rule along \(dV=0\); stability is used only to determine the sign of \(V_P<0\). I2 comes from the derivatives of a single Gibbs potential and the first-order homogeneity of the total quantities in the amounts of substance; the analytic Margules expressions are obtained by differentiating the input-defined \(g^E\) with respect to pressure, temperature and amounts of substance, rather than by repeatedly calling the production derivative functions. I3 comes from the definition of molality and the convention for dimensionless activities; the explicit \(m^\circ\) avoids taking the logarithm of a dimensional molality.

  Compared against the public source files at the merge commit: [Margules implementation](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/src/thermo/MargulesVPSSTP.cpp), [inherited density implementation](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/src/thermo/VPStandardStateTP.cpp), [Debye–Hückel implementation](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/src/thermo/DebyeHuckel.cpp), [molality public interface and conventions](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/include/cantera/thermo/MolalityVPSSTP.h). These sources fix the model and the API semantics; they are not an independent verification of experimental truth.

- Common operating domain (range of each input, with the basis):

  All tests vary only \(T,P,x\) or the equivalent molalities, with the species and model parameters fixed to the YAML at the merge commit. The ranges below are a conservative spot-check design and do not claim to cover the full validity domain of each equation of state; the state and the finite-difference neighborhood must lie within the species' temperature range and remain liquid.

  | Case and input range | Basis and restrictions |
  |---|---|
  | I1: `test/data/pdss_hkft.yaml`, phase `vpss_soln_pdss_hkft_phase`; \(T=300\)–350 K, \(P=1\)–10 atm; \(x_{Na+}=x_{Cl-}=q\), \(q=10^{-4}\)–0.02, \(x_{H+}=x_{OH-}=10^{-7}\), water making up the remainder to 1 | The water standard state declares 273.15–600 K; this window avoids boiling, the critical point and the HKFT 228 K singularity, and dilute aqueous solutions are chosen to avoid judging experiments with a high-concentration ideal-solution approximation. The heat-capacity-difference inequality is used only after checking that the actual \(V>0,V_P<0\). |
  | I1 constant-volume control: `IdealSolidSolnPhaseExample.yaml` / `VpssSolidSolutionExample`; \(T=300\)–400 K, \(P=1\)–5 atm, positive composition with each \(x_k\geq0.01\) | The same region as the existing consistency case; \(c_v=c_p\) is adopted only after confirming that all standard volumes are constant and there is no temperature-dependent excess volume. |
  | I2 Margules: `thermo-models.yaml` / `ethanol-water-rough`; \(T=298.15\)–308.15 K, \(P=1\)–20 atm, \(x_A=0.05\)–0.95 | The temperature range comes from the newly added case; the extended pressure range specifically activates the \(P-P_r\) term. The YAML is explicitly a fictitious rough test model, not an experimental fit for ethanol–water. The parameters are fixed to the file values; for example \(V_{H0}=-12.707\), \(V_{H1}=0.476\) cm³/mol, \(V_{S0}=-0.030\), \(V_{S1}=0\) cm³/mol/K, to be multiplied by \(10^{-3}\) when converted to m³/kmol. |
  | I2 SSVol: `Li_Liquid.yaml` / `Li(L)`; \(T=460\)–500 K, \(P=1\)–10 atm, \(x_{Li(L)}=0.1\)–0.9, the other species making up the remainder | Lies in the upper part of the newly added 300–500 K test interval, avoiding extrapolation to low-temperature liquid lithium; this is a model-consistency test and does not guarantee the measured accuracy of the model parameters. Points with \(P\ne P_r\) and \(V''\ne0\) are specifically taken. |
  | I2/I3 electrolyte: `debye-huckel-all.yaml` / `debye-huckel-dilute-IAPWS`, with `debye-huckel-dilute` as a constant-\(A\) control; \(T=300\)–350 K, \(P=1\)–10 atm; \(m_{Na+}=m_{Cl-}=10^{-4}\)–\(10^{-2}\) mol/kg, \(m_{H+}=m_{OH-}=10^{-7}\), other solutes zero | Electroneutral, with ionic strength \(I\lesssim0.011\) mol/kg, consistent with the low-concentration purpose of the dilute limit and far from the \(I=3\) truncation set in the file. The salt concentration of about 9 mol/kg in the original test can serve only as a regression point, and cannot serve as a basis for the physical validity range of dilute solutions. Logarithms are taken only for species with positive concentration. |
  | I2/I3 ideal molal: `thermo-models.yaml` / `ideal-molal-aqueous`; \(T=300\)–330 K, \(P=1\)–5 atm; CO2(aq), H2S(aq), CH4(aq) each \(10^{-3}\)–0.1 mol/kg | Around the 300/330 K of the original case; keeps water the overwhelming majority and avoids the low-solvent-fraction truncation. The fictitious species thermodynamic data here are used for verifying unit and basis consistency. |

  Cases and original state basis: [consistency-cases.yaml](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/test/data/consistency-cases.yaml), [newly added Margules parameters and fictitious-species description](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/test/data/thermo-models.yaml), [HKFT standard-state definition](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/test/data/pdss_hkft.yaml), [Debye–Hückel configuration](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/test/data/debye-huckel-all.yaml).

- Numerical check method (which public API of the merged commit, what to compute, what to compare against, suggested tolerance and why; if finite differences are needed, the step-size strategy):

  **Entry points.** The later verification uses the C++ API of the merge commit: `Cantera::newThermo(file, phaseName)`; set the composition through `setMoleFractionsByName`, or through `setMolalitiesByName` after casting to `MolalityVPSSTP&`, then set the state with `setState_TP(T,P)`. Record the composition actually returned, so that normalization and electroneutrality adjustments are not mistaken for inputs. Use the full path of the file, or `addDataDirectory`, to resolve references; reset \(T,P\) explicitly for every perturbation, restore the central state after each group of computations, and interleave the order of property queries to check whether the water state cache affects the results. A VPSS model that has an actual pressure response must not be skipped merely because `isCompressible()` returns false.

  **I1.** Read `cp_mole()`, `cv_mole()`, `thermalExpansionCoeff()`, `isothermalCompressibility()`, `meanMolecularWeight()`, `density()`, `molarVolume()`. Independently obtain \(V_T,V_P\) from \(\overline W/\rho(T\pm\delta T,P)\) and \(\overline W/\rho(T,P\pm\delta P)\), first check \(\alpha,\kappa_T\), and then compare with \(c_p-c_v\); merely substituting the public \(\alpha,\kappa_T\) back into the same formula as the production `cv` is not sufficient. Then, independently, at \(T\pm\delta T\) use a one-dimensional root finder to adjust \(P_\pm\) so that `molarVolume()` equals the central \(V\), and check \(c_v\) from the central difference of `intEnergy_mole()`. The initial pressure prediction is \(P_\pm=P\mp(V_T/V_P)\delta T\); use a bracketed root finder, and where necessary reduce \(\delta T\) to keep the pressure positive and the phase liquid; do not rely on `setState_TD`, which VPSS does not necessarily support. The volume-root residual target is relative \(10^{-11}\); if the water solver cannot reach it, record the error bound rather than claiming that the derivative has passed. The true constant-volume control checks only \(c_p=c_v\) and \(c_p=\partial_Tu\).

  **I2.** Use `getChemPotentials`, `getPartialMolarVolumes`, `getPartialMolarEntropies`, `getPartialMolarEnthalpies`, `getPartialMolarCp`, `getStandardVolumes`, `getStandardChemPotentials`, `getEnthalpy_RT`, `getEntropy_R`, `getCp_R` to obtain the arrays at the central and perturbed states; the last three must be multiplied by \(RT,R,R\) respectively to restore units. Compare, species by species, \(\mu_P\) with \(\bar V\), \(-\mu_T\) with \(\bar s\), and \(\partial_T\bar h\) with \(\bar c_p\), and at the same time check the phase-level Maxwell relation and the volume additivity. The Margules \(g^E,h^E,s^E,V^E,\bar V^E\) are computed independently from the input-parameter polynomials above, and compared at 1 atm and at no fewer than two non-reference pressures; one must not check only the total \(h\) and let the small excess terms be masked by the large standard-state enthalpy. For SSVol, take \(V''\) directly from the input polynomial and compare the difference of `getCp_R()*R` at two pressures. The two additional Debye equalities construct \(L_P,L_T,L_{TT}\) from central differences of the logarithm of `getMolalityActivityCoefficients`, and report the solvent residual specifically; compared with checking only the mixture \(c_p\), this locates solvent-term errors better.

  **I3.** At the same state, read `getActivities`, `getActivityCoefficients`, `MolalityVPSSTP::getMolalityActivityCoefficients`, `getMolalities`, `getMoleFractions`, `molecularWeight(0)`, and both chemical potentials. Verify all the conversions of I3 directly; for the chemical-potential comparison use \((\mu-\mu^\circ)/(RT)-\ln a\), to avoid exponent overflow. Repeat at no fewer than two dilutions; testing only the pure solvent is not enough, since an erroneous solute conversion factor is invisible there.

  **Differences and tolerances.** Start first-order temperature derivatives from \(\delta T/T=10^{-4}\) and first-order pressure derivatives from \(\delta P/P=10^{-3}\), trying \(h/2,h,2h\) for each; in the smooth case use \([f(q+h)-f(q-h)]/(2h)\), and for second derivatives \([f(q+h)-2f(q)+f(q-h)]/h^2\), the latter starting from \(\delta T=0.1\) K. Accept a result only when the neighborhood is valid and a convergence plateau appears; estimate the truncation error from the step-size differences, and if cancellation error appears when the step is too small, fall back, rather than mechanically shrinking the step. HKFT internally uses relative \(10^{-4}\) differences and WaterProps uses an internal temperature step of 0.1 K, so machine-precision derivatives are not required.

  For the difference-free multiplicative relations of I3, a relative error of \(10^{-10}\) plus a dimensionless absolute error of \(10^{-14}\) is suggested; for the chemical-potential logarithm residual, \(10^{-9}\) is suggested. For energy algebraic closure, \(|r|\le10^{-5}\) J/kmol plus \(10^{-12}\) times the sum of the absolute values of the terms involved is suggested, to accommodate the large-number cancellation of standard states. For the analytic excess volume and the density summation, an absolute \(10^{-11}\) m³/kmol plus a relative \(10^{-9}\) is suggested. First-order derivatives default to a relative error of \(10^{-6}\), relaxed to \(10^{-5}\) for water/HKFT; the absolute floor is \(10^{-3}\) J/kmol/K for heat capacities, \(10^{-8}\) m³/kmol for \(\mu_P\) obtained from volume derivatives, \(5\times10^{-14}\) Pa⁻¹ for \(\kappa_T\), and \(10^{-9}\) K⁻¹ for \(\alpha\). The heat-capacity correction obtained from second-order activity derivatives may start from a relative \(10^{-4}\) and an absolute \(10^{-2}\) J/kmol/K tolerance, but must also pass the first-order enthalpy-derivative check; for a near-zero correction a higher but still dilute concentration should be chosen so that the signal exceeds the numerical error. The tolerances are suggestions for later trials and have not been calibrated by running; when the step size does not converge or the signal is below the error floor, mark "numerically undecided", which must not be counted as a pass.

- Do the PR's own tests already check this invariant (yes/partial/no, naming test files and assertions):

  **I1: partial.** `test/thermo_consistency/consistency.cpp` adds `cv_eq_cp_minus_T_V_beta2_over_kappa_T`, and has finite-difference checks `betaT_eq_minus_dmv_dP_const_T_div_mv` and `alphaV_eq_dmv_dT_const_P_div_mv`; it adds the instances `MargulesSSVol` and `MargulesWithExcessVolume`, and HKFT sets `rtol_fd=4e-6`, `atol_c=5e-14`. But `cv_eq_dudT` and `cv_eq_dsdT_const_v_times_T` are skipped outright for `!isCompressible()`, the newly added heat-capacity-difference assertion shares its formula with the production code, and the zero-\(\kappa_T\) branch does not require \(V_T=0\). Hence it is not a complete verification of an independent constant-volume path within the domain restricted here.

  **I2: partial.** The same file has `gk_eq_hk_minus_T_sk`, `v_eq_sum_vk_Xk`, `cp_eq_dhdT`, `cp_eq_dsdT_const_p_times_T`, `cpk0_eq_dhk0dT` and `dsdP_const_T_eq_minus_dV_dT_const_P`. `consistency-cases.yaml` removes the relevant known-failure marks for IdealMolal, Debye with variable \(A\), and HMW; `test_problems/cathermo/ims/output_blessed.txt` updates the entropy regression values, and `test_problems/cathermo/HMW_graph_CpvT/output_blessed.txt` updates the heat-capacity results. The three central states of the newly added `ethanol-water-rough` are all at 1 atm (the Maxwell check still perturbs the pressure slightly), the analytic excess quantities over a wide pressure range are not explicitly covered, and the partial molar volumes are not verified species by species with \(\mu_P\); updating reference outputs does not by itself constitute independent physical proof.

  **I3: yes (core identities), partial (independent basis conversion).** `activity_coeffs` asserts \(a_k=\gamma_k^{(x)}x_k\), with tolerance \(10^{-10}\) relative plus \(10^{-20}\); `chem_potentials_to_activities` asserts \(\exp[(\mu_k-\mu_k^\circ)/(RT)]=a_k\), with tolerance \(10^{-9}\) relative plus \(10^{-14}\). The revised molal phase no longer marks this as a known failure, but these assertions do not independently reconstruct the full conversion through \(W_0/1000\) and \(\gamma^{(m)}\). The above is a judgment of source-code coverage, not the result of this execution. Test basis: [merged consistency.cpp](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/test/thermo_consistency/consistency.cpp).

- Do you expect the change to satisfy it (satisfies / may not satisfy / uncertain), with one or two sentences of reasoning (speculation, not a conclusion):

  **I1: satisfies (speculation, only for the stable compressible domain above and the true constant-volume domain)**; the new formula correctly uses the mixture volume derivatives; for models with zero compressibility but nonzero thermal expansion, the physical interpretation of `cv=cp` is **uncertain** and cannot be proved from this implementation. **I2: satisfies (speculation)**; adding \((P-P_r)V^E\), supplying the missing \(R\), and splitting out the solvent-specific derivatives all conform to the Gibbs derivatives and to dimensions; the remaining risks are the numerical error of the internal water properties and the insufficiently covered non-reference pressures. **I3: satisfies (speculation)**; dividing by \(x_0W_0/1000\) is exactly the factor required for the activity-basis conversion.

- Environment: merge commit SHA, language and dependencies needed, estimated setup difficulty (low/medium/high):

  Merge SHA: `d58a23265414e1789536cf5332d9fb29de095a43`. Preferably the C++20 Cantera library together with the same-commit YAML above (including the referenced species files), using the C++ public interface, with no need for the Python bindings; if building from this SHA later, a Python/SCons build driver and Cantera's dependencies such as Eigen, fmt, yaml-cpp and SUNDIALS are generally needed, and re-running the original GTest suite additionally needs GoogleTest/Boost headers, with exact versions resolved by the build checks of that commit. [SConstruct](https://github.com/Cantera/cantera/blob/d58a23265414e1789536cf5332d9fb29de095a43/SConstruct) specifies C++20 and provides the dependency configuration; estimated **medium**, the main cost being obtaining the development-version library matching the SHA and the complete test data, rather than the numerical workload.

  Per the merged source: the PR text claims to implement `cv_mole` for HKFT, but `PDSS_HKFT.h/.cpp` at this SHA has no such override; do not treat a direct call to `PDSS_HKFT::cv_mole()` as an implemented check entry point; this document uses `ThermoPhase::cv_mole()` of a phase containing HKFT species. For this document, only `input/pr.json`, `input/pr.diff` and the public files at the specified SHA on GitHub were read, and this document was written; no repository was cloned, nothing was installed, compiled or numerically executed, and no other audit track or file outside the working directory was read.

- Formalization (Lean) feasibility: can the invariant be written as a proposition over finite-dimensional reals, and what are the main obstacles:

  **Yes.** Represent the state by \((T,P,n_0,\ldots,n_{N-1})\in\mathbb R^{N+2}\), \(n_k>0\), \(x_k=n_k/\sum_jn_j\), and define \(G=(\sum n_k)g\); then the partial molar quantities and the Maxwell relations of I2 are propositions on partial derivatives of finite-dimensional real functions. The Margules special case is essentially differentiation of polynomial/rational functions plus first-order homogeneity, the SSVol pressure heat-capacity term is a second derivative with unit conversion, and I3 is an algebraic and logarithmic equality over positive reals; these are suitable for formalization first.

  I1 needs \(V_P\ne0\) to be made explicit and the implicit-function theorem/chain rule to be used to construct the constant-volume path, after which the non-negativity of the heat-capacity difference is proved under \(V_P<0\); the branch \(V_P=0,V_T\ne0\) should be treated as a model degeneracy that violates the premises, and the division identity must not be carried over. The main obstacles are the mathematical correspondence between the full HKFT/IAPWS parameterization and the iterative density solution, the piecewise smoothness caused by truncation, and the error bounds connecting the real-number model proof to C++ floating-point execution, cache state restoration and internal finite differences; a proof of the abstract Gibbs identities does not automatically prove that the production implementation satisfies them.
