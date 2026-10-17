## e2nIEE/pandapower#3016 Neutral earthing impedance (rn_ohm) and zigzag earthing transformers in the zero-sequence short circuit

- What changed (1–3 sentences, naming files and formulas):

  `pandapower/create/trafo_create.py` adds `rn_ohm` to the four single/bulk transformer creation APIs; `pandapower/pd2ppc_zero.py` writes the neutral impedance as (Z_N=R_N+jX_N), adds (3Z_N) to an ordinary earthed winding, and places it after the IEC transformer correction. A new zigzag branch represents the zero-sequence decoupling of the two sides by an HV-to-ground branch (α Z_{0k}+3Z_N) and (only for `ZNyn`) an LV-to-ground branch (Z_{0m}+(1-\alpha)Z_{0k}); the tests are in `pandapower/test/shortcircuit/test_1ph.py`, and the description is in `CHANGELOG.rst`. Note that the `ZN` in the input PR description is inconsistent with the merged source: the actual branch accepts the case-insensitive `znyn/znd/zny/zd`, and has no `zn`; the physical checks below use the explicit `ZNyn/ZNd/ZNy`, so the name-support issue is not mixed into the physical conclusions.

- Invariants (mathematical statements; write them as checkable equalities/inequalities, give variables, sign conventions, domain of validity):

  Unified conventions: linear, power-frequency, symmetrical-component short-circuit equivalent circuit, RMS phasors, (j^2=-1), passive impedance (Z=R+jX), inductive (X\ge0), port currents positive into the element. (U_b) is the rated line-to-line voltage of the bus (kV), (S_B) is an arbitrary calculation base (MVA), and (Z_{B,b}=U_b^2/S_B) (Ω); `rn_ohm/xn_ohm` are the actual Ω on the neutral side of the earthed winding, not per-unit values. (Z_{q,b}) denotes the q-sequence driving-point Thévenin impedance of bus b, and (I_f) is the magnitude of the single-phase-to-ground initial symmetrical short-circuit current (kA). All comparisons keep the topology, the positive-sequence equipment parameters and the voltage factor unchanged, with zero fault impedance, and no current-source inverters, generators or capacitive shunt branches to ground.

  **I1 — primary: the three-fold neutral impedance must sit at the correct winding port, and must not change the positive/negative sequence.**

  From (I_N=I_a+I_b+I_c=3I_0), the voltage drop at the neutral is (3Z_NI_0), so the series impedance added to the zero-sequence port of that winding is exactly (3Z_N), and the increments of the positive and negative sequences are zero. For an isolated fixture in which the fault point is at the bus of the earthed winding and that winding is the only zero-sequence return entry of that bus:

  \[
  Z_{0,b}(Z_N+\Delta Z_N)-Z_{0,b}(Z_N)=3\Delta Z_N,
  \quad Z_{1,b}(Z_N+\Delta Z_N)=Z_{1,b}(Z_N),
  \quad Z_{2,b}(Z_N+\Delta Z_N)=Z_{2,b}(Z_N).
  \]

  This means that a resistance perturbation changes only the zero-sequence real part, and a reactance perturbation changes only the imaginary part; the coefficient is not (3K_T), nor (3\alpha). Test first with the LV earthed end of `Dyn` and the HV earthed end of `YNd`; do not extend this driving-point equality to an arbitrary bus with several parallel earthing paths. Conversion across sides must obey (Z'=a^2Z), where (a=U_H/U_L); the per-unit neutral increment is (3\Delta Z_N/Z_{B,\mathrm{earth}}).

  The supplementary port check for `YNyn` also belongs to I1: following the HV-side selection semantics of the merged source, let `rn_ohm/xn_ohm` denote the HV neutral impedance, with the LV neutral solidly earthed. With all impedances referred to HV, let the original total zero-sequence leakage impedance be (Z_{k,r}), the corrected one (Z_k=K_TZ_{k,r}), the magnetizing impedance (Z_m), and (α=\texttt{si0_hv_partial}). With an LV external grid supplying and the HV being the fault end with no other earthing path, we should have

  \[
  Z_{0,H}=\alpha Z_k+3Z_N+
    \frac{Z_m\big[(1-\alpha)Z_k+Z_{g0,H}\big]}
         {Z_m+(1-\alpha)Z_k+Z_{g0,H}}.
  \]

  Hence for finite (Z_m,Z_{g0,H}) it must still hold that (\Delta Z_{0,H}=3\Delta Z_N). This formula expresses that the neutral impedance is added outside the HV branch, and must not first be added into (Z_k) and then split between HV/LV; this is the physical requirement after declaring the HV neutral semantics, and does not claim that this API can independently specify two neutrals.

  **I2 — primary: zero-sequence port decoupling of the zigzag and its respective earthing impedances.**

  For the ideal zigzag equivalent model declared by this PR, with a common power base and the respective winding voltage bases, let (z_k,z_m) be the corresponding per-unit impedances and (α\in(0,1)). We should have

  \[
  \begin{bmatrix}i_{0,H}\\i_{0,L}\end{bmatrix}=
  \begin{bmatrix}1/z_H&0\\0&y_L\end{bmatrix}
  \begin{bmatrix}v_{0,H}\\v_{0,L}\end{bmatrix},\qquad
  z_H=\alpha z_k+3Z_N/Z_{B,H},
  \]
  \[
  y_L=\begin{cases}
    1/[z_m+(1-\alpha)z_k],&\texttt{ZNyn},\\
    0,&\texttt{ZNd},\texttt{ZNy}.
  \end{cases}
  \]

  The zero-sequence cross admittance (Y_{0,HL}=Y_{0,LH}=0) is the core structural constraint. In a three-bus fixture in which the main transformer LV is delta, the system is earthed only by the zigzag, and the aux bus has no other equipment, (Z_{0,H}=Z_{B,H}z_H), and for `ZNyn` (Z_{0,L}=Z_{B,L}[z_m+(1-\alpha)z_k]); for `ZNd/ZNy` the ideal earth-fault current of the aux bus is zero. The HV (Z_0) does not change with the secondary `yn/d/y` type or the magnetizing impedance; the LV (Z_0) of `ZNyn` does not change with the HV neutral impedance. If an external grid is used to change the LV zero-sequence earthing strength while keeping the positive-sequence strength unchanged, the HV (Z_0) should still be unchanged. The total fault current after an arbitrary change of connection is not required to be the same here, because its positive-sequence network may also change.

  Here "the cross admittance is zero" is a structural requirement of the ideal zigzag ampere-turn cancellation model; how much of the input `vk0_percent` is interpreted as the zigzag's own leakage impedance depends on the parameterization convention of that PR. It cannot be claimed merely because the code uses (α z_k) that the actual equipment also conforms to this parameter mapping; its literature/measurement verification boundary must be kept separately.

  **I3 — secondary: three-sequence series fault formula and NER current limiting.**

  For the above passive fixture without rotating machines or current sources, (Z_2=Z_1), and a single-phase metallic earth fault satisfies

  \[
  I_0=I_1=I_2=I_f^{\mathrm{phasor}}/3,\qquad
  I_f=\frac{\sqrt3\,cU_b}{|2Z_{1,b}+Z_{0,b}|}.
  \]

  In the unique-return-path fixture of I1 or the zigzag HV fixture of I2, fix (X_N\ge0), define (D=2Z_{1,b}+Z_{0,b}(R_N=0)=A+jB), and assume (A,B\ge0); then

  \[
  I_f(R_N)=\frac{\sqrt3\,cU_b}{\sqrt{(A+3R_N)^2+B^2}},\quad
  \frac{dI_f}{dR_N}=-\frac{3\sqrt3\,cU_b(A+3R_N)}{[(A+3R_N)^2+B^2]^{3/2}}\le0.
  \]

  For (R_N>0), (I_f\le cU_b/(\sqrt3R_N)), and (R_NI_f\to cU_b/\sqrt3). This monotonicity is not extended to arbitrary networks containing resonant capacitances, parallel earthing loops or current sources.

- Sources of the invariants (first principles / literature formula / consistency with other implementations in the library / dimensions):

  I1 comes from Kirchhoff's current law, the definition of symmetrical components and Ohm's law at the neutral; the correct branch position comes from the port circuit above. The decoupling of I2 comes from the cancellation of zero-sequence ampere-turns in an ideal zigzag winding, and its numerical split of the leakage impedance is the equivalent model declared by this PR. I3 comes from the three-sequence series connection under the fault conditions (I_b=I_c=0,V_a=0) and (E_1=cU_b/\sqrt3); the current-limiting inequality and the derivative follow directly from this formula. The unit relations (\mathrm{kV}^2/\mathrm{MVA}=\Omega) and (\mathrm{kV}/\Omega=\mathrm{kA}) serve as independent dimensional checks.

  Merged-source evidence (all pinned to this merge SHA): [zero-sequence modeling and magnetizing/leakage impedance definitions](https://github.com/e2nIEE/pandapower/blob/0c6a1f0a65273af91dd1d5b854498cd89d457939/pandapower/pd2ppc_zero.py), [IEC correction factor implementation](https://github.com/e2nIEE/pandapower/blob/0c6a1f0a65273af91dd1d5b854498cd89d457939/pandapower/build_branch.py), [single-phase short-circuit formula](https://github.com/e2nIEE/pandapower/blob/0c6a1f0a65273af91dd1d5b854498cd89d457939/pandapower/shortcircuit/currents.py), [public result column extraction](https://github.com/e2nIEE/pandapower/blob/0c6a1f0a65273af91dd1d5b854498cd89d457939/pandapower/shortcircuit/results.py). These are evidence of API/parameter semantics and in-library consistency, and do not serve as independent proof of physical correctness.

  The text and source comments of the input [PR #3016](https://github.com/e2nIEE/pandapower/pull/3016) mention IEC 60909, Fig. 8.3 of the PowerFactory transformer technical reference, and Appendix 5 of the J&P Transformer Book; this time the corresponding versions and figures of these references were not independently obtained, so the exact zigzag parameter mapping, the conformity with IEC clauses, and the PowerFactory real-network verification claimed by the PR are all marked unverified. Comments, issues, forks or another audit track were not read.

- Common operating domain (range of each input, with the basis):

  A finite, linear, saturation-free short-circuit equivalent model is used; arcs, transient peaks, thermal limits, three-phase unbalanced power flow, and inference of the actual core structure are not covered. The ranges below are engineering scales and input-validity limits for stable random spot checks, and do not claim that IEC prescribes these upper and lower limits for every parameter.

  | Input | Spot-check range / fixed value | Basis |
  |---|---|---|
  | Voltage | (U_H=10\ldots110) kV, (U_L=0.4\ldots33) kV, (U_H>U_L); 110 kV usable at the main-transformer source end | Covers the PR's 20/0.4 and 110/20 kV fixtures |
  | Equipment rated power | (S_T=0.4\ldots40) MVA | From the PR's small earthing transformer to the main-transformer scale |
  | Arbitrary network power base | (S_B\in\{1,10,100\}) MVA | Changing the calculation base should not change the Ω/kA physical results |
  | Positive/zero-sequence short-circuit parameters | `vk_percent`, `vk0_percent` each 4…15; `vkr/vk`, `vkr0/vk0` each 0.02…0.4 | Nonzero, passive, and square roots strictly positive; avoids zero-parameter fallback logic |
  | Zero-sequence split / magnetizing | (α=0.1\ldots0.9), `mag0_percent`=10…1000, `mag0_rx`=0…1 | Internal split non-degenerate; values covered by the tests are 0.9, 100, 0 |
  | Neutral impedance | (R_N,X_N\ge0); baseline 0, perturbation such that (3R_N/|Z_{k,\mathrm{earth}}|) and (3X_N/|Z_{k,\mathrm{earth}}|\in[0.1,100]), including the PR's 0.01 Ω, 0.02 Ω, 5 Ω | Covers slight changes up to NER dominance, avoiding extremely ill-conditioned scales |
  | External grid | `s_sc_max_mva`=100…10000 MVA, `rx_max`=0.05…0.5, `x0x_max`=0.5…3, `r0x0_max`=0.05…0.5; give the corresponding positive values when min is needed | Finite passive source, covering the PR's 1000/5000 MVA |
  | Fixed modeling conditions | Rated winding voltage equal to the rated value of the connected bus; no tap, `parallel=1`, `in_service=True`, not a power-station unit; `pfe_kw=i0_percent=0`; 50 Hz | Excludes extra factors from the transformation ratio, parallel-neutral semantics and the power-station correction |
  | Solver options | `fault="1ph"`, `case="max"`, `lv_tol_percent=10`, `use_pre_fault_voltage=False`, `r_fault_ohm=x_fault_ohm=0`, `branch_results=False`; both `inverse_y` | Standard short-circuit path of the public API read; in this combination (c=1.1) at every voltage level |

  `shift_degree` is specified separately, and the vector group carries no clock-hour digits; the 330° of the newly added test can be reused, and the `YNyn` fixture uses 0°. This specification does not presume that `parallel>1` means a shared or independent NER, and does not cover the composite correction of power-station units. The ideal open-circuit check is handled separately, and its result is not required to be a finite complex impedance.

- Numerical check method (which public API of the merged commit, what to compute, what to compare against, suggested tolerance and why; if finite differences are needed, the step-size strategy):

  **Entry points and observables.** Build the fixtures below with `pandapower.create_empty_network`, `create_bus`, `create_ext_grid`, `create_transformer_from_parameters`, call `pandapower.shortcircuit.calc_sc`, and read only the public `rk0_ohm/xk0_ohm/rk_ohm/xk_ohm/ikss_ka` of `net.res_bus_sc`. Let (Z_0=\texttt{rk0_ohm}+j\texttt{xk0_ohm}) and (Z_1=\texttt{rk_ohm}+j\texttt{xk_ohm}); do not call private conversion functions, and do not use `_ppc0` as an acceptance entry. Save each result before running the next kind of fault, to prevent the result table from being overwritten. [API signature and default model](https://github.com/e2nIEE/pandapower/blob/0c6a1f0a65273af91dd1d5b854498cd89d457939/pandapower/shortcircuit/calc_sc.py).

  **I1.** In the two-bus `Dyn` fixture, put the ext_grid at HV and observe LV; in the reverse fixture `YNd`, put the ext_grid at LV and observe HV. Run (Z_N=0,R,jX,R+jX) respectively, and directly check the complex-increment equality, the unchanged positive sequence, and the unchanged currents of `fault="3ph"` and `fault="2ph"` (an indirect public observation of the negative-sequence influence). At the same time change the positive-sequence `vk/vkr` so that (K_T) changes, then compare the neutral increment within each fixed parameter group, verifying that the NER is not scaled by (K_T). The total (Z_0) must not be compared directly across different (K_T).

  `YNyn` supplementary fixture: two buses, only the LV has an ext_grid and the HV has no other element; compare (Z_0(R+jX)-Z_0(0)) of HV with (3(R+jX)), and use the T-circuit formula given in I1 as an independent absolute-value oracle. Let (p_0=\texttt{vk0_percent},r_0=\texttt{vkr0_percent},p=\texttt{vk_percent},r=\texttt{vkr_percent},q=\texttt{mag0_rx},m=\texttt{mag0_percent}/100); then in HV Ω

  \[
  Z_{k,r}=\frac{U_H^2}{S_T}\frac{r_0+j\sqrt{p_0^2-r_0^2}}{100},\quad
  K_T=\frac{0.95(1.1)}{1+0.6\sqrt{p^2-r^2}/100},\quad
  Z_m=m|Z_{k,r}|\frac{q+j}{\sqrt{1+q^2}}.
  \]

  To simplify the source-impedance oracle, take `x0x_max=1`, `r0x0_max=rx_max=ρ`, so that the LV grid referred to HV has (Z_{g0,H}=1.1U_H^2(\rho+j)/(S_g\sqrt{1+\rho^2})). Sweep finite magnetizing impedances and (α), to avoid using the infinite-magnetizing approximation to mask a neutral placed on the wrong branch. Here (Z_m) follows the library's existing definition without multiplication by (K_T); the audit target is the position of the neutral, not a separate ruling on the magnetizing correction rule.

  **I2.** Use the three-bus structure of the newly added test: 110 kV ext_grid → 110/20 kV `YNd` main transformer → 20 kV bus → 20/0.4 kV zigzag under test → aux bus; the aux has no load/source. Give all transformers complete positive/zero-sequence parameters; test `ZNyn/ZNd/ZNy`. Use the above (Z_{k,r},K_T,Z_m) formulas but with the earthing-transformer parameters to compute (z_k=K_TZ_{k,r}/Z_{B,H}) and (z_m=Z_m/Z_{B,H}), and compare with the HV/LV Ω oracles of I2; the main transformer delta blocks source-side zero sequence, so the HV oracle does not contain the zero-sequence impedance of the main grid. Note that LV Ω results must be multiplied by (Z_{B,L}), and the HV Ω must not be used directly.

  With the leakage impedance held fixed, perturb `rn_ohm/xn_ohm`, `mag0_percent/mag0_rx` and the secondary connection respectively, to check the independence of I2. Then add an ext_grid at the aux, keep its positive-sequence `s_sc_max_mva/rx_max` fixed and change only `x0x_max/r0x0_max`, and check that the HV (Z_0) is unchanged, testing the decoupling by actually changing the LV zero-sequence termination. For the passive aux open circuit of `ZNd/ZNy`, check only (I_f\approx0), and allow `rk0_ohm/xk0_ohm` to be infinity; the merged implementation simulates an open circuit with a large impedance, and result extraction also converts a very large impedance to infinity, so no `inf-inf` may be taken.

  **I3.** For each finite result, use the (Z_0) from the independent oracle of I1/I2 and the public result (Z_1) to compute (\sqrt3(1.1)U_b/|2Z_1+Z_0|), and compare with `ikss_ka`; another computation using the public (Z_0) serves only as an internal consistency check of the result, and cannot replace the independent oracle. For fixed (X_N) do an increasing (R_N) sweep to check monotonicity and the upper bound; the result of `fault="3ph"` in the same network can be used to cancel the voltage factor via (I_{1ph}/I_{3ph}=3|Z_1|/|2Z_1+Z_0|).

  The large-NER limit cannot be judged with an arbitrary fixed "sufficiently large" resistance: let (C=\sqrt3cU_b) and (t=|D|/(3R_N)<1); then we should have
  \(1/(1+t)\le I_f/(C/(3R_N))\le1\). Use (t\le0.01) to give a theoretical deviation of no more than about 1%, and ensure that the result has not entered the open-circuit truncation range.

  Use `inverse_y=True/False` for each fixture, and change (S_B=1,10,100) to check that the same physical Ω/kA results are obtained; these are auxiliary numerical consistency checks, not independent physical oracles. The core equalities use the complex-modulus residual ( |\hat Z-Z_*|\le10^{-8}\,\Omega+10^{-7}|Z_*|\), and currents use ( |\hat I-I_*|\le10^{-8}\,\mathrm{kA}+10^{-7}|I_*|\); the near-zero open-circuit current threshold can use the (10^{-6}\) kA of the newly added test. The basis is that a small-scale, non-degenerate double-precision linear circuit should be far better than an engineering percentage error; for ill-conditioned/open-circuit results, first record the solver failure and the scale separately, and do not blindly relax the threshold. Choose increments up to (3|\Delta Z_N|\ge0.1|Z_{k,\mathrm{earth}}|), to avoid the signal being drowned by subtraction of large numbers. Finite differences are not needed: the increment relations and the current-limiting formula can both be checked directly; an optional derivative check uses central differences (h\sim\epsilon^{1/3}\max(R_N,|D|/3)), requires (R_N>h), and compares with the analytic derivative after comparing the stability of (h,h/2,2h).

- Do the PR's own tests already check this invariant (yes/partial/no, naming test files and assertions):

  **I1: partial.** `test_trafo_neutral_earthing_impedance` in [test_1ph.py](https://github.com/e2nIEE/pandapower/blob/0c6a1f0a65273af91dd1d5b854498cd89d457939/pandapower/test/shortcircuit/test_1ph.py) asserts, in a fixture with `Dyn`, `parallel=1` and a 1 MVA network base, (\Delta R_0=0.03\) Ω and (\Delta X_0=0.06\) Ω, with the other component unchanged; in the resistance case the positive-sequence R/X are unchanged, with `atol=1e-6, rtol=0`. Both solution paths are parameterized; it does not cover the HV earthed group, the branch position of `YNyn`, mixed R+jX, different correction factors, or the calculation base.

  **I2: partial.** `test_zigzag_earthing_transformer` in the same file checks that the `ZNyn` HV current is positive, that its ratio to `YNd` is in 0.5…2, that the `ZNyn/ZNd` HV currents are equal, that the `ZNyn` LV current is positive and that the `ZNd` LV current is about zero. It does not verify the exact impedance of (α z_k) and of the LV magnetizing branch, and does not test `ZNy`, magnetizing perturbation, LV termination perturbation, or the independence of the LV from the HV NER; "the same order of magnitude" cannot confirm that the leakage-impedance formula is correct.

  **I3: partial.** The two newly added tests only assert that the fault current decreases after adding the specified resistance; they do not directly check the absolute-value three-sequence formula, continuous monotonicity, the analytic upper bound or the NER-dominated limit. The above is obtained from reading the assertions, and does not mean that the tests were executed this time, nor is the PR text's "all existing tests pass" taken as current run evidence.

- Do you expect the change to satisfy it (satisfies / may not satisfy / uncertain), with one or two sentences of reasoning (speculation, not a conclusion):

  **I1: may not satisfy.** Speculation: the standard short-circuit path of `Dyn/YNd` can satisfy the three-fold increment, but `YNyn`, which adds (3Z_N) into the total leakage impedance and then splits it by (α,1-\alpha), may violate the port formula that the HV neutral should be fully in series with the HV branch when the magnetizing branch is finite; not yet numerically verified.

  **I2: uncertain.** Speculation: a single, non-power-station fixture with rated ratio can satisfy the declared earthing branches and decoupling structure, because the implementation retains a very large HV–LV impedance and adds the bus admittances separately; whether the actual zigzag leakage impedance can be uniquely and correctly characterized by these inputs still lacks confirmation from independent literature/equipment data.

  **I3: satisfies (limited to fixtures that implement I1/I2 correctly; speculation).** The merged source explicitly computes the current with (2Z_1+Z_0), and a unique passive earthing path should have the current-limiting behavior above; if the upstream zero-sequence network is wrong, the internal consistency of the current formula cannot prove the physical result correct.

- Environment: merge commit SHA, language and dependencies needed, estimated setup difficulty (low/medium/high):

  Merge SHA: `0c6a1f0a65273af91dd1d5b854498cd89d457939`. Python; the [pyproject.toml](https://github.com/e2nIEE/pandapower/blob/0c6a1f0a65273af91dd1d5b854498cd89d457939/pyproject.toml) of this commit requires Python >=3.10, with runtime dependencies `numpy>=1.26,<2.4`, `scipy<1.17`, `pandas~=2.3`, `networkx~=3.4`, `packaging~=25.0`, `tqdm~=4.67`, `deepdiff~=8.6`, `geojson~=3.2`, `typing_extensions~=4.9`, `pandera~=0.26.1`. Repeating the original tests additionally needs pytest and the pytest-timeout used by its configuration; PowerFactory, a GPU or actual equipment is not needed. Estimated **low** (a common CPU environment, numerical-library wheels); specific platform compatibility is unverified. This time only the inputs and the pinned-commit source on GitHub were read and the specification written; no project code was installed, compiled or run.

- Formalization (Lean) feasibility: can the invariant be written as a proposition over finite-dimensional reals, and what are the main obstacles:

  Yes. Represent phasors by the equivalent two-dimensional real vectors of complex numbers, and nodal admittances by finite-dimensional real matrices; I1 is an algebraic identity of (I_N=3I_0) and series impedances, the supplementary `YNyn` formula is finite-circuit elimination; I2 is an identity of block-diagonal matrices and two-port admittances; the squared-current formula, the monotonicity inequality and the finite error bounds of I3 can be reduced to rational expressions, square roots and positive-denominator assumptions.

  Passivity, nonzero denominators, invertibility/earthing conditions, ideal-transformer voltage conversion and the zigzag ampere-turn assumption must be declared explicitly. The main obstacles are proving the mapping from real windings to the declared equivalent circuit, making the meaning of the IEC correction and the input parameters explicit, and relating this mathematical model to the floating-point and sparse-solver implementations of pandas/NumPy/SciPy; code comments are not sufficient to prove these premises. An open circuit should be represented by zero admittance or a limit, and IEEE infinity should not be treated as a real number; the absolute mathematical propositions and the error proofs for software tolerances/large-impedance approximations need to be separated.
