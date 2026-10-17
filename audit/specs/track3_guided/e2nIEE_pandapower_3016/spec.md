## 1. What changed

This specification reviews pandapower PR #3016, **“Neutral earthing impedance (rn_ohm) and zigzag earthing transformers in the zero-sequence short circuit”**, at merge commit `0c6a1f0a65273af91dd1d5b854498cd89d457939`, relative to base `be6b67a8562a448d538f385a43c46718956cbc31`. It is a source-based specification and proposed numerical check, not a report of executed simulations.

- `pandapower/create/trafo_create.py`: `create_transformer`, `create_transformers`, `create_transformer_from_parameters`, and `create_transformers_from_parameters` gain `rn_ohm`, defaulting to NaN. Single-element creators store non-NaN values; bulk creators propagate/store scalar or iterable values. The `xn_ohm` description now explicitly says reactance. The new argument is inserted before existing optional arguments, so positional callers can be affected; keyword-based calls are assumed below.
- `pandapower/pd2ppc_zero.py`, `_add_trafo_sc_impedance_zero`: combines the neutral parameters as (Z_N=R_N+jX_N\), treating missing/NaN resistance and NaN reactance as zero. The existing power-station-unit (PSU) expression changes from a purely reactive neutral to
  \[
  z_{0,\mathrm{PSU}}=\frac{k_{ST}Z_{0T,H}/p+3Z_N}{U_{bH}^2/S_b},
  \]
  where (p\) is `parallel`, (S_b\) is network MVA, and (U_b\) is line-to-line bus kV. (Z_{0T,H}\) is the uncorrected transformer zero-sequence impedance in ohms referred to its tap-adjusted HV voltage.
- For ordinary `Dyn`, `YNd`, `Yyn`, `YNyn`, `YNy`, and `Yzn`, the same function adds
  \[
  z_N=3Z_NS_b/U_{be}^2
  \]
  to the corrected leakage impedance **before** constructing the subsequent T-to-pi equivalent circuit. The selected earthing side (e\) is HV for `YNd`, `YNy`, and `YNyn`, and LV for the other three. PSU rows skip this ordinary addition. This does not by itself prove that the neutral ends up in the correct terminal arm of a two-port model.
- New case-insensitive branches accept **`ZNyn`, `ZNd`, `ZNy`, and `Zd`**, not `ZN`. This differs from the PR body and changelog, which advertise `ZN`. For these branches the ordinary neutral addition is deferred. The zero-sequence HV–LV series branch retains the pre-existing open-circuit approximation (R=X=10^{20}S_b\). Direct bus shunts are added with
  \[
  z_H=s z_\ell+z_N,\qquad
  z_L=z_m+(1-s)z_\ell\quad\text{(`ZNyn` only)},
  \]
  \[
  G_H+jB_H=S_b\,\mathbf1_{\rm service}/z_H,\qquad
  G_L+jB_L=S_b\,\mathbf1_{\rm service}/z_L.
  \]
  Here (s=\texttt{si0_hv_partial}\), (z_\ell\) is the previously computed `z0_k`, and (z_m\) is the existing zero-sequence magnetizing impedance. Other new groups have no LV ground shunt. Shunts accumulate with `np.add.at`, including multiple transformers at one bus.

For ordinary short-circuit transformers, the surrounding implementation computes
\[
z_{\ell}=\frac{k_T S_b}{p S_T}\left(\frac{U_{r,q}}{U_{b,q}}\right)^2
\left[\frac{v_{kr0}}{100}+j\sqrt{(v_{k0}/100)^2-(v_{kr0}/100)^2}\right],
\]
with (q=H\) only for `YNd`/`YNy`, otherwise (q=L\), **including the new HV zigzag groups**. In the ordinary positive-reactance domain,
\[
k_T=\frac{0.95c_{\max,L}}{1+0.6\sqrt{(v_k/100)^2-(v_{kr}/100)^2}}
\]
for `case="max"`, and (k_T=1\) for `case="min"`. Magnetizing impedance uses the uncorrected leakage magnitude: if (m=\texttt{mag0_percent}/100\), (\rho=\texttt{mag0_rx}\), and (z_{sc}\) is that magnitude in pu before parallel division, then (z_m=(z_{sc}m/p)(\rho+j)/\sqrt{1+\rho^2}\).

The neutral addition is outside the `mode == "sc"` guard: it also affects `runpp_3ph` for its supported groups (`Dyn`, `YNyn`, `Yzn`). That path multiplies the leakage base conversion by three but does not similarly change the newly added neutral conversion. New HV zigzag groups remain rejected by the earlier three-phase power-flow guard. Balanced positive/negative-sequence construction is not changed.

`pandapower/test/shortcircuit/test_1ph.py` adds two tests, each parametrized over `inverse_y=True/False`. `CHANGELOG.rst` describes the features and also contains an unrelated three-phase-generator fix note; the supplied diff contains no corresponding generator implementation change.

## 2. Invariants

All phasors below are RMS fundamental-frequency values. The default domain is a linear, balanced-parameter passive transformer and a metallic single-line-to-ground fault, with no converter current sources, saturation, capacitance, or arc dynamics. Sequence impedances in equations for current are referred to the fault bus in ohms.

**I1 — Primary: neutral-current law and correct placement.** For a grounded winding,
\[
I_N=I_a+I_b+I_c=3I_0,\quad V_N=Z_N I_N,
\quad\Delta V_0=3Z_N I_0.
\]
Thus the neutral contributes exactly (3Z_N\) to **that winding's zero-sequence arm**, once, after any correction of the transformer itself. On an isolated one-port ground-return path (ordinary `Dyn` or `YNd` in the default SC pi model, or the HV zigzag shunt),
\[
Z_{0,\rm th}(Z_N)-Z_{0,\rm th}(0)=3Z_N.
\]
This last terminal identity requires the transformer to be the sole ground-return path; it is not true for arbitrary meshed networks with other grounds.

For a finite-magnetizing `YNyn` T equivalent, with HV and LV referred to a common voltage base and HV-only neutral impedance, the correct arms are
\[
A=sZ_\ell+3Z_{NH},\quad B=(1-s)Z_\ell,\quad M=Z_m.
\]
For LV-only grounding impedance, add (3Z_{NL}\) to (B\) instead. If (a=1/A,b=1/B,m_0=1/M\), eliminating the internal node gives
\[
Y_0=\begin{pmatrix}a&0\\0&b\end{pmatrix}
-\frac{1}{a+b+m_0}\begin{pmatrix}a^2&ab\\ab&b^2\end{pmatrix}.
\]
This fixes placement more precisely than a claim about total leakage. The API has only one neutral impedance pair; it cannot specify two independent neutrals. Its docstrings do not define the side for `YNyn`; the changed code chooses HV. Use that explicit convention for checking, and report the documentation gap rather than assuming both neutrals have the same impedance. For groups without a grounded terminal (`Yy`, `Yd`, `Dy`, `Dd`), neutral parameters must have no effect.

**I2 — Primary: zigzag decoupling and the stated leakage-partition model.** For an ideal, symmetric, grounded HV zigzag with equal zig/zag sections and negligible interwinding capacitance,
\[
Y_0=\operatorname{diag}\left((sZ_\ell+3Z_N)^{-1},\,
[Z_m+(1-s)Z_\ell]^{-1}\right)
\]
for `ZNyn`; the second diagonal entry is zero for `ZNd`/`ZNy`. In particular, (Y_{0,HL}=Y_{0,LH}=0\), HV ground impedance is independent of LV termination and (Z_m\), and the LV zero-sequence path of `ZNyn` is independent of the HV neutral. `Zd` follows the same implementation but needs an explicit convention: a name without `N` does not itself establish an accessible grounded neutral. An unearthed zigzag must not acquire a galvanic ground merely from its winding connection.

Decoupling follows from ideal zigzag ampere-turn cancellation; the assignment of the measured zero-sequence impedance to (sZ_\ell\), and its IEC correction treatment, is a **model convention**, not a universal physical identity. The PR cites a particular T model but supplies no independently checkable transformer measurement or complete reference edition. This specification tests that stated convention conditionally; actual zigzag impedances must be supplied consistently with it. A separate physical standalone `ZN` device is not established by using a two-winding constructor with a dummy LV winding.

**I3 — Sequence selectivity and fault-current consistency.** Holding other parameters fixed, (\partial Z_1/\partial Z_N=\partial Z_2/\partial Z_N=0\), because balanced positive/negative-sequence currents sum to zero. Consequently balanced three-phase and line-to-line fault results are unchanged. With (Z_2=Z_1\), no current-source injection and zero fault impedance,
\[
I''_{k1}=\frac{\sqrt3\,cU_{LL}}{|2Z_1+Z_0|},\qquad
I''_{k3}=\frac{cU_{LL}}{\sqrt3|Z_1|}.
\]
Here kV/ohm gives kA. On a sole series return path with (2Z_1+Z_0(0)=R+jX\), (R,X\ge0\), increasing (R_N\ge0\) or inductive (X_N\ge0\) cannot increase (I''_{k1}\). For resistance domination, (I''_{k1}\sim cU_{LL}/(\sqrt3 R_N)\); without a zero-sequence return, (I''_{k1}\to0\). Monotonicity is not asserted for resonant, capacitive, or generally reconfigured networks.

**I4 — Primary: dimensional and reference-base covariance.** With physical equipment unchanged, changing arbitrary network base (S_b\to\lambda S_b\) must leave reported ohms, kA and physical voltages unchanged. A physical impedance on winding (e\) must convert as
\[
z^{pu}_e=Z_e S_b/U_{be}^2,\qquad
Z_H=(U_{rH}/U_{rL})^2Z_L
\]
when referred through the ideal transformer ratio (including applicable taps). In particular, an HV zigzag shunt must use its HV physical leakage base even if the auxiliary LV bus nominal voltage changes independently. The same physical neutral law must hold in `runpp_3ph`; changing a per-phase computational base cannot change (V_N=Z_N I_N\).

**I5 — Parallel composition, service status and data representation.** For (p\) identical independent transformer-plus-neutral units,
\[
Y_{0,\rm bank}=pY_{0,\rm unit},\qquad
Z_{0,\rm bank}=(Z_{0T}+3Z_N)/p
\]
for a one-port path. If the bank instead has one **shared** neutral impedance, the different formula (Z_{0T}/p+3Z_{N,\rm shared}\) applies. The creator documents `parallel` as the number of transformers but does not settle whether neutral impedance is per unit or shared; test both with explicitly mapped parameters. The implementation is consistent with a shared-neutral interpretation for simple paths, not automatically with per-unit neutral data.

An out-of-service transformer contributes zero admittance. Explicitly duplicated in-service shunts at a common bus must add, and construction order must not matter. Missing/NaN neutral components represent zero; all four creator APIs must produce the same physical results for equivalent data, including heterogeneous vectorized rows.

**I6 — Passivity and numerical realization.** For nonnegative resistances, nonnegative inductive reactances, (0<s<1\), positive leakage and magnetizing magnitudes, and positive correction factors, the assembled zero-sequence network is passive:
\[
\Re(V^*Y_0V)\ge0,
\]
and grounded driving-point impedances have (\Re Z_{0,\rm th}\ge0\). Finite grounded cases must yield finite, nonnegative current magnitudes. Both inverse-admittance and sparse-factorization solvers must agree within numerical tolerance for well-conditioned configurations. Ungrounded terminal impedances may correctly be infinite; such cases must not be mistaken for finite grounded calculations. Endpoint partitions that create zero denominators are limits, not part of the finite-domain assertion.

## 3. Sources of the invariants

- **I1, I3, I5, I6:** Kirchhoff's current/voltage laws, symmetrical-component definition (I_0=(I_a+I_b+I_c)/3\), passive circuit composition, and conservation of complex power. The factor three and the neutral's terminal placement follow directly from these laws. The T-network matrix above follows by a Schur complement of the internal-node nodal equation; it does not depend on a PR test oracle.
- **I4:** dimensional analysis and ideal-transformer voltage/current ratios. Network MVA is a normalization choice, not equipment data.
- **Library consistency:** [merged zero-sequence builder](https://github.com/e2nIEE/pandapower/blob/0c6a1f0a65273af91dd1d5b854498cd89d457939/pandapower/pd2ppc_zero.py), [base builder](https://github.com/e2nIEE/pandapower/blob/be6b67a8562a448d538f385a43c46718956cbc31/pandapower/pd2ppc_zero.py), [transformer correction implementation](https://github.com/e2nIEE/pandapower/blob/0c6a1f0a65273af91dd1d5b854498cd89d457939/pandapower/build_branch.py), and [single-phase fault-current implementation](https://github.com/e2nIEE/pandapower/blob/0c6a1f0a65273af91dd1d5b854498cd89d457939/pandapower/shortcircuit/currents.py). These establish what is actually calculated, including (Z_2=Z_1\), but are not independent evidence of physical correctness.
- **Standards cited by the code:** IEC 60909-0:2016 for equivalent-voltage-source short-circuit calculation and transformer correction; IEC 60909-4:2021, example 4.4.2, for the PSU example. The code cites section 6.3.3 for its max-only correction policy. The full standards were not independently inspected here; this specification does not certify that policy or derive zigzag leakage partitioning from an unverified clause.
- **I2's conditional model:** the changed builder cites DIgSILENT PowerFactory's *Two-Winding Transformer (3-Phase)* technical reference, Fig. 8.3, and the *J&P Transformer Book*, Appendix 5. The PR gives no edition/page-complete citation; the actual cited diagrams were not verified. Thus (sZ_\ell\) and the uncoupled LV magnetizing path are explicitly treated as the PR's chosen equivalent-circuit convention. The [merged transformer documentation](https://github.com/e2nIEE/pandapower/blob/0c6a1f0a65273af91dd1d5b854498cd89d457939/doc/elements/trafo.rst) identifies T and pi equivalents but does not resolve shared versus individual neutral impedances or `YNyn` neutral-side selection.

The PR's reported PowerFactory comparison supplies neither a reproducible network nor numerical records in the diff, so it is not used as validation evidence.

## 4. Operating domain

The following are representative engineering test ranges, not claimed API bounds or universal equipment limits.

| Input | Representative range and necessary restrictions |
|---|---|
| Bus and rated line-to-line voltages | LV 0.23–1 kV; collector/MV 3.3–33 kV; upstream 66–220 kV. Positive values; normally (U_{rH}>U_{rL}\). Test rated/bus ratios 0.9–1.1 and tap offsets up to ±10%. |
| Transformer `sn_mva` | 0.1–100 MVA; always positive. Earthing/auxiliary examples 0.4–0.63 MVA. |
| Network `sn_mva` | 1, 10, 100 MVA routinely; 0.1 and 123 MVA as normalization checks. |
| `vk_percent`, `vk0_percent` | 2–20%; zero-sequence data can differ substantially by construction. Require (v_k>|v_{kr}|\), (v_{k0}>|v_{kr0}|\) for finite inductive tests. |
| `vkr_percent`, `vkr0_percent` | 0.1–2%, always less than their associated magnitude. Use positive values to avoid the legacy all-zero fallback when constructing analytic oracles. |
| `rn_ohm`, `xn_ohm` | (R_N\ge0,X_N\ge0\); zero, LV examples 0.01/0.02 ohm, MV examples 1–100 ohm, including both simultaneously. Also sample (3|Z_N|/|Z_{0T,e}|=0,0.1,1,10,100\). NaN/absent are data sentinels. No capacitive compensation claims. |
| `si0_hv_partial` | Strictly between 0 and 1; sample 0.1, 0.5, 0.9, and 0.01/0.99 boundary-near cases. Values 0/1 require separate limiting/singularity handling. |
| `mag0_percent`, `mag0_rx` | 100–10,000% (impedance ratio 1–100) and R/X 0–1; include finite magnetizing impedance to expose arm-placement errors. Use 10 and 100,000% as sensitivity cases, not generic measured values. |
| `parallel`, service and PSU flags | Integers 1–4; service True/False; PSU missing/False/True. Valid PSU cases include a synchronous generator linked by `power_station_trafo`, with its required ratings and correction inputs. |
| External-grid short-circuit data | `s_sc_max_mva`, `s_sc_min_mva` 100–10,000 MVA, positive; `rx_max/min` 0.05–0.3, `x0x_max/min` 0.5–3, `r0x0_max/min` 0.05–0.3. Hold these physical data fixed in covariance tests. |
| Solver/model options | Both `case` values; both `inverse_y` values; `lv_tol_percent` 6/10; default SC pi model and the separately identified T path; `shift_degree` 0/±30/330 as appropriate to connection. `fault` 1ph, 2ph, 3ph. No clock digits in `vector_group`. |
| Additional PF inputs | 50/60 Hz with impedances specified at that frequency; asymmetric star-connected loads at 1–20% of rating, phase powers ≥0 and power factor 0.9–1; sufficiently small neutral impedance to maintain a convergent operating point. |

Use `pfe_kw=i0_percent=0` for the minimal SC oracles. Fault resistance/reactance are zero, `use_pre_fault_voltage=False`, and `ip=ith=branch_results=False` unless a specific secondary check below says otherwise. Thermal ratings, peak current, clearing time and arc resistance are outside these invariants.

## 5. Numerical check plan

Use public constructors, public table assignments, `pandapower.shortcircuit.calc_sc`, `pandapower.runpp`, and `pandapower.runpp_3ph`. Do not call `_add_trafo_sc_impedance_zero` to manufacture the expected answer. Extract (Z_0=\texttt{rk0_ohm}+j\texttt{xk0_ohm}\), (Z_1=\texttt{rk_ohm}+j\texttt{xk_ohm}\), and `ikss_ka` from `net.res_bus_sc`. Record every configuration and whether the solve completed.

Build two small fixtures with `create_empty_network`, `create_bus`, `create_ext_grid`, and `create_transformer_from_parameters`:

1. A 20/0.4 kV, 0.63 MVA `Dyn` supplied by a 1,000 MVA short-circuit grid, with (v_k=v_{k0}=6\), (v_{kr}=v_{kr0}=1\), (m=1\), (\rho=0\). This isolates the LV zero-sequence ground path.
2. A 110/20 kV `YNd` main transformer supplied from an external grid, and a 20/0.4 kV, 0.4 MVA grounding transformer on the delta-fed 20 kV bus. The main transformer blocks grid zero sequence at MV. Use (v_k=v_{k0}=4\), (v_{kr}=v_{kr0}=1\) for the grounding transformer. Its auxiliary LV bus has no load/source except when explicitly testing termination sensitivity.

For every new/affected connection, cross both solver choices and both SC cases. Use deterministic parameter values first; add 30 seeded log-uniform samples per connection over positive impedance/rating ranges, uniformly sampling (s\) and dimensionless R/X ratios. Pair each perturbed run with its own baseline. Do not sample invalid combinations with (v_{kr}\ge v_k\). Independently assemble the small analytical sequence circuits in physical units; solve their nodal equations with a dense complex solver and ideal-ratio referral.

- **I1:** On fixture 1, vary pure R, pure X, and combined (Z_N\); assert the exact complex (3Z_N\) increment and unchanged orthogonal component. In fixture 2 replace the grounding unit by `YNd` to obtain an isolated HV one-port. For `Yyn` and `YNy`, test the existing default pi one-port (Z_m+Z_\ell+3Z_N\); supply positive sequence from the opposite ungrounded winding and exclude other zero-sequence returns. For `YNyn`, compare driving-point impedances at **both** ends against the independently assembled T matrix with the full (3Z_N\) in its HV arm. Vary the finite source zero-sequence termination on the opposite side using external-grid zero-sequence data. This makes an incorrect split (s3Z_N,(1-s)3Z_N\) observable. Test `Yzn` with the reverse zigzag one-port convention (grounded LV leakage plus full (3Z_N\)); treat disagreement in the baseline leakage calibration separately from the neutral increment. Sweep (s\), (m\), and (k_T\) so multiplying the external neutral by a leakage partition or correction factor is detectable. For all four ungrounded groups, compare zero/nonzero neutral runs and require equality.
- **I2:** In fixture 2, run `ZNyn`, `ZNd`, `ZNy`, and the implemented `Zd` alias. Check exact (Z_{0,H}=sZ_{\ell,H}+3Z_N\), not just current positivity or order of magnitude. For `ZNyn`, check (Z_{0,L}=Z_{m,L}+(1-s)Z_{\ell,L}\). Change only grid zero-sequence strength at the remote side and verify no cross-side change; alternatively attach an LV external grid with fixed positive-sequence strength while varying only its zero-sequence strength. Compare impedances, since adding a source can alter fault current through (Z_1\). Vary HV neutral R/X and require unchanged LV (Z_0\); vary (m\) and require unchanged HV (Z_0\). With no auxiliary source, require LV `ikss_ka` below (10^{-6}\) kA for `ZNd`/`ZNy`/`Zd`. Attempt `ZN` separately and record its expected rejection as advertised-feature noncoverage; do not quietly substitute `Zd`. Repeat lowercase/mixed-case spellings and check that names ending in clock digits, e.g. `ZNyn11`, produce the explicit validation error.
- **I3:** Reconstruct the metallic earth-fault current from the two reported complex impedances and the chosen voltage factor. In this revision (c_{max}=1.1\) for MV/HV, (c_{min}=1\); below 1 kV, tolerance 10% gives 1.1/0.9 and tolerance 6% gives 1.05/0.95. Cross-check (Z_0\) independently using I1/I2, so a self-consistent but physically wrong impedance cannot pass by current reconstruction alone. On unchanged networks run `fault="3ph"` and `fault="2ph"` before/after neutral changes; require unchanged currents and (Z_1\). Sweep neutral R or X monotonically on isolated return paths; at large R compare with the asymptote, with relative error bounded by the remaining series impedance divided by (3R_N\), rather than demanding exact equality at finite R.
- **I4:** Repeat each grounded case for every network MVA base listed above. Sweep rated/bus voltage ratios independently on HV and LV and use physical HV leakage (Z_{\ell,H}=k_TU_{rH}^2(v_{kr0}+j\sqrt{v_{k0}^2-v_{kr0}^2})/(100S_Tp)\) for the HV shunt oracle. An especially discriminating test changes only the dummy LV bus nominal voltage for `ZNd`, holding winding ratings fixed: the physical HV zero-sequence impedance must not change. Separately exercise tap positions −2/0/+2 at 2.5% per step on each supported tap side, updating the independent referral calculation. This detects a wrong voltage base even when nominal-ratio fixtures pass.
- **I5:** Compare one row with `parallel=p` to (p\) explicit rows at common buses. First use identical per-unit neutrals and record the result under the independent-neutral convention; then set the grouped row's neutral to (Z_N/p\), which maps the independent bank to the implementation's shared-neutral parameterization. Repeat for zero neutrals and for ordinary/zigzag paths. Switch one duplicate out of service and compare against deleting precisely that unit. Shuffle row/index order, mix groups, and ensure multiple HV/LV shunts accumulate correctly. For each of the four creator APIs, use equivalent fully specified zero-sequence data; register a custom transformer standard type with `create_std_type` for the standard-type creators. Cover scalar bulk arguments, vectors, mixed zero/NaN/nonzero rows, omitted columns, and combined R/X. Compare results rather than only stored values.
- **I6:** On all finite grounded cases check finite nonnegative currents and \(\Re Z_0\ge-\epsilon_Z\), and compare solver outputs directly. In two-port cases the independently constructed admittance must have positive-semidefinite Hermitian part; terminal driving-point checks alone are only a necessary public-result check, not a full proof of passivity. Use a weak but finite grounding termination for conditioning tests before taking an open-circuit limit. Check out-of-service near-boundary cases separately; division by an exact zero before multiplication by service status can create NaNs.

**Additional changed branches:** Repeat a valid `YNd` PSU fixture from `iec_60909_4()` with `rn_ohm` omitted, zero and nonzero, and `xn_ohm` zero/nonzero; the expected HV neutral increment is still (3Z_N\), never (k_{ST}3Z_N\). Probe PSU=True with each accepted zigzag group using a correctly linked generator. If this combination is unsupported it should be explicitly rejected; if accepted, require the neutral to be counted once and the LV path to remain independent. Include rows with mixed PSU flags in one vector group.

The default public `calc_sc` path forces `trafo_model="pi"`. To reach its T-model path, first run `runpp(net, trafo_model="t")` and then `calc_sc(..., bus=target, use_pre_fault_voltage=True)`; repeat with `trafo_model="pi"`. Use these as arm-placement/sequence-impedance checks, not as evidence that the single-phase pre-fault-voltage method is otherwise correct. Do not pass a nonexistent `trafo_model` argument to `calc_sc`.

For the separately affected `pf_3ph` path, create asymmetric star loads with `create_asymmetric_load`, cycle the loaded phase through a/b/c, and call `runpp_3ph` for `Dyn`, `YNyn`, and `Yzn`. Reconstruct physical phase voltages from `res_bus_3ph` magnitudes/angles and phase currents from specified complex powers and voltages. Form (V_0=(V_a+V_b+V_c)/3\), (I_0=(I_a+I_b+I_c)/3\) and compare to a physical T-circuit oracle with full (3Z_N\) in the grounded arm and **no IEC SC correction**. For isolated Dyn loads, (V_0/I_0\) with consistent injection sign directly identifies the zero-sequence driving impedance. Repeat MVA-base changes and balanced loads; balanced-load results must be neutral-independent. Verify rejection of new HV zigzag groups in this PF mode. No balanced-network current monotonicity claim is extended to nonlinear load-flow states.

For well-conditioned SC cases, use `rtol=1e-7`, absolute tolerance (10^{-8}\) ohm for each impedance component and (10^{-8}\) kA for currents. Use these same tolerances for solver equality. Open-circuit current checks use (10^{-6}\) kA and allow infinite reported impedances: `results.py` explicitly converts very large pu impedances to infinity. PF comparisons use (10^{-5}\) relative tolerance and (10^{-6}\) pu absolute voltage residual, with converged power mismatch at most (10^{-8}\) MVA. Report conditioning failures separately rather than relaxing tolerances to conceal an order-one discrepancy.

## 6. Do the PR's own tests check these invariants

The two added tests are in [test_1ph.py](https://github.com/e2nIEE/pandapower/blob/0c6a1f0a65273af91dd1d5b854498cd89d457939/pandapower/test/shortcircuit/test_1ph.py). Both test inverse and factorized solves, but do not explicitly compare the two result arrays.

| Invariant | Coverage | Actual assertions and omissions |
|---|---|---|
| I1 | **Partial** | `test_trafo_neutral_earthing_impedance` checks a single ordinary Dyn at 20/0.4 kV, (S_b=1\), `case="max"`: R=0.01 adds 0.03 ohm to (R_0\), X=0.02 adds 0.06 ohm to (X_0\), and the other component is unchanged, all at absolute (10^{-6}\) ohm. No combined R/X, other grounded groups, finite-magnetizing YNyn placement, PSU resistance, or PF test. |
| I2 | **Partial** | `test_zigzag_earthing_transformer` checks positive MV current for ZNyn, MV current within a factor two of YNd, equal MV currents for ZNyn and ZNd within (10^{-6}\) kA, positive ZNyn LV current, and approximately zero ZNd LV current. It does not check exact leakage partition, complex impedance, LV termination sensitivity, magnetizing variation, ZNy/Zd/ZN, taps, or PSU zigzag. |
| I3 | **Partial** | The neutral test checks unchanged positive-sequence R/X for the resistance case and reduced earth-fault current. The zigzag test checks reduced MV current with 5 ohm resistance. Neither verifies the full current formula, the large-R asymptote, X monotonicity, or 2ph/3ph invariance. |
| I4 | **No** in added tests | Both new fixtures have matched rated/bus voltages and (S_b=1\). Existing `test_trafo` compares bases 1 and 123 and fixed reference currents for nine older groups, but has no nonzero neutral and no new zigzag group; `Yzn` is absent. It cannot validate the new base conversions. |
| I5 | **No** | Only single parameter-based creation and `parallel=1` are exercised. No bulk/standard-type creators, common-bus accumulation, service toggle, missing/NaN equivalence, or neutral-bank convention tests. |
| I6 | **Partial** | Both solver routes are subjected to the same limited physical assertions; positive/zero current assertions cover a few cases. There is no direct cross-solver equality, passivity matrix check, near-boundary test, or explicit nonfinite check. |

An existing [test_iec_60909_4_1ph](https://github.com/e2nIEE/pandapower/blob/0c6a1f0a65273af91dd1d5b854498cd89d457939/pandapower/test/shortcircuit/test_iec60909_4.py) compares ten initial currents with reference values (`atol=1e-4`, NumPy's default relative tolerance); its fixture includes a YNd PSU with `xn_ohm=22`. It is relevant regression coverage for the old reactive PSU behavior, not evidence for new resistance or zigzag behavior. Its peak-current assertion is commented out. The PR's statement that existing tests pass has not been independently rerun for this specification.

## 7. Expected verdict

| Invariant | Source-based expectation |
|---|---|
| I1 | **Expected pass** for ordinary Dyn/YNd default SC pi paths, the Yyn/YNy pi neutral increment, and valid YNd PSU R/X extension. **Expected failure** for the stated HV-only neutral convention in finite-magnetizing YNyn: adding (3Z_N\) to total leakage before splitting places (s3Z_N\) in one arm and ((1-s)3Z_N\) in the other. `Yzn` and T-model neutral placement are also suspect. A PSU zigzag appears to count neutral impedance inside `z0_k` and then again in `z0_hv`, giving ((s+1)3Z_N\) on the HV side and an unwanted LV dependence; numerical confirmation is required. |
| I2 | **Expected pass, conditional** for ordinary matched-voltage ZNyn/ZNd/ZNy with the stipulated leakage convention and finite impedances. Decoupling is approximate because of the (10^{20}S_b\) open branch. **Uncertain physical interpretation** for `Zd`. The advertised `ZN` case is expected to raise, not work. The reference-to-parameter mapping for real zigzag equipment remains unverified. |
| I3 | **Expected pass** for sequence selectivity, current reconstruction, and restricted monotonicity on the tested simple SC paths. This can pass even if (Z_0\) itself is physically wrong. |
| I4 | **Expected pass** for MVA-base invariance of ordinary simple paths and new zigzag shunts at matching ratios. **Expected failure** for HV zigzag leakage when HV and LV rated/bus ratios differ: `tap_corr` selects LV but the new shunt is stamped directly at HV. `Yzn` retains an extra `baseMVA` multiplier in its branch admittance, suggesting a base-dependence failure. **Likely failure** in the newly affected PF neutral scaling: leakage uses (3S_b/U^2\), while neutral is added using (S_b/U^2\); T-model placement is an additional issue. These require numerical confirmation. |
| I5 | **Expected pass** for creator propagation, NaN-to-zero behavior and finite-impedance service/accumulation cases. **Conditional** on the bank convention: the code leaves (3Z_N\) undivided by `parallel`, so an independent per-unit-neutral interpretation should fail the explicit-duplicate comparison; a shared-bank-neutral interpretation can pass with the stated mapping. |
| I6 | **Expected pass** in the finite passive interior, even for some incorrectly placed but still passive circuits. **Uncertain** near singular/open configurations and exact partition endpoints. Passing passivity is necessary but not sufficient for correct grounding physics. |

These are predictions from inspected code, not observed test outcomes. The additions appear sufficient for the narrow ordinary-Dyn and matched-voltage ordinary-zigzag examples; broader physical correctness is not established by the PR tests.

## 8. Environment

Run the check in an isolated environment with the **exact merge commit** above, not a moving branch or unspecified installed release. Use a separate installation of base commit `be6b67a8562a448d538f385a43c46718956cbc31` only for legacy regression comparisons; it is expected to reject new zigzag groups. No PowerFactory license is needed for the independent circuit checks.

The merged [pyproject.toml](https://github.com/e2nIEE/pandapower/blob/0c6a1f0a65273af91dd1d5b854498cd89d457939/pyproject.toml) requires Python ≥3.10 and declares NumPy ≥1.26,<2.4, SciPy <1.17, pandas ~=2.3, networkx ~=3.4, packaging ~=25.0, tqdm ~=4.67, deepdiff ~=8.6, geojson ~=3.2, typing_extensions ~=4.9 and pandera ~=0.26.1. Use a compatible Python 3.11 or 3.12 environment, install the checkout with its declared dependencies, and add pytest (test extra declares ~=9.0). NumPy/SciPy provide the independent complex circuit oracle and sparse factorization; optional numba acceleration is not a physical requirement. Record resolved versions and solver options.

The directly relevant regression command is `python -m pytest pandapower/test/shortcircuit/test_1ph.py pandapower/test/shortcircuit/test_iec60909_4.py`. Add a separate test module implementing section 5; executing only those existing files is insufficient. Small CPU-only networks suffice; neither remote compute nor large datasets are needed. This review read the supplied inputs and commit-pinned source/documentation, and did not install pandapower or execute this test plan.

## 9. Formalization feasibility

The core invariants are feasible over exact complex arithmetic in Lean. Define the symmetrical-component transform, prove (I_N=3I_0\), and apply Ohm's law to obtain the neutral contribution. A finite-dimensional nodal model permits proofs of the T-network Schur complement, zigzag block diagonality under its ideal cancellation assumptions, reference-base covariance, parallel-admittance addition, and service-zeroing. For positive real parameters, the stated single-path monotonicity reduces to inequalities for ((R+3R_N)^2+(X+3X_N)^2\). Passivity follows from summing nonnegative resistor losses and eliminating internal passive nodes.

Proofs require explicit nonzero-denominator and nonsingular-grounding hypotheses, positive voltage/rating bases, square-root domain constraints, and an explicit choice of neutral side, bank convention, and zigzag leakage model. Those modeling choices cannot be settled by proving transcribed code. In fact, the YNyn split and PSU-zigzag expressions permit algebraic counterexamples to the intended neutral-placement invariant, so a faithful transcription may prove failure rather than correctness.

The main implementation obstacles are pandas grouping/missing-value semantics, `np.where` branch evaluation, vectorized stamping with repeated indices, tap/base transformations, network connectivity, and the infinity/open-circuit approximation. IEEE floating-point error, conditioning, sparse factorization, and nonlinear three-phase load-flow convergence require separate numerical-analysis arguments; a proof over ideal complex numbers does not establish the section 5 tolerances. Real-transformer saturation, unequal zig/zag sections, and uncertain parameter identification remain external physical assumptions.
