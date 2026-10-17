# ct2107 implementation log (ΔQ1 unified spot check, Cantera#2107)

> Purpose: sentence-by-sentence mapping from the three tracks' check plans to code, and all implementation judgments, for the audit in §6 and the deviation register in §7 of the ΔQ1 protocol. Time slice: written on 2026-10-17 before the smoke run; after the smoke run only the section "Post-smoke corrections" is appended.
> Files: `spotcheck.py`, `harness.py` (identical copies in the four PR directories), `job_smoke.json`, `job_full.json`.
> Material read: the protocol; the sections of `audit/specs/track1.md` / `audit/specs/track2.md` for this PR; `audit/specs/track3_guided/Cantera_cantera_2107/spec.md` §2/§4/§5; from the source tree of merge commit d58a2326: `test/data/*.yaml` (including consistency-cases.yaml), `MargulesVPSSTP.h/.cpp`, `BinarySolutionTabulatedThermo.h/.cpp`; `env/build_cantera.sh`.

## 0. Environment and general conventions

- Native environment: compute node A `env/ctb`, built with scons from the source of merge commit d58a2326 by `env/build_cantera.sh` (`ct.__version__`=4.0.0a1, `__git_commit__`=d58a232). **A native environment is available** (built from merge-commit source, not a release). Phase files are read from `env/cantera-src/test/data` of the same source tree.
- Only the Python interface (`ct.ThermoPhase`) is used; the C++ harness required by the Track 3 plan was not written (see D1 in §6).
- Seed = 20261017 + 1×10 + track number; branch = phase (configuration), N split equally; the verdict, replay and failure-counting rules are as in §0 of the pp3016 log.
- Composition of molality-based phases: m_k is converted to X (X_o=1/(1+M̃_oΣm), X_k=m_kM̃_oX_o, M̃_o=molecular weight/1000 kg/mol) and then set with `TPX`; fixing X fixes m.
- When the plan gives only a relative tolerance the absolute term is 0; when it gives none, relative 1e-6 and absolute 1e-12.
- Determination of a "phase containing IAPWS": any species in the phase whose equation-of-state is liquid-water-IAPWS95 or HKFT (read from input_data at runtime).
- YAML rewriting (synthetic fixtures, cutoff models, standard-concentration basis, species order) uses ruamel.yaml to read source-tree files, change fields, and `ct.ThermoPhase(yaml=…)`.

## 1. Track 1 (`audit/specs/track1.md`)

Phase selection (judgment C1): the plan lists 7 files and no phase names → the phases in these files whose thermodynamic model belongs to the scope changed by this PR (Margules, ideal-molal-solution, Debye-Huckel, HMW-electrolyte, ideal-solution-VPSS, including HKFT/SSVol species) are taken, excluding `HMW-bogus`; 23 branches in all: the 10 phases of debye-huckel-all, debye-huckel-B-dot-ak/beta_ij of thermo-models, HMW-NaCl-electrolyte, HMW-NaCl-HKFT, ideal-molal-aqueous, ethanol-water-rough, molten-salt-Margules, IdealSolnGas-liquid, as well as HMW_NaCl, pdss_hkft, LiKCl_liquid, Li_Liquid, VpssSolidSolutionExample.

| Plan sentence | Code |
|---|---|
| Aqueous T∈[273.15,373.15] K (HKFT up to 573 K), p∈[1,100] atm, I∈[1e-4,6] mol/kg (HMW up to ~6.1); molten-salt Margules T∈[700,1400]; solid solution VPSS T∈[300,1000] | `t1_domain`, `t1_sampler` |
| A: read cv, cp, β, κ_T, v; compare β, κ_T with the central differences of v(T±h,p), v(T,p±δ) (h=1e-4·T, δ=1e-4·p; rtol 4e-6 for phases containing IAPWS, 1e-7 for others) | `A.beta=…`, `A.kappaT=…` of `t1_eval` |
| A: c_v = c_p − T v β²/κ_T (κ_T>0), c_v=c_p when κ_T=0 (library convention) | `A.Mayer…` (same rtol as A), `A.kappa=0 => cv=cp` (default tolerance) |
| A: isochoric path (T,p)→(T+h, p+(β/κ_T)h), 2–3 Newton steps to pull v back, c_v≈T·Δs/ΔT, rtol 1e-5 | `isochoric_newton` with 3 steps; one run on each of the positive and negative sides, taking the central difference (judgment C3); residuals go into records |
| B: compare chemical_potentials with standard_gibbs_RT·RT + RT·log(γX), rtol 1e-10 | `B.mu=mu0+RTln(gX)` |
| B: compare the solute term with the hand-computed γ^□m_k for m_k=X_k/(X_oM̃_o) (use the C++ getMolalityActivityCoefficients for γ^□ if Python has no interface) | Python has no γ^□ interface; `activities` (the molality-based activity a_k=γ^□m_k/m° of MolalityVPSSTP) is taken as γ^□m_k and compared with γ_kX_k, default tolerance (judgment C4) |
| C: central difference of μ_k(T±h) against −partial_molar_entropies (rtol 1e-6); ΣX_k s̄_k against entropy_mole (rtol 1e-12); central difference of h(T±h) against cp_mole; s(p±δ) and −v(T±h) compared with each other | `C.*`; the last two items have no tolerance in the plan → default |

Implementation judgments (Track 1):
- C1 Phase selection: see above.
- C2 Composition: the plan gives only the interval of I → for 1:1 Na⁺=Cl⁻, m=I (drawn uniformly), the other species are fixed at the values of state point 0 of consistency-cases.yaml (DH: H⁺ 1.05e-8, OH⁻ 1.3765e-6, NaCl(aq) 0.98492…; HMW: H⁺ 1.05e-5, OH⁻ 1e-5); the composition of ideal-molal, ethanol-water, molten salt, solid solution and HKFT (mole-fraction phases) is fixed at state point 0. The plan gives no pressure for molten salt and solid solution → fixed at 1 atm. ethanol-water is assigned to the interval of the "aqueous" class (not listed separately in the plan). The HKFT class (pdss_hkft, HMW-NaCl-HKFT) uses T up to 573.15 K.
- C3 Isochoric c_v: the truncation error of the one-sided T·Δs/ΔT is about h/(2T)·(T s_TT/s_T) ~5e-5, exceeding the plan's 1e-5; as a "central difference", the isochoric path is run once on each side.
- C4 γ^□m_k uses `activities` (see table above); this has one layer less independence than the C++ interface, recorded as a deviation candidate.
- C5 The absolute term for the κ_T and β comparisons is 0 (the plan gives only relative). When p is below the saturation pressure (e.g. 573 K, 1 atm) the PDSS of water may give a vapour branch or raise an error; it is computed anyway and an error counts as a failure.

## 2. Track 2 (`audit/specs/track2.md`)

Branches (17, the case phases of the models changed by this PR in consistency-cases.yaml): LiKCl_liquid, ethanol-water-rough, Li_Liquid, the 10 phases of debye-huckel-all, HMW_NaCl, ideal-molal-aqueous, VpssSolidSolutionExample, pdss_hkft.

| Plan sentence | Code |
|---|---|
| Usual operating domain: LiKCl T 900–1400, P 1–20 atm; ethanol–water T 298–308, X_eth 0.3–0.5; Li(L) T 300–500, P 1–10 atm; DH/HMW T 300–360, P 1–5 atm, m 0.01–9.4; HKFT T 300–360; "it is suggested to extend the pressure to 1–100 atm" | `t2_sampler`: T, X_eth, m as in the table; P is always 1–100 atm (judgment D1) |
| I1: set P±dP with TPX, compare the central difference of chemical_potentials with partial_molar_volumes; dP=max(1e-3P,1e3 Pa); rtol 1e-6, atol 1e-12·max\|v̄\|, 1e-5 for IAPWS phases | `I1.dmu/dP=vbar` (applies to the phases the plan lists, i.e. all 17) |
| I2: second central difference (v(T+h)−2v+v(T−h))/h², h=0.5–1 K; first central difference of cp with respect to P; rtol 1e-4 | h is 0.5 K (judgment D2); `I2.dcp/dP=-T d2v/dT2` |
| I2: SSVol species c_p,k°(T,P)−c_p,k°(T,P_ref)=−2c·T·(P−P_ref) | only for molar-volume-temperature-polynomial species (Li(L)2); the data is a cubic polynomial, so v°''=2c+6dT is used (judgment D3); P_ref=1 atm |
| I3: algebraic comparison of cv, cp, v, β, κ_T, rtol 1e-10; c_v=c_p when κ_T=0; c_v≤c_p | `I3.*` |
| I3: activities against activity_coefficients·X, rtol 1e-12 (MolalityVPSSTP subclasses) | DH, HMW, ideal-molal branches |
| I3: Gibbs–Duhem: composition perturbation δX (ΣX=1, \|δX\|~1e-6), compute ΣX_kΔμ_k, relative to max\|Δμ_k\| <1e-6 | δX_k=1e-6·X_k·(r_k−Σ_jX_jr_j) (r random per point); Δμ uses the central difference of X±δX (judgment D4) |

Implementation judgments (Track 2):
- D1 Pressure: the "suggested extension to 1–100 atm" in the plan's domain counts per the protocol ("suggestions count too") → P is uniform over 1–100 atm for all phases (including ethanol-water, HKFT, ideal-molal). For phases whose T is not listed in the table: ideal-molal uses its state-point range 300–340 K, VpssSolid 300–500 K.
- D2 The second-difference step is the lower end of the interval, 0.5 K.
- D3 The plan writes SSVol as a quadratic polynomial (v''=2c), but the actual data is a 4-term cubic polynomial → the exact second derivative of that polynomial is used; density-temperature-polynomial species (Li(L)) get no explicit formula and only go through the general I2.
- D4 Gibbs–Duhem: the perturbation is proportional to X_k, which guarantees positivity and ΣδX=0; the second-order term of a one-sided difference is of relative size ~1e-6·X, the same order as the 1e-6 criterion → central differences are used to cancel the second-order term (an implementation judgment about the plan's "Δμ").

## 3. Track 3 (spec.md §5)

Branches (26): the 10 phases of debye-huckel-all, HMW_NaCl, ideal-molal cutoff models none/poly/polyexp (3; polyexp is the original file), ideal-molal standard concentration unity/species-molar-volume (2; solvent-molar-volume is the original file), pdss_hkft, VpssSolidSolutionExample, Li_Liquid, LiKCl_liquid, ethanol-water-rough, margules_synthetic (ethanol–water interaction parameters randomised), Redlich-Kister-LiC6, Redlich-Kister-complex, anode (tabulated binary), tab_simplified (a simplified tabulated fixture with zero background and unity standard concentration). Mandatory points: stress molalities {0.1,1,6,9.4} for DH/HMW (298.15 K, 1 atm); all nodes, adjacent midpoints, x=0, 1, {1e-2,1e-6,1e-10,1e-14,1e-300} of the two tabulated fixtures, and symmetric trace states in which the trace of the reference species {1e-2,…,1e-300} is given explicitly, once for each of the two species orderings.

| Plan sentence | Code |
|---|---|
| Water/HKFT/DH/HMW: T 280–360 K, P 0.101325–10 MPa; NaCl 1e-6–1e-2 mol/kg (0.1,1,6,9.4 only as stress); interior composition 1e-6–1−1e-6; Margules/SSVol: ethanol–water 280–330, Li 300–500, molten salt 900–1400, P 1–20 atm; tabulated 280–350 K, 1–10 atm, x∈[0,1] | `t3_sampler` |
| Central difference h_T/T∈{1e-4,5e-5,2.5e-5}, h_P/P∈{1e-3,5e-4,2.5e-4}, first compare adjacent estimates; non-convergence = indeterminate | `deriv3`: a difference between the last two levels exceeding the tolerance is recorded as a failure, otherwise the finest level is used for the comparison |
| Algebraic identities \|a−b\| ≤ 1e-5 + 1e-12·max; derivatives relative 1e-6 (analytic)/1e-5 (water/HKFT), absolute floors cp 1e-3, volume 1e-10, Maxwell 1e-9, κ 5e-14, β 1e-9 | `alg`, `alg_v`, `deriv3`, `FLOOR` |
| I1: g=h−Ts, h=u+Pv, μ_k=h̄_k−Ts̄_k, f=ΣX_kf̄_k (v,h,s,g,c_p), c_p=h_T=Ts_T, s_P=−v_T, v̄_k=∂_Pμ_k | `t3_I1` |
| I1: SSVol with two kinds of polynomials, nonzero curvature, at P_r and far from P_r, compared with the analytic SSVol formula | `t3_ssvol` (the two species of Li_Liquid each have one kind of polynomial; P is the sampled value, P_r=1 atm) |
| I1: setTemperature then setPressure, reversed order, and setState_TP give consistent results | `t3_setorder`: approximated in Python with `TP=(T,None)`/`TP=(None,P)` (judgment E2) |
| I1: setCutoffModel none/poly/polyexp of ideal molal; solvent-enriched state and both sides of the transition | three branches by changing YAML `cutoff.model`; X_o drawn uniformly in [1e-6,1−1e-6] |
| I1: GibbsExcess subclasses including Margules (with/without excess volume) and Redlich–Kister | corresponding branches |
| I2: five PDSS (constant volume, two SSVol, IAPWS water, HKFT); dVdT/dVdP against molarVolume differenced externally; properties of the parent phase/water unchanged before and after the derivative calls | phase level: β, κ_T against the differenced volume (`deriv3`); h, s, v, μ unchanged after repeated derivative calls (`==`); PDSS-level dVdT/dVdP unreachable (D1) |
| I2: compare ρ, the sum of partial molar volumes, β, κ_T, c_p−c_v with the independent differenced volume | `I2.v=Mbar/rho`, v=ΣXv̄ of I1, `I2.cp-cv=-T vT^2/vP` |
| I2: when v_P≠0 solve v(T±h,P±)=v(T,P) in P with a bracketed scalar root, compare c_v with the isochoric derivatives of u and s; root residual <1e-11 relative volume; not bracketable = indeterminate | `_root_P` + `t3_I2` (compared after convergence of the two step levels) |
| I2: the fallback for v_P=0, v_T≠0 is only a convention check; the direct cv_mole of HKFT is recorded as unsupported | the convention is only recorded; the phase-level cv_mole of pdss_hkft is callable (PDSS level unreachable), recorded only on error |
| I3: a=Xγ (1e-10 relative + 1e-20 absolute); chemical potential difference (in RT units) 1e-9 absolute; ideal molal, HMW, 10 DH; standard concentration unity/species-molar-volume/solvent-molar-volume | `t3_I3`; the standard concentration changes YAML only for ideal-molal (DH/HMW have no such option) |
| I4: g^E,h^E,s^E,v^E computed from independent polynomials, compared with phase property minus standard state/ideal mixing; P=P*, 2P*, 20P*, 10 MPa; each coefficient includes zero values and both signs | `margules_excess` (formula from the MargulesVPSSTP.h documentation, parameters read from input_data as SI values), `t3_I4`; coefficient randomisation is in the margules_synthetic branch (each coefficient is 0 with probability 1/4, otherwise uniform over ±(data magnitude)) |
| I4: normalised composition path Σ_kX_kD_k=0; J_kj by differencing ln γ after normalising n_j↦n_je^{±ε} (ε=1e-4,1e-5,1e-6) | `I4.Gibbs-Duhem…`, `I4.sum_k X_k dlngamma_k/dln n_j=0` (tolerance 1e-6 relative + 1e-8 absolute) |
| I5: DH solvent and each solute partial molar c_p, volume compared with the chemical-potential/enthalpy derivatives | volume is in I1 (v̄=∂_Pμ), c_p is in `t3_I5` (c̄_p=∂_T h̄) |
| I6: anode and simplified fixtures; nodes, midpoints, endpoints, traces; two orderings; finiteness and the algebraic sum of I1; on the simplified fixture μ_tab−μ_ref=h_tab−Ts_tab (independent interpolation) and s=x s_tab−R ln(1−x); change the temperature and return to the same composition to check the cache | `t3_I6`; the interpolation uses a self-written linear interpolation (end values outside the table, same convention as the source) |

## 4. Unimplemented items (Track 3, unreachable through the Python interface; all deviation candidates)

- Comparison of `MargulesVPSSTP::getdlnActCoeffdT / getd2lnActCoeffdT2 / getdlnActCoeffds / getdlnActCoeffdlnN / _diag / dlnX_diag` with finite differences (I4).
- Derivative comparison for the four modes of `WaterProps::ADebye(T,P,ifunc)`, and the saturation-pressure clamping convention (I5).
- Direct comparison of dVdT/dVdP of `providePDSS(k)`, and the direct cv of PDSS_HKFT (I2).
- Separate calls of `setTemperature` / `setPressure` (approximated with the None form of TP), `setMolalities` (converted from X), `getMolalityActivityCoefficients`/`getMolalities` (m computed by hand; γ^□ can only be obtained from a/m, which would be a tautology, so no check is made).
- The synthetic mixtures required by the plan in I2 such as "water+rigid" and "water+HKFT into Margules" that are not in the test data; the ternary synthetic phases of I4 (zero/one/two interactions) were not constructed.
- The separate check in I6 of nonzero background/unequal standard concentration/the "explicit correction" at elevated pressure was not done (the plan says to flag only the dimensional problem).

## 5. Time and memory estimate (speculation)

A single TPX + properties: non-aqueous phases ~0.02–0.1 ms, IAPWS phases ~0.2 ms (empirical values). Track 1 has ~40 states per point → 5–20 ms; Track 2 ~30 states → the same order; Track 3 has ~200–1000 states per point (isochoric root finding up to ~800) → 0.05–0.3 s.

| Track | N=100 (including mandatory) | N=10⁴ |
|---|---|---|
| 1 | ~10 s | ~5 min |
| 2 | ~10 s | ~5 min |
| 3 | ~2 min (the ~350 mandatory tabulated points are fast) | ~30 min (upper estimate 1 h) |

Peak memory ~300 MB (Cantera + numpy; margules_synthetic builds a new phase per point, and old phases are freed with their references). The JobSpec declares 1 GiB, and the full run is serial.

## 6. Deviation candidates (for the orchestrating agent to rule on)

- D1 Track 3 requires a C++ harness; this script uses only Python, and the items listed in §4 were not tested. To make up for it, a small harness would have to be compiled on compute node A with `env/cantera-src/build/lib/libcantera_shared.dylib` + the source includes (building needs separate authorisation).
- D2 γ^□m_k of Track 1 B uses `activities` in place of the C++ interface (C4).
- D3 Track 2 Gibbs–Duhem uses central differences (D4).
