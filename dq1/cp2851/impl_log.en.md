# cp2851 implementation log (ΔQ1 unified spot check, CoolProp#2851)

> Purpose: sentence-by-sentence mapping from the three tracks' check plans to code, and all implementation judgments, for the audit in §6 and the deviation register in §7 of the ΔQ1 protocol. Time slice: written on 2026-10-17 before the smoke run; after the smoke run only the section "Post-smoke corrections" is appended.
> Files: `spotcheck.py`, `harness.py` (four identical copies), `job_smoke.json`, `job_full.json`.
> Material read: the protocol; the sections of `audit/specs/track1.md` / `audit/specs/track2.md` for this PR; `audit/specs/track3_guided/CoolProp_CoolProp_2851/spec.md` §2/§4/§5; `HelmholtzEOSMixtureBackend.cpp` of merge commit 54a2e327 and of v8.0.0 (ae81610e); `env/setup_pca_py.sh`.

## 0. Environment

- compute node A `env/py2`: CoolProp 8.0.0 wheel, `__gitrevision__`=ae81610e. GitHub compare: ae81610e is **291 commits newer and 0 older** than merge commit 54a2e327 (including this PR); the two `calc_first_saturation_deriv` functions (the one with SatL/SatV and the argument-less delegating version) are character-for-character identical in v8.0.0 and the merge commit (same md5). **Treated as having a native environment**, but it is not a build of the merge commit itself, and the rest of the code differs by 291 commits; recorded as deviation candidate D0.
- Seed = 20261017 + 6×10 + track number; verdict/replay/failure rules as in §0 of the pp3016 log; when the plan gives only a relative tolerance the absolute term is 0.
- A new AbstractState is created for every flash (Track 3 states "separate state instances" explicitly; Tracks 1 and 2 do the same, so that the previous state does not serve as an initial value and affect replay).

## 1. Track 1 (`audit/specs/track1.md`)

Branches (k=5): A_ME_Q0, A_ME_Q1, A_R32R125_Q0, A_R32R125_Q1, B_pure_limit_ME.

| Plan sentence | Code |
|---|---|
| Methane&Ethane (z_CH4∈[0.1,0.9], T∈[120,240] K, at least 5 K away from the critical point/cricondentherm); also one group of R32&R125 or R32&R1234yf (z∈[0.2,0.8], T∈[230,330]) | `T1_MIX`, `t1_sampler`; of the two alternatives R32&R125 (listed first) is taken |
| Exclusion of near-critical points | `t1_accept`: accepted only if the QT flashes of both Q=0 and Q=1 succeed at T−5, T and T+5, otherwise redrawn (≤60 times), judgment E1 |
| A: AS.update(QT_INPUTS,Q,T), first_saturation_deriv(iT,iP); also build states at T±h, flash to obtain p, central difference 2h/(p₊−p₋); Richardson extrapolation with h=0.2,0.1,0.05 K; do both Q=0 and Q=1; relative 1e-4 | `fd_slope`, `richardson` (S(h)=S₀+ah²+bh⁴ through three points to obtain S₀), `A.dT/dp_sat=…` |
| B: compare with the value of Methane at the same T, Q=0, for ε=1e-2,1e-3,1e-4, and check that the difference falls linearly with ε | `t1_eval_B`, `_linear_rate`: d/ε equal at adjacent levels (protocol default tolerance) + d decreasing; both the Q=0 and Q=1 mixture values against pure methane at Q=0 |
| C: compare first_saturation_deriv(iHmolar,iT) with the central difference of h(T±h) along the Q=0 line, relative 1e-4 | only in the Q=0 branch; h uses the three levels of A + Richardson (judgment E2) |
| Also check that Q=0.5 raises an exception | `expect_raise` at every point |

Implementation judgments (Track 1):
- E1 "At least 5 K from the critical point/cricondentherm" comes with no way to decide it → the operational definition is that the flashes of both saturation lines at ±5 K succeed; rejected samples are redrawn (equivalent to uniform on the allowed set).
- E2 The finite-difference step of C is not written in the plan → the 0.2/0.1/0.05 K + Richardson of A is reused; the bulk h uses `keyed_output(iHmolar)` of the QT state (at Q=0 the saturated liquid of composition z).
- E3 The "linear decrease" of B has no tolerance → by the protocol default (relative 1e-6), d/ε at adjacent levels is judged equal. The O(ε²) term makes d/ε differ at the order of ε, so this check **will probably be rejected widely because of the tolerance**: not certain, for the orchestrating agent to rule on whether to classify it as "plan wording" (see the report).
- E4 T of B: drawn uniformly in [120,240], requiring that both pure methane and the ε=1e-2 mixture can be flashed at T±5 K (same as E1); the critical temperature of methane Tc≈190.6 K, so the effective interval is about 120–185 K.

## 2. Track 2 (`audit/specs/track2.md`)

Branches (k=11): {ME, NM (Nitrogen&Methane), R32R125, PB (Propane&n-Butane), CM (CO2&Methane)} × {Q0,Q1} + I2_pure_limit_ME. Mandatory points: z₁∈{0.05,0.3,0.5,0.7,0.95} (discrete values of the plan) × every branch, T=0.75·T_c,mix.

| Plan sentence | Code |
|---|---|
| Five mixtures, z₁∈{0.05,…,0.95}, T taken as 0.6–0.9 times the mixture critical temperature | `t2_sampler`: z₁∈U[0.05,0.95], T=U[0.6,0.9]·T_c(z); T_c uses `T_critical()`, and when it fails (several critical points) falls back to `T_reducing()` and records it (judgment F1) |
| I1: central difference 2ΔT/(p(T+ΔT)−p(T−ΔT)), ΔT=0.01,0.02,0.05 K with Richardson; tolerance 1e-3; Q=0 and Q=1 | `I1.dT/dp_sat=Richardson FD` |
| Additionally: first_saturation_deriv(iP,iT)·first_saturation_deriv(iT,iP)=1 | `I1.dP/dT*dT/dP=1` (no tolerance → default) |
| I2: on the z₂→0 sequence (1e-2,1e-3,1e-4) compare with Methane, the deviation falls linearly with z₂ | `t2_eval_I2`, T=U[0.6,0.9]·T_c(methane), Q=0 and Q=1; linearity judged as in E3 of Track 1 |
| I3: compare first_saturation_deriv(iHmolar,iT) with (h(T+ΔT)−h(T−ΔT))/(2ΔT) (along Q=0), rtol 1e-3; an exception or a returned two-phase value is recorded as a finding | Q=0 branch; step and Richardson as in I1; exception → record a "finding" and count it as a failure per the protocol (judgment F2) |

Implementation judgments (Track 2): F1 see above; F2 the plan treats "raises an exception" as a finding, but the protocol provides that an exception is not a violation → counted as a failure and recorded separately; F3 I2 uses only Methane&Ethane (the plan only names Methane); F4 the three Richardson step ratios are unequal (2, 2.5), so an exact extrapolation through three points with a two-term polynomial in h², h⁴ is used.

## 3. Track 3 (spec.md §5)

Branches (k=10): {ME, MEP (Methane&Ethane&Propane), Water, nPropane, R410A} × {Q0,Q1}. Mandatory points: all discrete values of the table in §4 (ME 140–220 K every 10 K × methane fraction 0.1/0.3/0.5/0.7/0.9; MEP 160–220 K × two compositions; Water 300/400/500/550; n-Propane 230/270/300/330; R410A 250/280/310; both endpoints of each), plus an interior-quality rejection test point at 180 K, (0.7,0.3). Random points: ME T∈U[140,220], z∈U[0.1,0.9]; MEP T∈U[160,220], composition taken with equal probability from the two groups; pure fluids uniform over the T interval of the table.

| Plan sentence | Code |
|---|---|
| \|a−b\| ≤ r·max(\|a\|,\|b\|)+ε_abs; slope ε_abs=1e-8·T/p; property derivatives ε_abs=1e-8·F*/X* | `rel_ok`, `_close3` |
| I1: QT differences with δT=0.05,0.025,0.0125 K, adjacent reference values converge within relative 5e-4, then compare with r=1e-3; separately reproduce the PR's 0.05 K, 0.5% | `t3_I1`; pure fluids also with r=1e-3? The plan's I1 writes r=1e-3 only for mixtures; for pure fluids I2 uses 1e-4 (the code uses 1e-4 for I1 of pure fluids, judgment G1) |
| I1: PQ_INPUTS, [T(p+δp)−T(p−δp)]/(2δp) with δp/p=1e-4,5e-5,2.5e-5 | same function |
| I1 weighted residual: liquid/vapour compositions and saturated densities; build homogeneous states with specify_phase + PT_INPUTS; the central density agrees with the equilibrium branch; finite-difference ln φ_i to get a_{k,i}; \|A_T S+A_p\|/(\|A_T S\|+\|A_p\|)<2e-3; a root switch invalidates that level; composition normalisation 1e-10, max\|ln(f_L/f_V)\|<1e-6 | `t3_I1_weighted`; weights: the bubble line uses the vapour composition, the dew line the liquid composition (judgment G2); root-switch criterion: the density of the differenced state deviates from the equilibrium branch by >5% |
| I2: pure-fluid endpoint S against TΔv/Δh (saturation keyed output) r=1e-10; difference r=1e-4; bubble/dew slope at the same temperature r=1e-8; R410A only the expression check, base vs merge isolated builds compared | `t3_I2`; for the base vs merge nodes of R410A there is no isolated parent-commit build → recorded only, not run |
| I3: iDmolar, iHmolar, iSmolar, iUmolar and the mass versions, with respect to iT and iP; compared with differences of keyed_output at fixed q,z (mixtures 1e-3, pure 1e-4); chain rule with first_partial_deriv (1e-10); reciprocals and identities 1e-10 absolute | `t3_I3` |
| I4: internal quality 0.1/0.5/0.9 (180 K, 0.7/0.3), all four forms must raise; (iT,iDmolar) must raise at the endpoints | rejection tests at the mandatory points; (iT,iDmolar) checked at every point |
| I5: binary reversed order, 6 permutations of the ternary; p_sat, slope, rearranged phase composition r=1e-7, composition absolute 1e-9 | `t3_I5` (mixtures only) |
| I6: compare second_saturation_deriv(F,iP,iP) (F=T,Dmolar,Hmolar,Smolar,Umolar) with the pressure central difference of first_saturation_deriv(F,iP); reference convergence 0.1%; r=5e-3, ε_abs=1e-8F*/p²; unsupported calls reported separately | `t3_I6`; unsupported → records |

Implementation judgments (Track 3):
- G1 The pure-fluid endpoint I1 difference comparison uses 1e-4 (the plan's I2 writes "independently use the flash finite differences above at r=1e-4").
- G2 Weights of A_T, A_p: the definition in spec Section 1 (outside the readable range) cannot be seen by the implementer → the standard form derived from the phase-equilibrium total differential + Gibbs–Duhem (bubble line y, dew line x). Not certain.
- G3 The tolerance for the homogeneous-state density agreeing with the equilibrium branch is not given in the plan → default relative 1e-6.
- G4 F* is taken as max(1, \|F central value\|).

## 4. Time and memory estimate (speculation; compute node A probes: mixture QT flash ~0.7 ms, new AbstractState ~1–3 ms)

| Track | Per point | N=100 (including mandatory) | N=10⁴ |
|---|---|---|---|
| 1 | ~50–80 ms (including acceptance flashes) | ~10 s | ~15 min |
| 2 | ~30–50 ms | ~10 s | ~10 min |
| 3 | ~0.5 s (~250 state builds + flashes) | ~2 min (mandatory ~150 points) | ~1.5 h |

Peak memory ~200–300 MB per process. The full run has the three tracks in parallel (declared 2 GiB, 3 cores).

## 5. Deviation candidates

- D0 The environment is the v8.0.0 wheel (contains this PR, the changed functions are character-for-character identical), not a build of the merge commit.
- D1 The "linear decrease" of Track 1 B and Track 2 I2 is judged with the protocol default tolerance and is expected to be rejected widely (E3).
- D2 The R410A base-vs-merge check of Track 3 I2 was not run (no isolated parent-commit build).
