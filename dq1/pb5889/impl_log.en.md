# PyBaMM#5889 ΔQ1 spot-check implementation log

> Purpose: record the sentence-by-sentence mapping of the three tracks' check plans to `spotcheck.py`, and all implementation judgments, for review under §6/§7 of the ΔQ1 protocol.
> Structure: §0 environment → §1 Track 1 → §2 Track 2 → §3 Track 3 → §4 time and memory → §5 deviation candidates → §6 post-smoke corrections → §7 to be ruled on (classification of smoke rejections).
> Update rule: after the smoke run only script errors may be fixed; changing the plan mapping must be recorded in protocol §7 and appended to §5.
> References: upstream: the ΔQ1 protocol, the #5889 sections of `audit/specs/track1.md`/`track2.md`, `audit/specs/track3_guided/pybamm-team_PyBaMM_5889/spec.md` §2/4/5; downstream: `dq1/outputs/pb5889_t{1,2,3}.json`.

## §0 Environment and general conventions

- Node: compute node A, `env/py2` (pybamm 26.10.1.0, PyPI). Verified: `Chayambuka2022.py`, `data/D_p.csv`, `data/D_n.csv`, `lithium_ion_parameters.py`, `settings.py`, `expression_tree/interpolant.py` are byte-for-byte identical to merge commit eda03ed. Treated as the native environment of the merge commit.
- The `Chayambuka2022.py` of parent commit c235b90 differs from the merge commit in only one place at line 162 (c_max taken from `Initial concentration in electrolyte`); Track 2's "compare with the pre-fix version" uses this one difference to rebuild the parent-commit function in the same environment (`base_NVPF_diffusivity`, parameter override), without installing the parent commit separately.
- Seed 20261017 + 12×10 + track number. Each track has 1 branch (a single NVPF phase, a single linear interpolation, no options). Mandatory points are not counted in N.
- The data table is read independently by the script: `np.loadtxt(D_p.csv, delimiter=",", skiprows=1)` (the file inside the installed package, verified identical to the merge commit).
- For each (C_max, E) configuration a fresh `ParameterValues` is built, `update` is called once after building and it is not modified again; the internal cache is cleared every 2000 evaluations.

## §1 Track 1

| Plan text | Code |
|---|---|
| Usual operating domain sto∈[0.01,0.99], T∈[273,323] (the function does not depend on T) | `t1_sample`, linearly uniform |
| sto grid 0.01…0.99 in steps of 0.01, T=298.15 | mandatory points `grid_sto=…` (99 points, running all checks) |
| Compare `pv.evaluate(f(Scalar(sto), Scalar(T)))` with `numpy.interp(sto*15320, x_csv, y_csv)`, rtol 1e-10 (default linear) | `t1_api_vs_numpy_interp_rtol1e-10` |
| `pv.update({c_max: c})` changed to another c; at fixed c_s, f(c_s/c) is unchanged, rtol 1e-10 | `t1_cmax_invariance_c=…`: c ∈ {7660, 30640} |
| Secondary (magnitude): D_p ∈ [1.9e-17, 2.6e-15]; D_p ≤ 1e-16 for sto≥0.3 | `t1_magnitude_*` |
| Macroscopic corroboration: BasicDFN default parameters, 2C discharge to 2.0 V, capacity about 1.18 mAh | mandatory point `macro_BasicDFN_2C`: `Experiment(["Discharge at 2C until 2.0 V"])`, read the last value of `Discharge capacity [A.h]` |

Implementation judgments:
1. The value of "another c" is not given in the plan: 0.5× and 2× the default c_max (7660, 30640) are taken, with fixed c_s = sto×15320; the lookup argument c_s/c×c restores c_s, and c_s lies within the table.
2. "About 1.18 mAh" takes abs tolerance 0.005 mAh, half of the last place of 3 significant digits. The plan notes itself "not a physical invariant, only a regression reference", and it is still included per protocol §2.
3. The default interpolator is verified to be linear (`Interpolant(..., interpolator="linear", extrapolate=True)`), so the alternative comparison in the plan is not needed.

## §2 Track 2

| Plan text | Code |
|---|---|
| Usual operating domain: sto 0.01–0.99, T fixed at 298.15, rate C/50–3C, cutoff 2.0 V | `t2_sample`: sto and C_rate linearly uniform, T fixed |
| Operating domain "BasicDFN and the DFN/SPMe models" | only BasicDFN: `pybamm.sodium_ion` provides only `BasicDFN` (`__init__.py` exports only it); no sodium-ion DFN/SPMe is available |
| I1: for each node c_k of the table, sto_k=c_k/15320, compare with the table value y_k, rtol 1e-10 | mandatory points `node_k=0…33` |
| I1: random points between nodes compared with `numpy.interp`, rtol 1e-10 | sampled points `t2_I1_api_vs_numpy_interp_rtol1e-10` |
| I1: check positivity at sto∈{0,1} | mandatory points `endpoint_sto=0/1` |
| I1 (invariant statement): 1.7e-17 ≤ D_p ≤ 2.51e-15 on sto∈[0,1], and D_p>0 | `t2_I1_bounds_*`, `t2_I1_positive` checked at sampled points and the two endpoints |
| I1 (invariant statement): the negative electrode likewise D_n(sto)=L_n(sto·14540) | sampled points `t2_I1_negative_Dn_vs_numpy_interp_rtol1e-10` (`D_n.csv`) |
| I2: BasicDFN + `Experiment([f"Discharge at {C} until 2.0 V"])`, read the last value of `Discharge capacity [A.h]`; C/50 compared with the pre-fix version ≤2%; 2C compared with 1.18 mAh rtol 2%; monotonicity | mandatory point `anchors_C/50_1C_2C_parent`: four Experiment solves (C/50 merge, C/50 parent, 1C, 2C) |
| I2: capacity non-increasing with rate, Q(C/50) ≥ Q(1C) ≥ Q(2C) | the order of the three is checked at the anchors; each sampled rate C is then compared with the anchors: Q(C/50) ≥ Q(C); for C≤1 Q(C) ≥ Q(1C), otherwise ≤; for C≤2 Q(C) ≥ Q(2C), otherwise ≤ |
| I3: take `Positive particle radius [m]` and compute τ_p = R_p²/D_p(sto=0.5), which should be of the same order as 1800 s or larger | mandatory point `tau_p_sto=0.5` |

Implementation judgments:
1. The magnitude bound and the negative-electrode isomorphism written in invariant I1, although not listed word for word in the "numerical check method", are included per protocol §2 because that section refers to the invariant as "I1".
2. The capacity at sampled rates uses a fast route: the same BasicDFN is discretised only once, `Current function [A]` is set to `[input]`, `solve([0, 5·3600/C])`, and the model's own `Minimum voltage [V]` event terminates it at `Lower voltage cut-off [V]`=2.0 V (the parameter-set value, equal to the plan's 2.0 V). The Q(C/50), Q(1C), Q(2C) used in the anchor comparison are also computed separately by the same fast route (`anchor_fast`), to avoid mixing the two routes in a comparison; the anchor check named by the plan uses the Experiment route.
3. The monotonicity has no tolerance, and is compared loosely with the default rel 1e-6, abs 1e-12.
4. "Same order or larger" is taken as τ_p ≥ 180 s (not lower than 1/10 of 1800 s).
5. For the "pre-fix version (parent commit)" see §0: only this function is overridden.

## §3 Track 3

| Plan text | Code |
|---|---|
| Public route: build a new parameter object after each override; D with rtol 1e-10, atol 0, required finite | `pv(C,E)` caches a fresh object per configuration; `Dp()` |
| Independent reference: take the two enclosing rows and evaluate by the two-point formula, without using another Interpolant; end segments explicitly extrapolated | `L(c)` |
| I1: with the default C, at the 34 nodes and, for the 33 intervals, the division points 0.01/0.25/0.5/0.75/0.99 | mandatory points `knot_k` (34), `interval_k_frac_f` (165), running all the same checks as the sampled points |
| I1: draw 256 concentrations uniformly on the table (`default_rng(5889)`) | replaced by the protocol's N points and protocol seed: `t3_sample` c∈[x₀, x_last] linearly uniform |
| I1: independent values at θ=0.3, 0.5, 0.9: 5.6673280002e-17, 8.3317028568e-17, 3.5626719998e-17 | mandatory points `theta_*` |
| I1: samples passed in as `pybamm.Vector(theta_array)`, compared elementwise with the scalar results | mandatory points `vector_samples` (all sampled points), `vector_knots_intervals` |
| I1 model wiring: `BasicDFN().param.p.prim.D(Scalar(c), Scalar(298.15))` compared with L(c) (interior points) | sampled points `t3_I1_model_adapter_vs_L(c)_rtol1e-10` |
| I1 boundary: c/C = 0, τ/2, τ, 2τ, 1−2τ, 1−τ, 1−τ/2, 1 compared with L(C·clip(c/C,τ,1−τ)) | mandatory points `adapter_c/C=…` (τ = settings `D__c_s` = 1e-10) |
| Direct call at θ=0, x₀/(2C), (x_last/C+1)/2, 1 compared with the linear extrapolation of the end segments; θ=0, 1 about 2.7385714278e-15, 1.7387466637e-17 | mandatory points `direct_theta=…` |
| I2: repeat the default-C sampling with E=500/1000/1500, compare at fixed θ | sampled points `t3_I2_E=…_equals_default` |
| I2: C=12256/15320/18384 × the same E, θ=c/C, equal to L(c) at nodes and division points in the common admissible interval [x₀, 12256(1−τ)] | sampled points (and the node/division-point mandatory points) falling in this interval: `t3_I2_C=…_E=…_vs_L(c)` |
| I3: D>0, lying between the two enclosing table values (margin 1e-12×the larger) | `t3_I3_positive`, `t3_I3_within_bracketing_ordinates` |
| I3: at each interior node c=x_i±δ_i, δ_i=1e-6·min(adjacent spacing), compared with the one-sided linear formula, \|D−d_i\| ≤ \|slope\|δ_i+1e-12 d_i | mandatory points `knot_i_delta_±1` (32×2) |
| I3: flux corollary g=0, ±1e9, −Dg² ≤ 0 | `t3_I3_flux_sign_g=…` |
| I4: results at each listed temperature compared with the 298.15 K result, rtol 1e-12 | `t3_I4_T=273.15/323.15_equals_298.15` |

Implementation judgments:
1. Branches: the plan states "a single phase, a single linear interpolation, no charge/discharge/hysteresis/temperature branches"; the 33 segments, extrapolation and clamping cases are covered by mandatory points; the N points are no longer split equally by segment (k=1).
2. The model wiring uses the `param` (`LithiumIonParameters`) of BasicDFN, evaluated with the parameter object of the default configuration.
3. The vector check uses one Vector call over all sampled points (regenerating the same sequence with the same seed as the sampling), compares elementwise, and keeps the rtol 1e-10 of D as tolerance.
4. The negative-electrode HC function (which the plan says can serve as a control) was not done: the plan writes "can serve as a control" and gives no check.

## §4 Expected time and peak memory (speculation, not run)

| Track | Per point | 10⁴ points | Basis |
|---|---|---|---|
| t1 | ~10 ms (3 evaluate calls) | ~2 min + grid + 1 Experiment (~10–30 s) | |
| t2 | ~0.5–1 s (1 fast BasicDFN solve + 2 evaluate calls) | ~1.5–3 h; plus 4 Experiment anchors (C/50 is about 50 h of simulated time, the solve takes a dozen seconds or so) | BasicDFN is discretised only once |
| t3 | ~40 ms (about 13 evaluate calls) | ~7 min + about 450 mandatory points + vector comparison (10⁴ further scalar evaluations) | |

Peak memory: about 1–1.5 GB; the JobSpec declares 3 GiB. If t2 takes more than 1.4 s per point, N is reduced per protocol §3 and recorded in §7.

## §5 Deviation candidates

- The Track 2 operating domain mentions DFN/SPMe while the sodium-ion module has only BasicDFN, so only BasicDFN was done (candidate: record in §7 "the configuration named in the plan does not exist in the library").
- Track 2 "pre-fix version" is rebuilt by overriding a single function in the same environment instead of installing the parent commit (the parent commit differs in the source only in this one place, verified).

## §6 Post-smoke corrections (2026-10-17)

- None. The only rejection in the first round is a problem with how the plan wrote its bound (§7); the script is faithful to the plan and was not changed.

## §7 To be ruled on (classification of smoke rejections, script not changed)

1. **(P) Track 2 I1 "on sto ∈ [0,1], min(table) ≈ 1.7e-17 ≤ D_p ≤ max(table) = 2.51e-15 (the end extrapolation amount is small)" does not hold at the endpoint sto=0**. The only rejection in the first round: mandatory point `endpoint_sto=0.0`, observed D_p = 2.7386e-15 > 2.51e-15 (exceeding by 9.1%); replayable (recomputation with the same input is bitwise identical). Cause: the first two nodes of `D_p.csv` in the merge commit are (131.578947, 2.51e-15) and (592.105263, 1.71e-15), and PyBaMM's linear interpolation extrapolates linearly outside the table; at c=0, 2.51e-15 + 0.80e-15×131.58/460.53 = 2.7386e-15, which agrees with the observation to the last digit. At the other end, sto=1 (c=15320), extrapolation gives about 1.74e-17, and the plan's lower bound 1.7e-17 is written exactly from the extrapolated value; the upper bound, however, takes the table maximum, so the two ends are treated inconsistently. Where the plan is wrong: "≤ max(table) = 2.51e-15" contradicts "end extrapolation" in the same sentence and the library's default linear extrapolation (it should be ≈2.74e-15). The sampled points (sto∈[0.01,0.99], c ≥ 153 > 131.6) are unaffected, and only this one mandatory point will be rejected in the full run.
2. No (S) or (V) class was seen for this PR in the first round.
