# cp3339 implementation log (ΔQ1 unified spot check, CoolProp#3339)

> Purpose: sentence-by-sentence mapping from the three tracks' check plans to code, and all implementation judgments, for the audit in §6 and the deviation register in §7 of the ΔQ1 protocol. Time slice: written on 2026-10-17 before the smoke run; after the smoke run only the section "Post-smoke corrections" is appended.
> Files: `spotcheck.py`, `harness.py` (four identical copies), `job_smoke.json`, `job_full.json`.
> Material read: the protocol; the sections of `audit/specs/track1.md` / `audit/specs/track2.md` for this PR; `audit/specs/track3_guided/CoolProp_CoolProp_3339/spec.md` §2/§4/§5; `TransportRoutines.cpp::viscosity_rhosr` of merge commit 5418cec5, parent commit afcc7909 and v8.0.0 (ae81610e); the PR tests in `src/Tests/CoolProp-Tests.cpp` of the merge commit; the TRANSPORT.viscosity blocks of the 8 fluid JSON files under `dev/fluids/` of the merge commit; `env/setup_pca_py.sh`.

## 0. Environment (important: there is no native merge-commit environment)

- `env/py2` on compute node A/compute node B is CoolProp 8.0.0 (ae81610e). GitHub compare: ae81610e is an **ancestor of merge commit 5418cec5, 63 commits behind, and does not contain this PR**; its `viscosity_rhosr` is character-for-character identical to that of parent commit afcc7909 (centre hard-coded to 2). In addition, the output of `get_fluid_param_string(name,"JSON")` in 8.0.0 has no TRANSPORT section (probe: only ANCILLARIES/EOS/INFO/STATES).
- Conclusion: **there is no native environment for this PR on the nodes**; it can only be built from the merge-commit source. Per the common instructions, no build was attempted by the implementer. env/py2 can serve as a "parent-commit equivalent" baseline (`--base-python env/py2/bin/python`), but a cross-build comparison is not a bitwise comparison within the same build.
- The script probes `honours_x_crossover` on the current interpreter (clone x_c=2 and 12, compare at the PR test state), written to `versions.env_probe` in the JSON; if the current CoolProp is 8.0.0, Track 1 main 2 and Track 3 I1 will be "rejected" on a large scale, which is a wrong environment and not a measured result.
- The JobSpec assumes that a merge-commit build is installed at `env/cp3339` (does not exist; the orchestrating agent must decide whether to build it: compile the CoolProp Python extension with `uv pip install` from git 5418cec5 source, C++17 + CMake, about ten-odd minutes, speculation).
- Seed = 20261017 + 7×10 + track number; verdict/replay/failure rules as in §0 of the pp3016 log.
- Clone registration: following the PR test (change NAME/CAS, clear ALIASES, REFPROP_NAME="N/A", delete the INCHI/SMILES/CHEMSPIDER/2DPNG fields), `add_fluids_as_JSON("HEOS", …)`; `OVERWRITE_FLUIDS` is turned on and each purpose overwrites a fixed name repeatedly, so that the number of registrations within the process does not grow with the number of points.
- x is computed from public quantities: x = ρ·R·(τ·∂α^r/∂τ − α^r)/rhosr_critical (`rhomolar/gas_constant/tau/dalphar_dTau/alphar`); rhosr_critical, C, c_liq, c_vap are read from the fluid JSON.

## 1. Track 1 (`audit/specs/track1.md`)

Branches: 8 fluids (R1233zd(E), R1234yf, R1234ze(E), R124, R152A, R22, R245fa, R32), 1250 points each. Mandatory points: single-phase densities (found by bisection; skipped if out of reach) corresponding to x∈{x_c±0.01, ±0.05, ±0.5} (x_c∈{2,12}) at T=300 K.

| Plan sentence | Code |
|---|---|
| T∈[230,420] K, ρ from 1 mol/m³ to the saturated liquid density, single phase | `t1_sampler`, `sample_rho`: for T<Tc ρ∈U[1, ρ_L,sat(T)], redrawn if it falls in the two-phase region; for T≥Tc the upper limit is ρ_L,sat(230 K) (judgment H1) |
| Main 1: compute PropsSI("V","T",T,"Dmolar",ρ,fluid) once with the parent commit and once with the merge commit, compare float.hex | the baseline process (`--base-python`, `--serve-base` mode) returns the hex of PropsSI; `main1.parent==merge bitwise` |
| Main 2: register clones with x_c=a, b via add_fluids_as_JSON; compute x from smolar_residual() and rhomolar() and rhosr_critical; compute the theoretical ratio from C, c_liq, c_vap and compare with the PropsSI ratio, rtol 1e-12 | a=2 (the library value), b=12 (the PR test value); `main2.eta_a/eta_b=…` |
| Grid x∈{x_c±0.01, ±0.05, ±0.5} | as mandatory points |
| Minor: read c_vap[0] at ρ=1e-3 mol/m³ and see η/η_dilute→1 (η_dilute by extrapolation with decreasing ρ) | η₀=2η(1e-3)−η(2e-3); `minor.c_vap[0]==0`, `minor.eta/eta_dilute->1` (default tolerance) |

Implementation judgments (Track 1):
- H1 For T>Tc, "up to the saturated liquid density" is undefined → the upper limit is the largest saturated liquid density in the interval, ρ_L(230 K).
- H2 The plan writes x=−ρs^r/(ρs^r)_crit and in the same sentence x=ρR(τα_τ−α)/(ρs_r)_crit; whether the signs of the two expressions agree depends on the sign convention of smolar_residual → the latter expression (consistent with the source) is taken; at each point it is also computed with smolar_residual, and a disagreement is only recorded.
- H3 Main 2 takes two fixed clones a=2, b=12 (rather than generating a clone per point according to x), and the grid points cover the transition region through mandatory points.

## 2. Track 2 (`audit/specs/track2.md`)

Branches: 8 fluids, 1250 points each.

| Plan sentence | Code |
|---|---|
| T from the triple point to 1.5 T_c, ρ from 1e-6 ρ_c to the saturated liquid density, single phase | `t2_sampler` (for T≥Tc the upper limit is the saturated liquid density at the triple point, judgment H4); points beyond the EOS Tmax are computed anyway and an error counts as a failure |
| I1: the merge commit and the parent commit are bitwise equal in PropsSI("V",…) at (T,ρ) | same as main 1 of Track 1 |
| I2: clones a=2, b∈{1,3,6,12}; for x≤min(a,b)−δ or ≥max(a,b)+δ (δ=0.2), \|η_a/η_b−1\| ≤ \|f_L−f_V\|·e^{−100δ} in magnitude | "in magnitude" is judged as within a factor of 10: bound = 10·\|f_L−f_V\|·e^{−20} (judgment H5) |
| I2: for any x, L≡ln(1+(η/η₀−1)/C) ∈ [min(f_L,f_V), max(f_L,f_V)]; η₀ by ρ→0 extrapolation or computed from the JSON dilute parameters | η₀ is taken by linear ρ→0 extrapolation (two points 1e-6ρ_c, 2e-6ρ_c, judgment H6); both a=2 and each b are checked; boundary tolerance default |
| I3: take η at ρ=ρ_c·{1e-3,1e-4,1e-5,1e-6}, extrapolate linearly to 0 and compare with the smallest-ρ value, rtol 1e-6; read c_V[0] | the extrapolation uses the smallest two points; the slope sequence is only recorded ("approaches linearly with ρ" is qualitative, judgment H7) |

## 3. Track 3 (spec.md §5)

Branches: 8 fluids, 1250 points each. Mandatory points: for each fluid one I4 dilute-limit group (T∈{250,300,min(400,Tmax)} × ρ/ρ_c∈{1e-6,1e-8,1e-10,1e-12} × a∈{2,12,0}); Q=0/1 of QT at T/Tc∈{0.70,0.85,0.95}; supercritical T/Tc∈{1.01,1.05} × ρ/ρ_c∈{0.1,0.5,1,1.5,2,3} (within Tmax); the PR state R1233zd(E) 300 K, 1e4 mol/m³. These mandatory points are also cross-checked against PropsSI.

| Plan sentence | Code |
|---|---|
| 12 temperatures × 24 logarithmic pressures from 1 kPa to the pressure limit (PT_INPUTS), practical region 250 K–min(450 K,Tmax), 0<p≤min(20 MPa,pmax), ρ≲3ρ_c | random points: T uniform, p log-uniform (the plan states logarithmic explicitly); redrawn if ρ>3ρ_c (judgment H8); the plan's regular grid is replaced by uniform sampling (protocol §3) |
| I1: a C=0 clone gives η₀; a=x+d, d∈{−0.2,−0.05,−0.01,0,0.01,0.05,0.2}, plus a=2, 12; an independently written overflow-safe logistic; relative error ≤5e-12 | `t3_state` |
| I1: recover g=log[1+(η/η₀−1)/C], g(x+d)+g(x−d)=L+V, absolute 1e-10·max(1,\|L\|,\|V\|) | d∈{0.01,0.05,0.2} |
| I1: central difference of the a derivative with h=1e-5, 5e-6, relative 1e-5; in the degenerate case the absolute tolerance scales with η₀max(1,CE\|L−V\|) | at a=x; in the degenerate case (\|L−V\|<1e-6 or C=0) absolute tolerance 1e-10·η₀·max(1,CE\|L−V\|) (the coefficient 1e-10 is judgment H9) |
| I1: tail probe at a=x±8, inner exponential overflow should give weight 0 rather than a non-finite viscosity | `I1.tail` |
| I2: compare the baseline and merge builds for 8 fluids; bitwise for the same build, otherwise relative 1e-13; library fluids against the a=2 clone 1e-14 | cross-build judged at 1e-13, bitwise equality only recorded; clone comparison 1e-14 |
| I3: η finite and positive; bounds at the two endpoints (slack 5e-12); η≥η₀ for x≥0; negative x reported separately; with an equal-polynomial clone and the C=0 control, changing the crossover has no effect | the I3 section of `t3_state`; the equal-polynomial clone sets c_liq=c_vap |
| I4: 250/300/min(400,Tmax) K, ρ/ρ_c=1e-6…1e-12, a=2,12,0: compare with the finite-density formula at 5e-12; the last two points against the limit formula at 1e-8; for a=2,12 also check η/η₀≈1; for a=0 check the non-1 limit | `t3_dilute`; "check the non-1 limit" for a=0 = comparison with the limit formula, with the ratio also recorded |
| Use a fresh process for each fluid to isolate registration | not followed: single process, overwrite registration with fixed names (judgment H10, deviation candidate) |

## 4. Time and memory estimate (speculation)

PropsSI/AbstractState about 1 ms per call; clone registration (deep copy of ~70 KB JSON + parsing and building the fluid) ~5–20 ms.

| Track | Per point | N=100 (including mandatory) | N=10⁴ |
|---|---|---|---|
| 1 | ~15 ms + baseline round trip | ~5 s | ~3–5 min |
| 2 | ~30 ms | ~5 s | ~5–10 min |
| 3 | ~0.2–0.3 s (13 overwrite registrations per point) | ~1–2 min | ~45–60 min |

Peak memory: main process ~200–400 MB + baseline process ~150 MB; if overwrite registration leaks, it grows linearly (unverified) → declared 2 GiB, peak to be checked during the smoke run. The full run is serial.

## 5. Deviation candidates (for the orchestrating agent to rule on)

- D0 **There is no native merge-commit environment**; `env/cp3339` (merge commit) must be built before running; the baseline uses env/py2 (a viscosity_rhosr equivalent to the parent commit, but a different build). Track 1 main 1 and Track 2 I1 require "bitwise" equality, and a cross-build bitwise comparison may produce false rejections due to compilation differences.
- D1 Track 3: regular grid → uniform/log-uniform sampling; single-process overwrite registration instead of a fresh process per fluid (H10).
- D2 Track 2 I2 "in magnitude" is judged as within a factor of 10 (H5).
