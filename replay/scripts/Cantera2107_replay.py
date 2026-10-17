"""Cantera2107_replay.py -- native replay for Cantera PR #2107 (merge d58a23265414e1789536cf5332d9fb29de095a43).

Checks, against the Cantera Python module built from the merge commit (or any later commit),
the formulas and theorems of lean/PhysAudit/Cantera2107.lean:

  A. MargulesVPSSTP (phase "ethanol-water-rough" from test/data/thermo-models.yaml, embedded below),
     2000 random states T in [280, 360] K, P in [1, 500] atm, X_ethanol in [0.02, 0.98]:
       A1  ln(activity_coefficients) vs the Lean transcription `lnActCoeff`            (formula replay)
       A2  partial_molar_volumes vs v0 + Lean `vbarExcess`                              (formula replay)
       A3  d mu_k / dP (central difference) vs partial_molar_volumes                    (margules_dmu_dP_eq_vbar)
       A4  -d mu_k / dT (central difference) vs partial_molar_entropies                 (margules_dmu_dT_eq_neg_sbar)
       A5  partial_molar_enthalpies vs mu_k + T * s_k                                    (margules_hbar_eq_mu_plus_T_sbar)
       A6  thermal_expansion_coeff * v vs dv/dT (central difference); isothermal_compressibility == 0;
           cv_mole == cp_mole (kappa_T == 0 branch)                                      (margules_beta_correct,
                                                                                         margules_kappa_zero_branch)
       A7  pre-PR counterexample: by RT_lnActCoeff_affine, d mu_k^old / dP = v0_k, so the pre-PR
           violation equals -vbarExcess_k; reported at T = 298.15 K, X = (0.5, 0.5), P = 100 atm.
  B. IdealSolnGasVPSS with an IAPWS-95 water species + constant-volume ethanol (embedded YAML),
     2000 random liquid states T in [285, 360] K, P in [1, 200] atm:
       B1  cv_mole vs cp - T V beta^2 / kappa_T using the phase's own beta, kappa_T, V   (isg_mayer_test_tautology)
       B2  beta * V vs dV/dT, -kappa_T * V vs dV/dP (central differences)
       B3  cv_mole vs T * (ds/dT) along a numerically constructed isochore             (isg_cv_eq_isochoric_dudT)
  C. MargulesVPSSTP with IAPWS-95 water + ethanol and excess-volume parameters (kappa_T != 0 case):
       C1  cv_mole vs isochoric T * ds/dT                                               (margules_cv_eq_isochoric_dudT)
       C2  d mu_k / dP vs partial_molar_volumes
  D. PDSS_SSVol (molar-volume-temperature-polynomial) single-species ideal-solution-VPSS phase,
     2000 random states T in [300, 900] K, P in [2, 500] atm:
       D1  standard_cp_R vs cp0/R - T dP V''(T) / R  (post-PR formula)                  (ss_cp_eq_T_dsdT)
       D2  standard_cp_R vs T * d(standard_entropies_R)/dT (central difference)
       D3  pre-PR formula cp0/R - T dP V''(T) (missing /R) evaluated on the same states (ss_old_cp_fails)

Dependencies: Python 3, numpy, cantera built from d58a232 or later (env `ct2107`).
Run:  python Cantera2107_replay.py > Cantera2107_replay.json
If cantera cannot be imported, or lacks the PR's API, the script prints {"status": "SKIP", ...}.
"""

import json
import math
import sys

import numpy as np

RESULT = {"module": "Cantera2107", "merge_commit": "d58a23265414e1789536cf5332d9fb29de095a43"}

try:
    import cantera as ct
except Exception as exc:  # pragma: no cover
    print(json.dumps({**RESULT, "status": "SKIP", "reason": f"cannot import cantera: {exc!r}"}))
    sys.exit(0)

RESULT["cantera_version"] = getattr(ct, "__version__", None)
RESULT["cantera_git_commit"] = getattr(ct, "__git_commit__", None)

# ---------------------------------------------------------------- constants (ct_defs.h)
GAS_CONSTANT = 6.02214076e26 * 1.380649e-23
ONE_ATM = 1.01325e5
SMALL_NUMBER = 1e-300
N_POINTS = 2000
RNG = np.random.default_rng(2107)

# ---------------------------------------------------------------- embedded phases
ETHANOL_WATER_YAML = """
phases:
- name: ethanol-water-rough
  thermo: Margules
  elements: [H, C, O]
  species: [ethanol, water]
  interactions:
  - species: [ethanol, water]
    excess-enthalpy: [-2449 J/mol, -1613 J/mol]
    excess-entropy: [-14.45 J/mol/K, -11.64 J/mol/K]
    excess-volume-enthalpy: [-12.707 cm^3/mol, 0.476 cm^3/mol]
    excess-volume-entropy: [-0.030 cm^3/mol/K, 0.0 cm^3/mol/K]
species:
- name: ethanol
  composition: {C: 2, H: 6, O: 1}
  thermo: {model: constant-cp, T0: 298.15 K, h0: -277.7 kJ/mol, s0: 159.86 J/mol/K,
           cp0: 112.0 J/mol/K}
  equation-of-state: {model: constant-volume, molar-volume: 58.35 cm^3/mol}
- name: water
  composition: {H: 2, O: 1}
  thermo: {model: constant-cp, T0: 298.15 K, h0: -285.83 kJ/mol, s0: 69.91 J/mol/K,
           cp0: 75.3 J/mol/K}
  equation-of-state: {model: constant-volume, molar-volume: 18.07 cm^3/mol}
"""
# Lean `ethanolWater`, SI units (J/kmol, J/kmol/K, m3/kmol, m3/kmol/K); A = ethanol (index 0), B = water (1)
EW = dict(iA=0, iB=1, HEb=-2449000.0, HEc=-1613000.0, SEb=-14450.0, SEc=-11640.0,
          VHEb=-0.012707, VHEc=0.000476, VSEb=-0.00003, VSEc=0.0)
EW_V0 = np.array([0.05835, 0.01807])

IAPWS_PHASES_YAML = """
phases:
- name: isg-water-ethanol
  thermo: ideal-solution-VPSS
  elements: [H, C, O]
  species: [H2O(L), ethanol]
- name: margules-water-ethanol
  thermo: Margules
  elements: [H, C, O]
  species: [H2O(L), ethanol]
  interactions:
  - species: [ethanol, H2O(L)]
    excess-enthalpy: [-2449 J/mol, -1613 J/mol]
    excess-entropy: [-14.45 J/mol/K, -11.64 J/mol/K]
    excess-volume-enthalpy: [-12.707 cm^3/mol, 0.476 cm^3/mol]
    excess-volume-entropy: [-0.030 cm^3/mol/K, 0.002 cm^3/mol/K]
species:
- name: H2O(L)
  composition: {H: 2, O: 1}
  thermo: {model: constant-cp, T0: 298.15 K, h0: -285.83 kJ/mol, s0: 69.91 J/mol/K,
           cp0: 75.3 J/mol/K}
  equation-of-state: {model: liquid-water-IAPWS95}
- name: ethanol
  composition: {C: 2, H: 6, O: 1}
  thermo: {model: constant-cp, T0: 298.15 K, h0: -277.7 kJ/mol, s0: 159.86 J/mol/K,
           cp0: 112.0 J/mol/K}
  equation-of-state: {model: constant-volume, molar-volume: 58.35 cm^3/mol}
"""

SSVOL_COEFFS = [0.0125, 1.0e-6, 1.0e-8, 1.0e-11]  # m3/kmol, m3/kmol/K, m3/kmol/K2, m3/kmol/K3
SSVOL_CP0 = 30.0e3  # J/kmol/K
SSVOL_YAML = f"""
phases:
- name: ssvol-test
  thermo: ideal-solution-VPSS
  elements: [Li]
  species: [Li(L)]
species:
- name: Li(L)
  composition: {{Li: 1}}
  thermo: {{model: constant-cp, T0: 298.15 K, h0: 0.0 J/kmol, s0: 29.1e3 J/kmol/K,
            cp0: {SSVOL_CP0} J/kmol/K}}
  equation-of-state:
    model: molar-volume-temperature-polynomial
    data: {SSVOL_COEFFS}
"""


# ---------------------------------------------------------------- Lean transcriptions
def pair_sum(p, X, g0, g1, k):
    """Lean `pairSum` for a single binary interaction."""
    XA, XB = X[p["iA"]], X[p["iB"]]
    val = -1.0 * (XA * XB) * (g0 + g1 * XB) - XA * XB * XB * g1
    if k == p["iA"]:
        val += XB * (g0 + g1 * XB)
    if k == p["iB"]:
        val += XA * (g0 + g1 * XB) + XA * XB * g1
    return val


def ln_act_coeff(p, T, P, X, k):
    rt = T * GAS_CONSTANT
    g0 = (p["HEb"] - T * p["SEb"]) / rt + (P - ONE_ATM) / rt * (p["VHEb"] - T * p["VSEb"])
    g1 = (p["HEc"] - T * p["SEc"]) / rt + (P - ONE_ATM) / rt * (p["VHEc"] - T * p["VSEc"])
    return pair_sum(p, X, g0, g1, k)


def vbar_excess(p, T, X, k):
    return pair_sum(p, X, p["VHEb"] - T * p["VSEb"], p["VHEc"] - T * p["VSEc"], k)


def rel(a, b, floor=0.0):
    return abs(a - b) / max(abs(a), abs(b), floor, 1e-300)


class MaxTracker:
    def __init__(self):
        self.max_rel = 0.0
        self.worst = None
        self.n = 0

    def add(self, r, info):
        self.n += 1
        if not math.isfinite(r) or r > self.max_rel:
            self.max_rel = r if math.isfinite(r) else float("inf")
            self.worst = info

    def out(self, tol):
        return {"n": self.n, "max_rel_err": self.max_rel, "tol": tol,
                "pass": bool(self.n > 0 and self.max_rel <= tol), "worst": self.worst}


def set_tpx(ph, T, P, X):
    ph.TPX = T, P, X


# ---------------------------------------------------------------- A. ethanol-water Margules
def part_a():
    out = {}
    ph = ct.ThermoPhase(yaml=ETHANOL_WATER_YAML, name="ethanol-water-rough")
    trk = {k: MaxTracker() for k in ["A1", "A2", "A3", "A4", "A5", "A6_beta"]}
    kappa_max = 0.0
    cv_cp_max = 0.0
    for _ in range(N_POINTS):
        T = RNG.uniform(280.0, 360.0)
        P = RNG.uniform(1.0, 500.0) * ONE_ATM
        xe = RNG.uniform(0.02, 0.98)
        X = np.array([xe, 1.0 - xe])
        set_tpx(ph, T, P, X)
        Xn = ph.X
        lng = np.log(ph.activity_coefficients)
        vbar = ph.partial_molar_volumes
        mu = ph.chemical_potentials
        sbar = ph.partial_molar_entropies
        hbar = ph.partial_molar_enthalpies
        V = ph.volume_mole
        beta = ph.thermal_expansion_coeff
        kappa = ph.isothermal_compressibility
        cv, cp = ph.cv_mole, ph.cp_mole
        state = {"T": T, "P": P, "X_ethanol": xe}
        for k in range(2):
            lean_lng = ln_act_coeff(EW, T, P, Xn, k)
            trk["A1"].add(abs(lng[k] - lean_lng) / max(abs(lean_lng), 1e-3), {**state, "k": k})
            trk["A2"].add(rel(vbar[k], EW_V0[k] + vbar_excess(EW, T, Xn, k)), {**state, "k": k})
        dP = max(1e-3 * P, 1e4)
        set_tpx(ph, T, P + dP, X); mup = ph.chemical_potentials
        set_tpx(ph, T, P - dP, X); mum = ph.chemical_potentials
        dT = 1e-2
        set_tpx(ph, T + dT, P, X); muTp = ph.chemical_potentials; Vp = ph.volume_mole
        set_tpx(ph, T - dT, P, X); muTm = ph.chemical_potentials; Vm = ph.volume_mole
        for k in range(2):
            trk["A3"].add(rel((mup[k] - mum[k]) / (2 * dP), vbar[k]), {**state, "k": k})
            trk["A4"].add(rel(-(muTp[k] - muTm[k]) / (2 * dT), sbar[k], 1e3), {**state, "k": k})
            trk["A5"].add(rel(hbar[k], mu[k] + T * sbar[k], 1e3), {**state, "k": k})
        trk["A6_beta"].add(rel(beta * V, (Vp - Vm) / (2 * dT), 1e-12), state)
        kappa_max = max(kappa_max, abs(kappa))
        cv_cp_max = max(cv_cp_max, abs(cv - cp))
    out["A1_lnActCoeff_formula"] = trk["A1"].out(1e-10)
    out["A2_vbar_formula"] = trk["A2"].out(1e-12)
    out["A3_dmu_dP_eq_vbar"] = trk["A3"].out(1e-6)
    out["A4_neg_dmu_dT_eq_sbar"] = trk["A4"].out(1e-6)
    out["A5_hbar_eq_mu_plus_Tsbar"] = trk["A5"].out(1e-10)
    out["A6_beta_times_V_eq_dVdT"] = trk["A6_beta"].out(1e-6)
    out["A6_kappa_T_max_abs"] = kappa_max
    out["A6_cv_minus_cp_max_abs"] = cv_cp_max
    out["A6_note"] = ("kappa_T == 0 and beta != 0 for this phase: cv_mole returns cp_mole by "
                      "convention; no isochoric path exists (margules_no_isochore_when_kappa_zero)")

    # A7: pre-PR counterexample at the Lean witness state
    T, X = 298.15, np.array([0.5, 0.5])
    rows = []
    for P_atm in (1.0, 100.0, 500.0):
        P = P_atm * ONE_ATM
        set_tpx(ph, T, P, X)
        vbar = ph.partial_molar_volumes
        dP = 1e4
        set_tpx(ph, T, P + dP, X); mup = ph.chemical_potentials
        set_tpx(ph, T, P - dP, X); mum = ph.chemical_potentials
        dmu_new = (mup - mum) / (2 * dP)
        vE = np.array([vbar_excess(EW, T, X, k) for k in range(2)])
        # RT_lnActCoeff_affine: mu_old = mu_new - (P - P_ref) * vE  =>  d mu_old / dP = d mu_new / dP - vE
        dmu_old = dmu_new - vE
        rows.append({"P_atm": P_atm, "vbar_native": vbar.tolist(),
                     "dmu_dP_native": dmu_new.tolist(),
                     "dmu_dP_prePR": dmu_old.tolist(),
                     "prePR_violation_dmu_dP_minus_vbar": (dmu_old - vbar).tolist(),
                     "lean_vbarExcess": vE.tolist(),
                     "prePR_excess_term_in_mu_RT_lng": ((P - ONE_ATM) * vE).tolist()})
    out["A7_prePR_counterexample"] = {
        "input": {"T": T, "X": [0.5, 0.5]}, "rows": rows,
        "lean_vbarExcess_ethanol_expected": -9.40625e-4}
    return out


# ---------------------------------------------------------------- isochoric cv helper
def isochoric_cv(ph, T, P, X, h=0.05):
    """T * ds/dT along V = const, with P(T) solved by Newton using finite-difference dV/dP."""
    set_tpx(ph, T, P, X)
    V0 = ph.volume_mole
    beta, kappa = ph.thermal_expansion_coeff, ph.isothermal_compressibility
    slope = beta / kappa

    def s_on_isochore(Tn):
        Pn = P + slope * (Tn - T)
        for _ in range(12):
            set_tpx(ph, Tn, Pn, X)
            V = ph.volume_mole
            dp = max(1e-6 * Pn, 1.0)
            set_tpx(ph, Tn, Pn + dp, X)
            dVdP = (ph.volume_mole - V) / dp
            step = (V - V0) / dVdP
            Pn -= step
            if abs(step) < 1e-9 * Pn:
                break
        set_tpx(ph, Tn, Pn, X)
        return ph.entropy_mole, abs(ph.volume_mole - V0) / V0

    sp, ep = s_on_isochore(T + h)
    sm, em = s_on_isochore(T - h)
    set_tpx(ph, T, P, X)
    return T * (sp - sm) / (2 * h), max(ep, em)


def mayer_block(ph, name, Trange, Prange, n, check_mu):
    trk = {k: MaxTracker() for k in ["taut", "beta", "kappa", "iso", "mu"]}
    iso_resid = 0.0
    kappa_zero = 0
    for _ in range(n):
        T = RNG.uniform(*Trange)
        P = RNG.uniform(*Prange) * ONE_ATM
        xw = RNG.uniform(0.3, 0.98)
        X = np.array([xw, 1.0 - xw])
        set_tpx(ph, T, P, X)
        cv, cp = ph.cv_mole, ph.cp_mole
        V = ph.volume_mole
        beta, kappa = ph.thermal_expansion_coeff, ph.isothermal_compressibility
        state = {"T": T, "P": P, "X_water": xw}
        if kappa == 0.0:
            kappa_zero += 1
            continue
        trk["taut"].add(rel(cv, cp - T * V * beta * beta / kappa), state)
        dT = 1e-2
        set_tpx(ph, T + dT, P, X); Vp = ph.volume_mole
        set_tpx(ph, T - dT, P, X); Vm = ph.volume_mole
        trk["beta"].add(rel(beta * V, (Vp - Vm) / (2 * dT), 1e-12), state)
        dP = max(1e-4 * P, 100.0)
        set_tpx(ph, T, P + dP, X); Vp = ph.volume_mole; mup = ph.chemical_potentials
        set_tpx(ph, T, P - dP, X); Vm = ph.volume_mole; mum = ph.chemical_potentials
        trk["kappa"].add(rel(-kappa * V, (Vp - Vm) / (2 * dP), 1e-18), state)
        set_tpx(ph, T, P, X)
        if check_mu:
            vbar = ph.partial_molar_volumes
            for k in range(2):
                trk["mu"].add(rel((mup[k] - mum[k]) / (2 * dP), vbar[k]), {**state, "k": k})
        cv_iso, resid = isochoric_cv(ph, T, P, X)
        iso_resid = max(iso_resid, resid)
        trk["iso"].add(rel(cv, cv_iso), {**state, "cv_native": cv, "cv_isochoric": cv_iso,
                                          "cp_native": cp})
    res = {f"{name}_cv_eq_test_expression": trk["taut"].out(1e-10),
           f"{name}_beta_times_V_eq_dVdT": trk["beta"].out(1e-5),
           f"{name}_neg_kappa_times_V_eq_dVdP": trk["kappa"].out(1e-5),
           f"{name}_cv_eq_isochoric_T_dsdT": trk["iso"].out(1e-4),
           f"{name}_isochore_max_rel_volume_residual": iso_resid,
           f"{name}_kappa_zero_states": kappa_zero}
    if check_mu:
        res[f"{name}_dmu_dP_eq_vbar"] = trk["mu"].out(1e-5)
    return res


def part_b_c():
    out = {}
    try:
        ph = ct.ThermoPhase(yaml=IAPWS_PHASES_YAML, name="isg-water-ethanol")
        out.update(mayer_block(ph, "B", (285.0, 360.0), (1.0, 200.0), N_POINTS, False))
    except Exception as exc:
        out["B_status"] = f"SKIP: {exc!r}"
    try:
        ph = ct.ThermoPhase(yaml=IAPWS_PHASES_YAML, name="margules-water-ethanol")
        out.update(mayer_block(ph, "C", (285.0, 360.0), (1.0, 200.0), N_POINTS, True))
    except Exception as exc:
        out["C_status"] = f"SKIP: {exc!r}"
    return out


# ---------------------------------------------------------------- D. PDSS_SSVol
def part_d():
    out = {}
    try:
        ph = ct.ThermoPhase(yaml=SSVOL_YAML, name="ssvol-test")
    except Exception as exc:
        return {"D_status": f"SKIP: {exc!r}"}
    c = SSVOL_COEFFS
    p0 = ONE_ATM  # constant-cp reference pressure
    cp0R = SSVOL_CP0 / GAS_CONSTANT
    t1, t2 = MaxTracker(), MaxTracker()
    old_dev = 0.0
    for _ in range(N_POINTS):
        T = RNG.uniform(300.0, 900.0)
        P = RNG.uniform(2.0, 500.0) * ONE_ATM
        ph.TP = T, P
        cpR = ph.standard_cp_R[0]
        d2 = 2 * c[2] + 6 * T * c[3]
        lean_new = cp0R - T * (P - p0) * d2 / GAS_CONSTANT
        lean_old = cp0R - T * (P - p0) * d2
        t1.add(rel(cpR, lean_new), {"T": T, "P": P})
        dT = 1e-2
        ph.TP = T + dT, P; sp = ph.standard_entropies_R[0]
        ph.TP = T - dT, P; sm = ph.standard_entropies_R[0]
        t2.add(rel(cpR, T * (sp - sm) / (2 * dT)), {"T": T, "P": P})
        old_dev = max(old_dev, rel(cpR, lean_old))
    out["D1_cp_formula_postPR"] = t1.out(1e-12)
    out["D2_cp_eq_T_dsdT"] = t2.out(1e-6)
    out["D3_prePR_formula_max_rel_deviation"] = old_dev
    T, P = 600.0, 300.0 * ONE_ATM
    ph.TP = T, P
    d2 = 2 * c[2] + 6 * T * c[3]
    out["D3_prePR_counterexample"] = {
        "input": {"T": T, "P": P, "coeffs": c, "cp0_R": cp0R},
        "native_standard_cp_R": float(ph.standard_cp_R[0]),
        "postPR_formula": cp0R - T * (P - p0) * d2 / GAS_CONSTANT,
        "prePR_formula": cp0R - T * (P - p0) * d2}
    return out


def main():
    for attr in ("isothermal_compressibility", "thermal_expansion_coeff", "partial_molar_volumes"):
        if not hasattr(ct.ThermoPhase, attr):
            print(json.dumps({**RESULT, "status": "SKIP", "reason": f"missing ThermoPhase.{attr}"}))
            return
    RESULT["status"] = "RUN"
    for fn in (part_a, part_b_c, part_d):
        try:
            RESULT.update(fn())
        except Exception as exc:
            RESULT[f"{fn.__name__}_error"] = repr(exc)
    print(json.dumps(RESULT, indent=1, default=float))


if __name__ == "__main__":
    main()
