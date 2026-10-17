# spotcheck.py: unified ΔQ1 spot check for Cantera#2107 (merge commit d58a2326; env/ctb is a Python package built from that commit's source).
# Usage: python spotcheck.py --track {1,2,3} --n N --out PATH
# Track 1 = track1.md (A Mayer/β/κ finite differences and isochoric-path c_v, B chemical potential-activity, C partial molar entropy/entropy sum/c_p/Maxwell);
# Track 2 = track2.md (I1 ∂μ_k/∂P=v̄_k, I2 ∂c_p/∂P=−T∂²v/∂T² and the explicit SSVol formula, I3 Mayer/c_v≤c_p/activity/Gibbs–Duhem);
# Track 3 = audit/specs/track3_guided/Cantera_cantera_2107/spec.md §5 (I1–I6). The Track 3 plan asks for a C++ harness; this script uses only the Python interface,
# so interfaces that exist only in C++ (the Margules derivative methods, WaterProps::ADebye, providePDSS, separate setTemperature/setPressure calls,
# getMolalityActivityCoefficients) are not implemented, see "unimplemented items" in impl_log.md.
# Phase files are taken from the merge-commit source tree test/data (default ../../env/cantera-src/test/data, overridable via the environment variable CT_TEST_DATA).
import argparse
import io
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harness import U, choice, run  # noqa: E402

import cantera as ct  # noqa: E402

PR = "Cantera/cantera#2107"
PR_IDX = 1
R = ct.gas_constant
ATM = ct.one_atm
DATA = os.environ.get("CT_TEST_DATA") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "env", "cantera-src", "test", "data"))
ct.add_data_directory(DATA)

# ---------------- Phase registry ----------------
DH_STATE0 = {"Na+": 9.3549, "Cl-": 9.3549, "H+": 1.05e-08, "OH-": 1.3765e-06, "NaCl(aq)": 0.98492,
             "NaOH(aq)": 3.8836e-06, "NaH3SiO4(aq)": 6.8798e-05, "SiO2(aq)": 3.0179e-05, "H3SiO4-": 1.0231e-06}
HMW_STATE0 = {"Na+": 9.4, "Cl-": 9.4, "H+": 1.05e-05, "OH-": 1.0e-05}
DH_PHASES = ["debye-huckel-dilute", "debye-huckel-dilute-IAPWS", "debye-huckel-B-dot-ak", "debye-huckel-B-dot-ak-IAPWS",
             "debye-huckel-B-dot-a", "debye-huckel-B-dot-a-IAPWS", "debye-huckel-pitzer-beta_ij",
             "debye-huckel-pitzer-beta_ij-IAPWS", "debye-huckel-beta_ij", "debye-huckel-beta_ij-IAPWS"]
REG = {}
for _ph in DH_PHASES:
    REG["dh:" + _ph] = dict(file="debye-huckel-all.yaml", phase=_ph, molal=True, comp=DH_STATE0, kind="dh")
REG.update({
    "tm:debye-huckel-B-dot-ak": dict(file="thermo-models.yaml", phase="debye-huckel-B-dot-ak", molal=True, comp=DH_STATE0, kind="dh"),
    "tm:debye-huckel-beta_ij": dict(file="thermo-models.yaml", phase="debye-huckel-beta_ij", molal=True, comp=DH_STATE0, kind="dh"),
    "HMW_NaCl": dict(file="HMW_NaCl.yaml", phase="", molal=True, comp=HMW_STATE0, kind="hmw"),
    "tm:HMW-NaCl-electrolyte": dict(file="thermo-models.yaml", phase="HMW-NaCl-electrolyte", molal=True, comp=HMW_STATE0, kind="hmw"),
    "tm:HMW-NaCl-HKFT": dict(file="thermo-models.yaml", phase="HMW-NaCl-HKFT", molal=True, comp=HMW_STATE0, kind="hmw_hkft"),
    "tm:ideal-molal-aqueous": dict(file="thermo-models.yaml", phase="ideal-molal-aqueous", molal=True,
                                   comp={"CH4(aq)": 0.01, "H2S(aq)": 0.03, "CO2(aq)": 0.1}, kind="molal"),
    "tm:ethanol-water-rough": dict(file="thermo-models.yaml", phase="ethanol-water-rough", molal=False,
                                   comp={"ethanol": 0.3, "water": 0.7}, kind="margules_ev"),
    "pdss_hkft": dict(file="pdss_hkft.yaml", phase="", molal=False,
                      comp={"H2O(L)": 0.9, "Na+": 0.05, "Cl-": 0.05, "H+": 1e-7, "OH-": 1e-7}, kind="hkft_vpss"),
    "LiKCl_liquid": dict(file="LiKCl_liquid.yaml", phase="", molal=False, comp={"KCl(L)": 0.2, "LiCl(L)": 0.8}, kind="molten"),
    "tm:molten-salt-Margules": dict(file="thermo-models.yaml", phase="molten-salt-Margules", molal=False,
                                    comp={"KCl(l)": 0.3, "LiCl(l)": 0.7}, kind="molten"),
    "Li_Liquid": dict(file="Li_Liquid.yaml", phase="Li(L)", molal=False, comp={"Li(L)": 0.4, "Li(L)2": 0.6}, kind="ssvol"),
    "VpssSolidSolutionExample": dict(file="IdealSolidSolnPhaseExample.yaml", phase="VpssSolidSolutionExample", molal=False,
                                     comp={"C2H2-graph": 0.2, "C-graph": 0.5, "H2-solute": 0.3}, kind="vpss_solid"),
    "tm:IdealSolnGas-liquid": dict(file="thermo-models.yaml", phase="IdealSolnGas-liquid", molal=False, comp={"Li(l)": 1.0},
                                   kind="vpss_solid"),
    "tm:Redlich-Kister-LiC6": dict(file="thermo-models.yaml", phase="Redlich-Kister-LiC6", molal=False,
                                   comp={"Li(C6)": 0.85, "V(C6)": 0.15}, kind="rk"),
    "tm:Redlich-Kister-complex": dict(file="thermo-models.yaml", phase="Redlich-Kister-complex", molal=False,
                                      comp={"KCl(l)": 0.85, "LiCl(l)": 0.15}, kind="rk"),
    "anode": dict(file="BinarySolutionTabulatedThermo.yaml", phase="anode", molal=False,
                  comp={"Li[anode]": 0.3, "V[anode]": 0.7}, kind="tab"),
})

_PH = {}


def _yaml():
    from ruamel.yaml import YAML
    y = YAML()
    y.preserve_quotes = True
    return y


def _load_doc(fname):
    with open(os.path.join(DATA, fname)) as f:
        return _yaml().load(f)


def _dump(doc):
    s = io.StringIO()
    _yaml().dump(doc, s)
    return s.getvalue()


def _phase_entry(doc, name):
    for ph in doc["phases"]:
        if ph["name"] == name:
            return ph
    raise KeyError(name)


def variant(key):
    """Synthesized/rewritten phase (ct.ThermoPhase(yaml=...) after the YAML rewrite)."""
    if key.startswith("molal_cutoff:"):
        doc = _load_doc("thermo-models.yaml")
        ph = _phase_entry(doc, "ideal-molal-aqueous")
        ph["cutoff"]["model"] = key.split(":")[1]
        return ct.ThermoPhase(yaml=_dump(doc), name="ideal-molal-aqueous")
    if key.startswith("molal_stdconc:"):
        doc = _load_doc("thermo-models.yaml")
        ph = _phase_entry(doc, "ideal-molal-aqueous")
        ph["standard-concentration-basis"] = key.split(":")[1]
        return ct.ThermoPhase(yaml=_dump(doc), name="ideal-molal-aqueous")
    if key.startswith("tab_simplified") or key == "anode:swapped":
        doc = _load_doc("BinarySolutionTabulatedThermo.yaml")
        ph = _phase_entry(doc, "anode")
        if key.startswith("tab_simplified"):
            ph["standard-concentration-basis"] = "unity"
            for sp in doc["species"]:
                sp["thermo"] = {"model": "constant-cp", "T0": 298.15, "h0": 0.0, "s0": 0.0, "cp0": 0.0}
        if key.endswith("swapped"):
            ph["species"] = list(reversed(list(ph["species"])))
        return ct.ThermoPhase(yaml=_dump(doc), name="anode")
    raise KeyError(key)


def get_phase(key):
    if key not in _PH:
        if key in REG:
            r = REG[key]
            _PH[key] = ct.ThermoPhase(r["file"], r["phase"]) if r["phase"] else ct.ThermoPhase(r["file"])
        else:
            _PH[key] = variant(key)
    return _PH[key]


def margules_synthetic(coef):
    doc = _load_doc("thermo-models.yaml")
    ph = _phase_entry(doc, "ethanol-water-rough")
    it = ph["interactions"][0]
    it["excess-enthalpy"] = ["%r J/mol" % coef[0], "%r J/mol" % coef[1]]
    it["excess-entropy"] = ["%r J/mol/K" % coef[2], "%r J/mol/K" % coef[3]]
    it["excess-volume-enthalpy"] = ["%r cm^3/mol" % coef[4], "%r cm^3/mol" % coef[5]]
    it["excess-volume-entropy"] = ["%r cm^3/mol/K" % coef[6], "%r cm^3/mol/K" % coef[7]]
    return ct.ThermoPhase(yaml=_dump(doc), name="ethanol-water-rough")


def is_water_like(sol):
    """Phase containing IAPWS water or HKFT species (decides whether the tolerance relaxation for "phases containing IAPWS" applies in each track)."""
    for k in range(sol.n_species):
        eos = sol.species(k).input_data.get("equation-of-state", {})
        if isinstance(eos, dict) and eos.get("model") in ("liquid-water-IAPWS95", "HKFT"):
            return True
        if isinstance(eos, list) and any(e.get("model") in ("liquid-water-IAPWS95", "HKFT") for e in eos):
            return True
    return False


def comp_X(sol, comp, molal):
    """Composition: mole-fraction phases are used directly; molality phases are converted from m_k to X (X_o=1/(1+M̃_oΣm), X_k=m_k M̃_o X_o, with M̃_o in kg/mol)."""
    names = sol.species_names
    if not molal:
        X = np.array([float(comp.get(n, 0.0)) for n in names])
        return X / X.sum()
    Mo = sol.molecular_weights[0] / 1000.0
    m = np.array([0.0] + [float(comp.get(n, 0.0)) for n in names[1:]])
    Xo = 1.0 / (1.0 + Mo * m.sum())
    X = m * Mo * Xo
    X[0] = Xo
    return X


def molalities(sol):
    X = sol.X
    Mo = sol.molecular_weights[0] / 1000.0
    return X / (X[0] * Mo)


def setTPX(sol, T, P, X):
    sol.TPX = T, P, X


# ---------------- Common quantities ----------------
def props(sol):
    return dict(h=sol.enthalpy_mole, s=sol.entropy_mole, u=sol.int_energy_mole, g=sol.gibbs_mole, cp=sol.cp_mole,
                v=sol.volume_mole, mu=np.array(sol.chemical_potentials), vb=np.array(sol.partial_molar_volumes),
                sb=np.array(sol.partial_molar_entropies), hb=np.array(sol.partial_molar_enthalpies),
                cpb=np.array(sol.partial_molar_cp))


def at(sol, T, P, X, f):
    setTPX(sol, T, P, X)
    return f(sol)


def cd_T(sol, T, P, X, f, h):
    return (at(sol, T + h, P, X, f) - at(sol, T - h, P, X, f)) / (2 * h)


def cd_P(sol, T, P, X, f, d):
    return (at(sol, T, P + d, X, f) - at(sol, T, P - d, X, f)) / (2 * d)


def isochoric_newton(sol, T, P, X, v0, beta, kappa, h, steps=3):
    """Track 1 A: (T,p)→(T+h, p+(β/κ)h), then use Newton to pull v back to v0 (2–3 steps). Returns (p', s', u', residual)."""
    p = P + beta / kappa * h
    for _ in range(steps):
        setTPX(sol, T + h, p, X)
        v = sol.volume_mole
        dvdp = -v * sol.isothermal_compressibility
        if dvdp == 0:
            break
        p -= (v - v0) / dvdp
    setTPX(sol, T + h, p, X)
    return p, sol.entropy_mole, sol.int_energy_mole, abs(sol.volume_mole - v0) / abs(v0)


# ============================ Track 1 ============================
T1_PHASES = [k for k in REG if k not in ("tm:Redlich-Kister-LiC6", "tm:Redlich-Kister-complex", "anode")]


def t1_domain(kind):
    if kind in ("dh", "hmw", "molal", "margules_ev"):
        return (273.15, 373.15), (1 * ATM, 100 * ATM)
    if kind in ("hkft_vpss", "hmw_hkft"):
        return (273.15, 573.15), (1 * ATM, 100 * ATM)
    if kind == "molten":
        return (700.0, 1400.0), None
    return (300.0, 1000.0), None   # ssvol / vpss_solid: VPSS interval for solid solutions


def t1_sampler(key):
    r = REG[key]

    def s(rng):
        (t0, t1), prange = t1_domain(r["kind"])
        p = dict(key=key, T=U(rng, t0, t1), P=U(rng, *prange) if prange else ATM)
        if r["kind"] in ("dh", "hmw", "hmw_hkft"):
            p["I"] = U(rng, 1e-4, 6.1 if r["kind"] != "dh" else 6.0)
        return p
    return s


def _comp_for(key, p):
    r = REG[key]
    comp = dict(r["comp"])
    if "I" in p:
        comp["Na+"] = comp["Cl-"] = p["I"]
    return comp


def t1_eval(p, P):
    key = p["key"]
    r = REG[key]
    sol = get_phase(key)
    T, Pr = p["T"], p["P"]
    X = comp_X(sol, _comp_for(key, p), r["molal"])
    wl = is_water_like(sol)
    rt = 4e-6 if wl else 1e-7
    setTPX(sol, T, Pr, X)
    cv, cp, beta, kap, v = sol.cv_mole, sol.cp_mole, sol.thermal_expansion_coeff, sol.isothermal_compressibility, sol.volume_mole
    s0 = sol.entropy_mole
    h = 1e-4 * T
    dl = 1e-4 * Pr
    vf = lambda q: q.volume_mole  # noqa: E731
    beta_fd = cd_T(sol, T, Pr, X, vf, h) / v
    kap_fd = -cd_P(sol, T, Pr, X, vf, dl) / v
    P.close("A.beta=FD(v)/v", beta, beta_fd, rt, 0.0)
    P.close("A.kappaT=-FD(v)/v", kap, kap_fd, rt, 0.0)
    if kap > 0:
        P.close("A.Mayer cv=cp-Tvb^2/k", cv, cp - T * v * beta ** 2 / kap, rt, 0.0)
        sp_, _u, _r1 = isochoric_newton(sol, T, Pr, X, v, beta, kap, h)[1:]
        sm_, _u2, _r2 = isochoric_newton(sol, T, Pr, X, v, beta, kap, -h)[1:]
        P.record("A.isochoric_residual", plus=_r1, minus=_r2)
        P.close("A.cv=T*ds/dT|v(isochoric Newton)", cv, T * (sp_ - sm_) / (2 * h), 1e-5, 0.0)
    else:
        P.close("A.kappa=0 => cv=cp (library convention)", cv, cp)
    setTPX(sol, T, Pr, X)
    # B: μ_k = μ_k° + RT ln(γ_k X_k)
    mu = np.array(sol.chemical_potentials)
    g0 = np.array(sol.standard_gibbs_RT) * R * T
    gam = np.array(sol.activity_coefficients)
    for k in range(sol.n_species):
        if X[k] > 0:
            P.close("B.mu=mu0+RTln(gX)[%s]" % sol.species_name(k), mu[k], g0[k] + R * T * math.log(gam[k] * X[k]), 1e-10, 0.0)
    if r["molal"]:
        act = np.array(sol.activities)  # molality-basis activity a_k=γ^□ m_k/m°(MolalityVPSSTP)
        for k in range(1, sol.n_species):
            if X[k] > 0:
                P.close("B.gamma_k*X_k=gamma^m*m_k[%s]" % sol.species_name(k), gam[k] * X[k], act[k])
    # C
    sb = np.array(sol.partial_molar_entropies)
    dmudT = cd_T(sol, T, Pr, X, lambda q: np.array(q.chemical_potentials), h)
    for k in range(sol.n_species):
        P.close("C.dmu/dT=-sbar[%s]" % sol.species_name(k), dmudT[k], -sb[k], 1e-6, 0.0)
    P.close("C.sum Xk sk = s", float(np.dot(X, sb)), s0, 1e-12, 0.0)
    P.close("C.cp=dh/dT", cp, cd_T(sol, T, Pr, X, lambda q: q.enthalpy_mole, h))
    P.close("C.ds/dP=-dv/dT", cd_P(sol, T, Pr, X, lambda q: q.entropy_mole, dl), -cd_T(sol, T, Pr, X, vf, h))


def track1():
    return [(k, t1_sampler(k), t1_eval) for k in T1_PHASES], []


# ============================ Track 2 ============================
T2_PHASES = ["LiKCl_liquid", "tm:ethanol-water-rough", "Li_Liquid"] + ["dh:" + x for x in DH_PHASES] + [
    "HMW_NaCl", "tm:ideal-molal-aqueous", "VpssSolidSolutionExample", "pdss_hkft"]
T2_I1 = {"LiKCl_liquid", "tm:ethanol-water-rough", "Li_Liquid", "HMW_NaCl", "tm:ideal-molal-aqueous",
         "VpssSolidSolutionExample", "pdss_hkft"} | {"dh:" + x for x in DH_PHASES}


def t2_sampler(key):
    r = REG[key]
    k = r["kind"]

    def s(rng):
        if key == "LiKCl_liquid":
            T = U(rng, 900, 1400)
        elif key == "tm:ethanol-water-rough":
            T = U(rng, 298, 308)
        elif key in ("Li_Liquid", "VpssSolidSolutionExample"):
            T = U(rng, 300, 500)
        elif key == "tm:ideal-molal-aqueous":
            T = U(rng, 300, 340)
        else:
            T = U(rng, 300, 360)
        p = dict(key=key, T=T, P=U(rng, 1 * ATM, 100 * ATM))
        if k in ("dh", "hmw"):
            p["m"] = U(rng, 0.01, 9.4)
        if key == "tm:ethanol-water-rough":
            p["xe"] = U(rng, 0.3, 0.5)
        p["gd_seed"] = int(rng.integers(1 << 30))
        return p
    return s


def ssvol_poly(sol, k):
    """v°(T) polynomial coefficients (SI) of an SSVol species; returned only for molar-volume-temperature-polynomial."""
    eos = sol.species(k).input_data.get("equation-of-state", {})
    if isinstance(eos, dict) and eos.get("model") == "molar-volume-temperature-polynomial":
        return [float(c) for c in eos["data"]]
    return None


def ssvol_density_poly(sol, k):
    eos = sol.species(k).input_data.get("equation-of-state", {})
    if isinstance(eos, dict) and eos.get("model") == "density-temperature-polynomial":
        return [float(c) for c in eos["data"]]
    return None


def t2_eval(p, P):
    key = p["key"]
    r = REG[key]
    sol = get_phase(key)
    T, Pr = p["T"], p["P"]
    comp = dict(r["comp"])
    if "m" in p:
        comp["Na+"] = comp["Cl-"] = p["m"]
    if "xe" in p:
        comp = {"ethanol": p["xe"], "water": 1 - p["xe"]}
    X = comp_X(sol, comp, r["molal"])
    wl = is_water_like(sol)
    setTPX(sol, T, Pr, X)
    vb = np.array(sol.partial_molar_volumes)
    cp, cv, v, beta, kap = sol.cp_mole, sol.cv_mole, sol.volume_mole, sol.thermal_expansion_coeff, sol.isothermal_compressibility
    dP = max(1e-3 * Pr, 1e3)
    if key in T2_I1:
        dmu = cd_P(sol, T, Pr, X, lambda q: np.array(q.chemical_potentials), dP)
        atol = 1e-12 * float(np.max(np.abs(vb)))
        for k in range(sol.n_species):
            P.close("I1.dmu/dP=vbar[%s]" % sol.species_name(k), dmu[k], vb[k], 1e-5 if wl else 1e-6, atol)
    # I2
    hT = 0.5
    vf = lambda q: q.volume_mole  # noqa: E731
    v_tt = (at(sol, T + hT, Pr, X, vf) - 2 * at(sol, T, Pr, X, vf) + at(sol, T - hT, Pr, X, vf)) / hT ** 2
    dcp = cd_P(sol, T, Pr, X, lambda q: q.cp_mole, dP)
    P.close("I2.dcp/dP=-T d2v/dT2", dcp, -T * v_tt, 1e-4, 0.0)
    for k in range(sol.n_species):
        c = ssvol_poly(sol, k)
        if c is None:
            continue
        vtt = 2 * c[2] + (6 * c[3] * T if len(c) > 3 else 0.0)
        setTPX(sol, T, Pr, X)
        cpk = sol.standard_cp_R[k] * R
        setTPX(sol, T, ATM, X)
        cpk_ref = sol.standard_cp_R[k] * R
        P.close("I2.SSVol cp0(T,P)-cp0(T,Pref)=-T*v''*(P-Pref)[%s]" % sol.species_name(k), cpk - cpk_ref,
                -T * vtt * (Pr - ATM), 1e-4, 0.0)
    # I3
    setTPX(sol, T, Pr, X)
    if kap > 0:
        P.close("I3.Mayer cv=cp-Tvb^2/k", cv, cp - T * v * beta ** 2 / kap, 1e-10, 0.0)
    else:
        P.close("I3.kappa=0 => cv=cp", cv, cp, 1e-10, 0.0)
    P.cond("I3.cv<=cp", cv <= cp * (1 + 1e-12) + 1e-12, obs=cv, exp=cp)
    if r["molal"]:
        act = np.array(sol.activities)
        gam = np.array(sol.activity_coefficients)
        for k in range(sol.n_species):
            P.close("I3.activities=gamma*X[%s]" % sol.species_name(k), act[k], gam[k] * X[k], 1e-12, 0.0)
    rng = np.random.default_rng(p["gd_seed"])
    rr = rng.standard_normal(sol.n_species)
    d = 1e-6 * X * (rr - np.dot(X, rr))
    setTPX(sol, T, Pr, X - d)
    mu0 = np.array(sol.chemical_potentials)
    setTPX(sol, T, Pr, X + d)
    mu1 = np.array(sol.chemical_potentials)
    dmu = (mu1 - mu0) / 2  # central difference (implementation judgment: cancels the second-order term, see impl_log)
    rel = abs(float(np.dot(X, dmu))) / float(np.max(np.abs(dmu)))
    P.cond("I3.Gibbs-Duhem |sum X dmu|/max|dmu|<1e-6", rel < 1e-6, obs=rel, exp=0.0, resid=rel)


def track2():
    return [(k, t2_sampler(k), t2_eval) for k in T2_PHASES], []


# ============================ Track 3 ============================
HT = (1e-4, 5e-5, 2.5e-5)
HP = (1e-3, 5e-4, 2.5e-4)
FLOOR = dict(cp=1e-3, v=1e-10, maxwell=1e-9, kappa=5e-14, beta=1e-9)


def alg(P, name, a, b):
    """Track 3 algebraic identity: |a−b| ≤ 1e-5 + 1e-12·max(|a|,|b|) (energy units)."""
    a, b = float(a), float(b)
    if not (math.isfinite(a) and math.isfinite(b)):
        P.fail(name, "nonfinite %r %r" % (a, b))
        return
    d = abs(a - b)
    P.checks.append(dict(name=name, ok=bool(d <= 1e-5 + 1e-12 * max(abs(a), abs(b))), obs=a, exp=b, resid=d, tol=None))


def alg_v(P, name, a, b):
    """Volume-type algebraic identity: |a−b| ≤ 1e-10 + 1e-12·max (implementation judgment: use a volume floor instead of the energy floor)."""
    a, b = float(a), float(b)
    d = abs(a - b)
    P.checks.append(dict(name=name, ok=bool(d <= FLOOR["v"] + 1e-12 * max(abs(a), abs(b))), obs=a, exp=b, resid=d, tol=None))


def deriv3(P, name, analytic, fd_fun, steps, rtol, floor):
    """Differences at three step sizes; first check whether the last two have converged (not converged = undetermined, recorded as a fail), then compare with the analytic value."""
    est = [fd_fun(s) for s in steps]
    e2, e3 = est[1], est[2]
    tol_c = rtol * max(abs(e3), abs(analytic)) + floor
    if not (math.isfinite(e3) and math.isfinite(analytic)):
        P.fail(name, "nonfinite")
        return
    if abs(e2 - e3) > tol_c:
        P.fail(name, "stencil not converged e2=%r e3=%r" % (e2, e3))
        return
    d = abs(analytic - e3)
    P.checks.append(dict(name=name, ok=bool(d <= rtol * max(abs(e3), abs(analytic)) + floor), obs=analytic, exp=e3, resid=d, tol=tol_c))


T3_LIST = ["dh:" + x for x in DH_PHASES] + ["HMW_NaCl", "molal_cutoff:none", "molal_cutoff:poly", "molal_cutoff:polyexp",
                                            "molal_stdconc:unity", "molal_stdconc:species-molar-volume", "pdss_hkft",
                                            "VpssSolidSolutionExample", "Li_Liquid", "LiKCl_liquid", "tm:ethanol-water-rough",
                                            "margules_synthetic", "tm:Redlich-Kister-LiC6", "tm:Redlich-Kister-complex",
                                            "anode", "tab_simplified"]


def t3_kind(key):
    if key.startswith("molal_"):
        return "molal"
    if key == "margules_synthetic":
        return "margules_ev"
    if key == "tab_simplified":
        return "tab"
    return REG[key]["kind"]


def t3_sampler(key):
    kind = t3_kind(key)

    def s(rng):
        p = dict(key=key)
        if kind in ("dh", "hmw", "molal", "hkft_vpss"):
            p["T"] = U(rng, 280, 360)
            p["P"] = U(rng, 0.101325e6, 10e6)
            p["m"] = U(rng, 1e-6, 1e-2)
            if kind == "molal":
                p["Xo"] = U(rng, 1e-6, 1 - 1e-6)
        elif kind == "margules_ev":
            p["T"] = U(rng, 280, 330)
            p["P"] = U(rng, 1 * ATM, 20 * ATM)
            p["x1"] = U(rng, 1e-6, 1 - 1e-6)
        elif kind == "ssvol":
            p["T"] = U(rng, 300, 500)
            p["P"] = U(rng, 1 * ATM, 20 * ATM)
            p["x1"] = U(rng, 1e-6, 1 - 1e-6)
        elif kind == "molten":
            p["T"] = U(rng, 900, 1400)
            p["P"] = U(rng, 1 * ATM, 20 * ATM)
            p["x1"] = U(rng, 1e-6, 1 - 1e-6)
        elif kind == "rk":  # the plan gives no interval: take the range of the consistency-cases state points
            p["T"] = U(rng, 300, 440) if key.endswith("LiC6") else U(rng, 900, 950)
            p["P"] = U(rng, 0.1 * ATM, 10 * ATM)
            p["x1"] = U(rng, 1e-6, 1 - 1e-6)
        elif kind == "vpss_solid":
            p["T"] = U(rng, 300, 500)
            p["P"] = U(rng, 1 * ATM, 10 * ATM)
            w = rng.uniform(size=3)
            p["Xv"] = list(w / w.sum())
        elif kind == "tab":
            p["T"] = U(rng, 280, 350)
            p["P"] = U(rng, 1 * ATM, 10 * ATM)
            p["x"] = U(rng, 0, 1)
            p["order"] = choice(rng, [0, 1])
        if key == "margules_synthetic":
            base = [-2449, -1613, -14.45, -11.64, -12.707, 0.476, -0.030, 0.0]
            scale = [3000, 3000, 20, 20, 15, 15, 0.05, 0.05]
            p["coef"] = [0.0 if rng.uniform() < 0.25 else U(rng, -scale[i], scale[i]) for i in range(8)]
            p["coef_base"] = base
        p["gd_seed"] = int(rng.integers(1 << 30))
        return p
    return s


def t3_phase(p):
    if p["key"] == "margules_synthetic":
        return margules_synthetic(p["coef"])
    if t3_kind(p["key"]) == "tab" and p.get("order", 0) == 1:  # two species orderings (plan: both species orderings)
        return get_phase("anode:swapped" if p["key"] == "anode" else "tab_simplified:swapped")
    return get_phase(p["key"])


def t3_X(sol, p):
    key = p["key"]
    kind = t3_kind(key)
    if kind in ("dh", "hmw"):
        comp = dict(REG[key]["comp"])
        comp["Na+"] = comp["Cl-"] = p["m"]
        return comp_X(sol, comp, True)
    if kind == "hkft_vpss":
        Mo = 0.018015
        mH = 1e-7 / 0.9 / Mo
        Xo = 1.0 / (1.0 + Mo * (2 * p["m"] + 2 * mH))
        comp = {"H2O(L)": Xo, "Na+": p["m"] * Mo * Xo, "Cl-": p["m"] * Mo * Xo, "H+": mH * Mo * Xo, "OH-": mH * Mo * Xo}
        return comp_X(sol, comp, False)
    if kind == "molal":
        names = sol.species_names  # solutes in the proportions of state point 0: CO2 0.1, H2S 0.03, CH4 0.01
        order = {"CO2(aq)": 0.1, "H2S(aq)": 0.03, "CH4(aq)": 0.01}
        X = np.array([0.0] + [order[n] for n in names[1:]])
        X = X / X.sum() * (1 - p["Xo"])
        X[0] = p["Xo"]
        return X
    if kind == "vpss_solid":
        return np.array(p["Xv"])
    x1 = p["x1"]
    return np.array([x1, 1 - x1] + [0.0] * (sol.n_species - 2)) if sol.n_species > 2 else np.array([x1, 1 - x1])


def t3_I1(sol, P, T, Pr, X, wl):
    rt = 1e-5 if wl else 1e-6
    setTPX(sol, T, Pr, X)
    q = props(sol)
    alg(P, "I1.g=h-Ts", q["g"], q["h"] - T * q["s"])
    alg(P, "I1.h=u+Pv", q["h"], q["u"] + Pr * q["v"])
    for k in range(sol.n_species):
        alg(P, "I1.mu=hbar-T*sbar[%s]" % sol.species_name(k), q["mu"][k], q["hb"][k] - T * q["sb"][k])
    alg_v(P, "I1.v=sum X vbar", q["v"], float(np.dot(X, q["vb"])))
    alg(P, "I1.h=sum X hbar", q["h"], float(np.dot(X, q["hb"])))
    alg(P, "I1.s=sum X sbar", q["s"], float(np.dot(X, q["sb"])))
    alg(P, "I1.g=sum X mu", q["g"], float(np.dot(X, q["mu"])))
    alg(P, "I1.cp=sum X cpbar", q["cp"], float(np.dot(X, q["cpb"])))
    fT = lambda f: (lambda s: cd_T(sol, T, Pr, X, f, s * T))  # noqa: E731
    fP = lambda f: (lambda s: cd_P(sol, T, Pr, X, f, s * Pr))  # noqa: E731
    deriv3(P, "I1.cp=dh/dT", q["cp"], fT(lambda z: z.enthalpy_mole), HT, rt, FLOOR["cp"])
    deriv3(P, "I1.cp=T ds/dT", q["cp"], lambda s: T * fT(lambda z: z.entropy_mole)(s), HT, rt, FLOOR["cp"])
    vT = [fT(lambda z: z.volume_mole)(s) for s in HT]
    deriv3(P, "I1.ds/dP=-dv/dT", -vT[2], fP(lambda z: z.entropy_mole), HP, rt, FLOOR["maxwell"])
    for k in range(sol.n_species):
        if X[k] <= 0:
            continue
        deriv3(P, "I1.vbar=dmu/dP[%s]" % sol.species_name(k), q["vb"][k],
               fP(lambda z, k=k: z.chemical_potentials[k]), HP, rt, FLOOR["v"])
    setTPX(sol, T, Pr, X)
    return q


def t3_ssvol(sol, P, T, Pr, X):
    """SSVol analytic forms (two polynomials): g°−g_r°=ΔP V°, s°−s_r°=−ΔP V_T°, h°−h_r°=ΔP(V°−T V_T°), cp°−cp_r°=−T ΔP V_TT°."""
    setTPX(sol, T, Pr, X)
    g, s, h, c = (np.array(sol.standard_gibbs_RT) * R * T, np.array(sol.standard_entropies_R) * R,
                  np.array(sol.standard_enthalpies_RT) * R * T, np.array(sol.standard_cp_R) * R)
    setTPX(sol, T, ATM, X)
    gr, sr, hr, cr = (np.array(sol.standard_gibbs_RT) * R * T, np.array(sol.standard_entropies_R) * R,
                      np.array(sol.standard_enthalpies_RT) * R * T, np.array(sol.standard_cp_R) * R)
    dP = Pr - ATM
    for k in range(sol.n_species):
        cv_ = ssvol_poly(sol, k)
        cd_ = ssvol_density_poly(sol, k)
        if cv_ is not None:
            a = cv_ + [0.0] * (4 - len(cv_))
            V = a[0] + a[1] * T + a[2] * T ** 2 + a[3] * T ** 3
            VT = a[1] + 2 * a[2] * T + 3 * a[3] * T ** 2
            VTT = 2 * a[2] + 6 * a[3] * T
            tag = "volume-poly"
        elif cd_ is not None:
            d = cd_ + [0.0] * (4 - len(cd_))
            M = sol.molecular_weights[k]
            rho = d[0] + d[1] * T + d[2] * T ** 2 + d[3] * T ** 3
            rT = d[1] + 2 * d[2] * T + 3 * d[3] * T ** 2
            rTT = 2 * d[2] + 6 * d[3] * T
            V = M / rho
            VT = -M * rT / rho ** 2
            VTT = M * (2 * rT ** 2 / rho ** 3 - rTT / rho ** 2)
            tag = "density-poly"
        else:
            continue
        nm = "%s,%s" % (sol.species_name(k), tag)
        alg(P, "I1.SSVol g0-gr0=dP V[%s]" % nm, g[k] - gr[k], dP * V)
        alg(P, "I1.SSVol h0-hr0=dP(V-T VT)[%s]" % nm, h[k] - hr[k], dP * (V - T * VT))
        P.close("I1.SSVol s0-sr0=-dP VT[%s]" % nm, s[k] - sr[k], -dP * VT, 1e-6, FLOOR["cp"])
        P.close("I1.SSVol cp0-cpr0=-T dP VTT[%s]" % nm, c[k] - cr[k], -T * dP * VTT, 1e-6, FLOOR["cp"])


def t3_setorder(sol, P, T, Pr, X):
    """Track 3 I1: the results are consistent across three setting orders: T then P, P then T, and a single TP (in Python TP=(T,None)/(None,P) can only approximate setTemperature/setPressure)."""
    out = []
    for order in ("TP", "T_then_P", "P_then_T"):
        setTPX(sol, 0.9 * T, 0.5 * Pr, X)
        if order == "TP":
            sol.TP = T, Pr
        elif order == "T_then_P":
            sol.TP = T, None
            sol.TP = None, Pr
        else:
            sol.TP = None, Pr
            sol.TP = T, None
        out.append((sol.enthalpy_mole, sol.entropy_mole, sol.cp_mole, sol.volume_mole))
    for i in (1, 2):
        for j, nm in enumerate(("h", "s", "cp")):
            alg(P, "I1.set_order_equal[%d,%s]" % (i, nm), out[i][j], out[0][j])
        alg_v(P, "I1.set_order_equal[%d,v]" % i, out[i][3], out[0][3])


def _root_P(sol, T, X, v0, P0, rel=1e-11):
    """Isochoric root: find P such that v(T,P)=v0 (bracketing interval + bisection/secant). Returns P or None (cannot bracket = undetermined)."""
    def f(pp):
        setTPX(sol, T, pp, X)
        return sol.volume_mole - v0
    a, b = P0, P0
    fa = fb = f(P0)
    step = max(1e3, 1e-3 * P0)
    for _ in range(60):
        a, b = max(1.0, a - step), b + step
        fa, fb = f(a), f(b)
        if fa * fb <= 0:
            break
        step *= 2
    else:
        return None
    if fa * fb > 0:
        return None
    for _ in range(200):
        m = 0.5 * (a + b)
        fm = f(m)
        if abs(fm) <= rel * abs(v0):
            return m
        if fa * fm <= 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return 0.5 * (a + b) if abs(f(0.5 * (a + b))) <= rel * abs(v0) else None


def t3_I2(sol, P, T, Pr, X, wl, kind):
    rt = 1e-5 if wl else 1e-6
    setTPX(sol, T, Pr, X)
    snap0 = (sol.enthalpy_mole, sol.entropy_mole, sol.volume_mole, tuple(sol.chemical_potentials))
    v, beta, kap, cp = sol.volume_mole, sol.thermal_expansion_coeff, sol.isothermal_compressibility, sol.cp_mole
    try:
        cv = sol.cv_mole
    except Exception as e:
        if kind == "hkft_vpss":
            P.record("I2.HKFT direct cv_mole unsupported", err="%s: %s" % (type(e).__name__, str(e)[:100]))
        cv = None
    for _ in range(3):  # after repeated derivative calls, the state and other properties must not change
        sol.thermal_expansion_coeff, sol.isothermal_compressibility
    snap1 = (sol.enthalpy_mole, sol.entropy_mole, sol.volume_mole, tuple(sol.chemical_potentials))
    P.cond("I2.derivative_calls_leave_state_unchanged", snap0 == snap1)
    alg_v(P, "I2.v=Mbar/rho", v, sol.mean_molecular_weight / sol.density_mass)
    vf = lambda z: z.volume_mole  # noqa: E731
    vT = [cd_T(sol, T, Pr, X, vf, s * T) for s in HT]
    vP = [cd_P(sol, T, Pr, X, vf, s * Pr) for s in HP]
    deriv3(P, "I2.beta=v_T/v", beta, lambda s: cd_T(sol, T, Pr, X, vf, s * T) / v, HT, rt, FLOOR["beta"])
    deriv3(P, "I2.kappa=-v_P/v", kap, lambda s: -cd_P(sol, T, Pr, X, vf, s * Pr) / v, HP, rt, FLOOR["kappa"])
    setTPX(sol, T, Pr, X)
    if cv is None:
        return
    if vP[2] != 0.0:
        P.close("I2.cp-cv=-T vT^2/vP (FD volume)", cp - cv, -T * vT[2] ** 2 / vP[2], rt, FLOOR["cp"])
        res = {}
        for s in HT[1:]:
            hh = s * T
            vals = []
            for sg in (1, -1):
                pr = _root_P(sol, T + sg * hh, X, v, Pr)
                if pr is None:
                    vals = None
                    break
                setTPX(sol, T + sg * hh, pr, X)
                vals.append((sol.int_energy_mole, sol.entropy_mole))
            if vals is None:
                P.fail("I2.isochoric_root", "cannot bracket / root not converged to 1e-11 relative volume")
                return
            res[s] = ((vals[0][0] - vals[1][0]) / (2 * hh), T * (vals[0][1] - vals[1][1]) / (2 * hh))
        e2, e3 = res[HT[1]], res[HT[2]]
        for j, nm in enumerate(("du/dT|v", "T ds/dT|v")):
            if abs(e2[j] - e3[j]) > rt * abs(e3[j]) + FLOOR["cp"]:
                P.fail("I2.cv=%s" % nm, "stencil not converged")
            else:
                P.close("I2.cv=%s" % nm, cv, e3[j], rt, FLOOR["cp"])
    elif vT[2] != 0.0:
        P.record("I2.vP=0,vT!=0 convention cv=cp", cv=cv, cp=cp, equal=bool(abs(cv - cp) <= 1e-12 * abs(cp)))
    else:
        P.close("I2.rigid cv=cp", cv, cp, 1e-12, 0.0)
    setTPX(sol, T, Pr, X)


def t3_I3(sol, P, T, Pr, X):
    setTPX(sol, T, Pr, X)
    a = np.array(sol.activities)
    gam = np.array(sol.activity_coefficients)
    mu = np.array(sol.chemical_potentials)
    mu0 = np.array(sol.standard_gibbs_RT) * R * T
    for k in range(sol.n_species):
        if X[k] <= 0:
            continue
        P.close("I3.a=X*gamma[%s]" % sol.species_name(k), a[k], X[k] * gam[k], 1e-10, 1e-20)
        if a[k] > 0:
            P.close("I3.(mu-mu0)/RT=ln a[%s]" % sol.species_name(k), (mu[k] - mu0[k]) / (R * T), math.log(a[k]), 0.0, 1e-9)


def margules_excess(sol, T, Pr, X):
    """Independently compute g^E, h^E, s^E, v^E from the phase file's interaction parameters (input_data already converted to SI) (formulas as documented in the MargulesVPSSTP header)."""
    names = sol.species_names
    gE = hE = sE = vE = 0.0
    dP = Pr - ATM
    for it in sol.input_data.get("interactions", []):
        a, b = names.index(it["species"][0]), names.index(it["species"][1])
        xa, xb = X[a], X[b]
        h0, h1 = (list(it.get("excess-enthalpy", [0, 0])) + [0, 0])[:2]
        s0, s1 = (list(it.get("excess-entropy", [0, 0])) + [0, 0])[:2]
        vh0, vh1 = (list(it.get("excess-volume-enthalpy", [0, 0])) + [0, 0])[:2]
        vs0, vs1 = (list(it.get("excess-volume-entropy", [0, 0])) + [0, 0])[:2]
        w = xa * xb
        H = w * (h0 + h1 * xb)
        S = w * (s0 + s1 * xb)
        V = w * ((vh0 - T * vs0) + (vh1 - T * vs1) * xb)
        dVdT = -w * (vs0 + vs1 * xb)
        gE += H - T * S + dP * V
        sE += S - dP * dVdT
        vE += V
    hE = gE + T * sE
    return gE, hE, sE, vE


def t3_I4(sol, P, T, X):
    names = sol.species_names
    for Pr in (ATM, 2 * ATM, 20 * ATM, 10e6):
        setTPX(sol, T, Pr, X)
        g, h, s, v = sol.gibbs_mole, sol.enthalpy_mole, sol.entropy_mole, sol.volume_mole
        mu0 = np.array(sol.standard_gibbs_RT) * R * T
        h0 = np.array(sol.standard_enthalpies_RT) * R * T
        s0 = np.array(sol.standard_entropies_R) * R
        lnX = np.array([math.log(x) if x > 0 else 0.0 for x in X])
        v0 = []
        for k in range(len(names)):
            e = np.zeros(len(names))
            e[k] = 1.0
            setTPX(sol, T, Pr, e)
            v0.append(sol.volume_mole)
        gE, hE, sE, vE = margules_excess(sol, T, Pr, X)
        tag = "P=%.4g" % Pr
        alg(P, "I4.gE[%s]" % tag, g - float(np.dot(X, mu0)) - R * T * float(np.dot(X, lnX)), gE)
        alg(P, "I4.hE[%s]" % tag, h - float(np.dot(X, h0)), hE)
        P.close("I4.sE[%s]" % tag, s - float(np.dot(X, s0)) + R * float(np.dot(X, lnX)), sE, 1e-6, FLOOR["cp"])
        alg_v(P, "I4.vE[%s]" % tag, v - float(np.dot(X, np.array(v0))), vE)
    # Gibbs–Duhem: Σ X_k dlnγ_k/ds = 0 (isothermal-isobaric, normalized-composition path), and Σ_k X_k ∂lnγ_k/∂ln n_j = 0
    Pr = ATM

    def lng(Xc):
        setTPX(sol, T, Pr, Xc)
        return np.log(np.array(sol.activity_coefficients))
    dX = np.array([1.0, -1.0] + [0.0] * (len(names) - 2)) * min(X[0], X[1])
    est = []
    for eps in (1e-4, 1e-5, 1e-6):
        est.append((lng(X + eps * dX) - lng(X - eps * dX)) / (2 * eps))
    gd = float(np.dot(X, est[2]))
    scale = float(np.max(np.abs(est[2])))
    P.cond("I4.Gibbs-Duhem sum X dlngamma/ds=0", abs(gd) <= 1e-6 * scale + 1e-8, obs=gd, exp=0.0, resid=abs(gd))
    n = X.copy()
    for j in range(len(names)):
        if X[j] <= 0:
            continue
        cols = []
        for eps in (1e-4, 1e-5, 1e-6):
            npl, nmi = n.copy(), n.copy()
            npl[j] *= math.exp(eps)
            nmi[j] *= math.exp(-eps)
            cols.append((lng(npl / npl.sum()) - lng(nmi / nmi.sum())) / (2 * eps))
        gdj = float(np.dot(X, cols[2]))
        P.cond("I4.sum_k X_k dlngamma_k/dln n_j=0[%s]" % names[j], abs(gdj) <= 1e-6 * float(np.max(np.abs(cols[2]))) + 1e-8,
               obs=gdj, exp=0.0, resid=abs(gdj))


def t3_I5(sol, P, T, Pr, X, wl):
    """I5 (Python-reachable part): the partial molar c_p of the solvent and of each solute, and the volume and the chemical-potential/partial-molar-enthalpy derivatives, are consistent."""
    rt = 1e-5 if wl else 1e-6
    setTPX(sol, T, Pr, X)
    cpb = np.array(sol.partial_molar_cp)
    for k in range(sol.n_species):
        if X[k] <= 0:
            continue
        deriv3(P, "I5.cpbar=dhbar/dT[%s]" % sol.species_name(k), cpb[k],
               lambda s, k=k: cd_T(sol, T, Pr, X, lambda z: z.partial_molar_enthalpies[k], s * T), HT, rt, FLOOR["cp"])


def tab_table(sol):
    tt = sol.input_data["tabulated-thermo"]
    return np.array(tt["mole-fractions"], float), np.array(tt["enthalpy"], float), np.array(tt["entropy"], float)


def interp_clamped(x, xs, ys):
    if x > xs[-1]:
        return ys[-1]
    if x < xs[0]:
        return ys[0]
    i = int(np.searchsorted(xs, x, side="left"))
    if i == 0:
        return ys[0]
    return ys[i - 1] + (ys[i] - ys[i - 1]) * (x - xs[i - 1]) / (xs[i] - xs[i - 1])


def t3_I6(sol, P, T, Pr, x, simplified, xr=None):
    names = sol.species_names
    itab = names.index(sol.input_data["tabulated-species"])
    X = np.zeros(2)
    if xr is None:
        X[itab], X[1 - itab] = x, 1.0 - x
    else:  # give the trace reference species fraction explicitly, to avoid losing precision in the 1−x subtraction
        X[itab], X[1 - itab] = 1.0, xr
    sol.TPX = T, Pr, X
    xt = float(sol.X[itab])
    q = props(sol)
    for nm in ("h", "s", "g", "cp", "v"):
        P.cond("I6.finite[%s]" % nm, math.isfinite(q[nm]), obs=q[nm])
    P.cond("I6.finite[mu]", bool(np.all(np.isfinite(q["mu"]))), obs=q["mu"].tolist())
    if not all(math.isfinite(q[nm]) for nm in ("h", "s", "g")) or not np.all(np.isfinite(q["mu"])):
        return
    alg(P, "I6.g=h-Ts", q["g"], q["h"] - T * q["s"])
    alg(P, "I6.h=sum X hbar", q["h"], float(np.dot(sol.X, q["hb"])))
    alg(P, "I6.s=sum X sbar", q["s"], float(np.dot(sol.X, q["sb"])))
    alg(P, "I6.g=sum X mu", q["g"], float(np.dot(sol.X, q["mu"])))
    if simplified and abs(Pr - ATM) < 1e-9:
        xs, hs, ss = tab_table(sol)
        dg = q["mu"][itab] - q["mu"][1 - itab]
        alg(P, "I6.mu_tab-mu_ref=h_tab-T s_tab", dg, interp_clamped(xt, xs, hs) - T * interp_clamped(xt, xs, ss))
        if 0 < xt < 1:
            xo = float(sol.X[1 - itab])
            alg(P, "I6.s=x s_tab-R ln(1-x)", q["s"], xt * interp_clamped(xt, xs, ss) - R * math.log(xo))
    # change temperature and return to the same composition/temperature: caches consistent
    sol.TPX = T + 10.0, Pr, X
    sol.TPX = T, Pr, X
    q2 = props(sol)
    for nm in ("h", "s", "g", "cp"):
        P.cond("I6.cache_consistent_after_T_change[%s]" % nm, q2[nm] == q[nm], obs=q2[nm], exp=q[nm])


def t3_eval(p, P):
    key = p["key"]
    kind = t3_kind(key)
    sol = t3_phase(p)
    T, Pr = p["T"], p["P"]
    if kind == "tab":
        simplified = key == "tab_simplified"
        Pe = ATM if simplified else Pr
        t3_I6(sol, P, T, Pe, p["x"], simplified, p.get("xr"))
        return
    X = t3_X(sol, p)
    wl = is_water_like(sol)
    t3_I1(sol, P, T, Pr, X, wl)
    t3_setorder(sol, P, T, Pr, X)
    if kind == "ssvol":
        t3_ssvol(sol, P, T, Pr, X)
    t3_I2(sol, P, T, Pr, X, wl, kind)
    if kind in ("dh", "hmw", "molal"):
        t3_I3(sol, P, T, Pr, X)
    if kind in ("margules_ev", "molten", "ssvol"):
        t3_I4(sol, P, T, X)
    if kind in ("dh", "hmw"):
        t3_I5(sol, P, T, Pr, X, wl)


def track3():
    branches = [(k, t3_sampler(k), t3_eval) for k in T3_LIST]
    must = []
    # discrete values given by the plan: stress molalities 0.1/1/6/9.4 (DH, HMW), table nodes/midpoints/endpoints/trace
    for k in T3_LIST:
        kind = t3_kind(k)
        if kind in ("dh", "hmw"):
            for m in (0.1, 1.0, 6.0, 9.4):
                must.append((k, dict(key=k, T=298.15, P=ATM, m=m, gd_seed=1)))
    for k in ("anode", "tab_simplified"):
        xs = list(np.array(get_phase(k).input_data["tabulated-thermo"]["mole-fractions"], float))
        pts = xs + [0.5 * (a + b) for a, b in zip(xs[:-1], xs[1:])] + [0.0, 1.0, 1e-2, 1e-6, 1e-10, 1e-14, 1e-300]
        for order in (0, 1):
            for x in pts:
                must.append((k, dict(key=k, T=300.0, P=ATM, x=float(x), order=order, gd_seed=1)))
            for xr in (1e-2, 1e-6, 1e-10, 1e-14, 1e-300):  # symmetric trace reference state: give the reference species fraction explicitly
                must.append((k, dict(key=k, T=300.0, P=ATM, x=1.0, xr=xr, order=order, gd_seed=1)))
    return branches, must


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--track", type=int, choices=(1, 2, 3), required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    branches, must = {1: track1, 2: track2, 3: track3}[a.track]()
    versions = dict(cantera=ct.__version__, git_commit=getattr(ct, "__git_commit__", ""),
                    merge_commit="d58a23265414e1789536cf5332d9fb29de095a43", data_dir=DATA,
                    build="env/ctb: env/build_cantera.sh builds from the merge-commit source with scons")
    run(PR, PR_IDX, a.track, a.n, a.out, branches, must=must, versions=versions)


if __name__ == "__main__":
    main()
