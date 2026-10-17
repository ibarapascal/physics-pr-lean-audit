"""CoolProp2851 replay: native check of CoolProp PR #2851 (mixture first_saturation_deriv).

Merge commit 54a2e32764bfa0ae81877e835c33ccd13ba8f853; the patched function is identical in
tag v8.0.0 (HelmholtzEOSMixtureBackend.cpp:3832-3891), so the PyPI wheel CoolProp==8.0.0 is the
native target.

Dependencies: Python 3, CoolProp==8.0.0, numpy (only for medians).
Run:  python CoolProp2851_replay.py [--n 2000] [--seed 0] > CoolProp2851_replay.json

What it does (all output is one JSON object on stdout):
  For random states on the bubble (Q=0) and dew (Q=1) curves of binary mixtures
  (Methane&Ethane, R32&R125) it compares
  * native  = AS.first_saturation_deriv(iT, iP)
  * lean    = the Lean definition dTdP_sat_mixture (CoolProp2851.weightedDiff quotient),
              re-evaluated from fugacity coefficients of the coexisting phases; the partial
              derivatives (d ln phi_i/dp)_{T,n}, (d ln phi_i/dT)_{p,n} are obtained by central
              differences in (T, rho) of AS.fugacity_coefficient(i) on single-phase states with the
              saturated-phase composition and density, converted with analytic dp/drho, dp/dT
  * fd      = Richardson-extrapolated central difference of T along the same curve
              (theorems bubble_branch_is_tangent / dew_branch_is_tangent)
  * clapeyron_old = T*(1/rhoV - 1/rhoL)/(hV - hL) (pre-PR expression, for reference)
  and checks:
  * reciprocal: first_saturation_deriv(iP,iT) * first_saturation_deriv(iT,iP) == 1
  * chain:      first_saturation_deriv(iHmolar,iT) == first_saturation_deriv(iHmolar,iP)/dTdP
  * hmolar_along_curve: first_saturation_deriv(iHmolar,iT) vs central difference of the bulk
                hmolar along the curve (reported, no Lean theorem; bulk phase at Q=0 is liquid)
  * pure limit: z2 -> 0 sequence vs pure-fluid first_saturation_deriv (Clapeyron branch)
  * throw at Q=0.5 for a mixture
"""
import argparse
import json
import math
import random
import sys

try:
    import CoolProp
    import CoolProp.CoolProp as CP
except Exception as exc:  # pragma: no cover
    print(json.dumps({"status": "SKIP", "reason": "CoolProp import failed: %r" % (exc,)}))
    sys.exit(0)

try:
    import numpy as np
except Exception:  # pragma: no cover
    np = None

MIXTURES = [
    # name, fluids string, z1 range, T range [K]
    ("Methane&Ethane", "Methane&Ethane", (0.1, 0.9), (120.0, 240.0)),
    ("R32&R125", "R32&R125", (0.2, 0.8), (230.0, 330.0)),
]


def median(xs):
    if not xs:
        return None
    if np is not None:
        return float(np.median(xs))
    s = sorted(xs)
    return s[len(s) // 2]


def make_state(fluids, z):
    AS = CP.AbstractState("HEOS", fluids)
    AS.set_mole_fractions(list(z))
    return AS


def lnphi_partials(fluids, comp, rho, T, phase):
    """(d ln phi_i/dp)_{T,n} and (d ln phi_i/dT)_{p,n} at (T, rho, comp) for the given phase."""
    AS = make_state(fluids, comp)
    AS.specify_phase(phase)
    n = len(comp)
    hr = 1e-5 * rho
    hT = 1e-4 * T

    def lnphis(r, t):
        AS.update(CP.DmolarT_INPUTS, r, t)
        return [math.log(AS.fugacity_coefficient(i)) for i in range(n)]

    AS.update(CP.DmolarT_INPUTS, rho, T)
    dpdrho_T = AS.first_partial_deriv(CP.iP, CP.iDmolar, CP.iT)
    dpdT_rho = AS.first_partial_deriv(CP.iP, CP.iT, CP.iDmolar)
    lr_p = lnphis(rho + hr, T)
    lr_m = lnphis(rho - hr, T)
    lt_p = lnphis(rho, T + hT)
    lt_m = lnphis(rho, T - hT)
    dl_drho = [(a - b) / (2 * hr) for a, b in zip(lr_p, lr_m)]
    dl_dT_rho = [(a - b) / (2 * hT) for a, b in zip(lt_p, lt_m)]
    dl_dp_T = [d / dpdrho_T for d in dl_drho]
    dl_dT_p = [dt - dr * dpdT_rho / dpdrho_T for dt, dr in zip(dl_dT_rho, dl_drho)]
    return dl_dp_T, dl_dT_p


def weighted_diff(w, dL, dV):
    # CoolProp2851.weightedDiff
    return sum(wi * (a - b) for wi, a, b in zip(w, dL, dV))


def lean_dTdP_sat_mixture(Q, x, y, dLp, dVp, dLT, dVT):
    # CoolProp2851.dTdP_sat_mixture ; None models the throw
    if Q == 0:
        return -weighted_diff(y, dLp, dVp) / weighted_diff(y, dLT, dVT)
    if Q == 1:
        return -weighted_diff(x, dLp, dVp) / weighted_diff(x, dLT, dVT)
    return None


def psat(fluids, z, Q, T):
    AS = make_state(fluids, z)
    AS.update(CP.QT_INPUTS, Q, T)
    return AS.p(), AS.hmolar()


def fd_slope(fluids, z, Q, T, h=0.1):
    def D(hh):
        pp, _ = psat(fluids, z, Q, T + hh)
        pm, _ = psat(fluids, z, Q, T - hh)
        return 2 * hh / (pp - pm)

    d1, d2 = D(h), D(h / 2)
    return (4 * d2 - d1) / 3


def fd_h_along(fluids, z, Q, T, h=0.1):
    def D(hh):
        _, hp = psat(fluids, z, Q, T + hh)
        _, hm = psat(fluids, z, Q, T - hh)
        return (hp - hm) / (2 * hh)

    d1, d2 = D(h), D(h / 2)
    return (4 * d2 - d1) / 3


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def one_point(fluids, z, Q, T):
    AS = make_state(fluids, z)
    AS.update(CP.QT_INPUTS, Q, T)
    native = AS.first_saturation_deriv(CP.iT, CP.iP)
    native_dpdT = AS.first_saturation_deriv(CP.iP, CP.iT)
    x = list(AS.mole_fractions_liquid())
    y = list(AS.mole_fractions_vapor())
    rhoL = AS.saturated_liquid_keyed_output(CP.iDmolar)
    rhoV = AS.saturated_vapor_keyed_output(CP.iDmolar)
    hL = AS.saturated_liquid_keyed_output(CP.iHmolar)
    hV = AS.saturated_vapor_keyed_output(CP.iHmolar)
    if abs(rhoL - rhoV) / rhoL < 0.2:
        raise ValueError("near-critical state excluded")
    dLp, dLT = lnphi_partials(fluids, x, rhoL, T, CP.iphase_liquid)
    dVp, dVT = lnphi_partials(fluids, y, rhoV, T, CP.iphase_gas)
    lean = lean_dTdP_sat_mixture(Q, x, y, dLp, dVp, dLT, dVT)
    fd = fd_slope(fluids, z, Q, T)
    clap = T * (1 / rhoV - 1 / rhoL) / (hV - hL)
    rec = {
        "native": native,
        "lean": lean,
        "fd": fd,
        "clapeyron_old": clap,
        "rel_native_lean": rel(native, lean),
        "rel_native_fd": rel(native, fd),
        "rel_old_fd": rel(clap, fd),
        "reciprocal_err": abs(native * native_dpdT - 1.0),
    }
    try:
        dhdT = AS.first_saturation_deriv(CP.iHmolar, CP.iT)
        dhdp = AS.first_saturation_deriv(CP.iHmolar, CP.iP)
        rec["chain_err"] = rel(dhdT, dhdp / native)
        rec["rel_hmolar_along_curve"] = rel(dhdT, fd_h_along(fluids, z, Q, T))
    except Exception as exc:
        rec["chain_exception"] = str(exc)[:200]
    return rec


def summarize(vals):
    vals = [v for v in vals if v is not None and math.isfinite(v)]
    if not vals:
        return {"n": 0}
    return {"n": len(vals), "max": max(vals), "median": median(vals)}


def pure_limit(fluids, pure, T, Q):
    AP = CP.AbstractState("HEOS", pure)
    AP.update(CP.QT_INPUTS, Q, T)
    ref = AP.first_saturation_deriv(CP.iT, CP.iP)
    rows = []
    for eps in (1e-2, 1e-3, 1e-4):
        try:
            AS = make_state(fluids, [1 - eps, eps])
            AS.update(CP.QT_INPUTS, Q, T)
            v = AS.first_saturation_deriv(CP.iT, CP.iP)
            rows.append({"z2": eps, "value": v, "rel_diff_to_pure": rel(v, ref)})
        except Exception as exc:
            rows.append({"z2": eps, "exception": str(exc)[:200]})
    return {"pure": pure, "T": T, "Q": Q, "pure_value": ref, "sequence": rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    out = {
        "module": "CoolProp2851",
        "coolprop_version": CoolProp.__version__,
        "coolprop_gitrevision": getattr(CoolProp, "__gitrevision__", None),
        "n_requested": args.n,
        "mixtures": {},
    }
    per_mix = max(1, args.n // len(MIXTURES))
    examples = []
    for name, fluids, zr, Tr in MIXTURES:
        recs = {0: [], 1: []}
        failures = 0
        for k in range(per_mix):
            Q = k % 2
            z1 = rng.uniform(*zr)
            T = rng.uniform(*Tr)
            try:
                r = one_point(fluids, [z1, 1 - z1], Q, T)
                r.update({"z1": z1, "T": T, "Q": Q})
                recs[Q].append(r)
            except Exception:
                failures += 1
        mix = {"n_attempted": per_mix, "n_failed_or_excluded": failures}
        for Q in (0, 1):
            rs = recs[Q]
            key = "bubble_Q0" if Q == 0 else "dew_Q1"
            mix[key] = {
                "n": len(rs),
                "rel_native_vs_lean": summarize([r["rel_native_lean"] for r in rs]),
                "rel_native_vs_fd_tangent": summarize([r["rel_native_fd"] for r in rs]),
                "rel_old_clapeyron_vs_fd_tangent": summarize([r["rel_old_fd"] for r in rs]),
                "reciprocal_abs_err": summarize([r["reciprocal_err"] for r in rs]),
                "chain_rel_err": summarize([r.get("chain_err") for r in rs]),
                "rel_hmolar_satderiv_vs_fd": summarize([r.get("rel_hmolar_along_curve") for r in rs]),
                "n_chain_exceptions": sum(1 for r in rs if "chain_exception" in r),
                "worst_native_vs_fd": max(rs, key=lambda r: r["rel_native_fd"]) if rs else None,
            }
            if rs:
                examples.append(rs[0])
        out["mixtures"][name] = mix

    out["pure_limit"] = [
        pure_limit("Methane&Ethane", "Methane", 150.0, 0),
        pure_limit("Methane&Ethane", "Methane", 150.0, 1),
        pure_limit("R32&R125", "R32", 280.0, 0),
        pure_limit("R32&R125", "R32", 280.0, 1),
    ]

    # Q = 0.5 on a mixture must raise (theorem mixture_intermediate_Q_throws)
    try:
        AS = make_state("Methane&Ethane", [0.7, 0.3])
        AS.update(CP.QT_INPUTS, 0.5, 180.0)
        v = AS.first_saturation_deriv(CP.iT, CP.iP)
        out["intermediate_Q_throws"] = {"raised": False, "value": v}
    except Exception as exc:
        out["intermediate_Q_throws"] = {"raised": True, "message": str(exc)[:200]}

    # PR test point
    try:
        out["pr_test_point_methane_ethane_0.7_180K_Q0"] = one_point("Methane&Ethane", [0.7, 0.3], 0, 180.0)
        out["pr_test_point_methane_ethane_0.7_180K_Q1"] = one_point("Methane&Ethane", [0.7, 0.3], 1, 180.0)
    except Exception as exc:
        out["pr_test_point_error"] = str(exc)[:200]

    # Informational: two-phase derivatives that internally call SatL/SatV->calc_first_saturation_deriv
    info = {}
    for label, call in (
        ("first_two_phase_deriv_Dmolar_Hmolar_P", lambda A: A.first_two_phase_deriv(CP.iDmolar, CP.iHmolar, CP.iP)),
        ("second_saturation_deriv_T_P_P", lambda A: A.second_saturation_deriv(CP.iT, CP.iP, CP.iP)),
    ):
        for Q in (0.0, 0.5, 1.0):
            try:
                A = make_state("Methane&Ethane", [0.7, 0.3])
                A.update(CP.QT_INPUTS, Q, 180.0)
                info["%s_Q%g" % (label, Q)] = {"value": call(A)}
            except Exception as exc:
                info["%s_Q%g" % (label, Q)] = {"exception": str(exc)[:200]}
    out["informational_dependent_routines"] = info
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
