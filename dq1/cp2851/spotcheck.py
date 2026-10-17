# spotcheck.py: unified ΔQ1 spot check for CoolProp#2851 (mixture branch of first_saturation_deriv; the CoolProp 8.0.0 wheel in env/py2,
# version ae81610e is later than merge commit 54a2e327 and its two overloads of calc_first_saturation_deriv are verbatim identical to the merge commit, see impl_log.md).
# Usage: python spotcheck.py --track {1,2,3} --n N --out PATH
# Track 1 = track1.md (A tangent along the bubble/dew line + Richardson, B pure-fluid limit, C chain-rule derivative of h, Q=0.5 raises an exception);
# Track 2 = track2.md (I1 total derivative along the phase envelope + reciprocal, I2 z₂→0 limit, I3 saturation derivative of h);
# Track 3 = audit/specs/track3_guided/CoolProp_CoolProp_2851/spec.md §5 (I1 QT/PQ finite differences and fugacity-weighted residual, I2 pure-fluid Clapeyron, I3 chain rule,
# I4 interface convention, I5 invariance under component reordering, I6 second-order saturation derivative).
import argparse
import itertools
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harness import U, choice, run  # noqa: E402

import CoolProp as CP  # noqa: E402
from CoolProp.CoolProp import AbstractState  # noqa: E402

PR = "CoolProp/CoolProp#2851"
PR_IDX = 6


def AS(fluid, z=None):
    A = AbstractState("HEOS", fluid)
    if z is not None and "&" in fluid:
        A.set_mole_fractions(list(z))
    return A


def qt(fluid, z, q, T):
    A = AS(fluid, z)
    A.update(CP.QT_INPUTS, q, T)
    return A


def pq(fluid, z, q, p):
    A = AS(fluid, z)
    A.update(CP.PQ_INPUTS, p, q)
    return A


def richardson(hs, vals):
    """S(h)=S0+a h²+b h⁴ through three points; returns S0."""
    M = np.array([[1.0, h ** 2, h ** 4] for h in hs])
    return float(np.linalg.solve(M, np.array(vals))[0])


def fd_slope(fluid, z, q, T, h):
    return 2 * h / (qt(fluid, z, q, T + h).p() - qt(fluid, z, q, T - h).p())


def fd_keyed_T(fluid, z, q, T, h, key):
    return (qt(fluid, z, q, T + h).keyed_output(key) - qt(fluid, z, q, T - h).keyed_output(key)) / (2 * h)


def flash_ok(fluid, z, T, qs=(0, 1)):
    try:
        for q in qs:
            A = qt(fluid, z, q, T)
            if not (math.isfinite(A.p()) and A.p() > 0):
                return False
        return True
    except Exception:
        return False


# ============================ Track 1 ============================
T1_MIX = {"ME": ("Methane&Ethane", (0.1, 0.9), (120.0, 240.0)), "R32R125": ("R32&R125", (0.2, 0.8), (230.0, 330.0))}


def t1_accept(fluid, z, T):
    """Commonly used operating domain "at least 5 K from critical / cricondentherm": accept only if both the Q=0 and Q=1 flashes succeed at T−5, T and T+5."""
    return all(flash_ok(fluid, z, T + d) for d in (-5.0, 0.0, 5.0))


def t1_sampler(mk, q):
    fluid, zr, tr = T1_MIX[mk]

    def s(rng):
        for _ in range(60):
            z1, T = U(rng, *zr), U(rng, *tr)
            if t1_accept(fluid, [z1, 1 - z1], T):
                return dict(mix=mk, q=q, z1=z1, T=T)
        return dict(mix=mk, q=q, z1=None, T=None)
    return s


def t1_eval(p, P):
    if p["z1"] is None:
        P.fail("sampling", "no admissible state found in 60 draws")
        return
    fluid = T1_MIX[p["mix"]][0]
    z = [p["z1"], 1 - p["z1"]]
    q, T = p["q"], p["T"]
    A = qt(fluid, z, q, T)
    S = A.first_saturation_deriv(CP.iT, CP.iP)
    hs = (0.2, 0.1, 0.05)
    S0 = richardson(hs, [fd_slope(fluid, z, q, T, h) for h in hs])
    P.close("A.dT/dp_sat=Richardson(2h/dp)[Q=%d]" % q, S, S0, 1e-4, 0.0)
    if q == 0:
        D = qt(fluid, z, 0, T).first_saturation_deriv(CP.iHmolar, CP.iT)
        D0 = richardson(hs, [fd_keyed_T(fluid, z, 0, T, h, CP.iHmolar) for h in hs])
        P.close("C.dh/dT_sat(bubble)=Richardson FD", D, D0, 1e-4, 0.0)
    B = qt(fluid, z, 0.5, T)
    P.expect_raise("C.Q=0.5 raises", lambda: B.first_saturation_deriv(CP.iT, CP.iP))


def t1_sampler_B(rng):
    for _ in range(60):
        T = U(rng, 120.0, 240.0)
        if t1_accept("Methane", None, T) and t1_accept("Methane&Ethane", [1 - 1e-2, 1e-2], T):
            return dict(T=T)
    return dict(T=None)


def _linear_rate(P, tag, eps, d):
    for i in range(len(eps) - 1):
        P.close("%s slope d/eps const [%g vs %g]" % (tag, eps[i + 1], eps[i]), d[i + 1] / eps[i + 1], d[i] / eps[i])
    P.cond("%s d decreasing" % tag, all(d[i + 1] < d[i] for i in range(len(d) - 1)), obs=d)


def t1_eval_B(p, P):
    if p["T"] is None:
        P.fail("sampling", "no admissible T")
        return
    T = p["T"]
    Sp = qt("Methane", None, 0, T).first_saturation_deriv(CP.iT, CP.iP)
    eps = (1e-2, 1e-3, 1e-4)
    for q in (0, 1):
        d = [abs(qt("Methane&Ethane", [1 - e, e], q, T).first_saturation_deriv(CP.iT, CP.iP) - Sp) for e in eps]
        _linear_rate(P, "B.pure_limit[Q=%d]" % q, eps, d)


def track1():
    br = [("A_%s_Q%d" % (mk, q), t1_sampler(mk, q), t1_eval) for mk in T1_MIX for q in (0, 1)]
    br.append(("B_pure_limit_ME", t1_sampler_B, t1_eval_B))
    return br, []


# ============================ Track 2 ============================
T2_MIX = {"ME": "Methane&Ethane", "NM": "Nitrogen&Methane", "R32R125": "R32&R125", "PB": "Propane&n-Butane",
          "CM": "CO2&Methane"}


def tc_mix(fluid, z):
    A = AS(fluid, z)
    try:
        return A.T_critical(), "T_critical"
    except Exception:
        return A.T_reducing(), "T_reducing(fallback)"


def t2_sampler(mk, q):
    fluid = T2_MIX[mk]

    def s(rng):
        z1 = U(rng, 0.05, 0.95)
        tc, src = tc_mix(fluid, [z1, 1 - z1])
        return dict(mix=mk, q=q, z1=z1, T=U(rng, 0.6, 0.9) * tc, tc=tc, tc_src=src)
    return s


def t2_eval(p, P):
    fluid = T2_MIX[p["mix"]]
    z = [p["z1"], 1 - p["z1"]]
    q, T = p["q"], p["T"]
    P.record("tc_source", src=p["tc_src"], tc=p["tc"])
    A = qt(fluid, z, q, T)
    S = A.first_saturation_deriv(CP.iT, CP.iP)
    hs = (0.01, 0.02, 0.05)
    S0 = richardson(hs, [fd_slope(fluid, z, q, T, h) for h in hs])
    P.close("I1.dT/dp_sat=Richardson FD[Q=%d]" % q, S, S0, 1e-3, 0.0)
    P.close("I1.dP/dT*dT/dP=1[Q=%d]" % q, A.first_saturation_deriv(CP.iP, CP.iT) * S, 1.0)
    if q == 0:
        try:
            D = qt(fluid, z, 0, T).first_saturation_deriv(CP.iHmolar, CP.iT)
        except Exception as e:
            P.record("I3.finding: first_saturation_deriv(iHmolar,iT) raised at Q=0", err=str(e)[:150])
            P.fail("I3", "exception: %s" % e)
            return
        D0 = richardson(hs, [fd_keyed_T(fluid, z, 0, T, h, CP.iHmolar) for h in hs])
        P.close("I3.dh/dT_sat(bubble)=Richardson FD", D, D0, 1e-3, 0.0)


def t2_sampler_I2(rng):
    tc = AS("Methane").T_critical()
    return dict(T=U(rng, 0.6, 0.9) * tc)


def t2_eval_I2(p, P):
    T = p["T"]
    Sp = qt("Methane", None, 0, T).first_saturation_deriv(CP.iT, CP.iP)
    eps = (1e-2, 1e-3, 1e-4)
    for q in (0, 1):
        d = [abs(qt("Methane&Ethane", [1 - e, e], q, T).first_saturation_deriv(CP.iT, CP.iP) - Sp) for e in eps]
        _linear_rate(P, "I2.pure_limit[Q=%d]" % q, eps, d)


def track2():
    br = [("%s_Q%d" % (mk, q), t2_sampler(mk, q), t2_eval) for mk in T2_MIX for q in (0, 1)]
    br.append(("I2_pure_limit_ME", t2_sampler_I2, t2_eval_I2))
    must = []
    for mk in T2_MIX:
        for q in (0, 1):
            for z1 in (0.05, 0.3, 0.5, 0.7, 0.95):
                tc, src = tc_mix(T2_MIX[mk], [z1, 1 - z1])
                must.append(("%s_Q%d" % (mk, q), dict(mix=mk, q=q, z1=z1, T=0.75 * tc, tc=tc, tc_src=src)))
    return br, must


# ============================ Track 3 ============================
T3_CFG = {"ME": ("Methane&Ethane", "mix"), "MEP": ("Methane&Ethane&Propane", "mix"), "Water": ("Water", "pure"),
          "nPropane": ("n-Propane", "pure"), "R410A": ("R410A", "pseudo")}
MEP_Z = ([0.7, 0.2, 0.1], [0.3, 0.4, 0.3])
DT = (0.05, 0.025, 0.0125)
DPR = (1e-4, 5e-5, 2.5e-5)
FS = {CP.iDmolar: "Dmolar", CP.iHmolar: "Hmolar", CP.iSmolar: "Smolar", CP.iUmolar: "Umolar", CP.iDmass: "Dmass",
      CP.iHmass: "Hmass", CP.iSmass: "Smass", CP.iUmass: "Umass"}


def rel_ok(a, b, r, eabs):
    return abs(a - b) <= r * max(abs(a), abs(b)) + eabs


def t3_sampler(ck, q):
    fluid, kind = T3_CFG[ck]

    def s(rng):
        if ck == "ME":
            z1 = U(rng, 0.1, 0.9)
            return dict(cfg=ck, q=q, T=U(rng, 140, 220), z=[z1, 1 - z1])
        if ck == "MEP":
            return dict(cfg=ck, q=q, T=U(rng, 160, 220), z=list(choice(rng, MEP_Z)))
        rngT = {"Water": (300, 550), "nPropane": (230, 330), "R410A": (250, 310)}[ck]
        return dict(cfg=ck, q=q, T=U(rng, *rngT), z=None)
    return s


def _close3(P, name, a, b, r, eabs):
    if not (math.isfinite(a) and math.isfinite(b)):
        P.fail(name, "nonfinite %r %r" % (a, b))
        return
    P.checks.append(dict(name=name, ok=bool(rel_ok(a, b, r, eabs)), obs=a, exp=b, resid=abs(a - b), tol=None))


def _converged_ref(P, name, ests, r_conv):
    """Reference values at each step size must converge (relative difference between adjacent steps <= r_conv); otherwise "undetermined" is recorded as a fail; returns the finest step."""
    for i in range(len(ests) - 1):
        if not rel_ok(ests[i], ests[i + 1], r_conv, 0.0):
            P.fail(name, "reference not converged: %r" % (ests,))
            return None
    return ests[-1]


def t3_I1(P, fluid, z, q, T, kind):
    A = qt(fluid, z, q, T)
    S = A.first_saturation_deriv(CP.iT, CP.iP)
    p0 = A.p()
    eabs = 1e-8 * T / p0
    r = 1e-3 if kind == "mix" else 1e-4
    fdT = [fd_slope(fluid, z, q, T, h) for h in DT]
    ref = _converged_ref(P, "I1.QT reference", fdT, 5e-4)
    if ref is not None:
        _close3(P, "I1.S=QT central FD", S, ref, r, eabs)
    fd05 = fdT[0]
    P.cond("I1.PR test 0.05K 0.5%", abs(S - fd05) <= 5e-3 * abs(fd05), obs=S, exp=fd05, resid=abs(S - fd05))
    fdP = []
    for rr in DPR:
        dp = rr * p0
        fdP.append((pq(fluid, z, q, p0 + dp).T() - pq(fluid, z, q, p0 - dp).T()) / (2 * dp))
    refp = _converged_ref(P, "I1.PQ reference", fdP, 5e-4)
    if refp is not None:
        _close3(P, "I1.S=PQ central FD", S, refp, r, eabs)
    return A, S, p0


def _homog(fluid, comp, phase, T, p):
    B = AS(fluid, comp)
    B.specify_phase(phase)
    B.update(CP.PT_INPUTS, p, T)
    return B


def t3_I1_weighted(P, fluid, z, q, T, A, S, p0):
    x = list(A.mole_fractions_liquid())
    y = list(A.mole_fractions_vapor())
    P.cond("I1.precond normalized x,y", abs(sum(x) - 1) <= 1e-10 and abs(sum(y) - 1) <= 1e-10, obs=[sum(x), sum(y)])
    rhoL = A.saturated_liquid_keyed_output(CP.iDmolar)
    rhoV = A.saturated_vapor_keyed_output(CP.iDmolar)
    L = _homog(fluid, x, CP.iphase_liquid, T, p0)
    V = _homog(fluid, y, CP.iphase_gas, T, p0)
    P.close("I1.homog liquid density = equilibrium branch", L.rhomolar(), rhoL)
    P.close("I1.homog vapor density = equilibrium branch", V.rhomolar(), rhoV)
    n = len(x)
    lnf = max(abs(math.log(L.fugacity(i) / V.fugacity(i))) for i in range(n))
    P.cond("I1.precond max|ln fL/fV|<1e-6", lnf < 1e-6, obs=lnf)
    w = y if q == 0 else x  # weights: vapor composition on the bubble-point line, liquid composition on the dew-point line (from the total differential of Σ z_i K_i=1 + Gibbs–Duhem)
    res = []
    for k in range(3):
        dT, dpr = DT[k], DPR[k] * p0
        def dln(comp, ph, i, which):
            if which == "T":
                a, b = _homog(fluid, comp, ph, T + dT, p0), _homog(fluid, comp, ph, T - dT, p0)
                return (math.log(a.fugacity_coefficient(i)) - math.log(b.fugacity_coefficient(i))) / (2 * dT), a, b
            a, b = _homog(fluid, comp, ph, T, p0 + dpr), _homog(fluid, comp, ph, T, p0 - dpr)
            return (math.log(a.fugacity_coefficient(i)) - math.log(b.fugacity_coefficient(i))) / (2 * dpr), a, b
        AT = Ap = 0.0
        for i in range(n):
            lT, a1, b1 = dln(x, CP.iphase_liquid, i, "T")
            vT, a2, b2 = dln(y, CP.iphase_gas, i, "T")
            lp, a3, b3 = dln(x, CP.iphase_liquid, i, "p")
            vp, a4, b4 = dln(y, CP.iphase_gas, i, "p")
            for st, ref_rho in ((a1, rhoL), (b1, rhoL), (a3, rhoL), (b3, rhoL), (a2, rhoV), (b2, rhoV), (a4, rhoV), (b4, rhoV)):
                if abs(st.rhomolar() / ref_rho - 1) > 0.05:  # root switch: the finite difference at this step is invalid
                    P.fail("I1.weighted residual", "root switch in stencil")
                    return
            AT += w[i] * (lT - vT)
            Ap += w[i] * (lp - vp)
        den = abs(AT * S) + abs(Ap)
        res.append(abs(AT * S + Ap) / den if den > 0 else float("nan"))
    if not math.isfinite(res[-1]):
        P.fail("I1.weighted residual", "zero denominator")
        return
    P.record("I1.weighted residual by step", res=res)
    P.cond("I1.|A_T S + A_p|/(|A_T S|+|A_p|)<2e-3", res[-1] < 2e-3, obs=res[-1], exp=0.0, resid=res[-1])


def t3_I2(P, fluid, q, T, kind):
    A = qt(fluid, None, q, T)
    S = A.first_saturation_deriv(CP.iT, CP.iP)
    rL, rV = A.saturated_liquid_keyed_output(CP.iDmolar), A.saturated_vapor_keyed_output(CP.iDmolar)
    hL, hV = A.saturated_liquid_keyed_output(CP.iHmolar), A.saturated_vapor_keyed_output(CP.iHmolar)
    expr = T * (1 / rV - 1 / rL) / (hV - hL)
    _close3(P, "I2.S=T dv/dh (saturated keyed outputs)", S, expr, 1e-10, 1e-8 * T / A.p())
    if kind == "pure":
        Sother = qt(fluid, None, 1 - q, T).first_saturation_deriv(CP.iT, CP.iP)
        _close3(P, "I2.bubble slope = dew slope (pure)", S, Sother, 1e-8, 1e-8 * T / A.p())
    else:
        P.record("I2.R410A base-vs-merge comparison not run (no isolated base build on node)")


def t3_I3(P, fluid, z, q, T, kind, S, p0):
    r = 1e-3 if kind == "mix" else 1e-4
    A = qt(fluid, z, q, T)
    for F, nm in FS.items():
        try:
            _t3_I3_one(P, fluid, z, q, T, F, nm, r, S, p0)
        except Exception as e:
            P.fail("I3.%s" % nm, "%s: %s" % (type(e).__name__, str(e)[:120]))
    A = qt(fluid, z, q, T)
    _t3_I3_tail(P, A)


def _t3_I3_one(P, fluid, z, q, T, F, nm, r, S, p0):
    if True:
        A = qt(fluid, z, q, T)
        Fc = A.keyed_output(F)
        fdT, fdP = [], []
        for h in DT:
            a, b = qt(fluid, z, q, T + h), qt(fluid, z, q, T - h)
            fdT.append((a.keyed_output(F) - b.keyed_output(F)) / (2 * h))
        for rr in DPR:
            dp = rr * p0
            a, b = pq(fluid, z, q, p0 + dp), pq(fluid, z, q, p0 - dp)
            fdP.append((a.keyed_output(F) - b.keyed_output(F)) / (2 * dp))
        Fstar = max(1.0, abs(Fc))
        A = qt(fluid, z, q, T)
        DTF = A.first_saturation_deriv(F, CP.iT)
        DpF = A.first_saturation_deriv(F, CP.iP)
        refT = _converged_ref(P, "I3.%s T-ref" % nm, fdT, 5e-4)
        if refT is not None:
            _close3(P, "I3.d%s/dT_sat=FD" % nm, DTF, refT, r, 1e-8 * Fstar / T)
        refP = _converged_ref(P, "I3.%s p-ref" % nm, fdP, 5e-4)
        if refP is not None:
            _close3(P, "I3.d%s/dp_sat=FD" % nm, DpF, refP, r, 1e-8 * Fstar / p0)
        try:
            FT = A.first_partial_deriv(F, CP.iT, CP.iP)
            Fp = A.first_partial_deriv(F, CP.iP, CP.iT)
        except Exception as e:
            P.fail("I3.first_partial_deriv", str(e)[:120])
            return
        _close3(P, "I3.chain D_pF=F_p+F_T S[%s]" % nm, DpF, Fp + FT * S, 1e-10, 1e-8 * Fstar / p0)
        _close3(P, "I3.chain D_TF=F_T+F_p/S[%s]" % nm, DTF, FT + Fp / S, 1e-10, 1e-8 * Fstar / T)
        _close3(P, "I3.chain D_TF=D_pF/S[%s]" % nm, DTF, DpF / S, 1e-10, 1e-8 * Fstar / T)


def _t3_I3_tail(P, A):
    for nm, f in (("I3.D_pT*D_Tp=1", lambda: A.first_saturation_deriv(CP.iT, CP.iP) * A.first_saturation_deriv(CP.iP, CP.iT)),
                  ("I3.D_TT=1", lambda: A.first_saturation_deriv(CP.iT, CP.iT)),
                  ("I3.D_pp=1", lambda: A.first_saturation_deriv(CP.iP, CP.iP))):
        try:
            _close3(P, nm, f(), 1.0, 0.0, 1e-10)
        except Exception as e:
            P.fail(nm, "%s: %s" % (type(e).__name__, str(e)[:120]))
    P.expect_raise("I4.unsupported Wrt1 (iT,iDmolar) raises", lambda: A.first_saturation_deriv(CP.iT, CP.iDmolar))


def t3_I5(P, fluid, z, q, T, S, p0):
    names = fluid.split("&")
    A = qt(fluid, z, q, T)
    x0, y0 = list(A.mole_fractions_liquid()), list(A.mole_fractions_vapor())
    for perm in itertools.permutations(range(len(names))):
        if list(perm) == list(range(len(names))):
            continue
        f2 = "&".join(names[i] for i in perm)
        z2 = [z[i] for i in perm]
        B = qt(f2, z2, q, T)
        _close3(P, "I5.relabel p_sat[%s]" % f2, B.p(), p0, 1e-7, 0.0)
        _close3(P, "I5.relabel slope[%s]" % f2, B.first_saturation_deriv(CP.iT, CP.iP), S, 1e-7, 0.0)
        x2, y2 = list(B.mole_fractions_liquid()), list(B.mole_fractions_vapor())
        for j, i in enumerate(perm):
            P.cond("I5.relabel x[%s,%s]" % (f2, names[i]), abs(x2[j] - x0[i]) <= 1e-9, obs=x2[j], exp=x0[i])
            P.cond("I5.relabel y[%s,%s]" % (f2, names[i]), abs(y2[j] - y0[i]) <= 1e-9, obs=y2[j], exp=y0[i])


def t3_I6(P, fluid, z, q, T, p0):
    for F, nm in ((CP.iT, "T"), (CP.iDmolar, "Dmolar"), (CP.iHmolar, "Hmolar"), (CP.iSmolar, "Smolar"), (CP.iUmolar, "Umolar")):
        A = qt(fluid, z, q, T)
        try:
            an = A.second_saturation_deriv(F, CP.iP, CP.iP)
        except Exception as e:
            P.record("I6.second_saturation_deriv unsupported", F=nm, err=str(e)[:120])
            continue
        ests = []
        for rr in DPR:
            dp = rr * p0
            a, b = pq(fluid, z, q, p0 + dp), pq(fluid, z, q, p0 - dp)
            ests.append((a.first_saturation_deriv(F, CP.iP) - b.first_saturation_deriv(F, CP.iP)) / (2 * dp))
        Fstar = max(1.0, abs(A.keyed_output(F)))
        ref = _converged_ref(P, "I6.%s reference" % nm, ests, 1e-3)
        if ref is not None:
            _close3(P, "I6.d2%s/dp2_sat=FD of first_saturation_deriv" % nm, an, ref, 5e-3, 1e-8 * Fstar / p0 ** 2)


def t3_eval(p, P):
    fluid, kind = T3_CFG[p["cfg"]]
    z, q, T = p["z"], p["q"], p["T"]
    A, S, p0 = t3_I1(P, fluid, z, q, T, kind)
    if kind == "mix":
        t3_I1_weighted(P, fluid, z, q, T, A, S, p0)
    else:
        t3_I2(P, fluid, q, T, kind)
    t3_I3(P, fluid, z, q, T, kind, S, p0)
    if kind == "mix":
        t3_I5(P, fluid, z, q, T, S, p0)
    t3_I6(P, fluid, z, q, T, p0)
    if p.get("rejection_test"):
        for Qi in (0.1, 0.5, 0.9):
            B = qt("Methane&Ethane", [0.7, 0.3], Qi, 180.0)
            for (o, w), nm in (((CP.iT, CP.iP), "T,P"), ((CP.iP, CP.iT), "P,T"), ((CP.iHmolar, CP.iT), "H,T"),
                               ((CP.iHmolar, CP.iP), "H,P")):
                P.expect_raise("I4.interior Q=%g rejects (%s)" % (Qi, nm), lambda o=o, w=w, B=B: B.first_saturation_deriv(o, w))


def track3():
    br = [("%s_Q%d" % (ck, q), t3_sampler(ck, q), t3_eval) for ck in T3_CFG for q in (0, 1)]
    must = []
    for q in (0, 1):
        for T in range(140, 221, 10):
            for z1 in (0.1, 0.3, 0.5, 0.7, 0.9):
                must.append(("ME_Q%d" % q, dict(cfg="ME", q=q, T=float(T), z=[z1, 1 - z1])))
        for T in range(160, 221, 10):
            for zz in MEP_Z:
                must.append(("MEP_Q%d" % q, dict(cfg="MEP", q=q, T=float(T), z=list(zz))))
        for ck, Ts in (("Water", (300, 400, 500, 550)), ("nPropane", (230, 270, 300, 330)), ("R410A", (250, 280, 310))):
            for T in Ts:
                must.append(("%s_Q%d" % (ck, q), dict(cfg=ck, q=q, T=float(T), z=None)))
    must.append(("ME_Q0", dict(cfg="ME", q=0, T=180.0, z=[0.7, 0.3], rejection_test=True)))
    return br, must


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--track", type=int, choices=(1, 2, 3), required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    branches, must = {1: track1, 2: track2, 3: track3}[a.track]()
    versions = dict(CoolProp=CP.__version__, gitrevision=getattr(CP, "__gitrevision__", ""),
                    merge_commit="54a2e32764bfa0ae81877e835c33ccd13ba8f853",
                    note="CoolProp 8.0.0 wheel (ae81610e, 291 commits after merge); calc_first_saturation_deriv identical to merge")
    run(PR, PR_IDX, a.track, a.n, a.out, branches, must=must, versions=versions)


if __name__ == "__main__":
    main()
