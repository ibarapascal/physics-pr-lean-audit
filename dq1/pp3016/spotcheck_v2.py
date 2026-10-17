# spotcheck_v2.py: formal Track 3 (v2) ΔQ1 spot check of pandapower#3016 (merge commit 0c6a1f0a).
# Plan = "Numerical check method" in audit/specs/track3/e2nIEE_pandapower_3016/spec.md, together with Invariants (I1–I3) and Common operating domain of the same section.
# Branches (k=6, points split evenly): I1_Dyn, I1_YNd (reversed fixture), I1_YNyn (supplementary fixture, absolute T-circuit oracle), I2_ZNyn, I2_ZNd, I2_ZNy (three-bus zigzag fixture).
# Within each point: both inverse_y values are run and compared with each other, S_B∈{1,10,100} are compared with each other; the I3 current formula (once with the independent-oracle Z0, once with the public Z0), the I1ph/I3ph ratio, monotonicity and upper bound of the increasing R_N sweep,
# the large-NER bound at t=0.01, the central-difference derivative (h, h/2, 2h stabilized, then compared with the analytic derivative).
# Usage: python spotcheck_v2.py --n N --out PATH; seed 20261017+8×10+3; output JSON same as spotcheck.py, plus track="3", plan="codex_v2".
# The sentence-by-sentence plan-to-code mapping and implementation judgments are in impl_log.md §8 in the same directory. Dependencies: pandapower (env/pp3016), numpy; reuses the network-building / solving utilities of harness.py and spotcheck.py.
import argparse
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harness import U, run  # noqa: E402
import spotcheck as S  # noqa: E402  (importing sets warnings/logging to silent)

import numpy as np  # noqa: E402

PLAN = "codex_v2"
ZR, ZA = 1e-7, 1e-8      # Ω: |Ẑ−Z*| ≤ 1e-8 Ω + 1e-7|Z*|
IR, IA = 1e-7, 1e-8      # kA: |Î−I*| ≤ 1e-8 kA + 1e-7|I*|
IOPEN = 1e-6             # kA: open-circuit current threshold (the 1e-6 of the newly added test)
C11 = 1.1                # case=max, lv_tol 10%: c=1.1 at every voltage level
EPS3 = np.finfo(float).eps ** (1.0 / 3.0)
SBS = (10.0, 100.0)      # S_B=1 is the reference


# ============================ Common ============================
def zkr(p, u, pre=""):
    """Raw zero-sequence short-circuit impedance (Ω, referred to side u): U²/S_T·(r0+j√(p0²−r0²))/100."""
    return S.zleak(p[pre + "vk0"], p[pre + "vkr0"], u, p[pre + "sn"])


def kt(p, pre="", vk=None, vkr=None):
    """Plan formula: K_T=0.95·1.1/(1+0.6√(p²−r²)/100)."""
    vk = p[pre + "vk"] if vk is None else vk
    vkr = p[pre + "vkr"] if vkr is None else vkr
    return 0.95 * C11 / (1 + 0.6 * math.sqrt(vk ** 2 - vkr ** 2) / 100.0)


def zmag(p, zk_abs, pre=""):
    return S.zmag(p[pre + "mag0"], p[pre + "mag0_rx"], zk_abs)


def solve(net, fault="1ph", inv=True, sb=1.0):
    net.sn_mva = sb
    try:
        return S.sc(net, "max", inv, fault, 10, use_pre_fault_voltage=False, r_fault_ohm=0.0, x_fault_ohm=0.0,
                    branch_results=False)
    finally:
        net.sn_mva = 1.0


def fin_c(z):
    return math.isfinite(z.real) and math.isfinite(z.imag)


def zclose(P, name, obs, exp):
    P.cclose(name, obs, exp, ZR, ZA, comp=False)


def iclose(P, name, obs, exp):
    P.close(name, obs, exp, IR, IA)


def cmp_bus(P, tag, ra, rb):
    """Auxiliary consistency: Z0, Z1 compared by complex modulus, ikss by current tolerance; if both sides are infinity (ideal open circuit) no inf−inf is formed, and one side finite with the other infinite is recorded as a fail."""
    for nm, fz in (("Z0", S.z0), ("Z1", S.z1)):
        a, b = fz(ra), fz(rb)
        if fin_c(a) and fin_c(b):
            zclose(P, "%s.%s" % (tag, nm), a, b)
        elif fin_c(a) != fin_c(b):
            P.fail("%s.%s" % (tag, nm), "finite/infinite mismatch %r vs %r" % (a, b))
    iclose(P, tag + ".ikss", ra["ikss_ka"], rb["ikss_ka"])


def neutral(p, zk_abs):
    """3R_N/|Z_k,earth|, 3X_N/|Z_k,earth| ∈ [0.1,100] (uniform); mandatory points give RN/XN (Ω) directly."""
    if "RN" in p:
        return p["RN"], p["XN"]
    return p["kR"] * zk_abs / 3.0, p["kX"] * zk_abs / 3.0


def setzn(net, t, rn, xn):
    S.setzn(net, t, rn, xn)


# ---- I3 ----
def i3_formula(P, tag, r, u, z0pred=None):
    """ikss against √3·1.1·U_b/|2Z1+Z0|: once with Z0 from the independent oracle (if any), once with the public Z0 (internal consistency)."""
    z1 = S.z1(r)
    if z0pred is not None and fin_c(z1):
        iclose(P, "I3.ikss_vs_pred_Z0[%s]" % tag, r["ikss_ka"], S.SQ3 * C11 * u / abs(2 * z1 + z0pred))
    z0 = S.z0(r)
    if fin_c(z0) and fin_c(z1):
        iclose(P, "I3.ikss_internal[%s]" % tag, r["ikss_ka"], S.SQ3 * C11 * u / abs(2 * z1 + z0))


def i3_ratio(P, tag, r1, r3):
    """I_1ph/I_3ph = 3|Z1|/|2Z1+Z0| (the 3ph result on the same network cancels the voltage factor)."""
    z1, z0 = S.z1(r1), S.z0(r1)
    if fin_c(z0) and fin_c(z1) and r3["ikss_ka"] > 0:
        P.close("I3.I1ph/I3ph=3|Z1|/|2Z1+Z0|[%s]" % tag, r1["ikss_ka"] / r3["ikss_ka"], 3 * abs(z1) / abs(2 * z1 + z0), IR, IA)


def i3_sweep(P, tag, net, t, bus, u, R, X):
    """Increasing R_N sweep at fixed X_N: monotonically non-increasing, I_f ≤ cU/(√3R_N); at t=|D|/(3R_N)=0.01, 1/(1+t) ≤ I_f/(C/3R_N) ≤ 1; central-difference derivative."""
    C = S.SQ3 * C11 * u

    def ik_at(rn):
        setzn(net, t, rn, X)
        return solve(net)[bus]["ikss_ka"]

    setzn(net, t, 0.0, X)
    r0 = solve(net)[bus]
    D = 2 * S.z1(r0) + S.z0(r0)
    if not fin_c(D):
        P.fail("I3.sweep[%s]" % tag, "D nonfinite")
        return
    A, B = D.real, D.imag
    if A < 0 or B < 0:
        P.record("I3.precondition_A,B>=0_not_met", tag=tag, A=A, B=B)
        return
    rlim = abs(D) / (3 * 0.01)
    rs = sorted(set([0.0, R / 2, R, 2 * R, rlim]))
    iks = [r0["ikss_ka"] if rn == 0.0 else ik_at(rn) for rn in rs]
    P.cond("I3.monotone_RN[%s]" % tag, all(iks[i + 1] <= iks[i] + IA + IR * iks[i] for i in range(len(iks) - 1)),
           obs=list(zip(rs, iks)))
    for rn, ik in zip(rs, iks):
        if rn > 0:
            ub = C / (3 * rn)
            P.cond("I3.upper_bound_cU/(sqrt3RN)[%s,RN=%.6g]" % (tag, rn), ik <= ub + IA + IR * ub, obs=ik, exp=ub)
    ik_l = iks[rs.index(rlim)]
    if ik_l < IOPEN:
        P.fail("I3.limit[%s]" % tag, "open-circuit truncation range ik=%r" % ik_l)
    else:
        ub = C / (3 * rlim)
        lb = ub / 1.01
        P.cond("I3.large_NER_bound_t=0.01[%s]" % tag, lb - IA - IR * lb <= ik_l <= ub + IA + IR * ub, obs=ik_l,
               exp=[lb, ub], resid=ik_l * 3 * rlim / C)
    # central difference (optional derivative check): h ~ eps^(1/3)·max(R_N,|D|/3); requires R_N>h (2h is also used, hence requires R_N>2h)
    h = EPS3 * max(R, abs(D) / 3)
    if R > 2 * h:
        ds = []
        for hh in (h, h / 2, 2 * h):
            ds.append((ik_at(R + hh) - ik_at(R - hh)) / (2 * hh))
        spread = max(abs(d - ds[0]) for d in ds)
        if not spread <= 1e-6 * abs(ds[0]) + 1e-12:
            P.fail("I3.FD_unstable[%s]" % tag, "FD h,h/2,2h = %r" % ds)
        else:
            ana = -3 * C * (A + 3 * R) / ((A + 3 * R) ** 2 + B ** 2) ** 1.5
            P.close("I3.dIf/dRN_FD_vs_analytic[%s]" % tag, ds[0], ana)
    else:
        P.record("I3.FD_skipped_RN<=2h", tag=tag, R=R, h=h)
    setzn(net, t, R, X)


def aux_checks(P, tag, net, buses, ref_by_inv, states_nets):
    """Both inverse_y values compared with each other, S_B∈{1,10,100} compared with each other (two states: Z_N=0 and R+jX)."""
    for key, (rn, xn, t) in states_nets.items():
        setzn(net, t, rn, xn)
        for inv in (True, False):
            for sb in SBS:
                r = solve(net, inv=inv, sb=sb)
                for b in buses:
                    cmp_bus(P, "SB_cov[%s,%s,inv=%s,S_B=%g,b%d]" % (tag, key, inv, sb, b), r[b], ref_by_inv[(inv, key)][b])
        for b in buses:
            cmp_bus(P, "inverse_y_equal[%s,%s,b%d]" % (tag, key, b), ref_by_inv[(False, key)][b], ref_by_inv[(True, key)][b])


# ============================ Sampling (Common operating domain) ============================
def s_volts(rng):
    while True:
        uh, ul = U(rng, 10, 110), U(rng, 0.4, 33)
        if uh > ul:
            return uh, ul


def s_tr(rng, pre=""):
    vk, vk0 = U(rng, 4, 15), U(rng, 4, 15)
    return {pre + "sn": U(rng, 0.4, 40), pre + "vk": vk, pre + "vkr": U(rng, 0.02, 0.4) * vk, pre + "vk0": vk0,
            pre + "vkr0": U(rng, 0.02, 0.4) * vk0, pre + "mag0": U(rng, 10, 1000), pre + "mag0_rx": U(rng, 0, 1),
            pre + "si0": U(rng, 0.1, 0.9)}


def s_grid(rng, pre="g_"):
    return {pre + "ssc": U(rng, 100, 1e4), pre + "rx": U(rng, 0.05, 0.5), pre + "x0x": U(rng, 0.5, 3),
            pre + "r0x0": U(rng, 0.05, 0.5)}


def s_i1(rng):
    p = {}
    p["uh"], p["ul"] = s_volts(rng)
    p.update(s_tr(rng))
    vk2 = U(rng, 4, 15)
    p.update(vk_b=vk2, vkr_b=U(rng, 0.02, 0.4) * vk2)
    p.update(s_grid(rng))
    p.update(kR=U(rng, 0.1, 100), kX=U(rng, 0.1, 100))
    return p


def s_ynyn(rng):
    p = {}
    p["uh"], p["ul"] = s_volts(rng)
    p.update(s_tr(rng))
    p.update(rho=U(rng, 0.05, 0.5), g_ssc=U(rng, 100, 1e4))
    p.update(kR=U(rng, 0.1, 100), kX=U(rng, 0.1, 100))
    return p


def s_zz(rng):
    p = {}
    p["uh"], p["ul"] = s_volts(rng)
    p.update(s_tr(rng, "m_"))
    p.update(s_tr(rng))
    p.update(s_grid(rng))
    p.update(kR=U(rng, 0.1, 100), kX=U(rng, 0.1, 100), mag0_alt=U(rng, 10, 1000), mag0_rx_alt=U(rng, 0, 1),
             a_ssc=U(rng, 100, 1e4), a_rx=U(rng, 0.05, 0.5), a_x0x1=U(rng, 0.5, 3), a_r0x01=U(rng, 0.05, 0.5),
             a_x0x2=U(rng, 0.5, 3), a_r0x02=U(rng, 0.05, 0.5))
    return p


# ============================ Fixtures ============================
def _tr(net, hv, lv, p, vg, uh, ul, pre="", vk=None, vkr=None):
    return S.trafo(net, hv, lv, p[pre + "sn"], uh, ul, p[pre + "vk"] if vk is None else vk,
                   p[pre + "vkr"] if vkr is None else vkr, vg, vk0=p[pre + "vk0"], vkr0=p[pre + "vkr0"],
                   mag0=p[pre + "mag0"], mag0_rx=p[pre + "mag0_rx"], si0=p[pre + "si0"], rn=0.0, xn=0.0)


def _grid(net, bus, p, pre="g_", x0x=None, r0x0=None):
    return S.grid(net, bus, p[pre + "ssc"], rx=p[pre + "rx"], x0x=p[pre + "x0x"] if x0x is None else x0x,
                  r0x0=p[pre + "r0x0"] if r0x0 is None else r0x0)


def fx_i1(p, vg, vk, vkr):
    """Dyn: ext_grid on HV, observe LV; YNd: ext_grid on LV, observe HV."""
    net = S.pp.create_empty_network(sn_mva=1.0)
    bh = S.pp.create_bus(net, vn_kv=p["uh"])
    bl = S.pp.create_bus(net, vn_kv=p["ul"])
    _grid(net, bh if vg == "Dyn" else bl, p)
    t = _tr(net, bh, bl, p, vg, p["uh"], p["ul"], vk=vk, vkr=vkr)
    return net, t, (bl if vg == "Dyn" else bh), (p["ul"] if vg == "Dyn" else p["uh"])


def fx_ynyn(p):
    net = S.pp.create_empty_network(sn_mva=1.0)
    bh = S.pp.create_bus(net, vn_kv=p["uh"])
    bl = S.pp.create_bus(net, vn_kv=p["ul"])
    S.grid(net, bl, p["g_ssc"], rx=p["rho"], x0x=1.0, r0x0=p["rho"])
    t = _tr(net, bh, bl, p, "YNyn", p["uh"], p["ul"])
    return net, t, bh


def fx_zz(p, vg, aux=None):
    """110 kV ext_grid → 110/U_H YNd main transformer → U_H bus → U_H/U_L zigzag → aux; aux has no device by default, and an ext_grid is added when aux=(x0x,r0x0)."""
    net = S.pp.create_empty_network(sn_mva=1.0)
    bs = S.pp.create_bus(net, vn_kv=110.0)
    bh = S.pp.create_bus(net, vn_kv=p["uh"])
    ba = S.pp.create_bus(net, vn_kv=p["ul"])
    _grid(net, bs, p)
    _tr(net, bs, bh, p, "YNd", 110.0, p["uh"], pre="m_")
    t = _tr(net, bh, ba, p, vg, p["uh"], p["ul"])
    if aux is not None:
        S.grid(net, ba, p["a_ssc"], rx=p["a_rx"], x0x=aux[0], r0x0=aux[1])
    return net, t, bh, ba


# ============================ I1: Dyn / YNd ============================
def make_ev_i1(vg):
    def ev(p, P):
        u_e = p["ul"] if vg == "Dyn" else p["uh"]
        R, X = neutral(p, abs(zkr(p, u_e)))
        states = {"0": (0.0, 0.0), "R": (R, 0.0), "jX": (0.0, X), "RX": (R, X)}
        for gi, (vk, vkr) in enumerate(((p["vk"], p["vkr"]), (p["vk_b"], p["vkr_b"]))):
            net, t, fb, u = fx_i1(p, vg, vk, vkr)
            tg = "%s,g%d" % (vg, gi)
            res, r3 = {}, {}
            for inv in (True, False):
                for key, (rn, xn) in states.items():
                    setzn(net, t, rn, xn)
                    res[(inv, key)] = solve(net, inv=inv)
                z00 = S.z0(res[(inv, "0")][fb])
                for key in ("R", "jX", "RX"):
                    zclose(P, "I1.dZ0=3dZN[%s,%s,inv=%s]" % (tg, key, inv), S.z0(res[(inv, key)][fb]) - z00,
                           3 * complex(*states[key]))
                zclose(P, "I1.Z1_unchanged[%s,inv=%s]" % (tg, inv), S.z1(res[(inv, "RX")][fb]), S.z1(res[(inv, "0")][fb]))
                for fault in ("3ph", "2ph"):
                    rr = {}
                    for key in ("0", "RX"):
                        setzn(net, t, *states[key])
                        rr[key] = solve(net, fault, inv)[fb]
                        if inv and fault == "3ph":
                            r3[key] = rr[key]
                    iclose(P, "I1.%s_ikss_unchanged[%s,inv=%s]" % (fault, tg, inv), rr["RX"]["ikss_ka"], rr["0"]["ikss_ka"])
            z00 = S.z0(res[(True, "0")][fb])
            for key, (rn, xn) in states.items():
                i3_formula(P, "%s,%s" % (tg, key), res[(True, key)][fb], u,
                           None if key == "0" else z00 + 3 * complex(rn, xn))
            for key in ("0", "RX"):
                i3_ratio(P, "%s,%s" % (tg, key), res[(True, key)][fb], r3[key])
            if gi == 0:
                aux_checks(P, tg, net, [fb], {k: v for k, v in res.items() if k[1] in ("0", "RX")},
                           {"0": (0.0, 0.0, t), "RX": (R, X, t)})
                i3_sweep(P, tg, net, t, fb, u, R, X)
    return ev


# ============================ I1: YNyn supplementary fixture ============================
def ynyn_pred(p, zn):
    zk_r = zkr(p, p["uh"])
    zk = kt(p) * zk_r
    zm = zmag(p, abs(zk_r))
    a = p["si0"]
    rho = p["rho"]
    zg0 = C11 * p["uh"] ** 2 * (rho + 1j) / (p["g_ssc"] * math.sqrt(1 + rho ** 2))
    rest = (1 - a) * zk + zg0
    return a * zk + 3 * zn + zm * rest / (zm + rest)


def ev_ynyn(p, P):
    R, X = neutral(p, abs(zkr(p, p["uh"])))
    states = {"0": (0.0, 0.0), "R": (R, 0.0), "jX": (0.0, X), "RX": (R, X)}
    net, t, fb = fx_ynyn(p)
    u = p["uh"]
    tg = "YNyn"
    res, r3 = {}, {}
    for inv in (True, False):
        for key, (rn, xn) in states.items():
            setzn(net, t, rn, xn)
            res[(inv, key)] = solve(net, inv=inv)
            zclose(P, "I1.YNyn_T_pred_Z0H[%s,inv=%s]" % (key, inv), S.z0(res[(inv, key)][fb]), ynyn_pred(p, complex(rn, xn)))
        z00 = S.z0(res[(inv, "0")][fb])
        for key in ("R", "jX", "RX"):
            zclose(P, "I1.dZ0H=3dZN[%s,%s,inv=%s]" % (tg, key, inv), S.z0(res[(inv, key)][fb]) - z00, 3 * complex(*states[key]))
        zclose(P, "I1.Z1_unchanged[%s,inv=%s]" % (tg, inv), S.z1(res[(inv, "RX")][fb]), S.z1(res[(inv, "0")][fb]))
        for fault in ("3ph", "2ph"):
            rr = {}
            for key in ("0", "RX"):
                setzn(net, t, *states[key])
                rr[key] = solve(net, fault, inv)[fb]
                if inv and fault == "3ph":
                    r3[key] = rr[key]
            iclose(P, "I1.%s_ikss_unchanged[%s,inv=%s]" % (fault, tg, inv), rr["RX"]["ikss_ka"], rr["0"]["ikss_ka"])
    for key, (rn, xn) in states.items():
        i3_formula(P, "%s,%s" % (tg, key), res[(True, key)][fb], u, ynyn_pred(p, complex(rn, xn)))
    for key in ("0", "RX"):
        i3_ratio(P, "%s,%s" % (tg, key), res[(True, key)][fb], r3[key])
    aux_checks(P, tg, net, [fb], {k: v for k, v in res.items() if k[1] in ("0", "RX")},
               {"0": (0.0, 0.0, t), "RX": (R, X, t)})
    i3_sweep(P, tg, net, t, fb, u, R, X)


# ============================ I2: ZNyn / ZNd / ZNy ============================
def zz_pred(p, zn, mag=None):
    """HV Ω: Z0,H = α·K_T·Z_k,r + 3Z_N; LV Ω of ZNyn: Z0,L = (Z_m + (1−α)K_T Z_k,r)·(U_L/U_H)²."""
    zk_r = zkr(p, p["uh"])
    zk = kt(p) * zk_r
    a = p["si0"]
    q = dict(p) if mag is None else dict(p, mag0=mag[0], mag0_rx=mag[1])
    zm = zmag(q, abs(zk_r))
    return a * zk + 3 * zn, (zm + (1 - a) * zk) * (p["ul"] / p["uh"]) ** 2


def make_ev_zz(vg):
    others = [g for g in ("ZNyn", "ZNd", "ZNy") if g != vg]

    def lv_check(P, tag, r_ba, pred_l):
        if vg == "ZNyn":
            zclose(P, "I2.Z0L_pred[%s]" % tag, S.z0(r_ba), pred_l)
        else:
            ik = r_ba["ikss_ka"]
            if math.isnan(ik):
                P.fail("I2.aux_If[%s]" % tag, "ikss NaN")
            else:
                P.cond("I2.aux_If<=1e-6kA[%s]" % tag, abs(ik) <= IOPEN, obs=ik, exp=0.0)

    def ev(p, P):
        R, X = neutral(p, abs(zkr(p, p["uh"])))
        states = {"0": (0.0, 0.0), "R": (R, 0.0), "jX": (0.0, X), "RX": (R, X)}
        net, t, bh, ba = fx_zz(p, vg)
        u = p["uh"]
        res, r3 = {}, {}
        for inv in (True, False):
            for key, (rn, xn) in states.items():
                setzn(net, t, rn, xn)
                r = solve(net, inv=inv)
                res[(inv, key)] = r
                ph, pl = zz_pred(p, complex(rn, xn))
                zclose(P, "I2.Z0H_pred[%s,%s,inv=%s]" % (vg, key, inv), S.z0(r[bh]), ph)
                lv_check(P, "%s,%s,inv=%s" % (vg, key, inv), r[ba], pl)
            if vg == "ZNyn":
                zclose(P, "I2.ZNyn_Z0L_indep_HV_neutral[inv=%s]" % inv, S.z0(res[(inv, "RX")][ba]), S.z0(res[(inv, "0")][ba]))
        # magnetization perturbation: HV Z0 unchanged; ZNyn LV predicted with the new magnetization
        setzn(net, t, R, X)
        net.trafo.at[t, "mag0_percent"] = p["mag0_alt"]
        net.trafo.at[t, "mag0_rx"] = p["mag0_rx_alt"]
        rm = solve(net)
        net.trafo.at[t, "mag0_percent"] = p["mag0"]
        net.trafo.at[t, "mag0_rx"] = p["mag0_rx"]
        zclose(P, "I2.Z0H_indep_magnetizing[%s]" % vg, S.z0(rm[bh]), S.z0(res[(True, "RX")][bh]))
        lv_check(P, "%s,mag_alt" % vg, rm[ba], zz_pred(p, complex(R, X), (p["mag0_alt"], p["mag0_rx_alt"]))[1])
        # secondary-winding connection perturbation: HV Z0 does not change with yn/d/y
        for g2 in others:
            net.trafo.at[t, "vector_group"] = g2
            r2 = solve(net)
            zclose(P, "I2.Z0H_indep_secondary[%s->%s]" % (vg, g2), S.z0(r2[bh]), S.z0(res[(True, "RX")][bh]))
        net.trafo.at[t, "vector_group"] = vg
        # ext_grid added on aux, positive-sequence s_sc/rx fixed, only x0x/r0x0 changed: HV Z0 unchanged
        ra = []
        for aux in ((p["a_x0x1"], p["a_r0x01"]), (p["a_x0x2"], p["a_r0x02"])):
            na, ta, bha, _baa = fx_zz(p, vg, aux=aux)
            setzn(na, ta, R, X)
            ra.append(solve(na)[bha])
        zclose(P, "I2.Z0H_indep_LV_zero_seq_terminal[%s]" % vg, S.z0(ra[1]), S.z0(ra[0]))
        # I3
        for key, (rn, xn) in states.items():
            ph, pl = zz_pred(p, complex(rn, xn))
            i3_formula(P, "%s_HV,%s" % (vg, key), res[(True, key)][bh], u, ph)
            if vg == "ZNyn":
                i3_formula(P, "%s_LV,%s" % (vg, key), res[(True, key)][ba], p["ul"], pl)
        for key in ("0", "RX"):
            setzn(net, t, *states[key])
            r3[key] = solve(net, "3ph")[bh]
            i3_ratio(P, "%s_HV,%s" % (vg, key), res[(True, key)][bh], r3[key])
        aux_checks(P, vg, net, [bh, ba], {k: v for k, v in res.items() if k[1] in ("0", "RX")},
                   {"0": (0.0, 0.0, t), "RX": (R, X, t)})
        i3_sweep(P, vg + "_HV", net, t, bh, u, R, X)
    return ev


# ============================ Mandatory points (discrete values stated in the plan) ============================
def must_points():
    base = dict(sn=0.63, vk=6.0, vkr=1.0, vk0=6.0, vkr0=1.0, mag0=100.0, mag0_rx=0.0, si0=0.9)
    zns = ((0.01, 0.02), (5.0, 0.0))   # the PR's 0.01 Ω, 0.02 Ω, 5 Ω
    out = []
    for rn, xn in zns:
        out.append(("I1_Dyn", dict(base, uh=20.0, ul=0.4, vk_b=4.0, vkr_b=1.0, g_ssc=1000.0, g_rx=0.1, g_x0x=1.0,
                                   g_r0x0=0.1, RN=rn, XN=xn)))
        out.append(("I1_YNd", dict(uh=110.0, ul=20.0, sn=40.0, vk=12.0, vkr=0.4, vk0=12.0, vkr0=0.4, mag0=100.0,
                                   mag0_rx=0.0, si0=0.9, vk_b=6.0, vkr_b=1.0, g_ssc=5000.0, g_rx=0.1, g_x0x=1.0,
                                   g_r0x0=0.1, RN=rn, XN=xn)))
        out.append(("I1_YNyn", dict(base, uh=20.0, ul=0.4, rho=0.1, g_ssc=1000.0, RN=rn, XN=xn)))
        zz = dict(uh=20.0, ul=0.4, m_sn=40.0, m_vk=12.0, m_vkr=0.4, m_vk0=12.0, m_vkr0=0.4, m_mag0=100.0, m_mag0_rx=0.0,
                  m_si0=0.9, sn=0.4, vk=4.0, vkr=1.0, vk0=4.0, vkr0=1.0, mag0=100.0, mag0_rx=0.0, si0=0.9,
                  g_ssc=5000.0, g_rx=0.1, g_x0x=1.0, g_r0x0=0.1, mag0_alt=1000.0, mag0_rx_alt=0.5, a_ssc=1000.0,
                  a_rx=0.1, a_x0x1=1.0, a_r0x01=0.1, a_x0x2=3.0, a_r0x02=0.5, RN=rn, XN=xn)
        for vg in ("ZNyn", "ZNd", "ZNy"):
            out.append(("I2_" + vg, dict(zz)))
    return out


def branches():
    return [("I1_Dyn", s_i1, make_ev_i1("Dyn")), ("I1_YNd", s_i1, make_ev_i1("YNd")), ("I1_YNyn", s_ynyn, ev_ynyn),
            ("I2_ZNyn", s_zz, make_ev_zz("ZNyn")), ("I2_ZNd", s_zz, make_ev_zz("ZNd")), ("I2_ZNy", s_zz, make_ev_zz("ZNy"))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    import pandapower
    versions = dict(pandapower=pandapower.__version__, numpy=np.__version__,
                    merge_commit="0c6a1f0a65273af91dd1d5b854498cd89d457939",
                    # Chinese text kept verbatim: this string is written to the archived outputs/pp3016_t*.json; meaning: "env/pp3016: uv pip install pandapower from git at 0c6a1f0a (env/setup_pca_py.sh)"
                    install="env/pp3016：uv pip install pandapower@git+...@0c6a1f0a（env/setup_pca_py.sh）")
    run(S.PR, S.PR_IDX, 3, a.n, a.out, branches(), must=must_points(), versions=versions,
        info=dict(plan=PLAN, spec="specs_codex_v2/e2nIEE_pandapower_3016/spec.md"))
    with open(a.out) as f:
        d = json.load(f)
    d["track"] = "3"
    d["plan"] = PLAN
    with open(a.out + ".tmp", "w") as f:
        json.dump(d, f, ensure_ascii=False, indent=1)
    os.replace(a.out + ".tmp", a.out)


if __name__ == "__main__":
    main()
