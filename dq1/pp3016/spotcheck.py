# spotcheck.py: unified ΔQ1 spot check for pandapower#3016 (merge commit 0c6a1f0a) -- implemented according to each of the three track specifications' own "numerical verification method".
# Usage: python spotcheck.py --track {1,2,3} --n N --out PATH
# Track 1 = track1.md (A Dyn/YNd, B IEC current formula and monotonicity, C zigzag transformer, plus pf_3ph neutral-point displacement);
# Track 2 = track2.md (I1 3Z_N and Möbius fit, I2 current formula / monotonicity / asymptotics, I3 zigzag decoupling / Zd / ZNd low-voltage side);
# Track 3 = audit/specs/track3_guided/e2nIEE_pandapower_3016/spec.md §5 (I1–I6, PSU, T-model path, pf_3ph).
# The sentence-by-sentence plan-to-code comparison and all implementation judgments are in impl_log.md in the same directory. Dependencies: pandapower (env/pp3016), numpy;
# the network of the Track 3 PSU case is taken from the t1 + generator data of iec_60909_4() in the merge-commit test file test_iec60909_4.py (parameters copied verbatim, the test file is not imported).
import argparse
import cmath
import logging
import math
import os
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")
logging.disable(logging.CRITICAL)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harness import LU, U, choice, run  # noqa: E402

import pandapower as pp  # noqa: E402
from pandapower.shortcircuit import calc_sc  # noqa: E402

PR = "e2nIEE/pandapower#3016"
PR_IDX = 8
SQ3 = math.sqrt(3.0)
NAN = float("nan")
SHIFT = {"dyn": 330.0, "ynd": 330.0, "yzn": 330.0, "yd": 330.0, "dy": 330.0, "znyn": 330.0, "znd": 330.0,
         "zny": 330.0, "zd": 330.0}


# ============================ Common utilities ============================
def grid(net, bus, ssc=1000.0, ssc_min=None, rx=0.1, x0x=1.0, r0x0=0.1, rx_min=None, x0x_min=None, r0x0_min=None):
    return pp.create_ext_grid(net, bus, s_sc_max_mva=ssc, s_sc_min_mva=ssc if ssc_min is None else ssc_min,
                              rx_max=rx, rx_min=rx if rx_min is None else rx_min,
                              x0x_max=x0x, x0x_min=x0x if x0x_min is None else x0x_min,
                              r0x0_max=r0x0, r0x0_min=r0x0 if r0x0_min is None else r0x0_min)


def trafo(net, hv, lv, sn, vhv, vlv, vk, vkr, vg, vk0=None, vkr0=None, mag0=100.0, mag0_rx=0.0, si0=0.9,
          rn=0.0, xn=0.0, parallel=1, **kw):
    return pp.create_transformer_from_parameters(
        net, hv, lv, sn_mva=sn, vn_hv_kv=vhv, vn_lv_kv=vlv, vk_percent=vk, vkr_percent=vkr, pfe_kw=0.0,
        i0_percent=0.0, vector_group=vg, shift_degree=SHIFT.get(vg.lower(), 0.0),
        vk0_percent=vk if vk0 is None else vk0, vkr0_percent=vkr if vkr0 is None else vkr0,
        mag0_percent=mag0, mag0_rx=mag0_rx, si0_hv_partial=si0, rn_ohm=rn, xn_ohm=xn, parallel=parallel, **kw)


def sc(net, case, inv=True, fault="1ph", lvtol=10, **kw):
    calc_sc(net, fault=fault, case=case, inverse_y=inv, lv_tol_percent=lvtol, **kw)
    r = net.res_bus_sc
    return {int(b): {c: float(r.at[b, c]) for c in r.columns} for b in r.index}


def setzn(net, t, rn, xn):
    net.trafo.at[t, "rn_ohm"] = float(rn)
    net.trafo.at[t, "xn_ohm"] = float(xn)


def z0(r):
    return complex(r["rk0_ohm"], r["xk0_ohm"])


def z1(r):
    return complex(r["rk_ohm"], r["xk_ohm"])


def cfac(u_kv, case, lvtol=10):
    """Voltage factors of IEC 60909-0:2016 Table 1 (Track 3 §5 I3 states the same values): <1 kV: 10% tolerance 1.10/0.90, 6% tolerance 1.05/0.95; otherwise 1.10/1.00."""
    if u_kv < 1.0:
        cmax, cmin = (1.1, 0.9) if lvtol == 10 else (1.05, 0.95)
    else:
        cmax, cmin = 1.1, 1.0
    return cmax if case == "max" else cmin


def ik1(u_kv, c, zz1, zz0):
    return SQ3 * c * u_kv / abs(2 * zz1 + zz0)


def zleak(vk0, vkr0, u_kv, sn, p=1):
    """Transformer zero-sequence leakage reactance (Ω, referred to the side with voltage u_kv): (vkr0 + j·sqrt(vk0²−vkr0²))/100 · U²/S /p."""
    return (vkr0 + 1j * math.sqrt(vk0 ** 2 - vkr0 ** 2)) / 100.0 * u_kv ** 2 / sn / p


def zmag(mag0, rx, zl_abs):
    """Zero-sequence magnetizing impedance (Ω): per the create documentation, mag0_percent=(z_mag0/z0)·100%, mag0_rx=R/X. zl_abs is the uncorrected magnitude of the zero-sequence short-circuit impedance."""
    return mag0 / 100.0 * zl_abs * (rx + 1j) / math.sqrt(1 + rx ** 2)


def kt_iec(vk, vkr, cmax_lv):
    """IEC 60909-0 Eq. (12a): K_T = 0.95 c_max /(1 + 0.6 x_T), where x_T is the relative positive-sequence reactance."""
    return 0.95 * cmax_lv / (1 + 0.6 * math.sqrt(vk ** 2 - vkr ** 2) / 100.0)


def zgrid0(u_kv, ssc, rx, x0x, r0x0, c):
    """External-grid zero-sequence impedance (Ω): Z_Q=c·U²/S'', X1=Z_Q/sqrt(1+rx²), X0=x0x·X1, R0=r0x0·X0."""
    zq = c * u_kv ** 2 / ssc
    x1 = zq / math.sqrt(1 + rx ** 2)
    x0 = x0x * x1
    return complex(r0x0 * x0, x0)


def par(a, b):
    return a * b / (a + b)


def fin(x):
    return x is not None and math.isfinite(x)


# ============================ Track 1 ============================
def t1_tr(rng, pre=""):
    vk = U(rng, 4, 14)
    return {pre + "sn": U(rng, 0.1, 100), pre + "vk": vk, pre + "vkr": U(rng, 0.3, 1.5), pre + "si0": U(rng, 0.5, 0.95)}


def t1_s_dyn(rng):
    p = dict(mv=U(rng, 10, 35), RN=U(rng, 0, 1000), XN=U(rng, 0, 100), inv=bool(rng.integers(2)))
    p.update(t1_tr(rng))
    return p


def t1_s_ynd(rng):
    p = dict(hv=U(rng, 110, 132), mv=U(rng, 10, 35), RN=U(rng, 0, 1000), XN=U(rng, 0, 100), inv=bool(rng.integers(2)))
    p.update(t1_tr(rng))
    return p


def t1_s_zz(rng):
    p = dict(hv=U(rng, 110, 132), mv=U(rng, 10, 35), RN=U(rng, 0, 1000), XN=U(rng, 0, 100), inv=bool(rng.integers(2)))
    p.update(t1_tr(rng, "m_"))
    p.update(t1_tr(rng, "e_"))
    return p


def t1_s_pf(rng):
    p = dict(mv=U(rng, 10, 35), RN=U(rng, 0, 1000), XN=U(rng, 0, 100))
    p.update(t1_tr(rng))
    return p


def _scan(net, t, bus, case, inv, RN, XN):
    out = {}
    for key, (rn, xn) in (("0", (0.0, 0.0)), ("x", (0.0, XN)), ("h", (RN / 2, XN)), ("f", (RN, XN))):
        setzn(net, t, rn, xn)
        out[key] = sc(net, case, inv)[bus]
    return out


def t1_e_dyn(p, P):
    net = pp.create_empty_network(sn_mva=1.0)
    bm = pp.create_bus(net, vn_kv=p["mv"])
    bl = pp.create_bus(net, vn_kv=0.4)
    grid(net, bm, 1000.0)
    t = trafo(net, bm, bl, p["sn"], p["mv"], 0.4, p["vk"], p["vkr"], "Dyn", si0=p["si0"])
    zn = complex(p["RN"], p["XN"])
    for case in ("max", "min"):
        r = _scan(net, t, bl, case, p["inv"], p["RN"], p["XN"])
        P.cclose("A.dZ0=3ZN[%s]" % case, z0(r["f"]) - z0(r["0"]), 3 * zn, 1e-9, 1e-9, "max", comp=False)
        P.cclose("A.Z1_indep[%s]" % case, z1(r["f"]), z1(r["0"]), 1e-9, 1e-9, "max", comp=False)
        c = cfac(0.4, case)
        for key, rr in r.items():
            P.close("B.ikss=sqrt3cU/|2Z1+Z0|[%s,%s]" % (case, key), rr["ikss_ka"], ik1(0.4, c, z1(rr), z0(rr)), 1e-9, 0.0, "max")
        ik = [r[k]["ikss_ka"] for k in ("x", "h", "f")]
        P.cond("B.monotone_RN[%s]" % case, ik[1] <= ik[0] * (1 + 1e-9) and ik[2] <= ik[1] * (1 + 1e-9), obs=ik)


def t1_e_ynd(p, P):
    net = pp.create_empty_network(sn_mva=1.0)
    bh = pp.create_bus(net, vn_kv=p["hv"])
    bm = pp.create_bus(net, vn_kv=p["mv"])
    grid(net, bh, 5000.0)
    t = trafo(net, bh, bm, p["sn"], p["hv"], p["mv"], p["vk"], p["vkr"], "YNd", si0=p["si0"])
    for case in ("max", "min"):
        net.trafo.at[t, "in_service"] = False
        zrest = z0(sc(net, case, p["inv"])[bh])
        net.trafo.at[t, "in_service"] = True
        r = _scan(net, t, bh, case, p["inv"], p["RN"], p["XN"])
        zns = {"0": 0j, "x": complex(0, p["XN"]), "h": complex(p["RN"] / 2, p["XN"]), "f": complex(p["RN"], p["XN"])}
        const = {k: 1 / (1 / z0(r[k]) - 1 / zrest) - 3 * zns[k] for k in r}
        for k in ("x", "h", "f"):
            P.cclose("A.YNd_parallel_const[%s,%s]" % (case, k), const[k], const["0"], 1e-9, 0.0, "max", comp=False)
        P.cclose("A.Z1_indep[%s]" % case, z1(r["f"]), z1(r["0"]), 1e-9, 1e-9, "max", comp=False)
        c = cfac(p["hv"], case)
        for key, rr in r.items():
            P.close("B.ikss=sqrt3cU/|2Z1+Z0|[%s,%s]" % (case, key), rr["ikss_ka"], ik1(p["hv"], c, z1(rr), z0(rr)), 1e-9, 0.0, "max")
        ik = [r[k]["ikss_ka"] for k in ("x", "h", "f")]
        P.cond("B.monotone_RN[%s]" % case, ik[1] <= ik[0] * (1 + 1e-9) and ik[2] <= ik[1] * (1 + 1e-9), obs=ik)


def t1_e_zz(p, P):
    net = pp.create_empty_network(sn_mva=1.0)
    bh = pp.create_bus(net, vn_kv=p["hv"])
    bm = pp.create_bus(net, vn_kv=p["mv"])
    ba = pp.create_bus(net, vn_kv=0.4)
    g = grid(net, bh, 5000.0)
    trafo(net, bh, bm, p["m_sn"], p["hv"], p["mv"], p["m_vk"], p["m_vkr"], "YNd", si0=p["m_si0"])
    te = trafo(net, bm, ba, p["e_sn"], p["mv"], 0.4, p["e_vk"], p["e_vkr"], "ZNyn", si0=p["e_si0"])
    zn = complex(p["RN"], p["XN"])
    kt = kt_iec(p["e_vk"], p["e_vkr"], 1.1)
    zl = zleak(p["e_vk"], p["e_vkr"], p["mv"], p["e_sn"])
    res = {}
    for vg in ("ZNyn", "ZNd"):
        net.trafo.at[te, "vector_group"] = vg
        for case in ("max", "min"):
            for key, (rn, xn) in (("0", (0.0, 0.0)), ("f", (p["RN"], p["XN"]))):
                setzn(net, te, rn, xn)
                res[(vg, case, key)] = sc(net, case, p["inv"])
            if vg == "ZNyn":
                net.ext_grid.loc[g, ["x0x_max", "x0x_min", "r0x0_max", "r0x0_min"]] = [2.0, 2.0, 0.2, 0.2]
                res[(vg, case, "alt")] = sc(net, case, p["inv"])
                net.ext_grid.loc[g, ["x0x_max", "x0x_min", "r0x0_max", "r0x0_min"]] = [1.0, 1.0, 0.1, 0.1]
    for case in ("max", "min"):
        for key, zz in (("0", 0j), ("f", zn)):
            a, b = res[("ZNyn", case, key)][bm], res[("ZNd", case, key)][bm]
            P.cond("C.ZNyn_vs_ZNd_hv_bitwise[%s,%s]" % (case, key),
                   a["rk0_ohm"] == b["rk0_ohm"] and a["xk0_ohm"] == b["xk0_ohm"], obs=[a["rk0_ohm"], a["xk0_ohm"]],
                   exp=[b["rk0_ohm"], b["xk0_ohm"]])
            if case == "max":
                for vg in ("ZNyn", "ZNd"):
                    P.cclose("C.Z0hv=si0*KT*Z0leak+3ZN[%s,%s]" % (vg, key), z0(res[(vg, case, key)][bm]),
                             p["e_si0"] * kt * zl + 3 * zz, comp=False)
            P.cond("C.ZNd_lv_ikss=0[%s,%s]" % (case, key), abs(res[("ZNd", case, key)][ba]["ikss_ka"]) <= 1e-12,
                   obs=res[("ZNd", case, key)][ba]["ikss_ka"], exp=0.0)
        for vg in ("ZNyn", "ZNd"):
            P.cclose("C.dZ0hv=3ZN[%s,%s]" % (vg, case), z0(res[(vg, case, "f")][bm]) - z0(res[(vg, case, "0")][bm]),
                     3 * zn, comp=False)
        P.cclose("C.ZNyn_lv_Z0_indep_of_grid0[%s]" % case, z0(res[("ZNyn", case, "alt")][ba]),
                 z0(res[("ZNyn", case, "f")][ba]), comp=False)


def _v0_pu(net, bus):
    r = net.res_bus_3ph
    v = [r.at[bus, "vm_%s_pu" % ph] * cmath.exp(1j * math.radians(r.at[bus, "va_%s_degree" % ph])) for ph in "abc"]
    return sum(v) / 3.0


def t1_e_pf(p, P):
    # Post-smoke-test fix (impl_log "post-smoke-test fix" S1): the load is taken as 0.01·V_ph²/(|Z_N|+|Z_T|) (neutral-point displacement about 1%, a single-phase load can supply it),
    # the network base is 100 times the load (the outer convergence criterion 3e-8 of runpp_3ph is in per unit); the positive-sequence impedance of the external line is 1e-6·(|Z_N|+|Z_T|) Ω; the convergence tolerance is the library default.
    zt = abs(zleak(p["vk"], p["vkr"], 0.4, p["sn"]))
    zloop = abs(complex(p["RN"], p["XN"])) + zt
    s_load = 0.01 * (0.4 ** 2 / 3.0) / zloop
    pl, ql = s_load * 2 / math.sqrt(5), s_load / math.sqrt(5)
    z1_line = 1e-6 * zloop / math.sqrt(2)
    nets = []
    for variant in ("A", "B"):
        net = pp.create_empty_network(sn_mva=100.0 * s_load)
        bm = pp.create_bus(net, vn_kv=p["mv"])
        bl = pp.create_bus(net, vn_kv=0.4)
        grid(net, bm, 1000.0)
        if variant == "A":
            trafo(net, bm, bl, p["sn"], p["mv"], 0.4, p["vk"], p["vkr"], "Dyn", si0=p["si0"], rn=p["RN"], xn=p["XN"])
            bload = bl
        else:
            trafo(net, bm, bl, p["sn"], p["mv"], 0.4, p["vk"], p["vkr"], "Dyn", si0=p["si0"], rn=0.0, xn=0.0)
            bload = pp.create_bus(net, vn_kv=0.4)
            pp.create_line_from_parameters(net, bl, bload, length_km=1.0, r_ohm_per_km=z1_line, x_ohm_per_km=z1_line,
                                           c_nf_per_km=0.0, max_i_ka=10.0, r0_ohm_per_km=3 * p["RN"],
                                           x0_ohm_per_km=3 * p["XN"], c0_nf_per_km=0.0)
        pp.create_asymmetric_load(net, bload, p_a_mw=pl, q_a_mvar=ql)
        pp.runpp_3ph(net)
        if not net.converged:
            P.fail("pf3ph", "not converged " + variant)
            return
        nets.append(_v0_pu(net, bload))
    P.cclose("ADD.pf3ph_V0_Dyn+ZN_vs_external_ZN", nets[0], nets[1], comp=False)


def track1():
    mid_dyn = dict(mv=22.5, sn=50.05, vk=9.0, vkr=0.9, si0=0.725, inv=True)
    mid_ynd = dict(hv=121.0, mv=22.5, sn=50.05, vk=9.0, vkr=0.9, si0=0.725, inv=True)
    must = []
    for rn, xn in ((0.0, 0.0), (0.01, 0.0), (0.1, 0.05), (10.0, 0.0)):  # the discrete Z_N values stated in Track 1 A
        must.append(("A_Dyn", dict(mid_dyn, RN=rn, XN=xn)))
        must.append(("A_YNd", dict(mid_ynd, RN=rn, XN=xn)))
    branches = [("A_Dyn", t1_s_dyn, t1_e_dyn), ("A_YNd", t1_s_ynd, t1_e_ynd), ("C_zigzag", t1_s_zz, t1_e_zz),
                ("ADD_pf3ph_Dyn", t1_s_pf, t1_e_pf)]
    return branches, must


# ============================ Track 2 ============================
def t2_volts(rng, k):
    return sorted([U(rng, 0.4, 132) for _ in range(k)], reverse=True)


def t2_base(rng):
    vk0 = U(rng, 4, 14)
    return dict(vk0=vk0, vkr0=U(rng, 0.3, 1.5), si0=U(rng, 0.5, 0.95), RN=U(rng, 0, 200), XN=U(rng, 0, 50),
                ssc=U(rng, 100, 10000), x0x=U(rng, 0.5, 3), x0x2=U(rng, 0.5, 3))


def t2_s_radial(rng):
    p = t2_base(rng)
    p["hv"], p["lv"] = t2_volts(rng, 2)
    p["sn"] = U(rng, 0.1, 2.5)
    return p


def t2_s_zz(rng):
    p = t2_base(rng)
    p["hv"], p["mv"], p["lv"] = t2_volts(rng, 3)
    p["sn"] = U(rng, 0.1, 2.5)
    p["m_sn"], p["m_vk0"], p["m_vkr0"], p["m_si0"] = U(rng, 10, 100), U(rng, 4, 14), U(rng, 0.3, 1.5), U(rng, 0.5, 0.95)
    return p


def t2_s_par(rng):
    p = t2_base(rng)
    p["hv"], p["mv"] = t2_volts(rng, 2)
    p["sn"] = U(rng, 10, 100)
    p["Z"] = [complex(U(rng, 0, 200), U(rng, 0, 50)) for _ in range(3)]
    return p


def _t2_grid(net, b, p):
    return grid(net, b, p["ssc"], rx=0.1, x0x=p["x0x"], r0x0=0.1)


def _t2_I2(P, net, t, bus, u, case, inv, tag):
    c = cfac(u, case)
    ik = []
    for rn in (10.0, 100.0, 1000.0):
        setzn(net, t, rn, 0.0)
        r = sc(net, case, inv)[bus]
        ik.append(r["ikss_ka"])
        P.close("I2.ikss_formula[%s,RN=%g]" % (tag, rn), r["ikss_ka"], ik1(u, c, z1(r), z0(r)), 1e-9, 0.0, "max")
    P.cond("I2.strict_decrease_RN[%s]" % tag, ik[0] > ik[1] > ik[2], obs=ik)
    asym = c * u / (SQ3 * 1000.0)
    P.cond("I2.asymptote_RN1000<1%%[%s]" % tag, abs(ik[2] / asym - 1) < 0.01, obs=ik[2], exp=asym, resid=abs(ik[2] / asym - 1))


def _t2_I1_radial(P, r0, rf, p, tag):
    P.close("I1.dR0=3RN[%s]" % tag, rf["rk0_ohm"] - r0["rk0_ohm"], 3 * p["RN"], 0.0, 1e-9)
    P.close("I1.dX0=3XN[%s]" % tag, rf["xk0_ohm"] - r0["xk0_ohm"], 3 * p["XN"], 0.0, 1e-9)
    P.close("I1.R1_indep[%s]" % tag, rf["rk_ohm"], r0["rk_ohm"], 0.0, 1e-9)
    P.close("I1.X1_indep[%s]" % tag, rf["xk_ohm"], r0["xk_ohm"], 0.0, 1e-9)


def t2_make_e_dyn(case, inv):
    def ev(p, P):
        net = pp.create_empty_network(sn_mva=1.0)
        bh = pp.create_bus(net, vn_kv=p["hv"])
        bl = pp.create_bus(net, vn_kv=p["lv"])
        _t2_grid(net, bh, p)
        t = trafo(net, bh, bl, p["sn"], p["hv"], p["lv"], p["vk0"], p["vkr0"], "Dyn", si0=p["si0"])
        r0 = sc(net, case, inv)[bl]
        setzn(net, t, p["RN"], p["XN"])
        rf = sc(net, case, inv)[bl]
        _t2_I1_radial(P, r0, rf, p, "Dyn_lv")
        c = cfac(p["lv"], case)
        for key, r in (("0", r0), ("ZN", rf)):
            P.close("I2.ikss_formula[Dyn_lv,%s]" % key, r["ikss_ka"], ik1(p["lv"], c, z1(r), z0(r)), 1e-9, 0.0, "max")
        _t2_I2(P, net, t, bl, p["lv"], case, inv, "Dyn_lv")
    return ev


def _t2_zz_net(p, vg):
    net = pp.create_empty_network(sn_mva=1.0)
    bh = pp.create_bus(net, vn_kv=p["hv"])
    bm = pp.create_bus(net, vn_kv=p["mv"])
    ba = pp.create_bus(net, vn_kv=p["lv"])
    g = _t2_grid(net, bh, p)
    trafo(net, bh, bm, p["m_sn"], p["hv"], p["mv"], p["m_vk0"], p["m_vkr0"], "YNd", si0=p["m_si0"])
    te = trafo(net, bm, ba, p["sn"], p["mv"], p["lv"], p["vk0"], p["vkr0"], vg, si0=p["si0"])
    return net, bh, bm, ba, g, te


def t2_make_e_zz(vg, case, inv):
    def ev(p, P):
        net, bh, bm, ba, g, te = _t2_zz_net(p, vg)
        r0 = sc(net, case, inv)
        setzn(net, te, p["RN"], p["XN"])
        rf = sc(net, case, inv)
        if vg == "Zd":
            for key, r in (("0", r0), ("ZN", rf)):
                P.close("I3b.Zd_mv_ikss=0[%s]" % key, r[bm]["ikss_ka"], 0.0, 0.0, 1e-9)
            return
        _t2_I1_radial(P, r0[bm], rf[bm], p, vg + "_mv")
        c = cfac(p["mv"], case)
        for key, r in (("0", r0), ("ZN", rf)):
            P.close("I2.ikss_formula[%s_mv,%s]" % (vg, key), r[bm]["ikss_ka"], ik1(p["mv"], c, z1(r[bm]), z0(r[bm])), 1e-9, 0.0, "max")
        if vg == "ZNd":
            for key, r in (("0", r0), ("ZN", rf)):
                P.close("I3c.ZNd_lv_ikss=0[%s]" % key, r[ba]["ikss_ka"], 0.0)
        if vg == "ZNyn":
            net.ext_grid.loc[g, ["x0x_max", "x0x_min"]] = [p["x0x2"], p["x0x2"]]
            rx2 = sc(net, case, inv)
            net.ext_grid.loc[g, ["x0x_max", "x0x_min"]] = [p["x0x"], p["x0x"]]
            P.close("I3a.ZNyn_lv_R0_indep_x0x", rx2[ba]["rk0_ohm"], rf[ba]["rk0_ohm"], 0.0, 1e-9)
            P.close("I3a.ZNyn_lv_X0_indep_x0x", rx2[ba]["xk0_ohm"], rf[ba]["xk0_ohm"], 0.0, 1e-9)
            P.close("I3a.ZNyn_lv_R0_indep_hvRN", rf[ba]["rk0_ohm"], r0[ba]["rk0_ohm"], 0.0, 1e-9)
            P.close("I3a.ZNyn_lv_X0_indep_hvRN", rf[ba]["xk0_ohm"], r0[ba]["xk0_ohm"], 0.0, 1e-9)
        _t2_I2(P, net, te, bm, p["mv"], case, inv, vg + "_mv")
    return ev


def t2_make_e_par(case, inv):
    def ev(p, P):
        net = pp.create_empty_network(sn_mva=1.0)
        bh = pp.create_bus(net, vn_kv=p["hv"])
        bm = pp.create_bus(net, vn_kv=p["mv"])
        _t2_grid(net, bh, p)
        t = trafo(net, bh, bm, p["sn"], p["hv"], p["mv"], p["vk0"], p["vkr0"], "YNd", si0=p["si0"])
        zs = [0j] + [complex(z) for z in p["Z"]]
        rs = []
        for z in zs:
            setzn(net, t, z.real, z.imag)
            rs.append(sc(net, case, inv)[bh])
        Y = [1 / z0(r) for r in rs]
        K = (Y[1] - Y[0]) / (3 * (zs[0] - zs[1]))
        bsum, prod = Y[0] + Y[1], Y[0] * Y[1] - K
        disc = cmath.sqrt(bsum ** 2 - 4 * prod)
        best = None
        for u in ((bsum + disc) / 2, (bsum - disc) / 2):
            zb = 1 / (Y[0] - u) - 3 * zs[0]
            e3 = abs(1 / (Y[2] - u) - 3 * zs[2] - zb)
            if best is None or e3 < best[0]:
                best = (e3, u, zb)
        _e3, u, zb = best
        z4 = 1 / (u + 1 / (zb + 3 * zs[3]))
        P.record("I1.mobius_fit_point3_residual", value=_e3, Z0_p=1 / u, Z0_branch=zb)
        P.cclose("I1.mobius_predict_4th", z0(rs[3]), z4, 1e-9, 0.0, "max", comp=False)
        c = cfac(p["hv"], case)
        for i, r in enumerate(rs):
            P.close("I2.ikss_formula[YNd_hv,%d]" % i, r["ikss_ka"], ik1(p["hv"], c, z1(r), z0(r)), 1e-9, 0.0, "max")
            if i:
                P.close("I1.R1_indep[%d]" % i, r["rk_ohm"], rs[0]["rk_ohm"], 0.0, 1e-9)
                P.close("I1.X1_indep[%d]" % i, r["xk_ohm"], rs[0]["xk_ohm"], 0.0, 1e-9)
    return ev


def track2():
    branches = []
    for case in ("max", "min"):
        for inv in (True, False):
            tag = "[%s,inv=%s]" % (case, inv)
            branches += [("Dyn_radial" + tag, t2_s_radial, t2_make_e_dyn(case, inv)),
                         ("ZNd_radial" + tag, t2_s_zz, t2_make_e_zz("ZNd", case, inv)),
                         ("ZNyn_radial" + tag, t2_s_zz, t2_make_e_zz("ZNyn", case, inv)),
                         ("YNd_parallel" + tag, t2_s_par, t2_make_e_par(case, inv)),
                         ("Zd" + tag, t2_s_zz, t2_make_e_zz("Zd", case, inv))]
    return branches, []


# ============================ Track 3 ============================
CMP_Z = ("rk0_ohm", "xk0_ohm", "rk_ohm", "xk_ohm")
RT3, AZ3, AI3 = 1e-7, 1e-8, 1e-8
BASES = (10.0, 100.0, 0.1, 123.0)


def _vk_pair(rng):
    vk = LU(rng, 2, 20)
    while True:
        vkr = LU(rng, 0.1, 2)
        if vkr < vk:
            return vk, vkr


def t3_sample(rng):
    p = {}
    for pre in ("", "m_", "b_"):
        vk, vkr = _vk_pair(rng)
        vk0, vkr0 = _vk_pair(rng)
        p.update({pre + "sn": LU(rng, 0.1, 100), pre + "vk": vk, pre + "vkr": vkr, pre + "vk0": vk0, pre + "vkr0": vkr0,
                  pre + "mag0": LU(rng, 100, 1e4), pre + "mag0_rx": U(rng, 0, 1), pre + "s": U(rng, 0.01, 0.99)})
    ssc = LU(rng, 100, 1e4)
    p.update(dict(hv=LU(rng, 66, 220), mv=LU(rng, 3.3, 33), lv=LU(rng, 0.23, 1.0),
                  RN=LU(rng, 0.01, 100), XN=LU(rng, 0.01, 100), RN2=LU(rng, 0.01, 100), XN2=LU(rng, 0.01, 100),
                  ssc=ssc, ssc_min=LU(rng, 100, ssc), rx=U(rng, 0.05, 0.3), x0x=U(rng, 0.5, 3), r0x0=U(rng, 0.05, 0.3),
                  rx_min=U(rng, 0.05, 0.3), x0x_min=U(rng, 0.5, 3), r0x0_min=U(rng, 0.05, 0.3),
                  x0x_alt=U(rng, 0.5, 3), r0x0_alt=U(rng, 0.05, 0.3), mag0_alt=LU(rng, 100, 1e4),
                  case=choice(rng, ["max", "min"]), lvtol=choice(rng, [6, 10]),
                  ratio_h=U(rng, 0.9, 1.1), ratio_l=U(rng, 0.9, 1.1), tap_side=choice(rng, ["hv", "lv"]),
                  tap_pos=choice(rng, [-2, 0, 2]), lvbus_f=U(rng, 0.9, 1.1),
                  par=int(rng.integers(1, 5)), path=choice(rng, ["Dyn", "ZNd"]), s_oos=choice(rng, [0.01, 0.99]),
                  perm_seed=int(rng.integers(1 << 30)), tgroup=choice(rng, ["Dyn", "YNd", "Yyn", "YNy", "YNyn"]),
                  phase_load=U(rng, 0.01, 0.2), pf=U(rng, 0.9, 1.0), zn_ratio_pf=LU(rng, 0.1, 1.0),
                  zn_ang_pf=U(rng, 0, math.pi / 2), f_hz=choice(rng, [50.0, 60.0]), base_pf=choice(rng, list(BASES))))
    return p


def _g3(net, bus, p, alt=False):
    return grid(net, bus, p["ssc"], p["ssc_min"], p["rx"], p["x0x_alt"] if alt else p["x0x"],
                p["r0x0_alt"] if alt else p["r0x0"], p["rx_min"], p["x0x_alt"] if alt else p["x0x_min"],
                p["r0x0_alt"] if alt else p["r0x0_min"])


def _zg3(p, u, case, alt=False):
    if case == "max":
        return zgrid0(u, p["ssc"], p["rx"], p["x0x_alt"] if alt else p["x0x"], p["r0x0_alt"] if alt else p["r0x0"], cfac(u, "max"))
    return zgrid0(u, p["ssc_min"], p["rx_min"], p["x0x_alt"] if alt else p["x0x_min"],
                  p["r0x0_alt"] if alt else p["r0x0_min"], cfac(u, "min"))


def _tr3(net, hv, lv, p, vg, vhv, vlv, pre="", rn=0.0, xn=0.0, **kw):
    return trafo(net, hv, lv, p[pre + "sn"], vhv, vlv, p[pre + "vk"], p[pre + "vkr"], vg, vk0=p[pre + "vk0"],
                 vkr0=p[pre + "vkr0"], mag0=p[pre + "mag0"], mag0_rx=p[pre + "mag0_rx"], si0=p[pre + "s"], rn=rn, xn=xn, **kw)


def _kt3(p, case, u_lv, pre=""):
    return kt_iec(p[pre + "vk"], p[pre + "vkr"], cfac(u_lv, "max", p["lvtol"])) if case == "max" else 1.0


def _cmp(P, tag, ra, rb, buses, cols=None, rtol=RT3):
    for b in buses:
        for col in (cols or ra[b].keys()):
            if col not in ra[b] or col not in rb[b]:
                continue
            at = AI3 if col.endswith("_ka") else AZ3
            P.close("%s[b%d,%s]" % (tag, b, col), ra[b][col], rb[b][col], rtol, at)


def run3(net, P, p, tag, fault="1ph", **kw):
    """Track 3: every solve runs both inverse_y=True/False and compares them directly (I6), and checks that the currents are finite and non-negative and Re Z0≥−ε. Returns the result for inverse_y=True."""
    rT = sc(net, p["case"], True, fault, p["lvtol"], **kw)
    rF = sc(net, p["case"], False, fault, p["lvtol"], **kw)
    _cmp(P, "I6.solver_equal[%s]" % tag, rF, rT, rT.keys())
    for b, r in rT.items():
        ik = r["ikss_ka"]
        if math.isnan(ik):
            P.fail("I6.ikss_nan[%s,b%d]" % (tag, b), "ikss NaN")
        else:
            P.cond("I6.ikss_finite_nonneg[%s,b%d]" % (tag, b), math.isfinite(ik) and ik >= 0, obs=ik)
        if "rk0_ohm" in r and math.isfinite(r["rk0_ohm"]):
            P.cond("I6.ReZ0>=-eps[%s,b%d]" % (tag, b), r["rk0_ohm"] >= -AZ3, obs=r["rk0_ohm"])
    return rT


def bases3(net, P, p, ref, tag):
    for sb in BASES:
        net.sn_mva = sb
        r = sc(net, p["case"], True, "1ph", p["lvtol"])
        _cmp(P, "I4.base_cov[%s,Sb=%g]" % (tag, sb), r, ref, ref.keys())
    net.sn_mva = 1.0


def i3_3(net, P, p, t, bus, u, ref0, refz, tag, sweep=True):
    """I3: recompute with the current formula; 3ph/2ph unchanged before and after the neutral-point change; R, X sweeps on an isolated loop are monotonic and show the large-R asymptotics."""
    c = cfac(u, p["case"], p["lvtol"])
    for key, r in (("0", ref0[bus]), ("ZN", refz[bus])):
        if math.isfinite(r["rk0_ohm"]):
            P.close("I3.ikss_recon[%s,%s]" % (tag, key), r["ikss_ka"], ik1(u, c, z1(r), z0(r)), RT3, AI3)
    rn0, xn0 = float(net.trafo.at[t, "rn_ohm"]), float(net.trafo.at[t, "xn_ohm"])
    for fault in ("3ph", "2ph"):
        setzn(net, t, 0.0, 0.0)
        a = sc(net, p["case"], True, fault, p["lvtol"])
        setzn(net, t, p["RN"], p["XN"])
        b = sc(net, p["case"], True, fault, p["lvtol"])
        _cmp(P, "I3.%s_unchanged_by_ZN[%s]" % (fault, tag), b, a, a.keys(), ("ikss_ka", "rk_ohm", "xk_ohm"))
    if sweep:
        z0t = abs(z0(ref0[bus]))
        ikr = []
        for ratio in (0.1, 1.0, 10.0, 100.0):
            rn = ratio * z0t / 3
            setzn(net, t, rn, 0.0)
            ikr.append(sc(net, p["case"], True, "1ph", p["lvtol"])[bus]["ikss_ka"])
        P.cond("I3.monotone_R[%s]" % tag, all(ikr[i + 1] <= ikr[i] * (1 + RT3) + AI3 for i in range(3)), obs=ikr)
        rbig = 100.0 * z0t / 3
        rem = abs(2 * z1(ref0[bus]) + z0(ref0[bus]))
        asym = c * u / (SQ3 * rbig)
        P.cond("I3.asymptote_R[%s]" % tag, abs(ikr[3] / asym - 1) <= rem / (3 * rbig) + RT3, obs=ikr[3], exp=asym,
               resid=abs(ikr[3] / asym - 1))
        ikx = []
        for ratio in (0.1, 1.0, 10.0):
            setzn(net, t, 0.0, ratio * z0t / 3)
            ikx.append(sc(net, p["case"], True, "1ph", p["lvtol"])[bus]["ikss_ka"])
        P.cond("I3.monotone_X[%s]" % tag, all(ikx[i + 1] <= ikx[i] * (1 + RT3) + AI3 for i in range(2)), obs=ikx)
    setzn(net, t, rn0, xn0)


def _inc(P, tag, ra, rb, zn):
    P.cclose("I1.dZ0=3ZN[%s]" % tag, z0(rb) - z0(ra), 3 * zn, RT3, AZ3)


# ---- fixture 1: Dyn ----
def t3_e_dyn(p, P):
    net = pp.create_empty_network(sn_mva=1.0)
    bm = pp.create_bus(net, vn_kv=p["mv"])
    bl = pp.create_bus(net, vn_kv=p["lv"])
    _g3(net, bm, p)
    t = _tr3(net, bm, bl, p, "Dyn", p["mv"], p["lv"])
    r = {}
    for key, (rn, xn) in (("0", (0, 0)), ("R", (p["RN"], 0)), ("X", (0, p["XN"])), ("RX", (p["RN"], p["XN"]))):
        setzn(net, t, rn, xn)
        r[key] = run3(net, P, p, "Dyn," + key)
    _inc(P, "Dyn,pureR", r["0"][bl], r["R"][bl], complex(p["RN"], 0))
    _inc(P, "Dyn,pureX", r["0"][bl], r["X"][bl], complex(0, p["XN"]))
    _inc(P, "Dyn,RX", r["0"][bl], r["RX"][bl], complex(p["RN"], p["XN"]))
    P.cclose("I3.Z1_indep[Dyn]", z1(r["RX"][bl]), z1(r["0"][bl]), RT3, AZ3)
    bases3(net, P, p, r["RX"], "Dyn")
    i3_3(net, P, p, t, bl, p["lv"], r["0"], r["RX"], "Dyn")


# ---- fixture 2 ----
def _fx2(p, vg, alt=False, lv_scale=1.0):
    net = pp.create_empty_network(sn_mva=1.0)
    bh = pp.create_bus(net, vn_kv=p["hv"])
    bm = pp.create_bus(net, vn_kv=p["mv"])
    ba = pp.create_bus(net, vn_kv=p["lv"] * lv_scale)
    g = _g3(net, bh, p, alt)
    _tr3(net, bh, bm, p, "YNd", p["hv"], p["mv"], pre="m_")
    kw = {}
    if vg.lower() in ("znyn", "znd", "zny", "zd"):
        kw = dict(tap_side=p["tap_side"], tap_neutral=0, tap_min=-2, tap_max=2, tap_step_percent=2.5,
                  tap_pos=p["tap_pos"], tap_changer_type="Ratio")
        te = _tr3(net, bm, ba, p, vg, p["mv"] * p["ratio_h"], p["lv"] * p["ratio_l"], **kw)
    else:
        te = _tr3(net, bm, ba, p, vg, p["mv"], p["lv"])
    return net, bh, bm, ba, g, te


def t3_e_ynd(p, P):
    net, bh, bm, ba, g, te = _fx2(p, "YNd")
    r0 = run3(net, P, p, "YNd,0")
    setzn(net, te, p["RN"], p["XN"])
    rz = run3(net, P, p, "YNd,RX")
    _inc(P, "YNd_oneport_mv", r0[bm], rz[bm], complex(p["RN"], p["XN"]))
    P.cclose("I3.Z1_indep[YNd]", z1(rz[bm]), z1(r0[bm]), RT3, AZ3)
    bases3(net, P, p, rz, "YNd")
    i3_3(net, P, p, te, bm, p["mv"], r0, rz, "YNd")


def _urh_url(p):
    urh = p["mv"] * p["ratio_h"] * (1 + 0.025 * p["tap_pos"] if p["tap_side"] == "hv" else 1.0)
    url = p["lv"] * p["ratio_l"] * (1 + 0.025 * p["tap_pos"] if p["tap_side"] == "lv" else 1.0)
    return urh, url


def t3_make_e_zz(vg):
    def ev(p, P):
        net, bh, bm, ba, g, te = _fx2(p, vg)
        r0 = run3(net, P, p, vg + ",0")
        setzn(net, te, p["RN"], p["XN"])
        rz = run3(net, P, p, vg + ",RX")
        zn = complex(p["RN"], p["XN"])
        urh, url = _urh_url(p)
        kt = _kt3(p, p["case"], p["lv"])
        zlh = kt * zleak(p["vk0"], p["vkr0"], urh, p["sn"])
        for key, r, zz in (("0", r0, 0j), ("RX", rz, zn)):
            P.cclose("I2.Z0H=s*Zl,H+3ZN[%s,%s]" % (vg, key), z0(r[bm]), p["s"] * zlh + 3 * zz, RT3, AZ3)
        if vg == "ZNyn":
            zll_raw = zleak(p["vk0"], p["vkr0"], url, p["sn"])
            exp_l = zmag(p["mag0"], p["mag0_rx"], abs(zll_raw)) + (1 - p["s"]) * kt * zll_raw
            for key, r in (("0", r0), ("RX", rz)):
                P.cclose("I2.Z0L=Zm+(1-s)Zl,L[%s]" % key, z0(r[ba]), exp_l, RT3, AZ3)
            P.cclose("I2.ZNyn_lv_indep_hv_neutral", z0(rz[ba]), z0(r0[ba]), RT3, AZ3)
        else:
            for key, r in (("0", r0), ("RX", rz)):
                P.cond("I2.lv_ikss<1e-6[%s,%s]" % (vg, key), r[ba]["ikss_ka"] < 1e-6, obs=r[ba]["ikss_ka"], exp=0.0)
        # the zero-sequence strength of the remote grid changes: no cross-side change
        net2, *_x, te2 = _fx2(p, vg, alt=True)
        setzn(net2, te2, p["RN"], p["XN"])
        ra = sc(net2, p["case"], True, "1ph", p["lvtol"])
        P.cclose("I2.mv_Z0_indep_remote_grid0[%s]" % vg, z0(ra[bm]), z0(rz[bm]), RT3, AZ3)
        if vg == "ZNyn":
            P.cclose("I2.lv_Z0_indep_remote_grid0", z0(ra[ba]), z0(rz[ba]), RT3, AZ3)
        # change m: the high-voltage-side Z0 is unchanged
        net.trafo.at[te, "mag0_percent"] = p["mag0_alt"]
        rm = sc(net, p["case"], True, "1ph", p["lvtol"])
        net.trafo.at[te, "mag0_percent"] = p["mag0"]
        P.cclose("I2.hv_Z0_indep_m[%s]" % vg, z0(rm[bm]), z0(rz[bm]), RT3, AZ3)
        if vg == "ZNd":  # change only the rated voltage of the auxiliary low-voltage bus, winding ratings unchanged: the physical high-voltage-side Z0 is unchanged
            net3, *_y, te3 = _fx2(p, vg, lv_scale=p["lvbus_f"])
            setzn(net3, te3, p["RN"], p["XN"])
            rl = sc(net3, p["case"], True, "1ph", p["lvtol"])
            P.cclose("I4.ZNd_hv_Z0_indep_dummy_lv_vn", z0(rl[bm]), z0(rz[bm]), RT3, AZ3)
        bases3(net, P, p, rz, vg)
        i3_3(net, P, p, te, bm, p["mv"], r0, rz, vg)
    return ev


# ---- one-port Yyn / YNy, Yzn, YNyn, ungrounded groups ----
def t3_make_e_pi(vg):
    def ev(p, P):
        net = pp.create_empty_network(sn_mva=1.0)
        if vg == "Yyn":
            uh, ul = p["mv"], p["lv"]
        else:  # YNy: high-voltage side grounded, the grid is on the ungrounded low-voltage side
            uh, ul = p["hv"], p["mv"]
        bh = pp.create_bus(net, vn_kv=uh)
        bl = pp.create_bus(net, vn_kv=ul)
        if vg == "Yyn":
            _g3(net, bh, p)
            fb, uf = bl, ul
        else:
            _g3(net, bl, p)
            fb, uf = bh, uh
        t = _tr3(net, bh, bl, p, vg, uh, ul)
        kt = _kt3(p, p["case"], ul)
        zl = zleak(p["vk0"], p["vkr0"], uf, p["sn"])
        zm = zmag(p["mag0"], p["mag0_rx"], abs(zl))
        r0 = run3(net, P, p, vg + ",0")
        setzn(net, t, p["RN"], p["XN"])
        rz = run3(net, P, p, vg + ",RX")
        for key, r, zz in (("0", r0, 0j), ("RX", rz, complex(p["RN"], p["XN"]))):
            P.cclose("I1.pi_oneport_Zm+Zl+3ZN[%s,%s]" % (vg, key), z0(r[fb]), zm + kt * zl + 3 * zz, RT3, AZ3)
        _inc(P, vg, r0[fb], rz[fb], complex(p["RN"], p["XN"]))
        bases3(net, P, p, rz, vg)
        i3_3(net, P, p, t, fb, uf, r0, rz, vg)
    return ev


def t3_e_yzn(p, P):
    net = pp.create_empty_network(sn_mva=1.0)
    bh = pp.create_bus(net, vn_kv=p["mv"])
    bl = pp.create_bus(net, vn_kv=p["lv"])
    _g3(net, bh, p)
    t = _tr3(net, bh, bl, p, "Yzn", p["mv"], p["lv"])
    r0 = run3(net, P, p, "Yzn,0")
    setzn(net, t, p["RN"], p["XN"])
    rz = run3(net, P, p, "Yzn,RX")
    zn = complex(p["RN"], p["XN"])
    _inc(P, "Yzn_lv", r0[bl], rz[bl], zn)
    kt = _kt3(p, p["case"], p["lv"])
    base_exp = (1 - p["s"]) * kt * zleak(p["vk0"], p["vkr0"], p["lv"], p["sn"])
    d = z0(r0[bl]) - base_exp
    P.record("I1.Yzn_baseline_leakage_calibration(separate)", obs=z0(r0[bl]), exp=base_exp, resid=abs(d),
             within_tol=bool(abs(d.real) <= AZ3 + RT3 * abs(base_exp.real) and abs(d.imag) <= AZ3 + RT3 * abs(base_exp.imag)))
    bases3(net, P, p, rz, "Yzn")
    i3_3(net, P, p, t, bl, p["lv"], r0, rz, "Yzn")


def _ynyn_oracle(p, case, uh, ul, alt):
    kt = _kt3(p, case, ul)
    zlh = zleak(p["vk0"], p["vkr0"], uh, p["sn"])
    A = p["s"] * kt * zlh + 3 * complex(p["RN"], p["XN"])
    B = (1 - p["s"]) * kt * zlh
    M = zmag(p["mag0"], p["mag0_rx"], abs(zlh))
    zg = _zg3(p, uh, case, alt)
    zh = par(zg, A + M)
    zl_ = (B + par(M, A + zg)) * (ul / uh) ** 2
    a, b, m0 = 1 / A, 1 / B, 1 / M
    Y = np.array([[a, 0], [0, b]]) - np.array([[a * a, a * b], [a * b, b * b]]) / (a + b + m0)
    return zh, zl_, Y


def t3_e_ynyn(p, P):
    uh, ul = p["mv"], p["lv"]
    for alt in (False, True):
        net = pp.create_empty_network(sn_mva=1.0)
        bh = pp.create_bus(net, vn_kv=uh)
        bl = pp.create_bus(net, vn_kv=ul)
        _g3(net, bh, p, alt)
        t = _tr3(net, bh, bl, p, "YNyn", uh, ul, rn=p["RN"], xn=p["XN"])
        r = run3(net, P, p, "YNyn,alt=%s" % alt)
        zh, zl_, Y = _ynyn_oracle(p, p["case"], uh, ul, alt)
        P.cclose("I1.YNyn_T_drivingpoint_HV[alt=%s]" % alt, z0(r[bh]), zh, RT3, AZ3)
        P.cclose("I1.YNyn_T_drivingpoint_LV[alt=%s]" % alt, z0(r[bl]), zl_, RT3, AZ3)
        ev = np.linalg.eigvalsh((Y + Y.conj().T) / 2)
        P.cond("I6.oracle_Y0_hermitian_psd[alt=%s]" % alt, ev.min() >= -1e-12 * max(1.0, abs(ev).max()), obs=ev.tolist())
        if not alt:
            setzn(net, t, 0.0, 0.0)
            r0 = sc(net, p["case"], True, "1ph", p["lvtol"])
            setzn(net, t, p["RN"], p["XN"])
            bases3(net, P, p, r, "YNyn")
            i3_3(net, P, p, t, bl, ul, r0, r, "YNyn", sweep=False)


def t3_make_e_ungr(vg):
    def ev(p, P):
        net = pp.create_empty_network(sn_mva=1.0)
        bh = pp.create_bus(net, vn_kv=p["mv"])
        bl = pp.create_bus(net, vn_kv=p["lv"])
        _g3(net, bh, p)
        t = _tr3(net, bh, bl, p, vg, p["mv"], p["lv"])
        r0 = run3(net, P, p, vg + ",0")
        setzn(net, t, p["RN"], p["XN"])
        rz = run3(net, P, p, vg + ",RX")
        _cmp(P, "I1.ungrounded_neutral_no_effect[%s]" % vg, rz, r0, r0.keys())
    return ev


# ---- I5: parallel, in/out of service, ordering, the four create interfaces ----
def _i5_build(p, rows, path, perm=None):
    """rows: [(parallel, rn, xn, in_service, s override or None, parameter prefix)]; all attached to the same pair of buses."""
    net = pp.create_empty_network(sn_mva=1.0)
    if path == "Dyn":
        bm = pp.create_bus(net, vn_kv=p["mv"])
        bl = pp.create_bus(net, vn_kv=p["lv"])
        _g3(net, bm, p)
        hv, lv, uh, ul, vg, fb = bm, bl, p["mv"], p["lv"], "Dyn", bl
    else:
        bh = pp.create_bus(net, vn_kv=p["hv"])
        bm = pp.create_bus(net, vn_kv=p["mv"])
        bl = pp.create_bus(net, vn_kv=p["lv"])
        _g3(net, bh, p)
        _tr3(net, bh, bm, p, "YNd", p["hv"], p["mv"], pre="m_", index=1000)
        hv, lv, uh, ul, vg, fb = bm, bl, p["mv"], p["lv"], "ZNd", bm
    order = list(range(len(rows)))
    idx = [None] * len(rows)
    if perm is not None:
        order = perm
        idx = [int(10 + 7 * k) for k in range(len(rows))]
    for k in order:
        par_, rn, xn, ins, s_ov, pre, vgo = rows[k]
        kw = dict(parallel=par_, in_service=ins)
        if idx[k] is not None:
            kw["index"] = idx[k]
        q = dict(p)
        if s_ov is not None:
            q[pre + "s"] = s_ov
        _tr3(net, hv, lv, q, vgo or vg, uh, ul, pre=pre, rn=rn, xn=xn, **kw)
    return net, fb


def _std_data(p, vg, uh, ul):
    return dict(sn_mva=p["sn"], vn_hv_kv=uh, vn_lv_kv=ul, vk_percent=p["vk"], vkr_percent=p["vkr"], pfe_kw=0.0,
                i0_percent=0.0, shift_degree=SHIFT.get(vg.lower(), 0.0), vector_group=vg, vk0_percent=p["vk0"],
                vkr0_percent=p["vkr0"], mag0_percent=p["mag0"], mag0_rx=p["mag0_rx"], si0_hv_partial=p["s"])


def _i5_islands(p, path, api, neut):
    """Three identical islands (each with an external grid), neutral points by neut=[(rn,xn)|None]×3; api selects the create interface. None = parameter not passed."""
    net = pp.create_empty_network(sn_mva=1.0)
    uh, ul = (p["mv"], p["lv"]) if path == "Dyn" else (p["mv"], p["lv"])
    vg = "Dyn" if path == "Dyn" else "ZNd"
    hvs, lvs, fbs = [], [], []
    for _k in range(3):
        if path == "Dyn":
            bm = pp.create_bus(net, vn_kv=uh)
            bl = pp.create_bus(net, vn_kv=ul)
            _g3(net, bm, p)
            hvs.append(bm)
            lvs.append(bl)
            fbs.append(bl)
        else:
            bh = pp.create_bus(net, vn_kv=p["hv"])
            bm = pp.create_bus(net, vn_kv=uh)
            bl = pp.create_bus(net, vn_kv=ul)
            _g3(net, bh, p)
            _tr3(net, bh, bm, p, "YNd", p["hv"], p["mv"], pre="m_")
            hvs.append(bm)
            lvs.append(bl)
            fbs.append(bm)
    std = "dq1_std_%s" % vg
    if api in ("std", "stds"):
        pp.create_std_type(net, _std_data(p, vg, uh, ul), std, element="trafo")
    if api == "params":
        for k in range(3):
            kw = {} if neut[k] is None else dict(rn=neut[k][0], xn=neut[k][1])
            if neut[k] is None:
                pp.create_transformer_from_parameters(net, hvs[k], lvs[k], **{kk: v for kk, v in _std_data(p, vg, uh, ul).items()})
            else:
                _tr3(net, hvs[k], lvs[k], p, vg, uh, ul, **kw)
    elif api == "params_vec":
        d = _std_data(p, vg, uh, ul)
        kw = {}
        if any(x is not None for x in neut):
            kw = dict(rn_ohm=[NAN if x is None else x[0] for x in neut], xn_ohm=[NAN if x is None else x[1] for x in neut])
        pp.create_transformers_from_parameters(net, hvs, lvs, **d, **kw)
    elif api == "params_bulk_scalar":
        d = _std_data(p, vg, uh, ul)
        pp.create_transformers_from_parameters(net, hvs, lvs, **d, rn_ohm=neut[0][0], xn_ohm=neut[0][1])
    elif api == "std":
        for k in range(3):
            kw = {} if neut[k] is None else dict(rn_ohm=neut[k][0], xn_ohm=neut[k][1])
            pp.create_transformer(net, hvs[k], lvs[k], std, **kw)
    elif api == "stds":
        kw = {}
        if any(x is not None for x in neut):
            kw = dict(rn_ohm=[NAN if x is None else x[0] for x in neut], xn_ohm=[NAN if x is None else x[1] for x in neut])
        pp.create_transformers(net, hvs, lvs, std, **kw)
    return net, fbs


def t3_e_i5(p, P):
    path, n = p["path"], p["par"]
    zn = (p["RN"], p["XN"])
    tagp = "%s,p=%d" % (path, n)
    # (a) parallel=p with neutral-point Z_N/p  ==  p explicit rows each carrying Z_N
    ng, fb = _i5_build(p, [(n, zn[0] / n, zn[1] / n, True, None, "", None)], path)
    rg = run3(ng, P, p, "I5.grouped/p," + tagp)
    ne, _ = _i5_build(p, [(1, zn[0], zn[1], True, None, "", None)] * n, path)
    re_ = run3(ne, P, p, "I5.explicit," + tagp)
    _cmp(P, "I5.grouped(ZN/p)==explicit[%s]" % tagp, rg, re_, [fb])
    # (b) parallel=p with the same Z_N: record only (difference under the independent-neutral-point convention)
    ng2, _ = _i5_build(p, [(n, zn[0], zn[1], True, None, "", None)], path)
    rg2 = sc(ng2, p["case"], True, "1ph", p["lvtol"])
    P.record("I5.grouped(ZN)_vs_explicit(record)", grouped=z0(rg2[fb]), explicit=z0(re_[fb]), diff=z0(rg2[fb]) - z0(re_[fb]))
    # (c) zero neutral point
    ng0, _ = _i5_build(p, [(n, 0.0, 0.0, True, None, "", None)], path)
    ne0, _ = _i5_build(p, [(1, 0.0, 0.0, True, None, "", None)] * n, path)
    _cmp(P, "I5.zero_neutral_grouped==explicit[%s]" % tagp, sc(ng0, p["case"], True, "1ph", p["lvtol"]),
         sc(ne0, p["case"], True, "1ph", p["lvtol"]), [fb])
    # (d) one unit out of service (s near the boundary) == that unit deleted (the plan says "duplicate", so p=1 is not run)
    if n >= 2:
        _i5_oos(P, p, path, zn, n, tagp, fb)
    _i5_rest(P, p, path, zn, n, tagp, fb, re_)


def _i5_oos(P, p, path, zn, n, tagp, fb):
    rows = [(1, zn[0], zn[1], False, p["s_oos"], "", None)] + [(1, zn[0], zn[1], True, None, "", None)] * (n - 1)
    no, _ = _i5_build(p, rows, path)
    ro = run3(no, P, p, "I5.oos," + tagp)
    nd, _ = _i5_build(p, [(1, zn[0], zn[1], True, None, "", None)] * (n - 1), path)
    rd = sc(nd, p["case"], True, "1ph", p["lvtol"])
    _cmp(P, "I5.out_of_service==deleted[%s,s=%g]" % (tagp, p["s_oos"]), ro, rd, [fb])


def _i5_rest(P, p, path, zn, n, tagp, fb, re_):
    # (e) shuffle row / index order
    perm = list(np.random.default_rng(p["perm_seed"]).permutation(n))
    nsh, _ = _i5_build(p, [(1, zn[0] * (1 + 0.1 * k), zn[1], True, None, "", None) for k in range(n)], path, perm=perm)
    nun, _ = _i5_build(p, [(1, zn[0] * (1 + 0.1 * k), zn[1], True, None, "", None) for k in range(n)], path)
    _cmp(P, "I5.order_independent[%s]" % tagp, sc(nsh, p["case"], True, "1ph", p["lvtol"]),
         sc(nun, p["case"], True, "1ph", p["lvtol"]), [fb])
    # (f) two units with different parameters (and mixed vector groups) on the same bus: parallel shunt branches to ground accumulate
    vg2 = "Dyn" if path == "Dyn" else "ZNyn"
    nm, _ = _i5_build(p, [(1, zn[0], zn[1], True, None, "", None), (1, p["RN2"], p["XN2"], True, None, "b_", vg2)], path)
    rm = sc(nm, p["case"], True, "1ph", p["lvtol"])
    if path == "Dyn":
        za = _kt3(p, p["case"], p["lv"]) * zleak(p["vk0"], p["vkr0"], p["lv"], p["sn"]) + 3 * complex(*zn)
        zb = _kt3(p, p["case"], p["lv"], "b_") * zleak(p["b_vk0"], p["b_vkr0"], p["lv"], p["b_sn"]) + 3 * complex(p["RN2"], p["XN2"])
    else:
        za = p["s"] * _kt3(p, p["case"], p["lv"]) * zleak(p["vk0"], p["vkr0"], p["mv"], p["sn"]) + 3 * complex(*zn)
        zb = p["b_s"] * _kt3(p, p["case"], p["lv"], "b_") * zleak(p["b_vk0"], p["b_vkr0"], p["mv"], p["b_sn"]) + 3 * complex(p["RN2"], p["XN2"])
    P.cclose("I5.shunts_accumulate[%s+%s]" % (path, vg2), z0(rm[fb]), par(za, zb), RT3, AZ3)
    # (g) the four create interfaces, scalar batch, vector, mixed zero/NaN/non-zero, omitted columns
    neut = [zn, None, (0.0, 0.0)]
    res = {}
    for api in ("params", "params_vec", "std", "stds"):
        nn, fbs = _i5_islands(p, path, api, neut)
        res[api] = (sc(nn, p["case"], True, "1ph", p["lvtol"]), fbs)
    r_ref, fb_ref = res["params"]
    for api in ("params_vec", "std", "stds"):
        rr, fbs = res[api]
        for k in range(3):
            _cmp(P, "I5.creator_equal[%s vs params,island%d]" % (api, k), {0: rr[fbs[k]]}, {0: r_ref[fb_ref[k]]}, [0])
    _cmp(P, "I5.NaN_neutral==zero_neutral", {0: r_ref[fb_ref[1]]}, {0: r_ref[fb_ref[2]]}, [0])
    nb_, fbs = _i5_islands(p, path, "params_bulk_scalar", [zn] * 3)
    rb_ = sc(nb_, p["case"], True, "1ph", p["lvtol"])
    for k in range(3):
        _cmp(P, "I5.bulk_scalar==scalar[island%d]" % k, {0: rb_[fbs[k]]}, {0: r_ref[fb_ref[0]]}, [0])
    nom, fbs = _i5_islands(p, path, "params_vec", [None, None, None])
    rom = sc(nom, p["case"], True, "1ph", p["lvtol"])
    _cmp(P, "I5.omitted_columns==zero", {0: rom[fbs[0]]}, {0: r_ref[fb_ref[2]]}, [0])


# ---- PSU (t1 + G1 data of iec_60909_4()) ----
def _psu_net(p, vg="YNd", rn=None, xn=0.0, psu=True, with_plain=False, zn_plain=(0.0, 0.0)):
    net = pp.create_empty_network(sn_mva=34)
    b4 = pp.create_bus(net, vn_kv=110.0)
    hg1 = pp.create_bus(net, vn_kv=20.0)
    kw = dict(xn_ohm=xn)
    if rn is not None:
        kw["rn_ohm"] = rn
    t1 = pp.create_transformer_from_parameters(
        net, b4, hg1, sn_mva=150, pfe_kw=0, i0_percent=0, vn_hv_kv=115.0, vn_lv_kv=21, vk_percent=16,
        vkr_percent=0.5, pt_percent=12, oltc=True, vk0_percent=15.2, vkr0_percent=0.5, vector_group=vg,
        mag0_percent=10000, mag0_rx=0, si0_hv_partial=0.5, power_station_unit=psu, **kw)
    pp.create_gen(net, hg1, p_mw=0.85 * 150, vn_kv=21, xdss_pu=0.14, rdss_ohm=0.002, cos_phi=0.85, sn_mva=150,
                  pg_percent=0, power_station_trafo=t1, slack=True)
    extra = None
    if with_plain:
        bp = pp.create_bus(net, vn_kv=110.0)
        bq = pp.create_bus(net, vn_kv=20.0)
        _g3(net, bq, p)
        trafo(net, bp, bq, 40.0, 110.0, 20.0, 12.0, 0.4, "YNd", si0=0.5, rn=zn_plain[0], xn=zn_plain[1],
              power_station_unit=False)
        extra = bp
    return net, b4, hg1, extra


def t3_e_psu(p, P):
    RN, XN = p["RN"], p["XN"]
    cfgs = {"omit,0": (None, 0.0), "0,0": (0.0, 0.0), "R,0": (RN, 0.0), "omit,X": (None, XN), "R,X": (RN, XN)}
    r = {}
    for k, (rn, xn) in cfgs.items():
        net, b4, hg1, _ = _psu_net(p, rn=rn, xn=xn)
        r[k] = run3(net, P, p, "PSU," + k)
    _cmp(P, "PSU.rn_omitted==rn_zero", r["0,0"], r["omit,0"], r["omit,0"].keys())
    for k in ("R,0", "omit,X", "R,X"):
        rn, xn = cfgs[k]
        _inc(P, "PSU_YNd_hv,%s" % k, r["omit,0"][b4], r[k][b4], complex(rn or 0.0, xn))
    for vg in ("ZNyn", "ZNd", "ZNy", "Zd"):
        try:
            na, b4, hg1, _ = _psu_net(p, vg=vg, rn=0.0, xn=0.0)
            ra = sc(na, p["case"], True, "1ph", p["lvtol"])
            nb, _b, _h, _ = _psu_net(p, vg=vg, rn=RN, xn=XN)
            rb = sc(nb, p["case"], True, "1ph", p["lvtol"])
        except Exception as e:  # plan: if unsupported it should be rejected explicitly -- rejection counts as pass, recorded
            P.record("PSU.zigzag_rejected", vg=vg, err="%s: %s" % (type(e).__name__, str(e)[:120]))
            continue
        _inc(P, "PSU_%s_hv_counted_once" % vg, ra[b4], rb[b4], complex(RN, XN))
        if vg == "ZNyn":
            P.cclose("PSU.ZNyn_lv_indep_hv_neutral", z0(rb[hg1]), z0(ra[hg1]), RT3, AZ3)
    # PSU flags mixed within the same vector group
    na, b4, _h, bp = _psu_net(p, with_plain=True)
    ra = sc(na, p["case"], True, "1ph", p["lvtol"])
    nb, _b4, _h2, _bp = _psu_net(p, rn=RN, xn=XN, with_plain=True, zn_plain=(p["RN2"], p["XN2"]))
    rb = sc(nb, p["case"], True, "1ph", p["lvtol"])
    _inc(P, "PSU_mixed_flags,psu_unit", ra[b4], rb[b4], complex(RN, XN))
    _inc(P, "PSU_mixed_flags,plain_unit", ra[bp], rb[bp], complex(p["RN2"], p["XN2"]))


# ---- T-model path: runpp(trafo_model) → calc_sc(use_pre_fault_voltage=True) ----
def t3_e_tmodel(p, P):
    vg = p["tgroup"]
    if vg in ("Dyn", "Yyn", "YNyn"):
        uh, ul, gside, fside = p["mv"], p["lv"], "h", "l"
    else:  # YNd / YNy: the grid is on the ungrounded (low-voltage) side
        uh, ul, gside, fside = p["hv"], p["mv"], "l", "h"
    zn = complex(p["RN"], p["XN"])
    kt = _kt3(p, p["case"], ul)
    uf = ul if fside == "l" else uh
    zl = zleak(p["vk0"], p["vkr0"], uf, p["sn"])
    M = zmag(p["mag0"], p["mag0_rx"], abs(zl))
    s = p["s"]
    out = {}
    for model in ("t", "pi"):
        for key, zz in (("0", 0j), ("RX", zn)):
            net = pp.create_empty_network(sn_mva=1.0)
            bh = pp.create_bus(net, vn_kv=uh)
            bl = pp.create_bus(net, vn_kv=ul)
            _g3(net, bh if gside == "h" else bl, p)
            _tr3(net, bh, bl, p, vg, uh, ul, rn=zz.real, xn=zz.imag)
            fb = bl if fside == "l" else bh
            pp.runpp(net, trafo_model=model)
            calc_sc(net, fault="1ph", case=p["case"], bus=fb, use_pre_fault_voltage=True, lv_tol_percent=p["lvtol"])
            out[(model, key)] = z0({c: float(net.res_bus_sc.at[fb, c]) for c in net.res_bus_sc.columns})
            A, B = s * kt * zl, (1 - s) * kt * zl
            if vg == "YNyn":
                exp = _ynyn_oracle(dict(p, RN=zz.real, XN=zz.imag), p["case"], uh, ul, False)[1]
            elif model == "pi":
                exp = kt * zl + 3 * zz + (M if vg in ("Yyn", "YNy") else 0)
            elif vg == "Dyn":
                exp = B + 3 * zz + par(A, M)
            elif vg == "YNd":
                exp = A + 3 * zz + par(B, M)
            elif vg == "Yyn":
                exp = B + 3 * zz + M
            else:  # YNy
                exp = A + 3 * zz + M
            P.cclose("Tpath.Z0_vs_%s_oracle[%s,%s]" % (model, vg, key), out[(model, key)], exp, RT3, AZ3)
        if vg != "YNyn":
            P.cclose("Tpath.dZ0=3ZN[%s,%s]" % (model, vg), out[(model, "RX")] - out[(model, "0")], 3 * zn, RT3, AZ3)


# ---- pf_3ph ----
def _pf_oracle(p, vg, zn):
    """Low-voltage-side driving-point zero-sequence impedance (Ω, low-voltage side) of the physical T circuit (no IEC correction)."""
    uh, ul = p["mv"], p["lv"]
    s = p["s"]
    if vg == "Dyn":
        zl = zleak(p["vk0"], p["vkr0"], ul, p["sn"])
        M = zmag(p["mag0"], p["mag0_rx"], abs(zl))
        return (1 - s) * zl + 3 * zn + par(s * zl, M)
    if vg == "Yzn":
        zl = zleak(p["vk0"], p["vkr0"], ul, p["sn"])
        return (1 - s) * zl + 3 * zn
    zlh = zleak(p["vk0"], p["vkr0"], uh, p["sn"])  # YNyn: neutral point on the high-voltage arm (the convention chosen in the implementation), the external-grid zero sequence is the terminal
    A = s * zlh + 3 * zn
    B = (1 - s) * zlh
    M = zmag(p["mag0"], p["mag0_rx"], abs(zlh))
    zg = zgrid0(uh, p["ssc"], p["rx"], p["x0x"], p["r0x0"], 1.1)
    return (B + par(M, A + zg)) * (ul / uh) ** 2


def _pf_net(p, vg, zn, loads, sb=1.0):
    net = pp.create_empty_network(sn_mva=sb, f_hz=p["f_hz"])
    bh = pp.create_bus(net, vn_kv=p["mv"])
    bl = pp.create_bus(net, vn_kv=p["lv"])
    _g3(net, bh, p)
    _tr3(net, bh, bl, p, vg, p["mv"], p["lv"], rn=zn.real, xn=zn.imag)
    pp.create_asymmetric_load(net, bl, **loads)
    pp.runpp_3ph(net)
    if not net.converged:
        raise RuntimeError("runpp_3ph not converged")
    return net, bl


def _phasors(net, bus, u_kv):
    r = net.res_bus_3ph
    vb = u_kv / SQ3
    return [r.at[bus, "vm_%s_pu" % ph] * cmath.exp(1j * math.radians(r.at[bus, "va_%s_degree" % ph])) for ph in "abc"], vb


def t3_make_e_pf(vg):
    def ev(p, P):
        zl_ref = abs(zleak(p["vk0"], p["vkr0"], p["mv"] if vg == "YNyn" else p["lv"], p["sn"]))
        zmag_n = p["zn_ratio_pf"] * zl_ref / 3
        zn = complex(zmag_n * math.cos(p["zn_ang_pf"]), zmag_n * math.sin(p["zn_ang_pf"]))
        S = p["phase_load"] * p["sn"]
        pl, ql = S * p["pf"], S * math.sqrt(max(0.0, 1 - p["pf"] ** 2))
        zexp = _pf_oracle(p, vg, zn)
        va_ref = None
        for ph in "abc":
            loads = {"p_%s_mw" % ph: pl, "q_%s_mvar" % ph: ql}
            net, bl = _pf_net(p, vg, zn, loads)
            v, vb = _phasors(net, bl, p["lv"])
            i_load = [0j, 0j, 0j]
            k = "abc".index(ph)
            i_load[k] = (complex(pl, ql) / (v[k] * vb)).conjugate()  # kA (MVA/kV)
            v0 = sum(v) / 3.0
            i0 = sum(i_load) / 3.0
            v0_pred = -zexp * i0 / vb
            d = abs(v0 - v0_pred)
            P.cond("PF.V0=-Z0phys*I0[%s,phase %s]" % (vg, ph), d <= 1e-6 + 1e-5 * abs(v0), obs=v0, exp=v0_pred, resid=d)
            if ph == "a":
                va_ref = v
        net2, bl2 = _pf_net(p, vg, zn, {"p_a_mw": pl, "q_a_mvar": ql}, sb=p["base_pf"])
        v2, _ = _phasors(net2, bl2, p["lv"])
        for k in range(3):
            P.cond("PF.base_cov[%s,Sb=%g,ph%d]" % (vg, p["base_pf"], k), abs(v2[k] - va_ref[k]) <= 1e-6 + 1e-5 * abs(va_ref[k]),
                   obs=v2[k], exp=va_ref[k], resid=abs(v2[k] - va_ref[k]))
        bal = {"p_a_mw": pl / 3, "q_a_mvar": ql / 3, "p_b_mw": pl / 3, "q_b_mvar": ql / 3, "p_c_mw": pl / 3, "q_c_mvar": ql / 3}
        nz, blz = _pf_net(p, vg, zn, bal)
        n0, bl0 = _pf_net(p, vg, 0j, bal)
        vz, _ = _phasors(nz, blz, p["lv"])
        v0_, _ = _phasors(n0, bl0, p["lv"])
        for k in range(3):
            P.cond("PF.balanced_neutral_indep[%s,ph%d]" % (vg, k), abs(vz[k] - v0_[k]) <= 1e-6 + 1e-5 * abs(v0_[k]),
                   obs=vz[k], exp=v0_[k], resid=abs(vz[k] - v0_[k]))
    return ev


# ---- API conventions (mandatory points only) ----
def t3_e_api(p, P):
    net, bh, bm, ba, g, te = _fx2(p, "ZNyn")
    setzn(net, te, p["RN"], p["XN"])
    ref = sc(net, p["case"], True, "1ph", p["lvtol"])
    for sp in ("znyn", "ZNYN", "Znyn", "zNyN"):
        net.trafo.at[te, "vector_group"] = sp
        _cmp(P, "API.spelling_equal[%s]" % sp, sc(net, p["case"], True, "1ph", p["lvtol"]), ref, ref.keys())
    net.trafo.at[te, "vector_group"] = "ZNyn11"
    P.expect_raise("API.ZNyn11_explicit_validation_error", lambda: sc(net, p["case"], True, "1ph", p["lvtol"]), (ValueError,))
    net.trafo.at[te, "vector_group"] = "ZN"
    try:
        rzn = sc(net, p["case"], True, "1ph", p["lvtol"])
        P.record("API.ZN_attempt(noncoverage record)", outcome="accepted", z0_mv=z0(rzn[bm]))
    except Exception as e:
        P.record("API.ZN_attempt(noncoverage record)", outcome="rejected", err="%s: %s" % (type(e).__name__, str(e)[:120]))
    for vg in ("ZNyn", "ZNd", "ZNy", "Zd"):
        def f(vg=vg):
            _pf_net(p, vg, 0j, {"p_a_mw": 0.01 * p["sn"], "q_a_mvar": 0.0})
        P.expect_raise("API.pf3ph_rejects_hv_zigzag[%s]" % vg, f, (NotImplementedError, ValueError))


def t3_default(kind):
    """Deterministic fixture parameters given by the plan (fixture 1/2); those not given take the PR test values."""
    p = dict(sn=0.63, vk=6.0, vkr=1.0, vk0=6.0, vkr0=1.0, mag0=100.0, mag0_rx=0.0, s=0.9,
             m_sn=40.0, m_vk=12.0, m_vkr=0.4, m_vk0=12.0, m_vkr0=0.4, m_mag0=100.0, m_mag0_rx=0.0, m_s=0.9,
             b_sn=0.4, b_vk=4.0, b_vkr=1.0, b_vk0=4.0, b_vkr0=1.0, b_mag0=100.0, b_mag0_rx=0.0, b_s=0.9,
             hv=110.0, mv=20.0, lv=0.4, RN=0.01, XN=0.02, RN2=0.03, XN2=0.01, ssc=1000.0, ssc_min=1000.0, rx=0.1,
             x0x=1.0, r0x0=0.1, rx_min=0.1, x0x_min=1.0, r0x0_min=0.1, x0x_alt=2.0, r0x0_alt=0.2, mag0_alt=1000.0,
             case="max", lvtol=10, ratio_h=1.0, ratio_l=1.0, tap_side="hv", tap_pos=0, lvbus_f=1.0, par=2,
             path="Dyn", s_oos=0.01, perm_seed=1, tgroup="Dyn", phase_load=0.1, pf=0.95, zn_ratio_pf=0.5,
             zn_ang_pf=0.5, f_hz=50.0, base_pf=10.0)
    if kind == "fx2":
        p.update(sn=0.4, vk=4.0, vkr=1.0, vk0=4.0, vkr0=1.0, ssc=5000.0, ssc_min=5000.0)
    return p


def track3():
    branches = [("Dyn", t3_sample, t3_e_dyn), ("YNd_oneport", t3_sample, t3_e_ynd),
                ("Yyn", t3_sample, t3_make_e_pi("Yyn")), ("YNy", t3_sample, t3_make_e_pi("YNy")),
                ("YNyn", t3_sample, t3_e_ynyn), ("Yzn", t3_sample, t3_e_yzn)]
    branches += [(vg, t3_sample, t3_make_e_ungr(vg)) for vg in ("Yy", "Yd", "Dy", "Dd")]
    branches += [(vg, t3_sample, t3_make_e_zz(vg)) for vg in ("ZNyn", "ZNd", "ZNy", "Zd")]
    branches += [("I5_parallel_service_api", t3_sample, t3_e_i5), ("PSU", t3_sample, t3_e_psu),
                 ("Tmodel_path", t3_sample, t3_e_tmodel)]
    branches += [("PF3ph_" + vg, t3_sample, t3_make_e_pf(vg)) for vg in ("Dyn", "YNyn", "Yzn")]
    must = []
    s_list = (0.01, 0.1, 0.5, 0.9, 0.99)
    r_list = (0.1, 1.0, 10.0, 100.0, 1.0)
    for bname, _s, _e in branches:
        kind = "fx2" if bname in ("YNd_oneport", "ZNyn", "ZNd", "ZNy", "Zd") else "fx1"
        for s, ratio in zip(s_list, r_list):
            q = t3_default(kind)
            z0t = abs(zleak(q["vk0"], q["vkr0"], q["lv"] if kind == "fx1" else q["mv"], q["sn"]))
            zmag_ = ratio * z0t / 3
            q.update(s=s, RN=zmag_ / math.sqrt(2), XN=zmag_ / math.sqrt(2))
            must.append((bname, q))
    return branches, must


# ============================ Entry point ============================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--track", type=int, choices=(1, 2, 3), required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    import json
    import pandapower
    versions = dict(pandapower=pandapower.__version__, numpy=np.__version__,
                    merge_commit="0c6a1f0a65273af91dd1d5b854498cd89d457939",
                    # Chinese text kept verbatim: this string is written to the archived outputs/pp3016_t*.json; meaning: "env/pp3016: uv pip install pandapower from git at 0c6a1f0a (env/setup_pca_py.sh)"
                    install="env/pp3016：uv pip install pandapower@git+...@0c6a1f0a（env/setup_pca_py.sh）")
    branches, must = {1: track1, 2: track2, 3: track3}[a.track]()
    run(PR, PR_IDX, a.track, a.n, a.out, branches, must=must, versions=versions)
    if a.track == 3:  # API convention: only deterministic mandatory points, not part of the N-way split; results are merged into the same JSON
        api = run(PR, PR_IDX, 3, 0, a.out + ".api", [("API_contract", t3_sample, t3_e_api)],
                  must=[("API_contract", t3_default("fx2"))], versions=versions)
        with open(a.out) as f:
            d = json.load(f)
        d["api_contract"] = dict(branches=api["branches"], first_rejects=api["first_rejects"],
                                 records_sample=api["records_sample"])
        d["total_rejects"] += api["total_rejects"]
        d["total_fails"] += api["total_fails"]
        d["total_points"] += api["total_points"]
        d["first_rejects"] = (d["first_rejects"] + api["first_rejects"])[:20]
        with open(a.out + ".tmp", "w") as f:
            json.dump(d, f, ensure_ascii=False, indent=1)
        os.replace(a.out + ".tmp", a.out)
        os.remove(a.out + ".api")


if __name__ == "__main__":
    main()
