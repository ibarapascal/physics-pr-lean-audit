"""PyBaMM#5889 (Chayambuka2022 NVPF positive-electrode diffusivity looked up by solid-phase concentration) unified ΔQ1 spot check: following each of the three tracks' own check plans, sample uniformly on the native pybamm of the merge commit.
Usage: python spotcheck.py --track {1,2,3} --n N --out PATH
Environment: compute node A env/py2 (pybamm 26.10.1.0; Chayambuka2022.py, D_p.csv, D_n.csv, lithium_ion_parameters.py, settings.py,
interpolant.py are byte-identical to merge commit eda03ed, see impl_log.md). The plan mapping and implementation judgments are in impl_log.md.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spotlib as SL  # noqa: E402

import pybamm  # noqa: E402
from pybamm.input.parameters.sodium_ion import Chayambuka2022 as CHMOD  # noqa: E402

PR = "pybamm-team/PyBaMM#5889"
PR_IDX = 12
SET = "Chayambuka2022"
DP = "Positive particle diffusivity [m2.s-1]"
DN = "Negative particle diffusivity [m2.s-1]"
CK = "Maximum concentration in positive electrode [mol.m-3]"
EK = "Initial concentration in electrolyte [mol.m-3]"
C_DEF = 15320.0
CN_DEF = 14540.0
E_DEF = 1000.0

DATA = os.path.join(os.path.dirname(pybamm.__file__), "input", "parameters", "sodium_ion", "data")
XP, YP = np.loadtxt(os.path.join(DATA, "D_p.csv"), delimiter=",", skiprows=1).T
XN, YN = np.loadtxt(os.path.join(DATA, "D_n.csv"), delimiter=",", skiprows=1).T
TAU = pybamm.settings.tolerances["D__c_s"]

_PV = {}
_EVALS = [0]


def pv(C=C_DEF, E=E_DEF):
    """Build a brand-new ParameterValues for each (C,E) configuration; after building, call update only once and never modify it again (required by Track 3)."""
    key = (float(C), float(E))
    if key not in _PV:
        p = pybamm.ParameterValues(SET)
        p.update({CK: float(C), EK: float(E)})
        _PV[key] = p
    return _PV[key]


def _maybe_clear():
    _EVALS[0] += 1
    if _EVALS[0] % 2000 == 0:
        for p in _PV.values():
            proc = getattr(p, "_processor", None)
            if proc is not None and hasattr(proc, "clear_cache"):
                proc.clear_cache()


def _f(x):
    return float(np.asarray(x, dtype=float).reshape(-1)[0])


def Dp(theta, T=298.15, C=C_DEF, E=E_DEF):
    _maybe_clear()
    p = pv(C, E)
    return _f(p.evaluate(p[DP](pybamm.Scalar(float(theta)), pybamm.Scalar(float(T)))))


def Dn(theta, T=298.15):
    _maybe_clear()
    p = pv()
    return _f(p.evaluate(p[DN](pybamm.Scalar(float(theta)), pybamm.Scalar(float(T)))))


def versions():
    import scipy
    import casadi
    return dict(pybamm=pybamm.__version__, numpy=np.__version__, scipy=scipy.__version__, casadi=casadi.__version__,
                merge_commit="eda03ed6d1003b186cc2d70ad8b3240a7b6350f1",
                note="pybamm 26.10.1.0 PyPI; Chayambuka2022.py, D_p.csv, D_n.csv etc. are byte-identical to the merge commit",
                D__c_s=TAU)


# ---------------------------------------------------------------- Track 1
def t1_eval(inp):
    sto, T = inp["sto"], inp["T"]
    D = Dp(sto, T)
    cs = sto * C_DEF
    ch = [SL.close("t1_api_vs_numpy_interp_rtol1e-10", D, float(np.interp(cs, XP, YP)), rtol=1e-10, atol=0.0)]
    for c2 in (0.5 * C_DEF, 2.0 * C_DEF):
        ch.append(SL.close(f"t1_cmax_invariance_c={c2:g}_rtol1e-10", Dp(cs / c2, T, C=c2), D, rtol=1e-10, atol=0.0))
    ch.append(SL.cond("t1_magnitude_1.9e-17<=D<=2.6e-15", 1.9e-17 <= D <= 2.6e-15, observed=D,
                      expected=[1.9e-17, 2.6e-15]))
    if sto >= 0.3:
        ch.append(SL.cond("t1_magnitude_sto>=0.3_D<=1e-16", D <= 1e-16, observed=D, expected="<=1e-16"))
    return ch


def t1_sample(rng):
    return dict(sto=float(rng.uniform(0.01, 0.99)), T=float(rng.uniform(273, 323)))


def t1_macro_eval(inp):
    sim = pybamm.Simulation(pybamm.sodium_ion.BasicDFN(), parameter_values=pybamm.ParameterValues(SET),
                            experiment=pybamm.Experiment(["Discharge at 2C until 2.0 V"]))
    sol = sim.solve()
    Q = float(np.asarray(sol["Discharge capacity [A.h]"].entries).reshape(-1)[-1])
    return [SL.close("t1_macro_BasicDFN_2C_capacity_approx_1.18mAh", Q * 1e3, 1.18, rtol=0.0, atol=0.005)]


# ---------------------------------------------------------------- Track 2
def base_NVPF_diffusivity(sto, T):
    """Same-named function of parent commit c235b90 (the only difference: c_max is taken as the initial electrolyte concentration), used as the pre-fix control for I2."""
    name, (x, y) = CHMOD.D_p_data
    c_max = pybamm.Parameter("Initial concentration in electrolyte [mol.m-3]")
    return pybamm.Interpolant(x, y, sto * c_max, name)


class Capacity:
    def __init__(self):
        self.cache_exp = {}
        p = pybamm.ParameterValues(SET)
        self.Qnom = float(p["Nominal cell capacity [A.h]"])
        self.vcut = float(p["Lower voltage cut-off [V]"])
        p.update({"Current function [A]": "[input]"})
        self.sim = pybamm.Simulation(pybamm.sodium_ion.BasicDFN(), parameter_values=p)

    def fast(self, crate):
        sol = self.sim.solve([0, 5 * 3600.0 / crate], inputs={"Current function [A]": crate * self.Qnom})
        return float(np.asarray(sol["Discharge capacity [A.h]"].entries).reshape(-1)[-1])

    def experiment(self, crate_str, parent=False):
        key = (crate_str, parent)
        if key not in self.cache_exp:
            p = pybamm.ParameterValues(SET)
            if parent:
                p.update({DP: base_NVPF_diffusivity})
            sim = pybamm.Simulation(pybamm.sodium_ion.BasicDFN(), parameter_values=p,
                                    experiment=pybamm.Experiment([f"Discharge at {crate_str} until 2.0 V"]))
            sol = sim.solve()
            self.cache_exp[key] = float(np.asarray(sol["Discharge capacity [A.h]"].entries).reshape(-1)[-1])
        return self.cache_exp[key]


_CAP = None


def cap():
    global _CAP
    if _CAP is None:
        _CAP = Capacity()
        _CAP.anchor_fast = {r: _CAP.fast(r) for r in (1 / 50, 1.0, 2.0)}
    return _CAP


def t2_eval(inp):
    sto, crate = inp["sto"], inp["C_rate"]
    D = Dp(sto, 298.15)
    ch = [
        SL.close("t2_I1_api_vs_numpy_interp_rtol1e-10", D, float(np.interp(sto * C_DEF, XP, YP)), rtol=1e-10, atol=0.0),
        SL.cond("t2_I1_bounds_1.7e-17<=D<=2.51e-15", 1.7e-17 <= D <= 2.51e-15, observed=D, expected=[1.7e-17, 2.51e-15]),
        SL.cond("t2_I1_positive", D > 0, observed=D, expected=">0"),
        SL.close("t2_I1_negative_Dn_vs_numpy_interp_rtol1e-10", Dn(sto), float(np.interp(sto * CN_DEF, XN, YN)),
                 rtol=1e-10, atol=0.0),
    ]
    c = cap()
    Q = c.fast(crate)
    a = c.anchor_fast
    tolm = lambda q: SL.DEFAULT_ATOL + SL.DEFAULT_RTOL * abs(q)  # noqa: E731
    ch.append(SL.cond("t2_I2_Q(C/50)>=Q(C)", a[1 / 50] >= Q - tolm(Q), observed=Q, expected=f"<= {a[1 / 50]}"))
    if crate <= 1.0:
        ch.append(SL.cond("t2_I2_Q(C)>=Q(1C)_for_C<=1", Q >= a[1.0] - tolm(Q), observed=Q, expected=f">= {a[1.0]}"))
    else:
        ch.append(SL.cond("t2_I2_Q(C)<=Q(1C)_for_C>1", Q <= a[1.0] + tolm(Q), observed=Q, expected=f"<= {a[1.0]}"))
    if crate <= 2.0:
        ch.append(SL.cond("t2_I2_Q(C)>=Q(2C)_for_C<=2", Q >= a[2.0] - tolm(Q), observed=Q, expected=f">= {a[2.0]}"))
    else:
        ch.append(SL.cond("t2_I2_Q(C)<=Q(2C)_for_C>2", Q <= a[2.0] + tolm(Q), observed=Q, expected=f"<= {a[2.0]}"))
    return ch


def t2_sample(rng):
    return dict(sto=float(rng.uniform(0.01, 0.99)), C_rate=float(rng.uniform(1 / 50, 3.0)))


def t2_node_eval(inp):
    k = inp["k"]
    D = Dp(XP[k] / C_DEF, 298.15)
    return [SL.close("t2_I1_node_exact_rtol1e-10", D, YP[k], rtol=1e-10, atol=0.0)]


def t2_endpoint_eval(inp):
    D = Dp(inp["sto"], 298.15)
    return [SL.cond("t2_I1_endpoint_positive", D > 0, observed=D, expected=">0"),
            SL.cond("t2_I1_bounds_1.7e-17<=D<=2.51e-15", 1.7e-17 <= D <= 2.51e-15, observed=D,
                    expected=[1.7e-17, 2.51e-15])]


def t2_anchor_eval(inp):
    c = cap()
    q50 = c.experiment("C/50")
    q50p = c.experiment("C/50", parent=True)
    q1 = c.experiment("1C")
    q2 = c.experiment("2C")
    return [
        SL.close("t2_I2_C/50_merge_vs_parent_rel2pct", q50, q50p, rtol=0.02, atol=0.0),
        SL.close("t2_I2_2C_vs_PR_1.18mAh_rtol2pct", q2 * 1e3, 1.18, rtol=0.02, atol=0.0),
        SL.cond("t2_I2_monotone_Q(C/50)>=Q(1C)>=Q(2C)", q50 >= q1 >= q2, observed=[q50, q1, q2]),
    ]


def t2_tau_eval(inp):
    R = float(pybamm.ParameterValues(SET)["Positive particle radius [m]"])
    tau = R ** 2 / Dp(0.5, 298.15)
    return [SL.cond("t2_I3_tau_p_same_order_or_larger_than_1800s(>=180s)", tau >= 180.0, observed=tau,
                    expected=">=180 s (implementation judgment: same order of magnitude = not below 1/10)")]


# ---------------------------------------------------------------- Track 3
def L(c):
    """Independent reference: pick the two rows that bracket c and evaluate with the two-point formula; outside the table, extrapolate linearly from the end segment."""
    c = float(c)
    if c <= XP[0]:
        i = 0
    elif c >= XP[-1]:
        i = len(XP) - 2
    else:
        i = int(np.searchsorted(XP, c, side="right") - 1)
        i = min(max(i, 0), len(XP) - 2)
    lam = (c - XP[i]) / (XP[i + 1] - XP[i])
    return (1 - lam) * YP[i] + lam * YP[i + 1], i


_BDFN = None


def model_D(c, T=298.15):
    global _BDFN
    if _BDFN is None:
        _BDFN = pybamm.sodium_ion.BasicDFN()
    _maybe_clear()
    p = pv()
    return _f(p.evaluate(_BDFN.param.p.prim.D(pybamm.Scalar(float(c)), pybamm.Scalar(float(T)))))


C_ALT = (12256.0, 15320.0, 18384.0)
E_ALT = (500.0, 1000.0, 1500.0)


def t3_point_checks(c, include_model=True):
    theta = c / C_DEF
    ref, i = L(c)
    D = Dp(theta)
    ch = [SL.close("t3_I1_direct_vs_L(c)_rtol1e-10", D, ref, rtol=1e-10, atol=0.0)]
    if include_model:
        ch.append(SL.close("t3_I1_model_adapter_vs_L(c)_rtol1e-10", model_D(c), ref, rtol=1e-10, atol=0.0))
    for E in (500.0, 1500.0):
        ch.append(SL.close(f"t3_I2_E={E:g}_equals_default_rtol1e-10", Dp(theta, E=E), D, rtol=1e-10, atol=0.0))
    if XP[0] <= c <= 12256.0 * (1 - TAU):
        for C in C_ALT:
            for E in E_ALT:
                ch.append(SL.close(f"t3_I2_C={C:g}_E={E:g}_vs_L(c)_rtol1e-10", Dp(c / C, C=C, E=E), ref,
                                   rtol=1e-10, atol=0.0))
    lo, hi = min(YP[i], YP[i + 1]), max(YP[i], YP[i + 1])
    m = 1e-12 * hi
    ch.append(SL.cond("t3_I3_positive", D > 0, observed=D, expected=">0"))
    ch.append(SL.cond("t3_I3_within_bracketing_ordinates", lo - m <= D <= hi + m, observed=D, expected=[lo, hi]))
    for g in (0.0, 1e9, -1e9):
        ch.append(SL.cond(f"t3_I3_flux_sign_g={g:g}", -D * g * g <= 0, observed=-D * g * g, expected="<=0"))
    for T in (273.15, 323.15):
        ch.append(SL.close(f"t3_I4_T={T}_equals_298.15_rtol1e-12", Dp(theta, T=T), D, rtol=1e-12, atol=0.0))
    return ch


def t3_eval(inp):
    return t3_point_checks(inp["c"])


def t3_sample(rng):
    return dict(c=float(rng.uniform(XP[0], XP[-1])))


def t3_named_eval(inp):
    theta, expv = inp["theta"], inp["expected"]
    D = Dp(theta)
    ch = [SL.close("t3_I1_independent_value_rtol1e-10", D, expv, rtol=1e-10, atol=0.0)]
    ch += t3_point_checks(theta * C_DEF)
    return ch


def t3_adapter_boundary_eval(inp):
    r = inp["c_over_C"]
    c = r * C_DEF
    expv, _ = L(C_DEF * min(max(r, TAU), 1 - TAU))
    return [SL.close("t3_I1_model_adapter_clipped_vs_L_rtol1e-10", model_D(c), expv, rtol=1e-10, atol=0.0)]


def t3_direct_extrap_eval(inp):
    theta = inp["theta"]
    expv, _ = L(theta * C_DEF)
    ch = [SL.close("t3_I1_direct_endpoint_extension_rtol1e-10", Dp(theta), expv, rtol=1e-10, atol=0.0)]
    if "approx" in inp:
        ch.append(SL.close("t3_I1_direct_endpoint_approx_value_rtol1e-10", Dp(theta), inp["approx"],
                           rtol=1e-10, atol=0.0))
    return ch


def t3_knot_delta_eval(inp):
    i, sgn = inp["i"], inp["sign"]
    delta = 1e-6 * min(XP[i] - XP[i - 1], XP[i + 1] - XP[i])
    c = XP[i] + sgn * delta
    j = i if sgn > 0 else i - 1
    slope = (YP[j + 1] - YP[j]) / (XP[j + 1] - XP[j])
    one_sided = YP[j] + slope * (c - XP[j])
    D = Dp(c / C_DEF)
    return [SL.close("t3_I3_knot_delta_one_sided_rtol1e-10", D, one_sided, rtol=1e-10, atol=0.0),
            SL.cond("t3_I3_knot_delta_bound", abs(D - YP[i]) <= abs(slope) * delta + 1e-12 * YP[i],
                    observed=abs(D - YP[i]), expected=abs(slope) * delta + 1e-12 * YP[i])]


def t3_vector_eval(inp):
    cs = np.asarray(inp["cs"], float)
    theta = cs / C_DEF
    p = pv()
    vec = np.asarray(p.evaluate(p[DP](pybamm.Vector(theta), pybamm.Scalar(298.15))), float).reshape(-1)
    scal = np.array([Dp(t) for t in theta])
    return [SL.close("t3_I1_vector_vs_scalar_rtol1e-10", vec, scal, rtol=1e-10, atol=0.0)]


# ---------------------------------------------------------------- Assembly
def main():
    a = SL.parse_args()
    run = SL.Runner(PR, PR_IDX, a.track, a.n, a.out, versions())
    if a.track == 1:
        branches = [("NVPF_direct", t1_sample, t1_eval)]
        mand = [(f"grid_sto={k / 100:.2f}", dict(sto=k / 100, T=298.15), t1_eval) for k in range(1, 100)]
        mand.append(("macro_BasicDFN_2C", dict(case="BasicDFN default, 2C to 2.0 V"), t1_macro_eval))
    elif a.track == 2:
        branches = [("NVPF_direct+BasicDFN", t2_sample, t2_eval)]
        mand = [(f"node_k={k}", dict(k=k), t2_node_eval) for k in range(len(XP))]
        mand += [(f"endpoint_sto={s}", dict(sto=float(s)), t2_endpoint_eval) for s in (0.0, 1.0)]
        mand.append(("anchors_C/50_1C_2C_parent", dict(case="Experiment anchors"), t2_anchor_eval))
        mand.append(("tau_p_sto=0.5", dict(sto=0.5), t2_tau_eval))
    else:
        branches = [("NVPF_table_interior", t3_sample, t3_eval)]
        mand = [(f"knot_{k}", dict(c=float(XP[k])), t3_eval) for k in range(len(XP))]
        fr_cs = []
        for k in range(len(XP) - 1):
            for fr in (0.01, 0.25, 0.5, 0.75, 0.99):
                c = float(XP[k] + fr * (XP[k + 1] - XP[k]))
                fr_cs.append(c)
                mand.append((f"interval_{k}_frac_{fr}", dict(c=c), t3_eval))
        for th, ev in ((0.3, 5.6673280002e-17), (0.5, 8.3317028568e-17), (0.9, 3.5626719998e-17)):
            mand.append((f"theta_{th}", dict(theta=th, expected=ev), t3_named_eval))
        for r in (0.0, TAU / 2, TAU, 2 * TAU, 1 - 2 * TAU, 1 - TAU, 1 - TAU / 2, 1.0):
            mand.append((f"adapter_c/C={r!r}", dict(c_over_C=float(r)), t3_adapter_boundary_eval))
        mand.append(("direct_theta=0", dict(theta=0.0, approx=2.7385714278e-15), t3_direct_extrap_eval))
        mand.append(("direct_theta=x0/2C", dict(theta=float(XP[0] / (2 * C_DEF))), t3_direct_extrap_eval))
        mand.append(("direct_theta=(xlast/C+1)/2", dict(theta=float((XP[-1] / C_DEF + 1) / 2)), t3_direct_extrap_eval))
        mand.append(("direct_theta=1", dict(theta=1.0, approx=1.7387466637e-17), t3_direct_extrap_eval))
        for i in range(1, len(XP) - 1):
            for sgn in (+1, -1):
                mand.append((f"knot_{i}_delta_{sgn:+d}", dict(i=i, sign=sgn), t3_knot_delta_eval))
        # vector form: sampled points + nodes + interval-fraction points, one Vector call, compared element-wise with the scalar
        rng_v = np.random.default_rng(SL.seed_for(PR_IDX, 3))
        samp = [float(rng_v.uniform(XP[0], XP[-1])) for _ in range(a.n)]  # same seed and sequence as the main sampling
        mand.append(("vector_samples", dict(cs=samp), t3_vector_eval))
        mand.append(("vector_knots_intervals", dict(cs=[float(x) for x in XP] + fr_cs), t3_vector_eval))
    run.run(branches, mand, progress_every=100 if a.track == 2 else 1000)


if __name__ == "__main__":
    main()
