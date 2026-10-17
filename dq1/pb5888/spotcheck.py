"""PyBaMM#5888 (Ai2020 electrolyte diffusivity cm²/s→m²/s) unified ΔQ1 spot check: following each of the three tracks' own check plans, sample uniformly on the native pybamm of the merge commit.
Usage: python spotcheck.py --track {1,2,3} --n N --out PATH
Environment: compute node A env/py2 (pybamm 26.10.1.0; the relevant source files are byte-identical to merge commit 343ece2, see impl_log.md).
The plan mapping and all implementation judgments are in impl_log.md in the same directory.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spotlib as SL  # noqa: E402

import pybamm  # noqa: E402

PR = "pybamm-team/PyBaMM#5888"
PR_IDX = 11
DK = "Electrolyte diffusivity [m2.s-1]"
KK = "Electrolyte conductivity [S.m-1]"

AI = pybamm.ParameterValues("Ai2020")
XU = pybamm.ParameterValues("Xu2019")
LP = pybamm.LithiumIonParameters()
_EVALS = [0]


def _maybe_clear():
    # ParameterValues caches processed expressions by symbol internally; it would grow without bound over 10⁴ points, so it is cleared periodically (does not affect numerical values)
    _EVALS[0] += 1
    if _EVALS[0] % 2000 == 0:
        for p in (AI, XU):
            proc = getattr(p, "_processor", None)
            if proc is not None and hasattr(proc, "clear_cache"):
                proc.clear_cache()


def _f(x):
    return float(np.asarray(x, dtype=float).reshape(-1)[0])


def val(params, key, c, T):
    _maybe_clear()
    return _f(params.evaluate(params[key](pybamm.Scalar(float(c)), pybamm.Scalar(float(T)))))


def D_ai(c, T):
    return val(AI, DK, c, T)


def D_xu(c, T):
    return val(XU, DK, c, T)


def D_model(c, T):
    _maybe_clear()
    return _f(AI.evaluate(LP.D_e(pybamm.Scalar(float(c)), pybamm.Scalar(float(T)))))


def kappa_ai(c, T):
    return val(AI, KK, c, T)


def kappa_model(c, T):
    _maybe_clear()
    return _f(AI.evaluate(LP.kappa_e(pybamm.Scalar(float(c)), pybamm.Scalar(float(T)))))


def versions():
    import scipy
    import casadi
    return dict(pybamm=pybamm.__version__, numpy=np.__version__, scipy=scipy.__version__,
                casadi=casadi.__version__, merge_commit="343ece2062dbb39789dbfa977e0f2daf96e9e370",
                note="pybamm 26.10.1.0 PyPI; Ai2020.py/Xu2019.py/lithium_ion_parameters.py/settings.py are byte-identical to the merge commit",
                D_e__c_e=pybamm.settings.tolerances["D_e__c_e"],
                kappa_e__c_e=pybamm.settings.tolerances["kappa_e__c_e"])


# ---------------------------------------------------------------- Track 1
def t1_lit(c, T):
    return 1e-4 * 10 ** (-4.43 - 54 / (T - 229 - 5e-3 * c) - 0.22e-3 * c)


def t1_eval(inp):
    c, T = inp["c_e"], inp["T"]
    D = D_ai(c, T)
    ch = [
        SL.close("t1_api_vs_literature_rtol1e-12", D, t1_lit(c, T), rtol=1e-12, atol=0.0),
        SL.close("t1_ai_vs_xu2019_rtol1e-12", D, D_xu(c, T), rtol=1e-12, atol=0.0),
        SL.cond("t1_magnitude_1e-11<=D<=2e-9", 1e-11 <= D <= 2e-9, observed=D, expected=[1e-11, 2e-9]),
    ]
    dT = (D_ai(c, T + 0.01) - D_ai(c, T - 0.01)) / 0.02
    dc = (D_ai(c + 0.1, T) - D_ai(c - 0.1, T)) / 0.2
    ch.append(SL.cond("t1_fd_sign_dD_dT>0", dT > 0, observed=dT, expected=">0"))
    ch.append(SL.cond("t1_fd_sign_dD_dc<0", dc < 0, observed=dc, expected="<0"))
    return ch


def t1_sample(rng):
    return dict(c_e=float(rng.uniform(100, 3000)), T=float(rng.uniform(253, 333)))


def t1_dfn_eval(inp):
    model = pybamm.lithium_ion.DFN()
    sim = pybamm.Simulation(model, parameter_values=pybamm.ParameterValues("Ai2020"))
    sol = sim.solve([0, 3600])
    ce = np.asarray(sol["Electrolyte concentration [mol.m-3]"].entries)
    last = ce[:, -1] if ce.ndim == 2 else ce
    spread = float(np.max(last) - np.min(last))
    return [SL.cond("t1_optional_dfn_1C_final_ce_spread>=10", spread >= 10.0, observed=spread,
                    expected=">>1 mol/m3 (implementation judgment: >=10)")]


# ---------------------------------------------------------------- Track 2
def t2_formula(c, T):
    C = c / 1000.0  # mol/L
    return 1e-4 * 10 ** (-4.43 - 54 / (T - 229 - 5 * C) - 0.22 * C)


class DFNCase:
    """Track 2 I3: Ai2020 DFN, with the current taken as an input parameter so the built model can be reused."""

    def __init__(self):
        pv = pybamm.ParameterValues("Ai2020")
        self.Q = float(pv["Nominal cell capacity [A.h]"])
        self.eps = float(pv["Separator porosity"])
        self.b = float(pv["Separator Bruggeman coefficient (electrolyte)"])
        self.tplus = float(pv["Cation transference number"])
        self.F = float(pybamm.constants.F.value)
        try:
            self.T_amb = float(pv["Ambient temperature [K]"])
        except Exception:  # noqa: BLE001
            self.T_amb = 298.15
        pv.update({"Current function [A]": "[input]"})
        self.sim = pybamm.Simulation(pybamm.lithium_ion.DFN(), parameter_values=pv)

    def run(self, crate, use_cache=False):
        if use_cache and getattr(self, "_last", None) is not None and self._last[0] == crate:
            return self._last[1]
        out = self._run(crate)
        self._last = (crate, out)
        return out

    def _run(self, crate):
        I = crate * self.Q
        sol = self.sim.solve([0, 3600.0 / crate], inputs={"Current function [A]": I})
        t_end = float(sol.t[-1])
        c = np.asarray(sol["Separator electrolyte concentration [mol.m-3]"].entries)[:, -1]
        xs = np.asarray(sol["x_s [m]"].entries)
        xs = xs[:, -1] if xs.ndim == 2 else xs
        try:
            Ts = np.asarray(sol["Separator temperature [K]"].entries)[:, -1]
        except Exception:  # noqa: BLE001
            Ts = np.full_like(c, self.T_amb)
        i_cell = float(np.asarray(sol["Total current density [A.m-2]"].entries).reshape(-1)[-1])
        L = np.asarray(sol["Total lithium in electrolyte [mol]"].entries).reshape(-1)
        return dict(t_end=t_end, c=c, xs=xs, Ts=Ts, i=i_cell, L0=float(L[0]), L1=float(L[-1]))


_DFN = None


def dfn():
    global _DFN
    if _DFN is None:
        _DFN = DFNCase()
    return _DFN


def t2_eval(inp, use_cache=False):
    c, T, crate = inp["c_e"], inp["T"], inp["C_rate"]
    D = D_ai(c, T)
    ch = [
        SL.close("t2_I1_api_vs_formula_rtol1e-12", D, t2_formula(c, T), rtol=1e-12, atol=0.0),
        SL.cond("t2_I1_magnitude_1e-11<=D<=1e-9", 1e-11 <= D <= 1e-9, observed=D, expected=[1e-11, 1e-9]),
        SL.close("t2_I2_ai_vs_xu2019_rtol1e-12", D, D_xu(c, T), rtol=1e-12, atol=0.0),
    ]
    dT = (D_ai(c, T + 0.01) - D_ai(c, T - 0.01)) / 0.02
    dc = (D_ai(c + 0.1, T) - D_ai(c - 0.1, T)) / 0.2
    ch.append(SL.cond("t2_I2_fd_sign_dD_dT>0", dT > 0, observed=dT, expected=">0"))
    ch.append(SL.cond("t2_I2_fd_sign_dD_dc<0", dc < 0, observed=dc, expected="<0"))
    ch.extend(t2_dfn_checks(crate, use_cache=use_cache))
    return ch


def t2_replay(inp, bad):
    """Replay of a rejected point: if the rejected check involves the DFN (t2_I3_*), recompute the whole point including the DFN; if it involves only constitutive checks, recompute the constitutive part and reuse the DFN result just computed."""
    dfn_bad = any(c["name"].startswith("t2_I3") for c in bad)
    return t2_eval(inp, use_cache=not dfn_bad)


def t2_dfn_checks(crate, use_cache=False):
    d = dfn()
    r = d.run(crate, use_cache=use_cache)
    out = []
    if r["t_end"] >= 1000.0:
        dcdx = np.gradient(r["c"], r["xs"], edge_order=2)
        Deff = np.array([D_model(ci, Ti) for ci, Ti in zip(r["c"], r["Ts"])]) * d.eps ** d.b
        lhs = -Deff * dcdx
        rhs = (1 - d.tplus) * r["i"] / d.F
        out.append(SL.close("t2_I3_separator_flux_5pct", lhs, np.full_like(lhs, rhs), rtol=0.05, atol=0.0))
    else:
        out.append(SL.cond("t2_I3_quasi_steady_t_end>=1000s_applicable", True, observed=r["t_end"],
                           expected="applies only for t>=1000 s; recorded as pass when not applicable"))
    drift = abs(r["L1"] - r["L0"]) / abs(r["L0"])
    out.append(SL.cond("t2_I3_electrolyte_Li_drift<=1e-6", drift <= 1e-6, observed=drift, expected="<=1e-6"))
    return out


def t2_sample(rng):
    return dict(c_e=float(rng.uniform(0, 3000)), T=float(rng.uniform(263, 333)),
                C_rate=float(rng.uniform(0.1, 2.0)))


def t2_nominal_eval(inp):
    D = D_ai(1000.0, 298.15)
    return [SL.close("t2_I1_nominal_3.2227e-10_(5 sig)", D, 3.2227e-10, rtol=0.0, atol=0.00005e-10)]


# ---------------------------------------------------------------- Track 3
T3_C = [500, 750, 1000, 1250, 1500, 1750, 2000]
T3_T = [273.15, 278.15, 298.15, 318.15, 333.15]
T3_BC = [0, 5, 9.999, 10, 10.001, 20]
CMIN = 10.0


def t3_oracle_D(c, T):
    C = c / 1000.0
    Tg = 229 + 5 * C
    return float(np.exp(np.log(10) * (-4.43 - 54 / (T - Tg) - .22 * C)) * 1e-4)


def t3_oracle_kappa(c, T):
    P = ((-10.5 + 0.668e-3 * c + 0.494e-6 * c ** 2) + (0.074 - 1.78e-5 * c - 8.86e-10 * c ** 2) * T
         + (-6.96e-5 + 2.8e-8 * c) * T ** 2)
    return 1e-4 * c * P ** 2


def t3_dlnD(c, T):
    q = T - 229 - 0.005 * c
    return math.log(10) * 54 / q ** 2, -math.log(10) * (0.27 / q ** 2 + 0.00022)


def t3_eval(inp):
    c, T = inp["c_e"], inp["T"]
    D = D_ai(c, T)
    Dm = D_model(c, T)
    k = kappa_ai(c, T)
    km = kappa_model(c, T)
    ch = [
        SL.close("t3_I1_api_vs_oracle_rtol1e-12", D, t3_oracle_D(c, T), rtol=1e-12, atol=0.0),
        SL.close("t3_I2_ai_vs_xu2019_rtol1e-12", D, D_xu(c, T), rtol=1e-12, atol=0.0),
        SL.close("t3_I2_model_De_vs_raw_floor_rtol1e-12", Dm, D_ai(max(c, CMIN), T), rtol=1e-12, atol=0.0),
        SL.cond("t3_I3_raw_finite_positive", math.isfinite(D) and D > 0, observed=D, expected=">0"),
        SL.cond("t3_I3_model_finite_positive", math.isfinite(Dm) and Dm > 0, observed=Dm, expected=">0"),
        SL.close("t3_I5_kappa_vs_poly_rtol1e-12", k, t3_oracle_kappa(c, T), rtol=1e-12, atol=0.0),
        SL.cond("t3_I5_kappa_finite_nonneg", math.isfinite(k) and k >= 0, observed=k, expected=">=0"),
        SL.close("t3_I5_kappa_model_vs_poly_floor_rtol1e-12", km, t3_oracle_kappa(max(c, CMIN), T),
                 rtol=1e-12, atol=0.0),
    ]
    if 501 <= c <= 1999 and 273.16 <= T <= 333.14:
        hc, hT = 0.1, 0.001
        dlT = (math.log(D_ai(c, T + hT)) - math.log(D_ai(c, T - hT))) / (2 * hT)
        dlc = (math.log(D_ai(c + hc, T)) - math.log(D_ai(c - hc, T))) / (2 * hc)
        aT, ac = t3_dlnD(c, T)
        ch.append(SL.close("t3_I4_dlnD_dT_vs_analytic_rtol1e-6", dlT, aT, rtol=1e-6, atol=0.0))
        ch.append(SL.close("t3_I4_dlnD_dc_vs_analytic_rtol1e-6", dlc, ac, rtol=1e-6, atol=0.0))
    return ch


def t3_sample(rng):
    return dict(c_e=float(rng.uniform(500, 2000)), T=float(rng.uniform(273.15, 333.15)))


def t3_grid_monotone_eval(inp):
    vals = np.array([[D_ai(c, T) for T in T3_T] for c in T3_C])
    ch = []
    for i, c in enumerate(T3_C):
        ch.append(SL.cond(f"t3_I4_row_strict_increase_in_T_c={c}", bool(np.all(np.diff(vals[i]) > 0)),
                          observed=vals[i].tolist()))
    for j, T in enumerate(T3_T):
        ch.append(SL.cond(f"t3_I4_col_strict_decrease_in_c_T={T}", bool(np.all(np.diff(vals[:, j]) < 0)),
                          observed=vals[:, j].tolist()))
    return ch


def t3_nominal_eval(inp):
    return [
        SL.close("t3_I1_nominal_D(1000,298.15)_rtol1e-12", D_ai(1000.0, 298.15), 3.222722528605257e-10,
                 rtol=1e-12, atol=0.0),
        SL.close("t3_I5_nominal_kappa(1000,298.15)_rel1e-3", kappa_ai(1000.0, 298.15), 1.1943, rtol=1e-3, atol=0.0),
    ]


def t3_boundary_eval(inp):
    c, T = inp["c_e"], inp["T"]
    D = D_ai(c, T)
    Dm = D_model(c, T)
    k = kappa_ai(c, T)
    km = kappa_model(c, T)
    ch = [
        SL.cond("t3_I2_setting_D_e__c_e==10", pybamm.settings.tolerances["D_e__c_e"] == 10,
                observed=pybamm.settings.tolerances["D_e__c_e"], expected=10),
        SL.cond("t3_I5_setting_kappa_e__c_e==10", pybamm.settings.tolerances["kappa_e__c_e"] == 10,
                observed=pybamm.settings.tolerances["kappa_e__c_e"], expected=10),
        SL.close("t3_I2_model_De_vs_raw_floor_rtol1e-12", Dm, D_ai(max(c, CMIN), T), rtol=1e-12, atol=0.0),
        SL.close("t3_I1_api_vs_oracle_rtol1e-12", D, t3_oracle_D(c, T), rtol=1e-12, atol=0.0),
        SL.cond("t3_I3_raw_finite_positive", math.isfinite(D) and D > 0, observed=D, expected=">0"),
        SL.cond("t3_I3_model_finite_positive", math.isfinite(Dm) and Dm > 0, observed=Dm, expected=">0"),
        SL.close("t3_I5_kappa_vs_poly_rtol1e-12", k, t3_oracle_kappa(c, T), rtol=1e-12, atol=0.0),
        SL.cond("t3_I5_kappa_finite_nonneg", math.isfinite(k) and k >= 0, observed=k, expected=">=0"),
        SL.close("t3_I5_kappa_model_vs_poly_floor_rtol1e-12", km, t3_oracle_kappa(max(c, CMIN), T),
                 rtol=1e-12, atol=0.0),
    ]
    if c == 0:
        ch.append(SL.cond("t3_I5_raw_kappa(0,T)==0_exact", k == 0.0, observed=k, expected=0.0))
    if c < CMIN:
        ch.append(SL.close("t3_I4_wrapper_constant_below_floor_rtol1e-12", Dm, D_model(CMIN, T), rtol=1e-12, atol=0.0))
    if c > CMIN:
        prev = max(x for x in T3_BC if x < c)
        Dprev = D_model(prev, T)
        ch.append(SL.cond("t3_I4_wrapper_strict_decrease_above_floor", Dm < Dprev, observed=[Dprev, Dm],
                          expected=f"D_e({c}) < D_e({prev})"))
    return ch


def t3_vector_eval(inp):
    """On the grid, two combinations (paired Vector, scalar×Vector), compared element-wise with the scalar call (rtol 1e-12)."""
    cc, TT = np.meshgrid(np.array(T3_C, float), np.array(T3_T, float), indexing="ij")
    cf, Tf = cc.ravel(), TT.ravel()
    ch = []
    for key, fn in ((DK, D_ai), (KK, kappa_ai)):
        scal = np.array([fn(c, T) for c, T in zip(cf, Tf)])
        vec = np.asarray(AI.evaluate(AI[key](pybamm.Vector(cf), pybamm.Vector(Tf))), float).reshape(-1)
        ch.append(SL.close(f"t3_vector_paired_{key[:20]}", vec, scal, rtol=1e-12, atol=0.0))
        for T in T3_T:
            v = np.asarray(AI.evaluate(AI[key](pybamm.Vector(np.array(T3_C, float)), pybamm.Scalar(T))), float).reshape(-1)
            s = np.array([fn(c, T) for c in T3_C])
            ch.append(SL.close(f"t3_vector_c_scalar_T={T}_{key[:20]}", v, s, rtol=1e-12, atol=0.0))
        for c in T3_C:
            v = np.asarray(AI.evaluate(AI[key](pybamm.Scalar(c), pybamm.Vector(np.array(T3_T, float)))), float).reshape(-1)
            s = np.array([fn(c, T) for T in T3_T])
            ch.append(SL.close(f"t3_vector_scalar_c={c}_T_{key[:20]}", v, s, rtol=1e-12, atol=0.0))
    return ch


# ---------------------------------------------------------------- Assembly
def main():
    a = SL.parse_args()
    run = SL.Runner(PR, PR_IDX, a.track, a.n, a.out, versions())
    if a.track == 1:
        branches = [("Ai2020_constitutive", t1_sample, t1_eval)]
        mand = [(f"grid_c={c}_T={T}", dict(c_e=float(c), T=float(T)), t1_eval)
                for c in range(100, 3001, 100) for T in range(253, 334, 5)]
        mand.append(("optional_dfn_1C", dict(case="DFN 1C Ai2020 3600 s"), t1_dfn_eval))
    elif a.track == 2:
        branches = [("Ai2020_constitutive+DFN", t2_sample, t2_eval, t2_replay)]
        mand = [("nominal_1000_298.15", dict(c_e=1000.0, T=298.15), t2_nominal_eval),
                ("C_rate=1_c=1000_T=298.15", dict(c_e=1000.0, T=298.15, C_rate=1.0), t2_eval)]
    else:
        branches = [("Ai2020_raw_and_wrapper", t3_sample, t3_eval)]
        mand = [(f"grid_c={c}_T={T}", dict(c_e=float(c), T=float(T)), t3_eval) for c in T3_C for T in T3_T]
        mand.append(("grid_monotone", dict(grid="7x5"), t3_grid_monotone_eval))
        mand.append(("nominal", dict(c_e=1000.0, T=298.15), t3_nominal_eval))
        mand += [(f"boundary_c={c}_T={T}", dict(c_e=float(c), T=float(T)), t3_boundary_eval)
                 for c in T3_BC for T in T3_T]
        mand.append(("vector_grid", dict(grid="7x5"), t3_vector_eval))
    run.run(branches, mand, progress_every=100 if a.track == 2 else 1000)


if __name__ == "__main__":
    main()
