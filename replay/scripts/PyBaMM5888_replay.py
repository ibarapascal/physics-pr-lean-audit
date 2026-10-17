#!/usr/bin/env python3
"""
Replay / differential check for lean/PhysAudit/PyBaMM5888.lean
(pybamm-team/PyBaMM#5888, merge commit 343ece2062dbb39789dbfa977e0f2daf96e9e370).

Compares the native PyBaMM electrolyte diffusivity of the Ai2020 parameter set with the
formulas transcribed in the Lean module, and checks each Lean theorem numerically:
  (a) Ai2020 == Xu2019 (Valoen2005) on random (c_e, T)
  (b) Ai2020 == (1 cm)^2 * Valoen-Reimers eq. 14 (cgs, mol/L); parent version == 1e4 * merged
  (c) dD/dT > 0 and dD/dc_e < 0 on c_e in [0, 4000] mol/m3, T in [253, 333] K
      (analytic derivative from the Lean theorems vs native central differences),
      strict monotonicity on random pairs, and the two counterexamples across the pole
      T = 229 + 5e-3 * c_e.

Dependencies: pybamm==26.10.1.0 (the tag pybamm-v26.10.1.0 contains the merge commit),
numpy; mpmath is optional (used for a 50-digit evaluation of the Lean formula, float64 otherwise).
Run:
    <py2-venv>/bin/python PyBaMM5888_replay.py > PyBaMM5888_replay.json
Output: one JSON object on stdout.
"""
import inspect
import json
import math
import sys
import warnings

import numpy as np

try:
    import mpmath as mp

    mp.mp.dps = 50
except ImportError:  # pragma: no cover
    mp = None

try:
    import pybamm
except ImportError:
    print(json.dumps({"status": "SKIP", "reason": "pybamm not importable"}))
    sys.exit(0)

NAME = "Electrolyte diffusivity [m2.s-1]"
N_POINTS = 2000
SEED = 5888
C_RANGE = (0.0, 4000.0)
T_RANGE = (253.0, 333.0)


def ev(pv, sym):
    v = pv.evaluate(sym)
    return float(np.asarray(v, dtype=float).reshape(-1)[0])


def native(pv, c, T):
    return ev(pv, pv[NAME](pybamm.Scalar(c), pybamm.Scalar(T)))


# ---- Lean definitions (PyBaMM5888.aiExponent / diffAi2020 / diffAi2020Old / vrLog10D) ----
if mp is not None:
    M = mp.mpf

    def lean_exponent(c, T):
        c, T = M(c), M(T)
        return M("-4.43") - M(54) / (T - 229 - M("5e-3") * c) - M("0.22e-3") * c

    def lean_D(c, T):
        return M("1e-4") * mp.power(10, lean_exponent(c, T))

    def lean_D_old(c, T):
        return mp.power(10, lean_exponent(c, T))

    def lean_vr_SI(c, T):
        cl = M(c) / 1000
        e = M("-4.43") - M(54) / (M(T) - 229 - M("5.0") * cl) - M("0.22") * cl
        return (M(1) / 100) ** 2 * mp.power(10, e)

    def lean_xu(c, T):
        cm = M(c) / 1000
        Tg = 229 + 5 * cm
        return mp.power(10, (M("-4.43") - M(54) / (M(T) - Tg)) + M("-0.22") * cm) * M("1e-4")

    def lean_dDdT(c, T):
        u = M(T) - 229 - M("5e-3") * M(c)
        return lean_D(c, T) * mp.log(10) * (M(54) / u**2)

    def lean_dDdc(c, T):
        u = M(T) - 229 - M("5e-3") * M(c)
        return lean_D(c, T) * mp.log(10) * (-(M(54) * M("5e-3")) / u**2 - M("0.22e-3"))

    to_f = float
else:

    def lean_exponent(c, T):
        return -4.43 - 54.0 / (T - 229 - 5e-3 * c) - 0.22e-3 * c

    def lean_D(c, T):
        return 1e-4 * 10.0 ** lean_exponent(c, T)

    def lean_D_old(c, T):
        return 10.0 ** lean_exponent(c, T)

    def lean_vr_SI(c, T):
        cl = c / 1000
        return (1 / 100) ** 2 * 10.0 ** (-4.43 - 54.0 / (T - 229 - 5.0 * cl) - 0.22 * cl)

    def lean_xu(c, T):
        cm = c / 1000
        return 10.0 ** ((-4.43 - 54.0 / (T - (229 + 5 * cm))) - 0.22 * cm) * 1e-4

    def lean_dDdT(c, T):
        u = T - 229 - 5e-3 * c
        return lean_D(c, T) * math.log(10) * 54.0 / u**2

    def lean_dDdc(c, T):
        u = T - 229 - 5e-3 * c
        return lean_D(c, T) * math.log(10) * (-(54 * 5e-3) / u**2 - 0.22e-3)

    to_f = float


def rel(a, b):
    a, b = float(a), float(b)
    if a == b:
        return 0.0
    return abs(a - b) / max(abs(a), abs(b))


def main():
    out = {"module": "PyBaMM5888", "pybamm_version": pybamm.__version__,
           "lean_formula_precision": "mpmath-50-digits" if mp is not None else "float64"}

    pv_ai = pybamm.ParameterValues("Ai2020")
    pv_xu = pybamm.ParameterValues("Xu2019")
    f_ai = pv_ai[NAME]
    try:
        src = inspect.getsource(f_ai)
    except (OSError, TypeError):
        src = ""
    out["native_function"] = getattr(f_ai, "__qualname__", str(f_ai))
    out["source_contains_fix"] = "1e-4 * 10 ** (-4.43 - 54 / (T - 229 - 5e-3 * c_e) - 0.22e-3 * c_e)" in src

    rng = np.random.default_rng(SEED)
    cs = rng.uniform(*C_RANGE, N_POINTS)
    Ts = rng.uniform(*T_RANGE, N_POINTS)

    max_rel_native_vs_lean = 0.0
    max_rel_ai_vs_xu_native = 0.0
    max_rel_xu_native_vs_lean_xu = 0.0
    max_rel_lean_vs_literature = 0.0
    max_rel_old_ratio = 0.0
    max_rel_dT = 0.0
    max_rel_dc = 0.0
    sign_dT_ok = 0
    sign_dc_ok = 0
    min_u = math.inf
    D_min, D_max = math.inf, -math.inf
    worst = None

    dT, dc = 0.01, 0.1
    for c, T in zip(cs, Ts):
        c, T = float(c), float(T)
        d_nat = native(pv_ai, c, T)
        d_lean = lean_D(c, T)
        r = rel(d_nat, d_lean)
        if r > max_rel_native_vs_lean:
            max_rel_native_vs_lean = r
            worst = {"c_e": c, "T": T, "native": d_nat, "lean": to_f(d_lean)}
        d_xu = native(pv_xu, c, T)
        max_rel_ai_vs_xu_native = max(max_rel_ai_vs_xu_native, rel(d_nat, d_xu))
        max_rel_xu_native_vs_lean_xu = max(max_rel_xu_native_vs_lean_xu, rel(d_xu, lean_xu(c, T)))
        max_rel_lean_vs_literature = max(max_rel_lean_vs_literature, rel(d_lean, lean_vr_SI(c, T)))
        max_rel_old_ratio = max(max_rel_old_ratio, rel(lean_D_old(c, T) / d_nat, 1e4))
        D_min, D_max = min(D_min, d_nat), max(D_max, d_nat)
        min_u = min(min_u, T - 229 - 5e-3 * c)

        # native central differences (native function called directly on floats)
        fd_T = (f_ai(c, T + dT) - f_ai(c, T - dT)) / (2 * dT)
        fd_c = (f_ai(c + dc, T) - f_ai(c - dc, T)) / (2 * dc)
        an_T = lean_dDdT(c, T)
        an_c = lean_dDdc(c, T)
        sign_dT_ok += int(fd_T > 0 and an_T > 0)
        sign_dc_ok += int(fd_c < 0 and an_c < 0)
        max_rel_dT = max(max_rel_dT, rel(fd_T, an_T))
        max_rel_dc = max(max_rel_dc, rel(fd_c, an_c))

    out["n_points"] = N_POINTS
    out["domain"] = {"c_e_mol_m3": list(C_RANGE), "T_K": list(T_RANGE)}
    out["theorem_checks"] = {
        "native_vs_lean_diffAi2020_max_rel": max_rel_native_vs_lean,
        "native_vs_lean_worst_point": worst,
        "diffAi2020_eq_diffXu2019_native_max_rel": max_rel_ai_vs_xu_native,
        "xu2019_native_vs_lean_max_rel": max_rel_xu_native_vs_lean_xu,
        "diffAi2020_eq_SI_literature_max_rel": max_rel_lean_vs_literature,
        "diffAi2020Old_eq_1e4_mul_max_rel": max_rel_old_ratio,
        "box_subset_domain_min_denominator_K": min_u,
        "box_subset_domain_holds": min_u > 0,
        "deriv_T_pos_count": sign_dT_ok,
        "deriv_c_neg_count": sign_dc_ok,
        "hasDerivAt_T_formula_vs_native_fd_max_rel": max_rel_dT,
        "hasDerivAt_c_formula_vs_native_fd_max_rel": max_rel_dc,
        "native_D_range_m2_s": [D_min, D_max],
    }

    # strict monotonicity on random pairs (native ParameterValues path)
    mono_T_ok = 0
    anti_c_ok = 0
    n_pairs = 500
    for _ in range(n_pairs):
        c = float(rng.uniform(*C_RANGE))
        T1, T2 = sorted(rng.uniform(*T_RANGE, 2))
        mono_T_ok += int(native(pv_ai, c, float(T1)) < native(pv_ai, c, float(T2)))
        T = float(rng.uniform(*T_RANGE))
        c1, c2 = sorted(rng.uniform(*C_RANGE, 2))
        anti_c_ok += int(native(pv_ai, float(c2), T) < native(pv_ai, float(c1), T))
    out["theorem_checks"]["strictMonoOn_T_box_pairs_ok"] = f"{mono_T_ok}/{n_pairs}"
    out["theorem_checks"]["strictAntiOn_c_box_pairs_ok"] = f"{anti_c_ok}/{n_pairs}"

    # counterexamples across the pole
    ce = {}
    for (c, T) in [(0.0, 228.0), (0.0, 230.0), (3000.0, 250.0), (5000.0, 250.0)]:
        ce[f"c={c},T={T}"] = {"native": native(pv_ai, c, T), "lean": to_f(lean_D(c, T)),
                              "denominator": T - 229 - 5e-3 * c}
    out["counterexamples"] = {
        "not_monotone_T_across_pole": {
            "points": {k: ce[k] for k in ["c=0.0,T=228.0", "c=0.0,T=230.0"]},
            "native_D(0,230) < native_D(0,228)": ce["c=0.0,T=230.0"]["native"] < ce["c=0.0,T=228.0"]["native"],
        },
        "not_antitone_c_across_pole": {
            "points": {k: ce[k] for k in ["c=3000.0,T=250.0", "c=5000.0,T=250.0"]},
            "native_D(3000,250) < native_D(5000,250)": ce["c=3000.0,T=250.0"]["native"] < ce["c=5000.0,T=250.0"]["native"],
        },
    }

    # behaviour exactly at the pole (Lean uses x/0 = 0, so values differ there by construction)
    pole = {}
    try:
        pole["direct_float_call"] = float(f_ai(0.0, 229.0))
    except ZeroDivisionError as e:
        pole["direct_float_call"] = f"ZeroDivisionError: {e}"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            pole["parameter_values_evaluate"] = native(pv_ai, 0.0, 229.0)
        except Exception as e:  # noqa: BLE001
            pole["parameter_values_evaluate"] = f"{type(e).__name__}: {e}"
    pole["lean_value_with_x_div_0_eq_0"] = to_f(M("1e-4") * mp.power(10, M("-4.43"))) if mp is not None else 1e-4 * 10 ** -4.43
    out["pole_c0_T229"] = pole

    # regression values from the PR test
    out["pr_test_value_1000_298.15"] = {"native": native(pv_ai, 1000.0, 298.15), "expected": 3.2227e-10,
                                        "rel": rel(native(pv_ai, 1000.0, 298.15), 3.2227e-10)}

    tc = out["theorem_checks"]
    out["all_certificate_checks_pass"] = bool(
        out["source_contains_fix"]
        and tc["native_vs_lean_diffAi2020_max_rel"] < 1e-12
        and tc["diffAi2020_eq_diffXu2019_native_max_rel"] < 1e-12
        and tc["diffAi2020_eq_SI_literature_max_rel"] < 1e-12
        and tc["diffAi2020Old_eq_1e4_mul_max_rel"] < 1e-12
        and tc["box_subset_domain_holds"]
        and tc["deriv_T_pos_count"] == N_POINTS
        and tc["deriv_c_neg_count"] == N_POINTS
        and mono_T_ok == n_pairs
        and anti_c_ok == n_pairs
    )
    out["counterexamples_reproduced"] = bool(
        out["counterexamples"]["not_monotone_T_across_pole"]["native_D(0,230) < native_D(0,228)"]
        and out["counterexamples"]["not_antitone_c_across_pole"]["native_D(3000,250) < native_D(5000,250)"]
    )
    print(json.dumps(out, indent=2, default=str))


if __name__ == "__main__":
    main()
