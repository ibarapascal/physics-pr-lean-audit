#!/usr/bin/env python3
"""
Replay / differential check for lean/PhysAudit/PyBaMM5889.lean
(pybamm-team/PyBaMM#5889, merge commit eda03ed6d1003b186cc2d70ad8b3240a7b6350f1).

Compares the native PyBaMM positive particle diffusivity of the Chayambuka2022 parameter set
(NVPF_diffusivity_Chayambuka2022, pybamm.Interpolant -> scipy interp1d linear with extrapolation)
with the Lean transcription (searchsorted-left count, clip to [1, N-1], segment line), evaluated
exactly with fractions.Fraction on the float values of D_p.csv, and checks each Lean theorem:
  interp_at_node / nvpf_at_node, interp_pos_on_axis / interp_bounded_on_axis / nvpf_on_axis,
  lookup_indep_cmax, old_uses_first_three_nodes / old_ignores_upper_nodes,
  nvpf_pos_of_sto_le_one, and the Arithmetic constants (axis endpoints, extrapolation margins,
  old argument below the axis).
The parent-commit lookup (sto * initial electrolyte concentration) is rebuilt in this script
with the same pybamm.Interpolant call and installed through ParameterValues.update.

Dependencies: pybamm==26.10.1.0 (the tag pybamm-v26.10.1.0 contains the merge commit), numpy.
Run:
    <py2-venv>/bin/python PyBaMM5889_replay.py > PyBaMM5889_replay.json
Output: one JSON object on stdout.
"""
import inspect
import json
import sys
from fractions import Fraction as Fr

import numpy as np

try:
    import pybamm
except ImportError:
    print(json.dumps({"status": "SKIP", "reason": "pybamm not importable"}))
    sys.exit(0)

NAME = "Positive particle diffusivity [m2.s-1]"
CMAX_NAME = "Maximum concentration in positive electrode [mol.m-3]"
CE0_NAME = "Initial concentration in electrolyte [mol.m-3]"
N_POINTS = 2000
SEED = 5889

# Lean constants (PyBaMM5889 §2, §5)
C_POS_MAX = Fr(15320)
C_ELEC0 = Fr(1000)
DP_X0, DP_Y0 = Fr(131578947, 10**6), Fr(251, 10**17)
DP_X1, DP_Y1 = Fr(592105263, 10**6), Fr(171, 10**17)
DP_X2 = Fr(1052631579, 10**6)
DP_XM, DP_YM = Fr(15000), Fr(229, 10**19)
DP_XN, DP_YN = Fr(1519736842, 10**5), Fr(195, 10**19)


def ev(pv, sym):
    v = pv.evaluate(sym)
    return float(np.asarray(v, dtype=float).reshape(-1)[0])


def native(pv, sto, T=298.15):
    return ev(pv, pv[NAME](pybamm.Scalar(sto), pybamm.Scalar(T)))


# ---- Lean transcription (PyBaMM5889 §1) ----
def searchsorted_left(xs, x):
    return sum(1 for a in xs if a < x)


def interp_index(xs, x):
    return min(max(searchsorted_left(xs, x), 1), len(xs) - 1)


def segment(xs, ys, k, x):
    return (ys[k] - ys[k - 1]) / (xs[k] - xs[k - 1]) * (x - xs[k - 1]) + ys[k - 1]


def interp_lean(xs, ys, x):
    return segment(xs, ys, interp_index(xs, x), x)


def rel(a, b):
    a, b = float(a), float(b)
    if a == b:
        return 0.0
    return abs(a - b) / max(abs(a), abs(b))


def main():
    out = {"module": "PyBaMM5889", "pybamm_version": pybamm.__version__}
    pv = pybamm.ParameterValues("Chayambuka2022")
    f_new = pv[NAME]
    mod = sys.modules[f_new.__module__]
    try:
        src = inspect.getsource(f_new)
        src_hc = inspect.getsource(pv["Negative particle diffusivity [m2.s-1]"])
    except (OSError, TypeError):
        src, src_hc = "", ""
    out["source_contains_fix"] = 'pybamm.Parameter("Maximum concentration in positive electrode [mol.m-3]")' in src
    out["hc_uses_negative_cmax"] = 'pybamm.Parameter("Maximum concentration in negative electrode [mol.m-3]")' in src_hc
    out["c_max_positive"] = float(pv[CMAX_NAME])
    out["c_e0"] = float(pv[CE0_NAME])

    name, (x_list, y) = mod.D_p_data
    x = np.asarray(x_list[0], dtype=float)
    y = np.asarray(y, dtype=float)
    N = len(x)
    XF = [Fr(float(v)) for v in x]
    YF = [Fr(float(v)) for v in y]
    out["table"] = {
        "name": name, "N": N,
        "strictly_increasing": bool(np.all(np.diff(x) > 0)),
        "all_values_positive": bool(np.all(y > 0)),
        "x_first": float(x[0]), "x_last": float(x[-1]),
        "y_min": float(y.min()), "y_max": float(y.max()),
        "endpoint_constants_match_float": {
            "dpX0": float(DP_X0) == x[0], "dpY0": float(DP_Y0) == y[0],
            "dpX1": float(DP_X1) == x[1], "dpY1": float(DP_Y1) == y[1],
            "dpX2": float(DP_X2) == x[2],
            "dpXm": float(DP_XM) == x[-2], "dpYm": float(DP_YM) == y[-2],
            "dpXN": float(DP_XN) == x[-1], "dpYN": float(DP_YN) == y[-1],
        },
    }

    rng = np.random.default_rng(SEED)

    # 1. native vs Lean on random sto in [0, 1] (T random, unused by the function)
    stos = rng.uniform(0.0, 1.0, N_POINTS)
    Ts = rng.uniform(273.0, 323.0, N_POINTS)
    max_rel_exact = 0.0
    max_rel_float_arg = 0.0
    worst = None
    min_native = np.inf
    for s, T in zip(stos, Ts):
        s, T = float(s), float(T)
        d_nat = native(pv, s, T)
        d_exact = interp_lean(XF, YF, Fr(s) * C_POS_MAX)          # real-number argument
        d_farg = interp_lean(XF, YF, Fr(s * float(C_POS_MAX)))   # argument rounded as in float64
        r = rel(d_nat, d_exact)
        if r > max_rel_exact:
            max_rel_exact, worst = r, {"sto": s, "native": d_nat, "lean": float(d_exact)}
        max_rel_float_arg = max(max_rel_float_arg, rel(d_nat, d_farg))
        min_native = min(min_native, d_nat)
    out["native_vs_lean"] = {
        "n_points": N_POINTS, "sto_range": [0.0, 1.0],
        "max_rel_lean_exact_argument": max_rel_exact,
        "max_rel_lean_float64_argument": max_rel_float_arg,
        "worst_point": worst,
    }

    # 2. node reproduction (interp_at_node exactly in Fraction; nvpf_at_node natively)
    exact_nodes = all(interp_lean(XF, YF, XF[m]) == YF[m] for m in range(N))
    exact_nvpf_nodes = all(interp_lean(XF, YF, (XF[m] / C_POS_MAX) * C_POS_MAX) == YF[m] for m in range(N))
    node_rel = [rel(native(pv, float(x[m]) / 15320.0), y[m]) for m in range(N)]
    out["interp_at_node_exact_fraction"] = exact_nodes
    out["nvpf_at_node_exact_fraction"] = exact_nvpf_nodes
    out["nvpf_at_node_native_max_rel"] = max(node_rel)

    # 3. on-axis bounds and positivity (interp_bounded_on_axis / nvpf_on_axis)
    s_lo, s_hi = float(x[0]) / 15320.0, float(x[-1]) / 15320.0
    s_axis = np.linspace(s_lo, s_hi, 2001)
    vals = np.array([native(pv, float(s)) for s in s_axis])
    out["nvpf_on_axis"] = {
        "sto_interval": [s_lo, s_hi],
        "native_min": float(vals.min()), "native_max": float(vals.max()),
        "within_table_range": bool(vals.min() >= y.min() * (1 - 1e-12) and vals.max() <= y.max() * (1 + 1e-12)),
        "positive": bool(vals.min() > 0),
    }

    # 4. nvpf_pos_of_sto_le_one: dense sto grid in [-0.5, 1] including both extrapolation zones
    s_full = np.linspace(-0.5, 1.0, 3001)
    v_full = np.array([native(pv, float(s)) for s in s_full])
    v_unit = np.array([native(pv, float(s)) for s in np.linspace(0.0, 1.0, 2001)])
    out["nvpf_pos_of_sto_le_one"] = {
        "native_min_on_[-0.5,1]": float(v_full.min()),
        "native_min_on_[0,1]": float(v_unit.min()),
        "native_at_sto_0": native(pv, 0.0), "native_at_sto_1": native(pv, 1.0),
        "lean_at_sto_1_exact": float(interp_lean(XF, YF, C_POS_MAX)),
        "lean_right_segment_value_at_15320_from_constants": float(
            (DP_YN - DP_YM) / (DP_XN - DP_XM) * (Fr(15320) - DP_XM) + DP_YM),
        "positive": bool(v_full.min() > 0),
    }

    # 5. lookup_indep_cmax: fixed solid concentration, different c_max
    cmax_res = {}
    c_s = rng.uniform(float(x[0]), float(x[-1]), 200)
    base = np.array([native(pv, float(c) / 15320.0) for c in c_s])
    for cm in [10000.0, 20000.0, 22980.0]:
        pv2 = pv.copy()
        pv2.update({CMAX_NAME: cm})
        alt = np.array([native(pv2, float(c) / cm) for c in c_s])
        cmax_res[str(cm)] = float(np.max(np.abs(alt - base) / np.abs(base)))
    out["lookup_indep_cmax_max_rel"] = cmax_res

    # 6. parent-commit lookup rebuilt natively
    def make_old(xx, yy):
        def nvpf_old(sto, T):
            c_max = pybamm.Parameter(CE0_NAME)
            return pybamm.Interpolant(xx, yy, sto * c_max, name)
        return nvpf_old

    pv_old = pv.copy()
    pv_old.update({NAME: make_old([x], y)})
    y_mod = y.copy()
    y_mod[3:] *= 7.0
    pv_old_mod = pv.copy()
    pv_old_mod.update({NAME: make_old([x], y_mod)})

    def lean_old_formula(s):
        a = Fr(s) * C_ELEC0
        return segment(XF, YF, 1, a) if a <= XF[1] else segment(XF, YF, 2, a)

    s_old = rng.uniform(0.0, 1.0, N_POINTS)
    max_rel_old_vs_lean = 0.0
    max_rel_old_vs_three_node_formula = 0.0
    max_abs_old_vs_modified_table = 0.0
    below_axis = 0
    for s in s_old:
        s = float(s)
        d_old = native(pv_old, s)
        max_rel_old_vs_lean = max(max_rel_old_vs_lean, rel(d_old, interp_lean(XF, YF, Fr(s) * C_ELEC0)))
        max_rel_old_vs_three_node_formula = max(max_rel_old_vs_three_node_formula, rel(d_old, lean_old_formula(s)))
        max_abs_old_vs_modified_table = max(max_abs_old_vs_modified_table, abs(d_old - native(pv_old_mod, s)))
        below_axis += int(s * 1000.0 < x[0])
    grid = np.linspace(0.01, 0.99, 99)
    ratio = np.array([native(pv_old, float(s)) / native(pv, float(s)) for s in grid])
    out["parent_commit"] = {
        "native_old_vs_lean_max_rel": max_rel_old_vs_lean,
        "old_uses_first_three_nodes_max_rel": max_rel_old_vs_three_node_formula,
        "old_ignores_upper_nodes_max_abs_diff_after_scaling_nodes_3_onward_by_7": max_abs_old_vs_modified_table,
        "old_argument_range": [0.0, 1000.0],
        "old_argument_max_below_x2": bool(1000.0 < x[2]),
        "fraction_of_random_sto_with_old_argument_below_axis": below_axis / N_POINTS,
        "sto_threshold_old_below_axis": float(x[0]) / 1000.0,
        "old_over_new_ratio_on_sto_0.01_to_0.99": [float(ratio.min()), float(ratio.max())],
        "native_old_at_sto_0.5": native(pv_old, 0.5),
    }

    # 7. Arithmetic section constants
    out["arithmetic"] = {
        "extrapolation_margins_sto": [float(DP_X0 / C_POS_MAX), float(DP_XN / C_POS_MAX)],
        "new_arg_inside_axis_for_sto_0.01_0.99": bool(0.01 * 15320 >= float(DP_X0) and 0.99 * 15320 <= float(DP_XN)),
        "old_vs_new_example_sto_0.1": {"old_arg": 100.0, "new_arg": 1532.0, "x0": float(DP_X0)},
        "cElec0_lt_dpX2": bool(C_ELEC0 < DP_X2),
        "old_reach_fraction": float((C_ELEC0 - DP_X0) / (DP_XN - DP_X0)),
    }

    # 8. PR regression values
    out["pr_test_values"] = {
        str(s): {"native": native(pv, s), "expected": e, "rel": rel(native(pv, s), e)}
        for s, e in [(0.5, 8.3317e-17), (0.3, 5.6673e-17), (0.9, 3.5627e-17)]
    }

    t = out["table"]
    out["all_certificate_checks_pass"] = bool(
        out["source_contains_fix"]
        and t["strictly_increasing"] and t["all_values_positive"]
        and all(t["endpoint_constants_match_float"].values())
        and max_rel_float_arg < 1e-12
        and exact_nodes and exact_nvpf_nodes
        and out["nvpf_on_axis"]["within_table_range"] and out["nvpf_on_axis"]["positive"]
        and out["nvpf_pos_of_sto_le_one"]["positive"]
        and max(cmax_res.values()) < 1e-12
        and max_rel_old_vs_three_node_formula < 1e-12
        and max_abs_old_vs_modified_table == 0.0
    )
    print(json.dumps(out, indent=2, default=str))


if __name__ == "__main__":
    main()
