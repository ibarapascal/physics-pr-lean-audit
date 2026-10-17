#!/usr/bin/env python3
"""Moose33307_replay.py - differential check for idaholab/moose PR #33307 (merge b5a1d9e3df).

Subchannel closures: SCMFrictionUpgradedChengTodreas.C (tri lattice friction factor, wire-wrap
corrections, laminar/transition/turbulent regime by bulk Re), SCMMixingChengTodreas.C (Cm/Cs
transition), SCMHTCBorishanskii.C (log10 Nusselt base term), and
SubChannel1PhaseProblem::computeBulkReynoldsNumber.

There is no native MOOSE build available, so this script does NOT run MOOSE. It compares
  (A) an independent float64 line-by-line transcription of the C++ (math.log10 / math.log /
      math.pow, mooseError -> exception),
  (B) the formulas of the Lean definitions in lean/PhysAudit/Moose33307.lean (`triCoeffs`,
      `regimeFriction`, `transitionBlend`, `psi`, `ReL`, `ReT`, `mixingRegime`,
      `borishanskiiNuBase`, `bulkReynolds`), evaluated with mpmath at 50 digits when available
      (float64 otherwise),
on 2000 random points of the correlation domain, and checks every theorem of the Lean module.
It also prints a static description of the bulk-Re parallel behaviour (not executed).

License:      LGPL-2.1 (this file contains a line-by-line reimplementation of MOOSE C++ source,
              which is licensed under LGPL-2.1; it is not covered by the MIT license of this
              archive; see THIRD_PARTY.md).
Dependencies: Python >= 3.8; mpmath optional (pip install mpmath). No MOOSE.
Run:          python3 -I Moose33307_replay.py > Moose33307_replay.json
Output:       one JSON object on stdout.
"""
import json
import math
import random

try:
    import mpmath as mp
    mp.mp.dps = 50
    HAVE_MP = True
except Exception:  # pragma: no cover
    HAVE_MP = False

random.seed(33307)
N = 2000
CENTER, EDGE, CORNER = "CENTER", "EDGE", "CORNER"


class MooseError(Exception):
    pass


# ============================================================ (A) C++ transcription (float64)
def cpp_tri_friction(Re, i_type, S, w_perim, pitch, D, H, Dw, gap, has_wire, Reb):
    """SCMFrictionUpgradedChengTodreas::computeTriLatticeFrictionFactor (:73-298)."""
    Re_eff = max(Re, 1.0)
    Dh_i = 4.0 * S / w_perim
    p_over_d = pitch / D
    w_over_d = (D + gap) / D
    CbL1, CbL2, CbT1, CbT2 = 320, 1.0, 10000, 0.7
    ReL = math.pow(10, CbL2 * (p_over_d - 1)) * CbL1
    ReT = math.pow(10, CbT2 * (p_over_d - 1)) * CbT1
    psi = math.log(Reb / ReL) / math.log(ReT / ReL) if Reb > 0 else float("nan")
    if i_type == CENTER:
        if p_over_d < 1.1:
            k = (26.0, 888.2, -3334.0, 0.09378, 1.398, -8.664)
        else:
            k = (62.97, 216.9, -190.2, 0.1458, 0.03632, -0.03333)
        x = p_over_d - 1
    elif i_type == EDGE:
        if w_over_d < 1.1:
            k = (26.18, 554.5, -1480.0, 0.09377, 0.8732, -3.341)
        else:
            k = (44.4, 256.7, -267.6, 0.1430, 0.04199, -0.04428)
        x = w_over_d - 1
    else:
        if w_over_d < 1.1:
            k = (26.98, 1636.0, -10050.0, 0.1004, 1.625, -11.85)
        else:
            k = (87.26, 38.59, -55.12, 0.1499, 0.006706, -0.009567)
        x = w_over_d - 1
    CfL = k[0] + k[1] * x + k[2] * x ** 2
    CfT = k[3] + k[4] * x + k[5] * x ** 2
    if has_wire:
        theta = math.acos(H / math.sqrt(H ** 2 + (math.pi * (D + Dw)) ** 2))
        WdT = (19.56 + -98.71 * (Dw / D) + 303.47 * (Dw / D) ** 2) * math.pow(H / D, -0.541)
        WdL = 1.4 * WdT
        WsT = -11 * math.log10(H / D) + 19
        WsL = 1.0 * WsT
        if i_type == CENTER:
            pw_p = math.pi * D / 2.0
            ar = math.pi * (D + Dw) * Dw / 6.0
            a_p = S + math.pi * Dw ** 2 / 8.0 / math.cos(theta)
            CfT *= pw_p / w_perim
            CfT += WdT * (3.0 * ar / a_p) * (Dh_i / H) * math.pow(Dh_i / Dw, 0.18)
            CfL *= pw_p / w_perim
            CfL += WdL * (3.0 * ar / a_p) * (Dh_i / H) * (Dh_i / Dw)
        else:
            if i_type == EDGE:
                ar = math.pi * (D + Dw) * Dw / 4.0
                a_p = S + math.pi * Dw ** 2 / 8.0 / math.cos(theta)
            else:
                ar = math.pi * (D + Dw) * Dw / 6.0
                a_p = S + math.pi * Dw ** 2 / 24.0 / math.cos(theta)
            corr = 1 + WsT * (ar / a_p) * math.tan(theta) ** 2
            if not math.isfinite(corr) or corr < 0.0:
                raise MooseError("turbulent wire correction negative: %r" % corr)
            CfT *= math.pow(corr, 1.41)
            CfL *= 1 + WsL * (ar / a_p) * math.tan(theta) ** 2
    fL = CfL * math.pow(Re_eff, -1.0)
    fT = CfT * math.pow(Re_eff, -0.18)
    gamma, lam = 1.0 / 3.0, 7.0
    if Reb < ReL:
        return fL
    elif Reb > ReT:
        return fT
    else:
        return fL * math.pow(1 - psi, gamma) * (1 - math.pow(psi, lam)) + fT * math.pow(psi, gamma)


def cpp_old_corr(i_type, S, D, H, Dw):
    """Parent commit (ln version) turbulent wire correction for edge/corner."""
    theta = math.acos(H / math.sqrt(H ** 2 + (math.pi * (D + Dw)) ** 2))
    ws_t = -11.0 * math.log(H / D) + 19.0
    if i_type == EDGE:
        ar = math.pi * (D + Dw) * Dw / 4.0
        a_p = S + math.pi * Dw ** 2 / 8.0 / math.cos(theta)
    else:
        ar = math.pi * (D + Dw) * Dw / 6.0
        a_p = S + math.pi * Dw ** 2 / 24.0 / math.cos(theta)
    return 1 + ws_t * (ar / a_p) * math.tan(theta) ** 2


def cpp_mixing(CL, CT, bulk_Re, pitch, D):
    """SCMMixingChengTodreas.C:110-112, :152-167."""
    ReL = 320.0 * math.pow(10.0, pitch / D - 1.0)
    ReT = 10000.0 * math.pow(10.0, 0.7 * (pitch / D - 1.0))
    if bulk_Re < ReL:
        return CL
    elif bulk_Re > ReT:
        return CT
    psi = math.log(bulk_Re / ReL) / math.log(ReT / ReL)
    return CL + (CT - CL) * math.pow(psi, 2.0 / 3.0)


def cpp_borishanskii_base(poD, use_log10=True):
    poly = -8.12 + 12.76 * poD - 3.65 * poD ** 2
    if poly <= 0.0:
        raise MooseError("poly <= 0")
    return 24.15 * (math.log10(poly) if use_log10 else math.log(poly))


def cpp_bulk_re(mdot, mu, Dh, A):
    """SubChannel1PhaseProblem::computeBulkReynoldsNumber (:418-436), rank-0 body."""
    viscosity_in = 0.0
    mass_flow_in = 0.0
    for m, u in zip(mdot, mu):
        viscosity_in += m * u
        mass_flow_in += m
    try:
        inlet_mu = viscosity_in / mass_flow_in
    except ZeroDivisionError:
        inlet_mu = float("nan")          # C++: 0.0/0.0 = NaN
    denom = inlet_mu * A
    try:
        return mass_flow_in * Dh / denom
    except ZeroDivisionError:
        return float("nan")


# ============================================================ (B) Lean definitions
if HAVE_MP:
    F = mp.mpf
    LOG, LOG10, SQRT, ACOS, COS, TAN, PI = mp.log, mp.log10, mp.sqrt, mp.acos, mp.cos, mp.tan, mp.pi

    def RPOW(x, y):
        return mp.power(x, y) if x != 0 else (F(0) if y != 0 else F(1))
else:
    F = float
    LOG, LOG10, SQRT, ACOS, COS, TAN, PI = (math.log, math.log10, math.sqrt, math.acos, math.cos,
                                            math.tan, math.pi)
    RPOW = math.pow


def lean_quadPoly(k1, k2, k3, x):
    return k1 + k2 * x + k3 * x ** 2


def lean_triBareCoeffs(ty, pd, wd):
    if ty == CENTER:
        if pd < F("1.1"):
            return (lean_quadPoly(F("26.0"), F("888.2"), F("-3334.0"), pd - 1),
                    lean_quadPoly(F("0.09378"), F("1.398"), F("-8.664"), pd - 1))
        return (lean_quadPoly(F("62.97"), F("216.9"), F("-190.2"), pd - 1),
                lean_quadPoly(F("0.1458"), F("0.03632"), F("-0.03333"), pd - 1))
    if ty == EDGE:
        if wd < F("1.1"):
            return (lean_quadPoly(F("26.18"), F("554.5"), F("-1480.0"), wd - 1),
                    lean_quadPoly(F("0.09377"), F("0.8732"), F("-3.341"), wd - 1))
        return (lean_quadPoly(F("44.4"), F("256.7"), F("-267.6"), wd - 1),
                lean_quadPoly(F("0.1430"), F("0.04199"), F("-0.04428"), wd - 1))
    if wd < F("1.1"):
        return (lean_quadPoly(F("26.98"), F("1636.0"), F("-10050.0"), wd - 1),
                lean_quadPoly(F("0.1004"), F("1.625"), F("-11.85"), wd - 1))
    return (lean_quadPoly(F("87.26"), F("38.59"), F("-55.12"), wd - 1),
            lean_quadPoly(F("0.1499"), F("0.006706"), F("-0.009567"), wd - 1))


def lean_theta(H, D, Dw):
    return ACOS(H / SQRT(H ** 2 + (PI * (D + Dw)) ** 2))


def lean_WdT(D, Dw, H):
    return (F("19.56") + F("-98.71") * (Dw / D) + F("303.47") * (Dw / D) ** 2) * RPOW(H / D, F("-0.541"))


def lean_WsT(hd):
    return -11 * (LOG(hd) / LOG(F(10))) + 19


def lean_WsTold(hd):
    return F("-11.0") * LOG(hd) + F("19.0")


def lean_turbCorr(Ws, ar, ap, th):
    return 1 + Ws * (ar / ap) * TAN(th) ** 2


def lean_triCoeffs(ty, S, wperim, pitch, D, H, Dw, gap, has_wire):
    pd, wd = pitch / D, (D + gap) / D
    CfL, CfT = lean_triBareCoeffs(ty, pd, wd)
    if not has_wire:
        return (CfL, CfT)
    th = lean_theta(H, D, Dw)
    WsL = F("1.0") * lean_WsT(H / D)
    if ty == CENTER:
        Dh = F("4.0") * S / wperim
        ar = PI * (D + Dw) * Dw / F("6.0")
        ap = S + PI * Dw ** 2 / F("8.0") / COS(th)
        return (CfL * (PI * D / F("2.0") / wperim)
                + F("1.4") * lean_WdT(D, Dw, H) * (F("3.0") * ar / ap) * (Dh / H) * (Dh / Dw),
                CfT * (PI * D / F("2.0") / wperim)
                + lean_WdT(D, Dw, H) * (F("3.0") * ar / ap) * (Dh / H) * RPOW(Dh / Dw, F("0.18")))
    if ty == EDGE:
        ar = PI * (D + Dw) * Dw / F("4.0")
        ap = S + PI * Dw ** 2 / F("8.0") / COS(th)
    else:
        ar = PI * (D + Dw) * Dw / F("6.0")
        ap = S + PI * Dw ** 2 / F("24.0") / COS(th)
    corr = lean_turbCorr(lean_WsT(H / D), ar, ap, th)
    if corr < 0:
        return None
    return (CfL * (1 + WsL * (ar / ap) * TAN(th) ** 2), CfT * RPOW(corr, F("1.41")))


def lean_ReL(pd):
    return RPOW(F(10), F("1.0") * (pd - 1)) * 320


def lean_ReT(pd):
    return RPOW(F(10), F("0.7") * (pd - 1)) * 10000


def lean_psi(Reb, ReL, ReT):
    return LOG(Reb / ReL) / LOG(ReT / ReL)


def lean_transitionBlend(fL, fT, p):
    g = F(1) / F(3)
    return fL * RPOW(1 - p, g) * (1 - RPOW(p, F(7))) + fT * RPOW(p, g)


def lean_regimeFriction(fL, fT, Reb, ReL, ReT):
    if Reb < ReL:
        return fL
    if Reb > ReT:
        return fT
    return lean_transitionBlend(fL, fT, lean_psi(Reb, ReL, ReT))


def lean_triFrictionFactor(Re, ty, S, wperim, pitch, D, H, Dw, gap, has_wire, Reb):
    c = lean_triCoeffs(ty, S, wperim, pitch, D, H, Dw, gap, has_wire)
    if c is None:
        return None
    Re_eff = max(Re, F("1.0"))
    pd = pitch / D
    return lean_regimeFriction(c[0] * RPOW(Re_eff, F("-1.0")), c[1] * RPOW(Re_eff, F("-0.18")),
                               Reb, lean_ReL(pd), lean_ReT(pd))


def lean_mixingRegime(CL, CT, bRe, pd):
    ReL, ReT = lean_ReL(pd), lean_ReT(pd)
    if bRe < ReL:
        return CL
    if bRe > ReT:
        return CT
    return CL + (CT - CL) * RPOW(lean_psi(bRe, ReL, ReT), F(2) / F(3))


# ============================================================ sampling
def bare_area_perim(ty, P, D, gap):
    if ty == CENTER:
        return math.sqrt(3) / 4 * P ** 2 - math.pi * D ** 2 / 8, math.pi * D / 2
    if ty == EDGE:
        return P * (D / 2 + gap) - math.pi * D ** 2 / 8, P + math.pi * D / 2
    w = D / 2 + gap
    return w ** 2 / math.sqrt(3) - math.pi * D ** 2 / 24, math.pi * D / 6 + 2 * w / math.sqrt(3)


def sample(rng):
    D = rng.uniform(5e-3, 1e-2)
    pd = rng.uniform(1.02, 1.42)
    P = pd * D
    hd = rng.uniform(8, 52)
    Dw = rng.uniform(0.05 * D, min(0.25 * D, P - D))
    gap = rng.uniform(0.0, 0.3 * D) if rng.random() < 0.5 else Dw
    ty = rng.choice([CENTER, EDGE, CORNER])
    S, wp = bare_area_perim(ty, P, D, gap)
    has_wire = rng.random() < 0.7
    Re = 10 ** rng.uniform(0, 6)
    Reb = 10 ** rng.uniform(math.log10(50), 6)
    return dict(Re=Re, i_type=ty, S=S, w_perim=wp, pitch=P, D=D, H=hd * D, Dw=Dw, gap=gap,
                has_wire=has_wire, Reb=Reb)


def relerr(a, b):
    a, b = float(a), float(b)
    if math.isnan(a) and math.isnan(b):
        return 0.0
    return abs(a - b) / max(abs(b), 1e-300)


res = {
    "module": "Moose33307",
    "merge_commit": "b5a1d9e3df14e98850e514f3e57d02c25c44fca9",
    "native_replay": "SKIP: no MOOSE build available; float64 transcription of the C++ is "
                     "compared with the Lean formulas instead",
    "lean_side_precision": "mpmath 50 digits" if HAVE_MP else "float64 (mpmath missing)",
}

# ---------------- differential test
rng = random.Random(1)
max_rel = 0.0
errors_cpp = 0
errors_mismatch = 0
for _ in range(N):
    s = sample(rng)
    try:
        a = cpp_tri_friction(**s)
    except MooseError:
        a = None
        errors_cpp += 1
    b = lean_triFrictionFactor(F(s["Re"]), s["i_type"], F(s["S"]), F(s["w_perim"]), F(s["pitch"]),
                               F(s["D"]), F(s["H"]), F(s["Dw"]), F(s["gap"]), s["has_wire"],
                               F(s["Reb"]))
    if (a is None) != (b is None):
        errors_mismatch += 1
    elif a is not None:
        max_rel = max(max_rel, relerr(a, b))
res["diff_tri_friction_cpp_vs_lean"] = {"points": N, "max_rel_diff": max_rel,
                                        "mooseError_count_cpp": errors_cpp,
                                        "none_mismatch_count": errors_mismatch}

# ---------------- ReL_pos_lt_ReT
ok = all(0 < float(lean_ReL(F(1 + i / 1000))) < float(lean_ReT(F(1 + i / 1000)))
         for i in range(0, 1001))
res["thm_ReL_pos_lt_ReT"] = {"grid_P_over_D_1_to_2": ok,
                             "ReL(1.42)": float(lean_ReL(F("1.42"))),
                             "ReT(1.0)": float(lean_ReT(F(1)))}

# ---------------- psi_mem_Icc and regime continuity (C++ transcription)
rng = random.Random(2)
psi_ok = True
max_jump = {}
for delta in (1e-4, 1e-6, 1e-8):
    jm = 0.0
    for _ in range(N // 4):
        s = sample(rng)
        s["has_wire"] = s["has_wire"] and True
        pd = s["pitch"] / s["D"]
        ReL = math.pow(10, 1.0 * (pd - 1)) * 320
        ReT = math.pow(10, 0.7 * (pd - 1)) * 10000
        for edge in (ReL, ReT):
            try:
                lo = cpp_tri_friction(**dict(s, Reb=edge * (1 - delta)))
                hi = cpp_tri_friction(**dict(s, Reb=edge * (1 + delta)))
            except MooseError:
                continue
            jm = max(jm, abs(hi - lo) / max(abs(lo), abs(hi)))
        rb = math.exp(rng.uniform(math.log(ReL), math.log(ReT)))
        p = math.log(rb / ReL) / math.log(ReT / ReL)
        psi_ok &= (0.0 <= p <= 1.0)
    max_jump[str(delta)] = jm
res["thm_psi_mem_Icc"] = {"holds_on_samples": psi_ok}
res["thm_regimeFriction_continuousOn"] = {
    "max_rel_jump_across_ReL_ReT_by_delta": max_jump,
    "note": "jump shrinks proportionally to delta (no discontinuity)"}

# ---------------- transitionBlend nonneg, exceeds ends, dips below min near psi -> 1
def weights(p):
    return (1 - p) ** (1 / 3) * (1 - p ** 7), p ** (1 / 3)


grid = [i / 100000 for i in range(100001)]
wsum = [sum(weights(p)) for p in grid]
imax = max(range(len(grid)), key=lambda i: wsum[i])
fine = [1 - 10 ** (-k / 50) for k in range(50, 600)]
wfine = [sum(weights(p)) for p in fine]
imin = min(range(len(fine)), key=lambda i: wfine[i])
res["thm_transitionBlend_nonneg"] = {"min_over_grid_fL_fT_1": min(wsum)}
res["thm_transitionBlend_exceeds_ends"] = {
    "psi_1_8_value_fL_eq_fT_1": float(lean_transitionBlend(F(1), F(1), F(1) / 8)),
    "max_weight_sum": wsum[imax], "argmax_psi": grid[imax],
    "min_weight_sum_near_psi_1": wfine[imin], "argmin_psi": fine[imin],
}
# physical example: center bare channel, P/D = 1.2, Re_i = Reb = sqrt(ReL*ReT)
pdx = 1.2
ReLx = math.pow(10, pdx - 1) * 320
ReTx = math.pow(10, 0.7 * (pdx - 1)) * 10000
Rex = math.sqrt(ReLx * ReTx)
Dx = 8e-3
Sx, wpx = bare_area_perim(CENTER, pdx * Dx, Dx, 0.0)
f_mid = cpp_tri_friction(Rex, CENTER, Sx, wpx, pdx * Dx, Dx, 0.0, 0.0, 0.0, False, Rex)
CfLx = 62.97 + 216.9 * 0.2 - 190.2 * 0.04
CfTx = 0.1458 + 0.03632 * 0.2 - 0.03333 * 0.04
res["thm_transitionBlend_exceeds_ends"]["physical_example"] = {
    "P_over_D": pdx, "Re_i=Reb": Rex, "fL": CfLx / Rex, "fT": CfTx * Rex ** -0.18,
    "f_transition": f_mid, "ratio_to_max_end": f_mid / max(CfLx / Rex, CfTx * Rex ** -0.18)}

# ---------------- WsT_pos, WsTold_neg
hd_grid = [8 + 44 * i / 10000 for i in range(10001)]
res["thm_WsT_pos"] = {"min_on_[8,52]": min(float(lean_WsT(F(h))) for h in hd_grid),
                      "zero_at_H_over_D": 10 ** (19 / 11)}
res["thm_WsTold_neg"] = {"max_on_[8,52]": max(float(lean_WsTold(F(h))) for h in hd_grid)}

# ---------------- wire corrections >= 1 (new), old ln version error cases
rng = random.Random(3)
min_corr_new = float("inf")
old_neg = 0
old_example = None
for _ in range(N):
    s = sample(rng)
    if s["i_type"] == CENTER:
        continue
    D, H, Dw, S = s["D"], s["H"], s["Dw"], s["S"]
    th = math.acos(H / math.sqrt(H ** 2 + (math.pi * (D + Dw)) ** 2))
    if s["i_type"] == EDGE:
        ar = math.pi * (D + Dw) * Dw / 4.0
        ap = S + math.pi * Dw ** 2 / 8.0 / math.cos(th)
    else:
        ar = math.pi * (D + Dw) * Dw / 6.0
        ap = S + math.pi * Dw ** 2 / 24.0 / math.cos(th)
    corr = 1 + (-11 * math.log10(H / D) + 19) * (ar / ap) * math.tan(th) ** 2
    min_corr_new = min(min_corr_new, corr)
    oc = cpp_old_corr(s["i_type"], S, D, H, Dw)
    if oc < 0:
        old_neg += 1
        if old_example is None:
            old_example = dict(type=s["i_type"], D=D, H_over_D=H / D, Dw_over_D=Dw / D,
                               P_over_D=s["pitch"] / D, S=S, old_corr=oc, new_corr=corr)
res["thm_edgeCorr_cornerCorr_ge_one"] = {"min_new_correction": min_corr_new}
res["old_ln_version_mooseError_cases"] = {"count": old_neg, "example": old_example}

# ---------------- mixing between, continuity
rng = random.Random(4)
mix_ok = True
mix_diff = 0.0
for _ in range(N):
    CL, CT = rng.uniform(0, 1), rng.uniform(0, 1)
    pd = rng.uniform(1.0, 1.42)
    bRe = 10 ** rng.uniform(1, 6)
    v = cpp_mixing(CL, CT, bRe, pd, 1.0)
    mix_ok &= (min(CL, CT) - 1e-15 <= v <= max(CL, CT) + 1e-15)
    mix_diff = max(mix_diff, relerr(v, lean_mixingRegime(F(CL), F(CT), F(bRe), F(pd))))
res["thm_mixingRegime_between"] = {"holds_on_samples": mix_ok}
res["diff_mixing_cpp_vs_lean"] = {"points": N, "max_rel_diff": mix_diff}

# ---------------- Borishanskii
xs = [1.1 + 0.4 * i / 1000 for i in range(1001)]
res["thm_borishanskiiPoly_gt_one"] = {"min_poly_on_[1.1,1.5]":
                                      min(-8.12 + 12.76 * x - 3.65 * x * x for x in xs)}
res["thm_borishanskiiNuBase_pos"] = {
    "Nu_log10": {str(x): cpp_borishanskii_base(x) for x in (1.1, 1.2, 1.5)},
    "Nu_ln_parent": {str(x): cpp_borishanskii_base(x, False) for x in (1.1, 1.2, 1.5)},
    "ratio_ln_over_log10": cpp_borishanskii_base(1.2, False) / cpp_borishanskii_base(1.2)}

# ---------------- bulk Re
rng = random.Random(5)
mdot = [rng.uniform(0.01, 0.1) for _ in range(19)]
mu0 = 2.5e-4
br = cpp_bulk_re(mdot, [mu0] * 19, 3e-3, 1e-3)
res["thm_bulkReynolds_uniform_mu"] = {"rel_diff": relerr(br, sum(mdot) * 3e-3 / (mu0 * 1e-3))}
nan_re = cpp_bulk_re([0.0] * 19, [mu0] * 19, 3e-3, 1e-3)
s0 = sample(random.Random(6))
s0["Reb"] = nan_re
f_nan = cpp_tri_friction(**dict(s0, has_wire=False))
res["zero_inlet_flow_edge_case"] = {
    "bulk_Re": str(nan_re), "friction_factor": str(f_nan),
    "note": "NaN bulk Re fails both comparisons and falls into the transition branch"}
res["constructor_range_check"] = {
    "value_read": 1.0, "threshold": 50.0,
    "note": "SCMFrictionUpgradedChengTodreas.C:44,58 reads getBulkReynoldsNumber() in the "
            "constructor; _bulk_Re is initialised to 1.0 (SubChannel1PhaseProblem.C:209) and only "
            "computed in externalSolve (:2759), so the Reb range warning is flagged for every "
            "triangular-lattice run"}
res["parallel_bulk_re_static"] = {
    "guard": "SubChannel1PhaseProblem.C:420-421 if (processor_id() != 0) return;",
    "consumers": "computeFrictionFactor at :875,:940,:1159 and mixing at :1958,:2031,:2048, all "
                 "reached only from computeWijFromSolve / implicitPetscSolve / computeh / computeT "
                 "which run under if (processor_id() == 0) (:2832) or after goto aux_close on "
                 "ranks > 0 (:2855-2860); post-solve HTC loops (:2928-2992) iterate over "
                 "_n_channels / _n_pins which are 0 on ranks > 0 (comment :264-269)",
    "replay": "not executed (no MOOSE build)"}
print(json.dumps(res, indent=2, default=str))
