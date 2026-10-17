"""MOOSE#33307 (subchannel closure relations switched to log10, transition region switched to component Re_b, etc.) unified ΔQ1 spot check.
Non-native: there is no MOOSE build on the node; the modified functions are recomputed line by line in scm_port.py from the C++ source of merge commit b5a1d9e;
checks that need the full solver (component pressure drop, agreement of the three solution strategies, MPI, cross-flow estimate indexing) are listed in not_covered and not implemented.
Usage: python spotcheck.py --track {1,2,3} --n N --out PATH (pure Python + numpy, any Python 3.10+ environment)
Plan mapping and implementation judgments are in impl_log.md.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spotlib as SL  # noqa: E402
import scm_port as M  # noqa: E402
from scm_port import CENTER, EDGE, CORNER, MooseError  # noqa: E402

PR = "idaholab/moose#33307"
PR_IDX = 9
TYPES = (CENTER, EDGE, CORNER)
# template test19_implicit.i: nrings=3, pitch=0.012, D=0.01, dwire=0.002, hwire=0.0833, flat_to_flat=0.056
T_D = 0.01
T_PITCH = 0.012
T_GAP = 0.5 * (0.056 - (3 - 1) * 0.012 * math.sqrt(3.0) - 0.01)
NPIN_TO_NR = {1 + 3 * nr * (nr - 1): nr for nr in range(2, 10)}
MU_REF = 2.5e-4  # the local field enters only the warnings of the mixing closure; a sodium-order constant is used (implementation judgment)

NOT_COVERED = {
    1: ["A component pressure drop ΔP_wire≥ΔP_bare (needs the full solver; replaced by closure-level f_wire≥f_bare and correction factor ≥1)",
        "C back out Nu from wall temperature / bulk temperature (read the Nu returned by the closure directly)",
        "D explicit/implicit/monolithic solution strategies agree (needs the full solver)",
        "risk item: 1 process vs 4 processes agree (needs MPI)"],
    2: ["I1 component pressure drop and runtime-error cross-check (replaced by closure-level f_wire≥f_bare, correction factor ≥1, no mooseError triggered)",
        "I3 converged solutions of the three solution strategies agree (needs the full solver)", "I3 single process vs mpiexec -n 2 agree (needs MPI)",
        "water (Pr range not given; Borishanskii is used only for liquid metals)"],
    3: ["I1 whether bulk Re is stale under transient inlet changes, serial vs 2 processes agree (needs solver/MPI)",
        "I4 cross-flow estimate correlation matrix and mass conservation (needs the PETSc solver; the plan itself states that the public API cannot decide it)",
        "I6 single-subchannel zero-flow _max_sumWij/min_mdot in the implicit coupled solve (needs the solver)",
        "full mesh generation beyond I5 (only the three nrings/pitch/dwire checks are recomputed)"],
}


def versions():
    return dict(moose_merge_commit="b5a1d9e3df14e98850e514f3e57d02c25c44fca9", native=False,
                recompute="scm_port.py (line-by-line recomputation of SCMFrictionUpgradedChengTodreas / SCMMixingChengTodreas / "
                          "SCMHTCBorishanskii / SCMHTCGnielinski / SCMHTCClosureBase / computeBulkReynoldsNumber / "
                          "SCMTriAssemblyMeshGenerator checks / Tri,QuadSubChannelMesh area and wetted perimeter)",
                numpy=np.__version__, python=sys.version.split()[0])


def tri(pod, hod, dw, nr, D=T_D, gap=None, check=True):
    gap = T_GAP / T_D * D if gap is None else gap
    return M.TriGeom(pod * D, D, dw, hod * D if hod else 0.0, nr, gap, check_mesh=check)


def re_local(g, t, Reb):
    """Uniform mass-flux assumption: Re_i = Re_b·Dh_i/Dh_b (implementation judgment)."""
    _, _, Dhb = g.assembly()
    return Reb * g.Dh(t) / Dhb


def same_mdot_Reb(g_from, g_to, Reb):
    """Re_b after changing geometry at the same inlet mass flow and viscosity: Re = ṁ·Dh/(μ·A)."""
    A0, _, Dh0 = g_from.assembly()
    A1, _, Dh1 = g_to.assembly()
    return Reb * (A0 / Dh0) * (Dh1 / A1)


def fric(g, t, Re, Reb, info=None):
    return M.friction_tri(g, t, Re, g.flow_area(t), g.wetted_perimeter(t), Reb, info)


def gap_types(branch_type):
    return {CENTER: (CENTER, CENTER), EDGE: (EDGE, EDGE), CORNER: (CORNER, EDGE)}[branch_type]


def local_fields(g, ti, tj, Reb, scale_m=1.0, scale_mu=1.0):
    Si, Sj = g.flow_area(ti), g.flow_area(tj)
    _, _, Dhb = g.assembly()
    G = Reb * MU_REF / Dhb * scale_m
    mu = MU_REF * scale_mu
    return dict(S=(Si, Sj, Si, Sj), mu=(mu, mu, mu, mu), mdot=(G * Si, G * Sj, G * Si, G * Sj),
                wp=(g.wetted_perimeter(ti), g.wetted_perimeter(tj)))


def mix(g, branch_type, Reb, info=None, **kw):
    ti, tj = gap_types(branch_type)
    loc = local_fields(g, ti, tj, Reb, **kw)
    if branch_type == CENTER:
        return M.mixing_beta(g, ti, tj, Reb, loc, info)
    return M.sweep_beta(g, ti, tj, Reb, loc, info)


def jump(f_plus, f_minus):
    return abs(f_plus - f_minus) / abs(f_minus) if f_minus != 0 else abs(f_plus - f_minus)


def borishanskii_formula(pod):
    return 24.15 * math.log10(-8.12 + 12.76 * pod - 3.65 * pod ** 2)


# ================================================================ Track 1
def t1_sample(t):
    def s(rng):
        npin = int(rng.choice([7, 19, 37, 61, 217]))
        return dict(type=t, PoD=float(rng.uniform(1.05, 1.42)), HoD=float(rng.uniform(8, 52)),
                    DwoD=float(rng.uniform(0.05, 0.25)), Npin=npin, Reb=float(rng.uniform(50, 1e6)),
                    Pr=float(rng.uniform(0.003, 0.03)), Pe=float(rng.uniform(10, 2200)))
    return s


def t1_eval(inp):
    t = inp["type"]
    nr = NPIN_TO_NR[inp["Npin"]]
    g = tri(inp["PoD"], inp["HoD"], inp["DwoD"] * T_D, nr)  # mesh-generator check fails -> MooseError -> fail
    Reb = inp["Reb"]
    M.friction_constructor_tri(g, Reb)
    Rei = re_local(g, t, Reb)
    info = {}
    f = fric(g, t, Rei, Reb, info)
    ch = []
    # ---- A
    if t in (EDGE, CORNER):
        ch.append(SL.cond("t1_A_WsT>0", info["WsT"] > 0, observed=info["WsT"], expected=">0"))
        ch.append(SL.cond("t1_A_turbulent_wire_factor>=1", info["turbulent_wire_correction"] >= 1,
                          observed=info["turbulent_wire_correction"], expected=">=1"))
        ch.append(SL.cond("t1_A_laminar_wire_factor>=1", info["laminar_wire_correction"] >= 1,
                          observed=info["laminar_wire_correction"], expected=">=1"))
        gb = tri(inp["PoD"], inp["HoD"], 1e-6, nr)
        Rebb = same_mdot_Reb(g, gb, Reb)
        fb = fric(gb, t, re_local(gb, t, Rebb), Rebb)
        ch.append(SL.cond("t1_A_f_wire>=f_nearbare(Dw=1e-6m)_same_mdot", f >= fb, observed=f, expected=f">= {fb}"))
        fs = []
        for hod in (8.0, 20.0, 52.0):
            gm = tri(inp["PoD"], hod, inp["DwoD"] * T_D, nr)
            Rm = same_mdot_Reb(g, gm, Reb)
            fs.append(fric(gm, t, re_local(gm, t, Rm), Rm))
        ch.append(SL.cond("t1_A_f_monotone_increasing_as_HoD_decreases(8,20,52)", fs[0] >= fs[1] >= fs[2],
                          observed=fs))
    # ---- B
    for lab, R0, fend in (("ReL", info["ReL"], "fL"), ("ReT", info["ReT"], "fT")):
        # post-smoke-test fix: the fixed Re_i takes the uniform-mass-flux value at the Re_b=R0 under test (originally mistakenly used Re_i at the sampled Re_b)
        Rei0 = re_local(g, t, R0)
        ie = {}
        fe = fric(g, t, Rei0, R0, ie)
        ch.append(SL.close(f"t1_B_f(Reb={lab})=={fend}", fe, ie[fend]))
        for d in (1e-3, 1e-4):
            j = jump(fric(g, t, Rei0, R0 * (1 + d)), fric(g, t, Rei0, R0 * (1 - d)))
            ch.append(SL.cond(f"t1_B_f_jump_at_{lab}_delta={d:g}<10delta", j < 10 * d, observed=j, expected=f"<{10 * d:g}"))
        ic = {}
        ce = mix(g, t, R0, ic)
        Cname = "Cm" if t == CENTER else "Cs"
        ch.append(SL.close(f"t1_B_{Cname}(Reb={lab})==endpoint", ic["C"], ic["CL" if lab == "ReL" else "CT"]))
        for d in (1e-3, 1e-4):
            j = jump(mix(g, t, R0 * (1 + d)), mix(g, t, R0 * (1 - d)))
            ch.append(SL.cond(f"t1_B_{Cname}_jump_at_{lab}_delta={d:g}<10delta", j < 10 * d, observed=j,
                              expected=f"<{10 * d:g}"))
        del ce
    psis = []
    for tt in TYPES:
        ii = {}
        fric(g, tt, re_local(g, tt, Reb), Reb, ii)
        psis.append(ii["psi"])
    ch.append(SL.cond("t1_B_same_psi_all_channel_types", psis[0] == psis[1] == psis[2], observed=psis))
    # ---- C
    if inp["Pe"] <= 200.0:
        pre = M.nusselt_preinfo(g.pitch, g.D, t, inp["Pe"] / inp["Pr"], inp["Pr"])
        Nu = M.borishanskii_Nu(pre)
        ch.append(SL.close("t1_C_Borishanskii_Pe<=200_Nu==24.15log10(poly)", Nu, borishanskii_formula(inp["PoD"])))
        if 1.1 <= inp["PoD"] <= 1.5:
            ch.append(SL.cond("t1_C_Nu_in_[4.25,10.83]", 4.25 <= Nu <= 10.83, observed=Nu, expected=[4.25, 10.83]))
    # ---- E
    if Rei > 1e3:
        gi = {}
        pre = M.nusselt_preinfo(g.pitch, g.D, t, Rei, inp["Pr"])
        Nu = M.gnielinski_Nu(pre, f, gi)
        ch.append(SL.cond("t1_E_Gnielinski_NuT_finite>0", math.isfinite(gi["NuT"]) and gi["NuT"] > 0,
                          observed=gi["NuT"], expected=">0"))
        ch.append(SL.cond("t1_E_Gnielinski_Nu_finite>0", math.isfinite(Nu) and Nu > 0, observed=Nu, expected=">0"))
        for d in (1e-3, 1e-4):
            Nup = M.gnielinski_Nu(M.nusselt_preinfo(g.pitch, g.D, t, Rei, inp["Pr"] * (1 + d)), f)
            Num = M.gnielinski_Nu(M.nusselt_preinfo(g.pitch, g.D, t, Rei, inp["Pr"] * (1 - d)), f)
            j = jump(Nup, Num)
            ch.append(SL.cond(f"t1_E_Nu_continuous_in_Pr_delta={d:g}<10delta", j < 10 * d, observed=j,
                              expected=f"<{10 * d:g}"))
    return ch


# ================================================================ Track 2
def t2_sample(t):
    def s(rng):
        return dict(type=t, PoD=float(rng.uniform(1.0, 1.42)), HoD=float(rng.uniform(8, 52)),
                    Npin=int(rng.choice([7, 19, 37, 61, 91, 127, 169, 217])), Reb=float(rng.uniform(50, 1e6)),
                    Pr=float(rng.uniform(0.004, 0.03)), PoD_B=float(rng.uniform(1.1, 1.5)),
                    Pe=float(rng.uniform(200, 2000)))
    return s


def t2_eval(inp):
    t = inp["type"]
    nr = NPIN_TO_NR[inp["Npin"]]
    pod = inp["PoD"]
    g = tri(pod, inp["HoD"], (pod - 1.0) * T_D, nr)  # D_w = P − D (relation from template test19)
    Reb = inp["Reb"]
    M.friction_constructor_tri(g, Reb)
    Rei = re_local(g, t, Reb)
    ch = []
    info = {}
    try:
        f = fric(g, t, Rei, Reb, info)
        err = None
    except MooseError as e:
        f, err = float("nan"), str(e)
    ch.append(SL.cond("t2_I1_no_turbulent_wire_correction_error", err is None, expected="no mooseError",
                      observed=err or "none"))
    if err is not None:
        return ch
    if info["has_wire"]:
        ch.append(SL.cond("t2_I1_WsT_in_[0.12,9.07]", 0.12 <= info["WsT"] <= 9.07, observed=info["WsT"],
                          expected=[0.12, 9.07]))
        ch.append(SL.cond("t2_I1_WsT>0", info["WsT"] > 0, observed=info["WsT"], expected=">0"))
        gb = tri(pod, 0.0, 0.0, nr)  # no wire: dwire=hwire=0 (setting only dwire=0 would trigger a constructor error; implementation judgment)
        fb = fric(gb, t, Rei, Reb)
        if t in (EDGE, CORNER):
            ch.append(SL.cond("t2_I1_turbulent_wire_factor>=1", info["turbulent_wire_correction"] >= 1,
                              observed=info["turbulent_wire_correction"], expected=">=1"))
            ch.append(SL.cond("t2_I1_laminar_wire_factor>=1", info["laminar_wire_correction"] >= 1,
                              observed=info["laminar_wire_correction"], expected=">=1"))
            ch.append(SL.cond("t2_I1_f_wire>=f_bare_same_Re", f >= fb, observed=f, expected=f">= {fb}"))
        else:
            ch.append(SL.cond("t2_I1_center_drag_terms>0", info["dragT"] > 0 and info["dragL"] > 0,
                              observed=[info["dragT"], info["dragL"]], expected=">0"))
            ch.append(SL.close("t2_I1_center_CfT==CfT_bare*Pwp/Pw+drag", info["CfT"],
                               info["CfT_bare"] * info["pw_ratio"] + info["dragT"]))
    # ---- I2
    ch.append(SL.close("t2_I2_ReL==320*10^(P/D-1)", info["ReL"], 320 * 10 ** (pod - 1)))
    ch.append(SL.close("t2_I2_ReT==1e4*10^(0.7(P/D-1))", info["ReT"], 1e4 * 10 ** (0.7 * (pod - 1))))
    ch.append(SL.cond("t2_I2_f>=0", f >= 0, observed=f, expected=">=0"))
    for lab, R0, fend in (("ReL", info["ReL"], "fL"), ("ReT", info["ReT"], "fT")):
        # post-smoke-test fix: the fixed Re_i takes the uniform-mass-flux value at the Re_b=R0 under test (originally mistakenly used Re_i at the sampled Re_b)
        Rei0 = re_local(g, t, R0)
        ie = {}
        fe = fric(g, t, Rei0, R0, ie)
        ch.append(SL.close(f"t2_I2_f(Reb={lab})=={fend}", fe, ie[fend]))
        j = jump(fric(g, t, Rei0, R0 * (1 + 1e-3)), fric(g, t, Rei0, R0 * (1 - 1e-3)))
        ch.append(SL.cond(f"t2_I2_f_jump_at_{lab}_pm0.1pct<=1e-3", j <= 1e-3, observed=j, expected="<=1e-3"))
        if info["has_wire"]:
            ic = {}
            mix(g, t, R0, ic)
            Cname = "Cm" if t == CENTER else "Cs"
            ch.append(SL.close(f"t2_I2_{Cname}(Reb={lab})==endpoint", ic["C"], ic["CL" if lab == "ReL" else "CT"]))
            j = jump(mix(g, t, R0 * (1 + 1e-3)), mix(g, t, R0 * (1 - 1e-3)))
            ch.append(SL.cond(f"t2_I2_{Cname}_jump_at_{lab}_pm0.1pct<=1e-3", j <= 1e-3, observed=j, expected="<=1e-3"))
    psis = []
    for tt in TYPES:
        ii = {}
        fric(g, tt, re_local(g, tt, Reb), Reb, ii)
        psis.append(ii["psi"])
    ch.append(SL.cond("t2_I2_same_psi_all_channel_types", psis[0] == psis[1] == psis[2], observed=psis))
    # ---- I3 (Borishanskii; Pe≤200 applies only at the lower end of this track's Pe interval [200,2000])
    pre = M.nusselt_preinfo(inp["PoD_B"] * T_D, T_D, t, 200.0 / inp["Pr"], inp["Pr"])
    Nu = M.borishanskii_Nu(pre)
    ch.append(SL.close("t2_I3_Borishanskii_Pe=200_Nu==24.15log10(poly)", Nu, borishanskii_formula(inp["PoD_B"])))
    return ch


def t2_boris_value_eval(inp):
    pre = M.nusselt_preinfo(inp["PoD_B"] * T_D, T_D, inp["type"], 200.0 / inp["Pr"], inp["Pr"])
    Nu = M.borishanskii_Nu(pre)
    return [SL.close(f"t2_I3_Borishanskii_P/D={inp['PoD_B']}_approx", Nu, inp["approx"], rtol=0.0, atol=inp["atol"])]


# ================================================================ Track 3: independent oracle (written separately from the documented formulas + the source coefficient tables)
TRI_TAB = {  # (switch_var, [(kL1,kL2,kL3,kT1,kT2,kT3) for x<1.1, >=1.1])
    CENTER: ("pod", [(26.0, 888.2, -3334.0, 0.09378, 1.398, -8.664), (62.97, 216.9, -190.2, 0.1458, 0.03632, -0.03333)]),
    EDGE: ("wod", [(26.18, 554.5, -1480.0, 0.09377, 0.8732, -3.341), (44.4, 256.7, -267.6, 0.1430, 0.04199, -0.04428)]),
    CORNER: ("wod", [(26.98, 1636.0, -10050.0, 0.1004, 1.625, -11.85), (87.26, 38.59, -55.12, 0.1499, 0.006706, -0.009567)]),
}
QUAD_TAB = {
    CENTER: ("pod", [(26.37, 374.2, -493.9, 0.09423, 0.5806, -1.239), (35.55, 263.7, -190.2, 0.1339, 0.09059, -0.09926)]),
    EDGE: ("wod", [(26.18, 554.5, -1480, 0.09377, 0.8732, -3.341), (44.40, 256.7, -267.6, 0.1430, 0.04199, -0.04428)]),
    CORNER: ("wod", [(28.62, 715.9, -2807, 0.09755, 1.127, -6.304), (58.83, 160.7, -203.5, 0.1452, 0.02681, -0.03411)]),
}


def o_limits(pod):
    return 320.0 * 10.0 ** (pod - 1.0), 1e4 * 10.0 ** (0.7 * (pod - 1.0))


def o_interp(fl, ft, Reb, ReL, ReT):
    if Reb < ReL:
        return fl
    if Reb > ReT:
        return ft
    psi = math.log(Reb / ReL) / math.log(ReT / ReL)
    return fl * (1 - psi) ** (1 / 3) * (1 - psi ** 7) + ft * psi ** (1 / 3)


def o_poly(k, x):
    return np.polyval([k[2], k[1], k[0]], x), np.polyval([k[5], k[4], k[3]], x)


def o_friction_tri(P, D, Dw, H, gap, t, Re, S, Pw, Reb):
    pod, wod = P / D, (D + gap) / D
    sw, tabs = TRI_TAB[t]
    xv = pod if sw == "pod" else wod
    CL, CT = o_poly(tabs[0] if xv < 1.1 else tabs[1], xv - 1)
    if Dw != 0 and H != 0:
        th = math.atan(math.pi * (D + Dw) / H)  # equivalent form of acos(H/sqrt(H²+(π(D+Dw))²))
        dw = Dw / D
        WdT = (19.56 - 98.71 * dw + 303.47 * dw ** 2) * (H / D) ** -0.541
        WsT = 19.0 - 11.0 * math.log10(H / D)
        Dhi = 4 * S / Pw
        if t == CENTER:
            ar = math.pi * (D + Dw) * Dw / 6
            ap = S + math.pi * Dw ** 2 / 8 / math.cos(th)
            CT = CT * (math.pi * D / 2) / Pw + WdT * 3 * ar / ap * (Dhi / H) * (Dhi / Dw) ** 0.18
            CL = CL * (math.pi * D / 2) / Pw + 1.4 * WdT * 3 * ar / ap * (Dhi / H) * (Dhi / Dw)
        else:
            ar = math.pi * (D + Dw) * Dw / (4 if t == EDGE else 6)
            ap = S + math.pi * Dw ** 2 / (8 if t == EDGE else 24) / math.cos(th)
            fac = WsT * ar / ap * math.tan(th) ** 2
            CT = CT * (1 + fac) ** 1.41
            CL = CL * (1 + fac)
    ReL, ReT = o_limits(pod)
    Re = max(Re, 1.0)
    return o_interp(CL / Re, CT * Re ** -0.18, Reb, ReL, ReT)


def o_friction_quad(P, D, sg, t, Re, Reb):
    pod = P / D
    wod = (D / 2 + P / 2 + sg) / D
    sw, tabs = QUAD_TAB[t]
    xv = pod if sw == "pod" else wod
    CL, CT = o_poly(tabs[0] if pod < 1.1 else tabs[1], xv - 1)  # for the square lattice, all three types select the table by P/D
    ReL, ReT = o_limits(pod)
    Re = max(Re, 1.0)
    return o_interp(CL / Re, CT * Re ** -0.18, Reb, ReL, ReT)


def o_mixing(P, D, Dw, H, gap, Nr, kind, Reb):
    th = math.atan(math.pi * (D + Dw) / H)
    ReL, ReT = o_limits(P / D)
    if kind == "Cm":
        cl, ct = (0.055, 0.1) if Nr == 1 else (0.077, 0.14)
        CL, CT = cl * ((P - D) / D) ** -0.5, ct * ((P - D) / D) ** -0.5
        geo = math.sqrt((math.pi * (D + Dw) * Dw / 6) / (math.sqrt(3) / 4 * P ** 2 - math.pi * D ** 2 / 8)) * math.tan(th)
    else:
        cl, ct = (0.33, 0.6) if Nr == 1 else (0.413, 0.75)
        CL, CT = cl * (H / D) ** 0.3, ct * (H / D) ** 0.3
        geo = math.sqrt((math.pi * (D + Dw) * Dw / 4) / (P * (D + gap - D / 2) - math.pi * D ** 2 / 8)) * math.tan(th)
    if Reb < ReL:
        C = CL
    elif Reb > ReT:
        C = CT
    else:
        C = CL + (CT - CL) * (math.log(Reb / ReL) / math.log(ReT / ReL)) ** (2 / 3)
    return C * geo, C, CL, CT


def o_blend(Re, ReL, ReT, NuL, NuT):
    if Re <= ReL:
        return NuL
    if Re >= ReT:
        return NuT
    w = (Re - ReL) / (ReT - ReL)
    return w * NuT + (1 - w) * NuL


NU_LAM = {CENTER: 3.73, EDGE: 3.59, CORNER: 3.52}


def o_borishanskii(pod, t, Re, Pr):
    ReL, ReT = o_limits(pod)
    PeT = (ReT if ReL < Re < ReT else Re) * Pr
    NuT = 24.15 * math.log10(-8.12 + 12.76 * pod - 3.65 * pod ** 2)
    if PeT >= 200:
        NuT += 0.0174 * (1 - math.exp(6 - 6 * pod)) * (PeT - 200) ** 0.9
    return o_blend(Re, ReL, ReT, NU_LAM[t], NuT)


def o_gnielinski(pod, t, Re, Pr, f):
    ReL, ReT = o_limits(pod)
    p = Pr + 0.01
    d = max(1 + 12.7 * math.sqrt(f / 8) * (p ** (2 / 3) - 1), 0.1)
    NuT = (f / 8) * (Re - 1000) * p / d
    return o_blend(Re, ReL, ReT, NU_LAM[t], NuT)


TOL3 = dict(rtol=1e-10, atol=1e-12)


def o3(name, obs, exp):
    return SL.close(name, obs, exp, **TOL3)


def nr_sample(rng, choices=(2, 3, 4, 5, 6, 7, 8, 9)):
    return int(rng.choice(list(choices)))


def t3F_sample(lat, t):
    def s(rng):
        d = dict(branch="F", lattice=lat, type=t, PoD=float(rng.uniform(1.01, 1.4)), D=float(rng.uniform(0.005, 0.015)),
                 Reb=float(rng.uniform(50, 1e6)), Rei=float(rng.uniform(50, 1e6)))
        if lat != "quad":
            d.update(WoD=float(rng.uniform(1.01, 1.4)), Nr=nr_sample(rng))
        if lat == "tri_wire":
            d.update(HoD=float(rng.uniform(8, 52)), dw_frac=float(1.0 - rng.uniform(0, 1)))  # (0,1]
        if lat == "tri_bare":
            d.update(HoD=0.0, dw_frac=0.0)
        return d
    return s


def geom3(inp):
    D, P = inp["D"], inp["PoD"] * inp["D"]
    if inp["lattice"] == "quad":
        return M.QuadGeom(P, D, 0.1 * D)
    Dw = inp["dw_frac"] * (P - D)
    return M.TriGeom(P, D, Dw, inp["HoD"] * D, inp["Nr"], (inp["WoD"] - 1) * D)


def f3(g, inp, t, Re, Reb, info=None):
    if inp["lattice"] == "quad":
        return M.friction_quad(g, t, Re, Reb, info)
    return M.friction_tri(g, t, Re, g.flow_area(t), g.wetted_perimeter(t), Reb, info)


def of3(g, inp, t, Re, Reb):
    if inp["lattice"] == "quad":
        return o_friction_quad(g.pitch, g.D, g.side_gap, t, Re, Reb)
    return o_friction_tri(g.pitch, g.D, g.dwire, g.hwire, g.gap, t, Re, g.flow_area(t), g.wetted_perimeter(t), Reb)


def t3F_eval(inp):
    g = geom3(inp)
    t, Reb, Rei = inp["type"], inp["Reb"], inp["Rei"]
    info = {}
    f = f3(g, inp, t, Rei, Reb, info)
    ch = [o3("t3_I1_friction_vs_oracle", f, of3(g, inp, t, Rei, Reb))]
    ReL, ReT = info["ReL"], info["ReT"]
    ie = {}
    ch.append(o3("t3_I1_f(Reb=ReL)==fL", f3(g, inp, t, Rei, ReL, ie), ie["fL"]))
    ie = {}
    ch.append(o3("t3_I1_f(Reb=ReT)==fT", f3(g, inp, t, Rei, ReT, ie), ie["fT"]))
    i2 = {}
    f3(g, inp, t, Rei * 0.37 + 1.0, Reb, i2)
    ch.append(SL.cond("t3_I1_psi_independent_of_local_Re", i2["psi"] == info["psi"], observed=[info["psi"], i2["psi"]]))
    if ReL <= Reb <= ReT:
        ch.append(SL.cond("t3_I1_psi_in_[0,1]", 0.0 <= info["psi"] <= 1.0, observed=info["psi"]))
    f1 = f3(g, inp, t, 1.0, Reb)
    probes = [f3(g, inp, t, r, Reb) for r in (0.0, 1e-12, 0.1)]
    ch.append(SL.cond("t3_I6_friction_constant_for_0<=Re_i<=1", all(p == f1 for p in probes) and math.isfinite(f1),
                      observed=probes + [f1]))
    fp = f3(g, inp, t, 1.0 + 1e-6, Reb)
    ch.append(SL.cond("t3_I6_friction_finite_at_1+1e-6", math.isfinite(fp), observed=fp))
    if inp["lattice"] != "quad":
        A, _, Dhb = g.assembly()
        rng = np.random.default_rng(abs(hash((round(Reb, 6), round(Rei, 6)))) % (2 ** 32))
        counts = g.channel_counts()
        n = sum(counts.values())
        w = rng.uniform(0.2, 1.0, n)
        w /= w.sum()
        mu = rng.uniform(1e-4, 1e-3, n)
        mub = float(np.dot(w, mu))
        Mflow = Reb * mub * A / Dhb
        Rb = M.bulk_reynolds(list(w * Mflow), list(mu), A, Dhb)
        ch.append(o3("t3_I1_bulk_Re_vs_independent_sum", Rb, Mflow * Dhb / (mub * A)))
        w2 = rng.uniform(0.2, 1.0, n)
        w2 /= w2.sum()
        Rb2 = M.bulk_reynolds(list(w2 * Mflow), [mub] * n, A, Dhb)
        ch.append(o3("t3_I1_bulk_Re_same_M_mub_other_partition", Rb2, Rb))
    return ch


def t3M_sample(kind):
    def s(rng):
        return dict(branch="M", kind=kind, PoD=float(rng.uniform(1.067, 1.35)), HoD=float(rng.uniform(4, 52)),
                    D=float(rng.uniform(0.005, 0.015)), dw_frac=float(1.0 - rng.uniform(0, 1)),
                    WoD=float(rng.uniform(1.01, 1.4)), Nr=nr_sample(rng, (2, 3, 9)), Reb=float(rng.uniform(400, 1e6)),
                    pair=int(rng.integers(0, 2)))
    return s


PAIRS = {"Cm": [(CENTER, CENTER), (CENTER, EDGE)], "Cs": [(EDGE, EDGE), (EDGE, CORNER)]}


def t3M_eval(inp):
    D, P = inp["D"], inp["PoD"] * inp["D"]
    g = M.TriGeom(P, D, inp["dw_frac"] * (P - D), inp["HoD"] * D, inp["Nr"], (inp["WoD"] - 1) * D)
    M.mixing_constructor(g)
    ti, tj = PAIRS[inp["kind"]][inp["pair"]]
    Reb = inp["Reb"]
    loc = local_fields(g, ti, tj, Reb)
    fn, other = (M.mixing_beta, M.sweep_beta) if inp["kind"] == "Cm" else (M.sweep_beta, M.mixing_beta)
    info = {}
    beta = fn(g, ti, tj, Reb, loc, info)
    ob, C, CL, CT = o_mixing(P, D, g.dwire, g.hwire, g.gap, g.nrings, inp["kind"], Reb)
    ch = [o3("t3_I1_mixing_beta_vs_oracle", beta, ob),
          SL.cond("t3_I1_CL<=C<=CT", min(info["CL"], info["CT"]) - 1e-12 <= info["C"] <= max(info["CL"], info["CT"]) + 1e-12,
                  observed=[info["CL"], info["C"], info["CT"]])]
    for lab, R0, key in (("ReL", info["ReL"], "CL"), ("ReT", info["ReT"], "CT")):
        ie = {}
        fn(g, ti, tj, R0, loc, ie)
        ch.append(o3(f"t3_I1_C(Reb={lab})=={key}", ie["C"], ie[key]))
    i2 = {}
    fn(g, ti, tj, Reb * 1.001, loc, i2)
    ch.append(SL.cond("t3_I1_C_nondecreasing_in_Reb", i2["C"] >= info["C"] - 1e-12, observed=[info["C"], i2["C"]]))
    z = other(g, ti, tj, Reb, loc)
    ch.append(SL.cond("t3_I1_other_mixing_function_exact_zero", z == 0.0, observed=z, expected=0.0))
    loc2 = local_fields(g, ti, tj, Reb, scale_m=3.0, scale_mu=0.5)
    ch.append(SL.cond("t3_I1_beta_independent_of_local_fields", fn(g, ti, tj, Reb, loc2) == beta,
                      observed=[beta, fn(g, ti, tj, Reb, loc2)]))
    return ch


def t3B_sample(t):
    def s(rng):
        r = float(rng.uniform(1.1, 1.5))
        ReT = o_limits(r)[1]
        return dict(branch="B", type=t, r=r, D=float(rng.uniform(0.005, 0.015)),
                    Rei=float(rng.uniform(ReT * (1 + 1e-9), 1e6)), Pe=float(rng.uniform(0, 2200)))
    return s


def t3B_eval(inp):
    t, r, D, Re = inp["type"], inp["r"], inp["D"], inp["Rei"]
    Pr = inp["Pe"] / Re
    P = r * D
    g = M.TriGeom(P, D, 0.0, 0.0, 3, T_GAP / T_D * D, check_mesh=False)
    S, wp = g.flow_area(t), g.wetted_perimeter(t)
    bi = {}
    Nu = M.borishanskii_Nu(M.nusselt_preinfo(P, D, t, Re, Pr), bi)
    ch = [o3("t3_I2_Borishanskii_vs_oracle", Nu, o_borishanskii(r, t, Re, Pr)),
          SL.cond("t3_I2_Nu>0", Nu > 0, observed=Nu)]
    if inp["Pe"] <= 200:
        ch.append(o3("t3_I2_Pe<=200_plateau_B(r)", Nu, borishanskii_formula(r)))
    else:
        Nu2 = M.borishanskii_Nu(M.nusselt_preinfo(P, D, t, Re, Pr * 1.001))
        ch.append(SL.cond("t3_I2_Nu_increasing_in_Pe_above_200", Nu2 >= Nu, observed=[Nu, Nu2]))
    for k in (0.1, 1.0, 10.0, 100.0):
        h = M.compute_htc(Nu, k, S, wp)
        ch.append(o3(f"t3_I2_hDh/k==Nu_k={k:g}", h * (4 * S / wp) / k, Nu))
    for lam in (0.1, 10.0):
        gs = M.TriGeom(P * lam, D * lam, 0.0, 0.0, 3, T_GAP / T_D * D * lam, check_mesh=False)
        Nus = M.borishanskii_Nu(M.nusselt_preinfo(P * lam, D * lam, t, Re, Pr))
        ch.append(o3(f"t3_I2_length_scale_{lam:g}_Nu_unchanged", Nus, Nu))
        hs = M.compute_htc(Nus, 1.0, gs.flow_area(t), gs.wetted_perimeter(t))
        ch.append(o3(f"t3_I2_length_scale_{lam:g}_h_scales_inverse", hs, M.compute_htc(Nu, 1.0, S, wp) / lam))
    Dp = 0.98 * D
    Nud = M.borishanskii_Nu(M.nusselt_preinfo(P, D, t, Re, Pr, Dpin=Dp))
    ch.append(o3("t3_I2_deformed_Dpin_vs_oracle", Nud, o_borishanskii(P / Dp, t, Re, Pr)))
    ReL, ReT = o_limits(r)
    for lab, Rp in (("ReL", ReL), ("ReL(1+1e-6)", ReL * (1 + 1e-6)), ("ReT(1-1e-6)", ReT * (1 - 1e-6)), ("ReT", ReT)):
        Nb = M.borishanskii_Nu(M.nusselt_preinfo(P, D, t, Rp, Pr))
        ch.append(o3(f"t3_I2_blend_probe_{lab}_vs_oracle", Nb, o_borishanskii(r, t, Rp, Pr)))
        ch.append(SL.cond(f"t3_I2_blend_probe_{lab}_positive", Nb > 0, observed=Nb))
    return ch


def t3G_sample(lat, t):
    def s(rng):
        d = dict(branch="G", lattice=lat, type=t, PoD=float(rng.uniform(1.01, 1.4)), D=float(rng.uniform(0.005, 0.015)),
                 Reb=float(rng.uniform(50, 1e6)), Rei=float(rng.uniform(50, 1e6)), Pr=float(rng.uniform(1e-5, 2000)))
        if lat == "tri_wire":
            d.update(WoD=float(rng.uniform(1.01, 1.4)), Nr=nr_sample(rng), HoD=float(rng.uniform(8, 52)),
                     dw_frac=float(1.0 - rng.uniform(0, 1)))
        return d
    return s


def t3G_eval(inp):
    g = geom3(inp)
    t, Reb, Rei, Pr = inp["type"], inp["Reb"], inp["Rei"], inp["Pr"]
    f = f3(g, inp, t, Rei, Reb)
    gi = {}
    pre = M.nusselt_preinfo(g.pitch, g.D, t, Rei, Pr)
    Nu = M.gnielinski_Nu(pre, f, gi)
    pod = g.pitch / g.D
    ReL, ReT = o_limits(pod)
    ch = [o3("t3_I3_Gnielinski_vs_oracle", Nu, o_gnielinski(pod, t, Rei, Pr, f)),
          SL.cond("t3_I3_Nu_finite", math.isfinite(Nu), observed=Nu),
          SL.cond("t3_I3_Nu>=0_passive_(incl_blend)", Nu >= 0, observed=Nu, expected=">=0")]
    if Rei >= max(ReT, 1000.0):
        bound = 10 * (f / 8) * (Rei - 1000) * (Pr + 0.01)
        ch.append(SL.cond("t3_I3_Nu<=10(f/8)(Re-1000)(Pr+0.01)", Nu <= bound * (1 + 1e-12), observed=Nu, expected=bound))
    if gi["denominator"] >= 0.1:
        unreg = (f / 8) * (Rei - 1000) * (Pr + 0.01) / gi["denominator"]
        ch.append(o3("t3_I3_unfloored_equals_unregularized_formula", gi["NuT"], unreg))
    else:
        ch.append(SL.cond("t3_I3_floor_warning_flagged", "Gnielinski denominator floored" in gi["warnings"],
                          observed=gi["warnings"]))
    S = g.flow_area(t)
    wp = g.wetted_perimeter(t)
    for k in (0.1, 1.0, 10.0, 100.0):
        if Nu >= 0:
            h = M.compute_htc(Nu, k, S, wp)
            ch.append(o3(f"t3_I3_hDh/k==Nu_k={k:g}", h * (4 * S / wp) / k, Nu))
        else:
            try:
                M.compute_htc(Nu, k, S, wp)
                raised = False
            except MooseError:
                raised = True
            ch.append(SL.cond(f"t3_I3_negative_HTC_rejected_k={k:g}", raised, observed=raised, expected=True))
    Nud = M.gnielinski_Nu(M.nusselt_preinfo(g.pitch, g.D, t, Rei, Pr, is_duct=True), f)
    ch.append(o3("t3_I3_duct_sentinel_same_as_pin", Nud, Nu))
    return ch


# ---------------- Track 3 mandatory points
def t3_grid_wire_eval(inp):
    g = M.TriGeom(1.2 * T_D, T_D, inp["dw_frac"] * 0.2 * T_D, inp["HoD"] * T_D, 3, T_GAP)
    t = inp["type"]
    info = {}
    f3(g, dict(lattice="tri_wire"), t, 1e4, 1e4, info)
    ReL, ReT = info["ReL"], info["ReT"]
    Rb = {"0.5ReL": 0.5 * ReL, "ReL": ReL, "sqrt": math.sqrt(ReL * ReT), "ReT": ReT, "2ReT": 2 * ReT}[inp["Rb_kind"]]
    f = f3(g, dict(lattice="tri_wire"), t, inp["Rei"], Rb)
    return [o3("t3_I1_grid_friction_vs_oracle", f, of3(g, dict(lattice="tri_wire"), t, inp["Rei"], Rb)),
            SL.cond("t3_I1_grid_friction_finite", math.isfinite(f), observed=f)]


def t3_switch_eval(inp):
    lat, t, x = inp["lattice"], inp["type"], inp["x"]
    if lat == "quad":
        g = M.QuadGeom(x * T_D, T_D, 0.1 * T_D)
    else:
        dw = 0.5 * (x - 1) * T_D if lat == "tri_wire" else 0.0
        g = M.TriGeom(x * T_D, T_D, dw, 20 * T_D if lat == "tri_wire" else 0.0, 3, (x - 1) * T_D)
    f = f3(g, dict(lattice=lat), t, 1e4, 5e3)
    return [o3("t3_I1_table_switch_vs_oracle", f, of3(g, dict(lattice=lat), t, 1e4, 5e3))]


def t3_mix_grid_eval(inp):
    g = M.TriGeom(1.2 * T_D, T_D, 0.2 * T_D, 20 * T_D, inp["Nr"], T_GAP, check_mesh=inp["Nr"] >= 2)
    ReL, ReT = o_limits(1.2)
    Rb = {"0.5ReL": 0.5 * ReL, "ReL": ReL, "sqrt": math.sqrt(ReL * ReT), "ReT": ReT, "2ReT": 2 * ReT}[inp["Rb_kind"]]
    ti, tj = inp["pair"]
    loc = local_fields(g, ti, tj, Rb)
    fn = M.mixing_beta if inp["kind"] == "Cm" else M.sweep_beta
    beta = fn(g, ti, tj, Rb, loc)
    ob = o_mixing(g.pitch, g.D, g.dwire, g.hwire, g.gap, g.nrings, inp["kind"], Rb)[0]
    return [o3("t3_I1_mixing_grid_vs_oracle", beta, ob)]


def t3_boris_grid_eval(inp):
    t, r, Pe = inp["type"], inp["r"], inp["Pe"]
    Re = 2 * o_limits(r)[1]
    bi = {}
    Nu = M.borishanskii_Nu(M.nusselt_preinfo(r * T_D, T_D, t, Re, Pe / Re), bi)
    ch = [o3("t3_I2_Borishanskii_grid_vs_oracle", Nu, o_borishanskii(r, t, Re, Pe / Re))]
    if Pe > 2200:
        ch.append(SL.cond("t3_I2_Pe>2200_warning_flagged", any("Pe above" in w for w in bi["warnings"]),
                          observed=bi["warnings"]))
    return ch


def t3_boris_duct_eval(inp):
    try:
        M.borishanskii_Nu(M.nusselt_preinfo(1.2 * T_D, T_D, CENTER, 1e5, 0.005, is_duct=True), is_duct_closure=True)
        raised = False
    except MooseError:
        raised = True
    return [SL.cond("t3_I2_Borishanskii_as_duct_closure_errors", raised, observed=raised, expected=True)]


def t3_gn_grid_eval(inp):
    lat, t, Pr = inp["lattice"], inp["type"], inp["Pr"]
    g = M.QuadGeom(1.2 * T_D, T_D, 0.1 * T_D) if lat == "quad" else M.TriGeom(1.2 * T_D, T_D, 0.2 * T_D, 20 * T_D, 3, T_GAP)
    ReL, ReT = o_limits(1.2)
    Re = {"lam": 100.0, "trans": math.sqrt(ReL * ReT), "turb": 2e5}[inp["Re_kind"]]
    f = f3(g, dict(lattice=lat), t, Re, Re)
    Nu = M.gnielinski_Nu(M.nusselt_preinfo(g.pitch, g.D, t, Re, Pr), f)
    return [o3("t3_I3_Gnielinski_grid_vs_oracle", Nu, o_gnielinski(1.2, t, Re, Pr, f)),
            SL.cond("t3_I3_grid_Nu>=0", Nu >= 0, observed=Nu)]


def t3_gn_floor_eval(inp):
    g = M.TriGeom(1.2 * T_D, T_D, 0.2 * T_D, 20 * T_D, 3, T_GAP)
    t, Pr = CENTER, 0.005
    p = Pr + 0.01
    fc = 8 * (0.9 / (12.7 * (1 - p ** (2 / 3)))) ** 2
    Reb = 2 * o_limits(1.2)[1]
    info = {}
    f3(g, dict(lattice="tri_wire"), t, 1e4, Reb, info)
    ftarget = fc * (1 + inp["s"])
    Re = (info["CfT"] / ftarget) ** (1 / 0.18)
    f = f3(g, dict(lattice="tri_wire"), t, Re, Reb)
    gi = {}
    Nu = M.gnielinski_Nu(M.nusselt_preinfo(g.pitch, g.D, t, Re, Pr), f, gi)
    ch = [o3("t3_I3_floor_probe_vs_oracle", Nu, o_gnielinski(1.2, t, Re, Pr, f)),
          SL.cond("t3_I3_floor_probe_friction_straddles_fc", (f > fc) == (inp["s"] > 0) or inp["s"] == 0,
                  observed=[f, fc])]
    if gi["floored"]:
        ch.append(SL.cond("t3_I3_floor_probe_warning", "Gnielinski denominator floored" in gi["warnings"],
                          observed=gi["warnings"]))
    ch.append(SL.cond("t3_I3_floor_probe_Nu>=0", Nu >= 0, observed=Nu))
    return ch


def t3_geom_eval(inp):
    D = inp["D"]
    P = inp["PoD"] * D
    gapv = P - D
    dw = {"0": 0.0, "P-D": gapv, "P-D+0.5eps": gapv + 0.5e-12, "P-D+2eps": gapv + 2e-12}[inp["dw_kind"]]
    try:
        M.tri_mesh_generator_checks(P, D, dw, 3)
        outcome = "accept"
    except MooseError as e:
        outcome = "reject:" + ("pitch" if "pitch" in str(e) else "dwire")
    if P <= D:
        expected = "reject:pitch"
    elif dw > gapv + 1e-12:
        expected = "reject:dwire"
    else:
        expected = "accept"
    return [SL.cond("t3_I5_geometry_admissibility", outcome == expected, observed=outcome, expected=expected)]


def t3_bulk_scale_eval(inp):
    g = M.TriGeom(1.2 * T_D, T_D, 0.2 * T_D, 20 * T_D, 3, T_GAP)
    A, _, Dhb = g.assembly()
    n = sum(g.channel_counts().values())
    rng = np.random.default_rng(9)
    w = rng.uniform(0.2, 1.0, n)
    w /= w.sum()
    mu = rng.uniform(1e-4, 1e-3, n)
    M1 = 5e4 * float(np.dot(w, mu)) * A / Dhb
    R1 = M.bulk_reynolds(list(w * M1), list(mu), A, Dhb)
    tt = inp["t"]
    Rt = M.bulk_reynolds(list(w * M1 * tt), list(mu), A, Dhb)
    if tt == 0:
        return [SL.cond("t3_I6_zero_flow_bulk_Re_finite_or_diagnostic", math.isfinite(Rt), observed=Rt,
                        expected="finite value or explicit error")]
    return [o3(f"t3_I6_bulk_Re_scales_linearly_t={tt:g}", Rt, tt * R1)]


def t3_mixing_nr1_eval(inp):
    g = M.TriGeom(1.2 * T_D, T_D, 0.2 * T_D, 20 * T_D, 1, T_GAP, check_mesh=False)
    out = []
    for kind, pair in (("Cm", (CENTER, CENTER)), ("Cs", (EDGE, EDGE))):
        fn = M.mixing_beta if kind == "Cm" else M.sweep_beta
        b = fn(g, pair[0], pair[1], 5e4, local_fields(g, pair[0], pair[1], 5e4))
        out.append(o3(f"t3_I1_one_ring_probe_{kind}_vs_oracle", b, o_mixing(g.pitch, g.D, g.dwire, g.hwire, g.gap, 1, kind, 5e4)[0]))
    return out


# ================================================================ Assembly
def main():
    a = SL.parse_args()
    notes = dict(native=False, not_covered=NOT_COVERED[a.track])
    run = SL.Runner(PR, PR_IDX, a.track, a.n, a.out, versions(), notes)
    mand = []
    if a.track == 1:
        branches = [(f"tri_{t}", t1_sample(t), t1_eval) for t in TYPES]
        for t in TYPES:
            for hod in (8.0, 20.0, 52.0):
                mand.append((f"{t}_HoD={hod:g}_Reb=5e4", dict(type=t, PoD=1.2, HoD=hod, DwoD=0.2, Npin=19, Reb=5e4,
                                                               Pr=0.005, Pe=1000.0), t1_eval))
    elif a.track == 2:
        branches = [(f"tri_{t}", t2_sample(t), t2_eval) for t in TYPES]
        for t in TYPES:
            for hod in (8.0, 15.0, 25.0, 52.0):
                mand.append((f"{t}_HoD={hod:g}", dict(type=t, PoD=1.2, HoD=hod, Npin=19, Reb=5e4, Pr=0.005, PoD_B=1.2,
                                                      Pe=200.0), t2_eval))
            for pod, ap, at in ((1.1, 4.25, 0.005), (1.5, 10.8, 0.05)):
                for pr in (0.004, 0.03):
                    mand.append((f"{t}_Boris_P/D={pod}_Pr={pr}", dict(type=t, PoD_B=pod, Pr=pr, approx=ap, atol=at),
                                 t2_boris_value_eval))
    else:
        branches = []
        for lat in ("tri_bare", "tri_wire", "quad"):
            for t in TYPES:
                branches.append((f"F_{lat}_{t}", t3F_sample(lat, t), t3F_eval))
        for kind in ("Cm", "Cs"):
            branches.append((f"M_{kind}", t3M_sample(kind), t3M_eval))
        for t in TYPES:
            branches.append((f"B_{t}", t3B_sample(t), t3B_eval))
        for lat in ("tri_wire", "quad"):
            for t in TYPES:
                branches.append((f"G_{lat}_{t}", t3G_sample(lat, t), t3G_eval))
        for t in TYPES:
            for hod in (8.0, 20.0, 52.0):
                for fr in (0.1, 0.5, 1.0):
                    for rk in ("0.5ReL", "ReL", "sqrt", "ReT", "2ReT"):
                        for ri in (0.0, 0.1, 1.0, 50.0, 1e3, 1e5):
                            mand.append((f"grid_{t}_H{hod:g}_w{fr}_{rk}_Ri{ri:g}",
                                         dict(type=t, HoD=hod, dw_frac=fr, Rb_kind=rk, Rei=ri), t3_grid_wire_eval))
        for lat in ("tri_bare", "tri_wire", "quad"):
            for t in TYPES:
                for x in (1.1 - 1e-6, 1.1, 1.1 + 1e-6):
                    mand.append((f"switch_{lat}_{t}_{x!r}", dict(lattice=lat, type=t, x=x), t3_switch_eval))
        for nr in (2, 3, 9):
            for kind in ("Cm", "Cs"):
                for pair in PAIRS[kind]:
                    for rk in ("0.5ReL", "ReL", "sqrt", "ReT", "2ReT"):
                        mand.append((f"mix_{nr}_{kind}_{pair}_{rk}", dict(Nr=nr, kind=kind, pair=list(pair), Rb_kind=rk),
                                     t3_mix_grid_eval))
        mand.append(("mixing_one_ring_probe", dict(Nr=1), t3_mixing_nr1_eval))
        for t in TYPES:
            for r in (1.1, 1.2, 1.35, 1.5):
                for pe in (1.0, 199.999, 200.0, 200.001, 500.0, 2200.0, 2500.0):
                    mand.append((f"boris_{t}_{r}_{pe}", dict(type=t, r=r, Pe=pe), t3_boris_grid_eval))
        mand.append(("boris_as_duct_closure", dict(case="duct"), t3_boris_duct_eval))
        for lat in ("tri_wire", "quad"):
            for t in TYPES:
                for pr in (1e-5, 0.001, 0.005, 0.01, 0.5, 1.0, 10.0, 2000.0):
                    for rk in ("lam", "trans", "turb"):
                        mand.append((f"gn_{lat}_{t}_{pr}_{rk}", dict(lattice=lat, type=t, Pr=pr, Re_kind=rk), t3_gn_grid_eval))
        for s in (-1e-6, 0.0, 1e-6):
            mand.append((f"gn_floor_s={s}", dict(s=s), t3_gn_floor_eval))
        for D in (T_D, 10 * T_D):
            for pod in (1.0, 1.0 - 1e-6, 1.0 + 1e-6, 1.1):
                for dk in ("0", "P-D", "P-D+0.5eps", "P-D+2eps"):
                    mand.append((f"geom_D{D}_{pod!r}_{dk}", dict(D=D, PoD=pod, dw_kind=dk), t3_geom_eval))
        for tt in (1.0, 1e-3, 1e-6, 0.0):
            mand.append((f"bulk_scale_t={tt:g}", dict(t=tt), t3_bulk_scale_eval))
    run.run(branches, mand, progress_every=2000)


if __name__ == "__main__":
    main()
