"""Line-by-line Python recomputation of the subchannel closure relations of MOOSE#33307 merge commit b5a1d9e (not native: there is no MOOSE build on the node).
C++ functions covered (behavior corresponds line by line to the source; mooseError -> raise MooseError, flagSolutionWarning -> recorded in the warnings list):
- SCMFrictionUpgradedChengTodreas: checks and warnings in the constructor, computeTriLatticeFrictionFactor, computeQuadLatticeFrictionFactor
- SCMMixingChengTodreas: constructor warnings, computeMixingParameter, computeSweepFlowMixingParameter (geometric quantities are passed in directly; the solution field is not read)
- SCMHTCClosureBase: computeNusseltNumberPreInfo, computeHTC, turbulentReynoldsNumber, blendTurbulentNusseltNumber
- SCMHTCBorishanskii::computeNusseltNumber, SCMHTCGnielinski::computeNusseltNumber
- SubChannel1PhaseProblem::computeBulkReynoldsNumber (single-process semantics)
- the nrings / pitch / dwire checks in the SCMTriAssemblyMeshGenerator constructor (MooseUtils::absoluteFuzzyGreaterThan, default tolerance 1e-12)
- TriSubChannelMesh::getSubchannelFlowArea / getSubchannelWettedPerimeter / computeAssemblyHydraulicParameters (unchanged; used by the functions above to obtain geometry)
- QuadSubChannelMesh::getSubchannelFlowArea / getSubchannelWettedPerimeter (unchanged)
Not covered: the cross-flow estimate indexing of implicitPetscSolve (requires a PETSc solver context), the full externalSolve flow, MPI behavior.
Internal quantities (psi, ReL, ReT, WsT, correction factors, Cf, etc.) are returned through the info dictionary for use by the checks; the return values themselves are consistent with C++.
"""
import math

PI = math.pi
TOL_FUZZY = 1e-6 * 1e-6  # libMesh::TOLERANCE * libMesh::TOLERANCE


class MooseError(Exception):
    pass


CENTER, EDGE, CORNER = "CENTER", "EDGE", "CORNER"


class TriGeom:
    def __init__(self, pitch, D, dwire, hwire, nrings, duct_gap, check_mesh=True):
        self.pitch, self.D, self.dwire, self.hwire = float(pitch), float(D), float(dwire), float(hwire)
        self.nrings = int(nrings)
        self.gap = float(duct_gap)
        if check_mesh:
            tri_mesh_generator_checks(self.pitch, self.D, self.dwire, self.nrings)

    def theta(self):
        return math.acos(self.hwire / math.sqrt(self.hwire ** 2 + (PI * (self.D + self.dwire)) ** 2))

    def flow_area(self, t):
        th = self.theta()
        if t == CENTER:
            std = self.pitch ** 2 * math.sqrt(3.0) / 4.0
            rod = PI * self.D ** 2 / 8.0
            wire = PI * self.dwire ** 2 / 8.0 / math.cos(th)
        elif t == EDGE:
            std = self.pitch * (self.D / 2.0 + self.gap)
            rod = PI * self.D ** 2 / 8.0
            wire = PI * self.dwire ** 2 / 8.0 / math.cos(th)
        else:
            std = 1.0 / math.sqrt(3.0) * (self.D / 2.0 + self.gap) ** 2
            rod = PI * self.D ** 2 / 24.0
            wire = PI * self.dwire ** 2 / 24.0 / math.cos(th)
        return std - rod - wire

    def wetted_perimeter(self, t):
        th = self.theta()
        rc = PI * self.D
        wc = PI * self.dwire
        if t == CENTER:
            return 0.5 * rc + 0.5 * wc / math.cos(th)
        if t == EDGE:
            return 0.5 * rc + 0.5 * wc / math.cos(th) + self.pitch
        return (rc + wc / math.cos(th)) / 6.0 + 2.0 / math.sqrt(3.0) * (self.D / 2.0 + self.gap)

    def channel_counts(self):
        nr = self.nrings
        return {CENTER: 6 * (nr - 1) ** 2, EDGE: 6 * (nr - 1), CORNER: 6}

    def assembly(self):
        A = 0.0
        Pw = 0.0
        for t, n in self.channel_counts().items():
            A += n * self.flow_area(t)
            Pw += n * self.wetted_perimeter(t)
        if Pw == 0.0:
            raise MooseError("Assembly wetted perimeter is zero; cannot compute hydraulic diameter.")
        return A, Pw, 4.0 * A / Pw

    def Dh(self, t):
        return 4.0 * self.flow_area(t) / self.wetted_perimeter(t)


class QuadGeom:
    def __init__(self, pitch, D, side_gap):
        self.pitch, self.D, self.side_gap = float(pitch), float(D), float(side_gap)

    def flow_area(self, t):
        if t == CORNER:
            std = 0.25 * self.pitch ** 2
            rod = 0.25 * 0.25 * PI * self.D ** 2
            add = self.pitch * self.side_gap + self.side_gap ** 2
        elif t == EDGE:
            std = 0.5 * self.pitch ** 2
            rod = 0.5 * 0.25 * PI * self.D ** 2
            add = self.pitch * self.side_gap
        else:
            std = self.pitch ** 2
            rod = 0.25 * PI * self.D ** 2
            add = 0.0
        return std + add - rod

    def wetted_perimeter(self, t):
        rc = PI * self.D
        if t == CORNER:
            return 0.25 * rc + self.pitch + 2.0 * self.side_gap
        if t == EDGE:
            return 0.5 * rc + self.pitch
        return rc

    def Dh(self, t):
        return 4.0 * self.flow_area(t) / self.wetted_perimeter(t)


def tri_mesh_generator_checks(pitch, D, dwire, nrings):
    if nrings < 2:
        raise MooseError("'nrings' must be at least 2.")
    pin_gap = pitch - D
    if pin_gap <= 0.0:
        raise MooseError("The pin pitch must be greater than the pin diameter")
    if dwire > pin_gap + TOL_FUZZY:  # absoluteFuzzyGreaterThan(_dwire, pin_gap)
        raise MooseError("The wire diameter must not exceed the pin-to-pin gap")


def friction_constructor_tri(g, Reb):
    """SCMFrictionUpgradedChengTodreas constructor (triangular-lattice branch): errors and warnings."""
    warnings = []
    if (g.dwire != 0.0) != (g.hwire != 0.0):
        raise MooseError("Wire-wrapped bundle friction requires both wire diameter and wire lead length.")
    has_wire = g.dwire != 0.0 and g.hwire != 0.0
    pod = g.pitch / g.D
    hd = g.hwire / g.D
    npins = 1 + 3 * g.nrings * (g.nrings - 1)
    if pod < 1.0 or pod > 1.42:
        warnings.append("P/D outside friction data range")
    if has_wire and (hd < 8.0 or hd > 52.0):
        warnings.append("H/D outside friction data range")
    if npins < 7 or npins > 271:
        warnings.append("Npin outside friction data range")
    if Reb < 50.0 or Reb > 1.0e6:
        warnings.append("Reb outside friction data range")
    return warnings


def friction_tri(g, ch_type, Re, S, w_perim, Reb, info=None):
    """Line-by-line recomputation of computeTriLatticeFrictionFactor."""
    if info is None:
        info = {}
    Re_eff = max(Re, 1.0)
    Dh_i = 4.0 * S / w_perim
    pitch, pin_diameter = g.pitch, g.D
    wire_lead_length, wire_diameter = g.hwire, g.dwire
    gap = g.gap
    has_wire = g.dwire != 0.0 and g.hwire != 0.0
    p_over_d = pitch / pin_diameter
    w_over_d = (pin_diameter + gap) / pin_diameter
    CbL1, CbL2, CbT1, CbT2 = 320.0, 1.0, 10000.0, 0.7
    ReL = math.pow(10, CbL2 * (p_over_d - 1)) * CbL1
    ReT = math.pow(10, CbT2 * (p_over_d - 1)) * CbT1
    psi = math.log(Reb / ReL) / math.log(ReT / ReL)
    if ch_type == CENTER:
        if p_over_d < 1.1:
            k1L, k2L, k3L, k1T, k2T, k3T = 26.0, 888.2, -3334.0, 0.09378, 1.398, -8.664
        else:
            k1L, k2L, k3L, k1T, k2T, k3T = 62.97, 216.9, -190.2, 0.1458, 0.03632, -0.03333
        x = p_over_d - 1
    elif ch_type == EDGE:
        if w_over_d < 1.1:
            k1L, k2L, k3L, k1T, k2T, k3T = 26.18, 554.5, -1480.0, 0.09377, 0.8732, -3.341
        else:
            k1L, k2L, k3L, k1T, k2T, k3T = 44.4, 256.7, -267.6, 0.1430, 0.04199, -0.04428
        x = w_over_d - 1
    else:
        if w_over_d < 1.1:
            k1L, k2L, k3L, k1T, k2T, k3T = 26.98, 1636.0, -10050.0, 0.1004, 1.625, -11.85
        else:
            k1L, k2L, k3L, k1T, k2T, k3T = 87.26, 38.59, -55.12, 0.1499, 0.006706, -0.009567
        x = w_over_d - 1
    CfL = k1L + k2L * x + k3L * x ** 2
    CfT = k1T + k2T * x + k3T * x ** 2
    info.update(CfL_bare=CfL, CfT_bare=CfT, ReL=ReL, ReT=ReT, psi=psi, has_wire=has_wire)
    if has_wire:
        theta = math.acos(wire_lead_length / math.sqrt(wire_lead_length ** 2 + (PI * (pin_diameter + wire_diameter)) ** 2))
        CwT1, CwT2, CwT3, CwT4 = 19.56, -98.71, 303.47, -0.541
        CwL1, CwL2 = 1.4, 1.0
        WdT = (CwT1 + CwT2 * (wire_diameter / pin_diameter) + CwT3 * (wire_diameter / pin_diameter) ** 2) * \
            math.pow(wire_lead_length / pin_diameter, CwT4)
        WdL = CwL1 * WdT
        a, b = -11.0, 19.0
        WsT = a * math.log10(wire_lead_length / pin_diameter) + b
        WsL = CwL2 * WsT
        info.update(WsT=WsT, WsL=WsL, WdT=WdT, theta=theta)
        if ch_type == CENTER:
            pw_p = PI * pin_diameter / 2.0
            ar = PI * (pin_diameter + wire_diameter) * wire_diameter / 6.0
            a_p = S + PI * wire_diameter ** 2 / 8.0 / math.cos(theta)
            dragT = WdT * (3.0 * ar / a_p) * (Dh_i / wire_lead_length) * math.pow(Dh_i / wire_diameter, 0.18)
            dragL = WdL * (3.0 * ar / a_p) * (Dh_i / wire_lead_length) * (Dh_i / wire_diameter)
            info.update(pw_ratio=pw_p / w_perim, dragT=dragT, dragL=dragL)
            CfT *= (pw_p / w_perim)
            CfT += dragT
            CfL *= (pw_p / w_perim)
            CfL += dragL
        else:
            if ch_type == EDGE:
                ar = PI * (pin_diameter + wire_diameter) * wire_diameter / 4.0
                a_p = S + PI * wire_diameter ** 2 / 8.0 / math.cos(theta)
            else:
                ar = PI * (pin_diameter + wire_diameter) * wire_diameter / 6.0
                a_p = S + PI * wire_diameter ** 2 / 24.0 / math.cos(theta)
            twc = 1 + WsT * (ar / a_p) * math.tan(theta) ** 2
            lwc = 1 + WsL * (ar / a_p) * math.tan(theta) ** 2
            info.update(turbulent_wire_correction=twc, laminar_wire_correction=lwc)
            if (not math.isfinite(twc)) or twc < 0.0:
                info["moose_error"] = "turbulent wire correction must be non-negative and finite"
                raise MooseError("The exponentiated term in the Cheng-Todreas turbulent wire correction must be "
                                 f"non-negative and finite. Computed {twc}.")
            CfT *= math.pow(twc, 1.41)
            CfL *= lwc
    mL, mT = -1.0, -0.18
    fL = CfL * math.pow(Re_eff, mL)
    fT = CfT * math.pow(Re_eff, mT)
    gamma, lam = 1.0 / 3.0, 7.0
    info.update(CfL=CfL, CfT=CfT, fL=fL, fT=fT)
    if Reb < ReL:
        info["regime"] = "laminar"
        return fL
    if Reb > ReT:
        info["regime"] = "turbulent"
        return fT
    info["regime"] = "transition"
    return fL * math.pow((1 - psi), gamma) * (1 - math.pow(psi, lam)) + fT * math.pow(psi, gamma)


def friction_quad(q, ch_type, Re, Reb, info=None):
    """Line-by-line recomputation of computeQuadLatticeFrictionFactor."""
    if info is None:
        info = {}
    Re_eff = max(Re, 1.0)
    pitch, pin_diameter, side_gap = q.pitch, q.D, q.side_gap
    w = (pin_diameter / 2.0) + (pitch / 2.0) + side_gap
    p_over_d = pitch / pin_diameter
    w_over_d = w / pin_diameter
    ReL = math.pow(10, (p_over_d - 1)) * 320.0
    ReT = math.pow(10, 0.7 * (p_over_d - 1)) * 1.0E+4
    psi = math.log(Reb / ReL) / math.log(ReT / ReL)
    if ch_type == CENTER:
        if p_over_d < 1.1:
            k1L, k2L, k3L, k1T, k2T, k3T = 26.37, 374.2, -493.9, 0.09423, 0.5806, -1.239
        else:
            k1L, k2L, k3L, k1T, k2T, k3T = 35.55, 263.7, -190.2, 0.1339, 0.09059, -0.09926
        x = p_over_d - 1
    elif ch_type == EDGE:
        if p_over_d < 1.1:
            k1L, k2L, k3L, k1T, k2T, k3T = 26.18, 554.5, -1480, 0.09377, 0.8732, -3.341
        else:
            k1L, k2L, k3L, k1T, k2T, k3T = 44.40, 256.7, -267.6, 0.1430, 0.04199, -0.04428
        x = w_over_d - 1
    else:
        if p_over_d < 1.1:
            k1L, k2L, k3L, k1T, k2T, k3T = 28.62, 715.9, -2807, 0.09755, 1.127, -6.304
        else:
            k1L, k2L, k3L, k1T, k2T, k3T = 58.83, 160.7, -203.5, 0.1452, 0.02681, -0.03411
        x = w_over_d - 1
    CfL = k1L + k2L * x + k3L * x ** 2
    CfT = k1T + k2T * x + k3T * x ** 2
    fL = CfL * math.pow(Re_eff, -1.0)
    fT = CfT * math.pow(Re_eff, -0.18)
    info.update(CfL=CfL, CfT=CfT, fL=fL, fT=fT, ReL=ReL, ReT=ReT, psi=psi)
    if Reb < ReL:
        info["regime"] = "laminar"
        return fL
    if Reb > ReT:
        info["regime"] = "turbulent"
        return fT
    info["regime"] = "transition"
    return fL * math.pow((1 - psi), 1.0 / 3.0) * (1 - math.pow(psi, 7.0)) + fT * math.pow(psi, 1.0 / 3.0)


def mixing_constructor(g):
    warnings = []
    if g.hwire == 0 or g.dwire == 0:
        raise MooseError("This correlation applies only for wire-wrapped assemblies")
    pod = g.pitch / g.D
    hd = g.hwire / g.D
    npins = 1 + 3 * g.nrings * (g.nrings - 1)
    if pod < 1.067 or pod > 1.35:
        warnings.append("P/D outside mixing range")
    if hd < 4.0 or hd > 52.0:
        warnings.append("H/D outside mixing range")
    if npins < 7 or npins > 217:
        warnings.append("Npin outside mixing range")
    return warnings


def _local_gap_Re(local):
    """Local Re in computeMixingParameter (used only for warnings). local=(Si_in,Sj_in,Si_out,Sj_out,wpi,wpj,mu*4,mdot*4)."""
    Si_in, Sj_in, Si_out, Sj_out = local["S"]
    mu_i_in, mu_j_in, mu_i_out, mu_j_out = local["mu"]
    md_i_in, md_j_in, md_i_out, md_j_out = local["mdot"]
    wpi, wpj = local["wp"]
    S_total = Si_in + Sj_in + Si_out + Sj_out
    Si = 0.5 * (Si_in + Si_out)
    Sj = 0.5 * (Sj_in + Sj_out)
    avg_mu = (1.0 / S_total) * (mu_i_out * Si_out + mu_i_in * Si_in + mu_j_out * Sj_out + mu_j_in * Sj_in)
    avg_hD = 4.0 * (Si + Sj) / (wpi + wpj)
    avg_massflux = 0.5 * ((md_i_in + md_j_in) / (Si_in + Sj_in) + (md_i_out + md_j_out) / (Si_out + Sj_out))
    return avg_massflux * avg_hD / avg_mu


def mixing_beta(g, type_i, type_j, bulk_Re, local, info=None):
    """Line-by-line recomputation of computeMixingParameter (geometry / solution-field quantities are passed in via local)."""
    if info is None:
        info = {}
    warnings = []
    pitch, D, H, Dw, Nr = g.pitch, g.D, g.hwire, g.dwire, g.nrings
    Re = _local_gap_Re(local)
    if Re < 400.0 or Re > 1.0e6:
        warnings.append("local Re outside mixing range")
    ReL = 320.0 * math.pow(10.0, pitch / D - 1.0)
    ReT = 10000.0 * math.pow(10.0, 0.7 * (pitch / D - 1.0))
    beta = 0.0
    info.update(ReL=ReL, ReT=ReT, warnings=warnings, local_Re=Re)
    if type_i == CENTER or type_j == CENTER:
        theta = math.acos(H / math.sqrt(H ** 2 + (PI * (D + Dw)) ** 2))
        Ar1 = PI * (D + Dw) * Dw / 6.0
        A1prime = (math.sqrt(3.0) / 4.0) * pitch ** 2 - PI * D ** 2 / 8.0
        if Nr == 1:
            CmT_c, CmL_c = 0.1, 0.055
        else:
            CmT_c, CmL_c = 0.14, 0.077
        CmT = CmT_c * math.pow((pitch - D) / D, -0.5)
        CmL = CmL_c * math.pow((pitch - D) / D, -0.5)
        if bulk_Re < ReL:
            Cm = CmL
        elif bulk_Re > ReT:
            Cm = CmT
        else:
            psi = math.log(bulk_Re / ReL) / math.log(ReT / ReL)
            Cm = CmL + (CmT - CmL) * math.pow(psi, 2.0 / 3.0)
            info["psi"] = psi
        beta = Cm * math.sqrt(Ar1 / A1prime) * math.tan(theta)
        info.update(C=Cm, CL=CmL, CT=CmT, geom=math.sqrt(Ar1 / A1prime) * math.tan(theta))
    return beta


def sweep_beta(g, type_i, type_j, bulk_Re, local, info=None):
    """Line-by-line recomputation of computeSweepFlowMixingParameter."""
    if info is None:
        info = {}
    warnings = []
    pitch, D, H, Dw, Nr = g.pitch, g.D, g.hwire, g.dwire, g.nrings
    Re = _local_gap_Re(local)
    if Re < 400.0 or Re > 1.0e6:
        warnings.append("local Re outside mixing range")
    ReL = 320.0 * math.pow(10.0, pitch / D - 1.0)
    ReT = 10000.0 * math.pow(10.0, 0.7 * (pitch / D - 1.0))
    beta = 0.0
    info.update(ReL=ReL, ReT=ReT, warnings=warnings, local_Re=Re)
    if type_i in (CORNER, EDGE) and type_j in (CORNER, EDGE):
        theta = math.acos(H / math.sqrt(H ** 2 + (PI * (D + Dw)) ** 2))
        dpgap = g.gap
        w = D + dpgap
        Ar2 = PI * (D + Dw) * Dw / 4.0
        A2prime = pitch * (w - D / 2.0) - PI * D ** 2 / 8.0
        if Nr == 1:
            CsT_c, CsL_c = 0.6, 0.33
        else:
            CsT_c, CsL_c = 0.75, 0.413
        CsL = CsL_c * math.pow(H / D, 0.3)
        CsT = CsT_c * math.pow(H / D, 0.3)
        if bulk_Re < ReL:
            Cs = CsL
        elif bulk_Re > ReT:
            Cs = CsT
        else:
            psi = math.log(bulk_Re / ReL) / math.log(ReT / ReL)
            Cs = CsL + (CsT - CsL) * math.pow(psi, 2.0 / 3.0)
            info["psi"] = psi
        beta = Cs * math.sqrt(Ar2 / A2prime) * math.tan(theta)
        info.update(C=Cs, CL=CsL, CT=CsT, geom=math.sqrt(Ar2 / A2prime) * math.tan(theta))
    return beta


def nusselt_preinfo(pitch, D_mesh, ch_type, Re, Pr, Dpin=None, is_duct=False):
    if not is_duct and Dpin is not None:
        if Dpin > 0:
            D = Dpin
        else:
            raise MooseError("The diameter of the pin is equal or smaller than zero")
    else:
        D = D_mesh
    poD = pitch / D
    lam = 3.73 if ch_type == CENTER else (3.59 if ch_type == EDGE else 3.52)
    return dict(Re=Re, Pr=Pr, poD=poD, subch_type=ch_type, laminar_Nu=lam,
                ReL=320 * math.pow(10.0, (poD - 1.0)), ReT=1e4 * math.pow(10.0, 0.7 * (poD - 1.0)))


def turbulent_Re(pre):
    if pre["Re"] > pre["ReL"] and pre["Re"] < pre["ReT"]:
        return pre["ReT"]
    return pre["Re"]


def blend(pre, NuT):
    if pre["Re"] <= pre["ReL"]:
        return pre["laminar_Nu"]
    if pre["Re"] >= pre["ReT"]:
        return NuT
    w = (pre["Re"] - pre["ReL"]) / (pre["ReT"] - pre["ReL"])
    return w * NuT + (1.0 - w) * pre["laminar_Nu"]


def borishanskii_Nu(pre, info=None, is_duct_closure=False):
    if info is None:
        info = {}
    if is_duct_closure:
        raise MooseError("'Borishanskii' is not yet supported for the 'duct_htc_correlation'.")
    warnings = []
    Pe = pre["Re"] * pre["Pr"]
    tPe = turbulent_Re(pre) * pre["Pr"]
    if pre["poD"] < 1.1 or pre["poD"] > 1.5:
        warnings.append("P/D out of range for Borishanskii")
    poly = -8.12 + 12.76 * pre["poD"] - 3.65 * pre["poD"] ** 2
    if poly <= 0.0:
        raise MooseError("Logarithm argument non-positive in Borishanskii correlation")
    NuT = 24.15 * math.log10(poly)
    corr = 0.0174 * (1.0 - math.exp(6.0 - 6.0 * pre["poD"]))
    if Pe > 2200.0:
        warnings.append("Pe above recommended range for Borishanskii")
    if 200.0 <= tPe <= 2200.0:
        NuT += corr * math.pow(tPe - 200.0, 0.9)
    elif tPe < 200.0:
        pass
    else:
        NuT += corr * math.pow(tPe - 200.0, 0.9)
    info.update(NuT=NuT, poly=poly, Pe=Pe, turbulent_Pe=tPe, warnings=warnings)
    return blend(pre, NuT)


def gnielinski_Nu(pre, f_friction, info=None):
    """f_friction: return value of computeFrictionFactor of the configured friction closure (divided by 8 again in C++)."""
    if info is None:
        info = {}
    warnings = []
    if pre["Pr"] < 1e-5 or pre["Pr"] > 2e3:
        warnings.append("Pr out of range for Gnielinski")
    f_turb = f_friction / 8.0
    Pr_s = pre["Pr"] + 0.01
    den = 1 + 12.7 * math.sqrt(f_turb) * (math.pow(Pr_s, 2. / 3.) - 1.)
    den_raw = den
    if den < 0.1:
        warnings.append("Gnielinski denominator floored")
        den = 0.1
    NuT = f_turb * (pre["Re"] - 1e3) * Pr_s / den
    info.update(NuT=NuT, denominator=den_raw, floored=den_raw < 0.1, f_turb=f_turb, warnings=warnings)
    return blend(pre, NuT)


def compute_htc(Nu, k, S, w_perim):
    Dh_i = 4.0 * S / w_perim
    htc = Nu * k / Dh_i
    if (not math.isfinite(htc)) or htc < 0.0:
        raise MooseError("The heat transfer coefficient must be non-negative and finite.")
    return htc


def bulk_reynolds(mdot_in, mu_in, assembly_area, assembly_Dh):
    """computeBulkReynoldsNumber (process 0; other processes return directly and keep the old value -- not simulated here)."""
    viscosity_in = 0.0
    mass_flow_in = 0.0
    for m, mu in zip(mdot_in, mu_in):
        viscosity_in += m * mu
        mass_flow_in += m
    try:
        inlet_mu = viscosity_in / mass_flow_in
    except ZeroDivisionError:
        inlet_mu = float("nan")  # C++ floating-point 0/0 = NaN, no exception thrown
    try:
        return mass_flow_in * assembly_Dh / (inlet_mu * assembly_area)
    except ZeroDivisionError:
        return float("nan")
