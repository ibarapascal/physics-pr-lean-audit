-- Moose33307: idaholab/moose PR #33307 "Use log10 in relevant SCM closures and update transient regime
-- calculation" (merge commit b5a1d9e3df14e98850e514f3e57d02c25c44fca9) Lean certificate.
-- Transcribes SCMFrictionUpgradedChengTodreas.C of modules/subchannel (triangular-lattice friction factor: bare-rod coefficients, wire-wrap correction,
-- laminar/transition/turbulent piecewise definition), the Cm/Cs transition interpolation of SCMMixingChengTodreas.C, the base Nu term of SCMHTCBorishanskii.C,
-- and the formula of SubChannel1PhaseProblem::computeBulkReynoldsNumber; proves that the piecewise definition is continuous at Re_L and Re_T (for all geometries),
-- non-negative in the transition region, that the wire-wrap sweeping-flow coefficient Ws_T is positive for H/D ≤ 52 (hence the wire-wrap correction ≥ 1 and the mooseError branch is unreachable),
-- and gives a counterexample to "the transition-region friction factor lies between the two ends" and a counterexample for the old ln version.
import Mathlib

/-!
# MOOSE #33307 — log10 in SCM closures and bulk-Re transition regime

Transcription conventions
* `std::pow(x, y)` is transcribed as `Real.rpow`; used only where the base is ≥ 0 (in the transition branch ψ ∈ [0,1] is proved).
  For a negative base C++ gives NaN while Lean's rpow gives some other real; the two are not compared.
* `std::log10` → `Real.logb 10`, `std::log` → `Real.log`.
* The `mooseError` branch is transcribed as `Option.none`; `!std::isfinite(...)` is always false under real semantics and is omitted.
* `_scm_problem.getBulkReynoldsNumber()` enters as the real input `Reb`; the issue that outside rank 0 it stays at the initial value 1.0 is at the level of the parallel
  implementation and is not formalized in this file (see the report).
* The piecewise structure of the quadrilateral lattice `computeQuadLatticeFrictionFactor` (:300-431) is identical to that of the triangular lattice (:320-323,
  :415-430); this file transcribes only the coefficient table of the triangular lattice; the piecewise continuity theorem `regimeFriction_continuousOn` applies equally to both
  (`ReL_quad_eq` / `ReT_quad_eq` prove that the Re_L and Re_T formulas at the two places are equal).
* Floating-point rounding is compared by the diff script `replay/scripts/Moose33307_replay.py`.
-/

namespace PhysAudit.Moose33307

noncomputable section

/-! ## §1 Piecewise structure (SCMFrictionUpgradedChengTodreas.C:104-111, :276-297) -/

/-- ```
SCMFrictionUpgradedChengTodreas.C:104-107
  CbL1 = 320;
  CbL2 = 1.0;
  CbT1 = 10000;
  CbT2 = 0.7;
``` -/
def CbL1 : ℝ := 320
def CbL2 : ℝ := 1.0
def CbT1 : ℝ := 10000
def CbT2 : ℝ := 0.7

/-- ```
SCMFrictionUpgradedChengTodreas.C:108
  const auto ReL = std::pow(10, CbL2 * (p_over_d - 1)) * CbL1;
``` -/
def ReL (pOverD : ℝ) : ℝ := (10 : ℝ) ^ (CbL2 * (pOverD - 1)) * CbL1

/-- ```
SCMFrictionUpgradedChengTodreas.C:109
  const auto ReT = std::pow(10, CbT2 * (p_over_d - 1)) * CbT1;
``` -/
def ReT (pOverD : ℝ) : ℝ := (10 : ℝ) ^ (CbT2 * (pOverD - 1)) * CbT1

/-- ```
SCMFrictionUpgradedChengTodreas.C:111
  const auto psi = std::log(Reb / ReL) / std::log(ReT / ReL);
SCMMixingChengTodreas.C:164 / :283
      const Real psi = std::log(bulk_Re / ReL) / std::log(ReT / ReL);
``` -/
def psi (Reb ReL ReT : ℝ) : ℝ := Real.log (Reb / ReL) / Real.log (ReT / ReL)

/-- Transition-region interpolation (`gamma = 1.0 / 3.0`, `lambda = 7.0`).
```
SCMFrictionUpgradedChengTodreas.C:280-281
  gamma = 1.0 / 3.0;
  lambda = 7.0;
SCMFrictionUpgradedChengTodreas.C:295-296
    return fL * std::pow((1 - psi), gamma) * (1 - std::pow(psi, lambda)) +
           fT * std::pow(psi, gamma);
``` -/
def transitionBlend (fL fT ψ : ℝ) : ℝ :=
  fL * (1 - ψ) ^ (1.0 / 3.0 : ℝ) * (1 - ψ ^ (7.0 : ℝ)) + fT * ψ ^ (1.0 / 3.0 : ℝ)

/-- ```
SCMFrictionUpgradedChengTodreas.C:282-297
  if (Reb < ReL)
  { return fL; }
  else if (Reb > ReT)
  { return fT; }
  else
  { return fL * std::pow((1 - psi), gamma) * (1 - std::pow(psi, lambda)) +
           fT * std::pow(psi, gamma); }
``` -/
def regimeFriction (fL fT Reb ReL ReT : ℝ) : ℝ :=
  if Reb < ReL then fL else if Reb > ReT then fT else transitionBlend fL fT (psi Reb ReL ReT)

/-- The original Re_L, Re_T of the quadrilateral lattice (:320-321), to show consistency with the triangular lattice.
```
SCMFrictionUpgradedChengTodreas.C:320-321
  const auto ReL = std::pow(10, (p_over_d - 1)) * 320.0;
  const auto ReT = std::pow(10, 0.7 * (p_over_d - 1)) * 1.0E+4;
``` -/
def ReLquad (pOverD : ℝ) : ℝ := (10 : ℝ) ^ (pOverD - 1) * 320.0
def ReTquad (pOverD : ℝ) : ℝ := (10 : ℝ) ^ (0.7 * (pOverD - 1)) * 1.0E+4

/-- Generic three-segment structure (for use in proofs; `regimeFriction` and `mixingRegime` are both instances of it, see the `rfl` lemmas). -/
def regime (a b : ℝ) (T : ℝ → ℝ) (Reb ReL ReT : ℝ) : ℝ :=
  if Reb < ReL then a else if Reb > ReT then b else T (psi Reb ReL ReT)

theorem regimeFriction_eq_regime (fL fT Reb ReL ReT : ℝ) :
    regimeFriction fL fT Reb ReL ReT = regime fL fT (transitionBlend fL fT) Reb ReL ReT :=
  rfl

/-! ## §2 Triangular-lattice coefficients (SCMFrictionUpgradedChengTodreas.C:77-279) -/

/-- `EChannelType`. -/
inductive EChannelType
  | CENTER
  | EDGE
  | CORNER
  deriving DecidableEq

/-- `k1 + k2 * x + k3 * Utility::pow<2>(x)` (:137-139 etc.). -/
def quadPoly (k1 k2 k3 x : ℝ) : ℝ := k1 + k2 * x + k3 * x ^ 2

/-- Bare-rod coefficients `(CfL, CfT)`.
```
SCMFrictionUpgradedChengTodreas.C:116-190 (summary)
  if (subch_type == EChannelType::CENTER)
  { if (p_over_d < 1.1) { k1L = 26.0; k2L = 888.2; k3L = -3334.0; k1T = 0.09378; k2T = 1.398; k3T = -8.664; }
    else { k1L = 62.97; k2L = 216.9; k3L = -190.2; k1T = 0.1458; k2T = 0.03632; k3T = -0.03333; }
    CfL = k1L + k2L * (p_over_d - 1) + k3L * Utility::pow<2>((p_over_d - 1));
    CfT = k1T + k2T * (p_over_d - 1) + k3T * Utility::pow<2>((p_over_d - 1)); }
  else if (subch_type == EChannelType::EDGE)
  { if (w_over_d < 1.1) { k1L = 26.18; k2L = 554.5; k3L = -1480.0; k1T = 0.09377; k2T = 0.8732; k3T = -3.341; }
    else { k1L = 44.4; k2L = 256.7; k3L = -267.6; k1T = 0.1430; k2T = 0.04199; k3T = -0.04428; }
    CfL = k1L + k2L * (w_over_d - 1) + k3L * Utility::pow<2>((w_over_d - 1)); (CfT likewise) }
  else
  { if (w_over_d < 1.1) { k1L = 26.98; k2L = 1636.0; k3L = -10050.0; k1T = 0.1004; k2T = 1.625; k3T = -11.85; }
    else { k1L = 87.26; k2L = 38.59; k3L = -55.12; k1T = 0.1499; k2T = 0.006706; k3T = -0.009567; }
    CfL = k1L + k2L * (w_over_d - 1) + k3L * Utility::pow<2>((w_over_d - 1)); (CfT likewise) }
``` -/
def triBareCoeffs (ty : EChannelType) (pOverD wOverD : ℝ) : ℝ × ℝ :=
  match ty with
  | .CENTER =>
    if pOverD < 1.1 then
      (quadPoly 26.0 888.2 (-3334.0) (pOverD - 1), quadPoly 0.09378 1.398 (-8.664) (pOverD - 1))
    else
      (quadPoly 62.97 216.9 (-190.2) (pOverD - 1), quadPoly 0.1458 0.03632 (-0.03333) (pOverD - 1))
  | .EDGE =>
    if wOverD < 1.1 then
      (quadPoly 26.18 554.5 (-1480.0) (wOverD - 1), quadPoly 0.09377 0.8732 (-3.341) (wOverD - 1))
    else
      (quadPoly 44.4 256.7 (-267.6) (wOverD - 1), quadPoly 0.1430 0.04199 (-0.04428) (wOverD - 1))
  | .CORNER =>
    if wOverD < 1.1 then
      (quadPoly 26.98 1636.0 (-10050.0) (wOverD - 1), quadPoly 0.1004 1.625 (-11.85) (wOverD - 1))
    else
      (quadPoly 87.26 38.59 (-55.12) (wOverD - 1), quadPoly 0.1499 0.006706 (-0.009567) (wOverD - 1))

/-- ```
SCMFrictionUpgradedChengTodreas.C:197-200
    const auto theta =
        std::acos(wire_lead_length /
                  std::sqrt(Utility::pow<2>(wire_lead_length) +
                            Utility::pow<2>(libMesh::pi * (pin_diameter + wire_diameter))));
``` -/
def theta (H D Dw : ℝ) : ℝ := Real.arccos (H / Real.sqrt (H ^ 2 + (Real.pi * (D + Dw)) ^ 2))

/-- ```
SCMFrictionUpgradedChengTodreas.C:202-211
    const Real CwT1 = 19.56;  const Real CwT2 = -98.71;  const Real CwT3 = 303.47;
    const Real CwT4 = -0.541;
    const auto WdT = (CwT1 + CwT2 * (wire_diameter / pin_diameter) +
                      CwT3 * Utility::pow<2>((wire_diameter / pin_diameter))) *
                     std::pow((wire_lead_length / pin_diameter), CwT4);
``` -/
def WdT (D Dw H : ℝ) : ℝ :=
  (19.56 + (-98.71) * (Dw / D) + 303.47 * (Dw / D) ^ 2) * (H / D) ^ (-0.541 : ℝ)

/-- ```
SCMFrictionUpgradedChengTodreas.C:206, :212
    CwL1 = 1.4;
    const auto WdL = CwL1 * WdT;
``` -/
def WdL (D Dw H : ℝ) : ℝ := 1.4 * WdT D Dw H

/-- Wire-wrap sweeping-flow coefficient (changed to log10 by the PR).
```
SCMFrictionUpgradedChengTodreas.C:214-217
    const Real a = -11;
    const Real b = 19;
    // Sweep coefficient
    const auto WsT = a * std::log10(wire_lead_length / pin_diameter) + b;
``` -/
def WsT (HoverD : ℝ) : ℝ := -11 * Real.logb 10 HoverD + 19

/-- ```
SCMFrictionUpgradedChengTodreas.C:207, :218
    CwL2 = 1.0;
    const auto WsL = CwL2 * WsT;
``` -/
def WsL (HoverD : ℝ) : ℝ := 1.0 * WsT HoverD

/-- Old version (parent commit a72db3da46ef1191033c6733b2da54c40cee144e, SCMFrictionUpdatedChengTodreas.C:191):
```
    const auto ws_t = -11.0 * std::log(wire_lead_length / pin_diameter) + 19.0;
``` -/
def WsTold (HoverD : ℝ) : ℝ := -11.0 * Real.log HoverD + 19.0

/-- ```
SCMFrictionUpgradedChengTodreas.C:227  ar = libMesh::pi * (pin_diameter + wire_diameter) * wire_diameter / 6.0;
SCMFrictionUpgradedChengTodreas.C:259  ar = libMesh::pi * (pin_diameter + wire_diameter) * wire_diameter / 6.0;
``` -/
def arCenterCorner (D Dw : ℝ) : ℝ := Real.pi * (D + Dw) * Dw / 6.0

/-- ```
SCMFrictionUpgradedChengTodreas.C:241  ar = libMesh::pi * (pin_diameter + wire_diameter) * wire_diameter / 4.0;
``` -/
def arEdge (D Dw : ℝ) : ℝ := Real.pi * (D + Dw) * Dw / 4.0

/-- ```
SCMFrictionUpgradedChengTodreas.C:229  a_p = S + libMesh::pi * Utility::pow<2>(wire_diameter) / 8.0 / std::cos(theta);
SCMFrictionUpgradedChengTodreas.C:243  a_p = S + libMesh::pi * Utility::pow<2>(wire_diameter) / 8.0 / std::cos(theta);
``` -/
def apCenterEdge (S Dw θ : ℝ) : ℝ := S + Real.pi * Dw ^ 2 / 8.0 / Real.cos θ

/-- ```
SCMFrictionUpgradedChengTodreas.C:261  a_p = S + libMesh::pi * Utility::pow<2>(wire_diameter) / 24.0 / std::cos(theta);
``` -/
def apCorner (S Dw θ : ℝ) : ℝ := S + Real.pi * Dw ^ 2 / 24.0 / Real.cos θ

/-- ```
SCMFrictionUpgradedChengTodreas.C:245-246 / :263-264
      const Real turbulent_wire_correction =
          1 + WsT * (ar / a_p) * Utility::pow<2>(std::tan(theta));
``` -/
def turbulentWireCorrection (Ws ar ap θ : ℝ) : ℝ := 1 + Ws * (ar / ap) * Real.tan θ ^ 2

def edgeCorr (D Dw H S : ℝ) : ℝ :=
  turbulentWireCorrection (WsT (H / D)) (arEdge D Dw) (apCenterEdge S Dw (theta H D Dw)) (theta H D Dw)

def cornerCorr (D Dw H S : ℝ) : ℝ :=
  turbulentWireCorrection (WsT (H / D)) (arCenterCorner D Dw) (apCorner S Dw (theta H D Dw))
    (theta H D Dw)

/-- Wire-wrap correction for the central subchannel.
```
SCMFrictionUpgradedChengTodreas.C:225-236
      const Real pw_p = libMesh::pi * pin_diameter / 2.0;
      CfT *= (pw_p / w_perim);
      CfT += WdT * (3.0 * ar / a_p) * (Dh_i / wire_lead_length) *
             std::pow((Dh_i / wire_diameter), 0.18);
      CfL *= (pw_p / w_perim);
      CfL += WdL * (3.0 * ar / a_p) * (Dh_i / wire_lead_length) * (Dh_i / wire_diameter);
``` -/
def wireCenter (D Dw H S wperim Dh CfL CfT : ℝ) : ℝ × ℝ :=
  (CfL * (Real.pi * D / 2.0 / wperim) +
      WdL D Dw H * (3.0 * arCenterCorner D Dw / apCenterEdge S Dw (theta H D Dw)) * (Dh / H) *
        (Dh / Dw),
    CfT * (Real.pi * D / 2.0 / wperim) +
      WdT D Dw H * (3.0 * arCenterCorner D Dw / apCenterEdge S Dw (theta H D Dw)) * (Dh / H) *
        (Dh / Dw) ^ (0.18 : ℝ))

/-- Edge subchannel (`mooseError` → `none`).
```
SCMFrictionUpgradedChengTodreas.C:247-254
      if (!std::isfinite(turbulent_wire_correction) || turbulent_wire_correction < 0.0)
        mooseError("The exponentiated term in the Cheng-Todreas turbulent wire correction must be "...);
      CfT *= std::pow(turbulent_wire_correction, 1.41);
      CfL *= (1 + WsL * (ar / a_p) * Utility::pow<2>(std::tan(theta)));
``` -/
def wireEdge (D Dw H S CfL CfT : ℝ) : Option (ℝ × ℝ) :=
  if edgeCorr D Dw H S < 0.0 then none
  else some (CfL * (1 + WsL (H / D) * (arEdge D Dw / apCenterEdge S Dw (theta H D Dw)) *
      Real.tan (theta H D Dw) ^ 2), CfT * edgeCorr D Dw H S ^ (1.41 : ℝ))

/-- Corner subchannel (:265-272, isomorphic to the edge subchannel). -/
def wireCorner (D Dw H S CfL CfT : ℝ) : Option (ℝ × ℝ) :=
  if cornerCorr D Dw H S < 0.0 then none
  else some (CfL * (1 + WsL (H / D) * (arCenterCorner D Dw / apCorner S Dw (theta H D Dw)) *
      Real.tan (theta H D Dw) ^ 2), CfT * cornerCorr D Dw H S ^ (1.41 : ℝ))

/-- Geometry (the values of `_subchannel_mesh` / `_tri_sch_mesh`; `hasWireWrap` corresponds to `_has_wire_wrap`, :28-29). -/
structure TriGeom where
  pitch : ℝ
  pinDiameter : ℝ
  wireLeadLength : ℝ
  wireDiameter : ℝ
  ductToPinGap : ℝ
  hasWireWrap : Bool

/-- `p_over_d = pitch / pin_diameter` (:101). -/
def pOverD (g : TriGeom) : ℝ := g.pitch / g.pinDiameter

/-- `w_over_d = (pin_diameter + gap) / pin_diameter` (:102). -/
def wOverD (g : TriGeom) : ℝ := (g.pinDiameter + g.ductToPinGap) / g.pinDiameter

/-- Coefficients `(CfL, CfT)` (including the wire-wrap correction; `Dh_i = 4.0 * S / w_perim`, :84). -/
def triCoeffs (g : TriGeom) (ty : EChannelType) (S wperim : ℝ) : Option (ℝ × ℝ) :=
  if g.hasWireWrap then
    match ty with
    | .CENTER =>
      some (wireCenter g.pinDiameter g.wireDiameter g.wireLeadLength S wperim (4.0 * S / wperim)
        (triBareCoeffs .CENTER (pOverD g) (wOverD g)).1 (triBareCoeffs .CENTER (pOverD g) (wOverD g)).2)
    | .EDGE =>
      wireEdge g.pinDiameter g.wireDiameter g.wireLeadLength S
        (triBareCoeffs .EDGE (pOverD g) (wOverD g)).1 (triBareCoeffs .EDGE (pOverD g) (wOverD g)).2
    | .CORNER =>
      wireCorner g.pinDiameter g.wireDiameter g.wireLeadLength S
        (triBareCoeffs .CORNER (pOverD g) (wOverD g)).1 (triBareCoeffs .CORNER (pOverD g) (wOverD g)).2
  else some (triBareCoeffs ty (pOverD g) (wOverD g))

/-- Triangular-lattice friction factor.
```
SCMFrictionUpgradedChengTodreas.C:80   const Real Re_eff = std::max(Re, 1.0);
SCMFrictionUpgradedChengTodreas.C:276-279
  const Real mL = -1.0;
  const Real mT = -0.18;
  auto fL = CfL * std::pow(Re_eff, mL);
  auto fT = CfT * std::pow(Re_eff, mT);
(piecewise structure: see regimeFriction)
``` -/
def triFrictionFactor (g : TriGeom) (ty : EChannelType) (Re S wperim Reb : ℝ) : Option ℝ :=
  (triCoeffs g ty S wperim).map fun c =>
    regimeFriction (c.1 * (max Re 1.0) ^ (-1.0 : ℝ)) (c.2 * (max Re 1.0) ^ (-0.18 : ℝ)) Reb
      (ReL (pOverD g)) (ReT (pOverD g))

/-! ## §3 Mixing coefficients, Borishanskii, overall Re -/

/-- ```
SCMMixingChengTodreas.C:110 / :228
  const Real ReL = 320.0 * std::pow(10.0, pitch / pin_diameter - 1.0);
``` -/
def ReLmix (pOverD : ℝ) : ℝ := 320.0 * (10.0 : ℝ) ^ (pOverD - 1.0)

/-- ```
SCMMixingChengTodreas.C:111 / :229
  const Real ReT = 10000.0 * std::pow(10.0, 0.7 * (pitch / pin_diameter - 1.0));
``` -/
def ReTmix (pOverD : ℝ) : ℝ := 10000.0 * (10.0 : ℝ) ^ (0.7 * (pOverD - 1.0))

/-- ```
SCMMixingChengTodreas.C:165-166 (Cs same as :284-285)
      const Real gamma = 2.0 / 3.0;
      Cm = CmL + (CmT - CmL) * std::pow(psi, gamma);
``` -/
def mixingBlend (CL CT ψ : ℝ) : ℝ := CL + (CT - CL) * ψ ^ (2.0 / 3.0 : ℝ)

/-- ```
SCMMixingChengTodreas.C:152-167 (Cs same as :271-286)
    if (bulk_Re < ReL) { Cm = CmL; }
    else if (bulk_Re > ReT) { Cm = CmT; }
    else { const Real psi = std::log(bulk_Re / ReL) / std::log(ReT / ReL); ... }
``` -/
def mixingRegime (CL CT bulkRe ReL ReT : ℝ) : ℝ :=
  if bulkRe < ReL then CL else if bulkRe > ReT then CT else mixingBlend CL CT (psi bulkRe ReL ReT)

theorem mixingRegime_eq_regime (CL CT bulkRe ReL ReT : ℝ) :
    mixingRegime CL CT bulkRe ReL ReT = regime CL CT (mixingBlend CL CT) bulkRe ReL ReT :=
  rfl

/-- ```
SCMHTCBorishanskii.C:45
  const Real poly = -8.12 + 12.76 * pre.poD - 3.65 * Utility::pow<2>(pre.poD);
``` -/
def borishanskiiPoly (poD : ℝ) : ℝ := -8.12 + 12.76 * poD - 3.65 * poD ^ 2

/-- Changed to log10 by the PR.
```
SCMHTCBorishanskii.C:52
  auto NuT = 24.15 * std::log10(poly);
``` -/
def borishanskiiNuBase (poD : ℝ) : ℝ := 24.15 * Real.logb 10 (borishanskiiPoly poD)

/-- Overall Re (computed on rank 0; the process guard `if (processor_id() != 0) return;` is not modeled here).
```
SubChannel1PhaseProblem.C:424-435
  for (const auto i_ch : make_range(_n_channels))
  { ... const Real mdot_in = (*_mdot_soln)(node_in);
    viscosity_in += mdot_in * (*_mu_soln)(node_in);
    mass_flow_in += mdot_in; }
  const Real bulk_Dh = _subchannel_mesh.getAssemblyHydraulicDiameter();
  const Real inlet_mu = viscosity_in / mass_flow_in;
  _bulk_Re = mass_flow_in * bulk_Dh / (inlet_mu * _subchannel_mesh.getAssemblyFlowArea());
``` -/
def bulkReynolds {ι : Type*} (S : Finset ι) (mdot mu : ι → ℝ) (Dh A : ℝ) : ℝ :=
  (∑ i ∈ S, mdot i) * Dh / ((∑ i ∈ S, mdot i * mu i) / (∑ i ∈ S, mdot i) * A)

/-! ## §4 Certificates: Re_L < Re_T, ψ ∈ [0,1], piecewise continuity -/

/-- **Certificate**: for P/D ∈ [1, 2] (covering the correlation's data domain [1.0, 1.42]), `0 < Re_L < Re_T`. -/
theorem ReL_pos_lt_ReT (pd : ℝ) (h1 : 1 ≤ pd) (h2 : pd ≤ 2) : 0 < ReL pd ∧ ReL pd < ReT pd := by
  unfold ReL ReT CbL1 CbL2 CbT1 CbT2
  have e1 : (1.0 : ℝ) = 1 := by norm_num
  have e7 : (0.7 : ℝ) = 7 / 10 := by norm_num
  rw [e1, e7, one_mul]
  have hA : (10 : ℝ) ^ (pd - 1) ≤ 10 := by
    calc (10 : ℝ) ^ (pd - 1) ≤ (10 : ℝ) ^ (1 : ℝ) :=
          Real.rpow_le_rpow_of_exponent_le (by norm_num) (by linarith)
      _ = 10 := Real.rpow_one 10
  have hApos : 0 < (10 : ℝ) ^ (pd - 1) := by positivity
  have hB : 1 ≤ (10 : ℝ) ^ (7 / 10 * (pd - 1)) := Real.one_le_rpow (by norm_num) (by linarith)
  constructor <;> linarith

/-- The Re_L formulas of the quadrilateral and triangular lattices are the same. -/
theorem ReL_quad_eq (pd : ℝ) : ReLquad pd = ReL pd := by
  unfold ReLquad ReL CbL1 CbL2
  have e1 : (1.0 : ℝ) = 1 := by norm_num
  have e320 : (320.0 : ℝ) = 320 := by norm_num
  rw [e1, e320, one_mul]

/-- The Re_T formulas of the quadrilateral and triangular lattices are the same. -/
theorem ReT_quad_eq (pd : ℝ) : ReTquad pd = ReT pd := by
  unfold ReTquad ReT CbT1 CbT2
  have e4 : (1.0E+4 : ℝ) = 10000 := by norm_num
  rw [e4]

/-- The Re_L, Re_T of the mixing closure are the same as those of the friction closure. -/
theorem ReLmix_eq (pd : ℝ) : ReLmix pd = ReL pd := by
  unfold ReLmix ReL CbL1 CbL2
  have e1 : (1.0 : ℝ) = 1 := by norm_num
  have e10 : (10.0 : ℝ) = 10 := by norm_num
  have e320 : (320.0 : ℝ) = 320 := by norm_num
  rw [e1, e10, e320, one_mul, mul_comm]

theorem ReTmix_eq (pd : ℝ) : ReTmix pd = ReT pd := by
  unfold ReTmix ReT CbT1 CbT2
  have e1 : (1.0 : ℝ) = 1 := by norm_num
  have e10 : (10.0 : ℝ) = 10 := by norm_num
  have e4 : (10000.0 : ℝ) = 10000 := by norm_num
  rw [e1, e10, e4, mul_comm]

/-- **Certificate**: within the transition branch `ψ ∈ [0, 1]` (hence all bases of `std::pow` are non-negative). -/
theorem psi_mem_Icc (Reb ReL ReT : ℝ) (hL : 0 < ReL) (hLT : ReL < ReT) (h1 : ReL ≤ Reb)
    (h2 : Reb ≤ ReT) : 0 ≤ psi Reb ReL ReT ∧ psi Reb ReL ReT ≤ 1 := by
  have hden : 0 < Real.log (ReT / ReL) := Real.log_pos ((one_lt_div hL).mpr hLT)
  have hb : 0 < Reb := lt_of_lt_of_le hL h1
  have hnum0 : 0 ≤ Real.log (Reb / ReL) := Real.log_nonneg ((one_le_div hL).mpr h1)
  have hq : Reb / ReL ≤ ReT / ReL := by
    rw [div_eq_mul_inv, div_eq_mul_inv]
    exact mul_le_mul_of_nonneg_right h2 (inv_nonneg.mpr hL.le)
  have hnum1 : Real.log (Reb / ReL) ≤ Real.log (ReT / ReL) := Real.log_le_log (div_pos hb hL) hq
  exact ⟨div_nonneg hnum0 hden.le, (div_le_one hden).mpr hnum1⟩

/-- The three-segment structure = `T ∘ clamp₀₁ ∘ ψ` (as long as `T 0 = a`, `T 1 = b`). -/
theorem regime_eq_clamp (a b : ℝ) (T : ℝ → ℝ) (ReL ReT r : ℝ) (h0 : T 0 = a) (h1 : T 1 = b)
    (hL : 0 < ReL) (hLT : ReL < ReT) (hr : 0 < r) :
    regime a b T r ReL ReT = T (max 0 (min 1 (psi r ReL ReT))) := by
  have hden : 0 < Real.log (ReT / ReL) := Real.log_pos ((one_lt_div hL).mpr hLT)
  unfold regime
  split_ifs with hA hB
  · have hnum : Real.log (r / ReL) < 0 := Real.log_neg (div_pos hr hL) ((div_lt_one hL).mpr hA)
    have hψ : psi r ReL ReT < 0 := div_neg_of_neg_of_pos hnum hden
    rw [min_eq_right (show psi r ReL ReT ≤ 1 by linarith), max_eq_left hψ.le, h0]
  · have hq : ReT / ReL < r / ReL := by
      rw [div_eq_mul_inv, div_eq_mul_inv]
      exact mul_lt_mul_of_pos_right hB (inv_pos.mpr hL)
    have hnum : Real.log (ReT / ReL) < Real.log (r / ReL) :=
      Real.log_lt_log (div_pos (by linarith) hL) hq
    have hψ : 1 < psi r ReL ReT := (one_lt_div hden).mpr hnum
    rw [min_eq_left hψ.le, max_eq_right (zero_le_one : (0 : ℝ) ≤ 1), h1]
  · have hψ := psi_mem_Icc r ReL ReT hL hLT (not_lt.mp hA) (not_lt.mp hB)
    rw [min_eq_right hψ.2, max_eq_right hψ.1]

/-- The three-segment structure is continuous on `(0, ∞)` (for any continuous transition function `T` with `T 0 = a`, `T 1 = b`). -/
theorem regime_continuousOn (a b : ℝ) (T : ℝ → ℝ) (ReL ReT : ℝ) (hT : Continuous T)
    (h0 : T 0 = a) (h1 : T 1 = b) (hL : 0 < ReL) (hLT : ReL < ReT) :
    ContinuousOn (fun r => regime a b T r ReL ReT) (Set.Ioi 0) := by
  have hdiv : ContinuousOn (fun r : ℝ => r / ReL) (Set.Ioi 0) := continuousOn_id.div_const ReL
  have hlog : ContinuousOn (fun r : ℝ => Real.log (r / ReL)) (Set.Ioi 0) :=
    hdiv.log (fun x hx => div_ne_zero (Set.mem_Ioi.mp hx).ne' hL.ne')
  have hpsi : ContinuousOn (fun r => psi r ReL ReT) (Set.Ioi 0) := by
    show ContinuousOn (fun r => Real.log (r / ReL) / Real.log (ReT / ReL)) (Set.Ioi 0)
    exact hlog.div_const _
  have hcl : Continuous (fun x : ℝ => max 0 (min 1 x)) :=
    continuous_const.max (continuous_const.min continuous_id)
  have hclamp : ContinuousOn (fun r => max 0 (min 1 (psi r ReL ReT))) (Set.Ioi 0) :=
    hcl.comp_continuousOn hpsi
  refine (hT.comp_continuousOn hclamp).congr (fun r hr => ?_)
  show regime a b T r ReL ReT = T (max 0 (min 1 (psi r ReL ReT)))
  exact regime_eq_clamp a b T ReL ReT r h0 h1 hL hLT (Set.mem_Ioi.mp hr)

theorem transitionBlend_zero (fL fT : ℝ) : transitionBlend fL fT 0 = fL := by
  unfold transitionBlend
  simp only [sub_zero, Real.one_rpow]
  rw [Real.zero_rpow (by norm_num : (7.0 : ℝ) ≠ 0),
    Real.zero_rpow (by norm_num : (1.0 / 3.0 : ℝ) ≠ 0)]
  ring

theorem transitionBlend_one (fL fT : ℝ) : transitionBlend fL fT 1 = fT := by
  unfold transitionBlend
  simp only [sub_self, Real.one_rpow]
  ring

theorem transitionBlend_continuous (fL fT : ℝ) : Continuous (transitionBlend fL fT) := by
  have h3 : Continuous (fun x : ℝ => x ^ (1.0 / 3.0 : ℝ)) := Real.continuous_rpow_const (by norm_num)
  have h7 : Continuous (fun x : ℝ => x ^ (7.0 : ℝ)) := Real.continuous_rpow_const (by norm_num)
  have hA : Continuous (fun ψ : ℝ => (1 - ψ) ^ (1.0 / 3.0 : ℝ)) :=
    h3.comp (continuous_const.sub continuous_id)
  show Continuous (fun ψ : ℝ => fL * (1 - ψ) ^ (1.0 / 3.0 : ℝ) * (1 - ψ ^ (7.0 : ℝ)) +
    fT * ψ ^ (1.0 / 3.0 : ℝ))
  exact ((continuous_const.mul hA).mul (continuous_const.sub h7)).add (continuous_const.mul h3)

/-- **Certificate (primary invariant B)**: for any `fL, fT` (i.e. any geometry, any subchannel Re), the friction factor as a function of the overall Re_b
is continuous on `(0, ∞)`, in particular at `Re_b = Re_L` and `Re_b = Re_T`. -/
theorem regimeFriction_continuousOn (fL fT ReL ReT : ℝ) (hL : 0 < ReL) (hLT : ReL < ReT) :
    ContinuousOn (fun Reb => regimeFriction fL fT Reb ReL ReT) (Set.Ioi 0) :=
  regime_continuousOn fL fT (transitionBlend fL fT) ReL ReT (transitionBlend_continuous fL fT)
    (transitionBlend_zero fL fT) (transitionBlend_one fL fT) hL hLT

/-- **Certificate**: for any triangular-lattice geometry (P/D ∈ [1,2]), subchannel type, subchannel Re and flow area: as long as the coefficient computation does not
go through `mooseError` (`triCoeffs = some c`), the friction factor is `some (regimeFriction …)`, and as a function of Re_b it is continuous on
`(0, ∞)`. -/
theorem triFrictionFactor_continuousOn_Reb (g : TriGeom) (ty : EChannelType) (Re S wperim : ℝ)
    (c : ℝ × ℝ) (hc : triCoeffs g ty S wperim = some c)
    (hpd1 : 1 ≤ pOverD g) (hpd2 : pOverD g ≤ 2) :
    (∀ Reb, triFrictionFactor g ty Re S wperim Reb =
      some (regimeFriction (c.1 * (max Re 1.0) ^ (-1.0 : ℝ)) (c.2 * (max Re 1.0) ^ (-0.18 : ℝ)) Reb
        (ReL (pOverD g)) (ReT (pOverD g)))) ∧
    ContinuousOn (fun Reb => regimeFriction (c.1 * (max Re 1.0) ^ (-1.0 : ℝ))
      (c.2 * (max Re 1.0) ^ (-0.18 : ℝ)) Reb (ReL (pOverD g)) (ReT (pOverD g))) (Set.Ioi 0) := by
  obtain ⟨hL, hLT⟩ := ReL_pos_lt_ReT (pOverD g) hpd1 hpd2
  refine ⟨fun Reb => ?_, regimeFriction_continuousOn _ _ _ _ hL hLT⟩
  unfold triFrictionFactor
  rw [hc]
  rfl

/-- **Certificate**: for `fL, fT ≥ 0`, `ψ ∈ [0,1]` the transition interpolation is non-negative. -/
theorem transitionBlend_nonneg (fL fT ψ : ℝ) (hfL : 0 ≤ fL) (hfT : 0 ≤ fT) (h0 : 0 ≤ ψ)
    (h1 : ψ ≤ 1) : 0 ≤ transitionBlend fL fT ψ := by
  unfold transitionBlend
  have h7 : ψ ^ (7.0 : ℝ) ≤ 1 := Real.rpow_le_one h0 h1 (by norm_num)
  exact add_nonneg
    (mul_nonneg (mul_nonneg hfL (Real.rpow_nonneg (by linarith) _)) (by linarith))
    (mul_nonneg hfT (Real.rpow_nonneg h0 _))

/-- **Certificate**: for `fL, fT ≥ 0`, the friction factor is non-negative for all `Re_b > 0`. -/
theorem regimeFriction_nonneg (fL fT Reb ReL ReT : ℝ) (hfL : 0 ≤ fL) (hfT : 0 ≤ fT)
    (hL : 0 < ReL) (hLT : ReL < ReT) (hr : 0 < Reb) : 0 ≤ regimeFriction fL fT Reb ReL ReT := by
  rw [regimeFriction_eq_regime,
    regime_eq_clamp fL fT (transitionBlend fL fT) ReL ReT Reb (transitionBlend_zero fL fT)
      (transitionBlend_one fL fT) hL hLT hr]
  exact transitionBlend_nonneg fL fT _ hfL hfT (le_max_left _ _)
    (max_le zero_le_one (min_le_left _ _))

/-! ## §5 Counterexample: the transition-region friction factor does not lie between the two ends -/

theorem cube_rpow_third (x : ℝ) (hx : 0 ≤ x) : (x ^ 3) ^ ((1 : ℝ) / 3) = x := by
  rw [← Real.rpow_natCast, ← Real.rpow_mul hx]
  norm_num

/-- **Counterexample (to "lies between the two ends within the transition region")**: the sum of the interpolation weights `(1−ψ)^{1/3}(1−ψ⁷)` and `ψ^{1/3}` exceeds 1,
so when the two ends are equal (`fL = fT = f₀ > 0`), the value at ψ = 1/8 is strictly greater than the ends. This is the shape of the correlation itself
(the same expression already existed in the parent commit), not introduced by this PR. -/
theorem transitionBlend_exceeds_ends (f0 : ℝ) (hf0 : 0 < f0) :
    f0 < transitionBlend f0 f0 (1 / 8) := by
  unfold transitionBlend
  have e13 : (1.0 / 3.0 : ℝ) = 1 / 3 := by norm_num
  have e7 : ((1 : ℝ) / 8) ^ (7.0 : ℝ) = 1 / 2097152 := by
    rw [show (7.0 : ℝ) = ((7 : ℕ) : ℝ) by norm_num, Real.rpow_natCast]
    norm_num
  have ha : ((1 : ℝ) / 8) ^ ((1 : ℝ) / 3) = 1 / 2 := by
    rw [show ((1 : ℝ) / 8) = (1 / 2) ^ 3 by norm_num]
    exact cube_rpow_third _ (by norm_num)
  have hb : (9 / 10 : ℝ) ≤ (1 - (1 : ℝ) / 8) ^ ((1 : ℝ) / 3) := by
    rw [← cube_rpow_third (9 / 10 : ℝ) (by norm_num)]
    exact Real.rpow_le_rpow (by norm_num) (by norm_num) (by norm_num)
  rw [e7, e13, ha]
  nlinarith [mul_le_mul_of_nonneg_left hb hf0.le]

/-! ## §6 Certificate: wire-wrap sweeping-flow coefficient and wire-wrap correction (primary invariant A) -/

/-- **Certificate**: for `0 < H/D ≤ 52` (including the data domain [8, 52]), `Ws_T = 19 − 11·log₁₀(H/D) > 0`.
The core is `52¹¹ < 10¹⁹`. -/
theorem WsT_pos (h : ℝ) (h0 : 0 < h) (h52 : h ≤ 52) : 0 < WsT h := by
  unfold WsT Real.logb
  have hl10 : 0 < Real.log 10 := Real.log_pos (by norm_num)
  have h1 : Real.log h ≤ Real.log 52 := Real.log_le_log h0 h52
  have h2 : Real.log ((52 : ℝ) ^ 11) < Real.log ((10 : ℝ) ^ 19) :=
    Real.log_lt_log (by norm_num) (by norm_num)
  rw [Real.log_pow, Real.log_pow] at h2
  push_cast at h2
  have h3 : Real.log h / Real.log 10 < 19 / 11 := by
    rw [div_lt_iff₀ hl10]
    linarith
  linarith

/-- **Counterexample (old version, ln)**: on the data domain `H/D ≥ 8` the old expression `19 − 11·ln(H/D)` is always negative. -/
theorem WsTold_neg (h : ℝ) (h8 : 8 ≤ h) : WsTold h < 0 := by
  unfold WsTold
  have hl2 := Real.log_two_gt_d9
  have h8' : Real.log 8 = 3 * Real.log 2 := by
    have := Real.log_pow (n := 3) (x := (2 : ℝ))
    norm_num at this
    exact this
  have hle : Real.log 8 ≤ Real.log h := Real.log_le_log (by norm_num) h8
  norm_num at hl2 ⊢
  linarith

theorem cos_theta (H D Dw : ℝ) (hH : 0 < H) :
    Real.cos (theta H D Dw) = H / Real.sqrt (H ^ 2 + (Real.pi * (D + Dw)) ^ 2) := by
  have hspos : 0 < Real.sqrt (H ^ 2 + (Real.pi * (D + Dw)) ^ 2) :=
    Real.sqrt_pos.mpr (add_pos_of_pos_of_nonneg (pow_pos hH 2) (sq_nonneg _))
  have hHs : H ≤ Real.sqrt (H ^ 2 + (Real.pi * (D + Dw)) ^ 2) :=
    calc H = Real.sqrt (H ^ 2) := (Real.sqrt_sq hH.le).symm
      _ ≤ Real.sqrt (H ^ 2 + (Real.pi * (D + Dw)) ^ 2) :=
          Real.sqrt_le_sqrt (le_add_of_nonneg_right (sq_nonneg _))
  unfold theta
  exact Real.cos_arccos (by linarith [div_nonneg hH.le hspos.le]) ((div_le_one hspos).mpr hHs)

theorem cos_theta_nonneg (H D Dw : ℝ) (hH : 0 < H) : 0 ≤ Real.cos (theta H D Dw) := by
  rw [cos_theta H D Dw hH]
  exact div_nonneg hH.le (Real.sqrt_nonneg _)

theorem turbulentWireCorrection_ge_one (Ws ar ap θ : ℝ) (hWs : 0 ≤ Ws) (har : 0 ≤ ar)
    (hap : 0 ≤ ap) : 1 ≤ turbulentWireCorrection Ws ar ap θ := by
  unfold turbulentWireCorrection
  have : 0 ≤ Ws * (ar / ap) * Real.tan θ ^ 2 :=
    mul_nonneg (mul_nonneg hWs (div_nonneg har hap)) (sq_nonneg _)
  linarith

/-- **Certificate**: edge subchannel, for `D > 0, D_w ≥ 0, H > 0, H/D ≤ 52, S > 0` the turbulent wire-wrap correction is ≥ 1. -/
theorem edgeCorr_ge_one (D Dw H S : ℝ) (hD : 0 < D) (hDw : 0 ≤ Dw) (hH : 0 < H)
    (hHD : H / D ≤ 52) (hS : 0 < S) : 1 ≤ edgeCorr D Dw H S := by
  unfold edgeCorr
  apply turbulentWireCorrection_ge_one
  · exact (WsT_pos (H / D) (div_pos hH hD) hHD).le
  · unfold arEdge
    exact div_nonneg (mul_nonneg (mul_nonneg Real.pi_pos.le (by linarith)) hDw) (by norm_num)
  · unfold apCenterEdge
    have : 0 ≤ Real.pi * Dw ^ 2 / 8.0 / Real.cos (theta H D Dw) :=
      div_nonneg (div_nonneg (mul_nonneg Real.pi_pos.le (sq_nonneg Dw)) (by norm_num))
        (cos_theta_nonneg H D Dw hH)
    linarith

/-- **Certificate**: corner subchannel, same as above. -/
theorem cornerCorr_ge_one (D Dw H S : ℝ) (hD : 0 < D) (hDw : 0 ≤ Dw) (hH : 0 < H)
    (hHD : H / D ≤ 52) (hS : 0 < S) : 1 ≤ cornerCorr D Dw H S := by
  unfold cornerCorr
  apply turbulentWireCorrection_ge_one
  · exact (WsT_pos (H / D) (div_pos hH hD) hHD).le
  · unfold arCenterCorner
    exact div_nonneg (mul_nonneg (mul_nonneg Real.pi_pos.le (by linarith)) hDw) (by norm_num)
  · unfold apCorner
    have : 0 ≤ Real.pi * Dw ^ 2 / 24.0 / Real.cos (theta H D Dw) :=
      div_nonneg (div_nonneg (mul_nonneg Real.pi_pos.le (sq_nonneg Dw)) (by norm_num))
        (cos_theta_nonneg H D Dw hH)
    linarith

/-- **Certificate**: the `mooseError` branch of the edge subchannel is unreachable within the applicable domain, and the wire wrap only increases the drag, never decreases it
(for `CfL, CfT ≥ 0` the corrected coefficients are ≥ the bare-rod coefficients). -/
theorem wireEdge_spec (D Dw H S CfL CfT : ℝ) (hD : 0 < D) (hDw : 0 ≤ Dw) (hH : 0 < H)
    (hHD : H / D ≤ 52) (hS : 0 < S) (hCfL : 0 ≤ CfL) (hCfT : 0 ≤ CfT) :
    ∃ c : ℝ × ℝ, wireEdge D Dw H S CfL CfT = some c ∧ CfL ≤ c.1 ∧ CfT ≤ c.2 := by
  have hc := edgeCorr_ge_one D Dw H S hD hDw hH hHD hS
  have hlam : 1 ≤ 1 + WsL (H / D) * (arEdge D Dw / apCenterEdge S Dw (theta H D Dw)) *
      Real.tan (theta H D Dw) ^ 2 := by
    unfold edgeCorr turbulentWireCorrection at hc
    unfold WsL
    have e1 : (1.0 : ℝ) = 1 := by norm_num
    rw [e1, one_mul]
    linarith
  have hnl : ¬ (edgeCorr D Dw H S < 0.0) := by
    have e0 : (0.0 : ℝ) = 0 := by norm_num
    rw [e0, not_lt]
    linarith
  refine ⟨(CfL * (1 + WsL (H / D) * (arEdge D Dw / apCenterEdge S Dw (theta H D Dw)) *
      Real.tan (theta H D Dw) ^ 2), CfT * edgeCorr D Dw H S ^ (1.41 : ℝ)), ?_, ?_, ?_⟩
  · unfold wireEdge
    rw [if_neg hnl]
  · exact le_mul_of_one_le_right hCfL hlam
  · exact le_mul_of_one_le_right hCfT (Real.one_le_rpow hc (by norm_num))

/-- **Certificate**: corner subchannel, same as above. -/
theorem wireCorner_spec (D Dw H S CfL CfT : ℝ) (hD : 0 < D) (hDw : 0 ≤ Dw) (hH : 0 < H)
    (hHD : H / D ≤ 52) (hS : 0 < S) (hCfL : 0 ≤ CfL) (hCfT : 0 ≤ CfT) :
    ∃ c : ℝ × ℝ, wireCorner D Dw H S CfL CfT = some c ∧ CfL ≤ c.1 ∧ CfT ≤ c.2 := by
  have hc := cornerCorr_ge_one D Dw H S hD hDw hH hHD hS
  have hlam : 1 ≤ 1 + WsL (H / D) * (arCenterCorner D Dw / apCorner S Dw (theta H D Dw)) *
      Real.tan (theta H D Dw) ^ 2 := by
    unfold cornerCorr turbulentWireCorrection at hc
    unfold WsL
    have e1 : (1.0 : ℝ) = 1 := by norm_num
    rw [e1, one_mul]
    linarith
  have hnl : ¬ (cornerCorr D Dw H S < 0.0) := by
    have e0 : (0.0 : ℝ) = 0 := by norm_num
    rw [e0, not_lt]
    linarith
  refine ⟨(CfL * (1 + WsL (H / D) * (arCenterCorner D Dw / apCorner S Dw (theta H D Dw)) *
      Real.tan (theta H D Dw) ^ 2), CfT * cornerCorr D Dw H S ^ (1.41 : ℝ)), ?_, ?_, ?_⟩
  · unfold wireCorner
    rw [if_neg hnl]
  · exact le_mul_of_one_le_right hCfL hlam
  · exact le_mul_of_one_le_right hCfT (Real.one_le_rpow hc (by norm_num))

/-! ## §7 Certificates: mixing-coefficient transition interpolation, Borishanskii, overall Re -/

theorem mixingBlend_zero (CL CT : ℝ) : mixingBlend CL CT 0 = CL := by
  unfold mixingBlend
  rw [Real.zero_rpow (by norm_num : (2.0 / 3.0 : ℝ) ≠ 0)]
  ring

theorem mixingBlend_one (CL CT : ℝ) : mixingBlend CL CT 1 = CT := by
  unfold mixingBlend
  rw [Real.one_rpow]
  ring

theorem mixingBlend_continuous (CL CT : ℝ) : Continuous (mixingBlend CL CT) := by
  show Continuous (fun ψ : ℝ => CL + (CT - CL) * ψ ^ (2.0 / 3.0 : ℝ))
  exact continuous_const.add (continuous_const.mul (Real.continuous_rpow_const (by norm_num)))

/-- **Certificate**: Cm (and Cs) as a function of the overall Re_b is continuous on `(0, ∞)` (no jump at Re_L and Re_T). -/
theorem mixingRegime_continuousOn (CL CT ReL ReT : ℝ) (hL : 0 < ReL) (hLT : ReL < ReT) :
    ContinuousOn (fun r => mixingRegime CL CT r ReL ReT) (Set.Ioi 0) :=
  regime_continuousOn CL CT (mixingBlend CL CT) ReL ReT (mixingBlend_continuous CL CT)
    (mixingBlend_zero CL CT) (mixingBlend_one CL CT) hL hLT

theorem mixingBlend_between (CL CT ψ : ℝ) (h0 : 0 ≤ ψ) (h1 : ψ ≤ 1) :
    min CL CT ≤ mixingBlend CL CT ψ ∧ mixingBlend CL CT ψ ≤ max CL CT := by
  unfold mixingBlend
  have hw0 : 0 ≤ ψ ^ (2.0 / 3.0 : ℝ) := Real.rpow_nonneg h0 _
  have hw1 : ψ ^ (2.0 / 3.0 : ℝ) ≤ 1 := Real.rpow_le_one h0 h1 (by norm_num)
  constructor
  · rcases le_total CL CT with h | h
    · rw [min_eq_left h]; nlinarith [mul_nonneg (sub_nonneg.2 h) hw0]
    · rw [min_eq_right h]; nlinarith [mul_nonneg (sub_nonneg.2 h) (sub_nonneg.2 hw1)]
  · rcases le_total CL CT with h | h
    · rw [max_eq_right h]; nlinarith [mul_nonneg (sub_nonneg.2 h) (sub_nonneg.2 hw1)]
    · rw [max_eq_left h]; nlinarith [mul_nonneg (sub_nonneg.2 h) hw0]

/-- **Certificate**: for any Re_b, Cm (and Cs) lies between the laminar value and the turbulent value (unlike the friction factor, see §5). -/
theorem mixingRegime_between (CL CT bulkRe ReL ReT : ℝ) (hL : 0 < ReL) (hLT : ReL < ReT) :
    min CL CT ≤ mixingRegime CL CT bulkRe ReL ReT ∧ mixingRegime CL CT bulkRe ReL ReT ≤ max CL CT := by
  unfold mixingRegime
  split_ifs with hA hB
  · exact ⟨min_le_left _ _, le_max_left _ _⟩
  · exact ⟨min_le_right _ _, le_max_right _ _⟩
  · have hψ := psi_mem_Icc bulkRe ReL ReT hL hLT (not_lt.mp hA) (not_lt.mp hB)
    exact mixingBlend_between CL CT _ hψ.1 hψ.2

/-- **Certificate**: on the Borishanskii applicable domain `P/D ∈ [1.1, 1.5]`, `poly > 1` (the mooseError for `poly <= 0` is unreachable). -/
theorem borishanskiiPoly_gt_one (x : ℝ) (h1 : 1.1 ≤ x) (h2 : x ≤ 1.5) : 1 < borishanskiiPoly x := by
  unfold borishanskiiPoly
  norm_num at h1 h2 ⊢
  nlinarith [mul_nonneg (sub_nonneg.2 h1) (sub_nonneg.2 h2)]

/-- **Certificate**: on the applicable domain the log10 version of the base Nu term is positive. -/
theorem borishanskiiNuBase_pos (x : ℝ) (h1 : 1.1 ≤ x) (h2 : x ≤ 1.5) : 0 < borishanskiiNuBase x := by
  unfold borishanskiiNuBase
  exact mul_pos (by norm_num) (Real.logb_pos (by norm_num) (borishanskiiPoly_gt_one x h1 h2))

/-- The old ln version of the base Nu term = ln 10 × the new log10 version (for any P/D). -/
theorem borishanskii_old_eq_ln10_mul_new (x : ℝ) :
    24.15 * Real.log (borishanskiiPoly x) = Real.log 10 * borishanskiiNuBase x := by
  unfold borishanskiiNuBase Real.logb
  have hne : Real.log 10 ≠ 0 := (Real.log_pos (by norm_num : (1 : ℝ) < 10)).ne'
  rw [show Real.log 10 * (24.15 * (Real.log (borishanskiiPoly x) / Real.log 10)) =
      24.15 * Real.log (borishanskiiPoly x) * (Real.log 10 / Real.log 10) by ring,
    div_self hne, mul_one]

/-- **Certificate**: when the inlet viscosity is uniform (`μ_i ≡ μ₀`) and the total flow rate is nonzero, `Re_b = ṁ D_h / (μ₀ A)`. -/
theorem bulkReynolds_uniform_mu {ι : Type*} (S : Finset ι) (mdot mu : ι → ℝ) (Dh A mu0 : ℝ)
    (hmu : ∀ i ∈ S, mu i = mu0) (hm : ∑ i ∈ S, mdot i ≠ 0) :
    bulkReynolds S mdot mu Dh A = (∑ i ∈ S, mdot i) * Dh / (mu0 * A) := by
  unfold bulkReynolds
  have hv : ∑ i ∈ S, mdot i * mu i = (∑ i ∈ S, mdot i) * mu0 := by
    rw [Finset.sum_mul]
    exact Finset.sum_congr rfl (fun i hi => by rw [hmu i hi])
  have hq : (∑ i ∈ S, mdot i) * mu0 / (∑ i ∈ S, mdot i) = mu0 := by
    rw [mul_comm, mul_div_assoc, div_self hm, mul_one]
  rw [hv, hq]

/-! ## Arithmetic (side note: constant checks, not counted as certificates) -/
section Arithmetic

/-- Side note: constant check, not counted as a certificate. The endpoint condition `52¹¹ < 10¹⁹` for `Ws_T > 0`. -/
theorem ws_endpoint_arith : (52 : ℝ) ^ 11 < 10 ^ 19 := by norm_num

/-- Side note: constant check, not counted as a certificate. `53¹¹ < 10¹⁹ < 54¹¹`: the zero 10^{19/11} ≈ 53.4 of the log10 version of Ws_T lies
between 53 and 54, right at the upper end 52 of the data domain. -/
theorem ws_root_arith : (53 : ℝ) ^ 11 < 10 ^ 19 ∧ (10 : ℝ) ^ 19 < 54 ^ 11 := by norm_num

end Arithmetic

end

end PhysAudit.Moose33307
