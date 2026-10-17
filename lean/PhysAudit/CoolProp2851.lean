/-
CoolProp2851: Lean transcription and certificates for CoolProp/CoolProp PR #2851
"fix(HEOS): correct first_saturation_deriv for mixtures (#2091)",
merge commit 54a2e32764bfa0ae81877e835c33ccd13ba8f853 (identical code in tag v8.0.0,
there at HelmholtzEOSMixtureBackend.cpp:3832-3891).

Scope
* Transcribed: `HelmholtzEOSMixtureBackend::calc_first_saturation_deriv` (both overloads,
  src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp:3646-3705) and the pressure
  derivative `MixtureDerivatives::dln_fugacity_coefficient_dp__constT_n`
  (src/Backends/Helmholtz/MixtureDerivatives.cpp:40-47).
* NOT transcribed (taken as abstract parameters / hypotheses): the Helmholtz-energy EOS,
  the fugacity coefficients ln φ_i of each phase, their partial derivatives, the VLE flash
  solver and the multivariable chain rule for ln φ_i along a curve (stated as `HasDerivAt`
  hypotheses of chain-rule shape), the Gibbs–Duhem relation of the EOS, and the residual
  partial-molar-enthalpy identity for (∂ ln φ_i/∂T)_{p,n}.

Main certificates
* `tangent_relation`: on any C¹ curve of states satisfying phase equilibrium with one phase
  at fixed composition z, the tangent (T', P') satisfies A·T' + B·P' = 0 where A, B are the
  weighted sums computed by the code loop (weights = incipient-phase composition).
* `bubble_branch_is_tangent`, `dew_branch_is_tangent`: the code's Q=0 / Q=1 value equals dT/dp
  along the bubble / dew curve.
* `bubble_dPdT_branch_is_tangent`: the `(iP, iT)` dispatch value equals dp/dT along the curve.
* `gernert_eq_generalized_clapeyron`, `single_component_matches_pure_branch`,
  `pure_limit_weights`, `mixture_formula_continuous`: generalized Clapeyron form and
  reduction to the unchanged pure-fluid Clapeyron branch.
-/
import Mathlib

open Filter Topology

noncomputable section

namespace CoolProp2851

variable {N : ℕ}

/-! ## Transcription of the C++ code -/

/-- Pure / pseudo-pure branch.
C++ (HelmholtzEOSMixtureBackend.cpp:3650-3651):
```
if (is_pure_or_pseudopure) {
    dTdP_sat = T() * (1 / SatV.rhomolar() - 1 / SatL.rhomolar()) / (SatV.hmolar() - SatL.hmolar());
```
-/
def dTdP_sat_pure (T rhoL rhoV hL hV : ℝ) : ℝ :=
  T * (1 / rhoV - 1 / rhoL) / (hV - hL)

/-- Loop accumulator of the mixture branch (one of `dQ_dPsat`, `dQ_dTsat`).
C++ (HelmholtzEOSMixtureBackend.cpp:3658, 3661-3667; same shape 3670-3676 with `x[i]`):
```
CoolPropDbl dQ_dPsat = 0, dQ_dTsat = 0;
for (std::size_t i = 0; i < N; ++i) {
    dQ_dPsat += y[i]
                * (MixtureDerivatives::dln_fugacity_coefficient_dp__constT_n(SatL, i, xN_flag)
                   - MixtureDerivatives::dln_fugacity_coefficient_dp__constT_n(SatV, i, xN_flag));
    dQ_dTsat += y[i]
                * (MixtureDerivatives::dln_fugacity_coefficient_dT__constp_n(SatL, i, xN_flag)
                   - MixtureDerivatives::dln_fugacity_coefficient_dT__constp_n(SatV, i, xN_flag));
}
```
`w` = weights, `dL i` = derivative evaluated on SatL, `dV i` = on SatV. The floating-point
summation order is replaced by the (commutative) real sum. -/
def weightedDiff (w dL dV : Fin N → ℝ) : ℝ :=
  ∑ i, w i * (dL i - dV i)

/-- Mixture branch; `none` models the `throw`.
C++ (HelmholtzEOSMixtureBackend.cpp:3656-3657, 3660, 3669, 3678-3681):
```
std::vector<CoolPropDbl> x = SatL.get_mole_fractions();
std::vector<CoolPropDbl> y = SatV.get_mole_fractions();
if (_Q == 0) { ... weights y[i] ... }
else if (_Q == 1) { ... weights x[i] ... }
else { throw ValueError(...); }
dTdP_sat = -dQ_dPsat / dQ_dTsat;
```
`dLp/dVp` = `dln_fugacity_coefficient_dp__constT_n` on SatL/SatV,
`dLT/dVT` = `dln_fugacity_coefficient_dT__constp_n` on SatL/SatV.
Note: C++ division by zero gives ±inf/NaN, real division gives 0; every theorem below that
uses the quotient assumes a nonzero denominator. -/
def dTdP_sat_mixture (Q : ℝ) (x y dLp dVp dLT dVT : Fin N → ℝ) : Option ℝ :=
  if Q = 0 then some (-weightedDiff y dLp dVp / weightedDiff y dLT dVT)
  else if Q = 1 then some (-weightedDiff x dLp dVp / weightedDiff x dLT dVT)
  else none

/-- Branch selection on `is_pure_or_pseudopure` (HelmholtzEOSMixtureBackend.cpp:3650-3682). -/
def dTdP_sat (isPure : Bool) (Q T rhoL rhoV hL hV : ℝ) (x y dLp dVp dLT dVT : Fin N → ℝ) :
    Option ℝ :=
  if isPure then some (dTdP_sat_pure T rhoL rhoV hL hV)
  else dTdP_sat_mixture Q x y dLp dVp dLT dVT

/-- The `parameters` enum, restricted to what the dispatch distinguishes. -/
inductive Param
  | iT
  | iP
  | other (k : ℕ)
  deriving DecidableEq

/-- Output dispatch.
C++ (HelmholtzEOSMixtureBackend.cpp:3684-3700):
```
if (Of1 == iT && Wrt1 == iP) { return dTdP_sat; }
else if (Of1 == iP && Wrt1 == iT) { return 1 / dTdP_sat; }
else if (Wrt1 == iT) { return first_partial_deriv(Of1, iT, iP) + first_partial_deriv(Of1, iP, iT) / dTdP_sat; }
else if (Wrt1 == iP) { return first_partial_deriv(Of1, iP, iT) + first_partial_deriv(Of1, iT, iP) * dTdP_sat; }
else { throw ValueError(...); }
```
`dYdT_p` = `first_partial_deriv(Of1, iT, iP)`, `dYdp_T` = `first_partial_deriv(Of1, iP, iT)`
(evaluated on `this`; not transcribed, abstract inputs). -/
def firstSatDeriv (Of1 Wrt1 : Param) (dTdP dYdT_p dYdp_T : ℝ) : Option ℝ :=
  if Of1 = Param.iT ∧ Wrt1 = Param.iP then some dTdP
  else if Of1 = Param.iP ∧ Wrt1 = Param.iT then some (1 / dTdP)
  else if Wrt1 = Param.iT then some (dYdT_p + dYdp_T / dTdP)
  else if Wrt1 = Param.iP then some (dYdp_T + dYdT_p * dTdP)
  else none

/-- C++ (MixtureDerivatives.cpp:40-47, GERG eq. 7.30):
```
CoolPropDbl R_u = HEOS.gas_constant();
double partial_molar_volumeval = partial_molar_volume(HEOS, i, xN_flag);
double term1 = partial_molar_volumeval / (R_u * HEOS._T);
double term2 = 1.0 / HEOS.p();
return term1 - term2;
```
`partial_molar_volume` itself is not transcribed (abstract `vbar`). -/
def dlnphi_dp (vbar R T p : ℝ) : ℝ :=
  vbar / (R * T) - 1 / p

/-! ## Structural facts about the transcription -/

/-- Pure branch is selected regardless of Q (behaviour for pure fluids unchanged). -/
theorem pure_branch_independent_of_Q (Q T rhoL rhoV hL hV : ℝ) (x y dLp dVp dLT dVT : Fin N → ℝ) :
    dTdP_sat true Q T rhoL rhoV hL hV x y dLp dVp dLT dVT
      = some (dTdP_sat_pure T rhoL rhoV hL hV) := by
  simp [dTdP_sat]

/-- Mixtures at intermediate quality: the code throws. -/
theorem mixture_intermediate_Q_throws (Q T rhoL rhoV hL hV : ℝ) (x y dLp dVp dLT dVT : Fin N → ℝ)
    (h0 : Q ≠ 0) (h1 : Q ≠ 1) :
    dTdP_sat false Q T rhoL rhoV hL hV x y dLp dVp dLT dVT = none := by
  simp [dTdP_sat, dTdP_sat_mixture, h0, h1]

lemma weightedDiff_swap (w a b : Fin N → ℝ) : weightedDiff w a b = -weightedDiff w b a := by
  unfold weightedDiff
  rw [← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  ring

/-! ## Derivation: tangent of the phase envelope -/

/-- Gibbs–Duhem relation for the incipient phase in a form valid for both mole-number
(`c = 0`) and mole-fraction coordinates: `∑_i w_i ∂ ln φ_i / ∂ w_j` does not depend on `j`.
This is a property of the (not transcribed) Helmholtz EOS; used as a hypothesis. -/
def GibbsDuhem (w : Fin N → ℝ) (G : Fin N → Fin N → ℝ) : Prop :=
  ∃ c : ℝ, ∀ j, ∑ i, w i * G i j = c

/-- Linear algebra step: weight the linearized equilibrium equations with `w` and sum. -/
theorem linearized_elimination (w w' FT FP GT GP : Fin N → ℝ) (G : Fin N → Fin N → ℝ)
    (T' P' : ℝ) (hw : ∀ i, w i ≠ 0)
    (hlin : ∀ i, (FT i * T' + FP i * P') - w' i / w i
      - (GT i * T' + GP i * P' + ∑ j, G i j * w' j) = 0)
    (hsum : ∑ i, w' i = 0) (hGD : GibbsDuhem w G) :
    weightedDiff w FT GT * T' + weightedDiff w FP GP * P' = 0 := by
  obtain ⟨c, hc⟩ := hGD
  have key : ∑ i, w i * ((FT i * T' + FP i * P') - w' i / w i
      - (GT i * T' + GP i * P' + ∑ j, G i j * w' j)) = 0 :=
    Finset.sum_eq_zero (fun i _ => by rw [hlin i, mul_zero])
  have hterm : ∀ i, w i * ((FT i * T' + FP i * P') - w' i / w i
      - (GT i * T' + GP i * P' + ∑ j, G i j * w' j))
      = w i * (FT i - GT i) * T' + w i * (FP i - GP i) * P' - w' i
        - w i * ∑ j, G i j * w' j := by
    intro i
    have hwi := hw i
    have e : w i * (w' i / w i) = w' i := by
      field_simp
    linear_combination (-1 : ℝ) * e
  have hcross : ∑ i, w i * ∑ j, G i j * w' j = c * ∑ j, w' j := by
    calc ∑ i, w i * ∑ j, G i j * w' j
        = ∑ i, ∑ j, w i * (G i j * w' j) := by simp_rw [Finset.mul_sum]
      _ = ∑ j, ∑ i, w i * (G i j * w' j) := Finset.sum_comm
      _ = ∑ j, (∑ i, w i * G i j) * w' j := by
          refine Finset.sum_congr rfl (fun j _ => ?_)
          rw [Finset.sum_mul]
          refine Finset.sum_congr rfl (fun i _ => ?_)
          ring
      _ = ∑ j, c * w' j := by simp_rw [hc]
      _ = c * ∑ j, w' j := by rw [Finset.mul_sum]
  have hexpand : ∑ i, w i * ((FT i * T' + FP i * P') - w' i / w i
      - (GT i * T' + GP i * P' + ∑ j, G i j * w' j))
      = weightedDiff w FT GT * T' + weightedDiff w FP GP * P' - ∑ i, w' i
        - ∑ i, w i * ∑ j, G i j * w' j := by
    rw [Finset.sum_congr rfl (fun i _ => hterm i)]
    simp only [weightedDiff, Finset.sum_add_distrib, Finset.sum_sub_distrib, Finset.sum_mul]
  rw [hexpand, hcross, hsum] at key
  linarith

/-- **Tangent relation (certificate).** Curve `s ↦ (Tc s, Pc s, wc s)` of equilibrium states:
the phase with fixed composition `z` has `ln φ_i = F i T p`, the incipient phase (composition
`wc s`) has `ln φ_i = G i T p w`. Equilibrium `f_i` equal ⇔
`log z_i + F_i = log w_i + G_i` (the `log p` terms cancel). Hypotheses `hF`, `hG` are the
chain rule for the (not transcribed) EOS with partials `FT, FP, GT, GP, GW`. Conclusion:
`(∑ w_i ΔT_i) T' + (∑ w_i Δp_i) P' = 0` with Δ = fixed − incipient, i.e. the code's
accumulators define the normal of the envelope. -/
theorem tangent_relation
    (F : Fin N → ℝ → ℝ → ℝ) (G : Fin N → ℝ → ℝ → (Fin N → ℝ) → ℝ)
    (z : Fin N → ℝ) (Tc Pc : ℝ → ℝ) (wc : ℝ → Fin N → ℝ) (s0 T' P' : ℝ) (w' : Fin N → ℝ)
    (FT FP GT GP : Fin N → ℝ) (GW : Fin N → Fin N → ℝ)
    (hw : ∀ i, HasDerivAt (fun s => wc s i) (w' i) s0)
    (hne : ∀ i, wc s0 i ≠ 0)
    (hF : ∀ i, HasDerivAt (fun s => F i (Tc s) (Pc s)) (FT i * T' + FP i * P') s0)
    (hG : ∀ i, HasDerivAt (fun s => G i (Tc s) (Pc s) (wc s))
      (GT i * T' + GP i * P' + ∑ j, GW i j * w' j) s0)
    (heq : ∀ᶠ s in 𝓝 s0, ∀ i,
      Real.log (z i) + F i (Tc s) (Pc s) = Real.log (wc s i) + G i (Tc s) (Pc s) (wc s))
    (hnorm : ∀ᶠ s in 𝓝 s0, ∑ i, wc s i = 1)
    (hGD : GibbsDuhem (wc s0) GW) :
    weightedDiff (wc s0) FT GT * T' + weightedDiff (wc s0) FP GP * P' = 0 := by
  have hsum : ∑ i, w' i = 0 := by
    have h1 : HasDerivAt (fun s => ∑ i, wc s i) (∑ i, w' i) s0 := by
      have h := HasDerivAt.sum (u := Finset.univ) (fun i _ => hw i)
      convert h using 1
      funext s
      simp [Finset.sum_apply]
    have h2 : HasDerivAt (fun s => ∑ i, wc s i) 0 s0 :=
      (hasDerivAt_const s0 (1 : ℝ)).congr_of_eventuallyEq hnorm
    exact h1.unique h2
  have hlin : ∀ i, (FT i * T' + FP i * P') - w' i / wc s0 i
      - (GT i * T' + GP i * P' + ∑ j, GW i j * w' j) = 0 := by
    intro i
    have hlog : HasDerivAt (fun s => Real.log (wc s i)) (w' i / wc s0 i) s0 :=
      (hw i).log (hne i)
    have hD : HasDerivAt
        (fun s => (Real.log (z i) + F i (Tc s) (Pc s))
          - (Real.log (wc s i) + G i (Tc s) (Pc s) (wc s)))
        ((0 + (FT i * T' + FP i * P'))
          - (w' i / wc s0 i + (GT i * T' + GP i * P' + ∑ j, GW i j * w' j))) s0 :=
      ((hasDerivAt_const s0 (Real.log (z i))).add (hF i)).sub (hlog.add (hG i))
    have hD0 : HasDerivAt
        (fun s => (Real.log (z i) + F i (Tc s) (Pc s))
          - (Real.log (wc s i) + G i (Tc s) (Pc s) (wc s))) 0 s0 :=
      (hasDerivAt_const s0 (0 : ℝ)).congr_of_eventuallyEq
        (heq.mono fun s hs => by simp only [hs i, sub_self])
    have := hD.unique hD0
    linarith
  exact linearized_elimination (wc s0) w' FT FP GT GP GW T' P' hne hlin hsum hGD

/-- **Bubble branch = dT/dp along the bubble curve (certificate).**
Liquid at fixed overall composition `z` (so SatL composition `x = z`), incipient vapor
`y = wc s0`; curve parametrized by pressure (`Pc' = 1`). Then the code's Q=0 value equals the
curve slope `T'`. `dLp, dLT` = liquid partials, `dVp, dVT` = vapor partials. -/
theorem bubble_branch_is_tangent
    (lnphiL : Fin N → ℝ → ℝ → ℝ) (lnphiV : Fin N → ℝ → ℝ → (Fin N → ℝ) → ℝ)
    (z : Fin N → ℝ) (Tc Pc : ℝ → ℝ) (yc : ℝ → Fin N → ℝ) (s0 T' : ℝ) (y' : Fin N → ℝ)
    (dLT dLp dVT dVp : Fin N → ℝ) (GW : Fin N → Fin N → ℝ)
    (hy : ∀ i, HasDerivAt (fun s => yc s i) (y' i) s0)
    (hne : ∀ i, yc s0 i ≠ 0)
    (hL : ∀ i, HasDerivAt (fun s => lnphiL i (Tc s) (Pc s)) (dLT i * T' + dLp i * 1) s0)
    (hV : ∀ i, HasDerivAt (fun s => lnphiV i (Tc s) (Pc s) (yc s))
      (dVT i * T' + dVp i * 1 + ∑ j, GW i j * y' j) s0)
    (heq : ∀ᶠ s in 𝓝 s0, ∀ i,
      Real.log (z i) + lnphiL i (Tc s) (Pc s) = Real.log (yc s i) + lnphiV i (Tc s) (Pc s) (yc s))
    (hnorm : ∀ᶠ s in 𝓝 s0, ∑ i, yc s i = 1)
    (hGD : GibbsDuhem (yc s0) GW)
    (hA : weightedDiff (yc s0) dLT dVT ≠ 0) :
    dTdP_sat_mixture 0 z (yc s0) dLp dVp dLT dVT = some T' := by
  have h := tangent_relation lnphiL lnphiV z Tc Pc yc s0 T' 1 y' dLT dLp dVT dVp GW
    hy hne hL hV heq hnorm hGD
  have hT : T' = -weightedDiff (yc s0) dLp dVp / weightedDiff (yc s0) dLT dVT := by
    rw [eq_div_iff hA]
    linear_combination h
  rw [hT]
  simp [dTdP_sat_mixture]

/-- **Dew branch = dT/dp along the dew curve (certificate).**
Vapor at fixed composition `z` (SatV composition `y = z`), incipient liquid `x = xc s0`;
equilibrium written `log x_i + ln φ_i^L = log z_i + ln φ_i^V`. -/
theorem dew_branch_is_tangent
    (lnphiL : Fin N → ℝ → ℝ → (Fin N → ℝ) → ℝ) (lnphiV : Fin N → ℝ → ℝ → ℝ)
    (z : Fin N → ℝ) (Tc Pc : ℝ → ℝ) (xc : ℝ → Fin N → ℝ) (s0 T' : ℝ) (x' : Fin N → ℝ)
    (dLT dLp dVT dVp : Fin N → ℝ) (GW : Fin N → Fin N → ℝ)
    (hx : ∀ i, HasDerivAt (fun s => xc s i) (x' i) s0)
    (hne : ∀ i, xc s0 i ≠ 0)
    (hV : ∀ i, HasDerivAt (fun s => lnphiV i (Tc s) (Pc s)) (dVT i * T' + dVp i * 1) s0)
    (hL : ∀ i, HasDerivAt (fun s => lnphiL i (Tc s) (Pc s) (xc s))
      (dLT i * T' + dLp i * 1 + ∑ j, GW i j * x' j) s0)
    (heq : ∀ᶠ s in 𝓝 s0, ∀ i,
      Real.log (xc s i) + lnphiL i (Tc s) (Pc s) (xc s) = Real.log (z i) + lnphiV i (Tc s) (Pc s))
    (hnorm : ∀ᶠ s in 𝓝 s0, ∑ i, xc s i = 1)
    (hGD : GibbsDuhem (xc s0) GW)
    (hA : weightedDiff (xc s0) dLT dVT ≠ 0) :
    dTdP_sat_mixture 1 (xc s0) z dLp dVp dLT dVT = some T' := by
  have h := tangent_relation lnphiV lnphiL z Tc Pc xc s0 T' 1 x' dVT dVp dLT dLp GW
    hx hne hV hL (heq.mono fun s hs i => (hs i).symm) hnorm hGD
  rw [weightedDiff_swap (xc s0) dVT dLT, weightedDiff_swap (xc s0) dVp dLp] at h
  have hT : T' = -weightedDiff (xc s0) dLp dVp / weightedDiff (xc s0) dLT dVT := by
    rw [eq_div_iff hA]
    linear_combination (-1 : ℝ) * h
  rw [hT]
  simp [dTdP_sat_mixture]

/-- **`(iP, iT)` output = dp/dT along the bubble curve (certificate).** Curve parametrized by
temperature (`Tc' = 1`); the dispatch returns `1 / dTdP_sat` with the code's Q=0 value. -/
theorem bubble_dPdT_branch_is_tangent
    (lnphiL : Fin N → ℝ → ℝ → ℝ) (lnphiV : Fin N → ℝ → ℝ → (Fin N → ℝ) → ℝ)
    (z : Fin N → ℝ) (Tc Pc : ℝ → ℝ) (yc : ℝ → Fin N → ℝ) (s0 P' : ℝ) (y' : Fin N → ℝ)
    (dLT dLp dVT dVp : Fin N → ℝ) (GW : Fin N → Fin N → ℝ) (dYdT_p dYdp_T : ℝ)
    (hy : ∀ i, HasDerivAt (fun s => yc s i) (y' i) s0)
    (hne : ∀ i, yc s0 i ≠ 0)
    (hL : ∀ i, HasDerivAt (fun s => lnphiL i (Tc s) (Pc s)) (dLT i * 1 + dLp i * P') s0)
    (hV : ∀ i, HasDerivAt (fun s => lnphiV i (Tc s) (Pc s) (yc s))
      (dVT i * 1 + dVp i * P' + ∑ j, GW i j * y' j) s0)
    (heq : ∀ᶠ s in 𝓝 s0, ∀ i,
      Real.log (z i) + lnphiL i (Tc s) (Pc s) = Real.log (yc s i) + lnphiV i (Tc s) (Pc s) (yc s))
    (hnorm : ∀ᶠ s in 𝓝 s0, ∑ i, yc s i = 1)
    (hGD : GibbsDuhem (yc s0) GW)
    (hB : weightedDiff (yc s0) dLp dVp ≠ 0) :
    firstSatDeriv Param.iP Param.iT
      (-weightedDiff (yc s0) dLp dVp / weightedDiff (yc s0) dLT dVT) dYdT_p dYdp_T = some P' := by
  have h := tangent_relation lnphiL lnphiV z Tc Pc yc s0 1 P' y' dLT dLp dVT dVp GW
    hy hne hL hV heq hnorm hGD
  have hP : P' = 1 / (-weightedDiff (yc s0) dLp dVp / weightedDiff (yc s0) dLT dVT) := by
    rw [one_div_div, eq_div_iff (neg_ne_zero.mpr hB)]
    linear_combination (-1 : ℝ) * h
  rw [hP]
  simp [firstSatDeriv]

/-! ## Consistency of the output dispatch -/

/-- `(iP,iT)` and `(iT,iP)` outputs are reciprocal. -/
theorem firstSatDeriv_reciprocal (D a b c d : ℝ) (hD : D ≠ 0) :
    ∃ u v, firstSatDeriv Param.iT Param.iP D a b = some u ∧
      firstSatDeriv Param.iP Param.iT D c d = some v ∧ u * v = 1 := by
  refine ⟨D, 1 / D, by simp [firstSatDeriv], by simp [firstSatDeriv], ?_⟩
  field_simp

/-- For a generic output `Y`, the `Wrt1 = iT` value is the `Wrt1 = iP` value divided by
`dTdP_sat` (chain rule dY/dT|σ = (dY/dp|σ)/(dT/dp|σ)). -/
theorem firstSatDeriv_wrtT_eq_wrtP_div (k : ℕ) (D a b : ℝ) (hD : D ≠ 0) :
    ∃ u v, firstSatDeriv (Param.other k) Param.iT D a b = some u ∧
      firstSatDeriv (Param.other k) Param.iP D a b = some v ∧ u = v / D := by
  refine ⟨a + b / D, b + a * D, by simp [firstSatDeriv], by simp [firstSatDeriv], ?_⟩
  field_simp
  try ring

/-- For a generic output `Y(T,p)` with partials `a = (∂Y/∂T)_p`, `b = (∂Y/∂p)_T`, the
`Wrt1 = iP` value equals d/ds Y(Tc s, Pc s) on a pressure-parametrized curve with slope
`T' = dTdP_sat` (chain rule taken as hypothesis `hY`). -/
theorem firstSatDeriv_wrtP_is_total_derivative (k : ℕ) (Y : ℝ → ℝ → ℝ) (Tc Pc : ℝ → ℝ)
    (s0 T' a b : ℝ) (hY : HasDerivAt (fun s => Y (Tc s) (Pc s)) (a * T' + b * 1) s0) :
    ∃ v, firstSatDeriv (Param.other k) Param.iP T' a b = some v ∧
      HasDerivAt (fun s => Y (Tc s) (Pc s)) v s0 := by
  refine ⟨b + a * T', by simp [firstSatDeriv], ?_⟩
  convert hY using 1
  ring

/-! ## Generalized Clapeyron form and pure-fluid limit -/

/-- **Generalized Clapeyron (certificate).** With the code's GERG 7.30 pressure derivative
(`dlnphi_dp`, the `1/p` terms cancel) and the thermodynamic identity
`(∂ ln φ_i/∂T)_{p,n} = -(h̄_i - h_i^{ig})/(R T²)` (hypothesis layer: not transcribed; the
ideal-gas part `hig` is common to both phases at equal T), the code's mixture quotient equals
`T ∑ w_i (v̄_i^V - v̄_i^L) / ∑ w_i (h̄_i^V - h̄_i^L)`. -/
theorem gernert_eq_generalized_clapeyron (w vL vV hL hV hig : Fin N → ℝ) (R T p : ℝ)
    (hR : R ≠ 0) (hT : T ≠ 0) (hSH : ∑ i, w i * (hV i - hL i) ≠ 0) :
    -weightedDiff w (fun i => dlnphi_dp (vL i) R T p) (fun i => dlnphi_dp (vV i) R T p)
      / weightedDiff w (fun i => -(hL i - hig i) / (R * T ^ 2))
          (fun i => -(hV i - hig i) / (R * T ^ 2))
    = T * (∑ i, w i * (vV i - vL i)) / ∑ i, w i * (hV i - hL i) := by
  have hnum : weightedDiff w (fun i => dlnphi_dp (vL i) R T p) (fun i => dlnphi_dp (vV i) R T p)
      = -(∑ i, w i * (vV i - vL i)) / (R * T) := by
    simp only [weightedDiff, dlnphi_dp]
    rw [neg_div, Finset.sum_div, ← Finset.sum_neg_distrib]
    refine Finset.sum_congr rfl (fun i _ => ?_)
    ring
  have hden : weightedDiff w (fun i => -(hL i - hig i) / (R * T ^ 2))
      (fun i => -(hV i - hig i) / (R * T ^ 2))
      = (∑ i, w i * (hV i - hL i)) / (R * T ^ 2) := by
    simp only [weightedDiff]
    rw [Finset.sum_div]
    refine Finset.sum_congr rfl (fun i _ => ?_)
    ring
  rw [hnum, hden]
  generalize ∑ i, w i * (vV i - vL i) = SV
  generalize hsh : ∑ i, w i * (hV i - hL i) = SH at hSH ⊢
  field_simp
  try ring

/-- **Single component reproduces the unchanged pure branch (certificate).** For `N = 1`
(weight 1, partial molar = molar: `v̄ = 1/ρ`, `h̄ = h`) the Gernert quotient with the code's
`dlnphi_dp` equals the code's Clapeyron expression `dTdP_sat_pure`. -/
theorem single_component_matches_pure_branch (R T p rhoL rhoV hL hV hig : ℝ)
    (hR : R ≠ 0) (hT : T ≠ 0) (hh : hV ≠ hL) :
    -weightedDiff (fun _ : Fin 1 => (1 : ℝ))
        (fun _ => dlnphi_dp (1 / rhoL) R T p) (fun _ => dlnphi_dp (1 / rhoV) R T p)
      / weightedDiff (fun _ : Fin 1 => (1 : ℝ))
          (fun _ => -(hL - hig) / (R * T ^ 2)) (fun _ => -(hV - hig) / (R * T ^ 2))
    = dTdP_sat_pure T rhoL rhoV hL hV := by
  have hSH : ∑ i : Fin 1, (fun _ => (1 : ℝ)) i * ((fun _ => hV) i - (fun _ => hL) i) ≠ 0 := by
    simpa using sub_ne_zero.mpr hh
  have := gernert_eq_generalized_clapeyron (fun _ : Fin 1 => (1 : ℝ)) (fun _ => 1 / rhoL)
    (fun _ => 1 / rhoV) (fun _ => hL) (fun _ => hV) (fun _ => hig) R T p hR hT hSH
  rw [this]
  simp only [Fin.sum_univ_one, one_mul, dTdP_sat_pure]

/-- Weights concentrated on component 0 (the pure-fluid limit of the incipient-phase
composition) select component 0. -/
theorem pure_limit_weights {n : ℕ} (a b : Fin (n + 1) → ℝ) :
    weightedDiff (Pi.single (0 : Fin (n + 1)) (1 : ℝ)) a b = a 0 - b 0 := by
  unfold weightedDiff
  rw [Finset.sum_eq_single (0 : Fin (n + 1))]
  · simp
  · intro j _ hj
    simp [Pi.single_apply, hj]
  · intro h
    exact absurd (Finset.mem_univ _) h

/-- **Continuity of the mixture quotient (certificate, conditional on continuity of the
flash solution and EOS partials, hypotheses `h*`).** Used with `pure_limit_weights` and
`single_component_matches_pure_branch` for the limit z₂ → 0. -/
theorem mixture_formula_continuous {ι : Type*} {l : Filter ι}
    (w dLp dVp dLT dVT : ι → Fin N → ℝ) (w0 aLp aVp aLT aVT : Fin N → ℝ)
    (hw : ∀ i, Tendsto (fun k => w k i) l (𝓝 (w0 i)))
    (hLp : ∀ i, Tendsto (fun k => dLp k i) l (𝓝 (aLp i)))
    (hVp : ∀ i, Tendsto (fun k => dVp k i) l (𝓝 (aVp i)))
    (hLT : ∀ i, Tendsto (fun k => dLT k i) l (𝓝 (aLT i)))
    (hVT : ∀ i, Tendsto (fun k => dVT k i) l (𝓝 (aVT i)))
    (hden : weightedDiff w0 aLT aVT ≠ 0) :
    Tendsto (fun k => -weightedDiff (w k) (dLp k) (dVp k) / weightedDiff (w k) (dLT k) (dVT k))
      l (𝓝 (-weightedDiff w0 aLp aVp / weightedDiff w0 aLT aVT)) := by
  have tw : ∀ (a b : ι → Fin N → ℝ) (a0 b0 : Fin N → ℝ),
      (∀ i, Tendsto (fun k => a k i) l (𝓝 (a0 i))) →
      (∀ i, Tendsto (fun k => b k i) l (𝓝 (b0 i))) →
      Tendsto (fun k => weightedDiff (w k) (a k) (b k)) l (𝓝 (weightedDiff w0 a0 b0)) := by
    intro a b a0 b0 ha hb
    unfold weightedDiff
    exact tendsto_finset_sum _ (fun i _ => (hw i).mul ((ha i).sub (hb i)))
  exact ((tw dLp dVp aLp aVp hLp hVp).neg).div (tw dLT dVT aLT aVT hLT hVT) hden

end CoolProp2851

end
