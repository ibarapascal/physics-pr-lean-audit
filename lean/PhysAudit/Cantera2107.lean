-- Cantera2107: Cantera/cantera PR #2107 "Fix thermo consistency errors"
-- (merge commit d58a23265414e1789536cf5332d9fb29de095a43; PR base a6feb655151a3740f83703be87522322855f4a75) Lean certificate.
-- §1 General layer: Fréchet derivatives of abstract g(T,P), s, v + Maxwell symmetry ⇒ Mayer relation (∂u/∂T)_v = c_p − T v β²/κ_T;
--    when κ_T = 0 and β ≠ 0, no differentiable isochoric path exists.
-- §2 MargulesVPSSTP: transcription of ln γ, μ_k, v̄_k, dlnγ/dT, d²lnγ/dT², s̄_k, h̄_k, β, κ_T, cv;
--    ∂μ_k/∂P = v̄_k, −∂μ_k/∂T = s̄_k, temperature derivatives are true derivatives, cv = (∂u/∂T)_v; counterexample for the pre-PR version (without the (P−P_ref)V^E term).
-- §3 IdealSolnGasVPSS: transcription of cv/β/κ_T; the Mayer test expression is identical to the code expression for any V ≠ 0; with V = Σx_k v°_k, cv = (∂u/∂T)_v.
-- §4 PDSS_SSVol (tpoly): c_p,ss/R = T ∂(s_ss/R)/∂T; counterexample for the pre-PR version missing the /R.
import Mathlib

/-!
# Cantera #2107 — thermo consistency

Transcription conventions
* Floating point → reals; `std::max`, `if`, `fabs(·) < 1e-10` are transcribed as written. Accumulation loops `x[k] += …` are written as `Finset` sums
  (real addition is commutative and associative, so the order does not matter).
* The standard-state (PDSS) layer is not transcribed: μ°_k(T,P), v°_k(T,P), s°_k/R, h°_k/RT, PDSS::dVdT/dVdP enter as abstract function parameters,
  and the required relations (∂μ°/∂P = v°, ∂μ°/∂T = −s°, h° = μ° + T s°, dVdT/dVdP are true derivatives of v°) are written as explicit hypotheses.
* g, s, v of the general layer are abstract functions of (T,P); the Maxwell symmetry (∂s/∂P)_T = −(∂v/∂T)_P (which follows from g ∈ C² and Schwarz's theorem)
  is taken as the hypothesis `hMaxwell` and is not derived from g in Lean.
* `molarVolume()`: in VPStandardStateTP.cpp:199-204, density = meanMW / Σ_k X_k v̄_k and molarVolume = meanMW / density,
  which round-trip to an identity over the reals (meanMW ≠ 0), so Σ_k X_k v̄_k is used directly.
-/

namespace PhysAudit.Cantera2107

noncomputable section

open Filter Topology

/-! ## §0 Constants (include/cantera/base/ct_defs.h @ d58a232) -/

-- ct_defs.h:84   const double Avogadro = 6.02214076e26;
def Avogadro : ℝ := 6.02214076e26
-- ct_defs.h:87   const double Boltzmann = 1.380649e-23;
def Boltzmann : ℝ := 1.380649e-23
-- ct_defs.h:123  const double GasConstant = Avogadro * Boltzmann;
def GasConstant : ℝ := Avogadro * Boltzmann
-- ct_defs.h:99   const double OneAtm = 1.01325e5;
def OneAtm : ℝ := 1.01325e5
-- ct_defs.h:161  const double SmallNumber = 1.e-300;
def SmallNumber : ℝ := 1e-300

/-- ThermoPhase.h:1201-1203
```
double RT() const {
    return temperature() * GasConstant;
}
``` -/
def RT (T : ℝ) : ℝ := T * GasConstant

theorem GasConstant_pos : 0 < GasConstant := by
  unfold GasConstant Avogadro Boltzmann; norm_num

theorem GasConstant_ne_one : GasConstant ≠ 1 := by
  unfold GasConstant Avogadro Boltzmann; norm_num

/-- Helper: a / V * V = a. -/
theorem div_mul_self' {a V : ℝ} (hV : V ≠ 0) : a / V * V = a := by
  first
  | (field_simp; done)
  | exact div_mul_cancel₀ a hV
  | (field_simp <;> ring)

/-- Derivative of the identity function (fallback for the explicit/implicit argument differences of `hasDerivAt_id'`). -/
theorem hasDerivAt_idfun (x : ℝ) : HasDerivAt (fun y : ℝ => y) 1 x := by
  first
  | exact hasDerivAt_id' x
  | exact hasDerivAt_id x
  | exact hasDerivAt_id'

/-! ## §1 General layer (abstract Gibbs function; this layer transcribes no code and is taken as hypotheses) -/

/-- A continuous linear functional on ℝ×ℝ expanded in coordinates. -/
theorem clm_pair (L : ℝ × ℝ →L[ℝ] ℝ) (a b : ℝ) :
    L (a, b) = a * L (1, 0) + b * L (0, 1) := by
  have h : ((a, b) : ℝ × ℝ) = a • ((1 : ℝ), (0 : ℝ)) + b • ((0 : ℝ), (1 : ℝ)) := by
    ext <;> simp
  rw [h]
  first
  | (rw [map_add, map_smul, map_smul, smul_eq_mul, smul_eq_mul])
  | (simp only [map_add, map_smul, smul_eq_mul])
  | (simp only [map_add, ContinuousLinearMap.map_smul, smul_eq_mul])

/-- Derivative of the curve t ↦ (t, p t). -/
theorem curve_hasDerivAt {p : ℝ → ℝ} {p' T0 : ℝ} (hp : HasDerivAt p p' T0) :
    HasDerivAt (fun t => (t, p t)) ((1 : ℝ), p') T0 := by
  first
  | exact (hasDerivAt_idfun T0).prodMk hp
  | exact (hasDerivAt_idfun T0).prod hp
  | exact HasDerivAt.prodMk (hasDerivAt_idfun T0) hp

/-- Chain rule along the curve t ↦ (t, p t). -/
theorem comp_curve {F : ℝ × ℝ → ℝ} {L : ℝ × ℝ →L[ℝ] ℝ} {p : ℝ → ℝ} {p' T0 : ℝ}
    (hF : HasFDerivAt F L (T0, p T0)) (hp : HasDerivAt p p' T0) :
    HasDerivAt (fun t => F (t, p t)) (L ((1 : ℝ), p')) T0 := by
  have hγ := curve_hasDerivAt hp
  first
  | exact hF.comp_hasDerivAt T0 hγ
  | exact HasFDerivAt.comp_hasDerivAt T0 hF hγ

/-- Take the T partial derivative from the Fréchet derivative. -/
theorem partialT {F : ℝ × ℝ → ℝ} {L : ℝ × ℝ →L[ℝ] ℝ} {T P : ℝ}
    (h : HasFDerivAt F L (T, P)) : HasDerivAt (fun t => F (t, P)) (L (1, 0)) T :=
  comp_curve (p := fun _ => P) h (hasDerivAt_const T P)

/-- Take the P partial derivative from the Fréchet derivative. -/
theorem partialP {F : ℝ × ℝ → ℝ} {L : ℝ × ℝ →L[ℝ] ℝ} {T P : ℝ}
    (h : HasFDerivAt F L (T, P)) : HasDerivAt (fun q => F (T, q)) (L (0, 1)) P := by
  have hγ : HasDerivAt (fun q : ℝ => (T, q)) ((0 : ℝ), (1 : ℝ)) P := by
    first
    | exact (hasDerivAt_const P T).prodMk (hasDerivAt_idfun P)
    | exact (hasDerivAt_const P T).prod (hasDerivAt_idfun P)
    | exact HasDerivAt.prodMk (hasDerivAt_const P T) (hasDerivAt_idfun P)
  first
  | exact h.comp_hasDerivAt P hγ
  | exact HasFDerivAt.comp_hasDerivAt P h hγ

/-- Slope of the isochoric path: along t ↦ (t, p t), if v is constant then (∂v/∂T) + p' (∂v/∂P) = 0. -/
theorem isochore_slope (v : ℝ × ℝ → ℝ) (Lv : ℝ × ℝ →L[ℝ] ℝ) (T0 : ℝ) (p : ℝ → ℝ) (p' : ℝ)
    (hv : HasFDerivAt v Lv (T0, p T0)) (hp : HasDerivAt p p' T0)
    (hiso : ∀ᶠ t in 𝓝 T0, v (t, p t) = v (T0, p T0)) :
    Lv ((1 : ℝ), p') = 0 := by
  have dv := comp_curve hv hp
  have h0 : HasDerivAt (fun t => v (t, p t)) 0 T0 :=
    (hasDerivAt_const T0 (v (T0, p T0))).congr_of_eventuallyEq hiso
  exact dv.unique h0

/-- Constant-volume entropy derivative: (∂s/∂T)_v = (∂s/∂T)_P + (∂v/∂T)_P² / (∂v/∂P)_T, hence c_v = T(∂s/∂T)_v. -/
theorem isochoric_entropy (s v : ℝ × ℝ → ℝ) (Ls Lv : ℝ × ℝ →L[ℝ] ℝ) (T0 : ℝ)
    (p : ℝ → ℝ) (p' : ℝ)
    (hs : HasFDerivAt s Ls (T0, p T0)) (hv : HasFDerivAt v Lv (T0, p T0))
    (hMaxwell : Ls (0, 1) = -Lv (1, 0)) (he : Lv (0, 1) ≠ 0)
    (hp : HasDerivAt p p' T0) (hiso : ∀ᶠ t in 𝓝 T0, v (t, p t) = v (T0, p T0)) :
    HasDerivAt (fun t => s (t, p t)) (Ls (1, 0) + Lv (1, 0) ^ 2 / Lv (0, 1)) T0 := by
  have ds := comp_curve hs hp
  have hlin := isochore_slope v Lv T0 p p' hv hp hiso
  rw [clm_pair] at hlin
  have hp' : p' = -Lv (1, 0) / Lv (0, 1) := by
    rw [eq_div_iff he]; linarith
  refine ds.congr_deriv ?_
  rw [clm_pair Ls 1 p', hMaxwell, hp']
  field_simp <;> ring

/-- **Mayer relation (general theorem)**. With g the Gibbs function (Lg(1,0) = −s, Lg(0,1) = v) and u = g + T s − P v,
along a differentiable isochoric path t ↦ (t, p t) we have du/dT = c_p − T v β²/κ_T, where c_p = T(∂s/∂T)_P,
β = (∂v/∂T)_P / v, κ_T = −(∂v/∂P)_T / v. -/
theorem mayer_general (g s v : ℝ × ℝ → ℝ) (Lg Ls Lv : ℝ × ℝ →L[ℝ] ℝ) (T0 : ℝ)
    (p : ℝ → ℝ) (p' : ℝ)
    (hg : HasFDerivAt g Lg (T0, p T0)) (hs : HasFDerivAt s Ls (T0, p T0))
    (hv : HasFDerivAt v Lv (T0, p T0))
    (hgT : Lg (1, 0) = -s (T0, p T0)) (hgP : Lg (0, 1) = v (T0, p T0))
    (hMaxwell : Ls (0, 1) = -Lv (1, 0))
    (he : Lv (0, 1) ≠ 0) (hV : v (T0, p T0) ≠ 0)
    (hp : HasDerivAt p p' T0) (hiso : ∀ᶠ t in 𝓝 T0, v (t, p t) = v (T0, p T0)) :
    HasDerivAt (fun t => g (t, p t) + t * s (t, p t) - p t * v (t, p t))
      (T0 * Ls (1, 0) - T0 * v (T0, p T0) * (Lv (1, 0) / v (T0, p T0)) ^ 2
        / (-Lv (0, 1) / v (T0, p T0))) T0 := by
  have dg := comp_curve hg hp
  have ds := comp_curve hs hp
  have dv := comp_curve hv hp
  have hvz : Lv ((1 : ℝ), p') = 0 := isochore_slope v Lv T0 p p' hv hp hiso
  have hlin := hvz
  rw [clm_pair] at hlin
  have hp' : p' = -Lv (1, 0) / Lv (0, 1) := by
    rw [eq_div_iff he]; linarith
  have du := (dg.add ((hasDerivAt_idfun T0).mul ds)).sub (hp.mul dv)
  refine du.congr_deriv ?_
  rw [hvz, clm_pair Lg 1 p', clm_pair Ls 1 p', hgT, hgP, hMaxwell, hp']
  field_simp <;> ring

/-- When κ_T = 0 and β ≠ 0 no differentiable isochoric path exists: (∂u/∂T)_v is undefined, and the code returning c_p in this branch is only a convention. -/
theorem no_isochore (v : ℝ × ℝ → ℝ) (Lv : ℝ × ℝ →L[ℝ] ℝ) (T0 P0 : ℝ)
    (hv : HasFDerivAt v Lv (T0, P0)) (he : Lv (0, 1) = 0) (hc : Lv (1, 0) ≠ 0) :
    ¬ ∃ p : ℝ → ℝ, ∃ p' : ℝ, p T0 = P0 ∧ HasDerivAt p p' T0 ∧
      ∀ᶠ t in 𝓝 T0, v (t, p t) = v (T0, P0) := by
  rintro ⟨p, p', hp0, hp, hiso⟩
  subst hp0
  have hlin := isochore_slope v Lv T0 p p' hv hp hiso
  rw [clm_pair, he] at hlin
  exact hc (by linarith)

/-- Code-level bridge: if the code's β, κ_T equal the true partial derivatives of v (β v = ∂v/∂T, −κ_T v = ∂v/∂P), then
the formula of `cv_mole` (see `cvMole` in §2) equals du/dT along the isochoric path. -/
def cvMole (cp T V beta kappa : ℝ) : ℝ :=
  if kappa = 0 then cp else cp - T * V * beta * beta / kappa

theorem cv_bridge (g s v : ℝ × ℝ → ℝ) (Lg Ls Lv : ℝ × ℝ →L[ℝ] ℝ) (T0 : ℝ)
    (p : ℝ → ℝ) (p' : ℝ)
    (hg : HasFDerivAt g Lg (T0, p T0)) (hs : HasFDerivAt s Ls (T0, p T0))
    (hv : HasFDerivAt v Lv (T0, p T0))
    (hgT : Lg (1, 0) = -s (T0, p T0)) (hgP : Lg (0, 1) = v (T0, p T0))
    (hMaxwell : Ls (0, 1) = -Lv (1, 0))
    (hp : HasDerivAt p p' T0) (hiso : ∀ᶠ t in 𝓝 T0, v (t, p t) = v (T0, p T0))
    (beta kappa : ℝ) (hV : v (T0, p T0) ≠ 0) (hk : kappa ≠ 0)
    (hβ : HasDerivAt (fun t => v (t, p T0)) (beta * v (T0, p T0)) T0)
    (hκ : HasDerivAt (fun q => v (T0, q)) (-(kappa * v (T0, p T0))) (p T0)) :
    HasDerivAt (fun t => g (t, p t) + t * s (t, p t) - p t * v (t, p t))
      (cvMole (T0 * Ls (1, 0)) T0 (v (T0, p T0)) beta kappa) T0 := by
  have hc : Lv (1, 0) = beta * v (T0, p T0) := (partialT hv).unique hβ
  have he' : Lv (0, 1) = -(kappa * v (T0, p T0)) := (partialP hv).unique hκ
  have he : Lv (0, 1) ≠ 0 := by
    rw [he']; exact neg_ne_zero.mpr (mul_ne_zero hk hV)
  have hm := mayer_general g s v Lg Ls Lv T0 p p' hg hs hv hgT hgP hMaxwell he hV hp hiso
  refine hm.congr_deriv ?_
  unfold cvMole
  rw [if_neg hk, hc, he']
  field_simp <;> ring

/-! ## §2 MargulesVPSSTP -/

/-- Binary interaction parameter tables (m_pSpecies_A_ij … m_VSE_c_ij of MargulesVPSSTP.h; units J/kmol, J/kmol/K, m³/kmol, m³/kmol/K). -/
structure MargulesParams (K m : ℕ) where
  iA : Fin m → Fin K      -- m_pSpecies_A_ij
  iB : Fin m → Fin K      -- m_pSpecies_B_ij
  HEb : Fin m → ℝ         -- m_HE_b_ij
  HEc : Fin m → ℝ         -- m_HE_c_ij
  SEb : Fin m → ℝ         -- m_SE_b_ij
  SEc : Fin m → ℝ         -- m_SE_c_ij
  VHEb : Fin m → ℝ        -- m_VHE_b_ij
  VHEc : Fin m → ℝ        -- m_VHE_c_ij
  VSEb : Fin m → ℝ        -- m_VSE_b_ij
  VSEc : Fin m → ℝ        -- m_VSE_c_ij

variable {K m : ℕ}

/-- Accumulation structure shared by s_update_lnActCoeff / s_update_dlnActCoeff_dT / getPartialMolarVolumes
(g0, g1 are the respective coefficients):
```
MargulesVPSSTP.cpp:277-286 (d58a232)
        double XA = moleFractions_[iA];
        double XB = moleFractions_[iB];
        const double XAXB = XA * XB;
        const double g0g1XB = (g0 + g1 * XB);
        const double all = -1.0 * XAXB * g0g1XB - XAXB * XB * g1;
        for (size_t iK = 0; iK < m_kk; iK++) {
            lnActCoeff_Scaled_[iK] += all;
        }
        lnActCoeff_Scaled_[iA] += XB * g0g1XB;
        lnActCoeff_Scaled_[iB] += XA * g0g1XB + XAXB * g1;
MargulesVPSSTP.cpp:164-175
        double XA = moleFractions_[iA];
        double XB = moleFractions_[iB];
        double g0 = (m_VHE_b_ij[i] - T * m_VSE_b_ij[i]);
        double g1 = (m_VHE_c_ij[i] - T * m_VSE_c_ij[i]);
        const double temp1 = g0 + g1 * XB;
        const double all = -1.0*XA*XB*temp1 - XA*XB*XB*g1;
        for (size_t iK = 0; iK < m_kk; iK++) {
            vbar[iK] += all;
        }
        vbar[iA] += XB * temp1;
        vbar[iB] += XA * temp1 + XA*XB*g1;
``` -/
def pairSum (M : MargulesParams K m) (X : Fin K → ℝ) (g0 g1 : Fin m → ℝ) (k : Fin K) : ℝ :=
  ∑ i, ((-1 * (X (M.iA i) * X (M.iB i)) * (g0 i + g1 i * X (M.iB i))
          - X (M.iA i) * X (M.iB i) * X (M.iB i) * g1 i)
        + (if k = M.iA i then X (M.iB i) * (g0 i + g1 i * X (M.iB i)) else 0)
        + (if k = M.iB i then X (M.iA i) * (g0 i + g1 i * X (M.iB i))
              + X (M.iA i) * X (M.iB i) * g1 i else 0))

/-- The accumulation structure is linear in (g0, g1) (equivalent-rewriting lemma). -/
theorem pairSum_add_smul (M : MargulesParams K m) (X : Fin K → ℝ) (a a' b b' : Fin m → ℝ)
    (c : ℝ) (k : Fin K) :
    pairSum M X (fun i => a i + c * b i) (fun i => a' i + c * b' i) k
      = pairSum M X a a' k + c * pairSum M X b b' k := by
  unfold pairSum
  rw [Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  split_ifs <;> ring

theorem pairSum_smul (M : MargulesParams K m) (X : Fin K → ℝ) (a a' : Fin m → ℝ)
    (c : ℝ) (k : Fin K) :
    pairSum M X (fun i => c * a i) (fun i => c * a' i) k = c * pairSum M X a a' k := by
  unfold pairSum
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  split_ifs <;> ring

/-- Euler summation: Σ_k X_k · pairSum_k = Σ_i X_A X_B (g0 + g1 X_B) (requires Σ X = 1). -/
theorem pairSum_euler (M : MargulesParams K m) (X : Fin K → ℝ) (hX : ∑ k, X k = 1)
    (g0 g1 : Fin m → ℝ) :
    ∑ k, X k * pairSum M X g0 g1 k
      = ∑ i, X (M.iA i) * X (M.iB i) * (g0 i + g1 i * X (M.iB i)) := by
  unfold pairSum
  simp_rw [Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  first
  | ((simp only [mul_add, Finset.sum_add_distrib, mul_ite, mul_zero, Finset.sum_ite_eq',
      Finset.sum_ite_eq, Finset.mem_univ, ↓reduceIte, ← Finset.sum_mul, hX]) <;> ring)
  | ((simp [mul_add, Finset.sum_add_distrib, mul_ite, Finset.sum_ite_eq', ← Finset.sum_mul, hX])
      <;> ring)

/-- s_update_lnActCoeff (post-PR):
```
MargulesVPSSTP.cpp:265-276 (d58a232)
void MargulesVPSSTP::s_update_lnActCoeff() const
{
    double T = temperature();
    double deltaP_RT = (pressure() - OneAtm) / RT();
    ...
        double g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / RT()
                  + deltaP_RT * (m_VHE_b_ij[i] - T * m_VSE_b_ij[i]);
        double g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / RT()
                  + deltaP_RT * (m_VHE_c_ij[i] - T * m_VSE_c_ij[i]);
``` -/
def lnActCoeff (M : MargulesParams K m) (T P : ℝ) (X : Fin K → ℝ) (k : Fin K) : ℝ :=
  pairSum M X
    (fun i => (M.HEb i - T * M.SEb i) / RT T + (P - OneAtm) / RT T * (M.VHEb i - T * M.VSEb i))
    (fun i => (M.HEc i - T * M.SEc i) / RT T + (P - OneAtm) / RT T * (M.VHEc i - T * M.VSEc i))
    k

/-- s_update_lnActCoeff (pre-PR):
```
MargulesVPSSTP.cpp:237-238 (a6feb65, PR base)
        double g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / RT();
        double g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / RT();
``` -/
def lnActCoeffOld (M : MargulesParams K m) (T : ℝ) (X : Fin K → ℝ) (k : Fin K) : ℝ :=
  pairSum M X (fun i => (M.HEb i - T * M.SEb i) / RT T) (fun i => (M.HEc i - T * M.SEc i) / RT T) k

/-- getChemPotentials (same structure before and after the PR; mu0 = getStandardChemPotentials, an abstract parameter):
```
MargulesVPSSTP.cpp:40-52 (d58a232)
    getStandardChemPotentials(mu);
    s_update_lnActCoeff();
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(moleFractions_[k], SmallNumber);
        mu[k] += RT() * (log(xx) + lnActCoeff_Scaled_[k]);
    }
``` -/
def chemPotential (M : MargulesParams K m) (mu0 : Fin K → ℝ → ℝ → ℝ) (T P : ℝ)
    (X : Fin K → ℝ) (k : Fin K) : ℝ :=
  mu0 k T P + RT T * (Real.log (max (X k) SmallNumber) + lnActCoeff M T P X k)

/-- Chemical potential before the PR (MargulesVPSSTP.cpp:39-51 @ a6feb65, same structure, with the pre-PR version of ln γ). -/
def chemPotentialOld (M : MargulesParams K m) (mu0 : Fin K → ℝ → ℝ → ℝ) (T P : ℝ)
    (X : Fin K → ℝ) (k : Fin K) : ℝ :=
  mu0 k T P + RT T * (Real.log (max (X k) SmallNumber) + lnActCoeffOld M T X k)

/-- Excess part of getPartialMolarVolumes (MargulesVPSSTP.cpp:161-176 @ d58a232; unchanged by the PR,
identical to :126-141 of a6feb65). -/
def vbarExcess (M : MargulesParams K m) (T : ℝ) (X : Fin K → ℝ) (k : Fin K) : ℝ :=
  pairSum M X (fun i => M.VHEb i - T * M.VSEb i) (fun i => M.VHEc i - T * M.VSEc i) k

/-- getPartialMolarVolumes (MargulesVPSSTP.cpp:154-177; `getStandardVolumes(vbar)` gives v0, an abstract parameter). -/
def partialMolarVolume (M : MargulesParams K m) (v0 : Fin K → ℝ → ℝ → ℝ) (T P : ℝ)
    (X : Fin K → ℝ) (k : Fin K) : ℝ :=
  v0 k T P + vbarExcess M T X k

/-- RT·ln γ is affine in P: RT ln γ(P) = RT ln γ_old + (P − P_ref) v̄^E. -/
theorem RT_lnActCoeff_affine (M : MargulesParams K m) (T P : ℝ) (X : Fin K → ℝ) (k : Fin K)
    (hT : T ≠ 0) :
    RT T * lnActCoeff M T P X k
      = RT T * lnActCoeffOld M T X k + (P - OneAtm) * vbarExcess M T X k := by
  have hRT : RT T ≠ 0 := mul_ne_zero hT GasConstant_pos.ne'
  have h : lnActCoeff M T P X k
      = pairSum M X (fun i => (M.HEb i - T * M.SEb i) / RT T)
          (fun i => (M.HEc i - T * M.SEc i) / RT T) k
        + (P - OneAtm) / RT T
          * pairSum M X (fun i => M.VHEb i - T * M.VSEb i) (fun i => M.VHEc i - T * M.VSEc i) k :=
    pairSum_add_smul M X _ _ _ _ _ k
  have ho : lnActCoeffOld M T X k
      = pairSum M X (fun i => (M.HEb i - T * M.SEb i) / RT T)
          (fun i => (M.HEc i - T * M.SEc i) / RT T) k := rfl
  have hv : vbarExcess M T X k
      = pairSum M X (fun i => M.VHEb i - T * M.VSEb i) (fun i => M.VHEc i - T * M.VSEc i) k := rfl
  have hc : RT T * ((P - OneAtm) / RT T) = P - OneAtm := by
    rw [mul_comm (RT T)]; exact div_mul_self' hRT
  rw [h, ho, hv]
  linear_combination
    (pairSum M X (fun i => M.VHEb i - T * M.VSEb i) (fun i => M.VHEc i - T * M.VSEc i) k) * hc

/-- **Certificate (Maxwell, P)**: for all T ≠ 0, P, composition X (not required to be normalized) and species k,
∂μ_k/∂P |_{T,X} = v̄_k, with only the standard-state-layer ∂μ°_k/∂P = v°_k as premise. -/
theorem margules_dmu_dP_eq_vbar (M : MargulesParams K m) (mu0 v0 : Fin K → ℝ → ℝ → ℝ)
    (hstd : ∀ k T P, HasDerivAt (fun P => mu0 k T P) (v0 k T P) P)
    (T P : ℝ) (hT : T ≠ 0) (X : Fin K → ℝ) (k : Fin K) :
    HasDerivAt (fun P => chemPotential M mu0 T P X k) (partialMolarVolume M v0 T P X k) P := by
  have hfun : (fun q => chemPotential M mu0 T q X k) = fun q =>
      mu0 k T q + (RT T * Real.log (max (X k) SmallNumber) + RT T * lnActCoeffOld M T X k)
        + (q - OneAtm) * vbarExcess M T X k := by
    funext q
    unfold chemPotential
    rw [mul_add, RT_lnActCoeff_affine M T q X k hT]
    ring
  rw [hfun]
  have h2 : HasDerivAt (fun q => (q - OneAtm) * vbarExcess M T X k)
      (1 * vbarExcess M T X k) P :=
    ((hasDerivAt_idfun P).sub_const OneAtm).mul_const _
  have h3 := ((hstd k T P).add_const
    (RT T * Real.log (max (X k) SmallNumber) + RT T * lnActCoeffOld M T X k)).add h2
  refine h3.congr_deriv ?_
  unfold partialMolarVolume
  ring

/-- Pre-PR: ∂μ_k^old/∂P = v°_k (the excess part is independent of P). -/
theorem margules_old_dmu_dP (M : MargulesParams K m) (mu0 v0 : Fin K → ℝ → ℝ → ℝ)
    (hstd : ∀ k T P, HasDerivAt (fun P => mu0 k T P) (v0 k T P) P)
    (T P : ℝ) (X : Fin K → ℝ) (k : Fin K) :
    HasDerivAt (fun P => chemPotentialOld M mu0 T P X k) (v0 k T P) P :=
  (hstd k T P).add_const (RT T * (Real.log (max (X k) SmallNumber) + lnActCoeffOld M T X k))

/-- **Counterexample (pre-PR violates Maxwell)**: as long as v̄^E_k ≠ 0, the pre-PR μ_k does not satisfy ∂μ_k/∂P = v̄_k
(getPartialMolarVolumes already includes the excess volume while ln γ has no corresponding (P−P_ref) term). -/
theorem margules_old_maxwell_fails (M : MargulesParams K m) (mu0 v0 : Fin K → ℝ → ℝ → ℝ)
    (hstd : ∀ k T P, HasDerivAt (fun P => mu0 k T P) (v0 k T P) P)
    (T P : ℝ) (X : Fin K → ℝ) (k : Fin K) (hE : vbarExcess M T X k ≠ 0) :
    ¬ HasDerivAt (fun P => chemPotentialOld M mu0 T P X k) (partialMolarVolume M v0 T P X k) P := by
  intro h
  have h1 := h.unique (margules_old_dmu_dP M mu0 v0 hstd T P X k)
  unfold partialMolarVolume at h1
  exact hE (by linarith)

/-- Phase `ethanol-water-rough` of test/data/thermo-models.yaml (d58a232), converted to J/kmol, m³/kmol:
```
  - species: [ethanol, water]
    excess-enthalpy: [-2449 J/mol, -1613 J/mol]
    excess-entropy: [-14.45 J/mol/K, -11.64 J/mol/K]
    excess-volume-enthalpy: [-12.707 cm^3/mol, 0.476 cm^3/mol]
    excess-volume-entropy: [-0.030 cm^3/mol/K, 0.0 cm^3/mol/K]
``` -/
def ethanolWater : MargulesParams 2 1 where
  iA := fun _ => 0
  iB := fun _ => 1
  HEb := fun _ => -2449000
  HEc := fun _ => -1613000
  SEb := fun _ => -14450
  SEc := fun _ => -11640
  VHEb := fun _ => -0.012707
  VHEc := fun _ => 0.000476
  VSEb := fun _ => -0.00003
  VSEc := fun _ => 0

/-- Counterexample witness: the test phase ethanol-water-rough added by the PR has v̄^E_ethanol ≠ 0 at T = 298.15 K, X = (½, ½)
(value −9.40625e-4 m³/kmol). -/
theorem ethanolWater_vbarExcess_ne_zero :
    vbarExcess ethanolWater 298.15 ![1/2, 1/2] 0 ≠ 0 := by
  have h01 : (0 : Fin 2) ≠ 1 := by decide
  norm_num [vbarExcess, pairSum, ethanolWater, Fin.sum_univ_one, h01]

theorem margules_old_maxwell_fails_ethanolWater (mu0 v0 : Fin 2 → ℝ → ℝ → ℝ)
    (hstd : ∀ k T P, HasDerivAt (fun P => mu0 k T P) (v0 k T P) P) (P : ℝ) :
    ¬ HasDerivAt (fun P => chemPotentialOld ethanolWater mu0 298.15 P ![1/2, 1/2] 0)
        (partialMolarVolume ethanolWater v0 298.15 P ![1/2, 1/2] 0) P :=
  margules_old_maxwell_fails ethanolWater mu0 v0 hstd _ _ _ _ ethanolWater_vbarExcess_ne_zero

/-! ### Temperature derivatives -/

/-- Coefficients of s_update_dlnActCoeff_dT:
```
MargulesVPSSTP.cpp:292-303 (d58a232)
    double invT = 1.0 / temperature();
    double invRTT = 1.0 / GasConstant*invT*invT;
    double deltaP = pressure() - OneAtm;
    ...
        double g0 = -(m_HE_b_ij[i] + deltaP * m_VHE_b_ij[i]) * invRTT;
        double g1 = -(m_HE_c_ij[i] + deltaP * m_VHE_c_ij[i]) * invRTT;
``` -/
def gT0 (M : MargulesParams K m) (T P : ℝ) (i : Fin m) : ℝ :=
  -(M.HEb i + (P - OneAtm) * M.VHEb i) * (1 / GasConstant * (1 / T) * (1 / T))

def gT1 (M : MargulesParams K m) (T P : ℝ) (i : Fin m) : ℝ :=
  -(M.HEc i + (P - OneAtm) * M.VHEc i) * (1 / GasConstant * (1 / T) * (1 / T))

/-- dlnActCoeffdT_Scaled_ (MargulesVPSSTP.cpp:304-314, same accumulation structure as `pairSum`). -/
def dlnActCoeffdT (M : MargulesParams K m) (T P : ℝ) (X : Fin K → ℝ) (k : Fin K) : ℝ :=
  pairSum M X (gT0 M T P) (gT1 M T P) k

/-- d2lnActCoeffdT2_Scaled_:
```
MargulesVPSSTP.cpp:307-316 (d58a232)
        const double mult = 2.0 * invT;
        const double dT2all = mult * all;
        for (size_t iK = 0; iK < m_kk; iK++) {
            dlnActCoeffdT_Scaled_[iK] += all;
            d2lnActCoeffdT2_Scaled_[iK] -= dT2all;
        }
        ...
        d2lnActCoeffdT2_Scaled_[iA] -= mult * XB * g0g1XB;
        d2lnActCoeffdT2_Scaled_[iB] -= mult * (XA * g0g1XB + XAXB * g1);
``` -/
def d2lnActCoeffdT2 (M : MargulesParams K m) (T P : ℝ) (X : Fin K → ℝ) (k : Fin K) : ℝ :=
  ∑ i, (-(2 * (1 / T) * (-1 * (X (M.iA i) * X (M.iB i))
              * (gT0 M T P i + gT1 M T P i * X (M.iB i))
            - X (M.iA i) * X (M.iB i) * X (M.iB i) * gT1 M T P i))
        + (if k = M.iA i then
            -(2 * (1 / T) * X (M.iB i) * (gT0 M T P i + gT1 M T P i * X (M.iB i))) else 0)
        + (if k = M.iB i then
            -(2 * (1 / T) * (X (M.iA i) * (gT0 M T P i + gT1 M T P i * X (M.iB i))
              + X (M.iA i) * X (M.iB i) * gT1 M T P i)) else 0))

/-- Notation: ln γ = pairSum(dB) + T⁻¹ pairSum(dA). -/
def dA0 (M : MargulesParams K m) (P : ℝ) (i : Fin m) : ℝ :=
  (M.HEb i + (P - OneAtm) * M.VHEb i) / GasConstant
def dA1 (M : MargulesParams K m) (P : ℝ) (i : Fin m) : ℝ :=
  (M.HEc i + (P - OneAtm) * M.VHEc i) / GasConstant
def dB0 (M : MargulesParams K m) (P : ℝ) (i : Fin m) : ℝ :=
  -(M.SEb i + (P - OneAtm) * M.VSEb i) / GasConstant
def dB1 (M : MargulesParams K m) (P : ℝ) (i : Fin m) : ℝ :=
  -(M.SEc i + (P - OneAtm) * M.VSEc i) / GasConstant

theorem lnActCoeff_T_form (M : MargulesParams K m) (T P : ℝ) (X : Fin K → ℝ) (k : Fin K)
    (hT : T ≠ 0) :
    lnActCoeff M T P X k
      = pairSum M X (dB0 M P) (dB1 M P) k + T⁻¹ * pairSum M X (dA0 M P) (dA1 M P) k := by
  have hR := GasConstant_pos.ne'
  have h : lnActCoeff M T P X k
      = pairSum M X (fun i => dB0 M P i + T⁻¹ * dA0 M P i)
          (fun i => dB1 M P i + T⁻¹ * dA1 M P i) k := by
    unfold lnActCoeff
    congr 1 <;> funext i <;> simp only [dB0, dA0, dB1, dA1, RT] <;> field_simp <;> ring
  rw [h]
  exact pairSum_add_smul M X _ _ _ _ _ k

theorem dlnActCoeffdT_form (M : MargulesParams K m) (T P : ℝ) (X : Fin K → ℝ) (k : Fin K)
    (hT : T ≠ 0) :
    dlnActCoeffdT M T P X k = -(T ^ 2)⁻¹ * pairSum M X (dA0 M P) (dA1 M P) k := by
  have hR := GasConstant_pos.ne'
  have h : dlnActCoeffdT M T P X k
      = pairSum M X (fun i => -(T ^ 2)⁻¹ * dA0 M P i) (fun i => -(T ^ 2)⁻¹ * dA1 M P i) k := by
    unfold dlnActCoeffdT
    congr 1 <;> funext i <;> simp only [gT0, gT1, dA0, dA1] <;> field_simp <;> ring
  rw [h]
  exact pairSum_smul M X _ _ _ k

/-- Code structure: d²lnγ/dT² = −(2/T)·dlnγ/dT. -/
theorem d2lnActCoeffdT2_eq (M : MargulesParams K m) (T P : ℝ) (X : Fin K → ℝ) (k : Fin K) :
    d2lnActCoeffdT2 M T P X k = -(2 * (1 / T)) * dlnActCoeffdT M T P X k := by
  unfold d2lnActCoeffdT2 dlnActCoeffdT pairSum
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  split_ifs <;> ring

theorem eventually_pos_nhds {T : ℝ} (hT : 0 < T) : ∀ᶠ t in 𝓝 T, 0 < t := by
  first
  | exact eventually_gt_nhds hT
  | exact lt_mem_nhds hT
  | exact Ioi_mem_nhds hT

/-- **Certificate**: for all T > 0, P, X, k, the code's dlnγ/dT (including the deltaP·VHE term added by the PR) is the true derivative in T of ln γ (including the
(P−P_ref)V^E term added by the PR). -/
theorem margules_dlnActCoeff_dT_correct (M : MargulesParams K m) (P : ℝ) (X : Fin K → ℝ)
    (k : Fin K) (T : ℝ) (hT : 0 < T) :
    HasDerivAt (fun t => lnActCoeff M t P X k) (dlnActCoeffdT M T P X k) T := by
  have hev : (fun t => lnActCoeff M t P X k) =ᶠ[𝓝 T]
      (fun t => pairSum M X (dB0 M P) (dB1 M P) k + t⁻¹ * pairSum M X (dA0 M P) (dA1 M P) k) :=
    (eventually_pos_nhds hT).mono (fun t ht => lnActCoeff_T_form M t P X k ht.ne')
  have hd : HasDerivAt
      (fun t => pairSum M X (dB0 M P) (dB1 M P) k + t⁻¹ * pairSum M X (dA0 M P) (dA1 M P) k)
      (-(T ^ 2)⁻¹ * pairSum M X (dA0 M P) (dA1 M P) k) T :=
    ((hasDerivAt_inv hT.ne').mul_const _).const_add _
  refine (hd.congr_of_eventuallyEq hev).congr_deriv ?_
  rw [dlnActCoeffdT_form M T P X k hT.ne']

/-- **Certificate**: the code's d²lnγ/dT² is the true derivative in T of the code's dlnγ/dT (T > 0). -/
theorem margules_d2lnActCoeff_dT2_correct (M : MargulesParams K m) (P : ℝ) (X : Fin K → ℝ)
    (k : Fin K) (T : ℝ) (hT : 0 < T) :
    HasDerivAt (fun t => dlnActCoeffdT M t P X k) (d2lnActCoeffdT2 M T P X k) T := by
  have hev : (fun t => dlnActCoeffdT M t P X k) =ᶠ[𝓝 T]
      (fun t => -(t⁻¹ * t⁻¹) * pairSum M X (dA0 M P) (dA1 M P) k) :=
    (eventually_pos_nhds hT).mono (fun t ht => by
      try dsimp only
      rw [dlnActCoeffdT_form M t P X k ht.ne', pow_two, mul_inv])
  have hd : HasDerivAt (fun t => -(t⁻¹ * t⁻¹) * pairSum M X (dA0 M P) (dA1 M P) k)
      (-(-(T ^ 2)⁻¹ * T⁻¹ + T⁻¹ * -(T ^ 2)⁻¹) * pairSum M X (dA0 M P) (dA1 M P) k) T :=
    (((hasDerivAt_inv hT.ne').mul (hasDerivAt_inv hT.ne')).neg).mul_const _
  refine (hd.congr_of_eventuallyEq hev).congr_deriv ?_
  rw [d2lnActCoeffdT2_eq, dlnActCoeffdT_form M T P X k hT.ne']
  ring

/-- getPartialMolarEntropies (s0R = getEntropy_R, an abstract parameter):
```
MargulesVPSSTP.cpp:132-152 (d58a232)
    getEntropy_R(sbar);
    double T = temperature();
    s_update_lnActCoeff();
    s_update_dlnActCoeff_dT();
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(moleFractions_[k], SmallNumber);
        sbar[k] += - lnActCoeff_Scaled_[k] -log(xx) - T * dlnActCoeffdT_Scaled_[k];
    }
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] *= GasConstant;
    }
``` -/
def partialMolarEntropy (M : MargulesParams K m) (s0R : Fin K → ℝ → ℝ → ℝ) (T P : ℝ)
    (X : Fin K → ℝ) (k : Fin K) : ℝ :=
  (s0R k T P + (-(lnActCoeff M T P X k) - Real.log (max (X k) SmallNumber)
    - T * dlnActCoeffdT M T P X k)) * GasConstant

/-- getPartialMolarEnthalpies (h0RT = getEnthalpy_RT, an abstract parameter):
```
MargulesVPSSTP.cpp:93-110 (d58a232)
    getEnthalpy_RT(hbar);
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] *= RT();
    }
    s_update_lnActCoeff();
    s_update_dlnActCoeff_dT();
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] -= RT() * temperature() * dlnActCoeffdT_Scaled_[k];
    }
``` -/
def partialMolarEnthalpy (M : MargulesParams K m) (h0RT : Fin K → ℝ → ℝ → ℝ) (T P : ℝ)
    (X : Fin K → ℝ) (k : Fin K) : ℝ :=
  h0RT k T P * RT T - RT T * T * dlnActCoeffdT M T P X k

/-- **Certificate (Maxwell, T)**: for all T > 0, P, X, k, −∂μ_k/∂T = s̄_k (premise: standard-state ∂μ°/∂T = −s°). -/
theorem margules_dmu_dT_eq_neg_sbar (M : MargulesParams K m) (mu0 s0R : Fin K → ℝ → ℝ → ℝ)
    (hstdT : ∀ k T P, HasDerivAt (fun T => mu0 k T P) (-(s0R k T P * GasConstant)) T)
    (T P : ℝ) (hT : 0 < T) (X : Fin K → ℝ) (k : Fin K) :
    HasDerivAt (fun T => chemPotential M mu0 T P X k) (-(partialMolarEntropy M s0R T P X k)) T := by
  have hL := margules_dlnActCoeff_dT_correct M P X k T hT
  have hRT : HasDerivAt (fun t => RT t) GasConstant T :=
    ((hasDerivAt_idfun T).mul_const GasConstant).congr_deriv (one_mul _)
  have h := (hstdT k T P).add
    (hRT.mul ((hasDerivAt_const T (Real.log (max (X k) SmallNumber))).add hL))
  refine h.congr_deriv ?_
  unfold partialMolarEntropy RT
  simp only [Pi.add_apply]
  ring

/-- **Certificate (Gibbs–Helmholtz, algebraic)**: h̄_k = μ_k + T s̄_k (premise: standard-state h° = μ° + T s°). -/
theorem margules_hbar_eq_mu_plus_T_sbar (M : MargulesParams K m)
    (mu0 s0R h0RT : Fin K → ℝ → ℝ → ℝ)
    (hGH : ∀ k T P, h0RT k T P * RT T = mu0 k T P + T * (s0R k T P * GasConstant))
    (T P : ℝ) (X : Fin K → ℝ) (k : Fin K) :
    partialMolarEnthalpy M h0RT T P X k
      = chemPotential M mu0 T P X k + T * partialMolarEntropy M s0R T P X k := by
  unfold partialMolarEnthalpy chemPotential partialMolarEntropy
  rw [hGH]
  unfold RT
  ring

/-! ### β, κ_T, cv -/

/-- Molar volume (VPStandardStateTP.cpp:199-204; GibbsExcessVPSSTP::calcDensity was deleted by the PR):
```
void VPStandardStateTP::calcDensity()
{
    getPartialMolarVolumes(m_workS);
    double dd = meanMolecularWeight() / mean_X(m_workS);
    Phase::assignDensity(dd);
}
``` -/
def molarVolume (M : MargulesParams K m) (v0 : Fin K → ℝ → ℝ → ℝ) (T P : ℝ)
    (X : Fin K → ℝ) : ℝ :=
  ∑ k, X k * partialMolarVolume M v0 T P X k

/-- isothermalCompressibility (dVdP = PDSS::dVdP, an abstract parameter):
```
MargulesVPSSTP.cpp:64-73 (d58a232)
    double sum_xk_dVkdP = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        sum_xk_dVkdP += moleFractions_[k] * providePDSS(k)->dVdP();
    }
    return -sum_xk_dVkdP / molarVolume();
``` -/
def isothermalCompressibility (dVdP : Fin K → ℝ → ℝ → ℝ) (T P : ℝ) (X : Fin K → ℝ)
    (V : ℝ) : ℝ :=
  -(∑ k, X k * dVdP k T P) / V

/-- thermalExpansionCoeff (dVdT = PDSS::dVdT, an abstract parameter):
```
MargulesVPSSTP.cpp:75-91 (d58a232)
    double dVdT = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        dVdT += moleFractions_[k] * providePDSS(k)->dVdT();
    }
    for (size_t i = 0; i < numBinaryInteractions_; i++) {
        double XA = moleFractions_[m_pSpecies_A_ij[i]];
        double XB = moleFractions_[m_pSpecies_B_ij[i]];
        dVdT -= XA * XB * (m_VSE_b_ij[i] + m_VSE_c_ij[i] * XB);
    }
    return dVdT / molarVolume();
``` -/
def thermalExpansionCoeff (M : MargulesParams K m) (dVdT : Fin K → ℝ → ℝ → ℝ) (T P : ℝ)
    (X : Fin K → ℝ) (V : ℝ) : ℝ :=
  ((∑ k, X k * dVdT k T P)
    - ∑ i, X (M.iA i) * X (M.iB i) * (M.VSEb i + M.VSEc i * X (M.iB i))) / V

-- cv_mole (MargulesVPSSTP.cpp:54-62 (d58a232)) is `cvMole` of §1:
--     double kappa_T = isothermalCompressibility();
--     if (kappa_T == 0.0) { return cp_mole(); }
--     double beta = thermalExpansionCoeff();
--     return cp_mole() - temperature() * molarVolume() * beta * beta / kappa_T;

theorem vbarExcess_T_form (M : MargulesParams K m) (T : ℝ) (X : Fin K → ℝ) (k : Fin K) :
    vbarExcess M T X k
      = pairSum M X (fun i => M.VHEb i) (fun i => M.VHEc i) k
        + T * pairSum M X (fun i => -M.VSEb i) (fun i => -M.VSEc i) k := by
  have h : vbarExcess M T X k
      = pairSum M X (fun i => M.VHEb i + T * (-M.VSEb i)) (fun i => M.VHEc i + T * (-M.VSEc i)) k := by
    unfold vbarExcess
    congr 1 <;> funext i <;> ring
  rw [h]
  exact pairSum_add_smul M X _ _ _ _ T k

theorem vbarExcess_hasDerivAt_T (M : MargulesParams K m) (T : ℝ) (X : Fin K → ℝ) (k : Fin K) :
    HasDerivAt (fun t => vbarExcess M t X k)
      (pairSum M X (fun i => -M.VSEb i) (fun i => -M.VSEc i) k) T := by
  have hfun : (fun t => vbarExcess M t X k) = fun t =>
      pairSum M X (fun i => M.VHEb i) (fun i => M.VHEc i) k
        + t * pairSum M X (fun i => -M.VSEb i) (fun i => -M.VSEc i) k :=
    funext (fun t => vbarExcess_T_form M t X k)
  rw [hfun]
  exact (((hasDerivAt_idfun T).mul_const _).const_add _).congr_deriv (one_mul _)

/-- **Certificate**: the code's κ_T satisfies −κ_T·V = ∂V/∂P (V = molarVolume; premise: PDSS::dVdP is the true P derivative of v°, V ≠ 0).
That is, the comment that V^E is independent of P holds in this model. -/
theorem margules_kappa_correct (M : MargulesParams K m) (v0 dVdP : Fin K → ℝ → ℝ → ℝ)
    (T P : ℝ) (X : Fin K → ℝ)
    (hP : ∀ k, HasDerivAt (fun P => v0 k T P) (dVdP k T P) P)
    (hV : molarVolume M v0 T P X ≠ 0) :
    HasDerivAt (fun P => molarVolume M v0 T P X)
      (-(isothermalCompressibility dVdP T P X (molarVolume M v0 T P X)
        * molarVolume M v0 T P X)) P := by
  have h : HasDerivAt (fun P => ∑ k, X k * (v0 k T P + vbarExcess M T X k))
      (∑ k, X k * dVdP k T P) P := by
    first
    | (apply HasDerivAt.fun_sum; intro k _; exact ((hP k).add_const (vbarExcess M T X k)).const_mul (X k))
    | (apply HasDerivAt.sum; intro k _; exact ((hP k).add_const (vbarExcess M T X k)).const_mul (X k))
  refine h.congr_deriv ?_
  unfold isothermalCompressibility
  rw [div_mul_self' hV, neg_neg]

/-- **Certificate**: the code's β satisfies β·V = ∂V/∂T (including the excess-volume temperature term added by the PR; premise: PDSS::dVdT is the true T derivative of v°,
Σ X = 1, V ≠ 0). -/
theorem margules_beta_correct (M : MargulesParams K m) (v0 dVdT : Fin K → ℝ → ℝ → ℝ)
    (T P : ℝ) (X : Fin K → ℝ)
    (hT' : ∀ k, HasDerivAt (fun T => v0 k T P) (dVdT k T P) T)
    (hX : ∑ k, X k = 1) (hV : molarVolume M v0 T P X ≠ 0) :
    HasDerivAt (fun T => molarVolume M v0 T P X)
      (thermalExpansionCoeff M dVdT T P X (molarVolume M v0 T P X) * molarVolume M v0 T P X) T := by
  have h : HasDerivAt (fun T => ∑ k, X k * (v0 k T P + vbarExcess M T X k))
      (∑ k, X k * (dVdT k T P + pairSum M X (fun i => -M.VSEb i) (fun i => -M.VSEc i) k)) T := by
    first
    | (apply HasDerivAt.fun_sum; intro k _; exact ((hT' k).add (vbarExcess_hasDerivAt_T M T X k)).const_mul (X k))
    | (apply HasDerivAt.sum; intro k _; exact ((hT' k).add (vbarExcess_hasDerivAt_T M T X k)).const_mul (X k))
  refine h.congr_deriv ?_
  have he := pairSum_euler M X hX (fun i => -M.VSEb i) (fun i => -M.VSEc i)
  have hneg : ∑ i, X (M.iA i) * X (M.iB i) * (-M.VSEb i + -M.VSEc i * X (M.iB i))
      = -∑ i, X (M.iA i) * X (M.iB i) * (M.VSEb i + M.VSEc i * X (M.iB i)) := by
    rw [← Finset.sum_neg_distrib]
    exact Finset.sum_congr rfl (fun i _ => by ring)
  have hsplit : (∑ k, X k * (dVdT k T P + pairSum M X (fun i => -M.VSEb i) (fun i => -M.VSEc i) k))
      = (∑ k, X k * dVdT k T P)
        + ∑ k, X k * pairSum M X (fun i => -M.VSEb i) (fun i => -M.VSEc i) k := by
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl (fun k _ => by ring)
  rw [hsplit, he, hneg]
  unfold thermalExpansionCoeff
  rw [div_mul_self' hV]
  ring

/-- **Certificate (code-level Mayer, Margules)**: under the general-layer hypotheses, taking v = molarVolume,
`cv_mole` (κ_T ≠ 0 branch) equals (∂u/∂T)_v along the isochoric path, where c_p is taken as T(∂s/∂T)_P
(cp_mole itself is not transcribed and is absorbed into the general layer as a hypothesis). -/
theorem margules_cv_eq_isochoric_dudT (M : MargulesParams K m)
    (v0 dVdT dVdP : Fin K → ℝ → ℝ → ℝ) (X : Fin K → ℝ) (hX : ∑ k, X k = 1)
    (g s : ℝ × ℝ → ℝ) (Lg Ls Lv : ℝ × ℝ →L[ℝ] ℝ) (T0 : ℝ) (p : ℝ → ℝ) (p' : ℝ)
    (hg : HasFDerivAt g Lg (T0, p T0)) (hs : HasFDerivAt s Ls (T0, p T0))
    (hv : HasFDerivAt (fun q => molarVolume M v0 q.1 q.2 X) Lv (T0, p T0))
    (hgT : Lg (1, 0) = -s (T0, p T0)) (hgP : Lg (0, 1) = molarVolume M v0 T0 (p T0) X)
    (hMaxwell : Ls (0, 1) = -Lv (1, 0))
    (hT' : ∀ k, HasDerivAt (fun T => v0 k T (p T0)) (dVdT k T0 (p T0)) T0)
    (hP' : ∀ k, HasDerivAt (fun P => v0 k T0 P) (dVdP k T0 (p T0)) (p T0))
    (hV : molarVolume M v0 T0 (p T0) X ≠ 0)
    (hk : isothermalCompressibility dVdP T0 (p T0) X (molarVolume M v0 T0 (p T0) X) ≠ 0)
    (hp : HasDerivAt p p' T0)
    (hiso : ∀ᶠ t in 𝓝 T0, molarVolume M v0 t (p t) X = molarVolume M v0 T0 (p T0) X) :
    HasDerivAt (fun t => g (t, p t) + t * s (t, p t) - p t * molarVolume M v0 t (p t) X)
      (cvMole (T0 * Ls (1, 0)) T0 (molarVolume M v0 T0 (p T0) X)
        (thermalExpansionCoeff M dVdT T0 (p T0) X (molarVolume M v0 T0 (p T0) X))
        (isothermalCompressibility dVdP T0 (p T0) X (molarVolume M v0 T0 (p T0) X))) T0 :=
  cv_bridge g s (fun q => molarVolume M v0 q.1 q.2 X) Lg Ls Lv T0 p p' hg hs hv hgT hgP hMaxwell
    hp hiso _ _ hV hk
    (margules_beta_correct M v0 dVdT T0 (p T0) X hT' hX hV)
    (margules_kappa_correct M v0 dVdP T0 (p T0) X hP' hV)

/-- When PDSS::dVdP = 0 for all species (PDSS_ConstVol / PDSS_SSVol, such as the test phase ethanol-water-rough added by the PR),
the code's κ_T = 0 and cv_mole takes the `return cp_mole()` branch. -/
theorem margules_kappa_zero_branch (dVdP : Fin K → ℝ → ℝ → ℝ) (T P : ℝ) (X : Fin K → ℝ)
    (V cp beta : ℝ) (h0 : ∀ k, dVdP k T P = 0) :
    isothermalCompressibility dVdP T P X V = 0 ∧
      cvMole cp T V beta (isothermalCompressibility dVdP T P X V) = cp := by
  have hk : isothermalCompressibility dVdP T P X V = 0 := by
    unfold isothermalCompressibility; simp [h0]
  exact ⟨hk, by unfold cvMole; rw [if_pos hk]⟩

/-- **Counterexample (convention branch has no physical counterpart)**: when dVdP ≡ 0 while the code's β ≠ 0 (excess volume varies with T), no differentiable isochoric path exists,
(∂u/∂T)_v is undefined; the c_p returned by cv_mole is only a conventional value, and the Mayer test only checks the convention in this branch. -/
theorem margules_no_isochore_when_kappa_zero (M : MargulesParams K m)
    (v0 dVdT dVdP : Fin K → ℝ → ℝ → ℝ) (X : Fin K → ℝ) (hX : ∑ k, X k = 1)
    (T P : ℝ) (Lv : ℝ × ℝ →L[ℝ] ℝ)
    (hv : HasFDerivAt (fun q => molarVolume M v0 q.1 q.2 X) Lv (T, P))
    (hT' : ∀ k, HasDerivAt (fun T => v0 k T P) (dVdT k T P) T)
    (hP' : ∀ k, HasDerivAt (fun P => v0 k T P) (dVdP k T P) P)
    (h0 : ∀ k, dVdP k T P = 0) (hV : molarVolume M v0 T P X ≠ 0)
    (hβ : thermalExpansionCoeff M dVdT T P X (molarVolume M v0 T P X) ≠ 0) :
    ¬ ∃ p : ℝ → ℝ, ∃ p' : ℝ, p T = P ∧ HasDerivAt p p' T ∧
      ∀ᶠ t in 𝓝 T, molarVolume M v0 t (p t) X = molarVolume M v0 T P X := by
  have hc : Lv (1, 0)
      = thermalExpansionCoeff M dVdT T P X (molarVolume M v0 T P X) * molarVolume M v0 T P X :=
    (partialT hv).unique (margules_beta_correct M v0 dVdT T P X hT' hX hV)
  have he : Lv (0, 1) = 0 := by
    have h1 := (partialP hv).unique (margules_kappa_correct M v0 dVdP T P X hP' hV)
    rw [h1, (margules_kappa_zero_branch dVdP T P X (molarVolume M v0 T P X) 0 0 h0).1]
    simp
  exact no_isochore (fun q => molarVolume M v0 q.1 q.2 X) Lv T P hv he
    (by rw [hc]; exact mul_ne_zero hβ hV)

/-! ## §3 IdealSolnGasVPSS -/

/-- cv_mole (dVdT, dVdP are the values of each species' PDSS::dVdT/dVdP in the current state):
```
IdealSolnGasVPSS.cpp:62-77 (d58a232)
    double sum_xk_dVkdT = 0.0;
    double sum_xk_dVkdP = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        double Xk = moleFraction(k);
        sum_xk_dVkdT += Xk * providePDSS(k)->dVdT();
        sum_xk_dVkdP += Xk * providePDSS(k)->dVdP();
    }
    if (sum_xk_dVkdP == 0.0) {
        return cp_mole();
    }
    return cp_mole() - temperature() * sum_xk_dVkdT * sum_xk_dVkdT / (-sum_xk_dVkdP);
``` -/
def isgCvMole (cp T : ℝ) (X dVdT dVdP : Fin K → ℝ) : ℝ :=
  if (∑ k, X k * dVdP k) = 0 then cp
  else cp - T * (∑ k, X k * dVdT k) * (∑ k, X k * dVdT k) / (-(∑ k, X k * dVdP k))

/-- isothermalCompressibility (IdealSolnGasVPSS.cpp:79-87): `return -sum_xk_dVkdP / molarVolume();` -/
def isgKappa (X dVdP : Fin K → ℝ) (V : ℝ) : ℝ := -(∑ k, X k * dVdP k) / V

/-- thermalExpansionCoeff (IdealSolnGasVPSS.cpp:89-97): `return sum_xk_dVkdT / molarVolume();` -/
def isgBeta (X dVdT : Fin K → ℝ) (V : ℝ) : ℝ := (∑ k, X k * dVdT k) / V

/-- calcDensity (IdealSolnGasVPSS.cpp:106-111): `double v_mol = mean_X(getStandardVolumes());` -/
def isgMolarVolume (v0 : Fin K → ℝ → ℝ → ℝ) (T P : ℝ) (X : Fin K → ℝ) : ℝ :=
  ∑ k, X k * v0 k T P

-- The test cv_eq_cp_minus_T_V_beta2_over_kappa_T (test/thermo_consistency/consistency.cpp:423-441) compares:
--     if (kappa_T == 0.0) EXPECT_NEAR(cv, cp, atol);
--     else EXPECT_NEAR(cv, cp - T * V * beta * beta / kappa_T, ...);
-- i.e. cv_mole against `cvMole cp T V beta kappa_T` (β, κ_T, V taken from the same phase's thermalExpansionCoeff /
-- isothermalCompressibility / molarVolume).

/-- **Certificate (test tautology)**: for any V ≠ 0 (regardless of whether V equals Σ x_k v°_k, whether the PDSS derivatives are correct,
or whether cp is correct), the cv_mole of IdealSolnGasVPSS and the expected expression of the Mayer test are literally equal on both branches. -/
theorem isg_mayer_test_tautology (cp T : ℝ) (X dVdT dVdP : Fin K → ℝ) (V : ℝ) (hV : V ≠ 0) :
    isgCvMole cp T X dVdT dVdP = cvMole cp T V (isgBeta X dVdT V) (isgKappa X dVdP V) := by
  unfold isgCvMole cvMole isgBeta isgKappa
  by_cases h : (∑ k, X k * dVdP k) = 0
  · have h2 : -(∑ k, X k * dVdP k) / V = 0 := by rw [h]; simp
    rw [if_pos h, if_pos h2]
  · have h2 : -(∑ k, X k * dVdP k) / V ≠ 0 := div_ne_zero (neg_ne_zero.mpr h) hV
    rw [if_neg h, if_neg h2]
    field_simp <;> ring

/-- **Certificate (code-level Mayer, IdealSolnGasVPSS)**: if PDSS::dVdT/dVdP are the true derivatives of v°_k,
then cv_mole (κ_T ≠ 0 branch) equals the isochoric (∂u/∂T)_v with volume V = Σ x_k v°_k. This is something the test cannot provide:
it depends on the specific form of V and on the correctness of the PDSS derivatives. -/
theorem isg_cv_eq_isochoric_dudT (v0 dVdT dVdP : Fin K → ℝ → ℝ → ℝ) (X : Fin K → ℝ)
    (g s : ℝ × ℝ → ℝ) (Lg Ls Lv : ℝ × ℝ →L[ℝ] ℝ) (T0 : ℝ) (p : ℝ → ℝ) (p' : ℝ)
    (hg : HasFDerivAt g Lg (T0, p T0)) (hs : HasFDerivAt s Ls (T0, p T0))
    (hv : HasFDerivAt (fun q => isgMolarVolume v0 q.1 q.2 X) Lv (T0, p T0))
    (hgT : Lg (1, 0) = -s (T0, p T0)) (hgP : Lg (0, 1) = isgMolarVolume v0 T0 (p T0) X)
    (hMaxwell : Ls (0, 1) = -Lv (1, 0))
    (hT' : ∀ k, HasDerivAt (fun T => v0 k T (p T0)) (dVdT k T0 (p T0)) T0)
    (hP' : ∀ k, HasDerivAt (fun P => v0 k T0 P) (dVdP k T0 (p T0)) (p T0))
    (hV : isgMolarVolume v0 T0 (p T0) X ≠ 0)
    (hk : (∑ k, X k * dVdP k T0 (p T0)) ≠ 0)
    (hp : HasDerivAt p p' T0)
    (hiso : ∀ᶠ t in 𝓝 T0, isgMolarVolume v0 t (p t) X = isgMolarVolume v0 T0 (p T0) X) :
    HasDerivAt (fun t => g (t, p t) + t * s (t, p t) - p t * isgMolarVolume v0 t (p t) X)
      (isgCvMole (T0 * Ls (1, 0)) T0 X (fun k => dVdT k T0 (p T0))
        (fun k => dVdP k T0 (p T0))) T0 := by
  rw [isg_mayer_test_tautology _ _ X _ _ (isgMolarVolume v0 T0 (p T0) X) hV]
  have hβ : HasDerivAt (fun t => isgMolarVolume v0 t (p T0) X)
      (isgBeta X (fun k => dVdT k T0 (p T0)) (isgMolarVolume v0 T0 (p T0) X)
        * isgMolarVolume v0 T0 (p T0) X) T0 := by
    have h : HasDerivAt (fun t => ∑ k, X k * v0 k t (p T0)) (∑ k, X k * dVdT k T0 (p T0)) T0 := by
      first
      | (apply HasDerivAt.fun_sum; intro k _; exact (hT' k).const_mul (X k))
      | (apply HasDerivAt.sum; intro k _; exact (hT' k).const_mul (X k))
    unfold isgBeta
    rw [div_mul_self' hV]
    exact h
  have hκ : HasDerivAt (fun q => isgMolarVolume v0 T0 q X)
      (-(isgKappa X (fun k => dVdP k T0 (p T0)) (isgMolarVolume v0 T0 (p T0) X)
        * isgMolarVolume v0 T0 (p T0) X)) (p T0) := by
    have h : HasDerivAt (fun q => ∑ k, X k * v0 k T0 q) (∑ k, X k * dVdP k T0 (p T0)) (p T0) := by
      first
      | (apply HasDerivAt.fun_sum; intro k _; exact (hP' k).const_mul (X k))
      | (apply HasDerivAt.sum; intro k _; exact (hP' k).const_mul (X k))
    unfold isgKappa
    rw [div_mul_self' hV, neg_neg]
    exact h
  have hk' : isgKappa X (fun k => dVdP k T0 (p T0)) (isgMolarVolume v0 T0 (p T0) X) ≠ 0 :=
    div_ne_zero (neg_ne_zero.mpr hk) hV
  exact cv_bridge g s (fun q => isgMolarVolume v0 q.1 q.2 X) Lg Ls Lv T0 p p' hg hs hv hgT hgP
    hMaxwell hp hiso _ _ hV hk' hβ hκ

/-! ## §4 PDSS_SSVol (volumeModel_ == tpoly) -/

/-- PDSS_SSVol.cpp:110  `m_Vss = TCoeff_[0] + m_temp * (TCoeff_[1] + m_temp * (TCoeff_[2] + m_temp * TCoeff_[3]));` -/
def ssV (c : Fin 4 → ℝ) (T : ℝ) : ℝ := c 0 + T * (c 1 + T * (c 2 + T * c 3))

/-- PDSS_SSVol.cpp:112  `dVdT_ = TCoeff_[1] + 2.0 * m_temp * TCoeff_[2] + 3.0 * m_temp * m_temp * TCoeff_[3];` -/
def ssDVdT (c : Fin 4 → ℝ) (T : ℝ) : ℝ := c 1 + 2 * T * c 2 + 3 * T * T * c 3

/-- PDSS_SSVol.cpp:113  `d2VdT2_ = 2.0 * TCoeff_[2] + 6.0 * m_temp * TCoeff_[3];` -/
def ssD2VdT2 (c : Fin 4 → ℝ) (T : ℝ) : ℝ := 2 * c 2 + 6 * T * c 3

/-- m_sss_R of setTemperature (p0 = m_p0, s0R = m_s0_R(T)):
```
PDSS_SSVol.cpp:153-163 (d58a232)
    double deltaP = m_pres - m_p0;
    if (fabs(deltaP) < 1.0E-10) {
        ...
        m_sss_R = m_s0_R;
    } else {
        double sV_term = - deltaP / GasConstant * dVdT_;
        m_sss_R = m_s0_R + sV_term;
``` -/
def ssSssR (s0R : ℝ → ℝ) (c : Fin 4 → ℝ) (p0 P T : ℝ) : ℝ :=
  if |P - p0| < 1e-10 then s0R T else s0R T + -(P - p0) / GasConstant * ssDVdT c T

/-- m_cpss_R of setTemperature (post-PR):
```
PDSS_SSVol.cpp:158, 165 (d58a232)
        m_cpss_R = m_cp0_R;
        m_cpss_R = m_cp0_R - m_temp * deltaP * d2VdT2_ / GasConstant;
``` -/
def ssCpssR (cp0R : ℝ → ℝ) (c : Fin 4 → ℝ) (p0 P T : ℝ) : ℝ :=
  if |P - p0| < 1e-10 then cp0R T else cp0R T - T * (P - p0) * ssD2VdT2 c T / GasConstant

/-- m_cpss_R (pre-PR, PDSS_SSVol.cpp:155 @ a6feb65): `m_cpss_R = m_cp0_R - m_temp * deltaP * d2VdT2_;` -/
def ssCpssROld (cp0R : ℝ → ℝ) (c : Fin 4 → ℝ) (p0 P T : ℝ) : ℝ :=
  if |P - p0| < 1e-10 then cp0R T else cp0R T - T * (P - p0) * ssD2VdT2 c T

theorem ss_dVdT_correct (c : Fin 4 → ℝ) (T : ℝ) : HasDerivAt (ssV c) (ssDVdT c T) T := by
  have h3 : HasDerivAt (fun t => c 2 + t * c 3) (1 * c 3) T :=
    ((hasDerivAt_idfun T).mul_const (c 3)).const_add (c 2)
  have h2 : HasDerivAt (fun t => c 1 + t * (c 2 + t * c 3))
      (1 * (c 2 + T * c 3) + T * (1 * c 3)) T :=
    ((hasDerivAt_idfun T).mul h3).const_add (c 1)
  have h1 : HasDerivAt (fun t => c 0 + t * (c 1 + t * (c 2 + t * c 3)))
      (1 * (c 1 + T * (c 2 + T * c 3)) + T * (1 * (c 2 + T * c 3) + T * (1 * c 3))) T :=
    ((hasDerivAt_idfun T).mul h2).const_add (c 0)
  refine h1.congr_deriv ?_
  unfold ssDVdT
  ring

theorem ss_d2VdT2_correct (c : Fin 4 → ℝ) (T : ℝ) : HasDerivAt (ssDVdT c) (ssD2VdT2 c T) T := by
  have ha : HasDerivAt (fun t => 2 * t * c 2) (2 * 1 * c 2) T :=
    ((hasDerivAt_idfun T).const_mul 2).mul_const (c 2)
  have hb : HasDerivAt (fun t => 3 * t * t * c 3) ((3 * 1 * T + 3 * T * 1) * c 3) T :=
    (((hasDerivAt_idfun T).const_mul 3).mul (hasDerivAt_idfun T)).mul_const (c 3)
  have h := (ha.const_add (c 1)).add hb
  refine h.congr_deriv ?_
  unfold ssD2VdT2
  ring

/-- **Certificate**: for all T ≠ 0, P, p0 and coefficients c (valid in both |ΔP| branches), c_p,ss/R = T·∂(s_ss/R)/∂T,
premise: reference state cp0/R = T·d(s0/R)/dT. -/
theorem ss_cp_eq_T_dsdT (s0R cp0R : ℝ → ℝ) (c : Fin 4 → ℝ) (p0 P T : ℝ) (hT : T ≠ 0)
    (href : HasDerivAt s0R (cp0R T / T) T) :
    HasDerivAt (fun t => ssSssR s0R c p0 P t) (ssCpssR cp0R c p0 P T / T) T := by
  by_cases h : |P - p0| < 1e-10
  · have hf : (fun t => ssSssR s0R c p0 P t) = s0R := by
      funext t; simp [ssSssR, h]
    rw [hf]
    unfold ssCpssR
    rw [if_pos h]
    exact href
  · have hf : (fun t => ssSssR s0R c p0 P t)
        = fun t => s0R t + -(P - p0) / GasConstant * ssDVdT c t := by
      funext t; simp only [ssSssR, if_neg h]
    rw [hf]
    have hd := href.add ((ss_d2VdT2_correct c T).const_mul (-(P - p0) / GasConstant))
    refine hd.congr_deriv ?_
    unfold ssCpssR
    rw [if_neg h]
    have hR := GasConstant_pos.ne'
    field_simp <;> ring

/-- **Counterexample (pre-PR missing /R)**: when |ΔP| ≥ 1e-10 and ΔP·d²V/dT² ≠ 0, the pre-PR c_p,ss/R is not equal to T·∂(s_ss/R)/∂T. -/
theorem ss_old_cp_fails (s0R cp0R : ℝ → ℝ) (c : Fin 4 → ℝ) (p0 P T : ℝ) (hT : T ≠ 0)
    (href : HasDerivAt s0R (cp0R T / T) T) (h : ¬ |P - p0| < 1e-10)
    (hd : (P - p0) * ssD2VdT2 c T ≠ 0) :
    ¬ HasDerivAt (fun t => ssSssR s0R c p0 P t) (ssCpssROld cp0R c p0 P T / T) T := by
  intro hold
  have heq := hold.unique (ss_cp_eq_T_dsdT s0R cp0R c p0 P T hT href)
  unfold ssCpssROld ssCpssR at heq
  simp only [if_neg h] at heq
  rw [div_left_inj' hT] at heq
  have hR := GasConstant_pos.ne'
  have e : T * (P - p0) * ssD2VdT2 c T = T * (P - p0) * ssD2VdT2 c T / GasConstant := by
    linarith
  have e2 := (eq_div_iff hR).mp e
  have hy : T * ((P - p0) * ssD2VdT2 c T) * (GasConstant - 1) = 0 := by
    linear_combination e2
  rcases mul_eq_zero.mp hy with h1 | h1
  · rcases mul_eq_zero.mp h1 with h2 | h2
    · exact absurd h2 hT
    · exact hd h2
  · exact GasConstant_ne_one (by linarith)

end

end PhysAudit.Cantera2107
