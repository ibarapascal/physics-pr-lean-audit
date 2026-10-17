-- PyBaMM5888: pybamm-team/PyBaMM PR #5888 "fix: convert Ai2020 electrolyte diffusivity from cm2/s to m2/s"
-- (merge commit 343ece2062dbb39789dbfa977e0f2daf96e9e370) Lean certificate.
-- Transcribes `electrolyte_diffusivity_Ai2020` in Ai2020.py:426-456 (post-merge and parent commit versions) and
-- `electrolyte_diffusivity_Valoen2005` in Xu2019.py:121-152; proves
--   (a) the two are identical for all (c_e, T);
--   (b) the post-merge value = unit conversion factor (1 cm)² = 1e-4 m² times Valøen–Reimers (2005) eq. 14 (cm²/s, mol/L);
--       the parent commit's value is exactly the cgs number itself, 1e4 times the post-merge value, and nowhere equal to the SI value;
--   (c) on the domain T > 229 + 5e-3·c_e, D is strictly increasing in T and strictly decreasing in c_e (including the derivative signs in HasDerivAt form),
--       the block c_e ∈ [0,4000], T ∈ [253,333] lies entirely inside the domain; across the singularity T = 229 + 5e-3·c_e both monotonicity properties fail (counterexamples).
import Mathlib

/-!
# PyBaMM #5888 — Ai2020 electrolyte diffusivity

Transcription conventions
* Python's `10 ** x` (float or pybamm `Power`, with positive constant base 10) is transcribed as `Real.rpow`: `(10 : ℝ) ^ x`.
* Literals are copied verbatim (`5e-3`, `0.22e-3`, `-4.43`, `1e-4` are all written as Lean scientific-notation literals, i.e. exact decimal rationals);
    floating-point rounding is outside the scope of this file and is compared by the replay script `replay/scripts/PyBaMM5888_replay.py`.
* Division by zero: Python floats raise ZeroDivisionError at `T = 229 + 5e-3·c_e`, while pybamm/numpy give inf (D becomes 0);
    Lean has `x / 0 = 0`. So at the singularity the Lean definition differs from the native one. The identities (a)(b) hold in Lean semantics for all reals,
    and away from the singularity they are identities of the true numerical values; all of (c) carries the domain assumption `InDomain` and does not touch the singularity.
* The Valøen–Reimers eq. 14 is not in the repository; it is written as `vrLog10D` according to the paper (c in mol/L, D in cm²/s), as the reference object for (b).
-/

namespace PyBaMM5888

noncomputable section

/-! ## §1 Transcription -/

/-- Exponent shared by the post-merge version and the parent commit (Valøen–Reimers fit of log10 D[cm²/s], c_e in mol/m³).
```
Ai2020.py:454      D_c_e = 1e-4 * 10 ** (-4.43 - 54 / (T - 229 - 5e-3 * c_e) - 0.22e-3 * c_e)
``` -/
def aiExponent (c_e T : ℝ) : ℝ := -4.43 - 54 / (T - 229 - 5e-3 * c_e) - 0.22e-3 * c_e

/-- Post-merge (343ece2062).
```
Ai2020.py:426  def electrolyte_diffusivity_Ai2020(c_e, T):
Ai2020.py:453      # The Valoen and Reimers (2005) fit gives cm2/s; convert to m2/s
Ai2020.py:454      D_c_e = 1e-4 * 10 ** (-4.43 - 54 / (T - 229 - 5e-3 * c_e) - 0.22e-3 * c_e)
Ai2020.py:456      return D_c_e
``` -/
def diffAi2020 (c_e T : ℝ) : ℝ := 1e-4 * (10 : ℝ) ^ aiExponent c_e T

/-- Parent commit (before the change).
```
Ai2020.py:451 (parent)      D_c_e = 10 ** (-4.43 - 54 / (T - 229 - 5e-3 * c_e) - 0.22e-3 * c_e)
``` -/
def diffAi2020Old (c_e T : ℝ) : ℝ := (10 : ℝ) ^ aiExponent c_e T

/-- ```
Xu2019.py:144      c_e = c_e / 1000
``` -/
def xuMolar (c_e : ℝ) : ℝ := c_e / 1000

/-- ```
Xu2019.py:146      T_g = 229 + 5 * c_e
``` -/
def xuTg (c_e : ℝ) : ℝ := 229 + 5 * xuMolar c_e

/-- ```
Xu2019.py:147      D_0 = -4.43 - 54 / (T - T_g)
``` -/
def xuD0 (c_e T : ℝ) : ℝ := -4.43 - 54 / (T - xuTg c_e)

/-- ```
Xu2019.py:148      D_1 = -0.22
``` -/
def xuD1 : ℝ := -0.22

/-- ```
Xu2019.py:121  def electrolyte_diffusivity_Valoen2005(c_e, T):
Xu2019.py:152      return (10 ** (D_0 + D_1 * c_e)) * 1e-4
```
(`c_e` at line 152 has already been reassigned to mol/L at line 144.) -/
def diffXu2019 (c_e T : ℝ) : ℝ := (10 : ℝ) ^ (xuD0 c_e T + xuD1 * xuMolar c_e) * 1e-4

/-- Valøen & Reimers, J. Electrochem. Soc. 152(5) A882 (2005), eq. 14 (formula from the literature, outside the repository):
log10 D[cm²/s] = −4.43 − 54/(T − 229 − 5.0 c) − 0.22 c, with c in mol/L and T in K. -/
def vrLog10D (c T : ℝ) : ℝ := -4.43 - 54 / (T - 229 - 5.0 * c) - 0.22 * c

/-- cgs value D[cm²/s] given by the literature formula. -/
def vrDcgs (c T : ℝ) : ℝ := (10 : ℝ) ^ vrLog10D c T

/-- 1 cm = 1/100 m. -/
def cmInM : ℝ := 1 / 100

/-- Area unit conversion: 1 cm² = (1/100 m)². -/
def cm2ToM2 : ℝ := cmInM ^ 2

/-- 1 m³ = 1000 L. -/
def litrePerM3 : ℝ := 1000

/-- Concentration unit conversion: c[mol/L] = c[mol/m³] / 1000. -/
def molPerM3ToMolPerL (c : ℝ) : ℝ := c / litrePerM3

/-- Domain: the denominator in the exponent is positive (away from the singularity T = 229 + 5e-3·c_e). -/
def InDomain (c_e T : ℝ) : Prop := 0 < T - 229 - 5e-3 * c_e

/-! ## §2 Literals and unit factors -/

lemma lit_5e3 : (5e-3 : ℝ) = 1 / 200 := by norm_num

lemma lit_022e3 : (0.22e-3 : ℝ) = 11 / 50000 := by norm_num

lemma lit_022 : (0.22 : ℝ) = 11 / 50 := by norm_num

lemma lit_50 : (5.0 : ℝ) = 5 := by norm_num

lemma lit_443 : (4.43 : ℝ) = 443 / 100 := by norm_num

/-- Unit factor: (1 cm)² = 1e-4 m² (i.e. the `1e-4` multiplied in the code). -/
theorem cm2ToM2_eq : cm2ToM2 = 1e-4 := by
  unfold cm2ToM2 cmInM
  norm_num

/-! ## §3 (a) Ai2020 and Xu2019 are identical -/

theorem aiExponent_eq_xu (c_e T : ℝ) :
    aiExponent c_e T = xuD0 c_e T + xuD1 * xuMolar c_e := by
  have hden : T - 229 - 5e-3 * c_e = T - xuTg c_e := by
    unfold xuTg xuMolar
    rw [lit_5e3]
    ring
  unfold aiExponent xuD0 xuD1 xuMolar
  rw [hden, lit_022e3, lit_022]
  ring

/-- Certificate (a): the post-merge Ai2020 electrolyte diffusivity equals the Valøen2005 implementation of Xu2019 for all (c_e, T). -/
theorem diffAi2020_eq_diffXu2019 (c_e T : ℝ) : diffAi2020 c_e T = diffXu2019 c_e T := by
  unfold diffAi2020 diffXu2019
  rw [aiExponent_eq_xu]
  ring

/-! ## §4 (b) Unit conversion -/

theorem aiExponent_eq_vr (c_e T : ℝ) :
    aiExponent c_e T = vrLog10D (molPerM3ToMolPerL c_e) T := by
  have hden : T - 229 - 5e-3 * c_e = T - 229 - 5.0 * molPerM3ToMolPerL c_e := by
    unfold molPerM3ToMolPerL litrePerM3
    rw [lit_5e3, lit_50]
    ring
  unfold aiExponent vrLog10D
  rw [hden, lit_022e3]
  unfold molPerM3ToMolPerL litrePerM3
  rw [lit_022]
  -- Convert the scientific-notation literals to rationals first, then hand over to ring (ring fails kernel checking when it handles `5.0` directly)
  simp only [lit_50, lit_443]
  ring

/-- Certificate (b): post-merge D[m²/s] = (1 cm/1 m)² × literature formula D[cm²/s] (c is first converted from mol/m³ to mol/L). -/
theorem diffAi2020_eq_SI_literature (c_e T : ℝ) :
    diffAi2020 c_e T = cm2ToM2 * vrDcgs (molPerM3ToMolPerL c_e) T := by
  unfold diffAi2020 vrDcgs
  rw [aiExponent_eq_vr, cm2ToM2_eq]

/-- The parent commit returns exactly the cgs value (cm²/s) but labels it m²/s. -/
theorem diffAi2020Old_eq_cgs (c_e T : ℝ) :
    diffAi2020Old c_e T = vrDcgs (molPerM3ToMolPerL c_e) T := by
  unfold diffAi2020Old vrDcgs
  rw [aiExponent_eq_vr]

/-- Parent commit = post-merge × 1e4 (on the whole domain). -/
theorem diffAi2020Old_eq_1e4_mul (c_e T : ℝ) : diffAi2020Old c_e T = 1e4 * diffAi2020 c_e T := by
  unfold diffAi2020Old diffAi2020
  rw [← mul_assoc]
  norm_num

theorem diffAi2020_pos (c_e T : ℝ) : 0 < diffAi2020 c_e T := by
  unfold diffAi2020
  exact mul_pos (by norm_num) (Real.rpow_pos_of_pos (by norm_num) _)

/-- The parent commit never equals the post-merge value at any (c_e, T). -/
theorem diffAi2020Old_ne (c_e T : ℝ) : diffAi2020Old c_e T ≠ diffAi2020 c_e T := by
  rw [diffAi2020Old_eq_1e4_mul]
  have h := diffAi2020_pos c_e T
  have h4 : (1e4 : ℝ) = 10000 := by norm_num
  rw [h4]
  intro heq
  linarith

/-- Counterexample (before the change): the parent commit never equals the literature formula in SI at any (c_e, T). -/
theorem diffAi2020Old_ne_SI_literature (c_e T : ℝ) :
    diffAi2020Old c_e T ≠ cm2ToM2 * vrDcgs (molPerM3ToMolPerL c_e) T := by
  rw [← diffAi2020_eq_SI_literature]
  exact diffAi2020Old_ne c_e T

/-! ## §5 (c) Monotonicity -/

lemma div54_lt {u₁ u₂ : ℝ} (h₁ : 0 < u₁) (h₁₂ : u₁ < u₂) : 54 / u₂ < 54 / u₁ := by
  have h₂ : 0 < u₂ := lt_trans h₁ h₁₂
  have e : 54 / u₁ - 54 / u₂ = 54 * (u₂ - u₁) / (u₁ * u₂) := by
    rw [div_sub_div _ _ (ne_of_gt h₁) (ne_of_gt h₂)]
    ring
  have hp : 0 < 54 * (u₂ - u₁) / (u₁ * u₂) := div_pos (by linarith) (mul_pos h₁ h₂)
  linarith

lemma aiExponent_lt_of_T_lt {c_e T₁ T₂ : ℝ} (h₁ : InDomain c_e T₁) (h₁₂ : T₁ < T₂) :
    aiExponent c_e T₁ < aiExponent c_e T₂ := by
  unfold InDomain at h₁
  have h := div54_lt h₁ (by linarith : T₁ - 229 - 5e-3 * c_e < T₂ - 229 - 5e-3 * c_e)
  unfold aiExponent
  linarith

lemma aiExponent_lt_of_c_lt {T c₁ c₂ : ℝ} (h₂ : InDomain c₂ T) (h₁₂ : c₁ < c₂) :
    aiExponent c₂ T < aiExponent c₁ T := by
  unfold InDomain at h₂
  have hu : T - 229 - 5e-3 * c₂ < T - 229 - 5e-3 * c₁ := by
    rw [lit_5e3]
    linarith
  have h := div54_lt h₂ hu
  unfold aiExponent
  rw [lit_022e3]
  linarith

/-- Certificate (c): within the domain, D is strictly increasing in T (T₂ is automatically in the domain). -/
theorem diffAi2020_lt_of_T_lt {c_e T₁ T₂ : ℝ} (h₁ : InDomain c_e T₁) (h₁₂ : T₁ < T₂) :
    diffAi2020 c_e T₁ < diffAi2020 c_e T₂ := by
  unfold diffAi2020
  exact mul_lt_mul_of_pos_left
    (Real.rpow_lt_rpow_of_exponent_lt (by norm_num) (aiExponent_lt_of_T_lt h₁ h₁₂)) (by norm_num)

/-- Certificate (c): within the domain, D is strictly decreasing in c_e (only the larger c₂ needs to be in the domain). -/
theorem diffAi2020_lt_of_c_lt {T c₁ c₂ : ℝ} (h₂ : InDomain c₂ T) (h₁₂ : c₁ < c₂) :
    diffAi2020 c₂ T < diffAi2020 c₁ T := by
  unfold diffAi2020
  exact mul_lt_mul_of_pos_left
    (Real.rpow_lt_rpow_of_exponent_lt (by norm_num) (aiExponent_lt_of_c_lt h₂ h₁₂)) (by norm_num)

theorem diffAi2020_strictMonoOn_T (c_e : ℝ) :
    StrictMonoOn (fun T => diffAi2020 c_e T) {T | InDomain c_e T} :=
  fun _ h₁ _ _ h₁₂ => diffAi2020_lt_of_T_lt h₁ h₁₂

theorem diffAi2020_strictAntiOn_c (T : ℝ) :
    StrictAntiOn (fun c => diffAi2020 c T) {c | InDomain c T} :=
  fun _ _ _ h₂ h₁₂ => diffAi2020_lt_of_c_lt h₂ h₁₂

/-- The common operating block c_e ≤ 4000 mol/m³, T ≥ 253 K lies entirely in the domain (denominator ≥ 4 K). -/
theorem box_subset_domain {c_e T : ℝ} (hc : c_e ≤ 4000) (hT : 253 ≤ T) : InDomain c_e T := by
  unfold InDomain
  rw [lit_5e3]
  linarith

/-- Certificate (c): for c_e ∈ [0,4000], D is strictly increasing on T ∈ [253,333]. -/
theorem diffAi2020_strictMonoOn_T_box {c_e : ℝ} (hc : c_e ∈ Set.Icc (0 : ℝ) 4000) :
    StrictMonoOn (fun T => diffAi2020 c_e T) (Set.Icc 253 333) :=
  fun _ h₁ _ _ h₁₂ => diffAi2020_lt_of_T_lt (box_subset_domain hc.2 h₁.1) h₁₂

/-- Certificate (c): for T ∈ [253,333], D is strictly decreasing on c_e ∈ [0,4000]. -/
theorem diffAi2020_strictAntiOn_c_box {T : ℝ} (hT : T ∈ Set.Icc (253 : ℝ) 333) :
    StrictAntiOn (fun c => diffAi2020 c T) (Set.Icc 0 4000) :=
  fun _ _ _ h₂ h₁₂ => diffAi2020_lt_of_c_lt (box_subset_domain h₂.2 hT.1) h₁₂

/-! ### Derivative form -/

lemma diffAi2020_eq_exp (c_e T : ℝ) :
    diffAi2020 c_e T = 1e-4 * Real.exp (Real.log 10 * aiExponent c_e T) := by
  unfold diffAi2020
  rw [Real.rpow_def_of_pos (by norm_num : (0 : ℝ) < 10)]

/-- ∂D/∂T = D · ln10 · 54/(T − 229 − 5e-3·c_e)² (within the domain). -/
theorem hasDerivAt_diffAi2020_T (c_e T : ℝ) (h : InDomain c_e T) :
    HasDerivAt (fun T' => diffAi2020 c_e T')
      (diffAi2020 c_e T * Real.log 10 * (54 / (T - 229 - 5e-3 * c_e) ^ 2)) T := by
  unfold InDomain at h
  have hne : T - 229 - 5e-3 * c_e ≠ 0 := ne_of_gt h
  have hu : HasDerivAt (fun T' : ℝ => T' - 229 - 5e-3 * c_e) 1 T :=
    ((hasDerivAt_id T).sub_const 229).sub_const (5e-3 * c_e)
  have hq : HasDerivAt (fun T' : ℝ => 54 / (T' - 229 - 5e-3 * c_e))
      ((0 * (T - 229 - 5e-3 * c_e) - 54 * 1) / (T - 229 - 5e-3 * c_e) ^ 2) T :=
    (hasDerivAt_const T (54 : ℝ)).div hu hne
  have hg : HasDerivAt (fun T' : ℝ => -4.43 - 54 / (T' - 229 - 5e-3 * c_e) - 0.22e-3 * c_e)
      (0 - (0 * (T - 229 - 5e-3 * c_e) - 54 * 1) / (T - 229 - 5e-3 * c_e) ^ 2) T :=
    ((hasDerivAt_const T (-4.43 : ℝ)).sub hq).sub_const (0.22e-3 * c_e)
  have hE : HasDerivAt
      (fun T' : ℝ => 1e-4 * Real.exp (Real.log 10 *
        (-4.43 - 54 / (T' - 229 - 5e-3 * c_e) - 0.22e-3 * c_e)))
      (1e-4 * (Real.exp (Real.log 10 * (-4.43 - 54 / (T - 229 - 5e-3 * c_e) - 0.22e-3 * c_e)) *
        (Real.log 10 *
          (0 - (0 * (T - 229 - 5e-3 * c_e) - 54 * 1) / (T - 229 - 5e-3 * c_e) ^ 2)))) T :=
    ((hg.const_mul (Real.log 10)).exp).const_mul 1e-4
  have hE' : HasDerivAt (fun T' => diffAi2020 c_e T')
      (1e-4 * (Real.exp (Real.log 10 * (-4.43 - 54 / (T - 229 - 5e-3 * c_e) - 0.22e-3 * c_e)) *
        (Real.log 10 *
          (0 - (0 * (T - 229 - 5e-3 * c_e) - 54 * 1) / (T - 229 - 5e-3 * c_e) ^ 2)))) T := by
    simp only [diffAi2020_eq_exp, aiExponent]
    exact hE
  convert hE' using 1
  rw [diffAi2020_eq_exp, aiExponent]
  ring

/-- ∂D/∂c_e = D · ln10 · (−54·5e-3/(T − 229 − 5e-3·c_e)² − 0.22e-3) (within the domain). -/
theorem hasDerivAt_diffAi2020_c (c_e T : ℝ) (h : InDomain c_e T) :
    HasDerivAt (fun c' => diffAi2020 c' T)
      (diffAi2020 c_e T * Real.log 10 *
        (-(54 * 5e-3) / (T - 229 - 5e-3 * c_e) ^ 2 - 0.22e-3)) c_e := by
  unfold InDomain at h
  have hne : T - 229 - 5e-3 * c_e ≠ 0 := ne_of_gt h
  have hu : HasDerivAt (fun c' : ℝ => T - 229 - 5e-3 * c') (0 - 5e-3 * 1) c_e :=
    ((hasDerivAt_const c_e T).sub_const 229).sub ((hasDerivAt_id c_e).const_mul 5e-3)
  have hq : HasDerivAt (fun c' : ℝ => 54 / (T - 229 - 5e-3 * c'))
      ((0 * (T - 229 - 5e-3 * c_e) - 54 * (0 - 5e-3 * 1)) / (T - 229 - 5e-3 * c_e) ^ 2) c_e :=
    (hasDerivAt_const c_e (54 : ℝ)).div hu hne
  have hg : HasDerivAt (fun c' : ℝ => -4.43 - 54 / (T - 229 - 5e-3 * c') - 0.22e-3 * c')
      (0 - (0 * (T - 229 - 5e-3 * c_e) - 54 * (0 - 5e-3 * 1)) / (T - 229 - 5e-3 * c_e) ^ 2
        - 0.22e-3 * 1) c_e :=
    ((hasDerivAt_const c_e (-4.43 : ℝ)).sub hq).sub ((hasDerivAt_id c_e).const_mul 0.22e-3)
  have hE : HasDerivAt
      (fun c' : ℝ => 1e-4 * Real.exp (Real.log 10 *
        (-4.43 - 54 / (T - 229 - 5e-3 * c') - 0.22e-3 * c')))
      (1e-4 * (Real.exp (Real.log 10 * (-4.43 - 54 / (T - 229 - 5e-3 * c_e) - 0.22e-3 * c_e)) *
        (Real.log 10 *
          (0 - (0 * (T - 229 - 5e-3 * c_e) - 54 * (0 - 5e-3 * 1)) / (T - 229 - 5e-3 * c_e) ^ 2
            - 0.22e-3 * 1)))) c_e :=
    ((hg.const_mul (Real.log 10)).exp).const_mul 1e-4
  have hE' : HasDerivAt (fun c' => diffAi2020 c' T)
      (1e-4 * (Real.exp (Real.log 10 * (-4.43 - 54 / (T - 229 - 5e-3 * c_e) - 0.22e-3 * c_e)) *
        (Real.log 10 *
          (0 - (0 * (T - 229 - 5e-3 * c_e) - 54 * (0 - 5e-3 * 1)) / (T - 229 - 5e-3 * c_e) ^ 2
            - 0.22e-3 * 1)))) c_e := by
    simp only [diffAi2020_eq_exp, aiExponent]
    exact hE
  convert hE' using 1
  rw [diffAi2020_eq_exp, aiExponent]
  ring

/-- Certificate (c): within the domain, ∂D/∂T > 0. -/
theorem deriv_diffAi2020_T_pos (c_e T : ℝ) (h : InDomain c_e T) :
    0 < deriv (fun T' => diffAi2020 c_e T') T := by
  rw [(hasDerivAt_diffAi2020_T c_e T h).deriv]
  have hl : 0 < Real.log 10 := Real.log_pos (by norm_num)
  have hu : 0 < T - 229 - 5e-3 * c_e := h
  exact mul_pos (mul_pos (diffAi2020_pos c_e T) hl) (div_pos (by norm_num) (pow_pos hu 2))

/-- Certificate (c): within the domain, ∂D/∂c_e < 0. -/
theorem deriv_diffAi2020_c_neg (c_e T : ℝ) (h : InDomain c_e T) :
    deriv (fun c' => diffAi2020 c' T) c_e < 0 := by
  rw [(hasDerivAt_diffAi2020_c c_e T h).deriv]
  have hl : 0 < Real.log 10 := Real.log_pos (by norm_num)
  have hu : 0 < T - 229 - 5e-3 * c_e := h
  have h1 : 0 < 54 * 5e-3 / (T - 229 - 5e-3 * c_e) ^ 2 := div_pos (by norm_num) (pow_pos hu 2)
  have h2 : (0 : ℝ) < 0.22e-3 := by norm_num
  have h3 : -(54 * 5e-3) / (T - 229 - 5e-3 * c_e) ^ 2 - 0.22e-3 < 0 := by
    rw [neg_div]
    linarith
  exact mul_neg_of_pos_of_neg (mul_pos (diffAi2020_pos c_e T) hl) h3

/-! ### Counterexamples: the domain assumption cannot be dropped -/

/-- Counterexample: at c_e = 0, D at 228 K (below the singularity at 229 K) is larger than D at 230 K; monotonicity in T fails across the singularity. -/
theorem not_monotone_T_across_pole :
    ¬ InDomain 0 228 ∧ (228 : ℝ) < 230 ∧ diffAi2020 0 230 < diffAi2020 0 228 := by
  refine ⟨?_, by norm_num, ?_⟩
  · unfold InDomain
    norm_num
  · unfold diffAi2020
    refine mul_lt_mul_of_pos_left (Real.rpow_lt_rpow_of_exponent_lt (by norm_num) ?_) (by norm_num)
    unfold aiExponent
    norm_num

/-- Counterexample: at T = 250 K, D at c_e = 5000 (past the singularity at 4200 mol/m³) is larger than D at c_e = 3000; monotonicity in c fails across the singularity. -/
theorem not_antitone_c_across_pole :
    InDomain 3000 250 ∧ ¬ InDomain 5000 250 ∧ (3000 : ℝ) < 5000 ∧
      diffAi2020 3000 250 < diffAi2020 5000 250 := by
  refine ⟨?_, ?_, by norm_num, ?_⟩
  · unfold InDomain
    norm_num
  · unfold InDomain
    norm_num
  · unfold diffAi2020
    refine mul_lt_mul_of_pos_left (Real.rpow_lt_rpow_of_exponent_lt (by norm_num) ?_) (by norm_num)
    unfold aiExponent
    norm_num

/-! ## Arithmetic -/

section Arithmetic
-- Side note: check of constants, not counted as a certificate

/-- Denominator at the corner (c_e, T) = (4000, 253) of the common operating block = 4 K. -/
theorem box_corner_margin : (253 : ℝ) - 229 - 5e-3 * 4000 = 4 := by norm_num

/-- Exponent at (1000 mol/m³, 298.15 K) = −4.65 − 54/64.15 (≈ −5.4918, D ≈ 3.22e-10 m²/s). -/
theorem exponent_at_1000_298 : aiExponent 1000 298.15 = -4.65 - 54 / 64.15 := by
  unfold aiExponent
  norm_num

end Arithmetic

end

end PyBaMM5888
