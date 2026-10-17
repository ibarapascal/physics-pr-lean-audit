-- CoolProp3339: Lean certificate for CoolProp/CoolProp PR #3339 "fix(transport): rhosr-CS viscosity ignored x_crossover from the fluid file"
-- (merge commit 5418cec58c52183e3a010b2894e08b818ed4cc81, parent commit afcc79091533dfbd258cd59c065c0765747fd230).
-- Transcribes the crossover function and viscosity formula of both the pre-change and post-change versions of src/Backends/Helmholtz/TransportRoutines.cpp `viscosity_rhosr` (:1266-1295),
-- and proves that the two versions are formula-for-formula identical at x_crossover = 2, that the crossover function is continuous everywhere (logistic, no hard switch), and that outside the band the two versions differ only at the e^{-kδ} scale,
-- and gives two counterexamples: "identical formula-for-formula outside the band" fails over the reals, and x_crossover ≤ 0 breaks the dilute-gas limit.
import Mathlib

/-!
# CoolProp #3339 — rhosr-CS viscosity crossover

Transcription conventions
* Real-number semantics. The IEEE-754-level facts ("JSON integer 2 -> double 2.0 exactly") and `exp` overflow/underflow are out of scope for this file
  (the replay script `replay/scripts/CoolProp3339_replay.py` handles the floating-point comparison).
* The equation-of-state layer is not transcribed: `HEOS.rhomolar()`, `gas_constant()`, `tau()`, `dalphar_dTau()`, `alphar()`
  enter as arbitrary real inputs; for the dilute-gas part `viscosity_dilute_kinetic_theory` (:9-26) the transcription `viscosityDiluteKT` is given,
  but in the theorems `etaDil` is an arbitrary real parameter (it is independent of x_crossover).
* `c_liq` and `c_vap` are `std::vector<double>` in C++, and the code reads `[0..3]` unconditionally; here they are taken as `Fin 4 → ℝ`
  (i.e. the fluid file is assumed to supply at least 4 coefficients; extra entries are ignored).
-/

namespace PhysAudit.CoolProp3339

noncomputable section

/-! ## §1 Transcription -/

/-- ```
TransportRoutines.cpp:1280      const double crossover_sharpness = 100.0;
``` -/
def crossoverSharpness : ℝ := 100

/-- ```
TransportRoutines.cpp:1275      double x = HEOS.rhomolar() * HEOS.gas_constant() * (HEOS.tau() * HEOS.dalphar_dTau() - HEOS.alphar()) / data.rhosr_critical;
``` -/
def xRhosr (rhomolar R tau dalphar_dTau alphar rhosrCrit : ℝ) : ℝ :=
  rhomolar * R * (tau * dalphar_dTau - alphar) / rhosrCrit

/-- Crossover function of the merged version (`k` is `crossoverSharpness`, `xc` is `data.x_crossover`).
```
TransportRoutines.cpp:1281      double psi_liq = 1 / (1 + exp(-crossover_sharpness * (x - data.x_crossover)));
``` -/
def psiLiq (k x xc : ℝ) : ℝ := 1 / (1 + Real.exp (-k * (x - xc)))

/-- Crossover function before the change (parent commit afcc790):
```
TransportRoutines.cpp:1278 (afcc790)      double psi_liq = 1 / (1 + exp(-100.0 * (x - 2)));
``` -/
def psiLiqOld (x : ℝ) : ℝ := 1 / (1 + Real.exp (-100 * (x - 2)))

/-- ```
TransportRoutines.cpp:1285      double f_liq = cL[0] + x * (cL[1] + x * (cL[2] + x * (cL[3])));
TransportRoutines.cpp:1286      double f_vap = cV[0] + x * (cV[1] + x * (cV[2] + x * (cV[3])));
``` -/
def horner3 (c : Fin 4 → ℝ) (x : ℝ) : ℝ := c 0 + x * (c 1 + x * (c 2 + x * c 3))

/-- The exponent after crossover, `psi_liq * f_liq + (1 - psi_liq) * f_vap` (names the exp argument at :1289, structure unchanged). -/
def blendExponent (k xc fL fV x : ℝ) : ℝ := psiLiq k x xc * fL + (1 - psiLiq k x xc) * fV

/-- ```
TransportRoutines.cpp:1289      double etastar_ref = exp(psi_liq * f_liq + (1 - psi_liq) * f_vap);
``` -/
def etastarRef (k xc : ℝ) (cL cV : Fin 4 → ℝ) (x : ℝ) : ℝ :=
  Real.exp (psiLiq k x xc * horner3 cL x + (1 - psiLiq k x xc) * horner3 cV x)

/-- ```
TransportRoutines.cpp:1292      double etastar_fluid = 1 + data.C * (etastar_ref - 1);
``` -/
def etastarFluid (C k xc : ℝ) (cL cV : Fin 4 → ℝ) (x : ℝ) : ℝ :=
  1 + C * (etastarRef k xc cL cV x - 1)

/-- `viscosity_rhosr` of the merged version (`etaDil` = `viscosity_dilute_kinetic_theory(HEOS)`, :1272).
```
TransportRoutines.cpp:1294      return etastar_fluid * eta_dilute;
``` -/
def viscosityRhosr (etaDil C xc rhosrCrit : ℝ) (cL cV : Fin 4 → ℝ)
    (rhomolar R tau dalphar_dTau alphar : ℝ) : ℝ :=
  etastarFluid C crossoverSharpness xc cL cV (xRhosr rhomolar R tau dalphar_dTau alphar rhosrCrit)
    * etaDil

/-- `viscosity_rhosr` before the change (parent commit :1266-1292, line-for-line identical except for the crossover function). -/
def viscosityRhosrOld (etaDil C rhosrCrit : ℝ) (cL cV : Fin 4 → ℝ)
    (rhomolar R tau dalphar_dTau alphar : ℝ) : ℝ :=
  let x := xRhosr rhomolar R tau dalphar_dTau alphar rhosrCrit
  (1 + C * (Real.exp (psiLiqOld x * horner3 cL x + (1 - psiLiqOld x) * horner3 cV x) - 1)) * etaDil

/-- Dilute-gas part (independent of x_crossover, for comparison with the replay script).
```
TransportRoutines.cpp:11        CoolPropDbl Tstar = HEOS.T() / HEOS.components[0].transport.epsilon_over_k;
TransportRoutines.cpp:12        CoolPropDbl sigma_nm = HEOS.components[0].transport.sigma_eta * 1e9;
TransportRoutines.cpp:13        CoolPropDbl molar_mass_kgkmol = HEOS.molar_mass() * 1000;
TransportRoutines.cpp:18-19     CoolPropDbl OMEGA22 =
          1.16145 * pow(Tstar, static_cast<CoolPropDbl>(-0.14874)) + 0.52487 * exp(-0.77320 * Tstar) + 2.16178 * exp(-2.43787 * Tstar);
TransportRoutines.cpp:22        return 26.692e-9 * sqrt(molar_mass_kgkmol * HEOS.T()) / (pow(sigma_nm, 2) * OMEGA22);  // Pa-s
``` -/
def viscosityDiluteKT (T epsOverK sigmaEta molarMass : ℝ) : ℝ :=
  let Tstar := T / epsOverK
  let sigmaNm := sigmaEta * 1e9
  let mKmol := molarMass * 1000
  let omega22 := 1.16145 * Real.rpow Tstar (-0.14874) + 0.52487 * Real.exp (-0.77320 * Tstar)
    + 2.16178 * Real.exp (-2.43787 * Tstar)
  26.692e-9 * Real.sqrt (mKmol * T) / (sigmaNm ^ 2 * omega22)

/-! ## §2 (a) At x_crossover = 2 the two versions are formula-for-formula identical -/

/-- (certificate a, low content) For any state and any fluid coefficients, at x_crossover = 2 the merged version and the pre-change version are formula-for-formula identical (an identity at the definition level).
At the IEEE level one also needs "the JSON integer 2 is read as exactly the double 2.0"; that fact is not in Lean. -/
theorem viscosity_eq_old_of_xc_two (etaDil C rhosrCrit : ℝ) (cL cV : Fin 4 → ℝ)
    (rhomolar R tau dalphar_dTau alphar : ℝ) :
    viscosityRhosr etaDil C 2 rhosrCrit cL cV rhomolar R tau dalphar_dTau alphar
      = viscosityRhosrOld etaDil C rhosrCrit cL cV rhomolar R tau dalphar_dTau alphar := by
  simp only [viscosityRhosr, viscosityRhosrOld, etastarFluid, etastarRef, psiLiq, psiLiqOld,
    crossoverSharpness]

theorem psiLiq_eq_old (x : ℝ) : psiLiq crossoverSharpness x 2 = psiLiqOld x := by
  simp only [psiLiq, psiLiqOld, crossoverSharpness]

/-! ## §3 Properties of the crossover function -/

theorem psiLiq_pos (k x xc : ℝ) : 0 < psiLiq k x xc := by
  unfold psiLiq; positivity

theorem psiLiq_lt_one (k x xc : ℝ) : psiLiq k x xc < 1 := by
  unfold psiLiq
  rw [div_lt_one (by positivity)]
  linarith [Real.exp_pos (-k * (x - xc))]

/-- (certificate b) At the switching point ψ = 1/2 (the two branches have equal weight). -/
theorem psiLiq_at_crossover (k xc : ℝ) : psiLiq k xc xc = 1 / 2 := by
  norm_num [psiLiq]

/-- (certificate b) ψ is continuous in x (logistic switch, no hard switching point). -/
theorem continuous_psiLiq (k xc : ℝ) : Continuous (fun x => psiLiq k x xc) := by
  unfold psiLiq
  exact continuous_const.div (by fun_prop) (fun x => by positivity)

theorem continuous_horner3 (c : Fin 4 → ℝ) : Continuous (fun x => horner3 c x) := by
  unfold horner3; fun_prop

/-- (certificate b) The dimensionless viscosity η*(x) of the merged version is continuous in x everywhere (for any C, x_crossover, k, coefficients), including at the switching point. -/
theorem continuous_etastarFluid (C k xc : ℝ) (cL cV : Fin 4 → ℝ) :
    Continuous (fun x => etastarFluid C k xc cL cV x) := by
  have hψ := continuous_psiLiq k xc
  have hL := continuous_horner3 cL
  have hV := continuous_horner3 cV
  have hRef : Continuous (fun x => etastarRef k xc cL cV x) := by
    unfold etastarRef
    exact Real.continuous_exp.comp ((hψ.mul hL).add ((continuous_const.sub hψ).mul hV))
  unfold etastarFluid
  exact continuous_const.add (continuous_const.mul (hRef.sub continuous_const))

/-- (certificate b) At the switching point x = x_crossover, η* = 1 + C(exp((f_L + f_V)/2) − 1). -/
theorem etastarFluid_at_crossover (C k xc : ℝ) (cL cV : Fin 4 → ℝ) :
    etastarFluid C k xc cL cV xc
      = 1 + C * (Real.exp ((horner3 cL xc + horner3 cV xc) / 2) - 1) := by
  unfold etastarFluid etastarRef
  rw [psiLiq_at_crossover,
    show (1 : ℝ) / 2 * horner3 cL xc + (1 - 1 / 2) * horner3 cV xc
      = (horner3 cL xc + horner3 cV xc) / 2 by ring]

/-- For k > 0, ψ is strictly increasing in x. -/
theorem psiLiq_strictMono (k xc : ℝ) (hk : 0 < k) : StrictMono (fun x => psiLiq k x xc) := by
  intro x y hxy
  simp only [psiLiq]
  apply one_div_lt_one_div_of_lt
  · positivity
  · have : -k * (y - xc) < -k * (x - xc) := by nlinarith
    linarith [Real.exp_lt_exp.mpr this]

/-- For k > 0, ψ is strictly decreasing in x_crossover (a < b ⇒ ψ(x; b) < ψ(x; a)). -/
theorem psiLiq_strictAnti_xc (k : ℝ) (hk : 0 < k) (x : ℝ) {a b : ℝ} (hab : a < b) :
    psiLiq k x b < psiLiq k x a := by
  simp only [psiLiq]
  apply one_div_lt_one_div_of_lt
  · positivity
  · have : -k * (x - a) < -k * (x - b) := by nlinarith
    linarith [Real.exp_lt_exp.mpr this]

/-- (certificate b) The exponent is a convex combination of the two branches: it lies in [min f_L f_V, max f_L f_V]. -/
theorem blendExponent_mem_hull (k xc fL fV x : ℝ) :
    min fL fV ≤ blendExponent k xc fL fV x ∧ blendExponent k xc fL fV x ≤ max fL fV := by
  have hp := psiLiq_pos k x xc
  have hq := psiLiq_lt_one k x xc
  unfold blendExponent
  rcases le_total fL fV with h | h
  · rw [min_eq_left h, max_eq_right h]
    constructor
    · nlinarith [mul_nonneg (sub_nonneg.mpr h) (sub_nonneg.mpr hq.le)]
    · nlinarith [mul_nonneg (sub_nonneg.mpr h) hp.le]
  · rw [min_eq_right h, max_eq_left h]
    constructor
    · nlinarith [mul_nonneg (sub_nonneg.mpr h) hp.le]
    · nlinarith [mul_nonneg (sub_nonneg.mpr h) (sub_nonneg.mpr hq.le)]

theorem one_sub_psiLiq_le (k x xc : ℝ) : 1 - psiLiq k x xc ≤ Real.exp (-k * (x - xc)) := by
  have he := Real.exp_pos (-k * (x - xc))
  have h : 1 - Real.exp (-k * (x - xc)) ≤ psiLiq k x xc := by
    unfold psiLiq
    rw [le_div_iff₀ (by positivity)]
    nlinarith [sq_nonneg (Real.exp (-k * (x - xc)))]
  linarith

theorem psiLiq_le_exp (k x xc : ℝ) : psiLiq k x xc ≤ Real.exp (k * (x - xc)) := by
  have hE := Real.exp_pos (k * (x - xc))
  have hEe : Real.exp (k * (x - xc)) * Real.exp (-k * (x - xc)) = 1 := by
    rw [← Real.exp_add, show k * (x - xc) + -k * (x - xc) = 0 by ring, Real.exp_zero]
  unfold psiLiq
  rw [div_le_iff₀ (by positivity)]
  have : Real.exp (k * (x - xc)) * (1 + Real.exp (-k * (x - xc)))
      = Real.exp (k * (x - xc)) + 1 := by
    rw [mul_add, mul_one, hEe]
  linarith

/-! ## §4 (c) Outside the band: the two versions differ only at the e^{-kδ} scale; over the reals they are not formula-for-formula identical -/

/-- Above the switching band: x ≥ c + δ ⇒ 1 − ψ(x; c) ≤ e^{−kδ}. -/
theorem one_sub_psiLiq_le_above (k x c δ : ℝ) (hk : 0 ≤ k) (h : c + δ ≤ x) :
    1 - psiLiq k x c ≤ Real.exp (-k * δ) := by
  refine le_trans (one_sub_psiLiq_le k x c) (Real.exp_le_exp.mpr ?_)
  nlinarith [mul_nonneg hk (show (0 : ℝ) ≤ x - c - δ by linarith)]

/-- Below the switching band: x ≤ c − δ ⇒ ψ(x; c) ≤ e^{−kδ}. -/
theorem psiLiq_le_below (k x c δ : ℝ) (hk : 0 ≤ k) (h : x ≤ c - δ) :
    psiLiq k x c ≤ Real.exp (-k * δ) := by
  refine le_trans (psiLiq_le_exp k x c) (Real.exp_le_exp.mpr ?_)
  nlinarith [mul_nonneg hk (show (0 : ℝ) ≤ c - δ - x by linarith)]

/-- (certificate c) Two switching centers a, b, with x at least δ above both: |ψ_a − ψ_b| ≤ e^{−kδ}. -/
theorem psiLiq_close_above (k x a b δ : ℝ) (hk : 0 ≤ k) (ha : a + δ ≤ x) (hb : b + δ ≤ x) :
    |psiLiq k x a - psiLiq k x b| ≤ Real.exp (-k * δ) := by
  have h1 := one_sub_psiLiq_le_above k x a δ hk ha
  have h2 := one_sub_psiLiq_le_above k x b δ hk hb
  have h3 := psiLiq_lt_one k x a
  have h4 := psiLiq_lt_one k x b
  rw [abs_sub_le_iff]
  constructor <;> linarith

/-- (certificate c) x at least δ below both: |ψ_a − ψ_b| ≤ e^{−kδ}. -/
theorem psiLiq_close_below (k x a b δ : ℝ) (hk : 0 ≤ k) (ha : x ≤ a - δ) (hb : x ≤ b - δ) :
    |psiLiq k x a - psiLiq k x b| ≤ Real.exp (-k * δ) := by
  have h1 := psiLiq_le_below k x a δ hk ha
  have h2 := psiLiq_le_below k x b δ hk hb
  have h3 := psiLiq_pos k x a
  have h4 := psiLiq_pos k x b
  rw [abs_sub_le_iff]
  constructor <;> linarith

theorem blendExponent_sub (k a b fL fV x : ℝ) :
    blendExponent k a fL fV x - blendExponent k b fL fV x
      = (psiLiq k x a - psiLiq k x b) * (fL - fV) := by
  unfold blendExponent; ring

/-- (certificate c) Outside the band (x at distance ≥ δ from both centers and on the same side), the difference of the exponents of the two versions is ≤ e^{−kδ}·|f_L − f_V|.
With k = 100, δ = 0.2: e^{−20} ≈ 2.1e-9. -/
theorem blendExponent_close_outside (k a b fL fV x δ : ℝ) (hk : 0 ≤ k)
    (hside : (a + δ ≤ x ∧ b + δ ≤ x) ∨ (x ≤ a - δ ∧ x ≤ b - δ)) :
    |blendExponent k a fL fV x - blendExponent k b fL fV x|
      ≤ Real.exp (-k * δ) * |fL - fV| := by
  rw [blendExponent_sub, abs_mul]
  apply mul_le_mul_of_nonneg_right _ (abs_nonneg _)
  rcases hside with ⟨ha, hb⟩ | ⟨ha, hb⟩
  · exact psiLiq_close_above k x a b δ hk ha hb
  · exact psiLiq_close_below k x a b δ hk ha hb

/-- (counterexample c: "identical formula-for-formula outside the band" fails over the reals) As long as x_crossover differs (a < b), C ≠ 0,
f_L(x) ≠ f_V(x), the η* of the two versions differ at **every** x; outside the band they are only exponentially close. -/
theorem etastarFluid_ne_of_xc_lt (C a b x : ℝ) (cL cV : Fin 4 → ℝ) (hab : a < b) (hC : C ≠ 0)
    (hf : horner3 cL x ≠ horner3 cV x) :
    etastarFluid C crossoverSharpness a cL cV x ≠ etastarFluid C crossoverSharpness b cL cV x := by
  intro h
  unfold etastarFluid at h
  have hE : etastarRef crossoverSharpness a cL cV x = etastarRef crossoverSharpness b cL cV x := by
    have h' : C * (etastarRef crossoverSharpness a cL cV x
        - etastarRef crossoverSharpness b cL cV x) = 0 := by
      linear_combination h
    rcases mul_eq_zero.mp h' with h0 | h0
    · exact absurd h0 hC
    · linarith
  unfold etastarRef at hE
  have hg := Real.exp_injective hE
  have hpsi := psiLiq_strictAnti_xc crossoverSharpness (by norm_num [crossoverSharpness]) x hab
  have hprod : (psiLiq crossoverSharpness x a - psiLiq crossoverSharpness x b)
      * (horner3 cL x - horner3 cV x) = 0 := by
    linear_combination hg
  rcases mul_eq_zero.mp hprod with h0 | h0
  · linarith
  · exact hf (by linarith)

theorem viscosity_ne_of_xc_lt (etaDil C a b rhosrCrit : ℝ) (cL cV : Fin 4 → ℝ)
    (rhomolar R tau dalphar_dTau alphar : ℝ) (hab : a < b) (hC : C ≠ 0) (hη : etaDil ≠ 0)
    (hf : horner3 cL (xRhosr rhomolar R tau dalphar_dTau alphar rhosrCrit)
      ≠ horner3 cV (xRhosr rhomolar R tau dalphar_dTau alphar rhosrCrit)) :
    viscosityRhosr etaDil C a rhosrCrit cL cV rhomolar R tau dalphar_dTau alphar
      ≠ viscosityRhosr etaDil C b rhosrCrit cL cV rhomolar R tau dalphar_dTau alphar := by
  intro h
  unfold viscosityRhosr at h
  exact etastarFluid_ne_of_xc_lt C a b _ cL cV hab hC hf (mul_right_cancel₀ hη h)

/-! ## §5 Parameter semantics: the dilute-gas limit -/

/-- (certificate) At x = 0 (ρ → 0) with `c_vap[0] = 0`, the exponent of η* − 1 reduces to `ψ(0)·c_liq[0]`, and
|ψ(0)·c_liq[0]| ≤ e^{−k·x_c}·|c_liq[0]| (for x_c = 2, k = 100 this is ≤ e^{−200}|c_liq[0]|). -/
theorem dilute_exponent_bound (k xc : ℝ) (cL cV : Fin 4 → ℝ) (hcV : cV 0 = 0) :
    |psiLiq k 0 xc * horner3 cL 0 + (1 - psiLiq k 0 xc) * horner3 cV 0|
      ≤ Real.exp (k * (0 - xc)) * |cL 0| := by
  have hL : horner3 cL 0 = cL 0 := by simp [horner3]
  have hV : horner3 cV 0 = 0 := by simp [horner3, hcV]
  rw [hL, hV, mul_zero, add_zero, abs_mul, abs_of_pos (psiLiq_pos k 0 xc)]
  exact mul_le_mul_of_nonneg_right (psiLiq_le_exp k 0 xc) (abs_nonneg _)

/-- (counterexample: after the fix x_crossover becomes a live parameter, and the code does not validate its value) If the fluid file sets x_crossover to 0,
then at x = 0 we get ψ = 1/2, and the dilute-gas limit η* → 1 is broken (as long as C ≠ 0, c_liq[0] ≠ 0, c_vap[0] = 0). -/
theorem dilute_limit_broken_at_xc_zero (C : ℝ) (cL cV : Fin 4 → ℝ) (hC : C ≠ 0)
    (hcV : cV 0 = 0) (hcL : cL 0 ≠ 0) :
    etastarFluid C crossoverSharpness 0 cL cV 0 ≠ 1 := by
  intro h
  unfold etastarFluid etastarRef at h
  have hp : psiLiq crossoverSharpness 0 0 = 1 / 2 := psiLiq_at_crossover _ _
  have hL : horner3 cL 0 = cL 0 := by simp [horner3]
  have hV : horner3 cV 0 = 0 := by simp [horner3, hcV]
  rw [hp, hL, hV] at h
  have h2 : C * (Real.exp (1 / 2 * cL 0 + (1 - 1 / 2) * 0) - 1) = 0 := by linarith
  rcases mul_eq_zero.mp h2 with h0 | h0
  · exact hC h0
  · have h3 : Real.exp (1 / 2 * cL 0 + (1 - 1 / 2) * 0) = 1 := by linarith
    rw [Real.exp_eq_one_iff] at h3
    exact hcL (by linarith)

/-! ## section Arithmetic (side note: numerical bounds at the PR test point, not counted as a certificate) -/

section Arithmetic

/-- Side note: at the PR test point x = 9 (T = 300 K, ρ = 1e4 mol/m³): for the x_c = 12 clone ψ ≤ e^{−300},
and for the original fluid with x_c = 2, 1 − ψ ≤ e^{−700} (the two are almost entirely on the vapor branch and the liquid branch, respectively). -/
theorem pr_test_point_bounds :
    psiLiq crossoverSharpness 9 12 ≤ Real.exp (-300) ∧
      1 - psiLiq crossoverSharpness 9 2 ≤ Real.exp (-700) := by
  constructor
  · have h := psiLiq_le_exp crossoverSharpness 9 12
    have e : crossoverSharpness * (9 - 12) = -300 := by norm_num [crossoverSharpness]
    rwa [e] at h
  · have h := one_sub_psiLiq_le crossoverSharpness 9 2
    have e : -crossoverSharpness * (9 - 2) = -700 := by norm_num [crossoverSharpness]
    rwa [e] at h

end Arithmetic

end

end PhysAudit.CoolProp3339
