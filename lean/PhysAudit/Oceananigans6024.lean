-- Oceananigans6024: CliMA/Oceananigans.jl PR #6024 "Fix CATKE's bottom mixing length coefficient Cᵇ"
-- (merge commit 24de7026893af38df065920fbbf80f05024fd1b9) Lean certificate.
-- Transcribes CATKE's stable length scale, stability function, Richardson number, branch selection of the convective length, momentum mixing length, turbulent velocity, κu, dissipation length, and
-- dissipation rate; under the three model assumptions "neutral, local TKE equilibrium, constant stress", proves for arbitrary closure coefficients K_u = C·(Cˡᵒu³/CˡᵒD)^{1/4}·u★·d
-- (bottom boundary layer: C = Cᵇ, d = height above bottom); the numerical check of the default constants is placed in `section Arithmetic` (not counted as a certificate).
import Mathlib

/-!
# Oceananigans.jl #6024 — CATKE near-bottom mixing length coefficient

Transcription conventions
* Pointwise (single face / cell center) model: `stable_length_scaleᶜᶜᶠ` and `ᶜᶜᶜ` differ only in the evaluation location; here the same depth,
    height_above_bottom is taken (continuum-limit assumption). The face value `w★ = ℑz(√max(emin, e))` equals `√max(emin, e)` when e is locally uniform
    (assumption).
* `ℓᴺ = ifelse(N²⁺ == 0, Inf, …)`, `min(d, ℓᴺ)`: when N²⁺ = 0 the real-number transcription takes d directly (isomorphic to `min(d, Inf) = d`);
    the `isnan` branch does not occur over the reals.
* The formulas for the convective / entrainment lengths `ℓᶜ`, `ℓᵉ` are not transcribed and are taken as parameters; only the branches that decide whether to use them (convecting / entraining) are transcribed.
* `maximum_viscosity` defaults to `Inf` (catke_vertical_diffusivity.jl:51): transcribed as a real parameter `κmax`; the theorems only require K ≤ κmax.
* Three items are not transcribed and are taken as model assumptions: (M1) shear production = κu·S² (the code uses a discrete average of `shear_production`);
    (M2) local equilibrium: shear production = dissipation (at N² = 0 the buoyancy flux −κc N² = 0 is a fact of the code; vertical transport of TKE is dropped, steady state);
    (M3) constant-stress layer: κu·S = u★².
-/

namespace Oceananigans6024

noncomputable section

/-! ## §1 Pure-function transcription -/

/-- ```
TurbulenceClosures.jl:178 @inline clip(x) = max(zero(x), x)
``` -/
def clip (x : ℝ) : ℝ := max 0 x

/-- ```
catke_mixing_length.jl:208 @inline step(x, c, w) = max(zero(x), min(one(x), (x - c) / w))
``` -/
def stepFn (x c w : ℝ) : ℝ := max 0 (min 1 ((x - c) / w))

/-- ```
catke_mixing_length.jl:210 @inline function scale(Ri, σ⁻, σ⁰, σ∞, c, w)
catke_mixing_length.jl:211     σ⁺ = σ⁰ + (σ∞ - σ⁰) * step(Ri, c, w)
catke_mixing_length.jl:212     σ = σ⁻ * (Ri < 0) + σ⁺ * (Ri ≥ 0)
catke_mixing_length.jl:213     return σ
``` -/
def scaleFn (Ri σm σ0 σi c w : ℝ) : ℝ :=
  σm * (if Ri < 0 then 1 else 0) + (σ0 + (σi - σ0) * stepFn Ri c w) * (if Ri ≥ 0 then 1 else 0)

/-- ```
TKEBasedVerticalDiffusivities.jl:70 @inline function Riᶜᶜᶠ(i, j, k, grid, velocities, tracers, buoyancy)
TKEBasedVerticalDiffusivities.jl:73     S² = shearᶜᶜᶠ(i, j, k, grid, u, v)
TKEBasedVerticalDiffusivities.jl:74     N² = ∂z_b(i, j, k, grid, buoyancy, tracers)
TKEBasedVerticalDiffusivities.jl:75     Ri = N² / S²
TKEBasedVerticalDiffusivities.jl:76     return ifelse(N² == 0, zero(grid), Ri)
``` -/
def richardson (N2 S2 : ℝ) : ℝ := if N2 = 0 then 0 else N2 / S2

/-- ```
catke_mixing_length.jl:70 @inline function stable_length_scaleᶜᶜᶠ(i, j, k, grid, closure, e, velocities, tracers, buoyancy)
catke_mixing_length.jl:74     d_up   = Cˢ * depthᶜᶜᶠ(i, j, k, grid)
catke_mixing_length.jl:75     d_down = Cᵇ * height_above_bottomᶜᶜᶠ(i, j, k, grid)
catke_mixing_length.jl:76     d = min(d_up, d_down)
catke_mixing_length.jl:78     ℓᴺ = stratification_mixing_lengthᶜᶜᶠ(i, j, k, grid, closure, e, tracers, buoyancy)
catke_mixing_length.jl:80     ℓ = min(d, ℓᴺ)
catke_mixing_length.jl:81     ℓ = ifelse(isnan(ℓ), d, ℓ)
catke_mixing_length.jl:57     N²⁺ = clip(N²)
catke_mixing_length.jl:59     return ifelse(N²⁺ == 0, FT(Inf), w★ / sqrt(N²⁺))
``` -/
def stableLength (Cs Cb depth hab N2 wstar : ℝ) : ℝ :=
  if clip N2 = 0 then min (Cs * depth) (Cb * hab)
  else min (min (Cs * depth) (Cb * hab)) (wstar / Real.sqrt (clip N2))

/-- Branch selection of the convective length (ℓᶜ, ℓᵉ are parameters).
```
catke_mixing_length.jl:148     convecting = (Jᵇ > Jᵇᵋ) & (N² < 0)
catke_mixing_length.jl:149     entraining = (Jᵇ > Jᵇᵋ) & (N² > 0) & (N²_above < 0)
catke_mixing_length.jl:151     ℓ = ifelse(convecting, ℓᶜ,
catke_mixing_length.jl:152         ifelse(entraining, ℓᵉ, zero(grid)))
``` -/
def convSelect (Jb Jbe N2 N2above lconv lent : ℝ) : ℝ :=
  if Jb > Jbe ∧ N2 < 0 then lconv else if Jb > Jbe ∧ N2 > 0 ∧ N2above < 0 then lent else 0

/-- ```
catke_mixing_length.jl:241     ℓ★ = σ * stable_length_scaleᶜᶜᶠ(i, j, k, grid, closure, tracers.e, velocities, tracers, buoyancy)
catke_mixing_length.jl:245     ℓu = max(ℓ★, ℓʰ)
catke_mixing_length.jl:247     H = static_column_depthᶜᶜᵃ(i, j, grid)
catke_mixing_length.jl:248     return min(H, ℓu)
``` -/
def momentumMixingLength (H σ ℓstable ℓh : ℝ) : ℝ := min H (max (σ * ℓstable) ℓh)

/-- ```
TKEBasedVerticalDiffusivities.jl:139 @inline function turbulent_velocityᶜᶜᶜ(i, j, k, grid, closure, e)
TKEBasedVerticalDiffusivities.jl:140     eᵢ = @inbounds e[i, j, k]
TKEBasedVerticalDiffusivities.jl:141     eᵐⁱⁿ = closure.minimum_tke
TKEBasedVerticalDiffusivities.jl:142     return sqrt(max(eᵐⁱⁿ, eᵢ))
``` -/
def turbulentVelocity (emin e : ℝ) : ℝ := Real.sqrt (max emin e)

/-- ```
catke_vertical_diffusivity.jl:340 @inline function κuᶜᶜᶠ(i, j, k, grid, closure, velocities, tracers, buoyancy, surface_buoyancy_flux)
catke_vertical_diffusivity.jl:343     κu = ℓu * w★
catke_vertical_diffusivity.jl:345     κu★ = min(κu, κu_max)
``` -/
def kappaU (ℓu wstar κmax : ℝ) : ℝ := min (ℓu * wstar) κmax

/-- ```
catke_equation.jl:38     σᴰ = stability_functionᶜᶜᶜ(i, j, k, grid, closure, Cᵘⁿ, Cˡᵒ, Cʰⁱ, velocities, tracers, buoyancy)
catke_equation.jl:39     ℓ★ = stable_length_scaleᶜᶜᶜ(i, j, k, grid, closure, tracers.e, velocities, tracers, buoyancy)
catke_equation.jl:40     ℓ★ = ℓ★ / σᴰ
catke_equation.jl:45     ℓᴰ = max(ℓ★, ℓʰ)
catke_equation.jl:47     H = static_column_depthᶜᶜᵃ(i, j, grid)
catke_equation.jl:48     return min(H, ℓᴰ)
``` -/
def dissipationLength (H σD ℓstable ℓh : ℝ) : ℝ := min H (max (ℓstable / σD) ℓh)

/-- ```
catke_equation.jl:67     ω_numerical = 1 / closure.negative_tke_damping_time_scale
catke_equation.jl:68     ω_physical = sqrt(abs(eᵢ)) / ℓᴰ
catke_equation.jl:70     return ifelse(eᵢ < 0, ω_numerical, ω_physical)
``` -/
def dissipationRate (τ e ℓD : ℝ) : ℝ := if e < 0 then 1 / τ else Real.sqrt |e| / ℓD

/-- ```
catke_equation.jl:75     eᵢ = @inbounds tracers.e[i, j, k]
catke_equation.jl:76     ω = dissipation_rate(i, j, k, grid, closure, velocities, tracers, args...)
catke_equation.jl:77     return ω * eᵢ
``` -/
def dissipation (τ e ℓD : ℝ) : ℝ := dissipationRate τ e ℓD * e

/-- Pointwise κu: `κuᶜᶜᶠ ∘ momentum_mixing_lengthᶜᶜᶠ` (stability_function uses (Cᵘⁿu, Cˡᵒu, Cʰⁱu), catke_mixing_length.jl:230-249). -/
def kappaUPoint (Cs Cb Cuun Cu Cuhi CRi0 CRid H depth hab N2 N2above S2 emin e κmax Jb Jbe lconv lent : ℝ) :
    ℝ :=
  kappaU
    (momentumMixingLength H (scaleFn (richardson N2 S2) Cuun Cu Cuhi CRi0 CRid)
      (stableLength Cs Cb depth hab N2 (turbulentVelocity emin e))
      (convSelect Jb Jbe N2 N2above lconv lent))
    (turbulentVelocity emin e) κmax

/-- Pointwise dissipation: `dissipation ∘ dissipation_rate ∘ dissipation_length_scaleᶜᶜᶜ` (stability_function uses
(CᵘⁿD, CˡᵒD, CʰⁱD), catke_equation.jl:24-49). -/
def dissipationPoint (Cs Cb CDun CD CDhi CRi0 CRid H depth hab N2 N2above S2 emin e τ Jb Jbe lconvD lentD : ℝ) :
    ℝ :=
  dissipation τ e
    (dissipationLength H (scaleFn (richardson N2 S2) CDun CD CDhi CRi0 CRid)
      (stableLength Cs Cb depth hab N2 (turbulentVelocity emin e))
      (convSelect Jb Jbe N2 N2above lconvD lentD))

/-! ## §2 Simplification of each function under neutral stratification (facts of the code, not assumptions) -/

lemma richardson_neutral (S2 : ℝ) : richardson 0 S2 = 0 := by
  simp [richardson]

lemma stepFn_zero (c w : ℝ) (hc : 0 ≤ c) (hw : 0 < w) : stepFn 0 c w = 0 := by
  have h := div_nonneg hc hw.le
  have he : (0 - c) / w = -(c / w) := by ring
  unfold stepFn
  have h1 : -(c / w) ≤ 1 := by linarith
  have h2 : -(c / w) ≤ 0 := by linarith
  rw [he, min_eq_right h1, max_eq_left h2]

/-- When Ri = 0 (and CRi⁰ ≥ 0, CRiᵟ > 0) the stability function takes the low-Ri value σ⁰ (= Cˡᵒ). -/
lemma scaleFn_neutral (σm σ0 σi c w : ℝ) (hc : 0 ≤ c) (hw : 0 < w) :
    scaleFn 0 σm σ0 σi c w = σ0 := by
  unfold scaleFn
  rw [stepFn_zero c w hc hw]
  simp

lemma stableLength_neutral (Cs Cb depth hab wstar : ℝ) :
    stableLength Cs Cb depth hab 0 wstar = min (Cs * depth) (Cb * hab) := by
  simp [stableLength, clip]

lemma convSelect_neutral (Jb Jbe N2above lconv lent : ℝ) : convSelect Jb Jbe 0 N2above lconv lent = 0 := by
  simp [convSelect]

/-! ## §3 Algebra of local equilibrium -/

lemma eq_of_pow_four_eq {a b : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) (h : a ^ 4 = b ^ 4) : a = b := by
  have hfac : (a - b) * ((a + b) * (a ^ 2 + b ^ 2)) = 0 := by linear_combination h
  rcases mul_eq_zero.mp hfac with h1 | h1
  · linarith
  · rcases mul_eq_zero.mp h1 with h2 | h2
    · linarith
    · have ha2 : a ^ 2 = 0 := by nlinarith [sq_nonneg a, sq_nonneg b]
      have hb2 : b ^ 2 = 0 := by nlinarith [sq_nonneg a, sq_nonneg b]
      have ha0 : a = 0 := by nlinarith [sq_nonneg a]
      have hb0 : b = 0 := by nlinarith [sq_nonneg b]
      rw [ha0, hb0]

lemma rpow_quarter_pow_four {x : ℝ} (hx : 0 ≤ x) : (x ^ ((1 : ℝ) / 4)) ^ 4 = x := by
  rw [← Real.rpow_natCast, ← Real.rpow_mul hx]
  norm_num

/-- Local equilibrium `C_u d r S² = (r / (d / C_D)) e` (e = r²) and constant stress `C_u d r S = u★²` give the TKE level
`C_u C_D r⁴ = u★⁴`, i.e. `e = u★² / √(C_u C_D)`. -/
lemma neutral_balance_tke (Cu CD d r S ustar : ℝ) (hCD : 0 < CD) (hd : 0 < d) (hr : 0 < r)
    (hP : Cu * d * r * S ^ 2 = r / (d / CD) * r ^ 2) (hτ : Cu * d * r * S = ustar ^ 2) :
    Cu * CD * r ^ 4 = ustar ^ 4 := by
  have hd' : d ≠ 0 := hd.ne'
  have hr' : r ≠ 0 := hr.ne'
  have hrw : r / (d / CD) * r ^ 2 = CD * r ^ 3 / d := by
    rw [div_div_eq_mul_div]
    ring
  rw [hrw, eq_div_iff hd'] at hP
  have h1 : Cu * d ^ 2 * S ^ 2 = CD * r ^ 2 := by
    apply mul_left_cancel₀ hr'
    linear_combination hP
  linear_combination (-Cu * r ^ 2) * h1 + (Cu * d * r * S + ustar ^ 2) * hτ

/-- The same two equations give `K⁴ = (C_u³/C_D) d⁴ u★⁴`, with K = C_u d r. -/
lemma neutral_balance_K4 (Cu CD d r S ustar : ℝ) (hCD : 0 < CD) (hd : 0 < d) (hr : 0 < r)
    (hP : Cu * d * r * S ^ 2 = r / (d / CD) * r ^ 2) (hτ : Cu * d * r * S = ustar ^ 2) :
    (Cu * d * r) ^ 4 = Cu ^ 3 / CD * d ^ 4 * ustar ^ 4 := by
  have h2 := neutral_balance_tke Cu CD d r S ustar hCD hd hr hP hτ
  have hCD' : CD ≠ 0 := hCD.ne'
  rw [← h2, div_mul_eq_mul_div, div_mul_eq_mul_div, eq_div_iff hCD']
  ring

/-! ## §4 Certificate: neutral wall layer K_u = C (Cˡᵒu³/CˡᵒD)^{1/4} u★ d -/

/-- Certificate (general form): for arbitrary closure coefficients, at points with neutral stratification (N² = 0), no convection, length not truncated by the column depth H and κmax, e ≥ emin
and e > 0, if (M1)-(M3) hold then `K_u = (Cˡᵒu³/CˡᵒD)^{1/4} · d · u★`,
where `d = min(Cˢ·depth, Cᵇ·height_above_bottom)` is the wall-distance length in the code. -/
theorem neutral_wall_law
    (Cs Cb Cuun Cu Cuhi CDun CD CDhi CRi0 CRid H depth hab N2above S emin e κmax τ Jb Jbe lconv lent lconvD lentD
      ustar : ℝ)
    (hCu : 0 < Cu) (hCD : 0 < CD) (hCRi0 : 0 ≤ CRi0) (hCRid : 0 < CRid)
    (hd : 0 < min (Cs * depth) (Cb * hab))
    (hHu : Cu * min (Cs * depth) (Cb * hab) ≤ H) (hHD : min (Cs * depth) (Cb * hab) / CD ≤ H)
    (he : emin ≤ e) (he0 : 0 < e)
    (hκ : Cu * min (Cs * depth) (Cb * hab) * Real.sqrt e ≤ κmax)
    (hus : 0 < ustar)
    -- (M1)+(M2) local equilibrium: κu S² = dissipation
    (hP : kappaUPoint Cs Cb Cuun Cu Cuhi CRi0 CRid H depth hab 0 N2above (S ^ 2) emin e κmax Jb Jbe lconv lent *
        S ^ 2 =
      dissipationPoint Cs Cb CDun CD CDhi CRi0 CRid H depth hab 0 N2above (S ^ 2) emin e τ Jb Jbe lconvD lentD)
    -- (M3) constant stress
    (hτ : kappaUPoint Cs Cb Cuun Cu Cuhi CRi0 CRid H depth hab 0 N2above (S ^ 2) emin e κmax Jb Jbe lconv lent *
        S = ustar ^ 2) :
    kappaUPoint Cs Cb Cuun Cu Cuhi CRi0 CRid H depth hab 0 N2above (S ^ 2) emin e κmax Jb Jbe lconv lent =
      (Cu ^ 3 / CD) ^ ((1 : ℝ) / 4) * min (Cs * depth) (Cb * hab) * ustar := by
  obtain ⟨d, hdd⟩ : ∃ d, min (Cs * depth) (Cb * hab) = d := ⟨_, rfl⟩
  rw [hdd] at hd hHu hHD hκ ⊢
  have hK : kappaUPoint Cs Cb Cuun Cu Cuhi CRi0 CRid H depth hab 0 N2above (S ^ 2) emin e κmax Jb Jbe
      lconv lent = Cu * d * Real.sqrt e := by
    unfold kappaUPoint
    rw [richardson_neutral, scaleFn_neutral _ _ _ _ _ hCRi0 hCRid, stableLength_neutral,
      convSelect_neutral, hdd]
    unfold momentumMixingLength turbulentVelocity kappaU
    rw [max_eq_right he, max_eq_left (mul_pos hCu hd).le, min_eq_right hHu, min_eq_left hκ]
  have hE : dissipationPoint Cs Cb CDun CD CDhi CRi0 CRid H depth hab 0 N2above (S ^ 2) emin e τ Jb Jbe
      lconvD lentD = Real.sqrt e / (d / CD) * e := by
    unfold dissipationPoint
    rw [richardson_neutral, scaleFn_neutral _ _ _ _ _ hCRi0 hCRid, stableLength_neutral,
      convSelect_neutral, hdd]
    unfold dissipation dissipationLength dissipationRate
    rw [max_eq_left (div_pos hd hCD).le, min_eq_right hHD, if_neg (not_lt.mpr he0.le), abs_of_pos he0]
  have hr0 : 0 < Real.sqrt e := Real.sqrt_pos.mpr he0
  have her : e = Real.sqrt e ^ 2 := (Real.sq_sqrt he0.le).symm
  have hP' := hP
  rw [hK, hE] at hP'
  have hP'' : Cu * d * Real.sqrt e * S ^ 2 = Real.sqrt e / (d / CD) * Real.sqrt e ^ 2 := by
    rw [← her]
    exact hP'
  have hτ' := hτ
  rw [hK] at hτ'
  have hK4 := neutral_balance_K4 Cu CD d (Real.sqrt e) S ustar hCD hd hr0 hP'' hτ'
  have hq0 : 0 ≤ (Cu ^ 3 / CD) ^ ((1 : ℝ) / 4) := Real.rpow_nonneg (by positivity) _
  rw [hK]
  refine eq_of_pow_four_eq (mul_nonneg (mul_pos hCu hd).le (Real.sqrt_nonneg e))
    (mul_nonneg (mul_nonneg hq0 hd.le) hus.le) ?_
  rw [hK4, show ((Cu ^ 3 / CD) ^ ((1 : ℝ) / 4) * d * ustar) ^ 4 =
      ((Cu ^ 3 / CD) ^ ((1 : ℝ) / 4)) ^ 4 * d ^ 4 * ustar ^ 4 by ring,
    rpow_quarter_pow_four (div_nonneg (pow_nonneg hCu.le 3) hCD.le)]

/-- Certificate (bottom boundary layer): in the near-bottom region (Cᵇ·h ≤ Cˢ·depth) `K_u = Cᵇ (Cˡᵒu³/CˡᵒD)^{1/4} · u★ · h`, i.e. the
von Kármán constant implied by the code is exactly `κ_b = Cᵇ (Cˡᵒu³/CˡᵒD)^{1/4}` (for arbitrary Cᵇ, Cˡᵒu, CˡᵒD > 0). -/
theorem bottom_log_law
    (Cs Cb Cuun Cu Cuhi CDun CD CDhi CRi0 CRid H depth hab N2above S emin e κmax τ Jb Jbe lconv lent lconvD lentD
      ustar : ℝ)
    (hCb : 0 < Cb) (hCu : 0 < Cu) (hCD : 0 < CD) (hCRi0 : 0 ≤ CRi0) (hCRid : 0 < CRid)
    (hhab : 0 < hab) (hbottom : Cb * hab ≤ Cs * depth)
    (hHu : Cu * (Cb * hab) ≤ H) (hHD : Cb * hab / CD ≤ H)
    (he : emin ≤ e) (he0 : 0 < e) (hκ : Cu * (Cb * hab) * Real.sqrt e ≤ κmax) (hus : 0 < ustar)
    (hP : kappaUPoint Cs Cb Cuun Cu Cuhi CRi0 CRid H depth hab 0 N2above (S ^ 2) emin e κmax Jb Jbe lconv lent *
        S ^ 2 =
      dissipationPoint Cs Cb CDun CD CDhi CRi0 CRid H depth hab 0 N2above (S ^ 2) emin e τ Jb Jbe lconvD lentD)
    (hτ : kappaUPoint Cs Cb Cuun Cu Cuhi CRi0 CRid H depth hab 0 N2above (S ^ 2) emin e κmax Jb Jbe lconv lent *
        S = ustar ^ 2) :
    kappaUPoint Cs Cb Cuun Cu Cuhi CRi0 CRid H depth hab 0 N2above (S ^ 2) emin e κmax Jb Jbe lconv lent =
      Cb * (Cu ^ 3 / CD) ^ ((1 : ℝ) / 4) * ustar * hab := by
  have hmin : min (Cs * depth) (Cb * hab) = Cb * hab := min_eq_right hbottom
  have hd : 0 < min (Cs * depth) (Cb * hab) := by
    rw [hmin]
    exact mul_pos hCb hhab
  have h := neutral_wall_law Cs Cb Cuun Cu Cuhi CDun CD CDhi CRi0 CRid H depth hab N2above S emin e κmax τ
    Jb Jbe lconv lent lconvD lentD ustar hCu hCD hCRi0 hCRid hd (by rw [hmin]; exact hHu) (by rw [hmin]; exact hHD)
    he he0 (by rw [hmin]; exact hκ) hus hP hτ
  rw [h, hmin]
  ring

/-- Corollary: under the same conditions the near-bottom shear satisfies the log law `S = u★ / (κ_b h)`. -/
theorem bottom_log_law_shear
    (Cs Cb Cuun Cu Cuhi CDun CD CDhi CRi0 CRid H depth hab N2above S emin e κmax τ Jb Jbe lconv lent lconvD lentD
      ustar : ℝ)
    (hCb : 0 < Cb) (hCu : 0 < Cu) (hCD : 0 < CD) (hCRi0 : 0 ≤ CRi0) (hCRid : 0 < CRid)
    (hhab : 0 < hab) (hbottom : Cb * hab ≤ Cs * depth)
    (hHu : Cu * (Cb * hab) ≤ H) (hHD : Cb * hab / CD ≤ H)
    (he : emin ≤ e) (he0 : 0 < e) (hκ : Cu * (Cb * hab) * Real.sqrt e ≤ κmax) (hus : 0 < ustar)
    (hP : kappaUPoint Cs Cb Cuun Cu Cuhi CRi0 CRid H depth hab 0 N2above (S ^ 2) emin e κmax Jb Jbe lconv lent *
        S ^ 2 =
      dissipationPoint Cs Cb CDun CD CDhi CRi0 CRid H depth hab 0 N2above (S ^ 2) emin e τ Jb Jbe lconvD lentD)
    (hτ : kappaUPoint Cs Cb Cuun Cu Cuhi CRi0 CRid H depth hab 0 N2above (S ^ 2) emin e κmax Jb Jbe lconv lent *
        S = ustar ^ 2) :
    S = ustar / (Cb * (Cu ^ 3 / CD) ^ ((1 : ℝ) / 4) * hab) := by
  have hK := bottom_log_law Cs Cb Cuun Cu Cuhi CDun CD CDhi CRi0 CRid H depth hab N2above S emin e κmax τ
    Jb Jbe lconv lent lconvD lentD ustar hCb hCu hCD hCRi0 hCRid hhab hbottom hHu hHD he he0 hκ hus hP hτ
  rw [hK] at hτ
  have hq : 0 < (Cu ^ 3 / CD) ^ ((1 : ℝ) / 4) := Real.rpow_pos_of_pos (by positivity) _
  have hden : Cb * (Cu ^ 3 / CD) ^ ((1 : ℝ) / 4) * hab ≠ 0 := by positivity
  rw [eq_div_iff hden]
  apply mul_left_cancel₀ hus.ne'
  linear_combination hτ

/-! ## Arithmetic -/

section Arithmetic

-- Side note: check of constants, not counted as a certificate. Default values Cᵇ = 0.967, Cˢ = 1.131, Cˡᵒu = 0.361, CˡᵒD = 1.604
-- (catke_mixing_length.jl:28-29, :34; catke_equation.jl:9). The bounds on κ = C (Cˡᵒu³/CˡᵒD)^{1/4} are written with fourth powers to avoid rpow.

/-- Side note: κ⁴ = C⁴ Cˡᵒu³ / CˡᵒD (formalizes the rpow form as a fourth power). -/
theorem kappa_pow_four (C Cu CD : ℝ) (hCu : 0 < Cu) (hCD : 0 < CD) :
    (C * (Cu ^ 3 / CD) ^ ((1 : ℝ) / 4)) ^ 4 = C ^ 4 * (Cu ^ 3 / CD) := by
  rw [mul_pow, rpow_quarter_pow_four (div_nonneg (pow_nonneg hCu.le 3) hCD.le)]

/-- Side note: the new default Cᵇ = 0.967 gives κ_b ∈ [0.4001, 0.4003] (real value 0.40018). -/
theorem kappa_b_new_bounds :
    (0.4001 : ℝ) ^ 4 * 1.604 ≤ 0.967 ^ 4 * 0.361 ^ 3 ∧ (0.967 : ℝ) ^ 4 * 0.361 ^ 3 ≤ 0.4003 ^ 4 * 1.604 := by
  norm_num

/-- Side note: Cˢ = 1.131 gives κ_s ∈ [0.467, 0.469] (0.468 as stated in the PR). -/
theorem kappa_s_bounds :
    (0.467 : ℝ) ^ 4 * 1.604 ≤ 1.131 ^ 4 * 0.361 ^ 3 ∧ (1.131 : ℝ) ^ 4 * 0.361 ^ 3 ≤ 0.469 ^ 4 * 1.604 := by
  norm_num

/-- Side note: the old value Cᵇ = 0.28 gives κ_b ∈ [0.1158, 0.1159]. -/
theorem kappa_b_old_bounds :
    (0.1158 : ℝ) ^ 4 * 1.604 ≤ 0.28 ^ 4 * 0.361 ^ 3 ∧ (0.28 : ℝ) ^ 4 * 0.361 ^ 3 ≤ 0.1159 ^ 4 * 1.604 := by
  norm_num

/-- Side note: solving back for κ = 0.4 gives Cᵇ = 0.4 (CˡᵒD/Cˡᵒu³)^{1/4} ∈ [0.9665, 0.9675]; 0.967 is its three-significant-digit value. -/
theorem cb_inverse_bounds :
    (0.9665 : ℝ) ^ 4 * 0.361 ^ 3 ≤ 0.4 ^ 4 * 1.604 ∧ (0.4 : ℝ) ^ 4 * 1.604 ≤ 0.9675 ^ 4 * 0.361 ^ 3 := by
  norm_num

/-- Side note: the old formula stated in the PR description, Cᵇ = κ (CˡᵒD/Cˡᵒu)^{1/4} (Cˡᵒu not cubed), gives a value in [0.58, 0.59], inconsistent with the old value 0.28. -/
theorem pr_stated_old_formula_bounds :
    (0.58 : ℝ) ^ 4 * 0.361 ≤ 0.4 ^ 4 * 1.604 ∧ (0.4 : ℝ) ^ 4 * 1.604 ≤ 0.59 ^ 4 * 0.361 := by
  norm_num

/-- Side note: Cᵇ = κ (Cˡᵒu/CˡᵒD)^{1/4} with the two coefficients swapped and not cubed lies in [0.275, 0.276] and rounds to 0.28 at two digits,
but is not strictly equal to 0.28 (0.28⁴·CˡᵒD > 0.4⁴·Cˡᵒu). -/
theorem swapped_formula_bounds :
    (0.275 : ℝ) ^ 4 * 1.604 ≤ 0.4 ^ 4 * 0.361 ∧ (0.4 : ℝ) ^ 4 * 0.361 ≤ 0.276 ^ 4 * 1.604 ∧
      (0.4 : ℝ) ^ 4 * 0.361 < 0.28 ^ 4 * 1.604 := by
  norm_num

end Arithmetic

end

end Oceananigans6024
