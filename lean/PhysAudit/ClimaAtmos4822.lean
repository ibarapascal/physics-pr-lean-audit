-- ClimaAtmos4822: CliMA/ClimaAtmos.jl PR #4822 "Compute horizontal gradient along isentropes in SGS variance"
-- (merge commit ba2bf1225fc48058d709853c065b1db21a412204) Lean certificate.
-- Transcribes `isentropic_cell_cap` (:253-257), `isentropic_slope_ratio` (:270), `isentropic_gradient_invariant` (:279-280) of src/parameterized_tendencies/microphysics/cloud_fraction.jl,
-- and the combination of θ′θ′, q′q′, θ′q′ in `set_covariance_cache!` (:408-607),
-- and proves / refutes three groups of invariants over all reals (a: gauge invariance, b: no double counting, c: variance consistency).
import Mathlib

/-!
# ClimaAtmos.jl #4822 — isentropic SGS variance

Transcription conventions
* Real-number semantics: Julia `max(·, floor)` and `if` branches are transcribed as written as `max` / `if`; floating-point rounding is out of scope for this file (the diff script
  `replay/scripts/climaatmos4822_diff.jl` handles the comparison).
* The horizontal gradient invariants `inv_q = |∇_h q|²`, `inv_θq = ∇_h θ · ∇_h q`, `inv_θ = |∇_h θ|²` enter as three real inputs.
  In the code they are given by `hgrad_invariant!` / `hgrad_cross_invariant!` (:294-314, nodal product → `lumpedₕ` restriction → weighted DSS);
  the triple-form theorem of (a) only uses that "the triple transforms as a quadratic form under ∇q ↦ ∇q + b∇θ", which also holds for a triple after such linear filtering
  (the linearity of the filter is a speculation, not formalized in Lean); a pointwise vector version (`dot2`) is also given.
* Vertical gradient: in the code `isentropic_slope_ratio` / `isentropic_cell_cap` use `projected_vector_data(C3, ᶜgradᵥ_ψ, ᶜlg)`
  (:512-513, :524), and the vertical block uses the dot product of `Geometry.WVector(ᶜgradᵥ_ψ)` (:556-569). On a terrain-free grid both are the physical vertical component;
  this file substitutes the same real numbers `θz`, `qz` in both places (hypothesis, see the report). The horizontal block uses the gradient of `ᶜq_tot_nonneg` and the vertical block uses
  `ᶜgradᵥ_q_tot`, which are the same field when q_tot ≥ 0 (hypothesis).
* Julia numerical parameters: `θz_min = sgs_variance_isentropic_min_dtheta_dz`, `c_g`, `c_Δx`, `c_Δz`, `coeff`, `ρ_turb`,
  `sgs_variance_max_rel_std`, `eps(FT)` are all explicit parameters.
-/

namespace PhysAudit.ClimaAtmos4822

noncomputable section

/-! ## §1 Transcription of the pure functions -/

/-- Pointwise horizontal dot product (nodal-value model of `hgrad_invariant!` / `hgrad_cross_invariant!`).
```
cloud_fraction.jl:295  @. ᶜinv = lumpedₕ(norm_sqr(gradₕ(ᶜψ)))
cloud_fraction.jl:310  @. ᶜinv = lumpedₕ(dot(Geometry.Contravariant12Vector(gradₕ(ᶜψ1)), gradₕ(ᶜψ2)))
```
-/
def dot2 (u₁ u₂ v₁ v₂ : ℝ) : ℝ := u₁ * v₁ + u₂ * v₂

/-- ```
cloud_fraction.jl:254      s = sqrt(max(inv_θ, zero(inv_θ))) / max(∂θ∂z, θz_min)
``` -/
def cellSlope (invθ θz θzMin : ℝ) : ℝ := Real.sqrt (max invθ 0) / max θz θzMin

/-- ```
cloud_fraction.jl:255      x = s * Δx / Δz
``` -/
def cellX (invθ θz θzMin Δx Δz : ℝ) : ℝ := cellSlope invθ θz θzMin * Δx / Δz

/-- ```
cloud_fraction.jl:253  @inline function isentropic_cell_cap(inv_θ, ∂θ∂z, θz_min, Δx, Δz)
cloud_fraction.jl:254      s = sqrt(max(inv_θ, zero(inv_θ))) / max(∂θ∂z, θz_min)
cloud_fraction.jl:255      x = s * Δx / Δz
cloud_fraction.jl:256      return one(x) / (one(x) + x^2)
cloud_fraction.jl:257  end
``` -/
def isentropicCellCap (invθ θz θzMin Δx Δz : ℝ) : ℝ :=
  1 / (1 + cellX invθ θz θzMin Δx Δz ^ 2)

/-- ```
cloud_fraction.jl:270  @inline isentropic_slope_ratio(∂q∂z, ∂θ∂z, θz_min) = ∂q∂z / max(∂θ∂z, θz_min)
``` -/
def isentropicSlopeRatio (qz θz θzMin : ℝ) : ℝ := qz / max θz θzMin

/-- ```
cloud_fraction.jl:279  @inline isentropic_gradient_invariant(inv_q, inv_θq, inv_θ, r) =
cloud_fraction.jl:280      max(inv_q - 2 * r * inv_θq + r^2 * inv_θ, zero(inv_q))
``` -/
def isentropicGradientInvariant (invq invθq invθ r : ℝ) : ℝ :=
  max (invq - 2 * r * invθq + r ^ 2 * invθ) 0

/-- The expanded expression inside `max` (only names the first argument at :280, structure unchanged). -/
def invArg (invq invθq invθ r : ℝ) : ℝ := invq - 2 * r * invθq + r ^ 2 * invθ

theorem isentropicGradientInvariant_eq (invq invθq invθ r : ℝ) :
    isentropicGradientInvariant invq invθq invθ r = max (invArg invq invθq invθ r) 0 := rfl

/-! ## §2 The covariance combination of `set_covariance_cache!` (θ basis, steps 1–3; then to the T basis, the upper bound, and the diagnosed correlation) -/

/-- Parameters (`CAP.*`).
```
cloud_fraction.jl:426  c_g = CAP.sgs_variance_geometric_coeff(p.params)
cloud_fraction.jl:427  c_Δx = CAP.sgs_variance_horizontal_scale_factor(p.params)
cloud_fraction.jl:428  c_Δz = CAP.sgs_variance_vertical_scale_factor(p.params)
cloud_fraction.jl:443  coeff = CAP.diagnostic_covariance_coeff(p.params)
cloud_fraction.jl:467  ρ_turb = FT(correlation_Tq(p.params))
cloud_fraction.jl:505  θz_min = CAP.sgs_variance_isentropic_min_dtheta_dz(p.params)
``` -/
structure SGSParams where
  cg : ℝ
  cΔx : ℝ
  cΔz : ℝ
  coeff : ℝ
  ρturb : ℝ
  θzMin : ℝ

/-- Input at a single grid point: the horizontal invariant triple, the vertical gradients, the grid scales, the mixing length.
`Δx` = `Δx_h` (:432-436 `Spaces.node_horizontal_length_scale`), `Δz` = `ᶜΔz` (:491), `lmix` = `ᶜl_mix` (:442). -/
structure CellInputs where
  invq : ℝ
  invθq : ℝ
  invθ : ℝ
  qz : ℝ
  θz : ℝ
  Δx : ℝ
  Δz : ℝ
  lmix : ℝ

/-- θ-basis covariances (`ᶜq′q′`, `ᶜT′T′` (θ′θ′ is stored before the θ→T transform), `ᶜT′q′`). -/
structure Cov where
  qq : ℝ
  θθ : ℝ
  θq : ℝ

/-- ```
cloud_fraction.jl:444  cov_from_grad(C, L, ∇Φ, ∇Ψ) = 2 * C * L^2 * dot(∇Φ, ∇Ψ)
```
(`∇Φ`, `∇Ψ` are single-component `WVector`s, so the dot product = the product of components.) -/
def covFromGrad (C L a b : ℝ) : ℝ := 2 * C * L ^ 2 * (a * b)

/-- ```
cloud_fraction.jl:456  @. ᶜq′q′ = cov_from_grad(
cloud_fraction.jl:457      coeff, ᶜl_mix,
cloud_fraction.jl:458      Geometry.WVector(ᶜgradᵥ_q_tot),
cloud_fraction.jl:459      Geometry.WVector(ᶜgradᵥ_q_tot),
cloud_fraction.jl:460  )
``` -/
def turbQQ (P : SGSParams) (G : CellInputs) : ℝ := covFromGrad P.coeff G.lmix G.qz G.qz

/-- ```
cloud_fraction.jl:461  @. ᶜT′T′ = cov_from_grad(  # holds θ′θ′ until the ∂T/∂θ transform below.
cloud_fraction.jl:462      coeff, ᶜl_mix,
cloud_fraction.jl:463      Geometry.WVector(ᶜgradᵥ_θ_liq_ice),
cloud_fraction.jl:464      Geometry.WVector(ᶜgradᵥ_θ_liq_ice),
cloud_fraction.jl:465  )
``` -/
def turbθθ (P : SGSParams) (G : CellInputs) : ℝ := covFromGrad P.coeff G.lmix G.θz G.θz

/-- ```
cloud_fraction.jl:466  if diag_corr
cloud_fraction.jl:467      ρ_turb = FT(correlation_Tq(p.params))
cloud_fraction.jl:468      @. p.precomputed.ᶜT′q′ = ρ_turb * sqrt(ᶜT′T′ * ᶜq′q′)
cloud_fraction.jl:469  end
``` -/
def turbθq (P : SGSParams) (G : CellInputs) : ℝ :=
  P.ρturb * Real.sqrt (turbθθ P G * turbQQ P G)

/-- ```
cloud_fraction.jl:429  use_geometric = !iszero(c_Δx)
cloud_fraction.jl:431  geo_h = if use_geometric
cloud_fraction.jl:432      Δx_h = FT(...)
cloud_fraction.jl:437      c_g * (c_Δx * Δx_h)^2
cloud_fraction.jl:438  else
cloud_fraction.jl:439      zero(FT)
cloud_fraction.jl:440  end
``` -/
def geoH (P : SGSParams) (G : CellInputs) : ℝ :=
  if P.cΔx = 0 then 0 else P.cg * (P.cΔx * G.Δx) ^ 2

/-- ```
cloud_fraction.jl:507  @. ᶜinv_q = isentropic_gradient_invariant(
cloud_fraction.jl:508      ᶜinv_q,
cloud_fraction.jl:509      ᶜinv_θq,
cloud_fraction.jl:510      ᶜinv_θ,
cloud_fraction.jl:511      isentropic_slope_ratio(
cloud_fraction.jl:512          projected_vector_data(C3, ᶜgradᵥ_q_tot, ᶜlg),
cloud_fraction.jl:513          projected_vector_data(C3, ᶜgradᵥ_θ_liq_ice, ᶜlg),
cloud_fraction.jl:514          θz_min,
cloud_fraction.jl:515      ),
cloud_fraction.jl:516  )
``` -/
def invIsen (P : SGSParams) (G : CellInputs) : ℝ :=
  isentropicGradientInvariant G.invq G.invθq G.invθ (isentropicSlopeRatio G.qz G.θz P.θzMin)

/-- ```
cloud_fraction.jl:493  if use_geometric
cloud_fraction.jl:522      @. ᶜpart = isentropic_cell_cap(
cloud_fraction.jl:523          ᶜinv_θ,
cloud_fraction.jl:524          projected_vector_data(C3, ᶜgradᵥ_θ_liq_ice, ᶜlg),
cloud_fraction.jl:525          θz_min,
cloud_fraction.jl:526          Δx_h,
cloud_fraction.jl:527          ᶜΔz,
cloud_fraction.jl:528      )
cloud_fraction.jl:531  else
cloud_fraction.jl:532      @. ᶜpart = one(FT)
cloud_fraction.jl:533  end
```
Note: what is passed here is `Δx_h`, whereas the horizontal block scale is `c_Δx * Δx_h` (:437). -/
def part (P : SGSParams) (G : CellInputs) : ℝ :=
  if P.cΔx = 0 then 1 else isentropicCellCap G.invθ G.θz P.θzMin G.Δx G.Δz

/-- Increment of the horizontal geometric block to q′q′.
```
cloud_fraction.jl:493  if use_geometric
cloud_fraction.jl:529      @. ᶜinv_q *= ᶜpart
cloud_fraction.jl:530      @. ᶜq′q′ += geo_h * ᶜinv_q
``` -/
def horizQQ (P : SGSParams) (G : CellInputs) : ℝ :=
  if P.cΔx = 0 then 0 else geoH P G * (invIsen P G * part P G)

/-- ```
cloud_fraction.jl:552      geo_v(Δz) = c_g * (c_Δz * Δz)^2                 # per-cell coefficient
``` -/
def geoV (P : SGSParams) (Δz : ℝ) : ℝ := P.cg * (P.cΔz * Δz) ^ 2

/-- ```
cloud_fraction.jl:551  if use_vertical_grad
cloud_fraction.jl:553      @. ᶜq′q′ +=
cloud_fraction.jl:554          geo_v(ᶜΔz) * ᶜpart *
cloud_fraction.jl:555          dot(
cloud_fraction.jl:556              Geometry.WVector(ᶜgradᵥ_q_tot),
cloud_fraction.jl:557              Geometry.WVector(ᶜgradᵥ_q_tot),
cloud_fraction.jl:558          )
``` -/
def vertQQ (P : SGSParams) (G : CellInputs) : ℝ :=
  if P.cΔz = 0 then 0 else geoV P G.Δz * part P G * (G.qz * G.qz)

/-- ```
cloud_fraction.jl:559      @. ᶜT′T′ +=
cloud_fraction.jl:560          geo_v(ᶜΔz) * dot(
cloud_fraction.jl:561              Geometry.WVector(ᶜgradᵥ_θ_liq_ice),
cloud_fraction.jl:562              Geometry.WVector(ᶜgradᵥ_θ_liq_ice),
cloud_fraction.jl:563          )
``` -/
def vertθθ (P : SGSParams) (G : CellInputs) : ℝ :=
  if P.cΔz = 0 then 0 else geoV P G.Δz * (G.θz * G.θz)

/-- ```
cloud_fraction.jl:564      if diag_corr
cloud_fraction.jl:565          @. p.precomputed.ᶜT′q′ +=
cloud_fraction.jl:566              geo_v(ᶜΔz) * ᶜpart *
cloud_fraction.jl:567              dot(
cloud_fraction.jl:568                  Geometry.WVector(ᶜgradᵥ_θ_liq_ice),
cloud_fraction.jl:569                  Geometry.WVector(ᶜgradᵥ_q_tot),
cloud_fraction.jl:570              )
``` -/
def vertθq (P : SGSParams) (G : CellInputs) : ℝ :=
  if P.cΔz = 0 then 0 else geoV P G.Δz * part P G * (G.θz * G.qz)

/-- The θ-basis covariances at the end of steps 1–3 (before the θ→T transform, before :574): turbulence block + horizontal block + vertical block.
(The `θq` component corresponds to `ᶜT′q′` when `diag_corr = true`.) -/
def covTheta (P : SGSParams) (G : CellInputs) : Cov where
  qq := turbQQ P G + horizQQ P G + vertQQ P G
  θθ := turbθθ P G + vertθθ P G
  θq := turbθq P G + vertθq P G

/-- ```
cloud_fraction.jl:579  @. ᶜT′T′ = ᶜ∂T_∂θ^2 * ᶜT′T′
cloud_fraction.jl:580  if diag_corr
cloud_fraction.jl:581      @. p.precomputed.ᶜT′q′ = ᶜ∂T_∂θ * p.precomputed.ᶜT′q′
cloud_fraction.jl:582  end
``` -/
def toT (J : ℝ) (c : Cov) : Cov where
  qq := c.qq
  θθ := J ^ 2 * c.θθ
  θq := J * c.θq

/-- ```
cloud_fraction.jl:594  if diag_corr
cloud_fraction.jl:595      @. p.precomputed.ᶜT′q′ *= sqrt(min(
cloud_fraction.jl:596          one(FT),
cloud_fraction.jl:597          (r_max * ᶜq_tot_nonneg)^2 / max(ᶜq′q′, eps(FT)^2),
cloud_fraction.jl:598      ))
cloud_fraction.jl:599  end
cloud_fraction.jl:600  @. ᶜq′q′ = min(ᶜq′q′, (r_max * ᶜq_tot_nonneg)^2)
```
(`ε` = `eps(FT)`. `θq` uses the q′q′ before the upper bound, consistent with the line order in the code.) -/
def applyBound (rmax qtot ε : ℝ) (c : Cov) : Cov where
  qq := min c.qq ((rmax * qtot) ^ 2)
  θθ := c.θθ
  θq := c.θq * Real.sqrt (min 1 ((rmax * qtot) ^ 2 / max c.qq (ε ^ 2)))

/-- ```
cloud_fraction.jl:334      v = max(T′T′, zero(FT)) * max(q′q′, zero(FT))
``` -/
def corrV (TT qq : ℝ) : ℝ := max TT 0 * max qq 0

/-- Julia Base (base/math.jl): `clamp(x, lo, hi) = ifelse(x > hi, hi, ifelse(x < lo, lo, x))`. -/
def juliaClamp (x lo hi : ℝ) : ℝ := if x > hi then hi else if x < lo then lo else x

/-- ```
cloud_fraction.jl:332  @inline function sgs_correlation_Tq(T′T′, q′q′, T′q′, fallback)
cloud_fraction.jl:333      FT = typeof(T′T′)
cloud_fraction.jl:334      v = max(T′T′, zero(FT)) * max(q′q′, zero(FT))
cloud_fraction.jl:335      return ifelse(
cloud_fraction.jl:336          v > eps(FT)^2,
cloud_fraction.jl:337          clamp(T′q′ / sqrt(max(v, eps(FT)^2)), -one(FT), one(FT)),
cloud_fraction.jl:338          fallback,
cloud_fraction.jl:339      )
cloud_fraction.jl:340  end
``` -/
def sgsCorrelationTq (TT qq Tq fallback ε : ℝ) : ℝ :=
  if corrV TT qq > ε ^ 2 then
    juliaClamp (Tq / Real.sqrt (max (corrV TT qq) (ε ^ 2))) (-1) 1
  else fallback

/-- Covariances after steps 4–5 (`J` = `ᶜ∂T_∂θ`, `rmax` = `sgs_variance_max_rel_std`, `qtot` = `ᶜq_tot_nonneg`). -/
def finalCov (P : SGSParams) (G : CellInputs) (J rmax qtot ε : ℝ) : Cov :=
  applyBound rmax qtot ε (toT J (covTheta P G))

/-- Diagnosed correlation coefficient (:605 `set_tq_correlation!` → :360-364 `DiagnosedTqCorrelation`). -/
def corrFinal (P : SGSParams) (G : CellInputs) (J rmax qtot ε fallback : ℝ) : ℝ :=
  sgsCorrelationTq (finalCov P G J rmax qtot ε).θθ (finalCov P G J rmax qtot ε).qq
    (finalCov P G J rmax qtot ε).θq fallback ε

theorem covTheta_qq (P : SGSParams) (G : CellInputs) :
    (covTheta P G).qq = turbQQ P G + horizQQ P G + vertQQ P G := rfl
theorem covTheta_θθ (P : SGSParams) (G : CellInputs) :
    (covTheta P G).θθ = turbθθ P G + vertθθ P G := rfl
theorem covTheta_θq (P : SGSParams) (G : CellInputs) :
    (covTheta P G).θq = turbθq P G + vertθq P G := rfl
theorem toT_qq (J : ℝ) (c : Cov) : (toT J c).qq = c.qq := rfl
theorem toT_θθ (J : ℝ) (c : Cov) : (toT J c).θθ = J ^ 2 * c.θθ := rfl
theorem toT_θq (J : ℝ) (c : Cov) : (toT J c).θq = J * c.θq := rfl
theorem applyBound_qq (rmax qtot ε : ℝ) (c : Cov) :
    (applyBound rmax qtot ε c).qq = min c.qq ((rmax * qtot) ^ 2) := rfl
theorem applyBound_θθ (rmax qtot ε : ℝ) (c : Cov) : (applyBound rmax qtot ε c).θθ = c.θθ := rfl
theorem applyBound_θq (rmax qtot ε : ℝ) (c : Cov) :
    (applyBound rmax qtot ε c).θq
      = c.θq * Real.sqrt (min 1 ((rmax * qtot) ^ 2 / max c.qq (ε ^ 2))) := rfl

/-! ## §3 Basic properties of `isentropic_cell_cap` -/

theorem cap_mul_one_add_sq (invθ θz θzMin Δx Δz : ℝ) :
    isentropicCellCap invθ θz θzMin Δx Δz * (1 + cellX invθ θz θzMin Δx Δz ^ 2) = 1 := by
  have h : (0 : ℝ) < 1 + cellX invθ θz θzMin Δx Δz ^ 2 := by positivity
  unfold isentropicCellCap
  rw [one_div, inv_mul_cancel₀ h.ne']

/-- 0 < f (for all real inputs). -/
theorem cap_pos (invθ θz θzMin Δx Δz : ℝ) : 0 < isentropicCellCap invθ θz θzMin Δx Δz := by
  unfold isentropicCellCap
  positivity

/-- f ≤ 1 (for all real inputs). -/
theorem cap_le_one (invθ θz θzMin Δx Δz : ℝ) : isentropicCellCap invθ θz θzMin Δx Δz ≤ 1 := by
  have hf := cap_mul_one_add_sq invθ θz θzMin Δx Δz
  have h := mul_nonneg (cap_pos invθ θz θzMin Δx Δz).le (sq_nonneg (cellX invθ θz θzMin Δx Δz))
  nlinarith

theorem cellSlope_pos (invθ θz θzMin : ℝ) (hmin : 0 < θzMin) (hinv : 0 < invθ) :
    0 < cellSlope invθ θz θzMin := by
  unfold cellSlope
  exact div_pos (Real.sqrt_pos.mpr (lt_max_of_lt_left hinv)) (lt_max_of_lt_right hmin)

theorem cellX_pos (invθ θz θzMin Δx Δz : ℝ) (hmin : 0 < θzMin) (hinv : 0 < invθ)
    (hΔx : 0 < Δx) (hΔz : 0 < Δz) : 0 < cellX invθ θz θzMin Δx Δz := by
  unfold cellX
  exact div_pos (mul_pos (cellSlope_pos invθ θz θzMin hmin hinv) hΔx) hΔz

/-- f < 1 ⇔ the isentropic surface is tilted: strictly less than 1 when |∇_h θ|² > 0, Δx > 0 (other quantities positive). -/
theorem cap_lt_one (invθ θz θzMin Δx Δz : ℝ) (hmin : 0 < θzMin) (hinv : 0 < invθ)
    (hΔx : 0 < Δx) (hΔz : 0 < Δz) : isentropicCellCap invθ θz θzMin Δx Δz < 1 := by
  have hf := cap_mul_one_add_sq invθ θz θzMin Δx Δz
  have hX := cellX_pos invθ θz θzMin Δx Δz hmin hinv hΔx hΔz
  have h := mul_pos (cap_pos invθ θz θzMin Δx Δz) (pow_pos hX 2)
  nlinarith

theorem cellX_mul_Δz (invθ θz θzMin Δx Δz : ℝ) (hΔz : Δz ≠ 0) :
    cellX invθ θz θzMin Δx Δz * Δz = cellSlope invθ θz θzMin * Δx := by
  unfold cellX
  field_simp

theorem cellSlope_sq (invθ θz θzMin : ℝ) (hinv : 0 ≤ invθ) :
    cellSlope invθ θz θzMin ^ 2 = invθ / max θz θzMin ^ 2 := by
  unfold cellSlope
  rw [div_pow, max_eq_left hinv, Real.sq_sqrt hinv]

/-- Algebraic core: f(1+X²)=1, XΔz=SΔx ⇒ S²(fΔx²) + fΔz² = Δz². -/
theorem core_cap_identity (f S X Δx Δz : ℝ) (hf : f * (1 + X ^ 2) = 1)
    (hxz : X * Δz = S * Δx) : S ^ 2 * (f * Δx ^ 2) + f * Δz ^ 2 = Δz ^ 2 := by
  linear_combination Δz ^ 2 * hf - f * (S * Δx + X * Δz) * hxz

/-- **(b) Identity** `(sℓ)² + fΔz² = Δz²`, where `ℓ² = fΔx²` and `s` is the `s` in the code (:254). For any θz (including the floor branch). -/
theorem cap_identity (invθ θz θzMin Δx Δz : ℝ) (hΔz : 0 < Δz) :
    cellSlope invθ θz θzMin ^ 2 * (isentropicCellCap invθ θz θzMin Δx Δz * Δx ^ 2)
      + isentropicCellCap invθ θz θzMin Δx Δz * Δz ^ 2 = Δz ^ 2 :=
  core_cap_identity _ _ _ Δx Δz (cap_mul_one_add_sq invθ θz θzMin Δx Δz)
    (cellX_mul_Δz invθ θz θzMin Δx Δz hΔz.ne')

/-! ## §4 Theorem (a): gauge invariance and q = F(θ) ⇒ I = 0 -/

/-- The expanded expression equals the docstring's `|∇_h q − r ∇_h θ|²` (the max lower bound has no effect in exact arithmetic). -/
theorem isentropicGradientInvariant_eq_normSq (gq₁ gq₂ gθ₁ gθ₂ r : ℝ) :
    isentropicGradientInvariant (dot2 gq₁ gq₂ gq₁ gq₂) (dot2 gθ₁ gθ₂ gq₁ gq₂)
        (dot2 gθ₁ gθ₂ gθ₁ gθ₂) r
      = (gq₁ - r * gθ₁) ^ 2 + (gq₂ - r * gθ₂) ^ 2 := by
  have h : invArg (dot2 gq₁ gq₂ gq₁ gq₂) (dot2 gθ₁ gθ₂ gq₁ gq₂) (dot2 gθ₁ gθ₂ gθ₁ gθ₂) r
      = (gq₁ - r * gθ₁) ^ 2 + (gq₂ - r * gθ₂) ^ 2 := by
    unfold invArg dot2
    ring
  rw [isentropicGradientInvariant_eq, h]
  exact max_eq_left (by positivity)

/-- **(a) General shift formula** (any θz): under q ↦ q + bθ the quantity inside `max` changes exactly by
`(1 − k)·b·(2 inv_θq + b(1 − k) inv_θ − 2 r inv_θ)`, `k = θz / max(θz, θz_min)`. -/
theorem gauge_shift_general (invq invθq invθ qz θz θzMin b : ℝ) :
    invArg (invq + 2 * b * invθq + b ^ 2 * invθ) (invθq + b * invθ) invθ
        (isentropicSlopeRatio (qz + b * θz) θz θzMin)
      = invArg invq invθq invθ (isentropicSlopeRatio qz θz θzMin)
        + (1 - θz / max θz θzMin) * b
          * (2 * invθq + b * (1 - θz / max θz θzMin) * invθ
              - 2 * isentropicSlopeRatio qz θz θzMin * invθ) := by
  unfold invArg isentropicSlopeRatio
  ring

/-- **(a) Gauge invariance** (triple form): for θz ≥ θz_min, `I(∇q + b∇θ, q_z + bθ_z) = I(∇q, q_z)`. -/
theorem gauge_invariant (invq invθq invθ qz θz θzMin b : ℝ) (hmin : 0 < θzMin)
    (hθ : θzMin ≤ θz) :
    isentropicGradientInvariant (invq + 2 * b * invθq + b ^ 2 * invθ) (invθq + b * invθ) invθ
        (isentropicSlopeRatio (qz + b * θz) θz θzMin)
      = isentropicGradientInvariant invq invθq invθ (isentropicSlopeRatio qz θz θzMin) := by
  have hθ0 : θz ≠ 0 := (lt_of_lt_of_le hmin hθ).ne'
  simp only [isentropicGradientInvariant_eq]
  rw [gauge_shift_general, max_eq_left hθ, div_self hθ0]
  congr 1
  ring

/-- **(a) Gauge invariance** (pointwise vector form). -/
theorem gauge_invariant_vec (gq₁ gq₂ gθ₁ gθ₂ qz θz θzMin b : ℝ) (hmin : 0 < θzMin)
    (hθ : θzMin ≤ θz) :
    isentropicGradientInvariant
        (dot2 (gq₁ + b * gθ₁) (gq₂ + b * gθ₂) (gq₁ + b * gθ₁) (gq₂ + b * gθ₂))
        (dot2 gθ₁ gθ₂ (gq₁ + b * gθ₁) (gq₂ + b * gθ₂)) (dot2 gθ₁ gθ₂ gθ₁ gθ₂)
        (isentropicSlopeRatio (qz + b * θz) θz θzMin)
      = isentropicGradientInvariant (dot2 gq₁ gq₂ gq₁ gq₂) (dot2 gθ₁ gθ₂ gq₁ gq₂)
        (dot2 gθ₁ gθ₂ gθ₁ gθ₂) (isentropicSlopeRatio qz θz θzMin) := by
  have h1 : dot2 (gq₁ + b * gθ₁) (gq₂ + b * gθ₂) (gq₁ + b * gθ₁) (gq₂ + b * gθ₂)
      = dot2 gq₁ gq₂ gq₁ gq₂ + 2 * b * dot2 gθ₁ gθ₂ gq₁ gq₂ + b ^ 2 * dot2 gθ₁ gθ₂ gθ₁ gθ₂ := by
    unfold dot2
    ring
  have h2 : dot2 gθ₁ gθ₂ (gq₁ + b * gθ₁) (gq₂ + b * gθ₂)
      = dot2 gθ₁ gθ₂ gq₁ gq₂ + b * dot2 gθ₁ gθ₂ gθ₁ gθ₂ := by
    unfold dot2
    ring
  rw [h1, h2]
  exact gauge_invariant _ _ _ qz θz θzMin b hmin hθ

/-- **(a) q = F(θ) ⇒ I = 0** (triple form, θz ≥ θz_min):
`inv_q = F′² inv_θ`, `inv_θq = F′ inv_θ`, `q_z = F′ θ_z`. -/
theorem invariant_zero_of_func (invθ θz θzMin F : ℝ) (hmin : 0 < θzMin) (hθ : θzMin ≤ θz) :
    isentropicGradientInvariant (F ^ 2 * invθ) (F * invθ) invθ
      (isentropicSlopeRatio (F * θz) θz θzMin) = 0 := by
  have hθ0 : θz ≠ 0 := (lt_of_lt_of_le hmin hθ).ne'
  have hr : isentropicSlopeRatio (F * θz) θz θzMin = F := by
    unfold isentropicSlopeRatio
    rw [max_eq_left hθ, mul_div_assoc, div_self hθ0, mul_one]
  rw [isentropicGradientInvariant_eq, hr]
  exact max_eq_right (le_of_eq (by unfold invArg; ring))

/-- **(a) q = F(θ) ⇒ I = 0** (pointwise vector form: ∇q = F′∇θ, q_z = F′θ_z). -/
theorem invariant_zero_of_func_vec (gθ₁ gθ₂ θz θzMin F : ℝ) (hmin : 0 < θzMin)
    (hθ : θzMin ≤ θz) :
    isentropicGradientInvariant (dot2 (F * gθ₁) (F * gθ₂) (F * gθ₁) (F * gθ₂))
        (dot2 gθ₁ gθ₂ (F * gθ₁) (F * gθ₂)) (dot2 gθ₁ gθ₂ gθ₁ gθ₂)
        (isentropicSlopeRatio (F * θz) θz θzMin) = 0 := by
  have h1 : dot2 (F * gθ₁) (F * gθ₂) (F * gθ₁) (F * gθ₂) = F ^ 2 * dot2 gθ₁ gθ₂ gθ₁ gθ₂ := by
    unfold dot2
    ring
  have h2 : dot2 gθ₁ gθ₂ (F * gθ₁) (F * gθ₂) = F * dot2 gθ₁ gθ₂ gθ₁ gθ₂ := by
    unfold dot2
    ring
  rw [h1, h2]
  exact invariant_zero_of_func _ θz θzMin F hmin hθ

/-- **(a) Exact value on the floor branch** (θz < θz_min): for q = F(θ), `I = inv_θ · F′² · (1 − θz/θz_min)²`. -/
theorem invariant_of_func_below_floor (invθ θz θzMin F : ℝ) (hinv : 0 ≤ invθ)
    (hlt : θz < θzMin) :
    isentropicGradientInvariant (F ^ 2 * invθ) (F * invθ) invθ
        (isentropicSlopeRatio (F * θz) θz θzMin)
      = invθ * F ^ 2 * (1 - θz / θzMin) ^ 2 := by
  have hE : invArg (F ^ 2 * invθ) (F * invθ) invθ (isentropicSlopeRatio (F * θz) θz θzMin)
      = invθ * F ^ 2 * (1 - θz / θzMin) ^ 2 := by
    unfold invArg isentropicSlopeRatio
    rw [max_eq_right hlt.le]
    ring
  rw [isentropicGradientInvariant_eq, hE]
  exact max_eq_left (mul_nonneg (mul_nonneg hinv (sq_nonneg F)) (sq_nonneg _))

/-- **(a) Counterexample (universal form)**: for θz < θz_min, ∇θ ≠ 0, F′ ≠ 0, the I of q = F(θ) is strictly positive. -/
theorem invariant_pos_of_func_below_floor (invθ θz θzMin F : ℝ) (hmin : 0 < θzMin)
    (hinv : 0 < invθ) (hF : F ≠ 0) (hlt : θz < θzMin) :
    0 < isentropicGradientInvariant (F ^ 2 * invθ) (F * invθ) invθ
      (isentropicSlopeRatio (F * θz) θz θzMin) := by
  rw [invariant_of_func_below_floor invθ θz θzMin F hinv.le hlt]
  have hk : θz / θzMin < 1 := (div_lt_iff₀ hmin).mpr (by linarith)
  have h1 : (0 : ℝ) < (1 - θz / θzMin) ^ 2 := pow_pos (by linarith) 2
  have h2 : (0 : ℝ) < F ^ 2 := lt_of_le_of_ne (sq_nonneg F) (pow_ne_zero 2 hF).symm
  exact mul_pos (mul_pos hinv h2) h1

/-- **(a) Gauge invariance fails on the floor branch**: for θz < θz_min, shifting q = 0 to q = bθ (b ≠ 0) changes I. -/
theorem gauge_fails_below_floor (invθ θz θzMin b : ℝ) (hmin : 0 < θzMin) (hinv : 0 < invθ)
    (hb : b ≠ 0) (hlt : θz < θzMin) :
    isentropicGradientInvariant (0 + 2 * b * 0 + b ^ 2 * invθ) (0 + b * invθ) invθ
        (isentropicSlopeRatio (0 + b * θz) θz θzMin)
      ≠ isentropicGradientInvariant 0 0 invθ (isentropicSlopeRatio 0 θz θzMin) := by
  have h0 : isentropicGradientInvariant 0 0 invθ (isentropicSlopeRatio 0 θz θzMin) = 0 := by
    simp [isentropicGradientInvariant, isentropicSlopeRatio]
  rw [h0]
  simp only [mul_zero, zero_add]
  exact (invariant_pos_of_func_below_floor invθ θz θzMin b hmin hinv hb hlt).ne'

/-! ## §5 Theorem (b): no double counting -/

theorem part_eq (P : SGSParams) (G : CellInputs) (hcx : P.cΔx ≠ 0) :
    part P G = isentropicCellCap G.invθ G.θz P.θzMin G.Δx G.Δz := by
  unfold part
  rw [if_neg hcx]

theorem part_pos (P : SGSParams) (G : CellInputs) : 0 < part P G := by
  unfold part
  split_ifs
  · exact one_pos
  · exact cap_pos _ _ _ _ _

theorem part_le_one (P : SGSParams) (G : CellInputs) : part P G ≤ 1 := by
  unfold part
  split_ifs
  · exact le_rfl
  · exact cap_le_one _ _ _ _ _

theorem part_lt_one (P : SGSParams) (G : CellInputs) (hcx : P.cΔx ≠ 0) (hmin : 0 < P.θzMin)
    (hinv : 0 < G.invθ) (hΔx : 0 < G.Δx) (hΔz : 0 < G.Δz) : part P G < 1 := by
  rw [part_eq P G hcx]
  exact cap_lt_one _ _ _ _ _ hmin hinv hΔx hΔz

theorem horizQQ_eq (P : SGSParams) (G : CellInputs) (hcx : P.cΔx ≠ 0) :
    horizQQ P G = P.cg * (P.cΔx * G.Δx) ^ 2 * (invIsen P G * part P G) := by
  simp only [horizQQ, geoH, if_neg hcx]

theorem vertQQ_eq_any (P : SGSParams) (G : CellInputs) :
    vertQQ P G = geoV P G.Δz * part P G * (G.qz * G.qz) := by
  unfold vertQQ
  split_ifs with h
  · simp [geoV, h]
  · rfl

theorem vertθθ_eq_any (P : SGSParams) (G : CellInputs) :
    vertθθ P G = geoV P G.Δz * (G.θz * G.θz) := by
  unfold vertθθ
  split_ifs with h
  · simp [geoV, h]
  · rfl

theorem vertθq_eq_any (P : SGSParams) (G : CellInputs) :
    vertθq P G = geoV P G.Δz * part P G * (G.θz * G.qz) := by
  unfold vertθq
  split_ifs with h
  · simp [geoV, h]
  · rfl

/-- The isentropic invariant when ∇_h q = 0: `I = r² |∇_h θ|²`, `r = q_z / max(θz, θz_min)` (any θz). -/
theorem invIsen_of_flat (P : SGSParams) (G : CellInputs) (hinv : 0 ≤ G.invθ)
    (hq0 : G.invq = 0) (hqθ0 : G.invθq = 0) :
    invIsen P G = (G.qz / max G.θz P.θzMin) ^ 2 * G.invθ := by
  have hE : invArg G.invq G.invθq G.invθ (isentropicSlopeRatio G.qz G.θz P.θzMin)
      = (G.qz / max G.θz P.θzMin) ^ 2 * G.invθ := by
    unfold invArg isentropicSlopeRatio
    rw [hq0, hqθ0]
    ring
  unfold invIsen
  rw [isentropicGradientInvariant_eq, hE]
  exact max_eq_left (mul_nonneg (sq_nonneg _) hinv)

/-- Algebraic core: the sum of the horizontal block and the vertical block. -/
theorem core_flat_sum (cg cΔx cΔz f S X m qz t Δx Δz : ℝ) (hs2 : S ^ 2 = t / m ^ 2)
    (hxz : X * Δz = S * Δx) :
    cg * (cΔx * Δx) ^ 2 * ((qz / m) ^ 2 * t * f) + cg * (cΔz * Δz) ^ 2 * f * (qz * qz)
      = cg * qz ^ 2 * Δz ^ 2 * f * (cΔx ^ 2 * X ^ 2 + cΔz ^ 2) := by
  linear_combination (-(cg * cΔx ^ 2 * f * qz ^ 2 * Δx ^ 2)) * hs2
    - cg * cΔx ^ 2 * f * qz ^ 2 * (X * Δz + S * Δx) * hxz

/-- **(b) Exact sum for general c_Δx, c_Δz** (∇_h q = 0, any ∇θ, any θz including the floor branch):
horizontal block + vertical block = `c_g q_z² Δz² f (c_Δx² X² + c_Δz²)`, `X = s Δx_h / Δz` (the x at :255 of the code). -/
theorem no_double_count_general (P : SGSParams) (G : CellInputs) (hcx : P.cΔx ≠ 0)
    (hΔz : 0 < G.Δz) (hinv : 0 ≤ G.invθ) (hq0 : G.invq = 0) (hqθ0 : G.invθq = 0) :
    horizQQ P G + vertQQ P G
      = P.cg * G.qz ^ 2 * G.Δz ^ 2 * part P G
        * (P.cΔx ^ 2 * cellX G.invθ G.θz P.θzMin G.Δx G.Δz ^ 2 + P.cΔz ^ 2) := by
  have hs2 := cellSlope_sq G.invθ G.θz P.θzMin hinv
  have hxz := cellX_mul_Δz G.invθ G.θz P.θzMin G.Δx G.Δz hΔz.ne'
  rw [horizQQ_eq P G hcx, vertQQ_eq_any, invIsen_of_flat P G hinv hq0 hqθ0]
  unfold geoV
  linear_combination core_flat_sum P.cg P.cΔx P.cΔz (part P G) (cellSlope G.invθ G.θz P.θzMin)
    (cellX G.invθ G.θz P.θzMin G.Δx G.Δz) (max G.θz P.θzMin) G.qz G.invθ G.Δx G.Δz hs2 hxz

/-- **(b) No double counting**: for ∇_h q = 0 and c_Δx = c_Δz = c, horizontal block + vertical block = `c_g c² (Δz q_z)²`,
holding for any ∇θ and any θz (including the floor branch with θz < θz_min). -/
theorem no_double_count (P : SGSParams) (G : CellInputs) (hcx : P.cΔx ≠ 0)
    (hc : P.cΔx = P.cΔz) (hΔz : 0 < G.Δz) (hinv : 0 ≤ G.invθ) (hq0 : G.invq = 0)
    (hqθ0 : G.invθq = 0) :
    horizQQ P G + vertQQ P G = P.cg * P.cΔz ^ 2 * (G.Δz * G.qz) ^ 2 := by
  rw [no_double_count_general P G hcx hΔz hinv hq0 hqθ0, part_eq P G hcx, hc]
  linear_combination (P.cg * P.cΔz ^ 2 * G.Δz ^ 2 * G.qz ^ 2)
    * cap_mul_one_add_sq G.invθ G.θz P.θzMin G.Δx G.Δz

/-- **(b) Deviation when c_Δx ≠ c_Δz** (w1 suspect: f uses `Δx_h`, while the horizontal block scale uses `c_Δx Δx_h`):
sum − `c_g c_Δz² (Δz q_z)²` = `c_g q_z² Δz² f X² (c_Δx² − c_Δz²)`. -/
theorem double_count_defect (P : SGSParams) (G : CellInputs) (hcx : P.cΔx ≠ 0)
    (hΔz : 0 < G.Δz) (hinv : 0 ≤ G.invθ) (hq0 : G.invq = 0) (hqθ0 : G.invθq = 0) :
    horizQQ P G + vertQQ P G - P.cg * P.cΔz ^ 2 * (G.Δz * G.qz) ^ 2
      = P.cg * G.qz ^ 2 * G.Δz ^ 2 * isentropicCellCap G.invθ G.θz P.θzMin G.Δx G.Δz
        * cellX G.invθ G.θz P.θzMin G.Δx G.Δz ^ 2 * (P.cΔx ^ 2 - P.cΔz ^ 2) := by
  rw [no_double_count_general P G hcx hΔz hinv hq0 hqθ0, part_eq P G hcx]
  linear_combination (P.cg * G.qz ^ 2 * G.Δz ^ 2 * P.cΔz ^ 2)
    * cap_mul_one_add_sq G.invθ G.θz P.θzMin G.Δx G.Δz

/-- **(b) The deviation is nonzero**: when c_Δx² ≠ c_Δz², the isentropic surface is tilted, and q_z ≠ 0, the "no double counting" identity does not hold. -/
theorem double_count_defect_ne_zero (P : SGSParams) (G : CellInputs) (hcx : P.cΔx ≠ 0)
    (hcg : P.cg ≠ 0) (hmin : 0 < P.θzMin) (hΔx : 0 < G.Δx) (hΔz : 0 < G.Δz)
    (hinv : 0 < G.invθ) (hq0 : G.invq = 0) (hqθ0 : G.invθq = 0) (hqz : G.qz ≠ 0)
    (hneq : P.cΔx ^ 2 ≠ P.cΔz ^ 2) :
    horizQQ P G + vertQQ P G ≠ P.cg * P.cΔz ^ 2 * (G.Δz * G.qz) ^ 2 := by
  intro h
  have hd := double_count_defect P G hcx hΔz hinv.le hq0 hqθ0
  rw [h, sub_self] at hd
  have hX := cellX_pos G.invθ G.θz P.θzMin G.Δx G.Δz hmin hinv hΔx hΔz
  have hf := cap_pos G.invθ G.θz P.θzMin G.Δx G.Δz
  exact (mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero hcg
    (pow_ne_zero 2 hqz)) (pow_ne_zero 2 hΔz.ne')) hf.ne') (pow_ne_zero 2 hX.ne'))
    (sub_ne_zero.mpr hneq)) hd.symm

/-- Algebraic core: the projection part. -/
theorem core_proj (f S X m qz t Δx Δz : ℝ) (hs2 : S ^ 2 = t / m ^ 2) (hxz : X * Δz = S * Δx)
    (hf : f * (1 + X ^ 2) = 1) :
    (qz / m) ^ 2 * t * (f * Δx ^ 2) = (qz * Δz) ^ 2 * (1 - f) := by
  linear_combination (-(qz ^ 2 * f * Δx ^ 2)) * hs2 - qz ^ 2 * f * (X * Δz + S * Δx) * hxz
    + qz ^ 2 * Δz ^ 2 * hf

/-- **(b) Upper bound on the projection part** (w1 main B, second item): `r² |∇_h θ|² ℓ² ≤ (q_z Δz)²`, `ℓ² = fΔx²`; exactly
`= (q_z Δz)² (1 − f)`. -/
theorem projected_part_le (qz invθ θz θzMin Δx Δz : ℝ) (hΔz : 0 < Δz) (hinv : 0 ≤ invθ) :
    isentropicSlopeRatio qz θz θzMin ^ 2 * invθ * (isentropicCellCap invθ θz θzMin Δx Δz * Δx ^ 2)
      ≤ (qz * Δz) ^ 2 := by
  have hs2 := cellSlope_sq invθ θz θzMin hinv
  have hxz := cellX_mul_Δz invθ θz θzMin Δx Δz hΔz.ne'
  have hf := cap_mul_one_add_sq invθ θz θzMin Δx Δz
  have hf0 := cap_pos invθ θz θzMin Δx Δz
  have heq : isentropicSlopeRatio qz θz θzMin ^ 2 * invθ
        * (isentropicCellCap invθ θz θzMin Δx Δz * Δx ^ 2)
      = (qz * Δz) ^ 2 * (1 - isentropicCellCap invθ θz θzMin Δx Δz) := by
    unfold isentropicSlopeRatio
    linear_combination core_proj (isentropicCellCap invθ θz θzMin Δx Δz) (cellSlope invθ θz θzMin)
      (cellX invθ θz θzMin Δx Δz) (max θz θzMin) qz invθ Δx Δz hs2 hxz hf
  rw [heq]
  nlinarith [mul_nonneg (sq_nonneg (qz * Δz)) hf0.le]

/-! ## §6 Theorem (c): variance consistency (I3 of w2)

Assume q′ and θ′ are exactly functionally related (linearized): `inv_q = F′² inv_θ`, `inv_θq = F′ inv_θ` (derived from ∇q = F′∇θ),
`q_z = F′ θ_z`, and θz ≥ θz_min. -/

theorem turbθθ_zero (P : SGSParams) (G : CellInputs) (hturb : P.coeff = 0) : turbθθ P G = 0 := by
  unfold turbθθ covFromGrad
  rw [hturb]
  ring

theorem turbθq_zero (P : SGSParams) (G : CellInputs) (hturb : P.coeff = 0) : turbθq P G = 0 := by
  unfold turbθq
  rw [turbθθ_zero P G hturb, zero_mul, Real.sqrt_zero, mul_zero]

theorem invIsen_zero_of_func (P : SGSParams) (G : CellInputs) (F : ℝ) (hmin : 0 < P.θzMin)
    (hθ : P.θzMin ≤ G.θz) (hq : G.invq = F ^ 2 * G.invθ) (hqθ : G.invθq = F * G.invθ)
    (hqz : G.qz = F * G.θz) : invIsen P G = 0 := by
  unfold invIsen
  rw [hq, hqθ, hqz]
  exact invariant_zero_of_func G.invθ G.θz P.θzMin F hmin hθ

theorem horizQQ_zero_of_func (P : SGSParams) (G : CellInputs) (F : ℝ) (hmin : 0 < P.θzMin)
    (hθ : P.θzMin ≤ G.θz) (hq : G.invq = F ^ 2 * G.invθ) (hqθ : G.invθq = F * G.invθ)
    (hqz : G.qz = F * G.θz) : horizQQ P G = 0 := by
  unfold horizQQ
  split_ifs
  · rfl
  · rw [invIsen_zero_of_func P G F hmin hθ hq hqθ hqz]
    ring

theorem turbQQ_of_func (P : SGSParams) (G : CellInputs) (F : ℝ) (hqz : G.qz = F * G.θz) :
    turbQQ P G = F ^ 2 * turbθθ P G := by
  unfold turbQQ turbθθ covFromGrad
  rw [hqz]
  ring

theorem vertQQ_of_func (P : SGSParams) (G : CellInputs) (F : ℝ) (hqz : G.qz = F * G.θz) :
    vertQQ P G = part P G * F ^ 2 * vertθθ P G := by
  rw [vertQQ_eq_any, vertθθ_eq_any, hqz]
  ring

theorem vertθq_of_func (P : SGSParams) (G : CellInputs) (F : ℝ) (hqz : G.qz = F * G.θz) :
    vertθq P G = part P G * F * vertθθ P G := by
  rw [vertθq_eq_any, vertθθ_eq_any, hqz]
  ring

theorem vertθθ_pos (P : SGSParams) (G : CellInputs) (hcz : P.cΔz ≠ 0) (hcg : 0 < P.cg)
    (hΔz : 0 < G.Δz) (hθz : G.θz ≠ 0) : 0 < vertθθ P G := by
  rw [vertθθ_eq_any]
  unfold geoV
  have h1 : 0 < (P.cΔz * G.Δz) ^ 2 :=
    lt_of_le_of_ne (sq_nonneg _) (pow_ne_zero 2 (mul_ne_zero hcz hΔz.ne')).symm
  have h2 : 0 < G.θz * G.θz := mul_self_pos.mpr hθz
  exact mul_pos (mul_pos hcg h1) h2

/-- **(c) Exact formula (turbulence block kept as a parameter)**:
`q′q′ = F′² θ′θ′ − F′² (1 − f) · c_g (c_Δz Δz)² θ_z²` (f = `part`; when c_Δx ≠ 0, f = `isentropic_cell_cap`). -/
theorem qq_of_func (P : SGSParams) (G : CellInputs) (F : ℝ) (hmin : 0 < P.θzMin)
    (hθ : P.θzMin ≤ G.θz) (hq : G.invq = F ^ 2 * G.invθ) (hqθ : G.invθq = F * G.invθ)
    (hqz : G.qz = F * G.θz) :
    (covTheta P G).qq = F ^ 2 * (covTheta P G).θθ - F ^ 2 * (1 - part P G) * vertθθ P G := by
  rw [covTheta_qq, covTheta_θθ, turbQQ_of_func P G F hqz,
    horizQQ_zero_of_func P G F hmin hθ hq hqθ hqz, vertQQ_of_func P G F hqz]
  ring

/-- **(c) Turbulence block off**: the code gives `q′q′ = f · F′² · θ′θ′`. -/
theorem qq_eq_cap_mul (P : SGSParams) (G : CellInputs) (F : ℝ) (hturb : P.coeff = 0)
    (hmin : 0 < P.θzMin) (hθ : P.θzMin ≤ G.θz) (hq : G.invq = F ^ 2 * G.invθ)
    (hqθ : G.invθq = F * G.invθ) (hqz : G.qz = F * G.θz) :
    (covTheta P G).qq = part P G * F ^ 2 * (covTheta P G).θθ := by
  rw [qq_of_func P G F hmin hθ hq hqθ hqz, covTheta_θθ, turbθθ_zero P G hturb]
  ring

/-- **(c) Strict violation** (turbulence block arbitrary): when F′ ≠ 0, f < 1 (∇θ ≠ 0 and c_Δx, Δx, Δz, c_g, c_Δz nonzero/positive),
`q′q′ < F′² θ′θ′`. -/
theorem qq_lt_F2_θθ (P : SGSParams) (G : CellInputs) (F : ℝ) (hmin : 0 < P.θzMin)
    (hθ : P.θzMin ≤ G.θz) (hq : G.invq = F ^ 2 * G.invθ) (hqθ : G.invθq = F * G.invθ)
    (hqz : G.qz = F * G.θz) (hF : F ≠ 0) (hcx : P.cΔx ≠ 0) (hcz : P.cΔz ≠ 0)
    (hcg : 0 < P.cg) (hΔx : 0 < G.Δx) (hΔz : 0 < G.Δz) (hinv : 0 < G.invθ) :
    (covTheta P G).qq < F ^ 2 * (covTheta P G).θθ := by
  rw [qq_of_func P G F hmin hθ hq hqθ hqz]
  have hθz : G.θz ≠ 0 := (lt_of_lt_of_le hmin hθ).ne'
  have hp := part_lt_one P G hcx hmin hinv hΔx hΔz
  have hF2 : 0 < F ^ 2 := lt_of_le_of_ne (sq_nonneg F) (pow_ne_zero 2 hF).symm
  have hv := vertθθ_pos P G hcz hcg hΔz hθz
  have := mul_pos (mul_pos hF2 (sub_pos.mpr hp)) hv
  linarith

theorem θq_of_func (P : SGSParams) (G : CellInputs) (F : ℝ) (hturb : P.coeff = 0)
    (hqz : G.qz = F * G.θz) : (covTheta P G).θq = part P G * F * (covTheta P G).θθ := by
  rw [covTheta_θq, covTheta_θθ, turbθq_zero P G hturb, turbθθ_zero P G hturb,
    vertθq_of_func P G F hqz]
  ring

/-- Square of the θ-basis correlation coefficient: `(θ′q′)² = f · θ′θ′ · q′q′` (turbulence block off). -/
theorem corr_sq_theta (P : SGSParams) (G : CellInputs) (F : ℝ) (hturb : P.coeff = 0)
    (hmin : 0 < P.θzMin) (hθ : P.θzMin ≤ G.θz) (hq : G.invq = F ^ 2 * G.invθ)
    (hqθ : G.invθq = F * G.invθ) (hqz : G.qz = F * G.θz) :
    (covTheta P G).θq ^ 2 = part P G * ((covTheta P G).θθ * (covTheta P G).qq) := by
  rw [θq_of_func P G F hturb hqz, qq_eq_cap_mul P G F hturb hmin hθ hq hqθ hqz]
  ring

theorem θθ_pos_noturb (P : SGSParams) (G : CellInputs) (hturb : P.coeff = 0)
    (hmin : 0 < P.θzMin) (hθ : P.θzMin ≤ G.θz) (hcz : P.cΔz ≠ 0) (hcg : 0 < P.cg)
    (hΔz : 0 < G.Δz) : 0 < (covTheta P G).θθ := by
  rw [covTheta_θθ, turbθθ_zero P G hturb, zero_add]
  exact vertθθ_pos P G hcz hcg hΔz (lt_of_lt_of_le hmin hθ).ne'

theorem qq_pos_noturb (P : SGSParams) (G : CellInputs) (F : ℝ) (hturb : P.coeff = 0)
    (hmin : 0 < P.θzMin) (hθ : P.θzMin ≤ G.θz) (hq : G.invq = F ^ 2 * G.invθ)
    (hqθ : G.invθq = F * G.invθ) (hqz : G.qz = F * G.θz) (hF : F ≠ 0) (hcz : P.cΔz ≠ 0)
    (hcg : 0 < P.cg) (hΔz : 0 < G.Δz) : 0 < (covTheta P G).qq := by
  rw [qq_eq_cap_mul P G F hturb hmin hθ hq hqθ hqz]
  have hF2 : 0 < F ^ 2 := lt_of_le_of_ne (sq_nonneg F) (pow_ne_zero 2 hF).symm
  exact mul_pos (mul_pos (part_pos P G) hF2) (θθ_pos_noturb P G hturb hmin hθ hcz hcg hΔz)

/-- General: `a² = f·Q`, `Q > 0` ⇒ `|a / √Q| = √f`. -/
theorem abs_div_sqrt_eq (a Q f : ℝ) (hQ : 0 < Q) (h : a ^ 2 = f * Q) :
    |a / Real.sqrt Q| = Real.sqrt f := by
  have hsQ : 0 < Real.sqrt Q := Real.sqrt_pos.mpr hQ
  rw [abs_div, abs_of_pos hsQ, div_eq_iff hsQ.ne', ← Real.sqrt_mul' f hQ.le, ← h,
    Real.sqrt_sq_eq_abs]

theorem sqrt_lt_one_of_lt (f : ℝ) (h0 : 0 ≤ f) (h1 : f < 1) : Real.sqrt f < 1 := by
  have := Real.sqrt_lt_sqrt h0 h1
  rwa [Real.sqrt_one] at this

theorem sqrt_le_one_of_le (f : ℝ) (h1 : f ≤ 1) : Real.sqrt f ≤ 1 := by
  have := Real.sqrt_le_sqrt h1
  rwa [Real.sqrt_one] at this

/-- **(c) θ-basis diagnosed correlation coefficient** `|θ′q′ / √(θ′θ′ q′q′)| = √f` (turbulence block off). -/
theorem corr_theta_abs (P : SGSParams) (G : CellInputs) (F : ℝ) (hturb : P.coeff = 0)
    (hmin : 0 < P.θzMin) (hθ : P.θzMin ≤ G.θz) (hq : G.invq = F ^ 2 * G.invθ)
    (hqθ : G.invθq = F * G.invθ) (hqz : G.qz = F * G.θz) (hF : F ≠ 0) (hcz : P.cΔz ≠ 0)
    (hcg : 0 < P.cg) (hΔz : 0 < G.Δz) :
    |(covTheta P G).θq / Real.sqrt ((covTheta P G).θθ * (covTheta P G).qq)|
      = Real.sqrt (part P G) :=
  abs_div_sqrt_eq _ _ _
    (mul_pos (θθ_pos_noturb P G hturb hmin hθ hcz hcg hΔz)
      (qq_pos_noturb P G F hturb hmin hθ hq hqθ hqz hF hcz hcg hΔz))
    (corr_sq_theta P G F hturb hmin hθ hq hqθ hqz)

/-- **(c) Main theorem (universal form)**: when q′ = F′θ′ holds exactly, turbulence block off, θz ≥ θz_min, F′ ≠ 0, ∇θ ≠ 0
(other quantities positive), the merged code gives
`q′q′ = f·F′²·θ′θ′`, `q′q′ < F′²·θ′θ′`, and the θ-basis diagnosed correlation coefficient `|ρ| = √f < 1`, whereas a self-consistent distribution should have `|ρ| = 1`. -/
theorem isentropic_variance_inconsistency (P : SGSParams) (G : CellInputs) (F : ℝ)
    (hturb : P.coeff = 0) (hmin : 0 < P.θzMin) (hθ : P.θzMin ≤ G.θz)
    (hq : G.invq = F ^ 2 * G.invθ) (hqθ : G.invθq = F * G.invθ) (hqz : G.qz = F * G.θz)
    (hF : F ≠ 0) (hcx : P.cΔx ≠ 0) (hcz : P.cΔz ≠ 0) (hcg : 0 < P.cg) (hΔx : 0 < G.Δx)
    (hΔz : 0 < G.Δz) (hinv : 0 < G.invθ) :
    (covTheta P G).qq
        = isentropicCellCap G.invθ G.θz P.θzMin G.Δx G.Δz * F ^ 2 * (covTheta P G).θθ
      ∧ (covTheta P G).qq < F ^ 2 * (covTheta P G).θθ
      ∧ |(covTheta P G).θq / Real.sqrt ((covTheta P G).θθ * (covTheta P G).qq)|
          = Real.sqrt (isentropicCellCap G.invθ G.θz P.θzMin G.Δx G.Δz)
      ∧ Real.sqrt (isentropicCellCap G.invθ G.θz P.θzMin G.Δx G.Δz) < 1 := by
  refine ⟨?_, qq_lt_F2_θθ P G F hmin hθ hq hqθ hqz hF hcx hcz hcg hΔx hΔz hinv, ?_, ?_⟩
  · rw [← part_eq P G hcx]
    exact qq_eq_cap_mul P G F hturb hmin hθ hq hqθ hqz
  · rw [← part_eq P G hcx]
    exact corr_theta_abs P G F hturb hmin hθ hq hqθ hqz hF hcz hcg hΔz
  · exact sqrt_lt_one_of_lt _ (cap_pos _ _ _ _ _).le (cap_lt_one _ _ _ _ _ hmin hinv hΔx hΔz)

/-! ### (c) continued: after the θ→T transform, the σ_q upper bound and `sgs_correlation_Tq`, the correlation coefficient actually sampled by the code is still ±√f -/

theorem toT_sq (J f : ℝ) (c : Cov) (h : c.θq ^ 2 = f * (c.θθ * c.qq)) :
    (toT J c).θq ^ 2 = f * ((toT J c).θθ * (toT J c).qq) := by
  rw [toT_θq, toT_θθ, toT_qq]
  linear_combination J ^ 2 * h

theorem applyBound_sq (rmax qtot ε f : ℝ) (c : Cov) (h : c.θq ^ 2 = f * (c.θθ * c.qq))
    (hqq : 0 < c.qq) (hε : ε ^ 2 ≤ c.qq) (hQ : 0 < (rmax * qtot) ^ 2) :
    (applyBound rmax qtot ε c).θq ^ 2
      = f * ((applyBound rmax qtot ε c).θθ * (applyBound rmax qtot ε c).qq) := by
  rw [applyBound_θq, applyBound_θθ, applyBound_qq, max_eq_left hε, mul_pow,
    Real.sq_sqrt (le_min zero_le_one (div_nonneg hQ.le hqq.le))]
  rcases le_total c.qq ((rmax * qtot) ^ 2) with hle | hle
  · have h1 : 1 ≤ (rmax * qtot) ^ 2 / c.qq := (le_div_iff₀ hqq).mpr (by linarith)
    rw [min_eq_left h1, min_eq_left hle, mul_one]
    exact h
  · have h1 : (rmax * qtot) ^ 2 / c.qq ≤ 1 := (div_le_iff₀ hqq).mpr (by linarith)
    rw [min_eq_right h1, min_eq_right hle, h,
      show f * (c.θθ * c.qq) * ((rmax * qtot) ^ 2 / c.qq)
        = f * (c.θθ * (rmax * qtot) ^ 2) * (c.qq / c.qq) by ring,
      div_self hqq.ne', mul_one]

theorem juliaClamp_of_abs_le (x : ℝ) (h : |x| ≤ 1) : juliaClamp x (-1) 1 = x := by
  obtain ⟨h1, h2⟩ := abs_le.mp h
  unfold juliaClamp
  split_ifs with ha hb
  · exfalso
    linarith
  · exfalso
    linarith
  · rfl

theorem sgsCorrelationTq_abs (TT qq Tq fb ε f : ℝ) (hTT : 0 < TT) (hqq : 0 < qq)
    (hv : ε ^ 2 < TT * qq) (hsq : Tq ^ 2 = f * (TT * qq)) (hf1 : f ≤ 1) :
    |sgsCorrelationTq TT qq Tq fb ε| = Real.sqrt f := by
  have hV : corrV TT qq = TT * qq := by
    unfold corrV
    rw [max_eq_left hTT.le, max_eq_left hqq.le]
  have hv' : corrV TT qq > ε ^ 2 := by
    rw [hV]
    exact hv
  have habs := abs_div_sqrt_eq Tq (TT * qq) f (mul_pos hTT hqq) hsq
  have hcl : juliaClamp (Tq / Real.sqrt (TT * qq)) (-1) 1 = Tq / Real.sqrt (TT * qq) :=
    juliaClamp_of_abs_le _ (by rw [habs]; exact sqrt_le_one_of_le f hf1)
  unfold sgsCorrelationTq
  rw [if_pos hv', hV, max_eq_left hv.le, hcl]
  exact habs

/-- **(c) Correlation coefficient actually sampled by the code** (output of `set_tq_correlation!`): under the assumptions of (c), with ∂T/∂θ ≠ 0, the σ_q upper bound
`(r_max q_tot)² > 0`, and the variance above the `eps(FT)` threshold, `|ᶜcorr_Tq| = √f`. -/
theorem corrFinal_abs (P : SGSParams) (G : CellInputs) (F J rmax qtot ε fb : ℝ)
    (hturb : P.coeff = 0) (hmin : 0 < P.θzMin) (hθ : P.θzMin ≤ G.θz)
    (hq : G.invq = F ^ 2 * G.invθ) (hqθ : G.invθq = F * G.invθ) (hqz : G.qz = F * G.θz)
    (hF : F ≠ 0) (hcz : P.cΔz ≠ 0) (hcg : 0 < P.cg) (hΔz : 0 < G.Δz) (hJ : J ≠ 0)
    (hQ : 0 < (rmax * qtot) ^ 2) (hε : ε ^ 2 ≤ (covTheta P G).qq)
    (hv : ε ^ 2 < (finalCov P G J rmax qtot ε).θθ * (finalCov P G J rmax qtot ε).qq) :
    |corrFinal P G J rmax qtot ε fb| = Real.sqrt (part P G) := by
  have hθθ := θθ_pos_noturb P G hturb hmin hθ hcz hcg hΔz
  have hqq := qq_pos_noturb P G F hturb hmin hθ hq hqθ hqz hF hcz hcg hΔz
  have hsqθ := corr_sq_theta P G F hturb hmin hθ hq hqθ hqz
  have hsqT := toT_sq J (part P G) (covTheta P G) hsqθ
  have hsqB : (finalCov P G J rmax qtot ε).θq ^ 2
      = part P G * ((finalCov P G J rmax qtot ε).θθ * (finalCov P G J rmax qtot ε).qq) :=
    applyBound_sq rmax qtot ε (part P G) (toT J (covTheta P G)) hsqT hqq hε hQ
  have hJ2 : 0 < J ^ 2 := lt_of_le_of_ne (sq_nonneg J) (pow_ne_zero 2 hJ).symm
  have hTT : 0 < (finalCov P G J rmax qtot ε).θθ := mul_pos hJ2 hθθ
  have hQQ : 0 < (finalCov P G J rmax qtot ε).qq := lt_min hqq hQ
  unfold corrFinal
  exact sgsCorrelationTq_abs _ _ _ fb ε (part P G) hTT hQQ hv hsqB (part_le_one P G)

/-- **(c) The correlation coefficient actually sampled by the code is strictly less than 1** (adding conditions such as ∇θ ≠ 0, f < 1). -/
theorem corrFinal_lt_one (P : SGSParams) (G : CellInputs) (F J rmax qtot ε fb : ℝ)
    (hturb : P.coeff = 0) (hmin : 0 < P.θzMin) (hθ : P.θzMin ≤ G.θz)
    (hq : G.invq = F ^ 2 * G.invθ) (hqθ : G.invθq = F * G.invθ) (hqz : G.qz = F * G.θz)
    (hF : F ≠ 0) (hcz : P.cΔz ≠ 0) (hcg : 0 < P.cg) (hΔz : 0 < G.Δz) (hJ : J ≠ 0)
    (hQ : 0 < (rmax * qtot) ^ 2) (hε : ε ^ 2 ≤ (covTheta P G).qq)
    (hv : ε ^ 2 < (finalCov P G J rmax qtot ε).θθ * (finalCov P G J rmax qtot ε).qq)
    (hcx : P.cΔx ≠ 0) (hΔx : 0 < G.Δx) (hinv : 0 < G.invθ) :
    |corrFinal P G J rmax qtot ε fb| < 1 := by
  rw [corrFinal_abs P G F J rmax qtot ε fb hturb hmin hθ hq hqθ hqz hF hcz hcg hΔz hJ hQ hε hv]
  exact sqrt_lt_one_of_lt _ (part_pos P G).le (part_lt_one P G hcx hmin hinv hΔx hΔz)

/-- The assumptions of (c) are satisfiable (nonvacuous): at a concrete grid point `q′q′ < F′² θ′θ′`. -/
theorem isentropic_variance_inconsistency_witness :
    ∃ (P : SGSParams) (G : CellInputs) (F : ℝ),
      (covTheta P G).qq < F ^ 2 * (covTheta P G).θθ := by
  refine ⟨⟨1 / 12, 1, 1, 0, 0, 1 / 1000⟩, ⟨1, 1, 1, 1, 1, 1, 1, 1⟩, 1, ?_⟩
  exact qq_lt_F2_θθ _ _ 1 (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by norm_num) one_ne_zero (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by norm_num) (by norm_num)

end

end PhysAudit.ClimaAtmos4822
