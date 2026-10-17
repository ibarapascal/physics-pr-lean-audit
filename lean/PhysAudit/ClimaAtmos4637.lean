-- ClimaAtmos4637: CliMA/ClimaAtmos.jl PR #4637 "Fixes problems near upper boundary and tracer-mass continuity"
-- (merge commit a8efce5a0040924094a58739db91a0f3162f779c) Lean certificate.
-- Transcribes the top-boundary filter `top_velocity`/`set_velocity_at_top!`, the composition of the contravariant velocity, the single-column finite-volume vertical transport of mass and tracers
-- (`ᶜadvdivᵥ`/`ᶠinterp`/`vertical_transport`), the EDMFX difference-form SGS fluxes and the env ratio, and the q_tot SGS term newly added to the continuity equation;
-- proves (A) u³=0 at the top face and column mass conservation, (B) tracer-mass consistency for uniform χ, and gives the pre-change counterexamples and the model-boundary identity for "SGS water flux counted into ρ".
import Mathlib

/-!
# ClimaAtmos.jl #4637 — top boundary & tracer-mass consistency

Transcription conventions
* Real-number semantics: Julia `if`/`?:` branches are transcribed as written; floating-point rounding is out of scope for this file (the replay script
  `replay/scripts/ClimaAtmos4637_replay.jl` handles the comparison).
* Single-column discretization: cells `i ∈ {0,…,n-1}`, faces `k ∈ {0,…,n}`, and cell `i` lies between face `i` and face `i+1`
  (corresponding to ClimaCore's center `i+1` and face `k+½`). Fields are `ℕ → ℝ`, meaningful only on the indices above.
* ClimaCore operators (not part of this PR; transcribed as semantic hypotheses following the documented formulas in ClimaCore `src/Operators/finitedifference.jl`):
  `DivergenceF2C`: `D(v)[i] = (Jv³[i+½] − Jv³[i−½]) / J[i]`, and `SetValue(0)` sets v on the boundary faces to 0;
  `InterpolateC2F`: `I(x)[i] = ½(x[i+½] + x[i−½])`, and `Extrapolate` takes the nearest interior point at the boundary faces;
  `UpwindBiasedProductC2F`: `upwind_biased_product(v, a⁻, a⁺) = ((v + |v|) a⁻ + (v − |v|) a⁺) / 2`,
  with boundary ghost points filled by Extrapolate.
* The two reconstructions `:third_order` and `:vanleer_limiter` of `vertical_transport` are not transcribed; the theorems involving them state the required properties
  (exact for constants `ConstExact`, zero-preserving `ZeroPres`) as explicit hypotheses, while `:none` and `:first_order` are proved.
* Not formalized in this file: the implicit solver (Newton iteration and the hand-written sparse Jacobian `manual_sparse_jacobian.jl`), the sponge layer,
  the horizontal spectral-element operators and DSS, thermodynamics (`ᶜρ⁰` given by `TD.air_density` enters as an arbitrary real input),
  and `sgs_weight_function` (taken as an arbitrary function `W : ℝ → ℝ`).
-/

namespace PhysAudit.ClimaAtmos4637

open Finset

noncomputable section

/-! ## §A Top boundary: contravariant vertical velocity -/

/-- Composition of the contravariant vertical velocity.
```
precomputed_quantities.jl:438      @. ᶠu³ = ᶠuₕ³ + CT3(ᶠu₃)
```
`CT3(C3(u₃))` is the contravariant third component `g³³ u₃` of "a vector with only a covariant third component" (ClimaCore projection, taken as a semantic hypothesis).
`uh3` is `ᶠuₕ³` on the faces (`compute_ᶠuₕ³`, :366-369, obtained by weighted interpolation of center values), here an arbitrary real input. -/
def contravariantU3 (uh3 g33 u3cov : ℝ) : ℝ := uh3 + g33 * u3cov

/-- Full-metric form `u³ = g³¹u₁ + g³²u₂ + g³³u₃` (`uₕ³ = g³¹u₁ + g³²u₂` is `CT3(C12(uₕ))`). -/
def contravariantU3Metric (g31 g32 g33 u1 u2 u3 : ℝ) : ℝ := g31 * u1 + g32 * u2 + g33 * u3

/-- ```
precomputed_quantities.jl:398  function top_velocity(ᶠu₃, ᶠuₕ³)
precomputed_quantities.jl:399      top_level = Spaces.nlevels(axes(ᶠu₃)) - half
precomputed_quantities.jl:400      top_u₃ = Fields.level(ᶠu₃.components.data.:1, top_level)
precomputed_quantities.jl:401      top_uₕ³ = Fields.level(ᶠuₕ³.components.data.:1, top_level)
precomputed_quantities.jl:402      top_g³³ = g³³_field(axes(top_u₃))
precomputed_quantities.jl:403      return @. lazy(-top_uₕ³ / top_g³³) # u³ = uₕ³ + w³ = uₕ³ + w₃ * g³³
precomputed_quantities.jl:404  end
``` -/
def topVelocity (uh3 g33 : ℝ) : ℝ := -uh3 / g33

/-- ```
precomputed_quantities.jl:391  function surface_velocity(ᶠu₃, ᶠuₕ³)
precomputed_quantities.jl:395      return @. lazy(-sfc_uₕ³ / sfc_g³³) # u³ = uₕ³ + w³ = uₕ³ + w₃ * g³³
``` -/
def surfaceVelocity (uh3 g33 : ℝ) : ℝ := -uh3 / g33

/-- `set_velocity_at_top!` before the change (parent commit):
```
-    @. top_u₃ = 0
``` -/
def topVelocityOld (_uh3 _g33 : ℝ) : ℝ := 0

/-- (certificate A1) After the new filter, the contravariant vertical velocity at the model top is 0 (any `uₕ³`, any `g³³ ≠ 0`). -/
theorem top_u3_zero (uh3 g33 : ℝ) (hg : g33 ≠ 0) :
    contravariantU3 uh3 g33 (topVelocity uh3 g33) = 0 := by
  unfold contravariantU3 topVelocity
  first
    | (field_simp; done)
    | (field_simp; ring1)
    | (rw [mul_div_assoc', mul_div_cancel_left₀ _ hg]; ring1)

/-- (certificate A1′) Full-metric form: for any metric (including terrain slope, `g³¹, g³² ≠ 0`) and horizontal wind,
`u₃ = −(g³¹u₁ + g³²u₂)/g³³` gives `u³ = g³¹u₁ + g³²u₂ + g³³u₃ = 0`. -/
theorem top_u3_zero_metric (g31 g32 g33 u1 u2 : ℝ) (hg : g33 ≠ 0) :
    contravariantU3Metric g31 g32 g33 u1 u2 (topVelocity (g31 * u1 + g32 * u2) g33) = 0 := by
  have h := top_u3_zero (g31 * u1 + g32 * u2) g33 hg
  unfold contravariantU3 at h
  unfold contravariantU3Metric
  linarith

/-- The same holds at the ground (`set_velocity_at_surface!`, unchanged by this PR, for comparison). -/
theorem surface_u3_zero (uh3 g33 : ℝ) (hg : g33 ≠ 0) :
    contravariantU3 uh3 g33 (surfaceVelocity uh3 g33) = 0 :=
  top_u3_zero uh3 g33 hg

/-- (A2 counterexample basis) Before the change, the contravariant velocity at the top face equals `uₕ³`. -/
theorem top_u3_old (uh3 g33 : ℝ) : contravariantU3 uh3 g33 (topVelocityOld uh3 g33) = uh3 := by
  unfold contravariantU3 topVelocityOld; ring

/-- (A2 counterexample) There exist a sloped metric (`g³¹ = 1/10`) and horizontal wind (`u₁ = 10`) such that the pre-change top-face `u³ ≠ 0`. -/
theorem old_top_u3_nonzero_witness :
    ∃ g31 g32 g33 u1 u2 : ℝ, g33 ≠ 0 ∧
      contravariantU3Metric g31 g32 g33 u1 u2 (topVelocityOld (g31 * u1 + g32 * u2) g33) ≠ 0 := by
  refine ⟨1 / 10, 0, 1, 10, 0, one_ne_zero, ?_⟩
  unfold contravariantU3Metric topVelocityOld
  norm_num

/-- (certificate A3) EDMFX subdomains: `set_velocity_at_top!` sets `top_u₃ʲ = top_u₃` (:425-431), and the subdomains and the grid mean share `uₕ³`
(`set_velocity_quantities!`, :436-438), so the top-face SGS velocity difference `u³ʲ − u³ = g³³(u₃ʲ − u₃)` is 0. -/
theorem top_sgs_velocity_diff_zero (uh3 g33 u3 u3j : ℝ) (hj : u3j = u3) :
    contravariantU3 uh3 g33 u3j - contravariantU3 uh3 g33 u3 = 0 := by
  subst hj; ring

/-- The contravariant velocity difference between a subdomain and the grid mean comes only from the covariant component difference: `u³ʲ − u³ = g³³(u₃ʲ − u₃)`. -/
theorem sgs_contravariant_diff (uh3 g33 a b : ℝ) :
    contravariantU3 uh3 g33 a - contravariantU3 uh3 g33 b = g33 * (a - b) := by
  unfold contravariantU3; ring

/-! ## §B Single-column finite volume: mass conservation -/

/-- `SetValue(0)` boundary: the face values on the boundary faces (`k = 0` or `k = n`) are set to 0.
```
abbreviations.jl:35  const ᶜadvdivᵥ = Operators.DivergenceF2C(
abbreviations.jl:36      bottom = Operators.SetValue(CT3(0)),
abbreviations.jl:37      top = Operators.SetValue(CT3(0)),
abbreviations.jl:38  )
``` -/
def zeroBC (n : ℕ) (F : ℕ → ℝ) (k : ℕ) : ℝ := if k = 0 ∨ k = n then 0 else F k

/-- ClimaCore `DivergenceF2C` (no boundary condition, i.e. `abbreviations.jl:30  const ᶜdivᵥ = Operators.DivergenceF2C()`):
`D(v)[i] = (Jv³[i+½] − Jv³[i−½]) / J[i]`. `Jc` and `Jf` are the metric Jacobians `ᶜJ` and `ᶠJ` at centers and faces. -/
def divF2C (Jc Jf F : ℕ → ℝ) (i : ℕ) : ℝ := (Jf (i + 1) * F (i + 1) - Jf i * F i) / Jc i

/-- `ᶜadvdivᵥ`: divergence with the boundary-face fluxes set to zero. -/
def advdiv (n : ℕ) (Jc Jf F : ℕ → ℝ) (i : ℕ) : ℝ := divF2C Jc Jf (zeroBC n F) i

/-- ```
abbreviations.jl:57  const ᶠinterp = Operators.InterpolateC2F(
abbreviations.jl:58      bottom = Operators.Extrapolate(),
abbreviations.jl:59      top = Operators.Extrapolate(),
abbreviations.jl:60  )
``` -/
def interpC2F (n : ℕ) (x : ℕ → ℝ) (k : ℕ) : ℝ :=
  if k = 0 then x 0 else if k = n then x (n - 1) else (x (k - 1) + x k) / 2

/-- Mass flux on the faces, `ᶠinterp(Y.c.ρ * ᶜJ) / ᶠJ * ᶠu³`. -/
def massFlux (n : ℕ) (ρ Jc Jf u3 : ℕ → ℝ) (k : ℕ) : ℝ :=
  interpC2F n (fun j => ρ j * Jc j) k / Jf k * u3 k

/-- The merged continuity equation (vertical implicit part).
```
implicit_tendency.jl:127      @. Yₜ.c.ρ -= ᶜadvdivᵥ(ᶠinterp(Y.c.ρ * ᶜJ) / ᶠJ * ᶠu³)
``` -/
def rhoTendencyNew (n : ℕ) (ρ Jc Jf u3 : ℕ → ℝ) (i : ℕ) : ℝ :=
  -advdiv n Jc Jf (massFlux n ρ Jc Jf u3) i

/-- The continuity equation before the change:
```
-    @. Yₜ.c.ρ -= ᶜdivᵥ(ᶠinterp(Y.c.ρ * ᶜJ) / ᶠJ * ᶠu³)
``` -/
def rhoTendencyOld (n : ℕ) (ρ Jc Jf u3 : ℕ → ℝ) (i : ℕ) : ℝ :=
  -divF2C Jc Jf (massFlux n ρ Jc Jf u3) i

/-- After weighting by cell volume, `divF2C` is the difference of the adjacent face fluxes (term by term). -/
theorem sum_mul_divF2C (n : ℕ) (Jc Jf F : ℕ → ℝ) (hJ : ∀ i < n, Jc i ≠ 0) :
    ∑ i ∈ range n, Jc i * divF2C Jc Jf F i
      = ∑ i ∈ range n, (Jf (i + 1) * F (i + 1) - Jf i * F i) := by
  apply Finset.sum_congr rfl
  intro i hi
  have hi' := hJ i (Finset.mem_range.mp hi)
  unfold divF2C
  first
    | (field_simp; done)
    | (field_simp; ring1)
    | (rw [mul_div_assoc', mul_div_cancel_left₀ _ hi'])

/-- Telescoping sum: the volume integral over the column leaves only the fluxes at the top and bottom boundary faces. -/
theorem column_divF2C_sum (n : ℕ) (Jc Jf F : ℕ → ℝ) (hJ : ∀ i < n, Jc i ≠ 0) :
    ∑ i ∈ range n, Jc i * divF2C Jc Jf F i = Jf n * F n - Jf 0 * F 0 := by
  rw [sum_mul_divF2C n Jc Jf F hJ]
  exact Finset.sum_range_sub (fun k => Jf k * F k) n

/-- (certificate B1) `ᶜadvdivᵥ` gives a zero column integral for any face flux. -/
theorem column_advdiv_conserved (n : ℕ) (Jc Jf F : ℕ → ℝ) (hJ : ∀ i < n, Jc i ≠ 0) :
    ∑ i ∈ range n, Jc i * advdiv n Jc Jf F i = 0 := by
  unfold advdiv
  rw [column_divF2C_sum n Jc Jf (zeroBC n F) hJ]
  simp [zeroBC]

/-- (certificate B1′) The merged continuity equation conserves the column mass `Σᵢ Jᵢ ρᵢ` (vertical transport only, any ρ, u³, metric). -/
theorem new_column_mass_conserved (n : ℕ) (ρ Jc Jf u3 : ℕ → ℝ) (hJ : ∀ i < n, Jc i ≠ 0) :
    ∑ i ∈ range n, Jc i * rhoTendencyNew n ρ Jc Jf u3 i = 0 := by
  have h1 : ∀ i ∈ range n, Jc i * rhoTendencyNew n ρ Jc Jf u3 i
      = -(Jc i * advdiv n Jc Jf (massFlux n ρ Jc Jf u3) i) := by
    intro i _; unfold rhoTendencyNew; ring
  rw [Finset.sum_congr rfl h1, Finset.sum_neg_distrib,
    column_advdiv_conserved n Jc Jf _ hJ, neg_zero]

/-- (B2 counterexample) Before the change: the ground filter makes the bottom-face u³=0 while at the top face it only sets u₃=0 (so u³ = uₕ³), and the continuity equation uses `ᶜdivᵥ` without boundary conditions,
so the rate of change of the column mass is exactly `−ρ_{n-1} J_{n-1} uₕ³_top` (leakage through the top face). -/
theorem old_mass_leak (n : ℕ) (hn : 0 < n) (ρ Jc Jf u3 : ℕ → ℝ) (uh3 g33 : ℝ)
    (hJ : ∀ i < n, Jc i ≠ 0) (hJf : Jf n ≠ 0)
    (hbot : u3 0 = 0) (htop : u3 n = contravariantU3 uh3 g33 (topVelocityOld uh3 g33)) :
    ∑ i ∈ range n, Jc i * rhoTendencyOld n ρ Jc Jf u3 i = -(ρ (n - 1) * Jc (n - 1) * uh3) := by
  have h1 : ∀ i ∈ range n, Jc i * rhoTendencyOld n ρ Jc Jf u3 i
      = -(Jc i * divF2C Jc Jf (massFlux n ρ Jc Jf u3) i) := by
    intro i _; unfold rhoTendencyOld; ring
  rw [Finset.sum_congr rfl h1, Finset.sum_neg_distrib, column_divF2C_sum n Jc Jf _ hJ]
  have hI : interpC2F n (fun j => ρ j * Jc j) n = ρ (n - 1) * Jc (n - 1) := by
    simp [interpC2F, hn.ne']
  simp only [massFlux, hI, htop, hbot, contravariantU3, topVelocityOld, mul_zero, add_zero,
    sub_zero]
  first
    | (field_simp; done)
    | (field_simp; ring1)

/-- (B2′ counterexample) As long as the top-layer density, the Jacobian, and the top-face `uₕ³` are all nonzero, the pre-change column mass is not conserved. -/
theorem old_mass_not_conserved (n : ℕ) (hn : 0 < n) (ρ Jc Jf u3 : ℕ → ℝ) (uh3 g33 : ℝ)
    (hJ : ∀ i < n, Jc i ≠ 0) (hJf : Jf n ≠ 0)
    (hbot : u3 0 = 0) (htop : u3 n = contravariantU3 uh3 g33 (topVelocityOld uh3 g33))
    (hρ : ρ (n - 1) ≠ 0) (huh : uh3 ≠ 0) :
    ∑ i ∈ range n, Jc i * rhoTendencyOld n ρ Jc Jf u3 i ≠ 0 := by
  rw [old_mass_leak n hn ρ Jc Jf u3 uh3 g33 hJ hJf hbot htop]
  exact neg_ne_zero.mpr (mul_ne_zero (mul_ne_zero hρ (hJ _ (by omega))) huh)

/-- When the boundary-face flux is 0, `SetValue(0)` does not change the face values. -/
theorem zeroBC_eq_self (n : ℕ) (F : ℕ → ℝ) (hb : F 0 = 0) (ht : F n = 0) : zeroBC n F = F := by
  funext k
  unfold zeroBC
  by_cases h : k = 0 ∨ k = n
  · rw [if_pos h]
    rcases h with h | h
    · rw [h, hb]
    · rw [h, ht]
  · rw [if_neg h]

/-- (certificate B3) After the new filter (u³ at both top and bottom is set to zero by `top_velocity`/`surface_velocity`), `ᶜdivᵥ` and `ᶜadvdivᵥ`
give the same mass tendency in every cell, i.e. the PR's claim "the residual, boundary conditions, and the ρ row of the Jacobian are consistent" holds at the operator level. -/
theorem filtered_divv_eq_advdiv (n : ℕ) (ρ Jc Jf u3 : ℕ → ℝ) (uhb gb uht gt : ℝ)
    (hgb : gb ≠ 0) (hgt : gt ≠ 0)
    (hbot : u3 0 = contravariantU3 uhb gb (surfaceVelocity uhb gb))
    (htop : u3 n = contravariantU3 uht gt (topVelocity uht gt)) (i : ℕ) :
    rhoTendencyOld n ρ Jc Jf u3 i = rhoTendencyNew n ρ Jc Jf u3 i := by
  have h0 : u3 0 = 0 := by rw [hbot]; exact surface_u3_zero uhb gb hgb
  have hn : u3 n = 0 := by rw [htop]; exact top_u3_zero uht gt hgt
  have hF0 : massFlux n ρ Jc Jf u3 0 = 0 := by simp [massFlux, h0]
  have hFn : massFlux n ρ Jc Jf u3 n = 0 := by simp [massFlux, hn]
  unfold rhoTendencyOld rhoTendencyNew advdiv
  rw [zeroBC_eq_self n _ hF0 hFn]

/-- (certificate B3′) Corollary: under the new filter, even the pre-change `ᶜdivᵥ` conserves the column mass. -/
theorem filtered_old_operator_conserves (n : ℕ) (ρ Jc Jf u3 : ℕ → ℝ) (uhb gb uht gt : ℝ)
    (hJ : ∀ i < n, Jc i ≠ 0) (hgb : gb ≠ 0) (hgt : gt ≠ 0)
    (hbot : u3 0 = contravariantU3 uhb gb (surfaceVelocity uhb gb))
    (htop : u3 n = contravariantU3 uht gt (topVelocity uht gt)) :
    ∑ i ∈ range n, Jc i * rhoTendencyOld n ρ Jc Jf u3 i = 0 := by
  have h : ∀ i ∈ range n, Jc i * rhoTendencyOld n ρ Jc Jf u3 i
      = Jc i * rhoTendencyNew n ρ Jc Jf u3 i := by
    intro i _
    rw [filtered_divv_eq_advdiv n ρ Jc Jf u3 uhb gb uht gt hgb hgt hbot htop i]
  rw [Finset.sum_congr rfl h]
  exact new_column_mass_conserved n ρ Jc Jf u3 hJ

/-! ## §C Tracer transport and the "χ ≡ c check" -/

/-- Face reconstruction `R v χ k` (the face value of `ᶠu³ × χ̂`). -/
abbrev Recon := (ℕ → ℝ) → (ℕ → ℝ) → ℕ → ℝ

/-- `:none`: `ᶠu³ * ᶠinterp(ᶜχ)`. -/
def recNone (n : ℕ) : Recon := fun v χ k => v k * interpC2F n χ k

/-- ClimaCore `upwind_biased_product(v, a⁻, a⁺) = ((v + abs(v)) * a⁻ + (v - abs(v)) * a⁺) / 2`. -/
def ubp (v am ap : ℝ) : ℝ := ((v + |v|) * am + (v - |v|) * ap) / 2

/-- `abbreviations.jl:83  const ᶠupwind1 = Operators.UpwindBiasedProductC2F()`; the lower cell of face `k` is cell `k-1` and the upper cell is `k`,
boundary ghost points take the nearest interior point by Extrapolate (at `k = 0`, the ℕ truncation `0 - 1 = 0` is exactly Extrapolate). -/
def recUpwind1 (n : ℕ) : Recon := fun v χ k =>
  ubp (v k) (χ (k - 1)) (if k = n then χ (n - 1) else χ k)

/-- The generic `vertical_transport`: `-(ᶜadvdivᵥ(ᶠinterp(ᶜρ * ᶜJ) / ᶠJ * R(ᶠu³, ᶜχ)))`. -/
def vtGen (n : ℕ) (R : Recon) (ρ Jc Jf v χ : ℕ → ℝ) (i : ℕ) : ℝ :=
  -advdiv n Jc Jf (fun k => interpC2F n (fun j => ρ j * Jc j) k / Jf k * R v χ k) i

/-- ```
implicit_tendency.jl:77  function vertical_transport(ᶜρ, ᶠu³, ᶜχ, dt, ::Val{:none})
implicit_tendency.jl:78      ᶜJ = Fields.local_geometry_field(axes(ᶜρ)).J
implicit_tendency.jl:79      ᶠJ = Fields.local_geometry_field(axes(ᶠu³)).J
implicit_tendency.jl:80      return @. lazy(-(ᶜadvdivᵥ(ᶠinterp(ᶜρ * ᶜJ) / ᶠJ * ᶠu³ * ᶠinterp(ᶜχ))))
implicit_tendency.jl:81  end
``` -/
def vtNone (n : ℕ) (ρ Jc Jf u3 χ : ℕ → ℝ) (i : ℕ) : ℝ :=
  -advdiv n Jc Jf (fun k => interpC2F n (fun j => ρ j * Jc j) k / Jf k * u3 k * interpC2F n χ k) i

/-- ```
implicit_tendency.jl:82  function vertical_transport(ᶜρ, ᶠu³, ᶜχ, dt, ::Val{:first_order})
implicit_tendency.jl:85      return @. lazy(-(ᶜadvdivᵥ(ᶠinterp(ᶜρ * ᶜJ) / ᶠJ * ᶠupwind1(ᶠu³, ᶜχ))))
``` -/
def vtUpwind1 (n : ℕ) (ρ Jc Jf u3 χ : ℕ → ℝ) (i : ℕ) : ℝ :=
  -advdiv n Jc Jf (fun k => interpC2F n (fun j => ρ j * Jc j) k / Jf k * recUpwind1 n u3 χ k) i

/-- Associativity rewrite: `:none` is `vtGen` with `recNone`. -/
theorem vtNone_eq_vtGen (n : ℕ) (ρ Jc Jf u3 χ : ℕ → ℝ) (i : ℕ) :
    vtNone n ρ Jc Jf u3 χ i = vtGen n (recNone n) ρ Jc Jf u3 χ i := by
  simp only [vtNone, vtGen, recNone, mul_assoc]

theorem vtUpwind1_eq_vtGen (n : ℕ) (ρ Jc Jf u3 χ : ℕ → ℝ) (i : ℕ) :
    vtUpwind1 n ρ Jc Jf u3 χ i = vtGen n (recUpwind1 n) ρ Jc Jf u3 χ i := rfl

/-- The reconstruction is exact for constants: `R(v, c) = v·c`. -/
def ConstExact (R : Recon) : Prop := ∀ (v : ℕ → ℝ) (c : ℝ) (k : ℕ), R v (fun _ => c) k = v k * c

/-- The reconstruction maps the zero field to zero. -/
def ZeroPres (R : Recon) : Prop := ∀ (v : ℕ → ℝ) (k : ℕ), R v (fun _ => 0) k = 0

/-- The reconstruction is homogeneous in χ (linear reconstruction). -/
def Homog (R : Recon) : Prop :=
  ∀ (v ψ : ℕ → ℝ) (α : ℝ) (k : ℕ), R v (fun j => α * ψ j) k = α * R v ψ k

theorem interpC2F_const (n : ℕ) (c : ℝ) (k : ℕ) : interpC2F n (fun _ => c) k = c := by
  simp only [interpC2F]
  first | (split_ifs <;> ring) | ring

theorem interpC2F_smul (n : ℕ) (ψ : ℕ → ℝ) (α : ℝ) (k : ℕ) :
    interpC2F n (fun j => α * ψ j) k = α * interpC2F n ψ k := by
  simp only [interpC2F]
  first | (split_ifs <;> ring) | ring

theorem recNone_constExact (n : ℕ) : ConstExact (recNone n) := by
  intro v c k
  simp only [recNone, interpC2F_const]

theorem recUpwind1_constExact (n : ℕ) : ConstExact (recUpwind1 n) := by
  intro v c k
  simp only [recUpwind1, ubp]
  first | (split_ifs <;> ring) | ring

theorem recNone_homog (n : ℕ) : Homog (recNone n) := by
  intro v ψ α k
  simp only [recNone, interpC2F_smul]
  ring

theorem recUpwind1_homog (n : ℕ) : Homog (recUpwind1 n) := by
  intro v ψ α k
  simp only [recUpwind1, ubp]
  first | (split_ifs <;> ring) | ring

theorem zeroPres_of_constExact (R : Recon) (h : ConstExact R) : ZeroPres R := by
  intro v k
  rw [h v 0 k, mul_zero]

theorem advdiv_smul (n : ℕ) (Jc Jf F : ℕ → ℝ) (c : ℝ) (i : ℕ) :
    advdiv n Jc Jf (fun k => c * F k) i = c * advdiv n Jc Jf F i := by
  simp only [advdiv, divF2C, zeroBC]
  first | (split_ifs <;> ring) | ring

theorem advdiv_zero (n : ℕ) (Jc Jf : ℕ → ℝ) (i : ℕ) :
    advdiv n Jc Jf (fun _ => 0) i = 0 := by
  simp [advdiv, divF2C, zeroBC]

/-- (certificate C1) The grid-mean transport degenerates exactly to χ times the continuity equation for uniform χ:
`vertical_transport(ρ, u³, c) = c · (−ᶜadvdivᵥ(ᶠinterp(ρJ)/ᶠJ · u³))` (any ρ, u³, metric; the reconstruction only needs to be exact for constants). -/
theorem vtGen_const (n : ℕ) (R : Recon) (hR : ConstExact R) (ρ Jc Jf v : ℕ → ℝ) (c : ℝ) (i : ℕ) :
    vtGen n R ρ Jc Jf v (fun _ => c) i = c * rhoTendencyNew n ρ Jc Jf v i := by
  have hfun : (fun k => interpC2F n (fun j => ρ j * Jc j) k / Jf k * R v (fun _ => c) k)
      = (fun k => c * massFlux n ρ Jc Jf v k) := by
    funext k; rw [hR v c k]; simp only [massFlux]; ring
  unfold vtGen rhoTendencyNew
  rw [hfun, advdiv_smul]
  ring

theorem vtGen_homog (n : ℕ) (R : Recon) (hR : Homog R) (ρ Jc Jf v ψ : ℕ → ℝ) (α : ℝ) (i : ℕ) :
    vtGen n R ρ Jc Jf v (fun j => α * ψ j) i = α * vtGen n R ρ Jc Jf v ψ i := by
  have h2 : (fun k => interpC2F n (fun j => ρ j * Jc j) k / Jf k * R v (fun j => α * ψ j) k)
      = (fun k => α * (interpC2F n (fun j => ρ j * Jc j) k / Jf k * R v ψ k)) := by
    funext k; rw [hR v ψ α k]; ring
  unfold vtGen
  rw [h2, advdiv_smul]
  ring

/-- `variable_manipulations.jl:47  specific(ρχ, ρ) = ρχ / ρ` -/
def specific (ρχ ρ : ℝ) : ℝ := ρχ / ρ

/-- ```
mass_flux_closures.jl:73  function draft_area(ρa, ρ)
mass_flux_closures.jl:74      return ρa / ρ
``` -/
def draftArea (ρa ρ : ℝ) : ℝ := ρa / ρ

/-- The contribution of a single subdomain to the difference-form SGS flux.
```
edmfx_sgs_flux.jl:149                  @. ᶠu³_diff = ᶠu³ʲs.:($$j) - ᶠu³
edmfx_sgs_flux.jl:150                  @. ᶜa_scalar =
edmfx_sgs_flux.jl:151                      (ᶜχʲ - specific(ᶜρχ, Y.c.ρ)) *
edmfx_sgs_flux.jl:152                      draft_area(Y.c.sgsʲs.:($$j).ρa, ᶜρʲs.:($$j))
edmfx_sgs_flux.jl:153                  vtt = vertical_transport(
edmfx_sgs_flux.jl:154                      ᶜρʲs.:($j),
edmfx_sgs_flux.jl:155                      ᶠu³_diff,
edmfx_sgs_flux.jl:156                      ᶜa_scalar,
edmfx_sgs_flux.jl:157                      dt,
edmfx_sgs_flux.jl:158                      edmfx_sgsflux_upwinding,
edmfx_sgs_flux.jl:159                  )
```
(the env term :165-178 has the same form, with `ᶜρ⁰`, `ᶠu³⁰`, `ᶜχ⁰`, `draft_area(ᶜρa⁰, ᶜρ⁰)`; the q_tot terms :105-115 / :120-130 have the same form.) -/
def sgsTerm (n : ℕ) (R : Recon) (Jc Jf ρk vk u3 χk χ ak : ℕ → ℝ) (i : ℕ) : ℝ :=
  vtGen n R ρk Jc Jf (fun k => vk k - u3 k) (fun j => (χk j - χ j) * ak j) i

/-- (certificate C2) When the subdomain values equal the grid-mean value pointwise, the difference-form SGS flux is zero (any subdomain velocity, density, area; the reconstruction only needs to be zero-preserving). -/
theorem sgsTerm_zero_of_eq {n : ℕ} {R : Recon} (hR : ZeroPres R) {Jc Jf ρk vk u3 χk χ ak : ℕ → ℝ}
    {i : ℕ} (h : ∀ j, χk j = χ j) : sgsTerm n R Jc Jf ρk vk u3 χk χ ak i = 0 := by
  unfold sgsTerm vtGen
  have h0 : (fun j => (χk j - χ j) * ak j) = fun _ => (0 : ℝ) := by
    funext j; rw [h j]; ring
  have hR' : ∀ (v : ℕ → ℝ) (k : ℕ), R v (fun _ => 0) k = 0 := hR
  rw [h0]
  simp only [hR', mul_zero, advdiv_zero, neg_zero]

/-- (certificate C3) Translation invariance: χ ↦ χ + β (for the grid mean and the subdomains simultaneously) does not change the difference-form SGS flux (any reconstruction). -/
theorem sgsTerm_shift (n : ℕ) (R : Recon) (Jc Jf ρk vk u3 χk χ ak : ℕ → ℝ) (β : ℝ) (i : ℕ) :
    sgsTerm n R Jc Jf ρk vk u3 (fun j => χk j + β) (fun j => χ j + β) ak i
      = sgsTerm n R Jc Jf ρk vk u3 χk χ ak i := by
  simp only [sgsTerm]
  have h : (fun j => (χk j + β - (χ j + β)) * ak j) = (fun j => (χk j - χ j) * ak j) := by
    funext j; ring
  rw [h]

/-- (certificate C4) Scaling equivariance: under linear reconstruction (`:none`, `:first_order`), χ ↦ αχ multiplies the flux by α; together with C3 this gives affine equivariance in χ,
`F(αχ + β) = α F(χ)`. -/
theorem sgsTerm_homog (n : ℕ) (R : Recon) (hR : Homog R) (Jc Jf ρk vk u3 χk χ ak : ℕ → ℝ) (α : ℝ)
    (i : ℕ) :
    sgsTerm n R Jc Jf ρk vk u3 (fun j => α * χk j) (fun j => α * χ j) ak i
      = α * sgsTerm n R Jc Jf ρk vk u3 χk χ ak i := by
  simp only [sgsTerm]
  have h : (fun j => (α * χk j - α * χ j) * ak j) = (fun j => α * ((χk j - χ j) * ak j)) := by
    funext j; ring
  rw [h]
  exact vtGen_homog n R hR ρk Jc Jf (fun k => vk k - u3 k) (fun j => (χk j - χ j) * ak j) α i

theorem sgsTerm_affine (n : ℕ) (R : Recon) (hR : Homog R) (Jc Jf ρk vk u3 χk χ ak : ℕ → ℝ)
    (α β : ℝ) (i : ℕ) :
    sgsTerm n R Jc Jf ρk vk u3 (fun j => α * χk j + β) (fun j => α * χ j + β) ak i
      = α * sgsTerm n R Jc Jf ρk vk u3 χk χ ak i := by
  simp only [sgsTerm]
  have h : (fun j => (α * χk j + β - (α * χ j + β)) * ak j)
      = (fun j => α * ((χk j - χ j) * ak j)) := by
    funext j; ring
  rw [h]
  exact vtGen_homog n R hR ρk Jc Jf (fun k => vk k - u3 k) (fun j => (χk j - χ j) * ak j) α i

/-- Regularized ratio (`sgs_weight_function` is taken as an arbitrary function `W`, machine epsilon as the parameter `eps`).
```
variable_manipulations.jl:49  function specific(ρaχ, ρa, ρχ, ρ, turbconv_model)
variable_manipulations.jl:52      weight = sgs_weight_function(ρa / ρ, turbconv_model.a_half)
variable_manipulations.jl:56      return ρa < eps(typeof(ρ)) ? ρχ / ρ : weight * ρaχ / ρa + (1 - weight) * ρχ / ρ
``` -/
def specificReg (eps : ℝ) (W : ℝ → ℝ) (ρaχ ρa ρχ ρ : ℝ) : ℝ :=
  if ρa < eps then ρχ / ρ else W (ρa / ρ) * ρaχ / ρa + (1 - W (ρa / ρ)) * ρχ / ρ

/-- env ratio:
```
variable_manipulations.jl:296-297  ᶜenv_value(grid_scale_value, f_draft, gs) = grid_scale_value - draft_sum(f_draft, gs)
variable_manipulations.jl:341-347  ᶜρaχ⁰ = ᶜenv_value(ᶜρχ, sgsʲ -> sgsʲ.ρa * sgsʲ.χ, Y.c.sgsʲs)
variable_manipulations.jl:349      ᶜρa⁰ = @. lazy(ρa⁰(Y.c.ρ, Y.c.sgsʲs, turbconv_model))   # ρ - Σ ρaʲ (:435-441)
variable_manipulations.jl:355-361  return @. lazy(specific(ᶜρaχ⁰, ᶜρa⁰, ᶜρχ, Y.c.ρ, turbconv_model))
``` -/
def envSpecific (eps : ℝ) (W : ℝ → ℝ) {N : ℕ} (ρχ ρ : ℝ) (ρa χj : Fin N → ℝ) : ℝ :=
  specificReg eps W (ρχ - ∑ d, ρa d * χj d) (ρ - ∑ d, ρa d) ρχ ρ

theorem specific_uniform {ρ : ℝ} (c : ℝ) (hρ : ρ ≠ 0) : specific (c * ρ) ρ = c := by
  rw [specific, mul_div_assoc, div_self hρ, mul_one]

/-- (certificate C5) The env ratio of a uniform tracer is still the same constant, for both branches of the regularization, any weight function W, and any subdomain area
(including the non-physical state ρa⁰ ≤ 0, which takes the `ρa < eps` branch). -/
theorem envSpecific_uniform (eps : ℝ) (W : ℝ → ℝ) {N : ℕ} (ρχ ρ c : ℝ) (ρa χj : Fin N → ℝ)
    (hρ : ρ ≠ 0) (heps : 0 < eps) (hρχ : ρχ = c * ρ) (hχ : ∀ d, χj d = c) :
    envSpecific eps W ρχ ρ ρa χj = c := by
  have hs : ∑ d, ρa d * χj d = (∑ d, ρa d) * c := by
    rw [Finset.sum_mul]
    exact Finset.sum_congr rfl (fun d _ => by rw [hχ d])
  unfold envSpecific specificReg
  rw [hs, hρχ]
  by_cases h : ρ - ∑ d, ρa d < eps
  · rw [if_pos h, mul_div_assoc, div_self hρ, mul_one]
  · rw [if_neg h]
    have hA : ρ - ∑ d, ρa d ≠ 0 := ne_of_gt (lt_of_lt_of_le heps (not_lt.mp h))
    have e1 : c * ρ - (∑ d, ρa d) * c = c * (ρ - ∑ d, ρa d) := by ring
    rw [e1]
    first
      | (field_simp; done)
      | (field_simp; ring1)
      | (rw [mul_div_assoc, mul_div_assoc, mul_div_cancel_right₀ _ hA,
          mul_div_cancel_right₀ _ hρ]; ring1)

/-- Tendency contribution of all SGS difference fluxes (`n` drafts + env) to a given tracer. -/
def sgsAll (n : ℕ) (Rs : Recon) (Jc Jf ρ u3 ρχ : ℕ → ℝ) {N : ℕ}
    (ρj vj ρaj χj : Fin N → ℕ → ℝ) (ρ0 v0 : ℕ → ℝ) (eps : ℝ) (W : ℝ → ℝ) (i : ℕ) : ℝ :=
  (∑ d : Fin N, sgsTerm n Rs Jc Jf (ρj d) (vj d) u3 (χj d) (fun j => specific (ρχ j) (ρ j))
      (fun j => draftArea (ρaj d j) (ρj d j)) i)
  + sgsTerm n Rs Jc Jf ρ0 v0 u3
      (fun j => envSpecific eps W (ρχ j) (ρ j) (fun d => ρaj d j) (fun d => χj d j))
      (fun j => specific (ρχ j) (ρ j))
      (fun j => draftArea (ρ j - ∑ d, ρaj d j) (ρ0 j)) i

/-- Total vertical tendency of a tracer ρχ: grid-mean transport (`advection.jl:254-258`, condensates; q_tot has the same form in
`implicit_tendency.jl:135-137`) + all SGS difference fluxes (`edmfx_sgs_flux.jl:105-132`, `:144-181`). -/
def tracerTend (n : ℕ) (R Rs : Recon) (Jc Jf ρ u3 ρχ : ℕ → ℝ) {N : ℕ}
    (ρj vj ρaj χj : Fin N → ℕ → ℝ) (ρ0 v0 : ℕ → ℝ) (eps : ℝ) (W : ℝ → ℝ) (i : ℕ) : ℝ :=
  vtGen n R ρ Jc Jf u3 (fun j => specific (ρχ j) (ρ j)) i
  + sgsAll n Rs Jc Jf ρ u3 ρχ ρj vj ρaj χj ρ0 v0 eps W i

/-- Total vertical tendency of the merged ρ: continuity equation + the SGS flux of q_tot (newly added by the PR).
```
edmfx_sgs_flux.jl:117                  @. Yₜ.c.ρ += vtt  # Effect of SGS water flux on (moist) air mass
edmfx_sgs_flux.jl:132              @. Yₜ.c.ρ += vtt  # Effect of SGS water flux on (moist) air mass
``` -/
def rhoTend (n : ℕ) (Rs : Recon) (Jc Jf ρ u3 ρq : ℕ → ℝ) {N : ℕ}
    (ρj vj ρaj qj : Fin N → ℕ → ℝ) (ρ0 v0 : ℕ → ℝ) (eps : ℝ) (W : ℝ → ℝ) (i : ℕ) : ℝ :=
  rhoTendencyNew n ρ Jc Jf u3 i + sgsAll n Rs Jc Jf ρ u3 ρq ρj vj ρaj qj ρ0 v0 eps W i

/-- The total SGS contribution of a uniform tracer (grid mean and all drafts equal to c) is zero. -/
theorem sgsAll_zero_of_uniform {n : ℕ} {Rs : Recon} (hRz : ZeroPres Rs) {Jc Jf ρ u3 ρχ : ℕ → ℝ}
    {N : ℕ} {ρj vj ρaj χj : Fin N → ℕ → ℝ} {ρ0 v0 : ℕ → ℝ} {eps : ℝ} {W : ℝ → ℝ} {c : ℝ}
    (hρ : ∀ j, ρ j ≠ 0) (heps : 0 < eps) (hχ : ∀ j, ρχ j = c * ρ j) (hχj : ∀ d j, χj d j = c)
    (i : ℕ) : sgsAll n Rs Jc Jf ρ u3 ρχ ρj vj ρaj χj ρ0 v0 eps W i = 0 := by
  have hgm : ∀ j, specific (ρχ j) (ρ j) = c := fun j => by
    rw [hχ j]; exact specific_uniform c (hρ j)
  have hd : ∀ d : Fin N, sgsTerm n Rs Jc Jf (ρj d) (vj d) u3 (χj d)
      (fun j => specific (ρχ j) (ρ j)) (fun j => draftArea (ρaj d j) (ρj d j)) i = 0 := by
    intro d
    apply sgsTerm_zero_of_eq hRz
    intro j
    simp only [hχj d j, hgm j]
  have he : sgsTerm n Rs Jc Jf ρ0 v0 u3
      (fun j => envSpecific eps W (ρχ j) (ρ j) (fun d => ρaj d j) (fun d => χj d j))
      (fun j => specific (ρχ j) (ρ j))
      (fun j => draftArea (ρ j - ∑ d, ρaj d j) (ρ0 j)) i = 0 := by
    apply sgsTerm_zero_of_eq hRz
    intro j
    simp only [envSpecific_uniform eps W (ρχ j) (ρ j) c (fun d => ρaj d j) (fun d => χj d j)
      (hρ j) heps (hχ j) (fun d => hχj d j), hgm j]
  unfold sgsAll
  rw [Finset.sum_eq_zero (fun d _ => hd d), he, add_zero]

/-- (certificate C6, main theorem) Tracer-mass consistency: if a tracer χ and q_tot are both uniform in the grid mean and all drafts
(χ ≡ c, q_tot ≡ c_q), then on a single column `∂(ρχ)/∂t = c · ∂ρ/∂t` holds cell by cell, where ∂ρ/∂t includes the q_tot SGS term newly added by the PR.
It holds for any ρ, u³, subdomain velocity/density/area, metric, and weight function; the grid-mean reconstruction needs to be exact for constants, and the SGS reconstruction needs to be zero-preserving. -/
theorem uniform_tracer_mass_consistent (n : ℕ) (R Rs : Recon) (hRc : ConstExact R)
    (hRz : ZeroPres Rs) (Jc Jf ρ u3 ρχ ρq : ℕ → ℝ) {N : ℕ} (ρj vj ρaj χj qj : Fin N → ℕ → ℝ)
    (ρ0 v0 : ℕ → ℝ) (eps : ℝ) (W : ℝ → ℝ) (c cq : ℝ)
    (hρ : ∀ j, ρ j ≠ 0) (heps : 0 < eps)
    (hχ : ∀ j, ρχ j = c * ρ j) (hχj : ∀ d j, χj d j = c)
    (hq : ∀ j, ρq j = cq * ρ j) (hqj : ∀ d j, qj d j = cq) (i : ℕ) :
    tracerTend n R Rs Jc Jf ρ u3 ρχ ρj vj ρaj χj ρ0 v0 eps W i
      = c * rhoTend n Rs Jc Jf ρ u3 ρq ρj vj ρaj qj ρ0 v0 eps W i := by
  have hgm : (fun j => specific (ρχ j) (ρ j)) = fun _ => c := funext fun j => by
    show specific (ρχ j) (ρ j) = c
    rw [hχ j]; exact specific_uniform c (hρ j)
  unfold tracerTend rhoTend
  rw [sgsAll_zero_of_uniform hRz hρ heps hχ hχj i, sgsAll_zero_of_uniform hRz hρ heps hq hqj i,
    hgm, vtGen_const n R hRc]
  ring

/-- (C7 boundary identity) With χ uniform but q_tot arbitrary, the consistency gap is exactly `−c × (the SGS term of q_tot)`:
`∂(ρχ)/∂t − c ∂ρ/∂t = −c · S_qtot`. That is, "χ ≡ c is preserved" requires the SGS flux of q_tot to be zero as well; otherwise this is the dilution of specific-humidity-type tracers
by the SGS water-vapor flux (for the physical interpretation see the report). -/
theorem uniform_tracer_defect (n : ℕ) (R Rs : Recon) (hRc : ConstExact R)
    (hRz : ZeroPres Rs) (Jc Jf ρ u3 ρχ ρq : ℕ → ℝ) {N : ℕ} (ρj vj ρaj χj qj : Fin N → ℕ → ℝ)
    (ρ0 v0 : ℕ → ℝ) (eps : ℝ) (W : ℝ → ℝ) (c : ℝ)
    (hρ : ∀ j, ρ j ≠ 0) (heps : 0 < eps)
    (hχ : ∀ j, ρχ j = c * ρ j) (hχj : ∀ d j, χj d j = c) (i : ℕ) :
    tracerTend n R Rs Jc Jf ρ u3 ρχ ρj vj ρaj χj ρ0 v0 eps W i
        - c * rhoTend n Rs Jc Jf ρ u3 ρq ρj vj ρaj qj ρ0 v0 eps W i
      = -(c * sgsAll n Rs Jc Jf ρ u3 ρq ρj vj ρaj qj ρ0 v0 eps W i) := by
  have hgm : (fun j => specific (ρχ j) (ρ j)) = fun _ => c := funext fun j => by
    show specific (ρχ j) (ρ j) = c
    rw [hχ j]; exact specific_uniform c (hρ j)
  unfold tracerTend rhoTend
  rw [sgsAll_zero_of_uniform hRz hρ heps hχ hχj i, hgm, vtGen_const n R hRc]
  ring

/-- (certificate C8) A single explicit step: if `∂(ρχ)/∂t = c ∂ρ/∂t` and `ρχ = cρ`, the ratio is still c after the forward-Euler update (denominator nonzero).
(The same holds for each explicit stage in a multistage IMEX; the implicit stages depend on the Jacobian, see the report.) -/
theorem euler_preserves_uniform (ρ ρχ ρt dt c : ℝ) (hχ : ρχ = c * ρ) (hden : ρ + dt * ρt ≠ 0) :
    specific (ρχ + dt * (c * ρt)) (ρ + dt * ρt) = c := by
  unfold specific
  rw [hχ, div_eq_iff hden]
  ring

/-! ## §D Mass-weighted decomposition and the model boundary of "SGS water flux counted into ρ" -/

/-- env vertical velocity (mass-weighted closure + regularization):
```
variable_manipulations.jl:540  u₃⁰(ρaʲs, u₃ʲs, ρ, u₃, turbconv_model) = specific(
variable_manipulations.jl:541      ρ * u₃ - unrolled_dotproduct(ρaʲs, u₃ʲs),
variable_manipulations.jl:542      ρ - reduce(+, ρaʲs),
variable_manipulations.jl:543      ρ * u₃,
variable_manipulations.jl:544      ρ,
variable_manipulations.jl:545      turbconv_model,
variable_manipulations.jl:546  )
``` -/
def u3env (eps : ℝ) (W : ℝ → ℝ) {N : ℕ} (ρ u3 : ℝ) (m u : Fin N → ℝ) : ℝ :=
  specificReg eps W (ρ * u3 - ∑ d, m d * u d) (ρ - ∑ d, m d) (ρ * u3) ρ

/-- (certificate D1) Closure gap: `Σⱼ mⱼ(uⱼ − u) + m⁰(u⁰ − u) = λ · Σⱼ mⱼ(uⱼ − u)`, where `m⁰ = ρ − Σ mⱼ`,
`λ = 1` (`m⁰ < eps` branch, u⁰ taken as the grid mean) or `λ = 1 − W(m⁰/ρ)` (normal branch). Hence only when W = 1
is the grid-mean velocity strictly the mass-weighted (barycentric) velocity. -/
theorem u3env_closure_defect (eps : ℝ) (W : ℝ → ℝ) {N : ℕ} (ρ u3 : ℝ) (m u : Fin N → ℝ)
    (hρ : ρ ≠ 0) (heps : 0 < eps) :
    (∑ d, m d * (u d - u3)) + (ρ - ∑ d, m d) * (u3env eps W ρ u3 m u - u3)
      = (if ρ - ∑ d, m d < eps then 1 else 1 - W ((ρ - ∑ d, m d) / ρ))
          * ∑ d, m d * (u d - u3) := by
  have hS : ∑ d, m d * (u d - u3) = (∑ d, m d * u d) - (∑ d, m d) * u3 := by
    rw [Finset.sum_mul, ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl (fun d _ => by ring)
  rw [hS]
  unfold u3env specificReg
  split_ifs with h
  · have e : ρ * u3 / ρ = u3 := by rw [mul_comm, mul_div_assoc, div_self hρ, mul_one]
    rw [e]; ring
  · have hA : ρ - ∑ d, m d ≠ 0 := ne_of_gt (lt_of_lt_of_le heps (not_lt.mp h))
    first
      | (field_simp; done)
      | (field_simp; ring1)

/-- (certificate D2) Subdomain decomposition identity (face-by-face algebra): total subdomain flux minus the grid-mean flux =
difference-form flux + χ̄·(mass-flux closure gap) + v̄·(mass-weighted gap of χ). -/
theorem sgs_decomposition {m : ℕ} (M v χ : Fin m → ℝ) (vbar χbar : ℝ) :
    (∑ k, M k * v k * χ k) - (∑ k, M k) * vbar * χbar
      = (∑ k, M k * (v k - vbar) * (χ k - χbar)) + χbar * (∑ k, M k * (v k - vbar))
          + vbar * (∑ k, M k * (χ k - χbar)) := by
  have hR : (∑ k, M k * (v k - vbar) * (χ k - χbar)) + χbar * (∑ k, M k * (v k - vbar))
        + vbar * (∑ k, M k * (χ k - χbar))
      = ∑ k, (M k * (v k - vbar) * (χ k - χbar) + χbar * (M k * (v k - vbar))
          + vbar * (M k * (χ k - χbar))) := by
    rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
  have hL : (∑ k, M k * v k * χ k) - (∑ k, M k) * vbar * χbar
      = ∑ k, (M k * v k * χ k - M k * vbar * χbar) := by
    rw [Finset.sum_sub_distrib, Finset.sum_mul, Finset.sum_mul]
  rw [hL, hR]
  exact Finset.sum_congr rfl (fun k _ => by ring)

/-- (certificate D3) The difference-form SGS flux of dry air is exactly the negative of that of water: `F(1 − q) = −F(q)`. -/
theorem dry_flux_eq_neg_water {m : ℕ} (M v q : Fin m → ℝ) (vbar qbar : ℝ) :
    ∑ k, M k * (v k - vbar) * ((1 - q k) - (1 - qbar))
      = -(∑ k, M k * (v k - vbar) * (q k - qbar)) := by
  rw [← Finset.sum_neg_distrib]
  exact Finset.sum_congr rfl (fun k _ => by ring)

/-- (D4 counterexample) Under strictly mass-weighted closure (`Σ M(v − v̄) = 0`, q̄ the mass-weighted mean), the total SGS mass flux of moist air
(= water + dry = 0, by D2 and D3) is zero, whereas the term `F(q_tot)` that the PR adds to ρ can be nonzero, i.e. `Yₜ.c.ρ += vtt`
is equivalent to omitting the SGS flux of dry air `−F(q_tot)`. Two-subdomain instance. -/
theorem rho_sgs_term_nonzero_under_closure :
    ∃ (M v q : Fin 2 → ℝ) (vbar qbar : ℝ),
      (∑ k, M k * (v k - vbar)) = 0 ∧
      qbar * (∑ k, M k) = ∑ k, M k * q k ∧
      (∑ k, M k * (v k - vbar) * (q k - qbar))
        + (∑ k, M k * (v k - vbar) * ((1 - q k) - (1 - qbar))) = 0 ∧
      (∑ k, M k * (v k - vbar) * (q k - qbar)) ≠ 0 := by
  refine ⟨![1, 1], ![1, -1], ![1, 0], 0, 1 / 2, ?_, ?_, ?_, ?_⟩ <;>
    norm_num [Fin.sum_univ_two]

/-! ## section Arithmetic (side note: structure and constant checks, not counted as certificates) -/

section Arithmetic

/-- Side note: `top_velocity` uses the number of face-space levels `nlevels(axes(ᶠu₃)) − half`, and `set_velocity_at_top!` uses the number of center-space levels
`nlevels(axes(Y.c)) + half`; when the number of face levels = the number of center levels + 1 the two are the same level (ClimaCore convention, taken as a hypothesis). -/
theorem top_level_index_consistent (nc : ℕ) : ((nc + 1 : ℕ) : ℚ) - 1 / 2 = (nc : ℚ) + 1 / 2 := by
  push_cast; ring

end Arithmetic

end

end PhysAudit.ClimaAtmos4637
