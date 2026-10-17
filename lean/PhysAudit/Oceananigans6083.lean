-- Oceananigans6083: CliMA/Oceananigans.jl PR #6083 "Include drift velocities in the implicit part of adaptive
-- implicit tracer advection" (merge commit 1734ff9e4509198f945c5e7740a0bccacc69a720) Lean certificate.
-- Transcribes the explicit velocity scaling, implicit velocity, tridiagonal coefficients and explicit flux of AIVA (AdaptiveVerticallyImplicitDiscretization),
-- the new function `tracer_advecting_velocities`, and the GM advective-form slope tapering bundled into the diff (`advective_skew_diffusion.jl`);
-- over all reals, proves: velocities enter the operator only through the total velocity, conservation, slope upper bound, zeroing under unstable stratification, consistency with the diffusive form; and gives counterexamples for the pre-fix code and
-- for configurations with `minimum_bz > 0`.
import Mathlib

/-!
# Oceananigans.jl #6083 — AIVA total velocity + GM slope tapering

Transcription conventions
* Single column (Flat, Flat, Bounded), `density = nothing` (ρᶠ = ρᶜ = 1), `closure = nothing` (`ivd_diagonal = 1`, no implicit
    diffusion term), no implicit-explicit boundary flux (`boundary_flux_diagonal = 0`). The closure diffusion term is independent of velocity and does not affect the conclusions of §2,
    but this file does not transcribe it (speculation: the conclusions are unchanged if it is added).
* Subscripts are 0-based natural numbers: cells k = 0..N-1, faces f = 0..N; face k is the bottom face of cell k and face k+1 is its top face. Julia is 1-based,
    `k(Julia) = k(Lean) + 1`. Fields are `ℕ → ℝ`, and only values within the range are used.
* `!peripheral_node(...)` is a Bool multiplied onto a real (Julia `x * true = x`, `x * false = 0`), transcribed as 0/1-valued masks
    `act` (faces) and `actc` (cell centers).
* The explicit WENO reconstruction values are not transcribed and are taken as the parameter `recon W c f` (the reconstructed value at face f, which may depend on the whole velocity field and tracer field); this layer is an assumption.
* Floating point: `ifelse` → `if`; `max(·, zero)`/`min` copied verbatim. IEEE `x/0 = ±Inf` appears only in the GM part; see the comments in §4.
-/

namespace Oceananigans6083

noncomputable section

open Finset

/-! ## §1 Explicit velocity scaling and implicit velocity (single point) -/

/-- Vertical Courant number.
```
adaptive_implicit_vertical_advection.jl:23      α  = abs(w) * Δt / Δz
implicit_vertical_advection.jl:24               α  = abs(w) * Δt / Δz
``` -/
def courant (Δt Δz w : ℝ) : ℝ := |w| * Δt / Δz

/-- ```
adaptive_implicit_vertical_advection.jl:19  @inline function explicit_velocity_scaleᶜᶜᶠ(i, j, k, grid, scheme, td, W)
adaptive_implicit_vertical_advection.jl:20      Δt = _unwrap_for_gpu(td.Δt)
adaptive_implicit_vertical_advection.jl:21      Δz = Δzᶜᶜᶠ(i, j, k, grid)
adaptive_implicit_vertical_advection.jl:22      w  = @inbounds W[i, j, k]
adaptive_implicit_vertical_advection.jl:23      α  = abs(w) * Δt / Δz
adaptive_implicit_vertical_advection.jl:24      return ifelse(α > td.cfl, td.cfl / α, one(α))
``` -/
def explicitScale (Δt Δz cfl w : ℝ) : ℝ :=
  if courant Δt Δz w > cfl then cfl / courant Δt Δz w else 1

/-- ```
implicit_vertical_advection.jl:20  @inline function implicit_vertical_velocityᶜᶜᶠ(i, j, k, grid, scheme, td, W)
implicit_vertical_advection.jl:21      Δt = _unwrap_for_gpu(td.Δt)
implicit_vertical_advection.jl:22      Δz = previous_Δzᶜᶜᶠ(i, j, k, grid)
implicit_vertical_advection.jl:23      w  = @inbounds W[i, j, k]
implicit_vertical_advection.jl:24      α  = abs(w) * Δt / Δz
implicit_vertical_advection.jl:25      return w * (1 - ifelse(α > td.cfl, td.cfl / α, one(α)))
```
(On a static grid `previous_Δzᶜᶜᶠ = Δzᶜᶜᶠ`, :29; the σ⁻/σⁿ corrections for variable grids are outside the scope of this file.) -/
def implicitVel (Δt Δz cfl w : ℝ) : ℝ :=
  w * (1 - (if courant Δt Δz w > cfl then cfl / courant Δt Δz w else 1))

/-- Explicit-part velocity + implicit-part velocity = the velocity given, for any Courant number. -/
theorem explicit_add_implicit_eq (Δt Δz cfl w : ℝ) :
    explicitScale Δt Δz cfl w * w + implicitVel Δt Δz cfl w = w := by
  unfold explicitScale implicitVel
  ring

lemma courant_smul (Δt Δz s w : ℝ) : courant Δt Δz (s * w) = |s| * courant Δt Δz w := by
  unfold courant
  rw [abs_mul]
  ring

/-- The Courant number of the explicit part does not exceed cfl (any w). -/
theorem explicit_courant_le_cfl (Δt Δz cfl w : ℝ) (hcfl : 0 ≤ cfl) :
    courant Δt Δz (explicitScale Δt Δz cfl w * w) ≤ cfl := by
  unfold explicitScale
  split_ifs with h
  · have hα : 0 < courant Δt Δz w := lt_of_le_of_lt hcfl h
    have hα' : courant Δt Δz w ≠ 0 := hα.ne'
    have hs : 0 ≤ cfl / courant Δt Δz w := div_nonneg hcfl hα.le
    rw [courant_smul, abs_of_nonneg hs]
    have hcancel : cfl / courant Δt Δz w * courant Δt Δz w = cfl := by
      field_simp
    exact le_of_eq hcancel
  · rw [one_mul]
    exact not_lt.mp h

/-- When the Courant number does not exceed cfl, the implicit velocity is zero (fully explicit). -/
theorem implicitVel_eq_zero_of_le (Δt Δz cfl w : ℝ) (h : courant Δt Δz w ≤ cfl) :
    implicitVel Δt Δz cfl w = 0 := by
  unfold implicitVel
  rw [if_neg (show ¬ (courant Δt Δz w > cfl) from not_lt.mpr h)]
  ring

/-- Handing the velocity only to the explicit part while the implicit part receives 0 (how the pre-fix code treated drift velocities; the part masked out on boundary faces after the fix is analogous):
when the Courant number exceeds cfl, the velocity actually applied `s·w` is never equal to `w`. A counterexample on the parameter family. -/
theorem scaled_only_velocity_ne (Δt Δz cfl w : ℝ) (hcfl : 0 ≤ cfl) (h : courant Δt Δz w > cfl)
    (hw : w ≠ 0) : explicitScale Δt Δz cfl w * w + implicitVel Δt Δz cfl 0 ≠ w := by
  have hα : 0 < courant Δt Δz w := lt_of_le_of_lt hcfl h
  have hlt : cfl / courant Δt Δz w < 1 := (div_lt_one hα).mpr h
  have h0 : implicitVel Δt Δz cfl 0 = 0 := by
    unfold implicitVel
    ring
  rw [h0, add_zero]
  unfold explicitScale
  rw [if_pos h]
  intro heq
  have hz : (cfl / courant Δt Δz w - 1) * w = 0 := by linear_combination heq
  rcases mul_eq_zero.mp hz with h1 | h1
  · linarith
  · exact hw h1

/-! ## §2 Single-column tridiagonal coefficients and conservation -/

/-- Single-column geometry and AIVA parameters. `act f` = `!peripheral_node(i, j, f, grid, ℓx, ℓy, Face())`,
`actc k` = `!peripheral_node(i, j, k, grid, ℓx, ℓy, Center())`, taking values 0/1. -/
structure Column where
  N : ℕ
  Δt : ℝ
  cfl : ℝ
  V : ℕ → ℝ
  Az : ℕ → ℝ
  Δzf : ℕ → ℝ
  act : ℕ → ℝ
  actc : ℕ → ℝ

/-- Implicit velocity at face f `implicit_vertical_velocity(ℓx, ℓy, i, j, f, grid, scheme, td, w)` (implicit_vertical_advection.jl:90). -/
def wI (col : Column) (W : ℕ → ℝ) (f : ℕ) : ℝ := implicitVel col.Δt (col.Δzf f) col.cfl (W f)

/-- Upper diagonal (coefficient of q_{k+1} in row k).
```
implicit_vertical_advection.jl:102 @inline function implicit_advection_upper_diagonal(i, j, k, grid, advection::AIVA, w, Δt, ℓx, ℓy, ℓz::Center, density=nothing)
implicit_vertical_advection.jl:105     wⁱ  = implicit_vertical_velocity(ℓx, ℓy, i, j, k+1, grid, scheme, td, w)
implicit_vertical_advection.jl:106     Azᵢ = Az(i, j, k+1, grid, ℓx, ℓy, Face())
implicit_vertical_advection.jl:109     V⁻¹ = 1 / volume(i, j, k, grid, ℓx, ℓy, Center())
implicit_vertical_advection.jl:110     return Δt * V⁻¹ * Azᵢ * ρᶠ / ρᶜ * min(wⁱ, zero(wⁱ)) * !peripheral_node(i, j, k+1, grid, ℓx, ℓy, Face())
```
(ρᶠ = ρᶜ = 1.) -/
def upperAdv (col : Column) (W : ℕ → ℝ) (k : ℕ) : ℝ :=
  col.Δt * (1 / col.V k) * col.Az (k + 1) * min (wI col W (k + 1)) 0 * col.act (k + 1)

/-- Lower diagonal (Tridiagonal convention: the element with index k′ sits in row k′+1, multiplying q_{k′}).
```
implicit_vertical_advection.jl:115 @inline function implicit_advection_lower_diagonal(i, j, k′, grid, advection::AIVA, w, Δt, ℓx, ℓy, ℓz::Center, density=nothing)
implicit_vertical_advection.jl:118     k   = k′ + 1
implicit_vertical_advection.jl:119     wⁱ  = implicit_vertical_velocity(ℓx, ℓy, i, j, k, grid, scheme, td, w)
implicit_vertical_advection.jl:120     Azᵢ = Az(i, j, k, grid, ℓx, ℓy, Face())
implicit_vertical_advection.jl:123     V⁻¹ = 1 / volume(i, j, k, grid, ℓx, ℓy, Center())
implicit_vertical_advection.jl:124     return - Δt * V⁻¹ * Azᵢ * ρᶠ / ρᶜ * max(wⁱ, zero(wⁱ)) * !peripheral_node(i, j, k′, grid, ℓx, ℓy, Center())
``` -/
def lowerAdv (col : Column) (W : ℕ → ℝ) (k' : ℕ) : ℝ :=
  -col.Δt * (1 / col.V (k' + 1)) * col.Az (k' + 1) * max (wI col W (k' + 1)) 0 * col.actc k'

/-- Advective part of the main diagonal.
```
implicit_vertical_advection.jl:127 @inline function implicit_advection_diagonal(i, j, k, grid, advection::AIVA, w, Δt, ℓx, ℓy, ℓz::Center, density=nothing)
implicit_vertical_advection.jl:130     wⁱ⁺ = implicit_vertical_velocity(ℓx, ℓy, i, j, k+1, grid, scheme, td, w)
implicit_vertical_advection.jl:131     wⁱ⁻ = implicit_vertical_velocity(ℓx, ℓy, i, j, k,   grid, scheme, td, w)
implicit_vertical_advection.jl:140     active⁺ = !peripheral_node(i, j, k+1, grid, ℓx, ℓy, Face())
implicit_vertical_advection.jl:141     active⁻ = !peripheral_node(i, j, k,   grid, ℓx, ℓy, Face())
implicit_vertical_advection.jl:143     V⁻¹ = 1 / volume(i, j, k, grid, ℓx, ℓy, Center())
implicit_vertical_advection.jl:145     return Δt * V⁻¹ / ρᶜ * (Az⁺ * ρᶠ⁺ * max(wⁱ⁺, zero(wⁱ⁺)) * active⁺ -
implicit_vertical_advection.jl:146                             Az⁻ * ρᶠ⁻ * min(wⁱ⁻, zero(wⁱ⁻)) * active⁻)
``` -/
def diagAdv (col : Column) (W : ℕ → ℝ) (k : ℕ) : ℝ :=
  col.Δt * (1 / col.V k) *
    (col.Az (k + 1) * max (wI col W (k + 1)) 0 * col.act (k + 1) -
     col.Az k * min (wI col W k) 0 * col.act k)

/-- Main diagonal = `ivd_diagonal` (closure = nothing ⇒ 1) + advective part + boundary term (0).
```
vertically_implicit_diffusion_solver.jl:250     dκ  = ivd_diagonal(i, j, k, grid, clo, K, id, ℓx, ℓy, ℓz, Δt, clk, fields)
vertically_implicit_diffusion_solver.jl:251     dw  = implicit_advection_diagonal(i, j, k, grid, advection, w, Δt, ℓx, ℓy, ℓz, density)
vertically_implicit_diffusion_solver.jl:253     return dκ + dw + dbc
``` -/
def diagCoeff (col : Column) (W : ℕ → ℝ) (k : ℕ) : ℝ := 1 + diagAdv col W k

/-- Row k of the tridiagonal matrix acting on q (row 0 has no lower diagonal, row N-1 has no upper diagonal). -/
def matAct (col : Column) (W : ℕ → ℝ) (q : ℕ → ℝ) (k : ℕ) : ℝ :=
  (if 0 < k then lowerAdv col W (k - 1) * q (k - 1) else 0) + diagCoeff col W k * q k +
    (if k + 1 < col.N then upperAdv col W k * q (k + 1) else 0)

/-- Implicit first-order upwind face flux (as stated in the comment at implicit_vertical_advection.jl:79-80
`F_{k+1} = Az_{k+1} ρᶠ [max(wⁱ, 0) c_k + min(wⁱ, 0) c_{k+1}]`, multiplied by the face mask). -/
def implicitFlux (col : Column) (W : ℕ → ℝ) (q : ℕ → ℝ) (f : ℕ) : ℝ :=
  col.Az f * col.act f * (max (wI col W f) 0 * q (f - 1) + min (wI col W f) 0 * q f)

/-- Every row is a flux difference: `V_k (M q)_k = V_k q_k + Δt (F_{k+1} − F_k)`.
Mask assumption: the bottom and top faces are inactive (`act 0 = act N = 0`, Bounded column), and interior faces agree with the center mask of the cell below
(both are 1 on a grid without immersed boundary). -/
theorem row_identity (col : Column) (W q : ℕ → ℝ) (k : ℕ) (hk : k < col.N) (hV : col.V k ≠ 0)
    (h0 : col.act 0 = 0) (hN : col.act col.N = 0)
    (hmask : ∀ j, 0 < j → j < col.N → col.actc (j - 1) = col.act j) :
    col.V k * matAct col W q k =
      col.V k * q k + col.Δt * (implicitFlux col W q (k + 1) - implicitFlux col W q k) := by
  have hsub : k + 1 - 1 = k := by omega
  unfold matAct diagCoeff diagAdv upperAdv lowerAdv implicitFlux
  rw [hsub]
  by_cases hk0 : 0 < k
  · have hk1 : k - 1 + 1 = k := by omega
    rw [if_pos hk0, hk1, hmask k hk0 hk]
    by_cases hkN : k + 1 < col.N
    · rw [if_pos hkN]
      first
        | (field_simp; ring1)
        | field_simp
        | ring1
    · rw [if_neg hkN]
      have hkN' : k + 1 = col.N := by omega
      rw [hkN', hN]
      first
        | (field_simp; ring1)
        | field_simp
        | ring1
  · have hk00 : k = 0 := by omega
    subst hk00
    rw [if_neg (lt_irrefl (0 : ℕ)), h0]
    by_cases hkN : 0 + 1 < col.N
    · rw [if_pos hkN]
      first
        | (field_simp; ring1)
        | field_simp
        | ring1
    · rw [if_neg hkN]
      have hkN' : 0 + 1 = col.N := by omega
      rw [hkN', hN]
      first
        | (field_simp; ring1)
        | field_simp
        | ring1

/-- Conservation (implicit part): for any velocity field W and any Courant number, `Σ_k V_k (M q)_k = Σ_k V_k q_k`. -/
theorem implicit_conserves (col : Column) (W q : ℕ → ℝ)
    (hV : ∀ k, k < col.N → col.V k ≠ 0) (h0 : col.act 0 = 0) (hN : col.act col.N = 0)
    (hmask : ∀ j, 0 < j → j < col.N → col.actc (j - 1) = col.act j) :
    ∑ k ∈ range col.N, col.V k * matAct col W q k = ∑ k ∈ range col.N, col.V k * q k := by
  rw [Finset.sum_congr rfl (fun k hk => row_identity col W q k (Finset.mem_range.mp hk)
    (hV k (Finset.mem_range.mp hk)) h0 hN hmask)]
  rw [Finset.sum_add_distrib, ← Finset.mul_sum, Finset.sum_range_sub (implicitFlux col W q)]
  have hF0 : implicitFlux col W q 0 = 0 := by
    unfold implicitFlux
    rw [h0]
    ring
  have hFN : implicitFlux col W q col.N = 0 := by
    unfold implicitFlux
    rw [hN]
    ring
  rw [hF0, hFN]
  ring

/-- The implicit solve `M q' = q★` preserves the column total. -/
theorem implicit_solve_conserves (col : Column) (W q' qstar : ℕ → ℝ)
    (hV : ∀ k, k < col.N → col.V k ≠ 0) (h0 : col.act 0 = 0) (hN : col.act col.N = 0)
    (hmask : ∀ j, 0 < j → j < col.N → col.actc (j - 1) = col.act j)
    (hsolve : ∀ k, k < col.N → matAct col W q' k = qstar k) :
    ∑ k ∈ range col.N, col.V k * q' k = ∑ k ∈ range col.N, col.V k * qstar k := by
  rw [← implicit_conserves col W q' hV h0 hN hmask]
  exact Finset.sum_congr rfl (fun k hk => by rw [hsolve k (Finset.mem_range.mp hk)])

/-- Explicit vertical flux: scaling factor × explicit (ExplicitTimeDiscretization) flux `Az · W · reconstructed value`.
```
adaptive_implicit_vertical_advection.jl:59 @inline function advective_tracer_flux_z(i, j, k, grid, scheme, td::AVID, W, c)
adaptive_implicit_vertical_advection.jl:60     s = explicit_velocity_scaleᶜᶜᶠ(i, j, k, grid, scheme, td, W)
adaptive_implicit_vertical_advection.jl:61     return s * advective_tracer_flux_z(i, j, k, grid, scheme, ExplicitTimeDiscretization(), W, c)
```
(The WENO reconstructed value in the explicit flux = parameter `R f`, not transcribed.) -/
def explicitFlux (col : Column) (W R : ℕ → ℝ) (f : ℕ) : ℝ :=
  explicitScale col.Δt (col.Δzf f) col.cfl (W f) * (col.Az f * W f * R f)

/-- Explicit vertical tendency (z part of `div_Uc`, tracer_advection_operators.jl:31-34: `- V⁻¹ δz(F)`). -/
def explicitTend (col : Column) (W R : ℕ → ℝ) (k : ℕ) : ℝ :=
  -((1 / col.V k) * (explicitFlux col W R (k + 1) - explicitFlux col W R k))

/-- Conservation (explicit part): with W = 0 on the bottom and top faces (impermeable; the PR tests are set up this way), `Σ V_k G_k = 0`. -/
theorem explicit_conserves (col : Column) (W R : ℕ → ℝ)
    (hV : ∀ k, k < col.N → col.V k ≠ 0) (hW0 : W 0 = 0) (hWN : W col.N = 0) :
    ∑ k ∈ range col.N, col.V k * explicitTend col W R k = 0 := by
  have hrow : ∀ k ∈ range col.N, col.V k * explicitTend col W R k =
      -(explicitFlux col W R (k + 1) - explicitFlux col W R k) := by
    intro k hk
    have hVk := hV k (Finset.mem_range.mp hk)
    unfold explicitTend
    first
      | (field_simp; ring1)
      | field_simp
      | ring1
  rw [Finset.sum_congr rfl hrow, Finset.sum_neg_distrib, Finset.sum_range_sub (explicitFlux col W R)]
  unfold explicitFlux
  rw [hW0, hWN]
  ring

/-! ## §3 Total velocity: `tracer_advecting_velocities` and a single step -/

/-- New function (added by the PR): transport velocity + biogeochemical drift velocity + closure auxiliary velocity, with advective forcing added on top.
```
hydrostatic_free_surface_tendency_kernel_functions.jl:103 @inline function tracer_advecting_velocities(velocities, biogeochemistry, closure, closure_fields, forcing, val_tracer_name)
hydrostatic_free_surface_tendency_kernel_functions.jl:104     biogeochemical_velocities = biogeochemical_drift_velocity(biogeochemistry, val_tracer_name)
hydrostatic_free_surface_tendency_kernel_functions.jl:105     closure_velocities = closure_auxiliary_velocity(closure, closure_fields, val_tracer_name)
hydrostatic_free_surface_tendency_kernel_functions.jl:106     total_velocities = sum_of_velocities(velocities, biogeochemical_velocities, closure_velocities)
hydrostatic_free_surface_tendency_kernel_functions.jl:107     return with_advective_forcing(forcing, total_velocities)
```
(`sum_of_velocities` drops `nothing` terms = adds 0; `with_advective_forcing` = `sum_of_velocities(velocities(forcing), total)`,
advective_forcing.jl:78-79; `SumOfArrays` adds term by term by index, sum_of_arrays.jl:23-26. Only the w component is taken.) -/
def tracerAdvectingW (wt wb wc wF : ℕ → ℝ) : ℕ → ℝ := fun f => wF f + (wt f + wb f + wc f)

/-- One AIVA step (forward-Euler explicit part + implicit solve): `M(W) c_new = c + Δt G(W)`;
the same W enters both the explicit flux and the implicit coefficients. -/
def aivaStep (col : Column) (recon : (ℕ → ℝ) → (ℕ → ℝ) → ℕ → ℝ) (W c cnew : ℕ → ℝ) : Prop :=
  ∀ k, k < col.N → matAct col W cnew k = c k + col.Δt * explicitTend col W (recon W c) k

/-- After the fix (ab2_step.jl:244-256, rk_step.jl:230-242): `implicit_step!` receives `c_velocities = tracer_advecting_velocities(...)`,
using the same total velocity as the explicit tendency (tendency_kernel_functions.jl:143). -/
def newStep (col : Column) (recon : (ℕ → ℝ) → (ℕ → ℝ) → ℕ → ℝ) (wt wb wc wF c cnew : ℕ → ℝ) : Prop :=
  aivaStep col recon (tracerAdvectingW wt wb wc wF) c cnew

/-- Before the fix (parent commit 1bd69f4e): `implicit_step!(..., c_advection, model.transport_velocities)` — the implicit part receives only the transport velocity,
while the explicit part still uses the total velocity. -/
def oldStep (col : Column) (recon : (ℕ → ℝ) → (ℕ → ℝ) → ℕ → ℝ) (wt wb wc wF c cnew : ℕ → ℝ) : Prop :=
  ∀ k, k < col.N → matAct col wt cnew k =
    c k + col.Δt * explicitTend col (tracerAdvectingW wt wb wc wF)
      (recon (tracerAdvectingW wt wb wc wF) c) k

/-- Certificate A: the post-fix single step depends only on the total velocity. Two velocity decompositions with the same face-by-face total give exactly the same step relation (any Courant number,
any reconstruction `recon`). -/
theorem newStep_depends_only_on_total (col : Column) (recon : (ℕ → ℝ) → (ℕ → ℝ) → ℕ → ℝ)
    (wt wb wc wF wt' wb' wc' wF' c cnew : ℕ → ℝ)
    (h : ∀ f, tracerAdvectingW wt wb wc wF f = tracerAdvectingW wt' wb' wc' wF' f) :
    newStep col recon wt wb wc wF c cnew ↔ newStep col recon wt' wb' wc' wF' c cnew := by
  have hW : tracerAdvectingW wt wb wc wF = tracerAdvectingW wt' wb' wc' wF' := funext h
  unfold newStep
  rw [hW]

/-- Corollary: the same velocity field w, supplied as transport velocity, biogeochemical drift velocity or advective forcing, gives the same step relation. -/
theorem newStep_transport_iff_drift_iff_forcing (col : Column)
    (recon : (ℕ → ℝ) → (ℕ → ℝ) → ℕ → ℝ) (w c cnew : ℕ → ℝ) :
    (newStep col recon w 0 0 0 c cnew ↔ newStep col recon 0 w 0 0 c cnew) ∧
    (newStep col recon w 0 0 0 c cnew ↔ newStep col recon 0 0 0 w c cnew) := by
  constructor
  · apply newStep_depends_only_on_total
    intro f
    simp only [tracerAdvectingW, Pi.zero_apply]
    ring
  · apply newStep_depends_only_on_total
    intro f
    simp only [tracerAdvectingW, Pi.zero_apply]
    ring

/-- Certificate B: the post-fix single step is conservative (for any velocity decomposition and any Courant number). -/
theorem newStep_conserves (col : Column) (recon : (ℕ → ℝ) → (ℕ → ℝ) → ℕ → ℝ)
    (wt wb wc wF c cnew : ℕ → ℝ)
    (hV : ∀ k, k < col.N → col.V k ≠ 0) (h0 : col.act 0 = 0) (hN : col.act col.N = 0)
    (hmask : ∀ j, 0 < j → j < col.N → col.actc (j - 1) = col.act j)
    (hW0 : tracerAdvectingW wt wb wc wF 0 = 0) (hWN : tracerAdvectingW wt wb wc wF col.N = 0)
    (hstep : newStep col recon wt wb wc wF c cnew) :
    ∑ k ∈ range col.N, col.V k * cnew k = ∑ k ∈ range col.N, col.V k * c k := by
  unfold newStep aivaStep at hstep
  generalize tracerAdvectingW wt wb wc wF = W at hstep hW0 hWN
  generalize recon W c = R at hstep
  rw [implicit_solve_conserves col W cnew _ hV h0 hN hmask hstep]
  have hexp : ∀ k ∈ range col.N, col.V k * (c k + col.Δt * explicitTend col W R k) =
      col.V k * c k + col.Δt * (col.V k * explicitTend col W R k) := fun k _ => by ring
  rw [Finset.sum_congr rfl hexp, Finset.sum_add_distrib, ← Finset.mul_sum,
    explicit_conserves col W R hV hW0 hWN]
  ring

/-- Two-cell column used for the counterexample: Δt = 1, cfl = 1/2, V = Az = Δz = 1, only the interior face is active. -/
def exCol : Column where
  N := 2
  Δt := 1
  cfl := 1 / 2
  V := fun _ => 1
  Az := fun _ => 1
  Δzf := fun _ => 1
  act := fun f => if f = 1 then 1 else 0
  actc := fun _ => 1

/-- Velocity −1 (downward), Courant number 1 > cfl. -/
def exW : ℕ → ℝ := fun _ => -1

/-- Counterexample (pre-fix code): when the same total velocity is supplied as the two decompositions (transport = w, drift = 0) and (transport = 0, drift = w), the pre-fix implicit-part
coefficients (which see only the transport velocity) differ: one is −1/2, the other is 0. -/
theorem old_implicit_coeff_split_dependent :
    (∀ f, tracerAdvectingW exW 0 0 0 f = tracerAdvectingW 0 exW 0 0 f) ∧
    upperAdv exCol exW 0 ≠ upperAdv exCol 0 0 := by
  refine ⟨?_, ?_⟩
  · intro f
    simp only [tracerAdvectingW, Pi.zero_apply]
    ring
  · have h1 : upperAdv exCol exW 0 = -1 / 2 := by
      norm_num [upperAdv, wI, implicitVel, courant, exCol, exW, min_def]
    have h2 : upperAdv exCol 0 0 = 0 := by
      simp [upperAdv, wI, implicitVel]
    rw [h1, h2]
    norm_num

/-! ## §4 Bundled change: GM advective-form slope tapering (advective_skew_diffusion.jl) -/

/-- `bz = max(bz, convert(FT, slope_model.minimum_bz))` (isopycnal_skew_symmetric_diffusivity.jl:219,
advective_skew_diffusion.jl:71). -/
def bzFloor (bz minbz : ℝ) : ℝ := max bz minbz

/-- `Sx = - bx / bz + ∂x_z` (isopycnal_skew_symmetric_diffusivity.jl:221-222). -/
def rawSlope (bh bzf dhz : ℝ) : ℝ := -bh / bzf + dhz

/-- `Sx^2 + Sy^2`. -/
def taperQ (bx byv bzf dxz dyz : ℝ) : ℝ := rawSlope bx bzf dxz ^ 2 + rawSlope byv bzf dyz ^ 2

/-- ```
isopycnal_skew_symmetric_diffusivity.jl:216 @inline function calc_tapering(bx, by, bz, ∂x_z, ∂y_z, grid, slope_model, slope_limiter)
isopycnal_skew_symmetric_diffusivity.jl:217     FT = eltype(grid)
isopycnal_skew_symmetric_diffusivity.jl:218     Sₘ = convert(FT, slope_limiter.max_slope)
isopycnal_skew_symmetric_diffusivity.jl:219     bz = max(bz, convert(FT, slope_model.minimum_bz))
isopycnal_skew_symmetric_diffusivity.jl:221     Sx = - bx / bz + ∂x_z
isopycnal_skew_symmetric_diffusivity.jl:222     Sy = - by / bz + ∂y_z
isopycnal_skew_symmetric_diffusivity.jl:224     return ifelse(bz <= 0, zero(grid), min(one(grid), Sₘ^2 / (Sx^2 + Sy^2)))
```
IEEE isomorphism: when `Sx^2 + Sy^2 = 0`, Julia gives `Sₘ²/0 = +Inf` and `min(1, Inf) = 1` (when `Sₘ ≠ 0`); this is written here as an explicit branch;
when `Sₘ = 0` and `S² = 0`, Julia gives NaN, and this transcription does not cover that degenerate configuration. When `bz ≤ 0`, the `Sx` computed by Julia may be ±Inf/NaN, but
`ifelse` takes the 0 branch, so the result agrees with this transcription. -/
def calcTapering (bx byv bz dxz dyz minbz Sm : ℝ) : ℝ :=
  if bzFloor bz minbz ≤ 0 then 0
  else if taperQ bx byv (bzFloor bz minbz) dxz dyz = 0 then 1
  else min 1 (Sm ^ 2 / taperQ bx byv (bzFloor bz minbz) dxz dyz)

/-- ```
advective_skew_diffusion.jl:69 @inline function isopycnal_slope(bh, bz, ∂h_z, grid, slope_model)
advective_skew_diffusion.jl:70     FT = eltype(grid)
advective_skew_diffusion.jl:71     bz = max(bz, convert(FT, slope_model.minimum_bz))
advective_skew_diffusion.jl:72     return ifelse(bz <= 0, zero(grid), -bh / bz) + ∂h_z
``` -/
def isopycnalSlope (bh bz dhz minbz : ℝ) : ℝ :=
  (if bzFloor bz minbz ≤ 0 then 0 else -bh / bzFloor bz minbz) + dhz

/-- ```
advective_skew_diffusion.jl:86  @inline function ϵSxᶠᶜᶠ(i, j, k, grid, closure, b, C)
advective_skew_diffusion.jl:94      ϵ  = calc_tapering(bx, by, bz, ∂x_z, ∂y_z, grid, closure.isopycnal_tensor, closure.slope_limiter)
advective_skew_diffusion.jl:95      Sx = isopycnal_slope(bx, bz, ∂x_z, grid, closure.isopycnal_tensor)
advective_skew_diffusion.jl:98      inactive = peripheral_node(i, j, k, grid, Face(), Center(), Face())
advective_skew_diffusion.jl:100     return ifelse(inactive, zero(grid), ϵ * Sx)
```
(bx, by, bz are the interpolated values `bxᶠᶜᶠ/byᶠᶜᶠ/bzᶠᶜᶠ` (:76-78) at the (F,C,F) point, taken as inputs; `ϵSyᶜᶠᶠ` (:104-118) is symmetric.) -/
def epsSx (inactive : Bool) (bx byv bz dxz dyz minbz Sm : ℝ) : ℝ :=
  if inactive = true then 0
  else calcTapering bx byv bz dxz dyz minbz Sm * isopycnalSlope bx bz dxz minbz

/-- `ϵSyᶜᶠᶠ` (advective_skew_diffusion.jl:104-118). -/
def epsSy (inactive : Bool) (bx byv bz dxz dyz minbz Sm : ℝ) : ℝ :=
  if inactive = true then 0
  else calcTapering bx byv bz dxz dyz minbz Sm * isopycnalSlope byv bz dyz minbz

lemma slope_eq_raw (bh bz dhz minbz : ℝ) (h : ¬ bzFloor bz minbz ≤ 0) :
    isopycnalSlope bh bz dhz minbz = rawSlope bh (bzFloor bz minbz) dhz := by
  unfold isopycnalSlope rawSlope
  rw [if_neg h]

lemma taper_bound_core (Sm Q : ℝ) (hQ : 0 < Q) : (min 1 (Sm ^ 2 / Q)) ^ 2 * Q ≤ Sm ^ 2 := by
  have ht0 : 0 ≤ Sm ^ 2 / Q := div_nonneg (sq_nonneg _) hQ.le
  have hε0 : 0 ≤ min 1 (Sm ^ 2 / Q) := le_min zero_le_one ht0
  have hεt : min 1 (Sm ^ 2 / Q) ≤ Sm ^ 2 / Q := min_le_right _ _
  have hε1 : min 1 (Sm ^ 2 / Q) ≤ 1 := min_le_left _ _
  have hQ' : Q ≠ 0 := hQ.ne'
  have htQ : Sm ^ 2 / Q * Q = Sm ^ 2 := by field_simp
  by_cases hle : Q ≤ Sm ^ 2
  · have hε2 : (min 1 (Sm ^ 2 / Q)) ^ 2 ≤ 1 := by nlinarith
    have h3 := mul_le_mul_of_nonneg_right hε2 hQ.le
    linarith
  · push_neg at hle
    have ht1 : Sm ^ 2 / Q ≤ 1 := by nlinarith
    have h1 : (min 1 (Sm ^ 2 / Q)) ^ 2 ≤ (Sm ^ 2 / Q) ^ 2 := by nlinarith
    have h2 := mul_le_mul_of_nonneg_right h1 hQ.le
    have h3 : (Sm ^ 2 / Q) ^ 2 * Q = Sm ^ 2 / Q * Sm ^ 2 := by
      calc (Sm ^ 2 / Q) ^ 2 * Q = Sm ^ 2 / Q * (Sm ^ 2 / Q * Q) := by ring
        _ = Sm ^ 2 / Q * Sm ^ 2 := by rw [htQ]
    have h4 : Sm ^ 2 / Q * Sm ^ 2 ≤ Sm ^ 2 := by nlinarith [sq_nonneg Sm]
    linarith

/-- Certificate C1 (vector-form upper bound): the tapered full horizontal slope satisfies `ϵ²(Sx² + Sy²) ≤ Sₘ²`, for all inputs (including bz ≤ 0, any
minimum_bz). -/
theorem tapered_slope_vector_bound (bx byv bz dxz dyz minbz Sm : ℝ) :
    calcTapering bx byv bz dxz dyz minbz Sm ^ 2 *
      (isopycnalSlope bx bz dxz minbz ^ 2 + isopycnalSlope byv bz dyz minbz ^ 2) ≤ Sm ^ 2 := by
  unfold calcTapering
  split_ifs with h1 h2
  · rw [show (0 : ℝ) ^ 2 = 0 by norm_num, zero_mul]
    exact sq_nonneg _
  · rw [slope_eq_raw bx bz dxz minbz h1, slope_eq_raw byv bz dyz minbz h1]
    unfold taperQ at h2
    rw [h2, mul_zero]
    exact sq_nonneg _
  · rw [slope_eq_raw bx bz dxz minbz h1, slope_eq_raw byv bz dyz minbz h1]
    have hQnn : 0 ≤ taperQ bx byv (bzFloor bz minbz) dxz dyz := by
      unfold taperQ
      positivity
    have hQpos : 0 < taperQ bx byv (bzFloor bz minbz) dxz dyz := lt_of_le_of_ne hQnn (Ne.symm h2)
    have hcore := taper_bound_core Sm _ hQpos
    unfold taperQ at hcore ⊢
    exact hcore

/-- Certificate C2 (component upper bound): `|ϵ Sx| ≤ Sₘ` (including the peripheral-point zeroing branch). -/
theorem epsSx_abs_le (inactive : Bool) (bx byv bz dxz dyz minbz Sm : ℝ) (hSm : 0 ≤ Sm) :
    |epsSx inactive bx byv bz dxz dyz minbz Sm| ≤ Sm := by
  unfold epsSx
  split_ifs
  · simpa using hSm
  · have hv := tapered_slope_vector_bound bx byv bz dxz dyz minbz Sm
    generalize calcTapering bx byv bz dxz dyz minbz Sm = ε at hv ⊢
    generalize isopycnalSlope bx bz dxz minbz = a at hv ⊢
    generalize isopycnalSlope byv bz dyz minbz = b at hv ⊢
    have h : (ε * a) ^ 2 ≤ Sm ^ 2 := by nlinarith [sq_nonneg (ε * b)]
    rw [abs_le]
    constructor
    · nlinarith
    · nlinarith

/-- Certificate C2′: `|ϵ Sy| ≤ Sₘ`. -/
theorem epsSy_abs_le (inactive : Bool) (bx byv bz dxz dyz minbz Sm : ℝ) (hSm : 0 ≤ Sm) :
    |epsSy inactive bx byv bz dxz dyz minbz Sm| ≤ Sm := by
  unfold epsSy
  split_ifs
  · simpa using hSm
  · have hv := tapered_slope_vector_bound bx byv bz dxz dyz minbz Sm
    generalize calcTapering bx byv bz dxz dyz minbz Sm = ε at hv ⊢
    generalize isopycnalSlope bx bz dxz minbz = a at hv ⊢
    generalize isopycnalSlope byv bz dyz minbz = b at hv ⊢
    have h : (ε * b) ^ 2 ≤ Sm ^ 2 := by nlinarith [sq_nonneg (ε * a)]
    rw [abs_le]
    constructor
    · nlinarith
    · nlinarith

/-- Certificate C3: under unstable (including neutral) stratification `bz ≤ 0` with `minimum_bz ≤ 0` (default `minimum_bz = 0`,
isopycnal_rotation_tensor_components.jl:66), `ϵ Sx = 0`. -/
theorem epsSx_unstable_zero (inactive : Bool) (bx byv bz dxz dyz minbz Sm : ℝ) (hbz : bz ≤ 0)
    (hmin : minbz ≤ 0) : epsSx inactive bx byv bz dxz dyz minbz Sm = 0 := by
  unfold epsSx
  split_ifs
  · rfl
  · have hf : bzFloor bz minbz ≤ 0 := max_le hbz hmin
    unfold calcTapering
    rw [if_pos hf, zero_mul]

/-- Counterexample (domain of applicability): when `minimum_bz > 0`, unstable stratification is not zeroed: `max(bz, minimum_bz)` lifts bz to a positive value, and the slope is
computed as `−bx/minimum_bz` and then tapered. With bx = −1, by = 0, bz = −1, ∂z = 0, minimum_bz = 1, Sₘ = 1, `ϵ Sx = 1`. -/
theorem epsSx_floor_counterexample : epsSx false (-1) 0 (-1) 0 0 1 1 = 1 := by
  norm_num [epsSx, calcTapering, isopycnalSlope, bzFloor, taperQ, rawSlope, min_def, max_def]

/-- R₁₃ of the diffusive form (SmallSlopeIsopycnalTensor).
```
isopycnal_rotation_tensor_components.jl:68 @inline function isopycnal_rotation_tensor_xz_fcc(i, j, k, grid::AbstractGrid, buoyancy, tracers,
isopycnal_rotation_tensor_components.jl:70     bx   = ∂xᵣ_b(i, j, k, grid, buoyancy, tracers)
isopycnal_rotation_tensor_components.jl:71     bz   = ℑxzᶠᵃᶜ(i, j, k, grid, ∂z_b, buoyancy, tracers)
isopycnal_rotation_tensor_components.jl:72     bz   = max(bz, slope_model.minimum_bz)
isopycnal_rotation_tensor_components.jl:73     ∂x_z = ∂x_zᶠᶜᶜ(i, j, k, grid)
isopycnal_rotation_tensor_components.jl:74     Sx   = - bx / bz
isopycnal_rotation_tensor_components.jl:76     return ifelse(bz == 0, zero(grid), Sx) + ∂x_z
```
(The diffusive flux is multiplied by `ϵ = tapering_factorᶠᶜᶜ = calc_tapering(...)`, isopycnal_skew_symmetric_diffusivity.jl:183-191, :262-264.) -/
def rotationXZ (bx bz dxz minbz : ℝ) : ℝ :=
  (if bzFloor bz minbz = 0 then 0 else -bx / bzFloor bz minbz) + dxz

/-- Certificate C4 (consistent with the diffusive form): for the same inputs, `ϵ·isopycnal_slope` of the advective form and `ϵ·R₁₃` of the diffusive form are pointwise equal
(the same `calc_tapering`, the same bz lower bound), for any minimum_bz. The branch conditions `bz <= 0` and `bz == 0`
differ, but when `bz ≤ 0`, ϵ = 0 and the products agree. (The two forms take values at different grid locations: the advective form at (F,C,F), the diffusive form at (F,C,C);
this theorem is pointwise formula consistency.) -/
theorem advective_matches_diffusive (bx byv bz dxz dyz minbz Sm : ℝ) :
    calcTapering bx byv bz dxz dyz minbz Sm * isopycnalSlope bx bz dxz minbz =
      calcTapering bx byv bz dxz dyz minbz Sm * rotationXZ bx bz dxz minbz := by
  by_cases h : bzFloor bz minbz ≤ 0
  · unfold calcTapering
    rw [if_pos h, zero_mul, zero_mul]
  · have hne : bzFloor bz minbz ≠ 0 := fun h0 => h (le_of_eq h0)
    unfold isopycnalSlope rotationXZ
    rw [if_neg h, if_neg hne]

/-! ### The pre-fix (parent commit 1bd69f4e) advective form, used to characterize the behavior change -/

/-- ```
advective_skew_diffusion.jl(parent commit):44 @inline function Sxᶠᶜᶠ(i, j, k, grid, b, C)
advective_skew_diffusion.jl(parent commit):48     Sx   = ifelse(bz == 0, ∂x_z, - bx / bz + ∂x_z)
``` -/
def oldSx (bx bz dxz : ℝ) : ℝ := if bz = 0 then dxz else -bx / bz + dxz

/-- ```
advective_skew_diffusion.jl(parent commit):37 @inline function tapering_factor(Sx, Sy, slope_limiter::FluxTapering)
advective_skew_diffusion.jl(parent commit):38     S²  = Sx^2 + Sy^2
advective_skew_diffusion.jl(parent commit):39     Sₘ² = slope_limiter.max_slope^2
advective_skew_diffusion.jl(parent commit):40     return min(one(Sx), Sₘ² / S²)
```
Called with `Sy = 0` (parent commit :74); the IEEE branch for `S² = 0` is the same as in `calcTapering`. -/
def oldTaper (Sx Sm : ℝ) : ℝ := if Sx ^ 2 + 0 ^ 2 = 0 then 1 else min 1 (Sm ^ 2 / (Sx ^ 2 + 0 ^ 2))

/-- Parent commit `ϵSxᶠᶜᶠ` (:72-76): `ϵ = tapering_factor(Sx, zero(grid), slope_limiter); return ϵ * Sx`
(the peripheral-point Sx has already been zeroed in `Sxᶠᶜᶠ`). -/
def oldEpsSx (inactive : Bool) (bx bz dxz Sm : ℝ) : ℝ :=
  if inactive = true then 0 else oldTaper (oldSx bx bz dxz) Sm * oldSx bx bz dxz

lemma calcTapering_nonneg (bx byv bz dxz dyz minbz Sm : ℝ) :
    0 ≤ calcTapering bx byv bz dxz dyz minbz Sm := by
  unfold calcTapering
  split_ifs
  · exact le_refl 0
  · exact zero_le_one
  · exact le_min zero_le_one (div_nonneg (sq_nonneg _) (by unfold taperQ; positivity))

lemma calcTapering_le_one (bx byv bz dxz dyz minbz Sm : ℝ) :
    calcTapering bx byv bz dxz dyz minbz Sm ≤ 1 := by
  unfold calcTapering
  split_ifs
  · exact zero_le_one
  · exact le_refl 1
  · exact min_le_left _ _

/-- Certificate C5 (stable stratification, lower bound inactive: new ≤ old): when bz > 0 and minimum_bz ≤ bz, the new |ϵSx| does not exceed the pre-fix |ϵSx|
(tapering on the full horizontal slope is stronger than tapering on a single component). -/
theorem new_le_old (bx byv bz dxz dyz minbz Sm : ℝ) (hbz : 0 < bz) (hmin : minbz ≤ bz) :
    |epsSx false bx byv bz dxz dyz minbz Sm| ≤ |oldEpsSx false bx bz dxz Sm| := by
  have hf : bzFloor bz minbz = bz := max_eq_left hmin
  have hnot : ¬ bzFloor bz minbz ≤ 0 := by
    rw [hf]
    exact not_le.mpr hbz
  have hS : isopycnalSlope bx bz dxz minbz = rawSlope bx bz dxz := by
    rw [slope_eq_raw bx bz dxz minbz hnot, hf]
  have hSo : oldSx bx bz dxz = rawSlope bx bz dxz := by
    unfold oldSx rawSlope
    rw [if_neg hbz.ne']
  unfold epsSx oldEpsSx
  rw [if_neg Bool.false_ne_true, if_neg Bool.false_ne_true, hS, hSo, abs_mul, abs_mul]
  refine mul_le_mul_of_nonneg_right ?_ (abs_nonneg _)
  have hn := calcTapering_nonneg bx byv bz dxz dyz minbz Sm
  have hold0 : 0 ≤ oldTaper (rawSlope bx bz dxz) Sm := by
    unfold oldTaper
    split_ifs
    · exact zero_le_one
    · exact le_min zero_le_one (div_nonneg (sq_nonneg _) (by positivity))
  rw [abs_of_nonneg hn, abs_of_nonneg hold0]
  unfold oldTaper
  by_cases hS0 : rawSlope bx bz dxz ^ 2 + 0 ^ 2 = 0
  · rw [if_pos hS0]
    exact calcTapering_le_one bx byv bz dxz dyz minbz Sm
  · rw [if_neg hS0]
    unfold calcTapering
    rw [if_neg hnot, hf]
    have hpos : 0 < rawSlope bx bz dxz ^ 2 + 0 ^ 2 :=
      lt_of_le_of_ne (by positivity) (Ne.symm hS0)
    have hle : rawSlope bx bz dxz ^ 2 + 0 ^ 2 ≤ taperQ bx byv bz dxz dyz := by
      unfold taperQ
      nlinarith [sq_nonneg (rawSlope byv bz dyz)]
    have hQ : ¬ taperQ bx byv bz dxz dyz = 0 := by
      intro h0
      rw [h0] at hle
      linarith
    rw [if_neg hQ]
    apply min_le_min (le_refl 1)
    have h1 := one_div_le_one_div_of_le hpos hle
    calc Sm ^ 2 / taperQ bx byv bz dxz dyz = Sm ^ 2 * (1 / taperQ bx byv bz dxz dyz) := by ring
      _ ≤ Sm ^ 2 * (1 / (rawSlope bx bz dxz ^ 2 + 0 ^ 2)) :=
          mul_le_mul_of_nonneg_left h1 (sq_nonneg _)
      _ = Sm ^ 2 / (rawSlope bx bz dxz ^ 2 + 0 ^ 2) := by ring

/-- One counterexample (to "the bundled change does not alter numerical results"): with stable stratification, both horizontal slope components equal to 1, Sₘ = 1/2,
the new ϵSx = 1/8 and the pre-fix ϵSx = 1/4. -/
theorem tapering_change_witness :
    epsSx false (-1) (-1) 1 0 0 0 (1 / 2) = 1 / 8 ∧ oldEpsSx false (-1) 1 0 (1 / 2) = 1 / 4 := by
  constructor
  · norm_num [epsSx, calcTapering, isopycnalSlope, bzFloor, taperQ, rawSlope, min_def, max_def]
  · norm_num [oldEpsSx, oldTaper, oldSx, min_def]

/-- A second counterexample: with unstable stratification bz = −1 (minimum_bz = 0), the new ϵSx = 0 and the pre-fix ϵSx = −1/4 (the old code still has
GM eddy-induced transport under unstable stratification). -/
theorem unstable_change_witness :
    epsSx false (-1) 0 (-1) 0 0 0 (1 / 2) = 0 ∧ oldEpsSx false (-1) (-1) 0 (1 / 2) = -1 / 4 := by
  constructor
  · norm_num [epsSx, calcTapering, isopycnalSlope, bzFloor, taperQ, rawSlope, min_def, max_def]
  · norm_num [oldEpsSx, oldTaper, oldSx, min_def]

end

end Oceananigans6083
