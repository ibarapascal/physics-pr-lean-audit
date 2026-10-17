-- Moose34020: idaholab/moose PR #34020 "Pass the nodal slip rate to dynamic mortar function_friction"
-- (merge commit f90f508580f7611f7d75ee458cd177affc7288d1) Lean certificate.
-- Transcribes the nodal slip rate in modules/contact ComputeDynamicFrictionalForceLMMechanicalContact.C
-- (weighted tangential velocity / integral of the same test function) and computeFrictionValue, as well as
-- the accumulation of the normalization quantity ∫ψ_j in ComputeDynamicWeightedGapLMMechanicalContact.C;
-- for arbitrary finite families of quadrature points (any mesh, element order, quadrature rule, mortar segmentation, process partition), proves exact reproduction of constant fields,
-- exact reproduction of linear fields under the dual basis, (approximate) equivalence of the velocity-dependent friction law to a constant μ under uniform slip, and gives a counterexample for the old implementation.
import Mathlib

/-!
# MOOSE #34020 — nodal slip rate for dynamic mortar `function_friction`

Transcription conventions
* The accumulated terms at one node (dof) are represented by a finite index set `S : Finset ι`: each index `k` is one
    `computeQpIProperties` call (a mortar segment, a quadrature point `_qp`, the local index `_i` of the node).
    Cross-process `communicateGaps` / `communicateVelocities` merely add the per-process partial sums, so the value on the owning process
    is the sum over all of `S`; if terms are accumulated more than once between two `residualSetup` calls, that just enlarges `S`.
    (The communication layer is not transcribed; it is taken as the assumption "value on the owner = sum of all contributions".)
* `_test[_i][_qp]`, `_JxW_msm[_qp]`, `_coord[_qp]`, `_relative_velocity`, `nodal_tangents[d][_i]`
    are all treated as arbitrary real number / vector inputs at each index, with no sign assumption (test functions of the dual basis may be negative).
* `_function_friction->value(t, x, y, z)` is taken as an arbitrary real function `g : ℝ → ℝ → ℝ → ℝ → ℝ` (this layer is not transcribed).
* AD derivatives and `MetaPhysicL::raw_value` do not affect the numerical value and are transcribed as real values; floating-point rounding is checked by the differential script
    `replay/scripts/Moose34020_replay.py`.
-/

namespace PhysAudit.Moose34020

noncomputable section

open Finset

/-! ## §1 Transcription -/

/-- libMesh `RealVectorValue` (three components). -/
abbrev Vec3 := Fin 3 → ℝ

/-- libMesh vector `operator*` (dot product). -/
def dot3 (u v : Vec3) : ℝ := u 0 * v 0 + u 1 * v 1 + u 2 * v 2

/-- Scalar times vector (libMesh `Real * RealVectorValue`). -/
def smul3 (s : ℝ) (v : Vec3) : Vec3 := fun c => s * v c

/-- ```
ComputeDynamicFrictionalForceLMMechanicalContact.C:124
  _qp_tangential_velocity_nodal = -_relative_velocity * (_JxW_msm[_qp] * _coord[_qp]);
``` -/
def qpTangentialVelocityNodal (vrel : Vec3) (JxW coord : ℝ) : Vec3 :=
  fun c => -vrel c * (JxW * coord)

/-- ```
ComputeDynamicWeightedGapLMMechanicalContact.C:192
  _qp_factor = _JxW_msm[_qp] * _coord[_qp];
``` -/
def qpFactor (JxW coord : ℝ) : ℝ := JxW * coord

/-- Weighted tangential velocity (direction `tangent` = `nodal_tangents[0]` or `nodal_tangents[1]`).
```
ComputeDynamicFrictionalForceLMMechanicalContact.C:141-142
  _dof_to_weighted_tangential_velocity[dof][0] +=
      _test[_i][_qp] * _qp_tangential_velocity_nodal * nodal_tangents[0][_i];
ComputeDynamicFrictionalForceLMMechanicalContact.C:147-148
    _dof_to_weighted_tangential_velocity[dof][1] +=
        _test[_i][_qp] * _qp_tangential_velocity_nodal * nodal_tangents[1][_i];
```
The C++ evaluation order is `(test * vec) · t`; here it is written as `dot3 (smul3 test vec) t`. -/
def weightedTangentialVelocity {ι : Type*} (S : Finset ι) (test JxW coord : ι → ℝ)
    (vrel tangent : ι → Vec3) : ℝ :=
  ∑ k ∈ S, dot3 (smul3 (test k) (qpTangentialVelocityNodal (vrel k) (JxW k) (coord k))) (tangent k)

/-- Normalization quantity ∫ψ_j (accumulated unconditionally since the PR).
```
ComputeDynamicWeightedGapLMMechanicalContact.C:222-223
  // Always accumulated: also used to compute nodal slip rates for friction
  _dof_to_weighted_gap[dof].second += _test[_i][_qp] * _qp_factor;
``` -/
def normalization {ι : Type*} (S : Finset ι) (test JxW coord : ι → ℝ) : ℝ :=
  ∑ k ∈ S, test k * qpFactor (JxW k) (coord k)

/-- Nodal slip rate (one tangential component).
```
ComputeDynamicFrictionalForceLMMechanicalContact.C:183-186
    const Real normalization = libmesh_map_find(_dof_to_weighted_gap, dof).second;
    _dof_to_old_real_tangential_velocity[dof] = {
        {MetaPhysicL::raw_value(weighted_velocities[0]) / normalization,
         _3d ? MetaPhysicL::raw_value(weighted_velocities[1]) / normalization : 0.0}};
```
Note: under real-number semantics `x / 0 = 0`, whereas in C++ it is ±inf/NaN; all theorems assume the normalization quantity is nonzero. -/
def nodalSlipRate {ι : Type*} (S : Finset ι) (test JxW coord : ι → ℝ) (vrel tangent : ι → Vec3) :
    ℝ :=
  weightedTangentialVelocity S test JxW coord vrel tangent / normalization S test JxW coord

/-- Two components (`_3d ? ... : 0.0`), see the quotation in `nodalSlipRate` :184-186. -/
def nodalSlipRates {ι : Type*} (is3d : Bool) (S : Finset ι) (test JxW coord : ι → ℝ)
    (vrel t0 t1 : ι → Vec3) : ℝ × ℝ :=
  (nodalSlipRate S test JxW coord vrel t0,
    if is3d then nodalSlipRate S test JxW coord vrel t1 else 0)

/-- ```
ComputeDynamicFrictionalForceLMMechanicalContact.C:392-394
  const auto it = _dof_to_old_real_tangential_velocity.find(dof);
  return it != _dof_to_old_real_tangential_velocity.end() ? it->second
                                                          : std::array<Real, 2>{{0.0, 0.0}};
``` -/
def oldNodalSlipRate (history : Option (ℝ × ℝ)) : ℝ × ℝ := history.getD (0, 0)

/-- ```
ComputeDynamicFrictionalForceLMMechanicalContact.C:374-375
    ADReal tangential_vel_magnitude =
        sqrt(tangential_vel * tangential_vel + tangential_vel_dir * tangential_vel_dir + 1.0e-24);
``` -/
def slipMagnitude (v0 v1 : ℝ) : ℝ := Real.sqrt (v0 * v0 + v1 * v1 + 1e-24)

/-- Function branch.
```
ComputeDynamicFrictionalForceLMMechanicalContact.C:377
    mu_ad = _function_friction->value<ADReal>(0.0, contact_pressure, tangential_vel_magnitude, 0.0);
``` -/
def frictionValueFn (g : ℝ → ℝ → ℝ → ℝ → ℝ) (p v0 v1 : ℝ) : ℝ :=
  g 0 p (slipMagnitude v0 v1) 0

/-- ```
ComputeDynamicFrictionalForceLMMechanicalContact.C:368-378
  if (!_has_friction_function)
    mu_ad = _mu;
  else
  { ... (see slipMagnitude / frictionValueFn) }
```
Call sites: :288-289 (3D, `old_slip_rate[0], old_slip_rate[1]`) and :339 (2D, `oldNodalSlipRate(dof)[0], 0.0`). -/
def frictionValue (hasFunction : Bool) (mu : ℝ) (g : ℝ → ℝ → ℝ → ℝ → ℝ) (p v0 v1 : ℝ) : ℝ :=
  if hasFunction then frictionValueFn g p v0 v1 else mu

/-- Accumulation of the "true tangential velocity" in the old implementation (parent commit f3250a4703d405f2dc08acea0ec56b2760b79cd5), used for the counterexample.
```
ComputeDynamicFrictionalForceLMMechanicalContact.C (parent commit):124
  _qp_real_tangential_velocity_nodal = -_relative_velocity;
ComputeDynamicFrictionalForceLMMechanicalContact.C (parent commit):144-146
  _dof_to_real_tangential_velocity[dof][0] +=
      _test[_i][_qp] * MetaPhysicL::raw_value(_qp_real_tangential_velocity_nodal) *
      nodal_tangents[0][_i];
```
The old code uses this sum as-is as the old slip rate in timestepSetup (no JxW, no normalization). -/
def oldRealTangentialVelocity {ι : Type*} (S : Finset ι) (test : ι → ℝ) (vrel tangent : ι → Vec3) :
    ℝ :=
  ∑ k ∈ S, dot3 (smul3 (test k) (fun c => -vrel k c)) (tangent k)

/-! ## §2 Lemmas at the definition level -/

theorem dot3_smul3 (s : ℝ) (u v : Vec3) : dot3 (smul3 s u) v = s * dot3 u v := by
  simp only [dot3, smul3]; ring

theorem dot3_qpTangential (vrel : Vec3) (JxW coord : ℝ) (t : Vec3) :
    dot3 (qpTangentialVelocityNodal vrel JxW coord) t = -(qpFactor JxW coord) * dot3 vrel t := by
  simp only [dot3, qpTangentialVelocityNodal, qpFactor]; ring

/-- Each term = `−ψ_k · w_k · (v_k · t_k)`; `w_k = JxW·coord` is exactly the same weight as in the normalization quantity. -/
theorem weightedTangentialVelocity_eq {ι : Type*} (S : Finset ι) (test JxW coord : ι → ℝ)
    (vrel tangent : ι → Vec3) :
    weightedTangentialVelocity S test JxW coord vrel tangent =
      -∑ k ∈ S, test k * qpFactor (JxW k) (coord k) * dot3 (vrel k) (tangent k) := by
  unfold weightedTangentialVelocity
  rw [← Finset.sum_neg_distrib]
  apply Finset.sum_congr rfl
  intro k _
  rw [dot3_smul3, dot3_qpTangential]
  ring

theorem nodalSlipRates_fst {ι : Type*} (is3d : Bool) (S : Finset ι) (test JxW coord : ι → ℝ)
    (vrel t0 t1 : ι → Vec3) :
    (nodalSlipRates is3d S test JxW coord vrel t0 t1).1 = nodalSlipRate S test JxW coord vrel t0 :=
  rfl

theorem nodalSlipRates_false_snd {ι : Type*} (S : Finset ι) (test JxW coord : ι → ℝ)
    (vrel t0 t1 : ι → Vec3) : (nodalSlipRates false S test JxW coord vrel t0 t1).2 = 0 :=
  rfl

theorem nodalSlipRates_true_snd {ι : Type*} (S : Finset ι) (test JxW coord : ι → ℝ)
    (vrel t0 t1 : ι → Vec3) :
    (nodalSlipRates true S test JxW coord vrel t0 t1).2 = nodalSlipRate S test JxW coord vrel t1 :=
  rfl

theorem frictionValue_true (mu : ℝ) (g : ℝ → ℝ → ℝ → ℝ → ℝ) (p v0 v1 : ℝ) :
    frictionValue true mu g p v0 v1 = frictionValueFn g p v0 v1 :=
  rfl

/-! ## §3 Certificate: exact reproduction of constant fields (primary invariant A / I1) -/

/-- **Certificate**: for any finite family of contributions (any mesh, element order, quadrature rule, mortar segmentation, dual or standard basis,
weights of arbitrary sign), as long as the normalization quantity is nonzero and the tangential velocity component `v_k · t_k` of every term equals the same constant `a`,
the nodal slip rate is exactly `−a`. -/
theorem nodalSlipRate_const_field {ι : Type*} (S : Finset ι) (test JxW coord : ι → ℝ)
    (vrel tangent : ι → Vec3) (a : ℝ)
    (ha : ∀ k ∈ S, dot3 (vrel k) (tangent k) = a)
    (hN : normalization S test JxW coord ≠ 0) :
    nodalSlipRate S test JxW coord vrel tangent = -a := by
  have hW : weightedTangentialVelocity S test JxW coord vrel tangent =
      -a * normalization S test JxW coord := by
    rw [weightedTangentialVelocity_eq]
    unfold normalization
    rw [Finset.mul_sum, ← Finset.sum_neg_distrib]
    apply Finset.sum_congr rfl
    intro k hk
    rw [ha k hk]
    ring
  unfold nodalSlipRate
  rw [hW, mul_div_assoc, div_self hN, mul_one]

/-- **Certificate**: for rigid uniform slip (relative velocity constantly `V`, nodal tangent constantly `τ`), the nodal slip rate = `−V·τ`. -/
theorem nodalSlipRate_uniform {ι : Type*} (S : Finset ι) (test JxW coord : ι → ℝ)
    (vrel tangent : ι → Vec3) (V τ : Vec3)
    (hv : ∀ k ∈ S, vrel k = V) (ht : ∀ k ∈ S, tangent k = τ)
    (hN : normalization S test JxW coord ≠ 0) :
    nodalSlipRate S test JxW coord vrel tangent = -dot3 V τ :=
  nodalSlipRate_const_field S test JxW coord vrel tangent (dot3 V τ)
    (fun k hk => by rw [hv k hk, ht k hk]) hN

/-- **Certificate**: mesh independence. Two arbitrary discretizations (different meshes / quadrature rules / basis functions) under the same uniform slip field
give exactly the same nodal slip rate. -/
theorem nodalSlipRate_mesh_independent {ι κ : Type*}
    (S : Finset ι) (test JxW coord : ι → ℝ) (vrel tangent : ι → Vec3)
    (S' : Finset κ) (test' JxW' coord' : κ → ℝ) (vrel' tangent' : κ → Vec3) (V τ : Vec3)
    (hv : ∀ k ∈ S, vrel k = V) (ht : ∀ k ∈ S, tangent k = τ)
    (hv' : ∀ k ∈ S', vrel' k = V) (ht' : ∀ k ∈ S', tangent' k = τ)
    (hN : normalization S test JxW coord ≠ 0) (hN' : normalization S' test' JxW' coord' ≠ 0) :
    nodalSlipRate S test JxW coord vrel tangent = nodalSlipRate S' test' JxW' coord' vrel' tangent' := by
  rw [nodalSlipRate_uniform S test JxW coord vrel tangent V τ hv ht hN,
    nodalSlipRate_uniform S' test' JxW' coord' vrel' tangent' V τ hv' ht' hN']

/-! ## §4 Certificate: exact reproduction of finite-element fields under the dual basis (secondary invariant C / I2) -/

/-- **Certificate (conditional on the dual-basis biorthogonality assumption)**: given the node set `M` and the shape-function values at the quadrature points `phi m k`. If
* biorthogonality: `∑_k ψ_j(k) w_k φ_m(k) = δ_{mj} D` (`D ≠ 0`),
* partition of unity: `∑_m φ_m(k) = 1`,
* the tangential velocity component is a finite-element interpolant: `v_k · t_k = ∑_m φ_m(k) c_m`,
then the slip rate of node `j` is exactly `−c_j`. Biorthogonality and partition of unity are not transcribed from libMesh; they are taken as assumptions. -/
theorem nodalSlipRate_dual_exact {ι κ : Type*} [DecidableEq κ] (S : Finset ι) (M : Finset κ)
    (test JxW coord : ι → ℝ) (vrel tangent : ι → Vec3) (phi : κ → ι → ℝ) (cnode : κ → ℝ)
    (j : κ) (hj : j ∈ M) (D : ℝ) (hD : D ≠ 0)
    (hbiorth : ∀ m ∈ M,
      ∑ k ∈ S, test k * qpFactor (JxW k) (coord k) * phi m k = if m = j then D else 0)
    (hpu : ∀ k ∈ S, ∑ m ∈ M, phi m k = 1)
    (hinterp : ∀ k ∈ S, dot3 (vrel k) (tangent k) = ∑ m ∈ M, phi m k * cnode m) :
    nodalSlipRate S test JxW coord vrel tangent = -cnode j := by
  -- ∑_m c_m (if m = j then D else 0) = c_j D
  have hsel : ∀ c : κ → ℝ, ∑ m ∈ M, c m * (if m = j then D else 0) = c j * D := by
    intro c
    simp only [mul_ite, mul_zero]
    rw [Finset.sum_ite_eq', if_pos hj]
  have hW : weightedTangentialVelocity S test JxW coord vrel tangent = -(cnode j * D) := by
    rw [weightedTangentialVelocity_eq]
    congr 1
    calc ∑ k ∈ S, test k * qpFactor (JxW k) (coord k) * dot3 (vrel k) (tangent k)
        = ∑ k ∈ S, ∑ m ∈ M, cnode m * (test k * qpFactor (JxW k) (coord k) * phi m k) := by
          apply Finset.sum_congr rfl
          intro k hk
          rw [hinterp k hk, Finset.mul_sum]
          apply Finset.sum_congr rfl
          intro m _
          ring
      _ = ∑ m ∈ M, ∑ k ∈ S, cnode m * (test k * qpFactor (JxW k) (coord k) * phi m k) :=
          Finset.sum_comm
      _ = ∑ m ∈ M, cnode m * ∑ k ∈ S, test k * qpFactor (JxW k) (coord k) * phi m k := by
          apply Finset.sum_congr rfl
          intro m _
          rw [Finset.mul_sum]
      _ = ∑ m ∈ M, cnode m * (if m = j then D else 0) :=
          Finset.sum_congr rfl (fun m hm => by rw [hbiorth m hm])
      _ = cnode j * D := hsel cnode
  have hN : normalization S test JxW coord = D := by
    unfold normalization
    calc ∑ k ∈ S, test k * qpFactor (JxW k) (coord k)
        = ∑ k ∈ S, ∑ m ∈ M, test k * qpFactor (JxW k) (coord k) * phi m k := by
          apply Finset.sum_congr rfl
          intro k hk
          rw [← Finset.mul_sum, hpu k hk, mul_one]
      _ = ∑ m ∈ M, ∑ k ∈ S, test k * qpFactor (JxW k) (coord k) * phi m k := Finset.sum_comm
      _ = ∑ m ∈ M, (if m = j then D else 0) := Finset.sum_congr rfl (fun m hm => hbiorth m hm)
      _ = D := by rw [Finset.sum_ite_eq', if_pos hj]
  unfold nodalSlipRate
  rw [hW, hN, neg_div, mul_div_assoc, div_self hD, mul_one]

/-! ## §5 Certificates: value passed to the friction function under uniform slip (algebraic part of primary invariant B) -/

/-- **Certificate (2D)**: under uniform slip, on any mesh, the value passed to `function_friction` at every node is
`g(0, p, √((V·τ)² + 10⁻²⁴), 0)`, independent of the discretization. -/
theorem frictionValue_uniform_2d {ι : Type*} (S : Finset ι) (test JxW coord : ι → ℝ)
    (vrel t0 t1 : ι → Vec3) (V τ : Vec3)
    (hv : ∀ k ∈ S, vrel k = V) (ht : ∀ k ∈ S, t0 k = τ)
    (hN : normalization S test JxW coord ≠ 0) (mu : ℝ) (g : ℝ → ℝ → ℝ → ℝ → ℝ) (p : ℝ) :
    frictionValue true mu g p (nodalSlipRates false S test JxW coord vrel t0 t1).1
        (nodalSlipRates false S test JxW coord vrel t0 t1).2 =
      g 0 p (Real.sqrt (dot3 V τ ^ 2 + 1e-24)) 0 := by
  rw [frictionValue_true, nodalSlipRates_fst, nodalSlipRates_false_snd,
    nodalSlipRate_uniform S test JxW coord vrel t0 V τ hv ht hN]
  unfold frictionValueFn slipMagnitude
  have e : -dot3 V τ * -dot3 V τ + 0 * 0 + (1e-24 : ℝ) = dot3 V τ ^ 2 + 1e-24 := by ring
  rw [e]

/-- **Certificate (3D)**: when the two tangents `τ0, τ1` are orthonormal and the relative velocity `V = α τ0 + β τ1` lies in the tangent plane,
the value passed to `function_friction` at every node is `g(0, p, √(|V|² + 10⁻²⁴), 0)`. -/
theorem frictionValue_uniform_3d {ι : Type*} (S : Finset ι) (test JxW coord : ι → ℝ)
    (vrel t0 t1 : ι → Vec3) (τ0 τ1 : Vec3) (α β : ℝ)
    (h00 : dot3 τ0 τ0 = 1) (h01 : dot3 τ0 τ1 = 0) (h11 : dot3 τ1 τ1 = 1)
    (hv : ∀ k ∈ S, vrel k = fun c => α * τ0 c + β * τ1 c)
    (ht0 : ∀ k ∈ S, t0 k = τ0) (ht1 : ∀ k ∈ S, t1 k = τ1)
    (hN : normalization S test JxW coord ≠ 0) (mu : ℝ) (g : ℝ → ℝ → ℝ → ℝ → ℝ) (p : ℝ) :
    frictionValue true mu g p (nodalSlipRates true S test JxW coord vrel t0 t1).1
        (nodalSlipRates true S test JxW coord vrel t0 t1).2 =
      g 0 p (Real.sqrt (dot3 (fun c => α * τ0 c + β * τ1 c) (fun c => α * τ0 c + β * τ1 c)
        + 1e-24)) 0 := by
  rw [frictionValue_true, nodalSlipRates_fst, nodalSlipRates_true_snd,
    nodalSlipRate_uniform S test JxW coord vrel t0 _ τ0 hv ht0 hN,
    nodalSlipRate_uniform S test JxW coord vrel t1 _ τ1 hv ht1 hN]
  have hd0 : dot3 (fun c => α * τ0 c + β * τ1 c) τ0 = α := by
    simp only [dot3] at h00 h01 ⊢
    linear_combination α * h00 + β * h01
  have hd1 : dot3 (fun c => α * τ0 c + β * τ1 c) τ1 = β := by
    simp only [dot3] at h01 h11 ⊢
    linear_combination α * h01 + β * h11
  have hVV : dot3 (fun c => α * τ0 c + β * τ1 c) (fun c => α * τ0 c + β * τ1 c) = α ^ 2 + β ^ 2 := by
    simp only [dot3] at h00 h01 h11 ⊢
    linear_combination α ^ 2 * h00 + 2 * α * β * h01 + β ^ 2 * h11
  rw [hd0, hd1, hVV]
  unfold frictionValueFn slipMagnitude
  have e : -α * -α + -β * -β + (1e-24 : ℝ) = α ^ 2 + β ^ 2 + 1e-24 := by ring
  rw [e]

/-- Node with no history (first step or just entering contact): the friction function is evaluated at `y = √10⁻²⁴`. -/
theorem frictionValue_no_history (mu : ℝ) (g : ℝ → ℝ → ℝ → ℝ → ℝ) (p : ℝ) :
    frictionValue true mu g p (oldNodalSlipRate none).1 (oldNodalSlipRate none).2 =
      g 0 p (Real.sqrt 1e-24) 0 := by
  rw [frictionValue_true]
  unfold frictionValueFn slipMagnitude oldNodalSlipRate
  simp

/-- Upper and lower bounds of the regularization term `1e-24`: `|a| < √(a² + 10⁻²⁴) ≤ |a| + 10⁻¹²`. -/
theorem sqrt_sq_add_eps_bounds (a : ℝ) :
    |a| < Real.sqrt (a ^ 2 + 1e-24) ∧ Real.sqrt (a ^ 2 + 1e-24) ≤ |a| + 1e-12 := by
  have hpos : (0 : ℝ) < 1e-24 := by norm_num
  have hc : (0 : ℝ) ≤ 1e-12 := by norm_num
  have hcc : ((1e-12 : ℝ)) ^ 2 = 1e-24 := by norm_num
  constructor
  · rw [← Real.sqrt_sq_eq_abs]
    exact Real.sqrt_lt_sqrt (sq_nonneg a) (by linarith)
  · calc Real.sqrt (a ^ 2 + 1e-24) ≤ Real.sqrt ((|a| + 1e-12) ^ 2) :=
          Real.sqrt_le_sqrt (by nlinarith [sq_abs a, abs_nonneg a])
      _ = |a| + 1e-12 := Real.sqrt_sq (add_nonneg (abs_nonneg a) hc)

/-- **Certificate (quantitative equivalence)**: if the friction law `y ↦ g(0,p,y,0)` is `L`-Lipschitz, then under uniform slip, at every node
the difference between the friction coefficient and the constant `μ = g(0,p,|V·τ|,0)` is at most `L·10⁻¹²`, independent of the mesh. -/
theorem frictionValue_uniform_2d_close {ι : Type*} (S : Finset ι) (test JxW coord : ι → ℝ)
    (vrel t0 t1 : ι → Vec3) (V τ : Vec3)
    (hv : ∀ k ∈ S, vrel k = V) (ht : ∀ k ∈ S, t0 k = τ)
    (hN : normalization S test JxW coord ≠ 0) (mu : ℝ) (g : ℝ → ℝ → ℝ → ℝ → ℝ) (p : ℝ)
    (L : ℝ) (hL : 0 ≤ L) (hg : ∀ y y' : ℝ, |g 0 p y 0 - g 0 p y' 0| ≤ L * |y - y'|) :
    abs (frictionValue true mu g p (nodalSlipRates false S test JxW coord vrel t0 t1).1
        (nodalSlipRates false S test JxW coord vrel t0 t1).2 - g 0 p (abs (dot3 V τ)) 0) ≤ L * 1e-12 := by
  rw [frictionValue_uniform_2d S test JxW coord vrel t0 t1 V τ hv ht hN mu g p]
  obtain ⟨h1, h2⟩ := sqrt_sq_add_eps_bounds (dot3 V τ)
  have hb : abs (Real.sqrt (dot3 V τ ^ 2 + 1e-24) - abs (dot3 V τ)) ≤ 1e-12 := by
    rw [abs_of_nonneg (by linarith)]
    linarith
  calc abs (g 0 p (Real.sqrt (dot3 V τ ^ 2 + 1e-24)) 0 - g 0 p (abs (dot3 V τ)) 0)
      ≤ L * abs (Real.sqrt (dot3 V τ ^ 2 + 1e-24) - abs (dot3 V τ)) := hg _ _
    _ ≤ L * 1e-12 := mul_le_mul_of_nonneg_left hb hL

/-- **Counterexample (to "strictly equivalent to constant μ")**: if the friction law is strictly monotone in the slip rate, then under uniform slip the function friction coefficient
is never strictly equal to `g(0,p,|V·τ|,0)`; the difference comes from the `1e-24` regularization (magnitude given in the previous item, ≤ L·10⁻¹²). -/
theorem frictionValue_uniform_2d_ne_const {ι : Type*} (S : Finset ι) (test JxW coord : ι → ℝ)
    (vrel t0 t1 : ι → Vec3) (V τ : Vec3)
    (hv : ∀ k ∈ S, vrel k = V) (ht : ∀ k ∈ S, t0 k = τ)
    (hN : normalization S test JxW coord ≠ 0) (mu : ℝ) (g : ℝ → ℝ → ℝ → ℝ → ℝ) (p : ℝ)
    (hmono : StrictMono fun y => g 0 p y 0) :
    frictionValue true mu g p (nodalSlipRates false S test JxW coord vrel t0 t1).1
        (nodalSlipRates false S test JxW coord vrel t0 t1).2 ≠ g 0 p |dot3 V τ| 0 := by
  rw [frictionValue_uniform_2d S test JxW coord vrel t0 t1 V τ hv ht hN mu g p]
  intro h
  have hlt : g 0 p |dot3 V τ| 0 < g 0 p (Real.sqrt (dot3 V τ ^ 2 + 1e-24)) 0 :=
    hmono (sqrt_sq_add_eps_bounds (dot3 V τ)).1
  linarith

/-! ## §6 Counterexample: the old implementation is amplified by the number of quadrature points / mortar segments -/

/-- Under a constant field the old quantity = `−a · ∑ψ_k` (no JxW, no normalization). -/
theorem oldRealTangentialVelocity_const {ι : Type*} (S : Finset ι) (test : ι → ℝ)
    (vrel tangent : ι → Vec3) (a : ℝ) (ha : ∀ k ∈ S, dot3 (vrel k) (tangent k) = a) :
    oldRealTangentialVelocity S test vrel tangent = -a * ∑ k ∈ S, test k := by
  unfold oldRealTangentialVelocity
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro k hk
  rw [dot3_smul3]
  have : dot3 (fun c => -vrel k c) (tangent k) = -dot3 (vrel k) (tangent k) := by
    simp only [dot3]; ring
  rw [this, ha k hk]
  ring

/-- Values of the test function at the quadrature points for an interior node of a 1D interface: one linear element on each side, each with 2-point Gauss quadrature
(ξ = 1/2 ± √3/6); the linear basis at the node takes values 1/2 ± √3/6 at these 4 quadrature points (twice each; order does not affect the sum). -/
def twoElemTest (n : ℕ) : ℝ := 1 / 2 + (-1) ^ n * (Real.sqrt 3 / 6)

theorem twoElemTest_sum : ∑ n ∈ Finset.range 4, twoElemTest n = 2 := by
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, twoElemTest]
  ring

/-- **Counterexample (old implementation)**: for the interior node above (element length `h > 0`, quadrature weight `h/2`, `coord = 1`), uniform slip
`a = V·τ ≠ 0`: the new implementation gives `−a`, the old gives `−2a` (amplification factor = ∑ψ_k at the node = 2, which grows further with the number of mortar segments
and quadrature points). -/
theorem old_vs_new_interior_node (V τ : Vec3) (h : ℝ) (hh : 0 < h) (ha : dot3 V τ ≠ 0) :
    nodalSlipRate (Finset.range 4) twoElemTest (fun _ => h / 2) (fun _ => 1) (fun _ => V)
        (fun _ => τ) = -dot3 V τ ∧
      oldRealTangentialVelocity (Finset.range 4) twoElemTest (fun _ => V) (fun _ => τ) =
        -2 * dot3 V τ ∧
      oldRealTangentialVelocity (Finset.range 4) twoElemTest (fun _ => V) (fun _ => τ) ≠
        nodalSlipRate (Finset.range 4) twoElemTest (fun _ => h / 2) (fun _ => 1) (fun _ => V)
          (fun _ => τ) := by
  have hN : normalization (Finset.range 4) twoElemTest (fun _ => h / 2) (fun _ => 1) = h := by
    show ∑ n ∈ Finset.range 4, twoElemTest n * (h / 2 * 1) = h
    rw [← Finset.sum_mul, twoElemTest_sum]
    ring
  have hnew : nodalSlipRate (Finset.range 4) twoElemTest (fun _ => h / 2) (fun _ => 1)
      (fun _ => V) (fun _ => τ) = -dot3 V τ :=
    nodalSlipRate_uniform _ _ _ _ _ _ V τ (fun _ _ => rfl) (fun _ _ => rfl)
      (by rw [hN]; exact hh.ne')
  have hold : oldRealTangentialVelocity (Finset.range 4) twoElemTest (fun _ => V) (fun _ => τ) =
      -2 * dot3 V τ := by
    rw [oldRealTangentialVelocity_const _ _ _ _ (dot3 V τ) (fun _ _ => rfl), twoElemTest_sum]
    ring
  refine ⟨hnew, hold, ?_⟩
  rw [hnew, hold]
  intro hc
  apply ha
  linarith

end

end PhysAudit.Moose34020
