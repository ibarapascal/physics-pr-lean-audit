# MOOSE#34020 ΔQ1 spot-check implementation log

> Purpose: record the sentence-by-sentence mapping of the three tracks' check plans to `spotcheck.py` / `mortar_port.py`, and all implementation judgments, for review under §6/§7 of the ΔQ1 protocol.
> Structure: §0 environment (non-native) → §1 Track 1 → §2 Track 2 → §3 Track 3 → §4 time and memory → §5 deviation candidates.
> Update rule: after the smoke run only script errors may be fixed; changing the plan mapping must be recorded in protocol §7 and appended to §5.
> References: upstream: the protocol, the #34020 sections of `audit/specs/track1.md`/`track2.md`, `audit/specs/track3_guided/idaholab_moose_34020/spec.md` §2/4/5, the source of merge commit f90f508; downstream: `dq1/outputs/mo34020_t{1,2,3}.json`.

## §0 Environment: non-native recomputation (deviation candidate)

- Compute node A and compute node B have no MOOSE/libMesh/PETSc build (read-only check, same as #33307). No build was made.
- `mortar_port.py` recomputes, line by line from the merge-commit source, all numerical logic changed by the PR:
  - `ComputeDynamicWeightedGapLMMechanicalContact::computeQpProperties/computeQpIProperties`: `_qp_factor = JxW_msm·coord`, `weighted_gap[dof].second += test·_qp_factor` (accumulated unconditionally after the merge);
  - `ComputeDynamicFrictionalForceLMMechanicalContact::computeQpProperties/computeQpIProperties`: `qp_tv = −relative_velocity·JxW·coord`, `relative_velocity = primary side − secondary side` (the nodal velocities on each side interpolated with that side's element shape functions), `weighted_tv[dof][α] += test·qp_tv·t_α(node)`;
  - `timestepSetup` (the dofs owned by this process: s = weighted_tv / weighted_gap.second), `oldNodalSlipRate` (no history → 0), `computeFrictionValue` (F(0, p_old, √(s0²+s1²+1e-24), 0) or the constant mu);
  - constructor parameter validation (the `mu > 0` range check; error if both or neither are given).
- Geometric restriction: flat interfaces; 2D is a line-segment interface (1 tangent), 3D is an axis-aligned rectangular face (2 tangents); tensor-product meshes on both primary and secondary sides, a mortar segment = an interval/rectangle cut out by the union of the breakpoints of both sides, Gauss integration on each segment. The dual basis follows libMesh `compute_dual_shape_coeffs` (dual_coeff = A⁻¹D, A the element mass matrix, D=diag∫φ, ψ_i = Σ_j dual_coeff(j,i)φ_j; for rectangular elements the tensor product of the one-dimensional coefficients).
- Not recomputable: the nonlinear contact solve (normal/tangential LM, |λ_t|/λ_n, displacements), the Newmark velocity update, MPI communication itself (Track 3 replaces it by random partitioning of integration points and summing partial sums), node selection in edge dropping, restep, elemental LM.
- Wherever the plan backs out the slip rate from μ_eff = |λ_t|/λ_n, the actual argument passed to the friction function and its return value are read instead (in the slipping state, when the Coulomb law holds with equality, μ_eff equals this return value).
- The velocity field is given directly (not through Newmark); the slip field is recorded as "secondary velocity − primary velocity" = −v_rel in the code, with the primary side at rest.
- Seed 20261017 + 10×10 + track number.

## §1 Track 1 (`audit/specs/track1.md`, 2D plane problem)

Branches: `use_dual` ∈ {false, true} (k=2).

| Plan text | Code |
|---|---|
| Operating domain: \|V\|∈[1e-4,10] m/s; p∈[1e3,1e9] Pa; 4–64 elements on the interface; Newmark β=0.25, γ=0.5; μ=0.2+0.1·tanh(y/0.1) | `t1_sample`: \|V\| linearly uniform, random sign; p linearly uniform; element count integer-uniform 4–64; Δt, β, γ do not enter the recomputation |
| A: constant vector relative tangential velocity V, every active node s_j = −V·t_j, \|s_j\|=\|V\|, independent of mesh size, element order, quadrature order, and dual or not | `t1_A_s_j==-V_*` (relative velocity = primary − secondary = V) |
| A check: μ_eff compared with g(\|V\|) at step ≥2; read \|s_j\| directly when g(y)=y; relative 1e-6; μ_eff unchanged under 2 mesh refinements | `t1_A_mu_eff==g(|V|)_*`, `t1_A_identity_law_|s_j|==|V|_*`; every point is done for all 12 combinations of refinement ×1/×2/×4, element order 1/2, quadrature points p+1/p+2; `t1_A_mu_eff_mesh/order/quadrature_invariance` |
| B: with the same mesh compare `function_friction=g` with `mu=g(V)`, relative 1e-6 | `t1_B_function_vs_constant_mu_*` (only the friction coefficient is compared, not displacement/LM) |
| C: the upper block rotates uniformly, tangential velocity distributed linearly along the interface, use_dual=true, \|s_j\| inferred from μ_eff compared with the analytic linear field | dual branch: `t1_C_dual_linear_field_*` |
| At the first step and for nodes just entering the interface take s=0, excluded from the check | only the history established by the previous step's kinematics is examined |

Implementation judgments:
1. The interface length is not given: 1 is used (secondary face [0,1]); the primary face is [−0.5,1.5] with the same element count as the secondary face (element size twice, non-matching), imitating the different widths of the upper and lower blocks of the fixture.
2. Linear field: v(x) = V_lo + (V−V_lo)·x, V_lo = 1e-4·sign(V), so that the velocity over the whole interface falls in the operating domain [1e-4, \|V\|] (if V·x were taken, the velocity would be 0 at x=0, and the regularisation term of √(·+1e-24) would put the default tolerance of abs 1e-12 exactly at the boundary). C has no tolerance, so the default rel 1e-6, abs 1e-12 is used.
3. Invariant A states independence of "element order and quadrature order", while the check method mentions only mesh refinement; per protocol §2 they are included: order {1,2}, quadrature points per direction {p+1, p+2}.
4. "Rotation about the normal" does not lie in the plane in 2D, so, following the plan's parenthetical remark, a tangential velocity distributed linearly along the interface is used.

## §2 Track 2 (`audit/specs/track2.md`)

Branches (k=8): {2D, 3D} × {standard, dual} × {first order, second order}.

| Plan text | Code |
|---|---|
| Operating domain: slip velocity 1e-3–10 m/s; μ₀ 0.1–0.8; h at three refinement levels; first/second-order elements; 2D/3D; standard/dual | `t2_sample`: velocity linearly uniform, random sign in 2D, direction angle U[0,2π) in 3D; μ₀ linearly uniform |
| Template: `frictional-mortar-3d-dynamics-light-function` (upper block 0.5 wide, lower block 1 wide, 2×2 elements each) | secondary face [−0.25,0.25]^d, primary face [−0.5,0.5]^d, base 2 elements/direction, three refinement levels = 2, 4, 8 elements/direction |
| I1: constant tangential velocity v, `function_friction`=0.5−0.1·v, \|λ_t\|/λ_n compared with 0.5−0.1\|v\| rtol 1e-6; 3 mesh levels × 2 quadrature orders do not change it | `t2_I1_s_i==v.t_i_*`, `t2_I1_mu==0.5-0.1|v|_*` (3 levels × quadrature points {p+1,p+2} at every point) |
| I2: use_dual=true, the upper block rotates uniformly about the normal, the v inferred from nodal μ compared with ω r_i, relative 1e-8 | dual branch: `t2_I2_dual_rotation_s_i==omega_x_r_*` (components) and `t2_I2_dual_rotation_mu_inferred_speed_*` (y inferred from μ=0.5−0.1y) |
| I3: μ(v)=μ₀−b·v, when μ₀−b·v_max>0 the friction coefficient is in [μ₀−b·v_max, μ₀] and not negative; at the first step/for newly entering nodes μ=μ(p_old,0) | `t2_I3_mu_in_[mu0-b*vmax,mu0]_and_>0_*` (points not satisfying the condition are not checked), `t2_I3_first_step_mu==mu(p_old,0)` |

Implementation judgments:
1. The b of I3 is not given: 0.1 of the I1 law is used; v_max is the velocity drawn at that point.
2. The old pressure is not given in the operating domain, and neither friction law depends on pressure: p_old=1 is fixed.
3. Rotation: in 3D ω ẑ×r (centre at the midpoint of the interface), ω = v/(L/√2) so that the corner-point velocity equals the drawn v; 2D uses the linear field ω·x, ω = v/(L/2). The tolerance is relative 1e-8; the plan gives no absolute tolerance, so the protocol's 1e-12 is used. The theoretical slip at the centre node is 0, and the velocity inferred from μ is at least √1e-24 = 1e-12, exactly at the boundary of abs 1e-12 (rejections may appear; this is a regularisation convention, for ruling).
4. Slip direction convention: the plan's v is taken as the "−v_rel" (secondary minus primary) of the PR description, so the expected s_i = v·t_i.
5. Second-order dual basis: ∫φ is positive for both EDGE3/QUAD9, so libMesh does not take the QUAD8/TRI6 correction branch; computed as A⁻¹D.

## §3 Track 3 (spec.md §2/§4/§5)

Branches (k=16): {2D,3D} × normalize_c{false,true} × friction_projection_degree{ONE,TWO} × correct_edge_dropping{false,true} (the "core matrix" of plan §5).

| Plan text | Code |
|---|---|
| Planar patch of length 0.1–10, 1–8 elements per tangential direction, matching and non-matching; velocity components finite and signed (list 0…1 with both signs); pressure 0…1000; dual multipliers; "sufficient" quadrature order | `t3_sample`: L linearly uniform; the secondary/primary element counts are each integer-uniform 1–8 (equal means matching); components U[−1,1]; p_old U[0,1000]; quadrature points {2,3}; first-order elements |
| I1: W_i s_i = V_i | `t3_I1_W*s==V` |
| I1: uniform relative velocity → s=b, independent of support length ×0.1/1/10, refinement 1/2/4/8, non-matching primary-side subdivision, quadrature order; tolerance abs 1e-12 v* + rel 1e-9 | `t3_I1_constant_reproduction_s==b`, `t3_I1_support_scale_*`, `t3_I1_refine_x*` |
| I1: affine non-uniform slip field, independently integrate V_i/W_i with dual test functions and the mortar segment measure, compare with the norm after projection; including opposing slip within the support | `t3_I1_affine_…`, `t3_I1_opposing_…`: independent reference `oracle_affine` (high-order Gauss integration of the analytic field directly on the secondary elements, without primary-side interpolation or segmentation), compared against √(\|s\|²+1e-24) |
| I2: ≥3 accepted steps, accelerating slip and varying normal load; the callback argument at step n+1 = the projection of step n; at the first step velocity 1e-12, old pressure, t and z are 0; newly entering dofs take zero slip | `t3_I2_first_step_args_…`, `t3_I2_step{2,3}_speed_from_step…`, `t3_I2_step*_old_pressure` (tolerance 1e-9 p* + 1e-7\|p\|), `t3_I2_entering_dof_zero_slip_convention` (the previous step's primary face covers only the first half of the elements) |
| I2: `--test-restep` | cannot be recomputed, `not_covered` |
| I3: partition sum agrees with serial 1e-10 v* + 1e-8\|s\|; normalize_c, projection degree, edge dropping do not change the constitutive arguments | `t3_I3_partition_{2,4}_vs_serial` (integration points randomly assigned to 2/4 "processes"); the three switches do not enter s in the recomputation, so one trivially true check is recorded and explained in not_covered |
| I4: 2D ±b; 3D (b,0),(0,b),(−b,0),(0,−b),(3b/5,4b/5) and the reverses; geometry and field rotated rigidly together; speeds equal to 1e-10 relative + I1 absolute | `t3_I4_speed_isotropic_over_directions`, `t3_I4_rigid_rotation_speed_unchanged` (3D, the tangent basis rotates with it) |
| I5: constant mu=0.3 and a ParsedFunction that is constantly 0.3 give the same result; mu=0, −0.1, neither given, both given must raise errors, only a positive constant succeeds | `t3_I5_constant_mu_vs_constant_function`; mandatory points `valid_*` |
| I5: projection identity, disc bound, dissipation, epsilon threshold, merge/base regression | needs an LM solve, `not_covered` |
| Velocity component list 0,1e-12,1e-8,1e-4,0.01,0.02,0.1,1 (both signs), pressure 0, ε/2, ε, 2ε, 1, 10, 100, 1000 (ε=1e-7) | mandatory points `2D_b=…`, `3D_b=(…)` (96 cases, 8 pressures per case): constant reproduction, argument (0,p,velocity,0) exact, the given law positive for velocity ≤0.02 |
| wear_depth zero and nonzero | in the recomputation W contains only `_qp_factor` and is independent of wear; not listed separately (trivially true) |

Implementation judgments:
1. The velocity components are given only as a discrete list (0…1, both signs): the sampling interval is [−1, 1] and the list values are mandatory points.
2. Primary and secondary faces have the same outer frame (full support) and only the primary-face subdivision is changed; partial coverage is used only in the I2 "newly entering" check, with the coverage boundary aligned to the secondary-element boundary (ensuring W_i>0 at every node).
3. Each component of the affine-field gradient G is U[−1,1]/L; the opposing slip field is b·(x−L/2)/(L/2).
4. "Use prescribed kinematics": the recomputation is already given kinematics and involves no nonlinear trajectory.
5. First-order elements (fixture HEX8); dual multipliers (the plan calls this the established configuration).

## §4 Expected time and peak memory (speculation, not run)

| Track | Per point | 10⁴ points | Basis |
|---|---|---|---|
| t1 | ~20–50 ms (12 combinations, at most 256 elements × order 2) | ~5–8 min | numpy vectorised accumulation + per-node Python friction calls |
| t2 | ~10–30 ms (3D second order 8×8 elements is the largest) | ~3–5 min | |
| t3 | ~50–150 ms (about 23 accumulations, about 47 thousand integration points for 3D refinement ×8) | ~10–25 min | `np.add.at` |

Peak memory <300 MB; the JobSpec declares 1 GiB.

## §5 Deviation candidates

- **Non-native**: the nodes have no MOOSE build; the changed logic is recomputed in Python (coverage list in §0). μ_eff=|λ_t|/λ_n is always replaced by the argument/return value of the friction function; checks needing a nonlinear solve, restep, real MPI or edge dropping were not tested.
