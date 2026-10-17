#!/usr/bin/env python3
"""Moose34020_replay.py - differential check for idaholab/moose PR #34020 (merge f90f508580).

Nodal slip rate passed to the dynamic mortar `function_friction`
(ComputeDynamicFrictionalForceLMMechanicalContact.C:124,141-148,183-186,374-377 and
ComputeDynamicWeightedGapLMMechanicalContact.C:192,223).

There is no native MOOSE build available, so this script does NOT run MOOSE. It compares
  (A) an independent float64 transcription of the C++ accumulation loops
      (computeQpProperties / computeQpIProperties / timestepSetup / computeFrictionValue), driven by
      synthetic mortar segments, Gauss rules, standard and dual linear (2D) / bilinear (3D) bases,
      planar and axisymmetric (coord = 2*pi*r) weights,
  (B) the formulas of the Lean definitions in lean/PhysAudit/Moose34020.lean
      (`nodalSlipRate` = sum(dot3(smul3(test, qpTangentialVelocityNodal), t)) / sum(test*qpFactor)),
      evaluated in exact rational arithmetic (fractions.Fraction) on the same float inputs,
and checks every theorem of the Lean module numerically.

License:      LGPL-2.1 (this file contains a line-by-line reimplementation of MOOSE C++ source,
              which is licensed under LGPL-2.1; it is not covered by the MIT license of this
              archive; see THIRD_PARTY.md).
Dependencies: Python >= 3.8 standard library only (math, random, fractions, json).
Run:          python3 -I Moose34020_replay.py > Moose34020_replay.json
Output:       one JSON object on stdout.
"""
import json
import math
import random
from fractions import Fraction

random.seed(34020)
N_CASES = 2000          # random discretizations (float transcription)
N_EXACT = 200           # of which checked against the exact-rational Lean formula

# ---------------------------------------------------------------- quadrature / bases
GAUSS = {
    1: ([0.0], [2.0]),
    2: ([-1 / math.sqrt(3), 1 / math.sqrt(3)], [1.0, 1.0]),
    3: ([-math.sqrt(0.6), 0.0, math.sqrt(0.6)], [5 / 9, 8 / 9, 5 / 9]),
    4: ([-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526],
        [0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538]),
}


def lin_basis(xi):            # xi in [0,1] on the secondary element
    return [1.0 - xi, xi]


def dual_basis(xi):           # biorthogonal to lin_basis on [0,1]: psi1 = 2 phi1 - phi2
    p = lin_basis(xi)
    return [2 * p[0] - p[1], 2 * p[1] - p[0]]


def dot3(u, v):
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


# ---------------------------------------------------------------- (A) C++ transcription
def accumulate_2d(nodes, vel_fun, tangent, basis, nqp, max_seg, axisym, rng):
    """Emulate the mortar loop on a 1D secondary interface (2D problem).

    Returns per-node dicts: weighted tangential velocity, normalization, old (parent) quantity,
    and the raw contribution list (test, JxW, coord, vrel, t) for the exact Lean formula."""
    wtv, norm, old, contrib = {}, {}, {}, {}
    for e in range(len(nodes) - 1):
        xa, xb = nodes[e], nodes[e + 1]
        # mortar segments: the secondary element is cut by projected primary nodes
        cuts = sorted(rng.random() for _ in range(rng.randint(0, max_seg - 1)))
        segs = [0.0] + cuts + [1.0]
        for s in range(len(segs) - 1):
            ua, ub = segs[s], segs[s + 1]
            qx, qw = GAUSS[nqp]
            for x_ref, w_ref in zip(qx, qw):
                xi = ua + (ub - ua) * (x_ref + 1) / 2          # position on secondary element
                JxW = w_ref * (ub - ua) / 2 * (xb - xa)        # _JxW_msm[_qp]
                xphys = xa + xi * (xb - xa)
                coord = 2 * math.pi * (xphys + 0.5) if axisym else 1.0   # _coord[_qp]
                vrel = vel_fun(xphys)                           # _relative_velocity
                # :124  _qp_tangential_velocity_nodal = -_relative_velocity * (JxW * coord)
                qtvn = [-c * (JxW * coord) for c in vrel]
                # weighted gap :192  _qp_factor = JxW * coord
                qp_factor = JxW * coord
                test = basis(xi)
                for i in range(2):
                    dof = e + i
                    # :141-142  wtv[dof][0] += _test * _qp_tangential_velocity_nodal * t0
                    wtv[dof] = wtv.get(dof, 0.0) + dot3([test[i] * c for c in qtvn], tangent)
                    # weighted gap :223  .second += _test * _qp_factor
                    norm[dof] = norm.get(dof, 0.0) + test[i] * qp_factor
                    # parent commit :144-146  old += _test * raw(-_relative_velocity) * t0
                    old[dof] = old.get(dof, 0.0) + dot3([test[i] * (-c) for c in vrel], tangent)
                    contrib.setdefault(dof, []).append((test[i], JxW, coord, vrel, tangent))
    return wtv, norm, old, contrib


def cpp_slip(wtv, norm):
    # timestepSetup :183-186
    return {d: wtv[d] / norm[d] for d in wtv}


def cpp_friction_value(g, p, v0, v1):
    # computeFrictionValue :374-377
    mag = math.sqrt(v0 * v0 + v1 * v1 + 1.0e-24)
    return g(0.0, p, mag, 0.0)


# ---------------------------------------------------------------- (B) Lean definition, exact
def lean_nodal_slip_rate_exact(contribs):
    """nodalSlipRate S test JxW coord vrel tangent, exact rational arithmetic."""
    num = Fraction(0)
    den = Fraction(0)
    for test, JxW, coord, vrel, t in contribs:
        T, J, C = Fraction(test), Fraction(JxW), Fraction(coord)
        qtvn = [-Fraction(c) * (J * C) for c in vrel]          # qpTangentialVelocityNodal
        sm = [T * c for c in qtvn]                             # smul3
        num += sum(sm[k] * Fraction(t[k]) for k in range(3))   # dot3
        den += T * (J * C)                                     # test * qpFactor
    return num / den, den


# ---------------------------------------------------------------- friction laws
def g_test(t, x, y, z):          # mu_function of frictional-mortar-3d-dynamics-light-function.i
    return 0.3 + 0.5 * 2.17 ** (-x / 100.0) - 10.0 * y


def g_tanh(t, x, y, z):          # monotone velocity-strengthening law, Lipschitz L = 1
    return 0.2 + 0.1 * math.tanh(y / 0.1)


def relerr(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# ---------------------------------------------------------------- main checks
res = {
    "module": "Moose34020",
    "merge_commit": "f90f508580f7611f7d75ee458cd177affc7288d1",
    "native_replay": "SKIP: no MOOSE build available; float64 transcription of the C++ loops is "
                     "compared with the Lean formulas instead",
}

max_rel_cpp_vs_lean = 0.0
exact_const_ok = 0
exact_checked = 0
max_rel_const = 0.0
max_rel_old_over_new = 0.0
min_rel_old_over_new = float("inf")
max_mu_diff_vs_const = 0.0
dual_max_err = 0.0
std_lin_max_err = 0.0
cases_axisym = 0

for case in range(N_CASES):
    rng = random.Random(case)
    n_el = rng.randint(4, 16)
    h = [rng.uniform(0.5, 2.0) for _ in range(n_el)]
    nodes = [0.0]
    for hh in h:
        nodes.append(nodes[-1] + hh * 0.01)
    V = rng.choice([-1, 1]) * 10 ** rng.uniform(-4, 1)          # |V| in [1e-4, 10] m/s
    tangent = [1.0, 0.0, 0.0]
    basis = rng.choice([lin_basis, dual_basis])
    nqp = rng.randint(2, 4)
    max_seg = rng.randint(1, 4)
    axisym = rng.random() < 0.3
    cases_axisym += axisym
    wtv, norm, old, contrib = accumulate_2d(nodes, lambda x: [V, 0.0, 0.0], tangent, basis, nqp,
                                            max_seg, axisym, rng)
    slip = cpp_slip(wtv, norm)
    a = dot3([V, 0.0, 0.0], tangent)
    for d in slip:
        # Theorem nodalSlipRate_const_field / nodalSlipRate_uniform: slip == -a
        max_rel_const = max(max_rel_const, relerr(slip[d], -a))
        # Theorem old_vs_new: old == -a * sum(test)
        ratio = old[d] / slip[d]
        max_rel_old_over_new = max(max_rel_old_over_new, ratio)
        min_rel_old_over_new = min(min_rel_old_over_new, ratio)
        # frictionValue_uniform_2d / _close
        mu = cpp_friction_value(g_tanh, 1e6, slip[d], 0.0)
        max_mu_diff_vs_const = max(max_mu_diff_vs_const, abs(mu - g_tanh(0, 1e6, abs(a), 0)))
        if case < N_EXACT:
            lean_val, den = lean_nodal_slip_rate_exact(contrib[d])
            max_rel_cpp_vs_lean = max(max_rel_cpp_vs_lean, relerr(slip[d], float(lean_val)))
            exact_checked += 1
            exact_const_ok += int(den != 0 and lean_val == -Fraction(a))

    # nodalSlipRate_dual_exact: linear field, dual basis, planar, nqp >= 2 (exact for degree 2)
    if not axisym:
        c0, c1 = rng.uniform(-1, 1), rng.uniform(-50, 50)
        lin = lambda x: [c0 + c1 * x, 0.0, 0.0]
        for bas, tag in ((dual_basis, "dual"), (lin_basis, "std")):
            wtv2, norm2, _, _ = accumulate_2d(nodes, lin, tangent, bas, nqp, max_seg, False,
                                              random.Random(case + 10 ** 6))
            s2 = cpp_slip(wtv2, norm2)
            scale = max(abs(c0 + c1 * x) for x in nodes) + 1e-12
            err = max(abs(s2[d] + (c0 + c1 * nodes[d])) / scale for d in s2)
            if tag == "dual":
                dual_max_err = max(dual_max_err, err)
            else:
                std_lin_max_err = max(std_lin_max_err, err)

res["cases"] = N_CASES
res["cases_axisymmetric"] = cases_axisym
res["diff_cpp_float_vs_lean_exact"] = {"nodes_checked": exact_checked,
                                       "max_rel_diff": max_rel_cpp_vs_lean}
res["thm_nodalSlipRate_const_field"] = {
    "exact_rational_equal_count": exact_const_ok, "exact_rational_checked": exact_checked,
    "float_max_rel_err_vs_minus_a": max_rel_const, "holds": exact_const_ok == exact_checked}
res["thm_nodalSlipRate_mesh_independent"] = {
    "note": "every random discretization (4-16 elements, 1-4 mortar segments per element, "
            "2-4 Gauss points, standard/dual basis, planar/axisymmetric) gives -a",
    "float_max_rel_err": max_rel_const}
res["thm_nodalSlipRate_dual_exact"] = {
    "dual_basis_linear_field_max_rel_err": dual_max_err,
    "standard_basis_linear_field_max_rel_err": std_lin_max_err}

# 3D: bilinear quad face, random orthonormal tangents, tangential V
def rand_orthonormal(rng):
    n = [rng.gauss(0, 1) for _ in range(3)]
    nn = math.sqrt(dot3(n, n))
    n = [c / nn for c in n]
    a = [1.0, 0.0, 0.0] if abs(n[0]) < 0.9 else [0.0, 1.0, 0.0]
    t0 = [a[1] * n[2] - a[2] * n[1], a[2] * n[0] - a[0] * n[2], a[0] * n[1] - a[1] * n[0]]
    l0 = math.sqrt(dot3(t0, t0))
    t0 = [c / l0 for c in t0]
    t1 = [n[1] * t0[2] - n[2] * t0[1], n[2] * t0[0] - n[0] * t0[2], n[0] * t0[1] - n[1] * t0[0]]
    return t0, t1


max_rel_3d = 0.0
for case in range(500):
    rng = random.Random(7 * case + 1)
    t0, t1 = rand_orthonormal(rng)
    al, be = rng.uniform(-5, 5), rng.uniform(-5, 5)
    V = [al * t0[k] + be * t1[k] for k in range(3)]
    nq = rng.randint(2, 4)
    qx, qw = GAUSS[nq]
    w0 = w1 = nrm = 0.0
    # one node shared by 4 bilinear elements, each split into random mortar sub-rectangles
    for e in range(4):
        hx, hy = rng.uniform(0.5, 2), rng.uniform(0.5, 2)
        for _ in range(rng.randint(1, 3)):
            u0, u1 = sorted((rng.random(), rng.random()))
            for xa, wa in zip(qx, qw):
                for xb, wb in zip(qx, qw):
                    u = u0 + (u1 - u0) * (xa + 1) / 2
                    v = (xb + 1) / 2
                    JxW = wa * wb * (u1 - u0) / 2 * 0.5 * hx * hy
                    test = (1 - u) * (1 - v)            # node at local corner (0,0)
                    q = [-c * JxW for c in V]
                    w0 += dot3([test * c for c in q], t0)
                    w1 += dot3([test * c for c in q], t1)
                    nrm += test * JxW
    s0, s1 = w0 / nrm, w1 / nrm
    mu = cpp_friction_value(g_test, 0.0, s0, s1)
    expect = g_test(0.0, 0.0, math.sqrt(dot3(V, V) + 1e-24), 0.0)
    max_rel_3d = max(max_rel_3d, relerr(mu, expect))
res["thm_frictionValue_uniform_3d"] = {"cases": 500, "max_rel_err": max_rel_3d}

# frictionValue_uniform_2d_close and _ne_const: effect of the 1e-24 regularization in float64
absorb_threshold = None
for k in range(0, 2000):
    a = 10 ** (-12 + k * 0.005)
    if math.sqrt(a * a + 1e-24) == abs(a):
        absorb_threshold = a
        break
res["thm_frictionValue_uniform_2d_close"] = {
    "law": "0.2+0.1*tanh(y/0.1), L=1", "max_abs_diff_vs_constant_mu": max_mu_diff_vs_const,
    "bound_L_times_1e-12": 1e-12}
res["thm_frictionValue_uniform_2d_ne_const"] = {
    "real_number_statement": "sqrt(a^2+1e-24) > |a| strictly, so a strictly monotone law never "
                             "equals g(|a|) exactly",
    "float64_smallest_|a|_where_regularization_is_absorbed": absorb_threshold,
    "example_a_1e-6": {"sqrt": math.sqrt(1e-12 + 1e-24), "abs": 1e-6},
}
res["thm_frictionValue_no_history"] = {"slip_magnitude_used": math.sqrt(1e-24)}

# old_vs_new_interior_node (exact 2-element, 2-point Gauss configuration)
gp = [0.5 - math.sqrt(3) / 6, 0.5 + math.sqrt(3) / 6]
tests = [gp[0], gp[1], gp[1], gp[0]]  # right-end basis on elem 1, left-end basis on elem 2
Vx = 0.5
new = sum(t * (-Vx) * 0.005 for t in tests) / sum(t * 0.005 for t in tests)
oldv = sum(t * (-Vx) for t in tests)
res["thm_old_vs_new_interior_node"] = {
    "V": Vx, "new": new, "old": oldv, "old_over_new": oldv / new,
    "random_cases_old_over_new_range": [min_rel_old_over_new, max_rel_old_over_new],
    "test_input_law_example": {
        "law": "0.3 + 0.5*2.17^(-p/100) - 10*y, p=0",
        "mu_with_new_slip_V0.05": cpp_friction_value(g_test, 0.0, -0.05, 0.0),
        "mu_with_old_slip_4_segments_V0.05": cpp_friction_value(g_test, 0.0, -0.05 * 4, 0.0),
    },
}
print(json.dumps(res, indent=2))
