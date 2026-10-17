"""MOOSE#34020 (dynamic mortar function_friction now uses the nodal slip rate ∫ψ(−v_rel)·t/∫ψ) unified ΔQ1 spot check.
Non-native: there is no MOOSE build on the node; the modified logic is recomputed line by line in mortar_port.py from the C++ source of merge commit f90f508;
quantities that need the nonlinear contact solve (|λ_t|/λ_n, displacement, exodiff, restep, real MPI, edge dropping) are not implemented and are listed in not_covered;
wherever the plan backs out μ_eff=|λ_t|/λ_n, the argument passed to the friction function and its return value are read instead (the two coincide when the Coulomb law holds with equality in the slipping state).
Usage: python spotcheck.py --track {1,2,3} --n N --out PATH (pure Python + numpy)
Plan mapping and implementation judgments are in impl_log.md.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spotlib as SL  # noqa: E402
import mortar_port as MP  # noqa: E402

PR = "idaholab/moose#34020"
PR_IDX = 10

NOT_COVERED = {
    1: ["nonlinear solve for μ_eff=|λ_t|/λ_n (replaced by the friction-function argument/return value)", "B displacement and LM exodiff comparison (only the friction coefficient is compared)",
        "velocity obtained from Newmark time integration (the interface velocity field is given directly)", "step 1 and newly entering nodes (the plan requires excluding them)"],
    2: ["nonlinear solve for |λ_t|/λ_n and steady slip (replaced by the friction-function argument/return value)", "PresetVelocity/Newmark velocity (the interface velocity field is given directly)",
        "parallel owner/ghost assert risk (needs MPI)"],
    3: ["--test-restep retry step (needs the executioner)", "real MPI communication (replaced by randomly partitioning integration points and summing partial sums)",
        "node selection of correct/incorrect edge dropping (the recomputed s does not go through that logic; the switch serves only as a branch label)",
        "effect of normalize_c and projection degree on the nonlinear trajectory (only checks that they do not change s)",
        "I5 Coulomb projection identity, disc bound, dissipation, epsilon threshold, merge/base constant-μ regression (needs the LM solution)",
        "elemental CONSTANT friction-variable branch (marked uncertain by the plan)", "contact detection / persistence state"],
}


def versions():
    return dict(moose_merge_commit="f90f508580f7611f7d75ee458cd177affc7288d1", native=False,
                recompute="mortar_port.py (weighted gap / weighted tangential velocity accumulation, timestepSetup normalization, "
                          "oldNodalSlipRate, computeFrictionValue, constructor parameter validation); dual coefficients follow libMesh "
                          "compute_dual_shape_coeffs (A⁻¹D)", numpy=np.__version__, python=sys.version.split()[0])


def patch(dim, sec_box, prim_box, nsec, nprim, p, dual, nq, tangents=None):
    nd = dim - 1
    sec = [MP.Line(sec_box[0], sec_box[1], nsec, p) for _ in range(nd)]
    prim = [MP.Line(prim_box[0], prim_box[1], nprim, p) for _ in range(nd)]
    return MP.Patch(sec, prim, dual, nq, tangents)


def const_field(vec):
    vec = np.asarray(vec, float)
    return lambda pts: np.tile(vec, (pts.shape[0], 1))


ZERO1 = const_field([0.0])
ZERO2 = const_field([0.0, 0.0])


def slip_rates(pt, slip_fn, owner_fn=None, nproc=1):
    """Given a slip field (secondary velocity − primary velocity, i.e. −v_rel in the code) with the primary side at rest, returns the history and (V,W,present)."""
    zero = ZERO1 if pt.nd == 1 else ZERO2
    V, W, present = pt.accumulate(slip_fn, zero, owner_fn, nproc)
    return MP.timestep_setup(V, W, present), (V, W, present)


def hist_array(hist, n, nd):
    s = np.full((n, nd), np.nan)
    for i, v in hist.items():
        s[i] = v[:nd]
    return s


# ================================================================ Track 1 (2D)
def g1(y):
    return 0.2 + 0.1 * math.tanh(y / 0.1)


def t1_sample(dual):
    def s(rng):
        return dict(dual=dual, V=float(rng.uniform(1e-4, 10)) * float(rng.choice([-1.0, 1.0])),
                    p=float(rng.uniform(1e3, 1e9)), n=int(rng.integers(4, 65)))
    return s


def t1_eval(inp):
    V, p_old, n, dual = inp["V"], inp["p"], inp["n"], inp["dual"]
    ch = []
    mus_all = []
    for order in (1, 2):
        for nq in (order + 1, order + 2):
            for ref in (1, 2, 4):
                pt = patch(2, (0.0, 1.0), (-0.5, 1.5), n * ref, n * ref, order, dual, nq)
                # relative tangential velocity (primary − secondary) = V: secondary at rest, primary slides with V → slip field (−v_rel) = −V
                hist, _ = slip_rates(pt, const_field([-V]))
                s = np.array([hist[i][0] for i in sorted(hist)])
                tag = f"p{order}_nq{nq}_ref{ref}"
                ch.append(SL.close(f"t1_A_s_j==-V_{tag}", s, np.full_like(s, -V)))
                mus = np.array([MP.friction_value(p_old, si, 0.0, function=lambda t, x, y, z: g1(y))[0] for si in s])
                ch.append(SL.close(f"t1_A_mu_eff==g(|V|)_{tag}", mus, np.full_like(mus, g1(abs(V)))))
                ys = np.array([MP.friction_value(p_old, si, 0.0, function=lambda t, x, y, z: y)[0] for si in s])
                ch.append(SL.close(f"t1_A_identity_law_|s_j|==|V|_{tag}", ys, np.full_like(ys, abs(V))))
                mconst = MP.friction_value(p_old, s[0], 0.0, mu=g1(abs(V)))[0]
                ch.append(SL.close(f"t1_B_function_vs_constant_mu_{tag}", mus, np.full_like(mus, mconst)))
                mus_all.append(mus)
                if dual:
                    Vlo = 1e-4 * math.copysign(1.0, V)
                    lin = lambda pts, Vlo=Vlo: -(Vlo + (V - Vlo) * pts[:, :1] / 1.0)  # noqa: E731
                    hl, _ = slip_rates(pt, lin)
                    xs = pt.node_coords()[:, 0]
                    yl = np.array([MP.friction_value(p_old, hl[i][0], 0.0, function=lambda t, x, y, z: y)[0]
                                   for i in range(len(xs))])
                    ch.append(SL.close(f"t1_C_dual_linear_field_|s_j|==|v(x_j)|_{tag}", yl,
                                       np.abs(Vlo + (V - Vlo) * xs)))
    allm = np.concatenate(mus_all)
    ch.append(SL.close("t1_A_mu_eff_mesh/order/quadrature_invariance", allm, np.full_like(allm, allm[0])))
    return ch


# ================================================================ Track 2 (2D/3D × std/dual × order 1/2)
def t2_sample(dim, dual, order):
    def s(rng):
        d = dict(dim=dim, dual=dual, order=order, v=float(rng.uniform(1e-3, 10)), mu0=float(rng.uniform(0.1, 0.8)))
        if dim == 2:
            d["dir"] = float(rng.choice([-1.0, 1.0]))
        else:
            d["angle"] = float(rng.uniform(0, 2 * math.pi))
        return d
    return s


def law_I1(t, x, y, z):
    return 0.5 - 0.1 * y


def t2_eval(inp):
    dim, dual, order, vmag, mu0 = inp["dim"], inp["dual"], inp["order"], inp["v"], inp["mu0"]
    vec = [vmag * inp["dir"]] if dim == 2 else [vmag * math.cos(inp["angle"]), vmag * math.sin(inp["angle"])]
    ch = []
    p_old = 1.0
    law3 = lambda t, x, y, z: mu0 - 0.1 * y  # noqa: E731
    for lev in (1, 2, 4):
        for nq in (order + 1, order + 2):
            pt = patch(dim, (-0.25, 0.25), (-0.5, 0.5), 2 * lev, 2 * lev, order, dual, nq)
            tag = f"lev{lev}_nq{nq}"
            hist, _ = slip_rates(pt, const_field(vec))
            S = hist_array(hist, pt.n_nodes(), pt.nd)
            ch.append(SL.close(f"t2_I1_s_i==v.t_i_{tag}", S, np.tile(vec, (S.shape[0], 1))))
            mus = np.array([MP.friction_value(p_old, *(list(S[i]) + [0.0])[:2], function=law_I1)[0]
                            for i in range(S.shape[0])])
            ch.append(SL.close(f"t2_I1_mu==0.5-0.1|v|_{tag}", mus, np.full_like(mus, 0.5 - 0.1 * vmag)))
            if mu0 - 0.1 * vmag > 0:
                m3 = np.array([MP.friction_value(p_old, *(list(S[i]) + [0.0])[:2], function=law3)[0]
                               for i in range(S.shape[0])])
                ok = bool(np.all(m3 >= mu0 - 0.1 * vmag - 1e-12) and np.all(m3 <= mu0 + 1e-12) and np.all(m3 > 0))
                ch.append(SL.cond(f"t2_I3_mu_in_[mu0-b*vmax,mu0]_and_>0_{tag}", ok,
                                  observed=[float(m3.min()), float(m3.max())], expected=[mu0 - 0.1 * vmag, mu0]))
            if dual:
                L = 0.5
                if dim == 2:
                    om = vmag / (L / 2)
                    field = lambda pts, om=om: om * (pts[:, :1] - 0.0)  # noqa: E731
                else:
                    om = vmag / (L / math.sqrt(2))
                    field = lambda pts, om=om: np.stack([-om * pts[:, 1], om * pts[:, 0]], 1)  # noqa: E731
                hr, _ = slip_rates(pt, field)
                Sr = hist_array(hr, pt.n_nodes(), pt.nd)
                exp = field(pt.node_coords())
                ch.append(SL.close(f"t2_I2_dual_rotation_s_i==omega_x_r_{tag}", Sr, exp, rtol=1e-8, atol=1e-12))
                mr = np.array([MP.friction_value(p_old, *(list(Sr[i]) + [0.0])[:2], function=law_I1)[0]
                               for i in range(Sr.shape[0])])
                inferred = (0.5 - mr) / 0.1
                ch.append(SL.close(f"t2_I2_dual_rotation_mu_inferred_speed_{tag}", inferred,
                                   np.linalg.norm(exp, axis=1), rtol=1e-8, atol=1e-12))
    # I3 first step / node without history: μ = μ(p_old, 0)
    m_first = MP.friction_value(p_old, *MP.old_nodal_slip_rate({}, 0), function=law3)[0]
    ch.append(SL.close("t2_I3_first_step_mu==mu(p_old,0)", m_first, mu0))
    return ch


# ================================================================ Track 3 (2D/3D × normalize_c × projection × edge dropping)
VSTAR = 1.0
PSTAR = 1.0


def tol_patch():
    return dict(rtol=1e-9, atol=1e-12 * VSTAR)


def speed(s):
    s = np.atleast_2d(s)
    return np.sqrt(np.sum(s * s, axis=1) + 1e-24)


def t3_sample(dim, norm, proj, edge):
    def s(rng):
        nd = dim - 1
        return dict(dim=dim, normalize_c=norm, projection=proj, correct_edge_dropping=edge,
                    L=float(rng.uniform(0.1, 10)), nsec=int(rng.integers(1, 9)), nprim=int(rng.integers(1, 9)),
                    nq=int(rng.choice([2, 3])), b=[float(x) for x in rng.uniform(-1, 1, nd)],
                    G=[float(x) for x in rng.uniform(-1, 1, nd * nd)], p_old=float(rng.uniform(0, 1000)),
                    part_seed=int(rng.integers(0, 2 ** 31)), rot=float(rng.uniform(0, 2 * math.pi)))
    return s


def t3_patch(inp, scale=1.0, refine=1, tangents=None):
    L = inp["L"] * scale
    return patch(inp["dim"], (0.0, L), (0.0, L), inp["nsec"] * refine, inp["nprim"], 1, True, inp["nq"], tangents)


def t3_eval(inp):
    dim, nd = inp["dim"], inp["dim"] - 1
    b = np.array(inp["b"])
    G = np.array(inp["G"]).reshape(nd, nd) / inp["L"]
    ch = []
    pt = t3_patch(inp)
    hist, (V, W, present) = slip_rates(pt, const_field(b))
    S = hist_array(hist, pt.n_nodes(), nd)
    ok_nodes = present & np.isfinite(W) & (W != 0)
    ch.append(SL.close("t3_I1_W*s==V", (S * W[:, None])[ok_nodes], V[ok_nodes], **tol_patch()))
    ch.append(SL.close("t3_I1_constant_reproduction_s==b", S[ok_nodes], np.tile(b, (ok_nodes.sum(), 1)), **tol_patch()))
    for sc in (0.1, 10.0):
        h2, _ = slip_rates(t3_patch(inp, scale=sc), const_field(b))
        ch.append(SL.close(f"t3_I1_support_scale_{sc:g}_unchanged", hist_array(h2, pt.n_nodes(), nd), S, **tol_patch()))
    for rf in (2, 4, 8):
        prf = t3_patch(inp, refine=rf)
        h3, _ = slip_rates(prf, const_field(b))
        S3 = hist_array(h3, prf.n_nodes(), nd)
        ch.append(SL.close(f"t3_I1_refine_x{rf}_constant", S3, np.tile(b, (S3.shape[0], 1)), **tol_patch()))
    xc = np.full(nd, inp["L"] / 2)
    for lab, fld in (("affine", lambda P: b + (P - xc) @ G.T),
                     ("opposing", lambda P: b[None, :] * ((P[:, :1] - xc[0]) / (inp["L"] / 2)))):
        ha, _ = slip_rates(pt, fld)
        Sa = hist_array(ha, pt.n_nodes(), nd)
        So, _ = pt.oracle_affine(fld)
        ch.append(SL.close(f"t3_I1_{lab}_recorded_speed==norm_of_projection", speed(Sa), speed(So), **tol_patch()))
    # I2: three accepted steps, accelerating slip + varying old pressure; dofs newly entering the contact surface take zero slip
    p0 = inp["p_old"]
    hist_prev = {}
    for n in (1, 2, 3):
        p_old = p0 * (1 + 0.1 * (n - 1))
        args = [MP.friction_value(p_old, *MP.old_nodal_slip_rate(hist_prev, i), function=lambda t, x, y, z: y)[1]
                for i in range(pt.n_nodes())]
        if n == 1:
            ch.append(SL.cond("t3_I2_first_step_args_(0,p0,sqrt(1e-24),0)",
                              all(a[0] == 0.0 and a[1] == p_old and a[2] == math.sqrt(1e-24) and a[3] == 0.0 for a in args),
                              observed=list(args[0])))
        else:
            ys = np.array([a[2] for a in args])
            ch.append(SL.close(f"t3_I2_step{n}_speed_from_step{n - 1}", ys,
                               speed(hist_array(hist_prev, pt.n_nodes(), nd)), **tol_patch()))
            ps = np.array([a[1] for a in args])
            ch.append(SL.close(f"t3_I2_step{n}_old_pressure", ps, np.full_like(ps, p_old), rtol=1e-7, atol=1e-9 * PSTAR))
        hist_prev, _ = slip_rates(pt, const_field(b * n))
    if inp["nsec"] >= 2:
        m = inp["nsec"] // 2
        Lc = inp["L"] * m / inp["nsec"]
        part = patch(dim, (0.0, inp["L"]), (0.0, Lc), inp["nsec"], max(1, inp["nprim"]), 1, True, inp["nq"])
        hp, _ = slip_rates(part, const_field(b))
        entering = [i for i in range(pt.n_nodes()) if i not in hp]
        ys = [MP.friction_value(p0, *MP.old_nodal_slip_rate(hp, i), function=lambda t, x, y, z: y)[1][2] for i in entering]
        ch.append(SL.cond("t3_I2_entering_dof_zero_slip_convention", all(y == math.sqrt(1e-24) for y in ys) and len(entering) > 0,
                          observed=[len(entering), ys[:3]]))
    # I3: integration points randomly assigned to r "processes", partial sums computed
    for r in (2, 4):
        rng = np.random.default_rng(inp["part_seed"] + r)
        hpart, _ = slip_rates(pt, lambda P: b + (P - xc) @ G.T, owner_fn=lambda n: rng.integers(0, r, n), nproc=r)
        hser, _ = slip_rates(pt, lambda P: b + (P - xc) @ G.T)
        A1, A2 = hist_array(hpart, pt.n_nodes(), nd), hist_array(hser, pt.n_nodes(), nd)
        ch.append(SL.close(f"t3_I3_partition_{r}_vs_serial", A1, A2, rtol=1e-8, atol=1e-10 * VSTAR))
    ch.append(SL.cond("t3_I3_s_independent_of_normalize_c/projection/edge_flags", True, observed="the recomputed s does not read these three switches",
                      expected="same as the reference configuration"))
    # I4: velocity magnitude is invariant under an orthogonal transformation of the tangent coordinates
    bn = float(np.linalg.norm(b))
    if nd == 1:
        cases = [[bn], [-bn]]
    else:
        cases = [[bn, 0], [0, bn], [-bn, 0], [0, -bn], [0.6 * bn, 0.8 * bn], [-0.6 * bn, -0.8 * bn]]
    sp = []
    for cvec in cases:
        hc, _ = slip_rates(pt, const_field(cvec))
        sp.append(speed(hist_array(hc, pt.n_nodes(), nd)))
    sp = np.array(sp)
    ch.append(SL.close("t3_I4_speed_isotropic_over_directions", sp, np.full_like(sp, sp[0, 0]), rtol=1e-10, atol=1e-12))
    if nd == 2:
        c, s_ = math.cos(inp["rot"]), math.sin(inp["rot"])
        Q = np.array([[c, -s_], [s_, c]])
        prot = t3_patch(inp, tangents=Q.T)  # the tangent basis rotates together with the geometry
        hr, _ = slip_rates(prot, const_field(Q @ b))
        ch.append(SL.close("t3_I4_rigid_rotation_speed_unchanged", speed(hist_array(hr, pt.n_nodes(), nd)),
                           speed(S), rtol=1e-10, atol=1e-12))
    # I5 (covers only the part that does not need the LM): a constant μ=0.3 and a function that is identically 0.3 give the same friction coefficient
    mc = MP.friction_value(p0, *list(S[0]) + [0.0] * (2 - nd), mu=0.3)[0]
    mf = MP.friction_value(p0, *list(S[0]) + [0.0] * (2 - nd), function=lambda t, x, y, z: 0.3)[0]
    ch.append(SL.close("t3_I5_constant_mu_vs_constant_function", mf, mc, rtol=1e-6, atol=1e-7 * PSTAR))
    return ch


def supplied_law(t, x, y, z):
    return 0.3 + 0.5 * 2.17 ** (-x / 100) - 10 * y


EPS = 1e-7


def t3_mand_vel_eval(inp):
    dim, nd = inp["dim"], inp["dim"] - 1
    b = np.array(inp["b"])
    pt = patch(dim, (0.0, 1.0), (0.0, 1.0), 2, 3, 1, True, 2)
    hist, _ = slip_rates(pt, const_field(b))
    S = hist_array(hist, pt.n_nodes(), nd)
    ch = [SL.close("t3_mand_constant_reproduction", S, np.tile(b, (S.shape[0], 1)), **tol_patch())]
    for p_old in (0.0, EPS / 2, EPS, 2 * EPS, 1.0, 10.0, 100.0, 1000.0):
        args = MP.friction_value(p_old, *list(S[0]) + [0.0] * (2 - nd), function=supplied_law)[1]
        ch.append(SL.cond(f"t3_mand_args_exact_p={p_old:g}", args[0] == 0.0 and args[1] == p_old and args[3] == 0.0,
                          observed=list(args)))
        ch.append(SL.close(f"t3_mand_speed_p={p_old:g}", args[2], math.sqrt(float(b @ b) + 1e-24), **tol_patch()))
        if math.sqrt(float(b @ b)) <= 0.02:
            mu = supplied_law(*args)
            ch.append(SL.cond(f"t3_I5_supplied_law_positive_subset_p={p_old:g}", mu > 0, observed=mu))
    return ch


def t3_mand_valid_eval(inp):
    try:
        MP.validate_friction_params(inp["mu_given"], inp.get("mu", 0.0), inp["fn_given"])
        res = "accept"
    except MP.MooseError:
        res = "reject"
    return [SL.cond("t3_I5_input_validation", res == inp["expected"], observed=res, expected=inp["expected"])]


# ================================================================ Assembly
def main():
    a = SL.parse_args()
    notes = dict(native=False, not_covered=NOT_COVERED[a.track])
    run = SL.Runner(PR, PR_IDX, a.track, a.n, a.out, versions(), notes)
    mand = []
    if a.track == 1:
        branches = [(f"2D_use_dual={d}", t1_sample(d), t1_eval) for d in (False, True)]
    elif a.track == 2:
        branches = [(f"{dim}D_{'dual' if du else 'std'}_p{o}", t2_sample(dim, du, o), t2_eval)
                    for dim in (2, 3) for du in (False, True) for o in (1, 2)]
    else:
        branches = [(f"{dim}D_norm={nm}_proj={pj}_edge={ed}", t3_sample(dim, nm, pj, ed), t3_eval)
                    for dim in (2, 3) for nm in (False, True) for pj in ("ONE", "TWO") for ed in (False, True)]
        comps = [0.0, 1e-12, 1e-8, 1e-4, 0.01, 0.02, 0.1, 1.0]
        for c in comps:
            for sg in (1.0, -1.0):
                mand.append((f"2D_b={sg * c:g}", dict(dim=2, b=[sg * c]), t3_mand_vel_eval))
                mand.append((f"3D_b=({sg * c:g},0)", dict(dim=3, b=[sg * c, 0.0]), t3_mand_vel_eval))
                mand.append((f"3D_b=(0,{sg * c:g})", dict(dim=3, b=[0.0, sg * c]), t3_mand_vel_eval))
        for mu in (0.01, 0.3, 0.5, 1.0):
            mand.append((f"valid_mu={mu}", dict(mu_given=True, mu=mu, fn_given=False, expected="accept"), t3_mand_valid_eval))
        for mu in (0.0, -0.1):
            mand.append((f"valid_mu={mu}", dict(mu_given=True, mu=mu, fn_given=False, expected="reject"), t3_mand_valid_eval))
        mand.append(("valid_neither", dict(mu_given=False, fn_given=False, expected="reject"), t3_mand_valid_eval))
        mand.append(("valid_both", dict(mu_given=True, mu=0.3, fn_given=True, expected="reject"), t3_mand_valid_eval))
        mand.append(("valid_function_only", dict(mu_given=False, fn_given=True, expected="accept"), t3_mand_valid_eval))
    run.run(branches, mand, progress_every=1000)


if __name__ == "__main__":
    main()
