"""Line-by-line Python recomputation of the "nodal slip rate" of dynamic mortar friction in MOOSE#34020 merge commit f90f508 (not native: there is no MOOSE build on the node).
C++ logic covered:
- in ComputeDynamicWeightedGapLMMechanicalContact::computeQpProperties / computeQpIProperties,
  _qp_factor = JxW_msm·coord and _dof_to_weighted_gap[dof].second += test·_qp_factor (accumulated unconditionally after the merge);
- in ComputeDynamicFrictionalForceLMMechanicalContact::computeQpProperties / computeQpIProperties,
  _qp_tangential_velocity_nodal = −relative_velocity·JxW·coord, weighted_tangential_velocity[dof][α] += test·qp_tv·t_α(node);
  relative_velocity = primary-side velocity − secondary-side velocity (each interpolating nodal velocities with the shape functions of the element on its own side);
- timestepSetup: for dofs owned by this process, s = weighted_tangential_velocity / weighted_gap.second;
- oldNodalSlipRate: no history -> (0,0); computeFrictionValue: μ = F(t=0, x=p_old, y=sqrt(s0²+s1²+1e-24), z=0), or the constant mu;
- friction-parameter validation in the constructor (mu range check mu>0; error if both are given / neither is given).
Geometric restrictions: flat interface, axis-aligned rectangular (3D) or line-segment (2D) elements, tensor-product meshes on both the primary and secondary sides; a mortar segment = an interval/rectangle delimited by the union of the breakpoints of both sides;
the secondary-side test function uses Lagrange (standard) or the libMesh-style dual basis (ψ = D·A⁻¹·φ, with A the element mass matrix and D=diag∫φ).
Not covered: the LM obtained from the nonlinear solve, the normal constraint, the Newmark velocity update, MPI communication itself (replaced by partitioned summation), the edge-dropping logic, elemental LM.
"""
import math

import numpy as np


class MooseError(Exception):
    pass


def lagrange_1d(p, xi):
    """Order-p Lagrange on the reference element [−1,1] (nodes in ascending coordinate order). Returns (p+1, len(xi))."""
    xi = np.asarray(xi, float)
    nodes = np.linspace(-1.0, 1.0, p + 1)
    out = np.ones((p + 1, xi.size))
    for a in range(p + 1):
        for b in range(p + 1):
            if a != b:
                out[a] *= (xi - nodes[b]) / (nodes[a] - nodes[b])
    return out


def dual_coeff_1d(p):
    """libMesh compute_dual_shape_coeffs: A_ij=∫φiφj, D=diag∫φi, dual_coeff=A⁻¹D, ψ_i=Σ_j dual_coeff(j,i) φ_j."""
    xg, wg = np.polynomial.legendre.leggauss(p + 1)
    phi = lagrange_1d(p, xg)
    A = (phi * wg) @ phi.T
    Dm = np.diag(phi @ wg)
    C = np.linalg.solve(A, Dm)
    return C.T  # ψ = C.T @ φ


def tensor_dual(p):
    c1 = dual_coeff_1d(p)
    return c1  # rectangular affine element: tensor-product dual coefficients = Kronecker product of the 1D coefficients (applied direction by direction)


class Line:
    """Tensor-product mesh in one direction: interval [x0, x1], n elements of order p."""

    def __init__(self, x0, x1, n, p):
        self.x0, self.x1, self.n, self.p = float(x0), float(x1), int(n), int(p)
        self.h = (self.x1 - self.x0) / self.n
        self.nodes = self.x0 + np.arange(self.n * self.p + 1) * self.h / self.p

    def breakpoints(self):
        return self.x0 + np.arange(self.n + 1) * self.h

    def locate(self, x):
        e = np.clip(np.floor((x - self.x0) / self.h).astype(int), 0, self.n - 1)
        xi = 2.0 * (x - (self.x0 + e * self.h)) / self.h - 1.0
        return e, xi


def segments_1d(sec, prim, nq):
    """Mortar segments: sub-intervals obtained by partitioning the overlap region with the union of the secondary-side and primary-side breakpoints, each segment with nq Gauss points."""
    lo, hi = max(sec.x0, prim.x0), min(sec.x1, prim.x1)
    if hi <= lo:
        return None
    bp = np.unique(np.concatenate([sec.breakpoints(), prim.breakpoints(), [lo, hi]]))
    bp = bp[(bp >= lo) & (bp <= hi)]
    a, b = bp[:-1], bp[1:]
    keep = b - a > 0
    a, b = a[keep], b[keep]
    xg, wg = np.polynomial.legendre.leggauss(nq)
    mid = 0.5 * (a + b)
    half = 0.5 * (b - a)
    xq = (mid[:, None] + half[:, None] * xg[None, :]).ravel()
    wq = (half[:, None] * wg[None, :]).ravel()
    # the element a point within a segment belongs to is decided by the segment midpoint (avoids ambiguity of endpoint assignment)
    mids = np.repeat(mid, nq)
    es, _ = sec.locate(mids)
    ep, _ = prim.locate(mids)
    _, xis = sec_local(sec, es, xq)
    _, xip = sec_local(prim, ep, xq)
    return dict(xq=xq, wq=wq, es=es, ep=ep, xis=xis, xip=xip)


def sec_local(line, e, x):
    xi = 2.0 * (x - (line.x0 + e * line.h)) / line.h - 1.0
    return e, xi


class Patch:
    """dim=2 (line interface, 1 tangent) or dim=3 (surface interface, 2 tangents). sec/prim are lists of Line per direction."""

    def __init__(self, sec_lines, prim_lines, dual, nq, tangents=None):
        self.sec, self.prim = sec_lines, prim_lines
        self.nd = len(sec_lines)  # interface dimension 1 or 2
        self.p = sec_lines[0].p
        self.dual = dual
        self.nq = nq
        self.C = dual_coeff_1d(self.p) if dual else np.eye(self.p + 1)
        self.tangents = np.eye(self.nd) if tangents is None else np.asarray(tangents, float)
        self.segs = [segments_1d(s, m, nq) for s, m in zip(sec_lines, prim_lines)]

    def node_coords(self):
        if self.nd == 1:
            return self.sec[0].nodes[:, None]
        X, Y = np.meshgrid(self.sec[0].nodes, self.sec[1].nodes, indexing="ij")
        return np.stack([X.ravel(), Y.ravel()], axis=1)

    def n_nodes(self):
        return int(np.prod([len(s.nodes) for s in self.sec]))

    def accumulate(self, slip_sec_fn, slip_prim_fn, owner_fn=None, nproc=1):
        """Returns (V[α], W) for each secondary node and "whether it appears in the map". slip_*_fn(points)->(N, n_t) are the velocity fields in the tangent plane on the two sides.
        owner_fn: assigns each integration point to a "process"; partial sums are computed separately and then added (replaces the summation in communicateVelocities)."""
        if any(s is None for s in self.segs):
            n = self.n_nodes()
            return np.zeros((n, self.nd)), np.zeros(n), np.zeros(n, bool)
        p = self.p
        per_dir = []
        for d in range(self.nd):
            sg = self.segs[d]
            phis = lagrange_1d(p, sg["xis"])  # (p+1, Nq)
            psis = self.C @ phis
            phip = lagrange_1d(p, sg["xip"])
            per_dir.append((sg, phis, psis, phip))
        if self.nd == 1:
            sg, phis, psis, phip = per_dir[0]
            pts = sg["xq"][:, None]
            w = sg["wq"]
            gsec = sg["es"][None, :] * p + np.arange(p + 1)[:, None]  # (p+1, Nq) global node numbers
            gprim = sg["ep"][None, :] * p + np.arange(p + 1)[:, None]
            node_s = self.sec[0].nodes[:, None]
            node_p = self.prim[0].nodes[:, None]
            vs_nodes = slip_sec_fn(node_s)
            vp_nodes = slip_prim_fn(node_p)
            vs = np.einsum("aq,aqk->qk", phis, vs_nodes[gsec])
            vp = np.einsum("aq,aqk->qk", phip, vp_nodes[gprim])
            test = psis
            gtest = gsec
        else:
            (sx, phx, psx, ppx), (sy, phy, psy, ppy) = per_dir
            NX, NY = sx["xq"].size, sy["xq"].size
            X, Y = np.meshgrid(sx["xq"], sy["xq"], indexing="ij")
            pts = np.stack([X.ravel(), Y.ravel()], axis=1)
            w = (sx["wq"][:, None] * sy["wq"][None, :]).ravel()
            nnx_s, nnx_p = len(self.sec[0].nodes), len(self.prim[0].nodes)
            ns_nodes = self.node_coords()
            Xp, Yp = np.meshgrid(self.prim[0].nodes, self.prim[1].nodes, indexing="ij")
            np_nodes = np.stack([Xp.ravel(), Yp.ravel()], axis=1)
            vs_nodes = slip_sec_fn(ns_nodes)
            vp_nodes = slip_prim_fn(np_nodes)
            P1 = p + 1
            la = np.arange(P1)
            gx_s = sx["es"][None, :] * p + la[:, None]
            gy_s = sy["es"][None, :] * p + la[:, None]
            gx_p = sx["ep"][None, :] * p + la[:, None]
            gy_p = sy["ep"][None, :] * p + la[:, None]
            # local node (a,b) → global (gx*nny + gy); shape function = product over directions
            nny_s = len(self.sec[1].nodes)
            nny_p = len(self.prim[1].nodes)
            gsec = (gx_s[:, None, :, None] * nny_s + gy_s[None, :, None, :]).reshape(P1 * P1, NX * NY)
            gprim = (gx_p[:, None, :, None] * nny_p + gy_p[None, :, None, :]).reshape(P1 * P1, NX * NY)
            phis = (phx[:, None, :, None] * phy[None, :, None, :]).reshape(P1 * P1, NX * NY)
            psis = (psx[:, None, :, None] * psy[None, :, None, :]).reshape(P1 * P1, NX * NY)
            phip = (ppx[:, None, :, None] * ppy[None, :, None, :]).reshape(P1 * P1, NX * NY)
            del nnx_s, nnx_p
            vs = np.einsum("aq,aqk->qk", phis, vs_nodes[gsec])
            vp = np.einsum("aq,aqk->qk", phip, vp_nodes[gprim])
            test = psis
            gtest = gsec
        rel = vp - vs  # relative_velocity = primary − secondary
        qp_tv = -rel * w[:, None]  # −relative·JxW·coord
        proj = qp_tv @ self.tangents.T  # components along each tangent (node tangents are constant on a flat interface)
        n = self.n_nodes()
        if owner_fn is None or nproc == 1:
            owners = np.zeros(w.size, int)
            nproc = 1
        else:
            owners = owner_fn(w.size)
        V = np.zeros((nproc, n, self.nd))
        W = np.zeros((nproc, n))
        present = np.zeros(n, bool)
        for a in range(test.shape[0]):
            np.add.at(W, (owners, gtest[a]), test[a] * w)
            for al in range(self.nd):
                np.add.at(V[..., al], (owners, gtest[a]), test[a] * proj[:, al])
            present[gtest[a]] = True
        return V.sum(axis=0), W.sum(axis=0), present

    def oracle_affine(self, slip_fn, nq=6):
        """Independent reference: integrate the analytic slip field directly with a high-order Gauss rule on the secondary elements only (no interpolation from the primary side, no mortar segments), s=V/W."""
        p = self.p
        xg, wg = np.polynomial.legendre.leggauss(nq)
        n = self.n_nodes()
        V = np.zeros((n, self.nd))
        W = np.zeros(n)
        if self.nd == 1:
            L = self.sec[0]
            for e in range(L.n):
                x = L.x0 + e * L.h + 0.5 * L.h * (xg + 1)
                ww = 0.5 * L.h * wg
                psi = self.C @ lagrange_1d(p, xg)
                v = slip_fn(x[:, None]) @ self.tangents.T
                for a in range(p + 1):
                    gi = e * p + a
                    W[gi] += np.sum(psi[a] * ww)
                    V[gi] += (psi[a] * ww) @ v
        else:
            Lx, Ly = self.sec
            psi1 = self.C @ lagrange_1d(p, xg)
            for ex in range(Lx.n):
                for ey in range(Ly.n):
                    x = Lx.x0 + ex * Lx.h + 0.5 * Lx.h * (xg + 1)
                    y = Ly.x0 + ey * Ly.h + 0.5 * Ly.h * (xg + 1)
                    X, Y = np.meshgrid(x, y, indexing="ij")
                    ww = np.outer(0.5 * Lx.h * wg, 0.5 * Ly.h * wg).ravel()
                    v = slip_fn(np.stack([X.ravel(), Y.ravel()], 1)) @ self.tangents.T
                    for a in range(p + 1):
                        for b in range(p + 1):
                            ps = np.outer(psi1[a], psi1[b]).ravel()
                            gi = (ex * p + a) * len(Ly.nodes) + (ey * p + b)
                            W[gi] += np.sum(ps * ww)
                            V[gi] += (ps * ww) @ v
        return V / W[:, None], W


def timestep_setup(V, W, present):
    """timestepSetup after the merge: for dofs owned by this process, s = V/W. Returns {node: (s0, s1)}."""
    hist = {}
    for i in np.nonzero(present)[0]:
        norm = W[i]
        s = V[i] / norm
        hist[int(i)] = (float(s[0]), float(s[1]) if s.size > 1 else 0.0)
    return hist


def old_nodal_slip_rate(hist, node):
    return hist.get(int(node), (0.0, 0.0))


def friction_value(contact_pressure, tv, tv_dir, mu=None, function=None):
    """computeFrictionValue: returns (μ, the argument passed to the function, or None)."""
    if function is None:
        return mu, None
    y = math.sqrt(tv * tv + tv_dir * tv_dir + 1.0e-24)
    args = (0.0, contact_pressure, y, 0.0)
    return function(*args), args


def validate_friction_params(mu_given, mu_value, function_given):
    """Constructor validation: addRangeCheckedParam("mu", "mu > 0"), and the errors for both given / neither given."""
    if mu_given and not (mu_value > 0):
        raise MooseError("Range check failed for parameter mu: mu > 0")
    if not function_given and not mu_given:
        raise MooseError("A coefficient of friction needs to be provided as a constant value of via a function.")
    if function_given and mu_given:
        raise MooseError("Either provide a constant coefficient of friction or a function ... not both")
    return True
