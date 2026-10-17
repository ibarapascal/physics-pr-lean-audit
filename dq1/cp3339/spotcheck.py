# spotcheck.py: unified ΔQ1 spot check for CoolProp#3339 (rhosr-CS viscosity reads x_crossover from the fluid file; merge commit 5418cec5, parent commit afcc7909).
# Usage: python spotcheck.py --track {1,2,3} --n N --out PATH [--base-python PATH]
#   The current interpreter should have CoolProp built from the merge commit; --base-python points to a build of the parent commit (or a CoolProp
#   8.0.0 whose viscosity_rhosr is verbatim identical to the parent commit, env/py2), used for bitwise / cross-build regression comparison; if not given, those checks are not run and the JSON records this.
#   The env/py2 currently on the node is CoolProp 8.0.0 (ae81610e, 63 commits before the merge commit, does not contain this PR); there is no merge-commit build -- see impl_log.md.
# Track 1 = track1.md (primary 1: parent/merge bitwise equality, primary 2: clone ratio formula rtol 1e-12, secondary: dilute limit);
# Track 2 = track2.md (I1 bitwise regression, I2 switching locality and interpolation bounds, I3 dilute-limit extrapolation);
# Track 3 = audit/specs/track3_guided/CoolProp_CoolProp_3339/spec.md §5 (I1 blending weight and a derivative, I2 compatibility, I3 boundedness and positivity, I4 dilute limit).
# Internal mode --serve-base: runs as the baseline process, reads one {"fluid","T","rho"} per line from stdin and outputs {"hex"|"err"} (launched by the main process with --base-python).
import argparse
import copy
import json
import math
import os
import subprocess
import sys


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harness import LU, U, run  # noqa: E402

import CoolProp as CP  # noqa: E402
import CoolProp.CoolProp as CPC  # noqa: E402
from CoolProp.CoolProp import AbstractState  # noqa: E402

PR = "CoolProp/CoolProp#3339"
PR_IDX = 7
FLUIDS = ["R1233zd(E)", "R1234yf", "R1234ze(E)", "R124", "R152A", "R22", "R245fa", "R32"]
FLUIDS_DIR = os.environ.get("CP_FLUIDS_DIR") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "env", "coolprop-src", "dev", "fluids"))
K = 100.0  # switching steepness (crossover_sharpness in the merge-commit source)


def serve_base():
    for line in sys.stdin:
        q = json.loads(line)
        try:
            A = AbstractState("HEOS", q["fluid"])
            A.update(CP.DmolarT_INPUTS, q["rho"], q["T"])
            out = {"hex": float(A.viscosity()).hex(), "ps": float(CPC.PropsSI("V", "T", q["T"], "Dmolar", q["rho"], q["fluid"])).hex()}
        except Exception as e:
            out = {"err": str(e)[:200]}
        sys.stdout.write(json.dumps(out) + "\n")
        sys.stdout.flush()


class Base:
    proc = None
    python = None

    @classmethod
    def ask(cls, fluid, T, rho):
        if cls.python is None:
            return None
        if cls.proc is None:
            cls.proc = subprocess.Popen([cls.python, os.path.abspath(__file__), "--serve-base"], stdin=subprocess.PIPE,
                                        stdout=subprocess.PIPE, text=True, bufsize=1)
        cls.proc.stdin.write(json.dumps(dict(fluid=fluid, T=T, rho=rho)) + "\n")
        cls.proc.stdin.flush()
        return json.loads(cls.proc.stdout.readline())


# ---------------- Fluid JSON and clones ----------------
_JS = {}
_REG = set()


def fluid_json(name):
    if name not in _JS:
        d = json.loads(CPC.get_fluid_param_string(name, "JSON"))[0]
        src = "get_fluid_param_string"
        if "TRANSPORT" not in d:  # the 8.0.0 serialization does not contain TRANSPORT; fall back to dev/fluids in the merge-commit source tree
            with open(os.path.join(FLUIDS_DIR, name + ".json")) as f:
                d = json.load(f)
            src = "dev/fluids file"
        _JS[name] = (d, src)
    return _JS[name][0]


def vblock(d):
    v = d["TRANSPORT"]["viscosity"]
    return v[0] if isinstance(v, list) else v


def model(name):
    v = vblock(fluid_json(name))
    return dict(C=float(v["C"]), cL=[float(c) for c in v["c_liq"]], cV=[float(c) for c in v["c_vap"]],
                xc=float(v["x_crossover"]), rc=float(v["rhosr_critical"]))


def clone(name, tag, idx, **mods):
    """Mimics the PR test: change NAME/CAS, clear aliases and identifier metadata, then register; same-name override (OVERWRITE_FLUIDS=True)."""
    cname = "DQ1_%s_%s" % (name.replace("(", "").replace(")", ""), tag)
    key = (cname, json.dumps(mods, sort_keys=True))
    if key in _REG:
        return cname
    f = copy.deepcopy(fluid_json(name))
    f["INFO"]["NAME"] = cname
    f["INFO"]["CAS"] = "999-%02d-%02d" % (FLUIDS.index(name), idx)
    f["INFO"]["ALIASES"] = []
    f["INFO"]["REFPROP_NAME"] = "N/A"
    for k in ("INCHI_KEY", "INCHI_STRING", "SMILES", "CHEMSPIDER_ID", "2DPNG_URL"):
        f["INFO"].pop(k, None)
    v = vblock(f)
    for k, val in mods.items():
        v[k] = val
    CPC.add_fluids_as_JSON("HEOS", json.dumps([f]))
    _REG.difference_update({r for r in _REG if r[0] == cname})
    _REG.add(key)
    return cname


def eta_as(fluid, T, rho):
    A = AbstractState("HEOS", fluid)
    A.update(CP.DmolarT_INPUTS, rho, T)
    return float(A.viscosity()), A


def xval(A, rc):
    return A.rhomolar() * A.gas_constant() * (A.tau() * A.dalphar_dTau() - A.alphar()) / rc


def horner(c, x):
    return c[0] + x * (c[1] + x * (c[2] + x * c[3]))


def wlog(x, a):
    z = -K * (x - a)
    if z > 700:
        return 0.0
    return 1.0 / (1.0 + math.exp(z))


def eta_formula(eta0, m, x, a):
    L, V = horner(m["cL"], x), horner(m["cV"], x)
    w = wlog(x, a)
    return eta0 * (1 + m["C"] * (math.exp(w * L + (1 - w) * V) - 1)), L, V, w


def sat_rho(fluid, T):
    A = AbstractState("HEOS", fluid)
    A.update(CP.QT_INPUTS, 0, T)
    rl = A.rhomolar()
    A.update(CP.QT_INPUTS, 1, T)
    return rl, A.rhomolar()


def crit(fluid):
    A = AbstractState("HEOS", fluid)
    return A.T_critical(), A.rhomolar_critical(), A.Ttriple(), A.Tmax(), A.pmax()


def sample_rho(rng, fluid, T, lo, hi_super):
    Tc = crit(fluid)[0]
    if T < Tc:
        rl, rv = sat_rho(fluid, T)
        hi = rl
    else:
        rl = rv = None
        hi = hi_super
    for _ in range(200):
        r = U(rng, lo, hi)
        if rl is None or r <= rv or r >= rl:
            return r
    return None


def check_env_honours():
    """Probe at startup: whether the current CoolProp really reads x_crossover (PR test state R1233zd(E), 300 K, 1e4 mol/m³)."""
    try:
        a = clone("R1233zd(E)", "probe2", 90, x_crossover=2.0)
        b = clone("R1233zd(E)", "probe12", 91, x_crossover=12.0)
        e2 = CPC.PropsSI("V", "T", 300.0, "Dmolar", 1e4, a)
        e12 = CPC.PropsSI("V", "T", 300.0, "Dmolar", 1e4, b)
        return dict(honours_x_crossover=bool(abs(e12 / e2 - 1) > 0.01), eta_xc2=e2, eta_xc12=e12,
                    json_source=_JS["R1233zd(E)"][1])
    except Exception as e:
        return dict(probe_error=str(e)[:200])


# ============================ Track 1 ============================
def t1_sampler(fluid):
    def s(rng):
        T = U(rng, 230, 420)
        hi_super = sat_rho(fluid, 230.0)[0]
        return dict(fluid=fluid, T=T, rho=sample_rho(rng, fluid, T, 1.0, hi_super))
    return s


def t1_eval(p, P):
    fluid, T, rho = p["fluid"], p["T"], p["rho"]
    if rho is None:
        P.fail("sampling", "no single-phase density")
        return
    m = model(fluid)
    eta = float(CPC.PropsSI("V", "T", T, "Dmolar", rho, fluid))
    b = Base.ask(fluid, T, rho)
    if b is None:
        P.record("main1.not_run(no --base-python)")
    elif "err" in b:
        P.fail("main1.base", b["err"])
    else:
        P.cond("main1.parent==merge bitwise (float.hex)", b["ps"] == eta.hex(), obs=eta.hex(), exp=b["ps"])
    ca, cb = clone(fluid, "xc2", 2, x_crossover=2.0), clone(fluid, "xc12", 12, x_crossover=12.0)
    ea = float(CPC.PropsSI("V", "T", T, "Dmolar", rho, ca))
    eb = float(CPC.PropsSI("V", "T", T, "Dmolar", rho, cb))
    A = AbstractState("HEOS", fluid)
    A.update(CP.DmolarT_INPUTS, rho, T)
    x_alpha = xval(A, m["rc"])
    x_s = A.rhomolar() * A.smolar_residual() / m["rc"]
    if abs(x_s - x_alpha) > 1e-9 * max(1.0, abs(x_alpha)):
        P.record("main2.x sign/definition mismatch (used x from alphar form)", x_alpha=x_alpha, x_from_smolar_residual=x_s)
    ga = eta_formula(1.0, m, x_alpha, 2.0)[0]
    gb = eta_formula(1.0, m, x_alpha, 12.0)[0]
    P.close("main2.eta_a/eta_b=[1+C(e^ga-1)]/[1+C(e^gb-1)]", ea / eb, ga / gb, 1e-12, 0.0)
    e1 = float(CPC.PropsSI("V", "T", T, "Dmolar", 1e-3, fluid))
    e2 = float(CPC.PropsSI("V", "T", T, "Dmolar", 2e-3, fluid))
    eta0 = 2 * e1 - e2
    P.cond("minor.c_vap[0]==0 (file)", m["cV"][0] == 0.0, obs=m["cV"][0], exp=0.0)
    P.close("minor.eta(rho=1e-3)/eta_dilute(extrap) -> 1", e1 / eta0, 1.0)


def _rho_for_x(fluid, T, target, m):
    """must points: bisect on the single-phase segment to find x(ρ)=target. Returns None if not found."""
    Tc = crit(fluid)[0]
    segs = []
    if T < Tc:
        rl, rv = sat_rho(fluid, T)
        segs = [(1.0, rv), (rl, sat_rho(fluid, 230.0)[0] * 1.0)]
    for lo, hi in segs:
        try:
            def xf(r):
                A = AbstractState("HEOS", fluid)
                A.update(CP.DmolarT_INPUTS, r, T)
                return xval(A, m["rc"]) - target
            flo, fhi = xf(lo), xf(hi)
            if flo * fhi > 0:
                continue
            for _ in range(100):
                mid = 0.5 * (lo + hi)
                fm = xf(mid)
                if flo * fm <= 0:
                    hi, fhi = mid, fm
                else:
                    lo, flo = mid, fm
            return 0.5 * (lo + hi)
        except Exception:
            continue
    return None


def track1():
    br = [(f, t1_sampler(f), t1_eval) for f in FLUIDS]
    must = []
    for f in FLUIDS:
        try:
            m = model(f)
        except Exception as e:
            print("must-point generation skipped for %s: %s" % (f, e), flush=True)
            continue
        for xc in (2.0, 12.0):
            for d in (-0.5, -0.05, -0.01, 0.01, 0.05, 0.5):
                r = _rho_for_x(f, 300.0, xc + d, m)
                if r is not None:
                    must.append((f, dict(fluid=f, T=300.0, rho=r, x_target=xc + d)))
    return br, must


# ============================ Track 2 ============================
def t2_sampler(fluid):
    def s(rng):
        Tc, rc_, Tt, Tmax, pmax = crit(fluid)
        T = U(rng, Tt, 1.5 * Tc)
        return dict(fluid=fluid, T=T, rho=sample_rho(rng, fluid, T, 1e-6 * rc_, sat_rho(fluid, Tt)[0]), rhoc=rc_)
    return s


def t2_eval(p, P):
    fluid, T, rho = p["fluid"], p["T"], p["rho"]
    if rho is None:
        P.fail("sampling", "no single-phase density")
        return
    m = model(fluid)
    eta = float(CPC.PropsSI("V", "T", T, "Dmolar", rho, fluid))
    b = Base.ask(fluid, T, rho)
    if b is None:
        P.record("I1.not_run(no --base-python)")
    elif "err" in b:
        P.fail("I1.base", b["err"])
    else:
        P.cond("I1.parent==merge bitwise", b["ps"] == eta.hex(), obs=eta.hex(), exp=b["ps"])
    A = AbstractState("HEOS", fluid)
    A.update(CP.DmolarT_INPUTS, rho, T)
    x = xval(A, m["rc"])
    L, V = horner(m["cL"], x), horner(m["cV"], x)
    rhoc = p["rhoc"]
    r1, r2 = 1e-6 * rhoc, 2e-6 * rhoc
    eta0 = 2 * float(CPC.PropsSI("V", "T", T, "Dmolar", r1, fluid)) - float(CPC.PropsSI("V", "T", T, "Dmolar", r2, fluid))
    ca = clone(fluid, "xc2", 2, x_crossover=2.0)
    ea = float(CPC.PropsSI("V", "T", T, "Dmolar", rho, ca))
    for bb in (1.0, 3.0, 6.0, 12.0):
        cb = clone(fluid, "xc%g" % bb, int(bb) + 20, x_crossover=bb)
        eb = float(CPC.PropsSI("V", "T", T, "Dmolar", rho, cb))
        lo, hi = min(2.0, bb), max(2.0, bb)
        if x <= lo - 0.2 or x >= hi + 0.2:
            bound = 10 * abs(L - V) * math.exp(-K * 0.2)
            P.cond("I2.locality |eta_a/eta_b-1|<=|fL-fV|e^-20 (order) [b=%g]" % bb, abs(ea / eb - 1) <= bound,
                   obs=abs(ea / eb - 1), exp=bound)
        for tag, e in (("a=2", ea), ("b=%g" % bb, eb)):
            Lg = math.log(1 + (e / eta0 - 1) / m["C"])
            lo_, hi_ = min(L, V), max(L, V)
            tl = 1e-12 + 1e-6 * abs(lo_)
            th = 1e-12 + 1e-6 * abs(hi_)
            P.cond("I2.L in [min(fL,fV),max(fL,fV)] [%s]" % tag, lo_ - tl <= Lg <= hi_ + th, obs=Lg, exp=[lo_, hi_])
    # I3: dilute limit
    rs = [rhoc * f for f in (1e-3, 1e-4, 1e-5, 1e-6)]
    es = [float(CPC.PropsSI("V", "T", T, "Dmolar", r, fluid)) for r in rs]
    slope = (es[2] - es[3]) / (rs[2] - rs[3])
    e_lin0 = es[3] - slope * rs[3]
    P.close("I3.lim_{rho->0} eta (linear extrap) vs eta(1e-6 rho_c)", es[3], e_lin0, 1e-6, 0.0)
    P.record("I3.slopes (eta-eta0)/rho", slopes=[(e - e_lin0) / r for e, r in zip(es, rs)])
    P.cond("I3.c_V[0]==0 (file)", m["cV"][0] == 0.0, obs=m["cV"][0], exp=0.0)


def track2():
    return [(f, t2_sampler(f), t2_eval) for f in FLUIDS], []


# ============================ Track 3 ============================
D7 = (-0.2, -0.05, -0.01, 0.0, 0.01, 0.05, 0.2)


def t3_sampler(fluid):
    def s(rng):
        Tc, rc_, Tt, Tmax, pmax = crit(fluid)
        for _ in range(200):
            T = U(rng, 250.0, min(450.0, Tmax))
            P_ = LU(rng, 1e3, min(20e6, pmax))
            try:
                A = AbstractState("HEOS", fluid)
                A.update(CP.PT_INPUTS, P_, T)
                rho = A.rhomolar()
            except Exception:
                continue
            if rho <= 3 * rc_:
                return dict(fluid=fluid, T=T, rho=rho, P=P_, kind="state")
        return dict(fluid=fluid, T=None, rho=None, kind="state")
    return s


def _clone_eta(fluid, tag, idx, T, rho, **mods):
    return eta_as(clone(fluid, tag, idx, **mods), T, rho)[0]


def t3_state(p, P):
    fluid, T, rho = p["fluid"], p["T"], p["rho"]
    m = model(fluid)
    eta, A = eta_as(fluid, T, rho)
    x = xval(A, m["rc"])
    if p.get("propssi_xcheck"):
        P.close("xcheck PropsSI==AbstractState", float(CPC.PropsSI("V", "T", T, "Dmolar", rho, fluid)), eta, 1e-14, 0.0)
    eta0 = _clone_eta(fluid, "C0", 30, T, rho, C=0.0)
    L, V = horner(m["cL"], x), horner(m["cV"], x)
    C = m["C"]
    # I1: blending-weight form of a=x+d, fixed a=2/12, tail a=x±8
    gs = {}
    for i, d in enumerate(D7):
        e = _clone_eta(fluid, "d%d" % i, 40 + i, T, rho, x_crossover=x + d)
        ex = eta_formula(eta0, m, x, x + d)[0]
        P.close("I1.blend formula [d=%g]" % d, e, ex, 5e-12, 0.0)
        gs[d] = math.log(1 + (e / eta0 - 1) / C)
    for a, tag, idx in ((2.0, "xc2", 2), (12.0, "xc12", 12)):
        e = _clone_eta(fluid, tag, idx, T, rho, x_crossover=a)
        P.close("I1.blend formula [a=%g]" % a, e, eta_formula(eta0, m, x, a)[0], 5e-12, 0.0)
    for i, d in enumerate((8.0, -8.0)):
        e = _clone_eta(fluid, "tail%d" % i, 60 + i, T, rho, x_crossover=x + d)
        P.cond("I1.tail a=x%+g finite" % d, math.isfinite(e), obs=e)
        if math.isfinite(e):
            P.close("I1.tail a=x%+g formula" % d, e, eta_formula(eta0, m, x, x + d)[0], 5e-12, 0.0)
    tolc = 1e-10 * max(1.0, abs(L), abs(V))
    for d in (0.01, 0.05, 0.2):
        P.close("I1.complement g(x+d)+g(x-d)=L+V [d=%g]" % d, gs[-d] + gs[d], L + V, 0.0, tolc)
    E = math.exp(0.5 * (L + V))
    an = -K * eta0 * C * E * 0.25 * (L - V)
    for j, h in enumerate((1e-5, 5e-6)):
        ep = _clone_eta(fluid, "hp%d" % j, 70 + j, T, rho, x_crossover=x + h)
        em = _clone_eta(fluid, "hm%d" % j, 72 + j, T, rho, x_crossover=x - h)
        fd = (ep - em) / (2 * h)
        degenerate = abs(L - V) < 1e-6 or C == 0.0
        tol = (1e-10 * eta0 * max(1.0, C * E * abs(L - V))) if degenerate else 1e-5 * abs(an)
        P.cond("I1.d eta/da at a=x (h=%g)" % h, abs(fd - an) <= tol, obs=fd, exp=an, resid=abs(fd - an))
    # I2: library fluid vs a=2 JSON clone; cross-build
    e2 = _clone_eta(fluid, "xc2", 2, T, rho, x_crossover=2.0)
    P.close("I2.shipped == a=2 clone (1e-14)", eta, e2, 1e-14, 0.0)
    b = Base.ask(fluid, T, rho)
    if b is None:
        P.record("I2.base_vs_merge not_run(no --base-python)")
    elif "err" in b:
        P.fail("I2.base", b["err"])
    else:
        eb = float.fromhex(b["hex"])
        P.record("I2.bitwise base==merge", equal=(b["hex"] == eta.hex()))
        P.close("I2.base vs merge cross-build (1e-13)", eta, eb, 1e-13, 0.0)
    # I3: finite and positive, endpoint bounds, η≥η0 when x≥0, equal-polynomial and C=0 controls
    P.cond("I3.eta finite positive", math.isfinite(eta) and eta > 0, obs=eta)
    endL, endV = eta0 * (1 + C * (math.exp(L) - 1)), eta0 * (1 + C * (math.exp(V) - 1))
    sl = 5e-12 * max(abs(endL), abs(endV), abs(eta))
    P.cond("I3.eta between branch endpoints", min(endL, endV) - sl <= eta <= max(endL, endV) + sl, obs=eta, exp=[endL, endV])
    if x >= 0:
        P.cond("I3.x>=0 => eta>=eta0", eta >= eta0 - 5e-12 * eta0, obs=eta, exp=eta0)
    else:
        _e1, _L, _V, _w = eta_formula(1.0, m, x, m["xc"])
        Eg = math.exp(_w * _L + (1 - _w) * _V)
        P.record("I3.negative x state (reported separately)", x=x, eta=eta, eta0=eta0, C=C,
                 positivity_condition_holds=bool(C <= 1 or Eg > 1 - 1 / C))
    eq2 = _clone_eta(fluid, "eq2", 80, T, rho, x_crossover=2.0, c_liq=list(m["cV"]))
    eq12 = _clone_eta(fluid, "eq12", 81, T, rho, x_crossover=12.0, c_liq=list(m["cV"]))
    P.close("I3.equal polynomials => crossover inert", eq12, eq2, 5e-12, 0.0)
    c012 = _clone_eta(fluid, "C0_12", 82, T, rho, C=0.0, x_crossover=12.0)
    P.close("I3.C=0 => crossover inert", c012, eta0, 5e-12, 0.0)


def t3_dilute(p, P):
    fluid = p["fluid"]
    m = model(fluid)
    Tmax = crit(fluid)[3]
    rc_ = crit(fluid)[1]
    for T in (250.0, 300.0, min(400.0, Tmax)):
        for a, idx in ((2.0, 2), (12.0, 12), (0.0, 85)):
            vals = []
            for f in (1e-6, 1e-8, 1e-10, 1e-12):
                rho = f * rc_
                tag = "xc%g" % a if a != 0.0 else "xc0"
                e, A = eta_as(clone(fluid, tag, idx, x_crossover=a), T, rho)
                eta0 = _clone_eta(fluid, "C0", 30, T, rho, C=0.0)
                x = xval(A, m["rc"])
                P.close("I4.finite-density formula [T=%g,a=%g,rho/rhoc=%g]" % (T, a, f), e, eta_formula(eta0, m, x, a)[0], 5e-12, 0.0)
                vals.append(e / eta0)
            w0 = wlog(0.0, a)
            lim = 1 + m["C"] * (math.exp(w0 * m["cL"][0] + (1 - w0) * m["cV"][0]) - 1)
            for v in vals[2:]:
                P.close("I4.limit expression [T=%g,a=%g]" % (T, a), v, lim, 1e-8, 0.0)
                if a in (2.0, 12.0):
                    P.close("I4.eta/eta0~1 [T=%g,a=%g]" % (T, a), v, 1.0, 1e-8, 0.0)
                else:
                    P.record("I4.a=0 non-unit limit", T=T, ratio=v, limit=lim)


def t3_eval(p, P):
    if p.get("kind") == "dilute":
        return t3_dilute(p, P)
    if p["T"] is None:
        P.fail("sampling", "no admissible PT state")
        return
    return t3_state(p, P)


def track3():
    br = [(f, t3_sampler(f), t3_eval) for f in FLUIDS]
    must = []
    for f in FLUIDS:
        try:
            Tc, rc_, Tt, Tmax, pmax = crit(f)
        except Exception as e:
            print("must-point generation skipped for %s: %s" % (f, e), flush=True)
            continue
        must.append((f, dict(fluid=f, kind="dilute")))
        for tr in (0.70, 0.85, 0.95):
            T = tr * Tc
            if Tt <= T <= Tmax:
                rl, rv = sat_rho(f, T)
                for q, r in ((0, rl), (1, rv)):
                    must.append((f, dict(fluid=f, T=T, rho=r, kind="state", Q=q, propssi_xcheck=True)))
        for tr in (1.01, 1.05):
            T = tr * Tc
            if T < Tmax:
                for rr in (0.1, 0.5, 1.0, 1.5, 2.0, 3.0):
                    must.append((f, dict(fluid=f, T=T, rho=rr * rc_, kind="state", propssi_xcheck=True)))
    must.append(("R1233zd(E)", dict(fluid="R1233zd(E)", T=300.0, rho=10000.0, kind="state", propssi_xcheck=True)))
    return br, must


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--serve-base":
        return serve_base()
    ap = argparse.ArgumentParser()
    ap.add_argument("--track", type=int, choices=(1, 2, 3), required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--base-python", default=None)
    a = ap.parse_args()
    CPC.set_config_bool(CP.OVERWRITE_FLUIDS, True)
    Base.python = a.base_python
    probe = check_env_honours()
    branches, must = {1: track1, 2: track2, 3: track3}[a.track]()
    versions = dict(CoolProp=CP.__version__, gitrevision=getattr(CP, "__gitrevision__", ""),
                    merge_commit="5418cec58c52183e3a010b2894e08b818ed4cc81",
                    parent_commit="afcc79091533dfbd258cd59c065c0765747fd230", base_python=a.base_python,
                    env_probe=probe, fluids_json_dir_fallback=FLUIDS_DIR)
    run(PR, PR_IDX, a.track, a.n, a.out, branches, must=must, versions=versions)
    if Base.proc is not None:
        Base.proc.stdin.close()
        Base.proc.wait(timeout=30)


if __name__ == "__main__":
    main()
