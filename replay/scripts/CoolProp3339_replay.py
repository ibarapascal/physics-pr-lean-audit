#!/usr/bin/env python3
"""
CoolProp3339_replay.py - native replay for CoolProp PR #3339 against lean/PhysAudit/CoolProp3339.lean

Description
  For the eight rhosr-CS fluids (R1233zd(E), R1234yf, R1234ze(E), R124, R152A, R22, R245fa, R32) samples
  N single-phase states (T in [230, 420] K clipped to each fluid's range; vapor, compressed-liquid and
  supercritical branches) and compares the native viscosity with the Lean definitions
  (xRhosr, psiLiq, horner3, etastarRef, etastarFluid, viscosityRhosr, viscosityDiluteKT), evaluated in
  60-digit Decimal arithmetic from the native state's (rhomolar, R, tau, dalphar_dTau, alphar).
  Theorem checks:
    (a) viscosity_eq_old_of_xc_two : merged formula with x_c = 2 vs pre-PR formula (float64 bitwise and native).
    (b) psiLiq_at_crossover / continuous_etastarFluid : psi(x_c) = 1/2; |eta*(x_c +- h) - eta*(x_c)| -> 0.
    (c) blendExponent_close_outside : |g_a - g_b| <= exp(-k delta) |f_L - f_V| on random out-of-band x;
        etastarFluid_ne_of_xc_lt : in exact arithmetic the two versions differ at every x, while float64
        rounds the out-of-band difference away (reported example).
    dilute_limit_broken_at_xc_zero : predicted eta/eta_dilute at x = 0 for x_c = 0 vs native clone.
  Native clone checks (PR test pattern): clones with x_c in {2, 12, 3, 0} registered through
  add_fluids_as_JSON. If the clone with x_c = 12 returns the shipped value bit-for-bit, the native build
  does not contain PR #3339 and the clone values are compared against the pre-PR formula instead.

Environment
  PR #3339 merge commit 5418cec58c52183e3a010b2894e08b818ed4cc81 (2026-08-29). The PyPI wheel CoolProp==8.0.0
  (tag v8.0.0, 2026-06-27) is an ancestor of that commit and still has the literal `x - 2`
  (TransportRoutines.cpp:1278), so on the `py2` venv only the x_c = 2 identities can be checked against the
  merged formula; the x_c != 2 semantics need a source build of the merge commit
  (pip install "git+https://github.com/CoolProp/CoolProp@5418cec58c52183e3a010b2894e08b818ed4cc81#subdirectory=wrappers/Python").

Usage
  py2/bin/python -I CoolProp3339_replay.py [N=2000] [seed=3339] > CoolProp3339_replay.json

Dependencies
  Python >= 3.9, CoolProp (8.0.0 or a build of the merge commit); standard library only otherwise.

Change Log
  2026-10-17  first version
"""
import json
import math
import random
import sys
from decimal import Decimal, getcontext, localcontext

getcontext().prec = 60

N = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 3339

try:
    import CoolProp
    import CoolProp.CoolProp as CP
except Exception as exc:  # pragma: no cover
    print(json.dumps({"module": "CoolProp3339", "status": "SKIP", "reason": "CoolProp import failed: %s" % exc}))
    sys.exit(0)

FLUIDS = ["R1233zd(E)", "R1234yf", "R1234ze(E)", "R124", "R152A", "R22", "R245fa", "R32"]
K = 100.0  # crossover_sharpness (TransportRoutines.cpp:1280)

# ---------------------------------------------------------------------------
# Lean definitions, Decimal (high precision) and float64 (C++ operation order)
# ---------------------------------------------------------------------------
D = Decimal  # Decimal(float) is the exact binary value


def lean_x(rhomolar, R, tau, dadt, ar, rhosr_crit):  # xRhosr
    return D(rhomolar) * D(R) * (D(tau) * D(dadt) - D(ar)) / D(rhosr_crit)


def lean_psi(k, x, xc):  # psiLiq
    return D(1) / (D(1) + (-D(k) * (x - D(xc))).exp())


def lean_horner(c, x):  # horner3
    c = [D(v) for v in c]
    return c[0] + x * (c[1] + x * (c[2] + x * c[3]))


def lean_g(k, xc, cL, cV, x):  # blendExponent with f_L, f_V = horner3
    p = lean_psi(k, x, xc)
    return p * lean_horner(cL, x) + (D(1) - p) * lean_horner(cV, x)


def lean_etastar(C, k, xc, cL, cV, x):  # etastarFluid
    return D(1) + D(C) * (lean_g(k, xc, cL, cV, x).exp() - D(1))


def lean_dilute(T, eps_k, sigma, M):  # viscosityDiluteKT
    Tstar = D(T) / D(eps_k)
    sigma_nm = D(sigma) * D("1e9")
    mk = D(M) * D(1000)
    omega = (D("1.16145") * Tstar ** D("-0.14874") + D("0.52487") * (D("-0.77320") * Tstar).exp()
             + D("2.16178") * (D("-2.43787") * Tstar).exp())
    return D("26.692e-9") * (mk * D(T)).sqrt() / (sigma_nm ** 2 * omega)


def f64_etastar(C, xc, cL, cV, x, old=False):
    psi = 1 / (1 + math.exp(-100.0 * (x - 2))) if old else 1 / (1 + math.exp(-K * (x - xc)))
    f_liq = cL[0] + x * (cL[1] + x * (cL[2] + x * (cL[3])))
    f_vap = cV[0] + x * (cV[1] + x * (cV[2] + x * (cV[3])))
    return 1 + C * (math.exp(psi * f_liq + (1 - psi) * f_vap) - 1)


def rel(a, b):
    a, b = D(a), D(b)
    return float(abs(a - b) / max(abs(b), D("1e-300")))


# ---------------------------------------------------------------------------
# fluid data and native state helpers
# ---------------------------------------------------------------------------
def fluid_json(fluid):
    doc = json.loads(CP.get_fluid_param_string(fluid, "JSON"))
    visc = doc[0]["TRANSPORT"]["viscosity"]
    block = visc[0] if isinstance(visc, list) else visc
    return doc, block


def native_inputs(AS):
    try:
        return AS.rhomolar(), AS.gas_constant(), AS.tau(), AS.dalphar_dTau(), AS.alphar()
    except AttributeError:
        # x = rho * s_r / rhosr_crit with s_r = R (tau alphar_tau - alphar)
        return AS.rhomolar(), 1.0, 1.0, AS.smolar_residual(), 0.0


def dilute_params(AS, block):
    if "sigma_eta" in block and "epsilon_over_k" in block:
        return block["epsilon_over_k"], block["sigma_eta"], "fluid file"
    rho_red = AS.rhomolar_reducing() / 1000.0
    return AS.T_reducing() / 1.2593, 0.809 / rho_red ** (1.0 / 3.0) / 1e9, "Chung default (default_transport)"


def sample_state(rng, AS):
    Tc = AS.T_critical()
    Tlo = max(230.0, AS.Ttriple() + 1.0)
    Thi = min(420.0, AS.Tmax())
    rhoc = AS.rhomolar_critical()
    for _ in range(200):
        T = rng.uniform(Tlo, Thi)
        try:
            if T < Tc * 0.999:
                AS.update(CP.QT_INPUTS, 0, T)
                rl = AS.rhomolar()
                AS.update(CP.QT_INPUTS, 1, T)
                rv = AS.rhomolar()
                if rng.random() < 0.5:
                    rho = math.exp(rng.uniform(math.log(1.0), math.log(0.995 * rv)))
                else:
                    rho = rng.uniform(1.002 * rl, 1.08 * rl)
            else:
                rho = math.exp(rng.uniform(math.log(1.0), math.log(2.8 * rhoc)))
            AS.update(CP.DmolarT_INPUTS, rho, T)
            return T, rho
        except ValueError:
            continue
    raise RuntimeError("no single-phase state found")


def register_clone(doc, name, cas, xc):
    fluid = json.loads(json.dumps(doc[0]))
    fluid["INFO"]["NAME"] = name
    fluid["INFO"]["CAS"] = cas
    fluid["INFO"]["ALIASES"] = []
    fluid["INFO"]["REFPROP_NAME"] = "N/A"
    for key in ("INCHI_KEY", "INCHI_STRING", "SMILES", "CHEMSPIDER_ID", "2DPNG_URL"):
        fluid["INFO"].pop(key, None)
    visc = fluid["TRANSPORT"]["viscosity"]
    (visc[0] if isinstance(visc, list) else visc)["x_crossover"] = xc
    CP.add_fluids_as_JSON("HEOS", json.dumps([fluid]))


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
rng = random.Random(SEED)
out = {
    "module": "CoolProp3339",
    "pr": "CoolProp/CoolProp#3339",
    "merge_commit": "5418cec58c52183e3a010b2894e08b818ed4cc81",
    "coolprop_version": getattr(CoolProp, "__version__", "?"),
    "coolprop_gitrevision": getattr(CoolProp, "__gitrevision__", "?"),
    "n_points": N,
    "seed": SEED,
    "fluids": {},
}
per = max(1, N // len(FLUIDS))
honours = None

for fi, fluid in enumerate(FLUIDS):
    rec = {}
    try:
        doc, block = fluid_json(fluid)
        C, cL, cV = block["C"], block["c_liq"], block["c_vap"]
        rsc, xc = block["rhosr_critical"], float(block["x_crossover"])
        AS = CP.AbstractState("HEOS", fluid)
        eps_k, sigma, src = dilute_params(AS, block)
        M = AS.molar_mass()
        rec.update({"type": block.get("type"), "x_crossover": block["x_crossover"], "C": C,
                    "c_liq": cL, "c_vap": cV, "rhosr_critical": rsc, "dilute_params": src,
                    "c_vap0_is_zero": cV[0] == 0})
        max_rel_full = 0.0
        max_rel_bracket_ratio = 0.0
        max_rel_vs_old_f64 = 0.0
        bitwise_new_old_f64 = True
        xs = []
        for _ in range(per):
            T, rho = sample_state(rng, AS)
            eta_nat = AS.viscosity()
            rm, R, tau, dadt, ar = native_inputs(AS)
            x = lean_x(rm, R, tau, dadt, ar, rsc)
            xs.append(float(x))
            star = lean_etastar(C, K, xc, cL, cV, x)
            eta_lean = star * lean_dilute(T, eps_k, sigma, M)
            max_rel_full = max(max_rel_full, rel(eta_nat, eta_lean))
            xf = float(x)
            new64 = f64_etastar(C, xc, cL, cV, xf)
            old64 = f64_etastar(C, xc, cL, cV, xf, old=True)
            if xc == 2.0 and new64 != old64:
                bitwise_new_old_f64 = False
            max_rel_vs_old_f64 = max(max_rel_vs_old_f64, rel(new64, old64))
            # dilute-free check: eta(T, rho)/eta(T, rho') equals the bracket ratio at equal T
            AS.update(CP.DmolarT_INPUTS, max(1e-3, rho * 1e-6), T)
            eta_lo = AS.viscosity()
            rm2, R2, tau2, dadt2, ar2 = native_inputs(AS)
            x2 = lean_x(rm2, R2, tau2, dadt2, ar2, rsc)
            star2 = lean_etastar(C, K, xc, cL, cV, x2)
            max_rel_bracket_ratio = max(max_rel_bracket_ratio,
                                        rel(D(eta_nat) / D(eta_lo), star / star2))
        rec["native_vs_lean_full_max_rel"] = max_rel_full
        rec["native_ratio_vs_lean_bracket_ratio_max_rel"] = max_rel_bracket_ratio
        rec["thm_a_float64_bitwise_new_eq_old_at_xc2"] = bitwise_new_old_f64
        rec["x_range"] = [min(xs), max(xs)]

        # (b) value at the crossover and continuity
        psi_c = 1 / (1 + math.exp(-K * (xc - xc)))
        rec["thm_b_psi_at_crossover_float64"] = psi_c
        rec["thm_b_psi_at_crossover_decimal"] = str(lean_psi(K, D(xc), xc))
        s0 = lean_etastar(C, K, xc, cL, cV, D(xc))
        jumps = []
        for e in range(2, 13):
            h = D(10) ** (-e)
            jumps.append(float(max(abs(lean_etastar(C, K, xc, cL, cV, D(xc) + h) - s0),
                                   abs(lean_etastar(C, K, xc, cL, cV, D(xc) - h) - s0)) / h))
        rec["thm_b_difference_quotient_at_crossover_h=1e-2..1e-12"] = jumps  # bounded => continuous

        # (c) out-of-band bound and the exact-vs-float counterexample
        viol, checked, worst, skipped = 0, 0, 0.0, 0
        ctx = localcontext()
        ctx.__enter__()
        getcontext().prec = 250  # resolves exp(-400) against O(1) exponents
        for _ in range(500):
            xr = D(rng.uniform(0.0, 12.0))
            for b in (1.0, 3.0, 6.0, 12.0):
                a = 2.0
                if xr >= D(max(a, b)):
                    delta = xr - D(max(a, b))
                elif xr <= D(min(a, b)):
                    delta = D(min(a, b)) - xr
                else:
                    continue
                if D(K) * delta > 400:
                    skipped += 1
                    continue
                lhs = abs(lean_g(K, a, cL, cV, xr) - lean_g(K, b, cL, cV, xr))
                bound = (-D(K) * delta).exp() * abs(lean_horner(cL, xr) - lean_horner(cV, xr))
                checked += 1
                if lhs > bound * (1 + D("1e-40")):
                    viol += 1
                if bound > 0:
                    worst = max(worst, float(lhs / bound))
        rec["thm_c_bound_skipped_k_delta_gt_400"] = skipped
        rec["thm_c_bound_checked"] = checked
        rec["thm_c_bound_violations"] = viol
        rec["thm_c_max_lhs_over_bound"] = worst
        x0 = 0.5
        ga, gb = lean_g(K, 2.0, cL, cV, D(x0)), lean_g(K, 3.0, cL, cV, D(x0))
        rec["thm_c_counterexample_x=0.5_a=2_b=3"] = {
            "exact_g_difference": str(ga - gb),
            "exact_etastar_differs": lean_etastar(C, K, 2.0, cL, cV, D(x0)) != lean_etastar(C, K, 3.0, cL, cV, D(x0)),
            "float64_etastar_equal": f64_etastar(C, 2.0, cL, cV, x0) == f64_etastar(C, 3.0, cL, cV, x0),
        }
        ctx.__exit__(None, None, None)

        # dilute limit: x = 0
        rec["dilute_etastar_x0_shipped"] = str(lean_etastar(C, K, xc, cL, cV, D(0)))
        rec["dilute_etastar_x0_if_xc0"] = str(lean_etastar(C, K, 0.0, cL, cV, D(0)))

        # native clones (PR test pattern)
        tag = "%d%d" % (SEED % 1000, fi)
        clones = {}
        for j, xcl in enumerate((2.0, 12.0, 3.0, 0.0)):
            name = "XC%s_%s_%d" % (tag, fluid.replace("(", "").replace(")", ""), j)
            cas = "999-%s-%d" % (tag, j)
            try:
                register_clone(doc, name, cas, xcl)
                clones[xcl] = name
            except Exception as exc:
                rec["clone_error_%g" % xcl] = str(exc)
        pts = [(300.0, 10000.0), (300.0, 1e-3)]
        cl_out = {}
        for (T, rho) in pts:
            try:
                AS.update(CP.DmolarT_INPUTS, rho, T)
            except ValueError as exc:
                cl_out["T=%g,rho=%g" % (T, rho)] = {"error": str(exc)}
                continue
            shipped = AS.viscosity()
            rm, R, tau, dadt, ar = native_inputs(AS)
            x = lean_x(rm, R, tau, dadt, ar, rsc)
            row = {"x": float(x), "shipped": shipped}
            for xcl, name in clones.items():
                ASc = CP.AbstractState("HEOS", name)
                ASc.update(CP.DmolarT_INPUTS, rho, T)
                v = ASc.viscosity()
                pred_new = lean_etastar(C, K, xcl, cL, cV, x) / lean_etastar(C, K, xc, cL, cV, x)
                row["clone_xc=%g" % xcl] = {
                    "ratio_to_shipped": v / shipped,
                    "lean_merged_ratio": float(pred_new),
                    "bitwise_equal_to_shipped": v == shipped,
                }
            cl_out["T=%g,rho=%g" % (T, rho)] = row
        rec["clones"] = cl_out
        r12 = cl_out.get("T=300,rho=10000", {}).get("clone_xc=12")
        if r12 is not None and fluid == "R1233zd(E)":
            honours = not r12["bitwise_equal_to_shipped"]
    except Exception as exc:
        rec["error"] = "%s: %s" % (type(exc).__name__, exc)
    out["fluids"][fluid] = rec

out["native_honours_x_crossover"] = honours
out["note"] = ("native_honours_x_crossover == False means the build predates PR #3339; then clone ratios "
               "are expected to be 1 and only the x_c = 2 comparisons test the merged formula.")
print(json.dumps(out, indent=1, default=str))
