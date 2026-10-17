#!/usr/bin/env python3
"""judge_b.py: judge for the controlled comparison: (i) library regression and (ii) invariant spot checks.

Description:
    Runs inside an OpenModelica container with this directory (controlled/) mounted at /work. For every task and
    every variant (pre = library before the change, ctrl_* = known wrong/correct reference patches,
    sNN = generated samples):
    hard-link copy of the full pre-change library tree -> apply the unified diff -> deduplicate by patched content ->
    (i) simulate every regression example and record whether it completes (msl_1335 additionally compares all
        variables with pre, relative difference <= 1e-2) ->
    (ii) compile the test harness once and run N points with -override, judge each point against the task
        invariant and record the worst point.
    Results are cached per deduplication key in work/<task>/cache/<key>.json, so the run can be resumed.
Usage: python3 judge/judge_b.py <out.jsonl> --tasks a,b --variants pre,ctrl,gen --stages i,ii --n 10000 --par 3
Dependencies: Python standard library; omc 1.26.3 and MSL 4.1.0 in the container; libs/ prepared by prep_libs.sh
Change Log: 2026-10-16 first version
"""
import argparse, hashlib, json, math, os, pathlib, random, re, shutil, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

ROOT = pathlib.Path(__file__).resolve().parent.parent
LIBS, WORK, TASKS, GEN, J = ROOT / "libs", ROOT / "work", ROOT / "tasks", ROOT / "samples", ROOT / "judge"
SEED = 20261016
REG_ONLY = set()
# Environment-based exclusion: regression examples whose pre-change version fails to compile on OM 1.26.3 or times out, unrelated to the task, are skipped for all variants.
# Decided solely from the pre logs (environment check), fixed before any sample was examined; run-time failures caused by the original defect itself (e.g. msl_1530 SimpleAir) are not excluded
ENV_EXCLUDE = {
    "msl_1335": {
        "ModelicaTest.Fluid.Dissipation.Verifications.HeatTransfer.HeatExchanger.kc_flatTube",  # record contains a protected section, rejected by OM
        "ModelicaTest.Fluid.Dissipation.Verifications.HeatTransfer.HeatExchanger.kc_flatTube_KC",
        "ModelicaTest.Fluid.Dissipation.Verifications.HeatTransfer.StraightPipe.kc_overall",  # positional arguments of the record type do not match, rejected by OM
        "ModelicaTest.Fluid.Dissipation.Verifications.HeatTransfer.StraightPipe.kc_overall_KC",
        "ModelicaTest.Fluid.Dissipation.Verifications.PressureLoss.StraightPipe.dp_overall_DPMFLOW",
        "ModelicaTest.Fluid.Dissipation.Verifications.PressureLoss.StraightPipe.dp_twoPhaseOverall_DP",
        "ModelicaTest.Fluid.Dissipation.TestCases.PressureLoss.Bend",  # 1800 s timeout
    },
}

MSL323 = ['loadFile("{L}/msl_v323/ModelicaServices/package.mo", uses=false)', 'loadFile("{L}/msl_v323/Complex.mo", uses=false)',
          'loadFile("{L}/msl_v323/Modelica/package.mo", uses=false)']  # uses=false: prevents OM from automatically co-loading the MSL 4.1.0 bundled with the container according to the uses annotation
MSL41 = ['loadModel(Modelica, {{"4.1.0"}})']
CFG = {
    "buildings_2725": dict(lib="mb_309c814f", deps=MSL323, top=["Buildings"], harness="BoilerHarness",
                           ctrl={"wrong": "buildings_2725_wrong.diff", "correct": None}),
    "buildings_4775": dict(lib="mb_b83b2278", deps=MSL41, top=["Buildings"], harness="BoilerHarness",
                           ctrl={"wrong": "buildings_4775_wrong.diff"}),
    "buildings_3027": dict(lib="mb_dff251c2", deps=MSL41, top=["Buildings"], harness="CoilHarness",
                           ctrl={"wrong": "buildings_3027_wrong.diff", "correct": "buildings_3027_correct.diff"}),
    "msl_1335": dict(lib="msl_a3041258", deps=[], top=["ModelicaServices", "Complex.mo", "Modelica", "ModelicaTest"],
                     harness="StepHarness", ctrl={"wrong": "msl_1335_wrong.diff", "correct": "msl_1335_correct.diff"}),
    "msl_1530": dict(lib="msl_7c8ccafa", deps=[], top=["ModelicaServices", "Complex.mo", "Modelica", "ModelicaTest"],
                     harness="GasHarness", ctrl={"wrong": "msl_1530_wrong.diff", "correct": "msl_1530_correct.diff"}),
}


# ---------- unified diff ----------
def split_files(diff):
    files, cur, src = {}, None, None
    for line in diff.splitlines():
        if line.startswith("--- "):
            src = line[4:].strip().split("\t")[0]
            continue
        if line.startswith("+++ "):
            p = line[4:].strip().split("\t")[0]
            if p == "/dev/null":
                p, kind = src, "delete"
            else:
                kind = "new" if src == "/dev/null" else "modify"
            p = p[2:] if p.startswith(("a/", "b/")) else p
            cur = {"kind": kind, "lines": []}
            files[p] = cur
            continue
        if cur is not None:
            cur["lines"].append(line)
    return files


def hunks(lines):
    out, cur = [], None
    for l in lines:
        if l.startswith("@@"):
            cur = {"old": [], "new": []}
            out.append(cur)
        elif cur is not None:
            if l.startswith("-"):
                cur["old"].append(l[1:])
            elif l.startswith("+"):
                cur["new"].append(l[1:])
            elif l.startswith(" ") or l == "":
                cur["old"].append(l[1:] if l else "")
                cur["new"].append(l[1:] if l else "")
    return out


def norm(s):
    return re.sub(r"\s+", "", s)


def apply_hunks(text, hs):
    lines, method = text.split("\n"), "exact"
    for h in hs:
        old, new = h["old"], h["new"]
        while old and old[-1] == "" and new and new[-1] == "":
            old, new = old[:-1], new[:-1]
        k = len(old)
        idx = [i for i in range(len(lines) - k + 1) if [x.rstrip() for x in lines[i:i + k]] == [x.rstrip() for x in old]]
        if len(idx) != 1:
            idx = [i for i in range(len(lines) - k + 1) if [norm(x) for x in lines[i:i + k]] == [norm(x) for x in old]]
            method = "ws_normalized"
        if len(idx) != 1:
            raise ValueError(f"hunk location failed ({len(idx)} candidate positions)")
        lines[idx[0]:idx[0] + k] = new
    return "\n".join(lines), method


def snapshot_scope(task):
    snap = TASKS / task / "snapshot"
    files = {str(p.relative_to(snap)) for p in snap.rglob("*") if p.is_file()}
    dirs = {str(pathlib.Path(f).parent) for f in files}
    return files, dirs


def build_tree(task, vid, diff, check_scope):
    """Returns (tree root, dedup key, applied method, list of changed files); raises ValueError if the patch cannot be applied or goes out of bounds."""
    c = CFG[task]
    root = WORK / task / "trees" / vid
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)
    for t in c["top"]:
        src = LIBS / c["lib"] / t
        if subprocess.run(["cp", "-al", str(src), str(root / t)]).returncode != 0:  # fall back to a full copy when the mount does not support hard links
            shutil.rmtree(root / t, ignore_errors=True)
            subprocess.run(["cp", "-a", str(src), str(root / t)], check=True)
    extra = J / "extra" / task
    if extra.exists():
        for f in extra.rglob("*"):
            if f.is_file():
                dst = root / f.relative_to(extra)
                dst.parent.mkdir(parents=True, exist_ok=True)
                if dst.exists():
                    dst.unlink()
                shutil.copy(f, dst)
    if not diff:
        return root, "pre", "none", []
    files, dirs = snapshot_scope(task)
    methods, touched, digest = [], [], []
    for path, body in sorted(split_files(diff).items()):
        if check_scope and path not in files and not (body["kind"] == "new" and str(pathlib.Path(path).parent) in dirs):
            raise ValueError(f"modifies a file outside the snapshot: {path}")
        dst = root / path
        if body["kind"] == "delete":
            dst.unlink()
            txt = "<deleted>"
        elif body["kind"] == "new":
            txt = "\n".join(h for hk in hunks(body["lines"]) for h in hk["new"]) + "\n"
        else:
            txt, m = apply_hunks(dst.read_text(errors="replace"), hunks(body["lines"]))
            methods.append(m)
        if body["kind"] != "delete":
            if dst.exists():
                dst.unlink()  # break the hard link before writing, never modify libs/
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.write_text(txt)
        touched.append(path)
        digest.append(path + "\0" + txt)
    key = hashlib.sha256("\n\0".join(digest).encode()).hexdigest()[:16]
    return root, key, ("ws_normalized" if "ws_normalized" in methods else "exact"), touched


def load_lines(task, root):
    c = CFG[task]
    L = str(LIBS)
    out = [d.format(L=L) for d in c["deps"]]
    for t in c["top"]:
        p = root / t
        out.append(f'loadFile("{p / "package.mo" if p.is_dir() else p}", uses=false)')
    return "\n".join(x + "; getErrorString();" for x in out)


def omc(mos_text, d, timeout):
    d.mkdir(parents=True, exist_ok=True)
    (d / "run.mos").write_text(mos_text)
    t0 = time.time()
    try:
        r = subprocess.run(["omc", "run.mos"], cwd=d, capture_output=True, text=True, timeout=timeout)
        out, rc = r.stdout + r.stderr, r.returncode
    except subprocess.TimeoutExpired as e:
        out, rc = (e.stdout or b"").decode(errors="replace") if isinstance(e.stdout, bytes) else (e.stdout or ""), "timeout"
    (d / "omc.log").write_text(out)
    return out, rc, round(time.time() - t0, 1)


# ---------- (i) Regression ----------
def regression_models(task):
    m = json.loads((TASKS / task / "meta.json").read_text())
    names = []
    for e in m["regression_examples"] + m.get("regression_pr_added_tests", []):
        if REG_ONLY and e["model"].split(".")[-1] not in REG_ONLY:
            continue
        if e["model"] in ENV_EXCLUDE.get(task, ()):
            continue
        if e["model"] not in names:
            names.append(e["model"])
    return names


def read_csv(p):
    with open(p) as f:
        head = [h.strip().strip('"') for h in f.readline().split(",")]
        rows = {}
        for line in f:
            v = line.strip().split(",")
            try:
                t = float(v[0])
            except ValueError:
                continue
            rows.setdefault(t, v)
    return head, rows


def compare_csv(ref, new):
    hr, rr = read_csv(ref)
    hn, rn = read_csv(new)
    ix = {h: i for i, h in enumerate(hn)}
    worst = None
    for j, h in enumerate(hr):
        if j == 0 or h not in ix:
            continue
        col = []
        for t, v in rr.items():
            try:
                col.append(float(v[j]))
            except ValueError:
                pass
        if not col:
            continue
        rng = max(col) - min(col)
        for t, v in rr.items():
            if t not in rn:
                continue
            try:
                a, b = float(rn[t][ix[h]]), float(v[j])
            except ValueError:
                continue
            tol = max(1e-2 * abs(b), 1e-6 * rng, 1e-12)
            if not (abs(a - b) <= tol):
                ex = abs(a - b) / tol
                if worst is None or ex > worst[0]:
                    worst = (ex, h, t, a, b)
    missing = [h for h in hr[1:] if h not in ix]
    return worst, missing


def run_regression(task, key, root, par):
    models = regression_models(task)
    csv = task == "msl_1335"

    def one(mdl):
        d = WORK / task / "reg" / key / mdl.replace(".", "_")
        fmt = 'outputFormat="csv"' if csv else 'outputFormat="empty"'
        mos = (load_lines(task, root) + f'\nres := simulate({mdl}, {fmt});\nprint(res.messages);\n'
               f'print("\\nRESFILE=" + res.resultFile + "\\n");\ngetErrorString();\n')
        out, rc, wall = omc(mos, d, 1800)
        ok = "The simulation finished successfully" in out
        rec = {"model": mdl, "ok": ok, "rc": rc, "wall_s": wall}
        if not ok:
            rec["tail"] = out[-600:]
        if csv and ok:
            rec["csv"] = str(d / f"{mdl}_res.csv")
        return rec

    with ThreadPoolExecutor(max_workers=par) as ex:
        recs = list(ex.map(one, models))
    return recs


def regression_verdict(task, recs, pre_recs):
    ok = all(r["ok"] for r in recs)
    detail = {}
    if task == "msl_1335" and ok and pre_recs is not None:
        ref = {r["model"]: r.get("csv") for r in pre_recs}
        bad = []
        for r in recs:
            if ref.get(r["model"]) and r.get("csv"):
                worst, missing = compare_csv(ref[r["model"]], r["csv"])
                if worst:
                    bad.append({"model": r["model"], "excess": worst[0], "var": worst[1], "t": worst[2], "new": worst[3], "ref": worst[4]})
        detail["value_mismatch"] = bad
        ok = not bad
    return ok, detail


# ---------- (ii) Spot check ----------
def psat(T):
    return 611.657 * math.exp(17.2799 - 4102.99 / (T - 35.719))


def points(task, n):
    rng = random.Random(SEED)
    pts = []
    while len(pts) < n:
        if task.startswith("buildings_27") or task == "buildings_4775":
            pts.append(dict(y_in=rng.uniform(0, 1), eta=rng.uniform(0.8678, 0.9935), T_in=rng.uniform(293.15, 344.28),
                            T_amb=rng.uniform(263.15, 293.15), Q_nom=rng.uniform(3000.0, 1689847.79),
                            m_rel=1.0 - rng.random()))  # (0,1]
        elif task == "buildings_3027":
            T_a, X = rng.uniform(278.15, 303.15), rng.uniform(0.0035, 0.02)
            pw = psat(T_a)
            if X > 0.621964713 * pw / (101325 - (1 - 0.621964713) * pw):
                continue  # inlet relative humidity ≤1 (meta constraint)
            s = rng.randrange(3)
            m1n, m2n, UA = [(3.78, 2.646, 4748.0), (0.1, 0.14, None), (5.0, 5.0, 9999.0)][s]
            if UA is None:  # UA formula of PControl.mo: 2*m1*4200*(T_a1-T_b1)/lmtd(278.15,283.15,303.15,288.15)
                d1, d2 = 303.15 - 283.15, 288.15 - 278.15
                UA = 2 * 0.1 * 4200 * 5 / ((d1 - d2) / math.log(d1 / d2))
            pts.append(dict(T_w=rng.uniform(278.15, 293.15), T_a=T_a, X_w=X, m1_rel=rng.uniform(0.01, 1.0),
                            m2_rel=rng.uniform(-1.0, 1.0), m1_nom=m1n, m2_nom=m2n, UA=UA, sizing=s))
        elif task == "msl_1335":
            pairs = [(0.75, 0.25), (0.25, 0.75), (2200, 10000), (10000, 2200), (0.05, 0.0), (0.0, 0.05), (40000, 300000), (300000, 40000)]
            f, nf = pairs[rng.randrange(8)]
            w = abs(f - nf)
            x = rng.uniform(min(f, nf) - 0.5 * w, max(f, nf) + 0.5 * w)
            if x < 0:
                continue
            pts.append(dict(func=f, nofunc=nf, x=x))
        elif task == "msl_1530":
            pts.append(dict(TT=rng.uniform(273.15, 373.15), pp=rng.uniform(9e4, 1.1e5)))
    return pts


OUTVARS = {"BoilerHarness": "QFue,Qport,Hw,E", "CoilHarness": "H1,H2,mW2,mWat,Q1,Q2",
           "StepHarness": "r,d", "GasHarness": "hmu,RT,pd,du,u"}
CTRL_PARAMS = {"BoilerHarness": None}


def build_harness(task, key, root, qfue_fallback=False):
    c = CFG[task]
    h = c["harness"]
    txt = (J / "harness" / f"{h}.mo").read_text()
    if h == "BoilerHarness":
        txt = txt.replace("QFUE_EXPR", "y_in*Q_nom/boi.eta_nominal" if qfue_fallback else "boi.QFue_flow")
    d = WORK / task / "harn" / (key + ("_fb" if qfue_fallback else ""))
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{h}.mo").write_text(txt)
    mos = (load_lines(task, root) + f'\nloadFile("{h}.mo"); getErrorString();\n'
           f'res := buildModel({h}, stopTime=1, numberOfIntervals=2, tolerance=1e-8, outputFormat="empty");\n'
           f'print(res[1]); getErrorString();\n')
    out, rc, wall = omc(mos, d, 1800)
    exe = d / h
    return (exe if exe.exists() else None), d, out, wall


def run_point(exe, d, p, outvars):
    ov = ",".join(f"{k}={repr(float(v))}" for k, v in p.items() if k != "sizing")
    try:
        r = subprocess.run([str(exe), f"-override={ov}", f"-output={outvars}", "-lv=-LOG_SUCCESS"], cwd=d,
                           capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        return None
    m = re.search(r"time=[^\n]*", r.stdout)
    if r.returncode != 0 or not m:
        return {"__fail__": (r.stdout + r.stderr)[-300:]}
    vals = {}
    for kv in m.group(0).split(","):
        k, _, v = kv.partition("=")
        try:
            vals[k.strip()] = float(v)
        except ValueError:
            return None
    return vals


def violation(task, p, v):
    """Returns (whether violated, exceedance factor)."""
    if task in ("buildings_2725", "buildings_4775"):
        Q, E = v["QFue"], v["E"]
        tol = 1e-3 * max(abs(Q), abs(v["Qport"]), abs(v["Hw"]), 1.0)
        ex = max(-E, E - Q, -Q) / tol
        return ex > 1, ex
    if task == "buildings_3027":
        # steady-state energy closure: under the library convention simplify_mWat_flow=true, the net enthalpy inflow of the ports H1+H2 = 0;
        # if the sample instead drains condensate explicitly, H1+H2 is allowed to lie within the interval mW2*h_liq(T_c), with T_c between the inlet temperatures of the two fluid streams
        res = v["H1"] + v["H2"]
        hl = lambda T: 4184.0 * (T - 273.15)
        c = [v["mW2"] * hl(min(p["T_w"], p["T_a"])), v["mW2"] * hl(max(p["T_w"], p["T_a"]))]
        lo, hi = min(c), max(c)
        tol = 1e-3 * max(abs(v["H1"]), abs(v["H2"]), 1.0)
        ex = min(abs(res), max(lo - res, res - hi, 0.0)) / tol
        return ex > 1, ex
    if task == "msl_1530":
        s = max(abs(v["hmu"] + v["u"]), v["RT"])  # |h| = |hmu+u|
        ex = max(abs(v["hmu"] - v["RT"]), abs(v["pd"] - v["RT"])) / (1e-9 * s)
        ex = max(ex, abs(v["du"]) / (1e-9 * max(abs(v["u"]), v["RT"])))
        return ex > 1, ex
    if task == "msl_1335":
        r, f, nf, x = v["r"], p["func"], p["nofunc"], p["x"]
        ex = max(-r, r - 1) / 1e-12
        if (f > nf and x >= f) or (f < nf and x <= f):
            ex = max(ex, abs(r - 1) / 1e-9)
        if (f > nf and x <= nf) or (f < nf and x >= nf):
            ex = max(ex, abs(r) / 1e-9)
        return ex > 1, ex
    raise KeyError(task)


def step_global_checks(pts, vals):
    """msl_1335: for the same (func,nofunc), sorted by x, non-decreasing from nofunc toward func; derivative sign agrees with the finite difference."""
    groups = {}
    for p, v in zip(pts, vals):
        if v is not None and "__fail__" not in v:
            groups.setdefault((p["func"], p["nofunc"]), []).append((p["x"], v["r"], v["d"]))
    bad = []
    for (f, nf), g in groups.items():
        g.sort()
        sgn = 1 if f > nf else -1
        for (x0, r0, d0), (x1, r1, d1) in zip(g, g[1:]):
            dr = r1 - r0
            if sgn * dr < -1e-9:
                bad.append({"kind": "monotone", "func": f, "nofunc": nf, "x0": x0, "x1": x1, "r0": r0, "r1": r1})
            if abs(dr) > 1e-9 and (d0 * dr < 0 and d1 * dr < 0):
                bad.append({"kind": "der_sign", "func": f, "nofunc": nf, "x0": x0, "x1": x1, "dr": dr, "d0": d0, "d1": d1})
    return bad


def run_sampling(task, key, root, n, par):
    h = CFG[task]["harness"]
    exe, d, out, wall = build_harness(task, key, root)
    fb = False
    if exe is None and h == "BoilerHarness":
        exe, d, out, wall = build_harness(task, key, root, qfue_fallback=True)
        fb = True
    if exe is None:
        return {"harness_ok": False, "build_tail": out[-800:]}
    pts = points(task, n)
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=par) as ex:
        vals = list(ex.map(lambda p: run_point(exe, d, p, OUTVARS[h]), pts))
    nviol, worst, fails, viols, fail_ex = 0, None, 0, [], []
    for p, v in zip(pts, vals):
        if v is None or "__fail__" in v:
            fails += 1
            if len(fail_ex) < 3:
                fail_ex.append({"p": p, "err": (v or {}).get("__fail__", "timeout")})
            continue
        bad, exn = violation(task, p, v)
        if bad:
            nviol += 1
            if len(viols) < 5:
                viols.append({"p": p, "v": v, "excess": exn})
            if worst is None or exn > worst["excess"]:
                worst = {"p": p, "v": v, "excess": exn}
    rec = {"harness_ok": True, "qfue_fallback": fb, "build_wall_s": wall, "run_wall_s": round(time.time() - t0, 1),
           "n": n, "n_failed": fails, "failed_examples": fail_ex, "n_violation": nviol, "worst": worst, "first_violations": viols}
    if task == "msl_1335":
        g = step_global_checks(pts, vals)
        rec["global_violations"] = len(g)
        rec["global_first"] = g[:5]
        rec["n_violation"] += len(g)
    if task == "buildings_3027":  # save point-wise data to cross-check the energy-closure convention on pre
        with open(WORK / task / f"pts_{key}.jsonl", "w") as f:
            for p, v in zip(pts, vals):
                v = None if (v is None or "__fail__" in v) else v
                f.write(json.dumps({"p": p, "v": v}) + "\n")
    rec["rejected"] = rec["n_violation"] > 0
    return rec


# ---------- Variants ----------
def variants(task, which):
    out = []
    if "pre" in which:
        out.append(("pre", "", False))
    if "ctrl" in which:
        for k, f in CFG[task]["ctrl"].items():
            out.append((f"ctrl_{k}", (J / "controls" / f).read_text() if f else "", False))
    if "gen" in which:
        for p in sorted((GEN / task / "codex").glob("s*.json")):
            rec = json.loads(p.read_text())
            out.append((p.stem, rec.get("diff") or None, True))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--tasks", default=",".join(CFG))
    ap.add_argument("--variants", default="pre,ctrl,gen")
    ap.add_argument("--stages", default="i,ii")
    ap.add_argument("--n", type=int, default=10000)
    ap.add_argument("--par", type=int, default=3)
    ap.add_argument("--ppar", type=int, default=6)
    ap.add_argument("--reg-only", default="", help="run only these regression models (comma-separated, for debugging)")
    a = ap.parse_args()
    stages = a.stages.split(",")
    REG_ONLY.update(x for x in a.reg_only.split(",") if x)
    outf = open(a.out, "a")
    for task in a.tasks.split(","):
        cache = WORK / task / "cache"
        cache.mkdir(parents=True, exist_ok=True)
        pre_reg = None
        for vid, diff, is_gen in variants(task, a.variants.split(",")):
            rec = {"task": task, "variant": vid, "n_points": a.n}
            if is_gen and diff is None:
                rec.update(valid=False, invalid_reason="no diff")
                outf.write(json.dumps(rec, ensure_ascii=False) + "\n"); outf.flush(); print(task, vid, "invalid: no diff", flush=True)
                continue
            try:
                root, key, method, touched = build_tree(task, vid, diff, check_scope=is_gen)
            except (ValueError, FileNotFoundError) as e:
                rec.update(valid=False, invalid_reason=str(e))
                outf.write(json.dumps(rec, ensure_ascii=False) + "\n"); outf.flush(); print(task, vid, "invalid:", e, flush=True)
                continue
            rec.update(valid=True, key=key, apply_method=method, touched=touched)
            cf = cache / f"{key}_n{a.n}.json"
            res = json.loads(cf.read_text()) if cf.exists() else {}
            if "i" in stages and "reg" not in res:
                res["reg"] = run_regression(task, key, root, a.par)
            if "i" in stages:
                if vid == "pre":
                    pre_reg = res["reg"]
                ok, detail = regression_verdict(task, res["reg"], pre_reg)
                res["reg_pass"], res["reg_detail"] = ok, detail
            if "ii" in stages and "samp" not in res:
                res["samp"] = run_sampling(task, key, root, a.n, a.ppar)
            cf.write_text(json.dumps(res, ensure_ascii=False))
            rec.update(res)
            rec.pop("reg", None)
            rec["reg_failed"] = [r["model"] for r in res.get("reg", []) if not r["ok"]]
            outf.write(json.dumps(rec, ensure_ascii=False) + "\n"); outf.flush()
            s = res.get("samp", {})
            print(task, vid, key, "reg_pass", res.get("reg_pass"), "fail", rec["reg_failed"][:4],
                  "| samp rejected", s.get("rejected"), "viol", s.get("n_violation"), "failed_pts", s.get("n_failed"), flush=True)
            shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()
