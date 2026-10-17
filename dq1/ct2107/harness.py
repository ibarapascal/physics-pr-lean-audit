# harness.py: generic framework for the unified ΔQ1 spot check (an identical copy lives in each PR directory).
# Responsibilities: split the point budget evenly across branches, draw each point's random sample from [seed, branch index, point index],
# run the checks, replay rejected points unchanged to confirm, tally (points / rejects / fails / non-reproducible), keep the first 20 rejected points, write a .tmp file and rename it to produce the output JSON.
# Verdict: any check outside tolerance within a point = reject (must reproduce on replay); no out-of-tolerance check but an exception/NaN/solver failure = fail; otherwise pass.
import json
import math
import os
import sys
import time
import traceback

import numpy as np

DEFAULT_RTOL = 1e-6
DEFAULT_ATOL = 1e-12


def U(rng, a, b):
    return float(rng.uniform(a, b))


def LU(rng, a, b):
    return float(math.exp(rng.uniform(math.log(a), math.log(b))))


def choice(rng, seq):
    return seq[int(rng.integers(len(seq)))]


def _js(x):
    """Convert to a JSON-serializable object (complex -> [re, im], numpy scalar -> python, non-finite number -> string)."""
    if isinstance(x, dict):
        return {str(k): _js(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_js(v) for v in x]
    if isinstance(x, (complex, np.complexfloating)):
        return [_js(float(x.real)), _js(float(x.imag))]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, (float, np.floating)):
        x = float(x)
        return x if math.isfinite(x) else repr(x)
    if isinstance(x, np.ndarray):
        return _js(x.tolist())
    return x


class Point:
    """All checks at one sampled point."""

    def __init__(self):
        self.checks = []    # dict(name, ok, obs, exp, resid, tol)
        self.fails = []     # (name, msg): solver failure / exception / NaN, not counted as a violation
        self.records = []   # quantities that are only recorded, not judged (items the plan marks "record" or "report separately")

    def close(self, name, obs, exp, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL, mode="sum"):
        """Real-valued comparison. mode='sum': |d|<=atol+rtol|exp|; mode='max': |d|<=max(atol, rtol|exp|) (used when the plan says "absolute or relative")."""
        try:
            obs = float(obs)
            exp = float(exp)
        except Exception as e:  # noqa
            self.fails.append((name, "non-numeric: %r" % (e,)))
            return
        if not (math.isfinite(obs) and math.isfinite(exp)):
            if obs == exp:  # both +inf or both -inf
                self.checks.append(dict(name=name, ok=True, obs=obs, exp=exp, resid=0.0, tol=0.0))
            else:
                self.fails.append((name, "nonfinite obs=%r exp=%r" % (obs, exp)))
            return
        d = abs(obs - exp)
        tol = atol + rtol * abs(exp) if mode == "sum" else max(atol, rtol * abs(exp))
        self.checks.append(dict(name=name, ok=bool(d <= tol), obs=obs, exp=exp, resid=d, tol=tol))

    def cclose(self, name, obs, exp, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL, mode="sum", comp=True):
        """Complex-valued comparison. comp=True judges component-wise (real and imaginary parts separately); comp=False uses the complex modulus."""
        obs = complex(obs)
        exp = complex(exp)
        if comp:
            self.close(name + ".re", obs.real, exp.real, rtol, atol, mode)
            self.close(name + ".im", obs.imag, exp.imag, rtol, atol, mode)
            return
        if not all(math.isfinite(v) for v in (obs.real, obs.imag, exp.real, exp.imag)):
            if obs == exp:
                self.checks.append(dict(name=name, ok=True, obs=obs, exp=exp, resid=0.0, tol=0.0))
            else:
                self.fails.append((name, "nonfinite obs=%r exp=%r" % (obs, exp)))
            return
        d = abs(obs - exp)
        tol = atol + rtol * abs(exp) if mode == "sum" else max(atol, rtol * abs(exp))
        self.checks.append(dict(name=name, ok=bool(d <= tol), obs=obs, exp=exp, resid=d, tol=tol))

    def cond(self, name, ok, obs=None, exp=None, resid=None):
        self.checks.append(dict(name=name, ok=bool(ok), obs=obs, exp=exp, resid=resid, tol=None))

    def fail(self, name, msg):
        self.fails.append((name, str(msg)[:300]))

    def record(self, name, **kw):
        self.records.append(dict(name=name, **kw))

    def expect_raise(self, name, fn, exc=(Exception,)):
        """Check the plan requires to raise an exception: raised = pass, not raised = outside the criterion."""
        try:
            fn()
        except exc as e:
            self.checks.append(dict(name=name, ok=True, obs="raised %s" % type(e).__name__, exp="raise", resid=None, tol=None))
            return True
        self.checks.append(dict(name=name, ok=False, obs="no exception", exp="raise", resid=None, tol=None))
        return False

    def bad(self):
        return [c for c in self.checks if not c["ok"]]

    def status(self):
        if self.bad():
            return "reject"
        if self.fails:
            return "fail"
        return "pass"


def _evaluate(ev, params):
    P = Point()
    try:
        ev(params, P)
    except Exception as e:  # exception: not counted as a violation, recorded as a fail (still a reject if an out-of-tolerance check was already recorded)
        P.fail("exception", "%s: %s | %s" % (type(e).__name__, e, traceback.format_exc(limit=3).replace("\n", " / ")[-400:]))
    return P


def _sig(P):
    out = []
    for c in P.bad():
        o = c["obs"]
        if isinstance(o, float) and math.isnan(o):
            o = "nan"
        out.append((c["name"], repr(o)))
    return out


def run(pr, pr_idx, track, n, out, branches, must=None, versions=None, info=None):
    """branches: [(name, sampler(rng)->params, evaluator(params, Point))]; must: [(branch name, params)] mandatory test points for discrete values given by the plan (added on top, not counted in N)."""
    seed = 20261017 + pr_idx * 10 + track
    t0 = time.time()
    k = len(branches)
    base = n // k
    extra = n - base * k
    stats = {}
    first = []
    records = []
    nonrep = 0
    must = must or []
    bidx = {b[0]: i for i, b in enumerate(branches)}

    def handle(bname, bi, i, params, ismust):
        nonlocal nonrep
        if "sampler_error" in params:
            def ev(_p, P, _m=params["sampler_error"]):
                P.fail("sampler", _m)
        else:
            ev = branches[bi][2]
        st = stats[bname]
        tp = time.time()
        P = _evaluate(ev, params)
        st["seconds"] += time.time() - tp
        st["n"] += 1
        if ismust:
            st["n_must"] += 1
        s = P.status()
        if P.records and len(records) < 2000:
            records.append(dict(branch=bname, idx=i, must=ismust, records=_js(P.records)))
        if s == "reject":
            P2 = _evaluate(ev, params)  # replay: recompute on the same input; the bad checks and observed values must be identical
            if _sig(P2) == _sig(P):
                st["rejects"] += 1
                bad = P.bad()
                for c in bad:
                    st["reject_checks"][c["name"]] = st["reject_checks"].get(c["name"], 0) + 1
                if len(first) < 20:
                    first.append(_js(dict(branch=bname, idx=i, must=ismust, params=params,
                                          n_bad_checks=len(bad), bad_checks=bad[:6],
                                          fails=P.fails[:3])))
            else:
                nonrep += 1
                st["nonreproducible"] += 1
        elif s == "fail":
            st["fails"] += 1
            for name, msg in P.fails[:2]:
                key = name + ": " + msg[:120]
                if len(st["fail_examples"]) < 5 and key not in st["fail_examples"]:
                    st["fail_examples"].append(key)
        st["n_checks"] += len(P.checks)

    for bi, (bname, _s, _e) in enumerate(branches):
        stats[bname] = dict(n=0, n_must=0, rejects=0, fails=0, nonreproducible=0, n_checks=0,
                            reject_checks={}, fail_examples=[], seconds=0.0)
    for bname, params in must:
        handle(bname, bidx[bname], -1, params, True)
    for bi, (bname, sampler, _e) in enumerate(branches):
        nb = base + (1 if bi < extra else 0)
        for i in range(nb):
            rng = np.random.default_rng([seed, bi, i])
            try:
                params = sampler(rng)
            except Exception as e:  # the sampler itself errored (e.g. flash failure): record the point as a fail
                params = {"sampler_error": "%s: %s" % (type(e).__name__, str(e)[:200])}
            handle(bname, bi, i, params, False)
            if (i + 1) % 100 == 0:
                print("[%s t%d] %s %d/%d rej=%d fail=%d %.0fs" % (pr, track, bname, i + 1, nb, stats[bname]["rejects"],
                                                                   stats[bname]["fails"], time.time() - t0), flush=True)
    elapsed = time.time() - t0
    res = dict(pr=pr, track=track, n=n, seed=seed, n_branches=k,
               branches={b: dict(v, seconds=round(v["seconds"], 3)) for b, v in stats.items()},
               total_points=sum(v["n"] for v in stats.values()),
               total_rejects=sum(v["rejects"] for v in stats.values()),
               total_fails=sum(v["fails"] for v in stats.values()),
               total_nonreproducible=nonrep,
               first_rejects=first, records_sample=records[:200], n_records=len(records),
               versions=versions or {}, info=info or {}, elapsed_seconds=round(elapsed, 2),
               python=sys.version.split()[0])
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    tmp = out + ".tmp"
    with open(tmp, "w") as f:
        json.dump(_js(res), f, ensure_ascii=False, indent=1)
    os.replace(tmp, out)
    print("[%s t%d] done n=%d rejects=%d fails=%d nonrep=%d %.1fs -> %s" % (
        pr, track, res["total_points"], res["total_rejects"], res["total_fails"], nonrep, elapsed, out), flush=True)
    return res
