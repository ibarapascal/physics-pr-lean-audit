"""Shared skeleton of the dq1 spot checks (ΔQ1 unified spot-check protocol §3–§4): seeds, even split across branches, uniform-sampling driver, tolerance verdict,
replay of rejected points, JSON output (write PATH.tmp first, then rename). Each PR directory holds an identical copy; scripts import only the copy in their own directory.

Verdict convention:
- Each check at a point gives observed / expected / residual / tol / ok; if any check is finite and outside tolerance, the point is a reject candidate;
- A reject candidate is recomputed immediately on the same input and counts as a reject only if the result is bitwise identical (protocol §4 "can be replayed on its own"); otherwise it is counted in rejected_not_replayable;
- An exception, or a check whose observed is NaN/inf with no other true violation, counts as a fail (not a violation, listed separately).
"""
import argparse
import json
import math
import os
import platform
import sys
import time
import traceback

import numpy as np

DEFAULT_RTOL = 1e-6
DEFAULT_ATOL = 1e-12


def seed_for(pr_idx, track):
    return 20261017 + pr_idx * 10 + track


def split_counts(n, k):
    base, rem = divmod(n, k)
    return [base + (1 if i < rem else 0) for i in range(k)]


def _finite(x):
    try:
        arr = np.asarray(x, dtype=float)
    except (TypeError, ValueError):
        return True
    return bool(np.all(np.isfinite(arr)))


def close(name, observed, expected, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL):
    """|obs-exp| <= atol + rtol*|exp| (element-wise, worst case)."""
    obs = np.asarray(observed, dtype=float)
    exp = np.asarray(expected, dtype=float)
    if not (_finite(obs) and _finite(exp)):
        return dict(name=name, observed=observed, expected=expected, residual=None,
                    tol=dict(rtol=rtol, atol=atol), ok=False, nonfinite=True)
    diff = np.abs(obs - exp)
    lim = atol + rtol * np.abs(exp)
    excess = diff - lim
    worst = int(np.argmax(excess)) if excess.ndim else 0
    ok = bool(np.all(diff <= lim))
    res = float(diff.reshape(-1)[worst]) if diff.ndim else float(diff)
    rel = res / float(np.abs(exp).reshape(-1)[worst]) if float(np.abs(exp).reshape(-1)[worst]) != 0 else None
    return dict(name=name, observed=observed, expected=expected, residual=res, rel_residual=rel,
                tol=dict(rtol=rtol, atol=atol), ok=ok)


def cond(name, ok, observed=None, expected=None, residual=None, tol=None):
    """Boolean check (inequality, sign, whether an error is raised, etc.). Recorded as nonfinite when observed is non-finite."""
    nonfinite = observed is not None and not _finite(observed)
    return dict(name=name, observed=observed, expected=expected, residual=residual, tol=tol,
                ok=bool(ok) and not nonfinite, nonfinite=nonfinite)


def jsonable(x):
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return [jsonable(v) for v in x.tolist()]
    if isinstance(x, (np.floating, float)):
        v = float(x)
        if math.isnan(v):
            return "nan"
        if math.isinf(v):
            return "inf" if v > 0 else "-inf"
        return v
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    return x


def _signature(checks):
    sig = []
    for c in checks:
        sig.append((c["name"], repr(jsonable(c.get("observed"))), bool(c["ok"])))
    return sig


def classify(checks):
    bad = [c for c in checks if not c["ok"] and not c.get("nonfinite")]
    nonfinite = [c for c in checks if c.get("nonfinite")]
    if bad:
        return "rejected", bad
    if nonfinite:
        return "failed", nonfinite
    return "ok", []


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--track", type=int, choices=[1, 2, 3], required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--out", required=True)
    return ap.parse_args()


class Runner:
    """branches: [(name, sample_fn(rng)->inputs, eval_fn(inputs)->checks[, replay_fn(inputs, bad_checks)->checks])];
    mandatory: [(label, inputs, eval_fn)] (discrete values / grids / endpoints from the plan; mandatory test points added on top, not counted in N)."""

    def __init__(self, pr, pr_idx, track, n, out, versions, notes=None):
        self.pr, self.pr_idx, self.track, self.n, self.out = pr, pr_idx, track, n, out
        self.seed = seed_for(pr_idx, track)
        self.rng = np.random.default_rng(self.seed)
        self.versions = versions
        self.notes = notes or {}
        self.first_rejections = []
        self.first_failures = []
        self.check_stats = {}
        self.t0 = time.time()

    def _stat(self, checks, state, bad):
        for c in checks:
            s = self.check_stats.setdefault(c["name"], dict(n=0, rejected=0, nonfinite=0))
            s["n"] += 1
        if state == "rejected":
            for c in bad:
                self.check_stats[c["name"]]["rejected"] += 1
        elif state == "failed":
            for c in bad:
                self.check_stats[c["name"]]["nonfinite"] += 1

    def _eval_point(self, group, inputs, eval_fn, tally, replay_fn=None):
        tally["n"] += 1
        try:
            checks = eval_fn(inputs)
        except Exception as e:  # noqa: BLE001
            tally["failed"] += 1
            if len(self.first_failures) < 20:
                self.first_failures.append(dict(group=group, inputs=inputs, error=f"{type(e).__name__}: {e}",
                                                trace=traceback.format_exc(limit=3)))
            return
        state, bad = classify(checks)
        if state == "rejected":
            try:
                checks2 = replay_fn(inputs, bad) if replay_fn else eval_fn(inputs)
                replay = _signature(checks2) == _signature(checks)
            except Exception:  # noqa: BLE001
                replay = False
            if not replay:
                tally["rejected_not_replayable"] += 1
                state = "nonreplay"
            else:
                tally["rejected"] += 1
                if len(self.first_rejections) < 20:
                    self.first_rejections.append(dict(group=group, inputs=inputs, failed_checks=bad,
                                                      replay_identical=True))
        elif state == "failed":
            tally["failed"] += 1
            if len(self.first_failures) < 20:
                self.first_failures.append(dict(group=group, inputs=inputs, error="nonfinite observed",
                                                checks=bad))
        self._stat(checks, state, bad)

    def run(self, branches, mandatory, progress_every=500):
        counts = split_counts(self.n, len(branches))
        branch_tally = {}
        done = 0
        for br, nb in zip(branches, counts):
            name, sample_fn, eval_fn = br[:3]
            replay_fn = br[3] if len(br) > 3 else None
            tally = dict(n=0, rejected=0, failed=0, rejected_not_replayable=0)
            branch_tally[name] = tally
            for _ in range(nb):
                inputs = sample_fn(self.rng)
                self._eval_point(name, inputs, eval_fn, tally, replay_fn)
                done += 1
                if progress_every and done % progress_every == 0:
                    print(f"[{self.pr} t{self.track}] {done}/{self.n} sampled points, "
                          f"{time.time() - self.t0:.1f}s", flush=True)
        mand_tally = dict(n=0, rejected=0, failed=0, rejected_not_replayable=0)
        for label, inputs, eval_fn in mandatory:
            self._eval_point("mandatory:" + label, inputs, eval_fn, mand_tally)
        elapsed = time.time() - self.t0
        total_rej = sum(t["rejected"] for t in branch_tally.values()) + mand_tally["rejected"]
        total_fail = sum(t["failed"] for t in branch_tally.values()) + mand_tally["failed"]
        result = dict(
            pr=self.pr, track=self.track, n=self.n, seed=self.seed,
            branches=branch_tally, mandatory=mand_tally,
            total_rejected=total_rej, total_failed=total_fail,
            total_rejected_not_replayable=sum(t["rejected_not_replayable"] for t in branch_tally.values())
            + mand_tally["rejected_not_replayable"],
            check_stats=self.check_stats,
            first_rejections=self.first_rejections,
            first_failures=self.first_failures,
            versions=self.versions,
            host=dict(node=platform.node(), python=sys.version.split()[0], machine=platform.machine()),
            notes=self.notes,
            elapsed_s=elapsed,
            elapsed_per_point_s=elapsed / max(1, self.n + mand_tally["n"]),
        )
        os.makedirs(os.path.dirname(os.path.abspath(self.out)), exist_ok=True)
        tmp = self.out + ".tmp"
        with open(tmp, "w") as f:
            json.dump(jsonable(result), f, ensure_ascii=False, indent=1)
        os.replace(tmp, self.out)
        print(f"[{self.pr} t{self.track}] done: n={self.n} mandatory={mand_tally['n']} "
              f"rejected={total_rej} failed={total_fail} elapsed={elapsed:.1f}s -> {self.out}", flush=True)
        return result
