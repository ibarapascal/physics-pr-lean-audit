# zz_bitwise_probe.py: replay the first track's C_zigzag branch (full run, same seed and same indices); output only the points where "ZNyn and ZNd are bitwise identical on the high-voltage side" was rejected, together with the order of magnitude of their residuals
import sys, numpy as np
sys.path.insert(0, "dq1/pp3016")
import spotcheck as S, harness as H
branches, _ = S.track1()
bi = [b[0] for b in branches].index("C_zigzag"); _, samp, ev = branches[bi]
for i in range(2500):
    p = samp(np.random.default_rng([20261098, bi, i]))
    P = H._evaluate(ev, p)
    for c in P.checks:
        if c["name"].startswith("C.ZNyn_vs_ZNd_hv_bitwise") and not c["ok"]:
            a, b = complex(*c["obs"]), complex(*c["exp"])
            print(i, c["name"], "abs", abs(a - b), "rel", abs(a - b) / abs(b), flush=True)
print("DONE")
