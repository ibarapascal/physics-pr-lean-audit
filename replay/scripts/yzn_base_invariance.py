"""yzn_base_invariance.py: check whether the zero-sequence short-circuit impedance (in ohms) of Yzn (compared with Dyn) varies with the calculation base net.sn_mva,
for both Z_N = 0 and Z_N = 10 ohm; whichever pandapower environment it is run in determines which commit (merge commit or parent commit) is measured.
Usage: <env>/bin/python scripts/yzn_base_invariance.py <label> > outputs/yzn_base_invariance_<label>.json (reuses the network-building function from pandapower3016_replay.py)"""
import json, sys, importlib.metadata as md
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import pandapower3016_replay as R

P = R._pp()
p = {"vn_hv": 20.0, "vn_lv": 0.4, "sn_t": 0.63, "vk": 6.0, "vkr": 1.0, "mag0": 100.0, "mag0_rx": 0.0,
     "si": 0.9, "sn_net": 1.0, "s_sc": 1000.0, "x0x": 1.0}
try:
    du = json.loads(md.distribution("pandapower").read_text("direct_url.json") or "{}")
except Exception:
    du = {}
out = {"label": sys.argv[1] if len(sys.argv) > 1 else "", "pandapower_version": md.version("pandapower"),
       "direct_url": du, "input": p, "rows": []}
for vg in ("Yzn", "Dyn"):
    for zn in (0.0, 10.0):
        for sn in (1.0, 10.0, 100.0):
            row = {"vg": vg, "zn_ohm": zn, "net_sn_mva": sn}
            try:
                r = R.native_sc(P, vg, dict(p, sn_net=sn), zn, 0.0, "min")
                row["lv_z0_ohm"] = [r["lv"]["z0"].real, r["lv"]["z0"].imag]
                row["lv_ikss_ka"] = r["lv"]["ikss_ka"]
            except Exception as e:
                row["error"] = repr(e)[:300]
            out["rows"].append(row)
print(json.dumps(out, indent=1))
