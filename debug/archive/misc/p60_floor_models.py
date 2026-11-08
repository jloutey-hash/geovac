"""JOB A falsifier: is the floor demanded by the data, or an artifact of a
3-parameter fit applied to slow convergence?  Compare four models on the same
windows.  M1 has a free floor; M2/M3 force the limit to ZERO with progressively
slower decay.  If a zero-limit model fits comparably, the floor is not established.
"""
import json, sys
import numpy as np
from scipy.optimize import curve_fit

CHEM = 1.5936014616
lmax = int(sys.argv[1]) if len(sys.argv) > 1 else 3
rows = sorted(json.load(open("debug/data/p60_excited_l%d.json" % lmax)), key=lambda r: r["K"])
K = np.array([r["K"] for r in rows], float)
G = {0: np.array([r["gap0"] for r in rows]), 1: np.array([r["gap1"] for r in rows])}

MODELS = {
    "M1 c+bK^-q   (free floor, 3p)": (lambda k, c, b, q: c + b * k ** (-q), [1.0, 10.0, 0.8]),
    "M2 bK^-q     (zero limit, 2p)": (lambda k, b, q: b * k ** (-q), [10.0, 0.1]),
    "M3 b(lnK)^-s (zero limit, 2p)": (lambda k, b, s: b * np.log(k) ** (-s), [10.0, 0.5]),
    "M4 c+b/lnK   (free floor, 2p)": (lambda k, c, b: c + b / np.log(k), [1.0, 5.0]),
    "M5 b/(lnK)+d/(lnK)^2 (zero,2p)": (lambda k, b, d: b / np.log(k) + d / np.log(k) ** 2, [5.0, 5.0]),
}

out = {}
for root in (0, 1):
    g = G[root]
    print("=" * 86)
    print("ROOT %d %s" % (root, "ground (1^1S)" if root == 0 else "2^1S"))
    for wname, sl in (("ALL K=74..452 (n=10)", slice(0, 10)),
                      ("TAIL K=202..452 (n=6)", slice(4, 10)),
                      ("TAIL K=290..452 (n=4)", slice(6, 10))):
        kk, gg = K[sl], g[sl]
        print("  %s" % wname)
        for name, (f, p0) in MODELS.items():
            try:
                popt, _ = curve_fit(f, kk, gg, p0=p0, maxfev=400000)
                res = gg - f(kk, *popt)
                rms = float(np.sqrt(np.mean(res ** 2)))
                lim = float(popt[0]) if name.startswith(("M1", "M4")) else 0.0
                print("    %-32s rms=%.3e  limit=%+8.4f mHa (%.2fx chem)  params=%s"
                      % (name, rms, lim, lim / CHEM,
                         np.array2string(popt, precision=4, suppress_small=False)))
                out.setdefault(str(root), {}).setdefault(wname, {})[name] = dict(
                    rms=rms, limit=lim, params=[float(x) for x in popt])
            except Exception as e:
                print("    %-32s FAILED %s" % (name, e))
        print()
json.dump(out, open("debug/data/p60_floor_models_l%d.json" % lmax, "w"), indent=1)
print("wrote debug/data/p60_floor_models_l%d.json" % lmax)
