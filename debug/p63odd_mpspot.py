"""mp30 spot values of e^{2x} Phi_Lx[f_p E](x), L = 0, 2, x = 1, 3, 10 (referee for the float Laguerre route)."""
import json, sys, time
from pathlib import Path
import numpy as np, mpmath as mp
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import p63odd_engine as EN
mp.mp.dps = 30
out = []
for L in (0, 2):
    for x in (1.0, 3.0, 10.0):
        t0 = time.time()
        v = EN.phix_mp(x, L)
        f = float(EN.phix_f(np.array([x]), L)[0])
        out.append({"L": L, "x": x, "mp30": mp.nstr(v, 25), "float": f, "rel": float(abs(f / v - 1)), "seconds": time.time() - t0})
        print(out[-1], flush=True)
json.dump(out, open(HERE / "data" / "p63odd_mpspot.json", "w"), indent=1)
