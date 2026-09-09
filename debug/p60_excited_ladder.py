"""Excited-state accuracy at fixed encoding cost.

||M||_1 is a property of M, not of which root you extract -- so the ground state
and 2^1S cost the SAME to block-encode.  What differs is the accuracy delivered.
This measures that, on the full s+p+d+f family.
"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import debug.p60_engine as E
import geovac.sturmian_secular as S

GROUND = -2.903724377        # He 1^1S exact non-relativistic
EXC1 = -2.145974046          # He 2^1S exact non-relativistic
CHEM = 1.5936014616          # 1 kcal/mol in mHa

lmax = int(sys.argv[1])
OUT = "debug/data/p60_excited_l%d.json" % lmax
rows = json.load(open(OUT)) if os.path.exists(OUT) else []
done = {r["nmax"] for r in rows}
for nmax in [int(x) for x in sys.argv[2:]]:
    if nmax in done:
        continue
    E.set_grid(max(80.0, 5.0 * nmax * nmax), 24000, "grade", 2.0)
    t0 = time.time()
    cfgs = S.build_configs(E.family(nmax, lmax))
    M = S.build_M(cfgs, Z=2.0)
    p = np.sort(np.linalg.eigvalsh(M))[::-1]
    e0, e1 = -p[0] ** 2 / 2, -p[1] ** 2 / 2
    r = dict(nmax=nmax, K=len(cfgs), M1=float(np.abs(M).sum()),
             e0=e0, e1=e1, gap0=(e0 - GROUND) * 1000, gap1=(e1 - EXC1) * 1000,
             wall=time.time() - t0)
    rows.append(r)
    print("nmax=%2d K=%4d ||M||1=%9.2f | gnd gap=%7.3f mHa (%.2fx chem) | 2^1S gap=%7.3f mHa (%.2fx chem) [%.0fs]"
          % (nmax, r["K"], r["M1"], r["gap0"], r["gap0"] / CHEM, r["gap1"], r["gap1"] / CHEM, r["wall"]))
    sys.stdout.flush()
    json.dump(rows, open(OUT, "w"), indent=1)
