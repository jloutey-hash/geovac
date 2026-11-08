"""DIAGNOSTIC: does the fixed-l_max Goscinskian family converge?

Ground-truth calibration.  The l_max = 0 (s-only) sector has a KNOWN exact
answer -- the helium s-limit, -2.879028767 Ha -- so the same free-floor fit
that produced Paper 60's 6.44 mHa spdf floor can be run against a case where
the true asymptote is not in doubt.

Usage: python debug/p60_slimit_probe.py NMIN NMAX [NPTS] [BOXC]
"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import debug.p60_engine as E

S_LIMIT = -2.879028767          # He l=0-only exact (literature)
nmin, nmax_hi = int(sys.argv[1]), int(sys.argv[2])
npts = int(sys.argv[3]) if len(sys.argv) > 3 else 24000
boxc = float(sys.argv[4]) if len(sys.argv) > 4 else 5.0
OUT = f"debug/data/p60_slimit_N{npts}_c{boxc:g}.json"
rows = json.load(open(OUT)) if os.path.exists(OUT) else []
done = {r["nmax"] for r in rows}

for n in range(nmin, nmax_hi + 1):
    if n in done:
        continue
    cts = E.family(n, 0)
    rmax = boxc * n * n
    E.set_grid(rmax, npts, "grade", 2.0)
    t0 = time.time()
    m = E.norms(cts, Z=2.0)
    m.update(nmax=n, rmax=rmax, npts=npts, lmax=0, wall=time.time() - t0)
    rows.append(m)
    print("nmax=%2d K=%4d rmax=%7.1f E=%.7f  gap_to_s_limit=%8.4f mHa  M=%10.4f [%.0fs]"
          % (n, m["K"], rmax, m["E"], (m["E"] - S_LIMIT) * 1000, m["M_total"], m["wall"]))
    sys.stdout.flush()
    json.dump(rows, open(OUT, "w"), indent=1)
