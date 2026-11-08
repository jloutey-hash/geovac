"""Phase A -- K-dependent box rule for the Paper-60 secular 1-norm.

For each family nmax (K = 74..340), sweep R_MAX = c * nmax^2 and find the c at
which every 1-norm leg is stable.  Graded mesh (validated against the production
uniform mesh to 7 digits at 1/12 the points).
"""
import json, math, sys, time
import numpy as np
import debug.p60_engine as E

CS = [None, 1.0, 2.0, 3.0, 5.0]     # None -> the production R_MAX = 60
NPTS = int(sys.argv[2]) if len(sys.argv) > 2 else 24000
OUT = f"debug/data/p60_boxrule_n{sys.argv[1]}.json"

nmax = int(sys.argv[1])
cts = E.family(nmax)
rows = []
print(f"# nmax={nmax} K={len(cts)} npts={NPTS} graded p=2")
print(f"{'c':>6} {'Rmax':>8} {'M_total':>13} {'M_diag':>13} {'M_off':>13} "
      f"{'Tp_full':>13} {'Tp_diag':>12} {'T0':>12} {'E':>11} {'s':>5}")
for c in CS:
    rmax = 60.0 if c is None else c * nmax * nmax
    E.set_grid(rmax, NPTS, "grade", 2.0)
    t0 = time.time()
    m = E.norms(cts)
    m.update(c=c, rmax=rmax, npts=NPTS, nmax=nmax, wall=time.time() - t0)
    rows.append(m)
    print(f"{str(c):>6} {rmax:>8.1f} {m['M_total']:>13.6f} {m['M_diag']:>13.6f} "
          f"{m['M_off']:>13.6f} {m['Tp_full']:>13.6f} {m['Tp_diag']:>12.6f} "
          f"{m['T0']:>12.6f} {m['E']:>11.6f} {m['wall']:>5.0f}")
    sys.stdout.flush()
    with open(OUT, "w") as fh:
        json.dump(rows, fh, indent=1)
