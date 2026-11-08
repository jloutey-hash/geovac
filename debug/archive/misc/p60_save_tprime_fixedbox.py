"""Save T' on a FIXED box (e.g. the production R_MAX=60) for the same families,
so the converged vs production comparison is like-for-like in everything but the box."""
import os, sys, time
import numpy as np
import debug.p60_engine as E
import geovac.sturmian_secular as S

rmax = float(sys.argv[1]); nlo = int(sys.argv[2]); nhi = int(sys.argv[3])
npts = int(sys.argv[4]) if len(sys.argv) > 4 else 24000
lmax = int(sys.argv[5]) if len(sys.argv) > 5 else 3
kind = sys.argv[6] if len(sys.argv) > 6 else "grade"
for n in range(nlo, nhi + 1):
    out = f"debug/data/p60_Tp_l{lmax}_n{n}_R{rmax:g}_N{npts}_{kind}.npz"
    if os.path.exists(out):
        print("skip", out); continue
    cts = E.family(n, lmax)
    E.set_grid(rmax, npts, kind, 2.0)
    t0 = time.time()
    cfgs = S.build_configs(cts)
    Tp = S.build_Tprime(cfgs)
    np.savez_compressed(out, Tp=Tp, Rnu=np.array([c.Rnu for c in cfgs]),
                        lab=np.array([(c.l, c.na, c.nb) for c in cfgs]),
                        rmax=rmax, npts=npts, nmax=n, lmax=lmax)
    print(f"rmax={rmax:g} nmax={n} lmax={lmax} K={len(cfgs)} [{time.time()-t0:.0f}s]")
    sys.stdout.flush()
