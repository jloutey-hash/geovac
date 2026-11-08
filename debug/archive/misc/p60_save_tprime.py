"""Save the converged T' matrix (and R_nu) for each family, so every downstream
analysis -- Z sweep, T' vs offdiag split, local slopes -- is pure post-processing.

Box:  R_MAX = c * nmax^2 (default c=3, established converged to ~1e-7 by
debug/p60_boxrule.py).  Mesh: graded p=2.
"""
import os, sys, time
import numpy as np
import debug.p60_engine as E
import geovac.sturmian_secular as S

c = float(sys.argv[1]); nlo = int(sys.argv[2]); nhi = int(sys.argv[3])
npts = int(sys.argv[4]) if len(sys.argv) > 4 else 24000
lmax = int(sys.argv[5]) if len(sys.argv) > 5 else 3
for n in range(nlo, nhi + 1):
    out = f"debug/data/p60_Tp_l{lmax}_n{n}_c{c:g}_N{npts}.npz"
    if os.path.exists(out):
        print("skip", out); continue
    cts = E.family(n, lmax)
    E.set_grid(c * n * n, npts, "grade", 2.0)
    t0 = time.time()
    cfgs = S.build_configs(cts)
    Tp = S.build_Tprime(cfgs)
    Rnu = np.array([cf.Rnu for cf in cfgs])
    lab = np.array([(cf.l, cf.na, cf.nb) for cf in cfgs])
    np.savez_compressed(out, Tp=Tp, Rnu=Rnu, lab=lab, rmax=c * n * n, npts=npts,
                        nmax=n, lmax=lmax)
    print(f"nmax={n} lmax={lmax} K={len(cfgs)} rmax={c*n*n:.0f} N={npts} "
          f"[{time.time()-t0:.0f}s] -> {out}")
    sys.stdout.flush()
