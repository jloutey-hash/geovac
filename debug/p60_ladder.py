"""Phase B -- converged exponent ladder for the Paper-60 secular 1-norm.

Runs the s+p+d+f family nmax = 7..NMAX on a K-dependent box R_MAX = c * nmax^2
(graded mesh), recording every 1-norm leg.  Appends to a JSON so partial runs
are usable.
"""
import json, os, sys, time
import debug.p60_engine as E

c = float(sys.argv[1])
nmin, nmax_hi = int(sys.argv[2]), int(sys.argv[3])
npts = int(sys.argv[4]) if len(sys.argv) > 4 else 24000
Z = float(sys.argv[5]) if len(sys.argv) > 5 else 2.0
lmax = int(sys.argv[6]) if len(sys.argv) > 6 else 3
tag = sys.argv[7] if len(sys.argv) > 7 else ""
OUT = f"debug/data/p60_ladder_c{c:g}_N{npts}_Z{Z:g}_l{lmax}{tag}.json"

rows = json.load(open(OUT)) if os.path.exists(OUT) else []
done = {(r["nmax"], r["c"]) for r in rows}
for n in range(nmin, nmax_hi + 1):
    if (n, c) in done:
        continue
    cts = E.family(n, lmax)
    rmax = c * n * n
    E.set_grid(rmax, npts, "grade", 2.0)
    t0 = time.time()
    m = E.norms(cts, Z=Z)
    m.update(c=c, rmax=rmax, npts=npts, nmax=n, lmax=lmax, Z=Z,
             wall=time.time() - t0)
    rows.append(m)
    print(f"c={c:g} Z={Z:g} lmax={lmax} nmax={n:2d} K={m['K']:4d} Rmax={rmax:7.1f} "
          f"Mtot={m['M_total']:12.6f} Mdiag={m['M_diag']:12.6f} Moff={m['M_off']:12.6f} "
          f"Tpfull={m['Tp_full']:12.6f} Tpdiag={m['Tp_diag']:10.6f} T0={m['T0']:11.6f} "
          f"E={m['E']:10.6f} [{m['wall']:.0f}s]")
    sys.stdout.flush()
    json.dump(rows, open(OUT, "w"), indent=1)
