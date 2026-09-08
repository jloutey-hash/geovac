"""Resolution-axis convergence at the CONVERGED box (c=3): does N matter?"""
import sys, time
import debug.p60_engine as E

n = int(sys.argv[1])
cts = E.family(n)
print(f"# nmax={n} K={len(cts)} c=3 rmax={3*n*n}")
print(f"{'kind':>6} {'npts':>7} {'M_total':>13} {'M_diag':>13} {'M_off':>13} {'Tp_full':>13} {'s':>5}")
for kind, npts in [("grade", 12000), ("grade", 24000), ("grade", 48000),
                   ("grade", 96000), ("uni", 24000), ("uni", 96000)]:
    E.set_grid(3 * n * n, npts, kind, 2.0)
    t0 = time.time()
    m = E.norms(cts)
    print(f"{kind:>6} {npts:>7} {m['M_total']:>13.7f} {m['M_diag']:>13.7f} "
          f"{m['M_off']:>13.7f} {m['Tp_full']:>13.7f} {time.time()-t0:>5.0f}")
    sys.stdout.flush()
