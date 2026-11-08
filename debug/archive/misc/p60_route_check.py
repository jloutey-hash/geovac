"""Cross-check: graded mesh vs production uniform mesh at the same R_MAX."""
import time
import numpy as np
import debug.p60_engine as E

cts = E.family(7)
print(f"{'route':>7} {'rmax':>6} {'npts':>7} {'p':>4} {'M_total':>12} {'M_diag':>12} "
      f"{'M_off':>12} {'Tp_full':>12} {'E':>11} {'s':>5}")
for kind, rmax, npts, p in [
    ("uni", 60, 18000, 0), ("uni", 60, 72000, 0), ("uni", 60, 288000, 0),
    ("grade", 60, 6000, 2.0), ("grade", 60, 12000, 2.0), ("grade", 60, 24000, 2.0),
    ("grade", 60, 12000, 1.5), ("grade", 60, 24000, 1.5),
]:
    E.set_grid(rmax, npts, kind, p)
    t0 = time.time()
    m = E.norms(cts)
    print(f"{kind:>7} {rmax:>6} {npts:>7} {p:>4} {m['M_total']:>12.6f} {m['M_diag']:>12.6f} "
          f"{m['M_off']:>12.6f} {m['Tp_full']:>12.6f} {m['E']:>11.6f} {time.time()-t0:>5.0f}")
