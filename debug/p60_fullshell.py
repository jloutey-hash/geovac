"""Robustness: a DIFFERENT config family -- full hydrogenic shells (lmax = nmax-1)."""
import sys, time, math
import numpy as np
import debug.p60_engine as E
import geovac.sturmian_secular as S

rows = []
for n in range(5, int(sys.argv[1]) + 1):
    cts = S.gen_configs(n - 1, {l: n for l in range(n)})
    E.set_grid(3.0 * n * n, 24000, "grade", 2.0)
    t0 = time.time()
    m = E.norms(cts)
    rows.append(m)
    sl = ""
    if len(rows) > 1:
        a, b = rows[-2], rows[-1]
        lk = math.log(b["K"]) - math.log(a["K"])
        sl = " ".join(f"{k}={((math.log(b[k])-math.log(a[k]))/lk):.4f}"
                      for k in ("M_total", "M_off", "Tp_full", "T0"))
    print(f"nmax={n} lmax={n-1} K={m['K']:4d} Mtot={m['M_total']:11.5f} "
          f"Moff={m['M_off']:11.5f} Tpfull={m['Tp_full']:11.5f} T0={m['T0']:10.5f} "
          f"[{time.time()-t0:.0f}s]  local: {sl}", flush=True)
