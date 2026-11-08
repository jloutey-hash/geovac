"""Measure the uniform blockdiag(P,P) control at the PAPER's own N grid.

The paper table (L1310-1327) uses N = 12, 48, 192 with raw cond(A_1) =
183.0 / 2696.1 / 41699.7, and the sentence says the uniform control runs
"2766 -> 42008 over the same range, an exponent of N^1.95".

The claims reviewer flags that 2766 -> 42008 over N=12..192 (16x) gives N^0.98,
not N^1.95, so the displayed low value is stale.  Rather than accept its
proposed 729.2, measure the uniform column at the SAME three N points the raw
table uses, so the displayed pair is unambiguous over "the same range" and its
exponent actually matches what is printed.

The test's `_water_A1(n)` uses n = N/2 (raw at n=24 is 2696 = table N=48; raw
at n=96 is 41700 = table N=192).  So paper N=12/48/192 = n=6/24/96.
"""
from __future__ import annotations

import importlib.util as iu
import numpy as np

_s = iu.spec_from_file_location("_pre", "tests/test_paper60_preconditioner.py")
_m = iu.module_from_spec(_s)
_s.loader.exec_module(_m)
_water_A1 = _m._water_A1
tridiag = _m.tridiag
inv_sqrt = _m.inv_sqrt

# paper N -> test n = N/2
GRID = [(12, 6), (48, 24), (192, 96)]
rows = []
for Npaper, n in GRID:
    A = _water_A1(n)
    P = np.zeros_like(A)
    P[:n, :n] = tridiag(n)
    P[n:, n:] = tridiag(n)      # uniform blockdiag(P,P)
    Pis = inv_sqrt(P)
    raw = np.linalg.cond(A)
    uni = np.linalg.cond(Pis @ A @ Pis)
    rows.append((Npaper, raw, uni))

print(f"{'N(paper)':>9} {'raw cond':>14} {'uniform cond':>14}")
for Np, raw, uni in rows:
    print(f"{Np:>9} {raw:>14.4f} {uni:>14.4f}")

lnN = np.log([r[0] for r in rows])
e_raw = float(np.polyfit(lnN, np.log([r[1] for r in rows]), 1)[0])
e_uni = float(np.polyfit(lnN, np.log([r[2] for r in rows]), 1)[0])
print()
print(f"raw exponent over N=12..192   : N^{e_raw:.4f}")
print(f"uniform exponent over N=12..192: N^{e_uni:.4f}")
print()
lo, hi = rows[0][2], rows[-1][2]
print(f"uniform displayed pair over the full range: {lo:.1f} -> {hi:.1f}")
print(f"  two-point slope over N=12..192 (16x): "
      f"N^{np.log(hi/lo)/np.log(16):.4f}")
print(f"the STALE printed low value 2766 is the N=48 point: "
      f"{rows[1][2]:.1f} (matches? {abs(rows[1][2]-2766)<5})")
