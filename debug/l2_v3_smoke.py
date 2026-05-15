"""Quick smoke test of v3 pipeline using small n_max points only."""
from __future__ import annotations
import json
import sys
import time
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJ))

import mpmath
from mpmath import mp, mpf, log, pi
from geovac.central_fejer_su2 import T_n_via_sum_rule

# Use only small n for smoke test - should run in seconds
mp.dps = 200

def compute_h(n):
    Z = mpf(n * (n + 1)) / 2
    T = T_n_via_sum_rule(n, prec=200)
    g = pi - 4 * T / (pi * Z)
    return mpf(n) * g - (4 / pi) * log(mpf(n))

# Use Track 4's existing high-precision intermediate panel
with open(PROJ / "debug" / "data" / "l2_constant_c_precision_push.json") as f:
    data = json.load(f)
t4_panel = []
for p in data["panel_h_values"]:
    if len(p["h_250dps"]) > 200:
        t4_panel.append({"n": p["n"], "h_n": mpf(p["h_250dps"])})

print(f"Loaded {len(t4_panel)} Track 4 high-precision points")

# Also compute a few small high-precision points for smoke test
print("Computing small high-precision points...")
for n in [32, 64, 128, 256]:
    t0 = time.time()
    h = compute_h(n)
    t4_panel.append({"n": n, "h_n": h})
    print(f"  n={n}: {mpmath.nstr(h, 30)} ({time.time()-t0:.2f}s)")

t4_panel.sort(key=lambda d: d["n"])
print(f"Total panel: {len(t4_panel)} points, n values: {[p['n'] for p in t4_panel]}")

# Run Richardson at several K
def build_basis_matrix(n_values, K):
    rows = []
    for n in n_values:
        ln = log(mpf(n))
        inv_n = 1 / mpf(n)
        row = [mpf(1)]
        cur_inv = inv_n
        for k in range(1, K + 1):
            row.append(ln * cur_inv)
            row.append(cur_inv)
            cur_inv *= inv_n
        rows.append(row)
    return rows

def fit_subleading(panel, K):
    n_values = [d["n"] for d in panel]
    h_values = [d["h_n"] for d in panel]
    rows = build_basis_matrix(n_values, K)
    A = mpmath.matrix(rows)
    y = mpmath.matrix([[v] for v in h_values])
    n_pts = len(rows)
    p = len(rows[0])
    if n_pts == p:
        x = A ** (-1) * y
    else:
        AT = A.T
        x = (AT * A) ** (-1) * (AT * y)
    coeffs = [x[i, 0] for i in range(p)]
    fitted = [sum(c * r for c, r in zip(coeffs, row)) for row in rows]
    resid = [h - f for h, f in zip(h_values, fitted)]
    return {"K": K, "b_0": coeffs[0], "max_resid": max(abs(r) for r in resid)}

print()
print(f"Richardson tower with {len(t4_panel)} points:")
for K in range(3, 12):
    n_params = 2 * K + 1
    if n_params > len(t4_panel):
        print(f"  K={K} skipped (need {n_params}, have {len(t4_panel)})")
        continue
    try:
        fit = fit_subleading(t4_panel, K)
        print(f"  K={K:>2d} ({n_params:>2d} params): c={mpmath.nstr(fit['b_0'], 30)}, "
              f"max_resid={mpmath.nstr(fit['max_resid'], 6)}")
    except Exception as e:
        print(f"  K={K} fit failed: {e}")
