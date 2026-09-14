"""Measure every declared-debt literal from tracked code, so each can be
registered under C21 with a value I computed rather than copied from prose.

'measure or cite, never guess' (Sec.15 rule 3). If a measurement does NOT
reproduce the paper's printed literal, that is a finding to fix in the paper,
not a value to register.
"""
from __future__ import annotations

import importlib.util as iu
import numpy as np


def _load(path, name):
    s = iu.spec_from_file_location(name, path)
    m = iu.module_from_spec(s)
    s.loader.exec_module(m)
    return m


stur = _load("tests/test_paper60_sturmian.py", "_p60stur")
pre = _load("tests/test_paper60_preconditioner.py", "_p60pre")
from geovac.sturmian_l2_encoding import fit_lambda_exponent  # noqa: E402


def fit(Ns, ys):
    return float(np.polyfit(np.log(Ns), np.log(ys), 1)[0])


print("=" * 64)
print("1. eq:blowup  Q exponents (paper: hydrogenic 1.19, sturmian 3.33)")
Z = 2.0
p_hy = fit_lambda_exponent([1, 2, 3, 4, 5], "hydrogenic", Z=Z)
p_st = fit_lambda_exponent([1, 2, 3, 4, 5], "sturmian", Z=Z)
print(f"   hydrogenic Q^{p_hy:.4f}   sturmian Q^{p_st:.4f}")

print("=" * 64)
print("2. SW cond exponent (paper: N^1.85, window) and L2 overlap (N^1.70)")
R = 2.0
nmaxes = [2, 3, 4, 5, 6, 7, 8]
Ns = [2 * n for n in nmaxes]
sw = [float(np.linalg.cond(stur._sw_metric_mom(n, R))) for n in nmaxes]
print(f"   SW cond seq {[f'{x:.1f}' for x in sw]}")
print(f"   SW exponent full N={Ns[0]}..{Ns[-1]}: N^{fit(Ns, sw):.4f}")
# try to find the window that gives ~1.85
for lo in range(len(nmaxes) - 2):
    for hi in range(lo + 2, len(nmaxes)):
        e = fit(Ns[lo:hi + 1], sw[lo:hi + 1])
        if 1.82 < e < 1.88:
            print(f"     window N={Ns[lo]}..{Ns[hi]}: N^{e:.4f}  <- ~1.85")

# L2 overlap cond exponent (the 3.0/5.8/13.9/32.2 sequence)
r, _ = stur._grid()
l2 = []
Nl2 = [2, 3, 5, 8]
for N in Nl2:
    fs = [stur._sturmian_s(r, n) for n in range(1, N + 1)]
    S = np.array([[np.trapezoid(fs[i] * fs[j] * r * r, r) for j in range(N)]
                  for i in range(N)])
    l2.append(float(np.linalg.cond(S)))
print(f"   L2 overlap cond seq {[f'{x:.2f}' for x in l2]}  exponent N^{fit(Nl2, l2):.4f}")

print("=" * 64)
print("3. Gaussian ratio exponents (paper: ratio-1.6 -> N^6, ratio-3 -> N^1.4)")
for ratio in (1.6, 2.0, 3.0):
    ga = [float(np.linalg.cond(stur._gaussian_metric(n, R, ratio=ratio)))
          for n in nmaxes]
    print(f"   ratio={ratio}: cond seq {[f'{x:.3g}' for x in ga]}  "
          f"exponent N^{fit(Ns, ga):.4f}")

print("=" * 64)
print("4. water A_1 exponent (paper: N^1.97, cond 19.9->698 over N=6..36)")
wat_ns = [6, 12, 24, 48, 96]
wat = [float(np.linalg.cond(pre._water_A1(n))) for n in wat_ns]
Nw = [2 * n for n in wat_ns]     # paper N = 2n (raw water A1)
print(f"   water A1 cond seq {[f'{x:.1f}' for x in wat]}")
print(f"   water A1 exponent N={Nw[0]}..{Nw[-1]}: N^{fit(Nw, wat):.4f}")
print(f"   over N=12..192 only: N^{fit(Nw, wat):.4f}")
