"""Measure the THREE water A1 conditioning columns, so the paper's control
sentence can be written from data instead of from the blind comparison.

Column 1: raw cond(A).
Column 2: UNIFORM blockdiag(T,T) -- commutes with the rotation (I2 (x) T vs
          V (x) I), so it is a SELECTIVITY control, not a rotation control.
Column 3: SELECTIVE blockdiag(T,I) in the UNROTATED frame -- the tracked
          discriminating control.
Column 4: SELECTIVE blockdiag(T,I) in the ROTATED frame -- the lever.
"""
from __future__ import annotations
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import importlib.util as _iu
_s = _iu.spec_from_file_location("_p60pre", "tests/test_paper60_preconditioner.py")
_m = _iu.module_from_spec(_s); _s.loader.exec_module(_m)
_water_A1 = _m._water_A1
_null_direction_rotation = _m._null_direction_rotation
tridiag = _m.tridiag
inv_sqrt = _m.inv_sqrt

NS = (12, 24, 48, 96)
rows = []
for n in NS:
    A = _water_A1(n)
    Q = np.kron(_null_direction_rotation(), np.eye(n))
    Pu = np.zeros_like(A); Pu[:n, :n] = tridiag(n); Pu[n:, n:] = tridiag(n)
    Ps = np.zeros_like(A); Ps[:n, :n] = tridiag(n); Ps[n:, n:] = np.eye(n)
    Ui, Si = inv_sqrt(Pu), inv_sqrt(Ps)
    rows.append((n,
        np.linalg.cond(A),
        np.linalg.cond(Ui @ A @ Ui),
        np.linalg.cond(Si @ A @ Si),
        np.linalg.cond(Si @ (Q.T @ A @ Q) @ Si)))

print(f"{'n':>4} {'raw':>12} {'uniform':>12} {'sel/unrot':>12} {'sel/rot':>10}")
for r in rows:
    print(f"{r[0]:>4} {r[1]:>12.4g} {r[2]:>12.4g} {r[3]:>12.4g} {r[4]:>10.4g}")

ln = np.log(np.array(NS, float))
for j, name in ((1, "raw"), (2, "uniform"), (3, "sel/unrot"), (4, "sel/rot")):
    y = np.log(np.array([r[j] for r in rows], float))
    print(f"exponent {name:>10}: N^{np.polyfit(ln, y, 1)[0]:.3f}")
print()
print("uniform rotated == uniform unrotated?")
for n in NS:
    A = _water_A1(n); Q = np.kron(_null_direction_rotation(), np.eye(n))
    Pu = np.zeros_like(A); Pu[:n, :n] = tridiag(n); Pu[n:, n:] = tridiag(n)
    Ui = inv_sqrt(Pu)
    a = np.linalg.cond(Ui @ A @ Ui); b = np.linalg.cond(Ui @ (Q.T @ A @ Q) @ Ui)
    print(f"  n={n:>3}  {a:.6g} vs {b:.6g}   rel diff {abs(a-b)/a:.2e}")
