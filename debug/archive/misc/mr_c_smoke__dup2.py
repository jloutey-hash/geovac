"""Smoke test for MR-C: small n at moderate precision."""
import sys
from pathlib import Path
PROJ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJ))

import mpmath
from mpmath import mp, mpf, log, pi
from geovac.central_fejer_su2 import T_n_via_sum_rule

mp.dps = 60
N_VALUES = [16, 32, 64, 128, 256]

four_over_pi = 4 / pi
print(f"dps={mp.dps}  4/pi = {mpmath.nstr(four_over_pi, 20)}")
panel = []
for n in N_VALUES:
    Z = mpf(n*(n+1))/2
    T = T_n_via_sum_rule(n, prec=60)
    g = pi - 4*T/(pi*Z)
    h = mpf(n)*g - four_over_pi*log(mpf(n))
    panel.append((n, g, h))
    print(f"  n={n:4d}  gamma={mpmath.nstr(g, 18)}  h={mpmath.nstr(h, 18)}")

# Quick LS fit with K=4 -> 9 params, 5 points: underdetermined
# So use K=2: 5 params, 5 points, square system
import sys
sys.path.insert(0, str(PROJ / "debug"))
from mr_c_l2_subleading import fit_subleading

# Build a tiny panel-of-dicts
pdicts = [{"n": n, "gamma_n": g, "h_n": h} for n, g, h in panel]
fit = fit_subleading(pdicts, K=2)
print(f"\nK=2 fit:")
print(f"  c (b_0)   = {mpmath.nstr(fit['b_0_extracted'], 30)}")
for k, (a, b) in enumerate(fit['a_b_pairs'], 1):
    print(f"  a_{k}={mpmath.nstr(a, 14)}  b_{k}={mpmath.nstr(b, 14)}")
print(f"  max residual = {mpmath.nstr(max(abs(r) for r in fit['residuals']), 8)}")
