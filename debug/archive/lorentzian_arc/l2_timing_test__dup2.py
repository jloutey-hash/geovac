"""Timing test for L2 panel computation at high precision."""
from __future__ import annotations
import sys
import time
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJ))

import mpmath
from mpmath import mp, mpf, pi, log
from geovac.central_fejer_su2 import T_n_via_sum_rule


def compute_gamma_n(n: int, dps: int) -> mpmath.mpf:
    mp.dps = dps
    Z = mpf(n * (n + 1)) / 2
    T = T_n_via_sum_rule(n, prec=dps)
    return pi - 4 * T / (pi * Z)


def time_computation(n: int, dps: int):
    t0 = time.time()
    g = compute_gamma_n(n, dps=dps)
    h = mpf(n) * g - 4/pi * log(mpf(n))
    dt = time.time() - t0
    return dt, g, h


if __name__ == "__main__":
    # Time a few small n at dps=400 to extrapolate full panel cost
    for dps in [120, 250, 400]:
        print(f"\n--- dps = {dps} ---")
        for n in [256, 512, 1024]:
            dt, g, h = time_computation(n, dps)
            print(f"  n={n:>5d}: {dt:.2f}s  h={mpmath.nstr(h, 12)}")
