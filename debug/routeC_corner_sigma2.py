"""Diagnostic: does the substitution rho=sigma^2 in the Duffy corner make the
corner-triangle integral of J converge SPECTRALLY (no subtraction needed)?

Claim to test: at fixed Duffy angle alpha, J(rho*alpha, rho*(1-alpha)) is
ANALYTIC in sigma=sqrt(rho) -- a pure integer+half-integer power series in rho
with NO log terms (the k-integrand decays exponentially at large k for every
order in the rho-expansion, so no large-k cutoff / ln(rho) mechanism fires;
leading term rho^{3/2}, cf. routeC_corner_asymptotics.py).

If true, T1 = int_{s+t<=delta} J ds dt, written in Duffy (s=rho*alpha,
t=rho*(1-alpha), ds dt = rho drho dalpha) and then rho=sigma^2 (drho=2 sigma
dsigma), has integrand

    2 sigma^3 * J(sigma^2 alpha, sigma^2 (1-alpha))   on (sigma,alpha) in [0,sqrt(delta)]x[0,1]

which is ANALYTIC in sigma (starts at sigma^6), so plain Gauss-Legendre
converges spectrally -- unlike the rho^{3/2}-subtracted remainder (~rho^{5/2},
still non-analytic, GL only algebraic).

Compares three treatments of the SAME T1:
  (A) raw Duffy + GL in rho          -- baseline (non-analytic rho^{3/2}, slow)
  (B) rho^{3/2}-subtracted (corner_v2) -- one term removed (rho^{5/2}, faster-but-capped)
  (C) Duffy + rho=sigma^2 + GL        -- predicted spectral

Usage: python routeC_corner_sigma2.py [dps] [delta] [kpdeg]
"""
from __future__ import annotations
import sys
import time
import mpmath as mp

sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_probe9 import J_fixed
from routeC_probe8 import k_grid_sinh_paneled
from routeC_hp_evaluator import std_gl_nodes
from routeC_corner_v2 import corner_triangle_integral_subtracted


def T1_sigma2(delta, deg_sigma, deg_alpha, knodes):
    """(C) Duffy + rho=sigma^2, GL x GL. Integrand 2 sigma^3 J on
    (sigma,alpha) in [0,sqrt(delta)] x [0,1]. Uses s<->t (alpha<->1-alpha)
    symmetry: integrate alpha in [0,1/2], double."""
    std_s = std_gl_nodes(deg_sigma, mp.mp.prec)
    std_a = std_gl_nodes(deg_alpha, mp.mp.prec)
    sig_max = mp.sqrt(delta)
    half_s = sig_max / 2
    tot = mp.mpf(0)
    for xs, ws in std_s:
        sigma = half_s * (xs + 1)
        wsig = half_s * ws
        rho = sigma * sigma
        pref = wsig * 2 * sigma ** 3
        # alpha in [0,1/2] doubled by symmetry
        for xa, wa in std_a:
            alpha = (xa + 1) / 4            # maps [-1,1] -> [0,1/2]
            wal = wa / 4
            s = rho * alpha
            t = rho * (1 - alpha)
            tot += pref * wal * 2 * J_fixed(s, t, knodes)
    return tot


def T1_raw_duffy(delta, deg_rho, deg_alpha, knodes):
    """(A) raw Duffy + GL in rho (no substitution, no subtraction)."""
    std_r = std_gl_nodes(deg_rho, mp.mp.prec)
    std_a = std_gl_nodes(deg_alpha, mp.mp.prec)
    half_r = delta / 2
    tot = mp.mpf(0)
    for xr, wr in std_r:
        rho = half_r * (xr + 1)
        wrho = half_r * wr
        for xa, wa in std_a:
            alpha = (xa + 1) / 4
            wal = wa / 4
            s = rho * alpha
            t = rho * (1 - alpha)
            tot += wrho * wal * 2 * rho * J_fixed(s, t, knodes)
    return tot


def main():
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 45
    mp.mp.dps = dps
    delta = mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf('0.05')
    kpdeg = int(sys.argv[3]) if len(sys.argv) > 3 else 5

    knodes = k_grid_sinh_paneled(mp.mpf(16), mp.mpf(2), kpdeg)
    print(f"dps={dps} delta={delta} k_M={len(knodes)}", flush=True)

    print("\n(C) Duffy + rho=sigma^2 + GL  [predicted SPECTRAL]")
    prev = None
    for d in [2, 3, 4, 5, 6]:
        t0 = time.time()
        v = T1_sigma2(delta, d, d, knodes)
        diff = "" if prev is None else mp.nstr(abs(v - prev), 4)
        n = 3 * 2 ** (d - 1)
        print(f"  deg={d} (n={n}^2)  T1={mp.nstr(v, dps-5)}  diff={diff}  ({time.time()-t0:.1f}s)", flush=True)
        prev = v
    C_val = prev

    print("\n(A) raw Duffy + GL in rho  [baseline, non-analytic rho^{3/2}]")
    prev = None
    for d in [2, 3, 4, 5, 6]:
        v = T1_raw_duffy(delta, d, d, knodes)
        diff = "" if prev is None else mp.nstr(abs(v - prev), 4)
        print(f"  deg={d}  T1={mp.nstr(v, dps-5)}  diff={diff}", flush=True)
        prev = v

    print("\n(B) rho^{3/2}-subtracted (corner_v2)  [one term removed, rho^{5/2}]")
    prev = None
    for d in [2, 3, 4, 5, 6]:
        v, rem, ab = corner_triangle_integral_subtracted(delta, d, d, knodes)
        diff = "" if prev is None else mp.nstr(abs(v - prev), 4)
        print(f"  deg={d}  T1={mp.nstr(v, dps-5)}  diff={diff}", flush=True)
        prev = v
    B_val = prev

    print(f"\n(C) - (B) at top degree: {mp.nstr(abs(C_val - B_val), 6)}")


if __name__ == '__main__':
    main()
