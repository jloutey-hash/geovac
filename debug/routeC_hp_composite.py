"""High-precision COMPOSITE evaluator for T2 (Paper 59 collinear three-centre
observable), targeting >=40 cross-validated digits.

    T2 = (8/pi) int_{[0,1]^2} J(s,t) ds dt,   J symmetric in s,t
    J(s,t) = int_0^inf dk j0(k(s+t)) P(s,k) P(t,k)

State of play (from debug/sprint_hp_evaluator_memo.md, Continuations 1-2):
  * The ONLY non-analyticity of the outer integrand is the (0,0) corner, where
    J ~ rho^{3/2} (rho=s+t), degree-3/2, C^1-not-C^2 -- verified NOT a log
    (routeC_corner_asymptotics.py: (J-S)/J = O(rho) clean, factor-10 scaling).
  * Everywhere else J(s,t) is ANALYTIC (the small-c edges s,t->0,1 are smooth
    zeros P~c, not singularities), so plain Gauss-Legendre is SPECTRAL there --
    the earlier "RectB disagrees at 6 digits" was under-resolution at deg=4
    (N=24), not a deep obstruction.

Corner fix (this file): Duffy (s=rho*alpha, t=rho*(1-alpha)) THEN rho=sigma^2.
Because J is analytic in sigma=sqrt(rho) at fixed alpha (pure integer +
half-integer powers of rho, no logs -- the k-integrand decays exponentially at
large k for every order in the rho-expansion, so no ln(rho) mechanism fires),
the substituted radial integrand 2 sigma^3 J is ANALYTIC in sigma and GL is
SPECTRAL -- NO subtraction, NO asymptotic coefficients. This supersedes the
one-term rho^{3/2}-subtraction of routeC_corner_v2.py (whose remainder ~rho^{5/2}
is still non-analytic and only algebraically convergent).

Tiling of [0,1]^2 (delta small):
  T1    = {s+t <= delta}              -- corner triangle, sigma^2 Duffy [spectral]
  T2far = {s,t in[0,delta], s+t>delta}-- far triangle, reflected Duffy  [smooth]
  RectA = [0,delta] x [delta,1]                                         [smooth]
  RectB = [delta,1] x [0,1]                                            [smooth]

All smooth pieces: plain GL, degree pushed to plateau. Independent
cross-validation of the smooth pieces via tanh-sinh (mp.quad, a genuinely
different quadrature family) is in routeC_hp_composite_xcheck via --xcheck.

Usage: python routeC_hp_composite.py [dps] [delta] [deg] [kpdeg] [--conv]
  --conv : run a degree-refinement convergence study instead of one shot.
"""
from __future__ import annotations
import sys
import time
import mpmath as mp

sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_probe9 import J_fixed
from routeC_probe8 import k_grid_sinh_paneled
from routeC_hp_evaluator import std_gl_nodes


# --------------------------------------------------------------------------
# corner triangle T1 via Duffy + rho = sigma^2  (spectral)
# --------------------------------------------------------------------------
def T1_sigma2(delta, deg_sigma, deg_alpha, knodes):
    """int_{s+t<=delta} J ds dt.  Duffy s=rho*alpha, t=rho*(1-alpha),
    ds dt = rho drho dalpha; then rho=sigma^2, drho=2 sigma dsigma, giving
    integrand 2 sigma^3 J on (sigma,alpha) in [0,sqrt(delta)]x[0,1].
    alpha-integral folded to [0,1/2] via s<->t symmetry (x2)."""
    std_s = std_gl_nodes(deg_sigma, mp.mp.prec)
    std_a = std_gl_nodes(deg_alpha, mp.mp.prec)
    sig_max = mp.sqrt(delta)
    half_s = sig_max / 2
    tot = mp.mpf(0)
    for xs, ws in std_s:
        sigma = half_s * (xs + 1)
        rho = sigma * sigma
        pref = (half_s * ws) * 2 * sigma ** 3
        acc = mp.mpf(0)
        for xa, wa in std_a:
            alpha = (xa + 1) / 4          # [-1,1] -> [0,1/2]
            s = rho * alpha
            t = rho * (1 - alpha)
            acc += (wa / 4) * J_fixed(s, t, knodes)
        tot += pref * 2 * acc             # x2 for alpha in [1/2,1]
    return tot


# --------------------------------------------------------------------------
# far triangle T2far (reflected Duffy) -- smooth
# --------------------------------------------------------------------------
def T2far(delta, deg_rho, deg_alpha, knodes):
    """{s,t in [0,delta], s+t>delta} via reflected Duffy
    s=delta-rho*alpha, t=delta-rho*(1-alpha), rho in[0,delta], ds dt=rho drho dalpha."""
    std_r = std_gl_nodes(deg_rho, mp.mp.prec)
    std_a = std_gl_nodes(deg_alpha, mp.mp.prec)
    half_r = delta / 2
    tot = mp.mpf(0)
    for xr, wr in std_r:
        rho = half_r * (xr + 1)
        pref = (half_r * wr) * rho
        acc = mp.mpf(0)
        for xa, wa in std_a:
            alpha = (xa + 1) / 2
            s = delta - rho * alpha
            t = delta - rho * (1 - alpha)
            acc += (wa / 2) * J_fixed(s, t, knodes)
        tot += pref * acc
    return tot


# --------------------------------------------------------------------------
# rectangles -- smooth, plain GL
# --------------------------------------------------------------------------
def rect(deg_s, deg_t, slo, shi, tlo, thi, knodes):
    std_s = std_gl_nodes(deg_s, mp.mp.prec)
    std_t = std_gl_nodes(deg_t, mp.mp.prec)
    hs, ms = (shi - slo) / 2, (shi + slo) / 2
    ht, mt = (thi - tlo) / 2, (thi + tlo) / 2
    tot = mp.mpf(0)
    for xs, ws in std_s:
        s = ms + hs * xs
        acc = mp.mpf(0)
        for xt, wt in std_t:
            t = mt + ht * xt
            acc += wt * J_fixed(s, t, knodes)
        tot += ws * acc
    return (hs * ht) * tot


def assemble(delta, deg, kn, deg_rect=None, verbose=False):
    dr = deg_rect if deg_rect is not None else deg
    parts = {}
    t0 = time.time(); parts['T1'] = T1_sigma2(delta, deg, deg, kn)
    if verbose: print(f"    T1(sigma2)={mp.nstr(parts['T1'],20)} ({time.time()-t0:.1f}s)", flush=True)
    t0 = time.time(); parts['T2far'] = T2far(delta, deg, deg, kn)
    if verbose: print(f"    T2far   ={mp.nstr(parts['T2far'],20)} ({time.time()-t0:.1f}s)", flush=True)
    t0 = time.time(); parts['RectA'] = rect(dr, dr, mp.mpf(0), delta, delta, mp.mpf(1), kn)
    if verbose: print(f"    RectA   ={mp.nstr(parts['RectA'],20)} ({time.time()-t0:.1f}s)", flush=True)
    t0 = time.time(); parts['RectB'] = rect(dr, dr, delta, mp.mpf(1), mp.mpf(0), mp.mpf(1), kn)
    if verbose: print(f"    RectB   ={mp.nstr(parts['RectB'],20)} ({time.time()-t0:.1f}s)", flush=True)
    outer = parts['T1'] + parts['T2far'] + parts['RectA'] + parts['RectB']
    return (8 / mp.pi) * outer, parts


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    flags = [a for a in sys.argv[1:] if a.startswith('--')]
    dps = int(args[0]) if len(args) > 0 else 50
    delta = mp.mpf(args[1]) if len(args) > 1 else mp.mpf('0.05')
    deg = int(args[2]) if len(args) > 2 else 6
    kpdeg = int(args[3]) if len(args) > 3 else 6
    mp.mp.dps = dps
    kn = k_grid_sinh_paneled(mp.mpf(16), mp.mpf(2), kpdeg)
    print(f"dps={dps} delta={delta} k_M={len(kn)}", flush=True)

    if '--conv' in flags:
        prev = None
        for d in range(3, deg + 1):
            t0 = time.time()
            v, parts = assemble(delta, d, kn)
            diff = "" if prev is None else mp.nstr(abs(v - prev), 4)
            print(f"deg={d}: T2={mp.nstr(v,dps-4)}  diff={diff}  ({time.time()-t0:.1f}s)", flush=True)
            prev = v
    else:
        t0 = time.time()
        v, parts = assemble(delta, deg, kn, verbose=True)
        print(f"T2 = {mp.nstr(v, dps-4)}  ({time.time()-t0:.1f}s total)")


if __name__ == '__main__':
    main()
