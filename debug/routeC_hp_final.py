"""FINAL high-precision composite evaluator for T2 (Paper 59), targeting >=40
cross-validated digits.

Resolved facts feeding this (see sprint_hp_evaluator_memo addendum + this
session's routeC_rectb_resolve/arbiter + routeC_corner_sigma2):
  * The (0,0) corner is a rho^{3/2} non-analyticity (NOT a log). Handled by
    Duffy + rho=sigma^2, which is SPECTRAL (routeC_corner_sigma2.py: diffs
    1e-13->8e-16->8e-20->3e-25).
  * Everywhere else J(s,t) is analytic (small-c edges s,t->0,1 are smooth
    zeros P~c). GL converges there; the old "anchor vs composite digit-9
    discrepancy" was pure RectB GL under-resolution at deg=4 (GL deg=6 ==
    tanh-sinh to 6e-17), NOT an anchor bias. Anchor 0.3953557659017139
    CONFIRMED at ~12 digits.
  * k-grid pdeg7 (M=1536) gives J to 40-50 digits across the domain (only the
    doubly-weight-suppressed extreme edge s->1,t->0 drops to ~20 digits).

Symmetric tiling (cleaner than the wide RectB=[delta,1]x[0,1] which carried
both t->0 and t->1 edges):
    [0,1]^2 = T1 {s+t<=delta}              (corner triangle, sigma^2 Duffy)
            + T2far {[0,delta]^2, s+t>delta}(reflected Duffy)
            + 2 * RectA [0,delta]x[delta,1] (strip; x2 for the s<->t reflection)
            + BigSquare [delta,1]^2         (compact, only s=1,t=1 edges)
Smooth pieces use GL with a cluster-both-ends node map x=(1-cos(pi u))/2 (a
"sin^2" clustering matching the c=x(1-x) edge zeros) to accelerate convergence.

Usage: python routeC_hp_final.py [dps] [delta] [deg] [kpdeg] [--conv] [--parts]
"""
from __future__ import annotations
import sys
import time
import mpmath as mp

sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_probe9 import J_fixed
from routeC_probe8 import k_grid_sinh_paneled
from routeC_hp_evaluator import std_gl_nodes
from routeC_hp_composite import T1_sigma2, T2far


def clustered_nodes(deg, lo, hi):
    """GL in u in [0,1] pushed through x=(1-cos(pi u))/2 (clusters at BOTH
    ends), then affine to [lo,hi]. Returns (x, weight) with the dx/du and
    interval Jacobians folded in."""
    std = std_gl_nodes(deg, mp.mp.prec)
    span = hi - lo
    out = []
    for xu, wu in std:
        u = (xu + 1) / 2                 # [-1,1] -> [0,1]
        wU = wu / 2
        x = (1 - mp.cos(mp.pi * u)) / 2  # cluster both ends
        dxdu = (mp.pi / 2) * mp.sin(mp.pi * u)
        out.append((lo + span * x, wU * dxdu * span))
    return out


def clustered_rect(deg_s, deg_t, slo, shi, tlo, thi, kn):
    ns = clustered_nodes(deg_s, slo, shi)
    nt = clustered_nodes(deg_t, tlo, thi)
    tot = mp.mpf(0)
    for s, ws in ns:
        acc = mp.mpf(0)
        for t, wt in nt:
            acc += wt * J_fixed(s, t, kn)
        tot += ws * acc
    return tot


def assemble(delta, deg, kn, verbose=False):
    parts = {}
    t0 = time.time(); parts['T1'] = T1_sigma2(delta, deg, deg, kn)
    if verbose: print(f"    T1        ={mp.nstr(parts['T1'],22)} ({time.time()-t0:.1f}s)", flush=True)
    t0 = time.time(); parts['T2far'] = T2far(delta, deg, deg, kn)
    if verbose: print(f"    T2far     ={mp.nstr(parts['T2far'],22)} ({time.time()-t0:.1f}s)", flush=True)
    t0 = time.time(); parts['RectA'] = clustered_rect(deg, deg, mp.mpf(0), delta, delta, mp.mpf(1), kn)
    if verbose: print(f"    RectA     ={mp.nstr(parts['RectA'],22)} ({time.time()-t0:.1f}s)", flush=True)
    t0 = time.time(); parts['BigSq'] = clustered_rect(deg, deg, delta, mp.mpf(1), delta, mp.mpf(1), kn)
    if verbose: print(f"    BigSquare ={mp.nstr(parts['BigSq'],22)} ({time.time()-t0:.1f}s)", flush=True)
    outer = parts['T1'] + parts['T2far'] + 2 * parts['RectA'] + parts['BigSq']
    return (8 / mp.pi) * outer, parts


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    flags = [a for a in sys.argv[1:] if a.startswith('--')]
    dps = int(args[0]) if len(args) > 0 else 50
    delta = mp.mpf(args[1]) if len(args) > 1 else mp.mpf('0.05')
    deg = int(args[2]) if len(args) > 2 else 6
    kpdeg = int(args[3]) if len(args) > 3 else 7
    mp.mp.dps = dps
    kn = k_grid_sinh_paneled(mp.mpf(16), mp.mpf(2), kpdeg)
    print(f"dps={dps} delta={delta} k_M={len(kn)} (pdeg={kpdeg})", flush=True)
    anchor = mp.mpf('0.3953557659017139')

    if '--conv' in flags:
        prev = None
        for d in range(3, deg + 1):
            t0 = time.time()
            v, parts = assemble(delta, d, kn, verbose='--parts' in flags)
            diff = "" if prev is None else mp.nstr(abs(v - prev), 4)
            print(f"deg={d}: T2={mp.nstr(v,dps-4)}  diff={diff}  vsAnchor={mp.nstr(abs(v-anchor),4)}  ({time.time()-t0:.1f}s)", flush=True)
            prev = v
    else:
        t0 = time.time()
        v, parts = assemble(delta, deg, kn, verbose=True)
        print(f"T2 = {mp.nstr(v, dps-4)}", flush=True)
        print(f"vs anchor 0.3953557659017139 : {mp.nstr(abs(v-anchor),4)}  ({time.time()-t0:.1f}s total)")


if __name__ == '__main__':
    main()
