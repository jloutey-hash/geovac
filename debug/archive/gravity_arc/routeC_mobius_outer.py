"""Mobius-map quadrature for the OUTER (s,t) integral -- the same rational-map
trick already validated for the inner (semi-infinite) k-integral
(debug/routeC_probe2.py, probe3.py: k=K0(1-x)/(1+x)), applied here to a
FINITE interval [lo,hi] instead.

Composing x -> u=K0(1-x)/(1+x) in (0,inf) with u -> s=u/(1+u) in (0,1) gives
a single Mobius transformation of x directly onto (0,1):

    s(x) = K0(1-x) / [(1+K0) - (K0-1)x],   x in (-1,1)
    s(1)=0, s(-1)=1;  K0=1 reduces EXACTLY to the plain linear GL map
    s=(1-x)/2 (verified below) -- a genuinely different node family from
    both sin^2(phi) and plain-linear GL for K0 != 1 (rational vs
    trigonometric substitution; different endpoint clustering law).

Rescaled to a general interval via S(x) = lo + (hi-lo)*s(x), clustering
nodes near 'lo' for K0 far from 1 (either direction).
"""
from __future__ import annotations
import sys
import mpmath as mp

sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_hp_evaluator import std_gl_nodes  # noqa: E402


def mobius_interval_nodes(lo, hi, K0, degree):
    """Mobius-mapped GL nodes/weights on [lo,hi], clustering near lo when
    K0 != 1 (K0=1 reduces to plain linear GL, used as the correctness
    check)."""
    std_nodes = std_gl_nodes(degree, mp.mp.prec)
    span = hi - lo
    K0 = mp.mpf(K0)
    a = 1 + K0
    b = K0 - 1
    pts = []
    for x, w in std_nodes:
        denom = a - b * x
        s0 = K0 * (1 - x) / denom
        s = lo + span * s0
        jac = span * 2 * K0 / (denom * denom)
        pts.append((s, w * jac))
    return pts


def mobius_rect_pts(deg_s, deg_t, slo, shi, tlo, thi, K0s, K0t):
    ps = mobius_interval_nodes(slo, shi, K0s, deg_s)
    pt = mobius_interval_nodes(tlo, thi, K0t, deg_t)
    pts = []
    for s, ws in ps:
        for t, wt in pt:
            pts.append((s, t, ws * wt))
    return pts


if __name__ == '__main__':
    # Sanity check: K0=1 must reduce EXACTLY to plain linear GL.
    mp.mp.dps = 30
    nodes_k1 = mobius_interval_nodes(mp.mpf(0), mp.mpf(1), mp.mpf(1), 4)
    from routeC_hp_evaluator import outer_grid_plain
    s_plain, w_plain = outer_grid_plain(4)
    ok = all(abs(s1 - s2) < mp.mpf('1e-25') and abs(w1 - w2) < mp.mpf('1e-25')
             for (s1, w1), s2, w2 in zip(nodes_k1, s_plain, w_plain))
    print("K0=1 reduces to plain linear GL:", ok)
    # show clustering behavior at K0=5 and K0=0.2
    for K0 in [5, 0.2]:
        pts = mobius_interval_nodes(mp.mpf(0), mp.mpf(1), mp.mpf(K0), 3)
        smin = min(p[0] for p in pts)
        print(f"K0={K0}: s_min={mp.nstr(smin, 6)}  (n={len(pts)})")
