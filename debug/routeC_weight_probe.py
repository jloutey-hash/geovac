"""Route C weight probe: is the two-body 3-centre ERI weight-1 or weight-2?

SUPERSEDED (2026-08-16, same day): the transcendence CLASS was identified in
debug/routeC_momentum_poc.py -- it is ELLIPTIC (genus 1), NOT a polylogarithm, so the
dilog/zeta(2) PSLQ bases below are the wrong class and never land.  What survives from
this script and is still cited: (i) the high-precision collinear value 0.395355766...,
(ii) the weight-0 exclusion (the value is provably NOT a rational combination of
{1, e^-2, e^-4}, consistent with a genus-1 object).  Kept as the record of the
weight-vs-genus turn; do not re-run the PSLQ legs expecting a hit.

Original intent:

High-precision evaluation of the collinear symmetric case (X=0, Y=(0,0,1),
Z=(0,0,-1), all 1s zeta=1, so D1=D2=1 and |W(s,t)|=s+t) via the exact reduced
form (analytic d/dzeta, no finite difference):

    (XY|XZ) = (8/pi) int_0^1 ds int_0^1 dt int_0^inf dk
                j0(k(s+t)) P1(s,k) P2(t,k)
    Pi(x,k) = c e^{-D Delta}(D^2/Delta^3 + 3D/Delta^4 + 3/Delta^5),
              c = x(1-x), Delta = sqrt(c k^2 + 1), D = 1.

Then PSLQ the number against escalating bases:
  weight 0 : {1, e^-2, e^-4}                     (exp-polynomial, as the 2-centre
                                                   exchange J(R) turned out to be)
  weight 1 : + {gamma, ln2, E1(2), e^-2 * those} (the 2-centre exchange seed set)
  weight 2 : + {zeta(2)=pi^2/6, Li2(e^-2), Catalan}

A weight-<=1 relation => Route C closes at weight one (positive, matches the arc).
No weight-1 relation but a weight-2 one => three centres introduce weight two
(the structurally-informative negative).
"""

from __future__ import annotations

import sys
from pathlib import Path

import mpmath as mp

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

mp.mp.dps = 26


def P(x, k, D):
    c = x * (1 - x)
    Del = mp.sqrt(c * k * k + 1)
    return c * mp.e ** (-D * Del) * (D * D / Del ** 3 + 3 * D / Del ** 4 + 3 / Del ** 5)


def I_rad(s, t, D1=mp.mpf(1), D2=mp.mpf(1)):
    b = s + t
    # exponential envelope e^{-D sqrt(c k^2+1)} dominates the j0 oscillation;
    # plain tanh-sinh with breakpoints matches quadosc to full precision, ~10x faster.
    def f(k):
        j0 = mp.sin(k * b) / (k * b) if k * b > mp.mpf('1e-20') else mp.mpf(1)
        return j0 * P(s, k, D1) * P(t, k, D2)
    return mp.quad(f, [0, 1, 3, 8, 20, mp.inf])


def collinear_eri():
    """(s,t) integral by nested tanh-sinh (mp.quad) -- handles the s*ln(s) / sqrt(s)
    endpoint non-analyticity (from c=s(1-s)->0) that Gauss-Legendre resolves poorly.
    Exploits the s<->t symmetry of the collinear D1=D2 case is NOT done (kept direct
    for a clean error estimate)."""
    inner = lambda s: mp.quad(lambda t: I_rad(s, t), [0, mp.mpf('0.5'), 1])
    return (8 / mp.pi) * mp.quad(inner, [0, mp.mpf('0.5'), 1])


def main():
    val = collinear_eri()
    print("collinear (XY|XZ) via nested tanh-sinh:")
    print("  ", mp.nstr(val, 32))
    # cross-check against the double-precision reduced form (0.39535577)
    print("   (double-precision reference was 0.3953557770)")

    e1 = mp.e ** -1
    e2 = mp.e ** -2
    e4 = mp.e ** -4
    g = mp.euler
    ln2 = mp.log(2)
    E12 = mp.e1(2)
    E14 = mp.e1(4)
    z2 = mp.pi ** 2 / 6
    Li2_2 = mp.polylog(2, e2)
    Li2_4 = mp.polylog(2, e4)
    Li2_m2 = mp.polylog(2, -e2)
    Cat = mp.catalan

    def try_pslq(name, basis, tol=mp.mpf(10) ** -20):
        vec = [val] + basis
        rel = mp.pslq(vec, tol=tol, maxcoeff=10 ** 9, maxsteps=10 ** 5)
        print(f"  {name:10s}: {rel}")
        return rel

    print("PSLQ (first coeff multiplies the ERI value; small coeffs => genuine):")
    try_pslq("weight0", [mp.mpf(1), e2, e4])
    try_pslq("weight1", [mp.mpf(1), e2, e4, g, ln2, E12, E14])
    try_pslq("weight1b", [mp.mpf(1), e2, e4, g * e2, ln2 * e2, E12, E14])
    try_pslq("weight2", [mp.mpf(1), e2, e4, z2, Li2_2, Li2_4, z2 * e2, Li2_2 * e2])
    try_pslq("weight2b", [mp.mpf(1), e2, e4, g * e2, ln2 * e2, z2, Li2_2, Li2_m2, Cat])


if __name__ == "__main__":
    main()
