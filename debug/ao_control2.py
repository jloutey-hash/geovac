"""Full-size controls for the bond-angle operator search (record NA).
  pos K   : f = sum_{k=1..K} e^{c_k x} (xs - x)^(1/2 + k),  c_k = k/2 - 1 (k = 1..K): a known D-finite function of order K
            (LCLM of K first-order operators), same branch structure as the integral (half-integer power at xs).
  decoy   : f = Gamma(x + 2) (xs - x)^(3/2), not D-finite.
Legendre coefficients by Gauss-Legendre on [-1, 1] (3 * 2^(deg-1) nodes).
    python debug/ao_control2.py <pos|decoy> <K> <xs> <Lmax> <dps> <deg> <out.json>
"""
import json
import sys

import mpmath as mp
from mpmath.calculus.quadrature import GaussLegendre

kind, K, xs_s, Lmax, dps, deg, out = sys.argv[1], int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), \
    int(sys.argv[6]), sys.argv[7]
mp.mp.dps = dps
xs = mp.mpf(xs_s)
if kind == 'pos':
    f = lambda x: sum(mp.exp((mp.mpf(k) / 2 - 1) * x) * (xs - x) ** (mp.mpf(1) / 2 + k) for k in range(1, K + 1))
else:
    f = lambda x: mp.gamma(x + 2) * (xs - x) ** mp.mpf(1.5)
nodes = GaussLegendre(mp.mp).calc_nodes(deg, mp.mp.prec + 20)
F = [mp.mpf(0)] * (Lmax + 1)
for x, w in nodes:
    fx = w * f(x)
    p0, p1 = mp.mpf(1), x
    F[0] += fx
    F[1] += fx * p1
    for l in range(1, Lmax):
        p0, p1 = p1, ((2 * l + 1) * x * p1 - l * p0) / (l + 1)
        F[l + 1] += fx * p1
F = [(2 * l + 1) / mp.mpf(2) * v for l, v in enumerate(F)]
json.dump(dict(kind=kind, K=K, xs=xs_s, F=[mp.nstr(v, dps - 5) for v in F]), open(out, 'w'), indent=0)
print('done', len(nodes), mp.nstr(F[Lmax], 5))
