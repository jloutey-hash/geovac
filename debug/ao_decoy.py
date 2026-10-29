"""AO decoy: Legendre coefficients of a non-D-finite function with the same branch point,
f = Gamma(x + 2) (xs - x)^(7/2); the guesser must show no cliff within budget.  python debug/ao_decoy.py <xs> <Lmax> <dps> <out>"""
import json
import sys

import mpmath as mp
from mpmath.calculus.quadrature import GaussLegendre

xs_s, Lmax, dps, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
mp.mp.dps = dps
xs = mp.mpf(xs_s)
f = lambda x: mp.gamma(x + 2) * (xs - x) ** mp.mpf(3.5)
nodes = GaussLegendre(mp.mp).calc_nodes(9, mp.mp.prec + 20)
F = [mp.mpf(0)] * (Lmax + 1)
for x, w in nodes:
    fx = w * f(x)
    p0, p1 = mp.mpf(1), x
    F[0] += fx * p0
    F[1] += fx * p1
    for l in range(1, Lmax):
        p0, p1 = p1, ((2 * l + 1) * x * p1 - l * p0) / (l + 1)
        F[l + 1] += fx * p1
F = [(2 * l + 1) / mp.mpf(2) * v for l, v in enumerate(F)]
json.dump(dict(xs=xs_s, F=[mp.nstr(v, dps - 5) for v in F]), open(out, 'w'), indent=0)
print('done')
