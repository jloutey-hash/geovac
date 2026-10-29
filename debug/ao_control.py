"""AO guesser controls (record docs/predictions/2026-10-28_angle_operator.md).
  pos   : Legendre coefficients of a known D-finite function with an algebraic branch point outside [-1, 1]
          (the same singular exponent 7/2 as the Landau term), f = e^x (xs - x)^(7/2) + x e^(-2x) (xs - x)^(9/2);
          the guesser must find a recurrence that validates on held-out terms.
  decoy : the same coefficients with relative noise 10^-noise; the guesser must find none at the scored tolerance.
    python debug/ao_control.py <xs> <Lmax> <dps> <out.json>
"""
import json
import sys

import mpmath as mp
from mpmath.calculus.quadrature import GaussLegendre

xs_s, Lmax, dps, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
mp.mp.dps = dps
xs = mp.mpf(xs_s)
f = lambda x: mp.exp(x) * (xs - x) ** mp.mpf(3.5) + x * mp.exp(-2 * x) * (xs - x) ** mp.mpf(4.5)
nodes = GaussLegendre(mp.mp).calc_nodes(9, mp.mp.prec + 20)          # 768 points on [-1, 1]
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
print('done', mp.nstr(F[Lmax], 5))
