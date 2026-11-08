"""Exact local exponents over Q of an integer d/dz-form operator at rational points, with a logarithm check
for repeated or integer-spaced exponents left to the caller.  usage: python ib_exponents_Q.py <op json> z0 ..."""
import json, sys
from fractions import Fraction as F
from flint import fmpq_poly, fmpq, fmpz_poly
import sympy as sp

d = json.load(open(sys.argv[1]))
ops = [fmpq_poly([int(v) for v in p]) for p in d["polys"]]
order = len(ops) - 1
r = sp.Symbol("r")
for zs in sys.argv[2:]:
    z0 = F(zs)
    sh = fmpq_poly([fmpq(z0.numerator, z0.denominator), 1])
    Ps = [p(sh) for p in ops]
    v = [min((i for i, c in enumerate(p.coeffs()) if c != 0), default=10**9) for p in Ps]
    delta = min(v[k] - k for k in range(order + 1))
    fuchs = delta == v[-1] - order
    ind = 0
    for k in range(order + 1):
        if v[k] - k == delta:
            c = Ps[k].coeffs()[v[k]]
            ind += sp.Rational(int(c.p), int(c.q)) * sp.ff(r, k)
    roots = sp.roots(sp.Poly(sp.expand(ind), r))
    tot = sum(roots.values())
    lead_vanish = Ps[-1].coeffs()[0] == 0
    print(f"z0 = {zs}: {'regular singular' if fuchs else 'IRREGULAR'}; leading coefficient "
          f"{'vanishes' if lead_vanish else 'nonzero'}; exponents {dict(sorted(roots.items()))}"
          + ("" if tot == order else f" (found {tot} of {order})"))
