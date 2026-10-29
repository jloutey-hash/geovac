"""Record B2: the characteristic polynomial's irreducible quadratic factors (not handled by b2_asym), with
coefficients rationally reconstructed, and the x they map to (x = (lam + 1/lam)/2 solves a quadratic too).
    python b2_charpoly.py <recurrence.json>"""
import json
import sys

import flint

from b2_asym import ratrec

d = json.load(open(sys.argv[1]))
p = d['p']
P = [flint.nmod_poly(c, p) for c in d['polys']]
deg = max(q.degree() for q in P)
chi = flint.nmod_poly([int(q.coeffs()[deg]) if q.degree() == deg else 0 for q in P], p)
for f, mult in chi.factor()[1]:
    if f.degree() == 2:
        c0, c1, c2 = f.coeffs()
        a, b = ratrec(c1 / c2, p), ratrec(c0 / c2, p)
        # lam^2 + a lam + b = 0; lam*lam' = b, lam + lam' = -a
        print('lam^2 + (%s) lam + (%s), mult %d' % (a, b, mult))
