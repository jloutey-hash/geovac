"""Record NA, analysis of the found operator at every singular point, both primes (reported in the disposition).
For each irreducible factor f of the leading coefficient: degree, multiplicity, the full indicial polynomial's roots
(rationally reconstructed when f is linear; for higher-degree f the integer exponents in [0, 60] from
na_modanalyze.all_points).    python debug/na_points.py <case>   ->  data/na_points_<case>.json
"""
import json
import math
import random
import sys

import flint

import na_exact as E
import na_modanalyze as A
import na_modguess as G
from na_run import CASES, PRIMES


def ratrec(a, p):
    r0, r1 = p, int(a) % p
    s0, s1 = 0, 1
    B = math.isqrt(p // 2)
    while r1 > B:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    return f"{r1 * (1 if s1 > 0 else -1)}/{abs(s1)}"


case = sys.argv[1]
D1, D2, zX, zH = CASES[case]
out = {}
for p in PRIMES:
    Qs = E.coefficients(D1, D2, zX, zH, 1000, E.Field(p))
    keys = sorted(set(k for q in Qs for k in q), key=str)
    rnd = random.Random(2026)
    w = {k: flint.nmod(rnd.randrange(1, p), p) for k in keys}
    gen = [sum((w[k] * q.get(k, flint.nmod(0, p)) for k in keys), flint.nmod(0, p)) for q in Qs]
    B = G.Basis(gen, p)
    nul, ks, (N, nrow) = G.nullity(B, 8, 31, 30)
    c = A.standard_form(ks, [N[i, 0] for i in range(N.nrows())], p)
    pts = []
    for rec in A.all_points(c, p):
        f = rec['factor']
        d = dict(degree=f.degree(), multiplicity=rec['multiplicity'], regular=rec['regular'])
        if f.degree() == 1:
            x0 = -f.coeffs()[0] / f.coeffs()[1]
            d['x'] = ratrec(x0, p)
            loc = A.local(c, x0, p, {})
            fac = loc['indicial'].factor()[1] if loc['regular'] else []
            d['exponents'] = [(ratrec(-g.coeffs()[0] / g.coeffs()[1], p) if g.degree() == 1 else f'deg{g.degree()}', m)
                              for g, m in fac]
        else:
            d['integer_exponents'] = rec.get('integer_exponents')
        pts.append(d)
        print(case, p, d, flush=True)
    out[str(p)] = dict(nullity=nul, points=pts, infinity_regular=A.fuchs_at_infinity(c)[0],
                       degrees=[ck.degree() for ck in c])
json.dump(out, open(f'data/na_points_{case}.json', 'w'), indent=0)
