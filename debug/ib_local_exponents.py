"""Local exponents of a d/dz-form operator mod P at rational points (indicial polynomial; rational roots
found by testing r = j/den for |j| <= 60*den, den in 1, 2, 4).  usage: python ib_local_exponents.py <op json> z0 ..."""
import json, sys
from fractions import Fraction as F
from flint import nmod_poly
from ib_reduce import fr

d = json.load(open(sys.argv[1])); P = d["P"]
ops = [nmod_poly(p, P) for p in d["polys"]]
order = len(ops) - 1
for zs in sys.argv[2:]:
    z0 = F(zs)
    sh = nmod_poly([fr(z0), 1], P)
    Ps = [p(sh) for p in ops]           # p_k(z0 + s)
    v = [min((i for i, c in enumerate(p.coeffs()) if int(c)), default=10**9) for p in Ps]
    delta = min(v[k] - k for k in range(order + 1))
    fuchs = (delta == v[-1] - order)
    # indicial polynomial in r
    ind = nmod_poly([0], P)
    for k in range(order + 1):
        if v[k] - k == delta:
            ff = nmod_poly([1], P)
            for j in range(k):
                ff = ff * nmod_poly([(-j) % P, 1], P)
            ind = ind + ff * int(Ps[k].coeffs()[v[k]])
    roots = []
    rem = ind
    for den in (1, 2, 4):
        for j in range(-60 * den, 60 * den + 1):
            r = F(j, den)
            if r.denominator != den:
                continue
            lin = nmod_poly([(-fr(r)) % P, 1], P)
            while rem.degree() > 0 and rem(fr(r)) == 0:
                rem = rem // lin
                roots.append(r)
    print(f"z0 = {zs}: {'regular singular' if fuchs else 'IRREGULAR'}, indicial degree {ind.degree()}, "
          f"rational roots {sorted(roots)}" + (f", unfound degree {rem.degree()}" if rem.degree() > 0 else ""))
