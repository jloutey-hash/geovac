"""Record AB, P3a detail: the minimal operator of the generic sequence (A1/A2), every finite singular point with its
full exponent set, and whether its position moves when zeta_A moves (zeta_A + 1/997, same search).
    python debug/aabc_p3a_points.py <case> <prime> [shift]   -> data/aabc_p3a_points_<case>_<p>[_shift].json
Exponents: at a linear factor, the roots of the indicial polynomial (factored mod p, rational reconstruction); at a
factor f of degree > 1, the integer exponents in [0, 60] (na_modanalyze) plus the remaining one(s) from the sum of
the roots, computed in F_p[t]/f (reported as rational if it reduces to a constant).
"""
import json
import random
import sys
import time
from fractions import Fraction
from pathlib import Path

import flint

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import aabc_exact as X
import intclass_lib as I
import na_exact as E
import na_modanalyze as A
import na_modguess as G

case, p = sys.argv[1], int(sys.argv[2])
shift = len(sys.argv) > 3
D1, D2, zB, zC, zA = X.params(case)
if shift:
    zA = zA + Fraction(1, 997)
t0 = time.time()
K = E.Field(p)
pc = X.pieces(case, 1000, K, params_override=(D1, D2, zB, zC, zA))
f = X.combine(pc, zA, K)
keys = sorted(set(k for q in f for k in q), key=str)
rnd = random.Random(2026)
w = {k: flint.nmod(rnd.randrange(1, p), p) for k in keys}
gen = [sum((w[k] * q.get(k, flint.nmod(0, p)) for k in keys), flint.nmod(0, p)) for q in f]
print('coefficients', round(time.time() - t0), 's', len(keys), 'constants', flush=True)
B = G.Basis(gen, p)
Rm = None
for R in range(1, 25):
    Q0, nul = G.min_degree_bisect(B, R, 30, len(gen))
    if Q0 is not None:
        Rm = R
        break
nul, ks, (N, nrow) = G.nullity(B, Rm, Q0, 30)
c = A.standard_form(ks, [N[i, 0] for i in range(N.nrows())], p)
print('operator', Rm, Q0, nul, flush=True)
xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)


def rr(v):
    fr = I.ratrec(int(v), p, bound=10 ** 9)
    return str(fr) if fr is not None else f'mod-p:{int(v)}'


pts = []
for rec in A.all_points(c, p):
    fac = rec['factor']
    d = dict(degree=rec['degree'], mult=rec['multiplicity'], regular=rec['regular'],
             apparent_type=rec.get('apparent_type'))
    if rec['degree'] == 1:
        x0 = -fac[0]
        d['position'] = rr(x0)
        d['is_xstar'] = int(x0) == int(flint.nmod(xs.numerator, p) / flint.nmod(xs.denominator, p))
        loc = A.local(c, x0, p, {})
        if loc['regular']:
            ip = loc['indicial']
            roots = []
            for g, m in ip.factor()[1]:
                if g.degree() == 1:
                    roots += [rr(-g[0] / g[1])] * m
                else:
                    roots.append(f'deg{g.degree()}-factor x{m}')
            d['exponents'] = roots
    else:
        d['position'] = str(fac)[:80]
        ints = rec.get('integer_exponents') or []
        d['integer_exponents'] = ints
        R = len(c) - 1
        if rec['regular'] and len(ints) == R - 1:
            # sum of exponents = R(R-1)/2 - lead[R-1]/lead[R] in F_p[t]/f
            vals = [A.valuation(ck, fac) for ck in c]
            ords = [v[0] for v in vals]
            delta = ords[R] - R
            fp_ = fac.derivative()
            lead = []
            for k in (R - 1, R):
                if ords[k] == k + delta:
                    lead.append((vals[k][1] % fac) * (fp_ ** ords[k] % fac) % fac)
                else:
                    lead.append(flint.nmod_poly([], p))
            # inverse of lead[R] mod fac
            g, s_, t_ = lead[1].xgcd(fac)
            ratio = (lead[0] * s_ * flint.nmod_poly([int(1 / g[0])], p)) % fac
            tot = flint.nmod_poly([R * (R - 1) // 2], p) - ratio
            rem = (tot - flint.nmod_poly([sum(ints)], p)) % fac
            d['remaining_exponent'] = rr(rem[0]) if rem.degree() <= 0 else 'non-constant in F_p[t]/f'
    pts.append(d)
    print(d, flush=True)
out = dict(case=case, prime=p, zA=str(zA), R=Rm, Q0=Q0, nullity=nul, points=pts, seconds=round(time.time() - t0))
(HERE / 'data' / f"aabc_p3a_points_{case}_{p}{'_shift' if shift else ''}.json").write_text(json.dumps(out, indent=1))
