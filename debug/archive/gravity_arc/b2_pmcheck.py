"""Record B2, PM verification: the derived recurrence sum_k p_k(L) F_{L+k} = 0 applied to the AO3 generic sequence
built here directly from ao_exact (random.Random(2026) weights, as debug/ao3_run.py), L <= LMAX; plus a perturbed control.
    python debug/b2_pmcheck.py <p> <LMAX>
"""
import json
import random
import sys

import flint

import ao_exact as X

p, LMAX = int(sys.argv[1]), int(sys.argv[2])
d = json.load(open(f'data/b2_recurrence_{p}_4200.json'))
assert d['p'] == p
P = [flint.nmod_poly(c, p) for c in d['polys']]
m = len(P) - 1
Qs = X.coefficients('1.7', '2.5', '1.4', '1.1', LMAX, X.Field(p))
keys = sorted(set(k for q in Qs for k in q), key=repr)
rnd = random.Random(2026)
w = {k: flint.nmod(rnd.randrange(1, p), p) for k in keys}
gen = [sum((w[k] * q.get(k, flint.nmod(0, p)) for k in keys), flint.nmod(0, p)) for q in Qs]


def bad(P):
    return sum(1 for L in range(len(gen) - m)
               if sum((P[k](flint.nmod(L, p)) * gen[L + k] for k in range(m + 1)), flint.nmod(0, p)) != 0)


Pc = list(P)
Pc[7] = Pc[7] + flint.nmod_poly([1], p)
print(p, 'constants', len(keys), 'rows', len(gen) - m, 'nonzero', bad(P), 'control nonzero', bad(Pc), flush=True)
