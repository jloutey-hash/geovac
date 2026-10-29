"""Record BU, control for bu_compare.right_divide: the same exact division must FAIL when the guessed operator is
perturbed (one coefficient of c_0 changed by 1) and when the divisor is the derived operator's own order-8 truncation."""
import json
import random

import flint

import bu_compare as C
import na_modanalyze as A
import na_modguess as G

for p in C.PRIMES:
    seq = json.load(open(f'data/bu_check_seq_{p}.json'))['gen']
    B = G.Basis([flint.nmod(v, p) for v in seq], p)
    nul, keys, (N, nrow) = G.nullity(B, 8, 31, 30)
    Gs = A.standard_form(keys, [N[i, 0] for i in range(N.nrows())], p)
    d = json.load(open(f'data/bu_angle_v2_{p}_9_30.json'))
    Ds = A.standard_form([tuple(k) for k in d['keys']], [flint.nmod(v, p) for v in d['vecs'][0]], p)
    Gp = list(Gs)
    Gp[0] = Gp[0] + flint.nmod_poly([1], p)
    rnd = random.Random(5)
    Gr = [flint.nmod_poly([rnd.randrange(p) for _ in range(c.degree() + 1)], p) for c in Gs]
    print(p, 'true', C.right_divide(Ds, Gs)[0], 'perturbed', C.right_divide(Ds, Gp)[0],
          'random', C.right_divide(Ds, Gr)[0], 'truncation', C.right_divide(Ds, Ds[:-1])[0], flush=True)
