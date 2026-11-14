"""Positive control for record AB (docs/predictions/2026-11-14_aa_bc_three_centre.md): the NA engine, unchanged,
rebuilds the order-8, degree-31 operator of the one-electron integral at N1 (= the Coulomb piece of case A1).
Same weights and search as debug/na_run.py scan; orders 7 and 8 only.  Output debug/data/aabc_na_posctrl.json."""
import json
import random
import time

import flint

import na_exact as E
import na_modguess as G

out = {}
for p in [(1 << 61) - 1, 4611686018427387847]:
    t0 = time.time()
    Qs = E.coefficients('1.8', '2.6', '1.3', '1.0', 1000, E.Field(p))
    keys = sorted(set(k for q in Qs for k in q), key=str)
    rnd = random.Random(2026)
    w = {k: flint.nmod(rnd.randrange(1, p), p) for k in keys}
    gen = [sum((w[k] * q.get(k, flint.nmod(0, p)) for k in keys), flint.nmod(0, p)) for q in Qs]
    B = G.Basis(gen, p)
    r = {R: G.min_degree_bisect(B, R, 30, len(gen)) for R in (7, 8)}
    out[str(p)] = dict(constants=[str(k) for k in keys], trade={R: list(v) for R, v in r.items()},
                       seconds=round(time.time() - t0))
    print(p, out[str(p)], flush=True)
json.dump(out, open('data/aabc_na_posctrl.json', 'w'), indent=1)
