"""Record BU, CHECK DATA ONLY: the generic N1 sequence mod both primes, formed exactly as debug/na_run.py forms it
(random.Random(2026) weights over the sorted constant keys), plus each component, L = 0..1000.
Written to data/bu_check_seq_<p>.json.  Nothing here feeds the derivation (bu_ibp.py, bu_system.py, bu_derive.py)."""
import json
import random
import sys

import flint

import na_exact as E

PRIMES = [(1 << 61) - 1, 4611686018427387847]
for p in PRIMES:
    Qs = E.coefficients('1.8', '2.6', '1.3', '1.0', 1000, E.Field(p))
    keys = sorted(set(k for q in Qs for k in q), key=str)
    rnd = random.Random(2026)
    wts = {k: flint.nmod(rnd.randrange(1, p), p) for k in keys}
    gen = [sum((wts[k] * q.get(k, flint.nmod(0, p)) for k in keys), flint.nmod(0, p)) for q in Qs]
    comps = {str(k): [int(q.get(k, flint.nmod(0, p))) for q in Qs] for k in keys}
    json.dump(dict(p=p, gen=[int(v) for v in gen], components=comps, weights={str(k): int(v) for k, v in wts.items()}),
              open(f'data/bu_check_seq_{p}.json', 'w'))
    print('done', p, flush=True)
