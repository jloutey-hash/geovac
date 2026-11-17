"""Record B2, CHECK DATA ONLY: the exact W1 generic sequence (debug/ao_exact.py) mod both primes, AO3's weights
(keys sorted by repr, random.Random(2026).randrange(1, p) in that order, as debug/ao3_run.py scan), L = 0..LMAX,
plus each constant component.  -> data/b2_check_seq_<p>.json.  Nothing here feeds the derivation."""
import json
import random
import sys
import time

import flint

import ao_exact as X

PRIMES = [(1 << 61) - 1, 4611686018427387847]
LMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 300
for p in PRIMES:
    t0 = time.time()
    Qs = X.coefficients('1.7', '2.5', '1.4', '1.1', LMAX, X.Field(p))
    keys = sorted(set(k for q in Qs for k in q), key=repr)
    rnd = random.Random(2026)
    w = {k: flint.nmod(rnd.randrange(1, p), p) for k in keys}
    gen = [sum((w[k] * q.get(k, flint.nmod(0, p)) for k in keys), flint.nmod(0, p)) for q in Qs]
    comps = {repr(k): [int(q.get(k, flint.nmod(0, p))) for q in Qs] for k in keys}
    json.dump(dict(p=p, gen=[int(v) for v in gen], components=comps), open(f'data/b2_check_seq_{p}.json', 'w'))
    print('done', p, len(keys), 'constants', round(time.time() - t0), 's', flush=True)
