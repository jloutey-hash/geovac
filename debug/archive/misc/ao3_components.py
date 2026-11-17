"""Record AO3, reported only: the smallest operator of each of the 31 constant components of the two-electron
bond-angle coefficients, L = 0..1000, one prime.   python debug/ao3_components.py <case>  ->  data/ao3_comp_<case>.json
"""
import json
import sys
import time

import flint

import ao_exact as X
import na_modguess as G
CASES = {'W1': ('1.7', '2.5', '1.4', '1.1'), 'W2': ('1.6', '2.2', '1.0', '1.2')}  # as debug/ao3_run.py

case = sys.argv[1]
p = (1 << 61) - 1
D1, D2, zX, zH = CASES[case]
t0 = time.time()
Qs = X.coefficients(D1, D2, zX, zH, 1000, X.Field(p))
keys = sorted(set(k for q in Qs for k in q), key=repr)
print(case, 'coefficients', round(time.time() - t0), flush=True)
out = {}
for k in keys:
    seq = [q.get(k, flint.nmod(0, p)) for q in Qs]
    B = G.Basis(seq, p)
    best = None
    for R in range(1, 25):
        Q0, info = G.min_degree_bisect(B, R, 30, len(seq))
        if Q0 is not None:
            best = (R, Q0, info)
            break
    out[repr(k)] = best
    print(case, repr(k), best, round(time.time() - t0), flush=True)
json.dump(out, open(f'data/ao3_comp_{case}.json', 'w'), indent=0)
