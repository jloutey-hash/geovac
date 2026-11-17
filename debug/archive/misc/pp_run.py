"""Record PP (docs/predictions/2026-10-29_hydrogen_p_functions.md): one exact point.
    python debug/pp_run.py <name> <D1> <D2> <zX> <z1> <z2> [LMAX]      (exact fractions as 'a/b')
writes debug/data/pp/<name>.json with q_L, L = 0..LMAX (default 150), as 600-digit strings (evaluated at 6*LMAX+400 digits).
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp

import ao_exact
import na_exact
import pp_exact

name, D1, D2, zX, z1, z2 = sys.argv[1:7]
LMAX = int(sys.argv[7]) if len(sys.argv) > 7 else 150
t0 = time.time()
Q = pp_exact.coefficients2(D1, D2, zX, z1, z2, LMAX, na_exact.Field(None))
t1 = time.time()
mp.mp.dps = 6 * LMAX + 400
cv, vals = {}, []
for q in Q:
    s = mp.mpf(0)
    for k, v in q.items():
        if k not in cv:
            cv[k] = ao_exact.value(k, mp)
        s += mp.mpf(int(v.p)) / int(v.q) * cv[k]
    vals.append(mp.nstr(s, 600))
out = Path(__file__).parent / 'data' / 'pp'
out.mkdir(parents=True, exist_ok=True)
json.dump(dict(name=name, args=[D1, D2, zX, z1, z2], lmax=LMAX, q=vals,
               build_s=round(t1 - t0), eval_s=round(time.time() - t1)), open(out / f'{name}.json', 'w'))
print(name, 'build', round(t1 - t0), 's eval', round(time.time() - t1), 's', flush=True)
