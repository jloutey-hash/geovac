"""Build tests/data/paper59_intclass_2026_11_02.json, the stored copy the Paper 59 IC guard reads
(tests/test_paper59_integrated_class.py).  Recomputes the exact split of W1 (record IC, case W1 =
AO3's scored case) mod the second scoring prime to L = 600 with debug/ao_exact.py (unchanged since
AO3), and keeps the component sequences of the two ln constants, whose stuck residuals have only
degree-one singular factors (record IC, Caps erratum).  Run once; the guard never imports from debug/.
"""
import hashlib
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import ao_exact as X

P = 1073741833
L = 600
CASE = ('1.7', '2.5', '1.4', '1.1')            # W1: D1, D2, zX, zH
KEEP = ("(Fraction(-231, 50), (('ln', Fraction(5, 2)),))", "(Fraction(-231, 50), (('ln', Fraction(3, 10)),))")
OUT = HERE.parent / 'tests' / 'data' / 'paper59_intclass_2026_11_02.json'

t0 = time.time()
Qs = X.coefficients(*CASE, L, X.Field(P))
keys = sorted(set(k for q in Qs for k in q), key=repr)
comp = {}
for k in keys:
    if repr(k) in KEEP:
        comp[repr(k)] = [int(q[k]) if k in q else 0 for q in Qs]
assert sorted(comp) == sorted(KEEP), sorted(comp)
src = {f: hashlib.sha256((HERE / f).read_bytes()).hexdigest() for f in ('ao_exact.py', 'na_exact.py')}
OUT.write_text(json.dumps({'prime': P, 'L': L, 'case': {'D1': CASE[0], 'D2': CASE[1], 'zX': CASE[2], 'zH': CASE[3]},
                           'nconst': len(keys), 'source_sha256': src, 'components': comp}, indent=0))
print('wrote', OUT, len(keys), 'constants', round(time.time() - t0, 1), 's')
