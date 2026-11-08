"""Trace the master-basis closure under d/dz at one z, with a hard bound."""
import sys, time
from fractions import Fraction as F
from ib_reduce import Reducer, shifts_box, dz_terms, fr
from ib_kite import kite_family, TARGET
R_, S_ = int(sys.argv[1]), int(sys.argv[2])
fam = kite_family()
z = 987654321
fam.setup({("w", "w"): z}, 3)
seeds = [n for n in shifts_box(5, R_, S_, fam.integer) if not fam.is_zero(n)]
Rd = Reducer(fam); t0 = time.time()
Rd.seed(seeds, sym=lambda n: (n[1], n[0], n[3], n[2], n[4]))
print("seeded", round(time.time() - t0), "s", flush=True)
c = Rd.reduce({TARGET: 1})
masters = list(c); seen = set(masters); i = 0
while i < len(masters) and len(masters) < 40:
    m = masters[i]; i += 1
    r = Rd.reduce(dz_terms(fam, m, "w", z))
    new = [n for n in r if n not in seen]
    print(m, Rd.key(m)[:4], "->", len(r), "terms, new:", [(n, Rd.key(n)[:4], n in Rd.rows) for n in new], flush=True)
    for n in new:
        seen.add(n); masters.append(n)
print("closure size", len(masters), "closed" if i == len(masters) else "NOT closed")
