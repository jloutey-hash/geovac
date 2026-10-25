"""Probe: reduce the kite and its z-derivative at one point; count masters; time it."""
import sys, time
from fractions import Fraction as F
from ib_reduce import Reducer, shifts_box, dz_terms, fr
from ib_kite import kite_family, TARGET

R_, S_ = int(sys.argv[1]), int(sys.argv[2])
fam = kite_family()
z = F(7, 3)
fam.setup({("w", "w"): z}, 3)
seeds = shifts_box(5, R_, S_, fam.integer)
seeds = [n for n in seeds if not fam.is_zero(n)]
t0 = time.time()
R = Reducer(fam)
R.seed(seeds, sym=(lambda n: (n[1], n[0], n[3], n[2], n[4])) if len(sys.argv) > 3 else None)
t1 = time.time()
red0 = R.reduce({TARGET: 1})
red1 = R.reduce(dz_terms(fam, TARGET, "w", fr(z)))
ms = sorted(set(red0) | set(red1), key=R.key)
print(f"seeds {len(seeds)} rows {len(R.rows)} time {t1-t0:.1f}s")
print("integrals in target + derivative:", len(ms))
for n in ms:
    print("  ", n, R.key(n)[:4])
