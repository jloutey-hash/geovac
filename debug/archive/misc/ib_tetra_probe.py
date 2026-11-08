import sys, time
from fractions import Fraction as F
from ib_reduce import Reducer, shifts_box, fr
from ib_tetra import tetra_family, TARGET, dx_terms
R_, S_ = int(sys.argv[1]), int(sys.argv[2])
fam = tetra_family()
x = 987654321
m2 = [F(1), F(1), x, F(9, 4), F(9, 4), F(0)]
fam.setup({}, 3, masses2=[fr(v) if isinstance(v, F) else v for v in m2])
seeds = [n for n in shifts_box(6, R_, S_, fam.integer) if not fam.is_zero(n)]
print("seeds", len(seeds), flush=True)
t0 = time.time()
Rd = Reducer(fam)
Rd.seed(seeds, sym=lambda n: (n[1], n[0], n[2], n[4], n[3], n[5]))
print(f"seeded {time.time()-t0:.0f}s rows {len(Rd.rows)}", flush=True)
c = Rd.reduce({TARGET: 1})
ms, seen, i = list(c), set(c), 0
while i < len(ms) and len(ms) < 60:
    r = Rd.reduce(dx_terms(fam, ms[i])); i += 1
    for n in r:
        if n not in seen:
            seen.add(n); ms.append(n)
print("closure", len(ms), "closed" if i == len(ms) else "NOT closed", f"{time.time()-t0:.0f}s")
