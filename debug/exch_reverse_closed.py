"""exch_reverse_closed.py -- Question 1/2 for one case at one pair (a, b), a > b.

(i) fresh exact solve over Q at (a, b) (decaying solution, every solve re-checked, free unknowns counted);
(ii) its coordinates on the nine functions vs the frozen Q(t) closed form specialised at t = b/a (exact);
(iii) the closed form's value with the real branch (LC = ln((a+b)/|b-a|), PV Ei) at the given R, dps 60.
Run: python debug/exch_reverse_closed.py CASE A B R1,R2,...  -> debug/data/exch_reverse_closed_<CASE>_<A>_<B>.json"""
import json
import sys
import time
from fractions import Fraction

import mpmath as mp

sys.path.insert(0, "debug")
import exch_reverse_core as C  # noqa: E402

case, A, B, Rs = sys.argv[1], Fraction(sys.argv[2]), Fraction(sys.argv[3]), sys.argv[4].split(",")
T0 = time.time()
K, co, info = C.fresh_K_Q(case, A, B)
print(case, A, B, info, "%.1fs" % (time.time() - T0), flush=True)
st = C.store_coords_Q(case, A, B)
same = co is not None and all({p: c for p, c in d.items() if c != 0} == {p: c for p, c in e.items() if c != 0}
                              for d, e in zip(co, st))
present = [i for i, d in enumerate(co or []) if d]
print("fresh == stored Q(t) at t = b/a (rescaled):", same, "functions present:", present, flush=True)
mp.mp.dps = 60
vals = {}
for R in Rs:
    v1 = C.value(A, B, K, Fraction(R))
    v2 = C.value(A, B, C.coords_to_K(A, B, st), Fraction(R))
    vals[R] = {"fresh": mp.nstr(v1, 55), "stored": mp.nstr(v2, 55), "diff": mp.nstr(abs(v1 - v2), 3)}
    print(R, vals[R], flush=True)
out = {"case": case, "pair": [str(A), str(B)], "t": str(B / A), "fresh_info": info, "fresh_equals_stored": same,
       "nine_present": present, "laurent": [[min(d), max(d)] if d else None for d in co],
       "values_dps60": vals, "time_s": round(time.time() - T0, 1)}
tg = f"{case}_{A}_{B}".replace("/", "-")
json.dump(out, open(f"debug/data/exch_reverse_closed_{tg}.json", "w"), indent=1)
