"""exch_reverse_review_hp.py -- REVIEW (cold) of debug/exch_reverse_memo.md: exchange_hp reference values at NEW
a > b pairs.  Neumann sum over tau (geovac.qfd_core.exchange_hp), never sees the closed form.
Run: python debug/exch_reverse_review_hp.py CASE A B R TAU DPS
 -> debug/data/exch_reverse_review_hp_<CASE>_<A>_<B>_<R>_T<TAU>.json  (full per-tau list kept)"""
import json
import sys
import time

import mpmath as mp
import sympy as sp

sys.path.insert(0, ".")
from geovac import qfd_core as Q  # noqa: E402

LM = {"2psig": ((1, 0), (1, 0)), "3dsig": ((2, 0), (2, 0)), "2ppi": ((1, 1), (1, -1)), "3ddel": ((2, 2), (2, -2))}

if __name__ == "__main__":
    case, A, B, Rs, T, dps = (sys.argv[1], sp.Rational(sys.argv[2]), sp.Rational(sys.argv[3]), sys.argv[4],
                              int(sys.argv[5]), int(sys.argv[6]))
    (l1, m1), (l2, m2) = LM[case]
    orb = ((l1 + 1) * A, (l1 + 1, l1, m1), (1, 0, 0), B, (l2 + 1, l2, m2), (1, 0, 0))
    t0 = time.time()
    mp.mp.dps = dps
    v, per = Q.exchange_hp(*orb, sp.Rational(Rs), tau_max=T, dps=dps)
    out = {"case": case, "pair": [str(A), str(B)], "R": Rs, "tau_max": T, "dps": dps, "orbitals": str(orb),
           "value": mp.nstr(v, dps - 2), "per_tau": [mp.nstr(x, 5) for x in per],
           "time_s": round(time.time() - t0, 1)}
    tg = f"{case}_{A}_{B}_{Rs}_T{T}".replace("/", "-")
    json.dump(out, open(f"debug/data/exch_reverse_review_hp_{tg}.json", "w"), indent=1)
    print(out, flush=True)
