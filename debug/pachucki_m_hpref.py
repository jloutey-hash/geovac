"""Reference values from geovac.qfd_core.exchange_hp (Neumann sum; never sees the closed form) for m != 0 cases.
Run: python debug/pachucki_m_hpref.py CASE A B R TAU DPS -> debug/data/pachucki_m_hpref_<CASE>_<A>_<B>_<R>_T<TAU>.json
CASE: 2ppi  = (2p+1_A 1s_B | 2p-1_A 1s_B)   (orbital tuples as in tests/test_paper58_exchange_l1.py REF_ORBITALS)
      3ddel = (3d+2_A 1s_B | 3d-2_A 1s_B)   (as in tests/test_paper58_exchange_d_pp.py)
      ppdel = (2p+1_A 2p-1_B | 2p-1_A 2p+1_B)
      pppi0 = (2p+1_A 2p+1_B | 2p+1_A 2p+1_B)
      1s"""
import json
import sys
import time

import mpmath as mp
import sympy as sp

sys.path.insert(0, ".")
from geovac import qfd_core as Q  # noqa: E402

ORB = {
    "2ppi": lambda a, b: (2 * a, (2, 1, 1), (1, 0, 0), b, (2, 1, -1), (1, 0, 0)),
    "3ddel": lambda a, b: (3 * a, (3, 2, 2), (1, 0, 0), b, (3, 2, -2), (1, 0, 0)),
    "ppdel": lambda a, b: (2 * a, (2, 1, 1), (2, 1, -1), 2 * b, (2, 1, -1), (2, 1, 1)),
    "pppi0": lambda a, b: (2 * a, (2, 1, 1), (2, 1, 1), 2 * b, (2, 1, 1), (2, 1, 1)),
    "1s": lambda a, b: (a, (1, 0, 0), (1, 0, 0), b, (1, 0, 0), (1, 0, 0)),
}
case, A, B, Rs, T, dps = (sys.argv[1], sp.Rational(sys.argv[2]), sp.Rational(sys.argv[3]), sys.argv[4],
                          int(sys.argv[5]), int(sys.argv[6]))
t0 = time.time()
orb = ORB[case](A, B)
mp.mp.dps = dps
v, per = Q.exchange_hp(*orb, sp.Rational(Rs), tau_max=T, dps=dps)
out = {"case": case, "pair": [str(A), str(B)], "R": Rs, "tau_max": T, "dps": dps, "orbitals": str(orb),
       "value": mp.nstr(v, dps - 2), "last_terms": [mp.nstr(x, 3) for x in per[-3:]],
       "time_s": round(time.time() - t0, 1)}
tg = f"{case}_{A}_{B}_{Rs}_T{T}".replace("/", "-")
json.dump(out, open(f"debug/data/pachucki_m_hpref_{tg}.json", "w"), indent=1)
print(out, flush=True)
