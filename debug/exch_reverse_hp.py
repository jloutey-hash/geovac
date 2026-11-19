"""exch_reverse_hp.py -- reference values from geovac.qfd_core.exchange_hp (Neumann sum over tau; never sees the
closed form) for the four mixed exchange integrals with the l > 0 orbital on A carrying the exponent a, a > b.

Run: python debug/exch_reverse_hp.py CASE A B R TAU DPS
     -> debug/data/exch_reverse_hp_<CASE>_<A>_<B>_<R>_T<TAU>.json
CASE in 2psig, 3dsig (m = 0), 2ppi, 3ddel (m != 0); orbitals as in debug/eq_gen.py SPECS / ref_orbitals
(hydrogenic Z = n * rate)."""
import json
import sys
import time

import mpmath as mp
import sympy as sp

sys.path.insert(0, ".")
from geovac import qfd_core as Q  # noqa: E402

ORB = {
    "2psig": lambda a, b: (2 * a, (2, 1, 0), (1, 0, 0), b, (2, 1, 0), (1, 0, 0)),
    "3dsig": lambda a, b: (3 * a, (3, 2, 0), (1, 0, 0), b, (3, 2, 0), (1, 0, 0)),
    "2ppi": lambda a, b: (2 * a, (2, 1, 1), (1, 0, 0), b, (2, 1, -1), (1, 0, 0)),
    "3ddel": lambda a, b: (3 * a, (3, 2, 2), (1, 0, 0), b, (3, 2, -2), (1, 0, 0)),
}

if __name__ == "__main__":
    case, A, B, Rs, T, dps = (sys.argv[1], sp.Rational(sys.argv[2]), sp.Rational(sys.argv[3]), sys.argv[4],
                              int(sys.argv[5]), int(sys.argv[6]))
    t0 = time.time()
    orb = ORB[case](A, B)
    mp.mp.dps = dps
    v, per = Q.exchange_hp(*orb, sp.Rational(Rs), tau_max=T, dps=dps)
    out = {"case": case, "pair": [str(A), str(B)], "R": Rs, "tau_max": T, "dps": dps, "orbitals": str(orb),
           "value": mp.nstr(v, dps - 2), "last_terms": [mp.nstr(x, 3) for x in per[-4:]],
           "time_s": round(time.time() - t0, 1)}
    tg = f"{case}_{A}_{B}_{Rs}_T{T}".replace("/", "-")
    json.dump(out, open(f"debug/data/exch_reverse_hp_{tg}.json", "w"), indent=1)
    print(out, flush=True)
