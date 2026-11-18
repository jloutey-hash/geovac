"""Reference values: geovac.qfd_core.exchange_hp for an EQ/EP case at one pair and R.
Run: python debug/pachucki_ext_hpref.py CASE A B R TAU DPS  -> debug/data/pachucki_ext_hpref_<CASE>_<A>_<B>_<R>_T<TAU>.json"""
import json, sys, time
import mpmath as mp
import sympy as sp
sys.path.insert(0, "."); sys.path.insert(0, "debug")
from geovac import qfd_core as Q
import eq_gen as G
case, A, B, Rs, T, dps = sys.argv[1], sp.Rational(sys.argv[2]), sp.Rational(sys.argv[3]), sys.argv[4], int(sys.argv[5]), int(sys.argv[6])
t0 = time.time()
orb = G.ref_orbitals(case, A, B)
mp.mp.dps = dps
v, per = Q.exchange_hp(*orb, sp.Rational(Rs), tau_max=T, dps=dps)
out = {"case": case, "pair": [str(A), str(B)], "R": Rs, "tau_max": T, "dps": dps, "orbitals": str(orb),
       "value": mp.nstr(v, dps - 2), "last_terms": [mp.nstr(x, 3) for x in per[-3:]], "time_s": round(time.time() - t0, 1)}
tg = f"{case}_{A}_{B}_{Rs}_T{T}".replace("/", "-")
json.dump(out, open(f"debug/data/pachucki_ext_hpref_{tg}.json", "w"), indent=1)
print(out)
