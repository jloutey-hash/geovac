"""Checks on EVERY intermediate jet of a Q(t) run (state pickle of pachucki_m_symres.py): the proof 'for every t > 1'
uses uniqueness by decay at every multi-index, so every f_{n,alpha} must (i) have coefficient denominators without
zeros in t > 1 and (ii) consist of keys that decay for every t > 1.
Run: python debug/pachucki_m_symcheck.py CASE"""
import json
import pickle
import sys

import sympy as sp

case = sys.argv[1]
t = sp.Symbol("t")
raw = pickle.load(open(f"debug/data/pachucki_m_symres_{case}.pkl", "rb"))
den_factors, bad_keys, nkeys = set(), [], set()
cache = {}


def dfac(s):
    if s not in cache:
        d = sp.fraction(sp.together(sp.sympify(s)))[1]
        cache[s] = {str(f) for f, _ in sp.factor_list(d)[1]}
    return cache[s]


def decays(e):
    """e (linear in t) < 0 for all t > 1."""
    e = sp.sympify(e)
    p = sp.Poly(e, t)
    assert p.degree() <= 1
    return bool(e.subs(t, 1) <= 0 and p.coeff_monomial(t) <= 0 and not (e.subs(t, 1) == 0 and p.coeff_monomial(t) == 0))


count = 0
for k, d in raw.items():
    for al, Fn in d.items():
        count += 1
        for key, co in Fn.items():
            nkeys.add(key)
            if key[0] == "E":
                ok = decays(key[1])
            elif key[0] == "Ei":
                ok = decays(sp.sympify(key[1]) + sp.sympify(key[2]))
            else:
                ok = True
            if not ok:
                bad_keys.append((k, al, key))
            for p, c in co.items():
                den_factors |= dfac(c)
out = {"case": case, "jets_checked": count, "sets": {str(k): len(d) for k, d in raw.items()},
       "denominator_factors_all_intermediates": sorted(den_factors),
       "non_decaying_keys": [str(x) for x in bad_keys], "distinct_keys": len(nkeys)}
print(out)
json.dump(out, open(f"debug/data/pachucki_m_symcheck_{case}.json", "w"), indent=1)
