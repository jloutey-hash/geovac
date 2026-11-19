"""exch_reverse_decay.py -- key-by-key decay verdict on t in (0, 1) for every function key of every intermediate jet
of the frozen Q(t) solves (audit lists stored with the closed forms), and the denominator factors.
Run: python debug/exch_reverse_decay.py  -> debug/data/exch_reverse_decay.json"""
import json
import sys

import sympy as sp

sys.path.insert(0, "debug")
import exch_reverse_core as C  # noqa: E402

t = sp.Symbol("t")
out = {}
for case in C.CASES:
    au = C.stored(case)["audit"]
    rows = []
    allok = True
    for k in au["keys"]:
        if k["kind"] == "f":
            rows.append({"key": "f0", "rate_t": "-2*min(1,t)", "decays_on_(0,1)": True})
            continue
        lam = sp.sympify(k["lam"], locals={"t": t})
        mu = sp.sympify(k.get("mu", "0"), locals={"t": t})
        rate = sp.expand(lam + mu)
        r0, r1 = rate.subs(t, 0), rate.subs(t, 1)
        dec = bool(r0 <= 0 and r1 <= 0 and not (r0 == 0 and r1 == 0))     # linear in t
        zero_arg = [] if k["kind"] == "E" else [str(z) for z in sp.solve(mu, t)]
        sign_mu = None if k["kind"] == "E" else {"t<1": str(sp.sign(mu.subs(t, sp.Rational(1, 2)))),
                                                 "t>1": str(sp.sign(mu.subs(t, 2)))}
        rows.append({"key": k["kind"], "lam": str(lam), "mu": str(mu), "log_power": k.get("log_power"),
                     "rate": str(rate), "decays_on_(0,1)": dec, "Ei_argument_zero_at_t": zero_arg,
                     "sign_of_Ei_argument": sign_mu})
        allok &= dec
    dens = au["denominator_factors_max_multiplicity"]
    facs = sorted({f for d in dens.values() for f in d})
    out[case] = {"keys": rows, "all_decay_on_(0,1)": allok, "denominator_factors": facs}
    print(case, "all keys decay on (0,1):", allok, "denominator factors:", facs, flush=True)
    for r in rows:
        print("   ", r, flush=True)
json.dump(out, open("debug/data/exch_reverse_decay.json", "w"), indent=1)
