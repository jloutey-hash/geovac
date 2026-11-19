"""Review: spot-check Q(t) intermediate jets of the memo's resumable runs by specialising at a NEW t and comparing
with a direct Q computation at (1, t) (separate arithmetic path, separate solve).
Run: python debug/pachucki_m_review_qt.py CASE T k1:alpha k2:alpha ...   (alpha like 2011)"""
import json
import pickle
import sys
from fractions import Fraction

import sympy as sp

sys.path.insert(0, "debug")
import pachucki_ext_core as E  # noqa: E402
import pachucki_m_core as M  # noqa: E402

case, tv = sys.argv[1], Fraction(sys.argv[2])
picks = [(int(s.split(":")[0]), tuple(int(c) for c in s.split(":")[1])) for s in sys.argv[3:]]
t = sp.Symbol("t")
raw = pickle.load(open(f"debug/data/pachucki_m_symres_{case}.pkl", "rb"))
T = sp.Rational(tv.numerator, tv.denominator)


def q(s):
    v = sp.Rational(sp.sympify(s, locals={"t": t}).subs(t, T))
    return Fraction(int(v.p), int(v.q))


def spec(Fn):
    out = {}
    for key, d in Fn.items():
        if key[0] == "f":
            kk = key
        elif key[0] == "E":
            kk = ("E", q(key[1]), key[2])
        else:
            kk = ("Ei", q(key[1]), q(key[2]))
        dd = {p: q(c) for p, c in d.items()}
        out = E.fadd(out, {kk: dd})
    return E.clean(out)


ctx = E.Ctx("Q", 1, tv)
need = max(k for k, _ in picks)
gens = [al for _, al in picks]
D = E.downset(gens)
calc = M.MCalc(ctx)
_, f0, _, _, _ = E.solve_all(ctx, D, check=True)
fl = {0: f0}
for k in range(1, need + 1):
    RJ = M.RJets(ctx, D, calc)
    fl[k] = M.solve_fn(calc, RJ, D, 2 * k - 1, fl[k - 1], fl[k - 2] if k >= 2 else None, check=True)
res = []
for k, al in picks:
    sym = spec(raw[k][al])
    direct = fl[k][al]
    same = E.fadd(sym, direct, -1) == {}
    res.append({"k": k, "alpha": list(al), "n_keys": len(direct), "equal": same})
    print(case, "t =", tv, "f_%d" % (2 * k), al, "keys", len(direct), "Q(t)-specialised == direct Q:", same, flush=True)
out = {"case": case, "t": str(tv), "checks": res}
json.dump(out, open(f"debug/data/pachucki_m_review_qt_{case}_{str(tv).replace('/', '-')}.json", "w"), indent=1)
