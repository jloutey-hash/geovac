"""Controls (b) and (c)-specialisation for an m != 0 case.
 (b) closed form (dps 60; the only quadrature is the 1-D f0 integral) against exchange_hp reference JSONs;
 (c') the Q(t) closed form specialised at t = 2 equals the Q run at (1,2) EXACTLY (nine coordinates), and specialised
      at t = 15/8 with the scaling K(R;a,b) = a K(aR;1,b/a) equals the Q run at (4/5,3/2), R = 23/10, numerically.
Run: python debug/pachucki_m_verify.py CASE SYMTAG      (SYMTAG = <case>_sym or <case>_sym_res)"""
import glob
import json
import pickle
import sys
from fractions import Fraction

import mpmath as mp
import sympy as sp

sys.path.insert(0, "debug")
sys.path.insert(0, "tests")
import pachucki_ext_core as E  # noqa: E402
import pachucki_m_core as M  # noqa: E402
import _pachucki_exchange as PX  # noqa: E402

case = sys.argv[1]
symtag = sys.argv[2] if len(sys.argv) > 2 else None
mp.mp.dps = 60
out = {"case": case, "hp": [], "sym": {}}


def digits(x, y):
    d = abs(x - y)
    return 60 if d == 0 else int(mp.floor(-mp.log10(d / abs(y))))


Ks = {}
for fn in sorted(glob.glob(f"debug/data/pachucki_m_hpref_{case}_*_T24.json")):
    J = json.load(open(fn))
    a, b = Fraction(J["pair"][0]), Fraction(J["pair"][1])
    tg = f"{case}_{a}_{b}".replace("/", "-")
    if tg not in Ks:
        Ks[tg] = pickle.load(open(f"debug/data/pachucki_m_K_{tg}.pkl", "rb"))["K"]
    alg = PX.Algebra(a, b)
    x = mp.mpf(J["R"])
    cf = PX.evaluate(alg, Ks[tg], x, PX.f0_numeric(alg, x))
    hp = mp.mpf(J["value"])
    rec = {"pair": J["pair"], "R": J["R"], "tau_max": J["tau_max"], "hp_dps": J["dps"], "hp_last_terms": J["last_terms"],
           "closed_form": mp.nstr(cf, 52), "exchange_hp": J["value"], "abs_diff": mp.nstr(abs(cf - hp), 3),
           "agreeing_digits": digits(cf, hp)}
    out["hp"].append(rec)
    print(rec, flush=True)

if symtag:
    S = json.load(open(f"debug/data/pachucki_m_case_{symtag}.json"))
    t = sp.Symbol("t")
    coords = [{int(p): sp.sympify(c, locals={"t": t}) for p, c in d.items()} for d in S["nine_coords"]]

    def spec(tv):
        ctx = E.Ctx("Q", 1, tv)
        nine, _ = M.nine_span(ctx)
        K = {}
        for g, c in zip(nine, coords):
            for p, x_ in c.items():
                v = sp.Rational(x_.subs(t, sp.Rational(tv.numerator, tv.denominator)))
                K = E.fadd(K, E.fscale(g, Fraction(int(v.p), int(v.q)), p))
        return K
    K2 = spec(Fraction(2))
    Kq = pickle.load(open(f"debug/data/pachucki_m_K_{case}_1_2.pkl", "rb"))["K"]
    out["sym"]["t=2_equals_Q_run_exactly"] = E.fadd(K2, Kq, -1) == {}
    print("Q(t) closed form at t = 2 == Q run at (1,2), exactly:", out["sym"]["t=2_equals_Q_run_exactly"], flush=True)
    a, b, Rv = Fraction(4, 5), Fraction(3, 2), Fraction(23, 10)
    Kt = spec(b / a)
    alg1 = PX.Algebra(Fraction(1), b / a)
    x1 = mp.mpf((a * Rv).numerator) / (a * Rv).denominator
    v1 = mp.mpf(a.numerator) / a.denominator * PX.evaluate(alg1, Kt, x1, PX.f0_numeric(alg1, x1))
    Kd = pickle.load(open(f"debug/data/pachucki_m_K_{case}_4-5_3-2.pkl", "rb"))["K"]
    alg2 = PX.Algebra(a, b)
    x2 = mp.mpf(23) / 10
    v2 = PX.evaluate(alg2, Kd, x2, PX.f0_numeric(alg2, x2))
    out["sym"]["scaled_t=15/8_vs_Q_run_(4/5,3/2)_R=2.3"] = {"scaled": mp.nstr(v1, 52), "direct": mp.nstr(v2, 52),
                                                           "agreeing_digits": digits(v1, v2)}
    print(out["sym"], flush=True)
    for kk in ("denominator_factors_of_K_coefficients", "max_pole_order_at_t_1", "rank_at_t3_R7_5", "K_in_nine_span",
               "derivatives_in_nine_span"):
        out["sym"][kk] = S.get(kk)
json.dump(out, open(f"debug/data/pachucki_m_verify_{case}.json", "w"), indent=1)
