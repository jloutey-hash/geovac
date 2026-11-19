"""Review, attacks 3, 6, 7 on the UNPATCHED evaluator (debug/hybrid_ratio_engine.py is NOT imported).

  3. atom labels at the +a_d shift, read from hybrid_closed_form output at fresh points (H2 lt, H4 lt):
     which of E1((a_c + a_d)R) / E1((a_p + a_d)R) / E1((a_p + a_c + a_d)R) occur.
  6. exact minimal order (tests/_minimal_order.py: interval-certified minor + generic rank over Q(R)) at a
     coincidence a_d = a_c and a generic point, at exponents not used before (H7), plus H4 at a_d = a_p + a_c.
  7. scaling f(R; sZ_A, sZ_B) = s f(sR; Z_A, Z_B) (normalised orbitals), fam1 and H7, 40 digits.
  +  one double-precision comparison against geovac's independent 2-D quadrature (hybrid_quadrature).

Usage: python debug/hybrid_ratio_review_points.py <task>   task in {labels, order, scaling, quad}
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction as F
from pathlib import Path

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
import _minimal_order as M  # noqa: E402
from geovac.two_center_eri import hybrid_closed_form, hybrid_quadrature, R_s  # noqa: E402

O = {"1s": (1, 0, 0), "2s": (2, 0, 0), "2p0": (2, 1, 0), "2p1": (2, 1, 1), "3d0": (3, 2, 0)}
CASES = {"fam1": ("2p0", "2p0", "1s", "1s"), "H2": ("2p1", "2p1", "2p1", "2p1"),
         "H4": ("2p0", "2p0", "2p0", "3d0"), "H7": ("2p0", "3d0", "2s", "1s")}


def rates(case, za, zb):
    a, b, c, d = (O[x] for x in CASES[case])
    za, zb = sp.Rational(za), sp.Rational(zb)
    return za / a[0] + za / b[0], za / c[0], zb / d[0]          # a_p, a_c, a_d


def build(case, za, zb):
    a, b, c, d = (O[x] for x in CASES[case])
    t0 = time.time()
    e = hybrid_closed_form(F(za), a, b, c, F(zb), d)
    return e, time.time() - t0


def save(tag, obj):
    p = ROOT / "debug/data/hybrid_ratio_review_points.json"
    allc = json.loads(p.read_text()) if p.exists() else {}
    allc[tag] = obj
    p.write_text(json.dumps(allc, indent=1, default=str))


def task_labels():
    out = {}
    for case, za, zb in (("H2", "2", "1"), ("H4", "2", "3/2"), ("H7", "3", "1")):
        ap, ac, ad = rates(case, za, zb)
        e, tb = build(case, F(za), F(zb))
        f = M.Fn(e)
        labs = f.labels()
        plus = sorted({(str(c), str(k)) for (c, kind, k) in labs if kind == "y" and c == ad})
        minus = sorted({(str(c), str(k)) for (c, kind, k) in labs if kind == "y" and c == -ad})
        # E1(xR) = -Ei(-xR): kappa = -x
        out[f"{case} Z=({za},{zb})"] = {
            "a_p,a_c,a_d": [str(ap), str(ac), str(ad)], "plus_shift_kappas": plus, "minus_shift_kappas": minus,
            "E1((a_c+a_d)R) present": any(sp.Rational(k) == -(ac + ad) for _c, k in plus),
            "E1((a_p+a_c+a_d)R) present": any(sp.Rational(k) == -(ap + ac + ad) for _c, k in plus),
            "E1((a_p+a_d)R) present": any(sp.Rational(k) == -(ap + ad) for _c, k in plus),
            "rates": [str(r) for r in f.rates()], "build_s": round(tb, 1)}
        print(case, out[f"{case} Z=({za},{zb})"], flush=True)
    save("labels", out)


def task_order():
    out = {}
    pts = [("H7", "7", "7/2", 4),        # a_c = 7/2 = a_d            (first coincidence), fresh
           ("H7", "7", "13/5", 5),       # a_d < a_c, fresh generic
           ("H4", "4/3", "6", 4),        # a_p = 4/3, a_c = 2/3, a_d = 2 = a_p + a_c (second coincidence), fresh
           ("H2", "5", "17/3", 5)]       # a_p = 5, a_c = 5/2, a_d = 17/6 (middle), fresh
    for case, za, zb, expect in pts:
        ap, ac, ad = rates(case, za, zb)
        e, tb = build(case, F(za), F(zb))
        f = M.Fn(e)
        t0 = time.time()
        lo = M.order_at_least(f, expect)
        not_more = not M.order_at_least(f, expect + 1)
        up = M.order_at_most(f, expect)
        res = {"a_p,a_c,a_d": [str(ap), str(ac), str(ad)], "expected": expect, "order_at_least": lo,
               "order_at_most": up, "minor_of_size_expected+1_nonzero": not not_more,
               "n_rates": len(f.rates()), "rates": [str(r) for r in f.rates()],
               "labels": [str(x) for x in f.labels()], "build_s": round(tb, 1), "rank_s": round(time.time() - t0, 1)}
        out[f"{case} Z=({za},{zb})"] = res
        print(case, za, zb, res, flush=True)
    save("order", out)


def task_scaling():
    mp.mp.dps = 50
    out = {}
    for case, za, zb, s in (("fam1", "3", "5", "3/2"), ("H7", "2", "7/3", "5/4")):
        s_ = sp.Rational(s)
        e1, _ = build(case, F(za), F(zb))
        e2, _ = build(case, F(za) * F(s), F(zb) * F(s))
        devs = []
        for Rv in (sp.Rational(3, 4), sp.Rational(5, 2)):
            lhs = mp.mpf(sp.N(e2.subs(R_s, Rv), 45))
            rhs = mp.mpf(sp.N(s_ * e1.subs(R_s, s_ * Rv), 45))
            devs.append(mp.nstr(abs(lhs / rhs - 1), 5))
        out[case] = {"Z": [za, zb], "s": s, "rel_dev": devs}
        print(case, out[case], flush=True)
    save("scaling", out)


def task_quad():
    out = {}
    for case, za, zb, Rv in (("H7", "2", "7/3", 1.3), ("H4", "2", "3/2", 1.7)):
        a, b, c, d = (O[x] for x in CASES[case])
        e, _ = build(case, F(za), F(zb))
        cf = float(sp.N(e.subs(R_s, sp.nsimplify(Rv)), 30))
        q = hybrid_quadrature(F(za), a, b, c, F(zb), d, Rv)
        out[case] = {"Z": [za, zb], "R": Rv, "closed": cf, "quad": q, "rel": abs(cf / q - 1)}
        print(case, out[case], flush=True)
    save("quad", out)


if __name__ == "__main__":
    {"labels": task_labels, "order": task_order, "scaling": task_scaling, "quad": task_quad}[sys.argv[1]]()
