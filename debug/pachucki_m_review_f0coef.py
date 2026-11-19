"""Review: the memo's pinned f0 coefficient of 2ppi against the rebuilt closed form at (1,2) and (2/3,7/4)."""
import json, sys
from fractions import Fraction
import sympy as sp
sys.path.insert(0, "debug"); sys.path.insert(0, "tests")
import pachucki_ext_core as E
import pachucki_m_review_annihilate as A
R, a_, b_ = sp.symbols("R a b")
form = 2*a_**7*b_**5*(225 - 180*(a_**2+b_**2)*R**2 - (102*a_**4+336*a_**2*b_**2+102*b_**4)*R**4
       - (24*a_**6+176*a_**4*b_**2+176*a_**2*b_**4+24*b_**6)*R**6
       + (a_**8+16*a_**6*b_**2+36*a_**4*b_**4+16*a_**2*b_**6+b_**8)*R**8)/(R**3*(b_**2-a_**2)**6)
out = {}
for a, b in ((Fraction(1), Fraction(2)), (Fraction(2, 3), Fraction(7, 4))):
    ctx = E.Ctx("Q", a, b)
    calc, K = A.build_K(ctx, "2ppi")
    mine = sum(sp.Rational(c.numerator, c.denominator) * R**p for p, c in K[("f",)].items())
    ref = form.subs({a_: sp.Rational(a.numerator, a.denominator), b_: sp.Rational(b.numerator, b.denominator)})
    out[f"{a}_{b}"] = bool(sp.simplify(mine - ref) == 0)
    print(a, b, out[f"{a}_{b}"], flush=True)
json.dump(out, open("debug/data/pachucki_m_review_f0coef.json", "w"), indent=1)
