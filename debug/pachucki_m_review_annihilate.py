"""Review: apply the COMMITTED operators (OP_P_1_2, OP_3DDEL_1_2, read as text from the tests) to the closed form K
rebuilt by the review driver, as an exact identity in the function algebra (no operator search, no comparison with
the memo's derived operator).  Also: the operator for the WRONG family must not annihilate.
Run: python debug/pachucki_m_review_annihilate.py"""
import ast
import json
import re
import sys
from fractions import Fraction

import sympy as sp

sys.path.insert(0, "debug")
sys.path.insert(0, "tests")
import pachucki_ext_core as E  # noqa: E402
import pachucki_m_core as M  # noqa: E402
import pachucki_m_review_plants as RP  # noqa: E402

Rs = sp.Symbol("R")


def read_op(fn, name):
    src = open(fn, encoding="utf-8").read()
    m = re.search(r"^%s = (\[.*?^\])" % name, src, re.S | re.M)
    return [sp.Poly(sp.sympify(s_, locals={"R": Rs}), Rs) for s_ in ast.literal_eval(m.group(1))]


def build_K(ctx, case):
    # same assembly as RP.case_value, returning K
    import pachucki_m_review_plants as P
    A_, B_ = sp.Rational(ctx.a.numerator, ctx.a.denominator), sp.Rational(ctx.b.numerator, ctx.b.denominator)
    s1A, s1B, s2A, s2B, v, R = sp.symbols("s1A s1B s2A s2B v R", positive=True)
    z1 = (s1A + R**2 - s1B) / (2 * R)
    z2 = (s2A + R**2 - s2B) / (2 * R)
    X = (s1A + s2A - v) / 2 - z1 * z2
    SPEC = {"2ppi": (X, -A_**5 * B_**3 / 2),
            "3ddel": (2 * X**2 - (s1A - z1**2) * (s2A - z2**2), A_**7 * B_**3 / 12)}
    Pol, N2 = SPEC[case]
    poly = sp.Poly(sp.expand(Pol), s1A, s1B, s2A, s2B, v)
    terms = {}
    for (i1, j1, i2, j2, k), co in poly.terms():
        terms.setdefault(int(k), []).append(((2 * i1 + 1, 2 * j1 + 1, 2 * i2 + 1, 2 * j2 + 1), co))
    kmax = max(terms)
    Dk, gens = {}, []
    for k in range(kmax, -1, -1):
        gens = gens + [al for al, _ in terms.get(k, [])]
        Dk[k] = E.downset(gens)
    calc = M.MCalc(ctx)
    _, f0, _, _, _ = E.solve_all(ctx, Dk[0], check=True)
    fl = {0: f0}
    for k in range(1, kmax + 1):
        RJ = M.RJets(ctx, Dk[k], calc)
        fl[k] = P.solve_fn_plant(calc, RJ, Dk[k], 2 * k - 1, fl[k - 1], fl[k - 2] if k >= 2 else None, "none")
    pref = Fraction(str(16 * N2))
    K = {}
    for k, lst in terms.items():
        for al, cR in lst:
            for term in sp.Add.make_args(sp.expand(cR)):
                cc, e = term.as_coeff_exponent(R)
                cc = sp.Rational(cc)
                K = E.fadd(K, E.fscale(fl[k][al], pref * Fraction(int(cc.p), int(cc.q)), int(e) - 1))
    return calc, K


def apply(calc, K, op):
    ders = [K]
    for _ in range(len(op) - 1):
        ders.append(calc.diff(ders[-1]))
    tot = {}
    for i, p in enumerate(op):
        for (j,), c in p.terms():
            tot = E.fadd(tot, E.fscale(ders[i], Fraction(int(c)), j))
    return tot


ctx = E.Ctx("Q", 1, 2)
ops = {"2ppi": read_op("tests/test_paper58_exchange_l1.py", "OP_P_1_2"),
       "3ddel": read_op("tests/test_paper58_exchange_d_pp.py", "OP_3DDEL_1_2")}
out = {}
for case in ("2ppi", "3ddel"):
    calc, K = build_K(ctx, case)
    for opn, op in ops.items():
        res = apply(calc, K, op)
        out[f"{opn}_operator_on_{case}"] = {"annihilates_exactly": res == {}, "n_residual_keys": len(res),
                                            "order": len(op) - 1, "degree": max(p.degree() for p in op)}
        print(case, "<-", opn, out[f"{opn}_operator_on_{case}"], flush=True)
json.dump(out, open("debug/data/pachucki_m_review_annihilate.json", "w"), indent=1)
