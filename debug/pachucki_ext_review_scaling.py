"""Reviewer check (b): the stored Q(t) closed forms + the scaling K(R; a, b) = a K(aR; 1, b/a), evaluated at
(a, b, R) = (4/5, 3/2, 23/10) (t = 15/8, aR = 46/25) against the reviewer's Neumann integrals; also the (a, b)-restored
1s coefficients (closed_form_coefficients_ab), a wrong power of a (K = a^2 K(aR), K = K(aR)), and the a = b limit
formula (Sugiura form typed in debug/pachucki_ext_limit.py) against the Neumann value at a = b = 1, R = 7/5.
Reads only the author's JSON; f0 from the reviewer's Neumann code (no label-60 quadrature here).
"""
import json
import sys
from fractions import Fraction

import mpmath as mp
import sympy as sp

sys.path.insert(0, "debug")
from pachucki_ext_review_neumann import neumann  # noqa: E402

DPS = 40
a, b, R = Fraction(4, 5), Fraction(3, 2), Fraction(23, 10)
t, x = b / a, a * R


def funcs(A, B, X, f0):
    """the nine functions at exponents (A, B), distance X (names as in the author's JSON)."""
    LC = mp.log((A + B) / (B - A))
    return {"f0": f0, "-E1((a+b)R)": -mp.e1((A + B) * X), "e^{-(a+b)R}": mp.exp(-(A + B) * X),
            "e^{-2aR}": mp.exp(-2 * A * X), "e^{-2bR}": mp.exp(-2 * B * X),
            "-e^{2aR}E1(2(a+b)R)": -mp.exp(2 * A * X) * mp.e1(2 * (A + B) * X),
            "-e^{2bR}E1(2(a+b)R)": -mp.exp(2 * B * X) * mp.e1(2 * (A + B) * X),
            "psi_a": mp.exp(-2 * A * X) * (mp.ei(2 * (A - B) * X) - 2 * mp.ei((A - B) * X) - LC),
            "psi_b": mp.exp(-2 * B * X) * (mp.ei(2 * (B - A) * X) - 2 * mp.ei((B - A) * X) - LC)}


mp.mp.dps = DPS + 15
f0_1t, _ = neumann("f0", 1, t, x, DPS)            # f0(aR; 1, t)
f0_ab, _ = neumann("f0", a, b, R, DPS)            # f0(R; a, b)
print("f0 scaling  f0(R;a,b) / [a^-2 f0(aR;1,t)] - 1 =", mp.nstr(f0_ab / (f0_1t / (mp.mpf(a.numerator) / a.denominator) ** 2) - 1, 3))
mp.mp.dps = DPS + 15
T = sp.Symbol("t")
tm, xm, am, bm, Rm = [mp.mpf(v.numerator) / v.denominator for v in (t, x, a, b, R)]
g1 = funcs(mp.mpf(1), tm, xm, f0_1t)
for case, fn in (("1s", "pachucki_ext_symbolic_1s"), ("2psig", "pachucki_ext_sym_2psig"),
                 ("ppsig", "pachucki_ext_sym_ppsig"), ("pppi0", "pachucki_ext_sym_pppi0"),
                 ("3dsig", "pachucki_ext_sym_3dsig")):
    d = json.load(open("debug/data/%s.json" % fn))
    K1 = mp.mpf(0)
    for nm, c in d["closed_form_coefficients_t"].items():
        for p, v in c.items():
            cv = sp.sympify(v).subs(T, sp.Rational(t.numerator, t.denominator))
            K1 += mp.mpf(cv.p) / cv.q * xm ** int(p) * g1[nm]
    Knum, _ = neumann(case, a, b, R, DPS)
    mp.mp.dps = DPS + 15
    print("%-6s a K(aR;1,t) = %s  Neumann K(R;a,b) = %s  rel diff %s | wrong powers: a^2: %s, a^0: %s" % (
        case, mp.nstr(am * K1, 32), mp.nstr(Knum, 32), mp.nstr(abs(am * K1 / Knum - 1), 3),
        mp.nstr(abs(am ** 2 * K1 / Knum - 1), 3), mp.nstr(abs(K1 / Knum - 1), 3)))
    if case == "1s":
        A_, B_, R_ = sp.symbols("a b R", positive=True)
        gab = funcs(am, bm, Rm, f0_ab)
        Kab = mp.mpf(0)
        for nm, v in d["closed_form_coefficients_ab"].items():
            cv = sp.sympify(v, locals={"a": A_, "b": B_, "R": R_}).subs(
                {A_: sp.Rational(4, 5), B_: sp.Rational(3, 2), R_: sp.Rational(23, 10)})
            Kab += mp.mpf(cv.p) / cv.q * gab[nm]
        print("       (a,b)-restored coefficients: rel diff", mp.nstr(abs(Kab / Knum - 1), 3))

# a = b limit: Sugiura form as typed by the author, vs Neumann at a = b = 1
Rv = mp.mpf(7) / 5
e2, Lg, X2, X4 = mp.exp(-2 * Rv), mp.euler + mp.log(Rv), mp.e1(2 * Rv), mp.exp(2 * Rv) * mp.e1(4 * Rv)
S, Sp = 1 + Rv + Rv ** 2 / 3, 1 - Rv + Rv ** 2 / 3
sug = (-e2 * (-mp.mpf(25) / 8 + mp.mpf(23) / 4 * Rv + 3 * Rv ** 2 + Rv ** 3 / 3)
       + 6 / Rv * (e2 * S ** 2 * Lg - Sp ** 2 * X4 + 2 * S * Sp * X2)) / 5
K11, _ = neumann("1s", 1, 1, Fraction(7, 5), DPS)
mp.mp.dps = DPS + 15
print("a=b=1, R=7/5: Sugiura form", mp.nstr(sug, 35), " Neumann", mp.nstr(K11, 35), " rel diff", mp.nstr(abs(sug / K11 - 1), 3))
