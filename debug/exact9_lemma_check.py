"""exact9_lemma_check.py -- numerical checks of the local structure at R = 0 used by the independence lemma.

Claim used: on R > 0 each of the nine functions is  A_i(R) ln R + B_i(R)  with A_i, B_i entire, and
    A(f0) = -2c [Cinh(2bR) - Cinh(2aR)],  c = 1/(2(a^2-b^2)),  Cinh(x) = sum_{k>=1} x^{2k} / (2k (2k)!),
    A(Ei(-(a+b)R)) = 1, A(e^{2aR}Ei(-2(a+b)R)) = e^{2aR}, A(e^{2bR}Ei(..)) = e^{2bR},
    A(psi_a) = -e^{-2aR}, A(psi_b) = -e^{-2bR}, A(exponentials) = 0.
Checks (a = 1, b = t; t in 2, 3, 3/2, 5/4):
  C1  FM1 (label 60, helper's Func, mp.ei) = L ln R + M on R > 0, L = c(e^{2bR}-e^{2aR}+e^{-2bR}-e^{-2aR}),
      M built from the entire parts Ein; and M(0) = 0.
  C2  monodromy of f0: the contour integral of f0' = -(L log R + M)/R once around R = 0 (log continued) equals
      2 pi i A(f0)(R0).
  C3  B = f0 - A ln R (f0 by the helper's quadrature of label 59) satisfies B(R2) - B(R1) =
      -int_{R1}^{R2} (M + A)/s ds, an integral of an entire integrand (no log).
  C4  psi_a, psi_b, g6, g7, g2: g - A ln R equals the explicit entire series (definition of Ei on both signs).
Writes debug/data/exact9_lemma_check.json
"""
import json
import sys
from fractions import Fraction

import mpmath as mp

sys.path.insert(0, "tests")
import _pachucki_extension as PE  # noqa: E402

mp.mp.dps = 40
out = {}


def ein(z):          # Ei(x) = gamma + ln|x| + ein(x), entire
    s, term, k = mp.mpf(0), mp.mpf(1), 0
    while True:
        k += 1
        term = term * z / k
        add = term / k
        s += add
        if k > 20 and abs(add) < mp.mpf(10) ** (-mp.mp.dps - 5):
            return s


def cinh(x):
    return (ein(x) + ein(-x)) / 2


for tq in (Fraction(2), Fraction(3), Fraction(3, 2), Fraction(5, 4)):
    a, b = mp.mpf(1), mp.mpf(tq.numerator) / tq.denominator
    c = 1 / (2 * (a * a - b * b))
    LC = mp.log((a + b) / (b - a))
    g_ = mp.euler

    def Eent(mu, R):     # entire part of Ei(mu R): Ei(mu R) = ln R + Eent
        return g_ + mp.log(abs(mu)) + ein(mu * R)

    def L(R):
        return c * (mp.exp(2 * b * R) - mp.exp(2 * a * R) + mp.exp(-2 * b * R) - mp.exp(-2 * a * R))

    def M(R):
        return c * (mp.exp(2 * b * R) * Eent(-2 * (a + b), R) - mp.exp(2 * a * R) * Eent(-2 * (a + b), R)
                    - mp.exp(-2 * b * R) * (Eent(2 * (b - a), R) - 2 * Eent(b - a, R) - LC)
                    + mp.exp(-2 * a * R) * (Eent(2 * (a - b), R) - 2 * Eent(a - b, R) - LC))

    def A1(R):
        return -2 * c * (cinh(2 * b * R) - cinh(2 * a * R))

    calc = PE.Calc(PE.Ctx("Q", 1, tq))
    res = {}
    # C1
    d1 = max(abs(PE.evaluate(Fraction(1), tq, calc.FM1, x) - (L(x) * mp.log(x) + M(x)))
             for x in (mp.mpf("0.3"), mp.mpf("1.1"), mp.mpf("2.7")))
    res["C1_FM1_minus_(L lnR + M)"] = mp.nstr(d1, 5)
    res["C1_M(0)"] = mp.nstr(abs(M(mp.mpf(0))), 5)
    # C2
    R0 = mp.mpf("0.7")

    def integrand(th):
        z = R0 * mp.expj(th)
        logz = mp.log(R0) + 1j * th
        return -(L(z) * logz + M(z)) / z * (1j * z)

    mono = mp.quad(integrand, mp.linspace(0, 2 * mp.pi, 9))
    res["C2_monodromy_minus_2piiA"] = mp.nstr(abs(mono - 2j * mp.pi * A1(R0)), 5)
    res["C2_2piiA_size"] = mp.nstr(abs(2 * mp.pi * A1(R0)), 8)
    # C3
    R1, R2 = mp.mpf("0.05"), mp.mpf("1.3")
    f1, f2 = PE.f0_quadrature(Fraction(1), tq, R1), PE.f0_quadrature(Fraction(1), tq, R2)
    lhs = (f2 - A1(R2) * mp.log(R2)) - (f1 - A1(R1) * mp.log(R1))
    rhs = -mp.quad(lambda s: (M(s) + A1(s)) / s, mp.linspace(R1, R2, 6))
    res["C3_B_increment_diff"] = mp.nstr(abs(lhs - rhs), 5)
    res["C3_B_increment_size"] = mp.nstr(abs(lhs), 8)
    res["C3_ln_term_size"] = mp.nstr(abs(A1(R2) * mp.log(R2) - A1(R1) * mp.log(R1)), 8)
    # C4
    x = mp.mpf("0.9")
    nine = PE.nine_functions(PE.Ctx("Q", 1, tq))
    val = [None] + [PE.evaluate(Fraction(1), tq, g, x) for g, _ in nine[1:]]
    Avals = [None, 1, 0, 0, 0, mp.exp(2 * a * x), mp.exp(2 * b * x), -mp.exp(-2 * a * x), -mp.exp(-2 * b * x)]
    Bvals = [None, Eent(-(a + b), x), mp.exp(-(a + b) * x), mp.exp(-2 * a * x), mp.exp(-2 * b * x),
             mp.exp(2 * a * x) * Eent(-2 * (a + b), x), mp.exp(2 * b * x) * Eent(-2 * (a + b), x),
             mp.exp(-2 * a * x) * (Eent(2 * (a - b), x) - 2 * Eent(a - b, x) - LC),
             mp.exp(-2 * b * x) * (Eent(2 * (b - a), x) - 2 * Eent(b - a, x) - LC)]
    d4 = max(abs(val[i] - (Avals[i] * mp.log(x) + Bvals[i])) for i in range(1, 9))
    res["C4_max_diff_g2..g9"] = mp.nstr(d4, 5)
    out[str(tq)] = res
    print(tq, res, flush=True)

json.dump(out, open("debug/data/exact9_lemma_check.json", "w"), indent=1)
