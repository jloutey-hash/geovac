"""Closed form of Pachucki's two-centre master integral f0(R; a, b) at the exchange point (a < b) and checks.

f0(R) = int_R^oo FM1(s)/s ds,  FM1 = label 60 of arXiv:0905.1635 (u = a, w = b), as transcribed in
tests/_pachucki_extension.Calc.FM1.  FM1 = sum_i kap_i e^{lam_i s} Ei(mu_i s) + c LC (e^{-2bs} - e^{-2as}).

Termwise:  T(lam, mu; R) := int_R^oo e^{lam s} Ei(mu s) ds/s  (lam + mu < 0)
         = K(lam, mu) - g ln R - (1/2) ln^2 R - g Ein(lam R) - ln R Ein(lam R) + E2(lam R) - Phi(lam R, mu R)
   g = gamma + ln|mu|,  Ein(x) = sum_{k>=1} x^k/(k k!),  E2(x) = sum_{k>=1} x^k/(k^2 k!),
   Phi(x, y) = sum_{j>=0, k>=1} x^j y^k / (j! k k! (j+k)) = int_0^1 e^{xu} Ein(yu) du/u,
   K(lam, mu) = finite part at eps = 0 of int_0^oo s^{eps-1} e^{lam s} Ei(mu s) ds:
     mu < 0 (q = -mu, p = -lam):  K = -[(gamma + ln q)^2/2 + pi^2/12 + Li2(-p/q)]
     mu > 0 (p = -lam > mu):      K = -[(gamma + ln mu)^2/2 - 5 pi^2/12 + Re Li2(p/mu)]
Plus  int_R^oo e^{-ps}/s ds = E1(pR) for the LC terms.

Writes debug/data/master_weight_closed_form.json.  No paper/test edits.
"""
from __future__ import annotations

import json
import os
import sys
import time
from fractions import Fraction

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tests"))
import _pachucki_extension as PE  # noqa: E402


def terms(a, b):
    """[(kappa/c, lam, mu)] of label 60 at u = a, w = b (matches PE.Calc.FM1; checked below)."""
    return [(1, 2 * b, -2 * (a + b)), (-1, 2 * a, -2 * (a + b)), (-1, -2 * b, 2 * (b - a)),
            (2, -2 * b, b - a), (1, -2 * a, 2 * (a - b)), (-2, -2 * a, a - b)]


def Ein(x, N=None):
    return mp.nsum(lambda k: x ** k / (k * mp.factorial(k)), [1, mp.inf]) if x != 0 else mp.mpf(0)


def E2(x):
    return mp.nsum(lambda k: x ** k / (k ** 2 * mp.factorial(k)), [1, mp.inf]) if x != 0 else mp.mpf(0)


def Phi(x, y):
    """sum_{k>=1} y^k/(k k!) * m_k(x),  m_k(x) = int_0^1 e^{xu} u^{k-1} du = sum_j x^j/(j!(j+k))."""
    tot = mp.mpf(0)
    k = 1
    tiny = mp.mpf(10) ** (-mp.mp.dps - 5)
    while True:
        mk = mp.nsum(lambda j: x ** j / (mp.factorial(j) * (j + k)), [0, mp.inf])
        term = y ** k / (k * mp.factorial(k)) * mk
        tot += term
        if k > 5 and abs(term) < tiny * max(1, abs(tot)):
            break
        k += 1
    return tot


def Kconst(lam, mu):
    lam, mu = mp.mpf(lam), mp.mpf(mu)
    p = -lam
    if mu < 0:
        q = -mu
        return -((mp.euler + mp.log(q)) ** 2 / 2 + mp.pi ** 2 / 12 + mp.polylog(2, -p / q))
    return -((mp.euler + mp.log(mu)) ** 2 / 2 - 5 * mp.pi ** 2 / 12 + mp.re(mp.polylog(2, p / mu)))


def T_closed(lam, mu, R):
    lam, mu, R = mp.mpf(lam), mp.mpf(mu), mp.mpf(R)
    g = mp.euler + mp.log(abs(mu))
    L = mp.log(R)
    ein = Ein(lam * R)
    return (Kconst(lam, mu) - g * L - L ** 2 / 2 - g * ein - L * ein + E2(lam * R) - Phi(lam * R, mu * R))


def T_quad(lam, mu, R):
    lam, mu, R = mp.mpf(lam), mp.mpf(mu), mp.mpf(R)
    return mp.quad(lambda s: mp.exp(lam * s) * mp.ei(mu * s) / s, [R, R + 1, R + 3, R + 8, R + 25, mp.inf])


def f0_closed(a, b, R):
    a, b = mp.mpf(a), mp.mpf(b)
    c = 1 / (2 * (a * a - b * b))
    LC = mp.log((a + b) / (b - a))
    s = sum(k * T_closed(lam, mu, R) for k, lam, mu in terms(a, b))
    s += LC * (mp.e1(2 * b * R) - mp.e1(2 * a * R))
    return c * s


def f0_parts(a, b, R):
    """Split f0 = C(a,b) - Lam ln R + W1(R) + W2(R): W2 = -c sum kap Phi (the weight-two power series),
    W1 = the rest (Ein, E2 and the LC E1 tails minus their logs)."""
    a, b, R = mp.mpf(a), mp.mpf(b), mp.mpf(R)
    c = 1 / (2 * (a * a - b * b))
    LC = mp.log((a + b) / (b - a))
    C = c * (sum(k * Kconst(lam, mu) for k, lam, mu in terms(a, b)) + LC * mp.log(a / b))
    Lam = 2 * c * (mp.chi(2 * b * R) - mp.log(2 * b * R) - mp.euler - mp.chi(2 * a * R) + mp.log(2 * a * R) + mp.euler)
    W2 = -c * sum(k * Phi(lam * R, mu * R) for k, lam, mu in terms(a, b))
    W1 = c * sum(k * (-(mp.euler + mp.log(abs(mu))) * Ein(lam * R) + E2(lam * R)) for k, lam, mu in terms(a, b))
    W1 += c * LC * (Ein(-2 * a * R) - Ein(-2 * b * R))
    return C, Lam, W1, W2


def main():
    out = {"dps": 60, "checks": []}
    mp.mp.dps = 60
    t0 = time.time()
    # 0. the term list equals the helper's transcription of label 60
    for (a, b) in ((Fraction(1), Fraction(2)), (Fraction(1), Fraction(3, 2)), (Fraction(1), Fraction(3))):
        calc = PE.Calc(PE.Ctx("Q", a, b))
        c = 1 / (2 * (a * a - b * b))
        mine = {("Ei", lam, mu): {0: k * c} for k, lam, mu in terms(a, b)}
        mine[("E", -2 * b, 1)] = {0: c}
        mine[("E", -2 * a, 1)] = {0: -c}
        out.setdefault("terms_equal_helper_FM1", []).append([str(a), str(b), PE.clean(mine) == calc.FM1])
        # swapped-pair test used in the proof: for each Ei key (lam, mu), is (mu, lam) also a key?
        keys = {(lam, mu) for _, lam, mu in terms(a, b)}
        out.setdefault("swapped_pairs_present", []).append(
            [str(a), str(b), sorted([(str(l), str(m)) for (l, m) in keys if (m, l) in keys])])
    # 1. per-term closed form vs quadrature
    for (lam, mu) in ((4, -6), (2, -6), (-4, 2), (-4, 1), (-2, -2), (-2, -1), (-3, 1), (-3, 2)):
        for R in ("0.4", "1.7"):
            q, cf = T_quad(lam, mu, mp.mpf(R)), T_closed(lam, mu, mp.mpf(R))
            out["checks"].append({"what": "T(lam,mu;R)", "lam": lam, "mu": mu, "R": R,
                                  "quad": mp.nstr(q, 45), "closed": mp.nstr(cf, 45),
                                  "absdiff": mp.nstr(abs(q - cf), 5)})
    # 2. whole f0 vs the helper's quadrature (tests/_pachucki_extension.f0_quadrature)
    for (a, b) in ((Fraction(1), Fraction(2)), (Fraction(1), Fraction(3, 2)), (Fraction(2, 3), Fraction(5, 2))):
        for R in ("0.3", "1.4", "3.1"):
            Rq = mp.mpf(R)
            h = PE.f0_quadrature(a, b, Rq)
            cf = f0_closed(mp.mpf(a.numerator) / a.denominator, mp.mpf(b.numerator) / b.denominator, Rq)
            out["checks"].append({"what": "f0", "a": str(a), "b": str(b), "R": R, "helper_quad": mp.nstr(h, 45),
                                  "closed": mp.nstr(cf, 45), "absdiff": mp.nstr(abs(h - cf), 5),
                                  "reldiff": mp.nstr(abs(h - cf) / abs(h), 5)})
    # 3. decomposition C - Lam ln R + W1 + W2 reproduces f0; size of the weight-two part
    for (a, b) in ((mp.mpf(1), mp.mpf(2)),):
        for R in ("0.3", "1.4"):
            Rq = mp.mpf(R)
            C, Lam, W1, W2 = f0_parts(a, b, Rq)
            tot = C - Lam * mp.log(Rq) + W1 + W2
            out["checks"].append({"what": "parts", "a": "1", "b": "2", "R": R, "C": mp.nstr(C, 30),
                                  "Lam": mp.nstr(Lam, 30), "W1": mp.nstr(W1, 30), "W2": mp.nstr(W2, 30),
                                  "sum_minus_closed": mp.nstr(abs(tot - f0_closed(a, b, Rq)), 5)})
    # 4. R -> 0: f0(0) = C(a, b) (Lam ~ -R^2), check by quadrature at tiny R; ln R coefficient = -Lam
    for (a, b) in ((Fraction(1), Fraction(2)), (Fraction(1), Fraction(3, 2))):
        am, bm = mp.mpf(a.numerator) / a.denominator, mp.mpf(b.numerator) / b.denominator
        C = f0_parts(am, bm, mp.mpf(1))[0]
        Rs = mp.mpf("1e-12")
        h = PE.f0_quadrature(a, b, Rs)
        out["checks"].append({"what": "R->0", "a": str(a), "b": str(b), "C": mp.nstr(C, 40),
                              "helper_quad_at_1e-12": mp.nstr(h, 40), "diff": mp.nstr(abs(h - C), 5)})
    # 5. R -> oo: f0 ~ LC e^{-2aR} / (4a(b^2-a^2) R) (1 + O(1/R))
    for (a, b) in ((Fraction(1), Fraction(2)),):
        am, bm = mp.mpf(1), mp.mpf(2)
        LC = mp.log((am + bm) / (bm - am))
        for R in (10, 20, 40, 80):
            h = PE.f0_quadrature(a, b, mp.mpf(R))
            lead = LC * mp.exp(-2 * am * R) / (4 * am * (bm ** 2 - am ** 2) * R)
            out["checks"].append({"what": "R->oo", "R": R, "ratio_minus_1_times_R": mp.nstr((h / lead - 1) * R, 12)})
    out["seconds"] = round(time.time() - t0, 1)
    with open(os.path.join(HERE, "data", "master_weight_closed_form.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
