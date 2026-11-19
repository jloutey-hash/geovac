"""Helper for tests/test_paper58_master_weight.py -- closed form of Pachucki's two-centre master integral f0.

Source of the integral (prior art): K. Pachucki, Phys. Rev. A 80, 032520 (2009), arXiv:0905.1635, labels 58-60:
    f0(R; a, b) = int_R^oo F^(-1)(s) ds / s     (label 59),  so  f0' = -F^(-1)/R,
    F^(-1) = c [ sum_i k_i e^{lam_i s} Ei(mu_i s) + LC (e^{-2bs} - e^{-2as}) ]     (label 60, u = a < w = b),
    c = 1/(2(a^2 - b^2)),  LC = ln((a+b)/(b-a)).
The reference transcription of label 60 is tests/_pachucki_extension.Calc.FM1; TERMS below must equal it.

Closed form (derived 2026-11-19; this module only evaluates it):
    T(lam, mu; R) = int_R^oo e^{lam s} Ei(mu s) ds/s
                  = K(lam, mu) - g ln R - (1/2) ln^2 R - g Ein(lam R) - ln R Ein(lam R) + E2(lam R) - Phi(lam R, mu R),
    g = gamma + ln|mu|,  Ein(x) = sum_{k>=1} x^k/(k k!),  E2(x) = sum_{k>=1} x^k/(k^2 k!),
    Phi(x, y) = sum_{j>=0, k>=1} x^j y^k / (j! k k! (j+k)),
    K = Mellin finite part:  mu < 0 (q = -mu, p = -lam):  -[(gamma + ln q)^2/2 + pi^2/12 + Li2(-p/q)]
                             mu > 0 (p = -lam > mu):      -[(gamma + ln mu)^2/2 - 5 pi^2/12 + Re Li2(p/mu)]
    f0 = c sum_i k_i T(lam_i, mu_i; R) + c LC [E1(2bR) - E1(2aR)].
    C(a, b) = c [sum_i k_i K(lam_i, mu_i) + LC ln(a/b)]  (= f0(0)).
    Log coefficient at R = 0:  f0 = -Lam ln R + (entire),  Lam = 2c [Cinh(2bR) - Cinh(2aR)],
    Cinh(x) = sum_{k>=1} x^{2k}/(2k (2k)!) = Chi(x) - ln x - gamma.
    W2 = -c sum_i k_i Phi(lam_i R, mu_i R);  a -> b limit  R/(2a) d_x Phi(2aR, -4aR),
    d_x Phi(x, y) = int_0^1 e^{xu} Ein(yu) du = e^x Ein(y)/x - (Ein(x+y) - Ein(x))/x.

All numerics are mpmath at the ambient precision; callers set it with mp.workdps.  Nothing here reads or imports
debug/.  Each formula is one small function so a fire test can replace exactly one of them.
"""
from __future__ import annotations

from fractions import Fraction
from typing import List, Tuple

import mpmath as mp
import sympy as sp


def mpq(v) -> mp.mpf:
    """Fraction / int / mpf -> mpf at the ambient precision."""
    if isinstance(v, Fraction):
        return mp.mpf(v.numerator) / v.denominator
    return mp.mpf(v)


def terms(a, b) -> List[Tuple[int, object, object]]:
    """The six (k_i, lam_i, mu_i) of label 60 (generic arithmetic: Fraction, mpf or sympy)."""
    return [(1, 2 * b, -2 * (a + b)), (-1, 2 * a, -2 * (a + b)), (-1, -2 * b, 2 * (b - a)),
            (2, -2 * b, (b - a)), (1, -2 * a, 2 * (a - b)), (-2, -2 * a, (a - b))]


def _tol() -> mp.mpf:
    return mp.mpf(10) ** (-mp.mp.dps - 8)


def Ein(x) -> mp.mpf:
    """sum_{k>=1} x^k / (k k!)  (entire)."""
    x = mp.mpf(x)
    with mp.extradps(int(abs(x) / 2.3) + 10):
        tot, t, k = mp.mpf(0), mp.mpf(1), 0
        while True:
            k += 1
            t = t * x / k                      # x^k / k!
            term = t / k
            tot += term
            if k > abs(x) + 2 and abs(term) < _tol() * max(1, abs(tot)):
                break
    return +tot


def E2(x) -> mp.mpf:
    """sum_{k>=1} x^k / (k^2 k!)."""
    x = mp.mpf(x)
    with mp.extradps(int(abs(x) / 2.3) + 10):
        tot, t, k = mp.mpf(0), mp.mpf(1), 0
        while True:
            k += 1
            t = t * x / k
            term = t / (k * k)
            tot += term
            if k > abs(x) + 2 and abs(term) < _tol() * max(1, abs(tot)):
                break
    return +tot


def Phi(x, y) -> mp.mpf:
    """sum_{j>=0, k>=1} x^j y^k / (j! k k! (j + k))  (the double series, with its 1/(j+k))."""
    x, y = mp.mpf(x), mp.mpf(y)
    with mp.extradps(int((abs(x) + abs(y)) / 2.3) + 12):
        tol = _tol()
        tot, ty, k = mp.mpf(0), mp.mpf(1), 0
        while True:
            k += 1
            ty = ty * y / k                    # y^k / k!
            mk, tx, j = mp.mpf(0), mp.mpf(1), 0   # m_k = sum_j x^j / (j! (j + k))
            while True:
                term = tx / (j + k)
                mk += term
                if j > abs(x) + 2 and abs(term) < tol * max(1, abs(mk)):
                    break
                j += 1
                tx = tx * x / j
            term = ty / k * mk
            tot += term
            if k > abs(y) + 2 and abs(term) < tol * max(1, abs(tot)):
                break
    return +tot


def Kconst(lam, mu) -> mp.mpf:
    """Mellin finite part K(lam, mu) of int_0^oo s^{eps-1} e^{lam s} Ei(mu s) ds at eps = 0 (lam + mu < 0)."""
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
    return Kconst(lam, mu) - g * L - L ** 2 / 2 - g * ein - L * ein + E2(lam * R) - Phi(lam * R, mu * R)


def T_quad(lam, mu, R):
    """Direct quadrature of int_R^oo e^{lam s} Ei(mu s) ds/s (Ei = principal value for mu > 0)."""
    lam, mu, R = mp.mpf(lam), mp.mpf(mu), mp.mpf(R)
    pts = [R, R + mp.mpf(1) / 2, R + 2, R + 6, R + 20, R + 60, mp.inf]
    return mp.quad(lambda s: mp.exp(lam * s) * mp.ei(mu * s) / s, pts)


def consts(a, b):
    a, b = mpq(a), mpq(b)
    return a, b, 1 / (2 * (a * a - b * b)), mp.log((a + b) / (b - a))


def f0_closed(a, b, R):
    a, b, c, LC = consts(a, b)
    R = mp.mpf(R)
    s = sum(k * T_closed(lam, mu, R) for k, lam, mu in terms(a, b))
    s += LC * (mp.e1(2 * b * R) - mp.e1(2 * a * R))
    return c * s


def C_ab(a, b):
    """f0(0) as the six-term Li2 constant (should vanish)."""
    a, b, c, LC = consts(a, b)
    return c * (sum(k * Kconst(lam, mu) for k, lam, mu in terms(a, b)) + LC * mp.log(a / b))


def Cinh(x):
    """sum_{k>=1} x^{2k} / (2k (2k)!)."""
    x = mp.mpf(x)
    if x == 0:
        return mp.mpf(0)
    return mp.chi(x) - mp.log(x) - mp.euler


def Lam(a, b, R):
    a, b, c, _ = consts(a, b)
    R = mp.mpf(R)
    return 2 * c * (Cinh(2 * b * R) - Cinh(2 * a * R))


def Lam_sym(a, b, R):
    """Same Lam as a sympy expression (Cinh = Chi - log - EulerGamma); a, b, R sympy symbols."""
    c = 1 / (2 * (a ** 2 - b ** 2))

    def cinh(x):
        return sp.Chi(x) - sp.log(x) - sp.EulerGamma
    return 2 * c * (cinh(2 * b * R) - cinh(2 * a * R))


def W2(a, b, R):
    a, b, c, _ = consts(a, b)
    R = mp.mpf(R)
    return -c * sum(k * Phi(lam * R, mu * R) for k, lam, mu in terms(a, b))


def dxPhi(x, y):
    x, y = mp.mpf(x), mp.mpf(y)
    return mp.exp(x) * Ein(y) / x - (Ein(x + y) - Ein(x)) / x


def W2_limit(a, R):
    a, R = mpq(a), mp.mpf(R)
    return R / (2 * a) * dxPhi(2 * a * R, -4 * a * R)


def fm1_numeric(a, b, s):
    """F^(-1)(s) from TERMS (mpmath)."""
    a, b, c, LC = consts(a, b)
    s = mp.mpf(s)
    tot = sum(k * mp.exp(lam * s) * mp.ei(mu * s) for k, lam, mu in terms(a, b))
    return c * (tot + LC * (mp.exp(-2 * b * s) - mp.exp(-2 * a * s)))


# ------------------------------------------------------------------ Theorem A: rho' + 2b rho = rhs, rho rational
def rational_ansatz_solutions(rhs, b, R, deg: int, m0: int, m1: int):
    """Solutions rho = P(R) / (R^m0 (R - 1)^m1), deg P <= deg, of rho' + 2b rho = rhs (sympy), as a list of dicts
    (empty = no solution in the ansatz).  b is a nonzero rational."""
    cs = sp.symbols(f"p0:{deg + 1}")
    rho = sum(ci * R ** i for i, ci in enumerate(cs)) / (R ** m0 * (R - 1) ** m1)
    res = sp.together(sp.diff(rho, R) + 2 * b * rho - rhs)
    num = sp.Poly(sp.numer(res), R)
    eqs = [e for e in num.all_coeffs() if e != 0]
    return sp.solve(eqs, cs, dict=True) if eqs else [dict()]


def pole_order(expr, R, point) -> int:
    """Order of the pole of a rational function of R at `point` (0 if regular)."""
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    if num == 0:
        return 0
    d = sp.Poly(den, R)
    n = sp.Poly(num, R)
    od = on = 0
    lin = sp.Poly(R - point, R)
    while d.rem(lin).is_zero:
        d = d.quo(lin)
        od += 1
    while n.rem(lin).is_zero:
        n = n.quo(lin)
        on += 1
    return max(od - on, 0)
