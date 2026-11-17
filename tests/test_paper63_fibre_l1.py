"""Guard for record P63L1 (docs/predictions/2026-11-16_paper63_fibre_l1.md, scored GO): Paper 63 sec:dictionary,
"[SYMBOLIC + MEASURED] One l > 0 shape reaches its rung".

Object: chain of three densities on the slice (t; 1, 1); s-type root rho_t = t^3 K0(t r)/(2 pi^2); two p_z leaves
rho_p = d/dz [K0(r)/(2 pi^2)] = f_p cos(theta), f_p = -K1/(2 pi^2); Coulomb kernels; coincident centres.

Claims guarded:
  (a) the p_z leaf potential is g cos(theta), g = (2/pi)[K1 - chi/r^2], chi = int_0^r s K1 -- the z-derivative of
      the s-leaf potential (2/pi) chi/r, and Laplace's L = 1 multipole formula;
  (b) E = g + 1/r^2 ~ sqrt(2/pi) r^{-1/2} e^{-r} (1 + 11/(8r) + 97/(128 r^2) - 463/(1024 r^3) + ...), the
      r^{1/2} e^{-r} terms of the multipole pieces cancelling;
  (c) from that tail alone: eps^{3/2} S(eps) -> A = 1/(3 pi^2), eps^{3/2} S = A (1 + kappa eps) + O(eps^{3/2}),
      kappa = 47/8 (s-chain toy: 1/pi^2, 31/8);
  (d) the measured branch S(eps) = int_1^oo I0((2 - eps) x) H_sing(x) dx, H_sing = (2/(3 pi)) x^2 E^2, has
      log-slope -3/2 and extrapolated amplitude 1/(3 pi^2);
  (e) the Gaunt average <cos^2> = 1/3 at the s root (the 1/3 relative to the s chain).

Self-contained: the engine pieces of debug/p63l1_engine.py, p63l1_run.py and p63l1_controls.py (frozen in the record)
are ported, not imported.
"""
from __future__ import annotations

import math
from fractions import Fraction
from typing import Callable, List

import mpmath as mp
import numpy as np
import pytest
import sympy as sp
from scipy.integrate import cumulative_trapezoid
from scipy.special import i0e, k0, k1, k1e, roots_laguerre

# ----------------------------------------------------------------------------------------- record constants
COS2 = Fraction(1, 3)                       # <cos^2 theta> over the sphere at the s root (L = 0 part of Y10 Y10)
A_CLOSED = 1 / (3 * math.pi ** 2)           # record: eps^{3/2} S -> 1/(3 pi^2)
KAPPA = Fraction(47, 8)                     # record: first correction
A_S_CLOSED = 1 / math.pi ** 2               # s-chain toy
KAPPA_S = Fraction(31, 8)
TAIL = [Fraction(1), Fraction(11, 8), Fraction(97, 128), Fraction(-463, 1024)]   # record's E tail coefficients
BRANCH_EXP = 1.5                            # eps^{-3/2} (SA's n = 3 exponent)

# Record P63L1 scored run (debug/data/p63l1_run.json, mp dps 30): S(eps) at the first two eps.
S_REF = {0.01: 35.79763454389506258164132, 0.005623413251903491: 82.78546290866481031140301}


def pref() -> float:
    """H = PREF r^2 g^2 with PREF = 4 pi <cos^2> / (2 pi^2) (root density t^3 K0(t r)/(2 pi^2)) = 2/(3 pi)."""
    return 4 * math.pi * float(COS2) / (2 * math.pi ** 2)


# ----------------------------------------------------------------------------------------- mpmath leaf potential
def chi_mp(r) -> mp.mpf:
    return mp.quad(lambda s: s * mp.besselk(1, s), [0, r])


def g_closed_mp(r) -> mp.mpf:
    r = mp.mpf(r)
    return 2 / mp.pi * (mp.besselk(1, r) - chi_mp(r) / r ** 2)


def phi_s_closed_mp(r) -> mp.mpf:
    """s-leaf potential (2/pi) chi(r)/r."""
    r = mp.mpf(r)
    return 2 / mp.pi * chi_mp(r) / r


def phi_s_multipole_mp(r) -> mp.mpf:
    """Laplace L = 0 potential of rho_1 = K0/(2 pi^2): 4 pi [r^-1 int_0^r s^2 rho + int_r^oo s rho]."""
    r = mp.mpf(r)
    rho = lambda s: mp.besselk(0, s) / (2 * mp.pi ** 2)
    return 4 * mp.pi * (mp.quad(lambda s: s * s * rho(s), [0, r]) / r + mp.quad(lambda s: s * rho(s), [r, mp.inf]))


def fp_mp(s) -> mp.mpf:
    return -mp.besselk(1, s) / (2 * mp.pi ** 2)


def g_multipole_mp(r) -> mp.mpf:
    """Laplace L = 1: (4 pi/3) [r^-2 int_0^r s^3 f_p + r int_r^oo f_p]."""
    r = mp.mpf(r)
    a = mp.quad(lambda s: s ** 3 * fp_mp(s), [0, r])
    b = mp.quad(fp_mp, [r, r + 10, mp.inf])
    return 4 * mp.pi / 3 * (a / r ** 2 + r * b)


# ----------------------------------------------------------------------------------------- R = int_r^oo s K1, three routes
def Rscaled_quad_mp(x) -> mp.mpf:
    """e^x R(x) by direct quadrature."""
    x = mp.mpf(x)
    return mp.quad(lambda u: (x + u) * mp.besselk(1, x + u) * mp.exp(x + u) * mp.exp(-u), [0, 1, 10, mp.inf])


def Rscaled_struve_mp(x) -> mp.mpf:
    """e^x R(x), R = x K0 + pi/2 - int_0^x K0, int_0^x K0 = (pi x/2)[K0 L_{-1} + K1 L_0] (engine route below x = 80)."""
    dps = mp.mp.dps
    with mp.workdps(dps + int(mp.mpf(x) / 2.3) + 15):
        x = mp.mpf(x)
        v = x * mp.besselk(0, x) + mp.pi / 2 - mp.pi * x / 2 * (mp.besselk(0, x) * mp.struvel(-1, x)
                                                              + mp.besselk(1, x) * mp.struvel(0, x))
        v = v * mp.exp(x)
    return +v


def _R_series_coeffs(M: int = 400) -> List[Fraction]:
    a = [Fraction(1)]
    for k in range(1, M):
        a.append(a[-1] * Fraction(4 - (2 * k - 1) ** 2, 8 * k))
    d = [Fraction(0)] * M
    for k in range(M):
        b = Fraction(1, 2) - k
        ff = Fraction(1)
        for j in range(M - k):
            d[k + j] += a[k] * ff
            ff *= (b - j)
    return d


_DSER: List[Fraction] = []
X_SERIES = 80


def Rscaled_series_mp(x) -> mp.mpf:
    """e^x R(x) by its optimally truncated asymptotic series (engine route at x >= 80)."""
    x = mp.mpf(x)
    if not _DSER:
        _DSER.extend(_R_series_coeffs())
    tol = mp.mpf(10) ** (-mp.mp.dps - 4)
    s, p, prev = mp.mpf(0), mp.mpf(1), mp.inf
    for dm in _DSER:
        term = mp.mpf(dm.numerator) / dm.denominator * p
        if abs(term) > prev:
            raise RuntimeError("asymptotic series diverged before reaching tolerance")
        s += term
        if abs(term) < tol and dm != 0:
            break
        prev = abs(term) if dm != 0 else prev
        p /= x
    return mp.sqrt(mp.pi / 2) * mp.sqrt(x) * s


def Rscaled_mp(x) -> mp.mpf:
    return Rscaled_series_mp(x) if mp.mpf(x) >= X_SERIES else Rscaled_struve_mp(x)


def Escaled_mp(x) -> mp.mpf:
    """e^x E(x) = (2/pi) [e^x K1(x) + e^x R(x)/x^2]."""
    x = mp.mpf(x)
    return 2 / mp.pi * (mp.besselk(1, x) * mp.exp(x) + Rscaled_mp(x) / x ** 2)


def h_l1_mp(x) -> mp.mpf:
    return mp.mpf(pref()) * x ** 2 * Escaled_mp(x) ** 2


def h_s_mp(x) -> mp.mpf:
    return (2 / mp.pi) ** 3 * Rscaled_mp(x) ** 2


# ----------------------------------------------------------------------------------------- float64 branch (SA's scheme)
_UL, _WL = roots_laguerre(80)


def Rscaled_f(x: np.ndarray) -> np.ndarray:
    X = x[:, None] + _UL[None, :]
    return (X * k1e(X)) @ _WL


def h_l1_f(x: np.ndarray) -> np.ndarray:
    """e^{2x} H_sing(x), H_sing = PREF x^2 E^2."""
    Es = 2 / np.pi * (k1e(x) + Rscaled_f(x) / x ** 2)
    return pref() * x ** 2 * Es ** 2


def h_s_f(x: np.ndarray) -> np.ndarray:
    return (2 / np.pi) ** 3 * Rscaled_f(x) ** 2


def branch_float(eps: float, h: Callable[[np.ndarray], np.ndarray], npts: int = 400001) -> float:
    """S(eps) = int_1^oo I0((2-eps)x) H_sing dx, log grid x = e^y, y from 0 (x = 1)."""
    y = np.linspace(0.0, np.log(1.0 + 80.0 / eps), npts)
    x = np.exp(y)
    f = i0e((2 - eps) * x) * h(x) * np.exp(-eps * x)
    return float(np.trapezoid(f * x, y))


def _gl_nodes(n: int):
    xs, ws = [], []

    def pdp(z):
        p0, p1 = mp.mpf(1), z
        for j in range(2, n + 1):
            p0, p1 = p1, ((2 * j - 1) * z * p1 - (j - 1) * p0) / j
        return p1, n * (z * p1 - p0) / (z * z - 1)

    for k in range(1, n + 1):
        x = mp.cos(mp.pi * (k - mp.mpf(1) / 4) / (n + mp.mpf(1) / 2))
        for _ in range(100):
            p, dp = pdp(x)
            dx = p / dp
            x -= dx
            if abs(dx) < mp.mpf(10) ** (-mp.mp.dps - 5):
                break
        p, dp = pdp(x)
        xs.append(x)
        ws.append(2 / ((1 - x ** 2) * dp ** 2))
    return xs, ws


def branch_mp(eps, h: Callable, ngl: int = 40) -> mp.mpf:
    """S(eps) at mp.dps, composite Gauss-Legendre on [1, 100/eps] (engine scheme)."""
    eps = mp.mpf(eps)
    xs, ws = _gl_nodes(ngl)
    edges = [mp.mpf(1)]
    X = 100 / eps
    while edges[-1] < X:
        a = edges[-1]
        edges.append(min(a + min(a, 2 / eps), X))
    tot = mp.mpf(0)
    a2 = 2 - eps
    for a, b in zip(edges[:-1], edges[1:]):
        hm, hw = (a + b) / 2, (b - a) / 2
        for xi, wi in zip(xs, ws):
            x = hm + hw * xi
            tot += wi * hw * mp.besseli(0, a2 * x) * mp.exp(-a2 * x) * h(x) * mp.exp(-eps * x)
    return tot


def fit_branch(eps: np.ndarray, S: np.ndarray):
    """log-slope, and least squares of eps^{3/2} S on {1, u^2, u^3, u^4}, u = sqrt(eps) (run driver)."""
    slope = float(np.polyfit(np.log(eps), np.log(S), 1)[0])
    u = np.sqrt(eps)
    y = eps ** BRANCH_EXP * S
    M = np.stack([u ** 0, u ** 2, u ** 3, u ** 4], 1)
    c, *_ = np.linalg.lstsq(M, y, rcond=None)
    return slope, c


# ----------------------------------------------------------------------------------------- symbolic tail (C3)
_x = sp.symbols("x", positive=True)
NT = 7


def hankel(nu: int, n: int = NT):
    """K_nu(x) e^x sqrt(2x/pi) = sum a_k x^-k."""
    return [sp.Rational(sp.prod([4 * nu ** 2 - (2 * j - 1) ** 2 for j in range(1, k + 1)]), sp.factorial(k) * 8 ** k)
            for k in range(n)]


def tail_int(coeffs, a, n: int = NT):
    """int_r^oo s^(a - k) e^-s sum_k c_k ds ~ r^a e^-r sum_m d_m r^-m."""
    d = [sp.Integer(0)] * n
    for k, ck in enumerate(coeffs):
        b = a - k
        for j in range(n - k):
            d[k + j] += ck * sp.ff(b, j)
    return d


def _ser(cs):
    return sum(c * _x ** -k for k, c in enumerate(cs))


def E_tail_coeffs() -> list:
    """E ~ sqrt(2/pi) r^-1/2 e^-r sum e_k r^-k, from E = (2/pi)[K1 + R/r^2], R = int_r^oo s K1."""
    Rser = tail_int(hankel(1), sp.Rational(1, 2))
    Eser = sp.expand(_ser(hankel(1)) + _x ** -1 * _ser(Rser))
    return [Eser.coeff(_x, -k) for k in range(NT)]


def multipole_tail_pieces():
    """Exponential part of the Laplace L = 1 formula, times sqrt(2/pi)^-1 r^{1/2} e^r:
    piece a: (4pi/3) r^-2 int_r^oo s^3 f_p;  piece b: (4pi/3) r int_r^oo f_p = (4pi/3)(-1/(2pi^2)) r K0."""
    P3 = tail_int(hankel(1), sp.Rational(5, 2))
    piece_a = sp.expand(-_x ** -2 * _x ** sp.Rational(5, 2) * _ser(P3))
    piece_b = sp.expand(_x * _x ** sp.Rational(-1, 2) * _ser(hankel(0)))
    c = -sp.Rational(2, 3) / sp.pi * sp.sqrt(sp.pi / 2)
    return sp.expand(c * piece_a * sp.sqrt(_x)), sp.expand(c * piece_b * sp.sqrt(_x))


_eps = sp.symbols("epsilon", positive=True)


def amplitude(hser_coeffs, pre):
    """H_sing(x) ~ pre * x * e^-2x * sum h_k x^-k; I0(y) e^-y sqrt(2 pi y) ~ sum b_k y^-k, y = (2 - eps) x.
    Singular part of int_1^oo I0((2-eps)x) H_sing = sum_k C_k(eps) Gamma(3/2 - k) eps^(k - 3/2)."""
    b = [sp.Rational(sp.prod([(2 * j - 1) ** 2 for j in range(1, k + 1)]), sp.factorial(k) * 8 ** k) for k in range(4)]
    a2 = 2 - _eps
    total = 0
    for k in range(3):
        Ck = 0
        for i in range(k + 1):
            Ck += hser_coeffs[i] * b[k - i] * a2 ** (-(k - i))
        Ck *= pre / sp.sqrt(2 * sp.pi * a2)
        total += Ck * sp.gamma(sp.Rational(3, 2) - k) * _eps ** (k - sp.Rational(3, 2))
    s = sp.series(sp.expand(total * _eps ** sp.Rational(3, 2)), _eps, 0, 2).removeO()
    A = sp.simplify(s.subs(_eps, 0))
    return A, sp.simplify(s.coeff(_eps, 1) / A)


# ========================================================================================= (a) the leaf potential
R_PTS = ("0.5", "2", "7")


@pytest.mark.parametrize("r", R_PTS)
def test_g_is_z_derivative_of_s_potential(r):
    """g(r) = d/dr [(2/pi) chi/r] (on the z axis, d/dz of the s-leaf potential), to 1e-15 by mpmath numerical
    differentiation; and (2/pi) chi/r is the Laplace L = 0 potential of rho_1 = K0/(2 pi^2).  Rejects a K1 sign
    flip in g, chi replaced by pi/2, a wrong 2/pi, or a chi/r^2 -> chi/r slip."""
    with mp.workdps(30):
        d = mp.diff(phi_s_closed_mp, mp.mpf(r))
        g = g_closed_mp(r)
        assert abs(d / g - 1) < 1e-15, (r, d, g)
        p0, p1 = phi_s_closed_mp(r), phi_s_multipole_mp(r)
        assert abs(p0 / p1 - 1) < 1e-20, (r, p0, p1)


@pytest.mark.parametrize("r", R_PTS)
def test_g_matches_laplace_L1(r):
    """g(r) equals Laplace's L = 1 multipole formula for f_p = -K1/(2 pi^2) to 1e-20 (record C2: <= 1.1e-31 at dps
    30).  Independent of chi: rejects the same g errors as the derivative test, and a wrong f_p normalisation."""
    with mp.workdps(30):
        g, m = g_closed_mp(r), g_multipole_mp(r)
        assert abs(g / m - 1) < 1e-20, (r, g, m)


@pytest.mark.parametrize("x", ["1", "3.3", "50.3", "79.9"])
def test_R_struve_route_vs_quadrature(x):
    """e^x R(x), R = int_x^oo s K1 = pi/2 - chi: the modified-Struve route (mp branch integrand below x = 80) against
    direct quadrature to 1e-25, and against pi/2 - chi at x = 1, 3.3.  Rejects a sign flip of the pi/2 or of the
    Struve term, or a dropped precision boost."""
    with mp.workdps(30):
        s, q = Rscaled_struve_mp(x), Rscaled_quad_mp(x)
        assert abs(s / q - 1) < 1e-25, (x, s, q)
        if mp.mpf(x) < 4:
            c = (mp.pi / 2 - chi_mp(x)) * mp.exp(mp.mpf(x))
            assert abs(s / c - 1) < 1e-25, (x, s, c)


def test_R_routes_agree_at_switch_and_float():
    """At the x = 80 switch the asymptotic series and the Struve route agree to 1e-28; the float64 Laguerre route
    (float branch integrand) agrees with the mp routes to 1e-13 at x = 1, 7, 80, 300.  Rejects a broken series
    coefficient (a_k or the falling factorial), or a mis-scaled Laguerre route."""
    with mp.workdps(30):
        a, b = Rscaled_series_mp(80), Rscaled_struve_mp(80)
        assert abs(a / b - 1) < 1e-28, (a, b)
        xs = np.array([1.0, 7.0, 80.0, 300.0])
        f = Rscaled_f(xs)
        for xi, fi in zip(xs, f):
            ref = Rscaled_mp(xi)
            assert abs(fi / float(ref) - 1) < 1e-13, (xi, fi, ref)


# ========================================================================================= (b) the tail
def test_tail_coefficients_symbolic():
    """Derived (sympy) from K1's Hankel series and R = int s K1: E ~ sqrt(2/pi) r^-1/2 e^-r (1 + 11/(8r) +
    97/(128r^2) - 463/(1024r^3) + ...), equal to the record's literals; and in the multipole route each piece has a
    nonzero r^{1/2} e^{-r} term that cancels in their sum, leaving sqrt(2/pi) r^{-1/2} e^{-r}.  Rejects a changed
    literal (11/8 -> 7/8), a wrong Hankel coefficient, or an L = 1 formula without the cancellation."""
    e = E_tail_coeffs()
    assert [sp.Rational(c.numerator, c.denominator) for c in TAIL] == e[:4], e[:4]
    pa, pb = multipole_tail_pieces()
    la, lb = sp.limit(pa / _x, _x, sp.oo), sp.limit(pb / _x, _x, sp.oo)
    assert la != 0 and sp.simplify(la + lb) == 0, (la, lb)
    lead = sp.limit(sp.expand(pa + pb), _x, sp.oo)
    assert sp.simplify(lead - sp.sqrt(2 / sp.pi)) == 0, lead


def _E_over_lead(r: int) -> mp.mpf:
    """E / (sqrt(2/pi) r^-1/2 e^-r) through the Laplace L = 1 formula.  With int_0^oo s^3 K1 = 3 pi/2 the r^-2 part of
    g is exactly -1/r^2, so E = (4 pi/3)(-1/(2 pi^2)) [-r^-2 int_r^oo s^3 K1 + r int_r^oo K1]: the two multipole
    pieces, each O(r^{1/2} e^{-r}), whose leading terms cancel.  int_r^oo K1 = K0(r); the other by quadrature
    (scaled by e^r)."""
    rr = mp.mpf(r)
    a = mp.quad(lambda u: (rr + u) ** 3 * mp.besselk(1, rr + u) * mp.exp(rr + u) * mp.exp(-u), [0, 5, 20, 80])
    b = mp.besselk(0, rr) * mp.exp(rr)
    Es = 4 * mp.pi / 3 * (-1 / (2 * mp.pi ** 2)) * (-a / rr ** 2 + rr * b)
    return Es / (mp.sqrt(2 / mp.pi) * rr ** mp.mpf(-0.5))


def test_tail_series_against_full_E():
    """The full E = g + 1/r^2, through the multipole formula at dps 25 (exercising the cancellation), against the
    tail at r = 30, 40, 50.  With the record's 4 literal terms the relative residual is O(r^-4) (log-slope in
    [-5, -3]) and below 2e-6 at r = 30; with the 7 derived terms it falls as r^-7 (slope in [-8, -6], record C3:
    1.4e-9, 2.0e-10, 4.4e-11).  Rejects a changed 11/8 or 97/128 (residual O(r^-1) or O(r^-2)), a wrong lead
    sqrt(2/pi), and an E carrying r^{1/2} e^{-r}."""
    rs = (30, 40, 50)
    e7 = [mp.mpf(sp.Rational(c).p) / sp.Rational(c).q for c in E_tail_coeffs()]
    e4 = [mp.mpf(c.numerator) / c.denominator for c in TAIL]
    res4, res7 = [], []
    with mp.workdps(25):
        for r in rs:
            q = _E_over_lead(r)
            res4.append(abs(float(q / sum(c * mp.mpf(r) ** -k for k, c in enumerate(e4)) - 1)))
            res7.append(abs(float(q / sum(c * mp.mpf(r) ** -k for k, c in enumerate(e7)) - 1)))
    lr = np.log(rs)
    s4 = np.polyfit(lr, np.log(res4), 1)[0]
    s7 = np.polyfit(lr, np.log(res7), 1)[0]
    assert res4[0] < 2e-6 and -5 < s4 < -3, (res4, s4)
    assert res7[0] < 3e-9 and -8 < s7 < -6, (res7, s7)


# ========================================================================================= (c) closed forms
def test_amplitude_and_kappa_from_tail():
    """From the tail coefficients alone (sympy): H_sing = (2/(3pi)) x^2 E^2 = (2/(3pi))(2/pi) x e^{-2x} Eser^2 gives
    A = 1/(3 pi^2), kappa = 47/8; the s chain (H_sing = (2/pi)^3 R^2) gives 1/pi^2, 31/8.  Both compared with the
    module constants.  Rejects A = 1/pi^2 for the p chain, kappa = 31/8, <cos^2> = 1/2 (via PREF), and a changed
    tail literal."""
    ec = [sp.Rational(c.numerator, c.denominator) for c in TAIL]
    Eser = _ser(ec)
    h = sp.series(sp.expand(Eser ** 2), _x, sp.oo, 4).removeO()
    hc = [sp.expand(h).coeff(_x, -k) for k in range(3)]
    p = sp.nsimplify(4 * sp.pi * sp.Rational(COS2.numerator, COS2.denominator) / (2 * sp.pi ** 2))
    A, kap = amplitude(hc, p * 2 / sp.pi)
    assert sp.simplify(A - sp.Rational(1, 3) / sp.pi ** 2) == 0, A
    assert abs(float(A) / A_CLOSED - 1) < 1e-15, (A, A_CLOSED)
    assert kap == sp.Rational(KAPPA.numerator, KAPPA.denominator) == sp.Rational(47, 8), kap
    Rser = tail_int(hankel(1), sp.Rational(1, 2))
    hs = sp.series(sp.expand(_ser(Rser) ** 2), _x, sp.oo, 4).removeO()
    hsc = [sp.expand(hs).coeff(_x, -k) for k in range(3)]
    As, ks = amplitude(hsc, (2 / sp.pi) ** 3 * sp.pi / 2)
    assert sp.simplify(As - 1 / sp.pi ** 2) == 0 and abs(float(As) / A_S_CLOSED - 1) < 1e-15, As
    assert ks == sp.Rational(KAPPA_S.numerator, KAPPA_S.denominator) == sp.Rational(31, 8), ks
    assert sp.simplify(A / As - sp.Rational(1, 3)) == 0


def test_fibre_branch_prefactor():
    """Im FIBRE(-(2 - eps) + i0) = pi (2 - eps)^3 S + regular: from t^3 K0(t r) with K0(-y + i0) = K0(y) - i pi I0(y),
    checked numerically on mpmath's K0 at y = 1.7, and the singular amplitude 8 pi A with first correction
    kappa - 3/2 = 35/8.  Rejects a dropped (2 - eps)^3, a wrong sign of the i pi I0 jump, or t^2 for t^3."""
    y = mp.mpf("1.7")
    with mp.workdps(30):
        jump = mp.besselk(0, -y + mp.mpf(10) ** -25 * 1j) - mp.besselk(0, y)
        assert abs(jump - (-1j * mp.pi * mp.besseli(0, y))) < 1e-20, jump
    t = -(2 - _eps)
    im_fac = sp.simplify(t ** 3 * (-sp.pi))          # Im [t^3 K0(t r)] per unit I0 = pi (2 - eps)^3
    assert sp.simplify(im_fac - branch_prefactor_sym()) == 0
    kap = sp.Rational(KAPPA.numerator, KAPPA.denominator)
    s = sp.series(branch_prefactor_sym() * (1 + kap * _eps), _eps, 0, 2).removeO()
    assert s.subs(_eps, 0) == 8 * sp.pi and sp.simplify(s.coeff(_eps, 1) / (8 * sp.pi)) == sp.Rational(35, 8)


def branch_prefactor_sym():
    return sp.pi * (2 - _eps) ** 3


# ========================================================================================= (d) measured branch, reduced
EPS_FAST = np.logspace(-2, -3, 5)


@pytest.fixture(scope="module")
def scan_l1():
    return np.array([branch_float(e, h_l1_f) for e in EPS_FAST])


def test_branch_anchor_values(scan_l1):
    """float64 S at eps = 1e-2 and 10^-2.25 equals the record's mp30 values to 1e-11 (record FB0: 9.9e-14).  The only
    test sensitive to the lower limit x = 1 (it moves only the u^3 term, invisible to the amplitude fit), and to any
    integrand change.  Rejects a changed lower limit, a dropped e^{-eps x}, or a changed grid end."""
    for e, s in zip(EPS_FAST[:2], scan_l1[:2]):
        ref = S_REF[min(S_REF, key=lambda k: abs(k - e))]
        assert abs(s / ref - 1) < 1e-11, (e, s, ref)


def test_branch_slope_and_amplitude_fast(scan_l1):
    """eps = 1e-2 .. 1e-3 (5 points): log-slope within 0.05 of -3/2; the {1, u^2, u^3, u^4} fit of eps^{3/2} S gives
    A within 1e-4 of 1/(3 pi^2) and u^2 coefficient / (A kappa) in [0.95, 1.05].  Rejects a missing x^2 or
    <cos^2> = 1/2 in H_sing, a fit model in eps^1, A = 1/pi^2, kappa = 31/8, and a branch that does not reach -2."""
    slope, c = fit_branch(EPS_FAST, scan_l1)
    assert abs(slope + 1.5) < 0.05, slope
    assert abs(c[0] / A_CLOSED - 1) < 1e-4, (c[0], A_CLOSED)
    r = c[1] / (A_CLOSED * float(KAPPA))
    assert 0.95 < r < 1.05, r


def test_branch_toy_s_chain_fast():
    """The s-chain toy (H_sing = (2/pi)^3 R^2): same reduced scan, slope within 0.05 of -3/2, amplitude within 1e-4
    of 1/pi^2, u^2 / (A kappa_s) in [0.95, 1.05].  Checks the instrument on a known answer; the p-chain amplitude is
    then 1/3 of it.  Rejects a broken Laguerre R, I0 scaling, or fit."""
    S = np.array([branch_float(e, h_s_f) for e in EPS_FAST])
    slope, c = fit_branch(EPS_FAST, S)
    assert abs(slope + 1.5) < 0.05, slope
    assert abs(c[0] / A_S_CLOSED - 1) < 1e-4, (c[0], A_S_CLOSED)
    assert 0.95 < c[1] / (A_S_CLOSED * float(KAPPA_S)) < 1.05, c


# ========================================================================================= (e) Gaunt average
def test_cos2_average_is_one_third():
    """<cos^2 theta> over the unit sphere = (1/2) int_{-1}^{1} mu^2 dmu = 1/3, equal to COS2; and the L = 0
    component of Y10 Y10 is 1/sqrt(4 pi) (Gaunt), so the s root keeps exactly the 1/3.  Rejects COS2 = 1/2."""
    avg = mp.quad(lambda m: m * m, [-1, 1]) / 2
    assert abs(avg - mp.mpf(COS2.numerator) / COS2.denominator) < 1e-25, avg
    from sympy.physics.wigner import gaunt
    g = gaunt(1, 1, 0, 0, 0, 0)             # int Y10 Y10 Y00; cos^2 = (4 pi/3) Y10^2, Y00 = 1/sqrt(4 pi)
    avg_gaunt = sp.Rational(1, 3) * sp.sqrt(4 * sp.pi) * g
    assert sp.simplify(avg_gaunt - sp.Rational(COS2.numerator, COS2.denominator)) == 0, g


# ========================================================================================= extra: the reduction at real t
class _Grid:
    def __init__(self, npts: int = 200001, rmin: float = 1e-10, rmax: float = 90.0) -> None:
        self.U = np.linspace(np.log(rmin), np.log(rmax), npts)
        self.R = np.exp(self.U)

    def cum0(self, f):
        return cumulative_trapezoid(f * self.R, self.U, initial=0.0)

    def tail(self, f):
        c = self.cum0(f)
        return c[-1] - c

    def integral(self, f) -> float:
        return float(np.trapezoid(f * self.R, self.U))

    def rho_s(self, a):
        return a ** 3 * k0(a * self.R) / (2 * np.pi ** 2)

    def f_p(self):
        return -k1(self.R) / (2 * np.pi ** 2)

    def phi1(self, W):
        R = self.R
        return 4 * np.pi / 3 * (self.cum0(W * R ** 3) / R ** 2 + R * self.tail(W))


@pytest.mark.parametrize("t", [0.5, 1.6])
def test_reduction_three_routes(t):
    """FIBRE(t) rooted at the s root (multipole recursion), rooted at a leaf (an L = 1 density re-expanded), and the
    reduced form t^3 int K0(t r) H with H = PREF r^2 g^2, g closed: root vs leaf to 1e-10, closed vs recursion to
    1e-6 (record C4 at t = 1/2: 1.9e-12, 9e-9; measured 1.6e-7 at t = 1.6, grid error).  Rejects a wrong t^3, a wrong g in H, or a broken L = 1 recursion."""
    G = _Grid()
    ang = 4 * np.pi * float(COS2)
    Vp = G.phi1(G.f_p())
    root = ang * G.integral(G.rho_s(t) * Vp ** 2 * G.R ** 2)
    leaf = ang * G.integral(G.f_p() * G.phi1(G.rho_s(t) * Vp) * G.R ** 2)
    chi = G.cum0(G.R * k1(G.R))
    g = (2 / np.pi) * (k1(G.R) - chi / G.R ** 2)
    closed = t ** 3 * G.integral(k0(t * G.R) * pref() * G.R ** 2 * g ** 2)
    assert abs(root / leaf - 1) < 1e-10, (root, leaf)
    assert abs(closed / root - 1) < 1e-6, (closed, root)


# ========================================================================================= slow: the scored scan
EPS_FULL = np.logspace(-2, -4, 9)


@pytest.mark.slow
def test_full_scan_FB0_FB1_FB2_slow():
    """The record's 9-point scan, float64 and mp30 (engine schemes).  FB0: they agree to 1e-10 relative at every eps.
    FB1: last-5 log-slope in [-1.55, -1.45] (record -1.49775).  FB2: A_fit over the smallest 6 within 1e-5 of
    1/(3 pi^2) (record 1.2e-8).  FB3 (reported, not gating in the record; checked here): u^2 / (A kappa) in
    [0.95, 1.05] (record 1.00005).  Rejects every error the fast scan rejects at 10x tighter bands, plus an mp-path
    error (Struve route, series switch) through FB0."""
    Sf = np.array([branch_float(e, h_l1_f) for e in EPS_FULL])
    with mp.workdps(30):
        Sm = np.array([float(branch_mp(mp.mpf(e), h_l1_mp)) for e in EPS_FULL])
    assert np.max(np.abs(Sf / Sm - 1)) < 1e-10, Sf / Sm - 1
    slope5 = float(np.polyfit(np.log(EPS_FULL[-5:]), np.log(Sm[-5:]), 1)[0])
    assert -1.55 <= slope5 <= -1.45, slope5
    _, c = fit_branch(EPS_FULL[3:], Sm[3:])
    assert abs(c[0] / A_CLOSED - 1) <= 1e-5, (c[0], A_CLOSED)
    assert 0.95 <= c[1] / (A_CLOSED * float(KAPPA)) <= 1.05, c


# TIMINGS (2026-11-16, this machine): fast tests 57 s in all (21 tests; the reduced float64 scan ~8 s per object, the
# mp tail check at r = 30..50 11 s, the mp leaf-potential checks 1-5 s each).  Slow test 136 s (--slow; 9-point scan,
# float64 + mp30).  Fast scan measured: slope -1.4779, A/closed - 1 = -1.4e-6, u^2/(A kappa) 1.0008; toy slope
# -1.4874, A/closed - 1 = 8.6e-8.
#
# FIRE-TEST PLANTS (python debug/qa/fire_test.py tests/test_paper63_fibre_l1.py --plant "OLD=>NEW" -k <sel>),
# 2026-11-16, written and fired as a separate step after the guard; 28 runs, 27 fired; the 1 that did not is
# explained below.
#     <cos^2> = 1/2                      "COS2 = Fraction(1, 3)" -> Fraction(1, 2)           -k "cos2 or amplitude_and_kappa or slope_and_amplitude_fast"   FIRED
#     K1 sign in g (mp)                  "(mp.besselk(1, r) - chi_mp..." -> "(-mp.besselk..."  -k "g_is or g_matches"   FIRED
#     K1 sign in E (float branch)        "(k1e(x) + Rscaled_f..." -> "(-k1e(x) + ..."         -k "anchor or slope_and_amplitude_fast"   FIRED
#     K1 sign in g (grid reduction)      "(k1(G.R) - chi ..." -> "(-k1(G.R) - chi ..."         -k reduction   FIRED
#     chi -> pi/2 in g                   "- chi_mp(r) / r ** 2)" -> "- (mp.pi / 2) / r ** 2)"  -k "g_is or g_matches"   FIRED
#     A = 1/pi^2 (the s value)           "1 / (3 * math.pi ** 2)" -> "1 / (1 * math.pi ** 2)"  -k "amplitude_and_kappa or slope_and_amplitude_fast"   FIRED
#     kappa = 31/8                       "KAPPA = Fraction(47, 8)" -> Fraction(31, 8)          -k "amplitude_and_kappa or slope_and_amplitude_fast or prefactor"   FIRED
#     missing (2-eps)^3 (power 2)        "sp.pi * (2 - _eps) ** 3" -> "** 2"                   -k prefactor   FIRED
#     missing x^2 in H_sing              "pref() * x ** 2 * Es ** 2" -> "pref() * x * Es ** 2" -k "anchor or slope_and_amplitude_fast"   FIRED
#     tail 11/8 -> 7/8                   "Fraction(11, 8)" -> "Fraction(7, 8)"                 -k "tail or amplitude_and_kappa"   FIRED
#     tail 97/128 -> 95/128              "Fraction(97, 128)" -> "Fraction(95, 128)"            -k "tail or amplitude_and_kappa"   FIRED
#     lower limit 1 -> 2 (float64)       "np.linspace(0.0, ..." -> "np.linspace(np.log(2.0), ..."   -k anchor   FIRED
#     lower limit 1 -> 2 (mp)            "edges = [mp.mpf(1)]" -> "[mp.mpf(2)]"                -k full_scan (slow, FB0)   FIRED
#     Struve route pi/2 sign             "+ mp.pi / 2 - mp.pi * x / 2 *" -> "- mp.pi / 2 ..."  -k struve   FIRED
#     Struve term sign                   "+ mp.pi / 2 - mp.pi" -> "+ mp.pi / 2 + mp.pi"        -k "struve or routes_agree"   FIRED
#     fit exponent eps^-1                "BRANCH_EXP = 1.5" -> 1.0                             -k "slope_and_amplitude_fast or toy"   FIRED
#     fit basis u^2 -> u^1               "[u ** 0, u ** 2, ...]" -> "[u ** 0, u ** 1, ...]"    -k "slope_and_amplitude_fast or toy"   FIRED
#     Laguerre nodes 80 -> 8             "roots_laguerre(80)" -> "roots_laguerre(8)"           -k "routes_agree or anchor"   FIRED
#     Hankel 4nu^2 - -> +                "4 * nu ** 2 - (2 * j - 1)" -> "+"                    -k tail_coefficients   FIRED
#     R series falling factorial         "ff *= (b - j)" -> "ff *= (b + j)"                    -k routes_agree   FIRED
#     t^3 -> t^2 in the reduction        "closed = t ** 3 *" -> "t ** 2 *"                     -k reduction   FIRED
#     e^{-eps x} halved                  "np.exp(-eps * x)" -> "np.exp(-eps * x / 2)"          -k anchor   FIRED
#     Laplace L = 1 4pi/3 -> 4pi/5       in g_multipole_mp                                     -k g_matches   FIRED
#     multipole piece a sign             "piece_a = sp.expand(-_x" -> "sp.expand(_x"           -k tail_coefficients   FIRED
#     L = 1 recursion r^-2 -> r^-1       in _Grid.phi1                                         -k reduction   FIRED
#     I0 argument 2 - eps -> 2           "i0e((2 - eps) * x)" -> "i0e((2 - 0 * eps) * x)"      -k "anchor or slope_and_amplitude_fast"   FIRED
#     A off by 1e-4 (FB2 band check)     "1 / (3 * math.pi ** 2)" -> "1 / (3.0003 * math.pi ** 2)"   -k full_scan (slow)   FIRED
#   Did not fire, and why (not a guard defect):
#     lower limit 1 -> 2 with -k slope_and_amplitude_fast alone: a lower limit moves S by a function holomorphic at
#       eps = 0, i.e. only the u^3 term of eps^{3/2} S; the amplitude, slope and u^2 fit are blind to it by
#       construction.  The same plant fires test_branch_anchor_values (record mp30 values) and, in the mp path, FB0.
