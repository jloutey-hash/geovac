"""Guard for record P63L2 (docs/predictions/2026-11-17_paper63_attainment_l2.md, scored GO): Paper 63 sec:dictionary,
"[SYMBOLIC + MEASURED] Three l > 0 shapes reach the n = 3 branch at t = -2" -- the second chain.

Object: chain of three on the slice (t; 1, 1), fingerprint weights, Coulomb kernels.  Leaf p_z (f_p cos theta,
f_p = -K1/(2 pi^2)), leaf potential g cos theta with g = -1/r^2 + E (record P63L1); middle p_z, so the middle source
is f_p g cos^2 = f_p g [1/3 + (2/3) P2]; root of scale t:
  B1  d root  h_2(d) t^3 K0(t r)/(2 pi^2) = t^5 K2(t r) P2/(2 pi^2)            (L = 2 channel only)
  B2  mixed   (lap - (5/2) h_2(d)) t^3 K0/(2 pi^2) = t^5 [(K0 - K1/(t r)) - (5/2) K2 P2]/(2 pi^2)
Branch at t = -(2 - eps) + i0, y = (2 - eps) x, Phi_Lx = e^{-2x} part of the channel-L potential of f_p E:
  S1(eps) = int_1^oo I2(y) x^2 (2/3) Phi2x dx,      S2(eps) = int_1^oo x^2 [(I0 + I1/y) Phi0x - I2 Phi2x] dx.

Claims guarded:
  (a) angular bookkeeping: cos^2 = 1/3 + (2/3) P2; int cos^2 (1 - (5/2) P2) dOmega = 0, so the B2 channel weights
      (1 and 1, in S2's normalisation) cancel at leading order; the two root identities h_2 K0 = K2 P2,
      lap K0 = K0 - K1/r, and the root kernels' jumps at t = -y + i0 (I0 + I1/y, I2, same sign);
  (b) the local rule (e^{-kr} part of a channel-L potential = formal solution of (theta-L)(theta+L+1)Phi = -4pi r^2 s)
      reproduces the P63L1 controls (E tail 1, 11/8, 97/128, -463/1024; A = 1/(3 pi^2), kappa = 47/8; s chain
      1/pi^2, 31/8), and its 30-term series for e^{2x} Phi_L (L = 0, 2) equals a direct Laguerre quadrature of the
      multipole integral at x = 32, 40, 60;
  (c) closed forms derived here (sympy, from the tails, not read from any JSON): A1 = 1/(12 pi), kappa1 = 15/8;
      the eps^{-3/2} coefficient of B2 is exactly 0, A2 = 3/(8 pi), kappa2 = -3/8; and A1, A2 again by a hand
      leading-order route (leads of f_p and E, -4 pi s0/k^2, Hankel 1 - (4nu^2-1)/(8y), Gamma), independent of the
      series machinery;
  (d) fast: a fresh float64 evaluation of S1, S2 (direct integrand from x = 1, no series) at 7 eps in [1e-3, 1e-2]:
      anchors vs the record's referee values, eps^{3/2} S1/A1 - 1 - kappa1 eps = O(eps^{3/2}),
      eps^{1/2} S2/A2 - 1 = O(eps^{1/2}), log-slopes near -3/2 and -1/2, amplitude fits;
  (e) slow: the referee (frozen mp30 node cache on [1, 32] + this file's own 30-term local-rule series beyond 32)
      at the record's 9 eps, reproducing the record's referee S; FB0 float-vs-referee <= 1e-9; FB1/FB2 and FB4/FB5
      (smallest-6 fits) within the record's bars; FB3/FB6 checked.

Self-contained: the engine pieces of debug/p63odd_engine.py, p63odd_tail.py, p63l2_run.py and p63l1_engine.py (all
frozen in the record) are ported, not imported; the only debug/ input is the frozen cache debug/data/p63l2_cache.json
(slow tier, sha256-checked, read-only).
"""
from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Callable, Dict, List, Tuple

import mpmath as mp
import numpy as np
import pytest
import sympy as sp
from scipy.integrate import simpson
from scipy.special import i0e, i1e, ive, k1e, roots_laguerre

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "debug" / "data" / "p63l2_cache.json"
CACHE_SHA = "dc45a5b7be362e7c4e4b13a989dba885ed99b5188288a9a731782e4403fc9faf"

# ----------------------------------------------------------------------------------------- record constants
LEG0, LEG2 = Fraction(1, 3), Fraction(2, 3)   # cos^2 = LEG0 + LEG2 P2
B2_COEF = Fraction(5, 2)                       # mixed root: lap - B2_COEF h_2
A1 = 1 / (12 * math.pi)                        # record: eps^{3/2} S1 -> 1/(12 pi)
KAPPA1 = Fraction(15, 8)
A2 = 3 / (8 * math.pi)                         # record: eps^{1/2} S2 -> 3/(8 pi)
KAPPA2 = Fraction(-3, 8)
EXP1 = Fraction(3, 2)                          # S1 ~ eps^{-3/2}
EXP2 = Fraction(1, 2)                          # S2 ~ eps^{-1/2} (moved by one)

# Record P63L2 referee values (debug/data/p63l2_run.json, mp dps 30), eps = logspace(-2, -4, 9).
S1_REF = ["26.92167192862705313286", "63.46327167554553409146", "149.9466791131043847254", "354.8029164220489506879",
          "840.2884692639485361951", "1991.148645971359830073", "4719.724361185652976101", "11189.47150436674618028",
          "26530.69196144299731759"]
S2_REF = ["1.133564346633575467876", "1.532424286560809416271", "2.063938698123057483348", "2.77239796524999690392",
          "3.716870787380603366149", "4.976124433825436244293", "6.655190674282537615867", "8.894124284542608514497",
          "11.87968411926879836657"]
EPS_FULL = np.logspace(-2, -4, 9)


def channel_weights() -> Tuple[sp.Expr, sp.Expr]:
    """B2 channel weights in S2's normalisation (record: both 1): the L = 0 channel carries
    int (LEG0 + LEG2 P2) dOmega, the L = 2 channel B2_COEF int (LEG0 + LEG2 P2) P2 dOmega, each over 4 pi/3."""
    mu = sp.symbols("mu")
    P2 = (3 * mu ** 2 - 1) / 2
    src = sp.Rational(LEG0.numerator, LEG0.denominator) + sp.Rational(LEG2.numerator, LEG2.denominator) * P2
    w0 = 2 * sp.pi * sp.integrate(src, (mu, -1, 1))
    w2 = sp.Rational(B2_COEF.numerator, B2_COEF.denominator) * 2 * sp.pi * sp.integrate(src * P2, (mu, -1, 1))
    return sp.simplify(w0 / (4 * sp.pi / 3)), sp.simplify(w2 / (4 * sp.pi / 3))


W0, W2 = (float(w) for w in channel_weights())
B1_SCALE = float(LEG2)                         # S1's (2/3)


# ========================================================================================= local rule (sympy, exact)
NS = 8
_eps = sp.symbols("epsilon", positive=True)


def a_coef(nu: int, n: int = NS) -> List[sp.Rational]:
    """K_nu(x) ~ sqrt(pi/(2x)) e^-x sum a_k x^-k."""
    out = [sp.Integer(1)]
    for k in range(1, n):
        out.append(out[-1] * sp.Rational(4 * nu ** 2 - (2 * k - 1) ** 2, 8 * k))
    return out


class Ser:
    """pref * e^{-k x} x^a * sum_j c_j x^-j, c_j rational."""

    def __init__(self, pref, k, a, c) -> None:
        self.pref, self.k, self.a, self.c = pref, sp.Rational(k), sp.Rational(a), list(c)

    def __mul__(self, o: "Ser") -> "Ser":
        n = min(len(self.c), len(o.c))
        c = [sum(self.c[i] * o.c[j - i] for i in range(j + 1)) for j in range(n)]
        return Ser(self.pref * o.pref, self.k + o.k, self.a + o.a, c)

    def scale(self, s) -> "Ser":
        return Ser(self.pref * s, self.k, self.a, self.c)

    def xpow(self, p) -> "Ser":
        return Ser(self.pref, self.k, self.a + p, self.c)


def K_ser(nu: int, n: int = NS) -> Ser:
    return Ser(sp.sqrt(sp.pi / 2), 1, sp.Rational(-1, 2), a_coef(nu, n))


def local_potential(s: Ser, L: int) -> Ser:
    """e^{-k r} part of the channel-L potential of source s: the formal solution of
    (theta - L)(theta + L + 1) Phi = -4 pi r^2 s, (theta^2 + theta - L(L+1)) e^{-kr} r^b
    = [k^2 r^2 - 2k(b+1) r + b(b+1) - L(L+1)] e^{-kr} r^b."""
    k, a = s.k, s.a
    phi: List[sp.Rational] = []
    for j in range(len(s.c)):
        acc = s.c[j]
        if j >= 1:
            b = a - (j - 1)
            acc += 2 * k * (b + 1) * phi[j - 1]
        if j >= 2:
            b = a - (j - 2)
            acc -= (b * (b + 1) - L * (L + 1)) * phi[j - 2]
        phi.append(acc / k ** 2)
    return Ser(-4 * sp.pi * s.pref, k, a, phi)


def kernel_I(nu: int, extra_I1_over_y: bool = False, n: int = NS) -> List[sp.Rational]:
    """I_nu(y) ~ (2 pi y)^{-1/2} e^y sum (-1)^k a_k y^-k [+ I1(y)/y]."""
    c = [(-1) ** k * x for k, x in enumerate(a_coef(nu, n))]
    if extra_I1_over_y:
        c1 = [(-1) ** k * x for k, x in enumerate(a_coef(1, n))]
        c = [c[j] + (c1[j - 1] if j >= 1 else 0) for j in range(n)]
    return c


def branch_expansion(terms, order: int = 2) -> sp.Expr:
    """Singular part of sum sgn int_1^oo kernel(y) H dx, y = (2 - eps) x, H with e^{-2x};
    int_1^oo x^q e^{-eps x} = Gamma(q + 1) eps^{-q-1} + analytic."""
    tot = 0
    for sgn, c, H in terms:
        assert H.k == 2
        for j in range(order + 2):
            ser = sp.series((2 - _eps) ** (sp.Rational(-1, 2) - j), _eps, 0, order + 2).removeO()
            for kk in range(order + 2):
                q = H.a - sp.Rational(1, 2) - j - kk
                tot += (sgn * H.pref * H.c[kk] * c[j] * (2 * sp.pi) ** sp.Rational(-1, 2) * ser
                        * sp.gamma(q + 1) * _eps ** (-q - 1))
    return sp.expand(tot)


def powers(expr) -> Dict[sp.Rational, sp.Expr]:
    d: Dict = {}
    for term in sp.Add.make_args(sp.expand(expr)):
        c, p = term.as_coeff_exponent(_eps)
        d[p] = sp.simplify(d.get(p, 0) + c)
    return {p: v for p, v in sorted(d.items(), key=lambda kv: kv[0]) if v != 0}


def f_p_ser(n: int = NS) -> Ser:
    return K_ser(1, n).scale(-1 / (2 * sp.pi ** 2))


def E_full_ser(n: int = NS) -> Ser:
    """e^{-r} part of the L = 1 potential of f_p (the leaf's E)."""
    return local_potential(f_p_ser(n), 1)


def derive_B1():
    w = (f_p_ser() * E_full_ser()).scale(sp.Rational(LEG2.numerator, LEG2.denominator))
    H1 = local_potential(w, 2).xpow(2)
    return powers(branch_expansion([(1, kernel_I(2), H1)]))


def derive_B2():
    w0, w2 = channel_weights()
    src = f_p_ser() * E_full_ser()
    H0 = local_potential(src, 0).xpow(2)
    H2 = local_potential(src, 2).xpow(2)
    return powers(branch_expansion([(w0, kernel_I(0, True), H0), (-w2, kernel_I(2), H2)], order=3))


def _q(f: Fraction) -> sp.Rational:
    return sp.Rational(f.numerator, f.denominator)


# ========================================================================================= float64 route (direct)
_UR, _WR = roots_laguerre(80)
_UL, _WLG = roots_laguerre(120)


def Rscaled_f(x: np.ndarray) -> np.ndarray:
    """e^x R(x), R = int_x^oo s K1(s) ds."""
    X = x[:, None] + _UR[None, :]
    return (X * k1e(X)) @ _WR


def what_f(u: np.ndarray) -> np.ndarray:
    """e^{2u} f_p(u) E(u), E = (2/pi)[K1 + R/u^2]."""
    Es = 2 / np.pi * (k1e(u) + Rscaled_f(u) / u ** 2)
    return -k1e(u) * Es / (2 * np.pi ** 2)


def phix_f(x: np.ndarray, L: int) -> np.ndarray:
    """e^{2x} Phi_Lx[f_p E](x) = (4 pi/(2L+1)) x int_0^oo (1+v/x)^{L+2} expm1(-(2L+1) log1p(v/x)) w(x+v) dv
    (the multipole integral (4pi/(2L+1)) [x^L int_x u^{1-L} w - x^{-L-1} int_x u^{L+2} w], e^{-2v} by Laguerre)."""
    out = np.empty_like(x)
    for i in range(0, x.size, 4000):
        xs = x[i:i + 4000]
        V = _UL[None, :] / 2
        lp = np.log1p(V / xs[:, None])
        ker = np.exp((L + 2) * lp) * np.expm1(-(2 * L + 1) * lp)
        X = xs[:, None] + V
        vals = what_f(X.ravel()).reshape(X.shape)
        out[i:i + 4000] = 4 * np.pi / (2 * L + 1) * xs * ((ker * vals) @ _WLG) / 2
    return out


class FloatRoute:
    """S1, S2 on one log grid x = e^y, y in [0, log(1 + 80/eps_min)], Simpson; Phi_0, Phi_2 computed once."""

    def __init__(self, eps_min: float, npts: int) -> None:
        self.y = np.linspace(0.0, np.log(1.0 + 80.0 / eps_min), npts)
        self.x = np.exp(self.y)
        self.p0 = phix_f(self.x, 0)
        self.p2 = phix_f(self.x, 2)

    def S1(self, eps: float) -> float:
        x = self.x
        f = ive(2, (2 - eps) * x) * x ** 2 * B1_SCALE * self.p2 * np.exp(-eps * x)
        return float(simpson(f * x, x=self.y))

    def S2(self, eps: float) -> float:
        x = self.x
        yy = (2 - eps) * x
        f = (W0 * (i0e(yy) + i1e(yy) / yy) * self.p0 - W2 * ive(2, yy) * self.p2) * x ** 2 * np.exp(-eps * x)
        return float(simpson(f * x, x=self.y))


def fit(eps: np.ndarray, S: np.ndarray, pw: float, pows: List[int]) -> Dict[int, float]:
    u = np.sqrt(eps)
    M = np.stack([u ** p for p in pows], 1)
    c, *_ = np.linalg.lstsq(M, eps ** pw * S, rcond=None)
    return dict(zip(pows, c))


def slope(eps: np.ndarray, S: np.ndarray) -> float:
    return float(np.polyfit(np.log(eps), np.log(np.abs(S)), 1)[0])


# ========================================================================================= (a) angular + roots
def test_legendre_split_and_cancelling_weight():
    """cos^2 = LEG0 + LEG2 P2 = 1/3 + (2/3) P2 identically; int cos^2 (1 - (5/2) P2) dOmega = 0 with cos^2 literal;
    the separate channel integrals 4 pi/3 and (5/2)(8 pi/15) are each nonzero; S2's channel weights are (1, 1).
    Rejects the split 1/2 + 1/2 P2, a root coefficient other than 5/2 (a non-cancelling root), and a weight that
    cancels only because both channels vanish."""
    mu = sp.symbols("mu")
    P2 = (3 * mu ** 2 - 1) / 2
    assert sp.expand(mu ** 2 - (_q(LEG0) + _q(LEG2) * P2)) == 0
    i0 = 2 * sp.pi * sp.integrate(mu ** 2, (mu, -1, 1))
    i2 = 2 * sp.pi * sp.integrate(mu ** 2 * P2, (mu, -1, 1))
    assert i0 == 4 * sp.pi / 3 and i2 == 8 * sp.pi / 15
    assert sp.simplify(i0 - _q(B2_COEF) * i2) == 0
    assert channel_weights() == (1, 1)


def test_root_identities():
    """For f = K0(r): h_2(d) f = (f'' - f'/r) P2 with f'' - f'/r = K2, and lap f = f'' + 2 f'/r = K0 - K1/r (so the
    B1 root is K2 P2 and the B2 radial parts are (K0 - K1/r) - (5/2) K2 P2).  Checked by sympy's Bessel derivatives
    and numerically at r = 0.7, 3.  Rejects K1/r with the wrong sign or K2 replaced by K1."""
    r = sp.symbols("r", positive=True)
    f = sp.besselk(0, r)
    d1, d2 = sp.diff(f, r), sp.diff(f, r, 2)
    for rv in (sp.Rational(7, 10), sp.Integer(3)):
        hp = (d2 - d1 / r).subs(r, rv)
        lp = (d2 + 2 * d1 / r).subs(r, rv)
        assert abs(sp.N(hp - sp.besselk(2, rv), 30)) < 1e-25
        assert abs(sp.N(lp - (sp.besselk(0, rv) - sp.besselk(1, rv) / rv), 30)) < 1e-25


def test_root_kernel_jumps():
    """At z = -y + i0: Im[K0(z) - K1(z)/z] = -pi (I0(y) + I1(y)/y) and Im K2(z) = -pi I2(y) (same sign, so the B2
    kernels are (I0 + I1/y) and I2 with relative sign -1 from the -(5/2) P2).  Rejects dropping I1/y or flipping
    its sign."""
    with mp.workdps(30):
        for y in (mp.mpf("0.9"), mp.mpf("2.3")):
            z = -y + mp.mpf(10) ** -25 * 1j
            j0 = mp.im(mp.besselk(0, z) - mp.besselk(1, z) / z)
            j2 = mp.im(mp.besselk(2, z))
            assert abs(j0 + mp.pi * (mp.besseli(0, y) + mp.besseli(1, y) / y)) < 1e-18, j0
            assert abs(j2 + mp.pi * mp.besseli(2, y)) < 1e-18, j2


# ========================================================================================= (b) local rule controls
def test_local_rule_reproduces_P63L1_tail():
    """E = e^{-r} part of the L = 1 potential of f_p: prefactor sqrt(2/pi) and tail 1, 11/8, 97/128, -463/1024
    (record P63L1, derived there from K1's Hankel series and R = int s K1).  Rejects a sign or k^2 slip in the
    recurrence, and L entering at leading order."""
    E = E_full_ser()
    c0 = E.c[0]
    assert sp.simplify(E.pref * c0 - sp.sqrt(2 / sp.pi)) == 0
    assert [x / c0 for x in E.c[:4]] == [1, sp.Rational(11, 8), sp.Rational(97, 128), sp.Rational(-463, 1024)]


def _R_series(n: int = NS) -> List[sp.Rational]:
    """e^x R(x) ~ sqrt(pi/2) x^{1/2} sum d_m x^-m (P63L1 engine recurrence)."""
    a = a_coef(1, n)
    d = [sp.Integer(0)] * n
    for k in range(n):
        b = sp.Rational(1, 2) - k
        ff = sp.Integer(1)
        for j in range(n - k):
            d[k + j] += a[k] * ff
            ff *= (b - j)
    return d


def test_local_rule_reproduces_P63L1_and_s_chain_constants():
    """Through the same branch_expansion: P63L1's H = (2/(3pi)) x^2 E^2 with kernel I0 gives 1/(3 pi^2), 47/8; the s
    chain H = (2/pi)^3 R^2 gives 1/pi^2, 31/8.  Rejects a broken Gamma bookkeeping or (2 - eps) expansion."""
    E = E_full_ser()
    H = (E * E).xpow(2).scale(2 / (3 * sp.pi))
    P = powers(branch_expansion([(1, kernel_I(0), H)]))
    A = P[-sp.Rational(3, 2)]
    assert sp.simplify(A - 1 / (3 * sp.pi ** 2)) == 0 and sp.simplify(P[-sp.Rational(1, 2)] / A) == sp.Rational(47, 8)
    R = Ser(sp.sqrt(sp.pi / 2), 1, sp.Rational(1, 2), _R_series())
    Hs = (R * R).scale((2 / sp.pi) ** 3)
    P = powers(branch_expansion([(1, kernel_I(0), Hs)]))
    A = P[-sp.Rational(3, 2)]
    assert sp.simplify(A - 1 / sp.pi ** 2) == 0 and sp.simplify(P[-sp.Rational(1, 2)] / A) == sp.Rational(31, 8)


NSER = 30
_SER_CACHE: Dict[int, Tuple[float, List[mp.mpf]]] = {}


def phi_series_coeffs(L: int) -> Tuple[mp.mpf, List[mp.mpf]]:
    """e^{2x} Phi_Lx[f_p E] ~ x^a sum c_j x^-j, NSER terms (exact rationals, then mp at dps 40)."""
    if L not in _SER_CACHE:
        src = f_p_ser(NSER) * E_full_ser(NSER)
        Ph = local_potential(src, L)
        assert Ph.k == 2
        with mp.workdps(40):
            pref = mp.mpf(sp.N(Ph.pref, 45))
            _SER_CACHE[L] = (mp.mpf(sp.N(Ph.a, 45)), [pref * mp.mpf(c.p) / c.q for c in Ph.c])
    return _SER_CACHE[L]


def phi_series(x, L: int) -> mp.mpf:
    a, c = phi_series_coeffs(L)
    x = mp.mpf(x)
    return x ** a * mp.fsum(cj * x ** (-j) for j, cj in enumerate(c))


@pytest.mark.parametrize("L", [0, 2])
def test_series_vs_direct_quadrature(L):
    """The local rule's 30-term series of e^{2x} Phi_Lx (L = 0, 2) against the direct multipole quadrature (float,
    Laguerre) at x = 32, 40, 60: within 3e-14 (record D1: 6.8e-15 to 7.3e-15, the float route's known bias).  The
    only check of the local rule beyond its leading term that does not go through a branch fit; it pins every tail
    coefficient that kappa1, kappa2 are built from.  Rejects an L-dependence slip in the recurrence, a 2k(b+1) or
    b(b+1) error, or a wrong channel normalisation 4 pi/(2L+1)."""
    xs = np.array([32.0, 40.0, 60.0])
    f = phix_f(xs, L)
    with mp.workdps(30):
        for x, fv in zip(xs, f):
            s = phi_series(x, L)
            assert abs(fv / float(s) - 1) < 3e-14, (L, x, fv, s)


# ========================================================================================= (c) closed forms
def test_B1_closed_form():
    """From the tails (local rule, sympy): the leading power of S1 is eps^{-3/2} (the n = 3 exponent), its coefficient
    is A1 = 1/(12 pi), and eps^{3/2} S1 = A1 (1 + kappa1 eps) + O(eps^{3/2}) with kappa1 = 15/8.  Compared with the
    module constants.  Rejects A1 = 1/(3 pi^2) (P63L1's), 1/(6 pi), kappa1 = -15/8, and the split 1/2 + 1/2 P2."""
    P = derive_B1()
    lead = min(P)
    assert lead == -_q(EXP1) == -sp.Rational(3, 2), P
    A = P[lead]
    assert sp.simplify(A - 1 / (12 * sp.pi)) == 0, A
    assert abs(float(A) / A1 - 1) < 1e-15, (A, A1)
    kap = sp.simplify(P[lead + 1] / A)
    assert kap == _q(KAPPA1) == sp.Rational(15, 8), kap


def test_B2_closed_form():
    """From the tails: the eps^{-3/2} coefficient of S2 is exactly 0 (not merely small), the leading power is
    eps^{-1/2} with A2 = 3/(8 pi), and the eps^{1/2} coefficient is A2 kappa2, kappa2 = -3/8.  Rejects a
    non-cancelling root (weight != 0: eps^{-3/2} survives), the split 1/2 + 1/2 P2, dropping I1/y from the root
    kernel (A2 -> 1/(4 pi)), and an exponent -1."""
    P = derive_B2()
    lead = min(P)
    assert -sp.Rational(3, 2) not in P, P
    assert lead == -_q(EXP2) == -sp.Rational(1, 2), P
    A = P[lead]
    assert sp.simplify(A - 3 / (8 * sp.pi)) == 0, A
    assert abs(float(A) / A2 - 1) < 1e-15, (A, A2)
    kap = sp.simplify(P[lead + 1] / A)
    assert kap == _q(KAPPA2) == -sp.Rational(3, 8), kap


def test_B2_potentials_agree_through_one_over_x():
    """Phi0 and Phi2 of the same source agree through relative order 1/x (L enters at r^-2), so the B2 survivor comes
    only from the root kernels.  The two must differ at order x^-2 (else B2 would be the kernels alone at every
    order).  Rejects L entering the recurrence at order one."""
    src = f_p_ser() * E_full_ser()
    p0, p2 = local_potential(src, 0), local_potential(src, 2)
    assert p0.c[:2] == p2.c[:2] and p0.c[2] != p2.c[2]


def test_leading_amplitudes_by_hand():
    """Independent leading-order route, no Ser machinery: f_p ~ -(1/(2pi^2)) sqrt(pi/2) x^-1/2 e^-x, E ~ sqrt(2/pi)
    x^-1/2 e^-x, so f_p E ~ s0 x^-1 e^-2x with s0 = -1/(2 pi^2); Phi ~ -4 pi s0/k^2 x^-1 e^-2x (k = 2) = x^-1/(2 pi);
    H = x^2 Phi.  I_nu(y) ~ e^y (2 pi y)^-1/2 (1 - (4 nu^2 - 1)/(8y)).  Then
      A1 = LEG2 (1/(2pi)) Gamma(3/2)/sqrt(4 pi) = 1/(12 pi);
      B2: weights 1 - 1 cancel the lead; next, (I0 + I1/y) - I2 ~ [1/8 + 1 + 15/8]/y = 3/y, y = 2x at eps = 0:
      A2 = (1/(2pi)) (3/2) Gamma(1/2)/sqrt(4 pi) = 3/(8 pi).
    Rejects the same amplitude errors as the series derivation, by a different route."""
    s0 = -1 / (2 * sp.pi ** 2) * sp.sqrt(sp.pi / 2) * sp.sqrt(2 / sp.pi)
    phi0 = -4 * sp.pi * s0 / 4
    assert sp.simplify(phi0 - 1 / (2 * sp.pi)) == 0
    hank = lambda nu: -sp.Rational(4 * nu ** 2 - 1, 8)  # noqa: E731  (y^-1 coefficient of I_nu e^-y sqrt(2 pi y))
    A1h = _q(LEG2) * phi0 * sp.gamma(sp.Rational(3, 2)) / sp.sqrt(2 * sp.pi * 2)
    assert sp.simplify(A1h - 1 / (12 * sp.pi)) == 0, A1h
    w0, w2 = channel_weights()
    assert w0 - w2 == 0
    d1 = w0 * (hank(0) + 1) - w2 * hank(2)        # y^-1 coefficient of the kernel difference
    assert d1 == 3
    A2h = phi0 * d1 / 2 * sp.gamma(sp.Rational(1, 2)) / sp.sqrt(2 * sp.pi * 2)
    assert sp.simplify(A2h - 3 / (8 * sp.pi)) == 0, A2h
    assert abs(float(A1h) / A1 - 1) < 1e-15 and abs(float(A2h) / A2 - 1) < 1e-15


# ========================================================================================= (d) fast numeric
EPS_FAST = np.logspace(-2, -3, 7)


@pytest.fixture(scope="module")
def fast_scan():
    fr = FloatRoute(EPS_FAST[-1], 8001)
    return np.array([fr.S1(e) for e in EPS_FAST]), np.array([fr.S2(e) for e in EPS_FAST])


def test_fast_anchor_values(fast_scan):
    """float64 S1, S2 at eps = 1e-2 (the direct integrand from x = 1, no series) equal the record's referee values
    to 1e-11 (record FB0's bar is 1e-9; measured here 6.9e-15 and 3.8e-15).  The only fast test sensitive to the
    lower limit x = 1 (it moves only analytic terms of S) and to the overall normalisation of the integrand.
    Rejects a changed lower limit, a dropped e^{-eps x}, B2 channel weights that do not cancel, and I1/y dropped."""
    S1, S2 = fast_scan
    assert abs(S1[0] / float(S1_REF[0]) - 1) < 1e-11, (S1[0], S1_REF[0])
    assert abs(S2[0] / float(S2_REF[0]) - 1) < 1e-11, (S2[0], S2_REF[0])


def test_fast_B1_rate_and_amplitude(fast_scan):
    """eps in [1e-3, 1e-2] (7 points): r = eps^{3/2} S1/A1 -> 1 with r - 1 - kappa1 eps = O(eps^{3/2}) (|.| <= 6
    eps^{3/2} at every point; record: about -3.9 eps^{3/2}); log-slope within 0.02 of -3/2; the {1, u^2, u^3, u^4} fit
    gives A within 1e-6 of A1 (measured 1.7e-8) and u^2 coefficient / (A1 kappa1) in [0.95, 1.05].  Rejects A1 off by any factor
    (1/(3 pi^2), 1/(6 pi)), kappa1 sign flipped, the split 1/2 + 1/2 P2, and a fit exponent other than 3/2."""
    S1, _ = fast_scan
    r = EPS_FAST ** float(EXP1) * S1 / A1
    res = r - 1 - float(KAPPA1) * EPS_FAST
    assert np.all(np.abs(res) <= 6 * EPS_FAST ** 1.5), res / EPS_FAST ** 1.5
    assert abs(slope(EPS_FAST, S1) + 1.5) < 0.02, slope(EPS_FAST, S1)
    c = fit(EPS_FAST, S1, float(EXP1), [0, 2, 3, 4])
    assert abs(c[0] / A1 - 1) < 1e-6, c[0] / A1
    assert 0.95 < c[2] / (A1 * float(KAPPA1)) < 1.05, c[2] / (A1 * float(KAPPA1))


def test_fast_B2_rate_and_amplitude(fast_scan):
    """r = eps^{1/2} S2/A2 -> 1 with |r - 1| <= 1.0 eps^{1/2} (the u term, record about 0.48 u) and |r - 1| decreasing;
    log-slope of S2 in [-0.55, -0.45] (record local slopes -0.524 .. -0.509 on this range), in particular above -1
    (a non-cancelling root gives -3/2); the {1, u, u^2, u^3, u^4} fit gives A within 1e-5 of A2 (measured 4.7e-7) and
    u^2 / (A2 kappa2) in [0.95, 1.05] (measured 0.996).  Rejects a non-cancelling root, A2 off by a factor, I1/y dropped (A -> 1/(4pi)), and a fit
    exponent 1."""
    _, S2 = fast_scan
    r = EPS_FAST ** float(EXP2) * S2 / A2
    assert np.all(np.abs(r - 1) <= 1.0 * np.sqrt(EPS_FAST)), (r - 1) / np.sqrt(EPS_FAST)
    assert np.all(np.diff(np.abs(r - 1)) < 0), r
    s = slope(EPS_FAST, S2)
    assert -0.55 <= s <= -0.45, s
    c = fit(EPS_FAST, S2, float(EXP2), [0, 1, 2, 3, 4])
    assert abs(c[0] / A2 - 1) < 1e-5, c[0] / A2
    assert 0.95 < c[2] / (A2 * float(KAPPA2)) < 1.05, c[2] / (A2 * float(KAPPA2))


# ========================================================================================= (e) slow: referee + FB0-6
PANELS = [(1, 2), (2, 4), (4, 8), (8, 16), (16, 32)]
NGL = 40
XCUT = 32


def gl_nodes():
    t, w = np.polynomial.legendre.leggauss(NGL)
    xs, ws = [], []
    for a, b in PANELS:
        xs += list((b - a) / 2 * t + (a + b) / 2)
        ws += list((b - a) / 2 * w)
    return xs, ws


def kern_mp(y, which: int) -> mp.mpf:
    if which == 2:
        return mp.besseli(2, y) * mp.exp(-y)
    return (mp.besseli(0, y) + mp.besseli(1, y) / y) * mp.exp(-y)


def integrand_mp(shape: str, x, eps, p0, p2) -> mp.mpf:
    y = (2 - eps) * x
    damp = mp.exp(-eps * x) * x ** 2
    if shape == "B1":
        return kern_mp(y, 2) * damp * mp.mpf(LEG2.numerator) / LEG2.denominator * p2
    return damp * (W0 * kern_mp(y, 0) * p0 - W2 * kern_mp(y, 2) * p2)


def referee(shape: str, eps, cache: Dict[str, str]) -> mp.mpf:
    """GL panels on [1, 32] with the frozen mp30 node values; this file's 30-term local-rule series beyond 32."""
    xs, ws = gl_nodes()
    tot = mp.mpf(0)
    for x, w in zip(xs, ws):
        tot += w * integrand_mp(shape, mp.mpf(x), eps, mp.mpf(cache[f"0:{x!r}"]), mp.mpf(cache[f"2:{x!r}"]))
    f = lambda x: integrand_mp(shape, x, eps, phi_series(x, 0), phi_series(x, 2))  # noqa: E731
    tot += mp.quad(f, [XCUT, 64, 128, 256, 512, 1024, 4096, 16384, 65536, mp.inf])
    return tot


@pytest.mark.slow
def test_full_scan_slow():
    """The record's 9-point scan, both shapes.  Referee recomputed from the frozen node cache (sha256 checked) and this
    file's own local-rule series: equals the record's referee S to 1e-18 relative at every eps.  Float route (direct
    integrand from x = 1, Simpson, 20001 log nodes to 1 + 8e5): FB0, agrees with the referee to 1e-9.  On the referee:
    FB1 last-5 slope in [-1.55, -1.45]; FB2 |A_fit/A1 - 1| <= 1e-5 (smallest 6, {1, u^2, u^3, u^4}); FB3 u^2 ratio in
    [0.95, 1.05]; FB4 last-5 slope in [-0.75, -0.25]; FB5 |A_fit/A2 - 1| <= 1e-4 (smallest 6, {1, u, .., u^4}); FB6
    u^2 ratio in [0.9, 1.1]; S > 0 everywhere.  The same fits on the float route meet FB2/FB5 too (the
    series-independent amplitude).  Rejects any error the fast tier rejects at the record's bars, plus a broken
    series beyond x = 32 or a changed cache."""
    raw = CACHE.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == CACHE_SHA
    cache = json.loads(raw)
    fr = FloatRoute(EPS_FULL[-1], 20001)
    for shape, ref, A, kap, pw, pows, fbs, fba, fbk in (
            ("B1", S1_REF, A1, KAPPA1, EXP1, [0, 2, 3, 4], (-1.55, -1.45), 1e-5, (0.95, 1.05)),
            ("B2", S2_REF, A2, KAPPA2, EXP2, [0, 1, 2, 3, 4], (-0.75, -0.25), 1e-4, (0.9, 1.1))):
        with mp.workdps(30):
            Sr = [referee(shape, mp.mpf(e), cache) for e in EPS_FULL]
            for s, rv in zip(Sr, ref):
                assert abs(s / mp.mpf(rv) - 1) < 1e-18, (shape, s, rv)
        Sr = np.array([float(s) for s in Sr])
        Sf = np.array([fr.S1(e) if shape == "B1" else fr.S2(e) for e in EPS_FULL])
        assert np.max(np.abs(Sf / Sr - 1)) <= 1e-9, (shape, Sf / Sr - 1)
        assert np.all(Sr > 0)
        s5 = slope(EPS_FULL[-5:], Sr[-5:])
        assert fbs[0] <= s5 <= fbs[1], (shape, s5)
        for S in (Sr, Sf):
            c = fit(EPS_FULL[3:], S[3:], float(pw), pows)
            assert abs(c[0] / A - 1) <= fba, (shape, c[0] / A)
            assert fbk[0] <= c[2] / (A * float(kap)) <= fbk[1], (shape, c[2] / (A * float(kap)))


# TIMINGS (2026-11-17, this machine): fast tests 15 s in all (14 tests; the 8001-node float scan ~11 s, shared by the
# fast tier through a module fixture; the sympy derivations < 0.5 s each).  Slow test 72 s (--slow; 18 referee
# evaluations at mp30 from the frozen cache, plus the 20001-node float route).  Measured: fast anchors 6.9e-15 (B1),
# 3.8e-15 (B2) vs the record's referee; fast B1 slope -1.4945, A/A1 - 1 = 1.7e-8, u^2/(A1 kappa1) 0.99997,
# (r - 1 - kappa1 eps)/eps^{3/2} from -3.83 to -3.94; fast B2 slope -0.5155, A/A2 - 1 = 4.7e-7, u^2/(A2 kappa2) 0.9959,
# (r - 1)/u from -0.503 to -0.484.  Slow: referee vs record S 4.4e-22 (B1), 1.9e-22 (B2); FB0 7.1e-15, 2.9e-13;
# slopes -1.49934, -0.50459; smallest-6 A ratios 1 + 2.3e-10, 1 + 3.6e-9 (referee and float route alike); u^2 ratios
# 0.999997, 0.99975 -- the record's numbers reproduced.
#
# FIRE-TEST PLANTS (python debug/qa/fire_test.py tests/test_paper63_attainment_l2.py --plant "OLD=>NEW" -k <sel>),
# 2026-11-17, written and fired as a separate step after the guard; 41 runs, 38 fired; the 3 that did not are below.
#     A1 = 1/(3 pi^2) (P63L1's)          "A1 = 1 / (12 * math.pi)" -> "1 / (3 * math.pi ** 2)"   -k B1_closed; fast_B1    FIRED, FIRED
#     A1 = 1/(6 pi)                      "A1 = 1 / (12 * math.pi)" -> "1 / (6 * math.pi)"        -k B1_closed; fast_B1; by_hand   FIRED x3
#     A2 = 1/(4 pi)                      "A2 = 3 / (8 * math.pi)" -> "2 / (8 * math.pi)"         -k B2_closed; fast_B2    FIRED, FIRED
#     kappa1 sign flipped                "KAPPA1 = Fraction(15, 8)" -> Fraction(-15, 8)          -k B1_closed; fast_B1    FIRED, FIRED
#     kappa2 sign flipped                "KAPPA2 = Fraction(-3, 8)" -> Fraction(3, 8)            -k B2_closed; fast_B2    FIRED, FIRED
#     non-cancelling root (5/2 -> 2)     "B2_COEF = Fraction(5, 2)" -> Fraction(2, 1)            -k legendre_split; B2_closed; by_hand; fast_B2; anchor   FIRED x5
#     split 1/2 + 1/2 P2                 "LEG0, LEG2 = Fraction(1, 3), Fraction(2, 3)" -> (1/2, 1/2)   -k legendre_split; B1_closed; B2_closed; fast_B1; fast_B2   FIRED x5
#     B2 exponent -1                     "EXP2 = Fraction(1, 2)" -> Fraction(1, 1)               -k B2_closed; fast_B2    FIRED, FIRED
#     drop I1/y, derivation              "kernel_I(0, True), H0" -> "kernel_I(0, False), H0"      -k B2_closed             FIRED
#     drop I1/y, float route             "(i0e(yy) + i1e(yy) / yy)" -> "(i0e(yy) + 0 * i1e(yy) / yy)"   -k anchor; fast_B2   FIRED, FIRED
#     drop I1/y, hand route              "d1 = w0 * (hank(0) + 1)" -> "(hank(0) + 0)"             -k by_hand               FIRED
#     drop I1/y, referee (unique anchor) "return (mp.besseli(0, y) + mp.besseli(1, y) / y) * mp.exp(-y)" -> "+ 0 * ..."   -k slow   FIRED
#     recurrence 2k(b+1) -> 2k(b+2)      in local_potential                                       -k series_vs; P63L1_tail   FIRED, FIRED
#     recurrence L(L+1) -> L(L+2)        in local_potential                                       -k series_vs; B2_closed    FIRED, FIRED
#     lower limit 1 -> 1.5 (float)       "self.y = np.linspace(0.0," -> "np.linspace(np.log(1.5),"   -k anchor; slow (FB0)   FIRED, FIRED
#     e^{-eps x} dropped (float B1)      "... self.p2 * np.exp(-eps * x)" -> "np.exp(-0 * eps * x)"   -k anchor            FIRED
#     series 30 -> 6 terms (referee)     "NSER = 30" -> "NSER = 6"                                -k slow                  FIRED
#     A1 off by 1e-4 (FB2 bar)           "1 / (12 * math.pi)" -> "1 / (12.0012 * math.pi)"        -k slow                  FIRED
#     A2 off by 1e-3 (FB5 bar)           "3 / (8 * math.pi)" -> "3 / (8.008 * math.pi)"           -k slow                  FIRED
#   Did not fire, and why (not a guard defect):
#     lower limit 1 -> 1.5 with -k fast_B1 or -k fast_B2 alone: a lower limit moves S by a function holomorphic at
#       eps = 0, i.e. only integer powers of eps in S (the u^3 term of eps^{3/2} S1, the u term of eps^{1/2} S2, both
#       fitted); amplitude, slope and kappa are blind to it by construction.  The same plant fires the anchor test
#       and, in the slow tier, FB0.
#     drop I1/y in the referee, first attempt with the bare anchor "(mp.besseli(0, y) + mp.besseli(1, y) / y)": that
#       anchor occurs twice and fire_test.py mutates the first occurrence, which is in test_root_kernel_jumps, not in
#       kern_mp, so the slow test ran unmutated (fire_test prints a note).  Re-fired with the unique anchor above: FIRED.
