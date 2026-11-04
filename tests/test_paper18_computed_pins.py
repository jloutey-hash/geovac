"""Paper 18 computed pins (group3 guard pass, 2026-11-03).

1. eq:dirac_modular_residual.  The previous guard asserted dict entries of
   geovac.thermal_tensor_triple.modular_residual_dirac_S3 (plant 2 m^2 -> 3 m^2
   did not fire).  Here the module's epsilon_term_m is summed over m and
   compared with epsilon(t) = K(t) - [(sqrt(pi)/2) t^{-3/2} - (sqrt(pi)/4) t^{-1/2}],
   K(t) = sum_n 2(n+1)(n+2) e^{-t (n+3/2)^2} summed directly; > 100 digits at
   t = 1, 2 and to 1e-50 across t in [1/2, 10].

2. eq:fock_weight_values / eq:kappa_B_identity.  <n,l|(p^2+p0^2)^{-2}|n,l> at
   p0 = 1 derived from the hydrogenic momentum wavefunctions
   F_nl(p) ~ p^l (p^2+1)^{-(l+2)} C^{l+1}_{n-l-1}((p^2-1)/(p^2+1)), by exact
   Beta integrals (Fractions), cross-checked by mpmath quadrature; then the
   Hopf-weighted sum 63/4 = 6 B |kappa| with B = 42, kappa = -1/16.

3. alpha^2-weighted Ihara zeta, Rule A, n_max = 3 (Paper 18 "alpha^2-weighted
   Ihara zeta" paragraph): zeta^{-1} has alpha-degree 44 and its alpha^2
   coefficient is -35/648 x a primitive integer polynomial.  Recomputed
   independently of geovac.ihara_zeta_weighted: own f_SO weights from the
   Dirac labels, own Mizuno-Sato determinant det(I - s A_w + s^2 Q_w) by
   modular Gaussian elimination at two large primes, Lagrange interpolation
   in u = alpha^2 (degree) and in s (Jacobi d/du det = det M0 tr(M0^{-1} M1)).
   Only the 0/1 graph topology is taken from production
   (geovac.ihara_zeta_dirac.build_dirac_s3_graph).
"""
from __future__ import annotations

from fractions import Fraction
from math import factorial

import mpmath as mp
import pytest
import sympy as sp

from geovac.thermal_tensor_triple import modular_residual_dirac_S3

# ---------------------------------------------------------------------------
# 1. Dirac modular residual
# ---------------------------------------------------------------------------


def _K_direct(t: mp.mpf) -> mp.mpf:
    tot = mp.mpf(0)
    n = 0
    while True:
        term = 2 * (n + 1) * (n + 2) * mp.exp(-t * (n + mp.mpf(3) / 2) ** 2)
        tot += term
        if n > 5 and term < mp.mpf(10) ** (-mp.mp.dps - 10):
            return tot
        n += 1


def _eps_module(t: mp.mpf) -> mp.mpf:
    expr = modular_residual_dirac_S3()["epsilon_term_m"]
    ts = sp.Symbol("t", positive=True)
    ms = sp.Symbol("m", integer=True, positive=True)
    f = sp.lambdify((ms, ts), expr, modules="mpmath")
    tot = mp.mpf(0)
    m = 1
    while True:
        term = f(m, t)
        tot += term
        if abs(term) < mp.mpf(10) ** (-mp.mp.dps - 10):
            return tot
        m += 1


@pytest.mark.parametrize("t,digits", [("0.5", 50), ("1", 100), ("2", 100), ("5", 50), ("10", 50)])
def test_dirac_modular_residual_computed(t, digits):
    mp.mp.dps = 130
    t = mp.mpf(t)
    two_term = mp.sqrt(mp.pi) / 2 * t ** mp.mpf(-1.5) - mp.sqrt(mp.pi) / 4 * t ** mp.mpf(-0.5)
    eps_direct = _K_direct(t) - two_term
    assert abs(eps_direct) > mp.mpf(10) ** -8          # non-trivial residual
    assert abs(eps_direct - _eps_module(t)) < mp.mpf(10) ** -digits


# ---------------------------------------------------------------------------
# 2. Fock-weight values and the kappa-B identity
# ---------------------------------------------------------------------------

_q = sp.Symbol("q")
_x = sp.Symbol("x")
CELLS = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)]
PAPER_WEIGHTS = [Fraction(7, 16), Fraction(5, 8), Fraction(3, 8),
                 Fraction(5, 8), Fraction(17, 32), Fraction(11, 32)]


def _r(k: int) -> Fraction:
    # Gamma(k + 1/2) / sqrt(pi)
    return Fraction(factorial(2 * k), 4 ** k * factorial(k))


def _beta_int_over_pi(a: int, b: int) -> Fraction:
    # int_0^inf p^{2a} (1 + p^2)^{-b} dp / pi
    return _r(a) * _r(b - a - 1) / (2 * factorial(b - 1))


def _num_poly(n: int, l: int) -> sp.Expr:
    d = n - l - 1
    C = sp.gegenbauer(d, l + 1, _x)
    Ct = sp.expand(sp.cancel((_q + 1) ** d * C.subs(_x, (_q - 1) / (_q + 1))))
    return sp.expand(_q ** (l + 1) * Ct ** 2)   # p^2 F^2 (1+p^2)^{2n+2}, q = p^2


def _int(poly: sp.Expr, b: int) -> Fraction:
    tot = Fraction(0)
    for (a,), c in sp.Poly(poly, _q).terms():
        tot += Fraction(int(sp.numer(c)), int(sp.denom(c))) * _beta_int_over_pi(a, b)
    return tot


def _fock_weight(n: int, l: int) -> Fraction:
    P = _num_poly(n, l)
    return _int(P, 2 * n + 4) / _int(P, 2 * n + 2)


def test_fock_weight_values_derived():
    got = [_fock_weight(n, l) for (n, l) in CELLS]
    assert got == PAPER_WEIGHTS, got


def test_fock_weight_quadrature_cross_check():
    mp.mp.dps = 30
    for (n, l), w in zip(CELLS, PAPER_WEIGHTS):
        F = lambda p: p ** l / (p ** 2 + 1) ** (l + 2) * mp.gegenbauer(n - l - 1, l + 1, (p ** 2 - 1) / (p ** 2 + 1))
        N = mp.quad(lambda p: F(p) ** 2 * p ** 2, [0, 1, mp.inf])
        W = mp.quad(lambda p: F(p) ** 2 * p ** 2 / (p ** 2 + 1) ** 2, [0, 1, mp.inf])
        assert abs(W / N - mp.mpf(w.numerator) / w.denominator) < mp.mpf(10) ** -20


def test_kappa_B_identity():
    total = sum((2 * l + 1) * l * (l + 1) * _fock_weight(n, l) for (n, l) in CELLS)
    B, kappa = 42, Fraction(-1, 16)
    assert total == Fraction(63, 4) == 6 * B * abs(kappa)


# ---------------------------------------------------------------------------
# 3. alpha^2-weighted Ihara zeta, Rule A, n_max = 3
# ---------------------------------------------------------------------------

PRIMES = (2 ** 61 - 1, 2 ** 89 - 1)


@pytest.fixture(scope="module")
def graph():
    from geovac.ihara_zeta_dirac import build_dirac_s3_graph
    A01, labels, _, _ = build_dirac_s3_graph(3, "A")
    V = A01.shape[0]

    def f_so(lab):
        k = lab.kappa
        l = k if k > 0 else -k - 1
        if l == 0:
            return Fraction(0)
        return Fraction(-(k + 1)) / (4 * lab.n_fock ** 3 * l * Fraction(2 * l + 1, 2) * (l + 1))

    f = [f_so(x) for x in labels]
    nbr = [[j for j in range(V) if A01[i, j]] for i in range(V)]
    om = {(i, j): abs(f[i] + f[j]) / 2 for i in range(V) for j in nbr[i]}
    E = sum(len(x) for x in nbr) // 2
    # components via union-find
    par = list(range(V))

    def find(a):
        while par[a] != a:
            par[a] = par[par[a]]
            a = par[a]
        return a

    for i in range(V):
        for j in nbr[i]:
            par[find(i)] = find(j)
    c = len({find(i) for i in range(V)})
    return V, nbr, om, E, c


def _red(x: Fraction, p: int) -> int:
    x = Fraction(x)
    return x.numerator % p * pow(x.denominator, -1, p) % p


def _mats(graph, s: int, p: int):
    V, nbr, om, _, _ = graph
    M0 = [[0] * V for _ in range(V)]
    M1 = [[0] * V for _ in range(V)]
    for i in range(V):
        M0[i][i] = (1 + s * s * (len(nbr[i]) - 1)) % p
        M1[i][i] = s * s * _red(sum(om[(i, j)] for j in nbr[i]), p) % p
        for j in nbr[i]:
            M0[i][j] = (M0[i][j] - s) % p
            M1[i][j] = (M1[i][j] - s * _red(om[(i, j)], p)) % p
    return M0, M1


def _det_solve(M, B, p):
    n = len(M)
    A = [M[i][:] + B[i][:] for i in range(n)]
    det = 1
    for col in range(n):
        piv = next((r for r in range(col, n) if A[r][col]), None)
        if piv is None:
            return 0, None
        if piv != col:
            A[col], A[piv] = A[piv], A[col]
            det = -det
        det = det * A[col][col] % p
        inv = pow(A[col][col], -1, p)
        A[col] = [v * inv % p for v in A[col]]
        for r in range(n):
            if r != col and A[r][col]:
                fct = A[r][col]
                A[r] = [(v - fct * w) % p for v, w in zip(A[r], A[col])]
    return det % p, [row[n:] for row in A]


def _interp(pts, p):
    n = len(pts)
    coeffs = [0] * n
    for i, (xi, yi) in enumerate(pts):
        num, den = [1], 1
        for j, (xj, _) in enumerate(pts):
            if j != i:
                num = [(a - xj * b) % p for a, b in zip([0] + num, num + [0])]
                den = den * (xi - xj) % p
        fac = yi * pow(den, -1, p) % p
        coeffs = [(cc + fac * a) % p for cc, a in zip(coeffs, num)]
    return coeffs


@pytest.mark.parametrize("p", PRIMES)
def test_weighted_ihara_alpha_degree_44(graph, p):
    V = graph[0]
    for s in (2, 3):
        M0, M1 = _mats(graph, s, p)
        pts = []
        for u in range(V + 2):
            M = [[(a + u * b) % p for a, b in zip(r0, r1)] for r0, r1 in zip(M0, M1)]
            pts.append((u, _det_solve(M, [[0]] * V, p)[0]))
        cu = _interp(pts, p)
        u_degree = max(k for k, cc in enumerate(cu) if cc)
        assert 2 * u_degree == 44, (s, u_degree)


@pytest.mark.parametrize("p", PRIMES)
def test_weighted_ihara_alpha2_coefficient(graph, p):
    V, _, _, E, c = graph
    assert (V, E, c) == (28, 29, 5)
    r = E - V + c
    pts = []
    for s in range(2, 50):
        M0, M1 = _mats(graph, s, p)
        d, X = _det_solve(M0, M1, p)
        pts.append((s, d * sum(X[i][i] for i in range(V)) % p))
    cs = _interp(pts, p)
    # multiply by (1 - s^2)^(r - c)
    for _ in range(r - c):
        nxt = [0] * (len(cs) + 2)
        for k, a in enumerate(cs):
            nxt[k] = (nxt[k] + a) % p
            nxt[k + 2] = (nxt[k + 2] - a) % p
        cs = nxt
    s_ = sp.Symbol("s")
    ref = sp.Poly(-35 * s_**4 * (s_ - 1)**6 * (s_ + 1)**6 * (s_**2 + 1) * (s_**2 - s_ + 1)
                  * (s_**2 + s_ + 1) * (2*s_**3 - s_**2 + s_ - 1) * (2*s_**3 + s_**2 + s_ + 1)
                  * (60*s_**14 + 114*s_**12 + 80*s_**10 + 27*s_**8 - 6*s_**6 - 26*s_**4
                     - 20*s_**2 - 5) / 648, s_)
    # the reference (RH-K memo factorisation) has content 35/648, negative leading sign
    content, prim = ref.primitive()
    assert content == sp.Rational(35, 648) and prim.LC() < 0
    assert all(cs[k] == _red(Fraction(str(ref.coeff_monomial(s_ ** k))), p)
               for k in range(len(cs))), "alpha^2 coefficient differs from -35/648 x primitive"
