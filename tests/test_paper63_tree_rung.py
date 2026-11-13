"""Paper 63, Proposition prop:tree_rung: the three-dimensional Coulomb tree is bounded by the n-Bessel moments.

A = C[z]<theta>, z = r^2, theta z = z (theta + 2).  L_m = annihilator of K0^m (Sym^m of theta^2 - z, cyclic vector
K0^m).  phi: z -> u theta^2 (u = t^-2, theta u = u (theta - 2)), theta -> -1 - theta: the integration-by-parts
homomorphism of T(g) = int_0^oo K0(t r) g(r) dr.

Checked (each names the wrong answer it rejects):
- phi(L_2), phi(L_3) are F_3's and F_4's operators exactly (rejects phi(theta) = 1 + theta or theta - 1, and a
  wrong symmetric power);
- weight bound: every term z^k theta^j of L_m(theta - lam) has 2k + j <= m + 1, so phi gives order <= m + 1
  (rejects a piece operator of order above m + 1; the filtration step itself is backed by the proof and by the
  n = 3, 4 factor sets, not by this file);
- phi(L_m) annihilates F_{m+1}; phi(L_m(theta - 1)) annihilates int x K0(t x) K0(x)^m dx (numerical, 25 digits);
- the even-twist intertwiner phi(L(theta - lam - 2)) u theta^2 = u theta^2 phi(L(theta - lam)) (exact);
- the transfer T(theta g) = -(1 + theta_t) T(g), T(r^2 g) = t^-2 theta_t^2 T(g) on a tree function (rejects
  T(theta g) = +(1 + theta_t) T(g) and a missing boundary term);
- the radial Poisson equation theta (theta + 1) Phi W = -4 pi r^2 W of a Coulomb potential (rejects theta(theta-1));
- Clebsch-Gordan: Sym^a (x) Sym^b -> Sym^(a+b) has kernel of rank a b, hit injectively from Sym^(a-1) (x) Sym^(b-1).
"""
import itertools

import mpmath as mp
import numpy as np
import pytest
import sympy as sp
from scipy.integrate import cumulative_trapezoid
from scipy.special import k0, k1

z, th, y = sp.symbols("z theta y")


def _a_mul(P, Q, sgn):
    """normal form x^k theta^j with theta x = x (theta + 2 sgn)."""
    out = {}
    for (a, b), c1 in P.items():
        for (c, d), c2 in Q.items():
            for (e,), cf in sp.Poly(sp.expand((th + 2 * sgn * c) ** b * th**d), th).terms():
                out[(a + c, e)] = out.get((a + c, e), 0) + c1 * c2 * cf
    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}


def _sym_op(m):
    def act(vec):
        out = [sp.Integer(0)] * (m + 1)
        for a, p in enumerate(vec):
            if p == 0:
                continue
            out[a] += 2 * z * sp.diff(p, z)
            if a < m:
                out[a + 1] += (m - a) * p
            if a > 0:
                out[a - 1] += a * z * p
        return [sp.expand(q) for q in out]
    vs = [[sp.Integer(1)] + [sp.Integer(0)] * m]
    for _ in range(m + 1):
        vs.append(act(vs[-1]))
    M = sp.Matrix([[vs[j][a] for j in range(m + 1)] for a in range(m + 1)])
    q = M.LUsolve(sp.Matrix(vs[m + 1]))
    L = {(0, m + 1): sp.Integer(1)}
    for j, e in enumerate(q):
        for (k,), cf in sp.Poly(sp.expand(sp.simplify(e)), z).terms():
            L[(k, j)] = L.get((k, j), 0) - cf
    return {k: v for k, v in L.items() if v != 0}


def _shift(L, lam):
    out = {}
    for (k, j), c in L.items():
        for (e,), cf in sp.Poly(sp.expand((th - lam) ** j), th).terms():
            out[(k, e)] = out.get((k, e), 0) + c * cf
    return {k: v for k, v in out.items() if v != 0}


PHI_THETA = {(0, 0): sp.Integer(-1), (0, 1): sp.Integer(-1)}


def _phi(L, phi_theta=PHI_THETA):
    res = {}
    for (k, j), c in L.items():
        term = {(0, 0): sp.Integer(c)}
        for _ in range(k):
            term = _a_mul(term, {(1, 2): sp.Integer(1)}, -1)
        for _ in range(j):
            term = _a_mul(term, phi_theta, -1)
        for kk, v in term.items():
            res[kk] = res.get(kk, 0) + v
    return {k: sp.expand(v) for k, v in res.items() if sp.expand(v) != 0}


def _t_lists(P):
    K = max(k for k, _ in P)
    Q = [sp.Integer(0)] * (K + 1)
    for (k, j), c in P.items():
        Q[K - k] += c * th**j
    return [sp.expand(q) for q in Q]


def _proportional(A, B):
    if len(A) != len(B):
        return False
    r = sp.cancel(A[-1] / B[-1])
    return r.is_number and all(sp.expand(a - r * b) == 0 for a, b in zip(A, B))


F3 = [4 * th**3, -(th + 1) ** 3]                                                  # sec:density_ring
F4 = [9 * th**2 * (th - 2) ** 2, -2 * th**2 * (5 * th**2 + 1), (th + 1) ** 4]     # eq:four_density_op, in t


def test_phi_of_sym_powers_is_the_ladder():
    """phi(L_2) = F_3's operator and phi(L_3) = F_4's, exactly. Rejects phi(theta) = 1 + theta."""
    assert _proportional(_t_lists(_phi(_sym_op(2))), F3)
    assert _proportional(_t_lists(_phi(_sym_op(3))), F4)
    wrong = {(0, 0): sp.Integer(1), (0, 1): sp.Integer(1)}
    assert not _proportional(_t_lists(_phi(_sym_op(3), wrong)), F4)


@pytest.mark.parametrize("m", range(1, 7))
def test_weight_bound_gives_order_at_most_m_plus_1(m):
    """Every term of L_m(theta - lam) has 2k + j <= m + 1; phi(L_m(theta - lam)) has order <= m + 1, = m + 1 at
    lam = 0. Rejects a piece operator of order above m + 1."""
    L = _sym_op(m)
    assert max(j for (_, j) in L) == m + 1 and L[(0, m + 1)] == 1
    for lam in range(-3, 4):
        Ls = _shift(L, lam)
        assert max(2 * k + j for (k, j) in Ls) <= m + 1
        assert max(j for (_, j) in _phi(Ls)) <= m + 1
    assert max(j for (_, j) in _phi(L)) == m + 1


@pytest.mark.parametrize("m,lam", [(2, 0), (2, 1), (3, 0), (3, 2)])
def test_even_shift_intertwiner(m, lam):
    """phi(L(theta - lam - 2)) u theta^2 = u theta^2 phi(L(theta - lam)), exactly."""
    L = _sym_op(m)
    U = {(1, 2): sp.Integer(1)}
    lhs = _a_mul(_phi(_shift(L, lam + 2)), U, -1)
    rhs = _a_mul(U, _phi(_shift(L, lam)), -1)
    assert set(lhs) == set(rhs) and all(sp.expand(lhs[k] - rhs[k]) == 0 for k in lhs)


def _theta_derivs(J):
    a0, a1 = sp.Integer(1), sp.Integer(0)
    out = []
    for _ in range(J + 1):
        out.append((sp.lambdify(y, a0, "mpmath"), sp.lambdify(y, a1, "mpmath")))
        a0, a1 = sp.expand(y * sp.diff(a0, y) - a1 * y), sp.expand(-a0 * y + y * sp.diff(a1, y) - a1)
    return out


def _residual(P, m, w, tv):
    J = max(j for (_, j) in P)
    vals = []
    for g0, g1 in _theta_derivs(J):
        vals.append(mp.quad(lambda x: (g0(tv * x) * mp.besselk(0, tv * x) + g1(tv * x) * mp.besselk(1, tv * x))
                            * x**w * mp.besselk(0, x) ** m, [0, 1, 4, 12, mp.inf]))
    return sum(c * tv ** (-2 * k) * vals[j] for (k, j), c in P.items()) / max(abs(v) for v in vals)


@pytest.mark.parametrize("m", [1, 2, 3])
def test_phi_annihilates_even_and_odd_moments(m):
    """phi(L_m) kills F_{m+1} = int K0(tx) K0^m; phi(L_m(theta - 1)) kills int x K0(tx) K0^m (25 digits)."""
    mp.mp.dps = 25
    L = _sym_op(m)
    for tv in (mp.mpf("0.7"), mp.mpf("1.9")):
        assert abs(_residual(_phi(L), m, 0, tv)) < mp.mpf(10) ** -18
        assert abs(_residual(_phi(_shift(L, 1)), m, 1, tv)) < mp.mpf(10) ** -18
    # the odd moment is not killed by the even operator
    assert abs(_residual(_phi(L), m, 1, mp.mpf("1.9"))) > 1e-6


# ------------------------------------------------------------------ grid: Poisson and the transfer
U_ = np.linspace(np.log(1e-9), np.log(80.0), 400001)
R = np.exp(U_)


def _cum0(f):
    return cumulative_trapezoid(f * R, U_, initial=0.0)


def _coulomb(W):
    return 4 * np.pi * (_cum0(W * R**2) / R + (_cum0(W * R)[-1] - _cum0(W * R)))


def _theta(f):
    return np.gradient(f, U_)


@pytest.mark.parametrize("which", ["leaf", "inner"])
def test_radial_poisson(which):
    """theta (theta + 1) Phi W = -4 pi r^2 W. Rejects theta (theta - 1)."""
    W = k0(R) if which == "leaf" else k0(R) * _coulomb(k0(R))
    P = _coulomb(W)
    lhs = _theta(_theta(P) + P)
    sel = (R > 0.05) & (R < 15)
    rhs = -4 * np.pi * R[sel] ** 2 * W[sel]
    assert np.allclose(lhs[sel], rhs, rtol=1e-4, atol=1e-6 * np.abs(rhs).max())
    wrong = _theta(_theta(P) - P)
    assert not np.allclose(wrong[sel], -4 * np.pi * R[sel] ** 2 * W[sel], rtol=1e-2)


def test_transfer_homomorphism_on_a_tree_function():
    """T(theta g) = -(1 + theta_t) T(g) and T(r^2 g) = t^-2 theta_t^2 T(g) for g = r^2 Phi^3 (the star of four
    rooted at the centre). theta_t is taken under the integral, theta K0 = -x K1, theta^2 K0 = x^2 K0."""
    g = R**2 * _coulomb(k0(R)) ** 3
    for tv in (0.8, 1.6):
        def T(f, kern=lambda x: k0(x)):
            return np.trapezoid(kern(tv * R) * f * R, U_)
        Tg = T(g)
        th_Tg = T(g, lambda x: -x * k1(x))
        th2_Tg = T(g, lambda x: x**2 * k0(x))
        assert T(_theta(g)) == pytest.approx(-(Tg + th_Tg), rel=1e-6)
        assert T(R**2 * g) == pytest.approx(th2_Tg / tv**2, rel=1e-6)
        assert abs(T(_theta(g)) - (Tg + th_Tg)) > 1e-2 * abs(Tg)


@pytest.mark.parametrize("a,b", [(1, 1), (2, 1), (2, 2), (3, 2), (4, 3)])
def test_clebsch_gordan_kernel(a, b):
    """Multiplication Sym^a (x) Sym^b -> Sym^(a+b) (V = <x0, x1>) is onto with kernel of rank a b, and
    f (x) g (x) (x0^x1) -> x0 f (x) x1 g - x1 f (x) x0 g is injective into it."""
    x0, x1 = sp.symbols("x0 x1")
    ma = [x0 ** (a - i) * x1**i for i in range(a + 1)]
    mb = [x0 ** (b - i) * x1**i for i in range(b + 1)]
    mab = [x0 ** (a + b - i) * x1**i for i in range(a + b + 1)]
    pairs = list(itertools.product(range(a + 1), range(b + 1)))
    mult = sp.Matrix([[sp.Poly(ma[i] * mb[j], x0, x1).coeff_monomial(mm) for (i, j) in pairs] for mm in mab])
    assert mult.rank() == a + b + 1
    la = [x0 ** (a - 1 - i) * x1**i for i in range(a)]
    lb = [x0 ** (b - 1 - i) * x1**i for i in range(b)]
    cols = []
    for f, gg in itertools.product(la, lb):
        v = {}
        for (p, q, sgn) in ((x0 * f, x1 * gg, 1), (x1 * f, x0 * gg, -1)):
            v[(ma.index(sp.expand(p)), mb.index(sp.expand(q)))] = v.get((ma.index(sp.expand(p)), mb.index(sp.expand(q))), 0) + sgn
        cols.append([v.get(pq, 0) for pq in pairs])
    inc = sp.Matrix(cols).T
    assert inc.rank() == a * b
    assert (mult * inc).is_zero_matrix
