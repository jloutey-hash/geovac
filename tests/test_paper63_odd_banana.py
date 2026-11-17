"""Paper 63, Proposition prop:odd_banana: odd-lambda pieces carry two-dimensional banana operators.

Written as a separate guard pass after the proof (CLAUDE.md sec. 9); no debug/ import.  The operator algebra
(_a_mul, _sym_op, _shift, _phi) is copied from tests/test_paper63_tree_rung.py so this file stands alone.

Conventions.  An operator is a dict {(k, j): c} meaning sum c u^k theta^j, u = t^-2 on the LEFT, theta = t d/dt,
theta u = u (theta - 2); k may be negative (t^{2a} = u^{-a}).  L_m = annihilator of K0(r)^m in C[r^2]<theta_r>;
phi(r^2) = u theta^2, phi(theta_r) = -1 - theta.  Vanhove's operator (arXiv:1401.6438, Table 1, equal masses,
s = t^2) is V_m = sum_k q_k(s) d^k/ds^k with s^k d^k/ds^k = prod_{j<k}(theta_s - j), theta_s = theta / 2.

Each test names the wrong answer it rejects:
 1. test_odd_piece_is_theta2_times_vanhove: phi(L_m(theta - 1)) = (-2)^{m-1} t^{2 a_m} theta^2 V_m EXACTLY,
    a_m = 0, 1, 1, 2 (m = 2..5), with phi(L_m(theta - 1)) in the cleared t-list form sum_i t^{2i} Q_i(theta),
    i.e. t^{2 K_m} times the homomorphism image, K_m = floor((m+1)/2).  For the homomorphism itself the same
    identity reads phi(L_m(theta - 1)) = (-2)^{m-1} U V_m with U = t^-2 theta^2 for every m; both are checked.  Rejects a wrong sign/power of the constant, a wrong a_m, a rescaled or
    perturbed V_m, a shifted left factor (theta - r)^2, and phi(L_m(theta - 3)) in place of theta - 1.
 2. test_vanhove_does_not_divide_even_piece: V_m does not right-divide phi(L_m) (remainder != 0) in
    Q(t)<theta>, while it does right-divide phi(L_m(theta - 1)) (positive control through the SAME division
    routine, which is independent of test 1's normal-form comparison).  Rejects "the even piece also carries V_m".
 3. test_parity_intertwiner: phi(L_m(theta - lam - 2)) U = U phi(L_m(theta - lam)), U = t^-2 theta^2, exact, for
    m = 2..5 and odd and even lam (rejects U = t^-2 (theta - 1)^2, a shift by 1 instead of 2, and phi(theta) =
    +1 + theta); plus a non-commutation control so the identity is not U commuting with everything.
 4. test_m1_odd_moment_is_elementary: int_0^oo x I0(t x) K0(x) dx = 1/(1 - t^2) at two t (mpmath), and
    phi(L_1(theta - 1)) annihilates 1/(1 - t^2) exactly (sympy).  Rejects 1/(1 + t^2), 1/(1 - t^2)^2, and an
    operator that does not see the closed form.
 5. test_vanhove_constant / test_odd_operator_annihilates_banana: V_m applied to b_m(t) = int x I0(t x) K0(x)^m dx
    is the constant -m!/2^{m-1} at two t (rejects -m!/2^m, the opposite sign, and a rescaled V_m);
    phi(L_m(theta - 1)) annihilates b_m (rejects a wrong a_m-free operator, with an even-operator control that
    must NOT annihilate).  m = 2, 3 fast; m = 4 @slow.
 6. test_parity_bookkeeping_*: part (i).  See the docstrings; the rule is checked on the s-chain of three at the
    function level (the leaf potential = exterior tail + even Bessel part, from the shell-theorem integral), and
    its prediction for the four-density star is checked against the measured factor set of the order-8
    annihilator W8 (record docs/predictions/2026-10-21_four_density_star.md, guarded by
    tests/test_paper63_four_density_star.py): F_4 and the Legendre operator present, no order-3 (K3) piece.
"""
from __future__ import annotations

import math

import mpmath as mp
import pytest
import sympy as sp
from scipy.integrate import quad
from scipy.special import k0, k1

z, th, T = sp.symbols("z theta t")


# ------------------------------------------------------------------------------------------------ algebra
def _a_mul(P, Q, sgn=-1):
    """normal form u^k theta^j with theta u = u (theta + 2 sgn)."""
    out = {}
    for (a, b), c1 in P.items():
        for (c, d), c2 in Q.items():
            for (e,), cf in sp.Poly(sp.expand((th + 2 * sgn * c) ** b * th**d), th).terms():
                out[(a + c, e)] = out.get((a + c, e), 0) + c1 * c2 * cf
    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}


def _sym_op(m):
    """L_m in C[z]<theta>, z = r^2, monic of order m + 1 (Sym^m of theta^2 - z, cyclic vector K0^m)."""
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


PHI_THETA = {(0, 0): sp.Integer(-1), (0, 1): sp.Integer(-1)}   # phi(theta_r) = -1 - theta_t
PHI_Z = {(1, 2): sp.Integer(1)}                                 # phi(r^2) = t^-2 theta_t^2


def _phi(L):
    res = {}
    for (k, j), c in L.items():
        term = {(0, 0): sp.Integer(c)}
        for _ in range(k):
            term = _a_mul(term, PHI_Z)
        for _ in range(j):
            term = _a_mul(term, PHI_THETA)
        for kk, v in term.items():
            res[kk] = res.get(kk, 0) + v
    return {k: sp.expand(v) for k, v in res.items() if sp.expand(v) != 0}


def _P(m, lam):
    """phi(L_m(theta - lam))."""
    return _phi(_shift(_sym_op(m), lam))


def _scale(P, c):
    return {k: sp.expand(c * v) for k, v in P.items() if sp.expand(c * v) != 0}


def _equal(A, B):
    keys = set(A) | set(B)
    return all(sp.expand(A.get(k, 0) - B.get(k, 0)) == 0 for k in keys)


# ------------------------------------------------------------------------------------------------ Vanhove
# q_0, q_1, ..., q_{m-1} (leading coefficient last), arXiv:1401.6438 Table 1, equal masses, s = t^2.
VAN = {2: ["s-2", "s*(s-4)"],
       3: ["s-3", "3*s**2-20*s+9", "s*(s-1)*(s-9)"],
       4: ["s-4", "7*s**2-68*s+64", "6*s**3-90*s**2+192*s", "s**2*(s-4)*(s-16)"],
       5: ["s-5", "(3*s-5)*(5*s-57)", "25*s**3-518*s**2+1839*s-450",
           "10*s**4-280*s**3+1554*s**2-900*s", "s**2*(s-25)*(s-1)*(s-9)"]}
S_ = sp.Symbol("s")
THETA_S_PER_THETA_T = sp.Rational(1, 2)       # theta_s = theta_t / 2 for s = t^2


def _vanhove(m):
    """V_m as a u-theta dict: q_k(s) s^-k prod_{j<k} (theta_s - j), s^i = u^-i."""
    out = {}
    for k, qs in enumerate(VAN[m]):
        q = sp.Poly(sp.sympify(qs, locals={"s": S_}), S_)
        ff = sp.expand(sp.prod([THETA_S_PER_THETA_T * th - j for j in range(k)]))
        for (i,), c in q.terms():
            for (e,), cf in sp.Poly(ff, th).terms():
                out[(k - i, e)] = out.get((k - i, e), 0) + c * cf
    return {kk: sp.expand(v) for kk, v in out.items() if sp.expand(v) != 0}


# claimed data of the [SYMBOLIC] statement
U_OP = {(1, 2): sp.Integer(1)}       # t^-2 theta^2
A_M = {2: 0, 3: 1, 4: 1, 5: 2}
C_BASE = -2                 # constant (-2)^{m-1}
C_EXP_OFFSET = 1
LEFT_THETA_SHIFT = 0        # left factor t^{2a} (theta - LEFT_THETA_SHIFT)^2
ODD_LAM = 1


def _left_factor(m):
    sq = {(0, e): c for (e,), c in sp.Poly(sp.expand((th - LEFT_THETA_SHIFT) ** 2), th).terms()}
    return _a_mul({(-A_M[m], 0): sp.Integer(1)}, sq)


@pytest.mark.parametrize("m", [2, 3, 4, 5])
def test_odd_piece_is_theta2_times_vanhove(m):
    lhs = _P(m, ODD_LAM)
    c = sp.Integer(C_BASE) ** (m - C_EXP_OFFSET)
    # the paper's form is stated for the operator cleared of negative t-powers (the t-list form
    # sum_i t^{2i} Q_i(theta)), i.e. t^{2 K_m} phi(L_m(theta - 1)), K_m = largest u-power = floor((m+1)/2)
    Km = max(k for (k, _) in lhs)
    assert Km == (m + 1) // 2
    assert _equal(_a_mul({(-Km, 0): sp.Integer(1)}, lhs), _scale(_a_mul(_left_factor(m), _vanhove(m)), c))
    # equivalently, for the homomorphism phi itself: phi(L_m(theta - 1)) = (-2)^{m-1} U V_m, U = t^-2 theta^2,
    # the intertwiner of test 3 (a_m - K_m = -1 for every m)
    assert _equal(lhs, _scale(_a_mul(U_OP, _vanhove(m)), c))
    # order bookkeeping of the statement: V_m has order m - 1, the operator m + 1
    assert max(j for (_, j) in _vanhove(m)) == m - 1
    assert max(j for (_, j) in lhs) == m + 1


# ------------------------------------------------------------------------- right division in Q(t)<theta>
def _to_theta(P):
    """u-theta dict -> {j: c_j(t)} with the operator = sum_j c_j(t) theta^j (coefficients on the left)."""
    out = {}
    for (k, j), c in P.items():
        out[j] = out.get(j, 0) + c * T ** (-2 * k)
    return {j: sp.cancel(v) for j, v in out.items() if sp.cancel(v) != 0}


def _theta_times(coeffs, k):
    """theta^k * (sum_j c_j theta^j), using theta c = c theta + t c'."""
    cur = dict(coeffs)
    for _ in range(k):
        nxt = {}
        for j, c in cur.items():
            nxt[j + 1] = nxt.get(j + 1, 0) + c
            nxt[j] = nxt.get(j, 0) + T * sp.diff(c, T)
        cur = {j: sp.cancel(v) for j, v in nxt.items() if sp.cancel(v) != 0}
    return cur


def _rdiv_remainder(A, B):
    """remainder of right division A = Q B + R in Q(t)<theta> (A, B as {j: c_j(t)})."""
    A = dict(A)
    ob = max(B)
    while A and max(A) >= ob:
        oa = max(A)
        c = sp.cancel(A[oa] / B[ob])
        sub = _theta_times(B, oa - ob)
        for j, v in sub.items():
            A[j] = sp.cancel(A.get(j, 0) - c * v)
        A = {j: v for j, v in A.items() if v != 0}
        assert oa not in A
    return A


EVEN_LAM = 0


@pytest.mark.parametrize("m", [2, 3, 4, 5])
def test_vanhove_does_not_divide_even_piece(m):
    V = _to_theta(_vanhove(m))
    assert _rdiv_remainder(_to_theta(_P(m, ODD_LAM)), V) == {}      # positive control
    assert _rdiv_remainder(_to_theta(_P(m, EVEN_LAM)), V) != {}


# ------------------------------------------------------------------------------------------ intertwiner
INTERTWINE_STEP = 2


@pytest.mark.parametrize("m", [2, 3, 4, 5])
@pytest.mark.parametrize("lam", [-3, -1, 0, 1, 2, 3])
def test_parity_intertwiner(m, lam):
    L = _sym_op(m)
    lhs = _a_mul(_phi(_shift(L, lam + INTERTWINE_STEP)), U_OP)
    rhs = _a_mul(U_OP, _phi(_shift(L, lam)))
    assert _equal(lhs, rhs)
    # control: U does not commute with the operator, so the shift is doing the work
    assert not _equal(_a_mul(_phi(_shift(L, lam)), U_OP), rhs)


# ------------------------------------------------------------------------------------------ m <= 1
def M1_CLOSED(t):
    return 1 / (1 - t**2)


def test_m1_odd_moment_is_elementary():
    mp.mp.dps = 25
    for tv in (mp.mpf("0.3"), mp.mpf("0.75")):
        val = mp.quad(lambda x: x * mp.besseli(0, tv * x) * mp.besselk(0, x), [0, 1, 5, 20, mp.inf])
        assert abs(val - M1_CLOSED(tv)) < mp.mpf(10) ** -20 * abs(val)
    # phi(L_1(theta - 1)) kills the closed form exactly (order two, so its factors are of order one)
    f = M1_CLOSED(T)
    tot = 0
    for (k, j), c in _P(1, 1).items():
        g = f
        for _ in range(j):
            g = T * sp.diff(g, T)
        tot += c * T ** (-2 * k) * g
    assert sp.simplify(tot) == 0
    assert max(j for (_, j) in _P(1, 1)) == 2


# ------------------------------------------------------------------------------------------ numerics
_TH = {}


def _theta_I0(k):
    """theta^k I0(y) = a_k(y) I0(y) + b_k(y) y I1(y), exact polynomial a_k, b_k."""
    if k not in _TH:
        Y = sp.Symbol("Y")
        a, b = sp.Integer(1), sp.Integer(0)
        for _ in range(k):
            a, b = sp.expand(Y * sp.diff(a, Y) + Y**2 * b), sp.expand(a + Y * sp.diff(b, Y))
        _TH[k] = (sp.lambdify(Y, a, "mpmath"), sp.lambdify(Y, b, "mpmath"))
    return _TH[k]


def _theta_b(m, lam, tv, k):
    """theta_t^k b(t), b = int x^lam I0(t x) K0(x)^m dx, theta taken under the integral (exact in the kernel)."""
    a, b = _theta_I0(k)
    f = lambda x: x**lam * (a(tv * x) * mp.besseli(0, tv * x)  # noqa: E731
                           + b(tv * x) * tv * x * mp.besseli(1, tv * x)) * mp.besselk(0, x) ** m
    return mp.quad(f, [0, 1, 4, 12, 30, mp.inf])


def _apply(P, m, lam, tv):
    """(P b)(tv) and the largest term, P a u-theta dict."""
    K = max(j for (_, j) in P)
    ths = [_theta_b(m, lam, tv, k) for k in range(K + 1)]
    tot, scale = mp.mpf(0), mp.mpf(0)
    for (k, j), c in P.items():
        c = sp.Rational(c)
        term = mp.mpf(c.p) / c.q * tv ** (-2 * k) * ths[j]
        tot += term
        scale = max(scale, abs(term))
    return tot, scale


def VAN_CONST(m):
    return -sp.factorial(m) / sp.Integer(2) ** (m - 1)


_SLOW = pytest.mark.slow


@pytest.mark.parametrize("m", [2, 3, pytest.param(4, marks=_SLOW)])
def test_vanhove_constant(m):
    """V_m b_m = -m!/2^{m-1}, the same constant at t = 0.5 and 0.9."""
    mp.mp.dps = 20
    want = VAN_CONST(m)
    for tv in (mp.mpf("0.5"), mp.mpf("0.9")):
        val, _ = _apply(_vanhove(m), m, 1, tv)
        assert abs(val - mp.mpf(want.p) / want.q) < mp.mpf(10) ** -15


@pytest.mark.parametrize("m,tv", [(2, "0.6"), (2, "1.2"), (3, "0.6"),
                                  pytest.param(3, "2.2", marks=_SLOW), pytest.param(4, "0.6", marks=_SLOW),
                                  pytest.param(4, "3.2", marks=_SLOW)])
def test_odd_operator_annihilates_banana(m, tv):
    """phi(L_m(theta - 1)) b_m = 0 to 1e-17 of its largest term (dps 20)."""
    mp.mp.dps = 20
    val, scale = _apply(_P(m, ODD_LAM), m, 1, mp.mpf(tv))
    assert abs(val) < mp.mpf(10) ** -17 * scale


@pytest.mark.parametrize("m", [2, pytest.param(3, marks=_SLOW)])
def test_even_operator_does_not_annihilate_banana(m):
    """Control for the annihilation check: phi(L_m) leaves an O(1e-2) relative residual on the odd moment."""
    mp.mp.dps = 20
    val, scale = _apply(_P(m, EVEN_LAM), m, 1, mp.mpf("0.6"))
    assert abs(val) > mp.mpf(10) ** -3 * scale


# ------------------------------------------------------------------------------------------ part (i): parity
# The s-type leaf potential at the fingerprint weight W = K0/(2 pi^2), from the shell theorem
# Phi W(r) = 4 pi int W(s) s^2 / max(r, s) ds, splits as  TAIL_COEFF / r  -  (2/pi) R(r)/r,  R = int_r^oo s K1:
# an exterior tail r^{-1} (L = 0, twist -1: parity flipped, no Bessel factor) plus a Bessel part R/r whose
# twist is even, (theta + 1)(R/r) = theta K0 in B.  On the chain of three (root = centre, t on it),
# H = r^2 (Phi W)^2 = 1 - (4/pi) R + (4/pi^2) R^2: the R-linear piece comes with exactly ONE tail, and R itself
# has odd twist (theta R = r theta K0, in r B) -- lambda odd, m = 1 = n - 2.  The top piece R^2 has no tail.
TAIL_COEFF = 1.0


def _R(r):
    return quad(lambda s: s * k1(s), r, math.inf, epsabs=0.0, epsrel=1e-13, limit=200)[0]


def _leaf_potential(r):
    W = lambda s: k0(s) / (2 * math.pi**2)  # noqa: E731
    inner = quad(lambda s: W(s) * s**2, 0.0, r, epsabs=0.0, epsrel=1e-13, limit=200)[0] / r
    outer = quad(lambda s: W(s) * s, r, math.inf, epsabs=0.0, epsrel=1e-13, limit=200)[0]
    return 4 * math.pi * (inner + outer)


@pytest.mark.parametrize("r", [0.2, 1.0, 3.5])
def test_parity_leaf_potential_is_tail_plus_even_bessel_part(r):
    """Rejects a leaf potential without an exterior tail, a tail with a coefficient other than 1, and a Bessel
    part R (odd twist) in place of R/r (even twist); and an R-linear coefficient of H other than -4/pi (one tail)."""
    phi = _leaf_potential(r)
    assert phi == pytest.approx(TAIL_COEFF / r - (2 / math.pi) * _R(r) / r, rel=1e-11)
    # twists, by central differences: (theta + 1)(R/r) = theta K0 (even), theta R = r theta K0 (odd)
    h = 1e-4 * r
    dR = (_R(r + h) - _R(r - h)) / (2 * h)
    th_K0 = -r * k1(r)
    assert r * (dR / r - _R(r) / r**2) + _R(r) / r == pytest.approx(th_K0, rel=1e-7)
    assert r * dR == pytest.approx(r * th_K0, rel=1e-7)
    # H = r^2 (Phi W)^2 has R-linear coefficient 2 * TAIL * (-2/pi) = -4/pi
    R = _R(r)
    assert r**2 * phi**2 == pytest.approx(1 - (4 / math.pi) * R + (4 / math.pi**2) * R**2, rel=1e-10)


# The bookkeeping of part (i) as a rule: a piece of the s-type star of n densities (root carrying r^2, n - 1
# leaves) takes from each leaf either its Bessel part (twist BESSEL_TWIST, one K0) or its tail (twist TAIL_TWIST,
# no K0).  Part (ii)/(iii): a piece with m Bessel factors has, above order one, the factor of F_{m+1} (even lambda,
# m >= 1) or of V_m (odd lambda, order m - 1 > 1 only for m >= 3).
BESSEL_TWIST = 0
TAIL_TWIST = -1
ROOT_TWIST = 2


def _star_pieces(n):
    out = []
    for m in range(n):            # m leaves give their Bessel part, n - 1 - m their tail
        lam = ROOT_TWIST + m * BESSEL_TWIST + (n - 1 - m) * TAIL_TWIST
        out.append((m, lam))
    return out


def _factors_above_one(m, lam):
    if lam % 2 == 0:
        return {("F", m + 1)} if m >= 1 else set()
    return {("V", m)} if m - 1 > 1 else set()


def test_parity_bookkeeping_star_of_four_matches_measured_factors():
    """Prediction of parts (i)-(iii) for the four-density star vs its measured order-8 annihilator
    W8 = R4pp o X1 o E2 o (theta - 3) (tests/test_paper63_four_density_star.py): above order one exactly F_4's
    operator and the Legendre operator (F_2), no order-3 (K3, F_3) piece.  Rejects tails that keep parity
    (would predict an F_3 piece from the one-tail m = 2 term) and an odd Bessel twist (top piece odd: V_3 in
    place of F_4)."""
    predicted = set().union(*(_factors_above_one(m, lam) for m, lam in _star_pieces(4)))
    measured_w8 = {("F", 4), ("F", 2)}
    assert predicted == measured_w8


@pytest.mark.parametrize("n", range(3, 9))
def test_parity_bookkeeping_top_piece_even(n):
    """Every piece with m = n - 1 has even lambda; an odd piece has m <= n - 2 (one or more tails).  The chain
    of three's R-linear piece is the smallest odd one: m = 1, lambda = 1."""
    pieces = _star_pieces(n)
    assert all(lam % 2 == 0 for m, lam in pieces if m == n - 1)
    assert all(m <= n - 2 for m, lam in pieces if lam % 2)
    if n == 3:
        assert (1, 1) in pieces
