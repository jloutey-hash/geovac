"""Paper 59 sec:density_ring -- the chain at the Slater weights (eq:chain_slater_slice).

Record: docs/predictions/2026-10-19_p3e2_widened.md, P3e-3 (committed bb3ffe11 before computing).
Written as a separate step from the result.  Each test names the wrong answer it rejects.

  G(t) = int_0^inf x K1(t x) J(x)^2 dx,  J = 3 int_0^x K0 - 3 x K0 - x^2 K1,  theta = t d/dt.

  factorization     N2 o R3 o R1 o (theta + 2) equals the Mellin operator L_G exactly
  elementary pieces theta + 2 kills t^-2; R1 = (1 - t^2) theta - (1 + 7 t^2) kills t/(1 - t^2)^4
  K3 equivalence    R3 o A = 0 modulo F_3's operator 4 theta^3 - t^2 (theta + 1)^3, exactly, with A != 0
                    modulo it (so A carries F_3's solution space onto R3's)
  annihilation      W = R3 o R1 o (theta + 2) annihilates G, by quadrature of G and its theta-derivatives
                    in position space (independent of the Mellin route); L_G annihilates G too; a
                    one-part-in-198 change to one R3 coefficient does not; R1 (theta + 2) G != 0
  t = 1             the second solution of R1 o (theta + 2) has a pole of order 3 at t = 1, so the
                    operator's singularity there is not apparent

  order four        (record docs/predictions/2026-10-20_slice_artanh.md) R1 o S = 1 mod R3 has a unique
                    rational section within the exponent bound; W = Q1 o P with P of order four killing
                    t^-2; P annihilates G (so the artanh solution is absent); P has a log at t = 1 and a
                    half-integer exponent at t = 2, for the chain too, so analyticity does not force it
  parity proof      (record docs/predictions/2026-10-20_slice_parity_proof.md) R1 o P = -N o W exactly, so
                    P y is a multiple of the odd h = t/(1 - t^2)^4 for y in Sol(W); at infinity P keeps
                    exactly W's integer exponents; J and chi are odd at x = 0 (series), so G and C are
                    even at infinity and, given W G = 0 and W C = 0 (measured, not proved), P G = 0 and
                    P C = 0.  Slice and chain.

Not guarded: the no-split search (a mod-p search in the record) and the Mellin contour shift that
turns J's odd series into G's even expansion at infinity.  The paper states both caps.
"""
from __future__ import annotations

import mpmath as mp
import pytest
import sympy as sp

t, u, y = sp.symbols("t u y")

# The Mellin operator of G: L_G = sum_i t^(2i) Q_i(theta)  (record, step A).
QG = [
    -4 * (u - 5) * (u - 3) * (u - 1) ** 2 * (u + 1) * (u + 2) * (u + 5),
    (u - 3) * (u - 1) * (u + 1) * (u + 2) * (5 * u**3 + 88 * u**2 + 460 * u + 708),
    -(u - 1) * (u + 2) * (u**5 + 36 * u**4 + 463 * u**3 + 2718 * u**2 + 7396 * u + 7680),
    -6 * (u + 2) * (u + 4) * (u + 7) * (u + 8) * (u + 12),
]

# Factors, coefficient lists low -> high in theta (record, steps F and H).
TH2 = [sp.Integer(2), sp.Integer(1)]
R1 = [-(1 + 7 * t**2), 1 - t**2]
_d3 = (t - 2) * (t + 2) * (2 * t**2 + 7)
R3 = [60 * (4 * t**4 + 22 * t**2 + 7) / _d3, 4 * (46 * t**4 + 232 * t**2 + 91) / _d3,
      12 * (3 * t**4 + 11 * t**2 - 7) / _d3, sp.Integer(1)]
N2 = [(24 * t**8 + 12 * t**6 - 398 * t**4 - 1211 * t**2 - 980) / (2 * t**2 + 7) ** 2,
      6 * (3 * t**2 + 28) / (2 * t**2 + 7), (t - 2) * (t + 2)]
_dA = t**3 * (t**2 - 4) ** 5
A = [(4 * t**12 - 82 * t**10 + 745 * t**8 - 5924 * t**6 - 3328 * t**4 + 4416 * t**2 - 4096) / _dA,
     (8 * t**12 - 168 * t**10 + 1621 * t**8 - 14416 * t**6 - 18384 * t**4 + 16064 * t**2 - 10240) / _dA,
     2 * (2 * t**14 - 47 * t**12 + 517 * t**10 - 5012 * t**8 - 16400 * t**6 + 16160 * t**4
          - 21504 * t**2 + 16384) / (t**2 * _dA)]
L3 = [-t**2, -3 * t**2, -3 * t**2, 4 - t**2]          # 4 theta^3 - t^2 (theta + 1)^3


def _th(c):
    return sp.cancel(t * sp.diff(c, t))


def _compose(P, Q):
    """(P o Q) for P = sum p_k theta^k, Q = sum q_j theta^j."""
    out = [sp.Integer(0)] * (len(P) + len(Q) - 1)
    for k, p in enumerate(P):
        for j, q in enumerate(Q):
            qi = q
            for i in range(k + 1):
                out[k - i + j] += p * sp.binomial(k, i) * qi
                qi = _th(qi)
    return [sp.cancel(c) for c in out]


def _rrem(P, L):
    """Remainder of P on right division by L."""
    P = [sp.cancel(c) for c in P]
    m = len(L) - 1
    while len(P) - 1 >= m:
        n = len(P) - 1
        q = sp.cancel(P[n] / L[m])
        if q != 0:
            S = _compose([sp.Integer(0)] * (n - m) + [sp.Integer(1)], L)
            P = [sp.cancel(a - q * b) for a, b in zip(P, S)]
        P = P[:-1]
    return P


def _from_Q(Q):
    n = max(sp.Poly(q, u).degree() for q in Q) + 1
    out = [sp.Integer(0)] * n
    for i, q in enumerate(Q):
        for k, c in enumerate(sp.Poly(sp.expand(q), u).all_coeffs()[::-1]):
            out[k] += c * t ** (2 * i)
    return out


def _apply_sym(P, f):
    acc, g = sp.Integer(0), f
    for c in P:
        acc += c * g
        g = t * sp.diff(g, t)
    return sp.simplify(acc)


# ------------------------------------------------------------------ exact operator algebra
def test_factorization_is_exact():
    """Rejects: a wrong factor (any one), or a factorization that holds only up to a left multiplier."""
    W = _compose(N2, _compose(R3, _compose(R1, TH2)))
    LG = _from_Q(QG)
    assert len(W) == len(LG) == 8
    assert all(sp.cancel(a - b) == 0 for a, b in zip(W, LG))


def test_elementary_pieces():
    """Rejects: the printed solutions t^-2 and t/(1 - t^2)^4 not solving the printed factors."""
    assert _apply_sym(TH2, t**-2) == 0
    assert _apply_sym(R1, t / (1 - t**2) ** 4) == 0
    assert _apply_sym(R1, t / (1 - t**2) ** 3) != 0


def test_r3_is_equivalent_to_f3():
    """Rejects: R3 not equivalent to F_3's operator (R3 o A leaves a remainder mod L3), or a map A that
    kills F_3's solutions (A = 0 mod L3), which would make the remainder test empty."""
    assert all(c == 0 for c in _rrem(_compose(R3, A), L3))
    assert any(c != 0 for c in _rrem(A, L3))
    R3bad = [R3[0], R3[1], R3[2] * sp.Rational(199, 198), R3[3]]
    assert any(c != 0 for c in _rrem(_compose(R3bad, A), L3))


def test_t1_singularity_is_not_apparent():
    """Rejects: a removable singularity at t = 1 (the second solution of R1 o (theta + 2) regular there)."""
    s = sp.symbols("s")
    H = sp.integrate(s**2 / (1 - s**2) ** 4, s).subs(s, t) / t**2
    assert _apply_sym(_compose(R1, TH2), H) == 0
    lim = sp.limit((1 - t) ** 3 * H, t, 1)
    assert lim != 0 and lim.is_finite


# ------------------------------------------------------------------ numerical annihilation
def _theta_k1_powers(n):
    """theta_y^j K1(y) = a_j(y) K0(y) + b_j(y) K1(y), j = 0..n (theta K0 = -y K1, theta K1 = -y K0 - K1)."""
    a, b = sp.Integer(0), sp.Integer(1)
    out = [(a, b)]
    for _ in range(n):
        a, b = sp.expand(y * sp.diff(a, y) - y * b), sp.expand(y * sp.diff(b, y) - y * a - b)
        out.append((a, b))
    return [(sp.lambdify(y, a, "mpmath"), sp.lambdify(y, b, "mpmath")) for a, b in out]


def _J(x):
    # int_0^x K0 = (pi x / 2) (K0 L_{-1} + K1 L_0), Abramowitz-Stegun 11.1.8
    k0, k1 = mp.besselk(0, x), mp.besselk(1, x)
    intk0 = mp.pi * x / 2 * (k0 * mp.struvel(-1, x) + k1 * mp.struvel(0, x))
    return 3 * intk0 - 3 * x * k0 - x * x * k1


def _moments(tv, n):
    """I_j = theta_t^j G(t), j = 0..n, by quadrature (theta_t^j K1(t x) = (theta_y^j K1)(t x))."""
    pw = _theta_k1_powers(n)
    cache = {}

    def J2(x):
        if x not in cache:
            cache[x] = _J(x) ** 2
        return cache[x]

    out = []
    for a, b in pw:
        def h(x, a=a, b=b):
            yy = tv * x
            return x * J2(x) * (a(yy) * mp.besselk(0, yy) + b(yy) * mp.besselk(1, yy))
        out.append(mp.quad(h, [0, mp.mpf(1) / 8, 1, 4, 16, 64, mp.inf]))
    return out


def _apply_num(P, I, tv):
    terms = [mp.mpf(sp.N(c.subs(t, sp.Rational(tv)), mp.mp.dps + 5)) * I[k] for k, c in enumerate(P)]
    return mp.fsum(terms), max(abs(x) for x in terms)


def test_integrated_k_zero_closed_form():
    """Rejects: a wrong closed form for int_0^x K0, on which every quadrature below rests."""
    mp.mp.dps = 30
    for x in (mp.mpf("0.3"), mp.mpf(2), mp.mpf(7)):
        ref = mp.quad(lambda s: mp.besselk(0, s), [0, x])
        k0, k1 = mp.besselk(0, x), mp.besselk(1, x)
        assert abs(mp.pi * x / 2 * (k0 * mp.struvel(-1, x) + k1 * mp.struvel(0, x)) - ref) < mp.mpf(10) ** -25


G_RECORD = {"1/2": "63.71736086669894693826485", "13/10": "3.646420116881517536277732"}


@pytest.mark.slow
@pytest.mark.parametrize("tv", ["1/2", "13/10"])
def test_w_annihilates_g(tv):
    """Rejects: a wrong G (value against the record), G needing N2 (W G != 0), a Mellin operator that does not annihilate G, a guard too coarse to
    see a 1/198 change in R3, and a G without a K3 piece (R1 (theta + 2) G = 0)."""
    mp.mp.dps = 25
    tq = sp.Rational(tv)
    I = _moments(mp.mpf(tq.p) / tq.q, 7)
    # the record's value, from a different evaluator of int_0^x K0 (step A, termwise quadrature)
    assert abs(I[0] - mp.mpf(G_RECORD[tv])) < mp.mpf(10) ** -18
    W = _compose(R3, _compose(R1, TH2))
    r, s = _apply_num(W, I, tq)
    assert abs(r) / s < mp.mpf(10) ** -18
    rl, sl = _apply_num(_from_Q(QG), I, tq)
    assert abs(rl) / sl < mp.mpf(10) ** -18
    R3bad = [R3[0], R3[1], R3[2] * sp.Rational(199, 198), R3[3]]
    rb, sb = _apply_num(_compose(R3bad, _compose(R1, TH2)), I, tq)
    assert abs(rb) / sb > mp.mpf(10) ** -6
    e2, se = _apply_num(_compose(R1, TH2), I, tq)
    assert abs(e2) / se > mp.mpf(10) ** -3


# ------------------------------------------------------------------ the order-four factor (artanh absent)
# Record docs/predictions/2026-10-20_slice_artanh.md (SA1, SA2; committed cffafc1d before computing).
# Written as a separate step from the result.  Work in x = t^2 (the paper's w), delta = x d/dx (theta = 2 delta), with
# the chain guard's exact delta-form machinery; operators are lists [c_0, c_1, ...], functions on the left.
import test_paper59_density_ring as DR  # noqa: E402

X_ = DR.x_
Rat = sp.Rational
# the partial section S (R1 o S = 1 mod R3), on (1, delta, delta^2)
_dS = X_ * (X_ - 1) ** 4 * (2 * X_ + 7)
S_DELTA = [-(378 * X_**5 - 327 * X_**4 + 635 * X_**3 - 605 * X_**2 + 475 * X_ - 220) / (105 * _dS),
           -2 * (546 * X_**5 - 339 * X_**4 + 1309 * X_**3 - 2041 * X_**2 + 1893 * X_ - 792) / (315 * _dS),
           -4 * (X_ - 4) * (42 * X_**4 - 45 * X_**3 + 77 * X_**2 - 75 * X_ + 33) / (315 * _dS)]
# W = Q1 o P, P = (theta + 2) - S o R1 o (theta + 2)
_dQ = (X_ - 4) * (42 * X_**4 - 45 * X_**3 + 77 * X_**2 - 75 * X_ + 33)
Q1_SLATER = [315 * X_ * (X_ - 1) ** 3 * (2 * X_ + 7) * (7 * X_ + 1) / _dQ,
             630 * X_ * (X_ - 1) ** 4 * (2 * X_ + 7) / _dQ]


def _slater_delta():
    """R3, R1, theta + 2 and E = R1 o (theta + 2) of the slice, in delta form."""
    sub = {t: DR.t}
    R3d = DR._theta_to_delta([sp.sympify(c).subs(sub) for c in R3])
    R1d = DR._theta_to_delta([sp.sympify(c).subs(sub) for c in R1])
    T2d = [sp.Integer(2), sp.Integer(2)]
    return R3d, R1d, T2d, DR._dcompose(R1d, T2d)


def _P(S, T2d, Ed):
    return [sp.cancel(a - b) for a, b in zip(T2d + [0, 0, 0], DR._dcompose(S, Ed))]


def test_partial_section_and_order_four_factor():
    """Rejects a wrong section (s0 off by 1%), a W that is not a left multiple of P, a wrong Q1, and a P
    that does not kill t^-2 = x^-1 (delta^k x^-1 = (-1)^k x^-1)."""
    R3d, R1d, T2d, Ed = _slater_delta()
    assert DR._is_one_mod(R1d, S_DELTA, R3d)
    assert not DR._is_one_mod(R1d, [S_DELTA[0] * Rat(101, 100)] + S_DELTA[1:], R3d)
    P = _P(S_DELTA, T2d, Ed)
    assert len(P) == 5 and P[4] != 0
    Q, rem = DR._ddiv(DR._dcompose(R3d, Ed), P)
    assert all(c == 0 for c in rem)
    assert [sp.cancel(a - b) for a, b in zip(Q, Q1_SLATER)] == [0, 0]
    assert sp.cancel(sum(c * (-1) ** k for k, c in enumerate(P))) == 0


@pytest.mark.slow
def test_section_is_unique_within_the_bound():
    """Rejects 'no section exists' (SA1's falsifier), a section that is not unique, and a shrunken bound:
    the exponents at x = 0 and x = 1 must be the hand values, the exact system over Q must have 57
    unknowns at full rank, and its solution must be S_DELTA."""
    R3d, R1d, _, _ = _slater_delta()
    table = DR._split_bounds(R3d, R1d)
    assert list(table[sp.Integer(0)][1]) == [Rat(-5, 2), Rat(-1, 2), Rat(3, 2)]
    assert list(table[sp.Integer(1)][0]) == [-4, 0, 1, 2] and table[sp.Integer(1)][1] == "ord"
    nunk, rA, rAb, S = DR._split_search(R3d, R1d, table)
    assert (nunk, rA, rAb) == (57, 57, 57)
    assert [sp.cancel(a - b) for a, b in zip(S, S_DELTA)] == [0, 0, 0]


_r, _lg, _z = sp.symbols("r_ Lg z_")


def _layers(L, p):
    """z^n L = sum_i z^i P_i(theta_z) at x = p + z (d-form a_k d^k = a_k z^-k [theta_z]_k)."""
    a = DR._d_form(L)
    n = len(a) - 1
    terms = [(sp.cancel(c.subs(X_, p + _z) * _z ** (n - k)), sp.prod([_r - j for j in range(k)]))
             for k, c in enumerate(a)]
    den = sp.lcm([sp.denom(sp.together(c)) for c, _ in terms])
    out = {}
    for c, ff in terms:
        pc = sp.Poly(sp.expand(sp.cancel(c * den)), _z)
        for (e,), v in zip(pc.monoms(), pc.coeffs()):
            out[e] = out.get(e, 0) + v * ff
    lo = min(out)
    return {i - lo: sp.expand(q) for i, q in out.items()}


def _max_log(L, p, K=3):
    """Largest power of log(x - p) in the local solutions of L at x = p whose exponents are integers
    (generalised Frobenius y = sum_n z^(rho0 + n) q_n(log z), exact over Q)."""
    ex = [e for e in DR._exps_finite(L, X_ - p)[0] if e.is_integer]
    rho0, nmax = min(ex), int(max(ex) - min(ex))
    Pl = _layers(L, p)
    unk, ys = [], []
    for n in range(nmax + 1):
        cs = sp.symbols(f"a{n}_0:{K + 1}")
        unk += cs
        ys.append(sum(c * _lg**j for j, c in enumerate(cs)))
    eqs = []
    for n in range(nmax + 1):
        tot = 0
        for i, Pi in Pl.items():
            if n - i >= 0:
                sh = sp.Poly(sp.expand(Pi.subs(_r, _r + rho0 + n - i)), _r)
                tot += sum(c * sp.diff(ys[n - i], _lg, k) for (k,), c in zip(sh.monoms(), sh.coeffs()))
        tot = sp.expand(tot)
        eqs += [tot.coeff(_lg, j) for j in range(K + 1)]
    M = sp.Matrix([[sp.diff(e, v) for v in unk] for e in eqs])
    ns = M.nullspace()
    assert len(ns) == len(ex)
    deg = 0
    for v in ns:
        sub = dict(zip(unk, v))
        for y in ys:
            q = sp.expand(y.subs(sub))
            if q != 0:
                deg = max(deg, sp.Poly(q, _lg).degree())
    return deg


def _check_not_forced(P, Ed, R3d, half4):
    assert _max_log(Ed, 1) == 1 and _max_log(R3d, 1) == 0     # detector controls
    assert _max_log(P, 1) == 1
    e4 = DR._exps_finite(P, X_ - 4)[0]
    assert [e for e in e4 if not e.is_integer] == [half4]


def test_membership_is_not_forced_by_analyticity():
    """Rejects 'P is log-free at t = 1' and 'P has only integer exponents at t = 2': either would make
    ker P the solutions single-valued there, and analyticity of G on Re t > 0 alone would prove P G = 0.
    The paper says it does not, and derives P G = 0 by parity at infinity from the measured W G = 0
    instead (tests below).
    Detector controls: E = R1 (theta + 2) must show its artanh log at x = 1, and R3, ordinary there, none."""
    R3d, _, T2d, Ed = _slater_delta()
    _check_not_forced(_P(S_DELTA, T2d, Ed), Ed, R3d, Rat(-9, 2))


@pytest.mark.slow
def test_chain_membership_is_not_forced_by_analyticity():
    """The same for the chain's P: a log at t = 1 and the
    exponent -3/2 at t = 2."""
    R3d, _, T2d, Ed, _ = DR._chain_delta(DR._chain_pieces()[2])
    P = [sp.cancel(a - b) for a, b in zip(T2d + [0, 0, 0], DR._dcompose(DR.S1_DELTA, Ed))]
    _check_not_forced(P, Ed, R3d, Rat(-3, 2))


@pytest.mark.slow
@pytest.mark.parametrize("tv", ["1/2", "13/10"])
def test_g_lies_on_the_order_four_factor(tv):
    """Rejects 'G needs the artanh solution of R1 (theta + 2)' (minimal order five): P G must be below
    1e-22 of its largest term at 25 digits, while P from S with s0 times 1.01 must give more than 1e-5."""
    mp.mp.dps = 25
    _, _, T2d, Ed = _slater_delta()
    back = {DR.t: t}

    def theta_form(S):
        return [sp.sympify(c).subs(back) for c in DR._delta_to_theta(_P(S, T2d, Ed))]
    tq = sp.Rational(tv)
    I = _moments(mp.mpf(tq.p) / tq.q, 4)
    assert abs(I[0] - mp.mpf(G_RECORD[tv])) < mp.mpf(10) ** -18
    r, s = _apply_num(theta_form(S_DELTA), I, tq)
    assert abs(r) / s < mp.mpf(10) ** -22
    rb, sb = _apply_num(theta_form([S_DELTA[0] * Rat(101, 100)] + S_DELTA[1:]), I, tq)
    assert abs(rb) / sb > mp.mpf(10) ** -5


# ------------------------------------------------------------------ the parity proof that G and C lie on P
# Record docs/predictions/2026-10-20_slice_parity_proof.md (MP3; committed 5ec2b08d before computing).
# Written as a separate step from the result.
N_SLATER = (X_ - 4) * (42 * X_**4 - 45 * X_**3 + 77 * X_**2 - 75 * X_ + 33) / (315 * X_ * (X_ - 1) ** 3 * (2 * X_ + 7))
N_CHAIN = (X_ - 4) * (20 * X_**2 - 7 * X_ + 5) / (15 * X_ * (X_ - 1) * (4 * X_ + 11))


def _check_left_multiple(R3d, R1d, T2d, Ed, S, Nx):
    W = DR._dcompose(R3d, Ed)
    Q, rem = DR._ddiv(DR._dcompose(R1d, _P(S, T2d, Ed)), W)
    assert all(c == 0 for c in rem)
    assert len(Q) == 1 and sp.cancel(Q[0] + Nx) == 0
    _, remb = DR._ddiv(DR._dcompose(R1d, _P([S[0] * Rat(101, 100)] + list(S[1:]), T2d, Ed)), W)
    assert any(c != 0 for c in remb)


def test_r1_p_is_a_left_multiple_of_w():
    """Rejects 'P y need not solve R1' (the step that makes P y = alpha(y) h): R1 o P = -N o W with zero
    remainder and the recorded N, while P from s0 times 1.01 leaves a remainder."""
    R3d, R1d, T2d, Ed = _slater_delta()
    _check_left_multiple(R3d, R1d, T2d, Ed, S_DELTA, N_SLATER)
    # h = t/(1 - t^2)^4, the solution of R1, is odd in t: h(-t) = -h(t)
    h = t / (1 - t**2) ** 4
    assert _apply_sym(R1, h) == 0 and sp.cancel(h.subs(t, -t) + h) == 0


@pytest.mark.slow
def test_chain_r1_p_is_a_left_multiple_of_w():
    """The same for the chain, with h = t^3/(1 - t^2)^2, odd in t."""
    R3d, R1d, T2d, Ed, _ = DR._chain_delta(DR._chain_pieces()[2])
    _check_left_multiple(R3d, R1d, T2d, Ed, DR.S1_DELTA, N_CHAIN)
    h = DR.t**3 / (1 - DR.t**2) ** 2
    R1c = DR._delta_to_theta(R1d)
    assert sp.simplify(sum(c * hk for c, hk in zip(R1c, [h, DR.t * sp.diff(h, DR.t)]))) == 0


def test_p_keeps_the_integer_exponents_at_infinity():
    """Rejects a P that keeps the half-integer solution at infinity (then an even G could still have
    P G != 0): W at infinity is {1, 2, 7/2, 4, 6} in w = 1/x and P is {1, 2, 4, 6}."""
    R3d, _, T2d, Ed = _slater_delta()
    assert DR._exps_inf(DR._dcompose(R3d, Ed)) == [1, 2, Rat(7, 2), 4, 6]
    assert DR._exps_inf(_P(S_DELTA, T2d, Ed)) == [1, 2, 4, 6]


@pytest.mark.slow
def test_chain_p_keeps_the_integer_exponents_at_infinity():
    """The same for the chain: W {-1, 0, 1/2, 1, 2}, P {-1, 0, 1, 2}."""
    R3d, _, T2d, Ed, _ = DR._chain_delta(DR._chain_pieces()[2])
    assert DR._exps_inf(DR._dcompose(R3d, Ed)) == [-1, 0, Rat(1, 2), 1, 2]
    P = [sp.cancel(a - b) for a, b in zip(T2d + [0, 0, 0], DR._dcompose(DR.S1_DELTA, Ed))]
    assert DR._exps_inf(P) == [-1, 0, 1, 2]


_zs, _L, _c = sp.symbols("zs Ls cs")   # Ls = log z, cs = gamma - log 2


def _odd_series(NK=6):
    """J and chi near 0 from K0 = -(log(z/2) + gamma) I0 + sum_k H_k (z/2)^(2k)/(k!)^2 (A&S 9.6.13),
    K1 = -K0', termwise integration; returns (J, chi) as polynomials in z, Ls, cs."""
    I0 = sum((_zs / 2) ** (2 * k) / sp.factorial(k) ** 2 for k in range(NK))
    K0 = sp.expand(-(_L + _c) * I0 + sum(sp.harmonic(k) * (_zs / 2) ** (2 * k) / sp.factorial(k) ** 2
                                         for k in range(1, NK)))

    def d(e):
        return sp.expand(sp.diff(e, _zs) + sp.diff(e, _L) / _zs)

    def integ(e):
        out = 0
        for term in sp.Add.make_args(sp.expand(e)):
            n, j = sp.degree(term, _zs), sp.degree(term, _L)
            a = term / (_zs**n * _L**j)
            out += (a * (_zs ** (n + 1) * _L / (n + 1) - _zs ** (n + 1) / (n + 1) ** 2) if j == 1
                    else a * _zs ** (n + 1) / (n + 1))
        return sp.expand(out)
    K1 = sp.expand(-d(K0))
    J = sp.expand(3 * integ(K0) - 3 * _zs * K0 - _zs**2 * K1)
    chi = integ(sp.expand(_zs * K1))
    return J, chi


def _powers(e, top):
    out = set()
    for term in sp.Add.make_args(e):
        num, den = sp.fraction(sp.together(term))
        n = sp.degree(num, _zs) - sp.degree(den, _zs)
        if n < top:
            out.add(int(n))
    return sorted(out)


def test_j_and_chi_are_odd_at_zero():
    """Rejects an even power of x in J or chi near 0 (which would put odd powers of t in G's or C's
    expansion at infinity and break the parity step).  Control: the series reproduce J and chi at 0.3."""
    J, chi = _odd_series()
    assert _powers(J, 11) == [1, 3, 5, 7, 9]
    assert _powers(chi, 11) == [1, 3, 5, 7, 9]
    mp.mp.dps = 30
    zz = mp.mpf("0.3")
    sub = {_zs: zz, _L: mp.log(zz), _c: mp.euler - mp.log(2)}
    assert abs(mp.mpf(str(sp.N(J.subs(sub), 25))) - _J(zz)) < mp.mpf(10) ** -12
    chiref = mp.quad(lambda s: s * mp.besselk(1, s), [0, zz])
    assert abs(mp.mpf(str(sp.N(chi.subs(sub), 25))) - chiref) < mp.mpf(10) ** -12
