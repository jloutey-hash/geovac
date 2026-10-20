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

Not guarded: the no-split search (a mod-p search in the record).  The paper states it with that cap.
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
