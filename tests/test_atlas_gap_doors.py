"""Guard for three proposed atlas rows (docs/atlas.md doors 28-30, provisional numbering; derivations in the
2026-11-19 atlas-gap memo).  The computations are copied into tests/_atlas_gap_doors.py.

Door 28 (a uniform electric field on hydrogen, Stark).  The ground-state Rayleigh-Schroedinger coefficients are
exact rationals, odd orders vanish:
    E = -1/2 - (9/4) F^2 - (3555/64) F^4 - (2512779/512) F^6 - (13012777803/16384) F^8 - ...
The series diverges as a_k ~ -(6/pi)(3/2)^{2k}(2k)! (coefficient of F^{2k}; Oppenheimer width through the
dispersion relation).  First order: in every m block of shells n = 2, 3 the eigenvalues of z are
(3/2) n (n1 - n2), n1 + n2 = n - |m| - 1; the n = 2 shell splits 4 -> 1 + 2 + 1.
Door 29 (a uniform magnetic field, diamagnetic term (1/8) B^2 rho^2, beta = B^2):
    E = -1/2 + (1/4) beta - (53/192) beta^2 + (5581/4608) beta^3 - (21577397/2211840) beta^4 + ...
with large order a_k ~ (-1)^{k+1} (4/pi)^{5/2} pi^{-2k} Gamma(2k + 3/2).
Door 30 (lens spaces L(p,1) = S^3/Z_p): the multiplicity of n^2 - 1 is the number of the n^2 parabolic states
with L_z + A_z = 0 mod p; table for p = 2, 3, 5, n = 1..8 by three routes (weight count, character average,
Molien series) and the memo's generating functions; F_3 = sum m_n/n^4 = pi^2/18 + 8 sqrt3 pi^3/729;
sum_{k=1}^{p-1} cot^2(pi k/p) = (p-1)(p-2)/3, p = 2..11.

Each test's docstring names the wrong answer it rejects.  Slow tier: the two large-order constants (Richardson
on exact coefficients to F^60 and beta^30).
"""
from __future__ import annotations

from fractions import Fraction as Fr

import mpmath as mp
import pytest
import sympy as sp

import _atlas_gap_doors as H

STARK_EXPECTED = [Fr(-1, 2), Fr(0), Fr(-9, 4), Fr(0), Fr(-3555, 64), Fr(0), Fr(-2512779, 512), Fr(0),
                  Fr(-13012777803, 16384)]
DIAMAG_EXPECTED = [Fr(-1, 2), Fr(1, 4), Fr(-53, 192), Fr(5581, 4608), Fr(-21577397, 2211840)]

LENS_TABLE = {
    2: [1, 0, 9, 0, 25, 0, 49, 0],
    3: [1, 0, 3, 8, 5, 12, 21, 16],
    5: [1, 0, 3, 0, 5, 12, 7, 16],
}

# Richardson tolerance (relative).  Measured at k = 30 (F^60) and beta^30, order 8: Stark 1.5e-6,
# diamagnetic 2.2e-6 (to F^120 / beta^60 the memo reaches 8e-11 / 1e-10).  1e-5 leaves a 4-7x margin and
# still separates 6/pi from any nearby simple constant (2 is 4.7% away; 3/pi a factor 2).
RICH_TOL = mp.mpf("1e-5")


# ---------------------------------------------------------------------------------------------- door 28

def test_door28_stark_rs_coefficients_exact():
    """Rejects -3555/32 (or any wrong F^4 coefficient), a positive F^2 term (+9/4: the field would raise the
    ground state), -2512779/1024, and a non-vanishing odd order (parity broken)."""
    E = H.series(H.V_stark, 16)
    assert E[:9] == STARK_EXPECTED
    assert all(E[k] == 0 for k in range(1, 17, 2)), "odd orders must vanish (z is odd under parity)"
    # every even order is strictly negative (the field lowers the ground state at every order)
    assert all(E[k] < 0 for k in range(2, 17, 2))


def test_door28_linear_stark_parabolic_shifts():
    """Rejects shifts (3/2) n^2 (n1 - n2), 3 n (n1 - n2), or any l-diagonal first-order shift (degenerate
    perturbation theory skipped), and the n = 2 split 1 + 3 (or 2 + 2) instead of 1 + 2 + 1."""
    x = sp.symbols("x")
    n2_eigs = []
    for n in (2, 3):
        for m in range(-(n - 1), n):
            M = H.z_block(n, m)
            pred = H.parabolic_shifts(n, m)
            cp = sp.expand(M.charpoly(x).as_expr())
            assert sp.expand(cp - sp.prod([x - e for e in pred])) == 0, (n, m, pred)
            if n == 2:
                n2_eigs.extend(pred)
    vals = sorted(set(n2_eigs))
    assert [n2_eigs.count(v) for v in vals] == [1, 2, 1]
    assert vals == [-3, 0, 3]


@pytest.mark.slow
def test_door28_stark_large_order_constant():
    """Rejects a large-order constant 3/pi (ratio -> 2), 2 (ratio -> 1.047), 12/pi, or a wrong action
    (a (3/2)^{2k} replaced by 2^{2k} would send the ratio to 0 geometrically).  Richardson, orders 2..8, on
    exact coefficients to F^60 (k = 30)."""
    with mp.workdps(50):
        E = H.series(H.V_stark, 60)
        ks, xs = H.stark_ratios(E, 30)
        r2 = H.richardson(xs, ks, 2)
        r8 = H.richardson(xs, ks, 8)
        assert abs(r8 - 1) < RICH_TOL, mp.nstr(r8, 15)
        assert abs(r8 - 1) < abs(r2 - 1), "higher Richardson order must improve the estimate"
        assert abs(xs[-1] - 1) > 1e-3, "the raw ratio is still visibly away from 1; Richardson is doing the work"


# ---------------------------------------------------------------------------------------------- door 29

def test_door29_diamagnetic_coefficients_exact():
    """Rejects -53/96 at beta^2 (or any wrong coefficient), 1/8 at beta (the <rho^2> = 2 average dropped),
    and a non-alternating series (5581/4608 with the wrong sign)."""
    E = H.series(H.V_diamag, 4)
    assert E[:5] == DIAMAG_EXPECTED
    assert all((E[k] > 0) == (k % 2 == 1) for k in range(1, 5)), "the series alternates"


@pytest.mark.slow
def test_door29_diamagnetic_large_order_constant():
    """Rejects a large-order constant (4/pi)^2 or 4/pi^{5/2}, the sign (-1)^k instead of (-1)^{k+1} (ratio
    -> -C), and a growth pi^{-k} (ratio diverges).  Richardson order 8 on exact coefficients to beta^30."""
    with mp.workdps(50):
        E = H.series(H.V_diamag, 30)
        ks, xs = H.diamag_ratios(E, 30)
        C = (4 / mp.pi) ** (mp.mpf(5) / 2)
        r2 = H.richardson(xs, ks, 2)
        r8 = H.richardson(xs, ks, 8)
        assert abs(r8 / C - 1) < RICH_TOL, mp.nstr(r8, 15)
        assert abs(r8 / C - 1) < abs(r2 / C - 1)
        assert abs(xs[-1] / C - 1) > 1e-3


# ---------------------------------------------------------------------------------------------- door 30

@pytest.mark.parametrize("p", [2, 3, 5])
def test_door30_lens_multiplicities_two_routes(p):
    """Rejects the p = 3 row shifted by one n (0, 3, 8, 5, ...), multiplicity n^2/p (fractional), the
    condition L_z = 0 mod p instead of L_z + A_z (the p = 2 row would read 1, 2, 5, 8, ...), and the
    character average dropping the n from the SU(2)_R factor."""
    w = [H.mult_weights(p, n) for n in range(1, 9)]
    c = [H.mult_character(p, n) for n in range(1, 9)]
    assert w == LENS_TABLE[p]
    assert c == LENS_TABLE[p]


@pytest.mark.parametrize("p", [2, 3, 5])
def test_door30_lens_generating_function(p):
    """Rejects a generating function that matches only the first few terms (checked to t^30 against the
    exact Molien series and the weight count), and the p = 2 numerator 1 + 3 t^2 + t^4."""
    order = 30
    gen = H.series_coeffs(H.GEN_FUNCS[p], order)
    assert gen == H.molien_coeffs(p, order)
    assert gen == [H.mult_weights(p, n) for n in range(1, order + 2)]
    if p == 2:
        t = H.T
        assert sp.simplify(H.GEN_FUNCS[2] * (1 - t ** 2) ** 3 - (1 + 6 * t ** 2 + t ** 4)) == 0
        # the Molien series times (1 - t^2)^3 is the polynomial 1 + 6 t^2 + t^4, exactly, to order 30
        mc = H.molien_coeffs(2, order)
        prod = sp.expand(sum(mc[k] * t ** k for k in range(order + 1)) * (1 - t ** 2) ** 3)
        assert [prod.coeff(t, k) for k in range(order + 1)] == [1, 0, 6, 0, 1] + [0] * (order - 4)


def test_door30_F3_closed_form():
    """Rejects F_3 without the pi^3 term (pi^2/18 alone), with 16 sqrt3 pi^3/729, or with an added zeta(3)
    piece (the even part of e(n) must cancel).  Route: Hurwitz zeta on residue classes mod 6 at 70 digits,
    tolerance 1e-60; the class-linearity it relies on is checked to n = 600, and a plain partial sum to
    n = 4000 with a tail bound confirms the value independently at 1e-3."""
    p = 3
    q = 2 * p
    for a in range(1, q + 1):
        c1, c2 = H.c_count(p, a), H.c_count(p, a + q)
        for j in range(2, 100):
            assert H.c_count(p, a + j * q) == c1 + j * (c2 - c1), (a, j)
    with mp.workdps(70):
        closed = mp.pi ** 2 / 18 + 8 * mp.sqrt(3) * mp.pi ** 3 / 729
        assert abs(H.F_lens_numeric(3) - closed) < mp.mpf("1e-60")
        assert abs(H.F_lens_numeric(2) - mp.pi ** 2 / 8) < mp.mpf("1e-60")
    N = 4000
    partial = sum(Fr(H.c_count(3, n), n ** 3) for n in range(1, N + 1))
    # c(n) <= n/3 + 1, so the tail is below sum_{n>N} (1/(3n^2) + 1/n^3) < 1/(3N) + 1/(2N^2) < 1e-4
    assert 0 <= float(closed) - float(partial) < 1e-3


def test_door30_cotangent_sum():
    """Rejects (p-1)(p-2)/6, (p^2-1)/3 (the cosec^2 sum), and a non-rational value: exact route (Newton's
    identities on the polynomial whose roots are cot(pi k/p), roots checked) plus a 50-digit direct sum."""
    with mp.workdps(50):
        for p in range(2, 12):
            target = Fr((p - 1) * (p - 2), 3)
            assert H.cot_poly_power_sum(p) == target
            roots = [mp.cot(mp.pi * k / p) for k in range(1, p)]
            assert all(abs(H.cot_poly_value(p, x)) < mp.mpf("1e-40") * max(1, abs(x)) ** (p - 1) for x in roots)
            assert len({mp.nstr(x, 30) for x in roots}) == p - 1, "the p - 1 roots are distinct"
            assert abs(sum(x ** 2 for x in roots) - mp.mpf(target.numerator) / target.denominator) < mp.mpf("1e-40")


# Fire test (2026-11-19).  Each plant went into a fresh scratch copy of this file and its helper (never the
# tree), run with --slow and -k on the named test; FIRED = pytest exit 1, collection clean.  T = this file,
# H = tests/_atlas_gap_doors.py.
#
#   plant                                                        test                         result
#   F^4 expected -3555/32 (T)                                    stark_rs_coefficients        FIRED
#   F^2 expected +9/4 (T)                                        stark_rs_coefficients        FIRED
#   u P_j recursion (j+1)/(2j+1) -> (j+2)/(2j+1) (H)             stark_rs_coefficients        FIRED
#   Stark perturbation r^2 cos(theta) instead of r cos(theta) (H) stark_rs_coefficients       FIRED
#   shifts (3/2) n^2 (n1 - n2) (H)                               linear_stark                 FIRED
#   shifts 3 n (n1 - n2) (H)                                     linear_stark                 FIRED
#   n = 2 split expected 1 + 3 (T)                               linear_stark                 FIRED
#   wrong angular factor in the z matrix (H)                     linear_stark                 FIRED
#   large-order constant 3/pi (H)                                stark_large_order (slow)     FIRED
#   large-order constant 2 (H)                                   stark_large_order (slow)     FIRED
#   action (1.49)^{2k} instead of (3/2)^{2k} (H)                 stark_large_order (slow)     FIRED
#   beta^2 expected -53/96 (T)                                   diamagnetic_coefficients     FIRED
#   diamagnetic term (1/4) rho^2 instead of (1/8) (H)            diamagnetic_coefficients     FIRED
#   sign (-1)^k instead of (-1)^{k+1} (H)                        diamag_large_order (slow)    FIRED
#   constant (4/pi)^2 instead of (4/pi)^{5/2} (T)                diamag_large_order (slow)    FIRED
#   Gamma(2k + 1) instead of Gamma(2k + 3/2) (H)                 diamag_large_order (slow)    FIRED
#   p = 3 row shifted by one n (T)                               lens_multiplicities          FIRED 1/3 (p = 3)
#   weight condition L_z = 0 mod p (A_z dropped) (H)             lens_multiplicities          FIRED 3/3
#   character average without the SU(2)_R factor n (H)           lens_multiplicities          FIRED 3/3
#   p = 5 numerator 13 t^4 -> 12 t^4 (H)                         lens_generating_function     FIRED 1/3 (p = 5)
#   p = 2 numerator 1 + 3 t^2 + t^4 (H)                          lens_generating_function     FIRED 1/3 (p = 2)
#   Molien series without the (1 - t^2) harmonic factor (H)      lens_generating_function     FIRED 3/3
#   F_3 without the pi^3 term (T)                                F3_closed_form               FIRED
#   F_3 with 16 sqrt3 pi^3/729 (T)                               F3_closed_form               FIRED
#   F_3 plus 1e-40 zeta(3) (T)                                   F3_closed_form               FIRED
#   cotangent sum (p-1)(p-2)/6 (T)                               cotangent_sum                FIRED
#   cotangent sum (p^2-1)/3, the cosec^2 sum (T)                 cotangent_sum                FIRED
#   Newton identity e2 sign flipped (H)                          cotangent_sum                FIRED
#
# Clean runs after the plants: fast 11 passed, 2 skipped; --slow 13 passed.
