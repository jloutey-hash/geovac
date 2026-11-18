"""Guard for record L52 (docs/predictions/2026-11-17_aa_bc_contact_amplitudes.md): the amplitude of the contact
branch of

    f(x) = (1s_A 1s_A | 1s_B 1s_C) = sum_L F_L P_L(x),      x = cos(angle C B A),   |BC| = D1, |BA| = D2,

at unequal arms, and of the difference D = f_kappa - f when the Coulomb line is replaced by a Yukawa line
e^{-kappa r12}/r12.  Paper 59 sec:scope, the paragraph "The same expansion gives the amplitude at unequal arms".
Written 2026-11-17 as a separate step after the paper edit (CLAUDE.md sec. 9).

Nothing here imports debug/l52_local.py or reads debug/data/l52_predictions.json.  The routes are this file's own:

  * potential of the 1s_A^2 density (Coulomb and Yukawa): radial Green-function integral, done in sympy from an
    antiderivative that is checked by differentiation (the record used Fourier partial fractions);
  * pinch constant C_pq: Fourier transforms of |s|^p (c_p c_q / c_{p+q+3} divided by the Beta function), compared
    with the record's Gamma-function form;
  * segment geometry and Jacobian: sympy, with A on the complexified circle;
  * Legendre coefficients of (x* - x)^gamma: Gegenbauer generating function plus the Gegenbauer -> Legendre
    connection formula, a convergent series (the record used a Rodrigues integral), checked against quadrature;
  * the fit: own least squares.

Claims guarded, each with the wrong answer it rejects:

(1) 5/2 statement: the s^-1 and s^1 coefficients of the potential vanish identically for every kappa; the s^3
    coefficient is 2 zA^4/3 for every kappa; the s^5 coefficient moves by kappa^2 zA^4/45.  Rejects a point charge,
    an s^1 term, a kappa-dependent s^3 term.
(2) closed forms: Phi7 at A1, A2, A3 and DPhi9 at A1, A3 equal the values typed from the record to 1e-18.  Rejects a
    swapped segment weight, the weight taken at C only, a wrong C_pq, a wrong Jacobian.
(3) measurement: from the committed exact coefficients (sha256 pinned), F_L / (Phi7 T_L) and D_L / (DPhi9 T_L)
    extrapolate to 1 within the registered tolerances (four-term fit, L = 120..240 step 20), within the paper's
    4e-6 and 6e-6, and within 1e-6 for a five-term fit.  Rejects any amplitude off by more than 1e-6.
(4) the committed coefficients are not taken on trust at low L: F_L and the Yukawa F_L by direct two-centre
    projection quadrature.  (L = 120..240 cannot be reached by quadrature; the 2000/2600-digit agreement is checked.)
(5) the paper paragraph carries the measured figures, the tolerances, the geometry counts and the scope sentence.

Fire-test plants and outcomes: the comment block at the end of this file.
"""
from __future__ import annotations

import hashlib
import json
import re
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import pytest
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'debug' / 'data'
PAPER = ROOT / 'papers' / 'group2_quantum_chemistry' / 'paper_59_elliptic_bessel_moment.tex'

# record L52, table "Cases": D1, D2, zB, zC, zA, kappa
CASES = {'A1': ('9/5', '13/5', '13/10', '1', '9/10', '7/10'),
         'A2': ('3/2', '12/5', '1', '6/5', '4/5', None),
         'A3': ('2', '3', '11/10', '9/10', '6/5', '1')}
XSTAR = {'A1': '125/117', 'A2': '89/80', 'A3': '13/12'}

# typed from the record (not read from any debug/ file)
PHI7 = {'A1': '0.48406711402945841531', 'A2': '0.31156574196371525736', 'A3': '2.4426119167955476068'}
DPHI9 = {'A1': '0.02636189906870907033', 'A3': '0.34932701023832986156'}

# sha256 of the committed coefficient files, CRLF folded to LF (git may check text out with either)
SHA = {'l52_run_A1_2000.json': '4de380e85c19973eae9c09e3b355000218801a9b987dd6c26e8d7704543c4def',
       'l52_run_A2_2000.json': 'a39877e97bd44fb8f5b8680428cd22fd5dc1e20b63c09a5da48cb8819d492480',
       'l52_run_A3_2000.json': '90e30ccb7708e0064563789d1d71b287eb31def5277491449d73d56584ae0226'}

TOL_P1, TOL_P2A = mp.mpf('1e-3'), mp.mpf('3e-3')      # registered
PAPER_P1, PAPER_P2A = mp.mpf('4e-6'), mp.mpf('6e-6')  # Paper 59's printed figures
TOL_FIVE = mp.mpf('1e-6')                             # the record reports at most 1e-7
LS = list(range(120, 241, 20))                        # registered window
DPS = 40


def _mf(text):
    q = Fraction(text)
    return mp.mpf(q.numerator) / q.denominator


def geometry(case):
    D1, D2, zB, zC, zA, kap = CASES[case]
    return _mf(D1), _mf(D2), _mf(zB), _mf(zC), _mf(zA), (None if kap is None else _mf(kap))


# ---------------------------------------------------------------------------------------------------------------
# (1) potential of the 1s_A^2 density, Coulomb and Yukawa, by the radial Green function
# ---------------------------------------------------------------------------------------------------------------

_s, _r, _z, _k = sp.symbols('s r z kappa', positive=True)


def _moment(a, lo, hi):
    """int_lo^hi r e^{-a r} dr, from an antiderivative checked by differentiation."""
    G = lambda x: -sp.exp(-a * x) * (a * x + 1) / a ** 2
    assert sp.simplify(sp.diff(G(_r), _r) - _r * sp.exp(-a * _r)) == 0
    return (0 if hi is sp.oo else G(hi)) - G(lo)


def density(r):
    return _z ** 3 / sp.pi * sp.exp(-2 * _z * r)


def yukawa_potential():
    """V_kappa(s) = (2 pi/(kappa s)) int_0^inf r rho(r) [e^{-kappa |s - r|} - e^{-kappa (s + r)}] dr, kappa != 2 z."""
    c0 = _z ** 3 / sp.pi
    inner = c0 * sp.exp(-_k * _s) * (_moment(2 * _z - _k, 0, _s) - _moment(2 * _z + _k, 0, _s))
    outer = c0 * (sp.exp(_k * _s) - sp.exp(-_k * _s)) * _moment(2 * _z + _k, _s, sp.oo)
    V = sp.simplify(2 * sp.pi / (_k * _s) * (inner + outer))
    return V


def coulomb_potential():
    """(1/s) int_0^s 4 pi r^2 rho dr + int_s^inf 4 pi r rho dr, in closed form."""
    return 1 / _s - sp.exp(-2 * _z * _s) * (1 / _s + _z)


def odd_coefficients(V, upto=5):
    ser = sp.series(V, _s, 0, upto + 2).removeO()
    return {n: sp.simplify(ser.coeff(_s, n)) for n in range(-1, upto + 1, 2)}


_SYM = {}


def sym():
    if not _SYM:
        Vk, Vc = yukawa_potential(), coulomb_potential()
        _SYM.update(Vk=Vk, Vc=Vc, ok=odd_coefficients(Vk), oc=odd_coefficients(Vc))
    return _SYM


def test_coulomb_potential_solves_poisson_and_is_the_kappa_to_zero_limit():
    S = sym()
    lap = lambda V: sp.diff(_s * V, _s, 2) / _s
    assert sp.simplify(-lap(S['Vc']) - 4 * sp.pi * density(_s)) == 0
    assert sp.simplify(-lap(S['Vk']) + _k ** 2 * S['Vk'] - 4 * sp.pi * density(_s)) == 0
    assert sp.simplify(sp.limit(S['Vk'], _k, 0) - S['Vc']) == 0
    assert sp.limit(S['Vk'] * _s, _s, sp.oo) == 0 and sp.limit(S['Vc'] * _s, _s, sp.oo) == 1


def test_no_point_charge_and_no_linear_term_so_no_three_halves_or_five_halves_branch():
    """The 3/2 and 5/2 amplitudes carry v_-1 and v_1 (pairs (-1,1) and (1,1)): both vanish identically."""
    S = sym()
    for co in (S['ok'], S['oc']):
        assert co[-1] == 0 and co[1] == 0
    # the direct statement, not through the series routine: V is finite at s = 0 and V'(0) = 0
    for V in (S['Vk'], S['Vc']):
        assert sp.limit(V * _s, _s, 0) == 0
        assert sp.limit(sp.diff(V, _s), _s, 0) == 0


def test_cubic_term_is_kappa_free_and_quintic_term_moves_by_kappa_squared():
    S = sym()
    assert sp.simplify(S['oc'][3] - sp.Rational(2, 3) * _z ** 4) == 0
    assert sp.simplify(S['ok'][3] - sp.Rational(2, 3) * _z ** 4) == 0
    assert not S['ok'][3].has(_k)
    assert sp.simplify(S['ok'][5] - S['oc'][5] - _k ** 2 * _z ** 4 / 45) == 0
    assert sp.simplify(S['oc'][5] - sp.Rational(8, 45) * _z ** 6) == 0


# ---------------------------------------------------------------------------------------------------------------
# (2) the closed-form amplitudes
# ---------------------------------------------------------------------------------------------------------------


def _ft(p):
    """|s|^p has Fourier transform _ft(p) |k|^{-p-3} in three dimensions."""
    p = mp.mpf(p)
    return 2 ** (p + 3) * mp.pi ** 1.5 * mp.gamma((p + 3) / 2) / mp.gamma(-p / 2)


def C_pq(p, q):
    """Pinch constant: (|s|^p * |s|^q)(a) = [_ft(p) _ft(q) / _ft(p+q+3)] |a|^{p+q+3}, and the segment weight
    u^{(p+1)/2} (1-u)^{(q+1)/2} integrates to Beta((p+3)/2, (q+3)/2)."""
    return _ft(p) * _ft(q) / _ft(p + q + 3) / mp.beta(mp.mpf(p + 3) / 2, mp.mpf(q + 3) / 2)


def cusp_coefficient(q, zB, zC):
    """Coefficient of e^{-zB r_B} |r - C|^q in the normalized 1s_B 1s_C product: N (-zC)^q / q!."""
    t = sp.Symbol('t')
    c = sp.series(sp.exp(-t), t, 0, q + 1).removeO().coeff(t, q)          # (-1)^q / q!
    return mp.sqrt(zB ** 3 * zC ** 3) / mp.pi * mp.mpf(int(c.p)) / int(c.q) * zC ** q


def dist_to_B_squared(u, D1, D2):
    """|A + u (C - A) - B|^2 at the contact x = x*."""
    return (1 - u) * D2 ** 2 + u * D1 ** 2


def segment_average(p, q, D1, D2, zB):
    f = lambda u: (u ** (mp.mpf(p + 1) / 2) * (1 - u) ** (mp.mpf(q + 1) / 2)
                   * mp.exp(-zB * mp.sqrt(dist_to_B_squared(u, D1, D2))))
    return mp.quad(f, [0, mp.mpf(1) / 4, mp.mpf(1) / 2, mp.mpf(3) / 4, 1])


def amplitude(p, q, vp, D1, D2, zB, zC):
    """Coefficient of (x* - x)^{(p+q+3)/2}: |a|^2 = 2 D1 D2 (x* - x)."""
    gam = mp.mpf(p + q + 3) / 2
    return vp * cusp_coefficient(q, zB, zC) * C_pq(p, q) * (2 * D1 * D2) ** gam * segment_average(p, q, D1, D2, zB)


def _rational_values(case):
    zA, kap = CASES[case][4], CASES[case][5]
    return {_z: sp.Rational(zA), _k: sp.Rational(kap if kap is not None else 0)}


def _coefficient(expr, case):
    val = sp.Rational(expr.subs(_rational_values(case)))
    return mp.mpf(int(val.p)) / int(val.q)


def phi7(case):
    D1, D2, zB, zC, zA, _ = geometry(case)
    return amplitude(3, 1, _coefficient(sym()['oc'][3], case), D1, D2, zB, zC)


def dphi9(case):
    """Leading amplitude of f_kappa - f: the cubic terms must cancel, the quintic difference leads."""
    D1, D2, zB, zC, zA, kap = geometry(case)
    S = sym()
    assert sp.simplify(S['ok'][3] - S['oc'][3]) == 0, 'a kappa-dependent s^3 term would make D lead with 7/2'
    return amplitude(5, 1, _coefficient(S['ok'][5] - S['oc'][5], case), D1, D2, zB, zC)


def test_segment_geometry_and_jacobian():
    """With A = D2 (sin g, 0, cos g), C = (0, 0, D1), B = 0: |A - C|^2 = 2 D1 D2 (x* - x), and on the segment at
    x = x* the squared distance to B is (1-u) D2^2 + u D1^2."""
    x, u, D1, D2 = sp.symbols('x u D1 D2', positive=True)
    A = sp.Matrix([D2 * sp.sqrt(1 - x ** 2), 0, D2 * x])
    C = sp.Matrix([0, 0, D1])
    xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
    a2 = sp.expand(((A - C).T * (A - C))[0])
    assert sp.simplify(a2 - 2 * D1 * D2 * (xs - x)) == 0
    P = A + u * (C - A)
    assert sp.simplify(sp.expand((P.T * P)[0]).subs(x, xs) - dist_to_B_squared(u, D1, D2)) == 0
    for case, text in XSTAR.items():
        d1, d2 = Fraction(CASES[case][0]), Fraction(CASES[case][1])
        assert (d1 ** 2 + d2 ** 2) / (2 * d1 * d2) == Fraction(text)


def test_pinch_constant_matches_the_gamma_form_and_the_record_values():
    with mp.workdps(DPS):
        for p, q in ((-1, 1), (3, 1), (5, 1), (3, 3), (-1, 3)):
            gamma_form = mp.pi ** 1.5 * mp.gamma(-mp.mpf(p + q + 3) / 2) / (mp.gamma(-mp.mpf(p) / 2)
                                                                           * mp.gamma(-mp.mpf(q) / 2))
            assert abs(C_pq(p, q) / gamma_form - 1) < mp.mpf('1e-35'), (p, q)
        assert abs(C_pq(3, 1) / (-2 * mp.pi / 35) - 1) < mp.mpf('1e-35')
        assert abs(C_pq(5, 1) / (-2 * mp.pi / 63) - 1) < mp.mpf('1e-35')


def test_pinch_constant_against_a_numerical_convolution():
    """An independent number for the constant: int |s|^p |s - a|^q d^3s converges for -3 < p, q and p + q < -3 and
    equals [C_pq Beta] |a|^{p+q+3}; taken by quadrature in prolate spheroidal coordinates with foci 0 and a,
    |a| = 1: |s| = (xi + eta)/2, |s - a| = (xi - eta)/2, d^3s = (pi/4) (xi^2 - eta^2) dxi deta = pi |s| |s - a|."""
    with mp.workdps(20):
        for p, q, tol in ((-2, -2, '1e-12'), (mp.mpf(-9) / 4, mp.mpf(-7) / 4, '1e-8')):
            f = lambda xi, eta: mp.pi * ((xi + eta) / 2) ** (p + 1) * ((xi - eta) / 2) ** (q + 1)
            num = mp.quad(f, [1, 2, 8, mp.inf], [-1, 0, 1])
            assert abs(num / (C_pq(p, q) * mp.beta(mp.mpf(p + 3) / 2, mp.mpf(q + 3) / 2)) - 1) < mp.mpf(tol), (p, q)


@pytest.mark.parametrize('case', ['A1', 'A2', 'A3'])
def test_seven_halves_amplitude_equals_the_registered_value(case):
    with mp.workdps(DPS):
        assert abs(phi7(case) / mp.mpf(PHI7[case]) - 1) < mp.mpf('1e-18'), mp.nstr(phi7(case), 25)


@pytest.mark.parametrize('case', ['A1', 'A3'])
def test_nine_halves_amplitude_of_the_yukawa_difference_equals_the_registered_value(case):
    with mp.workdps(DPS):
        assert abs(dphi9(case) / mp.mpf(DPHI9[case]) - 1) < mp.mpf('1e-18'), mp.nstr(dphi9(case), 25)


def test_segment_formula_reduces_to_the_fourier_constant_at_equal_arms():
    """At D1 = D2 = D the weight is constant on the segment, so the amplitude is the plain convolution constant
    (the equal-arm closed form of tests/test_paper59_aa_bc_equal_arms.py, K(1, 3))."""
    with mp.workdps(DPS):
        D, zB, zC, zA = mp.mpf(2), mp.mpf(13) / 10, mp.mpf(1), mp.mpf(9) / 10
        v3 = 2 * zA ** 4 / 3
        direct = (v3 * mp.sqrt(zB ** 3 * zC ** 3) / mp.pi * (-zC) * mp.exp(-zB * D)
                  * _ft(3) * _ft(1) / _ft(7) * (2 * D * D) ** (mp.mpf(7) / 2))
        assert abs(amplitude(3, 1, v3, D, D, zB, zC) / direct - 1) < mp.mpf('1e-25')


# ---------------------------------------------------------------------------------------------------------------
# Legendre coefficients of (x* - x)^gamma
# ---------------------------------------------------------------------------------------------------------------


def legendre_coefficient(gam, xs, L):
    """(2L+1)/2 int_{-1}^{1} (xs - x)^gam P_L(x) dx, xs > 1.

    With t = xs - sqrt(xs^2 - 1): 1 - 2 x t + t^2 = 2 t (xs - x), so (xs - x)^gam = (2t)^{-gam} sum_n C_n^{(-gam)}(x) t^n,
    and C_n^{(lam)} = sum_k (2n - 4k + 1) (lam)_{n-k} (lam - 1/2)_k / ((3/2)_{n-k} k!) P_{n-2k}."""
    gam, lam = mp.mpf(gam), -mp.mpf(gam)
    t = xs - mp.sqrt(xs * xs - 1)
    term = mp.rf(lam, L) / mp.rf(mp.mpf(3) / 2, L)
    total, k = mp.mpf(0), 0
    while True:
        total += term
        term *= (lam + L + k) * (lam - mp.mpf(1) / 2 + k) / ((mp.mpf(3) / 2 + L + k) * (k + 1)) * t * t
        k += 1
        if k > 8 and abs(term) < mp.eps * abs(total):
            break
    return (2 * t) ** (-gam) * t ** L * (2 * L + 1) * total


@pytest.mark.parametrize('gam', [mp.mpf(7) / 2, mp.mpf(9) / 2])
def test_legendre_coefficients_of_the_branch_match_quadrature(gam):
    with mp.workdps(30):
        xs = _mf(XSTAR['A1'])
        for L in (0, 1, 2, 5, 9, 14):
            direct = (2 * L + 1) / mp.mpf(2) * mp.quad(lambda x: (xs - x) ** gam * mp.legendre(L, x),
                                                         mp.linspace(-1, 1, 9))
            assert abs(legendre_coefficient(gam, xs, L) / direct - 1) < mp.mpf('1e-20'), (float(gam), L)


def test_legendre_coefficients_resum_to_the_branch():
    """The series at large L, checked by summing it: sum_L T_L P_L(x) = (x* - x)^{7/2} to 1e-75 at A3's x*, where
    T_240 is about 1e-53, so an error of 1e-20 relative at L = 240 would show."""
    with mp.workdps(95):
        xs, x, gam = _mf(XSTAR['A3']), mp.mpf(4) / 5, mp.mpf(7) / 2
        P0, P1, tot = mp.mpf(1), x, mp.mpf(0)
        for L in range(431):
            tot += legendre_coefficient(gam, xs, L) * P0
            P0, P1 = P1, ((2 * L + 3) * x * P1 - (L + 1) * P0) / (L + 2)
        assert mp.mpf('1e-56') < abs(legendre_coefficient(gam, xs, 240)) < mp.mpf('1e-50')
        assert abs(tot - (xs - x) ** gam) < mp.mpf('1e-75')


# ---------------------------------------------------------------------------------------------------------------
# (3) the measurement, from the committed exact coefficients
# ---------------------------------------------------------------------------------------------------------------


def _bytes(name):
    return (DATA / name).read_bytes().replace(b'\r\n', b'\n')


def load(case):
    return json.loads(_bytes(f'l52_run_{case}_2000.json').decode())


@pytest.mark.parametrize('name', sorted(SHA))
def test_committed_coefficient_files_are_the_scored_ones(name):
    assert hashlib.sha256(_bytes(name)).hexdigest() == SHA[name]


def extrapolate(Ls, R, terms):
    """Least squares R_L = r0 + c1/L + ... + c_{terms-1}/L^{terms-1}; returns r0."""
    rows = [[mp.mpf(1) / mp.mpf(L) ** j for j in range(terms)] for L in Ls]
    n = terms
    M = mp.matrix(n, n)
    b = mp.matrix(n, 1)
    for row, val in zip(rows, R):
        for i in range(n):
            b[i] += row[i] * val
            for j in range(n):
                M[i, j] += row[i] * row[j]
    return mp.lu_solve(M, b)[0]


def ratios(case, key, Ls):
    D1, D2 = geometry(case)[:2]
    xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
    amp, gam = (phi7(case), mp.mpf(7) / 2) if key == 'F0' else (dphi9(case), mp.mpf(9) / 2)
    data = load(case)[key]
    return [mp.mpf(data[str(L)]) / (amp * legendre_coefficient(gam, xs, L)) for L in Ls]


CLAUSES = [('A1', 'F0'), ('A2', 'F0'), ('A3', 'F0'), ('A1', 'D'), ('A3', 'D')]


@pytest.mark.parametrize('case,key', CLAUSES)
def test_extrapolated_ratio_is_one_within_the_registered_tolerance(case, key):
    with mp.workdps(DPS):
        R = ratios(case, key, LS)
        r0 = extrapolate(LS, R, 4)
        assert abs(r0 - 1) <= (TOL_P1 if key == 'F0' else TOL_P2A), (case, key, mp.nstr(r0 - 1, 5))
        # the raw ratio is still 1-2% low at L = 240 and rising: the agreement is reached by extrapolation
        assert all(a < b for a, b in zip(R, R[1:])) and mp.mpf('0.97') < R[-1] < mp.mpf('0.995')


def test_four_term_residuals_are_the_papers_figures():
    """Paper 59 prints 'within 4e-6' (7/2, three geometries) and 'within 6e-6' (9/2, two geometries)."""
    with mp.workdps(DPS):
        dev = {ck: abs(extrapolate(LS, ratios(*ck, LS), 4) - 1) for ck in CLAUSES}
        p1 = max(v for (c, k), v in dev.items() if k == 'F0')
        p2 = max(v for (c, k), v in dev.items() if k == 'D')
        assert PAPER_P1 / 2 < p1 <= PAPER_P1, mp.nstr(p1, 5)
        assert PAPER_P2A / 2 < p2 <= PAPER_P2A, mp.nstr(p2, 5)


@pytest.mark.parametrize('case,key', CLAUSES)
def test_five_term_fit_closes_the_residual(case, key):
    """The four-term residual is extrapolation bias: one more correction term brings it below 1e-6, on the
    registered window and on every L from 120 to 240."""
    with mp.workdps(DPS):
        assert abs(extrapolate(LS, ratios(case, key, LS), 5) - 1) <= TOL_FIVE
        dense = list(range(120, 241))
        assert abs(extrapolate(dense, ratios(case, key, dense), 5) - 1) <= TOL_FIVE


@pytest.mark.parametrize('case', ['A1', 'A2', 'A3'])
def test_stored_difference_and_precision_twin(case):
    """D = Fk - F0 as stored, and the 2600-digit run agrees with the 2000-digit one at every L."""
    with mp.workdps(60):
        lo = load(case)
        hi = json.loads(_bytes(f'l52_run_{case}_2600.json').decode())
        for key in ('F0', 'Fk', 'D'):
            if key not in lo:
                assert case == 'A2'
                continue
            assert sorted(lo[key], key=int) == [str(L) for L in range(241)]
            for L in range(241):
                assert abs(mp.mpf(lo[key][str(L)]) / mp.mpf(hi[key][str(L)]) - 1) < mp.mpf('1e-35'), (key, L)
        if 'D' in lo:
            for L in range(241):
                d = mp.mpf(lo['Fk'][str(L)]) - mp.mpf(lo['F0'][str(L)])
                assert abs(d / mp.mpf(lo['D'][str(L)]) - 1) < mp.mpf('1e-30'), L


# ---------------------------------------------------------------------------------------------------------------
# (4) the committed coefficients against direct two-centre projection quadrature
# ---------------------------------------------------------------------------------------------------------------


def _legendre(L, x):
    P0, P1 = mp.mpf(1), x
    for n in range(L):
        P0, P1 = P1, ((2 * n + 3) * x * P1 - (n + 1) * P0) / (n + 2)
    return P0


def _projection(fun, r, Dc, L, nodes):
    """(2L+1)/2 int_{-1}^{1} fun(|r n - Dc e|) P_L(n.e) d(n.e), in the distance variable (smooth integrand)."""
    lo, hi = abs(r - Dc), r + Dc
    half, mid = (hi - lo) / 2, (hi + lo) / 2
    tot = mp.mpf(0)
    for xn, wn in nodes:
        sv = mid + half * xn
        tot += wn * fun(sv) * _legendre(L, (r * r + Dc * Dc - sv * sv) / (2 * r * Dc)) * sv
    return (2 * L + 1) / mp.mpf(2) * half * tot / (r * Dc)


def quadrature_FL(case, L, yukawa=False, dps=22, degree=5):
    """F_L = 4 pi/(2L+1) int r^2 rho_L(r) v_L(r) dr about B: rho_L of the 1s_B 1s_C product along B-C, v_L of the
    potential of the 1s_A^2 density along B-A."""
    with mp.workdps(dps):
        D1, D2, zB, zC, zA, kap = geometry(case)
        S = sym()
        if yukawa:
            assert kap is not None
            pot = sp.lambdify(_s, S['Vk'].subs(_rational_values(case)), 'mpmath')
        else:
            pot = sp.lambdify(_s, S['Vc'].subs(_rational_values(case)), 'mpmath')
        nodes = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(degree, mp.mp.prec)
        N = mp.sqrt(zB ** 3 * zC ** 3) / mp.pi

        def integrand(r):
            rhoL = N * mp.exp(-zB * r) * _projection(lambda t: mp.exp(-zC * t), r, D1, L, nodes)
            return 4 * mp.pi / (2 * L + 1) * r * r * rhoL * _projection(pot, r, D2, L, nodes)

        cuts = [mp.mpf(0), D1 / 2, D1, D2, D2 + 2, D2 + 6, D2 + 14, D2 + 30]
        tot = mp.mpf(0)
        for a, b in zip(cuts, cuts[1:]):
            half, mid = (b - a) / 2, (b + a) / 2
            tot += half * sum(wn * integrand(mid + half * xn) for xn, wn in nodes)
        return tot


@pytest.mark.parametrize('case,L', [('A1', 0), ('A1', 3), ('A2', 1), ('A3', 2)])
def test_coulomb_coefficient_matches_direct_quadrature(case, L):
    q = quadrature_FL(case, L)
    with mp.workdps(30):
        assert abs(mp.mpf(load(case)['F0'][str(L)]) / q - 1) < mp.mpf('1e-12'), (case, L)


@pytest.mark.parametrize('case,L', [('A1', 1), ('A3', 0)])
def test_yukawa_coefficient_matches_direct_quadrature(case, L):
    q = quadrature_FL(case, L, yukawa=True)
    with mp.workdps(30):
        assert abs(mp.mpf(load(case)['Fk'][str(L)]) / q - 1) < mp.mpf('1e-12'), (case, L)


@pytest.mark.slow
@pytest.mark.parametrize('case,key,L', [('A1', 'F0', 8), ('A1', 'Fk', 8), ('A3', 'F0', 12), ('A3', 'Fk', 12),
                                        ('A2', 'F0', 12)])
def test_higher_coefficients_match_direct_quadrature(case, key, L):
    q = quadrature_FL(case, L, yukawa=(key == 'Fk'), dps=40, degree=7)
    with mp.workdps(40):
        assert abs(mp.mpf(load(case)[key][str(L)]) / q - 1) < mp.mpf('1e-18'), (case, key, L)


# ---------------------------------------------------------------------------------------------------------------
# (5) the paper paragraph
# ---------------------------------------------------------------------------------------------------------------

START = r'The same expansion gives the amplitude at unequal arms'
END = r'2026-11-17\_aa\_bc\_contact\_amplitudes.md'


def paragraph():
    text = PAPER.read_text(encoding='utf-8')
    assert text.count(START) == 1 and text.count(END) == 1
    a, b = text.index(START), text.index(END)
    assert 0 < b - a < 2500
    return re.sub(r'\s+', ' ', text[a:b])


@pytest.mark.parametrize('needle', [
    r'for $1s$ functions',
    r'derived before computing',
    r'analyticity of the weight along the complexified segment assumed',
    r'At three geometries, one not computed before',
    r'exact Legendre coefficients up to $\ell=240$',
    r'within $4\times10^{-6}$, against a registered tolerance of $10^{-3}$',
    r'$|\mathbf{s}|^3$ term of the potential of the A density does not depend on $\kappa$',
    r'leads with $(x^{*}-x)^{9/2}$',
    r'within $6\times10^{-6}$ at two geometries (tolerance $3\times10^{-3}$)',
    r'This concerns the function, not the exponents of its minimal operator',
])
def test_paper_paragraph_states(needle):
    assert needle in paragraph(), needle


def test_paper_paragraph_tier_and_counts_agree_with_this_file():
    text = PAPER.read_text(encoding='utf-8')
    assert text[text.index(START) - 40:text.index(START)].rstrip().endswith(r'\textbf{[SYMBOLIC + MEASURED]}')
    assert sum(k == 'F0' for _, k in CLAUSES) == 3 and sum(k == 'D' for _, k in CLAUSES) == 2
    assert (PAPER_P1, PAPER_P2A, TOL_P1, TOL_P2A) == (mp.mpf('4e-6'), mp.mpf('6e-6'), mp.mpf('1e-3'), mp.mpf('3e-3'))


# FIRE-TEST PLANTS (python debug/qa/fire_test.py tests/test_paper59_contact_amplitudes.py --plant "OLD=>NEW" -k <sel>)
# 2026-11-17, written and fired as a separate step after the Paper 59 edit; 33/33 fired.  Where a plant is listed
# with two selectors it fired under each alone (closed-form test, then the measurement with the typed values unused).
#     segment weight swapped u <-> 1-u        in segment_average                       -k registered_value; -k registered_tolerance
#     weight taken at C only                  dist_to_B_squared(u, => dist_to_B_squared(1,          same two selectors
#     C_pq times 2; times sqrt(pi)            in C_pq's return                         -k pinch_constant; -k registered_tolerance
#     2^(p+3) -> 2^(p+2) in _ft               (wrong in both closed forms alike)       -k numerical_convolution
#     v3 kappa-dependent                      yukawa_potential: + _k ** 2 * _s ** 3 / 7     -k cubic_term; -k nine_halves
#     s^1 term added                          yukawa_potential: + _s / 7; coulomb_potential: + _z ** 2 * _s   -k no_point_charge
#     point charge left                       coulomb_potential: (1 / _s + _z) => (_z)  -k no_point_charge
#     Jacobian 2 D1 D2 -> 2 D1^2              in amplitude; in the symbolic assertion  -k registered_value; -k registered_tolerance; -k jacobian
#     typed Phi7 off 1e-3 (A1), 1e-15 (A3)    PHI7                                     -k seven_halves
#     typed DPhi9 off 1e-3 (A3)               DPHI9                                    -k nine_halves
#     D from the full v5, not the difference  in dphi9                                 -k nine_halves
#     cusp coefficient without zC^q           in cusp_coefficient                      -k seven_halves
#     connection formula lam - 1/2 -> lam     in legendre_coefficient                  -k legendre_coefficients; -k registered_tolerance
#     hash constant altered                   SHA                                      -k committed
#     one F_L (A1, L = 180) off 1e-6 in both precision files, D kept consistent (DATA pointed at a temp copy)
#                                                                                      -k committed; -k five_term; -k registered_tolerance; -k papers_figures
#     paper-figure constants loosened 10x     PAPER_P1, PAPER_P2A                      -k "papers_figures or tier_and_counts"
#     scope sentence deleted; 4e-6 -> 4e-7    PAPER pointed at a temp copy             -k paper_paragraph
#     potential off 1e-9 (zA factor)          coulomb_potential                        -k coulomb_coefficient
#     five-term bound read with three terms   in test_five_term_fit_closes_the_residual -k five_term
# Not reachable: an error in a stored coefficient at L >= 120 smaller than about 1e-7 relative (below the five-term
# bound; the quadrature referee stops at L = 12, slow tier).  The hash and the 2600-digit twin are the only cover.
