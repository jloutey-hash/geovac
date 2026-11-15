"""Guard for records AE, AE2 and AE2b (docs/predictions/2026-11-15_aa_bc_equal_arms*.md): the three-centre integrals

    S = (1s_A 1s_A | 1s_B 1s_C),  PC = (1s_A 1s_A | 1s_B 2pz_C),  PD = (1s_A 1s_A | 1s_B 3dsigma_C)

at equal arms, D1 = |BC| = |BA| = D, as Legendre series in x = cos(angle C B A).  Paper 59 sec:scope ("At equal
arms").  Written 2026-11-15 as a separate step after the paper edit (CLAUDE.md sec. 9).  The engine is the ported one
of tests/test_paper59_aa_bc_p.py, with the record-AE change (the middle radial region D1 < r < D2 is empty).

Claims guarded, each with the wrong answer it rejects:

(1) sign rule: the Legendre coefficients of (1 - x)^alpha have sign sign(Gamma(alpha+1)/Gamma(alpha+1-L)) (-1)^L
    (constant for L > alpha + 1), checked by quadrature; rejects AE's registered (-1)^L.
(2) the PC bracket 3 K33/(2D) + zB (K33 - (99/15) K53)/3 is proportional to (zB D - 9), and K53/K33 = 5/22; rejects a
    dropped term of the three (each moves the zero away from 9).
(3) per-channel values: the exact split at equal arms equals a direct-projection quadrature split at r = D; the summed
    series at E1, x = 1/2, equals the record's referee values.
(4) decay: beta_eff at L = 60 sits near the (1 - x)^{7/2} toy for S and the (1 - x)^{9/2} toy for PC, PD (E1), and
    near the (1 - x)^{11/2} toy for PC at E3 (zB D = 9), where PD stays at 9/2; rejects PC/PD keeping a 7/2 term
    (beta near -8) and an uncancelled PC at E3 (beta near -10).
(5) amplitudes: F_L / predicted at L = 59, 60, 61 positive and near 1 with the closed-form amplitudes and the corrected
    sign; PC at E3 against the registered (1 - x)^{11/2} amplitude; rejects the (-1)^L sign and a dropped B-gradient
    term in PC (about 40%).

Fire-test plants and outcomes: the comment block at the end of this file.
"""
from __future__ import annotations

import importlib.util
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import pytest

_spec = importlib.util.spec_from_file_location('_ap_guard_eq', Path(__file__).with_name('test_paper59_aa_bc_p.py'))
AP = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(AP)
AB = AP.AB
Field = AP.Field

CASES = {'E1': dict(D='2', zB='13/10', zC='1', zA='9/10'),
         'E3': dict(D='3', zB='3', zC='1', zA='9/10')}
DPS = 400          # the float conversion loses about 4.5 L digits (record AE2)


def params(case):
    g = {k: Fraction(v) for k, v in CASES[case].items()}
    return g['D'], g['zB'], g['zC'], g['zA']


# ---------------------------------------------------------------------------------------------------------------
# Engine at equal arms (record AE: two regions, [0, D] inner/near and [D, inf) outer/far)
# ---------------------------------------------------------------------------------------------------------------


def integrate_equal(dens, kers, D, zB, I, K):
    regions = [(Fraction(0), D, True, True), (D, None, False, False)]
    res = {}
    for (A, B, inner, near) in regions:
        for (e0d, bd), pd in dens[inner].items():
            for (e0k, bk), pk in kers[near].items():
                P = AB._conv(pd, pk, K)
                if not P:
                    continue
                bb = bd + bk - zB
                tab = I.table(bb, A, B, min(min(P) + 2, -1), max(max(P) + 2, -1))
                for n, v in P.items():
                    for (e1, kind, arg), u in tab[n + 2].items():
                        AB.add(res, (e0d + e0k + e1, kind, arg), K(-2) * v * u)
    return res


def coefficients(case, Lmax, integral):
    K = Field(None)
    D, zB, zC, zA = params(case)
    it, kt = AB.bessel_tables(K, Lmax + 3)
    I = AB.Integrals(K)
    dfun = {'S': AB.density, 'PC': AP.dens_pz, 'PD': AP.dens_dz2}[integral]
    co = {'NA': K(1), 'Y': K(-1), 'dY': K(zA)}
    out = []
    for L in range(Lmax + 1):
        dens = {inner: dfun(L, inner, D, zC, it, kt, K) for inner in (True, False)}
        kers = {nr: AP.kernel_factored(L, nr, co, it, kt, K, 2 * zA, D, 0) for nr in (True, False)}
        out.append({k: v for k, v in integrate_equal(dens, kers, D, zB, I, K).items() if v != 0})
    return out


def _mf(q):
    return mp.mpf(q.numerator) / q.denominator


def numeric(Q, case, integral):
    D, zB, zC, zA = (_mf(v) for v in params(case))
    tw = {'S': 2 * mp.sqrt(zB ** 3 * zC ** 3), 'PC': 2 * mp.sqrt(zB ** 3 * zC ** 5),
          'PD': 2 * mp.pi * mp.sqrt(zB ** 3 / mp.pi) * mp.sqrt(5 * (2 * zC) ** 7 / (4 * mp.pi * 720))}[integral]
    keys = set(k for q in Q for k in q)
    cv = {k: AB.const_value(k) for k in keys}
    return [tw * sum((mp.mpf(int(v.p)) / int(v.q) * cv[k] for k, v in q.items()), mp.mpf(0)) for q in Q]


_CACHE = {}


def F(case, integral, Lmax=61):
    key = (case, integral)
    if key not in _CACHE:
        Q = coefficients(case, Lmax, integral)
        with mp.workdps(DPS):
            _CACHE[key] = numeric(Q, case, integral)
    return _CACHE[key]


# ---------------------------------------------------------------------------------------------------------------
# Closed forms (record AE controls, corrected sign rule)
# ---------------------------------------------------------------------------------------------------------------


def c_L(a, L, sign_rule='corrected'):
    """(2L+1)/2 int_{-1}^{1} (1-x)^a P_L(x) dx in closed form."""
    mag = 2 ** (a + 1) * mp.gamma(a + 1) ** 2 / (mp.gamma(a + L + 2) * abs(mp.gamma(a + 1 - L)))
    sgn = mp.sign(mp.gamma(a + 1) / mp.gamma(a + 1 - L)) * (-1) ** L if sign_rule == 'corrected' else (-1) ** L
    return sgn * (2 * L + 1) / mp.mpf(2) * mag


def _c(q):
    return 2 ** (q + 3) * mp.pi ** 1.5 * mp.gamma((q + 3) / mp.mpf(2)) / mp.gamma(-q / mp.mpf(2))


def Kmn(m, n):
    return _c(m) * _c(n) / _c(m + n + 3)


def pc_bracket(D, zB):
    return 3 * Kmn(3, 3) / (2 * D) + zB * (Kmn(3, 3) - mp.mpf(99) / 15 * Kmn(5, 3)) / 3


def predicted(case, integral, L):
    D, zB, zC, zA = (_mf(v) for v in params(case))
    v3 = mp.mpf(2) / 3 * zA ** 4
    NB, f0 = mp.sqrt(zB ** 3 / mp.pi), mp.exp(-zB * D)
    if integral == 'S':
        A, a = NB * mp.sqrt(zC ** 3 / mp.pi) * f0 * (-zC) * v3 * Kmn(1, 3), mp.mpf(7) / 2
    elif integral == 'PC':
        A, a = NB * mp.sqrt(zC ** 5 / mp.pi) * f0 * zC * v3 * pc_bracket(D, zB), mp.mpf(9) / 2
    else:
        A, a = (NB * mp.sqrt(5 * (2 * zC) ** 7 / (4 * mp.pi * 720)) * f0 * zC * v3 * mp.mpf(99) / 30 * Kmn(5, 3),
                mp.mpf(9) / 2)
    return A * (2 * D * D) ** a * c_L(a, L)


A_E3_PC_11_2 = mp.mpf('-0.0048413')   # record AE2, control X: coefficient of (1 - x)^{11/2}, PC at E3


def beta(Fs, L):
    return float(mp.log(abs(Fs[L + 1] / Fs[L - 1])) / mp.log(mp.mpf(L + 1) / (L - 1)))


def beta_toy(a, L):
    return float(mp.log(abs(c_L(a, L + 1) / c_L(a, L - 1))) / mp.log(mp.mpf(L + 1) / (L - 1)))


# ---------------------------------------------------------------------------------------------------------------
# (1), (2) closed forms
# ---------------------------------------------------------------------------------------------------------------


@pytest.mark.parametrize('a', [mp.mpf(7) / 2, mp.mpf(9) / 2, mp.mpf(11) / 2])
def test_legendre_coefficient_sign_rule_matches_quadrature(a):
    """Corrected sign and magnitude against quadrature at L = 1..14; (-1)^L disagrees for L > alpha + 1."""
    with mp.workdps(30):
        wrong = 0
        for L in range(1, 15):
            q = (2 * L + 1) / mp.mpf(2) * mp.quad(lambda x: (1 - x) ** a * mp.legendre(L, x), [-1, 0, 0.9, 1])
            assert abs(c_L(a, L) / q - 1) < mp.mpf('1e-15'), (float(a), L)
            wrong += mp.sign(c_L(a, L, 'alternating')) != mp.sign(q)
        assert wrong >= 4


def test_pc_bracket_proportional_to_zB_D_minus_9():
    with mp.workdps(50):
        assert abs(Kmn(5, 3) / Kmn(3, 3) - mp.mpf(5) / 22) < mp.mpf('1e-45')
        for D in (mp.mpf(2), mp.mpf(3), mp.mpf(5) / 2):
            assert abs(pc_bracket(D, 9 / D)) < mp.mpf('1e-45') * abs(Kmn(3, 3))
            r1 = pc_bracket(D, mp.mpf(1) / 2) / (D / 2 - 9)
            r2 = pc_bracket(D, mp.mpf(7) / 3) / (D * 7 / 3 - 9)
            assert abs(r1 / r2 - 1) < mp.mpf('1e-45')


# ---------------------------------------------------------------------------------------------------------------
# (3) values
# ---------------------------------------------------------------------------------------------------------------


def quadrature_FL(case, integral, L, dps=30):
    with mp.workdps(dps):
        D, zB, zC, zA = (_mf(v) for v in params(case))
        nodes = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(5, mp.mp.prec)
        NB = mp.sqrt(zB ** 3 / mp.pi)
        NC = {'S': mp.sqrt(zC ** 3 / mp.pi), 'PC': mp.sqrt(zC ** 5 / mp.pi),
              'PD': mp.sqrt(5 * (2 * zC) ** 7 / (4 * mp.pi * 720))}[integral]
        ang = {'S': lambda s, t, r: 1, 'PC': lambda s, t, r: r * t - D,
               'PD': lambda s, t, r: (3 * (r * t - D) ** 2 - s * s) / 2}[integral]

        def integrand(r):
            rhoL = NB * NC * mp.exp(-zB * r) * AP._multipole(lambda s, t: ang(s, t, r) * mp.exp(-zC * s), r, D, L, nodes)
            vL = AP._multipole(lambda s, t: AB.potential_1s(s, zA), r, D, L, nodes)
            return 4 * mp.pi / (2 * L + 1) * r * r * rhoL * vL
        return mp.quad(integrand, [0, D, D + 4, mp.inf])


@pytest.mark.parametrize('integral', ['S', 'PC', 'PD'])
@pytest.mark.parametrize('L', [0, 3, 8])
def test_equal_arm_channel_matches_direct_projection(integral, L):
    q = quadrature_FL('E1', integral, L)
    with mp.workdps(30):
        assert abs(F('E1', integral)[L] / q - 1) < mp.mpf('1e-15'), (integral, L)


REFEREE_E1_HALF = {'S': '0.22145582331291747', 'PC': '-0.25281152601442053', 'PD': '0.14800720825235306'}  # R2 n=34


@pytest.mark.parametrize('integral', ['S', 'PC', 'PD'])
def test_summed_series_matches_referee_E1(integral):
    with mp.workdps(DPS):
        x = mp.mpf(1) / 2
        val = sum(F('E1', integral)[L] * mp.legendre(L, x) for L in range(62))
        assert abs(val - mp.mpf(REFEREE_E1_HALF[integral])) < 1e-10, (integral, float(val))


# ---------------------------------------------------------------------------------------------------------------
# (4) decay, (5) amplitudes
# ---------------------------------------------------------------------------------------------------------------


@pytest.mark.parametrize('case,integral,alpha', [('E1', 'S', 3.5), ('E1', 'PC', 4.5), ('E1', 'PD', 4.5),
                                                 ('E3', 'S', 3.5), ('E3', 'PC', 5.5), ('E3', 'PD', 4.5)])
def test_decay_exponent_at_L60(case, integral, alpha):
    with mp.workdps(DPS):
        b, bt = beta(F(case, integral), 60), beta_toy(mp.mpf(alpha), 60)
    assert abs(b - bt) < 0.3, (case, integral, b, bt)


@pytest.mark.parametrize('case,integral', [('E1', 'S'), ('E1', 'PC'), ('E1', 'PD'), ('E3', 'S'), ('E3', 'PD')])
def test_closed_form_amplitude_with_corrected_sign(case, integral):
    with mp.workdps(DPS):
        for L in (59, 60, 61):
            r = float(F(case, integral)[L] / predicted(case, integral, L))
            assert 0.95 < r < 1.05, (case, integral, L, r)


def test_E3_PC_next_order_amplitude():
    with mp.workdps(DPS):
        for L in (59, 60, 61):
            r = float(F('E3', 'PC')[L] / (A_E3_PC_11_2 * c_L(mp.mpf(11) / 2, L)))
            assert 0.9 < r < 1.1, (L, r)


# FIRE-TEST PLANTS (python debug/qa/fire_test.py tests/test_paper59_aa_bc_equal_arms.py --plant "OLD=>NEW" -k <sel>)
# 2026-11-15, written and fired as a separate step after the Paper 59 edit; 12/12 fired.
#     AE's (-1)^L sign rule     "sgn = mp.sign(mp.gamma(a + 1) / mp.gamma(a + 1 - L)) * (-1) ** L if=>sgn = (-1) ** L if"   -k amplitude; -k sign_rule
#     cusp term of PC dropped   "return 3 * Kmn(3, 3) / (2 * D) + zB=>return zB"                       -k bracket; -k amplitude
#     B-gradient h2 part lost   "zB * (Kmn(3, 3) - mp.mpf(99) / 15 * Kmn(5, 3)) / 3=>zB * Kmn(3, 3) / 3"   -k amplitude
#     wrong kernel outside D    "(D, None, False, False)]=>(D, None, False, True)]"                     -k projection
#     11/2 amplitude off 20%    "A_E3_PC_11_2 = mp.mpf('-0.0048413')=>A_E3_PC_11_2 = mp.mpf('-0.0058413')"   -k next_order
#     PC keeps 7/2 (E1)         "('E1', 'PC', 4.5)=>('E1', 'PC', 3.5)"                                 -k decay
#     no cancellation (E3)      "('E3', 'PC', 5.5)=>('E3', 'PC', 4.5)"                                 -k decay
#     referee off 2e-10         "'S': '0.22145582331291747'=>'S': '0.22145582351291747'"               -k summed
#     too few digits            "DPS = 400 =>DPS = 150 "                                               -k "decay or amplitude or next_order"
#     2pz node term dropped     --plant-in tests/test_paper59_aa_bc_p.py "(-K(Fraction(D1)), 0,=>(-K(Fraction(0)), 0,"   -k "projection or summed"
