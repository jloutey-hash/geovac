"""Frozen falsifier: S^(3) W10 identification -- full closed form.

W10-identification sprint 2026-06-12 (debug/sprint_s3_w10_identification_memo.md;
driver debug/s3_w10_ident.py).  Closes the v4.7.0 five-generator block:

  (i)   The two Charlton-Hoffman symmetry-pair sums (Thm 2.21 of
        arXiv:2204.14183 at depth 3, phi=0, explicit product terms):
          t3(4,4,2)+t3(2,4,4) and t3(4,3,3)+t3(3,3,4)
        are identified elements of Q[pi^2, zeta(3), zeta(5), zeta(7),
        zeta(5,3)] -- gated against independently computed 200-dps pins.
  (ii)  Each of the four weight-10 triples is individually identified:
        catalogued ring + rational combination of the antisymmetric
        doubles a(8,2), a(7,3), a(6,4).
  (iii) COLLAPSE THEOREM: in the assembly-normalized W10 corollary the
        three a-generators cancel identically (symbolic check) -- the
        S^(3) closed form needs NO level-2 depth-2 generator.
  (iv)  S^(3) full closed form in Q[pi^2, ln2, zeta(odd), zeta(5,3)]
        gated against the canonical 200-dps value.  zeta(5,3) (descending,
        sum_{m>n} m^-5 n^-3) is computed LIVE -- an error in its
        evaluation or its role in the closed form fails the gate.

Conventions: t-values/zetas descending (GeoVac/Hoffman), first index on
the largest summation variable.  All pinned values 195 digits.

Runtime: fast group < 30 s; the 200-dps gates are marked slow.
"""
from __future__ import annotations

from fractions import Fraction

import mpmath
from mpmath import mp, mpf, zeta as mpzeta
import pytest
import sympy as sp

# ---------------------------------------------------------------------------
# Pinned values (195 digits; sources: s3_pslq_cache.json @220,
# s3_w10_value_cache.json @200, s3_closure_recompute.json @200)
# ---------------------------------------------------------------------------

PIN = {
    't3(2,4,4)': "0.00168687050150783661538197068274966474752646087000638409735719295030498317883990145643893177061430175563597164878840729497833788916686833349979920079758076276833337700565032189163617113869258778458",
    't3(3,3,4)': "0.000618981623446194587378816197576214133861120816940211034673339437845041281870689752346479047346245158982876215444463792295895841732281123327185313957861024682928717288051009546282946480661879191946",
    't3(4,4,2)': "0.0000302890397680854182277765286132246505511347749262965464257786214990266466593991153262680790167942835009067627037141327053077952451803720817978339728289532832232944500387696470478771445079045100955",
    't3(4,3,3)': "0.0000937868838266138978874548513424225989145710475128365280310425547331356030733345887597378101308961115072081195311133503981307893675684485653794124586455843238510610942937864603578076366157614017736",
    'a(8,2)': "-0.233563963441682925406641612738528039559382035989258319411116771762948041694975930638551195263829120393802185271193800757529836206560688547135143735002734394786479200779928789561333332113950644457",
    'a(7,3)': "-0.0513345410541781283873714166580174676490035408002631678076462337464680455114927527221713871824133154375601636450493335137361083444680913423618965879443508030240661284118160909210277505566156651559",
    'a(6,4)': "-0.0132332547312598284231850248221570285930904276478509083453955261029051383132348039851350580884181382020593983853830567092718450072489514643113336221309244626377077543094472880071257000623730072294",
    'S3': "31.5725612075120227547619295450689833719575844695771941287005369766388345706854480909447330237578823749984324629677512219546838264469270715404509652942715095183123080504580446671013547025451659223",
}

# ---------------------------------------------------------------------------
# Frozen closed forms (exact rationals; driver s3_w10_ident.py `assemble`,
# data debug/data/s3_w10_ident_result.json).  Symbols: z53 = zeta(5,3)
# descending; a82/a73/a64 = antisymmetric double t-values.
# ---------------------------------------------------------------------------

_PI, _Z3, _Z5, _Z7, _Z9, _Z53, _LN2 = sp.symbols('pi z3 z5 z7 z9 z53 ln2')
_A82, _A73, _A64 = sp.symbols('a82 a73 a64')
R = sp.Rational

P1_FORM = (R(5, 6144) * _PI**4 * _Z3**2 - R(5, 256) * _PI**2 * _Z3 * _Z5
           - R(1, 256) * _PI**2 * _Z53 + R(137, 99532800) * _PI**10)
P2_FORM = (R(37, 12288) * _PI**4 * _Z3**2 - R(25, 1024) * _PI**2 * _Z3 * _Z5
           - R(889, 2048) * _Z3 * _Z7 - R(1, 256) * _PI**2 * _Z53
           + R(251, 58060800) * _PI**10)

T3_FORMS = {
    't3(2,4,4)': (-_A64 / 4 + _A82 / 4 - R(1, 6144) * _PI**4 * _Z3**2
                  - R(5, 512) * _PI**2 * _Z3 * _Z5 - R(1, 512) * _PI**2 * _Z53
                  + R(373, 174182400) * _PI**10),
    't3(3,3,4)': (-_A64 / 4 + _A73 / 4 + R(1, 12288) * _PI**4 * _Z3**2
                  + R(5, 1024) * _PI**2 * _Z3 * _Z5
                  - R(889, 4096) * _Z3 * _Z7 - R(1, 512) * _PI**2 * _Z53
                  + R(251, 116121600) * _PI**10),
    't3(4,4,2)': (_A64 / 4 - _A82 / 4 + R(1, 1024) * _PI**4 * _Z3**2
                  - R(5, 512) * _PI**2 * _Z3 * _Z5 - R(1, 512) * _PI**2 * _Z53
                  - R(533, 696729600) * _PI**10),
    't3(4,3,3)': (_A64 / 4 - _A73 / 4 + R(3, 1024) * _PI**4 * _Z3**2
                  - R(15, 512) * _PI**2 * _Z3 * _Z5
                  - R(889, 4096) * _Z3 * _Z7 - R(1, 512) * _PI**2 * _Z53
                  + R(251, 116121600) * _PI**10),
}

# Identified part of S^(3) (assembly sprint v4.7.0, 26 terms):
IDENT_FORM = (
    128 * _PI**2 * _LN2**2 - R(32, 3) * _PI**4 * _LN2**2
    + R(7232, 81) * _PI**2 * _LN2 - R(35504, 243) * _PI**4 * _LN2
    + R(352, 15) * _PI**6 * _LN2 - R(314, 315) * _PI**8 * _LN2
    - 96 * _PI**2 * _LN2 * _Z3 + 16 * _PI**4 * _LN2 * _Z3
    + 20 * _PI**2 * _Z3**2 - R(2200, 27) * _PI**2 * _Z3
    + R(4988, 81) * _PI**4 * _Z3 - R(66, 5) * _PI**6 * _Z3
    + R(157, 210) * _PI**8 * _Z3
    - 80 * _PI**2 * _LN2 * _Z5 + R(7460, 81) * _PI**2 * _Z5
    + R(85, 3) * _PI**4 * _Z5 - R(11, 3) * _PI**6 * _Z5
    - 7 * _PI**2 * _Z7 + R(91, 12) * _PI**4 * _Z7
    - 30 * _PI**2 * _Z9
    + R(4096, 6561) * _PI**2 - R(109240, 19683) * _PI**4
    + R(27427, 1215) * _PI**6 - R(288859, 51030) * _PI**8
    + R(1381, 2835) * _PI**10 - R(17851, 1247400) * _PI**12)


def lam(k):
    return (1 - R(1, 2**k)) * sp.zeta(k)


def w10_assembly_over_256(t3_244, t3_334):
    """The v4.7.0 assembly-normalized corollary with a(4,3), a(4,2)
    substituted by their catalogued closed forms."""
    # catalogued lower-weight antisymmetric doubles (Hoffman App-A rows):
    lam_ = {k: (1 - R(1, 2**k)) * sp.zeta(k) for k in range(2, 11)}
    lam_[3] = R(7, 8) * _Z3
    lam_[5] = R(31, 32) * _Z5
    lam_[7] = R(127, 128) * _Z7
    lam_[9] = R(511, 512) * _Z9
    t43 = (-R(1, 2) * lam_[7] + R(6, 7) * lam_[3] * lam_[4]
           - R(10, 31) * lam_[2] * lam_[5])
    t34 = lam_[3] * lam_[4] - lam_[7] - t43
    a43 = sp.expand(t43 - t34)
    t42 = -R(1, 7) * lam_[6] + R(1, 7) * lam_[3]**2
    t24 = R(11, 28) * lam_[6] - R(1, 7) * lam_[3]**2
    a42 = sp.expand(t42 - t24)
    return sp.expand(
        (6 * lam_[10] - lam_[2] * lam_[8] - 2 * lam_[3] * lam_[7]
         - 3 * lam_[4] * lam_[6] - 4 * lam_[3] * a43 - 2 * lam_[4] * a42
         + _A82 + 2 * _A73 - 3 * _A64 - 4 * t3_244 - 8 * t3_334
         ).subs(sp.pi, _PI))


# ---------------------------------------------------------------------------
# Numeric helpers
# ---------------------------------------------------------------------------

def zeta53_desc():
    """zeta(5,3) descending = sum_{m>n} m^-5 n^-3, Levin-safe (log-free)."""
    z3v = mpzeta(3)
    def term(n):
        return n ** mpf(-5) * (z3v - mpzeta(3, n))
    return mpmath.nsum(term, [2, mp.inf], method='levin')


def eval_form(expr, z53v):
    subs = {_PI: mp.pi, _Z3: mpzeta(3), _Z5: mpzeta(5), _Z7: mpzeta(7),
            _Z9: mpzeta(9), _Z53: z53v, _LN2: mp.log(2),
            _A82: mpf(PIN['a(8,2)']), _A73: mpf(PIN['a(7,3)']),
            _A64: mpf(PIN['a(6,4)'])}
    e = sp.expand(expr)
    poly = sp.Poly(e, *subs.keys())
    val = mpf(0)
    bases = list(subs.values())
    for monom, coeff in poly.terms():
        cr = R(coeff)
        t = mpf(int(cr.p)) / mpf(int(cr.q))
        for base, ee in zip(bases, monom):
            if ee:
                t *= base ** ee
        val += t
    return val


# ---------------------------------------------------------------------------
# (iii) Collapse theorem -- pure symbolic, no numerics
# ---------------------------------------------------------------------------

def test_a_generator_collapse_symbolic():
    w10a = w10_assembly_over_256(T3_FORMS['t3(2,4,4)'],
                                 T3_FORMS['t3(3,3,4)'])
    for gen in (_A82, _A73, _A64):
        coeff = sp.expand(w10a).coeff(gen)
        assert coeff == 0, f"a-generator {gen} survives: coeff {coeff}"


def test_pair_difference_stuffles_exact():
    """The two clean difference relations follow from four descending
    stuffles t(x)*t(y,z) by exact combinatorial expansion."""
    def stuffle(x, y, z):
        out = {}
        for key in [('T3', x, y, z), ('T3', y, x, z), ('T3', y, z, x)]:
            out[key] = out.get(key, 0) + 1
        out[('T2', x + y, z)] = out.get(('T2', x + y, z), 0) + 1
        out[('T2', y, x + z)] = out.get(('T2', y, x + z), 0) + 1
        return out
    RA, RB = stuffle(2, 4, 4), stuffle(4, 2, 4)
    diff = {k: RA.get(k, 0) - RB.get(k, 0)
            for k in set(RA) | set(RB)}
    # t3(4,2,4) must cancel; t3(4,4,2) - t3(2,4,4) survives with +1/-1
    assert diff.get(('T3', 4, 2, 4), 0) == 0
    assert diff[('T3', 4, 4, 2)] == 1 and diff[('T3', 2, 4, 4)] == -1
    RC, RD = stuffle(3, 3, 4), stuffle(3, 4, 3)
    diff2 = {k: RC.get(k, 0) - RD.get(k, 0) for k in set(RC) | set(RD)}
    assert diff2.get(('T3', 3, 4, 3), 0) == 0
    assert diff2[('T3', 3, 3, 4)] == 2 and diff2[('T3', 4, 3, 3)] == -2


# ---------------------------------------------------------------------------
# (i), (ii), (iv) -- fast 50-dps gates (always on)
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def z53_fast():
    mp.dps = 60
    return zeta53_desc()


def test_symmetry_pair_sums_fast(z53_fast):
    mp.dps = 60
    s1 = mpf(PIN['t3(4,4,2)']) + mpf(PIN['t3(2,4,4)'])
    s2 = mpf(PIN['t3(4,3,3)']) + mpf(PIN['t3(3,3,4)'])
    assert abs(eval_form(P1_FORM, z53_fast) - s1) < mpf(10) ** -50
    assert abs(eval_form(P2_FORM, z53_fast) - s2) < mpf(10) ** -50


def test_triple_closed_forms_fast(z53_fast):
    mp.dps = 60
    for nm, form in T3_FORMS.items():
        resid = abs(eval_form(form, z53_fast) - mpf(PIN[nm]))
        assert resid < mpf(10) ** -50, f"{nm}: {resid}"


def test_s3_full_closed_form_fast(z53_fast):
    mp.dps = 60
    s3 = IDENT_FORM + 256 * w10_assembly_over_256(
        T3_FORMS['t3(2,4,4)'], T3_FORMS['t3(3,3,4)'])
    resid = abs(eval_form(s3, z53_fast) - mpf(PIN['S3']))
    assert resid < mpf(10) ** -50, f"S3 closed form residual {resid}"


def test_zeta53_only_new_constant(z53_fast):
    """The closed form's generator set: pi, ln2, odd zetas, zeta(5,3),
    and NOTHING else -- in particular no a-doubles and no zeta(11)."""
    s3 = sp.expand(IDENT_FORM + 256 * w10_assembly_over_256(
        T3_FORMS['t3(2,4,4)'], T3_FORMS['t3(3,3,4)']))
    syms = s3.free_symbols
    assert syms <= {_PI, _Z3, _Z5, _Z7, _Z9, _Z53, _LN2}, syms
    assert s3.coeff(_Z53) == 6 * _PI**2


# ---------------------------------------------------------------------------
# slow 200-dps gates
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_triple_closed_forms_200dps():
    mp.dps = 205
    z53v = zeta53_desc()
    for nm, form in T3_FORMS.items():
        resid = abs(eval_form(form, z53v) - mpf(PIN[nm]))
        assert resid < mpf(10) ** -180, f"{nm}: {resid}"


@pytest.mark.slow
def test_s3_full_closed_form_200dps():
    mp.dps = 205
    z53v = zeta53_desc()
    s3 = IDENT_FORM + 256 * w10_assembly_over_256(
        T3_FORMS['t3(2,4,4)'], T3_FORMS['t3(3,3,4)'])
    resid = abs(eval_form(s3, z53v) - mpf(PIN['S3']))
    assert resid < mpf(10) ** -180, f"S3 closed form residual {resid}"
