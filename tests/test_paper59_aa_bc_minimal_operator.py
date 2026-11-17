"""Inventory row 6: the minimal bond-angle operator of (1s_A 1s_A | 1s_B 1s_C).  Regression guard for record MO
(docs/predictions/2026-11-17_aa_bc_minimal_operator.md, scored GO 2026-11-17), written 2026-11-17 as a separate step
after the record (CLAUDE.md Sec. 9, guard-writing rule).  Paper 59 sec:scope, the record-AB paragraph.

Claims guarded (record AB's cases A1, x* = 125/117, and A2, x* = 89/80; record AK's primes; all mod p):
(1) the minimal x-operator (gcrd of annihilators) of the generic sequence sum_c w_c q_{L,c} (21 constants,
    L <= 1000, weight seeds 2026 and 7) has order 14, degree 111, at both cases, both primes, both weights, at
    zeta_A and at zeta_A + 1/997;
(2) AB's search operator (first null vector, order 17, degree 48) is right-divisible by it;
(3) the only finite point with a non-integer exponent is x*, exponents {0..10} and 3/2, 7/2, 9/2, the same under both
    weights;
(4) every other genuine (non-apparent, weight-fixed) point is +-1 or a linear point with integer exponents
    {-3, 0..12} or {-2, 0..12} and a forced logarithm;
(5) the exponent -3 points move under zeta_A -> zeta_A + 1/997; the exponent -2 point does not; per record, the
    number of moving points is 4, 2 (A1, the two primes) and 2, 4 (A2); A2's fixed point reads 2509/1980.
Not guarded: minimality over Q (the record's cap: mod-p readings, not certified); MO3/MO4 (piece carriage; an
instrument check, reported); degree > 1 integer-exponent factors (class I, excluded by the record).

Instrument: tests/paper59_intclass_support/ (verbatim copies, sha256 checked below; s1_controls.py added 2026-11-17
from debug/s1_controls.py as frozen by record MO, commit 56f10c90).  Nothing is imported from debug/.
Embedded data: tests/data/paper59_mo_2026_11_17.json -- for each (case, prime, side, seed) the generic sequence
(s1_controls.generic_seq over AB's 21 constant sequences, which the scored run cached as
debug/data/s1_coeffs_<case>_<p>[_shift].json) and its minimal operator (intclass_lib.minimal_operator, content-reduced),
built 2026-11-17 with the frozen pipeline.  The data are not trusted on their own:
  * sequences: their first 41 terms are recomputed from the engine copy in the fast tier, all 1001 in the slow tier;
  * operators: each is checked to ANNIHILATE its sequence by an independent route (the D-form operator is rewritten
    as sum a_i(x) T^i, T = (1 - x^2) d/dx, after multiplying by (1 - x^2)^R, and applied to the Legendre
    coefficients with the local x and T rules; this code is not the search's), and the gcrd is recomputed from the
    sequence (fast: one entry; slow: all sixteen).
Fire-test record: the comment block at the end of this file.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import flint
import numpy as np
import pytest

HERE = Path(__file__).parent
SUP = HERE / 'paper59_intclass_support'
DATA = HERE / 'data' / 'paper59_mo_2026_11_17.json'
# sha256 after CRLF -> LF (as tests/test_paper59_aa_bc_class.py).  s1_controls.py: raw sha256 = record MO's frozen
# 1c50016620c5f6291ac735028d9c87a658d8a9c226f28f08a1c667eba9436687.
FROZEN = {'s1_controls.py': 'a150322d746cb445b3b47e8264cb51a53cfe1dd81b21f8f1ab00eb8f21049453',
          'intclass_lib.py': '2a142e5b5667b432c7de4be885f5b6767daf3d209c86261076dc6ae5e4f641d4',
          'na_modguess.py': '84691c87a4a280b6823438af574c0fff6af9bf3f7ba7b76d93f2bbe866847d7d',
          'na_modanalyze.py': 'de92dc3acdcc7082001e32b50ef44db8c3bd75d990f615603b5f97f6ea178005',
          'na_exact.py': '85fc39f18bfe9231f47925436b72258dca7afddc6f9619257d9fafc708004afb',
          'aabc_exact.py': '218b7283f525d25ed8a87df61f0749ca6f576b20abcc4bde1272de32436b2d30'}

PRIMES = (1073741441, 1073741381)        # record AK's
CASES = ('A1', 'A2')
SIDES = ('base', 'shift')
SEEDS = (2026, 7)
SHIFT = Fraction(1, 997)
XSTAR = {'A1': Fraction(125, 117), 'A2': Fraction(89, 80)}
ORDER, DEGREE = 14, 111
AB_ORDER, AB_DEGREE = 17, 48
XSTAR_EXPS = sorted([Fraction(k) for k in range(11)] + [Fraction(3, 2), Fraction(7, 2), Fraction(9, 2)])
EXTRA_EXPS = {-3: sorted([-3] + list(range(13))), -2: sorted([-2] + list(range(13)))}
N_MOVING = {('A1', 1073741441): 4, ('A1', 1073741381): 2, ('A2', 1073741441): 2, ('A2', 1073741381): 4}
A2_FIXED = '2509/1980'
L_PREFIX = 40


def _load():
    sys.path.insert(0, str(SUP))
    try:
        mods = []
        for name in ('s1_controls', 'aabc_exact'):
            spec = importlib.util.spec_from_file_location(f'p59mo_{name}', SUP / f'{name}.py')
            m = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(m)
            mods.append(m)
    finally:
        sys.path.remove(str(SUP))
    return mods


C, X = _load()
I, G = C.I, C.G
D = json.loads(DATA.read_text())['entries']


def key(case, p, side, seed):
    return f'{case}|{p}|{side}|{seed}'


def op_of(k, p):
    return [flint.nmod_poly(c, p) for c in D[k]['op']]


# --------------------------------------------------------------------------------------- independent annihilation
def to_T_form(c, p):
    """(1 - x^2)^R L = sum_i a_i(x) T^i, T = (1 - x^2) d/dx, for L = sum_k c_k(x) d^k.  T^k = sum_m t_{k,m} d^m with
    t_{k+1,m} = (1 - x^2)(t_{k,m}' + t_{k,m-1}) and (1 - x^2)^m | t_{k,m}; peeled from the top.  None if a division
    is inexact (cannot happen for a valid operator)."""
    R = len(c) - 1
    w = flint.nmod_poly([1, 0, p - 1], p)
    t = [[flint.nmod_poly([1], p)]]
    for k in range(R):
        t.append([w * ((t[k][m].derivative() if m <= k else flint.nmod_poly([], p)) +
                       (t[k][m - 1] if m >= 1 else flint.nmod_poly([], p))) for m in range(k + 2)])
    rem = [w ** R * ck for ck in c]
    a = [None] * (R + 1)
    for i in range(R, -1, -1):
        q, r = divmod(rem[i], w ** i)
        if r.degree() >= 0:
            return None
        a[i] = q
        for m in range(i + 1):
            rem[m] = rem[m] - q * t[i][m]
    return a if all(r.degree() < 0 for r in rem) else None


def apply_to_legendre(c, seq, p):
    """Legendre coefficients of (1 - x^2)^R L f, f = sum seq_L P_L, on the range where they are determined."""
    a = to_T_form(c, p)
    assert a is not None
    n = len(seq) + 2
    L = np.arange(n, dtype=np.int64)
    inv = lambda v: np.array([pow(int(z) % p, p - 2, p) for z in v], dtype=np.int64)   # noqa: E731
    xm, xp = L % p * inv(2 * L - 1) % p, (L + 1) % p * inv(2 * L + 3) % p                # x P_l recurrences
    tm, tp = (L - 1) * L % p * inv(2 * L - 1) % p, (L + 1) * (L + 2) % p * inv(2 * L + 3) % p

    def xop(s):
        m = len(s) - 1
        out = xp[:m] * s[1:] % p
        out[1:] = (out[1:] + xm[1:m] * s[:m - 1]) % p
        return out

    def top(s):
        m = len(s) - 1
        out = tp[:m] * s[1:] % p
        out[1:] = (out[1:] - tm[1:m] * s[:m - 1]) % p
        return out
    S = [np.array([int(v) % p for v in seq], dtype=np.int64)]
    for _ in range(len(a) - 1):
        S.append(top(S[-1]))
    total = None
    for i, ai in enumerate(a):
        co = [int(v) for v in ai.coeffs()]
        if not co:
            continue
        r = co[-1] * S[i] % p
        for cj in reversed(co[:-1]):
            r = xop(r)
            r = (r + cj * S[i][:len(r)]) % p
        m = len(r) if total is None else min(len(total), len(r))
        total = r[:m] if total is None else (total[:m] + r[:m]) % p
    return total


# --------------------------------------------------------------------------------------- points (frozen classifier)
@lru_cache(maxsize=None)
def pts(k):
    p = int(k.split('|')[1])
    return C.points(op_of(k, p), p)


def wf(case, p, side):
    return C.weight_fixed(*(pts(key(case, p, side, sd)) for sd in SEEDS))


def fk(x0, p):
    return ((-I.to_fp(x0, p)) % p, 1)


# ================================================================================================= fast tier
def test_support_copies_are_the_frozen_instrument():
    """Wrong answer rejected: a pipeline (classifier, gcrd, search, engine) edited after record MO was scored."""
    for name, h in FROZEN.items():
        assert hashlib.sha256((SUP / name).read_bytes().replace(b'\r\n', b'\n')).hexdigest() == h, name


def test_engine_cases_are_record_AB():
    """Wrong answer rejected: an edited geometry or exponent set (the guard would test another integral); x* not the
    root of (D1^2 + D2^2) / (2 D1 D2)."""
    assert X.CASES['A1'] == dict(D1='9/5', D2='13/5', zB='13/10', zC='1', zA='9/10')
    assert X.CASES['A2'] == dict(D1='3/2', D2='12/5', zB='1', zC='6/5', zA='4/5')
    for case in CASES:
        D1, D2 = X.params(case)[:2]
        assert (D1 ** 2 + D2 ** 2) / (2 * D1 * D2) == XSTAR[case]


@pytest.mark.parametrize('case', CASES)
@pytest.mark.parametrize('p', PRIMES)
def test_embedded_sequences_match_engine_prefix(case, p):
    """The embedded generic sequences are AB's integral: their first L_PREFIX + 1 terms equal the frozen engine's
    (both seeds, base and shift; the shift is zeta_A + 1/997).  Wrong answers rejected: data from another geometry
    or prime; the shift not applied, or applied to another exponent; seeds swapped; a constant set other than the 21."""
    D1, D2, zB, zC, zA = X.params(case)
    K = X.E.Field(p)
    for side in SIDES:
        z = zA if side == 'base' else zA + SHIFT
        pc = X.pieces(case, L_PREFIX, K, params_override=(D1, D2, zB, zC, z))
        f = X.combine(pc, z, K)
        ks = sorted(set(k for q in f for k in q), key=str)
        seqs = {k: [int(q.get(k, flint.nmod(0, p))) for q in f] for k in ks}
        assert len(seqs) == 21
        for sd in SEEDS:
            e = D[key(case, p, side, sd)]
            assert e['zA'] == str(z) and sorted(map(repr, ks)) == e['keys']
            assert len(e['seq']) == 1001
            assert C.generic_seq(seqs, p, sd) == e['seq'][:L_PREFIX + 1], (case, p, side, sd)


@pytest.mark.parametrize('case', CASES)
@pytest.mark.parametrize('p', PRIMES)
def test_embedded_operators_annihilate_and_have_order_14_degree_111(case, p):
    """Claim (1) on the embedded operators, through an independent route: each annihilates its own sequence over the
    whole determined range (>= 700 terms), and has order 14, degree 111.  Negative control: the operator of the other
    weight vector at the other side does not annihilate this sequence (the check can fail).  Wrong answers rejected:
    a corrupted or truncated operator; an operator of another sequence; order or degree other than 14 / 111."""
    for side in SIDES:
        for sd in SEEDS:
            k = key(case, p, side, sd)
            c = op_of(k, p)
            assert len(c) - 1 == ORDER and max(x.degree() for x in c) == DEGREE, k
            r = apply_to_legendre(c, D[k]['seq'], p)
            assert len(r) >= 700 and not np.any(r), (k, len(r), int(np.count_nonzero(r)))
    other = op_of(key(case, p, 'shift', 7), p)
    r = apply_to_legendre(other, D[key(case, p, 'base', 2026)]['seq'], p)
    assert np.count_nonzero(r) > 600


def test_minimal_operator_recomputed_and_divides_ab_search():
    """Claims (1)-(2), recomputed at A2, prime 1073741441, base, seed 2026: the gcrd of the sequence's annihilators
    equals the embedded operator; AB's search returns order 17, degree 48 and is right-divisible by it
    (pseudo-remainder zero), but not by A1's minimal operator.  Wrong answers rejected: a minimal operator other than
    the embedded one (e.g. AB's order 17 taken as minimal); AB's search not a left multiple; a divisibility test that
    passes any pair."""
    p = PRIMES[0]
    k = key('A2', p, 'base', 2026)
    seq = D[k]['seq']
    mn = C.generic_operator(seq, p)
    assert (mn['order'], mn['degree']) == (ORDER, DEGREE)
    assert [[int(a) for a in ck.coeffs()] for ck in mn['c']] == D[k]['op']
    ab = C.generic_operator_ab(seq, p)
    assert (ab['order'], ab['degree']) == (AB_ORDER, AB_DEGREE)
    rem = I.right_rem(ab['c'], op_of(k, p), p)
    assert len(rem) == 1 and rem[0].degree() < 0
    rem2 = I.right_rem(ab['c'], op_of(key('A1', p, 'base', 2026), p), p)
    assert not (len(rem2) == 1 and rem2[0].degree() < 0)


@pytest.mark.parametrize('case', CASES)
@pytest.mark.parametrize('p', PRIMES)
def test_xstar_is_the_only_noninteger_point(case, p):
    """Claim (3), base and shift: among weight-fixed finite points only x* has a non-integer exponent (class M); its
    exponents are 0..10 and 3/2, 7/2, 9/2 under each weight; no point is irregular or has unaccounted exponents.
    Wrong answers rejected: a second class-M point (AB's reading for the order-17 operator); a 5/2 or other
    half-integer at x*; x* integer-only; exponents that differ between weights."""
    kx = fk(XSTAR[case], p)
    for side in SIDES:
        W = wf(case, p, side)
        ms = [v['pos'] for v in W.values() if v['cls'] == 'M']
        assert ms == [str(XSTAR[case])], (side, ms)
        for sd in SEEDS:
            e = pts(key(case, p, side, sd))[kx]
            assert sorted(Fraction(z) for z in e['exps']) == XSTAR_EXPS, (side, sd, e['exps'])
            assert e['modz'] == ['0', '1/2']
        for v in C.genuine(W).values():
            assert v['cls'] != 'X' and v['unaccounted'] == 0, (side, v['pos'])


@pytest.mark.parametrize('case', CASES)
@pytest.mark.parametrize('p', PRIMES)
def test_other_genuine_points_integer_with_forced_log(case, p):
    """Claim (4): every genuine weight-fixed point other than x* is +-1 or a linear point of class L (integer
    exponents, a logarithm forced) with exponents {-3, 0..12} or {-2, 0..12}; +-1 are class L.  Wrong answers
    rejected: a pole without logarithm (class P); a lowest exponent other than -3/-2; a degree > 1 genuine factor;
    +-1 missing or log-free."""
    kx, k1, km1 = fk(XSTAR[case], p), fk(1, p), fk(-1, p)
    for side in SIDES:
        W = C.genuine(wf(case, p, side))
        assert W[k1]['cls'] == 'L' and W[km1]['cls'] == 'L', side
        extra = {k: v for k, v in W.items() if k not in (kx, k1, km1)}
        assert extra, side
        for v in extra.values():
            assert v['deg'] == 1 and v['cls'] == 'L', (side, v['pos'], v['cls'])
            ex = sorted(int(z) for z in v['exps'])
            assert ex in (EXTRA_EXPS[-3], EXTRA_EXPS[-2]), (side, v['pos'], ex)


@pytest.mark.parametrize('case', CASES)
@pytest.mark.parametrize('p', PRIMES)
def test_minus3_points_move_with_zeta_A_minus2_point_fixed(case, p):
    """Claim (5): read with the record's MO2 rule (C.mo2_points), every exponent -3 point moves under the zeta_A
    shift, the exponent -2 point (one per side) does not, the moving count is the record's, and MO2 holds (moving
    points integer-only).  Wrong answers rejected: no point moving (BORDERLINE); a -2 point that moves; a moving point
    with a non-integer exponent (STOP); a count other than the record's 4, 2, 2, 4; A2's fixed point elsewhere."""
    kx, k1, km1 = fk(XSTAR[case], p), fk(1, p), fk(-1, p)
    rows = C.mo2_points(wf(case, p, 'base'), wf(case, p, 'shift'), exclude=(kx, k1, km1))
    low = lambda r: min(int(z) for z in r['exps'])                            # noqa: E731
    assert all(r['moves'] == (low(r) == -3) for r in rows), rows
    fixed = [r for r in rows if not r['moves']]
    assert len(fixed) == 1 and low(fixed[0]) == -2
    assert sum(r['moves'] for r in rows) == N_MOVING[(case, p)]
    assert C.mo2_holds(rows)
    if case == 'A2':
        assert fixed[0]['pos'] == A2_FIXED


def test_mo2_reader_sees_a_moving_noninteger_point():
    """Instrument control (record MO's T7, rebuilt here): two toy Legendre sequences t^L and (2L+1) s^L; shifting
    t = 1/3 -> 3/10 moves x_t, which carries exponent 1/2 mod 1, and the MO2 reader must report it moving and NOT
    integer-only.  Plus T6: e_L = t^(L+1)/(L+1) has logarithmic (class L) points at +-1.  Wrong answers rejected: a
    reader that calls every moving point integer-only or never reports a move (the -3 test would then be vacuous);
    a log test that never fires."""
    p = PRIMES[0]
    n = 400
    t, s, t2 = (I.to_fp(Fraction(*v), p) for v in ((1, 3), (2, 5), (3, 10)))
    a = [pow(t, L, p) for L in range(n)]
    c = [(2 * L + 1) * pow(s, L, p) % p for L in range(n)]
    a2 = [pow(t2, L, p) for L in range(n)]

    def W(sq):
        return C.weight_fixed(*[C.points(C.generic_operator(C.generic_seq(sq, p, sd), p)['c'], p) for sd in SEEDS])
    rows = C.mo2_points(W({'a': a, 'c': c}), W({'a': a2, 'c': c}), exclude=())
    xt, xt2 = (1 + Fraction(1, 3) ** 2) / Fraction(2, 3), (1 + Fraction(3, 10) ** 2) / Fraction(6, 10)
    assert sorted((r['side'], r['pos'], r['integer_only']) for r in rows if r['moves']) == \
        [('base', str(xt), False), ('shift', str(xt2), False)]
    assert not C.mo2_holds(rows)
    e = [pow(t, L + 1, p) * pow(L + 1, p - 2, p) % p for L in range(n)]
    pe = C.genuine(C.points(C.piece_operator(e, p)['c'], p))
    assert sorted((v['pos'], v['cls']) for v in pe.values()) == [('-1', 'L'), ('1', 'L'), (str(xt), 'M')]


# ================================================================================================= slow tier
@pytest.mark.slow
@pytest.mark.parametrize('side', SIDES)
@pytest.mark.parametrize('case', CASES)
@pytest.mark.parametrize('p', PRIMES)
def test_embedded_sequences_full_recompute(case, p, side):
    """All 1001 terms of both seeds' sequences from the frozen engine (~11 min each).  Wrong answer rejected:
    embedded data that agree on a prefix only."""
    D1, D2, zB, zC, zA = X.params(case)
    z = zA if side == 'base' else zA + SHIFT
    K = X.E.Field(p)
    f = X.combine(X.pieces(case, 1000, K, params_override=(D1, D2, zB, zC, z)), z, K)
    ks = sorted(set(k for q in f for k in q), key=str)
    seqs = {k: [int(q.get(k, flint.nmod(0, p))) for q in f] for k in ks}
    for sd in SEEDS:
        assert C.generic_seq(seqs, p, sd) == D[key(case, p, side, sd)]['seq']


@pytest.mark.slow
@pytest.mark.parametrize('side', SIDES)
@pytest.mark.parametrize('case', CASES)
@pytest.mark.parametrize('p', PRIMES)
def test_minimal_operator_recomputed_all(case, p, side):
    """Claims (1)-(2) for every entry (~20 s per gcrd): gcrd recomputed equals the embedded operator; on the base side
    AB's search (order 17, degree 48) is right-divisible by it."""
    for sd in SEEDS:
        k = key(case, p, side, sd)
        mn = C.generic_operator(D[k]['seq'], p)
        assert [[int(a) for a in ck.coeffs()] for ck in mn['c']] == D[k]['op'], k
        ab = C.generic_operator_ab(D[k]['seq'], p)
        assert (ab['order'], ab['degree']) == (AB_ORDER, AB_DEGREE), k
        rem = I.right_rem(ab['c'], op_of(k, p), p)
        assert len(rem) == 1 and rem[0].degree() < 0, k


# Fire-test record (python debug/qa/fire_test.py tests/test_paper59_aa_bc_minimal_operator.py [--plant-in <file>]
# --plant "O=>N" -k <sel>): fire-tested 2026-11-17, every plant FIRED, every plant restored.  -k selections leave out
# the hash test, so each plant is caught by the behavioural test named, not by the sha256 check.
#   F01 this file  ORDER, DEGREE 14, 111 => 17, 48 (AB's search as minimal)   -k "annihilate or recomputed_and"  FIRED
#   F02 this file  x* exponent 9/2 => 5/2 (AB's 5/2 reading)                  -k xstar                           FIRED
#   F03 this file  extra-point lowest exponent -3 => -4                       -k other_genuine                   FIRED
#   F04 this file  n_moving (A1, 1073741441) 4 => 2                           -k minus3                          FIRED
#   F05 this file  A2 fixed point 2509/1980 => 2509/1981                      -k minus3                          FIRED
#   F06 this file  zeta_A shift 1/997 => 1/998                                -k prefix                          FIRED
#   F07 this file  A2 x* 89/80 => 89/81                                       -k xstar                           FIRED
#   F08 data json  one sequence term inside the prefix (+1)                   -k prefix                          FIRED
#   F09 data json  one sequence term beyond the prefix, L = 500 (+1)          -k annihilate                      FIRED
#   F10 data json  one operator coefficient (+1)                              -k annihilate                      FIRED
#   F11 s1_controls.py  mo2_points mv = k not in B => mv = False (nothing moves)   -k minus3                     FIRED
#   F12 s1_controls.py  integer_only => True (reader blind to non-integer moves)   -k mo2_reader                 FIRED
#   F13 s1_controls.py  forced-log class 'L' => 'P'                                -k other_genuine              FIRED
#   F14 intclass_lib.py  minimal_operator gcrd step dropped (search op returned)  -k recomputed_and             FIRED
#   F15 aabc_exact.py  A2 zA '4/5' => '4/6'                                       -k prefix                     FIRED
#   F16 s1_controls.py  docstring one-space edit                                  -k frozen                     FIRED
# Not planted: dropping the weight_fixed filter.  For the minimal operator the weight-moving factors are all apparent
# (class A; checked 2026-11-17 at all 8 (case, prime, side): genuine sets with and without the filter are equal), so
# that filter is not load-bearing for these claims and a plant on it would not fire.
