"""Inventory row 4: the one-electron three-centre integral NA, as a function of the bond angle, is NOT TRIANGULAR
(it leaves the two-centre-type class) at both geometries N1 and N2.  Regression guard, written 2026-11-16 as a
separate step after the records (CLAUDE.md Sec. 9, guard-writing rule).

Records: docs/predictions/2026-11-01_integrated_class.md (IC: N1 decided, its Caps erratum + addendum) and
docs/predictions/2026-11-15_na_n2_decided.md (IN: N2 decided on primes 1073741717 and 1073741689 under the
erratum check, GO).

Instrument: verbatim copies of the frozen instrument in tests/paper59_intclass_support/ (sha256 checked below;
na_exact.py and intclass_in_check.py were added 2026-11-16 from commit 247836bc, the IN registration).  Nothing
is imported from debug/.

What is recomputed, on primes NOT used by either record (both = 1 mod 4, with 11 and -69 squares, so N2's and
N1's lambda edges split):
  * fast: both geometries, all six constant components each, at L = 300 on 1073741857 (~40 s).  L = 300 is the
    smallest tried that reproduces the record's operator search exactly (R0 and the (Q0, nullity) table at
    R0..R0+3, as at L = 1200); at L = 240 the order-6 operators of the genus-one components are not found
    (R0 = 7), at L = 200 a component has no set-B operator.  Below L = 300 the guard would test a different object.
  * slow: the same at the records' L = 1200 on 1073741909 (~800 s).
Per component: minimal operator -> pullback to w = |Y - Z| -> order-one peel: both operator sets (A, B) STUCK
with no peeled step, clean (no flag), the same residual order as scored; the erratum check (IN) passes for both
sets; independently of that check, every indicial exponent at every F_p-irreducible factor of the leading
coefficient is an integer in F_p (so the search's blind spot, a degree >= 2 factor with a non-integer class or
exponents outside F_p, is empty), and at least one factor of degree >= 2 is examined (the check is not vacuous);
the rates at infinity are 0, (+-zH), and the pair +-sqrt(zH^2 - zX^2) found in F_p (the split that IC's primes
could not see at N2).  STUCK (no order-one right factor) excludes TRIANGULAR on its own: the bottom of a
composition series of order-one factors is an order-one right factor.  Many residuals here DO have an order-one
LEFT factor (adjoint_has_order1 True); that is not asserted either way and does not touch the verdict.
Controls: the erratum-check controls of record IN (planted Chebyshev with c a non-square: STUCK and the check
FAILS; with c a square: peels; Legendre K: STUCK and PASSES; Euler w^{+-sqrt c}: STUCK and FAILS; Bessel K0:
STUCK and PASSES), and a pipeline control: rho^L/(L+1) (a logarithm) through the same operator_sets + pullback
+ peel path returns TRIANGULAR, so the path is not STUCK-always.
Scope (IC's caps): one prime per run, function-level class at fixed rational lengths and exponents, 1s
functions, Qbar-independence of the six constants and minimality of the operators assumed as in the records.
"""
import hashlib
import importlib.util
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import flint
import pytest
from sympy import isprime, legendre_symbol

HERE = Path(__file__).parent
SUP = HERE / 'paper59_intclass_support'
# sha256 after CRLF -> LF.  The first four as in tests/test_paper59_integrated_class.py; na_exact.py and
# intclass_in_check.py equal debug/ at commit 247836bc (record IN's registration).
FROZEN = {'intclass_lib.py': '2a142e5b5667b432c7de4be885f5b6767daf3d209c86261076dc6ae5e4f641d4',
          'intclass_score.py': 'aaf700e482f53c0a204939c08e6e12c9786d544e48c77e176503f21e5add91ca',
          'na_modguess.py': '84691c87a4a280b6823438af574c0fff6af9bf3f7ba7b76d93f2bbe866847d7d',
          'na_modanalyze.py': 'de92dc3acdcc7082001e32b50ef44db8c3bd75d990f615603b5f97f6ea178005',
          'na_exact.py': '85fc39f18bfe9231f47925436b72258dca7afddc6f9619257d9fafc708004afb',
          'intclass_in_check.py': '50fb31871454014a1f2dee653985c778a8e276afbe478a66eecb6e2212fbf4c4'}

RECORD_PRIMES = {1073741789, 1073741833, 1073741717, 1073741689}     # IC's two, IN's two
PRIME_FAST, L_FAST = 1073741857, 300
PRIME_SLOW, L_SLOW = 1073741909, 1200


def _load():
    sys.path.insert(0, str(SUP))
    try:
        mods = []
        for name in ('intclass_in_check', 'na_exact'):
            spec = importlib.util.spec_from_file_location(f'p59na_{name}', SUP / f'{name}.py')
            m = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(m)
            mods.append(m)
    finally:
        sys.path.remove(str(SUP))
    C, NX = mods
    return C, C.I, C.S, NX


C, I, S, NX = _load()

# The records' scored table (debug/data/intclass_in_<case>_<p>.json, identical on both IN primes and at IC):
# component key -> (R0, {R: (Q0, nullity)} for R0..R0+3, residual order, residual type).
G1 = 'confluent-genus1-type'
EXPECT = {
    'N2': {"(Fraction(-177, 25), '1', None)": (6, {6: (26, 1), 7: (18, 2), 8: (15, 2), 9: (14, 4)}, 6, G1),
           "(Fraction(-3, 2), '1', None)": (5, {5: (18, 1), 6: (13, 2), 7: (11, 2), 8: (10, 2)}, 5, 'other'),
           "(Fraction(-87, 25), '1', None)": (6, {6: (26, 1), 7: (18, 2), 8: (15, 2), 9: (14, 4)}, 6, G1),
           "(Fraction(-9, 5), '1', None)": (5, {5: (15, 1), 6: (11, 2), 7: (9, 2), 8: (8, 2)}, 5, 'other'),
           "(Fraction(-9, 5), 'E1', Fraction(132, 25))": (4, {4: (14, 1), 5: (10, 2), 6: (8, 1), 7: (8, 4)}, 4,
                                                          'other'),
           "(Fraction(9, 5), 'E1', Fraction(132, 25))": (4, {4: (14, 1), 5: (10, 2), 6: (8, 1), 7: (8, 4)}, 4,
                                                         'other')},
    'N1': {"(Fraction(-117, 50), '1', None)": (5, {5: (18, 1), 6: (13, 2), 7: (11, 2), 8: (10, 2)}, 5, 'other'),
           "(Fraction(-209, 50), '1', None)": (6, {6: (26, 1), 7: (18, 2), 8: (15, 2), 9: (14, 4)}, 6, G1),
           "(Fraction(-389, 50), '1', None)": (6, {6: (26, 1), 7: (18, 2), 8: (15, 2), 9: (14, 4)}, 6, G1),
           "(Fraction(-9, 5), '1', None)": (5, {5: (15, 1), 6: (11, 2), 7: (9, 2), 8: (8, 2)}, 5, 'other'),
           "(Fraction(-9, 5), 'E1', Fraction(299, 50))": (4, {4: (14, 1), 5: (10, 2), 6: (8, 1), 7: (8, 4)}, 4,
                                                          'other'),
           "(Fraction(9, 5), 'E1', Fraction(299, 50))": (4, {4: (14, 1), 5: (10, 2), 6: (8, 1), 7: (8, 4)}, 4,
                                                         'other')}}
# The records' geometries (D1, D2, zX, zH), restated here so that the instrument's table is checked, not trusted.
GEOM = {'N1': ('1.8', '2.6', '1.3', '1.0'), 'N2': ('1.5', '2.4', '1.0', '1.2')}


def test_support_copies_are_the_frozen_instrument():
    """Wrong answer rejected: an instrument edited after the records were scored."""
    for name, h in FROZEN.items():
        assert hashlib.sha256((SUP / name).read_bytes().replace(b'\r\n', b'\n')).hexdigest() == h, name


def test_instrument_geometries_are_the_records():
    """Wrong answer rejected: a swapped or edited geometry (the guard would then test another integral)."""
    assert S.NACASES == GEOM


@pytest.mark.parametrize('p', [PRIME_FAST, PRIME_SLOW])
def test_prime_class(p):
    """Wrong answers rejected: a prime already used by a record (not an independent run); a prime where 11 or -69
    is a non-residue (the lambda edge then fails to split and every component is flagged, IC's N2 failure);
    p = 3 mod 4."""
    assert isprime(p) and p % 4 == 1 and p not in RECORD_PRIMES
    assert legendre_symbol(11, p) == 1 and legendre_symbol(-69 % p, p) == 1


@lru_cache(maxsize=None)
def coefficients(case, p, L):
    D1s, D2s, zX, zH = S.NACASES[case]
    Qs = NX.coefficients(D1s, D2s, zX, zH, L, NX.Field(p))
    keys = sorted(set(k for q in Qs for k in q), key=repr)
    return {repr(k): [int(q.get(k, flint.nmod(0, p))) for q in Qs] for k in keys}


def check_case(case, p, L):
    D1s, D2s, zXs, zHs = S.NACASES[case]
    D1, D2, zX, zH = (Fraction(v) for v in (D1s, D2s, zXs, zHs))
    comps = coefficients(case, p, L)
    exp = EXPECT[case]
    assert set(comps) == set(exp), (case, sorted(comps))
    a0 = I.to_fp((D1 ** 2 + D2 ** 2) / (2 * D1 * D2), p)
    a2 = I.to_fp(Fraction(-1) / (2 * D1 * D2), p)
    lam2 = I.to_fp(zH ** 2 - zX ** 2, p)
    small = {0, I.to_fp(zH, p), I.to_fp(-zH, p)}
    for key, seq in comps.items():
        R0, found, order, typ = exp[key]
        ops = S.operator_sets(seq, p, 24)
        assert ops is not None, (case, key, 'no operator')
        assert ops['R0'] == R0 and ops['found'] == found, (case, key, ops['R0'], ops['found'])
        rA = C.analyse_full(ops['A'], D1, D2, p)
        rB = C.analyse_full(ops['B'], D1, D2, p)
        # the record's verdict path
        assert S.component_verdict(rA, rB) == 'STUCK', (case, key)
        for tag, r in (('A', rA), ('B', rB)):
            assert r['verdict'] == 'STUCK' and r['steps'] == 0, (case, key, tag, r['verdict'], r['steps'])
            res = r['residual']
            assert res['clean'] and res['flags'] == [], (case, key, tag, res['flags'])
            assert res['order'] == order and res['type'] == typ, (case, key, tag, res['order'], res['type'])
            assert r['erratum']['passed'] is True and r['erratum']['bad'] == [], (case, key, tag, r['erratum'])
        # independent of erratum_check: the search's blind spot is empty, and the check had something to look at
        for tag in ('A', 'B'):
            pk = I.peel(I.pullback_quadratic(ops[tag], a0, a2, p), p, verbose=False)
            assert pk['verdict'] == 'STUCK'
            pts = I.finite_points(pk['residual'], p)
            assert pts, (case, key, tag, 'no finite singular factor')
            for pt in pts:
                assert pt['unaccounted'] == 0 and pt['croots'], (case, key, tag, pt['fdeg'])
                assert all(k == 'int' for (k, v), m in pt['croots']), (case, key, tag, pt['croots'])
            assert max(pt['fdeg'] for pt in pts) >= 2, (case, key, tag, 'no factor of degree >= 2')
            inf = pk['residual_diag']['inf']
            assert inf['rank_le_1'] and inf['lam_unsplit'] == 0, (case, key, tag, inf)
            lams = [r for r, m in inf['lam_roots']]
            pair = [x for x in lams if x not in small]
            assert 0 in lams and len(pair) == 2 and all(x * x % p == lam2 for x in pair), (case, key, tag, lams)
            assert (pair[0] + pair[1]) % p == 0
            assert (typ == G1) == (I.to_fp(zH, p) in lams and I.to_fp(-zH, p) in lams), (case, key, tag, lams)


@pytest.mark.parametrize('case', ['N2', 'N1'])
def test_na_not_triangular_fast(case):
    """Wrong answers rejected: any component that peels (TRIANGULAR, or STUCK after a peeled step); a flagged
    (unclean) STUCK, as at IC's N2 on primes where 11 is a non-residue; sets A and B disagreeing; a residual
    order, operator order or search table different from the records'; an erratum check that fails or finds
    anything; a residual with a non-integer or out-of-F_p exponent at any factor; a vacuous check (no factor of
    degree >= 2); rates at infinity not the records' (0, +-zH, +-sqrt(zH^2 - zX^2)); a component set different from
    the records' six (a geometry change or a dropped/added constant)."""
    check_case(case, PRIME_FAST, L_FAST)


@pytest.mark.slow
@pytest.mark.parametrize('case', ['N2', 'N1'])
def test_na_not_triangular_record_L(case):
    """As the fast test, at the records' L = 1200, on a second new prime."""
    check_case(case, PRIME_SLOW, L_SLOW)


def _nmod(lst, p):
    return [flint.nmod_poly([x % p for x in c], p) for c in lst]


def _peel_check(op, p):
    r = I.peel(op, p, verbose=False)
    return r['verdict'], (C.erratum_check(r['residual'], p)['passed'] if r['verdict'] == 'STUCK' else None)


def test_erratum_check_controls():
    """Record IN's erratum-check controls at PRIME_FAST.  Wrong answers rejected: an erratum check that passes
    everything (Chebyshev non-square and Euler would pass), or fails everything (Legendre K and Bessel K0 would
    fail); a search that misses the square-c Chebyshev peel.  The verdicts are asserted strictly: a control
    returning TRIANGULAR where STUCK is expected skips the check and fails here."""
    p = PRIME_FAST
    cn = next(c for c in range(2, 200) if legendre_symbol(c, p) == -1)
    cs = next(c for c in range(2, 200) if legendre_symbol(c, p) == 1 and int(c ** 0.5) ** 2 != c)
    cheb = lambda c: _nmod([[1], [0, -4], [4 * c, 0, -4]], p)          # 4(c - w^2) y'' - 4 w y' + y
    assert _peel_check(cheb(cn), p) == ('STUCK', False), cn
    assert _peel_check(cheb(cs), p)[0] == 'TRIANGULAR', cs
    q4 = (-pow(4, p - 2, p)) % p
    legendre_k = _nmod([[q4], [0, -2], [cn, 0, -1]], p)                 # (c - w^2) y'' - 2 w y' - y/4
    assert _peel_check(legendre_k, p) == ('STUCK', True)
    euler = _nmod([[-cn], [0, 1], [0, 0, 1]], p)                        # w^2 y'' + w y' - c y
    assert _peel_check(euler, p) == ('STUCK', False)
    bessel_k0 = _nmod([[0, -1], [1], [0, 1]], p)                        # w y'' + y' - w y
    assert _peel_check(bessel_k0, p) == ('STUCK', True)


def test_pipeline_can_return_triangular():
    """rho^L/(L+1) (a logarithm, triangular) through operator_sets -> pullback -> peel at N2's geometry.  Wrong
    answer rejected: a pipeline that returns STUCK for everything (the NA verdicts would then mean nothing)."""
    p = PRIME_FAST
    D1, D2 = (Fraction(v) for v in S.NACASES['N2'][:2])
    rho = D1 / D2
    seq = [I.to_fp(rho ** L / (L + 1), p) for L in range(300)]
    ops = S.operator_sets(seq, p, 24)
    assert ops is not None
    rA = C.analyse_full(ops['A'], D1, D2, p)
    rB = C.analyse_full(ops['B'], D1, D2, p)
    assert rA['verdict'] == 'TRIANGULAR' and S.component_verdict(rA, rB) == 'TRIANGULAR'


# Fire-test record (python debug/qa/fire_test.py tests/test_paper59_na_class.py --plant-in <file> --plant "O=>N"):
# fire-tested 2026-11-16, every plant FIRED, every plant restored (hash test green afterwards).  The -k selection
# leaves out the hash test, so each plant is caught by the behavioural test named, not by the sha256 check.
#   (a)  intclass_in_check.py  return dict(passed=not bad, ...) => passed=True     -k erratum_check_controls  FIRED
#   (a') intclass_in_check.py  same                             => passed=False    -k fast                    FIRED
#   (a2) intclass_in_check.py  clause (a) dropped (if nonint or pt['unaccounted'] > 0 => if nonint)
#                                                       -k erratum_check_controls (Euler control)  FIRED
#   (a3) intclass_in_check.py  clause (b) widened to fdeg >= 3  -k erratum_check_controls (Chebyshev)  FIRED
#   (b)  intclass_lib.py  while len(cur) - 1 >= 2 => >= 7 (peel returns TRIANGULAR at once)  -k fast   FIRED
#   (b') intclass_lib.py  res = dict(verdict='STUCK' => verdict='TRIANGULAR'                 -k fast   FIRED
#   (b") intclass_lib.py  if sol is None: => if True: (peel never peels)                    -k pipeline FIRED
#   (c)  this file  PRIME_FAST 1073741857 => 1073741621 (11 a non-residue): -k "fast and N2" FIRED (the flagged
#        residuals give component_verdict STUCK-UNCONFIRMED, not STUCK); -k prime_class FIRED
#   (d)  intclass_score.py  N2 geometry zX, zH swapped ('1.0', '1.2' => '1.2', '1.0')        -k "fast and N2" FIRED
#   (d') na_exact.py  the x i_{L+1} term dropped (K(2L+1) iL + K(x) iL1 => K(2L+1) iL)    -k fast   FIRED
#   (e)  na_exact.py  one-character edit (def zero => def zer0)                            -k frozen FIRED
