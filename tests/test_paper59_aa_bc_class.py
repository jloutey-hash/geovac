"""Inventory row 6: the screened (mu-carrying) part of (1s_A 1s_A | 1s_B 1s_C), as a function of the bond angle at
case A1, leaves the triangular (two-centre) class on its own.  Regression guard for record AK
(docs/predictions/2026-11-16_aa_bc_class_decided.md, scored GO), written 2026-11-16 as a separate step after the
record (CLAUDE.md Sec. 9, guard-writing rule).  Paper 59 sec:scope, the (AA|BC) paragraph.

Claims guarded:
(1) the four unsplit rates at infinity are the roots of lambda^4 - (51/20) lambda^2 + 81/25, namely
    (+-sqrt(615) +- i sqrt(105))/20; they all lie in F_p iff 615 and -105 are squares mod p; record AB's primes
    (1073741789, 1073741833) fail that test and record AK's (1073741441, 1073741381) pass it.
(2) on primes that split the rates, all 15 mu-carrying components of A1 are clean STUCK from both operator sets
    (IC's verdict), with IN's erratum check passing, residual orders 4 (12 components), 5 (e^{-162/25}),
    6 (e^{-189/50}, e^{-513/50}), and no lambda-edge flag.
(3) the rates the instrument finds at infinity satisfy the quartic; the remaining rates are 0 (the order-5
    component) and +-9/5 = +-mu (the order-6 pair); the order-4 components carry the quartic's four only.

Instrument: verbatim copies in tests/paper59_intclass_support/ (sha256 checked below); aabc_exact.py was added
2026-11-16 from debug/aabc_exact.py, equal to its version at record AB's commit 3a2fbfcf.  Nothing is imported from
debug/.  The scored table EXPECT below is copied from debug/data/aabc_ak_<p>.json (identical on both AK primes).

Recomputed, on primes NOT used by AB or AK (both = 1 mod 4 with 615 and -105 squares):
  * fast: all 15 components at L = 400 on 1073741477 (~2 min; coefficients ~110 s).  L = 300 does not reproduce
    the record's search for the order-6 pair (R0 = 7 instead of 6: too few rows for Q0 = 41 at R = 6); L = 378,
    380 and 400 reproduce all 15 tables exactly (measured 2026-11-16); 400 keeps a margin.
  * slow: all 15 at the record's L = 1200 on 1073740697 (~19 min; coefficients ~970 s).
Controls: IN's erratum-check controls (planted Chebyshev with c a non-square: STUCK and the check FAILS; c a square:
peels; Legendre K: STUCK, PASSES; Euler w^{+-sqrt c}: STUCK, FAILS; Bessel K0: STUCK, PASSES), and rho^L/(L+1)
through the same path at A1's geometry returns TRIANGULAR (the path is not STUCK-always).
Scope (IC's and AK's caps): one geometry (A1), one prime per run; Qbar-independence of the constants and minimality
of the operators assumed as in the records; a clean STUCK is evidence modulo the instrument's search, not a proof
over Qbar(x).  The rates' origin is not claimed and not tested.
"""
import hashlib
import importlib.util
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import flint
import pytest
import sympy as sp
from sympy import isprime, legendre_symbol, primerange

HERE = Path(__file__).parent
SUP = HERE / 'paper59_intclass_support'
# sha256 after CRLF -> LF.  First six as in tests/test_paper59_na_class.py; aabc_exact.py = debug/ at commit 3a2fbfcf.
FROZEN = {'intclass_lib.py': '2a142e5b5667b432c7de4be885f5b6767daf3d209c86261076dc6ae5e4f641d4',
          'intclass_score.py': 'aaf700e482f53c0a204939c08e6e12c9786d544e48c77e176503f21e5add91ca',
          'na_modguess.py': '84691c87a4a280b6823438af574c0fff6af9bf3f7ba7b76d93f2bbe866847d7d',
          'na_modanalyze.py': 'de92dc3acdcc7082001e32b50ef44db8c3bd75d990f615603b5f97f6ea178005',
          'na_exact.py': '85fc39f18bfe9231f47925436b72258dca7afddc6f9619257d9fafc708004afb',
          'intclass_in_check.py': '50fb31871454014a1f2dee653985c778a8e276afbe478a66eecb6e2212fbf4c4',
          'aabc_exact.py': '218b7283f525d25ed8a87df61f0749ca6f576b20abcc4bde1272de32436b2d30'}

AB_PRIMES = (1073741789, 1073741833)
AK_PRIMES = (1073741441, 1073741381)
PRIME_FAST, L_FAST = 1073741477, 400
PRIME_SLOW, L_SLOW = 1073740697, 1200

# The claimed quartic, monic, coefficients of lambda^4 .. lambda^0.
QUARTIC = (Fraction(1), Fraction(0), Fraction(-51, 20), Fraction(0), Fraction(81, 25))
MU = Fraction(9, 5)                                                     # 2 zeta_A at A1


def _load():
    sys.path.insert(0, str(SUP))
    try:
        mods = []
        for name in ('intclass_in_check', 'aabc_exact'):
            spec = importlib.util.spec_from_file_location(f'p59ak_{name}', SUP / f'{name}.py')
            m = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(m)
            mods.append(m)
    finally:
        sys.path.remove(str(SUP))   # one entry only: the instrument's own inserts stay (operator_sets imports late)
    C, X = mods
    return C, C.I, C.S, X


C, I, S, X = _load()

G1 = 'confluent-genus1-type'
_O4 = (4, {4: (28, 1), 5: (18, 2), 6: (14, 1), 7: (13, 4)}, 4)
_O6 = (6, {6: (41, 1), 7: (26, 1), 8: (21, 1), 9: (19, 3)}, 6)
# component key -> (R0, {R: (Q0, nullity)} for R0..R0+3, residual order); every residual type is G1.
EXPECT = {
    "(Fraction(-162, 25), '1', None)": (5, {5: (31, 1), 6: (20, 2), 7: (16, 3), 8: (14, 4)}, 5),
    "(Fraction(-162, 25), 'E1', Fraction(-27, 10))": _O4,
    "(Fraction(-162, 25), 'E1', Fraction(13, 10))": _O4,
    "(Fraction(-162, 25), 'E1', Fraction(189, 50))": _O4,
    "(Fraction(-162, 25), 'ln', Fraction(1, 2))": _O4,
    "(Fraction(-162, 25), 'ln', Fraction(21, 10))": _O4,
    "(Fraction(-162, 25), 'ln', Fraction(3, 2))": _O4,
    "(Fraction(-162, 25), 'ln', Fraction(41, 10))": _O4,
    "(Fraction(-189, 50), '1', None)": _O6,
    "(Fraction(-513, 50), '1', None)": _O6,
    "(Fraction(-72, 25), 'E1', Fraction(13, 10))": _O4,
    "(Fraction(-72, 25), 'E1', Fraction(369, 50))": _O4,
    "(Fraction(-72, 25), 'E1', Fraction(9, 10))": _O4,
    "(Fraction(162, 25), 'E1', Fraction(533, 50))": _O4,
    "(Fraction(72, 25), 'E1', Fraction(533, 50))": _O4,
}
# Rates at infinity besides the quartic's four, by residual order (record: 0 at order 5, +-mu at order 6).
EXTRA = {4: set(), 5: {Fraction(0)}, 6: {MU, -MU}}
A1 = {'D1': '9/5', 'D2': '13/5', 'zB': '13/10', 'zC': '1', 'zA': '9/10'}     # record AB's case A1, restated


def test_support_copies_are_the_frozen_instrument():
    """Wrong answer rejected: an instrument or engine edited after the records were scored."""
    for name, h in FROZEN.items():
        assert hashlib.sha256((SUP / name).read_bytes().replace(b'\r\n', b'\n')).hexdigest() == h, name


def test_engine_case_is_record_A1():
    """Wrong answer rejected: an edited A1 geometry or exponents (the guard would test another integral); an engine
    that picks up some other na_exact than the frozen copy."""
    assert X.CASES['A1'] == A1
    assert Path(X.E.__file__).resolve().parent == SUP.resolve()


def _quartic_roots_mod(p):
    poly = flint.nmod_poly([int(I.to_fp(c, p)) for c in reversed(QUARTIC)], p)
    return sorted(int(r) for r, m in poly.roots() for _ in range(m))


def test_quartic_roots_closed_form():
    """Claim (1), exact.  Wrong answers rejected: changed coefficients (51/20, 81/25, or a nonzero odd term); roots
    other than (+-sqrt 615 +- i sqrt 105)/20."""
    lam = sp.Symbol('lam')
    q = sum(sp.Rational(c.numerator, c.denominator) * lam ** (4 - i) for i, c in enumerate(QUARTIC))
    roots = [(s1 * sp.sqrt(615) + s2 * sp.I * sp.sqrt(105)) / 20 for s1 in (1, -1) for s2 in (1, -1)]
    assert sp.expand(sp.prod([lam - r for r in roots]) - q) == 0
    # the record notes |lambda|^2 = 9/5 for all four (a reading; its origin is not claimed)
    assert all(sp.simplify(r * sp.conjugate(r)) == sp.Rational(9, 5) for r in roots)


def test_split_iff_615_and_minus105_squares():
    """Claim (1), mod p.  The quartic has four roots in F_p iff 615 and -105 are both squares: checked on every prime
    in [1000, 30000) and on the record primes.  Wrong answers rejected: a criterion with another pair of numbers;
    AB's primes splitting (the flag at P4 would then be unexplained); AK's primes not splitting."""
    for p in list(primerange(1000, 30000)) + list(AB_PRIMES) + list(AK_PRIMES) + [PRIME_FAST, PRIME_SLOW]:
        crit = legendre_symbol(615, p) == 1 and legendre_symbol(-105 % p, p) == 1
        assert (len(_quartic_roots_mod(p)) == 4) == crit, p
    assert all(len(_quartic_roots_mod(p)) < 4 for p in AB_PRIMES)
    assert all(len(_quartic_roots_mod(p)) == 4 for p in AK_PRIMES)


@pytest.mark.parametrize('p', [PRIME_FAST, PRIME_SLOW])
def test_prime_class(p):
    """Wrong answers rejected: a prime already used by AB or AK (not an independent run); a prime on which the rates
    do not split (every component would carry the lambda-edge flag, AB's P4 failure); p = 3 mod 4."""
    assert isprime(p) and p % 4 == 1 and p not in AB_PRIMES + AK_PRIMES
    assert legendre_symbol(615, p) == 1 and legendre_symbol(-105 % p, p) == 1


@lru_cache(maxsize=None)
def components(p, L):
    """As debug/aabc_ak_check.components: the coefficients of A1's f minus NA's constants, mod p."""
    K = X.E.Field(p)
    pc = X.pieces('A1', L, K)
    f = X.combine(pc, X.params('A1')[4], K)
    na_keys = set(k for q in pc['NA'] for k in q)
    keys = sorted(set(k for q in f for k in q) - na_keys, key=repr)
    return {repr(k): [int(q.get(k, flint.nmod(0, p))) for q in f] for k in keys}


def check_all(p, L):
    D1, D2 = X.params('A1')[:2]
    comps = components(p, L)
    assert set(comps) == set(EXPECT), sorted(comps)
    a0 = I.to_fp((D1 ** 2 + D2 ** 2) / (2 * D1 * D2), p)
    a2 = I.to_fp(Fraction(-1) / (2 * D1 * D2), p)
    quart = set(_quartic_roots_mod(p))
    orders = []
    for key, seq in comps.items():
        R0, found, order = EXPECT[key]
        ops = S.operator_sets(seq, p, 24)
        assert ops is not None, (key, 'no operator')
        assert ops['R0'] == R0 and ops['found'] == found, (key, ops['R0'], ops['found'])
        rA = C.analyse_full(ops['A'], D1, D2, p)
        rB = C.analyse_full(ops['B'], D1, D2, p)
        v = S.component_verdict(rA, rB)                                 # the record's verdict path
        assert v == 'STUCK', (key, v, [((r.get('residual') or {}).get('flags')) for r in (rA, rB)])
        for tag, r in (('A', rA), ('B', rB)):
            assert r['verdict'] == 'STUCK' and r['steps'] == 0, (key, tag, r['verdict'], r['steps'])
            res = r['residual']
            assert res['clean'] and res['flags'] == [], (key, tag, res['flags'])
            assert res['order'] == order and res['type'] == G1, (key, tag, res['order'], res['type'])
            assert r['erratum']['passed'] is True and r['erratum']['bad'] == [], (key, tag, r['erratum'])
        orders.append(order)
        # claim (3): the rates at infinity, read from the peel's residual directly (not through the flag)
        for tag in ('A', 'B'):
            pk = I.peel(I.pullback_quadratic(ops[tag], a0, a2, p), p, verbose=False)
            assert pk['verdict'] == 'STUCK'
            inf = pk['residual_diag']['inf']
            assert inf['rank_le_1'] and inf['lam_unsplit'] == 0, (key, tag, inf)
            lams = [int(r) for r, m in inf['lam_roots']]
            assert len(lams) == len(set(lams)) == order, (key, tag, lams)
            assert quart <= set(lams), (key, tag, lams)
            assert set(lams) - quart == {int(I.to_fp(e, p)) for e in EXTRA[order]}, (key, tag, lams)
            for x in quart:                                            # the found four satisfy the quartic
                assert sum(int(I.to_fp(c, p)) * pow(x, 4 - i, p) for i, c in enumerate(QUARTIC)) % p == 0
    assert sorted(orders) == [4] * 12 + [5] + [6, 6]
    assert len(quart) == 4


def test_aabc_screened_part_not_triangular_fast():
    """Claims (2) and (3) at L = 400.  Wrong answers rejected: any component that peels (TRIANGULAR, or STUCK after a
    peeled step); a flagged STUCK (STUCK-UNCONFIRMED), as on AB's primes; sets A and B disagreeing; a residual order
    or search table different from the record's (orders 4 x 12, 5, 6, 6); an erratum check that fails or finds
    anything; an unsplit lambda edge; found rates that do not satisfy the quartic or extra rates other than 0 / +-mu;
    a component set different from the record's 15."""
    check_all(PRIME_FAST, L_FAST)


@pytest.mark.slow
def test_aabc_screened_part_not_triangular_record_L():
    """As the fast test, at the record's L = 1200, on a second new prime."""
    check_all(PRIME_SLOW, L_SLOW)


def _nmod(lst, p):
    return [flint.nmod_poly([x % p for x in c], p) for c in lst]


def _peel_check(op, p):
    r = I.peel(op, p, verbose=False)
    return r['verdict'], (C.erratum_check(r['residual'], p)['passed'] if r['verdict'] == 'STUCK' else None)


def test_erratum_check_controls():
    """Record IN's erratum-check controls at PRIME_FAST.  Wrong answers rejected: an erratum check that passes
    everything (Chebyshev non-square and Euler would pass) or fails everything (Legendre K and Bessel K0 would
    fail); a search that misses the square-c Chebyshev peel."""
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
    """rho^L/(L+1) (a logarithm, triangular) through operator_sets -> analyse_full at A1's geometry.  Wrong answer
    rejected: a pipeline that returns STUCK for everything."""
    p = PRIME_FAST
    D1, D2 = X.params('A1')[:2]
    rho = D1 / D2
    seq = [I.to_fp(rho ** L / (L + 1), p) for L in range(L_FAST)]
    ops = S.operator_sets(seq, p, 24)
    assert ops is not None
    rA = C.analyse_full(ops['A'], D1, D2, p)
    rB = C.analyse_full(ops['B'], D1, D2, p)
    assert rA['verdict'] == 'TRIANGULAR' and S.component_verdict(rA, rB) == 'TRIANGULAR'


# Fire-test record (python debug/qa/fire_test.py tests/test_paper59_aa_bc_class.py [--plant-in <file>] --plant "O=>N"):
# fire-tested 2026-11-16, every plant FIRED, every plant restored.  -k selections leave out the hash test, so each
# plant is caught by the behavioural test named, not by the sha256 check.
#   01  this file  QUARTIC -51/20 => -51/19                       -k "closed_form or split_iff"  FIRED
#   01b same                                                      -k not_triangular_fast          FIRED
#   02  this file  QUARTIC 81/25 => 81/16                         -k "closed_form or split_iff"  FIRED
#   03  this file  PRIME_FAST => 1073741237 (615 a non-residue)   -k prime_class                  FIRED
#   03b same                                                      -k not_triangular_fast          FIRED; checked by
#       hand that the first failure is the verdict, STUCK-UNCONFIRMED with 'lambda edge has 4 roots outside F_p'
#   04  intclass_score.py  component_verdict return 'STUCK' => 'TRIANGULAR'      -k not_triangular_fast  FIRED
#   05  intclass_in_check.py  passed=not bad => passed=True (check disabled)    -k erratum_check_controls FIRED
#   05b intclass_in_check.py  passed=not bad => passed=False (inverted)         -k not_triangular_fast  FIRED
#   06  this file  EXPECT order-5 component's residual order 5 => 4            -k not_triangular_fast  FIRED
#   06b intclass_lib.py  classify_residual order=order => order=order + 1     -k not_triangular_fast  FIRED
#   07  this file  MU 9/5 => 9/4 (extra rates of the order-6 pair)             -k not_triangular_fast  FIRED
#   08  aabc_exact.py  A1 zA '9/10' => '9/11'                                   -k engine_case          FIRED
#   08b aabc_exact.py  dY kernel's x i_{L+1} term dropped                       -k not_triangular_fast  FIRED
#   09  intclass_lib.py  if sol is None: => if True: (peel never peels)         -k pipeline             FIRED
#   10  aabc_exact.py  docstring one-word edit                                  -k frozen               FIRED
