"""Guard for record FX (docs/predictions/2026-11-17_four_centre_recurrence.md, scored STOP on FX3; FX0-FX2 held):
Paper 63 sec:four_centre, the paragraph "[MEASURED] On one slice of that family the coefficients are exact, and they
carry iterated exponential integrals of depth two that the known identities do not remove."

Object: the four-centre 1s integral (XY|ZW), X at the origin, Y = D1 z, Z = a n, W = b n, x = n.z, F_L its Legendre
coefficient in x; G_L = F_L/pi^2 = sum_c q_{L,c} c reduced to FX's canonical basis (debug/fx_indep.py) mod p.

Claims guarded (each against a named wrong answer):
  (a) geometry: the paper's and the record's (a, b, D1), (zZ, zW), zX = 13/10, zY = 1 equal the engine's cases, and
      both cases satisfy the matched condition zW/zZ = sqrt(a/b) with a/b a ratio of squares (rejects a typo or an
      off-slice exponent); the record's Landau contacts follow from x* = (D1^2 + c^2)/(2 D1 c);
  (b) on the slice the far-density multipole rho2_L(r) is elementary: debug/fx_engine.rho2_split returns only
      e^{-zZ u} and E1(zZ u) terms, and its value agrees with fc_engine's independent arb angular quadrature to
      1e-39 (the record's bar; observed ~1e-106) at radii in every panel; the same comparison against an off-slice
      exponent (FC's G2, zW = 7/10) misses by > 1e-3, and fx_engine refuses off-slice parameters;
  (c) the E1 coefficient R_L(r) is proportional to (b - a)^2: symbolically (a, b, zZ, r free), L = 0..2, R_L vanishes
      at b = a and R_L/((b - a)^2 (r^2 - ab)^2) does not; L = 0 equals -(1/(2ra)) zZ^2 gamma^2/8; numerically the
      split at a = b carries no E1 term;
  (d) an exact recomputation through FX's engine (debug/fx_exact.py + fx_indep.canonical, imported, not edited) for
      L <= 40 at M1 and M2r on BOTH primes equals the prefix of the saved sequences debug/data/fx_seq_<case>_<p>.json
      key by key, with a depth-2 key nonzero at every 2 <= L <= 40 (this file's classifier agrees with
      fx_indep.depth2_keys); the saved files' sha256 equal the hashes frozen in debug/fxk_run.py (SEQ_SHA);
  (e) from the saved sequences (FX1/FX2 as the paper states them): 67 keys at M1 and 76 at M2r, the same set on both
      primes, every key nonzero (on some prime) at every L = 2..1500, and some depth-2 key nonzero at every L; the
      measured depth-2 inventory (16 keys per case, each a single-E1 base) is pinned too;
  (f) the search budget the paper prints (order <= 32, degree <= 35, window l <= 1300, 2..1500) equals the frozen
      run files' settings and trade table, in which no run found an operator.
  The absence of a recurrence is NOT recomputed (record FX3, ~7000 s per run); (f) only reads the scored files.
  Slow tier: recomputation to L <= 150 on one prime.
"""
from __future__ import annotations

import ast
import hashlib
import json
import re
import sys
from fractions import Fraction
from math import isqrt
from pathlib import Path
from typing import Dict, List, Tuple

import pytest

flint = pytest.importorskip("flint")
mp = pytest.importorskip("mpmath")

ROOT = Path(__file__).resolve().parent.parent
DEBUG = ROOT / "debug"
SEQ_DIR = DEBUG / "data"
PAPER = ROOT / "papers" / "group3_foundations" / "paper_63_density_graph_ladder.tex"
RECORD = ROOT / "docs" / "predictions" / "2026-11-17_four_centre_recurrence.md"
sys.path.insert(0, str(DEBUG))

# ----------------------------------------------------------------------------------------- record constants
PRIMES = (1073741441, 1073741381)
SEQ_SHA = {('M1', 1073741441): 'edcf845b3186660750e59e595724967ddb97e6903906ecbedf8d9d2ec4cbc3e1',
           ('M1', 1073741381): 'faf2f135acc4f6b7166af44a24734360d36f39dd0448cb0fa58743b8acc0526e',
           ('M2r', 1073741441): '55b8dfdf882c2a06ff765eb395d95244491f5129063172fa56fe31caa00e8576',
           ('M2r', 1073741381): 'e6374838388e38404cdcf9af7a90a19237a4d50de39d418b06e7f557239d489c'}
ZX, ZY = Fraction(13, 10), Fraction(1)
GEOM = {'M1': dict(a=Fraction(1), b=Fraction(4), D1=Fraction(6, 5), zZ=Fraction(6, 5), zW=Fraction(3, 5)),
        'M2r': dict(a=Fraction(4, 9), b=Fraction(1), D1=Fraction(6, 5), zZ=Fraction(3, 2), zW=Fraction(1))}
CONTACTS = {'M1': ('61/60', '109/60'), 'M2r': ('829/540', '61/60')}     # record table: x_Z*, x_W*
N_KEYS = {'M1': 67, 'M2r': 76}
N_DEPTH2 = 16                    # measured from the saved sequences (not stated in the paper)
LMAX, WINDOW, RMAX, DEGCAP = 1500, 1300, 32, 35
L_RECOMP = 40
L_SLOW = 150


def _case_dict(case: str) -> Dict[str, str]:
    g = GEOM[case]
    return dict(a=str(g['a']), b=str(g['b']), D1=str(g['D1']), zX=str(ZX), zY=str(ZY), zZ=str(g['zZ']),
                zW=str(g['zW']))


def is_depth2(key: Tuple) -> bool:
    """A key carries a named base B (an integral of E1 factors against a rational kernel)."""
    return any(a[0] == 'B' for a in key[1])


# ----------------------------------------------------------------------------------------- loaders
_NS = {'Fraction': Fraction, '__builtins__': {}}
_CACHE: Dict[Tuple[str, int], Dict[str, list]] = {}


def _seq_file(case: str, p: int) -> Path:
    return SEQ_DIR / f'fx_seq_{case}_{p}.json'


def load_seqs(case: str, p: int) -> Dict[str, list]:
    if (case, p) not in _CACHE:
        fn = _seq_file(case, p)
        if not fn.exists():
            pytest.fail(f'missing saved sequence file {fn}')
        _CACHE[(case, p)] = json.loads(fn.read_text())
    return _CACHE[(case, p)]


def parse_key(k: str) -> Tuple:
    return eval(k, _NS)          # noqa: S307  (repo data, Fraction-only namespace)


def _paragraph() -> str:
    t = PAPER.read_text(encoding='utf-8')
    i = t.index('On one slice of that family the coefficients are')
    j = t.index('registered before computing).', i)
    return re.sub(r'\s+', '', t[i:j])


def _tex(q: Fraction) -> str:
    if q.denominator == 1:
        return str(q.numerator)
    if q.numerator < 10 and q.denominator < 10:
        return f'\\tfrac{q.numerator}{q.denominator}'
    return f'\\tfrac{{{q.numerator}}}{{{q.denominator}}}'


# ----------------------------------------------------------------------------------------- (a) geometry
@pytest.mark.parametrize('case', ['M1', 'M2r'])
def test_geometry_equals_engine_cases(case):
    import fx_exact as FX
    p = FX.params(case)
    g = GEOM[case]
    for k in ('a', 'b', 'D1', 'zZ', 'zW'):
        assert p[k] == g[k], (case, k, p[k], g[k])
    assert p['zX'] == ZX and p['zY'] == ZY


@pytest.mark.parametrize('case', ['M1', 'M2r'])
def test_matched_condition(case):
    g = GEOM[case]
    k2 = g['a'] / g['b']
    n, m = isqrt(k2.numerator), isqrt(k2.denominator)
    assert n * n == k2.numerator and m * m == k2.denominator, 'a/b not a ratio of squares'
    assert g['zW'] / g['zZ'] == Fraction(n, m), 'off the matched slice'
    assert g['a'] < g['b']


def test_paper_states_geometry():
    par = _paragraph()
    m1, m2 = GEOM['M1'], GEOM['M2r']
    need = [
        '(a,b,D_1)=(' + ','.join(_tex(m1[k]) for k in ('a', 'b', 'D1')) + ')',
        '(\\zeta_Z,\\zeta_W)=(' + ','.join(_tex(m1[k]) for k in ('zZ', 'zW')) + ')',
        '(' + ','.join(_tex(m2[k]) for k in ('a', 'b', 'D1')) + ')$with$(' + ','.join(_tex(m2[k]) for k in ('zZ', 'zW')) + ')$',
        '\\zeta_X=' + _tex(ZX), '\\zeta_Y=' + _tex(ZY),
        '\\zeta_W/\\zeta_Z=\\sqrt{a/b}',
        f'${N_KEYS["M1"]}$and${N_KEYS["M2r"]}$',
        'from$2$to$' + str(LMAX) + '$',
        f'orderatmost${RMAX}$', f'degreeatmost${DEGCAP}$', f'\\ell\\le{WINDOW}',
        'proportionalto$(b-a)^2$',
        'exactlymodulotwoprimes',
    ]
    for s in need:
        assert s in par, s


def test_record_states_geometry_and_contacts():
    t = RECORD.read_text(encoding='utf-8')
    for case in ('M1', 'M2r'):
        g = GEOM[case]
        row = (f"| {case} | ({g['a']}, {g['b']}, {g['D1']}) | ({g['zZ']}, {g['zW']}) |")
        assert row in t, row
        D1 = g['D1']
        xs = tuple(str((D1 ** 2 + c ** 2) / (2 * D1 * c)) for c in (g['a'], g['b']))
        assert xs == CONTACTS[case], (case, xs)
        assert f'x_Z* = {xs[0]}' in t and f'x_W* = {xs[1]}' in t
    assert f'primes {PRIMES[0]} and {PRIMES[1]}' in t


# ----------------------------------------------------------------------------------------- (b) rho2 elementary
def _fc_rho2(d: Dict[str, str], r: Fraction, Lmax: int):
    import fc_engine as FC
    FC.ctx.prec = 400            # rho2_point does not set it (design memo, control K1 note)
    return FC.rho2_point(dict(d), str(r), Lmax, 'A')


RADII = (Fraction(1, 3), Fraction(5, 7), Fraction(5, 2), Fraction(6, 1))


@pytest.mark.parametrize('case', ['M1', 'M2r'])
def test_rho2_elementary_matches_quadrature(case):
    import fx_engine as FE
    mp.mp.dps = 120
    d = _case_dict(case)
    worst = mp.mpf(0)
    for r in RADII:
        num = _fc_rho2(d, r, 6)
        for L in (0, 2, 6):
            sp_ = FE.rho2_split(d, r, L)
            assert {k[0] for k in sp_} <= {'e', 'E1'}, sp_.keys()
            assert any(k[0] == 'E1' for k in sp_), 'E1 part missing'
            ex = FE.rho2_value(sp_, d, mp)
            nv = mp.mpf(num[L].mid().str(110, radius=False))
            worst = max(worst, abs(ex - nv) / abs(nv))
    assert worst < mp.mpf('1e-39'), mp.nstr(worst, 3)


def test_rho2_offslice_control():
    """The comparison above discriminates: at FC's off-slice G2 exponent (zW = 7/10, same a, b, D1 as M1) the
    quadrature misses the matched closed form by > 1e-3, and the engine refuses the off-slice parameters."""
    import fx_engine as FE
    mp.mp.dps = 60
    on = _case_dict('M1')
    off = dict(on, zW='7/10')
    r = Fraction(5, 2)
    num = _fc_rho2(off, r, 2)
    ex = FE.rho2_value(FE.rho2_split(on, r, 2), on, mp)
    nv = mp.mpf(num[2].mid().str(50, radius=False))
    assert abs(ex - nv) / abs(nv) > mp.mpf('1e-3')
    with pytest.raises(ValueError):
        FE.params(off)


# ----------------------------------------------------------------------------------------- (c) E1 coefficient ~ (b-a)^2
def _residue(L: int):
    import sympy as sp
    a, b, z, r, u = sp.symbols('a b z r u', positive=True)
    g = (b - a) * (r ** 2 - a * b) / b
    s2 = (u ** 2 + 2 * g + g ** 2 / u ** 2) / 4
    mu = (r ** 2 + a ** 2 - s2) / (2 * r * a)
    h = sp.expand(sp.legendre(L, mu) * (u - g ** 2 / u ** 3) / 4)
    ser = sp.series(sp.exp(-z * u), u, 0, 2 * L + 4).removeO()
    res = sp.expand(h * ser).coeff(u, -1)
    return sp.Rational(2 * L + 1) / (2 * r * a) * res, (a, b, z, r, g)


@pytest.mark.parametrize('L', [0, 1, 2])
def test_e1_coefficient_is_order_two_in_b_minus_a(L):
    import sympy as sp
    R, (a, b, z, r, g) = _residue(L)
    assert sp.simplify(R.subs(b, a)) == 0
    q = sp.cancel(R / ((b - a) ** 2 * (r ** 2 - a * b) ** 2))
    assert sp.simplify(q.subs(b, a)) != 0, 'vanishes to higher order than (b - a)^2'
    num, den = sp.fraction(sp.together(q))
    assert sp.simplify(den.subs(b, a)) != 0, 'quotient singular at b = a'
    if L == 0:
        assert sp.simplify(R - (-(z ** 2) * g ** 2 / (8 * 2 * r * a))) == 0


def test_rho2_split_engine_agrees_with_symbolic_residue_and_a_eq_b():
    """The engine's E1 coefficient equals the symbolic residue (M1, L = 2, r = 5/2), and at a = b it is absent."""
    import fx_engine as FE
    import sympy as sp
    d = _case_dict('M1')
    g = GEOM['M1']
    r0 = Fraction(5, 2)
    R, (a, b, z, r, _) = _residue(2)
    val = R.subs({a: sp.Rational(1), b: sp.Rational(4), z: sp.Rational(6, 5), r: sp.Rational(5, 2)})
    sp_ = FE.rho2_split(d, r0, 2)
    kap = Fraction(1, 2)
    um = abs(r0 - g['a']) + kap * abs(r0 - g['b'])
    assert sp.Rational(str(sp_[('E1', um)])) == sp.nsimplify(val)
    eq = dict(a='1', b='1', D1='6/5', zX='13/10', zY='1', zZ='6/5', zW='6/5')
    assert not any(k[0] == 'E1' for k in FE.rho2_split(eq, r0, 2))


# ----------------------------------------------------------------------------------------- (d) recomputation
def _sha(fn: Path) -> str:
    return hashlib.sha256(fn.read_bytes()).hexdigest()


@pytest.mark.parametrize('case,p', [(c, p) for c in ('M1', 'M2r') for p in PRIMES])
def test_saved_sequence_sha(case, p):
    assert _sha(_seq_file(case, p)) == SEQ_SHA[(case, p)]


def test_sha_matches_fxk_run_if_present():
    fn = DEBUG / 'fxk_run.py'
    if not fn.exists():
        pytest.skip('debug/fxk_run.py absent')
    tree = ast.parse(fn.read_text(encoding='utf-8'))
    found = None
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(getattr(t, 'id', None) == 'SEQ_SHA' for t in node.targets):
            found = ast.literal_eval(node.value)
    if found is None:
        pytest.skip('SEQ_SHA not defined in fxk_run.py')
    assert found == SEQ_SHA


def _recompute(case: str, p: int, Lmax: int) -> List[dict]:
    import fx_exact as FX
    import fx_indep as C
    import na_exact as NE
    K = NE.Field(p)
    Fn, stable = C.natural_factors(case)
    assert stable, 'natural factor sets L-dependent'
    out = []
    for L, q in FX.Engine(case, K).run(Lmax):
        c = C.canonical(q, Fn, K)
        assert sorted(map(repr, C.depth2_keys(c))) == sorted(repr(k) for k in c if is_depth2(k))
        out.append(c)
    return out


def _check_prefix(case: str, p: int, Lmax: int) -> None:
    seqs = load_seqs(case, p)
    rec = _recompute(case, p, Lmax)
    for L, c in enumerate(rec):
        cr = {repr(k): int(v) for k, v in c.items()}
        assert set(cr) <= set(seqs), (case, p, L, sorted(set(cr) - set(seqs))[:3])
        for k, s in seqs.items():
            assert s[L] == cr.get(k, 0), (case, p, L, k[:80])
        if L >= 2:
            assert any(v for k, v in c.items() if is_depth2(k)), (case, p, L, 'no depth-2 key')


@pytest.mark.parametrize('case,p', [(c, p) for c in ('M1', 'M2r') for p in PRIMES])
def test_recompute_equals_saved_prefix(case, p):
    _check_prefix(case, p, L_RECOMP)


@pytest.mark.slow
@pytest.mark.parametrize('case', ['M1', 'M2r'])
def test_recompute_equals_saved_prefix_slow(case):
    _check_prefix(case, PRIMES[0], L_SLOW)


# ----------------------------------------------------------------------------------------- (e) FX1 / FX2
@pytest.mark.parametrize('case', ['M1', 'M2r'])
def test_key_count_and_same_set_on_both_primes(case):
    s1, s2 = load_seqs(case, PRIMES[0]), load_seqs(case, PRIMES[1])
    assert set(s1) == set(s2)
    assert len(s1) == N_KEYS[case]
    assert all(len(v) == LMAX + 1 for v in list(s1.values()) + list(s2.values()))


@pytest.mark.parametrize('case', ['M1', 'M2r'])
def test_key_set_constant_for_L_2_to_1500(case):
    s1, s2 = load_seqs(case, PRIMES[0]), load_seqs(case, PRIMES[1])
    for k in s1:
        a, b = s1[k], s2[k]
        miss = [L for L in range(2, LMAX + 1) if a[L] % PRIMES[0] == 0 and b[L] % PRIMES[1] == 0]
        assert not miss, (k[:80], miss[:5])


@pytest.mark.parametrize('case', ['M1', 'M2r'])
def test_depth2_key_nonzero_at_every_L(case):
    s1, s2 = load_seqs(case, PRIMES[0]), load_seqs(case, PRIMES[1])
    d2 = [k for k in s1 if is_depth2(parse_key(k))]
    assert len(d2) == N_DEPTH2
    for k in d2:                 # measured inventory: every canonical base has a single E1 factor
        bs = [a for a in parse_key(k)[1] if a[0] == 'B']
        assert len(bs) == 1 and len(bs[0][2]) == 1, k
    miss = [L for L in range(2, LMAX + 1) if not any(s1[k][L] or s2[k][L] for k in d2)]
    assert not miss, miss[:5]


# ----------------------------------------------------------------------------------------- (f) search budget
@pytest.mark.parametrize('case,p', [(c, p) for c in ('M1', 'M2r') for p in PRIMES])
def test_search_budget_and_no_operator_in_scored_runs(case, p):
    r = json.loads((SEQ_DIR / f'fx_run_{case}_{p}.json').read_text())
    assert r['instrument_error'] is None
    assert (r['LMAX'], r['LMAX'] - r['HOLD'], r['RMAX']) == (LMAX, WINDOW, RMAX)
    assert r['n_keys'] == N_KEYS[case] and r['natural_stable'] and r['no_depth2_L'] == []
    assert r['outside_alphabet'] == []
    ops = r['op']['ops']
    assert sorted(o['seed'] for o in ops) == [7, 2026]
    for o in ops:
        assert o['found'] is False
        assert sorted(map(int, o['trade'])) == list(range(1, RMAX + 1))
        assert o['trade'][str(RMAX)] == [None, DEGCAP]
