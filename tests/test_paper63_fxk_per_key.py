"""Guard for record FXK (docs/predictions/2026-11-17_fxk_per_key_recurrence.md, scored GO at M1 and M2r):
Paper 63 sec:four_centre, the sentences beginning "Searched constant by constant within the same budget".

Object: FX's matched-slice four-centre integral; G_L = F_L/pi^2 = sum_k c_k(L) key_k on FX's canonical keys; for each
key the minimal bond-angle operator of f_k(x) = sum_L c_k(L) P_L(x), searched on the window L <= 1300 (R <= 32), mod
the primes 1073741441 and 1073741381.  The input is FX's committed sequences debug/data/fx_seq_<case>_<p>.json.

Claims guarded (each against a named wrong answer):
  (a) inputs: the four sequence files have the sha256 frozen in record FX / fxk_run.SEQ_SHA (rejects an altered
      sequence); the contacts follow from the geometry, x* = (D1^2 + c^2)/(2 D1 c), and equal fx_run.landau (rejects
      a mistyped contact such as 62/60);
  (b) RECOMPUTED, not read back: for each contact, the named pure-exponential key's operator is rebuilt from the
      sequence (nullspace at the search shape (16, 52) on the window, nothing at (16, 51), nothing of order <= 15 at
      the window's largest degree, gcrd of the annihilators, content-reduced) and is pinned at (order, degree) =
      (15, 77) (rejects 15 -> 14, 77 -> 76); the operator, rewritten as sum q_i(x) T^i, annihilates the sequence on
      EVERY row the rewritten form reaches, including rows past the window (an independent numpy implementation of
      x and T on Legendre coefficients, not na_modguess); its leading coefficient vanishes at its contact and not at
      the other one (rejects the wrong key for a contact, and a contact carried twice); the point is regular singular
      and the indicial polynomial factors completely over Q with roots {0, ..., 13} and 7/2 (rejects 3/2 or 5/2, and
      any extra half-integer root); all of this identically on both primes;
  (c) every OTHER found key's operator, rebuilt at its search shape, has the pinned shape and a leading coefficient
      that does not vanish at either contact: each contact is a singular point of exactly one found operator;
  (d) from the scored run files (read, not recomputed; the absence of the other operators costs ~2000 s per run and
      is rechecked only in the slow tier): 12 of 67 keys found at M1 and 15 of 76 at M2r, none of depth two ('B'),
      the same keys with the same (order, degree) on both primes, the shape histogram the record prints;
  (e) the paper paragraph states 12, 67, 15, 76, order 15, 7/2, "the only half-integer one", "identically on both
      primes", and the scope disclaimer "it does not show that $F_\\ell$ itself carries"; the record's Results table
      names the same carriers.
  Slow tier: the full per-key search (fxk_run.per_key, ~2000 s per case and prime) reproduces the found set, shapes
  and fingerprints; and fx_run.minimal_windowed on each carrier key reproduces (15, 77) with its own R-bisection.

Paths can be redirected (fire tests on temp copies): GEOVAC_FXK_PAPER, GEOVAC_FXK_RECORD, GEOVAC_FXK_DATA.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import sys
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Tuple

import pytest

flint = pytest.importorskip("flint")
np = pytest.importorskip("numpy")

ROOT = Path(__file__).resolve().parent.parent
DEBUG = ROOT / "debug"
DATA = Path(os.environ.get("GEOVAC_FXK_DATA", DEBUG / "data"))
PAPER = Path(os.environ.get("GEOVAC_FXK_PAPER",
                            ROOT / "papers" / "group3_foundations" / "paper_63_density_graph_ladder.tex"))
RECORD = Path(os.environ.get("GEOVAC_FXK_RECORD",
                             ROOT / "docs" / "predictions" / "2026-11-17_fxk_per_key_recurrence.md"))
sys.path.insert(0, str(DEBUG))

# ----------------------------------------------------------------------------------------- record constants
PRIMES = (1073741441, 1073741381)
SEQ_SHA = {('M1', 1073741441): 'edcf845b3186660750e59e595724967ddb97e6903906ecbedf8d9d2ec4cbc3e1',
           ('M1', 1073741381): 'faf2f135acc4f6b7166af44a24734360d36f39dd0448cb0fa58743b8acc0526e',
           ('M2r', 1073741441): '55b8dfdf882c2a06ff765eb395d95244491f5129063172fa56fe31caa00e8576',
           ('M2r', 1073741381): 'e6374838388e38404cdcf9af7a90a19237a4d50de39d418b06e7f557239d489c'}
GEOM = {'M1': dict(a=Fraction(1), b=Fraction(4), D1=Fraction(6, 5)),
        'M2r': dict(a=Fraction(4, 9), b=Fraction(1), D1=Fraction(6, 5))}
CONTACTS = {'M1': {'xZ': Fraction(61, 60), 'xW': Fraction(109, 60)},
            'M2r': {'xZ': Fraction(829, 540), 'xW': Fraction(61, 60)}}
# carrier of each contact: a pure-exponential key (prefactor, no E1 factor)
CARRIER = {'M1': {'xZ': Fraction(-84, 25), 'xW': Fraction(-129, 25)},
           'M2r': {'xZ': Fraction(-476, 225), 'xW': Fraction(-359, 150)}}
CARRIER_SEARCH = (16, 52)            # (R, Q) at which the annihilators first appear on the window
CARRIER_SHAPE = (15, 77)             # (order, degree) of their gcrd
HALF_ROOTS = [Fraction(7, 2)]        # the only half-integer indicial root at each contact
INT_ROOTS = list(range(14))          # the remaining indicial roots there (record Reading: 0-13)
N_KEYS = {'M1': 67, 'M2r': 76}
N_FOUND = {'M1': 12, 'M2r': 15}
SHAPES = {'M1': {(6, 21): 2, (10, 37): 1, (10, 44): 6, (11, 44): 1, (15, 77): 2},
          'M2r': {(6, 21): 3, (10, 44): 8, (11, 44): 2, (15, 77): 2}}
LMAX, WINDOW, RMAX, MARGIN = 1500, 1300, 32, 30


def carrier_key(case: str, m: str) -> str:
    return repr((CARRIER[case][m], ()))


# ----------------------------------------------------------------------------------------- loaders
_NS = {'Fraction': Fraction, '__builtins__': {}}
_SEQ: Dict[Tuple[str, int], Dict[str, list]] = {}


def _seq_file(case: str, p: int) -> Path:
    return DATA / f'fx_seq_{case}_{p}.json'


def _sha(fn: Path) -> str:
    return hashlib.sha256(fn.read_bytes()).hexdigest()


def load_seqs(case: str, p: int) -> Dict[str, list]:
    if (case, p) not in _SEQ:
        fn = _seq_file(case, p)
        if not fn.exists():
            pytest.fail(f'missing saved sequence file {fn}')
        assert _sha(fn) == SEQ_SHA[(case, p)], f'{fn.name}: sha256 differs from the frozen value'
        d = json.loads(fn.read_text())
        _SEQ[(case, p)] = {k: [int(v) % p for v in s] for k, s in d.items()}
    return _SEQ[(case, p)]


def load_run(case: str, p: int) -> dict:
    fn = DATA / f'fxk_run_{case}_{p}.json'
    if not fn.exists():
        pytest.skip(f'{fn.name} absent')
    return json.loads(fn.read_text())


def parse_key(k: str) -> Tuple:
    return eval(k, _NS)          # noqa: S307  (repo data, Fraction-only namespace)


def is_depth2(key: Tuple) -> bool:
    return any(a[0] == 'B' for a in key[1])


# ----------------------------------------------------------------------------------------- operator machinery
def _operator(seq: List[int], p: int, R: int, Q: int):
    """Annihilators of order <= R, degree <= Q on the window (FX's ansatz, na_modguess), gcrd, content-reduced."""
    import intclass_lib as I
    import na_modanalyze as A
    import na_modguess as G
    Bg = G.Basis([flint.nmod(v, p) for v in seq[:WINDOW + 1]], p)
    nul, keys, info = G.nullity(Bg, R, Q, MARGIN)
    assert nul and nul >= 1, f'no annihilator at (R, Q) = ({R}, {Q})'
    N, nrow = info
    ops = [A.standard_form(keys, [N[i, c] for i in range(len(keys))], p) for c in range(nul)]
    g = ops[0]
    for o in ops[1:]:
        g = I.gcrd(g, o, p)
    return I.content_reduce(g), Bg, nrow


def _shape(c) -> Tuple[int, int]:
    return len(c) - 1, max(x.degree() for x in c)


class _Leg:
    """x and T = (1 - x^2) d/dx acting on Legendre coefficient sequences mod p, numpy int64 (written here, not
    taken from na_modguess, so the annihilation check does not share code with the search)."""

    def __init__(self, p: int, n: int):
        self.p = p
        l = np.arange(n, dtype=object)
        inv = lambda a: np.array([pow(int(v) % p, p - 2, p) for v in a], dtype=np.int64)  # noqa: E731
        self.xa = np.array([(int(v) * int(w)) % p for v, w in zip(l, inv(2 * l - 1))], dtype=np.int64)  # l/(2l-1)
        self.xb = np.array([(int(v + 1) * int(w)) % p for v, w in zip(l, inv(2 * l + 3))], dtype=np.int64)
        self.ta = np.array([(int((v - 1) * v) * int(w)) % p for v, w in zip(l, inv(2 * l - 1))], dtype=np.int64)
        self.tb = np.array([(int((v + 1) * (v + 2)) * int(w)) % p for v, w in zip(l, inv(2 * l + 3))],
                           dtype=np.int64)

    def X(self, s):
        n, p = len(s) - 1, self.p
        out = self.xb[:n] * s[1:n + 1] % p
        out[1:] = (out[1:] + self.xa[1:n] * s[:n - 1]) % p
        return out

    def T(self, s):
        n, p = len(s) - 1, self.p
        out = self.tb[:n] * s[1:n + 1] % p
        out[1:] = (out[1:] - self.ta[1:n] * s[:n - 1]) % p
        return out


def _t_form(c, p: int):
    """Rewrite sum_k c_k(x) D^k, left-multiplied by (1 - x^2)^R, as sum_k q_k(x) T^k; returns [q_k as int lists]."""
    import na_modanalyze as A
    R = len(c) - 1
    w = flint.nmod_poly([1, 0, -1], p)
    Tp = A.T_powers(R, p)
    Lp = [w ** R * ck for ck in c]
    q = [None] * (R + 1)
    for k in range(R, -1, -1):
        qk, r = divmod(Lp[k], w ** k)
        assert r.degree() < 0, 'T-form conversion: remainder'
        q[k] = qk
        for m, b in Tp[k].items():
            Lp[m] = Lp[m] - qk * b
    assert all(x.degree() < 0 for x in Lp), 'T-form conversion incomplete'
    return [[int(a) for a in qk.coeffs()] for qk in q]


def _apply_rows(c, seq: List[int], p: int):
    """(L f) on every Legendre row the T-form reaches; returns the residual array."""
    q = _t_form(c, p)
    lg = _Leg(p, len(seq) + 2)
    v = np.array(seq, dtype=np.int64)
    Tv = [v]
    for _ in range(len(q) - 1):
        Tv.append(lg.T(Tv[-1]))
    parts = []
    for i, qi in enumerate(q):
        if not qi:
            continue
        r = (qi[-1] * Tv[i]) % p
        for a in reversed(qi[:-1]):              # Horner in x
            r = lg.X(r)
            r = (r + a * Tv[i][:len(r)]) % p
        parts.append(r)
    n = min(len(r) for r in parts)
    tot = np.zeros(n, dtype=np.int64)
    for r in parts:
        tot = (tot + r[:n]) % p
    return tot


def _lc_at(c, x: Fraction, p: int) -> int:
    import intclass_lib as I
    return int(c[-1](flint.nmod(I.to_fp(x, p), p)))


def _indicial_roots(c, x: Fraction, p: int):
    import intclass_lib as I
    import na_modanalyze as A
    loc = A.local(c, flint.nmod(I.to_fp(x, p), p), p, {})
    assert loc['regular'], 'contact is not a regular singular point'
    ind = loc['indicial']
    _, facs = ind.factor()
    roots = []
    for f, mlt in facs:
        assert f.degree() == 1, f'indicial factor of degree {f.degree()} (irrational exponents)'
        r = (-int(f[0]) * pow(int(f[1]), p - 2, p)) % p
        fr = I.ratrec(r, p, bound=1000)
        assert fr is not None
        roots += [Fraction(fr)] * mlt
    assert len(roots) == ind.degree()
    return sorted(roots)


# ----------------------------------------------------------------------------------------- (a) inputs
@pytest.mark.parametrize('case,p', [(c, p) for c in ('M1', 'M2r') for p in PRIMES])
def test_saved_sequence_sha(case, p):
    assert _sha(_seq_file(case, p)) == SEQ_SHA[(case, p)]


def test_seq_sha_equals_fxk_run_pin():
    import ast
    tree = ast.parse((DEBUG / 'fxk_run.py').read_text(encoding='utf-8'))
    found = [ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign)
             and any(getattr(t, 'id', None) == 'SEQ_SHA' for t in n.targets)]
    assert found and found[0] == SEQ_SHA


@pytest.mark.parametrize('case', ['M1', 'M2r'])
def test_contacts_from_geometry(case):
    import fx_run as R
    g = GEOM[case]
    D1 = g['D1']
    xs = {m: (D1 ** 2 + c ** 2) / (2 * D1 * c) for m, c in (('xZ', g['a']), ('xW', g['b']))}
    assert xs == CONTACTS[case]
    lm = R.landau(case)
    assert {m: Fraction(lm[m]) for m in ('xZ', 'xW')} == CONTACTS[case]


# ----------------------------------------------------------------------------------------- (b) carriers, recomputed
_OPS: Dict[Tuple[str, int, str], object] = {}


def _carrier_op(case: str, p: int, m: str):
    if (case, p, m) not in _OPS:
        seq = load_seqs(case, p)[carrier_key(case, m)]
        c, Bg, nrow = _operator(seq, p, *CARRIER_SEARCH)
        _OPS[(case, p, m)] = (c, Bg, nrow)
    return _OPS[(case, p, m)]


@pytest.mark.parametrize('case,p,m', [(c, p, m) for c in ('M1', 'M2r') for p in PRIMES for m in ('xZ', 'xW')])
def test_carrier_operator_recomputed(case, p, m):
    import na_modguess as G
    seq = load_seqs(case, p)[carrier_key(case, m)]
    c, Bg, nrow = _carrier_op(case, p, m)
    R0, Q0 = CARRIER_SEARCH
    n, _, _ = G.nullity(Bg, R0, Q0 - 1, MARGIN)
    assert not n, 'annihilator below the search degree'
    # largest degree the window admits at order R0 - 1 (na_modguess.min_degree_bisect's Qtop)
    Qtop = 0
    while R0 * (Qtop + 2) + MARGIN <= WINDOW + 1 - (R0 - 1) - (Qtop + 1):
        Qtop += 1
    n, _, _ = G.nullity(Bg, R0 - 1, Qtop, MARGIN)
    assert not n, f'annihilator of order <= {R0 - 1} at degree <= {Qtop}'
    assert _shape(c) == CARRIER_SHAPE, _shape(c)
    # the gcrd itself annihilates the sequence on every row it reaches, past the window
    res = _apply_rows(c, seq, p)
    assert len(res) - 1 > WINDOW + 50, f'only {len(res)} rows checked'
    assert not res.any(), f'nonzero residual at rows {np.nonzero(res)[0][:5]}'
    # a perturbed sequence is NOT annihilated (the check discriminates)
    bad = list(seq)
    bad[WINDOW + 20] = (bad[WINDOW + 20] + 1) % p
    assert _apply_rows(c, bad, p).any()


@pytest.mark.parametrize('case,p,m', [(c, p, m) for c in ('M1', 'M2r') for p in PRIMES for m in ('xZ', 'xW')])
def test_carrier_singular_at_its_contact_with_seven_halves(case, p, m):
    c, _, _ = _carrier_op(case, p, m)
    other = 'xW' if m == 'xZ' else 'xZ'
    assert _lc_at(c, CONTACTS[case][m], p) == 0, 'leading coefficient does not vanish at the contact'
    assert _lc_at(c, CONTACTS[case][other], p) != 0, 'carrier also singular at the other contact'
    roots = _indicial_roots(c, CONTACTS[case][m], p)
    half = [r for r in roots if r.denominator == 2]
    assert half == HALF_ROOTS, half
    assert sorted(r for r in roots if r.denominator == 1) == INT_ROOTS
    assert all(r.denominator in (1, 2) for r in roots)


@pytest.mark.parametrize('case,m', [(c, m) for c in ('M1', 'M2r') for m in ('xZ', 'xW')])
def test_carrier_identical_on_both_primes(case, m):
    """Same shape, same contact roots on both primes; distinct carriers at the two contacts of a case."""
    shapes = {_shape(_carrier_op(case, p, m)[0]) for p in PRIMES}
    assert shapes == {CARRIER_SHAPE}
    roots = {tuple(_indicial_roots(_carrier_op(case, p, m)[0], CONTACTS[case][m], p)) for p in PRIMES}
    assert len(roots) == 1
    assert CARRIER[case]['xZ'] != CARRIER[case]['xW']


@pytest.mark.parametrize('case,p', [(c, p) for c in ('M1', 'M2r') for p in PRIMES])
def test_carrier_fingerprint_matches_scored_run(case, p):
    """Cross-check of the recomputation against the scored run (the scored file is not the source of (b))."""
    import fxk_run as K
    r = load_run(case, p)
    for m in ('xZ', 'xW'):
        pk = r['per_key'][carrier_key(case, m)]
        assert pk['found'] and (pk['order'], pk['degree']) == CARRIER_SHAPE
        assert K.fingerprint(_carrier_op(case, p, m)[0]) == pk['fp']
        assert r['contacts'][pk['fp']][m]['present']


# ----------------------------------------------------------------------------------------- (c) other found keys
@pytest.mark.parametrize('case', ['M1', 'M2r'])
def test_other_found_operators_not_singular_at_contacts(case):
    p = PRIMES[0]
    r = load_run(case, p)
    seqs = load_seqs(case, p)
    carriers = {carrier_key(case, m) for m in ('xZ', 'xW')}
    hist: Dict[Tuple[int, int], int] = {}
    for k, pk in r['per_key'].items():
        if not pk.get('found'):
            continue
        if k in carriers:
            hist[CARRIER_SHAPE] = hist.get(CARRIER_SHAPE, 0) + 1
            continue
        c, _, _ = _operator(seqs[k], p, *pk['search'])
        sh = _shape(c)
        hist[sh] = hist.get(sh, 0) + 1
        for m in ('xZ', 'xW'):
            assert _lc_at(c, CONTACTS[case][m], p) != 0, (k, m, 'a second found operator at the contact')
    assert hist == SHAPES[case], hist


# ----------------------------------------------------------------------------------------- (d) scored runs
@pytest.mark.parametrize('case', ['M1', 'M2r'])
def test_found_set_from_scored_runs(case):
    r1, r2 = (load_run(case, p) for p in PRIMES)
    for r in (r1, r2):
        assert r['instrument_error'] is None
        assert (r['LMAX'], r['LMAX'] - r['HOLD'], r['RMAX']) == (LMAX, WINDOW, RMAX)
        assert {m: Fraction(r['landau'][m]) for m in ('xZ', 'xW')} == CONTACTS[case]
        assert sum(v['in_scope'] for v in r['per_key'].values()) == N_KEYS[case]
    assert set(r1['per_key']) == set(load_seqs(case, PRIMES[0]))
    f1 = {k: (v['order'], v['degree']) for k, v in r1['per_key'].items() if v.get('found')}
    f2 = {k: (v['order'], v['degree']) for k, v in r2['per_key'].items() if v.get('found')}
    assert f1 == f2, 'found set or shapes differ between primes'
    assert len(f1) == N_FOUND[case]
    assert not [k for k in f1 if is_depth2(parse_key(k))], 'a depth-2 key has an operator'
    hist: Dict[Tuple[int, int], int] = {}
    for sh in f1.values():
        hist[sh] = hist.get(sh, 0) + 1
    assert hist == SHAPES[case]
    for r in (r1, r2):
        assert all(v['verified'] and v['min_held_out'] >= LMAX - WINDOW
                   for v in r['per_key'].values() if v.get('found'))
        for m in ('xZ', 'xW'):
            pres = [fp for fp, at in r['contacts'].items() if at[m]['present']]
            assert len(pres) == 1, (m, pres)
            assert pres[0] == r['per_key'][carrier_key(case, m)]['fp']


# ----------------------------------------------------------------------------------------- (e) paper + record
def _paragraph() -> str:
    t = PAPER.read_text(encoding='utf-8')
    i = t.index('Searched constant by constant')
    j = t.index('registered before computing).', i)
    return re.sub(r'\s+', '', t[i:j])


def test_paper_paragraph():
    par = _paragraph()
    need = [f'${N_FOUND["M1"]}$ofthe${N_KEYS["M1"]}$and${N_FOUND["M2r"]}$ofthe${N_KEYS["M2r"]}$',
            'noneofthemofdepthtwo',
            'eachcontactisasingularpointoftheoperatorintheangleofoneofthoseconstants',
            'apurelyexponentialone', 'adifferentoneateachcontact',
            f'oforder${CARRIER_SHAPE[0]}$ineverycase',
            'akeyoperatorhas$7/2$amongitsexponentsateachcontact,itsonlyhalf-integercandidate,readmodulobothprimes',
            'itdoesnotshowthat$F_\\ell$itselfcarries',
            'northatthebranchleads',
            '2026-11-17\\_fxk\\_per\\_key\\_recurrence.md']
    for s in need:
        assert s in par, s


def test_record_results_name_the_carriers():
    t = RECORD.read_text(encoding='utf-8')
    res = t[t.index('## Results'):]
    assert '**Verdict: GO** at both M1 and M2r' in res
    o, d = CARRIER_SHAPE
    for case in ('M1', 'M2r'):
        for m in ('xZ', 'xW'):
            assert f'key ({CARRIER[case][m]}, pure exp), ({o},{d}), class M' in res, (case, m)
    assert f'| 12: 3 pure-exponential, 9 single-E1 | 15: 4 pure-exponential, 11 single-E1 |' in res
    assert f'{N_FOUND["M1"]} of {N_KEYS["M1"]} keys' in res and f'{N_FOUND["M2r"]} of {N_KEYS["M2r"]} keys' in res
    assert 'x_Z* = 61/60, x_W* = 109/60' in res and 'x_W* = 61/60, x_Z* = 829/540' in res


# ----------------------------------------------------------------------------------------- slow tier
@pytest.mark.slow
@pytest.mark.parametrize('case,p', [(c, p) for c in ('M1', 'M2r') for p in PRIMES])
def test_carrier_minimal_windowed_slow(case, p):
    """FX's own search (R-bisection from R = 1) on the carrier keys: (15, 77), verified on >= 200 held-out rows."""
    import fx_run as R
    for m in ('xZ', 'xW'):
        r = R.minimal_windowed(load_seqs(case, p)[carrier_key(case, m)], p, WINDOW, RMAX)
        assert r['found'] and r['verified'] and r['min_held_out'] >= LMAX - WINDOW
        assert (r['order'], r['degree']) == CARRIER_SHAPE
        assert tuple(r['search'][:2]) == CARRIER_SEARCH


@pytest.mark.slow
@pytest.mark.parametrize('case,p', [(c, p) for c in ('M1', 'M2r') for p in PRIMES])
def test_full_per_key_search_slow(case, p):
    """~2000 s each: the whole per-key search reproduces the found set, shapes and fingerprints."""
    import fxk_run as K
    seqs = load_seqs(case, p)
    recs, ops = K.per_key(seqs, p, WINDOW, log=lambda *a, **k: None)
    found = {k: (v['order'], v['degree']) for k, v in recs.items() if v.get('found')}
    assert len(found) == N_FOUND[case]
    assert sum(v['in_scope'] for v in recs.values()) == N_KEYS[case]
    assert not [k for k in found if is_depth2(parse_key(k))]
    hist: Dict[Tuple[int, int], int] = {}
    for sh in found.values():
        hist[sh] = hist.get(sh, 0) + 1
    assert hist == SHAPES[case]
    r = load_run(case, p)
    assert {k: v['fp'] for k, v in recs.items() if v.get('found')} == \
        {k: v['fp'] for k, v in r['per_key'].items() if v.get('found')}
    for m in ('xZ', 'xW'):
        hits = [k for k, c in ops.items() if _lc_at(c, CONTACTS[case][m], p) == 0]
        assert hits == [carrier_key(case, m)], (m, hits)
