"""Guard for record OP52 (docs/predictions/2026-11-18_aa_bc_operator_exponents.md, scored GO 2026-11-18): the contact
exponents of the bond-angle operator of

    f'(x) = int V'(|r - A|) rho_BC(r) d^3r = NA - Y_mu,      V'(s) = (1 - e^{-mu s}) / s,   mu = 2 zeta_A,

the (AA|BC)-type integral with A's density replaced by mu^2 e^{-mu s} / (4 pi s), and the counting rule of Paper 59
sec:scope ("[ARGUMENT + MEASURED] A counting rule accounts for these exponents. ... registered before computing)").
Written 2026-11-18 as a separate step after the paper edit (CLAUDE.md sec. 9, guard-writing rule).

WHAT IS RECOMPUTED AND WHAT IS READ.  The scored run files debug/data/op52_run_<case>_<p>.json store READINGS only
(orders, degrees, exponent lists); they hold no operator and no sequence.  So this guard carries its own data,
tests/data/paper59_op52_2026_11_18.json, rebuilt 2026-11-18 from the frozen engine copies in
tests/paper59_intclass_support/ (nothing is imported from debug/): for each (case, prime) the 21 constant sequences of
f' (L <= 1000), the minimal operators of the generic sums (all 21; the three D1 endpoint constants; the two D2
endpoint constants; weight seeds 2026 and 7) and the minimal operator of every constant piece.  That file is not
trusted on its own:
  * sequences: first 31 terms recomputed from the engine copy (fast), all 1001 (slow);
  * operators: each is checked to ANNIHILATE its sequences by a route that is not the search's (T-form, applied to
    Legendre coefficients); the full operator must annihilate each of the five endpoint pieces, so its x* branches do not
    depend on the weights; one gcrd is recomputed in the fast tier, all in the slow tier;
  * exponents: read by this file's own indicial code (shift to x* + t, lowest valuation, chi(n), roots in F_p,
    rational reconstruction), not by s1_controls.points; all 14 roots of the full operator are accounted for, and
    chi(9/2) != 0 is asserted directly;
  * the committed run files (sha256 pinned, CRLF folded) must then agree with what was recomputed.
Shared with the scored run and NOT independent: the coefficient engine (aabc_exact / na_exact) and the operator
guesser (intclass_lib.minimal_operator).  The record's own cap stands: roots are read mod two primes, the operator is
not certified minimal and not certified over Q.

Claims guarded, each with the wrong answer it rejects:
(1) full operator: order 14, degree 109; at x* the indicial roots are exactly 0..10, 3/2, 5/2, 7/2 (rejects a 9/2,
    a missing 5/2, f's set {3/2, 7/2, 9/2}, an irregular point, a set that depends on the weights);
(2) blocks: D1 order 11 with half-integers {3/2, 5/2}; D2 order 10 with {5/2} (rejects D2 = {7/2}, a third D1 branch);
(3) pieces: exactly five are singular at x*; D1 pieces orders 5, 6, 6, each {3/2}; D2 pieces order 10, each {5/2};
(4) ingredients (sympy): V' solves Poisson for mu^2 e^{-mu s}/(4 pi s) and has s^1 coefficient -mu^2/2; the
    potential of the bounded 1s^2 density has no s^-1 and no s^1 term and first odd term s^3; point-charge and screened
    parts each lead with s^-1;
(5) the counting rule, as a pure function of (endpoint -> potentials seen): {3/2, 7/2, 9/2} for (AA|BC), {3/2, 5/2}
    for NA, {3/2, 5/2, 7/2} with blocks {3/2, 5/2} / {5/2} for f'; compared with the operators (f': recomputed here;
    (AA|BC): record MO's embedded operators, read by the same indicial code; NA and the (AA|BC) blocks: stored
    readings, sha256 pinned).  This is a consistency check of a rule matched on two objects and tested on one; the
    rule's two steps are not derived, and the guard does not make them so;
(6) instrument control: an order-2 operator with solutions (x_t - x)^{3/2}, (x_t - x)^k (x_u - x)^{1/2} is read as
    {3/2, k} for k = 5/2 and 7/2 (rejects a reader that reports a fixed set, or cannot see 5/2 absent);
(7) decay: slope of F_L / rho^L over L = 200..240 near -3 and on the 5/2 toy, about 1 away from the 3/2 and 7/2 toys;
    measured / predicted 5/2 amplitude extrapolates to 1 within 1e-3 (four terms, registered) and 1e-7 (five terms),
    with Phi5 typed from the record and recomputed from the segment formula;
(8) the paper paragraph and the record's Results carry the readings and the caveat.

Fire-test record: the comment block at the end of this file.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import random
import re
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import flint
import mpmath as mp
import numpy as np
import pytest
import sympy as sp

HERE = Path(__file__).parent
ROOT = HERE.parent
SUP = HERE / 'paper59_intclass_support'
EMBED = HERE / 'data' / 'paper59_op52_2026_11_18.json'
MO_EMBED = HERE / 'data' / 'paper59_mo_2026_11_17.json'
DATA = ROOT / 'debug' / 'data'
PAPER = ROOT / 'papers' / 'group2_quantum_chemistry' / 'paper_59_elliptic_bessel_moment.tex'
RECORD = ROOT / 'docs' / 'predictions' / '2026-11-18_aa_bc_operator_exponents.md'

# sha256 after CRLF -> LF, as tests/test_paper59_aa_bc_minimal_operator.py
FROZEN = {'s1_controls.py': 'a150322d746cb445b3b47e8264cb51a53cfe1dd81b21f8f1ab00eb8f21049453',
          'intclass_lib.py': '2a142e5b5667b432c7de4be885f5b6767daf3d209c86261076dc6ae5e4f641d4',
          'na_modguess.py': '84691c87a4a280b6823438af574c0fff6af9bf3f7ba7b76d93f2bbe866847d7d',
          'na_modanalyze.py': 'de92dc3acdcc7082001e32b50ef44db8c3bd75d990f615603b5f97f6ea178005',
          'na_exact.py': '85fc39f18bfe9231f47925436b72258dca7afddc6f9619257d9fafc708004afb',
          'aabc_exact.py': '218b7283f525d25ed8a87df61f0749ca6f576b20abcc4bde1272de32436b2d30'}
SHA = {'op52_run_A1_1073741441.json': '0c43355ec787e159d9142fcb17aac1507c8a4ab12946a2a1a6c1039f6e00a3c5',
       'op52_run_A1_1073741381.json': '3e89dd2157a00ba0ce49daf964cb6946baedc0bfe1ef01404f533b0c8caa2c74',
       'op52_run_A4_1073741441.json': 'cdc26a999c598a13e4137fc917004b8411343294fc2f4b493eaf70f1d00b5c13',
       'op52_run_A4_1073741381.json': '05d8d6a0c21bb115bdfae73aae32ee78c054b3f7f85e8fada480bbb294706a40',
       'op52_beta_A1_2000.json': 'e37d69763c43431ee98460236592fc4c7263c0bfceacb7d299409c223005362a',
       'op52_beta_A1_2600.json': '3816cd1119bc0e39cc59efe2a75e97f09ccd7d8bb09ad1aa77be6d469485a523',
       'op52_beta_A4_2000.json': '3476905d9bec48a2f1d21fcaafd435f1f3f8e144b552a81c6145226cc3a5a00b',
       'op52_beta_A4_2600.json': '30baa2d9fb3310ee69a7eaf678e32ad7cf8dbffcfb0f06c764e19875a5199cb8',
       'op52_na_1073741441.json': 'd6cdf15e453e3d1e189961fecbe0b956cbe4897f5f5bb911850c9b7280d29577',
       'op52_na_1073741381.json': '2c780ef8ece476eb62ec77c80e72474ad5e7dad590726683b689616bf2a1b8b0',
       'op52_subsets_A1_1073741441.json': 'c88a26b5f3db6d14e7b507b9cbe4f77e93329e17b3a34fd3e5f95dfc1d328f8e',
       'op52_subsets_A1_1073741381.json': '92155362e852ef5a6d6f7eaa6f314314e024bfdfa3c1d8b044ef30a3a10f1b73',
       'op52_subsets_A2_1073741441.json': 'f9e75b87f48aff3513c29c56448cd2433fed63d130444b6403dbf436f3c278d6',
       'op52_subsets_A2_1073741381.json': '14484c6b43c6b528452f6f836e3c62257c5d9b2383afe84d26027dd9b15a36ff'}

PRIMES = (1073741441, 1073741381)
SEEDS = (2026, 7)
# the record's table "Object": D1, D2, zB, zC, zA (mu = 2 zA)
CASES = {'A1': ('9/5', '13/5', '13/10', '1', '9/10'),
         'A4': ('8/5', '5/2', '13/10', '4/5', '7/10')}
XSTAR = {'A1': Fraction(125, 117), 'A4': Fraction(881, 800)}
CP = [(c, p) for c in CASES for p in PRIMES]

F32, F52, F72, F92 = (Fraction(n, 2) for n in (3, 5, 7, 9))
ORDER, DEGREE = 14, 109
FULL_HALVES = [F32, F52, F72]
FULL_ROOTS = sorted([Fraction(k) for k in range(11)] + FULL_HALVES)
D1_ORDER, D1_HALVES = 11, [F32, F52]
D2_ORDER, D2_HALVES = 10, [F52]
D1_PIECE_ORDERS, D2_PIECE_ORDERS = [5, 6, 6], [10, 10]
D1_PIECE_HALVES, D2_PIECE_HALVES = [F32], [F52]
MO_ROOTS = sorted([Fraction(k) for k in range(11)] + [F32, F72, F92])      # record MO: (AA|BC) at x*
L_PREFIX = 30

# typed from the record (OP5), not read from any file
PHI5 = {'A1': '-0.81014798223614941827', 'A4': '-0.2287660499486816037'}
LS = list(range(120, 241, 20))
DPS = 40


def _load_support():
    sys.path.insert(0, str(SUP))
    try:
        mods = []
        for name in ('s1_controls', 'aabc_exact'):
            spec = importlib.util.spec_from_file_location(f'p59op52_{name}', SUP / f'{name}.py')
            m = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(m)
            mods.append(m)
    finally:
        sys.path.remove(str(SUP))
    return mods


C, X = _load_support()
I = C.I


def _folded(path):
    return Path(path).read_bytes().replace(b'\r\n', b'\n')


@lru_cache(maxsize=None)
def _json(path):
    return json.loads(_folded(path).decode())


def embedded(case, p):
    return _json(str(EMBED))['entries'][f'{case}|{p}']


def stored(name):
    return _json(str(DATA / name))


def params(case):
    return tuple(Fraction(v) for v in CASES[case])


def endpoint_keys(case):
    """repr of the engine's keys for the endpoint constants: D1 = e^{-zB D1} {1, e^{-mu (D2 -+ D1)}},
    D2 = e^{+-zC D1 - (zB + zC) D2} (the record's 'Object' table)."""
    D1, D2, zB, zC, zA = params(case)
    mu = 2 * zA
    one = lambda e: repr((Fraction(e), '1', None))                                # noqa: E731
    return ([one(-zB * D1), one(-zB * D1 - mu * (D2 - D1)), one(-zB * D1 - mu * (D2 + D1))],
            [one(zC * D1 - (zB + zC) * D2), one(-zC * D1 - (zB + zC) * D2)])


def fprime_sequences(case, p, lmax):
    """The 21 constant sequences of f' = NA - Y_mu mod p from the frozen engine copy (kernels NA and Y; no dY)."""
    K = X.E.Field(p)
    pc = X.pieces('A1', lmax, K, which=('NA', 'Y'), params_override=params(case))
    out = []
    for qN, qY in zip(pc['NA'], pc['Y']):
        r = dict(qN)
        for k, v in qY.items():
            X.E.add(r, k, -v)
        out.append({k: v for k, v in r.items() if v != 0})
    keys = sorted(set(k for q in out for k in q), key=str)
    return {repr(k): [int(q.get(k, flint.nmod(0, p))) for q in out] for k in keys}


def polys(lists, p):
    return [flint.nmod_poly(c, p) for c in lists]


def op_lists(o):
    return [[int(a) for a in ck.coeffs()] for ck in o['c']]


# ================================================================================ own indicial reader (not the run's)
def fp(q, p):
    q = Fraction(q)
    return q.numerator % p * pow(q.denominator % p, p - 2, p) % p


def small_rational(a, p, dmax=64, nmax=4096):
    """The rational n/d = a mod p with d <= dmax, |n| <= nmax (unique: 2 dmax nmax < p), else None."""
    for d in range(1, dmax + 1):
        n = a * d % p
        if n > p // 2:
            n -= p
        if abs(n) <= nmax:
            return Fraction(n, d)
    return None


def falling(k, p):
    """n (n-1) ... (n-k+1) as a polynomial in n."""
    P = flint.nmod_poly([1], p)
    for j in range(k):
        P = P * flint.nmod_poly([(-j) % p, 1], p)
    return P


def indicial(c, x0, p):
    """Indicial polynomial at x0 of L = sum_k c_k(x) (d/dx)^k.  Each c_k is expanded at x0 + t; with v_k its
    t-valuation and a_k its lowest coefficient, L t^n = t^(n+m) [chi(n) + O(t)], m = min_k (v_k - k),
    chi(n) = sum_{k: v_k - k = m} a_k n(n-1)...(n-k+1).  Returns (chi, regular): regular iff the top order attains m
    (then deg chi = order).  None if x0 is an ordinary point of the leading coefficient."""
    shift = flint.nmod_poly([fp(x0, p), 1], p)
    low = []
    for k, ck in enumerate(c):
        if ck.degree() < 0:
            continue
        co = [int(v) for v in ck(shift).coeffs()]
        v = next(i for i, a in enumerate(co) if a)
        low.append((k, v, co[v]))
    R = len(c) - 1
    if dict((k, v) for k, v, _ in low)[R] == 0:
        return None
    m = min(v - k for k, v, _ in low)
    chi = flint.nmod_poly([], p)
    for k, v, a in low:
        if v - k == m:
            chi = chi + falling(k, p) * a
    return chi, chi.degree() == R


def roots_at(c, x0, p):
    """Sorted indicial roots at x0 with multiplicity, as Fractions; asserts the point is regular singular and that
    chi splits completely over F_p into small rationals (every root accounted for)."""
    got = indicial(c, x0, p)
    assert got is not None, 'not a singular point'
    chi, regular = got
    assert regular, 'irregular point'
    out = []
    for r, mult in chi.roots():
        q = small_rational(int(r), p)
        assert q is not None, int(r)
        out += [q] * int(mult)
    assert len(out) == len(c) - 1, ('roots unaccounted for', len(out), len(c) - 1)
    return sorted(out)


def is_root(c, x0, p, e):
    chi, _ = indicial(c, x0, p)
    return int(chi(flint.nmod(fp(e, p), p))) == 0


def halves(roots):
    return [r for r in roots if r.denominator == 2]


def singular_at(c, x0, p):
    """x0 is a root of the leading coefficient."""
    return int(c[-1](flint.nmod(fp(x0, p), p))) == 0


# ============================================================ independent annihilation (T-form; not the search's code)
def T_form(c, p):
    """(1 - x^2)^R L = sum_i a_i(x) T^i, T = (1 - x^2) d/dx.  T^k = sum_m t_{k,m} d^m with
    t_{k+1,m} = (1 - x^2)(t_{k,m}' + t_{k,m-1}); peeled from the top.  None if a division is inexact."""
    R = len(c) - 1
    w = flint.nmod_poly([1, 0, p - 1], p)
    zero = flint.nmod_poly([], p)
    t = [[flint.nmod_poly([1], p)]]
    for k in range(R):
        t.append([w * ((t[k][m].derivative() if m <= k else zero) + (t[k][m - 1] if m >= 1 else zero))
                  for m in range(k + 2)])
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


def apply_T(a, seq, p):
    """Legendre coefficients of sum_i a_i(x) T^i f, f = sum seq_L P_L, on the range where they are determined
    (x P_L and T P_L three-term rules)."""
    n = len(seq) + 2
    L = np.arange(n, dtype=np.int64)
    inv = lambda v: np.array([pow(int(z) % p, p - 2, p) for z in v], dtype=np.int64)     # noqa: E731
    xm, xp = L % p * inv(2 * L - 1) % p, (L + 1) % p * inv(2 * L + 3) % p
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


def annihilates(c, seqs, p, min_len=600):
    a = T_form(c, p)
    assert a is not None
    for sq in seqs:
        r = apply_T(a, sq, p)
        assert len(r) >= min_len, len(r)
        if np.any(r):
            return False
    return True


# ====================================================================================================== (4) potentials
_s, _mu, _z = sp.symbols('s mu zeta', positive=True)
POT = {'point': 1 / _s,                                        # 1/s: the point-charge part (NA's kernel)
       'yukawa': sp.exp(-_mu * _s) / _s,                       # Y_mu alone
       'screened': (1 / _s + _z) * sp.exp(-2 * _z * _s),       # S: the screened part of the 1s^2 potential
       'bounded': 1 / _s - (1 / _s + _z) * sp.exp(-2 * _z * _s),   # potential of the 1s^2 density (bounded)
       'vprime': (1 - sp.exp(-_mu * _s)) / _s}                 # potential of mu^2 e^{-mu s}/(4 pi s)
ODD = (-1, 1, 3, 5)


@lru_cache(maxsize=None)
def odd_coefficients(name):
    ser = sp.series(POT[name], _s, 0, ODD[-1] + 2).removeO()
    return tuple(sp.simplify(ser.coeff(_s, n)) for n in ODD)


def test_potential_of_the_singular_density_has_a_linear_term():
    """V' = (1 - e^{-mu s})/s is the potential of rho = mu^2 e^{-mu s}/(4 pi s) (Poisson, unit charge, finite at 0);
    its s^-1 coefficient is 0 and its s^1 coefficient is -mu^2/2, not zero.  Wrong answers rejected: a potential
    with no s^1 term (then no 5/2); a surviving point charge; a density that is not the one named."""
    V = POT['vprime']
    rho = _mu ** 2 * sp.exp(-_mu * _s) / (4 * sp.pi * _s)
    assert sp.simplify(-sp.diff(_s * V, _s, 2) / _s - 4 * sp.pi * rho) == 0
    assert sp.simplify(sp.integrate(4 * sp.pi * _s ** 2 * rho, (_s, 0, sp.oo)) - 1) == 0
    assert sp.limit(V * _s, _s, sp.oo) == 1 and sp.limit(V, _s, 0) == _mu
    v = odd_coefficients('vprime')
    assert v[0] == 0
    assert sp.simplify(v[1] + _mu ** 2 / 2) == 0 and v[1] != 0
    assert sp.simplify(v[2] + _mu ** 4 / 24) == 0
    assert sp.simplify(sp.limit(sp.diff(V, _s), _s, 0) + _mu ** 2 / 2) == 0       # not through the series routine
    assert sp.simplify(POT['point'] - POT['yukawa'] - V) == 0                     # V' = NA's kernel - Y's


def test_potential_of_the_bounded_density_has_no_linear_term_and_its_parts_are_point_charges():
    """The 1s^2 density zeta^3 e^{-2 zeta s}/pi: its potential has zero s^-1 and s^1 coefficients and first odd term
    s^3 (2 zeta^4/3); the point-charge part and the screened part each lead with s^-1, and the screened part has
    no s^1 term.  Wrong answers rejected: an s^1 term in the bounded potential (5/2 would arise in (AA|BC)); a
    screened part that is not a point charge at short range."""
    V = POT['bounded']
    rho = _z ** 3 / sp.pi * sp.exp(-2 * _z * _s)
    assert sp.simplify(-sp.diff(_s * V, _s, 2) / _s - 4 * sp.pi * rho) == 0
    b = odd_coefficients('bounded')
    assert b[0] == 0 and b[1] == 0
    assert sp.simplify(b[2] - sp.Rational(2, 3) * _z ** 4) == 0 and b[2] != 0
    assert sp.limit(sp.diff(V, _s), _s, 0) == 0
    assert odd_coefficients('point') == (1, 0, 0, 0)
    sc = odd_coefficients('screened')
    assert sc[0] == 1 and sc[1] == 0 and sc[2] != 0
    y = odd_coefficients('yukawa')
    assert y[0] == 1 and sp.simplify(y[1] - _mu ** 2 / 2) == 0
    assert sp.simplify(POT['point'] - POT['screened'] - V) == 0


# ================================================================================================ (5) the counting rule
def endpoint_exponents(potentials):
    """Half-integer exponents carried by the pieces at one endpoint.  One branch per potential seen there; the
    branch of s^p is (x* - x)^{(p+4)/2}; potentials with the same leading power are resolved inside the endpoint by
    their difference, i.e. the exponents are the leading odd powers of a valuation-adapted basis of the span
    (row echelon form of the odd-coefficient vectors)."""
    M = sp.Matrix([list(odd_coefficients(name)) for name in potentials])
    _, pivots = M.rref(simplify=True)
    return sorted(Fraction(ODD[j] + 4, 2) for j in pivots)


def resolve_collisions(exponents):
    """Branches of different endpoints with the same exponent differ at the next order."""
    out = []
    for e in sorted(exponents):
        while e in out:
            e += 1
        out.append(e)
    return out


def counting_rule(endpoints):
    """endpoints: {name: [potentials seen by the pieces there]} -> (blocks, full)."""
    blocks = {name: endpoint_exponents(pots) for name, pots in endpoints.items()}
    return blocks, resolve_collisions([e for b in blocks.values() for e in b])


# what each endpoint's pieces see (the paper paragraph; record steps 2-3 and K3)
OBJECTS = {'AABC': {'D1': ['point', 'screened'], 'D2': ['bounded']},
           'NA': {'D1': ['point'], 'D2': ['point']},
           'FPRIME': {'D1': ['point', 'yukawa'], 'D2': ['vprime']}}


def test_counting_rule_gives_the_three_exponent_sets():
    """The rule returns {3/2, 7/2, 9/2} for (AA|BC) with blocks {3/2, 7/2} / {7/2}; {3/2, 5/2} for NA; and
    {3/2, 5/2, 7/2} with blocks {3/2, 5/2} / {5/2} for f'.  Wrong answers rejected: a rule without the collision
    step (f' would read {3/2, 5/2}, (AA|BC) {3/2, 7/2}, NA {3/2}); a rule blind to the s^1 term (f' would read
    (AA|BC)'s set); two branches for one potential.  Not rejected, by construction: another rule that agrees on these
    three objects (the record: the two readings of the collision step coincide here)."""
    blocks, full = counting_rule(OBJECTS['AABC'])
    assert blocks == {'D1': [F32, F72], 'D2': [F72]} and full == [F32, F72, F92]
    blocks, full = counting_rule(OBJECTS['NA'])
    assert blocks == {'D1': [F32], 'D2': [F32]} and full == [F32, F52]
    blocks, full = counting_rule(OBJECTS['FPRIME'])
    assert blocks == {'D1': D1_HALVES, 'D2': D2_HALVES} and full == FULL_HALVES
    assert F92 not in full and F52 not in counting_rule(OBJECTS['AABC'])[1]


@pytest.mark.parametrize('case,p', CP)
def test_counting_rule_agrees_with_the_fprime_operators(case, p):
    """The rule's f' prediction against the embedded operators read by this file's indicial code (full and both
    blocks, both weight seeds).  Wrong answer rejected: a rule output that the operators do not show."""
    blocks, full = counting_rule(OBJECTS['FPRIME'])
    e, xs = embedded(case, p), XSTAR[case]
    for sd in SEEDS:
        assert halves(roots_at(polys(e['ops'][f'full|{sd}'], p), xs, p)) == full
        assert halves(roots_at(polys(e['ops'][f'D1|{sd}'], p), xs, p)) == blocks['D1']
        assert halves(roots_at(polys(e['ops'][f'D2|{sd}'], p), xs, p)) == blocks['D2']


def test_counting_rule_agrees_with_record_MO_and_the_stored_controls():
    """(AA|BC): record MO's embedded minimal operators (A1 and A2, both primes), read by this file's indicial code,
    have roots 0..10, 3/2, 7/2, 9/2 at x* (so the reader does not return f's set for every operator), which is the
    rule's set.  Blocks of (AA|BC) and the NA instance are STORED readings (controls K1, K2, K4; sha256 pinned).
    Wrong answers rejected: a rule set that differs from MO's operator; a reader that always returns {3/2, 5/2, 7/2}."""
    blocks, full = counting_rule(OBJECTS['AABC'])
    mo = _json(str(MO_EMBED))['entries']
    for case, xs in (('A1', Fraction(125, 117)), ('A2', Fraction(89, 80))):
        for p in PRIMES:
            c = polys(mo[f'{case}|{p}|base|2026']['op'], p)
            roots = roots_at(c, xs, p)
            assert roots == MO_ROOTS and halves(roots) == full
            assert not is_root(c, xs, p, F52)
            sub = stored(f'op52_subsets_{case}_{p}.json')['subsets']
            H = lambda k: [Fraction(v) for v in sub[k]['halves']]                 # noqa: E731
            assert H('trio_D1') == blocks['D1'] and H('pair_D2') == blocks['D2'] and H('all5') == full
    blocks, full = counting_rule(OBJECTS['NA'])
    for p in PRIMES:
        sub = stored(f'op52_na_{p}.json')['subsets']
        H = lambda k: [Fraction(v) for v in sub[k]['halves']]                     # noqa: E731
        assert H('D1') == blocks['D1'] and H('pair') == blocks['D2']
        assert H('D1_plus_pair') == full and H('all') == full


# ============================================================================================ (6) instrument control
def rank2_operator(k, p):
    """Order-2 operator y'' + p1 y' + p0 y with solutions y1 = (x_t - x)^{3/2}, y2 = (x_t - x)^k (x_u - x)^{1/2},
    x_t = 5/3, x_u = -3, from the logarithmic derivatives a = y1'/y1, b = y2'/y2:
    p1 = -[(b' + b^2) - (a' + a^2)] / (b - a),  p0 = [a (b' + b^2) - b (a' + a^2)] / (b - a); denominators cleared."""
    x = sp.Symbol('x')
    xt, xu = sp.Rational(5, 3), sp.Integer(-3)
    a = sp.Rational(3, 2) / (x - xt)
    b = sp.Rational(k) / (x - xt) + sp.Rational(1, 2) / (x - xu)
    A2, B2 = sp.diff(a, x) + a ** 2, sp.diff(b, x) + b ** 2
    p1 = sp.cancel(-(B2 - A2) / (b - a))
    p0 = sp.cancel((a * B2 - b * A2) / (b - a))
    for u, U2 in ((a, A2), (b, B2)):                               # y''/y + p1 y'/y + p0 = 0 for both solutions
        assert sp.cancel(U2 + p1 * u + p0) == 0
    den = sp.lcm(sp.denom(sp.together(p1)), sp.denom(sp.together(p0)))
    co = [sp.Poly(sp.cancel(v * den), x) for v in (p0, p1, sp.Integer(1))]
    return [flint.nmod_poly([fp(Fraction(int(q.p), int(q.q)), p) for q in reversed(P.all_coeffs())], p) for P in co]


@pytest.mark.parametrize('k', ['5/2', '7/2'])
@pytest.mark.parametrize('p', PRIMES)
def test_indicial_reader_on_a_rank_two_operator_with_known_exponents(k, p):
    """The reader returns exactly {3/2, k} at x_t = 5/3 (k = 7/2: 5/2 read absent) and {0, 1/2} at x_u = -3; and an
    ordinary point is reported as such.  Wrong answers rejected: a reader that reports {3/2, 5/2} or any fixed set;
    one that cannot tell 5/2 from 7/2; one that reads every point as singular."""
    c = rank2_operator(k, p)
    assert len(c) == 3
    assert roots_at(c, Fraction(5, 3), p) == sorted([F32, Fraction(k)])
    other = F72 if Fraction(k) == F52 else F52
    assert not is_root(c, Fraction(5, 3), p, other)
    assert roots_at(c, Fraction(-3), p) == [Fraction(0), Fraction(1, 2)]
    assert indicial(c, Fraction(1, 7), p) is None


# ================================================================================= (1)-(3) the operators, fast tier
def test_support_copies_are_the_frozen_instrument():
    """Wrong answer rejected: an engine or guesser edited after the record was scored."""
    for name, h in FROZEN.items():
        assert hashlib.sha256(_folded(SUP / name)).hexdigest() == h, name


@pytest.mark.parametrize('name', sorted(SHA))
def test_committed_data_files_are_the_scored_ones(name):
    """Wrong answer rejected: a run, beta or control file edited after scoring."""
    assert hashlib.sha256(_folded(DATA / name)).hexdigest() == SHA[name]


def test_geometry_and_endpoint_constants():
    """x* = (D1^2 + D2^2)/(2 D1 D2) is 125/117 and 881/800; the five endpoint constants are the record's, distinct,
    and among the 21 embedded constants.  Wrong answers rejected: another geometry; an endpoint constant that the
    engine does not produce."""
    want = {'A1': (['-117/50', '-189/50', '-513/50'], ['-209/50', '-389/50']),
            'A4': (['-52/25', '-167/50', '-391/50'], ['-397/100', '-653/100'])}
    for case in CASES:
        D1, D2 = params(case)[:2]
        assert (D1 ** 2 + D2 ** 2) / (2 * D1 * D2) == XSTAR[case]
        d1, d2 = endpoint_keys(case)
        assert d1 == [repr((Fraction(v), '1', None)) for v in want[case][0]]
        assert d2 == [repr((Fraction(v), '1', None)) for v in want[case][1]]
        for p in PRIMES:
            e = embedded(case, p)
            assert len(e['keys']) == 21 == len(set(e['keys'])) and sorted(e['seqs']) == sorted(e['keys'])
            assert set(d1 + d2) <= set(e['keys']) and len(set(d1 + d2)) == 5
            assert all(len(v) == 1001 for v in e['seqs'].values())


@pytest.mark.parametrize('case,p', CP)
def test_embedded_sequences_match_engine_prefix(case, p):
    """The embedded sequences are f' = NA - Y_mu at this case and prime: terms L <= 30 of all 21 constants from the
    frozen engine copy.  Wrong answers rejected: data of another geometry, prime or kernel combination (e.g. f with
    dY, or NA alone)."""
    e = embedded(case, p)
    seqs = fprime_sequences(case, p, L_PREFIX)
    assert sorted(seqs) == sorted(e['keys'])
    for k, sq in seqs.items():
        assert sq == e['seqs'][k][:L_PREFIX + 1], k
    assert any(any(sq[2:]) for sq in seqs.values())


@pytest.mark.parametrize('case,p', CP)
def test_full_operator_order_degree_and_annihilation(case, p):
    """Claim (1), structure: the embedded operators of the generic sums (both seeds) have order 14, degree 109; each
    annihilates its generic sequence, and each of the five endpoint pieces separately, so the x* branches belong to
    the operator whatever the weights.  (Measured here, not a paper claim: 19 of the 21 pieces are annihilated one by
    one; the two that are not are non-endpoint pure exponentials, regular at x*, so the operator of the generic sum
    depends on the weights only through them.)  Negative control: the other case's operator does not annihilate
    these sequences.  Wrong answers rejected: order 13 or 17, degree 111 (f's); a corrupted operator; an operator of
    another object; x* exponents that are an artifact of one weight vector."""
    e = embedded(case, p)
    d1, d2 = endpoint_keys(case)
    for sd in SEEDS:
        c = polys(e['ops'][f'full|{sd}'], p)
        assert len(c) - 1 == ORDER and max(ck.degree() for ck in c) == DEGREE, sd
        gen = C.generic_seq({k: e['seqs'][k] for k in e['keys']}, p, sd)
        assert annihilates(c, [gen] + [e['seqs'][k] for k in d1 + d2], p)
        missed = [k for k in e['keys'] if not annihilates(c, [e['seqs'][k]], p)]
        assert len(missed) == 2 and all(k.endswith("'1', None)") and k not in d1 + d2 for k in missed), missed
    other = next(cc for cc in CASES if cc != case)
    co = polys(embedded(other, p)['ops']['full|2026'], p)
    assert not annihilates(co, [gen], p) and not annihilates(co, [e['seqs'][d2[0]]], p)


@pytest.mark.parametrize('case,p', CP)
def test_full_operator_indicial_roots_at_xstar(case, p):
    """Claim (1): x* is a regular singular point of the full operator and its 14 indicial roots are exactly 0..10,
    3/2, 5/2, 7/2, each simple, under both weight seeds; chi(9/2) != 0 and chi(1/2) != 0 directly.  Wrong answers
    rejected: 9/2 present; 5/2 absent ((AA|BC)'s set); an unaccounted root; an irregular point."""
    e, xs = embedded(case, p), XSTAR[case]
    for sd in SEEDS:
        c = polys(e['ops'][f'full|{sd}'], p)
        roots = roots_at(c, xs, p)
        assert roots == FULL_ROOTS, (sd, roots)
        assert halves(roots) == FULL_HALVES
        assert not is_root(c, xs, p, F92) and not is_root(c, xs, p, Fraction(1, 2))
        assert all(is_root(c, xs, p, h) for h in FULL_HALVES)


@pytest.mark.parametrize('case,p', CP)
def test_blocks_orders_annihilation_and_exponents(case, p):
    """Claim (2): the D1 block operator has order 11, annihilates the three D1 pieces, and has half-integer roots
    {3/2, 5/2} at x*; the D2 block has order 10, annihilates the two D2 pieces, and has {5/2} alone; both seeds; the
    D2 block does not annihilate a D1 piece.  Wrong answers rejected: D2 = {7/2} (the jump at D2 not V' alone);
    a third D1 branch (Y's halves separating); 3/2 in the D2 block."""
    e, xs = embedded(case, p), XSTAR[case]
    d1, d2 = endpoint_keys(case)
    for sd in SEEDS:
        c1, c2 = polys(e['ops'][f'D1|{sd}'], p), polys(e['ops'][f'D2|{sd}'], p)
        assert len(c1) - 1 == D1_ORDER and len(c2) - 1 == D2_ORDER
        assert annihilates(c1, [e['seqs'][k] for k in d1], p)
        assert annihilates(c2, [e['seqs'][k] for k in d2], p)
        assert not annihilates(c2, [e['seqs'][d1[0]]], p)
        r1, r2 = roots_at(c1, xs, p), roots_at(c2, xs, p)
        assert halves(r1) == D1_HALVES and halves(r2) == D2_HALVES, (sd, r1, r2)
        assert all(r.denominator == 1 and r >= 0 for r in r1 + r2 if r.denominator != 2)
        assert not is_root(c2, xs, p, F32) and not is_root(c2, xs, p, F72) and not is_root(c1, xs, p, F72)


@pytest.mark.parametrize('case,p', CP)
def test_exactly_five_pieces_are_singular_at_xstar(case, p):
    """Claim (3): every embedded piece operator annihilates its own sequence; x* is a root of the leading
    coefficient for exactly the five endpoint pieces; the D1 pieces have orders 5, 6, 6 and {3/2}, the D2 pieces
    order 10 and {5/2}.  Wrong answers rejected: a sixth carrier; an endpoint piece regular at x*; a D2 piece with
    7/2 (f's reading) or 3/2."""
    e, xs = embedded(case, p), XSTAR[case]
    d1, d2 = endpoint_keys(case)
    ops = {k: polys(e['piece_ops'][k], p) for k in e['keys']}
    for k, c in ops.items():
        assert annihilates(c, [e['seqs'][k]], p, min_len=400), k
    assert sorted(k for k, c in ops.items() if singular_at(c, xs, p)) == sorted(d1 + d2)
    assert [len(ops[k]) - 1 for k in d1] == D1_PIECE_ORDERS
    assert [len(ops[k]) - 1 for k in d2] == D2_PIECE_ORDERS
    for k in d1:
        assert halves(roots_at(ops[k], xs, p)) == D1_PIECE_HALVES, k
    for k in d2:
        assert halves(roots_at(ops[k], xs, p)) == D2_PIECE_HALVES, k


@pytest.mark.parametrize('case,p', CP)
def test_committed_run_files_agree_with_the_recomputed_readings(case, p):
    """The scored run file (sha256 pinned above) stores the readings this guard recomputes: (14, 109) under both
    seeds, the 14 exponents, the blocks and their orders, and the five carriers with their orders and half-integers.
    Wrong answers rejected: a run file of another case or prime; a stored reading the operators do not support."""
    r = stored(f'op52_run_{case}_{p}.json')
    assert r['case'] == case and r['prime'] == p and Fraction(r['xstar']) == XSTAR[case]
    assert r['n_const'] == 21 and r['keys_ok'] is True and r['found'] is True and 'error' not in r
    assert r['full_od'] == [ORDER, DEGREE] == r['full_od_seed2']
    assert sorted(Fraction(v) for v in r['full_exps']) == FULL_ROOTS
    F = lambda v: [Fraction(z) for z in v]                                        # noqa: E731
    assert F(r['full'][0]) == FULL_HALVES == F(r['full'][1])
    assert F(r['D1']) == D1_HALVES == F(r['D1_s2']) and r['D1_order'] == D1_ORDER == r['D1_s2_order']
    assert F(r['D2']) == D2_HALVES == F(r['D2_s2']) and r['D2_order'] == D2_ORDER == r['D2_s2_order']
    d1, d2 = endpoint_keys(case)
    assert sorted(r['O3_carriers']) == sorted(d1 + d2) and r['O3_exactly_endpoint_five'] is True
    assert [r['O3_carriers'][k]['order'] for k in d1] == D1_PIECE_ORDERS
    assert [r['O3_carriers'][k]['order'] for k in d2] == D2_PIECE_ORDERS
    assert all(F(r['O3_carriers'][k]['halves']) == D1_PIECE_HALVES for k in d1)
    assert all(F(r['O3_carriers'][k]['halves']) == D2_PIECE_HALVES for k in d2)


def test_one_block_operator_recomputed_by_gcrd():
    """The guesser on embedded data, fast tier: at A4 (the fresh geometry), prime 1073741441, the gcrd of the
    annihilators of the D2 generic sequence equals the embedded operator.  (All operators: slow tier.)  Wrong answer
    rejected: an embedded operator that the frozen routine does not return."""
    case, p = 'A4', PRIMES[0]
    e = embedded(case, p)
    d2 = endpoint_keys(case)[1]
    o = C.generic_operator(C.generic_seq({k: e['seqs'][k] for k in d2}, p, SEEDS[0]), p)
    assert op_lists(o) == e['ops'][f'D2|{SEEDS[0]}']


# ====================================================================================================== (7) the decay
def _mf(text):
    q = Fraction(text)
    return mp.mpf(q.numerator) / q.denominator


def legendre_coefficient(gam, xs, L):
    """(2L+1)/2 int_{-1}^{1} (xs - x)^gam P_L(x) dx, xs > 1: with t = xs - sqrt(xs^2 - 1),
    (xs - x)^gam = (2t)^{-gam} sum_n C_n^{(-gam)}(x) t^n and the Gegenbauer -> Legendre connection formula
    (as tests/test_paper59_contact_amplitudes.py, where it is checked against quadrature)."""
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


def _ft(q):
    """|s|^q has Fourier transform _ft(q) |k|^{-q-3} in three dimensions."""
    q = mp.mpf(q)
    return 2 ** (q + 3) * mp.pi ** 1.5 * mp.gamma((q + 3) / 2) / mp.gamma(-q / 2)


def phi5(case):
    """The record's OP5 amplitude, coefficient of (x* - x)^{5/2}: the s^1 term of V' (v_1 = -mu^2/2) against the
    |r - C|^1 cusp of the normalized 1s_B 1s_C product (coefficient N (-zC), N = sqrt(zB^3 zC^3)/pi), pinch constant
    C_11 = ft(1)^2 / ft(5) / Beta(2, 2), |a|^2 = 2 D1 D2 (x* - x), weight e^{-zB r_B} averaged over the segment
    A -> C with u (1 - u), r_B(u)^2 = (1 - u) D2^2 + u D1^2."""
    D1, D2, zB, zC, zA = (_mf(v) for v in CASES[case])
    mu = 2 * zA
    q = sp.Rational(odd_coefficients('vprime')[1].subs(_mu, 2 * sp.Rational(CASES[case][4])))     # v_1 = -mu^2/2
    v1 = mp.mpf(int(q.p)) / int(q.q)
    assert abs(v1 / (-mu ** 2 / 2) - 1) < mp.mpf('1e-30')
    cusp = mp.sqrt(zB ** 3 * zC ** 3) / mp.pi * (-zC)
    c11 = _ft(1) ** 2 / _ft(5) / mp.beta(2, 2)
    seg = mp.quad(lambda u: u * (1 - u) * mp.exp(-zB * mp.sqrt((1 - u) * D2 ** 2 + u * D1 ** 2)),
                  [0, mp.mpf(1) / 4, mp.mpf(1) / 2, mp.mpf(3) / 4, 1])
    return v1 * cusp * c11 * (2 * D1 * D2) ** (mp.mpf(5) / 2) * seg


def extrapolate(Ls, R, terms):
    """Least squares R_L = r0 + c1/L + ... + c_{terms-1}/L^{terms-1}; returns r0."""
    rows = [[mp.mpf(1) / mp.mpf(L) ** j for j in range(terms)] for L in Ls]
    M, b = mp.matrix(terms, terms), mp.matrix(terms, 1)
    for row, val in zip(rows, R):
        for i in range(terms):
            b[i] += row[i] * val
            for j in range(terms):
                M[i, j] += row[i] * row[j]
    return mp.lu_solve(M, b)[0]


def slope(vals, xs, L1, L2):
    """d log(v_L / rho^L) / d log L between L1 and L2, rho = 1 / (x* + sqrt(x*^2 - 1))."""
    rho = 1 / (xs + mp.sqrt(xs * xs - 1))
    return mp.log(abs(vals[L2] / vals[L1]) / rho ** (L2 - L1)) / mp.log(mp.mpf(L2) / L1)


@pytest.mark.parametrize('case', sorted(CASES))
def test_fprime_decays_as_the_five_halves_branch(case):
    """From the committed exact F_L (2000 digits): the slope over L = 200..240 is within 0.02 of the 5/2 toy's (own
    Legendre coefficients of (x* - x)^{5/2}), between -3.10 and -3.02, and at least 0.9 from the 3/2 and 7/2 toys; the
    stored beta_eff(220) is within the registered 0.10 of the 5/2 toy at L = 220, and beta_eff rises toward -3; the
    2600-digit run agrees.  Wrong answers rejected: a leading 3/2 (surviving point charge, slope -2) or 7/2 (cancelled
    s^1 term, slope -4); a stored beta that is not this function's."""
    with mp.workdps(DPS):
        b = stored(f'op52_beta_{case}_2000.json')
        hi = stored(f'op52_beta_{case}_2600.json')
        assert b['case'] == case and sorted(b['F'], key=int) == [str(L) for L in LS] and 'error' not in b
        F = {L: mp.mpf(b['F'][str(L)]) for L in LS}
        for L in LS:
            assert abs(F[L] / mp.mpf(hi['F'][str(L)]) - 1) < mp.mpf('1e-12')
        assert all(F[L] > 0 for L in LS)
        xs = mp.mpf(XSTAR[case].numerator) / XSTAR[case].denominator
        toy = {g: {L: legendre_coefficient(mp.mpf(g) / 2, xs, L) for L in (200, 219, 221, 240)} for g in (3, 5, 7)}
        sF = slope(F, xs, 200, 240)
        assert mp.mpf('-3.10') < sF < mp.mpf('-3.02'), mp.nstr(sF, 6)
        assert abs(sF - slope(toy[5], xs, 200, 240)) < mp.mpf('0.02')
        assert abs(sF - slope(toy[3], xs, 200, 240)) > mp.mpf('0.9')
        assert abs(sF - slope(toy[7], xs, 200, 240)) > mp.mpf('0.9')
        # the registered gate OP3: beta_eff(220) from F_221 / F_219 against the 5/2 toy, tolerance 0.10
        beta = {int(k): v for k, v in b['beta'].items()}
        assert abs(beta[220] - float(slope(toy[5], xs, 219, 221))) < 0.10
        assert abs(beta[220] - float(sF)) < 0.005
        bl = [beta[L] for L in sorted(beta)]
        assert all(x < y < -3 for x, y in zip(bl, bl[1:]))
        assert {int(k): v for k, v in hi['beta'].items()} == beta


@pytest.mark.parametrize('case', sorted(CASES))
def test_five_halves_amplitude(case):
    """Phi5 typed from the record equals the segment formula (1e-18), and F_L / (Phi5 T^{5/2}_L), L = 120..240 step
    20, extrapolates to 1: |r0 - 1| < 1e-3 with the registered four-term fit and < 1e-7 with five terms; the raw
    ratio at L = 240 is still 0.6-0.8 % low.  Wrong answers rejected: an amplitude off by 1e-3 or more (a wrong v_1,
    pinch constant, cusp coefficient or segment weight); a stored r0 that the stored F_L do not give."""
    with mp.workdps(DPS):
        typed = mp.mpf(PHI5[case])
        assert abs(phi5(case) / typed - 1) < mp.mpf('1e-18'), mp.nstr(phi5(case), 25)
        b = stored(f'op52_beta_{case}_2000.json')
        xs = mp.mpf(XSTAR[case].numerator) / XSTAR[case].denominator
        R = [mp.mpf(b['F'][str(L)]) / (typed * legendre_coefficient(mp.mpf(5) / 2, xs, L)) for L in LS]
        r4, r5 = extrapolate(LS, R, 4), extrapolate(LS, R, 5)
        assert abs(r4 - 1) < mp.mpf('1e-3'), mp.nstr(r4 - 1, 5)
        assert abs(r5 - 1) < mp.mpf('1e-7'), mp.nstr(r5 - 1, 5)
        assert all(x < y for x, y in zip(R, R[1:])) and mp.mpf('0.990') < R[-1] < mp.mpf('0.995')
        op5 = b['OP5']
        assert abs(mp.mpf(op5['Phi5']) / typed - 1) < mp.mpf('1e-18') and op5['holds'] is True
        assert abs(mp.mpf(op5['r0']) - r4) < mp.mpf('1e-10') and abs(mp.mpf(op5['R_240']) - R[-1]) < mp.mpf('1e-10')


# ================================================================================================ (8) paper and record
def _flat(path):
    return re.sub(r'\s+', ' ', Path(path).read_text(encoding='utf-8'))


def test_paper_paragraph_states_the_readings_and_the_caveat():
    """Paper 59 sec:scope, from '[ARGUMENT + MEASURED]} A counting rule' to the record pointer.  Wrong answers
    rejected: the caveat sentence dropped (the rule's two steps presented as derived); a 9/2 claimed for f'; the
    blocks or the two-geometry / two-prime scope missing; the tier label changed."""
    text = _flat(PAPER)
    start = text.index(r'\textbf{[ARGUMENT + MEASURED]} A counting rule accounts for these exponents.')
    end = text.index(r'(record \texttt{docs/predictions/2026-11-18\_aa\_bc\_operator\_exponents.md}, registered before '
                     'computing', start)
    par = text[start:end]
    assert len(par) < 2500
    for phrase in (r'Only five of the $21$ pieces are singular at $x^{*}$',
                   r'$\mu^2e^{-\mu s}/(4\pi s)$, whose potential has an $|\mathbf{s}|$ term',
                   r'exponents $3/2$, $5/2$ and $7/2$ with no $9/2$',
                   r'the first group with $3/2$ and $5/2$, the second with $5/2$ alone',
                   r'a leading $(x^{*}-x)^{5/2}$ branch',
                   'held at two geometries (one not used in reading the rule), modulo two primes',
                   'On this reading no $5/2$ arises because the potential of a bounded density has no $|\\mathbf{s}|$ term',
                   'The rule was read from this integral and from the one-electron one ($3/2$ and $5/2$)',
                   'Two steps of the rule, one branch per group and potential and the resolution of coinciding '
                   'branches at the next order, are supported by that prediction and not derived'):
        assert phrase in par, phrase
    assert 'derived' not in par.replace('supported by that prediction and not derived', '')
    assert 'theorem' not in par.lower() and 'proved' not in par.lower()


def test_record_results_state_the_guarded_readings():
    """The record's Results section carries the figures pinned here.  Wrong answer rejected: a record whose scored
    figures differ from the guard's constants."""
    text = _flat(RECORD)
    res = text[text.index('## Results'):]
    for phrase in ('**Verdict: GO.**',
                   '| full operator (order, degree) | (reported) | (14, 109) | (14, 109) | (14, 109) | (14, 109) |',
                   '{3/2, 5/2, 7/2}, 9/2 absent | {3/2, 5/2, 7/2} | same | same | same |',
                   '| OP2: D1 block | {3/2, 5/2} | {3/2, 5/2} | same | same | same |',
                   '| OP2: D2 block | {5/2} | {5/2} | same | same | same |',
                   'the three D1 pieces have orders 5, 6, 6, each with {3/2}; the two D2 pieces have order 10, each '
                   'with {5/2}',
                   'all 14 indicial roots are accounted for: 0..10, 3/2, 5/2, 7/2',
                   'The D1 block has order 11 and the D2 block order 10',
                   '-3.0607 against -3.0697', '-3.0504 against -3.0577',
                   'are supported by this prediction, not derived'):
        assert phrase in res, phrase
    for case in CASES:
        assert f'{PHI5[case]}' in text


# ========================================================================================================= slow tier
@pytest.mark.slow
@pytest.mark.parametrize('case,p', CP)
def test_embedded_sequences_full_recompute(case, p):
    """All 1001 terms of the 21 constant sequences from the frozen engine (about 6 min each).  Wrong answer
    rejected: embedded data that agree with f' on a prefix only."""
    e = embedded(case, p)
    seqs = fprime_sequences(case, p, 1000)
    assert sorted(seqs) == sorted(e['keys'])
    for k, sq in seqs.items():
        assert sq == e['seqs'][k], k


@pytest.mark.slow
@pytest.mark.parametrize('case,p', CP)
def test_all_operators_recomputed_by_gcrd(case, p):
    """Every embedded operator is what the frozen guesser returns on the embedded sequences: full and both blocks
    under both seeds, and all 21 piece operators (about 8 min per case and prime)."""
    e = embedded(case, p)
    d1, d2 = endpoint_keys(case)
    for name, ks in (('full', e['keys']), ('D1', d1), ('D2', d2)):
        for sd in SEEDS:
            o = C.generic_operator(C.generic_seq({k: e['seqs'][k] for k in ks}, p, sd), p)
            assert op_lists(o) == e['ops'][f'{name}|{sd}'], (name, sd)
    for k in e['keys']:
        assert op_lists(C.piece_operator(e['seqs'][k], p)) == e['piece_ops'][k], k


@pytest.mark.slow
@pytest.mark.parametrize('case,p', CP)
def test_own_weights_give_the_same_exponents(case, p):
    """A third weight vector, drawn here (not s1_controls.generic_seq): the minimal operator of that generic sum has
    order 14 and the same 14 roots at x*; likewise the blocks.  Wrong answer rejected: exponents tied to the two
    scored weight vectors."""
    e, xs = embedded(case, p), XSTAR[case]
    d1, d2 = endpoint_keys(case)
    rnd = random.Random(31337)
    w = {k: rnd.randrange(1, p) for k in e['keys']}
    mix = lambda ks: [sum(w[k] * e['seqs'][k][L] for k in ks) % p for L in range(1001)]     # noqa: E731
    o = C.generic_operator(mix(e['keys']), p)
    assert (o['order'], o['degree']) == (ORDER, DEGREE)
    assert roots_at(o['c'], xs, p) == FULL_ROOTS
    o1, o2 = C.generic_operator(mix(d1), p), C.generic_operator(mix(d2), p)
    assert (o1['order'], o2['order']) == (D1_ORDER, D2_ORDER)
    assert halves(roots_at(o1['c'], xs, p)) == D1_HALVES and halves(roots_at(o2['c'], xs, p)) == D2_HALVES


# Fire-test record, 2026-11-18.  Every plant restored; the repo's paper, record and debug/ files were never written.
# (a) in this file, python debug/qa/fire_test.py tests/test_paper59_operator_exponents.py --plant "O=>N" -k <sel>
#     (each -k leaves out the hash tests, so the behavioural test named is what fires):
#   F01 FULL_HALVES + 9/2                                     -k full_operator_indicial           FIRED
#   F02 FULL_HALVES without 5/2                               -k full_operator_indicial           FIRED
#   F02b OBJECTS['FPRIME'] := (AA|BC)'s potentials            -k agrees_with_the_fprime           FIRED
#   F03 D2_HALVES {5/2} => {7/2}                              -k blocks_orders                    FIRED
#   F04 ORDER 14 => 13                                        -k order_degree_and_annihilation    FIRED
#   F05 DEGREE 109 => 111                                     -k order_degree_and_annihilation    FIRED
#   F06 D1_ORDER 11 => 12                                     -k blocks_orders                    FIRED
#   F07 D1 piece orders 5, 6, 6 => 6, 6, 6                    -k exactly_five                     FIRED
#   F08 D2 piece halves {5/2} => {7/2}                        -k exactly_five                     FIRED
#   F09 rule: collision step removed (while False)            -k counting_rule_gives              FIRED
#   F10 V' given no s^1 term (+ mu^2 s / 2)                   -k counting_rule_gives, singular_density   FIRED, FIRED
#   F11 expected s^1 coefficient of V' => 0                   -k singular_density                 FIRED
#   F12 bounded potential given an s^1 term                   -k bounded_density, counting_rule_gives    FIRED, FIRED
#   F13 rank-2 toy built with k + 1                           -k rank_two                         FIRED
#   F14 indicial reader: 'v - k == m' => 'v - k <= m + 1'     -k rank_two, full_operator_indicial        FIRED, FIRED
#   F15 typed Phi5(A1) x (1 + 1e-3)                           -k five_halves_amplitude            FIRED
#   F16 x*(A4) 881/800 => 881/801                             -k full_operator_indicial           FIRED
#   F17 decay compared with the 7/2 toy                       -k decays_as                        FIRED
#   F18 engine combination f' := NA alone                     -k engine_prefix                    FIRED
#   F19 MO_ROOTS 9/2 => 5/2                                   -k record_MO                        FIRED
#   F20 gcrd recomputed from the D1 keys                      -k one_block_operator               FIRED
#   F20x gcrd recomputed with weight seed 7 instead of 2026   -k one_block_operator               DID NOT FIRE: not a
#        wrong answer.  Both D2 pieces are annihilated by the block operator, so it is the same for every weight
#        vector (the embedded D2|2026 and D2|7 operators are equal).  The plant was wrong, not the guard.
# (b) temp copies (run file, beta file, embedded json, paper), module paths pointed at the copy, test functions
#     called directly:
#   T01 run file: full_exps 7/2 => 9/2                        hash FIRED, agreement FIRED
#   T02 run file: D2_order 10 => 11                           hash FIRED, agreement FIRED
#   T03 A1 and A4 run files swapped                           hash FIRED, agreement FIRED (both cases)
#   T04 A1 and A4 beta files swapped                          hash FIRED, decay FIRED, amplitude FIRED
#   T05 embedded full operator, one leading-coefficient entry +1    annihilation FIRED, indicial FIRED
#   T06 embedded full operator, one c_0 entry +1              annihilation FIRED; indicial DID NOT FIRE (expected:
#        c_0 has t-valuation 0 at x*, far above the lowest v_k - k, so it is not in chi; annihilation is the guard)
#   T07 embedded full operator, one c_12 entry +1             annihilation FIRED, indicial FIRED
#   T08 embedded A1 and A4 entries swapped                    prefix FIRED, indicial FIRED, pieces FIRED, blocks KeyError
#   T09 embedded sequence term L = 17 +1                      prefix FIRED, annihilation FIRED, pieces FIRED
#   T10 embedded sequence term L = 500 +1                     annihilation FIRED, pieces FIRED; prefix DID NOT FIRE
#        (expected: beyond L = 30; the slow tier recomputes all 1001 terms)
#   T11 embedded D2 block operator := D1 block operator       blocks FIRED
#   T12 embedded non-carrier piece operator := a D1 piece's   pieces FIRED
#   T13 paper: caveat sentence deleted                        FIRED
#   T14 paper: 'and not derived' => 'and derived'             FIRED
#   T15 paper: 'with no $9/2$' => 'and $9/2$'                 FIRED
#   T16 paper: 'held at two' => 'held at three'               FIRED
#   T17 Phi5(A1) x (1 + 1e-3), formula check bypassed         FIRED (the fit alone)
#   T18 Phi5(A4) x (1 - 1e-6), formula check bypassed         FIRED (five-term bound)
# Not planted: an error inside the shared engine or guesser (same code as the scored run; the hash test pins it, and
# the T-form annihilation check is the only route here that does not go through the guesser).
