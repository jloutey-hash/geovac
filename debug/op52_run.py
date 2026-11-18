"""Record OP52: the bond-angle operator of f' = NA - Y_mu at the contact x*.

f' is the (AA|BC)-type integral with A's 1s density replaced by rho_A = mu^2 e^{-mu s} / (4 pi s) (potential
V'(s) = (1 - e^{-mu s}) / s, odd terms s^1: -mu^2/2, s^3: -mu^4/24).  Built from record AB's frozen engine
(`aabc_exact.pieces`, kernels NA and Y only; no dY).  Record: docs/predictions/2026-11-18_aa_bc_operator_exponents.md.

CLI (Phase 1 = controls + fire; Phase 2 only after registration, with --registered):
  python debug/op52_run.py controls <p>          toy Legendre controls (known exponents) -> op52_ctrl_toy_<p>.json
  python debug/op52_run.py betatoys              Legendre toys (x*-x)^g, g = 3/2, 5/2, 7/2 -> op52_ctrl_betatoys.json
  python debug/op52_run.py fire                  gate fire test -> op52_ctrl_fire.json
  python debug/op52_run.py freeze                sha256 of frozen files -> op52_frozen.json (at registration)
  python debug/op52_run.py run <case> <p> --registered     f' mod p, operators       -> op52_run_<case>_<p>.json
  python debug/op52_run.py beta <case> <dps> --registered  f' over Q, beta_eff        -> op52_beta_<case>_<dps>.json
  python debug/op52_run.py score --registered              -> op52_score.json
"""
from __future__ import annotations

import copy
import hashlib
import json
import sys
import time
import traceback
from fractions import Fraction
from pathlib import Path

import flint
import mpmath as mp

HERE = Path(__file__).parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import intclass_lib as I         # noqa: E402
import s1_controls as C          # noqa: E402

DATA = HERE / 'data'
PRIMES = (1073741441, 1073741381)
SEEDS = (2026, 7)
LMAX = 1000
CASES = {'A1': dict(D1='9/5', D2='13/5', zB='13/10', zC='1', zA='9/10'),      # AB's A1; f' never formed
         'A4': dict(D1='8/5', D2='5/2', zB='13/10', zC='4/5', zA='7/10')}     # fresh; nothing ever computed
BETA_L = 220
BETA_LS = list(range(120, 241, 20))
TOL = dict(beta=0.10, dps=1e-12, OP5=1e-3)
WANT_FULL = ['3/2', '5/2', '7/2']
WANT_D1 = ['3/2', '5/2']
WANT_D2 = ['5/2']
FROZEN = ['debug/op52_run.py', 'debug/aabc_exact.py', 'debug/na_exact.py', 'debug/intclass_lib.py',
          'debug/s1_controls.py', 'debug/na_modguess.py', 'debug/na_modanalyze.py', 'debug/l52_local.py',
          'debug/l52_predict.py', 'debug/op52_subsets.py', 'debug/op52_na.py', 'debug/op52_pieces.py',
          'debug/data/op52_subsets_A1_1073741441.json', 'debug/data/op52_subsets_A1_1073741381.json',
          'debug/data/op52_subsets_A2_1073741441.json', 'debug/data/op52_subsets_A2_1073741381.json',
          'debug/data/op52_na_1073741441.json', 'debug/data/op52_na_1073741381.json',
          'debug/data/op52_ctrl_toy_1073741441.json', 'debug/data/op52_ctrl_toy_1073741381.json',
          'debug/data/op52_ctrl_betatoys.json', 'debug/data/op52_ctrl_fire.json',
          'debug/data/s1_coeffs_A1_1073741441.json', 'debug/data/s1_coeffs_A1_1073741381.json',
          'debug/data/s1_coeffs_A2_1073741441.json', 'debug/data/s1_coeffs_A2_1073741381.json']


def params(case):
    return tuple(Fraction(CASES[case][k]) for k in ('D1', 'D2', 'zB', 'zC', 'zA'))


def xstar(case):
    D1, D2 = params(case)[:2]
    return (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)


def endpoint_keys(case):
    D1, D2, zB, zC, zA = params(case)
    mu = 2 * zA
    one = lambda e: (Fraction(e), '1', None)                    # noqa: E731
    d1 = [one(-zB * D1), one(-zB * D1 - mu * (D2 - D1)), one(-zB * D1 - mu * (D2 + D1))]
    d2 = [one(zC * D1 - (zB + zC) * D2), one(-zC * D1 - (zB + zC) * D2)]
    return d1, d2


def halves(exps):
    return sorted((e for e in (exps or []) if Fraction(e).denominator == 2), key=Fraction)


def xs_exps(c, p, xs):
    kx = ((-I.to_fp(xs, p)) % p, 1)
    pts = C.points(c, p)
    return sorted(pts[kx]['exps'], key=Fraction) if kx in pts else None


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


# ------------------------------------------------------------------------------------------------ gate
def gate(s):
    """s = {'controls_pass': bool, 'hash_ok': bool, 'errors': [str],
            'cases': {case: {p: {'found': bool, 'keys_ok': bool, 'full_od': [o, d],
                                 'full': [halves seed1, halves seed2], 'D1': halves, 'D2': halves}}},
            'beta': {case: {'b': float, 'toy52': float, 'dps_rel': float}}}
    -> (verdict, clauses, void_flags).  VOID > STOP > BORDERLINE > GO.  A missing case is STOP."""
    void = []
    if not s.get('controls_pass'):
        void.append('control failed')
    if not s.get('hash_ok'):
        void.append('frozen hash mismatch')
    if s.get('errors'):
        void.append('exception')
    per = {}
    for case in CASES:
        cd = s.get('cases', {}).get(case)
        if cd is None or any(p not in cd for p in PRIMES):
            per[case] = None
            continue
        rows = [cd[p] for p in PRIMES]
        if not all(r['found'] for r in rows):
            void.append(f'{case}: minimal operator not found')
            per[case] = None
            continue
        if not all(r['keys_ok'] for r in rows):
            void.append(f'{case}: endpoint keys missing or coincident')
        if any(r['full'][0] != r['full'][1] or r.get('D1_s2', r['D1']) != r['D1'] or r.get('D2_s2', r['D2']) != r['D2']
               for r in rows):
            void.append(f'{case}: seeds disagree')
        if (rows[0]['full_od'] != rows[1]['full_od'] or rows[0]['full'][0] != rows[1]['full'][0]
                or rows[0]['D1'] != rows[1]['D1'] or rows[0]['D2'] != rows[1]['D2']):
            void.append(f'{case}: primes disagree')
        per[case] = rows[0]
    bet = s.get('beta', {})
    for case in CASES:
        if case in bet and bet[case]['dps_rel'] > TOL['dps']:
            void.append(f'{case}: beta precision')
    cl = {}
    complete = all(per.get(c) is not None for c in CASES) and all(c in bet for c in CASES)
    cl['O1_full'] = complete and all(per[c]['full'][0] == WANT_FULL for c in CASES)
    cl['O2_blocks'] = complete and all(per[c]['D1'] == WANT_D1 and per[c]['D2'] == WANT_D2 for c in CASES)
    cl['F1_beta'] = complete and all(abs(bet[c]['b'] - bet[c]['toy52']) <= TOL['beta'] for c in CASES)
    core = complete and all('3/2' in per[c]['full'][0] and '5/2' in per[c]['full'][0] for c in CASES)
    if void:
        v = 'VOID'
    elif not complete or not core or not cl['F1_beta']:
        v = 'STOP'
    elif cl['O1_full'] and cl['O2_blocks']:
        v = 'GO'
    else:
        v = 'BORDERLINE'
    return v, cl, void


def fire():
    row = lambda: dict(found=True, keys_ok=True, full_od=[9, 60], full=[list(WANT_FULL), list(WANT_FULL)],  # noqa
                       D1=list(WANT_D1), D2=list(WANT_D2), D1_s2=list(WANT_D1), D2_s2=list(WANT_D2))
    good = dict(controls_pass=True, hash_ok=True, errors=[],
                cases={c: {p: row() for p in PRIMES} for c in CASES},
                beta={'A1': dict(b=-3.08, toy52=-3.10, dps_rel=1e-30), 'A4': dict(b=-3.07, toy52=-3.09, dps_rel=1e-30)})
    P1, P2 = PRIMES
    plants = [('baseline', good, 'GO')]

    def plant(name, fn, want):
        x = copy.deepcopy(good)
        fn(x)
        plants.append((name, x, want))

    def both(case, key, val):
        def f(x):
            for p in PRIMES:
                x['cases'][case][p][key] = copy.deepcopy(val)
                x['cases'][case][p][key + '_s2'] = copy.deepcopy(val)
        return f

    def fullset(case, val):
        def f(x):
            for p in PRIMES:
                x['cases'][case][p]['full'] = [list(val), list(val)]
        return f
    # VOID branches
    plant('control failed', lambda x: x.__setitem__('controls_pass', False), 'VOID')
    plant('hash mismatch', lambda x: x.__setitem__('hash_ok', False), 'VOID')
    plant('exception', lambda x: x.__setitem__('errors', ['Traceback']), 'VOID')
    plant('operator not found, one prime', lambda x: x['cases']['A1'][P2].__setitem__('found', False), 'VOID')
    plant('endpoint keys missing', lambda x: x['cases']['A4'][P1].__setitem__('keys_ok', False), 'VOID')
    plant('seeds disagree', lambda x: x['cases']['A1'][P1].__setitem__('full', [WANT_FULL, ['3/2', '5/2']]), 'VOID')
    plant('primes disagree: full set', lambda x: x['cases']['A4'][P2].__setitem__('full', [['3/2', '5/2', '9/2']] * 2), 'VOID')
    plant('primes disagree: order', lambda x: x['cases']['A1'][P2].__setitem__('full_od', [10, 60]), 'VOID')
    plant('primes disagree: D1 block', lambda x: x['cases']['A1'][P1].__setitem__('D1', ['3/2']), 'VOID')
    plant('primes disagree: D2 block', lambda x: x['cases']['A4'][P1].__setitem__('D2', ['7/2']), 'VOID')
    plant('block seeds disagree (D1)', lambda x: x['cases']['A4'][P1].__setitem__('D1_s2', ['3/2']), 'VOID')
    plant('block seeds disagree (D2)', lambda x: x['cases']['A1'][P2].__setitem__('D2_s2', ['7/2']), 'VOID')
    plant('rank-2 control fails (controls_pass False)', lambda x: x.__setitem__('controls_pass', False), 'VOID')
    plant('beta precision', lambda x: x['beta']['A1'].__setitem__('dps_rel', 1e-6), 'VOID')
    plant('void beats stop', lambda x: (x.__setitem__('hash_ok', False), fullset('A1', ['3/2', '7/2', '9/2'])(x)), 'VOID')
    # STOP branches
    plant('5/2 absent (f-like set)', fullset('A1', ['3/2', '7/2', '9/2']), 'STOP')
    plant('3/2 absent', fullset('A4', ['5/2', '7/2']), 'STOP')
    plant('5/2 present only at A4', fullset('A1', ['3/2', '7/2']), 'STOP')
    plant('no half exponent', fullset('A4', []), 'STOP')
    plant('F1: beta at 3/2 toy', lambda x: x['beta']['A1'].__setitem__('b', -2.05), 'STOP')
    plant('F1: beta at 7/2 toy', lambda x: x['beta']['A4'].__setitem__('b', -4.10), 'STOP')
    plant('F1: beta 0.15 off', lambda x: x['beta']['A4'].__setitem__('b', -3.24), 'STOP')
    plant('case missing', lambda x: x['cases'].pop('A4'), 'STOP')
    plant('prime missing', lambda x: x['cases']['A1'].pop(P2), 'STOP')
    plant('beta missing', lambda x: x['beta'].pop('A1'), 'STOP')
    # BORDERLINE branches
    plant('9/2 instead of 7/2 (collision resolves late)', fullset('A1', ['3/2', '5/2', '9/2']), 'BORDERLINE')
    plant('extra 9/2 (four branches)', fullset('A4', ['3/2', '5/2', '7/2', '9/2']), 'BORDERLINE')
    plant('only {3/2, 5/2} (D1, D2 branches coincide)', fullset('A1', ['3/2', '5/2']), 'BORDERLINE')
    plant('D1 block one-dimensional', both('A1', 'D1', ['3/2']), 'BORDERLINE')
    plant('D1 block three branches', both('A4', 'D1', ['3/2', '5/2', '7/2']), 'BORDERLINE')
    plant('D2 block at 7/2', both('A4', 'D2', ['7/2']), 'BORDERLINE')
    plant('D2 block two branches', both('A1', 'D2', ['5/2', '7/2']), 'BORDERLINE')
    plant('OP5 amplitude off, OP4 non-monotone (reported only)',
          lambda x: x.__setitem__('reported', {'A1': dict(OP4_monotone=False, OP5=dict(r0='1.2', holds=False))}), 'GO')
    plant('beta 0.09 off (inside tol)', lambda x: x['beta']['A1'].__setitem__('b', -3.19), 'GO')
    res = []
    for name, x, want in plants:
        got = gate(x)[0]
        res.append(dict(plant=name, got=got, want=want, ok=got == want))
        print(f'{name:48s} got {got:10s} want {want:10s} {"ok" if got == want else "FAIL"}')
    ok = all(r['ok'] for r in res)
    print('FIRE TEST', 'PASS' if ok else 'FAIL', f'{sum(r["ok"] for r in res)}/{len(res)}')
    (DATA / 'op52_ctrl_fire.json').write_text(json.dumps(dict(all_pass=ok, plants=res), indent=1))
    return ok


# ------------------------------------------------------------------------------------------------ controls
def toy_seq(k, t, p, n, mult_x=0):
    """Legendre coefficients (mod p) of (1 - 2 x t + t^2)^(k - 1/2) * x^mult_x, i.e. of (x_t - x)^(k - 1/2) up to
    a constant and the polynomial x^mult_x; x_t = (1 + t^2)/(2t).  Exact recurrence from sum t^L P_L."""
    T = I.to_fp(t, p)
    N = n + k + mult_x + 2
    c = [pow(T, L, p) for L in range(N + k + mult_x + 2)]

    def inv(a):
        return pow(a % p, p - 2, p)

    def xmul(c):        # coefficients of x * sum c_m P_m
        out = []
        for L in range(len(c) - 1):
            v = c[L + 1] * (L + 1) * inv(2 * L + 3)
            if L >= 1:
                v += c[L - 1] * L * inv(2 * L - 1)
            out.append(v % p)
        return out
    for _ in range(k):
        xc = xmul(c)
        c = [((1 + T * T) * c[L] - 2 * T * xc[L]) % p for L in range(len(xc))]
    for _ in range(mult_x):
        c = xmul(c)
    return c[:n]


def int_seq(poly, t, p, n):
    """Legendre coefficients (mod p) of h_c(x) = int_0^t c(s) (1 - 2 x s + s^2)^(-1/2) ds, c(s) = sum_m poly[m] s^m:
    q_L = sum_m poly[m] t^(L+m+1) / (L+m+1).  Near x_t the half-integer branch of h_c is (x_t - x)^(k + 1/2) times a
    unit, k = ord_{s=t} c (endpoint expansion of the integral); a combination of pieces has the valuation of its
    combined integrand."""
    T = I.to_fp(t, p)
    poly = [I.to_fp(Fraction(c), p) for c in poly]
    return [sum(c * pow(T, L + m + 1, p) * pow(L + m + 1, p - 2, p) for m, c in enumerate(poly)) % p for L in range(n)]


def stp(k, t, extra=(1,)):
    """coefficients of (s - t)^k * extra(s), extra given low-to-high"""
    c = [Fraction(v) for v in extra]
    for _ in range(k):
        c = [(-t) * c[0]] + [c[i - 1] - t * c[i] for i in range(1, len(c))] + [c[-1]]
    return c


def rank2_operator(k):
    """Order-2 Fuchsian operator a2 y'' + a1 y' + a0 y (integer polynomial coefficients, sympy, exact) whose solutions
    are y1 = (x_t - x)^{3/2} and y2 = (x_t - x)^k (x_u - x)^{1/2}, x_t = 5/3, x_u = -3.  Known exponents at x_t:
    {3/2, k}; monodromy at x_u flips y2 only, so a generic solution spans both (rank 2).  Verified symbolically."""
    import sympy as sp
    x = sp.symbols('x')
    xt, xu = sp.Rational(5, 3), sp.Integer(-3)
    y1 = (xt - x) ** sp.Rational(3, 2)
    y2 = (xt - x) ** sp.Rational(k) * (xu - x) ** sp.Rational(1, 2)
    W = y1 * y2.diff(x) - y2 * y1.diff(x)
    q = sp.simplify((y1.diff(x) * y2.diff(x, 2) - y2.diff(x) * y1.diff(x, 2)) / W)
    pp = sp.simplify(-W.diff(x) / W)
    co = [q, pp, sp.Integer(1)]
    den = sp.lcm([sp.denom(sp.together(v)) for v in co])
    P = [sp.expand(sp.cancel(v * den)) for v in co]
    for yy in (y1, y2):
        assert sp.simplify(P[2] * yy.diff(x, 2) + P[1] * yy.diff(x) + P[0] * yy) == 0
    return [[Fraction(int(c.p), int(c.q)) for c in reversed(sp.Poly(v, x).all_coeffs())] for v in P]   # low-to-high


def rank2_seq(k, p, n, seed):
    """Legendre coefficients (mod p) of a solution of rank2_operator(k), from random initial values: the
    operator (1-x^2)^2 L acting on P_m, via (1-x^2) P_m' = m(m+1)/(2m+1) (P_{m-1} - P_{m+1}),
    (1-x^2)^2 P_m'' = 2x (1-x^2) P_m' - m(m+1)(1-x^2) P_m, and x P_m = ((m+1) P_{m+1} + m P_{m-1})/(2m+1),
    gives a banded linear recurrence on the coefficients, solved forward."""
    import random
    a0, a1, a2 = ([I.to_fp(c, p) for c in poly] for poly in rank2_operator(k))
    inv = lambda v: pow(v % p, p - 2, p)                  # noqa: E731

    def xmul(v):
        out = {}
        for m, c in v.items():
            out[m + 1] = (out.get(m + 1, 0) + c * (m + 1) * inv(2 * m + 1)) % p
            if m >= 1:
                out[m - 1] = (out.get(m - 1, 0) + c * m * inv(2 * m + 1)) % p
        return out

    def add(u, v, s=1):
        out = dict(u)
        for m, c in v.items():
            out[m] = (out.get(m, 0) + s * c) % p
        return out

    def pmul(poly, v):                                    # Horner
        out = {}
        for c in reversed(poly):
            out = add(xmul(out), {m: c * w % p for m, w in v.items()})
        return out

    def image(m):
        Pm = {m: 1}
        r = m * (m + 1) % p
        v1 = {m + 1: (-r * inv(2 * m + 1)) % p}
        if m >= 1:
            v1[m - 1] = r * inv(2 * m + 1) % p
        omx2 = lambda v: add(v, xmul(xmul(v)), -1)        # noqa: E731  (1 - x^2) v
        t2 = add({k_: 2 * c % p for k_, c in xmul(v1).items()}, {k_: r * c % p for k_, c in omx2(Pm).items()}, -1)
        return add(add(pmul(a2, t2), pmul(a1, omx2(v1))), pmul(a0, omx2(omx2(Pm))))
    N = n + 40
    img = [image(m) for m in range(N)]
    w = max(max(abs(r - m) for r in img[m]) for m in range(N))
    rows = {}                                             # rows[r] = {m: coefficient of c_m}
    for m in range(N):
        for r, c in img[m].items():
            rows.setdefault(r, {})[m] = c
    rnd = random.Random(seed)
    c = [rnd.randrange(p) for _ in range(w)]
    for r in range(0, N - w):
        row = rows.get(r, {})
        top = r + w
        lead = row.get(top, 0)
        if lead == 0:
            raise ValueError(f'zero leading coefficient at row {r}')
        s = sum(v * c[m] for m, v in row.items() if m < top) % p
        c.append((-s) * inv(lead) % p)
    assert all(sum(v * c[m] for m, v in rows.get(r, {}).items() if m < len(c)) % p == 0 for r in range(0, N - w))
    return c[:n]


def controls(p, n=500):
    t = Fraction(1, 3)
    xt = (1 + t * t) / (2 * t)
    g = lambda k, m=0: toy_seq(k, t, p, n, m)              # noqa: E731
    h = lambda k, ex=(1,): int_seq(stp(k, t, ex), t, p, n)  # noqa: E731
    sets = {  # rank one (rational multiples of one algebraic function): exactly one exponent, the lowest
            'alg_g3': ({'g3': g(3)}, ['5/2']),
            'alg_g4': ({'g4': g(4)}, ['7/2']),
            'alg_g2_g4_g5': ({'g2': g(2), 'g4': g(4), 'g5': g(5)}, ['3/2']),
            # integrals of s^m g0: they reduce (integration by parts in s) to two log-type periods plus
            # rational x sqrt(1 - 2xt + t^2), so their x_t half-branches are all rational multiples of ONE function:
            # rank one again, exactly the lowest exponent (measured first with a wrong rank-3 expectation; disclosed)
            'int_k1_k3_k4': ({'a': h(1), 'b': h(3), 'c': h(4)}, ['3/2']),
            'int_k1_k2_k3': ({'a': h(1), 'b': h(2), 'c': h(3)}, ['3/2']),
            'int_k1': ({'a': h(1)}, ['3/2'])}
    # rank two (independent of every engine): solutions of a known order-2 Fuchsian operator, random initial values;
    # exponents at x_t = 5/3: {3/2, 5/2} (variant: {3/2, 7/2}, 5/2 absent)
    sets['rank2_52'] = ({'s': rank2_seq(Fraction(5, 2), p, n, 11)}, ['3/2', '5/2'])
    sets['rank2_72'] = ({'s': rank2_seq(Fraction(7, 2), p, n, 12)}, ['3/2', '7/2'])
    res = dict(prime=p, xt=str(xt), sets={})
    for name, (sq, want) in sets.items():
        got = []
        for sd in SEEDS:
            op = C.generic_operator(C.generic_seq(sq, p, sd), p)
            got.append(halves(xs_exps(op['c'], p, xt)) if op else None)
        res['sets'][name] = dict(want=want, got=got, ok=all(gg == want for gg in got))
        print(name, got, want, flush=True)
    res['pass'] = all(v['ok'] for v in res['sets'].values())
    (DATA / f'op52_ctrl_toy_{p}.json').write_text(json.dumps(res, indent=1))
    return res


def betatoys():
    import l52_predict as LP
    mp.mp.dps = 40
    out = {}
    for case in CASES:
        xs = xstar(case)
        xsm = mp.mpf(xs.numerator) / xs.denominator
        out[case] = dict(xstar=str(xs), **{f'toy_{g}': float(LP.beta_toy(mp.mpf(Fraction(g).numerator) /
                                                                         Fraction(g).denominator, xsm, BETA_L))
                                            for g in ('3/2', '5/2', '7/2')})
        print(case, out[case], flush=True)
    (DATA / 'op52_ctrl_betatoys.json').write_text(json.dumps(out, indent=1))


def freeze():
    out = {f: sha(f) for f in FROZEN if f != 'debug/op52_run.py'}
    out['debug/op52_run.py'] = sha('debug/op52_run.py')
    (DATA / 'op52_frozen.json').write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


# ------------------------------------------------------------------------------------------------ Phase 2
def fprime_coeffs(case, K):
    import aabc_exact as X
    D1, D2, zB, zC, zA = params(case)
    pc = X.pieces(case, LMAX if K.p else 240, K, which=('NA', 'Y'), params_override=(D1, D2, zB, zC, zA))
    import na_exact as E
    out = []
    for qN, qY in zip(pc['NA'], pc['Y']):
        r = dict(qN)
        for k, v in qY.items():
            E.add(r, k, -v)
        out.append({k: v for k, v in r.items() if v != 0})
    return out


def run(case, p):
    import aabc_exact as X
    import na_exact as E
    X.CASES.setdefault(case, CASES[case])          # pieces() reads params via params_override only
    t0 = time.time()
    rec = dict(case=case, prime=p, xstar=str(xstar(case)))
    try:
        f = fprime_coeffs(case, E.Field(p))
        keys = sorted(set(k for q in f for k in q), key=str)
        seqs = {k: [int(q.get(k, flint.nmod(0, p))) for q in f] for k in keys}
        rec['n_const'] = len(keys)
        rec['t_coeff'] = round(time.time() - t0, 1)
        d1, d2 = endpoint_keys(case)
        rec['keys_ok'] = all(k in seqs for k in d1 + d2) and len(set(d1 + d2)) == 5
        xs = xstar(case)
        full = [C.generic_operator(C.generic_seq(seqs, p, sd), p) for sd in SEEDS]
        rec['found'] = all(o is not None for o in full)
        if rec['found']:
            rec['full_od'] = [full[0]['order'], full[0]['degree']]
            rec['full_od_seed2'] = [full[1]['order'], full[1]['degree']]
            rec['full_exps'] = xs_exps(full[0]['c'], p, xs)
            rec['full'] = [halves(xs_exps(o['c'], p, xs)) for o in full]
            for name, ks in (('D1', d1), ('D2', d2)):
                ks = [k for k in ks if k in seqs]
                for sd, tag in ((SEEDS[0], ''), (SEEDS[1], '_s2')):
                    op = C.generic_operator(C.generic_seq({k: seqs[k] for k in ks}, p, sd), p) if ks else None
                    rec[name + tag] = halves(xs_exps(op['c'], p, xs)) if op else None
                    rec[name + tag + '_order'] = op and op['order']
            # O3 (reported): which pieces carry x*
            car = {}
            for k, sq in seqs.items():
                po = C.piece_operator(sq, p)
                e = xs_exps(po['c'], p, xs) if po else 'not found'
                if e is not None:
                    car[repr(k)] = dict(order=po and po['order'], halves=halves(e) if e != 'not found' else e)
            rec['O3_carriers'] = car
            rec['O3_exactly_endpoint_five'] = set(car) == {repr(k) for k in d1 + d2}
    except Exception:
        rec['error'] = traceback.format_exc()
        rec['found'] = False
    rec['time'] = round(time.time() - t0, 1)
    (DATA / f'op52_run_{case}_{p}.json').write_text(json.dumps(rec, indent=1, default=str))
    print(json.dumps({k: v for k, v in rec.items() if k != 'O3_carriers'}, indent=1, default=str))


def beta(case, dps):
    import aabc_exact as X
    import na_exact as E
    X.CASES.setdefault(case, CASES[case])
    t0 = time.time()
    rec = dict(case=case, dps=dps)
    try:
        f = fprime_coeffs(case, E.Field(None))
        mp.mp.dps = dps
        keys = set(k for q in f for k in q)
        cv = {k: E.value(k, mp) for k in keys}
        D1, D2, zB, zC, zA = (mp.mpf(v.numerator) / v.denominator for v in params(case))
        tw = 2 * mp.sqrt(zB ** 3 * zC ** 3)              # L52's normalization of F_L
        F = [tw * sum((mp.mpf(int(v.p)) / int(v.q) * cv[k] for k, v in q.items()), mp.mpf(0)) for q in f]
        xs = xstar(case)
        xsm = mp.mpf(xs.numerator) / xs.denominator
        rho = 1 / (xsm + mp.sqrt(xsm * xsm - 1))
        G = [F[L] / rho ** L for L in range(len(F))]
        rec['beta'] = {L: float(mp.log(abs(G[L + 1] / G[L - 1])) / mp.log(mp.mpf(L + 1) / (L - 1)))
                       for L in (40, 80, 120, 160, 200, BETA_L)}
        rec['F'] = {str(L): mp.nstr(F[L], 40) for L in BETA_LS}
        bl = [rec['beta'][L] for L in sorted(rec['beta'])]
        rec['OP4_beta_monotone_toward_-3'] = all(abs(b + 3) > abs(c + 3) for b, c in zip(bl, bl[1:]))
        # OP5 (reported): R_L = F_L / (Phi5 T^{5/2}_L), Phi5 = L52 segment amplitude with p = q = 1, v_1 = -mu^2/2
        import l52_local as LL
        mu = 2 * zA
        phi5 = LL.amplitude(1, 1, -mu ** 2 / 2, D1, D2, zB, zC)
        R = [F[L] / (phi5 * LL.T(mp.mpf(5) / 2, xsm, L)) for L in BETA_LS]
        r0, c1 = LL.ratio_extrap(BETA_LS, R)
        rec['OP5'] = dict(Phi5=mp.nstr(phi5, 20), R_240=mp.nstr(R[-1], 12), r0=mp.nstr(r0, 12), c1=mp.nstr(c1, 8),
                          holds=bool(abs(r0 - 1) <= TOL['OP5']))
        rec['signs_L_ge_20'] = sorted(set(int(mp.sign(v)) for v in F[20:]))
    except Exception:
        rec['error'] = traceback.format_exc()
    rec['time'] = round(time.time() - t0, 1)
    (DATA / f'op52_beta_{case}_{dps}.json').write_text(json.dumps(rec, indent=1, default=str))
    print(json.dumps({k: v for k, v in rec.items() if k != 'F'}, indent=1, default=str))


def score():
    s = dict(controls_pass=True, hash_ok=True, errors=[], cases={}, beta={})
    try:
        frozen = json.loads((DATA / 'op52_frozen.json').read_text())
        s['hash_ok'] = all(sha(f) == h for f, h in frozen.items())
        ctl = [json.loads((DATA / f'op52_ctrl_toy_{p}.json').read_text())['pass'] for p in PRIMES]
        ctl.append(json.loads((DATA / 'op52_ctrl_fire.json').read_text())['all_pass'])
        for p in PRIMES:                                 # the measured (AA|BC) and NA controls, re-read
            for c in ('A1', 'A2'):
                d = json.loads((DATA / f'op52_subsets_{c}_{p}.json').read_text())['subsets']
                ctl += [d['trio_D1']['halves'] == ['3/2', '7/2'], d['pair_D2']['halves'] == ['7/2'],
                        d['all5']['halves'] == ['3/2', '7/2', '9/2']]
            d = json.loads((DATA / f'op52_na_{p}.json').read_text())['subsets']
            ctl += [d['D1']['halves'] == ['3/2'], d['pair']['halves'] == ['3/2'],
                    d['D1_plus_pair']['halves'] == ['3/2', '5/2']]
        s['controls_pass'] = all(ctl)
        toys = json.loads((DATA / 'op52_ctrl_betatoys.json').read_text())
        for case in CASES:
            s['cases'][case] = {}
            for p in PRIMES:
                fn = DATA / f'op52_run_{case}_{p}.json'
                if not fn.exists():
                    continue
                r = json.loads(fn.read_text())
                if 'error' in r:
                    s['errors'].append(f'{case} {p}')
                    continue
                if r['found'] and r['full_od'] != r['full_od_seed2']:
                    r['full'] = [r['full'][0], ['seed order differs']]
                s['cases'][case][p] = dict(found=r['found'], keys_ok=r.get('keys_ok', False),
                                           full_od=r.get('full_od'), full=r.get('full', [None, None]),
                                           D1=r.get('D1'), D2=r.get('D2'), D1_s2=r.get('D1_s2'),
                                           D2_s2=r.get('D2_s2'))
            lo, hi = (DATA / f'op52_beta_{case}_2000.json', DATA / f'op52_beta_{case}_2600.json')
            if lo.exists() and hi.exists():
                a, b = json.loads(lo.read_text()), json.loads(hi.read_text())
                if 'error' in a or 'error' in b:
                    s['errors'].append(f'{case} beta')
                    continue
                rel = max(float(abs(mp.mpf(a['F'][str(L)]) / mp.mpf(b['F'][str(L)]) - 1)) for L in BETA_LS)
                s['beta'][case] = dict(b=a['beta'][str(BETA_L)], toy52=toys[case]['toy_5/2'], dps_rel=rel)
                s.setdefault('reported', {})[case] = dict(OP4_monotone=a.get('OP4_beta_monotone_toward_-3'),
                                                          OP5=a.get('OP5'))
    except Exception:
        s['errors'].append(traceback.format_exc())
    v, cl, void = gate(s)
    out = dict(verdict=v, clauses=cl, void_flags=void, inputs=s)
    (DATA / 'op52_score.json').write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps(out, indent=1, default=str))


if __name__ == '__main__':
    a = sys.argv[1:]
    cmd = a[0]
    if cmd in ('run', 'beta', 'score') and '--registered' not in a:
        sys.exit('Phase 2 refuses to run without --registered (record OP52 must be committed first).')
    a = [x for x in a if x != '--registered']
    if cmd == 'controls':
        controls(int(a[1]))
    elif cmd == 'betatoys':
        betatoys()
    elif cmd == 'fire':
        sys.exit(0 if fire() else 1)
    elif cmd == 'freeze':
        freeze()
    elif cmd == 'run':
        run(a[1], int(a[2]))
    elif cmd == 'beta':
        beta(a[1], int(a[2]))
    elif cmd == 'score':
        score()
