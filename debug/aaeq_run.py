"""Record AE (docs/predictions/2026-11-15_aa_bc_equal_arms.md): (AA|BC) at equal arms D1 = D2.
Registers cases E1, E2 in THIS process only (aabc_exact.CASES is shared by aapc_exact).

    python debug/aaeq_run.py ap0             : AE0 (i) D1 < D2 output unchanged; (ii) per-channel direct-projection
                                               quadrature at E1, E2; (iii) continuity at D1 = D - 1e-40
    python debug/aaeq_run.py p1p4 [Lmax dps] : P1 values, P3 exponents, P4 amplitude ratios
    python debug/aaeq_run.py p2 <prime>      : P2 constant sets, L <= 1000
"""
import importlib.util
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import aabc_exact as X  # noqa: E402
import aapc_exact as P  # noqa: E402
import na_exact as E  # noqa: E402

DATA = HERE / 'data'
X.CASES['E1'] = dict(D1='2', D2='2', zB='13/10', zC='1', zA='9/10')
X.CASES['E2'] = dict(D1='12/5', D2='12/5', zB='1', zC='6/5', zA='4/5')
EQ = ('E1', 'E2')
INTS = {'S': 'AB', 'PC': 'PC', 'PD': 'PD'}
PTS = {'E1': ('1/2', '0'), 'E2': ('-1/3', '4/5')}
XKEY = {'1/2': 'x=0.5', '0': 'x=0', '-1/3': 'x=-0.333333', '4/5': 'x=0.8'}


def mfq(q):
    q = Fraction(q)
    return mp.mpf(q.numerator) / q.denominator


def numeric(f, case, integral):
    tw = P.twopiN(case, integral)
    keys = set(k for q in f for k in q)
    cv = {k: E.value(k, mp) for k in keys}
    return [tw * sum(mp.mpf(int(v.p)) / int(v.q) * cv[k] for k, v in q.items()) for q in f]


# ---------------------------------------------------------------- AE0
def legendre(L, t):
    p0, p1 = mp.mpf(1), t
    if L == 0:
        return p0
    for n in range(1, L):
        p0, p1 = p1, ((2 * n + 1) * t * p1 - n * p0) / (n + 1)
    return p1


def quad_FL(case, name, L, dps=22):
    """Direct projection; densities as written; split at r = D (A's and C's kinks coincide)."""
    with mp.workdps(dps):
        D1, D2, zB, zC, zA = (mfq(v) for v in P.params(case))
        D = D1
        NB = mp.sqrt(zB ** 3 / mp.pi)
        NC = {'S': mp.sqrt(zC ** 3 / mp.pi), 'PC': mp.sqrt(zC ** 5 / mp.pi),
              'PD': mp.sqrt(5 * (2 * zC) ** 7 / (4 * mp.pi * 720))}[name]
        ang = {'S': lambda s, t, r: 1, 'PC': lambda s, t, r: r * t - D,
               'PD': lambda s, t, r: (3 * (r * t - D) ** 2 - s * s) / 2}[name]
        nodes = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(5, mp.mp.prec)

        def mult(func, r):
            a, b = abs(r - D), r + D
            h, c = (b - a) / 2, (b + a) / 2
            tot = mp.mpf(0)
            for u, w in nodes:
                s = c + h * u
                t = (r * r + D * D - s * s) / (2 * r * D)
                tot += w * func(s, t) * s * legendre(L, t)
            return (2 * L + 1) / mp.mpf(2) * tot * h / (r * D)

        def VA(s):
            return 1 / s - (zA + 1 / s) * mp.exp(-2 * zA * s)

        def integrand(r):
            rho = NB * NC * mp.exp(-zB * r) * mult(lambda s, t: mp.exp(-zC * s) * ang(s, t, r), r)
            return 4 * mp.pi / (2 * L + 1) * r * r * rho * mult(lambda s, t: VA(s), r)
        return mp.quad(integrand, [0, D, D + 4, mp.inf])


def ap0():
    ok = True
    # (i) D1 < D2 unchanged
    spec = importlib.util.spec_from_file_location('engine_v6730', str(Path(sys.argv[2]) if len(sys.argv) > 2 else
                                                  HERE / 'aapc_exact.py'))
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    for pr in (None, (1 << 61) - 1):
        K = E.Field(pr)
        for case in ('A1', 'A2'):
            a = P.coefficients(case, 30, K, 'AB')
            b = X.combine(X.pieces(case, 30, K), P.params(case)[4], K)
            ok &= a == b
            print(case, 'Q' if pr is None else 'p', '(i) AB reproduced:', a == b, flush=True)
    Kp = E.Field((1 << 61) - 1)
    for integral in ('PC', 'PD'):
        same = P.coefficients('A1', 20, Kp, integral) == old.coefficients('A1', 20, Kp, integral)
        ok &= same
        print('A1', integral, '(i) identical to the committed v6.7.30 engine:', same, flush=True)
    # (ii) quadrature
    K = E.Field(None)
    mp.mp.dps = 80
    for case in EQ:
        for name, integral in INTS.items():
            f = P.coefficients(case, 8, K, integral)
            ex = numeric(f, case, integral)
            for L in (0, 3, 8):
                q = quad_FL(case, name, L)
                with mp.workdps(40):
                    rel = abs(q / ex[L] - 1)
                good = rel < mp.mpf('1e-18')
                ok &= bool(good)
                print(case, name, 'L', L, '(ii) exact', mp.nstr(ex[L], 18), 'quad', mp.nstr(q, 18), 'rel',
                      mp.nstr(rel, 3), good, flush=True)
    # (iii) continuity at D1 = D - 1e-40 (three-region path).  Evaluated at 400 digits: the first run at 80 lost
    # digits to the ~1700-digit rational coefficients (kept in data/aaeq_ap0_first_run_80digits.log).
    mp.mp.dps = 400
    eps = Fraction(1, 10 ** 40)
    for case in EQ:
        D1, D2, zB, zC, zA = P.params(case)
        for name, integral in INTS.items():
            f0 = P.coefficients(case, 20, K, integral)
            f1 = P.coefficients(case, 20, K, integral, params_override=(D1 - eps, D2, zB, zC, zA))
            a0, a1 = numeric(f0, case, integral), numeric(f1, case, integral)
            worst = max(abs(a1[L] / a0[L] - 1) for L in range(21))
            good = worst < mp.mpf('1e-35')
            ok &= bool(good)
            print(case, name, '(iii) max rel |F_L(D - 1e-40) / F_L(D) - 1|, L <= 20:', mp.nstr(worst, 3), good,
                  flush=True)
    print('AE0 self-checks pass:', ok)


# ---------------------------------------------------------------- P1, P3, P4
def p1p4(Lmax=240, dps=2000):
    ctrl = json.loads((DATA / 'aaeq_controls.json').read_text())
    out = {}
    K = E.Field(None)
    for case in EQ:
        for name, integral in INTS.items():
            t0 = time.time()
            f = P.coefficients(case, Lmax, K, integral)
            sec = round(time.time() - t0)
            mp.mp.dps = dps
            F = numeric(f, case, integral)
            vals, diffs = {}, {}
            for xs in PTS[case]:
                x = mfq(xs)
                val = sum(F[L] * legendre(L, x) for L in range(Lmax + 1))
                ref = ctrl['referee'][f'{case}_{XKEY[xs]}'][name]['34']
                r1 = ctrl['referee'][f'{case}_{XKEY[xs]}']['R1'][name]
                vals[xs] = mp.nstr(val, 30)
                diffs[xs] = dict(R2=float(val - mp.mpf(repr(ref))), R1=float(val - mp.mpf(repr(r1))))
            beta = {L: float(mp.log(abs(F[L + 1] / F[L - 1])) / mp.log(mp.mpf(L + 1) / (L - 1)))
                    for L in list(range(20, Lmax, 20)) + [Lmax - 20]}
            pred = ctrl['amplitude'][case][name]['pred_F']
            ratio = {L: float(F[int(L)] / mp.mpf(v)) for L, v in pred.items()}
            band = (-8.4, -7.6) if name == 'S' else (-10.4, -9.6)
            rec = dict(values=vals, diffs=diffs, P1=all(abs(d['R2']) <= 1e-10 for d in diffs.values()), beta=beta,
                       P3=band[0] <= beta[Lmax - 20] <= band[1], P4_ratio=ratio,
                       P4=0.9 <= ratio['220'] <= 1.1, signs_first12=[int(mp.sign(v)) for v in F[:12]],
                       signs_L200_205=[int(mp.sign(v)) for v in F[200:206]],
                       tail=mp.nstr(abs(F[Lmax]), 5), engine_s=sec, dps=dps, Lmax=Lmax)
            out[f'{case}_{name}'] = rec
            print(case, name, json.dumps(rec, indent=1), flush=True)
            (DATA / f'aaeq_p1p4_{dps}.json').write_text(json.dumps(out, indent=1))


def p2(pr, Lmax=1000):
    out = {}
    K = E.Field(pr)
    for case in EQ:
        for name, integral in INTS.items():
            t0 = time.time()
            f = P.coefficients(case, Lmax, K, integral)
            sets = [frozenset(q) for q in f]
            S0 = sets[2]
            changes = [L for L in range(3, len(sets)) if sets[L] != S0]
            rec = dict(seconds=round(time.time() - t0), n_constants_L2=len(S0), L_where_set_changes=changes[:20],
                       P2a=not changes, gamma_present=[repr(k) for k in S0 if k[1] == 'g'],
                       sizes_L0_1=[len(sets[0]), len(sets[1])], constants=sorted(map(repr, S0)))
            rec['P2c'] = not rec['gamma_present']
            out[f'{case}_{name}'] = rec
            print(case, name, pr, json.dumps({k: v for k, v in rec.items() if k != 'constants'}), flush=True)
            (DATA / f'aaeq_p2_{pr}.json').write_text(json.dumps(out, indent=1))
    for case in EQ:
        sets = {n: set(out[f'{case}_{n}']['constants']) for n in INTS}
        out[f'{case}_shared'] = dict(S_eq_PC=sets['S'] == sets['PC'], S_eq_PD=sets['S'] == sets['PD'],
                                     PC_eq_PD=sets['PC'] == sets['PD'])
    (DATA / f'aaeq_p2_{pr}.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'ap0':
        ap0()
    elif cmd == 'p1p4':
        p1p4(*(int(v) for v in sys.argv[2:4]))
    elif cmd == 'p2':
        p2(int(sys.argv[2]))
