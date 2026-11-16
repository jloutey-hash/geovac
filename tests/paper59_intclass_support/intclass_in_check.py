"""Record IN (docs/predictions/2026-11-15_na_n2_decided.md): the erratum check for record IC's frozen instrument,
and IN's scoring.  Verdict logic is IC's (intclass_score.operator_sets / component_verdict, intclass_lib.peel),
called unchanged; this file only adds the erratum check on each STUCK residual.

  python debug/intclass_in_check.py control <p>      -> erratum-check controls (no NA data)
  python debug/intclass_in_check.py run <N1|N2> <p>  -> debug/data/intclass_in_<case>_<p>.json
  python debug/intclass_in_check.py score            -> debug/data/intclass_in_score.json

Erratum check (IC's erratum + addendum): the hyperexponential search assigns one exponent per F_p-irreducible factor
of the leading coefficient, so it can miss a solution only at a factor of degree >= 2 whose exponents are not all
integers in F_p.  A STUCK residual PASSES the check iff every such factor has all indicial exponents classified
'int' and none unaccounted (exponents outside F_p).
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import flint

sys.path.insert(0, str(Path(__file__).parent))
import intclass_lib as I            # noqa: E402
import intclass_score as S          # noqa: E402

DATA = Path(__file__).parent / 'data'
NEW_PRIMES = (1073741717, 1073741689)
IC_PRIMES = S.PRIMES


def erratum_check(op, p):
    """PASS iff (a) at EVERY F_p-irreducible factor of the leading coefficient no indicial exponent lies outside F_p
    (unaccounted == 0; added on review: exponents like +-sqrt(2) at a rational point are hyperexponential over
    Qbar(w) and the search does not flag them), and (b) at every factor of degree >= 2 all exponents are integers."""
    bad = []
    for pt in I.finite_points(op, p):
        nonint = [str(v) for (k, v), m in pt['croots'] if k != 'int'] if pt['fdeg'] >= 2 else []
        if nonint or pt['unaccounted'] > 0:
            bad.append(dict(fdeg=pt['fdeg'], nonint=nonint, unaccounted=pt['unaccounted']))
    return dict(passed=not bad, bad=bad)


def analyse_full(c_x, D1, D2, p):
    a0 = I.to_fp((D1 ** 2 + D2 ** 2) / (2 * D1 * D2), p)
    a2 = I.to_fp(Fraction(-1) / (2 * D1 * D2), p)
    cw = I.pullback_quadratic(c_x, a0, a2, p)
    r = I.peel(cw, p, verbose=False)
    out = dict(verdict=r['verdict'], x_order=len(c_x) - 1, x_degree=max(x.degree() for x in c_x),
               w_degree=max(x.degree() for x in cw), steps=len(r['steps']),
               step_flags=sorted({f for s in r['steps'] for f in s['flags']}), time=round(r.get('time', 0), 1))
    if r['verdict'] == 'STUCK':
        cl = I.classify_residual(r, p)
        out['residual'] = dict(order=cl['order'], type=cl['type'], clean=cl['clean'], flags=cl['flags'],
                               lam=cl['lam'])
        out['erratum'] = erratum_check(r['residual'], p)
    return out


def run(case, p, L=1200, Rmax=24):
    import na_exact as NX
    T0 = time.time()
    D1s, D2s, zX, zH = S.NACASES[case]
    D1, D2 = Fraction(D1s), Fraction(D2s)
    Qs = NX.coefficients(D1s, D2s, zX, zH, L, NX.Field(p))
    keys = sorted(set(k for q in Qs for k in q), key=repr)
    res = dict(case=case, prime=p, L=L, components={})
    res['coeff_time'] = round(time.time() - T0, 1)
    print(case, p, 'coefficients to L =', L, len(keys), 'constants', res['coeff_time'], 's', flush=True)
    for k in keys:
        t1 = time.time()
        seq = [int(q.get(k, flint.nmod(0, p))) for q in Qs]
        ops = S.operator_sets(seq, p, Rmax)
        if ops is None:
            rec = dict(verdict='NOT-FOUND')
        else:
            rA = analyse_full(ops['A'], D1, D2, p)
            rB = analyse_full(ops['B'], D1, D2, p) if ops['B'] is not None else None
            rec = dict(R0=ops['R0'], found={str(a): b for a, b in ops['found'].items()}, A=rA, B=rB,
                       verdict=S.component_verdict(rA, rB))
            rec['erratum_pass'] = (rec['verdict'] != 'STUCK') or (rA['erratum']['passed'] and rB['erratum']['passed'])
        rec['time'] = round(time.time() - t1, 1)
        res['components'][repr(k)] = rec
        print(case, p, repr(k), rec['verdict'], rec.get('erratum_pass'), rec.get('found'),
              (rec.get('A') or {}).get('residual'), (rec.get('A') or {}).get('erratum'), rec['time'], 's', flush=True)
        (DATA / f'intclass_in_{case}_{p}.json').write_text(json.dumps(res, indent=1, default=str))
    res['time'] = round(time.time() - T0, 1)
    (DATA / f'intclass_in_{case}_{p}.json').write_text(json.dumps(res, indent=1, default=str))
    print('done', res['time'], 's')


def control(p):
    """Planted: 4(c - w^2) y'' - 4 w y' + y = 0, solutions (w -+ s)^(1/2), s^2 = c (IC's erratum).  TRIANGULAR.
    c a non-square mod p: the search returns STUCK and the erratum check must FAIL it (degree-2 factor, exponent
    1/2).  c a square: the search peels it (TRIANGULAR).  Clean negative: Bessel K0-type operator w y'' + y' - w y,
    whose only finite factor has degree 1: STUCK and the check must PASS it."""
    from sympy import legendre_symbol
    out = {}
    cn = next(c for c in range(2, 200) if legendre_symbol(c, p) == -1)
    cs = next(c for c in range(2, 200) if legendre_symbol(c, p) == 1 and int(c ** 0.5) ** 2 != c)
    for tag, c in (('nonsquare', cn), ('square', cs)):
        op = [flint.nmod_poly([1], p), flint.nmod_poly([0, p - 4], p), flint.nmod_poly([4 * c % p, 0, p - 4], p)]
        r = I.peel(op, p, verbose=False)
        rec = dict(c=c, verdict=r['verdict'])
        if r['verdict'] == 'STUCK':
            rec['erratum'] = erratum_check(r['residual'], p)
        out[tag] = rec
    # no-factor case (degree-1 leading coefficient only): STUCK, passes trivially
    k0 = [flint.nmod_poly([0, p - 1], p), flint.nmod_poly([1], p), flint.nmod_poly([0, 1], p)]
    r = I.peel(k0, p, verbose=False)
    out['k0_nofactor'] = dict(verdict=r['verdict'],
                              erratum=erratum_check(r['residual'], p) if r['verdict'] == 'STUCK' else None)
    # pass side WITH a degree-2 factor (review fix): Legendre nu = -1/2 (elliptic K type), (c - w^2) y'' - 2 w y'
    # - (1/4) y, c a non-square: irreducible, exponents 0, 0 (log) at w = +-sqrt(c).  Must be STUCK and PASS.
    q4 = (p - pow(4, p - 2, p)) % p
    leg = [flint.nmod_poly([q4], p), flint.nmod_poly([0, p - 2], p), flint.nmod_poly([cn % p, 0, p - 1], p)]
    r = I.peel(leg, p, verbose=False)
    out['legendre_K'] = dict(c=cn, verdict=r['verdict'],
                             erratum=erratum_check(r['residual'], p) if r['verdict'] == 'STUCK' else None)
    # exponents outside F_p at a degree-1 point (review fix): w^2 y'' + w y' - c y, solutions w^{+-sqrt(c)},
    # c a non-square: hyperexponential over Qbar(w) (TRIANGULAR), so a STUCK must FAIL the check.
    eul = [flint.nmod_poly([(p - cn) % p], p), flint.nmod_poly([0, 1], p), flint.nmod_poly([0, 0, 1], p)]
    r = I.peel(eul, p, verbose=False)
    out['euler_irrational'] = dict(c=cn, verdict=r['verdict'],
                                   erratum=erratum_check(r['residual'], p) if r['verdict'] == 'STUCK' else None)
    ok = (out['nonsquare']['verdict'] == 'STUCK' and not out['nonsquare']['erratum']['passed']
          and out['square']['verdict'] == 'TRIANGULAR'
          and out['k0_nofactor']['verdict'] == 'STUCK' and out['k0_nofactor']['erratum']['passed']
          and out['legendre_K']['verdict'] == 'STUCK' and out['legendre_K']['erratum']['passed']
          and (out['euler_irrational']['verdict'] == 'TRIANGULAR'
               or not out['euler_irrational']['erratum']['passed']))
    out['all_pass'] = ok
    (DATA / f'intclass_in_controls_{p}.json').write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps(out, default=str))


def score():
    out = {}

    def load(case, p):
        f = DATA / f'intclass_in_{case}_{p}.json'
        return json.loads(f.read_text())['components'] if f.exists() else None

    out['IN-0'] = 'HELD' if all(
        (DATA / f).exists() and json.loads((DATA / f).read_text())['all_pass'] for p in NEW_PRIMES
        for f in (f'intclass_controls_{p}.json', f'intclass_in_controls_{p}.json')) else 'MISSED'

    def case(case):
        per = [load(case, p) for p in NEW_PRIMES]
        if None in per:
            return 'MISSING', None
        keys = set(per[0]) | set(per[1])
        v = {k: tuple((d.get(k, {}).get('verdict', 'NOT-FOUND'), d.get(k, {}).get('erratum_pass')) for d in per)
             for k in keys}
        stuck_ok = [k for k, ((a, ea), (b, eb)) in v.items() if a == b == 'STUCK' and ea and eb]
        if stuck_ok:
            return 'NOT-TRIANGULAR', v
        if all(a == b == 'TRIANGULAR' for (a, _), (b, _) in v.values()):
            return 'TRIANGULAR', v
        return 'INCOMPLETE', v

    vN2, dN2 = case('N2')
    vN1, dN1 = case('N1')
    out['N2'], out['N1'] = vN2, vN1
    out['N2_components'], out['N1_components'] = dN2, dN1
    out['IN-1'] = 'HELD' if vN2 == 'NOT-TRIANGULAR' else ('MISSED' if vN2 == 'TRIANGULAR' else 'UNDECIDED')
    out['IN-2'] = 'HELD' if (vN1 == 'NOT-TRIANGULAR' and dN1 and len(dN1) == 6 and all(
        a == b == 'STUCK' and ea and eb for (a, ea), (b, eb) in dN1.values())) else 'MISSED'
    agree = True
    for case_, d in (('N2', dN2), ('N1', dN1)):
        if not d:
            continue
        old = [json.loads((DATA / f'intclass_NA_{case_}_{p}.json').read_text())['components'] for p in IC_PRIMES]
        for k, ((a, _), (b, _)) in d.items():
            dec = {x for x in (a, b) if x in ('STUCK', 'TRIANGULAR')}
            dec |= {o.get(k, {}).get('verdict') for o in old} & {'STUCK', 'TRIANGULAR'}
            if len(dec) > 1:
                agree = False
    out['IN-3'] = ('HELD' if agree else 'MISSED') if (dN1 and dN2) else 'UNDECIDED'
    if out['IN-0'] != 'HELD':
        gate = 'VOID'
    elif out['IN-1'] == 'HELD' and out['IN-2'] == 'HELD' and out['IN-3'] == 'HELD':
        gate = 'GO'
    elif out['IN-1'] == 'MISSED':
        gate = 'STOP'
    else:
        gate = 'BORDERLINE'
    out['gate'] = gate
    (DATA / 'intclass_in_score.json').write_text(json.dumps(out, indent=1, default=str))
    for k in ('IN-0', 'N1', 'N2', 'IN-1', 'IN-2', 'IN-3', 'gate'):
        print(k, out[k])


if __name__ == '__main__':
    m = sys.argv[1]
    if m == 'control':
        control(int(sys.argv[2]))
    elif m == 'run':
        run(sys.argv[2], int(sys.argv[3]))
    elif m == 'score':
        score()
