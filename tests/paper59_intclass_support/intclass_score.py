"""Record IC (docs/predictions/2026-11-01_integrated_class.md), Phase 2 scorer.  DO NOT RUN before the record is
committed on its own.

  python debug/intclass_score.py run <W1|W2> <prime> [L] [Rmax]   -> debug/data/intclass_<case>_<p>.json (+ .log)
  python debug/intclass_score.py na <N1|N2> <prime> [L]            -> debug/data/intclass_NA_<case>_<p>.json
  python debug/intclass_score.py score                             -> debug/data/intclass_score.json, prints IC-0..6

Per (case, prime): the exact split F_L / (4 zX^3 zH^3) = sum_c q_{L,c} c (debug/ao_exact.py, unchanged since AO3) to
L; for each of the constants c, its component sequence q_{.,c} -> minimal bond-angle operator, built twice
(set A: every null vector at the smallest findable order R0 plus the minimal-degree ones at R0+1, R0+2, gcrd;
set B: R0+1, R0+2, R0+3, gcrd) -> pullback to w = |Y - Z| -> order-one peel (debug/intclass_lib.py) -> classification.
Component verdicts:  TRIANGULAR | STUCK (clean, same residual order from sets A and B) | STUCK-UNCONFIRMED (sets
disagree or a flag) | NOT-FOUND (no operator within L, Rmax).
Case verdict:  NOT-TRIANGULAR if some component is STUCK on both primes;  TRIANGULAR if all components are found and
TRIANGULAR on both primes;  INCOMPLETE otherwise.
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import flint

sys.path.insert(0, str(Path(__file__).parent))
import intclass_lib as I

HERE = Path(__file__).parent
DATA = HERE / 'data'
CASES = {'W1': ('1.7', '2.5', '1.4', '1.1'), 'W2': ('1.6', '2.2', '1.0', '1.2')}      # AO3's scored cases
NACASES = {'N1': ('1.8', '2.6', '1.3', '1.0'), 'N2': ('1.5', '2.4', '1.0', '1.2')}    # NA's cases
PRIMES = (1073741789, 1073741833)            # both = 1 mod 4 (needed for P6's lambda = +-2i to split)
CONTROL_PRIMES = PRIMES


def operator_sets(seq, p, Rmax, margin=30):
    import na_modguess as G
    import na_modanalyze as A
    B = G.Basis([flint.nmod(int(v), p) for v in seq], p)
    R0 = None
    for R in range(1, Rmax + 1):
        Q0, nul = G.min_degree_bisect(B, R, margin, len(seq))
        if Q0 is not None:
            R0 = R
            break
    if R0 is None:
        return None
    found = {}
    for R in range(R0, R0 + 4):
        Q0, nul = G.min_degree_bisect(B, R, margin, len(seq))
        if Q0 is None:
            continue
        nv, keys, (N, nrow) = G.nullity(B, R, Q0, margin)
        found[R] = dict(Q0=Q0, nul=nv,
                        ops=[A.standard_form(keys, [N[i, col] for i in range(len(keys))], p) for col in range(nv)])

    def g_of(lst):
        if not lst:
            return None
        g = lst[0]
        for o in lst[1:]:
            g = I.gcrd(g, o, p)
        return I.content_reduce(g)
    setA = [o for R in (R0, R0 + 1, R0 + 2) if R in found for o in (found[R]['ops'] if R == R0 else found[R]['ops'][:1])]
    setB = [found[R]['ops'][0] for R in (R0 + 1, R0 + 2, R0 + 3) if R in found]
    gA, gB = g_of(setA), g_of(setB)
    return dict(R0=R0, found={R: (v['Q0'], v['nul']) for R, v in found.items()}, A=gA, B=gB)


def analyse(c_x, D1, D2, p):
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
                               lam=cl['lam'], logs=cl['logs'], adjoint_has_order1=cl['adjoint_has_order1'])
    return out


def component_verdict(rA, rB):
    if rA is None:
        return 'NOT-FOUND'
    if rA['verdict'] == 'TRIANGULAR' and (rB is None or rB['verdict'] == 'TRIANGULAR'):
        return 'TRIANGULAR'
    if rA['verdict'] == 'STUCK' and rB is not None and rB['verdict'] == 'STUCK' and rA['residual']['clean'] \
            and rB['residual']['clean'] and rA['residual']['order'] == rB['residual']['order']:
        return 'STUCK'
    if rA['verdict'] == 'TRIANGULAR' or (rB is not None and rB['verdict'] == 'TRIANGULAR'):
        return 'TRIANGULAR'     # a left multiplier can only ADD factors: one triangular annihilator suffices
    return 'STUCK-UNCONFIRMED'


def run(case, p, L, Rmax, na=False):
    import ao_exact as X
    import na_exact as NX
    T0 = time.time()
    D1s, D2s, zX, zH = (NACASES if na else CASES)[case]
    D1, D2 = Fraction(D1s), Fraction(D2s)
    tag = f"intclass_{'NA_' if na else ''}{case}_{p}"
    log = open(DATA / f'{tag}.log', 'w')

    def say(*a):
        s = ' '.join(str(x) for x in a)
        print(s, flush=True)
        log.write(s + '\n')
        log.flush()
    if na:
        Qs = NX.coefficients(D1s, D2s, zX, zH, L, NX.Field(p))
    else:
        Qs = X.coefficients(D1s, D2s, zX, zH, L, X.Field(p))
    keys = sorted(set(k for q in Qs for k in q), key=repr)
    say(case, p, 'coefficients to L =', L, len(keys), 'constants', round(time.time() - T0), 's')
    res = dict(case=case, prime=p, L=L, Rmax=Rmax, nconst=len(keys), components={})
    for k in keys:
        seq = [int(q.get(k, flint.nmod(0, p))) for q in Qs]
        t1 = time.time()
        ops = operator_sets(seq, p, Rmax)
        if ops is None:
            rec = dict(verdict='NOT-FOUND')
        else:
            rA = analyse(ops['A'], D1, D2, p)
            rB = analyse(ops['B'], D1, D2, p) if ops['B'] is not None else None
            rec = dict(R0=ops['R0'], found={str(a): b for a, b in ops['found'].items()}, A=rA, B=rB,
                       verdict=component_verdict(rA, rB))
        rec['time'] = round(time.time() - t1, 1)
        res['components'][repr(k)] = rec
        say(case, repr(k), rec['verdict'], {kk: rec.get(kk) for kk in ('R0', 'found')},
            (rec.get('A') or {}).get('residual'), rec['time'], 's')
        (DATA / f'{tag}.json').write_text(json.dumps(res, indent=1, default=str))
    res['time'] = round(time.time() - T0, 1)
    (DATA / f'{tag}.json').write_text(json.dumps(res, indent=1, default=str))
    say('done', res['time'], 's')


def case_verdict(case, na=False):
    per = []
    for p in PRIMES:
        f = DATA / f"intclass_{'NA_' if na else ''}{case}_{p}.json"
        if not f.exists():
            return 'MISSING', None
        per.append(json.loads(f.read_text())['components'])
    keys = set(per[0]) | set(per[1])
    v = {k: tuple(d.get(k, {}).get('verdict', 'NOT-FOUND') for d in per) for k in keys}
    if any(a == b == 'STUCK' for a, b in v.values()):
        return 'NOT-TRIANGULAR', v
    if all(a == b == 'TRIANGULAR' for a, b in v.values()):
        return 'TRIANGULAR', v
    return 'INCOMPLETE', v


def score():
    out = {}
    # IC-0: controls on both scoring primes (python debug/intclass_controls.py <p> must have been run first)
    ok = True
    for p in CONTROL_PRIMES:
        f = DATA / f'intclass_controls_{p}.json'
        ok = ok and f.exists() and json.loads(f.read_text())['all_pass']
    out['IC-0'] = 'HELD' if ok else 'MISSED'
    vW1, dW1 = case_verdict('W1')
    vW2, dW2 = case_verdict('W2')
    out['W1'], out['W2'] = vW1, vW2
    out['IC-1'] = 'HELD' if vW1 == 'NOT-TRIANGULAR' else ('MISSED' if vW1 == 'TRIANGULAR' else 'UNDECIDED')
    out['IC-2'] = ('HELD' if vW1 == vW2 else 'MISSED') if 'MISSING' not in (vW1, vW2) else 'UNDECIDED'
    # IC-3: elliptic reading of a stuck residual (conditional on IC-1)
    if vW1 == 'NOT-TRIANGULAR':
        types = set()
        for p in PRIMES:
            comps = json.loads((DATA / f'intclass_W1_{p}.json').read_text())['components']
            for k, r in comps.items():
                if r['verdict'] == 'STUCK':
                    types.add(r['A']['residual']['type'])
        out['W1_stuck_types'] = sorted(types)
        out['IC-3'] = 'HELD' if types & {'fuchsian-elliptic-type', 'confluent-genus1-type'} else 'MISSED'
    else:
        out['IC-3'] = 'UNSCORED'
    vN1, _ = case_verdict('N1', na=True)
    vN2, _ = case_verdict('N2', na=True)
    out['NA'] = (vN1, vN2)
    out['IC-4'] = 'HELD' if vN1 == vN2 == 'TRIANGULAR' else ('MISSED' if 'NOT-TRIANGULAR' in (vN1, vN2) else
                                                             'UNDECIDED')
    # IC-5: prime agreement of every decided component verdict, both cases
    agree = True
    for d in (dW1, dW2):
        if d:
            for a, b in d.values():
                if {a, b} <= {'TRIANGULAR', 'STUCK'} and a != b:
                    agree = False
    out['IC-5'] = 'HELD' if agree and dW1 and dW2 else ('MISSED' if not agree else 'UNDECIDED')
    # IC-6: completeness in W1
    out['IC-6'] = ('HELD' if dW1 and all('NOT-FOUND' not in ab for ab in dW1.values()) else 'MISSED') if dW1 \
        else 'UNDECIDED'
    held = lambda k: out.get(k) == 'HELD'
    if not held('IC-0'):
        gate = 'VOID'
    elif held('IC-1') and held('IC-2'):
        gate = 'GO'
    elif out['IC-1'] == 'MISSED':
        gate = 'STOP'
    else:
        gate = 'BORDERLINE'
    out['gate'] = gate
    (DATA / 'intclass_score.json').write_text(json.dumps(out, indent=1, default=str))
    for k, v in out.items():
        print(k, v)


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'run':
        run(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]) if len(sys.argv) > 4 else 2000,
            int(sys.argv[5]) if len(sys.argv) > 5 else 40)
    elif mode == 'na':
        run(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]) if len(sys.argv) > 4 else 1200, 24, na=True)
    elif mode == 'score':
        score()
