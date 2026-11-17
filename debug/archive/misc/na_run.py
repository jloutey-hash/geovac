"""Record NA scoring driver (docs/predictions/2026-10-28_nuclear_attraction_angle_operator.md).
    python debug/na_run.py controls      : NA0 (b), (c) on both primes
    python debug/na_run.py check <case>  : NA0 (a), exact vs data/na_<case>_num.json
    python debug/na_run.py scan <case>   : the search + analysis on both primes -> data/na_scan_<case>.json
"""
import json
import random
import sys
import time

import flint
import mpmath as mp

import na_ctrl_exact as C
import na_exact as E
import na_modanalyze as A
import na_modguess as G

PRIMES = [(1 << 61) - 1, 4611686018427387847]
CASES = {'N1': ('1.8', '2.6', '1.3', '1.0'), 'N2': ('1.5', '2.4', '1.0', '1.2')}
XSTAR = {'N1': (125, 117), 'N2': (89, 80)}
mode = sys.argv[1]

if mode == 'controls':
    for p in PRIMES:
        random.seed(7)
        gen, _ = C.control(4, 300, p, [random.randrange(1, p) for _ in range(8)])
        B = G.Basis(gen, p)
        res = {R: G.min_degree(B, R, 60, 30) for R in range(1, 7)}
        bad = list(gen)
        bad[150] = bad[150] + flint.nmod(1, p)
        Bb = G.Basis(bad, p)
        nb, _, _ = G.nullity(Bb, 4, 15, 30)
        print('p', p, {R: v for R, v in res.items()}, '| decoy nullity (4,15)', nb, 'decoy order-4',
              G.min_degree(Bb, 4, 60, 30), flush=True)

elif mode == 'check':
    case = sys.argv[2]
    D1, D2, zX, zH = CASES[case]
    mp.mp.dps = 200
    Q = E.coefficients(D1, D2, zX, zH, 40, E.Field(None))
    num = [mp.mpf(s.strip('[').split(' +/-')[0]) for s in json.load(open(f'data/na_{case}_num.json'))['F']]
    tw = 2 * mp.sqrt(mp.mpf(zX) ** 3 * mp.mpf(zH) ** 3)
    e = [abs(sum(mp.mpf(int(v.p)) / int(v.q) * E.value(k, mp) for k, v in Q[L].items()) * tw / num[L] - 1)
         for L in range(41)]
    print(case, 'NA0(a) max', mp.nstr(max(e), 3), 'at', e.index(max(e)), 'constants',
          sorted(set(k for q in Q for k in q), key=str), flush=True)

elif mode == 'scan':
    case = sys.argv[2]
    D1, D2, zX, zH = CASES[case]
    out = {}
    for p in PRIMES:
        t0 = time.time()
        Qs = E.coefficients(D1, D2, zX, zH, 1000, E.Field(p))
        keys = sorted(set(k for q in Qs for k in q), key=str)
        rnd = random.Random(2026)
        w = {k: flint.nmod(rnd.randrange(1, p), p) for k in keys}
        gen = [sum((w[k] * q.get(k, flint.nmod(0, p)) for k in keys), flint.nmod(0, p)) for q in Qs]
        comps = {str(k): [q.get(k, flint.nmod(0, p)) for q in Qs] for k in keys}
        print(case, p, 'coefficients', round(time.time() - t0), 's; constants', len(keys), flush=True)
        B = G.Basis(gen, p)
        trade = {}
        for R in range(1, 13):
            Q0, info = G.min_degree_bisect(B, R, 30, len(gen))
            trade[R] = (Q0, info)
            print(case, p, 'R', R, 'minQ', Q0, info, round(time.time() - t0), flush=True)
        found = [(R, v[0]) for R, v in trade.items() if v[0] is not None]
        rec = dict(trade={R: list(v) for R, v in trade.items()})
        if found:
            Rm, Q0 = min(found)
            nul, ks, (N, nrow) = G.nullity(B, Rm, Q0, 30)
            c = A.standard_form(ks, [N[i, 0] for i in range(N.nrows())], p)
            x0 = flint.nmod(XSTAR[case][0], p) / flint.nmod(XSTAR[case][1], p)
            inv2 = flint.nmod(1, p) / 2
            cand = {str(e): flint.nmod(e, p) for e in range(-3, 12)}
            cand.update({f'{2 * h + 1}/2': flint.nmod(2 * h + 1, p) * inv2 for h in range(-3, 10)})
            loc = A.local(c, x0, p, cand)
            pts = A.all_points(c, p)
            rec.update(R_min=Rm, Q0=Q0, nullity=nul, degrees=[ck.degree() for ck in c],
                       xstar=dict(orders=loc['orders'], regular=loc['regular'],
                                  roots=[k for k, v in loc.get('candidate_roots', {}).items() if v]),
                       infinity_regular=A.fuchs_at_infinity(c)[0],
                       points=[dict(factor=str(r['factor']), degree=r['degree'], mult=r['multiplicity'],
                                    regular=r['regular'], ints=r.get('integer_exponents'),
                                    apparent=r.get('apparent_type')) for r in pts])
            print(case, p, 'FOUND', Rm, Q0, 'nullity', nul, 'x*', rec['xstar'], 'inf regular', rec['infinity_regular'],
                  flush=True)
            for r in rec['points']:
                print('   ', r['degree'], r['mult'], r['regular'], r['ints'], r['apparent'],
                      r['factor'][:60], flush=True)
        # components, reported only
        comp = {}
        for name, seq in comps.items():
            Bc = G.Basis(seq, p)
            comp[name] = {R: G.min_degree_bisect(Bc, R, 30, len(seq))[0] for R in range(1, 7)}
        rec['components'] = comp
        print(case, p, 'components', comp, flush=True)
        out[str(p)] = rec
    json.dump(out, open(f'data/na_scan_{case}.json', 'w'), indent=0, default=str)
