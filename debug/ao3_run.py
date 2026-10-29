"""Record AO3 scoring driver (docs/predictions/2026-10-28_two_electron_angle_operator_exact.md).
    python debug/ao3_run.py check <case>          : AO3-0(a), exact vs data/ao3_<case>_num.json (ao_radial2, L <= 30)
    python debug/ao3_run.py scan <case> <prime>   : coefficients to L = 2000 mod p, search, analysis
                                                    -> data/ao3_scan_<case>_<p>.json
"""
import json
import random
import sys
import time

import flint
import mpmath as mp

import ao_exact as X
import na_modanalyze as A
import na_modguess as G

CASES = {'S1': ('1.8', '2.6', '1.3', '1.0'), 'W1': ('1.7', '2.5', '1.4', '1.1'), 'W2': ('1.6', '2.2', '1.0', '1.2')}
XSTAR = {'S1': (125, 117), 'W1': (457, 425), 'W2': (185, 176)}
LMAX, RMAX, MARGIN = 2000, 40, 30
mode, case = sys.argv[1], sys.argv[2]
D1, D2, zX, zH = CASES[case]

if mode == 'check':
    mp.mp.dps = 250
    Q = X.coefficients(D1, D2, zX, zH, 30, X.Field(None))
    num = [mp.mpf(s.strip('[').split(' +/-')[0]) for s in json.load(open(f'data/ao3_{case}_num.json'))['F']]
    pre = 4 * mp.mpf(zX) ** 3 * mp.mpf(zH) ** 3
    e = [abs(sum(mp.mpf(int(v.p)) / int(v.q) * X.value(k, mp) for k, v in Q[L].items()) * pre / num[L] - 1)
         for L in range(31)]
    print(case, 'AO3-0(a) max', mp.nstr(max(e), 3), 'at', e.index(max(e)), 'constants',
          len(set(k for q in Q for k in q)), flush=True)

elif mode == 'scan':
    p = int(sys.argv[3])
    t0 = time.time()
    Qs = X.coefficients(D1, D2, zX, zH, LMAX, X.Field(p))
    keys = sorted(set(k for q in Qs for k in q), key=repr)
    rnd = random.Random(2026)
    w = {k: flint.nmod(rnd.randrange(1, p), p) for k in keys}
    gen = [sum((w[k] * q.get(k, flint.nmod(0, p)) for k in keys), flint.nmod(0, p)) for q in Qs]
    print(case, p, 'coefficients', round(time.time() - t0), 's; constants', len(keys), flush=True)
    B = G.Basis(gen, p)
    trade, first = {}, None
    for R in range(1, RMAX + 1):
        Q0, info = G.min_degree_bisect(B, R, MARGIN, len(gen))
        trade[R] = (Q0, info)
        print(case, p, 'R', R, 'minQ', Q0, info, round(time.time() - t0), flush=True)
        if Q0 is not None and first is None:
            first = R
        if first is not None and R >= first + 4:
            break
    rec = dict(constants=[repr(k) for k in keys], trade={R: list(v) for R, v in trade.items()})
    if first is not None:
        Rm, Q0 = first, trade[first][0]
        nul, ks, (N, nrow) = G.nullity(B, Rm, Q0, MARGIN)
        c = A.standard_form(ks, [N[i, 0] for i in range(N.nrows())], p)
        x0 = flint.nmod(XSTAR[case][0], p) / flint.nmod(XSTAR[case][1], p)
        inv2 = flint.nmod(1, p) / 2
        cand = {str(e): flint.nmod(e, p) for e in range(-3, 31)}
        cand.update({f'{2 * h + 1}/2': flint.nmod(2 * h + 1, p) * inv2 for h in range(-3, 21)})
        loc = A.local(c, x0, p, cand)
        pts = []
        for r in A.all_points(c, p, Emax=80):
            f = r['factor']
            d = dict(degree=f.degree(), multiplicity=r['multiplicity'], regular=r['regular'],
                     ints=r.get('integer_exponents'))
            if f.degree() == 1 and r['regular']:
                z = -f.coeffs()[0] / f.coeffs()[1]
                lz = A.local(c, z, p, cand)
                d['x_mod_p'] = int(z)
                d['candidate_roots'] = [k for k, v in lz['candidate_roots'].items() if v]
            pts.append(d)
        rec.update(R_min=Rm, Q0=Q0, nullity=nul, degrees=[ck.degree() for ck in c],
                   xstar=dict(orders=loc['orders'], regular=loc['regular'],
                              roots=[k for k, v in loc.get('candidate_roots', {}).items() if v]),
                   infinity_regular=A.fuchs_at_infinity(c)[0], points=pts)
        print(case, p, 'FOUND', Rm, Q0, 'nullity', nul, 'x*', rec['xstar'], 'inf regular', rec['infinity_regular'],
              flush=True)
        for d in pts:
            print('   ', d, flush=True)
    rec['seconds'] = time.time() - t0
    json.dump(rec, open(f'data/ao3_scan_{case}_{p}.json', 'w'), indent=0, default=str)
