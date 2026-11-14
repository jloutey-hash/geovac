"""Record AB scoring driver (docs/predictions/2026-11-14_aa_bc_three_centre.md), Phase 2.
    python debug/aabc_run.py p1p3b [Lmax] [dps]   : P1 (value at the referee points) + P3b (decay) per case, exact
                                                    split at dps digits; -> data/aabc_p1p3b.json
    python debug/aabc_run.py perL                 : per-L check of the screened pieces vs independent quadrature
                                                    (Bessel functions + angular quadrature of the density)
    python debug/aabc_run.py p2 <prime>           : constant sets L = 0..1000 mod prime -> data/aabc_p2_<p>.json
    python debug/aabc_run.py p3a <case>           : NA's operator search on the generic sequence, both NA primes
    python debug/aabc_run.py p4 <prime>           : IC pipeline on A1's mu-carrying components
"""
import json
import random
import sys
import time
from fractions import Fraction
from pathlib import Path

import flint
import mpmath as mp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import aabc_exact as X
import na_exact as E

DATA = HERE / 'data'
REF = {'A1': dict(R2={'18': '0.21286449875185', '22': '0.21286449877733', '26': '0.21286449875979'},
                  R2_stated='0.2128644988', R1_finest='0.2128644972'),
       'A2': dict(R2={'18': '0.22313684936064', '22': '0.22313684950325', '26': '0.22313684949202'},
                  R2_stated='0.2231368495', R1_finest='0.2231368505')}
NA_PRIMES = [(1 << 61) - 1, 4611686018427387847]


def mpq(v):
    return mp.mpf(int(v.p)) / int(v.q)


def evalF(Q, consts):
    return [sum(mpq(v) * consts[k] for k, v in q.items()) for q in Q]


def p1p3b(Lmax=240, dps=2000):
    out = {}
    for case in X.CASES:
        D1, D2, zB, zC, zA = X.params(case)
        t0 = time.time()
        pc = X.pieces(case, Lmax, E.Field(None))
        f = X.combine(pc, zA, E.Field(None))
        print(case, 'exact split', round(time.time() - t0), 's', flush=True)
        mp.mp.dps = dps
        keys = set(k for P in list(pc.values()) + [f] for q in P for k in q)
        consts = {k: E.value(k, mp) for k in keys}
        tw = 2 * mp.sqrt(mp.mpf(zB.numerator) ** 3 / zB.denominator ** 3 * mp.mpf(zC.numerator) ** 3 / zC.denominator ** 3)
        F = [tw * v for v in evalF(f, consts)]
        FY = [tw * v for v in evalF(pc['Y'], consts)]
        FN = [tw * v for v in evalF(pc['NA'], consts)]
        FW = [-tw * v for v in evalF(pc['dY'], consts)]
        xr = mp.mpf(Fraction(X.REF_X[case]).numerator) / Fraction(X.REF_X[case]).denominator
        Pl = [mp.mpf(1), xr]
        for L in range(1, Lmax):
            Pl.append(((2 * L + 1) * xr * Pl[L] - L * Pl[L - 1]) / (L + 1))
        val = sum(F[L] * Pl[L] for L in range(Lmax + 1))
        parts = {n: mp.nstr(sum(G[L] * Pl[L] for L in range(Lmax + 1)), 25) for n, G in
                 (('NA', FN), ('Y', FY), ('W', FW))}
        tail = mp.nstr(abs(F[Lmax] * Pl[Lmax]), 5)
        diffs = {('R2_n' + k): float(val - mp.mpf(v)) for k, v in REF[case]['R2'].items()}
        diffs['R2_stated'] = float(val - mp.mpf(REF[case]['R2_stated']))
        diffs['R1_finest'] = float(val - mp.mpf(REF[case]['R1_finest']))
        # P3b: decay instrument
        xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
        xsm = mp.mpf(xs.numerator) / xs.denominator
        rho = 1 / (xsm + mp.sqrt(xsm * xsm - 1))

        def beta(Fs):
            G = [Fs[L] / rho ** L for L in range(Lmax + 1)]
            return {L: float(mp.log(abs(G[L + 1] / G[L - 1])) / mp.log(mp.mpf(L + 1) / (L - 1)))
                    for L in list(range(20, Lmax, 20)) + [Lmax - 20]}
        sg = sorted(set(int(mp.sign(v)) for v in F[20:]))
        sg5 = sorted(set(int(mp.sign(v)) for v in F[5:]))
        first_const_sign = None
        for L0 in range(Lmax, 0, -1):
            if len(set(int(mp.sign(v)) for v in F[L0 - 1:])) > 1:
                first_const_sign = L0
                break
        rec = dict(value=mp.nstr(val, 30), tail_term=tail, diffs=diffs, parts=parts, rho=mp.nstr(rho, 20),
                   beta_f=beta(F), beta_Y=beta(FY), beta_NA=beta(FN), beta_W=beta(FW),
                   signs_f_L_ge_20=sg, signs_f_L_ge_5=sg5, f_single_sign_from_L=first_const_sign,
                   n_constants=len(set(k for q in f for k in q)), dps=dps, Lmax=Lmax,
                   F_first=[mp.nstr(v, 15) for v in F[:6]])
        # precision floor check: the same beta at dps - 300 (K6' style)
        out[case] = rec
        print(case, json.dumps(rec, indent=1), round(time.time() - t0), 's', flush=True)
    (DATA / f'aabc_p1p3b_{dps}.json').write_text(json.dumps(out, indent=1))


def perL(Ls=(0, 1, 4, 12, 30), dps=50):
    """Independent route per L: rho_L by Gauss-Legendre angular quadrature of e^{-zB r - zC |r - D1 z|}, kernels from
    mpmath besseli/besselk, radial mp.quad split at D1, D2.  Compares 2 pi N q_L of each piece (NA, Y, W)."""
    mp.mp.dps = dps
    out = {}
    for case in X.CASES:
        D1, D2, zB, zC, zA = (mp.mpf(v.numerator) / v.denominator for v in X.params(case))
        mu = 2 * zA
        Lmax = max(Ls)
        pc = X.pieces(case, Lmax, E.Field(None))
        keys = set(k for P in pc.values() for q in P for k in q)
        with mp.workdps(1200):            # the exact split cancels digits: evaluate it far above the quadrature
            consts = {k: E.value(k, mp) for k in keys}
            zBh, zCh = (mp.mpf(v.numerator) / v.denominator for v in X.params(case)[2:4])
            twh = 2 * mp.sqrt(zBh ** 3 * zCh ** 3)
            exact = {L: dict(NA=twh * sum(mpq(v) * consts[k] for k, v in pc['NA'][L].items()),
                             Y=twh * sum(mpq(v) * consts[k] for k, v in pc['Y'][L].items()),
                             W=-twh * sum(mpq(v) * consts[k] for k, v in pc['dY'][L].items())) for L in Ls}
        N = mp.sqrt(zB ** 3 * zC ** 3) / mp.pi
        glx = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(7, mp.mp.prec)   # 192 nodes on [-1, 1]
        cache = {}

        def rhoL(r, L):
            # s = |r - D1 z|: t = (r^2 + D1^2 - s^2) / (2 r D1), dt = -s ds / (r D1); integrand analytic in s
            if (r, L) in cache:
                return cache[(r, L)]
            a, b = abs(r - D1), r + D1
            h, c = (b - a) / 2, (b + a) / 2
            tot = mp.mpf(0)
            for u, w in glx:
                s = c + h * u
                t = (r * r + D1 * D1 - s * s) / (2 * r * D1)
                p0, p1 = mp.mpf(1), t
                for l in range(1, L):
                    p0, p1 = p1, ((2 * l + 1) * t * p1 - l * p0) / (l + 1)
                PL = p0 if L == 0 else p1
                tot += w * mp.exp(-zC * s) * PL * s
            v = (2 * L + 1) / mp.mpf(2) * N * mp.exp(-zB * r) * tot * h / (r * D1)
            cache[(r, L)] = v
            return v

        def iL(L, z):
            return mp.sqrt(mp.pi / (2 * z)) * mp.besseli(L + mp.mpf(1) / 2, z)

        def kL(L, z):
            return mp.sqrt(2 / (mp.pi * z)) * mp.besselk(L + mp.mpf(1) / 2, z)

        res = {}
        for L in Ls:
            def kerY(r):
                a, b = (r, D2) if r < D2 else (D2, r)
                return mu * (2 * L + 1) * iL(L, mu * a) * kL(L, mu * b)

            def kerW(r):            # -dY/dmu by numerical differentiation (independent of the engine's identity)
                a, b = (r, D2) if r < D2 else (D2, r)
                return -mp.diff(lambda m: m * (2 * L + 1) * iL(L, m * a) * kL(L, m * b), mu)

            def kerN(r):
                a, b = (r, D2) if r < D2 else (D2, r)
                return a ** L / b ** (L + 1)
            vals = {}
            for name, ker in (('NA', kerN), ('Y', kerY), ('W', kerW)):
                q = mp.quad(lambda r: 4 * mp.pi / (2 * L + 1) * r * r * rhoL(r, L) * ker(r), [0, D1, D2, mp.inf])
                ex = exact[L][name]
                vals[name] = dict(quad=mp.nstr(q, 20), rel=mp.nstr(abs(q / ex - 1), 3))
            res[L] = vals
            print(case, L, vals, flush=True)
        out[case] = res
    (DATA / 'aabc_perL.json').write_text(json.dumps(out, indent=1))


def _gl(n):
    xs, ws = [], []
    for k, (x, w) in enumerate(_gauss_legendre(n)):
        yield x, w


def _gauss_legendre(n):
    # mpmath's Gauss-Legendre nodes via the quadrature class
    gl = mp.calculus.quadrature.GaussLegendre(mp.mp)
    deg = 1
    while 3 * 2 ** (deg - 1) < n:
        deg += 1
    nodes = gl.calc_nodes(deg, mp.mp.prec)
    return nodes


def const_sets(case, p, Lmax=1000):
    D1, D2, zB, zC, zA = X.params(case)
    K = E.Field(p)
    t0 = time.time()
    pc = X.pieces(case, Lmax, K)
    f = X.combine(pc, zA, K)
    return pc, f, round(time.time() - t0)


def candidate_set(case, na_keys):
    """C exactly as listed in the record's P2 (negative-power rule), as engine keys."""
    D1, D2, zB, zC, zA = X.params(case)
    mu = 2 * zA
    C = set(na_keys)
    P_in = -zC * D1 - mu * D2
    inner = [-zB + s * zC + t * mu for s in (1, -1) for t in (1, -1)]
    shared = [-zB - zC + t * mu for t in (1, -1)]
    for b in inner:
        C.add((P_in, 'E1', -b * D1))
        C.add((P_in, 'ln', abs(b)))
        C.add((P_in, '1', None))
        C.add((P_in + b * D1, '1', None))
    for b in shared:
        for s in (1, -1):
            C.add((s * zC * D1 - mu * D2, 'E1', -b * D1))
    for s in (1, -1):
        C.add((s * zC * D1 - mu * D2, 'E1', (zB + zC - mu) * D2))
        for t in (1, -1):
            C.add((s * zC * D1 + t * mu * D2, 'E1', (zB + zC + mu) * D2))
    # endpoint exponentials of the middle and outer components
    for s in (1, -1):
        P = s * zC * D1 - mu * D2
        for b in shared:
            C.add((P + b * D1, '1', None))
            C.add((P + b * D2, '1', None))
        for t in (1, -1):
            C.add((s * zC * D1 + t * mu * D2 + (-zB - zC - mu) * D2, '1', None))
    return C


def p2(p):
    out = {}
    for case in X.CASES:
        D1, D2, zB, zC, zA = X.params(case)
        pc, f, sec = const_sets(case, p)
        na_keys = set(k for q in pc['NA'] for k in q)
        C = candidate_set(case, na_keys)
        sets = [frozenset(q) for q in f]
        S2 = sets[2]
        changes = [L for L in range(3, len(sets)) if sets[L] != S2]
        union = set().union(*sets[2:])
        # mu-free: keys unchanged when zA moves by 1/997 (engine at L <= 40)
        pert = X.params(case)[:4] + (zA + Fraction(1, 997),)
        K = E.Field(p)
        pc2 = X.pieces(case, 40, K, params_override=pert)
        f2 = X.combine(pc2, pert[4], K)
        union2 = set(k for q in f2 for k in q)
        mufree = sorted((k for k in union if k in union2), key=repr)
        e1na = [k for k in na_keys if k[1] == 'E1']
        y_on_e1 = [L for L in range(len(f)) if any(k in pc['Y'][L] or k in pc['dY'][L] for k in e1na)]
        rec = dict(seconds=sec, n_constants_L2=len(S2), L_where_set_changes=changes[:20],
                   n_union=len(union), outside_C=sorted(map(repr, union - C)),
                   C_size=len(C), C_unused=sorted(map(repr, C - union)),
                   gamma_present=[repr(k) for k in union if k[1] == 'g'],
                   logs_present=[repr(k) for k in union if k[1] == 'ln'],
                   mu_free=[repr(k) for k in mufree], na_keys=sorted(map(repr, na_keys)),
                   mu_free_subset_NA=set(mufree) <= na_keys,
                   Y_or_W_on_NA_E1_at_L=y_on_e1[:20],
                   constants=sorted(map(repr, union)),
                   setsY=len(set(k for q in pc['Y'] for k in q)), setsdY=len(set(k for q in pc['dY'] for k in q)),
                   L01=[sorted(map(repr, sets[0])), sorted(map(repr, sets[1]))])
        rec['P2a'] = (not changes) and not (union - C)
        rec['P2b'] = not rec['gamma_present']
        rec['P2c'] = bool(rec['logs_present'])
        rec['P2d'] = rec['mu_free_subset_NA'] and not y_on_e1
        out[case] = rec
        print(case, p, json.dumps({k: v for k, v in rec.items() if k not in ('constants',)}, indent=1), flush=True)
        (DATA / f'aabc_p2_{p}.json').write_text(json.dumps(out, indent=1))
        p3a(case, p, f, sec)


def p3a(case, p, f, sec):
    import na_modanalyze as A
    import na_modguess as G
    D1, D2, zB, zC, zA = X.params(case)
    xs = (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)
    fn = DATA / f'aabc_p3a_{case}_{p}.json'
    out = json.loads(fn.read_text()) if fn.exists() else {}
    if True:
        keys = sorted(set(k for q in f for k in q), key=str)
        rnd = random.Random(2026)
        w = {k: flint.nmod(rnd.randrange(1, p), p) for k in keys}
        gen = [sum((w[k] * q.get(k, flint.nmod(0, p)) for k in keys), flint.nmod(0, p)) for q in f]
        print(case, p, 'coefficients', sec, 's; constants', len(keys), flush=True)
        B = G.Basis(gen, p)
        trade = {}
        t0 = time.time()
        for R in range(1, 25):
            Q0, info = G.min_degree_bisect(B, R, 30, len(gen))
            trade[R] = (Q0, info)
            print(case, p, 'R', R, 'minQ', Q0, info, round(time.time() - t0), flush=True)
            if Q0 is not None:
                break
        found = [(R, v[0]) for R, v in trade.items() if v[0] is not None]
        rec = dict(nconst=len(keys), trade={R: list(v) for R, v in trade.items()})
        if found:
            Rm, Q0 = min(found)
            nul, ks, (N, nrow) = G.nullity(B, Rm, Q0, 30)
            c = A.standard_form(ks, [N[i, 0] for i in range(N.nrows())], p)
            x0 = flint.nmod(xs.numerator, p) / flint.nmod(xs.denominator, p)
            inv2 = flint.nmod(1, p) / 2
            cand = {str(e): flint.nmod(e, p) for e in range(-3, 14)}
            cand.update({f'{2 * h + 1}/2': flint.nmod(2 * h + 1, p) * inv2 for h in range(-3, 12)})
            loc = A.local(c, x0, p, cand)
            pts = A.all_points(c, p)
            rec.update(R_min=Rm, Q0=Q0, nullity=nul, degrees=[ck.degree() for ck in c],
                       xstar=dict(orders=loc['orders'], regular=loc['regular'],
                                  roots=[k for k, v in loc.get('candidate_roots', {}).items() if v]),
                       infinity_regular=A.fuchs_at_infinity(c)[0],
                       points=[dict(factor=str(r['factor']), degree=r['degree'], mult=r['multiplicity'],
                                    regular=r['regular'], ints=r.get('integer_exponents'),
                                    apparent=r.get('apparent_type')) for r in pts])
            print(case, p, 'FOUND', json.dumps(rec, default=str)[:3000], flush=True)
        out[str(p)] = rec
        fn.write_text(json.dumps(out, indent=1, default=str))


def p4(p, L=1200, Rmax=24):
    import intclass_score as S
    case = 'A1'
    D1, D2, zB, zC, zA = X.params(case)
    T0 = time.time()
    K = E.Field(p)
    pc = X.pieces(case, L, K)
    f = X.combine(pc, zA, K)
    na_keys = set(k for q in pc['NA'] for k in q)
    keys = sorted(set(k for q in f for k in q) - na_keys, key=repr)
    log = open(DATA / f'aabc_p4_{p}.log', 'w')

    def say(*a):
        s = ' '.join(str(x) for x in a)
        print(s, flush=True)
        log.write(s + '\n')
        log.flush()
    say(case, p, 'coefficients to L =', L, len(keys), 'mu-carrying components', round(time.time() - T0), 's')
    res = dict(case=case, prime=p, L=L, Rmax=Rmax, components={})
    for k in keys:
        seq = [int(q.get(k, flint.nmod(0, p))) for q in f]
        t1 = time.time()
        ops = S.operator_sets(seq, p, Rmax)
        if ops is None:
            rec = dict(verdict='NOT-FOUND')
        else:
            rA = S.analyse(ops['A'], D1, D2, p)
            rB = S.analyse(ops['B'], D1, D2, p) if ops['B'] is not None else None
            rec = dict(R0=ops['R0'], found={str(a): b for a, b in ops['found'].items()}, A=rA, B=rB,
                       verdict=S.component_verdict(rA, rB))
        rec['time'] = round(time.time() - t1, 1)
        res['components'][repr(k)] = rec
        say(case, repr(k), rec['verdict'], {kk: rec.get(kk) for kk in ('R0', 'found')},
            (rec.get('A') or {}).get('residual'), rec['time'], 's')
        (DATA / f'aabc_p4_{p}.json').write_text(json.dumps(res, indent=1, default=str))
    res['time'] = round(time.time() - T0, 1)
    (DATA / f'aabc_p4_{p}.json').write_text(json.dumps(res, indent=1, default=str))
    say('done', res['time'], 's')


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'p1p3b':
        p1p3b(int(sys.argv[2]) if len(sys.argv) > 2 else 240, int(sys.argv[3]) if len(sys.argv) > 3 else 2000)
    elif mode == 'perL':
        perL()
    elif mode == 'p2':
        p2(int(sys.argv[2]))
    elif mode == 'p4':
        p4(int(sys.argv[2]))
