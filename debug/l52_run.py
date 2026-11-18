"""Record L52 Phase-2 driver (docs/predictions/2026-11-17_aa_bc_contact_amplitudes.md). NOT to be run before the record is committed.

    python debug/l52_run.py run <case> <dps>   : exact split (aabc_exact engine) of the Coulomb (AA|BC) and, where the case
                                                 has a kappa, of the Yukawa-line (AA|BC)_kappa; F_L at dps digits,
                                                 L <= 240 -> data/l52_run_<case>_<dps>.json
    python debug/l52_run.py perL <case>        : P0(iii) independent per-L quadrature of Y at mass kappa, L = 0, 1, 4, 12
    python debug/l52_run.py fk <case>          : P0(iv) assembled f_kappa per L vs quadrature of closed-form V_kappa
    python debug/l52_run.py ctrl_perL | ctrl_fk <Ls> <dps> : registration controls K7, K9 (A1 only)
    python debug/l52_run.py score              : reads the run/perL files, applies the gate -> data/l52_score.json
    python debug/l52_run.py fire               : gate fire test on synthetic inputs (no engine)

Yukawa line: for rho_A = (zA^3/pi) e^{-mu s}, mu = 2 zA, the potential through e^{-kappa r12}/r12 is
    V_kappa = mu^4 a (e^{-kappa s} - e^{-mu s})/s - mu^3/(2 (mu^2 - kappa^2)) e^{-mu s},  a = 1/(mu^2 - kappa^2)^2,
so (AA|BC)_kappa = mu^4 a (Y_kappa - Y_mu) - mu^3/(2(mu^2-kappa^2)) W_mu,  W = -dY/dmu  (engine pieces Y, dY).
"""
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import aabc_exact as X
import na_exact as E
import l52_local as LL
import l52_predict as PR

DATA = HERE / 'data'
LMAX = 240
LS = list(range(120, LMAX + 1, 20))
REF_VALUE = {'A1': ('1/2', '0.21286449876124374093484119906', 1e-25),
             'A2': ('-1/3', '0.223136849487580864384315474643', 1e-20)}   # AB P1, stored 2000-digit sums
TOL = dict(P0_dps=1e-12, P0_perL=1e-25, P0_fk=1e-20, P0_beta=0.1, P1=1e-3, P2a=3e-3, P2b=0.1)
PERL_LS = (0, 1, 4, 12)      # Y_kappa per-L referee: L measured at registration (control K7)
FK_LS = (0, 1, 4, 12)        # f_kappa assembly referee: L measured at registration (control K9, 30 digits)


def fr(s):
    return Fraction(s)


def mpq(v):
    return mp.mpf(int(v.p)) / int(v.q)


def evalF(Q, consts):
    return [sum(mpq(v) * consts[k] for k, v in q.items()) for q in Q]


def run(case, dps):
    g = PR.CASES[case]
    D1, D2, zB, zC, zA = (fr(g[k]) for k in ('D1', 'D2', 'zB', 'zC', 'zA'))
    t0 = time.time()
    K = E.Field(None)
    pc = X.pieces(case, LMAX, K, params_override=(D1, D2, zB, zC, zA))
    f0 = X.combine(pc, zA, K)
    pk = None
    if g['kappa']:
        kap = fr(g['kappa'])
        pk = X.pieces(case, LMAX, K, which=('Y',), params_override=(D1, D2, zB, zC, kap / 2))['Y']
    print(case, 'exact split', round(time.time() - t0), 's', flush=True)
    mp.mp.dps = dps
    keys = set(k for P in [f0] + ([pk, pc['Y'], pc['dY']] if pk else []) for q in P for k in q)
    consts = {k: E.value(k, mp) for k in keys}
    tw = 2 * mp.sqrt(mp.mpf(zB.numerator) ** 3 / zB.denominator ** 3 * mp.mpf(zC.numerator) ** 3 / zC.denominator ** 3)
    F0 = [tw * v for v in evalF(f0, consts)]
    out = dict(case=case, dps=dps, F0={str(L): mp.nstr(F0[L], 40) for L in range(LMAX + 1)})
    if pk:
        mu = 2 * mp.mpf(zA.numerator) / zA.denominator
        k = mp.mpf(kap.numerator) / kap.denominator
        a = 1 / (mu ** 2 - k ** 2) ** 2
        Yk, Ym, dYm = (evalF(P, consts) for P in (pk, pc['Y'], pc['dY']))
        Fk = [tw * (mu ** 4 * a * (Yk[L] - Ym[L]) + mu ** 3 / (2 * (mu ** 2 - k ** 2)) * dYm[L]) for L in range(LMAX + 1)]
        out['Fk'] = {str(L): mp.nstr(Fk[L], 40) for L in range(LMAX + 1)}
        out['D'] = {str(L): mp.nstr(Fk[L] - F0[L], 40) for L in range(LMAX + 1)}
    if case in REF_VALUE:
        xr = mp.mpf(fr(REF_VALUE[case][0]).numerator) / fr(REF_VALUE[case][0]).denominator
        Pl = [mp.mpf(1), xr]
        for L in range(1, LMAX):
            Pl.append(((2 * L + 1) * xr * Pl[L] - L * Pl[L - 1]) / (L + 1))
        out['value'] = mp.nstr(sum(F0[L] * Pl[L] for L in range(LMAX + 1)), 35)
    out['secs'] = round(time.time() - t0)
    (DATA / f'l52_run_{case}_{dps}.json').write_text(json.dumps(out, indent=1))
    print(case, dps, 'done', out['secs'], 's', flush=True)


def perL(case, Ls=PERL_LS, out_name=None):
    """Y at mass kappa per L: engine (exact split, 1200 digits) vs angular Gauss-Legendre + mpmath Bessel quadrature."""
    g = PR.CASES[case]
    D1, D2, zB, zC = (fr(g[k]) for k in ('D1', 'D2', 'zB', 'zC'))
    kap = fr(g['kappa'])
    pk = X.pieces(case, max(Ls), E.Field(None), which=('Y',), params_override=(D1, D2, zB, zC, kap / 2))['Y']
    with mp.workdps(1200):
        consts = {k: E.value(k, mp) for q in pk for k in q}
        ex = {L: 2 * mp.sqrt((mp.mpf(zB.numerator) / zB.denominator) ** 3 * (mp.mpf(zC.numerator) / zC.denominator) ** 3) * sum(mpq(v) * consts[k] for k, v in pk[L].items()) for L in Ls}
    mp.mp.dps = 50
    d1, d2, zb, zc, k = (mp.mpf(v.numerator) / v.denominator for v in (D1, D2, zB, zC, kap))
    N = mp.sqrt(zb ** 3 * zc ** 3) / mp.pi

    def rhoL(r, L):
        f = lambda s: mp.exp(-zc * s) * mp.legendre(L, (r * r + d1 * d1 - s * s) / (2 * r * d1)) * s
        return (2 * L + 1) / mp.mpf(2) * N * mp.exp(-zb * r) * mp.quad(f, [abs(r - d1), r + d1]) / (r * d1)

    def iL(L, z):
        return mp.sqrt(mp.pi / (2 * z)) * mp.besseli(L + mp.mpf(1) / 2, z)

    def kL(L, z):
        return mp.sqrt(2 / (mp.pi * z)) * mp.besselk(L + mp.mpf(1) / 2, z)
    res = {}
    for L in Ls:
        def ker(r):
            a, b = (r, d2) if r < d2 else (d2, r)
            return k * (2 * L + 1) * iL(L, k * a) * kL(L, k * b)
        q = mp.quad(lambda r: 4 * mp.pi / (2 * L + 1) * r * r * rhoL(r, L) * ker(r), [0, d1, d2, mp.inf])
        res[str(L)] = float(abs(q / ex[L] - 1))
        print(case, L, res[str(L)], flush=True)
    (DATA / (out_name or f'l52_perL_{case}.json')).write_text(json.dumps(res, indent=1))


def fkcheck(case, Ls=FK_LS, dps=30, out_name=None):
    """P0(iv): the assembled f_kappa per L (engine pieces Y_kappa, Y_mu, dY_mu combined as in run()) against direct
    quadrature of the closed-form V_kappa: F_L = 4 pi/(2L+1) int r^2 rho_L(r) V_L(r) dr, with
    V_L(r) = (2L+1)/2 int_{-1}^{1} V_kappa(|r - A|) P_L(t) dt done in s = |r - A| (t = (r^2 + D2^2 - s^2)/(2 r D2)).
    Shares nothing with the engine but rho_L's definition (also re-done here by quadrature)."""
    g = PR.CASES[case]
    D1, D2, zB, zC, zA = (fr(g[k]) for k in ('D1', 'D2', 'zB', 'zC', 'zA'))
    kap = fr(g['kappa'])
    t0 = time.time()
    K = E.Field(None)
    pc = X.pieces(case, max(Ls), K, which=('Y', 'dY'), params_override=(D1, D2, zB, zC, zA))
    pk = X.pieces(case, max(Ls), K, which=('Y',), params_override=(D1, D2, zB, zC, kap / 2))['Y']
    with mp.workdps(1200):
        consts = {k: E.value(k, mp) for P in (pk, pc['Y'], pc['dY']) for q in P for k in q}
        mu = 2 * mp.mpf(zA.numerator) / zA.denominator
        k = mp.mpf(kap.numerator) / kap.denominator
        a = 1 / (mu ** 2 - k ** 2) ** 2
        tw = 2 * mp.sqrt((mp.mpf(zB.numerator) / zB.denominator) ** 3 * (mp.mpf(zC.numerator) / zC.denominator) ** 3)
        Yk, Ym, dYm = (evalF(P, consts) for P in (pk, pc['Y'], pc['dY']))
        eng = {L: tw * (mu ** 4 * a * (Yk[L] - Ym[L]) + mu ** 3 / (2 * (mu ** 2 - k ** 2)) * dYm[L]) for L in Ls}
    mp.mp.dps = dps
    d1, d2, zb, zc, za, kk = (mp.mpf(v.numerator) / v.denominator for v in (D1, D2, zB, zC, zA, kap))
    m = 2 * za
    al = 1 / (m * m - kk * kk) ** 2
    N = mp.sqrt(zb ** 3 * zc ** 3) / mp.pi

    def V(s):
        return m ** 4 * al * (mp.exp(-kk * s) - mp.exp(-m * s)) / s - m ** 3 / (2 * (m * m - kk * kk)) * mp.exp(-m * s)

    def rhoL(r, L):
        f = lambda s: mp.exp(-zc * s) * mp.legendre(L, (r * r + d1 * d1 - s * s) / (2 * r * d1)) * s
        return (2 * L + 1) / mp.mpf(2) * N * mp.exp(-zb * r) * mp.quad(f, [abs(r - d1), r + d1]) / (r * d1)

    def VL(r, L):
        f = lambda s: V(s) * mp.legendre(L, (r * r + d2 * d2 - s * s) / (2 * r * d2)) * s
        return (2 * L + 1) / mp.mpf(2) * mp.quad(f, [abs(r - d2), r + d2]) / (r * d2)
    res = {}
    for L in Ls:
        q = mp.quad(lambda r: 4 * mp.pi / (2 * L + 1) * r * r * rhoL(r, L) * VL(r, L), [0, d1, d2, 2 * d2, mp.inf])
        res[str(L)] = dict(engine=mp.nstr(eng[L], 25), quad=mp.nstr(q, 25), rel=float(abs(q / eng[L] - 1)))
        print(case, L, res[str(L)], round(time.time() - t0), 's', flush=True)
    res['dps'] = dps
    res['secs'] = round(time.time() - t0)
    (DATA / (out_name or f'l52_fk_{case}.json')).write_text(json.dumps(res, indent=1))


def analyse(case, F, amp, gam, xs):
    R = [mp.mpf(F[str(L)]) / (amp * LL.T(gam, xs, L)) for L in LS]
    r0, c1 = LL.ratio_extrap(LS, R)
    rho = 1 / (xs + mp.sqrt(xs * xs - 1))
    G = lambda l: mp.mpf(F[str(l)]) / rho ** l
    beta = mp.log(abs(G(221) / G(219))) / mp.log(mp.mpf(221) / 219)
    return dict(R240=float(R[-1]), R_extrap=float(r0), c1=float(c1), beta220=float(beta))


def gate(s):
    """s: {'P0': {'value': [(err, tol)], 'dps': [rel], 'perL': [rel], 'fk': [rel], 'beta': {case: (beta_f, toy72)}},
           'P1': {case: R_extrap}, 'P2a': {case: R_extrap}, 'P2b': {case: (beta_D, toy92)}}  -> (verdict, clauses).
    P0 (instrument) failing is VOID; P2b is reported and never moves the verdict."""
    p0 = s['P0']
    cl = {}
    cl['P0'] = (len(p0['value']) == 2 and all(e <= t for e, t in p0['value'])
                and len(p0['dps']) >= 5 and all(v <= TOL['P0_dps'] for v in p0['dps'])
                and len(p0['perL']) == 2 * len(PERL_LS) and all(v <= TOL['P0_perL'] for v in p0['perL'])
                and len(p0['fk']) == 2 * len(FK_LS) and all(v <= TOL['P0_fk'] for v in p0['fk'])
                and 'A3' in p0['beta'] and all(abs(b - t) <= TOL['P0_beta'] for b, t in p0['beta'].values()))
    cl['P1'] = len(s['P1']) == 3 and all(abs(r - 1) <= TOL['P1'] for r in s['P1'].values())
    cl['P2a'] = len(s['P2a']) == 2 and all(abs(r - 1) <= TOL['P2a'] for r in s['P2a'].values())
    cl['P2b_reported'] = len(s['P2b']) == 2 and all(abs(b - t) <= TOL['P2b'] for b, t in s['P2b'].values())
    near = (len(s['P1']) == 3 and len(s['P2a']) == 2 and all(abs(r - 1) <= 5 * TOL['P1'] for r in s['P1'].values())
            and all(abs(r - 1) <= 5 * TOL['P2a'] for r in s['P2a'].values()))
    if not cl['P0']:
        v = 'VOID'
    elif cl['P1'] and cl['P2a']:
        v = 'GO'
    elif near:
        v = 'BORDERLINE'
    else:
        v = 'STOP'
    return v, cl


def score():
    pred = json.loads((DATA / 'l52_predictions.json').read_text())
    mp.mp.dps = 60
    s = {'P0': {'value': [], 'dps': [], 'perL': [], 'fk': [], 'beta': {}}, 'P1': {}, 'P2a': {}, 'P2b': {}}
    detail = {}
    for case, g in PR.CASES.items():
        lo = json.loads((DATA / f'l52_run_{case}_2000.json').read_text())
        hi = json.loads((DATA / f'l52_run_{case}_2600.json').read_text())
        D1, D2 = PR.mpf(g['D1']), PR.mpf(g['D2'])
        xs = LL.xstar(D1, D2)
        for key in [k for k in ('F0', 'D') if k in lo]:
            s['P0']['dps'].append(max(float(abs(mp.mpf(lo[key][str(L)]) / mp.mpf(hi[key][str(L)]) - 1)) for L in LS))
        if case in REF_VALUE:
            s['P0']['value'].append((float(abs(mp.mpf(lo['value']) - mp.mpf(REF_VALUE[case][1]))), REF_VALUE[case][2]))
        a1 = analyse(case, lo['F0'], mp.mpf(pred[case]['Phi7']), 3.5, xs)       # R240, c1 reported for f
        s['P1'][case] = a1['R_extrap']
        s['P0']['beta'][case] = (a1['beta220'], float(pred[case]['beta_toy_72_220']))
        detail[case] = dict(P1=a1)
        if 'D' in lo:
            a2 = analyse(case, lo['D'], mp.mpf(pred[case]['DeltaPhi9']), 4.5, xs)
            s['P2a'][case] = a2['R_extrap']
            s['P2b'][case] = (a2['beta220'], float(pred[case]['beta_toy_92_220']))
            detail[case]['P2'] = a2
            s['P0']['perL'] += list(json.loads((DATA / f'l52_perL_{case}.json').read_text()).values())
            fk = json.loads((DATA / f'l52_fk_{case}.json').read_text())
            s['P0']['fk'] += [fk[str(L)]['rel'] for L in FK_LS]
    v, cl = gate(s)
    out = dict(verdict=v, clauses=cl, inputs=s, detail=detail)
    (DATA / 'l52_score.json').write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps(out, indent=1, default=str))


def fire():
    """Each plant must move the verdict to the named one."""
    good = {'P0': {'value': [(1e-30, 1e-25), (1e-24, 1e-20)], 'dps': [1e-20] * 5, 'perL': [1e-40] * 8,
                   'fk': [1e-30] * 8, 'beta': {'A1': (-4.106, -4.119), 'A2': (-4.083, -4.094), 'A3': (-4.10, -4.108)}},
            'P1': {'A1': 1.0002, 'A2': 0.9998, 'A3': 1.0001}, 'P2a': {'A1': 1.001, 'A3': 0.999},
            'P2b': {'A1': (-5.17, -5.18), 'A3': (-5.15, -5.16)}}
    import copy
    plants = []

    def plant(name, fn, want):
        x = copy.deepcopy(good)
        fn(x)
        plants.append((name, gate(x)[0], want))
    plants.append(('baseline', gate(good)[0], 'GO'))
    plant('P0 A1 value off', lambda x: x['P0']['value'].__setitem__(0, (1e-10, 1e-25)), 'VOID')
    plant('P0 A2 value 1e-18 (tol 1e-20)', lambda x: x['P0']['value'].__setitem__(1, (1e-18, 1e-20)), 'VOID')
    plant('P0 A2 value missing', lambda x: x['P0']['value'].pop(1), 'VOID')
    plant('P0 dps off', lambda x: x['P0']['dps'].__setitem__(2, 1e-6), 'VOID')
    plant('P0 perL missing', lambda x: x['P0'].__setitem__('perL', []), 'VOID')
    plant('P0 perL one L off', lambda x: x['P0']['perL'].__setitem__(3, 1e-12), 'VOID')
    plant('P0 fk assembly sign error (rel 2)', lambda x: x['P0']['fk'].__setitem__(0, 2.0), 'VOID')
    plant('P0 fk W-term factor 2 (rel 0.3)', lambda x: x['P0']['fk'].__setitem__(5, 0.3), 'VOID')
    plant('P0 fk missing', lambda x: x['P0'].__setitem__('fk', []), 'VOID')
    plant('P0 A3 nearer singularity (beta -2.0)', lambda x: x['P0']['beta'].__setitem__('A3', (-2.0, -4.108)), 'VOID')
    plant('P0 A3 beta missing', lambda x: x['P0']['beta'].pop('A3'), 'VOID')
    plant('P0 A3 fk off with P2a STOP', lambda x: (x['P0']['fk'].__setitem__(7, 1e-3), x['P2a'].__setitem__('A3', 3.0)), 'VOID')
    plant('P1 3e-3', lambda x: x['P1'].__setitem__('A2', 1.003), 'BORDERLINE')
    plant('P1 6e-3 (past 5x)', lambda x: x['P1'].__setitem__('A2', 1.006), 'STOP')
    plant('P1 wrong amplitude (x2)', lambda x: x['P1'].__setitem__('A3', 2.0), 'STOP')
    plant('P1 7/2 coefficient zero (ratio 0)', lambda x: x['P1'].__setitem__('A1', 0.0), 'STOP')
    plant('P1 case missing', lambda x: x['P1'].pop('A3'), 'STOP')
    plant('P2a kappa-dependent 7/2 (ratio ~ L, huge)', lambda x: x['P2a'].__setitem__('A1', 37.0), 'STOP')
    plant('P2a case missing', lambda x: x['P2a'].pop('A1'), 'STOP')
    plant('P2a sign flip', lambda x: x['P2a'].__setitem__('A3', -1.0), 'STOP')
    plant('P2a 1e-2', lambda x: x['P2a'].__setitem__('A3', 1.01), 'BORDERLINE')
    plant('P2b only (reported)', lambda x: x['P2b'].__setitem__('A1', (-4.1, -5.18)), 'GO')
    ok = all(got == want for _, got, want in plants)
    for n, got, want in plants:
        print(f'{n:45s} got {got:10s} want {want:10s} {"ok" if got == want else "FAIL"}')
    print('FIRE TEST', 'PASS' if ok else 'FAIL', f'{sum(g == w for _, g, w in plants)}/{len(plants)}')
    return ok


if __name__ == '__main__':
    m = sys.argv[1]
    if m == 'run':
        run(sys.argv[2], int(sys.argv[3]))
    elif m == 'perL':
        perL(sys.argv[2])
    elif m == 'fk':
        fkcheck(sys.argv[2])
    elif m == 'ctrl_perL':          # registration control K7 (Y_kappa alone), A1 only
        perL('A1', out_name='l52_ctrl_perL_A1.json')
    elif m == 'ctrl_fk':            # registration control K9 (f_kappa at low L), A1 only
        fkcheck('A1', Ls=tuple(int(v) for v in sys.argv[2].split(',')), dps=int(sys.argv[3]),
                out_name=f'l52_ctrl_fk_A1_{sys.argv[2].replace(",", "_")}_{sys.argv[3]}.json')
    elif m == 'score':
        score()
    elif m == 'fire':
        sys.exit(0 if fire() else 1)
