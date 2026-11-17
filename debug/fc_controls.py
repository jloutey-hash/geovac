"""Controls for the four-centre (XY|ZW) Legendre-coefficient engine (debug/fc_engine.py).  No scored quantity is
computed here: coefficient values at L >= 40 are formed only inside precision self-differences and Legendre sums,
never printed or stored; no beta_eff of the engine's sequence; no amplitude at a branch point.

    python debug/fc_controls.py k1      multipoles rho2_L, rho1_L vs direct mpmath quadrature (60 digits)
    python debug/fc_controls.py k2      R2 referee: even-tempered n-Gaussian fits (n = 22, 26), closed-form s ERIs
    python debug/fc_controls.py main    K3 (sum to L = 260 vs R2) + K5 (level A vs B self-differences) + timings
    python debug/fc_controls.py k4      reduction b = a, Z = W: engine vs record AB's aabc_exact (case A1), L <= 30
    python debug/fc_controls.py k6      toy decay instrument on (x* - x)^nu, nu = 5/2, 7/2, 9/2
Output: debug/data/fc_controls.json (merged per command).
"""
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import mpmath as mp
from flint import arb, ctx

sys.path.insert(0, str(Path(__file__).parent))
import fc_engine as FE

OUT = Path(__file__).parent / "data" / "fc_controls.json"
XS = ('1/2', '-1/3')


def load():
    return json.loads(OUT.read_text()) if OUT.exists() else {}


def save_key(k, v):
    d = load()
    d[k] = v
    OUT.write_text(json.dumps(d, indent=1))


def fl(x):
    return float(x.mid()) if isinstance(x, arb) else float(x)


# ------------------------------------------------------------------ K1
def k1():
    mp.mp.dps = 60
    out = {}
    worst = {'rho2': 0.0, 'rho1': 0.0}
    Ls = (0, 1, 4, 12, 30)
    for case in ('G1', 'G2'):
        p = FE.params(case)
        a, b, D1 = p['a'], p['b'], p['D1']
        rs2 = [Fraction(1, 2), a - Fraction(1, 1000), a + Fraction(1, 1000), D1, b - Fraction(1, 1000),
               b + Fraction(1, 1000), Fraction(3) * b / 2]
        rs1 = [Fraction(1, 3), D1 - Fraction(1, 1000), D1 + Fraction(1, 1000), Fraction(5, 2)]
        zZ, zW, zX, zY = (mp.mpf(p[k].numerator) / p[k].denominator for k in ('zZ', 'zW', 'zX', 'zY'))
        am, bm, Dm = (mp.mpf(v.numerator) / v.denominator for v in (a, b, D1))
        rec = []
        ctx.prec = FE.LEVELS['A']['prec']
        for r in rs2:
            eng = FE.rho2_point(case, r, 30)
            rm = mp.mpf(r.numerator) / r.denominator
            f = lambda mu: mp.exp(-zZ * mp.sqrt(rm * rm + am * am - 2 * rm * am * mu)
                                  - zW * mp.sqrt(rm * rm + bm * bm - 2 * rm * bm * mu))
            pts = [-1, 0, mp.mpf('0.9'), mp.mpf('0.99'), mp.mpf('0.999'), mp.mpf('0.9999'), mp.mpf('0.99999'),
                   mp.mpf('0.999999'), mp.mpf('0.9999999'), 1]
            for L in Ls:
                ref = (2 * L + 1) * mp.quad(lambda mu: f(mu) * mp.legendre(L, mu), pts) / 2
                e = abs(mp.mpf(eng[L].mid().str(70, radius=False)) - ref) / abs(ref)
                worst['rho2'] = max(worst['rho2'], float(e))
                rec.append(dict(kind='rho2', r=str(r), L=L, rel_err=float(e)))
        for r in rs1:
            eng = FE.rho1_point(case, r, 30)
            rm = mp.mpf(r.numerator) / r.denominator
            f = lambda mu: mp.exp(-zX * rm - zY * mp.sqrt(rm * rm + Dm * Dm - 2 * rm * Dm * mu))
            pts = [-1, 0, mp.mpf('0.9'), mp.mpf('0.99'), mp.mpf('0.999'), mp.mpf('0.9999'), mp.mpf('0.99999'),
                   mp.mpf('0.999999'), 1]
            for L in Ls:
                ref = (2 * L + 1) * mp.quad(lambda mu: f(mu) * mp.legendre(L, mu), pts) / 2
                e = abs(mp.mpf(eng[L].mid().str(70, radius=False)) - ref) / abs(ref)
                worst['rho1'] = max(worst['rho1'], float(e))
                rec.append(dict(kind='rho1', r=str(r), L=L, rel_err=float(e)))
        out[case] = rec
        print(case, 'worst so far', worst, flush=True)
    save_key('K1', dict(worst=worst, points=out, note='engine level-A theta grid vs mp.quad (tanh-sinh, '
                        'mu variable, split toward mu = 1) at 60 digits'))


# ------------------------------------------------------------------ K2 (R2 referee)
def gauss_pair(al, c, z1, z2, C1, C2):
    a = z1 ** 2 * al[:, None]
    b = z2 ** 2 * al[None, :]
    p = a + b
    K = np.exp(-a * b / p * np.sum((C1 - C2) ** 2))
    P = (a[..., None] * C1 + b[..., None] * C2) / p[..., None]
    return (c[:, None] * c[None, :] * K).ravel(), p.ravel(), P.reshape(-1, 3)


def gauss_four(al, c, p, x):
    """(XY|ZW) unnormalised, every e^{-z r} replaced by the fit sum c_i e^{-al_i z^2 r^2}."""
    from aabc_controls import eri_charge
    g = {k: float(v) for k, v in p.items()}
    X = np.zeros(3)
    Y = np.array([0, 0, g['D1']])
    nv = np.array([np.sqrt(1 - x * x), 0, x])
    Z, W = g['a'] * nv, g['b'] * nv
    c1, p1, P1 = gauss_pair(al, c, g['zX'], g['zY'], X, Y)
    c2, p2, P2 = gauss_pair(al, c, g['zZ'], g['zW'], Z, W)
    tot = 0.0
    for i in range(len(c1)):
        tot += c1[i] * np.sum(c2 * eri_charge(p1[i], P1[i], p2, P2))
    return tot


def k2():
    from aabc_controls import fit
    t0 = time.time()
    fits = {n: fit(n) for n in (22, 26)}
    print('fits', {n: f[2] for n, f in fits.items()}, round(time.time() - t0), 's', flush=True)
    out = {'fit_res': {str(n): f[2] for n, f in fits.items()}}
    for case in ('G1', 'G2', 'K4A1'):
        p = FE.params(case)
        out[case] = {}
        for xs in XS:
            x = float(Fraction(xs))
            v = {str(n): gauss_four(f[0], f[1], p, x) for n, f in fits.items()}
            v['spread'] = abs(v['22'] - v['26'])
            out[case][xs] = v
            print(case, xs, v, flush=True)
    out['t'] = time.time() - t0
    save_key('K2', out)


# ------------------------------------------------------------------ K3 + K5 + timings
def main():
    d = load()
    R2 = d.get('K2')
    res = {}
    for case in ('G1', 'G2'):
        FA, ia = FE.coefficients(case, 260, 'A', verbose=True)
        FB, ib = FE.coefficients(case, 260, 'B', verbose=True)
        ctx.prec = FE.LEVELS['B']['prec']
        k5 = {}
        for L in list(range(219, 222)) + list(range(258, 261)):
            k5[str(L)] = float(abs((FA[L] - FB[L]) / FB[L]).mid())
        low = max(float(abs((FA[L] - FB[L]) / FB[L]).mid()) for L in range(40))
        k3 = {}
        for xs in XS:
            sA = FE.legendre_sum(FA, xs)
            sB = FE.legendre_sum(FB, xs)
            ent = dict(engine_A=fl(sA), engine_A_minus_B=fl(sA - sB))
            if R2 and case in R2:
                ent['R2_26'] = R2[case][xs]['26']
                ent['diff_vs_R2_26'] = fl(sA) - R2[case][xs]['26']
                ent['diff_vs_R2_22'] = fl(sA) - R2[case][xs]['22']
            k3[xs] = ent
        lan = FE.landau(case)
        res[case] = dict(K5_rel_selfdiff=k5, K5_max_rel_selfdiff_L_lt_40=low, K3=k3,
                         landau={k: (str(v) if isinstance(v, Fraction) else v) for k, v in lan.items()},
                         timing_A=ia, timing_B=ib)
        print(case, json.dumps(res[case], indent=1), flush=True)
    save_key('main', res)


# ------------------------------------------------------------------ K4
def k4():
    import na_exact as E
    import aabc_exact as AB
    t0 = time.time()
    K = E.Field(None)
    pc = AB.pieces('A1', 30, K)
    f = AB.combine(pc, Fraction(9, 10), K)
    t_ab = time.time() - t0
    mp.mp.dps = 500                                    # the exact split cancels heavily (AB K6')
    ab = [sum(mp.mpf(int(v.p)) / int(v.q) * E.value(k, mp) for k, v in q.items()) for q in f]
    F, info = FE.coefficients('K4A1', 30, 'A', verbose=True)
    ctx.prec = 340
    zA = mp.mpf(9) / 10
    pred = 2 * mp.pi ** 2 / zA ** 3                    # ours = (2 pi^2 / zA^3) * sum q c  if N = N_B N_C
    rows = []
    worst = 0.0
    for L in range(31):
        ours = mp.mpf(F[L].mid().str(90, radius=False))
        ratio = ours / ab[L]
        e = abs(ratio / pred - 1)
        worst = max(worst, float(e))
        rows.append(dict(L=L, rel_dev=float(e)))
    # value level vs record AB's R2 referee, normalised: N_B N_C zA^3/pi
    FF, _ = FE.coefficients('K4A1', 260, 'A')
    ctx.prec = 340
    s = FE.legendre_sum(FF, '1/2')
    zB, zC = mp.mpf(13) / 10, mp.mpf(1)
    norm = mp.sqrt(zB ** 3 / mp.pi) * mp.sqrt(zC ** 3 / mp.pi) * zA ** 3 / mp.pi
    val = mp.mpf(s.mid().str(60, radius=False)) * norm
    out = dict(worst_rel_dev_L_le_30=worst, per_L=rows, ratio_convention='ours / (sum q c) = 2 pi^2 / zA^3',
               value_x_half_normalised=float(val), AB_R2_referee='0.2128644988 (+-3e-11)',
               diff_vs_AB_R2=float(val - mp.mpf('0.2128644988')), t_aabc_exact=t_ab, timing=info)
    print(json.dumps({k: v for k, v in out.items() if k != 'per_L'}, indent=1), flush=True)
    save_key('K4', out)


# ------------------------------------------------------------------ K6 toys
def k6():
    ctx.prec = 640
    n = 700
    nodes = [arb.legendre_p_root(n, k, weight=True) for k in range(n)]
    nodes = [(arb(x.mid()), arb(w.mid())) for x, w in nodes]
    Lmax = 222
    out = {}
    for case in ('G1', 'G2'):
        lan = FE.landau(case)
        xs = lan['xZ'] if lan['dominant'] == 'Z' else lan['xW']
        rho = lan['rhoZ'] if lan['dominant'] == 'Z' else lan['rhoW']
        xa, ra = FE.A(xs), FE.A(rho)
        out[case] = dict(x_star=str(xs), rho=str(rho))
        for nu in ('5/2', '7/2', '9/2'):
            nua = FE.A(nu)
            c = [arb(0)] * (Lmax + 1)
            for x, w in nodes:
                fv = w * ((xa - x).log() * nua).exp()
                p0, p1 = arb(1), x
                c[0] += fv
                c[1] += fv * x
                for L in range(1, Lmax):
                    p0, p1 = p1, ((2 * L + 1) * x * p1 - L * p0) / (L + 1)
                    c[L + 1] += fv * p1
            G = [c[L] * (2 * L + 1) / 2 / ra ** L for L in range(Lmax + 1)]
            beta = {}
            for L in (60, 120, 220):
                beta[str(L)] = float((abs(G[L + 1] / G[L - 1]).log() / (arb(L + 1) / (L - 1)).log()).mid())
            out[case]['nu=' + nu] = beta
            print(case, nu, beta, flush=True)
    save_key('K6', out)


if __name__ == '__main__':
    {'k1': k1, 'k2': k2, 'main': main, 'k4': k4, 'k6': k6}[sys.argv[1]]()
