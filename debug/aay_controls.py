"""Controls for a pre-registration record AY: (AA|BC) with m != 0 3d pairs on A (DNX, DNY, DXY, DXX; engine
debug/aay_exact.py).  Nothing here computes a scored quantity: no summed exact series, no decay of the exact
coefficients, no constant set beyond L = 8.  Geometry and frame as aay_exact.

    python debug/aay_controls.py referee   : R2 (AB's even-tempered Gaussian fits n = 30, 34; generic quartic moment
                                             identity p(w) G = [p(d)/16b^4 + (Lap p)(d)/16b^3 + Lap^2 p/32b^2] G on
                                             A-derivatives of erf(wR)/R), R1 (prolate quadrature against the Gauss-law
                                             potentials, incomplete-gamma form), fit floors against geovac's closed
                                             forms (3d1 3d1|1s 1s), (3d2 3d2|1s 1s) (and 3d0) at R = 2
    python debug/aay_controls.py toy       : beta_eff(L) at L = 40, 220: P_L coefficients of (x* - x)^alpha, P_L^2 of
                                             (1 - x^2)(x* - x)^alpha, P_L^4 of (1 - x^2)^2 (x* - x)^alpha (extra),
                                             alpha = 5/2 .. 11/2
    python debug/aay_controls.py precision : exact coefficients over Q to L = 400 (both cases, all integrals), values
                                             at L = 210..221 and 395..400 at a ladder of working precisions; ONLY the
                                             relative self-differences between precisions are written (f0/f2/f4 at
                                             219..221 and 398..400; the total c at 219..221); engine timings
    python debug/aay_controls.py merge     : fold the part files into debug/data/aay_controls.json
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from math import comb, factorial
from pathlib import Path

import numpy as np
from scipy.special import erf as _erf, eval_hermite

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import aabc_controls as AB  # noqa: E402
import aaa_controls as AA  # noqa: E402

DATA = HERE / 'data'
OUT = DATA / 'aay_controls.json'
POINTS = {'A1': (Fraction(1, 2), Fraction(-1, 3)), 'A2': (Fraction(0), Fraction(4, 5))}
S3 = np.sqrt(3.0)


def save_part(key: str, val) -> None:
    (DATA / f'aay_part_{key}.json').write_text(json.dumps(val, indent=1))


def frame(A: np.ndarray):
    n = A / np.linalg.norm(A)
    x, sg = n[2], np.hypot(n[0], n[1])
    ep = np.array([x, 0.0, -sg]) if sg > 0 else np.array([1.0, 0.0, 0.0])
    ey = np.array([0.0, 1.0, 0.0])
    return n, ep, ey


def mats(A: np.ndarray) -> dict:
    n, ep, ey = frame(A)
    o = np.outer
    return {'DNX': S3 / 2 * (o(n, ep) + o(ep, n)), 'DNY': S3 / 2 * (o(n, ey) + o(ey, n)),
            'DXY': S3 / 2 * (o(ep, ey) + o(ey, ep)), 'DXX': S3 / 2 * (o(ep, ep) - o(ey, ey)),
            'DA': (3 * o(n, n) - np.eye(3)) / 2}


# ---------------- R2: generic quartic (w M w)^2 on A ----------------
def _Dk(R, w):
    """D^k f, k = 0..4, f = erf(w R)/R, D = R^-1 d/dR."""
    x = w * R
    g = np.exp(-x * x)
    Ej = [_erf(x)] + [w ** j * 2 / np.sqrt(np.pi) * (-1) ** (j - 1) * eval_hermite(j - 1, x) * g for j in range(1, 5)]
    f = [sum(comb(k, j) * Ej[j] * (-1) ** (k - j) * factorial(k - j) / R ** (k - j + 1) for j in range(k + 1))
         for k in range(5)]
    return [f[0], f[1] / R, f[2] / R ** 2 - f[1] / R ** 3, f[3] / R ** 3 - 3 * f[2] / R ** 4 + 3 * f[1] / R ** 5,
            f[4] / R ** 4 - 6 * f[3] / R ** 5 + 15 * f[2] / R ** 6 - 15 * f[1] / R ** 7]


def gauss_d(al, c, zA, zB, zC, A, B, C, M) -> float:
    """(d_M d_M | 1s_B 1s_C), d_M = N_d (w M w) e^{-zA s}, tr M = 0, tr M^2 = 3/2."""
    NB, NC = np.sqrt(zB ** 3 / np.pi), np.sqrt(zC ** 3 / np.pi)
    bA = (zA ** 2 * (al[:, None] + al[None, :])).ravel()
    cA = (c[:, None] * c[None, :]).ravel() * AA.Nd2(zA)
    p, cBC, P = AA.AP._bc_pairs(al, c, zB, zC, B, C)
    cBC = cBC * NB * NC
    d = A - P
    R = np.linalg.norm(d, axis=-1)
    M2 = M @ M
    t2 = np.trace(M2)
    dMd = np.einsum('ki,ij,kj->k', d, M, d)
    dM2d = np.einsum('ki,ij,kj->k', d, M2, d)
    tot = 0.0
    for i in range(len(bA)):
        be = bA[i]
        pref = (np.pi / p) ** 1.5 * (np.pi / be) ** 1.5
        D = _Dk(R, np.sqrt(p * be / (p + be)))
        q4 = 2 * t2 * D[2] + 4 * dM2d * D[3] + dMd ** 2 * D[4]          # (d M d)^2 f
        q2 = t2 * D[1] + dM2d * D[2]                                    # (d M^2 d) f
        val = q4 / (16 * be ** 4) + 8 * q2 / (16 * be ** 3) + 16 * t2 * D[0] / (32 * be ** 2)
        tot += cA[i] * np.sum(cBC * pref * val)
    return tot


# ---------------- R1: prolate quadrature, Gauss-law potential ----------------
def V_d(w: np.ndarray, M: np.ndarray, z: float) -> np.ndarray:
    s = np.linalg.norm(w, axis=-1)
    u = w / s[..., None]
    Vl = AA.V_multipoles('d', s, z)                       # = N_d^2 c_l g_l
    uMu = np.einsum('...i,ij,...j->...', u, M, u)
    uM2u = np.einsum('...i,ij,...j->...', u, M @ M, u)
    H2 = (8 * uM2u - 4) / 14
    H4 = uMu ** 2 - H2 - 0.2
    return Vl[0] + Vl[2] / (2 / 7) * H2 + Vl[4] / (18 / 35) * H4


def prolate(M, zA, zB, zC, D1, A, nxi, neta, nphi) -> float:
    NB, NC = np.sqrt(zB ** 3 / np.pi), np.sqrt(zC ** 3 / np.pi)
    h = D1 / 2
    lam = h * (zB + zC)
    u, wu = np.polynomial.laguerre.laggauss(nxi)
    xi = 1 + u / lam
    eta, we = np.polynomial.legendre.leggauss(neta)
    phi = 2 * np.pi * np.arange(nphi) / nphi
    XI, ET, PH = np.meshgrid(xi, eta, phi, indexing='ij')
    W = (wu / lam)[:, None, None] * we[None, :, None] * (2 * np.pi / nphi)
    rcyl = h * np.sqrt(np.clip((XI ** 2 - 1) * (1 - ET ** 2), 0, None))
    zz = h * (1 + XI * ET)
    pts = np.stack([rcyl * np.cos(PH), rcyl * np.sin(PH), zz], axis=-1)
    f = np.exp(-lam - h * (zB - zC) * ET) * V_d(pts - A, M, zA) * h ** 3 * (XI ** 2 - ET ** 2)
    return NB * NC * np.sum(W * f)


def floors_exact():
    import sympy as sp
    from geovac.two_center_eri import aabb_closed_form, R_s
    out = {}
    for name, nlm in (('d1', (3, 2, 1)), ('d2', (3, 2, 2)), ('d0', (3, 2, 0))):
        e = aabb_closed_form(Fraction(3), nlm, nlm, Fraction(1), (1, 0, 0), (1, 0, 0))
        out[name] = complex(sp.N(e.subs(R_s, sp.Integer(2)), 30))
    return out


def referee():
    out = {}
    t0 = time.time()
    fits = {nn: AB.fit(nn) for nn in (30, 34)}
    print('fits done', round(time.time() - t0), 's', flush=True)
    # self-consistency: generic gauss_d / V_d on DA == record AS's gauss_DA / V_multipoles at one point
    g = {k: float(Fraction(v)) for k, v in AB.CASES['A1'].items()}
    A = g['D2'] * np.array([np.sqrt(0.75), 0.0, 0.5])
    B, C = np.zeros(3), np.array([0.0, 0.0, g['D1']])
    al, c, _ = fits[34]
    mine = gauss_d(al, c, g['zA'], g['zB'], g['zC'], A, B, C, mats(A)['DA'])
    ref = AA.gauss_DA(al, c, g['zA'], g['zB'], g['zC'], A, B, C, n=A / np.linalg.norm(A))
    p1 = prolate(mats(A)['DA'], g['zA'], g['zB'], g['zC'], g['D1'], A, 60, 60, 72)
    p2 = AA.prolate_A('d', A / np.linalg.norm(A), g['zA'], g['zB'], g['zC'], g['D1'], A, 60, 60, 72)
    out['generic_vs_AS_DA'] = dict(R2_rel=abs(mine / ref - 1), R1_rel=abs(p1 / p2 - 1))
    print('generic vs AS code on DA:', out['generic_vs_AS_DA'], flush=True)
    # floors: B = C = origin, A = (0, 0, 2), zA = zB = zC = 1
    ex = floors_exact()
    o = np.zeros(3)
    Aax = np.array([0.0, 0.0, 2.0])
    Mx = mats(Aax)
    fl = {}
    for nn, (al, c, _) in fits.items():
        fl[str(nn)] = dict(DNX_minus_3d1=gauss_d(al, c, 1.0, 1.0, 1.0, Aax, o, o, Mx['DNX']) - ex['d1'].real,
                           DNY_minus_3d1=gauss_d(al, c, 1.0, 1.0, 1.0, Aax, o, o, Mx['DNY']) - ex['d1'].real,
                           DXX_minus_3d2=gauss_d(al, c, 1.0, 1.0, 1.0, Aax, o, o, Mx['DXX']) - ex['d2'].real,
                           DXY_minus_3d2=gauss_d(al, c, 1.0, 1.0, 1.0, Aax, o, o, Mx['DXY']) - ex['d2'].real,
                           DA_minus_3d0=gauss_d(al, c, 1.0, 1.0, 1.0, Aax, o, o, Mx['DA']) - ex['d0'].real)
        print('floor', nn, fl[str(nn)], flush=True)
    out['floors'] = dict(exact={k: [v.real, v.imag] for k, v in ex.items()}, by_n=fl,
                         note='B = C = origin, A = (0,0,2): DNX = DNY = (3d+1 3d+1|1s1s), DXX = DXY = (3d+2 3d+2|1s1s) '
                              '(complex |d_m|^2 = mean of the two real squares, equal here by symmetry); hydrogenic '
                              'Z = 3 3d = Slater 3d at 1')
    for case in AB.CASES:
        g = {k: float(Fraction(v)) for k, v in AB.CASES[case].items()}
        B, C = np.zeros(3), np.array([0.0, 0.0, g['D1']])
        for xq in POINTS[case]:
            x = float(xq)
            A = g['D2'] * np.array([np.sqrt(1 - x * x), 0.0, x])
            Ms = mats(A)
            for integral in ('DNX', 'DNY', 'DXY', 'DXX'):
                M = Ms[integral]
                r2 = {str(nn): gauss_d(al, c, g['zA'], g['zB'], g['zC'], A, B, C, M) for nn, (al, c, _) in fits.items()}
                r1 = {','.join(map(str, gr)): prolate(M, g['zA'], g['zB'], g['zC'], g['D1'], A, *gr)
                      for gr in ((60, 60, 72), (90, 90, 108), (120, 120, 144))}
                rec = dict(x=str(xq), R2=r2, R2_spread_30_34=r2['34'] - r2['30'], R1=r1,
                           R1_minus_R2_34=r1['120,120,144'] - r2['34'])
                out[f'{case}_{integral}_x={xq}'] = rec
                print(case, integral, str(xq), json.dumps(rec), round(time.time() - t0), 's', flush=True)
    save_part('referee', out)


# ---------------- toys ----------------
def toy():
    import mpmath as mp
    mp.mp.dps = 120
    nodes = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(9, mp.mp.prec)
    Lmax = 226
    t0 = time.time()

    def anti(cv):
        gv = [mp.mpf(0)] * len(cv)
        for L in range(1, len(cv) - 1):
            gv[L] = cv[L - 1] / (2 * L - 1) - cv[L + 1] / (2 * L + 3)
        return gv
    out = {}
    for case in AB.CASES:
        g = {k: mp.mpf(Fraction(v).numerator) / Fraction(v).denominator for k, v in AB.CASES[case].items()}
        xs = (g['D1'] ** 2 + g['D2'] ** 2) / (2 * g['D1'] * g['D2'])
        rho = 1 / (xs + mp.sqrt(xs * xs - 1))
        rec = dict(x_star=mp.nstr(xs, 20), rho=mp.nstr(rho, 20))
        for alpha in ('5/2', '7/2', '9/2', '11/2'):
            al = mp.mpf(Fraction(alpha).numerator) / Fraction(alpha).denominator
            cv = [mp.mpf(0)] * (Lmax + 2)
            for x, w in nodes:
                hv = w * (xs - x) ** al
                p0, p1 = mp.mpf(1), x
                cv[0] += hv
                cv[1] += hv * x
                for n in range(1, Lmax + 1):
                    p0, p1 = p1, ((2 * n + 1) * x * p1 - n * p0) / (n + 1)
                    cv[n + 1] += hv * p1
            cv = [(2 * n + 1) / mp.mpf(2) * v for n, v in enumerate(cv)]
            b2 = anti(anti(cv))                  # P_L^2 coefficients of (1 - x^2) h
            b4 = anti(anti(b2))                  # P_L^4 coefficients of (1 - x^2)^2 h  (extra)

            def beta(F, L):
                a_, b_ = F[L + 1] / rho ** (L + 1), F[L - 1] / rho ** (L - 1)
                return float(mp.log(abs(a_ / b_)) / mp.log(mp.mpf(L + 1) / (L - 1)))
            rec[alpha] = dict(P_L={L: beta(cv, L) for L in (40, 220)}, P_L2={L: beta(b2, L) for L in (40, 220)},
                              P_L4_extra={L: beta(b4, L) for L in (40, 220)})
            print(case, alpha, rec[alpha], round(time.time() - t0), 's', flush=True)
        # check of the P_L^4 map by direct projection (alpha = 7/2, L = 6, 9)
        al = mp.mpf(7) / 2
        cv = [mp.mpf(0)] * 16
        for x, w in nodes:
            for n in range(16):
                cv[n] += w * (xs - x) ** al * mp.legendre(n, x)
        cv = [(2 * n + 1) / mp.mpf(2) * v for n, v in enumerate(cv)]
        b4 = anti(anti(anti(anti(cv))))
        chk = {}
        for L in (6, 9):
            tot = mp.mpf(0)
            for x, w in nodes:
                P4 = (1 - x * x) ** 2 * mp.diff(lambda t: mp.legendre(L, t), x, 4)
                tot += w * (1 - x * x) ** 2 * (xs - x) ** al * P4
            tot *= (2 * L + 1) / mp.mpf(2) * mp.factorial(L - 4) / mp.factorial(L + 4)
            chk[L] = mp.nstr(abs(tot / b4[L] - 1), 3)
        rec['P_L4_map_vs_projection_rel'] = chk
        out[case] = rec
    out['definition'] = ('beta_eff(L) = log|F_{L+1}/F_{L-1}| / log((L+1)/(L-1)), F_L = coefficient / rho^L, rho = '
                         '1/(x* + sqrt(x*^2 - 1)); P_L^m = (1 - x^2)^{m/2} d^m P_L/dx^m (no CS phase); P_L^2 (P_L^4) '
                         'coefficients of (1 - x^2)^{m/2} h = Legendre coefficients of h integrated m times; 768-node '
                         'Gauss-Legendre, 120 digits (as record AS)')
    save_part('toy', out)


# ---------------- precision ladder + timings (values never written or printed) ----------------
def precision(Lmax: int = 400, cases=None, dpss=(1200, 1600, 2000, 2400, 3000, 3600)):
    import mpmath as mp
    import aay_exact as Y
    import na_exact as E
    W1 = (219, 221) if Lmax >= 230 else (Lmax - 10, Lmax - 8)
    W2 = (Lmax - 2, Lmax)
    keepL = set(range(W1[0] - 9, W1[1] + 1)) | set(range(W2[0] - 5, Lmax + 1))
    names = ('DNX', 'DNY', 'DXY', 'DXX', 'DA')
    res = {}
    for case in (cases or Y.CASES):
        K = E.Field(None)
        eng = Y.Engine(case, K, Lmax)
        D2 = eng.D2
        tw = {}
        cval = {}
        vals = {d: {nm: {} for nm in names} for d in dpss}      # vals[dps][name][(m, L)] = mpf
        tG = {'U4': 0.0, 'U2': 0.0, 'U0': 0.0}
        tasm = {nm: 0.0 for nm in names}
        tev = 0.0
        t_start = time.time()
        for L in range(Lmax + 1):
            Gs = {}
            for kern in ('U4', 'U2', 'U0'):
                t0 = time.time()
                Gs[kern] = eng.G_kernel(kern, L, Y.NMAX[kern])
                tG[kern] += time.time() - t0
            for nm in names:
                t0 = time.time()
                for m, rec in Y.channel_jets(nm, L).items():
                    acc = {}
                    for kern, jet in rec.items():
                        for (j, p), cc in jet.items():
                            a = K(cc) * K.pw(D2, -p)
                            for k, v in Gs[kern][j].items():
                                E.add(acc, k, a * v)
                    q = {k: v for k, v in acc.items() if v != 0}
                    tasm[nm] += time.time() - t0
                    if L in keepL:
                        t1 = time.time()
                        for d in dpss:
                            with mp.workdps(d):
                                if d not in tw:
                                    tw[d] = A_twopiN(case)
                                cs = cval.setdefault(d, {})
                                tot = mp.mpf(0)
                                for k, v in q.items():
                                    if k not in cs:
                                        cs[k] = E.value(k, mp)
                                    tot += mp.mpf(int(v.p)) / int(v.q) * cs[k]
                                vals[d][nm][(m, L)] = tw[d] * tot
                        tev += time.time() - t1
                    t0 = time.time()
            # free the exact cache for this L
            for k in [k for k in eng.cache if k[2] == L]:
                del eng.cache[k]
            eng.dens_cache.pop(L, None)
            if L % 50 == 0:
                print(case, 'L', L, 'elapsed', round(time.time() - t_start), 's', flush=True)
        # self-differences
        rec = dict(timing_s=dict(G_build_shared={k: round(v, 1) for k, v in tG.items()},
                                 assembly_per_integral={k: round(v, 1) for k, v in tasm.items()},
                                 evaluation_all_dps=round(tev, 1), total=round(time.time() - t_start, 1)))
        diffs = {}
        for i in range(len(dpss) - 1):
            d0, d1 = dpss[i], dpss[i + 1]
            pair = {}
            with mp.workdps(max(dpss)):
                for nm in names:
                    ms = sorted(set(m for (m, _) in vals[d1][nm]))
                    for m in ms:
                        for win in (W1, W2):
                            r = max(abs(vals[d1][nm][(m, L)] / vals[d0][nm][(m, L)] - 1)
                                    if vals[d0][nm][(m, L)] != 0 else mp.inf
                                    for L in range(win[0], win[1] + 1) if (m, L) in vals[d1][nm])
                            pair[f'{nm}_f{m}_L{win[0]}-{win[1]}'] = mp.nstr(r, 3)
                    if nm != 'DA':     # total c at 219..221 (needs L >= j only)
                        cs = {}
                        for d in (d0, d1):
                            f = {m: [mp.mpf(0)] * (Lmax + 1) for m in ms}
                            for (m, L), v in vals[d][nm].items():
                                f[m][L] = v
                            cs[d] = Y.total_coeffs(f)
                        r = max(abs(cs[d1][j] / cs[d0][j] - 1) for j in range(W1[0], W1[1] + 1))
                        pair[f'{nm}_c_L{W1[0]}-{W1[1]}'] = mp.nstr(r, 3)
            diffs[f'{d0}_vs_{d1}'] = pair
            worst = max(mp.mpf(v) for v in pair.values())
            print(case, f'{d0} vs {d1}: worst relative self-difference {mp.nstr(worst, 3)}', flush=True)
        rec['self_differences'] = diffs
        res[case] = rec
        print(case, 'timing', rec['timing_s'], flush=True)
        save_part(f'precision_{case}' + ('' if Lmax == 400 else f'_L{Lmax}'), rec)
    return res


def A_twopiN(case):
    import aay_exact as Y
    return Y.twopiN(case)


def merge():
    d = json.loads(OUT.read_text()) if OUT.exists() else {}
    for p in sorted(DATA.glob('aay_part_*.json')):
        key = p.stem[len('aay_part_'):]
        val = json.loads(p.read_text())
        if key.startswith('AY0_'):
            d.setdefault('AY0', {})[key[4:]] = val
        elif key.startswith('precision_'):
            d.setdefault('precision', {})[key[len('precision_'):]] = val
        else:
            d[key] = val
    OUT.write_text(json.dumps(d, indent=1))
    print('merged', sorted(d))


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'precision':   # precision <Lmax> <case> [dps,dps,...]
        dp = tuple(int(v) for v in sys.argv[4].split(',')) if len(sys.argv) > 4 else (1200, 1600, 2000, 2400, 3000, 3600)
        precision(int(sys.argv[2]), [sys.argv[3]], dp)
    else:
        {'referee': referee, 'toy': toy, 'merge': merge}[cmd]()
