"""Controls for a pre-registration record on PY, PXp, DA (non-s pairs on A against 1s_B 1s_C; engine debug/aaa_exact.py).
Nothing here computes a scored quantity: no exact constant set beyond L = 8, no decay of the exact coefficients,
no exact summed value.

Geometry of AB/AP/AQ/AD: B at the origin, C at (0, 0, D1), A at D2 (sin g, 0, cos g), x = cos g.
    PY  = (2p_y,A 2p_y,A | 1s_B 1s_C)          (p out of the ABC plane)
    PXp = (2p_x',A 2p_x',A | 1s_B 1s_C),  x' = (cos g, 0, -sin g)   (p in the plane, perpendicular to BA)
    DA  = (3ds_A 3ds_A | 1s_B 1s_C),  sigma along n = A/|A|, N_d^2 h2(w)^2 e^{-2 zA s}

    python debug/aaa_controls.py referee : R2 (AB's even-tempered Gaussian fits n = 30, 34; p pair as gauss_PA with
                                           the direction e = y or x'; d pair by the fourth-order moment identity
                                           p(s) G = [p(d)/(16 b^4) + (Lap p)(d)/(16 b^3) + Lap^2 p/(32 b^2)] G on
                                           A-derivatives of the s-type ERI), R1 (prolate quadrature against the
                                           multipole potentials, incomplete-gamma form), fit floors against geovac's
                                           closed forms (2p+1 2p+1 | 1s 1s) and (3d0 3d0 | 1s 1s) at R = 2
    python debug/aaa_controls.py toy     : beta_eff(L) at L = 40, 220 for (x* - x)^alpha, alpha = 3/2 .. 9/2,
                                           plain P_L series and the P_L^2 series of (1 - x^2)(x* - x)^alpha
    python debug/aaa_controls.py timing  : wall time of the exact coefficients (PY a+b, DA) to L = 30, 45, 60 mod p
                                           and over Q (results discarded unread), extrapolations
Output debug/data/aaa_controls.json (sections merged).
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.special import gammainc, gammaincc, gamma as Gam, erf as _erf

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import aabc_controls as AB  # noqa: E402
import aapc_controls as AP  # noqa: E402

OUT = HERE / 'data' / 'aaa_controls.json'
POINTS = {'A1': (Fraction(1, 2), Fraction(-1, 3)), 'A2': (Fraction(0), Fraction(4, 5))}


def load() -> dict:
    return json.loads(OUT.read_text()) if OUT.exists() else {}


def save(d: dict) -> None:
    OUT.write_text(json.dumps(d, indent=1))


def Nd2(z: float) -> float:
    return 5 * (2 * z) ** 7 / (4 * np.pi * 720)


# ---------------- R2 for the d pair on A ----------------
def _radial_ops():
    """callables (R, u, w) -> d_n^a Lap^b f,  f = erf(w R)/R, for the (a, b) the quartic needs.  Stable form:
    d_n^2 h = D h + R^2 u^2 D^2 h,  d_n^4 h = 3 D^2 h + 6 R^2 u^2 D^3 h + R^4 u^4 D^4 h  (D = R^-1 d/dR, from
    d_z D^k h = z D^{k+1} h), with f^(k) from erf and Hermite functions times e^{-w^2 R^2} (no exp overflow)."""
    from math import comb, factorial
    from scipy.special import eval_hermite
    Rs, us = sp.symbols('R u', positive=True)
    fs = sp.Function('f')(Rs)
    fk = sp.symbols('f0:5')

    def Lap(h):
        return sp.diff(h, Rs, 2) + 2 * sp.diff(h, Rs) / Rs

    def D(h):
        return sp.diff(h, Rs) / Rs
    exprs = {}
    for (a, b) in ((0, 0), (2, 0), (0, 1), (4, 0), (2, 1), (0, 2)):
        h = fs
        for _ in range(b):
            h = Lap(h)
        if a == 2:
            h = D(h) + Rs ** 2 * us ** 2 * D(D(h))
        elif a == 4:
            h = 3 * D(D(h)) + 6 * Rs ** 2 * us ** 2 * D(D(D(h))) + Rs ** 4 * us ** 4 * D(D(D(D(h))))
        h = sp.expand(h)
        for k in range(4, 0, -1):
            h = h.subs(sp.Derivative(fs, (Rs, k)), fk[k])
        h = h.subs(fs, fk[0])
        exprs[(a, b)] = sp.lambdify((Rs, us) + fk, h, 'numpy')

    def fders(R, w):
        x = w * R
        g = np.exp(-x * x)
        E = [_erf(x)] + [w ** j * 2 / np.sqrt(np.pi) * (-1) ** (j - 1) * eval_hermite(j - 1, x) * g for j in range(1, 5)]
        return [sum(comb(k, j) * E[j] * (-1) ** (k - j) * factorial(k - j) / R ** (k - j + 1) for j in range(k + 1))
                for k in range(5)]
    ops = {}
    for key, fn in exprs.items():
        ops[key] = (lambda fn: (lambda R, u, w: fn(R, u, *fders(R, w))))(fn)
    return ops


def _quartic_parts():
    """h = h2^2 = (9 u^4 - 6 u^2 S + S^2)/4 (n = z); Lap h, Lap^2 h as {(a, b): coef} in u^a S^b."""
    x, y, z = sp.symbols('x y z', real=True)
    S_, u_ = sp.symbols('S u')
    h = sp.expand(((3 * z ** 2 - (x ** 2 + y ** 2 + z ** 2)) / 2) ** 2)

    def lap(g):
        return sp.expand(sum(sp.diff(g, v, 2) for v in (x, y, z)))

    def to_uS(g):
        g = sp.expand(g.subs(y, 0).subs(x, sp.sqrt(S_ - z ** 2)))
        g = sp.expand(g.subs(z, u_))
        P = sp.Poly(g, u_, S_)
        return {(a, b): float(c) for (a, b), c in P.terms()}
    return to_uS(h), to_uS(lap(h)), to_uS(lap(lap(h)))


_OPS = None
_QP = None


def gauss_DA(al, c, zA, zB, zC, A, B, C, n=None) -> float:
    global _OPS, _QP
    if _OPS is None:
        _OPS, _QP = _radial_ops(), _quartic_parts()
    if n is None:
        n = A / np.linalg.norm(A)
    NB, NC = np.sqrt(zB ** 3 / np.pi), np.sqrt(zC ** 3 / np.pi)
    bA = (zA ** 2 * (al[:, None] + al[None, :])).ravel()
    cA = (c[:, None] * c[None, :]).ravel() * Nd2(zA)
    p, cBC, P = AP._bc_pairs(al, c, zB, zC, B, C)
    cBC = cBC * NB * NC
    d = A - P
    R = np.linalg.norm(d, axis=-1)
    un = d @ n / R
    h, lh, llh = _QP
    tot = 0.0
    for i in range(len(bA)):
        be = bA[i]
        pref = (np.pi / p) ** 1.5 * (np.pi / be) ** 1.5
        wv = np.sqrt(p * be / (p + be))
        D = {k: op(R, un, wv) for k, op in _OPS.items()}
        val = (sum(cf * D[k] for k, cf in h.items()) / (16 * be ** 4)
               + sum(cf * D[k] for k, cf in lh.items()) / (16 * be ** 3)
               + sum(cf * D[k] for k, cf in llh.items()) / (32 * be ** 2))
        tot += cA[i] * np.sum(cBC * pref * val)
    return tot


# ---------------- R1: multipole potentials in incomplete-gamma form (stable in float64) ----------------
def V_multipoles(kind: str, s, z: float) -> dict:
    """V_l(s) of the normalized p (density zA^5/pi wz^2 e^{-2 z s}) or d density; f_l = c_l t^k e^{-mu t}."""
    mu = 2 * z
    if kind == 'p':
        k, cl = 2, {0: z ** 5 / np.pi / 3, 2: 2 * z ** 5 / np.pi / 3}
    else:
        k, cl = 4, {l: Nd2(z) * f for l, f in ((0, 1 / 5), (2, 2 / 7), (4, 18 / 35))}
    out = {}
    for l, c in cl.items():
        ai, ao = l + 2 + k + 1, 1 - l + k + 1        # int_0^s t^{ai-1} e^{-mu t}, int_s^oo t^{ao-1} e^{-mu t}
        inner = gammainc(ai, mu * s) * Gam(ai) / mu ** ai
        outer = gammaincc(ao, mu * s) * Gam(ao) / mu ** ao
        out[l] = 4 * np.pi / (2 * l + 1) * c * (s ** (-l - 1) * inner + s ** l * outer)
    return out


def prolate_A(kind: str, e_fn, zA, zB, zC, D1, A, nxi, neta, nphi) -> float:
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
    z = h * (1 + XI * ET)
    pts = np.stack([rcyl * np.cos(PH), rcyl * np.sin(PH), z], axis=-1)
    sv = pts - A
    s = np.linalg.norm(sv, axis=-1)
    cosn = (sv @ e_fn) / s
    Vl = V_multipoles(kind, s, zA)
    V = sum(Vl[l] * np.polynomial.legendre.legval(cosn, [0] * l + [1]) for l in Vl)
    f = np.exp(-lam - h * (zB - zC) * ET) * V * h ** 3 * (XI ** 2 - ET ** 2)
    return NB * NC * np.sum(W * f)


def floors_exact():
    from geovac.two_center_eri import aabb_closed_form, R_s
    out = {}
    for name, args in (('p11', (Fraction(2), (2, 1, 1), (2, 1, 1), Fraction(1), (1, 0, 0), (1, 0, 0))),
                       ('p10', (Fraction(2), (2, 1, 0), (2, 1, 0), Fraction(1), (1, 0, 0), (1, 0, 0))),
                       ('d20', (Fraction(3), (3, 2, 0), (3, 2, 0), Fraction(1), (1, 0, 0), (1, 0, 0)))):
        e = aabb_closed_form(*args)
        out[name] = complex(sp.N(e.subs(R_s, sp.Integer(2)), 30))
    return out


def referee():
    d = load()
    out = {}
    t0 = time.time()
    fits = {nn: AB.fit(nn) for nn in (30, 34)}
    print('fits done', round(time.time() - t0), 's', flush=True)
    ex = floors_exact()
    o = np.zeros(3)
    Aax = np.array([0, 0, 2.0])
    fl = {}
    for nn, (al, c, _) in fits.items():
        py = AP.gauss_PA(al, c, 1.0, 1.0, 1.0, Aax, o, o, n=np.array([0, 1.0, 0]))
        pz = AP.gauss_PA(al, c, 1.0, 1.0, 1.0, Aax, o, o, n=np.array([0, 0, 1.0]))
        dz = gauss_DA(al, c, 1.0, 1.0, 1.0, Aax, o, o, n=np.array([0, 0, 1.0]))
        fl[str(nn)] = dict(py_minus_p11=py - ex['p11'].real, pz_minus_p10=pz - ex['p10'].real,
                           dsigma_minus_d20=dz - ex['d20'].real)
        print('floor', nn, fl[str(nn)], flush=True)
    out['floors'] = dict(exact={k: [v.real, v.imag] for k, v in ex.items()}, by_n=fl,
                         note='(2p+1 2p+1|1s1s) = (2py 2py|1s1s) for B, C on the A axis (m = +-2 cross term '
                              'integrates to zero); hydrogenic Z=2 2p = Slater 2p at 1, Z=3 3d = Slater 3d at 1')
    for case in AB.CASES:
        g = {k: float(Fraction(v)) for k, v in AB.CASES[case].items()}
        B, C = np.zeros(3), np.array([0.0, 0.0, g['D1']])
        for xq in POINTS[case]:
            x = float(xq)
            sg = np.sqrt(1 - x * x)
            A = g['D2'] * np.array([sg, 0.0, x])
            dirs = {'PY': np.array([0.0, 1.0, 0.0]), 'PXp': np.array([x, 0.0, -sg]), 'DA': A / np.linalg.norm(A)}
            for integral, e in dirs.items():
                fn = gauss_DA if integral == 'DA' else AP.gauss_PA
                r2 = {str(nn): fn(al, c, g['zA'], g['zB'], g['zC'], A, B, C, n=e) for nn, (al, c, _) in fits.items()}
                kind = 'd' if integral == 'DA' else 'p'
                r1 = {','.join(map(str, gr)): prolate_A(kind, e, g['zA'], g['zB'], g['zC'], g['D1'], A, *gr)
                      for gr in ((60, 60, 72), (90, 90, 108), (120, 120, 144))}
                rec = dict(x=str(xq), R2=r2, R2_spread_30_34=r2['34'] - r2['30'], R1=r1,
                           R1_minus_R2_34=r1['120,120,144'] - r2['34'])
                out[f'{case}_{integral}_x={xq}'] = rec
                print(case, integral, str(xq), json.dumps(rec), round(time.time() - t0), 's', flush=True)
    # orientation identity on the referee itself: PY + PXp + PA = 3 S, with S the spherical average (R2, n = 34)
    d['referee'] = out
    save(d)


# ---------------- toy controls for the decay instrument ----------------
def toy():
    import mpmath as mp
    mp.mp.dps = 120
    GL = mp.calculus.quadrature.GaussLegendre(mp.mp)
    t0 = time.time()
    nodes = GL.calc_nodes(9, mp.mp.prec)          # 768 nodes
    print('nodes', len(nodes), round(time.time() - t0), 's', flush=True)
    Lmax = 226
    out = {}

    def anti(c):
        """Legendre coefficients of an antiderivative of sum c_k P_k (index list)."""
        g = [mp.mpf(0)] * len(c)
        for L in range(1, len(c) - 1):
            g[L] = c[L - 1] / (2 * L - 1) - c[L + 1] / (2 * L + 3)
        return g
    for case in AB.CASES:
        g = {k: mp.mpf(Fraction(v).numerator) / Fraction(v).denominator for k, v in AB.CASES[case].items()}
        xs = (g['D1'] ** 2 + g['D2'] ** 2) / (2 * g['D1'] * g['D2'])
        rho = 1 / (xs + mp.sqrt(xs * xs - 1))
        rec = dict(x_star=mp.nstr(xs, 20), rho=mp.nstr(rho, 20))
        for alpha in ('3/2', '5/2', '7/2', '9/2'):
            al = mp.mpf(Fraction(alpha).numerator) / Fraction(alpha).denominator
            c = [mp.mpf(0)] * (Lmax + 2)
            for x, w in nodes:
                hv = w * (xs - x) ** al
                p0, p1 = mp.mpf(1), x
                c[0] += hv
                c[1] += hv * x
                for n in range(1, Lmax + 1):
                    p0, p1 = p1, ((2 * n + 1) * x * p1 - n * p0) / (n + 1)
                    c[n + 1] += hv * p1
            c = [(2 * n + 1) / mp.mpf(2) * v for n, v in enumerate(c)]
            b2 = anti(anti(c))            # P_L^2 coefficients of (1 - x^2) h  (L >= 2)

            def beta(F, L):
                a_, b_ = F[L + 1] / rho ** (L + 1), F[L - 1] / rho ** (L - 1)
                return float(mp.log(abs(a_ / b_)) / mp.log(mp.mpf(L + 1) / (L - 1)))
            rec[alpha] = dict(P_L={L: beta(c, L) for L in (40, 220)}, P_L2={L: beta(b2, L) for L in (40, 220)},
                              signs_P_L_ge_20=sorted(set(int(mp.sign(v)) for v in c[20:221])),
                              signs_P_L2_ge_20=sorted(set(int(mp.sign(v)) for v in b2[20:221])))
            if alpha == '7/2':   # check of the P_L^2 recurrence by direct projection, L = 5, 9
                chk = {}
                for L in (5, 9):
                    tot = mp.mpf(0)
                    for x, w in nodes:
                        p = mp.legendre(L, x)
                        pd2 = mp.diff(lambda t: mp.legendre(L, t), x, 2) if False else None
                        # P_L^2 = (1 - x^2) P_L'' via the ODE: (1 - x^2) P'' = 2 x P' - L(L+1) P
                        p1 = L * (mp.legendre(L - 1, x) - x * p) / (1 - x * x)
                        P2 = 2 * x * p1 - L * (L + 1) * p
                        tot += w * (1 - x * x) * (xs - x) ** al * P2
                    tot *= (2 * L + 1) / mp.mpf(2) * mp.factorial(L - 2) / mp.factorial(L + 2)
                    chk[L] = mp.nstr(abs(tot / b2[L] - 1), 3)
                rec['P_L2_recurrence_vs_projection_rel'] = chk
            print(case, alpha, rec[alpha], round(time.time() - t0), 's', flush=True)
        out[case] = rec
    out['definition'] = ('beta_eff(L) = log|F_{L+1}/F_{L-1}| / log((L+1)/(L-1)) with F_L = coefficient / rho^L; '
                         'P_L: Legendre coefficients of (x* - x)^alpha; P_L2: coefficients b_L of '
                         '(1 - x^2)(x* - x)^alpha = sum b_L P_L^2(x), P_L^2 = (1 - x^2) P_L\'\' (no CS phase), from '
                         'b = antiderivative twice of the Legendre coefficients; 768-node Gauss-Legendre, 120 digits')
    d = load()
    d['toy'] = out
    save(d)


# ---------------- timing (results discarded unread) ----------------
def timing():
    import aaa_exact as A
    import na_exact as E
    out = {}
    case = 'A1'
    for fld in ('p', 'Q'):
        K = E.Field((1 << 61) - 1) if fld == 'p' else E.Field(None)
        rec = {}
        for L in (30, 45, 60):
            for integral in ('PY', 'DA'):
                t0 = time.time()
                res = A.coefficients(case, L, K, integral)
                del res
                rec[f'{integral}_L{L}'] = round(time.time() - t0, 2)
                print(fld, integral, L, rec[f'{integral}_L{L}'], 's', flush=True)
        for integral in ('PY', 'DA'):
            t30, t60 = rec[f'{integral}_L30'], rec[f'{integral}_L60']
            k = float(np.log(t60 / t30) / np.log(2))
            rec[f'{integral}_exponent_30_60'] = round(k, 2)
            target = 1000 if fld == 'p' else 220
            rec[f'{integral}_extrapolated_L{target}_s'] = round(t60 * (target / 60) ** k)
        out[fld] = rec
    out['note'] = ('case A1; PY = a and b series (G of U2 to second D2-derivative, l <= L+2, plus W); DA = U4 to fourth '
                   'derivative + U2 + W; exponent from L = 30 -> 60, power-law extrapolation (an over-estimate if the '
                   'per-L cost saturates, an under-estimate if big-rational growth steepens over Q)')
    d = load()
    d['timing'] = out
    save(d)



# ---------------- amplitude toys for record AS (P4), added by the PM before registration ----------------
def amp() -> None:
    """Reference coefficients at L = 218..222: ta = P_L coefficients of (x* - x)^{7/2}; tb = P_L^2 coefficients of
    (1 - x^2)(x* - x)^{7/2}; tc = P_L coefficients of (x* - x)^{9/2}.  Same quadrature as toy()."""
    import mpmath as mp
    mp.mp.dps = 120
    nodes = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(9, mp.mp.prec)
    Lmax = 226

    def anti(c):
        g = [mp.mpf(0)] * len(c)
        for L in range(1, len(c) - 1):
            g[L] = c[L - 1] / (2 * L - 1) - c[L + 1] / (2 * L + 3)
        return g
    out = {}
    for case in AB.CASES:
        g = {k: mp.mpf(Fraction(v).numerator) / Fraction(v).denominator for k, v in AB.CASES[case].items()}
        xs = (g['D1'] ** 2 + g['D2'] ** 2) / (2 * g['D1'] * g['D2'])
        rec = dict(x_star=mp.nstr(xs, 30), x_star_sq_minus_1=mp.nstr(xs * xs - 1, 30))
        for alpha in ('7/2', '9/2'):
            al = mp.mpf(Fraction(alpha).numerator) / Fraction(alpha).denominator
            c = [mp.mpf(0)] * (Lmax + 2)
            for x, w in nodes:
                hv = w * (xs - x) ** al
                p0, p1 = mp.mpf(1), x
                c[0] += hv
                c[1] += hv * x
                for n in range(1, Lmax + 1):
                    p0, p1 = p1, ((2 * n + 1) * x * p1 - n * p0) / (n + 1)
                    c[n + 1] += hv * p1
            c = [(2 * n + 1) / mp.mpf(2) * v for n, v in enumerate(c)]
            if alpha == '7/2':
                b2 = anti(anti(c))
                rec['ta'] = {L: mp.nstr(c[L], 40) for L in range(218, 223)}
                rec['tb'] = {L: mp.nstr(b2[L], 40) for L in range(218, 223)}
            else:
                rec['tc'] = {L: mp.nstr(c[L], 40) for L in range(218, 223)}
        out[case] = rec
        print(case, rec, flush=True)
    d = load()
    d['amp_toys'] = out
    save(d)

if __name__ == '__main__':
    {'referee': referee, 'toy': toy, 'timing': timing, 'amp': amp}[sys.argv[1]]()
