"""Controls for a pre-registration record on (1s_A 1s_A | 1s_B 2px_C) (inventory row 8, part 2b-i; record AQ).
Nothing here computes a scored quantity: no exact constant split, no operator, no decay exponent.

Geometry of AB/AP: B at the origin, C at (0, 0, D1), A at D2 (sin g, 0, cos g), x = cos g, sin g >= 0.
    PX = (1s_A 1s_A | 1s_B 2px_C),   2px_C = sqrt(zC^5/pi) x e^{-zC rC}   (in the ABC plane, perpendicular to BC)

    python debug/aaqx_controls.py referee  : R2 (Gaussian fits of AB; x e^{-b rC^2} e^{-a rB^2} = K x e^{-p|r-P|^2},
                                             P on the z axis, so the ERI is (1/2p) d/dP_x of the s-type one) and
                                             R1 (prolate quadrature of rho_BC against the exact 1s potential V_A)
    python debug/aaqx_controls.py zero     : sanity: PX = 0 at g = 0 and g = pi (A on the BC axis), both routes
Output debug/data/aaqx_controls.json.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import aabc_controls as AB  # noqa: E402
import aapc_controls as AP  # noqa: E402

OUT = Path(__file__).parent / "data" / "aaqx_controls.json"


def gauss_PX(al, c, zA, zB, zC, A, B, C):
    NA2, NB, NC = zA ** 3 / np.pi, np.sqrt(zB ** 3 / np.pi), np.sqrt(zC ** 5 / np.pi)
    pA = (zA ** 2 * (al[:, None] + al[None, :])).ravel()
    cA = (c[:, None] * c[None, :]).ravel() * NA2
    p, cBC, P = AP._bc_pairs(al, c, zB, zC, B, C)
    cBC = cBC * NB * NC
    d = P - A
    R = np.linalg.norm(d, axis=-1)
    tot = 0.0
    for i in range(len(pA)):
        q = pA[i]
        pref = (np.pi / p) ** 1.5 * (np.pi / q) ** 1.5
        w = np.sqrt(p * q / (p + q))
        f0, f1, _ = AP._f(R, w)
        # x e^{-p|r-P|^2} = (P_x - C_x) e^{..} + (1/2p) d/dP_x e^{..};  P_x = C_x = 0
        val = pref * ((P[:, 0] - C[0]) * f0 + f1 * d[:, 0] / R / (2 * p))
        tot += cA[i] * np.sum(cBC * val)
    return tot


def prolate_PX(zA, zB, zC, D1, A, nxi, neta, nphi):
    NB, NC = np.sqrt(zB ** 3 / np.pi), np.sqrt(zC ** 5 / np.pi)
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
    s = np.linalg.norm(pts - A, axis=-1)
    V = AB.V_A(s, zA) * pts[..., 0]
    f = np.exp(-lam - h * (zB - zC) * ET) * V * h ** 3 * (XI ** 2 - ET ** 2)
    return NB * NC * np.sum(W * f)


def load():
    return json.loads(OUT.read_text()) if OUT.exists() else {}


def referee():
    d = load()
    out = {}
    t0 = time.time()
    fits = {nn: AB.fit(nn) for nn in (18, 22, 26)}
    for case in AB.CASES:
        g, A, B, C = AB.geometry(case)
        r2 = {str(nn): gauss_PX(al, c, g['zA'], g['zB'], g['zC'], A, B, C) for nn, (al, c, _) in fits.items()}
        r1 = {','.join(map(str, gr)): prolate_PX(g['zA'], g['zB'], g['zC'], g['D1'], A, *gr)
              for gr in ((60, 60, 72), (90, 90, 108), (120, 120, 144))}
        out[case] = dict(x=AB.REF_X[case], A=list(map(float, A)), R2=r2, R1=r1)
        print(case, 'R2', r2, '\n   R1', r1, round(time.time() - t0), 's', flush=True)
    d['referee'] = out
    OUT.write_text(json.dumps(d, indent=1))


def zero():
    d = load()
    al, c, _ = AB.fit(22)
    res = {}
    for case in AB.CASES:
        g, A, B, C = AB.geometry(case)
        D2 = float(np.linalg.norm(A))
        for name, Az in (('g0', D2), ('gpi', -D2)):
            Ax = np.array([0.0, 0.0, Az])
            res[f'{case}_{name}'] = dict(R2=gauss_PX(al, c, g['zA'], g['zB'], g['zC'], Ax, B, C),
                                         R1=prolate_PX(g['zA'], g['zB'], g['zC'], g['D1'], Ax, 60, 60, 72))
    print(res)
    d['zero'] = res
    OUT.write_text(json.dumps(d, indent=1))


EXTRA_X = {'A1': 0.0, 'A2': 0.8}   # third referee points (record AQ, P1c)


def extra():
    """R2 (n = 22, 26) and R1 at a second x per case: catches a reduction error that P1 at one point could miss."""
    d = load()
    fits = {nn: AB.fit(nn) for nn in (22, 26)}
    res = {}
    for case in AB.CASES:
        g, _, B, C = AB.geometry(case)
        x = EXTRA_X[case]
        A = g['D2'] * np.array([np.sqrt(1 - x * x), 0.0, x])
        r2 = {str(nn): gauss_PX(al, c, g['zA'], g['zB'], g['zC'], A, B, C) for nn, (al, c, _) in fits.items()}
        r1 = {','.join(map(str, gr)): prolate_PX(g['zA'], g['zB'], g['zC'], g['D1'], A, *gr)
              for gr in ((90, 90, 108), (120, 120, 144))}
        res[case] = dict(x=x, R2=r2, R1=r1)
        print(case, x, r2, r1, flush=True)
    d['extra'] = res
    OUT.write_text(json.dumps(d, indent=1))


def toy():
    """Positive control for P3's instrument: h = (x* - x)^{7/2} (AB's branch), PX-like series
    G_L = c_{L-1}/(2L-1) - c_{L+1}/(2L+3), c the Legendre coefficients of h; beta_eff as in AB.  Also the plain
    Legendre series of h (AB-type) for comparison.  No PX quantity is computed."""
    from fractions import Fraction as Fr
    N = 700
    import mpmath
    mpmath.mp.dps = 120
    X, W = [], []
    # Gauss-Legendre with N nodes at 120 digits (Golub-Welsch via mpmath polyroots is slow; use Newton on P_N)
    for k in range(1, N + 1):
        x = mpmath.cos(mpmath.pi * (k - mpmath.mpf(1) / 4) / (N + mpmath.mpf(1) / 2))
        for _ in range(100):
            p0, p1 = mpmath.mpf(1), x
            for n in range(1, N):
                p0, p1 = p1, ((2 * n + 1) * x * p1 - n * p0) / (n + 1)
            dp = N * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < mpmath.mpf(10) ** -110:
                break
        X.append(x)
        W.append(2 / ((1 - x * x) * dp * dp))
    out = {}
    for case in AB.CASES:
        g = {k: mpmath.mpf(Fr(v).numerator) / Fr(v).denominator for k, v in AB.CASES[case].items()}
        xs_ = (g['D1'] ** 2 + g['D2'] ** 2) / (2 * g['D1'] * g['D2'])
        rho = 1 / (xs_ + mpmath.sqrt(xs_ * xs_ - 1))
        Lmax = 242
        c = [mpmath.mpf(0)] * (Lmax + 2)
        for x, w in zip(X, W):
            hv = w * (xs_ - x) ** mpmath.mpf(3.5)
            p0, p1 = mpmath.mpf(1), x
            c[0] += hv
            c[1] += hv * x
            for n in range(1, Lmax + 1):
                p0, p1 = p1, ((2 * n + 1) * x * p1 - n * p0) / (n + 1)
                c[n + 1] += hv * p1
        c = [(2 * n + 1) / mpmath.mpf(2) * v for n, v in enumerate(c)]
        G = [None] + [c[L - 1] / (2 * L - 1) - c[L + 1] / (2 * L + 3) for L in range(1, Lmax + 1)]

        def beta(F, L):
            a, b = F[L + 1] / rho ** (L + 1), F[L - 1] / rho ** (L - 1)
            return float(mpmath.log(abs(a / b)) / mpmath.log(mpmath.mpf(L + 1) / (L - 1)))
        out[case] = dict(P_L={L: beta(c, L) for L in (40, 160, 220)},
                         P_L1={L: beta(G, L) for L in (40, 160, 220)})
        print(case, out[case], flush=True)
    d = load()
    d['toy_P3'] = out
    OUT.write_text(json.dumps(d, indent=1))


if __name__ == '__main__':
    {'referee': referee, 'zero': zero, 'extra': extra, 'toy': toy}[sys.argv[1]]()
