"""Controls for a pre-registration record on (1s_A 1s_A | 1s_B 3dsigma_C) (inventory row 8, part 2b; record AD).
Nothing here computes a scored quantity: no exact constant split, no operator, no decay exponent.

Geometry of AB/AP/AQ: B at the origin, C at (0, 0, D1), A at D2 (sin g, 0, cos g), x = cos g.
    PD = (1s_A 1s_A | 1s_B 3ds_C),   3ds_C = N_d h2(r - C) e^{-zC rC},   h2(v) = (3 v_z^2 - |v|^2)/2 = r^2 P2(cos),
    N_d = sqrt(5 (2 zC)^7 / (4 pi 720))   (normalized Slater 3d, n = 3, l = 2, m = 0)

    python debug/aadz_controls.py referee : R2 (AB's Gaussian fits of e^{-r}; h2 e^{-b rC^2} e^{-a rB^2} expanded about
                                            the product centre P on the axis; moments as P-derivatives of the s-type
                                            ERI, f = erf(wR)/R up to second derivatives) at two angles per case; R1
                                            (prolate quadrature against the exact 1s potential); fit floor against
                                            geovac's exact (1s1s|1s 3d0) at R = 2
Output debug/data/aadz_controls.json.
"""
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parents[1]))
import aabc_controls as AB  # noqa: E402
import aapc_controls as AP  # noqa: E402

OUT = Path(__file__).parent / "data" / "aadz_controls.json"
POINTS = {'A1': (0.5, 0.0), 'A2': (-1 / 3, 0.8)}


def Nd(z):
    return np.sqrt(5 * (2 * z) ** 7 / (4 * np.pi * 720))


def gauss_PD(al, c, zA, zB, zC, A, B, C):
    NA2, NB = zA ** 3 / np.pi, np.sqrt(zB ** 3 / np.pi)
    pA = (zA ** 2 * (al[:, None] + al[None, :])).ravel()
    cA = (c[:, None] * c[None, :]).ravel() * NA2
    p, cBC, P = AP._bc_pairs(al, c, zB, zC, B, C)
    cBC = cBC * NB * Nd(zC)
    delta = P[:, 2] - C[2]
    d = P - A
    R = np.linalg.norm(d, axis=-1)
    uz = d[:, 2] / R
    tot = 0.0
    for i in range(len(pA)):
        q = pA[i]
        pref = (np.pi / p) ** 1.5 * (np.pi / q) ** 1.5
        w = np.sqrt(p * q / (p + q))
        f0, f1, f2 = AP._f(R, w)
        dzz = f2 * uz ** 2 + f1 / R * (1 - uz ** 2)
        lap = f2 + 2 * f1 / R
        # 3 u_z^2 - |u|^2 + 4 delta u_z + 2 delta^2, u = r - P;  u_i u_j G = (d_i d_j / 4p^2 + delta_ij / 2p) G
        val = 0.5 * ((3 * dzz - lap) / (4 * p * p) + 4 * delta * f1 * uz / (2 * p) + 2 * delta ** 2 * f0)
        tot += cA[i] * np.sum(cBC * pref * val)
    return tot


def prolate_PD(zA, zB, zC, D1, A, nxi, neta, nphi):
    NB = np.sqrt(zB ** 3 / np.pi)
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
    v = pts - np.array([0.0, 0.0, D1])
    h2 = (3 * v[..., 2] ** 2 - np.sum(v * v, axis=-1)) / 2
    f = np.exp(-lam - h * (zB - zC) * ET) * AB.V_A(s, zA) * h2 * h ** 3 * (XI ** 2 - ET ** 2)
    return NB * Nd(zC) * np.sum(W * f)


def floor_exact():
    """geovac closed form (1s_A 1s_A | 1s_B 3d0_B): Z_A = 1, Z_B = 3 (1s rate 3, hydrogenic 3d = Slater 3d at rate 1)."""
    import sympy as sp
    from geovac.two_center_eri import aabb_closed_form, R_s
    e = aabb_closed_form(Fraction(1), (1, 0, 0), (1, 0, 0), Fraction(3), (1, 0, 0), (3, 2, 0))
    return float(sp.N(e.subs(R_s, sp.Integer(2)), 30))


def referee():
    out = {}
    t0 = time.time()
    fits = {nn: AB.fit(nn) for nn in (22, 26, 30, 34)}  # 34: floor 1.1e-11 (26 alone is 4e-10 for d)
    ex = floor_exact()
    o = np.zeros(3)
    out['floor'] = dict(exact=ex, by_n={str(nn): gauss_PD(al, c, 1.0, 3.0, 1.0, np.array([0, 0, -2.0]), o, o) - ex
                                        for nn, (al, c, _) in fits.items()})
    print('floor', out['floor'], flush=True)
    for case in AB.CASES:
        g = {k: float(Fraction(v)) for k, v in AB.CASES[case].items()}
        B, C = np.zeros(3), np.array([0.0, 0.0, g['D1']])
        for x in POINTS[case]:
            A = g['D2'] * np.array([np.sqrt(1 - x * x), 0.0, x])
            r2 = {str(nn): gauss_PD(al, c, g['zA'], g['zB'], g['zC'], A, B, C) for nn, (al, c, _) in fits.items()}
            r1 = {','.join(map(str, gr)): prolate_PD(g['zA'], g['zB'], g['zC'], g['D1'], A, *gr)
                  for gr in ((60, 60, 72), (90, 90, 108), (120, 120, 144))}
            out[f'{case}_x={x:.6g}'] = dict(x=x, R2=r2, R1=r1)
            print(case, x, 'R2', r2, '\n   R1', r1, round(time.time() - t0), 's', flush=True)
    OUT.write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    {'referee': referee}[sys.argv[1]]()
