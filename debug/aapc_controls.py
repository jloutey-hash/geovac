"""Controls for the pre-registration record docs/predictions/2026-11-14_aa_bc_p_orbitals.md (AP).
None of these computes a scored quantity: no exact constant split, no operator, no decay exponent of the new integrals.

Two integrals, Slater functions, the AB geometry (B at the origin, C at (0, 0, D1), A at D2 (sin g, 0, cos g),
x = cos g):
    PC = (1s_A 1s_A | 1s_B 2pz_C),      2pz_C = sqrt(zC^5/pi) (z - D1) e^{-zC rC}
    PA = (2ps_A 2ps_A | 1s_B 1s_C),     2ps_A = sqrt(zA^5/pi) (s . n) e^{-zA s},  s = r - A,  n = A / |A|

    python debug/aapc_controls.py referee   : R2 (Gaussian fits; p functions as centre derivatives of s Gaussians,
                                              closed-form s-Gaussian ERIs and their analytic derivatives) and
                                              R1 (prolate quadrature of rho_BC against the exact V_A), + fit floors
                                              against geovac closed forms
    python debug/aapc_controls.py identity  : the exact PA reduction V_A = (n.grad_A)^2 U + W0 (sympy)
Output debug/data/aapc_controls.json.
"""
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.special import erf

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parents[1]))
import aabc_controls as AB  # noqa: E402  (fit, eri_charge, V_A, geometry, CASES, REF_X)

OUT = Path(__file__).parent / "data" / "aapc_controls.json"
SQPI = np.sqrt(np.pi)


def load():
    return json.loads(OUT.read_text()) if OUT.exists() else {}


def save(d):
    OUT.write_text(json.dumps(d, indent=1))


# ---------------- R2 ----------------
def _f(R, w):
    """f = erf(wR)/R and its first two R-derivatives."""
    g = erf(w * R)
    g1 = 2 * w / SQPI * np.exp(-(w * R) ** 2)
    g2 = -2 * w * w * R * g1
    return g / R, g1 / R - g / R ** 2, g2 / R - 2 * g1 / R ** 2 + 2 * g / R ** 3


def _bc_pairs(al, c, zB, zC, B, C):
    a = zB ** 2 * al[:, None]
    b = zC ** 2 * al[None, :]
    p = a + b
    K = np.exp(-a * b / p * np.sum((B - C) ** 2))
    P = (a[..., None] * B + b[..., None] * C) / p[..., None]
    return p.ravel(), (c[:, None] * c[None, :] * K).ravel(), P.reshape(-1, 3)


def gauss_PC(al, c, zA, zB, zC, A, B, C):
    """(1s_A 1s_A | 1s_B pz_C): z_C e^{-a rC^2} e^{-b rB^2} = K [(z - Pz) + (Pz - Cz)] e^{-p |r - P|^2}."""
    NA2, NB, NC = zA ** 3 / np.pi, np.sqrt(zB ** 3 / np.pi), np.sqrt(zC ** 5 / np.pi)
    pA = (zA ** 2 * (al[:, None] + al[None, :])).ravel()
    cA = (c[:, None] * c[None, :]).ravel() * NA2
    p, cBC, P = _bc_pairs(al, c, zB, zC, B, C)
    cBC = cBC * NB * NC
    tot = 0.0
    for i in range(len(pA)):
        q = pA[i]
        pref = (np.pi / p) ** 1.5 * (np.pi / q) ** 1.5
        w = np.sqrt(p * q / (p + q))
        d = P - A
        R = np.linalg.norm(d, axis=-1)
        f0, f1, _ = _f(R, w)
        # (z - Pz) e^{-p|r-P|^2} = (1/2p) d/dPz e^{-p|r-P|^2};  d/dPz f(|P - A|) = f'(R) (Pz - Az)/R
        val = pref * ((P[:, 2] - C[2]) * f0 + f1 * d[:, 2] / R / (2 * p))
        tot += cA[i] * np.sum(cBC * val)
    return tot


def gauss_PA(al, c, zA, zB, zC, A, B, C, n=None):
    """(ps_A ps_A | 1s_B 1s_C): (s.n)^2 e^{-beta s^2} = (d_n^2 + 2 beta) e^{-beta |r - A|^2} / (4 beta^2)."""
    if n is None:
        n = A / np.linalg.norm(A)
    NA2, NB, NC = zA ** 5 / np.pi, np.sqrt(zB ** 3 / np.pi), np.sqrt(zC ** 3 / np.pi)
    bA = (zA ** 2 * (al[:, None] + al[None, :])).ravel()
    cA = (c[:, None] * c[None, :]).ravel() * NA2
    p, cBC, P = _bc_pairs(al, c, zB, zC, B, C)
    cBC = cBC * NB * NC
    d = A - P
    R = np.linalg.norm(d, axis=-1)
    un = d @ n / R
    tot = 0.0
    for i in range(len(bA)):
        be = bA[i]
        pref = (np.pi / p) ** 1.5 * (np.pi / be) ** 1.5
        w = np.sqrt(p * be / (p + be))
        f0, f1, f2 = _f(R, w)
        dnn = f2 * un ** 2 + f1 / R * (1 - un ** 2)
        tot += cA[i] * np.sum(cBC * pref * (dnn + 2 * be * f0) / (4 * be * be))
    return tot


# ---------------- R1 ----------------
def V_2ps(s, cosn, z):
    """Exact potential of (2ps)^2 at zeta = z: V0(s) + V2(s) P2(cos), from `identity`'s sympy integrals."""
    e = np.exp(-2 * z * s)
    V0 = 1 / s - e / s - 1.5 * z * e - z * z * s * e - z ** 3 * s * s * e / 3
    V2 = (3 / (s ** 3 * z * z)) * (-np.expm1(-2 * z * s)) - 6 * e / (s * s * z) - 6 * e / s - 4 * z * e \
        - 2 * z * z * s * e - 2 * z ** 3 * s * s * e / 3
    return V0 + V2 * (1.5 * cosn ** 2 - 0.5)


def prolate(kind, zA, zB, zC, D1, A, nxi, neta, nphi):
    NB = np.sqrt(zB ** 3 / np.pi)
    NC = np.sqrt((zC ** 5 if kind == 'PC' else zC ** 3) / np.pi)
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
    if kind == 'PC':
        V = AB.V_A(s, zA) * (z - D1)
    else:
        n = A / np.linalg.norm(A)
        V = V_2ps(s, (sv @ n) / s, zA)
    f = np.exp(-lam - h * (zB - zC) * ET) * V * h ** 3 * (XI ** 2 - ET ** 2)
    return NB * NC * np.sum(W * f)


def floors_exact():
    """geovac closed forms (hydrogenic 2p0 = Slater 2p at zeta = Z/2) for the fit floors."""
    import sympy as sp
    from geovac.two_center_eri import aabb_closed_form, R_s
    Rv = sp.Integer(2)
    # (1s_A 1s_A | 1s_B 2p0_B), Z_A = 1, Z_B = 2: 1s_B rate 2, 2p_B rate 1; A at -R z relative to B
    pc = aabb_closed_form(Fraction(1), (1, 0, 0), (1, 0, 0), Fraction(2), (1, 0, 0), (2, 1, 0))
    # (2p0_A 2p0_A | 1s_B 1s_B), Z_A = 2 (2p rate 1), Z_B = 1
    pa = aabb_closed_form(Fraction(2), (2, 1, 0), (2, 1, 0), Fraction(1), (1, 0, 0), (1, 0, 0))
    return float(sp.N(pc.subs(R_s, Rv), 30)), float(sp.N(pa.subs(R_s, Rv), 30))


def referee():
    d = load()
    out = {}
    t0 = time.time()
    fits = {}
    for nn in (18, 22, 26):
        fits[nn] = AB.fit(nn)
        print('fit', nn, fits[nn][2], round(time.time() - t0), 's', flush=True)
    ex_pc, ex_pa = floors_exact()
    o = np.zeros(3)
    fl = {}
    for nn, (al, c, _) in fits.items():
        e_pc = gauss_PC(al, c, 1.0, 2.0, 1.0, np.array([0, 0, -2.0]), o, o) - ex_pc
        e_pa = gauss_PA(al, c, 1.0, 1.0, 1.0, np.array([0, 0, 2.0]), o, o, n=np.array([0, 0, 1.0])) - ex_pa
        fl[str(nn)] = dict(PC_floor=e_pc, PA_floor=e_pa)
        print('floor', nn, fl[str(nn)], flush=True)
    out['floors'] = dict(exact_PC=ex_pc, exact_PA=ex_pa, by_n=fl)
    for case in AB.CASES:
        g, A, B, C = AB.geometry(case)
        for kind, fn in (('PC', gauss_PC), ('PA', gauss_PA)):
            r2 = {str(nn): fn(al, c, g['zA'], g['zB'], g['zC'], A, B, C) for nn, (al, c, _) in fits.items()}
            r1 = {}
            for grid in ((60, 60, 72), (90, 90, 108), (120, 120, 144)):
                r1[','.join(map(str, grid))] = prolate(kind, g['zA'], g['zB'], g['zC'], g['D1'], A, *grid)
            out[f'{case}_{kind}'] = dict(x=AB.REF_X[case], R2=r2, R1=r1)
            print(case, kind, 'R2', r2, '\n   R1', r1, round(time.time() - t0), 's', flush=True)
    d['referee'] = out
    save(d)


def identity():
    import sympy as sp
    s, t, z, c = sp.symbols('s t zeta c', positive=True)
    rho0 = z ** 5 / sp.pi * t ** 2 * sp.exp(-2 * z * t) / 3
    V = {}
    for l, r in ((0, rho0), (2, 2 * rho0)):
        inner = sp.integrate(t ** (l + 2) * r, (t, 0, s))
        outer = sp.integrate(t ** (1 - l) * r, (t, s, sp.oo))
        V[l] = sp.expand(sp.simplify(4 * sp.pi / (2 * l + 1) * (s ** (-l - 1) * inner + s ** l * outer)))
    e = sp.exp(-2 * z * s)
    U = 3 / (2 * s * z ** 2) - 3 * e / (2 * s * z ** 2) - 33 * e / (16 * z) - 9 * s * e / 8 - s ** 2 * z * e / 4
    W0 = 1 / s - e / s - 5 * z * e / 4 - s * z ** 2 * e / 2
    lhs = (U.diff(s, 2) - U.diff(s) / s) * c ** 2 + U.diff(s) / s + W0
    ok = sp.simplify(lhs - (V[0] + V[2] * (3 * c ** 2 - 1) / 2)) == 0
    charge = sp.limit(s * V[0], s, sp.oo)
    print('V0 =', V[0], '\nV2 =', V[2], '\nidentity holds:', ok, ' total charge', charge)
    dd = load()
    dd['identity'] = dict(V0=str(V[0]), V2=str(V[2]), U=str(U), W0=str(W0), holds=bool(ok), charge=str(charge))
    save(dd)


if __name__ == '__main__':
    {'referee': referee, 'identity': identity}[sys.argv[1]]()
