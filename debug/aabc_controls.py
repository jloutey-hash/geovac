"""Controls for the pre-registration record docs/predictions/2026-11-14_aa_bc_three_centre.md (AB).
None of these computes a scored quantity: no exact constant split of the screened pieces, no operator of (AA|BC).

    python debug/aabc_controls.py referee        : value of (1s_A 1s_A | 1s_B 1s_C) at the two referee points by
                                                   R2 (Gaussian fits, closed-form s-Gaussian ERIs; no use of V_A)
                                                   and R1 (prolate quadrature of rho_BC against V_A), + fit floors
    python debug/aabc_controls.py decay          : NA (case N1) exact coefficients, local decay exponent (positive
                                                   control for the P3b instrument; expected -> -2)
    python debug/aabc_controls.py candidates     : the P2 candidate constant list at the cases (arguments, degeneracy)

Geometry: B at the origin (1s, zeta_B = zX), C at (0, 0, D1) (1s, zeta_C = zH), A at D2 (sin g, 0, cos g)
(1s, zeta_A), x = cos g = cos(angle C B A).  Normalized 1s functions.  Output debug/data/aabc_controls.json.
"""
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import erf
import mpmath as mp

OUT = Path(__file__).parent / "data" / "aabc_controls.json"
CASES = {'A1': dict(D1='9/5', D2='13/5', zB='13/10', zC='1', zA='9/10'),
         'A2': dict(D1='3/2', D2='12/5', zB='1', zC='6/5', zA='4/5')}
REF_X = {'A1': '1/2', 'A2': '-1/3'}


def load():
    return json.loads(OUT.read_text()) if OUT.exists() else {}


def save(d):
    OUT.write_text(json.dumps(d, indent=1))


# ---------------- R2: Gaussian fits of e^{-r}, closed-form s-type ERIs ----------------
def fit(n):
    """Even-tempered least-squares fit e^{-r} ~ sum c_i e^{-a_i r^2} in L2(R^3); returns (alpha, c, residual^2/pi)."""
    mp.mp.dps = 50

    def solve(la, lb):
        a = [mp.e ** (la + i * lb) for i in range(n)]
        S = mp.matrix(n, n)
        b = mp.matrix(n, 1)
        for i in range(n):
            # int e^{-r} e^{-a r^2} d^3r = 4 pi int r^2 e^{-r - a r^2} dr
            sg = 1 / (2 * mp.sqrt(a[i]))
            gg = mp.e ** (sg * sg) * mp.erfc(sg)
            b[i] = 4 * mp.pi * mp.sqrt(mp.pi / a[i]) / 2 * ((2 + 4 * sg * sg) * gg - 4 * sg / mp.sqrt(mp.pi)) / (4 * a[i])
            for j in range(i + 1):
                S[i, j] = S[j, i] = (mp.pi / (a[i] + a[j])) ** 1.5
        c = mp.lu_solve(S, b)
        res = mp.pi - (b.T * c)[0]
        return a, c, res / mp.pi

    def obj(p):
        return float(mp.log10(abs(solve(mp.mpf(p[0]), mp.mpf(p[1]))[2]) + mp.mpf(10) ** -45))

    x0 = [np.log(0.03), np.log(2.6 if n < 16 else 1.9)]
    r = minimize(obj, x0, method='Nelder-Mead', options=dict(xatol=1e-4, fatol=1e-4, maxiter=400))
    a, c, res = solve(mp.mpf(r.x[0]), mp.mpf(r.x[1]))
    return np.array([float(v) for v in a]), np.array([float(v) for v in c]), float(res)


def eri_charge(p, P, q, Q):
    """int int e^{-p|r-P|^2} e^{-q|s-Q|^2} / |r - s|, arrays broadcast; P, Q shape (..., 3)."""
    R = np.linalg.norm(P - Q, axis=-1)
    pref = (np.pi / p) ** 1.5 * (np.pi / q) ** 1.5
    w = np.sqrt(p * q / (p + q))
    small = R < 1e-12
    val = np.where(small, 2 * w / np.sqrt(np.pi), erf(w * R) / np.where(small, 1.0, R))
    return pref * val


def gauss_aabc(al, c, zA, zB, zC, A, B, C):
    """(1s_A 1s_A | 1s_B 1s_C), every 1s Slater replaced by the fit (exponents scale as zeta^2)."""
    NA2, NB, NC = zA ** 3 / np.pi, np.sqrt(zB ** 3 / np.pi), np.sqrt(zC ** 3 / np.pi)
    # A pair: one-centre
    pA = (zA ** 2 * (al[:, None] + al[None, :])).ravel()
    cA = (c[:, None] * c[None, :]).ravel() * NA2
    # BC pair: Gaussian product theorem
    a = zB ** 2 * al[:, None]
    b = zC ** 2 * al[None, :]
    pBC = a + b
    K = np.exp(-a * b / pBC * np.sum((B - C) ** 2))
    P = (a[..., None] * B + b[..., None] * C) / pBC[..., None]
    cBC = (c[:, None] * c[None, :] * K).ravel() * NB * NC
    pBC = pBC.ravel()
    P = P.reshape(-1, 3)
    tot = 0.0
    for i in range(len(pA)):
        tot += cA[i] * np.sum(cBC * eri_charge(pA[i], A, pBC, P))
    return tot


# ---------------- R1: prolate quadrature of rho_BC against V_A (uses the reduction's V_A) ----------------
def V_A(s, z):
    return -np.expm1(-2 * z * s) / s - z * np.exp(-2 * z * s)


def prolate_aabc(zA, zB, zC, D1, A, nxi, neta, nphi):
    NB, NC = np.sqrt(zB ** 3 / np.pi), np.sqrt(zC ** 3 / np.pi)
    h = D1 / 2
    lam = h * (zB + zC)
    u, wu = np.polynomial.laguerre.laggauss(nxi)          # xi = 1 + u/lam, weight e^{-u}
    xi = 1 + u / lam
    eta, we = np.polynomial.legendre.leggauss(neta)
    phi = 2 * np.pi * np.arange(nphi) / nphi
    XI, ET, PH = np.meshgrid(xi, eta, phi, indexing='ij')
    W = (wu / lam)[:, None, None] * we[None, :, None] * (2 * np.pi / nphi)
    rcyl = h * np.sqrt(np.clip((XI ** 2 - 1) * (1 - ET ** 2), 0, None))
    z = h * (1 + XI * ET)
    pts = np.stack([rcyl * np.cos(PH), rcyl * np.sin(PH), z], axis=-1)
    s = np.linalg.norm(pts - A, axis=-1)
    # e^{-zB rB - zC rC} = e^{-lam xi - h (zB - zC) eta};  e^{-lam (xi-1)} = e^{-u} is in the weight
    f = np.exp(-lam - h * (zB - zC) * ET) * V_A(s, zA) * h ** 3 * (XI ** 2 - ET ** 2)
    return NB * NC * np.sum(W * f)


def geometry(case):
    g = {k: float(Fraction(v)) for k, v in CASES[case].items()}
    x = float(Fraction(REF_X[case]))
    B = np.zeros(3)
    C = np.array([0.0, 0.0, g['D1']])
    A = g['D2'] * np.array([np.sqrt(1 - x * x), 0.0, x])
    return g, A, B, C


def referee():
    d = load()
    out = {}
    t0 = time.time()
    fits = {}
    for n in (14, 18, 22, 26):
        al, c, res = fit(n)
        fits[n] = (al, c, res)
        print('fit n', n, 'L2 residual^2/pi', res, round(time.time() - t0), 's', flush=True)
    # fit floors on closed forms: (AA|AA) = 5 z/8; equal-exponent two-centre Coulomb at R = 2
    z, R = 1.0, 2.0
    J2 = 1 / R - np.exp(-2 * z * R) * (1 / R + 11 * z / 8 + 3 * z * z * R / 4 + z ** 3 * R * R / 6)
    floors = {}
    for n, (al, c, res) in fits.items():
        o = np.zeros(3)
        e1 = gauss_aabc(al, c, z, z, z, o, o, o) - 5 * z / 8
        e2 = gauss_aabc(al, c, z, z, z, o, np.array([0, 0, R]), np.array([0, 0, R])) - J2
        floors[n] = dict(fit_res=res, AAAA_err=e1, AABB_err=e2)
        print('floor n', n, floors[n], flush=True)
    out['fit_floors'] = {str(k): v for k, v in floors.items()}
    for case in CASES:
        g, A, B, C = geometry(case)
        r2 = {str(n): gauss_aabc(al, c, g['zA'], g['zB'], g['zC'], A, B, C) for n, (al, c, _) in fits.items()}
        r1 = {}
        for (nx, ne, nph) in ((40, 40, 48), (60, 60, 72), (90, 90, 108), (120, 120, 144)):
            r1[f'{nx},{ne},{nph}'] = prolate_aabc(g['zA'], g['zB'], g['zC'], g['D1'], A, nx, ne, nph)
        out[case] = dict(x=REF_X[case], R2_gauss=r2, R1_prolate=r1)
        print(case, 'R2', r2, '\n   R1', r1, flush=True)
    d['referee'] = out
    save(d)


# ---------------- decay instrument, positive control on NA (case N1) ----------------
def decay():
    sys.path.insert(0, str(Path(__file__).parent))
    import na_exact as E
    D1, D2, zX, zH = '1.8', '2.6', '1.3', '1.0'
    Lmax = 240
    t0 = time.time()
    Q = E.coefficients(D1, D2, zX, zH, Lmax, E.Field(None))
    print('exact coefficients', round(time.time() - t0), 's', flush=True)
    mp.mp.dps = 2000
    xs = (mp.mpf(D1) ** 2 + mp.mpf(D2) ** 2) / (2 * mp.mpf(D1) * mp.mpf(D2))
    rho = 1 / (xs + mp.sqrt(xs * xs - 1))
    F = [sum(mp.mpf(int(v.p)) / int(v.q) * E.value(k, mp) for k, v in q.items()) for q in Q]
    G = [F[L] / rho ** L for L in range(Lmax + 1)]
    beta = {L: float(mp.log(abs(G[L + 1] / G[L - 1])) / mp.log(mp.mpf(L + 1) / (L - 1))) for L in
            range(20, Lmax, 20)}
    signs = sorted(set(int(mp.sign(f)) for f in F[5:]))
    print('rho', mp.nstr(rho, 12), 'beta_eff', beta, 'signs', signs, flush=True)
    d = load()
    d['decay_NA_N1'] = dict(rho=mp.nstr(rho, 20), beta_eff=beta, signs_L_ge_5=signs, Lmax=Lmax)
    save(d)


# ---------------- candidate constant list (symbolic rule, evaluated) ----------------
def candidates():
    d = load()
    out = {}
    for case, g in CASES.items():
        D1, D2, zB, zC, zA = (Fraction(g[k]) for k in ('D1', 'D2', 'zB', 'zC', 'zA'))
        mu = 2 * zA
        assert D1 < D2
        rates_inner = sorted({-zB + s * zC + t * mu for s in (1, -1) for t in (1, -1)})
        rates_mid = sorted({-zB - zC + t * mu for t in (1, -1)})
        rate_out = -zB - zC - mu
        E1_D1 = sorted({-b * D1 for b in rates_inner + rates_mid})
        E1_D2 = sorted({-b * D2 for b in rates_mid + [rate_out]})
        logs = sorted({abs(b) for b in rates_inner})
        zero = [b for b in rates_inner + rates_mid + [rate_out] if b == 0]
        out[case] = dict(mu=str(mu), rates_inner=[str(v) for v in rates_inner], rates_mid=[str(v) for v in rates_mid],
                         rate_out=str(rate_out), E1_args_D1=[str(v) for v in E1_D1],
                         E1_args_D2=[str(v) for v in E1_D2], log_args=[str(v) for v in logs],
                         zero_rates=[str(v) for v in zero],
                         NA_E1_arg=str((zB + zC) * D2))
        print(case, out[case], flush=True)
    d['candidates'] = out
    save(d)


if __name__ == '__main__':
    {'referee': referee, 'decay': decay, 'candidates': candidates}[sys.argv[1]]()
