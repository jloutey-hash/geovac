"""Controls for a pre-registration record on (1s_A 1s_A | 1s_B 3d_{m != 0},C) (inventory row 8, part 2b; record AX).
Nothing here computes a scored quantity: no exact constant split, no operator, no decay exponent.

Geometry of AB/AD: B at the origin, C at (0, 0, D1), A at D2 (sin g, 0, cos g), x = cos g.
    DX = (1s_A 1s_A | 1s_B 3dxz_C),     3dxz_C   = N_xz  v_x v_z         e^{-zC rC},  N_xz  = sqrt(15/(4 pi)) R_norm
    DD = (1s_A 1s_A | 1s_B 3dx2y2_C),   3dx2y2_C = N_x2y2 (v_x^2 - v_y^2) e^{-zC rC}, N_x2y2 = sqrt(15/(16 pi)) R_norm
    v = r - C,  R_norm = sqrt((2 zC)^7 / 720)   (normalized real Slater 3d; N_d of AD = sqrt(5/(4 pi)) R_norm)
    3dyz and 3dxy vanish by y -> -y (checked numerically, `vanish`).

R2: AB's even-tempered Gaussian fits of e^{-r}; the BC product Gaussian centred at P (on the z axis); the
quadratic Q(r - C) = Q(u + delta), u = r - P, delta = P - C, expanded: u_i u_j G -> (d_i d_j / 4p^2 + delta_ij / 2p),
u_i G -> d_i / 2p acting on f(|P - A|) = erf(wR)/R.  General (any A, B, C, any traceless symmetric Q), so the
floors below can place the d function off the line of the other centres.

F (fit floor), two exact references, both through the zonal projection: everything except the d function is
axially symmetric about the line through C and the other centre(s), so only Q's zonal part about that line survives,
and for a degree-2 harmonic Q that part is Q(u_hat) |v|^2 P2(v_hat . u_hat).  Hence
    (..|..Q..) = (N_Q / N_d) Q(u_hat) x (..|..3d0 along u_hat..)  exactly.
    F1: B = C at the origin, A off axis  -> geovac aabb_closed_form (1s 1s | 1s 3d0), Z_A = 1, Z_B = 3, R = 2
    F2: A = B at the origin, C off axis  -> geovac hybrid_closed_form (1s 1s | 1s ; 3d0), Z_A = 1, Z_B = 3, R = 2
F2 runs the delta-linear and delta-constant terms (P != C) that F1 and AD's floor never run.

R1: prolate quadrature (foci B, C) against the exact 1s potential V_A, three grids, reported only.

    python debug/aadx_controls.py referee | vanish
Output debug/data/aadx_controls.json.
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

OUT = Path(__file__).parent / "data" / "aadx_controls.json"
POINTS = {'A1': (0.5, 0.0), 'A2': (-1 / 3, 0.8)}
NS = (22, 26, 30, 34)
GRIDS = ((60, 60, 72), (90, 90, 108), (120, 120, 144))

# Q(v) = v^T M v, with the real-harmonic normalization constant relative to R_norm
HARM = {
    'xz': (np.array([[0, 0, .5], [0, 0, 0], [.5, 0, 0]]), np.sqrt(15 / (4 * np.pi))),
    'x2y2': (np.diag([1.0, -1.0, 0.0]), np.sqrt(15 / (16 * np.pi))),
    'yz': (np.array([[0, 0, 0], [0, 0, .5], [0, .5, 0]]), np.sqrt(15 / (4 * np.pi))),
    'xy': (np.array([[0, .5, 0], [.5, 0, 0], [0, 0, 0]]), np.sqrt(15 / (4 * np.pi))),
}
ND_ANG = np.sqrt(5 / (4 * np.pi))


def Rnorm(z):
    return np.sqrt((2 * z) ** 7 / 720)


def gauss_Q(al, c, zA, zB, zC, A, B, C, key):
    M, ang = HARM[key]
    NA2, NB = zA ** 3 / np.pi, np.sqrt(zB ** 3 / np.pi)
    pA = (zA ** 2 * (al[:, None] + al[None, :])).ravel()
    cA = (c[:, None] * c[None, :]).ravel() * NA2
    p, cBC, P = AP._bc_pairs(al, c, zB, zC, B, C)
    cBC = cBC * NB * ang * Rnorm(zC)
    delta = P - C
    d = P - A
    R = np.linalg.norm(d, axis=-1)
    u = d / R[:, None]
    uMu = np.einsum('ki,ij,kj->k', u, M, u)
    trM = np.trace(M)
    dMd = np.einsum('ki,ij,kj->k', delta, M, delta)
    uMdel = np.einsum('ki,ij,kj->k', u, M, delta)
    tot = 0.0
    for i in range(len(pA)):
        q = pA[i]
        pref = (np.pi / p) ** 1.5 * (np.pi / q) ** 1.5
        w = np.sqrt(p * q / (p + q))
        f0, f1, f2 = AP._f(R, w)
        # sum M_ij H_ij, H_ij = f2 u_i u_j + f1/R (delta_ij - u_i u_j)
        MH = f2 * uMu + f1 / R * (trM - uMu)
        val = MH / (4 * p * p) + trM * f0 / (2 * p) + 2 * f1 * uMdel / (2 * p) + dMd * f0
        tot += cA[i] * np.sum(cBC * pref * val)
    return tot


def prolate_Q(zA, zB, zC, B, C, A, nxi, neta, nphi, key):
    """Foci B (origin) and C = (0,0,D1); phi in chunks to keep memory modest."""
    M, ang = HARM[key]
    D1 = C[2]
    NB = np.sqrt(zB ** 3 / np.pi)
    h = D1 / 2
    lam = h * (zB + zC)
    uu, wu = np.polynomial.laguerre.laggauss(nxi)
    xi = 1 + uu / lam
    eta, we = np.polynomial.legendre.leggauss(neta)
    XI, ET = np.meshgrid(xi, eta, indexing='ij')
    W2 = (wu / lam)[:, None] * we[None, :] * (2 * np.pi / nphi)
    rcyl = h * np.sqrt(np.clip((XI ** 2 - 1) * (1 - ET ** 2), 0, None))
    z = h * (1 + XI * ET)
    base = np.exp(-lam - h * (zB - zC) * ET) * h ** 3 * (XI ** 2 - ET ** 2) * W2
    tot = 0.0
    for ph in 2 * np.pi * np.arange(nphi) / nphi:
        pts = np.stack([rcyl * np.cos(ph), rcyl * np.sin(ph), z], axis=-1)
        s = np.linalg.norm(pts - A, axis=-1)
        v = pts - C
        Q = np.einsum('abi,ij,abj->ab', v, M, v)
        tot += np.sum(base * AB.V_A(s, zA) * Q)
    return NB * ang * Rnorm(zC) * tot


def floor_refs():
    import sympy as sp
    from geovac.two_center_eri import aabb_closed_form, hybrid_closed_form, R_s
    s1, d0 = (1, 0, 0), (3, 2, 0)
    e1 = aabb_closed_form(Fraction(1), s1, s1, Fraction(3), s1, d0)
    e2 = hybrid_closed_form(Fraction(1), s1, s1, s1, Fraction(3), d0)
    return (float(sp.N(e1.subs(R_s, sp.Integer(2)), 30)), float(sp.N(sp.re(e2.subs(R_s, sp.Integer(2))), 30)))


def floors(fits):
    E1, E2 = floor_refs()
    o = np.zeros(3)
    dirn = np.array([0.8, 0.0, 0.6])          # off-axis direction, Q(u_hat) != 0 for both harmonics
    out = {'E0_aabb_R2': E1, 'E0_hybrid_R2': E2}
    for key in ('xz', 'x2y2'):
        M, ang = HARM[key]
        # F1: B = C = 0, A = 2 dirn; zonal axis u_hat = dirn (Q even, sign irrelevant)
        ex1 = ang / ND_ANG * (dirn @ M @ dirn) * E1
        # F2: A = B = 0, C = 2 dirn; zonal axis through C toward A
        ex2 = ang / ND_ANG * (dirn @ M @ dirn) * E2
        out[f'F1_{key}'] = dict(exact=ex1, err={str(n): gauss_Q(al, c, 1.0, 3.0, 1.0, 2 * dirn, o, o, key) - ex1
                                                 for n, (al, c, _) in fits.items()})
        out[f'F2_{key}'] = dict(exact=ex2, err={str(n): gauss_Q(al, c, 1.0, 1.0, 1.0, o, o, 2 * dirn, key) - ex2
                                                 for n, (al, c, _) in fits.items()})
        print('floor', key, out[f'F1_{key}'], out[f'F2_{key}'], flush=True)
    return out


def load():
    return json.loads(OUT.read_text()) if OUT.exists() else {}


def referee():
    d = load()
    t0 = time.time()
    fits = {nn: AB.fit(nn) for nn in NS}
    d['fit_time_s'] = round(time.time() - t0, 1)
    d['floor'] = floors(fits)
    d['floor_time_s'] = round(time.time() - t0, 1)
    ref = {}
    for case in AB.CASES:
        g = {k: float(Fraction(v)) for k, v in AB.CASES[case].items()}
        B, C = np.zeros(3), np.array([0.0, 0.0, g['D1']])
        for x in POINTS[case]:
            A = g['D2'] * np.array([np.sqrt(1 - x * x), 0.0, x])
            for key in ('xz', 'x2y2'):
                t1 = time.time()
                r2 = {str(nn): gauss_Q(al, c, g['zA'], g['zB'], g['zC'], A, B, C, key) for nn, (al, c, _) in fits.items()}
                t2 = time.time()
                r1 = {','.join(map(str, gr)): prolate_Q(g['zA'], g['zB'], g['zC'], B, C, A, *gr, key) for gr in GRIDS}
                ref[f'{case}_x={x:.6g}_{key}'] = dict(x=x, R2=r2, R1=r1, t_R2=round(t2 - t1, 1),
                                                     t_R1=round(time.time() - t2, 1))
                print(case, x, key, 'R2', r2, '\n   R1', r1, round(time.time() - t0), 's', flush=True)
    d['referee'] = ref
    d['total_time_s'] = round(time.time() - t0, 1)
    OUT.write_text(json.dumps(d, indent=1))


def vanish():
    d = load()
    al, c, _ = AB.fit(26)
    res = {}
    case, x = 'A1', 0.5
    g = {k: float(Fraction(v)) for k, v in AB.CASES[case].items()}
    B, C = np.zeros(3), np.array([0.0, 0.0, g['D1']])
    A = g['D2'] * np.array([np.sqrt(1 - x * x), 0.0, x])
    for key in ('yz', 'xy'):
        res[key] = dict(R2_n26=gauss_Q(al, c, g['zA'], g['zB'], g['zC'], A, B, C, key),
                        R1_90=prolate_Q(g['zA'], g['zB'], g['zC'], B, C, A, 90, 90, 108, key))
    print(res, flush=True)
    d['vanish_A1_x=0.5'] = res
    OUT.write_text(json.dumps(d, indent=1))


def _gl_nodes(N, dps):
    import mpmath
    mpmath.mp.dps = dps
    X, W = [], []
    tol = mpmath.mpf(10) ** -(dps - 10)
    for k in range(1, N // 2 + 1):           # positive half; symmetric
        x = mpmath.cos(mpmath.pi * (k - mpmath.mpf(1) / 4) / (N + mpmath.mpf(1) / 2))
        for _ in range(100):
            p0, p1 = mpmath.mpf(1), x
            for n in range(1, N):
                p0, p1 = p1, ((2 * n + 1) * x * p1 - n * p0) / (n + 1)
            dp = N * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < tol:
                break
        w = 2 / ((1 - x * x) * dp * dp)
        X += [x, -x]
        W += [w, w]
    return X, W


def _toy_betas(N, dps, Lmax=222):
    """Associated-Legendre coefficients (no Condon-Shortley phase, AQ's convention P_L^m = (1-x^2)^{m/2} d^m P_L)
    of f_m = (1-x^2)^{m/2} (x* - x)^nu, m = 1, 2, nu = 7/2, 9/2, by Gauss-Legendre quadrature of
    int (1-x^2)^m (x*-x)^nu d^m P_L dx; norm int (P_L^m)^2 = 2/(2L+1) (L+m)!/(L-m)!.  beta_eff as in AQ
    (aaqx_controls.toy): log|F_{L+1} rho^-(L+1) / (F_{L-1} rho^-(L-1))| / log((L+1)/(L-1))."""
    import mpmath
    X, W = _gl_nodes(N, dps)
    mpmath.mp.dps = dps
    out = {}
    for case in AB.CASES:
        g = {k: mpmath.mpf(Fraction(v).numerator) / Fraction(v).denominator for k, v in AB.CASES[case].items()}
        xs_ = (g['D1'] ** 2 + g['D2'] ** 2) / (2 * g['D1'] * g['D2'])
        rho = 1 / (xs_ + mpmath.sqrt(xs_ * xs_ - 1))
        acc = {(m, nu): [mpmath.mpf(0)] * (Lmax + 2) for m in (1, 2) for nu in ('7/2', '9/2')}
        for x, w in zip(X, W):
            om = 1 - x * x
            hv = {'7/2': w * (xs_ - x) ** mpmath.mpf(3.5), '9/2': w * (xs_ - x) ** mpmath.mpf(4.5)}
            for m in (1, 2):
                # d^m P_L by the fixed-m recurrence (L-m+1) Q_{L+1} = (2L+1) x Q_L - (L+m) Q_{L-1}
                q0, q1 = mpmath.mpf(0), mpmath.mpf(1) if m == 1 else mpmath.mpf(3)   # Q_{m-1}=0, Q_m
                fac = om ** m
                for nu in ('7/2', '9/2'):
                    acc[(m, nu)][m] += hv[nu] * fac * q1
                a, b = q0, q1
                for L in range(m, Lmax + 1):
                    a, b = b, ((2 * L + 1) * x * b - (L + m) * a) / (L - m + 1)
                    for nu in ('7/2', '9/2'):
                        acc[(m, nu)][L + 1] += hv[nu] * fac * b
        res = {}
        for (m, nu), I in acc.items():
            F = [None] * (Lmax + 2)
            for L in range(m, Lmax + 2):
                F[L] = (2 * L + 1) * mpmath.factorial(L - m) / (2 * mpmath.factorial(L + m)) * I[L]

            def beta(L):
                a_, b_ = F[L + 1] / rho ** (L + 1), F[L - 1] / rho ** (L - 1)
                return float(mpmath.log(abs(a_ / b_)) / mpmath.log(mpmath.mpf(L + 1) / (L - 1)))
            res[f'T{m}_nu={nu}'] = {str(L): beta(L) for L in (40, 160, 220)}
            res[f'T{m}_nu={nu}_signs_L>=20'] = sorted({int(mpmath.sign(F[L])) for L in range(20, Lmax + 2)})
        out[case] = res
    return out


def toy():
    """Decay-instrument positive controls for P_L^1 (T1) and P_L^2 (T2), at two precisions."""
    d = load()
    t0 = time.time()
    lo = _toy_betas(700, 120)
    t1 = time.time()
    hi = _toy_betas(900, 160)
    t2 = time.time()
    aq = json.loads((Path(__file__).parent / 'data' / 'aaqx_controls.json').read_text()).get('toy_P3', {})
    cmp_aq = {c: {L: hi[c]['T1_nu=7/2'][L] - aq[c]['P_L1'][L] for L in ('40', '160', '220')} for c in aq}
    stab = {c: {k: max(abs(hi[c][k][L] - lo[c][k][L]) for L in ('40', '160', '220'))
                for k in hi[c] if not k.endswith('signs_L>=20')} for c in hi}
    d['toy_decay'] = dict(N700_dps120=lo, N900_dps160=hi, T1_minus_AQ_P_L1=cmp_aq, max_diff_two_precisions=stab,
                          time_s=[round(t1 - t0, 1), round(t2 - t1, 1)])
    print(json.dumps(d['toy_decay'], indent=1), flush=True)
    OUT.write_text(json.dumps(d, indent=1))


if __name__ == '__main__':
    {'referee': referee, 'vanish': vanish, 'toy': toy}[sys.argv[1]]()
