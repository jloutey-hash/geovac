"""WA0-WA4 of docs/predictions/2026-10-27_water_angle_singularities.md.

Channel coefficients F_l of the shared-centre three-centre integral (XY|XZ) over 1s Slater orbitals, which are the
Legendre coefficients of the integral as a function of x = cos(gamma). Formula and closed-form G from
tests/test_pairwise_three_center.py (re-implemented, exponent-general).

Usage: python debug/wa_channels.py [assign|control|cases]
"""
import sys
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import spherical_jn


def tanh_sinh(n, h):
    j = np.arange(-n, n + 1) * h
    x = np.tanh(0.5 * np.pi * np.sinh(j))
    w = 0.5 * np.pi * np.cosh(j) / np.cosh(0.5 * np.pi * np.sinh(j)) ** 2 * h
    return 0.5 * (x + 1), 0.5 * w


def G(k, t, D, zA, zB):
    N = np.sqrt(zA ** 3 / np.pi) * np.sqrt(zB ** 3 / np.pi)
    Dl = np.sqrt(k ** 2 * t * (1 - t) + t * zA ** 2 + (1 - t) * zB ** 2)
    return N * 2 * np.pi * t * (1 - t) * zA * zB * np.exp(-D * Dl) * (D ** 2 * Dl ** 2 + 3 * D * Dl + 3) / Dl ** 5


def grids(level):
    nk = 400 * 2 ** level
    nt, ht = 60 * 2 ** level, 0.07 / 2 ** level
    t, wt = tanh_sinh(nt, ht)
    keep = (t > 0) & (t < 1)
    t, wt = t[keep], wt[keep]
    u, wu = leggauss(nk)
    u = 0.5 * (u + 1); wu = 0.5 * wu
    return t, wt, u, wu


def F_coeffs(D1, D2, zX, zH, lmax, level=0, k_scale=3.0, assign='test'):
    """F_l, l = 0..lmax. assign='test': G(.., zA=zX, zB=zH) as the test file (equal zetas there);
    assign='swap': G(.., zA=zH, zB=zX)."""
    t, wt, u, wu = grids(level)
    k = k_scale * u / (1 - u)
    wk = wu * k_scale / (1 - u) ** 2
    zA, zB = (zX, zH) if assign == 'test' else (zH, zX)
    G1 = G(k[:, None], t[None, :], D1, zA, zB)
    G2 = G(k[:, None], t[None, :], D2, zA, zB)
    out = []
    for l in range(lmax + 1):
        A1 = (spherical_jn(l, np.outer(k, t) * D1) * G1) @ wt
        A2 = (spherical_jn(l, np.outer(k, t) * D2) * G2) @ wt
        out.append(2 * (2 * l + 1) / np.pi * np.sum(wk * A1 * A2))
    return np.array(out)


def rho_tilde_direct(kz, D, zX, zB, n=400):
    """On-axis Fourier transform of chi_X chi_B (B at (0,0,D)) by direct 2D quadrature in (r, cos) about X."""
    r, wr = leggauss(n); r = 0.5 * (r + 1); wr = 0.5 * wr
    R = 30.0 * r / (1 - r + 1e-300) ** 0 * 1.0
    R = 25.0 * r; wR = 25.0 * wr
    c, wc = leggauss(n)
    RR, CC = np.meshgrid(R, c, indexing='ij')
    rB = np.sqrt(RR ** 2 + D ** 2 - 2 * RR * D * CC)
    f = (np.sqrt(zX ** 3 / np.pi) * np.exp(-zX * RR)) * (np.sqrt(zB ** 3 / np.pi) * np.exp(-zB * rB))
    integrand = f * np.cos(kz * RR * CC) * RR ** 2 * 2 * np.pi
    return float(np.einsum('i,j,ij->', wR, wc, integrand))


def rho_tilde_feynman(kz, D, zA, zB, level=1):
    t, wt, _, _ = grids(level)
    # on-axis: e^{i t k.B} -> cos(t k D) for the real part
    return float(np.sum(wt * np.cos(t * kz * D) * G(kz, t, D, zA, zB)))


def fit_rate(F, lmin, lmax_use):
    ls = np.arange(lmin, lmax_use + 1)
    y = np.log(np.abs(F[ls]))
    A = np.vstack([np.ones_like(ls, dtype=float), -np.log(ls), ls.astype(float)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return {'a': coef[0], 'p': coef[1], 'r': float(np.exp(coef[2])), 'logr': float(coef[2])}


def alternation_triples(F, ls, fit):
    """Registered WA3 estimator: detrend by the smooth WA1 fit, then on adjacent triples
    eps_l = |q_l - (q_{l-1} + q_{l+1})/2| / 2, q = F / trend (trend fitted to |F|, positive)."""
    l_all = np.arange(len(F))
    trend = np.exp(fit['a'] - fit['p'] * np.log(np.maximum(l_all, 1)) + l_all * fit['logr'])
    q = F / trend
    vals = [abs(q[l] - 0.5 * (q[l - 1] + q[l + 1])) / 2 / abs(np.mean(q[l - 1:l + 2])) for l in ls
            if 1 <= l < len(F) - 1]
    return max(vals) if vals else float('nan')


def alternation(F, ls):
    """|T/S| from F_l = S_l + (-1)^l T_l with S, T locally linear in l: use 4 consecutive l."""
    vals = []
    for l in ls:
        if l + 3 >= len(F):
            break
        f = F[l:l + 4]
        sg = np.array([(-1) ** (l + j) for j in range(4)], dtype=float)
        j = np.arange(4, dtype=float)
        M = np.vstack([np.ones(4), j, sg, sg * j]).T
        c, *_ = np.linalg.lstsq(M, f, rcond=None)
        S, T = c[0] + 1.5 * c[1], c[2] + 1.5 * c[3]
        vals.append(abs(T / S))
    return max(vals) if vals else float('nan')


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'cases'
    if mode == 'assign':
        for (zX, zB) in [(2.25, 1.24)]:
            for kz in (0.7, 2.0):
                d = rho_tilde_direct(kz, 1.809, zX, zB)
                ft = rho_tilde_feynman(kz, 1.809, zX, zB)
                fs = rho_tilde_feynman(kz, 1.809, zB, zX)
                print(f"k={kz}: direct {d:.10f}  G(zA=zX) {ft:.10f}  G(zA=zH) {fs:.10f}")
    elif mode == 'control':
        X = np.zeros(3); Y = np.array([0, 0, 2.0]); Z = np.array([1.5, 0, 0.5])
        D1, D2 = np.linalg.norm(Y), np.linalg.norm(Z)
        x = float(Y @ Z / (D1 * D2))
        F = F_coeffs(D1, D2, 1.0, 1.0, 20)
        val = np.polynomial.legendre.legval(x, F)
        print(f"WA0(a) channel sum l<=20: {val:.12f}  vs REFEREE_NG16 0.204941652260 diff {val-0.204941652260:.2e}"
              f"  vs registered 0.2049417218 diff {val-0.2049417218:.2e}")
        for (a, b) in [(1.809, 2.5), (1.809, 1.0), (1.809, 1.809)]:
            lo, hi = min(a, b), max(a, b)
            Fe = np.array([lo ** l / hi ** (l + 1) for l in range(41)])
            fr = fit_rate(Fe, 15, 40)
            print(f"WA0(b) 1/|Y-Z| D=({a},{b}): r={fr['r']:.8f} (D</D> {lo/hi:.8f}) p={fr['p']:.2e}"
                  f"  alt={alternation(Fe, range(15, 38)):.2e}")
            Fa = np.array([(-1) ** l * lo ** l / hi ** (l + 1) for l in range(41)])
            print(f"WA0(c) 1/|Y+Z| D=({a},{b}): alt={alternation(Fa, range(15, 38)):.2e}")
            fa = fit_rate(Fa, 15, 40)
            print(f"   triples estimator: 1/|Y-Z| {alternation_triples(Fe, range(16, 40), fr):.2e}"
                  f"   1/|Y+Z| {alternation_triples(Fa, range(16, 40), fa):.2e}")
    else:
        res = {}
        cases = [(1.809, 1.809), (1.809, 2.5), (1.809, 1.0)]
        sets = {'E1': (1.0, 1.0), 'E2': (2.25, 1.24)}
        assign = sys.argv[2] if len(sys.argv) > 2 else 'test'
        for en, (zX, zH) in sets.items():
            for (a, b) in cases:
                F0 = F_coeffs(a, b, zX, zH, 40, level=0, assign=assign)
                F1 = F_coeffs(a, b, zX, zH, 40, level=1, assign=assign)
                floor = np.abs(F1 - F0)
                usable = [l for l in range(41) if abs(F1[l]) >= 100 * max(floor[l], 1e-300)]
                key = f"{en}_{a}_{b}"
                res[key] = {'F_level1': F1.tolist(), 'floor': floor.tolist(), 'usable': usable}
                print(key, "usable l:", usable[0] if usable else None, "..", usable[-1] if usable else None,
                      "count", len(usable))
        with open('debug/data/wa_channels.json', 'w') as fh:
            json.dump(res, fh)
