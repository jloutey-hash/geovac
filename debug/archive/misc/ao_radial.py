"""Position-space high-precision Legendre coefficients F_L of (XY|XZ) in the bond angle (record AO,
docs/predictions/2026-10-28_angle_operator.md).  Independent of the momentum route (debug/ao_engine.py).

X at the origin, 1s Slater functions, density rho_Y = phi_X(r) phi_H(|r - Y|), |Y| = D1, |Z| = D2, x = cos(gamma):
    (XY|XZ) = sum_L F_L P_L(x),   F_L = (4 pi/(2L+1))^2 R_L,
    R_L = int int r1^2 r2^2 rho_{Y,L}(r1) rho_{Z,L}(r2) r<^L / r>^(L+1),
    rho_{.,L}(r) = (2L+1)/2 N e^{-zX r} int_{-1}^{1} e^{-zH s} P_L(u) du,   s = |r - D|,
    int e^{-z s} P_L du = -2 d/dz [z i_L(z a) k_L(z b)],  a = min(r, D), b = max(r, D)
(i_0 = sinh x / x, k_0 = e^{-x}/x).  Radial integrals: Chebyshev-Lobatto on pieces split at D1, D2 (the cusp
kinks) and on a geometric tail, n + 2L nodes per piece; cumulative integrals by Chebyshev antiderivative, outer by Clenshaw-Curtis.

    python debug/ao_radial.py <D1> <D2> <zX> <zH> <Lmax> <n_per_piece> <dps> <out.json>
"""
import json
import sys
import time
from multiprocessing import Pool

import mpmath as mp


def il(L, x):
    return mp.sqrt(mp.pi / (2 * x)) * mp.besseli(L + mp.mpf(1) / 2, x)


def kl(L, x):
    return mp.sqrt(2 / (mp.pi * x)) * mp.besselk(L + mp.mpf(1) / 2, x)


def rho_L(L, r, D, zX, zH, N):
    a, b = (r, D) if r < D else (D, r)
    x, y = zH * a, zH * b
    iL, kL = il(L, x), kl(L, y)
    dI = il(L + 1, x) + L / x * iL          # i_L'
    dK = -kl(L + 1, y) + L / y * kL         # k_L'
    ang = -2 * (iL * kL + x * dI * kL + y * iL * dK)
    return (2 * L + 1) / mp.mpf(2) * N * mp.exp(-zX * r) * ang


def cheb_nodes(a, b, n):
    """Chebyshev-Lobatto nodes on [a, b], increasing."""
    return [a + (b - a) * (1 - mp.cos(mp.pi * j / n)) / 2 for j in range(n + 1)]


_COS = {}


def _cos_table(n):
    """cos(pi k j / n) for k = 0..n+1, j = 0..n, at the current precision."""
    key = (n, mp.mp.prec)
    if key not in _COS:
        _COS[key] = [[mp.cos(mp.pi * k * j / n) for j in range(n + 1)] for k in range(n + 2)]
    return _COS[key]


def cheb_cumint(vals, a, b):
    """Values at the Lobatto nodes (xi_j = -cos(pi j/n)) -> cumulative integral from a at each node (spectral)."""
    n = len(vals) - 1
    tab = _cos_table(n)
    wv = [v / 2 if j in (0, n) else v for j, v in enumerate(vals)]
    c = []
    for k in range(n + 1):
        s = mp.fsum(w * t for w, t in zip(wv, tab[k]))
        c.append((-1) ** k * 2 * s / n)
    c[0] /= 2
    c[n] /= 2
    # antiderivative in xi: int T_0 = T_1, int T_1 = T_2/4, int T_k = T_{k+1}/(2(k+1)) - T_{k-1}/(2(k-1))
    C = [mp.mpf(0)] * (n + 2)
    for k in range(n + 1):
        if k == 0:
            C[1] += c[0]
        elif k == 1:
            C[2] += c[1] / 4
        else:
            C[k + 1] += c[k] / (2 * (k + 1))
            C[k - 1] -= c[k] / (2 * (k - 1))
    # T_k(xi_j) = (-1)^k cos(pi k j / n)
    Cs = [(-1) ** k * Ck for k, Ck in enumerate(C)]
    out = [mp.fsum(Cs[k] * tab[k][j] for k in range(n + 2)) for j in range(n + 1)]
    return [(v - out[0]) * (b - a) / 2 for v in out]


def pieces(D1, D2, zsum):
    lo, hi = sorted([D1, D2])
    br = [mp.mpf(0), lo] + ([hi] if hi > lo else [])
    w = mp.mpf(2) / zsum
    end = br[-1]
    tail_end = end + mp.mpf(mp.mp.dps) * mp.log(10) / zsum * 2 + 10
    while end < tail_end:
        end = end + w
        br.append(end)
        w = w * 2
    return [(br[i], br[i + 1]) for i in range(len(br) - 1)]


_A = {}


def _init(a):
    _A.update(a)
    mp.mp.dps = a['dps']


def F_of_L(L):
    mp.mp.dps = _A['dps']
    D1, D2, zX, zH = (mp.mpf(_A[k]) for k in ('D1', 'D2', 'zX', 'zH'))
    n = _A['n'] + 2 * L
    N = mp.sqrt(zX ** 3 / mp.pi) * mp.sqrt(zH ** 3 / mp.pi)
    P = pieces(D1, D2, zX + zH)
    data = []
    for (a, b) in P:
        rs = cheb_nodes(a, b, n)
        if rs[0] == 0:
            rs[0] = mp.mpf(10) ** (-mp.mp.dps)
        g1 = [r * r * rho_L(L, r, D1, zX, zH, N) for r in rs]
        g2 = g1 if D1 == D2 else [r * r * rho_L(L, r, D2, zX, zH, N) for r in rs]
        data.append((a, b, rs, g1, g2))
    # C2(r) = int_0^r g2 r'^L ;  T2(r) = int_r^inf g2 r'^(-L-1)
    acc = mp.mpf(0)
    Cv = []
    for (a, b, rs, g1, g2) in data:
        cum = cheb_cumint([g * r ** L for g, r in zip(g2, rs)], a, b)
        Cv.append([acc + v for v in cum])
        acc += cum[-1]
    totT = mp.mpf(0)
    Tv = [None] * len(data)
    for i in range(len(data) - 1, -1, -1):
        a, b, rs, g1, g2 = data[i]
        cum = cheb_cumint([g / r ** (L + 1) for g, r in zip(g2, rs)], a, b)
        Tv[i] = [totT + cum[-1] - v for v in cum]
        totT += cum[-1]
    R = mp.mpf(0)
    for i, (a, b, rs, g1, g2) in enumerate(data):
        h = [g * (Cv[i][j] / r ** (L + 1) + Tv[i][j] * r ** L) for j, (g, r) in enumerate(zip(g1, rs))]
        R += cheb_cumint(h, a, b)[-1]
    F = (4 * mp.pi / (2 * L + 1)) ** 2 * R
    return mp.nstr(F, _A['dps'] - 5)


def run(D1, D2, zX, zH, Lmax, n, dps, procs=15):
    t0 = time.time()
    with Pool(procs, initializer=_init, initargs=(dict(D1=D1, D2=D2, zX=zX, zH=zH, n=n, dps=dps),)) as p:
        F = p.map(F_of_L, range(Lmax + 1), chunksize=1)
    return F, time.time() - t0


if __name__ == '__main__':
    D1, D2, zX, zH = sys.argv[1:5]
    Lmax, n, dps = int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
    F, sec = run(D1, D2, zX, zH, Lmax, n, dps)
    json.dump(dict(args=sys.argv[1:8], seconds=sec, F=F), open(sys.argv[8], 'w'), indent=0)
    print('seconds', round(sec, 1), flush=True)
