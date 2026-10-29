"""Position-space Legendre coefficients F_L of (XY|XZ) in the bond angle, ball arithmetic (python-flint arb), for the
large-budget operator search (record AO2).  Same formulas as debug/ao_radial.py:
    F_L = (4 pi/(2L+1))^2 R_L,  R_L = int int r1^2 r2^2 rho_{Y,L}(r1) rho_{Z,L}(r2) r<^L / r>^(L+1),
    rho_{.,L}(r) = (2L+1)/2 N e^{-zX r} (-2) [i_L k_L + x i_L'(x) k_L + y i_L k_L'(y)],  x = zH min(r,D), y = zH max(r,D).
Pieces: fixed n Chebyshev-Lobatto nodes each; around each kink D (the cusps at D1, D2) the pieces shrink
relative width 2/(L+2) everywhere (the r^(+-L) factors), cut where r^(2L) or the tail falls below working
precision; i_L, k_L elementary (debug/ao_bessel.py).

    python debug/ao_radial2.py <D1> <D2> <zX> <zH> <Lmax> <n> <prec_bits> <out.json>
"""
import json
import sys
import time
from multiprocessing import Pool

import flint
from flint import arb

from ao_bessel import il_pair, kl_pair

_COS = {}


def cos_table(n):
    key = (n, flint.ctx.prec)
    if key not in _COS:
        pi = arb.pi()
        _COS[key] = [[(pi * k * j / n).cos() for j in range(n + 1)] for k in range(n + 2)]
    return _COS[key]


def cheb_nodes(a, b, n):
    pi = arb.pi()
    return [a + (b - a) * (1 - (pi * j / n).cos()) / 2 for j in range(n + 1)]


def cheb_cumint(vals, a, b):
    n = len(vals) - 1
    tab = cos_table(n)
    wv = [v / 2 if j in (0, n) else v for j, v in enumerate(vals)]
    c = []
    for k in range(n + 1):
        s = arb(0)
        row = tab[k]
        for w, t in zip(wv, row):
            s += w * t
        c.append((-1) ** k * 2 * s / n)
    c[0] /= 2
    c[n] /= 2
    C = [arb(0)] * (n + 2)
    for k in range(n + 1):
        if k == 0:
            C[1] += c[0]
        elif k == 1:
            C[2] += c[1] / 4
        else:
            C[k + 1] += c[k] / (2 * (k + 1))
            C[k - 1] -= c[k] / (2 * (k - 1))
    Cs = [(-1) ** k * Ck for k, Ck in enumerate(C)]
    out = []
    for j in range(n + 1):
        s = arb(0)
        for k in range(n + 2):
            s += Cs[k] * tab[k][j]
        out.append(s)
    return [(v - out[0]) * (b - a) / 2 for v in out]


def breakpoints(L, D1, D2, zsum, prec):
    """Pieces of relative width w = 2/(L+2) (absolute width capped at 6/zsum), from r_min, where r^(2L+2) falls
    below 2^-(prec+20) relative to the inner kink, to the tail end, where (D>/r)^L e^{-zsum (r - D>)} does; both kinks
    are breakpoints.  A single piece [0, r_min]."""
    lo, hi = (D1, D2) if float(D1.mid()) < float(D2.mid()) else (D2, D1)
    bits = arb(prec + 20) * arb(2).log()
    r_min = lo * (-bits / (2 * L + 3)).exp()
    if float(r_min.mid()) < float((lo / 16).mid()):
        r_min = lo / 16
    w = arb(2) / (L + 2)
    cap = arb(6) / zsum
    pts = [arb(0), r_min]
    r = r_min
    kinks = [lo] + ([hi] if float((hi - lo).mid()) > 1e-12 else [])
    ki = 0
    while True:
        step = r * w
        if float(step.mid()) > float(cap.mid()):
            step = cap
        nxt = r + step
        if ki < len(kinks) and float(nxt.mid()) >= float(kinks[ki].mid()):
            nxt = kinks[ki]
            ki += 1
        pts.append(nxt)
        r = nxt
        if ki == len(kinks) and float((L * (r / hi).log() + zsum * (r - hi)).mid()) > float(bits.mid()) + 5:
            break
    return [(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]


def il(L, x):
    return (arb.pi() / (2 * x)).sqrt() * x.bessel_i(arb(2 * L + 1) / 2)


def kl(L, x):
    return (2 / (arb.pi() * x)).sqrt() * x.bessel_k(arb(2 * L + 1) / 2)


def rho_L(L, r, D, zX, zH, N):
    a, b = (r, D) if float(r.mid()) < float(D.mid()) else (D, r)
    x, y = zH * a, zH * b
    iL, iL1 = il_pair(L, x)
    kL, kL1 = kl_pair(L, y)
    dI = iL1 + L / x * iL
    dK = -kL1 + L / y * kL
    ang = -2 * (iL * kL + x * dI * kL + y * iL * dK)
    return arb(2 * L + 1) / 2 * N * (-zX * r).exp() * ang


_A = {}


def _init(a):
    _A.update(a)
    flint.ctx.prec = a['prec']


def F_of_L(L):
    flint.ctx.prec = _A['prec']
    D1, D2, zX, zH = (arb(_A[k]) for k in ('D1', 'D2', 'zX', 'zH'))
    n = _A['n']
    N = (zX ** 3 / arb.pi()).sqrt() * (zH ** 3 / arb.pi()).sqrt()
    P = breakpoints(L, D1, D2, zX + zH, _A['prec'])
    tiny = arb(2) ** (-_A['prec'])
    data = []
    for (a, b) in P:
        rs = cheb_nodes(a, b, n)
        if float(rs[0].mid()) == 0:
            rs[0] = tiny
        g1 = [r * r * rho_L(L, r, D1, zX, zH, N) for r in rs]
        g2 = g1 if D1 == D2 else [r * r * rho_L(L, r, D2, zX, zH, N) for r in rs]
        data.append((a, b, rs, g1, g2))
    acc = arb(0)
    Cv = []
    for (a, b, rs, g1, g2) in data:
        cum = cheb_cumint([g * r ** L for g, r in zip(g2, rs)], a, b)
        Cv.append([acc + v for v in cum])
        acc += cum[-1]
    totT = arb(0)
    Tv = [None] * len(data)
    for i in range(len(data) - 1, -1, -1):
        a, b, rs, g1, g2 = data[i]
        cum = cheb_cumint([g / r ** (L + 1) for g, r in zip(g2, rs)], a, b)
        Tv[i] = [totT + cum[-1] - v for v in cum]
        totT += cum[-1]
    R = arb(0)
    for i, (a, b, rs, g1, g2) in enumerate(data):
        h = [g * (Cv[i][j] / r ** (L + 1) + Tv[i][j] * r ** L) for j, (g, r) in enumerate(zip(g1, rs))]
        R += cheb_cumint(h, a, b)[-1]
    F = (4 * arb.pi() / (2 * L + 1)) ** 2 * R
    return F.str(int(_A['prec'] * 0.301) + 5, radius=True), len(P)


def run(D1, D2, zX, zH, Lmax, n, prec, procs=15):
    t0 = time.time()
    order = sorted(range(Lmax + 1), key=lambda L: -L)       # big L first: better load balance
    with Pool(procs, initializer=_init, initargs=(dict(D1=D1, D2=D2, zX=zX, zH=zH, n=n, prec=prec),)) as p:
        res = p.map(F_of_L, order, chunksize=1)
    F = [None] * (Lmax + 1)
    npieces = [None] * (Lmax + 1)
    for L, (f, m) in zip(order, res):
        F[L], npieces[L] = f, m
    return F, npieces, time.time() - t0


if __name__ == '__main__':
    D1, D2, zX, zH = sys.argv[1:5]
    Lmax, n, prec = int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
    F, npieces, sec = run(D1, D2, zX, zH, Lmax, n, prec)
    json.dump(dict(args=sys.argv[1:8], seconds=sec, npieces=npieces, F=F), open(sys.argv[8], 'w'), indent=0)
    print('seconds', round(sec, 1), 'pieces at Lmax', npieces[-1], flush=True)
