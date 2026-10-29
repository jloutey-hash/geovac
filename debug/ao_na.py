"""Legendre coefficients in the bond angle of the one-electron three-centre nuclear-attraction integral
V(x) = int phi_X(r) phi_H(|r - Y|) / |r - Z| d^3r,  |Y| = D1, |Z| = D2, x = cos(angle Y X Z)  (record NA, the stepping
stone for AO; same density multipoles rho_{Y,L} as debug/ao_radial2.py):
    V(x) = sum_L F_L P_L(x),   F_L = 4 pi/(2L+1) int_0^inf r^2 rho_{Y,L}(r) r<^L / r>^(L+1) dr,   r< , r> of (r, D2).
Same pieces as debug/ao_radial2.py (kinks at D1 and D2, relative width 2/(L+2)), Clenshaw-Curtis per piece.

    python debug/ao_na.py <D1> <D2> <zX> <zH> <Lmax> <n> <prec_bits> <out.json>
"""
import json
import sys
import time
from multiprocessing import Pool

import flint
from flint import arb

from ao_radial2 import breakpoints, cheb_cumint, cheb_nodes, rho_L

_A = {}


def _init(a):
    _A.update(a)
    flint.ctx.prec = a['prec']


def F_of_L(L):
    flint.ctx.prec = _A['prec']
    D1, D2, zX, zH = (arb(_A[k]) for k in ('D1', 'D2', 'zX', 'zH'))
    n = _A['n']
    N = (zX ** 3 / arb.pi()).sqrt() * (zH ** 3 / arb.pi()).sqrt()
    tiny = arb(2) ** (-_A['prec'])
    R = arb(0)
    for (a, b) in breakpoints(L, D1, D2, zX + zH, _A['prec']):
        rs = cheb_nodes(a, b, n)
        if float(rs[0].mid()) == 0:
            rs[0] = tiny
        h = []
        for r in rs:
            w = r ** L / D2 ** (L + 1) if float(r.mid()) < float(D2.mid()) else D2 ** L / r ** (L + 1)
            h.append(r * r * rho_L(L, r, D1, zX, zH, N) * w)
        R += cheb_cumint(h, a, b)[-1]
    return (4 * arb.pi() / (2 * L + 1) * R).str(int(_A['prec'] * 0.301) + 5, radius=True)


if __name__ == '__main__':
    D1, D2, zX, zH = sys.argv[1:5]
    Lmax, n, prec = int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
    t0 = time.time()
    order = sorted(range(Lmax + 1), key=lambda L: -L)
    with Pool(15, initializer=_init, initargs=(dict(D1=D1, D2=D2, zX=zX, zH=zH, n=n, prec=prec),)) as p:
        res = p.map(F_of_L, order, chunksize=1)
    F = [None] * (Lmax + 1)
    for L, f in zip(order, res):
        F[L] = f
    json.dump(dict(args=sys.argv[1:8], seconds=time.time() - t0, F=F), open(sys.argv[8], 'w'), indent=0)
    print('seconds', round(time.time() - t0, 1), flush=True)
