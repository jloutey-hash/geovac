"""High-precision Legendre coefficients F_l of the three-centre integral (XY|XZ) in the bond angle, for guessing the
angle operator (record AO, docs/predictions/2026-10-28_angle_operator.md).

Same channel formula and the same G as tests/test_paper59_water_bond_angle.py:
    F_l = (2(2l+1)/pi) int_0^inf dk A_l^{XY}(k) A_l^{XZ}(k),   A_l^{XB}(k) = int_0^1 dt j_l(t k D_B) G(k, t; D_B),
in ball arithmetic (python-flint arb).  Fixed double-exponential grids: tanh-sinh in t on [0, 1], exp-sinh in k on
(0, inf), step halved between levels.  j_l by downward recurrence from arb's j_L, j_{L+1} (stable).

    python debug/ao_engine.py <D1> <D2> <zX> <zH> <lmax> <level> <prec_bits> <out.json>
"""
import json
import sys
import time
from multiprocessing import Pool

import flint
from flint import arb

flint.ctx.prec = 64


def nodes_t(h, prec):
    """tanh-sinh on [0, 1]: t = s(pi sinh u), 1 - t = s(-pi sinh u), s = logistic; weight dt/du * h."""
    pi = arb.pi()
    out = []
    j = 0
    while True:
        for u in ((j * h,) if j == 0 else (j * h, -j * h)):
            u = arb(u)
            v = pi * u.sinh()
            t = 1 / (1 + (-v).exp())
            omt = 1 / (1 + v.exp())
            w = h * t * omt * pi * u.cosh()
            out.append((t, omt, w))
        if j > 0 and (float(h.mid()) * j > 1) and abs(float(out[-1][2].mid())) < 2.0 ** (-prec - 10) \
                and abs(float(out[-2][2].mid())) < 2.0 ** (-prec - 10):
            break
        j += 1
    return out


def nodes_k(h, umax):
    """exp-sinh on (0, inf): k = exp(pi/2 sinh u), dk/du = k pi/2 cosh u; u in [-umax, umax]."""
    pi = arb.pi()
    out = []
    n = int(umax / float(h.mid()))
    for j in range(-n, n + 1):
        u = arb(j) * h
        k = (pi / 2 * u.sinh()).exp()
        out.append((k, h * k * pi / 2 * u.cosh()))
    return out


def G(k, t, omt, D, zA, zB, N):
    Dl = (k * k * t * omt + t * zA * zA + omt * zB * zB).sqrt()
    return N * 2 * arb.pi() * t * omt * zA * zB * (-D * Dl).exp() * (D * D * Dl * Dl + 3 * D * Dl + 3) / Dl ** 5


_W = {}


def _setup(D1, D2, zX, zH, lmax, level, prec):
    flint.ctx.prec = prec
    D1, D2, zA, zB = arb(D1), arb(D2), arb(zX), arb(zH)
    N = (zA ** 3 / arb.pi()).sqrt() * (zB ** 3 / arb.pi()).sqrt()
    h = arb(1) / 2 ** (3 + level)
    Ds = [D1] if D1 == D2 else [D1, D2]
    return dict(prec=prec, lmax=lmax, tn=nodes_t(h, prec), kn=nodes_k(h, 4.2), Ds=Ds, zA=zA, zB=zB, N=N)


def _init(args):
    _W.update(_setup(*args))


def A_of_k(ik):
    """A_l(k) for l = 0..lmax at the ik-th k node, each bond length; returned as strings."""
    flint.ctx.prec = _W['prec']
    lmax, tn, Ds, zA, zB, N = _W['lmax'], _W['tn'], _W['Ds'], _W['zA'], _W['zB'], _W['N']
    k = _W['kn'][ik][0]
    nd = int(_W['prec'] * 0.302) + 15
    res = []
    for D in Ds:
        acc = [arb(0)] * (lmax + 1)
        L = lmax + 1
        for (t, omt, w) in tn:
            z = t * k * D
            g = G(k, t, omt, D, zA, zB, N) * w
            sq = (arb.pi() / (2 * z)).sqrt()
            jp = sq * z.bessel_j(arb(2 * L + 3) / 2)
            jc = sq * z.bessel_j(arb(2 * L + 1) / 2)
            js = [None] * (L + 1)
            js[L] = jc
            for l in range(L, 0, -1):
                jm = (2 * l + 1) / z * jc - jp
                js[l - 1] = jm
                jp, jc = jc, jm
            for l in range(lmax + 1):
                acc[l] += js[l] * g
        res.append([a.str(nd, radius=True) for a in acc])
    return res


def run(D1, D2, zX, zH, lmax, level, prec, procs=15):
    args = (D1, D2, zX, zH, lmax, level, prec)
    W = _setup(*args)
    t0 = time.time()
    with Pool(procs, initializer=_init, initargs=(args,)) as p:
        As = p.map(A_of_k, range(len(W['kn'])), chunksize=2)
    F = []
    for l in range(lmax + 1):
        s = arb(0)
        for (k, wk), A in zip(W['kn'], As):
            s += wk * arb(A[0][l]) * arb(A[-1][l])
        F.append(2 * (2 * l + 1) / arb.pi() * s)
    return F, dict(nt=len(W['tn']), nk=len(W['kn']), seconds=time.time() - t0)


if __name__ == '__main__':
    D1, D2, zX, zH = sys.argv[1:5]
    lmax, level, prec = int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
    F, info = run(D1, D2, zX, zH, lmax, level, prec)
    out = dict(args=sys.argv[1:8], info=info, F=[f.str(60, radius=True) for f in F])
    with open(sys.argv[8], 'w') as fh:
        json.dump(out, fh, indent=0)
    print(info, flush=True)
