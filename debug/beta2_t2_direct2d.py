"""INDEPENDENT check of the (KW) identity: evaluate
      F(k) = int_0^1 int_0^1 j0(k(s+t)) P(s,k) P(t,k) ds dt
directly as a 2D quadrature (graded y^6 maps on each half of [0,1] -> 4 blocks), with NO use of
the cos-representation j0(z)=int_0^1 cos(zw)dw, no w-integral, no s<->1-s folding.
Agreement with beta2_t2_kw_core.F_of_k tests the identity AND both quadratures.
"""
import sys; sys.path.insert(0, 'debug')
import mpmath as mp
from _fastgl import fast_gl
from beta2_t2_kw_core import P

_N = {}
def half_nodes(N, p):
    key = (N, p, mp.mp.dps)
    if key in _N: return _N[key]
    xs, ws = fast_gl(N)
    lo_s, lo_w, hi_s, hi_w = [], [], [], []
    for xg, wg in zip(xs, ws):
        y = (xg + 1) / 2
        s = y ** p / 2
        jac = p * y ** (p - 1) / 2
        lo_s.append(s);      lo_w.append((wg / 2) * jac)
        hi_s.append(1 - s);  hi_w.append((wg / 2) * jac)
    _N[key] = (lo_s + hi_s, lo_w + hi_w)
    return _N[key]

def F_direct(k, N, p=6):
    k = mp.mpf(k)
    S, W = half_nodes(N, p)
    Pv = [w * P(s, k) for s, w in zip(S, W)]
    tot = mp.mpf(0)
    for i, si in enumerate(S):
        if Pv[i] == 0: continue
        row = mp.mpf(0)
        for j, sj in enumerate(S):
            z = k * (si + sj)
            row += Pv[j] * (mp.sin(z) / z if z > mp.mpf('1e-40') else mp.mpf(1))
        tot += Pv[i] * row
    return tot

if __name__ == '__main__':
    import time, beta2_t2_kw_core as C
    mp.mp.dps = 60
    for k in ('0.7', '3', '12', '45', '110', '200'):
        t0 = time.time()
        a = F_direct(k, 170, 6)
        b = F_direct(k, 230, 6)
        t1 = time.time()
        kk = mp.mpf(k)
        c = C.F_of_k(kk, 160 + int(0.4*float(k)), 60 + int(0.7*float(k)), 6)
        print(f"k={k:>5}: F_direct={mp.nstr(b,30)}  selfconv {mp.nstr(abs(b-a)/abs(b),3)}"
              f"  vs F_KW rel {mp.nstr(abs(b-c)/abs(b),4)}  ({t1-t0:.0f}s)", flush=True)
