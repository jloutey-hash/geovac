"""Probe 9: isolate the k-integral convergence for a SINGLE (s,t) pair using the
paneled sinh grid vs. adaptive mp.quad (extended breakpoints for slow decay)."""
import sys
import time
import mpmath as mp

sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_probe5 import P
from routeC_probe8 import k_grid_sinh_paneled


def J_adaptive(s, t):
    b = s + t
    def f(k):
        j0 = mp.sin(k * b) / (k * b) if k * b > mp.mpf('1e-30') else mp.mpf(1)
        return j0 * P(s, k) * P(t, k)
    return mp.quad(f, [0, 1, 3, 8, 20, 100, 500, 3000, 20000, 200000, mp.inf])


def J_fixed(s, t, nodes):
    b = s + t
    tot = mp.mpf(0)
    for k, w in nodes:
        kb = k * b
        j0 = mp.sin(kb) / kb if kb > mp.mpf('1e-40') else mp.mpf(1) - kb * kb / 6
        tot += w * j0 * P(s, k) * P(t, k)
    return tot


def main():
    dps = 60
    mp.mp.dps = dps
    s, t = mp.mpf('3e-7'), mp.mpf('0.5')
    print(f"s={s} t={t}", flush=True)

    t0 = time.time()
    jad = J_adaptive(s, t)
    print(f"adaptive: {mp.nstr(jad, dps-3)}  ({time.time()-t0:.2f}s)", flush=True)

    for pdeg in [3, 4, 5, 6]:
        nodes = k_grid_sinh_paneled(mp.mpf(14), mp.mpf(2), pdeg)
        t0 = time.time()
        jf = J_fixed(s, t, nodes)
        dt = time.time() - t0
        print(f"  pdeg={pdeg} M={len(nodes):5d}  J={mp.nstr(jf, dps-3)}  diff={mp.nstr(abs(jf-jad),4)}  ({dt:.2f}s)", flush=True)


if __name__ == '__main__':
    main()
