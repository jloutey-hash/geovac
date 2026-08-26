"""Probe 10: verify the separable sin/cos-product formula reproduces J(s,t) exactly
for individual pairs, isolating whether the full-2D slow convergence is a bug in
the separable trick or a genuine accumulation effect over many outer pairs."""
import sys
import time
import mpmath as mp

sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_probe5 import P
from routeC_probe8 import k_grid_sinh_paneled
from routeC_probe9 import J_fixed as J_direct


def J_separable(s, t, nodes):
    b = s + t
    tot_A_st = mp.mpf(0)  # sum wp sin(k s)P(s,k) cos(k t)P(t,k)
    tot_A_ts = mp.mpf(0)  # sum wp sin(k t)P(t,k) cos(k s)P(s,k)
    for k, w in nodes:
        wp = w / k
        Ps = P(s, k)
        Pt = P(t, k)
        us = mp.sin(k * s) * Ps
        vs = mp.cos(k * s) * Ps
        ut = mp.sin(k * t) * Pt
        vt = mp.cos(k * t) * Pt
        tot_A_st += wp * us * vt
        tot_A_ts += wp * ut * vs
    return (tot_A_st + tot_A_ts) / b


def main():
    dps = 60
    mp.mp.dps = dps
    pairs = [
        (mp.mpf('3e-7'), mp.mpf('0.5')),
        (mp.mpf('0.5'), mp.mpf('0.5')),
        (mp.mpf('0.3'), mp.mpf('0.7')),
        (mp.mpf('0.01'), mp.mpf('0.02')),
    ]
    for pdeg in [5, 6]:
        nodes = k_grid_sinh_paneled(mp.mpf(14), mp.mpf(2), pdeg)
        print(f"pdeg={pdeg} M={len(nodes)}", flush=True)
        for s, t in pairs:
            jd = J_direct(s, t, nodes)
            js = J_separable(s, t, nodes)
            print(f"  s={mp.nstr(s,4)} t={mp.nstr(t,4)}  J_direct={mp.nstr(jd,dps-3)}  diff_sep={mp.nstr(abs(jd-js),4)}", flush=True)


if __name__ == '__main__':
    main()
