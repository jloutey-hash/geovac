"""Mobius-map outer assembly: T1 (rho^{3/2}-subtracted corner triangle,
UNCHANGED from routeC_corner_v2.py -- already cross-validated to ~14 digits
by raw-vs-subtracted, an independent check in its own right) + T2far/RectA/
RectB computed via the Mobius-mapped outer grid (routeC_mobius_outer.py)
instead of GL/dyadic-GL. Structurally different node family (rational Mobius
map vs polynomial Gauss-Legendre) for the "rest of domain" pieces, as
directed by the coordinator for cross-validating the digit-9 discrepancy
found with the GL-family (dyadic vs single-panel) assembly.
"""
import sys
import time
import mpmath as mp

sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_probe9 import J_fixed
from routeC_probe8 import k_grid_sinh_paneled
from routeC_mobius_outer import mobius_interval_nodes
from routeC_corner_v2 import (corner_triangle_integral_subtracted,
                               far_triangle_integral)


def mobius_rect_integral(deg_s, deg_t, slo, shi, tlo, thi, K0s, K0t, knodes):
    ps = mobius_interval_nodes(slo, shi, K0s, deg_s)
    pt = mobius_interval_nodes(tlo, thi, K0t, deg_t)
    tot = mp.mpf(0)
    for s, ws in ps:
        for t, wt in pt:
            tot += ws * wt * J_fixed(s, t, knodes)
    return tot


def main():
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    mp.mp.dps = dps
    delta = mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf('0.05')
    deg_tri = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    deg_rect = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    K0 = mp.mpf(sys.argv[5]) if len(sys.argv) > 5 else mp.mpf('4')

    knodes = k_grid_sinh_paneled(mp.mpf(16), mp.mpf(2), 5)
    print(f"dps={dps} delta={delta} deg_tri={deg_tri} deg_rect={deg_rect} K0={K0}", flush=True)

    t0 = time.time()
    T1, rem, ab = corner_triangle_integral_subtracted(delta, deg_tri, deg_tri, knodes)
    print(f"  T1 (unchanged, subtracted) = {mp.nstr(T1,20)}  ({time.time()-t0:.1f}s)", flush=True)

    t0 = time.time()
    T2v = far_triangle_integral(delta, deg_tri, deg_tri, knodes)
    print(f"  T2far (unchanged, polar GL) = {mp.nstr(T2v,20)}  ({time.time()-t0:.1f}s)", flush=True)

    t0 = time.time()
    # RectA=[0,delta]x[delta,1]: cluster s near 0 (K0), t plain (K0=1 near delta edge)
    RA = mobius_rect_integral(deg_rect, deg_rect, mp.mpf(0), delta, delta, mp.mpf(1), K0, mp.mpf(1), knodes)
    print(f"  RectA (Mobius) = {mp.nstr(RA,20)}  ({time.time()-t0:.1f}s)", flush=True)

    t0 = time.time()
    # RectB=[delta,1]x[0,1]: cluster s near delta (K0), cluster t near 0 (K0)
    RB = mobius_rect_integral(deg_rect, deg_rect, delta, mp.mpf(1), mp.mpf(0), mp.mpf(1), K0, K0, knodes)
    print(f"  RectB (Mobius) = {mp.nstr(RB,20)}  ({time.time()-t0:.1f}s)", flush=True)

    outer_sum = T1 + T2v + RA + RB
    v = (8 / mp.pi) * outer_sum
    print(f"T2 (Mobius-outer assembly) = {mp.nstr(v, dps-3)}")


if __name__ == '__main__':
    main()
