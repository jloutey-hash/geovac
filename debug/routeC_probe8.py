"""Probe 8: PANELED sinh-substitution k-grid. theta in [0, theta_max] split into
panels of width W, each with its own GL rule (order n_panel). Fixes the
resolution/reach tension: convergence rate per panel depends on (branch-point
distance pi/2) / (panel half-width), independent of how many panels we chain.
"""
import sys
import time
import mpmath as mp

sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_probe5 import P, std_gl_nodes, outer_grid_plain, outer_grid_sin2, compute_T2


def k_grid_sinh_paneled(theta_max, panel_width, n_panel_deg):
    std_nodes = std_gl_nodes(n_panel_deg, mp.mp.prec)  # reused per panel (same shape)
    nodes = []
    theta_lo = mp.mpf(0)
    while theta_lo < theta_max - mp.mpf('1e-30'):
        theta_hi = min(theta_lo + panel_width, theta_max)
        half = (theta_hi - theta_lo) / 2
        mid = (theta_hi + theta_lo) / 2
        for x, w in std_nodes:
            theta = mid + half * x
            k = 2 * mp.sinh(theta)
            jac = half * 2 * mp.cosh(theta)
            nodes.append((k, w * jac))
        theta_lo = theta_hi
    return nodes


def main():
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    mp.mp.dps = dps
    theta_max = mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf(14)
    panel_width = mp.mpf(sys.argv[3]) if len(sys.argv) > 3 else mp.mpf(2)
    N_deg = int(sys.argv[4]) if len(sys.argv) > 4 else 5
    pdegs = [int(a) for a in sys.argv[5:]] or [3, 4, 5]

    s_list, w_list = outer_grid_sin2(N_deg)
    print(f"dps={dps} theta_max={theta_max} panel_width={panel_width} N={len(s_list)}", flush=True)

    prev = None
    for pdeg in pdegs:
        knodes = k_grid_sinh_paneled(theta_max, panel_width, pdeg)
        kmax = max(kn[0] for kn in knodes)
        t0 = time.time()
        v = compute_T2(s_list, w_list, knodes)
        dt = time.time() - t0
        d = "" if prev is None else mp.nstr(abs(v - prev), 4)
        print(f"  pdeg={pdeg} M={len(knodes):5d} kmax={mp.nstr(kmax,6):>10s}  T2={mp.nstr(v, dps-3)}  diff_prev={d}  ({dt:.2f}s)", flush=True)
        prev = v


if __name__ == '__main__':
    main()
