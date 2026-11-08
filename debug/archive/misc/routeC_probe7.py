"""Probe 7: sinh-substitution k-grid (k = 2 sinh(theta)) on a FINITE theta domain.
Regularizes the sqrt(c k^2+1) branch point exactly at the reference scale c0=0.25
(the dominant/central s=t=0.5 case); the branch point sits at a fixed distance
pi/2 from the real theta-axis for ALL c in (0, 0.25], giving uniform geometric GL
convergence. k = 2 sinh(theta) grows exponentially, so theta only needs to range
to ~12-15 to reach k ~ 1e5-1e6 (full decay for even the slowest-decaying corners).
"""
import sys
import time
import mpmath as mp

sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_probe5 import P, std_gl_nodes, outer_grid_plain, outer_grid_sin2, compute_T2


def k_grid_sinh(degree, theta_max):
    std_nodes = std_gl_nodes(degree, mp.mp.prec)
    half = theta_max / 2
    nodes = []
    for x, w in std_nodes:
        theta = half * (x + 1)
        k = 2 * mp.sinh(theta)
        jac = theta_max * mp.cosh(theta)
        nodes.append((k, w * jac))
    return nodes


def main():
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    mp.mp.dps = dps
    theta_max = mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf(14)
    N_deg = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    kdegs = [int(a) for a in sys.argv[4:]] or [5, 6, 7, 8]

    s_list, w_list = outer_grid_sin2(N_deg)
    print(f"dps={dps} theta_max={theta_max} N={len(s_list)}", flush=True)

    prev = None
    for kdeg in kdegs:
        knodes = k_grid_sinh(kdeg, theta_max)
        kmax = max(kn[0] for kn in knodes)
        t0 = time.time()
        v = compute_T2(s_list, w_list, knodes)
        dt = time.time() - t0
        d = "" if prev is None else mp.nstr(abs(v - prev), 4)
        print(f"  kdeg={kdeg} M={len(knodes):5d} kmax={mp.nstr(kmax,6):>10s}  T2={mp.nstr(v, dps-3)}  diff_prev={d}  ({dt:.2f}s)", flush=True)
        prev = v


if __name__ == '__main__':
    main()
