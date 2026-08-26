"""RAW (no subtraction) J integrated on the SAME (rho,alpha) corner-triangle
grid as routeC_corner_v2.py's T1 -- for cross-validating the subtracted T1
against a treatment that uses NO analytic input at all."""
import sys
import time
import mpmath as mp

sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_probe9 import J_fixed
from routeC_probe8 import k_grid_sinh_paneled
from routeC_corner_v2 import triangle_pts


def raw_triangle_integral(delta, deg_rho, deg_alpha, knodes):
    pts = triangle_pts(delta, deg_rho, deg_alpha, reflect=False)
    tot = mp.mpf(0)
    for s, t, w, rho, alpha in pts:
        tot += w * J_fixed(s, t, knodes)
    return tot


def main():
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    mp.mp.dps = dps
    delta = mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf('0.05')
    degs = [int(a) for a in sys.argv[3:]] or [3, 4, 5]
    knodes = k_grid_sinh_paneled(mp.mpf(16), mp.mpf(2), 5)
    print(f"dps={dps} delta={delta}", flush=True)
    prev = None
    for deg in degs:
        t0 = time.time()
        v = raw_triangle_integral(delta, deg, deg, knodes)
        dt = time.time() - t0
        d = "" if prev is None else mp.nstr(abs(v - prev), 4)
        print(f"  deg={deg}  T1_raw={mp.nstr(v, dps-3)}  diff={d}  ({dt:.1f}s)", flush=True)
        prev = v


if __name__ == '__main__':
    main()
