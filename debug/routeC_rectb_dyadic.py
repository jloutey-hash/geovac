"""Test: does splitting RectB=[delta,1]x[0,1] into dyadic s-panels converge
faster than one big panel at matched total point count?"""
import sys
import time
import mpmath as mp

sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_probe9 import J_fixed
from routeC_probe8 import k_grid_sinh_paneled
from routeC_corner_v2 import rect_integral, rect_pts


def dyadic_rect_integral(delta, deg, knodes, tlo=None, thi=None, n_extra_panels=10):
    """s-panels: [delta,2d],[2d,4d],...,[2^{n-1}d, 1]; each x t=[tlo,thi]."""
    tlo = mp.mpf(0) if tlo is None else tlo
    thi = mp.mpf(1) if thi is None else thi
    bounds = [delta]
    b = delta
    for _ in range(n_extra_panels):
        if b >= 1:
            break
        b2 = min(2 * b, mp.mpf(1))
        bounds.append(b2)
        b = b2
    tot = mp.mpf(0)
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        tot += rect_integral(deg, deg, lo, hi, tlo, thi, knodes)
    return tot, bounds


def main():
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    mp.mp.dps = dps
    delta = mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf('0.05')
    degs = [int(a) for a in sys.argv[3:]] or [3, 4]
    knodes = k_grid_sinh_paneled(mp.mpf(16), mp.mpf(2), 5)

    print("Dyadic-panel RectB:", flush=True)
    prev = None
    for deg in degs:
        t0 = time.time()
        v, bounds = dyadic_rect_integral(delta, deg, knodes)
        dt = time.time() - t0
        d = "" if prev is None else mp.nstr(abs(v - prev), 4)
        print(f"  deg={deg} bounds={[mp.nstr(b,4) for b in bounds]}  v={mp.nstr(v,dps-3)}  diff={d}  ({dt:.1f}s)", flush=True)
        prev = v


if __name__ == '__main__':
    main()
