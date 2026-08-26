"""Probe 12: PANELED outer (s,t) grid -- split [0,1] into sub-intervals (denser
near 0,1) with plain GL on each panel, avoiding the global sin^2 substitution's
degenerate s_min~1/N^4 endpoint clustering. Panels far from the endpoints should
converge fast (no pathology); panels near 0/1 test whether narrow-panel GL beats
the global substitution for the accumulating-branch-point endpoint behavior."""
import sys
import time
import mpmath as mp

sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_probe5 import P, std_gl_nodes
from routeC_probe8 import k_grid_sinh_paneled
from routeC_probe11 import compute_T2_ts as compute_T2_generic


PANEL_BOUNDS = [mp.mpf(x) for x in
                ['0', '0.001', '0.01', '0.05', '0.15', '0.35', '0.5']]
# mirror to get the full [0,1] set: b, 1-b reversed
def full_panel_bounds():
    left = PANEL_BOUNDS
    right = [1 - b for b in reversed(left[:-1])]
    return left + right


def outer_grid_paneled(deg):
    bounds = full_panel_bounds()
    std_nodes = std_gl_nodes(deg, mp.mp.prec)
    s_list = []
    w_list = []
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        half = (hi - lo) / 2
        mid = (hi + lo) / 2
        for x, w in std_nodes:
            s_list.append(mid + half * x)
            w_list.append(half * w)
    return s_list, w_list


def main():
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    mp.mp.dps = dps
    pdeg = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    outer_degs = [int(a) for a in sys.argv[3:]] or [2, 3, 4]

    knodes = k_grid_sinh_paneled(mp.mpf(14), mp.mpf(2), pdeg)
    print(f"dps={dps} k: pdeg={pdeg} M={len(knodes)}  panels={len(full_panel_bounds())-1}", flush=True)

    prev = None
    for od in outer_degs:
        s_list, w_list = outer_grid_paneled(od)
        N = len(s_list)
        t0 = time.time()
        v = compute_T2_generic(s_list, w_list, knodes)
        dt = time.time() - t0
        d = "" if prev is None else mp.nstr(abs(v - prev), 4)
        print(f"  outer_deg={od} N={N:5d}  T2={mp.nstr(v, dps-3)}  diff_prev={d}  ({dt:.2f}s)", flush=True)
        prev = v


if __name__ == '__main__':
    main()
