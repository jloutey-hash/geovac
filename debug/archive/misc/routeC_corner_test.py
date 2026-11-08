"""Test: does (J - S) converge faster than J alone on the corner square?"""
import sys
import time
import mpmath as mp

sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_probe9 import J_fixed
from routeC_probe8 import k_grid_sinh_paneled
from routeC_corner_asymptotics import S_corner, corner_square_integral_of_S
from routeC_hp_evaluator import std_gl_nodes


def rect_grid(deg, lo, hi):
    std = std_gl_nodes(deg, mp.mp.prec)
    half = (hi - lo) / 2
    mid = (hi + lo) / 2
    xs = [mid + half * x for x, w in std]
    ws = [half * w for x, w in std]
    return xs, ws


def corner_remainder_integral(delta, deg, knodes, q_pdeg=6, q_theta_max=16):
    xs, ws = rect_grid(deg, mp.mpf(0), delta)
    tot = mp.mpf(0)
    for si, wi in zip(xs, ws):
        for ti, wj in zip(xs, ws):
            Jv = J_fixed(si, ti, knodes)
            Sv = S_corner(si, ti, n_panel_deg=q_pdeg, theta_max=mp.mpf(q_theta_max))
            tot += wi * wj * (Jv - Sv)
    return tot


def corner_raw_integral(delta, deg, knodes):
    xs, ws = rect_grid(deg, mp.mpf(0), delta)
    tot = mp.mpf(0)
    for si, wi in zip(xs, ws):
        for ti, wj in zip(xs, ws):
            tot += wi * wj * J_fixed(si, ti, knodes)
    return tot


def main():
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    mp.mp.dps = dps
    delta = mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf('0.05')
    degs = [int(a) for a in sys.argv[3:]] or [2, 3, 4]

    knodes = k_grid_sinh_paneled(mp.mpf(16), mp.mpf(2), 5)
    print(f"dps={dps} delta={delta}", flush=True)

    print("RAW J (no subtraction):")
    prev = None
    for deg in degs:
        t0 = time.time()
        v = corner_raw_integral(delta, deg, knodes)
        dt = time.time() - t0
        d = "" if prev is None else mp.nstr(abs(v - prev), 4)
        n = len(std_gl_nodes(deg, mp.mp.prec))
        print(f"  deg={deg} N={n}x{n}  v={mp.nstr(v, dps-3)}  diff={d}  ({dt:.2f}s)", flush=True)
        prev = v

    print("REMAINDER (J - S), plus closed-form add-back:")
    addback = corner_square_integral_of_S(delta, alpha_deg=6, n_panel_deg=6, theta_max=mp.mpf(16))
    print(f"  add-back = {mp.nstr(addback, dps-3)}", flush=True)
    prev = None
    for deg in degs:
        t0 = time.time()
        rem = corner_remainder_integral(delta, deg, knodes)
        v = rem + addback
        dt = time.time() - t0
        d = "" if prev is None else mp.nstr(abs(v - prev), 4)
        n = len(std_gl_nodes(deg, mp.mp.prec))
        print(f"  deg={deg} N={n}x{n}  remainder={mp.nstr(rem, dps-3)}  total={mp.nstr(v, dps-3)}  diff={d}  ({dt:.2f}s)", flush=True)
        prev = v


if __name__ == '__main__':
    main()
