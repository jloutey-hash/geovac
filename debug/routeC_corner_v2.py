"""Corrected corner-subtracted T2 evaluator.

BUG found and fixed vs the first attempt (routeC_corner_test.py): the
closed-form add-back int_{rho<=delta} S rho drho dalpha covers the TRIANGLE
{s+t<=delta}, not the SQUARE [0,delta]^2 -- comparing it against a remainder
computed on the square silently mismatched domains (raw J on the square gave
~8.6e-8, subtracted-total on the triangle gave ~1.78e-8: NOT the same
integral). Fixed by tiling [0,1]^2 EXACTLY as:
  T1 = {s,t>=0, s+t<=delta}            -- corner triangle, subtraction used
  T2 = {s,t in [0,delta], s+t>delta}    -- far triangle in the small square,
                                           reflected-Duffy grid, NO subtraction
                                           needed (bounded away from origin)
  RectA = [0,delta] x [delta,1]
  RectB = [delta,1] x [0,1]
"""
from __future__ import annotations
import sys
import time
import mpmath as mp

sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_probe9 import J_fixed
from routeC_probe8 import k_grid_sinh_paneled
from routeC_corner_asymptotics import W_of_alpha, H_of_alpha
from routeC_hp_evaluator import std_gl_nodes


def rect_pts(deg_s, deg_t, slo, shi, tlo, thi):
    std_s = std_gl_nodes(deg_s, mp.mp.prec)
    std_t = std_gl_nodes(deg_t, mp.mp.prec)
    halfs, mids = (shi - slo) / 2, (shi + slo) / 2
    halft, midt = (thi - tlo) / 2, (thi + tlo) / 2
    pts = []
    for xs, ws in std_s:
        s = mids + halfs * xs
        for xt, wt in std_t:
            t = midt + halft * xt
            pts.append((s, t, halfs * ws * halft * wt))
    return pts


def triangle_pts(delta, deg_rho, deg_alpha, reflect=False):
    """rho in [0,delta], alpha in [0,1]; s=rho*alpha, t=rho*(1-alpha)
    (reflect=False, the ORIGIN triangle) or s=delta-rho*alpha,
    t=delta-rho*(1-alpha) (reflect=True, the FAR triangle in [0,delta]^2)."""
    std_r = std_gl_nodes(deg_rho, mp.mp.prec)
    std_a = std_gl_nodes(deg_alpha, mp.mp.prec)
    halfr = delta / 2
    pts = []
    for xr, wr in std_r:
        rho = halfr * (xr + 1)
        wrho = halfr * wr
        for xa, wa in std_a:
            alpha = (xa + 1) / 2
            wal = wa / 2
            if reflect:
                s = delta - rho * alpha
                t = delta - rho * (1 - alpha)
            else:
                s = rho * alpha
                t = rho * (1 - alpha)
            jac = rho  # ds dt = rho drho dalpha (sign irrelevant, magnitude only)
            pts.append((s, t, wrho * wal * jac, rho, alpha))
    return pts


def S_val(s, t, w_cache, **kw):
    rho = s + t
    if rho <= 0:
        return mp.mpf(0)
    alpha = s / rho
    key = kw.get('_alpha_round')
    if key is not None:
        akey = mp.nstr(alpha, 25)
    else:
        akey = float(alpha)
    if akey in w_cache:
        Wv = w_cache[akey]
    else:
        Wv = W_of_alpha(alpha, n_panel_deg=kw.get('n_panel_deg', 6),
                         theta_max=kw.get('theta_max', mp.mpf(16)))
        w_cache[akey] = Wv
    return s * t * Wv / mp.sqrt(rho)


def corner_triangle_integral_subtracted(delta, deg_rho, deg_alpha, knodes, **kw):
    """T1: remainder (J-S) on the (rho,alpha) grid + closed-form add-back of S."""
    pts = triangle_pts(delta, deg_rho, deg_alpha, reflect=False)
    w_cache = {}
    rem = mp.mpf(0)
    for s, t, w, rho, alpha in pts:
        Jv = J_fixed(s, t, knodes)
        Sv = S_val(s, t, w_cache, **kw)
        rem += w * (Jv - Sv)
    # closed-form add-back: int_0^delta rho^{5/2} drho . int_0^1 H(alpha) dalpha
    std_a = std_gl_nodes(deg_alpha, mp.mp.prec)
    addback_alpha = mp.mpf(0)
    for xa, wa in std_a:
        alpha = (xa + 1) / 2
        addback_alpha += (wa / 2) * H_of_alpha(alpha, n_panel_deg=kw.get('n_panel_deg', 6),
                                                theta_max=kw.get('theta_max', mp.mpf(16)))
    addback = (mp.mpf(2) / 7) * delta ** mp.mpf('3.5') * addback_alpha
    return rem + addback, rem, addback


def far_triangle_integral(delta, deg_rho, deg_alpha, knodes):
    """T2: no subtraction needed (bounded away from origin)."""
    pts = triangle_pts(delta, deg_rho, deg_alpha, reflect=True)
    tot = mp.mpf(0)
    for s, t, w, rho, alpha in pts:
        tot += w * J_fixed(s, t, knodes)
    return tot


def rect_integral(deg_s, deg_t, slo, shi, tlo, thi, knodes):
    pts = rect_pts(deg_s, deg_t, slo, shi, tlo, thi)
    tot = mp.mpf(0)
    for s, t, w in pts:
        tot += w * J_fixed(s, t, knodes)
    return tot


def compute_T2_corner_subtracted(delta, deg_tri, deg_rect, knodes, verbose=False):
    t0 = time.time()
    T1, rem, ab = corner_triangle_integral_subtracted(delta, deg_tri, deg_tri, knodes)
    if verbose:
        print(f"    T1 (corner, subtracted): rem={mp.nstr(rem,15)} addback={mp.nstr(ab,15)} "
              f"T1={mp.nstr(T1,15)}  ({time.time()-t0:.1f}s)", flush=True)
    t0 = time.time()
    T2v = far_triangle_integral(delta, deg_tri, deg_tri, knodes)
    if verbose:
        print(f"    T2 (far tri, no sub): {mp.nstr(T2v,15)}  ({time.time()-t0:.1f}s)", flush=True)
    t0 = time.time()
    RA = rect_integral(deg_rect, deg_rect, mp.mpf(0), delta, delta, mp.mpf(1), knodes)
    if verbose:
        print(f"    RectA [0,d]x[d,1]: {mp.nstr(RA,15)}  ({time.time()-t0:.1f}s)", flush=True)
    t0 = time.time()
    RB = rect_integral(deg_rect, deg_rect, delta, mp.mpf(1), mp.mpf(0), mp.mpf(1), knodes)
    if verbose:
        print(f"    RectB [d,1]x[0,1]: {mp.nstr(RB,15)}  ({time.time()-t0:.1f}s)", flush=True)
    outer_sum = T1 + T2v + RA + RB
    return (8 / mp.pi) * outer_sum


def main():
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    mp.mp.dps = dps
    delta = mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf('0.05')
    deg_tri = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    deg_rect = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    kpdeg = int(sys.argv[5]) if len(sys.argv) > 5 else 6

    knodes = k_grid_sinh_paneled(mp.mpf(16), mp.mpf(2), kpdeg)
    print(f"dps={dps} delta={delta} deg_tri={deg_tri} deg_rect={deg_rect} k_M={len(knodes)}", flush=True)
    t0 = time.time()
    v = compute_T2_corner_subtracted(delta, deg_tri, deg_rect, knodes, verbose=True)
    dt = time.time() - t0
    print(f"T2 = {mp.nstr(v, dps-3)}  ({dt:.1f}s total)")


if __name__ == '__main__':
    main()
