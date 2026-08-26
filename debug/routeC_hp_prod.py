"""PRODUCTION high-precision evaluator for T2 (Paper 59), >=40 digits.

Same math as routeC_hp_final.py (symmetric tiling: T1 corner via sigma^2 Duffy
[spectral] + T2far + 2*RectA[0,delta]x[delta,1] + BigSquare[delta,1]^2, with
cluster-both-ends "sin^2" GL on the smooth pieces), but the two SMOOTH
rectangles use the RANK-1 SEPARABLE k-trick (product-to-sum on j0(k(s+t)))
instead of a fresh k-quadrature per (s,t) pair:

    J(s,t) = (A(s,t)+A(t,s))/(s+t),  A(s,t)=sum_m wp_m U(s,k_m) V(t,k_m)
    U(x,k)=sin(kx)P(x,k),  V(x,k)=cos(kx)P(x,k),  wp_m=w_m/k_m

so U,V are precomputed ONCE per node (O(N*M) transcendentals) and every pair
is O(M) cheap multiply-adds -- the ~1000x fewer exp/sqrt that makes deg>=7 at
M=1536 affordable.  s+t >= delta > 0 on both rectangles (corner excluded), so
the 1/(s+t) is safe.  BigSquare uses s<->t symmetry (diagonal + 2*upper-tri).

The corner T1 (sigma^2 Duffy) and T2far (reflected Duffy) are NOT product grids
in (s,t), so they keep the direct J_fixed evaluation -- they are small and
converge by deg 6-7.

Usage: python routeC_hp_prod.py [dps] [delta] [deg] [kpdeg] [--conv d0] [--parts]
"""
from __future__ import annotations
import sys
import time
import mpmath as mp

sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_probe9 import J_fixed
from routeC_probe8 import k_grid_sinh_paneled
from routeC_hp_evaluator import std_gl_nodes, P
from routeC_hp_composite import T1_sigma2, T2far


def clustered_nodes(deg, lo, hi):
    """GL in u in[0,1] through x=(1-cos(pi u))/2 (cluster BOTH ends), affine to
    [lo,hi]; weight folds dx/du and interval Jacobians."""
    std = std_gl_nodes(deg, mp.mp.prec)
    span = hi - lo
    out = []
    for xu, wu in std:
        u = (xu + 1) / 2
        x = (1 - mp.cos(mp.pi * u)) / 2
        dxdu = (mp.pi / 2) * mp.sin(mp.pi * u)
        out.append((lo + span * x, (wu / 2) * dxdu * span))
    return out


def build_UV(nodes, ks):
    """Per node x: U[m]=sin(k_m x)P(x,k_m), V[m]=cos(k_m x)P(x,k_m)."""
    U, V = [], []
    for x, _ in nodes:
        Ui = []
        Vi = []
        for k in ks:
            Pv = P(x, k)
            kx = k * x
            Ui.append(mp.sin(kx) * Pv)
            Vi.append(mp.cos(kx) * Pv)
        U.append(Ui)
        V.append(Vi)
    return U, V


def rect_rank1(sn, tn, ks, wp, symmetric=False):
    """Integral of J over the product grid sn x tn via the rank-1 k-trick.
    sn, tn: lists of (x, weight). symmetric=True: sn is tn (same grid), use
    s<->t symmetry."""
    Us, Vs = build_UV(sn, ks)
    if symmetric:
        Ut, Vt = Us, Vs
        tn = sn
    else:
        Ut, Vt = build_UV(tn, ks)
    M = len(ks)
    sx = [p[0] for p in sn]; sw = [p[1] for p in sn]
    tx = [p[0] for p in tn]; tw = [p[1] for p in tn]

    def Adot(Ui, Vj):
        acc = mp.mpf(0)
        for m in range(M):
            acc += wp[m] * Ui[m] * Vj[m]
        return acc

    tot = mp.mpf(0)
    if symmetric:
        n = len(sn)
        for i in range(n):
            # diagonal: J(s,s) = (A(s,s)+A(s,s))/(2s) = A(s,s)/s
            Aii = Adot(Us[i], Vs[i])
            tot += sw[i] * sw[i] * (Aii / sx[i])
            for j in range(i):
                Aij = Adot(Us[i], Vs[j])
                Aji = Adot(Us[j], Vs[i])
                tot += 2 * sw[i] * sw[j] * (Aij + Aji) / (sx[i] + sx[j])
    else:
        for i in range(len(sn)):
            for j in range(len(tn)):
                Aij = Adot(Us[i], Vt[j])
                Aji = Adot(Ut[j], Vs[i])
                tot += sw[i] * tw[j] * (Aij + Aji) / (sx[i] + tx[j])
    return tot


def compute_corner(delta, cdeg, kn, verbose=False):
    t0 = time.time(); T1 = T1_sigma2(delta, cdeg, cdeg, kn)
    if verbose: print(f"    T1(cdeg={cdeg})    ={mp.nstr(T1,26)} ({time.time()-t0:.1f}s)", flush=True)
    t0 = time.time(); T2f = T2far(delta, cdeg, cdeg, kn)
    if verbose: print(f"    T2far(cdeg={cdeg}) ={mp.nstr(T2f,26)} ({time.time()-t0:.1f}s)", flush=True)
    return T1 + T2f


def compute_rects(delta, rdeg, ks, wp, verbose=False):
    t0 = time.time()
    sn_A = clustered_nodes(rdeg, mp.mpf(0), delta)
    tn_A = clustered_nodes(rdeg, delta, mp.mpf(1))
    RA = rect_rank1(sn_A, tn_A, ks, wp)
    if verbose: print(f"    RectA(rdeg={rdeg})  ={mp.nstr(RA,26)} ({time.time()-t0:.1f}s)", flush=True)
    t0 = time.time()
    sn_B = clustered_nodes(rdeg, delta, mp.mpf(1))
    RB = rect_rank1(sn_B, sn_B, ks, wp, symmetric=True)
    if verbose: print(f"    BigSq(rdeg={rdeg})  ={mp.nstr(RB,26)} ({time.time()-t0:.1f}s)", flush=True)
    return RA, RB


def assemble(delta, deg, kn, verbose=False):
    ks = [kn_i[0] for kn_i in kn]
    wp = [kn_i[1] / kn_i[0] for kn_i in kn]
    corner = compute_corner(delta, deg, kn, verbose)
    RA, RB = compute_rects(delta, deg, ks, wp, verbose)
    outer = corner + 2 * RA + RB
    return (8 / mp.pi) * outer, {'corner': corner, 'RectA': RA, 'BigSq': RB}


def main():
    # final mode: corner ONCE at cdeg, sweep rectangle degree rd in [rd0..rd1]
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    dps = int(args[0]) if len(args) > 0 else 55
    delta = mp.mpf(args[1]) if len(args) > 1 else mp.mpf('0.05')
    cdeg = int(args[2]) if len(args) > 2 else 7          # corner degree (fixed)
    rd0 = int(args[3]) if len(args) > 3 else 7           # rect degree start
    rd1 = int(args[4]) if len(args) > 4 else 9           # rect degree end
    kpdeg = int(args[5]) if len(args) > 5 else 7         # rect k-grid pdeg (M)
    ckpdeg = int(args[6]) if len(args) > 6 else 6        # corner k-grid pdeg (smaller M)
    mp.mp.dps = dps
    kn = k_grid_sinh_paneled(mp.mpf(16), mp.mpf(2), kpdeg)          # rectangles
    kn_c = k_grid_sinh_paneled(mp.mpf(16), mp.mpf(2), ckpdeg)       # corner
    ks = [k for k, w in kn]; wp = [w / k for k, w in kn]
    print(f"dps={dps} delta={delta}  rect_M={len(kn)}(pdeg{kpdeg})  corner_M={len(kn_c)}(pdeg{ckpdeg})  cdeg={cdeg} rect={rd0}..{rd1}", flush=True)
    anchor = mp.mpf('0.3953557659017139')

    t0 = time.time()
    corner = compute_corner(delta, cdeg, kn_c, verbose=True)
    print(f"    corner(T1+T2far) = {mp.nstr(corner,26)}  ({time.time()-t0:.1f}s)", flush=True)

    prev = None
    for rd in range(rd0, rd1 + 1):
        t0 = time.time()
        RA, RB = compute_rects(delta, rd, ks, wp, verbose=True)
        v = (8 / mp.pi) * (corner + 2 * RA + RB)
        diff = "" if prev is None else mp.nstr(abs(v - prev), 4)
        print(f"  >> rect_deg={rd}: T2={mp.nstr(v,dps-4)}  selfdiff={diff}  vsA16={mp.nstr(abs(v-anchor),4)}  ({time.time()-t0:.1f}s)", flush=True)
        prev = v


if __name__ == '__main__':
    main()
