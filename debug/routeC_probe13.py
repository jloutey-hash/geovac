"""Probe 13: composite (Duffy-corner + 2 rectangles) outer quadrature.
Diagnosis (probes 9-12): individual (s,t) pairs and the domain INTERIOR converge
spectrally; only the (0,0) corner (where b=s+t -> 0, the separable-trick's 1/b
becomes 0/0) converges slowly (a homogeneous-degree-1 conical non-smoothness).
Fix: Duffy-transform the small square [0,delta]^2 into two triangles via s=rho,
t=rho*v (and its s<->t mirror); this makes the ds dt / b measure exactly
cancel to sin/smooth form. The rest of [0,1]^2 (away from the corner) is tiled
by two ordinary rectangles.
"""
import sys
import time
import mpmath as mp

sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_probe5 import P, std_gl_nodes
from routeC_probe8 import k_grid_sinh_paneled


def J_pair(s, t, ks, wp):
    """Direct per-pair k-integral via the separable sin/cos trick (no 1/b)."""
    b = s + t
    A = mp.mpf(0)
    Bv = mp.mpf(0)
    for k, wpm in zip(ks, wp):
        Ps = P(s, k)
        Pt = P(t, k)
        A += wpm * mp.sin(k * s) * Ps * mp.cos(k * t) * Pt
        Bv += wpm * mp.sin(k * t) * Pt * mp.cos(k * s) * Ps
    return (A + Bv) / b


def corner_triangle_points(delta, deg_rho, deg_v):
    """Lower triangle t<=s of [0,delta]^2 via s=rho, t=rho*v, rho in [0,delta], v in [0,1].
    ds dt = rho drho dv;  J = Asym(rho,rho v)/(rho(1+v)); combined measure:
        J ds dt = [Asym(rho,rho v)/(1+v)] drho dv   (smooth; rho-singularity cancels).
    Returns list of (s,t,weight) triples for this triangle's contribution ALONE
    (weight already includes the drho dv Jacobian; caller still must evaluate
    Asym/(1+v) using J_pair on the ACTUAL (s,t) -- simplest is to just call J_pair
    normally on (s,t) and multiply by the ORIGINAL ds dt weight = rho*w_rho*w_v;
    this avoids re-deriving Asym separately and stays bit-for-bit consistent with
    J_pair's own 1/b normalization (b=rho(1+v)).
    """
    std_rho = std_gl_nodes(deg_rho, mp.mp.prec)
    std_v = std_gl_nodes(deg_v, mp.mp.prec)
    halfr = delta / 2
    pts = []
    for xr, wr in std_rho:
        rho = halfr * (xr + 1)
        wrho = halfr * wr
        for xv, wv in std_v:
            v = (xv + 1) / 2
            wvv = wv / 2
            s = rho
            t = rho * v
            jac = rho  # ds dt = rho drho dv
            w = wrho * wvv * jac
            pts.append((s, t, w))
    return pts


def rect_points(deg_s, deg_t, slo, shi, tlo, thi):
    std_s = std_gl_nodes(deg_s, mp.mp.prec)
    std_t = std_gl_nodes(deg_t, mp.mp.prec)
    halfs = (shi - slo) / 2
    mids = (shi + slo) / 2
    halft = (thi - tlo) / 2
    midt = (thi + tlo) / 2
    pts = []
    for xs, ws in std_s:
        s = mids + halfs * xs
        wS = halfs * ws
        for xt, wt in std_t:
            t = midt + halft * xt
            wT = halft * wt
            pts.append((s, t, wS * wT))
    return pts


def build_all_points(delta, deg_corner, deg_rect):
    pts = []
    # corner: lower triangle (t<=s) doubled for upper (s<=t) by s<->t symmetry
    tri = corner_triangle_points(delta, deg_corner, deg_corner)
    for s, t, w in tri:
        pts.append((s, t, w))
        if s != t:
            pts.append((t, s, w))  # mirror (upper triangle), J symmetric
    # rectangle A: [0,delta] x [delta,1]
    pts.extend(rect_points(deg_rect, deg_rect, mp.mpf(0), delta, delta, mp.mpf(1)))
    # rectangle B: [delta,1] x [0,1] (mirror of A is included here for t in [0,delta] via full [0,1])
    pts.extend(rect_points(deg_rect, deg_rect, delta, mp.mpf(1), mp.mpf(0), mp.mpf(1)))
    return pts


def compute_T2_composite(pts, knodes):
    ks = [kn[0] for kn in knodes]
    wp = [kn[1] / kn[0] for kn in knodes]
    tot = mp.mpf(0)
    for s, t, w in pts:
        tot += w * J_pair(s, t, ks, wp)
    return (8 / mp.pi) * tot


def main():
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    mp.mp.dps = dps
    delta = mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf('0.1')
    pdeg = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    degs = [int(a) for a in sys.argv[4:]] or [3, 4, 5]

    knodes = k_grid_sinh_paneled(mp.mpf(14), mp.mpf(2), pdeg)
    print(f"dps={dps} delta={delta} k:pdeg={pdeg} M={len(knodes)}", flush=True)

    prev = None
    for deg in degs:
        pts = build_all_points(delta, deg, deg)
        t0 = time.time()
        v = compute_T2_composite(pts, knodes)
        dt = time.time() - t0
        d = "" if prev is None else mp.nstr(abs(v - prev), 4)
        print(f"  deg={deg} npts={len(pts):6d}  T2={mp.nstr(v, dps-3)}  diff_prev={d}  ({dt:.2f}s)", flush=True)
        prev = v


if __name__ == '__main__':
    main()
