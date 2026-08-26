"""High-precision evaluator for the Route C collinear integrated three-centre
observable T2 (Paper 59 frontier object).

    T2 = (8/pi) int_0^1 ds int_0^1 dt  J(s,t),   J symmetric in s,t
    J(s,t) = int_0^inf dk  j0(k(s+t)) P(s,k) P(t,k)
    P(x,k) = c e^{-Delta}(1/Delta^3 + 3/Delta^4 + 3/Delta^5),  c=x(1-x), Delta=sqrt(c k^2+1)

Collinear geometry: X=0, Y=(0,0,1), Z=(0,0,-1), 1s zeta=1 => D1=D2=1, |W|=s+t.

METHOD (replaces the slow adaptive-mp.quad-per-(s,t)-pair evaluator in
routeC_fast_evaluator.py):

1. k-INTEGRAL -- product-to-sum separation + a "sinh" branch-point-regularizing
   substitution.
   j0(kb)=sin(kb)/(kb), b=s+t, so
       J(s,t) = (1/b) int_0^inf dk/k [sin(ks)cos(kt) + cos(ks)sin(kt)] P(s,k)P(t,k)
   which SEPARATES the k-integral into two rank-1 pieces in (s,t):
       A(s,t) = int_0^inf (dk/k) sin(ks)P(s,k) . cos(kt)P(t,k)
       J(s,t) = (A(s,t) + A(t,s)) / b
   letting a SHARED k-grid be reused for every (s,t) pair (build U(s,k)=sin(ks)P(s,k)
   and V(s,k)=cos(ks)P(s,k) once per outer node, O(N*M); then every pair costs O(M)
   cheap multiply-adds instead of a fresh adaptive quadrature, O(N^2*M) total).

   For the semi-infinite k-domain: Delta(x,k)=sqrt(c k^2+1) has a branch point at
   k=i/sqrt(c). The substitution k=2*sinh(theta) (i.e. targeting the REFERENCE
   scale c0=1/4, the dominant s=t=1/2 fibre) makes Delta0=cosh(theta) EXACTLY at
   c=c0, and for every c in (0,1/4] the branch point sits at Im(theta)=pi/2 (a
   theta-independent vertical gap) -- this collapses the huge dynamic range needed
   in k (up to ~1e6 for the extreme small-c corners) into a small, FIXED theta
   range (~14), because k=2 sinh(theta) grows exponentially. theta in [0,14] is
   split into PANELS (width ~2) with a Gauss-Legendre rule per panel: each panel's
   convergence rate depends only on (branch-point gap)/(panel half-width), so
   panelling gives a uniformly good, fixed, REUSABLE k-grid. Validated to
   independently reproduce >=60 digits for individual (s,t) pairs against both a
   Mobius-rational-map GL grid and adaptive mp.quad with extended breakpoints
   (see debug/routeC_probe{7,8,9,10}.py).

2. OUTER (s,t) integral -- Gauss-Legendre with the sin^2 substitution
   (s=sin^2(phi)) used by the original evaluator, PLUS the separable-in-k trick
   above so the N^2/2 pair evaluations are cheap.

   DIAGNOSED LIMITATION (see debug/sprint_hp_evaluator_memo.md and
   debug/routeC_probe{11,12,13}.py): the ONLY slow-converging region of the
   (s,t) domain is a small neighbourhood of the (0,0) corner, where b=s+t -> 0
   makes the separable form's 1/b apparently singular. Away from that corner
   (tested: the interior [0.2,0.8]^2, and edge strips near s=0 with t away from
   0), the SAME grid converges to the dps floor by N~96 (<1e-50). At (0,0),
   Asym(s,t)/b is a homogeneous-degree-~1 ratio, and empirically the residual
   there converges only at a "few correct digits per doubling" rate consistent
   with a genuine (D ln D)-type non-analyticity of the kind already established
   for this Bessel-moment family (Paper 59 sprint memo, one-mass slice N(D)).
   A Duffy/polar re-parametrization of the small corner square SOFTENS but does
   not eliminate this (tested, see the memo) in the time available; a full fix
   needs an explicit log-subtraction at the corner and was not attempted here.

   Practical consequence: plain N-doubling of the GLOBAL sin^2+GL outer grid
   gives usable but SLOWLY improving digits (see the convergence table in the
   memo): ~13 digits at N=48, ~20 at N=96, ~23 at N=192, ~27 at N=384 (each
   compared to its predecessor; N=384 is this file's practical ceiling within a
   few-minute budget).

USAGE:
    python routeC_hp_evaluator.py [dps] [outer_deg] [k_pdeg] [theta_max] [panel_w]

    outer_deg: GL degree encoding for the outer sin^2 grid -> N = 3*2^(outer_deg-1)
               (outer_deg=8 -> N=384; outer_deg=7 -> N=192; ...)
    k_pdeg:    GL degree encoding PER PANEL for the k-grid -> n_panel=3*2^(k_pdeg-1)
               (k_pdeg=6 -> 96 pts/panel; with theta_max=14, panel_w=2 -> 7 panels
               -> M=672 total)

Reproduces the given anchor 0.39535576590171392 to all 17 quoted digits, and is
self-convergent (N=192 vs N=384) to ~27-28 digits; see the memo for the full
convergence table and the two independent cross-checks (a different outer
substitution -- plain GL, no sin^2 -- and a different k-grid family -- a
Mobius-rational-map GL grid).
"""
from __future__ import annotations
import sys
import time
import mpmath as mp


def P(x, k):
    c = x * (1 - x)
    Del = mp.sqrt(c * k * k + 1)
    return c * mp.e ** (-Del) * (1 / Del ** 3 + 3 / Del ** 4 + 3 / Del ** 5)


_gl_cache: dict = {}


def std_gl_nodes(degree: int, prec: int):
    """Standard [-1,1] Gauss-Legendre nodes/weights via mpmath's fast built-in
    recurrence-based root finder (n = 3*2^(degree-1) points), cached per
    (degree, prec)."""
    key = (degree, prec)
    if key in _gl_cache:
        return _gl_cache[key]
    rule = mp.calculus.quadrature.GaussLegendre(mp.mp)
    nodes = rule.calc_nodes(degree, prec)
    _gl_cache[key] = nodes
    return nodes


def k_grid_sinh_paneled(theta_max, panel_width, n_panel_deg):
    """Fixed k-grid via k=2*sinh(theta), theta in [0,theta_max] split into
    panels of width panel_width, each with its own GL rule of degree
    n_panel_deg (n=3*2^(n_panel_deg-1) points/panel). Regularizes the
    sqrt(c k^2+1) branch point (fixed Im(theta)=pi/2 gap for c<=1/4) and
    compresses the huge k-range needed for slow-decay (small-c) corners into a
    small theta-range, since k grows exponentially in theta."""
    std_nodes = std_gl_nodes(n_panel_deg, mp.mp.prec)
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


def outer_grid_sin2(degree: int):
    """s = sin^2(phi), phi in [0, pi/2] via GL (n=3*2^(degree-1) points)."""
    std_nodes = std_gl_nodes(degree, mp.mp.prec)
    half = mp.pi / 4
    s_list = []
    w_list = []
    for x, w in std_nodes:
        phi = half * (x + 1)
        s_list.append(mp.sin(phi) ** 2)
        w_list.append(half * w * mp.sin(2 * phi))
    return s_list, w_list


def outer_grid_plain(degree: int):
    """Plain GL directly on [0,1] (no substitution) -- used as an
    independent cross-check of the outer integral (different quadrature
    family from outer_grid_sin2)."""
    std_nodes = std_gl_nodes(degree, mp.mp.prec)
    s_list = [(x + 1) / 2 for x, w in std_nodes]
    w_list = [w / 2 for x, w in std_nodes]
    return s_list, w_list


def compute_T2(s_list, w_list, knodes, verbose: bool = False):
    """Assemble T2 from a shared outer (s,w) grid and a shared k-grid, using
    the sin/cos product-to-sum separation so each of the N(N+1)/2 outer pairs
    costs O(M) cheap multiply-adds (no repeated transcendental quadrature)."""
    N = len(s_list)
    M = len(knodes)
    ks = [kn[0] for kn in knodes]
    wp = [kn[1] / kn[0] for kn in knodes]  # w'_m = w_m / k_m (finite: sin(ks)/k -> s as k->0)

    t0 = time.time()
    U = [[mp.mpf(0)] * M for _ in range(N)]
    V = [[mp.mpf(0)] * M for _ in range(N)]
    for i, si in enumerate(s_list):
        Ui, Vi = U[i], V[i]
        for m in range(M):
            k = ks[m]
            Pv = P(si, k)
            Ui[m] = mp.sin(k * si) * Pv
            Vi[m] = mp.cos(k * si) * Pv
    if verbose:
        print(f"    build U,V ({N}x{M}): {time.time()-t0:.2f}s", flush=True)

    t0 = time.time()
    total = mp.mpf(0)
    for i in range(N):
        Ui, Vi, wi = U[i], V[i], w_list[i]
        Aii = mp.mpf(0)
        for m in range(M):
            Aii += wp[m] * Ui[m] * Vi[m]
        total += wi * wi * (Aii / s_list[i])  # J(s,s) = 2 A(s,s) / (2 s) = A(s,s)/s
        for j in range(i):
            Uj, Vj, wj = U[j], V[j], w_list[j]
            Aij = mp.mpf(0)
            Aji = mp.mpf(0)
            for m in range(M):
                Aij += wp[m] * Ui[m] * Vj[m]
                Aji += wp[m] * Uj[m] * Vi[m]
            total += 2 * wi * wj * (Aij + Aji) / (s_list[i] + s_list[j])
    if verbose:
        print(f"    bilinear sum: {time.time()-t0:.2f}s", flush=True)
    return (8 / mp.pi) * total


def main():
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 75
    outer_deg = int(sys.argv[2]) if len(sys.argv) > 2 else 7  # N=192
    k_pdeg = int(sys.argv[3]) if len(sys.argv) > 3 else 6      # 96 pts/panel
    theta_max = mp.mpf(sys.argv[4]) if len(sys.argv) > 4 else mp.mpf(14)
    panel_w = mp.mpf(sys.argv[5]) if len(sys.argv) > 5 else mp.mpf(2)

    mp.mp.dps = dps
    knodes = k_grid_sinh_paneled(theta_max, panel_w, k_pdeg)
    s_list, w_list = outer_grid_sin2(outer_deg)
    N, M = len(s_list), len(knodes)
    print(f"dps={dps}  outer N={N} (sin^2+GL)  k M={M} (paneled sinh, "
          f"theta_max={theta_max}, panel_w={panel_w})", flush=True)

    t0 = time.time()
    v = compute_T2(s_list, w_list, knodes, verbose=True)
    dt = time.time() - t0
    print(f"T2 = {mp.nstr(v, dps - 3)}")
    print(f"({dt:.2f}s total)")


if __name__ == '__main__':
    main()
