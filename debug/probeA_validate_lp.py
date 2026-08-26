"""Probe A -- validation of the LP minimax core against exact/brute-force answers."""
from __future__ import annotations

import numpy as np
from scipy.optimize import linprog

import probeA_specaware_lp as core


def abs_minimax(x, f, d):
    """min over deg<=d of max|p(x)-f(x)| (no boundedness) -- for exact checks."""
    A = core.cheb_design(x, d)
    n, nv = A.shape[0], d + 2
    A_ub = np.zeros((2 * n, nv)); b_ub = np.zeros(2 * n)
    A_ub[:n, :d + 1] = A;  A_ub[:n, -1] = -1.0;  b_ub[:n] = f
    A_ub[n:, :d + 1] = -A; A_ub[n:, -1] = -1.0;  b_ub[n:] = -f
    c = np.zeros(nv); c[-1] = 1.0
    res = linprog(c, A_ub=A_ub, b_ub=b_ub,
                  bounds=[(None, None)] * (d + 1) + [(0, None)], method="highs")
    return float(res.x[-1]), np.asarray(res.x[:d + 1])


print("=" * 74)
print("V1  Chebyshev exact minimax:  min_{deg<=k-1} ||x^k - q||_inf = 2^{1-k}")
m = 4001
xg = np.cos(np.pi * np.arange(m) / (m - 1))
for k in (3, 5, 8, 12, 20):
    r, _ = abs_minimax(xg, xg ** k, k - 1)
    exact = 2.0 ** (1 - k)
    print(f"   k={k:3d}   LP={r:.12e}   exact={exact:.12e}   rel.dev={abs(r/exact-1):.2e}")

print()
print("V2  exact interpolation: r*(d = #nodes-1) = 0 with boundedness relaxed")
rng = np.random.default_rng(7)
for nn in (6, 12, 20):
    xs = np.sort(rng.uniform(-0.9, 0.9, nn))
    ys = 1.0 + 0.5 * rng.uniform(size=nn)
    r_int, _ = core.min_rel_error(xs, ys, nn - 1, bound=1e6)
    r_low, _ = core.min_rel_error(xs, ys, nn - 2, bound=1e6)
    print(f"   nodes={nn:3d}  r*(d=n-1)={r_int:.3e}   r*(d=n-2)={r_low:.3e}")

print()
print("V3  grid independence of the bounded LP (BND_PER_DEG 3 / 6 / 12 / 24)")
lam = np.array([0.02, 0.08, 0.18, 0.32, 0.5, 0.72, 0.98, 1.28, 1.62, 2.0])
x, y, _ = core.make_problem(lam, "aware", "p60")
for bpd in (3, 6, 12, 24):
    core.BND_PER_DEG = bpd
    for d in (20, 40):
        r, cf = core.min_rel_error(x, y, d)
        print(f"   bpd={bpd:3d}  d={d:3d}  r*={r:.6e}  sup|p|(fine)={core.true_sup(cf):.6f}")
core.BND_PER_DEG = 6

print()
print("V4  Bernstein floor check.  For |p|<=1 on [-1,1], |p'(x)|<=d/sqrt(1-x^2).")
print("    Two bottom nodes x1<x2 near 0 with targets y1,y2 force d >= (y1-y2)/(x2-x1).")
for kappa in (20.0, 50.0, 120.0):
    x1, x2 = 1.0 / kappa, 4.0 / kappa          # quadratic clustering: lam_k ~ k^2 lam_1
    y1, y2 = 1.0, 0.5
    floor_aware = (y1 - y2) / (x2 - x1)
    x2g = 2.0 / kappa                           # continuum: nearest point at 2 lam_1
    floor_gen = (1.0 - 1 / np.sqrt(2)) / (x2g - x1)
    print(f"   kappa={kappa:6.1f}   d_aware >= {floor_aware:7.2f}   d_generic >= {floor_gen:7.2f}"
          f"   ratio={floor_gen/floor_aware:.3f}")


# ------------------------------------------------------------------ V5 / V6
def v5_v6():
    """Independent (non-LP) confirmation of the minimal aware degree, and a
    grid-refinement check of the generic arm."""
    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from geovac.sturmian_sigma_law import sw_cross_block, sigma_spectrum
    from scipy.optimize import minimize

    core.BND_PER_DEG = 6
    print("\nV5  aware/p60 minimal degree, LP vs an INDEPENDENT null-space minimisation")
    print("    (interpolate the 2n targets exactly -> affine family; minimise sup|p|")
    print("     over the family with a smooth-max + SLSQP, no linear programming)")
    for s, n in ((2.0, 4), (3.0, 4), (2.0, 6)):
        sig = sigma_spectrum(sw_cross_block(s, n))
        lam = np.sort(np.concatenate([1 - sig, 1 + sig]))
        x, y, _ = core.make_problem(lam, "aware", "p60", headroom=0.5)
        d_lp, _, _ = core.min_degree(x, y, 0.0, d_lo=2 * n, d_max=80)
        for d in (d_lp - 1, d_lp):
            A = core.cheb_design(x, d)
            # particular + null-space parametrisation of the interpolation constraints
            p0, *_ = np.linalg.lstsq(A, y, rcond=None)
            _, sv, Vt = np.linalg.svd(A)
            ns = Vt[np.sum(sv > 1e-10):].T
            m = 40 * (d + 1) + 1
            g = np.cos(np.pi * np.arange(m) / (m - 1))
            B = core.cheb_design(g, d)

            def sup(z, beta=400.0):
                v = np.abs(B @ (p0 + ns @ z))
                mx = v.max()
                return mx + np.log(np.mean(np.exp(beta * (v - mx)))) / beta

            best = np.inf
            for seed in range(4):
                z0 = np.zeros(ns.shape[1]) if seed == 0 else \
                    np.random.default_rng(seed).normal(0, 0.3, ns.shape[1])
                r = minimize(sup, z0, method="Nelder-Mead",
                             options=dict(maxiter=40000, maxfev=40000, xatol=1e-10,
                                          fatol=1e-12))
                best = min(best, float(np.max(np.abs(B @ (p0 + ns @ r.x)))))
            verdict = "feasible" if best <= 1.02 else "INFEASIBLE"
            print(f"    s={s} n={n}  d={d:3d}  min sup|p| over interpolants = {best:.4f}"
                  f"   -> {verdict}   (LP said d_min = {d_lp})")

    print("\nV6  generic/p60: accuracy-grid refinement (n_grid 60 / 130 / 400)")
    sig = sigma_spectrum(sw_cross_block(3.0, 6))
    lam = np.sort(np.concatenate([1 - sig, 1 + sig]))
    for ng in (60, 130, 400):
        x, y, _ = core.make_problem(lam, "generic", "p60", n_grid=ng, headroom=0.5)
        for d in (30, 60):
            r, _ = core.min_rel_error(x, y, d)
            print(f"    n_grid={ng:4d}  d={d:3d}  r*={r:.6e}")


if __name__ == "__main__" or True:
    pass
