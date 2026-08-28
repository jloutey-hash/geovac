"""Follow-up to ee_w_eigen_skeleton: identify the extra zero mode + test the 2n-1 law.

HYPOTHESIS (from the ns=3 surprise): the pair densities R_m R_n r^2 e^{-2r} are
polynomials of degree <= 2n-2 times a FIXED exponential, a space of dimension only
2n-1.  Consequences, all testable:
  H1  the function-space rank of the pair-density family is exactly 2n-1;
  H2  the matricization of the e-e tensor has exact rank <= 2n-1 for EVERY radial
      kernel (g, g_sep, W alike) -- linear in n, not quadratic;
  H3  the extra ns=3 null vector is a linear dependence of the DENSITIES themselves
      (so it is a null vector of every kernel matrix simultaneously), verifiable
      symbolically: sum_P c_P rho_P == 0 as a function.
Also: H4  the raw-Gram spectrum grows with ns (basis artifact); in the Loewdin-
      orthonormal orbital basis, does the W spectrum converge with ns?
      H5  resolution-identity check redone with the analytic tail (the earlier
      2.5e-3 discrepancy was exactly the 1/R_max cutoff).
"""
import io
import json
import os
import sys

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)
sys.path.insert(0, ROOT)

from geovac import transcorrelated_sturmian as TC  # noqa: E402

out = {}

print("=" * 78)
print("H1/H2 -- function-space rank of pair densities and kernel-independent 2n-1 law")
print("=" * 78)
print(f"{'ns':>3}{'2n-1':>6}{'dens rank':>11}{'rank g':>8}{'rank gsep':>10}"
      f"{'rank W':>8}{'rank(random kernel)':>21}")
rng = np.random.default_rng(1)
for ns in (2, 3, 4, 5, 6):
    r, wr = TC.make_grid(1.0, Ng=900)
    S, h1, Rtab, Wt = TC.build_one_body(ns, r, wr, 1.0, 2.0)
    prs = [(i, j) for i in range(1, ns + 1) for j in range(1, ns + 1)]
    Dm = np.array([Rtab[i][0] * Rtab[j][0] * Wt for (i, j) in prs])
    dens_rank = np.linalg.matrix_rank(Dm, tol=1e-10 * np.abs(Dm).max())
    R1g, R2g = np.meshgrid(r, r, indexing="ij")
    kernels = {
        "g": 1.0 / np.maximum(R1g, R2g),
        "gsep": 0.5 * (1.0 / R1g + 1.0 / R2g),
        "W": 0.5 * np.abs(1.0 / R1g - 1.0 / R2g),
    }
    # a RANDOM smooth radial kernel: rank must still be <= 2n-1 if H2 holds
    c = rng.normal(size=4)
    kernels["rand"] = (c[0] + c[1] * np.exp(-0.3 * (R1g + R2g))
                       + c[2] / (1 + R1g + R2g) + c[3] * np.exp(-0.1 * np.abs(R1g - R2g)))
    ranks = {}
    for nm, K in kernels.items():
        M = np.array([[float(Dm[a] @ K @ Dm[b]) for b in range(len(prs))]
                      for a in range(len(prs))])
        ranks[nm] = int(np.linalg.matrix_rank(M, tol=1e-10 * np.abs(M).max()))
    print(f"{ns:>3}{2*ns-1:>6}{dens_rank:>11}{ranks['g']:>8}{ranks['gsep']:>10}"
          f"{ranks['W']:>8}{ranks['rand']:>21}")
    out.setdefault("H1H2", {})[str(ns)] = dict(dens_rank=int(dens_rank), **ranks)

print()
print("=" * 78)
print("H3 -- the extra ns=3 null vector IS a density dependence (symbolic)")
print("=" * 78)
rv = sp.symbols("r", positive=True)
Rs = {n: sp.Rational(2, n) * sp.exp(-rv) * sp.assoc_laguerre(n - 1, 1, 2 * rv)
      for n in (1, 2, 3)}
prs3 = [(i, j) for i in range(1, 4) for j in range(1, 4)]
polys = [sp.expand(sp.cancel(Rs[i] * Rs[j] / sp.exp(-2 * rv)) * rv ** 2)
         for (i, j) in prs3]
# find the dependence among the 6 SYMMETRIC products (i<=j)
sym_idx = [a for a, (i, j) in enumerate(prs3) if i <= j]
basis = [rv ** p for p in range(2, 7)]                    # r^2..r^6: dim 5 (all densities carry r^2)
Mco = sp.Matrix([[sp.Poly(polys[a], rv).coeff_monomial(rv ** p) for p in range(2, 7)]
                 for a in sym_idx]).T
null = Mco.nullspace()
print(f"  symmetric products: {len(sym_idx)};  poly space dim: 5;  "
      f"dependencies found: {len(null)}")
if null:
    v = null[0]
    v = v / sp.gcd(tuple(v))
    combo = {str(prs3[sym_idx[a]]): str(sp.nsimplify(v[a])) for a in range(len(sym_idx))
             if v[a] != 0}
    residual = sp.expand(sum(v[a] * polys[sym_idx[a]] for a in range(len(sym_idx))))
    print(f"  exact dependence coefficients (per symmetric pair): {combo}")
    print(f"  symbolic residual of the combination: {residual}   (must be 0)")
    out["H3"] = dict(coeffs=combo, residual=str(residual))

print()
print("=" * 78)
print("H4 -- W spectrum in the Loewdin-ORTHONORMAL orbital basis vs ns")
print("=" * 78)
rows = {}
for ns in (2, 3, 4, 5, 6, 7, 8):
    r, wr = TC.make_grid(1.0, Ng=900)
    S, h1, Rtab, Wt = TC.build_one_body(ns, r, wr, 1.0, 2.0)
    R1g, R2g = np.meshgrid(r, r, indexing="ij")
    K_W = 0.5 * np.abs(1.0 / R1g - 1.0 / R2g)
    D = {(i, k): Rtab[i][0] * Rtab[k][0] * Wt for i in range(1, ns + 1)
         for k in range(1, ns + 1)}
    gw = np.zeros((ns,) * 4)
    for a, i in enumerate(range(1, ns + 1)):
        for b, j in enumerate(range(1, ns + 1)):
            for c2, kk in enumerate(range(1, ns + 1)):
                for d, ll in enumerate(range(1, ns + 1)):
                    gw[a, b, c2, d] = D[(i, kk)] @ K_W @ D[(j, ll)]
    X = TC.lowdin(S)
    gwo = TC.transform_2(gw, X)
    Wm = gwo.transpose(0, 2, 1, 3).reshape(ns * ns, ns * ns)
    lam = np.sort(np.abs(np.linalg.eigvalsh(Wm)))[::-1]
    rows[str(ns)] = [float(x) for x in lam[:4]]
    print(f"  ns={ns}: " + "  ".join(f"{x:.6f}" for x in lam[:4]))
out["H4"] = rows

print()
print("=" * 78)
print("H5 -- resolution identity with the analytic tail")
print("=" * 78)
rs = (0.3, 1.7)
Rmax = 400.0
Rg = np.linspace(1e-6, Rmax, 4_000_001)
res = np.trapezoid(((Rg > rs[0]) & (Rg > rs[1])) / Rg ** 2, Rg) + 1.0 / Rmax
print(f"  with tail: {res:.8f}   vs 1/max = {1/max(rs):.8f}   "
      f"dev {abs(res - 1/max(rs)):.1e}")
out["H5"] = dict(value=float(res), exact=1 / max(rs))

os.makedirs("debug/data", exist_ok=True)
with io.open("debug/data/ee_w_eigen_followup.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=2)
print("\nwrote debug/data/ee_w_eigen_followup.json")
