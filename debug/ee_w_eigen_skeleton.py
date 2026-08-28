"""PROBE: do W's leading eigenvectors have closed SKELETON forms?

Follow-on to rem:ee_partial_split (PI-directed).  W is the irreducible e-e remainder
(kernel |1/r1 - 1/r2|/2), scale-free after g(k) = k g(1).

PRE-REGISTERED
  P1  W(1), g(1), S are EXACTLY RATIONAL at k = 1 (the code's s-Sturmians are
      rational-coefficient polynomials x e^{-r}; the kernel integrates piecewise
      to rationals).  Computed exactly in sympy; numerics must agree.
  P2  the n^2 matricization has n(n-1)/2 EXACT zero modes (antisymmetric pair
      combinations: rho_ik = rho_ki), so the active block has dim n(n+1)/2.
  P3  THE DECIDABLE QUESTION: does the characteristic polynomial of the active
      block factor over Q into low-degree pieces?
        factors  -> closed forms exist; identify the leading eigenpair exactly.
        irreducible -> the finite-matrix eigenvectors are generic algebraic
                       numbers of degree n(n+1)/2; no simple closed form AT THE
                       MATRIX LEVEL -- the skeleton object is then P4.
  P4  in reciprocal radius u = 1/r the kernels are CLASSICAL:
        1/r_>       = min(u1, u2)   (Brownian covariance; Green's fn of -d^2/du^2,
                                     Dirichlet at u = 0)
        W's kernel  = |u1 - u2|/2   (d^2/du^2 of it = delta)
      so continuum eigenfunctions of the W kernel w.r.t. ANY weight w satisfy the
      ODE  lambda f''(u) = w(u) f(u).  Verified on a grid.  Also the resolution
      identity 1/max(r1,r2) = INT_0^inf 1[r1<R] 1[r2<R] dR/R^2 (the PSD/Gram
      reading behind Coulomb positivity).
  P5  leading eigenvalues of W(1) vs ns: measured (do they converge?).
"""
import io
import json
import os
import sys
import time

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)
sys.path.insert(0, ROOT)

from geovac import transcorrelated_sturmian as TC  # noqa: E402

out = {}
r1, r2 = sp.symbols("r1 r2", positive=True)


def R_sym(n):
    """code-convention s-Sturmian at k = 1: (2/n) e^{-r} L^1_{n-1}(2r)."""
    return sp.Rational(2, n) * sp.exp(-r1) * sp.assoc_laguerre(n - 1, 1, 2 * r1)


def pair_density(m, n, var):
    p = (sp.Rational(2, m) * sp.assoc_laguerre(m - 1, 1, 2 * var)
         * sp.Rational(2, n) * sp.assoc_laguerre(n - 1, 1, 2 * var))
    return sp.expand(p * var ** 2), sp.exp(-2 * var)


def exact_entry(kernel_case, m1, n1, m2, n2):
    """exact double integral of rho_{m1n1}(r1) K(r1,r2) rho_{m2n2}(r2)."""
    p1, e1 = pair_density(m1, n1, r1)
    p2raw, _ = pair_density(m2, n2, r1)
    p2 = p2raw.subs(r1, r2)
    e2 = sp.exp(-2 * r2)
    if kernel_case == "W":
        # region r1 < r2: (1/r1 - 1/r2)/2 ; mirror by symmetry of |.|
        inner = sp.integrate(p1 * e1 * (1 / r1 - 1 / r2) / 2, (r1, 0, r2))
        reg1 = sp.integrate(sp.expand(p2 * e2 * inner), (r2, 0, sp.oo))
        # region r1 > r2 == swap the two densities in region r1 < r2
        q1, _ = pair_density(m2, n2, r1)
        q2raw, _ = pair_density(m1, n1, r1)
        q2 = q2raw.subs(r1, r2)
        inner2 = sp.integrate(q1 * e1 * (1 / r1 - 1 / r2) / 2, (r1, 0, r2))
        reg2 = sp.integrate(sp.expand(q2 * e2 * inner2), (r2, 0, sp.oo))
        return sp.nsimplify(sp.simplify(reg1 + reg2), rational=True)
    if kernel_case == "g":
        inner = sp.integrate(p1 * e1 / r2, (r1, 0, r2)) \
              + sp.integrate(p1 * e1 / r1, (r1, r2, sp.oo))
        return sp.nsimplify(sp.simplify(
            sp.integrate(sp.expand(p2 * e2 * inner), (r2, 0, sp.oo))), rational=True)
    raise ValueError(kernel_case)


def exact_matrices(ns):
    pairs = [(i, j) for i in range(1, ns + 1) for j in range(1, ns + 1)]
    npair = len(pairs)
    Wx = sp.zeros(npair, npair)
    cache = {}
    for a, (m1, n1) in enumerate(pairs):
        for b, (m2, n2) in enumerate(pairs):
            if b < a:
                continue
            key = tuple(sorted([tuple(sorted((m1, n1))), tuple(sorted((m2, n2)))]))
            if key not in cache:
                cache[key] = exact_entry("W", m1, n1, m2, n2)
            Wx[a, b] = Wx[b, a] = cache[key]
    return pairs, Wx


print("=" * 78)
print("P1/P2/P3 -- exact rational W(1); zero modes; characteristic polynomial over Q")
print("=" * 78)
lam = sp.symbols("lam")
for ns in (2, 3):
    t0 = time.time()
    pairs, Wx = exact_matrices(ns)
    # P1: rationality
    all_rat = all(e.is_Rational for e in Wx)
    # numeric cross-check vs the grid engine
    r, wr = TC.make_grid(1.0, Ng=900)
    S, h1, Rtab, Wt = TC.build_one_body(ns, r, wr, 1.0, 2.0)
    R1g, R2g = np.meshgrid(r, r, indexing="ij")
    K_W = 0.5 * np.abs(1.0 / R1g - 1.0 / R2g)
    D = {(i, k): Rtab[i][0] * Rtab[k][0] * Wt for i in range(1, ns + 1)
         for k in range(1, ns + 1)}
    Wnum = np.array([[float(D[p] @ K_W @ D[q]) for q in pairs] for p in pairs])
    Wex = np.array([[float(Wx[a, b]) for b in range(len(pairs))]
                    for a in range(len(pairs))])
    num_dev = float(np.abs(Wnum - Wex).max())
    # P2: exact zero modes
    null = Wx.nullspace()
    # P3: char poly, factored
    cp = sp.factor(Wx.charpoly(lam).as_expr())
    print(f"\nns = {ns}  ({time.time()-t0:.0f}s exact)")
    print(f"  all entries rational: {all_rat};   grid-vs-exact max dev {num_dev:.1e}")
    print(f"  sample entries: W[1s^2,1s^2] = {Wx[0,0]},  W[1s^2,{pairs[1]}] = {Wx[0,1]}")
    print(f"  exact nullspace dim = {len(null)}   (predicted {ns*(ns-1)//2})")
    print(f"  char poly factored over Q:")
    print(f"    {cp}")
    out[f"ns={ns}"] = dict(all_rational=bool(all_rat), num_dev=num_dev,
                           W00=str(Wx[0, 0]), null_dim=len(null),
                           charpoly_factored=str(cp))

print()
print("=" * 78)
print("P4 -- reciprocal radius: Coulomb = min(u,u'); W-kernel = |u-u'|/2;")
print("      eigenfunctions solve  lambda f'' = w f   (verified on a grid)")
print("=" * 78)
# resolution identity, numeric spot check
rs = np.array([0.3, 1.7]); Rg = np.linspace(1e-6, 400.0, 4_000_001)
res = np.trapezoid(((Rg[None, :] > rs[0]) & (Rg[None, :] > rs[1])) / Rg ** 2,
                   Rg, axis=1)
print(f"  resolution identity: INT 1[r1<R]1[r2<R] dR/R^2 = {res[0]:.8f}"
      f"   vs 1/max = {1/max(rs):.8f}")
# min(u,u') check is the same identity in u; now the ODE:
U, Nu = 12.0, 3000
u = np.linspace(0, U, Nu)
du = u[1] - u[0]
w = np.exp(-u)                                   # any smooth weight
Kmat = 0.5 * np.abs(u[:, None] - u[None, :])
T = Kmat * (w[None, :] * du)                     # (Tf)(u) = INT K w f du'
ev, V = np.linalg.eig(T)
idx = np.argsort(-np.abs(ev))
ode_res = []
for j in idx[:3]:
    lam_j, f = float(np.real(ev[j])), np.real(V[:, j])
    fpp = np.gradient(np.gradient(f, du), du)
    lhs, rhs = lam_j * fpp[5:-5], (w * f)[5:-5]
    ode_res.append(float(np.linalg.norm(lhs - rhs) / np.linalg.norm(rhs)))
print("  ODE residual |lam f'' - w f| / |w f| for top-3 eigenfunctions: "
      + ", ".join(f"{x:.1e}" for x in ode_res))
out["P4"] = dict(resolution_dev=float(abs(res[0] - 1 / max(rs))), ode_residuals=ode_res)

print()
print("=" * 78)
print("P5 -- leading eigenvalues of W(1) vs ns (numeric)")
print("=" * 78)
rows = {}
for ns in (2, 3, 4, 5, 6, 7):
    r, wr = TC.make_grid(1.0, Ng=900)
    S, h1, Rtab, Wt = TC.build_one_body(ns, r, wr, 1.0, 2.0)
    R1g, R2g = np.meshgrid(r, r, indexing="ij")
    K_W = 0.5 * np.abs(1.0 / R1g - 1.0 / R2g)
    D = {(i, k): Rtab[i][0] * Rtab[k][0] * Wt for i in range(1, ns + 1)
         for k in range(1, ns + 1)}
    prs = [(i, j) for i in range(1, ns + 1) for j in range(1, ns + 1)]
    Wm = np.array([[float(D[p] @ K_W @ D[q]) for q in prs] for p in prs])
    lamv = np.sort(np.abs(np.linalg.eigvalsh(Wm)))[::-1]
    rows[ns] = [float(x) for x in lamv[:4]]
    print(f"  ns={ns}: " + "  ".join(f"{x:.6f}" for x in lamv[:4]))
out["P5"] = rows

os.makedirs("debug/data", exist_ok=True)
with io.open("debug/data/ee_w_eigen_skeleton.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=2)
print("\nwrote debug/data/ee_w_eigen_skeleton.json")
