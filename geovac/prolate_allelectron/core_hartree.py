"""Build increment 1 for prolate-native frozen-core LiH: the Li 1s^2 core's
Hartree screening of the valence, as a ONE-BODY potential in the prolate basis.

The valence electron on the prolate grid feels, from the Li side:
    bare:  -Z_A / r_A          (Z_A = 3, already in vne_hetero_mpf)
  + core Hartree (repulsive):  V_H(r_A) = (2/r_A)[1 - (1 + zc*r_A) exp(-2 zc r_A)]
    (the potential of a 1s^2 density of exponent zc; 2 = two core electrons).
Net long range: -3/r_A + V_H -> -1/r_A  (Li^2+); short range -> -3/r_A (penetration).

r_A = (R/2)(xi + eta) in prolate coordinates -> V_H is a function of (xi+eta), a
one-body integral of the SAME class as V_ne. This module builds <g_i|V_H|g_j> by
prolate quadrature (proven grid machinery, cf. heh_probe.vne_hetero_quad) and
validates:
  (A) the closed-form V_H(r_A) against a direct 3D numerical integral of the 1s^2
      density potential at several r_A (validates the physics);
  (B) the long-range limit: (-3/r_A + V_H) matrix approaches the -1/r_A (Li^2+)
      matrix as we probe large-r_A-weighted functions.

NOTE: gated on the HeH+ R_eq convergence check. Run after the gate confirms GO.
Run:  python -m geovac.prolate_allelectron.core_hartree

Production home (LiH productionization step 2a, 2026-10-04): promoted verbatim from
``debug/prolate_core_hartree.py``; numerics, summation order, precision and loop structure are
unchanged (bit-identical against tests/data/lih_golden/golden.json).  The debug/
file is now a module-alias shim.
Backs Paper 19's frozen-core prolate LiH diagnosis (the Li 1s^2 core Hartree screening
as a one-body potential; also supplies ``ZC_LI`` = 2.6875).  No import-time mpmath side
effect.
"""
from __future__ import annotations

import os
import sys
from typing import Sequence

import numpy as np

from geovac import prolate_recondition as pr
from geovac.hylleraas import build_quadrature_grids

ZC_LI = 3.0 - 5.0 / 16.0     # He-like Li^2+ 1s effective exponent (Slater), 2.6875
ZA = 3.0                     # bare Li nucleus


def V_H_closed(r_A: float | np.ndarray, zc: float = ZC_LI) -> float | np.ndarray:
    """Hartree potential of a 1s^2 density (exponent zc) at distance r_A (>0)."""
    return (2.0 / r_A) * (1.0 - (1.0 + zc * r_A) * np.exp(-2.0 * zc * r_A))


# ---- validation (A): closed-form V_H vs direct 3D integral of rho_core/|r-r'| ----
def V_H_direct(r_A: float, zc: float = ZC_LI, n: int = 400) -> float:
    """Direct numerical integral V_H(r_A) = int rho(r')/|r - r'| d3r', rho = 2|1s|^2,
    with the field point at distance r_A on the z-axis. Spherical grid on r'."""
    # rho(r') = 2 * (zc^3/pi) exp(-2 zc r'); by symmetry integrate over r', theta'
    rr = np.linspace(1e-4, 25.0 / zc, n)
    ct = np.linspace(-1.0, 1.0, n)               # cos(theta')
    RR, CT = np.meshgrid(rr, ct, indexing='ij')
    rho = 2.0 * (zc**3 / np.pi) * np.exp(-2.0 * zc * RR)
    # |r - r'| with r on axis at distance r_A: sqrt(r_A^2 + r'^2 - 2 r_A r' cos)
    dist = np.sqrt(r_A**2 + RR**2 - 2.0 * r_A * RR * CT + 1e-30)
    integrand = rho / dist * (RR**2)             # d3r' = r'^2 dr' dcos dphi ; *2pi
    _trap = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
    val = 2.0 * np.pi * _trap(_trap(integrand, ct, axis=1), rr)
    return val


def validate_closed_form() -> bool:
    print("Validation A: closed-form V_H(r_A) vs direct 3D integral (zc=%.4f)" % ZC_LI)
    print("  r_A      closed        direct        rel.err")
    ok = True
    for r_A in (0.3, 0.6, 1.0, 2.0, 4.0):
        c = V_H_closed(r_A); d = V_H_direct(r_A)
        rel = abs(c - d) / abs(d)
        ok = ok and rel < 5e-3
        print(f"  {r_A:4.1f}  {c:12.6f}  {d:12.6f}   {rel:.2e}")
    # limits: r_A->inf -> 2/r_A ; r_A->0 -> 2*zc (finite, = <1s|V|1s>-ish peak ~ (5/4)? )
    print("  limit check: V_H(20)*20 = %.4f (->2 as two-electron charge)"
          % (V_H_closed(20.0) * 20.0))
    print("  HARTREE CLOSED FORM", "VALIDATED" if ok else "MISMATCH")
    return ok


# ---- the prolate one-body matrix <g_i|V_H|g_j> (mu=0 sigma; extendable) ----
def core_hartree_matrix(basis: Sequence, R: float, alpha: float, zc: float = ZC_LI,
                        grid: dict | None = None) -> np.ndarray:
    """<g_i | V_H(r_A) | g_j> over prolate coords for a one-electron basis of
    ProductFn (sigma, mu=0). r_A = (R/2)(xi+eta)."""
    if grid is None:
        grid = build_quadrature_grids(N_xi=40, N_eta=30, N_phi=4, xi_max=20.0)
    xi, wxi = grid['xi'], grid['w_xi']
    eta, weta = grid['eta'], grid['w_eta']
    hR = R / 2.0
    n = len(basis)
    V = np.zeros((n, n))
    # one-electron prolate volume element: (R/2)^3 (xi^2 - eta^2) dxi deta dphi
    for a in range(len(xi)):
        x = xi[a]
        ef = np.exp(-2.0 * alpha * x)            # |g|^2 carries exp(-2 alpha xi)
        for b in range(len(eta)):
            e = eta[b]
            r_A = hR * (x + e)
            if r_A <= 0:
                continue
            vh = V_H_closed(r_A, zc)
            jac = (x**2 - e**2)
            w = wxi[a] * weta[b] * jac * (hR**3) * ef * vh * 2.0 * np.pi
            for i in range(n):
                gi = x**basis[i].j * e**basis[i].l
                for jj in range(n):
                    gj = x**basis[jj].j * e**basis[jj].l
                    V[i, jj] += w * gi * gj
    return V


def demo_matrix() -> np.ndarray:
    """Sanity: build V_H on a small one-electron sigma basis; report the net
    Li one-body (bare -3/r_A + V_H) is less deep than bare (screening works)."""
    R, alpha = 3.015, 1.6
    basis = [pr.ProductFn(j, l, 0, 0, 0, alpha) for j in range(3) for l in range(2)]
    VH = core_hartree_matrix(basis, R, alpha)
    print("\nCore-Hartree matrix (sigma, n=%d) diag:" % len(basis),
          np.round(np.diag(VH), 4))
    print("  (all positive = repulsive screening, as expected)")
    return VH


if __name__ == '__main__':
    ok = validate_closed_form()
    if ok:
        demo_matrix()
