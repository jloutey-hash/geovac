"""Bet 1, step 0: numerical ground truth for the 3-centre ONE-electron integral.

    I = <chi_Y | -Z_X / |r - X| | chi_Z>
      = integral  conj(chi_Y(r-Y)) chi_Z(r-Z) (-Z_X/|r-X|) d^3r

with Y, Z, X three DISTINCT centres (X off the Y-Z axis, so axial symmetry is
genuinely broken). This is the object Bet 1 asks a closed form for; any closed
form must reproduce this number. Nothing spheroidal here -- deliberately an
independent reference.

Singularity trick: integrate in spherical coordinates CENTRED AT X. Then
|r-X| = s and d^3r = s^2 ds dOmega, so the 1/s is cancelled by s^2, leaving
s ds dOmega -- a smooth integrand. Gauss-Legendre in s and cos(theta), spectral
(periodic-trapezoid) in phi. Convergence is checked by grid refinement.

Run:  python debug/bet1_three_center_1e_reference.py
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import sph_harm

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from geovac.two_center_eri import radial_norm, radial_poly  # noqa: E402


def chi_grid(Z, n, l, m, V):
    """chi at global points V (..,3) relative to its own centre already removed.

    V is the displacement r - centre, shape (...,3). Returns complex array.
    Matches the codebase chi convention (validated to 1e-18 against it).
    """
    x, y, z = V[..., 0], V[..., 1], V[..., 2]
    r = np.sqrt(x * x + y * y + z * z)
    rsafe = np.where(r > 0, r, 1.0)
    th = np.arccos(np.clip(z / rsafe, -1.0, 1.0))
    ph = np.arctan2(y, x)
    c, a = radial_poly(Z, n, l)
    N = radial_norm(Z, n, l)
    rad = np.zeros_like(r)
    for k, cc in c.items():
        rad = rad + float(N * cc) * r ** k
    rad = rad * np.exp(-float(a) * r)
    return rad * sph_harm(m, l, ph, th)


def three_center_1e(ZY, orbY, Yc, ZZ, orbZ, Zc, ZX, Xc,
                    n_s=400, n_u=64, n_phi=64, s_max=35.0):
    Yc, Zc, Xc = map(np.asarray, (Yc, Zc, Xc))
    # s in [0, s_max] via Gauss-Legendre
    xs, ws = leggauss(n_s)
    s = 0.5 * s_max * (xs + 1.0)
    ws_s = 0.5 * s_max * ws
    # cos(theta) in [-1,1] via Gauss-Legendre
    u, wu = leggauss(n_u)
    th = np.arccos(u)
    sin_th = np.sqrt(1.0 - u * u)
    # phi uniform, spectral for periodic integrand
    phi = 2.0 * np.pi * np.arange(n_phi) / n_phi
    w_phi = 2.0 * np.pi / n_phi

    # direction grid (n_u, n_phi, 3)
    nx = np.outer(sin_th, np.cos(phi))
    ny = np.outer(sin_th, np.sin(phi))
    nz = np.outer(np.ones_like(sin_th), np.ones_like(phi)) * u[:, None]
    ndir = np.stack([nx, ny, nz], axis=-1)  # (n_u, n_phi, 3)

    total = 0.0 + 0.0j
    for si, wsi in zip(s, ws_s):
        P = Xc[None, None, :] + si * ndir          # global points (n_u,n_phi,3)
        vY = P - Yc[None, None, :]
        vZ = P - Zc[None, None, :]
        rho = np.conj(chi_grid(ZY, *orbY, vY)) * chi_grid(ZZ, *orbZ, vZ)
        # angular quadrature: sum_u wu * sum_phi w_phi * rho
        ang = (wu[:, None] * rho).sum() * w_phi
        total += wsi * si * ang
    return complex(-ZX * total)


def converged(label, args, kwargs_lo, kwargs_hi):
    lo = three_center_1e(*args, **kwargs_lo)
    hi = three_center_1e(*args, **kwargs_hi)
    print(f"  {label}")
    print(f"    coarse : {lo.real:+.12f}  ({lo.imag:+.3e} i)")
    print(f"    fine   : {hi.real:+.12f}  ({hi.imag:+.3e} i)")
    print(f"    |delta|: {abs(hi-lo):.2e}   -> {'CONVERGED' if abs(hi-lo)<1e-8 else 'refine more'}")
    return hi


def main():
    Z1 = Fraction(1)
    Yc = (0.0, 0.0, 0.0)
    Zc = (0.0, 0.0, 2.0)
    Xc = (1.3, 0.0, 0.7)          # off the Y-Z (z) axis: axial symmetry broken
    ZX = 1.0
    print("3-centre one-electron reference  <chi_Y|-Z_X/|r-X||chi_Z>")
    print(f"  Y=1s@{Yc}  Z-centre@{Zc}  X(source)@{Xc}  Z_X={ZX}\n")

    lo = dict(n_s=300, n_u=48, n_phi=48, s_max=30.0)
    hi = dict(n_s=500, n_u=80, n_phi=80, s_max=40.0)

    vals = {}
    vals["A 1s_Y x 1s_Z"] = converged(
        "A  1s_Y x 1s_Z  (sigma_density=0)",
        (Z1, (1, 0, 0), Yc, Z1, (1, 0, 0), Zc, ZX, Xc), lo, hi)
    vals["B 2p0_Y x 1s_Z"] = converged(
        "B  2p0_Y x 1s_Z  (l>0 one side)",
        (Z1, (2, 1, 0), Yc, Z1, (1, 0, 0), Zc, ZX, Xc), lo, hi)
    vals["C 2p+1_Y x 2p0_Z"] = converged(
        "C  2p+1_Y x 2p0_Z  (sigma_density=-1, complex)",
        (Z1, (2, 1, 1), Yc, Z1, (2, 1, 0), Zc, ZX, Xc), lo, hi)

    print("\n  reference values (fine grid):")
    for k, v in vals.items():
        print(f"    {k:20s} = {v.real:+.12f}  {v.imag:+.3e} i")


if __name__ == "__main__":
    main()
