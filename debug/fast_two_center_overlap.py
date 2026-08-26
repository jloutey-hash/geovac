"""Fast float64 two-center overlap for GeoVac hydrogenic {1s, 2p} orbitals, matching the
topos3 (compute_topos3_two_center_meet.overlap_two_center) prolate-spheroidal convention
EXACTLY.  {1s,2p} are single-exponential STOs (zero radial nodes), so a numpy Gauss-Legendre
quadrature gives ~1e-9 in milliseconds -- the mpmath engine is correct but too slow for the
diffuse Z=1 2p m=1 (pi) integral needed by the bending study.

Convention (identical to topos3):
  phi_{nlm}(r,Z) = R_nl(Z,n,l,r) * theta_norm(l,|m|) * P_l^{|m|}(cos theta)   [azimuthal analytic]
  overlap = \int dxi \int deta  phi_A(r1,ct1) phi_B(r2,ct2) (R/2)^3 (xi^2-eta^2)
  r1=(R/2)(xi+eta), r2=(R/2)(xi-eta),
  ct1=(1+xi eta)/(xi+eta), ct2=(xi eta-1)/(xi-eta),  xi in [1,inf), eta in [-1,1].
"""
from __future__ import annotations
from math import factorial, comb, sqrt
import numpy as np


def _genlag(k, alpha, x):
    # generalized Laguerre L_k^alpha(x); k=0 -> 1 (all our states)
    if k == 0:
        return np.ones_like(x)
    return sum(((-1) ** j) * comb(k + alpha, k - j) / factorial(j) * x ** j
               for j in range(k + 1))


def _R_nl(Z, n, l, r):
    Z = float(Z)
    rho = 2 * Z * r / n
    norm = sqrt((2 * Z / n) ** 3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    return norm * np.exp(-rho / 2) * rho ** l * _genlag(n - l - 1, 2 * l + 1, rho)


def _theta_norm(l, m):
    return sqrt((2 * l + 1) / 2 * factorial(l - m) / factorial(l + m))


def _Plm(l, m, x):
    # Condon-Shortley associated Legendre for l<=1 (all we need)
    if l == 0:
        return np.ones_like(x)
    if l == 1:
        if m == 0:
            return x
        if m == 1:
            return -np.sqrt(np.clip(1 - x * x, 0, None))
    raise ValueError("only l<=1 supported")


def _angular(l, m, ct):
    return _theta_norm(l, abs(m)) * _Plm(l, abs(m), ct)


# cache GL nodes/weights
_GL = {}


def _gl(n):
    if n not in _GL:
        _GL[n] = np.polynomial.legendre.leggauss(n)
    return _GL[n]


def overlap_fast(Z1, n1, l1, Z2, n2, l2, m, R, n_xi=400, n_eta=80):
    """<phi_{n1 l1 m}(0,Z1) | phi_{n2 l2 m}(R zhat, Z2)>  (fast float64)."""
    R = float(R)
    half = R / 2.0
    zeta1, zeta2 = Z1 / n1, Z2 / n2                 # xi-decay rate = half*(zeta1+zeta2)
    rate = half * (zeta1 + zeta2)
    xi_max = 1.0 + 40.0 / rate                      # e^{-40} tail
    # eta in [-1,1]
    xe, we = _gl(n_eta)
    # xi in [1, xi_max]
    xx, wx = _gl(n_xi)
    xi = 0.5 * (xi_max - 1.0) * xx + 0.5 * (xi_max + 1.0)
    wxi = wx * 0.5 * (xi_max - 1.0)
    XI, ETA = np.meshgrid(xi, xe, indexing="ij")     # (n_xi, n_eta)
    Wgrid = np.outer(wxi, we)
    r1 = half * (XI + ETA)
    r2 = half * (XI - ETA)
    ct1 = (1 + XI * ETA) / (XI + ETA)
    ct2 = (XI * ETA - 1) / (XI - ETA)
    integ = (_R_nl(Z1, n1, l1, r1) * _angular(l1, m, ct1)
             * _R_nl(Z2, n2, l2, r2) * _angular(l2, m, ct2)
             * half ** 3 * (XI ** 2 - ETA ** 2))
    return float(np.sum(integ * Wgrid))


def funds_fast(R, Za, Zb):
    """(ss, SP=<s_A|pz_B>, PS=<pz_A|s_B>, ppsigma, pppi) for charges (Za,Zb) at separation R."""
    return (overlap_fast(Za, 1, 0, Zb, 1, 0, 0, R),
            overlap_fast(Za, 1, 0, Zb, 2, 1, 0, R),
            overlap_fast(Za, 2, 1, Zb, 1, 0, 0, R),
            overlap_fast(Za, 2, 1, Zb, 2, 1, 0, R),
            overlap_fast(Za, 2, 1, Zb, 2, 1, 1, R))


if __name__ == "__main__":
    import time
    # validate vs bargmann cache (m=0 fundamentals) + vs a known analytic 1s-1s overlap
    z = np.load("debug/data/beh2_overlap_cache.npz")
    RG, SBeH, SHH = z["RG"], z["SBeH"], z["SHH"]
    print("=== validate fast overlaps vs bargmann mpmath cache (m=0 block {1s,2p0}) ===")
    maxerr = 0.0
    for k in range(0, len(RG), 3):
        R = float(RG[k])
        for pair, Sc in [((2, 1), SBeH[k]), ((1, 1), SHH[k])]:
            ss = overlap_fast(pair[0], 1, 0, pair[1], 1, 0, 0, R)
            sp = overlap_fast(pair[0], 1, 0, pair[1], 2, 1, 0, R)
            ps = overlap_fast(pair[0], 2, 1, pair[1], 1, 0, 0, R)
            pp = overlap_fast(pair[0], 2, 1, pair[1], 2, 1, 0, R)
            got = np.array([[ss, sp], [ps, pp]])
            e = np.max(np.abs(got - Sc))
            maxerr = max(maxerr, e)
            if k % 6 == 0:
                print(f"  R={R:5.2f} pair{pair}: max|fast-cache|={e:.1e}")
    print(f"  MAX ERROR vs cache = {maxerr:.2e}  ({'OK' if maxerr < 1e-6 else 'FAIL'})")
    # analytic 1s-1s (Z=1) = (1+R+R^2/3)e^{-R}
    for R in [1.0, 2.0, 3.5]:
        an = (1 + R + R * R / 3) * np.exp(-R)
        fa = overlap_fast(1, 1, 0, 1, 1, 0, 0, R)
        print(f"  <1s|1s>(R={R}): analytic={an:.9f} fast={fa:.9f} err={abs(an-fa):.1e}")
    t = time.time(); v = overlap_fast(1, 2, 1, 1, 2, 1, 1, 2.4)
    print(f"  HH pppi R=2.4 = {v:.9f}  ({(time.time()-t)*1000:.1f} ms)  <- the slow-in-mpmath one")
