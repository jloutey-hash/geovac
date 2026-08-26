"""Independent numerical oracle for two-center ONE-electron integrals.

Purpose: validate a *separate* (production) two-center grid engine against a
numerical reference that shares no code with GeoVac's closed-form integral
machinery. This file builds that reference from scratch: direct 3D
quadrature of the defining integrals, done TWO independent ways, plus
hand-derived closed forms for the 1s-1s (l=0) case.

Orbital convention (nodeless hydrogenic = Slater-type orbital, "STO"):

    chi(r) = N * r^l * Y_lm(theta_c, phi_c) * exp(-zeta * r_c)

  - Y_lm are the COMPLEX, Condon-Shortley-phase, unit-normalized spherical
    harmonics (int |Y_lm|^2 dOmega = 1), quantized along a single shared lab
    z-axis (both centers sit on that axis: A at the origin, B at (0,0,R)).
  - "Nodeless" means principal quantum number n = l + 1 (no radial nodes):
    the radial factor is the bare power r^l, not a Laguerre polynomial.
  - N is fixed by int |chi|^2 d^3r = 1:

        N^2 * (2l+2)! / (2 zeta)^{2l+3} = 1   =>   N = sqrt((2 zeta)^{2l+3} / (2l+2)!)

    (l=0: N = 2 zeta^{3/2}; l=1: N = sqrt(4/3) zeta^{5/2}.)

Three reference integrals, orbitals on ANY centers, l in {0, 1}:

    s_ref(oi, oj, R) = <chi_i | chi_j>
    t_ref(oi, oj, R) = <chi_i | -1/2 grad^2 | chi_j>
    v_ref(oi, oj, R) = <chi_i | (-1/r_A - 1/r_B) | chi_j>     (Z_A = Z_B = 1)

Two independent quadrature engines
-----------------------------------
Engine 1 (primary) -- prolate spheroidal, foci at A and B:

    r_A = (R/2)(xi + eta), r_B = (R/2)(xi - eta), z = (R/2)(1 + xi eta)
    d^3r = (R/2)^3 (xi^2 - eta^2) dxi deta dphi,  xi >= 1, eta in [-1,1]

  The xi integral is done with Gauss-Laguerre quadrature EXACTLY matched to
  the pair's combined decay rate p = (zeta_i + zeta_j) R / 2 (a change of
  variables u = p(xi-1) turns int_1^inf f(xi) e^{-p xi} dxi into a standard
  Gauss-Laguerre sum); since every one-electron integrand here reduces,
  after multiplying by the (xi+eta)(xi-eta) Jacobian, to a FINITE-DEGREE
  POLYNOMIAL in xi times e^{-p xi} (the 1/r_A, 1/r_B, and kinetic 1/r
  singularities are all simple poles that the Jacobian's (xi+eta) or
  (xi-eta) factor cancels exactly), this quadrature is essentially EXACT
  (machine precision) at modest node counts. The eta integral (polynomial
  times e^{-q eta} on the finite interval [-1,1]) uses plain high-order
  Gauss-Legendre, which converges geometrically fast for this entire
  integrand. phi is done by the uniform trapezoid rule, which is spectrally
  exact for the periodic e^{i(m_j - m_i) phi} integrand.

Engine 2 (cross-check) -- spherical grid about the bond MIDPOINT:

    r in [0, inf) via an algebraic map to a Gauss-Legendre grid,
    cos(theta) via Gauss-Legendre, phi via uniform trapezoid.

  This grid shares no coordinate structure with Engine 1 (different origin,
  different radial variable, no exact polynomial cancellation to exploit),
  so agreement between the two engines is a genuine independent check.

Kinetic energy: -1/2 grad^2 acting on chi = N r^l Y_lm(theta_c,phi_c) e^{-zeta r_c}
is evaluated via the exact identity (Y_lm is an eigenfunction of the angular
part of the Laplacian, eigenvalue -l(l+1)):

    grad^2[f(r) Y_lm] = Y_lm * [f''(r) + (2/r) f'(r) - l(l+1)/r^2 f(r)]

For f(r) = N r^l e^{-zeta r}, direct differentiation gives (the r^{l-2} terms
cancel identically -- a standard STO Laplacian identity):

    f'' + (2/r) f' - l(l+1)/r^2 f = N e^{-zeta r} [zeta^2 r^l - 2 zeta(l+1) r^{l-1}]

so

    grad^2 chi = chi * zeta^2  -  2 zeta (l+1) * chi / r

This is evaluated ANALYTICALLY at each quadrature node -- no numerical
differentiation anywhere in the two production engines. A finite-difference
Laplacian check (5-point central stencil, O(h^4)) is included in __main__ as
an independent confirmation of the identity itself.

geovac.qfd_core is imported ONLY inside the __main__ validation block, to
cross-check the l=0 (s-only) numbers against GeoVac's production closed
forms. Nothing in the reference engines above imports from geovac.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache
from typing import Tuple, Union

import numpy as np

Center = Union[str, Tuple[float, float, float]]


# --------------------------------------------------------------- orbital type


@dataclass(frozen=True)
class Orbital:
    """(center, zeta, l, m). center is 'A', 'B', or an explicit (x,y,z)."""

    center: Center
    zeta: float
    l: int
    m: int

    def __post_init__(self) -> None:
        if self.l not in (0, 1):
            raise ValueError(f"only l in {{0,1}} supported, got l={self.l}")
        if abs(self.m) > self.l:
            raise ValueError(f"|m|={abs(self.m)} > l={self.l}")
        if self.zeta <= 0:
            raise ValueError("zeta must be positive")


def center_xyz(center: Center, R: float) -> np.ndarray:
    """A -> origin, B -> (0,0,R); anything else is taken as an explicit vector."""
    if isinstance(center, str):
        c = center.upper()
        if c == "A":
            return np.array([0.0, 0.0, 0.0])
        if c == "B":
            return np.array([0.0, 0.0, R])
        raise ValueError(f"unknown center label {center!r}")
    return np.asarray(center, dtype=float)


def _center_sign(center: Center) -> int:
    """+1 for A, -1 for B -- the sign of eta's coefficient in r_center(xi,eta)."""
    if isinstance(center, str):
        c = center.upper()
        if c == "A":
            return 1
        if c == "B":
            return -1
    raise ValueError("prolate-spheroidal engine requires center 'A' or 'B'")


# ---------------------------------------------------------------------- STO


def sto_norm(zeta: float, l: int) -> float:
    """N such that int |N r^l e^{-zeta r}|^2 r^2 dr = 1 (angular part |Y_lm|^2
    is separately normalized to 1 over solid angle): N^2 = (2 zeta)^{2l+3} / (2l+2)!."""
    return math.sqrt((2.0 * zeta) ** (2 * l + 3) / math.factorial(2 * l + 2))


def Ylm(l: int, m: int, theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """Complex, Condon-Shortley-phase, unit-normalized spherical harmonics,
    hard-coded closed forms for l in {0, 1} (no dependence on any external
    special-function library's harmonic-ordering convention)."""
    theta = np.asarray(theta, dtype=float)
    phi = np.asarray(phi, dtype=float)
    if l == 0:
        return np.full(np.broadcast(theta, phi).shape,
                        1.0 / math.sqrt(4.0 * math.pi), dtype=complex)
    if l == 1:
        if m == 0:
            return (math.sqrt(3.0 / (4.0 * math.pi)) * np.cos(theta)).astype(complex)
        if m == 1:
            return (-math.sqrt(3.0 / (8.0 * math.pi)) * np.sin(theta)
                    * np.exp(1j * phi))
        if m == -1:
            return (math.sqrt(3.0 / (8.0 * math.pi)) * np.sin(theta)
                    * np.exp(-1j * phi))
    raise ValueError(f"l={l}, m={m} not supported")


def chi_value(orb: Orbital, r_c: np.ndarray, theta_c: np.ndarray,
              phi_c: np.ndarray) -> np.ndarray:
    """chi_orb evaluated at points given by their (r,theta,phi) relative to
    orb's OWN center."""
    N = sto_norm(orb.zeta, orb.l)
    return (N * np.asarray(r_c, dtype=float) ** orb.l
            * Ylm(orb.l, orb.m, theta_c, phi_c)
            * np.exp(-orb.zeta * np.asarray(r_c, dtype=float)))


def lap_chi_value(orb: Orbital, r_c: np.ndarray, theta_c: np.ndarray,
                   phi_c: np.ndarray) -> np.ndarray:
    """grad^2 chi_orb, exact radial-Laplacian identity (see module docstring):
    grad^2 chi = chi * [zeta^2 - 2 zeta (l+1) / r]."""
    r_c = np.asarray(r_c, dtype=float)
    zeta, l = orb.zeta, orb.l
    r_safe = np.where(r_c < 1e-300, 1e-300, r_c)
    coeff = zeta ** 2 - 2.0 * zeta * (l + 1) / r_safe
    return coeff * chi_value(orb, r_c, theta_c, phi_c)


# ------------------------------------------------------------- Gauss caches


@lru_cache(maxsize=None)
def _gauss_legendre(n: int):
    x, w = np.polynomial.legendre.leggauss(n)
    return x, w


@lru_cache(maxsize=None)
def _gauss_laguerre(n: int):
    x, w = np.polynomial.laguerre.laggauss(n)
    return x, w


def _phi_grid(n_phi: int):
    phi = 2.0 * np.pi * np.arange(n_phi) / n_phi
    w = np.full(n_phi, 2.0 * np.pi / n_phi)
    return phi, w


# --------------------------------------------- Engine 1: prolate spheroidal


def _prolate_grid(oi: Orbital, oj: Orbital, R: float,
                   n_xi: int = 40, n_eta: int = 100, n_phi: int = 16):
    """Grid + weights for the two-center integral of chi_i and (chi_j or its
    Laplacian), matched to THIS orbital pair's combined decay rate p."""
    p = 0.5 * R * (oi.zeta + oj.zeta)
    q = 0.5 * R * (oi.zeta * _center_sign(oi.center)
                   + oj.zeta * _center_sign(oj.center))
    if p <= 0:
        raise ValueError("combined decay rate p must be positive")

    u, wu = _gauss_laguerre(n_xi)
    xi = 1.0 + u / p
    w_xi = wu * np.exp(u) / p          # exact re-weighting, see module docstring

    eta, w_eta = _gauss_legendre(n_eta)
    phi, w_phi = _phi_grid(n_phi)

    XI = xi[:, None, None]
    ETA = eta[None, :, None]
    PHI = phi[None, None, :]
    WXI = w_xi[:, None, None]
    WETA = w_eta[None, :, None]
    WPHI = w_phi[None, None, :]

    rA = (R / 2.0) * (XI + ETA)
    rB = (R / 2.0) * (XI - ETA)
    z = (R / 2.0) * (1.0 + XI * ETA)
    rho2 = np.clip((XI ** 2 - 1.0) * (1.0 - ETA ** 2), 0.0, None)
    rho = (R / 2.0) * np.sqrt(rho2)
    x = rho * np.cos(PHI)
    y = rho * np.sin(PHI)

    rA_safe = np.where(rA < 1e-300, 1e-300, rA)
    rB_safe = np.where(rB < 1e-300, 1e-300, rB)
    thetaA = np.arccos(np.clip(z / rA_safe, -1.0, 1.0))
    thetaB = np.arccos(np.clip((z - R) / rB_safe, -1.0, 1.0))
    # rA, rB, thetaA, thetaB carry no phi-dependence (axial symmetry); PHI
    # itself (shape (1,1,n_phi)) broadcasts against them wherever combined.
    phiA = PHI
    phiB = PHI

    jac = (R / 2.0) ** 3 * (XI ** 2 - ETA ** 2)
    weight = WXI * WETA * WPHI * jac
    return dict(rA=rA, rB=rB, thetaA=thetaA, thetaB=thetaB,
                phiA=phiA, phiB=phiB, weight=weight, p=p, q=q, x=x, y=y, z=z)


def _frame_prolate(orb: Orbital, grid: dict):
    c = orb.center.upper() if isinstance(orb.center, str) else None
    if c == "A":
        return grid["rA"], grid["thetaA"], grid["phiA"]
    if c == "B":
        return grid["rB"], grid["thetaB"], grid["phiB"]
    raise ValueError("prolate engine requires center 'A' or 'B'")


def _integral_prolate(kind: str, oi: Orbital, oj: Orbital, R: float,
                       ZA: float = 1.0, ZB: float = 1.0,
                       n_xi: int = 40, n_eta: int = 100, n_phi: int = 16):
    grid = _prolate_grid(oi, oj, R, n_xi, n_eta, n_phi)
    ri, ti, pi_ = _frame_prolate(oi, grid)
    rj, tj, pj = _frame_prolate(oj, grid)
    chi_i = chi_value(oi, ri, ti, pi_)

    if kind == "S":
        f = np.conj(chi_i) * chi_value(oj, rj, tj, pj)
    elif kind == "T":
        f = np.conj(chi_i) * (-0.5) * lap_chi_value(oj, rj, tj, pj)
    elif kind == "V":
        chi_j = chi_value(oj, rj, tj, pj)
        inv_r = ZA / np.where(grid["rA"] < 1e-300, 1e-300, grid["rA"]) \
            + ZB / np.where(grid["rB"] < 1e-300, 1e-300, grid["rB"])
        f = np.conj(chi_i) * (-inv_r) * chi_j
    elif kind in ("invA", "invB"):
        chi_j = chi_value(oj, rj, tj, pj)
        r_which = grid["rA"] if kind == "invA" else grid["rB"]
        f = np.conj(chi_i) * chi_value(oj, rj, tj, pj) / np.where(
            r_which < 1e-300, 1e-300, r_which)
    else:
        raise ValueError(kind)
    return complex(np.sum(f * grid["weight"]))


# ------------------------------------------- Engine 2: spherical @ midpoint


def _spherical_grid(R: float, n_r: int = 220, n_ct: int = 160, n_phi: int = 16,
                     L: float = 1.6, origin: np.ndarray | None = None):
    if origin is None:
        origin = np.array([0.0, 0.0, R / 2.0])
    t, wt = _gauss_legendre(n_r)
    r = L * (1.0 + t) / (1.0 - t)
    dr_dt = 2.0 * L / (1.0 - t) ** 2
    w_r = wt * dr_dt

    ct, w_ct = _gauss_legendre(n_ct)
    phi, w_phi = _phi_grid(n_phi)

    Rg = r[:, None, None]
    CT = ct[None, :, None]
    PH = phi[None, None, :]
    WR = w_r[:, None, None]
    WCT = w_ct[None, :, None]
    WPH = w_phi[None, None, :]

    ST = np.sqrt(np.clip(1.0 - CT ** 2, 0.0, None))
    x = Rg * ST * np.cos(PH) + origin[0]
    y = Rg * ST * np.sin(PH) + origin[1]
    z = Rg * CT + origin[2]
    weight = WR * WCT * WPH * (Rg ** 2)
    return dict(x=x, y=y, z=z, weight=weight, r=Rg)


def _frame_xyz(orb: Orbital, grid: dict, R: float):
    cx, cy, cz = center_xyz(orb.center, R)
    dx, dy, dz = grid["x"] - cx, grid["y"] - cy, grid["z"] - cz
    r = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
    r_safe = np.where(r < 1e-300, 1e-300, r)
    theta = np.arccos(np.clip(dz / r_safe, -1.0, 1.0))
    phi = np.arctan2(dy, dx)
    return r, theta, phi


def _integral_spherical(kind: str, oi: Orbital, oj: Orbital, R: float,
                         ZA: float = 1.0, ZB: float = 1.0,
                         n_r: int = 220, n_ct: int = 160, n_phi: int = 16,
                         L: float = 1.6):
    """S, T: single spherical grid about the bond MIDPOINT (integrand is
    smooth everywhere, the algebraic radial map resolves it easily).

    V: the -Z_C/r_C kernel has an (integrable, but locally sharp) 1/r
    singularity exactly AT nucleus C, off-center from the midpoint grid, so
    that grid under-resolves it. Instead V is split additively,

        <i|-ZA/rA - ZB/rB|j> = <i|-ZA/rA|j> + <i|-ZB/rB|j>,

    and each piece is evaluated on its OWN spherical grid centered exactly
    at the relevant nucleus (a two-center, i.e. Becke-style, superposition
    of atom-centered grids) -- the natural radial variable of that grid IS
    r_C, so r_C^2 (from d^3r) cancels the 1/r_C singularity exactly, and the
    algebraic radial map converges fast again.
    """
    if kind in ("S", "T"):
        grid = _spherical_grid(R, n_r, n_ct, n_phi, L)
        ri, ti, pi_ = _frame_xyz(oi, grid, R)
        rj, tj, pj = _frame_xyz(oj, grid, R)
        chi_i = chi_value(oi, ri, ti, pi_)
        if kind == "S":
            f = np.conj(chi_i) * chi_value(oj, rj, tj, pj)
        else:
            f = np.conj(chi_i) * (-0.5) * lap_chi_value(oj, rj, tj, pj)
        return complex(np.sum(f * grid["weight"]))

    if kind == "V":
        total = 0j
        for center_label, Z in (("A", ZA), ("B", ZB)):
            origin = center_xyz(center_label, R)
            grid = _spherical_grid(R, n_r, n_ct, n_phi, L, origin=origin)
            ri, ti, pi_ = _frame_xyz(oi, grid, R)
            rj, tj, pj = _frame_xyz(oj, grid, R)
            chi_i = chi_value(oi, ri, ti, pi_)
            chi_j = chi_value(oj, rj, tj, pj)
            r_C = grid["r"]                      # exact distance from `origin`
            f = np.conj(chi_i) * (-Z / r_C) * chi_j
            total += complex(np.sum(f * grid["weight"]))
        return total

    raise ValueError(kind)


# ------------------------------------------------------------- public API


def s_ref(oi: Orbital, oj: Orbital, R: float, engine: str = "prolate",
          **kwargs) -> complex:
    """<chi_i | chi_j>."""
    if engine == "prolate":
        return _integral_prolate("S", oi, oj, R, **kwargs)
    if engine == "spherical":
        return _integral_spherical("S", oi, oj, R, **kwargs)
    raise ValueError(engine)


def t_ref(oi: Orbital, oj: Orbital, R: float, engine: str = "prolate",
          **kwargs) -> complex:
    """<chi_i | -1/2 grad^2 | chi_j>."""
    if engine == "prolate":
        return _integral_prolate("T", oi, oj, R, **kwargs)
    if engine == "spherical":
        return _integral_spherical("T", oi, oj, R, **kwargs)
    raise ValueError(engine)


def v_ref(oi: Orbital, oj: Orbital, R: float, engine: str = "prolate",
          ZA: float = 1.0, ZB: float = 1.0, **kwargs) -> complex:
    """<chi_i | (-Z_A/r_A - Z_B/r_B) | chi_j>  (H2 default: Z_A = Z_B = 1)."""
    if engine == "prolate":
        return _integral_prolate("V", oi, oj, R, ZA=ZA, ZB=ZB, **kwargs)
    if engine == "spherical":
        return _integral_spherical("V", oi, oj, R, ZA=ZA, ZB=ZB, **kwargs)
    raise ValueError(engine)


def h_core_ref(oi: Orbital, oj: Orbital, R: float, engine: str = "prolate",
               ZA: float = 1.0, ZB: float = 1.0, **kwargs) -> complex:
    return (t_ref(oi, oj, R, engine=engine, **kwargs)
            + v_ref(oi, oj, R, engine=engine, ZA=ZA, ZB=ZB, **kwargs))


# ---------------------------------------------- 1s-1s (l=0) closed forms
#
# Standard Slater-orbital results, written out here so the numerical engines
# can be checked against them directly (not just against geovac.qfd_core).
# chi_a(r) = sqrt(a^3/pi) e^{-a r}  (== N_a Y_00 with N_a = 2 a^{3/2}).


def overlap_1s1s_same_center_closed(a: float, b: float) -> float:
    """<1s_a | 1s_b>, both on the SAME center: (2 sqrt(ab)/(a+b))^3."""
    return (2.0 * math.sqrt(a * b) / (a + b)) ** 3


def overlap_1s1s_two_center_equal_zeta_closed(zeta: float, R: float) -> float:
    """<1s_A(zeta) | 1s_B(zeta)>, equal exponents, two centers separated by R:
    S(zeta,R) = e^{-zeta R} (1 + zeta R + (zeta R)^2 / 3)."""
    x = zeta * R
    return math.exp(-x) * (1.0 + x + x * x / 3.0)


def kinetic_1s1s_same_center_closed(a: float, b: float) -> float:
    """<1s_a | -1/2 grad^2 | 1s_b>, same center: T = (a b / 2) * S,
    S = overlap_1s1s_same_center_closed(a, b).  Derivation: with
    N_x = 2 x^{3/2}, direct integration of
        -1/2 N_a N_b int_0^inf (b^2 r^2 - 2 b r) e^{-(a+b) r} dr
    gives  T = N_a N_b a b / (a+b)^3 = (a b / 2) * [N_a N_b * 2/(a+b)^3] = (ab/2) S.
    """
    S = overlap_1s1s_same_center_closed(a, b)
    return 0.5 * a * b * S


def v_1s1s_same_center_own_kernel_closed(a: float, b: float) -> float:
    """<1s_a | -1/r_A | 1s_b>, both 1s ORBITALS on A, kernel ALSO on A (the
    genuine one-center Coulomb integral): -N_a N_b / (a+b)^2, N_x = 2 x^{3/2}."""
    Na, Nb = 2.0 * a ** 1.5, 2.0 * b ** 1.5
    return -(Na * Nb) / (a + b) ** 2


__all__ = [
    "Orbital", "center_xyz", "sto_norm", "Ylm", "chi_value", "lap_chi_value",
    "s_ref", "t_ref", "v_ref", "h_core_ref",
    "overlap_1s1s_same_center_closed",
    "overlap_1s1s_two_center_equal_zeta_closed",
    "kinetic_1s1s_same_center_closed",
    "v_1s1s_same_center_own_kernel_closed",
]


# ===========================================================================
#                                VALIDATION
# ===========================================================================

if __name__ == "__main__":
    import sys
    import time

    t0 = time.time()
    rows = []          # (label, computed, reference, abs_err, tol, pass)

    def check(label, computed, reference, tol):
        err = abs(computed - reference)
        ok = err < tol
        rows.append((label, computed, reference, err, tol, ok))
        return ok

    # -------------------------------------------------------------- part 0
    # Finite-difference confirmation of the analytic Laplacian identity
    # (independent of BOTH quadrature engines -- pure calculus check).
    print("=" * 78)
    print("PART 0 -- finite-difference check of the analytic STO Laplacian identity")
    print("=" * 78)
    rng = np.random.default_rng(0)
    fd_h = 1e-3
    fd_pass = True
    for l, m in [(0, 0), (1, 0), (1, 1), (1, -1)]:
        orb = Orbital(center="A", zeta=1.15, l=l, m=m)
        for _ in range(4):
            # sample a point away from the origin so r, theta stay regular
            vec = rng.normal(size=3)
            vec = vec / np.linalg.norm(vec) * (0.6 + 1.2 * rng.random())
            x0, y0, z0 = vec

            def chi_cart(x, y, z):
                r = math.sqrt(x * x + y * y + z * z)
                theta = math.acos(max(-1.0, min(1.0, z / r)))
                phi = math.atan2(y, x)
                return complex(chi_value(orb, np.array(r), np.array(theta),
                                          np.array(phi)))

            # 5-point central stencil per axis, O(h^4)
            def d2(axis):
                base = list((x0, y0, z0))
                vals = []
                for k in (-2, -1, 0, 1, 2):
                    p = base.copy()
                    p[axis] += k * fd_h
                    vals.append(chi_cart(*p))
                return (-vals[4] + 16 * vals[3] - 30 * vals[2] + 16 * vals[1]
                        - vals[0]) / (12 * fd_h ** 2)

            lap_fd = d2(0) + d2(1) + d2(2)
            r0 = math.sqrt(x0 ** 2 + y0 ** 2 + z0 ** 2)
            th0 = math.acos(z0 / r0)
            ph0 = math.atan2(y0, x0)
            lap_an = complex(lap_chi_value(orb, np.array(r0), np.array(th0),
                                            np.array(ph0)))
            rel = abs(lap_fd - lap_an) / max(abs(lap_an), 1e-12)
            ok = rel < 1e-4
            fd_pass &= ok
            print(f"  l={l:d} m={m:2d}  r0={r0:.3f}  "
                  f"lap_analytic={lap_an.real:+.6f}{lap_an.imag:+.6f}j  "
                  f"lap_FD={lap_fd.real:+.6f}{lap_fd.imag:+.6f}j  "
                  f"rel_err={rel:.2e}  {'OK' if ok else 'FAIL'}")
    print(f"  -> Laplacian identity finite-difference check: "
          f"{'PASS' if fd_pass else 'FAIL'}\n")

    # -------------------------------------------------------------- part 1
    # 1s-1s (l=0) vs hand-derived closed forms
    print("=" * 78)
    print("PART 1 -- 1s-1s (l=0) vs closed forms")
    print("=" * 78)

    a, b = 1.0, 1.2
    oA_a = Orbital("A", a, 0, 0)
    oA_b = Orbital("A", b, 0, 0)
    S_same_ref = float(s_ref(oA_a, oA_b, 1.0).real)   # R irrelevant, same center
    S_same_cf = overlap_1s1s_same_center_closed(a, b)
    check(f"S same-center 1s(a={a})-1s(b={b})", S_same_ref, S_same_cf, 1e-10)

    T_same_ref = float(t_ref(oA_a, oA_b, 1.0).real)
    T_same_cf = kinetic_1s1s_same_center_closed(a, b)
    check(f"T same-center 1s(a={a})-1s(b={b})", T_same_ref, T_same_cf, 1e-10)

    V_ownA_ref = float(_integral_prolate("invA", oA_a, oA_b, 1.0).real)
    # closed form is <1s_a|-1/r_A|1s_b>; invA computes the +1/r_A kernel, so negate
    V_ownA_cf = -v_1s1s_same_center_own_kernel_closed(a, b)
    check(f"<1s_a|1/r_A|1s_b> same-center own-kernel (a={a},b={b})",
          V_ownA_ref, V_ownA_cf, 1e-10)

    zeta, R0 = 1.0, 1.4
    oA_z = Orbital("A", zeta, 0, 0)
    oB_z = Orbital("B", zeta, 0, 0)
    S_two_ref = float(s_ref(oA_z, oB_z, R0).real)
    S_two_cf = overlap_1s1s_two_center_equal_zeta_closed(zeta, R0)
    check(f"S two-center equal-zeta 1s-1s (zeta={zeta}, R={R0})",
          S_two_ref, S_two_cf, 1e-10)

    print(f"{'label':55s} {'computed':>16s} {'reference':>16s} "
          f"{'abs_err':>10s}  status")
    for label, comp, refv, err, tol, ok in rows:
        print(f"{label:55s} {comp:16.10f} {refv:16.10f} {err:10.2e}  "
              f"{'PASS' if ok else 'FAIL'}")
    part1_pass = all(r[5] for r in rows)
    print(f"  -> Part 1 (closed forms): {'PASS' if part1_pass else 'FAIL'}\n")

    # -------------------------------------------------------------- part 2
    # l=0 vs geovac.qfd_core (the ONLY place this file imports geovac)
    print("=" * 78)
    print("PART 2 -- l=0 vs geovac.qfd_core (overlap / kinetic / h_core)")
    print("=" * 78)
    from fractions import Fraction

    import geovac.qfd_core as qfd

    rows2 = []
    cases = [(1.0, 1.0, 1.4), (1.0, 1.3, 1.4), (Fraction(9, 10), Fraction(11, 10), 2.0)]
    for zA, zB, R in cases:
        zA_f, zB_f = float(zA), float(zB)
        oi = Orbital("A", zA_f, 0, 0)
        oj = Orbital("B", zB_f, 0, 0)
        oi_q = ("A", Fraction(zA) if not isinstance(zA, Fraction) else zA, 1)
        oj_q = ("B", Fraction(zB) if not isinstance(zB, Fraction) else zB, 1)

        S_num = float(s_ref(oi, oj, R).real)
        S_qfd = float(qfd.overlap(oi_q, oj_q, R))
        rows2.append((f"S(zA={zA_f},zB={zB_f},R={R})", S_num, S_qfd,
                      abs(S_num - S_qfd), 1e-8))

        T_num = float(t_ref(oi, oj, R).real)
        T_qfd = float(qfd.kinetic(oi_q, oj_q, R))
        rows2.append((f"T(zA={zA_f},zB={zB_f},R={R})", T_num, T_qfd,
                      abs(T_num - T_qfd), 1e-8))

        H_num = float(h_core_ref(oi, oj, R).real)
        H_qfd = float(qfd.h_core(oi_q, oj_q, 1, 1, R))
        rows2.append((f"h_core(zA={zA_f},zB={zB_f},R={R})", H_num, H_qfd,
                      abs(H_num - H_qfd), 1e-8))

        # same-center s-only too
        oiA = Orbital("A", zA_f, 0, 0)
        ojA = Orbital("A", zB_f, 0, 0)
        oiA_q, ojA_q = ("A", oi_q[1], 1), ("A", oj_q[1], 1)
        S_num_sc = float(s_ref(oiA, ojA, R).real)
        S_qfd_sc = float(qfd.overlap(oiA_q, ojA_q, R))
        rows2.append((f"S same-center (a={zA_f},b={zB_f})", S_num_sc, S_qfd_sc,
                      abs(S_num_sc - S_qfd_sc), 1e-8))

    print(f"{'label':40s} {'oracle':>16s} {'qfd_core':>16s} {'abs_err':>10s}  status")
    for label, num, ref, err, tol in rows2:
        ok = err < tol
        print(f"{label:40s} {num:16.10f} {ref:16.10f} {err:10.2e}  "
              f"{'PASS' if ok else 'FAIL'}")
    part2_pass = all(err < tol for _, _, _, err, tol in rows2)
    print(f"  -> Part 2 (vs qfd_core, l=0): {'PASS' if part2_pass else 'FAIL'}\n")

    # -------------------------------------------------------------- part 3
    # p functions: Engine 1 (prolate) vs Engine 2 (spherical @ midpoint)
    print("=" * 78)
    print("PART 3 -- p functions (l=1): Engine 1 (prolate) vs Engine 2 (spherical)")
    print("=" * 78)
    rows3 = []
    R1 = 1.4
    p_cases = [
        ("s(A)-pz(B)", Orbital("A", 1.0, 0, 0), Orbital("B", 1.05, 1, 0)),
        ("pz(A)-pz(B)", Orbital("A", 1.0, 1, 0), Orbital("B", 1.0, 1, 0)),
        ("pz(A)-pz(B) unequal zeta", Orbital("A", 0.9, 1, 0), Orbital("B", 1.2, 1, 0)),
        ("p+1(A)-p+1(B)", Orbital("A", 1.0, 1, 1), Orbital("B", 1.0, 1, 1)),
        ("p+1(A)-p-1(B)", Orbital("A", 1.0, 1, 1), Orbital("B", 1.0, 1, -1)),
        ("p0(A)-p+1(B) [m mismatch]", Orbital("A", 1.0, 1, 0), Orbital("B", 1.0, 1, 1)),
        ("s(A)-p+1(A) same-center", Orbital("A", 1.1, 0, 0), Orbital("A", 1.1, 1, 1)),
    ]
    for label, oi, oj in p_cases:
        for kind, fn in (("S", s_ref), ("T", t_ref), ("V", v_ref)):
            v1 = fn(oi, oj, R1, engine="prolate")
            v2 = fn(oi, oj, R1, engine="spherical")
            err = abs(v1 - v2)
            rows3.append((f"{kind} {label}", v1, v2, err, 1e-6))

    print(f"{'label':38s} {'prolate':>28s} {'spherical':>28s} {'abs_err':>10s}  status")
    for label, v1, v2, err, tol in rows3:
        ok = err < tol
        print(f"{label:38s} {v1.real:+.8f}{v1.imag:+.8f}j "
              f"{v2.real:+.8f}{v2.imag:+.8f}j {err:10.2e}  "
              f"{'PASS' if ok else 'FAIL'}")
    part3_pass = all(err < tol for _, _, _, err, tol in rows3)
    print(f"  -> Part 3 (dual-engine p-function agreement): "
          f"{'PASS' if part3_pass else 'FAIL'}\n")

    # -------------------------------------------------------------- part 4
    # structural sanity: normalization, positivity, hermiticity
    print("=" * 78)
    print("PART 4 -- structural sanity checks")
    print("=" * 78)
    rows4 = []
    sanity_orbs = [
        Orbital("A", 1.0, 0, 0), Orbital("A", 1.1, 1, 0),
        Orbital("A", 0.95, 1, 1), Orbital("B", 1.2, 1, -1),
    ]
    for orb in sanity_orbs:
        Sii = s_ref(orb, orb, R1)
        Tii = t_ref(orb, orb, R1)
        rows4.append((f"S_ii=1 ({orb.center},zeta={orb.zeta},l={orb.l},m={orb.m})",
                      Sii.real, 1.0, abs(Sii.real - 1.0), 1e-8))
        rows4.append((f"Im(S_ii)=0", Sii.imag, 0.0, abs(Sii.imag), 1e-10))
        rows4.append((f"T_ii>0 ({orb.center},zeta={orb.zeta},l={orb.l},m={orb.m})",
                      Tii.real, abs(Tii.real), 0.0, 0.0))

    # hermiticity: S(i,j,R) == conj(S(j,i,R))
    oi, oj = Orbital("A", 1.0, 1, 1), Orbital("B", 1.1, 1, -1)
    Sij = s_ref(oi, oj, R1)
    Sji = s_ref(oj, oi, R1)
    rows4.append(("Hermiticity S_ij = conj(S_ji)", Sij, np.conj(Sji),
                  abs(Sij - np.conj(Sji)), 1e-9))

    print(f"{'label':55s} {'value':>16s}  status")
    sanity_pass = True
    for row in rows4:
        label, val, refv, err, tol = row
        if label.startswith("T_ii>0"):
            ok = val > 0
            print(f"{label:55s} {val:16.10f}  {'PASS' if ok else 'FAIL'}")
        else:
            ok = err < tol
            vstr = f"{val:.3e}" if isinstance(val, complex) else f"{val:16.10f}"
            print(f"{label:55s} {str(val):>16s}  {'PASS' if ok else 'FAIL'}")
        sanity_pass &= ok
    print(f"  -> Part 4 (structural sanity): {'PASS' if sanity_pass else 'FAIL'}\n")

    # -------------------------------------------------------------- part 5
    # convergence check: double the prolate-engine grid density, confirm stable
    print("=" * 78)
    print("PART 5 -- grid-convergence check (prolate engine, doubled density)")
    print("=" * 78)
    oi, oj = Orbital("A", 1.0, 1, 0), Orbital("B", 1.0, 1, 0)
    lo = dict(n_xi=24, n_eta=60, n_phi=8)
    hi = dict(n_xi=48, n_eta=160, n_phi=16)
    conv_pass = True
    for kind, fn in (("S", s_ref), ("T", t_ref), ("V", v_ref)):
        v_lo = fn(oi, oj, R1, engine="prolate", **lo)
        v_hi = fn(oi, oj, R1, engine="prolate", **hi)
        err = abs(v_lo - v_hi)
        ok = err < 1e-8
        conv_pass &= ok
        print(f"  {kind}: lo-grid={v_lo.real:+.12f}  hi-grid={v_hi.real:+.12f}  "
              f"diff={err:.2e}  {'PASS' if ok else 'FAIL'}")
    print(f"  -> Part 5 (grid convergence): {'PASS' if conv_pass else 'FAIL'}\n")

    # -------------------------------------------------------------- summary
    overall = fd_pass and part1_pass and part2_pass and part3_pass and sanity_pass and conv_pass
    print("=" * 78)
    print(f"OVERALL: {'PASS' if overall else 'FAIL'}   "
          f"(elapsed {time.time() - t0:.1f}s)")
    print("=" * 78)
    sys.exit(0 if overall else 1)
