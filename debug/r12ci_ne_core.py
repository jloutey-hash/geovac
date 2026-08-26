"""N-electron s-only correlated-CI primitives (step 3): pair kernels + leg contractions.

An s-only N-electron matrix element of a product of pair functions reduces to a
contraction of per-electron radial densities d_k(r) = R_bra(r) R_ket(r) r^2 dr
against pair kernels.  Verified shapes (steps 1-2):

  0 legs                      trivial product of radial sums
  1 leg   (a,b)               d_a @ K0 @ d_b
  2 legs  same pair (a,b)     d_a @ K0[phi*psi] @ d_b   (kernel of the PRODUCT --
                                                         NOT the product of kernels)
  2 legs  shared vertex a     sum_ra d_a [K1 @ d_b][K2 @ d_c]     (RULE A / RULE B)
  3 legs  triangle            sum_L (2L+1)^-2 multipole contraction (step 1)

Self-tests at the bottom check every shape against direct 3-electron quadrature,
including a deliberate WRONG variant that must visibly fail.
"""
from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import leggauss, legval


# ---------------------------------------------------------------------------
# pair kernels
# ---------------------------------------------------------------------------
def _u_quad(r, nx):
    """Gauss-Legendre nodes/weights for the substitution x -> u = r_12.

    (1/2) INT_-1^1 g(r_12) dx  =  (1/(2 r_i r_j)) INT_{|ri-rj|}^{ri+rj} g(u) u du

    because dx = -u du/(r_i r_j).  The Jacobian's factor of u CANCELS a 1/r_12 in the
    integrand, so Coulomb-containing kernels become smooth integrals.  Gauss-Legendre
    in x converges only as ~1/nx on 1/r_12 (2.7e-3 rel at nx=160); in u it is exact
    for polynomial g and spectrally accurate otherwise.
    """
    R1, R2 = np.meshgrid(r, r, indexing="ij")
    a = np.abs(R1 - R2)
    b = R1 + R2
    t, w = leggauss(nx)
    half = (b - a) / 2.0
    u = a[None, :, :] + half[None, :, :] * (t[:, None, None] + 1.0)
    jac = half[None, :, :] * w[:, None, None]
    pref = 1.0 / (2.0 * R1 * R2)
    return R1, R2, u, jac, pref


def l0_kernel(fun, r, nx: int = 400) -> np.ndarray:
    """K0[i,j] = (1/2) INT_-1^1 fun(r_ij) dx, via the u = r_12 substitution."""
    R1, R2, u, jac, pref = _u_quad(r, nx)
    return pref * np.einsum("kij->ij", fun(u) * u * jac)


def l0_kernel_proj(fun, r, nx: int = 400) -> np.ndarray:
    """kA-type kernel: <fun(r_ab) * (rhat_ab . rhat_a)>, projected onto the FIRST index.

    rhat_ab . rhat_a = (r_a - r_b x)/r_ab, and in the u variable
    r_a - r_b x = (r_a^2 - r_b^2 + u^2) / (2 r_a), so the 1/u cancels the Jacobian's u.
    """
    R1, R2, u, jac, pref = _u_quad(r, nx)
    proj = (R1[None, :, :] ** 2 - R2[None, :, :] ** 2 + u * u) / (2.0 * R1[None, :, :])
    return pref * np.einsum("kij->ij", fun(u) * proj * jac)


def moments(fun, r, Lmax: int, nx: int = 400) -> np.ndarray:
    """a_L[i,j] with fun(r_ij) = sum_L a_L P_L(cos theta_ij), via the u substitution."""
    R1, R2, u, jac, pref = _u_quad(r, nx)
    x = (R1[None, :, :] ** 2 + R2[None, :, :] ** 2 - u * u) / (2.0 * R1[None, :, :] * R2[None, :, :])
    x = np.clip(x, -1.0, 1.0)
    fv = fun(u) * u * jac
    out = np.zeros((Lmax + 1, r.size, r.size))
    for L in range(Lmax + 1):
        cf = np.zeros(L + 1)
        cf[L] = 1.0
        out[L] = (2 * L + 1) / 2.0 * 2.0 * pref * np.einsum("kij->ij", fv * legval(x, cf))
    return out


def coul_moments(r, Lmax: int) -> np.ndarray:
    """EXACT Legendre moments of 1/r_12: a_L = r_<^L / r_>^(L+1). No quadrature."""
    R1, R2 = np.meshgrid(r, r, indexing="ij")
    lo = np.minimum(R1, R2)
    hi = np.maximum(R1, R2)
    return np.array([lo ** L / hi ** (L + 1) for L in range(Lmax + 1)])


# ---------------------------------------------------------------------------
# leg contractions. dens = list of per-electron radial density vectors.
# ---------------------------------------------------------------------------
def c_none(dens) -> float:
    out = 1.0
    for d in dens:
        out *= d.sum()
    return float(out)


def c_one(K, a: int, b: int, dens) -> float:
    val = dens[a] @ K @ dens[b]
    for k, d in enumerate(dens):
        if k not in (a, b):
            val *= d.sum()
    return float(val)


def c_shared(K1, K2, a: int, b: int, c: int, dens) -> float:
    """Legs (a,b) and (a,c) sharing vertex a. RULE A (scalar) / RULE B (vec-vec)."""
    val = dens[a] @ ((K1 @ dens[b]) * (K2 @ dens[c]))
    for k, d in enumerate(dens):
        if k not in (a, b, c):
            val *= d.sum()
    return float(val)


def c_triangle(mA, mB, mC, a: int, b: int, c: int, dens) -> float:
    """Legs (a,b)=mA, (a,c)=mB, (b,c)=mC. sum_L (2L+1)^-2 contraction."""
    tot = 0.0
    da, db, dc = dens[a], dens[b], dens[c]
    for L in range(mA.shape[0]):
        # P(r_a, r_b) = sum_rc  mB[L](ra,rc) * dc(rc) * mC[L](rb,rc)
        P = (mB[L] * dc[None, :]) @ mC[L].T
        tot += (da @ ((mA[L] * P) @ db)) / (2 * L + 1) ** 2
    for k, d in enumerate(dens):
        if k not in (a, b, c):
            tot *= d.sum()
    return float(tot)


# ---------------------------------------------------------------------------
# self-test reference: direct 3-electron orientation quadrature
# ---------------------------------------------------------------------------
def brute3(legs, r, dens, nang: int = 48) -> float:
    """legs = list of (a, b, callable). Fix r1hat = zhat, integrate Omega2, Omega3.

    Only the (1,2) leg depends on phi3, and it depends on (r2, r3) alone, so the phi
    sum collapses to an Ng x Ng array before any Ng^3 work.  Cost is
    nang^2 * (nang*Ng^2 + Ng^3) instead of nang^3 * Ng^3.
    """
    xt, wt = leggauss(nang)
    ph, wp = leggauss(nang)
    phi3 = np.pi * (ph + 1.0)
    wphi = np.pi * wp

    def pairfun(a, b, x):
        """product of all legs on pair (a,b) evaluated at cos = x (scalar or array)."""
        ra = r[:, None] if (a, b) != (1, 2) else r[:, None]
        rb = r[None, :]
        out = None
        for (p, q, fn) in legs:
            if (min(p, q), max(p, q)) != (a, b):
                continue
            rab = np.sqrt(np.maximum(ra ** 2 + rb ** 2 - 2 * ra * rb * x, 1e-30))
            out = fn(rab) if out is None else out * fn(rab)
        return out

    has01 = any((min(p, q), max(p, q)) == (0, 1) for p, q, _ in legs)
    has02 = any((min(p, q), max(p, q)) == (0, 2) for p, q, _ in legs)
    has12 = any((min(p, q), max(p, q)) == (1, 2) for p, q, _ in legs)
    d1, d2, d3 = dens
    ones = np.ones((r.size, r.size))
    tot = 0.0
    for c2, w2 in zip(xt, wt):
        s2 = np.sqrt(1 - c2 * c2)
        A01 = pairfun(0, 1, c2) if has01 else ones
        for c3, w3 in zip(xt, wt):
            s3 = np.sqrt(1 - c3 * c3)
            A02 = pairfun(0, 2, c3) if has02 else ones
            if has12:
                Cint = np.zeros((r.size, r.size))
                for p3, wpp in zip(phi3, wphi):
                    Cint += wpp * pairfun(1, 2, c2 * c3 + s2 * s3 * np.cos(p3))
            else:
                Cint = ones * (2 * np.pi)
            # M(r2,r3) = sum_r1 d1(r1) A01(r1,r2) A02(r1,r3)
            M = (A01 * d1[:, None]).T @ A02
            tot += w2 * w3 * float(d2 @ ((M * Cint) @ d3))
    return tot / (2.0 * 2.0 * 2 * np.pi)


if __name__ == "__main__":
    Ng = 22
    r = np.linspace(0.15, 4.5, Ng)
    dens = [np.exp(-1.3 * r) * r ** 2,
            np.exp(-0.9 * r) * r ** 2 * (1 + 0.3 * r),
            np.exp(-1.7 * r) * r ** 2]
    f = lambda x: np.exp(-0.7 * x)
    g = lambda x: np.exp(-1.3 * x)
    coul = lambda x: 1.0 / x

    print("=" * 74)
    print("leg-contraction primitives vs direct 3-electron quadrature")
    print("=" * 74)
    print(f"{'shape':<28}{'contraction':>16}{'brute':>16}{'rel.diff':>12}")

    v = c_one(l0_kernel(f, r), 0, 1, dens)
    b = brute3([(0, 1, f)], r, dens)
    print(f"{'1 leg (0,1)':<28}{v:>16.10f}{b:>16.10f}{abs(v - b) / abs(b):>12.2e}")

    v = c_one(l0_kernel(lambda x: f(x) * g(x), r), 0, 1, dens)
    b = brute3([(0, 1, f), (0, 1, g)], r, dens)
    print(f"{'2 legs same pair':<28}{v:>16.10f}{b:>16.10f}{abs(v - b) / abs(b):>12.2e}")

    vbad = c_one(l0_kernel(f, r) * l0_kernel(g, r), 0, 1, dens)
    print(f"{'  WRONG product-of-kernels':<28}{vbad:>16.10f}{b:>16.10f}"
          f"{abs(vbad - b) / abs(b):>12.2e}  <- must be LARGE")

    v = c_shared(l0_kernel(f, r), l0_kernel(g, r), 0, 1, 2, dens)
    b = brute3([(0, 1, f), (0, 2, g)], r, dens)
    print(f"{'2 legs shared (scalar)':<28}{v:>16.10f}{b:>16.10f}{abs(v - b) / abs(b):>12.2e}")

    Lmax = 16
    v = c_triangle(moments(f, r, Lmax), moments(g, r, Lmax), coul_moments(r, Lmax),
                   0, 1, 2, dens)
    b = brute3([(0, 1, f), (0, 2, g), (1, 2, coul)], r, dens, nang=64)
    print(f"{'3 legs triangle (f,g,1/r)':<28}{v:>16.10f}{b:>16.10f}"
          f"{abs(v - b) / abs(b):>12.2e}")

    h = lambda x: np.exp(-0.5 * x)
    v = c_triangle(moments(f, r, Lmax), moments(g, r, Lmax), moments(h, r, Lmax),
                   0, 1, 2, dens)
    b = brute3([(0, 1, f), (0, 2, g), (1, 2, h)], r, dens, nang=64)
    print(f"{'3 legs triangle (smooth)':<28}{v:>16.10f}{b:>16.10f}"
          f"{abs(v - b) / abs(b):>12.2e}")
