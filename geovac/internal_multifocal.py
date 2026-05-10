"""
Internal Multi-Focal Architecture for Atomic Excited States
============================================================

Multi-electron CI on a hydrogenic basis where each orbital carries its
own focal length (Sturmian exponent) lambda_p. This module is the
internal-atom analog of `geovac/cross_register_vne.py` (which implements
the cross-register multi-focal architecture for chemistry observables).

Background (Track 4 He oscillator strength autopsy, 2026-05-09)
---------------------------------------------------------------

Single-exponent Sturmian-style basis at k = Z = 2 reproduces the Helium
2^1P -> 1^1S oscillator strength at f ~ 0.444 vs Drake 0.276 (+61%
residual, plateaued at n_max = 7). The 2^1P excited state is
variationally violated by 0.08 Ha — same kappa-induced over-binding
mechanism as the small-Z graph-validity-boundary at He. The Lyman alpha
sanity check passes to machine precision, so the angular machinery is
exact; the failure is in the radial / CI sector (single common Sturmian
exponent cannot represent both compact 1s and diffuse 2p simultaneously).

This module addresses that limitation by allowing per-orbital focal
lengths lambda_p, with cross-exponent overlap, dipole, kinetic, and
two-electron Slater integrals all evaluated in closed form (sums of
incomplete-gamma terms via the same machinery
``shibuya_wulfman._split_integral_analytical`` already uses for cross-
center V_ne).

Architecture
------------

* **MultifocalSpec**: dataclass listing the orbitals as (n, l, m, lambda)
  tuples plus nuclear charge.
* **Cross-exponent radial primitives**: overlap, ``r^k`` matrix elements,
  kinetic, two-electron Slater R^k via Gauss-Laguerre on the outer radial
  + closed-form incomplete-gamma on the inner.
* **One-body assemblies**: ``h1_multifocal`` (kinetic + V_Ne),
  ``overlap_multifocal``, ``dipole_z_multifocal``.
* **CI on multi-focal basis**: generalized eigenproblem (H, S) since the
  basis is non-orthogonal at distinct lambda.
* **Driver**: oscillator strength assembly mirroring
  ``debug/calc_track_he_oscillator_v1.py`` but on the multi-focal basis.

Key design point: when all orbitals share lambda = Z/n_orb (the natural
Coulomb-Sturmian convention), the matrix elements reduce to the standard
single-focal hydrogenic results. This is a non-negotiable regression
test (``test_slater_matched_lambda_regression``).

Honest scope (Phase B first-pass)
---------------------------------

* The graph kappa = -1/16 adjacency from ``casimir_ci._build_graph_h1``
  is NOT applied at the multi-focal level. We evaluate kinetic + V_Ne
  directly from radial integrals, giving a standard variational
  hydrogenic CI. The graph-native property is sacrificed in this
  first-pass to isolate whether the basis-architecture fix alone closes
  the multi-electron oscillator-strength residual.
* Standard L2 hydrogenic normalization (∫|R|² r² dr = 1), NOT Coulomb-
  Sturmian normalization. We are NOT pursuing Sturmian closure here;
  we are addressing the orthogonal question of whether physical screening
  exponents per orbital close the multi-electron observable.
* Generalized eigenproblem H c = E S c via ``scipy.linalg.eigh`` with
  ``b=S_matrix``. S must be positive-definite; we monitor its condition
  number.

References
----------

* Track 4 memo: ``debug/he_oscillator_v1_memo.md`` (the diagnosis this
  module addresses).
* Architectural template: ``geovac/cross_register_vne.py``.
* Shibuya-Wulfman radial machinery: ``geovac/shibuya_wulfman.py``.
* Closed-form 1s × 1s overlap at distinct lambda: standard textbook
  ``S(lam_a, lam_b) = (2 sqrt(lam_a lam_b))^3 / (lam_a + lam_b)^3``.

Author: GeoVac Development Team (Phase A/B/C internal multi-focal sprint)
Date: 2026-05-09
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from math import factorial, sqrt
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from geovac.casimir_ci import _gaunt_ck
from geovac.shibuya_wulfman import (
    _hydrogenic_poly_coeffs_lam,
    _poly_product,
    _split_integral_analytical,
)


# ---------------------------------------------------------------------------
# Specification dataclass
# ---------------------------------------------------------------------------

@dataclass
class MultifocalOrbital:
    """A single orbital on the internal multi-focal basis.

    Attributes
    ----------
    n : int
        Principal quantum number (>= 1).
    l : int
        Angular momentum quantum number (0 <= l < n).
    m : int
        Magnetic quantum number (-l <= m <= l).
    lam : float
        Focal length / Sturmian exponent (positive, in atomic units of
        bohr^{-1}). For hydrogenic Z/n convention, lam = Z/n. For
        Slater-screened convention, lam = (Z - sigma)/n_eff with
        screening sigma per Slater rules.
    label : str
        Optional human-readable label (e.g. "1s_var" or "2p_screened").
    """
    n: int
    l: int
    m: int
    lam: float
    label: str = ""

    def __post_init__(self) -> None:
        if self.n < 1:
            raise ValueError(f"n must be >= 1, got {self.n}")
        if self.l < 0 or self.l >= self.n:
            raise ValueError(f"need 0 <= l < n, got l={self.l}, n={self.n}")
        if abs(self.m) > self.l:
            raise ValueError(f"need |m| <= l, got m={self.m}, l={self.l}")
        if self.lam <= 0:
            raise ValueError(f"lam must be positive, got {self.lam}")

    def key(self) -> Tuple[int, int, int, float]:
        """Cache key for this orbital."""
        return (self.n, self.l, self.m, float(self.lam))


@dataclass
class MultifocalSpec:
    """Specification of a multi-focal CI calculation.

    Attributes
    ----------
    orbitals : list of MultifocalOrbital
        The orbitals in the basis. May contain orbitals with the same
        (n, l, m) at different lambda — this is the whole point of the
        multi-focal extension.
    Z_nuc : float
        Nuclear charge (positive).
    label : str
        Optional human-readable label.
    """
    orbitals: List[MultifocalOrbital]
    Z_nuc: float = 1.0
    label: str = ""

    def __post_init__(self) -> None:
        if self.Z_nuc <= 0:
            raise ValueError(f"Z_nuc must be positive, got {self.Z_nuc}")
        if not self.orbitals:
            raise ValueError("orbitals list must be non-empty")

    @property
    def n_orbitals(self) -> int:
        return len(self.orbitals)


# ---------------------------------------------------------------------------
# Cross-exponent radial primitives
# ---------------------------------------------------------------------------

@lru_cache(maxsize=4096)
def _radial_polys(n: int, l: int, lam: float) -> Tuple[Tuple[float, ...], float]:
    """Cached radial polynomial coefficients R_{n,l}^{lam}(r) = e^{-lam r} sum c_k r^k.

    Returns (coefficients_tuple, alpha=lam). Tuple form for hashability.
    """
    coeffs, alpha = _hydrogenic_poly_coeffs_lam(lam, n, l)
    return tuple(coeffs), float(alpha)


def overlap_radial(
    n_p: int, l_p: int, lam_p: float,
    n_q: int, l_q: int, lam_q: float,
) -> float:
    """Cross-exponent radial overlap <R_p | R_q> = int R_p R_q r^2 dr.

    For l_p != l_q, returns 0 (orthogonality is angular, not radial).
    For l_p == l_q, the integral is closed form:
        sum_{i,j} c_p[i] c_q[j] (i+j+2)! / (lam_p + lam_q)^{i+j+3}.

    At lam_p = lam_q = Z/n with n_p == n_q this returns 1 (basis
    normalization). At lam_p = lam_q = Z/n with n_p != n_q this returns
    0 (standard hydrogenic orthogonality at matched exponent).

    Parameters
    ----------
    n_p, l_p : int
        Bra orbital quantum numbers.
    lam_p : float
        Bra orbital focal length.
    n_q, l_q : int
        Ket orbital quantum numbers.
    lam_q : float
        Ket orbital focal length.

    Returns
    -------
    float
        Radial overlap (positive on the diagonal).
    """
    if l_p != l_q:
        return 0.0

    c_p, alpha_p = _radial_polys(n_p, l_p, lam_p)
    c_q, alpha_q = _radial_polys(n_q, l_q, lam_q)
    alpha_total = alpha_p + alpha_q

    result = 0.0
    for i, ci in enumerate(c_p):
        if ci == 0.0:
            continue
        for j, cj in enumerate(c_q):
            if cj == 0.0:
                continue
            power = i + j + 2  # r^2 spherical-volume Jacobian
            result += ci * cj * factorial(power) / alpha_total ** (power + 1)
    return result


def matrix_element_rk(
    n_p: int, l_p: int, lam_p: float,
    n_q: int, l_q: int, lam_q: float,
    k_pow: int,
) -> float:
    """Cross-exponent radial matrix element <R_p | r^k | R_q>.

    = int R_p(r) R_q(r) r^k * r^2 dr
    = sum_{i,j} c_p[i] c_q[j] (i+j+2+k)! / (lam_p + lam_q)^{i+j+3+k}

    Convergence: requires i+j+2+k >= 0 for all (i, j) terms. For k >= 0
    always converges; for k < 0 (e.g. <1/r>) requires the smallest power
    to satisfy 0 + 0 + 2 + k >= 1 (i.e. k >= -1). For k = -2 (matrix
    element of 1/r^2) we need l_p + l_q >= 1 for the (0, 0) term to be
    well-defined, etc.

    For l_p != l_q with same convention, this still computes the radial
    factor; the angular factor must be applied externally if needed.

    Parameters
    ----------
    n_p, l_p, lam_p, n_q, l_q, lam_q : as before.
    k_pow : int
        Power of r in the matrix element. k=0: overlap. k=1: dipole.
        k=-1: <1/r> (Coulomb attraction). k=2: <r^2> (size).

    Returns
    -------
    float
        Radial matrix element.
    """
    c_p, alpha_p = _radial_polys(n_p, l_p, lam_p)
    c_q, alpha_q = _radial_polys(n_q, l_q, lam_q)
    alpha_total = alpha_p + alpha_q

    result = 0.0
    for i, ci in enumerate(c_p):
        if ci == 0.0:
            continue
        for j, cj in enumerate(c_q):
            if cj == 0.0:
                continue
            power = i + j + 2 + k_pow
            if power < 0:
                raise ValueError(
                    f"matrix_element_rk: integrand diverges for "
                    f"(i={i}, j={j}, k={k_pow}); need power >= 0, got {power}. "
                    f"Likely l_p + l_q too small for this k."
                )
            result += ci * cj * factorial(power) / alpha_total ** (power + 1)
    return result


def kinetic_radial(
    n_p: int, l_p: int, lam_p: float,
    n_q: int, l_q: int, lam_q: float,
) -> float:
    """Cross-exponent kinetic matrix element <R_p | -1/2 nabla^2 | R_q>.

    Diagonal in l (no angular mixing from kinetic). Returns 0 for l_p != l_q.

    The radial kinetic operator is
        T_l = -1/2 * (1/r^2) d/dr (r^2 d/dr) + l(l+1)/(2 r^2)
            = -1/2 d^2/dr^2 - (1/r) d/dr + l(l+1)/(2 r^2)

    Acting on R_q(r) = e^{-lam_q r} P_q(r) where P_q is a polynomial,
    we evaluate the action analytically:
        d/dr [e^{-lam r} P(r)] = e^{-lam r} (-lam P + P')
        d^2/dr^2 [e^{-lam r} P(r)] = e^{-lam r} (lam^2 P - 2 lam P' + P'')

    So
        -1/2 d^2/dr^2 R_q = e^{-lam_q r} (-1/2 lam_q^2 P_q + lam_q P_q' - 1/2 P_q'')
        -(1/r) d/dr R_q = e^{-lam_q r} (1/r)(lam_q P_q - P_q')

    Combining:
        T_l R_q = e^{-lam_q r} [
            (-1/2 lam_q^2) P_q
            + lam_q P_q'
            - 1/2 P_q''
            + (lam_q P_q - P_q')/r
            + (l(l+1)/2/r^2) P_q
        ]

    Then <R_p | T | R_q> = int R_p (T R_q) r^2 dr, evaluated as a sum of
    finite-power-times-exp integrals.

    Parameters
    ----------
    Same as overlap_radial.

    Returns
    -------
    float
        Kinetic matrix element.
    """
    if l_p != l_q:
        return 0.0
    l = l_p

    c_p, alpha_p = _radial_polys(n_p, l_p, lam_p)
    c_q, alpha_q = _radial_polys(n_q, l_q, lam_q)
    alpha_total = alpha_p + alpha_q

    # Build polynomial derivatives of P_q
    # P_q(r) = sum c_q[k] r^k
    # P_q'(r) = sum k c_q[k] r^{k-1}
    # P_q''(r) = sum k(k-1) c_q[k] r^{k-2}
    P_q = list(c_q) + [0.0]  # padded
    P_q_p = [0.0] * len(P_q)  # P_q'
    P_q_pp = [0.0] * len(P_q)  # P_q''
    for k in range(len(c_q)):
        if k >= 1:
            P_q_p[k - 1] += k * c_q[k]
        if k >= 2:
            P_q_pp[k - 2] += k * (k - 1) * c_q[k]

    # The action T R_q = e^{-lam_q r} * Q(r) where Q is rational in r^{-2}, r^{-1}, r^0, ...
    # We organize Q by power of r and integrate against R_p with r^2 Jacobian.
    #
    # Each piece contributes to <R_p | T | R_q> = int [sum c_p[i] r^i e^{-lam_p r}]
    #                                              * [piece] * r^2 dr
    # = sum_i c_p[i] * int r^{i+power+2} e^{-alpha_total r} dr
    # = sum_i c_p[i] * (i+power+2)! / alpha_total^{i+power+3}
    # provided power + 2 >= 0; cancellation with the 1/r and 1/r^2 terms is exact
    # for hydrogenic R because the centrifugal singularity matches the polynomial
    # leading behavior r^l.

    def _integrate_against_Rp(coeffs_q_local: List[float], extra_pow: int) -> float:
        """Integrate sum_k coeffs_q_local[k] r^{k + extra_pow} * R_p(r) * r^2 dr.

        Returns sum_{i,k} c_p[i] coeffs_q_local[k] (i+k+extra_pow+2)! / alpha_total^{i+k+extra_pow+3}.
        Requires i + k + extra_pow + 2 >= 0 for all contributing terms.
        """
        total = 0.0
        for i, ci in enumerate(c_p):
            if ci == 0.0:
                continue
            for k, ck in enumerate(coeffs_q_local):
                if ck == 0.0:
                    continue
                power = i + k + extra_pow + 2
                if power < 0:
                    # r^l leading behavior in P should ensure this cancels
                    # for hydrogenic functions; if it triggers, we have a bug.
                    raise ValueError(
                        f"kinetic_radial: divergent integrand at l={l}, "
                        f"(i={i}, k={k}, extra_pow={extra_pow})"
                    )
                total += ci * ck * factorial(power) / alpha_total ** (power + 1)
        return total

    # Term 1: -1/2 lam_q^2 * P_q   (extra_pow = 0)
    T1 = -0.5 * lam_q * lam_q * _integrate_against_Rp(list(c_q), 0)
    # Term 2: lam_q * P_q'   (extra_pow = 0)
    T2 = lam_q * _integrate_against_Rp(P_q_p, 0)
    # Term 3: -1/2 P_q''     (extra_pow = 0)
    T3 = -0.5 * _integrate_against_Rp(P_q_pp, 0)
    # Term 4: (lam_q P_q - P_q') / r   (extra_pow = -1)
    coeffs_4 = [lam_q * cq - p for cq, p in zip(c_q, P_q_p[:len(c_q)])]
    T4 = _integrate_against_Rp(coeffs_4, -1)
    # Term 5: l(l+1)/2 P_q / r^2   (extra_pow = -2)
    if l > 0:
        T5 = (l * (l + 1) / 2.0) * _integrate_against_Rp(list(c_q), -2)
    else:
        T5 = 0.0

    return T1 + T2 + T3 + T4 + T5


# ---------------------------------------------------------------------------
# Cross-exponent two-electron Slater integral R^k
# ---------------------------------------------------------------------------

def slater_rk_multifocal(
    n1: int, l1: int, lam1: float,
    n3: int, l3: int, lam3: float,
    n2: int, l2: int, lam2: float,
    n4: int, l4: int, lam4: float,
    k: int,
    n_quad: int = 100,
) -> float:
    """Cross-exponent two-electron Slater radial integral R^k.

    R^k(13;24) = int int R_{n1 l1}^{lam1}(r1) R_{n3 l3}^{lam3}(r1)
                          (r_<^k / r_>^{k+1})
                          R_{n2 l2}^{lam2}(r2) R_{n4 l4}^{lam4}(r2)
                          r1^2 r2^2 dr1 dr2

    The pair densities have decay rates lam1+lam3 (r1) and lam2+lam4 (r2),
    which can differ. The split-region representation:

        R^k = int_0^inf rho_13(r1) F_24(r1) dr1
        F_24(r1) = (1/r1^{k+1}) int_0^{r1} rho_24(r2) r2^{k+2} e^{-(lam2+lam4)r2} dr2
                  + r1^k int_{r1}^inf rho_24(r2) r2^{1-k} e^{-(lam2+lam4)r2} dr2

    Inner r2 integral closed-form via incomplete gamma at fixed r1
    (`_split_integral_analytical`); outer r1 by Gauss-Laguerre
    against weight e^{-(lam1+lam3)r1}.

    At all four lam equal, this reduces (modulo GL quadrature noise) to
    the standard ``compute_rk_float(n1,l1,n3,l3,n2,l2,n4,l4,k)`` of
    `geovac/hypergeometric_slater.py`. Quadrature noise typically <1e-9.

    Parameters
    ----------
    n1, l1, lam1, n3, l3, lam3 : the (1,3) electron-1 pair.
    n2, l2, lam2, n4, l4, lam4 : the (2,4) electron-2 pair.
    k : multipole order (>= 0).
    n_quad : Gauss-Laguerre quadrature order for the outer r1 integral.

    Returns
    -------
    float
        Slater R^k value at unit nuclear charge convention.
    """
    from scipy.special import roots_genlaguerre

    if k < 0:
        raise ValueError(f"k must be >= 0, got {k}")

    # Build product polynomials for each pair density (without r^2 weight;
    # that's added explicitly below to match the cross_register_vne pattern).
    c1, a1 = _radial_polys(n1, l1, lam1)
    c3, a3 = _radial_polys(n3, l3, lam3)
    coeffs_13 = _poly_product(np.array(c1), np.array(c3))
    alpha_13 = a1 + a3

    c2, a2 = _radial_polys(n2, l2, lam2)
    c4, a4 = _radial_polys(n4, l4, lam4)
    coeffs_24 = _poly_product(np.array(c2), np.array(c4))
    alpha_24 = a2 + a4

    # Multipole termination check: triangle inequality for k vs (l1, l3)
    # and (l2, l4) is already handled by the angular factor c_k externally.
    # Here we just compute the radial integral.

    # Gauss-Laguerre nodes/weights for the outer r1 integral against
    # weight e^{-t}, then substitute t = alpha_13 r1.
    nodes, weights = roots_genlaguerre(n_quad, 0.0)

    integral = 0.0
    for idx in range(n_quad):
        t = nodes[idx]
        w = weights[idx]
        r1 = t / alpha_13
        # Outer factor: rho_13(r1) * r1^2 = (sum c13[i] r1^i) e^{-alpha_13 r1} * r1^2
        # The exp is absorbed in the GL weight; the polynomial part:
        poly_13_with_r2 = 0.0
        for i in range(len(coeffs_13)):
            ci = coeffs_13[i]
            if ci != 0.0:
                poly_13_with_r2 += ci * r1 ** (i + 2)
        # Inner integral over r2 at fixed r1 = R_AB:
        #   F_24(r1) = (1/r1^{k+1}) int_0^{r1} rho_24(r2) r2^{k+2} e^{-alpha_24 r2} dr2
        #            + r1^k       int_{r1}^inf rho_24(r2) r2^{1-k} e^{-alpha_24 r2} dr2
        # _split_integral_analytical signature:
        #   (coeffs, alpha, p_inner, p_outer, L, R_AB)
        # with the integral
        #   I_inner * (1/R^{L+1}) + I_outer * R^L
        # where I_inner = int_0^R coeffs[k] r^{k + p_inner} e^{-alpha r} dr
        #       I_outer = int_R^inf coeffs[k] r^{k + p_outer} e^{-alpha r} dr
        # We need p_inner = k+2, p_outer = 1-k, L = k.
        F_24_r1 = _split_integral_analytical(
            coeffs_24, alpha_24, k + 2, 1 - k, k, r1
        )
        integral += w * poly_13_with_r2 * F_24_r1 / alpha_13

    return integral


def two_electron_integral_multifocal(
    p_a: MultifocalOrbital, p_b: MultifocalOrbital,
    p_c: MultifocalOrbital, p_d: MultifocalOrbital,
    n_quad: int = 100,
) -> float:
    """Two-electron physics-convention integral <ab|cd> on multi-focal basis.

    <ab|cd> = int int phi_a*(r1) phi_b*(r2) (1/r12) phi_c(r1) phi_d(r2) dr1 dr2

    With multi-focal orbitals, this decomposes via multipole expansion:

        <ab|cd> = sum_k c_k(a, c) c_k(b, d) R^k(ac; bd)

    where the angular factors c_k are unchanged from single-focal
    (Gaunt selection rules depend on (l, m) only, not on lambda).

    Selection rules (Paper 22 angular sparsity preserved):
        |l_a - l_c| <= k <= l_a + l_c
        |l_b - l_d| <= k <= l_b + l_d
        l_a + l_c + k even
        l_b + l_d + k even
        m_a - m_c = m_d - m_b  (cross-pair m-conservation)

    Parameters
    ----------
    p_a, p_b, p_c, p_d : MultifocalOrbital
        Four orbitals. (a, c) is the r_1 pair; (b, d) is the r_2 pair.
    n_quad : int
        Gauss-Laguerre order for the outer radial integral.

    Returns
    -------
    float
        Two-electron integral value.
    """
    # m-conservation
    if p_a.m + p_b.m != p_c.m + p_d.m:
        return 0.0

    k_min = max(abs(p_a.l - p_c.l), abs(p_b.l - p_d.l))
    k_max = min(p_a.l + p_c.l, p_b.l + p_d.l)

    result = 0.0
    for k in range(k_min, k_max + 1):
        if (p_a.l + p_c.l + k) % 2 != 0:
            continue
        if (p_b.l + p_d.l + k) % 2 != 0:
            continue

        ck_ac = _gaunt_ck(p_a.l, p_a.m, p_c.l, p_c.m, k)
        ck_bd = _gaunt_ck(p_b.l, p_b.m, p_d.l, p_d.m, k)
        if abs(ck_ac) < 1e-15 or abs(ck_bd) < 1e-15:
            continue

        rk = slater_rk_multifocal(
            p_a.n, p_a.l, p_a.lam,
            p_c.n, p_c.l, p_c.lam,
            p_b.n, p_b.l, p_b.lam,
            p_d.n, p_d.l, p_d.lam,
            k, n_quad=n_quad,
        )
        result += ck_ac * ck_bd * rk

    return result


# ---------------------------------------------------------------------------
# One-body matrices
# ---------------------------------------------------------------------------

def overlap_matrix(spec: MultifocalSpec) -> np.ndarray:
    """Build the (n_orb, n_orb) overlap matrix S_pq = <p|q>."""
    n = spec.n_orbitals
    S = np.zeros((n, n))
    for i, p in enumerate(spec.orbitals):
        for j, q in enumerate(spec.orbitals):
            if p.l != q.l or p.m != q.m:
                continue
            S[i, j] = overlap_radial(p.n, p.l, p.lam, q.n, q.l, q.lam)
    return S


def h1_matrix(spec: MultifocalSpec) -> np.ndarray:
    """Build the one-body Hamiltonian h1 = T - Z/r matrix.

    Diagonal in (l, m) — kinetic and Coulomb don't change angular character.

    Returns
    -------
    np.ndarray
        (n_orb, n_orb) symmetric matrix.
    """
    n = spec.n_orbitals
    h1 = np.zeros((n, n))
    Z = spec.Z_nuc
    for i, p in enumerate(spec.orbitals):
        for j, q in enumerate(spec.orbitals):
            if p.l != q.l or p.m != q.m:
                continue
            T = kinetic_radial(p.n, p.l, p.lam, q.n, q.l, q.lam)
            V = -Z * matrix_element_rk(p.n, p.l, p.lam, q.n, q.l, q.lam, k_pow=-1)
            h1[i, j] = T + V
    # Symmetrize numerically (tiny non-symmetry from independent T and V evals)
    h1 = 0.5 * (h1 + h1.T)
    return h1


def dipole_z_matrix(spec: MultifocalSpec) -> np.ndarray:
    """Build the (n_orb, n_orb) dipole z matrix <p|z|q>.

    z = r cos(theta) couples l = l_p +/- 1, m = m_p (no Delta_m for z).

    Returns
    -------
    np.ndarray
        (n_orb, n_orb) real matrix (anti-Hermitian-by-convention but real
        for real spherical harmonics; symmetric here).
    """
    n = spec.n_orbitals
    D = np.zeros((n, n))
    for i, p in enumerate(spec.orbitals):
        for j, q in enumerate(spec.orbitals):
            # Selection: l differs by 1, m equal
            if abs(p.l - q.l) != 1:
                continue
            if p.m != q.m:
                continue
            c1 = _gaunt_ck(p.l, p.m, q.l, q.m, 1)
            if abs(c1) < 1e-15:
                continue
            R = matrix_element_rk(p.n, p.l, p.lam, q.n, q.l, q.lam, k_pow=1)
            D[i, j] = c1 * R
    return D


# ---------------------------------------------------------------------------
# CI on multi-focal basis (singlet, two-electron)
# ---------------------------------------------------------------------------

def build_singlet_LM_subblock_multifocal(
    spec: MultifocalSpec,
    L_target: int,
    M_L_target: int,
    n_quad: int = 100,
) -> Tuple[np.ndarray, np.ndarray, List[Tuple[int, int]]]:
    """Build the multi-focal singlet two-electron CI matrices (H, S).

    Configurations are spatial pairs (i, j) where orbital i has l_i = 0
    (one electron in s) and orbital j has l_j = L_target (the other in
    the L_target shell), with m_i + m_j = M_L_target. For L_target = 0
    we take i <= j (avoid double-counting).

    The basis is non-orthogonal — distinct lambdas at the same (n, l, m)
    give nonzero overlap. The CI generalized eigenproblem is

        H c = E S c

    where S is the configuration-space overlap matrix (NOT the orbital
    overlap matrix; configurations are antisymmetrized products and have
    their own overlap structure).

    For singlet two-electron configurations |I> = (|i j> + |j i>) / N_IJ:

        S_{IJ} = (1/(N_I N_J)) sum_{(a,b)<-I, (c,d)<-J} sign * S_orb[a,c] S_orb[b,d]

    where sign accounts for the singlet (symmetric) symmetry and N_IJ = sqrt(2)
    for i != j, N_II = 1.

    Similarly for H = h1 + h1 + g (one-body diagonal + two-body Coulomb).

    Parameters
    ----------
    spec : MultifocalSpec
        Multi-focal orbital basis specification.
    L_target : int
        Target L for the second electron (0 for ss singlet, 1 for sp 1S/1P).
    M_L_target : int
        Total M_L. For ss: M_L = 0. For sp: M_L from 0 (m=0 sublevel).
    n_quad : int
        Gauss-Laguerre quadrature order for two-electron Slater integrals.

    Returns
    -------
    H : (n_configs, n_configs) Hamiltonian matrix
    S : (n_configs, n_configs) overlap matrix
    configs : list of (i, j) orbital index pairs
    """
    orbs = spec.orbitals

    # Identify s-orbitals (l=0, m=0)
    s_indices = [i for i, p in enumerate(orbs) if p.l == 0 and p.m == 0]
    target_indices = [
        i for i, p in enumerate(orbs)
        if p.l == L_target and p.m == M_L_target
    ]

    if L_target == 0:
        # ss singlet: i <= j (over s-orbitals only)
        configs: List[Tuple[int, int]] = [
            (i, j) for i in s_indices for j in s_indices if i <= j
        ]
    else:
        # sp singlet (L_target != 0): one s-orbital, one target-l-orbital
        configs = [(i, j) for i in s_indices for j in target_indices]

    n_configs = len(configs)
    if n_configs == 0:
        return np.zeros((0, 0)), np.zeros((0, 0)), []

    # Pre-compute one-body matrix and overlap
    S_orb = overlap_matrix(spec)
    h1_orb = h1_matrix(spec)

    H = np.zeros((n_configs, n_configs))
    S_ci = np.zeros((n_configs, n_configs))

    parity = 1.0  # singlet (symmetric spatial)

    # Cache for two-electron integrals
    g_cache: Dict[Tuple[int, int, int, int], float] = {}

    def g_int(a: int, b: int, c: int, d: int) -> float:
        key = (a, b, c, d)
        if key not in g_cache:
            g_cache[key] = two_electron_integral_multifocal(
                orbs[a], orbs[b], orbs[c], orbs[d], n_quad=n_quad,
            )
        return g_cache[key]

    for I in range(n_configs):
        i, j = configs[I]
        bra_perms = [(i, j, 1.0)]
        if i != j:
            bra_perms.append((j, i, parity))
        N_I = sqrt(float(len(bra_perms)))

        for J in range(I, n_configs):
            p, q = configs[J]
            ket_perms = [(p, q, 1.0)]
            if p != q:
                ket_perms.append((q, p, parity))
            N_J = sqrt(float(len(ket_perms)))

            ham_me = 0.0
            ovlp_me = 0.0
            for a, b, s_bra in bra_perms:
                for c, d, s_ket in ket_perms:
                    sign = s_bra * s_ket
                    # Overlap S = <a|c> <b|d>
                    ovlp_me += sign * S_orb[a, c] * S_orb[b, d]
                    # h1 contribution: h_ac <b|d> + <a|c> h_bd
                    ham_me += sign * h1_orb[a, c] * S_orb[b, d]
                    ham_me += sign * S_orb[a, c] * h1_orb[b, d]
                    # Two-electron g_ab,cd in physics convention <ab|cd>
                    ham_me += sign * g_int(a, b, c, d)

            ham_me /= (N_I * N_J)
            ovlp_me /= (N_I * N_J)
            H[I, J] = ham_me
            H[J, I] = ham_me
            S_ci[I, J] = ovlp_me
            S_ci[J, I] = ovlp_me

    return H, S_ci, configs


# ---------------------------------------------------------------------------
# Phase D extension: angular-CI extended subblocks
# ---------------------------------------------------------------------------
#
# The Phase B subblock builder restricts the L=0 (1S) sector to (s, s)
# configurations and the L=1 (1P) sector to (s, p) configurations. For
# He 2^1P -> 1^1S oscillator strength this misses TWO physically real
# contributions:
#
#   (a) 2p^2 1S configurations enter the 1^1S CI ground state (Drake
#       reports a ~1% admixture). They are NOT (s, s) and so are missed.
#       For multi-focal at finite n_max this admixture can be larger
#       than 1% because the basis is more flexible.
#
#   (b) (p, d) 1P configurations enter the 2^1P CI excited state with
#       similar admixture-level. They couple via the 2-electron Coulomb
#       to (s, p) dominantly.
#
# Both extensions are mechanical: build the configuration list as the
# union over all (l_a, l_b) compatible with the target L and singlet
# spatial parity (l_a + l_b + L even), then enumerate (orbital_a,
# orbital_b) at each m-distribution m_a + m_b = M_L. The Hamiltonian
# matrix elements have the same Slater-Condon form (the 2-electron
# Coulomb integrals already preserve all angular selection via Gaunt).
#
# We expose this through a NEW function that defaults to the same
# behaviour as Phase B when called with the default config_l_pairs:
#
#   - L=0: [(0, 0)] reproduces Phase B exactly.
#   - L=1: [(0, 1)] reproduces Phase B exactly.
#
# Calling with config_l_pairs=None enumerates ALL singlet-allowed pairs
# up to the maximum l in the spec. This is the Phase D production path.


def _enumerate_singlet_l_pairs(
    spec: MultifocalSpec, L_target: int,
) -> List[Tuple[int, int]]:
    """Enumerate (l_a, l_b) angular pairs for singlet two-electron CI at L_target.

    Singlet selection: spatial wavefunction is symmetric, so l_a + l_b + L
    must be EVEN (Pauli + spin singlet). Triangle inequality:
    |l_a - l_b| <= L_target <= l_a + l_b.

    For L=0 (1S): pairs (0,0), (1,1), (2,2), ... with l_a + l_b even.
    For L=1 (1P): pairs (0,1), (1,2), (2,3), ... with l_a + l_b + 1 even,
    i.e. l_a + l_b odd. So (0,1), (1,2), (2,3), ...
    For L=2 (1D): pairs (0,2), (1,1), (1,3), (2,2), ... with l_a + l_b even.

    We restrict to l_a <= l_b to avoid double-counting (the symmetric
    case is the singlet spatial wavefunction).

    Parameters
    ----------
    spec : MultifocalSpec
    L_target : int

    Returns
    -------
    list of (l_a, l_b) tuples with l_a <= l_b.
    """
    l_max = max(o.l for o in spec.orbitals)
    pairs: List[Tuple[int, int]] = []
    for la in range(l_max + 1):
        for lb in range(la, l_max + 1):  # la <= lb
            # Triangle inequality
            if abs(la - lb) > L_target or L_target > la + lb:
                continue
            # Singlet spatial parity: la + lb + L_target must be even
            if (la + lb + L_target) % 2 != 0:
                continue
            pairs.append((la, lb))
    return pairs


def _enumerate_m_distributions(
    la: int, lb: int, M_L: int,
) -> List[Tuple[int, int]]:
    """All (m_a, m_b) with |m_a| <= la, |m_b| <= lb, m_a + m_b = M_L."""
    out = []
    for ma in range(-la, la + 1):
        mb = M_L - ma
        if abs(mb) <= lb:
            out.append((ma, mb))
    return out


def build_singlet_LM_subblock_multifocal_extended(
    spec: MultifocalSpec,
    L_target: int,
    M_L_target: int,
    n_quad: int = 100,
    config_l_pairs: Optional[List[Tuple[int, int]]] = None,
) -> Tuple[np.ndarray, np.ndarray, List[Tuple[int, int]]]:
    """Build the multi-focal singlet two-electron CI matrices (H, S) — extended.

    Phase D enhancement of ``build_singlet_LM_subblock_multifocal``: includes
    ALL singlet-allowed angular pair channels (l_a, l_b), not just (0, 0)
    or (0, L_target).

    Configurations are spatial pairs (orb_i, orb_j) where orbital i has
    angular character l_a and orbital j has l_b, with l_a <= l_b drawn
    from ``config_l_pairs``, and m_a + m_b = M_L_target.

    For l_a == l_b: take orb_i <= orb_j (avoid double-counting); also
    require that the orbitals contain l_a, m_a sublevels.

    Parameters
    ----------
    spec : MultifocalSpec
    L_target : int
        Target total L for the CI sub-block.
    M_L_target : int
        Total M_L.
    n_quad : int
        GL order.
    config_l_pairs : list of (l_a, l_b), optional
        If provided, ONLY include these angular pairs. If None, include
        all singlet-allowed pairs up to the max l in the spec.

    Returns
    -------
    H, S_ci, configs : as ``build_singlet_LM_subblock_multifocal``.
    """
    orbs = spec.orbitals

    if config_l_pairs is None:
        config_l_pairs = _enumerate_singlet_l_pairs(spec, L_target)

    # Build configuration list across all (l_a, l_b) pairs and m-distributions
    configs: List[Tuple[int, int]] = []
    for la, lb in config_l_pairs:
        # Triangle / parity validation (defensive)
        if abs(la - lb) > L_target or L_target > la + lb:
            continue
        if (la + lb + L_target) % 2 != 0:
            continue

        m_dists = _enumerate_m_distributions(la, lb, M_L_target)

        for ma, mb in m_dists:
            a_indices = [i for i, p in enumerate(orbs)
                         if p.l == la and p.m == ma]
            b_indices = [i for i, p in enumerate(orbs)
                         if p.l == lb and p.m == mb]

            if la == lb and ma == mb:
                # Same (l, m) sublevel: use i <= j
                for i in a_indices:
                    for j in b_indices:
                        if i <= j:
                            configs.append((i, j))
            elif la == lb and ma != mb:
                # Same l but different m: i in (la, ma), j in (lb, mb).
                # To avoid double-counting we'd need to coordinate with
                # the (lb, mb) - (la, ma) pair which has the same (la, lb).
                # Use ordering: only include if (ma, i) < (mb, j) lex.
                for i in a_indices:
                    for j in b_indices:
                        # i in (la, ma), j in (la, mb), ma != mb.
                        # The complementary config would be i' in (la, mb)
                        # and j' in (la, ma); with our enumeration we get
                        # both (la, lb)=(la, la) entries, including (ma, mb)
                        # AND (mb, ma). To avoid duplicate, take only the
                        # case ma < mb (canonical ordering).
                        if ma < mb:
                            configs.append((i, j))
            else:
                # Different l: la < lb strictly. All combinations are
                # distinct (l-block disjoint).
                for i in a_indices:
                    for j in b_indices:
                        configs.append((i, j))

    n_configs = len(configs)
    if n_configs == 0:
        return np.zeros((0, 0)), np.zeros((0, 0)), []

    # Pre-compute one-body matrix and overlap (full, all l)
    S_orb = overlap_matrix(spec)
    h1_orb = h1_matrix(spec)

    H = np.zeros((n_configs, n_configs))
    S_ci = np.zeros((n_configs, n_configs))

    parity = 1.0  # singlet (symmetric spatial)

    # Cache for two-electron integrals
    g_cache: Dict[Tuple[int, int, int, int], float] = {}

    def g_int(a: int, b: int, c: int, d: int) -> float:
        key = (a, b, c, d)
        if key not in g_cache:
            g_cache[key] = two_electron_integral_multifocal(
                orbs[a], orbs[b], orbs[c], orbs[d], n_quad=n_quad,
            )
        return g_cache[key]

    for I in range(n_configs):
        i, j = configs[I]
        bra_perms = [(i, j, 1.0)]
        if i != j:
            bra_perms.append((j, i, parity))
        N_I = sqrt(float(len(bra_perms)))

        for J in range(I, n_configs):
            p, q = configs[J]
            ket_perms = [(p, q, 1.0)]
            if p != q:
                ket_perms.append((q, p, parity))
            N_J = sqrt(float(len(ket_perms)))

            ham_me = 0.0
            ovlp_me = 0.0
            for a, b, s_bra in bra_perms:
                for c, d, s_ket in ket_perms:
                    sign = s_bra * s_ket
                    # Overlap S = <a|c> <b|d>
                    ovlp_me += sign * S_orb[a, c] * S_orb[b, d]
                    # h1 contribution: h_ac <b|d> + <a|c> h_bd
                    ham_me += sign * h1_orb[a, c] * S_orb[b, d]
                    ham_me += sign * S_orb[a, c] * h1_orb[b, d]
                    # Two-electron g_ab,cd in physics convention <ab|cd>
                    ham_me += sign * g_int(a, b, c, d)

            ham_me /= (N_I * N_J)
            ovlp_me /= (N_I * N_J)
            H[I, J] = ham_me
            H[J, I] = ham_me
            S_ci[I, J] = ovlp_me
            S_ci[J, I] = ovlp_me

    return H, S_ci, configs


# ---------------------------------------------------------------------------
# CI eigensolve (generalized symmetric eigenproblem)
# ---------------------------------------------------------------------------

def solve_generalized_singlet(
    H: np.ndarray, S: np.ndarray,
    cond_threshold: float = 1e10,
) -> Tuple[np.ndarray, np.ndarray]:
    """Solve H c = E S c for the multi-focal singlet CI.

    Uses scipy.linalg.eigh in generalized mode. Logs a warning if
    the overlap matrix is ill-conditioned (cond > cond_threshold).

    Parameters
    ----------
    H : (n, n) symmetric Hamiltonian matrix.
    S : (n, n) symmetric positive-definite overlap matrix.
    cond_threshold : float
        Warn (raise RuntimeWarning) if cond(S) > threshold.

    Returns
    -------
    eigvals : (n,) sorted eigenvalues.
    eigvecs : (n, n) eigenvectors (each column normalized to c^T S c = 1).
    """
    from scipy.linalg import eigh

    if H.shape != S.shape or H.shape[0] != H.shape[1]:
        raise ValueError(f"Shape mismatch: H={H.shape}, S={S.shape}")

    if H.shape[0] == 0:
        return np.zeros(0), np.zeros((0, 0))

    cond_S = np.linalg.cond(S)
    if cond_S > cond_threshold:
        import warnings
        warnings.warn(
            f"Multi-focal CI: overlap matrix is ill-conditioned, "
            f"cond(S) = {cond_S:.3e}. Results may be unreliable.",
            RuntimeWarning,
        )

    # eigh with b=S solves the generalized symmetric eigenproblem
    eigvals, eigvecs = eigh(H, b=S)
    return eigvals, eigvecs


# ---------------------------------------------------------------------------
# Transition dipole on the multi-focal basis
# ---------------------------------------------------------------------------

def transition_dipole_multifocal(
    psi_init: np.ndarray, configs_init: List[Tuple[int, int]],
    psi_final: np.ndarray, configs_final: List[Tuple[int, int]],
    spec_init: MultifocalSpec,
    spec_final: MultifocalSpec,
    parity_init: float = 1.0,
    parity_final: float = 1.0,
) -> float:
    """Two-electron singlet-singlet transition dipole <Psi_init | z_1 + z_2 | Psi_final>.

    Each Psi is a CI vector in the (possibly different) multi-focal
    spec_init / spec_final basis. The cross-spec dipole matrix element
    <orb_init | z | orb_final> uses cross-exponent dipole.

    For the He 2^1P -> 1^1S case we use a single common spec (same orbital
    basis for both states), so spec_init == spec_final. The function
    supports distinct specs for generality (e.g. if final-state-optimized
    orbitals differ from initial-state-optimized).

    Parameters
    ----------
    psi_init, psi_final : np.ndarray
        CI eigenvectors (in their own configuration basis).
    configs_init, configs_final : list of (int, int)
        Configuration index pairs.
    spec_init, spec_final : MultifocalSpec
        Multi-focal orbital specs for the initial and final states.
        Must have identical orbital ordering (same indices refer to the
        same orbital). For Phase B we require spec_init is spec_final.
    parity_init, parity_final : float
        Spatial-symmetry parity (singlet = +1, triplet = -1).

    Returns
    -------
    float
        Transition dipole <Psi_init | z | Psi_final> (length form).
    """
    if spec_init is not spec_final:
        # Phase B: only support shared spec. Distinct specs would require
        # cross-spec one-body dipole (cross-orbital-set <orb_a|z|orb_c>),
        # which is mechanical but not in scope.
        raise NotImplementedError(
            "Phase B: spec_init must be the same object as spec_final; "
            "distinct specs not yet supported (would require cross-spec "
            "dipole machinery)."
        )

    spec = spec_init
    n_orb = spec.n_orbitals

    # Pre-compute single-orbital dipole matrix
    D_orb = dipole_z_matrix(spec)

    # Slater-Condon assembly for two-body singlet-to-singlet on a one-body operator
    matrix_element = 0.0
    for I, c_I in enumerate(psi_init):
        if abs(c_I) < 1e-15:
            continue
        i, j = configs_init[I]
        bra_perms = [(i, j, 1.0)]
        if i != j:
            bra_perms.append((j, i, parity_init))
        N_I = sqrt(float(len(bra_perms)))

        for J, c_J in enumerate(psi_final):
            if abs(c_J) < 1e-15:
                continue
            p, q = configs_final[J]
            ket_perms = [(p, q, 1.0)]
            if p != q:
                ket_perms.append((q, p, parity_final))
            N_J = sqrt(float(len(ket_perms)))

            me_IJ = 0.0
            for a, b, s_bra in bra_perms:
                for c, d, s_ket in ket_perms:
                    sign = s_bra * s_ket
                    # <a|z|c> <b|d> + <a|c> <b|z|d>
                    # In a non-orthogonal basis we must include the
                    # spectator overlap, NOT delta:
                    me_IJ += sign * D_orb[a, c] * overlap_radial(
                        spec.orbitals[b].n, spec.orbitals[b].l, spec.orbitals[b].lam,
                        spec.orbitals[d].n, spec.orbitals[d].l, spec.orbitals[d].lam,
                    ) * (1.0 if (spec.orbitals[b].l == spec.orbitals[d].l and
                                  spec.orbitals[b].m == spec.orbitals[d].m) else 0.0)
                    me_IJ += sign * D_orb[b, d] * overlap_radial(
                        spec.orbitals[a].n, spec.orbitals[a].l, spec.orbitals[a].lam,
                        spec.orbitals[c].n, spec.orbitals[c].l, spec.orbitals[c].lam,
                    ) * (1.0 if (spec.orbitals[a].l == spec.orbitals[c].l and
                                  spec.orbitals[a].m == spec.orbitals[c].m) else 0.0)

            me_IJ /= (N_I * N_J)
            matrix_element += c_I * c_J * me_IJ

    return matrix_element


# ---------------------------------------------------------------------------
# Driver: He 2^1P -> 1^1S oscillator strength on multi-focal basis
# ---------------------------------------------------------------------------

def compute_he_oscillator_strength_multifocal(
    spec: MultifocalSpec,
    n_quad: int = 100,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Compute He 2^1P -> 1^1S oscillator strength on a multi-focal basis.

    f = 2 * omega * |<1^1S | z | 2^1P_{m=0}>|^2

    The factor 2 (rather than 2/3) accounts for summing over the upper-
    state P sublevels (m = -1, 0, +1) which by rotational invariance all
    contribute equally; the m=0 component is computed and multiplied by 3,
    canceling the 1/3 in the (2/3) prefactor.

    Parameters
    ----------
    spec : MultifocalSpec
        The multi-focal basis. Must include enough s-orbitals for the
        1^1S CI (ideally with at least the variational 1s lambda),
        enough p-orbitals for the 2^1P CI (with appropriate screening
        lambda), and m=0 sublevels of p-orbitals.
    n_quad : int
        Gauss-Laguerre order for two-electron Slater integrals.
    verbose : bool
        Print diagnostic info.

    Returns
    -------
    Dict containing:
        f_length, omega_Ha, dipole_z_au, E_1S_Ha, E_2P_Ha,
        n_configs_1S, n_configs_2P, cond_S_1S, cond_S_2P.
    """
    if verbose:
        print(f"  Multi-focal spec: {spec.n_orbitals} orbitals, Z={spec.Z_nuc}")
        for o in spec.orbitals:
            print(f"    ({o.n}, {o.l}, m={o.m}, lam={o.lam:.4f})  [{o.label}]")

    # Build 1^1S sub-block (l_1 = l_2 = 0, M_L = 0)
    H_1S, S_1S, configs_1S = build_singlet_LM_subblock_multifocal(
        spec, L_target=0, M_L_target=0, n_quad=n_quad,
    )
    if verbose:
        print(f"  1S configs: {len(configs_1S)}")

    eigvals_1S, eigvecs_1S = solve_generalized_singlet(H_1S, S_1S)
    psi_1S = eigvecs_1S[:, 0]
    E_1S = eigvals_1S[0]

    # Build 2^1P sub-block (l_1 = 0, l_2 = 1, M_L = 0 -> m_2 = 0)
    H_2P, S_2P, configs_2P = build_singlet_LM_subblock_multifocal(
        spec, L_target=1, M_L_target=0, n_quad=n_quad,
    )
    if verbose:
        print(f"  2P configs: {len(configs_2P)}")

    eigvals_2P, eigvecs_2P = solve_generalized_singlet(H_2P, S_2P)
    psi_2P = eigvecs_2P[:, 0]
    E_2P = eigvals_2P[0]

    omega = E_2P - E_1S

    # Transition dipole z-component
    z_me = transition_dipole_multifocal(
        psi_1S, configs_1S, psi_2P, configs_2P,
        spec_init=spec, spec_final=spec,
    )

    f_length = 2.0 * omega * z_me ** 2

    cond_1S = float(np.linalg.cond(S_1S)) if S_1S.size > 0 else 1.0
    cond_2P = float(np.linalg.cond(S_2P)) if S_2P.size > 0 else 1.0

    return {
        "f_length": float(f_length),
        "omega_Ha": float(omega),
        "dipole_z_au": float(z_me),
        "dipole_squared_au": float(z_me ** 2),
        "E_1S_Ha": float(E_1S),
        "E_2P_Ha": float(E_2P),
        "n_configs_1S": len(configs_1S),
        "n_configs_2P": len(configs_2P),
        "n_orbitals": spec.n_orbitals,
        "Z_nuc": float(spec.Z_nuc),
        "cond_S_1S": cond_1S,
        "cond_S_2P": cond_2P,
        "n_quad": n_quad,
    }


# ---------------------------------------------------------------------------
# Convenience constructors for He
# ---------------------------------------------------------------------------

def compute_he_oscillator_strength_multifocal_extended(
    spec: MultifocalSpec,
    n_quad: int = 100,
    config_l_pairs_1S: Optional[List[Tuple[int, int]]] = None,
    config_l_pairs_2P: Optional[List[Tuple[int, int]]] = None,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Phase D extended driver: He 2^1P -> 1^1S oscillator strength on
    multi-focal basis with extended angular configuration set.

    Differs from ``compute_he_oscillator_strength_multifocal`` in that
    it uses ``build_singlet_LM_subblock_multifocal_extended``, allowing
    angular pairs (l_a, l_b) beyond just (0, 0) for 1S and (0, 1) for 1P.

    For physically meaningful gains we typically include:
      * 1S sector: [(0, 0), (1, 1), (2, 2), ...] up to spec's max l.
                   Each (l, l) pair contributes the (np, np') and similar
                   configurations summed over all m_a + m_b = 0.
      * 2P sector: [(0, 1), (1, 2), ...] singlet-allowed pairs.

    Parameters
    ----------
    spec : MultifocalSpec
        Multi-focal basis. Must include enough s, p, d sublevels for the
        target l-pairs to be populated. Typically you'd want all m
        sublevels for any l >= 1 you want to use in (l, l) configurations.
    n_quad : int
        Gauss-Laguerre order for two-electron integrals.
    config_l_pairs_1S : list of (l_a, l_b), optional
        Restrict 1S subblock to these angular pairs. Default: all
        singlet-allowed up to spec's max l.
    config_l_pairs_2P : list of (l_a, l_b), optional
        Restrict 2P subblock similarly. Default: all singlet-allowed.
    verbose : bool
        Print diagnostic info.

    Returns
    -------
    Dict with keys: f_length, omega_Ha, dipole_z_au, E_1S_Ha, E_2P_Ha,
    n_configs_1S, n_configs_2P, cond_S_1S, cond_S_2P, plus
    config_l_pairs_used_1S, config_l_pairs_used_2P (the actual pair lists).
    """
    if verbose:
        print(f"  Multi-focal spec (extended): {spec.n_orbitals} orbitals, "
              f"Z={spec.Z_nuc}")
        for o in spec.orbitals:
            print(f"    ({o.n}, {o.l}, m={o.m}, lam={o.lam:.4f})  [{o.label}]")

    # Build 1^1S sub-block (extended angular CI)
    H_1S, S_1S, configs_1S = build_singlet_LM_subblock_multifocal_extended(
        spec, L_target=0, M_L_target=0, n_quad=n_quad,
        config_l_pairs=config_l_pairs_1S,
    )
    if verbose:
        print(f"  1S configs: {len(configs_1S)}")

    eigvals_1S, eigvecs_1S = solve_generalized_singlet(H_1S, S_1S)
    psi_1S = eigvecs_1S[:, 0]
    E_1S = eigvals_1S[0]

    # Build 2^1P sub-block (extended angular CI)
    H_2P, S_2P, configs_2P = build_singlet_LM_subblock_multifocal_extended(
        spec, L_target=1, M_L_target=0, n_quad=n_quad,
        config_l_pairs=config_l_pairs_2P,
    )
    if verbose:
        print(f"  2P configs: {len(configs_2P)}")

    eigvals_2P, eigvecs_2P = solve_generalized_singlet(H_2P, S_2P)
    psi_2P = eigvecs_2P[:, 0]
    E_2P = eigvals_2P[0]

    omega = E_2P - E_1S

    z_me = transition_dipole_multifocal(
        psi_1S, configs_1S, psi_2P, configs_2P,
        spec_init=spec, spec_final=spec,
    )

    f_length = 2.0 * omega * z_me ** 2

    cond_1S = float(np.linalg.cond(S_1S)) if S_1S.size > 0 else 1.0
    cond_2P = float(np.linalg.cond(S_2P)) if S_2P.size > 0 else 1.0

    actual_pairs_1S = (
        _enumerate_singlet_l_pairs(spec, 0)
        if config_l_pairs_1S is None else list(config_l_pairs_1S)
    )
    actual_pairs_2P = (
        _enumerate_singlet_l_pairs(spec, 1)
        if config_l_pairs_2P is None else list(config_l_pairs_2P)
    )

    return {
        "f_length": float(f_length),
        "omega_Ha": float(omega),
        "dipole_z_au": float(z_me),
        "dipole_squared_au": float(z_me ** 2),
        "E_1S_Ha": float(E_1S),
        "E_2P_Ha": float(E_2P),
        "n_configs_1S": len(configs_1S),
        "n_configs_2P": len(configs_2P),
        "n_orbitals": spec.n_orbitals,
        "Z_nuc": float(spec.Z_nuc),
        "cond_S_1S": cond_1S,
        "cond_S_2P": cond_2P,
        "n_quad": n_quad,
        "config_l_pairs_1S": actual_pairs_1S,
        "config_l_pairs_2P": actual_pairs_2P,
    }


def he_extended_spec(
    n_max: int = 4,
    s_lams: Optional[List[float]] = None,
    p_lams: Optional[List[float]] = None,
    d_lams: Optional[List[float]] = None,
    include_p_all_m: bool = True,
    include_d_all_m: bool = True,
) -> MultifocalSpec:
    """Phase D: extended He spec with ALL m-sublevels for l>=1 orbitals.

    Necessary for the extended angular CI: (1, 1) configurations need
    p-orbitals at m = -1, 0, +1 to span the full M_L=0 sub-sector. (Drake
    2p^2 1S contains mostly the (m=+1, m=-1) and (m=0)^2 components.)

    Parameters
    ----------
    n_max : int
        Maximum principal quantum number.
    s_lams, p_lams, d_lams : list of float, optional
        Per-l focal lengths to include. If None, uses Slater rules with
        n=1 at variational 27/16, n=2 at 0.575, n=3 at 0.333.
    include_p_all_m : bool
        If True, include m = -1, 0, +1 sublevels for each p-orbital.
        Default True for Phase D extended CI.
    include_d_all_m : bool
        Similar for d-orbitals (m = -2, -1, 0, +1, +2).

    Returns
    -------
    MultifocalSpec
    """
    Z = 2.0
    orbitals: List[MultifocalOrbital] = []

    if s_lams is None:
        # Defaults: variational 1s + Slater n=2..n_max
        s_lams_default = [27.0 / 16.0]
        for n in range(2, n_max + 1):
            lam_n = (Z - 0.85) / n if n == 2 else (Z - 1.0) / n
            s_lams_default.append(lam_n)
        s_lams = s_lams_default
    if p_lams is None and n_max >= 2:
        p_lams_default = []
        for n in range(2, n_max + 1):
            lam_n = (Z - 0.85) / n if n == 2 else (Z - 1.0) / n
            p_lams_default.append(lam_n)
        p_lams = p_lams_default
    if d_lams is None and n_max >= 3:
        d_lams_default = []
        for n in range(3, n_max + 1):
            lam_n = (Z - 1.0) / n
            d_lams_default.append(lam_n)
        d_lams = d_lams_default

    # s-orbitals: assign n labels 1, 2, 3, ... in order
    if s_lams is None:
        s_lams = []
    for k, lam in enumerate(s_lams):
        n = k + 1
        orbitals.append(MultifocalOrbital(
            n=n, l=0, m=0, lam=lam, label=f"s{n}_lam{lam:.3f}",
        ))

    # p-orbitals: lowest n for p is 2; assign n=2, 3, ...
    if p_lams is None:
        p_lams = []
    for k, lam in enumerate(p_lams):
        n = k + 2  # p starts at n=2
        ms = [-1, 0, 1] if include_p_all_m else [0]
        for m in ms:
            orbitals.append(MultifocalOrbital(
                n=n, l=1, m=m, lam=lam, label=f"p{n}_m{m:+d}_lam{lam:.3f}",
            ))

    # d-orbitals: lowest n for d is 3
    if d_lams is None:
        d_lams = []
    for k, lam in enumerate(d_lams):
        n = k + 3  # d starts at n=3
        ms = [-2, -1, 0, 1, 2] if include_d_all_m else [0]
        for m in ms:
            orbitals.append(MultifocalOrbital(
                n=n, l=2, m=m, lam=lam, label=f"d{n}_m{m:+d}_lam{lam:.3f}",
            ))

    return MultifocalSpec(
        orbitals=orbitals,
        Z_nuc=Z,
        label=f"he_extended_nmax{n_max}_s{len(s_lams)}_p{len(p_lams)}_d{len(d_lams)}",
    )


def he_slater_spec(n_max: int = 2, m_p: int = 0) -> MultifocalSpec:
    """Build a Slater-rules multi-focal spec for He.

    Slater's rules for He at Z=2:
        1s: lam = Z - 0.30 = 1.70 (Slater for innermost, screening 0.30 per same-shell)
            (for 1s^2 the textbook variational is 27/16 = 1.6875)
        2s, 2p: lam = (Z - 0.85)/2 = 0.575 (n=2 valence, 1s shields by 0.85)
        3s, 3p, 3d: lam = (Z - 1.0)/3 ~= 0.333 (n=3 valence)

    For oscillator strength we use the variational 27/16 for 1s and
    Slater (Z-0.85)/n for n=2,3,...

    Parameters
    ----------
    n_max : int
        Maximum principal quantum number to include.
    m_p : int
        m sublevel for the p-orbital (typically 0 for the M_L=0 component
        of the 2^1P state).

    Returns
    -------
    MultifocalSpec
        Spec with orbitals (1s, 2s, 2p_{m_p}, 3s, 3p_{m_p}, 3d_{m_p}, ...)
        out to n_max.
    """
    Z = 2.0
    lam_1s = 27.0 / 16.0  # variational
    orbitals = [MultifocalOrbital(n=1, l=0, m=0, lam=lam_1s, label="1s_var")]

    for n in range(2, n_max + 1):
        lam_n = (Z - 0.85) / n if n == 2 else (Z - 1.0) / n
        for l in range(n):
            # Only include the m=0 sublevel for p-orbitals (or m=m_p for spec)
            # For l=0, m=0; for l >= 1 just include m=m_p in this driver
            # to keep the basis minimal
            if l == 0:
                orbitals.append(MultifocalOrbital(
                    n=n, l=l, m=0, lam=lam_n, label=f"{n}s_screen",
                ))
            else:
                orbitals.append(MultifocalOrbital(
                    n=n, l=l, m=m_p, lam=lam_n, label=f"{n}{'spdf'[l]}_m={m_p}",
                ))

    return MultifocalSpec(orbitals=orbitals, Z_nuc=Z, label=f"he_slater_nmax{n_max}")


def he_uniform_spec(n_max: int, lam: float, m_p: int = 0) -> MultifocalSpec:
    """Build a uniform-lambda He spec (for matched-lambda regression tests).

    All orbitals at the same lambda — should reduce to standard single-
    focal Sturmian basis.
    """
    Z = 2.0
    orbitals = []
    for n in range(1, n_max + 1):
        for l in range(n):
            if l == 0:
                orbitals.append(MultifocalOrbital(
                    n=n, l=l, m=0, lam=lam, label=f"{n}s",
                ))
            else:
                orbitals.append(MultifocalOrbital(
                    n=n, l=l, m=m_p, lam=lam, label=f"{n}{'spdf'[l]}",
                ))
    return MultifocalSpec(orbitals=orbitals, Z_nuc=Z, label=f"he_uniform_lam{lam}")
