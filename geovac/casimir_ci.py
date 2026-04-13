"""
Algebraic Casimir CI: Graph-native He ground-state calculation.

Computes the He ground-state energy entirely from graph-algebraic data:
- Casimir eigenvalues (pure integers from SO(6) representation theory)
- Analytical Slater integrals (exact rationals times k)
- Energy-shell self-consistency: k^2 = -2E

At each n_max, the FCI matrix H(k) = Ck^2 + Bk is a matrix polynomial
in the orbital exponent k. The self-consistent energy E* = -k*^2/2 is
an algebraic number (root of a polynomial with rational coefficients).

Reference: Track DI Sprint 2 design document.
           Paper 7 (S3 density overlap formula).
           Paper 18 (exchange constant taxonomy).
"""

import numpy as np
from typing import Tuple, List, Dict, Optional
from fractions import Fraction
from itertools import combinations
from functools import lru_cache


# ============================================================================
# Exact rational Slater integrals at k_orbital=1
# ============================================================================
#
# Two types of radial Slater integrals are stored:
#
# 1. F^k(a,b) = R^k(aa,bb):  Direct (density-density) integrals.
#    F^k = int int |R_a|^2 (r_<^k / r_>^{k+1}) |R_b|^2 r1^2 r2^2 dr1 dr2
#    Key: (na,la, na,la, nb,lb, nb,lb, k) — but stored with shorthand.
#
# 2. G^k(a,b) = R^k(ab,ab):  Exchange (transition-density) integrals.
#    G^k = int int R_a*R_b (r_<^k/r_>^{k+1}) R_a*R_b r1^2 r2^2 dr1 dr2
#
# All values verified by symbolic integration (sympy) against the
# hydrogenic wavefunctions R_nl(r) = N * (2r/n)^l * L_{n-l-1}^{2l+1}(2r/n)
# * exp(-r/n).
#
# General 4-index storage: _RK4_TABLE[(n1,l1,n3,l3, n2,l2,n4,l4, k)] = Fraction
# where R^k = int int R_{n1l1}*R_{n3l3} (r_<^k/r_>^{k+1}) R_{n2l2}*R_{n4l4} dr

_RK4_TABLE: Dict[Tuple[int, ...], Fraction] = {}


def _add_fk(n1: int, l1: int, n2: int, l2: int, k: int,
            num: int, den: int) -> None:
    """Register F^k(a,b) = R^k(aa,bb) — direct integral (both orderings)."""
    val = Fraction(num, den)
    _RK4_TABLE[(n1, l1, n1, l1, n2, l2, n2, l2, k)] = val
    _RK4_TABLE[(n2, l2, n2, l2, n1, l1, n1, l1, k)] = val


def _add_gk(n1: int, l1: int, n2: int, l2: int, k: int,
            num: int, den: int) -> None:
    """Register G^k(a,b) = R^k(ab,ab) — exchange integral.

    G^k is symmetric: R^k(ab,ab) = R^k(ba,ba) = R^k(ab,ba) = R^k(ba,ab).
    """
    val = Fraction(num, den)
    # All orderings of transition density pairs
    for (a1, b1, a2, b2) in [
        (n1, l1, n2, l2),
        (n2, l2, n1, l1),
    ]:
        for (c1, d1, c2, d2) in [
            (n1, l1, n2, l2),
            (n2, l2, n1, l1),
        ]:
            _RK4_TABLE[(a1, b1, a2, b2, c1, d1, c2, d2, k)] = val


# === F^0 (k=0): direct Coulomb integrals ===
_add_fk(1, 0, 1, 0, 0, 5, 8)
_add_fk(1, 0, 2, 0, 0, 17, 81)
_add_fk(1, 0, 2, 1, 0, 59, 243)
_add_fk(2, 0, 2, 0, 0, 77, 512)
_add_fk(2, 0, 2, 1, 0, 83, 512)
_add_fk(2, 1, 2, 1, 0, 93, 512)
_add_fk(1, 0, 3, 0, 0, 815, 8192)
_add_fk(1, 0, 3, 1, 0, 1783, 16384)
_add_fk(1, 0, 3, 2, 0, 1819, 16384)
_add_fk(2, 0, 3, 0, 0, 32857, 390625)
_add_fk(2, 0, 3, 1, 0, 35161, 390625)
_add_fk(2, 0, 3, 2, 0, 202301, 1953125)
_add_fk(2, 1, 3, 0, 0, 33889, 390625)
_add_fk(2, 1, 3, 1, 0, 36697, 390625)
_add_fk(2, 1, 3, 2, 0, 207677, 1953125)
_add_fk(3, 0, 3, 0, 0, 17, 256)
_add_fk(3, 0, 3, 1, 0, 317, 4608)
_add_fk(3, 0, 3, 2, 0, 337, 4608)
_add_fk(3, 1, 3, 1, 0, 1987, 27648)
_add_fk(3, 1, 3, 2, 0, 709, 9216)
_add_fk(3, 2, 3, 2, 0, 793, 9216)

# === R^1 (k=1): direct integrals ===
_add_fk(1, 0, 1, 0, 1, 3, 8)
_add_fk(1, 0, 2, 0, 1, 2, 27)
_add_fk(1, 0, 2, 1, 1, 8, 81)
_add_fk(2, 0, 2, 0, 1, 51, 512)
_add_fk(2, 0, 2, 1, 1, 53, 512)
_add_fk(2, 1, 2, 1, 1, 185, 1536)
_add_fk(1, 0, 3, 0, 1, 189, 8192)
_add_fk(1, 0, 3, 1, 1, 469, 16384)
_add_fk(1, 0, 3, 2, 1, 1797, 81920)
_add_fk(2, 0, 3, 0, 1, 16146, 390625)
_add_fk(2, 0, 3, 1, 1, 17708, 390625)
_add_fk(2, 0, 3, 2, 1, 122628, 1953125)
_add_fk(2, 1, 3, 0, 1, 15192, 390625)
_add_fk(2, 1, 3, 1, 1, 51998, 1171875)
_add_fk(2, 1, 3, 2, 1, 112606, 1953125)
_add_fk(3, 0, 3, 0, 1, 103, 2304)
_add_fk(3, 0, 3, 1, 1, 629, 13824)
_add_fk(3, 0, 3, 2, 1, 1079, 23040)
_add_fk(3, 1, 3, 1, 1, 3935, 82944)
_add_fk(3, 1, 3, 2, 1, 51, 1024)
_add_fk(3, 2, 3, 2, 1, 2779, 46080)

# === G^0 (k=0): exchange integrals ===
_add_gk(1, 0, 2, 0, 0, 16, 729)
_add_gk(1, 0, 2, 1, 0, 176, 2187)
_add_gk(2, 0, 2, 1, 0, 59, 512)
_add_gk(1, 0, 3, 0, 0, 189, 32768)
_add_gk(1, 0, 3, 1, 0, 1365, 65536)
_add_gk(1, 0, 3, 2, 0, 837, 327680)
_add_gk(2, 0, 3, 0, 0, 73008, 9765625)
_add_gk(2, 0, 3, 1, 0, 10752, 1953125)
_add_gk(2, 0, 3, 2, 0, 3290112, 48828125)
_add_gk(2, 1, 3, 0, 0, 69552, 9765625)
_add_gk(2, 1, 3, 1, 0, 96768, 9765625)
_add_gk(2, 1, 3, 2, 0, 2668032, 48828125)
_add_gk(3, 0, 3, 1, 0, 91, 1536)
_add_gk(3, 0, 3, 2, 0, 617, 23040)
_add_gk(3, 1, 3, 2, 0, 1087, 27648)

# === G^1 (k=1): exchange integrals ===
_add_gk(1, 0, 2, 0, 1, 16, 729)
_add_gk(1, 0, 2, 1, 1, 112, 2187)
_add_gk(2, 0, 2, 1, 1, 45, 512)

# === R^2 (k=2): direct integrals ===
_add_fk(2, 1, 2, 1, 2, 45, 512)  # corrected: was 43/512 (typo), verified algebraically + numerically


def get_rk4(n1: int, l1: int, n3: int, l3: int,
            n2: int, l2: int, n4: int, l4: int, k: int) -> Optional[Fraction]:
    """Look up the exact 4-index R^k at orbital exponent = 1.

    Returns None if not in the table (will fall back to numerical).
    """
    return _RK4_TABLE.get((n1, l1, n3, l3, n2, l2, n4, l4, k))


# ============================================================================
# Off-diagonal 1/r matrix elements: <nl|1/r|n'l> at k_orbital=1
# ============================================================================
#
# When k_orb != Z, the one-body Hamiltonian h = T - Z/r has off-diagonal
# elements in the hydrogenic basis with exponent k:
#   h_{ab} = delta_{ab} * (k^2/(2n^2) - Zk/n^2)  +  (k - Z) * <a|1/r|b>  [WRONG for a!=b]
#
# Actually: h_{ab} = <a|T|b> + <a|-Z/r|b>
# Since |a> are eigenstates of T - k/r with eigenvalue -k^2/(2n_a^2):
#   <a|T - k/r|b> = delta_{ab} * (-k^2/(2n^2))
#   => <a|T|b> = delta_{ab} * (-k^2/(2n^2)) + k * <a|1/r|b>
#   => h_{ab} = <a|T|b> - Z*<a|1/r|b>
#             = delta_{ab} * (-k^2/(2n^2)) + k * <a|1/r|b> - Z * <a|1/r|b>
#             = delta_{ab} * (-k^2/(2n^2)) + (k - Z) * <a|1/r|b>
#
# For k = Z: h_{ab} = delta_{ab} * (-Z^2/(2n^2)) — pure diagonal.
# For k != Z: off-diagonal elements from the (k-Z)/r perturbation.
#
# <nl|1/r|n'l> at k=1 (all verified by symbolic integration):
# These scale with k: <nl|1/r|n'l>_k = k * <nl|1/r|n'l>_{k=1}

_INV_R_TABLE: Dict[Tuple[int, int, int, int], float] = {}


def _add_inv_r(n1: int, l1: int, n2: int, l2: int, val: float) -> None:
    """Register <n1,l1|1/r|n2,l2> at k=1 (symmetric)."""
    _INV_R_TABLE[(n1, l1, n2, l2)] = val
    _INV_R_TABLE[(n2, l2, n1, l1)] = val


# Diagonal: <nl|1/r|nl> = 1/n^2
_add_inv_r(1, 0, 1, 0, 1.0)
_add_inv_r(2, 0, 2, 0, 0.25)
_add_inv_r(2, 1, 2, 1, 0.25)
_add_inv_r(3, 0, 3, 0, 1.0 / 9)
_add_inv_r(3, 1, 3, 1, 1.0 / 9)
_add_inv_r(3, 2, 3, 2, 1.0 / 9)

# Off-diagonal (same l, different n): verified by sympy
_add_inv_r(1, 0, 2, 0, 4 * np.sqrt(2) / 27)          # 4*sqrt(2)/27
_add_inv_r(1, 0, 3, 0, np.sqrt(3) / 16)                # sqrt(3)/16
_add_inv_r(2, 0, 3, 0, 92 * np.sqrt(6) / 3125)         # 92*sqrt(6)/3125
_add_inv_r(2, 1, 3, 1, 192 / 3125)                      # 192/3125


@lru_cache(maxsize=4096)
def _compute_inv_r_numerical(n1: int, l1: int, n2: int, l2: int) -> float:
    """Compute <n1,l1|1/r|n2,l2> at k=1 numerically.

    <nl|1/r|n'l> = integral R_{nl}(r) * (1/r) * R_{n'l}(r) * r^2 dr
                 = integral R_{nl}(r) * R_{n'l}(r) * r dr
    """
    if l1 != l2:
        return 0.0
    from numpy import trapezoid
    nr = 20000
    rmax = max(60.0, 15.0 * max(n1, n2) ** 2)
    r = np.linspace(1e-14, rmax, nr)
    integrand = _hydrogenic_R_numerical(n1, l1, r) * _hydrogenic_R_numerical(n2, l2, r) * r
    return trapezoid(integrand, r)


def _get_inv_r(n1: int, l1: int, n2: int, l2: int) -> float:
    """Get <n1,l1|1/r|n2,l2> at k=1 from table or numerical fallback."""
    val = _INV_R_TABLE.get((n1, l1, n2, l2))
    if val is not None:
        return val
    return _compute_inv_r_numerical(n1, l1, n2, l2)


def get_h1_element(
    n1: int, l1: int, m1: int,
    n2: int, l2: int, m2: int,
    Z: int, k_orb: float,
) -> float:
    """One-body matrix element <n1l1m1|T - Z/r|n2l2m2> at exponent k_orb.

    h_{ab} = delta_{ab} * (-k^2/(2n^2)) + (k - Z) * k * <a|1/r|b>_{k=1}

    The 1/r operator is diagonal in (l,m): nonzero only when l1==l2, m1==m2.
    """
    if l1 != l2 or m1 != m2:
        return 0.0

    result = 0.0
    # Diagonal term from kinetic + full nuclear attraction at exponent k
    if n1 == n2:
        result = -k_orb ** 2 / (2 * n1 ** 2)

    # Off-diagonal correction: (k - Z) * k * <n1l|1/r|n2l>_{k=1}
    inv_r = _get_inv_r(n1, l1, n2, l2)
    result += (k_orb - Z) * k_orb * inv_r

    return result


# ============================================================================
# Gaunt angular coefficients for two-electron integrals
# ============================================================================
#
# The two-electron integral <ab|cd> in the hydrogenic basis is:
#   <ab|cd> = sum_k c_k(la,ma,lc,mc) * c_k(lb,mb,ld,md) * R^k(na la, nb lb, nc lc, nd ld)
#
# where c_k(l1,m1,l2,m2) = (-1)^m1 * sqrt((2l1+1)(2l2+1)) *
#   (l1 k l2; 0 0 0) * (l1 k l2; -m1 m1-m2 m2)
#
# For the DIRECT integral <ab|ab>:
#   <ab|ab> = sum_k a_k(la,ma,lb,mb) * R^k(na la, nb lb)
#   a_k = c_k(la,ma,la,ma) * c_k(lb,mb,lb,mb)
#       = (2la+1)(2lb+1) * (la k la; 0 0 0)^2 * (lb k lb; 0 0 0)^2
#         * (la k la; -ma 0 ma)^2 * (lb k lb; -mb 0 mb)^2    [for direct, m1=m2]
#
# For the EXCHANGE integral <ab|ba>:
#   <ab|ba> = sum_k b_k(la,ma,lb,mb) * R^k(na la, nb lb)  [SAME R^k!]
#   b_k = c_k(la,ma,lb,mb) * c_k(lb,mb,la,ma)

@lru_cache(maxsize=4096)
def _wigner3j(j1: int, j2: int, j3: int, m1: int, m2: int, m3: int) -> float:
    """Wigner 3j symbol for integer arguments."""
    from math import factorial, sqrt
    if m1 + m2 + m3 != 0:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    if j3 > j1 + j2 or j3 < abs(j1 - j2):
        return 0.0

    a, b, c = j1, j2, j3
    triangle_num = factorial(a + b - c) * factorial(a - b + c) * factorial(-a + b + c)
    triangle_den = factorial(a + b + c + 1)

    prefactor_sq = (
        triangle_num
        * factorial(j1 + m1) * factorial(j1 - m1)
        * factorial(j2 + m2) * factorial(j2 - m2)
        * factorial(j3 + m3) * factorial(j3 - m3)
        / triangle_den
    )

    t_min = max(0, j2 - j3 - m1, j1 - j3 + m2)
    t_max = min(j1 + j2 - j3, j1 - m1, j2 + m2)

    s = 0.0
    for t in range(t_min, t_max + 1):
        denom = (
            factorial(t)
            * factorial(j3 - j2 + t + m1)
            * factorial(j3 - j1 + t - m2)
            * factorial(j1 + j2 - j3 - t)
            * factorial(j1 - t - m1)
            * factorial(j2 - t + m2)
        )
        s += (-1) ** t / denom

    phase = (-1) ** (j1 - j2 - m3)
    return phase * sqrt(prefactor_sq) * s


def _gaunt_ck(l1: int, m1: int, l2: int, m2: int, k: int) -> float:
    """Gaunt coefficient c_k(l1,m1,l2,m2).

    c_k = (-1)^m1 * sqrt((2l1+1)(2l2+1)) * (l1 k l2; 0 0 0) * (l1 k l2; -m1 q m2)
    where q = m1 - m2.
    """
    q = m1 - m2
    return ((-1) ** m1
            * np.sqrt((2 * l1 + 1) * (2 * l2 + 1))
            * _wigner3j(l1, k, l2, 0, 0, 0)
            * _wigner3j(l1, k, l2, -m1, q, m2))


def two_electron_integral(
    na: int, la: int, ma: int,
    nb: int, lb: int, mb: int,
    nc: int, lc: int, mc: int,
    nd: int, ld: int, md: int,
    k_orb: float = 1.0,
) -> float:
    """Compute <ab|cd> = sum_k c_k(a,c) * c_k(b,d) * R^k(ac,bd).

    The two-electron integral in the physics (Dirac) convention:
    <ab|cd> = int int phi_a*(r1) phi_b*(r2) (1/r12) phi_c(r1) phi_d(r2) dr1 dr2

    R^k(ac,bd) is the general 4-index radial Slater integral:
    R^k = int int [R_a(r1)*R_c(r1)] (r_<^k/r_>^{k+1}) [R_b(r2)*R_d(r2)] r1^2 r2^2 dr1 dr2

    All integrals scale linearly with the orbital exponent k_orb.

    Parameters
    ----------
    na, la, ma : quantum numbers of orbital a
    nb, lb, mb : quantum numbers of orbital b
    nc, lc, mc : quantum numbers of orbital c
    nd, ld, md : quantum numbers of orbital d
    k_orb : orbital exponent (default 1.0)

    Returns
    -------
    float : the two-electron integral value
    """
    result = 0.0
    # Sum over multipole order k
    # Triangle inequality: |la-lc| <= k <= la+lc AND |lb-ld| <= k <= lb+ld
    k_min = max(abs(la - lc), abs(lb - ld))
    k_max = min(la + lc, lb + ld)

    for k in range(k_min, k_max + 1):
        # Parity selection: la + lc + k must be even, lb + ld + k must be even
        if (la + lc + k) % 2 != 0:
            continue
        if (lb + ld + k) % 2 != 0:
            continue

        ck_ac = _gaunt_ck(la, ma, lc, mc, k)
        ck_bd = _gaunt_ck(lb, mb, ld, md, k)

        if abs(ck_ac) < 1e-15 or abs(ck_bd) < 1e-15:
            continue

        # R^k(ac, bd): 4-index radial Slater integral
        # Try exact table first, fall back to numerical
        rk_exact = get_rk4(na, la, nc, lc, nb, lb, nd, ld, k)
        if rk_exact is not None:
            rk_val = float(rk_exact)
        else:
            rk_val = _compute_rk_numerical(na, la, nc, lc, nb, lb, nd, ld, k)

        result += ck_ac * ck_bd * rk_val

    return result * k_orb


@lru_cache(maxsize=4096)
def _compute_rk_numerical(
    n1: int, l1: int, n3: int, l3: int,
    n2: int, l2: int, n4: int, l4: int,
    k: int,
) -> float:
    """Compute R^k numerically when not in the rational table.

    R^k(n1l1 n3l3, n2l2 n4l4) = int int P_13(r1) (r_<^k/r_>^{k+1}) P_24(r2) dr1 dr2
    where P_13(r) = R_{n1l1}(r) * R_{n3l3}(r) * r^2
    """
    from numpy import trapezoid

    nr = 20000
    rmax = max(60, 15 * max(n1, n2, n3, n4) ** 2)
    r = np.linspace(1e-14, rmax, nr)
    dr = r[1] - r[0]

    P13 = _hydrogenic_R_numerical(n1, l1, r) * _hydrogenic_R_numerical(n3, l3, r) * r ** 2
    P24 = _hydrogenic_R_numerical(n2, l2, r) * _hydrogenic_R_numerical(n4, l4, r) * r ** 2

    Q = np.cumsum(P24 * r ** k * dr)
    T = np.cumsum((P24 * r ** (-(k + 1)) * dr)[::-1])[::-1]
    integrand = P13 * (Q / r ** (k + 1) + T * r ** k)
    return trapezoid(integrand, r)


def _hydrogenic_R_numerical(n: int, l: int, r: np.ndarray) -> np.ndarray:
    """Normalized hydrogenic R_nl(r) at Z=1, numerically."""
    from scipy.special import eval_genlaguerre
    from math import factorial, sqrt

    rho = 2.0 * r / n
    N_sq = (2.0 / n) ** 3 * factorial(n - l - 1) / (2.0 * n * factorial(n + l))
    N = sqrt(N_sq)
    L = eval_genlaguerre(n - l - 1, 2 * l + 1, rho)
    return N * rho ** l * np.exp(-rho / 2) * L


# ============================================================================
# Orbital basis
# ============================================================================

def _build_orbital_basis(
    n_max: int, l_max: Optional[int] = None,
) -> List[Tuple[int, int, int]]:
    """Build the list of spatial orbitals (n, l, m) up to n_max.

    Ordered by n, then l, then m = -l, ..., +l.
    """
    if l_max is None:
        l_max = n_max - 1
    orbitals = []
    for n in range(1, n_max + 1):
        for l in range(min(n, l_max + 1)):
            for m in range(-l, l + 1):
                orbitals.append((n, l, m))
    return orbitals


# ============================================================================
# FCI matrix construction in the singlet spin-adapted basis
# ============================================================================
#
# For two electrons, the singlet (S=0) spatial wavefunction is symmetric:
#   |IJ> = [phi_i(1)*phi_j(2) + phi_j(1)*phi_i(2)] / N_IJ
# where N_IJ = sqrt(2) for i != j, and N_IJ = 1 for i == j.
#
# Matrix elements of H = h(1) + h(2) + 1/r_12:
#
# ONE-BODY (h diagonal in hydrogenic basis, h_ab = delta_ab * eps_a):
#   <IJ|h(1)+h(2)|PQ> = (1/N_IJ N_PQ) * [
#     delta(i,p)*delta(j,q) + delta(i,q)*delta(j,p)
#   ] * (eps_i + eps_j) * 2     ... but factor 2 overcounts, see below
#
# Actually: since h is diagonal, the only nonzero terms arise when
# {i,j} == {p,q} as SETS. Then:
#   <IJ|h|IJ> = eps_i + eps_j  (always, regardless of normalization)
#
# TWO-BODY (<IJ|g|PQ> with g = 1/r_12):
#   <IJ|g|PQ> = (1/N_IJ N_PQ) * [
#     <ij|pq> + <ij|qp> + <ji|pq> + <ji|qp>
#   ]
# where <ab|cd> is the physics-convention integral:
#   <ab|cd> = int phi_a*(1) phi_b*(2) (1/r12) phi_c(1) phi_d(2) d1 d2

def build_fci_matrix(
    Z: int,
    n_max: int,
    k_orb: float,
    l_max: Optional[int] = None,
    m_total: int = 0,
) -> np.ndarray:
    """Build the two-electron singlet FCI matrix at orbital exponent k_orb.

    Parameters
    ----------
    Z : int
        Nuclear charge.
    n_max : int
        Maximum principal quantum number.
    k_orb : float
        Orbital exponent parameter.
    l_max : int, optional
        Maximum angular momentum to include (default: n_max - 1).
    m_total : int
        Total M_L quantum number (default 0).

    Returns
    -------
    H : ndarray of shape (n_configs, n_configs)
    """
    orbitals = _build_orbital_basis(n_max, l_max)
    n_spatial = len(orbitals)

    # Build singlet configurations: pairs (i, j) with i <= j, M_L = m_total
    configs = []
    for i in range(n_spatial):
        for j in range(i, n_spatial):
            if orbitals[i][2] + orbitals[j][2] == m_total:
                configs.append((i, j))

    n_configs = len(configs)
    H = np.zeros((n_configs, n_configs))

    # One-body matrix elements: h_{ab} = <a|T - Z/r|b> at exponent k_orb
    # When k_orb != Z, includes off-diagonal terms from (k-Z)/r perturbation
    def h1(a: int, c: int) -> float:
        """One-body matrix element <a|T-Z/r|c>."""
        na, la, ma = orbitals[a]
        nc, lc, mc = orbitals[c]
        return get_h1_element(na, la, ma, nc, lc, mc, Z, k_orb)

    # Two-electron integrals (physics convention)
    def g_int(a: int, b: int, c: int, d: int) -> float:
        """<ab|cd> physics = int phi_a*(1) phi_b*(2) 1/r12 phi_c(1) phi_d(2)."""
        na, la, ma = orbitals[a]
        nb, lb, mb = orbitals[b]
        nc, lc, mc = orbitals[c]
        nd, ld, md = orbitals[d]
        return two_electron_integral(na, la, ma, nb, lb, mb,
                                     nc, lc, mc, nd, ld, md, k_orb)

    for I in range(n_configs):
        i, j = configs[I]

        for J in range(I, n_configs):
            p, q = configs[J]

            # Enumerate bra and ket product states explicitly.
            # |IJ> = sum_{(a,b) in bra_perms} |a(1)b(2)> / N_IJ
            bra_perms = [(i, j)]
            if i != j:
                bra_perms.append((j, i))
            ket_perms = [(p, q)]
            if p != q:
                ket_perms.append((q, p))

            N_IJ = np.sqrt(float(len(bra_perms)))
            N_PQ = np.sqrt(float(len(ket_perms)))

            me = 0.0
            for a, b in bra_perms:
                for c, d in ket_perms:
                    # One-body: <a(1)b(2)|h1+h2|c(1)d(2)>
                    #   = h1(a,c)*delta(b,d) + delta(a,c)*h1(b,d)
                    if b == d:
                        me += h1(a, c)
                    if a == c:
                        me += h1(b, d)

                    # Two-body: <a(1)b(2)|g|c(1)d(2)> = <ab|cd>
                    me += g_int(a, b, c, d)

            me /= (N_IJ * N_PQ)

            H[I, J] = me
            H[J, I] = me

    return H


# ============================================================================
# Self-consistent solver
# ============================================================================

def solve_self_consistent(
    Z: int,
    n_max: int,
    k_range: Tuple[float, float] = (0.5, 4.0),
    n_scan: int = 200,
    tol: float = 1e-12,
    l_max: Optional[int] = None,
) -> Tuple[float, float, np.ndarray]:
    """Find the self-consistent orbital exponent k* where k^2 = -2E(k).

    Scans k over k_range, computes E_0(k) at each point by diagonalizing
    the FCI matrix, then finds the zero of f(k) = k^2 + 2E_0(k) using
    Brent's method.

    Parameters
    ----------
    Z : int
        Nuclear charge.
    n_max : int
        Maximum principal quantum number.
    k_range : tuple
        Range of k values to scan.
    n_scan : int
        Number of scan points.
    tol : float
        Root-finding tolerance.
    l_max : int, optional
        Maximum angular momentum.

    Returns
    -------
    k_star : float
        Self-consistent orbital exponent.
    E_star : float
        Self-consistent energy = -k_star^2 / 2.
    scan_data : ndarray of shape (n_scan, 3)
        Columns: k, E_0(k), f(k) = k^2 + 2*E_0(k).
    """
    from scipy.optimize import brentq

    k_values = np.linspace(k_range[0], k_range[1], n_scan)
    E_values = np.zeros(n_scan)
    f_values = np.zeros(n_scan)

    for idx, k in enumerate(k_values):
        H = build_fci_matrix(Z, n_max, k, l_max=l_max)
        eigenvalues = np.linalg.eigvalsh(H)
        E_values[idx] = eigenvalues[0]
        f_values[idx] = k ** 2 + 2 * E_values[idx]

    scan_data = np.column_stack([k_values, E_values, f_values])

    # Find sign changes in f(k) — the self-consistency condition
    sign_changes = []
    for idx in range(len(f_values) - 1):
        if f_values[idx] * f_values[idx + 1] < 0:
            sign_changes.append(idx)

    if not sign_changes:
        # No zero crossing found — return the minimum |f| point
        best = np.argmin(np.abs(f_values))
        print(f"[CasimirCI] WARNING: No zero crossing found. Best |f| = {abs(f_values[best]):.2e} at k = {k_values[best]:.4f}")
        return k_values[best], -k_values[best] ** 2 / 2, scan_data

    # Use Brent's method on the first sign change
    idx = sign_changes[0]

    def f_func(k: float) -> float:
        H = build_fci_matrix(Z, n_max, k, l_max=l_max)
        eigenvalues = np.linalg.eigvalsh(H)
        return k ** 2 + 2 * eigenvalues[0]

    k_star = brentq(f_func, k_values[idx], k_values[idx + 1],
                    xtol=tol, rtol=tol)
    E_star = -k_star ** 2 / 2

    return k_star, E_star, scan_data


def solve_variational(
    Z: int,
    n_max: int,
    k_range: Tuple[float, float] = (0.5, 4.0),
    tol: float = 1e-10,
    l_max: Optional[int] = None,
) -> Tuple[float, float]:
    """Find the variational optimal k that minimizes E_0(k).

    Returns
    -------
    k_opt : float
        Optimal orbital exponent.
    E_opt : float
        Minimum energy E_0(k_opt).
    """
    from scipy.optimize import minimize_scalar

    def neg_E(k: float) -> float:
        H = build_fci_matrix(Z, n_max, k, l_max=l_max)
        return np.linalg.eigvalsh(H)[0]

    result = minimize_scalar(neg_E, bounds=k_range, method='bounded',
                             options={'xatol': tol})
    return result.x, result.fun


def convergence_table(
    Z: int,
    n_max_range: range,
    l_max: Optional[int] = None,
) -> List[Dict]:
    """Compute self-consistent energies at each n_max.

    Returns a list of dicts with keys:
      n_max, n_spatial, n_configs, k_star, E_star, E_exact, error_pct
    """
    # Exact He ground state energy
    exact_energies = {
        2: -2.903724,  # He
        1: -0.527751,  # H-
        3: -7.279913,  # Li+
    }
    E_exact = exact_energies.get(Z, None)

    results = []
    for n_max in n_max_range:
        print(f"\n[CasimirCI] n_max={n_max}, Z={Z}")

        # Count basis size
        orbitals = _build_orbital_basis(n_max)
        if l_max is not None:
            orbitals = [(n, l, m) for (n, l, m) in orbitals if l <= l_max]
        n_spatial = len(orbitals)

        # Count singlet M_L=0 configs
        n_configs = 0
        for i in range(n_spatial):
            for j in range(i, n_spatial):
                if orbitals[i][2] + orbitals[j][2] == 0:
                    n_configs += 1

        k_star, E_star, scan_data = solve_self_consistent(
            Z, n_max, l_max=l_max
        )

        error_pct = None
        if E_exact is not None:
            error_pct = abs((E_star - E_exact) / E_exact) * 100

        # Fixed-k energy for comparison
        H_fixed = build_fci_matrix(Z, n_max, float(Z), l_max=l_max)
        E_fixed = np.linalg.eigvalsh(H_fixed)[0]
        error_fixed = None
        if E_exact is not None:
            error_fixed = abs((E_fixed - E_exact) / E_exact) * 100

        # Variational optimal k
        k_var, E_var = solve_variational(Z, n_max, l_max=l_max)
        error_var = None
        if E_exact is not None:
            error_var = abs((E_var - E_exact) / E_exact) * 100

        entry = {
            'n_max': n_max,
            'n_spatial': n_spatial,
            'n_configs': n_configs,
            'k_star': k_star,
            'E_star': E_star,
            'k_var': k_var,
            'E_var': E_var,
            'error_var_pct': error_var,
            'E_fixed_kZ': E_fixed,
            'E_exact': E_exact,
            'error_pct': error_pct,
            'error_fixed_pct': error_fixed,
        }
        results.append(entry)

        print(f"  n_spatial={n_spatial}, n_configs={n_configs}")
        print(f"  k_sc = {k_star:.8f}, E_sc = {E_star:.8f} Ha")
        print(f"  k_var = {k_var:.8f}, E_var = {E_var:.8f} Ha")
        print(f"  E(k=Z) = {E_fixed:.8f} Ha")
        if error_pct is not None:
            print(f"  Error: SC={error_pct:.4f}%, var={error_var:.4f}%, fixed={error_fixed:.4f}%")

    return results


# ============================================================================
# Diagnostic: verify matrix polynomial structure
# ============================================================================

def verify_matrix_polynomial(Z: int, n_max: int, l_max: Optional[int] = None) -> Dict:
    """Verify that H(k) = k^2 * C + k * B by fitting at three k values.

    If the FCI matrix is truly a quadratic polynomial in k, then
    H(k) at any three distinct k values uniquely determines (A, B, C)
    where H(k) = A + Bk + Ck^2. For our case A=0 (no k-independent term).

    Returns
    -------
    dict with keys: 'B', 'C', 'max_residual', 'is_polynomial'
    """
    k_test = [1.0, 1.5, 2.0, 2.5]
    matrices = [build_fci_matrix(Z, n_max, k, l_max=l_max) for k in k_test]
    n = matrices[0].shape[0]

    # Fit H(k) = Ck^2 + Bk (no constant term, since T ~ k^2, V_ne ~ k, V_ee ~ k)
    # At k1, k2, k3: system is overdetermined, fit B and C
    B = np.zeros((n, n))
    C = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            h_vals = [matrices[m][i, j] for m in range(len(k_test))]
            # Fit h = c*k^2 + b*k using least squares
            A_mat = np.array([[k ** 2, k] for k in k_test])
            coeffs, residual, _, _ = np.linalg.lstsq(A_mat, h_vals, rcond=None)
            C[i, j] = coeffs[0]
            B[i, j] = coeffs[1]

    # Verify: reconstruct H at a test k value
    k_verify = 1.73
    H_verify = build_fci_matrix(Z, n_max, k_verify, l_max=l_max)
    H_reconstructed = C * k_verify ** 2 + B * k_verify
    max_residual = np.max(np.abs(H_verify - H_reconstructed))

    return {
        'B': B,
        'C': C,
        'max_residual': max_residual,
        'is_polynomial': max_residual < 1e-10,
        'n_configs': n,
    }


# ============================================================================
# Fast polynomial-based solvers for large n_max
# ============================================================================

def build_fci_polynomial(
    Z: int,
    n_max: int,
    l_max: Optional[int] = None,
    m_total: int = 0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Extract B, C matrices such that H(k) = B*k + C*k^2.

    Builds the FCI matrix at k=1 and k=2, then solves the linear system:
      H(1) = B + C
      H(2) = 2B + 4C
    => C = (H(2) - 2*H(1)) / 2,  B = H(1) - C

    Parameters
    ----------
    Z : int
        Nuclear charge.
    n_max : int
        Maximum principal quantum number.
    l_max : int, optional
        Maximum angular momentum.
    m_total : int
        Total M_L quantum number.

    Returns
    -------
    B : ndarray
        Coefficient of k (one-body nuclear + V_ee scaling).
    C : ndarray
        Coefficient of k^2 (kinetic energy scaling).
    """
    H1 = build_fci_matrix(Z, n_max, 1.0, l_max=l_max, m_total=m_total)
    H2 = build_fci_matrix(Z, n_max, 2.0, l_max=l_max, m_total=m_total)
    C = (H2 - 2.0 * H1) / 2.0
    B = H1 - C
    return B, C


def solve_variational_fast(
    Z: int,
    n_max: int,
    k_range: Tuple[float, float] = (0.5, 4.0),
    tol: float = 1e-10,
    l_max: Optional[int] = None,
    B: Optional[np.ndarray] = None,
    C: Optional[np.ndarray] = None,
) -> Tuple[float, float]:
    """Find the variational optimal k using precomputed polynomial matrices.

    Much faster than solve_variational() for large n_max because the FCI
    matrix is built only twice (at k=1 and k=2), then each k-evaluation
    is just H = B*k + C*k^2 followed by eigensolve.

    Parameters
    ----------
    Z : int
        Nuclear charge.
    n_max : int
        Maximum principal quantum number.
    k_range : tuple
        Search range for k.
    tol : float
        Optimization tolerance.
    l_max : int, optional
        Maximum angular momentum.
    B, C : ndarray, optional
        Precomputed polynomial matrices. If None, computed internally.

    Returns
    -------
    k_opt : float
        Optimal orbital exponent.
    E_opt : float
        Minimum energy E_0(k_opt).
    """
    from scipy.optimize import minimize_scalar

    if B is None or C is None:
        B, C = build_fci_polynomial(Z, n_max, l_max=l_max)

    def E_ground(k: float) -> float:
        H = B * k + C * k ** 2
        return np.linalg.eigvalsh(H)[0]

    result = minimize_scalar(E_ground, bounds=k_range, method='bounded',
                             options={'xatol': tol})
    return result.x, result.fun


def extended_convergence_table(
    Z: int,
    n_max_range: range,
    l_max: Optional[int] = None,
) -> List[Dict]:
    """Compute variational energies for extended n_max range.

    Uses the fast polynomial solver for efficiency at large n_max.

    Returns
    -------
    List[Dict] with keys: n_max, n_spatial, n_configs, k_var, E_var,
                          error_pct, E_exact
    """
    exact_energies = {
        2: -2.903724377,
        1: -0.527751,
        3: -7.279913,
    }
    E_exact = exact_energies.get(Z, None)

    results = []
    for n_max in n_max_range:
        orbitals = _build_orbital_basis(n_max, l_max)
        n_spatial = len(orbitals)

        # Count singlet M_L=0 configs
        n_configs = 0
        for i in range(n_spatial):
            for j in range(i, n_spatial):
                if orbitals[i][2] + orbitals[j][2] == 0:
                    n_configs += 1

        print(f"[CasimirCI] n_max={n_max}: {n_spatial} orbitals, "
              f"{n_configs} configs — building FCI polynomial...")

        import time
        t0 = time.time()
        B, C = build_fci_polynomial(Z, n_max, l_max=l_max)
        t_build = time.time() - t0
        print(f"  FCI polynomial built in {t_build:.1f}s")

        t0 = time.time()
        k_var, E_var = solve_variational_fast(
            Z, n_max, l_max=l_max, B=B, C=C
        )
        t_opt = time.time() - t0
        print(f"  k-optimization in {t_opt:.1f}s")

        error_pct = None
        error_ha = None
        if E_exact is not None:
            error_pct = abs((E_var - E_exact) / E_exact) * 100
            error_ha = E_var - E_exact

        entry = {
            'n_max': n_max,
            'n_spatial': n_spatial,
            'n_configs': n_configs,
            'k_var': float(k_var),
            'E_var': float(E_var),
            'error_pct': error_pct,
            'error_ha': error_ha,
            'E_exact': E_exact,
            'build_time_s': t_build,
            'opt_time_s': t_opt,
        }
        results.append(entry)

        print(f"  k_var={k_var:.6f}, E_var={E_var:.8f} Ha")
        if error_pct is not None:
            print(f"  Error: {error_pct:.4f}% ({error_ha:.6f} Ha)")

    return results


# ============================================================================
# Graph-native CI: hybrid h1 (exact diagonal + graph off-diagonal) + Slater V_ee
# ============================================================================
#
# Sprint 3D extension: build_graph_consistent_fci transforms V_ee into the
# graph eigenbasis, verifying that FCI is invariant under orbital rotation
# (the ~0.075% floor is NOT from basis mismatch).
# ============================================================================

def _build_graph_h1(Z: int, n_max: int) -> Tuple[np.ndarray, List[Tuple[int, int, int]]]:
    """Build the hybrid one-body matrix: exact diagonal + graph off-diagonal.

    Diagonal: -Z^2/(2n^2) for each orbital (exact hydrogenic eigenvalue).
    Off-diagonal: kappa * (-A[i,j]) where A is the graph adjacency matrix
    and kappa = -1/16 (universal topological constant).

    Parameters
    ----------
    Z : int
        Nuclear charge.
    n_max : int
        Maximum principal quantum number.

    Returns
    -------
    h1 : ndarray of shape (n_spatial, n_spatial)
        One-body Hamiltonian in the (n,l,m) basis.
    orbitals : list of (n,l,m) tuples
        Orbital ordering matching the matrix indices.
    """
    from .lattice import GeometricLattice

    lattice = GeometricLattice(max_n=n_max)
    orbitals = list(lattice.states)
    n_spatial = lattice.num_states

    # Exact diagonal: -Z^2/(2n^2)
    h1 = np.zeros((n_spatial, n_spatial))
    for i, (n, l, m) in enumerate(orbitals):
        h1[i, i] = -Z * Z / (2.0 * n * n)

    # Graph off-diagonal: kappa * (-A)
    # kappa = -1/16, so kappa * (-A[i,j]) = (1/16) * A[i,j]
    kappa = -1.0 / 16.0
    A = lattice.adjacency
    if hasattr(A, 'toarray'):
        A_dense = A.toarray()
    else:
        A_dense = np.array(A)

    for i in range(n_spatial):
        for j in range(n_spatial):
            if i != j and abs(A_dense[i, j]) > 1e-15:
                h1[i, j] = kappa * (-A_dense[i, j])

    return h1, orbitals


def build_graph_native_fci(
    Z: int,
    n_max: int,
    m_total: int = 0,
    spin: str = 'singlet',
) -> np.ndarray:
    """Build the graph-native two-electron FCI matrix.

    Uses hybrid h1 (exact diagonal + graph Laplacian off-diagonal)
    with analytical Slater V_ee integrals at k=Z. No free parameters.

    Parameters
    ----------
    Z : int
        Nuclear charge.
    n_max : int
        Maximum principal quantum number.
    m_total : int
        Total M_L quantum number (default 0).
    spin : str
        'singlet' (S=0, symmetric spatial) or 'triplet' (S=1, antisymmetric
        spatial). For triplet, both electrons must be in different orbitals.

    Returns
    -------
    H : ndarray of shape (n_configs, n_configs)
        FCI Hamiltonian matrix.
    """
    if spin not in ('singlet', 'triplet'):
        raise ValueError(f"spin must be 'singlet' or 'triplet', got '{spin}'")

    h1_mat, orbitals = _build_graph_h1(Z, n_max)
    n_spatial = len(orbitals)

    # Build configurations: pairs (i, j) with M_L = m_total
    # Singlet: i <= j (same orbital allowed), spatial symmetric
    # Triplet: i < j (different orbitals only), spatial antisymmetric
    configs = []
    if spin == 'singlet':
        for i in range(n_spatial):
            for j in range(i, n_spatial):
                if orbitals[i][2] + orbitals[j][2] == m_total:
                    configs.append((i, j))
    else:  # triplet
        for i in range(n_spatial):
            for j in range(i + 1, n_spatial):
                if orbitals[i][2] + orbitals[j][2] == m_total:
                    configs.append((i, j))

    n_configs = len(configs)
    if n_configs == 0:
        return np.zeros((0, 0))

    H = np.zeros((n_configs, n_configs))

    # Sign factor for spatial permutation: +1 singlet, -1 triplet
    parity = 1.0 if spin == 'singlet' else -1.0

    # One-body: use h1_mat directly (includes off-diagonal graph terms)
    def h1(a: int, c: int) -> float:
        return h1_mat[a, c]

    # Two-electron integrals at k=Z (physics convention)
    k_orb = float(Z)

    def g_int(a: int, b: int, c: int, d: int) -> float:
        na, la, ma = orbitals[a]
        nb, lb, mb = orbitals[b]
        nc, lc, mc = orbitals[c]
        nd, ld, md = orbitals[d]
        return two_electron_integral(na, la, ma, nb, lb, mb,
                                     nc, lc, mc, nd, ld, md, k_orb)

    for I in range(n_configs):
        i, j = configs[I]
        for J in range(I, n_configs):
            p, q = configs[J]

            # Spatial wavefunction:
            # Singlet: |IJ> = [|ij> + |ji>] / sqrt(2)  (i != j)
            #          |II> = |ii>                       (i == j)
            # Triplet: |IJ> = [|ij> - |ji>] / sqrt(2)  (always i < j)
            #
            # Enumerate product states with sign:
            # bra: (a,b,sign), ket: (c,d,sign)
            bra_perms = [(i, j, 1.0)]
            if i != j:
                bra_perms.append((j, i, parity))
            ket_perms = [(p, q, 1.0)]
            if p != q:
                ket_perms.append((q, p, parity))

            N_IJ = np.sqrt(float(len(bra_perms)))
            N_PQ = np.sqrt(float(len(ket_perms)))

            me = 0.0
            for a, b, s_bra in bra_perms:
                for c, d, s_ket in ket_perms:
                    sign = s_bra * s_ket
                    if b == d:
                        me += sign * h1(a, c)
                    if a == c:
                        me += sign * h1(b, d)
                    me += sign * g_int(a, b, c, d)

            me /= (N_IJ * N_PQ)
            H[I, J] = me
            H[J, I] = me

    return H


def build_graph_consistent_fci(
    Z: int,
    n_max: int,
    m_total: int = 0,
) -> np.ndarray:
    """Build FCI with both h1 and V_ee in the graph eigenbasis.

    Transforms V_ee into the graph eigenbasis so that both one-body
    and two-body operators are expressed consistently.  This verifies
    that FCI is invariant under unitary orbital rotations: the
    eigenvalues should be IDENTICAL to build_graph_native_fci (to
    machine precision), proving the ~0.075% floor is NOT from basis
    mismatch but from cusp/truncation.

    The h1 matrix is block-diagonalized by (l,m) sector to avoid
    mixing degenerate eigenstates across sectors (which numpy's eigh
    can do for exactly degenerate eigenvalues).

    Parameters
    ----------
    Z : int
        Nuclear charge.
    n_max : int
        Maximum principal quantum number (recommended <= 5 for speed;
        the dense integral tensor requires O(n_spatial^4) memory).
    m_total : int
        Total M_L quantum number (default 0).

    Returns
    -------
    H : ndarray of shape (n_configs, n_configs)
        FCI Hamiltonian matrix in graph eigenbasis.
    """
    h1_mat, orbitals = _build_graph_h1(Z, n_max)
    n_spatial = len(orbitals)

    # Block-diagonalize h1 by (l,m) sector to preserve quantum numbers.
    # The graph Laplacian only connects orbitals within the same (l,m)
    # sector (T± transitions change n but preserve l,m).  Full-matrix
    # eigh can mix degenerate states across sectors, corrupting the
    # m-assignment needed for singlet configuration selection.
    sectors: Dict[Tuple[int, int], List[int]] = {}
    for idx, (n, l, m) in enumerate(orbitals):
        key = (l, m)
        if key not in sectors:
            sectors[key] = []
        sectors[key].append(idx)

    U = np.eye(n_spatial)
    evals = np.zeros(n_spatial)
    eigen_m = [0] * n_spatial

    for (l, m), indices in sectors.items():
        if len(indices) == 1:
            idx = indices[0]
            evals[idx] = h1_mat[idx, idx]
            eigen_m[idx] = m
        else:
            block = h1_mat[np.ix_(indices, indices)]
            block_evals, block_U = np.linalg.eigh(block)
            for i_local, i_global in enumerate(indices):
                evals[i_global] = block_evals[i_local]
                eigen_m[i_global] = m
                for j_local, j_global in enumerate(indices):
                    U[j_global, i_global] = block_U[j_local, i_local]

    # Build the hydrogenic 2-electron integral tensor
    k_orb = float(Z)
    g_hydro = np.zeros((n_spatial, n_spatial, n_spatial, n_spatial))
    for i in range(n_spatial):
        ni, li, mi = orbitals[i]
        for j in range(n_spatial):
            nj, lj, mj = orbitals[j]
            for k in range(n_spatial):
                nk, lk, mk = orbitals[k]
                for l_idx in range(n_spatial):
                    nl, ll, ml = orbitals[l_idx]
                    g_hydro[i, j, k, l_idx] = two_electron_integral(
                        ni, li, mi, nj, lj, mj,
                        nk, lk, mk, nl, ll, ml, k_orb)

    # Transform V_ee to graph eigenbasis: 4-index rotation by U
    # g_graph[a,b,c,d] = sum_{ijkl} U[i,a]*U[j,b]*g[i,j,k,l]*U[k,c]*U[l,d]
    g_graph = np.einsum('ia,jb,ijkl,kc,ld->abcd', U, U, g_hydro, U, U,
                        optimize=True)

    # Build singlet configurations in graph eigenbasis
    configs = []
    for i in range(n_spatial):
        for j in range(i, n_spatial):
            if eigen_m[i] + eigen_m[j] == m_total:
                configs.append((i, j))

    n_configs = len(configs)
    H = np.zeros((n_configs, n_configs))

    for I in range(n_configs):
        i, j = configs[I]
        for J in range(I, n_configs):
            p, q = configs[J]

            bra_perms = [(i, j)]
            if i != j:
                bra_perms.append((j, i))
            ket_perms = [(p, q)]
            if p != q:
                ket_perms.append((q, p))

            N_IJ = np.sqrt(float(len(bra_perms)))
            N_PQ = np.sqrt(float(len(ket_perms)))

            me = 0.0
            for a, b in bra_perms:
                for c, d in ket_perms:
                    # One-body: diagonal in eigenbasis
                    if a == c and b == d:
                        me += evals[a] + evals[b]
                    # Two-body: transformed integrals
                    me += g_graph[a, b, c, d]

            me /= (N_IJ * N_PQ)
            H[I, J] = me
            H[J, I] = me

    return H


# ============================================================================
# He spectrum computation
# ============================================================================

# Nonrelativistic infinite-mass reference energies (Ha) from Drake's
# high-precision variational calculations.
# G.W.F. Drake, "High Precision Calculations for Helium" in Springer
# Handbook of Atomic, Molecular, and Optical Physics (2006).
HE_NR_REFERENCE = {
    '1_1S': -2.903724377034,
    '2_1S': -2.145974046,
    '2_3S': -2.175229378,
    '2_1P': -2.123843086,
    '2_3P': -2.133164190,
    '3_1S': -2.061272,
    '3_3S': -2.068690,
    '3_1P': -2.055146,
    '3_3P': -2.058081,
    '3_1D': -2.055636,
    '3_3D': -2.055637,
}


def compute_he_spectrum(
    n_max: int,
    Z: int = 2,
    n_states: int = 3,
    max_L: int = 2,
) -> Dict:
    """Compute the He atomic spectrum from graph-native CI.

    Builds the FCI matrix in each (spin, M_L) sector and extracts
    eigenvalues. Uses M_L = L to select states dominated by total
    angular momentum L. Note: the graph Laplacian's off-diagonal
    elements break SO(3) symmetry, so L is approximate — M_L is
    the exact quantum number.

    Parameters
    ----------
    n_max : int
        Maximum principal quantum number for the orbital basis.
    Z : int
        Nuclear charge (default 2 for He).
    n_states : int
        Number of eigenvalues to extract per sector.
    max_L : int
        Maximum total orbital angular momentum to compute (0=S, 1=P, 2=D).

    Returns
    -------
    result : dict
        'states': dict mapping state label to energy info
        'transitions': dict mapping transition label to energy difference info
        'reference': dict of Drake NR reference values used
        'n_max': int
        'sector_dims': dict mapping sector label to matrix dimension
    """
    states = {}
    sector_dims = {}

    for L in range(min(max_L + 1, n_max)):
        for spin in ['singlet', 'triplet']:
            H = build_graph_native_fci(Z=Z, n_max=n_max, m_total=L,
                                       spin=spin)
            dim = H.shape[0]
            if dim == 0:
                continue

            spin_mult = 1 if spin == 'singlet' else 3
            L_letter = 'SPDFGH'[L]
            sector_label = f'{spin_mult}{L_letter}'
            sector_dims[sector_label] = dim

            evals = np.sort(np.linalg.eigvalsh(H))
            n_extract = min(n_states, len(evals))

            for k in range(n_extract):
                if spin == 'triplet' and L == 0:
                    n_state = k + 2
                elif L > 0:
                    n_state = L + 1 + k
                else:
                    n_state = k + 1

                label = f'{n_state}_{spin_mult}{L_letter}'
                states[label] = evals[k]

    # Compute transition energies
    transitions = {}
    gs_label = '1_1S'
    if gs_label in states:
        transition_pairs = [
            ('1_1S', '2_1S'), ('1_1S', '2_3S'),
            ('1_1S', '2_1P'), ('1_1S', '2_3P'),
            ('2_3S', '2_1S'), ('2_3S', '2_3P'),
        ]
        for s1, s2 in transition_pairs:
            if s1 in states and s2 in states:
                dE = states[s2] - states[s1]
                ref1 = HE_NR_REFERENCE.get(s1)
                ref2 = HE_NR_REFERENCE.get(s2)
                if ref1 is not None and ref2 is not None:
                    dE_exact = ref2 - ref1
                    if abs(dE_exact) > 1e-10:
                        err = abs((dE - dE_exact) / dE_exact) * 100.0
                    else:
                        err = float('inf')
                    transitions[f'{s1}->{s2}'] = {
                        'dE_computed': dE,
                        'dE_exact': dE_exact,
                        'error_pct': err,
                    }

    # Add reference comparison for absolute energies
    state_info = {}
    for label, E in states.items():
        ref = HE_NR_REFERENCE.get(label)
        if ref is not None:
            state_info[label] = {
                'energy': E,
                'reference': ref,
                'error_pct': abs((E - ref) / ref) * 100.0,
            }
        else:
            state_info[label] = {'energy': E, 'reference': None, 'error_pct': None}

    return {
        'states': state_info,
        'transitions': transitions,
        'reference': HE_NR_REFERENCE,
        'n_max': n_max,
        'sector_dims': sector_dims,
    }
