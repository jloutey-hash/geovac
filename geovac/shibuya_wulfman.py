"""
Shibuya-Wulfman Two-Center Nuclear Attraction Integrals
=======================================================

Computes cross-center nuclear attraction matrix elements:

    I^{AB}_{nlm,n'l'm'} = <psi_{nlm}^A | (-Z_B / |r - R_B|) | psi_{n'l'm'}^A>

where psi_{nlm}^A are hydrogenic orbitals centered on nucleus A, and the
potential is from nucleus B located at distance R_AB along the z-axis.

The cross-center potential is expanded in a multipole series:

    1/|r - R_B| = sum_L (r_<^L / r_>^{L+1}) * P_L(cos theta)

For R_B along z-axis, the angular integral simplifies:
    - m = m' (diagonal in magnetic quantum number)
    - Selection rules: |l-l'| <= L <= l+l', l+L+l' even

Author: GeoVac Development Team (Track CD, Phase 1A)
Date: April 2026
"""

from math import factorial
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.special import gammainc, gammaincc, genlaguerre

from geovac.composed_qubit import (
    _enumerate_states,
    _radial_wf_grid,
    _wigner3j,
)


# ---------------------------------------------------------------------------
# Rotation helpers for non-collinear nuclei
# ---------------------------------------------------------------------------


def _build_rotation_matrix_l1(
    direction: Tuple[float, float, float],
) -> np.ndarray:
    """
    Build 3x3 real rotation matrix for l=1 real spherical harmonics.

    The real harmonics are ordered (m=-1, m=0, m=+1) = (y, z, x).
    The rotation maps the z-axis onto the given direction unit vector.

    Parameters
    ----------
    direction : (nx, ny, nz)
        Unit vector specifying nucleus direction from orbital center.

    Returns
    -------
    np.ndarray
        3x3 rotation matrix D^1 in the (y, z, x) basis.
    """
    nx, ny, nz = direction
    theta = np.arccos(np.clip(nz, -1.0, 1.0))
    phi = np.arctan2(ny, nx)

    ct, st = np.cos(theta), np.sin(theta)
    cp, sp = np.cos(phi), np.sin(phi)

    # Cartesian rotation R_xyz = R_z(phi) @ R_y(theta)
    # in (x, y, z) row/column order
    R_xyz = np.array([
        [cp * ct, -sp, cp * st],
        [sp * ct,  cp, sp * st],
        [-st,      0.0, ct],
    ])

    # Permute from (x, y, z) to (y, z, x) = (m=-1, m=0, m=+1)
    perm = [1, 2, 0]
    return R_xyz[np.ix_(perm, perm)]


def _build_block_rotation_matrix(
    states: List[Tuple[int, int, int]],
    direction: Tuple[float, float, float],
) -> np.ndarray:
    """
    Build block-diagonal rotation matrix for a set of orbitals.

    For l=0 orbitals the block is 1x1 identity.  For each group of l=1
    orbitals (m=-1, 0, +1 with the same n) the block is the 3x3 real
    spherical harmonic rotation matrix D^1.

    Parameters
    ----------
    states : list of (n, l, m)
        Orbital quantum numbers in canonical order.
    direction : (nx, ny, nz)
        Unit vector for the nucleus direction.

    Returns
    -------
    np.ndarray
        Block-diagonal rotation matrix of shape (N, N).

    Raises
    ------
    NotImplementedError
        If any orbital has l >= 2 (requires Ivanic-Ruedenberg recursion).
    """
    N = len(states)
    D = np.eye(N)

    # Pre-compute D^1 only if l=1 orbitals exist
    D1: Optional[np.ndarray] = None
    max_l = max(l for _, l, _ in states) if states else 0
    if max_l >= 2:
        raise NotImplementedError(
            f"Rotation for l={max_l} not implemented "
            "(need Ivanic-Ruedenberg recursion)"
        )
    if max_l >= 1:
        D1 = _build_rotation_matrix_l1(direction)

    i = 0
    while i < N:
        _, l, m = states[i]
        if l == 0:
            # D^0 = 1, already identity
            i += 1
        elif l == 1:
            assert i + 2 < N, "Incomplete l=1 block"
            assert (
                states[i][2] == -1
                and states[i + 1][2] == 0
                and states[i + 2][2] == 1
            ), f"Expected m=-1,0,+1 at index {i}, got {[states[i+k][2] for k in range(3)]}"
            D[i : i + 3, i : i + 3] = D1
            i += 3
        else:
            raise NotImplementedError(
                f"Rotation for l={l} not implemented "
                "(need Ivanic-Ruedenberg recursion)"
            )

    return D


def _angular_coefficient(l1: int, m1: int, l2: int, m2: int, L: int) -> float:
    """
    Angular coefficient for multipole expansion of cross-center V_ne
    with the second nucleus along the z-axis.

    The multipole expansion uses Legendre polynomials P_L(cos θ).
    Converting P_L = sqrt(4π/(2L+1)) Y_{L,0} and using the Gaunt
    integral for ∫ Y_{l1,m1}* Y_{L,0} Y_{l2,m2} dΩ gives:

        A_L = (-1)^m * sqrt((2l1+1)(2l2+1)) * (l1 L l2; 0 0 0) * (l1 L l2; -m 0 m)

    Returns 0 if m1 != m2 (M=0 selection rule for z-axis potential).
    """
    if m1 != m2:
        return 0.0

    # Triangle inequality
    if abs(l1 - l2) > L or L > l1 + l2:
        return 0.0

    # Parity: l1 + L + l2 must be even
    if (l1 + L + l2) % 2 != 0:
        return 0.0

    m = m1  # m1 == m2

    w1 = _wigner3j(l1, L, l2, 0, 0, 0)
    if abs(w1) < 1e-15:
        return 0.0

    w2 = _wigner3j(l1, L, l2, -m, 0, m)
    if abs(w2) < 1e-15:
        return 0.0

    prefactor = (-1) ** m * np.sqrt((2 * l1 + 1) * (2 * l2 + 1))
    return prefactor * w1 * w2


def _hydrogenic_poly_coeffs(Z: float, n: int, l: int) -> Tuple[np.ndarray, float]:
    """
    Decompose R_nl(r) = exp(-alpha*r) * sum_k c_k * r^k.

    Returns (coeffs, alpha) where coeffs[k] is the coefficient of r^k
    and alpha = Z/n is the exponential decay rate.

    Uses the same wavefunction as ``_radial_wf_grid`` but computes the
    exact normalization analytically via ∫₀^∞ r^j e^{-βr} dr = j!/β^{j+1}.
    """
    alpha = Z / n
    rho_scale = 2.0 * Z / n  # rho = rho_scale * r

    # Get Laguerre polynomial coefficients: L_{n-l-1}^{2l+1}(x)
    lag = genlaguerre(n - l - 1, 2 * l + 1)
    lag_coeffs = lag.coef[::-1]  # lag_coeffs[k] is coeff of x^k

    # Un-normalized wavefunction: rho^l * exp(-rho/2) * L(rho)
    # = (rho_scale*r)^l * exp(-alpha*r) * sum_k lag_coeffs[k] * (rho_scale*r)^k
    # = exp(-alpha*r) * sum_k [rho_scale^{l+k} * lag_coeffs[k]] * r^{l+k}

    degree = n - l - 1
    max_power = l + degree
    raw_coeffs = np.zeros(max_power + 1)
    for k in range(degree + 1):
        power = l + k
        raw_coeffs[power] += rho_scale ** (l + k) * lag_coeffs[k]

    # Compute exact normalization: ∫₀^∞ |R|² r² dr = 1
    # |R|² r² = [sum c_i r^i]² * e^{-2αr} * r²
    # = sum_{i,j} c_i c_j r^{i+j+2} e^{-2αr}
    beta = 2.0 * alpha
    norm_sq = 0.0
    for i in range(len(raw_coeffs)):
        if abs(raw_coeffs[i]) < 1e-30:
            continue
        for j in range(len(raw_coeffs)):
            if abs(raw_coeffs[j]) < 1e-30:
                continue
            p = i + j + 2
            norm_sq += raw_coeffs[i] * raw_coeffs[j] * factorial(p) / beta ** (p + 1)

    N = 1.0 / np.sqrt(norm_sq)
    return raw_coeffs * N, alpha


def _poly_product(c1: np.ndarray, c2: np.ndarray) -> np.ndarray:
    """Multiply two polynomial coefficient arrays (convolution)."""
    return np.convolve(c1, c2)


def _split_integral_analytical(
    coeffs: np.ndarray,
    alpha: float,
    p_inner: int,
    p_outer: int,
    L: int,
    R_AB: float,
) -> float:
    """
    Evaluate the split-region integral analytically:

    R_L = (1/R^{L+1}) * int_0^R [sum c_k r^k] * e^{-alpha*r} * r^{p_inner} dr
        + R^L * int_R^inf [sum c_k r^k] * e^{-alpha*r} * r^{p_outer} dr

    Uses the incomplete gamma function:
        int_0^R r^j e^{-alpha*r} dr = j! / alpha^{j+1} * gammainc(j+1, alpha*R)
        int_R^inf r^j e^{-alpha*r} dr = j! / alpha^{j+1} * gammaincc(j+1, alpha*R)

    where gammainc/gammaincc are REGULARIZED (divided by Gamma(j+1) = j!),
    so the actual formula is:
        int_0^R r^j e^{-alpha*r} dr = Gamma(j+1) / alpha^{j+1} * gammainc(j+1, alpha*R)
    """
    x = alpha * R_AB  # argument for incomplete gamma

    I_inner = 0.0
    I_outer = 0.0

    for k in range(len(coeffs)):
        if abs(coeffs[k]) < 1e-30:
            continue

        # Inner integral: power is k + p_inner
        j_in = k + p_inner
        # Gamma(j+1)/alpha^{j+1} * gammainc(j+1, x)
        # = j!/alpha^{j+1} * gammainc(j+1, x)  [since j is integer]
        scale_in = factorial(j_in) / alpha ** (j_in + 1)
        I_inner += coeffs[k] * scale_in * float(gammainc(j_in + 1, x))

        # Outer integral: power is k + p_outer
        j_out = k + p_outer
        scale_out = factorial(j_out) / alpha ** (j_out + 1)
        I_outer += coeffs[k] * scale_out * float(gammaincc(j_out + 1, x))

    return I_inner / R_AB ** (L + 1) + R_AB ** L * I_outer


def _radial_split_integral(
    Z_orb: float,
    n1: int, l1: int,
    n2: int, l2: int,
    L: int,
    R_AB: float,
    n_grid: int = 4000,
) -> float:
    """
    Compute the split-region radial integral analytically:

    R_L = (1/R^{L+1}) int_0^R R_nl(r) R_n'l'(r) r^{L+2} dr
        + R^L int_R^inf R_nl(r) R_n'l'(r) r^{1-L} dr

    Uses incomplete gamma functions for machine-precision evaluation.
    The n_grid parameter is ignored (kept for API compatibility).
    """
    c1, alpha1 = _hydrogenic_poly_coeffs(Z_orb, n1, l1)
    c2, alpha2 = _hydrogenic_poly_coeffs(Z_orb, n2, l2)

    # Product polynomial: R1*R2 = exp(-(alpha1+alpha2)*r) * sum c_k r^k
    prod_coeffs = _poly_product(c1, c2)
    alpha_total = alpha1 + alpha2

    # Inner integral has extra r^{L+2}, outer has r^{1-L}
    return _split_integral_analytical(
        prod_coeffs, alpha_total, L + 2, 1 - L, L, R_AB,
    )


def _radial_split_integral_grid(
    Z_orb: float,
    n1: int, l1: int,
    n2: int, l2: int,
    L: int,
    R_AB: float,
    n_grid: int = 4000,
) -> float:
    """
    Grid-based split-region radial integral (legacy, for cross-validation).

    Uses np.trapezoid on a uniform grid. O(1/n_grid) convergence.
    """
    r_max = max(80.0 / max(Z_orb, 0.5), 4.0 * R_AB)
    r_full = np.linspace(0, r_max, n_grid + 1)[1:]
    R1_full = _radial_wf_grid(Z_orb, n1, l1, r_full)
    R2_full = _radial_wf_grid(Z_orb, n2, l2, r_full)

    mask_inner = r_full <= R_AB
    mask_outer = r_full > R_AB

    r_in = r_full[mask_inner]
    if len(r_in) > 1:
        integrand_in = R1_full[mask_inner] * R2_full[mask_inner] * r_in ** (L + 2)
        I_inner = np.trapezoid(integrand_in, r_in)
    else:
        I_inner = 0.0

    r_out = r_full[mask_outer]
    if len(r_out) > 1:
        integrand_out = R1_full[mask_outer] * R2_full[mask_outer] * r_out ** (1 - L)
        I_outer = np.trapezoid(integrand_out, r_out)
    else:
        I_outer = 0.0

    return I_inner / R_AB ** (L + 1) + R_AB ** L * I_outer


def compute_cross_center_vne_element(
    Z_orb: float,
    n1: int, l1: int, m1: int,
    n2: int, l2: int, m2: int,
    Z_nuc: float,
    R_AB: float,
    L_max: int,
    n_grid: int = 4000,
    nuc_parity: int = 1,
) -> float:
    """
    Single matrix element of cross-center nuclear attraction:

    <psi_{n1,l1,m1}^A | (-Z_B / |r - R_B|) | psi_{n2,l2,m2}^A>

    Parameters
    ----------
    Z_orb : float
        Nuclear charge for the orbital wavefunctions.
    n1, l1, m1 : int
        Quantum numbers for the bra orbital.
    n2, l2, m2 : int
        Quantum numbers for the ket orbital.
    Z_nuc : float
        Charge of the off-center nucleus.
    R_AB : float
        Internuclear distance (bohr).
    L_max : int
        Maximum multipole order in the expansion.
    n_grid : int
        Number of radial grid points (ignored for analytical integrals).
    nuc_parity : int
        +1 if nucleus is in the +z direction from the orbital center,
        -1 if in the -z direction. Applies (-1)^L factor to each
        multipole term for nuclei at theta_B = pi. Default +1.

    Returns
    -------
    float
        Matrix element value in Hartree.
    """
    # m-diagonal selection rule
    if m1 != m2:
        return 0.0

    total = 0.0
    for L in range(0, L_max + 1):
        ang = _angular_coefficient(l1, m1, l2, m2, L)
        if abs(ang) < 1e-15:
            continue
        rad = _radial_split_integral(Z_orb, n1, l1, n2, l2, L, R_AB, n_grid)
        total += nuc_parity ** L * ang * rad

    return -Z_nuc * total


def compute_cross_center_vne(
    Z_orb: float,
    states: List[Tuple[int, int, int]],
    Z_nuc: float,
    R_AB: float,
    L_max: int,
    n_grid: int = 4000,
    nuc_parity: int = 1,
    direction: Optional[Tuple[float, float, float]] = None,
) -> np.ndarray:
    """
    Compute the full cross-center V_ne matrix for a set of orbitals.

    Parameters
    ----------
    Z_orb : float
        Nuclear charge for the orbital wavefunctions.
    states : list of (n, l, m) tuples
        Orbital quantum numbers.
    Z_nuc : float
        Charge of the off-center nucleus.
    R_AB : float
        Internuclear distance (bohr).
    L_max : int
        Maximum multipole order.
    n_grid : int
        Number of radial grid points (ignored for analytical integrals).
    nuc_parity : int
        +1 if nucleus is in +z direction, -1 for -z. Default +1.
        Ignored when ``direction`` is provided.
    direction : tuple of (nx, ny, nz), optional
        Unit vector from orbital center to off-center nucleus.
        When provided, takes precedence over ``nuc_parity``.
        The matrix is computed in the z-frame (nuc_parity=+1) and then
        rotated via a block-diagonal real spherical harmonic rotation
        D @ V_z @ D^T.

    Returns
    -------
    np.ndarray
        V_ne matrix of shape (N_orb, N_orb).
    """
    # --- Direction-based rotation approach ---
    if direction is not None:
        # Compute V_ne as if the nucleus is along +z
        vne_z = compute_cross_center_vne(
            Z_orb, states, Z_nuc, R_AB, L_max,
            n_grid=n_grid, nuc_parity=1, direction=None,
        )
        # Build block-diagonal rotation matrix and rotate
        D = _build_block_rotation_matrix(states, direction)
        return D @ vne_z @ D.T

    # --- Original z-axis code (nuc_parity path) ---
    N = len(states)
    vne = np.zeros((N, N))

    # Pre-compute unique (n,l) pairs
    unique_nl = sorted(set((n, l) for n, l, m in states))

    # Cache radial split integrals keyed by (n1,l1,n2,l2,L)
    radial_cache: Dict[Tuple[int, int, int, int, int], float] = {}
    for n1, l1 in unique_nl:
        for n2, l2 in unique_nl:
            for L in range(0, L_max + 1):
                if (l1 + L + l2) % 2 != 0:
                    continue
                if abs(l1 - l2) > L or L > l1 + l2:
                    continue
                key = (n1, l1, n2, l2, L)
                if key not in radial_cache:
                    radial_cache[key] = _radial_split_integral(
                        Z_orb, n1, l1, n2, l2, L, R_AB, n_grid
                    )

    # Build matrix
    for i, (n1, l1, m1) in enumerate(states):
        for j, (n2, l2, m2) in enumerate(states):
            if m1 != m2:
                continue

            val = 0.0
            for L in range(0, L_max + 1):
                if (l1 + L + l2) % 2 != 0:
                    continue
                if abs(l1 - l2) > L or L > l1 + l2:
                    continue
                ang = _angular_coefficient(l1, m1, l2, m2, L)
                if abs(ang) < 1e-15:
                    continue
                key = (n1, l1, n2, l2, L)
                rad = radial_cache.get(key, 0.0)
                val += nuc_parity ** L * ang * rad

            vne[i, j] = -Z_nuc * val

    return vne


def convergence_study(
    Z_orb: float,
    n1: int, l1: int, m1: int,
    n2: int, l2: int, m2: int,
    Z_nuc: float,
    R_AB: float,
    L_max_range: range,
    n_grid: int = 4000,
) -> Dict[int, float]:
    """
    Convergence study: compute cross-center V_ne element vs L_max.

    Parameters
    ----------
    Z_orb : float
        Nuclear charge for the orbital wavefunctions.
    n1, l1, m1 : int
        Bra orbital quantum numbers.
    n2, l2, m2 : int
        Ket orbital quantum numbers.
    Z_nuc : float
        Off-center nuclear charge.
    R_AB : float
        Internuclear distance (bohr).
    L_max_range : range
        Range of L_max values to test.
    n_grid : int
        Radial grid points.

    Returns
    -------
    dict
        Mapping L_max -> integral value.
    """
    results: Dict[int, float] = {}
    for L_max in L_max_range:
        val = compute_cross_center_vne_element(
            Z_orb, n1, l1, m1, n2, l2, m2, Z_nuc, R_AB, L_max, n_grid
        )
        results[L_max] = val
    return results
