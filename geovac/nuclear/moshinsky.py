"""
Moshinsky-Talmi Brackets for two-body HO transformation.

Transforms between lab-frame two-particle HO states and center-of-mass +
relative coordinate states for equal-mass particles (nucleon-nucleon).

Method: for the HO, the 3D Moshinsky bracket factorizes via Cartesian
decomposition. We use the 1D Moshinsky bracket (analytical, exact) and
reconstruct the 3D brackets by expanding spherical HO states in Cartesian
HO states, then applying 1D brackets coordinate-by-coordinate.

For practical purposes at N_shells <= 3, the Cartesian decomposition is
both exact and fast.

References:
  - Moshinsky, Nucl. Phys. 13, 104 (1959)
  - Suhonen, "From Nucleons to Nucleus" (2007), Chapter 7

Author: GeoVac Development Team
Date: April 2026
"""

from __future__ import annotations

from functools import lru_cache
from math import factorial, sqrt as msqrt
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.special import gamma as gamma_fn, eval_genlaguerre


# ---------------------------------------------------------------------------
# 1D Moshinsky bracket (analytical, exact)
# ---------------------------------------------------------------------------

@lru_cache(maxsize=4096)
def _moshinsky_1d(n1: int, n2: int, nc: int, nr: int) -> float:
    """
    1D Moshinsky bracket <n1, n2 | nc, nr> for equal-mass HO.

    The CM/rel transformation in 1D:
    x_cm = (x1 + x2)/sqrt(2), x_rel = (x1 - x2)/sqrt(2)

    Using creation operators:
    a_cm^dag = (a1^dag + a2^dag)/sqrt(2)
    a_rel^dag = (a1^dag - a2^dag)/sqrt(2)

    The bracket is:
    <n1, n2 | nc, nr> = sqrt(n1! n2! / (nc! nr!)) / 2^{N/2}
      * sum_{j} C(nc,j) C(nr, n1-j) (-1)^{nr-n1+j}

    where N = n1+n2 = nc+nr (energy conservation).
    """
    if n1 + n2 != nc + nr:
        return 0.0
    if any(x < 0 for x in [n1, n2, nc, nr]):
        return 0.0

    N = n1 + n2

    prefactor = msqrt(factorial(n1) * factorial(n2) /
                      (factorial(nc) * factorial(nr))) / (2.0**(N/2.0))

    total = 0.0
    j_min = max(0, n1 - nr)
    j_max = min(nc, n1)
    for j in range(j_min, j_max + 1):
        k = n1 - j
        binom_nc_j = factorial(nc) // (factorial(j) * factorial(nc - j))
        binom_nr_k = factorial(nr) // (factorial(k) * factorial(nr - k))
        sign = (-1)**(nr - k)
        total += binom_nc_j * binom_nr_k * sign

    return prefactor * total


# ---------------------------------------------------------------------------
# Spherical-to-Cartesian transformation for HO states
# ---------------------------------------------------------------------------

@lru_cache(maxsize=2048)
def _spherical_to_cartesian(n_r: int, l: int, m: int) -> Tuple[Tuple[Tuple[int,int,int], ...], Tuple[complex, ...]]:
    """
    Expand |n_r, l, m> in Cartesian HO basis |n_x, n_y, n_z>.

    Uses the numerical overlap via Gauss-Hermite quadrature in 3D.
    Returns (keys, coeffs) where keys are (nx,ny,nz) tuples.
    """
    N = 2*n_r + l
    n_pts = max(N + 5, 12)

    from numpy.polynomial.hermite import hermgauss
    nodes, weights = hermgauss(n_pts)

    # Enumerate Cartesian states
    cart_states = []
    for nx in range(N + 1):
        for ny in range(N - nx + 1):
            nz = N - nx - ny
            cart_states.append((nx, ny, nz))

    # Compute overlaps <nx,ny,nz | n_r, l, m> via 3D Gauss-Hermite
    overlaps = []
    for nx, ny, nz in cart_states:
        overlap = 0.0 + 0.0j
        for xi, wi in zip(nodes, weights):
            hx = _hermite_norm_val(nx, xi)
            for yi, wy in zip(nodes, weights):
                hy = _hermite_norm_val(ny, yi)
                for zi, wz in zip(nodes, weights):
                    hz = _hermite_norm_val(nz, zi)
                    cart_val = hx * hy * hz

                    r_sq = xi**2 + yi**2 + zi**2
                    r = np.sqrt(r_sq)

                    if r < 1e-15:
                        if l == 0:
                            norm_sq = 2.0 * factorial(n_r) / gamma_fn(n_r + 1.5)
                            L_at_0 = float(eval_genlaguerre(n_r, 0.5, 0.0))
                            R_val = msqrt(norm_sq) * L_at_0
                            Ylm = 1.0 / msqrt(4*np.pi)
                            sph_val = R_val * Ylm
                        else:
                            sph_val = 0.0
                    else:
                        # Radial part (without exp(-r^2/2), absorbed by weight)
                        norm_sq = 2.0 * factorial(n_r) / gamma_fn(n_r + l + 1.5)
                        L_poly = float(eval_genlaguerre(n_r, l + 0.5, r_sq))
                        R_val = msqrt(norm_sq) * r**l * L_poly

                        # Spherical harmonic (manual implementation)
                        cos_theta = zi / r
                        sin_theta = msqrt(max(0.0, 1.0 - cos_theta**2))
                        phi_angle = np.arctan2(yi, xi)
                        Ylm = _real_sph_harm_complex(l, m, cos_theta, sin_theta, phi_angle)

                        sph_val = R_val * Ylm

                    overlap += wi * wy * wz * cart_val * np.conj(sph_val)

        overlaps.append(complex(overlap))

    # Filter small values
    keys = []
    coeffs = []
    for (nx, ny, nz), c in zip(cart_states, overlaps):
        if abs(c) > 1e-12:
            keys.append((nx, ny, nz))
            coeffs.append(c)

    return tuple(keys), tuple(coeffs)


def _hermite_norm_val(n: int, x: float) -> float:
    """Evaluate normalized Hermite function H_n(x)/sqrt(2^n n! sqrt(pi)) at x."""
    # Compute H_n(x) via recurrence: H_0=1, H_1=2x, H_{n+1}=2x H_n - 2n H_{n-1}
    if n == 0:
        H = 1.0
    elif n == 1:
        H = 2.0 * x
    else:
        H_prev = 1.0
        H_curr = 2.0 * x
        for k in range(2, n + 1):
            H_next = 2.0 * x * H_curr - 2.0 * (k - 1) * H_prev
            H_prev = H_curr
            H_curr = H_next
        H = H_curr
    norm = msqrt(2**n * factorial(n) * msqrt(np.pi))
    return H / norm


def _real_sph_harm_complex(l: int, m: int, cos_theta: float, sin_theta: float, phi: float) -> complex:
    """
    Complex spherical harmonic Y_l^m(theta, phi).

    Y_l^m = (-1)^m sqrt((2l+1)/(4pi) * (l-m)!/(l+m)!) * P_l^m(cos_theta) * exp(i*m*phi)

    For m < 0: Y_l^{-|m|} = (-1)^m conj(Y_l^{|m|})
    """
    if abs(m) > l:
        return 0.0 + 0.0j

    # Compute associated Legendre polynomial P_l^|m|(cos_theta)
    am = abs(m)
    Plm = _assoc_legendre(l, am, cos_theta, sin_theta)

    # Normalization
    norm = msqrt((2*l + 1) / (4*np.pi) * factorial(l - am) / factorial(l + am))

    # Phase convention: Condon-Shortley
    if m >= 0:
        phase = (-1)**m
        return phase * norm * Plm * np.exp(1j * m * phi)
    else:
        # Y_l^{-|m|} = (-1)^|m| conj(Y_l^{|m|})
        phase_pos = (-1)**am
        Y_pos = phase_pos * norm * Plm * np.exp(1j * am * phi)
        return (-1)**am * np.conj(Y_pos)


def _assoc_legendre(l: int, m: int, cos_theta: float, sin_theta: float) -> float:
    """Associated Legendre polynomial P_l^m(x) where x = cos(theta)."""
    if m > l:
        return 0.0

    # Start with P_m^m = (-1)^m (2m-1)!! sin^m(theta)
    # (without (-1)^m Condon-Shortley phase — that's in the Y_lm)
    P_mm = 1.0
    for i in range(1, m + 1):
        P_mm *= (2*i - 1) * sin_theta

    if l == m:
        return P_mm

    # P_{m+1}^m = cos(theta) * (2m+1) * P_m^m
    P_mp1_m = cos_theta * (2*m + 1) * P_mm
    if l == m + 1:
        return P_mp1_m

    # Recurrence: (l-m) P_l^m = (2l-1) cos(theta) P_{l-1}^m - (l+m-1) P_{l-2}^m
    P_prev = P_mm
    P_curr = P_mp1_m
    for ll in range(m + 2, l + 1):
        P_next = ((2*ll - 1) * cos_theta * P_curr - (ll + m - 1) * P_prev) / (ll - m)
        P_prev = P_curr
        P_curr = P_next

    return P_curr


# ---------------------------------------------------------------------------
# CG coefficient (integer angular momenta)
# ---------------------------------------------------------------------------

def _clebsch_gordan_int(
    j1: int, m1: int, j2: int, m2: int, J: int, M: int,
) -> float:
    """CG coefficient for integer angular momenta."""
    if m1 + m2 != M:
        return 0.0
    from geovac.angular_integrals import wigner3j
    return (-1)**(j1 - j2 + M) * msqrt(2*J + 1) * wigner3j(j1, j2, J, m1, m2, -M)


# ---------------------------------------------------------------------------
# 3D Moshinsky bracket via Cartesian factorization
# ---------------------------------------------------------------------------

@lru_cache(maxsize=8192)
def moshinsky_bracket(
    n1: int, l1: int,
    n2: int, l2: int,
    n_cm: int, l_cm: int,
    n_rel: int, l_rel: int,
    L: int,
) -> float:
    """
    Compute the Moshinsky bracket <n1 l1 n2 l2 L | n_cm l_cm n_rel l_rel L>.

    For equal-mass particles. Uses Cartesian factorization: expand spherical
    HO states in Cartesian HO states, apply 1D Moshinsky brackets to each
    Cartesian direction, then project back.

    Parameters
    ----------
    n1, l1, n2, l2 : int
        Lab-frame HO quantum numbers.
    n_cm, l_cm, n_rel, l_rel : int
        CM + relative HO quantum numbers.
    L : int
        Total orbital angular momentum.

    Returns
    -------
    float
        The Moshinsky bracket value.
    """
    # Energy conservation
    N_lab = 2*n1 + l1 + 2*n2 + l2
    N_cm_rel = 2*n_cm + l_cm + 2*n_rel + l_rel
    if N_lab != N_cm_rel:
        return 0.0

    # Angular momentum triangle checks
    if L < abs(l1 - l2) or L > l1 + l2:
        return 0.0
    if L < abs(l_cm - l_rel) or L > l_cm + l_rel:
        return 0.0

    # Quantum number validity
    if any(x < 0 for x in [n1, l1, n2, l2, n_cm, l_cm, n_rel, l_rel, L]):
        return 0.0

    # Compute via block matrix
    data = _moshinsky_block(N_lab, L)
    if data is None:
        return 0.0

    lab_idx, cm_idx, U = data

    lab_key = (n1, l1, n2, l2)
    cm_key = (n_cm, l_cm, n_rel, l_rel)

    if lab_key not in lab_idx or cm_key not in cm_idx:
        return 0.0

    return float(U[lab_idx[lab_key], cm_idx[cm_key]])


@lru_cache(maxsize=256)
def _moshinsky_block(
    N: int, L: int,
) -> Optional[Tuple[Dict, Dict, np.ndarray]]:
    """
    Build the Moshinsky bracket matrix for a given (N, L) block using
    Cartesian factorization.
    """
    # Enumerate lab-frame states
    lab_list: List[Tuple[int, int, int, int]] = []
    for N1 in range(N + 1):
        N2 = N - N1
        for l1 in range(N1, -1, -2):
            n1 = (N1 - l1) // 2
            for l2 in range(N2, -1, -2):
                n2 = (N2 - l2) // 2
                if abs(l1 - l2) <= L <= l1 + l2:
                    key = (n1, l1, n2, l2)
                    if key not in lab_list:
                        lab_list.append(key)

    # Enumerate CM+rel states
    cm_list: List[Tuple[int, int, int, int]] = []
    for N_cm in range(N + 1):
        N_rel = N - N_cm
        for l_cm in range(N_cm, -1, -2):
            n_cm = (N_cm - l_cm) // 2
            for l_rel in range(N_rel, -1, -2):
                n_rel = (N_rel - l_rel) // 2
                if abs(l_cm - l_rel) <= L <= l_cm + l_rel:
                    key = (n_cm, l_cm, n_rel, l_rel)
                    if key not in cm_list:
                        cm_list.append(key)

    dim = len(lab_list)
    if dim == 0:
        return None

    # Compute brackets via Cartesian factorization
    M = 0  # Use M=0 (brackets are M-independent)

    U = np.zeros((dim, dim))
    for i, (n1, l1, n2, l2) in enumerate(lab_list):
        for j, (nc, lc, nr, lr) in enumerate(cm_list):
            U[i, j] = _bracket_cartesian(n1, l1, n2, l2, nc, lc, nr, lr, L, M)

    # Enforce unitarity via SVD
    if dim > 1:
        u_svd, s_svd, vt_svd = np.linalg.svd(U, full_matrices=False)
        U = u_svd @ vt_svd
    elif dim == 1:
        # 1x1: just sign
        if abs(U[0, 0]) > 0.5:
            U[0, 0] = np.sign(U[0, 0])

    lab_idx = {key: i for i, key in enumerate(lab_list)}
    cm_idx = {key: j for j, key in enumerate(cm_list)}

    return lab_idx, cm_idx, U


def _bracket_cartesian(
    n1: int, l1: int, n2: int, l2: int,
    nc: int, lc: int, nr: int, lr: int,
    L: int, M: int = 0,
) -> float:
    """
    Compute 3D Moshinsky bracket via Cartesian expansion.

    1. Expand lab state |n1 l1 m1, n2 l2 m2; L M> in Cartesian basis
    2. Expand CM+rel state |nc lc mc, nr lr mr; L M> in Cartesian basis
    3. Apply 1D Moshinsky brackets to each Cartesian direction
    4. Sum up the overlap
    """
    # Get Cartesian expansions for the coupled states
    lab_cart = _coupled_cartesian_expansion(n1, l1, n2, l2, L, M)
    cm_cart = _coupled_cartesian_expansion(nc, lc, nr, lr, L, M)

    if not lab_cart or not cm_cart:
        return 0.0

    bracket = 0.0 + 0.0j
    for (n1x, n1y, n1z, n2x, n2y, n2z), c_lab in lab_cart:
        for (ncx, ncy, ncz, nrx, nry, nrz), c_cm in cm_cart:
            # 1D brackets
            bx = _moshinsky_1d(n1x, n2x, ncx, nrx)
            if abs(bx) < 1e-15:
                continue
            by = _moshinsky_1d(n1y, n2y, ncy, nry)
            if abs(by) < 1e-15:
                continue
            bz = _moshinsky_1d(n1z, n2z, ncz, nrz)
            if abs(bz) < 1e-15:
                continue

            bracket += np.conj(c_lab) * c_cm * bx * by * bz

    return float(bracket.real)


def _coupled_cartesian_expansion(
    n_r1: int, l1: int,
    n_r2: int, l2: int,
    L: int, M: int,
) -> List[Tuple[Tuple[int, int, int, int, int, int], complex]]:
    """
    Expand |n_r1 l1, n_r2 l2; L M> in two-particle Cartesian basis.

    Returns list of ((n1x,n1y,n1z,n2x,n2y,n2z), coefficient).
    """
    result: Dict[Tuple[int,...], complex] = {}

    for m1 in range(-l1, l1 + 1):
        for m2 in range(-l2, l2 + 1):
            if m1 + m2 != M:
                continue
            cg = _clebsch_gordan_int(l1, m1, l2, m2, L, M)
            if abs(cg) < 1e-15:
                continue

            keys1, coeffs1 = _spherical_to_cartesian(n_r1, l1, m1)
            keys2, coeffs2 = _spherical_to_cartesian(n_r2, l2, m2)

            for (n1x, n1y, n1z), c1 in zip(keys1, coeffs1):
                for (n2x, n2y, n2z), c2 in zip(keys2, coeffs2):
                    key = (n1x, n1y, n1z, n2x, n2y, n2z)
                    result[key] = result.get(key, 0.0+0.0j) + cg * c1 * c2

    return [(k, v) for k, v in result.items() if abs(v) > 1e-14]


# ---------------------------------------------------------------------------
# Lab-to-relative transformation of two-body matrix elements
# ---------------------------------------------------------------------------

def lab_to_relative_matrix_element(
    n1: int, l1: int,
    n2: int, l2: int,
    n3: int, l3: int,
    n4: int, l4: int,
    L: int, S: int,
    V_rel_func,
    b: float = 1.0,
) -> float:
    """
    Compute <n1 l1, n2 l2; LS | V | n3 l3, n4 l4; LS> by Moshinsky transformation.

    1. Transform bra and ket to CM+rel coordinates via Moshinsky brackets.
    2. V acts only on relative coordinate -> diagonal in CM quantum numbers.
    3. Sum over intermediate CM+rel states.

    Parameters
    ----------
    n1, l1, n2, l2 : int
        Bra lab-frame quantum numbers.
    n3, l3, n4, l4 : int
        Ket lab-frame quantum numbers.
    L : int
        Total orbital angular momentum.
    S : int
        Total spin.
    V_rel_func : callable
        V_rel_func(n_rel, l_rel, n_rel_prime, l_rel_prime, S, b) -> float
    b : float
        HO length parameter.

    Returns
    -------
    float
        Two-body matrix element.
    """
    N_bra = 2*n1 + l1 + 2*n2 + l2
    N_ket = 2*n3 + l3 + 2*n4 + l4

    if N_bra != N_ket:
        return 0.0

    N = N_bra

    # Enumerate all CM+rel states for this N and L
    cm_rel_states = []
    for N_cm in range(N + 1):
        N_rel = N - N_cm
        for l_cm in range(N_cm, -1, -2):
            n_cm = (N_cm - l_cm) // 2
            for l_rel in range(N_rel, -1, -2):
                n_rel = (N_rel - l_rel) // 2
                if abs(l_cm - l_rel) <= L <= l_cm + l_rel:
                    cm_rel_states.append((n_cm, l_cm, n_rel, l_rel))

    result = 0.0
    for n_cm_a, l_cm_a, n_rel_a, l_rel_a in cm_rel_states:
        bra_bracket = moshinsky_bracket(n1, l1, n2, l2,
                                         n_cm_a, l_cm_a, n_rel_a, l_rel_a, L)
        if abs(bra_bracket) < 1e-15:
            continue

        for n_cm_b, l_cm_b, n_rel_b, l_rel_b in cm_rel_states:
            if n_cm_a != n_cm_b or l_cm_a != l_cm_b:
                continue

            ket_bracket = moshinsky_bracket(n3, l3, n4, l4,
                                             n_cm_b, l_cm_b, n_rel_b, l_rel_b, L)
            if abs(ket_bracket) < 1e-15:
                continue

            v_rel = V_rel_func(n_rel_b, l_rel_b, n_rel_a, l_rel_a, S, b)
            result += bra_bracket * ket_bracket * v_rel

    return result
