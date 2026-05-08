"""
Phillips-Kleinman Cross-Center Projection (W1b residual closure)
=================================================================

Implements Phillips-Kleinman-class projection on a valence orbital against
the frozen core of an OFF-CENTER nucleus. This is the W1b-residual
closure named in Phase C-W1c memo §6 (May 2026): the screened
cross-center V_ne (cross_center_screened_vne) reduces the bare cross-V_ne
attraction by 5-6x but leaves a residual ~0.36 Ha overattraction at NaH
because the H 1s orbital is non-orthogonal to the [Ne]-core orbitals on
Na, and Pauli exclusion (a kinetic-energy repulsion barrier in the strict
Schmidt-orthogonalization limit) is not enforced.

Mechanism: standard PK barrier in second-quantized matrix-element form:

    Delta H[p, q]^PK = sum_c (E_v - E_c) * S_pc * S_cq

where:
  - p, q index the valence orbitals on the off-center side (e.g., the H
    block in NaH)
  - c indexes the frozen-core orbitals on the heavy atom (e.g., 1s, 2s,
    2p_{x,y,z} of Na's [Ne] core)
  - E_c is the core orbital binding energy (Clementi-Roetti HF eigenvalue)
  - E_v is the valence reference energy (default 0 -> "absolute PK"
    barrier, purely repulsive)
  - S_pc is the cross-center overlap <phi_p^valence | phi_c^core>

S_pc is computed via 2D numerical quadrature in (r, u) where u = cos(theta)
in the valence-orbital-centered frame, with rho = sqrt(r^2 + R^2 - 2rRu)
the distance from the off-center core orbital. Azimuthal selection
(real-SH m must match for both orbitals) is exact.

The PK barrier is the standard operator approximation of strict Schmidt
orthogonalization (1 - P_core) in the limit |E_c| >> |E_v|. It is
repulsive whenever E_v - E_c > 0, exactly matching the kinetic-energy
repulsion that valence-on-core orthogonality enforces.

Note: this is the cross-center analog of the same-center PK already in
ab_initio_pk.py. The same-center PK barrier handles core electrons in the
SAME block as the valence (e.g., the explicit Li_core block of LiH); the
cross-center PK barrier handles frozen-core electrons in a DIFFERENT
block (e.g., the [Ne] core of Na sitting on a different atomic center).

Author: GeoVac Development Team (Phase C / Track 2 follow-up)
Date: May 2026
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.special import genlaguerre, lpmv

from geovac.cross_center_screened_vne import _detect_core_type
from geovac.shibuya_wulfman import _build_block_rotation_matrix


# ---------------------------------------------------------------------------
# Core orbital data: Clementi-Raimondi exponents and approximate eigenvalues
# ---------------------------------------------------------------------------

# Clementi-Roetti HF 1-electron eigenvalues for the orbitals in each frozen
# core (Hartree). Source: Clementi & Roetti, At. Data Nucl. Data Tables 14,
# 177 (1974), Table II (atomic Hartree-Fock orbital energies).
#
# These are the binding energies (-epsilon) used as the PK-barrier weight.
# A purely repulsive PK uses |E_c| as the weight (E_v = 0 reference).

# [Ne] core orbital energies (Hartree, neutral atom HF eigenvalues for the
# corresponding free atom). Used as approximate per-orbital eigenvalues;
# the (n,l)-resolved values are sufficient for PK weighting.
_CORE_ORBITAL_ENERGIES_NE = {
    # Z: {(n,l): epsilon}
    11: {(1, 0): -40.4787, (2, 0): -2.7967, (2, 1): -1.5181},   # Na (HF, neutral)
    12: {(1, 0): -49.0317, (2, 0): -3.7676, (2, 1): -2.2820},   # Mg
    13: {(1, 0): -58.5009, (2, 0): -4.9105, (2, 1): -3.2184},   # Al
    14: {(1, 0): -68.8123, (2, 0): -6.1567, (2, 1): -4.2557},   # Si
    15: {(1, 0): -79.9696, (2, 0): -7.5111, (2, 1): -5.4012},   # P
    16: {(1, 0): -92.0045, (2, 0): -8.9716, (2, 1): -6.6824},   # S
    17: {(1, 0): -104.8842, (2, 0): -10.6089, (2, 1): -8.0721}, # Cl
    18: {(1, 0): -118.6105, (2, 0): -12.3221, (2, 1): -9.5712}, # Ar
}


def _core_orbitals_for_Z(Z_nuc: int, core_type: str) -> List[Dict]:
    """Return list of frozen-core orbital descriptors for an off-center nucleus.

    Each descriptor is {'n', 'l', 'm', 'zeta', 'energy'} where:
      - n, l, m: quantum numbers (m in real-SH signed convention)
      - zeta: Clementi-Raimondi orbital exponent (so Z_eff = n * zeta for
        the hydrogenic radial wavefunction R_{nl}(r; Z_eff))
      - energy: orbital eigenvalue (Hartree, -|E_binding|)

    Currently supports [Ne] core only. [Ar], [Kr], [Xe] would extend
    similarly if needed.
    """
    Z_int = int(round(Z_nuc))

    if core_type != 'Ne':
        # Other cores: stub out with zero. PK barrier = 0 means the wall
        # is documented but not addressed for these systems in this sprint.
        return []

    if Z_int not in _CORE_ORBITAL_ENERGIES_NE:
        return []

    from geovac.neon_core import _CLEMENTI_ZETA_NE

    zeta_1s, zeta_2s, zeta_2p = _CLEMENTI_ZETA_NE[Z_int]
    energies = _CORE_ORBITAL_ENERGIES_NE[Z_int]

    orbs: List[Dict] = []
    # 1s (n=1, l=0, m=0)
    orbs.append({'n': 1, 'l': 0, 'm': 0,
                 'zeta': zeta_1s, 'energy': energies[(1, 0)]})
    # 2s (n=2, l=0, m=0)
    orbs.append({'n': 2, 'l': 0, 'm': 0,
                 'zeta': zeta_2s, 'energy': energies[(2, 0)]})
    # 2p (n=2, l=1, m=-1, 0, +1)
    for m in (-1, 0, 1):
        orbs.append({'n': 2, 'l': 1, 'm': m,
                     'zeta': zeta_2p, 'energy': energies[(2, 1)]})

    return orbs


# ---------------------------------------------------------------------------
# Hydrogenic radial wavefunction at arbitrary points (not on a fixed grid)
# ---------------------------------------------------------------------------


def _hydrogenic_radial_at_points(
    Z_eff: float, n: int, l: int, r_arr: np.ndarray,
) -> np.ndarray:
    """Hydrogenic R_{nl}(r; Z_eff) evaluated at arbitrary radial points.

    Returns the same normalization convention as
    geovac.composed_qubit._radial_wf_grid (normalized so that
    integral |R|^2 r^2 dr = 1 on a fine grid). For exact evaluation at
    arbitrary points, we use the analytical normalization constant for
    hydrogenic orbitals:

        R_{nl}(r) = N_{nl} * (2 Z_eff r / n)^l * exp(-Z_eff r / n) * L_{n-l-1}^{2l+1}(2 Z_eff r / n)

    with N_{nl} = sqrt[ (2 Z_eff / n)^3 * (n-l-1)! / (2n * (n+l)!) ]

    (Bethe-Salpeter, Eq. 3.18 conventions; consistent with
    composed_qubit._radial_wf_grid up to numerical normalization.)
    """
    from math import factorial

    rho = 2.0 * Z_eff * np.asarray(r_arr, dtype=float) / n
    L_poly = genlaguerre(n - l - 1, 2 * l + 1)(rho)
    radial = rho ** l * np.exp(-rho / 2.0) * L_poly

    # Analytical normalization (Bethe-Salpeter)
    norm = np.sqrt(
        (2.0 * Z_eff / n) ** 3
        * factorial(n - l - 1)
        / (2.0 * n * factorial(n + l))
    )
    return norm * radial


def _real_sh_normalization(l: int, m_signed: int) -> float:
    """Normalization N_{l,|m|} for real spherical harmonics in the convention

        Y_l^0(theta, phi) = N_{l,0} * P_l(cos theta)
        Y_l^|m|(theta, phi) = N_{l,|m|} * P_l^|m|(cos theta) * cos(|m| phi)  (m>0)
        Y_l^{-|m|}(theta, phi) = N_{l,|m|} * P_l^|m|(cos theta) * sin(|m| phi)  (m<0)

    where P_l^|m| follows scipy.special.lpmv convention (Condon-Shortley
    phase already absorbed: lpmv(m, l, x) = (-1)^m * P_l^m(x) standard).

        N_{l,0} = sqrt((2l+1) / (4 pi))
        N_{l,|m|} = sqrt(2 (2l+1) / (4 pi) * (l-|m|)! / (l+|m|)!)  (|m|>0)
    """
    from math import factorial
    abs_m = abs(m_signed)
    base = (2 * l + 1) / (4.0 * np.pi)
    if abs_m == 0:
        return float(np.sqrt(base))
    return float(np.sqrt(2.0 * base * factorial(l - abs_m) / factorial(l + abs_m)))


def _phi_integral_factor(m_p: int, m_c: int) -> float:
    """Integral over [0, 2pi] of the phi-dependent parts of two real SH.

    For real SH, the phi-dependence is:
      m=0: 1
      m>0: cos(m phi)
      m<0: sin(|m| phi)

    Phi-integral non-zero only when m_p == m_c (both signed). Returns:
      m_p = m_c = 0: 2 pi
      m_p = m_c > 0: pi (cos^2 integral)
      m_p = m_c < 0: pi (sin^2 integral)
      otherwise: 0
    """
    if m_p != m_c:
        return 0.0
    return 2.0 * np.pi if m_p == 0 else np.pi


# ---------------------------------------------------------------------------
# Cross-center overlap (the new computation)
# ---------------------------------------------------------------------------


def cross_center_overlap(
    Z_orb_p: float, n_p: int, l_p: int, m_p: int,
    Z_orb_c: float, n_c: int, l_c: int, m_c: int,
    R_AB: float,
    n_grid_r: int = 4000,
    n_grid_u: int = 64,
    r_max: Optional[float] = None,
) -> float:
    """Cross-center overlap <phi_p^A | phi_c^B> in z-axis frame.

    A is at origin, B is at +R_AB * z_hat. Both orbitals are hydrogenic
    with their respective effective charges. Real-SH convention with
    signed m index.

    Parameters
    ----------
    Z_orb_p : float
        Effective charge of the valence (A-side) orbital.
    n_p, l_p, m_p : int
        Quantum numbers of the valence orbital (real SH, signed m).
    Z_orb_c : float
        Effective charge of the core (B-side) orbital.
    n_c, l_c, m_c : int
        Quantum numbers of the core orbital.
    R_AB : float
        Internuclear distance (bohr), B at +R_AB * z_hat.
    n_grid_r : int
        Number of radial grid points.
    n_grid_u : int
        Number of u = cos(theta) Gauss-Legendre nodes.
    r_max : float, optional
        Maximum radial distance. Auto-chosen from orbital extent if None.

    Returns
    -------
    float
        Overlap value (dimensionless).
    """
    # Azimuthal selection: real-SH m_signed must match
    if m_p != m_c:
        return 0.0

    if r_max is None:
        # Choose to capture both orbital tails
        z_min = min(Z_orb_p, Z_orb_c)
        r_max = max(80.0 / max(z_min, 0.5), 4.0 * R_AB, 30.0)

    # Radial grid (np.linspace excluding 0)
    r = np.linspace(0, r_max, n_grid_r + 1)[1:]

    # Valence radial wavefunction R_p on r-grid
    R_p_r = _hydrogenic_radial_at_points(Z_orb_p, n_p, l_p, r)

    # Gauss-Legendre nodes/weights for u = cos(theta) on [-1, 1]
    u_nodes, u_weights = np.polynomial.legendre.leggauss(n_grid_u)

    # Vectorize over (r, u): rho(r, u), u_B(r, u)
    r_col = r[:, np.newaxis]
    u_row = u_nodes[np.newaxis, :]
    rho_2d = np.sqrt(r_col ** 2 + R_AB ** 2 - 2.0 * r_col * R_AB * u_row)
    rho_safe = np.where(rho_2d > 1e-12, rho_2d, 1e-12)
    u_B_2d = (r_col * u_row - R_AB) / rho_safe
    u_B_2d = np.clip(u_B_2d, -1.0, 1.0)

    # Core radial wavefunction R_c at the rho values (vectorized)
    rho_flat = rho_2d.reshape(-1)
    R_c_flat = _hydrogenic_radial_at_points(Z_orb_c, n_c, l_c, rho_flat)
    R_c_2d = R_c_flat.reshape(rho_2d.shape)

    # Associated Legendre P_l^|m|(u) for both orbitals
    abs_m = abs(m_p)
    P_lp_at_u = lpmv(abs_m, l_p, u_nodes)            # shape (n_u,)
    P_lc_at_uB = lpmv(abs_m, l_c, u_B_2d)            # shape (n_r, n_u)

    # Angular integrand at each (r, u)
    angular_2d = R_c_2d * P_lc_at_uB  # shape (n_r, n_u)
    weighted_2d = angular_2d * (u_weights * P_lp_at_u)[np.newaxis, :]

    u_integral = np.sum(weighted_2d, axis=1)  # shape (n_r,)

    # Radial integrand: r^2 * R_p(r) * u_integral(r)
    integrand_r = r ** 2 * R_p_r * u_integral
    radial_int = float(np.trapezoid(integrand_r, r))

    # Apply real-SH normalization and phi-integral factor
    norm_p = _real_sh_normalization(l_p, m_p)
    norm_c = _real_sh_normalization(l_c, m_c)
    phi_factor = _phi_integral_factor(m_p, m_c)

    return norm_p * norm_c * phi_factor * radial_int


# ---------------------------------------------------------------------------
# PK barrier matrix
# ---------------------------------------------------------------------------


def compute_pk_cross_center_barrier(
    Z_orb: float,
    valence_states: List[Tuple[int, int, int]],
    Z_nuc: float,
    R_AB: float,
    core_type: Optional[str] = None,
    E_valence_ref: float = 0.0,
    n_grid_r: int = 4000,
    n_grid_u: int = 64,
    direction: Optional[Tuple[float, float, float]] = None,
) -> np.ndarray:
    """Phillips-Kleinman barrier matrix Delta H[p, q] in the valence block.

    Delta H[p, q]^PK = sum_c (E_v - E_c) * S_pc * S_cq

    where p, q index valence_states (in the A-side orbital basis at Z_orb)
    and c indexes the frozen-core orbitals of the off-center nucleus B at
    distance R_AB. With E_v = 0 (default), Delta H is purely repulsive
    (since E_c < 0).

    Parameters
    ----------
    Z_orb : float
        Effective charge of the valence orbitals.
    valence_states : list of (n, l, m) tuples
        Real-SH labels of the valence orbitals (signed m).
    Z_nuc : float
        Bare nuclear charge of the off-center nucleus B.
    R_AB : float
        Internuclear distance (bohr).
    core_type : str, optional
        FrozenCore type for nucleus B. Auto-detected if None.
    E_valence_ref : float
        Valence reference energy E_v in (E_v - E_c). Default 0.0.
    n_grid_r : int
        Radial grid resolution for overlap.
    n_grid_u : int
        Gauss-Legendre nodes for angular projection.
    direction : tuple of (nx, ny, nz), optional
        Unit vector from valence-orbital center to off-center nucleus.
        When provided, the matrix is computed in z-frame and rotated.

    Returns
    -------
    np.ndarray
        Delta H matrix of shape (N_valence, N_valence). Returns zeros
        when no frozen core is detected for Z_nuc.
    """
    if core_type is None:
        core_type = _detect_core_type(Z_nuc)

    N = len(valence_states)
    if core_type is None:
        return np.zeros((N, N))

    core_orbs = _core_orbitals_for_Z(int(round(Z_nuc)), core_type)
    if not core_orbs:
        return np.zeros((N, N))

    # --- Direction-based: compute in z-frame, then rotate ---
    if direction is not None:
        delta_h_z = compute_pk_cross_center_barrier(
            Z_orb, valence_states, Z_nuc, R_AB,
            core_type=core_type, E_valence_ref=E_valence_ref,
            n_grid_r=n_grid_r, n_grid_u=n_grid_u, direction=None,
        )
        D = _build_block_rotation_matrix(valence_states, direction)
        return D @ delta_h_z @ D.T

    # --- z-axis frame ---
    # Build overlap matrix S[p, c] = <phi_p^A | phi_c^B>  (B at +R_AB z_hat)
    S = np.zeros((N, len(core_orbs)))
    for p, (n_p, l_p, m_p) in enumerate(valence_states):
        for c, core in enumerate(core_orbs):
            n_c, l_c, m_c = core['n'], core['l'], core['m']
            # Azimuthal selection: m must match (real SH)
            if m_p != m_c:
                continue
            zeta_c = core['zeta']
            # Z_eff for the core orbital wavefunction is n_c * zeta_c
            Z_eff_c = n_c * zeta_c
            S[p, c] = cross_center_overlap(
                Z_orb, n_p, l_p, m_p,
                Z_eff_c, n_c, l_c, m_c,
                R_AB,
                n_grid_r=n_grid_r, n_grid_u=n_grid_u,
            )

    # PK barrier: Delta H[p, q] = sum_c (E_v - E_c) * S_pc * S_cq
    energies = np.array([core['energy'] for core in core_orbs])
    weights = E_valence_ref - energies  # (E_v - E_c), positive for repulsion
    # Delta H = S @ diag(weights) @ S^T
    delta_h = S @ np.diag(weights) @ S.T
    return delta_h
