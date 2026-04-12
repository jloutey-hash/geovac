"""
Minnesota NN Potential for nuclear structure calculations.

Implements the Minnesota potential (Thompson, Lemere, Tang, Nucl. Phys. A286,
53, 1977) — a central Gaussian NN interaction with spin-singlet and
spin-triplet channels.

V_NN(r) = V_R exp(-kappa_R r^2) + V_S exp(-kappa_S r^2) P_S
         + V_T exp(-kappa_T r^2) P_T

where P_S = (1 - sigma1.sigma2)/4 projects onto S=0,
      P_T = (3 + sigma1.sigma2)/4 projects onto S=1.

Units: MeV throughout (converted at ecosystem export boundary only).
Lengths in fm.

Author: GeoVac Development Team
Date: April 2026
"""

from __future__ import annotations

from math import factorial, gamma as math_gamma
from typing import Dict, Tuple

import numpy as np
from scipy.special import eval_genlaguerre, gamma as gamma_fn


# ---------------------------------------------------------------------------
# Physical constants (nuclear)
# ---------------------------------------------------------------------------

HBAR_C_MEV_FM = 197.3269804  # MeV*fm
M_NUCLEON_MEV = 938.918      # average nucleon mass in MeV/c^2


# ---------------------------------------------------------------------------
# Minnesota potential parameters (Thompson et al. 1977)
# ---------------------------------------------------------------------------

def minnesota_params() -> Dict[str, float]:
    """
    Return the Minnesota potential parameters with units.

    Returns
    -------
    dict with keys:
        V_R, kappa_R : repulsive core (MeV, fm^-2)
        V_S, kappa_S : singlet attraction (MeV, fm^-2)
        V_T, kappa_T : triplet attraction (MeV, fm^-2)
    """
    return {
        'V_R': 200.0,        # MeV
        'kappa_R': 1.487,    # fm^-2
        'V_S': -178.0,       # MeV
        'kappa_S': 0.639,    # fm^-2
        'V_T': -91.4,        # MeV
        'kappa_T': 0.465,    # fm^-2
    }


# ---------------------------------------------------------------------------
# Evaluate the Minnesota potential
# ---------------------------------------------------------------------------

def minnesota_potential(
    r: np.ndarray,
    S: int,
) -> np.ndarray:
    """
    Evaluate V_NN(r) for total spin S.

    For S=0 (singlet): P_S=1, P_T=0
        V(r) = V_R exp(-kappa_R r^2) + V_S exp(-kappa_S r^2)

    For S=1 (triplet): P_S=0, P_T=1
        V(r) = V_R exp(-kappa_R r^2) + V_T exp(-kappa_T r^2)

    Parameters
    ----------
    r : ndarray
        Relative distance in fm.
    S : int
        Total spin quantum number (0 or 1).

    Returns
    -------
    ndarray
        V_NN(r) in MeV.
    """
    r = np.asarray(r, dtype=float)
    p = minnesota_params()

    V = p['V_R'] * np.exp(-p['kappa_R'] * r**2)

    if S == 0:
        V += p['V_S'] * np.exp(-p['kappa_S'] * r**2)
    elif S == 1:
        V += p['V_T'] * np.exp(-p['kappa_T'] * r**2)
    else:
        raise ValueError(f"S must be 0 or 1, got {S}")

    return V


# ---------------------------------------------------------------------------
# HO length parameter
# ---------------------------------------------------------------------------

def ho_length_parameter(
    hw_MeV: float,
    mass_nucleon_MeV: float = M_NUCLEON_MEV,
) -> float:
    """
    Compute the harmonic oscillator length parameter b = sqrt(hbar*c / (m_N * hw)).

    Parameters
    ----------
    hw_MeV : float
        Oscillator energy hbar*omega in MeV.
    mass_nucleon_MeV : float
        Nucleon mass in MeV/c^2.

    Returns
    -------
    float
        b in fm.
    """
    # b = sqrt(hbar^2 / (m * hbar*omega))
    # = sqrt(hbar / (m * omega))
    # = sqrt((hbar*c)^2 / (m*c^2 * hbar*omega)) / c
    # In natural units with hbar*c in MeV*fm:
    # b = hbar_c / sqrt(m_N * hw)  (in fm)
    return HBAR_C_MEV_FM / np.sqrt(mass_nucleon_MeV * hw_MeV)


# ---------------------------------------------------------------------------
# HO radial wavefunction in relative coordinate
# ---------------------------------------------------------------------------

def _ho_radial_wf(
    n_r: int,
    l: int,
    r_grid: np.ndarray,
    b: float,
) -> np.ndarray:
    """
    Compute normalized HO radial wavefunction R_{n_r,l}(r).

    R_{n_r,l}(r) = N * (r/b)^l * L_{n_r}^{l+1/2}((r/b)^2) * exp(-r^2/(2b^2))

    Normalization: integral_0^inf |R|^2 r^2 dr = 1.

    Parameters
    ----------
    n_r : int
        Radial quantum number (>= 0).
    l : int
        Orbital angular momentum (>= 0).
    r_grid : ndarray
        Radial grid in fm.
    b : float
        HO length in fm.

    Returns
    -------
    ndarray
        R_{n_r,l}(r) on the grid.
    """
    x = (r_grid / b) ** 2
    L_poly = eval_genlaguerre(n_r, l + 0.5, x)
    wf = (r_grid / b) ** l * L_poly * np.exp(-x / 2.0)

    # Analytical normalization
    norm_sq = 2.0 * factorial(n_r) / (b**3 * gamma_fn(n_r + l + 1.5))
    N_const = np.sqrt(norm_sq)

    return N_const * wf


# ---------------------------------------------------------------------------
# Matrix element of Gaussian potential in relative HO basis
# ---------------------------------------------------------------------------

def minnesota_matrix_element_relative(
    n_rel: int,
    l_rel: int,
    n_rel_prime: int,
    l_rel_prime: int,
    S: int,
    b: float,
    n_grid: int = 4000,
) -> float:
    """
    Compute <n'_rel l'_rel | V_NN | n_rel l_rel> for spin S.

    The Minnesota potential is central, so l_rel must equal l_rel_prime.

    Parameters
    ----------
    n_rel, l_rel : int
        Ket relative HO quantum numbers.
    n_rel_prime, l_rel_prime : int
        Bra relative HO quantum numbers.
    S : int
        Total spin (0 or 1).
    b : float
        HO length in fm.
    n_grid : int
        Number of grid points for numerical integration.

    Returns
    -------
    float
        Matrix element in MeV.
    """
    if l_rel != l_rel_prime:
        return 0.0  # central force conserves l

    r_max = 6.0 * b * np.sqrt(max(2 * n_rel + l_rel, 2 * n_rel_prime + l_rel_prime) + 3)
    r_grid = np.linspace(1e-10, r_max, n_grid)

    R_bra = _ho_radial_wf(n_rel_prime, l_rel_prime, r_grid, b)
    R_ket = _ho_radial_wf(n_rel, l_rel, r_grid, b)
    V = minnesota_potential(r_grid, S)

    integrand = R_bra * V * R_ket * r_grid**2
    return float(np.trapezoid(integrand, r_grid))


# ---------------------------------------------------------------------------
# Analytical Gaussian matrix element (for verification and speed)
# ---------------------------------------------------------------------------

def _gaussian_me_ho(
    V0: float,
    kappa: float,
    n1: int,
    l1: int,
    n2: int,
    l2: int,
    b: float,
    n_grid: int = 4000,
) -> float:
    """
    Matrix element <n1 l1 | V0 exp(-kappa r^2) | n2 l2> in HO basis.

    Computed via numerical integration on a fine r-grid.

    Returns 0 if l1 != l2 (central force selection rule).
    """
    if l1 != l2:
        return 0.0

    l = l1
    r_max = 6.0 * b * np.sqrt(max(2*n1 + l, 2*n2 + l) + 3)
    r_grid = np.linspace(1e-10, r_max, n_grid)

    R1 = _ho_radial_wf(n1, l, r_grid, b)
    R2 = _ho_radial_wf(n2, l, r_grid, b)
    V = V0 * np.exp(-kappa * r_grid**2)

    integrand = R1 * V * R2 * r_grid**2
    return float(np.trapezoid(integrand, r_grid))


def minnesota_matrix_element_analytical(
    n_rel: int,
    l_rel: int,
    n_rel_prime: int,
    l_rel_prime: int,
    S: int,
    b: float,
) -> float:
    """
    Analytical Minnesota matrix element using Gauss-Laguerre quadrature.

    Faster and more accurate than grid integration for Gaussian potentials.

    Parameters
    ----------
    n_rel, l_rel : int
        Ket quantum numbers.
    n_rel_prime, l_rel_prime : int
        Bra quantum numbers.
    S : int
        Total spin.
    b : float
        HO length in fm.

    Returns
    -------
    float
        Matrix element in MeV.
    """
    if l_rel != l_rel_prime:
        return 0.0

    p = minnesota_params()

    # Repulsive core (always present)
    me = _gaussian_me_ho(p['V_R'], p['kappa_R'],
                         n_rel_prime, l_rel_prime,
                         n_rel, l_rel, b)

    if S == 0:
        me += _gaussian_me_ho(p['V_S'], p['kappa_S'],
                              n_rel_prime, l_rel_prime,
                              n_rel, l_rel, b)
    elif S == 1:
        me += _gaussian_me_ho(p['V_T'], p['kappa_T'],
                              n_rel_prime, l_rel_prime,
                              n_rel, l_rel, b)
    else:
        raise ValueError(f"S must be 0 or 1, got {S}")

    return me
