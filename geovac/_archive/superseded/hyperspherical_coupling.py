"""
Non-adiabatic coupling matrix elements for the hyperspherical method.

Computes the first-derivative coupling P_μν(R) = ⟨Φ_μ|∂/∂R|Φ_ν⟩ using
the Hellmann-Feynman approach:

    P_μν(R) = ⟨Φ_μ|∂H_ang/∂R|Φ_ν⟩ / (μ_μ - μ_ν)   for μ ≠ ν

This avoids noisy numerical derivatives of eigenvectors and gives exact
coupling at each R point.

The diagonal Born-Oppenheimer correction (DBOC) is computed as:
    DBOC_μ(R) = ½ Σ_{ν≠μ} |P_μν|² = ½ Σ_{ν≠μ} |⟨Φ_μ|∂H/∂R|Φ_ν⟩|² / (μ_μ - μ_ν)²

References:
  - Macek, J. Phys. B 1, 831 (1968)
  - Lin, Phys. Rep. 257, 1 (1995)
  - Paper 13 (fiber bundle structure)
"""

import numpy as np
from math import sqrt
from typing import Tuple, Optional

from geovac.hyperspherical_angular import solve_angular, gaunt_integral, _precompute_gaunt


def _enforce_sign_consistency(
    vecs_prev: np.ndarray,
    vecs_curr: np.ndarray,
) -> np.ndarray:
    """
    Enforce sign consistency of eigenvectors between adjacent R points.

    Eigenvectors from eigh have arbitrary global sign. We fix this by
    requiring ⟨Φ_μ(R_i)|Φ_μ(R_{i-1})⟩ > 0 for each channel μ.

    Parameters
    ----------
    vecs_prev : ndarray of shape (n_channels, N)
        Eigenvectors at previous R point (already sign-fixed).
    vecs_curr : ndarray of shape (n_channels, N)
        Eigenvectors at current R point (to be sign-fixed).

    Returns
    -------
    vecs_fixed : ndarray of shape (n_channels, N)
        Sign-corrected eigenvectors.
    """
    n_ch = vecs_curr.shape[0]
    vecs_fixed = vecs_curr.copy()
    for mu in range(n_ch):
        overlap = np.dot(vecs_prev[mu], vecs_curr[mu])
        if overlap < 0:
            vecs_fixed[mu] *= -1
    return vecs_fixed


def _compute_dH_dR_matrix(
    R: float,
    Z: float,
    l_max: int,
    n_alpha: int,
) -> np.ndarray:
    """
    Compute ∂H_ang/∂R — the R-derivative of the angular Hamiltonian.

    H_ang = [kinetic + centrifugal] + R * C(Ω)

    where C(Ω) contains nuclear attraction (-Z/cos(α) - Z/sin(α)) and
    V_ee multipole coupling. Since kinetic and centrifugal terms are
    R-independent, ∂H/∂R = C(Ω).

    Parameters
    ----------
    R : float
        Hyperradius.
    Z : float
        Nuclear charge.
    l_max : int
        Maximum partial wave.
    n_alpha : int
        Angular FD grid points.

    Returns
    -------
    dHdR : ndarray of shape (N, N)
        R-derivative of the angular Hamiltonian matrix.
    """
    n_l = l_max + 1
    N = n_l * n_alpha

    h = (np.pi / 2) / (n_alpha + 1)
    alpha = (np.arange(n_alpha) + 1) * h

    sin_a = np.sin(alpha)
    cos_a = np.cos(alpha)

    dHdR = np.zeros((N, N))

    def idx(l: int, i: int) -> int:
        return l * n_alpha + i

    G = _precompute_gaunt(l_max)

    min_sc = np.minimum(sin_a, cos_a)
    max_sc = np.maximum(sin_a, cos_a)

    # Nuclear attraction: diagonal in l, point-wise in alpha
    for l in range(n_l):
        V_nuc = -Z / cos_a - Z / sin_a
        for i in range(n_alpha):
            ii = idx(l, i)
            dHdR[ii, ii] = V_nuc[i]

    # V_ee coupling (same structure as in solve_angular, but without R factor)
    for l in range(n_l):
        for lp in range(l, n_l):
            W = np.zeros(n_alpha)
            for k in range(abs(l - lp), l + lp + 1):
                if k > 2 * l_max:
                    continue
                g_val = G[l, k, lp]
                if abs(g_val) < 1e-15:
                    continue
                f_k = (min_sc / max_sc) ** k
                W += g_val * f_k / max_sc

            norm = sqrt((2 * l + 1) * (2 * lp + 1)) / 2.0

            for i in range(n_alpha):
                ii = idx(l, i)
                jj = idx(lp, i)
                val = norm * W[i]
                dHdR[ii, jj] += val
                if l != lp:
                    dHdR[jj, ii] += val

    return dHdR


def compute_coupling_matrices(
    R_grid: np.ndarray,
    Z: float = 2.0,
    l_max: int = 3,
    n_alpha: int = 100,
    n_channels: int = 3,
    return_vecs: bool = False,
) -> dict:
    """
    Compute non-adiabatic coupling P_μν(R) and DBOC on an R grid.

    Uses the Hellmann-Feynman approach:
        P_μν(R) = ⟨Φ_μ|∂H/∂R|Φ_ν⟩ / (μ_μ - μ_ν)   for μ ≠ ν

    This avoids numerical differentiation of eigenvectors.

    Parameters
    ----------
    R_grid : ndarray of shape (N_R,)
        Hyperradius grid.
    Z : float
        Nuclear charge.
    l_max : int
        Maximum partial wave.
    n_alpha : int
        Angular grid points.
    n_channels : int
        Number of adiabatic channels.
    return_vecs : bool
        If True, also return sign-consistent eigenvectors.

    Returns
    -------
    result : dict
        Keys:
        - 'P': ndarray (n_channels, n_channels, N_R) — first-derivative coupling
        - 'DBOC': ndarray (n_channels, N_R) — diagonal BO correction per channel
        - 'mu': ndarray (n_channels, N_R) — adiabatic eigenvalues
        - 'R_grid': the input R grid
    """
    N_R = len(R_grid)
    n_ch = n_channels

    mu = np.zeros((n_ch, N_R))
    P = np.zeros((n_ch, n_ch, N_R))
    DBOC = np.zeros((n_ch, N_R))

    all_vecs = []

    for i, R in enumerate(R_grid):
        evals, vecs = solve_angular(R, Z, l_max, n_alpha, n_channels)
        mu[:, i] = evals

        if i > 0:
            vecs = _enforce_sign_consistency(all_vecs[i - 1], vecs)
        all_vecs.append(vecs)

        # Compute dH/dR matrix
        dHdR = _compute_dH_dR_matrix(R, Z, l_max, n_alpha)

        # Hellmann-Feynman coupling
        for mu_idx in range(n_ch):
            for nu_idx in range(n_ch):
                if mu_idx == nu_idx:
                    P[mu_idx, nu_idx, i] = 0.0
                    continue

                gap = evals[mu_idx] - evals[nu_idx]
                if abs(gap) < 1e-10:
                    # Degenerate eigenvalues — coupling undefined, set to 0
                    P[mu_idx, nu_idx, i] = 0.0
                    continue

                # ⟨Φ_μ|dH/dR|Φ_ν⟩
                matrix_element = vecs[mu_idx] @ dHdR @ vecs[nu_idx]
                P[mu_idx, nu_idx, i] = matrix_element / gap

            # DBOC for channel mu: ½ Σ_{ν≠μ} P²_μν
            DBOC[mu_idx, i] = 0.5 * np.sum(P[mu_idx, :, i]**2)

    result = {
        'P': P,
        'DBOC': DBOC,
        'mu': mu,
        'R_grid': R_grid,
    }
    if return_vecs:
        result['vecs'] = all_vecs

    return result


def compute_coupling_on_grid(
    Z: float = 2.0,
    l_max: int = 3,
    n_alpha: int = 100,
    n_channels: int = 3,
    N_R: int = 200,
    R_max: float = 30.0,
) -> dict:
    """
    Compute coupling matrices on an adaptive R grid suitable for He.

    Denser near R = 1-4 bohr where avoided crossings produce large P_μν.

    Parameters
    ----------
    Z, l_max, n_alpha, n_channels : solver parameters
    N_R : int
        Approximate total number of R points.
    R_max : float
        Maximum hyperradius.

    Returns
    -------
    result : dict
        Same as compute_coupling_matrices.
    """
    R_grid = np.concatenate([
        np.linspace(0.1, 1.0, N_R // 3),
        np.linspace(1.0, 5.0, N_R // 3),
        np.linspace(5.0, R_max, N_R // 3 + 1),
    ])
    R_grid = np.unique(R_grid)

    return compute_coupling_matrices(
        R_grid, Z, l_max, n_alpha, n_channels
    )
