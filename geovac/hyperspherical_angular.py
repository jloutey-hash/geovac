"""
Hyperspherical angular eigenvalue solver for two-electron atoms.

At each fixed hyperradius R, solves the coupled angular eigenvalue problem
in the L=0 sector (He ground state) using FD discretization in alpha
and partial-wave expansion in theta_12.

The equation solved is:
    [Lambda^2/2 + R*C(alpha, theta_12)] Phi = mu(R) Phi

where C is the charge function and mu(R) relates to the adiabatic
potential via V_eff(R) = mu(R)/R^2 + 15/(8R^2).

Uses the Liouville substitution u_l(alpha) = sin(alpha)*cos(alpha)*f_l(alpha)
to transform the weighted Sturm-Liouville problem into a standard
Schrodinger equation with Dirichlet BCs u(0) = u(pi/2) = 0.

References:
  - Macek, J. Phys. B 1, 831 (1968)
  - Lin, Phys. Rep. 257, 1 (1995)
"""

import numpy as np
from scipy.linalg import eigh
from typing import Tuple
from math import factorial, sqrt


def gaunt_integral(l1: int, k: int, l2: int) -> float:
    """
    Compute the Gaunt integral: integral of P_l1(x) P_k(x) P_l2(x) dx
    over [-1, 1].

    Uses the Wigner 3j symbol relation.

    Parameters
    ----------
    l1, k, l2 : int
        Legendre polynomial orders.

    Returns
    -------
    float
        Value of the Gaunt integral.
    """
    s = l1 + k + l2
    if s % 2 != 0:
        return 0.0
    if l2 > l1 + k or l2 < abs(l1 - k):
        return 0.0

    g = s // 2
    if g < l1 or g < k or g < l2:
        return 0.0

    num = factorial(2 * (g - l1)) * factorial(2 * (g - k)) * factorial(2 * (g - l2))
    den = factorial(2 * g + 1)
    threej_sq = (factorial(g) ** 2 * num) / (
        factorial(g - l1) ** 2 * factorial(g - k) ** 2
        * factorial(g - l2) ** 2 * den
    )
    return 2.0 * threej_sq


def _precompute_gaunt(l_max: int) -> np.ndarray:
    """Precompute Gaunt integrals for all (l, k, l') up to l_max."""
    n_l = l_max + 1
    k_max = 2 * l_max
    G = np.zeros((n_l, k_max + 1, n_l))
    for l in range(n_l):
        for lp in range(n_l):
            for k in range(k_max + 1):
                G[l, k, lp] = gaunt_integral(l, k, lp)
    return G


def solve_angular(
    R: float,
    Z: float = 2.0,
    l_max: int = 3,
    n_alpha: int = 100,
    n_channels: int = 1,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Solve the hyperangular eigenvalue problem at fixed hyperradius R.

    Uses the Liouville substitution u_l = sin(a)*cos(a)*f_l to get:
        -1/2 u_l'' + V_l(alpha) u_l + sum_{l'} W_{ll'}(alpha) u_{l'} = mu u_l

    where V_l = -2 + l(l+1)(1/cos^2 + 1/sin^2)/2 + R(-Z/cos - Z/sin)
    and W_{ll'} encodes the V_ee multipole coupling.

    BC: u(0) = u(pi/2) = 0 (natural from u = sin*cos*f).

    Parameters
    ----------
    R : float
        Hyperradius (bohr).
    Z : float
        Nuclear charge.
    l_max : int
        Maximum partial wave (l1=l2=l for L=0).
    n_alpha : int
        Number of FD grid points in alpha.
    n_channels : int
        Number of eigenvalues to return.

    Returns
    -------
    mu : ndarray of shape (n_channels,)
        Eigenvalues mu(R) of the angular equation.
    vecs : ndarray of shape (n_channels, N)
        Eigenvectors in the FD basis.
    """
    n_l = l_max + 1
    N = n_l * n_alpha

    # --- Alpha FD grid (interior points, Dirichlet at boundaries) ---
    h = (np.pi / 2) / (n_alpha + 1)
    alpha = (np.arange(n_alpha) + 1) * h  # in (0, pi/2)

    sin_a = np.sin(alpha)
    cos_a = np.cos(alpha)

    # --- Build Hamiltonian in the u-basis ---
    # -1/2 u'' + V u = mu u (standard Schrodinger form)
    # FD: -1/2 (u_{i+1} - 2u_i + u_{i-1})/h^2 + V_i u_i = mu u_i
    # Tridiagonal: diag = 1/h^2 + V_i, off = -1/(2h^2)

    H = np.zeros((N, N))

    def idx(l: int, i: int) -> int:
        return l * n_alpha + i

    kinetic_diag = 1.0 / h**2
    kinetic_off = -0.5 / h**2

    # Precompute Gaunt integrals
    G = _precompute_gaunt(l_max)

    # V_ee coupling coefficients
    # C_ee(alpha) = 1/sqrt(1-sin2a*cos(theta)) = sum_k (min/max)^k/max P_k
    # In the R*C equation, this is multiplied by R.
    min_sc = np.minimum(sin_a, cos_a)
    max_sc = np.maximum(sin_a, cos_a)

    for l in range(n_l):
        # Potential for channel l (from Liouville substitution):
        # V_l(a) = -2 + l(l+1)(1/cos^2(a) + 1/sin^2(a))/2 + R(-Z/cos(a) - Z/sin(a))
        V_l = (-2.0
               + 0.5 * l * (l + 1) * (1.0 / cos_a**2 + 1.0 / sin_a**2)
               + R * (-Z / cos_a - Z / sin_a))

        # Diagonal: kinetic + potential
        for i in range(n_alpha):
            ii = idx(l, i)
            H[ii, ii] = kinetic_diag + V_l[i]

        # Off-diagonal: kinetic (tridiagonal within l-channel)
        for i in range(n_alpha - 1):
            ii = idx(l, i)
            jj = idx(l, i + 1)
            H[ii, jj] = kinetic_off
            H[jj, ii] = kinetic_off

    # V_ee coupling between channels (point-wise in alpha)
    for l in range(n_l):
        for lp in range(l, n_l):
            # Compute W_{ll'}(alpha) = norm * sum_k gaunt(l,k,lp) * f_k / max_sc
            W = np.zeros(n_alpha)
            for k in range(abs(l - lp), l + lp + 1):
                if k > 2 * l_max:
                    continue
                g_val = G[l, k, lp]
                if abs(g_val) < 1e-15:
                    continue
                f_k = (min_sc / max_sc) ** k
                W += g_val * f_k / max_sc

            # Normalization for orthonormal Legendre basis
            norm = sqrt((2 * l + 1) * (2 * lp + 1)) / 2.0

            # Add to Hamiltonian (diagonal in alpha grid index)
            # R*C_ee = R * sum_k (min/max)^k/max P_k
            for i in range(n_alpha):
                ii = idx(l, i)
                jj = idx(lp, i)
                val = R * norm * W[i]
                H[ii, jj] += val
                if l != lp:
                    H[jj, ii] += val

    # --- Diagonalize ---
    evals, evecs = eigh(H)

    return evals[:n_channels], evecs[:, :n_channels].T
