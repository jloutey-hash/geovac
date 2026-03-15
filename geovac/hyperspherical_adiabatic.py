"""
Adiabatic potential curves for the hyperspherical helium solver.

Computes mu(R) on a grid of hyperradius values. The angular solver
returns mu(R), eigenvalue of [Lambda^2/2 + R*C(Omega)].

The effective potential for the hyperradial equation is:
    V_eff(R) = mu(R)/R^2 + 15/(8R^2)

Asymptotically: mu(R) -> -Z^2 R^2 / 2, so V_eff -> -Z^2/2 (He+ threshold).

References:
  - Macek, J. Phys. B 1, 831 (1968)
  - Lin, Phys. Rep. 257, 1 (1995)
"""

import numpy as np
from typing import Optional

from geovac.hyperspherical_angular import solve_angular


def compute_adiabatic_curve(
    R_grid: np.ndarray,
    Z: float = 2.0,
    l_max: int = 3,
    n_alpha: int = 100,
    n_channels: int = 1,
) -> np.ndarray:
    """
    Compute adiabatic eigenvalues mu(R) on a grid.

    Parameters
    ----------
    R_grid : ndarray of shape (N_R,)
        Hyperradius values (bohr).
    Z : float
        Nuclear charge.
    l_max : int
        Maximum partial wave.
    n_alpha : int
        Number of FD grid points for alpha.
    n_channels : int
        Number of adiabatic channels to compute.

    Returns
    -------
    mu : ndarray of shape (n_channels, N_R)
        Angular eigenvalues mu(R) for each channel.
    """
    N_R = len(R_grid)
    mu = np.zeros((n_channels, N_R))

    for i, R in enumerate(R_grid):
        evals, _ = solve_angular(R, Z, l_max, n_alpha, n_channels)
        mu[:, i] = evals

    return mu


def effective_potential(
    R_grid: np.ndarray,
    mu_curve: np.ndarray,
) -> np.ndarray:
    """
    Compute the effective hyperradial potential.

    V_eff(R) = mu(R)/R^2 + 15/(8R^2)

    where mu(R) is the eigenvalue of [Lambda^2/2 + R*C(Omega)] and
    15/(8R^2) is the Jacobian centrifugal term from extracting R^{-5/2}.

    Parameters
    ----------
    R_grid : ndarray of shape (N_R,)
        Hyperradius values.
    mu_curve : ndarray of shape (N_R,)
        Angular eigenvalue curve (single channel).

    Returns
    -------
    V_eff : ndarray of shape (N_R,)
        Effective potential for the hyperradial equation.
    """
    return mu_curve / R_grid**2 + 15.0 / (8.0 * R_grid**2)


def plot_adiabatic_curves(
    R_grid: np.ndarray,
    mu: np.ndarray,
    Z: float = 2.0,
    save_path: Optional[str] = None,
) -> None:
    """
    Plot effective potential curves V_eff(R) = mu(R)/R^2 + 15/(8R^2).

    Parameters
    ----------
    R_grid : ndarray
        Hyperradius grid.
    mu : ndarray of shape (n_channels, N_R)
        Angular eigenvalues.
    Z : float
        Nuclear charge (for labeling thresholds).
    save_path : str, optional
        Path to save the plot.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    n_ch = mu.shape[0]
    for ch in range(n_ch):
        V_eff = effective_potential(R_grid, mu[ch])
        ax.plot(R_grid, V_eff, label=f'Channel {ch+1}')

    # He+ threshold
    ax.axhline(-Z**2 / 2, color='gray', linestyle='--', alpha=0.5,
               label=f'He+ threshold ({-Z**2/2:.1f} Ha)')

    ax.set_xlabel('Hyperradius R (bohr)')
    ax.set_ylabel('V_eff(R) (Ha)')
    ax.set_title(f'Adiabatic Potential Curves for Z={Z}')
    ax.set_xlim(0, R_grid[-1])
    ax.set_ylim(-5, 2)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"Saved plot to {save_path}")
    plt.close()
