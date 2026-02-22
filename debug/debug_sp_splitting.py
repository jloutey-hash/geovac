"""
Debug: Isolate 2s/2p Degeneracy Splitting as Discretization Artifact
====================================================================

Proves that the ~16% ΔE(2p-2s) splitting observed in the graph Laplacian
Hamiltonian is a discretization artifact of the finite lattice, NOT an
emergent physical force.

Strategy:
---------
For each max_n from 5 to 30, build the GeometricLattice and corresponding
Hamiltonian H = (-1/16) * (D - A), extract eigenvalues corresponding to
the 2s and 2p states, and compute the relative splitting:

    ΔE_rel = |E(2p) - E(2s)| / |E(2s)|

If this is a topological artifact of finite graph truncation, the splitting
should decay as a power law in the total vertex count V:

    ΔE_rel ~ V^(-α)

with a predictable scaling exponent α > 0.

Output:
-------
- Console table of n_max, V, E(2s), E(2p), ΔE_rel
- Log-log plot of ΔE_rel vs V with power-law fit → debug/plots/sp_splitting.png

Author: GeoVac Development Team
Date: February 2026
Status: Theoretical critique response (Paper 1 s/p degeneracy)
"""

import sys
import os
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac import GeometricLattice, UNIVERSAL_KINETIC_SCALE


def build_hamiltonian(lattice: GeometricLattice, scale: float) -> 'csr_matrix':
    """Build H = scale * (D - A) from a GeometricLattice."""
    adj = lattice.adjacency
    degree = np.array(adj.sum(axis=1)).flatten()
    D = diags(degree, 0, shape=(lattice.num_states, lattice.num_states), format='csr')
    return scale * (D - adj)


def identify_2s_2p_eigenvalues(
    lattice: GeometricLattice,
    eigenvalues: np.ndarray,
    eigenvectors: np.ndarray,
) -> tuple:
    """
    Identify the eigenvalues corresponding to the 2s and 2p states
    by projecting eigenvectors onto the known basis states.

    The 2s state is (n=2, l=0, m=0).
    The 2p states are (n=2, l=1, m=-1), (n=2, l=1, m=0), (n=2, l=1, m=+1).

    Returns (E_2s, E_2p_mean) where E_2p_mean averages over the three 2p states.
    """
    # Get indices of the 2s and 2p basis states
    idx_2s = lattice._state_index.get((2, 0, 0))
    idx_2p = [
        lattice._state_index.get((2, 1, -1)),
        lattice._state_index.get((2, 1, 0)),
        lattice._state_index.get((2, 1, 1)),
    ]

    if idx_2s is None or any(i is None for i in idx_2p):
        raise ValueError("Lattice too small to contain n=2 states")

    # For each eigenstate, compute overlap with the 2s and 2p basis vectors.
    # The eigenstate with maximum |<2s|psi_k>|^2 is the 2s eigenvalue.
    n_eig = eigenvalues.shape[0]

    # 2s: find eigenstate with max weight on the (2,0,0) basis vector
    overlap_2s = np.abs(eigenvectors[idx_2s, :])**2
    best_2s = np.argmax(overlap_2s)
    E_2s = eigenvalues[best_2s]

    # 2p: find eigenstates with max total weight on the three 2p basis vectors
    overlap_2p = np.zeros(n_eig)
    for idx in idx_2p:
        overlap_2p += np.abs(eigenvectors[idx, :])**2

    # Exclude the 2s eigenstate from 2p candidates
    overlap_2p[best_2s] = 0.0
    best_2p = np.argmax(overlap_2p)
    E_2p = eigenvalues[best_2p]

    return E_2s, E_2p


def main() -> None:
    print("=" * 78)
    print("PAPER 1 CRITIQUE: 2s/2p Degeneracy Splitting as Discretization Artifact")
    print("=" * 78)
    print()
    print("Exact hydrogen: E(2s) = E(2p) = -1/8 Ha (degenerate)")
    print("Graph Laplacian lifts this degeneracy due to finite lattice truncation.")
    print("If ΔE -> 0 as V -> inf with a clean power law, it's an artifact.\n")

    n_max_values = list(range(5, 31))
    results = []

    # We need enough eigenvalues to capture the n=2 shell.
    # The n=1 shell has 1 state, n=2 has 4 states => first 5 eigenvalues suffice.
    # Request a few extra for safety.
    k_eig = 10

    print(f"{'n_max':>6s}  {'V':>6s}  {'E(2s) [Ha]':>14s}  {'E(2p) [Ha]':>14s}  "
          f"{'dE_rel [%]':>12s}")
    print("-" * 62)

    for n_max in n_max_values:
        lattice = GeometricLattice(max_n=n_max)
        V = lattice.num_states
        H = build_hamiltonian(lattice, UNIVERSAL_KINETIC_SCALE)

        # Compute lowest eigenvalues
        k = min(k_eig, V - 1)
        eigenvalues, eigenvectors = eigsh(H, k=k, which='SA')
        sort_idx = np.argsort(eigenvalues)
        eigenvalues = eigenvalues[sort_idx]
        eigenvectors = eigenvectors[:, sort_idx]

        try:
            E_2s, E_2p = identify_2s_2p_eigenvalues(lattice, eigenvalues, eigenvectors)
            dE_rel = abs(E_2p - E_2s) / abs(E_2s) * 100  # percentage
            results.append((n_max, V, E_2s, E_2p, dE_rel))
            print(f"{n_max:>6d}  {V:>6d}  {E_2s:>14.8f}  {E_2p:>14.8f}  {dE_rel:>12.6f}")
        except Exception as e:
            print(f"{n_max:>6d}  {V:>6d}  {'ERROR':>14s}  {str(e)}")

    if len(results) < 3:
        print("\nInsufficient data points for power-law fit. Exiting.")
        return

    # --- Power-law fit: ΔE_rel ~ V^(-alpha) ---
    V_arr = np.array([r[1] for r in results], dtype=float)
    dE_arr = np.array([r[4] for r in results], dtype=float)

    # Filter out zero or negative values for log-log fit
    mask = dE_arr > 0
    log_V = np.log10(V_arr[mask])
    log_dE = np.log10(dE_arr[mask])

    # Linear fit in log-log space: log(ΔE) = -alpha * log(V) + const
    coeffs = np.polyfit(log_V, log_dE, 1)
    alpha = -coeffs[0]
    intercept = coeffs[1]

    # R-squared
    log_dE_pred = np.polyval(coeffs, log_V)
    ss_res = np.sum((log_dE - log_dE_pred)**2)
    ss_tot = np.sum((log_dE - np.mean(log_dE))**2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    print()
    print("=" * 78)
    print("POWER-LAW FIT: dE_rel ~ V^(-alpha)")
    print(f"  alpha    = {alpha:.4f}")
    print(f"  R^2      = {r_squared:.6f}")
    print(f"  Intercept = 10^{intercept:.4f}")
    print()
    if alpha > 0 and r_squared > 0.8:
        print("  CONCLUSION: The 2s/2p splitting is a DISCRETIZATION ARTIFACT.")
        print(f"  It decays as V^(-{alpha:.2f}) and vanishes in the continuum limit.")
    else:
        print("  WARNING: Power-law fit inconclusive. Check data.")
    print("=" * 78)

    # --- Plot ---
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Left: Log-log plot of ΔE_rel vs V
        ax1 = axes[0]
        ax1.loglog(V_arr[mask], dE_arr[mask], 'ko', markersize=6, label='Computed')
        V_fit = np.logspace(np.log10(V_arr[mask].min()), np.log10(V_arr[mask].max()), 100)
        dE_fit = 10**intercept * V_fit**(-alpha)
        ax1.loglog(V_fit, dE_fit, 'r-', linewidth=2,
                   label=rf'Fit: $\Delta E \propto V^{{-{alpha:.2f}}}$ ($R^2={r_squared:.4f}$)')
        ax1.set_xlabel('Total vertex count V', fontsize=13)
        ax1.set_ylabel(r'$|\Delta E_{2p-2s}| / |E_{2s}|$ (%)', fontsize=13)
        ax1.set_title('2s/2p Splitting vs. Lattice Size', fontsize=14)
        ax1.legend(fontsize=11)
        ax1.grid(True, which='both', alpha=0.3)

        # Right: Eigenvalues vs n_max
        ax2 = axes[1]
        nmax_arr = np.array([r[0] for r in results])
        E2s_arr = np.array([r[2] for r in results])
        E2p_arr = np.array([r[3] for r in results])
        ax2.plot(nmax_arr, E2s_arr, 'bs-', markersize=5, label='E(2s)')
        ax2.plot(nmax_arr, E2p_arr, 'r^-', markersize=5, label='E(2p)')
        ax2.axhline(-1/8, color='green', linestyle='--', linewidth=1.5,
                     label=r'Exact $-1/8$ Ha')
        ax2.set_xlabel(r'$n_{\max}$', fontsize=13)
        ax2.set_ylabel('Energy (Ha)', fontsize=13)
        ax2.set_title('2s and 2p Eigenvalue Convergence', fontsize=14)
        ax2.legend(fontsize=11)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()

        plot_dir = os.path.join(os.path.dirname(__file__), 'plots')
        os.makedirs(plot_dir, exist_ok=True)
        plot_path = os.path.join(plot_dir, 'sp_splitting.png')
        plt.savefig(plot_path, dpi=150)
        print(f"\nPlot saved to: {plot_path}")
        plt.close()

    except ImportError:
        print("\nmatplotlib not available - skipping plot generation.")


if __name__ == '__main__':
    main()
