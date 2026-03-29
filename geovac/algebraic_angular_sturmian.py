"""
Sturmian angular basis for Level 3 hyperspherical method.

Research prototype: replaces the free SO(6) eigenfunctions sin(2n*alpha) with
Sturmian functions that diagonalize H_0 + R_0 * V_nuclear_monopole.

The key insight: the Sturmian basis is constructed in a LARGER free basis
(n_construct >> n_basis), then the lowest n_basis Sturmian eigenvectors are
selected.  These n_basis vectors span a different (better-adapted) subspace
than the first n_basis free functions, because the Sturmian eigenvectors
are optimal linear combinations of many free functions that capture the
nuclear potential structure.

This is analogous to natural orbitals in quantum chemistry: the first n
natural orbitals (eigenvectors of the density matrix) span a better
n-dimensional subspace than the first n canonical orbitals.

Construction:
1. Build H_ref = diag(casimir) + R0 * V_nuclear in a large free basis (n_construct)
2. Diagonalize to get Sturmian eigenvectors U (columns of the rotation matrix)
3. Select the lowest n_basis Sturmian eigenvectors
4. Project the full Hamiltonian onto this truncated Sturmian subspace

This is NOT the same as the "Sturmian basis with shared p0" dead end in CLAUDE.md
Section 3.  That approach tried to use Sturmian functions for molecular encoding
(single S3 with shared momentum).  Here, the Sturmian basis is a rotation of the
angular eigenfunctions within the hyperspherical framework.

References:
  - Paper 13, Section XII (algebraic structure of angular equation)
  - Paper 12 (Neumann expansion precedent for algebraicization)
  - Abdouraman et al., J. Phys. B 49, 065004 (2016)
  - Aquilanti et al., Adv. Quantum Chem. 39, 103 (2001)
"""

import numpy as np
from scipy.linalg import eigh
from scipy.interpolate import CubicSpline
from typing import Tuple, Optional, Dict
from math import sqrt, pi

from geovac.algebraic_angular import AlgebraicAngularSolver


class SturmianAngularSolver:
    """Sturmian angular solver for Level 3 hyperspherical method.

    Constructs a Sturmian basis by diagonalizing H_free + R_0 * V_ref
    in a large free basis, then truncating to the lowest n_basis eigenvectors.

    Parameters
    ----------
    Z : float
        Nuclear charge.
    n_basis : int
        Number of Sturmian basis functions to keep per l-channel.
    l_max : int
        Maximum partial wave.
    R0 : float
        Reference hyperradius for Sturmian construction.
    n_construct : int
        Size of the free basis used to construct Sturmian functions.
        Must be >= n_basis.  Larger values give better Sturmian functions
        at the cost of a one-time diagonalization.
    ref_potential : str
        Which potential to include in the reference Hamiltonian:
        'nuclear': H_ref = H_free + R0 * V_nuclear (monopole Sturmian)
        'full': H_ref = H_free + R0 * V_coupling (full Sturmian, includes V_ee)
    symmetry : str
        'singlet' or 'triplet'.
    n_quad : int
        Quadrature points for underlying Gegenbauer solver.
    """

    def __init__(
        self,
        Z: float,
        n_basis: int = 10,
        l_max: int = 0,
        R0: float = 1.0,
        n_construct: int = 50,
        ref_potential: str = 'nuclear',
        symmetry: str = 'singlet',
        n_quad: int = 150,
    ) -> None:
        self.Z = Z
        self.n_basis = n_basis
        self.l_max = l_max
        self.R0 = R0
        self.n_construct = max(n_construct, n_basis)
        self.ref_potential = ref_potential
        self.symmetry = symmetry
        self.n_l = l_max + 1
        self._total_dim_kept = self.n_l * n_basis
        self._total_dim_full = self.n_l * self.n_construct

        # Build the large free-basis solver
        self._free_solver = AlgebraicAngularSolver(
            Z=Z, n_basis=self.n_construct, l_max=l_max,
            symmetry=symmetry, n_quad=n_quad,
        )

        # Extract matrices from the free solver (full size)
        self._casimir_free = np.concatenate(self._free_solver._channel_casimir)
        self._V_nuc_free = self._free_solver._nuclear_full.copy()
        self._V_ee_free = self._free_solver._vee_full.copy()
        self._V_coupling_free = self._free_solver._coupling_full.copy()

        # Decompose and build Sturmian basis
        self._decompose_coupling()
        self._build_sturmian_basis()
        self._project_to_sturmian()

    def _decompose_coupling(self) -> None:
        """Decompose V_coupling into nuclear monopole and remainder."""
        self._V_monopole = self._V_nuc_free.copy()
        self._V_remainder = self._V_ee_free.copy()

        nuc_norm = np.linalg.norm(self._V_monopole, 'fro')
        total_norm = np.linalg.norm(self._V_coupling_free, 'fro')
        vee_norm = np.linalg.norm(self._V_remainder, 'fro')

        self._monopole_fraction = nuc_norm / total_norm if total_norm > 0 else 0.0
        self._remainder_fraction = vee_norm / total_norm if total_norm > 0 else 0.0

    def _build_sturmian_basis(self) -> None:
        """Construct Sturmian basis in the large free basis, then truncate.

        Diagonalizes H_ref in the full n_construct-dimensional free basis.
        H_ref depends on ref_potential:
          'nuclear': H_ref = diag(casimir) + R0 * V_nuclear
          'full':    H_ref = diag(casimir) + R0 * V_coupling (nuclear + V_ee)

        For l_max > 0 with ref_potential='nuclear' (block-diagonal V_ref),
        uses per-channel selection: keeps the lowest n_basis eigenvectors
        within each l-channel independently, preventing l=0 bias.

        For ref_potential='full' (cross-channel coupling), uses global
        selection: keeps the lowest n_basis * n_l eigenvectors overall.
        """
        dim_full = self._total_dim_full
        dim_kept = self._total_dim_kept
        nc = self.n_construct
        nb = self.n_basis

        if self.ref_potential == 'nuclear':
            V_ref = self._V_monopole
        elif self.ref_potential == 'full':
            V_ref = self._V_coupling_free
        else:
            raise ValueError(f"Unknown ref_potential: {self.ref_potential}")

        self._V_ref = V_ref
        H_ref = np.diag(self._casimir_free) + self.R0 * V_ref

        if self.n_l == 1 or self.ref_potential == 'full':
            # Single channel or cross-channel ref: global selection
            mu_all, U_all = eigh(H_ref)
            self._mu_sturmian = mu_all[:dim_kept]
            self._U_truncated = U_all[:, :dim_kept]
        else:
            # Per-channel selection for block-diagonal V_ref
            # This ensures each l-channel gets exactly n_basis Sturmian functions
            mu_list = []
            U_blocks = []
            for l in range(self.n_l):
                i0, i1 = l * nc, (l + 1) * nc
                H_block = H_ref[i0:i1, i0:i1]
                mu_block, U_block = eigh(H_block)
                mu_list.append(mu_block[:nb])
                U_blocks.append(U_block[:, :nb])

            # Assemble the truncated rotation matrix (block-diagonal)
            self._mu_sturmian = np.concatenate(mu_list)
            self._U_truncated = np.zeros((dim_full, dim_kept))
            for l in range(self.n_l):
                i0_full, i1_full = l * nc, (l + 1) * nc
                j0_kept, j1_kept = l * nb, (l + 1) * nb
                self._U_truncated[i0_full:i1_full, j0_kept:j1_kept] = U_blocks[l]

    def _project_to_sturmian(self) -> None:
        """Project all Hamiltonian components onto the truncated Sturmian subspace.

        The projection matrix P = U_trunc (dim_full x dim_kept).
        Any operator A in the full space becomes A_proj = P^T @ A @ P
        in the truncated Sturmian subspace.
        """
        P = self._U_truncated

        # Project all operators
        self._V_monopole_proj = P.T @ self._V_monopole @ P
        self._V_remainder_proj = P.T @ self._V_remainder @ P
        self._V_coupling_proj = P.T @ self._V_coupling_free @ P
        self._V_ref_proj = P.T @ self._V_ref @ P
        self._casimir_proj = P.T @ np.diag(self._casimir_free) @ P

        # Residual coupling: part of V_coupling NOT in V_ref
        self._V_residual_proj = self._V_coupling_proj - self._V_ref_proj

        # Verify H_ref is diagonal in Sturmian basis
        H_ref_proj = self._casimir_proj + self.R0 * self._V_ref_proj
        off_diag = H_ref_proj - np.diag(np.diag(H_ref_proj))
        self._projection_error = np.max(np.abs(off_diag))

    def _build_H(self, R: float) -> np.ndarray:
        """Build projected Hamiltonian at hyperradius R.

        H(R) = diag(mu_sturm) + (R - R0) * V_ref_proj + R * V_residual_proj

        This is equivalent to casimir_proj + R * V_coupling_proj but
        exploits the fact that diag(mu_sturm) is diagonal.
        """
        return (np.diag(self._mu_sturmian)
                + (R - self.R0) * self._V_ref_proj
                + R * self._V_residual_proj)

    def solve(
        self, R: float, n_channels: int = 1
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Solve the angular eigenvalue problem at hyperradius R.

        Parameters
        ----------
        R : float
            Hyperradius (bohr).
        n_channels : int
            Number of eigenvalues/eigenvectors to return.

        Returns
        -------
        eigenvalues : ndarray of shape (n_channels,)
        eigenvectors : ndarray of shape (n_channels, dim_kept)
        """
        evals, evecs = eigh(self._build_H(R))
        return evals[:n_channels], evecs[:, :n_channels].T

    def solve_with_dboc(
        self, R: float, n_channels: int = 1
    ) -> Tuple[np.ndarray, np.ndarray, float]:
        """Solve angular problem and compute DBOC for the ground channel."""
        evals, evecs = eigh(self._build_H(R))

        phi_0 = evecs[:, 0]
        V_phi_0 = self._V_coupling_proj @ phi_0

        dboc = 0.0
        for mu in range(1, len(evals)):
            numerator = evecs[:, mu] @ V_phi_0
            denominator = evals[mu] - evals[0]
            if abs(denominator) > 1e-12:
                P_mu0 = numerator / denominator
                dboc += P_mu0 ** 2

        dboc *= 0.5
        return evals[:n_channels], evecs[:, :n_channels].T, dboc

    def monopole_decomposition(self) -> Dict[str, float]:
        """Return the monopole decomposition fractions."""
        return {
            'monopole_fraction': self._monopole_fraction,
            'remainder_fraction': self._remainder_fraction,
            'V_nuc_frobenius': np.linalg.norm(self._V_monopole, 'fro'),
            'V_ee_frobenius': np.linalg.norm(self._V_remainder, 'fro'),
            'V_total_frobenius': np.linalg.norm(self._V_coupling_free, 'fro'),
            'projection_error': self._projection_error,
            'n_construct': self.n_construct,
            'n_basis': self.n_basis,
        }

    def off_diagonal_analysis(self) -> Dict[str, float]:
        """Analyze off-diagonal coupling in the projected Sturmian basis."""
        V_mono = self._V_monopole_proj
        V_rem = self._V_remainder_proj

        mono_offdiag = V_mono - np.diag(np.diag(V_mono))
        rem_offdiag = V_rem - np.diag(np.diag(V_rem))

        return {
            'monopole_offdiag_norm': np.linalg.norm(mono_offdiag, 'fro'),
            'monopole_diag_norm': np.linalg.norm(np.diag(V_mono)),
            'remainder_offdiag_norm': np.linalg.norm(rem_offdiag, 'fro'),
            'remainder_diag_norm': np.linalg.norm(np.diag(V_rem)),
            'total_offdiag_norm': np.linalg.norm(
                mono_offdiag + rem_offdiag, 'fro'),
        }


def solve_hyperspherical_sturmian(
    Z: float = 2.0,
    n_electrons: int = 2,
    n_basis: int = 10,
    l_max: int = 0,
    R0: float = 1.0,
    n_construct: int = 50,
    ref_potential: str = 'nuclear',
    n_R: int = 200,
    R_min: float = 0.1,
    R_max: float = 30.0,
    N_R_radial: int = 2000,
    include_dboc: bool = False,
    verbose: bool = True,
) -> dict:
    """Full Level 3 solver using Sturmian angular basis.

    Drop-in replacement for solve_hyperspherical_algebraic, but using
    the SturmianAngularSolver with truncated Sturmian subspace.

    Parameters
    ----------
    Z : float
        Nuclear charge.
    n_electrons : int
        Number of electrons (must be 2).
    n_basis : int
        Number of Sturmian basis functions to keep per l-channel.
    l_max : int
        Maximum partial wave.
    R0 : float
        Reference hyperradius for Sturmian construction.
    n_construct : int
        Size of the free basis for constructing Sturmian functions.
    n_R : int
        Number of R points for computing mu(R).
    R_min, R_max : float
        Hyperradial grid boundaries.
    N_R_radial : int
        Number of grid points for the radial solve.
    include_dboc : bool
        If True, include DBOC in the effective potential.
    verbose : bool
        Print progress information.

    Returns
    -------
    result : dict
    """
    import time
    from geovac.hyperspherical_adiabatic import effective_potential
    from geovac.hyperspherical_radial import solve_radial

    if n_electrons != 2:
        raise ValueError("Only two-electron systems are supported.")

    t0 = time.time()

    solver = SturmianAngularSolver(
        Z, n_basis, l_max, R0=R0, n_construct=n_construct,
        ref_potential=ref_potential,
    )

    R_grid = np.concatenate([
        np.linspace(R_min, 1.0, n_R // 3),
        np.linspace(1.0, 5.0, n_R // 3),
        np.linspace(5.0, R_max, n_R // 3 + 1),
    ])
    R_grid = np.unique(R_grid)

    dboc_mode = "with DBOC" if include_dboc else "adiabatic"
    if verbose:
        decomp = solver.monopole_decomposition()
        print(f"Sturmian angular solver: Z={Z}, n_basis={n_basis}, "
              f"n_construct={n_construct}, l_max={l_max}, R0={R0} ({dboc_mode})")
        print(f"  Monopole fraction: {decomp['monopole_fraction']:.1%}")
        print(f"  Projection error: {decomp['projection_error']:.2e}")
        print(f"Computing mu(R) on {len(R_grid)} R points...")

    mu = np.zeros(len(R_grid))
    dboc_values = np.zeros(len(R_grid)) if include_dboc else None

    for i, R in enumerate(R_grid):
        if include_dboc:
            evals, _, dboc = solver.solve_with_dboc(R, n_channels=1)
            mu[i] = evals[0]
            dboc_values[i] = dboc
        else:
            evals, _ = solver.solve(R, n_channels=1)
            mu[i] = evals[0]

    t1 = time.time()
    if verbose:
        print(f"  Angular solve: {t1 - t0:.2f}s")
        V_eff_last = effective_potential(R_grid[-1:], mu[-1:])[0]
        print(
            f"  V_eff(R={R_grid[-1]:.1f}) = {V_eff_last:.4f} Ha "
            f"(should approach {-Z**2 / 2:.1f})"
        )

    V_eff = effective_potential(R_grid, mu)
    if include_dboc:
        V_eff_uncorrected = V_eff.copy()
        V_eff = V_eff + dboc_values
    V_eff_spline = CubicSpline(R_grid, V_eff, extrapolate=True)

    if verbose:
        print(f"Solving hyperradial equation (N_R={N_R_radial})...")

    E, F, R_grid_rad = solve_radial(
        V_eff_spline, R_min, R_max, N_R_radial, n_states=1
    )

    t2 = time.time()
    E_exact = -2.903724
    if verbose:
        print(f"  Radial solve: {t2 - t1:.2f}s")
        print(f"\n  Ground state energy: {E[0]:.6f} Ha")
        print(f"  Exact He:            {E_exact:.6f} Ha")
        err = abs(E[0] - E_exact) / abs(E_exact) * 100
        print(f"  Error:               {err:.4f}%")
        print(f"  Total time:          {t2 - t0:.2f}s")

    result = {
        'energy': E[0],
        'R_grid_angular': R_grid,
        'mu_curve': mu,
        'V_eff': V_eff,
        'R_grid_radial': R_grid_rad,
        'wavefunction': F[0],
        'solver': solver,
        'Z': Z,
        'n_basis': n_basis,
        'n_construct': n_construct,
        'l_max': l_max,
        'R0': R0,
        'include_dboc': include_dboc,
    }

    if include_dboc:
        V_eff_spline_nocorr = CubicSpline(
            R_grid, V_eff_uncorrected, extrapolate=True
        )
        E_no, _, _ = solve_radial(
            V_eff_spline_nocorr, R_min, R_max, N_R_radial, n_states=1
        )
        result['dboc_curve'] = dboc_values
        result['energy_no_dboc'] = E_no[0]

    return result


def optimize_R0(
    Z: float = 2.0,
    n_basis: int = 10,
    l_max: int = 0,
    n_construct: int = 50,
    R0_values: Optional[list] = None,
    n_R: int = 150,
    N_R_radial: int = 1500,
    verbose: bool = True,
) -> Dict[str, object]:
    """Find the optimal reference hyperradius R0 for the Sturmian basis.

    Scans R0 values and picks the one giving the lowest He energy.
    """
    if R0_values is None:
        R0_values = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]

    energies = []
    for R0 in R0_values:
        result = solve_hyperspherical_sturmian(
            Z=Z, n_basis=n_basis, l_max=l_max, R0=R0,
            n_construct=n_construct,
            n_R=n_R, N_R_radial=N_R_radial, verbose=False,
        )
        energies.append(result['energy'])
        if verbose:
            E_exact = -2.903724
            err = abs(result['energy'] - E_exact) / abs(E_exact) * 100
            print(f"  R0={R0:.1f}: E={result['energy']:.6f} Ha, err={err:.4f}%")

    best_idx = np.argmin(energies)
    return {
        'best_R0': R0_values[best_idx],
        'best_energy': energies[best_idx],
        'R0_values': R0_values,
        'energies': energies,
    }
