"""
Spectral angular solver for Level 4 molecule-frame hyperspherical method.

Replaces the FD-based solve_angular_multichannel() with a Jacobi polynomial
spectral basis, reducing matrix dimension from n_ch × n_alpha (~1000) to
n_ch × n_basis (~50).

Each (l1, l2) channel uses basis functions:
    φ_k^{(l1,l2)}(α) = N_k (cos α)^{l1+1} (sin α)^{l2+1} P_k^{(l2+½,l1+½)}(cos 2α)

where P_k^{(a,b)} is the Jacobi polynomial.  This reduces to the Level 3
Gegenbauer basis when l1 = l2 (since P_k^{(λ-½,λ-½)} ∝ C_k^λ).

The envelope (cos α)^{l1+1} (sin α)^{l2+1} absorbs the centrifugal
singularities l1(l1+1)/(2 cos²α) + l2(l2+1)/(2 sin²α), making all
integrands smooth.

Free eigenvalues are exact from SO(6) Casimir:
    μ_free(l1, l2, k) = ν(ν+4)/2  with  ν = l1 + l2 + 2k
    (= ½(l1+l2+2)² - 2 + 2k(k + l1 + l2 + 2))

The V_ee coupling is ρ-independent and precomputed once.  The nuclear
coupling depends on ρ = R/(2 R_e) and is recomputed per R_e point via
Gauss-Legendre quadrature on the spectral basis.

References:
  - Paper 15, Section V (multichannel expansion)
  - Paper 13, Section XII (algebraic angular structure, Level 3 precedent)
  - Track K plan (Level 4 angular speedup)
"""

import numpy as np
from scipy.linalg import eigh
from scipy.special import eval_jacobi
from typing import Tuple, List, Union, Optional
from math import sqrt

from geovac.hyperspherical_angular import gaunt_integral
from geovac.level4_multichannel import (
    _channel_list,
    _channel_list_extended,
    compute_nuclear_coupling,
    _ee_coupling,
    _ee_coupling_general,
    compute_core_screening_analytical,
    compute_pk_pseudopotential,
)


class Level4SpectralAngular:
    """Spectral Jacobi-polynomial solver for the Level 4 angular equation.

    Replaces the FD angular solver with a spectral basis, reducing matrix
    dimension from n_ch × n_alpha (~1000) to n_ch × n_basis (~50-100)
    while preserving eigenvalue accuracy to < 1e-4.

    The V_ee coupling matrix is precomputed once (ρ-independent).
    Nuclear coupling is recomputed per ρ-point via quadrature.

    Parameters
    ----------
    l_max : int
        Maximum angular momentum per electron.
    n_basis : int
        Number of Jacobi polynomial basis functions per channel.
    n_quad : int
        Number of Gauss-Legendre quadrature points per sub-interval.
        Total quadrature points = 2 × n_quad (split at π/4).
    Z_A, Z_B : float
        Nuclear charges.
    z0 : float
        Origin shift along internuclear axis.
    m_max : int
        Maximum |m| per electron.  0 = sigma-only.
    symmetry : str
        'singlet' restricts to exchange-symmetric k-indices for
        symmetric channels (l1=l2).  'all' includes all k.
    """

    def __init__(
        self,
        l_max: int = 2,
        n_basis: int = 10,
        n_quad: int = 100,
        Z_A: float = 1.0,
        Z_B: float = 1.0,
        z0: float = 0.0,
        m_max: int = 0,
        l_max_per_m: Union[dict, None] = None,
        symmetry: str = 'singlet',
    ) -> None:
        self.l_max = l_max
        self.n_basis = n_basis
        self.n_quad = n_quad
        self.Z_A = Z_A
        self.Z_B = Z_B
        self.z0 = z0
        self.m_max = m_max
        self.l_max_per_m = l_max_per_m
        self.symmetry = symmetry
        self.homonuclear = (Z_A == Z_B)

        # Build channel list
        if m_max == 0:
            channels_2 = _channel_list(l_max, homonuclear=self.homonuclear)
            self.channels_4 = [(l1, 0, l2, 0) for l1, l2 in channels_2]
            self.channels_2 = channels_2
        else:
            self.channels_4 = _channel_list_extended(
                l_max, m_max, l_max_per_m, homonuclear=self.homonuclear,
            )
            self.channels_2 = [(c[0], c[2]) for c in self.channels_4]
        self.n_ch = len(self.channels_4)

        self._setup_quadrature(n_quad)
        self._build_basis()
        self._precompute_vee()

    def _setup_quadrature(self, n_quad: int) -> None:
        """Setup Gauss-Legendre quadrature on [0,π/4] and [π/4,π/2].

        Split at π/4 handles the kink in max(cos,sin) and the
        split-region nuclear coupling boundary.
        """
        from numpy.polynomial.legendre import leggauss

        nodes, weights = leggauss(n_quad)

        a1, b1 = 0.0, np.pi / 4.0
        scale1 = (b1 - a1) / 2.0
        alpha1 = scale1 * nodes + (b1 + a1) / 2.0
        w1 = weights * scale1

        a2, b2 = np.pi / 4.0, np.pi / 2.0
        scale2 = (b2 - a2) / 2.0
        alpha2 = scale2 * nodes + (b2 + a2) / 2.0
        w2 = weights * scale2

        self._alpha = np.concatenate([alpha1, alpha2])
        self._weights = np.concatenate([w1, w2])
        self._sin_a = np.sin(self._alpha)
        self._cos_a = np.cos(self._alpha)
        self._cos_2a = np.cos(2.0 * self._alpha)

    def _build_basis(self) -> None:
        """Build orthonormalized Jacobi polynomial basis per channel.

        For channel (l1, l2):
            φ_k(α) = N_k (cos α)^{l1+1} (sin α)^{l2+1} P_k^{(a,b)}(cos 2α)
        with a = l2 + 0.5, b = l1 + 0.5.

        k-index selection:
          - Symmetric channels (l1=l2): even k for singlet, odd for triplet
          - Asymmetric channels (l1≠l2): all k values
        """
        self._channel_phi = []       # phi[ch][k, quad_point]
        self._channel_casimir = []   # free eigenvalues per channel
        self._channel_k_indices = [] # k values used per channel
        self._basis_offsets = []     # starting index in global basis

        offset = 0
        for ic, (l1, m1, l2, m2) in enumerate(self.channels_4):
            # Determine k-indices
            if l1 == l2 and self.symmetry == 'singlet':
                # Symmetric channel: even k for singlet
                k_indices = np.array([2 * j for j in range(self.n_basis)])
            elif l1 == l2 and self.symmetry == 'triplet':
                k_indices = np.array([2 * j + 1 for j in range(self.n_basis)])
            else:
                # Asymmetric channel or no symmetry filtering
                k_indices = np.arange(self.n_basis)

            # Jacobi parameters
            a_jac = l2 + 0.5  # associated with sin α
            b_jac = l1 + 0.5  # associated with cos α

            # Envelope: (cos α)^{l1+1} (sin α)^{l2+1}
            envelope = self._cos_a ** (l1 + 1) * self._sin_a ** (l2 + 1)

            # Free eigenvalues: μ = ½(l1+l2+2)² - 2 + 2k(k+l1+l2+2)
            L = l1 + l2 + 2
            casimir = 0.5 * L**2 - 2.0 + 2.0 * k_indices * (k_indices + L)

            # Build and normalize basis functions
            n_k = len(k_indices)
            phi = np.zeros((n_k, len(self._alpha)))
            for j, k in enumerate(k_indices):
                raw = envelope * eval_jacobi(k, a_jac, b_jac, self._cos_2a)
                norm_sq = np.dot(self._weights, raw * raw)
                if norm_sq > 1e-30:
                    phi[j] = raw / np.sqrt(norm_sq)
                else:
                    phi[j] = 0.0

            self._channel_phi.append(phi)
            self._channel_casimir.append(casimir)
            self._channel_k_indices.append(k_indices)
            self._basis_offsets.append(offset)
            offset += n_k

        self._total_dim = offset

    def _precompute_vee(self) -> None:
        """Precompute the V_ee coupling matrix (ρ-independent).

        Uses the existing _ee_coupling / _ee_coupling_general functions
        evaluated on the quadrature grid, then projected into the spectral
        basis.
        """
        dim = self._total_dim
        self._vee_matrix = np.zeros((dim, dim))

        for ic, (l1p, m1p, l2p, m2p) in enumerate(self.channels_4):
            phi_i = self._channel_phi[ic]
            n_i = phi_i.shape[0]
            i0 = self._basis_offsets[ic]

            for jc, (l1, m1, l2, m2) in enumerate(self.channels_4):
                if jc < ic:
                    continue

                phi_j = self._channel_phi[jc]
                n_j = phi_j.shape[0]
                j0 = self._basis_offsets[jc]

                # Evaluate V_ee on quadrature grid
                if self.m_max == 0:
                    W_ee = _ee_coupling(l1p, l2p, l1, l2, self._alpha, self.l_max)
                else:
                    W_ee = _ee_coupling_general(
                        l1p, m1p, l2p, m2p, l1, m1, l2, m2, self._alpha
                    )

                if np.max(np.abs(W_ee)) < 1e-15:
                    continue

                # Project: M[k',k] = Σ_q w_q φ_{k'}(α_q) W_ee(α_q) φ_k(α_q)
                weighted = phi_i * (self._weights * W_ee)[np.newaxis, :]
                block = weighted @ phi_j.T

                self._vee_matrix[i0:i0+n_i, j0:j0+n_j] = block
                if ic != jc:
                    self._vee_matrix[j0:j0+n_j, i0:i0+n_i] = block.T

    def _compute_nuclear_spectral(
        self,
        rho: float,
        R_e: float,
    ) -> np.ndarray:
        """Compute nuclear coupling matrix in spectral basis at given ρ.

        Evaluates the nuclear coupling V_nuc(α; ρ) on the quadrature grid
        using the existing algebraic multipole expansion, then projects
        into the spectral basis.

        Parameters
        ----------
        rho : float
            R / (2 R_e).
        R_e : float
            Electronic hyperradius.

        Returns
        -------
        V_nuc : ndarray of shape (total_dim, total_dim)
            Nuclear coupling in spectral basis.
        """
        dim = self._total_dim
        V_nuc = np.zeros((dim, dim))

        # Per-nucleus rho values for shifted origin
        rho_A = rho - self.z0 / R_e if self.z0 != 0.0 else rho
        rho_B = rho + self.z0 / R_e if self.z0 != 0.0 else rho

        for ic, (l1p, m1p, l2p, m2p) in enumerate(self.channels_4):
            phi_i = self._channel_phi[ic]
            n_i = phi_i.shape[0]
            i0 = self._basis_offsets[ic]

            for jc, (l1, m1, l2, m2) in enumerate(self.channels_4):
                if jc < ic:
                    continue

                phi_j = self._channel_phi[jc]
                n_j = phi_j.shape[0]
                j0 = self._basis_offsets[jc]

                # Selection rule: nuclear is diagonal in m
                if m1p != m1 or m2p != m2:
                    continue

                # Evaluate nuclear coupling on quadrature grid
                V_alpha = compute_nuclear_coupling(
                    l1p, l2p, l1, l2, m1, m2,
                    self._alpha, rho,
                    Z=self.Z_A,  # backward compat
                    Z_A=self.Z_A, Z_B=self.Z_B,
                    rho_A=rho_A, rho_B=rho_B,
                )

                if np.max(np.abs(V_alpha)) < 1e-15:
                    continue

                # Project into spectral basis
                weighted = phi_i * (self._weights * V_alpha)[np.newaxis, :]
                block = weighted @ phi_j.T

                V_nuc[i0:i0+n_i, j0:j0+n_j] = block
                if ic != jc:
                    V_nuc[j0:j0+n_j, i0:i0+n_i] = block.T

        return V_nuc

    def solve(
        self,
        rho: float,
        R_e: float,
        n_eig: int = 1,
        core_potentials: Union[List[dict], None] = None,
        pk_potentials: Union[List[dict], None] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Solve the angular eigenvalue problem at given (ρ, R_e).

        Assembles H = diag(μ_free) + R_e × (V_nuc(ρ) + V_ee)
        and diagonalizes.

        Parameters
        ----------
        rho : float
            R / (2 R_e).
        R_e : float
            Electronic hyperradius.
        n_eig : int
            Number of eigenvalues to return.
        core_potentials : list of dict or None
            Core screening specs.
        pk_potentials : list of dict or None
            Phillips-Kleinman pseudopotential specs.

        Returns
        -------
        mu : ndarray of shape (n_eig,)
            Lowest angular eigenvalues.
        channels : list
            Channel labels.
        """
        dim = self._total_dim

        # Free Hamiltonian: diagonal
        H = np.zeros((dim, dim))
        for ic in range(self.n_ch):
            n_k = len(self._channel_k_indices[ic])
            i0 = self._basis_offsets[ic]
            for j in range(n_k):
                H[i0 + j, i0 + j] = self._channel_casimir[ic][j]

        # Nuclear coupling (ρ-dependent)
        V_nuc = self._compute_nuclear_spectral(rho, R_e)
        H += R_e * V_nuc

        # V_ee coupling (precomputed)
        H += R_e * self._vee_matrix

        # Core screening (if applicable)
        if core_potentials:
            V_core = self._compute_core_spectral(R_e, rho, core_potentials)
            H += R_e * V_core

        # Phillips-Kleinman pseudopotential (if applicable)
        if pk_potentials:
            V_pk = self._compute_pk_spectral(R_e, rho, pk_potentials)
            H += R_e * V_pk

        # Eigensolve
        evals, evecs = eigh(H)

        if self.m_max == 0:
            channels = self.channels_2
        else:
            channels = self.channels_4

        return evals[:n_eig], channels

    def _compute_core_spectral(
        self,
        R_e: float,
        rho: float,
        core_potentials: List[dict],
    ) -> np.ndarray:
        """Core screening correction in spectral basis."""
        dim = self._total_dim
        V_core = np.zeros((dim, dim))

        rho_A = rho - self.z0 / R_e if self.z0 != 0.0 else rho
        rho_B = rho + self.z0 / R_e if self.z0 != 0.0 else rho

        V_pen_e1, V_pen_e2 = compute_core_screening_analytical(
            self._alpha, rho_A, R_e, core_potentials, rho_B=rho_B,
        )
        V_pen = V_pen_e1 + V_pen_e2

        # Core screening is diagonal in channels
        for ic in range(self.n_ch):
            phi = self._channel_phi[ic]
            n_k = phi.shape[0]
            i0 = self._basis_offsets[ic]

            weighted = phi * (self._weights * V_pen)[np.newaxis, :]
            block = weighted @ phi.T
            V_core[i0:i0+n_k, i0:i0+n_k] = block

        return V_core

    def _compute_pk_spectral(
        self,
        R_e: float,
        rho: float,
        pk_potentials: List[dict],
    ) -> np.ndarray:
        """Phillips-Kleinman pseudopotential in spectral basis."""
        dim = self._total_dim
        V_pk = np.zeros((dim, dim))

        rho_A = rho - self.z0 / R_e if self.z0 != 0.0 else rho
        rho_B = rho + self.z0 / R_e if self.z0 != 0.0 else rho

        V_pk_e1, V_pk_e2 = compute_pk_pseudopotential(
            self._alpha, rho_A, R_e, pk_potentials, rho_B=rho_B,
        )

        # Check for l-dependent PK mode
        pk_ch_mode = 'channel_blind'
        for p in pk_potentials:
            if p.get('channel_mode') == 'l_dependent':
                pk_ch_mode = 'l_dependent'
                break

        for ic, (l1, m1, l2, m2) in enumerate(self.channels_4):
            phi = self._channel_phi[ic]
            n_k = phi.shape[0]
            i0 = self._basis_offsets[ic]

            if pk_ch_mode == 'l_dependent':
                w1 = 1.0 if l1 == 0 else 0.0
                w2 = 1.0 if l2 == 0 else 0.0
            else:
                w1 = 1.0
                w2 = 1.0

            if w1 == 0.0 and w2 == 0.0:
                continue

            V_total = w1 * V_pk_e1 + w2 * V_pk_e2
            weighted = phi * (self._weights * V_total)[np.newaxis, :]
            block = weighted @ phi.T
            V_pk[i0:i0+n_k, i0:i0+n_k] = block

        return V_pk


def solve_angular_spectral(
    rho: float,
    R_e: float,
    l_max: int = 2,
    Z: float = 1.0,
    n_basis: int = 10,
    n_quad: int = 100,
    n_eig: int = 1,
    m_max: int = 0,
    l_max_per_m: Union[dict, None] = None,
    Z_A: float = None,
    Z_B: float = None,
    z0: float = 0.0,
    core_potentials: Union[List[dict], None] = None,
    pk_potentials: Union[List[dict], None] = None,
    symmetry: str = 'singlet',
    _solver: Optional['Level4SpectralAngular'] = None,
) -> Tuple[np.ndarray, list]:
    """
    Solve Level 4 angular equation using spectral Jacobi basis.

    Drop-in replacement for solve_angular_multichannel() but returns
    only eigenvalues and channel list (not eigenvectors or alpha grid).

    Parameters
    ----------
    _solver : Level4SpectralAngular or None
        Pre-built solver for reuse across R_e sweep.  If None, creates
        a new solver (slower, for one-off calls).

    Returns
    -------
    mu : ndarray of shape (n_eig,)
        Lowest angular eigenvalues.
    channels : list
        Channel labels.
    """
    if Z_A is None:
        Z_A = Z
    if Z_B is None:
        Z_B = Z

    if _solver is None:
        _solver = Level4SpectralAngular(
            l_max=l_max, n_basis=n_basis, n_quad=n_quad,
            Z_A=Z_A, Z_B=Z_B, z0=z0, m_max=m_max,
            l_max_per_m=l_max_per_m, symmetry=symmetry,
        )

    return _solver.solve(
        rho, R_e, n_eig=n_eig,
        core_potentials=core_potentials,
        pk_potentials=pk_potentials,
    )


def compute_adiabatic_curve_spectral(
    R: float,
    R_e_grid: np.ndarray,
    l_max: int = 2,
    Z: float = 1.0,
    n_basis: int = 10,
    n_quad: int = 100,
    m_max: int = 0,
    l_max_per_m: Union[dict, None] = None,
    Z_A: float = None,
    Z_B: float = None,
    z0: float = 0.0,
    core_potentials: Union[List[dict], None] = None,
    pk_potentials: Union[List[dict], None] = None,
    symmetry: str = 'singlet',
) -> np.ndarray:
    """
    Compute adiabatic potential U(R_e) using spectral angular solver.

    Drop-in replacement for compute_adiabatic_curve_mc() with spectral
    basis for the angular equation.

    Returns
    -------
    U : ndarray
        Effective potential U(R_e) in Ha.
    """
    if Z_A is None:
        Z_A = Z
    if Z_B is None:
        Z_B = Z

    # Build solver once, reuse for all R_e points
    solver = Level4SpectralAngular(
        l_max=l_max, n_basis=n_basis, n_quad=n_quad,
        Z_A=Z_A, Z_B=Z_B, z0=z0, m_max=m_max,
        l_max_per_m=l_max_per_m, symmetry=symmetry,
    )

    n_Re = len(R_e_grid)
    mu_vals = np.zeros(n_Re)

    for i, R_e in enumerate(R_e_grid):
        rho = R / (2.0 * R_e)
        evals, _ = solver.solve(
            rho, R_e, n_eig=1,
            core_potentials=core_potentials,
            pk_potentials=pk_potentials,
        )
        mu_vals[i] = evals[0]

    U = (mu_vals + 15.0 / 8.0) / R_e_grid**2
    return U
