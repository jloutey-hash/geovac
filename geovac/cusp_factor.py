"""
Kato cusp factor for Level 4 molecule-frame hyperspherical angular equation.

In Level 4 coordinates, r_12 = R_e * sin(2*alpha) * g(theta_12), so at the
coalescence point theta_12 = 0, the electron-electron cusp condition becomes
a function of the single coordinate alpha.  The multiplicative cusp factor

    f(alpha) = 1 + (R_e / 2) * sin(2*alpha)

satisfies the Kato cusp condition (1/psi)(dpsi/dr_12)|_{r_12=0} = 1/2 in
the alpha projection.  The similarity-transformed Hamiltonian

    H_tilde = f^{-1} H f

has no 1/r_12 divergence at alpha = pi/4, so its partial-wave expansion
should converge faster in l_max.

The transform preserves gerade/ungerade symmetry because
f(alpha) = f(pi/2 - alpha) (sin(2*alpha) is symmetric about pi/4).

The transform introduces two additional terms from the kinetic energy:
    -f'/f * d/d(alpha)        (first-derivative/advection term)
    -(1/2) * f''/f            (potential correction)

where f' = R_e * cos(2*alpha) and f'' = -2 * R_e * sin(2*alpha).

The first-derivative term breaks the Hermiticity of the kinetic operator
in the original (unweighted) inner product.  To restore Hermiticity, we
work in the f^2-weighted inner product: <u|v>_f = integral u(a) v(a) f(a)^2 da.
In this inner product, f^{-1} H f is self-adjoint.

Implementation: we modify the spectral basis functions by multiplying them
by f(alpha), forming phi_tilde_k = f * phi_k, and then project the ORIGINAL
Hamiltonian H into this modified basis.  The generalized eigenvalue problem
    H_tilde c = mu S_tilde c
where H_tilde_{ij} = <f*phi_i | H | f*phi_j> and S_tilde_{ij} = <f*phi_i | f*phi_j>
gives the same eigenvalues as the similarity-transformed problem but
maintains manifest Hermiticity.

References:
  - Paper 15, Section III (cusp structure in Level 4 coordinates)
  - Paper 12, Section VII (cusp diagnosis: 7.6% D_e gap)
  - Kato, Commun. Pure Appl. Math. 10, 151 (1957)
"""

import numpy as np
from scipy.linalg import eigh
from scipy.special import eval_jacobi
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh_tridiagonal
from typing import Tuple, List, Union, Optional

from geovac.level4_multichannel import (
    _channel_list,
    _channel_list_extended,
    compute_nuclear_coupling,
    _ee_coupling,
    _ee_coupling_general,
    compute_core_screening_analytical,
    compute_pk_pseudopotential,
)
from geovac.hyperspherical_angular import gaunt_integral


class CuspFactorSolver:
    """Spectral angular solver with Kato cusp factor for Level 4.

    Wraps the spectral Jacobi basis approach from Level4SpectralAngular,
    multiplying each basis function by the cusp factor f(alpha) to
    absorb the electron-electron coalescence singularity.

    The cusp factor is:
        f(alpha; R_e) = 1 + gamma * sin(2*alpha)

    where gamma = R_e / 2 (Kato cusp condition).  When gamma=0, this
    reduces exactly to the standard spectral solver.

    The solver uses a generalized eigenvalue problem:
        H_mod c = mu S_mod c
    where H_mod and S_mod are the Hamiltonian and overlap in the
    f-modified basis.

    Parameters
    ----------
    l_max : int
        Maximum angular momentum per electron.
    n_basis : int
        Number of Jacobi polynomial basis functions per channel.
    n_quad : int
        Gauss-Legendre quadrature points per sub-interval.
    Z_A, Z_B : float
        Nuclear charges.
    z0 : float
        Origin shift along internuclear axis.
    m_max : int
        Maximum |m| per electron.
    symmetry : str
        'singlet' or 'triplet' for k-index selection in symmetric channels.
    cusp_gamma : float or None
        Cusp factor strength.  If None, uses R_e/2 (Kato condition).
        If 0.0, disables cusp factor (reproduces baseline).
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
        cusp_gamma: Optional[float] = None,
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
        self._cusp_gamma = cusp_gamma  # None = use R_e/2; 0 = disabled

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
        self._build_bare_basis()

    def _setup_quadrature(self, n_quad: int) -> None:
        """Setup Gauss-Legendre quadrature on [0,pi/4] and [pi/4,pi/2]."""
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
        self._sin_2a = np.sin(2.0 * self._alpha)

    def _build_bare_basis(self) -> None:
        """Build un-normalized Jacobi polynomial basis per channel (without cusp factor).

        These are stored un-normalized; normalization and cusp factor are
        applied in _build_modified_basis() which is called per R_e point.
        """
        self._channel_phi_bare = []
        self._channel_casimir = []
        self._channel_k_indices = []
        self._basis_offsets = []

        offset = 0
        for ic, (l1, m1, l2, m2) in enumerate(self.channels_4):
            # k-index selection (same as Level4SpectralAngular)
            if l1 == l2 and self.symmetry == 'singlet':
                k_indices = np.array([2 * j for j in range(self.n_basis)])
            elif l1 == l2 and self.symmetry == 'triplet':
                k_indices = np.array([2 * j + 1 for j in range(self.n_basis)])
            else:
                k_indices = np.arange(self.n_basis)

            a_jac = l2 + 0.5
            b_jac = l1 + 0.5

            envelope = self._cos_a ** (l1 + 1) * self._sin_a ** (l2 + 1)

            L = l1 + l2 + 2
            casimir = 0.5 * L**2 - 2.0 + 2.0 * k_indices * (k_indices + L)

            n_k = len(k_indices)
            phi_bare = np.zeros((n_k, len(self._alpha)))
            for j, k in enumerate(k_indices):
                phi_bare[j] = envelope * eval_jacobi(k, a_jac, b_jac, self._cos_2a)

            self._channel_phi_bare.append(phi_bare)
            self._channel_casimir.append(casimir)
            self._channel_k_indices.append(k_indices)
            self._basis_offsets.append(offset)
            offset += n_k

        self._total_dim = offset

    def _cusp_factor(self, R_e: float) -> np.ndarray:
        """Compute the cusp factor f(alpha; R_e) on the quadrature grid.

        f(alpha) = 1 + gamma * sin(2*alpha)

        where gamma = R_e/2 (Kato) or a user-specified value.

        Returns
        -------
        f : ndarray of shape (n_quad_total,)
        """
        if self._cusp_gamma is not None:
            gamma = self._cusp_gamma
        else:
            gamma = R_e / 2.0

        return 1.0 + gamma * self._sin_2a

    def _cusp_factor_deriv(self, R_e: float) -> Tuple[np.ndarray, np.ndarray]:
        """Compute f'(alpha) and f''(alpha) on the quadrature grid.

        f'  = 2 * gamma * cos(2*alpha)
        f'' = -4 * gamma * sin(2*alpha)

        Returns
        -------
        fp : ndarray  -- f'(alpha)
        fpp : ndarray -- f''(alpha)
        """
        if self._cusp_gamma is not None:
            gamma = self._cusp_gamma
        else:
            gamma = R_e / 2.0

        fp = 2.0 * gamma * self._cos_2a
        fpp = -4.0 * gamma * self._sin_2a
        return fp, fpp

    def _build_modified_basis(self, R_e: float) -> List[np.ndarray]:
        """Build cusp-modified basis: phi_tilde_k = f(alpha; R_e) * phi_bare_k.

        These are NOT pre-normalized because the overlap matrix S_tilde
        will be used in the generalized eigenvalue problem.

        Returns
        -------
        channel_phi_mod : list of ndarray, each (n_k, n_quad)
        """
        f = self._cusp_factor(R_e)
        channel_phi_mod = []
        for ic in range(self.n_ch):
            phi_bare = self._channel_phi_bare[ic]
            phi_mod = phi_bare * f[np.newaxis, :]
            channel_phi_mod.append(phi_mod)
        return channel_phi_mod

    def _build_overlap_matrix(
        self,
        channel_phi: List[np.ndarray],
    ) -> np.ndarray:
        """Build overlap matrix S_{ij} = integral phi_i(a) phi_j(a) da.

        Parameters
        ----------
        channel_phi : list of ndarray
            Modified basis functions per channel.

        Returns
        -------
        S : ndarray of shape (total_dim, total_dim)
        """
        dim = self._total_dim
        S = np.zeros((dim, dim))

        for ic in range(self.n_ch):
            phi = channel_phi[ic]
            n_k = phi.shape[0]
            i0 = self._basis_offsets[ic]

            # Intra-channel overlap
            weighted = phi * self._weights[np.newaxis, :]
            block = weighted @ phi.T
            S[i0:i0+n_k, i0:i0+n_k] = block

        # Cross-channel overlap is zero because the Jacobi envelope
        # functions for different (l1,l2) channels are in different
        # function spaces (different powers of cos/sin). The cusp
        # factor f(alpha) preserves this because it's channel-independent.
        # However, the modified basis IS NOT orthogonal within a channel
        # when gamma != 0.

        return S

    def _build_kinetic_matrix(
        self,
        channel_phi: List[np.ndarray],
        R_e: float,
    ) -> np.ndarray:
        """Build kinetic + centrifugal + Liouville matrix in modified basis.

        Uses a manifestly symmetric formulation:

        1. Kinetic energy via integration by parts (symmetric by construction):
            K_{ij} = (1/2) integral phi_tilde_i'(a) phi_tilde_j'(a) da

           where phi_tilde = f * phi_bare and derivatives are analytic.

        2. Centrifugal + Liouville as multiplicative potential (symmetric):
            V_{ij} = integral phi_tilde_i(a) V_cent(a) phi_tilde_j(a) da

           with V_cent = l1(l1+1)/(2cos^2(a)) + l2(l2+1)/(2sin^2(a)) - 2.

        At gamma=0, phi_tilde = phi_bare and this reproduces the Casimir
        eigenvalues to quadrature precision.
        """
        dim = self._total_dim
        H_kin = np.zeros((dim, dim))
        f = self._cusp_factor(R_e)
        fp, fpp = self._cusp_factor_deriv(R_e)

        for ic, (l1, m1, l2, m2) in enumerate(self.channels_4):
            phi_mod = channel_phi[ic]  # f * phi_bare
            phi_bare = self._channel_phi_bare[ic]
            n_k = phi_mod.shape[0]
            i0 = self._basis_offsets[ic]

            # Centrifugal + Liouville (multiplicative, symmetric)
            V_cent = (0.5 * l1 * (l1 + 1) / self._cos_a**2
                      + 0.5 * l2 * (l2 + 1) / self._sin_a**2
                      - 2.0)

            weighted_V = phi_mod * (self._weights * V_cent)[np.newaxis, :]
            M_V = weighted_V @ phi_mod.T
            H_kin[i0:i0+n_k, i0:i0+n_k] += M_V

            # Kinetic energy via integration by parts (symmetric):
            # K_{ij} = (1/2) int phi_tilde_i' phi_tilde_j' da
            # phi_tilde' = f' * phi_bare + f * phi_bare'
            dphi_bare = self._analytic_derivative(ic, phi_bare)
            dphi_mod = fp[np.newaxis, :] * phi_bare + f[np.newaxis, :] * dphi_bare

            weighted_d = dphi_mod * self._weights[np.newaxis, :]
            M_K = 0.5 * (weighted_d @ dphi_mod.T)
            H_kin[i0:i0+n_k, i0:i0+n_k] += M_K

        return H_kin

    def _analytic_derivative(
        self, ic: int, phi_bare: np.ndarray,
    ) -> np.ndarray:
        """Compute d(phi_bare)/d(alpha) analytically at quadrature points.

        For channel (l1, l2) with Jacobi parameters a=l2+0.5, b=l1+0.5:
            phi_k(a) = env(a) * P_k^{(a_j,b_j)}(cos 2a)

        where env(a) = cos^{l1+1}(a) * sin^{l2+1}(a).

        Derivative:
            phi_k'(a) = env'(a) * P_k(x) + env(a) * dP_k/dx * (-2 sin 2a)

        env'(a) = cos^l1 * sin^l2 * [-(l1+1)*sin^2(a) + (l2+1)*cos^2(a)]

        dP_k^{(a,b)}/dx = (k+a+b+1)/2 * P_{k-1}^{(a+1,b+1)}(x)
        (DLMF 18.9.15, with dP_0/dx = 0)
        """
        l1, m1, l2, m2 = self.channels_4[ic]
        k_indices = self._channel_k_indices[ic]
        a_jac = l2 + 0.5
        b_jac = l1 + 0.5

        cos_a = self._cos_a
        sin_a = self._sin_a
        cos_2a = self._cos_2a
        sin_2a = self._sin_2a

        env = cos_a ** (l1 + 1) * sin_a ** (l2 + 1)
        env_deriv = (cos_a ** l1 * sin_a ** l2
                     * (-(l1 + 1) * sin_a**2 + (l2 + 1) * cos_a**2))

        n_k = len(k_indices)
        dphi = np.zeros_like(phi_bare)

        for j, k in enumerate(k_indices):
            P_k = eval_jacobi(k, a_jac, b_jac, cos_2a)

            if k == 0:
                dP_k_da = np.zeros_like(cos_2a)
            else:
                dP_k_dx = ((k + a_jac + b_jac + 1) / 2.0
                           * eval_jacobi(k - 1, a_jac + 1, b_jac + 1, cos_2a))
                dP_k_da = dP_k_dx * (-2.0 * sin_2a)

            dphi[j] = env_deriv * P_k + env * dP_k_da

        return dphi

    def _compute_vee_modified(
        self,
        channel_phi: List[np.ndarray],
    ) -> np.ndarray:
        """Compute V_ee in the cusp-modified basis.

        V_ee is rho-independent: we precompute it each time R_e changes
        (because the basis changes with R_e through the cusp factor).
        """
        dim = self._total_dim
        V_ee = np.zeros((dim, dim))

        for ic, (l1p, m1p, l2p, m2p) in enumerate(self.channels_4):
            phi_i = channel_phi[ic]
            n_i = phi_i.shape[0]
            i0 = self._basis_offsets[ic]

            for jc, (l1, m1, l2, m2) in enumerate(self.channels_4):
                if jc < ic:
                    continue

                phi_j = channel_phi[jc]
                n_j = phi_j.shape[0]
                j0 = self._basis_offsets[jc]

                if self.m_max == 0:
                    W_ee = _ee_coupling(l1p, l2p, l1, l2, self._alpha, self.l_max)
                else:
                    W_ee = _ee_coupling_general(
                        l1p, m1p, l2p, m2p, l1, m1, l2, m2, self._alpha,
                    )

                if np.max(np.abs(W_ee)) < 1e-15:
                    continue

                weighted = phi_i * (self._weights * W_ee)[np.newaxis, :]
                block = weighted @ phi_j.T

                V_ee[i0:i0+n_i, j0:j0+n_j] = block
                if ic != jc:
                    V_ee[j0:j0+n_j, i0:i0+n_i] = block.T

        return V_ee

    def _compute_nuclear_modified(
        self,
        channel_phi: List[np.ndarray],
        rho: float,
        R_e: float,
    ) -> np.ndarray:
        """Compute nuclear coupling in cusp-modified basis."""
        dim = self._total_dim
        V_nuc = np.zeros((dim, dim))

        rho_A = rho - self.z0 / R_e if self.z0 != 0.0 else rho
        rho_B = rho + self.z0 / R_e if self.z0 != 0.0 else rho

        for ic, (l1p, m1p, l2p, m2p) in enumerate(self.channels_4):
            phi_i = channel_phi[ic]
            n_i = phi_i.shape[0]
            i0 = self._basis_offsets[ic]

            for jc, (l1, m1, l2, m2) in enumerate(self.channels_4):
                if jc < ic:
                    continue

                phi_j = channel_phi[jc]
                n_j = phi_j.shape[0]
                j0 = self._basis_offsets[jc]

                if m1p != m1 or m2p != m2:
                    continue

                V_alpha = compute_nuclear_coupling(
                    l1p, l2p, l1, l2, m1, m2,
                    self._alpha, rho,
                    Z=self.Z_A,
                    Z_A=self.Z_A, Z_B=self.Z_B,
                    rho_A=rho_A, rho_B=rho_B,
                )

                if np.max(np.abs(V_alpha)) < 1e-15:
                    continue

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
    ) -> Tuple[np.ndarray, list]:
        """Solve the angular eigenvalue problem with cusp factor.

        Assembles H_mod and S_mod in the cusp-modified basis, then
        solves the generalized eigenvalue problem H_mod c = mu S_mod c.

        Parameters
        ----------
        rho : float
            R / (2 R_e).
        R_e : float
            Electronic hyperradius.
        n_eig : int
            Number of eigenvalues to return.

        Returns
        -------
        mu : ndarray of shape (n_eig,)
            Lowest angular eigenvalues.
        channels : list
            Channel labels.
        """
        # Build cusp-modified basis for this R_e
        channel_phi = self._build_modified_basis(R_e)

        # Overlap matrix
        S = self._build_overlap_matrix(channel_phi)

        # Kinetic + centrifugal + Liouville (includes cusp correction terms)
        H = self._build_kinetic_matrix(channel_phi, R_e)

        # Nuclear coupling
        V_nuc = self._compute_nuclear_modified(channel_phi, rho, R_e)
        H += R_e * V_nuc

        # V_ee coupling
        V_ee = self._compute_vee_modified(channel_phi)
        H += R_e * V_ee

        # Generalized eigenvalue problem
        evals, evecs = eigh(H, S)

        if self.m_max == 0:
            channels = self.channels_2
        else:
            channels = self.channels_4

        return evals[:n_eig], channels

    def solve_baseline(
        self,
        rho: float,
        R_e: float,
        n_eig: int = 1,
    ) -> Tuple[np.ndarray, list]:
        """Solve without cusp factor (gamma=0) for comparison.

        Temporarily disables cusp factor and solves, providing the
        baseline for A/B comparison.
        """
        saved_gamma = self._cusp_gamma
        self._cusp_gamma = 0.0
        try:
            result = self.solve(rho, R_e, n_eig=n_eig)
        finally:
            self._cusp_gamma = saved_gamma
        return result


def compute_adiabatic_curve_cusp(
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
    cusp_gamma: Optional[float] = None,
    symmetry: str = 'singlet',
) -> np.ndarray:
    """Compute adiabatic potential U(R_e) with cusp factor.

    U(R_e; R) = [mu(R_e; R) + 15/8] / R_e^2

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    R_e_grid : ndarray
        Electronic hyperradius grid points.
    cusp_gamma : float or None
        Cusp factor strength. None = R_e/2, 0 = disabled.

    Returns
    -------
    U : ndarray
        Effective potential U(R_e) in Ha.
    """
    if Z_A is None:
        Z_A = Z
    if Z_B is None:
        Z_B = Z

    solver = CuspFactorSolver(
        l_max=l_max, n_basis=n_basis, n_quad=n_quad,
        Z_A=Z_A, Z_B=Z_B, z0=z0, m_max=m_max,
        l_max_per_m=l_max_per_m, symmetry=symmetry,
        cusp_gamma=cusp_gamma,
    )

    n_Re = len(R_e_grid)
    mu_vals = np.zeros(n_Re)

    for i, R_e in enumerate(R_e_grid):
        rho = R / (2.0 * R_e)
        evals, _ = solver.solve(rho, R_e, n_eig=1)
        mu_vals[i] = evals[0]

    U = (mu_vals + 15.0 / 8.0) / R_e_grid**2
    return U


def solve_level4_h2_cusp(
    R: float,
    l_max: int = 2,
    Z: float = 1.0,
    n_Re: int = 400,
    R_e_min: float = 0.3,
    R_e_max: float = 15.0,
    verbose: bool = True,
    n_basis: int = 10,
    n_quad: int = 100,
    cusp_gamma: Optional[float] = None,
    radial_method: str = 'fd',
    n_basis_radial: int = 25,
    alpha_radial: float = 1.0,
    symmetry: str = 'singlet',
) -> dict:
    """Full Level 4 H2 solver with Kato cusp factor.

    Computes D_e for H2 using the cusp-modified spectral angular basis.
    Uses the same radial solver infrastructure as solve_level4_h2_multichannel.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    l_max : int
        Maximum angular momentum per electron.
    cusp_gamma : float or None
        Cusp factor strength. None = Kato (R_e/2), 0 = baseline.
    verbose : bool
        Print diagnostics.

    Returns
    -------
    result : dict
        Keys: E_elec, E_total, D_e, D_e_pct, etc.
    """
    import time
    t0 = time.time()

    Z_A = Z_B = Z

    # Non-uniform R_e grid
    R_e_angular = np.concatenate([
        np.linspace(R_e_min, 1.0, 40),
        np.linspace(1.0, 3.0, 40),
        np.linspace(3.0, 6.0, 30),
        np.linspace(6.0, R_e_max, 20),
    ])
    R_e_angular = np.unique(R_e_angular)

    channels = _channel_list(l_max, homonuclear=True)
    n_ch = len(channels)

    if verbose:
        gamma_str = f"R_e/2 (Kato)" if cusp_gamma is None else f"{cusp_gamma}"
        print(f"Level 4 cusp-factor solver for H2 at R = {R:.4f} bohr")
        print(f"  l_max={l_max}, n_ch={n_ch}, cusp_gamma={gamma_str}")
        print(f"  n_basis={n_basis}, n_quad={n_quad}")

    U_angular = compute_adiabatic_curve_cusp(
        R, R_e_angular, l_max, Z,
        n_basis=n_basis, n_quad=n_quad,
        cusp_gamma=cusp_gamma,
        symmetry=symmetry,
    )

    t1 = time.time()
    if verbose:
        i_min = np.argmin(U_angular)
        print(f"  Angular sweep: {t1 - t0:.2f}s")
        print(f"  U_min = {U_angular[i_min]:.6f} Ha at R_e = {R_e_angular[i_min]:.3f}")

    # Interpolate and solve radial equation
    U_spline = CubicSpline(R_e_angular, U_angular, extrapolate=True)

    if radial_method == 'spectral':
        from geovac.level4_multichannel import solve_adiabatic_radial_spectral
        E_elec, F, _ = solve_adiabatic_radial_spectral(
            U_spline, n_basis=n_basis_radial,
            alpha=alpha_radial, R_e_min=R_e_min,
            R_e_max=R_e_max,
        )
    else:
        h_Re = (R_e_max - R_e_min) / (n_Re + 1)
        R_e_radial = R_e_min + (np.arange(n_Re) + 1) * h_Re
        V_radial = U_spline(R_e_radial)

        diag = np.ones(n_Re) / h_Re**2 + V_radial
        off_diag = -0.5 * np.ones(n_Re - 1) / h_Re**2

        evals, evecs = eigh_tridiagonal(
            diag, off_diag,
            select='i', select_range=(0, 0),
        )

        E_elec = evals[0]
        F = evecs[:, 0]

    t2 = time.time()

    norm = np.sqrt(np.sum(F**2))
    if norm > 0:
        F /= norm

    V_NN = Z**2 / R
    E_total = E_elec + V_NN
    E_atoms = -Z**2  # 2 * (-Z^2/2)
    D_e = E_atoms - E_total

    E_exact = -1.17447
    D_e_exact = 0.17447
    D_e_pct = (D_e / D_e_exact * 100)

    if verbose:
        print(f"\n  === Results (l_max={l_max}, cusp factor) ===")
        print(f"  E_total     = {E_total:.6f} Ha  (exact: {E_exact:.6f})")
        print(f"  D_e         = {D_e:.6f} Ha  (exact: {D_e_exact:.6f})")
        print(f"  D_e / exact = {D_e_pct:.1f}%")
        print(f"  Total time: {t2 - t0:.2f}s")

    return {
        'E_elec': E_elec,
        'E_total': E_total,
        'D_e': D_e,
        'D_e_pct': D_e_pct,
        'R': R,
        'l_max': l_max,
        'n_ch': n_ch,
        'channels': channels,
        'cusp_gamma': cusp_gamma,
        'E_atoms': E_atoms,
        'V_NN': V_NN,
        'time_angular': t1 - t0,
        'time_total': t2 - t0,
    }
