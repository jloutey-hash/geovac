"""
Algebraic angular eigenvalue solver for Level 3 hyperspherical method.

Replaces the FD-based solve_angular() with a spectral basis expansion.
Matrix dimension ~ n_l * n_basis instead of ~ n_alpha * n_l (orders of
magnitude smaller for typical parameters).

Each l-channel uses basis functions (after Liouville substitution):
    u_{l,k}(alpha) = N (sin alpha cos alpha)^{l+1} C_k^{l+1}(cos 2 alpha)

where C_k^lambda is the Gegenbauer polynomial and k = 0,2,4,... (singlet)
or k = 1,3,5,... (triplet).  For l=0 this reduces to sin(2*n*alpha) with
n = k+1 via the Chebyshev U_k identity:
    sin(alpha) cos(alpha) * U_k(cos 2alpha) = sin(2(k+1) alpha) / 2

The free eigenvalues are mu_free(l,k) = 2(l+k+1)^2 - 2, from the SO(6)
Casimir nu(nu+4)/2 with nu = 2(l+k).

The centrifugal singularity l(l+1)/cos^2(alpha) that plagues the FD solver
is handled analytically by the basis functions (they are eigenfunctions of
the centrifugal operator).  This eliminates the l_max degradation in
Paper 13 Table I.

References:
  - Paper 13, Section III (angular equation), Section XII (algebraic structure)
  - Paper 12 (Neumann expansion precedent)
"""

import numpy as np
from scipy.linalg import eigh
from scipy.integrate import quad as scipy_quad
from typing import Tuple, Optional
from math import sqrt, pi

from geovac.hyperspherical_angular import gaunt_integral, _precompute_gaunt


def _gegenbauer(k: int, lam: float, x: np.ndarray) -> np.ndarray:
    """Evaluate Gegenbauer polynomial C_k^lambda(x) via stable recurrence.

    C_0^lam(x) = 1
    C_1^lam(x) = 2*lam*x
    C_{n+1}^lam(x) = (2(n+lam)x C_n - (n+2lam-1) C_{n-1}) / (n+1)

    Parameters
    ----------
    k : int
        Polynomial degree.
    lam : float
        Gegenbauer parameter (lambda = l+1 for channel l).
    x : ndarray
        Evaluation points.

    Returns
    -------
    ndarray
        C_k^lambda(x) at each point.
    """
    if k == 0:
        return np.ones_like(x, dtype=float)
    c_prev = np.ones_like(x, dtype=float)
    c_curr = 2.0 * lam * x
    if k == 1:
        return c_curr
    for n in range(1, k):
        c_next = (2.0 * (n + lam) * x * c_curr
                  - (n + 2 * lam - 1) * c_prev) / (n + 1)
        c_prev = c_curr
        c_curr = c_next
    return c_curr


class AlgebraicAngularSolver:
    """Spectral Gegenbauer solver for the Level 3 angular eigenvalue problem.

    Replaces the FD-based solve_angular() with algebraic matrix elements
    in the free SO(6) eigenbasis.  Matrix dimension ~ n_l * n_basis instead
    of ~ n_alpha * n_l.

    For channel l, the basis functions (after Liouville substitution) are:
        u_{l,k}(alpha) = N_{l,k} (sin alpha cos alpha)^{l+1} C_k^{l+1}(cos 2 alpha)
    where k = 0, 2, 4, ... (singlet) or 1, 3, 5, ... (triplet).

    The free eigenvalues are mu_free(l,k) = 2(l+k+1)^2 - 2, corresponding
    to SO(6) Casimir nu(nu+4)/2 with nu = 2(l+k).

    Parameters
    ----------
    Z : float
        Nuclear charge.
    n_basis : int
        Number of basis functions per l-channel.
    l_max : int
        Maximum partial wave (l1=l2=l for L=0 singlet).
    symmetry : str
        'singlet' selects even-k modes (exchange symmetric),
        'triplet' selects odd-k modes (exchange antisymmetric).
    n_quad : int
        Number of Gauss-Legendre quadrature points per sub-interval.
        The integration domain [0, pi/2] is split at pi/4 (where
        max(cos, sin) has a kink), so total quadrature points = 2 * n_quad.
    """

    def __init__(
        self,
        Z: float,
        n_basis: int = 10,
        l_max: int = 0,
        symmetry: str = 'singlet',
        n_quad: int = 100,
    ) -> None:
        self.Z = Z
        self.n_basis = n_basis
        self.l_max = l_max
        self.symmetry = symmetry
        self.n_quad = n_quad
        self.n_l = l_max + 1
        self._total_dim = self.n_l * n_basis

        # Backward-compatible basis_indices for l=0 channel
        if symmetry == 'singlet':
            self.basis_indices = np.array([2 * j + 1 for j in range(n_basis)])
        elif symmetry == 'triplet':
            self.basis_indices = np.array([2 * j + 2 for j in range(n_basis)])
        else:
            raise ValueError(f"Unknown symmetry: {symmetry}")

        # SO(6) Casimir eigenvalues for l=0 channel (backward compat)
        self._casimir_eigenvalues = (
            2.0 * self.basis_indices.astype(float) ** 2 - 2.0
        )

        self._setup_quadrature(n_quad)
        self._build_channels()
        self._precompute_coupling()

    def _setup_quadrature(self, n_quad: int) -> None:
        """Setup Gauss-Legendre quadrature on [0, pi/4] and [pi/4, pi/2].

        Splitting at pi/4 handles the kink in max(cos, sin) and clusters
        quadrature nodes near the boundary singularities (1/sin near 0,
        1/cos near pi/2), improving accuracy.
        """
        from numpy.polynomial.legendre import leggauss

        nodes, weights = leggauss(n_quad)

        # Map to [0, pi/4]
        a1, b1 = 0.0, np.pi / 4.0
        scale1 = (b1 - a1) / 2.0
        alpha1 = scale1 * nodes + (b1 + a1) / 2.0
        w1 = weights * scale1

        # Map to [pi/4, pi/2]
        a2, b2 = np.pi / 4.0, np.pi / 2.0
        scale2 = (b2 - a2) / 2.0
        alpha2 = scale2 * nodes + (b2 + a2) / 2.0
        w2 = weights * scale2

        # Concatenate
        self._alpha = np.concatenate([alpha1, alpha2])
        self._weights = np.concatenate([w1, w2])

        # Precompute trig values at quadrature points
        self._sin_a = np.sin(self._alpha)
        self._cos_a = np.cos(self._alpha)

    def _build_channels(self) -> None:
        """Build orthonormalized basis functions for each l-channel.

        Channel l uses basis functions:
            u_{l,k}(alpha) = N (sin a cos a)^{l+1} C_k^{l+1}(cos 2a)
        with k = 0,2,4,... (singlet) or 1,3,5,... (triplet).

        The envelope (sin a cos a)^{l+1} absorbs the centrifugal singularity
        l(l+1)/cos^2(a), ensuring all integrands are smooth.
        """
        self._channel_phi = []
        self._channel_casimir = []

        cos_2a = np.cos(2.0 * self._alpha)
        sin_cos = self._sin_a * self._cos_a

        for l in range(self.n_l):
            # Gegenbauer k-indices: even for singlet, odd for triplet
            if self.symmetry == 'singlet':
                k_indices = np.array([2 * j for j in range(self.n_basis)])
            else:
                k_indices = np.array([2 * j + 1 for j in range(self.n_basis)])

            # Free eigenvalues: mu = 2(l+k+1)^2 - 2
            casimir = 2.0 * (l + k_indices + 1.0) ** 2 - 2.0

            # Evaluate and normalize basis functions
            envelope = sin_cos ** (l + 1)
            lam = float(l + 1)

            phi = np.zeros((self.n_basis, len(self._alpha)))
            for j, kv in enumerate(k_indices):
                raw = envelope * _gegenbauer(kv, lam, cos_2a)
                norm_sq = np.dot(self._weights, raw * raw)
                phi[j] = raw / np.sqrt(norm_sq)

            self._channel_phi.append(phi)
            self._channel_casimir.append(casimir)

        # Store l=0 basis for backward compatibility
        self._phi = self._channel_phi[0]

    def _precompute_coupling(self) -> None:
        """Build R-independent nuclear and V_ee coupling matrices.

        Nuclear coupling is diagonal in l (isotropic potential).
        V_ee coupling mixes channels l <-> l' via Gaunt integrals:

            W_{ll'}(alpha) = sqrt((2l+1)(2l'+1))/2
                             * sum_k G(l,k,l') * (min/max)^k / max

        All integrands are smooth because the basis functions absorb
        the boundary singularities.
        """
        dim = self._total_dim
        nb = self.n_basis
        nuclear = np.zeros((dim, dim))
        vee = np.zeros((dim, dim))

        V_nuc_alpha = -self.Z * (1.0 / self._cos_a + 1.0 / self._sin_a)
        min_sc = np.minimum(self._sin_a, self._cos_a)
        max_sc = np.maximum(self._sin_a, self._cos_a)

        G = _precompute_gaunt(self.l_max)

        for l in range(self.n_l):
            phi_l = self._channel_phi[l]
            i0, i1 = l * nb, (l + 1) * nb

            # Nuclear: diagonal in l (potential is isotropic in theta_12)
            weighted_nuc = phi_l * (self._weights * V_nuc_alpha)[np.newaxis, :]
            nuclear[i0:i1, i0:i1] = weighted_nuc @ phi_l.T

            # V_ee: l <-> lp coupling via Gaunt integrals
            for lp in range(l, self.n_l):
                phi_lp = self._channel_phi[lp]
                j0, j1 = lp * nb, (lp + 1) * nb

                # Build W_{ll'}(alpha) = norm * sum_k G(l,k,l') * f_k/max
                W = np.zeros_like(self._alpha)
                for k in range(abs(l - lp), l + lp + 1):
                    if k > 2 * self.l_max:
                        continue
                    gv = G[l, k, lp]
                    if abs(gv) < 1e-15:
                        continue
                    f_k = (min_sc / max_sc) ** k
                    W += gv * f_k / max_sc

                norm_ll = sqrt((2 * l + 1) * (2 * lp + 1)) / 2.0
                W *= norm_ll

                weighted_vee = phi_l * (self._weights * W)[np.newaxis, :]
                block = weighted_vee @ phi_lp.T

                vee[i0:i1, j0:j1] = block
                if l != lp:
                    vee[j0:j1, i0:i1] = block.T

        self._nuclear_full = nuclear
        self._vee_full = vee
        self._coupling_full = nuclear + vee

        # Backward-compatible l=0 block views
        self._nuclear_matrix = nuclear[:nb, :nb].copy()
        self._vee_matrix = vee[:nb, :nb].copy()
        self._coupling_matrix = (nuclear[:nb, :nb] + vee[:nb, :nb]).copy()

    def solve_with_dboc(
        self, R: float, n_channels: int = 1
    ) -> Tuple[np.ndarray, np.ndarray, float]:
        """Solve angular problem and compute DBOC for the ground channel.

        The DBOC (diagonal Born-Oppenheimer correction) accounts for the
        kinetic energy cost of the angular eigenfunction changing shape
        as R varies.  Uses the Hellmann-Feynman theorem:

            P_{mu,0}(R) = <Phi_mu|V_coupling|Phi_0> / (mu_0 - mu_mu)

        where V_coupling = dH/dR is the precomputed R-independent coupling
        matrix.  The DBOC is then:

            DBOC(R) = (1/2) sum_{mu != 0} |P_{mu,0}(R)|^2

        This is a positive (repulsive) correction to V_eff(R), raising the
        energy toward the exact value.

        Parameters
        ----------
        R : float
            Hyperradius (bohr).
        n_channels : int
            Number of eigenvalues/eigenvectors to return.

        Returns
        -------
        eigenvalues : ndarray of shape (n_channels,)
            Angular eigenvalues mu(R), sorted ascending.
        eigenvectors : ndarray of shape (n_channels, total_dim)
            Expansion coefficients in the multichannel basis.
        dboc : float
            Diagonal Born-Oppenheimer correction at this R (Hartree).
        """
        casimir_all = np.concatenate(self._channel_casimir)
        H = np.diag(casimir_all) + R * self._coupling_full
        evals, evecs = eigh(H)

        # V_coupling applied to ground state eigenvector
        phi_0 = evecs[:, 0]
        V_phi_0 = self._coupling_full @ phi_0

        dboc = 0.0
        for mu in range(1, len(evals)):
            numerator = evecs[:, mu] @ V_phi_0
            denominator = evals[mu] - evals[0]
            if abs(denominator) > 1e-12:
                P_mu0 = numerator / denominator
                dboc += P_mu0 ** 2

        dboc *= 0.5

        return evals[:n_channels], evecs[:, :n_channels].T, dboc

    def solve(
        self, R: float, n_channels: int = 1
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Solve the angular eigenvalue problem at hyperradius R.

        H_ang(R) = diag(casimir_all) + R * V_coupling

        The coupling matrix V_coupling is precomputed once (R-independent);
        only the R scaling and small-matrix diagonalization happen at each R.

        Parameters
        ----------
        R : float
            Hyperradius (bohr).
        n_channels : int
            Number of eigenvalues/eigenvectors to return.

        Returns
        -------
        eigenvalues : ndarray of shape (n_channels,)
            Angular eigenvalues mu(R), sorted ascending.
        eigenvectors : ndarray of shape (n_channels, total_dim)
            Expansion coefficients in the multichannel basis.
        """
        casimir_all = np.concatenate(self._channel_casimir)
        H = np.diag(casimir_all) + R * self._coupling_full
        evals, evecs = eigh(H)
        return evals[:n_channels], evecs[:, :n_channels].T

    def coupling_integral_quad(self, ni: int, nj: int) -> Tuple[float, float]:
        """Compute coupling matrix element via scipy adaptive quadrature.

        For verification against the GL quadrature.  Uses the l=0
        sin(2*n*alpha) form directly.

        Parameters
        ----------
        ni, nj : int
            Basis function indices (the n in sin(2*n*alpha)).

        Returns
        -------
        nuc_val : float
            Nuclear coupling integral.
        vee_val : float
            V_ee monopole coupling integral.
        """
        Z = self.Z
        norm = 4.0 / np.pi

        def nuc_integrand(alpha: float) -> float:
            s_i = np.sin(2 * ni * alpha)
            s_j = np.sin(2 * nj * alpha)
            return norm * s_i * s_j * (-Z) * (
                1.0 / np.cos(alpha) + 1.0 / np.sin(alpha)
            )

        def vee_integrand_lo(alpha: float) -> float:
            """V_ee on [0, pi/4] where max = cos."""
            s_i = np.sin(2 * ni * alpha)
            s_j = np.sin(2 * nj * alpha)
            return norm * s_i * s_j / np.cos(alpha)

        def vee_integrand_hi(alpha: float) -> float:
            """V_ee on [pi/4, pi/2] where max = sin."""
            s_i = np.sin(2 * ni * alpha)
            s_j = np.sin(2 * nj * alpha)
            return norm * s_i * s_j / np.sin(alpha)

        nuc_val, nuc_err = scipy_quad(
            nuc_integrand, 0, np.pi / 2,
            limit=400, points=[np.pi / 4],
        )
        vee_lo, _ = scipy_quad(vee_integrand_lo, 0, np.pi / 4, limit=200)
        vee_hi, _ = scipy_quad(vee_integrand_hi, np.pi / 4, np.pi / 2, limit=200)
        vee_val = vee_lo + vee_hi

        return nuc_val, vee_val

    @property
    def nuclear_matrix(self) -> np.ndarray:
        """Nuclear coupling matrix for l=0 channel (backward compat)."""
        return self._nuclear_matrix.copy()

    @property
    def vee_matrix(self) -> np.ndarray:
        """V_ee coupling matrix for l=0 channel (backward compat)."""
        return self._vee_matrix.copy()

    @property
    def coupling_matrix(self) -> np.ndarray:
        """Total coupling matrix for l=0 channel (backward compat)."""
        return self._coupling_matrix.copy()

    @property
    def casimir_eigenvalues(self) -> np.ndarray:
        """SO(6) Casimir eigenvalues for l=0 channel (backward compat)."""
        return self._casimir_eigenvalues.copy()

    @property
    def nuclear_matrix_full(self) -> np.ndarray:
        """Full nuclear coupling matrix (all channels)."""
        return self._nuclear_full.copy()

    @property
    def vee_matrix_full(self) -> np.ndarray:
        """Full V_ee coupling matrix (all channels)."""
        return self._vee_full.copy()

    @property
    def coupling_matrix_full(self) -> np.ndarray:
        """Full coupling matrix (all channels)."""
        return self._coupling_full.copy()


def perturbation_series_mu(
    H0_diag: np.ndarray,
    V_C: np.ndarray,
    n_channels: int,
    max_order: int = 10,
) -> np.ndarray:
    """Compute Rayleigh-Schrödinger perturbation coefficients for mu(R).

    The angular Hamiltonian is H(R) = diag(H0) + R * V_C (linear matrix
    pencil).  The eigenvalues mu_nu(R) can be expanded as a power series
    in R:
        mu_nu(R) = a_0 + a_1 * R + a_2 * R^2 + ... + a_k * R^k

    where a_0 = H0_diag[nu] (SO(6) Casimir) and the higher-order
    coefficients are computed via standard RS perturbation theory.

    Parameters
    ----------
    H0_diag : ndarray of shape (dim,)
        Zeroth-order eigenvalues (SO(6) Casimir values).
    V_C : ndarray of shape (dim, dim)
        R-independent coupling matrix (nuclear + V_ee).
    n_channels : int
        Number of channels to compute perturbation series for.
    max_order : int
        Maximum perturbation order (default 10).

    Returns
    -------
    coeffs : ndarray of shape (n_channels, max_order + 1)
        coeffs[nu, k] = a_k for channel nu.
        coeffs[nu, 0] = H0_diag[nu] (zeroth order).
        coeffs[nu, 1] = <nu|V_C|nu> (first order).
    """
    dim = len(H0_diag)

    # Zeroth-order eigenvectors are identity (H0 is diagonal)
    # We compute the RS correction vectors |psi_nu^(k)> and energies a_k

    coeffs = np.zeros((n_channels, max_order + 1))

    for nu in range(n_channels):
        # a_0 = unperturbed eigenvalue
        E0 = H0_diag[nu]
        coeffs[nu, 0] = E0

        # First order: a_1 = V_nn
        coeffs[nu, 1] = V_C[nu, nu]

        # Build RS correction vectors iteratively
        # |psi^(0)> = e_nu (unit vector)
        # |psi^(k)> = sum_{m != nu} c^(k)_m |m>
        # where c^(k)_m = [sum_{j=1}^{k} a_j * c^(k-j)_m
        #                  - sum_{m'} V_{m,m'} c^(k-1)_{m'}] / (E0 - E_m)
        # but c^(k-1)_{m'} only involves m' != nu

        # Store correction vector coefficients (excluding nu component)
        # psi_k[m] = coefficient of |m> in |psi^(k)> for m != nu
        psi_vectors = []  # psi_vectors[k] is dict/array for order k

        # |psi^(0)> = e_nu → coefficients are all zero except nu
        psi_0 = np.zeros(dim)
        psi_0[nu] = 1.0
        psi_vectors.append(psi_0)

        # First-order wavefunction correction
        psi_1 = np.zeros(dim)
        for m in range(dim):
            if m == nu:
                continue
            gap = E0 - H0_diag[m]
            if abs(gap) < 1e-15:
                continue
            psi_1[m] = V_C[m, nu] / gap
        psi_vectors.append(psi_1)

        # Higher orders via RS recursion
        for k in range(2, max_order + 1):
            # a_k = <psi^(0)|V|psi^(k-1)> - sum_{j=1}^{k-1} a_j * <psi^(0)|psi^(k-j)>
            # Since |psi^(0)> = e_nu:
            # a_k = V[nu,:] @ psi^(k-1) - sum_{j=1}^{k-1} a_j * psi^(k-j)[nu]
            # But psi^(j)[nu] = 0 for j >= 1 (intermediate normalization)
            # So: a_k = sum_m V[nu,m] * psi^(k-1)[m] = V_C[nu,:] @ psi_{k-1}
            a_k = V_C[nu, :] @ psi_vectors[k - 1]
            coeffs[nu, k] = a_k

            # |psi^(k)> correction vector
            # For m != nu:
            # psi^(k)[m] = (1/(E0 - E_m)) * [V_C[m,:] @ psi^(k-1)
            #               - sum_{j=1}^{k} a_j * psi^(k-j)[m]]
            psi_k = np.zeros(dim)
            V_psi_prev = V_C @ psi_vectors[k - 1]

            for m in range(dim):
                if m == nu:
                    continue
                gap = E0 - H0_diag[m]
                if abs(gap) < 1e-15:
                    continue
                val = V_psi_prev[m]
                for j in range(1, k + 1):
                    val -= coeffs[nu, j] * psi_vectors[k - j][m]
                psi_k[m] = val / gap

            psi_vectors.append(psi_k)

    return coeffs


def evaluate_perturbation_series(
    coeffs: np.ndarray,
    R: np.ndarray,
) -> np.ndarray:
    """Evaluate mu(R) from perturbation series coefficients.

    Parameters
    ----------
    coeffs : ndarray of shape (n_channels, max_order + 1)
        Perturbation coefficients from perturbation_series_mu().
    R : ndarray of shape (N_R,)
        Hyperradius values.

    Returns
    -------
    mu : ndarray of shape (n_channels, N_R)
        mu_nu(R) evaluated via the truncated power series.
    """
    n_channels, n_terms = coeffs.shape
    R = np.asarray(R, dtype=float)
    mu = np.zeros((n_channels, len(R)))

    for nu in range(n_channels):
        # Horner's method for numerical stability
        result = np.full_like(R, coeffs[nu, n_terms - 1])
        for k in range(n_terms - 2, -1, -1):
            result = result * R + coeffs[nu, k]
        mu[nu] = result

    return mu


def pade_approximant(
    coeffs_1d: np.ndarray,
    p_order: int,
    q_order: int,
) -> tuple:
    """Construct [p/q] Padé approximant from power series coefficients.

    Given a_0 + a_1*x + a_2*x^2 + ..., find polynomials P(x) and Q(x)
    with deg(P) = p_order, deg(Q) = q_order, Q(0) = 1, such that
    P(x)/Q(x) matches the power series through order p + q.

    Parameters
    ----------
    coeffs_1d : ndarray of shape (n_terms,)
        Power series coefficients [a_0, a_1, ..., a_{n-1}].
    p_order : int
        Degree of numerator polynomial.
    q_order : int
        Degree of denominator polynomial.

    Returns
    -------
    p_coeffs : ndarray of shape (p_order + 1,)
        Numerator polynomial coefficients [p_0, p_1, ..., p_p].
    q_coeffs : ndarray of shape (q_order + 1,)
        Denominator polynomial coefficients [1, q_1, ..., q_q].

    Raises
    ------
    ValueError
        If not enough series coefficients or singular system.
    """
    n_needed = p_order + q_order + 1
    if len(coeffs_1d) < n_needed:
        raise ValueError(
            f"Need {n_needed} coefficients for [{p_order}/{q_order}] Padé, "
            f"have {len(coeffs_1d)}"
        )

    c = coeffs_1d[:n_needed]

    # Solve for q_coeffs from the linear system:
    # sum_{j=0}^{q} q_j * c_{p+1+i-j} = 0, i = 0, ..., q-1
    # with q_0 = 1.
    if q_order == 0:
        return c[:p_order + 1].copy(), np.array([1.0])

    # Build the linear system for q_1, ..., q_q
    A = np.zeros((q_order, q_order))
    b = np.zeros(q_order)
    for i in range(q_order):
        row_idx = p_order + 1 + i  # index into c
        b[i] = -c[row_idx] if row_idx < n_needed else 0.0
        for j in range(q_order):
            ci = row_idx - (j + 1)
            if 0 <= ci < n_needed:
                A[i, j] = c[ci]

    try:
        q_tail = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        raise ValueError("Singular Padé system")

    q_coeffs = np.zeros(q_order + 1)
    q_coeffs[0] = 1.0
    q_coeffs[1:] = q_tail

    # Compute p_coeffs: p_k = sum_{j=0}^{min(k, q)} q_j * c_{k-j}
    p_coeffs = np.zeros(p_order + 1)
    for k in range(p_order + 1):
        for j in range(min(k, q_order) + 1):
            p_coeffs[k] += q_coeffs[j] * c[k - j]

    return p_coeffs, q_coeffs


def evaluate_pade(
    p_coeffs: np.ndarray,
    q_coeffs: np.ndarray,
    R: np.ndarray,
) -> np.ndarray:
    """Evaluate a Padé approximant P(R)/Q(R).

    Parameters
    ----------
    p_coeffs, q_coeffs : ndarray
        Polynomial coefficients from pade_approximant().
    R : ndarray
        Evaluation points.

    Returns
    -------
    result : ndarray
        P(R)/Q(R) at each point.
    """
    R = np.asarray(R, dtype=float)
    # Horner evaluation
    p_val = np.full_like(R, p_coeffs[-1])
    for k in range(len(p_coeffs) - 2, -1, -1):
        p_val = p_val * R + p_coeffs[k]

    q_val = np.full_like(R, q_coeffs[-1])
    for k in range(len(q_coeffs) - 2, -1, -1):
        q_val = q_val * R + q_coeffs[k]

    return p_val / q_val


def perturbation_series_P_matrix(
    H0_diag: np.ndarray,
    V_C: np.ndarray,
    n_channels: int,
    max_order: int = 10,
) -> np.ndarray:
    """Compute the first-derivative coupling P_μν(R) as a series in R.

    The Hellmann-Feynman P-matrix is:
        P_μν(R) = <Φ_μ(R)|V_C|Φ_ν(R)> / (μ_ν(R) - μ_μ(R))

    Both numerator and denominator are power series in R (from the
    perturbation expansion of eigenvectors and eigenvalues).  The ratio
    is itself a power series, computable order by order.

    For the numerator N_μν(R) = <Φ_μ(R)|V_C|Φ_ν(R)>:
      N^(0) = V_C[μ,ν]  (zeroth-order eigenvectors are Kronecker deltas)
      N^(k) = sum over correction vectors

    For the denominator D_μν(R) = μ_ν(R) - μ_μ(R):
      D^(0) = E0_ν - E0_μ
      D^(k) = a_k(ν) - a_k(μ)

    P = N/D is expanded via the standard series division formula.

    Parameters
    ----------
    H0_diag : ndarray of shape (dim,)
        Zeroth-order eigenvalues.
    V_C : ndarray of shape (dim, dim)
        R-independent coupling matrix.
    n_channels : int
        Number of channels.
    max_order : int
        Maximum order for the P series.

    Returns
    -------
    P_coeffs : ndarray of shape (n_channels, n_channels, max_order + 1)
        P_coeffs[mu, nu, k] = coefficient of R^k in P_μν(R).
        Diagonal entries are zero by convention.
    """
    dim = len(H0_diag)

    # First get eigenvalue coefficients and eigenvector corrections
    mu_coeffs = perturbation_series_mu(H0_diag, V_C, dim, max_order)

    # Build eigenvector correction vectors for all channels
    # psi_vectors[nu][k] = correction vector at order k for channel nu
    all_psi = []
    for nu in range(dim):
        E0 = H0_diag[nu]
        psi_vecs = []

        psi_0 = np.zeros(dim)
        psi_0[nu] = 1.0
        psi_vecs.append(psi_0)

        psi_1 = np.zeros(dim)
        for m in range(dim):
            if m == nu:
                continue
            gap = E0 - H0_diag[m]
            if abs(gap) < 1e-15:
                continue
            psi_1[m] = V_C[m, nu] / gap
        psi_vecs.append(psi_1)

        for k in range(2, max_order + 1):
            V_psi_prev = V_C @ psi_vecs[k - 1]
            psi_k = np.zeros(dim)
            for m in range(dim):
                if m == nu:
                    continue
                gap = E0 - H0_diag[m]
                if abs(gap) < 1e-15:
                    continue
                val = V_psi_prev[m]
                for j in range(1, k + 1):
                    val -= mu_coeffs[nu, j] * psi_vecs[k - j][m]
                psi_k[m] = val / gap
            psi_vecs.append(psi_k)
        all_psi.append(psi_vecs)

    # Compute numerator N_μν^(k) = sum_{j=0}^{k} <psi_μ^(j)|V_C|psi_ν^(k-j)>
    N_coeffs = np.zeros((n_channels, n_channels, max_order + 1))
    for mu_idx in range(n_channels):
        for nu_idx in range(n_channels):
            if mu_idx == nu_idx:
                continue
            for k in range(max_order + 1):
                val = 0.0
                for j in range(k + 1):
                    val += all_psi[mu_idx][j] @ V_C @ all_psi[nu_idx][k - j]
                N_coeffs[mu_idx, nu_idx, k] = val

    # Denominator D_μν^(k) = mu_coeffs[nu, k] - mu_coeffs[mu, k]
    D_coeffs = np.zeros((n_channels, n_channels, max_order + 1))
    for mu_idx in range(n_channels):
        for nu_idx in range(n_channels):
            if mu_idx == nu_idx:
                continue
            D_coeffs[mu_idx, nu_idx, :] = (
                mu_coeffs[nu_idx, :max_order + 1]
                - mu_coeffs[mu_idx, :max_order + 1]
            )

    # P = N / D via series division: P^(k) = (N^(k) - sum_{j=0}^{k-1} P^(j) D^(k-j)) / D^(0)
    P_coeffs = np.zeros((n_channels, n_channels, max_order + 1))
    for mu_idx in range(n_channels):
        for nu_idx in range(n_channels):
            if mu_idx == nu_idx:
                continue
            D0 = D_coeffs[mu_idx, nu_idx, 0]
            if abs(D0) < 1e-15:
                continue
            for k in range(max_order + 1):
                val = N_coeffs[mu_idx, nu_idx, k]
                for j in range(k):
                    val -= P_coeffs[mu_idx, nu_idx, j] * D_coeffs[mu_idx, nu_idx, k - j]
                P_coeffs[mu_idx, nu_idx, k] = val / D0

    return P_coeffs


def solve_hyperspherical_algebraic(
    Z: float = 2.0,
    n_electrons: int = 2,
    n_basis: int = 10,
    l_max: int = 0,
    n_R: int = 200,
    R_min: float = 0.1,
    R_max: float = 30.0,
    N_R_radial: int = 2000,
    include_dboc: bool = False,
    verbose: bool = True,
) -> dict:
    """Full Level 3 solver using algebraic angular matrix elements.

    Drop-in replacement for the existing hyperspherical solver,
    but using the AlgebraicAngularSolver for the angular eigenvalues.

    Parameters
    ----------
    Z : float
        Nuclear charge.
    n_electrons : int
        Number of electrons (must be 2).
    n_basis : int
        Number of Gegenbauer basis functions for the angular solver.
    l_max : int
        Maximum partial wave.
    n_R : int
        Number of R points for computing mu(R).
    R_min, R_max : float
        Hyperradial grid boundaries.
    N_R_radial : int
        Number of grid points for the radial solve.
    include_dboc : bool
        If True, include the diagonal Born-Oppenheimer correction (DBOC)
        in the effective potential.  The DBOC is a positive (repulsive)
        correction that accounts for the kinetic energy cost of the angular
        eigenfunction changing shape as R varies.
    verbose : bool
        Print progress information.

    Returns
    -------
    result : dict
        Keys: 'energy', 'R_grid_angular', 'mu_curve', 'V_eff',
        'R_grid_radial', 'wavefunction', 'solver', 'Z', 'n_basis', 'l_max'.
        If include_dboc: also 'dboc_curve', 'energy_no_dboc'.
    """
    import time
    from scipy.interpolate import CubicSpline
    from geovac.hyperspherical_adiabatic import effective_potential
    from geovac.hyperspherical_radial import solve_radial

    if n_electrons != 2:
        raise ValueError("Only two-electron systems are supported.")

    t0 = time.time()

    solver = AlgebraicAngularSolver(Z, n_basis, l_max)

    # R grid: denser near origin where mu(R) changes rapidly
    R_grid = np.concatenate([
        np.linspace(R_min, 1.0, n_R // 3),
        np.linspace(1.0, 5.0, n_R // 3),
        np.linspace(5.0, R_max, n_R // 3 + 1),
    ])
    R_grid = np.unique(R_grid)

    dboc_mode = "with DBOC" if include_dboc else "adiabatic"
    if verbose:
        print(f"Algebraic angular solver: Z={Z}, n_basis={n_basis}, "
              f"l_max={l_max} ({dboc_mode})")
        print(f"Computing mu(R) on {len(R_grid)} R points...")

    # Solve angular problem at each R
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
        if include_dboc:
            print(f"  DBOC peak: {np.max(dboc_values):.6f} Ha")

    # Build V_eff and solve hyperradial equation
    V_eff = effective_potential(R_grid, mu)
    if include_dboc:
        # DBOC is added directly to V_eff (both in Hartree).
        # P_{mu0} = <Phi_mu|dPhi_0/dR> has units 1/bohr,
        # so |P|^2 has units 1/bohr^2 = Hartree in atomic units.
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
        'l_max': l_max,
        'include_dboc': include_dboc,
    }

    if include_dboc:
        # Also solve without DBOC to measure the correction
        V_eff_spline_nocorr = CubicSpline(
            R_grid, V_eff_uncorrected, extrapolate=True
        )
        E_no, _, _ = solve_radial(
            V_eff_spline_nocorr, R_min, R_max, N_R_radial, n_states=1
        )
        result['dboc_curve'] = dboc_values
        result['energy_no_dboc'] = E_no[0]

    return result
