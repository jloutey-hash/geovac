"""
Cusp-adapted angular basis for Level 4 molecule-frame hyperspherical method.

Investigates whether the Track J algebraicization pattern (absorbing singularities
into the basis weight function) can be applied to the theta_12 coordinate at Level 4
to accelerate partial-wave (l_max) convergence.

BACKGROUND:
Track U proved that an alpha-only cusp factor does not help because the cusp
is 2D in (alpha, theta_12). Track W proved that 1/r_12 cannot be absorbed into
a graph Laplacian on S^5 (dimensionality mismatch). This module investigates
whether the CHANNEL STRUCTURE (which encodes the theta_12 expansion) can be
modified to better represent the cusp.

THE KEY INSIGHT:
In Level 4, the theta_12 dependence is NOT a separate coordinate with its own
basis. It is encoded in the CHANNEL structure: each channel (l1, l2) corresponds
to a term in the coupled spherical harmonic expansion, and the inter-electron
angle theta_12 appears through the Gaunt integral coupling between channels.
Specifically:
  Psi(alpha, Omega1, Omega2) = sum_{l1,l2} F_{l1,l2}(alpha) Y_{l1}(Omega1) Y_{l2}(Omega2)
and the expansion in theta_12 = angle(r1, r2) is obtained by re-coupling the
Y_{l1} Y_{l2} product into Legendre polynomials P_l(cos theta_12).

THE TRACK J ANALOGY:
At Level 2, Track J absorbed the centrifugal 1/x singularity into an associated
Laguerre weight x^|m| e^{-x}. The analog here would be to absorb the 1/r_12
cusp into a modified weight in theta_12. But since theta_12 is encoded in the
channel structure, this means modifying how channels are combined.

APPROACH: Cusp-weighted channel extrapolation
Instead of modifying the basis functions, we use the KNOWN asymptotic form of the
channel expansion coefficients to extrapolate the l_max -> infinity limit. The
Schwartz (1962) result tells us that the energy contribution from channel l goes as:
  Delta_E_l ~ A / (l + 1/2)^4
This is the SAME physics as the cusp correction (Track X), but applied at the
angular eigenvalue level rather than the total energy.

We implement two approaches:
A. Channel coefficient analysis: extract and analyze the channel expansion
   coefficients of the angular eigenfunction to diagnose the cusp convergence.
B. Extrapolated angular eigenvalue: use the Schwartz partial-wave formula
   to extrapolate mu(rho) from finite l_max to the l_max -> infinity limit,
   applied per-rho-point to the angular eigenvalue rather than post-hoc
   to the total energy.

RESULT: This investigation demonstrates that the Track J pattern does NOT
directly transfer to the theta_12 channel structure because:
1. The channel index (l1, l2) is a COUPLED angular momentum, not a single
   coordinate with a weight function.
2. The cusp lives in the 2D coalescence manifold (alpha, theta_12), and
   no 1D weight modification can capture it.
3. The correct approach is the Schwartz extrapolation (Track X), which is
   a post-processing correction rather than a basis modification.

Author: Claude (Track Y investigation, v2.0.14 sprint)
"""

import numpy as np
from scipy.linalg import eigh
from scipy.special import eval_jacobi
from typing import Tuple, List, Dict, Optional, Union
from math import sqrt, pi

from geovac.hyperspherical_angular import gaunt_integral
from geovac.level4_multichannel import (
    _channel_list,
    _channel_list_extended,
    compute_nuclear_coupling,
    _ee_coupling,
    _ee_coupling_general,
)


# ============================================================================
# Part 1: Channel coefficient analysis
# ============================================================================

def extract_channel_coefficients(
    rho: float,
    R_e: float,
    l_max: int = 4,
    n_basis: int = 10,
    n_quad: int = 100,
    Z: float = 1.0,
    symmetry: str = 'singlet',
) -> Dict:
    """Extract the channel expansion coefficients of the ground-state eigenfunction.

    The ground-state angular eigenfunction is:
        Psi(alpha) = sum_{channels} sum_{k} c_{ch,k} phi_{ch,k}(alpha)

    The "channel weight" is the norm-squared of coefficients in each channel:
        W_{ch} = sum_k |c_{ch,k}|^2

    This weight distribution reveals how the cusp manifests in the channel
    expansion: slow l_max convergence means the high-l channels carry
    significant weight.

    Parameters
    ----------
    rho : float
        R / (2 R_e).
    R_e : float
        Electronic hyperradius.
    l_max : int
        Maximum angular momentum per electron.
    n_basis : int
        Number of Jacobi basis functions per channel.
    n_quad : int
        Quadrature points per sub-interval.
    Z : float
        Nuclear charge.
    symmetry : str
        'singlet' or 'triplet'.

    Returns
    -------
    dict
        Channel weights, eigenvalue, and convergence diagnostics.
    """
    from geovac.level4_spectral_angular import Level4SpectralAngular

    solver = Level4SpectralAngular(
        l_max=l_max, n_basis=n_basis, n_quad=n_quad,
        Z_A=Z, Z_B=Z, symmetry=symmetry,
    )

    # Solve and get full eigenvector
    dim = solver._total_dim
    H = np.zeros((dim, dim))

    # Free Hamiltonian
    for ic in range(solver.n_ch):
        n_k = len(solver._channel_k_indices[ic])
        i0 = solver._basis_offsets[ic]
        for j in range(n_k):
            H[i0 + j, i0 + j] = solver._channel_casimir[ic][j]

    # Nuclear + V_ee
    V_nuc = solver._compute_nuclear_spectral(rho, R_e)
    H += R_e * V_nuc
    H += R_e * solver._vee_matrix

    evals, evecs = eigh(H)
    mu_0 = evals[0]
    c_0 = evecs[:, 0]

    # Extract per-channel weights
    channels = solver.channels_2
    channel_weights = []
    channel_data = []

    for ic in range(solver.n_ch):
        l1, l2 = channels[ic]
        n_k = len(solver._channel_k_indices[ic])
        i0 = solver._basis_offsets[ic]
        c_ch = c_0[i0:i0 + n_k]
        weight = np.sum(c_ch**2)
        channel_weights.append(weight)
        channel_data.append({
            'l1': l1, 'l2': l2,
            'l_sum': l1 + l2,
            'weight': weight,
            'coefficients': c_ch.copy(),
            'n_basis': n_k,
        })

    # Sort by l_sum for convergence analysis
    channel_weights = np.array(channel_weights)

    # Group weights by l_sum = l1 + l2
    l_sum_weights = {}
    for cd in channel_data:
        ls = cd['l_sum']
        if ls not in l_sum_weights:
            l_sum_weights[ls] = 0.0
        l_sum_weights[ls] += cd['weight']

    return {
        'mu_0': mu_0,
        'channels': channels,
        'channel_data': channel_data,
        'channel_weights': channel_weights,
        'l_sum_weights': l_sum_weights,
        'total_weight': np.sum(channel_weights),
        'l_max': l_max,
        'rho': rho,
        'R_e': R_e,
    }


def analyze_channel_convergence(
    rho: float,
    R_e: float,
    l_max_values: List[int] = None,
    n_basis: int = 10,
    n_quad: int = 100,
    Z: float = 1.0,
) -> Dict:
    """Analyze how the angular eigenvalue converges with l_max.

    Computes mu_0(rho) at several l_max values and fits the convergence
    to the Schwartz 1/(l+1/2)^4 form.

    Parameters
    ----------
    rho : float
        R / (2 R_e).
    R_e : float
        Electronic hyperradius.
    l_max_values : list of int
        l_max values to test. Default [0, 2, 4, 6].
    n_basis : int
        Jacobi basis functions per channel.
    n_quad : int
        Quadrature points per sub-interval.
    Z : float
        Nuclear charge.

    Returns
    -------
    dict
        Convergence data and fitted parameters.
    """
    if l_max_values is None:
        l_max_values = [0, 2, 4]

    results = []
    for lm in l_max_values:
        data = extract_channel_coefficients(
            rho, R_e, l_max=lm, n_basis=n_basis, n_quad=n_quad, Z=Z,
        )
        results.append({
            'l_max': lm,
            'mu_0': data['mu_0'],
            'n_channels': len(data['channels']),
            'l_sum_weights': data['l_sum_weights'],
        })

    mu_values = np.array([r['mu_0'] for r in results])
    lmax_arr = np.array(l_max_values, dtype=float)

    # Fit incremental convergence: Delta_mu = mu(l_max) - mu(l_max-2)
    # Expected: Delta_mu ~ A / (l_max + 1/2)^4
    increments = []
    for i in range(1, len(results)):
        delta_mu = results[i]['mu_0'] - results[i-1]['mu_0']
        l_eff = (l_max_values[i] + l_max_values[i-1]) / 2.0
        increments.append({
            'l_max_from': l_max_values[i-1],
            'l_max_to': l_max_values[i],
            'delta_mu': delta_mu,
            'l_effective': l_eff,
        })

    return {
        'l_max_values': l_max_values,
        'mu_values': mu_values,
        'results': results,
        'increments': increments,
        'rho': rho,
        'R_e': R_e,
    }


# ============================================================================
# Part 2: Schwartz extrapolation of angular eigenvalue
# ============================================================================

def schwartz_extrapolate_mu(
    mu_values: np.ndarray,
    l_max_values: np.ndarray,
) -> Dict:
    """Extrapolate angular eigenvalue to l_max -> infinity using Schwartz formula.

    The partial-wave increment for the angular eigenvalue goes as:
        delta_mu_l ~ A_mu / (l + 1/2)^p

    where p is the Schwartz exponent (expected p=4 for Coulomb cusp).

    We fit A_mu and extrapolate mu_infinity = mu(l_max) + tail_correction.

    Parameters
    ----------
    mu_values : ndarray
        Angular eigenvalues at different l_max.
    l_max_values : ndarray
        Corresponding l_max values.

    Returns
    -------
    dict
        Extrapolated mu_infinity, fitted A_mu, and diagnostics.
    """
    if len(mu_values) < 2:
        return {
            'mu_extrapolated': mu_values[-1] if len(mu_values) > 0 else None,
            'converged': False,
            'reason': 'Need at least 2 l_max values for extrapolation',
        }

    # Compute increments
    delta_mu = np.diff(mu_values)
    l_mid = (l_max_values[:-1] + l_max_values[1:]) / 2.0

    if len(delta_mu) < 2:
        # Only 2 points: assume p=4 and fit A_mu
        p = 4.0
        A_mu = delta_mu[0] * (l_mid[0] + 0.5)**p
        tail = A_mu * _schwartz_tail(l_max_values[-1], p)
        return {
            'mu_extrapolated': mu_values[-1] + tail,
            'mu_last': mu_values[-1],
            'tail_correction': tail,
            'A_mu': A_mu,
            'p_assumed': p,
            'converged': True,
            'method': '2-point, p=4 assumed',
        }

    # 3+ points: fit both A_mu and p
    # log(|delta_mu|) = log(A_mu) - p * log(l_mid + 0.5)
    valid = np.abs(delta_mu) > 1e-15
    if np.sum(valid) < 2:
        return {
            'mu_extrapolated': mu_values[-1],
            'converged': False,
            'reason': 'Increments too small to fit',
        }

    log_delta = np.log(np.abs(delta_mu[valid]))
    log_l = np.log(l_mid[valid] + 0.5)
    coeffs = np.polyfit(log_l, log_delta, 1)
    p_fit = -coeffs[0]
    A_mu = np.exp(coeffs[1])

    # Sign of correction: delta_mu should be negative (cusp lowers energy)
    sign = np.sign(delta_mu[-1])

    tail = sign * A_mu * _schwartz_tail(l_max_values[-1], p_fit)

    return {
        'mu_extrapolated': mu_values[-1] + tail,
        'mu_last': mu_values[-1],
        'tail_correction': tail,
        'A_mu': A_mu,
        'p_fit': p_fit,
        'p_expected': 4.0,
        'sign': sign,
        'converged': True,
        'method': f'{len(mu_values)}-point fit',
    }


def _schwartz_tail(l_max: float, p: float) -> float:
    """Compute the Schwartz tail sum: sum_{l=l_max+1}^{infinity} 1/(l+1/2)^p.

    For p > 1, this converges. We compute it as zeta(p, l_max + 3/2) using
    the Hurwitz zeta function (approximated by direct summation for moderate p).

    Parameters
    ----------
    l_max : float
        Highest included l value.
    p : float
        Schwartz exponent.

    Returns
    -------
    float
        Tail sum value.
    """
    # Direct summation (converges quickly for p >= 4)
    total = 0.0
    for l in range(int(l_max) + 1, int(l_max) + 1000):
        term = 1.0 / (l + 0.5)**p
        total += term
        if term < 1e-15 * abs(total):
            break
    return total


# ============================================================================
# Part 3: Cusp-adapted spectral solver (Approach A: channel re-weighting)
# ============================================================================

class CuspAdaptedAngularSolver:
    """Spectral angular solver with cusp-adapted channel weighting.

    This solver explores whether modifying the channel combination
    can accelerate l_max convergence, following the Track J pattern
    of absorbing singularities into the basis.

    The idea: instead of the standard Jacobi basis
        phi_k^{(l1,l2)}(alpha) = N_k (cos a)^{l1+1} (sin a)^{l2+1} P_k(cos 2a)
    use a cusp-adapted basis where each basis function includes a
    coalescence weight:
        phi_k^{cusp}(alpha) = w(alpha) * phi_k^{(l1,l2)}(alpha)
    where w(alpha) captures the coalescence behavior.

    However, the critical insight is that the slow convergence is NOT
    in alpha (which is smooth) but in the CHANNEL sum (which encodes theta_12).
    The weight w(alpha) cannot accelerate the channel sum convergence because
    channels correspond to different angular momentum quantum numbers, not
    different alpha-space functions.

    What this solver DOES provide:
    1. Standard spectral eigenvalues for comparison
    2. Channel weight diagnostics showing where the cusp affects convergence
    3. Schwartz extrapolation of the angular eigenvalue per rho point

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
        'singlet' or 'triplet'.
    extrapolate : bool
        If True, apply Schwartz extrapolation to the angular eigenvalue.
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
        extrapolate: bool = False,
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
        self.extrapolate = extrapolate
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

        # Build sub-solvers for extrapolation if requested
        self._sub_solvers = {}
        if extrapolate and l_max >= 2:
            self._build_sub_solvers()

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

    def _build_basis(self) -> None:
        """Build orthonormalized Jacobi polynomial basis per channel."""
        self._channel_phi = []
        self._channel_casimir = []
        self._channel_k_indices = []
        self._basis_offsets = []

        offset = 0
        for ic, (l1, m1, l2, m2) in enumerate(self.channels_4):
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
        """Precompute V_ee coupling matrix."""
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

                self._vee_matrix[i0:i0+n_i, j0:j0+n_j] = block
                if ic != jc:
                    self._vee_matrix[j0:j0+n_j, i0:i0+n_i] = block.T

    def _compute_nuclear_spectral(
        self, rho: float, R_e: float,
    ) -> np.ndarray:
        """Compute nuclear coupling matrix in spectral basis."""
        dim = self._total_dim
        V_nuc = np.zeros((dim, dim))

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

    def _build_sub_solvers(self) -> None:
        """Build sub-solvers at lower l_max for Schwartz extrapolation."""
        # Build solvers at l_max=0 and intermediate values
        sub_lmax_values = []
        if self.l_max >= 4:
            sub_lmax_values = [0, 2]
        elif self.l_max >= 2:
            sub_lmax_values = [0]

        for lm in sub_lmax_values:
            self._sub_solvers[lm] = CuspAdaptedAngularSolver(
                l_max=lm, n_basis=self.n_basis, n_quad=self.n_quad,
                Z_A=self.Z_A, Z_B=self.Z_B, z0=self.z0,
                m_max=self.m_max, l_max_per_m=self.l_max_per_m,
                symmetry=self.symmetry, extrapolate=False,
            )

    def solve(
        self,
        rho: float,
        R_e: float,
        n_eig: int = 1,
    ) -> Tuple[np.ndarray, list]:
        """Solve the angular eigenvalue problem.

        If extrapolate=True, applies Schwartz extrapolation using
        sub-solvers at lower l_max.

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
        dim = self._total_dim
        H = np.zeros((dim, dim))

        for ic in range(self.n_ch):
            n_k = len(self._channel_k_indices[ic])
            i0 = self._basis_offsets[ic]
            for j in range(n_k):
                H[i0 + j, i0 + j] = self._channel_casimir[ic][j]

        V_nuc = self._compute_nuclear_spectral(rho, R_e)
        H += R_e * V_nuc
        H += R_e * self._vee_matrix

        evals, evecs = eigh(H)
        mu = evals[:n_eig].copy()

        if self.m_max == 0:
            channels = self.channels_2
        else:
            channels = self.channels_4

        if not self.extrapolate or not self._sub_solvers:
            return mu, channels

        # Schwartz extrapolation: collect mu_0 at multiple l_max values
        mu_series = []
        lmax_series = []

        for lm in sorted(self._sub_solvers.keys()):
            sub_mu, _ = self._sub_solvers[lm].solve(rho, R_e, n_eig=1)
            mu_series.append(sub_mu[0])
            lmax_series.append(lm)

        mu_series.append(mu[0])
        lmax_series.append(self.l_max)

        extrap = schwartz_extrapolate_mu(
            np.array(mu_series), np.array(lmax_series),
        )

        if extrap.get('converged', False):
            mu[0] = extrap['mu_extrapolated']

        return mu, channels

    def solve_with_diagnostics(
        self,
        rho: float,
        R_e: float,
    ) -> Dict:
        """Solve and return full diagnostic information.

        Returns eigenvalue, eigenvector, channel weights, and
        convergence analysis.
        """
        dim = self._total_dim
        H = np.zeros((dim, dim))

        for ic in range(self.n_ch):
            n_k = len(self._channel_k_indices[ic])
            i0 = self._basis_offsets[ic]
            for j in range(n_k):
                H[i0 + j, i0 + j] = self._channel_casimir[ic][j]

        V_nuc = self._compute_nuclear_spectral(rho, R_e)
        H += R_e * V_nuc
        H += R_e * self._vee_matrix

        evals, evecs = eigh(H)
        c_0 = evecs[:, 0]

        # Channel weight analysis
        channels = self.channels_2 if self.m_max == 0 else self.channels_4
        channel_weights = []
        for ic in range(self.n_ch):
            n_k = len(self._channel_k_indices[ic])
            i0 = self._basis_offsets[ic]
            c_ch = c_0[i0:i0 + n_k]
            weight = np.sum(c_ch**2)
            if self.m_max == 0:
                l1, l2 = self.channels_2[ic]
            else:
                l1, l2 = self.channels_4[ic][0], self.channels_4[ic][2]
            channel_weights.append({
                'l1': l1, 'l2': l2, 'l_sum': l1 + l2,
                'weight': weight,
            })

        # Group by l_sum
        l_sum_weights = {}
        for cw in channel_weights:
            ls = cw['l_sum']
            if ls not in l_sum_weights:
                l_sum_weights[ls] = 0.0
            l_sum_weights[ls] += cw['weight']

        return {
            'mu_0': evals[0],
            'eigenvalues': evals[:min(5, len(evals))],
            'channel_weights': channel_weights,
            'l_sum_weights': l_sum_weights,
            'H_symmetric': np.allclose(H, H.T, atol=1e-12),
            'total_dim': dim,
            'n_channels': self.n_ch,
        }


# ============================================================================
# Part 4: SO(6) Casimir verification
# ============================================================================

def verify_so6_casimir(l_max: int = 2, n_basis: int = 10) -> Dict:
    """Verify that the adapted solver reproduces SO(6) Casimir at rho=0.

    At rho -> 0 (no nuclear field), the angular eigenvalues should be
    exactly the SO(6) Casimir values: nu*(nu+4)/2 with nu = l1+l2+2k.

    This is the baseline sanity check.

    Returns
    -------
    dict
        Casimir verification results.
    """
    solver = CuspAdaptedAngularSolver(
        l_max=l_max, n_basis=n_basis, Z_A=1.0, Z_B=1.0,
        extrapolate=False,
    )

    # At rho=0, small R_e: free Hamiltonian dominates
    # Actually, rho=0 means R=0, so V_nuc is zero only if we set R_e large
    # Let's use R_e very large and rho very small
    R_e = 100.0
    rho = 0.001  # R = 2*R_e*rho = 0.2, but R_e*V_nuc ~ R_e * (small) -> negligible

    # Better: just check the free eigenvalues directly
    dim = solver._total_dim
    H_free = np.zeros((dim, dim))
    for ic in range(solver.n_ch):
        n_k = len(solver._channel_k_indices[ic])
        i0 = solver._basis_offsets[ic]
        for j in range(n_k):
            H_free[i0 + j, i0 + j] = solver._channel_casimir[ic][j]

    evals_free = np.sort(np.linalg.eigvalsh(H_free))

    # Expected: nu*(nu+4)/2 for nu = l1+l2+2k, k=0,1,...
    # For l_max=0, singlet: channel (0,0), k=0,2,4,...
    # nu = 0+0+2*0=0, 0+0+2*2=4, 0+0+2*4=8, ...
    # Casimir: 0, 4*8/2=16, 8*12/2=48, ...
    expected = []
    for ic in range(solver.n_ch):
        l1, l2 = solver.channels_2[ic]
        for k in solver._channel_k_indices[ic]:
            nu = l1 + l2 + 2 * k
            expected.append(nu * (nu + 4) / 2.0)
    expected = np.sort(expected)

    max_error = np.max(np.abs(evals_free[:len(expected)] - expected))

    return {
        'l_max': l_max,
        'n_eigenvalues': len(expected),
        'max_casimir_error': max_error,
        'casimir_exact': expected[:5],
        'casimir_computed': evals_free[:5],
        'pass': max_error < 1e-10,
    }


# ============================================================================
# Part 5: Comparison driver
# ============================================================================

def compare_standard_vs_adapted(
    rho: float,
    R_e: float,
    l_max_values: List[int] = None,
    n_basis: int = 10,
    n_quad: int = 100,
    Z: float = 1.0,
) -> Dict:
    """Compare standard spectral solver with cusp-adapted extrapolation.

    For each l_max, computes:
    1. Standard eigenvalue (from Level4SpectralAngular)
    2. Extrapolated eigenvalue (from CuspAdaptedAngularSolver with extrapolation)
    3. Channel weight distribution

    Parameters
    ----------
    rho, R_e : float
        Angular equation parameters.
    l_max_values : list of int
        l_max values to compare.
    n_basis, n_quad : int
        Basis and quadrature parameters.
    Z : float
        Nuclear charge.

    Returns
    -------
    dict
        Comparison results.
    """
    if l_max_values is None:
        l_max_values = [0, 2, 4]

    from geovac.level4_spectral_angular import Level4SpectralAngular

    results = []
    for lm in l_max_values:
        # Standard solver
        std_solver = Level4SpectralAngular(
            l_max=lm, n_basis=n_basis, n_quad=n_quad,
            Z_A=Z, Z_B=Z, symmetry='singlet',
        )
        std_mu, _ = std_solver.solve(rho, R_e, n_eig=1)

        # Adapted solver with extrapolation
        adapted = CuspAdaptedAngularSolver(
            l_max=lm, n_basis=n_basis, n_quad=n_quad,
            Z_A=Z, Z_B=Z, symmetry='singlet',
            extrapolate=(lm >= 2),
        )
        adp_mu, _ = adapted.solve(rho, R_e, n_eig=1)

        # Diagnostics
        diag = adapted.solve_with_diagnostics(rho, R_e)

        results.append({
            'l_max': lm,
            'mu_standard': std_mu[0],
            'mu_adapted': adp_mu[0],
            'mu_difference': adp_mu[0] - std_mu[0],
            'l_sum_weights': diag['l_sum_weights'],
            'n_channels': diag['n_channels'],
            'total_dim': diag['total_dim'],
        })

    return {
        'rho': rho,
        'R_e': R_e,
        'Z': Z,
        'results': results,
    }


# ============================================================================
# Part 6: Track J analogy assessment
# ============================================================================

def assess_track_j_analogy() -> Dict:
    """Assess whether the Track J pattern transfers to the theta_12 coordinate.

    Track J (Level 2, centrifugal singularity):
    - Singularity: m^2/x at x=0, a 1D singularity in a 1D coordinate
    - Solution: absorb into weight function x^|m| e^{-x}
    - Result: associated Laguerre L_n^{|m|}(x) basis
    - Effect: all matrix elements algebraic, faster convergence

    This investigation (Level 4, cusp singularity):
    - Singularity: 1/r_12 at r_12=0, a 2D singularity in (alpha, theta_12)
    - The theta_12 coordinate is NOT separable: it is encoded in the
      CHANNEL STRUCTURE (l1, l2) via the coupled spherical harmonic expansion
    - There is no weight function to modify because theta_12 is a derived
      quantity from the channel coupling, not an integration variable
    - The channel quantum numbers (l1, l2) are fixed by the Gaunt integral
      selection rules and cannot be continuously modified

    Returns
    -------
    dict
        Assessment of the analogy.
    """
    return {
        'track_j_singularity': {
            'type': 'centrifugal',
            'location': 'x=0 in radial coordinate',
            'dimension': '1D',
            'form': 'm^2/x',
            'coordinate': 'x (radial, continuous, integration variable)',
            'solution': 'Associated Laguerre basis L_n^{|m|}(x)',
            'mechanism': 'Absorb singularity into weight: x^|m| e^{-x}',
            'result': 'Algebraic matrix elements, faster convergence',
        },
        'level4_cusp': {
            'type': 'Coulomb cusp',
            'location': 'theta_12=0, alpha=pi/4 (coalescence)',
            'dimension': '2D',
            'form': '1/r_12 = 1/(R_e * sqrt(1 - sin(2a)cos(theta_12)))',
            'coordinate': 'theta_12 (encoded in channel structure, NOT integration variable)',
            'obstruction_1': (
                'theta_12 is not a separate coordinate with a weight function. '
                'It is encoded in the coupled spherical harmonic expansion '
                'through the channel indices (l1, l2) and the Gaunt integral coupling.'
            ),
            'obstruction_2': (
                'The cusp is 2D in (alpha, theta_12). Even if theta_12 were '
                'a separate coordinate, a 1D weight modification would not '
                'suffice (Track U proved alpha-only fails).'
            ),
            'obstruction_3': (
                'Channel quantum numbers (l1, l2) are discrete and constrained '
                'by Gaunt integral selection rules. Unlike the Laguerre weight '
                'exponent |m| which can be absorbed continuously, there is no '
                'continuous parameter to absorb the cusp into.'
            ),
        },
        'analogy_assessment': (
            'The Track J pattern does NOT transfer to the theta_12 cusp. '
            'At Level 2, the singularity lives in a 1D coordinate that is '
            'an explicit integration variable with a modifiable weight function. '
            'At Level 4, the theta_12 dependence is encoded in the DISCRETE '
            'channel structure through Gaunt integrals, and the cusp is 2D. '
            'There is no weight function to modify, and no continuous parameter '
            'to absorb the singularity into. The correct approach is the '
            'Schwartz extrapolation (Track X), which treats the slow channel '
            'convergence as a known asymptotic series and extrapolates, rather '
            'than trying to modify the basis.'
        ),
        'what_works_instead': (
            'Track X (perturbative cusp correction) via Schwartz partial-wave '
            'extrapolation. This treats the missing high-l channel contributions '
            'as a known series Delta_E_l ~ A/(l+1/2)^4 and sums the tail '
            'analytically. This is not a basis modification but a post-processing '
            'correction that leverages the KNOWN cusp physics. It breaks through '
            'the 0.19-0.20% convergence floor at Level 3 and provides ~1 pp D_e '
            'improvement at Level 4.'
        ),
        'classification': (
            'NEGATIVE RESULT: The Track J weight-absorption pattern requires '
            'a 1D singularity in a continuous coordinate with a modifiable weight. '
            'The Level 4 cusp violates all three conditions (2D, discrete channels, '
            'no weight function). The cusp is classified as an embedding exchange '
            'constant (Track W) that must be handled by extrapolation (Track X), '
            'not basis adaptation.'
        ),
    }
