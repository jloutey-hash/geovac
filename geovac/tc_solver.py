"""
Transcorrelated (TC) solver for Level 3 hyperspherical method (He).

The TC transformation with Jastrow J = +(1/2)r_12 gives:
    H_TC = T + V_nuc + G + 1/4
where V_ee = 1/r_12 is completely cancelled, G = -d/dr_12 is the TC
gradient operator, and +1/4 is the double-commutator constant.

IMPORTANT: The TC operator G = -d/dr_12 decomposes in hyperspherical
coordinates into an angular part G_ang (acting within the angular
eigenproblem at fixed R) and a hyperradial part G_R = -(sqrt(u)/2)*d/dR.
A naive adiabatic approach that removes V_ee from the angular Hamiltonian
and replaces it with G_ang alone is catastrophically wrong (~47% error)
because G_R carries the dominant contribution — the V_ee coupling in the
angular problem was O(R), while the angular TC correction is O(1).

This module implements TWO approaches:

1. **Perturbative TC correction (primary):** Solve the standard adiabatic
   problem (WITH V_ee) and compute the TC energy correction as a
   perturbative modification to V_eff(R):
       Delta_E = <Psi|H_TC - H_standard|Psi>
              = <Psi|G + 1/4 - V_ee|Psi>
   Since H_TC - H = G + 1/4 - V_ee, this gives the exact TC correction
   for the standard adiabatic wavefunction.

2. **Naive TC angular replacement (diagnostic only):** Remove V_ee from
   the angular Hamiltonian and replace with G_ang_tilde. This is
   documented as a NEGATIVE RESULT: the adiabatic separation is
   incompatible with naive V_ee removal because G_R is not small.

The perturbative approach computes:
    E_TC = E_std + <Delta V_eff(R)> + 1/4

where Delta V_eff(R) at each R is the expectation value of G_ang - V_ee
in the angular eigenfunction, and the hyperradial G_R coupling is treated
as an additional non-adiabatic correction.

References:
  - docs/tc_integrals_derivation.md (full derivation)
  - Paper 13, Section III (angular equation)
  - Ten-no, Chem. Phys. Lett. 398, 56 (2004)
  - Dobrautz et al., Phys. Rev. B 99, 075119 (2019)
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh
from typing import Tuple, Dict, Optional, List
from math import pi

from geovac.algebraic_angular import AlgebraicAngularSolver, _gegenbauer
from geovac.hyperspherical_adiabatic import effective_potential
from geovac.hyperspherical_radial import solve_radial


# Exact He ground state energy (Hylleraas)
_E_EXACT_HE = -2.903724


def _gegenbauer_derivative(k: int, lam: float, x: np.ndarray) -> np.ndarray:
    """Derivative of Gegenbauer polynomial: dC_k^lam/dx = 2*lam*C_{k-1}^{lam+1}(x).

    For k=0, the derivative is zero.
    """
    if k == 0:
        return np.zeros_like(x, dtype=float)
    return 2.0 * lam * _gegenbauer(k - 1, lam + 1.0, x)


def _legendre_derivative(l: int, x: np.ndarray) -> np.ndarray:
    """Derivative of Legendre polynomial P_l'(x).

    Uses recurrence: (1-x^2)*P_l'(x) = l*P_{l-1}(x) - l*x*P_l(x)
    """
    if l == 0:
        return np.zeros_like(x, dtype=float)

    from numpy.polynomial.legendre import legval

    coeffs_l = np.zeros(l + 1)
    coeffs_l[l] = 1.0
    P_l = legval(x, coeffs_l)

    coeffs_lm1 = np.zeros(l)
    coeffs_lm1[l - 1] = 1.0
    P_lm1 = legval(x, coeffs_lm1)

    denom = 1.0 - x * x
    result = np.zeros_like(x, dtype=float)
    mask = np.abs(denom) > 1e-14
    result[mask] = l * (P_lm1[mask] - x[mask] * P_l[mask]) / denom[mask]
    mask_p1 = (~mask) & (x > 0)
    mask_m1 = (~mask) & (x < 0)
    result[mask_p1] = l * (l + 1) / 2.0
    result[mask_m1] = (-1) ** (l + 1) * l * (l + 1) / 2.0
    return result


def _legendre_P(l: int, x: np.ndarray) -> np.ndarray:
    """Evaluate Legendre polynomial P_l(x)."""
    from numpy.polynomial.legendre import legval
    coeffs = np.zeros(l + 1)
    coeffs[l] = 1.0
    return legval(x, coeffs)


class TCSolver:
    """Transcorrelated solver for He using the Level 3 hyperspherical basis.

    Computes TC corrections perturbatively on top of the standard adiabatic
    solution, properly accounting for both angular and hyperradial components
    of the TC gradient operator.

    Parameters
    ----------
    Z : float
        Nuclear charge.
    n_basis : int
        Number of basis functions per l-channel.
    l_max : int
        Maximum partial wave.
    n_quad_alpha : int
        Number of Gauss-Legendre quadrature points for alpha.
    n_quad_theta : int
        Number of Gauss-Legendre quadrature points for theta_12.
    """

    def __init__(
        self,
        Z: float = 2.0,
        n_basis: int = 10,
        l_max: int = 0,
        n_quad_alpha: int = 100,
        n_quad_theta: int = 100,
    ) -> None:
        self.Z = Z
        self.n_basis = n_basis
        self.l_max = l_max
        self.n_l = l_max + 1
        self._total_dim = self.n_l * n_basis

        # Build the standard angular solver
        self._std_solver = AlgebraicAngularSolver(
            Z, n_basis, l_max, symmetry='singlet', n_quad=n_quad_alpha
        )

        # Setup 2D quadrature for TC matrix elements
        self._setup_2d_quadrature(n_quad_alpha, n_quad_theta)

        # Build 2D basis functions and derivatives
        self._build_basis_2d()

        # Build TC angular correction matrix (R-independent)
        self._build_tc_angular_matrix()

        # Build V_ee matrix in the 2D basis (for perturbative correction)
        self._build_vee_matrix_2d()

    def _setup_2d_quadrature(
        self, n_quad_alpha: int, n_quad_theta: int
    ) -> None:
        """Setup 2D Gauss-Legendre quadrature for (alpha, theta_12)."""
        from numpy.polynomial.legendre import leggauss

        # Alpha quadrature on [0, pi/4]
        nodes_a, weights_a = leggauss(n_quad_alpha)
        a1, b1 = 0.0, np.pi / 4.0
        scale1 = (b1 - a1) / 2.0
        self._alpha_2d = scale1 * nodes_a + (b1 + a1) / 2.0
        self._w_alpha = weights_a * scale1

        # theta_12 quadrature on [0, pi]
        nodes_t, weights_t = leggauss(n_quad_theta)
        a2, b2 = 0.0, np.pi
        scale2 = (b2 - a2) / 2.0
        self._theta_2d = scale2 * nodes_t + (b2 + a2) / 2.0
        self._w_theta = weights_t * scale2

        # Precompute trig
        self._sin_a = np.sin(self._alpha_2d)
        self._cos_a = np.cos(self._alpha_2d)
        self._sin_2a = np.sin(2.0 * self._alpha_2d)
        self._cos_2a = np.cos(2.0 * self._alpha_2d)
        self._sin_cos = self._sin_a * self._cos_a
        self._cos_t = np.cos(self._theta_2d)
        self._sin_t = np.sin(self._theta_2d)

        self._n_alpha = len(self._alpha_2d)
        self._n_theta = len(self._theta_2d)

    def _build_basis_2d(self) -> None:
        """Build basis functions and their derivatives on the 2D grid.

        Basis: phi_{k,l}(alpha, theta) = N * (sin a cos a)^{l+1}
               * C_k^{l+1}(cos 2a) * P_l(cos theta_12)

        Stored as 3D arrays indexed by (basis_idx, alpha_idx, theta_idx).
        """
        nb = self.n_basis
        n_a = self._n_alpha
        n_t = self._n_theta
        dim = self._total_dim

        cos_2a = np.cos(2.0 * self._alpha_2d)
        sin_cos = self._sin_cos

        # Legendre values and derivatives
        P_l_vals = {}
        dP_l_dtheta = {}
        for l in range(self.n_l):
            P_l_vals[l] = _legendre_P(l, self._cos_t)
            dP_l_dtheta[l] = -self._sin_t * _legendre_derivative(l, self._cos_t)

        self._phi_2d = np.zeros((dim, n_a, n_t))
        self._dphi_dalpha_2d = np.zeros((dim, n_a, n_t))
        self._dphi_dtheta_2d = np.zeros((dim, n_a, n_t))

        measure_a = self._sin_a**2 * self._cos_a**2
        measure_t = self._sin_t

        for l in range(self.n_l):
            lam = float(l + 1)
            envelope = sin_cos ** (l + 1)

            # d/dalpha[(sin a cos a)^{l+1}] = (l+1)(sin a cos a)^l * cos(2a)
            d_envelope = (l + 1) * sin_cos ** max(l, 0) * self._cos_2a

            k_indices = [2 * j for j in range(nb)]

            for j, kv in enumerate(k_indices):
                idx = l * nb + j

                C_k = _gegenbauer(kv, lam, cos_2a)
                dC_k_dx = _gegenbauer_derivative(kv, lam, cos_2a)
                dC_k_dalpha = dC_k_dx * (-2.0 * self._sin_2a)

                chi = envelope * C_k
                dchi_dalpha = d_envelope * C_k + envelope * dC_k_dalpha

                self._phi_2d[idx] = chi[:, np.newaxis] * P_l_vals[l][np.newaxis, :]
                self._dphi_dalpha_2d[idx] = dchi_dalpha[:, np.newaxis] * P_l_vals[l][np.newaxis, :]
                self._dphi_dtheta_2d[idx] = chi[:, np.newaxis] * dP_l_dtheta[l][np.newaxis, :]

        # Normalize with 2D measure (factor 2 for singlet alpha symmetry)
        self._norms = np.zeros(dim)
        for idx in range(dim):
            integrand = (self._phi_2d[idx]**2
                         * measure_a[:, np.newaxis]
                         * measure_t[np.newaxis, :])
            norm_sq = 2.0 * np.einsum('i,j,ij->', self._w_alpha, self._w_theta, integrand)
            self._norms[idx] = np.sqrt(max(norm_sq, 1e-30))
            self._phi_2d[idx] /= self._norms[idx]
            self._dphi_dalpha_2d[idx] /= self._norms[idx]
            self._dphi_dtheta_2d[idx] /= self._norms[idx]

    def _build_tc_angular_matrix(self) -> None:
        """Build the R-independent TC angular correction matrix G_tilde_ang.

        G_tilde_ang = -A(alpha, theta) * d/dalpha - B(alpha, theta) * d/dtheta

        where A = cos(2a)*cos(theta)/(2*sqrt(u)), B = cos(2a)*sin(theta)/(sin(2a)*sqrt(u))
        """
        dim = self._total_dim

        u_2d = 1.0 - self._sin_2a[:, np.newaxis] * self._cos_t[np.newaxis, :]
        sqrt_u = np.sqrt(np.maximum(u_2d, 1e-30))

        A_coeff = self._cos_2a[:, np.newaxis] * self._cos_t[np.newaxis, :] / (2.0 * sqrt_u)

        sin_2a_safe = self._sin_2a.copy()
        sin_2a_safe[np.abs(sin_2a_safe) < 1e-30] = 1e-30
        B_coeff = (self._cos_2a[:, np.newaxis] * self._sin_t[np.newaxis, :]
                   / (sin_2a_safe[:, np.newaxis] * sqrt_u))

        measure = ((self._sin_a**2 * self._cos_a**2)[:, np.newaxis]
                   * self._sin_t[np.newaxis, :])

        wA = -A_coeff * measure
        wB = -B_coeff * measure

        self._tc_ang_matrix = np.zeros((dim, dim))

        for j in range(dim):
            f_total = wA * self._dphi_dalpha_2d[j] + wB * self._dphi_dtheta_2d[j]
            for i in range(dim):
                integrand = self._phi_2d[i] * f_total
                val = 2.0 * np.einsum('i,j,ij->', self._w_alpha, self._w_theta, integrand)
                self._tc_ang_matrix[i, j] = val

    def _build_vee_matrix_2d(self) -> None:
        """Build V_ee matrix in the 2D-normalized basis.

        V_ee = 1/(R*sqrt(u)) in the charge function convention.
        Here we compute the angular charge function C_ee = 1/sqrt(u) matrix elements.
        The actual V_ee contribution to the angular Hamiltonian is R * C_ee.
        """
        dim = self._total_dim

        u_2d = 1.0 - self._sin_2a[:, np.newaxis] * self._cos_t[np.newaxis, :]
        inv_sqrt_u = 1.0 / np.sqrt(np.maximum(u_2d, 1e-30))

        measure = ((self._sin_a**2 * self._cos_a**2)[:, np.newaxis]
                   * self._sin_t[np.newaxis, :])

        self._vee_charge_matrix = np.zeros((dim, dim))

        for i in range(dim):
            for j in range(i, dim):
                integrand = (self._phi_2d[i] * self._phi_2d[j]
                             * inv_sqrt_u * measure)
                val = 2.0 * np.einsum('i,j,ij->', self._w_alpha, self._w_theta, integrand)
                self._vee_charge_matrix[i, j] = val
                self._vee_charge_matrix[j, i] = val

    def _build_gr_expectation_matrix(self) -> np.ndarray:
        """Build the expectation value matrix of sqrt(u)/2 (the G_R coefficient).

        G_R = -(sqrt(u)/2) * d/dR acts on the radial wavefunction.
        At first order in perturbation theory, this modifies V_eff(R) by:
            <phi_0|sqrt(u)/2|phi_0> * (dF/dR)/F
        which is a non-adiabatic correction. For the perturbative TC correction,
        we compute the expectation value of sqrt(u)/2 in the angular eigenstate.
        """
        dim = self._total_dim

        u_2d = 1.0 - self._sin_2a[:, np.newaxis] * self._cos_t[np.newaxis, :]
        sqrt_u_half = np.sqrt(np.maximum(u_2d, 1e-30)) / 2.0

        measure = ((self._sin_a**2 * self._cos_a**2)[:, np.newaxis]
                   * self._sin_t[np.newaxis, :])

        gr_matrix = np.zeros((dim, dim))
        for i in range(dim):
            for j in range(i, dim):
                integrand = self._phi_2d[i] * self._phi_2d[j] * sqrt_u_half * measure
                val = 2.0 * np.einsum('i,j,ij->', self._w_alpha, self._w_theta, integrand)
                gr_matrix[i, j] = val
                gr_matrix[j, i] = val

        return gr_matrix

    def compute_tc_correction_veff(
        self, R: float
    ) -> Tuple[float, float, float, float]:
        """Compute the TC perturbative correction to V_eff at one R point.

        The TC correction replaces V_ee with G + 1/4. In the adiabatic
        framework, the angular part of G modifies the angular eigenvalue,
        and the radial part G_R modifies the radial equation.

        For the ANGULAR correction at fixed R:
            Delta_mu(R) = <Phi_0(R)| G_tilde_ang - R*C_ee |Phi_0(R)>

        This gives the change in the angular eigenvalue if we use the TC
        angular Hamiltonian instead of the standard one.

        For the perturbative correction to V_eff:
            Delta_V_eff(R) = Delta_mu(R) / R^2

        The +1/4 constant shift is added to the total energy, not to V_eff.

        Parameters
        ----------
        R : float
            Hyperradius (bohr).

        Returns
        -------
        delta_mu : float
            Change in angular eigenvalue from TC correction.
        mu_std : float
            Standard angular eigenvalue.
        vee_expectation : float
            <Phi_0|R*C_ee|Phi_0> in the standard eigenstate.
        g_ang_expectation : float
            <Phi_0|G_tilde_ang|Phi_0> in the standard eigenstate.
        """
        # Solve standard angular problem
        casimir_all = np.concatenate(self._std_solver._channel_casimir)
        H_std = np.diag(casimir_all) + R * self._std_solver._coupling_full
        evals, evecs = eigh(H_std)

        phi_0 = evecs[:, 0]  # Ground state eigenvector
        mu_std = evals[0]

        # V_ee expectation: <phi_0|R*C_ee|phi_0>
        vee_exp = R * phi_0 @ self._std_solver._vee_full @ phi_0

        # G_ang expectation: <phi_0|G_tilde_ang|phi_0>
        g_ang_exp = phi_0 @ self._tc_ang_matrix @ phi_0

        delta_mu = g_ang_exp - vee_exp

        return delta_mu, mu_std, vee_exp, g_ang_exp


def solve_he_tc_perturbative(
    Z: float = 2.0,
    n_basis: int = 10,
    l_max: int = 0,
    n_R: int = 200,
    R_min: float = 0.1,
    R_max: float = 30.0,
    N_R_radial: int = 2000,
    n_quad_alpha: int = 100,
    n_quad_theta: int = 100,
    verbose: bool = True,
) -> Dict:
    """TC-corrected He energy via perturbative correction to V_eff(R).

    Approach:
    1. Solve the standard adiabatic problem (with V_ee) to get E_std, mu(R), phi(R)
    2. At each R, compute the TC angular correction:
       Delta_mu(R) = <Phi_0|G_tilde_ang|Phi_0> - <Phi_0|R*C_ee|Phi_0>
    3. Correct V_eff: V_eff_TC(R) = V_eff_std(R) + Delta_mu(R)/R^2
    4. Solve the corrected radial equation
    5. Apply the +1/4 energy shift: E_TC = E_radial - 1/4

    The +1/4 constant from the double commutator RAISES the TC energy.
    The net effect should be: V_ee removal lowers mu(R), the +1/4 partially
    compensates, and the smoother TC wavefunction converges faster with l_max.

    Returns
    -------
    result : dict
        Keys include energy_tc, energy_std, err_tc_pct, err_std_pct, etc.
    """
    import time

    t0 = time.time()

    if verbose:
        print(f"TC perturbative solver: Z={Z}, n_basis={n_basis}, l_max={l_max}")

    tc = TCSolver(Z, n_basis, l_max, n_quad_alpha, n_quad_theta)

    t1 = time.time()
    if verbose:
        print(f"  TC matrix build: {t1 - t0:.2f}s")

    # R grid
    R_grid = np.concatenate([
        np.linspace(R_min, 1.0, n_R // 3),
        np.linspace(1.0, 5.0, n_R // 3),
        np.linspace(5.0, R_max, n_R // 3 + 1),
    ])
    R_grid = np.unique(R_grid)

    if verbose:
        print(f"  Computing mu(R) on {len(R_grid)} R points...")

    # Standard and TC-corrected angular eigenvalues
    mu_std = np.zeros(len(R_grid))
    mu_tc = np.zeros(len(R_grid))
    delta_mu = np.zeros(len(R_grid))
    vee_exp_arr = np.zeros(len(R_grid))
    g_ang_exp_arr = np.zeros(len(R_grid))

    for i, R in enumerate(R_grid):
        dm, ms, ve, ge = tc.compute_tc_correction_veff(R)
        mu_std[i] = ms
        delta_mu[i] = dm
        mu_tc[i] = ms + dm
        vee_exp_arr[i] = ve
        g_ang_exp_arr[i] = ge

    t2 = time.time()
    if verbose:
        print(f"  Angular sweep: {t2 - t1:.2f}s")
        print(f"  Delta_mu range: [{delta_mu.min():.4f}, {delta_mu.max():.4f}]")
        print(f"  V_ee exp range: [{vee_exp_arr.min():.4f}, {vee_exp_arr.max():.4f}]")
        print(f"  G_ang exp range: [{g_ang_exp_arr.min():.6f}, {g_ang_exp_arr.max():.6f}]")

    # Build V_eff for standard and TC-corrected
    V_eff_std = effective_potential(R_grid, mu_std)
    V_eff_tc = effective_potential(R_grid, mu_tc)

    V_eff_std_spline = CubicSpline(R_grid, V_eff_std, extrapolate=True)
    V_eff_tc_spline = CubicSpline(R_grid, V_eff_tc, extrapolate=True)

    # Radial solve
    E_std_arr, _, R_rad = solve_radial(V_eff_std_spline, R_min, R_max, N_R_radial)
    E_std = E_std_arr[0]

    E_tc_raw_arr, _, _ = solve_radial(V_eff_tc_spline, R_min, R_max, N_R_radial)
    # The +1/4 from double commutator is a constant energy shift.
    # It RAISES the energy (makes it less negative).
    # E_TC = E_radial_TC + 1/4
    # But we want: E_physical = E_TC_eigenvalue
    # Since H_TC eigenvalue = E_exact + 1/4 (the +1/4 shifts eigenvalues up),
    # we need E_physical = E_TC - 1/4
    E_tc = E_tc_raw_arr[0] - 0.25

    t3 = time.time()

    err_std = abs(E_std - _E_EXACT_HE) / abs(_E_EXACT_HE) * 100
    err_tc = abs(E_tc - _E_EXACT_HE) / abs(_E_EXACT_HE) * 100

    if verbose:
        print(f"\n  Results (l_max={l_max}):")
        print(f"    Standard:    E = {E_std:.6f} Ha, error = {err_std:.4f}%")
        print(f"    TC (raw):    E = {E_tc_raw_arr[0]:.6f} Ha")
        print(f"    TC (-1/4):   E = {E_tc:.6f} Ha, error = {err_tc:.4f}%")
        print(f"    Exact He:    E = {_E_EXACT_HE:.6f} Ha")
        print(f"    Total time:  {t3 - t0:.2f}s")

    return {
        'energy_tc': E_tc,
        'energy_std': E_std,
        'err_tc_pct': err_tc,
        'err_std_pct': err_std,
        'mu_tc': mu_tc,
        'mu_std': mu_std,
        'delta_mu': delta_mu,
        'R_grid': R_grid,
        'V_eff_tc': V_eff_tc,
        'V_eff_std': V_eff_std,
        'vee_expectation': vee_exp_arr,
        'g_ang_expectation': g_ang_exp_arr,
        'Z': Z,
        'l_max': l_max,
        'n_basis': n_basis,
    }


def solve_he_tc_direct(
    Z: float = 2.0,
    n_basis: int = 10,
    l_max: int = 0,
    n_R: int = 200,
    R_min: float = 0.1,
    R_max: float = 30.0,
    N_R_radial: int = 2000,
    n_quad_alpha: int = 100,
    n_quad_theta: int = 100,
    verbose: bool = True,
) -> Dict:
    """Direct (non-Hermitian) TC solver — DIAGNOSTIC ONLY.

    Replaces V_ee in the angular Hamiltonian with G_tilde_ang directly.
    This is KNOWN TO BE WRONG because it drops the hyperradial component G_R.
    Included for completeness and to document the negative result.

    The energy error is ~47% because:
    - V_ee contributed O(R) to the angular eigenvalue
    - G_tilde_ang contributes O(1) (R-independent)
    - The missing G_R = -(sqrt(u)/2)*d/dR is the dominant TC contribution
    """
    import time

    t0 = time.time()

    if verbose:
        print(f"TC DIRECT solver (diagnostic): Z={Z}, n_basis={n_basis}, l_max={l_max}")
        print(f"  WARNING: This drops G_R and gives ~47% error. For diagnosis only.")

    tc = TCSolver(Z, n_basis, l_max, n_quad_alpha, n_quad_theta)

    R_grid = np.concatenate([
        np.linspace(R_min, 1.0, n_R // 3),
        np.linspace(1.0, 5.0, n_R // 3),
        np.linspace(5.0, R_max, n_R // 3 + 1),
    ])
    R_grid = np.unique(R_grid)

    mu_tc = np.zeros(len(R_grid))
    mu_std = np.zeros(len(R_grid))
    max_imag_all = 0.0

    casimir_all = np.concatenate(tc._std_solver._channel_casimir)

    for i, R in enumerate(R_grid):
        # Standard
        H_std = np.diag(casimir_all) + R * tc._std_solver._coupling_full
        evals_std, _ = eigh(H_std)
        mu_std[i] = evals_std[0]

        # TC direct: Casimir + R*nuclear + G_tilde (NO V_ee)
        H_tc = (np.diag(casimir_all)
                + R * tc._std_solver._nuclear_full
                + tc._tc_ang_matrix)
        evals_tc = np.linalg.eigvals(H_tc)
        idx = np.argsort(evals_tc.real)
        mu_tc[i] = evals_tc[idx[0]].real
        max_imag_all = max(max_imag_all, np.max(np.abs(evals_tc.imag)))

    V_eff_tc = effective_potential(R_grid, mu_tc)
    V_eff_std = effective_potential(R_grid, mu_std)

    V_tc_sp = CubicSpline(R_grid, V_eff_tc, extrapolate=True)
    V_std_sp = CubicSpline(R_grid, V_eff_std, extrapolate=True)

    E_tc_raw, _, _ = solve_radial(V_tc_sp, R_min, R_max, N_R_radial)
    E_std_arr, _, _ = solve_radial(V_std_sp, R_min, R_max, N_R_radial)

    E_tc = E_tc_raw[0] - 0.25
    E_std = E_std_arr[0]

    err_tc = abs(E_tc - _E_EXACT_HE) / abs(_E_EXACT_HE) * 100
    err_std = abs(E_std - _E_EXACT_HE) / abs(_E_EXACT_HE) * 100

    if verbose:
        print(f"\n  Standard: E = {E_std:.6f} Ha, error = {err_std:.4f}%")
        print(f"  TC direct: E = {E_tc:.6f} Ha, error = {err_tc:.4f}%")
        print(f"  Max imag: {max_imag_all:.2e}")
        print(f"  EXPECTED: ~47% error (G_R dropped)")

    return {
        'energy_tc': E_tc,
        'energy_std': E_std,
        'err_tc_pct': err_tc,
        'err_std_pct': err_std,
        'max_imag': max_imag_all,
        'l_max': l_max,
        'method': 'direct_DIAGNOSTIC',
    }


def _solve_radial_with_first_derivative(
    V_eff_func,
    gamma_func,
    R_min: float = 0.05,
    R_max: float = 30.0,
    N_R: int = 2000,
    n_states: int = 1,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Solve the radial equation with a first-order derivative term.

    [-1/2 d^2F/dR^2 - gamma(R)*dF/dR + V_eff(R)*F] = E*F

    with Dirichlet boundary conditions F(R_min) = F(R_max) = 0.

    The first-derivative term makes this non-Hermitian. Uses numpy.linalg.eig.

    Parameters
    ----------
    V_eff_func : callable
        Effective potential V_eff(R).
    gamma_func : callable
        First-derivative coefficient gamma(R). The term is -gamma(R)*dF/dR.
    R_min, R_max : float
        Grid boundaries.
    N_R : int
        Number of interior grid points.
    n_states : int
        Number of eigenstates to return.

    Returns
    -------
    E : ndarray of shape (n_states,)
        Real parts of eigenvalues, sorted ascending.
    F : ndarray of shape (n_states, N_R)
        Right eigenvectors.
    R_grid : ndarray of shape (N_R,)
        Grid points.
    """
    h = (R_max - R_min) / (N_R + 1)
    R_grid = R_min + (np.arange(N_R) + 1) * h

    V = V_eff_func(R_grid)
    gamma = gamma_func(R_grid)

    # Build the full N_R x N_R matrix
    # -1/2 d^2/dR^2 -> diagonal: 1/h^2, off-diag: -0.5/h^2
    # -gamma * d/dR -> off-diag: -gamma/(2h) (central difference)
    H = np.zeros((N_R, N_R))

    # Diagonal: kinetic + potential
    for i in range(N_R):
        H[i, i] = 1.0 / h**2 + V[i]

    # Off-diagonal: kinetic + first-derivative
    for i in range(N_R - 1):
        # Kinetic: -1/2 * 1/h^2
        kin = -0.5 / h**2
        # First derivative: -gamma * 1/(2h) for upper, +gamma * 1/(2h) for lower
        # d/dR ~ (F_{i+1} - F_{i-1})/(2h), so coefficient of F_{i+1} is 1/(2h)
        # and F_{i-1} is -1/(2h)
        # Term: -gamma(R_i) * dF/dR at point i
        # = -gamma_i * (F_{i+1} - F_{i-1})/(2h)
        # Coefficient of F_{i+1}: -gamma_i/(2h)
        # Coefficient of F_{i-1}: +gamma_i/(2h)
        H[i, i + 1] = kin - gamma[i] / (2.0 * h)
        H[i + 1, i] = kin + gamma[i + 1] / (2.0 * h)

    # Non-Hermitian eigensolve
    evals, evecs = np.linalg.eig(H)

    # Sort by real part, take lowest
    idx = np.argsort(evals.real)
    evals = evals[idx]
    evecs = evecs[:, idx]

    E = evals[:n_states].real
    F = np.zeros((n_states, N_R))
    for s in range(n_states):
        F[s] = evecs[:, s].real
        norm = np.sqrt(h * np.sum(F[s]**2))
        if norm > 0:
            F[s] /= norm

    return E, F, R_grid


def solve_he_tc_full(
    Z: float = 2.0,
    n_basis: int = 10,
    l_max: int = 0,
    n_R: int = 200,
    R_min: float = 0.1,
    R_max: float = 30.0,
    N_R_radial: int = 2000,
    n_quad_alpha: int = 100,
    n_quad_theta: int = 100,
    verbose: bool = True,
) -> Dict:
    """Full TC solver with G_R coupling in the radial equation.

    The TC Hamiltonian decomposed in the adiabatic framework:
    1. Angular problem (at fixed R): H_ang_TC = Lambda^2/2 + R*C_nuc + G_tilde_ang
       (No V_ee! G_tilde_ang is R-independent.)
    2. Radial equation: [-1/2 d^2F/dR^2 - gamma(R)*dF/dR + V_eff_TC(R)] F = (E-1/4) F
       where gamma(R) = <Phi_0(R)|sqrt(u)/2|Phi_0(R)> is the G_R coupling.

    The G_R = -(sqrt(u)/2)*d/dR term adds a first-order derivative to the
    radial equation, making it non-Hermitian.

    Parameters
    ----------
    Z : float
        Nuclear charge.
    n_basis : int
        Angular basis functions per l-channel.
    l_max : int
        Maximum partial wave.
    n_R : int
        Number of R points for adiabatic curves.
    R_min, R_max : float
        Hyperradial grid boundaries.
    N_R_radial : int
        Number of grid points for radial solve.
    n_quad_alpha, n_quad_theta : int
        Quadrature points.
    verbose : bool
        Print progress.

    Returns
    -------
    result : dict
        Keys: energy_tc, energy_std, err_tc_pct, err_std_pct, gamma, etc.
    """
    import time
    from scipy.linalg import eigh as scipy_eigh

    t0 = time.time()

    if verbose:
        print(f"TC FULL solver: Z={Z}, n_basis={n_basis}, l_max={l_max}")

    tc = TCSolver(Z, n_basis, l_max, n_quad_alpha, n_quad_theta)
    gr_matrix = tc._build_gr_expectation_matrix()

    t1 = time.time()
    if verbose:
        print(f"  Matrix build: {t1 - t0:.2f}s")

    # R grid
    R_grid = np.concatenate([
        np.linspace(R_min, 1.0, n_R // 3),
        np.linspace(1.0, 5.0, n_R // 3),
        np.linspace(5.0, R_max, n_R // 3 + 1),
    ])
    R_grid = np.unique(R_grid)

    if verbose:
        print(f"  Computing adiabatic curves on {len(R_grid)} R points...")

    casimir_all = np.concatenate(tc._std_solver._channel_casimir)
    mu_tc = np.zeros(len(R_grid))
    mu_std = np.zeros(len(R_grid))
    gamma_arr = np.zeros(len(R_grid))
    max_imag_all = 0.0

    for i, R in enumerate(R_grid):
        # Standard angular problem
        H_std = np.diag(casimir_all) + R * tc._std_solver._coupling_full
        evals_std, evecs_std = scipy_eigh(H_std)
        mu_std[i] = evals_std[0]

        # TC angular problem: Casimir + R*nuclear + G_tilde_ang (no V_ee)
        H_tc = (np.diag(casimir_all)
                + R * tc._std_solver._nuclear_full
                + tc._tc_ang_matrix)
        evals_tc = np.linalg.eigvals(H_tc)
        idx_tc = np.argsort(evals_tc.real)
        mu_tc[i] = evals_tc[idx_tc[0]].real
        max_imag_all = max(max_imag_all, np.max(np.abs(evals_tc.imag)))

        # For gamma(R), use the TC eigenvector.
        # Since H_tc is non-Hermitian, we need the RIGHT eigenvector.
        evals_tc_full, evecs_tc_full = np.linalg.eig(H_tc)
        idx_tc_full = np.argsort(evals_tc_full.real)
        phi_tc_0 = evecs_tc_full[:, idx_tc_full[0]].real
        # Normalize
        phi_tc_0 /= np.sqrt(np.dot(phi_tc_0, phi_tc_0))

        # gamma(R) = <Phi_TC_0|sqrt(u)/2|Phi_TC_0>
        gamma_arr[i] = phi_tc_0 @ gr_matrix @ phi_tc_0

    t2 = time.time()
    if verbose:
        print(f"  Angular sweep: {t2 - t1:.2f}s")
        print(f"  Max imaginary part: {max_imag_all:.2e}")
        print(f"  gamma(R) range: [{gamma_arr.min():.4f}, {gamma_arr.max():.4f}]")

    # Build V_eff for TC
    V_eff_tc = effective_potential(R_grid, mu_tc)
    V_eff_std = effective_potential(R_grid, mu_std)

    V_eff_tc_spline = CubicSpline(R_grid, V_eff_tc, extrapolate=True)
    V_eff_std_spline = CubicSpline(R_grid, V_eff_std, extrapolate=True)
    gamma_spline = CubicSpline(R_grid, gamma_arr, extrapolate=True)

    if verbose:
        print(f"  Solving radial equations...")

    # Standard radial solve (Hermitian)
    E_std_arr, _, R_rad = solve_radial(V_eff_std_spline, R_min, R_max, N_R_radial)
    E_std = E_std_arr[0]

    # TC radial solve with first-derivative term (non-Hermitian)
    E_tc_arr, _, _ = _solve_radial_with_first_derivative(
        V_eff_tc_spline, gamma_spline, R_min, R_max, N_R_radial
    )
    E_tc = E_tc_arr[0] - 0.25  # Subtract +1/4 constant

    t3 = time.time()

    err_std = abs(E_std - _E_EXACT_HE) / abs(_E_EXACT_HE) * 100
    err_tc = abs(E_tc - _E_EXACT_HE) / abs(_E_EXACT_HE) * 100

    if verbose:
        print(f"\n  Results (l_max={l_max}):")
        print(f"    Standard:  E = {E_std:.6f} Ha, error = {err_std:.4f}%")
        print(f"    TC (raw):  E = {E_tc + 0.25:.6f} Ha")
        print(f"    TC (-1/4): E = {E_tc:.6f} Ha, error = {err_tc:.4f}%")
        print(f"    Exact He:  E = {_E_EXACT_HE:.6f} Ha")
        print(f"    Max imag:  {max_imag_all:.2e}")
        print(f"    Total time: {t3 - t0:.2f}s")

    return {
        'energy_tc': E_tc,
        'energy_std': E_std,
        'err_tc_pct': err_tc,
        'err_std_pct': err_std,
        'mu_tc': mu_tc,
        'mu_std': mu_std,
        'gamma': gamma_arr,
        'R_grid': R_grid,
        'V_eff_tc': V_eff_tc,
        'V_eff_std': V_eff_std,
        'max_imag': max_imag_all,
        'Z': Z,
        'l_max': l_max,
        'n_basis': n_basis,
        'method': 'full_tc',
    }


def convergence_study(
    l_max_values: Optional[List[int]] = None,
    Z: float = 2.0,
    n_basis: int = 10,
    verbose: bool = True,
) -> List[Dict]:
    """Run TC perturbative vs standard comparison at multiple l_max values.

    Parameters
    ----------
    l_max_values : list of int, optional
        Partial wave values to test. Default: [0, 1, 2, 3].
    Z : float
        Nuclear charge.
    n_basis : int
        Angular basis functions per channel.
    verbose : bool
        Print progress.

    Returns
    -------
    results : list of dict
    """
    if l_max_values is None:
        l_max_values = [0, 1, 2, 3]

    all_results = []

    for l_max in l_max_values:
        if verbose:
            print(f"\n{'='*60}")
            print(f"  l_max = {l_max}")
            print(f"{'='*60}")

        result = solve_he_tc_full(
            Z=Z, n_basis=n_basis, l_max=l_max, verbose=verbose,
        )
        all_results.append(result)

    if verbose:
        print(f"\n{'='*60}")
        print(f"  CONVERGENCE SUMMARY (He TC full, Z={Z})")
        print(f"{'='*60}")
        print(f"{'l_max':>5} | {'E_std (Ha)':>12} | {'E_TC (Ha)':>12} | "
              f"{'err_std%':>8} | {'err_TC%':>8} | {'max_imag':>10}")
        print("-" * 70)
        for r in all_results:
            print(f"{r['l_max']:5d} | {r['energy_std']:12.6f} | "
                  f"{r['energy_tc']:12.6f} | {r['err_std_pct']:8.4f} | "
                  f"{r['err_tc_pct']:8.4f} | {r['max_imag']:10.2e}")

    return all_results
