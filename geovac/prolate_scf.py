"""
Prolate spheroidal SCF and CI for two-electron diatomics.

Implements three levels of orbital optimization:
1. Eckart variational SCF: single Z_eff parameter scaling
2. Relaxed-orbital CI: CI with Z_eff-optimized orbitals
3. Full grid-based SCF: Fock operator on 2D (xi,eta) grid

The grid SCF solves the Fock equation directly on the prolate spheroidal
grid, allowing the orbital to redistribute density beyond what uniform
Z_eff scaling can capture (polarization, shape flexibility).

References:
  - Eckart, Phys. Rev. 36, 878 (1930)
  - Bates, Ledsham & Stewart, Phil. Trans. A 246, 215 (1953)
"""

import numpy as np
from scipy.special import ellipk
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
from typing import Dict, Tuple
from math import factorial
import time

from geovac.prolate_spheroidal_lattice import (
    ProlateSpheroidalLattice,
    fit_spectroscopic_constants,
)
from geovac.molecular_sturmian import _angular_sep_const


# ============================================================
# Orbital generation on quadrature grid
# ============================================================

def get_orbital_on_grid(
    R: float,
    Z_A: float = 1.0,
    Z_B: float = 1.0,
    n_angular: int = 0,
    N_xi_solve: int = 5000,
    N_xi_grid: int = 80,
    N_eta_grid: int = 80,
    xi_max_grid: float = 12.0,
    m: int = 0,
) -> Dict:
    """Solve for a one-electron orbital and evaluate on quadrature grid.

    Supports fractional Z via internal parameter override. The solver
    only uses a = R*(Z_A+Z_B) and b = R*(Z_B-Z_A), so float Z works.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    Z_A, Z_B : float
        Nuclear charges (can be fractional for Z_eff optimization).
    n_angular : int
        Angular eigenvalue index (0=ground, 1=first excited, ...).
    N_xi_solve : int
        Radial grid points for eigenvalue solve.
    N_xi_grid, N_eta_grid : int
        Quadrature grid dimensions for wavefunction evaluation.
    xi_max_grid : float
        Maximum xi for quadrature grid.
    m : int
        Azimuthal quantum number (0 for sigma states).

    Returns
    -------
    dict with keys: E_elec, c2, A, xi, eta, w_xi, w_eta, psi, F_xi,
    G_eta, R, Z_A, Z_B, n_angular.
    """
    lat = ProlateSpheroidalLattice(
        R=R, Z_A=max(1, int(round(Z_A))), Z_B=max(1, int(round(Z_B))),
        N_xi=N_xi_solve, xi_max=max(25.0, R * 3),
        m=m, n_angular=n_angular,
    )
    # Override for fractional Z support
    lat._a = R * (Z_A + Z_B)
    lat._b = R * (Z_B - Z_A)

    E_elec, c2, A = lat.solve()
    c = np.sqrt(max(c2, 1e-15))

    # --- Radial wavefunction F(xi) on solve grid ---
    N = N_xi_solve
    xi_min_s = 1.0 + 5e-4
    xi_max_s = max(25.0, R * 3)
    h_s = (xi_max_s - xi_min_s) / (N + 1)
    xi_s = xi_min_s + (np.arange(N) + 1) * h_s

    a_param = R * (Z_A + Z_B)
    xi2_1 = xi_s ** 2 - 1
    q = A + a_param * xi_s - c2 * xi_s ** 2
    if m != 0:
        q -= m ** 2 / xi2_1

    p_plus = (xi_s + h_s / 2) ** 2 - 1
    p_minus = (xi_s - h_s / 2) ** 2 - 1

    diag = -(p_plus + p_minus) / h_s ** 2 + q
    off = p_plus[:-1] / h_s ** 2

    if m == 0:
        diag[0] = -p_plus[0] / h_s ** 2 + q[0]

    from scipy.linalg import eigh_tridiagonal
    evals, evecs = eigh_tridiagonal(diag, off)
    idx_sorted = np.argsort(evals)[::-1]
    F_raw = evecs[:, idx_sorted[0]]
    if F_raw[0] < 0:
        F_raw = -F_raw

    # --- Angular eigenvector (Legendre basis) ---
    m_abs = abs(m)
    n_basis = 50
    r_vals = np.arange(m_abs, m_abs + n_basis, dtype=float)
    norms = np.array(
        [2.0 / (2 * r + 1) * factorial(int(r + m_abs)) / factorial(int(r - m_abs))
         for r in r_vals]
    )

    nu_mat = np.zeros((n_basis, n_basis))
    for i in range(n_basis - 1):
        r = r_vals[i]
        val = ((r - m_abs + 1) * np.sqrt(norms[i + 1])
               / ((2 * r + 1) * np.sqrt(norms[i])))
        nu_mat[i + 1, i] = val
        nu_mat[i, i + 1] = val

    b_param = R * (Z_B - Z_A)
    H_ang = np.diag(-r_vals * (r_vals + 1))
    H_ang += c ** 2 * (nu_mat @ nu_mat)
    H_ang += b_param * nu_mat

    evals_ang, evecs_ang = np.linalg.eigh(H_ang)
    idx_ang = n_basis - 1 - n_angular
    coeffs_ang = evecs_ang[:, idx_ang]

    # --- Build quadrature grids ---
    eta, w_eta = np.polynomial.legendre.leggauss(N_eta_grid)
    u_gl, w_u_gl = np.polynomial.legendre.leggauss(N_xi_grid)
    t = (u_gl + 1) / 2
    xi = 1.0 + (xi_max_grid - 1.0) * t ** 2
    w_xi = w_u_gl * (xi_max_grid - 1.0) * t

    # --- Evaluate wavefunctions on quadrature grid ---
    F_xi = np.interp(xi, xi_s, F_raw, left=F_raw[0], right=0.0)

    from scipy.special import lpmv
    G_eta = np.zeros(N_eta_grid)
    for j, r in enumerate(r_vals):
        P_r = lpmv(m_abs, int(r), eta)
        G_eta += coeffs_ang[j] * P_r / np.sqrt(norms[j])

    # --- Build and normalize 2D wavefunction ---
    psi_2d = np.outer(F_xi, G_eta)

    XI, ETA = np.meshgrid(xi, eta, indexing='ij')
    J = (R / 2) ** 3 * (XI ** 2 - ETA ** 2)
    W_XI, W_ETA = np.meshgrid(w_xi, w_eta, indexing='ij')

    norm_sq = np.sum(psi_2d ** 2 * J * W_XI * W_ETA) * 2 * np.pi
    psi_2d /= np.sqrt(norm_sq)

    norm_check = np.sum(psi_2d ** 2 * J * W_XI * W_ETA) * 2 * np.pi
    assert abs(norm_check - 1.0) < 0.02, f"Normalization failed: {norm_check:.4f}"

    return {
        'E_elec': E_elec,
        'c2': c2,
        'A': A,
        'xi': xi,
        'eta': eta,
        'w_xi': w_xi,
        'w_eta': w_eta,
        'psi': psi_2d,
        'F_xi': F_xi,
        'G_eta': G_eta,
        'R': R,
        'Z_A': Z_A,
        'Z_B': Z_B,
        'n_angular': n_angular,
    }


# ============================================================
# Two-electron integrals via azimuthal averaging
# ============================================================

def compute_vee_integral(
    orb_a: Dict, orb_b: Dict, orb_c: Dict, orb_d: Dict,
) -> float:
    """Compute (ab|cd) two-electron integral via azimuthal averaging.

    Uses complete elliptic integral K for the phi-averaged 1/r12 kernel.
    Convention: chemists' notation (ab|cd) = [ac|bd] physicists'.
    """
    R = orb_a['R']
    xi = orb_a['xi']
    eta = orb_a['eta']
    w_xi = orb_a['w_xi']
    w_eta = orb_a['w_eta']

    rho_ac = orb_a['psi'] * orb_c['psi']
    rho_bd = orb_b['psi'] * orb_d['psi']

    XI, ETA = np.meshgrid(xi, eta, indexing='ij')
    rho_cyl = (R / 2) * np.sqrt(np.maximum((XI ** 2 - 1) * (1 - ETA ** 2), 0.0))
    z_cyl = (R / 2) * XI * ETA
    J = (R / 2) ** 3 * (XI ** 2 - ETA ** 2)
    W_XI, W_ETA = np.meshgrid(w_xi, w_eta, indexing='ij')
    W_2d = W_XI * W_ETA

    n_pts = len(xi) * len(eta)
    rho_ac_flat = (rho_ac * J * W_2d).flatten()
    rho_bd_flat = (rho_bd * J * W_2d).flatten()
    rho_flat = rho_cyl.flatten()
    z_flat = z_cyl.flatten()

    result = 0.0
    batch_size = 500
    n_batches = (n_pts + batch_size - 1) // batch_size

    for i_batch in range(n_batches):
        i_start = i_batch * batch_size
        i_end = min(i_start + batch_size, n_pts)

        rho1 = rho_flat[i_start:i_end, np.newaxis]
        z1 = z_flat[i_start:i_end, np.newaxis]
        rho2 = rho_flat[np.newaxis, :]
        z2 = z_flat[np.newaxis, :]

        dz = z1 - z2
        s = (rho1 + rho2) ** 2 + dz ** 2

        with np.errstate(divide='ignore', invalid='ignore'):
            k2 = np.where(s > 1e-30, 4 * rho1 * rho2 / s, 0.0)
            k2 = np.clip(k2, 0, 1 - 1e-15)

        s_safe = np.maximum(s, 1e-30)
        K_vals = ellipk(k2)
        coulomb_kernel = 8 * np.pi * K_vals / np.sqrt(s_safe)
        coulomb_kernel = np.where(s < 1e-30, 0.0, coulomb_kernel)

        w_ac = rho_ac_flat[i_start:i_end]
        w_bd = rho_bd_flat

        result += np.sum(w_ac[:, np.newaxis] * coulomb_kernel * w_bd[np.newaxis, :])

    return result


# ============================================================
# Nuclear attraction expectation value
# ============================================================

def compute_vnuc_expectation(orb: Dict) -> float:
    """Compute <phi|1/r_A + 1/r_B|phi> on the quadrature grid.

    Uses the identity: (1/r_A + 1/r_B) * J = R^2 * xi / 2
    where J = (R/2)^3 * (xi^2 - eta^2) is the volume-element Jacobian.
    """
    R = orb['R']
    xi = orb['xi']
    w_xi = orb['w_xi']
    w_eta = orb['w_eta']
    psi = orb['psi']

    XI = xi[:, np.newaxis]
    vnuc_J = R ** 2 * XI / 2

    W_XI = w_xi[:, np.newaxis]
    W_ETA = w_eta[np.newaxis, :]

    return 2 * np.pi * np.sum(psi ** 2 * vnuc_J * W_XI * W_ETA)


# ============================================================
# Coulomb potential on grid (for full grid SCF)
# ============================================================

def compute_coulomb_potential(
    psi_2d: np.ndarray,
    xi: np.ndarray,
    eta: np.ndarray,
    w_xi: np.ndarray,
    w_eta: np.ndarray,
    R: float,
) -> np.ndarray:
    """Compute V_J(xi,eta) = integral rho(r')/|r-r'| d3r' on grid.

    For closed-shell doubly-occupied orbital: rho = 2|psi|^2.
    The azimuthal integral is evaluated analytically using elliptic K.

    Parameters
    ----------
    psi_2d : array (N_xi, N_eta)
        Normalized orbital on (xi, eta) quadrature grid.
    xi, eta : arrays
        Grid coordinates.
    w_xi, w_eta : arrays
        Quadrature weights.
    R : float
        Internuclear distance.

    Returns
    -------
    V_J : array (N_xi, N_eta)
        Coulomb potential at each grid point.
    """
    N_xi = len(xi)
    N_eta = len(eta)

    # Density: rho = 2|psi|^2 (closed-shell)
    # In the integral, rho appears with volume element J * dxi * deta
    # (the 2*pi from phi is absorbed into the kernel)
    XI, ETA = np.meshgrid(xi, eta, indexing='ij')
    J = (R / 2) ** 3 * (XI ** 2 - ETA ** 2)
    W_XI, W_ETA = np.meshgrid(w_xi, w_eta, indexing='ij')

    # Weighted density: rho * J * w_xi * w_eta
    # rho = 2|psi|^2 / (2*pi) * (2*pi) = 2|psi|^2
    # (the 1/2pi from m=0 normalization × 2pi from azimuthal density integral = 1)
    # Actually: psi_2D is normalized so that int |psi_2D|^2 J dxi deta * 2pi = 1
    # So |psi(r)|^2 = |psi_2D|^2 / (2*pi)
    # rho = 2|psi|^2 = 2|psi_2D|^2/(2*pi) = |psi_2D|^2/pi
    # V_J = int rho(r') K_kernel(r,r') J' dxi' deta'
    #     = (1/pi) int |psi_2D|^2 * 4K(k)/sqrt(s) * J' dxi' deta'
    rho_weighted = (1.0 / np.pi) * psi_2d ** 2 * J * W_XI * W_ETA

    # Cylindrical coordinates for Coulomb kernel
    rho_cyl = (R / 2) * np.sqrt(np.maximum((XI ** 2 - 1) * (1 - ETA ** 2), 0.0))
    z_cyl = (R / 2) * XI * ETA

    rho_w_flat = rho_weighted.flatten()
    rho_flat = rho_cyl.flatten()
    z_flat = z_cyl.flatten()
    n_pts = N_xi * N_eta

    V_J = np.zeros(n_pts)
    batch_size = 500
    n_batches = (n_pts + batch_size - 1) // batch_size

    for i_batch in range(n_batches):
        i_start = i_batch * batch_size
        i_end = min(i_start + batch_size, n_pts)

        rho1 = rho_flat[i_start:i_end, np.newaxis]  # (batch, 1)
        z1 = z_flat[i_start:i_end, np.newaxis]

        rho2 = rho_flat[np.newaxis, :]  # (1, n_pts)
        z2 = z_flat[np.newaxis, :]

        dz = z1 - z2
        s = (rho1 + rho2) ** 2 + dz ** 2

        with np.errstate(divide='ignore', invalid='ignore'):
            k2 = np.where(s > 1e-30, 4 * rho1 * rho2 / s, 0.0)
            k2 = np.clip(k2, 0, 1 - 1e-15)

        s_safe = np.maximum(s, 1e-30)
        K_vals = ellipk(k2)
        kernel = 4.0 * K_vals / np.sqrt(s_safe)
        kernel = np.where(s < 1e-30, 0.0, kernel)

        # V_J[i] = sum_j rho_weighted[j] * kernel[i,j]
        V_J[i_start:i_end] = np.dot(kernel, rho_w_flat)

    return V_J.reshape(N_xi, N_eta)


# ============================================================
# 2D finite-difference Hamiltonian for grid SCF
# ============================================================

def _build_2d_hamiltonian(
    R: float,
    Z_A: float,
    Z_B: float,
    N_xi: int,
    N_eta: int,
    xi_max: float,
    V_ext: np.ndarray = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Build 2D FD Hamiltonian for the prolate spheroidal Fock equation.

    Uses a coordinate transformation xi = 1 + (xi_max-1)*t^2 to
    concentrate grid points near xi=1 where the wavefunction peaks.

    Solves the generalized eigenvalue problem A psi = c^2 B psi
    where c^2 = -R^2 E / 2.

    After multiplying the Schrodinger equation by -(R^2/2)(xi^2-eta^2)
    and by dxi/dt (the Jacobian of the xi mapping), we get a symmetric
    generalized eigenvalue problem in (t, eta) space.

    Parameters
    ----------
    R, Z_A, Z_B : float
        Internuclear distance and nuclear charges.
    N_xi, N_eta : int
        Grid dimensions.
    xi_max : float
        Maximum xi value.
    V_ext : array (N_xi, N_eta) or None
        External potential (e.g., V_J/2 for SCF). If None, solves H2+.

    Returns
    -------
    A : array (N_xi*N_eta, N_xi*N_eta)
        LHS matrix.
    B : array (N_xi*N_eta,)
        Diagonal of RHS weight matrix.
    xi_grid, eta_grid : arrays
        Grid coordinates (in physical xi, eta).
    """
    a_param = R * (Z_A + Z_B)
    b_param = R * (Z_B - Z_A)

    # --- Grids ---
    # t-space: uniform in [0, 1], maps to xi via xi = 1 + (xi_max-1)*t^2
    h_t = 1.0 / (N_xi + 1)
    t = (np.arange(N_xi) + 1) * h_t  # t in (0, 1)
    xi = 1.0 + (xi_max - 1.0) * t ** 2
    xi_prime = 2.0 * (xi_max - 1.0) * t  # dxi/dt

    # eta: uniform grid
    eta_min = -1.0 + 5e-4
    eta_max_val = 1.0 - 5e-4
    h_eta = (eta_max_val - eta_min) / (N_eta + 1)
    eta = eta_min + (np.arange(N_eta) + 1) * h_eta

    n_total = N_xi * N_eta

    # --- Build L_t: d/dt[P(t) d/dt] in t-space ---
    # P(t) = p(xi(t)) / xi'(t) = (xi^2 - 1) / xi'
    # At t->0: P -> t + O(t^3), so P(0) = 0 (natural Neumann BC)
    P = (xi ** 2 - 1.0) / xi_prime

    # P at half-points: P_{i+1/2} = (P_i + P_{i+1}) / 2
    P_half_plus = np.zeros(N_xi)
    P_half_minus = np.zeros(N_xi)
    for i in range(N_xi - 1):
        P_half_plus[i] = 0.5 * (P[i] + P[i + 1])
    P_half_plus[N_xi - 1] = 0.5 * P[N_xi - 1]  # Dirichlet at t=1

    # P at left boundary: P(0) = 0, so P_{-1/2} = P[0]/2
    P_half_minus[0] = 0.5 * P[0]  # BC: P(t=0) ≈ 0
    for i in range(1, N_xi):
        P_half_minus[i] = 0.5 * (P[i - 1] + P[i])

    diag_t = -(P_half_plus + P_half_minus) / h_t ** 2
    # Neumann at t=0: zero flux, P_{-1/2} ≈ 0
    diag_t[0] = -P_half_plus[0] / h_t ** 2

    upper_t = P_half_plus[:-1] / h_t ** 2  # connects i to i+1
    # lower_t = P_half_minus[1:] / h_t**2  # connects i to i-1 (= upper by symmetry)

    # --- Build L_eta: d/deta[(1-eta^2) d/deta] ---
    q_plus_eta = 1 - (eta + h_eta / 2) ** 2
    q_minus_eta = 1 - (eta - h_eta / 2) ** 2

    diag_eta = -(q_plus_eta + q_minus_eta) / h_eta ** 2
    upper_eta = q_plus_eta[:-1] / h_eta ** 2

    # --- Assemble full 2D matrix using Kronecker structure ---
    # Index: k = i * N_eta + j
    # The equation in (t, eta) space, multiplied by xi'(t):
    #   L_t psi + xi'(t) L_eta psi + xi'(t)[a*xi + b*eta - ...] psi
    #     = c^2 xi'(t)(xi^2-eta^2) psi
    #
    # So A = L_t (x) I + diag(xi') (x) L_eta + diag(V_total)
    # and B = diag(xi'(t) * (xi^2 - eta^2))

    A = np.zeros((n_total, n_total))

    for i in range(N_xi):
        xp = xi_prime[i]  # xi'(t_i)

        for j in range(N_eta):
            k = i * N_eta + j

            # Diagonal: L_t + xi' * L_eta + xi' * potential
            A[k, k] = (diag_t[i]
                        + xp * diag_eta[j]
                        + xp * (a_param * xi[i] + b_param * eta[j]))

            # V_ext contribution
            if V_ext is not None:
                w = xi[i] ** 2 - eta[j] ** 2
                A[k, k] -= xp * (R ** 2 / 2) * w * V_ext[i, j]

            # L_t off-diagonal (xi direction): i to i+1
            if i < N_xi - 1:
                k_right = (i + 1) * N_eta + j
                A[k, k_right] = upper_t[i]
                A[k_right, k] = upper_t[i]

            # xi' * L_eta off-diagonal (eta direction): j to j+1
            if j < N_eta - 1:
                k_up = i * N_eta + (j + 1)
                A[k, k_up] = xp * upper_eta[j]
                A[k_up, k] = xp * upper_eta[j]

    # Weight matrix B = xi'(t) * (xi^2 - eta^2)
    B_diag = np.zeros(n_total)
    for i in range(N_xi):
        for j in range(N_eta):
            k = i * N_eta + j
            B_diag[k] = xi_prime[i] * (xi[i] ** 2 - eta[j] ** 2)

    return A, B_diag, xi, eta


def _solve_2d_eigenvalue(
    A: np.ndarray,
    B_diag: np.ndarray,
    R: float,
) -> Tuple[float, np.ndarray]:
    """Solve generalized eigenvalue problem A psi = c^2 B psi.

    Returns the ground state energy and eigenvector.
    Uses B^{-1/2} transformation for symmetric standard form.
    """
    B_sqrt = np.sqrt(B_diag)
    B_inv_sqrt = 1.0 / B_sqrt

    # Transform: H_eff = B^{-1/2} A B^{-1/2}
    n = len(B_diag)
    H_eff = A * (B_inv_sqrt[:, np.newaxis] * B_inv_sqrt[np.newaxis, :])

    # Symmetrize (numerical noise)
    H_eff = 0.5 * (H_eff + H_eff.T)

    evals, evecs = eigh(H_eff)

    # Largest eigenvalue = largest c^2 = most bound state
    c2 = evals[-1]
    phi = evecs[:, -1]

    # Transform back: psi = B^{-1/2} phi
    psi = B_inv_sqrt * phi

    E_elec = -2.0 * c2 / R ** 2

    return E_elec, psi, c2


# ============================================================
# Full grid-based SCF (Part 1)
# ============================================================

def grid_scf(
    R: float,
    Z_A: float = 1.0,
    Z_B: float = 1.0,
    N_xi_fd: int = 40,
    N_eta_fd: int = 40,
    xi_max_fd: float = 15.0,
    N_xi_solve: int = 5000,
    max_iter: int = 30,
    tol: float = 1e-5,
    damping: float = 0.5,
    verbose: bool = True,
) -> Dict:
    """Full grid-based SCF for two-electron diatomics.

    Solves the Fock equation F phi = eps phi on a 2D (xi, eta) FD grid,
    where F = T + V_nuc + V_J/2 (closed-shell, one spatial orbital).

    The Coulomb potential V_J is computed on the FD grid via azimuthal
    averaging (elliptic K kernel). The orbital is updated iteratively
    with damping for stability.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    Z_A, Z_B : float
        Nuclear charges.
    N_xi_fd, N_eta_fd : int
        FD grid dimensions for eigenvalue solve.
    xi_max_fd : float
        Maximum xi for FD grid.
    N_xi_solve : int
        Radial grid for initial H2+ orbital.
    max_iter : int
        Maximum SCF iterations.
    tol : float
        Energy convergence threshold (Ha).
    damping : float
        Orbital mixing parameter (0 < alpha < 1). New orbital mixed as
        psi = damping * psi_new + (1-damping) * psi_old.
    verbose : bool
        Print iteration details.

    Returns
    -------
    dict with keys: E_total, E_elec, J_gg, eps, psi_2d, xi, eta,
    n_iter, converged, R, Z_A, Z_B.
    """
    t0 = time.time()
    a_param = R * (Z_A + Z_B)
    b_param = R * (Z_B - Z_A)
    V_NN = Z_A * Z_B / R

    # --- Step 1: Get initial orbital from H2+ solver ---
    # Solve H2+ on the fine grid, then interpolate onto the FD grid
    orb0 = get_orbital_on_grid(
        R=R, Z_A=Z_A, Z_B=Z_B, n_angular=0,
        N_xi_solve=N_xi_solve,
        N_xi_grid=60, N_eta_grid=60,
        xi_max_grid=xi_max_fd,
    )

    # --- Step 2: Set up FD grid (same mapping as _build_2d_hamiltonian) ---
    h_t = 1.0 / (N_xi_fd + 1)
    t_fd = (np.arange(N_xi_fd) + 1) * h_t
    xi_fd = 1.0 + (xi_max_fd - 1.0) * t_fd ** 2

    eta_min = -1.0 + 5e-4
    eta_max_val = 1.0 - 5e-4
    h_eta = (eta_max_val - eta_min) / (N_eta_fd + 1)
    eta_fd = eta_min + (np.arange(N_eta_fd) + 1) * h_eta

    # Interpolate initial orbital onto FD grid
    from scipy.interpolate import RegularGridInterpolator
    interp = RegularGridInterpolator(
        (orb0['xi'], orb0['eta']), orb0['psi'],
        method='linear', bounds_error=False, fill_value=0.0,
    )
    XI_fd, ETA_fd = np.meshgrid(xi_fd, eta_fd, indexing='ij')
    pts = np.column_stack([XI_fd.flatten(), ETA_fd.flatten()])
    psi_fd = interp(pts).reshape(N_xi_fd, N_eta_fd)

    # Normalize on FD grid (quadrature weights for non-uniform xi)
    J_fd = (R / 2) ** 3 * (XI_fd ** 2 - ETA_fd ** 2)
    # Quadrature weights: dxi = xi'(t) dt, so w_xi_i = xi'_i * h_t
    xi_prime_fd = 2.0 * (xi_max_fd - 1.0) * t_fd
    w_xi_fd = xi_prime_fd * h_t
    w_eta_fd = np.full(N_eta_fd, h_eta)
    W_XI_fd, W_ETA_fd = np.meshgrid(w_xi_fd, w_eta_fd, indexing='ij')

    norm_sq = np.sum(psi_fd ** 2 * J_fd * W_XI_fd * W_ETA_fd) * 2 * np.pi
    if norm_sq > 0:
        psi_fd /= np.sqrt(norm_sq)

    # --- SCF iteration ---
    eps_old = 0.0
    converged = False

    if verbose:
        print(f"\n  Grid SCF: R={R:.3f}, Z_A={Z_A:.2f}, Z_B={Z_B:.2f}")
        print(f"  FD grid: {N_xi_fd}x{N_eta_fd}, xi_max={xi_max_fd:.1f}")
        print(f"  {'Iter':>4s} {'E_elec':>12s} {'J_gg':>10s} {'dE':>12s}")

    for iteration in range(max_iter):
        # Compute V_J on FD grid
        V_J = compute_coulomb_potential(
            psi_fd, xi_fd, eta_fd, w_xi_fd, w_eta_fd, R,
        )

        # V_ext for Fock operator: V_J / 2
        V_ext = V_J / 2.0

        # Build and solve 2D Hamiltonian with V_ext
        A_mat, B_diag, _, _ = _build_2d_hamiltonian(
            R, Z_A, Z_B, N_xi_fd, N_eta_fd, xi_max_fd, V_ext,
        )

        E_elec, psi_new_flat, c2 = _solve_2d_eigenvalue(A_mat, B_diag, R)

        # Reshape and normalize
        psi_new = psi_new_flat.reshape(N_xi_fd, N_eta_fd)
        # Ensure positive at origin
        if psi_new[0, N_eta_fd // 2] < 0:
            psi_new = -psi_new
        norm_sq = np.sum(psi_new ** 2 * J_fd * W_XI_fd * W_ETA_fd) * 2 * np.pi
        if norm_sq > 0:
            psi_new /= np.sqrt(norm_sq)

        # Compute J_gg for energy
        # J_gg = (1/2) <phi|V_J|phi> where V_J uses rho=2|phi|^2
        # <phi|V_J|phi> = int |phi|^2 V_J J dxi deta * 2pi
        J_gg_val = np.sum(
            psi_new ** 2 * V_J * J_fd * W_XI_fd * W_ETA_fd
        ) * 2 * np.pi
        J_gg = J_gg_val / 2.0  # V_J is from rho=2|phi|^2, so <|V_J|> = 2*J_gg

        # Fock eigenvalue: eps = h + J_gg where h = <T+V_nuc>
        # E_elec from the eigenvalue is the one-electron energy WITH V_ext
        # eps = E_from_fock = h + <V_J/2> = h + J_gg (since <V_J> = 2*J_gg)
        eps = E_elec  # This is the Fock eigenvalue

        # Total HF energy: E = 2*h + J_gg + V_NN = 2*(eps - J_gg) + J_gg + V_NN
        # = 2*eps - J_gg + V_NN
        E_HF = 2 * eps - J_gg + V_NN

        dE = abs(eps - eps_old)

        if verbose:
            print(f"  {iteration:4d} {E_HF - V_NN:12.6f} {J_gg:10.6f} {dE:12.2e}")

        if dE < tol and iteration > 0:
            converged = True
            if verbose:
                print(f"  Converged after {iteration + 1} iterations")
            break

        eps_old = eps

        # Mix orbitals with damping
        if iteration > 0:
            psi_fd = damping * psi_new + (1 - damping) * psi_fd
            # Re-normalize
            norm_sq = np.sum(
                psi_fd ** 2 * J_fd * W_XI_fd * W_ETA_fd
            ) * 2 * np.pi
            if norm_sq > 0:
                psi_fd /= np.sqrt(norm_sq)
        else:
            psi_fd = psi_new.copy()

    dt = time.time() - t0

    if verbose:
        print(f"\n  E_HF(grid SCF) = {E_HF:.6f} Ha")
        print(f"  E_elec = {E_HF - V_NN:.6f}, J_gg = {J_gg:.6f}")
        print(f"  Time: {dt:.1f}s")

    return {
        'E_total': E_HF,
        'E_elec': E_HF - V_NN,
        'J_gg': J_gg,
        'eps': eps,
        'psi_2d': psi_fd,
        'xi': xi_fd,
        'eta': eta_fd,
        'n_iter': iteration + 1,
        'converged': converged,
        'R': R,
        'Z_A': Z_A,
        'Z_B': Z_B,
        'time': dt,
    }


# ============================================================
# Eckart variational SCF (Part 3 support)
# ============================================================

def eckart_scf_energy(
    R: float,
    Z_eff: float,
    N_xi_solve: int = 5000,
    N_grid: int = 50,
    xi_max_grid: float = 15.0,
) -> Dict:
    """Compute H2 HF energy using orbital generated at (Z_eff, Z_eff).

    The Eckart variational method: solve H2+ at R_eff = R*Z_eff with
    Z=1, then evaluate the physical energy using scaling relations.

    h_phys = Z_eff^2 * <T(R_eff)> - Z_eff * <V_nuc(R_eff)>
    E_HF = 2*h_phys + J_gg + 1/R
    """
    R_eff = R * Z_eff

    try:
        orb = get_orbital_on_grid(
            R=R_eff, Z_A=1.0, Z_B=1.0, n_angular=0,
            N_xi_solve=N_xi_solve,
            N_xi_grid=N_grid, N_eta_grid=N_grid,
            xi_max_grid=xi_max_grid,
        )
    except Exception as e:
        return {'R': R, 'Z_eff': Z_eff, 'E_HF': np.nan, 'error': str(e)}

    eps_zeff = orb['E_elec']

    # Re-normalize with physical R for integrals
    orb_phys = dict(orb)
    orb_phys['R'] = R
    XI, ETA = np.meshgrid(orb['xi'], orb['eta'], indexing='ij')
    J_phys = (R / 2) ** 3 * (XI ** 2 - ETA ** 2)
    W_XI, W_ETA = np.meshgrid(orb['w_xi'], orb['w_eta'], indexing='ij')
    norm_sq = np.sum(orb['psi'] ** 2 * J_phys * W_XI * W_ETA) * 2 * np.pi
    orb_phys['psi'] = orb['psi'] / np.sqrt(norm_sq)

    # Scaling relation: <T(R)> = Z_eff^2 * <T(R_eff)>, <V(R)> = Z_eff * <V(R_eff)>
    vnuc_reff = compute_vnuc_expectation(orb)
    T_reff = eps_zeff + vnuc_reff  # <T> = eps - <V_nuc> = eps + <1/r_A + 1/r_B>

    h_phys = Z_eff ** 2 * T_reff - Z_eff * vnuc_reff

    # V_ee with physical R
    J_gg = compute_vee_integral(orb_phys, orb_phys, orb_phys, orb_phys)

    E_HF = 2 * h_phys + J_gg + 1.0 / R

    return {
        'R': R,
        'Z_eff': Z_eff,
        'E_HF': E_HF,
        'h_gg': h_phys,
        'J_gg': J_gg,
        'eps_zeff': eps_zeff,
        'T_reff': T_reff,
        'vnuc_reff': vnuc_reff,
    }


def optimize_zeff_eckart(
    R: float,
    N_xi_solve: int = 5000,
    N_grid: int = 50,
    xi_max_grid: float = 15.0,
    verbose: bool = False,
) -> Dict:
    """Find optimal Z_eff that minimizes E_HF(Z_eff) at given R."""
    t0 = time.time()

    # Coarse scan
    zeff_vals = np.linspace(0.7, 1.5, 17)
    energies = []
    for z in zeff_vals:
        res = eckart_scf_energy(R, z, N_xi_solve, N_grid, xi_max_grid)
        energies.append(res['E_HF'])
    energies = np.array(energies)

    valid = ~np.isnan(energies)
    if np.sum(valid) < 3:
        return {'R': R, 'Z_eff_opt': np.nan, 'E_HF_opt': np.nan}

    idx_min = np.nanargmin(energies)

    # Refine
    if 0 < idx_min < len(zeff_vals) - 1:
        z_lo = zeff_vals[max(0, idx_min - 1)]
        z_hi = zeff_vals[min(len(zeff_vals) - 1, idx_min + 1)]
    else:
        z_lo, z_hi = zeff_vals[valid][0], zeff_vals[valid][-1]

    def obj(z):
        res = eckart_scf_energy(R, z, N_xi_solve, N_grid, xi_max_grid)
        return res['E_HF'] if not np.isnan(res['E_HF']) else 1e10

    opt = minimize_scalar(obj, bounds=(z_lo, z_hi), method='bounded',
                          options={'xatol': 0.001})
    z_opt = opt.x
    res_opt = eckart_scf_energy(R, z_opt, N_xi_solve, N_grid, xi_max_grid)
    res_frozen = eckart_scf_energy(R, 1.0, N_xi_solve, N_grid, xi_max_grid)

    dt = time.time() - t0

    if verbose:
        print(f"  Z_eff*={z_opt:.4f}, E_HF={res_opt['E_HF']:.6f} "
              f"(frozen={res_frozen['E_HF']:.6f}), J_gg={res_opt['J_gg']:.4f} [{dt:.0f}s]")

    return {
        'R': R,
        'Z_eff_opt': z_opt,
        'E_HF_opt': res_opt['E_HF'],
        'E_HF_frozen': res_frozen['E_HF'],
        'J_gg_opt': res_opt['J_gg'],
        'J_gg_frozen': res_frozen['J_gg'],
        'h_gg_opt': res_opt['h_gg'],
        'h_gg_frozen': res_frozen['h_gg'],
        'time': dt,
    }


# ============================================================
# Relaxed-orbital CI (Part 3)
# ============================================================

def relaxed_orbital_ci(
    R: float,
    N_xi_solve: int = 5000,
    N_grid: int = 50,
    xi_max_grid: float = 15.0,
    verbose: bool = True,
) -> Dict:
    """CI with Z_eff-optimized orbitals for H2.

    1. Optimize Z_eff via Eckart variational
    2. Generate sigma_g and sigma_u orbitals at Z_eff*
    3. Build and diagonalize 2x2 CI matrix

    Returns
    -------
    dict with E_total, E_HF_scf, E_HF_frozen, Z_eff_opt, J_gg, etc.
    """
    t0 = time.time()

    # Step 1: Optimize Z_eff
    opt = optimize_zeff_eckart(R, N_xi_solve, N_grid, xi_max_grid, verbose=verbose)
    if np.isnan(opt['Z_eff_opt']):
        return {'R': R, 'E_total': np.nan}

    z_opt = opt['Z_eff_opt']
    R_eff = R * z_opt

    if verbose:
        print(f"  Generating orbitals at Z_eff={z_opt:.4f} (R_eff={R_eff:.4f})")

    # Step 2: Generate both orbitals at Z_eff*
    orbs = {}
    for label, n_ang in [('g', 0), ('u', 1)]:
        orb = get_orbital_on_grid(
            R=R_eff, Z_A=1.0, Z_B=1.0, n_angular=n_ang,
            N_xi_solve=N_xi_solve,
            N_xi_grid=N_grid, N_eta_grid=N_grid,
            xi_max_grid=xi_max_grid,
        )
        # Re-normalize with physical R
        XI, ETA = np.meshgrid(orb['xi'], orb['eta'], indexing='ij')
        J_phys = (R / 2) ** 3 * (XI ** 2 - ETA ** 2)
        W_XI, W_ETA = np.meshgrid(orb['w_xi'], orb['w_eta'], indexing='ij')
        norm_sq = np.sum(orb['psi'] ** 2 * J_phys * W_XI * W_ETA) * 2 * np.pi
        orb['psi'] = orb['psi'] / np.sqrt(norm_sq)
        orb['R'] = R
        orbs[label] = orb

    # Step 3: Compute V_ee integrals
    J_gg = compute_vee_integral(orbs['g'], orbs['g'], orbs['g'], orbs['g'])
    J_uu = compute_vee_integral(orbs['u'], orbs['u'], orbs['u'], orbs['u'])
    V_gguu = compute_vee_integral(orbs['g'], orbs['g'], orbs['u'], orbs['u'])

    # Step 4: One-electron energies via scaling
    h_g = opt['h_gg_opt']

    # For sigma_u at Z_eff:
    orb_u_reff = get_orbital_on_grid(
        R=R_eff, Z_A=1.0, Z_B=1.0, n_angular=1,
        N_xi_solve=N_xi_solve,
        N_xi_grid=N_grid, N_eta_grid=N_grid,
        xi_max_grid=xi_max_grid,
    )
    eps_u_zeff = orb_u_reff['E_elec']
    vnuc_u_reff = compute_vnuc_expectation(orb_u_reff)
    T_u_reff = eps_u_zeff + vnuc_u_reff
    h_u = z_opt ** 2 * T_u_reff - z_opt * vnuc_u_reff

    # Step 5: Build CI matrix
    H_CI = np.array([
        [2 * h_g + J_gg, V_gguu],
        [V_gguu, 2 * h_u + J_uu],
    ])

    evals, evecs = eigh(H_CI)
    E_CI_elec = evals[0]
    V_NN = 1.0 / R
    E_total = E_CI_elec + V_NN
    E_HF_scf = 2 * h_g + J_gg + V_NN

    dt = time.time() - t0

    if verbose:
        c_g, c_u = evecs[0, 0], evecs[1, 0]
        print(f"  h_g={h_g:.6f}, h_u={h_u:.6f}")
        print(f"  J_gg={J_gg:.6f}, J_uu={J_uu:.6f}, V_gguu={V_gguu:.6f}")
        print(f"  CI: c_g={c_g:.4f}, c_u={c_u:.4f}")
        print(f"  E_HF(SCF) = {E_HF_scf:.6f}, E_CI(SCF) = {E_total:.6f}")
        print(f"  Time: {dt:.1f}s")

    return {
        'R': R,
        'E_total': E_total,
        'E_HF_scf': E_HF_scf,
        'E_HF_frozen': opt['E_HF_frozen'],
        'Z_eff_opt': z_opt,
        'J_gg': J_gg,
        'J_uu': J_uu,
        'V_gguu': V_gguu,
        'h_gg': h_g,
        'h_uu': h_u,
        'c_g': evecs[0, 0],
        'c_u': evecs[1, 0],
        'time': dt,
    }
