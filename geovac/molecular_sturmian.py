"""Molecular Sturmian betas via prolate spheroidal separation (v0.9.34).

Solves the two-center Coulomb problem in prolate spheroidal coordinates
(xi, eta) to find Sturmian eigenvalues beta for each molecular orbital.

The separated equations are:
  Angular: d/deta[(1-eta^2)dS/deta] + [A + c^2*eta^2 + b*eta - m^2/(1-eta^2)]S = 0
  Radial:  d/dxi[(xi^2-1)dR/dxi] + [-c^2*xi^2 + a*xi + A - m^2/(xi^2-1)]R = 0

Parameters (Sturmian with nuclear charge scaling beta):
  c^2 = p0^2*R^2/4    (fixed energy parameter)
  a = beta*(Z_A+Z_B)*R (nuclear attraction, beta-dependent)
  b = beta*(Z_B-Z_A)*R (heteronuclear splitting, beta-dependent)

Matching condition: A_angular(beta) + L_radial_top(beta) = 0.

v0.9.30: Added project_mo_betas_to_atom_centers() — dominant-overlap projection
of molecular orbital betas onto atom-centered hydrogenic orbitals.

v0.9.31: Added compute_h1_matrix() and compute_eri_matrix() for MO Sturmian FCI.
Numerical 2D quadrature on prolate spheroidal grid for one-electron integrals.
Population-weighted Slater F0 / Ohno-Klopman for two-electron integrals.
"""
from __future__ import annotations

import numpy as np
from math import factorial, sqrt as msqrt
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq
from scipy.special import assoc_laguerre, lpmv
from typing import Dict, List, Tuple


def _angular_sep_const(
    m_abs: int, n_sph: int, c: float, b: float = 0.0, n_basis: int = 50
) -> float:
    """Coulomb angular separation constant A_sep.

    Eigenvalue of H_eta = L_0 + c^2*eta^2 + b*eta from the TOP of
    the spectrum (the n_sph-th eigenvalue counting from the largest).

    For b=0 (homonuclear): equals -obl_cv(m, m+n_sph, c).
    """
    N = n_basis
    r_vals = np.arange(m_abs, m_abs + N, dtype=float)
    norms = np.array(
        [2.0 / (2*r+1) * factorial(int(r+m_abs)) / factorial(int(r-m_abs))
         for r in r_vals]
    )

    # nu matrix in orthonormal associated Legendre basis
    nu_mat = np.zeros((N, N))
    for i in range(N - 1):
        r = r_vals[i]
        val = (r - m_abs + 1) * np.sqrt(norms[i+1]) / ((2*r+1) * np.sqrt(norms[i]))
        nu_mat[i+1, i] = val
        nu_mat[i, i+1] = val

    H = np.diag(-r_vals * (r_vals + 1))
    H += c**2 * (nu_mat @ nu_mat)
    H += b * nu_mat

    evals = np.sort(np.linalg.eigvalsh(H))
    return evals[N - 1 - n_sph]


def _radial_top_evals(
    m_abs: int, c: float, a: float,
    n_grid: int = 1200, xi_max: float = 8.0, n_top: int = 5
) -> np.ndarray:
    """Top eigenvalues of the radial operator L_xi.

    L_xi = d/dxi[(xi^2-1)d/dxi] - c^2*xi^2 + a*xi - m^2/(xi^2-1)

    Uses self-adjoint FD with Neumann BC at xi=1 for m=0 (regular
    solution is nonzero at xi=1). Dirichlet R=0 at xi_max.

    Returns n_top largest eigenvalues in descending order.
    """
    N = n_grid
    xi_min = 1.0 + 0.0005
    h = (xi_max - xi_min) / (N + 1)
    xi = xi_min + (np.arange(N) + 1) * h

    xi2_1 = xi**2 - 1
    q = -c**2 * xi**2 + a * xi - m_abs**2 / xi2_1

    p_plus = (xi + h / 2)**2 - 1
    p_minus = (xi - h / 2)**2 - 1

    diag = -(p_plus + p_minus) / h**2 + q
    off = p_plus[:-1] / h**2

    if m_abs == 0:
        # Neumann BC at left: dR/dxi(xi_min) = 0
        # Ghost point: f_{-1} = f_0, so p_minus term vanishes
        diag[0] = -p_plus[0] / h**2 + q[0]

    evals = eigh_tridiagonal(diag, off, eigvals_only=True)
    return np.sort(evals)[::-1][:n_top]


def compute_molecular_sturmian_betas(
    Z_A: float, Z_B: float, R: float, p0: float, nmax: int,
    beta_min: float = 0.05, beta_max: float = 8.0, n_scan: int = 60,
    n_grid_radial: int = 1200
) -> List[Tuple[int, int, int, int, float]]:
    """Compute molecular Sturmian betas for all orbitals up to nmax.

    Each molecular orbital is labeled by (m, n_sph, n_rad) with
    principal quantum number n = 1 + m + n_sph + n_rad.

    Parameters
    ----------
    Z_A, Z_B : float
        Nuclear charges (Z_A >= Z_B by convention).
    R : float
        Internuclear distance in bohr.
    p0 : float
        Common momentum parameter (p0^2 = -2E).
    nmax : int
        Maximum principal quantum number.
    beta_min, beta_max : float
        Range for beta search.
    n_scan : int
        Number of scan points for root finding.
    n_grid_radial : int
        Grid points for radial FD solver.

    Returns
    -------
    results : list of (n, m, n_sph, n_rad, beta)
        Sorted by n, then m, then n_sph. Returns empty tuple for
        orbitals where no root was found.
    """
    c = p0 * R / 2
    Z_total = Z_A + Z_B
    Z_diff = Z_B - Z_A

    results: List[Tuple[int, int, int, int, float]] = []

    for n_princ in range(1, nmax + 1):
        for m in range(n_princ):
            for n_sph in range(n_princ - m):
                n_rad = n_princ - 1 - m - n_sph
                if n_rad < 0:
                    continue

                def mismatch(beta: float) -> float:
                    a = beta * Z_total * R
                    b = beta * Z_diff * R
                    A_ang = _angular_sep_const(m, n_sph, c, b)
                    L_ev = _radial_top_evals(
                        m, c, a, n_grid=n_grid_radial, n_top=n_rad + 3
                    )
                    if len(L_ev) <= n_rad:
                        return float('nan')
                    return A_ang + L_ev[n_rad]

                # Coarse scan
                betas = np.linspace(beta_min, beta_max, n_scan)
                vals = [mismatch(bt) for bt in betas]
                vals_arr = np.array(vals)

                # Find sign changes -> brentq
                beta_found = float('nan')
                for i in range(len(vals_arr) - 1):
                    v0, v1 = vals_arr[i], vals_arr[i+1]
                    if not np.isnan(v0) and not np.isnan(v1) and v0 * v1 < 0:
                        try:
                            root = brentq(mismatch, betas[i], betas[i+1],
                                          xtol=1e-8)
                            beta_found = root
                            break
                        except Exception:
                            pass

                results.append((n_princ, m, n_sph, n_rad, beta_found))

    return results


# ---------------------------------------------------------------------------
# Angular eigenvector (eta wavefunction) for MO-to-atom projection
# ---------------------------------------------------------------------------

def _angular_eigenvector(
    m_abs: int, n_sph: int, c: float, b: float = 0.0, n_basis: int = 50
) -> Tuple[float, np.ndarray, np.ndarray]:
    """Angular separation constant AND eigenvector.

    Returns the eigenvector in the orthonormal associated Legendre basis
    {P_r^m(eta) / norm_r}, r = m, m+1, ..., m+N-1.  The eigenvector
    coefficients allow reconstructing the angular wavefunction N(eta) on
    any eta grid.

    Returns
    -------
    A_sep : float
        Angular separation constant (same as _angular_sep_const).
    coeffs : ndarray, shape (n_basis,)
        Expansion coefficients in orthonormal P_r^m basis.
    r_vals : ndarray, shape (n_basis,)
        Legendre order values r = m, m+1, ...
    """
    N = n_basis
    r_vals = np.arange(m_abs, m_abs + N, dtype=float)
    norms = np.array(
        [2.0 / (2*r+1) * factorial(int(r+m_abs)) / factorial(int(r-m_abs))
         for r in r_vals]
    )

    # nu matrix in orthonormal associated Legendre basis
    nu_mat = np.zeros((N, N))
    for i in range(N - 1):
        r = r_vals[i]
        val = (r - m_abs + 1) * np.sqrt(norms[i+1]) / ((2*r+1) * np.sqrt(norms[i]))
        nu_mat[i+1, i] = val
        nu_mat[i, i+1] = val

    H = np.diag(-r_vals * (r_vals + 1))
    H += c**2 * (nu_mat @ nu_mat)
    H += b * nu_mat

    evals, evecs = np.linalg.eigh(H)
    # n_sph-th from top
    idx = N - 1 - n_sph
    return evals[idx], evecs[:, idx], r_vals


def _eval_angular_wavefunction(
    coeffs: np.ndarray, r_vals: np.ndarray, m_abs: int,
    eta_grid: np.ndarray
) -> np.ndarray:
    """Evaluate angular wavefunction N(eta) on a grid.

    N(eta) = sum_r c_r * P_r^m(eta) / sqrt(norm_r)

    where norm_r = 2/(2r+1) * (r+m)!/(r-m)!.

    Parameters
    ----------
    coeffs : ndarray
        Expansion coefficients from _angular_eigenvector.
    r_vals : ndarray
        Legendre orders.
    m_abs : int
        Azimuthal quantum number |m|.
    eta_grid : ndarray
        Grid of eta values in [-1, 1].

    Returns
    -------
    N_eta : ndarray
        Angular wavefunction values on the grid.
    """
    from scipy.special import lpmv

    N_eta = np.zeros_like(eta_grid)
    for i, r in enumerate(r_vals):
        ri = int(r)
        norm_r = 2.0 / (2*r+1) * factorial(ri + m_abs) / factorial(ri - m_abs)
        # P_r^m(eta) in scipy convention (includes Condon-Shortley phase)
        P_vals = lpmv(m_abs, ri, eta_grid)
        N_eta += coeffs[i] * P_vals / np.sqrt(norm_r)

    return N_eta


def _hydrogen_radial(n: int, l: int, Z: float, r: np.ndarray) -> np.ndarray:
    """Evaluate hydrogen-like radial wavefunction R_{nl}(r).

    R_{nl}(r) = N * (2Zr/n)^l * L_{n-l-1}^{2l+1}(2Zr/n) * exp(-Zr/n)
    """
    nr = n - l - 1
    N_coeff = (2.0 * Z / n) ** 1.5 * msqrt(
        factorial(nr) / (2.0 * n * factorial(n + l))
    )
    x = 2.0 * Z * r / n
    lag = assoc_laguerre(x, nr, 2 * l + 1)
    return N_coeff * x**l * lag * np.exp(-x / 2.0)


# ---------------------------------------------------------------------------
# MO-to-atom-center beta projection (v0.9.30)
# ---------------------------------------------------------------------------

def project_mo_betas_to_atom_centers(
    mo_results: List[Tuple[int, int, int, int, float]],
    Z_A: float, Z_B: float, R: float, p0: float, nmax: int,
    n_eta: int = 40, n_xi: int = 40, xi_max: float = 12.0,
) -> Tuple[Dict[Tuple[int, int, int], float], Dict[Tuple[int, int, int], float]]:
    """Project molecular orbital betas onto atom-centered hydrogenic orbitals.

    Uses dominant-overlap assignment: each atom-centered orbital gets the
    beta of the molecular orbital with which it has the largest squared
    overlap.

    The overlap is computed numerically on a 2D Gauss-Legendre grid in
    prolate spheroidal coordinates (xi, eta).

    Parameters
    ----------
    mo_results : list of (n, m, n_sph, n_rad, beta)
        Output from compute_molecular_sturmian_betas.
    Z_A, Z_B : float
        Nuclear charges (Z_A >= Z_B).
    R : float
        Internuclear distance in bohr.
    p0 : float
        Sturmian momentum parameter.
    nmax : int
        Maximum principal quantum number for atom-centered basis.
    n_eta, n_xi : int
        Number of Gauss-Legendre quadrature points in each dimension.
    xi_max : float
        Upper bound for xi integration.

    Returns
    -------
    beta_A : dict mapping (n, l, m) -> beta for A-center orbitals
    beta_B : dict mapping (n, l, m) -> beta for B-center orbitals
    """
    c = p0 * R / 2
    Z_total = Z_A + Z_B
    Z_diff = Z_B - Z_A

    # --- Build MO angular wavefunctions (eta-dependent) ---
    # Group MOs by |m| for block-diagonal projection
    mo_by_m: Dict[int, List[Tuple[int, int, int, float, np.ndarray, np.ndarray]]] = {}

    for n_princ, m_abs, n_sph, n_rad, beta in mo_results:
        if not np.isfinite(beta):
            continue
        b_param = beta * Z_diff * R
        _, coeffs, r_vals = _angular_eigenvector(m_abs, n_sph, c, b_param)
        mo_by_m.setdefault(m_abs, []).append(
            (n_princ, n_sph, n_rad, beta, coeffs, r_vals)
        )

    # --- Build atom-centered orbital list ---
    # GeoVac convention: states ordered (n, l, m) with m = -l, ..., +l
    atom_orbitals_A: List[Tuple[int, int, int]] = []
    atom_orbitals_B: List[Tuple[int, int, int]] = []
    for n in range(1, nmax + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                atom_orbitals_A.append((n, l, m))
                atom_orbitals_B.append((n, l, m))

    # --- Gauss-Legendre quadrature grids ---
    # eta in [-1, 1]
    eta_nodes, eta_weights = np.polynomial.legendre.leggauss(n_eta)
    # xi in [1, xi_max]: map from [-1, 1] to [1, xi_max]
    xi_01_nodes, xi_01_weights = np.polynomial.legendre.leggauss(n_xi)
    xi_nodes = (xi_max - 1.0) / 2.0 * xi_01_nodes + (xi_max + 1.0) / 2.0
    xi_weights = xi_01_weights * (xi_max - 1.0) / 2.0

    # --- Prolate spheroidal -> Cartesian coordinate maps ---
    # r_A = R*(xi + eta)/2, r_B = R*(xi - eta)/2
    # Volume element: (R/2)^3 * (xi^2 - eta^2) dxi deta dphi
    # For m=0 sigma orbitals: phi integral = 2*pi
    # For |m|>0: phi integral of |e^{im*phi}|^2 = 2*pi

    # Precompute 2D grids
    XI, ETA = np.meshgrid(xi_nodes, eta_nodes)  # shape (n_eta, n_xi)
    r_A = R * (XI + ETA) / 2.0  # distance from center A
    r_B = R * (XI - ETA) / 2.0  # distance from center B
    vol_elem = (R / 2.0)**3 * (XI**2 - ETA**2)  # no dphi factor yet

    # Weight matrix for 2D quadrature (outer product of weights)
    W2D = np.outer(eta_weights, xi_weights)  # shape (n_eta, n_xi)

    # cos(theta_A) and cos(theta_B) for spherical harmonics
    # For center A: cos(theta_A) = (xi*eta + 1) / (xi + eta)
    # For center B: cos(theta_B) = (xi*eta - 1) / (xi - eta)
    # Handle singularity at xi + eta = 0 or xi - eta = 0
    with np.errstate(divide='ignore', invalid='ignore'):
        cos_theta_A = np.where(
            np.abs(XI + ETA) > 1e-12,
            (XI * ETA + 1.0) / (XI + ETA),
            1.0
        )
        cos_theta_B = np.where(
            np.abs(XI - ETA) > 1e-12,
            (XI * ETA - 1.0) / (XI - ETA),
            1.0
        )

    # Clip to [-1, 1] for numerical safety
    cos_theta_A = np.clip(cos_theta_A, -1.0, 1.0)
    cos_theta_B = np.clip(cos_theta_B, -1.0, 1.0)

    # Ensure r_A, r_B non-negative
    r_A = np.maximum(r_A, 0.0)
    r_B = np.maximum(r_B, 0.0)

    # --- MO radial wavefunction evaluator (FD eigenvector on xi grid) ---
    def _mo_radial_on_xi(m_abs: int, n_rad: int, beta: float) -> np.ndarray:
        """Evaluate MO radial wavefunction M(xi) on xi_nodes via FD."""
        a_param = beta * Z_total * R
        n_grid_fd = 600
        xi_min_fd = 1.0 + 0.0005
        h = (xi_max - xi_min_fd) / (n_grid_fd + 1)
        xi_fd = xi_min_fd + (np.arange(n_grid_fd) + 1) * h

        xi2_1 = xi_fd**2 - 1
        q = -c**2 * xi_fd**2 + a_param * xi_fd - m_abs**2 / xi2_1

        p_plus = (xi_fd + h / 2)**2 - 1
        p_minus = (xi_fd - h / 2)**2 - 1

        diag_fd = -(p_plus + p_minus) / h**2 + q
        off_fd = p_plus[:-1] / h**2

        if m_abs == 0:
            diag_fd[0] = -p_plus[0] / h**2 + q[0]

        evals_fd, evecs_fd = eigh_tridiagonal(diag_fd, off_fd)
        # n_rad-th from top
        idx_fd = n_grid_fd - 1 - n_rad
        evec = evecs_fd[:, idx_fd]

        # Normalize on FD grid
        norm = np.sqrt(np.sum(evec**2) * h)
        if norm > 1e-30:
            evec /= norm

        # Interpolate to quadrature xi_nodes
        return np.interp(xi_nodes, xi_fd, evec, left=0.0, right=0.0)

    # --- Compute overlaps and assign betas ---
    beta_A: Dict[Tuple[int, int, int], float] = {}
    beta_B: Dict[Tuple[int, int, int], float] = {}

    for n_at, l_at, m_at in atom_orbitals_A:
        m_abs_at = abs(m_at)

        # Only MOs with same |m| can overlap
        if m_abs_at not in mo_by_m:
            # No matching MO — use atomic beta = 1/n (fallback)
            beta_A[(n_at, l_at, m_at)] = 1.0 / n_at
            beta_B[(n_at, l_at, m_at)] = 1.0 / n_at
            continue

        mos = mo_by_m[m_abs_at]

        # Evaluate atom-centered radial * angular on 2D grid
        # Center A: R_{nl}(r_A) * Y_l^m(theta_A, phi=0)
        # Y_l^m(theta, 0) = sqrt((2l+1)/(4pi) * (l-m)!/(l+m)!) * P_l^m(cos theta)
        # For overlap, the normalization constant cancels in the ratio;
        # we include it for proper magnitude.
        R_at_A = _hydrogen_radial(n_at, l_at, Z_A, r_A)
        Y_norm = np.sqrt(
            (2*l_at + 1) / (4.0 * np.pi) *
            factorial(l_at - m_abs_at) / factorial(l_at + m_abs_at)
        )
        P_at_A = lpmv(m_abs_at, l_at, cos_theta_A)
        chi_A = R_at_A * Y_norm * P_at_A

        R_at_B = _hydrogen_radial(n_at, l_at, Z_B, r_B)
        P_at_B = lpmv(m_abs_at, l_at, cos_theta_B)
        chi_B = R_at_B * Y_norm * P_at_B

        best_overlap_A = -1.0
        best_beta_A = 1.0 / n_at  # fallback
        best_overlap_B = -1.0
        best_beta_B = 1.0 / n_at

        for n_mo, n_sph, n_rad, beta, ang_coeffs, ang_r_vals in mos:
            # MO angular wavefunction on eta grid (broadcast to 2D)
            N_eta_1d = _eval_angular_wavefunction(
                ang_coeffs, ang_r_vals, m_abs_at, eta_nodes
            )
            # MO radial wavefunction on xi grid
            M_xi_1d = _mo_radial_on_xi(m_abs_at, n_rad, beta)

            # 2D MO wavefunction: M(xi) * N(eta)
            MO_2d = np.outer(N_eta_1d, M_xi_1d)  # shape (n_eta, n_xi)

            # Overlap with center A orbital: integral of chi_A * MO * vol * W
            # phi integral gives 2*pi for all m (normalization)
            overlap_A = np.sum(chi_A * MO_2d * vol_elem * W2D) * 2.0 * np.pi
            sq_A = overlap_A**2

            if sq_A > best_overlap_A:
                best_overlap_A = sq_A
                best_beta_A = beta

            # Overlap with center B orbital
            overlap_B = np.sum(chi_B * MO_2d * vol_elem * W2D) * 2.0 * np.pi
            sq_B = overlap_B**2

            if sq_B > best_overlap_B:
                best_overlap_B = sq_B
                best_beta_B = beta

        beta_A[(n_at, l_at, m_at)] = best_beta_A
        beta_B[(n_at, l_at, m_at)] = best_beta_B

    return beta_A, beta_B


# ---------------------------------------------------------------------------
# MO wavefunction evaluation on prolate spheroidal grid (v0.9.31)
# ---------------------------------------------------------------------------

def _eval_mo_wavefunction(
    mo_tuple: Tuple[int, int, int, int, float],
    xi_grid: np.ndarray,
    eta_grid: np.ndarray,
    R: float,
    p0: float,
    Z_A: float,
    Z_B: float,
    n_grid_fd: int = 600,
    xi_max_fd: float = 8.0,
) -> np.ndarray:
    """Evaluate MO Sturmian wavefunction on a 2D prolate spheroidal grid.

    Parameters
    ----------
    mo_tuple : (n, m, n_sph, n_rad, beta)
        MO quantum numbers and Sturmian eigenvalue.
    xi_grid : ndarray, shape (n_xi,)
        Radial coordinate values (xi >= 1).
    eta_grid : ndarray, shape (n_eta,)
        Angular coordinate values (eta in [-1, 1]).
    R : float
        Internuclear distance in bohr.
    p0 : float
        Sturmian momentum parameter.
    Z_A, Z_B : float
        Nuclear charges.
    n_grid_fd : int
        Grid points for FD radial solver.
    xi_max_fd : float
        Upper bound for FD radial grid.

    Returns
    -------
    psi : ndarray, shape (n_eta, n_xi)
        MO wavefunction values. Does NOT include phi-dependent factor.
    """
    n_princ, m_abs, n_sph, n_rad, beta = mo_tuple
    c = p0 * R / 2.0
    Z_total = Z_A + Z_B
    Z_diff = Z_B - Z_A

    # --- Angular part: N(eta) ---
    b_param = beta * Z_diff * R
    _, ang_coeffs, ang_r_vals = _angular_eigenvector(m_abs, n_sph, c, b_param)
    N_eta = _eval_angular_wavefunction(ang_coeffs, ang_r_vals, m_abs, eta_grid)

    # --- Radial part: M(xi) via FD eigenvector ---
    a_param = beta * Z_total * R
    xi_min_fd = 1.0 + 0.0005
    h = (xi_max_fd - xi_min_fd) / (n_grid_fd + 1)
    xi_fd = xi_min_fd + (np.arange(n_grid_fd) + 1) * h

    xi2_1 = xi_fd**2 - 1
    q = -c**2 * xi_fd**2 + a_param * xi_fd - m_abs**2 / xi2_1

    p_plus = (xi_fd + h / 2)**2 - 1
    p_minus = (xi_fd - h / 2)**2 - 1

    diag_fd = -(p_plus + p_minus) / h**2 + q
    off_fd = p_plus[:-1] / h**2

    if m_abs == 0:
        diag_fd[0] = -p_plus[0] / h**2 + q[0]

    evals_fd, evecs_fd = eigh_tridiagonal(diag_fd, off_fd)
    idx_fd = n_grid_fd - 1 - n_rad
    evec = evecs_fd[:, idx_fd]

    # Normalize on FD grid: integral of M(xi)^2 * (xi^2 - 1) dxi
    # Use weight (xi^2 - 1) for proper Sturm-Liouville normalization
    weight_xi = xi_fd**2 - 1
    norm = np.sqrt(np.sum(evec**2 * weight_xi) * h)
    if norm > 1e-30:
        evec /= norm

    # Interpolate to quadrature grid
    M_xi = np.interp(xi_grid, xi_fd, evec, left=0.0, right=0.0)

    # Also normalize angular part: integral of N(eta)^2 * (1 - eta^2)^... deta
    # The eigenvector from _angular_eigenvector is already in orthonormal basis

    # 2D product: psi(eta, xi) = N(eta) * M(xi)
    return np.outer(N_eta, M_xi)


# ---------------------------------------------------------------------------
# One-electron Hamiltonian matrix in MO Sturmian basis (v0.9.31)
# ---------------------------------------------------------------------------

def compute_h1_matrix(
    mo_results: List[Tuple[int, int, int, int, float]],
    Z_A: float,
    Z_B: float,
    R: float,
    p0: float,
    n_grid: int = 40,
    xi_max: float = 12.0,
    orth_method: str = 'canonical',
    orth_threshold: float = 1e-4,
) -> Tuple[np.ndarray, np.ndarray, int, np.ndarray, np.ndarray]:
    """Compute one-electron Hamiltonian in orthogonalized MO basis.

    Sturmian orbitals are orthogonal under the potential-weighted inner product
    <S_i|-V|S_j>, NOT the standard inner product <S_i|S_j>. The Slater-Condon
    rules require an orthonormal basis. This function:
      1. Computes raw H1 via the Sturmian identity
      2. Computes the overlap matrix S
      3. Orthogonalizes: canonical (default) or Löwdin

    Sturmian identity:
      H_ij = -p0^2/2 O_ij + (1 - beta_j) V_ij  (symmetrized)

    Parameters
    ----------
    mo_results : list of (n, m, n_sph, n_rad, beta)
        MO orbitals from compute_molecular_sturmian_betas.
    Z_A, Z_B : float
        Nuclear charges.
    R : float
        Internuclear distance in bohr.
    p0 : float
        Sturmian momentum parameter.
    n_grid : int
        Gauss-Legendre quadrature points per dimension.
    xi_max : float
        Upper radial integration bound.
    orth_method : str
        'canonical' (default) or 'lowdin'.
    orth_threshold : float
        Eigenvalue threshold for canonical orthogonalization.

    Returns
    -------
    H1_orth : ndarray, shape (n_orth, n_orth)
        Orthogonalized one-electron Hamiltonian in Hartree.
    X : ndarray, shape (n_mo, n_orth)
        Transformation matrix (n_orth <= n_mo after dropping near-LD vectors).
    n_dropped : int
        Number of linearly dependent basis functions removed.
    S : ndarray, shape (n_mo, n_mo)
        Raw overlap matrix (for diagnostics).
    H1_raw : ndarray, shape (n_mo, n_mo)
        Raw (non-orthogonalized) H1 matrix (for diagnostics).
    """
    # Filter valid MOs
    valid_mos = [(i, mo) for i, mo in enumerate(mo_results)
                 if np.isfinite(mo[4])]
    n_mo = len(valid_mos)
    H1 = np.zeros((n_mo, n_mo))

    if n_mo == 0:
        return H1, np.eye(n_mo)

    # --- Gauss-Legendre quadrature grids ---
    eta_nodes, eta_weights = np.polynomial.legendre.leggauss(n_grid)
    xi_01_nodes, xi_01_weights = np.polynomial.legendre.leggauss(n_grid)
    xi_nodes = (xi_max - 1.0) / 2.0 * xi_01_nodes + (xi_max + 1.0) / 2.0
    xi_weights = xi_01_weights * (xi_max - 1.0) / 2.0

    # 2D grid
    XI, ETA = np.meshgrid(xi_nodes, eta_nodes)  # (n_eta, n_xi)
    W2D = np.outer(eta_weights, xi_weights)

    # Prolate spheroidal volume element: (R/2)^3 (xi^2 - eta^2)
    vol_elem = (R / 2.0)**3 * (XI**2 - ETA**2)

    # Nuclear distances
    r_A = R * (XI + ETA) / 2.0
    r_B = R * (XI - ETA) / 2.0
    r_A = np.maximum(r_A, 1e-14)
    r_B = np.maximum(r_B, 1e-14)

    # Nuclear potential on grid: V = -Z_A/r_A - Z_B/r_B
    V_nuc = -Z_A / r_A - Z_B / r_B

    # --- Evaluate and normalize all MO wavefunctions on grid ---
    psi_mos = []
    mo_m_vals = []
    for idx_i, mo in valid_mos:
        psi = _eval_mo_wavefunction(mo, xi_nodes, eta_nodes, R, p0, Z_A, Z_B)
        # Normalize: <S_k|S_k> = integral |psi|^2 vol dV = 1
        phi_factor = 2.0 * np.pi
        norm_sq = np.sum(psi**2 * vol_elem * W2D) * phi_factor
        if norm_sq > 1e-30:
            psi /= np.sqrt(norm_sq)
        psi_mos.append(psi)
        mo_m_vals.append(mo[1])  # |m| value

    # --- Compute overlap matrix O_ij and potential matrix V_ij ---
    n_valid = len(valid_mos)
    O = np.zeros((n_valid, n_valid))
    V = np.zeros((n_valid, n_valid))

    phi_factor = 2.0 * np.pi
    for i in range(n_valid):
        m_i = mo_m_vals[i]
        psi_i = psi_mos[i]
        for j in range(i, n_valid):
            m_j = mo_m_vals[j]
            # phi orthogonality: <e^{im_i phi} | e^{im_j phi}> = 2pi delta_{m_i, m_j}
            if m_i != m_j:
                continue

            psi_j = psi_mos[j]

            # Overlap: <S_i|S_j>
            o_ij = np.sum(psi_i * psi_j * vol_elem * W2D) * phi_factor
            O[i, j] = o_ij
            O[j, i] = o_ij

            # Potential: <S_i|V_mol|S_j>
            v_ij = np.sum(psi_i * V_nuc * psi_j * vol_elem * W2D) * phi_factor
            V[i, j] = v_ij
            V[j, i] = v_ij

    # --- Assemble raw H1 using Sturmian identity ---
    # H_ij = -p0^2/2 * O_ij + (1 - beta_j) * V_ij
    # Symmetrize using average of column-j and column-i formulas:
    # H_ij^sym = -p0^2/2 O_ij + [(1-beta_i)+(1-beta_j)]/2 * V_ij
    betas = [valid_mos[k][1][4] for k in range(n_valid)]

    H1_raw = np.zeros((n_valid, n_valid))
    for i in range(n_valid):
        for j in range(n_valid):
            avg_factor = ((1.0 - betas[i]) + (1.0 - betas[j])) / 2.0
            H1_raw[i, j] = -p0**2 / 2.0 * O[i, j] + avg_factor * V[i, j]

    # --- Orthogonalization ---
    if orth_method == 'none_raw':
        # Return raw matrices without orthogonalization (for dual-p0 mode)
        return H1_raw, np.eye(n_valid), 0, O, H1_raw

    evals_S, evecs_S = np.linalg.eigh(O)

    if orth_method == 'canonical':
        # Canonical orthogonalization: drop near-LD eigenvectors
        keep = evals_S > orth_threshold
        n_dropped = int(np.sum(~keep))
        U_keep = evecs_S[:, keep]
        s_keep = evals_S[keep]
        X = U_keep @ np.diag(s_keep**(-0.5))  # shape (n_mo, n_orth)
        H1_orth = X.T @ H1_raw @ X
    else:
        # Löwdin: S^{-1/2} H1 S^{-1/2}
        n_dropped = 0
        evals_S = np.maximum(evals_S, 1e-10)
        X = evecs_S @ np.diag(evals_S**(-0.5)) @ evecs_S.T
        H1_orth = X @ H1_raw @ X

    return H1_orth, X, n_dropped, O, H1_raw


# ---------------------------------------------------------------------------
# Two-electron integrals in MO Sturmian basis (v0.9.31)
# ---------------------------------------------------------------------------

def compute_eri_matrix(
    mo_results: List[Tuple[int, int, int, int, float]],
    Z_A: float,
    Z_B: float,
    R: float,
    p0: float,
    nmax: int,
    U: np.ndarray = None,
    n_grid: int = 40,
    xi_max: float = 12.0,
) -> np.ndarray:
    """Compute two-electron integrals in MO Sturmian basis (approximate).

    Uses population-weighted decomposition:
    <S_k S_l | 1/r12 | S_m S_n> ≈ w_k^A w_l^A F0_AA + w_k^B w_l^B F0_BB
                                    + (w_k^A w_l^B + w_k^B w_l^A) J_AB

    If U (Löwdin S^{-1/2}) is provided, the ERI is transformed to the
    orthonormal basis: ERI'[p,q,r,s] = sum_{abcd} U[a,p]U[b,q]U[c,r]U[d,s] ERI[a,b,c,d]

    Parameters
    ----------
    mo_results : list of (n, m, n_sph, n_rad, beta)
        MO orbitals from compute_molecular_sturmian_betas.
    Z_A, Z_B : float
        Nuclear charges.
    R : float
        Internuclear distance in bohr.
    p0 : float
        Sturmian momentum parameter.
    nmax : int
        Maximum principal quantum number.
    n_grid : int
        Quadrature points for population weight computation.
    xi_max : float
        Upper radial integration bound.

    Returns
    -------
    ERI : ndarray, shape (n_mo, n_mo, n_mo, n_mo)
        Two-electron integral tensor (Coulomb-only approximation).
    """
    valid_mos = [(i, mo) for i, mo in enumerate(mo_results)
                 if np.isfinite(mo[4])]
    n_mo = len(valid_mos)
    ERI = np.zeros((n_mo, n_mo, n_mo, n_mo))

    if n_mo == 0:
        return ERI

    # --- Compute population weights for each MO ---
    weights_A, weights_B = _compute_population_weights(
        valid_mos, Z_A, Z_B, R, p0, nmax, n_grid, xi_max
    )

    # --- Compute same-center Slater F0 integrals (fast analytical table) ---
    # Use _slater_f0_fast() which avoids expensive numerical integration.
    mo_n_vals = [mo[0] for _, mo in valid_mos]

    F0_A = np.zeros((n_mo, n_mo))
    F0_B = np.zeros((n_mo, n_mo))

    for i in range(n_mo):
        ni = mo_n_vals[i]
        for j in range(i, n_mo):
            nj = mo_n_vals[j]
            f0_a = _slater_f0_fast(ni, nj, Z_A)
            f0_b = _slater_f0_fast(ni, nj, Z_B)
            F0_A[i, j] = F0_A[j, i] = f0_a
            F0_B[i, j] = F0_B[j, i] = f0_b

    # Cross-center Ohno-Klopman: J_AB = 1/sqrt(R^2 + (r_A + r_B)^2/4)
    J_AB = np.zeros((n_mo, n_mo))
    for i in range(n_mo):
        ni = mo_n_vals[i]
        r_i = ni**2 / max(Z_A, 0.5)
        for j in range(i, n_mo):
            nj = mo_n_vals[j]
            r_j = nj**2 / max(Z_B, 0.5)
            j_ab = 1.0 / np.sqrt(R**2 + ((r_i + r_j) / 2.0)**2)
            J_AB[i, j] = J_AB[j, i] = j_ab

    # --- Assemble ERI tensor ---
    # Only Coulomb-type integrals: <kl|kl>
    for k in range(n_mo):
        wkA, wkB = weights_A[k], weights_B[k]
        for l in range(n_mo):
            wlA, wlB = weights_A[l], weights_B[l]

            j_kl = (wkA * wlA * F0_A[k, l]
                    + wkB * wlB * F0_B[k, l]
                    + (wkA * wlB + wkB * wlA) * J_AB[k, l])

            ERI[k, l, k, l] = j_kl
            ERI[k, l, l, k] = j_kl

    # --- 4-index transform if U provided ---
    # U can be (n_mo, n_mo) for Löwdin or (n_mo, n_orth) for canonical
    if U is not None:
        n_orth = U.shape[1]
        # ERI'[p,q,r,s] = Σ_{abcd} U[a,p] U[b,q] U[c,r] U[d,s] ERI[a,b,c,d]
        tmp = np.einsum('ai,abcd->ibcd', U, ERI)
        tmp = np.einsum('bj,ibcd->ijcd', U, tmp)
        tmp = np.einsum('ck,ijcd->ijkd', U, tmp)
        ERI_orth = np.einsum('dl,ijkd->ijkl', U, tmp)
        return ERI_orth

    return ERI


def _slater_f0_fast(na: int, nb: int, Z: float) -> float:
    """Fast analytical Slater F0 direct Coulomb integral for s-orbitals.

    F0(ns, ns'; Z) = Z * f(na, nb) where f is a pure number.

    Known analytical values (from hydrogen radial integrals):
      F0(1s,1s) = 5Z/8
      F0(1s,2s) = 17Z/81
      F0(1s,3s) = 83Z/648
      F0(2s,2s) = 77Z/512
      F0(2s,3s) = 1793Z/19683  (approx 0.09109Z)
      F0(3s,3s) = 7129Z/131072 (approx 0.05440Z)

    For cases not in the table, uses the 1/n-scaling approximation:
      F0(ns, n's) ≈ Z * 5/(8*n*n')
    """
    # Ensure na <= nb for lookup
    n1, n2 = min(na, nb), max(na, nb)

    # Analytical table: f0_table[(na, nb)] = F0/Z
    f0_table = {
        (1, 1): 5.0 / 8.0,            # 0.625
        (1, 2): 17.0 / 81.0,          # 0.20988
        (1, 3): 83.0 / 648.0,         # 0.12809
        (2, 2): 77.0 / 512.0,         # 0.15039
        (2, 3): 1793.0 / 19683.0,     # 0.09109
        (3, 3): 7129.0 / 131072.0,    # 0.05440
    }

    f0_ratio = f0_table.get((n1, n2))
    if f0_ratio is not None:
        return Z * f0_ratio

    # Fallback: 1/n scaling
    return Z * 5.0 / (8.0 * n1 * n2)


def _compute_population_weights(
    valid_mos: List[Tuple[int, Tuple[int, int, int, int, float]]],
    Z_A: float,
    Z_B: float,
    R: float,
    p0: float,
    nmax: int,
    n_grid: int = 40,
    xi_max: float = 12.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute A-center and B-center population weights for each MO.

    w_k^A = integral of |S_k|^2 over A-half-space / total integral
    where A-half-space is eta > 0 (closer to nucleus A).

    Parameters
    ----------
    valid_mos : list of (index, (n, m, n_sph, n_rad, beta))
        Valid MO list.
    Z_A, Z_B : float
        Nuclear charges.
    R : float
        Internuclear distance.
    p0 : float
        Sturmian momentum.
    nmax : int
        Max principal quantum number.
    n_grid : int
        Quadrature points.
    xi_max : float
        Upper radial bound.

    Returns
    -------
    weights_A : ndarray, shape (n_mo,)
        A-center population weight for each MO.
    weights_B : ndarray, shape (n_mo,)
        B-center population weight (= 1 - weights_A).
    """
    n_mo = len(valid_mos)
    weights_A = np.zeros(n_mo)
    weights_B = np.zeros(n_mo)

    # Quadrature grids
    eta_nodes, eta_weights = np.polynomial.legendre.leggauss(n_grid)
    xi_01_nodes, xi_01_weights = np.polynomial.legendre.leggauss(n_grid)
    xi_nodes = (xi_max - 1.0) / 2.0 * xi_01_nodes + (xi_max + 1.0) / 2.0
    xi_weights = xi_01_weights * (xi_max - 1.0) / 2.0

    XI, ETA = np.meshgrid(xi_nodes, eta_nodes)
    W2D = np.outer(eta_weights, xi_weights)
    vol_elem = (R / 2.0)**3 * (XI**2 - ETA**2)

    for k, (idx_k, mo_k) in enumerate(valid_mos):
        psi = _eval_mo_wavefunction(mo_k, xi_nodes, eta_nodes, R, p0, Z_A, Z_B)
        density = psi**2 * vol_elem * W2D * 2.0 * np.pi  # includes phi factor

        total = np.sum(density)
        if total < 1e-30:
            weights_A[k] = 0.5
            weights_B[k] = 0.5
            continue

        # A-center: eta > 0 (nucleus A at eta=+1 in prolate spheroidal coords)
        mask_A = ETA > 0
        pop_A = np.sum(density[mask_A])

        weights_A[k] = pop_A / total
        weights_B[k] = 1.0 - weights_A[k]

    return weights_A, weights_B


# ---------------------------------------------------------------------------
# Exact direct Coulomb integrals via Poisson solve (v0.9.32)
# ---------------------------------------------------------------------------

def compute_exact_j_integrals(
    mo_results: List[Tuple[int, int, int, int, float]],
    Z_A: float,
    Z_B: float,
    R: float,
    p0: float,
    n_grid: int = 30,
    xi_max: float = 12.0,
    l_max_poisson: int = 8,
) -> np.ndarray:
    """Compute direct Coulomb integrals J_kl = <S_k S_l | 1/r12 | S_k S_l>.

    Uses Poisson solve via spherical harmonic expansion on prolate spheroidal
    grid. For each MO density rho_k(r) = |S_k(r)|^2:
      1. Expand rho_k in spherical harmonics about the bond midpoint
      2. Compute Hartree potential V_H[rho_k] via 1/r^{l+1} and r^l expansion
      3. J_kl = integral rho_l(r) V_H[rho_k](r) dV

    Parameters
    ----------
    mo_results : list of (n, m, n_sph, n_rad, beta)
        MO orbitals from compute_molecular_sturmian_betas.
    Z_A, Z_B : float
        Nuclear charges.
    R : float
        Internuclear distance in bohr.
    p0 : float
        Sturmian momentum parameter.
    n_grid : int
        Gauss-Legendre quadrature points per dimension.
    xi_max : float
        Upper radial integration bound.
    l_max_poisson : int
        Maximum l for spherical harmonic expansion of density.

    Returns
    -------
    J : ndarray, shape (n_mo, n_mo)
        Direct Coulomb integral matrix (symmetric).
    """
    valid_mos = [(i, mo) for i, mo in enumerate(mo_results)
                 if np.isfinite(mo[4])]
    n_mo = len(valid_mos)
    J = np.zeros((n_mo, n_mo))

    if n_mo == 0:
        return J

    # --- Build 3D quadrature grid in prolate spheroidal coords ---
    # xi in [1, xi_max], eta in [-1, 1], phi in [0, 2*pi]
    eta_nodes, eta_weights = np.polynomial.legendre.leggauss(n_grid)
    xi_01_nodes, xi_01_weights = np.polynomial.legendre.leggauss(n_grid)
    xi_nodes = (xi_max - 1.0) / 2.0 * xi_01_nodes + (xi_max + 1.0) / 2.0
    xi_weights = xi_01_weights * (xi_max - 1.0) / 2.0

    # For m=0 MOs (sigma): density is phi-independent, phi integral = 2*pi
    # For |m|>0: density |e^{im*phi}|^2 = 1, same result

    XI, ETA = np.meshgrid(xi_nodes, eta_nodes)  # (n_eta, n_xi)
    W2D = np.outer(eta_weights, xi_weights)
    vol_elem = (R / 2.0)**3 * (XI**2 - ETA**2)

    # Convert to spherical coords centered at bond midpoint
    # Prolate spheroidal -> Cartesian: z = R/2 * xi * eta, rho = R/2 * sqrt((xi^2-1)(1-eta^2))
    z_cart = R / 2.0 * XI * ETA
    rho_cart = R / 2.0 * np.sqrt(np.maximum((XI**2 - 1) * (1 - ETA**2), 0.0))
    r_sph = np.sqrt(z_cart**2 + rho_cart**2)
    r_sph = np.maximum(r_sph, 1e-14)
    cos_theta = z_cart / r_sph

    # --- Evaluate all MO densities on grid ---
    densities = []
    for k, (idx_k, mo_k) in enumerate(valid_mos):
        psi = _eval_mo_wavefunction(mo_k, xi_nodes, eta_nodes, R, p0, Z_A, Z_B)
        # Normalize
        norm_sq = np.sum(psi**2 * vol_elem * W2D) * 2.0 * np.pi
        if norm_sq > 1e-30:
            psi /= np.sqrt(norm_sq)
        densities.append(psi**2)  # |psi|^2

    # --- Compute Hartree potential for each density via multipole expansion ---
    # V_H[rho](r) = sum_l (4*pi/(2l+1)) * [r^l * B_l(r) / r + A_l(r) / r^{l+1}] * P_l(cos_theta)
    # where A_l(r) = integral_{r'<r} rho_l(r') r'^{l+2} dr'
    #       B_l(r) = integral_{r'>r} rho_l(r') / r'^{l-1} dr'
    # For azimuthally symmetric density (m=0 projection), only m=0 harmonics contribute.

    # Precompute Legendre polynomials P_l(cos_theta) on grid
    from scipy.special import eval_legendre
    P_l_grid = []
    for l in range(l_max_poisson + 1):
        P_l_grid.append(eval_legendre(l, cos_theta))

    # For each density, compute its Hartree potential on the grid
    def _hartree_potential(density_2d: np.ndarray) -> np.ndarray:
        """Compute Hartree potential from 2D density via multipole expansion.

        For each l, compute the radial integrals by sorting grid points by r.
        """
        n_eta, n_xi = density_2d.shape
        # Flatten for radial sorting
        r_flat = r_sph.ravel()
        cos_flat = cos_theta.ravel()
        rho_flat = density_2d.ravel()
        vol_flat = (vol_elem * W2D).ravel() * 2.0 * np.pi  # include phi factor

        # Sort by radius for cumulative integration
        sort_idx = np.argsort(r_flat)
        r_sorted = r_flat[sort_idx]
        rho_sorted = rho_flat[sort_idx]
        vol_sorted = vol_flat[sort_idx]

        n_pts = len(r_flat)
        V_H_flat = np.zeros(n_pts)

        for l in range(l_max_poisson + 1):
            # Project density onto P_l: rho_l(r_i) = rho(r_i) * P_l(cos_theta_i)
            P_l_flat = P_l_grid[l].ravel()
            rho_l_sorted = rho_sorted * P_l_flat[sort_idx]
            charge_sorted = rho_l_sorted * vol_sorted

            # Cumulative inner integral: A_l(r_i) = sum_{j: r_j < r_i} charge_j * r_j^l
            # Cumulative outer integral: B_l(r_i) = sum_{j: r_j > r_i} charge_j / r_j^{l+1}
            r_pow_l = r_sorted**l
            r_pow_lp1_inv = np.where(r_sorted > 1e-14, r_sorted**(-(l + 1)), 0.0)

            inner_terms = charge_sorted * r_pow_l
            outer_terms = charge_sorted * r_pow_lp1_inv

            A_l = np.cumsum(inner_terms)  # cumulative sum up to i
            B_l_rev = np.cumsum(outer_terms[::-1])[::-1]  # cumulative sum from i to end

            # V_H contribution from this l
            # V_l(r_i) = (4*pi/(2l+1)) * [A_l(r_i)/r_i^{l+1} + B_l(r_i)*r_i^l]
            # But A_l includes point i, need A_l(r_i) = sum_{j<i}, so shift
            A_l_shifted = np.zeros_like(A_l)
            A_l_shifted[1:] = A_l[:-1]
            B_l_shifted = np.zeros_like(B_l_rev)
            B_l_shifted[:-1] = B_l_rev[1:]

            r_inv_lp1 = np.where(r_sorted > 1e-14, r_sorted**(-(l + 1)), 0.0)

            # No extra prefactor: 1/|r-r'| = sum_l P_l(cos gamma) r_<^l/r_>^{l+1}
            # The 2*pi phi factor is already in the charge (vol includes dphi).
            V_l_sorted = A_l_shifted * r_inv_lp1 + B_l_shifted * r_pow_l

            # P_l on sorted points
            P_l_sorted = P_l_flat[sort_idx]
            V_H_sorted = V_l_sorted * P_l_sorted

            # Unsort back to original order
            unsort = np.argsort(sort_idx)
            V_H_flat += V_H_sorted[unsort]

        return V_H_flat.reshape(n_eta, n_xi)

    # --- Compute J_kl = integral rho_l * V_H[rho_k] dV ---
    hartree_potentials = []
    for k in range(n_mo):
        V_Hk = _hartree_potential(densities[k])
        hartree_potentials.append(V_Hk)

    for k in range(n_mo):
        for l in range(k, n_mo):
            # J_kl = integral rho_l(r) V_H[rho_k](r) dV
            j_kl = np.sum(densities[l] * hartree_potentials[k] * vol_elem * W2D) * 2.0 * np.pi
            J[k, l] = j_kl
            J[l, k] = j_kl

    return J


# ---------------------------------------------------------------------------
# Exact exchange integrals via Poisson solve (v0.9.33)
# ---------------------------------------------------------------------------

def compute_exact_k_integrals(
    mo_results: List[Tuple[int, int, int, int, float]],
    Z_A: float,
    Z_B: float,
    R: float,
    p0: float,
    n_grid: int = 30,
    xi_max: float = 12.0,
    l_max_poisson: int = 8,
) -> np.ndarray:
    """Compute exchange integrals K_kl = <S_k S_l|1/r12|S_l S_k>.

    Uses Poisson solve via multipole expansion on prolate spheroidal grid.
    Exchange density: rho_kl(r) = S_k(r) * S_l(r)
    Poisson solve: nabla^2 V_kl = -4*pi*rho_kl
    K_kl = integral rho_kl(r) * V_kl(r) dV

    m-selection rule: K_kl is nonzero only when m_k = m_l.
    This follows from the phi integral: both phi1 and phi2 integrals
    require m' = m_l - m_k AND m' = m_k - m_l, forcing m_k = m_l.

    K_kl = K_lk (symmetric). K_kk = J_kk by definition.

    Parameters
    ----------
    mo_results : list of (n, m, n_sph, n_rad, beta)
        MO orbitals from compute_molecular_sturmian_betas.
    Z_A, Z_B : float
        Nuclear charges.
    R : float
        Internuclear distance in bohr.
    p0 : float
        Sturmian momentum parameter.
    n_grid : int
        Gauss-Legendre quadrature points per dimension.
    xi_max : float
        Upper radial integration bound.
    l_max_poisson : int
        Maximum l for spherical harmonic expansion.

    Returns
    -------
    K : ndarray, shape (n_mo, n_mo)
        Exchange integral matrix (symmetric).
    """
    valid_mos = [(i, mo) for i, mo in enumerate(mo_results)
                 if np.isfinite(mo[4])]
    n_mo = len(valid_mos)
    K = np.zeros((n_mo, n_mo))

    if n_mo == 0:
        return K

    # --- Build 3D quadrature grid (same as compute_exact_j_integrals) ---
    eta_nodes, eta_weights = np.polynomial.legendre.leggauss(n_grid)
    xi_01_nodes, xi_01_weights = np.polynomial.legendre.leggauss(n_grid)
    xi_nodes = (xi_max - 1.0) / 2.0 * xi_01_nodes + (xi_max + 1.0) / 2.0
    xi_weights = xi_01_weights * (xi_max - 1.0) / 2.0

    XI, ETA = np.meshgrid(xi_nodes, eta_nodes)  # (n_eta, n_xi)
    W2D = np.outer(eta_weights, xi_weights)
    vol_elem = (R / 2.0)**3 * (XI**2 - ETA**2)

    # Convert to spherical coords centered at bond midpoint
    z_cart = R / 2.0 * XI * ETA
    rho_cart = R / 2.0 * np.sqrt(np.maximum((XI**2 - 1) * (1 - ETA**2), 0.0))
    r_sph = np.sqrt(z_cart**2 + rho_cart**2)
    r_sph = np.maximum(r_sph, 1e-14)
    cos_theta = z_cart / r_sph

    # --- Evaluate and normalize all MO wavefunctions ---
    psi_mos: List[np.ndarray] = []
    mo_m_vals: List[int] = []
    for k, (idx_k, mo_k) in enumerate(valid_mos):
        psi = _eval_mo_wavefunction(mo_k, xi_nodes, eta_nodes, R, p0, Z_A, Z_B)
        norm_sq = np.sum(psi**2 * vol_elem * W2D) * 2.0 * np.pi
        if norm_sq > 1e-30:
            psi /= np.sqrt(norm_sq)
        psi_mos.append(psi)
        mo_m_vals.append(mo_k[1])  # |m| value

    # --- Precompute Legendre polynomials P_l(cos_theta) ---
    from scipy.special import eval_legendre
    P_l_grid: List[np.ndarray] = []
    for l in range(l_max_poisson + 1):
        P_l_grid.append(eval_legendre(l, cos_theta))

    # --- Hartree potential via multipole expansion (same as J code) ---
    def _hartree_potential(density_2d: np.ndarray) -> np.ndarray:
        """Compute Hartree potential from 2D density via multipole expansion."""
        r_flat = r_sph.ravel()
        rho_flat = density_2d.ravel()
        vol_flat = (vol_elem * W2D).ravel() * 2.0 * np.pi

        sort_idx = np.argsort(r_flat)
        r_sorted = r_flat[sort_idx]
        rho_sorted = rho_flat[sort_idx]
        vol_sorted = vol_flat[sort_idx]

        n_pts = len(r_flat)
        V_H_flat = np.zeros(n_pts)

        for l in range(l_max_poisson + 1):
            P_l_flat = P_l_grid[l].ravel()
            rho_l_sorted = rho_sorted * P_l_flat[sort_idx]
            charge_sorted = rho_l_sorted * vol_sorted

            r_pow_l = r_sorted**l
            r_pow_lp1_inv = np.where(r_sorted > 1e-14, r_sorted**(-(l + 1)), 0.0)

            inner_terms = charge_sorted * r_pow_l
            outer_terms = charge_sorted * r_pow_lp1_inv

            A_l = np.cumsum(inner_terms)
            B_l_rev = np.cumsum(outer_terms[::-1])[::-1]

            A_l_shifted = np.zeros_like(A_l)
            A_l_shifted[1:] = A_l[:-1]
            B_l_shifted = np.zeros_like(B_l_rev)
            B_l_shifted[:-1] = B_l_rev[1:]

            r_inv_lp1 = np.where(r_sorted > 1e-14, r_sorted**(-(l + 1)), 0.0)

            V_l_sorted = A_l_shifted * r_inv_lp1 + B_l_shifted * r_pow_l

            P_l_sorted = P_l_flat[sort_idx]
            V_H_sorted = V_l_sorted * P_l_sorted

            unsort = np.argsort(sort_idx)
            V_H_flat += V_H_sorted[unsort]

        return V_H_flat.reshape(density_2d.shape)

    # --- Compute K_kl for pairs with m_k = m_l ---
    for k in range(n_mo):
        for l in range(k, n_mo):
            if mo_m_vals[k] != mo_m_vals[l]:
                continue  # K_kl = 0 by phi-orthogonality

            # Exchange density: psi_k * psi_l
            rho_kl = psi_mos[k] * psi_mos[l]

            # Hartree potential of exchange density
            V_Hkl = _hartree_potential(rho_kl)

            # K_kl = integral rho_kl * V_H[rho_kl] dV
            k_kl = np.sum(rho_kl * V_Hkl * vol_elem * W2D) * 2.0 * np.pi

            K[k, l] = k_kl
            K[l, k] = k_kl

    return K


# ---------------------------------------------------------------------------
# Cross-set overlap and H1 integrals for dual-p0 basis (v0.9.34)
# ---------------------------------------------------------------------------

def compute_cross_set_integrals(
    mo_results_A: List[Tuple[int, int, int, int, float]],
    mo_results_B: List[Tuple[int, int, int, int, float]],
    Z_A: float,
    Z_B: float,
    R: float,
    p0_A: float,
    p0_B: float,
    n_grid: int = 40,
    xi_max: float = 12.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute overlap O and H1 matrix elements between two MO sets.

    Used for dual-p0 basis: set A (e.g. Li-scale MOs at p0_A) and
    set B (e.g. H-scale MOs at p0_B).

    Sturmian identity applied with symmetrized p0:
      H1_{ij} = -(p0_A^2 + p0_B^2)/4 * O_{ij}
                + [(1 - beta_i^A) + (1 - beta_j^B)]/2 * V_{ij}

    Parameters
    ----------
    mo_results_A, mo_results_B : list of (n, m, n_sph, n_rad, beta)
        MO results from compute_molecular_sturmian_betas at p0_A and p0_B.
    Z_A, Z_B : float
        Nuclear charges (molecule, not MO set labels).
    R : float
        Internuclear distance in bohr.
    p0_A, p0_B : float
        Momentum parameters for each MO set.
    n_grid : int
        Gauss-Legendre quadrature points per dimension.
    xi_max : float
        Upper radial integration bound.

    Returns
    -------
    O_cross : ndarray, shape (nA, nB)
        Cross-set overlap matrix.
    H1_cross : ndarray, shape (nA, nB)
        Cross-set one-electron Hamiltonian matrix.
    """
    valid_A = [(i, mo) for i, mo in enumerate(mo_results_A)
               if np.isfinite(mo[4])]
    valid_B = [(i, mo) for i, mo in enumerate(mo_results_B)
               if np.isfinite(mo[4])]
    nA = len(valid_A)
    nB = len(valid_B)

    O_cross = np.zeros((nA, nB))
    H1_cross = np.zeros((nA, nB))

    if nA == 0 or nB == 0:
        return O_cross, H1_cross

    # --- Quadrature grid ---
    eta_nodes, eta_weights = np.polynomial.legendre.leggauss(n_grid)
    xi_01_nodes, xi_01_weights = np.polynomial.legendre.leggauss(n_grid)
    xi_nodes = (xi_max - 1.0) / 2.0 * xi_01_nodes + (xi_max + 1.0) / 2.0
    xi_weights = xi_01_weights * (xi_max - 1.0) / 2.0

    XI, ETA = np.meshgrid(xi_nodes, eta_nodes)
    W2D = np.outer(eta_weights, xi_weights)
    vol_elem = (R / 2.0)**3 * (XI**2 - ETA**2)

    # Nuclear potential
    r_A_grid = np.maximum(R * (XI + ETA) / 2.0, 1e-14)
    r_B_grid = np.maximum(R * (XI - ETA) / 2.0, 1e-14)
    V_nuc = -Z_A / r_A_grid - Z_B / r_B_grid

    phi_factor = 2.0 * np.pi

    # --- Evaluate and normalize A-set wavefunctions ---
    psi_A_list: List[np.ndarray] = []
    m_A: List[int] = []
    betas_A: List[float] = []
    for idx_i, mo in valid_A:
        psi = _eval_mo_wavefunction(
            mo, xi_nodes, eta_nodes, R, p0_A, Z_A, Z_B,
        )
        norm_sq = np.sum(psi**2 * vol_elem * W2D) * phi_factor
        if norm_sq > 1e-30:
            psi /= np.sqrt(norm_sq)
        psi_A_list.append(psi)
        m_A.append(mo[1])
        betas_A.append(mo[4])

    # --- Evaluate and normalize B-set wavefunctions ---
    psi_B_list: List[np.ndarray] = []
    m_B: List[int] = []
    betas_B: List[float] = []
    for idx_j, mo in valid_B:
        psi = _eval_mo_wavefunction(
            mo, xi_nodes, eta_nodes, R, p0_B, Z_A, Z_B,
        )
        norm_sq = np.sum(psi**2 * vol_elem * W2D) * phi_factor
        if norm_sq > 1e-30:
            psi /= np.sqrt(norm_sq)
        psi_B_list.append(psi)
        m_B.append(mo[1])
        betas_B.append(mo[4])

    # --- Compute cross-set overlap and potential ---
    V_cross = np.zeros((nA, nB))

    for i in range(nA):
        for j in range(nB):
            if m_A[i] != m_B[j]:
                continue  # phi orthogonality

            o_ij = np.sum(
                psi_A_list[i] * psi_B_list[j] * vol_elem * W2D
            ) * phi_factor
            O_cross[i, j] = o_ij

            v_ij = np.sum(
                psi_A_list[i] * V_nuc * psi_B_list[j] * vol_elem * W2D
            ) * phi_factor
            V_cross[i, j] = v_ij

    # --- H1 via symmetrized Sturmian identity ---
    avg_p0sq = (p0_A**2 + p0_B**2) / 2.0
    for i in range(nA):
        for j in range(nB):
            avg_beta_factor = (
                (1.0 - betas_A[i]) + (1.0 - betas_B[j])
            ) / 2.0
            H1_cross[i, j] = (
                -avg_p0sq / 2.0 * O_cross[i, j]
                + avg_beta_factor * V_cross[i, j]
            )

    return O_cross, H1_cross


def compute_combined_jk_integrals(
    mo_results_A: List[Tuple[int, int, int, int, float]],
    mo_results_B: List[Tuple[int, int, int, int, float]],
    Z_A: float,
    Z_B: float,
    R: float,
    p0_A: float,
    p0_B: float,
    n_grid: int = 30,
    xi_max: float = 12.0,
    l_max_poisson: int = 8,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute J and K integrals for combined dual-p0 basis.

    Evaluates all MO wavefunctions from both sets on a common grid,
    then computes the full (nA+nB) x (nA+nB) J and K matrices using
    Poisson multipole expansion.

    Parameters
    ----------
    mo_results_A, mo_results_B : list of (n, m, n_sph, n_rad, beta)
        MO results at p0_A and p0_B.
    Z_A, Z_B : float
        Nuclear charges.
    R : float
        Internuclear distance in bohr.
    p0_A, p0_B : float
        Momentum parameters for each MO set.
    n_grid : int
        Gauss-Legendre quadrature points per dimension.
    xi_max : float
        Upper radial integration bound.
    l_max_poisson : int
        Maximum l for Legendre expansion.

    Returns
    -------
    J : ndarray, shape (n_total, n_total)
        Direct Coulomb integral matrix.
    K : ndarray, shape (n_total, n_total)
        Exchange integral matrix.
    """
    valid_A = [(i, mo) for i, mo in enumerate(mo_results_A)
               if np.isfinite(mo[4])]
    valid_B = [(i, mo) for i, mo in enumerate(mo_results_B)
               if np.isfinite(mo[4])]
    nA = len(valid_A)
    nB = len(valid_B)
    n_total = nA + nB

    J = np.zeros((n_total, n_total))
    K = np.zeros((n_total, n_total))

    if n_total == 0:
        return J, K

    # --- Build quadrature grid ---
    eta_nodes, eta_weights = np.polynomial.legendre.leggauss(n_grid)
    xi_01_nodes, xi_01_weights = np.polynomial.legendre.leggauss(n_grid)
    xi_nodes = (xi_max - 1.0) / 2.0 * xi_01_nodes + (xi_max + 1.0) / 2.0
    xi_weights = xi_01_weights * (xi_max - 1.0) / 2.0

    XI, ETA = np.meshgrid(xi_nodes, eta_nodes)
    W2D = np.outer(eta_weights, xi_weights)
    vol_elem = (R / 2.0)**3 * (XI**2 - ETA**2)

    # Spherical coords for Poisson solve
    z_cart = R / 2.0 * XI * ETA
    rho_cart = R / 2.0 * np.sqrt(
        np.maximum((XI**2 - 1) * (1 - ETA**2), 0.0)
    )
    r_sph = np.sqrt(z_cart**2 + rho_cart**2)
    r_sph = np.maximum(r_sph, 1e-14)
    cos_theta = z_cart / r_sph

    phi_factor = 2.0 * np.pi

    # --- Evaluate all wavefunctions on common grid ---
    psi_all: List[np.ndarray] = []
    m_all: List[int] = []

    for idx_k, mo_k in valid_A:
        psi = _eval_mo_wavefunction(
            mo_k, xi_nodes, eta_nodes, R, p0_A, Z_A, Z_B,
        )
        norm_sq = np.sum(psi**2 * vol_elem * W2D) * phi_factor
        if norm_sq > 1e-30:
            psi /= np.sqrt(norm_sq)
        psi_all.append(psi)
        m_all.append(mo_k[1])

    for idx_k, mo_k in valid_B:
        psi = _eval_mo_wavefunction(
            mo_k, xi_nodes, eta_nodes, R, p0_B, Z_A, Z_B,
        )
        norm_sq = np.sum(psi**2 * vol_elem * W2D) * phi_factor
        if norm_sq > 1e-30:
            psi /= np.sqrt(norm_sq)
        psi_all.append(psi)
        m_all.append(mo_k[1])

    # --- Precompute Legendre polynomials ---
    from scipy.special import eval_legendre
    P_l_grid: List[np.ndarray] = []
    for l_val in range(l_max_poisson + 1):
        P_l_grid.append(eval_legendre(l_val, cos_theta))

    # --- Hartree potential via multipole expansion ---
    def _hartree_potential(density_2d: np.ndarray) -> np.ndarray:
        r_flat = r_sph.ravel()
        rho_flat = density_2d.ravel()
        vol_flat = (vol_elem * W2D).ravel() * phi_factor

        sort_idx = np.argsort(r_flat)
        r_sorted = r_flat[sort_idx]
        rho_sorted = rho_flat[sort_idx]
        vol_sorted = vol_flat[sort_idx]

        n_pts = len(r_flat)
        V_H_flat = np.zeros(n_pts)

        for l_val in range(l_max_poisson + 1):
            P_l_flat = P_l_grid[l_val].ravel()
            rho_l_sorted = rho_sorted * P_l_flat[sort_idx]
            charge_sorted = rho_l_sorted * vol_sorted

            r_pow_l = r_sorted**l_val
            r_pow_lp1_inv = np.where(
                r_sorted > 1e-14, r_sorted**(-(l_val + 1)), 0.0,
            )

            inner_terms = charge_sorted * r_pow_l
            outer_terms = charge_sorted * r_pow_lp1_inv

            A_l = np.cumsum(inner_terms)
            B_l_rev = np.cumsum(outer_terms[::-1])[::-1]

            A_l_shifted = np.zeros_like(A_l)
            A_l_shifted[1:] = A_l[:-1]
            B_l_shifted = np.zeros_like(B_l_rev)
            B_l_shifted[:-1] = B_l_rev[1:]

            r_inv_lp1 = np.where(
                r_sorted > 1e-14, r_sorted**(-(l_val + 1)), 0.0,
            )

            V_l_sorted = (
                A_l_shifted * r_inv_lp1 + B_l_shifted * r_pow_l
            )

            P_l_sorted = P_l_flat[sort_idx]
            V_H_sorted = V_l_sorted * P_l_sorted

            unsort = np.argsort(sort_idx)
            V_H_flat += V_H_sorted[unsort]

        return V_H_flat.reshape(density_2d.shape)

    # --- Compute J: J_kl = <rho_l | V_H[rho_k]> ---
    densities = [psi**2 for psi in psi_all]
    hartree_potentials = [_hartree_potential(d) for d in densities]

    for k in range(n_total):
        for l_idx in range(k, n_total):
            j_kl = np.sum(
                densities[l_idx] * hartree_potentials[k]
                * vol_elem * W2D
            ) * phi_factor
            J[k, l_idx] = j_kl
            J[l_idx, k] = j_kl

    # --- Compute K: K_kl = <rho_kl | V_H[rho_kl]> ---
    for k in range(n_total):
        for l_idx in range(k, n_total):
            if m_all[k] != m_all[l_idx]:
                continue  # K = 0 by phi-orthogonality

            rho_kl = psi_all[k] * psi_all[l_idx]
            V_Hkl = _hartree_potential(rho_kl)
            k_kl = np.sum(
                rho_kl * V_Hkl * vol_elem * W2D
            ) * phi_factor
            K[k, l_idx] = k_kl
            K[l_idx, k] = k_kl

    return J, K
