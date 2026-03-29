"""
Hyperradial solver for the adiabatic hyperspherical method.

Solves the 1D Schrodinger equation in the hyperradius R:
    [-1/2 d^2F/dR^2 + V_eff(R)] F(R) = E F(R)

where V_eff(R) = mu(R)/R^2 + 15/(8R^2), with mu(R) from the angular solve.

For coupled channels, solves the matrix equation:
    [-1/2 d^2/dR^2 δ_μν + V_eff_μ δ_μν - P_μν d/dR - 1/2 Q_μν] F_ν = E F_μ

Two radial methods are available:
  - 'fd' (default): Finite differences on a uniform grid (N_R ~ 2000-3000 points)
  - 'spectral': Spectral Laguerre basis (n_basis ~ 20-30 functions)
    Basis: phi_n(R) = R * exp(-alpha*R) * L_n(2*alpha*R)
    The R prefactor enforces F(0)=0; the Laguerre polynomial captures
    the radial oscillation; the exponential matches the asymptotic decay.
    Matrix elements are computed via Gauss-Laguerre quadrature.

References:
  - Macek, J. Phys. B 1, 831 (1968)
  - Lin, Phys. Rep. 257, 1 (1995)
  - Paper 11, v2.0.8 (spectral Laguerre pattern for prolate spheroidal)
"""

import numpy as np
from scipy.linalg import eigh_tridiagonal, eigh
from scipy.interpolate import CubicSpline
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import eigsh
from scipy.special import roots_laguerre, eval_laguerre
from typing import Tuple, Optional

from geovac.hyperspherical_adiabatic import (
    compute_adiabatic_curve,
    effective_potential,
)


def solve_radial(
    V_eff_func,
    R_min: float = 0.05,
    R_max: float = 30.0,
    N_R: int = 2000,
    n_states: int = 1,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Solve the hyperradial Schrodinger equation.

    [-1/2 d^2F/dR^2 + V_eff(R) F(R)] = E F(R)

    with Dirichlet boundary conditions F(R_min) = F(R_max) = 0.

    Parameters
    ----------
    V_eff_func : callable or ndarray
        Effective potential V_eff(R). If callable, evaluated on grid.
    R_min : float
        Left boundary (bohr). Must be > 0.
    R_max : float
        Right boundary (bohr).
    N_R : int
        Number of interior grid points.
    n_states : int
        Number of eigenstates to return.

    Returns
    -------
    E : ndarray of shape (n_states,)
        Energy eigenvalues (Ha).
    F : ndarray of shape (n_states, N_R)
        Radial wavefunctions on the grid.
    R_grid : ndarray of shape (N_R,)
        Grid points.
    """
    h = (R_max - R_min) / (N_R + 1)
    R_grid = R_min + (np.arange(N_R) + 1) * h

    if callable(V_eff_func):
        V = V_eff_func(R_grid)
    else:
        V = np.asarray(V_eff_func)

    # FD for -1/2 d^2/dR^2 + V(R)
    diag = np.ones(N_R) / h**2 + V
    off_diag = -0.5 * np.ones(N_R - 1) / h**2

    evals, evecs = eigh_tridiagonal(
        diag, off_diag,
        select='i', select_range=(0, n_states - 1)
    )

    for i in range(n_states):
        norm = np.sqrt(h * np.sum(evecs[:, i]**2))
        if norm > 0:
            evecs[:, i] /= norm

    return evals, evecs.T, R_grid


def _laguerre_moment_matrices(N: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build Laguerre moment matrices M_k[i,j] = int_0^inf x^k L_i(x) L_j(x) e^{-x} dx.

    Uses the three-term recurrence x*L_j = -(j+1)L_{j+1} + (2j+1)L_j - j*L_{j-1}
    to derive closed-form band structure. M0 = identity, M1 = tridiagonal, M2 = pentadiagonal.

    Same pattern as Level 2 (prolate_spheroidal_lattice.py).
    """
    M0 = np.eye(N)

    # M1: tridiagonal
    M1 = np.zeros((N, N))
    for i in range(N):
        M1[i, i] = 2 * i + 1
        if i + 1 < N:
            M1[i, i + 1] = -(i + 1)
            M1[i + 1, i] = -(i + 1)

    # M2: pentadiagonal (from two applications of the recurrence)
    M2 = np.zeros((N, N))
    for j in range(N):
        M2[j, j] = 6 * j * j + 6 * j + 2
        if j + 1 < N:
            val = -4 * (j + 1) ** 2
            M2[j, j + 1] = val
            M2[j + 1, j] = val
        if j + 2 < N:
            val = (j + 1) * (j + 2)
            M2[j, j + 2] = val
            M2[j + 2, j] = val

    return M0, M1, M2


def _build_laguerre_SK_algebraic(
    n_basis: int, alpha: float
) -> Tuple[np.ndarray, np.ndarray]:
    """Build overlap S and kinetic K matrices algebraically for Level 3.

    Basis: phi_n(R) = (R-R_min) * exp(-alpha*(R-R_min)) * L_n(2*alpha*(R-R_min))
    In x-space (x = 2*alpha*(R-R_min)): phi_n = (x/(2a)) * exp(-x/2) * L_n(x).

    Overlap: S_ij = 1/(8*alpha^3) * M2[i,j]
        where M2[i,j] = int x^2 L_i L_j exp(-x) dx (pentadiagonal).

    Kinetic: K_ij = 1/(4*alpha) * sum_k b[i,k] * b[j,k]
        where B_n(x) = (1 - x/2) L_n(x) + x L_n'(x) is the derivative kernel
        (dphi_n/dR = exp(-x/2) * B_n(x)), expanded in Laguerre polynomials as:
            B_n = -n/2 * L_{n-1} + 1/2 * L_n + (n+1)/2 * L_{n+1}
        This tridiagonal expansion follows from x*L_n'(x) = n*L_n - n*L_{n-1}
        (derived via the three-term recurrence applied to L_n' = -sum_{k<n} L_k).

    The potential V_eff(R) = mu(R)/R^2 + 15/(8R^2) is transcendental (Paper 13,
    Track G) and cannot be algebraicized — it stays as quadrature.

    Parameters
    ----------
    n_basis : int
        Number of Laguerre basis functions.
    alpha : float
        Exponential decay parameter.

    Returns
    -------
    S : ndarray of shape (n_basis, n_basis)
        Overlap matrix (pentadiagonal, exact).
    K : ndarray of shape (n_basis, n_basis)
        Kinetic energy matrix (exact).
    """
    N = n_basis
    _, _, M2 = _laguerre_moment_matrices(N)

    # Overlap: S_ij = 1/(8*alpha^3) * M2[i,j]
    S = M2 / (8.0 * alpha ** 3)

    # Kinetic: expand B_n in Laguerre basis
    # B_n = -n/2 * L_{n-1} + 1/2 * L_n + (n+1)/2 * L_{n+1}
    # b matrix is N x (N+1) since B_{N-1} has an L_N term
    b = np.zeros((N, N + 1))
    for n in range(N):
        if n > 0:
            b[n, n - 1] = -n / 2.0
        b[n, n] = 0.5
        b[n, n + 1] = (n + 1) / 2.0

    # K_ij = 1/(4*alpha) * sum_k b[i,k] b[j,k]  (orthogonality of L_k)
    K = (b @ b.T) / (4.0 * alpha)

    return S, K


def _build_laguerre_matrices_dirichlet(
    V_func,
    n_basis: int,
    alpha: float,
    R_min: float,
    matrix_method: str = 'quadrature',
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray,
           np.ndarray, np.ndarray, np.ndarray]:
    """Build Laguerre spectral matrices with Dirichlet BC at R_min.

    Domain: R in [R_min, inf), mapped via x = 2*alpha*(R - R_min).
    Basis: phi_n(R) = (R - R_min) * exp(-alpha*(R-R_min)) * L_n(x)

    The (R-R_min) prefactor enforces phi_n(R_min)=0, matching the FD
    solver's Dirichlet BC. In terms of x: phi_n = (x/2a) exp(-x/2) L_n(x).

    When matrix_method='algebraic', overlap S and kinetic K are computed
    from closed-form Laguerre recurrence relations (zero quadrature error).
    The potential V remains quadrature-based because V_eff(R) is transcendental.

    Returns S, K, V_mat, L_vals, B_vals, x_quad, w_quad, R_q.
    B_vals are the derivative kernel: dphi_n/dR = exp(-x/2) * B_n(x).
    """
    N = n_basis
    two_alpha = 2.0 * alpha

    n_quad = max(3 * N + 10, 80)
    x_quad, w_quad = roots_laguerre(n_quad)

    # Map to R: x = 2*alpha*(R - R_min)
    R_q = R_min + x_quad / two_alpha

    # Evaluate Laguerre polynomials: L[n, k]
    L_vals = np.zeros((N, n_quad))
    for n in range(N):
        L_vals[n] = eval_laguerre(n, x_quad)

    # Derivatives: L'_n(x) = L'_{n-1}(x) - L_{n-1}(x)
    dL_vals = np.zeros((N, n_quad))
    for n in range(1, N):
        dL_vals[n] = dL_vals[n - 1] - L_vals[n - 1]

    # Derivative of basis: dphi_n/dR = exp(-x/2) * B_n(x)
    # where B_n(x) = (1 - x/2) L_n(x) + x L'_n(x)
    B_vals = (1.0 - x_quad / 2.0)[np.newaxis, :] * L_vals + \
             x_quad[np.newaxis, :] * dL_vals

    inv_8a3 = 1.0 / (8.0 * alpha**3)
    inv_4a = 1.0 / (4.0 * alpha)
    x2_w = w_quad * x_quad**2

    if matrix_method == 'algebraic':
        # Overlap and kinetic from closed-form recurrence (exact)
        S, K = _build_laguerre_SK_algebraic(N, alpha)
    else:
        # Overlap: S_ij = 1/(8a^3) sum_k w_k x_k^2 L_i L_j
        x2_wL = x2_w[np.newaxis, :] * L_vals
        S = inv_8a3 * (x2_wL @ L_vals.T)

        # Kinetic: K_ij = 1/(4a) sum_k w_k B_i B_j
        wB = w_quad[np.newaxis, :] * B_vals
        K = inv_4a * (wB @ B_vals.T)

    # Potential: V_ij = 1/(8a^3) sum_k w_k x_k^2 V(R_k) L_i L_j
    # (always quadrature — V_eff(R) is transcendental)
    V_q = V_func(R_q)
    x2_wVL = (x2_w * V_q)[np.newaxis, :] * L_vals
    V_mat = inv_8a3 * (x2_wVL @ L_vals.T)

    return S, K, V_mat, L_vals, B_vals, x_quad, w_quad, R_q


def solve_radial_spectral(
    V_eff_func,
    n_basis: int = 25,
    alpha: float = 1.5,
    R_min: float = 0.05,
    n_states: int = 1,
    matrix_method: str = 'quadrature',
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Spectral Laguerre solver for the hyperradial Schrodinger equation.

    [-1/2 d^2F/dR^2 + V_eff(R) F(R)] = E F(R)

    on the domain [R_min, inf) with Dirichlet BC F(R_min)=0 (from the
    centrifugal barrier) and F(inf)=0 (exponential decay).

    Basis: phi_n(R) = (R-R_min) * exp(-alpha*(R-R_min)) * L_n(2*alpha*(R-R_min))

    The (R-R_min) prefactor enforces the Dirichlet BC, matching the FD
    solver. The domain shift avoids sampling the 1/R^2 singularity at R=0.

    Parameters
    ----------
    V_eff_func : callable
        Effective potential V_eff(R).
    n_basis : int
        Number of Laguerre basis functions.
    alpha : float
        Exponential decay parameter. Should approximate the wavefunction's
        asymptotic decay rate. For He, alpha ~ 1.0-2.0.
    R_min : float
        Left boundary of the domain (bohr). Must be > 0.
    n_states : int
        Number of eigenstates to return.
    matrix_method : str
        'quadrature' (default): all matrix elements via Gauss-Laguerre.
        'algebraic': overlap S and kinetic K from closed-form recurrence,
        potential V via quadrature (V_eff is transcendental).

    Returns
    -------
    E : ndarray of shape (n_states,)
        Energy eigenvalues (Ha).
    F : ndarray of shape (n_states, n_eval)
        Radial wavefunctions on evaluation grid.
    R_eval : ndarray of shape (n_eval,)
        Evaluation grid points.
    """
    S, K, V_mat, L_vals, B_vals, x_quad, w_quad, R_q = \
        _build_laguerre_matrices_dirichlet(
            V_eff_func, n_basis, alpha, R_min, matrix_method=matrix_method
        )
    N = n_basis
    two_alpha = 2.0 * alpha

    H_mat = K + V_mat

    # Solve generalized eigenvalue problem H c = E S c
    evals, evecs = eigh(H_mat, S)

    # Select lowest n_states
    idx = np.argsort(evals)
    evals = evals[idx[:n_states]]
    evecs = evecs[:, idx[:n_states]]

    # Reconstruct wavefunctions on an evaluation grid
    n_eval = 500
    R_eval = np.linspace(R_min, 25.0, n_eval)
    x_eval = two_alpha * (R_eval - R_min)
    L_eval = np.zeros((N, n_eval))
    for n in range(N):
        L_eval[n] = eval_laguerre(n, x_eval)

    # phi_n(R) = (R-R_min) * exp(-alpha*(R-R_min)) * L_n(x)
    s_eval = R_eval - R_min
    basis_eval = s_eval[np.newaxis, :] * np.exp(-alpha * s_eval)[np.newaxis, :] * L_eval

    F = np.zeros((n_states, n_eval))
    for s in range(n_states):
        F[s] = evecs[:, s] @ basis_eval
        dR = R_eval[1] - R_eval[0]
        norm = np.sqrt(dR * np.sum(F[s]**2))
        if norm > 0:
            F[s] /= norm

    return evals, F, R_eval


def solve_coupled_radial_spectral(
    V_eff_splines: list,
    P_splines: np.ndarray,
    Q_splines: Optional[np.ndarray] = None,
    n_channels: int = 2,
    n_basis: int = 25,
    alpha: float = 1.5,
    R_min: float = 0.05,
    n_states: int = 5,
    include_Q: bool = False,
    matrix_method: str = 'quadrature',
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Spectral Laguerre solver for the coupled-channel hyperradial equation.

    [-1/2 d^2/dR^2 delta_mu_nu + V_eff_mu delta_mu_nu
     - P_mu_nu d/dR - 1/2 Q_mu_nu] F_nu = E F_mu

    on domain [R_min, inf) with Dirichlet BC F(R_min)=0.
    Basis: phi_n(R) = (R-R_min) * exp(-alpha*(R-R_min)) * L_n(2*alpha*(R-R_min))

    Parameters
    ----------
    V_eff_splines : list of callables
        V_eff_mu(R) for each channel, length n_channels.
    P_splines : ndarray of callables, shape (n_channels, n_channels)
        P_mu_nu(R) interpolating functions.
    Q_splines : ndarray of callables, optional
        Q_mu_nu(R). Only used if include_Q=True.
    n_channels : int
        Number of coupled channels.
    n_basis : int
        Number of Laguerre basis functions per channel.
    alpha : float
        Exponential decay parameter.
    R_min : float
        Left boundary of the domain (bohr).
    n_states : int
        Number of eigenvalues to find.
    include_Q : bool
        Whether to include the Q (second-derivative) coupling term.
    matrix_method : str
        'quadrature' (default): all matrix elements via Gauss-Laguerre.
        'algebraic': overlap S and kinetic K from closed-form recurrence,
        potential V, P-coupling, Q-coupling via quadrature (R-dependent).

    Returns
    -------
    E : ndarray of shape (n_states,)
        Energy eigenvalues (Ha), sorted ascending.
    F : ndarray of shape (n_states, n_channels, n_eval)
        Channel wavefunctions on evaluation grid.
    R_eval : ndarray of shape (n_eval,)
        Evaluation grid points.
    """
    n_ch = n_channels
    N = n_basis
    two_alpha = 2.0 * alpha
    inv_8a3 = 1.0 / (8.0 * alpha**3)
    inv_4a = 1.0 / (4.0 * alpha)
    inv_4a2 = 1.0 / (4.0 * alpha**2)
    dim = n_ch * N

    # Gauss-Laguerre quadrature (still needed for V, P, Q which are R-dependent)
    n_quad = max(3 * N + 10, 80)
    x_quad, w_quad = roots_laguerre(n_quad)
    R_q = R_min + x_quad / two_alpha

    # Evaluate Laguerre polynomials and derivatives
    L_vals = np.zeros((N, n_quad))
    for n in range(N):
        L_vals[n] = eval_laguerre(n, x_quad)

    dL_vals = np.zeros((N, n_quad))
    for n in range(1, N):
        dL_vals[n] = dL_vals[n - 1] - L_vals[n - 1]

    # B_n(x) = (1 - x/2) L_n(x) + x L'_n(x)  [derivative kernel]
    B_vals = (1.0 - x_quad / 2.0)[np.newaxis, :] * L_vals + \
             x_quad[np.newaxis, :] * dL_vals

    # Common quadrature weights
    x2_w = w_quad * x_quad**2
    x_w = w_quad * x_quad

    if matrix_method == 'algebraic':
        # Overlap and kinetic from closed-form recurrence (exact)
        S_block, K_block = _build_laguerre_SK_algebraic(N, alpha)
    else:
        # Overlap block: S_ij = inv_8a3 * sum_k w_k x_k^2 L_i L_j
        x2_wL = x2_w[np.newaxis, :] * L_vals
        S_block = inv_8a3 * (x2_wL @ L_vals.T)

        # Kinetic block: K_ij = inv_4a * sum_k w_k B_i B_j
        wB = w_quad[np.newaxis, :] * B_vals
        K_block = inv_4a * (wB @ B_vals.T)

    # Build full block Hamiltonian and overlap
    H_full = np.zeros((dim, dim))
    S_full = np.zeros((dim, dim))

    for mu in range(n_ch):
        ms = mu * N
        me = (mu + 1) * N

        # Diagonal: kinetic + V_eff_mu
        V_q = V_eff_splines[mu](R_q)
        x2_wVL = (x2_w * V_q)[np.newaxis, :] * L_vals
        V_block = inv_8a3 * (x2_wVL @ L_vals.T)

        H_full[ms:me, ms:me] = K_block + V_block
        S_full[ms:me, ms:me] = S_block

        for nu in range(n_ch):
            ns = nu * N
            ne = (nu + 1) * N

            # P coupling: -int phi_i P_mn dphi_j/dR dR
            # phi_i = (x/2a) exp(-x/2) L_i, dphi_j/dR = exp(-x/2) B_j
            # = -inv_4a2 * sum_k w_k x_k P(R_k) L_i(x_k) B_j(x_k)
            P_q = P_splines[mu][nu](R_q)
            xP_wL = (x_w * P_q)[np.newaxis, :] * L_vals
            P_block = -inv_4a2 * (xP_wL @ B_vals.T)
            H_full[ms:me, ns:ne] += P_block

            # Q coupling: -1/2 int phi_i Q_mn phi_j dR
            # = -0.5 * inv_8a3 * sum_k w_k x_k^2 Q(R_k) L_i L_j
            if include_Q and Q_splines is not None:
                Q_q = Q_splines[mu][nu](R_q)
                x2_wQL = (x2_w * Q_q)[np.newaxis, :] * L_vals
                Q_block = -0.5 * inv_8a3 * (x2_wQL @ L_vals.T)
                H_full[ms:me, ns:ne] += Q_block

    # Solve generalized eigenvalue problem
    evals, evecs = eigh(H_full, S_full)

    idx = np.argsort(evals)
    n_ret = min(n_states, len(evals))
    evals = evals[idx[:n_ret]]
    evecs = evecs[:, idx[:n_ret]]

    # Reconstruct wavefunctions on evaluation grid
    n_eval = 500
    R_eval = np.linspace(R_min, 25.0, n_eval)
    x_eval = two_alpha * (R_eval - R_min)
    L_eval = np.zeros((N, n_eval))
    for n in range(N):
        L_eval[n] = eval_laguerre(n, x_eval)

    s_eval = R_eval - R_min
    basis_eval = s_eval[np.newaxis, :] * np.exp(-alpha * s_eval)[np.newaxis, :] * L_eval

    F = np.zeros((n_ret, n_ch, n_eval))
    for s in range(n_ret):
        for mu in range(n_ch):
            coeffs = evecs[mu * N:(mu + 1) * N, s]
            F[s, mu, :] = coeffs @ basis_eval

        dR = R_eval[1] - R_eval[0]
        norm_sq = dR * np.sum(F[s]**2)
        if norm_sq > 0:
            F[s] /= np.sqrt(norm_sq)

    return evals, F, R_eval


def solve_coupled_radial(
    V_eff_splines: list,
    P_splines: np.ndarray,
    Q_splines: Optional[np.ndarray] = None,
    n_channels: int = 2,
    R_min: float = 0.05,
    R_max: float = 30.0,
    N_R: int = 3000,
    n_states: int = 5,
    sigma: float = -3.0,
    include_Q: bool = False,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Solve the coupled-channel hyperradial equation.

    [-1/2 d²/dR² δ_μν + V_eff_μ(R) δ_μν - P_μν(R) d/dR (- 1/2 Q_μν)] F_ν = E F_μ

    Builds a block-tridiagonal sparse matrix and uses shift-invert eigsh.

    The Q term (second-derivative coupling) is optional. When n_channels is
    small, the closure approximation Q ≈ P² can be unreliable and is omitted
    by default. The P-only formulation captures the dominant non-adiabatic
    coupling.

    Parameters
    ----------
    V_eff_splines : list of callables
        V_eff_μ(R) for each channel μ, length n_channels.
    P_splines : ndarray of callables, shape (n_channels, n_channels)
        P_μν(R) interpolating functions.
    Q_splines : ndarray of callables, optional
        Q_μν(R) interpolating functions. Only used if include_Q=True.
    n_channels : int
        Number of coupled channels.
    R_min, R_max : float
        Radial boundaries.
    N_R : int
        Number of interior grid points.
    n_states : int
        Number of eigenvalues to find.
    sigma : float
        Shift for shift-invert mode (target energy region).
    include_Q : bool
        Whether to include the Q (second-derivative) coupling term.

    Returns
    -------
    E : ndarray of shape (n_states,)
        Energy eigenvalues (Ha), sorted ascending.
    F : ndarray of shape (n_states, n_channels, N_R)
        Channel wavefunctions on the grid.
    R_grid : ndarray of shape (N_R,)
        Grid points.
    """
    n_ch = n_channels
    h = (R_max - R_min) / (N_R + 1)
    R_grid = R_min + (np.arange(N_R) + 1) * h

    # Evaluate potentials and couplings on grid
    V_eff = np.zeros((n_ch, N_R))
    for mu in range(n_ch):
        V_eff[mu] = V_eff_splines[mu](R_grid)

    P_grid = np.zeros((n_ch, n_ch, N_R))
    for mu in range(n_ch):
        for nu in range(n_ch):
            P_grid[mu, nu] = P_splines[mu][nu](R_grid)

    Q_grid = np.zeros((n_ch, n_ch, N_R))
    if include_Q and Q_splines is not None:
        for mu in range(n_ch):
            for nu in range(n_ch):
                Q_grid[mu, nu] = Q_splines[mu][nu](R_grid)

    # Build block-tridiagonal sparse matrix of dimension (N_R * n_ch) x (N_R * n_ch)
    # Indexing: global index = i * n_ch + mu  (grid point i, channel mu)
    dim = N_R * n_ch
    H = lil_matrix((dim, dim))

    kinetic_diag = 1.0 / h**2
    kinetic_off = -0.5 / h**2

    for i in range(N_R):
        for mu in range(n_ch):
            row = i * n_ch + mu

            # Diagonal block: kinetic + V_eff (+ Q if included)
            H[row, row] = kinetic_diag + V_eff[mu, i] - 0.5 * Q_grid[mu, mu, i]

            # Off-diagonal channels at same grid point (from Q)
            if include_Q:
                for nu in range(n_ch):
                    if nu != mu:
                        col = i * n_ch + nu
                        H[row, col] += -0.5 * Q_grid[mu, nu, i]

        # Off-diagonal in R: coupling to i+1
        if i < N_R - 1:
            for mu in range(n_ch):
                row = i * n_ch + mu
                for nu in range(n_ch):
                    col = (i + 1) * n_ch + nu
                    # Kinetic: delta_mu_nu * kinetic_off
                    # P coupling: -P_mu_nu * 1/(2h)
                    val = 0.0
                    if mu == nu:
                        val += kinetic_off
                    val -= P_grid[mu, nu, i] / (2.0 * h)
                    if abs(val) > 1e-16:
                        H[row, col] = val

        # Off-diagonal in R: coupling to i-1
        if i > 0:
            for mu in range(n_ch):
                row = i * n_ch + mu
                for nu in range(n_ch):
                    col = (i - 1) * n_ch + nu
                    # Kinetic: delta_mu_nu * kinetic_off
                    # P coupling: +P_mu_nu * 1/(2h)
                    val = 0.0
                    if mu == nu:
                        val += kinetic_off
                    val += P_grid[mu, nu, i] / (2.0 * h)
                    if abs(val) > 1e-16:
                        H[row, col] = val

    H_csr = csr_matrix(H)

    # Shift-invert near sigma to find states near ground state
    evals, evecs = eigsh(H_csr, k=n_states, sigma=sigma, which='LM')

    # Sort by energy
    idx = np.argsort(evals)
    evals = evals[idx]
    evecs = evecs[:, idx]

    # Reshape eigenvectors: (N_R * n_ch, n_states) -> (n_states, n_ch, N_R)
    F = np.zeros((n_states, n_ch, N_R))
    for s in range(n_states):
        for mu in range(n_ch):
            F[s, mu, :] = evecs[mu::n_ch, s]

        # Normalize each state
        norm_sq = h * np.sum(F[s]**2)
        if norm_sq > 0:
            F[s] /= np.sqrt(norm_sq)

    return evals, F, R_grid


def solve_helium(
    Z: float = 2.0,
    l_max: int = 3,
    n_alpha: int = 100,
    N_R_angular: int = 200,
    R_min: float = 0.05,
    R_max: float = 30.0,
    N_R_radial: int = 2000,
    n_channels: int = 1,
    coupled: bool = False,
    verbose: bool = True,
    radial_method: str = 'fd',
    n_basis_radial: int = 25,
    alpha_radial: float = 2.0,
    matrix_method: str = 'quadrature',
) -> dict:
    """
    Full hyperspherical solver for He.

    Single-channel adiabatic (coupled=False) or coupled-channel (coupled=True).

    Parameters
    ----------
    Z : float
        Nuclear charge.
    l_max : int
        Maximum partial wave in angular expansion.
    n_alpha : int
        FD grid points for alpha.
    N_R_angular : int
        Number of R points for computing mu(R).
    R_min, R_max : float
        Hyperradial grid boundaries.
    N_R_radial : int
        Number of grid points for the radial solve.
    n_channels : int
        Number of adiabatic channels.
    coupled : bool
        If True, include non-adiabatic coupling between channels.
    verbose : bool
        Print progress information.
    radial_method : str
        'fd' (default) for finite-difference or 'spectral' for Laguerre basis.
    n_basis_radial : int
        Number of Laguerre basis functions (spectral method only).
    alpha_radial : float
        Exponential decay parameter for Laguerre basis (spectral method only).
    matrix_method : str
        'quadrature' (default) or 'algebraic' for spectral method.

    Returns
    -------
    result : dict
        Keys: 'energy', 'R_grid_angular', 'mu_curves', 'V_eff',
        'R_grid_radial', 'wavefunction', 'Z', 'l_max', 'n_channels',
        'coupled'.
        If coupled: also 'P_coupling', 'channel_weights'.
    """
    import time

    t0 = time.time()

    # --- Step 1: Compute adiabatic curves ---
    R_grid_ang = np.concatenate([
        np.linspace(0.1, 1.0, N_R_angular // 3),
        np.linspace(1.0, 5.0, N_R_angular // 3),
        np.linspace(5.0, R_max, N_R_angular // 3 + 1),
    ])
    R_grid_ang = np.unique(R_grid_ang)

    if verbose:
        mode = "coupled" if coupled else "single-channel"
        print(f"Computing adiabatic curves ({mode}, {n_channels} ch) "
              f"on {len(R_grid_ang)} R points...")
        print(f"  Z={Z}, l_max={l_max}, n_alpha={n_alpha}")

    if coupled and n_channels > 1:
        # Use coupling module with Hellmann-Feynman P and DBOC
        from geovac.hyperspherical_coupling import compute_coupling_matrices

        coupling = compute_coupling_matrices(
            R_grid_ang, Z, l_max, n_alpha, n_channels
        )
        mu_curves = coupling['mu']
        P = coupling['P']
        DBOC = coupling['DBOC']

        t1 = time.time()
        if verbose:
            print(f"  Angular + coupling solve: {t1 - t0:.2f}s")
            print(f"  DBOC peak (ch 0): {np.max(DBOC[0]):.6f} Ha")

        # --- Step 2: Build V_eff splines ---
        # V_eff_μ(R) = μ(R)/R² + 15/(8R²)
        # Note: DBOC is NOT added to diagonal. The off-diagonal P coupling
        # implicitly includes the DBOC through the coupled-channel equations.
        # Adding DBOC explicitly would double-count because the DBOC and
        # off-diagonal coupling nearly cancel (~97% cancellation for He).
        V_eff_splines = []
        for ch in range(n_channels):
            V_eff_ch = effective_potential(R_grid_ang, mu_curves[ch])
            V_eff_splines.append(
                CubicSpline(R_grid_ang, V_eff_ch, extrapolate=True)
            )

        P_splines = [[None] * n_channels for _ in range(n_channels)]
        for mu_idx in range(n_channels):
            for nu_idx in range(n_channels):
                P_splines[mu_idx][nu_idx] = CubicSpline(
                    R_grid_ang, P[mu_idx, nu_idx], extrapolate=True
                )

        # --- Step 3: Solve coupled radial equation ---
        if radial_method == 'spectral':
            if verbose:
                print(f"Solving coupled-channel radial equation "
                      f"(spectral, n_basis={n_basis_radial}, {n_channels} ch)...")
            E_all, F_all, R_grid_rad = solve_coupled_radial_spectral(
                V_eff_splines, P_splines,
                n_channels=n_channels,
                n_basis=n_basis_radial, alpha=alpha_radial,
                R_min=R_min,
                n_states=min(5, n_channels * 3),
                matrix_method=matrix_method,
            )
        else:
            if verbose:
                print(f"Solving coupled-channel radial equation "
                      f"(N_R={N_R_radial}, {n_channels} ch)...")
            E_all, F_all, R_grid_rad = solve_coupled_radial(
                V_eff_splines, P_splines,
                n_channels=n_channels,
                R_min=R_min, R_max=R_max, N_R=N_R_radial,
                n_states=min(5, n_channels * 3),
                sigma=-3.0,
            )

        t2 = time.time()
        E_exact = -2.903724

        # Channel weights for ground state
        h_rad = R_grid_rad[1] - R_grid_rad[0]
        weights = np.array([
            h_rad * np.sum(F_all[0, ch]**2)
            for ch in range(n_channels)
        ])
        weights /= weights.sum()

        if verbose:
            print(f"  Radial solve: {t2 - t1:.2f}s")
            print(f"\n  Ground state energy: {E_all[0]:.6f} Ha")
            print(f"  Exact He:            {E_exact:.6f} Ha")
            err = abs(E_all[0] - E_exact) / abs(E_exact) * 100
            print(f"  Error:               {err:.2f}%")
            print(f"  Channel weights:     {weights}")
            print(f"  Total time:          {t2 - t0:.2f}s")

        V_eff_angular = effective_potential(R_grid_ang, mu_curves[0])

        return {
            'energy': E_all[0],
            'energies': E_all,
            'R_grid_angular': R_grid_ang,
            'mu_curves': mu_curves,
            'V_eff': V_eff_angular,
            'DBOC': DBOC,
            'R_grid_radial': R_grid_rad,
            'wavefunction': F_all[0],
            'channel_weights': weights,
            'P_coupling': P,
            'Z': Z,
            'l_max': l_max,
            'n_channels': n_channels,
            'coupled': True,
        }

    else:
        # Single-channel adiabatic (original code path)
        mu = compute_adiabatic_curve(
            R_grid_ang, Z, l_max, n_alpha, n_channels
        )

        t1 = time.time()
        if verbose:
            print(f"  Angular solve: {t1 - t0:.2f}s")
            V_eff_last = effective_potential(
                R_grid_ang[-1:], mu[0, -1:]
            )[0]
            print(f"  V_eff(R={R_grid_ang[-1]:.1f}) = {V_eff_last:.4f} Ha "
                  f"(should approach {-Z**2/2:.1f})")

        # --- Step 2: Interpolate V_eff for radial solver ---
        V_eff_angular = effective_potential(R_grid_ang, mu[0])
        V_eff_spline = CubicSpline(
            R_grid_ang, V_eff_angular, extrapolate=True
        )

        # --- Step 3: Solve hyperradial equation ---
        if radial_method == 'spectral':
            if verbose:
                print(f"Solving hyperradial equation "
                      f"(spectral, n_basis={n_basis_radial})...")
            E, F, R_grid_rad = solve_radial_spectral(
                V_eff_spline, n_basis=n_basis_radial,
                alpha=alpha_radial, R_min=R_min, n_states=1,
                matrix_method=matrix_method,
            )
        else:
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
            print(f"  Error:               {err:.2f}%")
            print(f"  Total time:          {t2 - t0:.2f}s")

        return {
            'energy': E[0],
            'R_grid_angular': R_grid_ang,
            'mu_curves': mu,
            'V_eff': V_eff_angular,
            'R_grid_radial': R_grid_rad,
            'wavefunction': F[0],
            'Z': Z,
            'l_max': l_max,
            'n_channels': n_channels,
            'coupled': False,
        }
