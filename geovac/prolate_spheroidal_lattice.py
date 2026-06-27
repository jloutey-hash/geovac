"""
Prolate Spheroidal Lattice for H2+ and one-electron diatomics.

The molecular analog of Fock's S3 projection: the two-center Coulomb
problem separates in prolate spheroidal coordinates (xi, eta), and the
eigenvalues are obtained by matching separated 1D equations — no free
parameters beyond grid resolution.

Coordinate system (nuclei A, B at foci, distance R apart):
  xi  in [1, inf)  — confocal ellipsoidal (r_A + r_B = R*xi)
  eta in [-1, +1]  — confocal hyperboloidal (r_A - r_B = R*eta)
  phi in [0, 2*pi) — azimuthal

Separated equations for m=0 (sigma states):
  xi:  d/dxi[(xi^2-1)dF/dxi] + (A + a*xi - c^2*xi^2) F = 0
  eta: d/deta[(1-eta^2)dG/deta] + (-A + c^2*eta^2 + b*eta) G = 0

where a = R*(Z_A+Z_B), b = R*(Z_B-Z_A), c^2 = -R^2*E_elec/2.
Self-consistency: both equations share the separation constant A.

The angular equation is solved spectrally (Legendre basis, exact).
The radial equation is solved with self-adjoint FD (Neumann at xi=1).
Root-finding in c^2 matches the two equations.

Reference: Bates, Ledsham & Stewart, Phil. Trans. A 246, 215 (1953).
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal, eigh
from scipy.optimize import brentq
from scipy.special import (
    roots_laguerre, eval_laguerre, roots_genlaguerre, eval_genlaguerre,
    gammaln, exp1
)
from typing import Dict, Optional, Tuple
import time

from geovac.molecular_sturmian import _angular_sep_const


def _laguerre_moment_matrices(N: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build Laguerre moment matrices M_k[i,j] = int_0^inf x^k L_i(x) L_j(x) e^{-x} dx.

    Uses the three-term recurrence x*L_j = -(j+1)L_{j+1} + (2j+1)L_j - j*L_{j-1}
    to derive closed-form band structure. M0 = identity, M1 = tridiagonal, M2 = pentadiagonal.

    Parameters
    ----------
    N : int
        Matrix dimension (number of Laguerre basis functions).

    Returns
    -------
    M0, M1, M2 : np.ndarray
        Moment matrices of shape (N, N).
    """
    M0 = np.eye(N)

    # M1: tridiagonal
    # M1[i,i] = 2i+1, M1[i,i+1] = M1[i+1,i] = -(i+1)
    M1 = np.zeros((N, N))
    for i in range(N):
        M1[i, i] = 2 * i + 1
        if i + 1 < N:
            M1[i, i + 1] = -(i + 1)
            M1[i + 1, i] = -(i + 1)

    # M2: pentadiagonal
    # From x^2 * L_j expanded via two applications of the recurrence:
    # i=j:   6j^2 + 6j + 2
    # i=j+1: -4(j+1)^2      (and symmetric)
    # i=j+2: (j+1)(j+2)     (and symmetric)
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


def _associated_laguerre_moment_matrices(
    N: int, alpha: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build associated Laguerre moment matrices.

    M_k^(alpha)[i,j] = int_0^inf x^(alpha+k) L_i^(alpha)(x) L_j^(alpha)(x) e^{-x} dx

    where L_n^(alpha) are the associated Laguerre polynomials with weight x^alpha e^{-x}.

    Uses the three-term recurrence:
        x L_n^(a) = -(n+1) L_{n+1}^(a) + (2n+a+1) L_n^(a) - (n+a) L_{n-1}^(a)

    Parameters
    ----------
    N : int
        Matrix dimension (number of basis functions).
    alpha : int
        Associated Laguerre parameter (|m| for centrifugal applications).

    Returns
    -------
    M0, M1, M2 : np.ndarray
        Moment matrices of shape (N, N). M0 diagonal, M1 tridiagonal, M2 pentadiagonal.
    """
    if alpha == 0:
        return _laguerre_moment_matrices(N)

    # h_n^(alpha) = Gamma(n + alpha + 1) / n!
    h = np.array([np.exp(gammaln(n + alpha + 1) - gammaln(n + 1)) for n in range(N)])

    # M0: diagonal orthogonality
    M0 = np.diag(h)

    # M1: tridiagonal from single recurrence application
    # Diagonal: (2n + alpha + 1) h_n
    # Off-diagonal: -(n + alpha + 1) h_n  [for M1[n, n+1]]
    M1 = np.zeros((N, N))
    for n in range(N):
        M1[n, n] = (2 * n + alpha + 1) * h[n]
        if n + 1 < N:
            off = -(n + alpha + 1) * h[n]
            M1[n, n + 1] = off
            M1[n + 1, n] = off

    # M2: pentadiagonal from double recurrence application
    # x^2 L_j = (j+1)(j+2) L_{j+2}
    #         - 2(j+1)(2j+alpha+2) L_{j+1}
    #         + [6j^2 + 6(alpha+1)j + (alpha+1)(alpha+2)] L_j
    #         - 2(j+alpha)(2j+alpha) L_{j-1}
    #         + (j+alpha)(j+alpha-1) L_{j-2}
    M2 = np.zeros((N, N))
    a = alpha
    for j in range(N):
        M2[j, j] = (6 * j * j + 6 * (a + 1) * j + (a + 1) * (a + 2)) * h[j]
        if j + 1 < N:
            off1 = -2 * (j + a + 1) * (2 * j + a + 2) * h[j]
            M2[j, j + 1] = off1
            M2[j + 1, j] = off1
        if j + 2 < N:
            off2 = (j + a + 2) * (j + a + 1) * h[j]
            M2[j, j + 2] = off2
            M2[j + 2, j] = off2

    return M0, M1, M2


def _lowered_moment_matrix(N: int, alpha: int) -> np.ndarray:
    """Build the lowered moment matrix for associated Laguerre polynomials.

    M_{-1}^(alpha)[i,j] = int_0^inf x^(alpha-1) L_i^(alpha)(x) L_j^(alpha)(x) e^{-x} dx

    Requires alpha >= 1 (finite for |m| >= 1, exactly the centrifugal cases).
    Uses the identity L_n^(alpha)(x) = sum_{k=0}^{n} L_k^(alpha-1)(x) (DLMF 18.9.13)
    to reduce to orthogonality integrals of L^(alpha-1).

    M_{-1}[i,j] = sum_{p=0}^{min(i,j)} h_p^(alpha-1)

    where h_p^(alpha-1) = Gamma(p + alpha) / p!.

    Parameters
    ----------
    N : int
        Matrix dimension.
    alpha : int
        Associated Laguerre parameter (must be >= 1).

    Returns
    -------
    M_neg1 : np.ndarray
        Lowered moment matrix of shape (N, N).
    """
    if alpha < 1:
        raise ValueError(f"Lowered moment requires alpha >= 1, got {alpha}")

    # h_p^(alpha-1) = Gamma(p + alpha) / p!
    h_low = np.array([
        np.exp(gammaln(p + alpha) - gammaln(p + 1)) for p in range(N)
    ])

    # Cumulative sum: cum[k] = sum_{p=0}^{k} h_p^(alpha-1)
    cum_h = np.cumsum(h_low)

    # M_{-1}[i,j] = cum[min(i,j)]
    M = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            M[i, j] = cum_h[min(i, j)]

    return M


def _stieltjes_matrix(N: int, alpha: int, a: float) -> np.ndarray:
    """Build the Stieltjes integral matrix for associated Laguerre polynomials.

    J^(alpha)[i,j] = int_0^inf x^alpha e^{-x} L_i^(alpha)(x) L_j^(alpha)(x) / (x + a) dx

    where a > 0. Uses Miller's backward recurrence for the single-polynomial
    Stieltjes transforms S_n(a), then a hybrid forward/backward recurrence
    for the bilinear matrix J[i,j].

    The base case S_0^(alpha)(a) is computed from e^a E_1(a) via the
    recurrence S_0^(k)(a) = Gamma(k) - a * S_0^(k-1)(a).

    Accuracy: < 1e-12 for a >= 2 at any N <= 30. For a < 2, accuracy
    degrades gradually (~1e-8 at a=1, N=30) due to forward recurrence
    in the left segment. For centrifugal applications where a = 4*alpha_basis,
    typical a >= 2 and accuracy is machine-level.

    Parameters
    ----------
    N : int
        Matrix dimension.
    alpha : int
        Associated Laguerre parameter (>= 0).
    a : float
        Shift parameter (must be > 0). For centrifugal decomposition, a = 4*beta.

    Returns
    -------
    J : np.ndarray
        Stieltjes integral matrix of shape (N, N), real symmetric.
    """
    if a <= 0:
        raise ValueError(f"Shift parameter a must be > 0, got {a}")

    import math

    # h_n^(alpha) = Gamma(n + alpha + 1) / n!
    h = np.array([np.exp(gammaln(n + alpha + 1) - gammaln(n + 1)) for n in range(N)])

    # Step 1: Compute S_0^(alpha)(a) analytically
    # Base: S_0^(0)(a) = e^a * E_1(a). For large a, use asymptotic series
    # to avoid overflow: e^a E_1(a) = sum_{k=0}^K (-1)^k k!/a^{k+1}
    if a > 500:
        S0_val = 0.0
        term = 1.0 / a
        for k in range(30):
            S0_val += term
            term *= -(k + 1) / a
            if abs(term) < 1e-16 * abs(S0_val):
                break
    else:
        S0_val = float(np.exp(a) * exp1(a))
    for k in range(1, alpha + 1):
        S0_val = math.gamma(k) - a * S0_val

    # S_1 from the n=0 inhomogeneous equation
    S1_val = (alpha + 1 + a) * S0_val - math.gamma(alpha + 1)

    # Step 2: Miller's backward recurrence for S_n (stable for minimal solution)
    # For n >= 1, S_n satisfies the homogeneous recurrence:
    # (n+1) S_{n+1} = (2n+alpha+1+a) S_n - (n+alpha) S_{n-1}
    M_pad = N + 60
    f = np.zeros(M_pad + 1)
    f[M_pad] = 0.0
    f[M_pad - 1] = 1.0
    for n in range(M_pad - 1, 0, -1):
        f[n - 1] = ((2 * n + alpha + 1 + a) * f[n]
                    - (n + 1) * f[n + 1]) / (n + alpha)

    S = np.zeros(N)
    S[0] = S0_val
    if N > 1 and abs(f[1]) > 1e-300:
        scale = S1_val / f[1]
        for n in range(1, N):
            S[n] = f[n] * scale

    # Step 3: Hybrid forward/backward recurrence for J[i,j]
    # Left segment (j <= i): forward from j=0, accurate for i steps
    # Right segment (j > i): backward from j=M, stable for minimal solution
    J = np.zeros((N, N))
    M_J = N + 60

    for i in range(N):
        # Forward from j=0 to j=min(i+1, N-1)
        row_fwd = np.zeros(N)
        row_fwd[0] = S[i]
        if N > 1:
            row_fwd[1] = (alpha + 1 + a) * S[i] - h[i] * (1.0 if i == 0 else 0.0)
        for j in range(1, min(i + 1, N - 1)):
            row_fwd[j + 1] = ((2 * j + alpha + 1 + a) * row_fwd[j]
                              - (j + alpha) * row_fwd[j - 1]
                              - h[i] * (1.0 if i == j else 0.0)) / (j + 1)

        # Store left segment
        for j in range(min(i + 2, N)):
            J[i, j] = row_fwd[j]

        # Right segment: backward recurrence of homogeneous equation for j > i
        if i + 2 < N:
            # Bridge value from one forward step past i
            bridge = row_fwd[i + 1] if i + 1 < N else 0.0

            g = np.zeros(M_J + 1)
            g[M_J] = 0.0
            g[M_J - 1] = 1.0
            for j in range(M_J - 1, i + 1, -1):
                g[j - 1] = ((2 * j + alpha + 1 + a) * g[j]
                            - (j + 1) * g[j + 1]) / (j + alpha)

            if abs(g[i + 1]) > 1e-300:
                sc = bridge / g[i + 1]
                for j in range(i + 2, N):
                    J[i, j] = g[j] * sc

    # Enforce symmetry
    J = 0.5 * (J + J.T)

    return J


def _build_laguerre_matrices_algebraic(
    n_basis: int, alpha: float, A: float, a_param: float, c2: float, m: int
) -> Tuple[np.ndarray, np.ndarray]:
    """Build H and S matrices for the radial eigenvalue problem using algebraic
    Laguerre moment recurrence relations (no numerical quadrature).

    Valid for m=0 (sigma states). For m!=0, raises NotImplementedError because
    the m^2/(xi^2-1) term introduces a 1/x singularity in x-space that has no
    finite Laguerre moment expansion.

    Parameters
    ----------
    n_basis : int
        Number of Laguerre basis functions.
    alpha : float
        Exponential decay parameter for the basis.
    A : float
        Angular separation constant.
    a_param : float
        R*(Z_A + Z_B).
    c2 : float
        Separation parameter c^2 = -R^2*E/2.
    m : int
        Azimuthal quantum number.

    Returns
    -------
    H_mat, S : np.ndarray
        Hamiltonian and overlap matrices, shape (n_basis, n_basis).
    """
    if m != 0:
        return _build_laguerre_matrices_algebraic_associated(
            n_basis, alpha, A, a_param, c2, m
        )

    N = n_basis
    two_alpha = 2.0 * alpha
    inv_2a = 1.0 / two_alpha

    # Moment matrices
    M0, M1, M2 = _laguerre_moment_matrices(N)

    # --- Overlap: S_mn = delta_{mn} / (2*alpha) ---
    S = inv_2a * M0

    # --- Potential: V_mn = (1/(2a)) * [q0*M0 + q1*M1 + q2*M2]_{mn} ---
    # q(xi) = A + a_param*xi - c2*xi^2, with xi = 1 + x/(2*alpha)
    # q(x) = q0 + q1*x + q2*x^2
    q0 = A + a_param - c2
    q1 = (a_param - 2.0 * c2) / two_alpha
    q2 = -c2 / (two_alpha ** 2)
    V = inv_2a * (q0 * M0 + q1 * M1 + q2 * M2)

    # --- Kinetic: K_mn = -(1/2)*F*M1*F^T - (1/(8*alpha))*F*M2*F^T ---
    # where F[n,j] = 2 for j < n, 1 for j = n, 0 for j > n
    # (Laguerre expansion coefficients of f_n = L_n + 2*sum_{j<n} L_j)
    F = np.zeros((N, N))
    for n in range(N):
        F[n, :n] = 2.0
        F[n, n] = 1.0

    FM1 = F @ M1
    FM2 = F @ M2
    K = -0.5 * (FM1 @ F.T) - (1.0 / (8.0 * alpha)) * (FM2 @ F.T)

    H_mat = K + V
    return H_mat, S


def _build_laguerre_matrices_algebraic_associated(
    n_basis: int, alpha: float, A: float, a_param: float, c2: float, m: int
) -> Tuple[np.ndarray, np.ndarray]:
    """Build H and S matrices for m!=0 using associated Laguerre basis with algebraic moments.

    Uses basis functions phi_n(xi) = (xi-1)^{|m|/2} exp(-alpha*(xi-1)) L_n^{|m|}(2*alpha*(xi-1))
    which incorporate the correct (xi-1)^{|m|/2} behavior at the singular point xi=1.

    All matrix elements are computed from the algebraic moment matrices (M0, M1, M2),
    the lowered moment matrix (M_{-1}), and the Stieltjes matrix (J) — no spatial
    quadrature is used.

    The kinetic matrix is derived from integration by parts:
        K_{ij} = -int_1^inf (xi^2-1) phi_i'(xi) phi_j'(xi) dxi
    Using the identity x dL_n^(a)/dx = n L_n^(a) - (n+a) L_{n-1}^(a) and the
    product rule, the derivative decomposes into g_n = sum_k (G0_{nk} + G1_{nk} x) L_k
    with bidiagonal G0 and diagonal G1. The resulting integral reduces to:
        K = -beta^{-(alpha_L+1)} [G0 P0 G0^T + G0 P1 G1^T + G1 P1 G0^T + G1 P2 G1^T]
    where P_k = M_k + 2*beta * M_{k-1}.

    The centrifugal term -m^2/(xi^2-1) uses the "lowered Stieltjes" identity:
        int x^{a-1} e^{-x} L_i L_j / (x+a) dx = (M_{-1} - J(alpha, a)) / a

    Parameters
    ----------
    n_basis : int
        Number of associated Laguerre basis functions.
    alpha : float
        Exponential decay parameter for the basis.
    A : float
        Angular separation constant.
    a_param : float
        R*(Z_A + Z_B).
    c2 : float
        Separation parameter c^2 = -R^2*E/2.
    m : int
        Azimuthal quantum number (must be != 0).

    Returns
    -------
    H_mat, S : np.ndarray
        Hamiltonian and overlap matrices, shape (n_basis, n_basis).
    """
    N = n_basis
    alpha_L = abs(m)       # Associated Laguerre parameter
    sigma = abs(m) / 2.0   # Power-law exponent at xi=1
    beta = 2.0 * alpha     # Laguerre scale parameter

    # Moment matrices from Sub-agent 1 infrastructure
    M0, M1, M2 = _associated_laguerre_moment_matrices(N, alpha_L)
    M_neg1 = _lowered_moment_matrix(N, alpha_L)

    # --- Overlap: S = beta^{-(alpha_L+1)} * M0 ---
    pref = beta ** (-(alpha_L + 1))
    S = pref * M0

    # --- Regular potential: V_reg = pref * (q0*M0 + q1*M1 + q2*M2) ---
    # q(xi) = A + a_param*xi - c2*xi^2, with xi = 1 + x/beta
    q0 = A + a_param - c2
    q1 = (a_param - 2.0 * c2) / beta
    q2 = -c2 / beta ** 2
    V_reg = pref * (q0 * M0 + q1 * M1 + q2 * M2)

    # --- Centrifugal potential: -m^2/(xi^2-1) ---
    # V_cent = -m^2 * beta^{1-alpha_L} * J_{-1}(alpha_L, 2*beta)
    # where J_{-1}(a, shift) = (M_{-1} - J(a, shift)) / shift
    shift = 2.0 * beta
    J_stielt = _stieltjes_matrix(N, alpha_L, shift)
    V_cent = -m ** 2 * beta ** (1 - alpha_L) * (M_neg1 - J_stielt) / shift

    # --- Kinetic: K = -pref * [G0 P0 G0^T + G0 P1 G1^T + G1 P1 G0^T + G1 P2 G1^T] ---
    # Derivative decomposition: g_n(x) = sum_k (G0[n,k] + G1[n,k]*x) L_k^(alpha_L)(x)
    # G0: bidiagonal — G0[n,n] = sigma+n, G0[n,n-1] = -(n+alpha_L) for n>=1
    # G1: diagonal — G1[n,n] = -1/2
    G0 = np.zeros((N, N))
    for n in range(N):
        G0[n, n] = sigma + n
        if n >= 1:
            G0[n, n - 1] = -(n + alpha_L)
    G1 = np.diag(np.full(N, -0.5))

    # Weighted moment matrices: P_k = M_k + 2*beta * M_{k-1}
    P0 = M0 + 2.0 * beta * M_neg1
    P1 = M1 + 2.0 * beta * M0
    P2 = M2 + 2.0 * beta * M1

    I_mat = (G0 @ P0 @ G0.T + G0 @ P1 @ G1.T
             + G1 @ P1 @ G0.T + G1 @ P2 @ G1.T)
    K = -pref * I_mat

    H_mat = K + V_reg + V_cent
    return H_mat, S


class ProlateSpheroidalLattice:
    """
    Separated prolate spheroidal solver for one-electron diatomics.

    Solves H2+ (or HeH2+, etc.) by matching the separated xi and eta
    equations in prolate spheroidal coordinates. No free parameters.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    Z_A : int
        Charge of nucleus A.
    Z_B : int
        Charge of nucleus B.
    N_xi : int
        Number of radial grid points.
    xi_max : float
        Maximum xi value.
    m : int
        Azimuthal quantum number (0 for sigma, ±1 for pi, etc.).
    n_angular : int
        Angular excitation index (0 for gerade ground, 1 for ungerade, etc.).
        Controls which eigenvalue of the angular equation is used.
    n_radial : int
        Radial excitation index (0 for no radial nodes, 1 for one node, etc.).
        Controls which eigenvalue of the radial equation must vanish.
    """

    def __init__(
        self,
        R: float,
        Z_A: int = 1,
        Z_B: int = 1,
        N_xi: int = 5000,
        xi_max: float = 25.0,
        m: int = 0,
        n_angular: int = 0,
        n_radial: int = 0,
        radial_method: str = 'fd',
        n_basis: int = 20,
        matrix_method: str = 'quadrature',
    ):
        self.R = R
        self.Z_A = Z_A
        self.Z_B = Z_B
        self.N_xi = N_xi
        self.xi_max = xi_max
        self.m = m
        self.n_angular = n_angular
        self.n_radial = n_radial
        self.radial_method = radial_method
        self.n_basis = n_basis
        self.matrix_method = matrix_method

        # Derived
        self._a = R * (Z_A + Z_B)
        self._b = R * (Z_B - Z_A)

    def _radial_top_eigenvalue(
        self, c2: float, A: float
    ) -> float:
        """Top eigenvalue of the radial operator.

        L_xi = d/dxi[(xi^2-1)dF/dxi] + (A + a*xi - c^2*xi^2 - m^2/(xi^2-1))

        Returns the largest eigenvalue; should be 0 at the correct c^2.
        """
        N = self.N_xi
        xi_min = 1.0 + 5e-4
        h = (self.xi_max - xi_min) / (N + 1)
        xi = xi_min + (np.arange(N) + 1) * h

        xi2_1 = xi**2 - 1
        q = A + self._a * xi - c2 * xi**2
        if self.m != 0:
            q -= self.m**2 / xi2_1

        p_plus = (xi + h / 2)**2 - 1
        p_minus = (xi - h / 2)**2 - 1

        diag = -(p_plus + p_minus) / h**2 + q
        off = p_plus[:-1] / h**2

        # Neumann BC at xi=1 for m=0
        if self.m == 0:
            diag[0] = -p_plus[0] / h**2 + q[0]

        evals = eigh_tridiagonal(diag, off, eigvals_only=True)
        # n_radial=0 → top eigenvalue, n_radial=1 → second from top, etc.
        sorted_desc = np.sort(evals)[::-1]
        if self.n_radial >= len(sorted_desc):
            raise ValueError(
                f"n_radial={self.n_radial} exceeds available eigenvalues"
            )
        return sorted_desc[self.n_radial]

    def _radial_top_eigenvalue_spectral(
        self, c2: float, A: float
    ) -> float:
        """Spectral Laguerre basis solver for the radial xi-equation.

        Replaces the FD grid with a Galerkin expansion in Laguerre
        polynomials: phi_n(xi) = exp(-alpha*(xi-1)) * L_n(2*alpha*(xi-1)).

        The Laguerre basis is orthogonal with weight exp(-x) on [0, inf),
        so the overlap matrix S is diagonal.  Matrix elements can be evaluated
        via Gauss-Laguerre quadrature or algebraically via the three-term
        recurrence (matrix_method='algebraic').

        Returns the largest eigenvalue of the generalized eigenvalue
        problem H c = lambda S c.  This should be 0 at the correct c^2.
        """
        N = self.n_basis
        alpha = max(np.sqrt(max(c2, 0.01)), 0.5)

        if self.matrix_method == 'algebraic':
            H_mat, S = _build_laguerre_matrices_algebraic(
                N, alpha, A, self._a, c2, self.m
            )
        else:
            H_mat, S = self._build_laguerre_matrices_quadrature(
                N, alpha, A, c2
            )

        # Solve generalized eigenvalue problem H c = lambda S c
        evals = eigh(H_mat, S, eigvals_only=True)
        sorted_desc = np.sort(evals)[::-1]
        if self.n_radial >= len(sorted_desc):
            raise ValueError(
                f"n_radial={self.n_radial} exceeds available eigenvalues"
            )
        return sorted_desc[self.n_radial]

    def _build_laguerre_matrices_quadrature(
        self, N: int, alpha: float, A: float, c2: float
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Build H and S matrices via Gauss-Laguerre quadrature."""
        two_alpha = 2.0 * alpha

        # Gauss-Laguerre quadrature: integrates int_0^inf f(x) e^{-x} dx
        n_quad = max(3 * N, 60)
        x_quad, w_quad = roots_laguerre(n_quad)

        # Map to xi: x = 2*alpha*(xi - 1), xi = 1 + x/(2*alpha)
        xi_q = 1.0 + x_quad / two_alpha
        xi2m1 = xi_q**2 - 1.0  # p(xi) = xi^2 - 1

        # Evaluate Laguerre polynomials at quadrature points: L[n, k]
        L_vals = np.zeros((N, n_quad))
        for n in range(N):
            L_vals[n] = eval_laguerre(n, x_quad)

        # Derivatives of Laguerre polynomials: L'_n(x) = L'_{n-1}(x) - L_{n-1}(x)
        dL_vals = np.zeros((N, n_quad))
        for n in range(1, N):
            dL_vals[n] = dL_vals[n - 1] - L_vals[n - 1]

        # Derivative of basis function (without exponential factor):
        # d/dxi [e^{-alpha*s} L_n(x)] = e^{-alpha*s} * D_n(x)
        # where D_n(x) = -alpha * L_n(x) + 2*alpha * L'_n(x)
        D_vals = -alpha * L_vals + two_alpha * dL_vals

        # Build matrices via quadrature (factor 1/(2*alpha) from dxi = dx/(2*alpha))
        inv_2a = 1.0 / two_alpha

        # Overlap: S_ij = (1/2a) sum_k w_k L_i(x_k) L_j(x_k) = delta_ij / (2a)
        # Use quadrature for generality (handles m != 0 basis modifications)
        wL = w_quad[np.newaxis, :] * L_vals  # (N, n_quad)
        S = inv_2a * (wL @ L_vals.T)

        # Kinetic (integration by parts):
        # K_ij = -(1/2a) sum_k w_k (xi^2-1) D_i(x_k) D_j(x_k)
        wD = w_quad[np.newaxis, :] * D_vals * xi2m1[np.newaxis, :]
        K = -inv_2a * (wD @ D_vals.T)

        # Potential: V_ij = (1/2a) sum_k w_k q(xi_k) L_i(x_k) L_j(x_k)
        # where q(xi) = A + a*xi - c^2*xi^2
        q = A + self._a * xi_q - c2 * xi_q**2
        if self.m != 0:
            q -= self.m**2 / xi2m1
        wLq = w_quad[np.newaxis, :] * L_vals * q[np.newaxis, :]
        V = inv_2a * (wLq @ L_vals.T)

        H_mat = K + V
        return H_mat, S

    def _angular_separation_constant(self, c2: float) -> float:
        """Angular separation constant A for given c^2.

        Uses the spectral (Legendre basis) solver from molecular_sturmian.
        n_angular=0 → largest A (gerade), n_angular=1 → next (ungerade), etc.
        """
        c = np.sqrt(max(c2, 1e-15))
        return _angular_sep_const(
            abs(self.m), self.n_angular, c, b=self._b, n_basis=50
        )

    def _residual(self, c2: float) -> float:
        """Residual for root-finding: should be 0 at the correct c^2."""
        A = self._angular_separation_constant(c2)
        if self.radial_method == 'spectral':
            return self._radial_top_eigenvalue_spectral(c2, A)
        return self._radial_top_eigenvalue(c2, A)

    def solve(self) -> Tuple[float, float, float]:
        """Solve for the ground state electronic energy.

        Returns
        -------
        E_elec : float
            Electronic energy (Ha), not including V_NN.
        c2 : float
            Separation parameter c^2 = -R^2*E_elec/2.
        A : float
            Angular separation constant.
        """
        # Adaptive search range for c^2
        # c^2 = -R^2*E/2 > 0 for bound states
        # Upper bound: E can't be more negative than -2*(Z_A+Z_B)^2/R^2 (united atom)
        c2_max = self.R**2 * (self.Z_A + self.Z_B)**2
        c2_min = 0.01

        # Verify bracket
        f_lo = self._residual(c2_min)
        if f_lo < 0:
            raise ValueError(
                f"Residual negative at c2_min={c2_min}. "
                "No bound state found."
            )

        # Find upper bracket
        f_hi = self._residual(c2_max)
        if f_hi > 0:
            # Extend search
            for c2_try in np.linspace(c2_max, c2_max * 2, 20):
                if self._residual(c2_try) < 0:
                    c2_max = c2_try
                    break
            else:
                raise ValueError(
                    f"Could not bracket root up to c2={c2_max*2}"
                )

        c2_sol = brentq(self._residual, c2_min, c2_max, xtol=1e-12)
        E_elec = -2.0 * c2_sol / self.R**2
        A = self._angular_separation_constant(c2_sol)

        return E_elec, c2_sol, A

    def total_energy(self) -> float:
        """Electronic energy + nuclear repulsion."""
        E_elec, _, _ = self.solve()
        return E_elec + self.Z_A * self.Z_B / self.R

    def _radial_wavefunction_fd(
        self, c2: float, A: float
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Compute radial wavefunction F(xi) on the FD grid.

        Returns (xi_grid, F_xi) tuple.
        """
        N = self.N_xi
        xi_min = 1.0 + 5e-4
        h = (self.xi_max - xi_min) / (N + 1)
        xi = xi_min + (np.arange(N) + 1) * h

        xi2_1 = xi**2 - 1
        q = A + self._a * xi - c2 * xi**2
        if self.m != 0:
            q -= self.m**2 / xi2_1

        p_plus = (xi + h / 2)**2 - 1
        p_minus = (xi - h / 2)**2 - 1

        diag = -(p_plus + p_minus) / h**2 + q
        off = p_plus[:-1] / h**2

        if self.m == 0:
            diag[0] = -p_plus[0] / h**2 + q[0]

        from scipy.linalg import eigh_tridiagonal as _eigh_tri
        evals, evecs = _eigh_tri(diag, off)
        idx_sorted = np.argsort(evals)[::-1]
        F_xi = evecs[:, idx_sorted[self.n_radial]]
        if F_xi[0] < 0:
            F_xi = -F_xi
        F_xi /= np.sqrt(np.sum(F_xi**2) * h)
        return xi, F_xi

    def _radial_wavefunction_spectral(
        self, c2: float, A: float
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Compute radial wavefunction F(xi) from spectral Laguerre basis.

        Evaluates the eigenvector on a uniform xi grid for visualization
        and downstream use. Returns (xi_grid, F_xi) tuple.
        """
        N = self.n_basis
        alpha = max(np.sqrt(max(c2, 0.01)), 0.5)
        two_alpha = 2.0 * alpha

        if self.matrix_method == 'algebraic':
            H_mat, S = _build_laguerre_matrices_algebraic(
                N, alpha, A, self._a, c2, self.m
            )
        else:
            H_mat, S = self._build_laguerre_matrices_quadrature(
                N, alpha, A, c2
            )

        evals, evecs = eigh(H_mat, S)
        idx_sorted = np.argsort(evals)[::-1]
        coeffs = evecs[:, idx_sorted[self.n_radial]]

        # Evaluate on a uniform xi grid for output
        n_out = 500
        xi_grid = np.linspace(1.0 + 1e-4, self.xi_max, n_out)
        x_out = two_alpha * (xi_grid - 1.0)
        L_out = np.zeros((N, n_out))
        for n in range(N):
            L_out[n] = eval_laguerre(n, x_out)

        # F(xi) = sum_n c_n * exp(-alpha*(xi-1)) * L_n(2*alpha*(xi-1))
        F_xi = np.exp(-alpha * (xi_grid - 1.0)) * (coeffs @ L_out)

        # Normalize: positive at xi=1, unit norm
        if F_xi[0] < 0:
            F_xi = -F_xi
        dx = xi_grid[1] - xi_grid[0]
        norm = np.sqrt(np.sum(F_xi**2) * dx)
        if norm > 0:
            F_xi /= norm

        return xi_grid, F_xi

    def solve_with_wavefunction(self) -> Dict:
        """Solve and return energy + wavefunction on (xi, eta) grid.

        Returns dict with keys:
            E_elec, c2, A, xi_grid, eta_grid, F_xi, G_eta,
            radial_method (str indicating which solver produced F_xi)
        where F_xi is the radial wavefunction and G_eta is the angular.
        """
        E_elec, c2, A = self.solve()
        c = np.sqrt(max(c2, 1e-15))

        # --- Radial wavefunction F(xi) ---
        if self.radial_method == 'spectral':
            xi, F_xi = self._radial_wavefunction_spectral(c2, A)
        else:
            xi, F_xi = self._radial_wavefunction_fd(c2, A)

        # --- Angular wavefunction G(eta) ---
        n_eta = 200
        eta = np.linspace(-1 + 1e-6, 1 - 1e-6, n_eta)
        # Evaluate angular function using Legendre basis
        from geovac.molecular_sturmian import _angular_sep_const
        n_basis = 50
        from math import factorial
        m_abs = abs(self.m)
        r_vals = np.arange(m_abs, m_abs + n_basis, dtype=float)
        norms = np.array(
            [2.0 / (2*r+1) * factorial(int(r+m_abs)) / factorial(int(r-m_abs))
             for r in r_vals]
        )

        nu_mat = np.zeros((n_basis, n_basis))
        for i in range(n_basis - 1):
            r = r_vals[i]
            val = ((r - m_abs + 1) * np.sqrt(norms[i+1])
                   / ((2*r+1) * np.sqrt(norms[i])))
            nu_mat[i+1, i] = val
            nu_mat[i, i+1] = val

        H_ang = np.diag(-r_vals * (r_vals + 1))
        H_ang += c**2 * (nu_mat @ nu_mat)
        H_ang += self._b * nu_mat

        evals_ang, evecs_ang = np.linalg.eigh(H_ang)
        # n_angular-th from top
        idx_ang = n_basis - 1 - self.n_angular
        coeffs = evecs_ang[:, idx_ang]

        # Evaluate G(eta) = sum_r c_r * P_r^m(eta) / sqrt(norm_r)
        from scipy.special import lpmv
        G_eta = np.zeros(n_eta)
        for j, r in enumerate(r_vals):
            P_r = lpmv(m_abs, int(r), eta)
            G_eta += coeffs[j] * P_r / np.sqrt(norms[j])

        # Normalize
        d_eta = eta[1] - eta[0]
        norm_G = np.sqrt(np.sum(G_eta**2) * d_eta)
        if norm_G > 0:
            G_eta /= norm_G

        return {
            'E_elec': E_elec,
            'c2': c2,
            'A': A,
            'xi_grid': xi,
            'eta_grid': eta,
            'F_xi': F_xi,
            'G_eta': G_eta,
            'radial_method': self.radial_method,
        }


def scan_h2plus_pes(
    R_values: np.ndarray,
    Z_A: int = 1,
    Z_B: int = 1,
    N_xi: int = 5000,
    xi_max: float = 25.0,
    m: int = 0,
    n_angular: int = 0,
    n_radial: int = 0,
    verbose: bool = True,
    radial_method: str = 'fd',
    n_basis: int = 20,
    matrix_method: str = 'quadrature',
) -> Dict[str, np.ndarray]:
    """Scan H2+ PES over a range of R values.

    Parameters
    ----------
    radial_method : str
        'fd' (default) or 'spectral' for spectral Laguerre radial solver.
    n_basis : int
        Number of Laguerre basis functions (only used when radial_method='spectral').
    matrix_method : str
        'quadrature' (default) or 'algebraic' for algebraic Laguerre moments.
    """
    n = len(R_values)
    E_elec = np.zeros(n)
    E_total = np.zeros(n)

    t0 = time.time()
    for i, R in enumerate(R_values):
        lattice = ProlateSpheroidalLattice(
            R=R, Z_A=Z_A, Z_B=Z_B,
            N_xi=N_xi, xi_max=max(xi_max, R * 3), m=m,
            n_angular=n_angular, n_radial=n_radial,
            radial_method=radial_method, n_basis=n_basis,
            matrix_method=matrix_method,
        )
        try:
            E_el, _, _ = lattice.solve()
            E_elec[i] = E_el
            E_total[i] = E_el + Z_A * Z_B / R
        except ValueError:
            E_elec[i] = np.nan
            E_total[i] = np.nan
        if verbose:
            print(f"  R={R:6.3f}  E_elec={E_elec[i]:10.6f}  "
                  f"E_total={E_total[i]:10.6f}")

    dt = time.time() - t0
    if verbose:
        print(f"  PES scan: {n} points in {dt:.1f}s")

    return {
        'R': R_values,
        'E_elec': E_elec,
        'E_total': E_total,
        'V_NN': Z_A * Z_B / R_values,
    }


def fit_spectroscopic_constants(
    R_values: np.ndarray, E_values: np.ndarray
) -> Dict[str, float]:
    """Fit PES near minimum to extract R_eq, E_min, D_e, k."""
    valid = ~np.isnan(E_values)
    R_v = R_values[valid]
    E_v = E_values[valid]

    idx_min = np.argmin(E_v)

    if idx_min == 0 or idx_min == len(R_v) - 1:
        return {
            'R_eq': R_v[idx_min],
            'E_min': E_v[idx_min],
            'D_e': -0.5 - E_v[idx_min],
            'k': 0.0,
            'boundary': True,
        }

    lo = max(0, idx_min - 3)
    hi = min(len(R_v), idx_min + 4)
    coeffs = np.polyfit(R_v[lo:hi], E_v[lo:hi], 2)
    R_eq = -coeffs[1] / (2 * coeffs[0])
    E_min = np.polyval(coeffs, R_eq)

    return {
        'R_eq': R_eq,
        'E_min': E_min,
        'D_e': -0.5 - E_min,
        'k': 2 * coeffs[0],
        'boundary': False,
    }
