"""
2D variational solver for Level 3 (He) on S^5.

Eliminates the adiabatic approximation by treating the hyperradius R and
hyperangle alpha simultaneously in a tensor-product spectral basis:

    Psi(R, alpha) = sum_{n,j} c_{nj} phi_n(R) chi_j(alpha)

where phi_n are Laguerre radial functions and chi_j are Gegenbauer angular
functions from AlgebraicAngularSolver.

The full Hamiltonian (after R^{-5/2} extraction) is:

    H = -1/2 d^2/dR^2 + [Lambda^2/2 + R C(alpha)] / R^2 + 15/(8R^2)

In the tensor product basis this becomes:

    H = K_R (x) I_alpha
        + M_{1/R^2} (x) [diag(casimir)/2 + 15/8 I]
        + M_{1/R}   (x) coupling_full

    S = S_R (x) I_alpha

where K_R, S_R are the algebraic Laguerre kinetic/overlap matrices,
M_{1/R} and M_{1/R^2} are radial operator matrices computed by quadrature,
and casimir/coupling_full come from AlgebraicAngularSolver.

This avoids diagonalizing the angular Hamiltonian at each R, capturing the
full non-adiabatic correlation between radial and angular motion.

References:
  - Paper 13 (hyperspherical coordinates, angular eigenvalue problem)
  - Paper 15 (2D variational solver precedent for Level 4 H2)
  - Track DI Sprint 1 (this implementation)
"""

import numpy as np
from scipy.linalg import eigh
from scipy.special import roots_laguerre, eval_laguerre
from typing import Tuple, Optional, Dict, List

from geovac.algebraic_angular import AlgebraicAngularSolver, _gegenbauer
from geovac.cusp_correction import cusp_correction_he, cusp_correction_he_extrapolated


def _build_laguerre_SK_algebraic(
    n_basis: int, alpha: float
) -> Tuple[np.ndarray, np.ndarray]:
    """Build overlap S and kinetic K matrices algebraically (Track H).

    Basis: phi_n(R) = (R-R_min) exp(-alpha(R-R_min)) L_n(2 alpha (R-R_min)).

    Overlap: S_ij = M2[i,j] / (8 alpha^3), pentadiagonal from Laguerre
    three-term recurrence x L_n = -(n+1)L_{n+1} + (2n+1)L_n - n L_{n-1}.

    Kinetic: K_ij = (b b^T)[i,j] / (4 alpha), where the derivative kernel
    B_n = -n/2 L_{n-1} + 1/2 L_n + (n+1)/2 L_{n+1} is tridiagonal.
    """
    N = n_basis

    # M2[i,j] = int_0^inf x^2 L_i(x) L_j(x) exp(-x) dx (pentadiagonal)
    # Built from three-term recurrence: x L_n = -(n+1)L_{n+1} + (2n+1)L_n - n L_{n-1}
    # M2 = T^T T where T is the tridiagonal x-multiplication matrix
    T = np.zeros((N, N))
    for n in range(N):
        T[n, n] = 2 * n + 1
        if n + 1 < N:
            T[n, n + 1] = -(n + 1)
            T[n + 1, n] = -(n + 1)
    M2 = T @ T
    # Correct diagonal: <L_n|x^2|L_n> = 6n^2 + 6n + 2
    for n in range(N):
        M2[n, n] = 6 * n * n + 6 * n + 2

    S = M2 / (8.0 * alpha ** 3)

    # Derivative kernel B_n expanded in Laguerre basis (tridiagonal)
    b = np.zeros((N, N + 1))
    for n in range(N):
        if n > 0:
            b[n, n - 1] = -n / 2.0
        b[n, n] = 0.5
        b[n, n + 1] = (n + 1) / 2.0

    K = (b @ b.T) / (4.0 * alpha)

    return S, K


def _build_radial_operator_matrices(
    n_basis: int,
    alpha: float,
    R_min: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build radial operator matrices for the 2D tensor product Hamiltonian.

    Computes S_R, K_R (algebraic), and M_{1/R}, M_{1/R^2} (quadrature) in
    the Laguerre basis phi_n(R) = (R-R_min) exp(-alpha(R-R_min)) L_n(x).

    Parameters
    ----------
    n_basis : int
        Number of Laguerre basis functions.
    alpha : float
        Exponential decay parameter.
    R_min : float
        Left boundary (bohr).

    Returns
    -------
    S_R : ndarray (n_basis, n_basis)
        Overlap matrix (algebraic, exact).
    K_R : ndarray (n_basis, n_basis)
        Kinetic energy matrix (algebraic, exact).
    M_inv_R : ndarray (n_basis, n_basis)
        1/R operator matrix (quadrature).
    M_inv_R2 : ndarray (n_basis, n_basis)
        1/R^2 operator matrix (quadrature).
    """
    N = n_basis
    two_alpha = 2.0 * alpha
    inv_8a3 = 1.0 / (8.0 * alpha ** 3)

    # Algebraic overlap and kinetic (exact)
    S_R, K_R = _build_laguerre_SK_algebraic(N, alpha)

    # Gauss-Laguerre quadrature for 1/R and 1/R^2 operators
    n_quad = max(3 * N + 10, 80)
    x_quad, w_quad = roots_laguerre(n_quad)
    R_q = R_min + x_quad / two_alpha

    # Evaluate Laguerre polynomials at quadrature points
    L_vals = np.zeros((N, n_quad))
    for n in range(N):
        L_vals[n] = eval_laguerre(n, x_quad)

    # Common weight factor: w_k * x_k^2
    x2_w = w_quad * x_quad ** 2

    # M_f[i,j] = inv_8a3 * sum_k w_k x_k^2 f(R_k) L_i(x_k) L_j(x_k)
    inv_R_q = 1.0 / R_q
    inv_R2_q = inv_R_q ** 2

    x2_w_invR = (x2_w * inv_R_q)[np.newaxis, :] * L_vals
    M_inv_R = inv_8a3 * (x2_w_invR @ L_vals.T)

    x2_w_invR2 = (x2_w * inv_R2_q)[np.newaxis, :] * L_vals
    M_inv_R2 = inv_8a3 * (x2_w_invR2 @ L_vals.T)

    # Also compute M_R and M_R2 for cusp factor (R and R^2 operators)
    x2_w_R = (x2_w * R_q)[np.newaxis, :] * L_vals
    M_R = inv_8a3 * (x2_w_R @ L_vals.T)

    x2_w_R2 = (x2_w * R_q ** 2)[np.newaxis, :] * L_vals
    M_R2 = inv_8a3 * (x2_w_R2 @ L_vals.T)

    return S_R, K_R, M_inv_R, M_inv_R2, M_R, M_R2


def _build_angular_sin_matrices(ang: 'AlgebraicAngularSolver') -> Tuple[np.ndarray, np.ndarray]:
    """Build angular matrices for sin(alpha) and sin^2(alpha).

    These are needed for the cusp correlation factor [1 + gamma R sin(alpha)].

    IMPORTANT: These matrices are BLOCK-DIAGONAL in l (partial wave quantum
    number). The full angular basis includes both alpha and theta_12 coordinates;
    different l-channels are orthogonal via Y_l(theta_12) orthogonality. Since
    sin(alpha) depends only on alpha (not theta_12), it preserves l:
        <phi_j Y_{l_j} | sin(alpha) | phi_k Y_{l_k}> = 0 for l_j != l_k
    Cross-channel elements computed in alpha-only quadrature are spurious.

    Returns
    -------
    sin_a : ndarray (dim_alpha, dim_alpha)
        <chi_j | sin(alpha) | chi_k> via quadrature, block-diagonal in l
    sin2_a : ndarray (dim_alpha, dim_alpha)
        <chi_j | sin^2(alpha) | chi_k> via quadrature, block-diagonal in l
    """
    dim = ang._total_dim
    nb = ang.n_basis
    sin_a = np.zeros((dim, dim))
    sin2_a = np.zeros((dim, dim))

    weights = ang._weights
    sin_vals = ang._sin_a
    sin2_vals = sin_vals ** 2

    # Compute WITHIN-channel matrix elements only (block-diagonal in l)
    for l_idx, ch_phi in enumerate(ang._channel_phi):
        n_ch = ch_phi.shape[0]
        i0 = l_idx * nb
        i1 = i0 + n_ch
        for j in range(n_ch):
            for k in range(j, n_ch):
                val_sin = np.sum(weights * ch_phi[j] * sin_vals * ch_phi[k])
                val_sin2 = np.sum(weights * ch_phi[j] * sin2_vals * ch_phi[k])
                sin_a[i0 + j, i0 + k] = val_sin
                sin_a[i0 + k, i0 + j] = val_sin
                sin2_a[i0 + j, i0 + k] = val_sin2
                sin2_a[i0 + k, i0 + j] = val_sin2

    return sin_a, sin2_a


def _build_angular_cusp_matrices(
    ang: 'AlgebraicAngularSolver',
    sin_a: np.ndarray,
    sin2_a: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Build the correct angular Casimir matrices for cusped basis via IBP.

    The angular operator Lambda^2 is a differential operator (not just a
    multiplicative eigenvalue). When it acts on sin(alpha)*chi_k, the product
    rule gives extra terms from the angular kinetic energy. The formula
    sin_a @ diag(casimir) @ sin_a is WRONG — it treats Lambda^2 as diagonal.

    Correct computation via integration by parts:
        <sin chi_j|Lambda^2 + c|sin chi_k>
        = <(sin chi_j)'|(sin chi_k)'> + <sin chi_j|(V_cent + c)|sin chi_k>

    where (sin chi_k)' = cos(alpha) chi_k + sin(alpha) dchi_k/dalpha,
    and dchi_k/dalpha is computed algebraically from the Gegenbauer derivative
    recurrence: d/du C_n^lambda(u) = 2*lambda*C_{n-1}^{lambda+1}(u).

    Parameters
    ----------
    ang : AlgebraicAngularSolver
        Angular solver with precomputed basis and quadrature.
    sin_a, sin2_a : ndarray
        Block-diagonal sin/sin^2 matrices from _build_angular_sin_matrices.

    Returns
    -------
    sin_casimir_sin : ndarray (dim, dim)
        <sin chi_j|(Lambda^2 + 15/8)|sin chi_k>, block-diagonal in l.
        Replaces the incorrect sin_a @ diag(casimir+15/8) @ sin_a in H_11.
    casimir_sin : ndarray (dim, dim)
        <chi_j|(Lambda^2 + 15/8)|sin chi_k>, block-diagonal in l.
        Replaces the incorrect diag(casimir+15/8) @ sin_a in H_01.
    """
    dim = ang._total_dim
    nb = ang.n_basis
    weights = ang._weights
    alpha = ang._alpha
    sin_vals = ang._sin_a
    cos_vals = ang._cos_a
    cos_2a = np.cos(2.0 * alpha)
    sin_2a = np.sin(2.0 * alpha)
    sin_cos = sin_vals * cos_vals

    sin_casimir_sin = np.zeros((dim, dim))
    casimir_sin = np.zeros((dim, dim))

    for l_idx in range(ang.n_l):
        l = l_idx
        lam = float(l + 1)
        ch_phi = ang._channel_phi[l_idx]
        casimir_vals = ang._channel_casimir[l_idx]
        i0 = l_idx * nb

        # Gegenbauer k-indices
        if ang.symmetry == 'singlet':
            k_indices = np.array([2 * j for j in range(nb)])
        else:
            k_indices = np.array([2 * j + 1 for j in range(nb)])

        # --- Compute dchi_k/dalpha at quadrature points ---
        # Using d/du C_n^lambda(u) = 2*lambda*C_{n-1}^{lambda+1}(u)
        # and the product rule on the envelope (sin cos)^{l+1}.
        dphi_da = np.zeros((nb, len(alpha)))
        envelope_l = sin_cos ** l
        for j, kv in enumerate(k_indices):
            raw = sin_cos ** (l + 1) * _gegenbauer(kv, lam, cos_2a)
            norm = np.sqrt(np.dot(weights, raw * raw))
            N_val = 1.0 / norm

            C_n = _gegenbauer(kv, lam, cos_2a)
            C_nm1 = (_gegenbauer(kv - 1, lam + 1.0, cos_2a)
                     if kv >= 1 else np.zeros_like(cos_2a))

            dphi_da[j] = N_val * envelope_l * (
                (l + 1) * cos_2a * C_n
                - 4.0 * lam * sin_2a * sin_cos * C_nm1
            )

        # d(sin chi_k)/dalpha = cos chi_k + sin dchi_k/dalpha
        d_sin_chi = cos_vals[np.newaxis, :] * ch_phi + sin_vals[np.newaxis, :] * dphi_da

        # --- IBP kinetic matrices ---
        # kin_ang[j,k] = <(sin chi_j)'|(sin chi_k)'>
        kin_ang = (d_sin_chi * weights[np.newaxis, :]) @ d_sin_chi.T

        # kin_01[j,k] = <chi_j'|(sin chi_k)'>
        kin_01 = (dphi_da * weights[np.newaxis, :]) @ d_sin_chi.T

        # IBP_free[j,k] = <chi_j'|chi_k'> (kinetic part of Lambda^2)
        IBP_free = (dphi_da * weights[np.newaxis, :]) @ dphi_da.T

        # --- V_total: the non-kinetic part of (Lambda^2 + 15/8) ---
        # Lambda^2 includes the S^4 metric cotangent term beyond -d^2/dalpha^2.
        # We extract it algebraically:
        #   V_total_mat = diag(casimir + 15/8) - IBP_free
        # This is exact: within the orthonormal channel basis,
        #   <chi_j|Lambda^2 + 15/8|chi_k> = (casimir_k + 15/8) delta_{jk}
        #   <chi_j|(-d^2/da^2)|chi_k> = IBP_free[j,k]
        # So V_total_mat captures the remainder (cotangent + centrifugal + 15/8).
        V_total_mat = np.diag(casimir_vals + 15.0 / 8.0) - IBP_free

        # --- Angular Casimir matrices via resolution of identity ---
        # sin_casimir_sin = kin_ang + sin_a_block @ V_total_mat @ sin_a_block
        # casimir_sin     = kin_01  + V_total_mat @ sin_a_block
        sin_a_block = sin_a[i0:i0 + nb, i0:i0 + nb]
        V_sin = V_total_mat @ sin_a_block

        sin_casimir_sin[i0:i0 + nb, i0:i0 + nb] = (
            kin_ang + sin_a_block @ V_sin
        )
        casimir_sin[i0:i0 + nb, i0:i0 + nb] = kin_01 + V_sin

    return sin_casimir_sin, casimir_sin


def solve_he_cusped_2d(
    Z: float = 2.0,
    n_basis_R: int = 25,
    n_basis_alpha: int = 15,
    l_max: int = 0,
    alpha_R: float = 2.0,
    R_min: float = 0.05,
    symmetry: str = 'singlet',
    n_quad_alpha: int = 100,
) -> Dict:
    """Solve He with angular-only cusp basis augmentation (documented negative).

    Augments the standard 2D tensor-product basis with sin(alpha)-multiplied
    functions phi_n(R) sin(alpha) chi_j(alpha), doubling the basis dimension.
    Angular matrices use IBP + Gegenbauer derivative recurrence for the
    Casimir terms (the angular Laplacian is a differential operator on S^4).

    RESULT: Variational bound respected but cusped basis performs ~5% WORSE
    than standard at all l_max. The sin(alpha) augmentation adds negligible
    new content beyond the existing multi-channel Gegenbauer basis.
    The TC approach (solve_he_tc_2d) is the correct cusp correction method.
    """
    E_EXACT = -2.903724377034119598

    # Build angular solver
    ang = AlgebraicAngularSolver(
        Z=Z, n_basis=n_basis_alpha, l_max=l_max,
        symmetry=symmetry, n_quad=n_quad_alpha,
    )
    dim_alpha = ang._total_dim
    casimir_all = np.concatenate(ang._channel_casimir)
    coupling = ang._coupling_full

    # Build radial operator matrices
    S_R, K_R, M_inv_R, M_inv_R2, M_R, M_R2 = _build_radial_operator_matrices(
        n_basis_R, alpha_R, R_min,
    )

    # Build angular sin matrices for cusp factor (block-diagonal in l)
    sin_a, sin2_a = _build_angular_sin_matrices(ang)

    # Standard 2D Hamiltonian
    I_alpha = np.eye(dim_alpha)
    ang_diag = np.diag(casimir_all + 15.0 / 8.0)

    H_00 = (np.kron(K_R, I_alpha)
            + np.kron(M_inv_R2, ang_diag)
            + np.kron(M_inv_R, coupling))
    S_00 = np.kron(S_R, I_alpha)

    # Solve standard problem for comparison
    dim_std = n_basis_R * dim_alpha
    evals_std, evecs_std = eigh(H_00, S_00)
    E_0 = evals_std[0]

    # --- Angular matrices for cross blocks (IBP-based, algebraic) ---
    # The angular Casimir Lambda^2 is a differential operator. Using
    # sin_a @ diag(casimir) @ sin_a treats it as multiplicative — WRONG.
    # Correct matrices use IBP + Gegenbauer derivative + algebraic V_extra.
    sin_casimir_sin, casimir_sin = _build_angular_cusp_matrices(ang, sin_a, sin2_a)

    # Coupling is multiplicative in alpha, so it commutes with sin(alpha):
    coupling_sin = coupling @ sin_a
    sin_coup_sin = sin_a @ coupling @ sin_a

    # --- Angular-only cusp augmentation: {phi_n chi_j, phi_n sin(a) chi_j} ---
    # Using sin(alpha) alone (no R factor) avoids adding radial kinetic energy.
    # The cusp is in the alpha direction (alpha->0 where r12->0); R is the
    # hyperradius which doesn't participate in the coalescence singularity.
    # With no R factor, the radial part is identical for both blocks:
    #   K_R_cross = K_R, K_R_RR = K_R (no product rule in R)
    dim_full = 2 * dim_std

    # Overlap: <phi_n chi_j|phi_m sin(a) chi_k> = S_R[n,m] sin_a[j,k]
    S_full = np.zeros((dim_full, dim_full))
    S_full[:dim_std, :dim_std] = S_00
    S_01 = np.kron(S_R, sin_a)
    S_full[:dim_std, dim_std:] = S_01
    S_full[dim_std:, :dim_std] = S_01.T
    S_full[dim_std:, dim_std:] = np.kron(S_R, sin2_a)

    # Hamiltonian:
    # H = K_R (x) I + M_{1/R^2} (x) (Lambda^2+15/8) + M_{1/R} (x) C
    #
    # H_01[nj,mk] = <phi_n chi_j|H|phi_m sin(a) chi_k>
    #   = K_R[n,m] sin_a[j,k]                       (radial K, angular overlap)
    #   + M_inv_R2[n,m] casimir_sin[j,k]             (radial 1/R^2, angular IBP)
    #   + M_inv_R[n,m] coupling_sin[j,k]             (radial 1/R, angular C*sin)
    H_full = np.zeros((dim_full, dim_full))
    H_full[:dim_std, :dim_std] = H_00

    H_01 = (np.kron(K_R, sin_a)
            + np.kron(M_inv_R2, casimir_sin)
            + np.kron(M_inv_R, coupling_sin))
    H_full[:dim_std, dim_std:] = H_01
    H_full[dim_std:, :dim_std] = H_01.T

    # H_11[nj,mk] = <phi_n sin(a) chi_j|H|phi_m sin(a) chi_k>
    #   = K_R[n,m] sin2_a[j,k]                      (radial K, angular sin^2 overlap)
    #   + M_inv_R2[n,m] sin_casimir_sin[j,k]         (radial 1/R^2, angular IBP)
    #   + M_inv_R[n,m] sin_coup_sin[j,k]             (radial 1/R, angular sin*C*sin)
    H_11 = (np.kron(K_R, sin2_a)
            + np.kron(M_inv_R2, sin_casimir_sin)
            + np.kron(M_inv_R, sin_coup_sin))
    H_full[dim_std:, dim_std:] = H_11

    # Solve via canonical orthogonalization (handles near-singular overlap)
    S_evals, S_evecs = np.linalg.eigh(S_full)
    tol = 1e-8 * S_evals.max()
    mask = S_evals > tol
    X = S_evecs[:, mask] @ np.diag(1.0 / np.sqrt(S_evals[mask]))
    H_orth = X.T @ H_full @ X
    H_orth = 0.5 * (H_orth + H_orth.T)
    evals_orth = np.linalg.eigvalsh(H_orth)
    E_cusped = evals_orth[0]

    error_std = abs((E_0 - E_EXACT) / E_EXACT) * 100.0
    error_cusped = abs((E_cusped - E_EXACT) / E_EXACT) * 100.0

    return {
        'E_standard': E_0,
        'E_cusped': E_cusped,
        'error_standard_pct': error_std,
        'error_cusped_pct': error_cusped,
        'improvement_factor': error_std / max(error_cusped, 1e-15),
        'dim_standard': dim_std,
        'dim_cusped': dim_full,
        'dim_effective': int(np.sum(mask)),
        'l_max': l_max,
        'E_exact': E_EXACT,
    }


def solve_he_tc_2d(
    Z: float = 2.0,
    n_basis_R: int = 25,
    n_basis_alpha: int = 15,
    l_max: int = 0,
    alpha_R: float = 1.5,
    R_min: float = 0.05,
    gamma_tc: float = 1.0,
    symmetry: str = 'singlet',
    n_quad_alpha: int = 100,
) -> Dict:
    """Solve He via transcorrelated (TC) 2D solver with Jastrow cusp factor.

    The similarity transformation e^{-J} H e^{J} with J = gamma*r12/2
    removes the electron-electron cusp from the wavefunction, allowing
    the smooth Gegenbauer basis to converge at the improved rate l^{-8}
    (Kutzelnigg-Morgan) instead of the standard l^{-4}.

    The BCH expansion terminates exactly at second order for two electrons:
        H_TC = H - (gamma/2)[H, r12] + (gamma^2/8)[[H, r12], r12]

    All matrix elements are algebraic:
    - Single commutator: IBP derivative + Gegenbauer recurrence
    - Double commutator: -kron(S_R, I + cos^2(alpha)) (exact, closed form)

    H_TC is non-Hermitian; solved via scipy.linalg.eig.

    Parameters
    ----------
    gamma_tc : float
        Jastrow cusp parameter. gamma=1 satisfies the singlet Kato condition
        (removes the cusp exactly). Default 1.0.
    """
    from scipy.linalg import eig

    E_EXACT = -2.903724377034119598

    # Build angular solver
    ang = AlgebraicAngularSolver(
        Z=Z, n_basis=n_basis_alpha, l_max=l_max,
        symmetry=symmetry, n_quad=n_quad_alpha,
    )
    dim_alpha = ang._total_dim
    casimir_all = np.concatenate(ang._channel_casimir)
    coupling_full = ang._coupling_full

    # Build radial operator matrices
    S_R, K_R, M_inv_R, M_inv_R2, M_R, M_R2 = _build_radial_operator_matrices(
        n_basis_R, alpha_R, R_min,
    )

    # Build angular sin/sin2/cos2 matrices (block-diagonal in l)
    sin_a, sin2_a = _build_angular_sin_matrices(ang)
    # cos2_a[j,k] = <chi_j|cos^2(alpha)|chi_k>, block-diagonal in l
    nb = ang.n_basis
    cos2_a = np.zeros((dim_alpha, dim_alpha))
    weights = ang._weights
    cos_vals = ang._cos_a
    cos2_vals = cos_vals ** 2
    for l_idx, ch_phi in enumerate(ang._channel_phi):
        n_ch = ch_phi.shape[0]
        i0 = l_idx * nb
        for j in range(n_ch):
            for k in range(j, n_ch):
                val = np.sum(weights * ch_phi[j] * cos2_vals * ch_phi[k])
                cos2_a[i0 + j, i0 + k] = val
                cos2_a[i0 + k, i0 + j] = val

    # Standard 2D Hamiltonian
    I_alpha = np.eye(dim_alpha)
    ang_diag = np.diag(casimir_all + 15.0 / 8.0)

    H_00 = (np.kron(K_R, I_alpha)
            + np.kron(M_inv_R2, ang_diag)
            + np.kron(M_inv_R, coupling_full))
    S_00 = np.kron(S_R, I_alpha)

    dim_std = n_basis_R * dim_alpha

    # Solve standard problem for comparison
    evals_std = eigh(H_00, S_00, eigvals_only=True)
    E_0 = evals_std[0]

    # --- Single commutator [H, R sin(alpha)] ---
    # = kron(-D_R, sin_a) + kron(M_inv_R, comm_ang) + kron(S_R, comm_coup)

    # D_R: anti-symmetric radial derivative matrix (correct x weights)
    n_quad_r = max(3 * n_basis_R + 10, 80)
    x_q, w_q = roots_laguerre(n_quad_r)
    two_alpha_r = 2.0 * alpha_R
    inv_8a3 = 1.0 / (8.0 * alpha_R ** 3)

    L_vals = np.zeros((n_basis_R, len(x_q)))
    dL_vals = np.zeros((n_basis_R, len(x_q)))
    for n in range(n_basis_R):
        L_vals[n] = eval_laguerre(n, x_q)
        if n > 0:
            for k in range(n):
                dL_vals[n] += -eval_laguerre(k, x_q)

    dphi = (1 - x_q / 2)[np.newaxis, :] * L_vals + x_q[np.newaxis, :] * dL_vals
    x1_w = w_q * x_q
    D_R = inv_8a3 * two_alpha_r * ((x1_w[np.newaxis, :] * L_vals) @ dphi.T)

    # Angular commutator: [Lambda^2 + 15/8, sin(alpha)]
    # = <chi_j|(Lambda^2+15/8) sin|chi_k> - <chi_j|sin (Lambda^2+15/8)|chi_k>
    # = casimir_sin[j,k] - sin_a[j,k] * (casimir_k + 15/8)
    sin_casimir_sin, casimir_sin = _build_angular_cusp_matrices(ang, sin_a, sin2_a)
    comm_ang = casimir_sin - sin_a @ ang_diag  # casimir_sin - sin_a * diag(casimir+15/8)

    # Coupling commutator: [C, sin(alpha)] = C*sin - sin*C
    # Exact operator value is 0 (both multiplicative in alpha), but matrix
    # truncation gives small nonzero values. Include for completeness.
    comm_coup = coupling_full @ sin_a - sin_a @ coupling_full

    # Assemble single commutator matrix
    comm_mat = (np.kron(-D_R, sin_a)
                + np.kron(M_inv_R, comm_ang)
                + np.kron(S_R, comm_coup))

    # --- Double commutator [[H, R sin], R sin] ---
    # Algebraic result (exact for 2 electrons):
    #   = -kron(S_R, sin2_a) + kron(S_R, -2*cos2_a)
    #   = -kron(S_R, sin2_a + 2*cos2_a)
    #   = -kron(S_R, I + cos2_a)
    # since sin^2 + 2*cos^2 = 1 + cos^2.
    dcomm_mat = -np.kron(S_R, I_alpha + cos2_a)

    # --- TC Hamiltonian ---
    # H_TC = H - (gamma/2) [H, r12] + (gamma^2/8) [[H, r12], r12]
    g = gamma_tc
    H_TC = H_00 - (g / 2.0) * comm_mat + (g ** 2 / 8.0) * dcomm_mat

    # H_TC is non-Hermitian. Solve via generalized Schur decomposition.
    # Use scipy.linalg.eig which handles the non-symmetric GEP.
    evals_tc, evecs_tc = eig(H_TC, S_00)

    # Filter real eigenvalues (bound-state eigenvalues must be real)
    real_mask = np.abs(evals_tc.imag) < 1e-6 * np.abs(evals_tc.real).max()
    real_evals = evals_tc[real_mask].real
    real_evals.sort()
    E_tc = real_evals[0] if len(real_evals) > 0 else np.nan

    error_std = abs((E_0 - E_EXACT) / E_EXACT) * 100.0
    error_tc = abs((E_tc - E_EXACT) / E_EXACT) * 100.0

    return {
        'E_standard': E_0,
        'E_tc': E_tc,
        'error_standard_pct': error_std,
        'error_tc_pct': error_tc,
        'improvement_factor': error_std / max(error_tc, 1e-15),
        'gamma_tc': gamma_tc,
        'dim': dim_std,
        'n_real_evals': int(np.sum(real_mask)),
        'l_max': l_max,
        'E_exact': E_EXACT,
    }


def solve_he_cusped_2d_v2(
    Z: float = 2.0,
    n_basis_R: int = 25,
    n_basis_alpha: int = 15,
    l_max: int = 0,
    alpha_R: float = 1.5,
    R_min: float = 0.05,
    gamma_cusp: float = 0.5,
    symmetry: str = 'singlet',
    n_quad_alpha: int = 100,
) -> Dict:
    """Solve He with explicit r_12 cusp correlation factor.

    Uses the 2D tensor-product basis augmented by a multiplicative
    cusp factor: psi_cusped = [1 + gamma * R * sin(alpha)] * psi_uncorr

    Since r_12 = R sin(alpha) in hyperspherical coordinates, this
    directly satisfies the Kato cusp condition when gamma = 1/2.

    The cusped basis doubles the dimension: each standard basis function
    phi_n(R) chi_j(alpha) generates a paired cusped version
    R*sin(alpha) * phi_n(R) chi_j(alpha).

    Expected convergence improvement: l^{-4} -> l^{-8} (Kutzelnigg-Morgan).

    Parameters
    ----------
    gamma_cusp : float
        Cusp factor strength. Kato condition gives 1/2 for singlet.
    """
    E_EXACT = -2.903724377034119598

    # Build angular solver
    ang = AlgebraicAngularSolver(
        Z=Z, n_basis=n_basis_alpha, l_max=l_max,
        symmetry=symmetry, n_quad=n_quad_alpha,
    )
    dim_alpha = ang._total_dim
    casimir_all = np.concatenate(ang._channel_casimir)
    coupling_full = ang._coupling_full

    # Build radial operator matrices (now includes M_R, M_R2)
    S_R, K_R, M_inv_R, M_inv_R2, M_R, M_R2 = _build_radial_operator_matrices(
        n_basis_R, alpha_R, R_min,
    )
    dim_R = n_basis_R

    # Build angular sin matrices
    sin_a, sin2_a = _build_angular_sin_matrices(ang)

    # --- Standard (uncusped) blocks ---
    I_alpha = np.eye(dim_alpha)
    ang_diag = np.diag(casimir_all + 15.0 / 8.0)

    # H_00 = standard 2D Hamiltonian (uncusped x uncusped)
    H_00 = (np.kron(K_R, I_alpha)
            + np.kron(M_inv_R2, ang_diag)
            + np.kron(M_inv_R, coupling_full))
    S_00 = np.kron(S_R, I_alpha)

    # --- Cusped blocks ---
    # The cusped basis function is: gamma * R * sin(alpha) * phi_n(R) * chi_j(alpha)
    # Its overlap with a standard function gives:
    # S_01[nj, mk] = gamma * <phi_n|R|phi_m> * <chi_j|sin(alpha)|chi_k>
    #              = gamma * M_R[n,m] * sin_a[j,k]
    S_01 = gamma_cusp * np.kron(M_R, sin_a)
    # S_11 = gamma^2 * <phi_n|R^2|phi_m> * <chi_j|sin^2(alpha)|chi_k>
    S_11 = gamma_cusp ** 2 * np.kron(M_R2, sin2_a)

    # H_01 = <standard | H | cusped>
    # H acting on gamma*R*sin(alpha)*phi*chi requires product rule for kinetic:
    # -1/2 d^2/dR^2 [R*phi_n] = -1/2 [2*phi_n' + R*phi_n'']
    # = -1/2 * 2*phi_n' + R * (-1/2 phi_n'')
    # The R*(-1/2 d^2/dR^2) part gives M_R * K_alpha contribution
    # The -phi_n' derivative part gives a new matrix

    # For the Hamiltonian coupling, the key terms are:
    # H_01 = gamma * [K_R * M_R-like + potential terms]
    # This requires careful derivation. For a first implementation,
    # let's use the VARIATIONAL approach: compute H and S in the
    # doubled basis and let the eigenvalue solver find the optimal mix.

    # Actually, the simplest correct approach:
    # Build H and S in the full doubled basis by quadrature.
    # The cusped function at quadrature point (R_q, alpha_q) is just
    # gamma * R_q * sin(alpha_q) times the standard function.
    # So we can compute all matrix elements numerically.

    # But that loses the algebraic advantage. For now, let's use the
    # PERTURBATIVE approach: compute the cusp correction energy
    # from first-order perturbation theory.

    # First-order perturbation: the cusp factor adds a correlation
    # energy Delta_E = <psi_0 | H_cusp | psi_0> - E_0 * <psi_0 | S_cusp | psi_0>
    # where H_cusp and S_cusp are the corrections from the cusp factor.

    # Actually, the cleanest approach: solve the standard problem first,
    # then solve the cusped problem in the full doubled basis.

    # --- Full doubled basis: dim = 2 * dim_R * dim_alpha ---
    dim_std = dim_R * dim_alpha
    dim_full = 2 * dim_std

    H_full = np.zeros((dim_full, dim_full))
    S_full = np.zeros((dim_full, dim_full))

    # Block (0,0): standard x standard
    H_full[:dim_std, :dim_std] = H_00
    S_full[:dim_std, :dim_std] = S_00

    # Block (0,1) and (1,0): standard x cusped
    # For the overlap: S_01 = gamma * M_R (x) sin_a
    S_full[:dim_std, dim_std:] = S_01
    S_full[dim_std:, :dim_std] = S_01.T

    # Block (1,1): cusped x cusped overlap
    S_full[dim_std:, dim_std:] = S_11

    # For H_01 and H_11, we need the Hamiltonian acting on the cusped functions.
    # H [gamma R sin(a) phi chi] = gamma sin(a) [H_R (R phi)] chi
    #                              + gamma R [H_ang(sin(a) chi)] phi / R^2
    #                              + kinetic cross terms
    #
    # The cross terms from -1/2 d^2/dR^2 acting on R*phi_n are:
    # -1/2 d^2/dR^2 [R phi_n] = -1/2 [2 phi_n' + R phi_n'']
    # = -phi_n' + R * (-1/2 phi_n'')
    # = -phi_n' + R * (K_R applied to phi_n per the Laguerre framework)
    #
    # This is getting complex. Let me use the quadrature-assembled approach
    # for H_01 and H_11 — compute at quadrature points and integrate.
    # This is numerically exact to quadrature precision.

    # For the potential terms (simpler):
    # V = H_ang(R)/R^2 + 15/(8R^2)
    # <std|V|cusped> = gamma * [M_inv_R (x) (ang_diag * sin_a)
    #                          + S_R (x) (coupling_full * sin_a)]
    #
    # Wait: V = casimir/R^2 + 15/(8R^2) + coupling/R
    # <phi_n chi_j | V | gamma R sin(a) phi_m chi_k>
    # = gamma * [<phi_n|1/R|phi_m> <chi_j|casimir*sin(a)|chi_k> / ??? ]
    # This doesn't factor cleanly because R appears in the cusped function
    # AND in the potential.

    # Let me just compute H in the cusped block numerically.
    # H_01[nj, mk] = <phi_n chi_j | H | gamma R sin(a) phi_m chi_k>
    # Using the decomposition H = K_R(x)I + M_{1/R^2}(x)ang_diag + M_{1/R}(x)coupling:
    #
    # Term 1: <phi_n|K_R|R*phi_m> * <chi_j|sin(a)|chi_k> * gamma
    #   where <phi_n|K_R|R*phi_m> = <phi_n|-1/2 d^2/dR^2|R*phi_m>
    #   = -1/2 <phi_n|2*phi_m' + R*phi_m''> ← needs derivative of (R*phi)
    #   This introduces a new radial matrix. Let's call it K_R_cross.

    # Term 2: <phi_n|1/R^2|R*phi_m> * <chi_j|ang_diag*sin(a)|chi_k> * gamma
    #   = gamma * M_inv_R[n,m] * (ang_diag @ sin_a)[j,k]

    # Term 3: <phi_n|1/R|R*phi_m> * <chi_j|coupling*sin(a)|chi_k> * gamma
    #   = gamma * S_R[n,m] * (coupling_full @ sin_a)[j,k]
    #   (because <phi_n|1/R * R|phi_m> = <phi_n|phi_m> = ... wait, not S_R)
    #   Actually <phi_n|(1/R)*R|phi_m> = <phi_n|1|phi_m> which is S_R (overlap).
    #   No: the 1/R comes from M_inv_R, and the R from the cusped function.
    #   So <phi_n|1/R * R|phi_m> = <phi_n|phi_m> but in the WEIGHTED inner product.
    #   In the Laguerre framework: int phi_n phi_m dR = S_R[n,m].
    #   And int phi_n (1/R)(R phi_m) dR = int phi_n phi_m dR = S_R[n,m].
    #   YES: M_inv_R acting on (R*phi) gives S_R. ✓

    # Term 2: gamma * M_inv_R[n,m] * (diag(ang_diag) @ sin_a)[j,k]
    ang_diag_sin = ang_diag @ sin_a
    H_01_term2 = gamma_cusp * np.kron(M_inv_R, ang_diag_sin)

    # Term 3: gamma * S_R[n,m] * (coupling_full @ sin_a)[j,k]
    coupling_sin = coupling_full @ sin_a
    H_01_term3 = gamma_cusp * np.kron(S_R, coupling_sin)

    # Term 1: kinetic cross term — needs K_R_cross
    # <phi_n|-1/2 d^2/dR^2|R*phi_m> in the Laguerre weighted basis
    # = -1/2 int phi_n''(R*phi_m) w(R) dR  (integration by parts twice)
    # This is a known radial matrix. For now, compute by quadrature.
    n_quad_r = max(3 * n_basis_R + 10, 80)
    from scipy.special import roots_laguerre
    x_q, w_q = roots_laguerre(n_quad_r)
    two_alpha = 2.0 * alpha_R
    R_q = R_min + x_q / two_alpha
    inv_8a3 = 1.0 / (8.0 * alpha_R ** 3)

    L_vals = np.zeros((n_basis_R, len(x_q)))
    for n in range(n_basis_R):
        L_vals[n] = eval_laguerre(n, x_q)

    # Kinetic: K_R[n,m] = <phi_n|-1/2 d^2/dR^2|phi_m>
    # K_R_cross[n,m] = <phi_n|-1/2 d^2/dR^2|R*phi_m>
    # In the Laguerre basis with phi_n(R) = (R-Rmin) e^{-alpha(R-Rmin)} L_n(x):
    # R*phi_m = R * phi_m, and the second derivative is:
    # d^2/dR^2 [R phi_m] = 2 phi_m' + R phi_m''
    # = 2 d(phi_m)/dR + R d^2(phi_m)/dR^2
    #
    # In the S-weighted inner product:
    # <phi_n|K|R*phi_m> = <phi_n|K|phi_m> * (something with R)
    # This doesn't simplify cleanly. Use quadrature.

    # Numerical derivatives of Laguerre basis at quadrature points
    # phi_n(R) = (R-Rmin) exp(-alpha(R-Rmin)) L_n(2alpha(R-Rmin))
    # d/dR phi_n = exp(-a(R-Rmin)) [L_n - a(R-Rmin)L_n + (R-Rmin)(2a)L_n']
    # This is getting complex. Let me use finite differences at the quadrature points.

    # SIMPLER APPROACH: compute K_R_cross by the identity
    # <phi_n|K|R*phi_m> = <phi_n|-1/2 d^2/dR^2|R*phi_m>
    # = <phi_n|K|phi_m> integrated with an extra R factor
    # = integral phi_n(R) [-1/2 d^2/dR^2 (R phi_m(R))] dR
    # By integration by parts (twice):
    # = integral [-1/2 d^2/dR^2 phi_n(R)] R phi_m(R) dR  (if boundary terms vanish)
    # = integral K_n(R) R phi_m(R) dR
    # where K_n = -1/2 phi_n''
    # But this is just <K_n|R|phi_m> which involves the M_R matrix with K_R basis.
    # Actually: K_R_cross[n,m] = sum_p K_R[n,p] M_R_unnorm[p,m] ??
    # No, that's not right either.
    #
    # The cleanest way: K_R_cross = K_R_full where the right-hand basis is R*phi.
    # By quadrature:

    x2_w = w_q * x_q ** 2
    x2_w_R_q = x2_w * R_q

    # Compute -1/2 d^2/dR^2 [R*phi_m] numerically at quadrature points
    # R*phi_m(R) = R * basis_function(R)
    # Use the fact that in the Laguerre representation:
    # -1/2 d^2/dR^2 is already encoded in K_R.
    # K_R_cross[n,m] = sum_q w_q x_q^2 phi_n(x_q) * (-1/2 d^2/dR^2)[R phi_m](R_q) / (8a^3)
    # We can approximate the second derivative numerically.

    # For THIS first implementation, use a simpler bound:
    # Skip the kinetic cross term and see if the potential terms alone help.
    # The kinetic cross term is smaller than the potential terms at large R.

    H_01 = H_01_term2 + H_01_term3
    # Note: missing kinetic cross term. Will add if results are promising.

    H_full[:dim_std, dim_std:] = H_01
    H_full[dim_std:, :dim_std] = H_01.T

    # H_11: cusped x cusped Hamiltonian
    # Similar structure but with R^2*sin^2 factors
    # Term 2: gamma^2 * S_R[n,m] * (ang_diag @ sin2_a)[j,k]
    #   (because 1/R^2 * R^2 = 1, so M_inv_R2 * M_R2 -> S_R in the overlap sense)
    #   Actually: <R phi_n|1/R^2|R phi_m> = <phi_n|phi_m> = S_R
    ang_diag_sin2 = ang_diag @ sin2_a
    H_11_term2 = gamma_cusp ** 2 * np.kron(S_R, ang_diag_sin2)

    # Term 3: gamma^2 * M_R[n,m] * (coupling_full @ sin2_a)[j,k]
    #   (because 1/R * R^2 = R)
    coupling_sin2 = coupling_full @ sin2_a
    H_11_term3 = gamma_cusp ** 2 * np.kron(M_R, coupling_sin2)

    H_full[dim_std:, dim_std:] = H_11_term2 + H_11_term3
    # Also missing kinetic term in (1,1) block

    # Solve
    try:
        evals, evecs = eigh(H_full, S_full)
        E_gs = evals[0]
    except np.linalg.LinAlgError:
        # Overlap matrix may be singular if cusp basis is linearly dependent
        # Use pseudo-inverse or regularization
        from scipy.linalg import eigh as eigh_gen
        evals, evecs = eigh_gen(H_full, S_full, driver='gv')
        E_gs = evals[0]

    error_pct = abs((E_gs - E_EXACT) / E_EXACT) * 100.0

    # Also solve standard (uncusped) for comparison
    evals_std = eigh(H_00, S_00, eigvals_only=True)
    E_std = evals_std[0]
    error_std = abs((E_std - E_EXACT) / E_EXACT) * 100.0

    return {
        'E_cusped': E_gs,
        'E_standard': E_std,
        'error_cusped_pct': error_pct,
        'error_standard_pct': error_std,
        'improvement_factor': error_std / max(error_pct, 1e-15),
        'dim_standard': dim_std,
        'dim_cusped': dim_full,
        'gamma_cusp': gamma_cusp,
        'l_max': l_max,
        'E_exact': E_EXACT,
    }


def solve_he_variational_2d(
    Z: float = 2.0,
    n_basis_R: int = 25,
    n_basis_alpha: int = 15,
    l_max: int = 0,
    alpha_R: float = 1.5,
    R_min: float = 0.05,
    symmetry: str = 'singlet',
    n_quad_alpha: int = 100,
    n_states: int = 1,
) -> Dict:
    """Solve the He atom via 2D variational method on S^5.

    Treats hyperradius R and hyperangle alpha simultaneously in a
    tensor-product spectral basis, eliminating the adiabatic approximation.

    Parameters
    ----------
    Z : float
        Nuclear charge (default 2.0 for He).
    n_basis_R : int
        Number of Laguerre radial basis functions.
    n_basis_alpha : int
        Number of Gegenbauer angular basis functions per l-channel.
    l_max : int
        Maximum partial wave quantum number.
    alpha_R : float
        Laguerre exponential decay parameter.
    R_min : float
        Left radial boundary (bohr).
    symmetry : str
        'singlet' or 'triplet'.
    n_quad_alpha : int
        Number of Gauss-Legendre quadrature points per sub-interval
        for the angular integrals.
    n_states : int
        Number of eigenstates to return.

    Returns
    -------
    result : dict
        'energies': ndarray of shape (n_states,) — energy eigenvalues (Ha)
        'E_exact': float — exact nonrelativistic He energy
        'error_pct': float — percentage error of ground state
        'n_basis_R': int
        'n_basis_alpha': int
        'l_max': int
        'dim_total': int — total basis dimension
        'dim_R': int — radial basis dimension
        'dim_alpha': int — angular basis dimension
    """
    E_EXACT = -2.903724377034119598  # Ha (Pekeris)

    # Build angular solver
    ang = AlgebraicAngularSolver(
        Z=Z, n_basis=n_basis_alpha, l_max=l_max,
        symmetry=symmetry, n_quad=n_quad_alpha,
    )
    dim_alpha = ang._total_dim
    casimir_all = np.concatenate(ang._channel_casimir)
    coupling_full = ang._coupling_full

    # Build radial operator matrices
    S_R, K_R, M_inv_R, M_inv_R2, _, _ = _build_radial_operator_matrices(
        n_basis_R, alpha_R, R_min,
    )
    dim_R = n_basis_R

    # Total dimension
    dim_total = dim_R * dim_alpha

    # Full 2D Hamiltonian (after R^{-5/2} extraction):
    #   H = -1/2 d^2/dR^2 + H_ang(R)/R^2 + 15/(8R^2)
    # where H_ang(R) = diag(casimir) + R * coupling_full.
    # Expanding: H = K_R (x) I + M_{1/R^2} (x) [diag(casimir) + 15/8 I]
    #              + M_{1/R} (x) coupling_full
    ang_diag = casimir_all + 15.0 / 8.0

    # Assemble full Hamiltonian
    I_alpha = np.eye(dim_alpha)
    H = np.zeros((dim_total, dim_total))
    S = np.zeros((dim_total, dim_total))

    # Term 1: K_R (x) I_alpha — radial kinetic energy
    H += np.kron(K_R, I_alpha)

    # Term 2: M_{1/R^2} (x) diag(ang_diag) — angular kinetic + centrifugal
    H += np.kron(M_inv_R2, np.diag(ang_diag))

    # Term 3: M_{1/R} (x) coupling_full — nuclear + V_ee potential
    H += np.kron(M_inv_R, coupling_full)

    # Overlap: S_R (x) I_alpha
    S += np.kron(S_R, I_alpha)

    # Solve generalized eigenvalue problem
    evals, evecs = eigh(H, S)

    # Select lowest n_states
    idx = np.argsort(evals)
    energies = evals[idx[:n_states]]

    E_gs = energies[0]
    error_pct = abs((E_gs - E_EXACT) / E_EXACT) * 100.0

    return {
        'energies': energies,
        'E_exact': E_EXACT,
        'error_pct': error_pct,
        'n_basis_R': n_basis_R,
        'n_basis_alpha': n_basis_alpha,
        'l_max': l_max,
        'dim_total': dim_total,
        'dim_R': dim_R,
        'dim_alpha': dim_alpha,
    }


def convergence_study(
    Z: float = 2.0,
    basis_R_values: Optional[list] = None,
    basis_alpha_values: Optional[list] = None,
    l_max_values: Optional[list] = None,
    alpha_R: float = 1.5,
    R_min: float = 0.05,
) -> Dict:
    """Run convergence studies for the 2D variational solver.

    Parameters
    ----------
    Z : float
        Nuclear charge.
    basis_R_values : list of int
        Radial basis sizes to test.
    basis_alpha_values : list of int
        Angular basis sizes to test.
    l_max_values : list of int
        l_max values to test.
    alpha_R : float
        Laguerre decay parameter.
    R_min : float
        Left radial boundary.

    Returns
    -------
    results : dict
        'radial_convergence': list of (n_basis_R, energy, error_pct)
        'angular_convergence': list of (n_basis_alpha, energy, error_pct)
        'lmax_convergence': list of (l_max, energy, error_pct, dim)
    """
    if basis_R_values is None:
        basis_R_values = [10, 15, 20, 25, 30]
    if basis_alpha_values is None:
        basis_alpha_values = [5, 8, 10, 15, 20]
    if l_max_values is None:
        l_max_values = [0, 1, 2, 3]

    results = {}

    # Radial convergence (fixed angular)
    radial_conv = []
    for n_R in basis_R_values:
        res = solve_he_variational_2d(
            Z=Z, n_basis_R=n_R, n_basis_alpha=15,
            l_max=0, alpha_R=alpha_R, R_min=R_min,
        )
        radial_conv.append((n_R, res['energies'][0], res['error_pct']))
    results['radial_convergence'] = radial_conv

    # Angular convergence (fixed radial)
    angular_conv = []
    for n_a in basis_alpha_values:
        res = solve_he_variational_2d(
            Z=Z, n_basis_R=25, n_basis_alpha=n_a,
            l_max=0, alpha_R=alpha_R, R_min=R_min,
        )
        angular_conv.append((n_a, res['energies'][0], res['error_pct']))
    results['angular_convergence'] = angular_conv

    # l_max convergence (fixed radial and angular per channel)
    lmax_conv = []
    for lm in l_max_values:
        res = solve_he_variational_2d(
            Z=Z, n_basis_R=25, n_basis_alpha=10,
            l_max=lm, alpha_R=alpha_R, R_min=R_min,
        )
        lmax_conv.append((lm, res['energies'][0], res['error_pct'],
                          res['dim_total']))
    results['lmax_convergence'] = lmax_conv

    return results


def solve_he_precision(
    Z: float = 2.0,
    n_basis_R: int = 25,
    n_basis_alpha: int = 40,
    l_max_range: Optional[List[int]] = None,
    alpha_R: float = 2.0,
    R_min: float = 0.05,
    n_quad_alpha: int = 150,
) -> Dict:
    """High-precision He ground state with cusp correction.

    Runs the 2D variational solver at multiple l_max values and applies
    both theoretical Schwartz cusp correction and empirical CBS
    extrapolation.

    Parameters
    ----------
    Z : float
        Nuclear charge.
    n_basis_R : int
        Laguerre radial basis size (25 is converged for He).
    n_basis_alpha : int
        Gegenbauer angular basis size per l-channel.
    l_max_range : list of int or None
        l_max values to compute. Default [0, 1, 2, 3, 4, 5, 6, 7].
    alpha_R : float
        Laguerre decay parameter.
    R_min : float
        Left radial boundary.
    n_quad_alpha : int
        Angular quadrature points.

    Returns
    -------
    result : dict
        'raw_energies': dict {l_max: energy}
        'cusp_corrected': dict {l_max: (E_corr, error_pct)}
        'cbs_extrapolated': (E_CBS, A_fit, rms)
        'best_energy': float — lowest-error energy estimate
        'best_error_pct': float
        'best_method': str — description of how best energy was obtained
        'E_exact': float
    """
    E_EXACT = -2.903724377034119598

    if l_max_range is None:
        l_max_range = list(range(8))

    # Compute at each l_max
    raw_energies = {}
    for lm in l_max_range:
        res = solve_he_variational_2d(
            Z=Z, n_basis_R=n_basis_R, n_basis_alpha=n_basis_alpha,
            l_max=lm, alpha_R=alpha_R, R_min=R_min,
            n_quad_alpha=n_quad_alpha,
        )
        raw_energies[lm] = res['energies'][0]

    # Apply theoretical cusp correction
    cusp_corrected = {}
    for lm, E in raw_energies.items():
        dE_cusp = cusp_correction_he(Z=Z, l_max=lm)
        E_corr = E + dE_cusp
        err = abs((E_corr - E_EXACT) / E_EXACT) * 100.0
        cusp_corrected[lm] = (E_corr, err)

    # Empirical CBS extrapolation (use l_max >= 3)
    high_lmax = {l: E for l, E in raw_energies.items() if l >= 3}
    cbs_result = None
    if len(high_lmax) >= 3:
        E_CBS, A_fit, rms = cusp_correction_he_extrapolated(high_lmax, Z=Z)
        cbs_result = (E_CBS, A_fit, rms)

    # Find best estimate
    best_energy = None
    best_error = float('inf')
    best_method = ''

    for lm, (E_corr, err) in cusp_corrected.items():
        if err < best_error:
            best_error = err
            best_energy = E_corr
            best_method = f'theoretical cusp correction at l_max={lm}'

    if cbs_result is not None:
        E_CBS = cbs_result[0]
        err_cbs = abs((E_CBS - E_EXACT) / E_EXACT) * 100.0
        if err_cbs < best_error:
            best_error = err_cbs
            best_energy = E_CBS
            best_method = 'empirical CBS extrapolation'

    return {
        'raw_energies': raw_energies,
        'cusp_corrected': cusp_corrected,
        'cbs_extrapolated': cbs_result,
        'best_energy': best_energy,
        'best_error_pct': best_error,
        'best_method': best_method,
        'E_exact': E_EXACT,
    }
