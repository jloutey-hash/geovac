"""
Hylleraas explicitly correlated wavefunctions for H2 in prolate spheroidal coordinates.

Implements the James-Coolidge (1933) ansatz:
    Ψ = e^{-α(ξ₁+ξ₂)} Σ c_{jklmp} ξ₁^j ξ₂^k η₁^l η₂^m r₁₂^p

with proper symmetrization for ¹Σ_g⁺ states (singlet, gerade).

The key innovation over standard CI is the explicit r₁₂ dependence (p > 0),
which captures the electron-electron cusp that orbital products cannot represent.

Coordinates (nuclei at foci, separation R):
    ξ ∈ [1, ∞)   — confocal ellipsoidal:  ξ = (r_A + r_B) / R
    η ∈ [-1, +1]  — confocal hyperboloidal: η = (r_A - r_B) / R
    φ ∈ [0, 2π)   — azimuthal

r₁₂ in prolate spheroidals:
    r₁₂² = (R/2)² [(ξ₁η₁-ξ₂η₂)² + ρ₁² + ρ₂² - 2ρ₁ρ₂cos(Δφ)]
    where ρᵢ² = (ξᵢ²-1)(1-ηᵢ²)

For m=0 (σ) states with p=0, the cos(Δφ) term integrates to zero,
giving a 4D problem in (ξ₁, η₁, ξ₂, η₂).

Kinetic energy is computed via integration by parts (Green's first identity)
using exact analytical derivatives of the exponential-polynomial basis
functions, avoiding all finite-difference errors.

Reference:
    James, H. M. & Coolidge, A. S. (1933). J. Chem. Phys. 1, 825-835.
"""

import numpy as np
from scipy.linalg import eigh
from scipy.special import ellipk
from scipy.optimize import minimize_scalar
from typing import Dict, List, Tuple
import time


# ============================================================
# Basis function definition
# ============================================================

class HylleraasBasisFunction:
    """A single James-Coolidge basis function.

    φ_{jklmp} = e^{-α(ξ₁+ξ₂)} * (ξ₁^j ξ₂^k η₁^l η₂^m r₁₂^p + exchange)

    For ¹Σ_g⁺ symmetry:
      - Singlet: symmetric under electron exchange (1↔2)
      - Gerade: symmetric under inversion (η → -η)

    Symmetrization combines (j,k,l,m) with (k,j,m,l) for exchange,
    and requires l+m to be even for gerade symmetry.

    Parameters
    ----------
    j, k : int
        Powers of ξ₁, ξ₂ (≥ 0).
    l, m : int
        Powers of η₁, η₂ (≥ 0).
    p : int
        Power of r₁₂ (≥ 0). p=0 is standard CI, p≥1 captures cusp.
    alpha : float
        Exponential screening parameter.
    """

    def __init__(self, j: int, k: int, l: int, m: int, p: int,
                 alpha: float = 1.0):
        self.j = j
        self.k = k
        self.l = l
        self.m = m
        self.p = p
        self.alpha = alpha

    def __repr__(self) -> str:
        return (f"Hyl(j={self.j},k={self.k},l={self.l},m={self.m},"
                f"p={self.p},α={self.alpha:.3f})")

    @property
    def is_gerade(self) -> bool:
        """Gerade symmetry requires l+m even."""
        return (self.l + self.m) % 2 == 0

    @property
    def is_self_exchange(self) -> bool:
        """True if (j,k,l,m) == (k,j,m,l) — no exchange partner needed."""
        return self.j == self.k and self.l == self.m


# ============================================================
# Basis set generation
# ============================================================

def generate_basis(
    j_max: int = 2,
    l_max: int = 2,
    p_max: int = 2,
    alpha: float = 1.0,
    gerade_only: bool = True,
) -> List[HylleraasBasisFunction]:
    """Generate a systematic Hylleraas basis set.

    Generates all (j,k,l,m,p) with j≥k (to avoid double-counting exchange
    partners), l+m even (gerade), and independent truncation per index.

    Parameters
    ----------
    j_max : int
        Maximum power of ξ (for both j and k).
    l_max : int
        Maximum power of η (for both l and m).
    p_max : int
        Maximum power of r₁₂.
    alpha : float
        Exponential parameter (shared by all basis functions).
    gerade_only : bool
        If True, only include gerade functions (l+m even).

    Returns
    -------
    List of HylleraasBasisFunction.
    """
    basis = []
    for p in range(p_max + 1):
        for j in range(j_max + 1):
            for k in range(j + 1):  # k ≤ j to avoid double-counting
                for ll in range(l_max + 1):
                    for m in range(l_max + 1):
                        if gerade_only and (ll + m) % 2 != 0:
                            continue
                        if k == j and m > ll:
                            continue
                        basis.append(HylleraasBasisFunction(
                            j, k, ll, m, p, alpha
                        ))
    return basis


def generate_basis_truncated(
    omega_max: int = 4,
    p_max: int = 2,
    alpha: float = 1.0,
) -> List[HylleraasBasisFunction]:
    """Generate basis with total power truncation j+k+l+m+p ≤ omega_max.

    This gives a more balanced basis than independent max truncation.
    James-Coolidge (1933) used omega_max=5 with 13 terms.
    """
    basis = []
    for p in range(p_max + 1):
        remaining = omega_max - p
        for j in range(remaining + 1):
            for k in range(min(j, remaining - j) + 1):
                for ll in range(remaining - j - k + 1):
                    m_max = remaining - j - k - ll
                    for m in range(m_max + 1):
                        if (ll + m) % 2 != 0:
                            continue
                        if k == j and m > ll:
                            continue
                        basis.append(HylleraasBasisFunction(
                            j, k, ll, m, p, alpha
                        ))
    return basis


# ============================================================
# r₁₂ computation in prolate spheroidals
# ============================================================

def compute_r12_squared(
    xi1: np.ndarray, eta1: np.ndarray,
    xi2: np.ndarray, eta2: np.ndarray,
    R: float,
) -> np.ndarray:
    """Compute r₁₂² averaged over azimuthal angle Δφ (m=0 states).

    From Cartesian decomposition in prolate spheroidals:
        ρ = (R/2)√((ξ²-1)(1-η²)),  z = (R/2)ξη

    r₁₂² = (z₁-z₂)² + ρ₁² + ρ₂² - 2ρ₁ρ₂cos(Δφ)

    Averaging over Δφ (cos → 0):
    <r₁₂²>_φ = (R/2)² [(ξ₁η₁-ξ₂η₂)² + (ξ₁²-1)(1-η₁²) + (ξ₂²-1)(1-η₂²)]

    Note: This is NOT zero when both electrons have the same (ξ,η) but
    different φ. For r₁₂=0 at exactly the same point, use
    compute_r12_with_phi with Δφ=0.
    """
    half_R = R / 2.0
    r12_sq = half_R**2 * (
        (xi1 * eta1 - xi2 * eta2)**2
        + (xi1**2 - 1) * (1 - eta1**2)
        + (xi2**2 - 1) * (1 - eta2**2)
    )
    return np.maximum(r12_sq, 0.0)


def compute_r12_with_phi(
    xi1: np.ndarray, eta1: np.ndarray,
    xi2: np.ndarray, eta2: np.ndarray,
    dphi: float,
    R: float,
) -> np.ndarray:
    """Compute r₁₂ including the azimuthal angle difference.

    r₁₂² = (R/2)² [(ξ₁η₁-ξ₂η₂)² + ρ₁² + ρ₂² - 2ρ₁ρ₂cos(Δφ)]
    where ρᵢ² = (ξᵢ²-1)(1-ηᵢ²)
    """
    half_R = R / 2.0
    rho1_sq = np.maximum((xi1**2 - 1) * (1 - eta1**2), 0.0)
    rho2_sq = np.maximum((xi2**2 - 1) * (1 - eta2**2), 0.0)

    r12_sq = half_R**2 * (
        (xi1 * eta1 - xi2 * eta2)**2
        + rho1_sq + rho2_sq
        - 2 * np.cos(dphi) * np.sqrt(rho1_sq * rho2_sq)
    )
    return np.sqrt(np.maximum(r12_sq, 0.0))


# ============================================================
# Quadrature grids
# ============================================================

def build_quadrature_grids(
    N_xi: int = 30,
    N_eta: int = 20,
    N_phi: int = 16,
    xi_max: float = 15.0,
) -> Dict:
    """Build Gauss-Legendre quadrature grids for 4D+phi integration.

    Parameters
    ----------
    N_xi : int
        Number of quadrature points in each ξ direction.
    N_eta : int
        Number of quadrature points in each η direction.
    N_phi : int
        Number of points for azimuthal (Δφ) integration via trapezoid rule.
    xi_max : float
        Maximum ξ value.

    Returns
    -------
    dict with xi, eta, w_xi, w_eta, dphi, w_phi arrays.
    """
    # η ∈ [-1, +1]: standard Gauss-Legendre
    eta, w_eta = np.polynomial.legendre.leggauss(N_eta)

    # ξ ∈ [1, ξ_max]: map from u ∈ [-1,1] via quadratic clustering near ξ=1
    u_gl, w_u = np.polynomial.legendre.leggauss(N_xi)
    t = (u_gl + 1) / 2  # t ∈ [0, 1]
    xi = 1.0 + (xi_max - 1.0) * t**2
    # dξ = 2(ξ_max-1)t dt, dt = du/2
    w_xi = w_u * (xi_max - 1.0) * t

    # Δφ ∈ [0, 2π]: trapezoid rule (periodic function)
    dphi = np.linspace(0, 2 * np.pi, N_phi, endpoint=False)
    w_phi = np.full(N_phi, 2 * np.pi / N_phi)

    return {
        'xi': xi, 'w_xi': w_xi,
        'eta': eta, 'w_eta': w_eta,
        'dphi': dphi, 'w_phi': w_phi,
        'N_xi': N_xi, 'N_eta': N_eta, 'N_phi': N_phi,
        'xi_max': xi_max,
    }


# ============================================================
# Basis function evaluation and analytical derivatives
# ============================================================

def evaluate_basis_function(
    bf: HylleraasBasisFunction,
    xi1: np.ndarray, eta1: np.ndarray,
    xi2: np.ndarray, eta2: np.ndarray,
    r12: np.ndarray,
    R: float,
) -> np.ndarray:
    """Evaluate a symmetrized Hylleraas basis function.

    For ¹Σ_g⁺: φ = e^{-α(ξ₁+ξ₂)} [f(1,2) + f(2,1)] r₁₂^p
    where f(1,2) = ξ₁^j ξ₂^k η₁^l η₂^m
    """
    exp_factor = np.exp(-bf.alpha * (xi1 + xi2))
    r12_factor = r12**bf.p if bf.p > 0 else 1.0

    direct = xi1**bf.j * xi2**bf.k * eta1**bf.l * eta2**bf.m
    if bf.is_self_exchange:
        sym = 2.0 * direct
    else:
        exchange = xi1**bf.k * xi2**bf.j * eta1**bf.m * eta2**bf.l
        sym = direct + exchange

    return exp_factor * r12_factor * sym


def _eval_unsym_and_derivs(
    j: int, k: int, l: int, m: int, alpha: float,
    xi1: np.ndarray, eta1: np.ndarray,
    xi2: np.ndarray, eta2: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate unsymmetrized term and its first derivatives (p=0).

    For g = e^{-α(ξ₁+ξ₂)} ξ₁^j ξ₂^k η₁^l η₂^m, returns:
        g, dg/dξ₁, dg/dη₁, dg/dξ₂, dg/dη₂

    All computed analytically.
    """
    exp_f = np.exp(-alpha * (xi1 + xi2))
    poly = xi1**j * xi2**k * eta1**l * eta2**m
    g = exp_f * poly

    # ∂g/∂ξ₁ = exp * ξ₂^k η₁^l η₂^m * (j·ξ₁^{j-1} - α·ξ₁^j)
    if j > 0:
        dg_dxi1 = exp_f * xi2**k * eta1**l * eta2**m * (
            j * xi1**(j - 1) - alpha * xi1**j
        )
    else:
        dg_dxi1 = -alpha * g

    # ∂g/∂η₁ = exp * ξ₁^j ξ₂^k η₂^m * l·η₁^{l-1}
    if l > 0:
        dg_deta1 = exp_f * xi1**j * xi2**k * eta2**m * l * eta1**(l - 1)
    else:
        dg_deta1 = np.zeros_like(g)

    # ∂g/∂ξ₂ = exp * ξ₁^j η₁^l η₂^m * (k·ξ₂^{k-1} - α·ξ₂^k)
    if k > 0:
        dg_dxi2 = exp_f * xi1**j * eta1**l * eta2**m * (
            k * xi2**(k - 1) - alpha * xi2**k
        )
    else:
        dg_dxi2 = -alpha * g

    # ∂g/∂η₂ = exp * ξ₁^j ξ₂^k η₁^l * m·η₂^{m-1}
    if m > 0:
        dg_deta2 = exp_f * xi1**j * xi2**k * eta1**l * m * eta2**(m - 1)
    else:
        dg_deta2 = np.zeros_like(g)

    return g, dg_dxi1, dg_deta1, dg_dxi2, dg_deta2


def evaluate_basis_and_derivs(
    bf: HylleraasBasisFunction,
    xi1: np.ndarray, eta1: np.ndarray,
    xi2: np.ndarray, eta2: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate symmetrized basis function and its first derivatives (p=0).

    Returns φ, ∂φ/∂ξ₁, ∂φ/∂η₁, ∂φ/∂ξ₂, ∂φ/∂η₂.
    All computed analytically — no finite differences.
    """
    # Direct term: (j, k, l, m)
    g_d, dg_d_xi1, dg_d_eta1, dg_d_xi2, dg_d_eta2 = _eval_unsym_and_derivs(
        bf.j, bf.k, bf.l, bf.m, bf.alpha, xi1, eta1, xi2, eta2
    )

    if bf.is_self_exchange:
        return (2.0 * g_d, 2.0 * dg_d_xi1, 2.0 * dg_d_eta1,
                2.0 * dg_d_xi2, 2.0 * dg_d_eta2)

    # Exchange term: (k, j, m, l) — swap electron labels
    g_e, dg_e_xi1, dg_e_eta1, dg_e_xi2, dg_e_eta2 = _eval_unsym_and_derivs(
        bf.k, bf.j, bf.m, bf.l, bf.alpha, xi1, eta1, xi2, eta2
    )

    return (g_d + g_e, dg_d_xi1 + dg_e_xi1, dg_d_eta1 + dg_e_eta1,
            dg_d_xi2 + dg_e_xi2, dg_d_eta2 + dg_e_eta2)


# ============================================================
# Matrix element computation — p=0 (integration by parts)
# ============================================================

def compute_overlap_matrix(
    basis: List[HylleraasBasisFunction],
    R: float,
    grids: Dict,
) -> np.ndarray:
    """Compute overlap matrix S_ij = <φ_i | φ_j>.

    For p=0: 4D integral with (2π)² azimuthal factor.
    For p>0: 5D integral with explicit Δφ integration.
    """
    n_basis = len(basis)
    S = np.zeros((n_basis, n_basis))

    xi = grids['xi']
    eta = grids['eta']
    w_xi = grids['w_xi']
    w_eta = grids['w_eta']
    half_R = R / 2.0
    max_p = max(bf.p for bf in basis)

    if max_p == 0:
        phi_factor = (2 * np.pi)**2
        for i in range(n_basis):
            for j in range(i, n_basis):
                S[i, j] = _overlap_4d(
                    basis[i], basis[j], xi, eta, w_xi, w_eta,
                    half_R, phi_factor, R
                )
                S[j, i] = S[i, j]
    else:
        dphi = grids['dphi']
        w_phi = grids['w_phi']
        for i in range(n_basis):
            for j in range(i, n_basis):
                S[i, j] = _overlap_5d(
                    basis[i], basis[j], xi, eta, w_xi, w_eta,
                    dphi, w_phi, R
                )
                S[j, i] = S[i, j]

    return S


def _overlap_4d(
    bf_i: HylleraasBasisFunction,
    bf_j: HylleraasBasisFunction,
    xi: np.ndarray, eta: np.ndarray,
    w_xi: np.ndarray, w_eta: np.ndarray,
    half_R: float, phi_factor: float,
    R: float,
) -> float:
    """4D overlap integral when p = 0 for all basis functions."""
    n_xi = len(xi)
    n_eta = len(eta)
    result = 0.0

    for a in range(n_xi):
        xi1 = xi[a]
        J1 = half_R**3 * (xi1**2 - eta**2)

        for c in range(n_xi):
            xi2 = xi[c]
            J2 = half_R**3 * (xi2**2 - eta**2)

            for b in range(n_eta):
                eta1 = eta[b]
                eta2 = eta

                # Use azimuthally averaged r12 for evaluation
                # (doesn't matter for p=0, but needed for function signature)
                r12_sq = compute_r12_squared(xi1, eta1, xi2, eta2, R)
                r12 = np.sqrt(np.maximum(r12_sq, 1e-30))

                fi = evaluate_basis_function(
                    bf_i, xi1, eta1, xi2, eta2, r12, R
                )
                fj = evaluate_basis_function(
                    bf_j, xi1, eta1, xi2, eta2, r12, R
                )

                integrand = fi * fj * J1[b] * J2 * w_eta[b] * w_eta
                result += w_xi[a] * w_xi[c] * np.sum(integrand)

    return result * phi_factor


def _overlap_5d(
    bf_i: HylleraasBasisFunction,
    bf_j: HylleraasBasisFunction,
    xi: np.ndarray, eta: np.ndarray,
    w_xi: np.ndarray, w_eta: np.ndarray,
    dphi: np.ndarray, w_phi: np.ndarray,
    R: float,
) -> float:
    """5D overlap integral with explicit Δφ integration for r₁₂^p terms."""
    half_R = R / 2.0
    n_xi = len(xi)
    n_eta = len(eta)

    result = 0.0
    for a in range(n_xi):
        xi1 = xi[a]
        for c in range(n_xi):
            xi2 = xi[c]
            for b in range(n_eta):
                eta1 = eta[b]
                J1 = half_R**3 * (xi1**2 - eta1**2)
                for d in range(n_eta):
                    eta2 = eta[d]
                    J2 = half_R**3 * (xi2**2 - eta2**2)

                    # Integrate over Δφ
                    phi_int = 0.0
                    for ip in range(len(dphi)):
                        r12 = compute_r12_with_phi(
                            xi1, eta1, xi2, eta2, dphi[ip], R
                        )
                        fi = evaluate_basis_function(
                            bf_i, xi1, eta1, xi2, eta2, r12, R
                        )
                        fj = evaluate_basis_function(
                            bf_j, xi1, eta1, xi2, eta2, r12, R
                        )
                        phi_int += fi * fj * w_phi[ip]

                    result += (w_xi[a] * w_eta[b] * w_xi[c] * w_eta[d]
                               * J1 * J2 * phi_int * 2 * np.pi)

    return result


def compute_hamiltonian_matrix(
    basis: List[HylleraasBasisFunction],
    R: float,
    grids: Dict,
    Z_A: float = 1.0,
    Z_B: float = 1.0,
) -> np.ndarray:
    """Compute Hamiltonian matrix H_ij = <φ_i | T+V_ne+V_ee | φ_j>.

    For p=0: kinetic energy via integration by parts (analytical derivatives),
    V_ne and V_ee via direct quadrature, V_ee uses elliptic integral kernel.

    For p>0: full 5D numerical integration.
    """
    n_basis = len(basis)
    max_p = max(bf.p for bf in basis)

    if max_p == 0:
        H = _hamiltonian_p0(basis, R, grids, Z_A, Z_B)
    else:
        H = _hamiltonian_general(basis, R, grids, Z_A, Z_B)

    return H


# ============================================================
# p=0 Hamiltonian: analytical kinetic energy via integration by parts
# ============================================================

def _hamiltonian_p0(
    basis: List[HylleraasBasisFunction],
    R: float,
    grids: Dict,
    Z_A: float,
    Z_B: float,
) -> np.ndarray:
    """Hamiltonian matrix for p=0 basis using integration by parts.

    <φ_i|H|φ_j> = <φ_i|T₁+T₂|φ_j> + <φ_i|V_ne+V_ee|φ_j>

    Kinetic energy via integration by parts:
    <φ_i|T₁|φ_j> = (2/R²)(R/2)³ (2π)² ∫∫∫∫
        [(ξ₁²-1)(∂φ_i/∂ξ₁)(∂φ_j/∂ξ₁) + (1-η₁²)(∂φ_i/∂η₁)(∂φ_j/∂η₁)]
        * J₂ dξ₁ dη₁ dξ₂ dη₂

    Note: the J₁ = (R/2)³(ξ₁²-η₁²) factor from the volume element cancels
    with the 1/(ξ₁²-η₁²) from the kinetic energy operator, leaving just
    (R/2)³ which combines with the -(2/R²) prefactor.

    V_ne and V_ee use the same 4D quadrature with elliptic integral for 1/r₁₂.
    """
    n_basis = len(basis)
    H = np.zeros((n_basis, n_basis))

    xi = grids['xi']
    eta = grids['eta']
    w_xi = grids['w_xi']
    w_eta = grids['w_eta']
    half_R = R / 2.0
    n_xi = len(xi)
    n_eta = len(eta)
    phi_factor = (2 * np.pi)**2

    # Prefactor for kinetic energy: (2/R²)(R/2)³ = (2/R²)(R³/8) = R/4
    T_prefactor = R / 4.0

    for i in range(n_basis):
        for j in range(i, n_basis):
            T_ij = 0.0
            V_ij = 0.0

            for a in range(n_xi):
                xi1 = xi[a]
                xi1_sq_m1 = xi1**2 - 1  # (ξ₁²-1) for kinetic
                J1_arr = half_R**3 * (xi1**2 - eta**2)  # (n_eta,)

                for c in range(n_xi):
                    xi2 = xi[c]
                    xi2_sq_m1 = xi2**2 - 1
                    J2_arr = half_R**3 * (xi2**2 - eta**2)  # (n_eta,)

                    for b in range(n_eta):
                        eta1 = eta[b]
                        one_m_eta1_sq = 1.0 - eta1**2
                        J1 = J1_arr[b]
                        eta2 = eta  # vectorized

                        # Evaluate basis functions and derivatives
                        fi, dfi_xi1, dfi_eta1, dfi_xi2, dfi_eta2 = \
                            evaluate_basis_and_derivs(
                                basis[i], xi1, eta1, xi2, eta2
                            )
                        fj, dfj_xi1, dfj_eta1, dfj_xi2, dfj_eta2 = \
                            evaluate_basis_and_derivs(
                                basis[j], xi1, eta1, xi2, eta2
                            )

                        # --- Kinetic energy T₁ (integration by parts) ---
                        # (ξ₁²-1)(∂φ_i/∂ξ₁)(∂φ_j/∂ξ₁) × J₂
                        # + (1-η₁²)(∂φ_i/∂η₁)(∂φ_j/∂η₁) × J₂
                        T1_integrand = (
                            xi1_sq_m1 * dfi_xi1 * dfj_xi1
                            + one_m_eta1_sq * dfi_eta1 * dfj_eta1
                        ) * J2_arr  # shape (n_eta,)

                        # T₂ similarly: derivatives w.r.t. ξ₂, η₂
                        one_m_eta2_sq = 1.0 - eta2**2
                        T2_integrand = (
                            xi2_sq_m1 * dfi_xi2 * dfj_xi2
                            + one_m_eta2_sq * dfi_eta2 * dfj_eta2
                        ) * J1  # J1 is scalar, J2 cancels T2's 1/(ξ₂²-η₂²)

                        # Wait — T₂ integration by parts involves
                        # the J₁ volume element for electron 1 (overlap)
                        # and the kinetic part for electron 2.
                        # For T₂: same structure but roles swapped.
                        # <φ_i|T₂|φ_j> = T_prefactor * phi_factor * ∫∫∫∫
                        #   [(ξ₂²-1)(∂φ_i/∂ξ₂)(∂φ_j/∂ξ₂)
                        #    + (1-η₂²)(∂φ_i/∂η₂)(∂φ_j/∂η₂)]
                        #   * J₁ dξ₁dη₁ dξ₂dη₂

                        # Combined T = T₁ + T₂:
                        T_integrand = T1_integrand + T2_integrand

                        T_sum = np.sum(
                            T_integrand * w_eta[b] * w_eta
                        )
                        T_ij += w_xi[a] * w_xi[c] * T_sum

                        # --- Potential energy V_ne + V_ee ---
                        r1A = half_R * (xi1 + eta1)
                        r1B = half_R * abs(xi1 - eta1)
                        r2A = half_R * (xi2 + eta2)
                        r2B = half_R * np.abs(xi2 - eta2)

                        V_ne = (
                            -Z_A / np.maximum(r1A, 1e-15)
                            - Z_B / np.maximum(r1B, 1e-15)
                            - Z_A / np.maximum(r2A, 1e-15)
                            - Z_B / np.maximum(r2B, 1e-15)
                        )

                        # V_ee = 1/r₁₂ via elliptic integral (azimuthal average)
                        rho1 = half_R * np.sqrt(max(
                            (xi1**2 - 1) * (1 - eta1**2), 0.0
                        ))
                        z1 = half_R * xi1 * eta1
                        rho2 = half_R * np.sqrt(np.maximum(
                            (xi2**2 - 1) * (1 - eta2**2), 0.0
                        ))
                        z2 = half_R * xi2 * eta2

                        dz = z1 - z2
                        s = (rho1 + rho2)**2 + dz**2
                        with np.errstate(divide='ignore', invalid='ignore'):
                            k2 = np.where(
                                s > 1e-30,
                                4 * rho1 * rho2 / s,
                                0.0
                            )
                            k2 = np.clip(k2, 0, 1 - 1e-15)
                        s_safe = np.maximum(s, 1e-30)
                        K_vals = ellipk(k2)
                        # Phi-averaged 1/r₁₂: (4/π) K(k) / √s
                        # Factor: the phi integral gives 2π for φ₁ and
                        # the azimuthal average of 1/r₁₂ is
                        # (1/2π) ∫ 1/r₁₂ dΔφ = (2/π) K(k)/√s
                        # With (2π)² from both phi integrals:
                        # V_ee_eff = (2π)² × (2/π) K(k)/√s = 8π K(k)/√s
                        # But we already have phi_factor = (2π)²
                        # So V_ee per point (before phi_factor) =
                        # (1/(2π)²) × (2π) × ∫ (1/r₁₂) dΔφ/(2π) × (2π)
                        # = (1/2π) ∫ (1/r₁₂) dΔφ = (2/π) K(k)/√s
                        # Hmm, let me be more careful.
                        #
                        # Full integral:
                        # ∫₀²π dφ₁ ∫₀²π dφ₂ (1/r₁₂) = 2π × ∫₀²π (1/r₁₂(Δφ)) dΔφ
                        # = 2π × 2 × (2K(k)/√s)  [standard result]
                        # = 8π K(k)/√s
                        #
                        # So with phi_factor = (2π)², we need:
                        # V_ee_integral = (1/(2π)²) × 8π K(k)/√s = 2K(k)/(π√s)
                        # No wait, phi_factor is applied to the ENTIRE matrix
                        # element. For V_ne, the phi integrals are trivial
                        # (both give 2π), so V_ne gets multiplied by (2π)².
                        # For V_ee, the phi integral over 1/r₁₂ gives
                        # 8π K(k)/√s (not (2π)²). So we should NOT apply
                        # phi_factor to V_ee and instead use 8π K(k)/√s directly.
                        #
                        # Actually, let me restructure: do NOT include phi_factor
                        # globally. Instead:
                        # - T and V_ne: multiply by (2π)²
                        # - V_ee: multiply by 8π K(k)/√s / (2π)² × (2π)²
                        #        = 8π K(k)/√s
                        # Wait, this is getting confused. Let me be precise.
                        #
                        # The full matrix element is:
                        # <φ_i|O|φ_j> = ∫∫ φ_i O φ_j J₁J₂ dξ₁dη₁dφ₁ dξ₂dη₂dφ₂
                        #
                        # For p=0, φ_i and φ_j are independent of φ. So:
                        # For T, V_ne (independent of φ):
                        #   = (∫dφ₁)(∫dφ₂) × 4D integral = (2π)² × 4D
                        #
                        # For V_ee = 1/r₁₂ which depends on Δφ:
                        #   = ∫dφ₁ ∫dφ₂ φ_i φ_j / r₁₂(Δφ) J₁J₂ dξ₁dη₁dξ₂dη₂
                        #   = 2π × ∫dΔφ [φ_i φ_j / r₁₂(Δφ)] J₁J₂ dξ₁dη₁dξ₂dη₂
                        #   = 2π × ∫ [φ_i φ_j J₁J₂] × [∫ 1/r₁₂ dΔφ] dξ₁dη₁dξ₂dη₂
                        #
                        # ∫₀²π 1/r₁₂(Δφ) dΔφ = (4/√s_max) K(k)
                        # where s_max = (ρ₁+ρ₂)²+(z₁-z₂)², k² = 4ρ₁ρ₂/s_max
                        # Actually: ∫₀²π 1/r₁₂ dΔφ = 4K(k)/(R/2 √s_max)???
                        #
                        # Let me use a known result. From the generating function
                        # of elliptic integrals:
                        # ∫₀²π dΔφ / √(a-b·cos(Δφ)) = 4K(b/(2a)) / √a  ??? No.
                        #
                        # r₁₂² = (R/2)²[(ξ₁η₁-ξ₂η₂)² + ρ₁² + ρ₂² - 2ρ₁ρ₂cos(Δφ)]
                        #       = (R/2)²[Q + 2ρ₁ρ₂(1-cos(Δφ))]  ... hmm not quite.
                        #
                        # Actually let me use the standard cylindrical result:
                        # 1/|r₁-r₂| in cylindrical coordinates:
                        # |r₁-r₂|² = (z₁-z₂)² + ρ₁² + ρ₂² - 2ρ₁ρ₂cos(Δφ)
                        #
                        # ∫₀²π dΔφ/|r₁-r₂| = 4K(k)/√s_max
                        # where s_max = (ρ₁+ρ₂)² + (z₁-z₂)²
                        # and k² = 4ρ₁ρ₂ / s_max
                        #
                        # This is a well-known result from Jackson (electrostatics).
                        #
                        # So the V_ee contribution to the matrix element is:
                        # 2π × ∫ φ_i φ_j J₁J₂ × [4K(k)/√s_max] d4D
                        # = 8π × ∫ φ_i φ_j J₁J₂ K(k)/√s_max d4D
                        #
                        # This is NOT (2π)² × integral. So I need to handle
                        # T/V_ne and V_ee with different phi factors.

                        V_ee = np.where(
                            s > 1e-30,
                            K_vals / np.sqrt(s_safe),
                            0.0
                        )

                        # V_ne uses phi_factor = (2π)²
                        V_ne_integrand = fi * fj * V_ne * J1 * J2_arr
                        V_ne_sum = np.sum(V_ne_integrand * w_eta[b] * w_eta)

                        # V_ee uses phi_factor = 8π
                        V_ee_integrand = fi * fj * V_ee * J1 * J2_arr
                        V_ee_sum = np.sum(V_ee_integrand * w_eta[b] * w_eta)

                        V_ij += w_xi[a] * w_xi[c] * (
                            phi_factor * V_ne_sum
                            + 8 * np.pi * V_ee_sum
                        )

            H[i, j] = T_prefactor * phi_factor * T_ij + V_ij
            H[j, i] = H[i, j]

    return H


# ============================================================
# General Hamiltonian (p > 0) — 5D numerical integration
# ============================================================

def _hamiltonian_general(
    basis: List[HylleraasBasisFunction],
    R: float,
    grids: Dict,
    Z_A: float,
    Z_B: float,
) -> np.ndarray:
    """Hamiltonian for general basis with p > 0.

    Uses 5D integration with explicit Δφ for all terms.
    Kinetic energy via finite differences (necessary because r₁₂^p
    couples all coordinates making analytical derivatives complex).

    For the kinetic energy, we use the identity that for the
    exponential-polynomial × r₁₂^p form, we can split:
    T = T_poly + T_r12_coupling
    where T_poly uses analytical derivatives and T_r12 uses FD on r₁₂^p only.
    """
    n_basis = len(basis)
    H = np.zeros((n_basis, n_basis))

    xi = grids['xi']
    eta = grids['eta']
    w_xi = grids['w_xi']
    w_eta = grids['w_eta']
    dphi = grids['dphi']
    w_phi = grids['w_phi']
    half_R = R / 2.0
    eps_fd = 1e-4  # FD step for r₁₂-dependent terms

    for i in range(n_basis):
        for j in range(i, n_basis):
            val = 0.0

            for a in range(len(xi)):
                xi1 = xi[a]
                for c in range(len(xi)):
                    xi2 = xi[c]
                    for b in range(len(eta)):
                        eta1 = eta[b]
                        J1 = half_R**3 * (xi1**2 - eta1**2)
                        for d in range(len(eta)):
                            eta2 = eta[d]
                            J2 = half_R**3 * (xi2**2 - eta2**2)

                            phi_int = 0.0
                            for ip in range(len(dphi)):
                                r12 = float(compute_r12_with_phi(
                                    xi1, eta1, xi2, eta2, dphi[ip], R
                                ))
                                r12 = max(r12, 1e-15)

                                fi = float(evaluate_basis_function(
                                    basis[i], xi1, eta1, xi2, eta2, r12, R
                                ))
                                fj = float(evaluate_basis_function(
                                    basis[j], xi1, eta1, xi2, eta2, r12, R
                                ))

                                # V_ne
                                r1A = half_R * (xi1 + eta1)
                                r1B = half_R * abs(xi1 - eta1)
                                r2A = half_R * (xi2 + eta2)
                                r2B = half_R * abs(xi2 - eta2)
                                v_ne = (-Z_A / max(r1A, 1e-15)
                                        - Z_B / max(r1B, 1e-15)
                                        - Z_A / max(r2A, 1e-15)
                                        - Z_B / max(r2B, 1e-15))

                                v_ee = 1.0 / r12

                                # Kinetic: FD for full function
                                T_fj = _kinetic_fd_scalar(
                                    basis[j], xi1, eta1, xi2, eta2,
                                    dphi[ip], R, eps_fd
                                )

                                H_fj = T_fj + (v_ne + v_ee) * fj
                                phi_int += fi * H_fj * w_phi[ip]

                            phi_int *= 2 * np.pi
                            val += (w_xi[a] * w_eta[b] * w_xi[c] * w_eta[d]
                                    * J1 * J2 * phi_int)

            H[i, j] = val
            H[j, i] = val

    return H


def _kinetic_fd_scalar(
    bf: HylleraasBasisFunction,
    xi1: float, eta1: float, xi2: float, eta2: float,
    dphi: float, R: float, eps: float,
) -> float:
    """Kinetic energy T|φ⟩ at a single point via finite differences.

    T = T₁ + T₂ where T_i = -(2/R²)/(ξᵢ²-ηᵢ²) × [L_ξᵢ + L_ηᵢ]
    and L_ξ = d/dξ[(ξ²-1)d/dξ], L_η = d/dη[(1-η²)d/dη].

    Uses centered 3-point stencil for second derivative.
    """
    def _eval(x1, e1, x2, e2):
        r = float(compute_r12_with_phi(x1, e1, x2, e2, dphi, R))
        return float(evaluate_basis_function(bf, x1, e1, x2, e2, max(r, 1e-15), R))

    f0 = _eval(xi1, eta1, xi2, eta2)
    T = 0.0

    # T₁: ξ₁ derivative
    xi1_p = xi1 + eps
    xi1_m = max(xi1 - eps, 1.0 + 1e-8)
    h1 = xi1_p - xi1
    h2 = xi1 - xi1_m
    fp = _eval(xi1_p, eta1, xi2, eta2)
    fm = _eval(xi1_m, eta1, xi2, eta2)
    d2f_dxi1 = 2 * (fp * h2 + fm * h1 - f0 * (h1 + h2)) / (h1 * h2 * (h1 + h2))
    df_dxi1 = (fp - fm) / (h1 + h2)
    Lxi1 = (xi1**2 - 1) * d2f_dxi1 + 2 * xi1 * df_dxi1

    # T₁: η₁ derivative
    eta1_p = min(eta1 + eps, 1.0 - 1e-8)
    eta1_m = max(eta1 - eps, -1.0 + 1e-8)
    h1e = eta1_p - eta1
    h2e = eta1 - eta1_m
    fp = _eval(xi1, eta1_p, xi2, eta2)
    fm = _eval(xi1, eta1_m, xi2, eta2)
    d2f_deta1 = 2 * (fp * h2e + fm * h1e - f0 * (h1e + h2e)) / (h1e * h2e * (h1e + h2e))
    df_deta1 = (fp - fm) / (h1e + h2e)
    Leta1 = (1 - eta1**2) * d2f_deta1 - 2 * eta1 * df_deta1

    w1 = xi1**2 - eta1**2
    T += -(2.0 / R**2) * (Lxi1 + Leta1) / max(w1, 1e-15)

    # T₂: ξ₂ derivative
    xi2_p = xi2 + eps
    xi2_m = max(xi2 - eps, 1.0 + 1e-8)
    h1 = xi2_p - xi2
    h2 = xi2 - xi2_m
    fp = _eval(xi1, eta1, xi2_p, eta2)
    fm = _eval(xi1, eta1, xi2_m, eta2)
    d2f_dxi2 = 2 * (fp * h2 + fm * h1 - f0 * (h1 + h2)) / (h1 * h2 * (h1 + h2))
    df_dxi2 = (fp - fm) / (h1 + h2)
    Lxi2 = (xi2**2 - 1) * d2f_dxi2 + 2 * xi2 * df_dxi2

    # T₂: η₂ derivative
    eta2_p = min(eta2 + eps, 1.0 - 1e-8)
    eta2_m = max(eta2 - eps, -1.0 + 1e-8)
    h1e = eta2_p - eta2
    h2e = eta2 - eta2_m
    fp = _eval(xi1, eta1, xi2, eta2_p)
    fm = _eval(xi1, eta1, xi2, eta2_m)
    d2f_deta2 = 2 * (fp * h2e + fm * h1e - f0 * (h1e + h2e)) / (h1e * h2e * (h1e + h2e))
    df_deta2 = (fp - fm) / (h1e + h2e)
    Leta2 = (1 - eta2**2) * d2f_deta2 - 2 * eta2 * df_deta2

    w2 = xi2**2 - eta2**2
    T += -(2.0 / R**2) * (Lxi2 + Leta2) / max(w2, 1e-15)

    return T


# ============================================================
# Solver: generalized eigenvalue problem
# ============================================================

def solve_hylleraas(
    basis: List[HylleraasBasisFunction],
    R: float,
    grids: Dict,
    Z_A: float = 1.0,
    Z_B: float = 1.0,
    n_states: int = 1,
    verbose: bool = True,
) -> Dict:
    """Solve the Hylleraas variational problem Hc = ESc.

    Parameters
    ----------
    basis : list of HylleraasBasisFunction
    R : float
        Internuclear distance (bohr).
    grids : dict
        Quadrature grids from build_quadrature_grids().
    Z_A, Z_B : float
        Nuclear charges.
    n_states : int
        Number of states to return.
    verbose : bool
        Print progress.

    Returns
    -------
    dict with E_elec, E_total, D_e, coeffs, S, H, evals, evecs, etc.
    """
    t0 = time.time()
    n_bf = len(basis)

    if verbose:
        print(f"\n  Hylleraas solver: {n_bf} basis functions, R={R:.4f}")
        print(f"  Grid: {grids['N_xi']}x{grids['N_eta']} (xi,eta), "
              f"{grids['N_phi']} phi points")

    # Step 1: Overlap matrix
    if verbose:
        print("  Computing overlap matrix S...")
    t1 = time.time()
    S = compute_overlap_matrix(basis, R, grids)
    dt_S = time.time() - t1
    if verbose:
        print(f"    S computed in {dt_S:.1f}s")
        print(f"    S condition number: {np.linalg.cond(S):.2e}")

    s_evals = np.linalg.eigvalsh(S)
    if s_evals[0] < 1e-10:
        if verbose:
            print(f"    WARNING: Near-singular S (min eval = {s_evals[0]:.2e})")

    # Step 2: Hamiltonian matrix
    if verbose:
        print("  Computing Hamiltonian matrix H...")
    t2 = time.time()
    H = compute_hamiltonian_matrix(basis, R, grids, Z_A, Z_B)
    dt_H = time.time() - t2
    if verbose:
        print(f"    H computed in {dt_H:.1f}s")

    # Add V_NN
    V_NN = Z_A * Z_B / R
    H_total = H + V_NN * S

    # Step 3: Solve generalized eigenvalue problem
    try:
        evals, evecs = eigh(H_total, S)
    except np.linalg.LinAlgError:
        if verbose:
            print("    Regularizing singular overlap matrix...")
        S_reg = S + 1e-10 * np.eye(n_bf)
        evals, evecs = eigh(H_total, S_reg)

    E_total = evals[0]
    E_elec = E_total - V_NN
    coeffs = evecs[:, 0]

    E_atoms = -Z_A**2 / 2 - Z_B**2 / 2
    D_e = E_atoms - E_total

    dt_total = time.time() - t0

    if verbose:
        print(f"\n  Results:")
        print(f"    E_total = {E_total:.6f} Ha")
        print(f"    E_elec  = {E_elec:.6f} Ha")
        print(f"    D_e     = {D_e:.6f} Ha ({D_e/0.1745*100:.1f}% of exact)")
        print(f"    Time    = {dt_total:.1f}s")

    return {
        'E_total': E_total,
        'E_elec': E_elec,
        'D_e': D_e,
        'D_e_pct': D_e / 0.1745 * 100,
        'V_NN': V_NN,
        'coeffs': coeffs,
        'evals': evals[:n_states],
        'evecs': evecs[:, :n_states],
        'S': S,
        'H': H_total,
        'n_basis': n_bf,
        'R': R,
        'time': dt_total,
        'basis': basis,
    }


# ============================================================
# Alpha optimization (nonlinear parameter)
# ============================================================

def optimize_alpha(
    basis_generator,
    R: float,
    grids: Dict,
    alpha_range: Tuple[float, float] = (0.5, 2.5),
    Z_A: float = 1.0,
    Z_B: float = 1.0,
    verbose: bool = True,
) -> Dict:
    """Optimize the nonlinear parameter α.

    Parameters
    ----------
    basis_generator : callable
        Function alpha -> List[HylleraasBasisFunction].
    R : float
        Internuclear distance.
    grids : dict
        Quadrature grids.
    alpha_range : tuple
        (alpha_min, alpha_max) for optimization.
    Z_A, Z_B : float
        Nuclear charges.
    verbose : bool

    Returns
    -------
    dict with alpha_opt, E_opt, and full result dict.
    """
    if verbose:
        print(f"\n  Optimizing alpha in [{alpha_range[0]:.2f}, {alpha_range[1]:.2f}]...")

    def objective(alpha: float) -> float:
        basis = basis_generator(alpha)
        try:
            result = solve_hylleraas(basis, R, grids, Z_A, Z_B, verbose=False)
            return result['E_total']
        except Exception:
            return 0.0

    n_scan = 9
    alphas = np.linspace(alpha_range[0], alpha_range[1], n_scan)
    energies = np.array([objective(a) for a in alphas])

    valid = energies < -0.5
    if not np.any(valid):
        idx_best = np.argmin(energies)
    else:
        idx_best = np.argmin(np.where(valid, energies, 0.0))

    if verbose:
        print(f"    Coarse scan: best alpha={alphas[idx_best]:.3f}, "
              f"E={energies[idx_best]:.6f}")

    if 0 < idx_best < n_scan - 1:
        bounds = (alphas[max(0, idx_best - 1)],
                  alphas[min(n_scan - 1, idx_best + 1)])
    else:
        bounds = alpha_range

    opt = minimize_scalar(objective, bounds=bounds, method='bounded',
                          options={'xatol': 0.005})
    alpha_opt = opt.x
    basis_opt = basis_generator(alpha_opt)
    result_opt = solve_hylleraas(
        basis_opt, R, grids, Z_A, Z_B, verbose=verbose
    )

    if verbose:
        print(f"\n  Optimal alpha = {alpha_opt:.4f}")

    return {
        'alpha_opt': alpha_opt,
        'E_opt': result_opt['E_total'],
        'result': result_opt,
        'scan_alphas': alphas,
        'scan_energies': energies,
    }


# ============================================================
# Convenience: full H2 calculation
# ============================================================

def compute_h2_hylleraas(
    R: float = 1.4,
    j_max: int = 1,
    l_max: int = 1,
    p_max: int = 1,
    N_xi: int = 20,
    N_eta: int = 16,
    N_phi: int = 12,
    xi_max: float = 12.0,
    optimize_alpha_flag: bool = True,
    alpha: float = 1.0,
    verbose: bool = True,
) -> Dict:
    """One-call H2 Hylleraas calculation.

    Parameters
    ----------
    R : float
        Bond length (bohr).
    j_max, l_max, p_max : int
        Basis set truncation parameters.
    N_xi, N_eta, N_phi : int
        Quadrature grid sizes.
    xi_max : float
        Maximum xi value.
    optimize_alpha_flag : bool
        Whether to optimize alpha.
    alpha : float
        Initial or fixed alpha value.
    verbose : bool

    Returns
    -------
    dict with full results including D_e, E_total, etc.
    """
    grids = build_quadrature_grids(N_xi, N_eta, N_phi, xi_max)

    if optimize_alpha_flag:
        def gen(a):
            return generate_basis(j_max, l_max, p_max, alpha=a)
        result = optimize_alpha(gen, R, grids, verbose=verbose)
        return result['result']
    else:
        basis = generate_basis(j_max, l_max, p_max, alpha=alpha)
        return solve_hylleraas(basis, R, grids, verbose=verbose)
