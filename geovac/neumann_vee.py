"""
Neumann expansion for V_ee (1/r₁₂) in prolate spheroidal coordinates.

The Coulomb kernel expands as (for σ states, m=0 only):

    <1/|r₁-r₂|>_φ = (2/R) Σ_l (2l+1) P_l(ξ_<) Q_l(ξ_>) P_l(η₁) P_l(η₂)

where ξ_< = min(ξ₁,ξ₂), ξ_> = max(ξ₁,ξ₂).

This converts V_ee matrix elements from 5D numerical quadrature to
sums of 1D algebraic integrals — exact within the Neumann truncation.

References:
    Roothaan (1951) J. Chem. Phys. 19, 1445
    Shavitt (1963) Methods in Computational Physics, Vol. 2
    Harris & Michels (1966) Adv. Chem. Phys. 13, 205
"""

import numpy as np
from scipy.special import exp1  # E₁(x) = exponential integral
from typing import Tuple, List, Dict, Optional
from functools import lru_cache

# Module-level cache for B_l tables (keyed on l_max, p_max, alpha)
_Bl_cache: Dict[Tuple[int, int, float], np.ndarray] = {}


# ============================================================
# Auxiliary integral C_l(p): η moments
# ============================================================

def compute_Cl_table(p_max: int, l_max: int) -> np.ndarray:
    """Compute table of C_l(p) = ∫₋₁¹ η^p P_l(η) dη.

    Uses the Legendre recurrence:
        (l+1) C_{l+1}(p) = (2l+1) C_l(p+1) - l C_{l-1}(p)

    Selection rule: C_l(p) = 0 if p < l or (p-l) is odd.

    Parameters
    ----------
    p_max : int
        Maximum η power needed.
    l_max : int
        Maximum Legendre order.

    Returns
    -------
    C : ndarray, shape (l_max+1, p_max+3)
        C[l, p] = C_l(p). Extra columns for recurrence workspace.
    """
    # Need p up to p_max + l_max for the recurrence chain:
    # C_{l+1}(p) uses C_l(p+1), so C_{l_max}(p) needs workspace to p + l_max
    p_ext = p_max + l_max
    C = np.zeros((l_max + 1, p_ext + 1))

    # Base cases: l = 0
    for p in range(p_ext + 1):
        if p % 2 == 0:
            C[0, p] = 2.0 / (p + 1)
        # else: C[0, p] = 0 (odd powers integrate to zero with P_0=1)

    # Base cases: l = 1 (P_1(η) = η)
    if l_max >= 1:
        for p in range(p_ext + 1):
            if p % 2 == 1:  # p odd → p+1 even → ∫η^{p+1} dη = 2/(p+2)
                C[1, p] = 2.0 / (p + 2)

    # Recurrence upward in l
    for l in range(1, l_max):
        for p in range(p_ext + 1):
            # (l+1) C_{l+1}(p) = (2l+1) C_l(p+1) - l C_{l-1}(p)
            if p + 1 <= p_ext:
                C[l + 1, p] = ((2 * l + 1) * C[l, p + 1] - l * C[l - 1, p]) / (l + 1)

    return C[:, :p_max + 1]


# ============================================================
# Auxiliary integral A_n(α): incomplete gamma base functions
# ============================================================

def compute_An_table(n_max: int, alpha: float) -> np.ndarray:
    """Compute A_n(α) = ∫₁^∞ ξ^n e^{-αξ} dξ for n = 0, 1, ..., n_max.

    Uses upward recurrence: A_n(α) = [n·A_{n-1}(α) + e^{-α}] / α
    Starting from A_0(α) = e^{-α} / α.

    Parameters
    ----------
    n_max : int
        Maximum power.
    alpha : float
        Exponential parameter (must be > 0).

    Returns
    -------
    A : ndarray, shape (n_max+1,)
        A[n] = A_n(α).
    """
    A = np.zeros(n_max + 1)
    e_neg_a = np.exp(-alpha)

    A[0] = e_neg_a / alpha
    for n in range(1, n_max + 1):
        A[n] = (n * A[n - 1] + e_neg_a) / alpha

    return A


def compute_An_negative(n_min: int, alpha: float) -> Dict[int, float]:
    """Compute A_n(α) for negative n using downward recurrence.

    A_{-1}(α) = E₁(α) (exponential integral)
    A_{-n-1}(α) = [e^{-α} - α·A_{-n}(α)] / n   for n ≥ 1

    Parameters
    ----------
    n_min : int
        Most negative index needed (e.g., -5).
    alpha : float
        Exponential parameter.

    Returns
    -------
    A_neg : dict mapping int -> float
        A_neg[n] = A_n(α) for n = n_min, ..., -1.
    """
    A_neg = {}
    e_neg_a = np.exp(-alpha)

    # A_{-1} = E₁(α)
    A_neg[-1] = float(exp1(alpha))

    # Downward: A_{-n-1} = [e^{-α} - α·A_{-n}] / n
    for n in range(1, abs(n_min)):
        idx = -(n + 1)
        A_neg[idx] = (e_neg_a - alpha * A_neg[-n]) / n

    return A_neg


# ============================================================
# Auxiliary integral A_l(p, α): first-kind Legendre moments
# ============================================================

def compute_Al_table(
    l_max: int, p_max: int, alpha: float,
) -> np.ndarray:
    """Compute A_l(p, α) = ∫₁^∞ ξ^p e^{-αξ} P_l(ξ) dξ.

    Uses the Legendre recurrence in l:
        (l+1) A_{l+1}(p, α) = (2l+1) A_l(p+1, α) - l A_{l-1}(p, α)

    Parameters
    ----------
    l_max : int
        Maximum Legendre order.
    p_max : int
        Maximum ξ power.
    alpha : float
        Exponential parameter.

    Returns
    -------
    Al : ndarray, shape (l_max+1, p_max+1)
        Al[l, p] = A_l(p, α).
    """
    # Need A_n up to p_max + l_max + 1 for the recurrence
    An = compute_An_table(p_max + l_max + 1, alpha)

    # Workspace: need p up to p_max + l_max for recurrence
    p_ext = p_max + l_max + 1
    Al = np.zeros((l_max + 1, p_ext + 1))

    # l = 0: A_0(p, α) = A_p(α)
    for p in range(min(p_ext + 1, len(An))):
        Al[0, p] = An[p]

    # l = 1: A_1(p, α) = A_{p+1}(α)  [since P_1(ξ) = ξ]
    if l_max >= 1:
        for p in range(min(p_ext, len(An) - 1)):
            Al[1, p] = An[p + 1]

    # Recurrence upward in l
    for l in range(1, l_max):
        for p in range(p_ext):
            Al[l + 1, p] = ((2 * l + 1) * Al[l, p + 1] - l * Al[l - 1, p]) / (l + 1)

    return Al[:, :p_max + 1]


# ============================================================
# Auxiliary integral B_l(p, α): second-kind Legendre moments
# ============================================================

def compute_B0_table(p_max: int, alpha: float) -> np.ndarray:
    """Compute B_0(p, α) = ∫₁^∞ ξ^p e^{-αξ} Q_0(ξ) dξ.

    Q_0(ξ) = ½ ln[(ξ+1)/(ξ-1)], which has a logarithmic singularity
    at ξ=1. We use the identity:

        B_0(p, α) = ½[L⁺_p(α) - L⁻_p(α)]

    where L⁺_p = ∫₁^∞ ξ^p e^{-αξ} ln(ξ+1) dξ (smooth integrand)
          L⁻_p = ∫₁^∞ ξ^p e^{-αξ} ln(ξ-1) dξ (log singularity at ξ=1)

    Both L⁺ and L⁻ satisfy recurrences in p via integration by parts:
        L±_p(α) = [p·L±_{p-1}(α) + A_p(α)] / α   (wrong sign — use direct)

    We compute B_0(0, α) and B_0(1, α) directly using the formula:
        B_0(0, α) = ½ ∫₁^∞ e^{-αξ} ln((ξ+1)/(ξ-1)) dξ

    Then use the recurrence in p:
        B_0(p, α) = [p·B_0(p-1, α) + A_p(α)·½·ln(2)] / α  -- NO, IBP gives:

    Actually, use direct upward recurrence from IBP on the full integrand:
        ∫₁^∞ ξ^p e^{-αξ} Q_0(ξ) dξ
    Let u = ξ^p Q_0(ξ), dv = e^{-αξ} dξ, then:
        IBP gives: B_0(p) = (1/α)[e^{-α}·Q_0(1)·1 + p·B_0(p-1) + G_0(p)]
    where G_0(p) = ∫₁^∞ ξ^p e^{-αξ} (-1/(ξ²-1)) dξ  [from Q_0'(ξ) = -1/(ξ²-1)]

    This is getting complex. Use direct numerical integration (scipy.integrate.quad
    handles the logarithmic singularity efficiently).

    Parameters
    ----------
    p_max : int
        Maximum ξ power.
    alpha : float
        Exponential parameter.

    Returns
    -------
    B0 : ndarray, shape (p_max+1,)
        B0[p] = B_0(p, α).
    """
    from scipy.integrate import quad

    B0 = np.zeros(p_max + 1)

    def Q0(xi: float) -> float:
        if xi <= 1.0:
            return 0.0
        return 0.5 * np.log((xi + 1.0) / (xi - 1.0))

    for p in range(p_max + 1):
        def integrand(xi: float, p_: int = p) -> float:
            return xi**p_ * np.exp(-alpha * xi) * Q0(xi)

        # Split integral: [1, 1+eps] has log singularity, [1+eps, ∞) smooth
        # quad with weight='alg-loga' handles algebraic-log singularities
        # Or just use points= to help quad find the singularity
        val, err = quad(integrand, 1.0, np.inf, limit=200, epsabs=1e-14)
        B0[p] = val

    return B0


def compute_Bl_table(
    l_max: int, p_max: int, alpha: float,
) -> np.ndarray:
    """Compute B_l(p, alpha) = integral_1^inf xi^p e^{-alpha*xi} Q_l(xi) dxi.

    Uses scipy.integrate.quad with vectorized Q_l evaluation via lqmn.
    All l values are computed simultaneously for each quadrature point
    by using a shared integration grid.

    Parameters
    ----------
    l_max : int
        Maximum Legendre order.
    p_max : int
        Maximum xi power.
    alpha : float
        Exponential parameter.

    Returns
    -------
    Bl : ndarray, shape (l_max+1, p_max+1)
        Bl[l, p] = B_l(p, alpha).
    """
    from scipy.integrate import quad
    from scipy.special import lqmn

    # Check cache
    cache_key = (l_max, p_max, round(alpha, 10))
    if cache_key in _Bl_cache:
        return _Bl_cache[cache_key].copy()

    # Check if a larger cached table exists
    for (cl, cp, ca), cached in _Bl_cache.items():
        if ca == round(alpha, 10) and cl >= l_max and cp >= p_max:
            return cached[:l_max + 1, :p_max + 1].copy()

    Bl = np.zeros((l_max + 1, p_max + 1))

    for l in range(l_max + 1):
        for p in range(p_max + 1):
            def integrand(xi: float, l_: int = l, p_: int = p) -> float:
                Qvals, _ = lqmn(0, l_, xi)
                Ql = Qvals[0, l_]
                return xi**p_ * np.exp(-alpha * xi) * Ql

            val, err = quad(integrand, 1.0, np.inf, limit=200, epsabs=1e-14)
            Bl[l, p] = val

    _Bl_cache[cache_key] = Bl.copy()
    return Bl


# ============================================================
# Legendre polynomial coefficients
# ============================================================

def legendre_poly_coeffs(l: int) -> np.ndarray:
    """Return coefficients of P_l(x) as a polynomial.

    P_l(x) = Σ_k c_k x^k, returns array c of length l+1.

    Uses the recurrence: (l+1)P_{l+1} = (2l+1)x P_l - l P_{l-1}
    """
    if l == 0:
        return np.array([1.0])
    if l == 1:
        return np.array([0.0, 1.0])

    # Build up using recurrence
    P_prev = np.zeros(l + 1)
    P_prev[0] = 1.0  # P_0
    P_curr = np.zeros(l + 1)
    P_curr[1] = 1.0  # P_1

    for n in range(1, l):
        P_next = np.zeros(l + 1)
        # (n+1) P_{n+1}(x) = (2n+1) x P_n(x) - n P_{n-1}(x)
        # x P_n(x) shifts coefficients up by 1
        for k in range(n + 2):
            if k > 0:
                P_next[k] += (2 * n + 1) * P_curr[k - 1]
            P_next[k] -= n * P_prev[k]
        P_next /= (n + 1)
        P_prev = P_curr.copy()
        P_curr = P_next.copy()

    return P_curr


def poly_product_coeffs(p: int, l: int) -> np.ndarray:
    """Coefficients of ξ^p · P_l(ξ) as a polynomial.

    Returns array c of length p+l+1 where ξ^p P_l(ξ) = Σ_k c[k] ξ^k.
    """
    Pl = legendre_poly_coeffs(l)
    # Multiply by ξ^p: shift all coefficients up by p
    result = np.zeros(p + l + 1)
    for k in range(l + 1):
        result[p + k] = Pl[k]
    return result


# ============================================================
# 2D ξ-integral X_l(p₁, p₂, α) via integration by parts
# ============================================================

def compute_Xl(
    l_max: int,
    p1_max: int,
    p2_max: int,
    alpha: float,
) -> np.ndarray:
    """Compute X_l(p₁, p₂, α) = ∫∫ ξ₁^{p₁} ξ₂^{p₂} e^{-α(ξ₁+ξ₂)} P_l(ξ_<) Q_l(ξ_>) dξ₁dξ₂.

    Split into two ordered regions (ξ₁<ξ₂ and ξ₁>ξ₂):

    X_l = ∫₁^∞ ξ₂^{p₂} e^{-αξ₂} Q_l(ξ₂) S_l(p₁,α,ξ₂) dξ₂ + (1↔2)

    where S_l(p,α,ξ₀) = ∫₁^{ξ₀} ξ^p e^{-αξ} P_l(ξ) dξ is computed by
    integration by parts, yielding:

    S_l(p,α,ξ₀) = A_l(p,α) - e^{-αξ₀} Σ_j a_j Σ_{k=0}^j [j!/(j-k)!/α^{k+1}] ξ₀^{j-k}
                   + e^{-α} Σ_j a_j Σ_{k=0}^j [j!/(j-k)!/α^{k+1}]

    where {a_j} are coefficients of ξ^p P_l(ξ).

    Substituting into the outer integral gives correction terms involving
    B_l(·, 2α) (the factor 2α from combining e^{-αξ₂} and e^{-αξ₂}).

    Parameters
    ----------
    l_max : int
        Maximum Neumann order.
    p1_max, p2_max : int
        Maximum ξ powers.
    alpha : float
        Exponential parameter.

    Returns
    -------
    Xl : ndarray, shape (l_max+1, p1_max+1, p2_max+1)
        Xl[l, p1, p2] = X_l(p₁, p₂, α).
    """
    pmax = max(p1_max, p2_max)

    # Precompute tables at both α and 2α
    # Need large enough p range for the correction terms
    p_need = pmax + l_max + 2
    Al_a = compute_Al_table(l_max, p_need, alpha)
    Bl_a = compute_Bl_table(l_max, p_need, alpha)
    Bl_2a = compute_Bl_table(l_max, p_need + pmax + l_max, 2 * alpha)

    # Precompute polynomial coefficients for ξ^p P_l(ξ)
    # and the IBP correction coefficients
    e_neg_a = np.exp(-alpha)

    Xl = np.zeros((l_max + 1, p1_max + 1, p2_max + 1))

    for l in range(l_max + 1):
        Pl_coeffs = legendre_poly_coeffs(l)

        for p1 in range(p1_max + 1):
            for p2 in range(p2_max + 1):
                # Region I: ξ₁ < ξ₂
                # = A_l(p1,α) · B_l(p2,α) - correction_I
                I1 = Al_a[l, p1] * Bl_a[l, p2]

                # Correction from IBP of S_l(p1, α, ξ₂):
                # For each monomial ξ^j in ξ^{p1} P_l(ξ):
                #   ∫₁^{ξ₂} ξ^j e^{-αξ} dξ = A_j(α) - e^{-αξ₂} Σ_{k=0}^j j!/((j-k)! α^{k+1}) ξ₂^{j-k}
                # The correction involves the e^{-αξ₂} terms, which when
                # multiplied by the outer e^{-αξ₂} give e^{-2αξ₂} terms.
                corr1 = 0.0
                poly1 = poly_product_coeffs(p1, l)
                for j in range(len(poly1)):
                    if abs(poly1[j]) < 1e-30:
                        continue
                    for k in range(j + 1):
                        # coefficient: poly1[j] * j! / ((j-k)! * α^{k+1})
                        coeff = poly1[j]
                        for f in range(j - k + 1, j + 1):  # j!/(j-k)!
                            coeff *= f
                        coeff /= alpha ** (k + 1)

                        # The outer integral becomes:
                        # ∫₁^∞ ξ₂^{p2+j-k} e^{-2αξ₂} Q_l(ξ₂) dξ₂ = B_l(p2+j-k, 2α)
                        p2_shifted = p2 + j - k
                        if p2_shifted >= 0 and p2_shifted < Bl_2a.shape[1]:
                            corr1 += coeff * Bl_2a[l, p2_shifted]

                I1 -= corr1

                # Region II: ξ₂ < ξ₁ (swap p1 ↔ p2)
                I2 = Al_a[l, p2] * Bl_a[l, p1]

                corr2 = 0.0
                poly2 = poly_product_coeffs(p2, l)
                for j in range(len(poly2)):
                    if abs(poly2[j]) < 1e-30:
                        continue
                    for k in range(j + 1):
                        coeff = poly2[j]
                        for f in range(j - k + 1, j + 1):
                            coeff *= f
                        coeff /= alpha ** (k + 1)

                        p1_shifted = p1 + j - k
                        if p1_shifted >= 0 and p1_shifted < Bl_2a.shape[1]:
                            corr2 += coeff * Bl_2a[l, p1_shifted]

                I2 -= corr2

                Xl[l, p1, p2] = I1 + I2

    return Xl


# ============================================================
# V_ee matrix element via Neumann expansion
# ============================================================

def compute_vee_matrix_neumann(
    basis: list,
    R: float,
    l_max: int = 20,
    verbose: bool = False,
) -> np.ndarray:
    """Compute V_ee matrix using Neumann expansion of 1/r₁₂.

    For σ states (m=0), the Neumann expansion gives:
        1/r₁₂ = (4/R) Σ_l P_l(ξ_<) Q_l(ξ_>) P_l(η₁) P_l(η₂)

    The V_ee matrix element ⟨φ_i|1/r₁₂|φ_j⟩ factorizes (after handling
    the Jacobian) into products of η-integrals C_l and ξ-integrals X_l.

    Parameters
    ----------
    basis : list of HylleraasBasisFunction
        Basis functions (p=0 only for now).
    R : float
        Internuclear distance (bohr).
    l_max : int
        Truncation of Neumann series.
    verbose : bool
        Print diagnostic info.

    Returns
    -------
    Vee : ndarray, shape (n_basis, n_basis)
        V_ee matrix in the basis.
    """
    n_bf = len(basis)

    # Determine maximum quantum numbers
    j_max = max(bf.j for bf in basis)
    k_max = max(bf.k for bf in basis)
    l_max_basis = max(bf.l for bf in basis)
    m_max = max(bf.m for bf in basis)

    # Combined powers from bra × ket × Jacobian:
    # ξ₁ power: j_i + j_j + 2 (from Jacobian ξ₁²), max = 2*j_max + 2
    # η₁ power: l_i + l_j + 2 (from Jacobian η₁²), max = 2*l_max_basis + 2
    p_xi_max = 2 * max(j_max, k_max) + 2
    p_eta_max = 2 * max(l_max_basis, m_max) + 2

    # Selection rule: C_l(q) = 0 if l > q or (q-l) odd
    # So effective l_max is min(l_max, p_eta_max)
    l_eff = min(l_max, p_eta_max)

    if verbose:
        print(f"  Neumann V_ee: {n_bf} basis functions, l_max={l_eff}")
        print(f"  xi powers up to {p_xi_max}, eta powers up to {p_eta_max}")

    # Get all unique alpha values
    alphas = sorted(set(bf.alpha for bf in basis))
    if len(alphas) > 1:
        raise NotImplementedError("Multiple alpha values not yet supported")
    alpha = alphas[0]
    alpha_tot = 2 * alpha  # bra + ket exponentials

    # Precompute tables
    Cl = compute_Cl_table(p_eta_max, l_eff)

    # Precompute A_l and B_l tables ONCE for all elements
    p_need = p_xi_max + l_eff + 2
    if verbose:
        import time as _time
        _t0 = _time.time()
    Al_a = compute_Al_table(l_eff, p_need, alpha_tot)
    Bl_a = compute_Bl_table(l_eff, p_need, alpha_tot)
    Bl_2a = compute_Bl_table(l_eff, p_need + p_xi_max + l_eff + 1, 2 * alpha_tot)
    if verbose:
        print(f"  Precomputed auxiliary tables in {_time.time()-_t0:.1f}s")

    # Precompute X_l table for all needed (l, p1, p2) combinations
    Xl_cache = _precompute_Xl_table(
        l_eff, p_xi_max, alpha_tot, Al_a, Bl_a, Bl_2a
    )

    # Prefactor: (2/R) × (R/2)⁶ × (2π)² = (2/R)(R⁶/64)(4π²) = π²R⁵/8
    # The (2l+1) factor is applied inside the Neumann sum loop.
    prefactor = np.pi**2 * R**5 / 8.0

    # For each pair (i,j), compute V_ee
    Vee = np.zeros((n_bf, n_bf))

    for i in range(n_bf):
        for j in range(i, n_bf):
            val = _compute_vee_element_cached(
                basis[i], basis[j], l_eff, Cl, Xl_cache, prefactor
            )
            Vee[i, j] = val
            Vee[j, i] = val

    return Vee


def _precompute_Xl_table(
    l_max: int, p_max: int, alpha: float,
    Al_a: np.ndarray, Bl_a: np.ndarray, Bl_2a: np.ndarray,
) -> np.ndarray:
    """Precompute X_l(p₁, p₂, α) for all needed combinations.

    Uses precomputed A_l and B_l tables for efficiency.

    Returns
    -------
    Xl : ndarray, shape (l_max+1, p_max+1, p_max+1)
    """
    Xl = np.zeros((l_max + 1, p_max + 1, p_max + 1))

    for l in range(l_max + 1):
        for p1 in range(p_max + 1):
            poly1 = poly_product_coeffs(p1, l)
            for p2 in range(p_max + 1):
                # Region I: ξ₁ < ξ₂
                I1 = Al_a[l, p1] * Bl_a[l, p2]

                corr1 = 0.0
                for j in range(len(poly1)):
                    if abs(poly1[j]) < 1e-30:
                        continue
                    for k in range(j + 1):
                        coeff = poly1[j]
                        for f in range(j - k + 1, j + 1):
                            coeff *= f
                        coeff /= alpha ** (k + 1)
                        p2s = p2 + j - k
                        if 0 <= p2s < Bl_2a.shape[1]:
                            corr1 += coeff * Bl_2a[l, p2s]
                I1 -= corr1

                # Region II: ξ₂ < ξ₁
                I2 = Al_a[l, p2] * Bl_a[l, p1]

                poly2 = poly_product_coeffs(p2, l)
                corr2 = 0.0
                for j in range(len(poly2)):
                    if abs(poly2[j]) < 1e-30:
                        continue
                    for k in range(j + 1):
                        coeff = poly2[j]
                        for f in range(j - k + 1, j + 1):
                            coeff *= f
                        coeff /= alpha ** (k + 1)
                        p1s = p1 + j - k
                        if 0 <= p1s < Bl_2a.shape[1]:
                            corr2 += coeff * Bl_2a[l, p1s]
                I2 -= corr2

                Xl[l, p1, p2] = I1 + I2

    return Xl


def _compute_vee_element(
    bf_i, bf_j, alpha: float, R: float,
    l_max: int, Cl: np.ndarray, prefactor: float,
) -> float:
    """Compute a single V_ee matrix element ⟨φ_i|1/r₁₂|φ_j⟩.

    Each symmetrized basis function φ = [g(1,2) + g(2,1)] where
    g(1,2) = e^{-α(ξ₁+ξ₂)} ξ₁^j ξ₂^k η₁^l η₂^m.

    The product φ_i φ_j has up to 4 unsymmetrized terms.
    """
    alpha_tot = 2 * alpha  # bra + ket exponentials

    # Enumerate unsymmetrized term pairs from bra and ket
    terms_i = _get_unsym_terms(bf_i)
    terms_j = _get_unsym_terms(bf_j)

    result = 0.0
    for (ji, ki, li, mi) in terms_i:
        for (jj, kj, lj, mj) in terms_j:
            # Combined powers from bra × ket
            p1 = ji + jj   # ξ₁ power
            p2 = ki + kj   # ξ₂ power
            q1 = li + lj   # η₁ power
            q2 = mi + mj   # η₂ power

            result += _vee_unsym(p1, p2, q1, q2, alpha_tot, R, l_max, Cl, prefactor)

    return result


def _compute_vee_element_cached(
    bf_i, bf_j,
    l_max: int, Cl: np.ndarray, Xl_cache: np.ndarray, prefactor: float,
) -> float:
    """Compute V_ee element using precomputed X_l table."""
    terms_i = _get_unsym_terms(bf_i)
    terms_j = _get_unsym_terms(bf_j)

    result = 0.0
    for (ji, ki, li, mi) in terms_i:
        for (jj, kj, lj, mj) in terms_j:
            p1 = ji + jj
            p2 = ki + kj
            q1 = li + lj
            q2 = mi + mj
            result += _vee_unsym_cached(
                p1, p2, q1, q2, l_max, Cl, Xl_cache, prefactor
            )

    return result


def _get_unsym_terms(bf) -> list:
    """Get unsymmetrized (j, k, l, m) terms from a symmetrized basis function.

    Returns list of (j, k, l, m) tuples. For self-exchange (j=k, l=m),
    returns one term with factor 2 (handled by doubling).
    """
    if bf.is_self_exchange:
        # φ = 2 · g(1,2), but the factor 2 is in evaluate_basis_function
        # For V_ee, ⟨2g_i|V|2g_j⟩ involves 4 terms, all same → 4× one term
        # But we handle the symmetrization factor elsewhere
        return [(bf.j, bf.k, bf.l, bf.m), (bf.k, bf.j, bf.m, bf.l)]
    else:
        return [(bf.j, bf.k, bf.l, bf.m), (bf.k, bf.j, bf.m, bf.l)]


def _vee_unsym(
    p1: int, p2: int, q1: int, q2: int,
    alpha: float, R: float, l_max: int,
    Cl: np.ndarray, prefactor: float,
) -> float:
    """V_ee for a single unsymmetrized term with given combined powers.

    The Jacobian (ξ₁²-η₁²)(ξ₂²-η₂²) expands into 4 sub-terms:
    (+) ξ₁²ξ₂²:  P₁=p1+2, P₂=p2+2, Q₁=q1,   Q₂=q2
    (-) ξ₁²η₂²:  P₁=p1+2, P₂=p2,   Q₁=q1,   Q₂=q2+2
    (-) η₁²ξ₂²:  P₁=p1,   P₂=p2+2, Q₁=q1+2, Q₂=q2
    (+) η₁²η₂²:  P₁=p1,   P₂=p2,   Q₁=q1+2, Q₂=q2+2
    """
    val = 0.0

    jacobian_terms = [
        (+1, p1 + 2, p2 + 2, q1, q2),
        (-1, p1 + 2, p2, q1, q2 + 2),
        (-1, p1, p2 + 2, q1 + 2, q2),
        (+1, p1, p2, q1 + 2, q2 + 2),
    ]

    for sign, P1, P2, Q1, Q2 in jacobian_terms:
        # Compute X_l and C_l for each l, sum Neumann series
        neumann_sum = 0.0
        for l in range(l_max + 1):
            # Selection rule: C_l(Q) = 0 if Q < l or (Q-l) odd
            if Q1 < l or Q2 < l:
                continue
            if (Q1 - l) % 2 != 0 or (Q2 - l) % 2 != 0:
                continue

            c1 = Cl[l, Q1] if Q1 < Cl.shape[1] else 0.0
            c2 = Cl[l, Q2] if Q2 < Cl.shape[1] else 0.0
            if abs(c1) < 1e-30 or abs(c2) < 1e-30:
                continue

            # Compute X_l(P1, P2, α) on the fly
            xl = _compute_Xl_single(l, P1, P2, alpha)
            neumann_sum += (2 * l + 1) * xl * c1 * c2

        val += sign * neumann_sum

    return prefactor * val


def _vee_unsym_cached(
    p1: int, p2: int, q1: int, q2: int,
    l_max: int, Cl: np.ndarray, Xl_cache: np.ndarray, prefactor: float,
) -> float:
    """V_ee for unsymmetrized term using precomputed X_l table."""
    val = 0.0

    jacobian_terms = [
        (+1, p1 + 2, p2 + 2, q1, q2),
        (-1, p1 + 2, p2, q1, q2 + 2),
        (-1, p1, p2 + 2, q1 + 2, q2),
        (+1, p1, p2, q1 + 2, q2 + 2),
    ]

    for sign, P1, P2, Q1, Q2 in jacobian_terms:
        neumann_sum = 0.0
        for l in range(l_max + 1):
            if Q1 < l or Q2 < l:
                continue
            if (Q1 - l) % 2 != 0 or (Q2 - l) % 2 != 0:
                continue

            c1 = Cl[l, Q1] if Q1 < Cl.shape[1] else 0.0
            c2 = Cl[l, Q2] if Q2 < Cl.shape[1] else 0.0
            if abs(c1) < 1e-30 or abs(c2) < 1e-30:
                continue

            xl = Xl_cache[l, P1, P2]
            neumann_sum += (2 * l + 1) * xl * c1 * c2

        val += sign * neumann_sum

    return prefactor * val


def _compute_Xl_single(l: int, p1: int, p2: int, alpha: float) -> float:
    """Compute a single X_l(p₁, p₂, α) value.

    This is less efficient than batch computation but simpler for debugging.
    """
    # Need A_l and B_l at both α and 2α
    p_need = max(p1, p2) + l + 2
    Al_a = compute_Al_table(l, p_need, alpha)
    Bl_a = compute_Bl_table(l, p_need, alpha)
    Bl_2a = compute_Bl_table(l, p_need + max(p1, p2) + l + 1, 2 * alpha)

    # Region I: ξ₁ < ξ₂
    I1 = Al_a[l, p1] * Bl_a[l, p2]

    # Correction from IBP
    poly1 = poly_product_coeffs(p1, l)
    corr1 = 0.0
    for j in range(len(poly1)):
        if abs(poly1[j]) < 1e-30:
            continue
        for k in range(j + 1):
            coeff = poly1[j]
            for f in range(j - k + 1, j + 1):
                coeff *= f
            coeff /= alpha ** (k + 1)
            p2s = p2 + j - k
            if 0 <= p2s < Bl_2a.shape[1]:
                corr1 += coeff * Bl_2a[l, p2s]
    I1 -= corr1

    # Region II: ξ₂ < ξ₁ (swap roles)
    I2 = Al_a[l, p2] * Bl_a[l, p1]

    poly2 = poly_product_coeffs(p2, l)
    corr2 = 0.0
    for j in range(len(poly2)):
        if abs(poly2[j]) < 1e-30:
            continue
        for k in range(j + 1):
            coeff = poly2[j]
            for f in range(j - k + 1, j + 1):
                coeff *= f
            coeff /= alpha ** (k + 1)
            p1s = p1 + j - k
            if 0 <= p1s < Bl_2a.shape[1]:
                corr2 += coeff * Bl_2a[l, p1s]
    I2 -= corr2

    return I1 + I2


# ============================================================
# Convenience: compute full Hamiltonian with Neumann V_ee
# ============================================================

def compute_hamiltonian_neumann(
    basis: list,
    R: float,
    grids: dict,
    Z_A: float = 1.0,
    Z_B: float = 1.0,
    l_max: int = 20,
    verbose: bool = False,
) -> np.ndarray:
    """Compute H = T + V_ne (numerical) + V_ee (Neumann).

    T and V_ne use the existing IBP quadrature (exact for polynomials).
    V_ee uses the Neumann expansion (algebraic, exact within truncation).

    Parameters
    ----------
    basis : list of HylleraasBasisFunction
    R : float
        Internuclear distance.
    grids : dict
        Quadrature grids (for T + V_ne only).
    Z_A, Z_B : float
        Nuclear charges.
    l_max : int
        Neumann series truncation.
    verbose : bool

    Returns
    -------
    H : ndarray, shape (n_basis, n_basis)
    """
    from geovac.hylleraas import compute_hamiltonian_matrix

    # T + V_ne via Numba kernel (fast, no V_ee)
    H_tv = compute_hamiltonian_matrix(
        basis, R, grids, Z_A, Z_B,
        use_numba=True, include_vee=False,
    )

    # Exact V_ee via Neumann expansion
    V_ee_neumann = compute_vee_matrix_neumann(basis, R, l_max, verbose)

    if verbose:
        print(f"  H(T+V_ne) range: [{H_tv.min():.6f}, {H_tv.max():.6f}]")
        print(f"  V_ee(Neumann) range: [{V_ee_neumann.min():.6f}, {V_ee_neumann.max():.6f}]")

    return H_tv + V_ee_neumann


def _compute_vee_numerical_p0(
    basis: list, R: float, grids: dict,
) -> np.ndarray:
    """Compute V_ee matrix via numerical elliptic integral quadrature.

    Vectorized over eta2 dimension for speed.
    """
    from geovac.hylleraas import evaluate_basis_function
    from scipy.special import ellipk

    n_bf = len(basis)
    xi = grids['xi']
    eta = grids['eta']
    w_xi = grids['w_xi']
    w_eta = grids['w_eta']
    half_R = R / 2.0
    n_xi = len(xi)
    n_eta = len(eta)

    # Precompute basis function values at all grid points
    # Shape: (n_bf, n_xi, n_xi, n_eta, n_eta)
    # Too much memory for large grids. Instead, precompute 1D factors.
    # For p=0 basis: phi = e^{-alpha*(xi1+xi2)} * xi1^j * xi2^k * eta1^l * eta2^m
    # So fi(xi1,eta1,xi2,eta2) = f1(xi1,eta1) * f2(xi2,eta2) ... not quite,
    # because of symmetrization.

    # Precompute V_ee kernel on the grid: K(k)/sqrt(s) for all (xi1,xi2,eta1,eta2)
    # V_ee depends only on geometry, not basis. Precompute once.
    rho_arr = half_R * np.sqrt(np.maximum((xi[:, None]**2 - 1) * (1 - eta[None, :]**2), 0.0))
    z_arr = half_R * xi[:, None] * eta[None, :]  # (n_xi, n_eta)

    Vee = np.zeros((n_bf, n_bf))

    for i in range(n_bf):
        for j in range(i, n_bf):
            val = 0.0
            for a in range(n_xi):
                xi1 = xi[a]
                rho1 = rho_arr[a, :]  # (n_eta,) but we need scalar for each eta1
                z1_arr = z_arr[a, :]  # (n_eta,)
                J1_arr = half_R**3 * (xi1**2 - eta**2)

                for c in range(n_xi):
                    xi2 = xi[c]
                    J2_arr = half_R**3 * (xi2**2 - eta**2)

                    for b in range(n_eta):
                        eta1 = eta[b]
                        J1 = J1_arr[b]
                        r1 = rho_arr[a, b]  # scalar
                        zz1 = z_arr[a, b]  # scalar

                        eta2 = eta  # (n_eta,)
                        r2 = rho_arr[c, :]  # (n_eta,)
                        zz2 = z_arr[c, :]  # (n_eta,)

                        fi = evaluate_basis_function(basis[i], xi1, eta1, xi2, eta2)
                        fj = evaluate_basis_function(basis[j], xi1, eta1, xi2, eta2)

                        dz = zz1 - zz2
                        s = (r1 + r2)**2 + dz**2
                        with np.errstate(divide='ignore', invalid='ignore'):
                            k2 = np.where(s > 1e-30, 4 * r1 * r2 / s, 0.0)
                            k2 = np.clip(k2, 0, 1 - 1e-15)
                        s_safe = np.maximum(s, 1e-30)
                        K_vals = ellipk(k2)
                        vee = np.where(s > 1e-30, K_vals / np.sqrt(s_safe), 0.0)

                        integrand = fi * fj * vee * J1 * J2_arr
                        vee_sum = np.sum(integrand * w_eta[b] * w_eta)
                        val += w_xi[a] * w_xi[c] * 8 * np.pi * vee_sum

            Vee[i, j] = val
            Vee[j, i] = val

    return Vee
