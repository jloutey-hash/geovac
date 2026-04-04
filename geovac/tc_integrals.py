"""
Transcorrelated (TC) Two-Body Integrals for Composed Qubit Pipeline
====================================================================

Computes TC-modified two-body integrals that replace the standard ERI
in the composed qubit Hamiltonian.  The TC transformation with Jastrow
factor J = (1/2) r_12 removes the 1/r_12 singularity:

    V_ee + (-1/2 sum_i nabla_i^2 J) = 1/r_12 - 1/r_12 = 0

The remaining gradient operator G = -(1/2) d/dr_12 is a smooth,
first-order differential operator whose matrix elements converge
exponentially in l_max (vs algebraically for the cusp).

The TC Hamiltonian in second quantization:
    H_TC = sum h1_pq a+_p a_q + (1/2) sum g^TC_pqrs a+_p a+_q a_s a_r - 1/4

where h1 is UNCHANGED and g^TC replaces the standard ERI.
The constant -1/4 arises from (1/2)|nabla J|^2 = (1/2)(1/4 + 1/4) = 1/4,
entering with a minus sign from the BCH second-order term.

The gradient operator G = -(1/2) r_hat_12 . (nabla_1 - nabla_2) decomposes
into radial and angular gradient contributions:
  - Radial:  (r_hat_12 . r_hat_i) * dR/dr * Y_{lm}
  - Angular: (R/r) * (r_hat_12)_perp . nabla_Omega Y_{lm}

The angular gradient uses the identity:
  sin(theta_12) * (e_12 . nabla_Omega_1 Y_{lc,mc})
  = -(4pi/3) sum_q [ (lc+1) C_{lc-1} Y_{lc-1,mc-q} - lc C_{lc+1} Y_{lc+1,mc-q} ] Y_{1q}
where C_L are Gaunt coefficients, derived via the Laplacian identity
nabla_Omega f . nabla_Omega g = (1/2)[Delta(fg) - f Delta g - g Delta f].

References:
    - docs/tc_integrals_derivation.md (full derivation)
    - docs/tc_geovac_scoping.md (feasibility analysis)
    - Ten-no (2004), Luo (2010), Dobrautz et al. (2019) for TC methods

Author: GeoVac Development Team
Date: April 2026
"""

from math import factorial, sqrt
from typing import Dict, List, Tuple

import numpy as np
from scipy.special import genlaguerre


# ---------------------------------------------------------------------------
# Radial wavefunction and derivative
# ---------------------------------------------------------------------------

def _radial_wf_grid(
    Z: float, n: int, l: int, r_grid: np.ndarray,
) -> np.ndarray:
    """Normalized hydrogenic radial wavefunction R_{nl}(r) on grid."""
    rho = 2.0 * Z * r_grid / n
    L_poly = genlaguerre(n - l - 1, 2 * l + 1)(rho)
    wf = rho ** l * np.exp(-rho / 2.0) * L_poly
    norm_sq = np.trapezoid(wf ** 2 * r_grid ** 2, r_grid)
    if norm_sq < 1e-30:
        return np.zeros_like(r_grid)
    return wf / np.sqrt(norm_sq)


def _radial_wf_derivative_grid(
    Z: float, n: int, l: int, r_grid: np.ndarray,
) -> np.ndarray:
    """
    Derivative dR_{nl}/dr on the grid, computed analytically.

    Uses the product rule on R_{nl}(r) = N * rho^l * exp(-rho/2) * L_n^{2l+1}(rho)
    where rho = 2Zr/n, so d/dr = (2Z/n) d/drho.
    """
    zeta = 2.0 * Z / n  # rho = zeta * r
    rho = zeta * r_grid

    alpha_L = 2 * l + 1
    nr = n - l - 1  # radial quantum number

    L_poly = genlaguerre(nr, alpha_L)(rho)

    # Derivative of L_nr^{alpha}(rho) w.r.t. rho:
    # d/drho L_nr^alpha = -L_{nr-1}^{alpha+1}(rho) for nr >= 1, else 0
    if nr >= 1:
        dL_poly = -genlaguerre(nr - 1, alpha_L + 1)(rho)
    else:
        dL_poly = np.zeros_like(rho)

    # R_{nl} = N * rho^l * exp(-rho/2) * L(rho)
    # dR/dr = N * zeta * d/drho [rho^l * exp(-rho/2) * L(rho)]
    # d/drho = l*rho^{l-1}*exp(-rho/2)*L - (1/2)*rho^l*exp(-rho/2)*L + rho^l*exp(-rho/2)*dL

    envelope = np.exp(-rho / 2.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        if l > 0:
            term1 = l * rho ** (l - 1) * envelope * L_poly
        else:
            term1 = np.zeros_like(rho)
        term2 = -0.5 * rho ** l * envelope * L_poly
        term3 = rho ** l * envelope * dL_poly

    deriv_rho = term1 + term2 + term3  # d/drho of the unnormalized wf
    deriv_r = zeta * deriv_rho  # chain rule: dR/dr = zeta * dR/drho

    # Normalize using same normalization as _radial_wf_grid
    wf_unnorm = rho ** l * envelope * L_poly
    norm_sq = np.trapezoid(wf_unnorm ** 2 * r_grid ** 2, r_grid)
    if norm_sq < 1e-30:
        return np.zeros_like(r_grid)
    return deriv_r / np.sqrt(norm_sq)


# ---------------------------------------------------------------------------
# Direction cosine multipole kernels
# ---------------------------------------------------------------------------

def _direction_cosine_kernel_K_vectorized(
    r1: np.ndarray, r2: np.ndarray, K: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Vectorized direction cosine kernels F^{(1)}_K and F^{(2)}_K.

    F^{(1)}_K encodes the K-th Legendre component of (r_hat_12 . r_hat_1).
    F^{(2)}_K encodes the K-th Legendre component of (r_hat_12 . r_hat_2).

    r1: shape (N,) broadcast to (N,1)
    r2: shape (M,) broadcast to (1,M)
    Returns F1_K, F2_K each shape (N, M).
    """
    r1_2d = r1[:, np.newaxis]  # (N, 1)
    r2_2d = r2[np.newaxis, :]  # (1, M)

    r_min = np.minimum(r1_2d, r2_2d)
    r_max = np.maximum(r1_2d, r2_2d)

    # Avoid division by zero
    r_max_safe = np.where(r_max > 1e-30, r_max, 1.0)

    ratio = r_min / r_max_safe
    neum_K = ratio ** K / r_max_safe  # r_<^K / r_>^{K+1}

    if K >= 1:
        neum_Km1 = r_min ** (K - 1) / r_max_safe ** K
        coeff_Km1 = K / (2 * K - 1)
    else:
        neum_Km1 = np.zeros_like(r_min)
        coeff_Km1 = 0.0

    neum_Kp1 = r_min ** (K + 1) / r_max_safe ** (K + 2)
    coeff_Kp1 = (K + 1) / (2 * K + 3)

    p1_coupling = coeff_Km1 * neum_Km1 + coeff_Kp1 * neum_Kp1

    F1_K = r1_2d * neum_K - r2_2d * p1_coupling
    F2_K = r1_2d * p1_coupling - r2_2d * neum_K

    # Zero out where r_max was tiny
    mask = r_max < 1e-30
    F1_K = np.where(mask, 0.0, F1_K)
    F2_K = np.where(mask, 0.0, F2_K)

    return F1_K, F2_K


# ---------------------------------------------------------------------------
# TC Y^K potential (radial gradient contribution)
# ---------------------------------------------------------------------------

def _compute_tc_yk_potentials(
    R_on_grid: Dict[Tuple[int, int], np.ndarray],
    R_deriv_on_grid: Dict[Tuple[int, int], np.ndarray],
    unique_nl: List[Tuple[int, int]],
    r_grid: np.ndarray,
    K: int,
) -> Tuple[Dict, Dict]:
    """
    Compute TC Y^K potentials for all (n_r, l_r, n_s, l_s) pairs.

    Returns two dicts:
      yk_tc1[(nr,lr,ns,ls)] = array, for the (r_hat_12.r_hat_1) kernel with R_r*R_s
      yk_tc2[(nr,lr,ns,ls)] = array, for the (r_hat_12.r_hat_2) kernel with R_r*R'_s
    """
    dr = r_grid[1] - r_grid[0]

    # Pre-compute the direction cosine kernel matrices (N x N)
    F1_K, F2_K = _direction_cosine_kernel_K_vectorized(r_grid, r_grid, K)

    yk_tc1: Dict[Tuple[int, int, int, int], np.ndarray] = {}
    yk_tc2: Dict[Tuple[int, int, int, int], np.ndarray] = {}

    for nr, lr in unique_nl:
        for ns, ls in unique_nl:
            key = (nr, lr, ns, ls)

            # Density-like products for electron 2
            f_rs = R_on_grid[(nr, lr)] * R_on_grid[(ns, ls)] * r_grid ** 2
            f_rs_deriv = R_on_grid[(nr, lr)] * R_deriv_on_grid[(ns, ls)] * r_grid ** 2

            # Vectorized integration: yk1[i] = sum_j F1_K[i,j] * f_rs[j] * dr
            yk_tc1[key] = F1_K @ f_rs * dr
            yk_tc2[key] = F2_K @ f_rs_deriv * dr

    return yk_tc1, yk_tc2


# ---------------------------------------------------------------------------
# Wigner 3j and Gaunt coefficients
# ---------------------------------------------------------------------------

def _wigner3j(j1: int, j2: int, j3: int,
              m1: int, m2: int, m3: int) -> float:
    """Wigner 3j symbol for integer arguments."""
    if m1 + m2 + m3 != 0:
        return 0.0
    if abs(j1 - j2) > j3 or j3 > j1 + j2:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0

    def _tri(a: int, b: int, c: int) -> float:
        return (factorial(a + b - c) * factorial(a - b + c)
                * factorial(-a + b + c)
                / factorial(a + b + c + 1))

    pre = ((-1) ** (j1 - j2 - m3)
           * np.sqrt(_tri(j1, j2, j3)
                     * factorial(j1 + m1) * factorial(j1 - m1)
                     * factorial(j2 + m2) * factorial(j2 - m2)
                     * factorial(j3 + m3) * factorial(j3 - m3)))

    t_min = max(0, j2 - j3 - m1, j1 - j3 + m2)
    t_max = min(j1 + j2 - j3, j1 - m1, j2 + m2)

    s = 0.0
    for t in range(t_min, t_max + 1):
        s += ((-1) ** t
              / (factorial(t)
                 * factorial(j3 - j2 + m1 + t)
                 * factorial(j3 - j1 - m2 + t)
                 * factorial(j1 + j2 - j3 - t)
                 * factorial(j1 - m1 - t)
                 * factorial(j2 + m2 - t)))
    return pre * s


def _ck_coefficient(la: int, ma: int, lc: int, mc: int, k: int) -> float:
    """Gaunt angular coupling coefficient c^k(l,m,l',m')."""
    q = mc - ma
    pre = ((-1) ** ma * np.sqrt((2 * la + 1) * (2 * lc + 1)))
    w1 = _wigner3j(la, k, lc, 0, 0, 0)
    if abs(w1) < 1e-15:
        return 0.0
    w2 = _wigner3j(la, k, lc, -ma, q, mc)
    return pre * w1 * w2


def _gaunt_integral(l1: int, m1: int, l2: int, m2: int,
                    l3: int, m3: int) -> float:
    """
    Real Gaunt integral: int Y_{l1,m1} Y_{l2,m2} Y*_{l3,m3} dOmega.

    Uses: (-1)^{m3} sqrt((2l1+1)(2l2+1)(2l3+1)/(4pi))
          * 3j(l1,l2,l3;0,0,0) * 3j(l1,l2,l3;m1,m2,-m3).

    Selection rules: m1 + m2 = m3, l1+l2+l3 even, triangle inequality.
    """
    if m1 + m2 != m3:
        return 0.0
    if (l1 + l2 + l3) % 2 != 0:
        return 0.0
    if abs(l1 - l2) > l3 or l3 > l1 + l2:
        return 0.0
    w1 = _wigner3j(l1, l2, l3, 0, 0, 0)
    if abs(w1) < 1e-15:
        return 0.0
    w2 = _wigner3j(l1, l2, l3, m1, m2, -m3)
    return ((-1) ** m3
            * sqrt((2 * l1 + 1) * (2 * l2 + 1) * (2 * l3 + 1) / (4 * np.pi))
            * w1 * w2)


# ---------------------------------------------------------------------------
# Angular gradient coupling coefficients
# ---------------------------------------------------------------------------

def _angular_gradient_coefficients(
    lc: int, mc: int,
) -> List[Tuple[int, int, int, float]]:
    """
    Angular gradient coupling for Y_{lc,mc}.

    Using the identity:
        sin(theta_12) * (e_12 . nabla_Omega_1 Y_{lc,mc})
        = -(4pi/3) sum_q [ (lc+1) C_{lc-1} Y_{lc-1,mc-q}(Omega_1)
                          - lc   C_{lc+1} Y_{lc+1,mc-q}(Omega_1) ] * Y_{1,q}(Omega_2)

    Derived from:
        nabla_Omega Y*_{1,q} . nabla_Omega Y_{lc,mc}
        = (1/2) sum_L C_L [-L(L+1) + 2 + lc(lc+1)] Y_{L,mc-q}

    For L = lc-1: factor = 2(lc+1)/2 = (lc+1)
    For L = lc+1: factor = -2lc/2     = -lc

    The C_L Gaunt coefficient is: int Y_{L,mc-q} Y_{1,q} Y_{lc,mc} dOmega
    (decomposition of Y_{1,q} * Y_{lc,mc} into Y_L basis).

    Returns list of (l_eff, m_eff, q, coefficient) where:
      - l_eff, m_eff = effective orbital quantum numbers on Omega_1
      - q            = Y_{1,q}(Omega_2) index
      - coefficient  = -(4pi/3) * angular_factor * C_L
    """
    if lc == 0:
        return []  # Angular gradient of Y_{00} is zero

    result: List[Tuple[int, int, int, float]] = []
    for q in [-1, 0, 1]:
        m_eff = mc - q

        for L in [lc - 1, lc + 1]:
            if L < 0:
                continue
            if abs(m_eff) > L:
                continue
            # Parity: L + 1 + lc must be even
            if (L + 1 + lc) % 2 != 0:
                continue

            # Gaunt coefficient: int Y_{L,m_eff} Y_{1,q} Y_{lc,mc} dOmega
            # with m_eff + q + (-mc) = mc-q + q - mc = 0 ... wait
            # Actually we need: int Y*_{L,m_eff} Y_{1,q} Y_{lc,mc} dOmega
            # For real spherical harmonics Y*=Y, so this is
            # _gaunt_integral(1, q, lc, mc, L, m_eff)
            # which requires q + mc = m_eff, i.e., m_eff = mc + q.
            #
            # But we defined m_eff = mc - q above!
            # Let me re-derive: the product Y*_{1,q} Y_{lc,mc} decomposes as
            # sum_L C_L Y_{L,...}. For real Y, Y*_{1,q} = Y_{1,q}.
            # So: Y_{1,q} Y_{lc,mc} = sum_L G(1,q; lc,mc; L,q+mc) Y_{L,q+mc}
            # This means the L-component has m = q + mc, not mc - q.
            #
            # Going back to the derivation:
            # -sin(theta_12) e_12 = nabla_Omega_1 cos(theta_12)
            # = (4pi/3) sum_q [nabla_Omega_1 Y*_{1,q}(Omega_1)] Y_{1,q}(Omega_2)
            #
            # For REAL spherical harmonics, Y*_{1,q} = Y_{1,q}, so:
            # nabla_Omega Y_{1,q} . nabla_Omega Y_{lc,mc}
            # = (1/2) sum_L decomp_coeff(L) * [-L(L+1) + 2 + lc(lc+1)] Y_{L,m_L}
            #
            # where Y_{1,q} Y_{lc,mc} = sum_L decomp_coeff(L) Y_{L, q+mc}
            # So m_L = q + mc.
            m_eff_correct = q + mc
            if abs(m_eff_correct) > L:
                continue

            # Gaunt: int Y*_{L,m_eff_correct} Y_{1,q} Y_{lc,mc} dOmega
            C_L = _gaunt_integral(1, q, lc, mc, L, m_eff_correct)
            if abs(C_L) < 1e-15:
                continue

            # Angular gradient factor from Laplacian identity
            if L == lc - 1:
                ang_factor = (lc + 1)
            else:  # L == lc + 1
                ang_factor = -lc

            # Total: -(4pi/3) * ang_factor * C_L
            total = -(4 * np.pi / 3) * ang_factor * C_L

            if abs(total) > 1e-15:
                result.append((L, m_eff_correct, q, total))

    return result


def _neumann_kernel_K_vectorized(
    r1: np.ndarray, r2: np.ndarray, K: int,
) -> np.ndarray:
    """
    Standard Neumann kernel: r_<^K / r_>^{K+1}.

    Returns shape (len(r1), len(r2)).
    """
    r1_2d = r1[:, np.newaxis]
    r2_2d = r2[np.newaxis, :]
    r_min = np.minimum(r1_2d, r2_2d)
    r_max = np.maximum(r1_2d, r2_2d)
    r_max_safe = np.where(r_max > 1e-30, r_max, 1.0)
    ratio = r_min / r_max_safe
    result = ratio ** K / r_max_safe
    mask = r_max < 1e-30
    return np.where(mask, 0.0, result)


# ---------------------------------------------------------------------------
# Main entry: compute TC two-body integrals for a block
# ---------------------------------------------------------------------------

def compute_tc_integrals_block(
    Z: float,
    states: List[Tuple[int, int, int]],
    n_grid: int = 2000,
    include_angular: bool = True,
) -> Dict[Tuple[int, int, int, int], float]:
    """
    Compute TC two-body integrals g^TC_{pqrs} for a single-center block.

    These replace the standard ERI in the composed qubit Hamiltonian.
    Implements BOTH radial and angular gradient contributions.
    The angular gradient is exact for all orbital angular momenta.

    The TC transformation with J = (1/2) r_12:
    - Cancels V_ee completely (V_ee + nabla^2 J = 0)
    - Replaces it with G = -(1/2) d/dr_12 (smooth, non-Hermitian)
    - Adds -1/4 constant energy shift per electron pair (handled by caller)

    Parameters
    ----------
    Z : float
        Nuclear charge for the block.
    states : list of (n, l, m)
        Spatial orbital quantum numbers.
    n_grid : int
        Radial grid points.
    include_angular : bool
        Include angular gradient contribution (default True).
        Set False to reproduce radial-only results from BX-3.

    Returns
    -------
    dict mapping (p, q, r, s) -> float
        TC two-body integrals in PHYSICIST notation <pr|G|qs>.
        Spatial orbital indices relative to the block.
    """
    unique_nl = sorted(set((n, l) for n, l, m in states))
    n_sp = len(states)

    r_max = 80.0 / max(Z, 0.5)
    r_grid = np.linspace(0, r_max, n_grid + 1)[1:]
    dr = r_grid[1] - r_grid[0]

    # Pre-compute radial wavefunctions and their derivatives
    R_on_grid: Dict[Tuple[int, int], np.ndarray] = {}
    R_deriv_on_grid: Dict[Tuple[int, int], np.ndarray] = {}
    # R/r for angular gradient
    R_over_r_on_grid: Dict[Tuple[int, int], np.ndarray] = {}
    for n, l in unique_nl:
        R_on_grid[(n, l)] = _radial_wf_grid(Z, n, l, r_grid)
        R_deriv_on_grid[(n, l)] = _radial_wf_derivative_grid(Z, n, l, r_grid)
        with np.errstate(divide='ignore', invalid='ignore'):
            ror = np.where(r_grid > 1e-30, R_on_grid[(n, l)] / r_grid, 0.0)
        R_over_r_on_grid[(n, l)] = ror

    # Determine which K values are needed (from Gaunt selection rules)
    max_l = max(l for _, l in unique_nl)
    K_max = 2 * max_l  # Same as standard V_ee: triangle inequality

    # Pre-compute standard Gaunt coefficient table
    ck_table: Dict[Tuple[int, int, int], float] = {}
    for a in range(n_sp):
        la, ma = states[a][1], states[a][2]
        for c in range(n_sp):
            lc, mc = states[c][1], states[c][2]
            for K in range(0, K_max + 1):
                if (la + lc + K) % 2 != 0:
                    continue
                val = _ck_coefficient(la, ma, lc, mc, K)
                if abs(val) > 1e-15:
                    ck_table[(a, c, K)] = val

    # ===================================================================
    # PART 1: Radial gradient contribution (existing)
    # ===================================================================
    tc_eri: Dict[Tuple[int, int, int, int], float] = {}

    for K in range(0, K_max + 1):
        has_K = any(k == K for (_, _, k) in ck_table.keys())
        if not has_K:
            continue

        yk_tc1, yk_tc2 = _compute_tc_yk_potentials(
            R_on_grid, R_deriv_on_grid, unique_nl, r_grid, K
        )

        ac_K_map: Dict[Tuple[int, int], float] = {}
        for (a, c, k), val in ck_table.items():
            if k == K:
                ac_K_map[(a, c)] = val

        for (a, c), ck_ac in ac_K_map.items():
            na, la, ma = states[a]
            nc, lc, mc = states[c]
            for (b, d), ck_bd in ac_K_map.items():
                nb, lb, mb = states[b]
                nd, ld, md = states[d]

                # Term 1: electron 1 radial gradient
                yk1_key = (nb, lb, nd, ld)
                if yk1_key in yk_tc1:
                    integrand1 = (R_on_grid[(na, la)]
                                  * R_deriv_on_grid[(nc, lc)]
                                  * yk_tc1[yk1_key]
                                  * r_grid ** 2)
                    val1 = -0.5 * ck_ac * ck_bd * np.trapezoid(integrand1, r_grid)
                else:
                    val1 = 0.0

                # Term 2: electron 2 radial gradient
                yk2_key = (nb, lb, nd, ld)
                if yk2_key in yk_tc2:
                    integrand2 = (R_on_grid[(na, la)]
                                  * R_on_grid[(nc, lc)]
                                  * yk_tc2[yk2_key]
                                  * r_grid ** 2)
                    val2 = 0.5 * ck_ac * ck_bd * np.trapezoid(integrand2, r_grid)
                else:
                    val2 = 0.0

                total = val1 + val2
                if abs(total) > 1e-15:
                    key = (a, b, c, d)
                    tc_eri[key] = tc_eri.get(key, 0.0) + total

    # ===================================================================
    # PART 2: Angular gradient contribution
    # ===================================================================
    if not include_angular:
        return tc_eri
    # The angular gradient of orbital c produces effective orbitals at
    # (lc+1, mc+q) or (lc-1, mc+q), coupled with Y_{1,q} on the Omega_2 side.
    #
    # The full matrix element for the angular gradient (electron 1 term):
    # -(1/2) * int R_a(r1) [R_c(r1)/r1] * [r2/(r1 * r12)] *
    #          [sin(theta_12) * e_12.grad_Omega Y_{lc,mc}] *
    #          R_b(r2) Y_{ld,md}(Omega_2) * r1^2 dr1 * r2^2 dr2 * dOmega_1 * dOmega_2
    #
    # The angular part factorizes (see _angular_gradient_coefficients):
    # On Omega_1: Y*_{la,ma} * Y_{L,mc+q}  -> delta(la,L) delta(ma, mc+q)
    # On Omega_2: Y*_{lb,mb} * Y_{K,Q} * Y_{1,q} * Y_{ld,md}
    #           -> double Gaunt integral
    #
    # The radial part uses the standard Neumann kernel multiplied by r2:
    # r2 * r_<^K / r_>^{K+1} (from r2/r12 in the transverse projection)

    # Pre-compute angular gradient coefficients for all orbital lc values
    ang_grad_cache: Dict[Tuple[int, int], List] = {}
    for _, lc, mc in states:
        if (lc, mc) not in ang_grad_cache:
            ang_grad_cache[(lc, mc)] = _angular_gradient_coefficients(lc, mc)

    # K_max for angular gradient: may need higher K due to additional coupling
    K_max_ang = K_max + 1  # angular gradient introduces K=1 coupling

    # Pre-compute Neumann kernels for each K value
    neumann_cache: Dict[int, np.ndarray] = {}
    for K in range(0, K_max_ang + 1):
        neumann_cache[K] = _neumann_kernel_K_vectorized(r_grid, r_grid, K)

    # Build state lookup: map (l, m) -> list of state indices
    lm_to_idx: Dict[Tuple[int, int], List[int]] = {}
    for idx, (n, l, m) in enumerate(states):
        key = (l, m)
        if key not in lm_to_idx:
            lm_to_idx[key] = []
        lm_to_idx[key].append(idx)

    # For each orbital c that has angular gradient terms:
    for c_idx in range(n_sp):
        nc, lc, mc = states[c_idx]
        ang_terms = ang_grad_cache.get((lc, mc), [])
        if not ang_terms:
            continue

        for (L_eff, m_eff, q_val, ang_coeff) in ang_terms:
            # Omega_1 integral: int Y*_{la,ma} Y_{L_eff,m_eff} dOmega = delta(la,L_eff) delta(ma,m_eff)
            # So only states a with (la, ma) = (L_eff, m_eff) contribute
            a_indices = lm_to_idx.get((L_eff, m_eff), [])
            if not a_indices:
                continue

            for a_idx in a_indices:
                na, la, ma = states[a_idx]

                # Now handle the Omega_2 side and the Neumann expansion.
                # The Omega_2 integral involves:
                # int Y*_{lb,mb} Y_{K,Q} Y_{1,q_val} Y_{ld,md} dOmega_2
                #
                # Y_{K,Q} comes from the Neumann expansion of 1/r12.
                # Y_{1,q_val} comes from the angular gradient decomposition.
                #
                # The product Y_{K,Q} Y_{1,q_val} = sum_J Gaunt(K,Q; 1,q_val; J,Q+q_val) Y_{J,Q+q_val}
                # Then: int Y*_{lb,mb} Y_{J,Q+q_val} Y_{ld,md} dOmega
                #     = Gaunt(J, Q+q_val; ld, md; lb, mb)  [nonzero if Q+q_val+md=mb]
                #
                # But the Neumann addition theorem gives:
                # 1/r12 = sum_K (4pi/(2K+1)) sum_Q Y*_{KQ}(Omega_1) Y_{KQ}(Omega_2) * neum_K(r1,r2)
                #
                # The Omega_1 integral is already evaluated as delta(la, L_eff)
                # from the angular gradient factorization. The Y*_{KQ}(Omega_1)
                # is SEPARATE from the angular gradient Y_{L_eff}(Omega_1).
                #
                # Wait -- I need to reconsider the structure. The full integral is:
                #
                # int dOmega_1 dOmega_2  Y*_{la,ma}(O1) * [angular_grad_factor(O1,O2)]
                #                        * Y*_{KQ}(O1) Y_{KQ}(O2) * Y*_{lb,mb}(O2) Y_{ld,md}(O2)
                #
                # The angular_grad_factor puts Y_{L,m_eff}(O1) * Y_{1,q}(O2).
                # So Omega_1 integral: int Y*_{la,ma} Y_{L,m_eff} Y*_{KQ} dOmega_1
                #                    = Gaunt(la, L, K; ma, m_eff, -Q) type integral
                #                    (product of three Y's)
                #
                # This is NOT simply a delta function! I was wrong earlier.
                # The Y*_{KQ} from the Neumann expansion also lives on Omega_1.
                #
                # So the Omega_1 integral is a triple product:
                # int Y*_{la,ma} Y_{L_eff,m_eff} Y*_{K,Q} dOmega_1
                # = _gaunt_integral(L_eff, m_eff, la, ma, K, Q)
                # where Q = m_eff - ma (from m conservation)
                #
                # And the Omega_2 integral is also a triple product:
                # int Y_{K,Q} Y_{1,q_val} Y*_{lb,mb} Y_{ld,md} dOmega_2
                # = sum_J [Gaunt(K,Q; 1,q; J,Q+q)] * [Gaunt(J,Q+q; ld,md; lb,mb)]
                #
                # Hmm, the Omega_2 has FOUR Y's. Let me handle this properly.

                # The Omega_2 integral:
                # int Y_{K,Q}(O2) * Y_{1,q_val}(O2) * Y*_{lb,mb}(O2) * Y_{ld,md}(O2) dO2
                #
                # Using completeness: Y_K * Y_1 = sum_J C^J Y_J, then
                # int Y_J * Y*_lb * Y_ld = Gaunt.
                # But this is expensive for a sum over all K and Q.
                #
                # Alternatively: int Y_K * Y_1 * Y*_lb * Y_ld dO
                # = sum_J Gaunt(K,1,J) * Gaunt(J,lb,ld)
                # where the Gaunt integrals enforce selection rules.

                # For each K value in the Neumann expansion:
                for K in range(0, K_max_ang + 1):
                    Q_val = m_eff - ma  # from Omega_1 integral

                    # Check if Q is valid for Y_{K,Q}
                    if abs(Q_val) > K:
                        continue

                    # Omega_1 triple integral:
                    # int Y*_{la,ma}(O1) Y_{L_eff,m_eff}(O1) Y*_{K,Q}(O1) dO1
                    # = int Y_{la,ma} Y_{L_eff,m_eff} Y_{K,-Q} dO1  (for real Y)
                    # Wait, Y*_{K,Q} = Y_{K,Q} for real spherical harmonics.
                    # But the Neumann theorem uses Y*_{KQ}(O1) Y_{KQ}(O2).
                    # For real Y, Y*_{KQ} = Y_{KQ}, so:
                    # Omega_1 integral = int Y_{la,ma} Y_{L_eff,m_eff} Y_{K,Q} dO1
                    #
                    # For this to be nonzero: ma + m_eff + Q = ... no, the Gaunt
                    # integral int Y_a Y_b Y_c dO requires specific selection.
                    #
                    # Using _gaunt_integral(la, ma, L_eff, m_eff, K, Q_val)?
                    # _gaunt_integral(l1,m1,l2,m2,l3,m3) = int Y_{l1} Y_{l2} Y*_{l3} dO
                    # with selection m1+m2=m3.
                    #
                    # But I want int Y_{la,ma} Y_{L,m_eff} Y_{K,Q} dO
                    # = _gaunt_integral(la, ma, L_eff, m_eff, K, -Q_val)?
                    # No, that would require ma + m_eff = -Q_val.
                    #
                    # Actually, the three-Y integral for REAL spherical harmonics:
                    # int Y_{l1,m1} Y_{l2,m2} Y_{l3,m3} dO
                    # This is NOT the same as _gaunt_integral which has Y*_{l3}.
                    #
                    # For real Y: Y*_{l,m} = Y_{l,m}, so _gaunt_integral gives:
                    # int Y_{l1,m1} Y_{l2,m2} Y_{l3,m3} dO with m1+m2=m3
                    #
                    # I need: int Y_{la,ma} Y_{L,m_eff} Y_{K,Q} dO
                    # = _gaunt_integral(la, ma, L_eff, m_eff, K, Q_val)
                    # if ma + m_eff = Q_val.
                    #
                    # But Q_val = m_eff - ma, so ma + m_eff = 2*m_eff - Q_val.
                    # That doesn't match Q_val unless ma = 0 and m_eff = Q_val.
                    #
                    # Let me reconsider: for the Gaunt integral
                    # int Y_a Y_b Y*_c dO with m_a + m_b = m_c,
                    # if I want int Y_a Y_b Y_c dO, I use:
                    # int Y_a Y_b Y_c dO = int Y_a Y_b Y*_c dO (since Y_c = Y*_c for real)
                    # This requires m_a + m_b = m_c.
                    #
                    # So: int Y_{la,ma} Y_{L,m_eff} Y_{K,Q} dO
                    # requires ma + m_eff = Q.
                    # Since Q = m_eff - ma, this gives ma + m_eff = m_eff - ma -> 2*ma = 0 -> ma = 0.
                    # That can't be right in general!

                    # I think the issue is in my definition of Q_val.
                    # Let me reconsider.
                    #
                    # The Neumann expansion: 1/r12 = sum_K sum_Q (4pi/(2K+1)) Y*_KQ(O1) Y_KQ(O2) neum_K(r1,r2)
                    # For real Y: Y*_KQ = Y_KQ.
                    #
                    # The angular gradient gives: ang_factor * Y_{L,m_eff}(O1) * Y_{1,q}(O2)
                    #
                    # Combined Omega_1 integral:
                    # int Y*_{la,ma}(O1) Y_{L,m_eff}(O1) Y_{K,Q}(O1) dO1
                    # = _gaunt_integral(L, m_eff, K, Q, la, ma)
                    # requires m_eff + Q = ma.
                    # So Q = ma - m_eff.

                    Q_val = ma - m_eff  # CORRECTED

                    if abs(Q_val) > K:
                        continue

                    # Parity: la + L_eff + K must be even
                    if (la + L_eff + K) % 2 != 0:
                        continue
                    # Triangle: |la - L_eff| <= K <= la + L_eff
                    if K < abs(la - L_eff) or K > la + L_eff:
                        continue

                    gaunt_O1 = _gaunt_integral(L_eff, m_eff, K, Q_val, la, ma)
                    if abs(gaunt_O1) < 1e-15:
                        continue

                    # Omega_2 integral:
                    # int Y*_{lb,mb}(O2) Y_{K,Q}(O2) Y_{1,q_val}(O2) Y_{ld,md}(O2) dO2
                    #
                    # Decompose: Y_{K,Q} Y_{1,q_val} = sum_J Gaunt(...) Y_{J, Q+q_val}
                    # Then: int Y*_{lb,mb} Y_{J,Q+q_val} Y_{ld,md} dO2 = _gaunt_integral(J,Q+q_val, ld,md, lb,mb)
                    #   requires Q+q_val+md = mb

                    m_J = Q_val + q_val
                    mb_needed = m_J  # from the Gaunt on Omega_2, we need md + m_J = mb
                    # Wait: the second Gaunt: int Y*_{lb,mb} Y_{J,m_J} Y_{ld,md} dO2
                    # = _gaunt_integral(J, m_J, ld, md, lb, mb) requires m_J + md = mb

                    # Only states b with mb = m_J + md can contribute.
                    # We need to loop over valid (b,d) pairs.

                    for d_idx in range(n_sp):
                        nd, ld, md = states[d_idx]
                        mb_req = m_J + md

                        b_indices = lm_to_idx.get((None, mb_req), None)
                        # Actually, we need to loop over all b with mb = mb_req
                        for b_idx in range(n_sp):
                            nb, lb, mb = states[b_idx]
                            if mb != mb_req:
                                continue

                            # Compute Omega_2 integral as sum over J
                            gaunt_O2_total = 0.0
                            for J in range(abs(K - 1), K + 2):  # J = K-1, K, K+1
                                if J < 0:
                                    continue
                                if abs(m_J) > J:
                                    continue
                                # Gaunt 1: int Y*_{J,m_J} Y_{K,Q} Y_{1,q_val} dO
                                # = _gaunt_integral(K, Q_val, 1, q_val, J, m_J)
                                if (K + 1 + J) % 2 != 0:
                                    continue
                                if J < abs(K - 1) or J > K + 1:
                                    continue
                                g1 = _gaunt_integral(K, Q_val, 1, q_val, J, m_J)
                                if abs(g1) < 1e-15:
                                    continue

                                # Gaunt 2: int Y*_{lb,mb} Y_{J,m_J} Y_{ld,md} dO
                                # = _gaunt_integral(J, m_J, ld, md, lb, mb)
                                if (J + ld + lb) % 2 != 0:
                                    continue
                                if lb < abs(J - ld) or lb > J + ld:
                                    continue
                                g2 = _gaunt_integral(J, m_J, ld, md, lb, mb)
                                if abs(g2) < 1e-15:
                                    continue

                                gaunt_O2_total += g1 * g2

                            if abs(gaunt_O2_total) < 1e-15:
                                continue

                            # Radial integral for angular gradient (electron 1):
                            # -(1/2) * R_a(r1) * [R_c(r1)/r1] * r1^2 * dr1
                            #        * r2 * neum_K(r1,r2) * R_b(r2) * R_d(r2) * r2^2 * dr2
                            #        * (4pi/(2K+1))
                            #        * ang_coeff * gaunt_O1 * gaunt_O2_total
                            #
                            # The factor R_c/r1 * r1^2 = R_c * r1.
                            # The factor r2 * neum_K * R_b * R_d * r2^2 = yk_ang
                            #
                            # yk_ang[i] = sum_j r2[j] * neum_K[i,j] * R_b(r2[j]) * R_d(r2[j]) * r2[j]^2 * dr

                            neum_K_mat = neumann_cache[K]
                            density_2 = (r_grid * R_on_grid[(nb, lb)]
                                         * R_on_grid[(nd, ld)] * r_grid ** 2)
                            yk_ang = neum_K_mat @ density_2 * dr

                            # r1 integral: R_a(r1) * R_c(r1)/r1 * yk_ang(r1) * r1^2 dr1
                            # = R_a * (R_c/r) * yk_ang * r1^2
                            integrand = (R_on_grid[(na, la)]
                                         * R_over_r_on_grid[(nc, lc)]
                                         * yk_ang
                                         * r_grid ** 2)
                            rad_integral = np.trapezoid(integrand, r_grid)

                            # Prefactor: -(1/2) * (4pi/(2K+1)) * ang_coeff * gaunt_O1 * gaunt_O2
                            prefactor = (-0.5 * (4 * np.pi / (2 * K + 1))
                                         * ang_coeff * gaunt_O1 * gaunt_O2_total)
                            val_ang = prefactor * rad_integral

                            if abs(val_ang) > 1e-15:
                                key = (a_idx, b_idx, c_idx, d_idx)
                                tc_eri[key] = tc_eri.get(key, 0.0) + val_ang

    # Similarly for electron 2 angular gradient (swap roles of electrons 1 and 2)
    # By the antisymmetry structure of G = -(1/2)(r_hat_12.nabla_1 - r_hat_12.nabla_2),
    # the electron 2 angular gradient has the opposite sign and swapped electron indices.
    #
    # For the electron 2 term: +(1/2) r_hat_12 . nabla_2 phi_d
    # The angular gradient of phi_d produces effective orbitals at (ld+1,md+q) or (ld-1,md+q)
    # coupled with Y_{1,q} on the Omega_1 side (swapped from electron 1).
    #
    # By the structure of the two-body integral, the electron 2 angular gradient
    # gives the same contribution with (c,d) swapped and sign flipped.
    # This means: val_ang_2(a,b,c,d) = -val_ang_1(b,a,d,c) ... but the exact
    # relation depends on the index conventions.
    #
    # For a complete and correct implementation, we compute the electron 2 term
    # by swapping the roles: the angular gradient acts on orbital d (electron 2),
    # and the Y_{1,q} factor appears on the Omega_1 side.

    # Pre-compute angular gradient for orbital d
    for d_idx in range(n_sp):
        nd, ld, md = states[d_idx]
        ang_terms = ang_grad_cache.get((ld, md), [])
        if not ang_terms:
            continue

        for (L_eff, m_eff, q_val, ang_coeff) in ang_terms:
            # For electron 2: the angular gradient gives Y_{L,m_eff}(O2), Y_{1,q}(O1)
            # Omega_2 integral: int Y*_{lb,mb} Y_{L,m_eff} Y_{K,Q} dO2
            b_indices_for_L = lm_to_idx.get((L_eff, m_eff), None)
            # Wait, the Omega_2 integral is:
            # int Y*_{lb,mb}(O2) Y_{L,m_eff}(O2) Y_{K,Q}(O2) dO2
            # which requires mb = m_eff + Q (if Q from Neumann).
            # This is not simply delta(lb, L_eff)!

            # The Omega_1 integral for electron 2 term:
            # int Y*_{la,ma}(O1) Y_{K,Q}(O1) Y_{1,q_val}(O1) Y_{lc,mc}(O1) dO1
            # = four-Y integral on Omega_1

            # This has the same structure as the electron 1 case but with
            # (a <-> b) and (c <-> d) swapped, and the sign flipped.

            for K in range(0, K_max_ang + 1):
                for b_idx in range(n_sp):
                    nb, lb, mb = states[b_idx]
                    Q_val = mb - m_eff
                    if abs(Q_val) > K:
                        continue
                    if (lb + L_eff + K) % 2 != 0:
                        continue
                    if K < abs(lb - L_eff) or K > lb + L_eff:
                        continue

                    gaunt_O2 = _gaunt_integral(L_eff, m_eff, K, Q_val, lb, mb)
                    if abs(gaunt_O2) < 1e-15:
                        continue

                    m_J = Q_val + q_val

                    for a_idx in range(n_sp):
                        na, la, ma = states[a_idx]
                        for c_idx_inner in range(n_sp):
                            nc, lc, mc = states[c_idx_inner]
                            ma_req = m_J + mc
                            if ma != ma_req:
                                continue

                            gaunt_O1_total = 0.0
                            for J in range(abs(K - 1), K + 2):
                                if J < 0:
                                    continue
                                if abs(m_J) > J:
                                    continue
                                if (K + 1 + J) % 2 != 0:
                                    continue
                                if J < abs(K - 1) or J > K + 1:
                                    continue
                                g1 = _gaunt_integral(K, Q_val, 1, q_val, J, m_J)
                                if abs(g1) < 1e-15:
                                    continue
                                if (J + lc + la) % 2 != 0:
                                    continue
                                if la < abs(J - lc) or la > J + lc:
                                    continue
                                g2 = _gaunt_integral(J, m_J, lc, mc, la, ma)
                                if abs(g2) < 1e-15:
                                    continue
                                gaunt_O1_total += g1 * g2

                            if abs(gaunt_O1_total) < 1e-15:
                                continue

                            # Radial integral for electron 2 angular gradient
                            neum_K_mat = neumann_cache[K]
                            density_1 = (r_grid * R_on_grid[(na, la)]
                                         * R_on_grid[(nc, lc)] * r_grid ** 2)
                            yk_ang = neum_K_mat.T @ density_1 * dr
                            # Note: for electron 2, the r factor is r1 not r2

                            integrand = (R_on_grid[(nb, lb)]
                                         * R_over_r_on_grid[(nd, ld)]
                                         * yk_ang
                                         * r_grid ** 2)
                            rad_integral = np.trapezoid(integrand, r_grid)

                            # Sign is +1/2 for electron 2 (opposite to electron 1)
                            prefactor = (+0.5 * (4 * np.pi / (2 * K + 1))
                                         * ang_coeff * gaunt_O2 * gaunt_O1_total)
                            val_ang = prefactor * rad_integral

                            if abs(val_ang) > 1e-15:
                                key = (a_idx, b_idx, c_idx_inner, d_idx)
                                tc_eri[key] = tc_eri.get(key, 0.0) + val_ang

    return tc_eri


def tc_eri_to_chemist(
    tc_eri_phys: Dict[Tuple[int, int, int, int], float],
    n_spatial: int,
) -> np.ndarray:
    """
    Convert TC ERI from physicist notation dict to dense chemist notation array.

    Physicist: <pr|G|qs> -> Chemist: (pq|rs) = <pr|G|qs>

    Parameters
    ----------
    tc_eri_phys : dict
        TC integrals in physicist notation {(p,r,q,s): value}.
    n_spatial : int
        Number of spatial orbitals.

    Returns
    -------
    np.ndarray, shape (n_spatial, n_spatial, n_spatial, n_spatial)
        TC ERI in chemist notation.
    """
    eri = np.zeros((n_spatial, n_spatial, n_spatial, n_spatial))
    for (p, r, q, s), val in tc_eri_phys.items():
        # Physicist <pr|G|qs> = Chemist (pq|rs)
        eri[p, q, r, s] = val
    return eri


# ---------------------------------------------------------------------------
# Diagnostic: compare standard and TC integrals
# ---------------------------------------------------------------------------

def compare_standard_vs_tc(
    Z: float,
    states: List[Tuple[int, int, int]],
    n_grid: int = 2000,
) -> Dict:
    """
    Compute both standard V_ee and TC gradient integrals for comparison.

    Returns a dict with diagnostic information.
    """
    from geovac.composed_qubit import _compute_rk_integrals_block, _build_eri_block

    # Standard ERI
    rk_cache = _compute_rk_integrals_block(Z, states, n_grid)
    std_eri = _build_eri_block(Z, states, rk_cache)

    # TC integrals
    tc_eri = compute_tc_integrals_block(Z, states, n_grid)

    # Count nonzero
    n_std = sum(1 for v in std_eri.values() if abs(v) > 1e-12)
    n_tc = sum(1 for v in tc_eri.values() if abs(v) > 1e-12)

    # Sum of absolute values (proxy for 1-norm contribution)
    sum_std = sum(abs(v) for v in std_eri.values())
    sum_tc = sum(abs(v) for v in tc_eri.values())

    return {
        'n_std_nonzero': n_std,
        'n_tc_nonzero': n_tc,
        'sum_abs_std': sum_std,
        'sum_abs_tc': sum_tc,
        'std_eri': std_eri,
        'tc_eri': tc_eri,
    }


# ---------------------------------------------------------------------------
# TC-modified composed Hamiltonian builder
# ---------------------------------------------------------------------------

def build_tc_composed_hamiltonian(
    spec: 'MolecularSpec',
    pk_in_hamiltonian: bool = True,
    n_grid: int = 1500,
    include_angular: bool = False,
    verbose: bool = False,
) -> Dict:
    """
    Build a TC-modified composed qubit Hamiltonian.

    Replaces standard V_ee (ERI) with TC gradient integrals in each block,
    then assembles the full Hamiltonian via JW transform.

    The TC constant shift (-1/4 per 2-electron block) is added to the
    nuclear repulsion constant.

    Parameters
    ----------
    spec : MolecularSpec
        Molecular specification (same as standard builder).
    pk_in_hamiltonian : bool
        Include PK in h1 (default True).
    n_grid : int
        Radial grid points for TC integral evaluation.
    include_angular : bool
        Include angular gradient for l>0 orbitals (default False).
        Angular gradient adds ~2.7x Pauli terms for <0.02 pp accuracy
        improvement at max_n=2 — a net negative for quantum efficiency.
    verbose : bool
        Print progress.

    Returns
    -------
    dict
        Same structure as build_composed_hamiltonian, plus TC metadata.
    """
    import time
    from geovac.composed_qubit import (
        _enumerate_states, _compute_pk_matrix_elements,
    )
    from geovac.qubit_encoding import build_fermion_op_from_integrals
    from openfermion import jordan_wigner

    t0 = time.perf_counter()

    # ------------------------------------------------------------------
    # 1. Enumerate states per block and compute offsets
    # ------------------------------------------------------------------
    block_info: List[Dict] = []
    offset = 0
    n_pairs = 0  # count 2-electron blocks for constant shift

    for blk in spec.blocks:
        info: Dict = {
            'label': blk.label,
            'block_type': blk.block_type,
            'Z_center': blk.Z_center,
            'n_electrons': blk.n_electrons,
        }
        center_states = _enumerate_states(blk.max_n)
        info['center_states'] = center_states
        info['center_offset'] = offset
        info['center_M'] = len(center_states)
        offset += len(center_states)

        if blk.has_h_partner:
            partner_n = blk.max_n_partner if blk.max_n_partner > 0 else blk.max_n
            partner_states = _enumerate_states(partner_n)
            info['partner_states'] = partner_states
            info['partner_offset'] = offset
            info['partner_M'] = len(partner_states)
            info['Z_partner'] = blk.Z_partner
            offset += len(partner_states)
        else:
            info['partner_states'] = None
            info['partner_offset'] = None
            info['partner_M'] = 0
            info['Z_partner'] = None

        info['pk_A'] = blk.pk_A
        info['pk_B'] = blk.pk_B
        block_info.append(info)

        if blk.n_electrons == 2:
            n_pairs += 1

    M = offset
    Q = 2 * M

    if verbose:
        print(f"[TC build] {spec.name}: M={M}, Q={Q}, {n_pairs} electron pairs")

    # ------------------------------------------------------------------
    # 2. Build h1 (same as standard)
    # ------------------------------------------------------------------
    h1 = np.zeros((M, M))
    for bi in block_info:
        off_c = bi['center_offset']
        Z_c = bi['Z_center']
        for i, (n, l, m) in enumerate(bi['center_states']):
            h1[off_c + i, off_c + i] = -Z_c ** 2 / (2.0 * n ** 2)

        if bi['partner_states'] is not None:
            off_p = bi['partner_offset']
            Z_p = bi['Z_partner']
            for i, (n, l, m) in enumerate(bi['partner_states']):
                h1[off_p + i, off_p + i] = -Z_p ** 2 / (2.0 * n ** 2)

    # ------------------------------------------------------------------
    # 3. PK pseudopotential (same as standard)
    # ------------------------------------------------------------------
    h1_pk_full = np.zeros((M, M))
    for bi in block_info:
        if bi['pk_A'] > 0.0:
            h1_pk_block = _compute_pk_matrix_elements(
                bi['Z_center'], bi['center_states'],
                bi['pk_A'], bi['pk_B'],
            )
            off_c = bi['center_offset']
            M_c = bi['center_M']
            for i in range(M_c):
                for j in range(M_c):
                    if abs(h1_pk_block[i, j]) > 1e-15:
                        h1_pk_full[off_c + i, off_c + j] = h1_pk_block[i, j]

    if pk_in_hamiltonian:
        h1 = h1 + h1_pk_full

    # ------------------------------------------------------------------
    # 4. Build TC ERI (block-diagonal, using TC gradient integrals)
    # ------------------------------------------------------------------
    eri = np.zeros((M, M, M, M))
    tc_cache: Dict[Tuple[float, int], Dict] = {}

    def _get_tc_eri(Z: float, max_n: int) -> Dict:
        key = (Z, max_n)
        if key not in tc_cache:
            st = _enumerate_states(max_n)
            tc_cache[key] = compute_tc_integrals_block(
                Z, st, n_grid, include_angular=include_angular,
            )
        return tc_cache[key]

    for bi in block_info:
        # Center sub-block
        Z_c = bi['Z_center']
        off_c = bi['center_offset']
        max_n_c = max(n for n, l, m in bi['center_states'])
        tc_eri_c = _get_tc_eri(Z_c, max_n_c)

        for (a, b, c, d), val in tc_eri_c.items():
            # physicist <ab|G|cd> -> chemist (ac|bd)
            eri[a + off_c, c + off_c, b + off_c, d + off_c] = val

        # Partner sub-block
        if bi['partner_states'] is not None:
            Z_p = bi['Z_partner']
            off_p = bi['partner_offset']
            max_n_p = max(n for n, l, m in bi['partner_states'])
            tc_eri_p = _get_tc_eri(Z_p, max_n_p)

            for (a, b, c, d), val in tc_eri_p.items():
                eri[a + off_p, c + off_p, b + off_p, d + off_p] = val

    # ------------------------------------------------------------------
    # 5. JW transform with TC constant shift
    # ------------------------------------------------------------------
    tc_constant = -0.25 * n_pairs
    nuclear_repulsion = spec.nuclear_repulsion_constant + tc_constant

    fermion_op = build_fermion_op_from_integrals(h1, eri, nuclear_repulsion)
    qubit_op = jordan_wigner(fermion_op)
    N_pauli = len(qubit_op.terms)

    elapsed = time.perf_counter() - t0

    if verbose:
        print(f"[TC build] {spec.name}: Q={Q}, N_pauli={N_pauli}, "
              f"TC_shift={tc_constant:.3f}, wall_time={elapsed:.1f}s")

    n_eri_total = int(np.count_nonzero(np.abs(eri) > 1e-15))

    return {
        'M': M,
        'Q': Q,
        'N_pauli': N_pauli,
        'nuclear_repulsion': nuclear_repulsion,
        'tc_constant': tc_constant,
        'n_pairs_tc': n_pairs,
        'wall_time_s': elapsed,
        'n_eri_total': n_eri_total,
        'pk_in_hamiltonian': pk_in_hamiltonian,
        'h1': h1,
        'h1_pk': h1_pk_full,
        'eri': eri,
        'qubit_op': qubit_op,
        'fermion_op': fermion_op,
    }
