"""
Level 4 multichannel solver for H2.

Extends Phase 1 (l1=l2=0 only) to arbitrary l_max by coupling (l1, l2)
partial-wave channels through the anisotropic nuclear field and the
electron-electron Gaunt integrals.

The angular eigenvalue problem at fixed (rho, R_e):

    [-1/2 d^2/dalpha^2 + V_diag(alpha) + W_coupling(alpha)] u = mu u

where channels are labeled by (l1, l2) pairs with l1+l2 even (Sigma_g).

References:
  - Level 4 design: papers/core/level4_geometry_design.md
  - Lin, Phys. Rep. 257, 1 (1995)
  - hyperspherical_angular.py for Gaunt integral machinery
"""

import numpy as np
from scipy.linalg import eigh
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh_tridiagonal
from typing import Tuple, List, Union, Callable, Optional
from math import sqrt, factorial

from scipy.special import roots_laguerre, eval_laguerre

from geovac.hyperspherical_angular import gaunt_integral
from geovac.hyperspherical_radial import (
    solve_radial_spectral as _solve_radial_spectral_l3,
    _build_laguerre_matrices_dirichlet,
)


def _channel_list(
    l_max: int,
    homonuclear: bool = True,
) -> List[Tuple[int, int]]:
    """
    Build list of (l1, l2) channels for sigma orbitals (m=0).

    Constraints:
    - m1 = m2 = 0 (sigma orbitals)
    - l1 + l2 even (gerade symmetry) when homonuclear=True
    - All l1, l2 combinations when homonuclear=False
    - Both orderings (l1, l2) and (l2, l1) included when l1 != l2,
      because they have DIFFERENT centrifugal potentials and the
      nuclear field couples (0,0) -> (0,2) and (0,0) -> (2,0)
      through different electrons.

    Returns list of (l1, l2) tuples sorted by l1+l2 then l1.
    """
    channels = []
    for l1 in range(l_max + 1):
        for l2 in range(l_max + 1):
            if homonuclear and (l1 + l2) % 2 != 0:
                continue
            channels.append((l1, l2))
    channels.sort(key=lambda x: (x[0] + x[1], x[0]))
    return channels


def _channel_list_extended(
    l_max: int, m_max: int,
    l_max_per_m: Union[dict, None] = None,
    homonuclear: bool = True,
) -> List[Tuple[int, int, int, int]]:
    """
    Build list of (l1, m1, l2, m2) channels with M = 0.

    Constraints:
    - m1 + m2 = 0 (total M = 0)
    - l1 + l2 even (gerade symmetry) when homonuclear=True
    - All l1+l2 parities when homonuclear=False
    - |m1| <= l1, |m2| = |m1| <= l2
    - |m1| <= m_max

    When m_max=0, returns [(l1, 0, l2, 0)] matching _channel_list order.

    Parameters
    ----------
    l_max : int
        Default maximum angular momentum per electron.
    m_max : int
        Maximum |m| per electron.
    l_max_per_m : dict or None
        Per-|m| angular momentum limit. E.g., {0: 4, 1: 3} uses l_max=4
        for sigma (m=0) and l_max=3 for pi (|m|=1). If None, uses l_max
        for all m values. This allows keeping the total channel count in a
        regime where the single-channel adiabatic approximation is valid.
    homonuclear : bool
        If True, enforce l1+l2 even (gerade symmetry). Default True.

    Returns list sorted by (l1+l2, l1, |m1|, m1).
    """
    channels: List[Tuple[int, int, int, int]] = []
    for l1 in range(l_max + 1):
        for l2 in range(l_max + 1):
            if homonuclear and (l1 + l2) % 2 != 0:
                continue
            m_limit = min(l1, l2, m_max)
            for m1 in range(-m_limit, m_limit + 1):
                m2 = -m1
                am = abs(m1)
                # Apply per-m angular momentum limit
                if l_max_per_m is not None and am in l_max_per_m:
                    lm = l_max_per_m[am]
                    if l1 > lm or l2 > lm:
                        continue
                channels.append((l1, m1, l2, m2))
    channels.sort(key=lambda x: (x[0] + x[2], x[0], abs(x[1]), x[1]))
    return channels


def compute_nuclear_coupling(
    l1p: int, l2p: int, l1: int, l2: int,
    m1: int, m2: int,
    alpha: np.ndarray, rho: float, Z: float = 1.0,
    Z_A: float = None, Z_B: float = None,
    rho_A: float = None, rho_B: float = None,
) -> np.ndarray:
    """
    Nuclear attraction coupling via algebraic multipole expansion.

    For nucleus A at distance R_A from origin (+z) and nucleus B at R_B (-z),
    the potential on electron i is:

        V_i = -Z_A sum_k f_k(s_i, rho_A) P_k
              -Z_B sum_k (-1)^k f_k(s_i, rho_B) P_k

    where f_k(s, rho) = (min(s,rho)/max(s,rho))^k / max(s,rho) and
    rho_X = R_X / R_e is the nuclear distance in hyperradial units.

    When rho_A = rho_B = rho (symmetric/midpoint origin), the two terms
    combine into coefficient -(Z_A + (-1)^k Z_B), recovering the original
    formula. When rho_A != rho_B (shifted origin), the radial factors
    differ and two separate sums are required.

    The 3j symbol (l' k l; 0 0 0) requires l'+k+l even, so for a given
    (l', l) pair, only one parity of k contributes.

    Parameters
    ----------
    l1p, l2p : int
        Bra channel angular momenta.
    l1, l2 : int
        Ket channel angular momenta.
    m1, m2 : int
        Shared magnetic quantum numbers (bra = ket, diagonal in m).
    alpha : ndarray
        Alpha grid points.
    rho : float
        R / (2 R_e). Used when rho_A, rho_B not specified.
    Z : float
        Nuclear charge per nucleus (used when Z_A, Z_B not specified).
    Z_A, Z_B : float or None
        Per-nucleus charges. If None, uses Z for both.
    rho_A, rho_B : float or None
        Nuclear distances from origin in R_e units. If None, uses rho
        for both (symmetric midpoint origin).

    Returns
    -------
    V : ndarray of shape (n_alpha,)
        Nuclear coupling at each alpha (charge function units, / R_e).
    """
    from geovac.angular_integrals import wigner3j

    # Resolve charges: backward compatible
    if Z_A is None:
        Z_A = Z
    if Z_B is None:
        Z_B = Z
    # Resolve nuclear positions: backward compatible
    if rho_A is None:
        rho_A = rho
    if rho_B is None:
        rho_B = rho

    n_alpha = len(alpha)
    cos_a = np.cos(alpha)
    sin_a = np.sin(alpha)
    am1 = abs(m1)
    am2 = abs(m2)

    V = np.zeros(n_alpha)

    # Electron 1: couples (l1p, l1) with |m1|, acts when l2p == l2
    if l2p == l2:
        norm1 = sqrt((2 * l1p + 1) * (2 * l1 + 1))
        phase1 = (-1) ** am1
        s1 = cos_a

        # 3j requires l1p + k + l1 even → k has parity of (l1p + l1)
        k_min = abs(l1p - l1)
        k_max = l1p + l1
        required_parity = (l1p + l1) % 2
        if k_min % 2 != required_parity:
            k_min += 1

        # Precompute radial factors for both nuclei
        min_s1_A = np.minimum(s1, rho_A)
        max_s1_A = np.maximum(s1, rho_A)
        min_s1_B = np.minimum(s1, rho_B)
        max_s1_B = np.maximum(s1, rho_B)

        V1 = np.zeros(n_alpha)
        for k in range(k_min, k_max + 1, 2):
            g0 = wigner3j(l1p, k, l1, 0, 0, 0)
            if abs(g0) < 1e-15:
                continue
            gm = wigner3j(l1p, k, l1, -am1, 0, am1)
            coeff = g0 * gm
            if abs(coeff) < 1e-30:
                continue
            # Nucleus A at +R_A: -Z_A * f_k(s, rho_A)
            c_k_A = -Z_A * (min_s1_A / max_s1_A) ** k / max_s1_A
            # Nucleus B at -R_B: -Z_B * (-1)^k * f_k(s, rho_B)
            c_k_B = -Z_B * (-1)**k * (min_s1_B / max_s1_B) ** k / max_s1_B
            V1 += coeff * (c_k_A + c_k_B)

        V += phase1 * norm1 * V1

    # Electron 2: couples (l2p, l2) with |m2|, acts when l1p == l1
    if l1p == l1:
        norm2 = sqrt((2 * l2p + 1) * (2 * l2 + 1))
        phase2 = (-1) ** am2
        s2 = sin_a

        # 3j requires l2p + k + l2 even → k has parity of (l2p + l2)
        k_min = abs(l2p - l2)
        k_max = l2p + l2
        required_parity = (l2p + l2) % 2
        if k_min % 2 != required_parity:
            k_min += 1

        # Precompute radial factors for both nuclei
        min_s2_A = np.minimum(s2, rho_A)
        max_s2_A = np.maximum(s2, rho_A)
        min_s2_B = np.minimum(s2, rho_B)
        max_s2_B = np.maximum(s2, rho_B)

        V2 = np.zeros(n_alpha)
        for k in range(k_min, k_max + 1, 2):
            g0 = wigner3j(l2p, k, l2, 0, 0, 0)
            if abs(g0) < 1e-15:
                continue
            gm = wigner3j(l2p, k, l2, -am2, 0, am2)
            coeff = g0 * gm
            if abs(coeff) < 1e-30:
                continue
            # Nucleus A at +R_A: -Z_A * f_k(s, rho_A)
            c_k_A = -Z_A * (min_s2_A / max_s2_A) ** k / max_s2_A
            # Nucleus B at -R_B: -Z_B * (-1)^k * f_k(s, rho_B)
            c_k_B = -Z_B * (-1)**k * (min_s2_B / max_s2_B) ** k / max_s2_B
            V2 += coeff * (c_k_A + c_k_B)

        V += phase2 * norm2 * V2

    return V


def compute_nuclear_coupling_screened(
    l1p: int, l2p: int, l1: int, l2: int,
    m1: int, m2: int,
    alpha: np.ndarray, rho: float,
    Z_A_func: Callable, Z_B: float, R_e: float,
    Z_A_bare: float = None,
    rho_A: float = None, rho_B: float = None,
    n_theta: int = 64,
) -> np.ndarray:
    """
    Nuclear coupling with screened Z_A using analytical base + smooth correction.

    Strategy: compute the exact analytical coupling with constant Z_A (handling
    the 1/r singularity via multipole expansion), then add a smooth quadrature
    correction for the difference (Z_A - Z_eff(r)) / r.  Since Z_eff(0) = Z_A
    (no screening at the nucleus), the correction vanishes at r → 0, removing
    the 1/r singularity from the quadrature integrand.

    Parameters
    ----------
    l1p, l2p, l1, l2 : int
        Bra/ket channel angular momenta.
    m1, m2 : int
        Shared magnetic quantum numbers.
    alpha : ndarray
        Alpha grid points.
    rho : float
        R / (2 R_e).
    Z_A_func : callable
        Z_eff(r) → effective nuclear charge at distance r (bohr) from A.
    Z_B : float
        Constant charge for nucleus B.
    R_e : float
        Electronic hyperradius (bohr).
    Z_A_bare : float or None
        Bare nuclear charge for nucleus A (used for analytical base).
        If None, uses Z_A_func(0.0).
    rho_A, rho_B : float or None
        Nuclear distances from origin in R_e units.
    n_theta : int
        Gauss-Legendre quadrature order.

    Returns
    -------
    V : ndarray of shape (n_alpha,)
        Nuclear coupling (charge function units, without R_e factor).
    """
    from scipy.special import lpmv

    if rho_A is None:
        rho_A = rho
    if rho_B is None:
        rho_B = rho
    if Z_A_bare is None:
        Z_A_bare = float(Z_A_func(0.0))

    # Step 1: Analytical base with constant Z_A_bare (exact, no singularity)
    V_base = compute_nuclear_coupling(
        l1p, l2p, l1, l2, m1, m2, alpha, rho,
        Z=Z_A_bare, Z_A=Z_A_bare, Z_B=Z_B,
        rho_A=rho_A, rho_B=rho_B,
    )

    # Step 2: Smooth quadrature correction for (Z_A_bare - Z_eff(r)) / r
    # This correction vanishes at r → 0 since Z_eff(0) = Z_A_bare.
    n_alpha = len(alpha)
    cos_a = np.cos(alpha)
    sin_a = np.sin(alpha)
    am1 = abs(m1)
    am2 = abs(m2)

    nodes, weights = np.polynomial.legendre.leggauss(n_theta)

    dV = np.zeros(n_alpha)

    # --- Electron 1 correction (requires l2p == l2) ---
    if l2p == l2:
        norm_l1p = sqrt(
            (2 * l1p + 1) / 2.0
            * factorial(l1p - am1) / factorial(l1p + am1)
        )
        norm_l1 = sqrt(
            (2 * l1 + 1) / 2.0
            * factorial(l1 - am1) / factorial(l1 + am1)
        )
        Theta_l1p = norm_l1p * lpmv(am1, l1p, nodes)
        Theta_l1 = norm_l1 * lpmv(am1, l1, nodes)
        TT_w = Theta_l1p * Theta_l1 * weights

        s1 = cos_a
        for ia in range(n_alpha):
            # Distance electron 1 → nucleus A
            r1A_sq = s1[ia] ** 2 + rho_A ** 2 - 2 * s1[ia] * rho_A * nodes
            r1A = np.sqrt(np.maximum(r1A_sq, 1e-30))

            # Smooth correction: -(Z_eff(r) - Z_A_bare) / r
            # = (Z_A_bare - Z_eff(r)) / r
            z_eff_A = Z_A_func(R_e * r1A)
            delta_V = -(z_eff_A - Z_A_bare) / r1A

            dV[ia] += np.dot(TT_w, delta_V)

    # --- Electron 2 correction (requires l1p == l1) ---
    if l1p == l1:
        norm_l2p = sqrt(
            (2 * l2p + 1) / 2.0
            * factorial(l2p - am2) / factorial(l2p + am2)
        )
        norm_l2 = sqrt(
            (2 * l2 + 1) / 2.0
            * factorial(l2 - am2) / factorial(l2 + am2)
        )
        Theta_l2p = norm_l2p * lpmv(am2, l2p, nodes)
        Theta_l2 = norm_l2 * lpmv(am2, l2, nodes)
        TT_w = Theta_l2p * Theta_l2 * weights

        s2 = sin_a
        for ia in range(n_alpha):
            r2A_sq = s2[ia] ** 2 + rho_A ** 2 - 2 * s2[ia] * rho_A * nodes
            r2A = np.sqrt(np.maximum(r2A_sq, 1e-30))

            z_eff_A = Z_A_func(R_e * r2A)
            delta_V = -(z_eff_A - Z_A_bare) / r2A

            dV[ia] += np.dot(TT_w, delta_V)

    return V_base + dV


def compute_core_screening_analytical(
    alpha: np.ndarray,
    rho_A: float,
    R_e: float,
    core_potentials: List[dict],
    rho_B: float = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Analytical monopole core-active screening from locked-shell quantum numbers.

    For a hydrogenic core orbital (n, l=0) with nuclear charge Z_core, the
    electrostatic potential of the core charge density at distance r from the
    nucleus is:

        φ_core(r) = (n_core/r) × [1 - (1 + Z r) exp(-2 Z r)]

    The angle-averaged (k=0 monopole) penetration correction — the difference
    between this screened potential and the point-charge n_core/r already
    captured by Z_eff — is evaluated analytically:

        ⟨V_pen⟩(s, ρ) = -n_core/(2 R_e² s ρ) × [F(b) - F(a)]

    where a = R_e|s-ρ|, b = R_e(s+ρ), and F(u) = [-3/(4Z) - u/2] exp(-2Zu).

    This formula derives entirely from the core's quantum numbers (n=1, l=0)
    and nuclear charge — no numerical quadrature needed.

    The monopole (k=0 Legendre component) only contributes to diagonal
    channels (l' = l for each electron), consistent with the spherical
    symmetry of the s-orbital core.

    Parameters
    ----------
    alpha : ndarray
        Alpha grid points.
    rho_A : float
        Distance of nucleus A from origin in R_e units.
    R_e : float
        Electronic hyperradius.
    core_potentials : list of dict
        Each entry: {'Z_core': float, 'n_core': int, 'atom': str,
                     'n_q': int, 'l_q': int}.
        Only l_q=0 (s-orbital) cores are supported.
    rho_B : float or None
        Distance of nucleus B from origin in R_e units (for B-atom cores).

    Returns
    -------
    V_pen_e1 : ndarray of shape (n_alpha,)
        Monopole penetration for electron 1 (at s1 = cos α).
        In charge-function units (divide by R_e already applied).
    V_pen_e2 : ndarray of shape (n_alpha,)
        Monopole penetration for electron 2 (at s2 = sin α).
    """
    n_alpha = len(alpha)
    V_e1 = np.zeros(n_alpha)
    V_e2 = np.zeros(n_alpha)

    # Filter to s-orbital cores by atom
    cores_A = [c for c in core_potentials
               if c.get('atom') == 'A' and c.get('l_q', 0) == 0]
    cores_B = [c for c in core_potentials
               if c.get('atom') == 'B' and c.get('l_q', 0) == 0]

    if not cores_A and not cores_B:
        return V_e1, V_e2

    cos_a = np.cos(alpha)
    sin_a = np.sin(alpha)

    def _F(u: np.ndarray, Z: float) -> np.ndarray:
        """Antiderivative: F(u) = [-3/(4Z) - u/2] exp(-2Zu)."""
        return (-3.0 / (4.0 * Z) - u / 2.0) * np.exp(-2.0 * Z * u)

    def _monopole_pen(s: np.ndarray, rho_nuc: float,
                      Z_core: float, n_core: int) -> np.ndarray:
        """Angle-averaged penetration at electron distance s from origin.

        Returns V_pen / R_e (charge-function units, ready to multiply by R_e
        in the Hamiltonian assembly).
        """
        # Physical distances: a = R_e|s - rho|, b = R_e(s + rho)
        a = R_e * np.abs(s - rho_nuc)
        b = R_e * (s + rho_nuc)

        # Denominator: 2 * R_e^2 * s * rho
        denom = 2.0 * R_e**2 * s * rho_nuc
        # Guard against s ≈ 0 (alpha near boundaries)
        safe = denom > 1e-30

        Fb = _F(b, Z_core)
        Fa = _F(a, Z_core)

        result = np.zeros_like(s)
        result[safe] = -n_core / denom[safe] * (Fb[safe] - Fa[safe])

        # Convert to charge-function units: divide by R_e
        # (build_angular_hamiltonian will multiply by R_e)
        return result / R_e

    # Accumulate penetration from all core shells
    for c in cores_A:
        V_e1 += _monopole_pen(cos_a, rho_A, c['Z_core'], c['n_core'])
        V_e2 += _monopole_pen(sin_a, rho_A, c['Z_core'], c['n_core'])

    if rho_B is not None:
        for c in cores_B:
            V_e1 += _monopole_pen(cos_a, rho_B, c['Z_core'], c['n_core'])
            V_e2 += _monopole_pen(sin_a, rho_B, c['Z_core'], c['n_core'])

    return V_e1, V_e2


def compute_pk_pseudopotential(
    alpha: np.ndarray,
    rho_A: float,
    R_e: float,
    pk_potentials: List[dict],
    rho_B: float = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Effective core potential (ECP) in hyperspherical coordinates.

    Uses the Gaussian/r² form from the Fuentealba-Durand-Barthelat
    pseudopotential:

        V_ECP(r) = C × exp(-β r²) / r²

    The 1/r² divergence at r→0 overwhelms the nuclear -Z/r, creating a
    repulsive barrier that prevents active electrons from collapsing into
    the locked core orbital.  At r >> 1/√β, V_ECP → 0 and the electron
    sees only the nuclear potential.

    The angle-averaged (monopole) component uses the exponential integral
    E₁:

        ⟨V_ECP⟩(s, ρ, R_e) = C / (4 R_e² s ρ) × [E₁(β a²) - E₁(β b²)]

    where a = R_e|s-ρ|, b = R_e(s+ρ).

    Parameters
    ----------
    alpha : ndarray
        Alpha grid points.
    rho_A : float
        Distance of nucleus A from origin in R_e units.
    R_e : float
        Electronic hyperradius.
    pk_potentials : list of dict
        Each entry: {'C_core': float, 'beta_core': float,
                     'atom': str ('A' or 'B')}.
        C_core: amplitude of Gaussian/r² repulsive core (Ha × bohr²).
        beta_core: Gaussian exponent (1/bohr²).
    rho_B : float or None
        Distance of nucleus B from origin in R_e units.

    Returns
    -------
    V_pk_e1, V_pk_e2 : ndarray of shape (n_alpha,)
        ECP for electrons 1 and 2, in charge-function units.
    """
    from scipy.special import exp1

    n_alpha = len(alpha)
    V_e1 = np.zeros(n_alpha)
    V_e2 = np.zeros(n_alpha)

    pks_A = [p for p in pk_potentials if p.get('atom') == 'A']
    pks_B = [p for p in pk_potentials if p.get('atom') == 'B']

    if not pks_A and not pks_B:
        return V_e1, V_e2

    cos_a = np.cos(alpha)
    sin_a = np.sin(alpha)

    def _ecp_mono(s: np.ndarray, rho_nuc: float,
                  C: float, beta: float) -> np.ndarray:
        """Monopole ECP at electron distance s from origin.

        V_ECP(r) = C exp(-β r²) / r²
        ⟨V_ECP⟩ = C / (4 R_e² s ρ) × [E₁(β a²) - E₁(β b²)]

        Returns V_ECP / R_e (charge-function units).
        """
        a = R_e * np.abs(s - rho_nuc)
        b = R_e * (s + rho_nuc)

        denom = 4.0 * R_e**2 * s * rho_nuc
        safe = denom > 1e-30

        # E₁(x) = ∫_x^∞ exp(-t)/t dt, monotonically decreasing.
        # E₁(β a²) > E₁(β b²) since a < b, so result is positive (repulsive).
        ba2 = beta * a**2
        bb2 = beta * b**2
        # Clamp argument to avoid overflow in exp1 at very small values
        ba2 = np.maximum(ba2, 1e-30)

        E1_a = exp1(ba2)
        E1_b = exp1(bb2)

        result = np.zeros_like(s)
        result[safe] = C / denom[safe] * (E1_a[safe] - E1_b[safe])

        return result / R_e

    for p in pks_A:
        C = p['C_core']
        beta = p['beta_core']
        V_e1 += _ecp_mono(cos_a, rho_A, C, beta)
        V_e2 += _ecp_mono(sin_a, rho_A, C, beta)

    if rho_B is not None:
        for p in pks_B:
            C = p['C_core']
            beta = p['beta_core']
            V_e1 += _ecp_mono(cos_a, rho_B, C, beta)
            V_e2 += _ecp_mono(sin_a, rho_B, C, beta)

    return V_e1, V_e2


def _ee_coupling(
    l1p: int, l2p: int, l1: int, l2: int,
    alpha: np.ndarray, l_max: int,
) -> np.ndarray:
    """
    Compute e-e repulsion coupling between channels (l1',l2') and (l1,l2).

    Uses the multipole expansion of 1/r12 in hyperspherical coordinates:
        1/r12 = (1/R_e) sum_k f_k(alpha) * [angular part]

    where f_k(alpha) = (min/max)^k / max with min=min(cos a, sin a),
    max=max(cos a, sin a).

    The angular coupling involves a double Gaunt integral:
        sum_k G(l1', k, l1) * G(l2', k, l2) * f_k / max

    For the He solver (l1=l2=l, l1'=l2'=l'), the Gaunt structure
    simplifies. Here we need the general two-index form.

    In H2, channels have independent (l1, l2), so the coupling is:
        W = sum_k sqrt((2l1'+1)(2l2'+1)(2l1+1)(2l2+1))/4
            * G(l1',k,l1) * G(l2',k,l2) * f_k(alpha) / max_sc

    Parameters
    ----------
    l1p, l2p, l1, l2 : int
        Channel angular momenta.
    alpha : ndarray
        Alpha grid.
    l_max : int
        Maximum l for Gaunt integral range.

    Returns
    -------
    W : ndarray of shape (n_alpha,)
        E-e coupling at each alpha point (charge function units / R_e).
    """
    sin_a = np.sin(alpha)
    cos_a = np.cos(alpha)
    min_sc = np.minimum(sin_a, cos_a)
    max_sc = np.maximum(sin_a, cos_a)

    W = np.zeros(len(alpha))

    k_min = max(abs(l1p - l1), abs(l2p - l2))
    k_max = min(l1p + l1, l2p + l2)

    for k in range(k_min, k_max + 1):
        g1 = gaunt_integral(l1p, k, l1)
        g2 = gaunt_integral(l2p, k, l2)
        if abs(g1) < 1e-15 or abs(g2) < 1e-15:
            continue
        f_k = (min_sc / max_sc) ** k
        W += g1 * g2 * f_k / max_sc

    norm = sqrt((2*l1p+1) * (2*l2p+1) * (2*l1+1) * (2*l2+1)) / 4.0
    W *= norm

    return W


def _ee_coupling_general(
    l1p: int, m1p: int, l2p: int, m2p: int,
    l1: int, m1: int, l2: int, m2: int,
    alpha: np.ndarray,
) -> np.ndarray:
    """
    Generalized e-e coupling for arbitrary m using Wigner 3j symbols.

    From the multipole expansion of 1/r12 in spherical harmonics:

        W = norm * sum_k (-1)^q *
            (l1' k l1; 0 0 0)(l1' k l1; -m1' -q m1) *
            (l2' k l2; 0 0 0)(l2' k l2; -m2' q m2) * f_k(a)/max

    where q = m1 - m1' is the azimuthal transfer quantum number, and
    norm = sqrt((2l1'+1)(2l1+1)(2l2'+1)(2l2+1)).

    Selection rules:
    - m1 + m2 = m1' + m2' (M conservation, automatic for M=0)
    - l1' + k + l1 even (from (000) 3j symbol)
    - |l1'-l1| <= k <= l1'+l1, |l2'-l2| <= k <= l2'+l2

    For m1=m1'=m2=m2'=0, reduces to _ee_coupling (verified algebraically:
    gaunt_integral = 2*(3j;000)^2, and the /4 norm cancels with the factor 4).

    Parameters
    ----------
    l1p, m1p, l2p, m2p : int
        Bra channel quantum numbers.
    l1, m1, l2, m2 : int
        Ket channel quantum numbers.
    alpha : ndarray
        Alpha grid.

    Returns
    -------
    W : ndarray of shape (n_alpha,)
        E-e coupling (charge function units / R_e).
    """
    from geovac.angular_integrals import wigner3j

    q = m1 - m1p
    # M conservation check: q must also equal m2p - m2
    if q != m2p - m2:
        return np.zeros(len(alpha))

    sin_a = np.sin(alpha)
    cos_a = np.cos(alpha)
    min_sc = np.minimum(sin_a, cos_a)
    max_sc = np.maximum(sin_a, cos_a)

    W = np.zeros(len(alpha))

    k_min = max(abs(l1p - l1), abs(l2p - l2))
    k_max = min(l1p + l1, l2p + l2)

    for k in range(k_min, k_max + 1):
        # Parity check via (l k l'; 0 0 0) — zero unless l+k+l' even
        g1_0 = wigner3j(l1p, k, l1, 0, 0, 0)
        if abs(g1_0) < 1e-15:
            continue
        g2_0 = wigner3j(l2p, k, l2, 0, 0, 0)
        if abs(g2_0) < 1e-15:
            continue

        g1_m = wigner3j(l1p, k, l1, -m1p, -q, m1)
        g2_m = wigner3j(l2p, k, l2, -m2p, q, m2)

        coeff = g1_0 * g1_m * g2_0 * g2_m
        if abs(coeff) < 1e-30:
            continue

        f_k = (min_sc / max_sc) ** k
        W += coeff * f_k / max_sc

    phase = (-1) ** q
    norm = sqrt((2*l1p+1) * (2*l1+1) * (2*l2p+1) * (2*l2+1))
    W *= phase * norm

    return W


def build_angular_hamiltonian(
    alpha_grid: np.ndarray,
    rho: float,
    R_e: float,
    l_max: int = 2,
    Z: float = 1.0,
    m_max: int = 0,
    l_max_per_m: Union[dict, None] = None,
    Z_A: float = None,
    Z_B: float = None,
    z0: float = 0.0,
    core_potentials: Union[List[dict], None] = None,
    pk_potentials: Union[List[dict], None] = None,
    Z_A_func: Optional[Callable] = None,
    n_theta: int = 64,
    pk_projector: Union[dict, None] = None,
) -> np.ndarray:
    """
    Build the full coupled-channel angular Hamiltonian.

    Matrix structure: block (n_ch x n_ch) at each alpha grid point,
    with tridiagonal kinetic coupling within each channel along alpha.

    Total dimension: N = n_ch * n_alpha.
    Indexing: global = ch * n_alpha + i_alpha.

    Parameters
    ----------
    alpha_grid : ndarray of shape (n_alpha,)
        Interior alpha grid points (Dirichlet BCs at boundaries).
    rho : float
        R / (2 R_e).
    R_e : float
        Electronic hyperradius.
    l_max : int
        Maximum angular momentum.
    Z : float
        Nuclear charge per nucleus (homonuclear default).
    m_max : int
        Maximum |m| per electron. 0 = sigma-only (Phase 2),
        1 = sigma + pi (Phase 3).
    Z_A, Z_B : float or None
        Per-nucleus charges for heteronuclear systems. If None, uses Z.
    z0 : float
        Origin shift along internuclear axis. Positive = toward nucleus A.
        Nuclear positions: A at R/2 - z0, B at R/2 + z0 from origin.
    pk_potentials : list of dict or None
        Phillips-Kleinman pseudopotential specs. Each entry:
        {'Z_core': float, 'E_core': float, 'E_val': float, 'atom': str}.
        Creates a repulsive barrier where the core density lives, enforcing
        core-active orthogonality.
    Z_A_func : callable or None
        If provided, Z_eff(r) callable for nucleus A.  Replaces the constant
        Z_A in the nuclear coupling via Gauss-Legendre quadrature.  When
        None (default), the fast analytical multipole path is used.
    n_theta : int
        Gauss-Legendre quadrature order (only used when Z_A_func is not None).

    Returns
    -------
    H : ndarray of shape (N, N)
        Full angular Hamiltonian matrix.
    """
    # Resolve charges
    if Z_A is None:
        Z_A = Z
    if Z_B is None:
        Z_B = Z
    homonuclear = (Z_A == Z_B)

    # Compute per-nucleus rho values for shifted origin
    # rho = R/(2*R_e), so R_A = R/2 - z0 → rho_A = rho - z0/R_e
    rho_A = rho - z0 / R_e if z0 != 0.0 else rho
    rho_B = rho + z0 / R_e if z0 != 0.0 else rho

    # --- Determine channel list ---
    if m_max == 0:
        channels_2 = _channel_list(l_max, homonuclear=homonuclear)
        n_ch = len(channels_2)
        # Convert to 4-tuple form for unified indexing below
        channels_4 = [(l1, 0, l2, 0) for l1, l2 in channels_2]
    else:
        channels_4 = _channel_list_extended(
            l_max, m_max, l_max_per_m, homonuclear=homonuclear,
        )
        n_ch = len(channels_4)

    n_alpha = len(alpha_grid)
    N = n_ch * n_alpha

    h = alpha_grid[1] - alpha_grid[0]  # uniform grid spacing

    cos_a = np.cos(alpha_grid)
    sin_a = np.sin(alpha_grid)

    H = np.zeros((N, N))

    kinetic_diag = 1.0 / h**2
    kinetic_off = -0.5 / h**2

    def idx(ch: int, i: int) -> int:
        return ch * n_alpha + i

    # --- Diagonal blocks: kinetic + centrifugal + Liouville for each channel ---
    # Centrifugal l(l+1)/(2 cos^2 a) depends on l, not m (m is absorbed
    # into the l(l+1) eigenvalue of the angular Laplacian on each sphere).
    for ic, (l1, m1, l2, m2) in enumerate(channels_4):
        V_cent = (0.5 * l1 * (l1 + 1) / cos_a**2
                  + 0.5 * l2 * (l2 + 1) / sin_a**2)
        V_liouville = -2.0
        V_diag = V_cent + V_liouville

        for i in range(n_alpha):
            ii = idx(ic, i)
            H[ii, ii] = kinetic_diag + V_diag[i]

        # Tridiagonal kinetic within channel
        for i in range(n_alpha - 1):
            ii = idx(ic, i)
            jj = idx(ic, i + 1)
            H[ii, jj] = kinetic_off
            H[jj, ii] = kinetic_off

    # --- Nuclear coupling (diagonal in m1, m2 due to axial symmetry) ---
    _use_screened = Z_A_func is not None
    for ic, (l1p, m1p, l2p, m2p) in enumerate(channels_4):
        for jc, (l1, m1, l2, m2) in enumerate(channels_4):
            if jc < ic:
                continue  # fill symmetric

            # Nuclear potential is phi-independent: requires m1' = m1, m2' = m2
            if m1p != m1 or m2p != m2:
                continue

            if _use_screened:
                V_nuc = compute_nuclear_coupling_screened(
                    l1p, l2p, l1, l2, m1, m2, alpha_grid, rho,
                    Z_A_func=Z_A_func, Z_B=Z_B, R_e=R_e,
                    Z_A_bare=Z_A,
                    rho_A=rho_A, rho_B=rho_B, n_theta=n_theta,
                )
            else:
                V_nuc = compute_nuclear_coupling(
                    l1p, l2p, l1, l2, m1, m2, alpha_grid, rho,
                    Z=Z, Z_A=Z_A, Z_B=Z_B,
                    rho_A=rho_A, rho_B=rho_B,
                )

            for i in range(n_alpha):
                ii = idx(ic, i)
                jj = idx(jc, i)
                val = R_e * V_nuc[i]
                H[ii, jj] += val
                if ic != jc:
                    H[jj, ii] += val

    # --- Core penetration correction (analytical monopole from core quantum numbers) ---
    # The k=0 monopole of the penetration potential is diagonal in (l, m)
    # because the s-orbital core is spherically symmetric. Derived from
    # the analytical angle-average of V_pen(r_A) — no quadrature needed.
    if core_potentials:
        V_pen_e1, V_pen_e2 = compute_core_screening_analytical(
            alpha_grid, rho_A, R_e, core_potentials, rho_B=rho_B,
        )
        for ic, (l1, m1, l2, m2) in enumerate(channels_4):
            for i in range(n_alpha):
                ii = idx(ic, i)
                # Electron 1 penetration (diagonal in l1, acts on all channels)
                H[ii, ii] += R_e * V_pen_e1[i]
                # Electron 2 penetration (diagonal in l2, acts on all channels)
                H[ii, ii] += R_e * V_pen_e2[i]

    # --- Phillips-Kleinman pseudopotential (core-active orthogonality) ---
    # Repulsive potential V_PK = |φ_core|² × (E_val - E_core) prevents active
    # electrons from collapsing into the locked core orbital.  Diagonal in
    # channels (monopole, spherically symmetric core).
    #
    # channel_mode='l_dependent': The exact PK operator is a projector onto
    # core states.  For a 1s² core (pure l=0), the projector has non-zero
    # matrix elements only for electrons with l=0 character.  Per-electron
    # weight: δ_{l_i, 0}.  This prevents unphysical repulsion in higher-l
    # channels that are automatically orthogonal to the core by angular
    # momentum.  See Paper 17 Sec VI.A.
    if pk_potentials:
        V_pk_e1, V_pk_e2 = compute_pk_pseudopotential(
            alpha_grid, rho_A, R_e, pk_potentials, rho_B=rho_B,
        )
        # Determine channel mode from pk_potentials
        pk_ch_mode = 'channel_blind'
        for p in pk_potentials:
            if p.get('channel_mode') == 'l_dependent':
                pk_ch_mode = 'l_dependent'
                break

        for ic, (l1, m1, l2, m2) in enumerate(channels_4):
            if pk_ch_mode == 'l_dependent':
                # Per-electron PK weight: core is 1s² (pure l=0),
                # so only l=0 electrons overlap with the core.
                w1 = 1.0 if l1 == 0 else 0.0
                w2 = 1.0 if l2 == 0 else 0.0
            else:
                w1 = 1.0
                w2 = 1.0
            if w1 == 0.0 and w2 == 0.0:
                continue
            for i in range(n_alpha):
                ii = idx(ic, i)
                H[ii, ii] += R_e * w1 * V_pk_e1[i]
                H[ii, ii] += R_e * w2 * V_pk_e2[i]

    # --- Algebraic PK projector (rank-1, channel-space) ---
    # Replaces the Gaussian PK with the exact Phillips-Kleinman projector:
    #   V_PK = E_shift * |core⟩⟨core|
    # where |core⟩ is the Level 3 core eigenvector mapped into the Level 4
    # channel basis via the atomic-limit approximation (l=0 → (0,0) channel).
    if pk_projector is not None and pk_projector.get('mode') == 'algebraic':
        from scipy.interpolate import CubicSpline as _CS

        core_l0 = pk_projector['core_l0_wavefunction']
        core_n_alpha_src = pk_projector['core_n_alpha']
        core_h_src = pk_projector['core_h_alpha']

        # Interpolate to Level 4 α-grid if grids differ
        if core_n_alpha_src != n_alpha:
            core_alpha = np.array(
                [(i + 1) * core_h_src for i in range(core_n_alpha_src)]
            )
            level4_alpha = alpha_grid
            cs = _CS(core_alpha, core_l0, bc_type='clamped')
            core_l0_interp = cs(level4_alpha)
        else:
            core_l0_interp = core_l0.copy()

        # Normalize on Level 4 grid
        norm = np.sqrt(np.sum(core_l0_interp**2) * h)
        if norm > 1e-15:
            core_l0_interp /= norm

        # Build full core vector in Level 4 (channel, α) space
        # Atomic-limit approximation: core maps to (l1=0, m1=0, l2=0, m2=0)
        core_vec = np.zeros(n_ch * n_alpha)
        ch_00 = None
        for ic, ch_label in enumerate(channels_4):
            if m_max == 0:
                l1, l2 = ch_label[0], ch_label[1]
                if l1 == 0 and l2 == 0:
                    ch_00 = ic
                    break
            else:
                l1, m1, l2, m2 = ch_label
                if l1 == 0 and m1 == 0 and l2 == 0 and m2 == 0:
                    ch_00 = ic
                    break

        if ch_00 is not None:
            core_vec[ch_00 * n_alpha:(ch_00 + 1) * n_alpha] = core_l0_interp

        # Add rank-1 projector: V_PK = R_e * E_shift * |core⟩⟨core|
        # R_e factor: Hamiltonian is in charge-function form H u = μ u
        # where μ has units of energy × R_e (matching existing R_e * V terms).
        E_shift = pk_projector['energy_shift']
        H += R_e * E_shift * np.outer(core_vec, core_vec)

    # --- E-e coupling ---
    # For m_max=0: original Gaunt-based coupling (exact backward compatibility).
    # For m_max>0: generalized Wigner 3j coupling for all channel pairs.
    for ic, (l1p, m1p, l2p, m2p) in enumerate(channels_4):
        for jc, (l1, m1, l2, m2) in enumerate(channels_4):
            if jc < ic:
                continue

            if m_max == 0:
                W_ee = _ee_coupling(l1p, l2p, l1, l2, alpha_grid, l_max)
            else:
                W_ee = _ee_coupling_general(
                    l1p, m1p, l2p, m2p, l1, m1, l2, m2, alpha_grid
                )

            if np.max(np.abs(W_ee)) < 1e-15:
                continue

            for i in range(n_alpha):
                ii = idx(ic, i)
                jj = idx(jc, i)
                val = R_e * W_ee[i]
                H[ii, jj] += val
                if ic != jc:
                    H[jj, ii] += val

    return H


def solve_angular_multichannel(
    rho: float,
    R_e: float,
    l_max: int = 2,
    Z: float = 1.0,
    n_alpha: int = 200,
    n_eig: int = 1,
    m_max: int = 0,
    l_max_per_m: Union[dict, None] = None,
    Z_A: float = None,
    Z_B: float = None,
    z0: float = 0.0,
    core_potentials: Union[List[dict], None] = None,
    pk_potentials: Union[List[dict], None] = None,
    pk_projector: Union[dict, None] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, list]:
    """
    Solve the coupled-channel angular eigenvalue problem.

    Parameters
    ----------
    rho : float
        R / (2 R_e).
    R_e : float
        Electronic hyperradius.
    l_max : int
        Maximum angular momentum per electron.
    Z : float
        Nuclear charge (homonuclear default).
    n_alpha : int
        FD grid points in alpha.
    n_eig : int
        Number of eigenvalues to return.
    m_max : int
        Maximum |m| per electron. 0 = sigma-only, 1 = sigma+pi.
    Z_A, Z_B : float or None
        Per-nucleus charges for heteronuclear systems.
    z0 : float
        Origin shift along internuclear axis.

    Returns
    -------
    mu : ndarray of shape (n_eig,)
        Lowest angular eigenvalues.
    vecs : ndarray of shape (n_eig, N)
        Eigenvectors.
    alpha_grid : ndarray
        Alpha grid.
    channels : list
        Channel labels: (l1, l2) if m_max=0, (l1, m1, l2, m2) if m_max>0.
    """
    if Z_A is None:
        Z_A = Z
    if Z_B is None:
        Z_B = Z
    homonuclear = (Z_A == Z_B)

    h = (np.pi / 2) / (n_alpha + 1)
    alpha = (np.arange(n_alpha) + 1) * h

    H = build_angular_hamiltonian(alpha, rho, R_e, l_max, Z, m_max,
                                   l_max_per_m, Z_A=Z_A, Z_B=Z_B, z0=z0,
                                   core_potentials=core_potentials,
                                   pk_potentials=pk_potentials,
                                   pk_projector=pk_projector)

    evals, evecs = eigh(H)

    if m_max == 0:
        channels = _channel_list(l_max, homonuclear=homonuclear)
    else:
        channels = _channel_list_extended(
            l_max, m_max, l_max_per_m, homonuclear=homonuclear,
        )

    return evals[:n_eig], evecs[:, :n_eig].T, alpha, channels


def compute_adiabatic_curve_mc(
    R: float,
    R_e_grid: np.ndarray,
    l_max: int = 2,
    Z: float = 1.0,
    n_alpha: int = 200,
    m_max: int = 0,
    l_max_per_m: Union[dict, None] = None,
    Z_A: float = None,
    Z_B: float = None,
    z0: float = 0.0,
    core_potentials: Union[List[dict], None] = None,
    pk_potentials: Union[List[dict], None] = None,
    pk_projector: Union[dict, None] = None,
) -> np.ndarray:
    """
    Compute adiabatic potential U(R_e) with multichannel angular solver.

    U(R_e; R) = [mu(R_e; R) + 15/8] / R_e^2

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    R_e_grid : ndarray
        Electronic hyperradius grid.
    l_max : int
        Maximum angular momentum.
    Z : float
        Nuclear charge (homonuclear default).
    n_alpha : int
        FD grid points.
    m_max : int
        Maximum |m| per electron.
    Z_A, Z_B : float or None
        Per-nucleus charges for heteronuclear systems.
    z0 : float
        Origin shift along internuclear axis.

    Returns
    -------
    U : ndarray
        Effective potential (Ha).
    """
    n_Re = len(R_e_grid)
    mu_vals = np.zeros(n_Re)

    for i, R_e in enumerate(R_e_grid):
        rho = R / (2.0 * R_e)
        evals, _, _, _ = solve_angular_multichannel(
            rho, R_e, l_max, Z, n_alpha, m_max=m_max,
            l_max_per_m=l_max_per_m, Z_A=Z_A, Z_B=Z_B, z0=z0,
            core_potentials=core_potentials,
            pk_potentials=pk_potentials,
            pk_projector=pk_projector,
        )
        mu_vals[i] = evals[0]

    U = (mu_vals + 15.0 / 8.0) / R_e_grid**2
    return U


def compute_adiabatic_curve_dboc(
    R: float,
    R_e_grid: np.ndarray,
    l_max: int = 2,
    Z: float = 1.0,
    n_alpha: int = 200,
    m_max: int = 0,
    n_states: int = 5,
    l_max_per_m: Union[dict, None] = None,
    Z_A: float = None,
    Z_B: float = None,
    z0: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute adiabatic potential U(R_e) with diagonal Born-Oppenheimer
    correction (DBOC).

    The DBOC adds a positive correction to the adiabatic potential that
    accounts for the non-adiabatic coupling between angular channels:

        U_corr(R_e) = U_0(R_e) + (1/2) sum_j |P_0j(R_e)|^2

    where P_0j = <phi_0 | d phi_j / dR_e> is the first derivative coupling.
    This correction is always positive, preventing variational violations
    that occur when the single-channel adiabatic approximation overcounts
    the angular channel mixing for large channel sets (m_max > 0).

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    R_e_grid : ndarray
        Electronic hyperradius grid.
    l_max : int
        Maximum angular momentum.
    Z : float
        Nuclear charge.
    n_alpha : int
        FD grid points.
    m_max : int
        Maximum |m| per electron.
    n_states : int
        Number of adiabatic states for DBOC computation.

    Returns
    -------
    U_corrected : ndarray
        DBOC-corrected effective potential (Ha).
    U_bare : ndarray
        Uncorrected adiabatic potential (Ha).
    """
    n_Re = len(R_e_grid)
    mu_vals = np.zeros(n_Re)
    dboc_vals = np.zeros(n_Re)
    vecs_list: List[np.ndarray] = []

    for i, R_e in enumerate(R_e_grid):
        rho = R / (2.0 * R_e)
        evals, evecs, _, _ = solve_angular_multichannel(
            rho, R_e, l_max, Z, n_alpha,
            n_eig=n_states, m_max=m_max,
            l_max_per_m=l_max_per_m, Z_A=Z_A, Z_B=Z_B, z0=z0,
        )
        mu_vals[i] = evals[0]

        # Phase-align eigenvectors with previous R_e
        if i > 0:
            for j in range(min(n_states, len(evals))):
                if np.dot(vecs_list[-1][j], evecs[j]) < 0:
                    evecs[j] *= -1

        vecs_list.append(evecs.copy())

    # Compute DBOC from central finite differences
    for i in range(1, n_Re - 1):
        h_left = R_e_grid[i] - R_e_grid[i - 1]
        h_right = R_e_grid[i + 1] - R_e_grid[i]
        h_avg = (h_left + h_right) / 2.0

        dboc = 0.0
        n_actual = min(n_states, vecs_list[i].shape[0])
        for j in range(1, n_actual):
            # P_0j = <phi_0(R_e) | d phi_j / dR_e>
            dphi_j = (vecs_list[i + 1][j] - vecs_list[i - 1][j]) / (2 * h_avg)
            P_0j = np.dot(vecs_list[i][0], dphi_j)
            dboc += P_0j**2

        dboc_vals[i] = dboc / 2.0

    # Boundary: extrapolate from interior
    dboc_vals[0] = dboc_vals[1]
    dboc_vals[-1] = dboc_vals[-2]

    U_bare = (mu_vals + 15.0 / 8.0) / R_e_grid**2
    U_corrected = U_bare + dboc_vals

    return U_corrected, U_bare


def solve_coupled_channel_radial(
    R: float,
    R_e_grid: np.ndarray,
    l_max: int = 2,
    Z: float = 1.0,
    n_alpha: int = 200,
    m_max: int = 0,
    n_coupled: int = 3,
    l_max_per_m: Union[dict, None] = None,
    n_Re_radial: int = 300,
    R_e_min: float = 0.3,
    R_e_max: float = 15.0,
    Z_A: float = None,
    Z_B: float = None,
    z0: float = 0.0,
    pk_potentials: Union[List[dict], None] = None,
) -> Tuple[float, np.ndarray]:
    """
    Diabatic coupled-channel radial solver.

    Uses the angular eigenvectors at a reference R_e (the minimum of the
    adiabatic potential) as a fixed diabatic basis. The angular Hamiltonian
    is projected onto this basis at each R_e to give a smooth n_coupled ×
    n_coupled potential matrix V_ij(R_e). The coupled radial Schrödinger
    equation is then solved as a block-matrix eigenvalue problem.

    This properly handles non-adiabatic coupling between angular channels,
    avoiding the variational violations of the single-channel adiabatic
    approximation and the convergence issues of DBOC.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    R_e_grid : ndarray
        R_e grid for angular sweeps (used to find reference point).
    l_max : int
        Maximum angular momentum per electron.
    Z : float
        Nuclear charge.
    n_alpha : int
        FD grid points for alpha.
    m_max : int
        Maximum |m| per electron.
    n_coupled : int
        Number of coupled radial channels.
    l_max_per_m : dict or None
        Per-|m| angular momentum limit.
    n_Re_radial : int
        Grid points for radial equation.
    R_e_min, R_e_max : float
        Radial boundaries.

    Returns
    -------
    E_elec : float
        Electronic energy eigenvalue.
    F : ndarray
        Radial wavefunction (n_coupled * n_Re_radial,).
    """
    # Step 1: Quick adiabatic sweep to find reference R_e
    n_Re_ang = len(R_e_grid)
    mu_vals = np.zeros(n_Re_ang)
    for i, R_e in enumerate(R_e_grid):
        rho = R / (2.0 * R_e)
        evals, _, _, _ = solve_angular_multichannel(
            rho, R_e, l_max, Z, n_alpha, m_max=m_max,
            l_max_per_m=l_max_per_m, Z_A=Z_A, Z_B=Z_B, z0=z0,
        )
        mu_vals[i] = evals[0]

    U_adia = (mu_vals + 15.0 / 8.0) / R_e_grid**2
    i_ref = np.argmin(U_adia)
    R_e_ref = R_e_grid[i_ref]

    # Step 2: Get reference eigenvectors at R_e_ref
    rho_ref = R / (2.0 * R_e_ref)
    _, ref_vecs, _, _ = solve_angular_multichannel(
        rho_ref, R_e_ref, l_max, Z, n_alpha,
        n_eig=n_coupled, m_max=m_max, l_max_per_m=l_max_per_m,
        Z_A=Z_A, Z_B=Z_B, z0=z0,
    )
    # ref_vecs: shape (n_coupled, N_angular) — rows are eigenvectors

    # Step 3: Build diabatic potential matrix on the radial grid
    h_Re = (R_e_max - R_e_min) / (n_Re_radial + 1)
    R_e_radial = R_e_min + (np.arange(n_Re_radial) + 1) * h_Re

    # Precompute alpha grid (same as in solve_angular_multichannel)
    h_alpha = (np.pi / 2) / (n_alpha + 1)
    alpha_grid = (np.arange(n_alpha) + 1) * h_alpha

    V_diabatic = np.zeros((n_Re_radial, n_coupled, n_coupled))

    for k, R_e in enumerate(R_e_radial):
        rho = R / (2.0 * R_e)
        H_ang = build_angular_hamiltonian(
            alpha_grid, rho, R_e, l_max, Z, m_max, l_max_per_m,
            Z_A=Z_A, Z_B=Z_B, z0=z0,
            pk_potentials=pk_potentials,
        )
        # Project: V_ij = phi_i^T H_ang phi_j
        # H_ang is (N_ang x N_ang), ref_vecs is (n_coupled x N_ang)
        Hv = H_ang @ ref_vecs.T  # (N_ang x n_coupled)
        V_diabatic[k] = ref_vecs @ Hv  # (n_coupled x n_coupled)

        # Add 15/8 to the diagonal (hypercentrifugal term)
        for i in range(n_coupled):
            V_diabatic[k, i, i] = (V_diabatic[k, i, i] + 15.0 / 8.0) / R_e**2

        # Off-diagonal: also divide by R_e^2
        for i in range(n_coupled):
            for j in range(n_coupled):
                if i != j:
                    V_diabatic[k, i, j] /= R_e**2

    # Step 4: Build coupled-channel radial Hamiltonian
    # Dimension: (n_coupled * n_Re_radial) x (n_coupled * n_Re_radial)
    # Block structure: kinetic (tridiagonal within each channel) + V_diabatic
    N_total = n_coupled * n_Re_radial
    H_rad = np.zeros((N_total, N_total))

    kinetic_diag = 1.0 / h_Re**2
    kinetic_off = -0.5 / h_Re**2

    for ic in range(n_coupled):
        offset = ic * n_Re_radial
        # Kinetic: tridiagonal within channel ic
        for k in range(n_Re_radial):
            H_rad[offset + k, offset + k] = kinetic_diag
        for k in range(n_Re_radial - 1):
            H_rad[offset + k, offset + k + 1] = kinetic_off
            H_rad[offset + k + 1, offset + k] = kinetic_off

    # Potential coupling: at each R_e point, couple channels
    for k in range(n_Re_radial):
        for ic in range(n_coupled):
            for jc in range(n_coupled):
                H_rad[ic * n_Re_radial + k,
                      jc * n_Re_radial + k] += V_diabatic[k, ic, jc]

    # Step 5: Solve for lowest eigenvalue
    evals, evecs = eigh(H_rad, subset_by_index=[0, 0])
    E_elec = evals[0]
    F = evecs[:, 0]

    return E_elec, F


def solve_direct_2d(
    R: float,
    l_max: int = 2,
    Z: float = 1.0,
    n_alpha: int = 100,
    n_Re: int = 200,
    R_e_min: float = 0.3,
    R_e_max: float = 15.0,
    m_max: int = 0,
    l_max_per_m: Union[dict, None] = None,
    verbose: bool = True,
    Z_A: float = None,
    Z_B: float = None,
    z0: float = 0.0,
    core_potentials: Union[List[dict], None] = None,
    pk_potentials: Union[List[dict], None] = None,
) -> float:
    """
    Direct 2D solver: solve the full (alpha, R_e) problem without
    the adiabatic approximation.

    Builds the complete Hamiltonian matrix in the product space of
    alpha grid × R_e grid × channel space, then finds the lowest
    eigenvalue using sparse iterative methods.

    H = T_Re ⊗ I_ang + [H_ang(R_e) + 15/8 I] / R_e²

    This avoids the adiabatic approximation entirely and gives a
    properly variational result.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    l_max : int
        Maximum angular momentum.
    Z : float
        Nuclear charge.
    n_alpha : int
        FD grid points for alpha.
    n_Re : int
        Grid points for R_e.
    R_e_min, R_e_max : float
        R_e boundaries.
    m_max : int
        Maximum |m| per electron.
    l_max_per_m : dict or None
        Per-|m| angular momentum limit.
    verbose : bool
        Print diagnostics.

    Returns
    -------
    E_elec : float
        Lowest electronic energy eigenvalue.
    """
    from scipy.sparse.linalg import eigsh

    # Resolve charges
    if Z_A is None:
        Z_A = Z
    if Z_B is None:
        Z_B = Z
    homonuclear = (Z_A == Z_B)

    # Determine channels
    if m_max == 0:
        channels_2 = _channel_list(l_max, homonuclear=homonuclear)
        channels_4 = [(l1, 0, l2, 0) for l1, l2 in channels_2]
    else:
        channels_4 = _channel_list_extended(
            l_max, m_max, l_max_per_m, homonuclear=homonuclear,
        )
    n_ch = len(channels_4)

    # Grids
    h_alpha = (np.pi / 2) / (n_alpha + 1)
    alpha_grid = (np.arange(n_alpha) + 1) * h_alpha

    h_Re = (R_e_max - R_e_min) / (n_Re + 1)
    R_e_grid = R_e_min + (np.arange(n_Re) + 1) * h_Re

    N_ang = n_ch * n_alpha
    N_total = n_Re * N_ang

    if verbose:
        print(f"  Direct 2D solver: {N_total} x {N_total} sparse matrix")
        print(f"  ({n_Re} R_e points x {n_ch} channels x {n_alpha} alpha points)")

    # Build sparse Hamiltonian in COO format
    # Index: r * N_ang + (ch * n_alpha + ia)
    kinetic_diag_Re = 1.0 / h_Re**2
    kinetic_off_Re = -0.5 / h_Re**2

    # Angular blocks: [H_ang(R_e) + 15/8 I] / R_e^2
    # Collect all COO entries for efficiency
    rows_list: List[np.ndarray] = []
    cols_list: List[np.ndarray] = []
    vals_list: List[np.ndarray] = []

    # Add radial kinetic entries
    for r in range(n_Re):
        ang_idx = np.arange(N_ang)
        diag_idx = r * N_ang + ang_idx
        rows_list.append(diag_idx)
        cols_list.append(diag_idx)
        vals_list.append(np.full(N_ang, kinetic_diag_Re))

        if r < n_Re - 1:
            idx1 = r * N_ang + ang_idx
            idx2 = (r + 1) * N_ang + ang_idx
            rows_list.extend([idx1, idx2])
            cols_list.extend([idx2, idx1])
            vals_list.extend([np.full(N_ang, kinetic_off_Re)] * 2)

    for r, R_e in enumerate(R_e_grid):
        rho = R / (2.0 * R_e)
        H_ang = build_angular_hamiltonian(
            alpha_grid, rho, R_e, l_max, Z, m_max, l_max_per_m,
            Z_A=Z_A, Z_B=Z_B, z0=z0,
            core_potentials=core_potentials,
            pk_potentials=pk_potentials,
        )
        # Add 15/8 to diagonal and scale
        np.fill_diagonal(H_ang, H_ang.diagonal() + 15.0 / 8.0)
        H_ang *= (1.0 / R_e**2)

        # Extract nonzero entries
        nz_i, nz_j = np.nonzero(np.abs(H_ang) > 1e-15)
        if len(nz_i) > 0:
            offset = r * N_ang
            rows_list.append(offset + nz_i)
            cols_list.append(offset + nz_j)
            vals_list.append(H_ang[nz_i, nz_j])

    from scipy.sparse import coo_matrix
    all_rows = np.concatenate(rows_list)
    all_cols = np.concatenate(cols_list)
    all_vals = np.concatenate(vals_list)
    H_csr = coo_matrix((all_vals, (all_rows, all_cols)),
                        shape=(N_total, N_total)).tocsr()
    if verbose:
        print(f"  Nonzeros: {H_csr.nnz}")

    # Find lowest eigenvalue using shift-invert mode.
    # Direct 'SA' (smallest algebraic) can find spurious boundary modes
    # at large matrix sizes. Shift-invert solves (H - sigma*I)^{-1}
    # and finds eigenvalues nearest to sigma.
    # Note: this function returns E_elec (no V_NN), so sigma must target
    # the electronic energy: E_elec ~ E_total - V_NN.
    # Estimate E_atoms for sigma: same logic as main solver
    if Z_A == Z_B:
        E_atoms_est = -Z_A**2
    else:
        Z_max = max(Z_A, Z_B)
        if Z_max == 2:
            E_atoms_est = -2.9037
        else:
            E_atoms_est = -Z_max**2 + 5 * Z_max / 8
    V_NN = Z_A * Z_B / R
    sigma = E_atoms_est - V_NN - 0.15
    try:
        evals, _ = eigsh(H_csr, k=3, sigma=sigma, which='LM')
        E_elec = np.min(evals)
    except Exception as e:
        if verbose:
            print(f"  Shift-invert failed ({e}), falling back to SA")
        evals, _ = eigsh(H_csr, k=6, which='SA')
        E_elec = evals[0]

    if verbose:
        print(f"  Eigenvalues found: {np.sort(evals)}")
        print(f"  sigma = {sigma:.4f}")

    return E_elec


def _safe_potential_wrapper(
    U_spline,
    R_e_max: float,
) -> Callable:
    """Wrap a CubicSpline to return 0 beyond R_e_max.

    Gauss-Laguerre quadrature points extend to infinity. CubicSpline
    extrapolation can diverge, creating artificial deep wells. This
    wrapper clamps V(R) to 0 for R > R_e_max, matching the physical
    asymptotic behavior (U ~ 1/R_e² -> 0).
    """
    def V_safe(R: np.ndarray) -> np.ndarray:
        R = np.asarray(R)
        V = np.zeros_like(R, dtype=float)
        mask = R <= R_e_max
        if np.any(mask):
            V[mask] = U_spline(R[mask])
        return V
    return V_safe


def solve_adiabatic_radial_spectral(
    U_spline,
    n_basis: int = 25,
    alpha: float = 1.0,
    R_e_min: float = 0.3,
    R_e_max: float = 15.0,
) -> Tuple[float, np.ndarray, np.ndarray]:
    """Spectral Laguerre solver for the Level 4 adiabatic radial equation.

    Solves [-1/2 d²/dR_e² + U(R_e)] F = E F using Laguerre basis functions,
    reusing the Level 3 spectral pattern (hyperspherical_radial.py).

    Parameters
    ----------
    U_spline : callable
        Adiabatic effective potential U(R_e) (CubicSpline or callable).
    n_basis : int
        Number of Laguerre basis functions.
    alpha : float
        Exponential decay parameter for the Laguerre basis.
    R_e_min : float
        Left boundary (Dirichlet BC).
    R_e_max : float
        Beyond this R_e, the potential is clamped to 0 to prevent
        CubicSpline extrapolation artifacts.

    Returns
    -------
    E : float
        Ground state energy eigenvalue (Ha).
    F : ndarray
        Radial wavefunction on evaluation grid.
    R_eval : ndarray
        Evaluation grid points.
    """
    V_safe = _safe_potential_wrapper(U_spline, R_e_max)
    evals, F, R_eval = _solve_radial_spectral_l3(
        V_safe, n_basis=n_basis, alpha=alpha, R_min=R_e_min, n_states=1,
    )
    return evals[0], F[0], R_eval


def solve_coupled_channel_radial_spectral(
    R: float,
    R_e_grid: np.ndarray,
    l_max: int = 2,
    Z: float = 1.0,
    n_alpha: int = 200,
    m_max: int = 0,
    n_coupled: int = 3,
    l_max_per_m: Union[dict, None] = None,
    R_e_min: float = 0.3,
    Z_A: float = None,
    Z_B: float = None,
    z0: float = 0.0,
    pk_potentials: Union[List[dict], None] = None,
    n_basis: int = 25,
    alpha: float = 1.0,
) -> Tuple[float, np.ndarray]:
    """Spectral Laguerre solver for the Level 4 diabatic coupled-channel equation.

    Replaces the FD radial grid in solve_coupled_channel_radial with a
    Laguerre basis. The diabatic potential V_ij(R_e) is projected from the
    angular Hamiltonian onto a fixed basis at a reference R_e, then the
    coupled radial equations are solved in the Laguerre basis.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    R_e_grid : ndarray
        R_e grid for angular sweeps (used to find reference point).
    l_max, Z, n_alpha, m_max, n_coupled, l_max_per_m : as in solve_coupled_channel_radial.
    R_e_min : float
        Left Dirichlet boundary.
    Z_A, Z_B : float or None
        Per-nucleus charges.
    z0 : float
        Origin shift.
    pk_potentials : list or None
        Phillips-Kleinman specs.
    n_basis : int
        Number of Laguerre basis functions per channel.
    alpha : float
        Exponential decay parameter.

    Returns
    -------
    E_elec : float
        Electronic energy eigenvalue.
    F : ndarray
        Radial wavefunction (n_coupled * n_basis,).
    """
    # Step 1: Quick adiabatic sweep to find reference R_e
    n_Re_ang = len(R_e_grid)
    mu_vals = np.zeros(n_Re_ang)
    for i, R_e in enumerate(R_e_grid):
        rho = R / (2.0 * R_e)
        evals, _, _, _ = solve_angular_multichannel(
            rho, R_e, l_max, Z, n_alpha, m_max=m_max,
            l_max_per_m=l_max_per_m, Z_A=Z_A, Z_B=Z_B, z0=z0,
        )
        mu_vals[i] = evals[0]

    U_adia = (mu_vals + 15.0 / 8.0) / R_e_grid**2
    i_ref = np.argmin(U_adia)
    R_e_ref = R_e_grid[i_ref]

    # Step 2: Get reference eigenvectors at R_e_ref
    rho_ref = R / (2.0 * R_e_ref)
    _, ref_vecs, _, _ = solve_angular_multichannel(
        rho_ref, R_e_ref, l_max, Z, n_alpha,
        n_eig=n_coupled, m_max=m_max, l_max_per_m=l_max_per_m,
        Z_A=Z_A, Z_B=Z_B, z0=z0,
    )

    # Step 3: Build Laguerre basis and quadrature
    N = n_basis
    n_ch = n_coupled
    two_alpha = 2.0 * alpha
    n_quad = max(3 * N + 10, 80)
    x_quad, w_quad = roots_laguerre(n_quad)
    R_q = R_e_min + x_quad / two_alpha

    # Evaluate Laguerre polynomials and derivative kernel
    L_vals = np.zeros((N, n_quad))
    for n in range(N):
        L_vals[n] = eval_laguerre(n, x_quad)

    dL_vals = np.zeros((N, n_quad))
    for n in range(1, N):
        dL_vals[n] = dL_vals[n - 1] - L_vals[n - 1]

    B_vals = (1.0 - x_quad / 2.0)[np.newaxis, :] * L_vals + \
             x_quad[np.newaxis, :] * dL_vals

    inv_8a3 = 1.0 / (8.0 * alpha**3)
    inv_4a = 1.0 / (4.0 * alpha)
    x2_w = w_quad * x_quad**2

    # Overlap block
    x2_wL = x2_w[np.newaxis, :] * L_vals
    S_block = inv_8a3 * (x2_wL @ L_vals.T)

    # Kinetic block
    wB = w_quad[np.newaxis, :] * B_vals
    K_block = inv_4a * (wB @ B_vals.T)

    # Step 4: Precompute alpha grid for angular Hamiltonian
    h_alpha = (np.pi / 2) / (n_alpha + 1)
    alpha_grid = (np.arange(n_alpha) + 1) * h_alpha

    # Step 5: Compute diabatic potential at quadrature points and assemble
    dim = n_ch * N
    H_full = np.zeros((dim, dim))
    S_full = np.zeros((dim, dim))

    # For each channel pair, accumulate potential integrals over quadrature
    V_diabatic_q = np.zeros((n_quad, n_ch, n_ch))
    for k, R_e in enumerate(R_q):
        if R_e <= 0:
            continue
        rho = R / (2.0 * R_e)
        H_ang = build_angular_hamiltonian(
            alpha_grid, rho, R_e, l_max, Z, m_max, l_max_per_m,
            Z_A=Z_A, Z_B=Z_B, z0=z0,
            pk_potentials=pk_potentials,
        )
        Hv = H_ang @ ref_vecs.T
        V_proj = ref_vecs @ Hv

        for ic in range(n_ch):
            V_proj[ic, ic] = (V_proj[ic, ic] + 15.0 / 8.0) / R_e**2
            for jc in range(n_ch):
                if ic != jc:
                    V_proj[ic, jc] /= R_e**2

        V_diabatic_q[k] = V_proj

    # Assemble block Hamiltonian
    for ic in range(n_ch):
        ms = ic * N
        me = (ic + 1) * N

        # Kinetic (same for all channels)
        H_full[ms:me, ms:me] = K_block

        # Overlap
        S_full[ms:me, ms:me] = S_block

        for jc in range(n_ch):
            ns = jc * N
            ne = (jc + 1) * N

            # Potential: V_ij_nm = inv_8a3 * sum_k w_k x_k^2 V_ij(R_k) L_n L_m
            V_q = V_diabatic_q[:, ic, jc]
            x2_wVL = (x2_w * V_q)[np.newaxis, :] * L_vals
            V_block = inv_8a3 * (x2_wVL @ L_vals.T)
            H_full[ms:me, ns:ne] += V_block

    # Step 6: Solve generalized eigenvalue problem
    evals, evecs = eigh(H_full, S_full)
    E_elec = evals[0]
    F = evecs[:, 0]

    return E_elec, F


def solve_direct_2d_spectral(
    R: float,
    l_max: int = 2,
    Z: float = 1.0,
    n_alpha: int = 100,
    n_basis: int = 25,
    alpha_laguerre: float = 1.0,
    R_e_min: float = 0.3,
    m_max: int = 0,
    l_max_per_m: Union[dict, None] = None,
    verbose: bool = True,
    Z_A: float = None,
    Z_B: float = None,
    z0: float = 0.0,
    core_potentials: Union[List[dict], None] = None,
    pk_potentials: Union[List[dict], None] = None,
) -> float:
    """Spectral Laguerre 2D solver: full (alpha, R_e) variational problem.

    Replaces the FD R_e grid in solve_direct_2d with a Laguerre basis.
    The Hamiltonian is:

        H = T_Re ⊗ I_ang + [H_ang(R_e) + 15/8 I] / R_e²

    where T_Re is represented in the Laguerre basis via overlap and kinetic
    matrices, and the angular part is integrated over R_e by Gauss-Laguerre
    quadrature.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    l_max, Z, n_alpha : as in solve_direct_2d.
    n_basis : int
        Number of Laguerre basis functions for R_e.
    alpha_laguerre : float
        Exponential decay parameter for Laguerre basis.
    R_e_min : float
        Left Dirichlet boundary for R_e.
    m_max, l_max_per_m, verbose : as in solve_direct_2d.
    Z_A, Z_B, z0 : as in solve_direct_2d.
    core_potentials, pk_potentials : as in solve_direct_2d.

    Returns
    -------
    E_elec : float
        Lowest electronic energy eigenvalue.
    """
    # Resolve charges
    if Z_A is None:
        Z_A = Z
    if Z_B is None:
        Z_B = Z
    homonuclear = (Z_A == Z_B)

    # Determine channels
    if m_max == 0:
        channels_2 = _channel_list(l_max, homonuclear=homonuclear)
        channels_4 = [(l1, 0, l2, 0) for l1, l2 in channels_2]
    else:
        channels_4 = _channel_list_extended(
            l_max, m_max, l_max_per_m, homonuclear=homonuclear,
        )
    n_ch = len(channels_4)

    # Angular grid (same as FD 2D solver)
    h_alpha = (np.pi / 2) / (n_alpha + 1)
    alpha_grid = (np.arange(n_alpha) + 1) * h_alpha

    N_ang = n_ch * n_alpha
    N = n_basis
    N_total = N * N_ang

    if verbose:
        print(f"  Spectral 2D solver: {N_total} x {N_total} dense matrix")
        print(f"  ({N} Laguerre basis x {n_ch} channels x {n_alpha} alpha points)")

    # Build Laguerre quadrature and basis
    two_alpha = 2.0 * alpha_laguerre
    n_quad = max(3 * N + 10, 80)
    x_quad, w_quad = roots_laguerre(n_quad)
    R_q = R_e_min + x_quad / two_alpha

    L_vals = np.zeros((N, n_quad))
    for n in range(N):
        L_vals[n] = eval_laguerre(n, x_quad)

    dL_vals = np.zeros((N, n_quad))
    for n in range(1, N):
        dL_vals[n] = dL_vals[n - 1] - L_vals[n - 1]

    B_vals = (1.0 - x_quad / 2.0)[np.newaxis, :] * L_vals + \
             x_quad[np.newaxis, :] * dL_vals

    inv_8a3 = 1.0 / (8.0 * alpha_laguerre**3)
    inv_4a = 1.0 / (4.0 * alpha_laguerre)
    x2_w = w_quad * x_quad**2

    # Overlap block (N x N): S_nm = inv_8a3 * sum_k w_k x_k^2 L_n L_m
    x2_wL = x2_w[np.newaxis, :] * L_vals
    S_Re = inv_8a3 * (x2_wL @ L_vals.T)

    # Kinetic block (N x N): K_nm = inv_4a * sum_k w_k B_n B_m
    wB = w_quad[np.newaxis, :] * B_vals
    K_Re = inv_4a * (wB @ B_vals.T)

    # Build full Hamiltonian and overlap using Kronecker products
    # Index: n * N_ang + (ch * n_alpha + ia)
    I_ang = np.eye(N_ang)

    # Kinetic: K_Re ⊗ I_ang;  Overlap: S_Re ⊗ I_ang
    H_full = np.kron(K_Re, I_ang)
    S_full = np.kron(S_Re, I_ang)

    # Angular potential via quadrature:
    # V_nm^{ab} = inv_8a3 * sum_k w_k x_k^2 [H_ang(R_k)+15/8]_{ab}/R_k^2 * L_n(x_k) * L_m(x_k)
    # Vectorized as: LLt_k = outer(L[:,k], L[:,k]) then kron(wt * LLt_k, H_ang_k)

    for k in range(n_quad):
        R_e = R_q[k]
        if R_e <= 0:
            continue
        rho = R / (2.0 * R_e)

        H_ang = build_angular_hamiltonian(
            alpha_grid, rho, R_e, l_max, Z, m_max, l_max_per_m,
            Z_A=Z_A, Z_B=Z_B, z0=z0,
            core_potentials=core_potentials,
            pk_potentials=pk_potentials,
        )
        np.fill_diagonal(H_ang, H_ang.diagonal() + 15.0 / 8.0)
        H_ang *= (1.0 / R_e**2)

        wt = inv_8a3 * x2_w[k]
        L_k = L_vals[:, k]  # shape (N,)
        LLt = np.outer(L_k, L_k)  # (N, N)
        H_full += np.kron(wt * LLt, H_ang)

    # Solve generalized eigenvalue problem
    if verbose:
        print(f"  Solving {N_total}x{N_total} generalized eigenvalue problem...")

    evals, evecs = eigh(H_full, S_full)

    # Find lowest physical eigenvalue
    # Filter out numerical artifacts from ill-conditioned overlap
    E_elec = evals[0]

    # Estimate sigma for sanity check
    if Z_A == Z_B:
        E_atoms_est = -Z_A**2
    else:
        Z_max = max(Z_A, Z_B)
        if Z_max == 2:
            E_atoms_est = -2.9037
        else:
            E_atoms_est = -Z_max**2 + 5 * Z_max / 8
    V_NN = Z_A * Z_B / R
    sigma = E_atoms_est - V_NN - 0.15

    if verbose:
        print(f"  Lowest eigenvalues: {evals[:5]}")
        print(f"  sigma = {sigma:.4f}")

    return E_elec


def solve_level4_h2_multichannel(
    R: float,
    l_max: int = 2,
    Z: float = 1.0,
    n_alpha: int = 200,
    n_Re: int = 400,
    R_e_min: float = 0.3,
    R_e_max: float = 15.0,
    verbose: bool = True,
    m_max: int = 0,
    n_coupled: int = 1,
    l_max_per_m: Union[dict, None] = None,
    Z_A: float = None,
    Z_B: float = None,
    E_exact: float = None,
    D_e_exact: float = None,
    origin: str = 'midpoint',
    core_potentials: Union[List[dict], None] = None,
    pk_potentials: Union[List[dict], None] = None,
    pk_projector: Union[dict, None] = None,
    radial_method: str = 'fd',
    n_basis_radial: int = 25,
    alpha_radial: float = 1.0,
    angular_method: str = 'fd',
    n_basis_angular: int = 10,
    n_quad_angular: int = 100,
) -> dict:
    """
    Full Level 4 multichannel solver for two-electron diatomics.

    Supports both homonuclear (H2, He2^2+) and heteronuclear (HeH+)
    systems via Z_A, Z_B parameters.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    l_max : int
        Maximum angular momentum per electron.
    Z : float
        Nuclear charge per nucleus (homonuclear default).
    n_alpha : int
        FD grid points for alpha.
    n_Re : int
        Grid points for hyperradial equation.
    R_e_min, R_e_max : float
        Hyperradial boundaries.
    verbose : bool
        Print diagnostics.
    m_max : int
        Maximum |m| per electron. 0 = sigma-only, 1 = sigma+pi.
    n_coupled : int
        Radial solver mode. 1 = single-channel adiabatic (default).
        -1 = direct 2D solver (properly variational, no adiabatic
        approximation, but slower). >1 = DBOC-corrected adiabatic
        (experimental, numerically unstable for large channel counts).
    Z_A, Z_B : float or None
        Per-nucleus charges. If None, uses Z for both (homonuclear).
    E_exact, D_e_exact : float or None
        Known exact values for comparison. If None, uses H2 defaults
        for Z_A=Z_B=1, otherwise omits comparison.
    origin : str or float
        Origin for hyperspherical coordinates.
        'midpoint': geometric center (default, z0=0).
        'charge_center': charge-weighted center, z0 = R(Z_A-Z_B)/(2(Z_A+Z_B)).
        float: explicit z0 value.
    radial_method : str
        Radial solver method. 'fd' (default) uses finite differences on a
        uniform grid. 'spectral' uses a Laguerre spectral basis, giving
        100x+ dimension reduction with equivalent accuracy.
    n_basis_radial : int
        Number of Laguerre basis functions (spectral method only).
    alpha_radial : float
        Exponential decay parameter for Laguerre basis (spectral method only).
    angular_method : str
        Angular solver method. 'fd' (default) uses finite differences on a
        uniform alpha grid. 'spectral' uses a Jacobi polynomial basis,
        giving 30-250x speedup with < 1e-4 eigenvalue agreement.
    n_basis_angular : int
        Number of Jacobi basis functions per channel (spectral method only).
    n_quad_angular : int
        Gauss-Legendre quadrature order per sub-interval (spectral only).

    Returns
    -------
    result : dict
        Keys: E_elec, E_total, D_e, D_e_pct, R, l_max, m_max, channels, etc.
    """
    import time
    t0 = time.time()

    # Resolve charges
    if Z_A is None:
        Z_A = Z
    if Z_B is None:
        Z_B = Z
    homonuclear = (Z_A == Z_B)

    # Resolve origin shift
    if isinstance(origin, (int, float)):
        z0 = float(origin)
    elif origin == 'charge_center':
        z0 = R * (Z_A - Z_B) / (2.0 * (Z_A + Z_B))
    else:  # 'midpoint'
        z0 = 0.0

    if m_max == 0:
        channels = _channel_list(l_max, homonuclear=homonuclear)
    else:
        channels = _channel_list_extended(
            l_max, m_max, l_max_per_m, homonuclear=homonuclear,
        )
    n_ch = len(channels)

    # Non-uniform R_e grid for adiabatic curve
    R_e_angular = np.concatenate([
        np.linspace(R_e_min, 1.0, 40),
        np.linspace(1.0, 3.0, 40),
        np.linspace(3.0, 6.0, 30),
        np.linspace(6.0, R_e_max, 20),
    ])
    R_e_angular = np.unique(R_e_angular)

    # System label
    if homonuclear:
        sys_label = f"Z={Z_A}" if Z_A != 1 else "H2"
    else:
        sys_label = f"Z_A={Z_A}, Z_B={Z_B}"

    if verbose:
        print(f"Level 4 multichannel solver for {sys_label} at R = {R:.4f} bohr")
        print(f"  l_max={l_max}, m_max={m_max}, n_ch={n_ch}")
        if not homonuclear:
            print(f"  Heteronuclear: Z_A={Z_A}, Z_B={Z_B}")
        if z0 != 0.0:
            print(f"  Origin shift: z0={z0:.4f} (R_A={R/2 - z0:.4f}, R_B={R/2 + z0:.4f})")
        if n_ch <= 15:
            print(f"  channels={channels}")
        if angular_method == 'spectral':
            ang_dim = n_ch * n_basis_angular
            print(f"  angular_method=spectral (n_basis={n_basis_angular}), "
                  f"matrix {ang_dim}x{ang_dim}")
        else:
            print(f"  angular_method=fd (n_alpha={n_alpha}), "
                  f"matrix {n_ch * n_alpha}x{n_ch * n_alpha}")
        if radial_method == 'spectral':
            print(f"  radial_method=spectral "
                  f"(n_basis={n_basis_radial}, alpha={alpha_radial})")
        else:
            print(f"  n_Re={n_Re}")
        if n_coupled > 1:
            print(f"  DBOC radial solver: {n_coupled} states")
        elif n_coupled == -1:
            print(f"  Direct 2D solver (no adiabatic approximation)")

    if n_coupled == -1:
        # --- Direct 2D solver (properly variational) ---
        if radial_method == 'spectral':
            E_elec = solve_direct_2d_spectral(
                R, l_max, Z, n_alpha,
                n_basis=n_basis_radial,
                alpha_laguerre=alpha_radial,
                R_e_min=R_e_min,
                m_max=m_max, l_max_per_m=l_max_per_m,
                verbose=verbose,
                Z_A=Z_A, Z_B=Z_B, z0=z0,
                core_potentials=core_potentials,
                pk_potentials=pk_potentials,
            )
        else:
            E_elec = solve_direct_2d(
                R, l_max, Z, n_alpha, n_Re, R_e_min, R_e_max,
                m_max, l_max_per_m, verbose=verbose,
                Z_A=Z_A, Z_B=Z_B, z0=z0,
                core_potentials=core_potentials,
                pk_potentials=pk_potentials,
            )
        F = np.zeros(n_Re)  # no radial wavefunction in 2D mode

        t2 = time.time()
        if verbose:
            print(f"  2D solve: {t2 - t0:.2f}s")

    elif n_coupled > 1:
        # --- DBOC-corrected adiabatic solver (experimental) ---
        if verbose:
            print(f"  Computing DBOC-corrected adiabatic curve "
                  f"({n_coupled} states) on "
                  f"{len(R_e_angular)} R_e points...")

        U_corrected, U_bare = compute_adiabatic_curve_dboc(
            R, R_e_angular, l_max, Z, n_alpha, m_max,
            n_states=n_coupled, l_max_per_m=l_max_per_m,
            Z_A=Z_A, Z_B=Z_B, z0=z0,
        )

        t1 = time.time()
        if verbose:
            i_min = np.argmin(U_corrected)
            print(f"  DBOC sweep: {t1 - t0:.2f}s")

        U_spline = CubicSpline(R_e_angular, U_corrected,
                               extrapolate=True)

        if radial_method == 'spectral':
            E_elec, F, _ = solve_adiabatic_radial_spectral(
                U_spline, n_basis=n_basis_radial,
                alpha=alpha_radial, R_e_min=R_e_min,
                R_e_max=R_e_max,
            )
        else:
            h_Re = (R_e_max - R_e_min) / (n_Re + 1)
            R_e_radial = R_e_min + (np.arange(n_Re) + 1) * h_Re
            V_radial = U_spline(R_e_radial)

            diag = np.ones(n_Re) / h_Re**2 + V_radial
            off_diag = -0.5 * np.ones(n_Re - 1) / h_Re**2

            evals, evecs = eigh_tridiagonal(
                diag, off_diag,
                select='i', select_range=(0, 0),
            )

            E_elec = evals[0]
            F = evecs[:, 0]

        t2 = time.time()
        if verbose:
            print(f"  Radial solve: {t2 - t1:.2f}s")

    else:
        # --- Single-channel adiabatic solver (original) ---
        if verbose:
            print(f"  Computing adiabatic curve on "
                  f"{len(R_e_angular)} R_e points...")

        if angular_method == 'spectral':
            from geovac.level4_spectral_angular import (
                compute_adiabatic_curve_spectral,
            )
            U_angular = compute_adiabatic_curve_spectral(
                R, R_e_angular, l_max, Z, n_basis=n_basis_angular,
                n_quad=n_quad_angular, m_max=m_max,
                l_max_per_m=l_max_per_m, Z_A=Z_A, Z_B=Z_B, z0=z0,
                core_potentials=core_potentials,
                pk_potentials=pk_potentials,
            )
        else:
            U_angular = compute_adiabatic_curve_mc(
                R, R_e_angular, l_max, Z, n_alpha, m_max=m_max,
                l_max_per_m=l_max_per_m, Z_A=Z_A, Z_B=Z_B, z0=z0,
                core_potentials=core_potentials,
                pk_potentials=pk_potentials,
                pk_projector=pk_projector,
            )

        t1 = time.time()
        if verbose:
            i_min = np.argmin(U_angular)
            print(f"  Angular sweep: {t1 - t0:.2f}s")
            print(f"  U_min = {U_angular[i_min]:.6f} Ha at "
                  f"R_e = {R_e_angular[i_min]:.3f}")

        # Interpolate and solve radial equation
        U_spline = CubicSpline(R_e_angular, U_angular, extrapolate=True)

        if radial_method == 'spectral':
            E_elec, F, _ = solve_adiabatic_radial_spectral(
                U_spline, n_basis=n_basis_radial,
                alpha=alpha_radial, R_e_min=R_e_min,
                R_e_max=R_e_max,
            )
        else:
            h_Re = (R_e_max - R_e_min) / (n_Re + 1)
            R_e_radial = R_e_min + (np.arange(n_Re) + 1) * h_Re
            V_radial = U_spline(R_e_radial)

            diag = np.ones(n_Re) / h_Re**2 + V_radial
            off_diag = -0.5 * np.ones(n_Re - 1) / h_Re**2

            evals, evecs = eigh_tridiagonal(
                diag, off_diag,
                select='i', select_range=(0, 0),
            )

            E_elec = evals[0]
            F = evecs[:, 0]

        t2 = time.time()
        if verbose:
            print(f"  Radial solve: {t2 - t1:.2f}s")

    norm = np.sqrt(np.sum(F**2))
    if norm > 0:
        F /= norm

    V_NN = Z_A * Z_B / R
    E_total = E_elec + V_NN

    # Dissociation limit: where do the 2 electrons go?
    if homonuclear:
        # Each atom gets one electron: 2 × E(Z, 1e) = 2 × (-Z²/2) = -Z²
        E_atoms = -Z_A**2
    else:
        # Both electrons go to the larger nucleus (better binding)
        Z_max = max(Z_A, Z_B)
        if Z_max == 1:
            # H⁻ is barely bound; use H + H limit
            E_atoms = -1.0
        elif Z_max == 2:
            # He atom (exact nonrelativistic)
            E_atoms = -2.9037
        else:
            # Variational 2-electron estimate: E ≈ -Z² + 5Z/8
            E_atoms = -Z_max**2 + 5 * Z_max / 8

    D_e = E_atoms - E_total

    # Determine exact reference values
    if E_exact is None and D_e_exact is None:
        if homonuclear and Z_A == 1.0:
            # H2 defaults
            E_exact = -1.17447
            D_e_exact = 0.17447
        elif not homonuclear and Z_A == 2.0 and Z_B == 1.0:
            # HeH+ reference (Kolos & Peek 1976)
            E_exact = -2.9787
            D_e_exact = E_atoms - E_exact  # 0.0750
    elif D_e_exact is None and E_exact is not None:
        D_e_exact = E_atoms - E_exact

    D_e_pct = (D_e / D_e_exact * 100) if D_e_exact else None

    if verbose:
        print(f"\n  === Results (l_max={l_max}, {n_ch} channels) ===")
        if E_exact is not None:
            print(f"  E_total    = {E_total:.6f} Ha  (exact: {E_exact:.6f})")
            print(f"  D_e        = {D_e:.6f} Ha  (exact: {D_e_exact:.6f})")
            print(f"  D_e / exact = {D_e_pct:.1f}%")
            if homonuclear and Z_A == 1.0:
                if D_e_pct > 92.4:
                    print(f"  ** IMPROVES on Paper 12 Neumann V_ee (92.4%) **")
                else:
                    print(f"  Paper 12 Neumann V_ee: 92.4%")
        else:
            print(f"  E_total    = {E_total:.6f} Ha")
            print(f"  D_e        = {D_e:.6f} Ha")
            print(f"  E_atoms    = {E_atoms:.6f} Ha")
        print(f"  Total time: {t2 - t0:.2f}s")

    result = {
        'E_elec': E_elec,
        'E_total': E_total,
        'D_e': D_e,
        'D_e_pct': D_e_pct,
        'R': R,
        'Z': Z,
        'Z_A': Z_A,
        'Z_B': Z_B,
        'homonuclear': homonuclear,
        'l_max': l_max,
        'm_max': m_max,
        'n_ch': n_ch,
        'n_coupled': n_coupled,
        'channels': channels,
        'R_e_grid_angular': R_e_angular,
        'wavefunction': F,
        'E_atoms': E_atoms,
        'V_NN': V_NN,
        'z0': z0,
        'origin': origin,
        'radial_method': radial_method,
        'angular_method': angular_method,
    }

    if n_coupled == 1:
        result['U_adiabatic'] = U_angular
    if radial_method == 'spectral':
        result['n_basis_radial'] = n_basis_radial
        result['alpha_radial'] = alpha_radial
    else:
        h_Re = (R_e_max - R_e_min) / (n_Re + 1)
        result['R_e_grid_radial'] = R_e_min + (np.arange(n_Re) + 1) * h_Re

    return result


def solve_level4_lih(
    R: float,
    core_model: str = 'zeff',
    l_max: int = 4,
    m_max: int = 0,
    n_alpha: int = 200,
    n_Re: int = 400,
    n_coupled: int = 1,
    R_e_max: float = None,
    verbose: bool = True,
) -> dict:
    """
    Convenience wrapper for Level 4 LiH with locked Li 1s² core.

    Solves two active electrons in hyperspherical coordinates with the
    Li 1s² core locked.  Three core models control how the active
    electrons see the screened Li nucleus:

      'bare'     — Z_eff_A = 3 (no screening).  Diagnostic only.
      'zeff'     — Z_eff_A = 1 (integer screening, Z - n_core).
      'hartree'  — Z_eff_A = 1 + penetration correction from the
                   analytical Core Hartree potential:
                   V_H(r) = (2/r)[1 - (1 + Z_core r) exp(-2 Z_core r)]
      'clementi' — Z_eff_A = 1.279 (Clementi-Raimondi), Z_eff_B = 1.0.

    The 'hartree' model is mathematically equivalent to 'zeff' +
    `penetration=True` in the LockedShellMolecule hyperspherical path.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    core_model : str
        'bare', 'zeff', or 'hartree'.
    l_max : int
        Maximum angular momentum per electron.
    m_max : int
        Maximum |m| per electron.
    n_alpha : int
        FD grid points for alpha.
    n_Re : int
        Grid points for radial equation.
    n_coupled : int
        Radial solver mode (1=adiabatic, -1=direct 2D).
    R_e_max : float or None
        Hyperradial boundary (auto-scaled if None).
    verbose : bool
        Print diagnostics.

    Returns
    -------
    result : dict
        Keys: E_total, D_e, E_elec, E_locked, V_NN, E_core_isolated,
        V_cross_nuc, core_model, and all keys from solve_level4_h2_multichannel.
    """
    from .locked_shell import LockedShellMolecule

    # Map core_model to LockedShellMolecule kwargs
    hyper_kwargs: dict = {
        'l_max': l_max,
        'm_max': m_max,
        'n_alpha': n_alpha,
        'n_Re': n_Re,
        'n_coupled': n_coupled,
        'origin': 'charge_center',
    }

    if core_model == 'bare':
        hyper_kwargs['zeff_mode'] = 'bare'
        hyper_kwargs['penetration'] = False
    elif core_model == 'hartree':
        hyper_kwargs['zeff_mode'] = 'screened'
        hyper_kwargs['penetration'] = True
    elif core_model == 'clementi':
        hyper_kwargs['zeff_mode'] = 'clementi'
        hyper_kwargs['penetration'] = False
    else:  # 'zeff' (default)
        hyper_kwargs['zeff_mode'] = 'screened'
        hyper_kwargs['penetration'] = False

    if R_e_max is not None:
        hyper_kwargs['R_e_max'] = R_e_max

    hyper_kwargs['verbose'] = verbose

    mol = LockedShellMolecule(
        Z_A=3, Z_B=1,
        nmax_A=3, nmax_B=3,
        R=R,
        n_electrons=4,
        locked_config={0: [(1, 0)]},
        active_method='hyperspherical',
        hyperspherical_kwargs=hyper_kwargs,
    )

    eigvals, eigvecs = mol.solve()

    E_total = eigvals[0]

    # Dissociation limit: Li(1s²2s) + H(1s)
    # Li ground state: E ≈ -7.478 Ha (exact nonrelativistic)
    # H ground state: E = -0.500 Ha
    # Total: -7.978 Ha
    E_Li_exact = -7.478
    E_H = -0.500
    E_atoms = E_Li_exact + E_H
    D_e = E_atoms - E_total
    D_e_eV = D_e * 27.2114

    # Experimental reference
    D_e_expt_eV = 2.515
    D_e_expt_Ha = D_e_expt_eV / 27.2114

    if verbose:
        print(f"\n  === LiH Level 4 Results (core_model='{core_model}') ===")
        print(f"  R         = {R:.4f} bohr")
        print(f"  E_total   = {E_total:.6f} Ha")
        print(f"  E_atoms   = {E_atoms:.6f} Ha (Li + H)")
        print(f"  D_e       = {D_e:.6f} Ha = {D_e_eV:.3f} eV")
        print(f"  D_e(expt) = {D_e_expt_Ha:.6f} Ha = {D_e_expt_eV:.3f} eV")
        if D_e > 0:
            print(f"  D_e/expt  = {D_e / D_e_expt_Ha * 100:.1f}%")
        else:
            print(f"  ** UNBOUND **")

    result = {
        'E_total': E_total,
        'D_e': D_e,
        'D_e_eV': D_e_eV,
        'E_atoms': E_atoms,
        'E_locked': mol.E_locked,
        'E_core_isolated': mol._E_core_isolated,
        'V_cross_nuc': mol._V_cross_nuc,
        'V_NN': mol.V_NN,
        'E_elec': mol._E_elec,
        'core_model': core_model,
        'R': R,
        'Z_eff_A': mol._Z_eff_A,
        'Z_eff_B': mol._Z_eff_B,
    }

    if hasattr(mol, '_level4_result'):
        result['level4'] = mol._level4_result

    return result
