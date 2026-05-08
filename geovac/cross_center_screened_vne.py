"""
Screened Cross-Center Nuclear Attraction Integrals (Phase C-W1c)
================================================================

Computes screened cross-center nuclear attraction matrix elements for
frozen-core species:

    I^{AB}_{nlm,n'l'm'} = <psi_{nlm}^A | -Z_eff^B(|r-R_B|) / |r-R_B| | psi_{n'l'm'}^A>

where psi_{nlm}^A are hydrogenic orbitals centered on nucleus A, and
Z_eff^B(rho) is the FrozenCore-screened effective charge of nucleus B
seen at radial distance rho from B.

This is the screened analog of ``shibuya_wulfman.compute_cross_center_vne``,
which uses the bare nuclear charge Z_B regardless of frozen core. For
second-row hydrides (NaH, MgH2, HCl, ...) the bare-charge cross-center V_ne
overattracts the H-side orbitals by ~10x at typical bond lengths because
the [Ne] core is fully internalized at R > 1.5 bohr (Phase B-W1c-diag,
``debug/multifocal_b_w1c_diag_memo.md``).

Algebraic decomposition (Phase B-W1c-diag, Section 4):

    V_screened(r, R_B) = -Z_eff^B(rho) / rho - tilde_Phi^B(rho)
                      = -Z_B / rho + N_core^B(rho) / rho + Phi_outer^B(rho)
                      = bare Coulomb - Hartree potential of core electrons

where rho = |r - R_B|, both terms are spherically symmetric about B (Newton's
shell theorem on a spherical core density), and tilde_Phi^B is the smooth
"outer" Hartree integral that is finite at the origin.

The multipole expansion in A-centered spherical coordinates terminates at
L_max = 2 * l_max by the Gaunt 3j triangle inequality (which acts on the
orbital angular labels, NOT on the radial form of the potential). The
angular machinery (``_angular_coefficient``) is reused verbatim from the
bare-Coulomb implementation.

The radial multipole pieces f_L(r, R_B) are computed by direct Legendre
projection:

    f_L(r, R_B) = (2L + 1) / 2 * integral_{-1}^{1} f(rho(u)) * P_L(u) du

with rho(u) = sqrt(r^2 + R_B^2 - 2*r*R_B*u). For the bare Coulomb part,
this reduces analytically to r_<^L / r_>^{L+1}; for the smooth outer part,
it is computed by Gauss-Legendre quadrature. We split the screened
potential into:

    f(rho) = -Z_eff^B(rho) / rho - tilde_Phi^B(rho)

and treat the two pieces separately:

  - The bare Coulomb -Z_B/rho is handled analytically (existing path).
  - The screening correction +N_core^B(rho)/rho - tilde_Phi^B(rho) is
    smooth, exponentially decaying at large rho, and finite at rho=0.
    This piece is integrated numerically against P_L(u) by Gauss-Legendre.

Honest scope: the radial integrals are grid-quadrature, not closed-form
incomplete-gamma sums. The B-W1c-diag memo identifies a fully algebraic
path via Clementi-Raimondi exponential-shell decomposition; the present
module uses the simpler numerical-quadrature approach as a first-pass
implementation that exactly reproduces the bare result in the limit
Z_eff^B(rho) -> Z_B (no core) and converges with grid resolution.

Author: GeoVac Development Team (Phase C-W1c)
Date: May 2026
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np

from geovac.composed_qubit import _radial_wf_grid
from geovac.neon_core import FrozenCore
from geovac.shibuya_wulfman import (
    _angular_coefficient,
    _build_block_rotation_matrix,
    _radial_split_integral,
    compute_cross_center_vne,
)


# ---------------------------------------------------------------------------
# Frozen-core auto-detection
# ---------------------------------------------------------------------------


def _detect_core_type(Z_nuc: int) -> Optional[str]:
    """Return the FrozenCore core_type for Z, or None if no frozen core.

    Same auto-detection rules as ``FrozenCore.__init__``.
    Returns None for first-row Z = 1..10 (no frozen core).
    """
    Z_int = int(round(Z_nuc))
    if Z_int < 11:
        return None
    # Defer to FrozenCore's registry
    try:
        fc = FrozenCore(Z=Z_int)
    except ValueError:
        return None
    return fc.core_type


# ---------------------------------------------------------------------------
# Screening profile evaluation
# ---------------------------------------------------------------------------


def _screening_correction_grid(
    Z_nuc: float,
    fc: FrozenCore,
    rho_grid: np.ndarray,
) -> np.ndarray:
    """Evaluate the screening correction f_screen(rho) = -V_screened + V_bare.

    V_screened(rho) = -Z_eff^B(rho) / rho - tilde_Phi^B(rho)
    V_bare(rho)     = -Z_B / rho

    The screening correction added to V_bare to get V_screened is:

        f_screen(rho) = V_screened(rho) - V_bare(rho)
                      = -Z_eff^B(rho)/rho + Z_B/rho - tilde_Phi^B(rho)
                      = (Z_B - Z_eff^B(rho))/rho - tilde_Phi^B(rho)
                      = N_core^B(rho)/rho - tilde_Phi^B(rho)

    where tilde_Phi^B(rho) = integral_rho^infty n_r^B(rho')/rho' drho'
    is the smooth outer Hartree integral (n_r is FrozenCore's internal
    radial density, already including the 4*pi*r^2 volume element so that
    integral_0^infty n_r dr = N_core).

    This function is finite everywhere:
      - at rho=0: N_core(0)=0 and tilde_Phi(0) is finite
      - at rho=infty: both terms decay (N_core->N_core_total, divided
                      by rho gives 1/rho decay, but THIS goes back into
                      Z_eff_total/rho — note the trick is f_screen is
                      finite at the ORIGIN, not at infinity, where it
                      tends to (Z_B - Z_eff_asymptotic)/rho which is
                      the cumulative core contribution)

    For frozen-core species at typical bond R = 2-5 bohr, rho ranges widely
    over the orbital extent on side A; we evaluate on the supplied grid.
    """
    rho_grid = np.asarray(rho_grid, dtype=float)

    # Cumulative core charge inside rho:
    # FrozenCore's _N_core_spline returns the integral of n_r from 0 to rho.
    # n_r is the radial number density (already includes 4*pi*r^2).
    # Therefore N_core(rho) = integral_0^rho n_r(rho') drho' is exactly
    # the count of core electrons inside radius rho.
    N_core_at_rho = np.clip(
        fc._N_core_spline(rho_grid), 0.0, float(fc.n_core_electrons),
    )

    # Outer Hartree integral: tilde_Phi(rho) = integral_rho^infty n_r(rho') / rho' drho'
    # Compute on FrozenCore's internal r_grid where n_r is known, then interpolate.
    r_int = fc.r_grid
    n_r_int = fc.n_r
    # Avoid division by zero at the smallest grid point
    integrand = np.zeros_like(n_r_int)
    mask = r_int > 1e-12
    integrand[mask] = n_r_int[mask] / r_int[mask]

    # Cumulative integral from rho to infinity:
    # tilde_Phi(rho_i) = sum_{j>=i} integrand[j] * dr
    # Use cumulative_trapezoid in REVERSE direction.
    from scipy.integrate import cumulative_trapezoid
    cum_rev = cumulative_trapezoid(integrand[::-1], r_int[::-1], initial=0)
    # cum_rev is integrated against -dr (reversed grid); flip back
    tilde_Phi_int = -cum_rev[::-1]
    # tilde_Phi_int[i] = integral_{r_int[i]}^{r_int[-1]} n_r(rho')/rho' drho'
    # The truncation at r_max may underestimate; add a tail correction
    # assuming exponential decay (safe for well-internalized cores at r_max=20)
    # In practice, n_r at r_int[-1] is essentially zero for [Ne]/[Ar]/etc., so
    # tail correction is negligible at default r_max=20 bohr.

    from scipy.interpolate import CubicSpline
    tilde_Phi_spline = CubicSpline(r_int, tilde_Phi_int, extrapolate=True)
    tilde_Phi_at_rho = np.where(
        rho_grid > r_int[-1],
        0.0,  # extrapolate to zero beyond grid (n_core fully internalized)
        tilde_Phi_spline(rho_grid),
    )
    # Ensure non-negative (tilde_Phi is a positive Hartree-like integral)
    tilde_Phi_at_rho = np.clip(tilde_Phi_at_rho, 0.0, None)

    # f_screen(rho) = N_core(rho)/rho - tilde_Phi(rho)
    # At rho=0: N_core(0)=0, but N_core(rho)/rho -> 4*pi*n_core(0)*rho/rho
    # which is finite. We can use the limit form for very small rho:
    # near origin, N_core(rho) ~ rho * n_r(rho_min)/rho (linear).
    # Practical: clip rho_grid to a small min and the result is finite.
    rho_safe = np.where(rho_grid > 1e-10, rho_grid, 1e-10)
    f_screen = N_core_at_rho / rho_safe - tilde_Phi_at_rho

    return f_screen


# ---------------------------------------------------------------------------
# Multipole projection of the screening correction
# ---------------------------------------------------------------------------


def _f_L_screening_grid(
    fc: FrozenCore,
    Z_nuc: float,
    r_grid: np.ndarray,
    R_AB: float,
    L: int,
    n_legendre: int = 64,
) -> np.ndarray:
    """Compute f_screen_L(r, R_B) on the radial grid.

        f_screen_L(r, R_B) = (2L+1)/2 * integral_{-1}^{1}
                                f_screen(rho(u)) * P_L(u) du

    where rho(u) = sqrt(r^2 + R_B^2 - 2*r*R_B*u).

    Uses Gauss-Legendre quadrature in u with n_legendre nodes.

    Returns array of same shape as r_grid.
    """
    r_grid = np.asarray(r_grid, dtype=float)

    # Gauss-Legendre nodes and weights on [-1, 1]
    u_nodes, u_weights = np.polynomial.legendre.leggauss(n_legendre)

    # Legendre P_L evaluated at nodes
    P_L_nodes = np.zeros(n_legendre)
    for i, u in enumerate(u_nodes):
        # Use scipy or direct evaluation. For small L, direct is fine.
        # numpy.polynomial.legendre uses [c_0, c_1, ...] expansion.
        coefs = np.zeros(L + 1)
        coefs[L] = 1.0
        P_L_nodes[i] = np.polynomial.legendre.legval(u, coefs)

    # For each r, compute rho(u) at all nodes, evaluate f_screen, sum.
    # Vectorize: shape (Nr, n_legendre)
    r_col = r_grid[:, np.newaxis]
    u_row = u_nodes[np.newaxis, :]
    rho_grid_2d = np.sqrt(r_col ** 2 + R_AB ** 2 - 2.0 * r_col * R_AB * u_row)

    # Flatten, evaluate, reshape
    rho_flat = rho_grid_2d.reshape(-1)
    f_screen_flat = _screening_correction_grid(Z_nuc, fc, rho_flat)
    f_screen_2d = f_screen_flat.reshape(rho_grid_2d.shape)

    # Quadrature sum
    integrand = f_screen_2d * (u_weights * P_L_nodes)[np.newaxis, :]
    f_L_screening = (2.0 * L + 1.0) / 2.0 * np.sum(integrand, axis=1)

    return f_L_screening


# ---------------------------------------------------------------------------
# Radial integral with screening
# ---------------------------------------------------------------------------


def _screened_radial_integral(
    Z_orb: float,
    n1: int, l1: int,
    n2: int, l2: int,
    fc: FrozenCore,
    Z_nuc: float,
    L: int,
    R_AB: float,
    n_grid: int = 4000,
    r_max: Optional[float] = None,
    n_legendre: int = 64,
) -> Tuple[float, float]:
    """Compute the L-th radial integral pieces of the screened cross-center V_ne.

    Returns (rad_bare, rad_screen):
      rad_bare = bare Coulomb radial integral (existing analytical path)
      rad_screen = screening correction radial integral (numerical)

    The total radial integral for V_screened = -Z_eff/rho - tilde_Phi is:
      rad_total = -Z_nuc * rad_bare + rad_screen
                = -Z_nuc * rad_bare
                + integral_0^infty R_a(r) R_b(r) f_screen_L(r, R_B) r^2 dr

    NOTE: rad_bare here is the radial part WITHOUT the -Z_nuc factor;
    f_screen_L is the L-th projection of f_screen(rho), which is
    (V_screened - V_bare). So total V_ne = -Z_nuc * rad_bare + rad_screen.
    """
    # 1. Bare Coulomb radial integral (analytical, exact) — same as
    #    shibuya_wulfman._radial_split_integral
    rad_bare = _radial_split_integral(
        Z_orb, n1, l1, n2, l2, L, R_AB, n_grid=n_grid,
    )

    # 2. Screening correction radial integral (numerical quadrature)
    if r_max is None:
        # Choose r_max from orbital extent: 80/Z is canonical for hydrogenic
        # orbitals, but make sure to extend several bohr past R_AB.
        r_max = max(80.0 / max(Z_orb, 0.5), 4.0 * R_AB, 30.0)

    # Use uniform radial grid (np.linspace excludes 0)
    r_full = np.linspace(0, r_max, n_grid + 1)[1:]

    R1 = _radial_wf_grid(Z_orb, n1, l1, r_full)
    R2 = _radial_wf_grid(Z_orb, n2, l2, r_full)

    # f_screen_L on r-grid
    f_L_grid = _f_L_screening_grid(
        fc, Z_nuc, r_full, R_AB, L, n_legendre=n_legendre,
    )

    integrand = R1 * R2 * f_L_grid * r_full ** 2
    rad_screen = float(np.trapezoid(integrand, r_full))

    return rad_bare, rad_screen


# ---------------------------------------------------------------------------
# Public API: matrix element and full matrix
# ---------------------------------------------------------------------------


def compute_screened_cross_center_vne_element(
    Z_orb: float,
    n1: int, l1: int, m1: int,
    n2: int, l2: int, m2: int,
    Z_nuc: float,
    R_AB: float,
    L_max: int,
    core_type: Optional[str] = None,
    n_grid: int = 4000,
    n_legendre: int = 64,
    nuc_parity: int = 1,
) -> float:
    """Single matrix element of the screened cross-center nuclear attraction.

    <psi_{n1,l1,m1}^A | -Z_eff^B(|r-R_B|) / |r-R_B| | psi_{n2,l2,m2}^A>

    Drop-in analog of ``shibuya_wulfman.compute_cross_center_vne_element``.
    For first-row Z_nuc (no frozen core), exactly reproduces the bare result.

    Parameters
    ----------
    Z_orb : float
        Nuclear charge for the orbital wavefunctions (on side A).
    n1, l1, m1 : int
        Quantum numbers for the bra orbital.
    n2, l2, m2 : int
        Quantum numbers for the ket orbital.
    Z_nuc : float
        Bare nuclear charge of the off-center nucleus (B).
    R_AB : float
        Internuclear distance (bohr).
    L_max : int
        Maximum multipole order in the expansion.
    core_type : str, optional
        FrozenCore type for nucleus B: 'Ne', 'Ar', 'Ar3d10', 'Kr', 'Xe'.
        If None, auto-detected from int(Z_nuc). If no frozen core exists
        for this Z (e.g., Z=1..10), reverts to bare-Coulomb behavior.
    n_grid : int
        Radial grid points for screened-correction integration.
    n_legendre : int
        Gauss-Legendre quadrature nodes for the multipole angular projection.
    nuc_parity : int
        +1 if nucleus is in +z direction, -1 if -z.

    Returns
    -------
    float
        Matrix element value in Hartree.
    """
    if m1 != m2:
        return 0.0

    # Auto-detect core type
    if core_type is None:
        core_type = _detect_core_type(Z_nuc)

    # No frozen core: revert to bare Coulomb (exact compatibility)
    if core_type is None:
        from geovac.shibuya_wulfman import compute_cross_center_vne_element
        return compute_cross_center_vne_element(
            Z_orb, n1, l1, m1, n2, l2, m2, Z_nuc, R_AB, L_max,
            n_grid=n_grid, nuc_parity=nuc_parity,
        )

    # Build FrozenCore for B
    fc = FrozenCore(Z=int(round(Z_nuc)), core_type=core_type)
    fc.solve()

    total = 0.0
    for L in range(0, L_max + 1):
        ang = _angular_coefficient(l1, m1, l2, m2, L)
        if abs(ang) < 1e-15:
            continue
        rad_bare, rad_screen = _screened_radial_integral(
            Z_orb, n1, l1, n2, l2, fc, Z_nuc, L, R_AB,
            n_grid=n_grid, n_legendre=n_legendre,
        )
        # V_screened = -Z_nuc/rho + f_screen(rho) gives matrix element:
        # -Z_nuc * rad_bare + rad_screen
        total += nuc_parity ** L * ang * (-Z_nuc * rad_bare + rad_screen)

    return total


def compute_screened_cross_center_vne(
    Z_orb: float,
    states: List[Tuple[int, int, int]],
    Z_nuc: float,
    R_AB: float,
    L_max: int,
    core_type: Optional[str] = None,
    n_grid: int = 4000,
    n_legendre: int = 64,
    nuc_parity: int = 1,
    direction: Optional[Tuple[float, float, float]] = None,
) -> np.ndarray:
    """Full screened cross-center V_ne matrix.

    Drop-in analog of ``shibuya_wulfman.compute_cross_center_vne``. For
    first-row Z_nuc (no frozen core), exactly reproduces the bare result.

    See ``compute_screened_cross_center_vne_element`` for parameter docs.

    Returns
    -------
    np.ndarray
        V_ne matrix of shape (N_orb, N_orb).
    """
    # Auto-detect core type
    if core_type is None:
        core_type = _detect_core_type(Z_nuc)

    # No frozen core: revert to bare path (exact compatibility)
    if core_type is None:
        return compute_cross_center_vne(
            Z_orb, states, Z_nuc, R_AB, L_max,
            n_grid=n_grid, nuc_parity=nuc_parity, direction=direction,
        )

    # --- Direction-based rotation: compute in z-frame, then rotate ---
    if direction is not None:
        vne_z = compute_screened_cross_center_vne(
            Z_orb, states, Z_nuc, R_AB, L_max,
            core_type=core_type,
            n_grid=n_grid, n_legendre=n_legendre,
            nuc_parity=1, direction=None,
        )
        D = _build_block_rotation_matrix(states, direction)
        return D @ vne_z @ D.T

    # --- z-axis frame ---
    fc = FrozenCore(Z=int(round(Z_nuc)), core_type=core_type)
    fc.solve()

    N = len(states)
    vne = np.zeros((N, N))

    # Cache radial integrals keyed by (n1, l1, n2, l2, L)
    unique_nl = sorted(set((n, l) for n, l, m in states))
    radial_cache: Dict[Tuple[int, int, int, int, int], Tuple[float, float]] = {}
    for n1, l1 in unique_nl:
        for n2, l2 in unique_nl:
            for L in range(0, L_max + 1):
                if (l1 + L + l2) % 2 != 0:
                    continue
                if abs(l1 - l2) > L or L > l1 + l2:
                    continue
                key = (n1, l1, n2, l2, L)
                if key not in radial_cache:
                    radial_cache[key] = _screened_radial_integral(
                        Z_orb, n1, l1, n2, l2, fc, Z_nuc, L, R_AB,
                        n_grid=n_grid, n_legendre=n_legendre,
                    )

    # Build matrix
    for i, (n1, l1, m1) in enumerate(states):
        for j, (n2, l2, m2) in enumerate(states):
            if m1 != m2:
                continue

            val = 0.0
            for L in range(0, L_max + 1):
                if (l1 + L + l2) % 2 != 0:
                    continue
                if abs(l1 - l2) > L or L > l1 + l2:
                    continue
                ang = _angular_coefficient(l1, m1, l2, m2, L)
                if abs(ang) < 1e-15:
                    continue
                key = (n1, l1, n2, l2, L)
                rad_bare, rad_screen = radial_cache.get(key, (0.0, 0.0))
                val += nuc_parity ** L * ang * (
                    -Z_nuc * rad_bare + rad_screen
                )

            vne[i, j] = val

    return vne
