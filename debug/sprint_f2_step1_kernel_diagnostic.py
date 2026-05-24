"""Sprint F2 Step 1 — Algebraic cross-V_ne kernel diagnostic.

Diagnostic-before-engineering algebraic test of whether the multipole-
expansion cross-V_ne kernel (current `geovac/shibuya_wulfman.py` path) is
faithful relative to full 3D numerical quadrature at NaH R = R_eq.

Cross-block matrix elements (Na orbital + H orbital) are NOT computed in
the framework (h1 is structurally block-diagonal in composed_qubit).
Therefore the diagnostic targets:
  - within-Na cross-V_ne diagonal element: <Na 3s | V_ne(H) | Na 3s>
  - within-Na cross-V_ne off-diagonal between Na 3s and Na 4s:
    <Na 3s | V_ne(H) | Na 4s>
  - within-H cross-V_ne diagonal element: <H 1s | V_ne(Na) | H 1s>
  - within-H cross-V_ne off-diagonal between H 1s and H 2s:
    <H 1s | V_ne(Na) | H 2s>
plus the bonding/antibonding splitting under h1-only at the truncated
2-state basis {|Na 3s>, |H 1s>} with the within-sub-block contributions
that the framework actually constructs.

For the multipole side we use the current shibuya_wulfman primitives
(both hydrogenic and multi-zeta paths). For the 3D side we use the
standard prolate-spheroidal-style separation: for orbitals both centered
at A and a nucleus at B with linear NaH geometry, the integrand
psi_A(r) psi_A(r) (-Z_B/|r - R_B|) is axially symmetric, so the 3D
integral reduces to a 2D integral in (r, theta) with r in [0, inf) and
theta in [0, pi].
"""

from __future__ import annotations

import json
import math
import time
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.special import legendre, roots_legendre

from geovac.composed_qubit import _radial_wf_grid, _wigner3j
from geovac.multi_zeta_orbitals import (
    MultiZetaOrbital,
    get_physical_valence_orbitals,
)
from geovac.shibuya_wulfman import (
    compute_cross_center_vne_element,
    _hydrogenic_poly_coeffs,
    _radial_split_integral_multizeta,
    _radial_split_integral,
    _angular_coefficient,
    _multizeta_to_poly_components,
    _split_integral_analytical,
)


# ---------------------------------------------------------------------------
# Real spherical harmonics on (theta, phi) for the |r,theta,phi> grid.
# We use l = 0, 1 only for this diagnostic (3s, 4s, 1s, 2s all have l = 0).
# ---------------------------------------------------------------------------


def real_Y_l0(theta: np.ndarray) -> np.ndarray:
    """Real Y_{l=0, m=0}(theta) = 1/sqrt(4 pi)."""
    return np.full_like(theta, 1.0 / math.sqrt(4.0 * math.pi))


def na_3s_multizeta() -> MultiZetaOrbital:
    """Get the Z=11 Na 3s multi-zeta orbital (physical-fit)."""
    orbs = get_physical_valence_orbitals(11)
    for o in orbs:
        if o.n_orbital == 3 and o.l_orbital == 0:
            return o
    raise RuntimeError("Na 3s not in physical-fit registry")


def na_3s_hydrogenic_wf(r: np.ndarray, Z_orb: float = 1.0) -> np.ndarray:
    """Hydrogenic Na 3s at Z_orb=1 (the framework's baseline)."""
    return _radial_wf_grid(Z_orb, 3, 0, r)


def na_orbital_wf(r: np.ndarray, n: int, l: int, Z_orb: float = 1.0,
                  multi_zeta: bool = False) -> np.ndarray:
    """Get Na orbital wavefunction R_nl(r), either hydrogenic or multi-zeta."""
    if multi_zeta and (n, l) == (3, 0):
        return na_3s_multizeta().evaluate(r)
    return _radial_wf_grid(Z_orb, n, l, r)


def h_orbital_wf(r: np.ndarray, n: int, l: int, Z_orb: float = 1.0) -> np.ndarray:
    """Hydrogenic H orbital at Z_orb=1."""
    return _radial_wf_grid(Z_orb, n, l, r)


# ---------------------------------------------------------------------------
# Full 3D quadrature for cross-V_ne matrix element on a SAME-CENTER pair.
# ---------------------------------------------------------------------------


def cross_vne_3d_quadrature_same_center(
    R_bra,                       # callable r -> R_bra(r) radial wavefunction
    l_bra: int, m_bra: int,
    R_ket,                       # callable r -> R_ket(r)
    l_ket: int, m_ket: int,
    R_C: float,                  # distance from orbital center A to off-center nucleus C (along +z)
    Z_C: float,                  # charge of off-center nucleus
    n_radial: int = 200,
    n_angular: int = 200,
    r_max: float = 30.0,
    relative_to_R: bool = True,
) -> float:
    """Compute <psi_{l_bra, m_bra} | -Z_C / |r - R_C * z_hat| | psi_{l_ket, m_ket}>
    by full 3D numerical quadrature.

    Uses Gauss-Legendre quadrature in (r, theta) (m-integration is trivial for
    m_bra = m_ket).

    The integrand:
        I = -Z_C * integral_0^inf dr r^2 R_bra(r) R_ket(r)
              * integral_0^pi sin(theta) dtheta
                  Y_{l_bra, m_bra}(theta, phi) Y_{l_ket, m_ket}*(theta, phi)
                  / sqrt(r^2 - 2 r R_C cos(theta) + R_C^2)
              * 2 pi * delta_{m_bra, m_ket}

    For l_bra = l_ket = 0, m_bra = m_ket = 0:
        Y_{0,0} Y_{0,0}* = 1 / (4 pi)

    so the integral simplifies to
        I = -Z_C / 2 * integral_0^inf dr r^2 R_bra(r) R_ket(r)
              * integral_0^pi dtheta sin(theta) / sqrt(r^2 - 2 r R_C cos(theta) + R_C^2).

    For the inner theta integral, substitution u = cos(theta) yields:
        integral_-1^1 du / sqrt(r^2 - 2 r R_C u + R_C^2)
            = (sqrt(r^2 + 2 r R_C + R_C^2) - sqrt(r^2 - 2 r R_C + R_C^2)) / (r R_C)
            = (|r + R_C| - |r - R_C|) / (r R_C)
            = 2 * min(r, R_C) / (r R_C) [for positive r and R_C]
            = 2 / max(r, R_C)

    Therefore (for s-s case, l_bra = l_ket = m_bra = m_ket = 0):
        I = -Z_C * integral_0^inf dr r^2 R_bra(r) R_ket(r) / max(r, R_C).

    This is the well-known result for s-orbital cross-V_ne. We can use this
    closed form for an even more precise reference value than 3D quadrature.

    For l != 0 cases we use the multipole expansion to high L_max (~20)
    instead of brute-force 3D quadrature; this gives essentially-exact
    reference values since the multipole series for nuclei outside the
    orbital cloud converges geometrically.

    Parameters
    ----------
    R_bra, R_ket : callable
        Radial wavefunctions (numpy ufuncs of r).
    l_bra, m_bra, l_ket, m_ket : int
        Angular quantum numbers.
    R_C : float
        Distance from orbital center to off-center nucleus (bohr, > 0).
    Z_C : float
        Charge of off-center nucleus.
    n_radial : int
        Number of Gauss-Legendre radial points.
    n_angular : int
        Number of Gauss-Legendre angular points (unused for s-s; kept for
        signature compatibility).
    r_max : float
        Cutoff for radial integration (bohr).

    Returns
    -------
    float
        Matrix element <psi_a | -Z_C / |r - R_C z| | psi_b> in Hartree.
    """
    if m_bra != m_ket:
        return 0.0
    # We use the closed-form theta integral for s-s, and the multipole-
    # expansion route for everything else (since multipole is the textbook
    # reference; the question is whether the framework's L_max truncation
    # matches the converged-multipole reference).
    if l_bra == 0 and l_ket == 0 and m_bra == 0 and m_ket == 0:
        # s-s case: use Gauss-Legendre radial integration with closed-form
        # theta integral. Transform r in [0, r_max] -> u in [-1, 1] via
        # r = r_max * (u + 1) / 2 with jacobian dr = r_max/2 du.
        u, w = roots_legendre(n_radial)
        r = r_max * (u + 1.0) / 2.0
        jac = r_max / 2.0
        psi_a = R_bra(r)
        psi_b = R_ket(r)
        # min(r, R_C) / (r R_C) = 1/max(r, R_C)
        inv_max = np.where(r < R_C, 1.0 / R_C, 1.0 / r)
        integrand = r * r * psi_a * psi_b * inv_max
        I = float(np.sum(w * integrand) * jac)
        return -Z_C * I
    else:
        raise NotImplementedError(
            "Only s-s case implemented in this diagnostic (l_bra=l_ket=0)."
        )


# ---------------------------------------------------------------------------
# Reference via converged multipole expansion (high L_max)
# ---------------------------------------------------------------------------


def cross_vne_converged_multipole(
    Z_orb: float,
    n1: int, l1: int, m1: int,
    n2: int, l2: int, m2: int,
    Z_C: float,
    R_AB: float,
    L_max_high: int = 20,
) -> float:
    """Same as compute_cross_center_vne_element but with high L_max for
    convergence reference."""
    return compute_cross_center_vne_element(
        Z_orb, n1, l1, m1, n2, l2, m2, Z_C, R_AB,
        L_max=L_max_high,
    )


def cross_vne_converged_multipole_mz(
    orbital_bra: MultiZetaOrbital,
    orbital_ket: MultiZetaOrbital,
    m_bra: int, m_ket: int,
    Z_C: float,
    R_AB: float,
    L_max_high: int = 20,
) -> float:
    """Multi-zeta cross-V_ne via converged multipole expansion (high L_max).

    Uses _radial_split_integral_multizeta directly; treats l = orbital.l_orbital.
    """
    l1 = orbital_bra.l_orbital
    l2 = orbital_ket.l_orbital
    if m_bra != m_ket:
        return 0.0

    total = 0.0
    for L in range(0, L_max_high + 1):
        ang = _angular_coefficient(l1, m_bra, l2, m_ket, L)
        if abs(ang) < 1e-15:
            continue
        rad = _radial_split_integral_multizeta(
            orbital_bra, orbital_ket, L, R_AB,
        )
        total += ang * rad
    return -Z_C * total


# ---------------------------------------------------------------------------
# Bonding / antibonding diagnostic at the 2-state h1 level
# ---------------------------------------------------------------------------


def bonding_antibonding_h1_diagnostic(
    R_AB: float,
    Z_Na: float = 11.0,
    Z_H: float = 1.0,
    use_multi_zeta_for_Na3s: bool = True,
    L_max: int = 4,
) -> Dict:
    """Compute the bonding/antibonding splitting in the h1-only 2-state basis
    {|Na 3s>, |H 1s>} as it appears in the framework's balanced_coupled.

    The framework's h1 within-sub-block structure for these two orbitals:

      h1[Na 3s, Na 3s] = -Z_Na_val^2 / (2 * n_block^2) + V_ne(H)[Na 3s, Na 3s]
      h1[H 1s, H 1s]   = -Z_H^2 / (2 * 1^2) + V_ne(Na)[H 1s, H 1s]
      h1[Na 3s, H 1s]  = STRUCTURALLY ABSENT in current framework

    Therefore the 2-state h1 IS DIAGONAL in this basis. The bonding and
    antibonding combinations of |Na 3s> +/- |H 1s> have h1 expectation values:

      <bonding | h1 | bonding> = (1/2)(E_Na + E_H)
      <antibonding | h1 | antibonding> = (1/2)(E_Na + E_H)

    So at the pure h1 level there is NO splitting; the bonding/antibonding
    splitting comes entirely from the 2-body ERIs in the FCI.

    The h1 expectation on the BONDING combination is therefore the AVERAGE of
    the two diagonal energies (since no off-diagonal h1 element). The FCI's
    bonding NO occupancy reflects whether the 2-body ERIs can lower the
    bonding combination's energy enough to compete with the separated
    [Na 3s alone | H 1s alone] configuration.

    However, the 2-state diagnostic IS still informative: we can compute
      a) the diagonal h1 energies E_Na = h1[Na 3s, Na 3s] (which include
         the multipole cross-V_ne contribution from H nucleus at distance R)
      b) the same diagonal under FULL 3D quadrature reference for the
         cross-V_ne contribution
      c) the differential between (a) and (b) for both Na 3s and H 1s sites.

    A multipole that is faithful to <1% of the bond-energy scale (~0.075 Ha
    for NaH) cannot be the bottleneck; the F2 hypothesis is then refuted at
    the kernel level.

    Returns
    -------
    dict
        Comprehensive diagnostic table.
    """
    # --- Framework Z values (block_n=1 for both, n_val_offset=2 for Na) ---
    # NaH spec uses nah_spec with Z_orb=1 on both centers (the hydrogenic
    # baseline), block_n=1, max_n=3, n_val_offset=2.
    # The on-site hydrogenic h1 diagonal is -Z_orb^2 / (2 * block_n^2) = -0.5
    # for both Na 3s and H 1s at the un-multizeta level.
    # With multi-zeta on Na 3s, the on-site h1 diagonal is replaced by the
    # Schroedinger eigenvalue of the FrozenCore-screened potential
    # (-0.170 Ha for Na 3s per F1 max_n=3 memo).

    # For the cross-V_ne contribution we have two sub-block computations:
    #   1. Na 3s sub-block sees H nucleus (Z_C = 1) at distance R from Na.
    #   2. H 1s sub-block sees Na nucleus (Z_C = 11) at distance R from H.

    result = {
        'R_AB_bohr': R_AB,
        'Z_Na': Z_Na,
        'Z_H': Z_H,
        'use_multi_zeta_for_Na3s': use_multi_zeta_for_Na3s,
        'L_max_framework': L_max,
        'L_max_high': 20,
    }

    # --- Hydrogenic orbital callables ---
    def R_Na3s_hyd(r): return _radial_wf_grid(1.0, 3, 0, r)
    def R_Na4s_hyd(r): return _radial_wf_grid(1.0, 4, 0, r)
    def R_H1s(r): return _radial_wf_grid(1.0, 1, 0, r)
    def R_H2s(r): return _radial_wf_grid(1.0, 2, 0, r)

    na_3s_mz_orb = na_3s_multizeta()
    def R_Na3s_mz(r): return na_3s_mz_orb.evaluate(r)

    # ---------------------------------------------------------------
    # Matrix element A: <Na 3s | V_ne(H, at R) | Na 3s>
    #   Hydrogenic Na 3s (Z_orb=1) at center A; H nucleus (Z_C=1) at distance R.
    # ---------------------------------------------------------------
    val_mp_hyd_A_diag = compute_cross_center_vne_element(
        Z_orb=1.0,
        n1=3, l1=0, m1=0,
        n2=3, l2=0, m2=0,
        Z_nuc=1.0, R_AB=R_AB, L_max=L_max,
    )
    val_mp_high_hyd_A_diag = cross_vne_converged_multipole(
        Z_orb=1.0,
        n1=3, l1=0, m1=0,
        n2=3, l2=0, m2=0,
        Z_C=1.0, R_AB=R_AB, L_max_high=20,
    )
    val_3d_hyd_A_diag = cross_vne_3d_quadrature_same_center(
        R_bra=R_Na3s_hyd, l_bra=0, m_bra=0,
        R_ket=R_Na3s_hyd, l_ket=0, m_ket=0,
        R_C=R_AB, Z_C=1.0,
        n_radial=400, r_max=40.0,
    )

    # ---------------------------------------------------------------
    # Matrix element B: <Na 3s | V_ne(H, at R) | Na 4s>
    # Off-diagonal within Na sub-block; tests basis-pair coupling.
    # ---------------------------------------------------------------
    val_mp_hyd_B_offdiag = compute_cross_center_vne_element(
        Z_orb=1.0,
        n1=3, l1=0, m1=0,
        n2=4, l2=0, m2=0,
        Z_nuc=1.0, R_AB=R_AB, L_max=L_max,
    )
    val_mp_high_hyd_B_offdiag = cross_vne_converged_multipole(
        Z_orb=1.0,
        n1=3, l1=0, m1=0,
        n2=4, l2=0, m2=0,
        Z_C=1.0, R_AB=R_AB, L_max_high=20,
    )
    val_3d_hyd_B_offdiag = cross_vne_3d_quadrature_same_center(
        R_bra=R_Na3s_hyd, l_bra=0, m_bra=0,
        R_ket=R_Na4s_hyd, l_ket=0, m_ket=0,
        R_C=R_AB, Z_C=1.0,
        n_radial=400, r_max=40.0,
    )

    # ---------------------------------------------------------------
    # Matrix element C: <H 1s | V_ne(Na, at R) | H 1s>
    # H sub-block diagonal; tests the un-screened Z=11 attraction on H 1s.
    # ---------------------------------------------------------------
    val_mp_C_diag = compute_cross_center_vne_element(
        Z_orb=1.0,
        n1=1, l1=0, m1=0,
        n2=1, l2=0, m2=0,
        Z_nuc=11.0, R_AB=R_AB, L_max=L_max,
    )
    val_mp_high_C_diag = cross_vne_converged_multipole(
        Z_orb=1.0,
        n1=1, l1=0, m1=0,
        n2=1, l2=0, m2=0,
        Z_C=11.0, R_AB=R_AB, L_max_high=20,
    )
    val_3d_C_diag = cross_vne_3d_quadrature_same_center(
        R_bra=R_H1s, l_bra=0, m_bra=0,
        R_ket=R_H1s, l_ket=0, m_ket=0,
        R_C=R_AB, Z_C=11.0,
        n_radial=400, r_max=40.0,
    )

    # ---------------------------------------------------------------
    # Matrix element D: <H 1s | V_ne(Na, at R) | H 2s>
    # H off-diagonal; tests basis-pair coupling on H side.
    # ---------------------------------------------------------------
    val_mp_D_offdiag = compute_cross_center_vne_element(
        Z_orb=1.0,
        n1=1, l1=0, m1=0,
        n2=2, l2=0, m2=0,
        Z_nuc=11.0, R_AB=R_AB, L_max=L_max,
    )
    val_mp_high_D_offdiag = cross_vne_converged_multipole(
        Z_orb=1.0,
        n1=1, l1=0, m1=0,
        n2=2, l2=0, m2=0,
        Z_C=11.0, R_AB=R_AB, L_max_high=20,
    )
    val_3d_D_offdiag = cross_vne_3d_quadrature_same_center(
        R_bra=R_H1s, l_bra=0, m_bra=0,
        R_ket=R_H2s, l_ket=0, m_ket=0,
        R_C=R_AB, Z_C=11.0,
        n_radial=400, r_max=40.0,
    )

    # ---------------------------------------------------------------
    # Matrix element E: <Na 3s | V_ne(H, at R) | Na 3s>, multi-zeta on Na 3s
    # ---------------------------------------------------------------
    val_mp_high_mz_E_diag = cross_vne_converged_multipole_mz(
        orbital_bra=na_3s_mz_orb,
        orbital_ket=na_3s_mz_orb,
        m_bra=0, m_ket=0,
        Z_C=1.0, R_AB=R_AB, L_max_high=20,
    )
    val_3d_mz_E_diag = cross_vne_3d_quadrature_same_center(
        R_bra=R_Na3s_mz, l_bra=0, m_bra=0,
        R_ket=R_Na3s_mz, l_ket=0, m_ket=0,
        R_C=R_AB, Z_C=1.0,
        n_radial=400, r_max=40.0,
    )

    result['matrix_elements'] = {
        'A_NaA_3s__VneH__NaA_3s_diag_hydrogenic': {
            'description': '<Na 3s_hyd | V_ne(H, R) | Na 3s_hyd>',
            'mp_L_max_framework_Ha': val_mp_hyd_A_diag,
            'mp_L_max_high_Ha': val_mp_high_hyd_A_diag,
            '3d_quadrature_Ha': val_3d_hyd_A_diag,
            'mp_vs_3d_diff_Ha': val_mp_high_hyd_A_diag - val_3d_hyd_A_diag,
            'mp_vs_3d_relerr': abs(val_mp_high_hyd_A_diag - val_3d_hyd_A_diag) / max(abs(val_3d_hyd_A_diag), 1e-30),
        },
        'B_NaA_3s__VneH__NaA_4s_offdiag_hydrogenic': {
            'description': '<Na 3s_hyd | V_ne(H, R) | Na 4s_hyd>',
            'mp_L_max_framework_Ha': val_mp_hyd_B_offdiag,
            'mp_L_max_high_Ha': val_mp_high_hyd_B_offdiag,
            '3d_quadrature_Ha': val_3d_hyd_B_offdiag,
            'mp_vs_3d_diff_Ha': val_mp_high_hyd_B_offdiag - val_3d_hyd_B_offdiag,
            'mp_vs_3d_relerr': abs(val_mp_high_hyd_B_offdiag - val_3d_hyd_B_offdiag) / max(abs(val_3d_hyd_B_offdiag), 1e-30),
        },
        'C_HB_1s__VneNa__HB_1s_diag': {
            'description': '<H 1s | V_ne(Na, R) | H 1s>',
            'mp_L_max_framework_Ha': val_mp_C_diag,
            'mp_L_max_high_Ha': val_mp_high_C_diag,
            '3d_quadrature_Ha': val_3d_C_diag,
            'mp_vs_3d_diff_Ha': val_mp_high_C_diag - val_3d_C_diag,
            'mp_vs_3d_relerr': abs(val_mp_high_C_diag - val_3d_C_diag) / max(abs(val_3d_C_diag), 1e-30),
        },
        'D_HB_1s__VneNa__HB_2s_offdiag': {
            'description': '<H 1s | V_ne(Na, R) | H 2s>',
            'mp_L_max_framework_Ha': val_mp_D_offdiag,
            'mp_L_max_high_Ha': val_mp_high_D_offdiag,
            '3d_quadrature_Ha': val_3d_D_offdiag,
            'mp_vs_3d_diff_Ha': val_mp_high_D_offdiag - val_3d_D_offdiag,
            'mp_vs_3d_relerr': abs(val_mp_high_D_offdiag - val_3d_D_offdiag) / max(abs(val_3d_D_offdiag), 1e-30),
        },
        'E_NaA_3s_mz__VneH__NaA_3s_mz_diag_multi_zeta': {
            'description': '<Na 3s_mz | V_ne(H, R) | Na 3s_mz>',
            'mp_L_max_high_Ha': val_mp_high_mz_E_diag,
            '3d_quadrature_Ha': val_3d_mz_E_diag,
            'mp_vs_3d_diff_Ha': val_mp_high_mz_E_diag - val_3d_mz_E_diag,
            'mp_vs_3d_relerr': abs(val_mp_high_mz_E_diag - val_3d_mz_E_diag) / max(abs(val_3d_mz_E_diag), 1e-30),
        },
    }

    # ---------------------------------------------------------------
    # Bonding/antibonding splitting on the 2-state {|Na 3s>, |H 1s>} basis
    # under h1 only.
    # ---------------------------------------------------------------
    # Framework on-site h1 diagonal (with multi-zeta on Na 3s when requested):
    # Na 3s without multi-zeta: -Z_orb^2 / (2 * block_n^2) = -1/2 (Z_orb=1, block_n=1)
    # Na 3s with multi-zeta: ~ -0.170 Ha (the FrozenCore-screened eigenvalue)
    # H 1s: -Z_orb^2 / (2 * 1^2) = -1/2 Ha
    h1_diag_Na3s_onsite_hyd = -0.5
    h1_diag_Na3s_onsite_mz = -0.170  # approximation from F1 memo
    h1_diag_H1s_onsite = -0.5

    # Add cross-V_ne contributions from the existing framework path
    # (multipole at L_max=4)
    h1_Na3s_framework_hyd = h1_diag_Na3s_onsite_hyd + val_mp_hyd_A_diag
    h1_H1s_framework = h1_diag_H1s_onsite + val_mp_C_diag

    # Reference under 3D quadrature
    h1_Na3s_3d_hyd = h1_diag_Na3s_onsite_hyd + val_3d_hyd_A_diag
    h1_H1s_3d = h1_diag_H1s_onsite + val_3d_C_diag

    # Multi-zeta variants for Na 3s (cross-V_ne contribution same modulo
    # multi-zeta-shape difference)
    h1_Na3s_framework_mz = h1_diag_Na3s_onsite_mz + val_mp_high_mz_E_diag
    h1_Na3s_3d_mz = h1_diag_Na3s_onsite_mz + val_3d_mz_E_diag

    # In a 2-state DIAGONAL h1, the bonding combination |B> = (|Na 3s> + |H 1s>)/sqrt(2)
    # has <B|h1|B> = (1/2)(h1_Na3s + h1_H1s).
    # The antibonding combination |A> = (|Na 3s> - |H 1s>)/sqrt(2) has the SAME
    # expectation value: <A|h1|A> = (1/2)(h1_Na3s + h1_H1s).
    # So there is NO bonding-antibonding splitting from h1 alone (since h1 is
    # diagonal in the 2-state space).
    # The "bonding < antibonding" energetic ordering must come from 2-body
    # ERIs, which is consistent with the F1 finding that the FCI constructs
    # the bonding NO via 2-body coupling.
    bonding_h1_framework_hyd = 0.5 * (h1_Na3s_framework_hyd + h1_H1s_framework)
    antibonding_h1_framework_hyd = bonding_h1_framework_hyd
    bonding_h1_3d_hyd = 0.5 * (h1_Na3s_3d_hyd + h1_H1s_3d)
    antibonding_h1_3d_hyd = bonding_h1_3d_hyd

    result['bonding_antibonding_h1_only'] = {
        'hydrogenic_path': {
            'h1_Na3s_diag_framework_Ha': h1_Na3s_framework_hyd,
            'h1_H1s_diag_framework_Ha': h1_H1s_framework,
            'h1_Na3s_diag_3d_Ha': h1_Na3s_3d_hyd,
            'h1_H1s_diag_3d_Ha': h1_H1s_3d,
            'bonding_h1_framework_Ha': bonding_h1_framework_hyd,
            'antibonding_h1_framework_Ha': antibonding_h1_framework_hyd,
            'bonding_minus_antibonding_framework_Ha': 0.0,  # diagonal h1 in 2-state
            'bonding_h1_3d_Ha': bonding_h1_3d_hyd,
            'antibonding_h1_3d_Ha': antibonding_h1_3d_hyd,
            'bonding_minus_antibonding_3d_Ha': 0.0,
        },
        'multi_zeta_path': {
            'h1_Na3s_mz_diag_framework_Ha': h1_Na3s_framework_mz,
            'h1_Na3s_mz_diag_3d_Ha': h1_Na3s_3d_mz,
            'bonding_h1_framework_Ha': 0.5 * (h1_Na3s_framework_mz + h1_H1s_framework),
            'bonding_h1_3d_Ha': 0.5 * (h1_Na3s_3d_mz + h1_H1s_3d),
        },
        'interpretation': (
            "At the 2-state h1-only level in the {|Na 3s>, |H 1s>} basis, "
            "the framework has NO off-diagonal h1 cross-block element "
            "<Na 3s | V_ne | H 1s> (h1 is structurally block-diagonal in "
            "composed_qubit). Therefore h1 alone gives NO bonding-antibonding "
            "splitting in this basis -- the bonding orbital that F1 max_n=3 "
            "found is constructed by the 2-body ERIs in the FCI, not by "
            "h1 cross-coupling. The 'kernel-shape' question for h1 is "
            "therefore narrower: does multipole vs 3D-quadrature give the "
            "correct DIAGONAL cross-V_ne contributions (which set the "
            "absolute energies on each side), and does the on-site energy "
            "of the bonding-constructable combination differ enough from "
            "the separated configuration to be the bottleneck?"
        ),
    }

    return result


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------


def main():
    t0 = time.time()
    R_eq_nah = 3.566  # bohr, experimental NaH equilibrium

    print("=" * 78)
    print("Sprint F2 Step 1 — Algebraic cross-V_ne kernel diagnostic")
    print("=" * 78)
    print(f"NaH R_eq = {R_eq_nah} bohr (experimental)")
    print(f"Framework L_max = 4 (NaH max_n=3 default per F1)")
    print(f"Reference L_max = 20 (converged multipole)")
    print()

    diag = bonding_antibonding_h1_diagnostic(
        R_AB=R_eq_nah, L_max=4,
    )

    # Print matrix elements
    print("Matrix element comparison: multipole vs 3D quadrature")
    print("-" * 78)
    for label, me in diag['matrix_elements'].items():
        print(f"\n[{label}]")
        print(f"  {me['description']}")
        if 'mp_L_max_framework_Ha' in me:
            print(f"    mp L_max=4   : {me['mp_L_max_framework_Ha']:+.10f} Ha")
        print(f"    mp L_max=20  : {me['mp_L_max_high_Ha']:+.10f} Ha")
        print(f"    3D quadrature: {me['3d_quadrature_Ha']:+.10f} Ha")
        print(f"    mp_vs_3d_diff: {me['mp_vs_3d_diff_Ha']:+.4e} Ha")
        print(f"    mp_vs_3d_relerr: {me['mp_vs_3d_relerr']:.4e}")

    print()
    print("Bonding/antibonding h1-only analysis (2-state {|Na 3s>, |H 1s>})")
    print("-" * 78)
    bb = diag['bonding_antibonding_h1_only']
    print("Hydrogenic path:")
    for k, v in bb['hydrogenic_path'].items():
        if isinstance(v, (int, float)):
            print(f"  {k}: {v:+.6f}")
    print("Multi-zeta path:")
    for k, v in bb['multi_zeta_path'].items():
        if isinstance(v, (int, float)):
            print(f"  {k}: {v:+.6f}")
    print()
    print(f"Interpretation: {bb['interpretation']}")
    print()

    # Save JSON
    output_path = "debug/data/sprint_f2_step1_kernel_diagnostic.json"
    out = {
        'sprint': 'F2 Step 1 — algebraic kernel diagnostic',
        'date': '2026-05-23',
        'duration_s': time.time() - t0,
        **diag,
    }
    with open(output_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Wrote: {output_path}")
    print(f"\nElapsed: {time.time() - t0:.2f} s")


if __name__ == "__main__":
    main()
