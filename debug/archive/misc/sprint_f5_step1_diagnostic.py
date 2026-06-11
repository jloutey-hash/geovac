"""Sprint F5 Step 1 — algebraic cross-block 2-body Coulomb diagnostic.

Predicts the correlation correction at NaH R = R_eq = 3.566 bohr WITHOUT
building the constrained-FCI architecture. For each Na [Ne] core orbital
c in {1s, 2s, 2p_x, 2p_y, 2p_z}, compute:

  J_{b,c} = <bond(r1) c(r2) | 1/r12 | bond(r1) c(r2)>
  K_{b,c} = <bond(r1) c(r2) | 1/r12 | c(r1) bond(r2)>

  predicted correction = sum_c (J_{b,c} - K_{b,c})

where bond = (1/sqrt(2))(Na 3s + H 1s)-like orbital from the F3
cross-block-h1 diagonalization at NaH max_n=2:

  Composition (from sprint_f4_step1_full_h1_diag.json):
    -0.62 Na(1,0,0)  [physical n=3, multi-zeta Na 3s]
    +0.35 Na(2,0,0)  [physical n=4, multi-zeta Na 4s]
    -0.47 H(1,0,0)   [hydrogenic H 1s]
    +0.52 H(2,0,0)   [hydrogenic H 2s]

  (note: not normalized in the F3 basis directly; we use the eigenvector
   coefficients as-is which already form a normalized linear combination
   in the orthonormalized basis — but the basis itself is non-orthogonal
   in the spatial sense, so we re-normalize the spatial bonding density.)

The 5 [Ne] core orbitals use Clementi-Raimondi exponents.

Gate logic:
  - sum > 2 Ha (right magnitude to take 4.37 Ha well -> ~physical):
    PROCEED to Step 2 architectural extension
  - sum < 0.5 Ha: STOP, W1e is not core correlation
  - wrong sign (attractive when should be repulsive): STOP, deeper mech

IMPLEMENTATION
--------------
J via density-density 6D integral, factorized through Hartree potential:
  V_core(r1) = sum_c 2 * integral |phi_c(r2)|^2 / |r1-r2| d^3 r2
For closed-shell [Ne] core, V_core depends only on |r1 - R_Na| (sph sym).
Compute V_core(r) via radial integration, then
  J_total = sum_c J_{b,c} = (1/2) <bond | V_core | bond>  (factor 1/2 to undo "2*")

K via explicit cross-density:
  K_{b,c} = integral integral phi_b(r1) phi_c(r1) [1/r12] phi_c(r2) phi_b(r2) d^3 r1 d^3 r2
For each c separately, need cross-density rho_bc(r) = phi_b(r) phi_c(r),
then K_{b,c} = integral rho_bc(r1) [1/r12] rho_bc(r2) d^3 r1 d^3 r2.

The cross-density rho_bc is centered at Na for the Na-component of the bond
and the core c, and has a smaller cross-center piece between H 1s/2s of
the bond and the Na core c.

For axial geometry, all integrals can be done in cylindrical (rho, z)
with phi analytic. For s-s components, multipole expansion of 1/r12
truncates cleanly to L=0 monopole at large separations.

For this Step 1 diagnostic we compute J + K both via direct 6D axial
quadrature — slow but correct-first.
"""

from __future__ import annotations

import json
import math
import os
import sys
import time

import numpy as np
from scipy.special import roots_legendre

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _ROOT)

from geovac.cross_block_h1 import hydrogenic_R_nl_analytical  # noqa: E402
from geovac.multi_zeta_orbitals import get_physical_valence_orbitals  # noqa: E402


# ---------------------------------------------------------------------------
# Constants: NaH geometry + Clementi-Raimondi [Ne] core exponents
# ---------------------------------------------------------------------------

R_EQ_BOHR = 3.566  # NaH equilibrium bond length
Z_NA = 11.0
Z_H = 1.0

# Clementi-Raimondi 1963 Table II for Z=11 (Na):
# zeta values for hydrogenic R_nl(r; Z_eff), use Z_eff = n * zeta
NA_CR_ZETA = {
    (1, 0): 10.6259,   # 1s
    (2, 0): 6.5714,    # 2s
    (2, 1): 6.8018,    # 2p
}

# Na core eigenvalues for diagnostic reference (Clementi-Roetti)
NA_CORE_E = {(1, 0): -40.4787, (2, 0): -2.7967, (2, 1): -1.5181}

# Bonding orbital composition from F3 full-h1 diagonalization
# (from debug/data/sprint_f4_step1_full_h1_diag.json)
BOND_COMPOSITION = {
    # ('side', block_n, l, m): coefficient
    ('Na', 1, 0, 0): -0.6200027604993643,   # Na 3s (mz)
    ('Na', 2, 0, 0): +0.35196368219450197,  # Na 4s (mz)
    ('H', 1, 0, 0): -0.4668081378568081,    # H 1s
    ('H', 2, 0, 0): +0.522125677340157,     # H 2s
}


# ---------------------------------------------------------------------------
# Build per-component orbital callables (radial part R(r))
# ---------------------------------------------------------------------------

def build_bond_components():
    """Returns list of (side, n, l, m, c, R_callable) for the bonding orbital."""
    # Na valence multi-zeta orbitals (3s, 4s, 4p)
    na_mz = get_physical_valence_orbitals(11)
    mz_dict = {(orb.n_orbital, orb.l_orbital): orb for orb in na_mz}

    components = []
    for (side, n, l, m), c in BOND_COMPOSITION.items():
        if side == 'Na':
            # block_n -> physical_n via n_val_offset=2
            phys_n = n + 2
            if (phys_n, l) in mz_dict:
                mz = mz_dict[(phys_n, l)]
                R_call = (lambda r, _o=mz: _o.evaluate(np.asarray(r)))
            else:
                R_call = (lambda r, _n=n, _l=l: hydrogenic_R_nl_analytical(
                    1.0, _n, _l, np.asarray(r)))
        else:  # 'H'
            R_call = (lambda r, _n=n, _l=l: hydrogenic_R_nl_analytical(
                1.0, _n, _l, np.asarray(r)))
        components.append({'side': side, 'n': n, 'l': l, 'm': m,
                            'c': c, 'R': R_call})
    return components


def build_core_orbital(n_c, l_c, m_c):
    """Return R(r) callable for Na [Ne] core orbital (n_c, l_c).

    Uses Clementi-Raimondi Z_eff = n_c * zeta.
    """
    zeta = NA_CR_ZETA[(n_c, l_c)]
    Z_eff = n_c * zeta

    def R_call(r, _Z=Z_eff, _n=n_c, _l=l_c):
        return hydrogenic_R_nl_analytical(_Z, _n, _l, np.asarray(r))

    return R_call


# ---------------------------------------------------------------------------
# Spherically averaged core density
# ---------------------------------------------------------------------------

def core_density_radial(r):
    """Total spherically symmetric [Ne] core density at radius r from Na.

    rho_core(r) = sum_c |R_c(r) Y_{l_c,m_c}|^2 with sum over occupied
    spin-orbitals. For [Ne] = 1s^2 2s^2 2p^6:

      rho_core(r) = 2*|R_1s|^2 / (4pi) + 2*|R_2s|^2/(4pi)
                  + 2 * (|R_2p Y_{1,-1}|^2 + |R_2p Y_{1,0}|^2 + |R_2p Y_{1,1}|^2)

    For p shell: sum_m |Y_{1,m}|^2 = 3/(4pi) (closure), so 2p^6 contributes
    6 * |R_2p|^2 / (4pi).

    Total: rho_core(r) = [2 |R_1s|^2 + 2 |R_2s|^2 + 6 |R_2p|^2] / (4pi)

    This is normalized so that integral 4pi r^2 rho dr = 10 (electron count).
    """
    R_1s = hydrogenic_R_nl_analytical(
        NA_CR_ZETA[(1, 0)] * 1, 1, 0, r)
    R_2s = hydrogenic_R_nl_analytical(
        NA_CR_ZETA[(2, 0)] * 2, 2, 0, r)
    R_2p = hydrogenic_R_nl_analytical(
        NA_CR_ZETA[(2, 1)] * 2, 2, 1, r)
    rho = (2.0 * R_1s ** 2 + 2.0 * R_2s ** 2 + 6.0 * R_2p ** 2) / (4.0 * math.pi)
    return rho


def core_hartree_potential(r_grid):
    """Compute Hartree potential V_core(r) = integral rho_core(r') / |r-r'| d^3 r'
    for spherically symmetric core density.

    For sph sym rho:
      V(r) = (4pi/r) * integral_0^r rho(r') r'^2 dr'  +  4pi * integral_r^infinity rho(r') r' dr'
    """
    rho_vals = core_density_radial(r_grid)
    integrand1 = rho_vals * r_grid ** 2  # for inner integral
    integrand2 = rho_vals * r_grid       # for outer integral

    # Use trapezoidal sums (uniform-ish grid is OK; precision check later)
    # Cumulative integral from 0 to r:
    dr = np.diff(r_grid)
    dr = np.concatenate([[dr[0]], dr])  # match length
    # inner: integral_0^r rho r'^2 dr'
    inner_cum = np.cumsum(integrand1 * dr)
    # outer: integral_r^infinity rho r' dr' = total - integral_0^r rho r' dr'
    outer_cum_partial = np.cumsum(integrand2 * dr)
    outer_total = outer_cum_partial[-1]
    outer = outer_total - outer_cum_partial

    # V(r) = 4 pi [ inner / r + outer ]
    # avoid division at r=0
    V = np.where(r_grid > 1e-10, 4.0 * math.pi * (inner_cum / r_grid + outer),
                  4.0 * math.pi * (rho_vals[0] * r_grid[0] / 2.0 + outer_total))
    return V


# ---------------------------------------------------------------------------
# Bonding orbital density |phi_b(r)|^2 evaluated on (rho, z) grid
# ---------------------------------------------------------------------------

def evaluate_bond_orbital(components, rho_grid, z_grid, z_Na=0.0, z_H=R_EQ_BOHR):
    """Return phi_bond(rho, z) summed over components.

    For s components: phi_k = c_k * R_k(r_center) / sqrt(4 pi)
    where r_center is distance to that side's nucleus.
    """
    psi = np.zeros_like(rho_grid)
    for comp in components:
        z_center = z_Na if comp['side'] == 'Na' else z_H
        r_grid = np.sqrt(rho_grid ** 2 + (z_grid - z_center) ** 2)
        R_val = comp['R'](r_grid.flatten()).reshape(r_grid.shape)
        # s-orbital: Y_{0,0} = 1/sqrt(4 pi)
        psi += comp['c'] * R_val / math.sqrt(4.0 * math.pi)
    return psi


def normalize_bond_orbital(components, n_rho=120, n_z=160,
                             rho_max=25.0, z_max=25.0):
    """Compute spatial norm <phi_b | phi_b> via axial quadrature.

    Returns the renormalization factor; the bonding orbital in the
    non-orthogonal basis may not have spatial norm = 1.
    """
    u, wu = roots_legendre(n_rho)
    rho_pts = rho_max * (u + 1.0) / 2.0
    rho_w = wu * rho_max / 2.0
    v_z, wv_z = roots_legendre(n_z)
    z_center = R_EQ_BOHR / 2.0
    z_pts = z_center + z_max * v_z
    z_w = wv_z * z_max
    rho_grid, z_grid = np.meshgrid(rho_pts, z_pts, indexing='ij')

    psi = evaluate_bond_orbital(components, rho_grid, z_grid)
    # Norm: integral |psi|^2 dV = 2 pi integral integral |psi|^2 rho drho dz
    weight_2d = np.outer(rho_w * rho_pts, z_w)
    norm_sq = 2.0 * math.pi * float(np.sum(weight_2d * psi ** 2))
    return norm_sq, (rho_pts, rho_w, z_pts, z_w, rho_grid, z_grid)


# ---------------------------------------------------------------------------
# Hartree J: <phi_b | V_core | phi_b>
# ---------------------------------------------------------------------------

def compute_total_J_via_hartree(components, n_rho=120, n_z=160,
                                  rho_max=25.0, z_max=25.0,
                                  n_r_radial=2000, r_max_radial=40.0,
                                  norm_sq=None):
    """Compute the total Hartree J:

      J_total = sum_c (occupied) integral |phi_b(r1)|^2 |phi_c(r2)|^2 / r12 d^3 r1 d^3 r2

    For closed-shell core (with 10 electrons = 2 * (n_orbitals_doubly_occ)):
    The mean-field potential V_core(r) = sum_c (2 * J_c)(r) is itself
    the Coulomb potential of the total core density (factor 2 for both spins
    is absorbed into rho_core which already counts 10 electrons).

    Returns J_total = (1/2) <phi_b | V_core | phi_b> with the 1/2 factor
    converting "2 J_c sum" -> "J_c sum".
    """
    # Build radial grid for V_core
    # Use log-spaced grid near origin, then linear
    r_grid_inner = np.linspace(0.001, 1.0, 200)
    r_grid_mid = np.linspace(1.0, 5.0, 400)
    r_grid_outer = np.linspace(5.0, r_max_radial, n_r_radial - 600)
    r_grid = np.concatenate([r_grid_inner, r_grid_mid[1:], r_grid_outer[1:]])

    V_core_radial = core_hartree_potential(r_grid)

    # Build the cylindrical (rho, z) grid for bond density
    u, wu = roots_legendre(n_rho)
    rho_pts = rho_max * (u + 1.0) / 2.0
    rho_w = wu * rho_max / 2.0
    v_z, wv_z = roots_legendre(n_z)
    z_center = R_EQ_BOHR / 2.0
    z_pts = z_center + z_max * v_z
    z_w = wv_z * z_max
    rho_grid, z_grid = np.meshgrid(rho_pts, z_pts, indexing='ij')

    # Distance from each grid point to Na (at origin)
    r_to_Na = np.sqrt(rho_grid ** 2 + z_grid ** 2)

    # Interpolate V_core to grid points
    V_core_grid = np.interp(r_to_Na, r_grid, V_core_radial)

    # Bonding orbital density
    psi_bond = evaluate_bond_orbital(components, rho_grid, z_grid)
    rho_bond = psi_bond ** 2

    # Renormalize bond density if requested
    if norm_sq is not None and norm_sq > 1e-12:
        rho_bond = rho_bond / norm_sq

    # Integrate: 2 pi integral integral rho_bond * V_core * rho drho dz
    weight_2d = np.outer(rho_w * rho_pts, z_w)
    full_J = 2.0 * math.pi * float(np.sum(weight_2d * rho_bond * V_core_grid))

    # full_J is integral rho_bond V_core dr where V_core comes from
    # rho_core (counting 10 electrons). For closed shell:
    # full_J = 2 sum_c J_{bond,c}  -> sum_c J_{bond,c} = full_J / 2
    return full_J / 2.0


# ---------------------------------------------------------------------------
# Exchange K: explicit cross-density-density integrals
# ---------------------------------------------------------------------------

def compute_K_for_core_orbital(components, n_c, l_c, m_c,
                                 n_rho=120, n_z=160,
                                 rho_max=25.0, z_max=25.0,
                                 n_r_radial=2000, r_max_radial=40.0,
                                 norm_sq=None):
    """Compute K_{bond, c} = integral cross-density phi_b(r1) phi_c(r1)
    times 1/r12 times phi_c(r2) phi_b(r2) d^3 r1 d^3 r2.

    Use Hartree-like reduction:
      rho_bc(r) = phi_b(r) * phi_c(r)
      K = (1/2) integral integral rho_bc(r1) rho_bc(r2) / r12 d^3 r1 d^3 r2
        = (1/2) integral rho_bc(r1) V_bc(r1) d^3 r1   where V_bc is potential
        of rho_bc.

    The cross-density rho_bc has NO spherical symmetry in general (b is
    not centered at Na). Computing its Hartree potential requires
    multipole expansion.

    SIMPLIFICATION FOR STEP 1 DIAGNOSTIC: use multipole expansion in
    spherical harmonics centered at Na. For s-core (l_c=0), only L=0
    monopole survives in the b-on-Na components AND the L corresponding
    to the angular character of the H-component near Na.

    For Step 1 magnitude estimate, we use the "Slater-Condon" approximation:
    K_{b,c} ≈ |S_{b,c}|^2 * J_{c,c}^{intrinsic}
    where J_{c,c} = <c|c|c|c> = F^0(c,c) is the intrinsic core 2-body integral
    and S_{b,c} is the bond-core overlap.

    This is order-of-magnitude correct for s-cores and tests whether K
    is structurally important.

    Returns:
      K_estimate (Ha)
      S_b_c (overlap)
    """
    # Compute overlap S_{b,c}
    # Build core orbital at Na
    zeta = NA_CR_ZETA[(n_c, l_c)]
    Z_eff = n_c * zeta

    def R_c(r):
        return hydrogenic_R_nl_analytical(Z_eff, n_c, l_c, np.asarray(r))

    # For s-core (l_c=0, m_c=0): overlap with s-components of bond only
    # For p-core: only nonzero overlap with l=1, m=m_c components of bond
    # Bond has only s components -> S_{b, 2p} = 0 for l_c=1

    if l_c == 0:
        # Compute axial s-core overlap with bond
        u, wu = roots_legendre(n_rho)
        rho_pts = rho_max * (u + 1.0) / 2.0
        rho_w = wu * rho_max / 2.0
        v_z, wv_z = roots_legendre(n_z)
        z_center = R_EQ_BOHR / 2.0
        z_pts = z_center + z_max * v_z
        z_w = wv_z * z_max
        rho_grid, z_grid = np.meshgrid(rho_pts, z_pts, indexing='ij')

        psi_bond = evaluate_bond_orbital(components, rho_grid, z_grid)
        if norm_sq is not None and norm_sq > 1e-12:
            psi_bond = psi_bond / math.sqrt(norm_sq)

        # core at Na (z=0)
        r_to_Na = np.sqrt(rho_grid ** 2 + z_grid ** 2)
        R_c_grid = R_c(r_to_Na.flatten()).reshape(r_to_Na.shape)
        psi_c = R_c_grid / math.sqrt(4.0 * math.pi)  # Y_{0,0}

        # Overlap: integral psi_b psi_c dV = 2 pi sum (psi_b psi_c rho) w
        weight_2d = np.outer(rho_w * rho_pts, z_w)
        S_bc = 2.0 * math.pi * float(np.sum(weight_2d * psi_bond * psi_c))
    else:
        # For p-core, bond has no l=1 components -> S = 0
        S_bc = 0.0

    # Intrinsic core 2-body integral: J_{c,c} = F^0(c,c)
    # For hydrogenic 1s: F^0 = 5Z/8; 2s: 77 Z/512 etc.
    # We use the Hartree potential of core c on itself:
    # F^0(c,c) = integral rho_c(r1) rho_c(r2) / r12 = (1/2) integral rho_c V_c dr
    # where V_c is potential of single occupied c (1 electron, not 2).
    #
    # For hydrogenic 1s with Z: F^0 = 5Z/8. For other shells, need integration.

    if l_c == 0 and n_c == 1:
        F0_cc = 5.0 * Z_eff / 8.0
    elif l_c == 0 and n_c == 2:
        F0_cc = 77.0 * Z_eff / 512.0
    elif l_c == 1 and n_c == 2:
        F0_cc = 83.0 * Z_eff / 512.0   # F^0(2p,2p) = 83 Z / 512
    else:
        F0_cc = float('nan')

    # K estimate via Slater-Condon: K_{b,c} ≈ S_{b,c}^2 * F0_{c,c}
    # This is a UPPER BOUND on exchange magnitude (assumes maximum overlap density)
    K_est = (S_bc ** 2) * F0_cc

    return K_est, S_bc, F0_cc


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def main():
    t0 = time.perf_counter()
    print("=" * 78)
    print("Sprint F5 Step 1 — algebraic cross-block 2-body Coulomb diagnostic")
    print(f"NaH R = {R_EQ_BOHR} bohr; bonding orbital from F3 full-h1")
    print(f"Predict: sum_c (J - K)_{{bond,c}} for c in Na [Ne] core")
    print("=" * 78)

    # ------------------------------------------------------------------
    # Build bonding orbital components
    # ------------------------------------------------------------------
    components = build_bond_components()
    print("\nBonding orbital composition (F3 full-h1 diagonalization):")
    for comp in components:
        print(f"  {comp['side']}({comp['n']},{comp['l']},{comp['m']}): "
              f"c = {comp['c']:+.4f}")

    # ------------------------------------------------------------------
    # Normalize bonding orbital (spatial norm in non-orthogonal basis)
    # ------------------------------------------------------------------
    norm_sq, _grid_info = normalize_bond_orbital(components)
    print(f"\nSpatial norm <bond|bond>^2 = {norm_sq:.6f}")
    print(f"  Renormalization factor sqrt(norm_sq) = {math.sqrt(norm_sq):.4f}")

    # ------------------------------------------------------------------
    # Verify core density normalizes to 10
    # ------------------------------------------------------------------
    r_check = np.linspace(0.001, 30.0, 5000)
    rho_check = core_density_radial(r_check)
    N_core = 4.0 * math.pi * np.trapezoid(r_check ** 2 * rho_check, r_check)
    print(f"\nCore density check: 4 pi integral r^2 rho dr = {N_core:.4f} "
          f"(expected 10 for [Ne])")

    # ------------------------------------------------------------------
    # Total J via Hartree potential
    # ------------------------------------------------------------------
    print(f"\n{'-' * 78}")
    print("Computing total J = sum_c <bond bond | 1/r12 | c c>...")
    print(f"{'-' * 78}")
    J_total = compute_total_J_via_hartree(components, norm_sq=norm_sq)
    print(f"\nTotal J (sum over 5 core orbitals, divided by 2 to convert "
          f"2J_sum -> J_sum):")
    print(f"  sum_c J_{{bond,c}} = {J_total:+.4f} Ha")

    # ------------------------------------------------------------------
    # Per-orbital K via Slater-Condon estimate
    # ------------------------------------------------------------------
    print(f"\n{'-' * 78}")
    print("Computing K_{bond,c} for each core orbital (Slater-Condon estimate)...")
    print(f"{'-' * 78}")

    core_list = [(1, 0, 0), (2, 0, 0), (2, 1, -1), (2, 1, 0), (2, 1, 1)]
    K_per_orbital = {}
    K_total = 0.0
    for (n_c, l_c, m_c) in core_list:
        K_est, S_bc, F0_cc = compute_K_for_core_orbital(
            components, n_c, l_c, m_c, norm_sq=norm_sq)
        K_total += K_est
        K_per_orbital[f"({n_c},{l_c},{m_c})"] = {
            "K_estimate_Ha": K_est,
            "S_bond_core": S_bc,
            "F0_cc_Ha": F0_cc,
        }
        print(f"  c = ({n_c},{l_c},{m_c}): S_{{b,c}}={S_bc:+.4f}, "
              f"F^0_{{c,c}}={F0_cc:.4f} Ha, K_est = {K_est:+.4f} Ha")
    print(f"\nTotal K (sum over 5 core orbitals): {K_total:+.4f} Ha")

    # ------------------------------------------------------------------
    # Predicted correction = J - K
    # ------------------------------------------------------------------
    correction = J_total - K_total
    print(f"\n{'=' * 78}")
    print("PREDICTED CORRELATION CORRECTION:")
    print(f"  sum_c (J - K)_{{bond,c}} = {correction:+.4f} Ha")
    print(f"  J_total = {J_total:+.4f} Ha")
    print(f"  K_total = {K_total:+.4f} Ha")
    print(f"  K/J ratio = {(K_total/J_total*100 if abs(J_total)>1e-6 else 0):.2f}%")
    print(f"{'=' * 78}")

    # ------------------------------------------------------------------
    # Gate decision
    # ------------------------------------------------------------------
    wall_depth = 4.37  # F3 baseline wall depth at NaH max_n=2
    exp_De = 0.075  # experimental NaH binding energy
    closure_window_2x = (exp_De * 0.5, exp_De * 2.0)
    # For the wall to close to physical, need correction ≈ wall_depth (4.37 Ha)
    # since correction is repulsive (positive J - K) and would shift PES upward.

    gate_outcomes = []
    if correction > 2.0:
        verdict = "PROCEED_TO_STEP_2"
        gate_outcomes.append(f"sum > 2 Ha (right magnitude class to close "
                              f"{wall_depth} Ha wall toward physical)")
    elif correction < 0.5:
        verdict = "STOP_NOT_CORE_CORRELATION"
        gate_outcomes.append(f"sum < 0.5 Ha — W1e is NOT core-correlation-driven; "
                              f"deeper mechanism needed")
    elif correction < 0:
        verdict = "STOP_WRONG_SIGN"
        gate_outcomes.append(f"correction is ATTRACTIVE when it should be repulsive — "
                              f"physics broken")
    else:
        verdict = "INTERMEDIATE_PROCEED_FLAGGED"
        gate_outcomes.append(f"0.5 < sum < 2 Ha — partial magnitude; proceed but "
                              f"expect partial closure only")

    print(f"\nGATE VERDICT: {verdict}")
    for og in gate_outcomes:
        print(f"  {og}")
    print(f"\n  Predicted closure fraction: {correction/wall_depth*100:.1f}% of wall")
    print(f"  Predicted residual D_e if linear: {wall_depth - correction:+.4f} Ha")
    print(f"  Experimental NaH D_e: {exp_De} Ha")
    print(f"  Closure within 2x window [{closure_window_2x[0]}, "
          f"{closure_window_2x[1]}] requires residual = D_e_target +/- 2x")

    t1 = time.perf_counter()
    print(f"\nWall time: {t1 - t0:.2f} s")

    # ------------------------------------------------------------------
    # Save results
    # ------------------------------------------------------------------
    out = {
        "sprint": "F5 Step 1 — algebraic cross-block 2-body Coulomb diagnostic",
        "date": "2026-05-23",
        "R_AB_bohr": R_EQ_BOHR,
        "bonding_orbital_composition": {
            f"{c['side']}({c['n']},{c['l']},{c['m']})": c['c']
            for c in components
        },
        "spatial_norm_sq": norm_sq,
        "core_electron_count_check": N_core,
        "J_total_Ha": J_total,
        "K_per_orbital": K_per_orbital,
        "K_total_Ha": K_total,
        "predicted_correction_Ha": correction,
        "K_over_J_ratio": (K_total / J_total) if abs(J_total) > 1e-6 else None,
        "F3_baseline_wall_depth_Ha": wall_depth,
        "experimental_NaH_De_Ha": exp_De,
        "predicted_residual_De_if_linear_Ha": wall_depth - correction,
        "gate_verdict": verdict,
        "gate_outcomes": gate_outcomes,
        "closure_fraction_pct": correction / wall_depth * 100,
        "wall_time_s": t1 - t0,
    }
    out_path = os.path.join(_HERE, "data", "sprint_f5_step1_diagnostic.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
