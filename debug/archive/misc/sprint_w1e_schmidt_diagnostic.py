"""Sprint W1e Schmidt-orthogonalization Day-1 diagnostic
========================================================

Question (decision gate): does basis-level Schmidt orthogonalization of H 1s
against the Na [Ne] core (single-zeta Clementi-Raimondi) produce a TOTAL
diagonal-element differential |Delta h1 + Delta V_ne_screened| greater than
~100 mHa at NaH R = 3.5 bohr? If so, GO on a multi-week production-Schmidt
sprint. If 30-100 mHa, BORDERLINE (extend to Na 3s, 3p). If < 30 mHa, STOP.

This is a diagnostic-only computation. NO production code is modified.
NO FCI runs. Pure matrix-element / wavefunction-on-grid computation.

Setup
-----
- Na nucleus at origin (Z = 11).
- H nucleus at (R_AB, 0, 0) bohr, R_AB = 3.5.
- Hydrogenic H 1s (Z_orb = 1) centered on H.
- Na [Ne] core orbitals (1s, 2s, 2p_{-1, 0, +1}) centered on Na, single-zeta
  Clementi-Raimondi (1963), with Z_eff = n * zeta per Bethe-Salpeter conv.
- All orbitals in REAL spherical harmonics (framework convention).

Outputs
-------
The 5 cross-center overlaps S_c = <core_c | H 1s>_{cross-center}, the
total projection mass Sigma |S_c|^2, the diagonal-element differentials
on h1 (kinetic + V_H_nuc) and V_Na_screened, and the gate verdict.

Decision gate (frozen before computation):
  |Delta_total| > 100 mHa   => GO
  30 mHa <= |Delta_total| <= 100 mHa => BORDERLINE
  |Delta_total| < 30 mHa    => STOP
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy.special import genlaguerre, lpmv

# Reuse framework machinery
from geovac.phillips_kleinman_cross_center import (
    cross_center_overlap,
    _hydrogenic_radial_at_points,
    _real_sh_normalization,
    _core_orbitals_for_Z,
)
from geovac.neon_core import FrozenCore, _CLEMENTI_ZETA_NE


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

R_AB = 3.5  # NaH internuclear separation in bohr
Z_NA = 11
Z_H = 1

# Hydrogenic H 1s: (n=1, l=0, m=0), Z_eff = 1
H_1s_quantum = (1, 0, 0)
Z_orb_H = 1.0

# Output files
OUT_JSON = Path("debug/data/sprint_w1e_schmidt_diagnostic.json")
OUT_MEMO = Path("debug/sprint_w1e_schmidt_diagnostic_memo.md")


# ---------------------------------------------------------------------------
# Step 1: Load Na [Ne] core orbital descriptors
# ---------------------------------------------------------------------------

def get_core_orbitals():
    """Return list of (label, n, l, m, Z_eff_core, energy) for Na [Ne] core."""
    core_data = _core_orbitals_for_Z(Z_NA, 'Ne')
    # Each entry has 'n', 'l', 'm', 'zeta', 'energy'
    # Convert to (label, n, l, m, Z_eff, energy)
    label_map = {
        (1, 0, 0): '1s',
        (2, 0, 0): '2s',
        (2, 1, -1): '2p_-1',
        (2, 1, 0): '2p_0',
        (2, 1, +1): '2p_+1',
    }
    orbs = []
    for d in core_data:
        n = d['n']
        l = d['l']
        m = d['m']
        zeta = d['zeta']
        Z_eff = n * zeta  # framework convention (Bethe-Salpeter)
        label = label_map[(n, l, m)]
        orbs.append({
            'label': label,
            'n': n, 'l': l, 'm': m,
            'zeta': zeta, 'Z_eff': Z_eff,
            'energy': d['energy'],
        })
    return orbs


# ---------------------------------------------------------------------------
# Step 2: Compute 5 cross-center overlaps S_c = <core_c | H 1s>
# ---------------------------------------------------------------------------
# IMPORTANT: cross_center_overlap as currently written places "A" at origin
# and "B" at +R_AB * z_hat. So if we want H at +R_AB * x_hat and Na at origin
# (the physical NaH setup), we have an azimuthal-axis mismatch.
#
# Approach: use a z-aligned configuration (Na at origin, H at +R_AB * z_hat).
# This is identical physics; real-SH azimuthal selection is preserved because
# the axial symmetry is now along z. We compute overlaps with
#   "A" orbital = Na core orbital (centered on A at origin)
#   "B" orbital = H 1s (centered on B at +R_AB * z_hat).
# The signature of cross_center_overlap is:
#   <phi_p^A | phi_c^B>  where A is at origin, B at +R_AB z_hat.
# So we feed:
#   p = Na core orbital, c = H 1s.
# Real-SH azimuthal selection (m_p == m_c) gives nonzero only when the Na
# orbital has m=0, i.e. the 1s, 2s, 2p_0 orbitals. The 2p_+1 and 2p_-1
# give zero by symmetry (z-axial geometry).
# ---------------------------------------------------------------------------

def compute_overlaps(core_orbs, n_grid_r=8000, n_grid_u=128):
    """Compute the 5 cross-center overlaps S_c."""
    results = {}
    for c, orb in enumerate(core_orbs):
        S_c = cross_center_overlap(
            Z_orb_p=orb['Z_eff'],
            n_p=orb['n'], l_p=orb['l'], m_p=orb['m'],
            Z_orb_c=Z_orb_H,
            n_c=1, l_c=0, m_c=0,
            R_AB=R_AB,
            n_grid_r=n_grid_r,
            n_grid_u=n_grid_u,
        )
        results[orb['label']] = {
            'n': orb['n'], 'l': orb['l'], 'm': orb['m'],
            'Z_eff_core': orb['Z_eff'],
            'energy': orb['energy'],
            'S_c': float(S_c),
            'S_c_squared': float(S_c ** 2),
        }
    return results


# ---------------------------------------------------------------------------
# 3D-grid wavefunction construction for h1 and V_Na_screened diagonals
# ---------------------------------------------------------------------------
#
# For the diagonal differentials we need:
#   Delta h1 = <H1s_orth | T + V_H | H1s_orth> - <H1s | T + V_H | H1s>
#   Delta V_ne = <H1s_orth | V_Na_scr | H1s_orth> - <H1s | V_Na_scr | H1s>
#
# Since H 1s (Z=1) is an EXACT eigenstate of T + V_H with eigenvalue -1/2,
# we have <H1s | T + V_H | H1s> = -1/2 exactly. (No wavefunction tail beyond
# the radial integration matters here -- this is analytic.)
#
# For the orthogonalized state, we use the algebraic identity:
#   |H1s_orth> = (|H1s> - sum_c S_c |phi_c>) / sqrt(N)
# where N = 1 - sum_c |S_c|^2 (assuming the core orbitals are themselves
# orthonormal to each other -- which they are in the single-zeta Clementi-
# Raimondi approximation since they're separated by orthogonal (n,l,m)
# quantum numbers and use scaled hydrogenic radial functions for each n).
#
# Wait: in the C-R single-zeta convention, each shell uses a DIFFERENT
# Z_eff (n * zeta), so the 1s and 2s ARE NOT orthogonal in the
# strict hydrogenic sense; they would be orthogonal only in a self-
# consistent HF wavefunction. For the diagnostic, we'll proceed using
# THE FRAMEWORK'S [Ne] core (which is what production uses), check
# what <core_c | core_c'> looks like (should be ~0 for c != c'), and
# use the framework's actual <H1s|...> overlaps for the projection.
#
# The simplest correct approach: build a real-space grid representation of
# |H1s> on a 3D grid covering both centers; build |phi_c> for each Na core;
# build |H1s_orth> = |H1s> - sum_c S_c |phi_c> (unnormalized); then renormalize.
# Then compute matrix elements by numerical integration on the grid.
#
# 3D grid is expensive. Better: build wavefunctions on a 1D radial grid
# centered on Na and use spherical harmonic basis expansion. The H 1s is
# off-center; we expand it in spherical harmonics around Na via:
#   psi_H1s(r) = (1/sqrt(pi)) * exp(-|r - R_H|)
# where R_H = R_AB * z_hat. Use real-SH multipole expansion of this around
# the Na center.
#
# Easiest robust approach: use direct 3D Cartesian quadrature for the
# diagonal matrix elements. Coarse-but-converged 3D grid is feasible for
# a single diagnostic.
# ---------------------------------------------------------------------------


def hydrogenic_wavefunction_at_points(Z_eff, n, l, m, x, y, z, x0=0.0, y0=0.0, z0=0.0):
    """Hydrogenic wavefunction (with real-SH) evaluated at Cartesian points (x,y,z).

    Centered at (x0, y0, z0). Convention matches framework:
       psi_nlm(r) = R_nl(r; Z_eff) * Y_lm(theta, phi)
    where Y_lm is REAL spherical harmonic.

    R_nl normalization: int |R|^2 r^2 dr = 1.

    Returns array same shape as x.
    """
    dx = x - x0
    dy = y - y0
    dz = z - z0
    r = np.sqrt(dx * dx + dy * dy + dz * dz)
    # Avoid r=0 singularity for l=0 (R is finite but cos/lpmv stable)
    r_safe = np.where(r < 1e-12, 1e-12, r)

    # Radial part
    R_nl = _hydrogenic_radial_at_points(Z_eff, n, l, r_safe)

    # Angular part (real SH)
    # theta = arccos(dz/r), phi = arctan2(dy, dx)
    cos_theta = np.clip(dz / r_safe, -1.0, 1.0)
    phi = np.arctan2(dy, dx)

    abs_m = abs(m)
    # lpmv(abs_m, l, cos_theta): associated Legendre with Condon-Shortley
    P_lm = lpmv(abs_m, l, cos_theta)

    # Real-SH angular factor
    norm_sh = _real_sh_normalization(l, m)
    if m == 0:
        Y_lm = norm_sh * P_lm
    elif m > 0:
        Y_lm = norm_sh * P_lm * np.cos(abs_m * phi)
    else:  # m < 0
        Y_lm = norm_sh * P_lm * np.sin(abs_m * phi)

    return R_nl * Y_lm


def build_grid_3d(n_per_axis=80, extent=15.0):
    """Build a 3D Cartesian grid with a center near the H atom.

    Returns x, y, z arrays of shape (N,N,N) and dV (volume element).
    We center the grid box between Na (at origin) and H (at R_AB z_hat),
    so the box runs from z = -extent to z = R_AB + extent, with the
    transverse box symmetric.
    """
    # Box runs symmetrically around the midpoint of Na-H
    z_mid = R_AB / 2.0
    z_half = max(extent, R_AB / 2.0 + 4.0)
    # Use the same extent for transverse for simplicity
    x_lin = np.linspace(-extent, extent, n_per_axis)
    y_lin = np.linspace(-extent, extent, n_per_axis)
    z_lin = np.linspace(z_mid - z_half, z_mid + z_half, n_per_axis)
    dx = x_lin[1] - x_lin[0]
    dy = y_lin[1] - y_lin[0]
    dz = z_lin[1] - z_lin[0]
    dV = dx * dy * dz
    X, Y, Z = np.meshgrid(x_lin, y_lin, z_lin, indexing='ij')
    return X, Y, Z, dV


def build_grid_3d_two_scale(n_inner=60, r_inner=0.6, n_outer=70, extent=12.0):
    """Build a 3D Cartesian grid that has finer resolution near Na origin
    (to resolve the Na 1s/2s/2p high-Z core) AND covers the H atom at R_AB.

    Two grids are constructed and concatenated: an "inner" cube of size
    2*r_inner centered at Na with n_inner points/axis, and an "outer" 1D
    layout that fills the rest of the box up to extent. For simplicity and
    correctness of the FD Laplacian, we instead use a SINGLE uniform grid
    but choose dimensions to put many points inside the Na core.

    With n=140 per axis and extent=12, dx = 24/140 = 0.171 bohr. The Na 1s
    length scale is n/Z = 1/10.6 = 0.094 bohr, so we have ~1 point per
    core radius -- still under-resolved.

    A radial integration of psi_core itself wouldn't suffer this because
    it uses an adaptive r-grid in [0, infinity). For our diagnostic, the
    issue is that we need <psi_orth | O | psi_orth> on the same grid for
    BOTH psi_orth and psi_H. Since psi_orth includes a tiny admixture of
    the high-Z Na cores, the 3D grid is the natural shared space.

    Better strategy: increase n and use a tightly-cropped grid focused
    on the chemistry-relevant region (a few bohr around H), accept that
    the inner Na 1s peak is under-resolved, but check that the answer
    is stable.
    """
    return build_grid_3d(n_inner, extent)


def compute_diagonal_matrix_elements(core_orbs, S_dict, n_per_axis=80, extent=15.0):
    """Compute diagonal matrix elements on H1s and H1s_orth.

    Returns dict with:
      - h1_H1s, h1_H1s_orth   (kinetic + V_H_nuc)
      - vne_H1s, vne_H1s_orth (V_Na_screened)
    Plus norms as sanity check.
    """
    print(f"  Building {n_per_axis}^3 = {n_per_axis**3} Cartesian grid (extent={extent}, R_AB={R_AB})...")
    X, Y, Z, dV = build_grid_3d(n_per_axis, extent)

    # H 1s on H center
    print("  Building psi_H1s on grid...")
    psi_H = hydrogenic_wavefunction_at_points(
        Z_eff=Z_orb_H, n=1, l=0, m=0,
        x=X, y=Y, z=Z,
        x0=0.0, y0=0.0, z0=R_AB,  # H is at +R_AB z_hat
    )

    # Norm check for H1s
    norm_H = np.sum(psi_H * psi_H) * dV
    print(f"    <H1s | H1s>_grid = {norm_H:.6f} (should be ~1.0)")

    # Build all 5 Na core orbitals on grid
    print("  Building Na core orbitals on grid (1s, 2s, 2p_-1, 2p_0, 2p_+1)...")
    psi_cores = {}
    for orb in core_orbs:
        label = orb['label']
        psi_c = hydrogenic_wavefunction_at_points(
            Z_eff=orb['Z_eff'],
            n=orb['n'], l=orb['l'], m=orb['m'],
            x=X, y=Y, z=Z,
            x0=0.0, y0=0.0, z0=0.0,  # Na at origin
        )
        norm_c = np.sum(psi_c * psi_c) * dV
        print(f"    <{label} | {label}>_grid = {norm_c:.6f} (should be ~1.0)")
        psi_cores[label] = psi_c

    # Check overlap consistency: <core_c | H1s>_grid vs analytical S_c
    print("  Checking grid overlaps vs analytical:")
    for orb in core_orbs:
        label = orb['label']
        S_grid = np.sum(psi_cores[label] * psi_H) * dV
        S_anal = S_dict[label]['S_c']
        diff = abs(S_grid - S_anal)
        print(f"    {label}: S_grid={S_grid:+.4e}  S_anal={S_anal:+.4e}  |diff|={diff:.1e}")

    # Build H1s_orth (unnormalized): |H1s> - sum_c S_c |phi_c>, using the
    # analytical S_c values (more accurate than grid for the projection step)
    print("  Building H1s_orth = H1s - sum_c S_c phi_c ...")
    psi_orth = psi_H.copy()
    for orb in core_orbs:
        label = orb['label']
        S_c = S_dict[label]['S_c']
        psi_orth -= S_c * psi_cores[label]

    # Norm of orth (unnormalized): should be sqrt(1 - sum |S_c|^2)
    norm_orth_unnorm = np.sum(psi_orth * psi_orth) * dV
    sum_Sc_sq = sum(d['S_c_squared'] for d in S_dict.values())
    expected_norm_sq = 1.0 - sum_Sc_sq
    print(f"    <H1s_orth | H1s_orth>_unnorm = {norm_orth_unnorm:.6f}")
    print(f"    1 - sum |S_c|^2              = {expected_norm_sq:.6f}")

    # Normalize
    psi_orth_n = psi_orth / math.sqrt(norm_orth_unnorm)

    # ----- Compute h1 = T + V_H matrix elements -----
    # For H 1s (exact eigenstate of T - 1/r_H with eigenvalue -1/2):
    #   <H1s | T + V_H | H1s> = -1/2 exactly
    # Use grid-based numerical kinetic via psi*(-1/2 nabla^2) psi.
    # For numerical stability, compute via the action of -1/2 Laplacian
    # using 6-point stencil. Then add V_H = -1/|r - R_H|.

    print("  Computing T + V_H matrix elements...")
    # V_H = -1/|r - r_H|
    dist_H = np.sqrt(X**2 + Y**2 + (Z - R_AB)**2)
    dist_H_safe = np.where(dist_H < 1e-3, 1e-3, dist_H)  # avoid singularity
    V_H = -1.0 / dist_H_safe

    # Kinetic via finite difference Laplacian. Grid spacing:
    dx = X[1, 0, 0] - X[0, 0, 0]
    dy = Y[0, 1, 0] - Y[0, 0, 0]
    dz = Z[0, 0, 1] - Z[0, 0, 0]
    # We need -1/2 nabla^2 psi acting on psi.
    # Use second-order centered FD; boundary tapered to zero.

    def laplacian(psi):
        lap = np.zeros_like(psi)
        lap[1:-1, :, :] += (psi[2:, :, :] - 2 * psi[1:-1, :, :] + psi[:-2, :, :]) / dx**2
        lap[:, 1:-1, :] += (psi[:, 2:, :] - 2 * psi[:, 1:-1, :] + psi[:, :-2, :]) / dy**2
        lap[:, :, 1:-1] += (psi[:, :, 2:] - 2 * psi[:, :, 1:-1] + psi[:, :, :-2]) / dz**2
        return lap

    # h1 on H1s
    T_psi_H = -0.5 * laplacian(psi_H)
    h1_H = np.sum(psi_H * T_psi_H) * dV + np.sum(psi_H * V_H * psi_H) * dV

    T_psi_orth = -0.5 * laplacian(psi_orth_n)
    h1_orth = np.sum(psi_orth_n * T_psi_orth) * dV + np.sum(psi_orth_n * V_H * psi_orth_n) * dV

    print(f"    <H1s | h1 | H1s>            = {h1_H:.6f} Ha (exact: -0.500000)")
    print(f"    <H1s_orth | h1 | H1s_orth>  = {h1_orth:.6f} Ha")

    # ----- Compute V_Na_screened matrix elements -----
    print("  Computing V_Na_screened matrix elements...")
    # V_Na_screened(r) = -Z_eff(|r|) / |r|
    fc = FrozenCore(Z=Z_NA, core_type='Ne')
    fc.solve()
    r_Na = np.sqrt(X**2 + Y**2 + Z**2)
    r_Na_safe = np.where(r_Na < 1e-3, 1e-3, r_Na)
    # z_eff works elementwise on arrays
    z_eff_arr = fc.z_eff(r_Na_safe)
    V_Na_scr = -z_eff_arr / r_Na_safe
    # At r_Na -> 0, Z_eff -> Z = 11, so V_Na_scr -> -11/r (singular, contained by mask)

    vne_H = np.sum(psi_H * V_Na_scr * psi_H) * dV
    vne_orth = np.sum(psi_orth_n * V_Na_scr * psi_orth_n) * dV
    print(f"    <H1s | V_Na_scr | H1s>            = {vne_H:+.6f} Ha")
    print(f"    <H1s_orth | V_Na_scr | H1s_orth>  = {vne_orth:+.6f} Ha")

    return {
        'norm_H1s_grid': float(norm_H),
        'norm_H1s_orth_unnorm': float(norm_orth_unnorm),
        'sum_Sc_squared': float(sum_Sc_sq),
        'expected_orth_norm_sq': float(expected_norm_sq),
        'h1_H1s': float(h1_H),
        'h1_H1s_orth': float(h1_orth),
        'vne_H1s': float(vne_H),
        'vne_H1s_orth': float(vne_orth),
    }


# ---------------------------------------------------------------------------
# Decision gate
# ---------------------------------------------------------------------------

def gate_verdict(delta_total_mHa):
    """Apply the decision gate frozen before computation."""
    abs_delta = abs(delta_total_mHa)
    if abs_delta > 100.0:
        return 'GO', 'Schmidt sprint worth multi-week commitment.'
    elif abs_delta >= 30.0:
        return 'BORDERLINE', 'Extend diagnostic to Na 3s + 3p before commitment.'
    else:
        return 'STOP', 'Wall is not basis-level Schmidt; report honestly.'


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    t0 = time.time()
    print("=" * 78)
    print("Sprint W1e Schmidt-orthogonalization Day-1 diagnostic")
    print("=" * 78)
    print(f"  Geometry: Na (Z={Z_NA}) at origin, H (Z={Z_H}) at +{R_AB}*z_hat bohr")
    print(f"  Convention: z-aligned (computational; physically equivalent to x-aligned)")
    print(f"  [Ne] core: single-zeta Clementi-Raimondi (1963)")
    print()

    # Step 1: Load core orbital descriptors
    print("Step 1: Loading Na [Ne] core orbital descriptors...")
    core_orbs = get_core_orbitals()
    for orb in core_orbs:
        print(f"  {orb['label']}: n={orb['n']}, l={orb['l']}, m={orb['m']}, "
              f"zeta={orb['zeta']:.4f}, Z_eff={orb['Z_eff']:.4f}, E={orb['energy']:+.4f} Ha")
    print()

    # Step 2: Compute 5 cross-center overlaps
    print("Step 2: Computing cross-center overlaps S_c = <core_c | H 1s>...")
    S_dict = compute_overlaps(core_orbs)
    for label, d in S_dict.items():
        m_signed = d['m']
        if m_signed == 0:
            tag = '(nonzero by symmetry)'
        else:
            tag = '(zero by m-mismatch symmetry)'
        print(f"  S({label}) = {d['S_c']:+.6e}   |S|^2 = {d['S_c_squared']:.6e}   {tag}")
    sum_Sc_sq = sum(d['S_c_squared'] for d in S_dict.values())
    print(f"  Projection mass: sum |S_c|^2 = {sum_Sc_sq:.6e}")
    print(f"  Orthogonalization norm:        sqrt(1 - sum |S_c|^2) = {math.sqrt(1.0 - sum_Sc_sq):.6f}")
    print()

    # Step 3-4: Compute diagonal matrix elements via 3D grid quadrature
    # Convergence sweep over grid resolution
    print("Step 3-4: Computing diagonal matrix elements on 3D Cartesian grid...")
    print("  Convergence sweep over (n_per_axis, extent)...")
    diag_runs = []
    for (n_pa, ext) in [(80, 12.0), (100, 10.0), (120, 8.0), (140, 8.0)]:
        print(f"\n  === Grid ({n_pa}, ext={ext}) ===")
        try:
            d = compute_diagonal_matrix_elements(core_orbs, S_dict,
                                                 n_per_axis=n_pa, extent=ext)
        except MemoryError:
            print(f"    MemoryError at {n_pa}; skipping.")
            continue
        d_h1 = d['h1_H1s_orth'] - d['h1_H1s']
        d_vne = d['vne_H1s_orth'] - d['vne_H1s']
        d_tot = d_h1 + d_vne
        print(f"    => Delta h1 = {d_h1*1000:+.3f} mHa  "
              f"Delta V_ne = {d_vne*1000:+.3f} mHa  "
              f"Delta total = {d_tot*1000:+.3f} mHa")
        d['grid_n'] = n_pa
        d['grid_extent'] = ext
        diag_runs.append(d)

    # Use the finest grid for the headline number, but report all runs
    diag = diag_runs[-1] if diag_runs else None
    print()
    print("  Convergence summary:")
    print(f"  {'n':>4} {'ext':>6} {'<H|h1|H>':>11} {'<orth|h1|orth>':>15} "
          f"{'<H|Vne|H>':>11} {'<orth|Vne|orth>':>15} {'Dh1(mHa)':>10} {'DVne(mHa)':>11} {'Dtot(mHa)':>11}")
    for d in diag_runs:
        d_h1 = d['h1_H1s_orth'] - d['h1_H1s']
        d_vne = d['vne_H1s_orth'] - d['vne_H1s']
        d_tot = d_h1 + d_vne
        print(f"  {d['grid_n']:>4} {d['grid_extent']:>6.1f} "
              f"{d['h1_H1s']:>11.5f} {d['h1_H1s_orth']:>15.5f} "
              f"{d['vne_H1s']:>11.5f} {d['vne_H1s_orth']:>15.5f} "
              f"{d_h1*1000:>10.3f} {d_vne*1000:>11.3f} {d_tot*1000:>11.3f}")
    print()

    # Differentials
    delta_h1 = diag['h1_H1s_orth'] - diag['h1_H1s']
    delta_vne = diag['vne_H1s_orth'] - diag['vne_H1s']
    delta_total = delta_h1 + delta_vne
    delta_h1_mHa = delta_h1 * 1000.0
    delta_vne_mHa = delta_vne * 1000.0
    delta_total_mHa = delta_total * 1000.0

    print("Step 5: Differentials")
    print(f"  Delta h1   = {delta_h1:+.6f} Ha  =  {delta_h1_mHa:+.3f} mHa")
    print(f"  Delta V_ne = {delta_vne:+.6f} Ha  =  {delta_vne_mHa:+.3f} mHa")
    print(f"  Delta total= {delta_total:+.6f} Ha  =  {delta_total_mHa:+.3f} mHa")
    print()

    verdict, comment = gate_verdict(delta_total_mHa)
    print(f"DECISION GATE: {verdict}  ({comment})")
    print(f"  Threshold: |Delta_total| > 100 mHa => GO; 30-100 => BORDERLINE; < 30 => STOP")
    print()

    elapsed = time.time() - t0
    print(f"Total elapsed: {elapsed:.1f} s")
    print()

    # ----- Save JSON output -----
    out = {
        'sprint': 'W1e Schmidt orthogonalization Day-1 diagnostic',
        'date': '2026-05-23',
        'geometry': {
            'Na_at': [0.0, 0.0, 0.0],
            'H_at':  [0.0, 0.0, R_AB],
            'R_AB_bohr': R_AB,
            'z_aligned': True,
            'note': 'z-aligned for simplicity; physically equivalent to x-aligned NaH',
        },
        'core_orbitals': {orb['label']: {k: orb[k] for k in ('n', 'l', 'm', 'zeta', 'Z_eff', 'energy')}
                          for orb in core_orbs},
        'overlap_signed': {label: d['S_c'] for label, d in S_dict.items()},
        'overlap_squared': {label: d['S_c_squared'] for label, d in S_dict.items()},
        'projection_norm_squared': float(sum_Sc_sq),
        'expected_orth_norm_unnorm': float(1.0 - sum_Sc_sq),
        'grid_convergence_sweep': [
            {**{k: float(v) if isinstance(v, (int, float, np.floating)) else v
                for k, v in d.items()},
             'delta_h1_Ha': float(d['h1_H1s_orth'] - d['h1_H1s']),
             'delta_v_ne_Ha': float(d['vne_H1s_orth'] - d['vne_H1s']),
             'delta_total_Ha': float((d['h1_H1s_orth'] - d['h1_H1s']) + (d['vne_H1s_orth'] - d['vne_H1s'])),
             }
            for d in diag_runs
        ],
        'finest_grid': diag,
        'delta_h1_Ha': float(delta_h1),
        'delta_v_ne_Ha': float(delta_vne),
        'delta_total_Ha': float(delta_total),
        'delta_h1_mHa': float(delta_h1_mHa),
        'delta_v_ne_mHa': float(delta_vne_mHa),
        'delta_total_mHa': float(delta_total_mHa),
        'gate_verdict': verdict,
        'gate_comment': comment,
        'elapsed_seconds': float(elapsed),
    }
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with OUT_JSON.open('w') as f:
        json.dump(out, f, indent=2)
    print(f"  Saved JSON: {OUT_JSON}")

    return out


if __name__ == '__main__':
    main()
