"""
Phase 4: H2 SCF via Orbital Exponent Optimization
===================================================
The frozen H2+ orbitals are too compact for H2 (J_gg = 0.78 vs optimal ~0.63).
Fix: vary the effective nuclear charge Z_eff used to generate the orbital,
then evaluate the energy with the physical Z=1 Hamiltonian.

This is the Eckart variational approach:
  phi(Z_eff) = 1sigma_g orbital of (Z_eff, Z_eff) problem at distance R
  h(Z_eff) = eps(Z_eff) + (Z_eff - 1) * <1/r_A + 1/r_B>
  E_HF(Z_eff) = 2*h + J_gg + V_NN

Minimize E_HF(Z_eff) over Z_eff to get the SCF energy.

Then run CI on top of the optimized orbital for correlation energy.

Date: 2026-03-13
"""

import warnings
warnings.filterwarnings('ignore')

import sys
import os
if sys.stdout.encoding and sys.stdout.encoding.lower().startswith('cp'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.linalg import eigh
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.prolate_spheroidal_lattice import (
    ProlateSpheroidalLattice,
    fit_spectroscopic_constants,
)
from debug.stress_test_prolate_heh_plus import get_orbital_on_grid_general
from debug.stress_test_prolate_h2 import compute_vee_integral


# ============================================================
# Nuclear attraction expectation value
# ============================================================
def compute_vnuc_expectation(orb: dict) -> float:
    """Compute <phi|1/r_A + 1/r_B|phi> on the quadrature grid.

    Uses the identity:
      (1/r_A + 1/r_B) * J = R^2 * xi / 2
    where J = (R/2)^3 * (xi^2 - eta^2) is the Jacobian.

    This cancellation makes the integral very stable numerically.
    """
    R = orb['R']
    xi = orb['xi']
    eta = orb['eta']
    w_xi = orb['w_xi']
    w_eta = orb['w_eta']
    psi = orb['psi']

    # V_nuc * J = R^2 * xi / 2 (the (xi^2 - eta^2) cancels!)
    XI = xi[:, np.newaxis]  # (N_xi, 1)
    vnuc_J = R**2 * XI / 2  # (N_xi, 1) broadcasts over eta

    W_XI = w_xi[:, np.newaxis]
    W_ETA = w_eta[np.newaxis, :]

    # <V_nuc> = 2*pi * sum(psi^2 * vnuc_J * w_xi * w_eta)
    result = 2 * np.pi * np.sum(psi**2 * vnuc_J * W_XI * W_ETA)
    return result


# ============================================================
# SCF energy at given Z_eff
# ============================================================
def compute_scf_energy(
    R: float,
    Z_eff: float,
    N_xi_solve: int = 5000,
    N_grid: int = 50,
    xi_max_grid: float = 15.0,
) -> dict:
    """Compute H2 HF energy using orbital generated at (Z_eff, Z_eff).

    h_gg = eps(Z_eff) + (Z_eff - 1) * <1/r_A + 1/r_B>
    E_HF = 2*h_gg + J_gg + 1/R

    For Z_eff=1: recovers the frozen H2+ result.
    For Z_eff<1: more diffuse orbital, lower J_gg.
    """
    # Generate orbital at (Z_eff, Z_eff)
    # Z_A, Z_B must be int for the solver, but we can use float Z
    # by scaling: solve at Z=1 with R_eff = R*Z_eff, then rescale
    # Actually, the solver takes int Z. For fractional Z_eff we need
    # a workaround.
    #
    # Workaround: the prolate spheroidal equation depends on
    #   a = R*(Z_A + Z_B) and c^2 = -R^2*E/2
    # If we set Z_A = Z_B = 1 but use R_eff = R*Z_eff, then
    #   a_eff = R_eff * 2 = 2*R*Z_eff = R*(Z_eff + Z_eff)
    # This is equivalent to Z_A = Z_B = Z_eff at distance R.
    #
    # The energy scales as: E(R_eff) = E_elec at R_eff
    # And h_gg = E(R_eff) + (Z_eff - 1) * <1/r_A + 1/r_B>
    #
    # But wait — the orbital shape depends on R_eff, and the
    # grid coordinates are in terms of R_eff too. The V_ee integral
    # and V_nuc expectation need to use the physical R.
    #
    # Cleaner approach: directly modify the solver parameters.
    # Since ProlateSpheroidalLattice uses _a = R*(Z_A+Z_B) internally,
    # and we want _a = R*(Z_eff + Z_eff) = 2*R*Z_eff, we can
    # pass R_scaled = R * Z_eff with Z_A=Z_B=1.

    R_eff = R * Z_eff

    try:
        orb = get_orbital_on_grid_general(
            R=R_eff, Z_A=1, Z_B=1, n_angular=0,
            N_xi_solve=N_xi_solve,
            N_xi_grid=N_grid, N_eta_grid=N_grid,
            xi_max_grid=xi_max_grid,
        )
    except Exception as e:
        return {'R': R, 'Z_eff': Z_eff, 'E_HF': np.nan, 'error': str(e)}

    # The orbital is defined in (xi, eta) coordinates of R_eff.
    # But the physical system has internuclear distance R.
    # In prolate spheroidal coords:
    #   r_A = R*xi*eta ... no, r_A = R*(xi - eta)/2, r_B = R*(xi + eta)/2
    #
    # If the orbital was generated at R_eff, its (xi,eta) grid maps to
    # physical coordinates via R_eff. To evaluate integrals at physical R,
    # we need the orbital on the R grid.
    #
    # Actually, the scaling is: if we solve at (Z_eff, Z_eff, R), the
    # equation parameter a = R*(2*Z_eff). This is the same as solving
    # at (1, 1, R*Z_eff) with a = R_eff*2. The ORBITAL SHAPE in (xi,eta)
    # is identical, but the physical coordinates are scaled:
    #   r_phys = (R/2) * xi * eta  vs  r_eff = (R_eff/2) * xi * eta
    #
    # For computing V_ee integrals, we need the orbital on the PHYSICAL
    # grid (R), not the R_eff grid. But the (xi,eta) shape is the same.
    # We just need to adjust R in the Jacobian and cylindrical coordinates.
    #
    # Simplest correct approach: overwrite orb['R'] = R (physical distance)
    # and re-normalize. The psi values on the grid are the same shape,
    # the Jacobian changes with R.

    # Re-normalize with physical R
    orb_phys = dict(orb)
    orb_phys['R'] = R

    XI, ETA = np.meshgrid(orb['xi'], orb['eta'], indexing='ij')
    J_phys = (R / 2)**3 * (XI**2 - ETA**2)
    W_XI, W_ETA = np.meshgrid(orb['w_xi'], orb['w_eta'], indexing='ij')
    norm_sq = np.sum(orb['psi']**2 * J_phys * W_XI * W_ETA) * 2 * np.pi
    orb_phys['psi'] = orb['psi'] / np.sqrt(norm_sq)

    # eps(Z_eff) from the solver (electronic energy at R_eff)
    eps_zeff = orb['E_elec']

    # <1/r_A + 1/r_B> with PHYSICAL R
    vnuc_exp = compute_vnuc_expectation(orb_phys)

    # One-electron energy with physical Z=1:
    # h = T + V_nuc(Z=1) = [T + V_nuc(Z_eff)] + (1 - Z_eff)*V_nuc/Z_eff
    # But eps(Z_eff) = T + Z_eff * V_nuc_exp_eff
    #
    # Actually, this is trickier with the R_eff scaling. Let me think...
    #
    # The orbital phi satisfies: [T(R_eff) + V_nuc(Z=1, R_eff)] phi = eps phi
    # where T(R_eff) uses R_eff in the kinetic energy.
    #
    # But the physical Hamiltonian has T(R) and V_nuc(Z=1, R).
    #
    # For the Eckart approach, we use phi as a trial function:
    # <phi|H_phys|phi> = <phi|T(R)|phi> + <phi|V_nuc(Z=1,R)|phi>
    #
    # Since phi was generated at R_eff, the kinetic energy involves
    # R_eff. The kinetic energy in prolate coords scales as 1/R^2.
    #
    # Let me use a different, more direct approach.

    # Direct approach: compute h_gg = <phi|T+V_nuc|phi> by quadrature
    # T = total - V_nuc, but that's circular.
    #
    # Better: use the virial-like relation. For the H2+ eigenstate at R_eff:
    #   eps(R_eff) = <T(R_eff)> + <V_nuc(Z=1, R_eff)>
    #
    # The physical one-electron Hamiltonian:
    #   h_phys = T(R) + V_nuc(Z=1, R)
    #
    # In prolate spheroidal, T depends on R through the 2/R^2 prefactor.
    # <T(R)> = (R_eff/R)^2 * <T(R_eff)> [since T ~ 1/R^2 and same (xi,eta) grid]
    #
    # And V_nuc(Z=1, R) = -(1/r_A + 1/r_B) with r_A = R(xi-eta)/2.
    # <V_nuc(Z=1,R)> = -<1/r_A(R) + 1/r_B(R)> = -(R_eff/R) * <1/r_A(R_eff) + 1/r_B(R_eff)>
    # Wait, 1/r_A = 2/(R(xi-eta)), so <1/r_A> scales as 1/R.
    # <V_nuc(R)> = (R_eff/R) * <V_nuc(R_eff)>
    #
    # And <T(R_eff)> = eps(R_eff) - <V_nuc(R_eff)>
    #
    # So: h_phys = (R_eff/R)^2 * <T(R_eff)> + (R_eff/R) * <V_nuc(R_eff)>

    # <V_nuc(R_eff)> = -<1/r_A(R_eff) + 1/r_B(R_eff)>
    # We need this at R_eff, which is just vnuc from the orb (not orb_phys)
    orb_reff = dict(orb)
    orb_reff['R'] = R_eff
    # orb is already at R_eff and normalized there
    vnuc_reff = compute_vnuc_expectation(orb)  # <1/r_A + 1/r_B> at R_eff

    T_reff = eps_zeff + vnuc_reff  # <T> = eps - <V_nuc> = eps - (-vnuc) = eps + vnuc

    ratio = Z_eff  # R_eff / R = Z_eff
    h_phys = ratio**2 * T_reff - ratio * vnuc_reff
    # = Z_eff^2 * <T(R_eff)> - Z_eff * <1/r_A(R_eff) + 1/r_B(R_eff)>

    # Compute J_gg with physical R
    J_gg = compute_vee_integral(orb_phys, orb_phys, orb_phys, orb_phys)

    E_HF = 2 * h_phys + J_gg + 1.0 / R

    return {
        'R': R,
        'Z_eff': Z_eff,
        'E_HF': E_HF,
        'h_gg': h_phys,
        'J_gg': J_gg,
        'eps_zeff': eps_zeff,
        'T_reff': T_reff,
        'vnuc_reff': vnuc_reff,
        'vnuc_phys': vnuc_exp,
    }


# ============================================================
# Optimize Z_eff at given R
# ============================================================
def optimize_zeff(
    R: float,
    N_xi_solve: int = 5000,
    N_grid: int = 50,
    xi_max_grid: float = 15.0,
    verbose: bool = True,
) -> dict:
    """Find optimal Z_eff that minimizes E_HF(Z_eff) at given R."""
    t0 = time.time()

    if verbose:
        print(f"\n  Optimizing Z_eff at R={R:.3f} bohr")

    # Coarse scan to find bracket
    zeff_vals = np.linspace(0.7, 1.5, 17)
    energies = []
    for z in zeff_vals:
        res = compute_scf_energy(R, z, N_xi_solve, N_grid, xi_max_grid)
        energies.append(res['E_HF'])
        if verbose:
            print(f"    Z_eff={z:.3f}: E_HF={res['E_HF']:.6f} Ha  "
                  f"J_gg={res.get('J_gg', np.nan):.4f}")

    energies = np.array(energies)
    valid = ~np.isnan(energies)
    if np.sum(valid) < 3:
        return {'R': R, 'Z_eff_opt': np.nan, 'E_HF_opt': np.nan}

    idx_min = np.nanargmin(energies)

    # Refine with golden section
    if 0 < idx_min < len(zeff_vals) - 1:
        z_lo = zeff_vals[max(0, idx_min - 1)]
        z_hi = zeff_vals[min(len(zeff_vals)-1, idx_min + 1)]
    else:
        z_lo, z_hi = zeff_vals[valid][0], zeff_vals[valid][-1]

    def obj(z):
        res = compute_scf_energy(R, z, N_xi_solve, N_grid, xi_max_grid)
        return res['E_HF'] if not np.isnan(res['E_HF']) else 1e10

    opt = minimize_scalar(obj, bounds=(z_lo, z_hi), method='bounded',
                          options={'xatol': 0.001})
    z_opt = opt.x
    res_opt = compute_scf_energy(R, z_opt, N_xi_solve, N_grid, xi_max_grid)

    # Also get frozen (Z_eff=1) result for comparison
    res_frozen = compute_scf_energy(R, 1.0, N_xi_solve, N_grid, xi_max_grid)

    dt = time.time() - t0

    if verbose:
        print(f"\n  Z_eff* = {z_opt:.4f}")
        print(f"  E_HF(Z_eff*) = {res_opt['E_HF']:.6f} Ha")
        print(f"  E_HF(Z=1)    = {res_frozen['E_HF']:.6f} Ha")
        print(f"  Improvement:   {res_frozen['E_HF'] - res_opt['E_HF']:.6f} Ha")
        print(f"  J_gg(Z_eff*) = {res_opt['J_gg']:.6f} Ha")
        print(f"  J_gg(Z=1)    = {res_frozen['J_gg']:.6f} Ha")
        print(f"  Time: {dt:.1f}s")

    return {
        'R': R,
        'Z_eff_opt': z_opt,
        'E_HF_opt': res_opt['E_HF'],
        'E_HF_frozen': res_frozen['E_HF'],
        'J_gg_opt': res_opt['J_gg'],
        'J_gg_frozen': res_frozen['J_gg'],
        'h_gg_opt': res_opt['h_gg'],
        'h_gg_frozen': res_frozen['h_gg'],
        'time': dt,
    }


# ============================================================
# SCF + CI at optimized Z_eff
# ============================================================
def h2_scf_ci(
    R: float,
    N_xi_solve: int = 5000,
    N_grid: int = 50,
    xi_max_grid: float = 15.0,
    verbose: bool = True,
) -> dict:
    """Full SCF + CI calculation at given R.

    1. Optimize Z_eff for the bonding orbital
    2. Generate bonding + antibonding orbitals at Z_eff*
    3. Run 2x2 CI on top
    """
    t0 = time.time()

    # Step 1: Optimize Z_eff
    opt = optimize_zeff(R, N_xi_solve, N_grid, xi_max_grid, verbose=verbose)
    if np.isnan(opt['Z_eff_opt']):
        return {'R': R, 'E_total': np.nan}

    z_opt = opt['Z_eff_opt']
    R_eff = R * z_opt

    # Step 2: Generate both orbitals at Z_eff*
    if verbose:
        print(f"\n  Generating orbitals at Z_eff={z_opt:.4f} (R_eff={R_eff:.4f})")

    orbs = {}
    for label, n_ang in [('g', 0), ('u', 1)]:
        orb = get_orbital_on_grid_general(
            R=R_eff, Z_A=1, Z_B=1, n_angular=n_ang,
            N_xi_solve=N_xi_solve,
            N_xi_grid=N_grid, N_eta_grid=N_grid,
            xi_max_grid=xi_max_grid,
        )
        # Re-normalize with physical R
        XI, ETA = np.meshgrid(orb['xi'], orb['eta'], indexing='ij')
        J_phys = (R / 2)**3 * (XI**2 - ETA**2)
        W_XI, W_ETA = np.meshgrid(orb['w_xi'], orb['w_eta'], indexing='ij')
        norm_sq = np.sum(orb['psi']**2 * J_phys * W_XI * W_ETA) * 2 * np.pi
        orb['psi'] = orb['psi'] / np.sqrt(norm_sq)
        orb['R'] = R  # Use physical R for integrals
        orbs[label] = orb

    # Step 3: Compute integrals
    J_gg = compute_vee_integral(orbs['g'], orbs['g'], orbs['g'], orbs['g'])
    J_uu = compute_vee_integral(orbs['u'], orbs['u'], orbs['u'], orbs['u'])
    V_gguu = compute_vee_integral(orbs['g'], orbs['g'], orbs['u'], orbs['u'])

    if verbose:
        print(f"    J_gg = {J_gg:.6f}, J_uu = {J_uu:.6f}, V_gguu = {V_gguu:.6f}")

    # Step 4: Build CI matrix
    # One-electron energies: use the scaling relation
    for label in ['g', 'u']:
        orb = orbs[label]
        # Recompute vnuc at R_eff for this orbital
        orb_reff = dict(orb)
        orb_reff['R'] = R_eff
        # Need to re-normalize at R_eff
        J_reff = (R_eff / 2)**3 * (XI**2 - ETA**2)
        norm_reff = np.sum(orb['psi']**2 * J_reff * W_XI * W_ETA) * 2 * np.pi
        psi_reff = orb['psi'] * np.sqrt(
            np.sum(orb['psi']**2 * J_phys * W_XI * W_ETA) * 2 * np.pi / norm_reff
        )
        orb_reff['psi'] = psi_reff

        # Actually, simpler: the same psi on grid, different R for vnuc
        vnuc_tmp = dict(orb)
        vnuc_tmp['R'] = R_eff
        # Renormalize for R_eff
        norm_reff = np.sum(orb['psi']**2 * (R_eff/2)**3 * (XI**2 - ETA**2) * W_XI * W_ETA) * 2 * np.pi
        vnuc_tmp['psi'] = orb['psi'] / np.sqrt(norm_reff) * np.sqrt(
            np.sum(orb['psi']**2 * J_phys * W_XI * W_ETA) * 2 * np.pi
        )

    # Simpler: just use h from the optimize step for g, and compute for u similarly
    # For the antibonding orbital:
    orb_u_reff = get_orbital_on_grid_general(
        R=R_eff, Z_A=1, Z_B=1, n_angular=1,
        N_xi_solve=N_xi_solve,
        N_xi_grid=N_grid, N_eta_grid=N_grid,
        xi_max_grid=xi_max_grid,
    )
    eps_u_zeff = orb_u_reff['E_elec']
    vnuc_u_reff = compute_vnuc_expectation(orb_u_reff)
    T_u_reff = eps_u_zeff + vnuc_u_reff
    h_u_phys = z_opt**2 * T_u_reff - z_opt * vnuc_u_reff

    h_g_phys = opt['h_gg_opt']

    H_CI = np.array([
        [2*h_g_phys + J_gg, V_gguu],
        [V_gguu, 2*h_u_phys + J_uu],
    ])

    evals, evecs = eigh(H_CI)
    E_CI_elec = evals[0]
    V_NN = 1.0 / R
    E_total = E_CI_elec + V_NN
    E_HF = 2*h_g_phys + J_gg + V_NN

    dt = time.time() - t0

    if verbose:
        print(f"\n  h_g = {h_g_phys:.6f}, h_u = {h_u_phys:.6f}")
        print(f"  CI matrix:")
        print(f"    H_11 = {H_CI[0,0]:.6f}, H_22 = {H_CI[1,1]:.6f}, H_12 = {H_CI[0,1]:.6f}")
        print(f"  E_HF(SCF)  = {E_HF:.6f} Ha")
        print(f"  E_CI(SCF)  = {E_total:.6f} Ha")
        print(f"  E_HF(froz) = {opt['E_HF_frozen']:.6f} Ha")
        print(f"  Total time: {dt:.1f}s")

    return {
        'R': R,
        'E_total': E_total,
        'E_HF_scf': E_HF,
        'E_HF_frozen': opt['E_HF_frozen'],
        'Z_eff_opt': z_opt,
        'J_gg': J_gg,
        'h_gg': h_g_phys,
        'time': dt,
    }


# ============================================================
# Main
# ============================================================
def main():
    print("=" * 70)
    print("  H2 SCF VIA ORBITAL EXPONENT OPTIMIZATION")
    print("  Eckart variational: minimize E_HF(Z_eff)")
    print("  Date: 2026-03-13")
    print("=" * 70)

    E_H = -0.5
    E_H2_exact = -1.1745
    R_eq_exact = 1.401
    D_e_exact = 0.1745

    # --- Test 1: Z_eff optimization at R=1.4 ---
    print("\n" + "="*60)
    print("  TEST 1: Z_eff optimization at R=1.4")
    print("="*60)
    opt = optimize_zeff(R=1.4, N_xi_solve=5000, N_grid=50, verbose=True)

    # --- Test 2: SCF + CI at R=1.4 ---
    print("\n" + "="*60)
    print("  TEST 2: SCF + CI at R=1.4")
    print("="*60)
    res = h2_scf_ci(R=1.4, N_xi_solve=5000, N_grid=50, verbose=True)

    if not np.isnan(res['E_total']):
        D_e_scf = 2*E_H - res['E_total']
        print(f"\n  D_e(SCF+CI) = {D_e_scf:.6f} Ha ({D_e_scf/D_e_exact*100:.1f}% of exact)")

    # --- Test 3: PES scan ---
    print("\n" + "="*60)
    print("  TEST 3: H2 SCF PES scan")
    print("="*60)
    R_vals = np.array([0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0])
    results = []
    for R in R_vals:
        print(f"\n  --- R = {R:.2f} ---")
        try:
            res = h2_scf_ci(R=R, N_xi_solve=5000, N_grid=40, verbose=False)
            results.append(res)
            De = 2*E_H - res['E_total']
            print(f"  Z_eff*={res.get('Z_eff_opt',np.nan):.3f}  "
                  f"E_CI={res['E_total']:.6f}  D_e={De:.6f}  "
                  f"J_gg={res.get('J_gg',np.nan):.4f}  [{res['time']:.0f}s]")
        except Exception as e:
            print(f"  FAILED: {e}")
            results.append({'R': R, 'E_total': np.nan, 'E_HF_frozen': np.nan})

    R_arr = np.array([r['R'] for r in results])
    E_scf = np.array([r['E_total'] for r in results])
    E_froz = np.array([r.get('E_HF_frozen', np.nan) for r in results])

    print(f"\n  {'R':>6s} {'E_SCF+CI':>12s} {'E_frozen':>12s} {'D_e_SCF':>10s} {'Z_eff*':>8s}")
    print("-" * 55)
    for i, r in enumerate(results):
        R = r['R']
        E = r['E_total']
        Ef = r.get('E_HF_frozen', np.nan)
        De = 2*E_H - E
        z = r.get('Z_eff_opt', np.nan)
        print(f"  {R:6.3f} {E:12.6f} {Ef:12.6f} {De:10.6f} {z:8.4f}")

    valid = ~np.isnan(E_scf)
    if np.sum(valid) >= 5:
        fit_scf = fit_spectroscopic_constants(R_arr[valid], E_scf[valid])

        D_scf = 2*E_H - fit_scf['E_min']
        R_err = abs(fit_scf['R_eq'] - R_eq_exact) / R_eq_exact * 100
        D_err = abs(D_scf - D_e_exact) / D_e_exact * 100

        print(f"\n  SCF+CI spectroscopic constants:")
        print(f"    R_eq = {fit_scf['R_eq']:.4f} bohr (exact {R_eq_exact}, err {R_err:.1f}%)")
        print(f"    D_e  = {D_scf:.6f} Ha (exact {D_e_exact}, err {D_err:.1f}%)")

        print(f"\n  {'='*55}")
        print(f"  SCF IMPROVEMENT OVER FROZEN ORBITALS:")
        print(f"    D_e: {D_scf:.4f} Ha vs 0.072 Ha (frozen)")
        print(f"    R_eq: {fit_scf['R_eq']:.3f} vs 1.787 (frozen)")
        print(f"  {'='*55}")

    # Save results
    outpath = os.path.join(os.path.dirname(__file__), 'data', 'prolate_h2_scf.txt')
    with open(outpath, 'w', encoding='utf-8') as f:
        f.write("# H2 SCF via Z_eff Optimization\n")
        f.write(f"# Date: 2026-03-13\n")
        f.write(f"# Method: Eckart variational (optimize Z_eff in H2+ orbital)\n")
        f.write(f"#\n")
        f.write(f"# {'R':>8s} {'E_SCF_CI':>14s} {'E_frozen':>14s} "
                f"{'D_e_SCF':>14s} {'Z_eff':>10s} {'J_gg':>10s}\n")
        for r in results:
            R = r['R']
            E = r['E_total']
            Ef = r.get('E_HF_frozen', np.nan)
            z = r.get('Z_eff_opt', np.nan)
            J = r.get('J_gg', np.nan)
            f.write(f"  {R:8.4f} {E:14.8f} {Ef:14.8f} "
                    f"{2*E_H-E:14.8f} {z:10.4f} {J:10.6f}\n")
    print(f"\n  Results saved to {outpath}")


if __name__ == '__main__':
    main()
