"""
Heteronuclear SCF for prolate spheroidal diatomics.

Implements per-atom Z_eff optimization for heteronuclear systems like HeH+,
where independent orbital exponents Z_eff_A and Z_eff_B allow the orbital
to redistribute density between the two centers.

For HeH+ with frozen HeH2+ orbitals, J_11 ~ 1.3 Ha makes the system
unbound. Independent Z_eff optimization should:
  - Screen Z_He: Z_eff_A < 2.0
  - Attract toward H: Z_eff_B > 1.0 (partial delocalization)
  - Reduce J_11 enough to bind

The effective nuclear charges modify the one-electron Hamiltonian:
  V_nuc^eff = -Z_eff_A / r_A - Z_eff_B / r_B
"""

import numpy as np
from scipy.optimize import minimize
from scipy.linalg import eigh
from typing import Dict
import time

from geovac.prolate_scf import (
    get_orbital_on_grid,
    compute_vee_integral,
    compute_vnuc_expectation,
)
from geovac.prolate_spheroidal_lattice import fit_spectroscopic_constants


# ============================================================
# Heteronuclear one-electron energy at given Z_eff_A, Z_eff_B
# ============================================================

def heteronuclear_scf_energy(
    R: float,
    Z_A_phys: float,
    Z_B_phys: float,
    Z_eff_A: float,
    Z_eff_B: float,
    N_xi_solve: int = 5000,
    N_grid: int = 50,
    xi_max_grid: float = 15.0,
) -> Dict:
    """Compute HF energy with per-atom effective charges.

    Generates the 1sigma orbital from the (Z_eff_A, Z_eff_B, R) problem,
    then evaluates the total energy with physical charges (Z_A, Z_B).

    The one-electron energy with physical charges:
      h_phys = <phi|T|phi> + <phi|-Z_A/r_A - Z_B/r_B|phi>

    We compute this via the identity:
      <phi|T + V_nuc(Z_eff)|phi> = eps(Z_eff)  [eigenvalue]
      <phi|V_nuc(Z)|phi> = -(Z_A/Z_eff_A + Z_B/Z_eff_B) * <phi|V_nuc(Z_eff)|phi>
      ... which is not quite right for heteronuclear.

    Instead, compute h_phys directly on the grid.

    Parameters
    ----------
    R : float
        Internuclear distance.
    Z_A_phys, Z_B_phys : float
        Physical nuclear charges.
    Z_eff_A, Z_eff_B : float
        Effective charges for orbital generation.
    """
    try:
        orb = get_orbital_on_grid(
            R=R, Z_A=Z_eff_A, Z_B=Z_eff_B, n_angular=0,
            N_xi_solve=N_xi_solve,
            N_xi_grid=N_grid, N_eta_grid=N_grid,
            xi_max_grid=xi_max_grid,
        )
    except Exception as e:
        return {'R': R, 'E_HF': np.nan, 'error': str(e)}

    eps_eff = orb['E_elec']

    # Compute <1/r_A> and <1/r_B> separately
    xi = orb['xi']
    eta = orb['eta']
    w_xi = orb['w_xi']
    w_eta = orb['w_eta']
    psi = orb['psi']

    XI = xi[:, np.newaxis]
    ETA = eta[np.newaxis, :]
    W_XI = w_xi[:, np.newaxis]
    W_ETA = w_eta[np.newaxis, :]

    # r_A = R*(xi+eta)/2, r_B = R*(xi-eta)/2
    # 1/r_A * J = (R/2)^3 (xi^2-eta^2) * 2/(R(xi+eta)) = (R/2)^2 (xi-eta)
    # 1/r_B * J = (R/2)^3 (xi^2-eta^2) * 2/(R(xi-eta)) = (R/2)^2 (xi+eta)
    inv_rA_J = (R / 2) ** 2 * (XI - ETA)
    inv_rB_J = (R / 2) ** 2 * (XI + ETA)

    expect_inv_rA = 2 * np.pi * np.sum(psi ** 2 * inv_rA_J * W_XI * W_ETA)
    expect_inv_rB = 2 * np.pi * np.sum(psi ** 2 * inv_rB_J * W_XI * W_ETA)

    # One-electron energy with physical charges:
    # h_phys = <T> - Z_A_phys * <1/r_A> - Z_B_phys * <1/r_B>
    # <T> = eps_eff - <V_nuc(Z_eff)> = eps_eff + Z_eff_A*<1/r_A> + Z_eff_B*<1/r_B>
    T_expect = eps_eff + Z_eff_A * expect_inv_rA + Z_eff_B * expect_inv_rB
    h_phys = T_expect - Z_A_phys * expect_inv_rA - Z_B_phys * expect_inv_rB

    # J_gg
    J_gg = compute_vee_integral(orb, orb, orb, orb)

    # Total HF energy
    V_NN = Z_A_phys * Z_B_phys / R
    E_HF = 2 * h_phys + J_gg + V_NN

    return {
        'R': R,
        'Z_eff_A': Z_eff_A,
        'Z_eff_B': Z_eff_B,
        'E_HF': E_HF,
        'h_gg': h_phys,
        'J_gg': J_gg,
        'T_expect': T_expect,
        'inv_rA': expect_inv_rA,
        'inv_rB': expect_inv_rB,
        'eps_eff': eps_eff,
        'V_NN': V_NN,
    }


# ============================================================
# Per-atom Z_eff optimization
# ============================================================

def optimize_heteronuclear_zeff(
    R: float,
    Z_A_phys: float = 2.0,
    Z_B_phys: float = 1.0,
    N_xi_solve: int = 5000,
    N_grid: int = 50,
    xi_max_grid: float = 15.0,
    verbose: bool = True,
) -> Dict:
    """Optimize per-atom Z_eff to minimize E_HF for heteronuclear systems.

    Uses Nelder-Mead optimization over (Z_eff_A, Z_eff_B) with physical
    constraints: 0.5 <= Z_eff <= Z_bare + 0.5.

    Parameters
    ----------
    R : float
        Internuclear distance.
    Z_A_phys, Z_B_phys : float
        Physical nuclear charges (e.g., 2, 1 for HeH+).
    """
    t0 = time.time()

    if verbose:
        print(f"\n  Optimizing Z_eff at R={R:.3f} (Z_A={Z_A_phys}, Z_B={Z_B_phys})")

    def objective(params: np.ndarray) -> float:
        z_a, z_b = params
        # Bounds enforcement via penalty
        if z_a < 0.3 or z_a > Z_A_phys + 1.0 or z_b < 0.3 or z_b > Z_B_phys + 1.0:
            return 1e10
        res = heteronuclear_scf_energy(
            R, Z_A_phys, Z_B_phys, z_a, z_b,
            N_xi_solve, N_grid, xi_max_grid,
        )
        return res['E_HF'] if not np.isnan(res['E_HF']) else 1e10

    # Initial guess: physical charges
    x0 = np.array([Z_A_phys, Z_B_phys])

    # Coarse grid search first
    best_E = 1e10
    best_z = x0.copy()
    z_A_vals = np.linspace(max(0.5, Z_A_phys - 0.8), Z_A_phys + 0.3, 7)
    z_B_vals = np.linspace(max(0.5, Z_B_phys - 0.3), Z_B_phys + 0.5, 7)

    if verbose:
        print(f"  Coarse grid: Z_A in [{z_A_vals[0]:.1f}, {z_A_vals[-1]:.1f}], "
              f"Z_B in [{z_B_vals[0]:.1f}, {z_B_vals[-1]:.1f}]")

    for za in z_A_vals:
        for zb in z_B_vals:
            E = objective(np.array([za, zb]))
            if E < best_E:
                best_E = E
                best_z = np.array([za, zb])

    if verbose:
        print(f"  Coarse best: Z_eff_A={best_z[0]:.3f}, Z_eff_B={best_z[1]:.3f}, "
              f"E={best_E:.6f}")

    # Refine with Nelder-Mead
    result = minimize(
        objective, best_z, method='Nelder-Mead',
        options={'xatol': 0.005, 'fatol': 1e-5, 'maxiter': 100},
    )

    z_opt = result.x
    res_opt = heteronuclear_scf_energy(
        R, Z_A_phys, Z_B_phys, z_opt[0], z_opt[1],
        N_xi_solve, N_grid, xi_max_grid,
    )

    # Frozen result (physical charges)
    res_frozen = heteronuclear_scf_energy(
        R, Z_A_phys, Z_B_phys, Z_A_phys, Z_B_phys,
        N_xi_solve, N_grid, xi_max_grid,
    )

    dt = time.time() - t0

    if verbose:
        print(f"\n  Optimized: Z_eff_A={z_opt[0]:.4f}, Z_eff_B={z_opt[1]:.4f}")
        print(f"  E_HF(opt)    = {res_opt['E_HF']:.6f} Ha")
        print(f"  E_HF(frozen) = {res_frozen['E_HF']:.6f} Ha")
        print(f"  Improvement:   {res_frozen['E_HF'] - res_opt['E_HF']:.6f} Ha")
        print(f"  J_gg(opt) = {res_opt['J_gg']:.6f}, J_gg(frozen) = {res_frozen['J_gg']:.6f}")
        print(f"  Time: {dt:.1f}s")

    return {
        'R': R,
        'Z_eff_A': z_opt[0],
        'Z_eff_B': z_opt[1],
        'E_HF_opt': res_opt['E_HF'],
        'E_HF_frozen': res_frozen['E_HF'],
        'J_gg_opt': res_opt['J_gg'],
        'J_gg_frozen': res_frozen['J_gg'],
        'h_gg_opt': res_opt['h_gg'],
        'h_gg_frozen': res_frozen['h_gg'],
        'inv_rA_opt': res_opt['inv_rA'],
        'inv_rB_opt': res_opt['inv_rB'],
        'time': dt,
    }


# ============================================================
# Heteronuclear SCF + CI
# ============================================================

def heteronuclear_scf_ci(
    R: float,
    Z_A_phys: float = 2.0,
    Z_B_phys: float = 1.0,
    N_xi_solve: int = 5000,
    N_grid: int = 50,
    xi_max_grid: float = 15.0,
    verbose: bool = True,
) -> Dict:
    """SCF + CI for heteronuclear diatomic (e.g., HeH+).

    1. Optimize (Z_eff_A, Z_eff_B) for the 1sigma orbital
    2. Generate 1sigma and 2sigma orbitals at optimized Z_eff
    3. Run 2x2 CI for correlation energy

    Returns
    -------
    dict with E_total, E_HF_scf, E_HF_frozen, Z_eff_A/B, etc.
    """
    t0 = time.time()
    V_NN = Z_A_phys * Z_B_phys / R

    # Step 1: Optimize Z_eff
    opt = optimize_heteronuclear_zeff(
        R, Z_A_phys, Z_B_phys, N_xi_solve, N_grid, xi_max_grid,
        verbose=verbose,
    )

    z_a = opt['Z_eff_A']
    z_b = opt['Z_eff_B']

    # Step 2: Generate both orbitals
    orbitals = {}
    for label, n_ang in [('1', 0), ('2', 1)]:
        try:
            orb = get_orbital_on_grid(
                R=R, Z_A=z_a, Z_B=z_b, n_angular=n_ang,
                N_xi_solve=N_xi_solve,
                N_xi_grid=N_grid, N_eta_grid=N_grid,
                xi_max_grid=xi_max_grid,
            )
            orbitals[label] = orb
        except Exception as e:
            if verbose:
                print(f"  Failed to generate orbital {label}: {e}")
            return {'R': R, 'E_total': np.nan, 'error': str(e)}

    # Step 3: One-electron energies with physical charges
    h_vals = {}
    for label in ['1', '2']:
        orb = orbitals[label]
        xi = orb['xi']
        eta = orb['eta']
        w_xi = orb['w_xi']
        w_eta = orb['w_eta']
        psi = orb['psi']

        XI = xi[:, np.newaxis]
        ETA = eta[np.newaxis, :]
        W_XI = w_xi[:, np.newaxis]
        W_ETA = w_eta[np.newaxis, :]

        inv_rA_J = (R / 2) ** 2 * (XI - ETA)
        inv_rB_J = (R / 2) ** 2 * (XI + ETA)

        inv_rA = 2 * np.pi * np.sum(psi ** 2 * inv_rA_J * W_XI * W_ETA)
        inv_rB = 2 * np.pi * np.sum(psi ** 2 * inv_rB_J * W_XI * W_ETA)

        eps_eff = orb['E_elec']
        T_exp = eps_eff + z_a * inv_rA + z_b * inv_rB
        h_phys = T_exp - Z_A_phys * inv_rA - Z_B_phys * inv_rB
        h_vals[label] = h_phys

    # Step 4: V_ee integrals
    o1, o2 = orbitals['1'], orbitals['2']
    J_11 = compute_vee_integral(o1, o1, o1, o1)
    J_22 = compute_vee_integral(o2, o2, o2, o2)
    V_1122 = compute_vee_integral(o1, o1, o2, o2)

    if verbose:
        print(f"  J_11={J_11:.6f}, J_22={J_22:.6f}, V_1122={V_1122:.6f}")

    # Step 5: Build CI matrix
    H_CI = np.array([
        [2 * h_vals['1'] + J_11, V_1122],
        [V_1122, 2 * h_vals['2'] + J_22],
    ])

    evals, evecs = eigh(H_CI)
    E_CI_elec = evals[0]
    E_total = E_CI_elec + V_NN
    E_HF = 2 * h_vals['1'] + J_11 + V_NN

    dt = time.time() - t0

    if verbose:
        print(f"\n  h_1={h_vals['1']:.6f}, h_2={h_vals['2']:.6f}")
        print(f"  E_HF(SCF) = {E_HF:.6f}, E_CI(SCF) = {E_total:.6f}")
        print(f"  CI coeffs: c_1={evecs[0, 0]:.4f}, c_2={evecs[1, 0]:.4f}")
        print(f"  Time: {dt:.1f}s")

    return {
        'R': R,
        'E_total': E_total,
        'E_HF_scf': E_HF,
        'E_HF_frozen': opt['E_HF_frozen'],
        'Z_eff_A': z_a,
        'Z_eff_B': z_b,
        'J_11': J_11,
        'J_22': J_22,
        'V_1122': V_1122,
        'h_1': h_vals['1'],
        'h_2': h_vals['2'],
        'c_1': evecs[0, 0],
        'c_2': evecs[1, 0],
        'time': dt,
    }
