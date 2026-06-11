"""
Phase 5: HeH+ Two-Electron CI on Prolate Spheroidal Lattice
============================================================
HeH+ is the simplest heteronuclear two-electron system (Z_A=2, Z_B=1).
Uses the same minimal-basis CI as H2 but with asymmetric nuclear charges.

Key differences from H2:
  - Orbitals are NOT gerade/ungerade (no inversion symmetry)
  - 1sigma orbital localizes on He at large R
  - Dissociation: HeH+ -> He(1s^2) + H+  [E_ref = -2.9037 Ha]

Date: 2026-03-13
"""

import warnings
warnings.filterwarnings('ignore')

import sys
import os
if sys.stdout.encoding and sys.stdout.encoding.lower().startswith('cp'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import numpy as np
from scipy.special import ellipk
from scipy.linalg import eigh
from math import factorial
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.prolate_spheroidal_lattice import (
    ProlateSpheroidalLattice,
    fit_spectroscopic_constants,
)
from debug.stress_test_prolate_h2 import compute_vee_integral


# ============================================================
# Generalized orbital generation with Z_A, Z_B support
# ============================================================
def get_orbital_on_grid_general(
    R: float,
    Z_A: int = 1,
    Z_B: int = 1,
    n_angular: int = 0,
    N_xi_solve: int = 5000,
    N_xi_grid: int = 80,
    N_eta_grid: int = 80,
    xi_max_grid: float = 12.0,
    m: int = 0,
) -> dict:
    """Solve for an orbital and evaluate it on a (xi,eta) quadrature grid.

    Generalized version supporting heteronuclear diatomics (Z_A != Z_B).
    """
    # --- Solve for energy ---
    lat = ProlateSpheroidalLattice(
        R=R, Z_A=Z_A, Z_B=Z_B, N_xi=N_xi_solve,
        xi_max=max(25.0, R * 3),
        m=m, n_angular=n_angular,
    )
    E_elec, c2, A = lat.solve()
    c = np.sqrt(max(c2, 1e-15))

    # --- Radial wavefunction F(xi) on solve grid ---
    N = N_xi_solve
    xi_min_s = 1.0 + 5e-4
    xi_max_s = max(25.0, R * 3)
    h_s = (xi_max_s - xi_min_s) / (N + 1)
    xi_s = xi_min_s + (np.arange(N) + 1) * h_s

    a_param = R * (Z_A + Z_B)
    xi2_1 = xi_s**2 - 1
    q = A + a_param * xi_s - c2 * xi_s**2
    if m != 0:
        q -= m**2 / xi2_1

    p_plus = (xi_s + h_s / 2)**2 - 1
    p_minus = (xi_s - h_s / 2)**2 - 1

    diag = -(p_plus + p_minus) / h_s**2 + q
    off = p_plus[:-1] / h_s**2

    if m == 0:
        diag[0] = -p_plus[0] / h_s**2 + q[0]

    from scipy.linalg import eigh_tridiagonal
    evals, evecs = eigh_tridiagonal(diag, off)
    idx_sorted = np.argsort(evals)[::-1]
    F_raw = evecs[:, idx_sorted[0]]
    if F_raw[0] < 0:
        F_raw = -F_raw

    # --- Angular eigenvector (Legendre basis) ---
    m_abs = abs(m)
    n_basis = 50
    r_vals = np.arange(m_abs, m_abs + n_basis, dtype=float)
    norms = np.array(
        [2.0 / (2*r+1) * factorial(int(r+m_abs)) / factorial(int(r-m_abs))
         for r in r_vals]
    )

    nu_mat = np.zeros((n_basis, n_basis))
    for i in range(n_basis - 1):
        r = r_vals[i]
        val = ((r - m_abs + 1) * np.sqrt(norms[i+1])
               / ((2*r+1) * np.sqrt(norms[i])))
        nu_mat[i+1, i] = val
        nu_mat[i, i+1] = val

    b_param = R * (Z_B - Z_A)
    H_ang = np.diag(-r_vals * (r_vals + 1))
    H_ang += c**2 * (nu_mat @ nu_mat)
    H_ang += b_param * nu_mat

    evals_ang, evecs_ang = np.linalg.eigh(H_ang)
    idx_ang = n_basis - 1 - n_angular
    coeffs_ang = evecs_ang[:, idx_ang]

    # --- Build quadrature grids ---
    eta_gl, w_eta_gl = np.polynomial.legendre.leggauss(N_eta_grid)
    eta = eta_gl
    w_eta = w_eta_gl

    u_gl, w_u_gl = np.polynomial.legendre.leggauss(N_xi_grid)
    t = (u_gl + 1) / 2
    xi = 1.0 + (xi_max_grid - 1.0) * t**2
    w_xi = w_u_gl * (xi_max_grid - 1.0) * t

    # --- Evaluate wavefunctions on quadrature grid ---
    F_xi = np.interp(xi, xi_s, F_raw, left=F_raw[0], right=0.0)

    from scipy.special import lpmv
    G_eta = np.zeros(N_eta_grid)
    for j, r in enumerate(r_vals):
        P_r = lpmv(m_abs, int(r), eta)
        G_eta += coeffs_ang[j] * P_r / np.sqrt(norms[j])

    # --- Build 2D wavefunction ---
    psi_2d = np.outer(F_xi, G_eta)

    XI, ETA = np.meshgrid(xi, eta, indexing='ij')
    J = (R / 2)**3 * (XI**2 - ETA**2)
    W_XI, W_ETA = np.meshgrid(w_xi, w_eta, indexing='ij')

    norm_sq = np.sum(psi_2d**2 * J * W_XI * W_ETA) * 2 * np.pi
    psi_2d /= np.sqrt(norm_sq)

    norm_check = np.sum(psi_2d**2 * J * W_XI * W_ETA) * 2 * np.pi
    assert abs(norm_check - 1.0) < 0.02, f"Normalization failed: {norm_check:.4f}"

    return {
        'E_elec': E_elec,
        'c2': c2,
        'A': A,
        'xi': xi,
        'eta': eta,
        'w_xi': w_xi,
        'w_eta': w_eta,
        'psi': psi_2d,
        'F_xi': F_xi,
        'G_eta': G_eta,
        'R': R,
        'Z_A': Z_A,
        'Z_B': Z_B,
        'n_angular': n_angular,
    }


# ============================================================
# HeH+ CI driver
# ============================================================
def heh_plus_ci(
    R: float,
    N_xi_solve: int = 5000,
    N_grid: int = 60,
    xi_max_grid: float = 15.0,
    verbose: bool = True,
) -> dict:
    """Minimal-basis (1sigma, 2sigma) CI for HeH+.

    Same structure as H2 CI but with Z_A=2, Z_B=1.
    No gerade/ungerade symmetry — both orbitals couple freely.
    """
    t0 = time.time()
    Z_A, Z_B = 2, 1

    if verbose:
        print(f"\n  HeH+ CI at R={R:.3f} bohr (Z_A={Z_A}, Z_B={Z_B})")
        print(f"  Grid: {N_grid}x{N_grid}, N_xi_solve={N_xi_solve}")

    # Step 1: Generate orbitals
    orbitals = {}
    for label, n_ang, name in [('1', 0, '1sigma'), ('2', 1, '2sigma')]:
        if verbose:
            print(f"  Solving for {name} (n_angular={n_ang})...")
        try:
            orb = get_orbital_on_grid_general(
                R=R, Z_A=Z_A, Z_B=Z_B, n_angular=n_ang,
                N_xi_solve=N_xi_solve,
                N_xi_grid=N_grid, N_eta_grid=N_grid,
                xi_max_grid=xi_max_grid,
            )
            orbitals[label] = orb
            if verbose:
                print(f"    E_elec({name}) = {orb['E_elec']:.6f} Ha")
        except Exception as e:
            if verbose:
                print(f"    FAILED: {e}")
            return {'R': R, 'E_total': np.nan, 'error': str(e)}

    # Step 2: V_ee integrals
    if verbose:
        print("  Computing V_ee integrals...")
    t1 = time.time()

    o1, o2 = orbitals['1'], orbitals['2']

    J_11 = compute_vee_integral(o1, o1, o1, o1)
    J_22 = compute_vee_integral(o2, o2, o2, o2)
    J_12 = compute_vee_integral(o1, o2, o1, o2)
    K_12 = compute_vee_integral(o1, o2, o2, o1)
    V_1122 = compute_vee_integral(o1, o1, o2, o2)  # double excitation coupling

    t_vee = time.time() - t1
    if verbose:
        print(f"    J_11 = {J_11:.6f}, J_22 = {J_22:.6f}")
        print(f"    J_12 = {J_12:.6f}, K_12 = {K_12:.6f}")
        print(f"    V_1122 = {V_1122:.6f}")
        print(f"    V_ee time: {t_vee:.1f}s")

    # Step 3: Build 2x2 CI matrix
    eps_1 = o1['E_elec']
    eps_2 = o2['E_elec']

    H_CI = np.array([
        [2*eps_1 + J_11, V_1122],
        [V_1122, 2*eps_2 + J_22],
    ])

    V_NN = Z_A * Z_B / R
    evals, evecs = eigh(H_CI)
    E_CI_elec = evals[0]
    E_total = E_CI_elec + V_NN
    E_HF = 2*eps_1 + J_11 + V_NN

    dt = time.time() - t0

    if verbose:
        c1 = evecs[0, 0]
        c2_coeff = evecs[1, 0]
        print(f"\n  CI matrix:")
        print(f"    H_11 = {H_CI[0,0]:.6f}, H_22 = {H_CI[1,1]:.6f}")
        print(f"    H_12 = {H_CI[0,1]:.6f}")
        print(f"  Eigenvalues: {evals}")
        print(f"  GS coefficients: c_1={c1:.4f}, c_2={c2_coeff:.4f}")
        print(f"\n  V_NN = {V_NN:.6f}")
        print(f"  E_HF    = {E_HF:.6f} Ha")
        print(f"  E_CI    = {E_total:.6f} Ha")
        print(f"  E_corr  = {E_total - E_HF:.6f} Ha")
        print(f"  Time: {dt:.1f}s")

    return {
        'R': R,
        'E_total': E_total,
        'E_HF': E_HF,
        'E_CI_elec': E_CI_elec,
        'V_NN': V_NN,
        'eps_1': eps_1,
        'eps_2': eps_2,
        'J_11': J_11,
        'J_22': J_22,
        'J_12': J_12,
        'K_12': K_12,
        'V_1122': V_1122,
        'c_1': evecs[0, 0],
        'c_2': evecs[1, 0],
        'time': dt,
    }


# ============================================================
# PES scan
# ============================================================
def scan_heh_pes(
    R_values: np.ndarray,
    N_xi_solve: int = 5000,
    N_grid: int = 50,
    xi_max_grid: float = 15.0,
    verbose: bool = True,
) -> dict:
    """Scan HeH+ PES."""
    results = []
    for R in R_values:
        try:
            res = heh_plus_ci(
                R=R, N_xi_solve=N_xi_solve, N_grid=N_grid,
                xi_max_grid=xi_max_grid, verbose=verbose,
            )
            results.append(res)
        except Exception as e:
            print(f"  R={R:.3f}: FAILED - {e}")
            results.append({'R': R, 'E_total': np.nan, 'E_HF': np.nan, 'time': 0})

    return {
        'R': np.array([r['R'] for r in results]),
        'E_total': np.array([r['E_total'] for r in results]),
        'E_HF': np.array([r.get('E_HF', np.nan) for r in results]),
        'results': results,
    }


# ============================================================
# Main
# ============================================================
def main():
    print("=" * 70)
    print("  HeH+ TWO-ELECTRON CI ON PROLATE SPHEROIDAL LATTICE")
    print("  Z_A=2 (He), Z_B=1 (H)  |  CI: 2x2 minimal basis")
    print("  Date: 2026-03-13")
    print("=" * 70)

    # Reference values
    E_He_exact = -2.9037  # He atom exact
    E_Hplus = 0.0         # H+ energy
    E_ref = E_He_exact + E_Hplus  # Dissociation limit
    R_eq_exact = 1.46     # bohr (approximate)
    D_e_exact = 0.075     # Ha (approximate)

    # --- Test 1: Single point near R_eq ---
    print("\n" + "="*60)
    print("  TEST 1: Single point at R=1.5 bohr")
    print("="*60)
    res = heh_plus_ci(R=1.5, N_xi_solve=5000, N_grid=50, verbose=True)

    if not np.isnan(res['E_total']):
        D_e = E_ref - res['E_total']
        print(f"\n  D_e = E_ref - E = {E_ref:.4f} - ({res['E_total']:.4f}) = {D_e:.4f} Ha")
        print(f"  Bound: {res['E_total'] < E_ref}")

    # --- Test 2: Grid convergence ---
    print("\n" + "="*60)
    print("  TEST 2: Grid convergence at R=1.5")
    print("="*60)
    for N_grid in [30, 40, 50, 60]:
        rc = heh_plus_ci(R=1.5, N_xi_solve=5000, N_grid=N_grid, verbose=False)
        if not np.isnan(rc['E_total']):
            De = E_ref - rc['E_total']
            print(f"  N_grid={N_grid:3d}: E={rc['E_total']:.6f}  D_e={De:.6f}  [{rc['time']:.1f}s]")

    # --- Test 3: PES scan ---
    print("\n" + "="*60)
    print("  TEST 3: HeH+ PES scan")
    print("="*60)
    R_vals = np.array([0.8, 1.0, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0])
    pes = scan_heh_pes(R_vals, N_xi_solve=5000, N_grid=50, xi_max_grid=15.0, verbose=False)

    print(f"\n  {'R':>6s} {'E_CI':>12s} {'E_HF':>12s} {'D_e_CI':>10s} {'D_e_HF':>10s}")
    print("-" * 60)
    for i, R in enumerate(pes['R']):
        Eci = pes['E_total'][i]
        Ehf = pes['E_HF'][i]
        De_ci = E_ref - Eci
        De_hf = E_ref - Ehf
        print(f"  {R:6.3f} {Eci:12.6f} {Ehf:12.6f} {De_ci:10.6f} {De_hf:10.6f}")

    # Fit spectroscopic constants
    valid = ~np.isnan(pes['E_total'])
    if np.sum(valid) >= 5:
        # Use He+H+ dissociation limit for D_e
        fit_ci = fit_spectroscopic_constants(pes['R'][valid], pes['E_total'][valid])
        fit_hf = fit_spectroscopic_constants(pes['R'][valid], pes['E_HF'][valid])
        # Override D_e to use proper reference
        fit_ci['D_e'] = E_ref - fit_ci['E_min']
        fit_hf['D_e'] = E_ref - fit_hf['E_min']

        print(f"\n  Spectroscopic constants:")
        print(f"  {'':20s} {'CI':>12s} {'HF':>12s} {'Exact':>12s}")
        print(f"  {'R_eq (bohr)':20s} {fit_ci['R_eq']:12.4f} {fit_hf['R_eq']:12.4f} {R_eq_exact:12.4f}")
        print(f"  {'E_min (Ha)':20s} {fit_ci['E_min']:12.6f} {fit_hf['E_min']:12.6f} {'':>12s}")
        print(f"  {'D_e (Ha)':20s} {fit_ci['D_e']:12.6f} {fit_hf['D_e']:12.6f} {D_e_exact:12.6f}")

        R_err = abs(fit_ci['R_eq'] - R_eq_exact) / R_eq_exact * 100
        D_err = abs(fit_ci['D_e'] - D_e_exact) / D_e_exact * 100

        print(f"\n  CI errors:")
        print(f"    R_eq: {R_err:.1f}%")
        print(f"    D_e:  {D_err:.1f}%")

        print(f"\n  {'='*55}")
        if fit_ci['D_e'] > 0:
            print(f"  HeH+ IS BOUND (D_e = {fit_ci['D_e']:.4f} Ha)")
            if R_err < 30:
                print(f"  Error pattern {'MATCHES' if R_err > 15 else 'BETTER THAN'} H2 (27% R_eq error)")
        else:
            print(f"  HeH+ NOT BOUND in our framework")
        print(f"  {'='*55}")

    # Save results
    outpath = os.path.join(os.path.dirname(__file__), 'data', 'prolate_heh_ci.txt')
    with open(outpath, 'w', encoding='utf-8') as f:
        f.write("# HeH+ Minimal-Basis CI on Prolate Spheroidal Lattice\n")
        f.write(f"# Date: 2026-03-13\n")
        f.write(f"# Z_A=2 (He), Z_B=1 (H), N_xi_solve=5000, N_grid=50\n")
        f.write(f"# Dissociation limit: He + H+ = {E_ref:.6f} Ha\n")
        f.write(f"#\n")
        f.write(f"# {'R':>8s} {'E_CI':>14s} {'E_HF':>14s} {'D_e_CI':>14s} {'D_e_HF':>14s}\n")
        for i, R in enumerate(pes['R']):
            Eci = pes['E_total'][i]
            Ehf = pes['E_HF'][i]
            f.write(f"  {R:8.4f} {Eci:14.8f} {Ehf:14.8f} "
                    f"{E_ref-Eci:14.8f} {E_ref-Ehf:14.8f}\n")
    print(f"\n  Results saved to {outpath}")


if __name__ == '__main__':
    main()
