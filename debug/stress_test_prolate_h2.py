"""
Phase 2: H2 Minimal-Basis CI on Prolate Spheroidal Lattice
============================================================
Build a 2-electron CI using 1sigma_g and 1sigma_u orbitals from the
one-electron H2+ solver. Compute V_ee integrals via azimuthal averaging
(complete elliptic integral K) + numerical quadrature on (xi,eta) grid.

The critical test: can we get bound H2 with correct R_eq and D_e?

Date: 2026-03-13
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
from scipy.special import ellipk
from scipy.linalg import eigh
import time
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.prolate_spheroidal_lattice import (
    ProlateSpheroidalLattice,
    fit_spectroscopic_constants,
)


# ============================================================
# Wavefunction extraction on a common (xi, eta) grid
# ============================================================
def get_orbital_on_grid(
    R: float,
    n_angular: int,
    N_xi_solve: int = 5000,
    N_xi_grid: int = 80,
    N_eta_grid: int = 80,
    xi_max_grid: float = 12.0,
    m: int = 0,
) -> dict:
    """Solve for an orbital and evaluate it on a (xi,eta) grid.

    Uses Gauss-Legendre quadrature for eta and a sinh-mapped grid for xi
    (concentrating points near xi=1 where wavefunctions peak).

    Returns dict with E_elec, c2, A, xi, eta, w_xi, w_eta, psi(xi,eta).
    """
    # --- Solve for energy ---
    lat = ProlateSpheroidalLattice(
        R=R, N_xi=N_xi_solve, xi_max=25.0,
        m=m, n_angular=n_angular,
    )
    E_elec, c2, A = lat.solve()
    c = np.sqrt(max(c2, 1e-15))

    # --- Radial wavefunction F(xi) on solve grid (for eigenvector) ---
    N = N_xi_solve
    xi_min_s = 1.0 + 5e-4
    h_s = (25.0 - xi_min_s) / (N + 1)
    xi_s = xi_min_s + (np.arange(N) + 1) * h_s

    a_param = R * (1 + 1)  # Z_A + Z_B for H2+
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
    F_raw = evecs[:, idx_sorted[0]]  # Top eigenvalue
    if F_raw[0] < 0:
        F_raw = -F_raw

    # --- Angular eigenvector (Legendre basis) ---
    from math import factorial
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

    H_ang = np.diag(-r_vals * (r_vals + 1))
    H_ang += c**2 * (nu_mat @ nu_mat)

    evals_ang, evecs_ang = np.linalg.eigh(H_ang)
    idx_ang = n_basis - 1 - n_angular
    coeffs_ang = evecs_ang[:, idx_ang]

    # --- Build quadrature grids ---

    # eta: Gauss-Legendre on [-1, 1]
    eta_gl, w_eta_gl = np.polynomial.legendre.leggauss(N_eta_grid)
    eta = eta_gl
    w_eta = w_eta_gl

    # xi: sinh-mapped grid to cluster points near xi=1
    # Map u in [-1,1] (Gauss-Legendre) to xi via xi = 1 + alpha*sinh(beta*u)
    # Choose alpha, beta so that xi goes from ~1 to xi_max
    u_gl, w_u_gl = np.polynomial.legendre.leggauss(N_xi_grid)
    # Linear map: u in [-1,1] -> t in [0, t_max]
    # Then xi = 1 + t^2/t_max * (xi_max - 1) concentrates near xi=1
    # Actually, simpler: use Gauss-Laguerre-like mapping
    # xi = 1 + (xi_max-1) * ((u+1)/2)^2  (quadratic concentration)
    t = (u_gl + 1) / 2  # t in [0, 1]
    xi = 1.0 + (xi_max_grid - 1.0) * t**2
    dxi_du = (xi_max_grid - 1.0) * 2 * t * 0.5  # d(xi)/d(u) * d(u_raw)
    # Weight: w_xi = w_u * dxi/du where dxi/du = (xi_max-1)*2*t * (1/2)
    w_xi = w_u_gl * (xi_max_grid - 1.0) * t

    # --- Evaluate wavefunctions on quadrature grid ---
    # Interpolate F(xi) from solve grid
    F_xi = np.interp(xi, xi_s, F_raw, left=F_raw[0], right=0.0)

    # Evaluate G(eta) using Legendre expansion
    from scipy.special import lpmv
    G_eta = np.zeros(N_eta_grid)
    for j, r in enumerate(r_vals):
        P_r = lpmv(m_abs, int(r), eta)
        G_eta += coeffs_ang[j] * P_r / np.sqrt(norms[j])

    # --- Build 2D wavefunction psi(xi, eta) = F(xi)*G(eta) ---
    psi_2d = np.outer(F_xi, G_eta)

    # --- Normalize using quadrature weights ---
    # int |psi|^2 dV = sum_ij |psi_ij|^2 * J_ij * w_xi_i * w_eta_j * 2*pi
    XI, ETA = np.meshgrid(xi, eta, indexing='ij')
    J = (R / 2)**3 * (XI**2 - ETA**2)

    W_XI, W_ETA = np.meshgrid(w_xi, w_eta, indexing='ij')

    norm_sq = np.sum(psi_2d**2 * J * W_XI * W_ETA) * 2 * np.pi
    psi_2d /= np.sqrt(norm_sq)

    # Verify normalization
    norm_check = np.sum(psi_2d**2 * J * W_XI * W_ETA) * 2 * np.pi
    assert abs(norm_check - 1.0) < 0.01, f"Normalization failed: {norm_check:.4f}"

    return {
        'E_elec': E_elec,
        'c2': c2,
        'A': A,
        'xi': xi,
        'eta': eta,
        'w_xi': w_xi,
        'w_eta': w_eta,
        'psi': psi_2d,  # shape (N_xi, N_eta)
        'F_xi': F_xi,
        'G_eta': G_eta,
        'R': R,
        'n_angular': n_angular,
    }


# ============================================================
# V_ee two-electron integrals via azimuthal averaging
# ============================================================
def compute_vee_integral(
    orb_a: dict, orb_b: dict, orb_c: dict, orb_d: dict,
) -> float:
    """Compute (ab|cd) = int psi_a(1)*psi_c(1) * (1/r12) * psi_b(2)*psi_d(2) dV1 dV2.

    Uses azimuthal averaging: for m=0 orbitals, integrate out phi1-phi2
    analytically using complete elliptic integral K.

    Convention: chemists' notation (ab|cd) = [ac|bd] in physicists' notation.

    The double-phi integral for m=0 orbitals:
      int_0^{2pi} dphi1 int_0^{2pi} dphi2 / r12
      = 2pi * int_0^{2pi} 1/r12(dphi) d(dphi)
      = 2pi * 4*K(k) / sqrt(s) = 8pi * K(k) / sqrt(s)

    where s = (rho1+rho2)^2 + dz^2, k^2 = 4*rho1*rho2/s.
    """
    R = orb_a['R']
    xi = orb_a['xi']
    eta = orb_a['eta']
    w_xi = orb_a['w_xi']
    w_eta = orb_a['w_eta']

    N_xi = len(xi)
    N_eta = len(eta)

    # Density products on grid
    rho_ac = orb_a['psi'] * orb_c['psi']  # (N_xi, N_eta)
    rho_bd = orb_b['psi'] * orb_d['psi']  # (N_xi, N_eta)

    # Coordinate arrays
    XI, ETA = np.meshgrid(xi, eta, indexing='ij')

    # Cylindrical coordinates
    rho_cyl = (R / 2) * np.sqrt(np.maximum((XI**2 - 1) * (1 - ETA**2), 0.0))
    z_cyl = (R / 2) * XI * ETA

    # Jacobian: (R/2)^3 * (xi^2 - eta^2)
    J = (R / 2)**3 * (XI**2 - ETA**2)

    # Quadrature weights (2D product)
    W_XI, W_ETA = np.meshgrid(w_xi, w_eta, indexing='ij')
    W_2d = W_XI * W_ETA  # (N_xi, N_eta)

    # Flatten for pairwise computation
    n_pts = N_xi * N_eta
    # Weighted density = psi_a * psi_c * J * quadrature_weight
    rho_ac_flat = (rho_ac * J * W_2d).flatten()
    rho_bd_flat = (rho_bd * J * W_2d).flatten()
    rho_flat = rho_cyl.flatten()
    z_flat = z_cyl.flatten()

    # O(n_pts^2) computation in batches
    result = 0.0
    batch_size = 500
    n_batches = (n_pts + batch_size - 1) // batch_size

    for i_batch in range(n_batches):
        i_start = i_batch * batch_size
        i_end = min(i_start + batch_size, n_pts)

        rho1 = rho_flat[i_start:i_end, np.newaxis]
        z1 = z_flat[i_start:i_end, np.newaxis]

        rho2 = rho_flat[np.newaxis, :]
        z2 = z_flat[np.newaxis, :]

        dz = z1 - z2
        s = (rho1 + rho2)**2 + dz**2

        with np.errstate(divide='ignore', invalid='ignore'):
            k2 = np.where(s > 1e-30, 4 * rho1 * rho2 / s, 0.0)
            k2 = np.clip(k2, 0, 1 - 1e-15)

        s_safe = np.maximum(s, 1e-30)

        K_vals = ellipk(k2)
        coulomb_kernel = 8 * np.pi * K_vals / np.sqrt(s_safe)
        coulomb_kernel = np.where(s < 1e-30, 0.0, coulomb_kernel)

        w_ac = rho_ac_flat[i_start:i_end]
        w_bd = rho_bd_flat

        result += np.sum(w_ac[:, np.newaxis] * coulomb_kernel * w_bd[np.newaxis, :])

    return result


def compute_h1_integral(
    orb_a: dict, orb_b: dict, R: float,
) -> float:
    """Compute <a|h1|b> = <a|-nabla^2/2 - Z_A/r_A - Z_B/r_B|b>.

    For orbitals that are eigenstates of the one-electron Hamiltonian,
    <a|h1|b> = E_a * delta_{ab}. We use this for efficiency.
    """
    if orb_a['n_angular'] == orb_b['n_angular']:
        return orb_a['E_elec']
    else:
        # Off-diagonal: should be zero for eigenstates of same Hamiltonian
        return 0.0


# ============================================================
# H2 Minimal CI
# ============================================================
def h2_minimal_ci(
    R: float,
    N_xi_solve: int = 5000,
    N_grid: int = 60,
    xi_max_grid: float = 15.0,
    verbose: bool = True,
) -> dict:
    """Build minimal-basis (1sigma_g, 1sigma_u) CI for H2.

    The CI basis:
    |1> = |g_up, g_dn>  (both in sigma_g)
    |2> = |u_up, u_dn>  (both in sigma_u)
    |3> = (|g_up, u_dn> - |u_up, g_dn>) / sqrt(2)  (singlet open-shell)

    In the singlet sector, the CI matrix is 2x2:
    H = [[E_gg + J_gg, K_gu],
         [K_gu, E_uu + J_uu]]

    where:
    E_gg = 2*eps_g (one-electron energies)
    E_uu = 2*eps_u
    J_gg = (gg|gg) (Coulomb self-interaction)
    J_uu = (uu|uu)
    K_gu = (gu|ug) (exchange, equals (gg|uu) for real orbitals)

    The open-shell singlet |3> has energy E_gu = eps_g + eps_u + J_gu - K_gu
    but only couples to |1> and |2> through (gg|uu).

    Full 2x2 CI in {|gg>, |uu>} basis:
    H_11 = 2*eps_g + (gg|gg)
    H_22 = 2*eps_u + (uu|uu)
    H_12 = (gg|uu)  [double excitation coupling]
    """
    t0 = time.time()

    if verbose:
        print(f"\n  H2 CI at R={R:.3f} bohr")
        print(f"  Grid: {N_grid}x{N_grid}, N_xi_solve={N_xi_solve}")

    # Step 1: Get 1sigma_g and 1sigma_u orbitals
    if verbose:
        print("  Solving for 1sigma_g...")
    orb_g = get_orbital_on_grid(
        R=R, n_angular=0, N_xi_solve=N_xi_solve,
        N_xi_grid=N_grid, N_eta_grid=N_grid,
        xi_max_grid=xi_max_grid,
    )

    if verbose:
        print(f"    E_elec(g) = {orb_g['E_elec']:.6f} Ha")
        print("  Solving for 1sigma_u...")
    orb_u = get_orbital_on_grid(
        R=R, n_angular=1, N_xi_solve=N_xi_solve,
        N_xi_grid=N_grid, N_eta_grid=N_grid,
        xi_max_grid=xi_max_grid,
    )

    if verbose:
        print(f"    E_elec(u) = {orb_u['E_elec']:.6f} Ha")

    # Step 2: Compute two-electron integrals
    if verbose:
        print("  Computing V_ee integrals...")

    t1 = time.time()

    # (gg|gg): both electrons in sigma_g
    J_gg = compute_vee_integral(orb_g, orb_g, orb_g, orb_g)
    if verbose:
        print(f"    (gg|gg) = {J_gg:.6f} Ha")

    # (uu|uu): both electrons in sigma_u
    J_uu = compute_vee_integral(orb_u, orb_u, orb_u, orb_u)
    if verbose:
        print(f"    (uu|uu) = {J_uu:.6f} Ha")

    # (gu|gu): Coulomb integral between g and u
    J_gu = compute_vee_integral(orb_g, orb_u, orb_g, orb_u)
    if verbose:
        print(f"    (gu|gu) = {J_gu:.6f} Ha")

    # (gu|ug): Exchange integral (= (gg|uu) for real orbitals)
    K_gu = compute_vee_integral(orb_g, orb_u, orb_u, orb_g)
    if verbose:
        print(f"    (gu|ug) = {K_gu:.6f} Ha")

    # (gg|uu): Double excitation coupling
    V_gguu = compute_vee_integral(orb_g, orb_g, orb_u, orb_u)
    if verbose:
        print(f"    (gg|uu) = {V_gguu:.6f} Ha")

    t_vee = time.time() - t1
    if verbose:
        print(f"    V_ee integrals: {t_vee:.1f}s")

    # Step 3: Build CI matrix in {|gg>, |uu>} basis
    eps_g = orb_g['E_elec']
    eps_u = orb_u['E_elec']

    H_CI = np.array([
        [2*eps_g + J_gg, V_gguu],
        [V_gguu, 2*eps_u + J_uu],
    ])

    if verbose:
        print(f"\n  CI matrix:")
        print(f"    H_11 = 2*eps_g + J_gg = {2*eps_g:.4f} + {J_gg:.4f} = {H_CI[0,0]:.6f}")
        print(f"    H_22 = 2*eps_u + J_uu = {2*eps_u:.4f} + {J_uu:.4f} = {H_CI[1,1]:.6f}")
        print(f"    H_12 = (gg|uu) = {V_gguu:.6f}")

    # Step 4: Diagonalize
    evals, evecs = eigh(H_CI)
    E_CI = evals[0]  # Ground state electronic energy

    # Add nuclear repulsion
    V_NN = 1.0 / R  # Z_A * Z_B / R = 1/R for H2
    E_total = E_CI + V_NN

    if verbose:
        c_g = evecs[0, 0]
        c_u = evecs[1, 0]
        print(f"\n  CI eigenvalues: {evals}")
        print(f"  Ground state: E_CI = {E_CI:.6f} Ha")
        print(f"    Coefficients: c_g = {c_g:.4f}, c_u = {c_u:.4f}")
        print(f"    Weight: |c_g|^2 = {c_g**2:.4f}, |c_u|^2 = {c_u**2:.4f}")
        print(f"  V_NN = {V_NN:.6f} Ha")
        print(f"  E_total = {E_total:.6f} Ha")

    # Also compute Hartree-Fock (no CI mixing)
    E_HF = 2*eps_g + J_gg + V_NN

    # Open-shell singlet energy (for reference)
    E_os = eps_g + eps_u + J_gu - K_gu + V_NN

    dt = time.time() - t0
    if verbose:
        print(f"\n  E_HF = {E_HF:.6f} Ha")
        print(f"  E_CI = {E_total:.6f} Ha")
        print(f"  Correlation energy = {E_total - E_HF:.6f} Ha")
        print(f"  Total time: {dt:.1f}s")

    return {
        'R': R,
        'E_total': E_total,
        'E_CI_elec': E_CI,
        'E_HF': E_HF,
        'E_open_shell': E_os,
        'V_NN': V_NN,
        'eps_g': eps_g,
        'eps_u': eps_u,
        'J_gg': J_gg,
        'J_uu': J_uu,
        'J_gu': J_gu,
        'K_gu': K_gu,
        'V_gguu': V_gguu,
        'c_g': evecs[0, 0],
        'c_u': evecs[1, 0],
        'time': dt,
    }


# ============================================================
# H2 PES scan
# ============================================================
def scan_h2_pes(
    R_values: np.ndarray,
    N_xi_solve: int = 5000,
    N_grid: int = 50,
    xi_max_grid: float = 15.0,
    verbose: bool = True,
) -> dict:
    """Scan H2 PES with minimal CI."""
    results = []

    for R in R_values:
        try:
            res = h2_minimal_ci(
                R=R, N_xi_solve=N_xi_solve, N_grid=N_grid,
                xi_max_grid=xi_max_grid, verbose=verbose,
            )
            results.append(res)
        except Exception as e:
            print(f"  R={R:.3f}: FAILED - {e}")
            results.append({
                'R': R, 'E_total': np.nan, 'E_HF': np.nan,
                'time': 0, 'error': str(e),
            })

    R_arr = np.array([r['R'] for r in results])
    E_arr = np.array([r['E_total'] for r in results])
    E_HF_arr = np.array([r['E_HF'] for r in results])

    return {
        'R': R_arr,
        'E_total': E_arr,
        'E_HF': E_HF_arr,
        'results': results,
    }


# ============================================================
# Main
# ============================================================
def main():
    print("=" * 70)
    print("  H2 MINIMAL-BASIS CI ON PROLATE SPHEROIDAL LATTICE")
    print("  Date: 2026-03-13")
    print("=" * 70)

    # --- Reference values ---
    E_H = -0.5  # H atom energy
    E_H2_exact = -1.1745  # Exact H2 total energy at R_eq
    R_eq_exact = 1.401    # bohr
    D_e_exact = 0.1745    # Ha

    # --- Test 1: Single point at R=1.4 (near equilibrium) ---
    print("\n--- Test 1: Single point at R=1.4 bohr ---")
    res = h2_minimal_ci(R=1.4, N_xi_solve=5000, N_grid=50, verbose=True)

    E_err = abs(res['E_total'] - E_H2_exact) / abs(E_H2_exact) * 100
    D_e_ours = 2*E_H - res['E_total']
    print(f"\n  Results:")
    print(f"    E_total = {res['E_total']:.6f} Ha  (exact {E_H2_exact:.4f}, err {E_err:.1f}%)")
    print(f"    D_e = {D_e_ours:.6f} Ha  (exact {D_e_exact:.4f})")
    print(f"    Bound: {res['E_total'] < 2*E_H}")

    # --- Test 2: Grid convergence of V_ee at R=1.4 ---
    print("\n--- Test 2: V_ee grid convergence at R=1.4 ---")
    for N_grid in [30, 40, 50, 60]:
        res_conv = h2_minimal_ci(R=1.4, N_xi_solve=5000, N_grid=N_grid, verbose=False)
        print(f"  N_grid={N_grid:3d}: E_total={res_conv['E_total']:.6f}  "
              f"J_gg={res_conv['J_gg']:.6f}  V_gguu={res_conv['V_gguu']:.6f}  "
              f"[{res_conv['time']:.1f}s]")

    # --- Test 3: PES scan ---
    print("\n--- Test 3: H2 PES scan ---")
    R_vals = np.array([0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0])
    pes = scan_h2_pes(
        R_vals, N_xi_solve=5000, N_grid=50,
        xi_max_grid=15.0, verbose=False,
    )

    print(f"\n  {'R':>6s} {'E_CI':>12s} {'E_HF':>12s} {'D_e_CI':>10s} {'D_e_HF':>10s}")
    print("-" * 60)
    for i, R in enumerate(pes['R']):
        E_ci = pes['E_total'][i]
        E_hf = pes['E_HF'][i]
        De_ci = 2*E_H - E_ci
        De_hf = 2*E_H - E_hf
        print(f"  {R:6.3f} {E_ci:12.6f} {E_hf:12.6f} {De_ci:10.6f} {De_hf:10.6f}")

    # Fit spectroscopic constants
    valid = ~np.isnan(pes['E_total'])
    if np.sum(valid) >= 5:
        fit_ci = fit_spectroscopic_constants(pes['R'][valid], pes['E_total'][valid])
        fit_hf = fit_spectroscopic_constants(pes['R'][valid], pes['E_HF'][valid])

        print(f"\n  Spectroscopic constants:")
        print(f"  {'':20s} {'CI':>12s} {'HF':>12s} {'Exact':>12s}")
        print(f"  {'R_eq (bohr)':20s} {fit_ci['R_eq']:12.4f} {fit_hf['R_eq']:12.4f} {R_eq_exact:12.4f}")
        E_min_ci = fit_ci['E_min']
        E_min_hf = fit_hf['E_min']
        D_e_ci = 2*E_H - E_min_ci
        D_e_hf = 2*E_H - E_min_hf
        print(f"  {'E_min (Ha)':20s} {E_min_ci:12.6f} {E_min_hf:12.6f} {E_H2_exact:12.6f}")
        print(f"  {'D_e (Ha)':20s} {D_e_ci:12.6f} {D_e_hf:12.6f} {D_e_exact:12.6f}")

        R_err = abs(fit_ci['R_eq'] - R_eq_exact) / R_eq_exact * 100
        E_err = abs(E_min_ci - E_H2_exact) / abs(E_H2_exact) * 100
        D_err = abs(D_e_ci - D_e_exact) / D_e_exact * 100

        print(f"\n  CI errors:")
        print(f"    R_eq: {R_err:.1f}%")
        print(f"    E_min: {E_err:.1f}%")
        print(f"    D_e: {D_err:.1f}%")

        # Verdict
        print(f"\n  {'='*50}")
        if D_e_ci > 0:
            print(f"  H2 IS BOUND (D_e = {D_e_ci:.4f} Ha)")
            if R_err < 10:
                print(f"  R_eq within 10% -- SUCCESS")
            else:
                print(f"  R_eq off by {R_err:.0f}% -- NEEDS WORK")
        else:
            print(f"  H2 IS NOT BOUND -- FAIL")
        print(f"  {'='*50}")

    # Save results
    outpath = os.path.join(os.path.dirname(__file__), 'data', 'prolate_h2_ci.txt')
    with open(outpath, 'w') as f:
        f.write("# H2 Minimal-Basis CI on Prolate Spheroidal Lattice\n")
        f.write(f"# Date: 2026-03-13\n")
        f.write(f"# N_xi_solve=5000, N_grid=50\n")
        f.write(f"# E(H) = {E_H:.6f} Ha\n")
        f.write(f"#\n")
        f.write(f"# {'R':>8s} {'E_CI':>14s} {'E_HF':>14s} {'D_e_CI':>14s} {'D_e_HF':>14s}\n")
        for i, R in enumerate(pes['R']):
            E_ci = pes['E_total'][i]
            E_hf = pes['E_HF'][i]
            f.write(f"  {R:8.4f} {E_ci:14.8f} {E_hf:14.8f} "
                    f"{2*E_H-E_ci:14.8f} {2*E_H-E_hf:14.8f}\n")
    print(f"\n  Results saved to {outpath}")


if __name__ == '__main__':
    main()
