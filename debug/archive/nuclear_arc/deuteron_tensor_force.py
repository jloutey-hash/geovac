"""
Deuteron with Minnesota + Gaussian tensor force (coupled ³S₁ + ³D₁)
====================================================================
Adds a Gaussian tensor term V_T(r) * S₁₂ to the Minnesota central potential.
The tensor operator S₁₂ = 3(sigma1·r_hat)(sigma2·r_hat) - sigma1·sigma2
couples l=0 (S-wave) to l=2 (D-wave), giving:
  - Nonzero quadrupole moment Q_d (experiment: 0.286 fm²)
  - Corrected magnetic moment mu_d (experiment: 0.857 n.m.)
  - D-state probability P_D ~ 4-7%

Uses the Sturmian (exponential) basis from the previous diagnostic,
extended to a two-channel coupled system.

In the coupled (L, S=1, J=1) basis:
  <³D₁ | S₁₂ | ³S₁> = sqrt(8)  (standard result)
  <³S₁ | S₁₂ | ³S₁> = -2
  <³D₁ | S₁₂ | ³D₁> = -2/5  (diagonal D-wave tensor)

The radial Schrodinger equation becomes a 2x2 coupled system:
  [T_S + V_c(S=1) + V_T*(-2)]   [V_T * sqrt(8)]      [u_S]       [u_S]
  [V_T * sqrt(8)]                [T_D + V_c(S=1) + V_T*(-2/5)]  [u_D]  = E [u_D]

where V_c is the central Minnesota, V_T = V_T0 * exp(-kappa_T * r²),
T_S has l=0 centrifugal, T_D has l=2 centrifugal.
"""
import sys, json, time
import numpy as np
from scipy.special import genlaguerre
from scipy.integrate import quad
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
sys.path.insert(0, '.')

from geovac.nuclear.minnesota import minnesota_potential

HBAR_C = 197.3269804
M_NUCLEON = 938.918
M_REDUCED = M_NUCLEON / 2
MU_P = 2.7928473
MU_N = -1.9130427

# Tensor S₁₂ matrix elements in the (L, S=1, J=1) basis
# Standard results: Bohr & Mottelson Vol. I; Ring & Schuck (1980)
S12_SS = 0.0         # <³S₁|S₁₂|³S₁> = 0 (angular average of rank-2 tensor in L=0)
S12_SD = np.sqrt(8)  # <³D₁|S₁₂|³S₁> = <³S₁|S₁₂|³D₁>
S12_DD = -2.0        # <³D₁|S₁₂|³D₁>


def sturmian_radial_on_grid(n, l, alpha, r_grid):
    """Normalized Sturmian basis function on a pre-computed grid."""
    rho = alpha * r_grid
    L_poly = genlaguerre(n, 2*l + 1)(rho)
    psi = rho**l * np.exp(-rho / 2) * L_poly
    norm_sq = np.trapezoid(psi**2 * r_grid**2, r_grid)
    if norm_sq < 1e-30:
        return np.zeros_like(r_grid)
    return psi / np.sqrt(norm_sq)


def build_coupled_channel_matrices(n_basis_S, n_basis_D, alpha_S, alpha_D,
                                     V_T0, kappa_T, n_grid=8000):
    """Build the coupled ³S₁ + ³D₁ Hamiltonian matrices.

    Returns (H, S_overlap) as block matrices:
      H = [[H_SS, H_SD],    S = [[S_SS, 0  ],
           [H_DS, H_DD]]         [0,   S_DD]]

    Blocks:
      H_SS = T(l=0) + V_central(S=1) + V_T0 * S12_SS * V_T_radial
      H_DD = T(l=2) + V_central(S=1) + V_T0 * S12_DD * V_T_radial
      H_SD = V_T0 * S12_SD * V_T_radial (off-diagonal tensor coupling)
      S_SS, S_DD = overlap matrices within each channel
    """
    r_max = max(50.0/alpha_S, 50.0/alpha_D, 30.0)
    r_grid = np.linspace(1e-6, r_max, n_grid)
    dr = r_grid[1] - r_grid[0]
    hbar2_2mu = HBAR_C**2 / (2 * M_REDUCED)

    # Build basis functions
    phi_S = np.zeros((n_basis_S, n_grid))
    for n in range(n_basis_S):
        phi_S[n] = sturmian_radial_on_grid(n, 0, alpha_S, r_grid)

    phi_D = np.zeros((n_basis_D, n_grid))
    for n in range(n_basis_D):
        phi_D[n] = sturmian_radial_on_grid(n, 2, alpha_D, r_grid)

    # Potentials on grid
    V_central = minnesota_potential(r_grid, S=1)
    V_tensor_radial = V_T0 * np.exp(-kappa_T * r_grid**2)

    dim = n_basis_S + n_basis_D

    def overlap_block(phi_a, phi_b):
        na, nb = phi_a.shape[0], phi_b.shape[0]
        S = np.zeros((na, nb))
        for i in range(na):
            for j in range(nb):
                S[i, j] = np.trapezoid(phi_a[i] * phi_b[j] * r_grid**2, r_grid)
        return S

    def potential_block(phi_a, phi_b, V):
        na, nb = phi_a.shape[0], phi_b.shape[0]
        M = np.zeros((na, nb))
        for i in range(na):
            for j in range(nb):
                M[i, j] = np.trapezoid(phi_a[i] * V * phi_b[j] * r_grid**2, r_grid)
        return M

    def kinetic_block(phi, l_val):
        n_b = phi.shape[0]
        T = np.zeros((n_b, n_b))
        for j in range(n_b):
            d2phi = np.zeros(n_grid)
            d2phi[1:-1] = (phi[j, 2:] - 2*phi[j, 1:-1] + phi[j, :-2]) / dr**2
            dphi = np.zeros(n_grid)
            dphi[1:-1] = (phi[j, 2:] - phi[j, :-2]) / (2*dr)
            cent = np.zeros(n_grid)
            mask = r_grid > 1e-10
            cent[mask] = l_val*(l_val+1) / r_grid[mask]**2 * phi[j, mask]
            first_deriv = np.zeros(n_grid)
            first_deriv[mask] = 2.0 / r_grid[mask] * dphi[mask]
            neg_lap = -(d2phi + first_deriv) + cent
            for i in range(n_b):
                T[i, j] = np.trapezoid(phi[i] * hbar2_2mu * neg_lap * r_grid**2, r_grid)
        return (T + T.T) / 2

    # Overlap
    S_SS = overlap_block(phi_S, phi_S)
    S_DD = overlap_block(phi_D, phi_D)
    S_full = np.block([
        [S_SS, np.zeros((n_basis_S, n_basis_D))],
        [np.zeros((n_basis_D, n_basis_S)), S_DD],
    ])

    # Kinetic
    T_SS = kinetic_block(phi_S, 0)
    T_DD = kinetic_block(phi_D, 2)

    # Central potential (same for S and D since both are S=1 triplet)
    Vc_SS = potential_block(phi_S, phi_S, V_central)
    Vc_DD = potential_block(phi_D, phi_D, V_central)

    # Tensor potential
    VT_SS = potential_block(phi_S, phi_S, V_tensor_radial) * S12_SS
    VT_DD = potential_block(phi_D, phi_D, V_tensor_radial) * S12_DD
    VT_SD = potential_block(phi_S, phi_D, V_tensor_radial) * S12_SD
    VT_DS = VT_SD.T

    # Assemble full H
    H_SS = T_SS + Vc_SS + VT_SS
    H_DD = T_DD + Vc_DD + VT_DD
    H_full = np.block([
        [H_SS, VT_SD],
        [VT_DS, H_DD],
    ])

    return H_full, S_full, phi_S, phi_D, r_grid


def solve_and_analyze(n_basis_S, n_basis_D, alpha_S, alpha_D, V_T0, kappa_T,
                       n_grid=8000, verbose=True):
    """Solve the coupled-channel problem and compute observables."""

    H, S_ov, phi_S, phi_D, r_grid = build_coupled_channel_matrices(
        n_basis_S, n_basis_D, alpha_S, alpha_D, V_T0, kappa_T, n_grid
    )

    cond = np.linalg.cond(S_ov)
    if cond > 1e14:
        return None

    evals, evecs = eigh(H, S_ov)
    E_gs = evals[0]
    B_d = -E_gs

    gs = evecs[:, 0]
    c_S = gs[:n_basis_S]
    c_D = gs[n_basis_S:]

    # Reconstruct wavefunctions on grid
    psi_S = np.zeros_like(r_grid)
    for n in range(n_basis_S):
        psi_S += c_S[n] * phi_S[n]

    psi_D = np.zeros_like(r_grid)
    for n in range(n_basis_D):
        psi_D += c_D[n] * phi_D[n]

    # Normalize total wavefunction
    norm_S = np.trapezoid(psi_S**2 * r_grid**2, r_grid)
    norm_D = np.trapezoid(psi_D**2 * r_grid**2, r_grid)
    norm_total = norm_S + norm_D
    psi_S /= np.sqrt(norm_total)
    psi_D /= np.sqrt(norm_total)
    norm_S /= norm_total
    norm_D /= norm_total

    P_D = norm_D  # D-state probability

    # <r²> for each channel
    r2_S = np.trapezoid(psi_S**2 * r_grid**4, r_grid)
    r2_D = np.trapezoid(psi_D**2 * r_grid**4, r_grid)
    r2_rel = r2_S + r2_D

    # Quadrupole moment
    # Q_d = (1/20) * <r²_D> * sqrt(8) / ...
    # More precisely, for the deuteron:
    # Q_d = (1/sqrt(50)) * integral u_S(r) * u_D(r) * r² dr
    #      + (1/20) * integral u_D(r)² * r² dr
    # where u = r * R (reduced radial wavefunction)
    # With our normalization (R functions, not u):
    # Q_d = (1/sqrt(50)) * integral psi_S * psi_D * r^4 dr
    #      + (1/20) * integral psi_D^2 * r^4 dr

    # Standard formula: Q_d = (1/sqrt(50)) * <S|r²|D> + (1/20) * <D|r²|D>
    # where |S> and |D> are the S and D components of the wavefunction
    SD_r2 = np.trapezoid(psi_S * psi_D * r_grid**4, r_grid)
    Q_d = (1.0/np.sqrt(50)) * SD_r2 + (1.0/20.0) * r2_D

    # Magnetic moment
    # mu_d = mu_S * (1 - P_D) * (3/2) + mu_D * P_D * (...)
    # Standard result for ³S₁-³D₁ mixture:
    # mu_d = (mu_p + mu_n) * (1 - 3/2 * P_D) + 3/4 * P_D
    # = (mu_p + mu_n) - (3/2 * (mu_p + mu_n) - 3/4) * P_D
    mu_S_wave = MU_P + MU_N
    mu_d = mu_S_wave * (1 - 1.5 * P_D) + 0.75 * P_D

    # Charge radius
    r2_pp = r2_rel / 4
    r_pp = np.sqrt(max(r2_pp, 0))
    r_p_intr = 0.8414
    r_n2_intr = -0.1155
    M_d = 938.272 + 939.565
    darwin_foldy = 3 * HBAR_C**2 / (4 * M_d**2)
    r2_d = r_p_intr**2 + r_n2_intr + r2_pp + darwin_foldy
    r_d = np.sqrt(max(r2_d, 0))

    result = {
        'E_gs': E_gs, 'B_d': B_d,
        'P_D': P_D, 'P_S': 1 - P_D,
        'Q_d': Q_d, 'mu_d': mu_d,
        'r_d': r_d, 'r_pp': r_pp,
        'r2_rel': r2_rel,
        'V_T0': V_T0, 'kappa_T': kappa_T,
        'alpha_S': alpha_S, 'alpha_D': alpha_D,
        'n_basis_S': n_basis_S, 'n_basis_D': n_basis_D,
        'cond_S': cond,
    }

    if verbose:
        print(f"  B_d = {B_d:.4f} MeV, P_D = {P_D*100:.2f}%")
        print(f"  Q_d = {Q_d:.4f} fm^2, mu_d = {mu_d:.4f} n.m.")
        print(f"  r_d = {r_d:.4f} fm, r_pp = {r_pp:.3f} fm")

    return result


def fit_tensor_to_Qd(n_basis_S=12, n_basis_D=8, alpha_S=1.1, alpha_D=0.8,
                      kappa_T=0.5, Q_target=0.2860):
    """Fit the tensor strength V_T0 to reproduce Q_d = 0.286 fm²."""

    print(f"Fitting V_T0 at kappa_T={kappa_T:.2f} fm^-2 to Q_d={Q_target:.4f} fm^2")
    print(f"  Basis: n_S={n_basis_S}, n_D={n_basis_D}, alpha_S={alpha_S}, alpha_D={alpha_D}")

    def objective(V_T0):
        r = solve_and_analyze(n_basis_S, n_basis_D, alpha_S, alpha_D,
                               V_T0, kappa_T, verbose=False)
        if r is None:
            return 1e6
        return (r['Q_d'] - Q_target)**2

    # Bracket search: tensor must be attractive (V_T0 < 0)
    print(f"\n  Scanning V_T0...")
    best_V = None
    best_cost = float('inf')
    for V_T0 in np.arange(-5, -80, -5):
        r = solve_and_analyze(n_basis_S, n_basis_D, alpha_S, alpha_D,
                               V_T0, kappa_T, verbose=False)
        if r is None:
            continue
        cost = (r['Q_d'] - Q_target)**2
        print(f"    V_T0={V_T0:>6.1f}: Q_d={r['Q_d']:.4f}, B_d={r['B_d']:.4f}, "
              f"P_D={r['P_D']*100:.2f}%")
        if cost < best_cost:
            best_cost = cost
            best_V = V_T0

    if best_V is None:
        print("  No valid solution found!")
        return None

    # Refine with golden section
    print(f"\n  Refining around V_T0 = {best_V:.1f} ...")
    result = minimize_scalar(objective, bounds=(best_V - 10, best_V + 10),
                              method='bounded')
    V_T0_opt = result.x

    print(f"\n  Optimal V_T0 = {V_T0_opt:.2f} MeV")
    final = solve_and_analyze(n_basis_S, n_basis_D, alpha_S, alpha_D,
                               V_T0_opt, kappa_T, verbose=True)
    return final


def main():
    print("=" * 72)
    print("DEUTERON WITH MINNESOTA + GAUSSIAN TENSOR FORCE")
    print("Coupled ³S₁ + ³D₁ channels")
    print("=" * 72)

    # Experimental values
    B_exp = 2.224575
    Q_exp = 0.2860
    mu_exp = 0.8574382
    r_d_exp = 2.12799

    print(f"\n  Experimental targets:")
    print(f"    B_d = {B_exp:.6f} MeV")
    print(f"    Q_d = {Q_exp:.4f} fm^2")
    print(f"    mu_d = {mu_exp:.7f} n.m.")
    print(f"    r_d = {r_d_exp:.5f} fm")

    # First: verify the pure central (V_T0=0) reproduces the Sturmian result
    print(f"\n{'='*60}")
    print("Step 1: Verify pure central (V_T0 = 0)")
    print(f"{'='*60}")

    r_central = solve_and_analyze(
        n_basis_S=14, n_basis_D=8, alpha_S=1.1, alpha_D=0.8,
        V_T0=0.0, kappa_T=0.5, verbose=True
    )
    if r_central:
        print(f"  P_D should be 0%: {r_central['P_D']*100:.4f}% (tensor=0)")

    # Step 2: Fit tensor at kappa_T = 0.5 fm^-2 (range ~ 1.4 fm, typical NN)
    print(f"\n{'='*60}")
    print("Step 2: Fit tensor strength at kappa_T = 0.5 fm^-2")
    print(f"{'='*60}")

    result_05 = fit_tensor_to_Qd(
        n_basis_S=14, n_basis_D=8, alpha_S=1.1, alpha_D=0.8,
        kappa_T=0.5, Q_target=Q_exp
    )

    # Step 3: Try different kappa_T values
    print(f"\n{'='*60}")
    print("Step 3: Scan kappa_T (tensor range)")
    print(f"{'='*60}")

    results_by_kappa = []
    for kappa_T in [0.3, 0.5, 0.7]:
        print(f"\n--- kappa_T = {kappa_T:.1f} fm^-2 ---")
        r = fit_tensor_to_Qd(
            n_basis_S=14, n_basis_D=8, alpha_S=1.1, alpha_D=0.8,
            kappa_T=kappa_T, Q_target=Q_exp
        )
        if r:
            results_by_kappa.append(r)

    # Summary
    print(f"\n{'='*72}")
    print("SUMMARY")
    print(f"{'='*72}")
    print(f"\n  {'kappa_T':>7} {'V_T0':>7} {'B_d':>7} {'B%':>7} {'Q_d':>7} {'mu_d':>7} "
          f"{'mu%':>7} {'P_D':>6} {'r_d':>7} {'r%':>7}")
    print(f"  " + "-" * 80)

    for r in results_by_kappa:
        B_err = (r['B_d'] - B_exp) / B_exp * 100
        mu_err = (r['mu_d'] - mu_exp) / mu_exp * 100
        r_err = (r['r_d'] - r_d_exp) / r_d_exp * 100
        print(f"  {r['kappa_T']:>7.2f} {r['V_T0']:>7.1f} {r['B_d']:>7.3f} {B_err:>+6.1f}% "
              f"{r['Q_d']:>7.4f} {r['mu_d']:>7.4f} {mu_err:>+6.1f}% "
              f"{r['P_D']*100:>5.1f}% {r['r_d']:>7.4f} {r_err:>+6.1f}%")

    print(f"\n  Experiment:")
    print(f"  {'':>7} {'':>7} {B_exp:>7.3f} {'':>7} {Q_exp:>7.4f} {mu_exp:>7.4f} "
          f"{'':>7} {'4-7%':>6} {r_d_exp:>7.4f}")

    # Save
    all_data = {
        'central_only': {k: float(v) if isinstance(v, (int, float, np.floating)) else v
                          for k, v in (r_central or {}).items()},
        'tensor_fits': [{k: float(v) if isinstance(v, (int, float, np.floating)) else v
                          for k, v in r.items()} for r in results_by_kappa],
        'experimental': {'B_d': B_exp, 'Q_d': Q_exp, 'mu_d': mu_exp, 'r_d': r_d_exp},
    }
    with open('debug/data/deuteron_tensor_force.json', 'w') as f:
        json.dump(all_data, f, indent=2, default=float)
    print(f"\n  Results saved to debug/data/deuteron_tensor_force.json")


if __name__ == '__main__':
    main()
