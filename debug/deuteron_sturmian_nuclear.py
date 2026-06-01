"""
Deuteron with Sturmian-like exponential basis for the relative coordinate
=========================================================================
The HO basis gives Gaussian tails that are too fat (Zemach +14%, Friar +51%).
Physical bound states have exponential tails exp(-gamma*r) where
gamma = sqrt(2*mu*B) ~ 0.23 fm^-1 for the deuteron.

This script replaces the HO relative-coordinate basis with exponential
(Laguerre-type) basis functions that have the correct asymptotic behavior,
while keeping the Minnesota NN potential.

Basis: phi_n(r; alpha) = N_n * r^l * exp(-alpha*r/2) * L_n^{2l+1}(alpha*r)
(Coulomb-Sturmian-like, with alpha as the scale parameter)

For the deuteron (l_rel = 0, S=1 channel only with Minnesota):
phi_n(r; alpha) = N_n * exp(-alpha*r/2) * L_n^1(alpha*r)

The Minnesota matrix elements are computed via numerical integration
on a radial grid (the Sturmian basis functions are known analytically).
"""
import sys, json
import numpy as np
from scipy.special import genlaguerre
from scipy.integrate import quad
from scipy.linalg import eigh, solve

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
sys.path.insert(0, '.')

from geovac.nuclear.minnesota import minnesota_potential, minnesota_params

HBAR_C = 197.3269804   # MeV*fm
M_NUCLEON = 938.918     # MeV/c^2 (average)
M_REDUCED = M_NUCLEON / 2  # reduced mass for equal-mass pn pair
E2_MEV_FM = 1.4399764  # e^2 in MeV*fm


def sturmian_radial(r, n, l, alpha):
    """Coulomb-Sturmian radial basis function.

    phi_n(r) = N_n * (alpha*r)^l * exp(-alpha*r/2) * L_n^{2l+1}(alpha*r)

    Normalized: integral |phi_n|^2 r^2 dr = 1
    """
    rho = alpha * np.asarray(r, dtype=float)
    L_poly = genlaguerre(n, 2*l + 1)(rho)
    psi = rho**l * np.exp(-rho / 2) * L_poly

    # Normalization via numerical integration
    def integrand(r_val):
        rho_val = alpha * r_val
        L_val = genlaguerre(n, 2*l + 1)(rho_val)
        f = rho_val**l * np.exp(-rho_val / 2) * L_val
        return f**2 * r_val**2

    norm_sq, _ = quad(integrand, 0, 50.0/alpha, limit=200)
    if norm_sq < 1e-30:
        return np.zeros_like(r, dtype=float)
    return psi / np.sqrt(norm_sq)


def build_sturmian_matrices(n_basis, l, alpha, S=1, n_grid=4000):
    """Build overlap, kinetic, and Minnesota potential matrices in Sturmian basis.

    H_{ij} = T_{ij} + V_{ij}
    S_{ij} = <phi_i | phi_j>
    T_{ij} = (hbar^2 / 2*mu) * <phi_i | (-d^2/dr^2 + l(l+1)/r^2 - 2/r * d/dr) | phi_j>

    In practice: T = integral phi_i(r) * [-hbar^2/(2*mu) * nabla^2] * phi_j(r) * r^2 dr
    Using integration by parts or the kinetic-energy matrix element formula.

    For simplicity, compute everything via numerical integration on a fine grid.
    """
    r_max = max(50.0 / alpha, 30.0)  # extend to capture the tail
    r_grid = np.linspace(1e-6, r_max, n_grid)
    dr = r_grid[1] - r_grid[0]

    # Pre-compute basis functions on grid
    phi = np.zeros((n_basis, n_grid))
    for n in range(n_basis):
        phi[n, :] = sturmian_radial(r_grid, n, l, alpha)

    # Overlap matrix
    S_mat = np.zeros((n_basis, n_basis))
    for i in range(n_basis):
        for j in range(i, n_basis):
            integrand = phi[i] * phi[j] * r_grid**2
            S_mat[i, j] = np.trapezoid(integrand, r_grid)
            S_mat[j, i] = S_mat[i, j]

    # Potential energy matrix (Minnesota in relative coordinate)
    V_nn = minnesota_potential(r_grid, S)
    V_mat = np.zeros((n_basis, n_basis))
    for i in range(n_basis):
        for j in range(i, n_basis):
            integrand = phi[i] * V_nn * phi[j] * r_grid**2
            V_mat[i, j] = np.trapezoid(integrand, r_grid)
            V_mat[j, i] = V_mat[i, j]

    # Kinetic energy: T = -hbar^2/(2*mu) * (d^2/dr^2 + 2/r * d/dr - l(l+1)/r^2)
    # Using the second-derivative via finite differences on the basis functions
    hbar2_over_2mu = HBAR_C**2 / (2 * M_REDUCED)

    T_mat = np.zeros((n_basis, n_basis))
    for j in range(n_basis):
        # Second derivative via central differences
        d2phi = np.zeros(n_grid)
        d2phi[1:-1] = (phi[j, 2:] - 2*phi[j, 1:-1] + phi[j, :-2]) / dr**2
        # First derivative
        dphi = np.zeros(n_grid)
        dphi[1:-1] = (phi[j, 2:] - phi[j, :-2]) / (2*dr)
        # Centrifugal
        cent = np.zeros(n_grid)
        mask = r_grid > 1e-10
        cent[mask] = l*(l+1) / r_grid[mask]**2 * phi[j, mask]
        # First derivative term: 2/r * dphi
        first_deriv = np.zeros(n_grid)
        first_deriv[mask] = 2.0 / r_grid[mask] * dphi[mask]

        # -nabla^2 = -(d^2/dr^2 + 2/r*d/dr) + l(l+1)/r^2
        neg_laplacian = -(d2phi + first_deriv) + cent

        for i in range(n_basis):
            integrand = phi[i] * hbar2_over_2mu * neg_laplacian * r_grid**2
            T_mat[i, j] = np.trapezoid(integrand, r_grid)

    # Symmetrize T (should be symmetric, small numerical asymmetry from FD)
    T_mat = (T_mat + T_mat.T) / 2

    H_mat = T_mat + V_mat

    return S_mat, T_mat, V_mat, H_mat, phi, r_grid


def solve_generalized_eigenvalue(H, S):
    """Solve the generalized eigenvalue problem H*c = E*S*c."""
    evals, evecs = eigh(H, S)
    return evals, evecs


def compute_observables(evals, evecs, phi, r_grid, S_mat, alpha, l=0):
    """Compute observables from the ground state."""
    gs_coeffs = evecs[:, 0]
    E_gs = evals[0]

    n_basis = len(gs_coeffs)
    n_grid = len(r_grid)

    # Ground state wavefunction on grid
    psi_gs = np.zeros(n_grid)
    for n in range(n_basis):
        psi_gs += gs_coeffs[n] * phi[n]

    # Normalize
    norm = np.trapezoid(psi_gs**2 * r_grid**2, r_grid)
    psi_gs /= np.sqrt(norm)

    # <r^2> for the relative coordinate
    r2_rel = np.trapezoid(psi_gs**2 * r_grid**4, r_grid)

    # <r> for Zemach
    r1_rel = np.trapezoid(psi_gs**2 * r_grid**3, r_grid)

    # <r^3> for Friar
    r3_rel = np.trapezoid(psi_gs**2 * r_grid**5, r_grid)

    # <r^4>
    r4_rel = np.trapezoid(psi_gs**2 * r_grid**6, r_grid)

    # Point-proton: r_p = r_rel/2
    r2_pp = r2_rel / 4
    r_pp = np.sqrt(r2_pp)

    # Charge radius
    r_p_intr = 0.8414  # fm
    r_n2_intr = -0.1155  # fm^2
    M_d = 938.272 + 939.565
    darwin_foldy = 3 * HBAR_C**2 / (4 * M_d**2)
    r2_d = r_p_intr**2 + r_n2_intr + r2_pp + darwin_foldy
    r_d = np.sqrt(max(r2_d, 0))

    # Zemach radius (point nucleon)
    # For equal-mass, point-nucleon: r_Z = <|r_rel|> / 2 * convolution factor
    # Actually, Zemach = integral integral rho_ch(r1) rho_mag(r2) |r1-r2| d3r1 d3r2
    # For point proton with density |psi(r_rel)|^2:
    # rho_p(r_p) = integral |psi(2*r_p)|^2 * 8 d(angles)  [change of variable]
    # This is complex. Simpler: for point nucleons and S-wave:
    # r_Z ≈ <r_rel>/2 * geometric_factor
    # The exact relation involves the Fourier transform.
    # For now, use the direct convolution via the charge form factor.

    # Charge form factor: F_ch(q) = <psi| exp(i*q*r_p) |psi> = <psi| exp(i*q*r_rel/2) |psi>
    # For S-wave: F_ch(q) = integral |psi(r)|^2 * sin(q*r/2)/(q*r/2) * r^2 dr
    q_grid = np.linspace(0.01, 20.0, 200)
    Fch = np.zeros(len(q_grid))
    for iq, q in enumerate(q_grid):
        # sin(q*r/2)/(q*r/2) with r = r_rel
        sinc_factor = np.sinc(q * r_grid / (2 * np.pi))  # np.sinc(x) = sin(pi*x)/(pi*x)
        # Actually np.sinc(x) = sin(pi*x)/(pi*x), so for sin(qr/2)/(qr/2):
        arg = q * r_grid / 2
        sinc_factor = np.ones_like(arg)
        mask = arg > 1e-10
        sinc_factor[mask] = np.sin(arg[mask]) / arg[mask]

        Fch[iq] = np.trapezoid(psi_gs**2 * sinc_factor * r_grid**2, r_grid)

    # Zemach radius from form factor: r_Z = -4 * dF_ch/d(q^2) |_{q^2=0} * ...
    # Actually: r_Z = -(48/pi) * integral_0^inf [F_ch(q)*F_mag(q) - 1] / q^4 dq
    # For point nucleons: F_mag = F_ch (same spatial distribution)
    # This integral is UV-divergent for point nucleons. Use the rms expansion instead.
    # r_Z = 2 * r_pp  (for Gaussian)
    # For exponential: r_Z ~ (4/3) * <r_p>

    # Just report the moments directly
    r1_pp = r1_rel / 2  # <r_p> = <r_rel>/2

    return {
        'E_gs': E_gs,
        'r2_rel': r2_rel,
        'r_rel': np.sqrt(r2_rel),
        'r1_rel': r1_rel,
        'r3_rel': r3_rel,
        'r2_pp': r2_pp,
        'r_pp': r_pp,
        'r1_pp': r1_pp,
        'r_d': r_d,
        'psi_gs': psi_gs,
        'Fch': Fch,
        'q_grid': q_grid,
    }


def main():
    print("=" * 72)
    print("DEUTERON WITH STURMIAN (EXPONENTIAL) BASIS")
    print("=" * 72)

    # Experimental values
    r_d_exp = 2.12799  # fm
    B_d_exp = 2.224575  # MeV (deuteron binding energy)
    gamma_exp = np.sqrt(2 * M_REDUCED * B_d_exp) / HBAR_C  # fm^-1
    print(f"\n  Experimental: B_d = {B_d_exp:.6f} MeV")
    print(f"  Asymptotic decay: gamma = sqrt(2*mu*B)/hbar_c = {gamma_exp:.4f} fm^-1")
    print(f"  Characteristic length: 1/gamma = {1/gamma_exp:.2f} fm")

    # Scan alpha (Sturmian scale parameter) and n_basis
    print(f"\n  Scanning alpha and n_basis for S=1, l=0 channel...")
    print(f"\n  {'alpha':>6} {'n':>3} {'E_gs':>8} {'B_d':>8} {'r_pp':>7} {'r_d':>7} "
          f"{'r_d%':>7} {'<r>_pp':>7} {'<r^3>':>8}")
    print(f"  " + "-" * 75)

    best_result = None
    best_err = float('inf')

    for n_basis in [3, 5, 8, 12]:
        for alpha in [0.2, 0.3, 0.4, 0.5, 0.7, 1.0]:
            try:
                S_mat, T_mat, V_mat, H_mat, phi, r_grid = build_sturmian_matrices(
                    n_basis, l=0, alpha=alpha, S=1, n_grid=6000
                )

                # Check overlap matrix conditioning
                cond = np.linalg.cond(S_mat)
                if cond > 1e12:
                    continue

                evals, evecs = solve_generalized_eigenvalue(H_mat, S_mat)
                E_gs = evals[0]

                # Binding energy (kinetic zero-point is included)
                # For the Sturmian basis, E_gs is the total energy of the relative motion
                # B_d = -E_gs (binding energy is positive for bound state)
                B_d = -E_gs

                obs = compute_observables(evals, evecs, phi, r_grid, S_mat, alpha)

                r_d_err = (obs['r_d'] - r_d_exp) / r_d_exp * 100
                B_err = (B_d - B_d_exp) / B_d_exp * 100

                print(f"  {alpha:>6.2f} {n_basis:>3} {E_gs:>8.3f} {B_d:>8.3f} "
                      f"{obs['r_pp']:>7.3f} {obs['r_d']:>7.4f} {r_d_err:>+6.1f}% "
                      f"{obs['r1_pp']:>7.3f} {obs['r3_rel']:>8.2f}")

                if abs(B_err) < abs(best_err) and B_d > 0:
                    best_err = B_err
                    best_result = {
                        'alpha': alpha, 'n_basis': n_basis,
                        'E_gs': E_gs, 'B_d': B_d,
                        **obs,
                    }
            except Exception as e:
                pass

    # Report best result
    if best_result:
        print(f"\n{'='*72}")
        print(f"BEST RESULT: alpha={best_result['alpha']:.2f}, n_basis={best_result['n_basis']}")
        print(f"{'='*72}")
        print(f"  B_d = {best_result['B_d']:.4f} MeV (exp: {B_d_exp:.4f}, "
              f"err: {(best_result['B_d']-B_d_exp)/B_d_exp*100:+.1f}%)")
        print(f"  r_pp = {best_result['r_pp']:.4f} fm")
        print(f"  r_d = {best_result['r_d']:.5f} fm (exp: {r_d_exp:.5f}, "
              f"err: {(best_result['r_d']-r_d_exp)/r_d_exp*100:+.1f}%)")
        print(f"  <r_rel> = {best_result['r1_rel']:.4f} fm")
        print(f"  <r^3_rel> = {best_result['r3_rel']:.2f} fm^3")

        # Compare with HO basis
        print(f"\n  Comparison with HO basis (hw=8):")
        print(f"  {'':>20} {'Sturmian':>10} {'HO(hw=8)':>10} {'Experiment':>10}")
        print(f"  {'r_d (charge)':>20} {best_result['r_d']:>10.4f} {'2.1186':>10} {r_d_exp:>10.5f}")
        print(f"  {'r_pp (point-p)':>20} {best_result['r_pp']:>10.4f} {'1.972':>10} {'~1.975':>10}")

    # Save
    if best_result:
        save_data = {k: float(v) if isinstance(v, (int, float, np.floating)) else None
                     for k, v in best_result.items()
                     if not isinstance(v, np.ndarray)}
        save_data['B_d_exp'] = B_d_exp
        save_data['r_d_exp'] = r_d_exp
        with open('debug/data/deuteron_sturmian_nuclear.json', 'w') as f:
            json.dump(save_data, f, indent=2, default=float)
        print(f"\n  Results saved to debug/data/deuteron_sturmian_nuclear.json")


if __name__ == '__main__':
    main()
