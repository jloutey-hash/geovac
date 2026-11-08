"""
Deuteron Polarizability via Coulomb-Sturmian (Laguerre) Basis
=============================================================
Solves the deuteron in RELATIVE coordinates as a one-body problem:
    [-hbar^2/(2*mu) * (d^2/dr^2 - l(l+1)/r^2) + V_Minnesota(r)] * u(r) = E * u(r)

Uses Laguerre basis functions with exponential tails:
    phi_{n,l}(r) = N_{n,l} * (2*alpha*r)^{l+1} * L_n^{2l+1}(2*alpha*r) * exp(-alpha*r)

where alpha is a scale parameter. When alpha ~ kappa (binding momentum),
the basis matches the deuteron's exponential tail exactly.

This should eliminate the HO hw-gap problem and converge the polarizability.
"""
import sys, io, json
import numpy as np
from scipy.special import gamma as gamma_fn, factorial
from scipy.integrate import quad
from scipy.linalg import eigh as scipy_eigh

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# Compatibility
if not hasattr(np, 'trapz'):
    np.trapz = np.trapezoid

# Constants
HBAR_C = 197.3269804  # MeV*fm
M_N = 938.918         # MeV (average nucleon)
MU_D = M_N / 2        # reduced mass n-p
E2_MEV_FM = 1.4399764 # e^2 in MeV*fm
KAPPA_D = np.sqrt(2 * MU_D * 2.2246) / HBAR_C  # 0.2316 fm^-1
ALPHA_E_EXP = 0.6328  # fm^3

# Minnesota triplet channel (S=1, deuteron)
V_R, K_R = 200.0, 1.487    # MeV, fm^-2
V_T, K_T = -91.4, 0.465    # MeV, fm^-2

def V_minnesota_triplet(r):
    """Minnesota potential in the triplet (S=1) channel."""
    return V_R * np.exp(-K_R * r**2) + V_T * np.exp(-K_T * r**2)


def laguerre_basis(r, n, l, alpha):
    """Normalized Laguerre basis function for radial Schrodinger.
    u_{n,l}(r) = r * R_{n,l}(r), where u satisfies the reduced eqn.

    phi_n(r; alpha, l) = N * (2*alpha*r)^(l+1) * L_n^(2l+1)(2*alpha*r) * exp(-alpha*r)

    Normalized: integral_0^inf |phi_n|^2 dr = 1
    """
    x = 2 * alpha * r
    # Normalization: int_0^inf x^{2l+2} [L_n^{2l+1}(x)]^2 e^{-x} dx/(2*alpha)
    # = Gamma(n + 2l + 2) / (n! * (2*alpha))
    norm_sq = gamma_fn(n + 2*l + 2) / (factorial(n, exact=True) * 2 * alpha)
    norm = 1.0 / np.sqrt(norm_sq)

    # Generalized Laguerre via recursion
    L = _laguerre(n, 2*l + 1, x)
    return norm * x**(l + 1) * L * np.exp(-x / 2)


def _laguerre(n, alpha_lag, x):
    """Generalized Laguerre polynomial L_n^alpha(x) via recursion."""
    if n == 0:
        return np.ones_like(x)
    elif n == 1:
        return 1 + alpha_lag - x
    else:
        L0 = np.ones_like(x)
        L1 = 1 + alpha_lag - x
        for k in range(2, n + 1):
            L2 = ((2*k - 1 + alpha_lag - x) * L1 - (k - 1 + alpha_lag) * L0) / k
            L0, L1 = L1, L2
        return L1


def build_matrices(n_max, l, alpha, n_quad=200):
    """Build overlap S, kinetic T, and potential V matrices in Laguerre basis.

    Uses Gauss-Laguerre-like quadrature on a transformed grid.
    """
    dim = n_max + 1

    # Quadrature grid (map to semi-infinite interval)
    # Use simple composite Simpson on [0, R_max] with R_max large enough
    R_max = 40.0 / alpha  # cover ~40 e-folding lengths
    r_grid = np.linspace(1e-10, R_max, n_quad)
    dr = r_grid[1] - r_grid[0]

    # Precompute basis functions on grid
    phi = np.zeros((dim, n_quad))
    dphi = np.zeros((dim, n_quad))
    for n in range(dim):
        phi[n, :] = laguerre_basis(r_grid, n, l, alpha)
        # Numerical derivative for kinetic energy
        h = 1e-5
        phi_p = laguerre_basis(r_grid + h, n, l, alpha)
        phi_m = laguerre_basis(r_grid - h, n, l, alpha)
        dphi[n, :] = (phi_p - phi_m) / (2 * h)

    # Overlap matrix
    S = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(i, dim):
            integrand = phi[i, :] * phi[j, :]
            S[i, j] = np.trapz(integrand, r_grid)
            S[j, i] = S[i, j]

    # Kinetic energy: T = -hbar^2/(2*mu) * [d^2/dr^2 - l(l+1)/r^2]
    # Using integration by parts: <i|T|j> = hbar^2/(2*mu) * [<i'|j'> + l(l+1) <i|1/r^2|j>]
    hbar2_over_2mu = HBAR_C**2 / (2 * MU_D)  # MeV * fm^2

    T = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(i, dim):
            # Derivative term
            integrand_d = dphi[i, :] * dphi[j, :]
            # Centrifugal term
            integrand_c = phi[i, :] * phi[j, :] * l * (l + 1) / (r_grid**2 + 1e-30)
            T[i, j] = hbar2_over_2mu * (np.trapz(integrand_d, r_grid) +
                                          np.trapz(integrand_c, r_grid))
            T[j, i] = T[i, j]

    # Potential matrix
    V_grid = V_minnesota_triplet(r_grid)
    V = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(i, dim):
            integrand = phi[i, :] * V_grid * phi[j, :]
            V[i, j] = np.trapz(integrand, r_grid)
            V[j, i] = V[i, j]

    return S, T, V, r_grid, phi


def build_dipole_matrix(n_max_s, n_max_p, alpha, n_quad=200):
    """Build dipole matrix elements <n',l=1|r|n,l=0> in Laguerre basis.

    The E1 operator r*Y_{10} connects l=0 to l=1.
    Matrix element: <phi_{n',l=1}| r | phi_{n,l=0}>
    (the angular part gives a Clebsch-Gordan factor of 1/sqrt(3) for m=0)
    """
    R_max = 40.0 / alpha
    r_grid = np.linspace(1e-10, R_max, n_quad)

    phi_s = np.zeros((n_max_s + 1, n_quad))
    phi_p = np.zeros((n_max_p + 1, n_quad))

    for n in range(n_max_s + 1):
        phi_s[n, :] = laguerre_basis(r_grid, n, 0, alpha)
    for n in range(n_max_p + 1):
        phi_p[n, :] = laguerre_basis(r_grid, n, 1, alpha)

    # Dipole matrix: D[n',n] = <n',l=1| r | n,l=0>
    D = np.zeros((n_max_p + 1, n_max_s + 1))
    for np_ in range(n_max_p + 1):
        for ns in range(n_max_s + 1):
            integrand = phi_p[np_, :] * r_grid * phi_s[ns, :]
            D[np_, ns] = np.trapz(integrand, r_grid)

    return D


def solve_deuteron_sturmian(alpha, n_max_s=10, n_max_p=10, n_quad=300):
    """Solve the deuteron in Sturmian basis and compute polarizability."""

    # S-wave (l=0) sector — contains the ground state
    S_s, T_s, V_s, r_grid, phi_s_grid = build_matrices(n_max_s, 0, alpha, n_quad)
    H_s = T_s + V_s

    # Generalized eigenvalue problem H c = E S c
    evals_s, evecs_s = scipy_eigh(H_s, S_s)
    E_gs = evals_s[0]
    c_gs = evecs_s[:, 0]

    # P-wave (l=1) sector — E1 intermediate states
    S_p, T_p, V_p, _, phi_p_grid = build_matrices(n_max_p, 1, alpha, n_quad)
    H_p = T_p + V_p

    evals_p, evecs_p = scipy_eigh(H_p, S_p)

    # Dipole matrix elements in the basis
    D_basis = build_dipole_matrix(n_max_s, n_max_p, alpha, n_quad)

    # Transform to eigenbasis
    # |gs> = sum_n c_gs[n] |phi_n^s>
    # |p_m> = sum_n c_p[n,m] |phi_n^p>
    # <p_m|r|gs> = sum_{n,n'} c_p[n',m] * D_basis[n',n] * c_gs[n]
    # = c_p[:,m]^T . D_basis . c_gs

    D_eigen = evecs_p.T @ D_basis @ c_gs  # array of <p_m|r|gs> for each p-wave eigenstate

    # Polarizability: alpha_E = 2 * e^2 * sum_m |<p_m|r|gs>|^2 / (E_m - E_gs)
    # Factor 2 accounts for the 3 Cartesian components (isotropic: 3 × 2/3 = 2)
    # Actually: for m=0 component of E1, the angular factor is 1/3.
    # alpha_E = (2*e^2/3) * sum_m |<m||r||0>|^2 / (E_m - E_0)
    # where the reduced matrix element |<m||r||0>|^2 = 3 * |<m,l=1|r|0,l=0>|^2
    # So alpha_E = 2 * e^2 * sum |<m|r_z|0>|^2 / (E_m - E_0) for z-component
    # with the 1/sqrt(3) angular factor already in D_basis? No — D_basis is
    # the RADIAL part only. The angular part gives C.G. coefficient.
    # For <l=1,m=0|cos(theta)|l=0,m=0> = 1/sqrt(3)
    # So full D_z matrix element = D_radial * (1/sqrt(3))
    # And alpha_E = 2*e^2 * sum |D_radial/sqrt(3)|^2 / (E_m - E_0)
    #            = (2*e^2/3) * sum |D_radial|^2 / (E_m - E_0)
    # But for isotropic system, summing over all m_l of the intermediate
    # state gives factor 3 (m_l = -1, 0, +1 all contribute equally)
    # So alpha_E = 2*e^2 * sum_m |D_radial|^2 / (3 * (E_m - E_0)) * 3
    #            = 2*e^2 * sum |D_radial|^2 / (E_m - E_0)
    # The factors cancel! alpha_E = 2*e^2 * sum |<p_m|r|gs>|^2 / (E_m - E_gs)

    # BUT: for the deuteron (n-p system with charge e only on proton),
    # the effective E1 operator in relative coordinates is:
    #   D_z = e * z_p - e*(Z/A)*z_cm = e * (1/2) * z_rel  (for Z=1, A=2)
    # So D_eff = (e/2) * r * cos(theta)
    # alpha_E = 2 * e^2 * (1/2)^2 * sum |<m|r|0>|^2 / (E_m - E_0)
    #         = (e^2/2) * sum |<m|r|0>|^2 / (E_m - E_0)

    alpha_E = 0.0
    n_contributing = 0
    for m in range(len(evals_p)):
        dE = evals_p[m] - E_gs
        if dE < 1e-6:
            continue
        alpha_E += 0.5 * E2_MEV_FM * D_eigen[m]**2 / dE
        if abs(D_eigen[m]) > 1e-10:
            n_contributing += 1

    return {
        'alpha': alpha,
        'E_gs': E_gs,
        'alpha_E': alpha_E,
        'n_max_s': n_max_s,
        'n_max_p': n_max_p,
        'n_contributing': n_contributing,
        'evals_s': evals_s[:5].tolist(),
        'evals_p': evals_p[:5].tolist(),
    }


# =========================================================================
print("=" * 72)
print("DEUTERON POLARIZABILITY: COULOMB-STURMIAN BASIS")
print("(exponential tails, relative-coordinate one-body problem)")
print("=" * 72)

print(f"\nDeuteron: kappa = {KAPPA_D:.4f} fm^-1, 1/kappa = {1/KAPPA_D:.2f} fm")
print(f"Minnesota triplet: V(r) = {V_R} exp(-{K_R}r^2) {V_T} exp(-{K_T}r^2) MeV")
print(f"Experimental alpha_E = {ALPHA_E_EXP:.4f} fm^3")

# Step 1: Scan over alpha (Sturmian scale parameter)
print(f"\n--- Step 1: alpha scan at n_max=10 ---")
alpha_values = np.arange(0.10, 1.01, 0.05)
scan_results = []

for alpha in alpha_values:
    res = solve_deuteron_sturmian(alpha, n_max_s=10, n_max_p=10, n_quad=400)
    scan_results.append(res)
    marker = " <-- kappa_d" if abs(alpha - KAPPA_D) < 0.03 else ""
    print(f"  alpha={alpha:.2f} fm^-1: E_gs={res['E_gs']:+.4f} MeV, "
          f"alpha_E={res['alpha_E']:.4f} fm^3 "
          f"({(res['alpha_E']-ALPHA_E_EXP)/ALPHA_E_EXP*100:+.1f}%){marker}")

# Step 2: Convergence in basis size at optimal alpha
print(f"\n--- Step 2: basis convergence at alpha near kappa ---")

# Find alpha that minimizes E_gs
best_alpha = min(scan_results, key=lambda x: x['E_gs'])['alpha']
print(f"  Variational-optimal alpha = {best_alpha:.2f} fm^-1")
print(f"  (deuteron kappa = {KAPPA_D:.2f} fm^-1)")

for n_max in [2, 4, 6, 8, 10, 15, 20]:
    res = solve_deuteron_sturmian(best_alpha, n_max_s=n_max, n_max_p=n_max, n_quad=500)
    print(f"  n_max={n_max:>2}: E_gs={res['E_gs']:+.4f} MeV, "
          f"alpha_E={res['alpha_E']:.4f} fm^3 ({(res['alpha_E']-ALPHA_E_EXP)/ALPHA_E_EXP*100:+.1f}%), "
          f"n_contrib={res['n_contributing']}")

# Step 3: Self-consistency check
print(f"\n--- Step 3: self-consistency (does variational alpha give right alpha_E?) ---")
res_opt = solve_deuteron_sturmian(best_alpha, n_max_s=15, n_max_p=15, n_quad=500)

# Find alpha that gives closest alpha_E to experiment
best_ae = min(scan_results, key=lambda x: abs(x['alpha_E'] - ALPHA_E_EXP))
alpha_ae = best_ae['alpha']

print(f"  Variational alpha = {best_alpha:.2f}: E_gs = {res_opt['E_gs']:+.4f} MeV, "
      f"alpha_E = {res_opt['alpha_E']:.4f} fm^3 ({(res_opt['alpha_E']-ALPHA_E_EXP)/ALPHA_E_EXP*100:+.1f}%)")
print(f"  Alpha-matched alpha = {alpha_ae:.2f}: E_gs = {best_ae['E_gs']:+.4f} MeV, "
      f"alpha_E = {best_ae['alpha_E']:.4f} fm^3 ({(best_ae['alpha_E']-ALPHA_E_EXP)/ALPHA_E_EXP*100:+.1f}%)")
print(f"  Gap: {best_alpha - alpha_ae:+.2f} fm^-1")

if abs(best_alpha - alpha_ae) < 0.1:
    print(f"  *** SELF-CONSISTENT: variational alpha gives correct alpha_E! ***")
else:
    print(f"  Gap exists but smaller than HO ({abs(best_alpha - alpha_ae):.2f} fm^-1 in alpha "
          f"vs 5 MeV in hw)")

# Step 4: Compare to HO
print(f"\n{'='*72}")
print(f"COMPARISON: STURMIAN vs HO")
print(f"{'='*72}")
print(f"\n  {'Method':<35} {'alpha_E':>10} {'vs exp':>10} {'Self-consistent?':>18}")
print(f"  {'-'*75}")
print(f"  {'HO N_shells=2, hw=8 (tuned)':35} {'0.6944':>10} {'+9.7%':>10} {'No (hw gap 5 MeV)':>18}")
print(f"  {'HO N_shells=2, hw=3 (variational)':35} {'4.9025':>10} {'+675%':>10} {'N/A':>18}")
res_sc = solve_deuteron_sturmian(best_alpha, n_max_s=15, n_max_p=15, n_quad=500)
print(f"  {'Sturmian n=15, alpha=opt':35} {res_sc['alpha_E']:>10.4f} "
      f"{(res_sc['alpha_E']-ALPHA_E_EXP)/ALPHA_E_EXP*100:>+9.1f}% "
      f"{'gap=' + f'{abs(best_alpha-alpha_ae):.2f}':>18}")
print(f"  {'Pionless EFT NNLO':35} {'0.63':>10} {'-0.4%':>10} {'(analytic)':>18}")
print(f"  {'SLEGS experiment 2026':35} {'0.6328':>10} {'---':>10} {'---':>18}")

# Save
output = {
    'scan': [{'alpha': r['alpha'], 'E_gs': r['E_gs'], 'alpha_E': r['alpha_E']}
             for r in scan_results],
    'best_alpha_var': best_alpha,
    'best_alpha_ae': alpha_ae,
    'kappa_d': KAPPA_D,
    'alpha_E_exp': ALPHA_E_EXP,
}
with open('debug/data/deuteron_sturmian.json', 'w') as f:
    json.dump(output, f, indent=2, default=float)
print(f"\nSaved to debug/data/deuteron_sturmian.json")
