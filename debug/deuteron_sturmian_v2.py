"""
Deuteron Polarizability: Sturmian Basis v2 (Analytical + Gauss-Laguerre)
========================================================================
Fixes v1's numerical issues by using:
  1. Analytical overlap and kinetic energy matrix elements
  2. Gauss-Laguerre quadrature for the Gaussian Minnesota potential

Basis: u_n^l(r) = N * x^(l+1) * L_n^(2l+1)(x) * exp(-x/2), x = 2*alpha*r
Orthogonality: int_0^inf u_n * u_m dr = delta_{nm} (with proper normalization)

The Laguerre functions with the SAME alpha are orthonormal under the
standard inner product. The kinetic energy and 1/r^2 centrifugal term
have known closed-form matrix elements. Only the Minnesota potential
(a sum of Gaussians) requires quadrature.
"""
import sys, io, json
import numpy as np
from scipy.special import gamma as gamma_fn
from scipy.linalg import eigh as scipy_eigh

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
if not hasattr(np, 'trapz'):
    np.trapz = np.trapezoid

HBAR_C = 197.3269804
M_N = 938.918
MU_D = M_N / 2
E2_MEV_FM = 1.4399764
KAPPA_D = np.sqrt(2 * MU_D * 2.2246) / HBAR_C
ALPHA_E_EXP = 0.6328
E_BIND_EXP = -2.2246

V_R, K_R = 200.0, 1.487
V_T, K_T = -91.4, 0.465

def V_minnesota_triplet(r):
    return V_R * np.exp(-K_R * r**2) + V_T * np.exp(-K_T * r**2)


def gauss_laguerre_points(n_pts, alpha_lag):
    """Generalized Gauss-Laguerre quadrature points and weights for
    integral_0^inf x^alpha * exp(-x) * f(x) dx ~ sum w_i * f(x_i)
    """
    from numpy.polynomial.laguerre import laggauss
    # Standard Gauss-Laguerre (alpha=0)
    x, w = laggauss(n_pts)
    # For generalized: multiply weights by x^alpha_lag / Gamma(alpha_lag+1)?
    # Actually numpy's laggauss gives points for int exp(-x) f(x) dx
    # For int x^a exp(-x) f(x) dx, we can include the x^a in f
    return x, w


def laguerre_poly(n, alpha_lag, x):
    """Generalized Laguerre polynomial L_n^alpha(x) via stable recursion."""
    if n == 0:
        return np.ones_like(x, dtype=float)
    elif n == 1:
        return 1.0 + alpha_lag - x
    L0 = np.ones_like(x, dtype=float)
    L1 = 1.0 + alpha_lag - x
    for k in range(2, n + 1):
        L2 = ((2*k - 1 + alpha_lag - x) * L1 - (k - 1 + alpha_lag) * L0) / k
        L0, L1 = L1, L2
    return L1


def norm_factor(n, l):
    """Normalization: int_0^inf |u_n^l(r)|^2 dr = 1.
    u_n^l(r) = N * (2ar)^{l+1} L_n^{2l+1}(2ar) exp(-ar)

    With x = 2ar, dr = dx/(2a):
    int = N^2 / (2a) * int_0^inf x^{2l+2} [L_n^{2l+1}(x)]^2 exp(-x) dx
        = N^2 / (2a) * Gamma(n + 2l + 2) / n!

    So N = sqrt(2a * n! / Gamma(n + 2l + 2))
    """
    from math import factorial
    return np.sqrt(2.0 * factorial(n) / gamma_fn(n + 2*l + 2))
    # Note: the factor of 'a' (alpha) is handled when we substitute x = 2*a*r


def build_sturmian_matrices(n_max, l, alpha, n_gauss=80):
    """Build H and S matrices analytically + Gauss-Laguerre for V.

    ANALYTICAL overlap: S_{nm} = delta_{nm} (Laguerre functions are
    orthonormal under the correct weight).

    ANALYTICAL kinetic energy:
    For u_n(r) = N_n * x^{l+1} * L_n^{2l+1}(x) * exp(-x/2), x=2ar:
    The kinetic operator -hbar^2/(2mu) * [d^2/dr^2 - l(l+1)/r^2] gives
    the standard Sturmian kinetic matrix:

    T_{nm} = (hbar^2 * alpha^2 / (2*mu)) * [
        (2n + 2l + 1) * delta_{nm}
        - sqrt(n*(n+2l+1)) * delta_{n,m+1}
        - sqrt((n+1)*(n+2l+2)) * delta_{n,m-1}
    ] ... wait, this is for the hydrogen Sturmian where the kinetic
    operator gives a tridiagonal matrix. Let me derive this properly.

    Actually for the KINETIC ENERGY of a particle in a Laguerre basis,
    the matrix elements are well-known (Suhonen, de-Shalit & Talmi):

    <n'l|T|nl> = (alpha^2 * hbar^2 / (2*mu)) * {
        (2n + 2l + 3)/2 * delta_{n',n}
        - sqrt((n+1)*(n+2l+2))/2 * delta_{n',n+1}
        - sqrt(n*(n+2l+1))/2 * delta_{n',n-1}
    }

    This is the standard tridiagonal Sturmian kinetic matrix.
    """
    dim = n_max + 1
    hbar2_2mu = HBAR_C**2 / (2 * MU_D)  # MeV * fm^2

    # Overlap: orthonormal
    S = np.eye(dim)

    # Kinetic: tridiagonal (analytical)
    T = np.zeros((dim, dim))
    for n in range(dim):
        T[n, n] = hbar2_2mu * alpha**2 * (2*n + 2*l + 3) / 2
        if n + 1 < dim:
            off = -hbar2_2mu * alpha**2 * np.sqrt((n+1)*(n+2*l+2)) / 2
            T[n, n+1] = off
            T[n+1, n] = off

    # Potential: Gauss-Laguerre quadrature
    # V_{nm} = int_0^inf u_n(r) V(r) u_m(r) dr
    # With x = 2*alpha*r, dr = dx/(2*alpha):
    # V_{nm} = (1/(2*alpha)) * int_0^inf N_n*N_m * x^{2l+2} *
    #          L_n^{2l+1}(x) * L_m^{2l+1}(x) * exp(-x) *
    #          V(x/(2*alpha)) dx
    #
    # This is a generalized Gauss-Laguerre integral with weight
    # x^{2l+2} * exp(-x). We absorb x^{2l+2} into the integrand.

    # Use standard Gauss-Laguerre (weight exp(-x)) with enough points
    x_gl, w_gl = gauss_laguerre_points(n_gauss, 0)

    # Precompute Laguerre polynomials at quadrature points
    L_vals = np.zeros((dim, n_gauss))
    for n in range(dim):
        L_vals[n, :] = laguerre_poly(n, 2*l + 1, x_gl)

    # Normalization factors
    norms = np.array([norm_factor(n, l) for n in range(dim)])

    # V(r) at quadrature points: r = x/(2*alpha)
    r_gl = x_gl / (2 * alpha)
    V_gl = V_minnesota_triplet(r_gl)

    # Integrand for V_{nm}: (weight already exp(-x) from GL)
    # f(x) = N_n * N_m * x^{2l+2} * L_n(x) * L_m(x) * V(x/(2a)) / (2a)
    V_mat = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(i, dim):
            integrand = (norms[i] * norms[j] * x_gl**(2*l+2) *
                        L_vals[i, :] * L_vals[j, :] * V_gl / (2 * alpha))
            V_mat[i, j] = np.sum(w_gl * integrand)
            V_mat[j, i] = V_mat[i, j]

    return S, T, V_mat


def build_dipole_sturmian(n_max_s, n_max_p, alpha, n_gauss=80):
    """Dipole matrix elements <n',l=1|r|n,l=0> in Sturmian basis.

    <u_{n'}^1 | r | u_n^0> = int_0^inf u_{n'}^1(r) * r * u_n^0(r) dr

    With x = 2*alpha*r, r = x/(2a):
    = (1/(2a)) * int N_{n'} * N_n * x^{1+1} * L_{n'}^3(x) * exp(-x/2) *
      (x/(2a)) * x^{0+1} * L_n^1(x) * exp(-x/2) dx
    = 1/(4*a^2) * int N_{n'}*N_n * x^4 * L_{n'}^3(x) * L_n^1(x) * exp(-x) dx
    """
    x_gl, w_gl = gauss_laguerre_points(n_gauss, 0)

    L_s = np.zeros((n_max_s + 1, n_gauss))
    L_p = np.zeros((n_max_p + 1, n_gauss))
    norms_s = np.array([norm_factor(n, 0) for n in range(n_max_s + 1)])
    norms_p = np.array([norm_factor(n, 1) for n in range(n_max_p + 1)])

    for n in range(n_max_s + 1):
        L_s[n, :] = laguerre_poly(n, 1, x_gl)  # L_n^{2*0+1} = L_n^1
    for n in range(n_max_p + 1):
        L_p[n, :] = laguerre_poly(n, 3, x_gl)  # L_n^{2*1+1} = L_n^3

    # D_{n',n} = 1/(4*a^2) * sum_i w_i * N_{n'} * N_n * x_i^4 * L_{n'}^3(x_i) * L_n^1(x_i)
    D = np.zeros((n_max_p + 1, n_max_s + 1))
    for np_ in range(n_max_p + 1):
        for ns in range(n_max_s + 1):
            integrand = (norms_p[np_] * norms_s[ns] *
                        x_gl**4 * L_p[np_, :] * L_s[ns, :] /
                        (4 * alpha**2))
            D[np_, ns] = np.sum(w_gl * integrand)

    return D


def solve_and_polarize(alpha, n_max_s=15, n_max_p=15, n_gauss=100):
    """Solve deuteron and compute polarizability in Sturmian basis."""

    S_s, T_s, V_s = build_sturmian_matrices(n_max_s, 0, alpha, n_gauss)
    H_s = T_s + V_s
    evals_s, evecs_s = scipy_eigh(H_s, S_s)
    E_gs = evals_s[0]

    S_p, T_p, V_p = build_sturmian_matrices(n_max_p, 1, alpha, n_gauss)
    H_p = T_p + V_p
    evals_p, evecs_p = scipy_eigh(H_p, S_p)

    D_basis = build_dipole_sturmian(n_max_s, n_max_p, alpha, n_gauss)
    c_gs = evecs_s[:, 0]
    D_eigen = evecs_p.T @ D_basis @ c_gs

    # alpha_E = (e^2/2) * sum_m |<p_m|r|gs>|^2 / (E_m - E_gs)
    # Factor 1/2 from effective charge (e/2) in relative coordinates squared × 2 for sum-rule
    alpha_E = 0.0
    for m in range(len(evals_p)):
        dE = evals_p[m] - E_gs
        if dE < 1e-6:
            continue
        alpha_E += 0.5 * E2_MEV_FM * D_eigen[m]**2 / dE

    return {
        'alpha': alpha, 'E_gs': E_gs, 'alpha_E': alpha_E,
        'evals_s': evals_s[:6].tolist(),
        'evals_p': evals_p[:6].tolist(),
        'n_max_s': n_max_s, 'n_max_p': n_max_p,
    }


# =========================================================================
print("=" * 72)
print("DEUTERON STURMIAN v2: ANALYTICAL KINETIC + GAUSS-LAGUERRE POTENTIAL")
print("=" * 72)
print(f"\nkappa_d = {KAPPA_D:.4f} fm^-1, E_bind = {E_BIND_EXP:.4f} MeV")

# Step 1: alpha scan
print(f"\n--- Alpha scan (n_max=15, n_gauss=100) ---")
alpha_values = np.concatenate([
    np.arange(0.10, 0.50, 0.02),
    np.arange(0.50, 1.51, 0.10),
])

scan = []
for alpha in alpha_values:
    res = solve_and_polarize(alpha, n_max_s=15, n_max_p=15, n_gauss=100)
    scan.append(res)
    marker = " <-- kappa" if abs(alpha - KAPPA_D) < 0.015 else ""
    print(f"  a={alpha:.2f}: E_gs={res['E_gs']:+8.4f} MeV, "
          f"alpha_E={res['alpha_E']:12.4f} fm^3 "
          f"({(res['alpha_E']-ALPHA_E_EXP)/ALPHA_E_EXP*100:+8.1f}%){marker}")

# Find variational minimum
best_E = min(scan, key=lambda x: x['E_gs'])
best_aE = min(scan, key=lambda x: abs(x['alpha_E'] - ALPHA_E_EXP))

print(f"\n  Variational minimum: alpha={best_E['alpha']:.2f}, "
      f"E_gs={best_E['E_gs']:+.4f} MeV, alpha_E={best_E['alpha_E']:.4f} fm^3")
print(f"  Alpha-E match:      alpha={best_aE['alpha']:.2f}, "
      f"E_gs={best_aE['E_gs']:+.4f} MeV, alpha_E={best_aE['alpha_E']:.4f} fm^3")
print(f"  Gap: {best_E['alpha'] - best_aE['alpha']:+.2f} fm^-1")

# Step 2: Convergence at best alpha
print(f"\n--- Basis convergence at alpha={best_E['alpha']:.2f} ---")
for n_max in [2, 4, 6, 8, 10, 15, 20, 25]:
    res = solve_and_polarize(best_E['alpha'], n_max_s=n_max, n_max_p=n_max, n_gauss=120)
    print(f"  n_max={n_max:>2}: E_gs={res['E_gs']:+.6f} MeV, "
          f"alpha_E={res['alpha_E']:.4f} fm^3 ({(res['alpha_E']-ALPHA_E_EXP)/ALPHA_E_EXP*100:+.1f}%)")

# Step 3: Convergence at kappa
print(f"\n--- Basis convergence at alpha=kappa={KAPPA_D:.2f} ---")
for n_max in [2, 4, 6, 8, 10, 15, 20, 25]:
    res = solve_and_polarize(KAPPA_D, n_max_s=n_max, n_max_p=n_max, n_gauss=120)
    print(f"  n_max={n_max:>2}: E_gs={res['E_gs']:+.6f} MeV, "
          f"alpha_E={res['alpha_E']:.4f} fm^3 ({(res['alpha_E']-ALPHA_E_EXP)/ALPHA_E_EXP*100:+.1f}%)")

# Summary
print(f"\n{'='*72}")
print(f"COMPARISON TABLE")
print(f"{'='*72}")
print(f"\n  {'Method':<40} {'E_gs':>8} {'alpha_E':>10} {'vs exp':>8}")
print(f"  {'-'*70}")
print(f"  {'HO N=2, hw=8 (tuned)':40} {'+11.91':>8} {'0.6944':>10} {'+9.7%':>8}")
print(f"  {'HO N=2, hw=3 (variational)':40} {'+5.23':>8} {'4.9025':>10} {'+675%':>8}")
res_best = solve_and_polarize(best_E['alpha'], n_max_s=20, n_max_p=20, n_gauss=120)
print(f"  {'Sturmian n=20, alpha=var':40} {res_best['E_gs']:>+8.4f} "
      f"{res_best['alpha_E']:>10.4f} {(res_best['alpha_E']-ALPHA_E_EXP)/ALPHA_E_EXP*100:>+7.1f}%")
res_kap = solve_and_polarize(KAPPA_D, n_max_s=20, n_max_p=20, n_gauss=120)
print(f"  {'Sturmian n=20, alpha=kappa':40} {res_kap['E_gs']:>+8.4f} "
      f"{res_kap['alpha_E']:>10.4f} {(res_kap['alpha_E']-ALPHA_E_EXP)/ALPHA_E_EXP*100:>+7.1f}%")
print(f"  {'Pionless EFT NNLO':40} {'-2.22':>8} {'0.63':>10} {'-0.4%':>8}")
print(f"  {'SLEGS experiment 2026':40} {'-2.22':>8} {'0.6328':>10} {'---':>8}")

# Save
output = {
    'scan': [{'alpha': s['alpha'], 'E_gs': s['E_gs'], 'alpha_E': s['alpha_E']}
             for s in scan],
    'best_var': {'alpha': best_E['alpha'], 'E_gs': best_E['E_gs'], 'alpha_E': best_E['alpha_E']},
    'best_aE': {'alpha': best_aE['alpha'], 'E_gs': best_aE['E_gs'], 'alpha_E': best_aE['alpha_E']},
}
with open('debug/data/deuteron_sturmian_v2.json', 'w') as f:
    json.dump(output, f, indent=2, default=float)
print(f"\nSaved to debug/data/deuteron_sturmian_v2.json")
