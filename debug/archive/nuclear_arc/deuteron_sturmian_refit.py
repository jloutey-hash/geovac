"""
Deuteron Sturmian v3: Refit Minnesota V_T to give correct binding energy,
then compute the polarizability. Tests whether the Sturmian basis produces
the right alpha_E when the interaction is properly calibrated.
"""
import sys, io, json
import numpy as np
from scipy.special import gamma as gamma_fn
from scipy.linalg import eigh as scipy_eigh
from scipy.optimize import brentq

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

HBAR_C = 197.3269804
MU_D = 938.918 / 2
E2_MEV_FM = 1.4399764
KAPPA_D = np.sqrt(2 * MU_D * 2.2246) / HBAR_C
ALPHA_E_EXP = 0.6328
E_BIND_TARGET = -2.2246

V_R, K_R = 200.0, 1.487
K_T = 0.465

def laguerre_poly(n, a, x):
    if n == 0: return np.ones_like(x, dtype=float)
    if n == 1: return 1.0 + a - x
    L0, L1 = np.ones_like(x, dtype=float), 1.0 + a - x
    for k in range(2, n+1):
        L2 = ((2*k-1+a-x)*L1 - (k-1+a)*L0) / k
        L0, L1 = L1, L2
    return L1

def norm_factor(n, l):
    from math import factorial
    return np.sqrt(2.0 * factorial(n) / gamma_fn(n + 2*l + 2))

def gauss_laguerre(n_pts):
    from numpy.polynomial.laguerre import laggauss
    return laggauss(n_pts)

def build_H(n_max, l, alpha, V_T_val, n_gauss=100):
    """Build Hamiltonian (H = T + V) and overlap S."""
    dim = n_max + 1
    hbar2_2mu = HBAR_C**2 / (2 * MU_D)

    S = np.eye(dim)
    T = np.zeros((dim, dim))
    for n in range(dim):
        T[n,n] = hbar2_2mu * alpha**2 * (2*n + 2*l + 3) / 2
        if n+1 < dim:
            off = -hbar2_2mu * alpha**2 * np.sqrt((n+1)*(n+2*l+2)) / 2
            T[n,n+1] = off; T[n+1,n] = off

    x_gl, w_gl = gauss_laguerre(n_gauss)
    L_vals = np.zeros((dim, n_gauss))
    norms = np.array([norm_factor(n, l) for n in range(dim)])
    for n in range(dim):
        L_vals[n,:] = laguerre_poly(n, 2*l+1, x_gl)

    r_gl = x_gl / (2*alpha)
    V_gl = V_R * np.exp(-K_R * r_gl**2) + V_T_val * np.exp(-K_T * r_gl**2)

    V_mat = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(i, dim):
            integrand = norms[i]*norms[j] * x_gl**(2*l+2) * L_vals[i,:]*L_vals[j,:] * V_gl / (2*alpha)
            V_mat[i,j] = np.sum(w_gl * integrand)
            V_mat[j,i] = V_mat[i,j]

    return S, T + V_mat

def ground_state_energy(V_T_val, alpha, n_max, n_gauss=100):
    S, H = build_H(n_max, 0, alpha, V_T_val, n_gauss)
    evals, _ = scipy_eigh(H, S)
    return evals[0]

def solve_full(V_T_val, alpha, n_max_s=20, n_max_p=20, n_gauss=100):
    """Full solve: s-wave ground state + p-wave intermediates + polarizability."""
    S_s, H_s = build_H(n_max_s, 0, alpha, V_T_val, n_gauss)
    evals_s, evecs_s = scipy_eigh(H_s, S_s)
    E_gs = evals_s[0]
    c_gs = evecs_s[:,0]

    S_p, H_p = build_H(n_max_p, 1, alpha, V_T_val, n_gauss)
    evals_p, evecs_p = scipy_eigh(H_p, S_p)

    # Dipole matrix
    x_gl, w_gl = gauss_laguerre(n_gauss)
    norms_s = np.array([norm_factor(n, 0) for n in range(n_max_s+1)])
    norms_p = np.array([norm_factor(n, 1) for n in range(n_max_p+1)])
    L_s = np.zeros((n_max_s+1, n_gauss))
    L_p = np.zeros((n_max_p+1, n_gauss))
    for n in range(n_max_s+1): L_s[n,:] = laguerre_poly(n, 1, x_gl)
    for n in range(n_max_p+1): L_p[n,:] = laguerre_poly(n, 3, x_gl)

    D = np.zeros((n_max_p+1, n_max_s+1))
    for i in range(n_max_p+1):
        for j in range(n_max_s+1):
            D[i,j] = np.sum(w_gl * norms_p[i]*norms_s[j] * x_gl**4 * L_p[i,:]*L_s[j,:] / (4*alpha**2))

    D_eigen = evecs_p.T @ D @ c_gs
    alpha_E = sum(0.5*E2_MEV_FM*D_eigen[m]**2/(evals_p[m]-E_gs)
                  for m in range(len(evals_p)) if evals_p[m]-E_gs > 1e-6)

    return E_gs, alpha_E, evals_s[:6], evals_p[:6]


print("=" * 72)
print("DEUTERON STURMIAN: REFITTED MINNESOTA")
print("=" * 72)

# For each alpha, find V_T that gives E_gs = -2.2246 MeV
alpha_values = [0.20, 0.23, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.80, 1.00]
n_max = 20

print(f"\n--- Step 1: Refit V_T at each alpha to give E_gs = {E_BIND_TARGET} MeV ---")
print(f"   (Original V_T = -91.4 MeV)")

results = []
for alpha in alpha_values:
    # Find V_T such that E_gs(V_T) = E_BIND_TARGET
    def objective(V_T_val):
        return ground_state_energy(V_T_val, alpha, n_max) - E_BIND_TARGET

    # Search: V_T near 0 gives unbound; V_T very negative gives deeply bound
    try:
        # Check brackets
        e_weak = ground_state_energy(-10.0, alpha, n_max)
        e_strong = ground_state_energy(-200.0, alpha, n_max)
        if (e_weak - E_BIND_TARGET) * (e_strong - E_BIND_TARGET) > 0:
            # Try wider range
            e_strong = ground_state_energy(-500.0, alpha, n_max)

        V_T_fit = brentq(objective, -500.0, -1.0, xtol=1e-6)
    except Exception as e:
        print(f"  alpha={alpha:.2f}: refit FAILED ({e})")
        continue

    E_gs_check, alpha_E, evals_s, evals_p = solve_full(V_T_fit, alpha, n_max, n_max)
    residual = (alpha_E - ALPHA_E_EXP) / ALPHA_E_EXP * 100

    print(f"  alpha={alpha:.2f}: V_T={V_T_fit:+8.3f} MeV, E_gs={E_gs_check:+.4f} MeV, "
          f"alpha_E={alpha_E:.4f} fm^3 ({residual:+.1f}%)")

    results.append({
        'alpha': alpha, 'V_T_fit': V_T_fit,
        'E_gs': E_gs_check, 'alpha_E': alpha_E,
        'residual_pct': residual,
    })

# Step 2: alpha stability of the refitted polarizability
print(f"\n--- Step 2: Alpha-stability of refitted alpha_E ---")
aE_values = [r['alpha_E'] for r in results]
aE_mean = np.mean(aE_values)
aE_std = np.std(aE_values)
print(f"  alpha_E across alpha values: {aE_mean:.4f} +/- {aE_std:.4f} fm^3")
print(f"  Spread: {aE_std/aE_mean*100:.1f}%")
print(f"  vs experiment: {(aE_mean-ALPHA_E_EXP)/ALPHA_E_EXP*100:+.1f}%")

if aE_std / aE_mean < 0.1:
    print(f"  *** ALPHA-STABLE: polarizability converged across alpha values! ***")
elif aE_std / aE_mean < 0.3:
    print(f"  *** PARTIALLY STABLE: moderate alpha dependence ***")
else:
    print(f"  *** NOT STABLE: strong alpha dependence ***")

# Step 3: Basis convergence at best alpha
best = min(results, key=lambda r: abs(r['alpha_E'] - ALPHA_E_EXP))
print(f"\n--- Step 3: Basis convergence at alpha={best['alpha']:.2f}, V_T={best['V_T_fit']:.3f} ---")
for nm in [4, 6, 8, 10, 15, 20, 25, 30]:
    E, aE, _, _ = solve_full(best['V_T_fit'], best['alpha'], nm, nm)
    print(f"  n_max={nm:>2}: E_gs={E:+.6f} MeV, alpha_E={aE:.4f} fm^3 "
          f"({(aE-ALPHA_E_EXP)/ALPHA_E_EXP*100:+.1f}%)")

# Summary
print(f"\n{'='*72}")
print(f"FINAL COMPARISON")
print(f"{'='*72}")
print(f"\n  {'Method':<45} {'alpha_E':>10} {'vs exp':>8}")
print(f"  {'-'*65}")
print(f"  {'HO N=2 hw=8 (tuned hw, orig Minnesota)':45} {'0.6944':>10} {'+9.7%':>8}")
print(f"  {'Sturmian n=20 (refit V_T, alpha-averaged)':45} {aE_mean:>10.4f} {(aE_mean-ALPHA_E_EXP)/ALPHA_E_EXP*100:>+7.1f}%")
print(f"  {'Sturmian n=20 (refit V_T, best alpha)':45} {best['alpha_E']:>10.4f} {best['residual_pct']:>+7.1f}%")
print(f"  {'Pionless EFT NNLO':45} {'0.63':>10} {'-0.4%':>8}")
print(f"  {'SLEGS experiment 2026':45} {'0.6328':>10} {'---':>8}")

with open('debug/data/deuteron_sturmian_refit.json', 'w') as f:
    json.dump(results, f, indent=2, default=float)
print(f"\nSaved to debug/data/deuteron_sturmian_refit.json")
