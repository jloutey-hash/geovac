"""Root-finding for exact Gaussian spectral action match to K/pi."""
import sys, numpy as np
from scipy.optimize import brentq
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent))
from geovac.spectral_triple import FockSpectralTriple
from sympy import Rational

K_over_pi = 43.619934066848226
B = 42.0

print("Gaussian spectral action root-finding: Tr exp(-D^2/L^2) = target")
print("="*80)

results = []

for n in range(4, 9):
    print(f"\nn_max={n}:", flush=True)
    st = FockSpectralTriple(n_max=n, kappa=Rational(-1,16),
                            j_type='kramers', adjacency_weights='cg')
    D = np.array(st.dirac_operator.tolist(), dtype=float)
    evals = np.linalg.eigvalsh(D)
    evals_sq = evals**2
    dim_H = len(evals)
    max_abs = max(abs(evals))

    def gauss_action(L):
        return float(np.sum(np.exp(-evals_sq / L**2)))

    # For K/pi target
    # The Gaussian action is monotonically increasing in L (from 0 to dim_H)
    # Find the root of gauss_action(L) - K_over_pi = 0
    f_Kpi = lambda L: gauss_action(L) - K_over_pi

    # Bracket: at L=0.01, action~0; at L=100, action~dim_H
    L_lo, L_hi = 0.1, 100.0
    if f_Kpi(L_lo) * f_Kpi(L_hi) > 0:
        print(f"  K/pi: cannot bracket (f(0.1)={f_Kpi(0.1):.4f}, f(100)={f_Kpi(100.0):.4f})")
        L_Kpi = None
    else:
        L_Kpi = brentq(f_Kpi, L_lo, L_hi, xtol=1e-15, rtol=1e-15)
        val_Kpi = gauss_action(L_Kpi)
        residual_Kpi = val_Kpi - K_over_pi
        rel_err_Kpi = abs(residual_Kpi) / K_over_pi
        print(f"  K/pi match: Lambda = {L_Kpi:.15f}")
        print(f"    value = {val_Kpi:.15f}")
        print(f"    target = {K_over_pi:.15f}")
        print(f"    residual = {residual_Kpi:.2e}")
        print(f"    rel_err = {rel_err_Kpi:.2e}")
        print(f"    Lambda / max|lambda_D| = {L_Kpi/max_abs:.6f}")
        print(f"    Lambda / (n_max + 0.5) = {L_Kpi/(n+0.5):.6f}")

    # For B=42 target
    f_B = lambda L: gauss_action(L) - B
    if f_B(L_lo) * f_B(L_hi) > 0:
        print(f"  B: cannot bracket")
        L_B = None
    else:
        L_B = brentq(f_B, L_lo, L_hi, xtol=1e-15, rtol=1e-15)
        val_B = gauss_action(L_B)
        residual_B = val_B - B
        rel_err_B = abs(residual_B) / B
        print(f"  B match: Lambda = {L_B:.15f}")
        print(f"    value = {val_B:.15f}")
        print(f"    residual = {residual_B:.2e}")
        print(f"    rel_err = {rel_err_B:.2e}")

    if L_Kpi and L_B:
        ratio = L_Kpi / L_B
        print(f"  Lambda(K/pi) / Lambda(B) = {ratio:.10f}")
        print(f"  (K/pi) / B = {K_over_pi/B:.10f}")
        print(f"  ratio of ratios = {ratio / (K_over_pi/B):.10f}")

    results.append({
        'n_max': n, 'dim_H': dim_H, 'max_abs': max_abs,
        'L_Kpi': L_Kpi, 'L_B': L_B,
    })

print("\n" + "="*80)
print("Lambda convergence:")
Ls = [r['L_Kpi'] for r in results if r['L_Kpi']]
ns = [r['n_max'] for r in results if r['L_Kpi']]
for i, r in enumerate(results):
    if r['L_Kpi']:
        print(f"  n_max={r['n_max']}: Lambda(K/pi) = {r['L_Kpi']:.10f}, "
              f"Lambda(B) = {r['L_B']:.10f}")

if len(Ls) >= 3:
    print(f"\nRichardson extrapolation (last 3 points):")
    L1, L2, L3 = Ls[-3], Ls[-2], Ls[-1]
    # Assume L(n) = L_inf + c/n^p
    # With 3 points, estimate L_inf
    # Simple linear extrapolation in 1/n:
    n1, n2, n3 = ns[-3], ns[-2], ns[-1]
    # Aitken delta-squared
    denom = (L3 - L2)**2 - (L3 - L2)*(L2 - L1)
    if abs(L3 - 2*L2 + L1) > 1e-15:
        L_aitken = L3 - (L3 - L2)**2 / (L3 - 2*L2 + L1)
        print(f"  Aitken delta^2: Lambda_inf ~ {L_aitken:.10f}")

    # Also try 1/n extrapolation
    # L = a + b/n
    # L2 = a + b/n2, L3 = a + b/n3
    b = (L3 - L2) / (1.0/n3 - 1.0/n2)
    a = L3 - b/n3
    print(f"  Linear 1/n extrapolation: Lambda_inf ~ {a:.10f}")

    # 1/n^2 extrapolation
    b2 = (L3 - L2) / (1.0/n3**2 - 1.0/n2**2)
    a2 = L3 - b2/n3**2
    print(f"  Linear 1/n^2 extrapolation: Lambda_inf ~ {a2:.10f}")

# Check Lambda against structural constants
print(f"\nStructural comparison:")
print(f"  pi = {np.pi:.10f}")
print(f"  sqrt(pi) = {np.sqrt(np.pi):.10f}")
print(f"  2*sqrt(pi) = {2*np.sqrt(np.pi):.10f}")
print(f"  e = {np.e:.10f}")
print(f"  sqrt(2*pi) = {np.sqrt(2*np.pi):.10f}")
for r in results:
    if r['L_Kpi']:
        L = r['L_Kpi']
        print(f"\n  n_max={r['n_max']}: Lambda = {L:.10f}")
        print(f"    Lambda^2 = {L**2:.10f}")
        print(f"    Lambda^2 / pi = {L**2/np.pi:.10f}")
        print(f"    Lambda^2 / (2*pi) = {L**2/(2*np.pi):.10f}")
        print(f"    Lambda / sqrt(pi) = {L/np.sqrt(np.pi):.10f}")
