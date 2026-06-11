"""Continuum limit of the Gaussian spectral action on S^3 Dirac.

In the continuum limit (n_max -> inf), the Dirac eigenvalues on unit S^3
are |lambda_n| = n + 3/2 with degeneracy g_n = 2(n+1)(n+2) (Camporesi-Higuchi).

Tr exp(-D^2/L^2) = sum_{n=0}^{inf} 2(n+1)(n+2) exp(-(n+3/2)^2/L^2)

Find L such that this equals K/pi = 43.619934066848226.
"""
import numpy as np
from scipy.optimize import brentq
from mpmath import mp, mpf, exp as mpexp, fsum, nstr

mp.dps = 50  # 50 decimal places

K_over_pi = mpf('43.619934066848226')
B = mpf('42')

def continuum_gaussian_action(L, n_terms=10000):
    """Compute sum_{n=0}^{n_terms} 2(n+1)(n+2) exp(-(n+3/2)^2/L^2)."""
    L2 = L * L
    terms = []
    for n in range(n_terms):
        g = 2 * (n+1) * (n+2)
        lam_sq = (n + mpf('1.5'))**2
        terms.append(g * mpexp(-lam_sq / L2))
    return fsum(terms)

def continuum_gaussian_action_np(L, n_terms=10000):
    """numpy version for root-finding."""
    ns = np.arange(n_terms)
    gs = 2.0 * (ns+1) * (ns+2)
    lam_sqs = (ns + 1.5)**2
    return float(np.sum(gs * np.exp(-lam_sqs / L**2)))

# Find the root
print("Continuum Gaussian spectral action: Tr exp(-D^2/L^2) on S^3 Dirac")
print("="*70)

# Check the function at various L
print("\nFunction evaluation:")
for L in [2.0, 3.0, 3.5, 3.7, 3.8, 4.0, 5.0, 10.0]:
    val = continuum_gaussian_action_np(L)
    print(f"  L = {L:.1f}: Tr exp(-D^2/L^2) = {val:.6f}")

# Root-find for K/pi
f_Kpi = lambda L: continuum_gaussian_action_np(L) - float(K_over_pi)
L_Kpi = brentq(f_Kpi, 1.0, 50.0, xtol=1e-15, rtol=1e-15)
val = continuum_gaussian_action_np(L_Kpi)
print(f"\nK/pi match (float):")
print(f"  Lambda = {L_Kpi:.15f}")
print(f"  value = {val:.15f}")
print(f"  target = {float(K_over_pi):.15f}")
print(f"  residual = {val - float(K_over_pi):.2e}")

# High-precision evaluation at the found Lambda
print(f"\nHigh-precision (50 digits):")
L_mp = mpf(L_Kpi)
val_mp = continuum_gaussian_action(L_mp)
print(f"  Lambda = {nstr(L_mp, 20)}")
print(f"  value = {nstr(val_mp, 30)}")
print(f"  target = {nstr(K_over_pi, 30)}")
print(f"  residual = {nstr(val_mp - K_over_pi, 10)}")

# Root-find for B=42
f_B = lambda L: continuum_gaussian_action_np(L) - float(B)
L_B = brentq(f_B, 1.0, 50.0, xtol=1e-15, rtol=1e-15)
val_B = continuum_gaussian_action_np(L_B)
print(f"\nB match (float):")
print(f"  Lambda = {L_B:.15f}")
print(f"  value = {val_B:.15f}")

# Structural analysis of Lambda
print(f"\n{'='*70}")
print(f"Structural analysis of Lambda_inf:")
print(f"  Lambda(K/pi) = {L_Kpi:.15f}")
print(f"  Lambda(B) = {L_B:.15f}")
print(f"  ratio = {L_Kpi/L_B:.15f}")
print(f"  K_pi/B = {float(K_over_pi/B):.15f}")

import math
pi = math.pi
sqrt_pi = math.sqrt(pi)

candidates = {
    'pi': pi,
    'sqrt(pi)': sqrt_pi,
    '2*sqrt(pi)': 2*sqrt_pi,
    'e': math.e,
    'sqrt(2*pi)': math.sqrt(2*pi),
    'sqrt(e*pi)': math.sqrt(math.e * pi),
    'pi/sqrt(e)': pi/math.sqrt(math.e),
    'sqrt(14)': math.sqrt(14),
    'sqrt(13)': math.sqrt(13),
    'sqrt(13.5)': math.sqrt(13.5),
    'sqrt(13.75)': math.sqrt(13.75),
    '2*sqrt(pi/e)': 2*math.sqrt(pi/math.e),
    '(3/2)*sqrt(pi)': 1.5*sqrt_pi,
    'sqrt(3)*sqrt(pi)': math.sqrt(3)*sqrt_pi,
    'cbrt(42)': 42**(1/3),
    'cbrt(K/pi)': float(K_over_pi)**(1/3),
    'pi**(2/3)': pi**(2/3),
}

print(f"\nCandidate identifications for Lambda(K/pi) = {L_Kpi:.10f}:")
sorted_candidates = sorted(candidates.items(), key=lambda x: abs(x[1] - L_Kpi))
for name, val in sorted_candidates[:10]:
    err = (val - L_Kpi) / L_Kpi
    print(f"  {name:25s} = {val:.10f}  rel_err = {err:+.6e}")

# Also check Lambda^2
L2 = L_Kpi**2
print(f"\nCandidate identifications for Lambda^2 = {L2:.10f}:")
candidates_sq = {
    'pi^2': pi**2,
    '2*pi': 2*pi,
    'e*pi': math.e*pi,
    '4*pi': 4*pi,
    'pi^2/e': pi**2/math.e,
    '14': 14.0,
    '13': 13.0,
    '13.5': 13.5,
    '42/pi': 42/pi,
    'K_pi/pi': float(K_over_pi)/pi,
    '3*pi': 3*pi,
    '4*e': 4*math.e,
    'e^2': math.e**2,
    'pi^2 - 3/2': pi**2 - 1.5,
    'pi + e': pi + math.e,
    'B/pi': float(B)/pi,
    '4*pi/3': 4*pi/3,
}
sorted_sq = sorted(candidates_sq.items(), key=lambda x: abs(x[1] - L2))
for name, val in sorted_sq[:10]:
    err = (val - L2) / L2
    print(f"  {name:25s} = {val:.10f}  rel_err = {err:+.6e}")

# Convergence from finite graphs to continuum
print(f"\n{'='*70}")
print("Convergence of finite-graph Lambda to continuum limit:")
finite_Ls = [
    (4, 4.812758589542093),
    (5, 4.023982804269919),
    (6, 3.813999560904021),
    (7, 3.743862595436180),
    (8, 3.720202039135583),
]
print(f"{'n_max':>5} {'Lambda_graph':>15} {'Lambda_cont':>15} {'diff':>12}")
for n, L_graph in finite_Ls:
    diff = L_graph - L_Kpi
    print(f"{n:>5} {L_graph:>15.10f} {L_Kpi:>15.10f} {diff:>+12.6e}")

# Check: does the graph-to-continuum convergence scale as 1/n_max?
print(f"\nGraph-continuum differences:")
diffs = [(n, L_graph - L_Kpi) for n, L_graph in finite_Ls]
for i in range(len(diffs)-1):
    n1, d1 = diffs[i]
    n2, d2 = diffs[i+1]
    if abs(d2) > 1e-15 and abs(d1) > 1e-15:
        ratio = d2/d1
        print(f"  n={n1}->{n2}: diff ratio = {ratio:.4f}, "
              f"n-ratio = {n1/n2:.4f}, n^2-ratio = {(n1/n2)**2:.4f}")
