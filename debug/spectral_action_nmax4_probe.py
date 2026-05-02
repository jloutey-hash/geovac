"""Quick n_max=4 probe: does dim_H=80 allow K/pi match?"""
import sys, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent))
from geovac.spectral_triple import FockSpectralTriple
from sympy import Rational

st = FockSpectralTriple(n_max=4, kappa=Rational(-1,16),
                        j_type='kramers', adjacency_weights='cg')
D = np.array(st.dirac_operator.tolist(), dtype=float)
evals = np.linalg.eigvalsh(D)
evals_sq = evals**2
dim_H = len(evals)

K_over_pi = 43.619934066848226
B = 42.0

print(f"n_max=4 Kramers+CG: dim_H = {dim_H}")
print(f"  |lambda| range: [{min(abs(evals)):.6f}, {max(abs(evals)):.6f}]")
print(f"  CH max eigenvalue: 4.5")
print(f"  K/pi = {K_over_pi:.4f}, B = {B}")
print(f"  dim_H > K/pi: {dim_H > K_over_pi} (ceiling allows match)")

Lambda = max(abs(evals))
Lambda_sq = Lambda**2
sharp = dim_H
gauss = float(np.sum(np.exp(-evals_sq / Lambda_sq)))
sig01 = float(np.sum(1.0 / (1.0 + np.exp((evals_sq/Lambda_sq - 1.0) / 0.1))))

print(f"\nAt Lambda = max|lambda| = {Lambda:.6f}:")
print(f"  sharp  = {sharp}")
print(f"  gauss  = {gauss:.4f}")
print(f"  sig0.1 = {sig01:.4f}")

# Search for Lambda where various test functions hit K/pi
print(f"\nSearching for Lambda giving K/pi = {K_over_pi:.4f}:")
Ls = np.linspace(0.5, 30, 5000)

for name, func in [
    ("gaussian", lambda L: float(np.sum(np.exp(-evals_sq / L**2)))),
    ("sigmoid_0.1", lambda L: float(np.sum(1.0 / (1.0 + np.exp((evals_sq/L**2 - 1.0) / 0.1))))),
    ("sigmoid_0.3", lambda L: float(np.sum(1.0 / (1.0 + np.exp((evals_sq/L**2 - 1.0) / 0.3))))),
    ("sigmoid_0.5", lambda L: float(np.sum(1.0 / (1.0 + np.exp((evals_sq/L**2 - 1.0) / 0.5))))),
]:
    best_L, best_val, best_err = None, None, 1e10
    for L in Ls:
        v = func(L)
        err = abs(v - K_over_pi)
        if err < best_err:
            best_L, best_val, best_err = L, v, err
    rel = best_err / K_over_pi
    print(f"  {name:15s}: Lambda={best_L:.4f}, value={best_val:.6f}, "
          f"rel_err={rel:.6f} ({rel*100:.3f}%)")

# Also search for B=42
print(f"\nSearching for Lambda giving B = {B:.1f}:")
for name, func in [
    ("gaussian", lambda L: float(np.sum(np.exp(-evals_sq / L**2)))),
    ("sigmoid_0.1", lambda L: float(np.sum(1.0 / (1.0 + np.exp((evals_sq/L**2 - 1.0) / 0.1))))),
]:
    best_L, best_val, best_err = None, None, 1e10
    for L in Ls:
        v = func(L)
        err = abs(v - B)
        if err < best_err:
            best_L, best_val, best_err = L, v, err
    rel = best_err / B
    print(f"  {name:15s}: Lambda={best_L:.4f}, value={best_val:.6f}, "
          f"rel_err={rel:.6f}")

# Sharp cutoff for B=42 would need exactly 42 eigenvalues below Lambda
# Check if any Lambda gives exactly 42
unique_abs = sorted(set(np.round(np.abs(evals), 10)))
cumcount = 0
print(f"\nSharp cutoff staircase (unique |lambda| values):")
for a in unique_abs[:15]:
    c = int(np.sum(np.abs(evals) <= a + 1e-10))
    if cumcount != c:
        print(f"  |lambda| <= {a:.8f}: N = {c}")
        if c >= 42 and cumcount < 42:
            print(f"    *** Jumps past 42: {cumcount} -> {c}, no Lambda gives sharp=42 ***")
        cumcount = c

# Supertrace checks
gamma = np.array(st.grading.tolist(), dtype=float)
gamma_diag = np.diag(gamma)
str_I = float(gamma_diag.sum())
str_D2 = float(np.sum(gamma_diag * evals_sq))
casimir_cont = sum(2*(n+1)*(n+2)*(n+1.5)**2 for n in range(4))

print(f"\nStructural checks:")
print(f"  Str(I) = {str_I:.1f}")
print(f"  Str(D^2) = {str_D2:.4f}  (vs -dim_H = {-dim_H})")
print(f"  Tr(D^2) = {sum(evals_sq):.4f}")
print(f"  Casimir(cont) = {casimir_cont:.4f}")
print(f"  Tr(D^2)/Casimir = {sum(evals_sq)/casimir_cont:.6f}")
