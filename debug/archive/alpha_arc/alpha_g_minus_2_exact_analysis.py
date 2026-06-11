"""
Deeper analysis of the exact forms from alpha_g_minus_2_exact_forms.py.

Questions:
1. What is the exact power-law exponent for F2 vs lambda?
2. Does lam^2 * F2 converge, or is the exponent slightly different from -2?
3. Is there a clean expansion lam^2 * F2 = c0 + c1/lam^2 + ...?
4. Do the V_magnetic expressions have a universal form?
"""

import numpy as np
from sympy import (Rational, sqrt, simplify, S, nsimplify,
                   symbols, radsimp, cancel, together, factor)
import json

# Load the exact data (updated with n_ext=1..7)
with open('debug/data/alpha_g_minus_2_exact_forms.json') as f:
    data = json.load(f)

results = data['results']
print(f"Loaded {len(results)} data points (n_ext=1..{results[-1]['n_ext']})")

print("=" * 70)
print("ANALYSIS OF EXACT ALGEBRAIC FORMS")
print("=" * 70)

# 1. Power-law exponent from successive pairs
print("\n--- 1. Effective power-law exponent ---")
print("  log(F2[i+1]/F2[i]) / log(lam[i+1]/lam[i]) for successive pairs:")
for i in range(len(results) - 1):
    f2_i = results[i]['F2_nint1_float']
    f2_ip1 = results[i + 1]['F2_nint1_float']
    lam_i = results[i]['lambda_float']
    lam_ip1 = results[i + 1]['lambda_float']
    if f2_i > 0 and f2_ip1 > 0:
        eff_exp = np.log(f2_ip1 / f2_i) / np.log(lam_ip1 / lam_i)
        print(f"  n_ext=({results[i]['n_ext']},{results[i+1]['n_ext']}): "
              f"exponent = {eff_exp:.6f}")

# Overall fit
lam_arr = np.array([r['lambda_float'] for r in results])
f2_arr = np.array([r['F2_nint1_float'] for r in results])
log_lam = np.log(lam_arr)
log_f2 = np.log(f2_arr)
coeffs = np.polyfit(log_lam, log_f2, 1)
print(f"\n  Overall log-log fit: F2 ~ lam^{coeffs[0]:.4f}")
print(f"  (coefficient: exp({coeffs[1]:.4f}) = {np.exp(coeffs[1]):.6f})")

# 2. Check if lam^2 * F2 has a convergent expansion
print("\n--- 2. lam^2 * F2 expansion ---")
lam2_f2 = np.array([r['lam2_F2_float'] for r in results])
inv_lam2 = 1.0 / lam_arr**2

print("  Data:")
for i, r in enumerate(results):
    print(f"    n_ext={r['n_ext']}: 1/lam^2={inv_lam2[i]:.6f}, "
          f"lam^2*F2={lam2_f2[i]:.10f}")

# Linear fit: lam^2 * F2 = a + b/lam^2
A = np.column_stack([np.ones_like(inv_lam2), inv_lam2])
fit_coeffs, residuals, _, _ = np.linalg.lstsq(A, lam2_f2, rcond=None)
a, b = fit_coeffs
print(f"\n  Linear fit: lam^2 * F2 = {a:.8f} + {b:.6f}/lam^2")
print(f"  Extrapolated limit (lam->inf): {a:.8f}")

# Try to identify the limit
try:
    limit_guess = nsimplify(a, rational=False, tolerance=1e-4)
    print(f"  nsimplify: {a:.8f} ~ {limit_guess} = {float(limit_guess):.8f}")
except Exception:
    pass

# Quadratic fit: lam^2 * F2 = a + b/lam^2 + c/lam^4
A2 = np.column_stack([np.ones_like(inv_lam2), inv_lam2, inv_lam2**2])
fit_coeffs2, _, _, _ = np.linalg.lstsq(A2, lam2_f2, rcond=None)
a2, b2, c2 = fit_coeffs2
print(f"\n  Quadratic fit: lam^2 * F2 = {a2:.8f} + {b2:.6f}/lam^2 + {c2:.4f}/lam^4")
print(f"  Extrapolated limit (lam->inf): {a2:.8f}")

# Residuals for both fits
print(f"\n  Residuals (linear):")
for i, r in enumerate(results):
    pred = a + b * inv_lam2[i]
    print(f"    n_ext={r['n_ext']}: data={lam2_f2[i]:.10f}, "
          f"fit={pred:.10f}, resid={lam2_f2[i]-pred:.4e}")

print(f"\n  Residuals (quadratic):")
for i, r in enumerate(results):
    pred2 = a2 + b2 * inv_lam2[i] + c2 * inv_lam2[i]**2
    print(f"    n_ext={r['n_ext']}: data={lam2_f2[i]:.10f}, "
          f"fit={pred2:.10f}, resid={lam2_f2[i]-pred2:.4e}")

# 3. V_magnetic structure
print("\n--- 3. V_magnetic structure ---")
print("  Check: V_magnetic = A*sqrt(n*(n+1)) - B*sqrt((n-1)*n) or similar")
for r in results:
    n = r['n_ext']
    v = r['V_magnetic_float']
    lam = r['lambda_float']
    # Try various normalizations
    print(f"  n_ext={n}: V_mag={v:.8f}, "
          f"V_mag*lam={v*lam:.6f}, "
          f"V_mag*lam^2={v*lam**2:.6f}, "
          f"V_mag*(n+1)={v*(n+1):.6f}")

# 4. Alternative: check if the product lam^alpha * F2 is constant for alpha != 2
print("\n--- 4. Find optimal power alpha such that lam^alpha * F2 = const ---")
# Richardson extrapolation using last 3 points
# If F2 ~ C * lam^{-alpha}, then alpha = -log(F2_i/F2_{i+1}) / log(lam_i/lam_{i+1})
# Already done above. Let's try fitting more carefully.
from scipy.optimize import minimize_scalar


def var_of_product(alpha):
    vals = f2_arr * lam_arr**alpha
    return np.var(vals) / np.mean(vals)**2  # coefficient of variation squared


result = minimize_scalar(var_of_product, bounds=(1.5, 3.0), method='bounded')
alpha_opt = result.x
print(f"  Optimal alpha: {alpha_opt:.6f}")
vals_opt = f2_arr * lam_arr**alpha_opt
print(f"  lam^{alpha_opt:.4f} * F2 values:")
for i, r in enumerate(results):
    print(f"    n_ext={r['n_ext']}: {vals_opt[i]:.10f}")
print(f"  Mean: {np.mean(vals_opt):.10f}")
print(f"  CV: {np.std(vals_opt)/np.mean(vals_opt)*100:.2f}%")

# 5. Compute B(n_int=0) for parity check
print("\n--- 5. n_int=0 always zero (parity check) ---")
for r in results:
    print(f"  n_ext={r['n_ext']}: B(n_int=0) = {r['B_nint0_float']}")

# 6. Factor structure of denominators
print("\n--- 6. Denominator analysis ---")
# The denominators of B(n_int=1) have a clear pattern related to
# (2*n_int+3)^4 * q*(q+2) with the CG-induced factors
for r in results:
    n = r['n_ext']
    # From the exact form, extract denominator pattern
    # The internal propagator gives factor 1/(lam_int^4 * mu_q)
    # lam_int = 5/2 for n_int=1, so lam_int^4 = 625/16
    # mu_q = q*(q+2) for allowed q values
    lam_int = Rational(5, 2)
    lam4 = lam_int**4  # = 625/16
    print(f"  n_ext={n}: exact B = {r['B_nint1_exact']}")
    # Check if denominators are multiples of 5^k * 3^j * (2n+3)^m
    # Parse the expression to understand denominator structure

# 7. Check V_magnetic exact against known CG patterns
print("\n--- 7. V_magnetic exact expressions ---")
# For n_ext, j=1/2: the electron has (jL, jR) = ((n+1)/2, n/2)
# The probe photon is q=1, so (jpL, jpR) can be (1, 0) or (0, 1)
# For (1, 0): probe has angular momentum only in L sector
# For (0, 1): probe has angular momentum only in R sector
for r in results:
    n = r['n_ext']
    jL = Rational(n + 1, 2)
    jR = Rational(n, 2)
    print(f"  n_ext={n}: (jL, jR) = ({jL}, {jR}), exact V_mag = {r['V_magnetic_exact']}")

# 8. Rationalize V_magnetic denominators
print("\n--- 8. Rationalized V_magnetic ---")
# sqrt(a)/sqrt(b) patterns
V_mag_exprs = {
    1: -2*sqrt(3)/9 + 2*sqrt(2)/3,
    2: -sqrt(2)/3 + 2*sqrt(15)/9,
    3: -2*sqrt(15)/15 + sqrt(6)/3,
    4: -2*sqrt(6)/9 + 2*sqrt(35)/15,
    5: -2*sqrt(35)/21 + 4*sqrt(3)/9,
}

# Look at the pattern: V_magnetic = a_n * sqrt(prod_n) + b_n * sqrt(prod_n')
# For n=1: 2*sqrt(2)/3 - 2*sqrt(3)/9
# For n=2: -sqrt(2)/3 + 2*sqrt(15)/9
# The arguments under sqrt:
# n=1: sqrt(2), sqrt(3)  -> 1*2, 1*3
# n=2: sqrt(2), sqrt(15) -> 1*2, 3*5
# n=3: sqrt(15), sqrt(6) -> 3*5, 2*3
# n=4: sqrt(6), sqrt(35) -> 2*3, 5*7
# n=5: sqrt(35), sqrt(3) -> 5*7, 1*3
# Pattern in sqrt arguments: n(n+1)/2 and (n+1)(n+2)/2 ?
print("  sqrt arguments and their structure:")
sqrt_args = {
    1: (2, 3),       # = (1*2, 1*3) or (2, 3)
    2: (2, 15),      # = (2, 3*5)
    3: (15, 6),      # = (3*5, 2*3)
    4: (6, 35),      # = (2*3, 5*7)
    5: (35, 3),      # = (5*7, 3) ... hmm
}

for n, (a, b) in sqrt_args.items():
    # Check: is this related to n*(n+1), (n+1)*(n+2) ?
    jL = (n + 1) / 2
    jR = n / 2
    # jL*(jL+1) = (n+1)/2 * (n+3)/2 = (n+1)(n+3)/4
    # jR*(jR+1) = n/2 * (n+2)/2 = n(n+2)/4
    # (2jL+1)*(2jR+1) = (n+2)*(n+1)
    dim_prod = (n + 2) * (n + 1)
    print(f"  n={n}: sqrt({a}, {b}), dim_prod=(n+1)(n+2)={dim_prod}")

# 9. Check the relation between V_magnetic^2 and a rational number
print("\n--- 9. V_magnetic^2 (should rationalize) ---")
for n, expr in V_mag_exprs.items():
    v2 = simplify(expr**2)
    print(f"  n_ext={n}: V_mag^2 = {v2} = {float(v2):.10f}")

# 10. Richardson extrapolation for the limit of lam^2 * F2
print("\n--- 10. Richardson extrapolation ---")
# If lam^2 * F2 = c0 + c1/lam^2 + c2/lam^4 + ...
# Then we can use successive elimination
# Let f(x) = c0 + c1*x + c2*x^2 + ... where x = 1/lam^2
x_vals = [1.0 / r['lambda_float']**2 for r in results]
y_vals = [r['lam2_F2_float'] for r in results]

# Neville's method for Richardson extrapolation
# Build the triangular table
n_pts = len(x_vals)
table = [[0.0] * n_pts for _ in range(n_pts)]
for i in range(n_pts):
    table[i][0] = y_vals[i]

for j in range(1, n_pts):
    for i in range(n_pts - j):
        table[i][j] = (x_vals[i + j] * table[i][j - 1] - x_vals[i] * table[i + 1][j - 1]) / (x_vals[i + j] - x_vals[i])

print("  Neville table (extrapolation to x=0, i.e., lam->inf):")
for j in range(n_pts):
    vals_at_j = [table[i][j] for i in range(n_pts - j)]
    print(f"  Order {j}: {[f'{v:.10f}' for v in vals_at_j]}")

c0_estimate = table[0][n_pts - 1]
print(f"\n  Best estimate of c0 = lim(lam->inf) lam^2 * F2: {c0_estimate:.10f}")

# Try to identify c0
try:
    c0_guess = nsimplify(c0_estimate, rational=False, tolerance=1e-5)
    print(f"  nsimplify: {c0_estimate:.10f} ~ {c0_guess} = {float(c0_guess):.10f}")
except Exception:
    pass

# Try specific targets
targets = {
    '1/180': 1.0/180,
    '1/175': 1.0/175,
    '1/(6*pi^2)': 1.0/(6*np.pi**2),
    '1/(32*pi)': 1.0/(32*np.pi),
    'pi/600': np.pi/600,
    '1/200': 1.0/200,
    '1/192': 1.0/192,
}
print(f"\n  Comparison with targets:")
for name, val in targets.items():
    print(f"    {name} = {val:.10f}, ratio = {c0_estimate/val:.6f}")


print("\n\nDone.")
