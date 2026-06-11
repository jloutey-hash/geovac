"""
g-2 curvature analysis: extract the curvature expansion coefficients
from the S^3 vertex correction data.

Key question: F2(n_ext) / Schwinger = 1 + c1/lam^2 + c2/lam^4 + ...
What are c1, c2? Are they clean algebraic numbers?

Also: Can we identify F2(n_ext=1) itself (the converged sum) via PSLQ
against a purely algebraic basis (since F2 is algebraic, not transcendental)?
"""

import json
import sys
import numpy as np

sys.path.insert(0, '.')

# Load the existing data
with open('debug/data/alpha_g_minus_2_ratio_investigation.json') as f:
    ratio_data = json.load(f)

with open('debug/data/alpha_g_minus_2_exact_forms.json') as f:
    exact_data = json.load(f)

print("=" * 70)
print("  g-2 ON S^3: CURVATURE ANALYSIS")
print("=" * 70)

# === Part 1: Full sum convergence at n_ext=1 ===
print("\n--- 1. Convergence of F2/Schwinger at n_ext=1 ---")
levels = ratio_data['per_level']
schwinger = ratio_data['schwinger']
V_mag = ratio_data['V_magnetic_float']

for lev in levels:
    n = lev['n_int']
    if n == 0:
        continue
    r = lev['F2_over_schwinger']
    b = lev['B_level']
    print(f"  n_int={n:2d}: B_level={b:12.4e}  F2/S = {r:.10f}")

# Richardson extrapolation on F2/Schwinger
print("\n--- 2. Richardson extrapolation on F2/Schwinger ---")
# Use the last 6 points for extrapolation
ratios = [lev['F2_over_schwinger'] for lev in levels if lev['n_int'] >= 1]
n_vals = [lev['n_int'] for lev in levels if lev['n_int'] >= 1]

# The convergence rate: each pair of levels drops by ~factor 10
# Power law: B_level ~ A * n^{-p}
# Fit p from log-log of B_level vs n_int
b_levels = [lev['B_level'] for lev in levels if lev['n_int'] >= 2]
n_b = [lev['n_int'] for lev in levels if lev['n_int'] >= 2]
log_n = np.log(np.array(n_b, dtype=float))
log_b = np.log(np.array(b_levels))
# Fit log_b = a - p * log_n
A = np.column_stack([np.ones_like(log_n), -log_n])
fit, _, _, _ = np.linalg.lstsq(A, log_b, rcond=None)
p_fit = fit[1]
print(f"  Power-law exponent of B_level decay: p = {p_fit:.2f}")

# Aitken's delta-squared for the ratio sequence
print("\n  Aitken's delta-squared acceleration:")
for i in range(len(ratios) - 2):
    s0, s1, s2 = ratios[i], ratios[i+1], ratios[i+2]
    denom = s2 - 2*s1 + s0
    if abs(denom) > 1e-20:
        s_aitken = s0 - (s1 - s0)**2 / denom
        print(f"    n_int={n_vals[i]:2d},{n_vals[i+1]:2d},{n_vals[i+2]:2d}: "
              f"Aitken = {s_aitken:.12f}")

# Wynn epsilon acceleration
print("\n  Wynn epsilon (last 8 points):")
eps_data = ratios[-8:]
n_eps = len(eps_data)
# Build the epsilon table
eps = [[0.0] * (n_eps + 1) for _ in range(n_eps)]
for i in range(n_eps):
    eps[i][0] = 0.0
    eps[i][1] = eps_data[i]
for j in range(2, n_eps + 1):
    for i in range(n_eps - j + 1):
        diff = eps[i+1][j-1] - eps[i][j-1]
        if abs(diff) > 1e-30:
            eps[i][j] = eps[i+1][j-2] + 1.0 / diff
        else:
            eps[i][j] = 1e30
# Even columns are the accelerated estimates
for j in range(2, min(n_eps + 1, 9), 2):
    print(f"    eps[0][{j}] = {eps[0][j]:.12f}")

# Simple tail estimate: sum remaining B_level terms
# B_level ~ C * n^{-p} with p ~ 7
C = np.exp(fit[0])
tail = 0.0
for n in range(13, 10000):
    tail += C * n**(-p_fit)
B_tail = tail
F2_tail = B_tail / V_mag
ratio_tail = F2_tail / schwinger
final_ratio = ratio_data['final_ratio']
print(f"\n  Tail correction (n_int=13..10000): {ratio_tail:.2e}")
print(f"  Extrapolated F2/Schwinger: {final_ratio + ratio_tail:.10f}")

# === Part 2: PSLQ on F2 itself (algebraic basis) ===
print("\n--- 3. PSLQ on F2(n_ext=1) against algebraic basis ---")
try:
    import mpmath
    mpmath.mp.dps = 50

    F2_val = mpmath.mpf(str(ratio_data['final_F2']))

    # F2 should be algebraic -- try a basis of sqrt combinations
    sqrt2 = mpmath.sqrt(2)
    sqrt3 = mpmath.sqrt(3)
    sqrt5 = mpmath.sqrt(5)
    sqrt6 = mpmath.sqrt(6)
    sqrt7 = mpmath.sqrt(7)
    sqrt10 = mpmath.sqrt(10)
    sqrt15 = mpmath.sqrt(15)

    # Basic algebraic basis
    basis = [F2_val, mpmath.mpf(1), sqrt2, sqrt3, sqrt6]
    labels = ["F2", "1", "sqrt(2)", "sqrt(3)", "sqrt(6)"]

    try:
        rel = mpmath.pslq(basis, tol=1e-15, maxcoeff=1000000)
        if rel is not None and rel[0] != 0:
            print(f"  PSLQ hit: {' + '.join(f'{r}*{l}' for r, l in zip(rel, labels) if r != 0)}")
        else:
            print(f"  PSLQ: no relation found in basic sqrt basis")
    except Exception as e:
        print(f"  PSLQ error: {e}")

    # Try with more surds (from the exact n_int=1 form)
    basis2 = [F2_val, mpmath.mpf(1), sqrt2, sqrt3, sqrt5, sqrt6, sqrt7,
              sqrt10, sqrt15, mpmath.sqrt(21), mpmath.sqrt(42)]
    labels2 = ["F2", "1", "rt2", "rt3", "rt5", "rt6", "rt7", "rt10", "rt15", "rt21", "rt42"]

    try:
        rel2 = mpmath.pslq(basis2, tol=1e-12, maxcoeff=100000)
        if rel2 is not None and rel2[0] != 0:
            print(f"  Extended PSLQ: ", end="")
            terms = [f"{r}*{l}" for r, l in zip(rel2, labels2) if r != 0]
            print(" + ".join(terms))
        else:
            print(f"  Extended PSLQ: no relation found")
    except Exception as e:
        print(f"  Extended PSLQ error: {e}")

except ImportError:
    print("  mpmath not available")

# === Part 3: Curvature expansion from n_ext dependence ===
print("\n--- 4. Curvature expansion from n_ext dependence ---")
print("  Using n_int=1 dominant term (proven closed forms)")

# Data from exact_forms
for r in exact_data['results']:
    n = r['n_ext']
    lam = r['lambda_float']
    f2 = r['F2_nint1_float']
    lam2_f2 = r['lam2_F2_float']
    print(f"  n_ext={n}: lam={lam:.1f}  F2(n_int=1)={f2:.6e}  lam^2*F2={lam2_f2:.8f}")

# Fit: F2(lam) = a0/lam^2 + a1/lam^4 + a2/lam^6
lam_vals = np.array([r['lambda_float'] for r in exact_data['results']])
f2_vals = np.array([r['F2_nint1_float'] for r in exact_data['results']])

# Working variable: x = 1/lam^2
x = 1.0 / lam_vals**2
y = f2_vals

# Fit y = a0*x + a1*x^2 + a2*x^3 (no constant -- F2 -> 0 as lam -> inf)
X = np.column_stack([x, x**2, x**3])
coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
a0, a1, a2 = coeffs

print(f"\n  Fit: F2(lam) = {a0:.6f}/lam^2 + {a1:.6f}/lam^4 + {a2:.6f}/lam^6")
print(f"  a0 = {a0:.10f}")
print(f"  a0 * 2pi = {a0 * 2 * np.pi:.10f}")

# Check: does a0 = Schwinger = alpha/(2*pi)?
print(f"\n  Schwinger = alpha/(2pi) = {schwinger:.10f}")
print(f"  a0 = {a0:.10f}")
print(f"  a0 / Schwinger = {a0 / schwinger:.6f}")
print(f"  If a0 = alpha/(2pi), then Schwinger IS the leading flat-space limit.")

# Check: what is a1/a0?
print(f"\n  a1/a0 = {a1/a0:.6f}")
print(f"  This is the curvature correction coefficient c1")
print(f"  Parker-Toms predicts R_scalar/12 = 6/12 = 1/2 = 0.5")
print(f"  Measured c1 = {a1/a0:.6f}")

# More careful: fit with just a0 and a1 (2-parameter)
X2 = np.column_stack([x, x**2])
coeffs2, _, _, _ = np.linalg.lstsq(X2, y, rcond=None)
a0_2, a1_2 = coeffs2
print(f"\n  2-parameter fit: F2 = {a0_2:.8f}/lam^2 + {a1_2:.8f}/lam^4")
print(f"  a0_2 = {a0_2:.10f}")
print(f"  a0_2 / Schwinger = {a0_2 / schwinger:.6f}")
print(f"  c1 = a1/a0 = {a1_2/a0_2:.6f}")

# Use Neville/Richardson to extrapolate lam^2*F2 -> a0
print("\n--- 5. Richardson extrapolation: lam^2*F2 -> a0 as lam->inf ---")
lam2_f2 = lam_vals**2 * f2_vals
# Neville table in 1/lam^2
x_rich = 1.0 / lam_vals**2
n_pts = len(x_rich)
table = np.zeros((n_pts, n_pts))
table[:, 0] = lam2_f2
for j in range(1, n_pts):
    for i in range(n_pts - j):
        table[i, j] = (x_rich[i+j] * table[i, j-1] - x_rich[i] * table[i+1, j-1]) / (x_rich[i+j] - x_rich[i])
print(f"  Richardson table diagonal:")
for j in range(n_pts):
    print(f"    order {j}: a0 ~ {table[0, j]:.12f}")

a0_rich = table[0, n_pts - 1]
print(f"\n  Best Richardson estimate: a0 = {a0_rich:.12f}")
print(f"  a0 / Schwinger = {a0_rich / schwinger:.8f}")

# === Part 4: PSLQ on a0 ===
print("\n--- 6. PSLQ on a0 (the flat-space limit) ---")
try:
    a0_mp = mpmath.mpf(str(a0_rich))
    pi_mp = mpmath.pi
    alpha = mpmath.mpf('7.2973525693e-3')  # CODATA 2018
    schwinger_mp = alpha / (2 * pi_mp)

    # Is a0 = alpha/(2*pi)?
    diff = abs(a0_mp - schwinger_mp)
    print(f"  a0 = {a0_mp}")
    print(f"  alpha/(2pi) = {schwinger_mp}")
    print(f"  |a0 - alpha/(2pi)| = {diff}")
    print(f"  Relative diff = {float(diff / schwinger_mp):.6e}")

    # PSLQ: is a0 a rational multiple of alpha/(2pi)?
    basis3 = [a0_mp, schwinger_mp, mpmath.mpf(1)]
    try:
        rel3 = mpmath.pslq(basis3, tol=1e-8, maxcoeff=10000)
        if rel3 is not None:
            print(f"  PSLQ {rel3}: {rel3[0]}*a0 + {rel3[1]}*Schwinger + {rel3[2]} = 0")
            if rel3[0] != 0 and rel3[2] == 0:
                ratio = -mpmath.mpf(rel3[1]) / rel3[0]
                print(f"  => a0 = {ratio} x Schwinger")
        else:
            print(f"  PSLQ: no integer relation found")
    except Exception as e:
        print(f"  PSLQ error: {e}")

except Exception as e:
    print(f"  Error: {e}")

# === Part 5: Summary ===
print("\n" + "=" * 70)
print("  SUMMARY")
print("=" * 70)
print(f"""
  At n_ext=1 (ground state external electron on unit S^3):
    F2 / Schwinger = {ratio_data['final_ratio']:.8f}  (converged at n_int=12)
    Deviation from flat space: {(ratio_data['final_ratio'] - 1)*100:.2f}%

  Parker-Toms leading curvature correction:
    R/(12*lam^2) = 6/(12 * 25/4) = 2/25 = 0.08 = {0.08*100:.0f}%
    Accounts for {0.08/0.08445*100:.0f}% of total deviation

  n_ext -> inf limit (flat-space recovery):
    a0 (Richardson) = {a0_rich:.10f}
    Schwinger       = {schwinger:.10f}
    Ratio           = {a0_rich/schwinger:.8f}
""")
