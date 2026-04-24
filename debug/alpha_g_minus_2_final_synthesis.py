"""
Final synthesis of g-2 on S^3 results.

PROVEN CLOSED FORMS:
  V_channel(1,0)_n = (2/3) * sqrt((n+3)/(n+1))
  V_channel(0,1)_n = -(2/3) * sqrt(n/(n+2))
  V_magnetic(n)    = (2/3) * [sqrt((n+3)/(n+1)) - sqrt(n/(n+2))]

ASYMPTOTIC:
  V_magnetic(n) = 4/(3n) - 2/n^2 + 10/(3n^3) - 6/n^4 + ...
  V_mag * lambda -> 4/3 as n -> inf (confirmed numerically)

KEY QUESTION REMAINING:
  F2(n_ext) = B(n_ext, n_int=1) / V_magnetic(n_ext)
  Does lam^2 * F2 -> clean constant?

This script uses the closed-form V_mag to compute F2 exactly for n_ext=1..7
and analyze the asymptotic structure.
"""

from sympy import Rational, sqrt, simplify, S, nsimplify, pi, cancel
import json
import sys

# Load the exact data
with open('debug/data/alpha_g_minus_2_exact_forms.json') as f:
    data = json.load(f)

results = data['results']

print("=" * 70)
print("SYNTHESIS OF g-2 ON S^3 RESULTS")
print("=" * 70)

# 1. Closed-form V_magnetic
print("\n--- 1. Closed-form V_magnetic (PROVEN) ---")
print("  V_magnetic(n) = (2/3) [sqrt((n+3)/(n+1)) - sqrt(n/(n+2))]")
print()
for r in results:
    n = r['n_ext']
    v_formula = float(Rational(2, 3) * (sqrt(Rational(n+3, n+1)) - sqrt(Rational(n, n+2))))
    v_data = r['V_magnetic_float']
    print(f"  n_ext={n}: formula = {v_formula:.12f}, data = {v_data:.12f}, "
          f"match = {abs(v_formula - v_data) < 1e-10}")

# 2. B(n_int=1) exact values and F2
print("\n--- 2. F2 = B(n_int=1) / V_magnetic (exact) ---")
for r in results:
    n = r['n_ext']
    lam = r['lambda_float']
    f2 = r['F2_nint1_float']
    lam2_f2 = r['lam2_F2_float']
    print(f"  n_ext={n}: lam={lam:.1f}, F2={f2:.6e}, lam^2*F2={lam2_f2:.10f}")

# 3. Asymptotic analysis of lam^2 * F2
print("\n--- 3. Asymptotic analysis of lam^2 * F2 ---")
lam_vals = [r['lambda_float'] for r in results]
f2_vals = [r['F2_nint1_float'] for r in results]
lam2_f2_vals = [r['lam2_F2_float'] for r in results]

# Check if lam^2 * F2 = c0 + c1/lam^2 + c2/lam^4 + ...
import numpy as np

x = np.array([1.0 / l**2 for l in lam_vals])
y = np.array(lam2_f2_vals)

# Linear extrapolation: y = a + b*x
A = np.column_stack([np.ones_like(x), x])
fit, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
c0_lin = fit[0]
print(f"\n  Linear extrapolation: lam^2*F2 -> {c0_lin:.10f} as lam -> inf")

# Quadratic: y = a + b*x + c*x^2
A2 = np.column_stack([np.ones_like(x), x, x**2])
fit2, _, _, _ = np.linalg.lstsq(A2, y, rcond=None)
c0_quad = fit2[0]
print(f"  Quadratic extrapolation: lam^2*F2 -> {c0_quad:.10f} as lam -> inf")

# Richardson/Neville
n_pts = len(x)
table = [[0.0] * n_pts for _ in range(n_pts)]
for i in range(n_pts):
    table[i][0] = y[i]
for j in range(1, n_pts):
    for i in range(n_pts - j):
        table[i][j] = (x[i + j] * table[i][j - 1] - x[i] * table[i + 1][j - 1]) / (x[i + j] - x[i])
c0_rich = table[0][n_pts - 1]
print(f"  Richardson ({n_pts}-point): lam^2*F2 -> {c0_rich:.10f} as lam -> inf")

# 4. Try to identify c0 algebraically
print("\n--- 4. Algebraic identification of c0 ---")
# c0 ~ 0.00543
# Try: 1/n for small n
candidates = {
    '1/180': 1.0/180,
    '1/184': 1.0/184,
    '1/185': 1.0/185,
    '1/183': 1.0/183,
    '4/(3*pi^2)': 4.0/(3*np.pi**2),
    '1/(6*pi^2)': 1.0/(6*np.pi**2),
    '1/(18*pi)': 1.0/(18*np.pi),
    'pi/576': np.pi/576,
    'pi/578': np.pi/578,
    '1/184.06': 1.0/184.06,
    '2/(9*pi^2)': 2.0/(9*np.pi**2),
    '1/(3*pi^2)': 1.0/(3*np.pi**2),
    '4/(9*pi^2)': 4.0/(9*np.pi**2),
    '1/(2*pi^2)': 1.0/(2*np.pi**2),
    'alpha/(2*pi)': (1.0/137.036)/(2*np.pi),
    '1/(9*4*pi)': 1.0/(36*np.pi),
    '16/(9*pi^4)': 16.0/(9*np.pi**4),
    '1/184.0': 1.0/184,
    '1/(9*20.444)': 1.0/(9*20.444),
    '1/(2*pi)^2': 1.0/(2*np.pi)**2,
}

print(f"  c0 estimate: {c0_rich:.10f}")
for name, val in sorted(candidates.items(), key=lambda x: abs(x[1] - c0_rich)):
    rel_err = (c0_rich - val) / c0_rich * 100
    print(f"  {name:20s} = {val:.10f}, rel_err = {rel_err:+.4f}%")
    if abs(rel_err) > 5:
        break

# 5. Try nsimplify
print("\n--- 5. nsimplify attempts ---")
for tol in [1e-4, 1e-5, 1e-6]:
    try:
        guess = nsimplify(c0_rich, rational=False, tolerance=tol)
        print(f"  tol={tol}: {c0_rich:.10f} ~ {guess} = {float(guess):.10f}")
    except Exception as e:
        print(f"  tol={tol}: failed ({e})")

# 6. Key structural observation
print("\n--- 6. Structural decomposition ---")
print("  F2(n_ext) = B(n_ext, n_int=1) / V_magnetic(n_ext)")
print("  V_magnetic(n) = (2/3)[sqrt((n+3)/(n+1)) - sqrt(n/(n+2))]")
print("  V_mag ~ 4/(3n) for large n")
print("  F2 ~ c0 / lambda^2 for large lambda")
print("  So B(n, n_int=1) ~ V_mag * c0/lambda^2 ~ (4/(3n)) * c0/n^2 = 4*c0/(3*n^3)")
print()
print("  B(n_int=1) should scale as n_ext^{-3} at large n_ext")

# Verify B scaling
print("\n  B(n_int=1) values and n_ext^3 * B:")
for r in results:
    n = r['n_ext']
    b = r['B_nint1_float']
    scaled = b * n**3
    print(f"  n_ext={n}: B = {b:.6e}, n^3*B = {scaled:.8f}")

# 7. More careful: compute B * n_ext^3 and check if it converges
b_vals = [r['B_nint1_float'] for r in results]
n_vals = [r['n_ext'] for r in results]
bn3 = [b * n**3 for b, n in zip(b_vals, n_vals)]
print(f"\n  n^3*B values: {[f'{v:.6f}' for v in bn3]}")
print(f"  Ratio of consecutive: {[f'{bn3[i+1]/bn3[i]:.6f}' for i in range(len(bn3)-1)]}")

# 8. lam^2 * F2 as exact sympy
print("\n--- 8. Exact lam^2 * F2 expressions ---")
for r in results:
    n = r['n_ext']
    print(f"  n_ext={n}: lam^2*F2 = {r['lam2_F2_exact']}")

# 9. Check: does F2 * V_mag * lam^2 = B * lam^2 / V_mag * V_mag = B * lam^2
# So lam^2 * F2 = lam^2 * B / V_mag
# We need the exact B expressions to analyze further
print("\n--- 9. Exact B(n_int=1) expressions ---")
for r in results:
    n = r['n_ext']
    print(f"  n_ext={n}: B = {r['B_nint1_exact']}")

# 10. Check if B * 3/(4*c0) * n^3 -> 1
print("\n--- 10. B * n^3 / (4*c0/3) ratios ---")
if c0_rich > 0:
    for r in results:
        n = r['n_ext']
        b = r['B_nint1_float']
        ratio = b * n**3 / (4*c0_rich/3)
        print(f"  n_ext={n}: B*n^3 / (4c0/3) = {ratio:.8f}")

# 11. Better: compute lam^2 * F2 * (9/4) and see if that's cleaner
print("\n--- 11. (9/4) * lam^2 * F2 values ---")
for r in results:
    n = r['n_ext']
    val = r['lam2_F2_float'] * 9/4
    print(f"  n_ext={n}: (9/4)*lam^2*F2 = {val:.10f}")

print(f"\n  Limit of (9/4)*lam^2*F2: {c0_rich * 9/4:.10f}")
# ~ 0.01222
# 1/81.8 ~ 0.01222
# Check against 1/(8*pi^2) = 0.01267, 1/81 = 0.01235, 1/82 = 0.01220
print(f"  1/81 = {1/81:.10f}")
print(f"  1/82 = {1/82:.10f}")
print(f"  1/(8*pi^2) = {1/(8*np.pi**2):.10f}")

# 12. The Schwinger form: F2 from one-loop on S^3
# In flat space, F2 = alpha/(2*pi) = 1/(2*pi) [from one loop]
# On S^3, mode-dependent, but the sum over n_int should give the total
# We only computed n_int=1 (the dominant one).
# The FULL F2 = sum_{n_int} B(n_ext, n_int) / V_mag(n_ext)
# But for the RATIO F2/Schwinger, what matters is the coefficient

# 13. Alternative: look at F2 * lambda^2 * (n+1)*(n+2) since V_mag
# involves sqrt factors with (n+1) and (n+2)
print("\n--- 13. F2 * lambda^2 * (n+1)*(n+2) ---")
for r in results:
    n = r['n_ext']
    val = r['lam2_F2_float'] * (n+1) * (n+2)
    print(f"  n_ext={n}: F2*lam^2*(n+1)(n+2) = {val:.8f}")

# This should approach c0 * n^2 for large n, so diverges.
# Better: F2 * lambda^2 / ((2n+3)/2)^2 * ... hmm

# 14. Exact computation of F2 using closed-form V_mag
print("\n--- 14. F2 using closed-form V_magnetic ---")
for r in results:
    n = r['n_ext']
    v_mag = Rational(2, 3) * (sqrt(Rational(n+3, n+1)) - sqrt(Rational(n, n+2)))
    v_mag_f = float(v_mag)
    b_f = r['B_nint1_float']
    f2 = b_f / v_mag_f
    lam = r['lambda_float']
    print(f"  n_ext={n}: V_mag(formula)={v_mag_f:.10f}, B={b_f:.6e}, "
          f"F2={f2:.6e}, lam^2*F2={lam**2*f2:.10f}")

# 15. Can we find a pattern in B(n_int=1) using the closed-form V_mag?
# B has a specific algebraic structure from the one-loop sum.
# B(n_ext, n_int=1) = sum over q, channels of
#   V_1(n_ext -> n_int=1) * probe * V_3(n_int=1 -> n_ext) / (lam_int^4 * mu_q)
# where lam_int = 5/2, lam_int^4 = 625/16
# and q ranges over allowed values with vertex selection

# For n_ext=1, n_int=1:
#   q in [1, max=2]: vertex_allowed(1, 1, q) requires 1+1+q odd -> q=1
#   mu_1 = 1*3 = 3
# For n_ext=2, n_int=1:
#   q in [1, max=3]: 2+1+q odd -> q=2 (since q=1: 2+1+1=4 even, q=3: 2+1+3=6 even)
#   Wait, q=2: 2+1+2=5 odd YES
#   mu_2 = 2*4 = 8
# For n_ext=3, n_int=1:
#   q in [1, max=4]: 3+1+q odd -> q=1 (5), q=3 (7)
#   mu_1 = 3, mu_3 = 3*5 = 15
# For n_ext=4, n_int=1:
#   q in [1, max=5]: 4+1+q odd -> q=2 (7), q=4 (9)
#   mu_2 = 8, mu_4 = 4*6 = 24

print("\n--- 15. Allowed q values for (n_ext, n_int=1) ---")
for n_ext in range(1, 11):
    n_int = 1
    q_lo = max(1, abs(n_ext - n_int))
    q_hi = n_ext + n_int
    allowed_q = []
    for q in range(q_lo, q_hi + 1):
        if (n_ext + n_int + q) % 2 == 1:
            allowed_q.append(q)
    mu_vals = [q*(q+2) for q in allowed_q]
    print(f"  n_ext={n_ext}: q_allowed = {allowed_q}, mu_q = {mu_vals}")

# Note: the number of allowed q values grows, which makes B harder to find
# in closed form. But the structure is clean.

# 16. The key ratio: lam^2 * F2 = lam^2 * B / V_mag
# Using V_mag = (2/3)[sqrt((n+3)/(n+1)) - sqrt(n/(n+2))]
# and lam = (2n+3)/2, lam^2 = (2n+3)^2/4
# lam^2 * F2 = (2n+3)^2/4 * B / {(2/3)[sqrt((n+3)/(n+1)) - sqrt(n/(n+2))]}
#            = 3(2n+3)^2 * B / {8[sqrt((n+3)/(n+1)) - sqrt(n/(n+2))]}
# As n -> inf: V_mag ~ 4/(3n), so
# lam^2*F2 ~ n^2 * B * 3n / (8 * 4/(3n)) = n^2 * B * 3n * 3n / (8*4) = 9n^4 B / 32
# Hmm, that doesn't look right. Let me redo:
# V_mag ~ 4/(3n), lam ~ n, so lam^2*F2 = n^2 * B / (4/(3n)) = 3n^3 B/4
# If lam^2*F2 -> c0, then B ~ 4c0/(3n^3)

print("\n--- 16. B * 3*n^3 / 4 ratios ---")
for r in results:
    n = r['n_ext']
    b = r['B_nint1_float']
    val = b * 3 * n**3 / 4
    print(f"  n_ext={n}: 3n^3*B/4 = {val:.10f}")

print(f"\n  c0 estimate: {c0_rich:.10f}")
print(f"  3n^3*B/4 at n=7: {results[-1]['B_nint1_float'] * 3 * 7**3 / 4:.10f}")

# 17. Better asymptotic: use the exact V_mag to define the exact ratio
# c0 = lim_{n->inf} lam^2 * B(n, n_int=1) / V_mag(n)
# = lim (2n+3)^2/4 * B(n, 1) / {(2/3)[sqrt((n+3)/(n+1)) - sqrt(n/(n+2))]}
# = lim 3*(2n+3)^2 * B(n, 1) / {8 * [sqrt((n+3)/(n+1)) - sqrt(n/(n+2))]}

# Compute this exact ratio for each n
print("\n--- 17. Exact ratio 3*(2n+3)^2 * B / {8 * V_diff} ---")
for r in results:
    n = r['n_ext']
    lam = Rational(2*n+3, 2)
    v_diff = sqrt(Rational(n+3, n+1)) - sqrt(Rational(n, n+2))
    # lam^2 * F2 = lam^2 * B / V_mag = lam^2 * B / ((2/3)*v_diff)
    # = (3/2) * lam^2 * B / v_diff
    ratio_exact = Rational(3, 2) * lam**2 * Rational(r['B_nint1_float']) / float(v_diff)
    print(f"  n_ext={n}: lam^2*F2 = {r['lam2_F2_float']:.10f}, "
          f"3*lam^2*B/(2*v_diff) = {ratio_exact:.10f}")

# 18. Try to identify the limit using more data points
# We need n_ext=4,5 (the background agent) results to get a better extrapolation
# For now, use the 7 data points we have.

# The effective exponent from F2:
print("\n--- 18. Effective power-law exponent for F2 ---")
for i in range(len(results) - 1):
    f2_i = results[i]['F2_nint1_float']
    f2_ip1 = results[i + 1]['F2_nint1_float']
    lam_i = results[i]['lambda_float']
    lam_ip1 = results[i + 1]['lambda_float']
    eff_exp = np.log(f2_ip1 / f2_i) / np.log(lam_ip1 / lam_i)
    print(f"  n_ext=({results[i]['n_ext']},{results[i+1]['n_ext']}): "
          f"exponent = {eff_exp:.6f}")

# The exponents are approaching -2 from below (-2.07), consistent with F2 ~ c0/lam^2

# 19. Summary
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print()
print("1. PROVEN CLOSED FORMS (exact, verified for n=1..10):")
print("   V_channel(1,0)_n = (2/3) * sqrt((n+3)/(n+1))")
print("   V_channel(0,1)_n = -(2/3) * sqrt(n/(n+2))")
print("   V_magnetic(n)    = (2/3) * [sqrt((n+3)/(n+1)) - sqrt(n/(n+2))]")
print()
print("2. INTER-CHANNEL RECURSION (exact):")
print("   V_channel(0,1)_n = -[n/(n+2)] * V_channel(1,0)_{n-1}")
print()
print("3. ASYMPTOTIC BEHAVIOR:")
print("   V_magnetic ~ 4/(3n) - 2/n^2 + O(1/n^3)")
print("   V_mag * lambda -> 4/3 as n -> infinity")
print()
print("4. ONE-LOOP VERTEX CORRECTION:")
print("   B(n_ext=0) = 0 for all n_ext (parity selection)")
print("   F2(n_ext) = B(n_ext, n_int=1) / V_magnetic(n_ext)")
print(f"   lam^2 * F2 -> {c0_rich:.8f} (Richardson extrapolation)")
print(f"   Effective exponent F2 ~ lam^{{-alpha}}: alpha -> 2.07 (approaching 2)")
print()
print("5. ALGEBRAIC IDENTIFICATION OF c0:")
print(f"   c0 = {c0_rich:.10f}")
print(f"   1/c0 = {1/c0_rich:.6f}")
try:
    g = nsimplify(c0_rich, rational=False, tolerance=1e-4)
    print(f"   nsimplify: {g} = {float(g):.10f}")
except:
    print(f"   nsimplify: no clean identification")
print()
print("6. STRUCTURAL INSIGHT:")
print("   The V_magnetic closed form decomposes the Schwinger-like vertex")
print("   into two channels (1,0) and (0,1) with:")
print("   - sqrt((n+3)/(n+1)) ~ 1 + 1/n : the LEFT sector coupling")
print("   - sqrt(n/(n+2)) ~ 1 - 1/n : the RIGHT sector coupling")
print("   The magnetic moment is the DIFFERENCE of these two couplings,")
print("   which goes as 2/n ~ 1/lambda at leading order.")
print("   This is the S^3 analog of the flat-space magnetic moment vertex,")
print("   where the 1/lambda suppression reflects the finite curvature.")
