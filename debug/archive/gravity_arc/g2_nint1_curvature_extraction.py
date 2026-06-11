"""
g-2 curvature expansion: high-precision extraction from n_int=1 exact forms.

Strategy: We have 7 exact algebraic forms for F2(n_int=1, n_ext) at n_ext=1..7.
Evaluate these at 80+ digit precision using mpmath, then fit the curvature
expansion lam^2 * F2(nint1) = A0 + A1/lam^2 + A2/lam^4 + ...

The n_int=1 term captures the DOMINANT contribution (99% at n_ext=1).
Its curvature expansion structure reveals the geometry.
"""

import json
import sys
import os

sys.path.insert(0, '.')

try:
    import mpmath
    mpmath.mp.dps = 80
except ImportError:
    print("mpmath required")
    sys.exit(1)

with open('debug/data/alpha_g_minus_2_exact_forms.json') as f:
    exact = json.load(f)

alpha_mp = mpmath.mpf('7.2973525693e-3')
schwinger_mp = alpha_mp / (2 * mpmath.pi)

print("=" * 70)
print("  g-2 n_int=1 CURVATURE EXTRACTION (80-digit precision)")
print("=" * 70)

# === 1. Evaluate exact forms at high precision ===
print("\n--- 1. High-precision n_int=1 data ---")
n_ext_vals = []
lam_vals = []
f2_nint1_vals = []
lam2_f2_vals = []

for r in exact['results']:
    n = r['n_ext']
    lam = mpmath.mpf(r['lambda_float'])

    # Reconstruct exact expression from string
    expr_str = r['F2_nint1_exact']

    # Evaluate using mpmath (parse the exact expressions)
    # These involve sqrt(2), sqrt(3), sqrt(5), sqrt(6), sqrt(7), sqrt(10),
    # sqrt(15), sqrt(21), sqrt(35), sqrt(42), sqrt(210)
    f2_val = mpmath.mpf(str(r['F2_nint1_float']))

    # For higher precision, evaluate the lam2_F2 exact form
    lam2_f2_str = r['lam2_F2_exact']
    lam2_f2_val = mpmath.mpf(str(r['lam2_F2_float']))

    n_ext_vals.append(n)
    lam_vals.append(lam)
    f2_nint1_vals.append(f2_val)
    lam2_f2_vals.append(lam2_f2_val)

    ratio = f2_val / schwinger_mp
    lam2_ratio = lam**2 * ratio

    print(f"  n_ext={n}: lam={float(lam):.1f}  "
          f"F2(nint1)/S = {float(ratio):.10f}  "
          f"lam^2*F2/S = {float(lam2_ratio):.10f}")

# === 2. Better approach: use sympy to evaluate EXACTLY ===
print("\n--- 2. Exact symbolic evaluation ---")
try:
    import sympy
    from sympy import sqrt, Rational, nsimplify, N

    # Re-evaluate using the exact algebraic expressions
    # Parse from the JSON strings
    exact_f2_hp = []
    exact_lam2_f2_hp = []

    for r in exact['results']:
        n = r['n_ext']
        lam_exact = sympy.Rational(2*n + 3, 2)

        # Evaluate the lam2_F2 exact expression using sympy
        lam2_f2_expr_str = r['lam2_F2_exact']
        # Clean up the expression for sympy eval
        ctx = {
            'sqrt': sympy.sqrt,
        }
        try:
            lam2_f2_sym = sympy.sympify(lam2_f2_expr_str)
            lam2_f2_hp = mpmath.mpf(str(sympy.N(lam2_f2_sym, 80)))
            exact_lam2_f2_hp.append(lam2_f2_hp)
        except Exception as e:
            print(f"  n_ext={n}: sympy parse failed ({e}), using float")
            exact_lam2_f2_hp.append(mpmath.mpf(str(r['lam2_F2_float'])))

        try:
            f2_sym = sympy.sympify(r['F2_nint1_exact'])
            f2_hp = mpmath.mpf(str(sympy.N(f2_sym, 80)))
            exact_f2_hp.append(f2_hp)
        except Exception as e:
            print(f"  n_ext={n}: F2 parse failed ({e})")
            exact_f2_hp.append(mpmath.mpf(str(r['F2_nint1_float'])))

    print("\n  High-precision lam^2 * F2(nint1) values:")
    for i, n in enumerate(n_ext_vals):
        print(f"    n_ext={n}: {exact_lam2_f2_hp[i]}")

    # Compute lam^2 * F2 / Schwinger at high precision
    lam2_ratio_hp = []
    for i, n in enumerate(n_ext_vals):
        lam = mpmath.mpf(2*n + 3) / 2
        r = lam**2 * exact_f2_hp[i] / schwinger_mp
        lam2_ratio_hp.append(r)
        print(f"    n_ext={n}: lam^2*F2/S = {r}")

except ImportError:
    print("  sympy not available, using float precision")
    lam2_ratio_hp = [mpmath.mpf(str(v)) for v in lam2_f2_vals]
    for i, n in enumerate(n_ext_vals):
        lam = mpmath.mpf(2*n + 3) / 2
        lam2_ratio_hp[i] = lam**2 * f2_nint1_vals[i] / schwinger_mp

# === 3. Polynomial fit in 1/lam^2 ===
print("\n--- 3. Polynomial fit: lam^2*F2/S = A0 + A1*x + A2*x^2 + ... (x = 1/lam^2) ---")

# Build the x values: x_i = 1/lam_i^2
x_hp = []
y_hp = []
for i, n in enumerate(n_ext_vals):
    lam = mpmath.mpf(2*n + 3) / 2
    x_i = 1 / lam**2
    x_hp.append(x_i)
    y_hp.append(lam2_ratio_hp[i])

# Vandermonde matrix for polynomial fit
# y = A0 + A1*x + A2*x^2 + A3*x^3 + A4*x^4 + A5*x^5 + A6*x^6
# With 7 data points, we can fit exactly up to degree 6
n_points = len(x_hp)
for degree in range(2, min(n_points, 7)):
    # Build system
    mat = mpmath.matrix(n_points, degree + 1)
    rhs = mpmath.matrix(n_points, 1)
    for i in range(n_points):
        for j in range(degree + 1):
            mat[i, j] = x_hp[i]**j
        rhs[i] = y_hp[i]

    # Solve least squares (QR)
    try:
        # Use normal equations for simplicity
        A = mat
        b = rhs
        ATA = A.T * A
        ATb = A.T * b
        coeffs = mpmath.lu_solve(ATA, ATb)

        # Compute residual
        residuals = A * coeffs - b
        max_res = max(abs(residuals[i]) for i in range(n_points))

        print(f"\n  Degree {degree} fit (max residual = {float(max_res):.2e}):")
        for j in range(degree + 1):
            c = float(coeffs[j])
            print(f"    A{j} = {c:+.12f}", end="")
            # Try to identify as rational
            if abs(c) > 1e-10:
                frac = mpmath.identify(mpmath.mpf(str(c)), tol=1e-6)
                if frac:
                    print(f"  (~{frac})", end="")
            print()

        # A0 is the flat-space limit
        if degree >= 2:
            A0 = float(coeffs[0])
            A1 = float(coeffs[1])
            print(f"    Flat-space limit A0 = {A0:.10f}")
            print(f"    Leading correction A1/A0 = {A1/A0:.10f}")
    except Exception as e:
        print(f"  Degree {degree}: solve failed ({e})")

# === 4. Richardson extrapolation at high precision ===
print("\n--- 4. Richardson extrapolation (high precision) ---")
n_pts = len(x_hp)
table = [[mpmath.mpf(0)] * n_pts for _ in range(n_pts)]
for i in range(n_pts):
    table[i][0] = y_hp[i]

for j in range(1, n_pts):
    for i in range(n_pts - j):
        num = x_hp[i+j] * table[i][j-1] - x_hp[i] * table[i+1][j-1]
        den = x_hp[i+j] - x_hp[i]
        table[i][j] = num / den

print("  Richardson diagonal:")
for j in range(n_pts):
    val = float(table[0][j])
    print(f"    order {j}: {val:.15f}")

A0_rich = table[0][n_pts - 1]
print(f"\n  Best estimate A0 = {float(A0_rich):.15f}")

# Try to identify A0
print(f"\n  PSLQ identification of A0:")
basis_A0 = [A0_rich, mpmath.mpf(1), mpmath.pi, mpmath.pi**2,
            mpmath.sqrt(2), mpmath.sqrt(3)]
labels_A0 = ["A0", "1", "pi", "pi^2", "sqrt(2)", "sqrt(3)"]
try:
    rel = mpmath.pslq(basis_A0, tol=1e-8, maxcoeff=10000)
    if rel is not None and rel[0] != 0:
        terms = [f"{r}*{l}" for r, l in zip(rel, labels_A0) if r != 0]
        print(f"    {' + '.join(terms)} = 0")
        A0_form = sum(-mpmath.mpf(r)/rel[0] * v for r, l, v in
                      zip(rel[1:], labels_A0[1:],
                          [mpmath.mpf(1), mpmath.pi, mpmath.pi**2,
                           mpmath.sqrt(2), mpmath.sqrt(3)]))
        print(f"    A0 = {float(A0_form):.15f}")
    else:
        print(f"    No relation found")
except Exception as e:
    print(f"    Error: {e}")

# === 5. Curvature expansion of the FULL F2/Schwinger ===
print("\n--- 5. Full F2/Schwinger curvature expansion at n_ext=1 ---")
# At n_ext=1 we have the converged full sum
# F2/Schwinger = 1 + delta
# delta = c1/lam^2 + c2/lam^4 + ...
# where lam = 5/2, so 1/lam^2 = 4/25

with open('debug/data/alpha_g_minus_2_ratio_investigation.json') as f:
    ratio_data = json.load(f)

final_ratio = mpmath.mpf(str(ratio_data['final_ratio']))
delta_full = final_ratio - 1

print(f"  Full delta = {float(delta_full):.12f}")
print(f"  Parker-Toms c1 = R/12 = 1/2 = 0.5")
print(f"  delta * lam^2 = {float(delta_full * mpmath.mpf('25')/4):.12f}")
print(f"  (Should approach R/12 = 0.5 as lam -> inf)")

# From the n_int=1 data, get the n_int=1 contribution to delta at n_ext=1
lam1 = mpmath.mpf('5') / 2
f2_nint1_at_1 = exact_f2_hp[0]
delta_nint1 = f2_nint1_at_1 / schwinger_mp
delta_rest = delta_full - delta_nint1
ratio_nint1 = delta_nint1 / delta_full

print(f"\n  delta from n_int=1: {float(delta_nint1):.12f}")
print(f"  delta from n_int>=2: {float(delta_rest):.12f}")
print(f"  n_int=1 fraction of delta: {float(ratio_nint1)*100:.2f}%")

# === 6. Direct extraction of c1, c2 using two-point method ===
print("\n--- 6. Two-point extraction (n_ext=1 and n_ext=2) ---")
# The n_ext=2 data is NOT well converged for the full sum.
# But we can estimate what it SHOULD be.
#
# Alternative: use the n_int=1 exact data to get the n_int=1 contribution
# to the curvature coefficients, then estimate the n_int>=2 correction.

# n_int=1 curvature expansion from polynomial fit
# F2(nint1,n_ext) / S = A0/lam^2 + A1/lam^4 + ...
# So lam^2 * F2(nint1)/S = A0 + A1/lam^2 + ...

# From the 7-point fit, extract the nint1-contribution curvature coefficients
# The FULL curvature expansion is:
# F2(all)/S = [A0^(sum)]/lam^2 + [A1^(sum)]/lam^4 + ...
# where A0^(sum) = sum over n_int of A0^(nint)
# and we want A0^(sum) = Schwinger (= 1 in F2/S normalization)

# For the n_int=1 part, use degree-4 polynomial
degree = 4
mat = mpmath.matrix(n_points, degree + 1)
rhs = mpmath.matrix(n_points, 1)
for i in range(n_points):
    for j in range(degree + 1):
        mat[i, j] = x_hp[i]**j
    rhs[i] = y_hp[i]
ATA = mat.T * mat
ATb = mat.T * rhs
coeffs_4 = mpmath.lu_solve(ATA, ATb)

print(f"  n_int=1 curvature expansion (degree-4 fit):")
print(f"    lam^2 * F2(nint1)/S = A0 + A1/lam^2 + A2/lam^4 + ...")
for j in range(degree + 1):
    print(f"    A{j} = {float(coeffs_4[j]):+.10f}")

A0_nint1 = float(coeffs_4[0])
A1_nint1 = float(coeffs_4[1])
A2_nint1 = float(coeffs_4[2])

# The curvature expansion of F2(nint1)/S is:
# F2(nint1)/S = A0/lam^2 + A1/lam^4 + A2/lam^6 + ...
# The curvature expansion of F2(all)/S is:
# 1 + c1/lam^2 + c2/lam^4 + ...
# So: sum_nint A0(nint) = 1 (Schwinger normalization)
#     c1 = sum_nint A1(nint)
#     c2 = sum_nint A2(nint)

print(f"\n  Flat-space structure:")
print(f"    A0(nint=1) = {A0_nint1:.10f}")
print(f"    Need: sum A0(nint) = Schwinger * lam^2 / F2")
print(f"    = 1.0 (in F2/S normalization)")
print(f"    Missing from n_int>=2: {1.0 - A0_nint1:.10f}")
print(f"    n_int=1 is {A0_nint1*100:.1f}% of the flat-space limit")

# === 7. Per-level lam^2*F2 analysis ===
print("\n--- 7. Per-level contributions to curvature coefficients ---")
levels = ratio_data['per_level']
V_mag = mpmath.mpf(str(ratio_data['V_magnetic_float']))

# At n_ext=1, the per-level contribution to lam^2 * F2/S is:
# lam^2 * B_level / (V_mag * S)
lam1_sq = mpmath.mpf('25') / 4

print(f"  n_int  lam^2*B/(V*S)    cumulative")
cum = mpmath.mpf(0)
for lev in levels:
    ni = lev['n_int']
    if ni == 0:
        continue
    B = mpmath.mpf(str(lev['B_level']))
    contrib = lam1_sq * B / (V_mag * schwinger_mp)
    cum += contrib
    f2_s_cum = cum / lam1_sq  # F2/S cumulative
    print(f"  {ni:4d}  {float(contrib):14.8f}    {float(cum):14.8f}  "
          f"(cum F2/S = {float(1 + f2_s_cum):.10f})")

# === 8. Asymptotic ratio of consecutive levels ===
print("\n--- 8. Ratio of consecutive B_levels ---")
B_all = [mpmath.mpf(str(lev['B_level'])) for lev in levels if lev['n_int'] >= 1]
n_all = [lev['n_int'] for lev in levels if lev['n_int'] >= 1]

for i in range(1, len(B_all)):
    ratio_b = B_all[i] / B_all[i-1]
    print(f"  B({n_all[i]})/B({n_all[i-1]}) = {float(ratio_b):.6f}  "
          f"(n_ratio = {n_all[i]/n_all[i-1]:.2f})")

# === 9. PSLQ on delta with curvature basis ===
print("\n--- 9. High-precision PSLQ on delta ---")
# Test: delta = (1/2)/lam^2 + c2/lam^4 + c3/lam^6
# Rearranging: delta - (1/2)/lam^2 = c2/lam^4 + c3/lam^6

x1 = mpmath.mpf(4) / 25  # 1/lam^2 at n_ext=1
residual_after_c1 = delta_full - mpmath.mpf(1)/2 * x1

print(f"  delta = {float(delta_full):.15f}")
print(f"  (1/2)/lam^2 = {float(mpmath.mpf(1)/2 * x1):.15f}")
print(f"  residual = {float(residual_after_c1):.15f}")
print(f"  residual / (1/lam^4) = {float(residual_after_c1 / x1**2):.15f}")

# What if c1 is NOT exactly 1/2?
# The Parker-Toms formula for massless Dirac on S^3 may have a different coefficient
# since the mass regulator is replaced by the discrete spectrum.
# Let's be more careful.

# delta * lam^2 = c1 + c2/lam^2 + c3/lam^4 + ...
delta_times_lam2 = delta_full * lam1_sq
print(f"\n  delta * lam^2 = {float(delta_times_lam2):.15f}")
print(f"  (Should be c1 + corrections)")
print(f"  c1 = R/12 = 1/2 gives correction = {float(delta_times_lam2 - 0.5):.15f}")

# PSLQ: is delta_times_lam2 a simple function of known constants?
basis_d = [delta_times_lam2, mpmath.mpf(1), mpmath.pi, mpmath.pi**2,
           mpmath.sqrt(2), mpmath.sqrt(3), mpmath.log(2)]
labels_d = ["d*lam^2", "1", "pi", "pi^2", "sqrt(2)", "sqrt(3)", "log(2)"]

try:
    rel_d = mpmath.pslq(basis_d, tol=1e-4, maxcoeff=1000)
    if rel_d is not None and rel_d[0] != 0:
        terms = [f"{r}*{l}" for r, l in zip(rel_d, labels_d) if r != 0]
        print(f"\n  PSLQ: {' + '.join(terms)} = 0")
    else:
        print(f"\n  PSLQ: no relation found at tol=1e-4")
except Exception as e:
    print(f"\n  PSLQ error: {e}")

# === 10. Summary table ===
print("\n" + "=" * 70)
print("  SUMMARY: n_int=1 curvature expansion")
print("=" * 70)
print(f"""
  n_int=1 contribution (exact algebraic forms, 7 n_ext values):
    Flat-space limit A0 = {A0_nint1:.8f}
    Leading curvature A1 = {A1_nint1:.8f}
    Second curvature A2 = {A2_nint1:.8f}

  Full F2/Schwinger at n_ext=1:
    delta = {float(delta_full):.10f}
    Parker-Toms leading = 2/25 = 0.080000
    Residual = {float(delta_full) - 0.08:.10f}

  Key insight: n_int=1 gives {A0_nint1*100:.0f}% of the flat-space limit.
  The remaining ~{(1-A0_nint1)*100:.0f}% comes from n_int >= 2.
  Curvature coefficients require the FULL sum.
""")
