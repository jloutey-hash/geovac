"""
Investigate the strikingly linear CG_sum / V_mag^2 ratio.

From the previous run:
  n=1: CG_sum/V_mag^2 = 0.2620
  n=2: CG_sum/V_mag^2 = 0.4218
  n=3: CG_sum/V_mag^2 = 0.5809
  n=4: CG_sum/V_mag^2 = 0.7396
  n=5: CG_sum/V_mag^2 = 0.8981
  n=6: CG_sum/V_mag^2 = 1.0564
  n=7: CG_sum/V_mag^2 = 1.2146

Differences ~0.1584-0.1587, almost exactly constant.
Linear fit: ratio ~ a + b*n

Also: CG_sum / V_mag approaches ~0.2 (slowly growing).
And: CG_sum * (n+1)*(n+2) appears EXACTLY linear in n!
"""

from sympy import (Rational, sqrt, simplify, S, sympify, nsimplify,
                   cancel, factor, collect, radsimp, together)
from sympy.physics.wigner import clebsch_gordan
import json

# Load exact results
with open('debug/data/alpha_g_minus_2_exact_forms.json') as f:
    data = json.load(f)
results = data['results']


def V_magnetic_closed(n):
    return Rational(2, 3) * (sqrt(Rational(n+3, n+1)) - sqrt(Rational(n, n+2)))


def V10(n):
    return Rational(2, 3) * sqrt(Rational(n+3, n+1))


def V01(n):
    return -Rational(2, 3) * sqrt(Rational(n, n+2))


# 1. Analyze the CG_sum / V_mag^2 ratio more carefully
print("=== CG_sum / V_mag^2: checking linearity ===")
ratios = []
for r in results:
    n = r['n_ext']
    b_exact = sympify(r['B_nint1_exact'])
    denom_exact = Rational(625, 16) * n * (n + 2)
    cg_exact = simplify(b_exact * denom_exact)
    v_mag = V_magnetic_closed(n)
    v_mag_sq = simplify(v_mag**2)
    ratio = simplify(cg_exact / v_mag_sq)
    ratio_f = float(ratio)
    ratios.append((n, ratio_f, ratio))
    print(f"  n={n}: ratio = {ratio_f:.12f}")

print("\n  Successive differences:")
for i in range(1, len(ratios)):
    diff = ratios[i][1] - ratios[i-1][1]
    print(f"  n={ratios[i][0]}: diff = {diff:.12f}")

# Linear fit: r = a + b*n
# Using n=6,7: b = (r7 - r6) = diff, a = r7 - 7*b
b_est = ratios[-1][1] - ratios[-2][1]
a_est = ratios[-1][1] - 7 * b_est
print(f"\n  Linear fit: a = {a_est:.12f}, b = {b_est:.12f}")
print(f"  Check: 9*a/b = {9*a_est/b_est:.10f}")
print(f"  b = {b_est:.12f}")
print(f"  9*b = {9*b_est:.12f}")
print(f"  1/(9*b) = {1/(9*b_est):.12f}")

# Try to identify b as a rational number
from fractions import Fraction
b_frac = Fraction(b_est).limit_denominator(10000)
print(f"  b as fraction: {b_frac} = {float(b_frac):.12f}")
a_frac = Fraction(a_est).limit_denominator(10000)
print(f"  a as fraction: {a_frac} = {float(a_frac):.12f}")

# Check fit quality
print("\n  Linear fit residuals:")
for n, ratio_f, _ in ratios:
    pred = a_est + b_est * n
    resid = ratio_f - pred
    print(f"  n={n}: actual={ratio_f:.12f}, pred={pred:.12f}, resid={resid:.2e}")


# 2. CG_sum * (n+1)*(n+2) also looks linear - check
print("\n\n=== CG_sum * (n+1)*(n+2): checking linearity ===")
prod_vals = []
for r in results:
    n = r['n_ext']
    b_exact = sympify(r['B_nint1_exact'])
    denom_exact = Rational(625, 16) * n * (n + 2)
    cg_exact = simplify(b_exact * denom_exact)
    val = float(cg_exact) * (n+1) * (n+2)
    prod_vals.append((n, val))
    print(f"  n={n}: CG_sum*(n+1)*(n+2) = {val:.12f}")

print("\n  Successive differences:")
for i in range(1, len(prod_vals)):
    diff = prod_vals[i][1] - prod_vals[i-1][1]
    print(f"  n={prod_vals[i][0]}: diff = {diff:.12f}")

# Check: differences are ~0.280. Is 0.280 = 7/25 = 0.28? Or 9/32 = 0.28125?
diffs_prod = [prod_vals[i][1] - prod_vals[i-1][1] for i in range(1, len(prod_vals))]
avg_diff = sum(diffs_prod) / len(diffs_prod)
print(f"\n  Average diff = {avg_diff:.12f}")
print(f"  7/25 = {7/25:.12f}")
print(f"  9/32 = {9/32:.12f}")
print(f"  2/7  = {2/7:.12f}")
print(f"  14/50 = {14/50:.12f}")

diff_frac = Fraction(avg_diff).limit_denominator(1000)
print(f"  Best fraction for avg_diff: {diff_frac} = {float(diff_frac):.12f}")


# 3. Let's look at the EXACT sympy ratio CG_sum / V_mag^2
# V_mag^2 = (4/9)[(n+3)/(n+1) + n/(n+2) - 2*sqrt(n(n+3)/((n+1)(n+2)))]
# If CG_sum / V_mag^2 is close to a + b*n, then:
# CG_sum ~ (a + b*n) * V_mag^2
# But V_mag^2 has an irrational cross-term. So (a+b*n) can't be exactly rational.
# Unless CG_sum itself has the same irrational structure...

print("\n\n=== Exact sympy: CG_sum / V_mag^2 ===")
for r in results:
    n = r['n_ext']
    b_exact = sympify(r['B_nint1_exact'])
    denom_exact = Rational(625, 16) * n * (n + 2)
    cg_exact = simplify(b_exact * denom_exact)
    v_mag = V_magnetic_closed(n)
    v_mag_sq = simplify(v_mag**2)

    # Rationalize the ratio
    ratio = simplify(cg_exact / v_mag_sq)
    ratio_r = radsimp(ratio)
    print(f"  n={n}: CG_sum/V_mag^2 = {ratio_r}")


# 4. Instead of V_mag^2, let's try V_mag * something rational
# V_mag = (2/3)[sqrt((n+3)/(n+1)) - sqrt(n/(n+2))]
# V_mag * sqrt((n+1)(n+2)) = (2/3)[sqrt((n+3)(n+2)) - sqrt(n(n+1))]
# Let f(n) = sqrt(n(n+1)). Then:
# V_mag * sqrt((n+1)(n+2)) = (2/3)[f(n+2) - f(n)]
# This is a second-difference of f!

print("\n\n=== V_mag * sqrt((n+1)*(n+2)) ===")
for n in range(1, 11):
    v = V_magnetic_closed(n)
    factor_val = sqrt(Rational((n+1)*(n+2)))
    product = simplify(v * factor_val)
    print(f"  n={n}: V_mag*sqrt((n+1)(n+2)) = {product} = {float(product):.10f}")
    # Check: is this (2/3)*(sqrt((n+2)(n+3)) - sqrt(n(n+1)))?
    expected = Rational(2, 3) * (sqrt(Rational((n+2)*(n+3))) - sqrt(Rational(n*(n+1))))
    diff = simplify(product - expected)
    print(f"         check: diff = {diff}")


# 5. Let's try: CG_sum * sqrt((n+1)*(n+2)) -- might clean up
print("\n\n=== CG_sum * sqrt((n+1)*(n+2)) ===")
for r in results:
    n = r['n_ext']
    b_exact = sympify(r['B_nint1_exact'])
    denom_exact = Rational(625, 16) * n * (n + 2)
    cg_exact = simplify(b_exact * denom_exact)
    factor_val = sqrt(Rational((n+1)*(n+2)))
    product = simplify(cg_exact * factor_val)
    print(f"  n={n}: {product} = {float(product):.12f}")


# 6. What about F2 = B / V_mag = CG_sum / (denom * V_mag)?
# F2 * denom = CG_sum / V_mag
# This approached ~0.19 and is growing slowly
# F2 * lam^2 -> c0

# F2 * lam^2 = CG_sum / (denom * V_mag) * lam^2
# = CG_sum * lam^2 / (625/16 * n*(n+2) * V_mag)
# ~ CG_sum * n^2 / (625/16 * n^3 * 4/(3n))
# = CG_sum * n^2 * 3n / (625/16 * n^3 * 4)
# = CG_sum * 3 / (625/16 * 4)
# = CG_sum * 3*16 / (625*4)
# = CG_sum * 12/625
# Wait, this doesn't converge because CG_sum ~ 1/n.

# Let me recalculate:
# V_mag ~ 4/(3n)
# denom = 625/16 * n*(n+2) ~ 625/16 * n^2
# F2 = CG_sum / (denom * V_mag) ~ CG_sum * 16 / (625 * n^2 * 4/(3n))
# = CG_sum * 16 * 3n / (625 * n^2 * 4)
# = CG_sum * 12 / (625 * n)
# F2 * lam^2 ~ F2 * n^2 = CG_sum * 12n / 625

# So c0 = lim CG_sum * 12n / 625
# CG_sum * n grows but slowly...

print("\n\n=== CG_sum * 12*n / 625 (should approach c0) ===")
for r in results:
    n = r['n_ext']
    b_exact = sympify(r['B_nint1_exact'])
    denom_exact = Rational(625, 16) * n * (n + 2)
    cg_exact = simplify(b_exact * denom_exact)
    val = float(cg_exact) * 12 * n / 625
    print(f"  n={n}: {val:.12f}")
print(f"  c0 (Richardson) = 0.005431")

# CG_sum * n is NOT converging fast. Let's try CG_sum * n * (n+2) / (n+1)?
# Or other power corrections.

print("\n\n=== CG_sum * n with power corrections ===")
for r in results:
    n = r['n_ext']
    b_exact = sympify(r['B_nint1_exact'])
    denom_exact = Rational(625, 16) * n * (n + 2)
    cg_exact = simplify(b_exact * denom_exact)
    cg_f = float(cg_exact)
    v1 = cg_f * n
    v2 = cg_f * n * (1 + 1/(2*n))
    v3 = cg_f * n * (n+1)/(n+0.5)
    v4 = cg_f * (n + 0.5)
    print(f"  n={n}: cg*n={v1:.8f}, cg*(n+0.5)={v4:.8f}, "
          f"cg*n*(1+1/2n)={v2:.8f}, cg*n*(n+1)/(n+.5)={v3:.8f}")


# 7. KEY INSIGHT: CG_sum * (n+1)*(n+2) is nearly exactly linear: a + b*n
# Let's check if it's exactly: CG_sum * (n+1)*(n+2) = p + q*n for rational p, q
# This would mean CG_sum = (p + q*n) / ((n+1)*(n+2))

# If CG_sum*(n+1)*(n+2) = p + q*n:
# At n=6: p + 6q = 1.8865654807
# At n=7: p + 7q = 2.1668939812
# => q = 0.2803285005
# => p = 1.8865654807 - 6*0.2803285005 = 0.2045944780

print("\n\n=== CG_sum*(n+1)*(n+2) = p + q*n test ===")
q_val = prod_vals[-1][1] - prod_vals[-2][1]
p_val = prod_vals[-1][1] - 7 * q_val
print(f"  From n=6,7: p = {p_val:.12f}, q = {q_val:.12f}")

# Check these
for n, val in prod_vals:
    pred = p_val + q_val * n
    resid = val - pred
    print(f"  n={n}: actual={val:.12f}, pred={pred:.12f}, resid={resid:.6e}")

# The residuals are O(1/n), not exactly linear. But let's try
# CG_sum * (n+1)*(n+2) = p + q*n + r/n for some r
# Or: CG_sum * (n+1)*(n+2) * n = p*n + q*n^2 + r
# Using n=5,6,7:
# 5*(p+5q+r/5) = 5p + 25q + r
# 6*(p+6q+r/6) = 6p + 36q + r
# 7*(p+7q+r/7) = 7p + 49q + r

v5 = prod_vals[4][1]  # n=5
v6 = prod_vals[5][1]  # n=6
v7 = prod_vals[6][1]  # n=7

# n*v(n) = p*n + q*n^2 + r
# 5*v5 = 5p + 25q + r
# 6*v6 = 6p + 36q + r
# 7*v7 = 7p + 49q + r

# From differences:
# 6*v6 - 5*v5 = p + 11q
# 7*v7 - 6*v6 = p + 13q
# => 2q = (7*v7-6*v6) - (6*v6-5*v5) = 7*v7 - 12*v6 + 5*v5
q2 = (7*v7 - 12*v6 + 5*v5) / 2
p2 = (6*v6 - 5*v5) - 11*q2
r2 = 5*v5 - 5*p2 - 25*q2

print(f"\n  Three-parameter fit CG*(n+1)(n+2) = p + q*n + r/n:")
print(f"  p = {p2:.12f}, q = {q2:.12f}, r = {r2:.12f}")

for n, val in prod_vals:
    pred = p2 + q2*n + r2/n
    resid = val - pred
    print(f"  n={n}: actual={val:.12f}, pred={pred:.12f}, resid={resid:.6e}")


# 8. Try nsimplify on the key parameters
print("\n\n=== nsimplify on q (slope of CG*(n+1)(n+2)) ===")
from sympy import nsimplify as ns, pi, Rational as R
q_sym = ns(q2, rational=False, tolerance=1e-6)
print(f"  q = {q2:.12f}")
print(f"  nsimplify(q) = {q_sym} = {float(q_sym):.12f}")

# Try specific candidates for q
candidates = [R(7,25), R(9,32), R(14,50), R(11,39), R(28,100),
              R(7,25), R(31,110), R(14,49), R(2,7), R(4,14+1)]
for c in candidates:
    err = abs(float(c) - q2) / q2
    print(f"  {c} = {float(c):.12f}, rel err = {err:.6e}")

# Also check q against 4/9 * 7/10 = 28/90 = 14/45
print(f"  14/45 = {14/45:.12f}, rel err = {abs(14/45 - q2)/q2:.6e}")
print(f"  4/9 * V_mag(1) = {4/9 * float(V_magnetic_closed(1)):.12f}")
print(f"  4/(9*sqrt(6)) = {4/(9*6**0.5):.12f}")


# 9. Let's try the exact sympy approach for n=6 where CG_sum is simple
print("\n\n=== n=6 deep dive (simplest CG_sum) ===")
r6 = results[5]  # n_ext=6
b6_exact = sympify(r6['B_nint1_exact'])
denom6 = Rational(625, 16) * 6 * 8
cg6 = simplify(b6_exact * denom6)
print(f"  CG_sum(6) = {cg6}")
print(f"  CG_sum(6) = {float(cg6):.12f}")

# CG_sum(6) = -11*sqrt(7)/5292 + 8*sqrt(42)/1323
# = (-11*sqrt(7) + 32*sqrt(42)) / 5292
# = (-11*sqrt(7) + 32*sqrt(42)) / 5292
# Note: sqrt(42) = sqrt(6*7) = sqrt(6)*sqrt(7)
# So = sqrt(7)*(-11 + 32*sqrt(6)) / 5292
cg6_alt = sqrt(7) * (-11 + 32*sqrt(6)) / 5292
print(f"  Factored: sqrt(7)*(-11 + 32*sqrt(6))/5292 = {float(cg6_alt):.12f}")
print(f"  Check: {float(simplify(cg6 - cg6_alt))}")

# V_mag(6) = (2/3)*(sqrt(9/7) - sqrt(6/8)) = (2/3)*(3/sqrt(7) - sqrt(3/4))
# = (2/3)*(3/sqrt(7) - sqrt(3)/2)
v_mag_6 = V_magnetic_closed(6)
print(f"  V_mag(6) = {v_mag_6} = {float(v_mag_6):.12f}")

# CG_sum(6) / V_mag(6)^2
ratio_6 = simplify(cg6 / v_mag_6**2)
print(f"  CG_sum(6)/V_mag(6)^2 = {ratio_6} = {float(ratio_6):.12f}")
ratio_6_r = radsimp(ratio_6)
print(f"  radsimp: {ratio_6_r}")

# CG_sum(6) * (7)*(8) = CG_sum(6) * 56
prod_6 = simplify(cg6 * 56)
print(f"  CG_sum(6)*56 = {prod_6} = {float(prod_6):.12f}")

# Now let's try: is CG_sum(6) / V_mag(6) rational?
ratio_6b = simplify(cg6 / v_mag_6)
ratio_6b_r = radsimp(ratio_6b)
print(f"  CG_sum(6)/V_mag(6) = {ratio_6b_r} = {float(ratio_6b_r):.12f}")


# 10. V(1,0)_6 = (2/3)*sqrt(9/7) = 2/(sqrt(7)) and V(0,1)_6 = -(2/3)*sqrt(6/8) = -(2/3)*sqrt(3/4)
# = -sqrt(3)/3
# CG_sum has sqrt(7) and sqrt(42) = sqrt(6)*sqrt(7)
# V(1,0) has sqrt(7) only, V(0,1) has sqrt(3) only
# So CG can be written as: a*sqrt(7) + b*sqrt(6)*sqrt(7) = sqrt(7)*(a + b*sqrt(6))
# And V(1,0) = 2*sqrt(7)/7 (proportional to sqrt(7))
# So CG / V(1,0) = (a + b*sqrt(6)) * 7 / (2*sqrt(7)/7) ...

# Let me just compute CG / V(1,0) and CG / V(0,1) directly
print(f"\n  CG(6) / V(1,0)_6 = {simplify(cg6 / V10(6))} = {float(simplify(cg6 / V10(6))):.12f}")
print(f"  CG(6) / V(0,1)_6 = {simplify(cg6 / V01(6))} = {float(simplify(cg6 / V01(6))):.12f}")

# CG(6) has sqrt(7) as a common factor. V(1,0)_6 also has sqrt(7).
# CG/V10 = sqrt(7)*(-11+32*sqrt(6))/5292 / (2*sqrt(7)/7)
# = (-11+32*sqrt(6)) * 7 / (5292 * 2)
# = (-11+32*sqrt(6)) / 1512
cg_over_v10 = (-11 + 32*sqrt(6)) / 1512
print(f"\n  Check: (-11+32*sqrt(6))/1512 = {float(cg_over_v10):.12f}")

# Can we express this as p + q*V01(6)/V10(6)?
# V01/V10 = -sqrt(n/(n+2)) / sqrt((n+3)/(n+1)) = -sqrt(n*(n+1)/((n+2)*(n+3)))
# At n=6: -sqrt(42/72) = -sqrt(7/12) = -sqrt(7)/(2*sqrt(3))
v_ratio_6 = simplify(V01(6) / V10(6))
print(f"  V01(6)/V10(6) = {v_ratio_6} = {float(v_ratio_6):.12f}")


print("\n\n=== SUMMARY ===")
print("CG_sum / V_mag^2 is almost exactly linear in n, with slope ~0.1582")
print("CG_sum * (n+1)*(n+2) is nearly linear: ~0.205 + 0.280*n + 0.028/n")
print("The CG_sum expressions have complex multi-sqrt structure")
print("Even the simplest (n=6: two sqrt terms) doesn't rationalize cleanly against V_mag")
print()
print("Next: need more data points (larger n) or a fundamentally different approach")
print("(e.g., Racah algebra / 6j symbols to evaluate the CG sums in closed form)")
