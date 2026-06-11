"""
Identify the exact closed form for V(1,0) per channel.

Discovered: V(1,0)_n / [2*sqrt(jL*(jL+1))] = r_n where:
  n=1: 1/3     n=2: 2/9     n=3: 1/6     n=4: 2/15
  n=5: 1/9     n=6: 2/21    n=7: 1/12    n=8: 2/27
  n=9: 1/15    n=10: 2/33

And V(0,1)_n = -[n/(n+2)] * V(1,0)_{n-1}  (exact recursion).

Goal: find the universal formula for r_n.
"""

from sympy import Rational, sqrt, simplify, S, factor

# The ratio values
ratios = {
    1: Rational(1, 3),
    2: Rational(2, 9),
    3: Rational(1, 6),
    4: Rational(2, 15),
    5: Rational(1, 9),
    6: Rational(2, 21),
    7: Rational(1, 12),
    8: Rational(2, 27),
    9: Rational(1, 15),
    10: Rational(2, 33),
}

print("=== Looking for universal formula ===")
for n, r in ratios.items():
    jL = Rational(n + 1, 2)
    jR = Rational(n, 2)
    dimL = 2*jL + 1  # = n+2
    dimR = 2*jR + 1  # = n+1

    # Candidate formulas
    # The denominators are: 3, 9, 6, 15, 9, 21, 12, 27, 15, 33
    # = 3, 9, 6, 15, 9, 21, 12, 27, 15, 33
    # Factor: 3, 3^2, 2*3, 3*5, 3^2, 3*7, 3*4, 3^3, 3*5, 3*11
    # All divisible by 3. d/3 = 1, 3, 2, 5, 3, 7, 4, 9, 5, 11
    # That's: 1, 3, 2, 5, 3, 7, 4, 9, 5, 11
    # Odd n: d/3 = 1, 2, 3, 4, 5 -> (n+1)/2
    # Even n: d/3 = 3, 5, 7, 9, 11 -> n+1

    # For odd n: r = 1/(3*(n+1)/2) = 2/(3*(n+1))
    # For even n: r = 2/(3*(n+1)) as well? n=2: 2/(3*3)=2/9 YES!
    # n=4: 2/(3*5)=2/15 YES!
    # n=6: 2/(3*7)=2/21 YES!
    # n=8: 2/(3*9)=2/27 YES!
    # n=10: 2/(3*11)=2/33 YES!

    # For odd n:
    # n=1: 1/3 vs 2/(3*2)=1/3 YES!
    # n=3: 1/6 vs 2/(3*4)=1/6 YES!
    # n=5: 1/9 vs 2/(3*6)=1/9 YES!
    # n=7: 1/12 vs 2/(3*8)=1/12 YES!
    # n=9: 1/15 vs 2/(3*10)=1/15 YES!

    # UNIVERSAL: r_n = 2/(3*(n+1)) for ALL n!

    candidate = Rational(2, 3*(n+1))
    match = (r == candidate)
    print(f"  n={n:2d}: r = {r}, 2/(3*(n+1)) = {candidate}, match = {match}")

print("\n=== RESULT: UNIVERSAL CLOSED FORM ===")
print("  V_channel(1,0)_n = 2*sqrt(jL*(jL+1)) * 2/(3*(n+1))")
print("                   = 4*sqrt(jL*(jL+1)) / (3*(n+1))")
print("  where jL = (n+1)/2, so jL*(jL+1) = (n+1)(n+3)/4")
print("  V(1,0)_n = 4*sqrt((n+1)(n+3)/4) / (3*(n+1))")
print("           = 2*sqrt((n+1)(n+3)) / (3*(n+1))")
print("           = 2*sqrt((n+3)/(n+1)) / 3")
print("           = 2/(3*sqrt((n+1)/(n+3)))")

print("\n  Simplified: V(1,0)_n = (2/3) * sqrt((n+3)/(n+1))")

# Verify
print("\n=== Verification ===")
v10_exact = {
    1: 2*sqrt(Rational(2))/3,
    2: 2*sqrt(Rational(15))/9,
    3: sqrt(Rational(6))/3,
    4: 2*sqrt(Rational(35))/15,
    5: 4*sqrt(Rational(3))/9,
    6: 2*sqrt(Rational(7))/7,
    7: sqrt(Rational(5))/3,
    8: 2*sqrt(Rational(11))/9,
    9: 2*sqrt(Rational(30))/15,
    10: 2*sqrt(Rational(143))/33,
}

for n, v_exact in v10_exact.items():
    v_formula = Rational(2, 3) * sqrt(Rational(n+3, n+1))
    diff = simplify(v_exact - v_formula)
    print(f"  n={n:2d}: exact = {float(v_exact):.10f}, formula = {float(v_formula):.10f}, "
          f"diff = {diff}")

print("\n\n=== CLOSED FORM FOR V(0,1) ===")
print("  V(0,1)_n = -[n/(n+2)] * V(1,0)_{n-1}")
print("           = -[n/(n+2)] * (2/3) * sqrt(n+2)/sqrt(n)")
print("           = -(2/3) * n/(n+2) * sqrt((n+2)/n)")
print("           = -(2/3) * sqrt(n/(n+2))")

v01_exact = {
    1: -2*sqrt(Rational(3))/9,
    2: -sqrt(Rational(2))/3,
    3: -2*sqrt(Rational(15))/15,
    4: -2*sqrt(Rational(6))/9,
    5: -2*sqrt(Rational(35))/21,
    6: -sqrt(Rational(3))/3,
    7: -2*sqrt(Rational(7))/9,
    8: -4*sqrt(Rational(5))/15,
    9: -2*sqrt(Rational(11))/11,
    10: -sqrt(Rational(30))/9,
}

print("\n  Verification of V(0,1) formula:")
for n, v_exact in v01_exact.items():
    # V(0,1)_n = -n/(n+2) * V(1,0)_{n-1}
    # V(1,0)_{n-1} = (2/3) * sqrt((n+2)/n)  [substituting n -> n-1 in the formula: (n-1+3)/(n-1+1) = (n+2)/n]
    # So V(0,1)_n = -n/(n+2) * (2/3) * sqrt((n+2)/n)
    #            = -(2/3) * n/(n+2) * sqrt((n+2)/n)
    #            = -(2/3) * sqrt(n*(n+2)) / (n+2)
    #            = -(2/3) * sqrt(n) / sqrt(n+2)
    #            = -(2/3) * sqrt(n/(n+2))
    v_formula = -Rational(2, 3) * sqrt(Rational(n, n+2))
    diff = simplify(v_exact - v_formula)
    print(f"  n={n:2d}: exact = {float(v_exact):.10f}, formula = {float(v_formula):.10f}, "
          f"diff = {diff}")

print("\n\n=== FULL V_MAGNETIC CLOSED FORM ===")
print("  V_magnetic(n) = V(1,0)_n + V(0,1)_n")
print("                = (2/3)*sqrt((n+3)/(n+1)) - (2/3)*sqrt(n/(n+2))")
print("                = (2/3) * [sqrt((n+3)/(n+1)) - sqrt(n/(n+2))]")

print("\n  Verification of total V_magnetic:")
for n in range(1, 11):
    jL = Rational(n + 1, 2)
    jR = Rational(n, 2)

    v_mag_formula = Rational(2, 3) * (sqrt(Rational(n+3, n+1)) - sqrt(Rational(n, n+2)))

    # Compare with exact computed values
    v_exact = v10_exact.get(n, S.Zero) + v01_exact.get(n, S.Zero)

    diff = simplify(v_mag_formula - v_exact)
    print(f"  n={n:2d}: formula = {float(v_mag_formula):.10f}, "
          f"exact = {float(v_exact):.10f}, diff = {diff}")

print("\n\n=== ASYMPTOTIC ANALYSIS ===")
# V_magnetic(n) = (2/3) * [sqrt((n+3)/(n+1)) - sqrt(n/(n+2))]
# For large n: sqrt((n+3)/(n+1)) ≈ 1 + 1/n - 1/(2n²) + ...
#              sqrt(n/(n+2)) ≈ 1 - 1/n + 1/(2n²) + ...
# Difference ≈ 2/n - 1/n² + ...
# So V_mag ≈ (2/3) * 2/n = 4/(3n) for large n

# More precisely: let f(x) = sqrt((1+3/x)/(1+1/x)) - sqrt(1/(1+2/x))
# At large x:
# (1+3/x)/(1+1/x) = 1 + 2/x - 2/x² + ...
# sqrt(...) = 1 + 1/x - 3/(2x²) + ...
# 1/(1+2/x) = 1 - 2/x + 4/x² - ...
# sqrt(...) = 1 - 1/x + 3/(2x²) - ...
# diff = 2/x - 3/x² + ...
# V_mag = (2/3)(2/x - 3/x² + ...) = 4/(3x) - 2/x² + ...
# with x = n

from sympy import series, Symbol, oo as symoo
x = Symbol('x', positive=True)
f1 = sqrt((x+3)/(x+1))
f2 = sqrt(x/(x+2))
v_asymp = Rational(2,3) * (f1 - f2)

# Series around x = infinity
v_series = series(v_asymp, x, symoo, 6)
print(f"  V_magnetic(n) = {v_series}")

# So V_mag * lambda where lambda = n + 3/2
# ≈ (4/(3n)) * n = 4/3 at leading order -- this explains the apparent 4/3 limit at moderate n
# But V_mag = 4/(3n) - 2/n² + ... and lambda = n + 3/2
# V_mag * lambda = (4/(3n) - 2/n² + ...) * (n + 3/2) = 4/3 + 2/n - 2/n + ...
# Need careful expansion

lam = x + Rational(3, 2)
v_lam = simplify(v_asymp * lam)
v_lam_series = series(v_lam, x, symoo, 6)
print(f"\n  V_magnetic * lambda = {v_lam_series}")

# V_mag * lambda^2
v_lam2 = simplify(v_asymp * lam**2)
v_lam2_series = series(v_lam2, x, symoo, 6)
print(f"\n  V_magnetic * lambda^2 = {v_lam2_series}")

print("\n\n=== NUMERICAL CHECK AT LARGE n ===")
for n in [10, 20, 50, 100, 200, 500, 1000]:
    v = float(Rational(2, 3) * (sqrt(Rational(n+3, n+1)) - sqrt(Rational(n, n+2))))
    lam = n + 1.5
    print(f"  n={n:5d}: V_mag = {v:.10e}, V*lam = {v*lam:.8f}, V*lam^2 = {v*lam**2:.6f}")

print("\n\n=== KEY INSIGHT FOR F2 ===")
print("  V_magnetic(n) = (2/3)[sqrt((n+3)/(n+1)) - sqrt(n/(n+2))]")
print("  Leading asymptotic: V_mag ~ 4/(3n) - 2/n^2 + O(1/n^3)")
print("  V_mag * lambda ~ 4/3 + constant/n + O(1/n^2)")
print("  V_mag * lambda^2 ~ 4n/3 + ... (diverges)")
print("  So V_mag ~ lambda^{-1} at leading order, NOT lambda^{-2}")
print("  The lambda^{-2} scaling is for F2 = B/V_mag, not for V_mag itself")

# Now let's see what V_mag^2 simplifies to
print("\n=== V_magnetic^2 simplified ===")
for n in range(1, 11):
    v = Rational(2, 3) * (sqrt(Rational(n+3, n+1)) - sqrt(Rational(n, n+2)))
    v2 = simplify(v**2)
    print(f"  n={n:2d}: V_mag^2 = {v2}")

# And V_mag^2 * something rational?
print("\n=== V_magnetic^2 * (n+1)(n+2) ===")
for n in range(1, 11):
    v = Rational(2, 3) * (sqrt(Rational(n+3, n+1)) - sqrt(Rational(n, n+2)))
    v2 = simplify(v**2 * (n+1) * (n+2))
    print(f"  n={n:2d}: V_mag^2*(n+1)(n+2) = {v2}")

print("\n=== V_magnetic^2 * 9/4 ===")
for n in range(1, 11):
    v = Rational(2, 3) * (sqrt(Rational(n+3, n+1)) - sqrt(Rational(n, n+2)))
    v2 = simplify(v**2 * Rational(9, 4))
    # V^2 * 9/4 = [(n+3)/(n+1) + n/(n+2) - 2*sqrt(n(n+3)/((n+1)(n+2)))]
    print(f"  n={n:2d}: V_mag^2*9/4 = {v2} = {float(v2):.10f}")

# Actually the full expansion:
# V_mag = (2/3)[sqrt((n+3)/(n+1)) - sqrt(n/(n+2))]
# V_mag^2 = (4/9)[(n+3)/(n+1) + n/(n+2) - 2*sqrt(n(n+3))/sqrt((n+1)(n+2))]
# The cross term involves sqrt(n(n+3)/(n+1)/(n+2)), which is irrational.
# So V_mag^2 is NOT rational in general. Let me verify:
print("\n=== Is V_mag^2 rational? ===")
for n in range(1, 11):
    v = Rational(2, 3) * (sqrt(Rational(n+3, n+1)) - sqrt(Rational(n, n+2)))
    v2 = simplify(v**2)
    # Check if it's rational
    v2_expanded = simplify(v2)
    print(f"  n={n:2d}: V_mag^2 = {v2_expanded}")
    # Try to check if this is rational by testing float vs round
    v2f = float(v2)
    # Check if any simple fraction matches
    from fractions import Fraction
    frac = Fraction(v2f).limit_denominator(10000)
    print(f"          ≈ {frac} ({float(frac):.10f}), "
          f"match = {abs(v2f - float(frac)) < 1e-10}")
