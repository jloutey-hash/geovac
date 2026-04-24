"""
Find the exact closed form for V_channel(1,0) and V_channel(0,1).

Known exact values:
n=1:  V(1,0) = 2*sqrt(2)/3        V(0,1) = -2*sqrt(3)/9
n=2:  V(1,0) = 2*sqrt(15)/9       V(0,1) = -sqrt(2)/3
n=3:  V(1,0) = sqrt(6)/3          V(0,1) = -2*sqrt(15)/15
n=4:  V(1,0) = 2*sqrt(35)/15      V(0,1) = -2*sqrt(6)/9
n=5:  V(1,0) = 4*sqrt(3)/9        V(0,1) = -2*sqrt(35)/21
n=6:  V(1,0) = 2*sqrt(7)/7        V(0,1) = -sqrt(3)/3
n=7:  V(1,0) = sqrt(5)/3          V(0,1) = -2*sqrt(7)/9
n=8:  V(1,0) = 2*sqrt(11)/9       V(0,1) = -4*sqrt(5)/15
n=9:  V(1,0) = 2*sqrt(30)/15      V(0,1) = -2*sqrt(11)/11
n=10: V(1,0) = 2*sqrt(143)/33     V(0,1) = -sqrt(30)/9
"""

from sympy import Rational, sqrt, simplify, S, factor, nsimplify

# Channel (1,0) values
# jL = (n+1)/2 for the electron, q=1 probe with (jpL, jpR) = (1, 0)
# The CG structure should give a formula in terms of jL, jR

# Let's square these and look at the rational numbers
print("=== V(1,0)^2 ===")
v10_data = {
    1: (2, 2, 3),    # 2*sqrt(2)/3 -> 8/9
    2: (2, 15, 9),   # 2*sqrt(15)/9 -> 60/81 = 20/27
    3: (1, 6, 3),    # sqrt(6)/3 -> 6/9 = 2/3
    4: (2, 35, 15),  # 2*sqrt(35)/15 -> 140/225 = 28/45
    5: (4, 3, 9),    # 4*sqrt(3)/9 -> 48/81 = 16/27
    6: (2, 7, 7),    # 2*sqrt(7)/7 -> 28/49 = 4/7
    7: (1, 5, 3),    # sqrt(5)/3 -> 5/9
    8: (2, 11, 9),   # 2*sqrt(11)/9 -> 44/81
    9: (2, 30, 15),  # 2*sqrt(30)/15 -> 120/225 = 8/15
    10: (2, 143, 33),# 2*sqrt(143)/33 -> 572/1089 = 52/99
}

# Compute V(1,0)^2 = coeff^2 * arg_under_sqrt / denom^2
for n, (c, arg, d) in v10_data.items():
    v_sq = Rational(c**2 * arg, d**2)
    jL = Rational(n + 1, 2)
    jR = Rational(n, 2)
    dimL = 2*jL + 1  # = n+2
    dimR = 2*jR + 1  # = n+1

    # Try: V^2 = 4*jL*(jL+1) / (2*jL+1)^2
    formula1 = 4 * jL * (jL + 1) / dimL**2
    # = (n+1)(n+3) / (n+2)^2
    formula2 = Rational((n+1)*(n+3), (n+2)**2)

    print(f"  n={n:2d}: V^2 = {v_sq} = {float(v_sq):.8f}, "
          f"(n+1)(n+3)/(n+2)^2 = {formula2} = {float(formula2):.8f}, "
          f"match = {simplify(v_sq - formula2) == 0}")

print("\n=== V(0,1)^2 ===")
v01_data = {
    1: (2, 3, 9),      # 2*sqrt(3)/9 -> 12/81 = 4/27
    2: (1, 2, 3),       # sqrt(2)/3 -> 2/9
    3: (2, 15, 15),     # 2*sqrt(15)/15 -> 60/225 = 4/15
    4: (2, 6, 9),       # 2*sqrt(6)/9 -> 24/81 = 8/27
    5: (2, 35, 21),     # 2*sqrt(35)/21 -> 140/441 = 20/63
    6: (1, 3, 3),       # sqrt(3)/3 -> 3/9 = 1/3
    7: (2, 7, 9),       # 2*sqrt(7)/9 -> 28/81
    8: (4, 5, 15),      # 4*sqrt(5)/15 -> 80/225 = 16/45
    9: (2, 11, 11),     # 2*sqrt(11)/11 -> 44/121 = 4/11
    10: (1, 30, 9),     # sqrt(30)/9 -> 30/81 = 10/27
}

for n, (c, arg, d) in v01_data.items():
    v_sq = Rational(c**2 * arg, d**2)
    jL = Rational(n + 1, 2)
    jR = Rational(n, 2)
    dimL = 2*jL + 1  # = n+2
    dimR = 2*jR + 1  # = n+1

    # Try: V^2 = 4*jR*(jR+1) / (2*jR+1)^2
    formula1 = 4 * jR * (jR + 1) / dimR**2
    # = n*(n+2) / (n+1)^2
    formula2 = Rational(n*(n+2), (n+1)**2)

    print(f"  n={n:2d}: V^2 = {v_sq} = {float(v_sq):.8f}, "
          f"n(n+2)/(n+1)^2 = {formula2} = {float(formula2):.8f}, "
          f"match = {simplify(v_sq - formula2) == 0}")

# These do NOT match. Let's see what the actual pattern is.
print("\n=== Finding the pattern for V(1,0)^2 ===")
for n, (c, arg, d) in v10_data.items():
    v_sq = Rational(c**2 * arg, d**2)
    # Factor numerator and denominator
    print(f"  n={n:2d}: V(1,0)^2 = {v_sq}, "
          f"num = {v_sq.p} = {int(v_sq.p)}, "
          f"den = {v_sq.q} = {int(v_sq.q)}")

print("\n=== Finding the pattern for V(0,1)^2 ===")
for n, (c, arg, d) in v01_data.items():
    v_sq = Rational(c**2 * arg, d**2)
    print(f"  n={n:2d}: V(0,1)^2 = {v_sq}, "
          f"num = {v_sq.p} = {int(v_sq.p)}, "
          f"den = {v_sq.q} = {int(v_sq.q)}")

# Look for the formula by checking simple rationals
print("\n=== Trying formulas for V(1,0)^2 ===")
for n, (c, arg, d) in v10_data.items():
    v_sq = Rational(c**2 * arg, d**2)
    # Check 4*(n+1)/(3*(n+2)) ?
    for a in range(1, 20):
        for b in range(1, 20):
            candidate = Rational(a, b)
            # Try various forms
            for p, q in [(n+1, n+2), (n+2, n+3), (n, n+1), (n+1, n+3),
                          ((n+1)*(n+2), (n+3)*(n+2)), (n*(n+1), (n+2)*(n+3)),
                          ((n+1)**2, (n+2)**2)]:
                if q == 0:
                    continue
                trial = candidate * Rational(p, q)
                if trial == v_sq:
                    print(f"  n={n}: V^2 = {candidate} * {p}/{q}")

print("\n=== V(1,0)^2 * (n+2)^2 ===")
for n, (c, arg, d) in v10_data.items():
    v_sq = Rational(c**2 * arg, d**2)
    val = v_sq * (n+2)**2
    print(f"  n={n:2d}: V(1,0)^2 * (n+2)^2 = {val}")

print("\n=== V(0,1)^2 * (n+1)^2 ===")
for n, (c, arg, d) in v01_data.items():
    v_sq = Rational(c**2 * arg, d**2)
    val = v_sq * (n+1)**2
    print(f"  n={n:2d}: V(0,1)^2 * (n+1)^2 = {val}")

# V(1,0)^2 * (n+2)^2 should give us the numerator pattern
# V(0,1)^2 * (n+1)^2 should give us the numerator pattern

print("\n=== Pattern search: V(1,0)^2 * denominator ===")
for n, (c, arg, d) in v10_data.items():
    v_sq = Rational(c**2 * arg, d**2)
    # jL = (n+1)/2
    # 2jL+1 = n+2, (2jL+1)^2 = (n+2)^2
    # V^2 * (n+2)^2 should be an integer or simple fraction
    val = v_sq * Rational((n+2)**2, 1)
    # This is c^2 * arg / d^2 * (n+2)^2
    # But d is already (n+2) or 3 or ... let's just check
    print(f"  n={n:2d}: c={c}, arg={arg}, d={d}, "
          f"c^2*arg={c**2*arg}, d^2={d**2}, "
          f"(n+2)^2={(n+2)**2}, "
          f"V^2*(n+2)^2 = {val}")

# The denominators d are: 3, 9, 3, 15, 9, 7, 3, 9, 15, 33
# These are: (n+2): 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
# d values:    3, 9, 3, 15, 9, 7, 3, 9, 15, 33
# Hmm. d is NOT simply (n+2). Let me look at d more carefully.
# d = {3, 9, 3, 15, 9, 7, 3, 9, 15, 33}
# = {3, 9, 3, 15, 9, 7, 3, 9, 15, 33}
# Factored: 3, 3^2, 3, 3*5, 3^2, 7, 3, 3^2, 3*5, 3*11

# Actually, looking at V(1,0) directly:
# n=1: 2*sqrt(2)/3     = 2*sqrt(2*1)/3
# n=2: 2*sqrt(15)/9    = 2*sqrt(3*5)/(3^2)
# n=3: sqrt(6)/3       = sqrt(2*3)/3
# n=4: 2*sqrt(35)/15   = 2*sqrt(5*7)/(3*5)
# n=5: 4*sqrt(3)/9     = 4*sqrt(3)/(3^2)
# n=6: 2*sqrt(7)/7     = 2*sqrt(7)/7
# n=7: sqrt(5)/3       = sqrt(5)/3
# n=8: 2*sqrt(11)/9    = 2*sqrt(11)/(3^2)
# n=9: 2*sqrt(30)/15   = 2*sqrt(2*3*5)/(3*5)
# n=10: 2*sqrt(143)/33 = 2*sqrt(11*13)/(3*11)

# The arguments under sqrt:
# n=1: 2 = 1*2
# n=2: 15 = 3*5
# n=3: 6 = 2*3
# n=4: 35 = 5*7
# n=5: 3 = 1*3
# n=6: 7 = 1*7
# n=7: 5 = 1*5
# n=8: 11 = 1*11
# n=9: 30 = 2*3*5
# n=10: 143 = 11*13

# Let me check: does (n+1)/2 * (n+3)/2 match the arg?
# jL = (n+1)/2, jL+1 = (n+3)/2
# jL*(jL+1) = (n+1)(n+3)/4
# For n=1: (2*4)/4 = 2 YES
# For n=2: (3*5)/4 = 15/4, but arg=15 -> sqrt(15/4) * 2/(something)?
#   V = 2*sqrt(15)/9. If V = 2*sqrt(jL*(jL+1))/(something), then
#   2*sqrt(15/4)/(something) = 2*sqrt(15)/9 -> sqrt(15/4)*2/9 = sqrt(15)/2 * 2/9 = sqrt(15)/9
#   But we need 2*sqrt(15)/9, so 2*sqrt(15/4)*2/(something)?
#   Hmm, let me just check V^2 directly.

# V(1,0)^2:
# n=1: 8/9
# n=2: 60/81 = 20/27
# n=3: 6/9 = 2/3
# n=4: 140/225 = 28/45
# n=5: 48/81 = 16/27
# n=6: 28/49 = 4/7
# n=7: 5/9
# n=8: 44/81
# n=9: 120/225 = 8/15
# n=10: 572/1089 = 52/99

# Let's check: V^2 = 4*jL*(jL+1) / (2*jL+1)^2
# = (n+1)(n+3) / (n+2)^2
# n=1: 2*4/9 = 8/9 YES!
# n=2: 3*5/16 = 15/16 != 20/27 NO

# So the formula works for n=1 but not n=2.
# The CG structure is more complex than a single Racah coefficient.

# Let me look at the V(0,1) values and check if there's a shifted version
# V(0,1):
# n=1: -2*sqrt(3)/9      V(0,1)^2 = 12/81 = 4/27
# n=2: -sqrt(2)/3         V(0,1)^2 = 2/9
# Note: V(0,1) at n=2 = -sqrt(2)/3 = V(1,0) at n=1 * (-1) * sqrt(3/2)/3
# Actually V(0,1) at n=2 has sqrt(2), and V(1,0) at n=1 has sqrt(2).
# V(0,1) at n=3 has sqrt(15), and V(1,0) at n=2 has sqrt(15).
# V(0,1) at n=4 has sqrt(6), and V(1,0) at n=3 has sqrt(6).
# V(0,1) at n=5 has sqrt(35), and V(1,0) at n=4 has sqrt(35).

print("\n\n=== KEY OBSERVATION: V(0,1) at n uses the SAME sqrt as V(1,0) at n-1 ===")
for n in range(2, 11):
    v10_prev = v10_data.get(n-1)
    v01_curr = v01_data.get(n)
    if v10_prev and v01_curr:
        print(f"  V(1,0) at n={n-1}: {v10_prev[0]}*sqrt({v10_prev[1]})/{v10_prev[2]}")
        print(f"  V(0,1) at n={n}:   {v01_curr[0]}*sqrt({v01_curr[1]})/{v01_curr[2]}")
        print(f"    sqrt args match: {v10_prev[1] == v01_curr[1]}")
        # Ratio of V(0,1,n)/V(1,0,n-1)
        v10_f = v10_prev[0] * v10_prev[1]**0.5 / v10_prev[2]
        v01_f = v01_curr[0] * v01_curr[1]**0.5 / v01_curr[2]
        print(f"    V(0,1,n)/V(1,0,n-1) = {-v01_f/v10_f:.8f}")
        print()

# So V(0,1) at n = -r(n) * V(1,0) at n-1
# where r(n) is some coefficient.
# V(1,0)_1 = 2*sqrt(2)/3,  V(0,1)_2 = -sqrt(2)/3,    ratio = 1/2
# V(1,0)_2 = 2*sqrt(15)/9, V(0,1)_3 = -2*sqrt(15)/15, ratio = 2*9/(2*15) = 9/15 = 3/5
# V(1,0)_3 = sqrt(6)/3,    V(0,1)_4 = -2*sqrt(6)/9,   ratio = 2*3/(1*9) = 2/3
# V(1,0)_4 = 2*sqrt(35)/15,V(0,1)_5 = -2*sqrt(35)/21, ratio = 2*15/(2*21) = 15/21 = 5/7
# V(1,0)_5 = 4*sqrt(3)/9,  V(0,1)_6 = -sqrt(3)/3,     ratio = 1*9/(4*3) = 3/4
# V(1,0)_6 = 2*sqrt(7)/7,  V(0,1)_7 = -2*sqrt(7)/9,   ratio = 2*7/(2*9) = 7/9
# V(1,0)_7 = sqrt(5)/3,    V(0,1)_8 = -4*sqrt(5)/15,  ratio = 4*3/(1*15) = 4/5
# V(1,0)_8 = 2*sqrt(11)/9, V(0,1)_9 = -2*sqrt(11)/11, ratio = 2*9/(2*11) = 9/11
# V(1,0)_9 = 2*sqrt(30)/15,V(0,1)_10 = -sqrt(30)/9,   ratio = 1*15/(2*9) = 5/6

# Ratios: 1/2, 3/5, 2/3, 5/7, 3/4, 7/9, 4/5, 9/11, 5/6
# = n/(2n-1): 1/2, 2/3 no...
# n=2: 1/2 = 1/2
# n=3: 3/5
# n=4: 2/3 = 4/6
# n=5: 5/7
# n=6: 3/4 = 6/8
# n=7: 7/9
# n=8: 4/5 = 8/10
# n=9: 9/11
# n=10: 5/6 = 10/12

# Pattern: n/(n+1) for even n, n/(n+2) for odd n?
# n=2: 1/2 = (n/2)/(n/2+1) = 1/2 YES if pattern is floor((n+1)/2) / (floor((n+1)/2)+1)
# Actually: the pattern is (n-1)/(n+1) for odd n, n/(n+1) for even n? Let me check:
# Ratios are: n/(2n-1) for ... no.
# Let me just try: r(n) = (n-1)/(n+1)? n=2: 1/3 NO
# r(n) = n/(n+2)? n=2: 2/4=1/2 YES, n=3: 3/5 YES, n=4: 4/6=2/3 YES!!
print("\n=== Ratio pattern ===")
print("  V(0,1) at n = -[n/(n+2)] * V(1,0) at n-1")
for n in range(2, 11):
    v10_prev = v10_data.get(n-1)
    v01_curr = v01_data.get(n)
    if v10_prev and v01_curr:
        v10_f = v10_prev[0] * v10_prev[1]**0.5 / v10_prev[2]
        v01_f = v01_curr[0] * v01_curr[1]**0.5 / v01_curr[2]
        actual_ratio = -v01_f / v10_f
        predicted_ratio = n / (n + 2)
        print(f"  n={n}: actual ratio = {actual_ratio:.8f}, "
              f"n/(n+2) = {predicted_ratio:.8f}, "
              f"match = {abs(actual_ratio - predicted_ratio) < 1e-10}")

# PERFECT MATCH!
# V(0,1)_n = -[n/(n+2)] * V(1,0)_{n-1}
# This means V_magnetic = V(1,0)_n - [n/(n+2)] * V(1,0)_{n-1}
# This is a recursion. Now we need V(1,0)_n formula.

# V(1,0)^2:
# n=1: 8/9
# n=2: 20/27
# n=3: 2/3
# n=4: 28/45
# n=5: 16/27
# n=6: 4/7
# n=7: 5/9
# n=8: 44/81
# n=9: 8/15
# n=10: 52/99

# Can simplify: write V(1,0)^2 as f(n) and look for pattern
# Let me compute V(1,0)^2 * 9 for odd n and see:
# n=1: 8/9 * 9 = 8
# n=3: 2/3 * 9 = 6
# n=5: 16/27 * 9 = 16/3
# n=7: 5/9 * 9 = 5
# n=9: 8/15 * 9 = 24/5

# Not clean. Let me try V(1,0)^2 * (2jL+1)^2 = V(1,0)^2 * (n+2)^2:
print("\n=== V(1,0)^2 * (n+2)^2 ===")
for n, (c, arg, d) in v10_data.items():
    v_sq_scaled = Rational(c**2 * arg * (n+2)**2, d**2)
    print(f"  n={n:2d}: V(1,0)^2 * (n+2)^2 = {v_sq_scaled}")

# Results:
# n=1: 8, n=2: 80/3, n=3: 50/3, n=4: 168/5, n=5: 784/27, n=6: 256/7,
# n=7: 405/9=45, n=8: 4400/81, n=9: 968/15, n=10: 7488/99

# Not clean. Let me try a different approach.
# Look at the sqrt arguments under V(1,0):
# n=1: sqrt(2), n=2: sqrt(15), n=3: sqrt(6), n=4: sqrt(35),
# n=5: sqrt(3), n=6: sqrt(7), n=7: sqrt(5), n=8: sqrt(11),
# n=9: sqrt(30), n=10: sqrt(143)

# Product of neighbors: 2*3=6=jL(jL+1) at n=1 -> (1*2)(1*3)/4 = 6/4... no
# Actually: jR = n/2, so jR*(jR+1) = n(n+2)/4 -> n(n+2) = 3, 8, 15, 24, 35, 48, 63, 80, 99, 120
# jL = (n+1)/2, jL*(jL+1) = (n+1)(n+3)/4 -> (n+1)(n+3) = 8, 15, 24, 35, 48, 63, 80, 99, 120, 143

# sqrt args for V(1,0): 2, 15, 6, 35, 3, 7, 5, 11, 30, 143
# (n+1)(n+3)/4: 2, 15/4, 6, 35/4, 12, 63/4, 20, 99/4, 30, 143/4

# AH! For odd n: (n+1)(n+3)/4 IS the sqrt arg!
# n=1: 2*4/4 = 2 YES
# n=3: 4*6/4 = 6 YES
# n=5: 6*8/4 = 12 but arg is 3. NO.

# For n=5: arg=3, (n+1)(n+3)/4 = 6*8/4=12. sqrt(3) != sqrt(12).
# But sqrt(12) = 2*sqrt(3). So V(1,0)_5 = 4*sqrt(3)/9 = 2*sqrt(12)/9. Hmm.

# Actually: 4*sqrt(3)/9 = 2*2*sqrt(3)/9. And 2*sqrt(12)/9 = 2*2*sqrt(3)/9 = 4*sqrt(3)/9. YES!
# So V(1,0)_5 = 2*sqrt(12)/9 where 12 = (n+1)(n+3)/4 * 4/(something).
# 12 = 6*8/4 = 12. So sqrt arg = (n+1)(n+3)/4 = 12, and V = 2*sqrt(12)/9 = 4*sqrt(3)/9.

# Let me rewrite all V(1,0) with the arg = jL*(jL+1) = (n+1)(n+3)/4:
print("\n\n=== Rewriting V(1,0) using jL*(jL+1) ===")
for n, (c, arg, d) in v10_data.items():
    jL_sq = Rational((n+1)*(n+3), 4)
    v_exact = c * sqrt(arg) / d
    v_from_jL = float(2 * sqrt(float(jL_sq)) / d)  # trying 2*sqrt(jL*(jL+1))/d
    print(f"  n={n:2d}: jL*(jL+1)={jL_sq}, sqrt(jL*(jL+1))={float(sqrt(jL_sq)):.6f}, "
          f"exact V(1,0)={float(v_exact):.8f}, "
          f"2*sqrt(jL*(jL+1))/d = {v_from_jL:.8f}")

# That doesn't match because d varies. Let me check V(1,0) = 2*sqrt(jL*(jL+1)) / D(n)
# where D(n) is the denominator.
print("\n=== V(1,0) / [2*sqrt(jL*(jL+1))] ===")
for n, (c, arg, d) in v10_data.items():
    jL_sq = Rational((n+1)*(n+3), 4)
    v_exact_sym = Rational(c,1) * sqrt(Rational(arg,1)) / Rational(d,1)
    ratio = simplify(v_exact_sym / (2 * sqrt(jL_sq)))
    print(f"  n={n:2d}: V(1,0) / [2*sqrt(jL*(jL+1))] = {ratio} = {float(ratio):.8f}")
