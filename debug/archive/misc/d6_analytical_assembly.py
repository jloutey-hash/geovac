"""
D6 analytical assembly: substitute all weight-11 Euler sum decompositions
into D6 = (-2289/512)z(10) + (567/128)z(11) + sum_r beta_r*(S_{r,11-r} - z(11))
and collect coefficients of each MZV basis element.

All Euler sum decompositions found via sub-basis PSLQ at 400 dps.
"""
import sys, io, functools
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
print = functools.partial(print, flush=True)

from sympy import Rational, S, simplify

# ============================================================
# Euler sum decompositions (weight 11)
# S_{r,s} = a*z(11) + b*z(2)z(9) + c*z(3)z(8) + d*z(4)z(7) + e*z(5)z(6)
# Format: {r: (a, b, c, d, e)}
# ============================================================
euler_sums = {
    1: (S(6), S(-1), S(-1), S(-1), S(-1)),          # Euler's formula
    2: (S(-27), S(9), S(2), S(6), S(4)),             # PSLQ 400 dps
    # 3: not needed (beta_3 = 0)
    4: (Rational(-329, 2), S(84), S(0), S(21), S(4)),  # PSLQ 400 dps
    5: (Rational(463, 2), S(-126), S(0), S(-21), S(0)), # PSLQ sub-basis
    6: (Rational(-461, 2), S(126), S(0), S(21), S(1)),  # PSLQ sub-basis
    7: (Rational(331, 2), S(-84), S(0), S(-20), S(-4)), # PSLQ sub-basis
    8: (S(-82), S(36), S(1), S(15), S(6)),             # stuffle from S_{3,8}
}

# Beta coefficients from D6 Euler sum expansion
betas = {
    1: Rational(315, 32),
    2: Rational(-245, 32),
    4: Rational(63, 32),
    5: Rational(35, 64),
    6: Rational(-21, 64),
    7: Rational(-21, 64),
    8: Rational(-7, 64),
}

# Rational part of D6
rat_z10 = Rational(-2289, 512)
rat_z11 = Rational(567, 128)

# ============================================================
# Assemble: D6 = rat_z10*z(10) + rat_z11*z(11) + sum_r beta_r*(S_r - z(11))
# S_r - z(11) = (a-1)*z(11) + b*z(2)z(9) + c*z(3)z(8) + d*z(4)z(7) + e*z(5)z(6)
# ============================================================

# Accumulate coefficients: [z(10), z(11), z(2)z(9), z(3)z(8), z(4)z(7), z(5)z(6)]
names = ["z(10)", "z(11)", "z(2)z(9)", "z(3)z(8)", "z(4)z(7)", "z(5)z(6)"]
coeffs = [rat_z10, rat_z11, S(0), S(0), S(0), S(0)]

print("=" * 70)
print("D6 ANALYTICAL ASSEMBLY")
print("=" * 70)

# Rational part
print(f"\nRational part:")
print(f"  z(10): {rat_z10}")
print(f"  z(11): {rat_z11}")

# Add each beta_r * (S_r - z(11))
for r, beta in sorted(betas.items()):
    a, b, c, d, e = euler_sums[r]
    # S_r - z(11) has z(11) coefficient (a - 1)
    z11_contrib = beta * (a - 1)
    z2z9_contrib = beta * b
    z3z8_contrib = beta * c
    z4z7_contrib = beta * d
    z5z6_contrib = beta * e

    coeffs[1] += z11_contrib
    coeffs[2] += z2z9_contrib
    coeffs[3] += z3z8_contrib
    coeffs[4] += z4z7_contrib
    coeffs[5] += z5z6_contrib

    print(f"\n  beta_{r} = {beta}, S_{{{r},{11-r}}} - z(11):")
    print(f"    z(11): {beta}*({a}-1) = {z11_contrib}")
    print(f"    z(2)z(9): {beta}*{b} = {z2z9_contrib}")
    print(f"    z(3)z(8): {beta}*{c} = {z3z8_contrib}")
    print(f"    z(4)z(7): {beta}*{d} = {z4z7_contrib}")
    print(f"    z(5)z(6): {beta}*{e} = {z5z6_contrib}")

# ============================================================
# Final result
# ============================================================
print("\n" + "=" * 70)
print("D6 DECOMPOSITION")
print("=" * 70)

for name, coeff in zip(names, coeffs):
    print(f"  {name}: {coeff}")

# Pretty print
terms = []
for name, coeff in zip(names, coeffs):
    if coeff != 0:
        terms.append(f"({coeff}){name}")
print(f"\n  D6 = {' + '.join(terms)}")

# ============================================================
# Product survival rule verification
# ============================================================
print("\n" + "=" * 70)
print("PRODUCT SURVIVAL RULE CHECK (p=6, weight 11)")
print("=" * 70)
print(f"  z(2)z(9) coeff = {coeffs[2]}")
if coeffs[2] == 0:
    print(f"  *** z(2)z(9) = 0 CONFIRMED ***")
    print(f"  Product survival rule: at p=6, weight=2p-1=11,")
    print(f"    z(2)*z(odd) cancels; surviving products = floor((2*6-5)/2) = 3")
    surv = sum(1 for c in coeffs[3:] if c != 0)
    print(f"    Actual surviving products: {surv}")
    print(f"    Expected: 3 (z(3)z(8), z(4)z(7), z(5)z(6))")
else:
    print(f"  *** z(2)z(9) != 0 — RULE VIOLATED ***")

# ============================================================
# Numerical verification at 200 dps
# ============================================================
print("\n" + "=" * 70)
print("NUMERICAL VERIFICATION")
print("=" * 70)

from mpmath import mp, mpf, nstr, zeta as mzeta
mp.dps = 200

D6_ref = mpf("-0.08420102765654912920293685898880703062501967433561437590179792")

z = lambda s: mzeta(s)

D6_assembled = (mpf(str(coeffs[0])) * z(10)
              + mpf(str(coeffs[1])) * z(11)
              + mpf(str(coeffs[2])) * z(2) * z(9)
              + mpf(str(coeffs[3])) * z(3) * z(8)
              + mpf(str(coeffs[4])) * z(4) * z(7)
              + mpf(str(coeffs[5])) * z(5) * z(6))

diff = abs(D6_assembled - D6_ref)
import mpmath
match = -int(mpmath.log10(diff)) if diff > 0 else 200
print(f"  D6 (assembled)  = {nstr(D6_assembled, 60)}")
print(f"  D6 (reference)  = {nstr(D6_ref, 60)}")
print(f"  Difference      = {nstr(diff, 5)}")
print(f"  Matching digits = {match}")

if match >= 50:
    print(f"  *** VERIFIED to {match} digits ***")
else:
    print(f"  *** FAILED — only {match} matching digits ***")

# ============================================================
# Compact summary for Paper 28
# ============================================================
print("\n" + "=" * 70)
print("PAPER 28 SUMMARY")
print("=" * 70)
print(f"  D6 = ({coeffs[0]})z(10) + ({coeffs[1]})z(11)")
for name, coeff in zip(names[2:], coeffs[2:]):
    if coeff != 0:
        print(f"       + ({coeff}){name}")
print(f"\n  Weight: 11 = 2*6 - 1")
print(f"  z(2)z(9): ABSENT (product survival rule)")
print(f"  Surviving products: {sum(1 for c in coeffs[3:] if c != 0)}")
print(f"  Expected by rule: max(0, floor((2*6-5)/2)) = 3")
print("Done.")
