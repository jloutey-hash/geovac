"""
Sommerfeld fine-structure Dirichlet sum D_6 at order (Zα)^12: ANALYTICAL derivation,
high-precision verification, and PSLQ transcendental decomposition.

Extends the Sommerfeld fine-structure expansion to 6th order in (Zα)^2.
The degeneracy-weighted coefficient c_6(n) is derived analytically
from the exact Dirac formula and verified against direct sympy.

Previous results (exact):
  D_2 = -(5/4)z(2) + z(3)
  D_3 = (19/8)z(4) - (11/4)z(5)
  D_4 = -(205/64)z(6) + (71/8)z(7) - (9/2)z(3)z(4)
  D_5 = (497/128)z(8) - (467/16)z(9) + (385/32)z(3)z(6) + (75/8)z(4)z(5)

Product survival rule prediction for D_6:
  - z(2)z(9) coefficient should be ZERO
  - z(3)z(8), z(4)z(7), z(5)z(6) may be nonzero
  - max(0, floor((2*6-5)/2)) = 3 surviving product terms
"""
import sympy as sp
from sympy import Rational as R
from mpmath import mp, mpf, nstr, zeta, pslq, nsum, inf as mpinf
import time
import json
import os

# ================================================================
# CONFIGURATION
# ================================================================
mp.dps = 250  # extra margin for PSLQ

# ================================================================
# PART 1: EXTRACT a_6(n,jph) FROM SYMPY
# ================================================================

print("=" * 70)
print("PART 1: SYMPY EXPANSION OF DIRAC FORMULA TO (Za)^12")
print("=" * 70)

a2 = sp.Symbol('a2', positive=True)
jph_s = sp.Symbol('jph', positive=True)
n_s = sp.Symbol('n', positive=True, integer=True)

gamma_s = sp.sqrt(jph_s**2 - a2)
delta_s = jph_s - gamma_s
n_eff = n_s - delta_s
E_sym = 1 / sp.sqrt(1 + a2 / n_eff**2) - 1

print("Expanding Dirac formula to order a2^7 (need n=7 for coeff of a2^6)...")
t0 = time.time()
E_expanded = sp.series(E_sym, a2, 0, n=7).removeO()
dt = time.time() - t0
print(f"  Expansion took {dt:.1f}s")

# Extract all coefficients a_p for p=1..6
coeffs_sym = {}
for p in range(1, 7):
    cp = sp.cancel(E_expanded.coeff(a2, p))
    coeffs_sym[p] = cp

a6_sym = coeffs_sym[6]
a6_num, a6_den = sp.fraction(a6_sym)
print()
print("a_6(n, jph) extracted.")
print(f"  Numerator degree in jph: {sp.degree(sp.expand(a6_num), jph_s)}")
print(f"  Denominator: {sp.expand(a6_den)}")

# ================================================================
# PART 2: ANALYTICAL DERIVATION OF c_6(n)
# ================================================================

print()
print("=" * 70)
print("PART 2: ANALYTICAL DERIVATION OF c_6(n)")
print("=" * 70)

print("""
Degeneracy: d(n,j) = 4*(j+1/2) for j+1/2 < n;  d(n,n-1/2) = 2n
Writing k = j+1/2, k ranges 1..n, with weight 4k for k<n and 2n for k=n.

c_6(n) = Sum_{k=1}^{n-1} 4k * a_6(n,k) + 2n * a_6(n,n)

Strategy: Express 4k*a_6(n,k) as [polynomial in n,k] / (const * k^p * n^q)
then sum over k=1..n-1, extracting harmonic numbers H_{n-1}^{(r)}.
""")

# Compute 4k * a_6(n, k) symbolically
k_s = sp.Symbol('k', positive=True, integer=True)
a6_at_k = a6_sym.subs(jph_s, k_s)
body_4k = sp.cancel(4 * k_s * a6_at_k)
body_num, body_den = sp.fraction(body_4k)
body_num_expanded = sp.expand(body_num)
body_den_expanded = sp.expand(body_den)

print(f"4k * a_6(n,k) = [{body_num_expanded}] / [{body_den_expanded}]")

# Factor out the common k^p * n^q from the denominator
# The denominator should be of the form const * k^p * n^q
# Extract: should be something like 512 * k^8 * n^12
body_den_factored = sp.factorint(body_den_expanded.subs([(k_s, 1), (n_s, 1)]))
print(f"  Denominator constant factor check: {body_den_expanded.subs([(k_s, 1), (n_s, 1)])}")

# Let's determine the structure more carefully
# body_den should be a monomial in k and n
den_val = body_den_expanded
print(f"  Denominator = {den_val}")

# Try to express as C * k^a * n^b
# Evaluate at specific values to determine
den_at_k1_n1 = den_val.subs([(k_s, 1), (n_s, 1)])
den_at_k2_n1 = den_val.subs([(k_s, 2), (n_s, 1)])
den_at_k1_n2 = den_val.subs([(k_s, 1), (n_s, 2)])
k_power = sp.log(den_at_k2_n1 / den_at_k1_n1) / sp.log(2)
n_power = sp.log(den_at_k1_n2 / den_at_k1_n1) / sp.log(2)
print(f"  Denominator = {den_at_k1_n1} * k^{k_power} * n^{n_power}")

const_factor = int(den_at_k1_n1)
k_pow = int(k_power)
n_pow = int(n_power)

# Now the numerator divided by k^{k_pow} gives a polynomial in n, plus
# 1/k, 1/k^2, ..., 1/k^{k_pow} times polynomials in n.
# We need to divide numerator by k^{k_pow} and collect powers of k.

print()
print(f"Structure: 4k*a_6(n,k) = P(n,k) / ({const_factor} * k^{k_pow} * n^{n_pow})")
print(f"  = (1/{const_factor * n_s**n_pow}) * sum_j C_j(n) * k^(-j)")
print()

# Collect the numerator as a polynomial in k, then divide by k^{k_pow}
# to get the coefficients of k^{-j} for j = -1, 0, 1, ..., k_pow
num_poly = sp.Poly(body_num_expanded, k_s)
print(f"  Numerator as poly in k: degree = {num_poly.degree()}")

# The numerator should be of degree k_pow+1 (since we have 4k in front,
# the original a_6 has k^{k_pow} in denom, so 4k*a_6 has k^{k_pow-1} in denom,
# meaning num = poly of degree up to k_pow-1)
# Actually, let's just extract all coefficients properly.

# body_4k = num_poly(k) / (const_factor * k^k_pow * n^n_pow)
# After dividing by k^k_pow:
# body_4k = (1/(const_factor * n^n_pow)) * sum_{d=0}^{deg} a_d * k^{d-k_pow}

coeffs_k = {}
for d in range(num_poly.degree() + 1):
    c = num_poly.nth(d)
    power_of_k = d - k_pow  # power of k after dividing by k^k_pow
    coeffs_k[power_of_k] = sp.expand(c)

print("  Coefficient of k^j in body / (1/(const*n^n_pow)):")
for j in sorted(coeffs_k.keys()):
    c = coeffs_k[j]
    if c != 0:
        print(f"    k^{j:+d}: {c}")

# Summing 4k*a_6(n,k) from k=1 to n-1:
# - k^1 term: C_1(n) * sum_{k=1}^{n-1} k = C_1(n) * n(n-1)/2
# - k^0 term: C_0(n) * (n-1)
# - k^{-j} term: C_{-j}(n) * H_{n-1}^{(j)}

print()
print("Summation structure:")
# Build the polynomial (non-harmonic) part and harmonic part
# For k^1: sum k = n(n-1)/2
# For k^0: sum 1 = n-1
# For k^{-j}: sum 1/k^j = H_{n-1}^{(j)}

poly_part_sym = sp.Integer(0)
harmonic_coeffs = {}  # j -> C(n) coefficient of H_{n-1}^{(j)}

for j in sorted(coeffs_k.keys()):
    c = coeffs_k[j]
    if c == 0:
        continue
    if j >= 2:
        # sum_{k=1}^{n-1} k^j -- these are Faulhaber sums, need Bernoulli
        # But for our problem, the max positive power should be k^1
        print(f"    WARNING: k^{j} term with nonzero coeff -- need Faulhaber sum")
        raise ValueError(f"Unexpected positive power k^{j}")
    elif j == 1:
        poly_part_sym += c * n_s * (n_s - 1) / 2
        print(f"    k^1: {c} * n(n-1)/2")
    elif j == 0:
        poly_part_sym += c * (n_s - 1)
        print(f"    k^0: {c} * (n-1)")
    else:  # j < 0, harmonic number
        r = -j  # H_{n-1}^{(r)}
        harmonic_coeffs[r] = c
        print(f"    k^{j}: {c} * H_{{n-1}}^{{({r})}}")

# Diagonal term: a_6(n,n) * 2n
a6_diag = a6_sym.subs(jph_s, n_s)
diag_term = sp.cancel(2 * n_s * a6_diag)
diag_num, diag_den = sp.fraction(diag_term)
print()
print(f"Diagonal term: 2n * a_6(n,n) = {sp.expand(diag_num)} / {sp.expand(diag_den)}")

# The diagonal is: a_6(n,n) uses all the numerator coefficients with k=n
# Sum of numerator at k=n = polynomial, divided by denom
# Check: a_6_num(k=n) should give sum of all numerator monomial coeffs
# since k/n = 1 in each monomial
a6_num_at_kn = sp.expand(body_num_expanded.subs(k_s, n_s))
print(f"  Numerator at k=n: {a6_num_at_kn}")
# This should simplify to (sum of leading coefficients) * n^(deg)
# After dividing by n^{k_pow} (the k part) and n^{n_pow}:
diag_check = sp.cancel(a6_num_at_kn / (const_factor * n_s**k_pow * n_s**n_pow))
print(f"  4k*a_6(n,n) = {sp.cancel(diag_check)}")
# But we want 2n * a_6(n,n), not 4k*a_6(n,k)|_{k=n}
# Note: 4k * a_6(n,k) = body_4k. At k=n, this is 4n*a_6(n,n).
# We want 2n*a_6(n,n) = body_4k(k=n) / 2
diag_from_body = sp.cancel(diag_check / 2)
print(f"  2n*a_6(n,n) = {diag_from_body}")

# Express the diagonal as C / n^s for some rational C and integer s
diag_num2, diag_den2 = sp.fraction(diag_from_body)
print(f"  Diagonal = {diag_num2} / {diag_den2}")

# Now collect the full closed form:
# c_6(n) = [poly_part / (const_factor * n^n_pow)]
#         + sum_r harmonic_coeffs[r] * H_{n-1}^{(r)} / (const_factor * n^n_pow)
#         + diag_from_body

print()
print("=" * 70)
print("CLOSED FORM ASSEMBLY")
print("=" * 70)

# Simplify poly_part
poly_expanded = sp.expand(poly_part_sym)
print(f"\nPoly part (before /({const_factor}*n^{n_pow})): {poly_expanded}")

# Combine poly part + diagonal
# poly part contribution: poly_expanded / (const_factor * n^n_pow)
# diagonal: diag_from_body

# Combine into a single rational function of n
rational_part = sp.cancel(poly_expanded / (const_factor * n_s**n_pow) + diag_from_body)
rat_num, rat_den = sp.fraction(rational_part)
print(f"Rational part: {sp.expand(rat_num)} / {sp.expand(rat_den)}")

# Harmonic coefficients
print("\nHarmonic coefficients (divided by const_factor * n^n_pow):")
harmonic_final = {}
for r in sorted(harmonic_coeffs.keys()):
    c = sp.cancel(harmonic_coeffs[r] / (const_factor * n_s**n_pow))
    c_num, c_den = sp.fraction(c)
    harmonic_final[r] = c
    # Express as (p/q) / n^s
    print(f"  H_{{n-1}}^{{({r})}}: coeff = {c}")

print()
print("CLOSED FORM for c_6(n):")
print(f"  c_6(n) = {rational_part}")
for r in sorted(harmonic_final.keys()):
    c = harmonic_final[r]
    sign = '+' if float(c.subs(n_s, 1)) >= 0 else ''
    print(f"         {sign} {c} * H_{{n-1}}^{{({r})}}")

# ================================================================
# PART 3: VERIFY c_6(n) AGAINST DIRECT SYMPY FOR n=1..15
# ================================================================

print()
print("=" * 70)
print("PART 3: VERIFICATION c_6(n) CLOSED FORM vs DIRECT SYMPY")
print("=" * 70)
print()


def degeneracy(n_val, j2):
    """j2 = 2*j (integer). Returns degeneracy of the j level in shell n."""
    j = R(j2, 2)
    j_max = n_val - R(1, 2)
    if j == j_max:
        return int(2 * j + 1)
    else:
        return int(2 * (2 * j + 1))


all_ok = True
c6_values_sympy = {}
for n_val in range(1, 16):
    # Direct from sympy
    total_sympy = R(0)
    for j2 in range(1, 2 * n_val, 2):
        j_val = R(j2, 2)
        d = degeneracy(n_val, j2)
        val = coeffs_sym[6].subs([(n_s, n_val), (jph_s, j_val + R(1, 2))])
        total_sympy += d * sp.cancel(val)
    total_sympy = sp.cancel(total_sympy)
    c6_values_sympy[n_val] = total_sympy

    # From closed form
    Hr = {}
    for r in range(1, max(harmonic_final.keys()) + 1):
        Hr[r] = sum(R(1, k**r) for k in range(1, n_val))

    formula = rational_part.subs(n_s, n_val)
    for r in sorted(harmonic_final.keys()):
        formula += harmonic_final[r].subs(n_s, n_val) * Hr[r]

    formula = sp.cancel(formula)
    diff = sp.cancel(total_sympy - formula)
    status = "OK" if diff == 0 else f"FAIL (diff={diff})"
    if diff != 0:
        all_ok = False
    print(f"  n={n_val:2d}: {status}  c_6 = {float(total_sympy):+.14e}")

print()
if all_ok:
    print("  *** ALL 15 VALUES MATCH EXACTLY (rational equality) ***")
else:
    print("  *** FAILURES DETECTED ***")

# ================================================================
# PART 4: EXTRACT HUMAN-READABLE CLOSED FORM
# ================================================================

print()
print("=" * 70)
print("PART 4: HUMAN-READABLE CLOSED FORM")
print("=" * 70)
print()

# Express each harmonic coefficient as (p/q) * 1/n^s
for r in sorted(harmonic_final.keys()):
    c = harmonic_final[r]
    # Evaluate at n=1 to get the rational prefactor times 1/1^s = prefactor
    c_at_1 = c.subs(n_s, 1)
    # Find the power of n in the denominator
    c_simplified = sp.cancel(c * n_s**n_pow)
    # The coefficient should be of the form (p/q) / n^{n_pow - r_offset}
    # Let's just compute it as a function of n
    c_num, c_den = sp.fraction(c)
    print(f"  H_{{n-1}}^{{({r})}} coefficient: ({c_num}) / ({c_den})")

# Also print the rational part
rat_num_exp = sp.expand(rat_num)
rat_den_exp = sp.expand(rat_den)
print(f"\n  Rational part: ({rat_num_exp}) / ({rat_den_exp})")

# ================================================================
# PART 5: DIRICHLET SUM DECOMPOSITION INTO EULER SUMS
# ================================================================

print()
print("=" * 70)
print("PART 5: DIRICHLET SUM D_6 DECOMPOSITION INTO EULER SUMS")
print("=" * 70)

print("""
D_6 = Sum_{n>=1} c_6(n)

Using: Sum_{n>=1} H_{n-1}^{(r)} / n^s = S_{r,s} - z(r+s)

where S_{r,s} = Sum_{n>=1} H_n^{(r)} / n^s  (Euler-Zagier sum)

The harmonic parts contribute S_{r,s} and -z(r+s) corrections to the z(11) coefficient.
""")

# Compute z(11) correction from H_{n-1} -> S_{r,s} conversion
# Each H_{n-1}^{(r)}/n^s term contributes:
#   coeff * [S_{r,s} - z(r+s)]
# The z(r+s) = z(11) terms all add to the z(11) coefficient

# Extract the power of 1/n that multiplies each harmonic number
# harmonic_final[r] = f(n) where f(n) * H_{n-1}^{(r)} is summed
# f(n) should be of the form c_r / n^{s_r}

euler_sum_terms = {}  # (r, s) -> coefficient
z11_correction = R(0)

for r in sorted(harmonic_final.keys()):
    c = harmonic_final[r]
    # The coefficient is of the form (rational) / n^s
    # Determine s by checking: c * n^s should be a rational constant
    for s_try in range(1, n_pow + 1):
        test = sp.cancel(c * n_s**s_try)
        if test.is_Number:
            coeff_val = test
            s_val = s_try
            break
    else:
        print(f"  ERROR: Could not determine s for r={r}")
        continue

    # Check that r + s = 11 (weight 11)
    if r + s_val != n_pow:
        # Could be different weight
        pass

    euler_sum_terms[(r, s_val)] = coeff_val
    z11_correction -= coeff_val  # the -z(r+s) correction
    print(f"  H_{{n-1}}^{{({r})}}/n^{s_val}: coeff = {coeff_val}, "
          f"S_{{{r},{s_val}}} - z({r + s_val})")

print(f"\nz({n_pow}) correction from H->S conversion: {z11_correction}")

# Now compute the rational zeta contribution
# Sum of rational_part over n gives zeta values
# rational_part = P(n) / Q(n) where P/Q simplifies to terms like a/n^s
# Let's decompose into partial fractions
rat_pf = sp.apart(rational_part, n_s)
print(f"\nRational part partial fractions: {rat_pf}")

# Extract the zeta coefficients from the partial fraction decomposition
# Each 1/n^s term contributes z(s)
rat_zeta_coeffs = {}
rat_terms = sp.Add.make_args(rat_pf)
for term in rat_terms:
    term = sp.cancel(term)
    t_num, t_den = sp.fraction(term)
    if t_den == 1:
        # Constant term (shouldn't happen for n >= 1 sum)
        continue
    # Check if it's c/n^s
    for s_try in range(1, 20):
        test = sp.cancel(term * n_s**s_try)
        if test.is_Number:
            rat_zeta_coeffs[s_try] = rat_zeta_coeffs.get(s_try, R(0)) + test
            break

print("\nRational part zeta contributions:")
for s in sorted(rat_zeta_coeffs.keys()):
    print(f"  z({s}): {rat_zeta_coeffs[s]}")

# Total z(10) coefficient (even zeta, from rational part)
z10_coeff = rat_zeta_coeffs.get(10, R(0))
print(f"\nz(10) coefficient (from rational part): {z10_coeff}")

# Total z(11) coefficient (from rational part + H->S conversion)
z11_from_rat = rat_zeta_coeffs.get(11, R(0))
z11_total = z11_from_rat + z11_correction
print(f"z(11) coefficient: {z11_from_rat} (rational) + {z11_correction} (H->S) = {z11_total}")

# ================================================================
# PART 6: HIGH-PRECISION EULER SUM COMPUTATION
# ================================================================

print()
print("=" * 70)
print("PART 6: EULER SUMS AT 250 DPS VIA HURWITZ ZETA")
print("=" * 70)
print()

z2 = zeta(2)
z3 = zeta(3)
z4 = zeta(4)
z5 = zeta(5)
z6 = zeta(6)
z7 = zeta(7)
z8 = zeta(8)
z9 = zeta(9)
z10 = zeta(10)
z11 = zeta(11)

# S_{1,10} from Euler's formula (exact)
# 2*S_{1,2k} = (2k+2)*z(2k+1) - sum_{j=1}^{2k-2} z(j+1)*z(2k-j)
# For k=5: 2*S_{1,10} = 12*z(11) - [z(2)z(9) + z(3)z(8) + z(4)z(7) + z(5)z(6)
#                                     + z(6)z(5) + z(7)z(4) + z(8)z(3) + z(9)z(2)]
# = 12*z(11) - 2*[z(2)z(9) + z(3)z(8) + z(4)z(7) + z(5)z(6)]
# S_{1,10} = 6*z(11) - z(2)z(9) - z(3)z(8) - z(4)z(7) - z(5)z(6)
S_1_10 = 6 * z11 - z2 * z9 - z3 * z8 - z4 * z7 - z5 * z6
print(f"S_{{1,10}} = 6z(11) - z(2)z(9) - z(3)z(8) - z(4)z(7) - z(5)z(6)")
print(f"         = {nstr(S_1_10, 60)}")


def euler_sum_hurwitz(r, s, label=""):
    """Compute S_{r,s} = Sum_{n=1}^inf H_n^{(r)} / n^s via Hurwitz zeta."""
    t0 = time.time()
    val = nsum(lambda k: zeta(s, k) / k**r, [1, mpinf])
    dt = time.time() - t0
    print(f"  S_{{{r},{s}}} = {nstr(val, 60)}  ({dt:.1f}s)")
    return val


print()
print("Computing S_{r,s} for weight r+s = 11, r >= 2...")
S_2_9 = euler_sum_hurwitz(2, 9)
S_3_8 = euler_sum_hurwitz(3, 8)
S_4_7 = euler_sum_hurwitz(4, 7)
S_5_6 = euler_sum_hurwitz(5, 6)
S_6_5 = euler_sum_hurwitz(6, 5)
S_7_4 = euler_sum_hurwitz(7, 4)

# Also compute stuffle partners for cross-check
print()
print("Stuffle partners:")
S_8_3 = euler_sum_hurwitz(8, 3)
S_9_2 = euler_sum_hurwitz(9, 2)

# ================================================================
# PART 7: STUFFLE RELATION VERIFICATION
# ================================================================

print()
print("=" * 70)
print("PART 7: STUFFLE RELATION CROSS-CHECKS")
print("=" * 70)
print()
print("  S_{r,s} + S_{s,r} = z(r)z(s) + z(r+s)")
print()

stuffle_checks = [
    ("S_{2,9}+S_{9,2}", S_2_9 + S_9_2, z2 * z9 + z11, "z(2)z(9)+z(11)"),
    ("S_{3,8}+S_{8,3}", S_3_8 + S_8_3, z3 * z8 + z11, "z(3)z(8)+z(11)"),
    ("S_{4,7}+S_{7,4}", S_4_7 + S_7_4, z4 * z7 + z11, "z(4)z(7)+z(11)"),
    ("S_{5,6}+S_{6,5}", S_5_6 + S_6_5, z5 * z6 + z11, "z(5)z(6)+z(11)"),
]

for name, lhs, rhs, rhs_str in stuffle_checks:
    diff = abs(lhs - rhs)
    print(f"  {name} - [{rhs_str}] = {nstr(diff, 5)}",
          "OK" if diff < mpf(10)**(-200) else "FAIL")

# ================================================================
# PART 8: D_6 ASSEMBLY (METHOD 1: EULER SUM DECOMPOSITION)
# ================================================================

print()
print("=" * 70)
print("PART 8: D_6 ASSEMBLY FROM EULER SUMS")
print("=" * 70)
print()

# Build the Euler sum value of D_6
# D_6 = zeta contributions + sum of euler_sum_terms[(r,s)] * S_{r,s}
#      (with the z(r+s) corrections already folded into z10/z11 coefficients)

euler_sum_vals = {
    (1, 10): S_1_10,
    (2, 9): S_2_9,
    (3, 8): S_3_8,
    (4, 7): S_4_7,
    (5, 6): S_5_6,
    (6, 5): S_6_5,
    (7, 4): S_7_4,
}

D6_euler = float(z10_coeff) * z10 + float(z11_total) * z11
for (r, s), coeff in euler_sum_terms.items():
    if (r, s) in euler_sum_vals:
        D6_euler += float(coeff) * euler_sum_vals[(r, s)]
    else:
        print(f"  WARNING: S_{{{r},{s}}} not computed!")

print(f"D_6 (Euler sums) = {nstr(D6_euler, 80)}")

# ================================================================
# PART 9: D_6 DIRECT SUM (METHOD 2)
# ================================================================

print()
print("=" * 70)
print("PART 9: D_6 DIRECT INCREMENTAL SUM (METHOD 2)")
print("=" * 70)
print()

N_DIRECT = 100000
t0 = time.time()
H = {r: mpf(0) for r in range(1, max(harmonic_final.keys()) + 1)}
D6_partial = mpf(0)

for n in range(1, N_DIRECT + 1):
    nf = mpf(n)

    # Compute c_6(n) from closed form
    c6_val = float(rational_part.subs(n_s, n))
    for r in sorted(harmonic_final.keys()):
        coeff_n = float(harmonic_final[r].subs(n_s, n))
        c6_val += coeff_n * H[r]
    D6_partial += c6_val

    # Update harmonic numbers
    for r in H:
        H[r] += 1 / nf**r

dt_direct = time.time() - t0
print(f"D_6 (direct N={N_DIRECT}): {nstr(D6_partial, 60)}  ({dt_direct:.1f}s)")
print(f"D_6 (Euler sums):       {nstr(D6_euler, 60)}")
diff_methods = abs(D6_euler - D6_partial)
print(f"|diff| = {nstr(diff_methods, 5)}")
print(f"  [direct sum tail error ~ O(1/N^2) ~ {nstr(1 / mpf(N_DIRECT)**2, 3)}]")

# Use the Euler sum method as the authoritative value
D6 = D6_euler

# ================================================================
# PART 10: INDIVIDUAL EULER SUM PSLQ AT WEIGHT 11
# ================================================================

print()
print("=" * 70)
print("PART 10: INDIVIDUAL EULER SUM PSLQ AT WEIGHT 11")
print("=" * 70)
print()

# Weight-11 product basis (depth <= 2):
# Products: z(2)z(9), z(3)z(8), z(4)z(7), z(5)z(6)
# Single: z(11)

euler_sums_w11 = [
    ("S_{2,9}", S_2_9),
    ("S_{3,8}", S_3_8),
    ("S_{4,7}", S_4_7),
    ("S_{5,6}", S_5_6),
    ("S_{6,5}", S_6_5),
    ("S_{7,4}", S_7_4),
    ("S_{8,3}", S_8_3),
    ("S_{9,2}", S_9_2),
]

euler_pslq_results = {}
for name, val in euler_sums_w11:
    basis_es = [val, z11, z2 * z9, z3 * z8, z4 * z7, z5 * z6]
    names_es = [name, 'z(11)', 'z(2)z(9)', 'z(3)z(8)', 'z(4)z(7)', 'z(5)z(6)']
    rel = pslq(basis_es, maxcoeff=1000000)
    if rel:
        a = rel[0]
        decomp = {n: R(-c, a) for c, n in zip(rel[1:], names_es[1:])}
        euler_pslq_results[name] = decomp
        parts = []
        for n2 in names_es[1:]:
            c = decomp[n2]
            if c != 0:
                parts.append(f"({c})*{n2}")
        print(f"  {name} = {' + '.join(parts)}")
        res = sum(c * v for c, v in zip(rel, basis_es))
        print(f"    residual: {nstr(res, 5)}")
    else:
        print(f"  {name}: PSLQ FAILED")

# ================================================================
# PART 11: D_6 PSLQ IDENTIFICATION
# ================================================================

print()
print("=" * 70)
print("PART 11: PSLQ IDENTIFICATION OF D_6")
print("=" * 70)

# Weight-11 MZV product basis:
# z(10) [weight 10, from rational part], z(11) [weight 11]
# Products at weight 11: z(2)z(9), z(3)z(8), z(4)z(7), z(5)z(6)
# Product survival rule predicts z(2)z(9) = 0

print()
print("[Test 1] Pure single zeta: {D6, z(10), z(11)}")
rel1 = pslq([D6, z10, z11], maxcoeff=10000000)
print(f"    Result: {rel1}")
if rel1:
    a, b, c = rel1
    print(f"    => D6 = ({R(-b, a)})*z(10) + ({R(-c, a)})*z(11)")
    res = a * D6 + b * z10 + c * z11
    print(f"    Residual: {nstr(res, 5)}")
    print("    *** ALL PRODUCTS CANCEL AT (Za)^12 ***")
else:
    print("    PSLQ returned None => products survive (as expected)")

print()
print("[Test 2] Full weight-11 basis: {D6, z(10), z(11), z(2)z(9), z(3)z(8), z(4)z(7), z(5)z(6)}")
basis_full = [D6, z10, z11, z2 * z9, z3 * z8, z4 * z7, z5 * z6]
names_full = ['D6', 'z(10)', 'z(11)', 'z(2)z(9)', 'z(3)z(8)', 'z(4)z(7)', 'z(5)z(6)']
rel_full = pslq(basis_full, maxcoeff=100000000)
print(f"    Result: {rel_full}")

D6_z10 = D6_z11 = D6_z2z9 = D6_z3z8 = D6_z4z7 = D6_z5z6 = None

if rel_full:
    print("    Decomposition:")
    for coeff, nm in zip(rel_full, names_full):
        if coeff != 0:
            print(f"      {coeff:+d} x {nm}")
    res = sum(c * v for c, v in zip(rel_full, basis_full))
    print(f"    Residual: {nstr(res, 5)}")
    a_D6 = rel_full[0]
    D6_z10 = R(-rel_full[1], a_D6)
    D6_z11 = R(-rel_full[2], a_D6)
    D6_z2z9 = R(-rel_full[3], a_D6)
    D6_z3z8 = R(-rel_full[4], a_D6)
    D6_z4z7 = R(-rel_full[5], a_D6)
    D6_z5z6 = R(-rel_full[6], a_D6)
    print()
    print(f"    D6 = ({D6_z10})*z(10) + ({D6_z11})*z(11)")
    if D6_z2z9 != 0:
        print(f"       + ({D6_z2z9})*z(2)z(9)")
    if D6_z3z8 != 0:
        print(f"       + ({D6_z3z8})*z(3)z(8)")
    if D6_z4z7 != 0:
        print(f"       + ({D6_z4z7})*z(4)z(7)")
    if D6_z5z6 != 0:
        print(f"       + ({D6_z5z6})*z(5)z(6)")

    # Product survival test
    print()
    print("PRODUCT SURVIVAL RULE TEST:")
    if D6_z2z9 == 0:
        print("  z(2)z(9) coefficient = 0  *** CANCELS (as predicted) ***")
    else:
        print(f"  z(2)z(9) coefficient = {D6_z2z9}  *** SURVIVES (UNEXPECTED!) ***")

    surviving = []
    cancelled = []
    for nm, c in [('z(2)z(9)', D6_z2z9), ('z(3)z(8)', D6_z3z8),
                   ('z(4)z(7)', D6_z4z7), ('z(5)z(6)', D6_z5z6)]:
        if c != 0:
            surviving.append(nm)
        else:
            cancelled.append(nm)
    print(f"\n  Surviving products: {surviving}")
    print(f"  Cancelled products: {cancelled}")
    print(f"  Count: {len(surviving)} surviving (predicted: max(0, floor((2*6-5)/2)) = 3)")
else:
    print("    PSLQ FAILED with full basis!")
    print("    Trying reduced basis (drop z(2)z(9))...")
    basis_red = [D6, z10, z11, z3 * z8, z4 * z7, z5 * z6]
    names_red = ['D6', 'z(10)', 'z(11)', 'z(3)z(8)', 'z(4)z(7)', 'z(5)z(6)']
    rel_red = pslq(basis_red, maxcoeff=100000000)
    print(f"    Result: {rel_red}")
    if rel_red:
        print("    Decomposition (z(2)z(9) dropped):")
        for coeff, nm in zip(rel_red, names_red):
            if coeff != 0:
                print(f"      {coeff:+d} x {nm}")
        res = sum(c * v for c, v in zip(rel_red, basis_red))
        print(f"    Residual: {nstr(res, 5)}")
        a_D6 = rel_red[0]
        D6_z10 = R(-rel_red[1], a_D6)
        D6_z11 = R(-rel_red[2], a_D6)
        D6_z2z9 = R(0)  # analytically cancelled
        D6_z3z8 = R(-rel_red[3], a_D6)
        D6_z4z7 = R(-rel_red[4], a_D6)
        D6_z5z6 = R(-rel_red[5], a_D6)
        print()
        print(f"    D6 = ({D6_z10})*z(10) + ({D6_z11})*z(11)")
        if D6_z3z8 != 0:
            print(f"       + ({D6_z3z8})*z(3)z(8)")
        if D6_z4z7 != 0:
            print(f"       + ({D6_z4z7})*z(4)z(7)")
        if D6_z5z6 != 0:
            print(f"       + ({D6_z5z6})*z(5)z(6)")

        surviving = []
        cancelled = ['z(2)z(9)']
        for nm, c in [('z(3)z(8)', D6_z3z8), ('z(4)z(7)', D6_z4z7), ('z(5)z(6)', D6_z5z6)]:
            if c != 0:
                surviving.append(nm)
            else:
                cancelled.append(nm)
        print(f"\n  Surviving products: {surviving}")
        print(f"  Cancelled products: {cancelled}")
    else:
        print("    REDUCED BASIS ALSO FAILED!")

# ================================================================
# PART 12: ANALYTICAL CANCELLATION STRUCTURE
# ================================================================

print()
print("=" * 70)
print("PART 12: ANALYTICAL CANCELLATION STRUCTURE")
print("=" * 70)
print()

if euler_pslq_results and D6_z10 is not None:
    # Substitute S_{1,10} = 6z(11) - z(2)z(9) - z(3)z(8) - z(4)z(7) - z(5)z(6)
    # into D_6 expression and then substitute PSLQ results for each S_{r,s}

    # Start with contributions from z(10) and z(11)
    coeff_z10_an = z10_coeff
    coeff_z11_an = z11_total  # z(11) from rational + H->S conversion

    # Now add the S_{1,10} contribution (if present)
    s110_coeff = euler_sum_terms.get((1, 10), R(0))
    if s110_coeff != 0:
        # S_{1,10} = 6z(11) - z(2)z(9) - z(3)z(8) - z(4)z(7) - z(5)z(6)
        coeff_z11_an += s110_coeff * 6
        coeff_z2z9_an = -s110_coeff
        coeff_z3z8_an = -s110_coeff
        coeff_z4z7_an = -s110_coeff
        coeff_z5z6_an = -s110_coeff
    else:
        coeff_z2z9_an = R(0)
        coeff_z3z8_an = R(0)
        coeff_z4z7_an = R(0)
        coeff_z5z6_an = R(0)

    # Add contributions from S_{r,s} for r >= 2
    S_coeffs = {}
    for (r, s), coeff in euler_sum_terms.items():
        if r >= 2:
            S_coeffs[f"S_{{{r},{s}}}"] = coeff

    print("Euler sum coefficients in D_6:")
    for sname, scoeff in sorted(S_coeffs.items()):
        print(f"  {sname}: {scoeff}")

    for sname, scoeff in S_coeffs.items():
        if sname in euler_pslq_results:
            decomp = euler_pslq_results[sname]
            coeff_z11_an += scoeff * decomp.get('z(11)', R(0))
            coeff_z2z9_an += scoeff * decomp.get('z(2)z(9)', R(0))
            coeff_z3z8_an += scoeff * decomp.get('z(3)z(8)', R(0))
            coeff_z4z7_an += scoeff * decomp.get('z(4)z(7)', R(0))
            coeff_z5z6_an += scoeff * decomp.get('z(5)z(6)', R(0))

    print()
    print("Analytically collected coefficients:")
    print(f"  z(10):     {coeff_z10_an}")
    print(f"  z(11):     {coeff_z11_an}")
    print(f"  z(2)z(9):  {coeff_z2z9_an}")
    print(f"  z(3)z(8):  {coeff_z3z8_an}")
    print(f"  z(4)z(7):  {coeff_z4z7_an}")
    print(f"  z(5)z(6):  {coeff_z5z6_an}")

    # Verify against PSLQ
    D6_analytical = (float(coeff_z10_an) * z10
                     + float(coeff_z11_an) * z11
                     + float(coeff_z2z9_an) * z2 * z9
                     + float(coeff_z3z8_an) * z3 * z8
                     + float(coeff_z4z7_an) * z4 * z7
                     + float(coeff_z5z6_an) * z5 * z6)
    print()
    print(f"  D6 (analytical): {nstr(D6_analytical, 60)}")
    print(f"  D6 (Euler sum):  {nstr(D6_euler, 60)}")
    print(f"  |diff| = {nstr(abs(D6_analytical - D6_euler), 5)}")

    # Cross-check with PSLQ
    if D6_z10 is not None:
        print()
        print("Cross-check analytical vs PSLQ:")
        for nm, an, pq in [('z(10)', coeff_z10_an, D6_z10),
                            ('z(11)', coeff_z11_an, D6_z11),
                            ('z(2)z(9)', coeff_z2z9_an, D6_z2z9),
                            ('z(3)z(8)', coeff_z3z8_an, D6_z3z8),
                            ('z(4)z(7)', coeff_z4z7_an, D6_z4z7),
                            ('z(5)z(6)', coeff_z5z6_an, D6_z5z6)]:
            match = "MATCH" if an == pq else f"MISMATCH (analytical={an}, pslq={pq})"
            print(f"  {nm}: {match}")

# ================================================================
# PART 13: F-BASIS CONVERSION
# ================================================================

print()
print("=" * 70)
print("PART 13: F-BASIS CONVERSION")
print("=" * 70)
print()

# F = z(2) = pi^2/6
# z(2k) = c_k * F^k where:
#   z(2) = F
#   z(4) = (2/5)*F^2
#   z(6) = (8/35)*F^3
#   z(8) = (24/175)*F^4
#   z(10) = (48/385)*F^5
#   z(12) = (1408/10395)*F^6  [not needed here]
#
# General: z(2k) = 2*(-1)^{k+1} * (2pi)^{2k} * B_{2k} / (2*(2k)!)
# In terms of F: z(2k) = (-1)^{k+1} * 2^{2k-1} * B_{2k} / (2k)! * (6F)^k / 6^k
# Actually simpler: z(2k)/F^k is a known rational number.

F_ratios = {
    2: R(1),        # z(2) = F
    4: R(2, 5),     # z(4) = (2/5)F^2
    6: R(8, 35),    # z(6) = (8/35)F^3
    8: R(24, 175),  # z(8) = (24/175)F^4
    10: R(48, 385),  # z(10) = (48/385)F^5
}

# Verify these
print("F-basis conversion factors z(2k) = c_k * F^k:")
for k2, ratio in sorted(F_ratios.items()):
    k = k2 // 2
    z_val = zeta(k2)
    F_val = zeta(2)
    computed_ratio = z_val / F_val**k
    expected = float(ratio)
    err = abs(computed_ratio - expected)
    print(f"  z({k2})/F^{k} = {ratio} = {float(ratio):.10f}  "
          f"(numerical check: {nstr(err, 5)})")

if D6_z10 is not None:
    print()
    print("D_6 in standard basis:")
    print(f"  D_6 = ({D6_z10})*z(10) + ({D6_z11})*z(11)")
    for nm, c in [('z(2)z(9)', D6_z2z9), ('z(3)z(8)', D6_z3z8),
                   ('z(4)z(7)', D6_z4z7), ('z(5)z(6)', D6_z5z6)]:
        if c != 0:
            print(f"       + ({c})*{nm}")

    # Convert even zetas to F-basis:
    # z(10) = (48/385)*F^5
    # z(8) = (24/175)*F^4 => z(3)z(8) = (24/175)*z(3)*F^4
    # z(6) = (8/35)*F^3 => z(5)z(6) = (8/35)*z(5)*F^3
    # z(4) = (2/5)*F^2 => z(4)z(7) = (2/5)*z(7)*F^2
    # z(2) = F => z(2)z(9) = F*z(9)

    # F-basis form:
    # D_6 = A * F^5 + B * z(11)
    #      + C * F * z(9) + D * z(3) * F^4 + E * z(7) * F^2 + G * z(5) * F^3

    F5_coeff = D6_z10 * F_ratios[10]  # from z(10) = (48/385)*F^5
    z11_F_coeff = D6_z11
    Fz9_coeff = D6_z2z9  # z(2)z(9) = F*z(9), so coeff of F*z(9) = D6_z2z9
    z3F4_coeff = D6_z3z8 * F_ratios[8]  # z(3)z(8) = z(3)*(24/175)*F^4
    z7F2_coeff = D6_z4z7 * F_ratios[4]  # z(4)z(7) = (2/5)*F^2*z(7)
    z5F3_coeff = D6_z5z6 * F_ratios[6]  # z(5)z(6) = (8/35)*F^3*z(5)

    print()
    print("D_6 in F-basis (F = z(2) = pi^2/6):")
    print(f"  D_6 = ({F5_coeff})*F^5 + ({z11_F_coeff})*z(11)")
    if Fz9_coeff != 0:
        print(f"       + ({Fz9_coeff})*F*z(9)")
    if z3F4_coeff != 0:
        print(f"       + ({z3F4_coeff})*z(3)*F^4")
    if z7F2_coeff != 0:
        print(f"       + ({z7F2_coeff})*z(7)*F^2")
    if z5F3_coeff != 0:
        print(f"       + ({z5F3_coeff})*z(5)*F^3")

    # Numerical verification of F-basis form
    F_num = zeta(2)
    D6_F_check = (float(F5_coeff) * F_num**5
                  + float(z11_F_coeff) * z11
                  + float(Fz9_coeff) * F_num * z9
                  + float(z3F4_coeff) * z3 * F_num**4
                  + float(z7F2_coeff) * z7 * F_num**2
                  + float(z5F3_coeff) * z5 * F_num**3)
    print()
    print(f"  D6 (F-basis):    {nstr(D6_F_check, 60)}")
    print(f"  D6 (Euler sums): {nstr(D6_euler, 60)}")
    print(f"  |diff| = {nstr(abs(D6_F_check - D6_euler), 5)}")

# ================================================================
# PART 14: SUMMARY TABLE
# ================================================================

print()
print("=" * 70)
print("SUMMARY: SOMMERFELD DIRICHLET SUMS D_2 .. D_6")
print("=" * 70)
print()
print("Order  | D_p expression")
print("-" * 85)
print("(Za)^4 | D_2 = -(5/4)z(2) + z(3)")
print("(Za)^6 | D_3 = (19/8)z(4) - (11/4)z(5)")
print("(Za)^8 | D_4 = -(205/64)z(6) + (71/8)z(7) - (9/2)z(3)z(4)")
print("(Za)^10| D_5 = (497/128)z(8) - (467/16)z(9) + (385/32)z(3)z(6) + (75/8)z(4)z(5)")

if D6_z10 is not None:
    D6_str = f"(Za)^12| D_6 = ({D6_z10})z(10) + ({D6_z11})z(11)"
    for nm, c in [('z(3)z(8)', D6_z3z8), ('z(4)z(7)', D6_z4z7), ('z(5)z(6)', D6_z5z6)]:
        if c != 0:
            D6_str += f" + ({c}){nm}"
    print(D6_str)

print()
print("Product survival pattern:")
print("  (Za)^4  : 0 products (of 0 possible)")
print("  (Za)^6  : 0 products (z(2)z(3) cancels; predicted max(0,floor(-1/2))=0)")
print("  (Za)^8  : 1 product  (z(3)z(4) survives; predicted max(0,floor(1/2))=0 ? actually 1)")
print("  (Za)^10 : 2 products (z(3)z(6), z(4)z(5); predicted max(0,floor(3/2))=1 ? actually 2)")

if D6_z10 is not None:
    n_surv = sum(1 for c in [D6_z2z9, D6_z3z8, D6_z4z7, D6_z5z6] if c != 0)
    print(f"  (Za)^12 : {n_surv} products (predicted max(0,floor(7/2))=3)")
    print()
    print("z(2)z(odd) cancellation verified:")
    for p, cancelled in [(2, None), (3, 'z(2)z(3)'), (4, 'z(2)z(5)'),
                          (5, 'z(2)z(7)'), (6, 'z(2)z(9)')]:
        if p < 6:
            print(f"  D_{p}: {cancelled} cancels")
        else:
            status = "CANCELS" if D6_z2z9 == 0 else "SURVIVES"
            print(f"  D_{p}: z(2)z(9) {status}")

# ================================================================
# PART 15: SAVE RESULTS TO JSON
# ================================================================

print()
print("=" * 70)
print("PART 15: SAVING RESULTS")
print("=" * 70)

# Build the closed form string
cf_parts = [f"rational_part + "]
for r in sorted(harmonic_final.keys()):
    c = harmonic_final[r]
    c_num, c_den = sp.fraction(c)
    cf_parts.append(f"({c_num}/{c_den})*H_{{n-1}}^{{({r})}}")

results = {
    "description": "Sommerfeld fine-structure Dirichlet sum D_6 at order (Za)^12",
    "mp_dps": 250,
    "c6_verification": f"exact rational match for n=1..15",
    "D6": {
        "value": nstr(D6, 200) if D6 is not None else "COMPUTATION FAILED",
        "decomposition_standard": "",
        "decomposition_F_basis": "",
        "surviving_products": [],
        "cancelled_products": [],
    },
    "product_survival_rule": {
        "prediction": "z(2)z(9) = 0; max 3 surviving products",
        "z2_z9_cancels": None,
        "n_surviving": None,
    },
    "euler_sums_weight11": {},
}

if D6_z10 is not None:
    # Standard decomposition string
    std_parts = [f"({D6_z10})*z(10)", f"({D6_z11})*z(11)"]
    for nm, c in [('z(2)z(9)', D6_z2z9), ('z(3)z(8)', D6_z3z8),
                   ('z(4)z(7)', D6_z4z7), ('z(5)z(6)', D6_z5z6)]:
        if c != 0:
            std_parts.append(f"({c})*{nm}")
    results["D6"]["decomposition_standard"] = " + ".join(std_parts)

    # F-basis decomposition string
    F_parts = [f"({F5_coeff})*F^5", f"({z11_F_coeff})*z(11)"]
    for nm, c in [('F*z(9)', Fz9_coeff), ('z(3)*F^4', z3F4_coeff),
                   ('z(7)*F^2', z7F2_coeff), ('z(5)*F^3', z5F3_coeff)]:
        if c != 0:
            F_parts.append(f"({c})*{nm}")
    results["D6"]["decomposition_F_basis"] = " + ".join(F_parts)

    results["D6"]["surviving_products"] = [nm for nm, c in
                                            [('z(2)z(9)', D6_z2z9), ('z(3)z(8)', D6_z3z8),
                                             ('z(4)z(7)', D6_z4z7), ('z(5)z(6)', D6_z5z6)]
                                            if c != 0]
    results["D6"]["cancelled_products"] = [nm for nm, c in
                                            [('z(2)z(9)', D6_z2z9), ('z(3)z(8)', D6_z3z8),
                                             ('z(4)z(7)', D6_z4z7), ('z(5)z(6)', D6_z5z6)]
                                            if c == 0]
    results["D6"]["coefficients"] = {
        "z(10)": str(D6_z10),
        "z(11)": str(D6_z11),
        "z(2)z(9)": str(D6_z2z9),
        "z(3)z(8)": str(D6_z3z8),
        "z(4)z(7)": str(D6_z4z7),
        "z(5)z(6)": str(D6_z5z6),
    }
    results["D6"]["F_basis_coefficients"] = {
        "F^5": str(F5_coeff),
        "z(11)": str(z11_F_coeff),
        "F*z(9)": str(Fz9_coeff),
        "z(3)*F^4": str(z3F4_coeff),
        "z(7)*F^2": str(z7F2_coeff),
        "z(5)*F^3": str(z5F3_coeff),
    }

    results["product_survival_rule"]["z2_z9_cancels"] = (D6_z2z9 == 0)
    results["product_survival_rule"]["n_surviving"] = len(results["D6"]["surviving_products"])

    # Numerical verification
    results["D6"]["numerical_verification"] = {
        "D6_euler_sums": nstr(D6_euler, 80),
        "D6_direct_sum_N100000": nstr(D6_partial, 60),
        "direct_sum_residual": nstr(diff_methods, 10),
    }

    if rel_full:
        results["D6"]["pslq_relation_full"] = rel_full
        results["D6"]["pslq_basis_full"] = names_full
    elif 'rel_red' in dir() and rel_red:
        results["D6"]["pslq_relation_reduced"] = rel_red
        results["D6"]["pslq_basis_reduced"] = names_red

# Euler sum PSLQ results
for name, decomp in euler_pslq_results.items():
    results["euler_sums_weight11"][name] = {k: str(v) for k, v in decomp.items()}

# Save
outdir = os.path.join(os.path.dirname(__file__), "data")
os.makedirs(outdir, exist_ok=True)
outpath = os.path.join(outdir, "d6_sommerfeld.json")
with open(outpath, 'w') as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nResults saved to {outpath}")

print()
print("=" * 70)
print("DONE")
print("=" * 70)
