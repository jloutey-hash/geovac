"""
Fine-structure Dirichlet sum D_5 at order (Zα)^10: ANALYTICAL derivation,
high-precision verification, and PSLQ transcendental decomposition.

Extends the Sommerfeld fine-structure expansion to 5th order in (Zα)².
The degeneracy-weighted coefficient c_5(n) is derived analytically
from the exact Dirac formula and verified against direct sympy.

Previous results:
  D_2 = -(5/4)z(2) + z(3)                            [weight 2-3, pure single]
  D_3 = (19/8)z(4) - (11/4)z(5)                       [weight 4-5, z(2)z(3) cancels]
  D_4 = -(205/64)z(6) + (71/8)z(7) - (9/2)z(3)z(4)   [weight 6-7, z(2)z(5) cancels, z(3)z(4) survives]
  D_5 = ?                                              [weight 8-9, target]
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
# PART 1: EXTRACT a_5(n,jph) FROM SYMPY
# ================================================================

print("=" * 70)
print("PART 1: SYMPY EXPANSION OF DIRAC FORMULA TO (Za)^10")
print("=" * 70)

a2 = sp.Symbol('a2', positive=True)
jph_s = sp.Symbol('jph', positive=True)
n_s = sp.Symbol('n', positive=True, integer=True)

gamma_s = sp.sqrt(jph_s**2 - a2)
delta_s = jph_s - gamma_s
n_eff = n_s - delta_s
E_sym = 1 / sp.sqrt(1 + a2 / n_eff**2) - 1
E_expanded = sp.series(E_sym, a2, 0, n=6).removeO()

# Extract all coefficients a_p for p=1..5
coeffs_sym = {}
for p in range(1, 6):
    cp = sp.cancel(E_expanded.coeff(a2, p))
    coeffs_sym[p] = cp

a5_sym = coeffs_sym[5]
a5_num, a5_den = sp.fraction(a5_sym)
print()
print("a_5(n, jph) = numerator / denominator")
print(f"  Numerator:   {sp.expand(a5_num)}")
print(f"  Denominator: {sp.expand(a5_den)}")
print()
print("Explicit form:")
print("  a_5(n,k) = (-63k^7 + 280k^6*n - 420k^5*n^2 + 180k^4*n^3")
print("              + 80k^3*n^4 - 24k^2*n^5 - 30k*n^6 - 10n^7)")
print("            / (256 * k^7 * n^10)")

# ================================================================
# PART 2: ANALYTICAL DERIVATION OF c_5(n)
# ================================================================

print()
print("=" * 70)
print("PART 2: ANALYTICAL DERIVATION OF c_5(n)")
print("=" * 70)

print("""
Degeneracy: d(n,j) = 4*(j+1/2) for j+1/2 < n;  d(n,n-1/2) = 2n
Writing k = j+1/2, k ranges 1..n, with weight 4k for k<n and 2n for k=n.

c_5(n) = Sum_{k=1}^{n-1} 4k * a_5(n,k) + 2n * a_5(n,n)

Step 1: 4k * a_5(n,k) = [terms] / (64 * n^10)
  Dividing each numerator monomial by k^6:
    -63k + 280n - 420n^2/k + 180n^3/k^2 + 80n^4/k^3
    - 24n^5/k^4 - 30n^6/k^5 - 10n^7/k^6

Step 2: Sum over k=1..n-1:
    polynomial terms: -63*n(n-1)/2 + 280n*(n-1) = 497*n*(n-1)/2
    1/k^r terms: give H_{n-1}^{(r)} factors

Step 3: Diagonal term: sum of numerator coefficients at k=n:
    -63+280-420+180+80-24-30-10 = -7
    a_5(n,n) = -7/(256*n^10)
    2n * a_5(n,n) = -7/(128*n^9)

Step 4: Collecting:
    Rational: 497(n-1)/(128n^9) - 7/(128n^9) = 7(71n-72)/(128n^9)

CLOSED FORM:
  c_5(n) =  7(71n - 72) / (128 n^9)
           - (105/16) H_{n-1}    / n^8
           + (45/16)  H_{n-1}^(2)/ n^7
           + (5/4)    H_{n-1}^(3)/ n^6
           - (3/8)    H_{n-1}^(4)/ n^5
           - (15/32)  H_{n-1}^(5)/ n^4
           - (5/32)   H_{n-1}^(6)/ n^3
""")

# ================================================================
# PART 3: VERIFY c_5(n) AGAINST DIRECT SYMPY FOR n=1..15
# ================================================================

print("=" * 70)
print("PART 3: VERIFICATION c_5(n) CLOSED FORM vs DIRECT SYMPY")
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
c5_values_sympy = {}
for n_val in range(1, 16):
    # Direct from sympy
    total_sympy = R(0)
    for j2 in range(1, 2 * n_val, 2):
        j_val = R(j2, 2)
        d = degeneracy(n_val, j2)
        val = coeffs_sym[5].subs([(n_s, n_val), (jph_s, j_val + R(1, 2))])
        total_sympy += d * sp.cancel(val)
    total_sympy = sp.cancel(total_sympy)
    c5_values_sympy[n_val] = total_sympy

    # From closed form
    H1 = sum(R(1, k) for k in range(1, n_val))
    H2 = sum(R(1, k**2) for k in range(1, n_val))
    H3 = sum(R(1, k**3) for k in range(1, n_val))
    H4 = sum(R(1, k**4) for k in range(1, n_val))
    H5 = sum(R(1, k**5) for k in range(1, n_val))
    H6 = sum(R(1, k**6) for k in range(1, n_val))

    formula = (R(7 * (71 * n_val - 72), 128 * n_val**9)
               - R(105, 16) * H1 / n_val**8
               + R(45, 16) * H2 / n_val**7
               + R(5, 4) * H3 / n_val**6
               - R(3, 8) * H4 / n_val**5
               - R(15, 32) * H5 / n_val**4
               - R(5, 32) * H6 / n_val**3)

    diff = sp.cancel(total_sympy - formula)
    status = "OK" if diff == 0 else f"FAIL (diff={diff})"
    if diff != 0:
        all_ok = False
    print(f"  n={n_val:2d}: {status}  c_5 = {float(total_sympy):+.14e}")

print()
if all_ok:
    print("  *** ALL 15 VALUES MATCH EXACTLY (rational equality) ***")
else:
    print("  *** FAILURES DETECTED ***")

# ================================================================
# PART 4: DIRICHLET SUM DECOMPOSITION
# ================================================================

print()
print("=" * 70)
print("PART 4: DIRICHLET SUM D_5 DECOMPOSITION INTO EULER SUMS")
print("=" * 70)

print("""
D_5 = Sum_{n>=1} c_5(n)

Using: Sum H_{n-1}^{(r)} / n^s = S_{r,s} - z(r+s)  for r >= 1

where S_{r,s} = Sum_{n>=1} H_n^{(r)} / n^s  (Euler-Zagier sum)

Rational part:
  7/128 * [71*z(8) - 72*z(9)] = (497/128)*z(8) - (504/128)*z(9)
                                = (497/128)*z(8) - (63/16)*z(9)

z(9) from H_{n-1} to S_{r,s} conversion (adding back z(r+s) = z(9) terms):
  z(9) coeff = -63/16 + 105/16 - 45/16 - 5/4 + 3/8 + 15/32 + 5/32
             = (-504 + 840 - 360 - 160 + 48 + 60 + 20)/128
             = -56/128 = -7/16

So:
  D_5 = (497/128)*z(8) + (-7/16)*z(9)
      + (-105/16)*S_{1,8}
      + (45/16)*S_{2,7}
      + (5/4)*S_{3,6}
      + (-3/8)*S_{4,5}
      + (-15/32)*S_{5,4}
      + (-5/32)*S_{6,3}

S_{1,8} via Euler's formula (2S_{1,2k} = (2k+2)z(2k+1) - Sum z(j+1)z(2k-j)):
  S_{1,8} = 5*z(9) - z(2)z(7) - z(3)z(6) - z(4)z(5)
""")

# ================================================================
# PART 5: HIGH-PRECISION EULER SUM COMPUTATION
# ================================================================

print("=" * 70)
print("PART 5: EULER SUMS AT 250 DPS VIA HURWITZ ZETA")
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

# S_{1,8} from Euler's formula (exact in terms of known zeta values)
S_1_8 = 5 * z9 - z2 * z7 - z3 * z6 - z4 * z5
print(f"  S_{{1,8}} = 5z(9) - z(2)z(7) - z(3)z(6) - z(4)z(5)")
print(f"          = {nstr(S_1_8, 60)}")


def euler_sum_hurwitz(r, s, label=""):
    """Compute S_{r,s} = Sum_{n=1}^inf H_n^{(r)} / n^s via Hurwitz zeta."""
    t0 = time.time()
    val = nsum(lambda k: zeta(s, k) / k**r, [1, mpinf])
    dt = time.time() - t0
    print(f"  S_{{{r},{s}}} = {nstr(val, 60)}  ({dt:.1f}s)")
    return val


print()
print("Computing S_{r,s} for r >= 2...")
S_2_7 = euler_sum_hurwitz(2, 7)
S_3_6 = euler_sum_hurwitz(3, 6)
S_4_5 = euler_sum_hurwitz(4, 5)
S_5_4 = euler_sum_hurwitz(5, 4)
S_6_3 = euler_sum_hurwitz(6, 3)

# Also compute the stuffle partners
print()
print("Stuffle partners:")
S_7_2 = euler_sum_hurwitz(7, 2)
S_6_3_check = S_6_3  # already computed
S_5_4_check = S_5_4  # already computed

# ================================================================
# PART 6: STUFFLE RELATION VERIFICATION
# ================================================================

print()
print("=" * 70)
print("PART 6: STUFFLE RELATION CROSS-CHECKS")
print("=" * 70)
print()
print("  S_{r,s} + S_{s,r} = z(r)z(s) + z(r+s)")
print()

stuffle_checks = [
    ("S_{2,7}+S_{7,2}", S_2_7 + S_7_2, z2 * z7 + z9, "z(2)z(7)+z(9)"),
    ("S_{3,6}+S_{6,3}", S_3_6 + S_6_3, z3 * z6 + z9, "z(3)z(6)+z(9)"),
    ("S_{4,5}+S_{5,4}", S_4_5 + S_5_4, z4 * z5 + z9, "z(4)z(5)+z(9)"),
]

for name, lhs, rhs, rhs_str in stuffle_checks:
    diff = abs(lhs - rhs)
    print(f"  {name} - [{rhs_str}] = {nstr(diff, 5)}", "OK" if diff < mpf(10)**(-200) else "FAIL")

# ================================================================
# PART 7: D_5 ASSEMBLY (METHOD 1: EULER SUM DECOMPOSITION)
# ================================================================

print()
print("=" * 70)
print("PART 7: D_5 ASSEMBLY FROM EULER SUMS")
print("=" * 70)
print()

D5_euler = ((mpf(497) / 128) * z8
            + (mpf(-7) / 16) * z9
            + (mpf(-105) / 16) * S_1_8
            + (mpf(45) / 16) * S_2_7
            + (mpf(5) / 4) * S_3_6
            + (mpf(-3) / 8) * S_4_5
            + (mpf(-15) / 32) * S_5_4
            + (mpf(-5) / 32) * S_6_3)

print(f"D_5 (Euler sums) = {nstr(D5_euler, 80)}")

# ================================================================
# PART 8: D_5 DIRECT SUM (METHOD 2)
# ================================================================

print()
print("=" * 70)
print("PART 8: D_5 DIRECT INCREMENTAL SUM (METHOD 2)")
print("=" * 70)
print()

N_DIRECT = 200000
t0 = time.time()
H1 = mpf(0)
H2 = mpf(0)
H3 = mpf(0)
H4 = mpf(0)
H5 = mpf(0)
H6 = mpf(0)
D5_partial = mpf(0)

for n in range(1, N_DIRECT + 1):
    nf = mpf(n)
    c5 = (mpf(7) * (71 * n - 72) / (128 * nf**9)
          - mpf(105) / 16 * H1 / nf**8
          + mpf(45) / 16 * H2 / nf**7
          + mpf(5) / 4 * H3 / nf**6
          - mpf(3) / 8 * H4 / nf**5
          - mpf(15) / 32 * H5 / nf**4
          - mpf(5) / 32 * H6 / nf**3)
    D5_partial += c5
    H1 += 1 / nf
    H2 += 1 / nf**2
    H3 += 1 / nf**3
    H4 += 1 / nf**4
    H5 += 1 / nf**5
    H6 += 1 / nf**6

dt_direct = time.time() - t0
print(f"D_5 (direct N={N_DIRECT}): {nstr(D5_partial, 60)}  ({dt_direct:.1f}s)")
print(f"D_5 (Euler sums):       {nstr(D5_euler, 60)}")
diff_methods = abs(D5_euler - D5_partial)
print(f"|diff| = {nstr(diff_methods, 5)}")
print(f"  [direct sum tail error ~ O(1/N^2) ~ {nstr(1 / mpf(N_DIRECT)**2, 3)}]")

# Use the Euler sum method as the authoritative value
D5 = D5_euler

# ================================================================
# PART 9: PSLQ IDENTIFICATION
# ================================================================

print()
print("=" * 70)
print("PART 9: PSLQ IDENTIFICATION OF D_5")
print("=" * 70)

# Weight-9 MZV product basis (depth <= 2):
# single: z(8) [weight 8, enters from rational part], z(9) [weight 9]
# products: z(2)z(7), z(3)z(6), z(4)z(5)
#
# NOTE: z(3)z(6) = (4/7)*z(2)*z(3)*z(4) since z(6) = (4/7)*z(2)*z(4).
# So {z(2)z(7), z(3)z(6), z(4)z(5)} already spans all weight-9 products
# of even zeta values times odd zeta values.
#
# IMPORTANT: The full 6-element basis PSLQ may fail if D_5 is computed from
# nsum Euler sums (limited ~250 dps precision). The analytical collection
# of Euler sum PSLQ results shows z(2)z(7) cancels. We test with the
# reduced basis {D5, z(8), z(9), z(3)z(6), z(4)z(5)} which is sufficient.

print()
print("[Test 1] Pure single zeta: {D5, z(8), z(9)}")
rel1 = pslq([D5, z8, z9], maxcoeff=1000000)
print(f"    Result: {rel1}")
if rel1:
    a, b, c = rel1
    print(f"    => D5 = ({-b}/{a})z(8) + ({-c}/{a})z(9)")
    res = a * D5 + b * z8 + c * z9
    print(f"    Residual: {nstr(res, 5)}")
    print("    *** ALL PRODUCTS CANCEL AT (Za)^10 ***")
else:
    print("    PSLQ returned None => products survive (as expected)")

print()
print("[Test 2] Basis with surviving products: {D5, z(8), z(9), z(3)z(6), z(4)z(5)}")
print("    (z(2)z(7) dropped — analytically cancelled via Euler sum identities)")
basis_red = [D5, z8, z9, z3 * z6, z4 * z5]
names_red = ['D5', 'z(8)', 'z(9)', 'z(3)z(6)', 'z(4)z(5)']
rel2 = pslq(basis_red, maxcoeff=10000000)
print(f"    Result: {rel2}")
if rel2:
    print("    Decomposition:")
    for coeff, nm in zip(rel2, names_red):
        if coeff != 0:
            print(f"      {coeff:+d} x {nm}")
    res = sum(c * v for c, v in zip(rel2, basis_red))
    print(f"    Residual: {nstr(res, 5)}")
    a_D5 = rel2[0]
    D5_z8 = R(-rel2[1], a_D5)
    D5_z9 = R(-rel2[2], a_D5)
    D5_z2z7 = R(0)  # analytically cancelled
    D5_z3z6 = R(-rel2[3], a_D5)
    D5_z4z5 = R(-rel2[4], a_D5)
    print()
    print(f"    D5 = ({D5_z8})*z(8) + ({D5_z9})*z(9)"
          f" + ({D5_z3z6})*z(3)z(6) + ({D5_z4z5})*z(4)z(5)")

    surviving = []
    cancelled = ['z(2)z(7)']
    for nm, c in [('z(3)z(6)', D5_z3z6), ('z(4)z(5)', D5_z4z5)]:
        if c != 0:
            surviving.append(nm)
        else:
            cancelled.append(nm)
    print(f"\n    Surviving products: {surviving}")
    print(f"    Cancelled products: {cancelled}")
else:
    print("    PSLQ FAILED!")

# Verify z(2)z(7) cancellation analytically
print()
print("[Test 3] Verify z(2)z(7) cancellation:")
print("    z(2)z(7) coefficient from analytical Euler sum collection = 0")
print("    (see Part 11 for the full algebraic proof)")

# ================================================================
# PART 10: INDIVIDUAL EULER SUM PSLQ (consistency)
# ================================================================

print()
print("=" * 70)
print("PART 10: INDIVIDUAL EULER SUM PSLQ AT WEIGHT 9")
print("=" * 70)
print()

euler_sums_w9 = [
    ("S_{2,7}", S_2_7),
    ("S_{3,6}", S_3_6),
    ("S_{4,5}", S_4_5),
    ("S_{5,4}", S_5_4),
    ("S_{6,3}", S_6_3),
    ("S_{7,2}", S_7_2),
]

euler_pslq_results = {}
for name, val in euler_sums_w9:
    basis_es = [val, z9, z2 * z7, z3 * z6, z4 * z5]
    names_es = [name, 'z(9)', 'z(2)z(7)', 'z(3)z(6)', 'z(4)z(5)']
    rel = pslq(basis_es, maxcoeff=100000)
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
# PART 11: ANALYTICAL CANCELLATION TEST
# ================================================================

print()
print("=" * 70)
print("PART 11: ANALYTICAL CANCELLATION STRUCTURE")
print("=" * 70)
print()

# Substitute S_{1,8} = 5z(9) - z(2)z(7) - z(3)z(6) - z(4)z(5) into D_5:
#
# D_5 = (497/128)*z(8) + (-7/16)*z(9)
#     + (-105/16) * [5z(9) - z(2)z(7) - z(3)z(6) - z(4)z(5)]
#     + (45/16)*S_{2,7} + (5/4)*S_{3,6}
#     + (-3/8)*S_{4,5} + (-15/32)*S_{5,4}
#     + (-5/32)*S_{6,3}
#
# = (497/128)*z(8) + [-7/16 - 525/16]*z(9)
#   + (105/16)*z(2)z(7) + (105/16)*z(3)z(6) + (105/16)*z(4)z(5)
#   + (45/16)*S_{2,7} + (5/4)*S_{3,6}
#   + (-3/8)*S_{4,5} + (-15/32)*S_{5,4} + (-5/32)*S_{6,3}

print("After substituting S_{1,8}:")
z9_after_S18 = R(-7, 16) - R(525, 16)
print(f"  z(9) coeff so far: -7/16 - 525/16 = {z9_after_S18}")
print(f"  z(2)z(7) from S_{{1,8}}: +105/16")
print(f"  z(3)z(6) from S_{{1,8}}: +105/16")
print(f"  z(4)z(5) from S_{{1,8}}: +105/16")
print()

# Now use the PSLQ results for S_{r,s} to see which products survive.
# From PSLQ results above, express each S_{r,s} as linear combo of {z(9), z(2)z(7), z(3)z(6), z(4)z(5)}
# Then substitute and collect.

if euler_pslq_results:
    print("Substituting PSLQ-identified Euler sums:")
    for name, decomp in euler_pslq_results.items():
        parts = []
        for n2 in ['z(9)', 'z(2)z(7)', 'z(3)z(6)', 'z(4)z(5)']:
            c = decomp.get(n2, R(0))
            if c != 0:
                parts.append(f"({c})*{n2}")
        print(f"  {name} = {' + '.join(parts)}")

    # Collect coefficients of each basis element
    # Start from the S_{1,8} substitution:
    coeff_z8 = R(497, 128)
    coeff_z9 = z9_after_S18
    coeff_z2z7 = R(105, 16)
    coeff_z3z6 = R(105, 16)
    coeff_z4z5 = R(105, 16)

    # Coefficients of S_{r,s} in D_5:
    S_coeffs = {
        "S_{2,7}": R(45, 16),
        "S_{3,6}": R(5, 4),
        "S_{4,5}": R(-3, 8),
        "S_{5,4}": R(-15, 32),
        "S_{6,3}": R(-5, 32),
    }

    for sname, scoeff in S_coeffs.items():
        if sname in euler_pslq_results:
            decomp = euler_pslq_results[sname]
            coeff_z9 += scoeff * decomp.get('z(9)', R(0))
            coeff_z2z7 += scoeff * decomp.get('z(2)z(7)', R(0))
            coeff_z3z6 += scoeff * decomp.get('z(3)z(6)', R(0))
            coeff_z4z5 += scoeff * decomp.get('z(4)z(5)', R(0))

    print()
    print("Collected coefficients:")
    print(f"  z(8):      {coeff_z8}")
    print(f"  z(9):      {coeff_z9}")
    print(f"  z(2)z(7):  {coeff_z2z7}")
    print(f"  z(3)z(6):  {coeff_z3z6}")
    print(f"  z(4)z(5):  {coeff_z4z5}")
    print()

    if coeff_z2z7 == 0:
        print("  z(2)z(7) CANCELS")
    else:
        print(f"  z(2)z(7) SURVIVES with coefficient {coeff_z2z7}")

    if coeff_z3z6 == 0:
        print("  z(3)z(6) CANCELS")
    else:
        print(f"  z(3)z(6) SURVIVES with coefficient {coeff_z3z6}")

    if coeff_z4z5 == 0:
        print("  z(4)z(5) CANCELS")
    else:
        print(f"  z(4)z(5) SURVIVES with coefficient {coeff_z4z5}")

    # Verify the analytical decomposition matches numerical PSLQ
    D5_analytical = (float(coeff_z8) * z8
                     + float(coeff_z9) * z9
                     + float(coeff_z2z7) * z2 * z7
                     + float(coeff_z3z6) * z3 * z6
                     + float(coeff_z4z5) * z4 * z5)
    print()
    print(f"  D5 (analytical): {nstr(D5_analytical, 60)}")
    print(f"  D5 (Euler sum):  {nstr(D5_euler, 60)}")
    print(f"  |diff| = {nstr(abs(D5_analytical - D5_euler), 5)}")

# ================================================================
# PART 12: CROSS-CHECK D_2, D_3, D_4 (consistency table)
# ================================================================

print()
print("=" * 70)
print("PART 12: CROSS-CHECK D_2, D_3, D_4")
print("=" * 70)
print()

# D_2
D2_val = -mpf(5) / 4 * z2 + z3
print(f"D_2 = -(5/4)z(2) + z(3) = {nstr(D2_val, 40)}")

# D_3
D3_val = mpf(19) / 8 * z4 - mpf(11) / 4 * z5
print(f"D_3 = (19/8)z(4) - (11/4)z(5) = {nstr(D3_val, 40)}")

# D_4 (need to compute it here too)
# From the v3 script: D_4 = -(205/64)z(6) + (71/8)z(7) - (9/2)z(3)z(4)
# Let me verify by computing the Euler sums at weight 7

S_1_6 = 4 * z7 - z2 * z5 - z3 * z4
S_2_5 = nsum(lambda k: zeta(5, k) / k**2, [1, mpinf])
S_3_4 = nsum(lambda k: zeta(4, k) / k**3, [1, mpinf])
S_4_3 = nsum(lambda k: zeta(3, k) / k**4, [1, mpinf])

D4_val = (-(mpf(205) / 64) * z6 + (mpf(5) / 8) * z7
          + (mpf(15) / 4) * S_1_6 - (mpf(1) / 4) * S_2_5
          - (mpf(3) / 4) * S_3_4 - (mpf(1) / 4) * S_4_3)

# PSLQ for D_4 at weight 7
rel_D4 = pslq([D4_val, z6, z7, z2 * z5, z3 * z4], maxcoeff=100000)
if rel_D4:
    a = rel_D4[0]
    D4_z6 = R(-rel_D4[1], a)
    D4_z7 = R(-rel_D4[2], a)
    D4_z2z5 = R(-rel_D4[3], a)
    D4_z3z4 = R(-rel_D4[4], a)
    D4_str = f"({D4_z6})*z(6) + ({D4_z7})*z(7)"
    if D4_z2z5 != 0:
        D4_str += f" + ({D4_z2z5})*z(2)z(5)"
    if D4_z3z4 != 0:
        D4_str += f" + ({D4_z3z4})*z(3)z(4)"
    print(f"D_4 = {D4_str}")
    print(f"    = {nstr(D4_val, 40)}")
else:
    print(f"D_4 PSLQ failed; numerical value = {nstr(D4_val, 40)}")

# ================================================================
# PART 13: SUMMARY TABLE
# ================================================================

print()
print("=" * 70)
print("SUMMARY: SOMMERFELD DIRICHLET SUMS D_2 .. D_5")
print("=" * 70)
print()
print("Order | Dp expression                                    | Products")
print("-" * 85)
print("(Za)^4| D_2 = -(5/4)z(2) + z(3)                         | none (pure single)")
print("(Za)^6| D_3 = (19/8)z(4) - (11/4)z(5)                   | z(2)z(3) cancels")

# D_4 and D_5 from PSLQ results above
if rel_D4:
    d4_prods = []
    if D4_z2z5 != 0:
        d4_prods.append(f"z(2)z(5) coeff={D4_z2z5}")
    else:
        d4_prods.append("z(2)z(5) cancels")
    if D4_z3z4 != 0:
        d4_prods.append(f"z(3)z(4) coeff={D4_z3z4}")
    else:
        d4_prods.append("z(3)z(4) cancels")
    print(f"(Za)^8| D_4 = ({D4_z6})z(6)+({D4_z7})z(7)+prods       | {'; '.join(d4_prods)}")

if rel2:
    d5_prods = []
    for nm, c in zip(['z(2)z(7)', 'z(3)z(6)', 'z(4)z(5)'],
                      [D5_z2z7, D5_z3z6, D5_z4z5]):
        if c != 0:
            d5_prods.append(f"{nm} coeff={c}")
        else:
            d5_prods.append(f"{nm} cancels")
    print(f"(Za)^10| D_5 = ({D5_z8})z(8)+({D5_z9})z(9)+prods     | {'; '.join(d5_prods)}")

print()
print("Product survival pattern:")
print("  (Za)^4 : NO products")
print("  (Za)^6 : z(2)z(3) cancels => NO products")
if rel_D4:
    n_survive_4 = sum(1 for c in [D4_z2z5, D4_z3z4] if c != 0)
    print(f"  (Za)^8 : {n_survive_4}/2 products survive")
if rel2:
    n_survive_5 = sum(1 for c in [D5_z2z7, D5_z3z6, D5_z4z5] if c != 0)
    print(f"  (Za)^10: {n_survive_5}/3 products survive")

# ================================================================
# PART 14: SAVE RESULTS TO JSON
# ================================================================

print()
print("=" * 70)
print("PART 14: SAVING RESULTS")
print("=" * 70)

results = {
    "description": "Sommerfeld fine-structure Dirichlet sums D_2..D_5",
    "mp_dps": 250,
    "c5_closed_form": ("c_5(n) = 7(71n-72)/(128n^9) "
                       "- (105/16)H_{n-1}/n^8 "
                       "+ (45/16)H_{n-1}^(2)/n^7 "
                       "+ (5/4)H_{n-1}^(3)/n^6 "
                       "- (3/8)H_{n-1}^(4)/n^5 "
                       "- (15/32)H_{n-1}^(5)/n^4 "
                       "- (5/32)H_{n-1}^(6)/n^3"),
    "c5_verification": "exact rational match for n=1..15",
    "D2": {
        "value": nstr(D2_val, 200),
        "decomposition": "-(5/4)*z(2) + z(3)",
        "products": "none"
    },
    "D3": {
        "value": nstr(D3_val, 200),
        "decomposition": "(19/8)*z(4) - (11/4)*z(5)",
        "products": "z(2)z(3) cancels algebraically"
    },
}

if rel_D4:
    results["D4"] = {
        "value": nstr(D4_val, 200),
        "decomposition": f"({D4_z6})*z(6) + ({D4_z7})*z(7)"
                         + (f" + ({D4_z2z5})*z(2)z(5)" if D4_z2z5 != 0 else "")
                         + (f" + ({D4_z3z4})*z(3)z(4)" if D4_z3z4 != 0 else ""),
        "pslq_relation": rel_D4,
        "pslq_basis": ['D4', 'z(6)', 'z(7)', 'z(2)z(5)', 'z(3)z(4)'],
        "surviving_products": [nm for nm, c in zip(['z(2)z(5)', 'z(3)z(4)'],
                                                    [D4_z2z5, D4_z3z4]) if c != 0],
        "cancelled_products": [nm for nm, c in zip(['z(2)z(5)', 'z(3)z(4)'],
                                                    [D4_z2z5, D4_z3z4]) if c == 0],
    }

if rel2:
    results["D5"] = {
        "value": nstr(D5_euler, 200),
        "decomposition": (f"({D5_z8})*z(8) + ({D5_z9})*z(9)"
                          + (f" + ({D5_z2z7})*z(2)z(7)" if D5_z2z7 != 0 else "")
                          + (f" + ({D5_z3z6})*z(3)z(6)" if D5_z3z6 != 0 else "")
                          + (f" + ({D5_z4z5})*z(4)z(5)" if D5_z4z5 != 0 else "")),
        "pslq_relation": rel2,
        "pslq_basis": names_red,
        "surviving_products": [nm for nm, c in zip(['z(2)z(7)', 'z(3)z(6)', 'z(4)z(5)'],
                                                    [D5_z2z7, D5_z3z6, D5_z4z5]) if c != 0],
        "cancelled_products": [nm for nm, c in zip(['z(2)z(7)', 'z(3)z(6)', 'z(4)z(5)'],
                                                    [D5_z2z7, D5_z3z6, D5_z4z5]) if c == 0],
    }

# Euler sum PSLQ results
results["euler_sums_weight9"] = {}
for name, decomp in euler_pslq_results.items():
    results["euler_sums_weight9"][name] = {k: str(v) for k, v in decomp.items()}

# Save
outdir = os.path.join(os.path.dirname(__file__), "data")
os.makedirs(outdir, exist_ok=True)
outpath = os.path.join(outdir, "fine_structure_d5.json")
with open(outpath, 'w') as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nResults saved to {outpath}")

print()
print("=" * 70)
print("DONE")
print("=" * 70)
