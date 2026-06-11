"""
Fine-structure Euler sum cancellation: ANALYTICAL derivation at (Za)^8.

The per-state coefficient a_p(n,jph) has been derived from sympy. The degeneracy
weighted sum c_p(n) = Sum_j d(n,j) * a_p(n,j) can be computed ANALYTICALLY because
j+1/2 ranges over 1,2,...,n and the degeneracy is 4*jph for jph<n, 2n for jph=n.

This gives c_p(n) in terms of polynomial sums (closed-form) and harmonic numbers.
"""
import sympy as sp
from sympy import Rational as R
from mpmath import mp, mpf, nstr, zeta, pslq
import time

mp.dps = 120

# ================================================================
# ANALYTICAL DERIVATION of c_4(n)
# ================================================================

# From sympy:
# a_4(n,k) = (35k^5 - 120k^4*n + 120k^3*n^2 - 8k^2*n^3 - 24k*n^4 - 8n^5)
#           / (128 * k^5 * n^8)
#
# Degeneracy: d(n,j) = 4*jph for jph < n;  d(n,n) = 2n   [jph = j+1/2 = 1..n]
#
# c_4(n) = Sum_{k=1}^{n-1} 4k * a_4(n,k) + 2n * a_4(n,n)
#
# 4k * a_4(n,k) = [35k - 120n + 120n^2/k - 8n^3/k^2 - 24n^4/k^3 - 8n^5/k^4]
#                 / (32 * n^8)
#
# Summing over k=1..n-1:
#   Sum 35k = 35*n*(n-1)/2
#   Sum -120n = -120n*(n-1)
#   Sum 120n^2/k = 120n^2 * H_{n-1}
#   Sum -8n^3/k^2 = -8n^3 * H_{n-1}^{(2)}
#   Sum -24n^4/k^3 = -24n^4 * H_{n-1}^{(3)}
#   Sum -8n^5/k^4 = -8n^5 * H_{n-1}^{(4)}
#
# Diagonal: a_4(n,n) = (35-120+120-8-24-8)/(128*n^8) = -5/(128*n^8)
# So 2n * a_4(n,n) = -5/(64*n^7)
#
# Collecting:
# c_4(n) = [35n(n-1)/2 - 120n(n-1)]/(32n^8) - 5/(64n^7)
#         + (15/4) H_{n-1}/n^6
#         - (1/4) H_{n-1}^{(2)}/n^5
#         - (3/4) H_{n-1}^{(3)}/n^4
#         - (1/4) H_{n-1}^{(4)}/n^3
#
# Rational part:
#   [35(n-1)/2 - 120(n-1)]/(32n^7) - 5/(64n^7)
#   = (n-1)(-205/2)/(32n^7) - 5/(64n^7)
#   = [-205(n-1)]/(64n^7) - 5/(64n^7)
#   = [-205n + 205 - 5]/(64n^7)
#   = -5(41n - 40)/(64n^7)

print("="*70)
print("VERIFYING c_4(n) CLOSED FORM")
print("="*70)
print()
print("c_4(n) = -5(41n-40)/(64n^7) + (15/4)H_{n-1}/n^6")
print("         - (1/4)H_{n-1}^{(2)}/n^5 - (3/4)H_{n-1}^{(3)}/n^4")
print("         - (1/4)H_{n-1}^{(4)}/n^3")
print()

# Verify against sympy for n=1..10
a2 = sp.Symbol('a2', positive=True)
jph_s = sp.Symbol('jph', positive=True)
n_s = sp.Symbol('n', positive=True, integer=True)

gamma_s = sp.sqrt(jph_s**2 - a2)
delta_s = jph_s - gamma_s
n_eff = n_s - delta_s
E_sym = 1/sp.sqrt(1 + a2/n_eff**2) - 1
E_expanded = sp.series(E_sym, a2, 0, n=6).removeO()
a4_sym = sp.cancel(E_expanded.coeff(a2, 4))

def degeneracy(n_val, j2):
    j = R(j2, 2)
    j_max = n_val - R(1, 2)
    if j == j_max:
        return int(2*j + 1)
    else:
        return int(2*(2*j + 1))

for n_val in range(1, 11):
    # Direct from sympy
    total_sympy = R(0)
    for j2 in range(1, 2*n_val, 2):
        j_val = R(j2, 2)
        d = degeneracy(n_val, j2)
        val = a4_sym.subs([(n_s, n_val), (jph_s, j_val + R(1,2))])
        total_sympy += d * sp.cancel(val)
    total_sympy = sp.cancel(total_sympy)

    # From closed form
    H1 = sum(R(1,k) for k in range(1, n_val))
    H2 = sum(R(1,k**2) for k in range(1, n_val))
    H3 = sum(R(1,k**3) for k in range(1, n_val))
    H4 = sum(R(1,k**4) for k in range(1, n_val))

    formula = (R(-5*(41*n_val-40), 64*n_val**7)
               + R(15,4)*H1/n_val**6
               - R(1,4)*H2/n_val**5
               - R(3,4)*H3/n_val**4
               - R(1,4)*H4/n_val**3)

    diff = sp.cancel(total_sympy - formula)
    status = "OK" if diff == 0 else f"FAIL (diff={diff})"
    print(f"  n={n_val}: sympy={total_sympy}, formula={float(formula):.15e}, {status}")

# ================================================================
# DIRICHLET SUM D_4 decomposition
# ================================================================

print()
print("="*70)
print("DIRICHLET SUM D_4 = Sum c_4(n)")
print("="*70)

# D_4 = Sum c_4(n) = R_part + H1_part + H2_part + H3_part + H4_part
#
# R_part = -5/64 * Sum (41n-40)/n^7 = -5/64 * [41*z(6) - 40*z(7)]
#
# H1_part = 15/4 * Sum H_{n-1}/n^6 = 15/4 * [S_{1,6} - z(7)]
#   (because Sum H_{n-1}/n^s = Sum [H_n - 1/n]/n^s = S_{1,s} - z(s+1))
#   Wait: H_{n-1} = H_n - 1/n, so Sum H_{n-1}/n^s = S_{1,s} - z(s+1).
#   Actually: Sum_{n>=1} H_{n-1}/n^s. At n=1, H_0 = 0. At n>=2, H_{n-1} = H_{n} - 1/n.
#   Sum_{n>=1} H_{n-1}/n^s = Sum_{n>=2} (H_n - 1/n)/n^s = Sum_{n>=1} H_n/n^s - H_1/1^s - Sum_{n>=2} 1/n^{s+1}
#   Wait, this is getting confused. Let me be careful.
#
#   Sum_{n=1}^inf H_{n-1}/n^s where H_0 = 0, H_k = 1 + 1/2 + ... + 1/k
#   = Sum_{n=2}^inf H_{n-1}/n^s  (since H_0 = 0)
#   = Sum_{m=1}^inf H_m / (m+1)^s  [substitute m = n-1]
#
#   Alternative: H_{n-1} = H_n - 1/n for n >= 1 (with H_0 = 0 so H_{1-1} = H_1 - 1/1 = 0, OK)
#   Sum_{n=1}^inf H_{n-1}/n^s = Sum_{n=1}^inf (H_n - 1/n)/n^s = S_{1,s} - z(s+1)
#
#   S_{1,s} = Sum_{n=1}^inf H_n / n^s
#   So Sum H_{n-1}/n^s = S_{1,s} - z(s+1). Verified.
#
# Similarly: Sum H_{n-1}^{(r)}/n^s = S_{r,s} - z(r+s) for r >= 1
#   where S_{r,s} = Sum_{n=1}^inf H_n^{(r)} / n^s
#   and H_{n-1}^{(r)} = H_n^{(r)} - 1/n^r
#   Sum H_{n-1}^{(r)}/n^s = S_{r,s} - z(r+s). Verified for r >= 1.
#
# So:
# D_4 = -(205/64)z(6) + (200/64)z(7)
#      + 15/4 * [S_{1,6} - z(7)]
#      - 1/4 * [S_{2,5} - z(7)]
#      - 3/4 * [S_{3,4} - z(7)]
#      - 1/4 * [S_{4,3} - z(7)]
#
# Collecting z(7) terms:
# (200/64) - 15/4 + 1/4 + 3/4 + 1/4 = 25/8 - 15/4 + 5/4 = 25/8 - 10/4 = 25/8 - 20/8 = 5/8

print()
print("D_4 = -(205/64)z(6) + (5/8)z(7)")
print("    + (15/4)S_{1,6} - (1/4)S_{2,5} - (3/4)S_{3,4} - (1/4)S_{4,3}")
print()

# ================================================================
# EULER SUM EVALUATIONS at weight 7
# ================================================================

# S_{1,6}: Euler's formula for even s:
# 2*S_{1,2k} = (2k+2)z(2k+1) - Sum_{j=1}^{2k-2} z(j+1)*z(2k-j)
# For k=3: 2*S_{1,6} = 8*z(7) - [z(2)z(5) + z(3)z(4) + z(4)z(3) + z(5)z(2)]
# = 8*z(7) - 2*z(2)z(5) - 2*z(3)z(4)
# S_{1,6} = 4*z(7) - z(2)z(5) - z(3)z(4)

print("Known evaluations:")
print("  S_{1,6} = 4z(7) - z(2)z(5) - z(3)z(4)  [Euler's formula]")

# For S_{2,5}, S_{3,4}, S_{4,3} we need additional identities.
# Constraint: S_{r,s} + S_{s,r} = z(r)z(s) + z(r+s)
#   S_{2,5} + S_{5,2} = z(2)z(5) + z(7)
#   S_{3,4} + S_{4,3} = z(3)z(4) + z(7)

# ================================================================
# HIGH-PRECISION NUMERICAL COMPUTATION of Euler sums
# ================================================================

print()
print("Computing Euler sums S_{r,s} numerically at 120 dps...")

from mpmath import nsum, inf as mpinf

def euler_sum_hurwitz(r, s):
    """Compute S_{r,s} = Sum_{n=1}^inf H_n^{(r)} / n^s via Hurwitz zeta.
    Uses identity: S_{r,s} = Sum_{k=1}^inf zeta(s, k) / k^r
    where zeta(s, k) = Sum_{n>=k} 1/n^s is the Hurwitz zeta.
    """
    t0 = time.time()
    val = nsum(lambda k: zeta(s, k) / k**r, [1, mpinf])
    dt = time.time() - t0
    return val, dt

z2 = zeta(2); z3 = zeta(3); z4 = zeta(4); z5 = zeta(5); z6 = zeta(6); z7 = zeta(7)

# S_{1,6} from Euler's formula (exact)
S_1_6 = 4*z7 - z2*z5 - z3*z4
print(f"  S_{{1,6}} = {nstr(S_1_6, 60)} [Euler formula]")

S_2_5, dt = euler_sum_hurwitz(2, 5)
print(f"  S_{{2,5}} = {nstr(S_2_5, 60)} ({dt:.1f}s)")

S_3_4, dt = euler_sum_hurwitz(3, 4)
print(f"  S_{{3,4}} = {nstr(S_3_4, 60)} ({dt:.1f}s)")

S_4_3, dt = euler_sum_hurwitz(4, 3)
print(f"  S_{{4,3}} = {nstr(S_4_3, 60)} ({dt:.1f}s)")

S_5_2, dt = euler_sum_hurwitz(5, 2)
print(f"  S_{{5,2}} = {nstr(S_5_2, 60)} ({dt:.1f}s)")

# Verify stuffle identities (z2..z7 defined above)
print(f"\nStuffle checks:")
print(f"  S_{{2,5}}+S_{{5,2}} - z(2)z(5) - z(7) = {nstr(S_2_5+S_5_2-z2*z5-z7, 5)}")
print(f"  S_{{3,4}}+S_{{4,3}} - z(3)z(4) - z(7) = {nstr(S_3_4+S_4_3-z3*z4-z7, 5)}")

# Verify S_{1,6} against Euler formula
S_1_6_euler = 4*z7 - z2*z5 - z3*z4
print(f"  S_{{1,6}} - [4z(7)-z(2)z(5)-z(3)z(4)] = {nstr(S_1_6-S_1_6_euler, 5)}")

# ================================================================
# PSLQ identification of S_{r,s} in terms of {z(7), z(2)z(5), z(3)z(4)}
# ================================================================

print()
print("="*70)
print("PSLQ IDENTIFICATION of Euler sums")
print("="*70)

for name, val in [("S_{2,5}", S_2_5), ("S_{3,4}", S_3_4), ("S_{4,3}", S_4_3), ("S_{5,2}", S_5_2)]:
    rel = pslq([val, z7, z2*z5, z3*z4])
    if rel:
        a, b, c, d = rel
        print(f"  {name} = ({-b}/{a})z(7) + ({-c}/{a})z(2)z(5) + ({-d}/{a})z(3)z(4)")
        res = a*val + b*z7 + c*z2*z5 + d*z3*z4
        print(f"    Residual: {nstr(res, 5)}")
    else:
        print(f"  {name}: PSLQ failed!")

# ================================================================
# D_4 ASSEMBLY
# ================================================================

print()
print("="*70)
print("D_4 ASSEMBLY AND CANCELLATION TEST")
print("="*70)

D4 = -(mpf(205)/64)*z6 + (mpf(5)/8)*z7 + (mpf(15)/4)*S_1_6 - (mpf(1)/4)*S_2_5 - (mpf(3)/4)*S_3_4 - (mpf(1)/4)*S_4_3

print(f"\nD_4 = {nstr(D4, 60)}")

# PSLQ: reduced basis (cancellation hypothesis)
print(f"\n[1] Reduced basis {{D4, z(6), z(7)}}:")
rel1 = pslq([D4, z6, z7], maxcoeff=10000)
print(f"    Result: {rel1}")
if rel1:
    a, b, c = rel1
    print(f"    => D4 = ({-b}/{a})z(6) + ({-c}/{a})z(7)")
    res = a*D4 + b*z6 + c*z7
    print(f"    Residual: {nstr(res, 5)}")
    print(f"    *** PRODUCTS CANCEL AT (Za)^8 ***")

# PSLQ: full weight-7 basis
print(f"\n[2] Full basis {{D4, z(6), z(7), z(2)z(5), z(3)z(4)}}:")
rel2 = pslq([D4, z6, z7, z2*z5, z3*z4], maxcoeff=100000)
print(f"    Result: {rel2}")
if rel2:
    names = ['D4', 'z(6)', 'z(7)', 'z(2)z(5)', 'z(3)z(4)']
    print("    Decomposition:")
    for coeff, nm in zip(rel2, names):
        if coeff != 0:
            print(f"      {coeff:+d} x {nm}")
    res = sum(c*v for c,v in zip(rel2, [D4, z6, z7, z2*z5, z3*z4]))
    print(f"    Residual: {nstr(res, 5)}")
    prod_coeffs = rel2[3:]
    if all(c == 0 for c in prod_coeffs):
        print("    *** ALL PRODUCT COEFFICIENTS ZERO => CANCELLATION CONFIRMED ***")
    else:
        print("    *** PRODUCTS PRESENT => CANCELLATION FAILS ***")

# ================================================================
# CROSS-CHECK D_3 at (Za)^6
# ================================================================

print()
print("="*70)
print("CROSS-CHECK: D_3 at (Za)^6")
print("="*70)

# D_3 = -(19/8)(5/4) wait...
# c_3(n) = (19n-20)/(8n^5) - (3/2)H_{n-1}/n^4 - (1/2)H_{n-1}^{(2)}/n^3
# D_3 = Sum c_3(n) = (19/8)z(4) - (20/8)z(5) - (3/2)[S_{1,4} - z(5)] - (1/2)[S_{2,3} - z(5)]
# = (19/8)z(4) - (5/2)z(5) - (3/2)S_{1,4} + (3/2)z(5) - (1/2)S_{2,3} + (1/2)z(5)
# = (19/8)z(4) + [-5/2+3/2+1/2]z(5) - (3/2)S_{1,4} - (1/2)S_{2,3}
# = (19/8)z(4) - (1/2)z(5) - (3/2)S_{1,4} - (1/2)S_{2,3}

# S_{1,4}: Euler formula for k=2: 2S_{1,4} = 6z(5) - [z(2)z(3) + z(3)z(2)] = 6z(5) - 2z(2)z(3)
# S_{1,4} = 3z(5) - z(2)z(3)

# S_{2,3}: stuffle: S_{2,3} + S_{3,2} = z(2)z(3) + z(5)
# From the literature: S_{2,3} = 3z(2)z(3) - (9/2)z(5)
# Check stuffle: S_{3,2} = z(2)z(3) + z(5) - S_{2,3} = z(2)z(3) + z(5) - 3z(2)z(3) + (9/2)z(5)
#              = -2z(2)z(3) + (11/2)z(5)
# This can be verified independently.

# Note: S_{2,3} here uses H_n convention. But we need H_{n-1} convention.
# Sum H_{n-1}^{(2)}/n^3 = S_{2,3} - z(5) = 3z(2)z(3) - (9/2)z(5) - z(5) = 3z(2)z(3) - (11/2)z(5)

S_1_4 = mpf(3)*z5 - z2*z3
S_2_3 = mpf(3)*z2*z3 - (mpf(9)/2)*z5

D3_analytical = (mpf(19)/8)*z4 - (mpf(1)/2)*z5 - (mpf(3)/2)*S_1_4 - (mpf(1)/2)*S_2_3
print(f"D_3 (analytical) = {nstr(D3_analytical, 50)}")

# Expand:
# = (19/8)z4 - (1/2)z5 - (3/2)[3z5 - z2z3] - (1/2)[3z2z3 - (9/2)z5]
# = (19/8)z4 - (1/2)z5 - (9/2)z5 + (3/2)z2z3 - (3/2)z2z3 + (9/4)z5
# = (19/8)z4 + [-1/2 - 9/2 + 9/4]z5 + [3/2 - 3/2]z2z3
# = (19/8)z4 + [-2/4 - 18/4 + 9/4]z5
# = (19/8)z4 - (11/4)z5

D3_simplified = (mpf(19)/8)*z4 - (mpf(11)/4)*z5
print(f"D_3 (simplified) = (19/8)z(4) - (11/4)z(5) = {nstr(D3_simplified, 50)}")
print(f"|diff| = {nstr(abs(D3_analytical - D3_simplified), 5)}")

# Direct verification: sum c_3(n) to N=100000
D3_direct = mpf(0)
H1_inc = mpf(0)
H2_inc = mpf(0)
for n_val in range(1, 100001):
    nf = mpf(n_val)
    c3 = mpf(19*n_val - 20)/(8*nf**5) - mpf(3)/2 * H1_inc/nf**4 - mpf(1)/2 * H2_inc/nf**3
    D3_direct += c3
    H1_inc += mpf(1)/nf
    H2_inc += mpf(1)/nf**2

# Tail correction for D3: leading term is ~ -(205/8)/n^4 for large n
# c_3(n) ~ (19/8)/n^4 - (3/2)ln(n)/n^4 - ... for large n
# Actually the asymptotic is messy. Just use the Euler sum formula.
print(f"D_3 (direct N=100000) = {nstr(D3_direct, 30)}")
print(f"|simplified - direct| = {nstr(abs(D3_simplified - D3_direct), 5)}")
print(f"[Note: direct sum convergence is O(1/N^2), ~10 digits at N=100000]")

print(f"\n*** z(2)z(3) product coefficient in D_3: (3/2) - (3/2) = 0 ***")
print(f"*** CANCELLATION CONFIRMED ALGEBRAICALLY ***")

# ================================================================
# ANALYTICAL CANCELLATION TEST at (Za)^8
# ================================================================

print()
print("="*70)
print("ANALYTICAL CANCELLATION AT (Za)^8")
print("="*70)

# D_4 = -(205/64)z(6) + (5/8)z(7) + (15/4)S_{1,6} - (1/4)S_{2,5} - (3/4)S_{3,4} - (1/4)S_{4,3}
#
# Substituting S_{1,6} = 4z(7) - z(2)z(5) - z(3)z(4):
#
# D_4 = -(205/64)z(6) + (5/8)z(7) + 15z(7) - (15/4)z(2)z(5) - (15/4)z(3)z(4)
#       - (1/4)S_{2,5} - (3/4)S_{3,4} - (1/4)S_{4,3}
#
# = -(205/64)z(6) + (5/8 + 15)z(7) - (15/4)z(2)z(5) - (15/4)z(3)z(4)
#   - (1/4)S_{2,5} - (3/4)S_{3,4} - (1/4)S_{4,3}
#
# = -(205/64)z(6) + (125/8)z(7) - (15/4)z(2)z(5) - (15/4)z(3)z(4)
#   - (1/4)S_{2,5} - (3/4)S_{3,4} - (1/4)S_{4,3}

# Now use PSLQ results for S_{r,s} to express everything in terms of z(7), z(2)z(5), z(3)z(4).
# Then check if the z(2)z(5) and z(3)z(4) coefficients vanish.

print("\nSubstituting PSLQ-identified Euler sums and collecting product terms...")
print()
print("From PSLQ we need: S_{2,5}, S_{3,4}, S_{4,3} as Q-linear combos of {z(7), z(2)z(5), z(3)z(4)}")
print("[See PSLQ results above]")
print()
print("If the product coefficients cancel, we get:")
print("  D_4 = A*z(6) + B*z(7)  [pure single zetas]")
print()
print("Otherwise products persist.")

# ================================================================
# SUMMARY
# ================================================================

print()
print("="*70)
print("SUMMARY")
print("="*70)
print()
print("(Za)^4: D_2 = -(5/4)z(2) + z(3)  [pure single, PROVEN]")
print("(Za)^6: D_3 = (19/8)z(4) - (11/4)z(5)")
print("         z(2)z(3) coefficient = 3/2 - 3/2 = 0  [ALGEBRAIC PROOF]")
print("(Za)^8: D_4 = see PSLQ result above")
