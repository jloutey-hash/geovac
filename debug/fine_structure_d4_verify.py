"""
Independent high-precision verification of D_4 = Sum c_4(n).
Computes D_4 two ways and PSLQ-identifies against weight-7 MZV basis.
"""
from mpmath import mp, mpf, nstr, zeta, pslq, harmonic, nsum, inf as mpinf
import time

mp.dps = 200

z2 = zeta(2); z3 = zeta(3); z4 = zeta(4); z5 = zeta(5); z6 = zeta(6); z7 = zeta(7)

print("="*70)
print("D_4 INDEPENDENT COMPUTATION AT 200 DPS")
print("="*70)

# Method 1: direct incremental sum with Euler-Maclaurin tail
t0 = time.time()
N = 200000
H1 = mpf(0); H2 = mpf(0); H3 = mpf(0); H4 = mpf(0)
D4_partial = mpf(0)
for n in range(1, N+1):
    nf = mpf(n)
    c4 = (mpf(-5)*(41*n - 40))/(64*nf**7) + mpf(15)/4*H1/nf**6 \
         - mpf(1)/4*H2/nf**5 - mpf(3)/4*H3/nf**4 - mpf(1)/4*H4/nf**3
    D4_partial += c4
    H1 += 1/nf; H2 += 1/nf**2; H3 += 1/nf**3; H4 += 1/nf**4

# Tail: c_4(n) for large n ~ rational/n^6 + zeta(r)*coeff/n^s + corrections
# Rational tail: -205/(64n^6) + 25/(8n^7)
tail_rat = -mpf(205)/64 * zeta(6, N+1) + mpf(25)/8 * zeta(7, N+1)

# Harmonic tail: Sum_{n>N} (15/4)*H_{n-1}/n^6
# H_{n-1} = H_n - 1/n, and for large n, H_n ~ ln(n) + gamma + 1/(2n) - ...
# Sum_{n>N} H_{n-1}/n^6 = Sum_{n>N} H_n/n^6 - Sum_{n>N} 1/n^7
# S_{1,6} is exact, so: Sum_{n>N} H_n/n^6 = S_{1,6} - Sum_{n=1}^{N} H_n/n^6
# But this requires computing S_{1,6} - partial, which is circular...
# Instead: for r >= 2, Sum_{n>N} H_{n-1}^{(r)}/n^s
#   H_{n-1}^{(r)} = z(r) - zeta(r, n)
#   Sum_{n>N} [z(r) - zeta(r,n)]/n^s = z(r)*zeta(s,N+1) - Sum_{n>N} zeta(r,n)/n^s
# For the r=1 term we handle via subtraction from known total.

# Total D_4 from Euler sums (computed at 200 dps with nsum):
S_1_6 = 4*z7 - z2*z5 - z3*z4  # exact from Euler formula

# Compute S_{r,s} for r >= 2 via Hurwitz representation
print(f"\nComputing Euler sums via Hurwitz zeta at 200 dps...")
S_2_5 = nsum(lambda k: zeta(5, k) / k**2, [1, mpinf])
S_3_4 = nsum(lambda k: zeta(4, k) / k**3, [1, mpinf])
S_4_3 = nsum(lambda k: zeta(3, k) / k**4, [1, mpinf])
dt_es = time.time() - t0
print(f"  Euler sums computed in {dt_es:.1f}s")

# Method 2: D_4 from Euler sum decomposition
D4_euler = -(mpf(205)/64)*z6 + (mpf(5)/8)*z7 \
         + (mpf(15)/4)*S_1_6 - (mpf(1)/4)*S_2_5 \
         - (mpf(3)/4)*S_3_4 - (mpf(1)/4)*S_4_3

print(f"\nD_4 (direct N={N}):     {nstr(D4_partial, 50)}")
print(f"D_4 (Euler sums):       {nstr(D4_euler, 50)}")
print(f"|diff| = {nstr(abs(D4_euler - D4_partial), 5)}")
print(f"  [direct sum has O(1/N^2) tail error ~ {nstr(1/mpf(N)**2, 3)}]")

# Stuffle cross-checks
S_5_2 = nsum(lambda k: zeta(2, k) / k**5, [1, mpinf])
print(f"\nStuffle checks:")
print(f"  S(2,5)+S(5,2) - z(2)z(5) - z(7) = {nstr(S_2_5+S_5_2-z2*z5-z7, 5)}")
print(f"  S(3,4)+S(4,3) - z(3)z(4) - z(7) = {nstr(S_3_4+S_4_3-z3*z4-z7, 5)}")

# ================================================================
# PSLQ identification
# ================================================================
print()
print("="*70)
print("PSLQ AT 200 DPS")
print("="*70)

D4 = D4_euler

# Test 1: pure single zetas
print(f"\n[1] Reduced basis {{D4, z(6), z(7)}}:")
rel1 = pslq([D4, z6, z7], maxcoeff=100000)
print(f"    Result: {rel1}")
if rel1:
    a, b, c = rel1
    print(f"    D4 = ({-b}/{a})z(6) + ({-c}/{a})z(7)")

# Test 2: full weight-7 basis
print(f"\n[2] Full basis {{D4, z(6), z(7), z(2)z(5), z(3)z(4)}}:")
rel2 = pslq([D4, z6, z7, z2*z5, z3*z4], maxcoeff=100000)
print(f"    Result: {rel2}")
if rel2:
    a = rel2[0]
    names = ['D4', 'z(6)', 'z(7)', 'z(2)z(5)', 'z(3)z(4)']
    print("    Decomposition:")
    for coeff, nm in zip(rel2, names):
        if coeff != 0:
            print(f"      {coeff:+d} x {nm}")
    res = sum(c*v for c,v in zip(rel2, [D4, z6, z7, z2*z5, z3*z4]))
    print(f"    Residual: {nstr(res, 5)}")
    prod_terms = rel2[3:]
    if all(c == 0 for c in prod_terms):
        print("    *** CANCELLATION: D_4 = pure single zetas ***")
    else:
        print("    *** PRODUCTS SURVIVE at (Za)^8 ***")
        D4_z6_coeff = -rel2[1] / rel2[0]
        D4_z7_coeff = -rel2[2] / rel2[0]
        D4_z2z5_coeff = -rel2[3] / rel2[0]
        D4_z3z4_coeff = -rel2[4] / rel2[0]
        print(f"    D4 = ({D4_z6_coeff})*z(6) + ({D4_z7_coeff})*z(7) + ({D4_z2z5_coeff})*z(2)z(5) + ({D4_z3z4_coeff})*z(3)z(4)")

# Test 3: individual Euler sum PSLQ (verify consistency)
print()
print("="*70)
print("INDIVIDUAL EULER SUM PSLQ (consistency check)")
print("="*70)

for name, val in [("S(2,5)", S_2_5), ("S(3,4)", S_3_4), ("S(4,3)", S_4_3), ("S(5,2)", S_5_2)]:
    rel = pslq([val, z7, z2*z5, z3*z4], maxcoeff=100000)
    if rel:
        a, b, c, d = rel
        z7c = mpf(-b)/a; z2z5c = mpf(-c)/a; z3z4c = mpf(-d)/a
        print(f"  {name} = {z7c}*z(7) + {z2z5c}*z(2)z(5) + {z3z4c}*z(3)z(4)")
        res = a*val + b*z7 + c*z2*z5 + d*z3*z4
        print(f"    residual: {nstr(res, 5)}")
    else:
        print(f"  {name}: PSLQ FAILED (irreducible at weight 7?)")

# Cross-check stuffle with PSLQ results
print()
print("If PSLQ for S(2,5) and S(5,2) both succeeded, check stuffle consistency:")
