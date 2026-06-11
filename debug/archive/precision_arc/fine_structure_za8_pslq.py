"""
Fine structure (Zα)^8 Euler sum cancellation test.

Determines whether products ζ(even)×ζ(odd) cancel in D_{S^(8)}(8),
extending the pattern from (Zα)^6. Uses high-precision depth-2 MZV
computation via cumulative harmonic numbers + Euler-Maclaurin tail.
"""
from mpmath import mp, mpf, zeta, pi, pslq, nstr, fraction
import time

mp.dps = 120
N_MZV = 50000

def depth2_mzv(s, r, N):
    """
    Compute ζ(s,r) = Σ_{n>k≥1} 1/(n^s k^r)
    = Σ_{k=1}^∞ [ζ(s) - H_k^{(s)}] / k^r

    Uses partial sum to N + 2-term Euler-Maclaurin tail correction.
    Error: O(N^{-(s+r+1)}).
    """
    zs = zeta(s)
    H = mpf(0)
    partial = mpf(0)
    for k in range(1, N+1):
        H += mpf(1) / mpf(k)**s
        partial += (zs - H) / mpf(k)**r
    # Tail: Σ_{k>N} ζ(s,k+1)/k^r, with m=k+1 substitution
    # ≈ 1/(s-1)*ζ(s+r-1, N+2) + [r/(s-1)+1/2]*ζ(s+r, N+2)
    t1 = mpf(1)/(s-1) * zeta(s+r-1, N+2)
    t2 = (mpf(r)/(s-1) + mpf(1)/2) * zeta(s+r, N+2)
    return partial + t1 + t2

# ============================================================
# Sanity check: ζ(2,1) = ζ(3)
# ============================================================
print("Sanity check: ζ(2,1) should equal ζ(3)")
t0 = time.time()
z21 = depth2_mzv(2, 1, 5000)
z3 = zeta(3)
print(f"  ζ(2,1) = {nstr(z21, 40)}")
print(f"  ζ(3)   = {nstr(z3, 40)}")
print(f"  diff   = {nstr(z21 - z3, 5)} ({time.time()-t0:.1f}s)")
print()

# ============================================================
# (Zα)^6 verification: D_{S^(6)}(6) products should cancel
# ============================================================
print("="*60)
print("(Zα)^6 verification (weight 11)")
print("="*60)

# Need ζ(9,2) for S_{2,9}
t0 = time.time()
z92 = depth2_mzv(9, 2, N_MZV)
print(f"ζ(9,2) = {nstr(z92, 50)} ({time.time()-t0:.1f}s)")

S_2_9 = z92 + zeta(11)

# S_{1,10} from Euler's formula
P6 = (zeta(2)*zeta(9) + zeta(3)*zeta(8) + zeta(4)*zeta(7) + zeta(5)*zeta(6))
S_1_10 = 6*zeta(11) - P6

# D_{S^(6)}(6) = -27/16 ζ(10) + 1/2 ζ(11) + 5/4 S_{1,10} - 1/4 S_{2,9}
D6 = mpf(-27)/16*zeta(10) + mpf(1)/2*zeta(11) + mpf(5)/4*S_1_10 - mpf(1)/4*S_2_9

print(f"\nD_{{S^(6)}}(6) = {nstr(D6, 60)}")

# PSLQ: reduced basis (cancellation hypothesis)
print("\nPSLQ {D6, ζ(10), ζ(11)}:")
rel6_red = pslq([D6, zeta(10), zeta(11)])
print(f"  Result: {rel6_red}")
if rel6_red:
    a, b, c = rel6_red
    print(f"  => D6 = ({-b}/{a})ζ(10) + ({-c}/{a})ζ(11)")
    res6 = a*D6 + b*zeta(10) + c*zeta(11)
    print(f"  Residual: {nstr(res6, 5)}")

# PSLQ: full basis (with products)
print("\nPSLQ {D6, ζ(10), ζ(11), π²ζ(9), π⁴ζ(7), π⁶ζ(5), π⁸ζ(3)}:")
basis6f = [D6, zeta(10), zeta(11), pi**2*zeta(9), pi**4*zeta(7),
           pi**6*zeta(5), pi**8*zeta(3)]
rel6_full = pslq(basis6f)
print(f"  Result: {rel6_full}")
if rel6_full:
    names6 = ['D6', 'ζ(10)', 'ζ(11)', 'π²ζ(9)', 'π⁴ζ(7)', 'π⁶ζ(5)', 'π⁸ζ(3)']
    for coeff, nm in zip(rel6_full, names6):
        if coeff != 0:
            print(f"    {coeff:+d} × {nm}")

print()

# ============================================================
# (Zα)^8 main computation (weight 15)
# ============================================================
print("="*60)
print("(Zα)^8 analysis (weight 15)")
print("="*60)

# Compute depth-2 MZVs
print("\nComputing depth-2 MZVs...")
t0 = time.time()
z132 = depth2_mzv(13, 2, N_MZV)
dt1 = time.time()-t0
print(f"  ζ(13,2) = {nstr(z132, 50)} ({dt1:.1f}s)")

t0 = time.time()
z123 = depth2_mzv(12, 3, N_MZV)
dt2 = time.time()-t0
print(f"  ζ(12,3) = {nstr(z123, 50)} ({dt2:.1f}s)")

t0 = time.time()
z114 = depth2_mzv(11, 4, N_MZV)
dt3 = time.time()-t0
print(f"  ζ(11,4) = {nstr(z114, 50)} ({dt3:.1f}s)")

# Euler sums
S_2_13 = z132 + zeta(15)
S_3_12 = z123 + zeta(15)
S_4_11 = z114 + zeta(15)

# S_{1,14} from Euler's formula
P8 = (zeta(2)*zeta(13) + zeta(3)*zeta(12) + zeta(4)*zeta(11) +
      zeta(5)*zeta(10) + zeta(6)*zeta(9) + zeta(7)*zeta(8))
S_1_14 = 8*zeta(15) - P8

print(f"\nEuler sums:")
print(f"  S_{{1,14}} = {nstr(S_1_14, 50)}")
print(f"  S_{{2,13}} = {nstr(S_2_13, 50)}")
print(f"  S_{{3,12}} = {nstr(S_3_12, 50)}")
print(f"  S_{{4,11}} = {nstr(S_4_11, 50)}")

# Assemble D_{S^(8)}(8) from decomposition:
# D = -205/64 ζ(14) + 5/8 ζ(15) + 15/4 S_{1,14} - 1/4 S_{2,13}
#     - 3/4 S_{3,12} - 1/4 S_{4,11}
D8 = (mpf(-205)/64*zeta(14) + mpf(5)/8*zeta(15) +
      mpf(15)/4*S_1_14 - mpf(1)/4*S_2_13 -
      mpf(3)/4*S_3_12 - mpf(1)/4*S_4_11)

print(f"\nD_{{S^(8)}}(8) = {nstr(D8, 80)}")

# Cross-check against previously computed value
D8_prev = mpf('-0.07833462075101528980838792233787228')
print(f"Previous:      {nstr(D8_prev, 35)}")
print(f"Agreement to:  {nstr(abs(D8 - D8_prev), 5)}")

# Direct computation cross-check (double sum, fewer digits)
print("\nDirect double-sum cross-check (N=500)...")
D8_direct = mpf(0)
for n in range(1, 501):
    Sn = mpf(0)
    for k in range(1, n):
        B4 = (35*mpf(k)**5 - 120*mpf(k)**4*n + 120*mpf(k)**3*n**2
              - 8*mpf(k)**2*n**3 - 24*mpf(k)*n**4 - 8*mpf(n)**5) / (128*mpf(k)**5*mpf(n)**8)
        Sn += 4*k*B4
    B4_diag = (35*mpf(n)**5 - 120*mpf(n)**5 + 120*mpf(n)**5
               - 8*mpf(n)**5 - 24*mpf(n)**5 - 8*mpf(n)**5) / (128*mpf(n)**13)
    Sn += 2*n*B4_diag
    D8_direct += Sn / mpf(n)**8
print(f"D8 (direct, N=500) = {nstr(D8_direct, 35)}")
print(f"D8 (Euler sums)    = {nstr(D8, 35)}")
print(f"Agreement to:      {nstr(abs(D8 - D8_direct), 5)}")

# ============================================================
# PSLQ analysis
# ============================================================
print("\n" + "="*60)
print("PSLQ Analysis for (Zα)^8")
print("="*60)

# Reduced basis: cancellation hypothesis
print("\n[1] Reduced basis {D8, ζ(14), ζ(15)}:")
rel8_red = pslq([D8, zeta(14), zeta(15)])
print(f"    Result: {rel8_red}")
if rel8_red:
    a, b, c = rel8_red
    print(f"    => D8 = ({-b}/{a})ζ(14) + ({-c}/{a})ζ(15)")
    res = a*D8 + b*zeta(14) + c*zeta(15)
    print(f"    Residual: {nstr(res, 5)}")
    print("    *** PRODUCTS CANCEL at (Zα)^8 ***")
else:
    print("    No relation found with reduced basis")

# Full basis: with products
print("\n[2] Full 9-element basis {D8, ζ(14), ζ(15), π²ζ(13), ..., π¹²ζ(3)}:")
basis8f = [D8, zeta(14), zeta(15),
           pi**2*zeta(13), pi**4*zeta(11), pi**6*zeta(9),
           pi**8*zeta(7), pi**10*zeta(5), pi**12*zeta(3)]
rel8_full = pslq(basis8f)
print(f"    Result: {rel8_full}")
if rel8_full:
    names8 = ['D8', 'ζ(14)', 'ζ(15)', 'π²ζ(13)', 'π⁴ζ(11)',
              'π⁶ζ(9)', 'π⁸ζ(7)', 'π¹⁰ζ(5)', 'π¹²ζ(3)']
    print("    Decomposition:")
    for coeff, nm in zip(rel8_full, names8):
        if coeff != 0:
            print(f"      {coeff:+d} × {nm}")
    res8 = sum(c*v for c,v in zip(rel8_full, basis8f))
    print(f"    Residual: {nstr(res8, 5)}")
    # Check if product coefficients are zero
    prod_coeffs = rel8_full[3:]
    if all(c == 0 for c in prod_coeffs):
        print("    *** ALL PRODUCT COEFFICIENTS ZERO => CANCELLATION CONFIRMED ***")
    else:
        print("    Product coefficients nonzero => cancellation FAILS at (Zα)^8")

# Summary
print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"(Zα)^4: complementarity D_S(s)+D_g(s)=(9/2)ζ(s-2) [exact, all s]")
print(f"(Zα)^6: Euler sum cancellation -> {'YES' if rel6_red else 'UNDETERMINED'}")
print(f"(Zα)^8: Euler sum cancellation -> {'YES' if rel8_red else 'NO/UNDETERMINED'}")
