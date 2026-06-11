"""
Fine structure (Zα)^8 Euler sum cancellation — v2 with analytical tail corrections.

Strategy: partial sums with cumulative H_n^{(r)} + Hurwitz zeta tail corrections.
The tail of Σ H_{n-1}^{(r)}/n^s is expressed exactly via ζ(s,a), ζ'(s,a), and
Bernoulli-number corrections, bypassing the flawed Euler-Maclaurin tail in v1.
"""
from mpmath import (mp, mpf, zeta, pi, pslq, nstr, euler as GAMMA,
                    bernoulli, factorial)
import time

mp.dps = 200
N_CUT = 50000
M_BERN = 8

def pochhammer(x, m):
    result = mpf(1)
    for j in range(m):
        result *= (x + j)
    return result

# ── Analytical tail corrections ──────────────────────────────────

def tail_rational(coeffs_powers, N_cut):
    """Σ_{n>N} Σ_i c_i / n^{p_i}  =  Σ_i c_i · ζ(p_i, N+1)"""
    a = N_cut + 1
    return sum(c * zeta(p, a) for c, p in coeffs_powers)

def tail_H1(s, N_cut, M=M_BERN):
    """
    Σ_{n>N} H_{n-1} / n^s

    H_{n-1} = γ + ψ(n), ψ(n) ~ ln n - 1/(2n) - Σ_{k≥1} B_{2k}/(2k · n^{2k})

    = γ·ζ(s,a) - ζ'(s,a) - ζ(s+1,a)/2 - Σ_{k=1}^M B_{2k}/(2k)·ζ(s+2k,a)
    """
    a = N_cut + 1
    result  = GAMMA * zeta(s, a)
    result -= zeta(s, a, derivative=1)        # -ζ'(s,a) = Σ ln(n+a)/(n+a)^s
    result -= zeta(s + 1, a) / 2
    for k in range(1, M + 1):
        result -= bernoulli(2*k) / (2*k) * zeta(s + 2*k, a)
    return result

def tail_Hr(r, s, N_cut, M=M_BERN):
    """
    Σ_{n>N} H_{n-1}^{(r)} / n^s   for r ≥ 2

    H_{n-1}^{(r)} = ζ(r) - ζ(r, n)

    Expansion of ζ(r, n):
      n^{1-r}/(r-1) + n^{-r}/2 + Σ_{k≥1} B_{2k}·(r)_{2k-1}/(2k)! · n^{1-r-2k}

    Result = ζ(r)·ζ(s,a) − 1/(r-1)·ζ(r+s−1,a) − ζ(r+s,a)/2
             − Σ_{k=1}^M B_{2k}·(r)_{2k-1}/(2k)! · ζ(r+s+2k−1,a)
    """
    a = N_cut + 1
    result  = zeta(r) * zeta(s, a)
    result -= zeta(r + s - 1, a) / (r - 1)
    result -= zeta(r + s, a) / 2
    for k in range(1, M + 1):
        coeff = bernoulli(2*k) * pochhammer(mpf(r), 2*k - 1) / factorial(2*k)
        result -= coeff * zeta(r + s + 2*k - 1, a)
    return result

# ── Partial sum + tail for D_{S^(2p)}(2p) ───────────────────────

def compute_D6(N_cut):
    """D_{S^(6)}(6) = Σ S^(6)(n)/n^6"""
    t0 = time.time()
    H1 = mpf(0)
    H2 = mpf(0)
    partial = mpf(0)
    for n in range(1, N_cut + 1):
        inv = mpf(1) / n
        partial += -3*(9*n - 8) / (16 * mpf(n)**11)
        partial += mpf(3)/2 * H1 * inv**10
        partial -= mpf(1)/2 * H2 * inv**9
        H1 += inv
        H2 += inv**2

    tail  = tail_rational([(-mpf(27)/16, 10), (mpf(3)/2, 11)], N_cut)
    tail += mpf(3)/2 * tail_H1(10, N_cut)
    tail -= mpf(1)/2 * tail_Hr(2, 9, N_cut)

    D = partial + tail
    print(f"  D6 computed in {time.time()-t0:.1f}s")
    return D

def compute_D8(N_cut):
    """D_{S^(8)}(8) = Σ S^(8)(n)/n^8"""
    t0 = time.time()
    H1 = mpf(0)
    H2 = mpf(0)
    H3 = mpf(0)
    H4 = mpf(0)
    partial = mpf(0)
    for n in range(1, N_cut + 1):
        inv = mpf(1) / n
        partial += -5*(41*n - 40) / (64 * mpf(n)**15)
        partial += mpf(15)/4 * H1 * inv**14
        partial -= mpf(1)/4  * H2 * inv**13
        partial -= mpf(3)/4  * H3 * inv**12
        partial -= mpf(1)/4  * H4 * inv**11
        H1 += inv
        H2 += inv**2
        H3 += inv**3
        H4 += inv**4

    tail  = tail_rational([(-mpf(205)/64, 14), (mpf(25)/8, 15)], N_cut)
    tail += mpf(15)/4 * tail_H1(14, N_cut)
    tail -= mpf(1)/4  * tail_Hr(2, 13, N_cut)
    tail -= mpf(3)/4  * tail_Hr(3, 12, N_cut)
    tail -= mpf(1)/4  * tail_Hr(4, 11, N_cut)

    D = partial + tail
    print(f"  D8 computed in {time.time()-t0:.1f}s")
    return D

# ── Main ─────────────────────────────────────────────────────────

print("="*60)
print("Sanity: ζ(2,1) = ζ(3)")
print("="*60)
# Direct check of tail_Hr(2, 1, ...) + partial
H2_partial = mpf(0)
partial_z21 = mpf(0)
for n in range(1, 1001):
    inv = mpf(1)/n
    partial_z21 += H2_partial * inv
    H2_partial += inv**2
z21_tail = tail_Hr(2, 1, 1000)
z21 = partial_z21 + z21_tail + zeta(3)  # S_{2,1} = ζ(2,1) + ζ(3), and ζ(2,1) = partial + tail
# Actually: Σ H_{n-1}^{(2)}/n = S_{2,1} = ζ(1,2) + ζ(3). And ζ(2,1) is the MZV.
# S_{2,1} = ζ(2)ζ(1)?? No, ζ(1) diverges.
# Let me just check: Σ_{n≥1} H_{n-1}^{(2)}/n^1 doesn't converge! H_{n-1}^{(2)} → ζ(2).
# Bad sanity check. Let me use S_{2,3} = Σ H_n^{(2)}/n^3 instead.
print("  Using S_{2,3} = Σ H_n^{(2)}/n^3")
H2p = mpf(0)
partial_s23 = mpf(0)
for n in range(1, 1001):
    inv = mpf(1)/n
    H2p += inv**2
    partial_s23 += H2p * inv**3
tail_s23 = tail_Hr(2, 3, 1000) + zeta(5)  # Σ H_{n-1}^{(2)}/n^3 + ζ(5) = S_{2,3}
# Wait: Σ H_n^{(2)}/n^3 = Σ [H_{n-1}^{(2)} + 1/n^2]/n^3 = Σ H_{n-1}^{(2)}/n^3 + ζ(5)
# partial_s23 computes Σ_{n=1}^{1000} H_n^{(2)}/n^3 (where H_n includes 1/n^2)
# So S_{2,3} ≈ partial_s23 + tail for Σ_{n>1000} H_n^{(2)}/n^3
# Σ_{n>N} H_n^{(2)}/n^3 = Σ_{n>N} H_{n-1}^{(2)}/n^3 + ζ(5, N+1) = tail_Hr(2,3,N) + ζ(5,N+1)
s23 = partial_s23 + tail_Hr(2, 3, 1000) + zeta(5, 1001)
# Known: S_{2,3} = 11/2 ζ(5) - 3 ζ(2)ζ(3)  [Euler, various sources]
# Actually the known value is: S_{2,3} = (11/2)ζ(5) - 3ζ(2)ζ(3)
# Hmm, let me verify: stuffle gives S_{2,3} + S_{3,2} = ζ(2)ζ(3) + ζ(5)
# And Euler's formula gives S_{2,3} = ?
# From the MZV literature: ζ(3,2) + ζ(2,3) = ζ(2)ζ(3) - ζ(5) [stuffle for MZVs]
# And S_{r,s} = ζ(s,r) + ζ(s+r), so S_{2,3} = ζ(3,2) + ζ(5)
# Also S_{3,2} = ζ(2,3) + ζ(5)
# S_{2,3} + S_{3,2} = ζ(3,2) + ζ(2,3) + 2ζ(5) = ζ(2)ζ(3) - ζ(5) + 2ζ(5) = ζ(2)ζ(3) + ζ(5) ✓
# Specific value: S_{2,3} = 3ζ(2)ζ(3) - 11/2 ζ(5) ... OR the other way?
# From Flajolet-Salvy: S_{2,3} = (11/2)ζ(5) - 3ζ(2)ζ(3)? That gives a negative number since ζ(2)ζ(3) ≈ 1.977.
# 11/2 * 1.0369 ≈ 5.703, 3*1.977 ≈ 5.932... so 5.703 - 5.932 = -0.229? But S_{2,3} should be positive!
# Let me compute: Σ H_n^{(2)}/n^3 for small n: n=1: 1/1 = 1; n=2: (1+1/4)/8 = 5/32 ≈ 0.156; sum ≈ 1.156 + ...
# So S_{2,3} > 1. The formula S_{2,3} = (11/2)ζ(5) - 3ζ(2)ζ(3) gives ≈ -0.229, which is WRONG.
# Let me use: S_{2,3} = 3ζ(2)ζ(3) - (11/2)ζ(5) ≈ 5.932 - 5.703 = 0.229...
# But we said S_{2,3} > 1. Hmm, 0.229 < 1. Let me recheck with n=1: H_1^{(2)} = 1, 1/1^3 = 1. So S_{2,3} ≥ 1. Contradiction!
# Wait: H_1^{(2)} = 1 and H_1^{(2)}/1^3 = 1. But then for n=2: H_2^{(2)} = 1 + 1/4 = 5/4, and 5/4/8 = 5/32 ≈ 0.156.
# For n=3: H_3^{(2)} = 1 + 1/4 + 1/9 = 49/36, and 49/(36·27) ≈ 0.050.
# So S_{2,3} ≈ 1 + 0.156 + 0.050 + ... ≈ 1.23. NOT 0.229.
# The formula S_{2,3} = 3ζ(2)ζ(3) - (11/2)ζ(5) ≈ 0.229 is clearly wrong.
# CORRECT formula must be: S_{2,3} = ζ(2)ζ(3) - ζ(5)/2  ??? No...
# Let me just skip the sanity check with a specific closed form and just compare at two N values.
print("  (Skipping closed-form sanity; will verify D6 against known result)")
print()

# ── (Zα)^6 verification ─────────────────────────────────────────

print("="*60)
print("(Zα)^6 verification")
print("="*60)

D6 = compute_D6(N_CUT)
print(f"\n  D6 = {nstr(D6, 100)}")

# Cross-check at N=25000
D6_check = compute_D6(25000)
print(f"  D6 (N=25000) = {nstr(D6_check, 100)}")
diff6 = abs(D6 - D6_check)
print(f"  |diff| = {nstr(diff6, 5)}")

# Known result from previous session (direct double-sum, 48 digits):
D6_known = mpf('-0.187687779614796351979684043203780765')
print(f"  D6 (known)   = {nstr(D6_known, 40)}")
print(f"  |D6 - known| = {nstr(abs(D6 - D6_known), 5)}")

print("\n  PSLQ {D6, ζ(10), ζ(15)}:")
rel = pslq([D6, zeta(10), zeta(11)], maxcoeff=10000)
print(f"    Result: {rel}")
if rel:
    a, b, c = rel
    print(f"    => D6 = ({-b}/{a})ζ(10) + ({-c}/{a})ζ(11)")
    res = a*D6 + b*zeta(10) + c*zeta(11)
    print(f"    Residual: {nstr(res, 5)}")

print("\n  PSLQ full {D6, ζ(10), ζ(11), π²ζ(9), π⁴ζ(7), π⁶ζ(5), π⁸ζ(3)}:")
basis6 = [D6, zeta(10), zeta(11), pi**2*zeta(9), pi**4*zeta(7),
          pi**6*zeta(5), pi**8*zeta(3)]
rel6f = pslq(basis6, maxcoeff=100000)
print(f"    Result: {rel6f}")
if rel6f:
    names6 = ['D6', 'ζ(10)', 'ζ(11)', 'π²ζ(9)', 'π⁴ζ(7)', 'π⁶ζ(5)', 'π⁸ζ(3)']
    for coeff, nm in zip(rel6f, names6):
        if coeff != 0:
            print(f"      {coeff:+d} × {nm}")
    prod_coeffs = rel6f[3:]
    if all(c == 0 for c in prod_coeffs):
        print("    *** PRODUCTS CANCEL ***")
    res6 = sum(c*v for c,v in zip(rel6f, basis6))
    print(f"    Residual: {nstr(res6, 5)}")

print()

# ── (Zα)^8 ──────────────────────────────────────────────────────

print("="*60)
print("(Zα)^8 analysis")
print("="*60)

D8 = compute_D8(N_CUT)
print(f"\n  D8 = {nstr(D8, 100)}")

D8_check = compute_D8(25000)
print(f"  D8 (N=25000) = {nstr(D8_check, 100)}")
diff8 = abs(D8 - D8_check)
print(f"  |diff| = {nstr(diff8, 5)}")

# Previous value
D8_prev = mpf('-0.07833462075101528980838792233787228')
print(f"  D8 (prev)    = {nstr(D8_prev, 40)}")
print(f"  |D8 - prev|  = {nstr(abs(D8 - D8_prev), 5)}")

print("\n  [1] PSLQ {D8, ζ(14), ζ(15)}:")
rel8r = pslq([D8, zeta(14), zeta(15)], maxcoeff=100000)
print(f"      Result: {rel8r}")
if rel8r:
    a, b, c = rel8r
    print(f"      => D8 = ({-b}/{a})ζ(14) + ({-c}/{a})ζ(15)")
    res = a*D8 + b*zeta(14) + c*zeta(15)
    print(f"      Residual: {nstr(res, 5)}")
    print("      *** PRODUCTS CANCEL at (Zα)^8 ***")
else:
    print("      No relation — cancellation FAILS")

print("\n  [2] PSLQ 9-element basis {D8, ζ(14), ζ(15), π²ζ(13), ..., π¹²ζ(3)}:")
basis8 = [D8, zeta(14), zeta(15),
          pi**2*zeta(13), pi**4*zeta(11), pi**6*zeta(9),
          pi**8*zeta(7), pi**10*zeta(5), pi**12*zeta(3)]
rel8f = pslq(basis8, maxcoeff=1000000)
print(f"      Result: {rel8f}")
if rel8f:
    names8 = ['D8', 'ζ(14)', 'ζ(15)', 'π²ζ(13)', 'π⁴ζ(11)',
              'π⁶ζ(9)', 'π⁸ζ(7)', 'π¹⁰ζ(5)', 'π¹²ζ(3)']
    print("      Decomposition:")
    for coeff, nm in zip(rel8f, names8):
        if coeff != 0:
            print(f"        {coeff:+d} × {nm}")
    res8 = sum(c*v for c,v in zip(rel8f, basis8))
    print(f"      Residual: {nstr(res8, 5)}")
    prod_coeffs = rel8f[3:]
    if all(c == 0 for c in prod_coeffs):
        print("      *** ALL PRODUCT COEFFICIENTS ZERO ***")
    else:
        print("      Product coefficients NONZERO")

# ── Summary ──────────────────────────────────────────────────────
print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"(Zα)^4: complementarity D_S(s)+D_g(s)=(9/2)ζ(s-2)  [exact, all s]")
print(f"(Zα)^6: Euler sum cancellation -> {'YES' if (rel and all(c==0 for c in rel[3:])) else ('PARTIAL' if rel6f else 'UNDETERMINED')}")
print(f"(Zα)^8: reduced PSLQ -> {'CANCEL' if rel8r else 'NO CANCEL'}")
print(f"(Zα)^8: full PSLQ   -> {'FOUND' if rel8f else 'NO RELATION'}")
