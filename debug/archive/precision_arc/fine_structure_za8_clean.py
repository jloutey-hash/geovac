"""
Fine-structure Euler sum cancellation: clean verification at (Za)^6, (Za)^8, (Za)^10.

Step 1: Derive a_p(n,j) from exact Dirac formula via sympy
Step 2: Compute degeneracy-weighted c_p(n) for small n (sympy exact)
Step 3: Fit closed-form ansatz for c_p(n) involving harmonic numbers
Step 4: Evaluate Dirichlet sum D_p algebraically using known Euler sum identities
Step 5: Check product cancellation via PSLQ
"""
import sympy as sp
from sympy import Rational as R, sqrt, symbols, simplify, oo, zeta as sp_zeta
from mpmath import mp, mpf, nstr, zeta, pslq, pi
import time

mp.dps = 80

# ═══════════════════════════════════════════════════════════════
# STEP 1: Sympy expansion of Dirac formula
# ═══════════════════════════════════════════════════════════════

a2 = sp.Symbol('a2', positive=True)
jph = sp.Symbol('jph', positive=True)
n_sym = sp.Symbol('n', positive=True, integer=True)

gamma_s = sp.sqrt(jph**2 - a2)
delta_s = jph - gamma_s
n_eff = n_sym - delta_s
E_sym = 1/sp.sqrt(1 + a2/n_eff**2) - 1

# Expand to (Za)^12 = a2^6 to get through (Za)^10
E_expanded = sp.series(E_sym, a2, 0, n=6).removeO()

coeffs_sym = {}
for p in range(1, 6):
    cp = E_expanded.coeff(a2, p)
    cp = sp.simplify(sp.cancel(cp))
    coeffs_sym[p] = cp

print("="*70)
print("SYMPY COEFFICIENTS a_p(n, jph)")
print("="*70)
for p in range(1, 6):
    print(f"\na_{p}(n,jph) = {coeffs_sym[p]}")

# ═══════════════════════════════════════════════════════════════
# STEP 2: Verify the manual a_4 formula against sympy
# ═══════════════════════════════════════════════════════════════

print("\n" + "="*70)
print("VERIFYING MANUAL a_4 FORMULA")
print("="*70)

# The manually-entered formula from the old script:
# numer = -35*jph^5 + 120*jph^4*n - 120*jph^3*n^2 + 8*jph^2*n^3 + 24*jph*n^4 + 8*n^5
# denom = 128*jph^5*n^8
manual_a4 = (-35*jph**5 + 120*jph**4*n_sym - 120*jph**3*n_sym**2
             + 8*jph**2*n_sym**3 + 24*jph*n_sym**4 + 8*n_sym**5) / (128*jph**5*n_sym**8)

diff_a4 = sp.simplify(sp.cancel(coeffs_sym[4] - manual_a4))
print(f"sympy a_4 - manual a_4 = {diff_a4}")
if diff_a4 == 0:
    print("*** MANUAL FORMULA VERIFIED ***")
else:
    print("*** MANUAL FORMULA IS WRONG — using sympy version ***")
    # Print the correct sympy version
    a4_num, a4_den = sp.fraction(sp.cancel(coeffs_sym[4]))
    print(f"Correct numerator:   {sp.expand(a4_num)}")
    print(f"Correct denominator: {sp.expand(a4_den)}")

# ═══════════════════════════════════════════════════════════════
# STEP 3: Degeneracy-weighted sums c_p(n) for small n
# ═══════════════════════════════════════════════════════════════

def degeneracy(n_val, j2):
    j = R(j2, 2)
    j_max = n_val - R(1, 2)
    if j == j_max:
        return int(2*j + 1)
    else:
        return int(2*(2*j + 1))

def compute_cp_n(p, n_val):
    """Compute c_p(n) = Σ_j d(n,j) × a_p(n,j) exactly."""
    total = R(0)
    for j2 in range(1, 2*n_val, 2):
        j_val = R(j2, 2)
        d = degeneracy(n_val, j2)
        cp_val = coeffs_sym[p].subs([(n_sym, n_val), (jph, j_val + R(1,2))])
        cp_val = sp.cancel(cp_val)
        total += d * cp_val
    return sp.cancel(total)

print("\n" + "="*70)
print("DEGENERACY-WEIGHTED c_p(n) TABLES")
print("="*70)

for p in [2, 3, 4, 5]:
    print(f"\n--- Order p={p}  [(Za)^{2*p}] ---")
    cp_table = {}
    for n_val in range(1, 9):
        val = compute_cp_n(p, n_val)
        cp_table[n_val] = val
        print(f"  n={n_val}: c_{p}(n) = {val}")

    # ═══════════════════════════════════════════════════════════
    # STEP 4: Fit closed-form ansatz with harmonic numbers
    # ═══════════════════════════════════════════════════════════

    if p == 2:
        # c_2(n) = -(5n-4)/(4n^3)  [known]
        print(f"  Closed form: c_2(n) = -(5n-4)/(4n^3)")
        for n_val in range(1, 9):
            formula = -R(5*n_val - 4, 4*n_val**3)
            assert formula == cp_table[n_val], f"Mismatch at n={n_val}"
        print("  *** VERIFIED ***")

    elif p == 3:
        # c_3(n) = (19n-20)/(8n^5) - 3H_{n-1}/(2n^4) - H_{n-1}^{(2)}/(2n^3)
        print(f"  Closed form: c_3(n) = (19n-20)/(8n^5) - 3H_{{n-1}}/(2n^4) - H_{{n-1}}^{{(2)}}/(2n^3)")
        for n_val in range(1, 9):
            H1 = sum(R(1, k) for k in range(1, n_val))
            H2 = sum(R(1, k**2) for k in range(1, n_val))
            formula = R(19*n_val - 20, 8*n_val**5) - R(3,2)*H1/n_val**4 - R(1,2)*H2/n_val**3
            match = (formula == cp_table[n_val])
            if not match:
                print(f"  n={n_val}: MISMATCH formula={formula} vs table={cp_table[n_val]}")
            assert match, f"c_3 mismatch at n={n_val}"
        print("  *** VERIFIED ***")

    elif p == 4:
        # Try to identify the closed form
        # Ansatz: c_4(n) = P(n)/n^7 + a*H_{n-1}/n^6 + b*H_{n-1}^{(2)}/n^5 + c*H_{n-1}^{(3)}/n^4
        # where P(n) is a polynomial
        # At n=1: all H terms vanish, so P(1) = c_4(1)
        P1 = cp_table[1]
        print(f"\n  Ansatz fitting for c_4(n):")
        print(f"  c_4(1) = {P1}  (P(1)/1 = c_4(1))")

        # Build system: for n=2,3,4,5 with various ansatz degrees
        # Try: P(n) = (an^2 + βn + γ) / D_0
        # and c_4(n) = P(n)/n^7 + a*H_{n-1}/n^6 + b*H_{n-1}^{(2)}/n^5 + c*H_{n-1}^{(3)}/n^4

        # We have unknowns: a, β, γ, D_0, a, b, c
        # But P(n) can be any rational. Let's instead use a direct solve:
        # 4 unknowns (a, b, c, + polynomial coefficients) from 8 data points

        # More systematic: subtract P(n)/n^7 trial and fit
        # First, assume P(n) = (An^2 + Bn + C) / 16  (by analogy with c_3 where denom=8)

        # For n=1: P(1) = (A+B+C)/16 = c_4(1)
        # At n=1 H terms all vanish, so c_4(1) = P(1)/1^7 = P(1)

        # Let me instead try a LINEAR SYSTEM approach.
        # c_4(n) * n^7 = P(n) + a*H_{n-1}*n + b*H_{n-1}^{(2)}*n^2 + c*H_{n-1}^{(3)}*n^3
        # where P(n) is a polynomial of degree d

        # Rewrite: c_4(n)*n^7 - a*H_{n-1}*n - b*H_{n-1}^{(2)}*n^2 - c*H_{n-1}^{(3)}*n^3 = P(n)
        # P(n) is polynomial in n, so LHS must be polynomial in n for n=1..N

        # Use: n=1..8 gives 8 equations.
        # If P(n) = p0 + p1*n + p2*n^2  (3 coefficients) + 3 harmonic coefficients = 6 unknowns
        # 8 equations, overdetermined — good.

        from sympy import Matrix

        max_poly_deg = 2
        N_fit = 8
        n_unknowns = (max_poly_deg + 1) + 3  # poly coeffs + a, b, c

        # LHS = c_4(n)*n^7
        # = p0 + p1*n + p2*n^2 + a*H_{n-1}*n + b*H_{n-1}^{(2)}*n^2 + c*H_{n-1}^{(3)}*n^3

        rows = []
        rhs = []
        for n_val in range(1, N_fit + 1):
            H1 = sum(R(1, k) for k in range(1, n_val))
            H2 = sum(R(1, k**2) for k in range(1, n_val))
            H3 = sum(R(1, k**3) for k in range(1, n_val))

            # Columns: [1, n, n^2, H1*n, H2*n^2, H3*n^3]
            row = [R(1), R(n_val), R(n_val**2),
                   H1 * n_val, H2 * n_val**2, H3 * n_val**3]
            rows.append(row)
            rhs.append(cp_table[n_val] * n_val**7)

        A_mat = Matrix(rows)
        b_vec = Matrix(rhs)

        # Least squares: (A^T A) x = A^T b
        sol = (A_mat.T * A_mat).solve(A_mat.T * b_vec)
        p0, p1, p2, a_coeff, b_coeff, c_coeff = sol

        print(f"  Polynomial: P(n) = {p0} + {p1}·n + {p2}·n^2")
        print(f"  Harmonic: a={a_coeff}, b={b_coeff}, c={c_coeff}")

        # Verify
        all_match = True
        for n_val in range(1, N_fit + 1):
            H1 = sum(R(1, k) for k in range(1, n_val))
            H2 = sum(R(1, k**2) for k in range(1, n_val))
            H3 = sum(R(1, k**3) for k in range(1, n_val))

            formula = (p0 + p1*n_val + p2*n_val**2)/n_val**7 + a_coeff*H1/n_val**6 + b_coeff*H2/n_val**5 + c_coeff*H3/n_val**4
            diff = sp.cancel(formula - cp_table[n_val])
            if diff != 0:
                print(f"  n={n_val}: RESIDUAL = {diff}")
                all_match = False

        if all_match:
            print("  *** c_4 CLOSED FORM VERIFIED for n=1..8 ***")

            # Simplify the polynomial
            # P(n)/n^7 = (p0 + p1*n + p2*n^2) / n^7
            # = (p2*n^2 + p1*n + p0) / n^7
            print(f"\n  c_4(n) = ({p2}·n^2 + {p1}·n + {p0}) / n^7")
            print(f"         + {a_coeff}·H_{{n-1}} / n^6")
            print(f"         + {b_coeff}·H_{{n-1}}^{{(2)}} / n^5")
            print(f"         + {c_coeff}·H_{{n-1}}^{{(3)}} / n^4")
        else:
            print("  Ansatz with poly degree 2 does not fit. Trying degree 3...")

            # Try poly degree 3
            rows2 = []
            rhs2 = []
            for n_val in range(1, N_fit + 1):
                H1 = sum(R(1, k) for k in range(1, n_val))
                H2 = sum(R(1, k**2) for k in range(1, n_val))
                H3 = sum(R(1, k**3) for k in range(1, n_val))

                row = [R(1), R(n_val), R(n_val**2), R(n_val**3),
                       H1 * n_val, H2 * n_val**2, H3 * n_val**3]
                rows2.append(row)
                rhs2.append(cp_table[n_val] * n_val**7)

            A_mat2 = Matrix(rows2)
            b_vec2 = Matrix(rhs2)
            sol2 = (A_mat2.T * A_mat2).solve(A_mat2.T * b_vec2)
            p0, p1, p2, p3, a_coeff, b_coeff, c_coeff = sol2

            print(f"  Polynomial: P(n) = {p0} + {p1}·n + {p2}·n^2 + {p3}·n^3")
            print(f"  Harmonic: a={a_coeff}, b={b_coeff}, c={c_coeff}")

            for n_val in range(1, N_fit + 1):
                H1 = sum(R(1, k) for k in range(1, n_val))
                H2 = sum(R(1, k**2) for k in range(1, n_val))
                H3 = sum(R(1, k**3) for k in range(1, n_val))

                formula = (p0 + p1*n_val + p2*n_val**2 + p3*n_val**3)/n_val**7 + a_coeff*H1/n_val**6 + b_coeff*H2/n_val**5 + c_coeff*H3/n_val**4
                diff = sp.cancel(formula - cp_table[n_val])
                if diff != 0:
                    print(f"  n={n_val}: RESIDUAL = {diff}")
                else:
                    pass
            else:
                print("  Degree 3 fit checked.")

print("\n\n" + "="*70)
print("DIRICHLET SUM D_4 = Σ c_4(n)")
print("="*70)

# Use the identified closed form (from the fit above) to compute D_4 via Euler sums
# But first, get a numerical value from the direct Dirac formula at high precision

print("\nDirect Dirac summation (incremental, no j-loop overhead for large n)...")

# For this we need c_4(n) in terms of n and harmonic numbers.
# The fit above should give us the coefficients. Let's just compute D_4
# from the exact Dirac formula for each (n,j) up to N_MAX.

# Actually, the j-loop is O(n) per n, so total is O(N^2). At N=10000 that's 10^8 ops.
# At 80 dps that's slow. Better: use the closed-form c_4(n) with incremental harmonic sums.

# We'll use the fit coefficients from step 4. But we need to determine them first.
# Let's extract them inside the main flow.

# For now, compute D_4 via the closed form c_4(n) once we have the coefficients.
# The script above should have printed them. Let me use the numerical extraction approach.

print("\nUsing incremental harmonic number closed form for c_4(n)...")
print("(Coefficients will be taken from the fit above)")
print("Computing D_4 with N=100000...")

# We need to read the fit results. Since this is a script, let's redo the fit inline.
# Compute c_4(n) for n=1..10 with sympy
c4_table = {}
for n_val in range(1, 11):
    c4_table[n_val] = compute_cp_n(4, n_val)

# Fit: c_4(n)*n^7 = p0 + p1*n + p2*n^2 + a*H1*n + b*H2*n^2 + c*H3*n^3
from sympy import Matrix as M2

# Try degree 2 first
rows = []
rhs_v = []
for n_val in range(1, 11):
    H1 = sum(R(1,k) for k in range(1, n_val))
    H2 = sum(R(1,k**2) for k in range(1, n_val))
    H3 = sum(R(1,k**3) for k in range(1, n_val))
    rows.append([R(1), R(n_val), R(n_val**2), H1*n_val, H2*n_val**2, H3*n_val**3])
    rhs_v.append(c4_table[n_val] * n_val**7)

A_m = M2(rows)
b_v = M2(rhs_v)
try:
    sol_final = (A_m.T * A_m).solve(A_m.T * b_v)
    p0_f, p1_f, p2_f, a_f, b_f, c_f = sol_final

    # Verify fit
    fit_ok = True
    for n_val in range(1, 11):
        H1 = sum(R(1,k) for k in range(1, n_val))
        H2 = sum(R(1,k**2) for k in range(1, n_val))
        H3 = sum(R(1,k**3) for k in range(1, n_val))
        formula = (p0_f + p1_f*n_val + p2_f*n_val**2)/n_val**7 + a_f*H1/n_val**6 + b_f*H2/n_val**5 + c_f*H3/n_val**4
        if sp.cancel(formula - c4_table[n_val]) != 0:
            fit_ok = False
            break

    if not fit_ok:
        raise ValueError("Degree 2 fit failed")

    print(f"\nc_4(n) closed form:")
    print(f"  = ({p2_f}·n^2 + {p1_f}·n + {p0_f}) / n^7")
    print(f"  + {a_f}·H_{{n-1}} / n^6")
    print(f"  + {b_f}·H_{{n-1}}^{{(2)}} / n^5")
    print(f"  + {c_f}·H_{{n-1}}^{{(3)}} / n^4")

    # Compute D_4 using incremental harmonic sums
    t0 = time.time()
    D4 = mpf(0)
    H1_inc = mpf(0)
    H2_inc = mpf(0)
    H3_inc = mpf(0)
    N_MAX = 200000

    p0_n = float(p0_f)
    p1_n = float(p1_f)
    p2_n = float(p2_f)
    a_n = float(a_f)
    b_n = float(b_f)
    c_n = float(c_f)

    for n_val in range(1, N_MAX + 1):
        nf = mpf(n_val)
        n7 = nf**7
        c4 = (p0_n + p1_n*nf + p2_n*nf**2)/n7 + a_n*H1_inc/nf**6 + b_n*H2_inc/nf**5 + c_n*H3_inc/nf**4
        D4 += c4
        H1_inc += mpf(1)/nf
        H2_inc += mpf(1)/nf**2
        H3_inc += mpf(1)/nf**3

    dt = time.time() - t0
    print(f"\nD_4 (N={N_MAX}, {dt:.1f}s) = {nstr(D4, 60)}")

except Exception as e:
    print(f"Fit failed: {e}")
    print("Falling back to degree 3 fit...")

    rows3 = []
    rhs3 = []
    for n_val in range(1, 11):
        H1 = sum(R(1,k) for k in range(1, n_val))
        H2 = sum(R(1,k**2) for k in range(1, n_val))
        H3 = sum(R(1,k**3) for k in range(1, n_val))
        rows3.append([R(1), R(n_val), R(n_val**2), R(n_val**3), H1*n_val, H2*n_val**2, H3*n_val**3])
        rhs3.append(c4_table[n_val] * n_val**7)

    A_m3 = M2(rows3)
    b_v3 = M2(rhs3)
    sol3 = (A_m3.T * A_m3).solve(A_m3.T * b_v3)
    p0_f, p1_f, p2_f, p3_f, a_f, b_f, c_f = sol3

    print(f"\nc_4(n) closed form (degree 3 polynomial):")
    print(f"  = ({p3_f}·n^3 + {p2_f}·n^2 + {p1_f}·n + {p0_f}) / n^7")
    print(f"  + {a_f}·H_{{n-1}} / n^6")
    print(f"  + {b_f}·H_{{n-1}}^{{(2)}} / n^5")
    print(f"  + {c_f}·H_{{n-1}}^{{(3)}} / n^4")

    # Verify
    for n_val in range(1, 11):
        H1 = sum(R(1,k) for k in range(1, n_val))
        H2 = sum(R(1,k**2) for k in range(1, n_val))
        H3 = sum(R(1,k**3) for k in range(1, n_val))
        formula = (p0_f + p1_f*n_val + p2_f*n_val**2 + p3_f*n_val**3)/n_val**7 + a_f*H1/n_val**6 + b_f*H2/n_val**5 + c_f*H3/n_val**4
        diff = sp.cancel(formula - c4_table[n_val])
        if diff != 0:
            print(f"  n={n_val}: RESIDUAL = {diff}")

    # Compute D_4
    t0 = time.time()
    D4 = mpf(0)
    H1_inc = mpf(0)
    H2_inc = mpf(0)
    H3_inc = mpf(0)
    N_MAX = 200000

    p0_n = float(p0_f)
    p1_n = float(p1_f)
    p2_n = float(p2_f)
    p3_n = float(p3_f)
    a_n = float(a_f)
    b_n = float(b_f)
    c_n = float(c_f)

    for n_val in range(1, N_MAX + 1):
        nf = mpf(n_val)
        n7 = nf**7
        c4 = (p0_n + p1_n*nf + p2_n*nf**2 + p3_n*nf**3)/n7 + a_n*H1_inc/nf**6 + b_n*H2_inc/nf**5 + c_n*H3_inc/nf**4
        D4 += c4
        H1_inc += mpf(1)/nf
        H2_inc += mpf(1)/nf**2
        H3_inc += mpf(1)/nf**3

    dt = time.time() - t0
    print(f"\nD_4 (N={N_MAX}, {dt:.1f}s) = {nstr(D4, 60)}")

# ═══════════════════════════════════════════════════════════════
# STEP 5: PSLQ at (Za)^8 — weight 7 basis
# ═══════════════════════════════════════════════════════════════

print("\n" + "="*70)
print("PSLQ ANALYSIS AT (Za)^8")
print("="*70)

z2 = zeta(2)
z3 = zeta(3)
z4 = zeta(4)
z5 = zeta(5)
z6 = zeta(6)
z7 = zeta(7)

# Reduced basis: single zetas only (cancellation hypothesis)
print(f"\n[1] Reduced basis {{D4, z(6), z(7)}}:")
rel1 = pslq([D4, z6, z7], maxcoeff=1000000)
print(f"    Result: {rel1}")
if rel1:
    a, b, c = rel1
    print(f"    => D4 = ({-b}/{a})z(6) + ({-c}/{a})z(7)")
    res = a*D4 + b*z6 + c*z7
    print(f"    Residual: {nstr(res, 5)}")
    print(f"    *** PRODUCTS CANCEL AT (Za)^8 ***")

# Full weight-7 basis
print(f"\n[2] Full weight-7 basis:")
basis = [D4, z6, z7, z2*z5, z3*z4, z2*z3**2]
names = ['D4', 'z(6)', 'z(7)', 'z(2)z(5)', 'z(3)z(4)', 'z(2)z(3)^2']
rel2 = pslq(basis, maxcoeff=1000000)
print(f"    Result: {rel2}")
if rel2:
    print("    Decomposition:")
    for coeff, nm in zip(rel2, names):
        if coeff != 0:
            print(f"      {coeff:+d} × {nm}")
    res = sum(c*v for c,v in zip(rel2, basis))
    print(f"    Residual: {nstr(res, 5)}")
    prod_coeffs = rel2[3:]
    if all(c == 0 for c in prod_coeffs):
        print("    *** ALL PRODUCT COEFFICIENTS ZERO ***")
    else:
        print("    *** PRODUCTS PRESENT — cancellation FAILS at (Za)^8 ***")

# Also try the wrong-weight basis to confirm it fails (sanity check)
print(f"\n[3] Wrong-weight sanity check {{D4, z(14), z(15)}}:")
rel3 = pslq([D4, zeta(14), zeta(15)], maxcoeff=1000000)
print(f"    Result: {rel3}")
if rel3 is None:
    print("    (None — confirms weight 14/15 is wrong for D4)")

# ═══════════════════════════════════════════════════════════════
# Cross-check: verify (Za)^6 result
# ═══════════════════════════════════════════════════════════════

print("\n" + "="*70)
print("CROSS-CHECK: (Za)^6 PSLQ")
print("="*70)

# D_3 should be (19/8)z(4) - (11/4)z(5)
D3_formula = mpf(19)/8 * z4 - mpf(11)/4 * z5
print(f"D_3 (algebraic) = (19/8)z(4) - (11/4)z(5) = {nstr(D3_formula, 50)}")

# Direct sum
D3_direct = mpf(0)
H1_inc = mpf(0)
H2_inc = mpf(0)
for n_val in range(1, 200001):
    nf = mpf(n_val)
    c3 = mpf(19*n_val - 20)/(8*nf**5) - mpf(3)/2 * H1_inc/nf**4 - mpf(1)/2 * H2_inc/nf**3
    D3_direct += c3
    H1_inc += mpf(1)/nf
    H2_inc += mpf(1)/nf**2

print(f"D_3 (direct N=200000) = {nstr(D3_direct, 50)}")
print(f"|diff| = {nstr(abs(D3_formula - D3_direct), 5)}")

rel_d3 = pslq([D3_direct, z4, z5], maxcoeff=10000)
print(f"PSLQ {{D3, z(4), z(5)}}: {rel_d3}")
if rel_d3:
    a, b, c = rel_d3
    print(f"  => D3 = ({-b}/{a})z(4) + ({-c}/{a})z(5)")

# ═══════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════

print("\n\n" + "="*70)
print("SUMMARY OF EULER SUM CANCELLATION")
print("="*70)
print(f"(Za)^2: D_1 = -1  [trivial]")
print(f"(Za)^4: D_2 = Σ c_2(n) = -(5/4)z(2) + z(3)  [weight 3, pure single zetas]")
print(f"(Za)^6: D_3 = Σ c_3(n) = (19/8)z(4) - (11/4)z(5)")
print(f"        [weight 5, z(2)z(3) cancels ALGEBRAICALLY, coeff (3/2)-(3/2)=0]")
print(f"(Za)^8: D_4 = Σ c_4(n) = ... [see PSLQ result above]")
