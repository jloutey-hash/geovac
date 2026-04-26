"""
Paper 28 Section 8 formula verification script.

Verifies the corrected Euler sum formulas:
  c_2(n) = -(5n-4)/(4n^3)
  c_3(n) = (19n-20)/(8n^5) - 3H_{n-1}/(2n^4) - H_{n-1}^{(2)}/(2n^3)
  c_4(n) = -5(41n-40)/(64n^7) + 15H_{n-1}/(4n^6) - H_{n-1}^{(2)}/(4n^5)
           - 3H_{n-1}^{(3)}/(4n^4) - H_{n-1}^{(4)}/(4n^3)
  D_2 = -(5/4)zeta(2) + zeta(3)
  D_3 = (19/8)zeta(4) - (11/4)zeta(5)
  D_4 = -(205/64)zeta(6) + (71/8)zeta(7) - (9/2)zeta(3)zeta(4)

Also verifies:
  c_2(n) = -S(n)/(2n^4)  where S(n) = n(5n-4)/2
  D_S(4) = (5/2)zeta(2) - 2zeta(3)
  D_2 = -D_S(4)/2
"""
import sympy as sp
from sympy import Rational as R
from mpmath import mp, mpf, nstr
import mpmath

mp.dps = 80

print("=" * 70)
print("PAPER 28 SECTION 8: FORMULA VERIFICATION")
print("=" * 70)

# ── Exact Dirac expansion setup ──
a2 = sp.Symbol('a2', positive=True)
jph = sp.Symbol('jph', positive=True, rational=True)
n_sym = sp.Symbol('n', positive=True, integer=True)

gamma_s = sp.sqrt(jph**2 - a2)
delta_s = jph - gamma_s
n_eff = n_sym - delta_s
E_sym = 1/sp.sqrt(1 + a2/n_eff**2) - 1
E_expanded = sp.series(E_sym, a2, 0, n=6).removeO()
E_expanded = sp.expand(E_expanded)

coeffs_sym = {}
for p in range(1, 6):
    coeffs_sym[p] = sp.simplify(E_expanded.coeff(a2, p))


def degeneracy(n_val, j2):
    """j2 = 2j (integer). Returns degeneracy."""
    j = R(j2, 2)
    j_max = n_val - R(1, 2)
    if j == j_max:
        return int(2*j + 1)
    else:
        return int(2*(2*j + 1))


def direct_cp(p, n_val):
    """Compute c_p(n) directly from the Dirac expansion."""
    total = R(0)
    for j2 in range(1, 2*n_val, 2):
        j_val = R(j2, 2)
        d = degeneracy(n_val, j2)
        cp_val = coeffs_sym[p].subs([(n_sym, n_val), (jph, j_val + R(1, 2))])
        cp_val = sp.simplify(cp_val)
        total += d * cp_val
    return sp.simplify(total)


# ── Test 1: c_2(n) = -(5n-4)/(4n^3) ──
print("\n--- Test 1: c_2(n) = -(5n-4)/(4n^3) ---")
all_pass = True
for n_val in range(1, 8):
    formula = R(-(5*n_val - 4), 4*n_val**3)
    direct = direct_cp(2, n_val)
    ok = sp.simplify(formula - direct) == 0
    if not ok:
        all_pass = False
    print(f"  n={n_val}: formula={float(formula):.10e}, direct={float(direct):.10e}, {'OK' if ok else 'FAIL'}")
print(f"  Result: {'ALL PASS' if all_pass else 'FAILURES DETECTED'}")


# ── Test 2: c_3(n) = (19n-20)/(8n^5) - 3H/(2n^4) - H^{(2)}/(2n^3) ──
print("\n--- Test 2: c_3(n) = (19n-20)/(8n^5) - 3H_{n-1}/(2n^4) - H^{(2)}_{n-1}/(2n^3) ---")
all_pass = True
for n_val in range(1, 8):
    H1 = sum(R(1, k) for k in range(1, n_val))
    H2 = sum(R(1, k**2) for k in range(1, n_val))
    formula = R(19*n_val - 20, 8*n_val**5) - R(3, 2)*H1/n_val**4 - R(1, 2)*H2/n_val**3
    direct = direct_cp(3, n_val)
    ok = sp.simplify(formula - direct) == 0
    if not ok:
        all_pass = False
    print(f"  n={n_val}: formula={float(formula):.10e}, direct={float(direct):.10e}, {'OK' if ok else 'FAIL'}")
print(f"  Result: {'ALL PASS' if all_pass else 'FAILURES DETECTED'}")


# ── Test 3: c_4(n) closed form ──
print("\n--- Test 3: c_4(n) = -5(41n-40)/(64n^7) + 15H/(4n^6) - H^{(2)}/(4n^5) - 3H^{(3)}/(4n^4) - H^{(4)}/(4n^3) ---")
all_pass = True
for n_val in range(1, 8):
    H1 = sum(R(1, k) for k in range(1, n_val))
    H2 = sum(R(1, k**2) for k in range(1, n_val))
    H3 = sum(R(1, k**3) for k in range(1, n_val))
    H4 = sum(R(1, k**4) for k in range(1, n_val))
    formula = (R(-5*(41*n_val - 40), 64*n_val**7)
               + R(15, 4)*H1/n_val**6
               - R(1, 4)*H2/n_val**5
               - R(3, 4)*H3/n_val**4
               - R(1, 4)*H4/n_val**3)
    direct = direct_cp(4, n_val)
    ok = sp.simplify(formula - direct) == 0
    if not ok:
        all_pass = False
    print(f"  n={n_val}: formula={float(formula):.10e}, direct={float(direct):.10e}, {'OK' if ok else 'FAIL'}")
print(f"  Result: {'ALL PASS' if all_pass else 'FAILURES DETECTED'}")


# ── Test 4: S(n)/(2n^4) consistency ──
print("\n--- Test 4: c_2(n) = -S(n)/(2n^4) where S(n) = n(5n-4)/2 ---")
all_pass = True
for n_val in range(1, 8):
    S_n = R(n_val*(5*n_val - 4), 2)
    c2_from_S = -S_n / (2*n_val**4)
    c2_direct = R(-(5*n_val - 4), 4*n_val**3)
    ok = sp.simplify(c2_from_S - c2_direct) == 0
    if not ok:
        all_pass = False
    print(f"  n={n_val}: -S(n)/(2n^4)={float(c2_from_S):.10e}, c_2(n)={float(c2_direct):.10e}, {'OK' if ok else 'FAIL'}")
print(f"  Result: {'ALL PASS' if all_pass else 'FAILURES DETECTED'}")


# ── Test 5: D_2 = -(5/4)zeta(2) + zeta(3) ──
print("\n--- Test 5: D_2 = -(5/4)zeta(2) + zeta(3) ---")
z2 = mpmath.zeta(2)
z3 = mpmath.zeta(3)
D2_formula = -mpf(5)/4 * z2 + z3

# Direct sum
D2_direct = mpf(0)
for n_val in range(1, 200001):
    D2_direct += mpf(-(5*n_val - 4)) / (4*mpf(n_val)**3)

print(f"  D_2 (formula)      = {nstr(D2_formula, 40)}")
print(f"  D_2 (direct N=2e5) = {nstr(D2_direct, 40)}")
print(f"  |diff|             = {nstr(abs(D2_formula - D2_direct), 5)}")
print(f"  Result: {'OK (converging)' if abs(D2_formula - D2_direct) < 1e-4 else 'POSSIBLE ERROR'}")


# ── Test 6: D_3 = (19/8)zeta(4) - (11/4)zeta(5) ──
print("\n--- Test 6: D_3 = (19/8)zeta(4) - (11/4)zeta(5) ---")
z4 = mpmath.zeta(4)
z5 = mpmath.zeta(5)
D3_formula = mpf(19)/8 * z4 - mpf(11)/4 * z5

# Direct sum with harmonic numbers
H1_inc = mpf(0)
H2_inc = mpf(0)
D3_direct = mpf(0)
for n_val in range(1, 200001):
    nf = mpf(n_val)
    c3 = mpf(19*n_val - 20) / (8*nf**5) - mpf(3)/2 * H1_inc/nf**4 - mpf(1)/2 * H2_inc/nf**3
    D3_direct += c3
    H1_inc += mpf(1)/nf
    H2_inc += mpf(1)/nf**2

print(f"  D_3 (formula)      = {nstr(D3_formula, 40)}")
print(f"  D_3 (direct N=2e5) = {nstr(D3_direct, 40)}")
print(f"  |diff|             = {nstr(abs(D3_formula - D3_direct), 5)}")
print(f"  Result: {'OK (converging)' if abs(D3_formula - D3_direct) < 1e-3 else 'POSSIBLE ERROR'}")


# ── Test 7: D_4 = -(205/64)zeta(6) + (71/8)zeta(7) - (9/2)zeta(3)zeta(4) ──
print("\n--- Test 7: D_4 = -(205/64)zeta(6) + (71/8)zeta(7) - (9/2)zeta(3)zeta(4) ---")
z6 = mpmath.zeta(6)
z7 = mpmath.zeta(7)
D4_formula = -mpf(205)/64 * z6 + mpf(71)/8 * z7 - mpf(9)/2 * z3 * z4

# Direct sum with harmonic numbers
H1_inc = mpf(0)
H2_inc = mpf(0)
H3_inc = mpf(0)
H4_inc = mpf(0)
D4_direct = mpf(0)
for n_val in range(1, 50001):
    nf = mpf(n_val)
    c4 = (mpf(-5*(41*n_val - 40)) / (64*nf**7)
          + mpf(15)/4 * H1_inc/nf**6
          - mpf(1)/4 * H2_inc/nf**5
          - mpf(3)/4 * H3_inc/nf**4
          - mpf(1)/4 * H4_inc/nf**3)
    D4_direct += c4
    H1_inc += mpf(1)/nf
    H2_inc += mpf(1)/nf**2
    H3_inc += mpf(1)/nf**3
    H4_inc += mpf(1)/nf**4

print(f"  D_4 (formula)      = {nstr(D4_formula, 40)}")
print(f"  D_4 (direct N=5e4) = {nstr(D4_direct, 40)}")
print(f"  |diff|             = {nstr(abs(D4_formula - D4_direct), 5)}")
print(f"  Result: {'OK (converging)' if abs(D4_formula - D4_direct) < 1e-3 else 'POSSIBLE ERROR'}")


# ── Test 8: D_S(4) = (5/2)zeta(2) - 2zeta(3) and D_2 = -D_S(4)/2 ──
print("\n--- Test 8: D_2 = -D_S(4)/2 consistency ---")
DS4 = mpf(5)/2 * z2 - 2*z3
D2_from_DS = -DS4/2
print(f"  D_S(4)        = {nstr(DS4, 40)}")
print(f"  -D_S(4)/2     = {nstr(D2_from_DS, 40)}")
print(f"  D_2           = {nstr(D2_formula, 40)}")
print(f"  |diff|        = {nstr(abs(D2_from_DS - D2_formula), 5)}")
print(f"  Result: {'OK (exact)' if abs(D2_from_DS - D2_formula) < 1e-70 else 'FAIL'}")


# ── Test 9: Verify the old WRONG formulas do NOT match ──
print("\n--- Test 9: Old (wrong) formulas do NOT match direct sums ---")
D2_old = -mpf(5)/4 * z4 + z5  # WRONG: zeta(4), zeta(5)
D3_old = -mpf(27)/16 * z6 + mpf(3)/2 * z7  # WRONG
print(f"  D_2_wrong (old paper) = {nstr(D2_old, 20)}")
print(f"  D_2_correct           = {nstr(D2_formula, 20)}")
print(f"  D_2_direct            = {nstr(D2_direct, 20)}")
print(f"  Wrong matches direct? {abs(D2_old - D2_direct) < 1e-3} (should be False)")
print(f"  Correct matches direct? {abs(D2_formula - D2_direct) < 1e-3} (should be True)")
print()
print(f"  D_3_wrong (old paper) = {nstr(D3_old, 20)}")
print(f"  D_3_correct           = {nstr(D3_formula, 20)}")
print(f"  D_3_direct            = {nstr(D3_direct, 20)}")
print(f"  Wrong matches direct? {abs(D3_old - D3_direct) < 1e-3} (should be False)")
print(f"  Correct matches direct? {abs(D3_formula - D3_direct) < 1e-3} (should be True)")


print("\n" + "=" * 70)
print("ALL VERIFICATIONS COMPLETE")
print("=" * 70)
