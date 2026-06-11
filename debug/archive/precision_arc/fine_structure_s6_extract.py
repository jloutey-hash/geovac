"""
Extract S^(2p)(n) coefficients from the EXACT Dirac formula by polynomial fitting.
Independent verification of the closed forms for (Zα)^4, (Zα)^6, (Zα)^8.
"""
from mpmath import mp, mpf, sqrt, nstr, matrix, lu_solve, zeta, pi, pslq
import sympy as sp

mp.dps = 120

def dirac_energy_exact(n, j_half_int, x):
    """
    E_{nj} / mc² from the exact Dirac formula.
    j_half_int: j as a half-integer (1 for j=1/2, 3 for j=3/2, etc.)
    x = (Zα)²
    Returns E/mc² - 1 (without rest mass).
    """
    j = mpf(j_half_int) / 2
    kappa = j + mpf(1)/2
    n_r = mpf(n) - kappa
    gamma = sqrt(kappa**2 - x)
    N = n_r + gamma
    return mpf(1) / sqrt(1 + x / N**2) - 1

def degeneracy(n, j_half_int):
    """Degeneracy of (n, j) level. j_half_int = 2j (so 1 for j=1/2, etc.)."""
    j = j_half_int / 2
    j_max = n - 0.5
    if abs(j - j_max) < 0.01:
        return int(2*j + 1)
    else:
        return int(2*(2*j + 1))

def total_energy_coeff(n, x):
    """
    Degeneracy-weighted sum of E_{nj}/mc² at given x = (Zα)².
    Returns Σ_j d(n,j) × E_{nj}/mc²
    """
    total = mpf(0)
    for j2 in range(1, 2*n, 2):  # j = 1/2, 3/2, ..., n-1/2
        d = degeneracy(n, j2)
        E = dirac_energy_exact(n, j2, x)
        total += d * E
    return total

def extract_coefficients(n, max_order=5):
    """
    Extract coefficients c_p such that Σ_j d(n,j) E_{nj}/mc² = Σ_p c_p x^p
    using polynomial fitting at multiple x values.
    """
    # Use several small x values
    num_pts = max_order + 4  # over-determined for stability
    x_vals = [mpf(10)**(-4 - k*2) for k in range(num_pts)]

    # Compute total energy at each x
    y_vals = [total_energy_coeff(n, x) for x in x_vals]

    # Solve Vandermonde system: Σ_p c_p x_i^p = y_i
    # For better conditioning, solve the divided-out system
    # c_1 + c_2 x + c_3 x² + ... = y/x
    # But more robust: just use the Vandermonde directly

    A = matrix(num_pts, max_order)
    b = matrix(num_pts, 1)
    for i in range(num_pts):
        for p in range(max_order):
            A[i, p] = x_vals[i]**(p+1)  # x^1, x^2, ..., x^max_order
        b[i] = y_vals[i]

    # Least squares: A^T A c = A^T b
    ATA = A.T * A
    ATb = A.T * b
    c = lu_solve(ATA, ATb)

    return [c[p] for p in range(max_order)]

print("="*70)
print("EXTRACTING FINE-STRUCTURE COEFFICIENTS FROM EXACT DIRAC FORMULA")
print("="*70)

# Known S^(4)(n) = n(5n-4)/2 from previous session
print("\n--- (Zα)^4 verification ---")
for n_val in range(1, 7):
    coeffs = extract_coefficients(n_val, max_order=5)
    # coeffs[0] = coeff of x^1 in Σ d E
    # coeffs[1] = coeff of x^2 in Σ d E  (this is the (Zα)^4 correction)
    # coeffs[2] = coeff of x^3 in Σ d E  (this is the (Zα)^6 correction)

    s4_extracted = coeffs[1]
    s4_formula = mpf(n_val) * (5*n_val - 4) / 2

    # The convention: E_total at order x² = -x/(2n²) × S^(4)(n)/(2n²) × x
    # Actually the x² coefficient of Σ_j d(n,j) E_{nj}/mc² should be:
    # -(1/(2n²)) × (Zα)^4 correction / x² = ???
    # Let me just extract and compare

    # Leading: coeff of x = Σ d × (-1/(2n²)) = -2n² × 1/(2n²) = -1
    print(f"\nn={n_val}:")
    print(f"  c_1 (leading):  {nstr(coeffs[0], 20)}  (expect -1)")
    print(f"  c_2 (α⁴ corr):  {nstr(coeffs[1], 20)}")
    print(f"  c_3 (α⁶ corr):  {nstr(coeffs[2], 20)}")
    print(f"  c_4 (α⁸ corr):  {nstr(coeffs[3], 20)}")

    # Verify c_2: it should be -1/(2n^4) × S^(4)(n) where S(n) convention varies
    # From the formula: E/mc² at x² = -x²/(2n⁴) × [n/(j+1/2) - 3/4]
    # Degeneracy-weighted sum: Σ d × [-1/(2n⁴) × (n/(j+1/2) - 3/4)]
    # = -1/(2n⁴) × Σ d × (n/(j+1/2) - 3/4)
    # = -1/(2n⁴) × S(n)  where S(n) = n(5n-4)/2

    s4_check = -s4_formula / (2*mpf(n_val)**4)
    print(f"  c_2 expected:    {nstr(s4_check, 20)}  (-S(n)/(2n⁴))")

print("\n\n--- S^(6)(n) comparison ---")
print("Convention: c_3 = Σ_j d(n,j) × [coeff of x³ in E_{nj}/mc²]")
print("Closed form claims: S^(6)(n) = -3(9n-8)/(16n^5) + 3H_{n-1}/(2n^4) - H_{n-1}^{(2)}/(2n^3)")
print("But with normalization: c_3 = -S^(6)(n)/(2n^6) ... or something else")
print()

for n_val in range(1, 7):
    coeffs = extract_coefficients(n_val, max_order=5)
    c3 = coeffs[2]

    # Formula for S^(6)(n)
    H1 = sum(mpf(1)/k for k in range(1, n_val))
    H2 = sum(mpf(1)/k**2 for k in range(1, n_val))
    s6_raw = -3*(9*n_val - 8)/(16*mpf(n_val)**5) + mpf(3)/2 * H1/mpf(n_val)**4 - mpf(1)/2 * H2/mpf(n_val)**3

    # What normalization relates c_3 to S^(6)(n)?
    # By analogy with c_2: c_2 = -S(n)/(2n⁴)
    # Try: c_3 = -S6(n)/(2n⁶)?
    c3_guess = -s6_raw / (2*mpf(n_val)**6)

    # Or: c_3 = S6(n)/n^6?
    c3_guess2 = s6_raw / mpf(n_val)**6

    # Or just divide c_3 by s6_raw to find the ratio
    if abs(s6_raw) > 1e-50:
        ratio = c3 / s6_raw * mpf(n_val)**6
    else:
        ratio = mpf(0)

    print(f"n={n_val}: c_3 = {nstr(c3, 25)}")
    print(f"       S6/n^6 = {nstr(s6_raw/mpf(n_val)**6, 25)}")
    print(f"       ratio (c3 × n^6 / S6_raw) = {nstr(ratio, 15)}")
    print()

print("\n--- Direct degeneracy-weighted (Zα)^6 coefficient ---")
# For sympy verification
a2 = sp.Symbol('a2', positive=True)
jph = sp.Symbol('jph', positive=True)
n_sym = sp.Symbol('n', positive=True, integer=True)

# delta = jph - sqrt(jph^2 - a2)
gamma_s = sp.sqrt(jph**2 - a2)
delta_s = jph - gamma_s
n_eff = n_sym - delta_s

# E/mc^2 = 1/sqrt(1 + a2/n_eff^2) - 1
# Expand to order a2^4
E_sym = 1/sp.sqrt(1 + a2/n_eff**2) - 1
E_expanded = sp.series(E_sym, a2, 0, n=5).removeO()

# Coefficient of a2^3
c3_sym = E_expanded.coeff(a2, 3)
c3_simplified = sp.simplify(c3_sym)
print(f"Symbolic c_3(n, jph) = {c3_simplified}")

# Now compute S^(6)(n) = Σ_j d(n,j) × c_3(n,j) for n=1..5
print("\nDegeneracy-weighted sum (symbolic):")
for n_val in range(1, 6):
    total = sp.Rational(0)
    for j2 in range(1, 2*n_val, 2):
        j_val = sp.Rational(j2, 2)
        d = degeneracy(n_val, j2)
        c3_val = c3_simplified.subs([(n_sym, n_val), (jph, j_val + sp.Rational(1,2))])
        c3_val = sp.simplify(c3_val)
        total += d * c3_val
    total = sp.simplify(total)
    print(f"  n={n_val}: Σ d·c_3 = {total} = {float(total):.15e}")

    # Compare with S6_raw / n^6
    H1_val = sum(sp.Rational(1, k) for k in range(1, n_val))
    H2_val = sum(sp.Rational(1, k**2) for k in range(1, n_val))
    s6_formula_val = sp.Rational(-3*(9*n_val-8), 16*n_val**5) + sp.Rational(3,2)*H1_val/n_val**4 - sp.Rational(1,2)*H2_val/n_val**3
    s6_over_n6 = s6_formula_val / n_val**6
    print(f"          S6/n^6 = {s6_over_n6} = {float(s6_over_n6):.15e}")
