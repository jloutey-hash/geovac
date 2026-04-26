"""
Verify S^(6)(n) closed form against exact Dirac energy expansion.

Computes the (Zα)^6 correction coefficient for each (n,j) state from
the exact Dirac formula, weights by degeneracy, and compares with the
claimed closed form:
  S^(6)(n) = -3(9n-8)/(16n^5) + 3H_{n-1}/(2n^4) - H_{n-1}^{(2)}/(2n^3)
"""
from mpmath import mp, mpf, sqrt, nstr

mp.dps = 80

def dirac_energy(n, j, alpha):
    """Exact Dirac energy E_{nj}/mc^2 (without rest mass)."""
    x = alpha**2
    kappa_tilde = j + mpf(1)/2
    n_r = n - kappa_tilde
    gamma = sqrt(kappa_tilde**2 - x)
    denom = n_r + gamma
    return mpf(1)/sqrt(1 + x/denom**2) - 1

def degeneracy(n, j):
    """Degeneracy of (n,j) level in Dirac hydrogen."""
    j_max = n - mpf(1)/2
    if abs(j - j_max) < 1e-10:
        return int(2*j + 1)
    else:
        return int(2*(2*j + 1))

def extract_coefficient(n, j, order, alpha_vals=None):
    """
    Extract the coefficient of (Zα)^{2*order} in E_{nj}/mc^2
    using Richardson extrapolation.
    """
    if alpha_vals is None:
        alpha_vals = [mpf(10)**(-k) for k in range(8, 15)]

    vals = []
    for a in alpha_vals:
        E = dirac_energy(n, j, a)
        # Subtract known lower orders
        x = a**2
        # Leading: -x/(2n^2)
        E -= -x / (2*mpf(n)**2)
        if order >= 2:
            # (Zα)^4 coefficient: -(x^2)/(2n^4) * (n/(j+1/2) - 3/4)
            c1 = mpf(n)/(j + mpf(1)/2) - mpf(3)/4
            E -= -x**2 / (2*mpf(n)**4) * c1
        if order >= 3:
            # (Zα)^6 coefficient: extract
            # E at this point ~ c_3 * x^3 + higher
            vals.append(E / x**3)
        elif order == 2:
            vals.append(E / x**2)

    # Richardson: take the average of nearby extrapolation values
    # At very small alpha, the leading term dominates
    return vals[-1]  # smallest alpha, most accurate

def S6_formula(n):
    """Closed form: S^(6)(n) = -3(9n-8)/(16n^5) + 3H_{n-1}/(2n^4) - H_{n-1}^{(2)}/(2n^3)"""
    H1 = sum(mpf(1)/k for k in range(1, n))
    H2 = sum(mpf(1)/k**2 for k in range(1, n))
    return -3*(9*n - 8)/(16*mpf(n)**5) + mpf(3)/2 * H1/mpf(n)**4 - mpf(1)/2 * H2/mpf(n)**3

def S6_from_dirac(n):
    """Compute S^(6)(n) from exact Dirac formula."""
    total = mpf(0)
    j_vals = [mpf(k)/2 for k in range(1, 2*n, 2)]  # j = 1/2, 3/2, ..., n-1/2
    for j in j_vals:
        d = degeneracy(n, j)
        # The (Zα)^6 correction: E_{nj} ~ ... + c_3(n,j) * (Zα)^6 + ...
        # We need to extract c_3(n,j) and weight by d
        c3 = extract_coefficient(n, j, order=3)
        total += d * c3
    return total

print("Verifying S^(6)(n) closed form vs exact Dirac expansion")
print("="*60)

for n in range(1, 6):
    s6_exact = S6_from_dirac(n)
    s6_form = S6_formula(n)
    diff = abs(s6_exact - s6_form)
    print(f"\nn={n}:")
    print(f"  From Dirac:    {nstr(s6_exact, 30)}")
    print(f"  Closed form:   {nstr(s6_form, 30)}")
    print(f"  |diff|:        {nstr(diff, 5)}")

    # Also show the (Zα)^4 coefficient for cross-check
    total_s4 = mpf(0)
    j_vals = [mpf(k)/2 for k in range(1, 2*n, 2)]
    for j in j_vals:
        d = degeneracy(n, j)
        c1 = mpf(n)/(j + mpf(1)/2) - mpf(3)/4
        c1_term = -c1 / (2*mpf(n)**4)
        total_s4 += d * c1_term
    # S^(4)(n) in convention where E ~ S^(4)(n) * x^2
    print(f"  S^(4)(n) check: {nstr(total_s4, 20)}")
    s4_expected = mpf(n) * (5*n - 4) / (2 * 2*mpf(n)**4)  # from S(n) = n(5n-4)/2
    # Actually S(n)/(2n^4) where S(n) = n(5n-4)/2
    # Hmm, convention issue. Let me just print both.

print("\n\nNow checking D6 decomposition directly")
print("="*60)

# Compute D6 = Σ S^(6)(n)/n^6 two ways:
# Way 1: direct sum of S^(6)(n)/n^6 with the closed form
N = 500
D6_direct = mpf(0)
for n in range(1, N+1):
    D6_direct += S6_formula(n) / mpf(n)**6

# Way 2: from Euler sum decomposition
from mpmath import zeta, pi
# S_{1,10} from Euler's formula
P6 = (zeta(2)*zeta(9) + zeta(3)*zeta(8) + zeta(4)*zeta(7) + zeta(5)*zeta(6))
S_1_10 = 6*zeta(11) - P6

# S_{2,9} from depth-2 MZV (using mpmath zeta at high precision, no tail needed)
# ζ(s,r) = Σ_{n>k≥1} 1/(n^s k^r) = Σ_{k≥1} [ζ(s) - H_k^{(s)}] / k^r
z92_partial = mpf(0)
zs = zeta(9)
H = mpf(0)
for k in range(1, N+1):
    H += mpf(1)/mpf(k)**9
    z92_partial += (zs - H) / mpf(k)**2
# Tail is small since we sum to 500, but let's estimate
# The tail term is Σ_{k>N} [ζ(9) - H_k^{(9)}] / k^2
# ζ(9) - H_k^{(9)} ~ 1/(8k^8) for large k, so tail ~ Σ 1/(8k^10) ~ 1/(8) * ζ(10,501)
z92_tail = zeta(10, N+1) / 8  # leading term only
z92 = z92_partial + z92_tail
S_2_9 = z92 + zeta(11)

# D6 from v1 decomposition
D6_v1 = -mpf(27)/16 * zeta(10) + mpf(1)/2 * zeta(11) + mpf(5)/4 * S_1_10 - mpf(1)/4 * S_2_9

# D6 from v2 decomposition
D6_v2 = -mpf(27)/16 * zeta(10) + mpf(1)/2 * zeta(11) + mpf(3)/2 * S_1_10 - mpf(1)/2 * S_2_9

print(f"D6 (direct sum N={N}):     {nstr(D6_direct, 40)}")
print(f"D6 (v1 coeff 5/4, -1/4):  {nstr(D6_v1, 40)}")
print(f"D6 (v2 coeff 3/2, -1/2):  {nstr(D6_v2, 40)}")
print(f"|v1 - direct|: {nstr(abs(D6_v1 - D6_direct), 5)}")
print(f"|v2 - direct|: {nstr(abs(D6_v2 - D6_direct), 5)}")
