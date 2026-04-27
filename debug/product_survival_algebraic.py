"""
Algebraic analysis of the product survival rule in the Sommerfeld expansion.

Goal: understand WHY a_{p,1} = 0 (the ζ(2)×ζ(odd) coefficient vanishes) at every
order p in the Dirichlet sums D_p = Σ_n c_p(n).

Approach:
1. Extract the general closed form of c_p(n) for p=2..6
   c_p(n) = R_p(n) + Σ_r α_{p,r} H_{n-1}^{(r)} / n^{s_{p,r}}
   where R_p(n) is a rational function and α_{p,r} are rational coefficients.

2. The Dirichlet sum D_p decomposes as:
   D_p = [rational zeta terms from R_p] + Σ_r α_{p,r} × [S_{r, s_{p,r}} - ζ correction]
   where S_{r,s} = Σ H_n^{(r)}/n^s are Euler sums with known evaluations.

3. The ζ(2)ζ(odd) content comes from specific Euler sum evaluations.
   Show that the linear combination always zeroes out the ζ(2)×ζ(odd) coefficient.

4. Look for a pattern in the coefficients α_{p,r} that makes the cancellation manifest.
"""
import sympy as sp
from sympy import Rational as R, symbols, sqrt, series, simplify, cancel, expand
from fractions import Fraction
from mpmath import mp, mpf, nstr, zeta as mpzeta
import json
import os

mp.dps = 100

# ================================================================
# PART 1: Extract c_p(n) closed forms for p=2..6
# ================================================================

print("=" * 72)
print("PART 1: SYSTEMATIC EXTRACTION OF c_p(n) CLOSED FORMS")
print("=" * 72)

a2 = sp.Symbol('a2', positive=True)
jph = sp.Symbol('jph', positive=True)
n_sym = sp.Symbol('n', positive=True, integer=True)

gamma_s = sp.sqrt(jph**2 - a2)
delta_s = jph - gamma_s
n_eff = n_sym - delta_s
E_sym = 1 / sp.sqrt(1 + a2 / n_eff**2) - 1

# Expand to (Zα)^12 = a2^6
E_expanded = sp.series(E_sym, a2, 0, n=7).removeO()

coeffs_sym = {}
for p in range(1, 7):
    coeffs_sym[p] = sp.cancel(E_expanded.coeff(a2, p))


def degeneracy(n_val, j2):
    j = R(j2, 2)
    j_max = n_val - R(1, 2)
    if j == j_max:
        return int(2 * j + 1)
    else:
        return int(2 * (2 * j + 1))


# For each order p, compute c_p(n) for n=1..12 and fit the harmonic number structure
# c_p(n) = rational_part(n) + Σ_r coeff_r(n) × H_{n-1}^{(r)}
#
# The structure is: c_p(n) = R(n)/n^{2p+1} + Σ_{r=1}^{2p-3} β_{p,r} × H_{n-1}^{(r)} / n^{2p-1-r}
# where R(n) is a polynomial of degree 1 in n and β_{p,r} are rational constants.

# First, compute c_p(n) directly for enough values to determine the structure
print()
for p in range(2, 7):
    print(f"\n--- Order p={p} [(Zα)^{2*p}] ---")
    cp_values = {}
    for n_val in range(1, 13):
        total = R(0)
        for j2 in range(1, 2 * n_val, 2):
            j_val = R(j2, 2)
            d = degeneracy(n_val, j2)
            cp_val = coeffs_sym[p].subs([(n_sym, n_val), (jph, j_val + R(1, 2))])
            cp_val = sp.cancel(cp_val)
            total += d * cp_val
        total = sp.cancel(total)
        cp_values[n_val] = total

    # Now determine the harmonic number structure.
    # At n=1, H_{0}^{(r)} = 0 for all r, so c_p(1) = rational_part(1).
    # The rational part is a polynomial in 1/n with leading term ~1/n^{2p-1}.

    print(f"  c_{p}(1) = {cp_values[1]}  [= rational part at n=1]")

    # Compute harmonic numbers
    def H(n_val, r):
        return sum(R(1, k**r) for k in range(1, n_val))

    # For n >= 2, subtract the rational part (determined by c_p(1) structure)
    # and look for harmonic number coefficients.
    #
    # Strategy: the H^{(r)} coefficients can be extracted by noting that
    # c_p(n) - c_p(1)/n^{2p-1} should be expressible in terms of H's.
    # But this is tricky. Let me use a different approach:
    # solve a linear system for the coefficients.

    # The general form for order p has harmonic number terms H^{(r)} for r=1..2p-3
    # (based on the pattern from p=2,3,4,5) with the H^{(r)} term divided by n^{2p-1-r}.
    # Plus a rational part of the form (An+B)/(C*n^{2p-1}).
    max_r = 2 * p - 3  # maximum harmonic number order

    # Build the system: for each n_val >= 2, we have
    # c_p(n) = a/(n^{2p-1}) + b/(n^{2p-2}) + Σ_{r=1}^{max_r} β_r × H_{n-1}^{(r)} / n^{2p-1-r}
    # where a, b are the rational part coefficients.
    # Wait, the rational part might be more complex. Let me parametrize:
    # rational_part = Σ_{k=0}^{1} a_k / n^{2p-1-k}  (degree 1 polynomial in 1/n)
    # Actually from the known forms:
    # p=2: c_2(n) = (-5n+4)/(4n^3) = -5/(4n^2) + 1/n^3
    #   ... hmm, the rational part has terms in 1/n^2 and 1/n^3.
    # p=3: c_3(n) = (19n-20)/(8n^5) = 19/(8n^4) - 20/(8n^5) = 19/(8n^4) - 5/(2n^5)
    # p=4: from d5 computation, c_4(n) has rational part (An-B)/(Cn^7)
    # p=5: c_5(n) has rational part 7(71n-72)/(128n^9)
    #
    # Pattern: rational part = (αn - β)/(γ n^{2p-1})
    # where α, β, γ are specific to each p.

    # Use the exact values at small n to extract the harmonic number coefficients.
    # At n=1: all H = 0, so c_p(1) = rational_part(1)
    # At n=2: H_{1}^{(r)} = 1 for all r, so c_p(2) = rational_part(2) + Σ β_r / 2^{2p-1-r}
    # At n=3: H_{2}^{(r)} = 1 + 1/2^r

    # Let's set up a system. The unknowns are:
    # rational_coeff_a (coeff of 1/n^{2p-2}), rational_coeff_b (coeff of 1/n^{2p-1})
    # and β_r for r=1..max_r

    num_unknowns = 2 + max_r  # a, b, β_1, ..., β_{max_r}

    # We need num_unknowns equations, so use n=1..num_unknowns+1
    A_matrix = []
    b_vector = []

    for n_val in range(1, num_unknowns + 2):
        row = []
        # rational part: a/n^{2p-2} + b/n^{2p-1}
        row.append(R(1, n_val**(2*p - 2)))
        row.append(R(1, n_val**(2*p - 1)))
        # harmonic number terms: β_r × H_{n-1}^{(r)} / n^{2p-1-r}
        for r in range(1, max_r + 1):
            Hnr = H(n_val, r)
            row.append(Hnr / R(n_val**(2*p - 1 - r)))
        A_matrix.append(row)
        b_vector.append(cp_values[n_val])

    # Solve using sympy
    A_sp = sp.Matrix(A_matrix[:num_unknowns])
    b_sp = sp.Matrix(b_vector[:num_unknowns])
    solution = A_sp.solve(b_sp)

    # Verify with extra equations
    all_ok = True
    for i in range(num_unknowns, len(A_matrix)):
        pred = sum(A_matrix[i][j] * solution[j] for j in range(num_unknowns))
        actual = b_vector[i]
        diff = sp.cancel(pred - actual)
        if diff != 0:
            all_ok = False
            print(f"  VERIFICATION FAIL at n={i+1}: diff = {diff}")

    a_coeff = sp.cancel(solution[0])
    b_coeff = sp.cancel(solution[1])
    beta_coeffs = [sp.cancel(solution[2 + r]) for r in range(max_r)]

    print(f"  Rational part: ({a_coeff})/n^{2*p-2} + ({b_coeff})/n^{2*p-1}")
    print(f"  Harmonic number coefficients:")
    for r in range(max_r):
        print(f"    β_{r+1} (H^({r+1})/n^{2*p-2-r}): {beta_coeffs[r]}")
    print(f"  Verified against extra n-values: {all_ok}")

    # Store for later analysis
    if p == 2:
        # Known: c_2(n) = (-5n+4)/(4n^3) = -5/(4n^2) + 1/n^3
        # No harmonic number terms (max_r = 1, β_1 should be 0... wait)
        # Actually for p=2, max_r = 2*2-3 = 1. So there IS one harmonic number term.
        # Let me check...
        pass

# ================================================================
# PART 2: EXTRACT THE ζ(2)×ζ(odd) COEFFICIENT AT EACH ORDER
# ================================================================

print()
print("=" * 72)
print("PART 2: ζ(2)×ζ(odd) COEFFICIENT EXTRACTION")
print("=" * 72)
print()

# For each order p, the Dirichlet sum D_p = Σ c_p(n) decomposes as:
# D_p = [rational zetas from R_p] + Σ_r β_r × Σ_n H_{n-1}^{(r)} / n^{s_r}
#
# The harmonic number sums are:
# Σ_{n≥1} H_{n-1}^{(r)} / n^s = S_{r,s} - ζ(r+s)
#   where S_{r,s} = Σ_{n≥1} H_n^{(r)} / n^s  (standard Euler sum)
#
# The Euler sums S_{r,s} evaluate to combinations of zeta products.
# The ζ(2)×ζ(w-2) content at weight w = r+s comes from specific evaluations.
#
# KEY IDENTITY (Euler's formula):
# S_{1,s} = (s/2 + 1)ζ(s+1) - (1/2)Σ_{j=2}^{s-2} ζ(j)ζ(s+1-j)   [for s >= 3]
#
# The ζ(2)ζ(s-1) content of S_{1,s} is: -(1/2)ζ(2)ζ(s-1)  [the j=2 term]
#
# For S_{r,s} with r ≥ 2, the evaluations are more complex but known.

# Rather than implement all Euler sum evaluations symbolically, let me use
# a different approach: compute D_p both as a high-precision number AND
# as the theoretical expression, then use PSLQ to check that the ζ(2)×ζ(odd)
# coefficient is zero.

# But actually, the most illuminating approach for proving the cancellation
# algebraically is to track the ζ(2)×ζ(odd) coefficient through each step.

# At order p, the ζ(2)×ζ(2p-3) coefficient in D_p comes from:
# (i) The Euler sums Σ H_{n-1}^{(r)}/n^{2p-1-r} for each r
# (ii) Each such sum = S_{r, 2p-1-r} - ζ(2p-1)

# The ζ(2)ζ(2p-3) content of each Euler sum S_{r, s} at weight r+s = 2p-1:
# This requires the Euler sum evaluation for each (r, s) pair.

# Let me use a numerical approach: for each order p, compute the coefficient
# of ζ(2)×ζ(2p-3) in D_p by numerical extraction.

print("Strategy: compute the ζ(2)×ζ(2p-3) coefficient numerically via PSLQ")
print("at each order p=3..5, then look for the cancellation mechanism.")
print()

# For p=3: weight-5 products. ζ(2)ζ(3) coefficient.
# S_{1,4}(H_{n-1}) = Σ H_{n-1}^{(1)}/n^4 = S_{1,4} - ζ(5) = [3ζ(5) - ζ(2)ζ(3)] - ζ(5) = 2ζ(5) - ζ(2)ζ(3)
# S_{2,3}(H_{n-1}) = Σ H_{n-1}^{(2)}/n^3 = S_{2,3} - ζ(5) = [3ζ(2)ζ(3) - 9/2·ζ(5)] - ζ(5) = 3ζ(2)ζ(3) - 11/2·ζ(5)
#
# ζ(2)ζ(3) coefficient from S_{1,4}(H_{n-1}): -1
# ζ(2)ζ(3) coefficient from S_{2,3}(H_{n-1}): +3
#
# Total ζ(2)ζ(3) coefficient = β₁ × (-1) + β₂ × (+3)
# For p=3: β₁ = -3/2, β₂ = -1/2
# = (-3/2)(-1) + (-1/2)(3) = 3/2 - 3/2 = 0  ✓

# The pattern is: the ζ(2)ζ(odd) content of S_{r,s}(H_{n-1}) at weight w=r+s
# has a coefficient that depends on r. The β_r coefficients from the c_p(n)
# formula combine to give exactly zero.

# Let me extract the ζ(2)ζ(w-2) coefficient of S_{r,s}(H_{n-1}) for all relevant (r,s).
# This requires knowledge of Euler sum evaluations.

# The cleanest approach: use mpmath to COMPUTE each Euler sum numerically
# and extract the ζ(2)ζ(odd) coefficient via PSLQ.

from mpmath import pslq

def euler_sum_numerical(r, s, N=200000):
    """Compute S_{r,s}(H_{n-1}) = Σ_{n=1}^∞ H_{n-1}^{(r)} / n^s numerically."""
    mp.dps = 100
    total = mpf(0)
    H = mpf(0)
    for n_val in range(1, N + 1):
        total += H / mpf(n_val)**s
        H += mpf(1) / mpf(n_val)**r
    return total


def extract_z2_zodd_coeff(r, s):
    """Extract the coefficient of ζ(2)×ζ(r+s-2) in S_{r,s}(H_{n-1})."""
    w = r + s  # total weight
    odd_z = w - 2  # the odd zeta in the product ζ(2)×ζ(odd)

    # Compute the Euler sum numerically
    S_val = euler_sum_numerical(r, s)

    # PSLQ basis: [S_{r,s}, ζ(2)ζ(odd), ζ(w), 1]
    # For weight w = r+s, the terms that appear are:
    # ζ(w) (single zeta), ζ(2)ζ(w-2), and possibly other products
    z2 = mpzeta(2)
    z_odd = mpzeta(odd_z)
    z_w = mpzeta(w)

    # For small weights, we can do exact PSLQ
    # At weight 5: basis = [S, ζ(2)ζ(3), ζ(5)]
    # At weight 7: basis = [S, ζ(2)ζ(5), ζ(3)ζ(4), ζ(7)]
    # At weight 9: basis = [S, ζ(2)ζ(7), ζ(3)ζ(6), ζ(4)ζ(5), ζ(9)]

    if w == 5:
        basis = [S_val, z2 * z_odd, z_w]
        rel = pslq(basis)
        if rel:
            z2_coeff = -Fraction(rel[1], rel[0])
            return z2_coeff, rel
    elif w == 7:
        basis = [S_val, z2 * mpzeta(5), mpzeta(3) * mpzeta(4), z_w]
        rel = pslq(basis)
        if rel:
            z2_coeff = -Fraction(rel[1], rel[0])
            return z2_coeff, rel
    elif w == 9:
        basis = [S_val, z2 * mpzeta(7), mpzeta(3) * mpzeta(6), mpzeta(4) * mpzeta(5), z_w]
        rel = pslq(basis)
        if rel:
            z2_coeff = -Fraction(rel[1], rel[0])
            return z2_coeff, rel
    elif w == 11:
        basis = [S_val, z2 * mpzeta(9), mpzeta(3) * mpzeta(8), mpzeta(4) * mpzeta(7),
                 mpzeta(5) * mpzeta(6), z_w]
        rel = pslq(basis)
        if rel:
            z2_coeff = -Fraction(rel[1], rel[0])
            return z2_coeff, rel

    return None, None


# ================================================================
# PART 3: Build the cancellation table
# ================================================================

print()
print("=" * 72)
print("PART 3: ζ(2)×ζ(odd) CANCELLATION TABLE")
print("=" * 72)

# For each order p, extract the β_r coefficients, compute the ζ(2)ζ(odd) content
# of each Euler sum, and show that the weighted sum is zero.

# Re-extract β coefficients for p=2..6 using sympy (from Part 1 above, but let me
# redo it more carefully in a clean loop)

all_betas = {}
all_rational_parts = {}

for p in range(2, 7):
    max_r = 2 * p - 3
    cp_values = {}
    for n_val in range(1, max(13, 2 + max_r + 3)):
        total = R(0)
        for j2 in range(1, 2 * n_val, 2):
            j_val = R(j2, 2)
            d = degeneracy(n_val, j2)
            cp_val = coeffs_sym[p].subs([(n_sym, n_val), (jph, j_val + R(1, 2))])
            cp_val = sp.cancel(cp_val)
            total += d * cp_val
        cp_values[n_val] = sp.cancel(total)

    def H_sympy(n_val, r):
        return sum(R(1, k**r) for k in range(1, n_val))

    num_unknowns = 2 + max_r
    A_matrix = []
    b_vector = []

    for n_val in range(1, num_unknowns + 2):
        row = []
        row.append(R(1, n_val**(2*p - 2)))
        row.append(R(1, n_val**(2*p - 1)))
        for r in range(1, max_r + 1):
            Hnr = H_sympy(n_val, r)
            row.append(Hnr / R(n_val**(2*p - 1 - r)))
        A_matrix.append(row)
        b_vector.append(cp_values[n_val])

    A_sp = sp.Matrix(A_matrix[:num_unknowns])
    b_sp = sp.Matrix(b_vector[:num_unknowns])
    solution = A_sp.solve(b_sp)

    betas = [sp.cancel(solution[2 + r]) for r in range(max_r)]
    all_betas[p] = betas
    all_rational_parts[p] = (sp.cancel(solution[0]), sp.cancel(solution[1]))

    print(f"\np={p}: β coefficients (H^(r)/n^(2p-1-r) for r=1..{max_r}):")
    for r in range(max_r):
        print(f"  β_{r+1} = {betas[r]}")

# ================================================================
# PART 4: Compute ζ(2)×ζ(odd) content of each Euler sum S_{r,s}
# ================================================================

print()
print("=" * 72)
print("PART 4: ζ(2)×ζ(odd) CONTENT OF EULER SUMS")
print("=" * 72)

# For each (r, s) pair that appears at order p, extract the ζ(2)×ζ(odd) coeff

z2_coeffs_by_weight = {}

for p in range(3, 7):
    max_r = 2 * p - 3
    w = 2 * p - 1  # weight of the odd-zeta terms
    print(f"\np={p}: weight {w} terms")

    z2_content = []
    for r in range(1, max_r + 1):
        s = 2 * p - 1 - r
        if s < 2:
            z2_content.append((r, s, None, "s too small"))
            continue

        z2_coeff, rel = extract_z2_zodd_coeff(r, s)
        z2_content.append((r, s, z2_coeff, rel))
        print(f"  S_{{{r},{s}}}(H_{{n-1}}): ζ(2)ζ({w-2}) coeff = {z2_coeff}   [PSLQ: {rel}]")

    # Now compute the weighted sum: Σ_r β_r × (ζ(2)ζ(odd) coeff of S_{r,s})
    betas = all_betas[p]
    total_z2_coeff = Fraction(0)
    print(f"\n  Weighted sum for ζ(2)ζ({w-2}) coefficient:")
    for r_idx in range(max_r):
        r = r_idx + 1
        s = 2 * p - 1 - r
        if s < 2:
            continue
        z2c = z2_content[r_idx][2]
        if z2c is None:
            continue
        beta_frac = Fraction(betas[r_idx].p, betas[r_idx].q)
        contribution = beta_frac * z2c
        total_z2_coeff += contribution
        print(f"    β_{r} × ζ₂-coeff(S_{{{r},{s}}}) = ({beta_frac}) × ({z2c}) = {contribution}")

    print(f"\n  TOTAL ζ(2)ζ({w-2}) coefficient in D_{p}: {total_z2_coeff}")
    if total_z2_coeff == 0:
        print(f"  *** CANCELLATION CONFIRMED at p={p} ***")
    else:
        print(f"  *** WARNING: NON-ZERO! ***")

    z2_coeffs_by_weight[p] = {
        'weight': w,
        'z2_zodd': f"ζ(2)ζ({w-2})",
        'euler_sum_contributions': [(r+1, z2_content[r][2]) for r in range(max_r) if z2_content[r][2] is not None],
        'beta_coefficients': [str(b) for b in betas],
        'total_z2_coeff': str(total_z2_coeff),
        'cancellation': total_z2_coeff == 0
    }

# ================================================================
# PART 5: LOOK FOR PATTERN IN THE CANCELLATION MECHANISM
# ================================================================

print()
print("=" * 72)
print("PART 5: CANCELLATION MECHANISM PATTERN ANALYSIS")
print("=" * 72)
print()

# For each p, the cancellation is:
# Σ_r β_{p,r} × z2_coeff(S_{r, 2p-1-r}) = 0
#
# Let's denote z2_coeff(S_{r,s}) = μ_r (the ζ(2)ζ(s-2) coefficient of S_{r,s})
# Then the cancellation is: Σ β_r × μ_r = 0
#
# This is a single linear constraint on the β_r vector. Why does it hold?
# The β_r are determined by the Dirac formula. The μ_r are determined by Euler sum theory.
# The claim is that these two independent structures conspire to satisfy this constraint.

# Let me tabulate β_r and μ_r side by side for each p:

for p in range(3, 7):
    max_r = 2 * p - 3
    betas = all_betas[p]
    print(f"\np={p}: β_r and μ_r (ζ(2)ζ({2*p-3}) coefficients)")
    print(f"  {'r':>3s}  {'β_r':>15s}  {'μ_r':>15s}  {'β_r × μ_r':>15s}")
    running = Fraction(0)
    for r_idx in range(max_r):
        r = r_idx + 1
        s = 2 * p - 1 - r
        if s < 2:
            print(f"  {r:3d}  {str(betas[r_idx]):>15s}  {'(s<2)':>15s}  {'---':>15s}")
            continue
        z2c = z2_coeffs_by_weight[p]['euler_sum_contributions']
        mu_r = None
        for rr, cc in z2c:
            if rr == r:
                mu_r = cc
                break
        if mu_r is None:
            print(f"  {r:3d}  {str(betas[r_idx]):>15s}  {'N/A':>15s}  {'---':>15s}")
            continue
        beta_frac = Fraction(betas[r_idx].p, betas[r_idx].q)
        prod = beta_frac * mu_r
        running += prod
        print(f"  {r:3d}  {str(beta_frac):>15s}  {str(mu_r):>15s}  {str(prod):>15s}")
    print(f"  {'SUM':>3s}  {'':>15s}  {'':>15s}  {str(running):>15s}")

# ================================================================
# SAVE RESULTS
# ================================================================

results = {
    'description': 'Product survival rule algebraic analysis',
    'beta_coefficients': {
        str(p): [str(b) for b in betas]
        for p, betas in all_betas.items()
    },
    'z2_coefficients': {
        str(p): data for p, data in z2_coeffs_by_weight.items()
    }
}

outdir = os.path.join(os.path.dirname(__file__), 'data')
os.makedirs(outdir, exist_ok=True)
outpath = os.path.join(outdir, 'product_survival_algebraic.json')
with open(outpath, 'w') as f:
    json.dump(results, f, indent=2, default=str)

print(f"\nResults saved to {outpath}")
