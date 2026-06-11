"""
D6 Sommerfeld: analytical partial-fraction approach.

Strategy:
  1. Sympy Taylor expansion of Dirac formula → a_p(n, k)
  2. Partial fraction decomposition in k → exact beta_r coefficients
  3. Assemble D_p from Euler sums (Hurwitz, 250 dps)
  4. PSLQ decompose against MZV basis
"""
import sys, io, functools
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
print = functools.partial(print, flush=True)

import sympy as sp
from sympy import Rational as Rat, apart, cancel, expand, simplify, oo
from mpmath import mp, mpf, nstr, nsum, inf as mpinf
import mpmath
import time
import json
import os

mp.dps = 300

# ================================================================
# PART 1: SYMPY EXPANSION
# ================================================================
print("=" * 70)
print("PART 1: SYMPY DIRAC EXPANSION")
print("=" * 70)

a2 = sp.Symbol('a2', positive=True)
k = sp.Symbol('k', positive=True)  # k = j + 1/2
n = sp.Symbol('n', positive=True, integer=True)

gamma_s = sp.sqrt(k**2 - a2)
delta_s = k - gamma_s
n_eff = n - delta_s
E_sym = 1 / sp.sqrt(1 + a2 / n_eff**2) - 1

t0 = time.time()
E_expanded = sp.series(E_sym, a2, 0, n=7).removeO()
E_expanded = expand(E_expanded)
print(f"  Expansion took {time.time()-t0:.1f}s")

# Extract a_p(n, k) for p=2..6
coeffs_sym = {}
for p in range(1, 7):
    coeffs_sym[p] = cancel(E_expanded.coeff(a2, p))
    # Verify it's rational in k
    numer, denom = sp.fraction(coeffs_sym[p])
    print(f"  a_{p}(n,k): num deg(k)={sp.degree(numer, k)}, den deg(k)={sp.degree(denom, k)}")

# ================================================================
# PART 2: PARTIAL FRACTION DECOMPOSITION IN k
# ================================================================
print("\n" + "=" * 70)
print("PART 2: PARTIAL FRACTION IN k")
print("=" * 70)

def extract_harmonic_structure(p):
    """
    Decompose the interior sum contribution via partial fractions in k.

    c_p(n) = 2n * a_p(n,n) + sum_{k=1}^{n-1} 4k * a_p(n,k)

    For the interior sum, we need the partial fraction of 4k * a_p(n,k) in k.
    Each term c(n)/k^b contributes:
      b >= 2: c(n) * H_{n-1}^{(b-1)}  (harmonic number)
      b = 1:  c(n) * (n-1)             (polynomial)
      b = 0:  c(n) * n(n-1)/2          (polynomial)
      b < 0:  c(n) * power_sum         (polynomial)
    """
    ap = coeffs_sym[p]

    # Interior integrand: 4*k*a_p(n,k)
    integrand = 4 * k * ap
    integrand = cancel(integrand)

    # Partial fraction in k
    t0 = time.time()
    pf = apart(integrand, k)
    dt = time.time() - t0
    print(f"  p={p}: partial fraction took {dt:.1f}s")

    # Collect terms by power of k
    # pf is a sum of terms like c(n) * k^m or c(n) / k^b
    terms_by_power = {}  # power -> coefficient (function of n)

    # Use sp.Add.make_args to split the sum
    for term in sp.Add.make_args(pf):
        # Extract the power of k in this term
        # term = coeff_n * k^power
        coeff_n, power_k = sp.Wild('c'), sp.Wild('p')

        # Try to identify the k-dependence
        if term.is_Mul or term.is_Pow or term.has(k):
            # Collect terms as c(n) * k^m
            # Factor out the k-dependence
            c_of_n = term.subs(k, 1)
            ratio = cancel(term / c_of_n)

            # ratio should be k^m for some m
            if ratio == 1:
                power = 0
            elif ratio == k:
                power = 1
            elif ratio.is_Pow and ratio.base == k:
                power = ratio.exp
            elif (1/ratio).is_Pow and (1/ratio).base == k:
                power = -(1/ratio).exp
            elif ratio == 1/k:
                power = -1
            else:
                # More complex k-dependence; try to express as k^m
                # log(ratio) / log(k) should be constant
                try:
                    power = sp.simplify(sp.log(ratio) / sp.log(k))
                    if not power.is_Number:
                        print(f"    WARNING: non-power term: {term}")
                        continue
                except:
                    print(f"    WARNING: unrecognized term: {term}")
                    continue

            c_of_n = cancel(term * k**(-power))
            power = int(power)
        else:
            # Constant in k
            power = 0
            c_of_n = term

        if power not in terms_by_power:
            terms_by_power[power] = sp.Rational(0)
        terms_by_power[power] += c_of_n

    # Simplify coefficients
    for pw in terms_by_power:
        terms_by_power[pw] = cancel(terms_by_power[pw])

    return terms_by_power

# Actually, the above approach is fragile. Let me use a more robust method.
# Instead of parsing the partial fraction output, I'll extract coefficients
# by evaluating at specific k values (Laurent expansion around k=0).

def extract_harmonic_structure_v2(p):
    """
    Extract the k-Laurent coefficients of 4*k*a_p(n,k) around k=0.

    4*k*a_p(n,k) = sum_{m} c_m(n) * k^m

    The coefficient of k^{-b} (b >= 0) contributes to the sum over k=1..n-1:
      b >= 1 (i.e., 1/k^b in 4*k*a_p): sum k * (1/k^b) = sum 1/k^{b-1}

    Wait, I already have 4*k*a_p, so the term is c(n)/k^b in 4*k*a_p.
    The SUM from k=1 to n-1 of c(n)/k^b is:
      b >= 1: c(n) * H_{n-1}^{(b)}
      b = 0:  c(n) * (n-1)
      b < 0:  c(n) * sum_{k=1}^{n-1} k^{|b|} = c(n) * Faulhaber(|b|, n-1)
    """
    ap = coeffs_sym[p]
    integrand = 4 * k * ap
    integrand = cancel(integrand)

    # Laurent expansion around k = 0
    # First, find the order of the pole
    numer, denom = sp.fraction(integrand)
    denom_expanded = sp.Poly(denom, k)
    pole_order = 0
    while denom_expanded.eval(0) == 0:
        denom_expanded = sp.Poly(sp.quo(denom_expanded.as_expr(), k, k), k)
        pole_order += 1

    # Now expand as Laurent series
    # integrand = sum_{m=-pole_order}^{max_m} c_m(n) * k^m
    max_m = 5  # positive powers of k (contributes power sums)

    coefficients = {}  # m -> c_m(n)

    # Use series expansion
    series_expr = sp.series(integrand, k, 0, n=max_m + 1)

    for m in range(-pole_order - 2, max_m + 1):
        c_m = series_expr.coeff(k, m)
        if c_m != 0:
            c_m = cancel(c_m)
            coefficients[m] = c_m

    return coefficients, pole_order

# Even simpler: just use sympy's series expansion directly
def extract_structure_series(p):
    """
    Expand 4*k*a_p(n,k) as a Laurent series in k around k=0.
    """
    ap = coeffs_sym[p]
    integrand = 4 * k * ap
    integrand = cancel(integrand)

    t0 = time.time()

    # Determine pole order by checking denominator
    numer, denom = sp.fraction(integrand)

    # Factor out k from denominator
    pole_order = 0
    d = denom
    while d.subs(k, 0) == 0:
        d = cancel(d / k)
        pole_order += 1

    # Laurent series: multiply by k^pole_order, expand, divide back
    shifted = expand(integrand * k**pole_order)

    # Taylor expand the shifted expression
    max_terms = pole_order + 8  # enough to get harmonic + polynomial terms
    taylor = sp.series(shifted, k, 0, n=max_terms).removeO()

    coefficients = {}
    for m in range(max_terms):
        c_m = taylor.coeff(k, m)
        if c_m != 0:
            actual_power = m - pole_order  # actual power of k
            coefficients[actual_power] = cancel(c_m)

    dt = time.time() - t0
    print(f"  p={p}: series expansion took {dt:.1f}s, pole_order={pole_order}")

    return coefficients, pole_order

# Test with p=3 first (known answer: beta_1 = -3/2, beta_2 = -1/2)
print("\n--- Testing p=3 (known: beta_1=-3/2, beta_2=-1/2) ---")

coefficients_3, pole_3 = extract_structure_series(3)
print(f"  Laurent coefficients of 4k*a_3(n,k):")
for m in sorted(coefficients_3.keys()):
    print(f"    k^{m}: {coefficients_3[m]}")

# The coefficient of k^{-b} in 4k*a_p gives:
# sum_{k=1}^{n-1} c_{-b}(n)/k^b = c_{-b}(n) * H_{n-1}^{(b)} for b >= 1
# The "beta_r" in our original notation is: the coefficient of H^{(r)}/n^{2p-1-r}
# which comes from: c_{-r}(n) = beta_r / n^{2p-1-r}

# Check: for p=3, the H^{(1)} coefficient should give beta_1 = -3/2
# From the k^{-1} term: this coefficient times 1/k gives H^{(1)} sum
# Actually: 4k*a_p = ... + c_{-1}(n)/k + ... and sum 1/k = H_{n-1}^{(1)}
# Wait, the term in the integrand (4k*a_p) with k^{-1} gives sum_{k=1}^{n-1} of that term
# = c_{-1}(n) * sum_{k=1}^{n-1} 1/k = c_{-1}(n) * H_{n-1}^{(1)}

# For p=3: the H^{(1)} coefficient should be c_{-1}(n) and when we write
# c_3(n) = ... + c_{-1}(n) * H_{n-1}^{(1)} + ...
# Comparing with c_3(n) = ... + (-3/2)/n^4 * H_{n-1}^{(1)} + ...
# So c_{-1}(n) should be -3/2 / n^4 = -3/(2n^4)

# Similarly for H^{(2)}: from k^{-2} term:
# c_{-2}(n) * H_{n-1}^{(2)}
# Should give (-1/2)/n^3, so c_{-2}(n) = -1/(2n^3)

# ================================================================
# PART 3: FULL EXTRACTION FOR ALL p
# ================================================================
print("\n" + "=" * 70)
print("PART 3: EXTRACT BETAS AND RATIONAL PARTS")
print("=" * 70)

def assemble_cp_structure(p):
    """
    Fully decompose c_p(n) into harmonic and rational parts.

    Returns:
      betas: dict r -> beta_r (rational, function of n simplified)
      rational_terms: dict d -> a_d (coefficient of 1/n^d in rational part)
      endpoint_rational: rational function of n from 2n*a_p(n,n)
    """
    coeffs, pole_order = extract_structure_series(p)

    # Interior sum contributions
    # k^{-b} term in 4k*a_p → contributes c_{-b}(n) * H_{n-1}^{(b)} for b >= 1
    # k^0 term → c_0(n) * (n-1)
    # k^m term (m>0) → c_m(n) * power_sum(m, n-1)

    betas = {}  # r -> c_{-r}(n) coefficient of H_{n-1}^{(r)}
    poly_from_interior = sp.Rational(0)  # rational part from interior sum

    for m in sorted(coeffs.keys()):
        c_m = coeffs[m]
        if m < 0:
            # Harmonic contribution: H_{n-1}^{(-m)}
            r = -m
            betas[r] = c_m
        elif m == 0:
            # Constant in k: contributes c_0(n) * (n-1)
            poly_from_interior += c_m * (n - 1)
        else:
            # Positive power of k: contributes c_m(n) * sum_{k=1}^{n-1} k^m
            # Use Faulhaber formula or Bernoulli
            power_sum = sp.summation(sp.Symbol('kk')**m, (sp.Symbol('kk'), 1, n - 1))
            poly_from_interior += c_m * power_sum

    # Endpoint contribution: 2n * a_p(n, n)
    ap_nn = coeffs_sym[p].subs(k, n)
    endpoint = cancel(2 * n * ap_nn)

    # Total rational part = endpoint + poly_from_interior
    rational_total = cancel(endpoint + poly_from_interior)

    # Express as sum of a_d / n^d
    # First, expand as a rational function of n
    rational_total = cancel(rational_total)

    return betas, rational_total, endpoint

# Verify with p=3
print("\n--- p=3 ---")
betas_3, rat_3, endpt_3 = assemble_cp_structure(3)
print(f"  Betas:")
for r in sorted(betas_3):
    print(f"    beta_{r} = coeff of H^({r}): {betas_3[r]}")
print(f"  Rational part: {rat_3}")

# Verify: compute c_3(n) from this decomposition for specific n values
print("\n  Verification:")
for n_val in [1, 2, 3, 5, 10, 20]:
    # Computed from decomposition
    predicted = rat_3.subs(n, n_val)
    for r in betas_3:
        H_val = sum(sp.Rational(1, kk**r) for kk in range(1, n_val))
        predicted += betas_3[r].subs(n, n_val) * H_val

    # Exact from sympy
    exact = sp.Rational(0)
    for j2 in range(1, 2*n_val, 2):
        j_val = sp.Rational(j2, 2)
        k_val = j_val + sp.Rational(1, 2)
        if j2 == 2*n_val - 1:
            d_val = int(2*j_val + 1)
        else:
            d_val = int(2*(2*j_val + 1))
        ap_val = coeffs_sym[3].subs([(n, n_val), (k, k_val)])
        exact += d_val * simplify(ap_val)

    err = simplify(predicted - exact)
    print(f"    n={n_val}: err = {err}")

# Now do all p=2..6
print("\n" + "=" * 70)
print("FULL EXTRACTION p=2..6")
print("=" * 70)

all_betas = {}
all_rationals = {}

for p in range(2, 7):
    print(f"\n{'='*50}")
    print(f"p={p}")
    print(f"{'='*50}")

    betas_p, rat_p, endpt_p = assemble_cp_structure(p)

    print(f"  Harmonic betas (as functions of n):")
    for r in sorted(betas_p):
        print(f"    beta_{r}(n) = {betas_p[r]}")

    # Verify for n=1..25
    max_err = 0
    ok = True
    for n_val in range(1, 26):
        predicted = rat_p.subs(n, n_val)
        for r in betas_p:
            H_val = sum(sp.Rational(1, kk**r) for kk in range(1, n_val))
            predicted += betas_p[r].subs(n, n_val) * H_val

        # Exact
        exact = sp.Rational(0)
        for j2 in range(1, 2*n_val, 2):
            j_val = sp.Rational(j2, 2)
            k_val = j_val + sp.Rational(1, 2)
            if j2 == 2*n_val - 1:
                d_val = int(2*j_val + 1)
            else:
                d_val = int(2*(2*j_val + 1))
            ap_val = coeffs_sym[p].subs([(n, n_val), (k, k_val)])
            exact += d_val * simplify(ap_val)

        err = simplify(predicted - exact)
        if err != 0:
            ok = False
            if n_val <= 5 or n_val in [10, 20, 25]:
                print(f"    n={n_val}: MISMATCH err={err}")

    print(f"  Verified n=1..25: {'PASS' if ok else 'FAIL'}")

    all_betas[p] = betas_p
    all_rationals[p] = rat_p

# ================================================================
# PART 4: COMPUTE D_p FROM EULER SUMS
# ================================================================
print("\n" + "=" * 70)
print("PART 4: ASSEMBLE D_p FROM EULER SUMS (250 dps)")
print("=" * 70)

z = lambda s: mpmath.zeta(s)

def euler_sum_hurwitz(r, s):
    """S_{r,s} = sum_{k>=1} H_k^{(r)} / k^s via Hurwitz zeta."""
    old_dps = mp.dps
    mp.dps = 310
    def term(kk):
        return mpmath.hurwitz(s, kk) / mpf(kk)**r
    result = nsum(term, [1, mpinf], method='richardson')
    mp.dps = old_dps
    return result

def compute_Dp_from_structure(p, betas_p, rat_p):
    """
    Compute D_p = sum_{n>=1} c_p(n) analytically.

    c_p(n) = R(n) + sum_r beta_r(n) * H_{n-1}^{(r)}

    D_p = sum_{n>=1} R(n) + sum_r [sum_{n>=1} beta_r(n) * H_{n-1}^{(r)}]

    For the harmonic part:
    sum_{n>=1} f(n) * H_{n-1}^{(r)} = sum_{n>=1} f(n) * [H_n^{(r)} - 1/n^r]
    = sum_{n>=1} f(n)*H_n^{(r)} - sum_{n>=1} f(n)/n^r

    If f(n) = c/n^s, then:
    sum f(n)*H_n^{(r)} = c * S_{r,s}
    sum f(n)/n^r = c * zeta(r+s)

    So the harmonic contribution of beta_r(n) = c/n^s is:
    c * [S_{r,s} - zeta(r+s)]
    """
    w = 2*p - 1  # total weight

    # First, decompose beta_r(n) into powers of 1/n
    harmonic_contributions = {}  # (r, s) -> coefficient

    for r in betas_p:
        beta_expr = betas_p[r]
        # beta_r(n) should be c/n^s for some rational c and integer s
        # Express as a sum of terms a_j / n^j
        beta_expanded = sp.series(beta_expr, n, oo, n=10)

        for j in range(1, 20):
            coeff_j = beta_expanded.coeff(n, -j)
            if coeff_j != 0:
                harmonic_contributions[(r, j)] = sp.nsimplify(coeff_j)

    # Similarly, decompose R(n) into powers of 1/n
    rational_contributions = {}  # j -> coefficient of 1/n^j
    rat_expanded = sp.series(rat_p, n, oo, n=15)
    for j in range(1, 20):
        coeff_j = rat_expanded.coeff(n, -j)
        if coeff_j != 0:
            rational_contributions[j] = sp.nsimplify(coeff_j)

    print(f"\n  p={p}:")
    print(f"  Harmonic structure:")
    for (r, s), c in sorted(harmonic_contributions.items()):
        print(f"    ({c}) * H_{{n-1}}^{{({r})}} / n^{s}  -> ({c}) * [S_{{{r},{s}}} - z({r+s})]")

    print(f"  Rational structure:")
    for j, c in sorted(rational_contributions.items()):
        print(f"    ({c}) / n^{j}  -> ({c}) * z({j})")

    # Now compute D_p numerically at 300 dps
    D_val = mpf(0)

    # Rational sum: sum c_j * zeta(j)
    for j, c in rational_contributions.items():
        c_f = mpf(sp.Rational(c))
        if j == 1:
            print(f"  WARNING: 1/n term in rational part (divergent!)")
            continue
        D_val += c_f * z(j)

    # Harmonic sum: c * [S_{r,s} - zeta(r+s)]
    for (r, s), c in harmonic_contributions.items():
        c_f = mpf(sp.Rational(c))
        t0 = time.time()
        S_rs = euler_sum_hurwitz(r, s)
        dt = time.time() - t0
        contribution = c_f * (S_rs - z(r + s))
        print(f"    S_{{{r},{s}}} computed in {dt:.1f}s, contribution = {nstr(contribution, 30)}")
        D_val += contribution

    return D_val

# Known D_p values for cross-check
known = {
    2: lambda: mpf(-5)/4*z(2) + z(3),
    3: lambda: mpf(19)/8*z(4) - mpf(11)/4*z(5),
    4: lambda: mpf(-205)/64*z(6) + mpf(71)/8*z(7) - mpf(9)/2*z(3)*z(4),
    5: lambda: mpf(497)/128*z(8) - mpf(467)/16*z(9) + mpf(385)/32*z(3)*z(6) + mpf(75)/8*z(4)*z(5),
}

Dp_values = {}
for p in range(2, 7):
    print(f"\n{'='*50}")

    if p not in all_betas:
        print(f"  p={p}: no structure extracted, skipping")
        continue

    D_val = compute_Dp_from_structure(p, all_betas[p], all_rationals[p])
    Dp_values[p] = D_val
    print(f"\n  D_{p} = {nstr(D_val, 60)}")

    # Cross-check against known values
    if p in known:
        known_val = known[p]()
        diff = abs(D_val - known_val)
        if diff > 0:
            digits = -int(mpmath.log10(diff))
        else:
            digits = 300
        print(f"  Cross-check: {digits} matching digits {'OK' if digits >= 200 else 'WARN'}")

# ================================================================
# PART 5: PSLQ DECOMPOSITION
# ================================================================
print("\n" + "=" * 70)
print("PART 5: PSLQ DECOMPOSITION")
print("=" * 70)

pslq_results = {}

for p in range(2, 7):
    if p not in Dp_values:
        continue

    w = 2*p - 1
    D_val = Dp_values[p]

    # Build weight-w MZV basis
    basis_vals = [D_val]
    basis_names = [f"D_{p}"]

    # Single zetas
    for s in range(w, 1, -1):
        basis_vals.append(z(s))
        basis_names.append(f"z({s})")

    # Products z(a)*z(b) with a+b=w, 2 <= a <= b
    for a in range(2, w//2 + 1):
        b = w - a
        basis_vals.append(z(a) * z(b))
        basis_names.append(f"z({a})z({b})")

    basis_vals.append(mpf(1))
    basis_names.append("1")

    print(f"\n--- D_{p} (weight {w}, {len(basis_vals)} basis) ---")

    rel = None
    for mc in [1000, 10000, 100000, 1000000]:
        rel = mpmath.pslq(basis_vals, maxcoeff=mc)
        if rel is not None:
            break

    if rel is None:
        print(f"  PSLQ FAILED")
        pslq_results[p] = None
        continue

    d_coeff = rel[0]
    terms = []
    for c, nm in zip(rel[1:], basis_names[1:]):
        if c != 0:
            terms.append(f"({sp.Rational(-c, d_coeff)}){nm}")
    print(f"  D_{p} = {' + '.join(terms)}")

    res = sum(c*v for c, v in zip(rel, basis_vals))
    print(f"  Residual: {nstr(res, 5)}")

    # Check z(2)*z(2p-3)
    z2_target = f"z(2)z({w-2})"
    z2_coeff = 0
    for c, nm in zip(rel, basis_names):
        if nm == z2_target:
            z2_coeff = c

    if z2_coeff == 0:
        print(f"  *** z(2)*z({w-2}) = 0  CANCELLED ***")
    else:
        print(f"  *** z(2)*z({w-2}) PRESENT (coeff {z2_coeff}) ***")

    pslq_results[p] = {
        'relation': [int(x) for x in rel],
        'basis': basis_names,
        'z2_odd_cancelled': (z2_coeff == 0),
    }

# ================================================================
# PART 6: CANCELLATION TABLE
# ================================================================
print("\n" + "=" * 70)
print("PART 6: CANCELLATION TABLE (analytical)")
print("=" * 70)

def pslq_euler_sum(r, s, w):
    """PSLQ decompose S_{r,s} and extract z(2)*z(w-2) coefficient."""
    S_val = euler_sum_hurwitz(r, s)

    basis_vals = [S_val]
    basis_names = [f"S_{{{r},{s}}}"]

    for sv in range(w, 1, -1):
        basis_vals.append(z(sv))
        basis_names.append(f"z({sv})")

    for a in range(2, w//2 + 1):
        b = w - a
        basis_vals.append(z(a) * z(b))
        basis_names.append(f"z({a})z({b})")

    basis_vals.append(mpf(1))
    basis_names.append("1")

    rel = None
    for mc in [1000, 10000, 100000, 1000000]:
        rel = mpmath.pslq(basis_vals, maxcoeff=mc)
        if rel is not None:
            break

    if rel is None:
        return None, None

    z2_target = f"z(2)z({w-2})"
    mu = sp.Rational(0)
    if rel[0] != 0:
        for c, nm in zip(rel, basis_names):
            if nm == z2_target:
                mu = sp.Rational(-c, rel[0])
                break

    d_coeff = rel[0]
    terms = []
    for c, nm in zip(rel[1:], basis_names[1:]):
        if c != 0:
            terms.append(f"({sp.Rational(-c, d_coeff)}){nm}")
    decomp_str = ' + '.join(terms)

    return mu, decomp_str

for p in range(3, 7):
    if p not in all_betas:
        continue

    w = 2*p - 1
    print(f"\n=== p={p}, weight {w}, target: z(2)*z({w-2}) ===")

    # Get the simple beta_r values (coefficient in H^{(r)}/n^{2p-1-r} form)
    # From the harmonic structure, beta_r(n) = c_r / n^{s_r}
    # The "beta" in the cancellation table is just c_r, and the Euler sum is S_{r, s_r}

    # Extract (r, s, c) triples from the harmonic contributions
    beta_expr = all_betas[p]

    cancel_sum = sp.Rational(0)
    for r in sorted(beta_expr):
        # beta_r(n) as a function of n — expand to find c/n^s
        expr = beta_expr[r]
        series_r = sp.series(expr, n, oo, n=10)

        # Find the leading term c/n^s
        for s in range(1, 20):
            c_val = series_r.coeff(n, -s)
            if c_val != 0:
                c_val = sp.nsimplify(c_val)
                # This contributes c_val * [S_{r,s} - zeta(r+s)]
                # The z(2)*z(w-2) content comes from S_{r,s}

                print(f"\n  r={r}, s={s}: beta coeff = {c_val}")
                mu, decomp = pslq_euler_sum(r, s, w)
                if mu is None:
                    print(f"    S_{{{r},{s}}} PSLQ failed")
                    # Try to compute anyway
                    continue
                else:
                    print(f"    S_{{{r},{s}}} = {decomp}")
                    print(f"    mu_{r} (z(2)*z({w-2}) coeff) = {mu}")
                    contrib = c_val * mu
                    cancel_sum += contrib
                    print(f"    contribution = {c_val} * {mu} = {contrib}")

        # Check if there are additional terms in the expansion
        # (multiple powers of 1/n for the same r)

    print(f"\n  TOTAL z(2)*z({w-2}) coefficient = {cancel_sum}")
    if cancel_sum == 0:
        print(f"  *** CANCELLATION CONFIRMED ***")
    else:
        print(f"  *** NO CANCELLATION: residual = {cancel_sum} ***")

# ================================================================
# PART 7: SAVE
# ================================================================
print("\n" + "=" * 70)
print("PART 7: SAVE + SUMMARY")
print("=" * 70)

results = {
    'description': 'D6 analytical: partial-fraction harmonic extraction + Euler sum assembly',
    'Dp_values': {f'D{p}': nstr(Dp_values[p], 60) for p in Dp_values},
    'pslq': {f'D{p}': pslq_results[p] for p in pslq_results if pslq_results[p]},
    'betas': {},
}
for p in all_betas:
    results['betas'][f'p{p}'] = {str(r): str(all_betas[p][r]) for r in all_betas[p]}

out_path = os.path.join(os.path.dirname(__file__), 'data', 'product_survival_d6.json')
with open(out_path, 'w') as f:
    json.dump(results, f, indent=2, default=str)
print(f"Saved to {out_path}")

print("\n=== SUMMARY ===")
for p in range(2, 7):
    if p in pslq_results and pslq_results[p]:
        rel = pslq_results[p]['relation']
        names = pslq_results[p]['basis']
        d_coeff = rel[0]
        terms = []
        for c, nm in zip(rel[1:], names[1:]):
            if c != 0:
                terms.append(f"({sp.Rational(-c, d_coeff)}){nm}")
        z2c = pslq_results[p]['z2_odd_cancelled']
        print(f"  D_{p} = {' + '.join(terms)}  [z(2)z(odd) cancelled: {z2c}]")
    elif p in Dp_values:
        print(f"  D_{p} = {nstr(Dp_values[p], 40)} (PSLQ failed)")
    else:
        print(f"  D_{p}: not computed")

print("\nDone.")
