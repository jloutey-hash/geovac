"""
D6 Sommerfeld computation + product survival cancellation proof.

Strategy:
  1. Sympy: expand Dirac formula to a2^7, get a_p(n, jph) for p=2..6
  2. Lambdify with mpmath for fast 250-dps evaluation
  3. Compute D_p = sum c_p(n) via mpmath.nsum (Richardson extrapolation)
  4. PSLQ decompose against MZV basis at each weight
  5. Euler sums via Hurwitz zeta for the cancellation table
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import sympy as sp
from sympy import Rational as Rat, cancel, expand, simplify
from mpmath import mp, mpf, nstr, nsum, inf as mpinf
import mpmath
import time
import json
import os

mp.dps = 300  # extra margin; PSLQ needs ~250

# ================================================================
# PART 1: SYMPY EXPANSION
# ================================================================
print("=" * 70)
print("PART 1: SYMPY DIRAC EXPANSION")
print("=" * 70)

a2 = sp.Symbol('a2', positive=True)
jph_s = sp.Symbol('jph', positive=True)
n_s = sp.Symbol('n', positive=True, integer=True)

gamma_s = sp.sqrt(jph_s**2 - a2)
delta_s = jph_s - gamma_s
n_eff = n_s - delta_s
E_sym = 1 / sp.sqrt(1 + a2 / n_eff**2) - 1

t0 = time.time()
E_expanded = sp.series(E_sym, a2, 0, n=7).removeO()
E_expanded = expand(E_expanded)
print(f"  Expansion took {time.time()-t0:.1f}s")

coeffs_sym = {}
for p in range(1, 7):
    coeffs_sym[p] = cancel(E_expanded.coeff(a2, p))

# Create lambdified evaluators
evaluators = {}
for p in range(2, 7):
    evaluators[p] = sp.lambdify((n_s, jph_s), coeffs_sym[p], modules=['mpmath'])
    print(f"  a_{p}(n, jph) lambdified")

# ================================================================
# PART 2: VERIFY LAMBDIFIED EVALUATORS
# ================================================================
print("\n" + "=" * 70)
print("PART 2: CROSS-CHECK LAMBDIFIED VS SYMPY EXACT")
print("=" * 70)

def cp_sympy_exact(p, n_val):
    total = sp.Rational(0)
    for j2 in range(1, 2*n_val, 2):
        j_val = sp.Rational(j2, 2)
        k_val = j_val + sp.Rational(1, 2)
        if j2 == 2*n_val - 1:
            d = int(2*j_val + 1)
        else:
            d = int(2*(2*j_val + 1))
        ap_val = coeffs_sym[p].subs([(n_s, n_val), (jph_s, k_val)])
        total += d * simplify(ap_val)
    return simplify(total)

def cp_mpmath(p, n_int):
    n = mpf(n_int)
    f = evaluators[p]
    total = mpf(0)
    for k_int in range(1, n_int + 1):
        k = mpf(k_int)
        if k_int < n_int:
            d = 4 * k
        else:
            d = 2 * n
        total += d * f(n, k)
    return total

for p in [2, 3, 6]:
    for n_val in [1, 3, 5]:
        exact = float(cp_sympy_exact(p, n_val))
        numerical = float(cp_mpmath(p, n_val))
        reldiff = abs(exact - numerical) / (abs(exact) + 1e-300)
        print(f"  c_{p}({n_val}): exact={exact:.15e}, mpmath={numerical:.15e}, reldiff={reldiff:.2e}")

# ================================================================
# PART 3: COMPUTE D_p VIA mpmath.nsum
# ================================================================
print("\n" + "=" * 70)
print("PART 3: COMPUTE D_p VIA mpmath.nsum (250+ dps)")
print("=" * 70)

Dp_values = {}
z = lambda s: mpmath.zeta(s)

for p in range(2, 7):
    t0 = time.time()
    f = evaluators[p]

    def make_term(p_val, f_val):
        def term(n_int_mpf):
            n_int = int(n_int_mpf)
            n = mpf(n_int)
            total = mpf(0)
            for k_int in range(1, n_int + 1):
                k = mpf(k_int)
                if k_int < n_int:
                    d = 4 * k
                else:
                    d = 2 * n
                total += d * f_val(n, k)
            return total
        return term

    term_func = make_term(p, f)

    # Use nsum with enough working precision
    # For slowly converging sums (p=2), we need more terms
    try:
        D_val = nsum(term_func, [1, mpinf], method='richardson')
    except Exception as e:
        print(f"  D_{p}: nsum failed ({e}), trying direct sum + tail")
        # Fallback: direct sum to N + Euler-Maclaurin tail
        N = 2000
        D_val = mpf(0)
        for nn in range(1, N + 1):
            D_val += term_func(mpf(nn))
        # Rough tail correction: c_p(n) ~ C/n^{2p-2} for large n
        # Tail ~ C * [zeta(2p-2) - H_N^{(2p-2)}]
        # Estimate C from the last term
        c_last = term_func(mpf(N))
        C_est = c_last * mpf(N)**(2*p - 2)
        tail = C_est * (z(2*p - 2) - sum(mpf(1)/mpf(k)**(2*p-2) for k in range(1, N+1)))
        D_val += tail

    dt = time.time() - t0
    Dp_values[p] = D_val
    print(f"  D_{p} = {nstr(D_val, 60)}  ({dt:.1f}s)")

# ================================================================
# PART 4: CROSS-CHECK D_2..D_5 AGAINST KNOWN DECOMPOSITIONS
# ================================================================
print("\n" + "=" * 70)
print("PART 4: CROSS-CHECK AGAINST KNOWN D_2..D_5")
print("=" * 70)

known = {
    2: mpf(-5)/4*z(2) + z(3),
    3: mpf(19)/8*z(4) - mpf(11)/4*z(5),
    4: mpf(-205)/64*z(6) + mpf(71)/8*z(7) - mpf(9)/2*z(3)*z(4),
    5: mpf(497)/128*z(8) - mpf(467)/16*z(9) + mpf(385)/32*z(3)*z(6) + mpf(75)/8*z(4)*z(5),
}

for p in range(2, 6):
    diff = abs(Dp_values[p] - known[p])
    if diff > 0:
        digits = -int(mpmath.log10(diff))
    else:
        digits = 300
    status = "OK" if digits >= 200 else f"WARN ({digits} digits)"
    print(f"  D_{p}: |computed - known| = {nstr(diff, 5)} -> {digits} matching digits [{status}]")

# ================================================================
# PART 5: PSLQ DECOMPOSITION
# ================================================================
print("\n" + "=" * 70)
print("PART 5: PSLQ DECOMPOSITION")
print("=" * 70)

pslq_results = {}

for p in range(2, 7):
    w = 2*p - 1
    D_val = Dp_values[p]

    # Build weight-w MZV basis
    basis_vals = [D_val]
    basis_names = [f"D_{p}"]

    # Single zetas: only those of the correct parity contribute
    # At weight w=2p-1 (odd), the single zetas are zeta(w) and zeta(w-1)=zeta(2p-2)
    # But include all for safety
    for s in range(w, 1, -1):
        basis_vals.append(z(s))
        basis_names.append(f"z({s})")

    # Products of weight w: z(a)*z(b) with a+b=w, 2<=a<=b
    for a in range(2, w//2 + 1):
        b = w - a
        basis_vals.append(z(a) * z(b))
        basis_names.append(f"z({a})z({b})")

    basis_vals.append(mpf(1))
    basis_names.append("1")

    print(f"\n--- D_{p} (weight {w}, {len(basis_vals)} basis elements) ---")

    # Try PSLQ with increasing maxcoeff
    rel = None
    for mc in [1000, 10000, 100000, 1000000]:
        rel = mpmath.pslq(basis_vals, maxcoeff=mc)
        if rel is not None:
            break

    if rel is None:
        print(f"  PSLQ FAILED")
        pslq_results[p] = None
        continue

    # Display decomposition
    d_coeff = rel[0]
    print(f"  Relation: {rel}")
    print(f"  D_{p} = ", end="")
    terms_str = []
    for c, nm in zip(rel[1:], basis_names[1:]):
        if c != 0:
            coeff = sp.Rational(-c, d_coeff)
            terms_str.append(f"({coeff}){nm}")
    print(" + ".join(terms_str))

    res = sum(c*v for c, v in zip(rel, basis_vals))
    print(f"  Residual: {nstr(res, 5)}")

    # Check z(2)*z(2p-3) coefficient
    z2_target = f"z(2)z({w-2})"
    z2_coeff = 0
    for c, nm in zip(rel, basis_names):
        if nm == z2_target:
            z2_coeff = c
    if z2_coeff == 0:
        print(f"  *** z(2)*z({w-2}) = 0 -- CANCELLED ***")
    else:
        print(f"  *** z(2)*z({w-2}) PRESENT (coeff {z2_coeff}) ***")

    pslq_results[p] = {
        'relation': [int(x) for x in rel],
        'basis': basis_names,
        'residual': float(abs(res)),
        'z2_odd_cancelled': (z2_coeff == 0),
    }

# ================================================================
# PART 6: EULER SUMS + CANCELLATION TABLE
# ================================================================
print("\n" + "=" * 70)
print("PART 6: EULER SUM CANCELLATION TABLE")
print("=" * 70)

def euler_sum_hurwitz(r, s):
    """S_{r,s} = sum_{k>=1} hurwitz(s, k) / k^r"""
    old_dps = mp.dps
    mp.dps = 310
    def term(k):
        return mpmath.hurwitz(s, k) / mpf(k)**r
    result = nsum(term, [1, mpinf], method='richardson')
    mp.dps = old_dps
    return result

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
        return None, S_val

    # Extract z(2)*z(w-2) coefficient
    z2_target = f"z(2)z({w-2})"
    mu = sp.Rational(0)
    if rel[0] != 0:
        for c, nm in zip(rel, basis_names):
            if nm == z2_target:
                mu = sp.Rational(-c, rel[0])
                break

    # Display
    d_coeff = rel[0]
    terms = []
    for c, nm in zip(rel[1:], basis_names[1:]):
        if c != 0:
            terms.append(f"({sp.Rational(-c, d_coeff)}){nm}")
    print(f"    S_{{{r},{s}}} = {' + '.join(terms)}")
    res = sum(c*v for c, v in zip(rel, basis_vals))
    print(f"    Residual: {nstr(res, 5)}")

    return mu, S_val

# For the cancellation proof, I need the harmonic structure of c_p(n).
# Instead of symbolic extraction, I'll use the KNOWN structures for p=3..5
# and derive p=6 from the pattern.

# Known harmonic structures (from the literature on Sommerfeld expansion):
# c_p(n) = R(n)/n^{2p-1} + sum_{r=1}^{p-1} beta_r * H_{n-1}^{(r)} / n^{2p-1-r}

# For p=3: beta_1 = -3/2, beta_2 = -1/2
# Verified: D_3 = (19/8)z(4) - (11/4)z(5) with z(2)z(3) cancelled

# I need to extract the betas for p=4,5,6 from the sympy formula.
# Approach: compute c_p(n) exactly for many n, then solve for betas.

print("\nExtracting harmonic structures from sympy exact values...")

def extract_betas(p, n_check=25):
    """
    Extract beta_r coefficients from exact c_p(n) values.
    Ansatz: c_p(n) = sum_{d=0}^{D} a_d/n^{2p-1-d} + sum_{r=1}^{R} beta_r * H_{n-1}^{(r)} / n^{2p-1-r}
    """
    max_r = p - 1  # harmonic depth
    max_d = p      # polynomial degree in R(n)

    # Compute exact c_p(n)
    cp_exact = []
    for n_val in range(1, n_check + 1):
        cp_exact.append(cp_sympy_exact(p, n_val))

    # Build linear system
    # Unknowns: a_0, ..., a_{max_d} (rational part), beta_1, ..., beta_{max_r} (harmonic)
    n_unknowns = (max_d + 1) + max_r
    unknowns = sp.symbols(f'x0:{n_unknowns}')

    a_syms = unknowns[:max_d + 1]
    beta_syms = unknowns[max_d + 1:]

    equations = []
    for idx, n_val in enumerate(range(1, n_check + 1)):
        # Predicted value
        pred = sp.Rational(0)
        # Rational part: sum a_d / n^{2p-1-d}
        for d in range(max_d + 1):
            pred += a_syms[d] / sp.Integer(n_val)**(2*p - 1 - d)
        # Harmonic part: sum beta_r * H_{n-1}^{(r)} / n^{2p-1-r}
        for ir in range(max_r):
            r = ir + 1
            H_val = sum(sp.Rational(1, k**r) for k in range(1, n_val))
            pred += beta_syms[ir] * H_val / sp.Integer(n_val)**(2*p - 1 - r)

        equations.append(sp.Eq(pred, cp_exact[idx]))

    # Solve (overdetermined -> use first n_unknowns equations)
    sol = sp.solve(equations[:n_unknowns], unknowns)

    if not sol:
        print(f"  p={p}: FAILED to solve (trying more unknowns)")
        # Try with more polynomial terms
        max_d2 = max_d + 2
        n_unknowns2 = (max_d2 + 1) + max_r
        unknowns2 = sp.symbols(f'y0:{n_unknowns2}')
        a_syms2 = unknowns2[:max_d2 + 1]
        beta_syms2 = unknowns2[max_d2 + 1:]

        equations2 = []
        for idx, n_val in enumerate(range(1, min(n_check + 1, n_unknowns2 + 3))):
            pred = sp.Rational(0)
            for d in range(max_d2 + 1):
                pred += a_syms2[d] / sp.Integer(n_val)**(2*p - 1 - d)
            for ir in range(max_r):
                r = ir + 1
                H_val = sum(sp.Rational(1, k**r) for k in range(1, n_val))
                pred += beta_syms2[ir] * H_val / sp.Integer(n_val)**(2*p - 1 - r)
            equations2.append(sp.Eq(pred, cp_exact[idx]))

        sol = sp.solve(equations2[:n_unknowns2], unknowns2)
        if not sol:
            print(f"  p={p}: FAILED with extended ansatz too")
            return None
        a_syms = a_syms2
        beta_syms = beta_syms2
        max_d = max_d2

    # Extract betas
    betas = {}
    for ir in range(len(beta_syms)):
        if beta_syms[ir] in sol:
            betas[ir + 1] = sol[beta_syms[ir]]

    # Verify against all n_check values
    ok = True
    for idx, n_val in enumerate(range(1, n_check + 1)):
        pred = sp.Rational(0)
        for d in range(len(a_syms)):
            if a_syms[d] in sol:
                pred += sol[a_syms[d]] / sp.Integer(n_val)**(2*p - 1 - d)
        for ir in range(len(beta_syms)):
            r = ir + 1
            if beta_syms[ir] in sol:
                H_val = sum(sp.Rational(1, k**r) for k in range(1, n_val))
                pred += sol[beta_syms[ir]] * H_val / sp.Integer(n_val)**(2*p - 1 - r)

        err = sp.simplify(pred - cp_exact[idx])
        if err != 0:
            ok = False
            print(f"  p={p}, n={n_val}: MISMATCH err={err}")
            if idx > n_check - 3:
                break

    print(f"  p={p}: max_r={len(betas)}, verified={ok}")
    for r in sorted(betas):
        print(f"    beta_{r} = {betas[r]}")

    return betas

# Extract betas for each p
all_betas = {}
for p in range(3, 7):
    print(f"\n--- p={p} ---")
    betas = extract_betas(p, n_check=20)
    if betas is not None:
        all_betas[p] = betas

# Build cancellation table
print("\n" + "-" * 50)
print("CANCELLATION TABLE: sum beta_r * mu_r = 0 ?")
print("-" * 50)

for p in sorted(all_betas):
    betas = all_betas[p]
    w = 2*p - 1
    print(f"\n=== p={p}, weight {w}, target: z(2)*z({w-2}) ===")

    cancel_sum = sp.Rational(0)
    for r in sorted(betas):
        s = 2*p - 1 - r
        print(f"\n  r={r}, s={s}: beta_{r} = {betas[r]}")
        mu, S_val = pslq_euler_sum(r, s, w)
        if mu is None:
            print(f"    PSLQ failed for S_{{{r},{s}}}")
            continue
        print(f"    mu_{r} (coeff of z(2)z({w-2})) = {mu}")
        contrib = betas[r] * mu
        cancel_sum += contrib
        print(f"    beta_{r} * mu_{r} = {contrib}")

    print(f"\n  SUM = {cancel_sum}")
    if cancel_sum == 0:
        print(f"  *** CANCELLATION CONFIRMED ***")
    else:
        print(f"  *** RESIDUAL: {cancel_sum} ***")

# ================================================================
# PART 7: SAVE RESULTS
# ================================================================
print("\n" + "=" * 70)
print("PART 7: SAVE + SUMMARY")
print("=" * 70)

results = {
    'description': 'D6 product survival: Sommerfeld Dirichlet sums D2-D6 at 250 dps',
    'Dp_values': {},
    'pslq_results': {},
    'harmonic_betas': {},
    'product_survival': {},
}

for p in Dp_values:
    results['Dp_values'][f'D{p}'] = nstr(Dp_values[p], 60)

for p in pslq_results:
    if pslq_results[p]:
        results['pslq_results'][f'D{p}'] = pslq_results[p]

for p in all_betas:
    results['harmonic_betas'][f'p{p}'] = {str(r): str(v) for r, v in all_betas[p].items()}

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
        z2_cancelled = pslq_results[p]['z2_odd_cancelled']
        print(f"  D_{p} = {' + '.join(terms)}  [z(2)*z(odd) cancelled: {z2_cancelled}]")
    elif p in Dp_values:
        print(f"  D_{p} = {nstr(Dp_values[p], 30)} (PSLQ failed)")

print("\nDone.")
