"""
D6 fix: corrected series truncation + reduced PSLQ basis.

The d6_analytical.py run had two bugs:
1. sp.series(expr, n, oo, n=10) only captures terms to 1/n^9,
   missing beta_1 = 315/(32*n^10) for p=6 and the 1/n^{10},1/n^{11} rational terms.
2. Full MZV PSLQ basis is ill-conditioned (even zetas are π-dependent).
   Fix: use reduced basis {z(2p-2), z(2p-1), weight-(2p-1) products only}.
"""
import sys, io, functools, time, json, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
print = functools.partial(print, flush=True)

import sympy as sp
from sympy import Rational as Rat, cancel, expand, simplify
from mpmath import mp, mpf, nstr, nsum, inf as mpinf
import mpmath

mp.dps = 300
z = lambda s: mpmath.zeta(s)

# ============================================================
# PART 1: Sympy expansion (fast, reused from d6_analytical.py)
# ============================================================
print("=" * 70)
print("PART 1: SYMPY EXPANSION + PARTIAL FRACTION (FIXED)")
print("=" * 70)

a2 = sp.Symbol('a2', positive=True)
k = sp.Symbol('k', positive=True)
n = sp.Symbol('n', positive=True, integer=True)

gamma_s = sp.sqrt(k**2 - a2)
delta_s = k - gamma_s
n_eff = n - delta_s
E_sym = 1 / sp.sqrt(1 + a2 / n_eff**2) - 1

t0 = time.time()
E_expanded = sp.series(E_sym, a2, 0, n=7).removeO()
E_expanded = expand(E_expanded)
print(f"  Expansion took {time.time()-t0:.1f}s")

coeffs_sym = {}
for p in range(1, 7):
    coeffs_sym[p] = cancel(E_expanded.coeff(a2, p))

# ============================================================
# PART 2: Extract harmonic structure with FIXED series order
# ============================================================
print("\n" + "=" * 70)
print("PART 2: EXTRACT BETAS (series order n=20, fixes D6)")
print("=" * 70)

def extract_structure_series(p, max_order=20):
    ap = coeffs_sym[p]
    integrand = 4 * k * ap
    integrand = cancel(integrand)

    # Determine pole order
    pole_order = 0
    d = sp.fraction(integrand)[1]
    while d.subs(k, 0) == 0:
        d = cancel(d / k)
        pole_order += 1

    # Laurent series
    shifted = expand(integrand * k**pole_order)
    max_terms = pole_order + max_order
    taylor = sp.series(shifted, k, 0, n=max_terms).removeO()

    coefficients = {}
    for m in range(max_terms):
        c_m = taylor.coeff(k, m)
        if c_m != 0:
            actual_power = m - pole_order
            coefficients[actual_power] = cancel(c_m)

    return coefficients, pole_order


def assemble_cp_structure(p, max_order=20):
    coeffs, pole_order = extract_structure_series(p, max_order)

    betas = {}
    poly_from_interior = sp.Rational(0)

    for m in sorted(coeffs.keys()):
        c_m = coeffs[m]
        if m < 0:
            r = -m
            betas[r] = c_m
        elif m == 0:
            poly_from_interior += c_m * (n - 1)
        else:
            kk = sp.Symbol('kk')
            power_sum = sp.summation(kk**m, (kk, 1, n - 1))
            poly_from_interior += c_m * power_sum

    ap_nn = coeffs_sym[p].subs(k, n)
    endpoint = cancel(2 * n * ap_nn)
    rational_total = cancel(endpoint + poly_from_interior)

    return betas, rational_total


all_betas = {}
all_rationals = {}

for p in range(2, 7):
    print(f"\n--- p={p} ---")
    t0 = time.time()
    betas_p, rat_p = assemble_cp_structure(p)
    dt = time.time() - t0
    print(f"  Extraction took {dt:.1f}s")

    for r in sorted(betas_p):
        print(f"  beta_{r}(n) = {betas_p[r]}")

    # Verify for n=1..25
    ok = True
    for n_val in range(1, 26):
        predicted = rat_p.subs(n, n_val)
        for r in betas_p:
            H_val = sum(sp.Rational(1, kk**r) for kk in range(1, n_val))
            predicted += betas_p[r].subs(n, n_val) * H_val

        exact = sp.Rational(0)
        for j2 in range(1, 2*n_val, 2):
            j_val = sp.Rational(j2, 2)
            k_val = j_val + sp.Rational(1, 2)
            d_val = int(2*j_val + 1) if j2 == 2*n_val - 1 else int(2*(2*j_val + 1))
            ap_val = coeffs_sym[p].subs([(n, n_val), (k, k_val)])
            exact += d_val * simplify(ap_val)

        if simplify(predicted - exact) != 0:
            ok = False
            print(f"  MISMATCH at n={n_val}")

    print(f"  Verified n=1..25: {'PASS' if ok else 'FAIL'}")
    all_betas[p] = betas_p
    all_rationals[p] = rat_p

# ============================================================
# PART 3: Assemble D_p from Euler sums (FIXED series order)
# ============================================================
print("\n" + "=" * 70)
print("PART 3: COMPUTE D_p (FIXED)")
print("=" * 70)

def euler_sum_hurwitz(r, s):
    old_dps = mp.dps
    mp.dps = 310
    def term(kk):
        return mpmath.hurwitz(s, kk) / mpf(kk)**r
    result = nsum(term, [1, mpinf], method='richardson')
    mp.dps = old_dps
    return result


def compute_Dp(p, betas_p, rat_p):
    """Compute D_p with INCREASED series order for rational/harmonic decomposition."""
    SERIES_ORDER = 20  # enough for 1/n^{19}

    harmonic_contribs = {}
    for r in betas_p:
        expr = betas_p[r]
        series_r = sp.series(expr, n, sp.oo, n=SERIES_ORDER)
        for j in range(1, 2*SERIES_ORDER):
            c = series_r.coeff(n, -j)
            if c != 0:
                harmonic_contribs[(r, j)] = sp.nsimplify(c)

    rational_contribs = {}
    rat_series = sp.series(rat_p, n, sp.oo, n=SERIES_ORDER)
    for j in range(1, 2*SERIES_ORDER):
        c = rat_series.coeff(n, -j)
        if c != 0:
            rational_contribs[j] = sp.nsimplify(c)

    print(f"\n  p={p}:")
    print(f"  Harmonic terms: {len(harmonic_contribs)}")
    for (r, s), c in sorted(harmonic_contribs.items()):
        print(f"    ({c}) * H^{{({r})}} / n^{s}  -> ({c}) * [S_{{{r},{s}}} - z({r+s})]")
    print(f"  Rational terms: {len(rational_contribs)}")
    for j, c in sorted(rational_contribs.items()):
        print(f"    ({c}) / n^{j}  -> ({c}) * z({j})")

    D_val = mpf(0)

    for j, c in rational_contribs.items():
        if j == 1:
            print(f"  WARNING: divergent 1/n term!")
            continue
        D_val += mpf(sp.Rational(c)) * z(j)

    for (r, s), c in sorted(harmonic_contribs.items()):
        c_f = mpf(sp.Rational(c))
        t0 = time.time()
        S_rs = euler_sum_hurwitz(r, s)
        dt = time.time() - t0
        contrib = c_f * (S_rs - z(r + s))
        print(f"    S_{{{r},{s}}} in {dt:.1f}s, contrib = {nstr(contrib, 30)}")
        D_val += contrib

    return D_val


known = {
    2: lambda: mpf(-5)/4*z(2) + z(3),
    3: lambda: mpf(19)/8*z(4) - mpf(11)/4*z(5),
    4: lambda: mpf(-205)/64*z(6) + mpf(71)/8*z(7) - mpf(9)/2*z(3)*z(4),
    5: lambda: mpf(497)/128*z(8) - mpf(467)/16*z(9) + mpf(385)/32*z(3)*z(6) + mpf(75)/8*z(4)*z(5),
}

Dp_values = {}
for p in range(2, 7):
    print(f"\n{'='*50}")
    D_val = compute_Dp(p, all_betas[p], all_rationals[p])
    Dp_values[p] = D_val
    print(f"\n  D_{p} = {nstr(D_val, 60)}")

    if p in known:
        kv = known[p]()
        diff = abs(D_val - kv)
        digits = -int(mpmath.log10(diff)) if diff > 0 else 300
        print(f"  Cross-check vs known: {digits} matching digits")

# ============================================================
# PART 4: PSLQ with REDUCED basis
# ============================================================
print("\n" + "=" * 70)
print("PART 4: PSLQ (REDUCED BASIS)")
print("=" * 70)

pslq_results = {}

for p in range(2, 7):
    w = 2*p - 1
    D_val = Dp_values[p]

    # Reduced basis: z(2p-2), z(2p-1), products z(a)z(b) with a+b=w, 3 <= a <= b
    basis_vals = [D_val]
    basis_names = [f"D_{p}"]

    basis_vals.append(z(w-1))  # z(2p-2), even zeta
    basis_names.append(f"z({w-1})")

    basis_vals.append(z(w))    # z(2p-1), odd zeta
    basis_names.append(f"z({w})")

    # Products with a+b = w, a >= 2
    for a in range(2, w//2 + 1):
        b = w - a
        if b > a or (b == a and a >= 2):
            basis_vals.append(z(a) * z(b))
            basis_names.append(f"z({a})z({b})")

    print(f"\n--- D_{p} (weight {w}, {len(basis_vals)} basis: {', '.join(basis_names[1:])}) ---")

    rel = None
    for mc in [1000, 10000, 100000, 1000000, 10000000]:
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
    decomp = ' + '.join(terms)
    print(f"  D_{p} = {decomp}")

    res = sum(c*v for c, v in zip(rel, basis_vals))
    print(f"  Residual: {nstr(res, 5)}")

    z2_target = f"z(2)z({w-2})"
    z2_found = any(nm == z2_target for nm in basis_names[1:])
    z2_coeff = 0
    if z2_found:
        for c, nm in zip(rel[1:], basis_names[1:]):
            if nm == z2_target:
                z2_coeff = c
    if z2_coeff == 0:
        print(f"  *** z(2)*z({w-2}) = 0  CANCELLED ***")
    else:
        print(f"  *** z(2)*z({w-2}) PRESENT ***")

    pslq_results[p] = {
        'relation': [int(x) for x in rel],
        'basis': basis_names,
        'decomposition': decomp,
        'z2_odd_cancelled': (z2_coeff == 0),
    }

# ============================================================
# PART 5: CANCELLATION TABLE with reduced Euler sum PSLQ
# ============================================================
print("\n" + "=" * 70)
print("PART 5: CANCELLATION TABLE")
print("=" * 70)

def pslq_euler_sum_reduced(r, s):
    """Decompose S_{r,s} using reduced basis."""
    w = r + s  # total weight of the Euler sum

    S_val = euler_sum_hurwitz(r, s)

    # Reduced basis: S_val, z(w), products z(a)z(b) with a+b=w
    basis_vals = [S_val]
    basis_names = ["S"]

    basis_vals.append(z(w))
    basis_names.append(f"z({w})")

    for a in range(2, w//2 + 1):
        b = w - a
        basis_vals.append(z(a) * z(b))
        basis_names.append(f"z({a})z({b})")

    rel = None
    for mc in [1000, 10000, 100000, 1000000, 10000000]:
        rel = mpmath.pslq(basis_vals, maxcoeff=mc)
        if rel is not None:
            break

    if rel is None:
        return None, None, None

    d_coeff = rel[0]
    terms = []
    for c, nm in zip(rel[1:], basis_names[1:]):
        if c != 0:
            terms.append(f"({sp.Rational(-c, d_coeff)}){nm}")
    decomp = ' + '.join(terms)

    # Extract z(2)*z(w-2) coefficient
    z2_target = f"z(2)z({w-2})"
    mu = sp.Rational(0)
    for c, nm in zip(rel[1:], basis_names[1:]):
        if nm == z2_target and d_coeff != 0:
            mu = sp.Rational(-c, d_coeff)

    return mu, decomp, S_val


cancel_results = {}

for p in range(3, 7):
    w = 2*p - 1
    print(f"\n=== p={p}, weight {w}, target: z(2)*z({w-2}) ===")

    cancel_sum = sp.Rational(0)
    all_resolved = True
    terms_info = []

    for r in sorted(all_betas[p]):
        expr = all_betas[p][r]
        # Extract leading coefficient c/n^s
        series_r = sp.series(expr, n, sp.oo, n=20)
        for s in range(1, 30):
            c_val = series_r.coeff(n, -s)
            if c_val != 0:
                c_val = sp.nsimplify(c_val)
                euler_s = s  # S_{r,s} weight = r + s

                print(f"\n  r={r}, s={euler_s}: beta = {c_val}")
                mu, decomp, S_val = pslq_euler_sum_reduced(r, euler_s)

                if mu is None:
                    print(f"    S_{{{r},{euler_s}}} PSLQ failed (weight {r+euler_s})")
                    all_resolved = False
                else:
                    print(f"    S_{{{r},{euler_s}}} = {decomp}")
                    print(f"    mu_{r} (z(2)*z({w-2}) coeff) = {mu}")
                    contrib = c_val * mu
                    cancel_sum += contrib
                    print(f"    beta * mu = {c_val} * {mu} = {contrib}")

                terms_info.append({
                    'r': r, 's': euler_s, 'beta': str(c_val),
                    'mu': str(mu) if mu is not None else None,
                    'decomp': decomp,
                })
                break

    print(f"\n  TOTAL z(2)*z({w-2}) coefficient = {cancel_sum}")
    if cancel_sum == 0 and all_resolved:
        print(f"  *** CANCELLATION CONFIRMED (genuine) ***")
    elif cancel_sum == 0 and not all_resolved:
        print(f"  *** CANCELLATION UNRESOLVED (some PSLQ failed) ***")
    else:
        print(f"  *** NO CANCELLATION ***")

    cancel_results[p] = {
        'target': f'z(2)*z({w-2})',
        'total': str(cancel_sum),
        'all_resolved': all_resolved,
        'terms': terms_info,
    }

# ============================================================
# PART 6: SAVE
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

results = {
    'description': 'D6 analytical (FIXED): corrected series order + reduced PSLQ basis',
    'Dp_values': {},
    'pslq': {},
    'betas': {},
    'cancellation': cancel_results,
}

for p in Dp_values:
    results['Dp_values'][f'D{p}'] = nstr(Dp_values[p], 60)

for p in pslq_results:
    if pslq_results[p]:
        results['pslq'][f'D{p}'] = pslq_results[p]
        print(f"  D_{p} = {pslq_results[p]['decomposition']}  [z(2)*z(odd) cancelled: {pslq_results[p]['z2_odd_cancelled']}]")
    else:
        print(f"  D_{p} = {nstr(Dp_values[p], 40)} (PSLQ failed)")

for p in all_betas:
    results['betas'][f'p{p}'] = {str(r): str(all_betas[p][r]) for r in all_betas[p]}

out = os.path.join(os.path.dirname(__file__), 'data', 'd6_fix.json')
with open(out, 'w') as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nSaved to {out}")
print("Done.")
