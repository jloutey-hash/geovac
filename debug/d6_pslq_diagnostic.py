"""
D6 PSLQ diagnostic: try FULL basis (including z(2)z(9)) at 400+ dps.
This tests whether the failure is:
  (a) z(2)z(9) is actually needed (product survival rule fails)
  (b) precision issue (need more dps)
  (c) coefficient size issue
"""
import sys, io, functools, json, time
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
print = functools.partial(print, flush=True)

from mpmath import mp, mpf, nstr, nsum, inf as mpinf
import mpmath
import sympy as sp

def euler_sum(r, s, dps):
    old = mp.dps
    mp.dps = dps + 40
    def term(kk):
        return mpmath.hurwitz(s, kk) / mpf(kk)**r
    result = nsum(term, [1, mpinf], method='richardson')
    mp.dps = old
    return result

def compute_D6(dps):
    mp.dps = dps
    z = lambda s: mpmath.zeta(s)
    D6 = mpf(-2289)/512 * z(10) + mpf(567)/128 * z(11)
    pairs = [
        (mpf(315)/32, 1, 10), (mpf(-245)/32, 2, 9),
        (mpf(63)/32, 4, 7), (mpf(35)/64, 5, 6),
        (mpf(-21)/64, 6, 5), (mpf(-21)/64, 7, 4),
        (mpf(-7)/64, 8, 3),
    ]
    for coeff, r, s in pairs:
        t0 = time.time()
        S = euler_sum(r, s, dps)
        c = coeff * (S - z(r + s))
        D6 += c
        print(f"  S_{{{r},{s}}} ({time.time()-t0:.0f}s)")
    return D6

def run_pslq_test(name, D_val, weight, include_z2_odd, dps):
    mp.dps = dps
    z = lambda s: mpmath.zeta(s)
    w = weight
    basis_vals = [D_val]
    basis_names = [name]

    basis_vals.append(z(w-1))
    basis_names.append(f"z({w-1})")
    basis_vals.append(z(w))
    basis_names.append(f"z({w})")

    start_a = 2 if include_z2_odd else 3
    for a in range(start_a, w//2 + 1):
        b = w - a
        if b >= a:
            basis_vals.append(z(a) * z(b))
            basis_names.append(f"z({a})z({b})")

    print(f"\n  Basis ({len(basis_vals)}): {', '.join(basis_names)}")
    if include_z2_odd:
        print(f"  z(2)z({w-2}) INCLUDED")
    else:
        print(f"  z(2)z({w-2}) EXCLUDED")

    rel = None
    for mc in [1000, 10000, 100000, 1000000, 10000000, 100000000]:
        rel = mpmath.pslq(basis_vals, maxcoeff=mc)
        if rel is not None:
            print(f"  PSLQ at maxcoeff={mc}")
            break

    if rel is None:
        print(f"  PSLQ FAILED at {dps} dps")
        return None, None

    d = rel[0]
    terms = []
    decomp = {}
    for c, nm in zip(rel[1:], basis_names[1:]):
        if c != 0:
            r = sp.Rational(-c, d)
            terms.append(f"({r}){nm}")
            decomp[nm] = str(r)

    res = sum(c*v for c, v in zip(rel, basis_vals))
    print(f"  {name} = {' + '.join(terms)}")
    print(f"  Residual: {nstr(res, 5)}")
    print(f"  Relation: {[int(x) for x in rel]}")

    z2_key = f"z(2)z({w-2})"
    z2_coeff = decomp.get(z2_key, "0")
    print(f"  z(2)z({w-2}) coeff = {z2_coeff}")

    return [int(x) for x in rel], decomp

# ============================================================
# TEST 1: D6 with FULL basis at 400 dps
# ============================================================
print("=" * 70)
print("TEST 1: D6 at 400 dps, FULL basis (z(2)z(9) INCLUDED)")
print("=" * 70)

DPS = 400
print(f"\nComputing D6 from Euler sums at {DPS} dps...")
D6_400 = compute_D6(DPS)
print(f"\n  D6 = {nstr(D6_400, 70)}")

print("\n--- PSLQ with z(2)z(9) INCLUDED ---")
rel_full, decomp_full = run_pslq_test("D6", D6_400, 11, True, DPS)

print("\n--- PSLQ with z(2)z(9) EXCLUDED ---")
rel_excl, decomp_excl = run_pslq_test("D6", D6_400, 11, False, DPS)

# ============================================================
# TEST 2: D6 at 600 dps if 400 failed
# ============================================================
if rel_full is None and rel_excl is None:
    print("\n" + "=" * 70)
    print("TEST 2: D6 at 600 dps (both bases failed at 400)")
    print("=" * 70)

    DPS2 = 600
    print(f"\nComputing D6 from Euler sums at {DPS2} dps...")
    D6_600 = compute_D6(DPS2)
    print(f"\n  D6 = {nstr(D6_600, 70)}")

    print("\n--- PSLQ with z(2)z(9) INCLUDED ---")
    rel_full2, decomp_full2 = run_pslq_test("D6", D6_600, 11, True, DPS2)

    print("\n--- PSLQ with z(2)z(9) EXCLUDED ---")
    rel_excl2, decomp_excl2 = run_pslq_test("D6", D6_600, 11, False, DPS2)
else:
    rel_full2, decomp_full2 = None, None
    rel_excl2, decomp_excl2 = None, None

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

results = {
    'D6_400dps': nstr(D6_400, 80),
    'pslq_400_full': {'relation': rel_full, 'decomp': decomp_full} if rel_full else 'FAILED',
    'pslq_400_excl': {'relation': rel_excl, 'decomp': decomp_excl} if rel_excl else 'FAILED',
}

if rel_full:
    z2_coeff = decomp_full.get('z(2)z(9)', '0')
    print(f"  400 dps FULL: SUCCEEDED, z(2)z(9) coeff = {z2_coeff}")
    if z2_coeff == '0':
        print(f"  *** z(2)z(9) = 0 PROVEN (present in basis, coefficient zero) ***")
    else:
        print(f"  *** PRODUCT SURVIVAL RULE VIOLATED: z(2)z(9) != 0 ***")
elif rel_excl:
    print(f"  400 dps EXCLUDED: SUCCEEDED")
    print(f"  *** z(2)z(9) = 0 PROVEN (absent from basis, PSLQ succeeded) ***")
else:
    print(f"  400 dps: BOTH FAILED")

if rel_full2 or rel_excl2:
    results['pslq_600_full'] = {'relation': rel_full2, 'decomp': decomp_full2} if rel_full2 else 'FAILED'
    results['pslq_600_excl'] = {'relation': rel_excl2, 'decomp': decomp_excl2} if rel_excl2 else 'FAILED'
    if rel_full2:
        z2_coeff = decomp_full2.get('z(2)z(9)', '0')
        print(f"  600 dps FULL: SUCCEEDED, z(2)z(9) coeff = {z2_coeff}")
    elif rel_excl2:
        print(f"  600 dps EXCLUDED: SUCCEEDED")
    else:
        print(f"  600 dps: BOTH FAILED")

with open('debug/data/d6_pslq_diagnostic.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nSaved to debug/data/d6_pslq_diagnostic.json")
print("Done.")
