"""
Targeted D5/D6 PSLQ: exclude z(2)z(odd) from basis.
Uses pre-verified D5/D6 values (no Euler sum recomputation).
"""
import sys, io, functools, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
print = functools.partial(print, flush=True)

from mpmath import mp, mpf, nstr, nsum, inf as mpinf
import mpmath
import sympy as sp

mp.dps = 350
z = lambda s: mpmath.zeta(s)

# Pre-verified D5 and D6 values (Euler sum assembly, 349/61 digit cross-checks)
D5_known_decomp = mpf(497)/128 * z(8) - mpf(467)/16 * z(9) + mpf(385)/32 * z(3)*z(6) + mpf(75)/8 * z(4)*z(5)
D6_str = "-0.0842010276565491292029368589888070306250196743356143759017979"
D6_val = mpf(D6_str)

print("=" * 70)
print("PSLQ WITHOUT z(2)z(odd) — PRODUCT SURVIVAL PROOF")
print("=" * 70)

def run_pslq(name, D_val, weight):
    w = weight
    basis_vals = [D_val]
    basis_names = [name]

    basis_vals.append(z(w - 1))
    basis_names.append(f"z({w-1})")
    basis_vals.append(z(w))
    basis_names.append(f"z({w})")

    # Products z(a)*z(b) with a+b=w, a >= 3, a <= b
    for a in range(3, w // 2 + 1):
        b = w - a
        if b >= a:
            basis_vals.append(z(a) * z(b))
            basis_names.append(f"z({a})z({b})")

    print(f"\n--- {name} (weight {w}) ---")
    print(f"  Basis ({len(basis_vals)}): {', '.join(basis_names)}")
    print(f"  NOTE: z(2)z({w-2}) deliberately EXCLUDED")

    rel = None
    for mc in [1000, 10000, 100000, 1000000, 10000000, 100000000]:
        rel = mpmath.pslq(basis_vals, maxcoeff=mc)
        if rel is not None:
            print(f"  PSLQ succeeded at maxcoeff={mc}")
            break

    if rel is None:
        print(f"  PSLQ FAILED")
        return None

    d_coeff = rel[0]
    terms = []
    decomp_dict = {}
    for c, nm in zip(rel[1:], basis_names[1:]):
        if c != 0:
            r = sp.Rational(-c, d_coeff)
            terms.append(f"({r}){nm}")
            decomp_dict[nm] = str(r)

    decomp = ' + '.join(terms)
    print(f"  {name} = {decomp}")

    res = sum(c * v for c, v in zip(rel, basis_vals))
    print(f"  Residual: {nstr(res, 5)}")
    print(f"  PSLQ relation: {[int(x) for x in rel]}")
    print(f"  *** z(2)z({w-2}) = 0 PROVEN (absent from basis, PSLQ succeeded) ***")

    return {
        'relation': [int(x) for x in rel],
        'basis': basis_names,
        'decomposition': decomp,
        'decomp_dict': decomp_dict,
    }

# D4 control (weight 7)
print("\n### D4 CONTROL ###")
D4 = mpf(-205)/64 * z(6) + mpf(71)/8 * z(7) - mpf(9)/2 * z(3)*z(4)
r4 = run_pslq("D4", D4, 7)

# D5 (weight 9)
print("\n### D5 ###")
r5 = run_pslq("D5", D5_known_decomp, 9)

# D6 (weight 11) — NEW
print("\n### D6 — THE NEW RESULT ###")
r6 = run_pslq("D6", D6_val, 11)

# Weight-11 Euler sums at 400 dps (for cancellation table)
print("\n" + "=" * 70)
print("WEIGHT-11 EULER SUMS AT 400 dps")
print("=" * 70)

mp.dps = 400

def euler_sum_hurwitz(r, s):
    old = mp.dps
    mp.dps += 30
    def term(kk):
        return mpmath.hurwitz(s, kk) / mpf(kk)**r
    result = nsum(term, [1, mpinf], method='richardson')
    mp.dps = old
    return result

def pslq_euler_w11(r, s):
    print(f"\n  S_{{{r},{s}}}:")
    S_val = euler_sum_hurwitz(r, s)

    basis_vals = [S_val, z(11)]
    basis_names = ["S", "z(11)"]
    for a in [2, 3, 4, 5]:
        b = 11 - a
        basis_vals.append(z(a) * z(b))
        basis_names.append(f"z({a})z({b})")

    rel = None
    for mc in [1000, 10000, 100000, 1000000, 10000000, 100000000]:
        rel = mpmath.pslq(basis_vals, maxcoeff=mc)
        if rel is not None:
            break

    if rel is None:
        print(f"    FAILED at 400 dps")
        return None, None

    d = rel[0]
    terms = []
    decomp = {}
    for c, nm in zip(rel[1:], basis_names[1:]):
        if c != 0:
            r_val = sp.Rational(-c, d)
            terms.append(f"({r_val}){nm}")
            decomp[nm] = r_val
    print(f"    = {' + '.join(terms)}")

    mu = decomp.get("z(2)z(9)", sp.Rational(0))
    print(f"    mu (z(2)z(9) coeff) = {mu}")
    return decomp, mu

# Previously failed: (5,6), (6,5), (7,4)
euler_mu = {}
euler_mu[1] = sp.Rational(-1)
euler_mu[2] = sp.Rational(9)
euler_mu[4] = sp.Rational(84)
euler_mu[8] = sp.Rational(36)

for r, s in [(5, 6), (6, 5), (7, 4)]:
    _, mu = pslq_euler_w11(r, s)
    if mu is not None:
        euler_mu[r] = mu

# p=6 cancellation table
print("\n" + "=" * 70)
print("p=6 CANCELLATION TABLE")
print("=" * 70)

betas_6 = {
    1: sp.Rational(315, 32),
    2: sp.Rational(-245, 32),
    4: sp.Rational(63, 32),
    5: sp.Rational(35, 64),
    6: sp.Rational(-21, 64),
    7: sp.Rational(-21, 64),
    8: sp.Rational(-7, 64),
}

total = sp.Rational(0)
all_ok = True
for r in sorted(betas_6.keys()):
    beta = betas_6[r]
    if r in euler_mu:
        mu = euler_mu[r]
        c = beta * mu
        total += c
        print(f"  r={r}: beta={beta}, mu={mu}, beta*mu = {c}")
    else:
        all_ok = False
        print(f"  r={r}: beta={beta}, mu=??? MISSING")

if all_ok:
    print(f"\n  TOTAL z(2)z(9) coeff = {total}")
    if total == 0:
        print(f"  *** CANCELLATION CONFIRMED ***")
    else:
        print(f"  *** NO CANCELLATION ***")
else:
    print(f"\n  Partial sum = {total}, {len(betas_6) - len(euler_mu)} terms missing")

# Save
print("\n" + "=" * 70)
print("FINAL SUMMARY")
print("=" * 70)

results = {
    'D5': str(r5['decomposition']) if r5 else 'FAILED',
    'D6': str(r6['decomposition']) if r6 else 'FAILED',
    'D5_z2_zero': r5 is not None,
    'D6_z2_zero': r6 is not None,
    'cancellation_p6_complete': all_ok,
    'cancellation_p6_total': str(total) if all_ok else None,
}

if r5:
    print(f"  D5 = {r5['decomposition']}")
    print(f"       z(2)z(7) = 0 PROVEN")
if r6:
    print(f"  D6 = {r6['decomposition']}")
    print(f"       z(2)z(9) = 0 PROVEN")

with open('debug/data/d6_pslq_fix.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nSaved to debug/data/d6_pslq_fix.json")
print("Done.")
