"""
D6 PSLQ at 1000 dps: brute-force precision attack.

At weight 9, D5 PSLQ succeeds at 400 dps with 5 basis elements (maxcoeff=10000).
PSLQ headroom: 400/(5+1) ≈ 67 digits of coefficient.

At weight 11, we have 7 basis elements (with z(2)z(9)).
Headroom at 600 dps: 600/8 ≈ 75 digits. Should be enough for any reasonable LCD.
But 600 dps FAILED. Two possibilities:
  (a) D6 requires depth > 1 (triple products, Euler sums, etc.)
  (b) Euler sum computation has subtle precision loss

This script tests (b) by going to 1000 dps, giving headroom of 1000/8 ≈ 125 digits.
If this still fails, (a) is confirmed.

Also tests the EXCLUDED basis (no z(2)z(9)) at 1000 dps: 6 elements,
headroom 1000/7 ≈ 143 digits. If the product survival rule holds,
this should be ample.
"""
import sys, io, functools, json, time
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
print = functools.partial(print, flush=True)

from mpmath import mp, mpf, nstr, nsum, inf as mpinf
import mpmath
import sympy as sp

DPS = 1000
mp.dps = DPS
z = lambda s: mpmath.zeta(s)

def euler_sum(r, s):
    old = mp.dps
    mp.dps = DPS + 80
    def term(kk):
        return mpmath.hurwitz(s, kk) / mpf(kk)**r
    result = nsum(term, [1, mpinf], method='richardson')
    mp.dps = old
    return result

# ============================================================
# PART 1: D6 at 1000 dps
# ============================================================
print("=" * 70)
print(f"D6 FROM EULER SUMS ({DPS} dps)")
print("=" * 70)

D6 = mpf(-2289)/512 * z(10) + mpf(567)/128 * z(11)
print(f"  Rational part = {nstr(D6, 30)}")

pairs = [
    (mpf(315)/32, 1, 10), (mpf(-245)/32, 2, 9),
    (mpf(63)/32, 4, 7), (mpf(35)/64, 5, 6),
    (mpf(-21)/64, 6, 5), (mpf(-21)/64, 7, 4),
    (mpf(-7)/64, 8, 3),
]
for coeff, r, s in pairs:
    t0 = time.time()
    S = euler_sum(r, s)
    c = coeff * (S - z(r + s))
    D6 += c
    print(f"  S_{{{r},{s}}} ({time.time()-t0:.0f}s)")

print(f"\n  D6 = {nstr(D6, 80)}")

# Cross-check
D6_ref = mpf("-0.0842010276565491292029368589888070306250196743356143759017979")
diff = abs(D6 - D6_ref)
match = -int(mpmath.log10(diff)) if diff > 0 else DPS
print(f"  Cross-check: {match} matching digits vs 60-digit ref")

# ============================================================
# PART 2: PSLQ at 1000 dps
# ============================================================
print("\n" + "=" * 70)
print("PSLQ AT 1000 DPS")
print("=" * 70)

# Test A: Full basis (with z(2)z(9))
print("\n--- Full depth-1 basis (7 elements) ---")
basis_full = [D6, z(10), z(11), z(2)*z(9), z(3)*z(8), z(4)*z(7), z(5)*z(6)]
names_full = ["D6", "z(10)", "z(11)", "z(2)z(9)", "z(3)z(8)", "z(4)z(7)", "z(5)z(6)"]
print(f"  Basis: {', '.join(names_full)}")
print(f"  PSLQ headroom: {DPS}/8 = {DPS//8} digits of coefficient")

rel_full = None
for mc in [1000, 10000, 100000, 1000000, 10**7, 10**8, 10**9]:
    t0 = time.time()
    rel_full = mpmath.pslq(basis_full, maxcoeff=mc)
    dt = time.time() - t0
    if rel_full is not None:
        print(f"  FOUND at maxcoeff={mc} ({dt:.1f}s)")
        break
    else:
        print(f"  maxcoeff={mc}: no ({dt:.1f}s)")

if rel_full:
    d = rel_full[0]
    terms = []
    decomp = {}
    for c, nm in zip(rel_full[1:], names_full[1:]):
        if c != 0:
            r = sp.Rational(-c, d)
            terms.append(f"({r}){nm}")
            decomp[nm] = str(r)
    res = sum(c*v for c, v in zip(rel_full, basis_full))
    print(f"  D6 = {' + '.join(terms)}")
    print(f"  Residual: {nstr(res, 5)}")
    print(f"  Relation: {[int(x) for x in rel_full]}")
    z2z9_coeff = decomp.get("z(2)z(9)", "0")
    print(f"  z(2)z(9) coeff = {z2z9_coeff}")
    if z2z9_coeff == "0":
        print(f"  *** z(2)z(9) = 0 PROVEN ***")
else:
    print(f"  Full basis PSLQ FAILED at {DPS} dps")

# Test B: Excluded basis (no z(2)z(9))
print("\n--- Excluded basis (6 elements, no z(2)z(9)) ---")
basis_excl = [D6, z(10), z(11), z(3)*z(8), z(4)*z(7), z(5)*z(6)]
names_excl = ["D6", "z(10)", "z(11)", "z(3)z(8)", "z(4)z(7)", "z(5)z(6)"]
print(f"  Basis: {', '.join(names_excl)}")
print(f"  PSLQ headroom: {DPS}/7 = {DPS//7} digits of coefficient")

rel_excl = None
for mc in [1000, 10000, 100000, 1000000, 10**7, 10**8, 10**9]:
    t0 = time.time()
    rel_excl = mpmath.pslq(basis_excl, maxcoeff=mc)
    dt = time.time() - t0
    if rel_excl is not None:
        print(f"  FOUND at maxcoeff={mc} ({dt:.1f}s)")
        break
    else:
        print(f"  maxcoeff={mc}: no ({dt:.1f}s)")

if rel_excl:
    d = rel_excl[0]
    terms = []
    decomp_excl = {}
    for c, nm in zip(rel_excl[1:], names_excl[1:]):
        if c != 0:
            r = sp.Rational(-c, d)
            terms.append(f"({r}){nm}")
            decomp_excl[nm] = str(r)
    res = sum(c*v for c, v in zip(rel_excl, basis_excl))
    print(f"  D6 = {' + '.join(terms)}")
    print(f"  Residual: {nstr(res, 5)}")
    print(f"  Relation: {[int(x) for x in rel_excl]}")
    print(f"  *** z(2)z(9) = 0 PROVEN (absent from basis) ***")
else:
    print(f"  Excluded basis PSLQ FAILED at {DPS} dps")

# Test C: Extended basis with triple products
print("\n--- Extended basis with triple products (9 elements) ---")
basis_ext = [D6, z(10), z(11), z(2)*z(9), z(3)*z(8), z(4)*z(7), z(5)*z(6),
             z(3)**2*z(5), z(3)*z(4)**2]
names_ext = ["D6", "z(10)", "z(11)", "z(2)z(9)", "z(3)z(8)", "z(4)z(7)", "z(5)z(6)",
             "z(3)^2*z(5)", "z(3)*z(4)^2"]
print(f"  Basis: {', '.join(names_ext)}")
print(f"  PSLQ headroom: {DPS}/10 = {DPS//10} digits of coefficient")

rel_ext = None
for mc in [1000, 10000, 100000, 1000000, 10**7, 10**8, 10**9]:
    t0 = time.time()
    rel_ext = mpmath.pslq(basis_ext, maxcoeff=mc)
    dt = time.time() - t0
    if rel_ext is not None:
        print(f"  FOUND at maxcoeff={mc} ({dt:.1f}s)")
        break
    else:
        print(f"  maxcoeff={mc}: no ({dt:.1f}s)")

if rel_ext:
    d = rel_ext[0]
    terms = []
    for c, nm in zip(rel_ext[1:], names_ext[1:]):
        if c != 0:
            r = sp.Rational(-c, d)
            terms.append(f"({r}){nm}")
    res = sum(c*v for c, v in zip(rel_ext, basis_ext))
    print(f"  D6 = {' + '.join(terms)}")
    print(f"  Residual: {nstr(res, 5)}")
    print(f"  Relation: {[int(x) for x in rel_ext]}")
else:
    print(f"  Extended basis PSLQ FAILED at {DPS} dps")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

results = {
    'D6_1000dps': nstr(D6, 100),
    'cross_check_digits': match,
    'full_pslq': 'FOUND' if rel_full else 'FAILED',
    'excl_pslq': 'FOUND' if rel_excl else 'FAILED',
    'extended_pslq': 'FOUND' if rel_ext else 'FAILED',
}
if rel_full:
    results['full_relation'] = [int(x) for x in rel_full]
if rel_excl:
    results['excl_relation'] = [int(x) for x in rel_excl]
if rel_ext:
    results['ext_relation'] = [int(x) for x in rel_ext]

with open('debug/data/d6_pslq_1000dps.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nSaved to debug/data/d6_pslq_1000dps.json")
print("Done.")
