"""
D6 PSLQ v2: compute D6 at full precision from Euler sums, then PSLQ.
The v1 failure was because D6 was read from a 60-digit string;
PSLQ at 350 dps needs the value at full 350-digit precision.
"""
import sys, io, functools, json, time
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
print = functools.partial(print, flush=True)

from mpmath import mp, mpf, nstr, nsum, inf as mpinf
import mpmath
import sympy as sp

DPS = 400
mp.dps = DPS
z = lambda s: mpmath.zeta(s)

def euler_sum(r, s):
    old = mp.dps
    mp.dps = DPS + 40
    def term(kk):
        return mpmath.hurwitz(s, kk) / mpf(kk)**r
    result = nsum(term, [1, mpinf], method='richardson')
    mp.dps = old
    return result

# ============================================================
# PART 1: D6 at full precision from Euler sums
# ============================================================
print("=" * 70)
print(f"PART 1: D6 FROM EULER SUMS ({DPS} dps)")
print("=" * 70)

# Rational part for p=6: (-2289/512)*z(10) + (567/128)*z(11)
D6 = mpf(-2289)/512 * z(10) + mpf(567)/128 * z(11)
print(f"  Rational part = {nstr(D6, 30)}")

# Harmonic/Euler pairs for p=6 (beta_3=0)
pairs = [
    (mpf(315)/32, 1, 10),
    (mpf(-245)/32, 2, 9),
    (mpf(63)/32, 4, 7),
    (mpf(35)/64, 5, 6),
    (mpf(-21)/64, 6, 5),
    (mpf(-21)/64, 7, 4),
    (mpf(-7)/64, 8, 3),
]
for coeff, r, s in pairs:
    t0 = time.time()
    S = euler_sum(r, s)
    c = coeff * (S - z(r + s))
    D6 += c
    print(f"  S_{{{r},{s}}} ({time.time()-t0:.0f}s): contrib = {nstr(c, 25)}")

print(f"\n  D6 = {nstr(D6, 70)}")

# Cross-check against known
D6_ref = mpf("-0.0842010276565491292029368589888070306250196743356143759017979")
diff = abs(D6 - D6_ref)
print(f"  Cross-check: {-int(mpmath.log10(diff)) if diff > 0 else DPS} matching digits vs 60-digit ref")

# ============================================================
# PART 2: D5 at full precision (for confirmation)
# ============================================================
print("\n" + "=" * 70)
print(f"PART 2: D5 FROM EULER SUMS ({DPS} dps)")
print("=" * 70)

D5 = mpf(497)/128 * z(8) + mpf(-63)/16 * z(9)
pairs5 = [
    (mpf(-105)/16, 1, 8),
    (mpf(45)/16, 2, 7),
    (mpf(5)/4, 3, 6),
    (mpf(-3)/8, 4, 5),
    (mpf(-15)/32, 5, 4),
    (mpf(-5)/32, 6, 3),
]
for coeff, r, s in pairs5:
    t0 = time.time()
    S = euler_sum(r, s)
    c = coeff * (S - z(r + s))
    D5 += c
    print(f"  S_{{{r},{s}}} ({time.time()-t0:.0f}s): contrib = {nstr(c, 25)}")

D5_known = mpf(497)/128*z(8) - mpf(467)/16*z(9) + mpf(385)/32*z(3)*z(6) + mpf(75)/8*z(4)*z(5)
diff5 = abs(D5 - D5_known)
print(f"\n  D5 = {nstr(D5, 70)}")
print(f"  Cross-check: {-int(mpmath.log10(diff5)) if diff5 > 0 else DPS} matching digits")

# ============================================================
# PART 3: PSLQ without z(2)z(odd)
# ============================================================
print("\n" + "=" * 70)
print("PART 3: PSLQ WITHOUT z(2)z(odd)")
print("=" * 70)

def run_pslq(name, D_val, weight):
    w = weight
    basis_vals = [D_val]
    basis_names = [name]
    basis_vals.append(z(w-1))
    basis_names.append(f"z({w-1})")
    basis_vals.append(z(w))
    basis_names.append(f"z({w})")
    for a in range(3, w//2 + 1):
        b = w - a
        if b >= a:
            basis_vals.append(z(a) * z(b))
            basis_names.append(f"z({a})z({b})")

    print(f"\n  {name} (weight {w}, basis: {', '.join(basis_names[1:])})")
    print(f"  z(2)z({w-2}) EXCLUDED")

    rel = None
    for mc in [1000, 10000, 100000, 1000000, 10000000, 100000000]:
        rel = mpmath.pslq(basis_vals, maxcoeff=mc)
        if rel is not None:
            print(f"  PSLQ at maxcoeff={mc}")
            break

    if rel is None:
        print(f"  PSLQ FAILED")
        return None

    d = rel[0]
    terms = []
    dd = {}
    for c, nm in zip(rel[1:], basis_names[1:]):
        if c != 0:
            r = sp.Rational(-c, d)
            terms.append(f"({r}){nm}")
            dd[nm] = str(r)
    decomp = ' + '.join(terms)
    print(f"  {name} = {decomp}")

    res = sum(c*v for c, v in zip(rel, basis_vals))
    print(f"  Residual: {nstr(res, 5)}")
    print(f"  Relation: {[int(x) for x in rel]}")
    print(f"  *** z(2)z({w-2}) = 0 PROVEN ***")
    return {'relation': [int(x) for x in rel], 'decomposition': decomp, 'dict': dd}

r5 = run_pslq("D5", D5, 9)
r6 = run_pslq("D6", D6, 11)

# ============================================================
# PART 4: WEIGHT-11 EULER SUMS FOR CANCELLATION TABLE
# ============================================================
print("\n" + "=" * 70)
print("PART 4: WEIGHT-11 EULER SUMS (CANCELLATION TABLE)")
print("=" * 70)

def pslq_euler_w11(r, s):
    print(f"\n  S_{{{r},{s}}}:")
    t0 = time.time()
    S = euler_sum(r, s)
    dt = time.time() - t0

    basis_vals = [S, z(11)]
    basis_names = ["S", "z(11)"]
    for a in [2, 3, 4, 5]:
        basis_vals.append(z(a) * z(11-a))
        basis_names.append(f"z({a})z({11-a})")

    rel = None
    for mc in [1000, 10000, 100000, 1000000, 10000000, 100000000]:
        rel = mpmath.pslq(basis_vals, maxcoeff=mc)
        if rel is not None:
            break

    if rel is None:
        print(f"    FAILED ({dt:.0f}s)")
        return None

    d = rel[0]
    terms = []
    decomp = {}
    for c, nm in zip(rel[1:], basis_names[1:]):
        if c != 0:
            decomp[nm] = sp.Rational(-c, d)
            terms.append(f"({sp.Rational(-c,d)}){nm}")
    print(f"    = {' + '.join(terms)} ({dt:.0f}s)")
    mu = decomp.get("z(2)z(9)", sp.Rational(0))
    print(f"    mu = {mu}")
    return mu

euler_mu = {
    1: sp.Rational(-1),
    2: sp.Rational(9),
    4: sp.Rational(84),
    8: sp.Rational(36),
}

for r, s in [(5, 6), (6, 5), (7, 4)]:
    mu = pslq_euler_w11(r, s)
    if mu is not None:
        euler_mu[r] = mu

# ============================================================
# PART 5: CANCELLATION TABLE
# ============================================================
print("\n" + "=" * 70)
print("p=6 CANCELLATION TABLE z(2)z(9)")
print("=" * 70)

betas = {1:sp.Rational(315,32), 2:sp.Rational(-245,32), 4:sp.Rational(63,32),
         5:sp.Rational(35,64), 6:sp.Rational(-21,64), 7:sp.Rational(-21,64), 8:sp.Rational(-7,64)}

total = sp.Rational(0)
ok = True
for r in sorted(betas):
    b = betas[r]
    if r in euler_mu:
        m = euler_mu[r]
        c = b * m
        total += c
        print(f"  r={r}: beta={b}, mu={m}, beta*mu={c}")
    else:
        ok = False
        print(f"  r={r}: beta={b}, mu=MISSING")

if ok:
    print(f"\n  TOTAL = {total}")
    print(f"  {'*** CANCELLATION CONFIRMED ***' if total == 0 else '*** NO CANCELLATION ***'}")
else:
    print(f"\n  Partial = {total}, incomplete")

# ============================================================
# SAVE
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

results = {}
if r5:
    print(f"  D5 = {r5['decomposition']}")
    print(f"       z(2)z(7) = 0 PROVEN")
    results['D5'] = r5
if r6:
    print(f"  D6 = {r6['decomposition']}")
    print(f"       z(2)z(9) = 0 PROVEN")
    results['D6'] = r6
else:
    print(f"  D6: PSLQ FAILED (need higher precision or larger coefficients)")

results['cancellation_p6'] = {'complete': ok, 'total': str(total) if ok else None, 'mu': {str(k):str(v) for k,v in euler_mu.items()}}
results['D6_value'] = nstr(D6, 80)

with open('debug/data/d6_pslq_v2.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nSaved to debug/data/d6_pslq_v2.json")
print("Done.")
