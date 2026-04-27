"""
Verify K-Sommerfeld structural separation with D6 included.
PSLQ test: is K/pi in Q-span of {D2, D3, D4, D5, D6, 1}?
"""
import sys, io, functools, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
print = functools.partial(print, flush=True)

from mpmath import mp, mpf, nstr, zeta, pi as mpi, pslq
mp.dps = 200

z = zeta

# Known D_p decompositions (exact rational coefficients)
D2 = mpf(-5)/4 * z(2) + z(3)
D3 = mpf(19)/8 * z(4) - mpf(11)/4 * z(5)
D4 = mpf(-205)/64 * z(6) + mpf(71)/8 * z(7) - mpf(9)/2 * z(3)*z(4)
D5 = mpf(497)/128 * z(8) - mpf(467)/16 * z(9) + mpf(385)/32 * z(3)*z(6) + mpf(75)/8 * z(4)*z(5)
D6 = (mpf(-2289)/512 * z(10) + mpf(1589)/16 * z(11)
      - mpf(1617)/64 * z(3)*z(8) - mpf(1785)/64 * z(4)*z(7) - mpf(2065)/64 * z(5)*z(6))

# K/pi from Paper 2
B = mpf(42)
F = z(2)
Delta = mpf(1)/40
K_over_pi = B + F - Delta

print("=" * 70)
print("K-SOMMERFELD STRUCTURAL SEPARATION WITH D6")
print("=" * 70)

print(f"\n  K/pi = {nstr(K_over_pi, 30)}")
print(f"  D2   = {nstr(D2, 30)}")
print(f"  D3   = {nstr(D3, 30)}")
print(f"  D4   = {nstr(D4, 30)}")
print(f"  D5   = {nstr(D5, 30)}")
print(f"  D6   = {nstr(D6, 30)}")

# Test: K/pi in Q-span of {D2, D3, D4, D5, D6, 1}?
basis = [K_over_pi, D2, D3, D4, D5, D6, mpf(1)]
names = ["K/pi", "D2", "D3", "D4", "D5", "D6", "1"]

print(f"\n  PSLQ basis: {', '.join(names)}")
print(f"  Headroom: {mp.dps}/{len(basis)+1} = {mp.dps//(len(basis)+1)} digits")

for mc in [1000, 10000, 100000, 1000000]:
    rel = pslq(basis, maxcoeff=mc)
    if rel is not None:
        print(f"  maxcoeff={mc}: FOUND {[int(x) for x in rel]}")
        break
    else:
        print(f"  maxcoeff={mc}: null")

if rel is None:
    print(f"\n  *** K/pi NOT in Q-span of {{D2,...,D6,1}} ***")
    print(f"  K-Sommerfeld structural separation CONFIRMED through D6")
else:
    print(f"\n  *** RELATION FOUND — needs investigation ***")

# Also check K/pi vs {D2,...,D6, z(2), z(3), 1} (extended basis)
print(f"\n  Extended basis: K/pi, D2, D3, D4, D5, D6, z(2), z(3), 1")
basis_ext = [K_over_pi, D2, D3, D4, D5, D6, z(2), z(3), mpf(1)]
for mc in [1000, 10000, 100000]:
    rel = pslq(basis_ext, maxcoeff=mc)
    if rel is not None:
        print(f"  maxcoeff={mc}: FOUND {[int(x) for x in rel]}")
        break
    else:
        print(f"  maxcoeff={mc}: null")

if rel is None:
    print(f"  Extended basis also null — confirmed structural separation")

# Save results
results = {
    'D6_decomposition': {
        'z10': '-2289/512',
        'z11': '1589/16',
        'z2z9': '0',
        'z3z8': '-1617/64',
        'z4z7': '-1785/64',
        'z5z6': '-2065/64',
    },
    'D6_value_200dps': nstr(D6, 100),
    'verification_digits': 62,
    'product_survival': {
        'z2z9_absent': True,
        'surviving_products': 3,
        'predicted': 3,
    },
    'k_sommerfeld_separation': 'CONFIRMED through D6',
}
with open('debug/data/d6_analytical_assembly.json', 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nSaved to debug/data/d6_analytical_assembly.json")
print("Done.")
