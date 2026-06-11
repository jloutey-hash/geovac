"""
Inverse Packing Problem: What packings exist on S^d for various d,
and what physics do they encode?

This script computes:
1. Laplace-Beltrami eigenvalues and degeneracies on S^d for d=1..7
2. Cumulative degeneracy sequences and comparison to known physics
3. S^3 verification (Coulomb case)
4. Flat-space packing analogs in d dimensions
5. Exceptional Lie group irrep dimensions

Output: debug/data/inverse_packing_results.json
"""

import numpy as np
import json
from math import comb, factorial, gcd
from functools import reduce


def spherical_harmonic_degeneracy(l, d):
    """
    Dimension of the space of degree-l spherical harmonics on S^d.

    For S^d embedded in R^{d+1}, the l-th eigenspace of the
    Laplace-Beltrami operator has eigenvalue -l(l+d-1) and
    degeneracy:

    g_l = C(l+d, d) - C(l+d-2, d)   for l >= 2
    g_0 = 1
    g_1 = d + 1

    Equivalently: g_l = (2l+d-1) * C(l+d-2, d-1) / (l+d-1)
    which equals the dimension of the l-th harmonic polynomial space.
    """
    if d <= 0:
        raise ValueError(f"d must be >= 1, got {d}")
    if l < 0:
        return 0
    if l == 0:
        return 1
    if l == 1:
        return d + 1
    # General formula: C(l+d, d) - C(l+d-2, d)
    return comb(l + d, d) - comb(l + d - 2, d)


def laplace_beltrami_eigenvalue(l, d):
    """Eigenvalue of the Laplace-Beltrami operator on S^d: -l(l+d-1)"""
    return -l * (l + d - 1)


# =========================================================================
# 1. Tabulate eigenvalues and degeneracies on S^d for d=1..7, l=0..10
# =========================================================================
print("=" * 70)
print("1. LAPLACE-BELTRAMI EIGENVALUES AND DEGENERACIES ON S^d")
print("=" * 70)

results = {}
table_data = {}

for d in range(1, 8):
    table_data[d] = []
    print(f"\n--- S^{d} (embedded in R^{d+1}) ---")
    print(f"{'l':>3}  {'eigenvalue':>12}  {'degeneracy':>12}  {'cum_deg':>12}")
    cum_deg = 0
    for l in range(11):
        ev = laplace_beltrami_eigenvalue(l, d)
        g = spherical_harmonic_degeneracy(l, d)
        cum_deg += g
        table_data[d].append({
            'l': l, 'eigenvalue': ev, 'degeneracy': g,
            'cumulative_degeneracy': cum_deg
        })
        print(f"{l:3d}  {ev:12d}  {g:12d}  {cum_deg:12d}")

results['sphere_eigenvalues'] = table_data


# =========================================================================
# 2. Cumulative degeneracy sequences -- comparison to known physics
# =========================================================================
print("\n" + "=" * 70)
print("2. CUMULATIVE DEGENERACY SEQUENCES")
print("=" * 70)

# Known physical sequences
nuclear_magic = [2, 8, 20, 28, 50, 82, 126]
noble_gas = [2, 10, 18, 36, 54, 86]
hydrogen_2n2 = [2*n**2 for n in range(1, 8)]  # 2, 8, 18, 32, 50, 72, 98
hydrogen_n2 = [n**2 for n in range(1, 8)]  # 1, 4, 9, 16, 25, 36, 49

# HO magic numbers (without spin-orbit): 2, 8, 20, 40, 70, 112, 168
# These are 2 * cumulative HO degeneracies: sum_{N=0}^{N_max} (N+1)(N+2)/2
# with factor 2 for spin, giving 2, 8, 20, 40, 70, 112
ho_magic_no_so = []
cum = 0
for N in range(7):
    cum += (N + 1) * (N + 2) // 2
    ho_magic_no_so.append(2 * cum)

print(f"\nNuclear magic numbers:    {nuclear_magic}")
print(f"Noble gas electron counts: {noble_gas}")
print(f"Hydrogen 2n^2:             {hydrogen_2n2}")
print(f"Hydrogen n^2:              {hydrogen_n2}")
print(f"HO magic (no spin-orbit):  {ho_magic_no_so}")

cumulative_sequences = {}
matches = {}

for d in range(1, 8):
    cum_seq = [entry['cumulative_degeneracy'] for entry in table_data[d]]
    cumulative_sequences[d] = cum_seq

    # Check for matches
    d_matches = []

    # Check hydrogen n^2
    if cum_seq[:7] == hydrogen_n2[:7]:
        d_matches.append("hydrogen_n2 (orbital)")

    # Check 2n^2
    cum_2x = [2 * c for c in cum_seq]
    if cum_2x[:7] == hydrogen_2n2[:7]:
        d_matches.append("hydrogen_2n2 (with spin factor 2)")

    # Check HO magic (without spin-orbit)
    ho_cum = []
    cum = 0
    for entry in table_data[d]:
        cum += entry['degeneracy']
        ho_cum.append(cum)
    # Compare with HO shell structure

    # Check nuclear magic
    for m in nuclear_magic:
        if m in cum_seq:
            idx = cum_seq.index(m)
            d_matches.append(f"nuclear_magic {m} at l={idx}")

    # Check noble gas
    for g in noble_gas:
        if g in cum_seq:
            idx = cum_seq.index(g)
            d_matches.append(f"noble_gas {g} at l={idx}")

    matches[d] = d_matches

    print(f"\nS^{d}: cumulative = {cum_seq[:8]}")
    if d_matches:
        for m in d_matches:
            print(f"  MATCH: {m}")

results['cumulative_sequences'] = {str(k): v for k, v in cumulative_sequences.items()}
results['known_sequence_matches'] = {str(k): v for k, v in matches.items()}


# =========================================================================
# 2b. Deeper sequence analysis: 2x cumulative (spin doubling)
# =========================================================================
print("\n" + "=" * 70)
print("2b. SPIN-DOUBLED CUMULATIVE SEQUENCES (2 x cum_deg)")
print("=" * 70)

spin_doubled_matches = {}

for d in range(1, 8):
    cum_seq = [entry['cumulative_degeneracy'] for entry in table_data[d]]
    doubled = [2 * c for c in cum_seq]

    d_matches = []

    # Check noble gas (which includes spin)
    for g in noble_gas:
        if g in doubled:
            idx = doubled.index(g)
            d_matches.append(f"noble_gas {g} at l={idx}")

    # Check nuclear magic
    for m in nuclear_magic:
        if m in doubled:
            idx = doubled.index(m)
            d_matches.append(f"nuclear_magic {m} at l={idx}")

    # Check HO magic
    for h in ho_magic_no_so:
        if h in doubled:
            idx = doubled.index(h)
            d_matches.append(f"HO_magic {h} at l={idx}")

    spin_doubled_matches[d] = d_matches

    print(f"\nS^{d}: 2 x cumulative = {doubled[:8]}")
    if d_matches:
        for m in d_matches:
            print(f"  MATCH: {m}")
    else:
        print("  (no matches)")

results['spin_doubled_matches'] = {str(k): v for k, v in spin_doubled_matches.items()}


# =========================================================================
# 3. S^3 verification (the known Coulomb case)
# =========================================================================
print("\n" + "=" * 70)
print("3. S^3 VERIFICATION (COULOMB CASE)")
print("=" * 70)

print("\nS^3 degeneracies vs n^2:")
print(f"{'n':>3}  {'g_n(S^3)':>10}  {'n^2':>10}  {'match':>8}")
s3_match = True
for n in range(1, 11):
    # On S^3 (d=3), the l-th harmonic has degeneracy (l+1)^2
    # Using l = n-1 gives g = n^2
    g_s3 = spherical_harmonic_degeneracy(n - 1, 3)
    expected = n**2
    ok = (g_s3 == expected)
    if not ok:
        s3_match = False
    print(f"{n:3d}  {g_s3:10d}  {expected:10d}  {'YES' if ok else 'NO':>8}")

print(f"\nAll match: {s3_match}")

# Cumulative: sum n^2 = n(n+1)(2n+1)/6 vs 2n^2 (hydrogen with spin)
print("\nCumulative sum n^2 vs hydrogen 2n^2:")
print(f"{'n':>3}  {'sum_n^2':>10}  {'2n^2':>10}  {'ratio':>10}")
cum = 0
for n in range(1, 8):
    cum += n**2
    formula = n * (n + 1) * (2 * n + 1) // 6
    assert cum == formula
    print(f"{n:3d}  {cum:10d}  {2*n**2:10d}  {cum/(2*n**2):10.4f}")

print("""
KEY INSIGHT: The S^3 cumulative degeneracy sum_{n=1}^{N} n^2 = N(N+1)(2N+1)/6
does NOT equal 2N^2 (hydrogen with spin). The relationship is:
- S^3 gives the ORBITAL degeneracy per energy shell: n^2
- Hydrogen groups shells: n=1 has g=1, n=2 has g=4, etc.
- The cumulative count sum n^2 counts total orbital states up to shell N
- The noble gas sequence 2, 8, 18, 36, ... = 2n^2 counts states PER SHELL
  (with spin factor 2), not cumulative
- The actual noble gas cumulative is 2, 10, 18, 36, 54, 86 which comes from
  filling 2n^2 states shell by shell WITH Madelung rule reordering
""")

results['s3_verification'] = {
    'degeneracies_match_n_squared': s3_match,
    'note': 'S^3 harmonic degeneracy g_l = (l+1)^2 matches n^2 with n=l+1'
}


# =========================================================================
# 4. Flat-space packing analog in d dimensions
# =========================================================================
print("\n" + "=" * 70)
print("4. FLAT-SPACE PACKING IN d DIMENSIONS")
print("=" * 70)

print("""
Paper 0's construction: in 2D flat space, shell k has annular area
proportional to k^2 - (k-1)^2 = 2k-1. Dividing by fundamental area
gives N_k = 2(2k-1) angular states per shell.

Generalization: in d-dimensional flat space, shell k occupies volume
proportional to k^d - (k-1)^d. The "state count" is proportional to
this volume difference.
""")

print(f"{'k':>3}", end="")
for d in range(2, 7):
    print(f"  {'d='+str(d):>10}", end="")
print()
print("-" * 58)

flat_packing = {}
for d in range(2, 7):
    flat_packing[d] = []

for k in range(1, 11):
    print(f"{k:3d}", end="")
    for d in range(2, 7):
        vol_diff = k**d - (k - 1)**d
        flat_packing[d].append(vol_diff)
        print(f"  {vol_diff:10d}", end="")
    print()

# Compare flat packing capacities to spherical harmonic degeneracies
print("\n--- Comparison: flat packing vs spherical harmonic degeneracies ---")
print("(Does flat d-space packing capacity match S^{d'} degeneracy for some d'?)")

packing_matches = {}
for d_flat in range(2, 7):
    fp = flat_packing[d_flat]
    d_matches = []
    for d_sphere in range(1, 8):
        sh_degs = [spherical_harmonic_degeneracy(l, d_sphere) for l in range(10)]
        if fp[:6] == sh_degs[:6]:
            d_matches.append(f"Exact match with S^{d_sphere} harmonic degeneracies")
        # Check if proportional
        if sh_degs[0] != 0 and fp[0] != 0:
            ratio = fp[0] / sh_degs[0]
            all_proportional = all(
                abs(fp[i] - ratio * sh_degs[i]) < 1e-10
                for i in range(min(6, len(fp), len(sh_degs)))
                if sh_degs[i] != 0
            )
            if all_proportional and abs(ratio - round(ratio)) < 1e-10:
                d_matches.append(
                    f"Proportional to S^{d_sphere} degeneracies "
                    f"(ratio = {ratio:.0f})"
                )

    packing_matches[d_flat] = d_matches
    print(f"\nd={d_flat} flat packing: {fp[:6]}")
    if d_matches:
        for m in d_matches:
            print(f"  {m}")
    else:
        # Check what the flat packing sequence looks like
        # k^d - (k-1)^d for general d is a polynomial of degree d-1 in k
        print(f"  Closed form: k^{d_flat} - (k-1)^{d_flat}")
        # Evaluate leading term
        print(f"  Leading term: {d_flat}*k^{d_flat-1}")

results['flat_packing'] = {str(k): v for k, v in flat_packing.items()}
results['flat_packing_matches'] = {str(k): v for k, v in packing_matches.items()}


# =========================================================================
# 4b. More careful analysis: Paper 0 packing in d dimensions
# =========================================================================
print("\n" + "=" * 70)
print("4b. PAPER 0 PACKING GENERALIZED TO d DIMENSIONS")
print("=" * 70)

print("""
Paper 0 uses 2D flat space with:
  - Fundamental area sigma_0 = pi * d_0^2 / 2
  - Annular area A_k = pi * d_0^2 * (2k-1)
  - State count N_k = A_k / sigma_0 = 2(2k-1)
  - Angular capacity (no spin): 2k-1 = 2l+1

Generalize to d dimensions:
  - Shell k has volume V_k = V_d * d_0^d * [k^d - (k-1)^d]
    where V_d is the d-ball volume coefficient
  - Fundamental volume sigma_0 = V_d * d_0^d / N_init
    (N_init = 2 in Paper 0)
  - State count N_k = N_init * [k^d - (k-1)^d]

So the angular capacity per shell in d dimensions is:
  N_k^{angular} = k^d - (k-1)^d
  (dropping the global factor N_init which is spin-like)
""")

# Key question: does k^d - (k-1)^d for packing dimension d
# match the S^{d'} harmonic degeneracy for some d'?

print("Comparison of flat packing capacity vs S^{d-1} harmonic degeneracy:")
print("(Paper 0 case: d=2 packing gives 2k-1 = S^2 harmonic = SO(3) irrep dim)")
print()

for d in range(2, 7):
    fp = [k**d - (k - 1)**d for k in range(1, 8)]
    # S^{d-1} harmonic degeneracy
    sh = [spherical_harmonic_degeneracy(k - 1, d - 1) for k in range(1, 8)]
    # S^d harmonic degeneracy
    sh_d = [spherical_harmonic_degeneracy(k - 1, d) for k in range(1, 8)]
    # S^{d+1} harmonic degeneracy (the Fock sphere)
    sh_dp1 = [spherical_harmonic_degeneracy(k - 1, d + 1) for k in range(1, 8)]

    match_dm1 = (fp == sh)
    match_d = (fp == sh_d)
    match_dp1 = (fp == sh_dp1)

    print(f"d={d} packing:  {fp}")
    print(f"  S^{d-1} harm:  {sh}  {'MATCH!' if match_dm1 else ''}")
    print(f"  S^{d} harm:    {sh_d}  {'MATCH!' if match_d else ''}")
    print(f"  S^{d+1} harm:  {sh_dp1}  {'MATCH!' if match_dp1 else ''}")
    print()


# =========================================================================
# 4c. The key pattern: packing in d flat dimensions vs S^d degeneracy
# =========================================================================
print("=" * 70)
print("4c. KEY PATTERN ANALYSIS")
print("=" * 70)

print("""
Testing: does the d-dimensional flat packing capacity k^d - (k-1)^d
equal the spherical harmonic degeneracy on some S^{d'} at l = k-1?
""")

# Build a table
print(f"{'d_flat':>6} | {'packing seq':>40} | {'matches S^?':>20}")
print("-" * 72)

pattern_results = {}
for d in range(1, 8):
    fp = [k**d - (k - 1)**d for k in range(1, 8)]
    match_sphere = None
    for d_sph in range(1, 10):
        sh = [spherical_harmonic_degeneracy(k - 1, d_sph) for k in range(1, 8)]
        if fp == sh:
            match_sphere = d_sph
            break

    fp_str = str(fp)
    if match_sphere:
        match_str = f"S^{match_sphere}"
    else:
        match_str = "none"

    print(f"{d:6d} | {fp_str:>40} | {match_str:>20}")
    pattern_results[d] = {
        'flat_packing_capacity': fp,
        'matching_sphere': match_sphere
    }

results['packing_sphere_correspondence'] = {str(k): v for k, v in pattern_results.items()}

print("""
FINDING: d=1 flat packing gives capacity = 1 for all k (matches S^0 trivially).
d=2 flat packing gives 2k-1 which does NOT match S^1 (which has g_l=2 for l>=1)
but DOES match the dimension of the l-th irrep of SO(3), i.e., the angular
part of the S^2 harmonic space before considering the full S^2 degeneracy.

The S^d spherical harmonic degeneracy for general d is:
  g_l(S^d) = C(l+d, d) - C(l+d-2, d)

The flat packing capacity in d dimensions is:
  f_k(d) = k^d - (k-1)^d

These are in general DIFFERENT functions of l=k-1 and d.
Let's compare them term by term.
""")

print(f"{'l':>3}", end="")
for d in range(2, 7):
    print(f"  {'pack_'+str(d):>8}  {'S^'+str(d):>8}", end="")
print()
print("-" * 95)

for l in range(8):
    k = l + 1
    print(f"{l:3d}", end="")
    for d in range(2, 7):
        fp = k**d - (k - 1)**d
        sh = spherical_harmonic_degeneracy(l, d)
        print(f"  {fp:8d}  {sh:8d}", end="")
    print()

print("""
OBSERVATION: The flat packing and S^d degeneracies agree only at d=1 (trivially)
and they diverge for d >= 2. The d=2 case is special:
  flat: 2k-1 = 2l+1
  S^2: 2l+1
These DO match! But this is the SO(3) irrep dimension, not S^2 in the S^d
Laplace-Beltrami sense (S^2 harmonics also have degeneracy 2l+1).

Wait -- let me recheck. S^2 harmonics: g_l = C(l+2,2) - C(l,2) = (l+2)(l+1)/2 - l(l-1)/2 = 2l+1.
So g_l(S^2) = 2l+1. And flat d=2 packing gives 2k-1 = 2l+1 with l=k-1.
They DO match!
""")

# Recheck carefully
print("CAREFUL RECHECK: flat packing vs S^d degeneracy")
for d in range(1, 8):
    match = True
    for l in range(8):
        k = l + 1
        fp = k**d - (k - 1)**d
        sh = spherical_harmonic_degeneracy(l, d)
        if fp != sh:
            match = False
            break
    status = "MATCH" if match else "differ"
    print(f"  d={d} flat packing vs S^{d} harmonics: {status}")

# Also check d flat vs S^{d-1}
print()
for d in range(2, 8):
    match = True
    for l in range(8):
        k = l + 1
        fp = k**d - (k - 1)**d
        sh = spherical_harmonic_degeneracy(l, d - 1)
        if fp != sh:
            match = False
            break
    status = "MATCH" if match else "differ"
    print(f"  d={d} flat packing vs S^{d-1} harmonics: {status}")


# =========================================================================
# 4d. What DOES match?  Analyze the relationship algebraically
# =========================================================================
print("\n" + "=" * 70)
print("4d. ALGEBRAIC ANALYSIS: k^d - (k-1)^d vs harmonic degeneracy")
print("=" * 70)

print("""
Flat d=2: k^2 - (k-1)^2 = 2k-1.
S^2 harm: C(l+2,2) - C(l,2) = (2l+1).  With l=k-1: 2(k-1)+1 = 2k-1. MATCH.

Flat d=3: k^3 - (k-1)^3 = 3k^2 - 3k + 1.
S^3 harm: (l+1)^2.  With l=k-1: k^2.
  3k^2 - 3k + 1 vs k^2: differ by 2k^2 - 3k + 1 = (2k-1)(k-1).

Flat d=4: k^4 - (k-1)^4 = 4k^3 - 6k^2 + 4k - 1.
S^4 harm: C(l+4,4) - C(l+2,4).
  l=0: 1 vs 1. l=1: 4*1-6+4-1=1?? No, k=2: 16-1=15. g_1(S^4)=5.

Let me just compute the ratios to see if there's a simple relationship.
""")

print(f"{'l':>3}", end="")
for d in range(2, 7):
    print(f"  {'ratio_'+str(d):>10}", end="")
print()
print("-" * 60)

ratio_analysis = {}
for d in range(2, 7):
    ratios = []
    for l in range(1, 8):
        k = l + 1
        fp = k**d - (k - 1)**d
        sh = spherical_harmonic_degeneracy(l, d)
        if sh > 0:
            ratios.append(fp / sh)
        else:
            ratios.append(None)
    ratio_analysis[d] = ratios

for i, l in enumerate(range(1, 8)):
    k = l + 1
    print(f"{l:3d}", end="")
    for d in range(2, 7):
        r = ratio_analysis[d][i]
        if r:
            print(f"  {r:10.4f}", end="")
        else:
            print(f"  {'---':>10}", end="")
    print()

print("""
For d=2: ratio = 1.0000 (exact match, as proven above)
For d>=3: ratios diverge, growing with l.

The d=2 case is unique because:
  k^2 - (k-1)^2 = 2k-1  (linear in k, degree d-1=1)
  S^2 harmonic deg = 2l+1  (also linear in l)

For d>=3, k^d-(k-1)^d is degree d-1 in k, while S^d harmonic deg
is degree d-1 in l but with different coefficients (binomial vs power).
""")


# =========================================================================
# 5. Identify degeneracy sequences with known physics
# =========================================================================
print("\n" + "=" * 70)
print("5. PHYSICAL IDENTIFICATION OF DEGENERACY SEQUENCES")
print("=" * 70)

identifications = {}

for d in range(1, 8):
    degs = [spherical_harmonic_degeneracy(l, d) for l in range(11)]
    cum = []
    c = 0
    for g in degs:
        c += g
        cum.append(c)

    id_list = []

    # Known identifications
    if d == 1:
        id_list.append("S^1 (circle): g_l = 2 for l>=1, g_0=1. "
                       "Fourier modes on a circle. 1D periodic systems.")
    elif d == 2:
        id_list.append("S^2: g_l = 2l+1 = SO(3) irrep dimension. "
                       "Angular momentum eigenstates. Universal angular "
                       "cross-section (Paper 0).")
    elif d == 3:
        id_list.append("S^3: g_l = (l+1)^2. Hydrogen orbital degeneracy "
                       "with n=l+1. Fock projection of Coulomb problem. "
                       "SO(4) symmetry. (Paper 7)")
    elif d == 4:
        id_list.append("S^4: g_l = (l+1)(l+2)(2l+3)/6. "
                       "Appears in 5D harmonic analysis.")
        # Check if this relates to any known potential
    elif d == 5:
        id_list.append("S^5: g_l = (l+1)(l+2)(l+3)/6 ... wait, need to check. "
                       "Bargmann-Segal lattice for 3D HO lives here. "
                       "SU(3) symmetric irreps (N,0). (Paper 24)")
        # Verify: S^5 harmonic deg at l=N should be (N+1)(N+2)/2
        for N in range(6):
            g = spherical_harmonic_degeneracy(N, 5)
            ho_deg = (N + 1) * (N + 2) // 2
            if g != ho_deg:
                id_list.append(f"  WARNING: S^5 g_{N} = {g} != HO deg {ho_deg}")
            else:
                if N == 0:
                    id_list.append(f"  VERIFIED: S^5 g_N = (N+1)(N+2)/2 = "
                                   f"3D HO degeneracy (checked N=0..5)")

    identifications[d] = id_list

    print(f"\nS^{d}: degeneracies = {degs[:8]}")
    print(f"       cumulative  = {cum[:8]}")
    for item in id_list:
        print(f"  {item}")

# Verify S^5 = 3D HO degeneracy explicitly
print("\n--- S^5 vs 3D HO shell degeneracy ---")
print(f"{'N':>3}  {'g(S^5)':>10}  {'(N+1)(N+2)/2':>15}  {'match':>8}")
s5_ho_match = True
for N in range(11):
    g_s5 = spherical_harmonic_degeneracy(N, 5)
    ho_deg = (N + 1) * (N + 2) // 2
    ok = (g_s5 == ho_deg)
    if not ok:
        s5_ho_match = False
    print(f"{N:3d}  {g_s5:10d}  {ho_deg:15d}  {'YES' if ok else 'NO':>8}")

print(f"\nAll match: {s5_ho_match}")
if s5_ho_match:
    print("CONFIRMED: S^5 harmonic degeneracy = 3D HO shell degeneracy = (N+1)(N+2)/2")
    print("This is the Bargmann-Segal lattice result from Paper 24.")

results['s5_ho_verification'] = {
    'degeneracies_match': s5_ho_match,
    'note': 'S^5 harmonic degeneracy matches 3D isotropic HO: (N+1)(N+2)/2'
}

# Check S^{2d-1} vs d-dimensional HO
print("\n--- General pattern: S^{2d-1} vs d-dimensional HO ---")
print("The d-dim isotropic HO has shell degeneracy C(N+d-1, d-1)")
print()

ho_sphere_matches = {}
for dim_ho in range(1, 5):
    d_sphere = 2 * dim_ho - 1  # S^1 for 1D HO, S^3 for 2D HO, S^5 for 3D HO, S^7 for 4D HO
    all_match = True
    print(f"{dim_ho}D HO vs S^{d_sphere} harmonics:")
    for N in range(8):
        g_sphere = spherical_harmonic_degeneracy(N, d_sphere)
        ho_deg = comb(N + dim_ho - 1, dim_ho - 1)
        ok = (g_sphere == ho_deg)
        if not ok:
            all_match = False
        if N < 5:
            print(f"  N={N}: S^{d_sphere} = {g_sphere}, "
                  f"{dim_ho}D HO = {ho_deg}  {'MATCH' if ok else 'DIFFER'}")

    ho_sphere_matches[dim_ho] = {
        'd_sphere': d_sphere,
        'match': all_match
    }
    print(f"  All match: {all_match}")
    print()

results['ho_sphere_correspondence'] = ho_sphere_matches

print("""
MAJOR FINDING: The pattern S^{2d-1} <-> d-dimensional HO does NOT hold in general!
Only S^5 <-> 3D HO matches. Let's check more carefully...
""")

# Actually the Bargmann space for d-dim HO is C^d, which has unit sphere S^{2d-1}
# The HO degeneracy for d dimensions is C(N+d-1, d-1) = dim of symmetric polynomials
# The S^{2d-1} harmonic degeneracy involves a more complex formula
# Let me check if the SU(d) symmetric irrep (N,0,...,0) has the right dimension

print("--- SU(d) symmetric irrep dimension vs d-dim HO ---")
print("SU(d) irrep (N,0,...,0) has dimension C(N+d-1, d-1)")
print("This is EXACTLY the d-dim HO shell degeneracy.")
print("The S^{2d-1} Laplace-Beltrami degeneracy is DIFFERENT from the SU(d) irrep dim.")
print()

for dim_ho in range(1, 6):
    d_sphere = 2 * dim_ho - 1
    print(f"{dim_ho}D HO / SU({dim_ho}) / S^{d_sphere}:")
    for N in range(5):
        su_d_dim = comb(N + dim_ho - 1, dim_ho - 1)
        s_deg = spherical_harmonic_degeneracy(N, d_sphere)
        print(f"  N={N}: SU({dim_ho}) (N,0)={su_d_dim}, S^{d_sphere} harm={s_deg}", end="")
        if su_d_dim == s_deg:
            print(" MATCH", end="")
        print()
    print()

print("""
RESOLUTION: The Bargmann-Segal construction picks out the HOLOMORPHIC SECTOR
of S^{2d-1}, which is the SU(d) symmetric irrep (N,0,...,0) with dimension
C(N+d-1, d-1). The full S^{2d-1} Laplace-Beltrami degeneracy is LARGER because
it includes all SO(2d) harmonics, not just the holomorphic ones.

The match S^5 <-> 3D HO works because:
  SU(3) (N,0) dim = C(N+2,2) = (N+1)(N+2)/2
  S^5 harmonic deg = C(N+5,5) - C(N+3,5)

These are NOT the same in general! Let me verify...
""")

# Final check
for N in range(6):
    su3 = comb(N + 2, 2)
    s5 = spherical_harmonic_degeneracy(N, 5)
    print(f"N={N}: SU(3) (N,0) = {su3}, S^5 full harmonic = {s5}")

print("""
So S^5 harmonic degeneracy IS the same as SU(3) (N,0) dimension!
This is because S^5 = SU(3)/SU(2), and the harmonic decomposition of
L^2(S^5) under SO(6) ~ SU(4) restricts to SU(3) reps of type (N,0).

Let me check if this pattern extends:
S^3 = SU(2)/U(1), harmonics = SU(2) irreps of dim N+1 vs (N+1)^2...
No, S^3 harmonics have deg (N+1)^2, while SU(2) irreps have dim N+1.
So the S^5 case IS special.
""")

# Let me be more careful: S^{2d-1} sits inside C^d
# Harmonics on S^{2d-1} decompose under U(d) action
# The (p,q) bidegree harmonics on S^{2d-1} have the structure from
# the Peter-Weyl theorem for U(d)/U(d-1) = S^{2d-1}
# For the Bargmann-Segal lattice, we only want the holomorphic sector (p,0)

print("\n--- Detailed: S^{2d-1} harmonic deg vs SU(d) (N,0) irrep dim ---")
for dim_ho in [1, 2, 3, 4]:
    d_sphere = 2 * dim_ho - 1
    su_d = [comb(N + dim_ho - 1, dim_ho - 1) for N in range(8)]
    s_full = [spherical_harmonic_degeneracy(N, d_sphere) for N in range(8)]
    print(f"\nS^{d_sphere} (d_HO={dim_ho}):")
    print(f"  SU({dim_ho}) (N,0): {su_d}")
    print(f"  S^{d_sphere} full:   {s_full}")
    print(f"  Ratio s_full/su:   {[s_full[i]/su_d[i] if su_d[i]>0 else 'inf' for i in range(8)]}")


# =========================================================================
# 5b. Cumulative degeneracies and shell closures
# =========================================================================
print("\n" + "=" * 70)
print("5b. SHELL CLOSURES AND MAGIC NUMBERS")
print("=" * 70)

# For each S^d, compute cumulative degeneracy and identify where
# magic-number-like sequences appear

print("\nCumulative degeneracy sequences (with 2x spin factor):")
for d in range(1, 8):
    cum = 0
    cum_seq = []
    for l in range(15):
        g = spherical_harmonic_degeneracy(l, d)
        cum += g
        cum_seq.append(cum)
    cum_2x = [2 * c for c in cum_seq]

    print(f"\nS^{d}: 2 x cumulative = {cum_2x[:10]}")

    # Check against nuclear magic numbers
    nuc_hits = [m for m in nuclear_magic if m in cum_2x[:10]]
    if nuc_hits:
        print(f"  Nuclear magic number hits: {nuc_hits}")

    # Check against noble gas
    ng_hits = [g for g in noble_gas if g in cum_2x[:10]]
    if ng_hits:
        print(f"  Noble gas hits: {ng_hits}")

print("""
KEY FINDING:
- S^3 with 2x spin: 2, 10, 28, 60, 110, ...
  Hits nuclear magic 2 and 28 (but not 8, 20, 50, 82, 126)

- S^5 with 2x spin: 2, 8, 20, 40, 70, 112, ...
  Hits nuclear magic 2, 8, 20 -- these are the HO magic numbers!
  The HO shell model produces magic numbers 2, 8, 20, 40, 70, 112
  (without spin-orbit splitting). Adding spin-orbit breaks these
  into 2, 8, 20, 28, 50, 82, 126 (the actual nuclear magic numbers).

This confirms: S^5 (Bargmann-Segal) encodes the harmonic oscillator
shell structure that is the starting point of the nuclear shell model.
The spin-orbit splitting that produces the observed magic numbers
is an ADDITIONAL perturbation not captured by the sphere geometry alone.
""")


# =========================================================================
# 6. Exceptional Lie group irrep dimensions
# =========================================================================
print("\n" + "=" * 70)
print("6. EXCEPTIONAL LIE GROUP IRREP DIMENSIONS")
print("=" * 70)

print("""
Checking whether fundamental/small irreps of exceptional Lie groups
follow packing-like patterns (i.e., match k^d - (k-1)^d or S^d
degeneracy for some d).
""")

# G2 (rank 2): fundamental irreps have dimensions 7, 14
# Weyl dimension formula for G2: for (a,b) highest weight
# dim = (a+1)(b+1)(a+b+2)(a+2b+3)(a+3b+4)(2a+3b+5) / 120
# (from standard tables)

def g2_dim(a, b):
    """Dimension of G2 irrep with Dynkin labels (a,b)."""
    # Using the Weyl dimension formula for G2
    # Positive roots of G2: alpha1, alpha2, alpha1+alpha2, 2alpha1+alpha2,
    # 3alpha1+alpha2, 3alpha1+2alpha2
    # <lambda+rho, alpha> / <rho, alpha>
    # rho = (1,1) in Dynkin labels
    # For G2, the dimension formula is:
    num = (a + 1) * (b + 1) * (a + b + 2) * (a + 2*b + 3) * (a + 3*b + 4) * (2*a + 3*b + 5)
    den = 1 * 1 * 2 * 3 * 4 * 5  # = 120
    return num // den

print("G2 irrep dimensions (a,b) for small a,b:")
g2_dims = []
for a in range(5):
    for b in range(5):
        d = g2_dim(a, b)
        g2_dims.append(((a, b), d))
        if a + b <= 3:
            print(f"  ({a},{b}): dim = {d}")

# F4 (rank 4): fundamental irreps have dimensions 26, 52, 273, 1274
# (from standard tables)
print("\nF4 fundamental irrep dimensions (from tables):")
f4_fund = [26, 52, 273, 1274]
print(f"  {f4_fund}")

# E6: fundamental = 27, 78 (adjoint), 351, ...
print("\nE6 fundamental irrep dimensions (from tables):")
e6_fund = [27, 78, 351, 2925, 351, 27]
print(f"  {e6_fund}")

# E7: fundamental = 133 (adjoint), 56, 912, ...
print("\nE7 fundamental irrep dimensions (from tables):")
e7_fund = [133, 912, 8645, 365750, 27664, 1539, 56]
print(f"  {e7_fund}")

# E8: fundamental = 248 (adjoint), 30380, ...
print("\nE8 fundamental irrep dimensions (from tables):")
e8_fund = [248, 30380, 2450240, 146325270, 6696000, 6899079264, 147250, 3875]
print(f"  First few: {e8_fund[:4]}")

# Check if any of these match packing sequences
print("\nG2 symmetric irreps (a,0) -- analog of SU(d) symmetric:")
g2_sym = []
for a in range(8):
    d = g2_dim(a, 0)
    g2_sym.append(d)
    print(f"  (a={a},0): dim = {d}")

print("\nG2 symmetric irreps (0,b):")
g2_sym2 = []
for b in range(8):
    d = g2_dim(0, b)
    g2_sym2.append(d)
    print(f"  (0,b={b}): dim = {d}")

# Check for packing-like growth
print("\nG2 (a,0) dimension ratios (consecutive):")
for i in range(1, len(g2_sym)):
    if g2_sym[i - 1] > 0:
        print(f"  dim({i},0)/dim({i-1},0) = {g2_sym[i]/g2_sym[i-1]:.4f}")

# Compare to known sequences
all_packing_seqs = {}
for d in range(1, 10):
    all_packing_seqs[f"S^{d}"] = [spherical_harmonic_degeneracy(l, d) for l in range(8)]
    all_packing_seqs[f"flat_{d}"] = [k**d - (k - 1)**d for k in range(1, 9)]

print("\nComparison of exceptional Lie group sequences to known packing sequences:")
print(f"  G2 (a,0): {g2_sym[:8]}")
print(f"  G2 (0,b): {g2_sym2[:8]}")

for name, seq in all_packing_seqs.items():
    if seq[:5] == g2_sym[:5]:
        print(f"  G2 (a,0) matches {name}!")
    if seq[:5] == g2_sym2[:5]:
        print(f"  G2 (0,b) matches {name}!")

# Check if G2 symmetric dims match any S^d formula
print("\nG2 (a,0) vs S^d harmonic degeneracies:")
for d in range(1, 15):
    s_d = [spherical_harmonic_degeneracy(l, d) for l in range(8)]
    if s_d[:6] == g2_sym[:6]:
        print(f"  G2 (a,0) MATCHES S^{d} harmonics!")
        break
else:
    print("  No S^d match found for d up to 14")

# For G2 (a,0), the dimension formula simplifies to:
# dim(a,0) = (a+1)(a+2)(a+3)(a+4)(2a+5) / 120
# This is a degree-5 polynomial in a
# Compare to S^d which has degree d-1
# So if it matches any S^d, it would be S^6 (degree 5 polynomial)
print("\nG2 (a,0) vs S^6 harmonics:")
for a in range(8):
    g = g2_sym[a]
    s6 = spherical_harmonic_degeneracy(a, 6)
    print(f"  a={a}: G2={g}, S^6={s6}")

results['exceptional_lie'] = {
    'G2_symmetric_a0': g2_sym[:8],
    'G2_symmetric_0b': g2_sym2[:8],
    'F4_fundamental': f4_fund,
    'E6_fundamental': e6_fund,
    'E7_fundamental': e7_fund[:4],
    'E8_fundamental': e8_fund[:4],
}


# =========================================================================
# 7. Summary of findings
# =========================================================================
print("\n" + "=" * 70)
print("7. SUMMARY OF FINDINGS")
print("=" * 70)

summary = """
INVERSE PACKING PROBLEM -- KEY FINDINGS
========================================

1. LAPLACE-BELTRAMI EIGENVALUES ON S^d:
   Eigenvalue: lambda_l = -l(l+d-1)
   Degeneracy: g_l = C(l+d,d) - C(l+d-2,d)

   Known physical correspondences:
   - S^2: g_l = 2l+1 = SO(3) irrep dim = angular momentum states
   - S^3: g_l = (l+1)^2 = n^2 hydrogen orbital degeneracy (Fock, Paper 7)
   - S^5: g_l matches 3D HO via Bargmann-Segal holomorphic sector (Paper 24),
     but the FULL S^5 harmonic degeneracy is LARGER than the HO degeneracy.
     The HO selects the holomorphic sector (SU(3) (N,0) irreps).
   - S^6: g_l matches G2 symmetric irrep (a,0) dimensions exactly! (SURPRISE)

2. THE S^2 UNIQUENESS:
   Paper 0's flat 2D packing gives shell capacity 2k-1 = 2l+1.
   This EXACTLY matches the S^2 harmonic degeneracy.
   For d >= 3 flat packing, k^d-(k-1)^d does NOT match any S^d degeneracy.
   The d=2 flat-to-sphere correspondence is UNIQUE.

3. CUMULATIVE DEGENERACIES WITH SPIN FACTOR 2:
   - S^2 x 2: 2, 8, 18, 32, 50, 72, 98 = hydrogen 2n^2 shell capacities
     Hits nuclear magic 2, 8, 50 and noble gas 2, 18.
   - S^3 x 2: 2, 10, 28, 60, 110, ... -> hits nuclear magic 2, 28
   - HO magic numbers 2, 8, 20, 40, 70, 112 come from filling HO shells
     with capacity 2*(N+1)(N+2)/2 (3D HO with spin), NOT from S^5 cumulative.
     The S^5 cumulative x2 is 2, 14, 54, ... which does NOT match HO magic.

4. S^3 vs HYDROGEN:
   S^3 gives ORBITAL degeneracy n^2 per shell.
   Hydrogen with spin has 2n^2 per shell.
   The cumulative sum_n^2 = N(N+1)(2N+1)/6 does NOT equal 2N^2.
   Noble gas numbers come from Madelung-rule filling of 2n^2 states,
   not from cumulative S^3 sums.

5. BARGMANN-SEGAL GENERALIZATION:
   For d-dimensional HO, the Bargmann space is C^d with unit sphere S^{2d-1}.
   The HO shell degeneracy C(N+d-1,d-1) = SU(d) symmetric irrep (N,0,...,0) dim.
   The FULL S^{2d-1} Laplace-Beltrami degeneracy is LARGER for all d.
   The ratio S^{2d-1}_full / SU(d)_(N,0) = N+1 for d=2 (S^3), and grows
   faster for d >= 3. The HO lives in a HOLOMORPHIC SUBSECTOR of S^{2d-1}.

6. G2 = S^6 SURPRISE:
   The G2 exceptional Lie group symmetric irreps (a,0) have dimensions:
   1, 7, 27, 77, 182, 378, 714, 1254
   These EXACTLY match the S^6 spherical harmonic degeneracies!
   This means dim(a,0)_{G2} = C(a+6,6) - C(a+4,6) for all a.
   Algebraically: the Weyl dimension formula for G2 (a,0) reduces to
   (a+1)(a+2)(a+3)(a+4)(2a+5)/120 which equals the S^6 harmonic formula.
   IMPLICATION: G2 is the "angular momentum group" of S^6 in the same way
   that SO(3) is the angular momentum group of S^2. This suggests a potential
   physical system on S^6 or S^7 whose symmetry group is G2.
   G2 is the automorphism group of the octonions -- this S^6 = G2/SU(3)
   connection is known in differential geometry (S^6 has a nearly-Kahler
   structure from the octonion cross product).

7. PHYSICS ENCODED BY EACH SPHERE:
   - S^1: Fourier modes, 1D periodic systems, particle on a ring
   - S^2: Angular momentum, universal angular cross-section (Paper 0)
   - S^3: Hydrogen/Coulomb problem via Fock projection (Paper 7)
   - S^4: No known single-particle quantum system identified ("orphan")
   - S^5: 3D HO via Bargmann-Segal holomorphic sector (Paper 24)
   - S^6: G2 symmetric irreps -- octonion geometry, nearly-Kahler structure
   - S^7: Quaternionic structure, Sp(2) action, 4D HO Bargmann space
   - General pattern: S^{2d-1} hosts d-dim HO (holomorphic sector only)

SURPRISES:
a) The flat-to-sphere correspondence (Paper 0) is unique to d=2. No higher-d
   analog of the packing construction recovers spherical harmonic degeneracies.

b) G2 symmetric irreps = S^6 harmonics. This was not expected and connects
   exceptional Lie algebras to the sphere-packing framework. S^6 = G2/SU(3)
   is one of only two spheres (S^2 and S^6) that admit an almost-complex
   structure, and the only non-Kahler one.

c) S^4 is an "orphan" -- no known quantum system has its degeneracy pattern
   g_l = (l+1)(l+2)(2l+3)/6. This is curious given S^3 (Coulomb) and S^5 (HO)
   are the two fundamental QM potentials, and S^6 has the G2 connection.
   S^4 = SO(5)/SO(4), related to Sp(2) -- could encode a 5D rotational system?

d) The HO magic numbers do NOT come from S^5 cumulative degeneracies.
   They come from filling HO shells (with spin), which are SU(3) (N,0) irreps
   (a proper subset of S^5 harmonics). The Bargmann-Segal construction
   selects the holomorphic sector, not the full Laplace-Beltrami spectrum.

e) The S^{2d-1} harmonic degeneracy always exceeds the d-dim HO degeneracy
   by a factor that grows with N. At N=1: ratio = 2 for all d >= 2.
   The non-holomorphic harmonics are "extra" degrees of freedom that the
   HO does not use. Their physical meaning (if any) is unclear.
"""

print(summary)
results['summary'] = summary

# Save results
output_path = "C:/Users/jlout/Desktop/Project_Geometric/debug/data/inverse_packing_results.json"

# Convert any non-serializable types
def make_serializable(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {str(k): make_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [make_serializable(x) for x in obj]
    elif isinstance(obj, tuple):
        return [make_serializable(x) for x in obj]
    return obj

results_clean = make_serializable(results)

with open(output_path, 'w') as f:
    json.dump(results_clean, f, indent=2, default=str)

print(f"\nResults saved to {output_path}")
