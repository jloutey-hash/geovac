#!/usr/bin/env python3
"""
Phase 2: Combinatorial Numerology Test for Paper 2
===================================================
Answers: "Given natural spectral ingredients of S1, S2, S3, how likely
is it to hit alpha^-1 ~ 137.036 to 8.8e-8 precision by chance?"

Search space: formulas of the form  alpha^-1 = T * (A + B + C)
where T is a transcendental prefactor and A, B, C are ratios of
spectral integers or transcendental constants. Each formula is
tested both as a linear prediction (alpha^-1 = T*S) and as the
bare coupling in a cubic self-consistency equation (alpha^3 - K*alpha + 1 = 0).
"""

import numpy as np
import bisect
import time
import os
import sys
from math import gcd, comb

# ================================================================
# Constants
# ================================================================
CODATA = 137.035999084  # CODATA 2018 alpha^-1
alpha_codata = 1.0 / CODATA
K_EXACT = 1.0 / alpha_codata + alpha_codata**2  # K that gives CODATA via cubic

output_dir = os.path.dirname(os.path.abspath(__file__))

# ================================================================
# STEP 1: Build ingredient pools
# ================================================================
print("=" * 72)
print("STEP 1: Building ingredient pools")
print("=" * 72)

# --- Base integer pool (task specification) ---
base_integers = set()

# S3 eigenvalues |lambda_n| for n=1..4: {0, 3, 8, 15}
for n in range(1, 5):
    base_integers.add(n * n - 1)

# S3 degeneracies g_n for n=1..4: {1, 4, 9, 16}
for n in range(1, 5):
    base_integers.add(n * n)

# S2 Casimir eigenvalues l(l+1) for l=0..3: {0, 2, 6, 12}
for l in range(4):
    base_integers.add(l * (l + 1))

# S2 degeneracies 2l+1 for l=0..3: {1, 3, 5, 7}
for l in range(4):
    base_integers.add(2 * l + 1)

# Cumulative state counts N(n) for n=1..4: {1, 5, 14, 30}
c = 0
for n in range(1, 5):
    c += n * n
    base_integers.add(c)

# Small integers: {1, 2, 3, 4}
for k in range(1, 5):
    base_integers.add(k)

# Special Casimir values
base_integers.add(42)  # DW Casimir sum at nmax=3
base_integers.add(10)  # Unweighted Casimir sum at nmax=3

print(f"Base integers ({len(base_integers)}): {sorted(base_integers)}")

# --- Extended pool: add products p*q <= 50 ---
# This captures the paper's boundary term 8*5=40, as well as
# alternative constructions like 8*4=32. Capped at 50 to keep
# the term pool manageable while covering all relevant products.
extended = set(base_integers)
base_nz = sorted(x for x in base_integers if x > 0)
products_added = set()
for p in base_nz:
    for q in base_nz:
        prod = p * q
        if prod <= 50 and prod not in base_integers:
            extended.add(prod)
            products_added.add(prod)

integers = sorted(extended)
integers_nz = sorted(x for x in integers if x > 0)

print(f"Products added ({len(products_added)}): {sorted(products_added)}")
print(f"Extended pool ({len(integers)}): {integers}")

# --- Transcendental ingredients ---
zeta2 = np.pi**2 / 6       # 1.6449340668...
zeta3 = 1.2020569031595942  # Apery's constant
zeta4 = np.pi**4 / 90       # 1.0823232337...
ln2 = np.log(2)             # 0.6931471806...

transcendentals = [
    ('pi',      np.pi),
    ('pi^2',    np.pi**2),
    ('zeta(2)', zeta2),
    ('zeta(3)', zeta3),
    ('zeta(4)', zeta4),
    ('ln2',     ln2),
]
print(f"Transcendentals ({len(transcendentals)}): {[t[0] for t in transcendentals]}")

# --- Build unique term values with labels ---
# Terms = {0} ∪ {±p/q : p,q ∈ pool} ∪ {±transcendental}
term_dict = {}  # round(value, 12) -> (value, label)


def add_term(val, label):
    key = round(val, 12)
    if key not in term_dict or len(label) < len(term_dict[key][1]):
        term_dict[key] = (val, label)


# Zero
add_term(0.0, '0')

# All ratios p/q from the extended integer pool
for p in integers:
    for q in integers_nz:
        val = p / q
        g = gcd(p, q)
        pn, qn = p // g, q // g
        lbl = str(pn) if qn == 1 else f"{pn}/{qn}"
        add_term(val, lbl)
        if val != 0:
            add_term(-val, f"-{lbl}")

# Transcendentals (positive and negative)
for name, val in transcendentals:
    add_term(val, name)
    add_term(-val, f"-{name}")

# Sort by value
sorted_items = sorted(term_dict.values(), key=lambda x: x[0])
terms = np.array([x[0] for x in sorted_items])
labels = [x[1] for x in sorted_items]
N = len(terms)

print(f"\nTerm pool: {N} unique values in [{terms[0]:.4f}, {terms[-1]:.4f}]")

# Verify paper's ingredients are present
paper_terms = {'42': 42.0, 'zeta(2)': zeta2, '-1/40': -1/40}
for name, val in paper_terms.items():
    key = round(val, 12)
    found = key in term_dict
    print(f"  Paper term {name:>10s} = {val:>12.6f}: {'FOUND' if found else 'MISSING'}"
          f" as '{term_dict[key][1]}'" if found else "")

# ================================================================
# STEP 2: Define formula space
# ================================================================
print("\n" + "=" * 72)
print("STEP 2: Formula space")
print("=" * 72)

prefactors = [
    (np.pi,           'pi'),
    (2 * np.pi,       '2pi'),
    (np.pi**2,        'pi^2'),
    (np.pi / 2,       'pi/2'),
    (1 / (2 * np.pi), '1/(2pi)'),
    (4 * np.pi,       '4pi'),
    (4 * np.pi**2,    '4pi^2'),
    (1.0,             '1'),
    (2.0,             '2'),
    (0.5,             '1/2'),
    (np.pi / 4,       'pi/4'),
    (np.pi / 6,       'pi/6'),
    (1 / np.pi,       '1/pi'),
    (np.sqrt(np.pi),  'sqrt(pi)'),
]

N_triples = comb(N + 2, 3)
N_per_mode = len(prefactors) * N_triples
N_total = 2 * N_per_mode  # linear + cubic

print(f"Terms: {N}")
print(f"Unordered triples (with replacement): C({N}+2, 3) = {N_triples:,}")
print(f"Prefactors: {len(prefactors)}")
print(f"Formulas per mode: {N_per_mode:,}")
print(f"Total (linear + cubic): {N_total:,}")

# Feasibility check
min_sum = 3 * terms[0]
max_sum = 3 * terms[-1]
print(f"\nAchievable sum range: [{min_sum:.1f}, {max_sum:.1f}]")
print(f"Targets by prefactor:")
for T_val, T_label in prefactors:
    tgt = CODATA / T_val
    reachable = min_sum <= tgt <= max_sum
    print(f"  {T_label:>10s}: target_sum = {tgt:10.3f}  "
          f"{'reachable' if reachable else 'OUT OF RANGE'}")

# ================================================================
# STEP 3: Search for hits
# ================================================================
print("\n" + "=" * 72)
print("STEP 3: Searching (this may take a few minutes)...")
print("=" * 72)

thresholds = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 8.8e-8]
PAPER_ERR = 8.8e-8

hits_lin = {t: 0 for t in thresholds}
hits_cub = {t: 0 for t in thresholds}

TOP_CAP = 50
top_matches = []  # (rel_err, formula_str, alpha_inv, mode)


def clean_formula(T_label, idx_i, idx_j, idx_k):
    """Build human-readable formula string."""
    parts = []
    for idx in [idx_i, idx_j, idx_k]:
        lbl = labels[idx]
        if lbl != '0':
            parts.append(lbl)
    if not parts:
        return f"{T_label} * 0"

    expr = parts[0]
    for p in parts[1:]:
        if p.startswith('-'):
            expr += f" - {p[1:]}"
        else:
            expr += f" + {p}"

    if T_label == '1':
        return expr
    return f"{T_label} * ({expr})"


def try_update_top(rel_err, idx_i, idx_j, idx_k, T_label, alpha_inv, mode):
    """Insert into top-N list if good enough."""
    global top_matches
    if len(top_matches) >= TOP_CAP and rel_err >= top_matches[-1][0]:
        return
    formula = clean_formula(T_label, idx_i, idx_j, idx_k)
    if mode == 'cubic':
        formula = f"[cubic] {formula}"
    top_matches.append((rel_err, formula, alpha_inv, mode))
    top_matches.sort(key=lambda x: x[0])
    if len(top_matches) > TOP_CAP:
        top_matches = top_matches[:TOP_CAP]


t_start = time.time()

for pf_idx, (T_val, T_label) in enumerate(prefactors):
    target_lin = CODATA / T_val
    target_cub = K_EXACT / T_val

    # Skip prefactors where target is unreachable
    min_s = 3 * terms[0]
    max_s = 3 * terms[-1]
    margin = thresholds[0] * CODATA / abs(T_val)

    if target_lin < min_s - margin and target_cub < min_s - margin:
        elapsed = time.time() - t_start
        print(f"  [{pf_idx+1:2d}/{len(prefactors)}] {T_label:>10s}  "
              f"SKIPPED (target {target_lin:.1f} below range)  {elapsed:.1f}s")
        continue
    if target_lin > max_s + margin and target_cub > max_s + margin:
        elapsed = time.time() - t_start
        print(f"  [{pf_idx+1:2d}/{len(prefactors)}] {T_label:>10s}  "
              f"SKIPPED (target {target_lin:.1f} above range)  {elapsed:.1f}s")
        continue

    # 1% tolerance in sum-space
    max_tol = thresholds[0] * CODATA / abs(T_val)

    lin_hits_this = 0
    cub_hits_this = 0

    for i in range(N):
        for j in range(i, N):
            sij = terms[i] + terms[j]

            # --- LINEAR ---
            needed = target_lin - sij
            lo = bisect.bisect_left(terms, needed - max_tol, j, N)
            hi = bisect.bisect_right(terms, needed + max_tol, j, N)

            for k in range(lo, hi):
                ainv = T_val * (sij + terms[k])
                re = abs(ainv - CODATA) / CODATA
                for t in thresholds:
                    if re < t:
                        hits_lin[t] += 1
                    else:
                        break
                if re < thresholds[0]:
                    lin_hits_this += 1
                try_update_top(re, i, j, k, T_label, ainv, 'linear')

            # --- CUBIC ---
            needed = target_cub - sij
            lo = bisect.bisect_left(terms, needed - max_tol, j, N)
            hi = bisect.bisect_right(terms, needed + max_tol, j, N)

            for k in range(lo, hi):
                K_val = T_val * (sij + terms[k])
                if K_val < 3:
                    continue
                # Approximation: alpha^-1 ≈ K - 1/K^2
                # Accurate to O(1/K^5) ~ 4e-11, well below finest threshold
                ainv = K_val - 1.0 / (K_val * K_val)
                re = abs(ainv - CODATA) / CODATA
                for t in thresholds:
                    if re < t:
                        hits_cub[t] += 1
                    else:
                        break
                if re < thresholds[0]:
                    cub_hits_this += 1
                try_update_top(re, i, j, k, T_label, ainv, 'cubic')

    elapsed = time.time() - t_start
    print(f"  [{pf_idx+1:2d}/{len(prefactors)}] {T_label:>10s}  "
          f"target={target_lin:10.3f}  hits(1%): lin={lin_hits_this:,} "
          f"cub={cub_hits_this:,}  {elapsed:.1f}s")

total_time = time.time() - t_start
print(f"\nSearch complete in {total_time:.1f}s")

# Refine top cubic matches with exact cubic solve
print("Refining top matches with exact cubic roots...")
for idx in range(len(top_matches)):
    re, formula, ainv, mode = top_matches[idx]
    if mode == 'cubic':
        # Recover K from the approximation: K ≈ ainv + 1/ainv^2
        K_approx = ainv + 1.0 / (ainv * ainv)
        roots = np.roots([1, 0, -K_approx, 1])
        real_pos = [r.real for r in roots
                    if abs(r.imag) < 1e-10 and r.real > 0]
        if real_pos:
            ainv_exact = 1.0 / min(real_pos)
            re_exact = abs(ainv_exact - CODATA) / CODATA
            top_matches[idx] = (re_exact, formula, ainv_exact, mode)

top_matches.sort(key=lambda x: x[0])

# ================================================================
# STEP 4: Results
# ================================================================
print("\n" + "=" * 72)
print("RESULTS")
print("=" * 72)

hits_total = {t: hits_lin[t] + hits_cub[t] for t in thresholds}

print(f"\n{'Threshold':>14s}  {'Linear':>10s}  {'Cubic':>10s}  "
      f"{'Total':>10s}  {'Fraction':>14s}")
print("-" * 66)
for t in thresholds:
    frac = hits_total[t] / N_total if N_total > 0 else 0
    marker = " *" if t == PAPER_ERR else ""
    print(f"  {t:>12.1e}  {hits_lin[t]:>10,}  {hits_cub[t]:>10,}  "
          f"{hits_total[t]:>10,}  {frac:>14.4e}{marker}")

# P-value
paper_hits = hits_total[PAPER_ERR]
p_value = paper_hits / N_total if N_total > 0 else 0

print(f"\n{'=' * 66}")
print("P-VALUE CALCULATION")
print(f"{'=' * 66}")
print(f"Total formulas searched:                {N_total:>14,}")
print(f"Formulas matching paper precision (8.8e-8): {paper_hits:>10,}")
if paper_hits > 0:
    print(f"p-value = {paper_hits} / {N_total:,} = {p_value:.4e}")
else:
    upper = 1.0 / N_total
    print(f"p-value < {upper:.2e} (zero hits; upper bound)")
    p_value = upper  # for assessment

# ================================================================
# TOP 20
# ================================================================
print(f"\n{'=' * 72}")
print("TOP 20 BEST MATCHES")
print(f"{'=' * 72}")
print(f"{'#':>3}  {'|Rel Err|':>12}  {'alpha^-1':>14}  Formula")
print("-" * 95)
for rank, (re, formula, ainv, mode) in enumerate(top_matches[:20], 1):
    # Flag the paper's formula
    is_paper = ('42' in formula and 'zeta(2)' in formula and '1/40' in formula
                and 'pi' in formula and 'cubic' in formula)
    marker = "  <-- PAPER" if is_paper else ""
    print(f" {rank:2d}  {re:12.2e}  {ainv:14.8f}  {formula}{marker}")

# Show more if interesting
if len(top_matches) > 20:
    print(f"\n  (... {len(top_matches) - 20} more in saved results file)")

# ================================================================
# STEP 5: Interpretation
# ================================================================
print(f"\n{'=' * 72}")
print("INTERPRETATION")
print(f"{'=' * 72}")

# Assessment
if p_value < 1e-4:
    verdict = "VERY UNLIKELY TO BE ACCIDENTAL"
    detail = (f"Only {paper_hits} out of {N_total:,} formulas ({p_value:.2e}) achieve "
              f"comparable or better precision than the paper's 8.8e-8.")
elif p_value < 1e-2:
    verdict = "SUGGESTIVE BUT NOT CONCLUSIVE"
    detail = (f"{paper_hits} out of {N_total:,} formulas ({p_value:.2e}) achieve "
              f"comparable precision. The result is notable but not uniquely precise.")
else:
    verdict = "LIKELY NUMEROLOGICAL"
    detail = (f"{paper_hits} out of {N_total:,} formulas ({p_value:.2e}) match "
              f"as well or better. Many formulas achieve comparable precision.")

print(f"\n  Verdict: {verdict}")
print(f"  {detail}")

# Context and analysis
print(f"\n  Context:")
print(f"  - The formula space uses {N} term values from spectral quantities")
print(f"    of S1, S2, S3 (eigenvalues, degeneracies, Casimir traces)")
print(f"  - Combined with {len(prefactors)} transcendental prefactors")
print(f"  - Both linear (alpha^-1 = T*S) and cubic (alpha^3 - T*S*alpha + 1 = 0)")
print(f"    interpretations tested: {N_total:,} total formulas")
print(f"  - The cubic correction shifts alpha^-1 by ~{K_EXACT - CODATA:.1e},")
print(f"    providing ~5x precision improvement over the linear interpretation")

# Competitor analysis
n_cubic_top = sum(1 for _, _, _, m in top_matches[:20] if m == 'cubic')
n_linear_top = sum(1 for _, _, _, m in top_matches[:20] if m == 'linear')
print(f"\n  Top-20 composition: {n_cubic_top} cubic, {n_linear_top} linear")

# How many distinct prefactors appear in top 20?
top20_prefactors = set()
for _, formula, _, _ in top_matches[:20]:
    for _, T_label in prefactors:
        if formula.startswith(f"[cubic] {T_label} ") or \
           formula.startswith(f"{T_label} ") or \
           (T_label == '1' and not any(formula.startswith(f"{t[1]} ") for t in prefactors if t[1] != '1')):
            top20_prefactors.add(T_label)
            break
print(f"  Distinct prefactors in top-20: {top20_prefactors}")

# Is paper's formula in top results?
paper_found = False
for rank, (re, formula, ainv, mode) in enumerate(top_matches, 1):
    if '42' in formula and 'zeta(2)' in formula and '1/40' in formula:
        paper_found = True
        print(f"\n  Paper's formula found at rank {rank} with error {re:.2e}")
        break
if not paper_found:
    print(f"\n  Paper's specific formula (pi*(42 + zeta(2) - 1/40)) was NOT found")
    print(f"  in the top-{TOP_CAP} results.")

# Scaling analysis: how does hit rate scale with threshold?
print(f"\n  Scaling analysis (log-log slope of hits vs threshold):")
import warnings
warnings.filterwarnings('ignore')
log_thresh = [np.log10(t) for t in thresholds if hits_total[t] > 0]
log_hits = [np.log10(hits_total[t]) for t in thresholds if hits_total[t] > 0]
if len(log_thresh) >= 2:
    # Fit line to log-log
    coeffs = np.polyfit(log_thresh, log_hits, 1)
    print(f"  hits ~ threshold^{coeffs[0]:.2f}")
    if abs(coeffs[0] - 1.0) < 0.3:
        print(f"  (Slope ~1: consistent with uniform distribution in alpha^-1 space)")
    else:
        print(f"  (Slope differs from 1: non-trivial clustering near CODATA)")

# ================================================================
# Save results
# ================================================================
results_path = os.path.join(output_dir, 'phase2_results.txt')
with open(results_path, 'w') as f:
    f.write("Phase 2: Combinatorial Numerology Test — Results\n")
    f.write("=" * 72 + "\n\n")

    f.write(f"CODATA alpha^-1 = {CODATA}\n")
    f.write(f"K_exact (cubic) = {K_EXACT:.10f}\n")
    f.write(f"Base integers: {sorted(base_integers)}\n")
    f.write(f"Extended (products <= 50): {integers}\n")
    f.write(f"Transcendentals: {[t[0] for t in transcendentals]}\n")
    f.write(f"Prefactors: {[p[1] for p in prefactors]}\n")
    f.write(f"Unique terms: {N}\n")
    f.write(f"Total formulas: {N_total:,}\n\n")

    f.write("HIT COUNTS\n" + "-" * 66 + "\n")
    f.write(f"{'Threshold':>14s}  {'Linear':>10s}  {'Cubic':>10s}  "
            f"{'Total':>10s}  {'Fraction':>14s}\n")
    for t in thresholds:
        frac = hits_total[t] / N_total
        f.write(f"  {t:>12.1e}  {hits_lin[t]:>10,}  {hits_cub[t]:>10,}  "
                f"{hits_total[t]:>10,}  {frac:>14.4e}\n")

    f.write(f"\np-value = {p_value:.4e}")
    if paper_hits == 0:
        f.write(f" (upper bound; zero hits)\n")
    else:
        f.write(f"\n")
    f.write(f"({paper_hits} formulas of {N_total:,} match paper precision)\n\n")

    f.write(f"TOP {TOP_CAP} MATCHES\n" + "-" * 95 + "\n")
    f.write(f"{'#':>3}  {'|Rel Err|':>12}  {'alpha^-1':>14}  Formula\n")
    for rank, (re, formula, ainv, mode) in enumerate(top_matches, 1):
        f.write(f" {rank:2d}  {re:12.2e}  {ainv:14.8f}  {formula}\n")

    f.write(f"\nVERDICT: {verdict}\n")
    f.write(f"{detail}\n")

print(f"\nResults saved to: {results_path}")
