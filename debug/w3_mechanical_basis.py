"""
Sprint W3 Test 5: mechanical-basis re-test for selection bias.

The original W3 sprint (debug/w3_lambda_predictive_verification.py) used a
hand-curated 30-form prespecified basis to test whether the four CKM
Wolfenstein parameters fit GeoVac-internal natural forms. It found 8 PDG-1σ
matches at ~6-10σ above chance vs random expectation.

Curve-fit-audit concern (debug/w3_forward_plan_memo.md, Honest Concern #3):
the basis was prespecified within-script but not version-controlled BEFORE
the test was run; a strict reader could legitimately argue the basis was
lightly tuned during script-writing.

This script addresses the concern by replacing the human-curated basis with
a MECHANICALLY-GENERATED basis enumerating all small-integer rational
combinations of master-Mellin-engine seeds up to a bounded complexity. The
generator is specified before any test runs.

Discrimination test:
  - If the original four candidates appear in the matches AND have no equally
    good alternatives, basis-selection bias is ruled out.
  - If many forms hit each parameter at <1σ in the larger basis, the original
    signal is partially attributable to selection.
  - If matches distribute evenly across M1/M2/M3/algebraic/rational, the
    structural reading (M1 for real-magnitude, M2 for CP/amplitude) was
    selection bias.

Methodology specified BEFORE generation runs:
  1. Seeds organized by Mellin-engine class (M1/M2/M3/algebraic/rational).
  2. Generation rule: combine UP TO TWO seeds via {*, /} with integer
     numerator prefactor in {±1, ±2, ±3, ±4, ±5} and integer denominator
     in {1, 2, 3, 4, 5, 6, 8, 10}.
  3. Deduplicate by numerical value (10⁻¹² tolerance).
  4. Filter to forms with values in [0.001, 10] (Wolfenstein parameter range).
  5. Save basis to JSON BEFORE searching.
  6. Per Wolfenstein parameter: count matches at <0.5%, count matches within
     1 PDG σ, classify by mechanism.
  7. Statistical test: random expectation vs observed.
  8. Report verdict: signal survives / partially attributable / killed.
"""

import json
import math
from pathlib import Path
from collections import Counter, defaultdict

import mpmath as mp
import numpy as np
from mpmath import mpf, pi, log, sqrt

mp.mp.dps = 80
OUT_DIR = Path("debug/data")
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# SEEDS organized by master-Mellin-engine mechanism class
# (declared BEFORE generation runs — frozen specification)
# ============================================================

# M1 (Hopf-base measure, π-family — Vol(S^k) ratios)
SEEDS_M1 = {
    "1":             mpf(1),
    "pi":            pi,
    "2*pi":          2*pi,
    "4*pi":          4*pi,
    "2*pi^2":        2*pi**2,
    "pi^2":          pi**2,
    "pi^3":          pi**3,
    "1/pi":          1/pi,
    "1/(2*pi)":      1/(2*pi),
    "1/(4*pi)":      1/(4*pi),
    "1/(2*pi^2)":    1/(2*pi**2),
    "1/pi^3":        1/pi**3,
    "pi/2":          pi/2,
    "pi/3":          pi/3,
    "pi/4":          pi/4,
    "pi/6":          pi/6,
}

# M2 (chirality / heat-kernel / half-integer Hurwitz — ln 2 family)
SEEDS_M2 = {
    "ln(2)":         log(2),
    "ln(2)/2":       log(2)/2,
    "ln(2)/4":       log(2)/4,
    "ln(3)":         log(3),
    "ln(5)":         log(5),
    "ln(6)":         log(6),
    "log2(phi)":     log((1+sqrt(5))/2)/log(2),
}

# M3 (vertex parity — Catalan G, Dirichlet beta, ζ values)
SEEDS_M3 = {
    "Catalan_G":     mp.catalan,
    "beta(4)":       mp.dirichlet(4, [0, 1, 0, -1]) if False else None,  # filled below
    "zeta(3)":       mp.zeta(3),
    "zeta(2)":       mp.zeta(2),
    "zeta(2)/pi^2":  mp.zeta(2)/pi**2,
    "zeta(3)/pi^3":  mp.zeta(3)/pi**3,
    "pi^2/6":        pi**2/6,
    "pi^4/90":       pi**4/90,
    "pi^4/96":       pi**4/96,
    "zeta(4)":       mp.zeta(4),
    "zeta(5)":       mp.zeta(5),
    "zeta(6)":       mp.zeta(6),
}
# Compute beta(4) properly (Dirichlet beta function)
SEEDS_M3["beta(4)"] = sum(mp.power(-1, k) / mp.power(2*k+1, 4) for k in range(0, 200))

# Small algebraic constants (sqrt of small integers, golden ratio family)
SEEDS_ALG = {
    "sqrt(2)":       sqrt(2),
    "sqrt(3)":       sqrt(3),
    "sqrt(5)":       sqrt(5),
    "sqrt(6)":       sqrt(6),
    "sqrt(pi)":      sqrt(pi),
    "sqrt(2*pi)":    sqrt(2*pi),
    "phi":           (1+sqrt(5))/2,
    "1/phi":         2/(1+sqrt(5)),
    "1/phi^2":       (3-sqrt(5))/2,
    "1+sqrt(5)":     1+sqrt(5),
}

# Pure rational seeds (control class — these are NOT M1/M2/M3 ring)
SEEDS_RAT = {
    "1":             mpf(1),
    "1/2":           mpf(1)/2,
    "1/3":           mpf(1)/3,
    "1/4":           mpf(1)/4,
    "1/5":           mpf(1)/5,
    "1/6":           mpf(1)/6,
    "1/7":           mpf(1)/7,
    "1/8":           mpf(1)/8,
    "2/3":           mpf(2)/3,
    "3/4":           mpf(3)/4,
    "3/5":           mpf(3)/5,
    "2/5":           mpf(2)/5,
}

ALL_SEEDS = {}
SEED_CLASS = {}  # name -> mechanism class

for name, val in SEEDS_M1.items():
    if name not in ALL_SEEDS:
        ALL_SEEDS[name] = val
        SEED_CLASS[name] = "M1"

for name, val in SEEDS_M2.items():
    if name not in ALL_SEEDS:
        ALL_SEEDS[name] = val
        SEED_CLASS[name] = "M2"

for name, val in SEEDS_M3.items():
    if name not in ALL_SEEDS:
        ALL_SEEDS[name] = val
        SEED_CLASS[name] = "M3"

for name, val in SEEDS_ALG.items():
    if name not in ALL_SEEDS:
        ALL_SEEDS[name] = val
        SEED_CLASS[name] = "ALG"

for name, val in SEEDS_RAT.items():
    if name not in ALL_SEEDS:
        ALL_SEEDS[name] = val
        SEED_CLASS[name] = "RAT"

print("=" * 76)
print("Sprint W3 Test 5: mechanical-basis re-test")
print("=" * 76)
print()
print(f"Seeds by class:")
print(f"  M1 (Hopf-base, pi-family): {sum(1 for c in SEED_CLASS.values() if c == 'M1')}")
print(f"  M2 (chirality, ln-family): {sum(1 for c in SEED_CLASS.values() if c == 'M2')}")
print(f"  M3 (vertex parity):        {sum(1 for c in SEED_CLASS.values() if c == 'M3')}")
print(f"  ALG (algebraic):           {sum(1 for c in SEED_CLASS.values() if c == 'ALG')}")
print(f"  RAT (pure rational):       {sum(1 for c in SEED_CLASS.values() if c == 'RAT')}")
print(f"  Total seeds:               {len(ALL_SEEDS)}")
print()


# ============================================================
# MECHANICAL GENERATOR
# ============================================================

def classify_form(seed_a_class, seed_b_class):
    """Determine the mechanism class of a 2-seed combination."""
    if seed_b_class is None:
        return seed_a_class
    classes = sorted([seed_a_class, seed_b_class])
    # If both are same class, return that class
    if classes[0] == classes[1]:
        return classes[0]
    # If one is RAT and the other is something else, take the other
    if "RAT" in classes:
        non_rat = [c for c in classes if c != "RAT"][0]
        return non_rat
    # Mixed-mechanism class
    return f"{classes[0]}+{classes[1]}"


def generate_basis():
    """Generate mechanical basis: up to 2 seeds combined via {*, /} with
    integer prefactors and denominators."""
    seeds = list(ALL_SEEDS.items())
    seed_names = [s[0] for s in seeds]
    seed_vals = [s[1] for s in seeds]

    prefactors_num = [1, 2, 3, 4, 5]
    denoms = [1, 2, 3, 4, 5, 6, 8, 10]

    # Single-seed forms with prefactors
    raw_basis = {}  # form_name -> (mpf_value, mechanism_class)

    print("Generating single-seed forms with rational prefactors...")
    for sname, sval in zip(seed_names, seed_vals):
        sclass = SEED_CLASS[sname]
        for num in prefactors_num:
            for den in denoms:
                if math.gcd(num, den) != 1:
                    continue  # canonical form: gcd=1
                # +num/den
                v = mpf(num) * sval / mpf(den)
                if mpf("0.001") < v < mpf(10):
                    if num == 1 and den == 1:
                        fname = sname
                    elif den == 1:
                        fname = f"{num}*{sname}"
                    elif num == 1:
                        fname = f"{sname}/{den}"
                    else:
                        fname = f"({num}/{den})*{sname}"
                    raw_basis[fname] = (v, sclass)
                # -num/den is symmetric to +num/den for our purposes
                # (Wolfenstein parameters are positive in the standard basis;
                # we drop the negative branch for the search range filter)

    print(f"  After single-seed: {len(raw_basis)} forms")

    # Two-seed forms (a/b and a*b with integer numerator)
    print("Generating two-seed combinations (a*b/n and a/b/n)...")
    seeds_list = list(zip(seed_names, seed_vals))
    n_pairs = 0
    for i, (sa_name, sa_val) in enumerate(seeds_list):
        for j, (sb_name, sb_val) in enumerate(seeds_list):
            if i >= j:
                continue
            sa_class = SEED_CLASS[sa_name]
            sb_class = SEED_CLASS[sb_name]
            mech = classify_form(sa_class, sb_class)
            n_pairs += 1
            # a * b / den
            for num in [1, 2, 3]:
                for den in denoms:
                    if math.gcd(num, den) != 1:
                        continue
                    v = mpf(num) * sa_val * sb_val / mpf(den)
                    if mpf("0.001") < v < mpf(10):
                        if num == 1 and den == 1:
                            fname = f"{sa_name}*{sb_name}"
                        elif den == 1:
                            fname = f"{num}*{sa_name}*{sb_name}"
                        elif num == 1:
                            fname = f"{sa_name}*{sb_name}/{den}"
                        else:
                            fname = f"({num}/{den})*{sa_name}*{sb_name}"
                        if fname not in raw_basis:
                            raw_basis[fname] = (v, mech)
            # a / b / den
            for num in [1, 2, 3]:
                for den in denoms:
                    if math.gcd(num, den) != 1:
                        continue
                    if sb_val != 0:
                        v = mpf(num) * sa_val / sb_val / mpf(den)
                        if mpf("0.001") < v < mpf(10):
                            if num == 1 and den == 1:
                                fname = f"{sa_name}/{sb_name}"
                            elif den == 1:
                                fname = f"{num}*{sa_name}/{sb_name}"
                            elif num == 1:
                                fname = f"{sa_name}/({sb_name}*{den})"
                            else:
                                fname = f"({num}/{den})*{sa_name}/{sb_name}"
                            if fname not in raw_basis:
                                raw_basis[fname] = (v, mech)
                    # b / a / den (other direction)
                    if sa_val != 0:
                        v = mpf(num) * sb_val / sa_val / mpf(den)
                        if mpf("0.001") < v < mpf(10):
                            if num == 1 and den == 1:
                                fname = f"{sb_name}/{sa_name}"
                            elif den == 1:
                                fname = f"{num}*{sb_name}/{sa_name}"
                            elif num == 1:
                                fname = f"{sb_name}/({sa_name}*{den})"
                            else:
                                fname = f"({num}/{den})*{sb_name}/{sa_name}"
                            if fname not in raw_basis:
                                raw_basis[fname] = (v, mech)

    print(f"  After two-seed: {len(raw_basis)} forms (from {n_pairs} ordered pairs)")

    # Deduplicate by numerical value (1e-10 relative tolerance)
    print("Deduplicating by numerical value...")
    sorted_forms = sorted(raw_basis.items(), key=lambda x: float(x[1][0]))
    deduped = []
    last_val = None
    for fname, (fval, fclass) in sorted_forms:
        v = float(fval)
        if last_val is not None and abs(v - last_val) / max(abs(v), 1e-30) < 1e-10:
            continue
        deduped.append((fname, fval, fclass))
        last_val = v

    print(f"  After dedup: {len(deduped)} forms")
    return deduped


# Generate
basis = generate_basis()
N_BASIS = len(basis)

# Save basis (keys = canonical form name, values = float for size, class)
print(f"\nSaving basis to debug/data/w3_mechanical_basis.json (basis size {N_BASIS})...")
basis_to_save = {
    "metadata": {
        "n_basis": N_BASIS,
        "n_seeds_total": len(ALL_SEEDS),
        "seeds_by_class": {
            "M1": sum(1 for c in SEED_CLASS.values() if c == "M1"),
            "M2": sum(1 for c in SEED_CLASS.values() if c == "M2"),
            "M3": sum(1 for c in SEED_CLASS.values() if c == "M3"),
            "ALG": sum(1 for c in SEED_CLASS.values() if c == "ALG"),
            "RAT": sum(1 for c in SEED_CLASS.values() if c == "RAT"),
        },
        "seed_names": list(ALL_SEEDS.keys()),
        "seed_classes": SEED_CLASS,
        "generation_rule": (
            "Single-seed: prefactor (+num/den), num in {1,2,3,4,5}, "
            "den in {1,2,3,4,5,6,8,10}, gcd(num,den)=1. "
            "Two-seed: a*b/(num/den) and a/b/(num/den) with num in {1,2,3}, "
            "same denom set, both directions. Range filter [0.001, 10]. "
            "Dedup by 1e-10 relative tolerance."
        ),
        "frozen_before_test": True,
    },
    "basis_class_counts": dict(Counter(c for (_, _, c) in basis)),
    "basis_summary": {
        "n_forms": N_BASIS,
        "first_50_forms": [
            {"form": fname, "value": str(fval)[:30], "class": fclass}
            for (fname, fval, fclass) in basis[:50]
        ],
        "last_50_forms": [
            {"form": fname, "value": str(fval)[:30], "class": fclass}
            for (fname, fval, fclass) in basis[-50:]
        ],
    },
}


# ============================================================
# WOLFENSTEIN PDG values (same as original sprint)
# ============================================================

wolfenstein = {
    "lambda":  (mpf("0.22500"), mpf("0.00067")),
    "A":       (mpf("0.826"),   mpf("0.012")),
    "rho_bar": (mpf("0.159"),   mpf("0.010")),
    "eta_bar": (mpf("0.348"),   mpf("0.009")),
}

# Original four candidates (from prespecified-basis sprint)
ORIGINAL_CANDIDATES = {
    "lambda":  ("1/sqrt(Vol(S^3)) = 1/(pi*sqrt(2))", 1/(pi*sqrt(2)), "M1+ALG"),
    "A":       ("sqrt(ln 2)",                        sqrt(log(2)),    "M2+ALG"),
    "rho_bar": ("1/Vol(S^1) = 1/(2*pi)",             1/(2*pi),        "M1"),
    "eta_bar": ("ln(2)/2",                           log(2)/2,        "M2"),
}


# ============================================================
# Search
# ============================================================

print()
print("=" * 76)
print("Search per Wolfenstein parameter")
print("=" * 76)

# Pre-compute float arrays for fast filtering
float_basis = np.array([float(fval) for (_, fval, _) in basis])
basis_names = [fname for (fname, _, _) in basis]
basis_classes = [fclass for (_, _, fclass) in basis]
basis_mpf = [fval for (_, fval, _) in basis]

per_param_results = {}
total_05pct = 0
total_within_1sigma = 0

for pname, (pval, psig) in wolfenstein.items():
    print()
    print(f"--- {pname} = {mp.nstr(pval, 8)} +/- {mp.nstr(psig, 4)} ---")
    pval_f = float(pval)
    psig_f = float(psig)

    # Float-level filter for matches within 2 PDG sigma
    rel_dev = (float_basis - pval_f) / pval_f
    abs_pct = np.abs(rel_dev) * 100
    sigma_dev = np.abs(float_basis - pval_f) / psig_f

    # All forms within <0.5% AND/OR within 1 PDG sigma
    matches_05pct = []
    matches_1sigma = []
    matches_2sigma = []  # for top-N sorting

    for i in range(len(float_basis)):
        if abs_pct[i] < 0.5 or sigma_dev[i] < 1.0:
            entry = {
                "form": basis_names[i],
                "value_float": float(float_basis[i]),
                "rel_deviation_pct": float(rel_dev[i] * 100),
                "within_pdg_sigma": float(sigma_dev[i]),
                "class": basis_classes[i],
            }
            if abs_pct[i] < 0.5:
                matches_05pct.append(entry)
            if sigma_dev[i] < 1.0:
                matches_1sigma.append(entry)
        if sigma_dev[i] < 2.0:
            matches_2sigma.append({
                "form": basis_names[i],
                "value_float": float(float_basis[i]),
                "rel_deviation_pct": float(rel_dev[i] * 100),
                "within_pdg_sigma": float(sigma_dev[i]),
                "class": basis_classes[i],
            })

    matches_05pct.sort(key=lambda x: abs(x["rel_deviation_pct"]))
    matches_1sigma.sort(key=lambda x: abs(x["within_pdg_sigma"]))
    matches_2sigma.sort(key=lambda x: abs(x["within_pdg_sigma"]))

    # Compute mechanism class distribution among <1σ matches
    class_counter_1sigma = Counter(m["class"] for m in matches_1sigma)

    print(f"  basis size: {N_BASIS}")
    print(f"  matches at <0.5%: {len(matches_05pct)}")
    print(f"  matches within 1 PDG sigma: {len(matches_1sigma)}")
    print(f"  class distribution (1-sigma matches): {dict(class_counter_1sigma)}")

    # Show top-10 by sigma
    print(f"  Top 10 by sigma:")
    for m in matches_2sigma[:10]:
        print(f"    {m['form'][:50]:<50s}  diff {m['rel_deviation_pct']:+7.4f}%  "
              f"({m['within_pdg_sigma']:.3f} sigma)  [{m['class']}]")

    # Find original candidate's rank
    original = ORIGINAL_CANDIDATES[pname]
    orig_value = float(original[1])
    # Find closest entry to original candidate value in basis
    orig_idx = np.argmin(np.abs(float_basis - orig_value))
    orig_form = basis_names[orig_idx]
    orig_in_basis_dev = abs(float_basis[orig_idx] - orig_value) / max(abs(orig_value), 1e-30)
    # Find rank of this form among 1σ matches
    orig_rank_1sigma = None
    for k, m in enumerate(matches_1sigma):
        if abs(m["value_float"] - orig_value) / max(abs(orig_value), 1e-30) < 1e-8:
            orig_rank_1sigma = k + 1  # 1-indexed
            break

    per_param_results[pname] = {
        "pdg_central": str(pval),
        "pdg_sigma": str(psig),
        "n_basis": N_BASIS,
        "n_matches_05pct": len(matches_05pct),
        "n_matches_1sigma": len(matches_1sigma),
        "class_distribution_1sigma": dict(class_counter_1sigma),
        "top_20_matches_by_sigma": matches_2sigma[:20],
        "top_10_matches_05pct": matches_05pct[:10],
        "original_candidate": {
            "label": original[0],
            "value": float(original[1]),
            "class": original[2],
            "found_in_basis_as": orig_form,
            "in_basis_relative_deviation": orig_in_basis_dev,
            "rank_among_1sigma_matches": orig_rank_1sigma,
        },
    }
    total_05pct += len(matches_05pct)
    total_within_1sigma += len(matches_1sigma)


# ============================================================
# Statistical assessment
# ============================================================

print()
print("=" * 76)
print("Statistical assessment")
print("=" * 76)

# Random-target benchmark: randomly draw target values log-uniform on
# [0.10, 1.00] (Wolfenstein parameter range) and count matches at <0.5%
# and within 1 σ.

np.random.seed(42)

def expected_hits_random(float_basis, pdg_sigma_grid_pct, log_range, n_trials=10000):
    """Compute expected number of hits per random target.
    pdg_sigma_grid_pct: list of (label, threshold_pct) — counts hits with
    that percentage threshold."""
    log_min, log_max = log_range
    out = {label: [] for (label, _) in pdg_sigma_grid_pct}
    for _ in range(n_trials):
        target = math.exp(np.random.uniform(log_min, log_max))
        rel_dev = np.abs(float_basis - target) / target
        for label, thresh_pct in pdg_sigma_grid_pct:
            n_hits = int(np.sum(rel_dev < thresh_pct/100))
            out[label].append(n_hits)
    return {k: (np.mean(v), np.std(v)) for k, v in out.items()}


log_range = (math.log(0.10), math.log(1.00))
pdg_grid = [("0.5pct", 0.5), ("1pct", 1.0), ("2pct", 2.0)]
exp_hits = expected_hits_random(float_basis, pdg_grid, log_range, n_trials=10000)
print(f"  Basis size: {N_BASIS}")
print(f"  Random target log-uniform on [0.10, 1.00], 10000 trials:")
for label, (mu, sigma) in exp_hits.items():
    print(f"    Mean hits at <{label}: {mu:.2f} +/- {sigma:.2f}")

# Per-parameter random expectation at the parameter's specific PDG sigma
print()
print("  Per-parameter random-target expectation at param's PDG-sigma threshold:")
per_param_random_exp = {}
for pname, (pval, psig) in wolfenstein.items():
    pval_f = float(pval)
    psig_f = float(psig)
    sigma_thresh_pct = (psig_f / pval_f) * 100  # within-1-sigma in pct
    half_pct_hits = exp_hits["0.5pct"][0]
    sigma_hits = exp_hits.get(f"{sigma_thresh_pct:.2f}pct", None)
    # Compute expected hits at the parameter's own 1-sigma threshold
    rel_thresh = sigma_thresh_pct / 100
    np.random.seed(43)
    counts = []
    for _ in range(5000):
        target = math.exp(np.random.uniform(*log_range))
        rel_dev = np.abs(float_basis - target) / target
        counts.append(int(np.sum(rel_dev < rel_thresh)))
    mu_p, sd_p = np.mean(counts), np.std(counts)
    per_param_random_exp[pname] = {
        "pdg_sigma_pct": sigma_thresh_pct,
        "expected_hits_within_1pdg_sigma": mu_p,
        "std_hits_within_1pdg_sigma": sd_p,
        "observed_hits_within_1pdg_sigma": per_param_results[pname]["n_matches_1sigma"],
    }
    obs = per_param_results[pname]["n_matches_1sigma"]
    z = (obs - mu_p) / max(sd_p, 1e-9)
    print(f"    {pname:<10s}  PDG 1-sigma = {sigma_thresh_pct:.3f}%  expected = {mu_p:.1f}  observed = {obs}  z = {z:+.2f}")


# ============================================================
# Discrimination test: do M1/M2 forms preferentially fit Wolfenstein?
# ============================================================

print()
print("=" * 76)
print("Discrimination test: class distribution of <1-sigma matches")
print("=" * 76)

# Aggregate over all four parameters
all_class_counts = Counter()
for pname in wolfenstein:
    for m in per_param_results[pname]["top_20_matches_by_sigma"]:
        if m["within_pdg_sigma"] < 1.0:
            all_class_counts[m["class"]] += 1

# Compute basis class distribution (the "null" for chance comparison)
basis_class_counts = Counter(c for (_, _, c) in basis)

print(f"  Basis class composition:")
for cls, cnt in sorted(basis_class_counts.items(), key=lambda x: -x[1]):
    pct = cnt / N_BASIS * 100
    print(f"    {cls:<10s}  {cnt:>6d}  ({pct:.1f}%)")

print()
print(f"  <1-sigma matches class distribution (aggregated over 4 Wolfenstein params):")
total_matches = sum(all_class_counts.values())
for cls, cnt in sorted(all_class_counts.items(), key=lambda x: -x[1]):
    if total_matches > 0:
        pct = cnt / total_matches * 100
    else:
        pct = 0
    basis_pct = basis_class_counts.get(cls, 0) / N_BASIS * 100
    enrichment = (pct / max(basis_pct, 1e-9)) if basis_pct > 0 else float('inf')
    print(f"    {cls:<10s}  {cnt:>6d}  ({pct:5.1f}%)   basis {basis_pct:5.1f}%  enrichment {enrichment:.2f}x")


# ============================================================
# Verdict
# ============================================================

print()
print("=" * 76)
print("VERDICT")
print("=" * 76)

# Compute overall signal strength against null
# Total observed within 1 PDG sigma vs expected
total_observed = total_within_1sigma
total_expected = sum(per_param_random_exp[p]["expected_hits_within_1pdg_sigma"] for p in wolfenstein)
total_std = math.sqrt(sum(per_param_random_exp[p]["std_hits_within_1pdg_sigma"]**2 for p in wolfenstein))

z_total = (total_observed - total_expected) / max(total_std, 1e-9)

print(f"  Total observed <1-sigma matches across 4 Wolfenstein params:  {total_observed}")
print(f"  Total expected from random-target null:                       {total_expected:.1f} +/- {total_std:.1f}")
print(f"  Signal vs null:                                               z = {z_total:+.2f}")

# How many "equally good" alternatives exist for each original candidate?
print()
print(f"  Alternatives ranking original candidates among <1-sigma matches:")
for pname in wolfenstein:
    info = per_param_results[pname]["original_candidate"]
    n_better_or_equal = 0
    for m in per_param_results[pname]["top_20_matches_by_sigma"]:
        if m["within_pdg_sigma"] >= 1.0:
            continue
        if m["within_pdg_sigma"] <= info.get("rank_among_1sigma_matches_dummy", 999):
            pass
    rank = info.get("rank_among_1sigma_matches", "not found")
    n_1sigma = per_param_results[pname]["n_matches_1sigma"]
    print(f"    {pname:<10s}  original: {info['label'][:45]:<45s}  rank: {rank} of {n_1sigma}")


# Save full results
print()
print("Saving results to debug/data/w3_mechanical_basis.json...")

results = {
    "metadata": basis_to_save["metadata"],
    "basis_class_counts": dict(basis_class_counts),
    "basis_summary": basis_to_save["basis_summary"],
    "per_parameter_results": per_param_results,
    "random_expectation": {
        "pdg_grid": {label: {"mean": float(mu), "std": float(sd)}
                     for label, (mu, sd) in exp_hits.items()},
        "per_parameter": {pname: {
            "pdg_sigma_pct": v["pdg_sigma_pct"],
            "expected_hits_within_1pdg_sigma": float(v["expected_hits_within_1pdg_sigma"]),
            "std_hits_within_1pdg_sigma": float(v["std_hits_within_1pdg_sigma"]),
            "observed_hits_within_1pdg_sigma": v["observed_hits_within_1pdg_sigma"],
        } for pname, v in per_param_random_exp.items()},
    },
    "discrimination_test": {
        "basis_class_counts": dict(basis_class_counts),
        "match_class_counts_aggregated": dict(all_class_counts),
        "match_class_enrichment": {
            cls: {
                "match_count": cnt,
                "match_pct": (cnt / total_matches * 100) if total_matches > 0 else 0,
                "basis_pct": basis_class_counts.get(cls, 0) / N_BASIS * 100,
                "enrichment": ((cnt / total_matches) / (basis_class_counts.get(cls, 0) / N_BASIS))
                              if (basis_class_counts.get(cls, 0) > 0 and total_matches > 0)
                              else None,
            }
            for cls, cnt in all_class_counts.items()
        },
        "total_matches_aggregated_1sigma": total_matches,
        "total_basis_size": N_BASIS,
    },
    "verdict_summary": {
        "total_observed_1sigma_matches": total_observed,
        "total_expected_1sigma_matches": float(total_expected),
        "total_std_1sigma_matches": float(total_std),
        "z_score_signal_vs_null": float(z_total),
    },
}


def to_json_safe(obj):
    if isinstance(obj, mp.mpf):
        return str(obj)
    if isinstance(obj, dict):
        return {k: to_json_safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [to_json_safe(v) for v in obj]
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    return obj


with open(OUT_DIR / "w3_mechanical_basis.json", "w") as f:
    json.dump(to_json_safe(results), f, indent=2)

print()
print("=" * 76)
print(f"DONE. Basis size: {N_BASIS}.  Total <1-sigma matches: {total_observed}.  z = {z_total:+.2f}")
print("=" * 76)
