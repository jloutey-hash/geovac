"""
Sprint W3 follow-up: predictive verification of Candidate A.

Candidate A: lambda_Wolfenstein = 1/sqrt(Vol(S^3)) = 1/(pi sqrt 2).
Single-data-point match at 0.035 percent.

Curve-fit-audit discipline says: a single coincidence-class match cannot be
elevated to "structural fact" unless the same hypothesis predicts independent
observables that we can test against PDG.

Hypothesis: the four Wolfenstein CKM parameters (lambda, A, rho_bar, eta_bar)
are each set by a clean GeoVac-internal natural form. We have lambda's
candidate; we test the other three for natural-form matches in a prespecified
basis, then check whether the combined four-candidate hypothesis predicts the
remaining CKM observables within PDG uncertainty.

The key disciplinary question: how many "natural form" fits are expected by
chance for the parameter density we have, and how many do we find?

Structural design:
  Step 1: Prespecify the GeoVac-internal natural-form basis (frozen before
          looking at the CKM data).
  Step 2: For each Wolfenstein parameter, search the basis for natural-form
          matches at <0.5 percent.
  Step 3: For matches found, report the GeoVac-internal interpretation.
  Step 4: Compute combined predictions for derived CKM observables (|V_cb|,
          |V_ub|, Cabibbo angle, Jarlskog) and compare to PDG.
  Step 5: Honest assessment of expected vs found false-positive rate.
  Step 6: Cross-sector implications: do the same natural forms predict
          anything in the lepton sector?
"""

import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from mpmath import mpf, pi, log, sqrt, exp, sin, cos, atan, atan2

mp.mp.dps = 80
OUT_DIR = Path("debug/data")

# ============================================================
# Wolfenstein parameters from PDG 2024 (central value, sigma)
# ============================================================

wolfenstein = {
    "lambda":  (mpf("0.22500"), mpf("0.00067")),   # PDG ~0.3% relative
    "A":       (mpf("0.826"),   mpf("0.012")),    # PDG ~1.5% relative
    "rho_bar": (mpf("0.159"),   mpf("0.010")),    # PDG ~6% relative
    "eta_bar": (mpf("0.348"),   mpf("0.009")),    # PDG ~3% relative
}

# ============================================================
# Step 1: PRESPECIFIED natural-form basis (frozen)
# ============================================================
#
# Each "natural form" is a closed-form expression in the GeoVac-internal
# vocabulary: framework manifold volumes, master Mellin engine constants,
# and small-integer combinations.

natural_forms = {
    # Framework manifold amplitudes -- "1/sqrt(Vol)"
    "1/sqrt(Vol(S^1)) = 1/sqrt(2 pi)":    1/sqrt(2*pi),
    "1/sqrt(Vol(S^2)) = 1/sqrt(4 pi)":    1/sqrt(4*pi),
    "1/sqrt(Vol(S^3)) = 1/(pi sqrt 2)":   1/(pi*sqrt(2)),
    "1/sqrt(Vol(S^5)) = 1/sqrt(pi^3)":    1/sqrt(pi**3),

    # Framework manifold volume reciprocals
    "1/Vol(S^1) = 1/(2 pi)":              1/(2*pi),
    "1/Vol(S^2) = 1/(4 pi)":              1/(4*pi),
    "1/Vol(S^3) = 1/(2 pi^2)":            1/(2*pi*pi),

    # Mellin engine M2 (heat-kernel / bit-vs-nat) constants and simple combos
    "ln(2)":                               log(2),
    "ln(2)/2":                             log(2)/2,
    "ln(2)/4":                             log(2)/4,
    "sqrt(ln 2)":                          sqrt(log(2)),
    "1/ln(2)":                             1/log(2),

    # Mellin engine M1 (Hopf) -- pi-type
    "pi/14":                               pi/14,
    "pi/16":                               pi/16,
    "pi/(2 phi)":                          pi/(2*(1+sqrt(5))/2),

    # GeoVac-internal Paper 2 / Paper 18 constants
    "kappa^2 = 1/256":                     mpf(1)/256,
    "Delta = 1/40":                        mpf(1)/40,
    "1/B = 1/42":                          mpf(1)/42,
    "F = pi^2/6":                          pi*pi/6,
    "alpha = 1/137.036":                   1/mpf("137.035999084"),

    # Golden / phi-family
    "1/phi = (sqrt 5 - 1)/2":              (sqrt(5)-1)/2,
    "1/phi^2 = (3-sqrt 5)/2":              (3-sqrt(5))/2,
    "phi - 1":                             (sqrt(5)-1)/2,
    "phi^2 - phi = 1":                     mpf(1),

    # Trace-natural simple combos
    "1/sqrt(20) = 1/(2 sqrt 5)":           1/sqrt(20),
    "log(2)/pi":                           log(2)/pi,
    "1/pi":                                1/pi,
    "1/(2 pi^2)":                          1/(2*pi*pi),
    "sqrt(2)/(2 pi)":                      sqrt(2)/(2*pi),
    "1/(pi sqrt(pi))":                     1/(pi*sqrt(pi)),
}


def find_matches(target_value, target_sigma, target_label,
                 forms, threshold_pct=1.0):
    """Search forms for matches within threshold_pct of target.
    Returns matches sorted by relative deviation."""
    matches = []
    for fname, fval in forms.items():
        if fval == 0 or target_value == 0:
            continue
        diff = (fval - target_value) / target_value
        within_sigma = abs(fval - target_value) / target_sigma
        if abs(float(diff)*100) < threshold_pct:
            matches.append({
                "form": fname,
                "value": str(fval),
                "rel_deviation_pct": float(diff)*100,
                "within_pdg_sigma": float(within_sigma),
            })
    matches.sort(key=lambda x: abs(x["rel_deviation_pct"]))
    return matches


# ============================================================
# Step 2: Search natural-form basis for each Wolfenstein parameter
# ============================================================

print("=" * 76)
print("Step 2: natural-form match search per Wolfenstein parameter")
print("=" * 76)

results = {"step_2_natural_form_matches": {}}
total_tests = 0
total_hits_1pct = 0
total_hits_05pct = 0
total_hits_within_1sigma = 0

for pname, (pval, psig) in wolfenstein.items():
    print()
    print(f"--- {pname} = {mp.nstr(pval, 8)} +/- {mp.nstr(psig, 4)} ---")
    matches = find_matches(pval, psig, pname, natural_forms, threshold_pct=2.0)
    n_tests = len(natural_forms)
    total_tests += n_tests
    n_1pct = sum(1 for m in matches if abs(m["rel_deviation_pct"]) < 1.0)
    n_05pct = sum(1 for m in matches if abs(m["rel_deviation_pct"]) < 0.5)
    n_1sigma = sum(1 for m in matches if abs(m["within_pdg_sigma"]) < 1.0)
    total_hits_1pct += n_1pct
    total_hits_05pct += n_05pct
    total_hits_within_1sigma += n_1sigma
    print(f"  forms tested: {n_tests}; matches at <1%: {n_1pct}; at <0.5%: {n_05pct}; "
          f"within 1 sigma PDG: {n_1sigma}")
    results["step_2_natural_form_matches"][pname] = {
        "value":  str(pval),
        "sigma":  str(psig),
        "matches": matches,
        "n_tested": n_tests,
        "n_1pct":   n_1pct,
        "n_05pct":  n_05pct,
        "n_1sigma": n_1sigma,
    }
    for m in matches[:5]:
        marker = " *" if abs(m["within_pdg_sigma"]) < 1 else ""
        print(f"     {m['form']:<40s}  diff {m['rel_deviation_pct']:+7.3f}%  "
              f"({m['within_pdg_sigma']:.2f} sigma){marker}")

print()
print(f"Aggregate: {total_tests} tests, {total_hits_1pct} hits <1%, "
      f"{total_hits_05pct} hits <0.5%, {total_hits_within_1sigma} within 1 PDG sigma.")

# ============================================================
# Step 3: Pick best match per parameter (the four-candidate hypothesis)
# ============================================================

print()
print("=" * 76)
print("Step 3: Best natural-form fit per Wolfenstein parameter (the 4-candidate hypothesis)")
print("=" * 76)

best = {}
for pname, info in results["step_2_natural_form_matches"].items():
    if info["matches"]:
        best[pname] = info["matches"][0]
    else:
        best[pname] = None

print()
for pname in ["lambda", "A", "rho_bar", "eta_bar"]:
    b = best.get(pname)
    if b:
        print(f"  {pname:<10s}  ~  {b['form']:<40s}  diff {b['rel_deviation_pct']:+7.3f}%")
    else:
        print(f"  {pname:<10s}  no <2% natural form match")

results["four_candidate_hypothesis"] = {pname: best[pname] for pname in best}

# ============================================================
# Step 4: Combined predictions for derived CKM observables
# ============================================================

print()
print("=" * 76)
print("Step 4: derived CKM observables under the 4-candidate hypothesis")
print("=" * 76)

# Use best-fit candidates (or fallback to PDG central if no fit)
def get_best_value(pname):
    b = best.get(pname)
    if b is None:
        return wolfenstein[pname][0]
    return mpf(b["value"])

L = get_best_value("lambda")
A_w = get_best_value("A")
R = get_best_value("rho_bar")
H = get_best_value("eta_bar")

print(f"  Using candidate values:")
print(f"    lambda  = {mp.nstr(L, 8)}")
print(f"    A       = {mp.nstr(A_w, 8)}")
print(f"    rho_bar = {mp.nstr(R, 8)}")
print(f"    eta_bar = {mp.nstr(H, 8)}")

# Predicted CKM matrix elements (Wolfenstein expansion at relevant order)
sin_theta_12 = L
sin_theta_23 = A_w * L * L
sin_theta_13_complex_mag = A_w * L**3 * sqrt(R*R + H*H)

theta_12_pred = mp.asin(sin_theta_12) * 180/pi
V_us_pred = L
V_cb_pred = A_w * L * L
V_ub_pred = A_w * L**3 * sqrt(R*R + H*H)
J_pred   = A_w**2 * L**6 * H
delta_CP_pred = atan2(H, R) * 180/pi

# PDG values for comparison
ckm_pdg = {
    "V_us":            (mpf("0.22500"), mpf("0.00067")),
    "V_cb":            (mpf("0.04182"), mpf("0.00085")),
    "V_ub":            (mpf("0.00369"), mpf("0.00011")),
    "Cabibbo (deg)":   (mpf("13.04"),   mpf("0.05")),
    "delta_CP (deg)":  (mpf("65.4"),    mpf("3.0")),
    "Jarlskog J":      (mpf("3.18e-5"), mpf("1.5e-7")),
}

predictions = {
    "V_us":           V_us_pred,
    "V_cb":           V_cb_pred,
    "V_ub":           V_ub_pred,
    "Cabibbo (deg)":  theta_12_pred,
    "delta_CP (deg)": delta_CP_pred,
    "Jarlskog J":     J_pred,
}

print()
print(f"  {'Observable':<18s} {'Predicted':>16s} {'PDG central':>16s} {'Diff':>10s} {'PDG sigma':>10s}")
print(f"  {'-'*72}")
predictions_table = {}
for obs in ckm_pdg:
    pred = predictions[obs]
    pdg_c, pdg_s = ckm_pdg[obs]
    diff_pct = float((pred - pdg_c)/pdg_c * 100)
    n_sigma = float(abs(pred - pdg_c)/pdg_s)
    flag = ""
    if n_sigma < 1:
        flag = " ok (within 1 sigma)"
    elif n_sigma < 2:
        flag = " marginal"
    else:
        flag = " *** OUTSIDE PDG 2 sigma ***"
    print(f"  {obs:<18s} {mp.nstr(pred,6):>16s} {mp.nstr(pdg_c,6):>16s} "
          f"{diff_pct:>+9.3f}% {n_sigma:>9.2f}{flag}")
    predictions_table[obs] = {
        "predicted": str(pred),
        "pdg_central": str(pdg_c),
        "pdg_sigma": str(pdg_s),
        "diff_pct": diff_pct,
        "n_sigma": n_sigma,
    }

results["step_4_derived_observables"] = predictions_table

# ============================================================
# Step 5: Statistical assessment of false-positive rate
# ============================================================

print()
print("=" * 76)
print("Step 5: how many natural-form matches expected by chance?")
print("=" * 76)

# For each parameter, estimate the density of natural-form values within
# 1 percent of the parameter. A 1-percent-wide band in the natural-form
# space contains about (1/100) * (range covered by basis) / (range
# covered by natural-form values).

# Crude estimator: count basis values within 1 percent of UNIFORMLY DISTRIBUTED
# random points on the same range as the Wolfenstein parameters,
# averaged over many trials.

import numpy as np
np.random.seed(11)

def expected_hits_random(forms, target_range_log, n_trials=10000, threshold_pct=1.0):
    fvals = np.array([float(f) for f in forms.values()])
    log_min, log_max = target_range_log
    hits = []
    for _ in range(n_trials):
        target = np.exp(np.random.uniform(log_min, log_max))
        n_hit = sum(1 for f in fvals if f != 0 and abs(f-target)/target < threshold_pct/100)
        hits.append(n_hit)
    return np.mean(hits), np.std(hits)

# Range covered by Wolfenstein parameters: ~[0.16, 0.83] => log range ~[-1.83, -0.19]
log_range = (math.log(0.10), math.log(1.00))
mean_hits, std_hits = expected_hits_random(natural_forms, log_range,
                                            n_trials=10000, threshold_pct=1.0)
print(f"  Random target in range [0.10, 1.00] (log-uniform):")
print(f"  Mean natural-form matches per random target at 1%: {mean_hits:.2f} +/- {std_hits:.2f}")
print(f"  Expected total over 4 Wolfenstein parameters: {4*mean_hits:.2f} +/- {math.sqrt(4)*std_hits:.2f}")
print(f"  Observed total: {total_hits_1pct} matches at <1%")

# Tighter discipline: how many at 0.5%?
mean_hits_05, std_hits_05 = expected_hits_random(natural_forms, log_range,
                                                  n_trials=10000, threshold_pct=0.5)
print()
print(f"  Mean natural-form matches per random target at 0.5%: {mean_hits_05:.2f} +/- {std_hits_05:.2f}")
print(f"  Expected total over 4 parameters: {4*mean_hits_05:.2f}")
print(f"  Observed: {total_hits_05pct}")

results["step_5_false_positive_assessment"] = {
    "n_natural_forms": len(natural_forms),
    "range_log_uniform": [math.log(0.10), math.log(1.00)],
    "expected_hits_per_param_1pct": float(mean_hits),
    "expected_hits_per_param_05pct": float(mean_hits_05),
    "observed_total_hits_1pct": total_hits_1pct,
    "observed_total_hits_05pct": total_hits_05pct,
    "observed_within_1sigma": total_hits_within_1sigma,
}

# ============================================================
# Step 6: Cross-sector implications
# ============================================================
# Question: do the same natural forms appear in lepton observables?

print()
print("=" * 76)
print("Step 6: cross-sector check -- do candidate forms predict anything in lepton sector?")
print("=" * 76)

# Lepton observables
m_e   = mpf("0.51099895069")
m_mu  = mpf("105.6583755")
m_tau = mpf("1776.86")
masses = [m_e, m_mu, m_tau]
sum_m = sum(masses)
sum_sqrt = sum(sqrt(m) for m in masses)
K_lep = sum_m / sum_sqrt**2

# Lepton observables to test
lep_observables = {
    "Koide K":                            K_lep,
    "K - 2/3":                             K_lep - mpf(2)/3,
    "m_mu / m_e":                          m_mu/m_e,
    "m_tau / m_mu":                        m_tau/m_mu,
    "log(m_mu/m_e)":                       log(m_mu/m_e),
    "log(m_tau/m_mu)":                     log(m_tau/m_mu),
    "ratio of log spacings":               log(m_tau/m_mu) / log(m_mu/m_e),
    "sqrt(m_e/m_tau)":                     sqrt(m_e/m_tau),
}

print(f"  Lepton observable                vs natural form")
print(f"  {'-'*70}")
cross_matches = []
for lname, lval in lep_observables.items():
    for fname, fval in natural_forms.items():
        if fval == 0 or lval == 0:
            continue
        ratio = lval / fval
        # Test if ratio is close to a small simple rational
        for fac_name, fac_val in [
            ("1", mpf(1)), ("2", mpf(2)), ("3", mpf(3)), ("4", mpf(4)), ("5", mpf(5)),
            ("1/2", mpf(1)/2), ("1/3", mpf(1)/3), ("1/4", mpf(1)/4),
            ("2/3", mpf(2)/3), ("3/2", mpf(3)/2),
        ]:
            test = ratio / fac_val
            diff = abs(test - 1)
            if diff < mpf("0.005"):  # 0.5% threshold for cross-sector
                pct = float(diff)*100
                cross_matches.append({
                    "lepton_obs": lname,
                    "lepton_value": str(lval),
                    "natural_form": fname,
                    "factor": fac_name,
                    "rel_diff_pct": pct,
                })

# Print only sub-0.5% (already filtered)
shown = 0
for cm in sorted(cross_matches, key=lambda x: x["rel_diff_pct"]):
    if cm["rel_diff_pct"] < 0.3:
        print(f"  {cm['lepton_obs']:<28s}  ~  {cm['factor']:>6s} * {cm['natural_form']:<35s}  ({cm['rel_diff_pct']:.4f}%)")
        shown += 1
if shown == 0:
    print("  (no sub-0.3% non-trivial cross-sector matches)")

results["step_6_cross_sector"] = {"matches": cross_matches}

# ============================================================
# Step 7: Specific predictions of Candidate A in isolation
# ============================================================

print()
print("=" * 76)
print("Step 7: specific predictions of Candidate A (lambda = 1/sqrt(Vol S^3)) only")
print("=" * 76)

# If lambda = 1/sqrt(Vol(S^3)), what are the implied amplitude scales for
# OTHER framework manifolds?

print()
print("  Hypothesis: lambda is the 'natural amplitude on the framework manifold S^3'.")
print("  Then 1/sqrt(Vol(S^k)) gives natural amplitudes on other manifolds:")
for k in [1, 2, 3, 5]:
    if k == 1: vol = 2*pi
    elif k == 2: vol = 4*pi
    elif k == 3: vol = 2*pi*pi
    elif k == 5: vol = pi**3
    amp = 1/sqrt(vol)
    print(f"    1/sqrt(Vol(S^{k})) = {mp.nstr(amp, 8)}")

# Check: do any framework observables sit at 1/sqrt(Vol(S^k)) for k = 1, 2, 5?
# Useful framework observables to test:
print()
print("  Other Wolfenstein parameters vs 1/sqrt(Vol(S^k)) family:")
print(f"    A       = {mp.nstr(A_w, 8)} (PDG)")
print(f"    rho_bar = {mp.nstr(R, 8)}")
print(f"    eta_bar = {mp.nstr(H, 8)}")
print()
print(f"    Closest 1/sqrt(Vol(S^k)) for each:")
for pname, pval in [("A", get_best_value("A")), ("rho_bar", get_best_value("rho_bar")),
                     ("eta_bar", get_best_value("eta_bar"))]:
    best_k_match = None
    for k in [1, 2, 3, 5]:
        if k == 1: vol = 2*pi
        elif k == 2: vol = 4*pi
        elif k == 3: vol = 2*pi*pi
        elif k == 5: vol = pi**3
        amp = 1/sqrt(vol)
        diff = (pval - amp) / pval
        if best_k_match is None or abs(float(diff)) < abs(float(best_k_match[1])):
            best_k_match = (k, diff, amp)
    print(f"    {pname:<10s}  closest k={best_k_match[0]}, "
          f"value {mp.nstr(best_k_match[2], 6)}, diff {float(best_k_match[1])*100:+7.3f}%")

# Save
print()
print("=" * 76)
print("Saving full results to debug/data/w3_lambda_predictive_verification.json")
print("=" * 76)


def to_json_safe(obj):
    if isinstance(obj, mp.mpf):
        return str(obj)
    if isinstance(obj, dict):
        return {k: to_json_safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [to_json_safe(v) for v in obj]
    if isinstance(obj, np.floating):
        return float(obj)
    return obj


with open(OUT_DIR / "w3_lambda_predictive_verification.json", "w") as f:
    json.dump(to_json_safe(results), f, indent=2)
