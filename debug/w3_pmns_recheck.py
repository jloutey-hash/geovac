"""
Sprint W3 Test 1 (forward plan): PMNS recheck with the SAME prespecified 30-form
natural-form basis used for the CKM Wolfenstein verification.

Question being tested:
  Is the M1/M2 spectral identification of CKM Wolfenstein parameters
  (lambda, A, rho_bar, eta_bar) sector-universal in SM mixing matrices?

If PMNS angles fit M1/M2 forms within 1 PDG sigma using the same basis,
the M1/M2 split is sector-universal and W3 strengthens dramatically.
If PMNS does not fit (more likely outcome based on the first-pass probe),
W3 is restricted to the non-Majorana / quark-mixing sector.

Methodology:
  1. PRESPECIFIED basis frozen from w3_lambda_predictive_verification.py.
     (Do NOT modify; that would invalidate the universality test.)
  2. PMNS parameters from PDG 2024 / NuFIT 5.3 (Normal Hierarchy):
     sin^2 theta_12 = 0.307 +- 0.013
     sin^2 theta_23 = 0.546 +- 0.021  (octant ambiguity flagged)
     sin^2 theta_13 = 0.0220 +- 0.0007
     delta_CP_PMNS  = 230 +- 25 deg   (loose PDG)
     Plus derived: J_PMNS, sin theta_ij, mixing angles in deg/rad.
  3. For each parameter, search basis for matches at <1%, <0.5%, within 1 sigma PDG.
  4. Statistical accounting using random-target log-uniform reference (matches CKM script).
  5. ALSO compare to TBM benchmark (sin^2 theta_12 = 1/3, sin^2 theta_23 = 1/2,
     sin^2 theta_13 = 0). TBM is pure-rational, NOT in the M1/M2 rings.
     Report which (M1/M2 vs TBM-rational) is closer.
  6. If anything fits cleanly, derive predictions analogous to A^2 = 2 eta_bar
     and tan(delta_CP) = pi * ln 2; test against PDG.
  7. Honest verdict.

Run time: < 5 min at dps=80.
"""

import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from mpmath import mpf, pi, log, sqrt, exp, sin, cos, atan, atan2

mp.mp.dps = 80
OUT_DIR = Path("debug/data")
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# PMNS parameters from PDG 2024 / NuFIT 5.3 Normal Hierarchy
# ============================================================
# Source: PDG 2024 Review (Particle Data Group), with 1-sigma uncertainties.
# NuFIT 5.3 (Esteban et al. 2024) values are nearly identical at 1 sigma.
# Octant ambiguity in sin^2 theta_23 (lower vs upper octant) is flagged but
# not folded in -- best-fit value used.

# Primary PMNS parameters
sin2_t12 = (mpf("0.307"),  mpf("0.013"))   # solar; ~4% relative
sin2_t23 = (mpf("0.546"),  mpf("0.021"))   # atmospheric; ~4% relative; octant ambiguity
sin2_t13 = (mpf("0.0220"), mpf("0.0007"))  # reactor; ~3% relative
delta_cp = (mpf("230.0"),  mpf("25.0"))    # CP phase in deg; very loose PDG

# Derived (sin theta, sqrt(sin^2)), mixing angles in deg/rad
def derive_pmns():
    s12_2, _ = sin2_t12
    s23_2, _ = sin2_t23
    s13_2, _ = sin2_t13
    dcp,    _ = delta_cp

    s12 = sqrt(s12_2);  c12 = sqrt(1 - s12_2)
    s23 = sqrt(s23_2);  c23 = sqrt(1 - s23_2)
    s13 = sqrt(s13_2);  c13 = sqrt(1 - s13_2)
    dcp_rad = dcp * pi / 180

    # Jarlskog (CP-violation amplitude)
    J_pmns = c12 * c23 * c13**2 * s12 * s23 * s13 * sin(dcp_rad)
    # PDG sigma propagation (rough quadrature; dominated by delta_CP)
    # |J| has rough relative uncertainty ~25 deg / 230 deg ~ 11% -> 1 sigma in J ~ 11%
    # But also angle uncertainties contribute. Use 15% as conservative bound.
    sig_J = abs(J_pmns) * mpf("0.15")

    # Mixing angles in radians and degrees
    th12 = mp.asin(s12); th23 = mp.asin(s23); th13 = mp.asin(s13)
    th12_deg = th12 * 180/pi; th23_deg = th23 * 180/pi; th13_deg = th13 * 180/pi
    # Sigma in angle (deg) from sigma in sin^2:
    # sin^2 theta = X => theta = asin(sqrt X); dtheta/dX = 1/(2*sin theta * cos theta) = 1/sin(2 theta)
    # so sigma_theta = sigma_X / sin(2 theta).
    sig_th12_rad = sin2_t12[1] / sin(2*th12)
    sig_th23_rad = sin2_t23[1] / sin(2*th23)
    sig_th13_rad = sin2_t13[1] / sin(2*th13)
    sig_th12_deg = sig_th12_rad * 180/pi
    sig_th23_deg = sig_th23_rad * 180/pi
    sig_th13_deg = sig_th13_rad * 180/pi

    return {
        "s12": (s12, sin2_t12[1] / (2*s12)),  # sin theta with sigma
        "c12": (c12, sin2_t12[1] / (2*c12)),
        "s23": (s23, sin2_t23[1] / (2*s23)),
        "c23": (c23, sin2_t23[1] / (2*c23)),
        "s13": (s13, sin2_t13[1] / (2*s13)),
        "c13": (c13, sin2_t13[1] / (2*c13)),
        "th12_deg": (th12_deg, sig_th12_deg),
        "th23_deg": (th23_deg, sig_th23_deg),
        "th13_deg": (th13_deg, sig_th13_deg),
        "th12_rad": (th12, sig_th12_rad),
        "th23_rad": (th23, sig_th23_rad),
        "th13_rad": (th13, sig_th13_rad),
        "J_pmns":   (J_pmns, sig_J),
        "delta_cp_rad": (delta_cp[0]*pi/180, delta_cp[1]*pi/180),
    }


derived = derive_pmns()

# Test set: primary + selected derived observables
# Use a meaningful subset (omit redundancies); each is an independent test
pmns_observables = {
    "sin2_theta_12":   sin2_t12,
    "sin2_theta_23":   sin2_t23,
    "sin2_theta_13":   sin2_t13,
    "delta_CP_deg":    delta_cp,
    # Derived (informative additions; same dataset, different parameterizations)
    "sin_theta_12":    derived["s12"],
    "sin_theta_23":    derived["s23"],
    "sin_theta_13":    derived["s13"],
    "theta_12_deg":    derived["th12_deg"],
    "theta_23_deg":    derived["th23_deg"],
    "theta_13_deg":    derived["th13_deg"],
    "J_PMNS":          derived["J_pmns"],
}


# ============================================================
# Step 1: PRESPECIFIED natural-form basis (FROZEN, identical to CKM script)
# ============================================================
# DO NOT MODIFY. This basis is copied verbatim from
# debug/w3_lambda_predictive_verification.py.

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

# Tag each form by Mellin family for the M1/M2/Other split
# M1 = Hopf-base / pi-family (volumes of S^k, pi powers)
# M2 = chirality / ln(2)-family (Hurwitz at half-integer shift)
# Other = framework Paper 2/18 numbers, golden ratio, alpha
M1_forms = {
    "1/sqrt(Vol(S^1)) = 1/sqrt(2 pi)",
    "1/sqrt(Vol(S^2)) = 1/sqrt(4 pi)",
    "1/sqrt(Vol(S^3)) = 1/(pi sqrt 2)",
    "1/sqrt(Vol(S^5)) = 1/sqrt(pi^3)",
    "1/Vol(S^1) = 1/(2 pi)",
    "1/Vol(S^2) = 1/(4 pi)",
    "1/Vol(S^3) = 1/(2 pi^2)",
    "pi/14",
    "pi/16",
    "pi/(2 phi)",
    "1/pi",
    "1/(2 pi^2)",
    "sqrt(2)/(2 pi)",
    "1/(pi sqrt(pi))",
}
M2_forms = {
    "ln(2)",
    "ln(2)/2",
    "ln(2)/4",
    "sqrt(ln 2)",
    "1/ln(2)",
    "log(2)/pi",  # mixed M1/M2 -- count as M2 because ln(2) is the operative scale
}


def family(form_name):
    if form_name in M1_forms:
        return "M1"
    if form_name in M2_forms:
        return "M2"
    return "Other"


# ============================================================
# Step 2: Search natural-form basis for each PMNS parameter
# ============================================================

def find_matches(target_value, target_sigma, target_label,
                 forms, threshold_pct=2.0):
    matches = []
    for fname, fval in forms.items():
        if fval == 0 or target_value == 0:
            continue
        diff = (fval - target_value) / target_value
        within_sigma = abs(fval - target_value) / target_sigma if target_sigma > 0 else mpf("inf")
        if abs(float(diff)*100) < threshold_pct:
            matches.append({
                "form": fname,
                "family": family(fname),
                "value": str(fval),
                "rel_deviation_pct": float(diff)*100,
                "within_pdg_sigma": float(within_sigma),
            })
    matches.sort(key=lambda x: abs(x["rel_deviation_pct"]))
    return matches


print("=" * 76)
print("Step 2: natural-form match search per PMNS observable")
print("=" * 76)

results = {"step_2_natural_form_matches": {}}
total_tests = 0
total_hits_1pct = 0
total_hits_05pct = 0
total_hits_within_1sigma = 0
hits_M1 = 0
hits_M2 = 0

# Restrict the count of "primary" parameters for fair statistical comparison
# to the CKM script (which used 4 primary params). The 4 primary PMNS params:
#   sin^2 theta_12, sin^2 theta_23, sin^2 theta_13, delta_CP
PRIMARY_PARAMS = ["sin2_theta_12", "sin2_theta_23", "sin2_theta_13", "delta_CP_deg"]

for pname, (pval, psig) in pmns_observables.items():
    is_primary = pname in PRIMARY_PARAMS
    print()
    label = "primary" if is_primary else "derived"
    print(f"--- {pname} = {mp.nstr(pval, 8)} +/- {mp.nstr(psig, 4)} ({label}) ---")
    matches = find_matches(pval, psig, pname, natural_forms, threshold_pct=2.0)
    n_tests = len(natural_forms)
    if is_primary:
        total_tests += n_tests
    n_1pct  = sum(1 for m in matches if abs(m["rel_deviation_pct"]) < 1.0)
    n_05pct = sum(1 for m in matches if abs(m["rel_deviation_pct"]) < 0.5)
    n_1sig  = sum(1 for m in matches if abs(m["within_pdg_sigma"]) < 1.0)
    if is_primary:
        total_hits_1pct += n_1pct
        total_hits_05pct += n_05pct
        total_hits_within_1sigma += n_1sig
        for m in matches:
            if abs(m["within_pdg_sigma"]) < 1.0:
                if m["family"] == "M1":
                    hits_M1 += 1
                elif m["family"] == "M2":
                    hits_M2 += 1
    print(f"  forms tested: {n_tests}; <1%: {n_1pct}; <0.5%: {n_05pct}; "
          f"within 1 sigma PDG: {n_1sig}")
    results["step_2_natural_form_matches"][pname] = {
        "is_primary": is_primary,
        "value":  str(pval),
        "sigma":  str(psig),
        "matches": matches,
        "n_tested": n_tests,
        "n_1pct":   n_1pct,
        "n_05pct":  n_05pct,
        "n_1sigma": n_1sig,
    }
    for m in matches[:5]:
        marker = " *" if abs(m["within_pdg_sigma"]) < 1 else ""
        print(f"     [{m['family']:<2s}] {m['form']:<40s}  "
              f"diff {m['rel_deviation_pct']:+7.3f}%  "
              f"({m['within_pdg_sigma']:.2f} sigma){marker}")

print()
print(f"Aggregate over {len(PRIMARY_PARAMS)} primary PMNS parameters:")
print(f"  total tests = {total_tests}")
print(f"  hits <1%    = {total_hits_1pct}")
print(f"  hits <0.5%  = {total_hits_05pct}")
print(f"  within 1 PDG sigma = {total_hits_within_1sigma}")
print(f"  M1 family (Hopf-base / pi)        = {hits_M1}")
print(f"  M2 family (chirality / ln 2)      = {hits_M2}")


# ============================================================
# Step 3: Best match per primary PMNS parameter
# ============================================================

print()
print("=" * 76)
print("Step 3: Best natural-form fit per primary PMNS parameter")
print("=" * 76)

best = {}
for pname in PRIMARY_PARAMS:
    info = results["step_2_natural_form_matches"][pname]
    best[pname] = info["matches"][0] if info["matches"] else None

print()
for pname in PRIMARY_PARAMS:
    b = best.get(pname)
    if b:
        print(f"  {pname:<18s} ~ [{b['family']:<2s}] {b['form']:<40s}  "
              f"diff {b['rel_deviation_pct']:+7.3f}%  "
              f"({b['within_pdg_sigma']:.2f} sigma PDG)")
    else:
        print(f"  {pname:<18s} no <2% natural form match")

results["best_per_primary"] = {p: best[p] for p in PRIMARY_PARAMS}


# ============================================================
# Step 4: TBM (tri-bimaximal) benchmark comparison
# ============================================================
# TBM is the structural alternative for PMNS that is widely studied:
#   sin^2 theta_12 = 1/3, sin^2 theta_23 = 1/2, sin^2 theta_13 = 0
# These are PURE RATIONALS, not in the M1 (pi-family) or M2 (ln 2-family) rings.
# If PMNS sits closer to TBM rationals than to M1/M2 forms, that signals the
# leptonic mixing sector follows a different mechanism (group-theoretic
# discrete-symmetry breaking) rather than the spectral-zeta mechanism.

print()
print("=" * 76)
print("Step 4: TBM (tri-bimaximal) benchmark comparison")
print("=" * 76)

TBM_targets = {
    "sin2_theta_12": (mpf(1)/3,  "1/3 (TBM)"),
    "sin2_theta_23": (mpf(1)/2,  "1/2 (TBM)"),
    "sin2_theta_13": (mpf(0),    "0 (TBM)"),
}

print()
print(f"  {'Parameter':<18s} {'PDG':>10s} {'TBM':>10s} {'PDG sigma':>10s} {'Best M1/M2':>20s} {'M1/M2 sigma':>12s}")
print(f"  {'-'*100}")

tbm_table = {}
for pname in ["sin2_theta_12", "sin2_theta_23", "sin2_theta_13"]:
    pval, psig = pmns_observables[pname]
    tbm_val, tbm_label = TBM_targets[pname]
    # TBM distance in PDG sigma
    if psig > 0:
        tbm_sigma_dist = abs(tbm_val - pval) / psig
    else:
        tbm_sigma_dist = mpf("inf")
    # Best M1/M2 match
    matches = results["step_2_natural_form_matches"][pname]["matches"]
    best_m1m2 = None
    best_m1m2_sigma = mpf("inf")
    for m in matches:
        if m["family"] in ("M1", "M2"):
            sd = abs(mpf(m["within_pdg_sigma"]))
            if sd < best_m1m2_sigma:
                best_m1m2_sigma = sd
                best_m1m2 = m
    if best_m1m2 is None:
        m1m2_str = "(no M1/M2 in basis < 2%)"
        m1m2_sigma_str = "-"
    else:
        m1m2_str = f"[{best_m1m2['family']}] {best_m1m2['form'][:18]}"
        m1m2_sigma_str = f"{best_m1m2['within_pdg_sigma']:.2f}"
    print(f"  {pname:<18s} {mp.nstr(pval,5):>10s} {tbm_label:>10s} "
          f"{float(tbm_sigma_dist):>10.2f} {m1m2_str:>20s} {m1m2_sigma_str:>12s}")
    tbm_table[pname] = {
        "pdg_central": str(pval),
        "pdg_sigma":   str(psig),
        "tbm_value":   str(tbm_val),
        "tbm_sigma_distance": float(tbm_sigma_dist),
        "best_m1m2_form":  best_m1m2["form"]   if best_m1m2 else None,
        "best_m1m2_family": best_m1m2["family"] if best_m1m2 else None,
        "best_m1m2_sigma": best_m1m2["within_pdg_sigma"] if best_m1m2 else None,
    }

results["tbm_comparison"] = tbm_table


# ============================================================
# Step 5: Statistical assessment (false-positive rate)
# ============================================================

print()
print("=" * 76)
print("Step 5: false-positive rate by random-target log-uniform reference")
print("=" * 76)

# Same methodology as CKM script: how many natural-form matches are expected
# by chance for a random target in a comparable parameter range?

# PMNS primary parameters span:
#   sin^2 theta_13 ~ 0.022 (smallest)
#   sin^2 theta_23 ~ 0.546 (largest of the sin^2)
#   delta_CP ~ 230 deg (different scale)
# We compare like-to-like: log-uniform range matching the parameter under test.
# For sin^2 angles the natural log range is [0.01, 1.0]; for delta_CP it's
# [10, 360] but we test it as the angle in degrees so range [10, 400].
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


# For three sin^2 parameters, the relevant range covers [0.02, 0.6]
# matching the PMNS sin^2 values. delta_CP is scaled separately
# but for fair comparison with CKM we use the same log-uniform reference
# range used in the CKM verification: [0.10, 1.00].
log_range_sin2 = (math.log(0.02), math.log(1.00))   # natural for sin^2 (PMNS-specific)
log_range_ckm  = (math.log(0.10), math.log(1.00))   # used in CKM script (cross-comparable)

mean_sin2, std_sin2 = expected_hits_random(natural_forms, log_range_sin2,
                                            n_trials=10000, threshold_pct=1.0)
mean_ckm,  std_ckm  = expected_hits_random(natural_forms, log_range_ckm,
                                            n_trials=10000, threshold_pct=1.0)

# Match what the CKM script did: same 4 primary params, same log range [0.10, 1.00]
# even though delta_CP is at a different scale. This is most directly comparable.
expected_4 = 4 * mean_ckm
expected_4_sigma = math.sqrt(4) * std_ckm

print(f"  basis size (forms tested) = {len(natural_forms)}")
print(f"  PMNS primary params       = {len(PRIMARY_PARAMS)}")
print(f"  observed hits <1%         = {total_hits_1pct}")
print(f"  observed hits <0.5%       = {total_hits_05pct}")
print(f"  observed within 1 sigma   = {total_hits_within_1sigma}")
print()
print(f"  Random target log-uniform [0.10, 1.00] (matches CKM script):")
print(f"    mean hits per param at <1%  = {mean_ckm:.3f} +/- {std_ckm:.3f}")
print(f"    expected total over 4 prim  = {expected_4:.2f} +/- {expected_4_sigma:.2f}")
print(f"    observed                     = {total_hits_1pct}")
deviation_pmns_chance = (total_hits_1pct - expected_4) / max(expected_4_sigma, 1e-9)
print(f"    sigma above chance           = {deviation_pmns_chance:+.2f}")
print()
print(f"  Random target log-uniform [0.02, 1.00] (PMNS-specific sin^2 range):")
print(f"    mean hits per param at <1%  = {mean_sin2:.3f} +/- {std_sin2:.3f}")

mean_05_ckm,  _ = expected_hits_random(natural_forms, log_range_ckm,
                                       n_trials=10000, threshold_pct=0.5)
print()
print(f"  At <0.5% threshold, random log-uniform [0.10, 1.00]:")
print(f"    mean hits per param at <0.5% = {mean_05_ckm:.3f}")
print(f"    expected over 4 params       = {4*mean_05_ckm:.2f}")
print(f"    observed                      = {total_hits_05pct}")

results["step_5_false_positive_assessment"] = {
    "n_natural_forms": len(natural_forms),
    "n_primary_pmns_params": len(PRIMARY_PARAMS),
    "log_range_ckm": [math.log(0.10), math.log(1.00)],
    "log_range_pmns_specific": [math.log(0.02), math.log(1.00)],
    "expected_per_param_1pct_ckm_range":   float(mean_ckm),
    "expected_per_param_05pct_ckm_range":  float(mean_05_ckm),
    "expected_per_param_1pct_pmns_range":  float(mean_sin2),
    "expected_total_4params_1pct_ckm_range":  float(expected_4),
    "expected_total_4params_1pct_sigma":       float(expected_4_sigma),
    "observed_total_hits_1pct":      int(total_hits_1pct),
    "observed_total_hits_05pct":     int(total_hits_05pct),
    "observed_within_1sigma":        int(total_hits_within_1sigma),
    "observed_M1_family_hits":       int(hits_M1),
    "observed_M2_family_hits":       int(hits_M2),
    "deviation_above_chance_sigma":  float(deviation_pmns_chance),
}


# ============================================================
# Step 6: CKM comparison (side-by-side numbers from JSON if present)
# ============================================================

print()
print("=" * 76)
print("Step 6: side-by-side CKM vs PMNS")
print("=" * 76)

ckm_json = OUT_DIR / "w3_lambda_predictive_verification.json"
ckm_data = None
if ckm_json.exists():
    with open(ckm_json) as f:
        ckm_data = json.load(f)
    ckm_obs    = ckm_data["step_5_false_positive_assessment"]["observed_total_hits_1pct"]
    ckm_obs_05 = ckm_data["step_5_false_positive_assessment"]["observed_total_hits_05pct"]
    ckm_obs_1s = ckm_data["step_5_false_positive_assessment"]["observed_within_1sigma"]
    ckm_exp    = ckm_data["step_5_false_positive_assessment"]["expected_hits_per_param_1pct"]
    print(f"  {'Metric':<35s} {'CKM':>10s} {'PMNS':>10s}")
    print(f"  {'-'*60}")
    print(f"  {'observed hits at <1%':<35s} {ckm_obs:>10d} {total_hits_1pct:>10d}")
    print(f"  {'observed hits at <0.5%':<35s} {ckm_obs_05:>10d} {total_hits_05pct:>10d}")
    print(f"  {'within 1 PDG sigma':<35s} {ckm_obs_1s:>10d} {total_hits_within_1sigma:>10d}")
    print(f"  {'expected per param at <1% (chance)':<35s} {ckm_exp:>10.3f} {mean_ckm:>10.3f}")
else:
    print("  (CKM JSON not found; comparison skipped)")
    ckm_obs = ckm_obs_05 = ckm_obs_1s = None

results["ckm_pmns_comparison"] = {
    "ckm_observed_1pct":  ckm_obs,
    "ckm_observed_05pct": ckm_obs_05,
    "ckm_observed_1sigma": ckm_obs_1s,
    "pmns_observed_1pct":  total_hits_1pct,
    "pmns_observed_05pct": total_hits_05pct,
    "pmns_observed_1sigma": total_hits_within_1sigma,
}


# ============================================================
# Step 7: Derived predictions analogous to A^2 = 2 eta_bar
# ============================================================
# In the CKM result, A^2 = 2 eta_bar followed from
#   A    = sqrt(ln 2)         => A^2 = ln 2
#   eta  = (ln 2)/2           => 2 eta = ln 2
# i.e. A^2 = 2 eta_bar is a derived relation among the four candidates that
# survives PDG at 0.52 sigma. It reduces 4 free parameters to 3 + 1 constraint.
#
# Question: do the best-fit M1/M2 candidates for any subset of PMNS params
# imply derived relations that survive PDG?

print()
print("=" * 76)
print("Step 7: derived relations in PMNS sector (analogous to A^2 = 2 eta_bar)")
print("=" * 76)

derived_relations = []

# Loop over pairs of primary params with M1/M2 best fits within 1 sigma PDG
within_sigma_fits = {}
for pname in PRIMARY_PARAMS:
    info = results["step_2_natural_form_matches"][pname]
    if info["matches"] and abs(info["matches"][0]["within_pdg_sigma"]) < 1.0:
        b = info["matches"][0]
        if b["family"] in ("M1", "M2"):
            within_sigma_fits[pname] = b

print(f"  PMNS primary params with within-1-sigma M1/M2 fit:")
for p, b in within_sigma_fits.items():
    print(f"    {p:<18s} ~ [{b['family']}] {b['form']}")

if len(within_sigma_fits) >= 2:
    print()
    print(f"  Looking for low-integer-coefficient relations among the {len(within_sigma_fits)} fits...")
    # Test x_i / x_j against simple rationals 1, 2, 1/2, 3, 1/3, ...
    items = list(within_sigma_fits.items())
    for i in range(len(items)):
        for j in range(i+1, len(items)):
            pi_name, fi = items[i]
            pj_name, fj = items[j]
            xi = mpf(fi["value"])
            xj = mpf(fj["value"])
            ratio = xi / xj
            # Test against simple rationals and pi powers
            test_targets = [
                (mpf(1),    "1"),
                (mpf(2),    "2"),
                (mpf(1)/2,  "1/2"),
                (mpf(3),    "3"),
                (mpf(1)/3,  "1/3"),
                (mpf(4),    "4"),
                (mpf(1)/4,  "1/4"),
                (pi,        "pi"),
                (1/pi,      "1/pi"),
                (mpf(2)*pi, "2 pi"),
                (1/(2*pi),  "1/(2 pi)"),
            ]
            for tval, tname in test_targets:
                if abs(ratio - tval) / tval < mpf("0.05"):
                    rel_pct = float(abs(ratio - tval) / tval * 100)
                    derived_relations.append({
                        "lhs": pi_name,
                        "rhs": pj_name,
                        "candidate_relation": f"{pi_name} / {pj_name} = {tname}",
                        "fit_lhs": fi["form"],
                        "fit_rhs": fj["form"],
                        "rel_dev_pct": rel_pct,
                    })

if not derived_relations:
    print(f"    (no clean derived relations found at <5% threshold)")
else:
    for r in derived_relations:
        print(f"    {r['candidate_relation']:<40s}  rel dev {r['rel_dev_pct']:.3f}%")
        print(f"      ({r['fit_lhs']}) / ({r['fit_rhs']})")

results["derived_relations"] = derived_relations


# ============================================================
# Step 8: tan(delta_CP) probe analogous to CKM's tan = pi * ln 2
# ============================================================
# In the CKM result, tan(delta_CP_CKM) = eta_bar / rho_bar = (ln 2 / 2) / (1/(2 pi))
#                                      = pi * ln 2  (M2 prototype * M1 prototype)
# Test whether tan(delta_CP_PMNS) also lands on a clean M1*M2 product.

print()
print("=" * 76)
print("Step 8: tan(delta_CP_PMNS) vs natural M1*M2 / M2/M1 products")
print("=" * 76)

dcp_rad = derived["delta_cp_rad"][0]
tan_dcp = mp.tan(dcp_rad)
print(f"  delta_CP_PMNS = {float(delta_cp[0])} +/- {float(delta_cp[1])} deg")
print(f"  tan(delta_CP_PMNS) = {mp.nstr(tan_dcp, 8)}")
print()

# Build M1*M2 and M2/M1 product candidates from the basis
M1_vals = {n: natural_forms[n] for n in M1_forms}
M2_vals = {n: natural_forms[n] for n in M2_forms}

tan_dcp_matches = []
for m1_name, m1_v in M1_vals.items():
    for m2_name, m2_v in M2_vals.items():
        for op, op_str in [("*", m1_v * m2_v), ("M2/M1", m2_v / m1_v), ("M1/M2", m1_v / m2_v)]:
            test_val = op_str
            if test_val == 0:
                continue
            diff = (test_val - tan_dcp) / tan_dcp
            if abs(float(diff)) < 0.05:
                tan_dcp_matches.append({
                    "operation": op,
                    "M1_form": m1_name,
                    "M2_form": m2_name,
                    "value":   str(test_val),
                    "rel_dev_pct": float(diff)*100,
                })

if tan_dcp_matches:
    print(f"  Matches at <5%:")
    for m in sorted(tan_dcp_matches, key=lambda x: abs(x["rel_dev_pct"]))[:6]:
        print(f"    {m['operation']:<8s} ({m['M1_form'][:25]}, {m['M2_form'][:15]})  "
              f"diff {m['rel_dev_pct']:+7.3f}%")
else:
    print(f"  (no M1*M2 or M2/M1 products within 5% of tan(delta_CP_PMNS))")

results["tan_delta_cp_probe"] = {
    "delta_cp_deg": str(delta_cp[0]),
    "delta_cp_sigma_deg": str(delta_cp[1]),
    "tan_delta_cp": str(tan_dcp),
    "matches": tan_dcp_matches,
}


# ============================================================
# Step 9: Cross-rationals (PMNS angles vs simple fractions of pi)
# ============================================================
# Test theta_13 in deg vs simple p/q*deg, theta_23 vs 45 deg, etc.

print()
print("=" * 76)
print("Step 9: PMNS angles vs simple rational fractions of pi (in radians) and degrees")
print("=" * 76)

simple_targets = [
    (mpf(45), "45 deg = pi/4"),
    (mpf(33.56), "asin(1/sqrt 3) (TBM theta_12)"),
    (mpf(0), "0 deg (TBM theta_13)"),
    (mpf(35.264), "asin(1/sqrt 3) other"),
    (pi/6 * 180/pi, "pi/6 = 30 deg"),
    (pi/4 * 180/pi, "pi/4 = 45 deg"),
    (pi/8 * 180/pi, "pi/8 = 22.5 deg"),
    (pi/12 * 180/pi, "pi/12 = 15 deg"),
    (mp.acos(1/sqrt(3))*180/pi, "acos(1/sqrt 3) ~ 54.74 deg"),
]

print(f"  {'Param':<20s} {'PDG':>14s} {'sigma':>10s} {'Best simple fit':>30s} {'sigma':>10s}")
print(f"  {'-'*90}")
simple_table = {}
for pname in ["theta_12_deg", "theta_23_deg", "theta_13_deg"]:
    pval, psig = pmns_observables[pname]
    best_simple = None
    best_simple_sigma = mpf("inf")
    for tval, tlabel in simple_targets:
        if psig > 0:
            sd = abs(tval - pval) / psig
            if sd < best_simple_sigma:
                best_simple_sigma = sd
                best_simple = (tval, tlabel)
    print(f"  {pname:<20s} {mp.nstr(pval,5):>14s} {mp.nstr(psig,3):>10s} "
          f"{best_simple[1]:>30s} {float(best_simple_sigma):>10.2f}")
    simple_table[pname] = {
        "pdg":  str(pval),
        "sigma": str(psig),
        "best_simple_target": best_simple[1],
        "best_simple_value":  str(best_simple[0]),
        "best_simple_sigma":  float(best_simple_sigma),
    }
results["simple_rationals_in_deg"] = simple_table


# ============================================================
# VERDICT
# ============================================================

print()
print("=" * 76)
print("VERDICT")
print("=" * 76)
print()

# Build the verdict string based on the numerical evidence
# Criterion 1: how many primary PMNS params have within-1-sigma M1/M2 fit?
n_primary_within_sigma_m1m2 = sum(1 for p in PRIMARY_PARAMS
                                   if results["step_2_natural_form_matches"][p]["matches"]
                                   and abs(results["step_2_natural_form_matches"][p]["matches"][0]["within_pdg_sigma"]) < 1.0
                                   and results["step_2_natural_form_matches"][p]["matches"][0]["family"] in ("M1", "M2"))

# Criterion 2: TBM closer than M1/M2 for sin^2 angles?
tbm_wins = 0
m1m2_wins = 0
ties = 0
for pname in ["sin2_theta_12", "sin2_theta_23", "sin2_theta_13"]:
    tbm_sd = tbm_table[pname]["tbm_sigma_distance"]
    m1m2_sd = tbm_table[pname]["best_m1m2_sigma"]
    if m1m2_sd is None:
        tbm_wins += 1
    elif tbm_sd < m1m2_sd:
        tbm_wins += 1
    elif m1m2_sd < tbm_sd:
        m1m2_wins += 1
    else:
        ties += 1

print(f"  Primary PMNS params with within-1-sigma M1/M2 fit: {n_primary_within_sigma_m1m2} of 4")
print(f"  TBM-rational vs best-M1/M2 (sin^2 angles, lower sigma wins):")
print(f"    TBM wins: {tbm_wins}, M1/M2 wins: {m1m2_wins}, tie: {ties}")
print(f"  Aggregate hits at <1%: {total_hits_1pct} (CKM had 7 of similar 4-param test)")
print(f"  Sigma above chance: {deviation_pmns_chance:+.2f}")
print()

if n_primary_within_sigma_m1m2 >= 3 and m1m2_wins >= tbm_wins:
    verdict_str = "SECTOR-UNIVERSAL: PMNS fits M1/M2 like CKM. W3 strengthens dramatically."
elif n_primary_within_sigma_m1m2 >= 2:
    verdict_str = "PARTIAL: some PMNS params fit M1/M2; M1/M2 split is incompletely sector-universal."
elif tbm_wins > m1m2_wins:
    verdict_str = "TBM-INSTEAD: PMNS sits closer to TBM rationals than M1/M2 forms; leptonic mixing follows a different mechanism (group-theoretic discrete-symmetry breaking) rather than spectral-zeta. W3 is sector-restricted to CKM."
else:
    verdict_str = "SECTOR-RESTRICTED: PMNS does not fit M1/M2 cleanly; CKM result does not generalize. W3 universality hypothesis is weakened."

print(f"  VERDICT: {verdict_str}")

results["verdict"] = {
    "n_primary_within_sigma_m1m2": n_primary_within_sigma_m1m2,
    "tbm_wins": tbm_wins,
    "m1m2_wins": m1m2_wins,
    "ties": ties,
    "total_hits_1pct": total_hits_1pct,
    "sigma_above_chance": float(deviation_pmns_chance),
    "verdict_string": verdict_str,
}


# ============================================================
# Save results
# ============================================================

print()
print("=" * 76)
print("Saving full results to debug/data/w3_pmns_recheck.json")
print("=" * 76)


def to_json_safe(obj):
    if isinstance(obj, mp.mpf):
        return str(obj)
    if isinstance(obj, dict):
        return {k: to_json_safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [to_json_safe(v) for v in obj]
    if isinstance(obj, tuple):
        return [to_json_safe(v) for v in obj]
    if isinstance(obj, (np.floating, np.integer)):
        return float(obj) if isinstance(obj, np.floating) else int(obj)
    return obj


with open(OUT_DIR / "w3_pmns_recheck.json", "w") as f:
    json.dump(to_json_safe(results), f, indent=2)

print()
print("DONE.")
