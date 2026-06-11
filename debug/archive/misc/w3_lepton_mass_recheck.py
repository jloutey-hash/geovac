"""
Sprint W3 follow-up, Test 2 of forward plan: Lepton mass spectrum recheck
with prespecified-basis methodology.

Background:
  On 2026-05-08 (evening), w3_lambda_predictive_verification.py applied a
  30-form prespecified GeoVac-internal natural-form basis to the four
  Wolfenstein CKM parameters and found 8 PDG-1-sigma matches (~6-10 sigma
  above chance). The lepton sector independently has the Koide cone fact
  (K = 2/3 at 1 arcsec). But m_e, m_mu, m_tau as standalone observables --
  ratios, log-spacings -- have not been searched with the same methodology.

Hypothesis under test:
  If W3 (calibration data lives in master Mellin engine M1/M2 rings) is
  universal, the dimensionless lepton-mass observables should also fit
  M1/M2/M3 forms within experimental uncertainty. If they don't, W3 is
  restricted to the mixing-matrix sector.

Methodology (locked):
  Step 0: Prespecify the observable panel BEFORE looking at data
          (this script header).
  Step 1: Use the EXACT same 30-form basis as
          w3_lambda_predictive_verification.py. Do not extend.
  Step 2: For each observable, find matches at <1%, <0.5%, and within 1
          measurement-precision sigma.
  Step 3: Statistical accounting.
  Step 4: Cross-sector check: do any lepton observables share a form
          identified for a CKM Wolfenstein parameter?
  Step 5: Honest verdict.

Prespecified observable panel (LOCKED, before searching the basis):
  Direct mass ratios:
    1. r_eμ      = m_μ / m_e            (~206.77)
    2. r_μτ      = m_τ / m_μ            (~16.82)
    3. r_eτ      = m_τ / m_e            (~3477.2)
    4. 1/r_eμ    = m_e / m_μ            (~0.004836)
    5. 1/r_μτ    = m_μ / m_τ            (~0.05946)
    6. 1/r_eτ    = m_e / m_τ            (~0.000288)
  Logarithmic spacings:
    7. log(r_eμ)                        (~5.331)
    8. log(r_μτ)                        (~2.823)
    9. log(r_eτ)                        (~8.154)
   10. log spacing ratio = log(r_μτ)/log(r_eμ)   (~0.530)
  Square-root mass ratios (Koide-natural):
   11. sqrt(m_e/m_μ)                    (~0.06955)
   12. sqrt(m_μ/m_τ)                    (~0.2438)
   13. sqrt(m_e/m_τ)                    (~0.01696)
   14. sqrt(m_μ/m_e)                    (~14.378)
   15. sqrt(m_τ/m_μ)                    (~4.101)
  Koide observables:
   16. Koide K = (Σm)/(Σ√m)²            (~0.66666)
   17. Koide deviation K - 2/3
   18. Cone half-angle (degrees)        (~44.99974)
   19. Cone half-angle from democratic in radians
   20. Cone deviation from 45° (radians)
  Mass-fraction observables (democratic-normalized):
   21. m_e / (m_e + m_μ + m_τ)
   22. m_μ / (m_e + m_μ + m_τ)
   23. m_τ / (m_e + m_μ + m_τ)
   24. m_μ - m_e fraction:
       (m_μ - m_e) / (m_τ - m_e)        (~0.05884)
   25. (m_τ - m_μ) / m_τ                (~0.9405)
  sqrt-mass-fraction observables:
   26. √m_e / (√m_e + √m_μ + √m_τ)
   27. √m_μ / (√m_e + √m_μ + √m_τ)
   28. √m_τ / (√m_e + √m_μ + √m_τ)
  Geometric-series tests:
   29. m_μ² / (m_e * m_τ)               (geometric-series test, ≠ 1)
   30. log(m_μ²) - log(m_e * m_τ)       (= 0 if perfect geometric)

This is the panel. It is committed to the script header. The basis below
matches w3_lambda_predictive_verification.py verbatim. Only THEN do we run
the search.
"""

import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from mpmath import mpf, pi, log, sqrt, atan2

mp.mp.dps = 80
OUT_DIR = Path("debug/data")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# PDG 2024 charged-lepton masses (MeV / c^2) and uncertainties
# These match the values used in w3_lambda_predictive_verification.py
# (CODATA-2022 for m_e, m_mu; PDG 2024 for m_tau)
# ============================================================

m_e   = mpf("0.51099895069")
m_mu  = mpf("105.6583755")
m_tau = mpf("1776.86")

# CODATA-2022 / PDG-2024 1-sigma uncertainties
sigma_e   = mpf("0.00000000016")
sigma_mu  = mpf("0.0000023")
sigma_tau = mpf("0.12")

# ============================================================
# PRESPECIFIED natural-form basis (frozen, identical to
# w3_lambda_predictive_verification.py, lines 63-107 verbatim)
# ============================================================

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

assert len(natural_forms) == 30, f"basis must be 30 forms, got {len(natural_forms)}"

# ============================================================
# Compute the prespecified observable panel
# ============================================================

# Convenience
masses = [m_e, m_mu, m_tau]
sum_m  = sum(masses)
sqrt_m = [sqrt(m) for m in masses]
sum_sqrt_m = sum(sqrt_m)
K_lep  = sum_m / sum_sqrt_m**2

# Cone geometry: cos^2(angle) = (sum sqrt m_i)^2 / (3 sum m_i) = 1/(3K)
cos_sq_angle = (sum_sqrt_m * sum_sqrt_m) / (3 * sum_m)
angle_rad    = mp.acos(sqrt(cos_sq_angle))
angle_deg    = angle_rad * 180 / pi
cone_dev_deg = angle_deg - 45
cone_dev_rad = angle_rad - pi/4

observables = {
    # Direct mass ratios
    "1. r_eμ = m_μ/m_e":                m_mu / m_e,
    "2. r_μτ = m_τ/m_μ":                m_tau / m_mu,
    "3. r_eτ = m_τ/m_e":                m_tau / m_e,
    "4. 1/r_eμ = m_e/m_μ":              m_e / m_mu,
    "5. 1/r_μτ = m_μ/m_τ":              m_mu / m_tau,
    "6. 1/r_eτ = m_e/m_τ":              m_e / m_tau,
    # Logarithmic spacings
    "7. log(r_eμ)":                     log(m_mu / m_e),
    "8. log(r_μτ)":                     log(m_tau / m_mu),
    "9. log(r_eτ)":                     log(m_tau / m_e),
    "10. log spacing ratio":            log(m_tau / m_mu) / log(m_mu / m_e),
    # Square-root mass ratios
    "11. sqrt(m_e/m_μ)":                sqrt(m_e / m_mu),
    "12. sqrt(m_μ/m_τ)":                sqrt(m_mu / m_tau),
    "13. sqrt(m_e/m_τ)":                sqrt(m_e / m_tau),
    "14. sqrt(m_μ/m_e)":                sqrt(m_mu / m_e),
    "15. sqrt(m_τ/m_μ)":                sqrt(m_tau / m_mu),
    # Koide observables
    "16. Koide K = (Σm)/(Σ√m)²":        K_lep,
    "17. Koide K - 2/3":                K_lep - mpf(2)/3,
    "18. cone half-angle (deg)":        angle_deg,
    "19. cone half-angle (rad)":        angle_rad,
    "20. cone deviation from 45° (rad)":cone_dev_rad,
    # Mass-fraction observables
    "21. m_e / Σm":                     m_e / sum_m,
    "22. m_μ / Σm":                     m_mu / sum_m,
    "23. m_τ / Σm":                     m_tau / sum_m,
    "24. (m_μ - m_e)/(m_τ - m_e)":      (m_mu - m_e)/(m_tau - m_e),
    "25. (m_τ - m_μ)/m_τ":              (m_tau - m_mu) / m_tau,
    # sqrt-mass-fraction observables
    "26. √m_e / Σ√m":                   sqrt_m[0] / sum_sqrt_m,
    "27. √m_μ / Σ√m":                   sqrt_m[1] / sum_sqrt_m,
    "28. √m_τ / Σ√m":                   sqrt_m[2] / sum_sqrt_m,
    # Geometric-series tests
    "29. m_μ² / (m_e m_τ)":             (m_mu * m_mu) / (m_e * m_tau),
    "30. log(m_μ²) - log(m_e m_τ)":     log(m_mu*m_mu) - log(m_e*m_tau),
}

assert len(observables) == 30, f"panel must have 30 observables, got {len(observables)}"

# ============================================================
# Per-observable measurement-precision sigma propagation
# These are the "experimental sigma" for each observable, propagated
# from input mass uncertainties via partial derivatives.
# ============================================================

# Helper: relative uncertainty of each mass
rel_sigma = {
    "e":   sigma_e   / m_e,   # ~3.1e-10
    "mu":  sigma_mu  / m_mu,  # ~2.2e-8
    "tau": sigma_tau / m_tau, # ~6.8e-5
}

def sigma_ratio(num_label, den_label, value):
    """Relative sigma of m_a / m_b is sqrt((σ_a/m_a)² + (σ_b/m_b)²)."""
    rel = sqrt(rel_sigma[num_label]**2 + rel_sigma[den_label]**2)
    return value * rel

def sigma_log_ratio(num_label, den_label):
    """sigma of log(m_a/m_b) is sqrt((σ_a/m_a)² + (σ_b/m_b)²)."""
    return sqrt(rel_sigma[num_label]**2 + rel_sigma[den_label]**2)

def sigma_sqrt_ratio(num_label, den_label, value):
    """sigma of sqrt(m_a/m_b) is (1/2) sqrt(...)*value."""
    rel = mpf("0.5") * sqrt(rel_sigma[num_label]**2 + rel_sigma[den_label]**2)
    return value * rel

# Use Monte Carlo for the more involved observables (Koide K, cone, fractions).
rng = np.random.default_rng(11)
N_MC = 80_000
me_f, mmu_f, mtau_f = float(m_e), float(m_mu), float(m_tau)
se_f, smu_f, stau_f = float(sigma_e), float(sigma_mu), float(sigma_tau)

mc_samples = []
for _ in range(N_MC):
    pe = me_f + se_f * rng.normal()
    pmu = mmu_f + smu_f * rng.normal()
    ptau = mtau_f + stau_f * rng.normal()
    sm = pe + pmu + ptau
    ssr = math.sqrt(pe) + math.sqrt(pmu) + math.sqrt(ptau)
    K = sm / (ssr * ssr)
    cos_sq = (ssr*ssr) / (3 * sm)
    angle = math.acos(math.sqrt(cos_sq))
    angle_d = angle * 180 / math.pi
    mc_samples.append({
        "K": K,
        "K_minus_2_3": K - 2/3,
        "angle_deg": angle_d,
        "angle_rad": angle,
        "cone_dev_rad": angle - math.pi/4,
        "m_e_frac": pe / sm,
        "m_mu_frac": pmu / sm,
        "m_tau_frac": ptau / sm,
        "diff_frac": (pmu - pe) / (ptau - pe),
        "tau_minus_mu_over_tau": (ptau - pmu) / ptau,
        "sqrt_e_frac": math.sqrt(pe) / ssr,
        "sqrt_mu_frac": math.sqrt(pmu) / ssr,
        "sqrt_tau_frac": math.sqrt(ptau) / ssr,
        "geom_test": (pmu*pmu) / (pe*ptau),
        "log_geom": math.log(pmu*pmu) - math.log(pe*ptau),
    })

def mc_sigma(field_name):
    arr = np.array([s[field_name] for s in mc_samples])
    return arr.std()

# Build sigma map per observable
sigmas = {
    "1. r_eμ = m_μ/m_e":                sigma_ratio("mu", "e", m_mu/m_e),
    "2. r_μτ = m_τ/m_μ":                sigma_ratio("tau", "mu", m_tau/m_mu),
    "3. r_eτ = m_τ/m_e":                sigma_ratio("tau", "e", m_tau/m_e),
    "4. 1/r_eμ = m_e/m_μ":              sigma_ratio("e", "mu", m_e/m_mu),
    "5. 1/r_μτ = m_μ/m_τ":              sigma_ratio("mu", "tau", m_mu/m_tau),
    "6. 1/r_eτ = m_e/m_τ":              sigma_ratio("e", "tau", m_e/m_tau),
    "7. log(r_eμ)":                     sigma_log_ratio("mu", "e"),
    "8. log(r_μτ)":                     sigma_log_ratio("tau", "mu"),
    "9. log(r_eτ)":                     sigma_log_ratio("tau", "e"),
    # log ratio: closed form too messy; use MC
    "10. log spacing ratio":            mpf(str(np.std([math.log(s["geom_test"]/(np.exp(0))) for s in mc_samples])))*0  # placeholder, computed below
    if False else mpf("0.0"),
    "11. sqrt(m_e/m_μ)":                sigma_sqrt_ratio("e", "mu", sqrt(m_e/m_mu)),
    "12. sqrt(m_μ/m_τ)":                sigma_sqrt_ratio("mu", "tau", sqrt(m_mu/m_tau)),
    "13. sqrt(m_e/m_τ)":                sigma_sqrt_ratio("e", "tau", sqrt(m_e/m_tau)),
    "14. sqrt(m_μ/m_e)":                sigma_sqrt_ratio("mu", "e", sqrt(m_mu/m_e)),
    "15. sqrt(m_τ/m_μ)":                sigma_sqrt_ratio("tau", "mu", sqrt(m_tau/m_mu)),
    "16. Koide K = (Σm)/(Σ√m)²":        mpf(str(mc_sigma("K"))),
    "17. Koide K - 2/3":                mpf(str(mc_sigma("K_minus_2_3"))),
    "18. cone half-angle (deg)":        mpf(str(mc_sigma("angle_deg"))),
    "19. cone half-angle (rad)":        mpf(str(mc_sigma("angle_rad"))),
    "20. cone deviation from 45° (rad)":mpf(str(mc_sigma("cone_dev_rad"))),
    "21. m_e / Σm":                     mpf(str(mc_sigma("m_e_frac"))),
    "22. m_μ / Σm":                     mpf(str(mc_sigma("m_mu_frac"))),
    "23. m_τ / Σm":                     mpf(str(mc_sigma("m_tau_frac"))),
    "24. (m_μ - m_e)/(m_τ - m_e)":      mpf(str(mc_sigma("diff_frac"))),
    "25. (m_τ - m_μ)/m_τ":              mpf(str(mc_sigma("tau_minus_mu_over_tau"))),
    "26. √m_e / Σ√m":                   mpf(str(mc_sigma("sqrt_e_frac"))),
    "27. √m_μ / Σ√m":                   mpf(str(mc_sigma("sqrt_mu_frac"))),
    "28. √m_τ / Σ√m":                   mpf(str(mc_sigma("sqrt_tau_frac"))),
    "29. m_μ² / (m_e m_τ)":             mpf(str(mc_sigma("geom_test"))),
    "30. log(m_μ²) - log(m_e m_τ)":     mpf(str(mc_sigma("log_geom"))),
}

# Special case: log spacing ratio sigma via MC
log_ratio_samples = []
for _ in range(N_MC):
    pe = me_f + se_f * rng.normal()
    pmu = mmu_f + smu_f * rng.normal()
    ptau = mtau_f + stau_f * rng.normal()
    log_ratio_samples.append(math.log(ptau/pmu) / math.log(pmu/pe))
sigmas["10. log spacing ratio"] = mpf(str(np.std(log_ratio_samples)))


# ============================================================
# Search the basis
# ============================================================

def find_matches(target_value, target_sigma, forms,
                 deviation_threshold_pct=2.0):
    """Search forms for matches within deviation_threshold_pct percent of target.
    Returns matches sorted by relative deviation magnitude.
    Also tries scaled versions (factor in {1, 2, 3, 4, 5, 1/2, 1/3, 1/4, 2/3, 3/2})
    since lepton observables span O(1e-4) to O(1e3) and the basis sits in O(0.01)
    to O(2). The factor list is fixed and noted as part of the prespecification.
    """
    factors = [
        ("1",       mpf(1)),
        ("2",       mpf(2)),
        ("3",       mpf(3)),
        ("4",       mpf(4)),
        ("5",       mpf(5)),
        ("1/2",     mpf(1)/2),
        ("1/3",     mpf(1)/3),
        ("1/4",     mpf(1)/4),
        ("2/3",     mpf(2)/3),
        ("3/2",     mpf(3)/2),
    ]
    matches = []
    for fname, fval in forms.items():
        for facname, facval in factors:
            tval = facval * fval
            if tval == 0:
                continue
            if target_value == 0:
                # absolute deviation only
                continue
            diff = (tval - target_value) / target_value
            if abs(float(diff)*100) < deviation_threshold_pct:
                if target_sigma > 0:
                    n_sigma = abs(tval - target_value) / target_sigma
                else:
                    n_sigma = mp.inf
                form_label = fname if facname == "1" else f"{facname} * ({fname})"
                matches.append({
                    "form": form_label,
                    "form_base": fname,
                    "factor": facname,
                    "value": str(tval),
                    "rel_deviation_pct": float(diff)*100,
                    "within_pdg_sigma": float(n_sigma) if n_sigma != mp.inf else float('inf'),
                })
    matches.sort(key=lambda x: abs(x["rel_deviation_pct"]))
    # de-duplicate near-identical forms (e.g. "(sqrt 5 - 1)/2" appearing twice
    # under different labels) by keeping smallest deviation per (factor, base).
    seen = set()
    deduped = []
    for m in matches:
        key = m["form"]
        if key not in seen:
            seen.add(key)
            deduped.append(m)
    return deduped


print("=" * 80)
print("W3 Test 2: Lepton mass spectrum recheck with prespecified-basis methodology")
print("=" * 80)
print()
print(f"  PDG/CODATA inputs (MeV):")
print(f"    m_e   = {mp.nstr(m_e, 12)} ± {mp.nstr(sigma_e, 4)}")
print(f"    m_μ   = {mp.nstr(m_mu, 12)} ± {mp.nstr(sigma_mu, 4)}")
print(f"    m_τ   = {mp.nstr(m_tau, 12)} ± {mp.nstr(sigma_tau, 4)}")
print()
print(f"  Prespecified panel: {len(observables)} observables")
print(f"  Prespecified basis: {len(natural_forms)} forms × 10 factors = "
      f"{len(natural_forms)*10} effective tests per observable")
print()

results = {
    "input": {
        "m_e_MeV":   {"value": str(m_e),   "sigma": str(sigma_e)},
        "m_mu_MeV":  {"value": str(m_mu),  "sigma": str(sigma_mu)},
        "m_tau_MeV": {"value": str(m_tau), "sigma": str(sigma_tau)},
    },
    "n_basis_forms": len(natural_forms),
    "n_factors": 10,
    "n_observables": len(observables),
    "step_2_per_observable": {},
}

per_obs_summary = []
total_hits_1pct = 0
total_hits_05pct = 0
total_hits_within_1sigma = 0

for oname, oval in observables.items():
    osig = sigmas[oname]
    matches = find_matches(oval, osig, natural_forms,
                           deviation_threshold_pct=2.0)
    n_1pct = sum(1 for m in matches if abs(m["rel_deviation_pct"]) < 1.0)
    n_05pct = sum(1 for m in matches if abs(m["rel_deviation_pct"]) < 0.5)
    n_1sigma = sum(1 for m in matches if m["within_pdg_sigma"] < 1.0)
    total_hits_1pct += n_1pct
    total_hits_05pct += n_05pct
    total_hits_within_1sigma += n_1sigma

    print(f"--- {oname} ---")
    print(f"    value = {mp.nstr(oval, 10)}, sigma = {mp.nstr(osig, 4)}")
    print(f"    matches: <1%={n_1pct}, <0.5%={n_05pct}, within 1σ={n_1sigma}")
    for m in matches[:5]:
        flag = " *" if m["within_pdg_sigma"] < 1.0 else ""
        sig_str = (f"{m['within_pdg_sigma']:.2f}σ" if m["within_pdg_sigma"] != float('inf')
                   else "∞σ")
        print(f"       {m['form']:<55s}  diff {m['rel_deviation_pct']:+8.4f}%  "
              f"({sig_str}){flag}")
    print()
    per_obs_summary.append({
        "observable": oname,
        "value": str(oval),
        "sigma": str(osig),
        "n_1pct": n_1pct,
        "n_05pct": n_05pct,
        "n_1sigma": n_1sigma,
        "best": matches[:5],
    })
    results["step_2_per_observable"][oname] = {
        "value": str(oval),
        "sigma": str(osig),
        "n_1pct": n_1pct,
        "n_05pct": n_05pct,
        "n_1sigma": n_1sigma,
        "matches": matches,
    }

print()
print("=" * 80)
print("Aggregate statistical accounting")
print("=" * 80)

n_observables = len(observables)
n_forms = len(natural_forms)
n_factors = 10
total_tests = n_observables * n_forms * n_factors

print(f"  Total tests: {n_observables} obs × {n_forms} forms × {n_factors} factors "
      f"= {total_tests}")
print(f"  Hits at <1%   threshold: {total_hits_1pct}")
print(f"  Hits at <0.5% threshold: {total_hits_05pct}")
print(f"  Hits within 1σ measurement: {total_hits_within_1sigma}")

# ============================================================
# Statistical accounting: random expectation
# ============================================================

# For each observable, estimate the random expectation as follows:
# the "natural form basis × factors" creates 300 candidate values.
# A random target T in the same order of magnitude as the observable
# expects (300 * 2*threshold/100) ~ N_basis * 2*threshold/100 hits at threshold.

# More carefully: empirical density. For each observable, count basis
# values within an *order-of-magnitude* match window. Use log-uniform.

print()
print("Random-baseline expectation:")
np.random.seed(11)

def expected_hits_random(forms, target_order, n_trials=10000,
                         threshold_pct=1.0, factors=None):
    if factors is None:
        factors = [1.0, 2.0, 3.0, 4.0, 5.0, 0.5, 1.0/3, 0.25, 2.0/3, 1.5]
    fvals = []
    for fv in forms.values():
        ff = float(fv)
        for fac in factors:
            fvals.append(ff * fac)
    fvals = np.array(fvals)
    log_target_min = math.log(target_order * 0.1)
    log_target_max = math.log(target_order * 10)
    hits = []
    for _ in range(n_trials):
        target = math.exp(np.random.uniform(log_target_min, log_target_max))
        n_hit = sum(1 for f in fvals
                    if f != 0 and abs(f-target)/target < threshold_pct/100)
        hits.append(n_hit)
    return float(np.mean(hits)), float(np.std(hits))

# Group observables by order of magnitude; estimate hits
orders = {}
for oname, oval in observables.items():
    if abs(float(oval)) < 1e-30:
        continue
    om = abs(float(oval))
    log_om = math.log10(om)
    bucket = round(log_om)
    orders.setdefault(bucket, []).append(oname)

# Aggregate random-expectation across buckets
total_expected_1pct = 0.0
total_expected_05pct = 0.0
total_expected_var_1pct = 0.0
total_expected_var_05pct = 0.0
for bucket, names in orders.items():
    target_order = 10**bucket
    mh1, sh1 = expected_hits_random(natural_forms, target_order,
                                     n_trials=2000, threshold_pct=1.0)
    mh05, sh05 = expected_hits_random(natural_forms, target_order,
                                       n_trials=2000, threshold_pct=0.5)
    n = len(names)
    total_expected_1pct  += n * mh1
    total_expected_05pct += n * mh05
    total_expected_var_1pct  += n * (sh1**2)
    total_expected_var_05pct += n * (sh05**2)

print(f"  Observables grouped into {len(orders)} order-of-magnitude buckets")
print(f"  Expected hits @<1%   (sum over panel): {total_expected_1pct:6.2f} "
      f"± {math.sqrt(total_expected_var_1pct):.2f}")
print(f"  Observed hits @<1%:                  {total_hits_1pct}")
print(f"  Expected hits @<0.5% (sum over panel): {total_expected_05pct:6.2f} "
      f"± {math.sqrt(total_expected_var_05pct):.2f}")
print(f"  Observed hits @<0.5%:                {total_hits_05pct}")

# z-score for the deviation from random expectation
def z_score(observed, expected, std):
    if std == 0:
        return float('inf')
    return (observed - expected) / std

z_1pct = z_score(total_hits_1pct, total_expected_1pct,
                 math.sqrt(total_expected_var_1pct))
z_05pct = z_score(total_hits_05pct, total_expected_05pct,
                  math.sqrt(total_expected_var_05pct))

print(f"  Signal vs chance @<1%:   z = {z_1pct:+.2f} "
      f"({'above chance' if z_1pct > 0 else 'below or at chance'})")
print(f"  Signal vs chance @<0.5%: z = {z_05pct:+.2f}")

results["statistical_accounting"] = {
    "total_tests": total_tests,
    "total_hits_1pct": total_hits_1pct,
    "total_hits_05pct": total_hits_05pct,
    "total_hits_within_1sigma": total_hits_within_1sigma,
    "expected_hits_1pct": total_expected_1pct,
    "expected_hits_1pct_std": math.sqrt(total_expected_var_1pct),
    "expected_hits_05pct": total_expected_05pct,
    "expected_hits_05pct_std": math.sqrt(total_expected_var_05pct),
    "z_score_1pct": z_1pct,
    "z_score_05pct": z_05pct,
}

# ============================================================
# Comparison with CKM result (8 hits / 120 tests under same basis)
# ============================================================

print()
print("=" * 80)
print("Cross-comparison: lepton vs CKM result (same basis)")
print("=" * 80)

# w3_lambda_predictive_verification.py reported aggregate at the end.
# 4 parameters × 30 forms = 120 tests; 8 within-1-sigma hits
# (per CLAUDE.md §2 W3 entry).

ckm_n_tests = 4 * 30
ckm_hits_1sigma = 8

# Lepton equivalent normalized
lep_n_tests = n_observables * n_forms  # without factor multiplier — match CKM script
lep_hit_density_1sigma = total_hits_within_1sigma / lep_n_tests if lep_n_tests else 0
ckm_hit_density_1sigma = ckm_hits_1sigma / ckm_n_tests

print(f"  CKM (Test 0):    {ckm_hits_1sigma} hits within 1σ / {ckm_n_tests} tests "
      f"(density {ckm_hit_density_1sigma*100:.2f}%)")
print(f"  Leptons (Test 2): {total_hits_within_1sigma} hits within 1σ / "
      f"{n_observables * n_forms * n_factors} tests "
      f"(density {total_hits_within_1sigma/(n_observables*n_forms*n_factors)*100:.4f}%)")
print()
print(f"  Note: the lepton sigmas are very tight (CODATA m_e to 1 part in 10^9),")
print(f"  so 'within 1σ' is a much stricter bar than for CKM (where σ is ~1%).")

results["ckm_comparison"] = {
    "ckm_n_tests": ckm_n_tests,
    "ckm_hits_1sigma": ckm_hits_1sigma,
    "ckm_hit_density_1sigma_pct": ckm_hit_density_1sigma*100,
    "lepton_n_tests_w_factors": n_observables * n_forms * n_factors,
    "lepton_hits_1sigma": total_hits_within_1sigma,
    "lepton_hit_density_1sigma_pct": total_hits_within_1sigma/(n_observables*n_forms*n_factors)*100,
    "note": "lepton sigmas are O(1e-9), CKM sigmas are O(1e-2); 1-sigma bars are not directly comparable",
}

# ============================================================
# Cross-sector check: do lepton observables share forms with CKM Wolfenstein?
# ============================================================

print()
print("=" * 80)
print("Cross-sector check: lepton observable ↔ CKM Wolfenstein-form sharing")
print("=" * 80)

# CKM candidate forms identified in w3_lambda_predictive_verification.py:
ckm_candidate_forms = {
    "lambda^2 = 1/Vol(S^3)":        "1/Vol(S^3) = 1/(2 pi^2)",
    "rho_bar = 1/Vol(S^1)":         "1/Vol(S^1) = 1/(2 pi)",
    "eta_bar = ln(2)/2":            "ln(2)/2",
    "A^2 = ln(2)":                  "ln(2)",
}

ckm_form_keys = set(ckm_candidate_forms.values())

cross_sector_hits = []
for oname, info in results["step_2_per_observable"].items():
    for m in info["matches"]:
        if abs(m["rel_deviation_pct"]) < 1.0 and m["form_base"] in ckm_form_keys:
            ckm_param = next(k for k, v in ckm_candidate_forms.items()
                             if v == m["form_base"])
            cross_sector_hits.append({
                "lepton_observable": oname,
                "lepton_value": info["value"],
                "shared_form": m["form"],
                "shared_form_base": m["form_base"],
                "ckm_parameter": ckm_param,
                "rel_deviation_pct": m["rel_deviation_pct"],
                "within_lepton_sigma": m["within_pdg_sigma"],
            })

if cross_sector_hits:
    print(f"  Found {len(cross_sector_hits)} cross-sector form-shares at <1%:")
    print()
    for cs in cross_sector_hits[:30]:
        sig_str = (f"{cs['within_lepton_sigma']:.2f}σ"
                   if cs['within_lepton_sigma'] != float('inf') else "∞σ")
        print(f"    {cs['lepton_observable']:<40s}")
        print(f"      ↔ {cs['ckm_parameter']:<22s}  via {cs['shared_form']:<40s}")
        print(f"      rel diff {cs['rel_deviation_pct']:+8.4f}%  ({sig_str})")
        print()
else:
    print("  No cross-sector form-shares at <1%.")

results["cross_sector_hits"] = cross_sector_hits

# ============================================================
# Save
# ============================================================

print()
print("=" * 80)
print("Saving results to debug/data/w3_lepton_mass_recheck.json")
print("=" * 80)


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
    if isinstance(obj, float):
        if obj == float('inf'):
            return "inf"
        if obj == float('-inf'):
            return "-inf"
    return obj


with open(OUT_DIR / "w3_lepton_mass_recheck.json", "w") as f:
    json.dump(to_json_safe(results), f, indent=2)

print()
print("Done. Summary:")
print(f"  Observables tested: {n_observables}")
print(f"  Basis forms × factors: {n_forms} × {n_factors} = {n_forms*n_factors}")
print(f"  Hits at <1%:    {total_hits_1pct}")
print(f"  Hits at <0.5%:  {total_hits_05pct}")
print(f"  Hits within 1σ: {total_hits_within_1sigma}")
print(f"  Cross-sector form-sharing hits: {len(cross_sector_hits)}")
print(f"  z-score vs random @<1%:   {z_1pct:+.2f}")
print(f"  z-score vs random @<0.5%: {z_05pct:+.2f}")
