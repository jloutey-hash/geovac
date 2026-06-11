"""
Sprint W3-CKM-PMNS: information-theoretic / packing-style probe of the CKM
and PMNS mixing matrices.

Question: do quark mixing (CKM) and neutrino mixing (PMNS) admit
packing-axiom-style structural constraints, the way the charged lepton
mass triple admits the Koide 45-deg cone constraint?

Universality logic:
  - Lepton masses: Koide K = 2/3 holds essentially exactly (1 arcsec on cone).
  - Quark masses: Koide K is 10-27% off; the 45-deg cone breaks for quarks.
  - Asymmetry maps to Connes inner factor split: (C + H) lepton vs M_3(C) quark.

If the quark sector has its OWN packing constraint -- different from Koide
because it lives in M_3(C) instead of (C + H) -- it should show up as a
geometric fact in the CKM matrix. Same for PMNS in the lepton sector.

Probes:
  1. Wolfenstein hierarchy: do |V_us|, |V_cb|, |V_ub| sit at geometric series?
  2. Cabibbo and other CKM angles: closed-form / natural-constant matches.
  3. Koide-style cone constraints on three CKM-derived triples.
  4. Same probes for PMNS magnitudes and angles.
  5. Jarlskog invariant: structural form?
  6. GeoVac-internal cross-references.
  7. PSLQ on key invariants.
"""

import json
import math
from pathlib import Path

import numpy as np
import mpmath as mp
from mpmath import mpf, pi, log, sqrt, atan, atan2, sin, cos, exp


mp.mp.dps = 80
OUT_DIR = Path("debug/data")

# ============================================================
# Data: CKM magnitudes (PDG 2024 central values)
# ============================================================

V_CKM_mag = np.array([
    [0.97435, 0.22500, 0.00369],   # u row
    [0.22486, 0.97349, 0.04182],   # c row
    [0.00857, 0.04110, 0.999118],  # t row
])

# Wolfenstein parameters (PDG 2024)
lam     = mpf("0.22500")
A_w     = mpf("0.826")
rho_bar = mpf("0.159")
eta_bar = mpf("0.348")

# CKM CP phase (rad), from rho_bar, eta_bar
delta_ckm = atan2(eta_bar, rho_bar)

# ============================================================
# Data: PMNS magnitudes from PDG 2024 best-fit angles
# ============================================================

sin2_t12 = mpf("0.307")
sin2_t23 = mpf("0.546")
sin2_t13 = mpf("0.0220")
delta_cp_pmns_deg = mpf("230")

s12 = sqrt(sin2_t12); c12 = sqrt(1 - sin2_t12)
s23 = sqrt(sin2_t23); c23 = sqrt(1 - sin2_t23)
s13 = sqrt(sin2_t13); c13 = sqrt(1 - sin2_t13)
delta_cp_pmns_rad = delta_cp_pmns_deg * pi / 180

# Standard PDG parameterization:
# V = R23(theta23) * U_delta * R13(theta13) * U_delta^dagger * R12(theta12)
# Magnitudes:
V_PMNS_mag = np.array([
    [float(c12*c13),                                         float(s12*c13),                                         float(s13)],
    [float(sqrt((s12*c23)**2 + (c12*s23*s13)**2 - 2*s12*c23*c12*s23*s13*cos(delta_cp_pmns_rad))),
     float(sqrt((c12*c23)**2 + (s12*s23*s13)**2 + 2*c12*c23*s12*s23*s13*cos(delta_cp_pmns_rad))),
     float(s23*c13)],
    [float(sqrt((s12*s23)**2 + (c12*c23*s13)**2 + 2*s12*s23*c12*c23*s13*cos(delta_cp_pmns_rad))),
     float(sqrt((c12*s23)**2 + (s12*c23*s13)**2 - 2*c12*s23*s12*c23*s13*cos(delta_cp_pmns_rad))),
     float(c23*c13)],
])

# Sanity: rows and columns should sum to 1 (squared)
for i in range(3):
    row_sq = sum(V_PMNS_mag[i, j]**2 for j in range(3))
    col_sq = sum(V_PMNS_mag[j, i]**2 for j in range(3))
    assert abs(row_sq - 1) < 5e-3, f"PMNS row {i}: {row_sq}"
    assert abs(col_sq - 1) < 5e-3, f"PMNS col {i}: {col_sq}"

results = {}


def pretty(x, n=12):
    return mp.nstr(x, n)


# ============================================================
# PROBE 1: Wolfenstein hierarchy
# ============================================================

print("=" * 72)
print("PROBE 1: Wolfenstein hierarchy on CKM off-diagonals")
print("=" * 72)

V_us = mpf(str(V_CKM_mag[0, 1]))
V_cb = mpf(str(V_CKM_mag[1, 2]))
V_ub = mpf(str(V_CKM_mag[0, 2]))

log_V_us = log(V_us)
log_V_cb = log(V_cb)
log_V_ub = log(V_ub)

# log spacings
spacing_1 = log_V_us - log_V_cb   # log(|V_us| / |V_cb|)
spacing_2 = log_V_cb - log_V_ub   # log(|V_cb| / |V_ub|)
ratio_of_spacings = spacing_2 / spacing_1

print(f"  |V_us| = {pretty(V_us, 8)}    (Wolfenstein order: lambda^1)")
print(f"  |V_cb| = {pretty(V_cb, 8)}    (lambda^2 expected: {pretty(lam**2, 6)})")
print(f"  |V_ub| = {pretty(V_ub, 8)}    (lambda^3 expected: {pretty(lam**3, 6)})")
print()
print(f"  log(|V_us|/|V_cb|) = {pretty(spacing_1, 8)}")
print(f"  log(|V_cb|/|V_ub|) = {pretty(spacing_2, 8)}")
print(f"  ratio of spacings = {pretty(ratio_of_spacings, 8)}")
print()
# If purely geometric series: ratio = 1
geom_test = abs(ratio_of_spacings - 1)
print(f"  Pure geometric series (ratio = 1)? deviation: {pretty(geom_test, 4)} ({float(geom_test*100):.2f}%)")
# Check ratio against 1 (geometric), 2 (each step doubles), small rationals
for q in range(1, 8):
    for p in range(1, 4*q):
        if math.gcd(p, q) > 1:
            continue
        diff = abs(ratio_of_spacings - mpf(p)/q)
        if float(diff) < 0.05:
            print(f"     Close to {p}/{q} = {pretty(mpf(p)/q, 6)}: diff {pretty(diff, 4)}")

results["wolfenstein_hierarchy"] = {
    "V_us": str(V_us),
    "V_cb": str(V_cb),
    "V_ub": str(V_ub),
    "lambda_squared": str(lam**2),
    "lambda_cubed": str(lam**3),
    "log_spacing_1": str(spacing_1),
    "log_spacing_2": str(spacing_2),
    "ratio_of_spacings": str(ratio_of_spacings),
    "deviation_from_pure_geom": str(geom_test),
}

# ============================================================
# PROBE 2: Cabibbo angle and other quark mixing angles
# ============================================================

print()
print("=" * 72)
print("PROBE 2: CKM angles vs natural / framework constants")
print("=" * 72)

theta_C_rad = atan(V_us / mpf(str(V_CKM_mag[0, 0])))
theta_C_deg = theta_C_rad * 180 / pi

# Known relation: tan(theta_C) ~ sqrt(m_d / m_s)  [Gatto-Sartori-Tonin]
m_d = mpf("4.70")
m_s = mpf("93.5")
sqrt_md_ms = sqrt(m_d / m_s)
gst_diff = (V_us / mpf(str(V_CKM_mag[0,0])) - sqrt_md_ms) / sqrt_md_ms

# CKM angles via standard parameterization
# theta_12 ~ Cabibbo, theta_23 ~ |V_cb|, theta_13 ~ |V_ub|
theta12_ckm = theta_C_deg
theta23_ckm = mp.asin(V_cb) * 180 / pi
theta13_ckm = mp.asin(V_ub) * 180 / pi

print(f"  theta_12 (Cabibbo)    = {pretty(theta_C_deg, 8)} deg")
print(f"  theta_23              = {pretty(theta23_ckm, 8)} deg")
print(f"  theta_13              = {pretty(theta13_ckm, 8)} deg")
print(f"  delta_CP (Wolfenstein) = {pretty(delta_ckm * 180/pi, 6)} deg")
print()
print(f"  Gatto-Sartori-Tonin: tan(theta_C) ?= sqrt(m_d/m_s)")
print(f"    tan(theta_C)    = {pretty(V_us/mpf(str(V_CKM_mag[0,0])), 8)}")
print(f"    sqrt(m_d/m_s)   = {pretty(sqrt_md_ms, 8)}")
print(f"    relative diff   = {pretty(gst_diff*100, 4)}%")
print()

# Try to identify theta_C in degrees as a simple combination
# Tested constants: 13 (closest integer), 12+something, sqrt-related forms
print(f"  Closest small rationals for theta_C / 1 deg = {pretty(theta_C_deg, 6)}:")
for q in range(1, 30):
    for p in range(8*q, 16*q):
        if math.gcd(p, q) > 1:
            continue
        diff = abs(theta_C_deg - mpf(p)/q)
        if float(diff) < 0.1:
            print(f"     {p}/{q} = {pretty(mpf(p)/q, 6)}, diff {pretty(diff, 4)}")

# theta_C in rad against natural constants:
# 1/137 ~ 0.0073 (alpha)
# sqrt(alpha) ~ 0.0854
# pi/alpha ?
# alpha * 137 = 1 (definition)
print()
print(f"  theta_C / pi = {pretty(theta_C_rad/pi, 8)}")
print(f"  theta_C * 137 / pi (ie theta_C / (pi * alpha)) = {pretty(theta_C_rad / (pi * mpf('1')/mpf('137.036')), 6)}")

results["ckm_angles"] = {
    "theta_C_deg":       str(theta_C_deg),
    "theta_23_deg":      str(theta23_ckm),
    "theta_13_deg":      str(theta13_ckm),
    "delta_CP_deg":      str(delta_ckm*180/pi),
    "gst_relation": {
        "tan_theta_C":   str(V_us/mpf(str(V_CKM_mag[0,0]))),
        "sqrt_md_ms":    str(sqrt_md_ms),
        "rel_diff":      str(gst_diff),
    },
}

# ============================================================
# PROBE 3: Koide-style cone tests on CKM-derived triples
# ============================================================

print()
print("=" * 72)
print("PROBE 3: Koide-style cone constraints on CKM triples")
print("=" * 72)


def koide_cone(triple, label):
    """For three positive numbers x_i, compute the Koide K and the angle
    of the sqrt-vector from the democratic direction (1,1,1)/sqrt(3)."""
    x = [mpf(str(v)) for v in triple]
    sx = sum(x)
    rsx = sum(sqrt(xi) for xi in x)
    K = sx / (rsx * rsx)
    z = [sqrt(xi) for xi in x]
    z_norm_sq = sum(zi*zi for zi in z)
    z_dot_d = sum(z) / sqrt(3)
    cos_sq = (z_dot_d * z_dot_d) / z_norm_sq
    # cos_sq must be in [1/3, 1] always for positive triples
    if cos_sq > 1:
        cos_sq = mpf(1)
    angle_deg = mp.acos(sqrt(cos_sq)) * 180 / pi
    print(f"  {label}")
    print(f"    triple = {[pretty(xi, 6) for xi in x]}")
    print(f"    K (Koide-like)        = {pretty(K, 8)}")
    print(f"    angle from (1,1,1)    = {pretty(angle_deg, 8)} deg")
    print(f"    deviation from 45 deg = {pretty(angle_deg - 45, 4)} deg")
    return {
        "triple": [str(xi) for xi in x],
        "K": str(K),
        "angle_deg": str(angle_deg),
        "dev_from_45_deg": str(angle_deg - 45),
    }


# Triple A: diagonal CKM magnitudes (mostly ~1)
triple_A = [V_CKM_mag[i, i] for i in range(3)]
res_A = koide_cone(triple_A, "Triple A: diagonal CKM magnitudes (|V_ud|, |V_cs|, |V_tb|)")

# Triple B: off-diagonal CKM hierarchy
triple_B = [V_CKM_mag[0, 1], V_CKM_mag[1, 2], V_CKM_mag[0, 2]]
res_B = koide_cone(triple_B, "Triple B: |V_us|, |V_cb|, |V_ub| (the Wolfenstein hierarchy)")

# Triple C: lower-triangular off-diagonals
triple_C = [V_CKM_mag[1, 0], V_CKM_mag[2, 1], V_CKM_mag[2, 0]]
res_C = koide_cone(triple_C, "Triple C: |V_cd|, |V_ts|, |V_td|")

# Triple D: 1 - diagonal magnitudes (deviation from identity)
triple_D = [1 - V_CKM_mag[i, i] for i in range(3)]
res_D = koide_cone(triple_D, "Triple D: 1 - |V_ii| (mixing strength per row)")

# Triple E: square magnitudes of diagonal -- "stay-prob" probabilities
triple_E = [V_CKM_mag[i, i]**2 for i in range(3)]
res_E = koide_cone(triple_E, "Triple E: |V_ii|^2  (no-flavor-change probabilities)")

# Triple F: sin^2 of the three CKM mixing angles (PDG-like)
triple_F = [
    (V_us / mpf(str(V_CKM_mag[0,0])))**2 / (1 + (V_us / mpf(str(V_CKM_mag[0,0])))**2),  # sin^2(theta_12)
    V_cb**2,                                                                              # ~sin^2(theta_23)
    V_ub**2,                                                                              # ~sin^2(theta_13)
]
res_F = koide_cone(triple_F, "Triple F: sin^2 of CKM angles")

results["ckm_koide_cones"] = {
    "diagonal_magnitudes":      res_A,
    "Wolfenstein_hierarchy":    res_B,
    "lower_offdiag":            res_C,
    "mixing_strength":          res_D,
    "stay_probabilities":       res_E,
    "sin_squared_angles":       res_F,
}

# ============================================================
# PROBE 4: Same Koide-style probes on PMNS
# ============================================================

print()
print("=" * 72)
print("PROBE 4: Koide-style probes on PMNS")
print("=" * 72)

# Print PMNS for sanity
print(f"  PMNS magnitude matrix:")
for i in range(3):
    print(f"    {V_PMNS_mag[i,0]:.4f}  {V_PMNS_mag[i,1]:.4f}  {V_PMNS_mag[i,2]:.4f}")
print()

# Same triples as for CKM
pmns_results = {}
pmns_results["diagonal"] = koide_cone(
    [V_PMNS_mag[i, i] for i in range(3)],
    "Triple A: diagonal PMNS")
pmns_results["hierarchy_top"] = koide_cone(
    [V_PMNS_mag[0, 1], V_PMNS_mag[1, 2], V_PMNS_mag[0, 2]],
    "Triple B: |U_e2|, |U_mu3|, |U_e3|")
pmns_results["lower_offdiag"] = koide_cone(
    [V_PMNS_mag[1, 0], V_PMNS_mag[2, 1], V_PMNS_mag[2, 0]],
    "Triple C: |U_mu1|, |U_tau2|, |U_tau1|")
pmns_results["sin_squared_angles"] = koide_cone(
    [float(sin2_t12), float(sin2_t23), float(sin2_t13)],
    "Triple F: sin^2 of PMNS angles")
pmns_results["mixing_angles_radians"] = koide_cone(
    [float(mp.asin(s12)), float(mp.asin(s23)), float(mp.asin(s13))],
    "Triple G: PMNS mixing angles in radians")

results["pmns_koide_cones"] = pmns_results

# Highlight: PMNS has nearly tri-bimaximal pattern
# Tri-bimaximal: sin^2(theta_12) = 1/3, sin^2(theta_23) = 1/2, sin^2(theta_13) = 0
print()
print("  Tri-bimaximal benchmark:")
print(f"    sin^2(theta_12) = {sin2_t12} (TB: 1/3 = 0.333) diff = {float(sin2_t12)-1/3:+.4f}")
print(f"    sin^2(theta_23) = {sin2_t23} (TB: 1/2 = 0.500) diff = {float(sin2_t23)-1/2:+.4f}")
print(f"    sin^2(theta_13) = {sin2_t13} (TB: 0)         diff = {float(sin2_t13):+.4f}")

# ============================================================
# PROBE 5: Jarlskog invariant
# ============================================================

print()
print("=" * 72)
print("PROBE 5: Jarlskog invariants (CP-violation measure)")
print("=" * 72)

# CKM Jarlskog from Wolfenstein: J = A^2 lambda^6 eta_bar
J_ckm_wolf = A_w**2 * lam**6 * eta_bar
print(f"  J_CKM (Wolfenstein) = A^2 lambda^6 eta_bar = {pretty(J_ckm_wolf, 6)}")
print(f"  PDG measured J_CKM ~ 3.18e-5")

# PMNS Jarlskog
# J_PMNS = c12 c23 c13^2 s12 s23 s13 sin(delta_cp)
J_pmns = c12 * c23 * c13**2 * s12 * s23 * s13 * abs(sin(delta_cp_pmns_rad))
print(f"  J_PMNS (PDG best fit, NH) = {pretty(J_pmns, 6)}")
print(f"  Ratio J_PMNS / J_CKM = {pretty(J_pmns/J_ckm_wolf, 6)}")

# Is J_PMNS / J_CKM close to a simple ratio?
ratio_J = J_pmns / J_ckm_wolf
print(f"  log10(J_PMNS / J_CKM) = {pretty(mp.log10(ratio_J), 4)}")

results["jarlskog"] = {
    "J_CKM_wolfenstein": str(J_ckm_wolf),
    "J_PMNS_PDG":         str(J_pmns),
    "ratio":              str(J_pmns / J_ckm_wolf),
    "log10_ratio":        str(mp.log10(ratio_J)),
}

# ============================================================
# PROBE 6: Wolfenstein lambda vs natural / GeoVac constants
# ============================================================

print()
print("=" * 72)
print("PROBE 6: Wolfenstein lambda and Cabibbo angle vs constants")
print("=" * 72)

# lambda has been speculated to be tied to:
#   - sqrt(2)/sqrt(40)
#   - 1 / (something)
#   - alpha-related forms
constants = {
    "1": mpf(1),
    "1/(4 pi)":             1/(4*pi),
    "1/sqrt(20)":           1/sqrt(20),
    "1/sqrt(Delta_inv)":    1/sqrt(40),
    "alpha^(1/2)":          sqrt(1/mpf("137.036")),
    "1/sqrt(2*pi^2)":       1/sqrt(2*pi*pi),
    "(B Delta)^(1/2)":      sqrt(mpf(42)/mpf(40)),
    "log(2)/pi":            log(2)/pi,
    "1/pi":                 1/pi,
    "1/sqrt(20)":           1/sqrt(20),
    "1/(2 sqrt(5))":        1/(2*sqrt(5)),
}

print(f"  Wolfenstein lambda = {pretty(lam, 8)}")
print(f"  {'Candidate':<28s} {'Value':>14s} {'Diff %':>10s}")
print(f"  {'-'*52}")
matches = []
for name, val in constants.items():
    diff = (lam - val) / lam
    print(f"  {name:<28s} {pretty(val, 8):>14s} {float(diff)*100:>9.3f}%")
    if abs(float(diff)) < 0.05:
        matches.append({"name": name, "value": str(val), "rel_diff": str(diff)})

# Check: theta_C = arctan(lambda) deg vs simple forms
# theta_C ~ 13.04 deg
print()
print(f"  Cabibbo angle theta_C = {pretty(theta_C_deg, 8)} deg")
print(f"  theta_C in radians: {pretty(theta_C_rad, 8)}")
print(f"  pi/24 = 7.5 deg, pi/14 ~ 12.86 deg, pi/13.8 ~ 13.04 deg")
print(f"  theta_C * 14 / pi = {pretty(theta_C_rad * 14 / pi, 8)} (=1 if theta_C = pi/14)")
print(f"  theta_C * 27 / (2 pi) = {pretty(theta_C_rad * 27 / (2*pi), 8)} (=1 if theta_C = 2pi/27)")

results["lambda_constants_test"] = {
    "lambda": str(lam),
    "theta_C_deg": str(theta_C_deg),
    "theta_C_rad": str(theta_C_rad),
    "matches_within_5pct": matches,
}

# ============================================================
# PROBE 7: GeoVac internal cross-references
# ============================================================

print()
print("=" * 72)
print("PROBE 7: GeoVac internal cross-references (CKM/PMNS invariants)")
print("=" * 72)

gv_numbers = {
    "kappa = -1/16":            mpf("-1")/16,
    "B = 42":                   mpf(42),
    "F = pi^2/6":               pi*pi/6,
    "Delta = 1/40":             mpf(1)/40,
    "1/Delta = 40":             mpf(40),
    "K_alpha":                  pi*(mpf(42) + pi*pi/6 - mpf(1)/40),
    "alpha_inv = 137.036":      mpf("137.035999084"),
    "alpha":                    1/mpf("137.035999084"),
    "Vol(S^3) = 2 pi^2":        2*pi*pi,
    "Vol(S^2) = 4 pi":          4*pi,
    "Hopf measure pi":          pi,
    "d_max = 4":                mpf(4),
    "1/4pi (S^2 Weyl)":         1/(4*pi),
    "Casimir m=3 (B):42":       mpf(42),
    "g_3^Dirac = 40":           mpf(40),
}

invariants = {
    "lambda_Wolfenstein":               lam,
    "lambda^2":                          lam**2,
    "lambda^3":                          lam**3,
    "Cabibbo theta_C (rad)":             theta_C_rad,
    "Cabibbo theta_C (deg)":             theta_C_deg,
    "|V_cb|":                            V_cb,
    "|V_ub|":                            V_ub,
    "Jarlskog CKM":                      J_ckm_wolf,
    "Jarlskog PMNS":                     J_pmns,
    "PMNS theta_12 (deg)":               mp.asin(s12) * 180/pi,
    "PMNS theta_23 (deg)":               mp.asin(s23) * 180/pi,
    "PMNS theta_13 (deg)":               mp.asin(s13) * 180/pi,
    "ratio log spacings (CKM)":          ratio_of_spacings,
}

print(f"  Looking for sub-1% non-trivial matches...")
print()
near_misses = []
for iname, ival in invariants.items():
    for gname, gval in gv_numbers.items():
        if gval == 0 or ival == 0:
            continue
        ratio = ival / gval
        for fac_name, fac_val in [
            ("1", mpf(1)), ("2", mpf(2)), ("3", mpf(3)), ("4", mpf(4)), ("5", mpf(5)),
            ("1/2", mpf(1)/2), ("1/3", mpf(1)/3), ("1/4", mpf(1)/4), ("1/5", mpf(1)/5),
            ("2/3", mpf(2)/3), ("3/2", mpf(3)/2), ("3/4", mpf(3)/4), ("4/3", mpf(4)/3),
            ("pi", pi), ("1/pi", 1/pi), ("sqrt(2)", sqrt(2)), ("sqrt(3)", sqrt(3)),
        ]:
            test = ratio / fac_val
            diff = abs(test - 1)
            if diff < mpf("0.01"):
                rel_pct = float(diff)*100
                near_misses.append({
                    "invariant": iname,
                    "invariant_value": str(ival),
                    "geovac_number": gname,
                    "geovac_value": str(gval),
                    "factor": fac_name,
                    "rel_diff_pct": rel_pct,
                })

# Print only sub-0.5% non-trivial
shown = 0
for nm in sorted(near_misses, key=lambda x: x["rel_diff_pct"]):
    if nm["rel_diff_pct"] < 0.5:
        print(f"  {nm['invariant']:<28s} ~ {nm['factor']:>10s} * {nm['geovac_number']:<24s} (off {nm['rel_diff_pct']:.4f}%)")
        shown += 1
if shown == 0:
    print("  (no sub-0.5% non-trivial matches)")

results["geovac_cross_refs"] = {
    "near_misses": near_misses,
    "criterion": "sub-1% match for invariant = factor * GeoVac_number",
}

# ============================================================
# PROBE 8: PSLQ on Wolfenstein lambda and Cabibbo
# ============================================================

print()
print("=" * 72)
print("PROBE 8: PSLQ -- is lambda a closed form in natural constants?")
print("=" * 72)

basis = [pi, 1/pi, log(2), sqrt(2), sqrt(3), sqrt(5), mp.euler, mpf(1)]
basis_names = ["pi", "1/pi", "log(2)", "sqrt(2)", "sqrt(3)", "sqrt(5)", "Euler", "1"]

for tname, tval in [("lambda", lam), ("theta_C_rad", theta_C_rad), ("ratio_spacings", ratio_of_spacings)]:
    print(f"  PSLQ on {tname} = {pretty(tval, 30)}")
    try:
        rel = mp.pslq([tval] + basis, tol=mpf("1e-15"), maxcoeff=200)
        if rel and rel[0] != 0:
            print(f"    relation: target_coef={rel[0]}, basis_coefs={list(rel[1:])}")
            terms = [f"{rel[i+1]:+d}*{basis_names[i]}" for i in range(len(basis)) if rel[i+1] != 0]
            print(f"    => {rel[0]} * {tname} {' '.join(terms)} = 0")
        else:
            print(f"    no relation found within tol/maxcoeff")
    except Exception as ex:
        print(f"    PSLQ failed: {ex}")

# ============================================================
# Save
# ============================================================

print()
print("=" * 72)
print("Done. Saving results to debug/data/w3_ckm_pmns_probe.json")
print("=" * 72)

def to_json_safe(obj):
    if isinstance(obj, mp.mpf):
        return str(obj)
    if isinstance(obj, dict):
        return {k: to_json_safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [to_json_safe(v) for v in obj]
    if isinstance(obj, tuple):
        return [to_json_safe(v) for v in obj]
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    return obj

with open(OUT_DIR / "w3_ckm_pmns_probe.json", "w") as f:
    json.dump(to_json_safe(results), f, indent=2)
