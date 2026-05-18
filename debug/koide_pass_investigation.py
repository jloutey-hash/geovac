"""
Koide cone focused investigation (2026-05-18).

Read-only, no production code modifications. Builds on Sprint W3 (May 2026)
but is Koide-specific (W3 tested general lepton mass spectrum, not the
Koide cone structure directly).

Part 1: Numerical verification of Koide at PDG 2024 precision.
Part 2: Focused PSLQ on Koide-specific quantities against a Koide-motivated
        basis (master Mellin engine M1/M2/M3 + SU(3) Casimir rationals +
        framework algebraics).
Part 3: Structural angles — k ∈ {0,1,2} mapping, SU(3) 45° check, inner-factor
        Mellin engine theorem implications.
Part 4: Honest verdict.

Output:
    debug/data/koide_pass_investigation.json
"""

import json
from pathlib import Path

import mpmath as mp
from mpmath import mpf, mpc, pi, log, sqrt, exp, atan2, acos, cos, sin, zeta


mp.mp.dps = 100

OUT_DIR = Path("debug/data")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# Part 1: PDG 2024 lepton masses
# ============================================================

# PDG 2024 charged lepton masses (MeV)
m_e = mpf("0.51099895000")
m_mu = mpf("105.6583755")
m_tau = mpf("1776.86")

sigma_e = mpf("0.00000000015")
sigma_mu = mpf("0.0000023")
sigma_tau = mpf("0.12")  # dominant uncertainty

masses = [m_e, m_mu, m_tau]
sqrt_masses = [sqrt(m) for m in masses]

# Koide ratio
K = sum(masses) / (sum(sqrt_masses)) ** 2
K_target = mpf("2") / 3
K_dev_abs = K - K_target
# Propagate tau uncertainty (dominant)
# K = S1 / S_half^2, dK/dm_tau computed via finite difference
m_tau_plus = m_tau + sigma_tau
masses_plus = [m_e, m_mu, m_tau_plus]
sqrt_masses_plus = [sqrt(m) for m in masses_plus]
K_plus = sum(masses_plus) / (sum(sqrt_masses_plus)) ** 2
dK_dmtau = (K_plus - K) / sigma_tau
# Linear error propagation
sigma_K_from_tau = abs(dK_dmtau * sigma_tau)
# Total sigma — tau dominates by orders of magnitude
sigma_K_total = sigma_K_from_tau
K_dev_sigma = K_dev_abs / sigma_K_total

# Angle from democratic (1,1,1) direction
democratic = [mpf(1), mpf(1), mpf(1)]
dot = sum(sqrt_masses[i] * democratic[i] for i in range(3))
norm_sqrt = sqrt(sum(s * s for s in sqrt_masses))
norm_dem = sqrt(mpf(3))
cos_theta = dot / (norm_sqrt * norm_dem)
theta_rad = acos(cos_theta)
theta_deg = theta_rad * 180 / pi
theta_45 = pi / 4
theta_45_deg = mpf(45)
deviation_arcsec = (theta_deg - mpf(45)) * 3600
cos2_theta = cos_theta * cos_theta

# Also compute cos^2 propagated uncertainty (for the 1/2 test)
masses_plus2 = [m_e, m_mu, m_tau_plus]
sqrt_masses_plus2 = [sqrt(m) for m in masses_plus2]
dot_p = sum(sqrt_masses_plus2[i] for i in range(3))
norm_p = sqrt(sum(s * s for s in sqrt_masses_plus2))
cos_theta_p = dot_p / (norm_p * sqrt(mpf(3)))
sigma_cos2 = abs(cos_theta_p * cos_theta_p - cos2_theta)

part1 = {
    "PDG_2024": {
        "m_e_MeV": str(m_e),
        "m_mu_MeV": str(m_mu),
        "m_tau_MeV": str(m_tau),
        "sigma_tau_MeV": str(sigma_tau),
    },
    "Koide_ratio_K": mp.nstr(K, 30),
    "Koide_target_2over3": mp.nstr(K_target, 30),
    "K_minus_2over3": mp.nstr(K_dev_abs, 15),
    "sigma_K_from_tau": mp.nstr(sigma_K_total, 15),
    "K_deviation_sigma": mp.nstr(K_dev_sigma, 6),
    "cos_theta": mp.nstr(cos_theta, 30),
    "cos2_theta": mp.nstr(cos2_theta, 30),
    "cos2_minus_half": mp.nstr(cos2_theta - mpf("0.5"), 15),
    "theta_deg": mp.nstr(theta_deg, 15),
    "theta_minus_45_arcsec": mp.nstr(deviation_arcsec, 8),
    "sigma_cos2": mp.nstr(sigma_cos2, 15),
}

print("=" * 60)
print("Part 1: Numerical verification")
print("=" * 60)
for k, v in part1.items():
    print(f"  {k}: {v}")

# ============================================================
# Part 2: Focused PSLQ on Koide-specific quantities
# ============================================================

# Koide-specific targets
r1 = sqrt(m_e / m_mu)
r2 = sqrt(m_e / m_tau)
r3 = sqrt(m_mu / m_tau)

# Also the angle and cos^2
target_panel = {
    "r1_sqrt_me_over_mmu": r1,
    "r2_sqrt_me_over_mtau": r2,
    "r3_sqrt_mmu_over_mtau": r3,
    "theta_radians": theta_rad,
    "cos2_theta": cos2_theta,
    "K_minus_2over3": K_dev_abs,
    "theta_minus_pi_over_4": theta_rad - pi / 4,
}

# ============================================================
# Build Koide-motivated basis
# ============================================================
# CRITICAL: This basis is HAND-CURATED but motivated by structural candidates,
# not by what fits the data. We document each entry's source.

# (a) Master Mellin engine M1 (Hopf-base measure family)
M1_basis = {
    "1/(2*pi^2)": 1 / (2 * pi * pi),   # Vol(S^3) reciprocal
    "1/pi": 1 / pi,
    "4/pi": 4 / pi,
    "1/(4*pi)": 1 / (4 * pi),          # Vol(S^2) reciprocal
    "pi/16": pi / 16,
    "pi/12": pi / 12,
    "pi/8": pi / 8,
    "pi/6": pi / 6,
    "pi/4": pi / 4,
    "pi/3": pi / 3,
    "pi/2": pi / 2,
    "pi": pi,
    "pi^2/6": pi ** 2 / 6,
    "pi^2": pi ** 2,
}

# (b) Master Mellin engine M2 (Seeley-DeWitt / heat-kernel)
M2_basis = {
    "sqrt(pi)": sqrt(pi),
    "sqrt(pi)/2": sqrt(pi) / 2,
    "1/sqrt(pi)": 1 / sqrt(pi),
    "zeta(2)": zeta(2),
    "zeta(4)": zeta(4),
    "zeta(3)": zeta(3),  # included as control
    "log(2)": log(2),
    "log(2)^2": log(2) ** 2,
}

# (c) Master Mellin engine M3 (vertex parity Hurwitz / Dirichlet-L)
# Catalan G = beta(2), beta(4)
G = mp.catalan
M3_basis = {
    "Catalan_G": G,
    "G^2": G * G,
    "1/G": 1 / G,
    "beta(2)_equals_G": G,  # same as Catalan
}

# (d) Three-generation structural candidates: SU(3) Casimir rationals
SU3_basis = {
    "1/3": mpf(1) / 3,
    "2/3": mpf(2) / 3,
    "1/2": mpf(1) / 2,
    "1/6": mpf(1) / 6,
    "1/8": mpf(1) / 8,   # octet dim
    "1/9": mpf(1) / 9,
    "1/10": mpf(1) / 10,  # decuplet dim
    "8/9": mpf(8) / 9,
    "3/8": mpf(3) / 8,
    "1/27": mpf(1) / 27,
    "27": mpf(27),
}

# (e) Framework-internal small algebraics (no Yukawas — those are inner factor)
FRAMEWORK_basis = {
    "kappa_squared_1/256": mpf(1) / 256,
    "Delta_1/40": mpf(1) / 40,
    "1/B_1/42": mpf(1) / 42,
    "F_pi2_6": pi ** 2 / 6,
    "alpha": mpf("0.0072973525693"),  # fine structure
    "1/alpha": 1 / mpf("0.0072973525693"),
}

# (f) Golden ratio (as control for selection-bias check)
phi = (1 + sqrt(5)) / 2
GOLDEN_basis = {
    "phi": phi,
    "1/phi": 1 / phi,
    "phi^2": phi * phi,
    "phi-1": phi - 1,
}

# All bases combined
ALL_BASES = {
    "M1_Hopf": M1_basis,
    "M2_HeatKernel": M2_basis,
    "M3_Hurwitz": M3_basis,
    "SU3_Casimir": SU3_basis,
    "Framework": FRAMEWORK_basis,
    "Golden": GOLDEN_basis,
}

# Small integer factors (consistent with W3 methodology)
factors = [mpf(1), mpf(2), mpf(3), mpf(4), mpf(5),
           mpf(1) / 2, mpf(1) / 3, mpf(1) / 4, mpf(1) / 5, mpf(1) / 6,
           mpf(2) / 3, mpf(3) / 2, mpf(3) / 4, mpf(4) / 3]

# Total basis size for selection-bias accounting
total_basis_size = sum(len(b) for b in ALL_BASES.values())
total_tests = total_basis_size * len(factors)


def search_target(target, target_name, tolerance_relative=1e-3, ceiling=10000):
    """Search Koide-motivated basis for matches to target."""
    target_abs = abs(target)
    hits = []
    for basis_name, basis in ALL_BASES.items():
        for form_name, form_val in basis.items():
            for f in factors:
                candidate = f * form_val
                if abs(candidate) < mpf("1e-30"):
                    continue
                rel_err = abs(candidate - target) / max(target_abs, abs(candidate))
                if rel_err < tolerance_relative:
                    hits.append({
                        "basis": basis_name,
                        "form": form_name,
                        "factor": mp.nstr(f, 10),
                        "candidate_value": mp.nstr(candidate, 15),
                        "relative_error": mp.nstr(rel_err, 8),
                        "log10_rel_err": float(mp.log10(rel_err)) if rel_err > 0 else -100,
                    })
    return sorted(hits, key=lambda h: float(h["relative_error"]))


# Run PSLQ-style searches at three tolerance levels
print("\n" + "=" * 60)
print("Part 2: Focused PSLQ on Koide-specific quantities")
print(f"Total basis: {total_basis_size} forms x {len(factors)} factors = {total_tests} candidates per target")
print("=" * 60)

part2 = {
    "basis_size_total_forms": total_basis_size,
    "factor_count": len(factors),
    "tests_per_target": total_tests,
    "n_targets": len(target_panel),
    "selection_bias_note": (
        "HAND-CURATED basis motivated by Koide-cone structural candidates. "
        "Each entry's source documented. Selection bias non-zero but bounded "
        f"by basis size {total_basis_size}. Hits within target sigma reported "
        "honestly with rel_err and context."
    ),
    "search_results": {},
}

# Tolerance levels (1%, 0.1%, 0.01%)
for tol_label, tol in [("1pct", 1e-2), ("0p1pct", 1e-3), ("0p01pct", 1e-4)]:
    part2["search_results"][tol_label] = {}
    for tname, tval in target_panel.items():
        hits = search_target(tval, tname, tolerance_relative=tol)
        part2["search_results"][tol_label][tname] = {
            "target_value": mp.nstr(tval, 20),
            "n_hits": len(hits),
            "hits": hits[:5],  # top 5
        }

# Print summary
for tol_label in ["1pct", "0p1pct", "0p01pct"]:
    print(f"\nTolerance {tol_label}:")
    for tname, result in part2["search_results"][tol_label].items():
        if result["n_hits"] > 0:
            print(f"  {tname}: {result['n_hits']} hits")
            for h in result["hits"][:2]:
                print(f"    -> {h['basis']}/{h['form']} x {h['factor']} (rel_err={h['relative_error']})")

# ============================================================
# Part 3: Structural angles
# ============================================================

print("\n" + "=" * 60)
print("Part 3: Structural angles")
print("=" * 60)

part3 = {}

# (A) Can (m_e, m_mu, m_tau) come from one Mellin transform at k=0, 1, 2?
# This would require the masses to be M-engine outputs.
# Test: do the masses (or sqrt-masses) fit any natural k-graded pattern?

# Normalize: divide by sqrt(m_tau) to get dimensionless ratios
v = [sqrt(m_e / m_tau), sqrt(m_mu / m_tau), mpf(1)]
v_str = [mp.nstr(x, 15) for x in v]

# Test: does v[i] match M_k engine output at k=0,1,2 for any natural D, t?
# This is purely structural — the engine produces continuous functions,
# but a triple of values would require the engine to be evaluated at 3 specific
# parameter values that "select" the lepton triple. No natural mechanism is known.
part3["mellin_k_mapping_test"] = {
    "v_dimensionless": v_str,
    "verdict": (
        "NO NATURAL MAPPING. Master Mellin engine produces a continuous family "
        "M[Tr(D^k e^{-tD^2})](s) of values. Three discrete values (m_e, m_mu, m_tau) "
        "would require 3 specific (k, t, s) selections with no internal mechanism. "
        "Inner-factor Mellin engine theorem (CLAUDE.md 2026-05-07) says Yukawas "
        "live in inner factor Dirichlet ring, categorically disjoint from outer M1/M2/M3. "
        "Koide cone is in inner factor; not addressable by master Mellin engine."
    ),
}

# (B) SU(3) 45° check
# Does any embedding of (1,1) adjoint, or (3,0)/(0,3) decuplet, give 45° naturally?
# The 45° (= cos^2 = 1/2) angle is the constraint that the sqrt-mass vector
# lies on a circular cone of half-opening pi/4 around the democratic direction.
#
# In SU(3) rep theory, weights of (1,1) adjoint lie on a hexagonal pattern in
# the weight lattice; democratic direction (1,1,1) under SU(3) is the singlet,
# 45° is not a natural angle in this lattice (natural angles are 60°, 120°).
#
# Sanity check: compute the angle between SU(3) (1,1) weight vector and (1,1,1).
# The fundamental weights of SU(3) form 60° angles, not 45°. Verify:

# Standard SU(3) simple roots in 3-space (with sum-to-zero constraint):
alpha1 = [mpf(1), -mpf(1), mpf(0)]
alpha2 = [mpf(0), mpf(1), -mpf(1)]
# Angle between alpha1 and alpha2:
dot12 = sum(alpha1[i] * alpha2[i] for i in range(3))
n1 = sqrt(sum(a * a for a in alpha1))
n2 = sqrt(sum(a * a for a in alpha2))
cos_alpha12 = dot12 / (n1 * n2)
angle_alpha12_deg = acos(cos_alpha12) * 180 / pi

part3["su3_45deg_check"] = {
    "su3_simple_root_angle_deg": mp.nstr(angle_alpha12_deg, 10),
    "expected_60_or_120": "120 deg (Cartan matrix of A_2)",
    "verdict": (
        "NO NATURAL 45° IN SU(3). Simple roots of SU(3) (= A_2 Lie algebra) "
        "have 120° between them. Fundamental weights have 60° to each other. "
        "45° is not a natural angle in the SU(3) weight lattice. The Koide cone "
        "is not an SU(3) representation-theoretic fact."
    ),
}

# (C) Inner-factor Mellin engine theorem
# CLAUDE.md §2 (2026-05-07): η-trivialization + AC factorization
# D^2 = D_GV^2 ⊗ 1 + 1 ⊗ D_F^2 (cross term vanishes by outer chirality anticomm)
# Combined Mellin: (outer M_i ring) × (inner Yukawa Dirichlet ring)
# Inner Dirichlet ring is ℚ[y_1^{-2s}, ..., y_n^{-2s}]
# Koide depends on (sqrt(y1), sqrt(y2), sqrt(y3)), which are NOT in the natural
# variables {y_i^{-2s}} of the inner Dirichlet ring at any integer s.
#
# Question: does any analytic property of the inner sum at specific s-values
# reduce to the Koide constraint?

# Test: Σ y_i^{-2s} at s = -1/2 is Σ y_i (sum of masses, linear in masses)
# at s = -1/4 is Σ y_i^{1/2} (sum of sqrt-masses, what Koide uses)
# at s = -1/2: y_i^1 = y_i, gives Σ y_i = m_e + m_mu + m_tau
# at s = -1/4: y_i^{1/2} = sqrt(y_i), gives Σ sqrt(y_i)
# Koide ratio is then [Σ y_i^{-2*(-1/2)}] / [Σ y_i^{-2*(-1/4)}]^2 = (Σ y_i) / (Σ √y_i)^2

# So Koide IS expressible as a ratio of inner Dirichlet values at s=-1/2 and s=-1/4.
# This is a STRUCTURAL OBSERVATION — but it doesn't derive Koide; it just
# expresses it. The constraint K = 2/3 requires the specific Yukawa values to
# satisfy a non-trivial polynomial relation.

inner_s1 = -mpf(1) / 2  # gives Σ y_i (sum of masses)
inner_s2 = -mpf(1) / 4  # gives Σ √y_i (sum of sqrt-masses)
inner_sum_s1 = m_e + m_mu + m_tau
inner_sum_s2 = sqrt(m_e) + sqrt(m_mu) + sqrt(m_tau)
koide_via_inner = inner_sum_s1 / (inner_sum_s2 ** 2)

part3["inner_factor_mellin_check"] = {
    "inner_sum_s_neg_half_equals_sum_of_masses": mp.nstr(inner_sum_s1, 20),
    "inner_sum_s_neg_quarter_equals_sum_sqrt_masses": mp.nstr(inner_sum_s2, 20),
    "koide_K_via_inner_ratio": mp.nstr(koide_via_inner, 30),
    "matches_direct_K": mp.nstr(koide_via_inner - K, 10),
    "verdict": (
        "Koide IS expressible as ratio of inner Dirichlet values at half-integer s "
        "(s=-1/2 and s=-1/4). This is a TAUTOLOGICAL reformulation, not a derivation. "
        "The K = 2/3 constraint becomes a polynomial relation among the three Yukawas "
        "y_i that does NOT follow from any axiom of the inner Dirichlet ring. "
        "Inner-factor Mellin engine theorem says the inner ring is CONSTRAINED but "
        "not DETERMINED by the framework — Koide's K=2/3 is a 1-parameter constraint "
        "the framework cannot autonomously generate."
    ),
    "half_integer_s_caveat": (
        "Half-integer s sits outside the integer-s regime where the inner Dirichlet "
        "ring is rigorously characterized. Continuation to half-integer s opens "
        "analytic-continuation freedom that doesn't help: the values are just the "
        "sums Σ y_i and Σ √y_i, which are not constrained by any framework axiom."
    ),
}

# ============================================================
# Part 4: Honest verdict
# ============================================================

# Tally PSLQ hits at most stringent tolerance
critical_tolerance = "0p01pct"
non_trivial_hits_at_critical = 0
trivial_hits_at_critical = 0  # restatements like 2/3 * 1 = 2/3
for tname, result in part2["search_results"][critical_tolerance].items():
    for h in result["hits"]:
        # Is this a structural new identification or a trivial restatement?
        if tname == "K_minus_2over3":
            # This target is by construction near zero; hits are tautological
            trivial_hits_at_critical += 1
        elif "1/" in h["form"] or "2/3" in h["form"] or h["form"] == "K_target":
            # Need finer judgment — log for review
            non_trivial_hits_at_critical += 1
        else:
            non_trivial_hits_at_critical += 1

# Total hits at 1% (the soft tolerance used in W3)
total_hits_1pct = sum(r["n_hits"] for r in part2["search_results"]["1pct"].values())
total_hits_0p01pct = sum(r["n_hits"] for r in part2["search_results"]["0p01pct"].values())

# Random-coincidence expectation: tolerance × N_targets × N_tests
# For uniform random log-scale, expected hits at tolerance ε is ε × N_targets × N_basis × N_factors
n_targets = len(target_panel)
expected_random_1pct = 0.01 * 2 * n_targets * total_basis_size * len(factors)  # factor of 2 for two-sided
expected_random_0p01pct = 1e-4 * 2 * n_targets * total_basis_size * len(factors)

part4_verdict = {
    "PSLQ_summary": {
        "total_hits_1pct": total_hits_1pct,
        "total_hits_0p01pct": total_hits_0p01pct,
        "expected_random_1pct": float(expected_random_1pct),
        "expected_random_0p01pct": float(expected_random_0p01pct),
    },
    "Koide_precision": {
        "K_vs_2over3_sigma": mp.nstr(K_dev_sigma, 4),
        "theta_vs_45deg_arcsec": mp.nstr(deviation_arcsec, 6),
    },
    "Structural_findings": {
        "mellin_k_mapping": "NO MAPPING",
        "su3_45deg": "NO — SU(3) gives 60°/120°, not 45°",
        "inner_factor_check": "Koide expressible as ratio of inner Dirichlet at half-integer s but constraint NOT derivable",
    },
}

print("\n" + "=" * 60)
print("Part 4: Summary")
print("=" * 60)
for k, v in part4_verdict.items():
    print(f"\n  {k}:")
    if isinstance(v, dict):
        for kk, vv in v.items():
            print(f"    {kk}: {vv}")

# ============================================================
# Write everything to JSON
# ============================================================

final_output = {
    "sprint": "koide_pass_investigation",
    "date": "2026-05-18",
    "context": "Focused Koide-specific re-test after Sprint W3 (May 2026) tested general lepton mass spectrum",
    "framework_state": "Papers 38-45 (Lorentzian propinquity, modular Hamiltonian, inner-factor Mellin engine theorem)",
    "precision_dps": 100,
    "part1_numerical": part1,
    "part2_pslq": part2,
    "part3_structural": part3,
    "part4_verdict": part4_verdict,
}

with open(OUT_DIR / "koide_pass_investigation.json", "w") as f:
    json.dump(final_output, f, indent=2)

print(f"\n[OK] Written to debug/data/koide_pass_investigation.json")
