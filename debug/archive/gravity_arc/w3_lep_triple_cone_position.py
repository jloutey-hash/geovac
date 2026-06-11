"""
Sprint W3-LepTriple, follow-up: where on the Koide 45-deg cone do the leptons sit?

Result of the first probe: the sqrt-mass vector z = (sqrt m_e, sqrt m_mu, sqrt m_tau)
lies essentially exactly at 45 deg from the democratic direction d = (1,1,1)/sqrt(3),
to within 1 arcsecond. This means z lives on a 2D cone of half-opening 45 deg
around d -- a one-parameter constraint that collapses three mass values to two.

But the cone is 2D. The remaining two parameters are:
  (a) the overall radius |z| (sets mass scale -- presumably tied to Higgs VEV)
  (b) the azimuthal phase phi on the cone

This script computes phi for each charged-fermion triple (leptons, up-quarks,
down-quarks) and asks whether phi shows further structure.

If phi sits at a small simple angle (0, pi/n, etc.) or matches a GeoVac internal
number, that's a second compression -- evidence that the masses sit at a
structurally-specific point, not just on a structurally-specific cone.
"""

import json
import math
from pathlib import Path

import mpmath as mp
from mpmath import mpf, pi, log, sqrt

mp.mp.dps = 80

OUT_DIR = Path("debug/data")

# ============================================================
# Charged fermion masses (PDG / CODATA 2022, MeV/c^2)
# ============================================================

# Charged leptons -- highly precise
m_e   = mpf("0.51099895069")
m_mu  = mpf("105.6583755")
m_tau = mpf("1776.86")

# Up-type quarks -- MS-bar at appropriate scale, PDG 2024 central values
# u: 2.16 MeV, c: 1.27 GeV, t: 172.69 GeV (pole, PDG)
m_u = mpf("2.16")
m_c = mpf("1270.0")
m_t = mpf("172690.0")

# Down-type quarks
# d: 4.70 MeV, s: 93.5 MeV, b: 4180 MeV
m_d = mpf("4.70")
m_s = mpf("93.5")
m_b = mpf("4180.0")

triples = {
    "leptons (e, mu, tau)":      [m_e, m_mu, m_tau],
    "up-quarks (u, c, t)":       [m_u, m_c, m_t],
    "down-quarks (d, s, b)":     [m_d, m_s, m_b],
}

# ============================================================
# Cone-parameterization helpers
# ============================================================

# Decompose z = (sqrt m_1, sqrt m_2, sqrt m_3) into (parallel, perp1, perp2)
# where parallel = projection onto (1,1,1)/sqrt(3),
# and perp1, perp2 span the orthogonal plane.
#
# Natural orthonormal basis for the perpendicular plane:
#   e1 = (2, -1, -1) / sqrt(6)
#   e2 = (0,  1, -1) / sqrt(2)
# (e1 picks out "first mass vs the other two";
#  e2 picks out "second vs third".)
#
# Phi = atan2(z . e2, z . e1) is the azimuthal phase on the cone.

def koide_K(M):
    s = sum(M)
    rs = sum(sqrt(m) for m in M)
    return s / (rs * rs)

def cone_decomposition(M):
    z = [sqrt(m) for m in M]
    z_norm_sq = sum(zi*zi for zi in z)
    z_dot_d = sum(z) / sqrt(3)
    par_sq = z_dot_d * z_dot_d
    perp_sq = z_norm_sq - par_sq
    cos_sq = par_sq / z_norm_sq
    angle_rad = mp.acos(sqrt(cos_sq))
    angle_deg = angle_rad * 180 / pi

    # Project z onto (e1, e2) in perp plane
    e1 = [mpf(2), mpf(-1), mpf(-1)]
    e1_norm = sqrt(sum(e*e for e in e1))
    e1 = [e/e1_norm for e in e1]
    e2 = [mpf(0), mpf(1), mpf(-1)]
    e2_norm = sqrt(sum(e*e for e in e2))
    e2 = [e/e2_norm for e in e2]
    proj_e1 = sum(z[i] * e1[i] for i in range(3))
    proj_e2 = sum(z[i] * e2[i] for i in range(3))
    phi_rad = mp.atan2(proj_e2, proj_e1)
    phi_deg = phi_rad * 180 / pi

    return {
        "K": koide_K(M),
        "cos_sq": cos_sq,
        "angle_from_democratic_deg": angle_deg,
        "deviation_from_45_deg": angle_deg - 45,
        "parallel_norm_sq": par_sq,
        "perp_norm_sq": perp_sq,
        "z_norm_sq": z_norm_sq,
        "phi_rad": phi_rad,
        "phi_deg": phi_deg,
        "proj_e1": proj_e1,
        "proj_e2": proj_e2,
    }

def pretty(x, n=12):
    return mp.nstr(x, n)

# ============================================================
# Run on all three families
# ============================================================

print("=" * 70)
print("Cone parameterization for charged fermion triples")
print("=" * 70)
print()
print("z = (sqrt m_1, sqrt m_2, sqrt m_3) decomposed as:")
print("  parallel: along (1,1,1)/sqrt(3) [democratic direction]")
print("  perp1:    along (2,-1,-1)/sqrt(6)  [first mass vs others]")
print("  perp2:    along (0, 1,-1)/sqrt(2)  [second vs third]")
print()
print("phi = atan2(<z, perp2>, <z, perp1>)  is the azimuthal phase on the cone.")
print("If Koide held exactly, the half-opening angle would be exactly 45 deg")
print("for all three families.")
print()

results = {}

for name, M in triples.items():
    print("-" * 70)
    print(f"  Family: {name}")
    print("-" * 70)
    decomp = cone_decomposition(M)
    print(f"  Koide K          = {pretty(decomp['K'], 12)}")
    print(f"  Angle from (1,1,1) = {pretty(decomp['angle_from_democratic_deg'], 10)} deg")
    print(f"  Deviation from 45 deg = {pretty(decomp['deviation_from_45_deg'], 6)} deg")
    print()
    print(f"  ||z||^2          = sum m_i = {pretty(decomp['z_norm_sq'], 10)} MeV")
    print(f"  ||z_par||^2      = {pretty(decomp['parallel_norm_sq'], 8)} MeV")
    print(f"  ||z_perp||^2     = {pretty(decomp['perp_norm_sq'], 8)} MeV")
    print()
    print(f"  proj on perp1  = {pretty(decomp['proj_e1'], 8)} (sqrt MeV)")
    print(f"  proj on perp2  = {pretty(decomp['proj_e2'], 8)} (sqrt MeV)")
    print()
    print(f"  azimuthal phi  = {pretty(decomp['phi_rad'], 12)} rad")
    print(f"               = {pretty(decomp['phi_deg'], 12)} deg")
    print()

    # Compare phi to natural angles
    natural_angles_deg = {
        "0":           mpf(0),
        "30":          mpf(30),
        "45":          mpf(45),
        "60":          mpf(60),
        "90":          mpf(90),
        "180/pi":      mpf(180)/pi,
        "180/pi^2":    mpf(180)/(pi*pi),
        "60 - delta":  mpf(60),  # placeholder
    }
    print(f"  Distance to natural angles:")
    for label, val in natural_angles_deg.items():
        diff = decomp['phi_deg'] - val
        print(f"    phi - {label} = {pretty(diff, 6)} deg")

    # Compare phi/(2pi) -- maybe phi is a rational multiple of 2pi
    phi_over_2pi = decomp['phi_rad'] / (2*pi)
    print(f"  phi / (2 pi)   = {pretty(phi_over_2pi, 12)}")
    # Try to find phi as fraction of pi
    phi_over_pi = decomp['phi_rad'] / pi
    print(f"  phi / pi       = {pretty(phi_over_pi, 12)}")

    # Test: is phi/pi close to a small rational?
    print(f"  Closest small rational p/q for phi/pi:")
    best = None
    for q in range(1, 20):
        for p in range(-2*q, 2*q+1):
            if math.gcd(abs(p), q) > 1:
                continue
            test = mpf(p)/q
            diff = phi_over_pi - test
            if best is None or abs(diff) < abs(best[2]):
                best = (p, q, diff)
    print(f"    {best[0]}/{best[1]} (diff = {pretty(best[2], 6)})")

    results[name] = {
        "Koide_K":               str(decomp['K']),
        "angle_from_democratic": str(decomp['angle_from_democratic_deg']),
        "dev_from_45":           str(decomp['deviation_from_45_deg']),
        "phi_rad":               str(decomp['phi_rad']),
        "phi_deg":               str(decomp['phi_deg']),
        "phi_over_pi":           str(phi_over_pi),
        "phi_over_2pi":          str(phi_over_2pi),
        "best_rational_phi_over_pi": {
            "p": int(best[0]), "q": int(best[1]),
            "value": float(mpf(best[0])/best[1]),
            "diff": float(best[2]),
        },
        "z_norm_sq_MeV":         str(decomp['z_norm_sq']),
        "proj_e1":               str(decomp['proj_e1']),
        "proj_e2":               str(decomp['proj_e2']),
    }
    print()

# ============================================================
# Test: do up-quarks AND down-quarks AND leptons all live on
# the same cone (45 deg)?  If so, this is a stronger structural fact:
# the 45-deg cone is universal across charged-fermion families.
# ============================================================

print("=" * 70)
print("UNIVERSALITY CHECK: do all three families satisfy Koide's 45 deg?")
print("=" * 70)
print()
print(f"{'Family':<30s} {'Koide K':>15s} {'angle (deg)':>15s} {'(angle - 45)':>15s}")
for name, M in triples.items():
    d = cone_decomposition(M)
    k_str = pretty(d['K'], 8)
    ang_str = pretty(d['angle_from_democratic_deg'], 8)
    dev_str = pretty(d['deviation_from_45_deg'], 4)
    print(f"  {name:<28s} {k_str:>15s} {ang_str:>15s} {dev_str:>15s}")

print()

# Save
with open(OUT_DIR / "w3_lep_triple_cone_position.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print(f"Saved to {OUT_DIR / 'w3_lep_triple_cone_position.json'}")
