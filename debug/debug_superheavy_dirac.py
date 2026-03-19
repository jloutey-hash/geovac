#!/usr/bin/env python3
"""
debug_superheavy_dirac.py
=========================
Investigation: SO(3N) structure at large N and the Dirac Z ~ 137 barrier.

Questions:
  1. Does anything special happen at N = 137 in the representation theory?
  2. Is there a maximum N beyond which mu_free formula breaks?
  3. Can we see the Dirac instability as a geometric/topological feature?
  4. What is the GeoVac prediction for the "end of the periodic table"?

Date: 2026-03-18
Status: Exploratory analysis
"""

import numpy as np
from fractions import Fraction
from math import comb, factorial, sqrt, log, pi

# =========================================================================
# PHYSICAL CONSTANTS
# =========================================================================

ALPHA = 1.0 / 137.035999084   # fine structure constant
ALPHA_INV = 137.035999084
MC2_HA = 1.0 / ALPHA**2        # mc^2 in Hartree = 1/alpha^2 ~ 18778.86

# Reference data for superheavy elements
SUPERHEAVY_DATA = {
    # name: (Z, N_neutral, config_note, predicted_chemistry)
    'Ne':  (10,  10,  '[He]2s2 2p6',        'Noble gas'),
    'Ar':  (18,  18,  '[Ne]3s2 3p6',        'Noble gas'),
    'Kr':  (36,  36,  '[Ar]3d10 4s2 4p6',   'Noble gas'),
    'Xe':  (54,  54,  '[Kr]4d10 5s2 5p6',   'Noble gas'),
    'Rn':  (86,  86,  '[Xe]4f14 5d10 6s2 6p6', 'Noble gas'),
    'Og':  (118, 118, '[Rn]5f14 6d10 7s2 7p6', 'Noble gas? Maybe not'),
    'Cn':  (112, 112, '[Rn]5f14 6d10 7s2',  '"Relativistic noble liquid"'),
    'Fl':  (114, 114, '[Rn]5f14 6d10 7s2 7p2', 'More like Hg than Pb'),
    'Nh':  (113, 113, '[Rn]5f14 6d10 7s2 7p1', 'Thallium analog'),
    'Mc':  (115, 115, '[Rn]5f14 6d10 7s2 7p3', 'Bismuth analog'),
    'Lv':  (116, 116, '[Rn]5f14 6d10 7s2 7p4', 'Polonium analog'),
    'Ts':  (117, 117, '[Rn]5f14 6d10 7s2 7p5', 'Astatine analog'),
}

# Exact total energies (Hartree) for lighter atoms (Chakravorty et al.)
EXACT_ENERGIES = {
    1: -0.5,           # H
    2: -2.9037,        # He
    3: -7.4781,        # Li
    4: -14.6674,       # Be
    10: -128.9376,     # Ne
    18: -527.5400,     # Ar  (approx)
    36: -2752.055,     # Kr  (approx)
    54: -7232.14,      # Xe  (approx)
}

# Thomas-Fermi total energy: E_TF = -0.7687 * Z^(7/3) Hartree
TF_COEFF = -0.7687


# =========================================================================
# STEP 1: Extrapolate mu_free formula across the periodic table
# =========================================================================

print("=" * 90)
print("STEP 1: mu_free = 2(N-2)^2 EXTRAPOLATION")
print("=" * 90)

print(f"""
  Ground state formula (all atoms with S < N/2):
    nu = N - 2
    mu_free = nu * (nu + 3N - 2) / 2 = (N-2)(4N-4)/2 = 2(N-2)^2

  This is the SO(3N) Casimir eigenvalue for the ground-state spatial irrep.
""")

# Generate table
N_values = [1, 2, 3, 4, 10, 18, 36, 54, 86, 100, 112, 114, 118,
            120, 126, 137, 150, 170, 184, 200]

print(f"  {'N':>5}  {'Element':>8}  {'S^(3N-1)':>10}  {'dim(SO(3N))':>12}  "
      f"{'nu':>6}  {'mu_free':>10}  {'l_eff':>8}  {'n_eff':>8}")
print("-" * 90)

element_names = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 10: 'Ne', 18: 'Ar',
    36: 'Kr', 54: 'Xe', 86: 'Rn', 100: 'Fm', 112: 'Cn',
    114: 'Fl', 118: 'Og', 120: 'Ubn', 126: 'Ubh',
    137: 'Z=137', 150: 'Z=150', 170: 'Z=170', 184: 'Z=184', 200: 'Z=200'
}

results = {}
for N in N_values:
    d = 3 * N              # dimension of config space R^(3N)
    sphere_dim = d - 1     # S^(3N-1) dimension
    so_dim = d * (d - 1) // 2  # dim of SO(3N) Lie algebra
    nu = max(N - 2, 0)
    mu_free = 2 * (N - 2)**2 if N >= 2 else 0
    l_eff = Fraction(3*N - 3, 2)
    n_eff = Fraction(3*N - 1, 2)
    elem = element_names.get(N, f'Z={N}')

    results[N] = {
        'elem': elem, 'sphere_dim': sphere_dim, 'so_dim': so_dim,
        'nu': nu, 'mu_free': mu_free, 'l_eff': float(l_eff),
        'n_eff': float(n_eff)
    }

    print(f"  {N:>5}  {elem:>8}  S^{sphere_dim:<6d}    {so_dim:>10,d}  "
          f"{nu:>6,d}  {mu_free:>10,d}  {float(l_eff):>8.1f}  {float(n_eff):>8.1f}")

print(f"""
  KEY OBSERVATIONS:
  - mu_free grows as ~2N^2 (quadratic). No singularity, no divergence.
  - At N = 137: mu_free = {2*(137-2)**2:,d}, living on S^{3*137-1}
  - At N = 170: mu_free = {2*(170-2)**2:,d}, living on S^{3*170-1}
  - Nothing special happens at N = 137 in the mu_free formula.
  - The formula is SMOOTH through Z = 1/alpha.
""")


# =========================================================================
# STEP 2: Casimir eigenvalue density -- ground state vs maximum
# =========================================================================

print("=" * 90)
print("STEP 2: CASIMIR FILLING FRACTION -- HOW CLOSE IS THE GROUND STATE TO THE EDGE?")
print("=" * 90)

print(f"""
  On S^(d-1) with d = 3N, the Laplace-Beltrami eigenvalues are:
    lambda_nu = -nu(nu + d - 2)  for nu = 0, 1, 2, ...

  The CONTINUOUS sphere has no maximum nu (nu -> infinity).
  But for a DISCRETE graph with V vertices, nu_max ~ V^(1/(d-1)).

  For the GeoVac lattice, V = sum_n (2n-1) from n=1 to nmax.
  At nmax, V = nmax^2. So nu_max ~ nmax.

  But physically, the question is: how does nu_gs = N-2 relate to
  fundamental scales of the geometry?

  KEY RATIO: nu / (d - 2) = (N-2) / (3N-2)
  This is the "excitation fraction" -- how far up the spectrum the
  ground state sits relative to the natural scale of the sphere.
""")

print(f"  {'N':>5}  {'nu':>6}  {'d-2':>6}  {'nu/(d-2)':>10}  {'-> limit':>10}")
print("-" * 55)

for N in N_values:
    nu = max(N - 2, 0)
    d_minus_2 = 3*N - 2
    ratio = nu / d_minus_2 if d_minus_2 > 0 else 0
    limit_val = 1/3  # asymptotic limit
    print(f"  {N:>5}  {nu:>6,d}  {d_minus_2:>6,d}  {ratio:>10.6f}  {'1/3':>10}")

print(f"""
  RESULT: nu/(d-2) -> 1/3 monotonically from below.
  The ground state always sits at ~1/3 of the natural angular scale.
  This ratio is SMOOTH through N = 137. No critical behavior.

  Physical interpretation: the Pauli exclusion principle forces the
  ground state to have nu = N-2, but the sphere dimension grows as 3N,
  so the "excitation fraction" stabilizes at 1/3.

  This means the ground state is NEVER close to the edge of the
  representation space. There is NO topological squeeze.
""")


# =========================================================================
# STEP 3: The Dirac energy in hyperspherical coordinates
# =========================================================================

print("=" * 90)
print("STEP 3: DIRAC ENERGY vs QUASI-COULOMB IN HYPERSPHERICAL COORDINATES")
print("=" * 90)

print(f"""
  NON-RELATIVISTIC (Schrodinger):
    E_1s = -Z^2 / 2  (hydrogen-like, Hartree)

  RELATIVISTIC (Dirac):
    E_1s = mc^2 * [sqrt(1 - (Z*alpha)^2) - 1]    (point nucleus)
         = mc^2 * [sqrt(1 - Z^2/137.036^2) - 1]

  At Z = 137: (Z*alpha)^2 ~ 1, so E_1s -> -mc^2 = -{MC2_HA:.1f} Ha
  This is the "dive into the continuum."

  QUASI-COULOMB (hyperspherical, Paper 13):
    E ~ -Z_eff^2 / (2 * n_eff^2)
    where n_eff = (3N-1)/2

  For a ONE-electron atom (N=1): n_eff = 1, E = -Z^2/2 (exact match).
  For MULTI-electron neutral atoms: n_eff = (3N-1)/2, Z_eff < Z due to screening.

  Let's compare these three models for hydrogen-like ions (N=1, vary Z):
""")

print(f"  {'Z':>5}  {'E_NR (Ha)':>12}  {'E_Dirac (Ha)':>14}  {'Ratio E_D/E_NR':>16}  "
      f"{'v/c = Z*alpha':>14}  {'Dive?':>6}")
print("-" * 80)

Z_scan = [1, 10, 20, 30, 50, 70, 80, 90, 100, 110, 118, 130, 137, 140, 150, 170]

for Z in Z_scan:
    E_NR = -Z**2 / 2.0
    Za = Z * ALPHA
    Za2 = Za**2

    if Za2 < 1.0:
        E_Dirac = MC2_HA * (sqrt(1.0 - Za2) - 1.0)
        ratio = E_Dirac / E_NR
        dive = ""
    else:
        E_Dirac = float('-inf')
        ratio = float('inf')
        dive = " <-- DIVE"

    if Za2 < 1.0:
        print(f"  {Z:>5}  {E_NR:>12.2f}  {E_Dirac:>14.2f}  {ratio:>16.6f}  "
              f"{Za:>14.6f}{dive}")
    else:
        print(f"  {Z:>5}  {E_NR:>12.2f}  {'IMAGINARY':>14}  {'---':>16}  "
              f"{Za:>14.6f}{dive}")

print(f"""
  KEY INSIGHT:
  The Dirac energy diverges at Z*alpha = 1 (Z ~ 137), but the
  non-relativistic energy is perfectly smooth: E_NR = -Z^2/2.

  The hyperspherical formalism is NON-RELATIVISTIC. It recovers
  E = -Z^2/(2n^2) exactly via the Fock projection (Paper 7).

  Therefore: the Dirac limit Z ~ 137 is INVISIBLE to the
  non-relativistic SO(3N) representation theory.

  The mu_free formula does not break. The topology does not change.
  The Dirac instability is a RELATIVISTIC effect that lies outside
  the current GeoVac framework.
""")


# =========================================================================
# STEP 4: Relativistic corrections as geometric perturbations
# =========================================================================

print("=" * 90)
print("STEP 4: RELATIVISTIC CORRECTIONS AS GEOMETRIC PERTURBATIONS ON S^(3N-1)")
print("=" * 90)

print(f"""
  The Dirac equation for hydrogen-like atoms gives:
    E_nj = mc^2 / sqrt(1 + (Z*alpha/(n - delta_j))^2) - mc^2

  where delta_j = j + 1/2 - sqrt((j+1/2)^2 - (Z*alpha)^2)

  Expanding to O((Z*alpha)^4):
    E_nj ~ -Z^2/(2n^2) * [1 + (Z*alpha)^2/n * (1/(j+1/2) - 3/(4n)) + ...]

  The relativistic correction to the non-relativistic energy is:
    Delta_E / E_NR = (Z*alpha)^2/n * (1/(j+1/2) - 3/(4n))

  IN HYPERSPHERICAL TERMS:
  The Fock projection maps the graph Laplacian eigenvalues to
  E_n = -Z^2/(2n^2). The relativistic correction would modify
  the METRIC on S^(3N-1):

    g_ij -> g_ij + (Z*alpha)^2 * delta_g_ij

  This is a CONFORMAL PERTURBATION of the unit S^3 geometry.
  The Casimir eigenvalue acquires a correction:

    mu -> mu + delta_mu,  where delta_mu ~ (Z*alpha)^2 * mu

  Let's quantify this for the 1s electron at various Z:
""")

print(f"  {'Z':>5}  {'(Za)^2':>10}  {'Delta_E/E_NR':>14}  {'Status':>30}")
print("-" * 70)

for Z in [1, 10, 30, 50, 79, 92, 100, 118, 130, 137]:
    Za = Z * ALPHA
    Za2 = Za**2
    # 1s correction: n=1, j=1/2 -> 1/(j+1/2) - 3/4 = 1 - 3/4 = 1/4
    # Delta_E/E = (Za)^2 * (1/4) for 1s
    # Actually more precisely: (Za)^2 * (1 - 3/4) = (Za)^2/4
    # But let's use the exact Dirac vs NR comparison
    E_NR = -Z**2 / 2.0
    if Za2 < 1.0:
        E_D = MC2_HA * (sqrt(1.0 - Za2) - 1.0)
        correction = (E_D - E_NR) / abs(E_NR)
        if abs(correction) < 0.01:
            status = "Non-relativistic regime"
        elif abs(correction) < 0.1:
            status = "Perturbative corrections"
        elif abs(correction) < 0.5:
            status = "STRONGLY relativistic"
        else:
            status = "APPROACHING DIVE"
        print(f"  {Z:>5}  {Za2:>10.6f}  {correction:>14.6f}  {status:>30}")
    else:
        print(f"  {Z:>5}  {Za2:>10.6f}  {'DIVERGENT':>14}  {'BEYOND DIRAC LIMIT':>30}")

print(f"""
  GEOMETRIC INTERPRETATION:
  =========================
  The non-relativistic S^3 (or S^(3N-1)) has a FLAT conformal factor.
  Relativistic corrections WARP this sphere:

    ds^2 = (1 + (Za)^2 * f(theta)) * ds^2_flat

  where f(theta) concentrates near the "poles" (= nucleus).

  At Z*alpha << 1: small warping, perturbation theory works.
  At Z*alpha ~ 0.5 (Z ~ 70): significant warping, QED corrections needed.
  At Z*alpha -> 1 (Z ~ 137): the sphere PINCHES at the nuclear pole.

  The "pinch" is the geometric avatar of the Dirac dive:
  the conformal factor goes to zero at the nuclear contact point,
  creating a conical singularity on S^(3N-1).

  This is NOT captured by the non-relativistic mu_free formula,
  which assumes a ROUND (constant curvature) sphere.
""")


# =========================================================================
# STEP 5: N = 137 -- special structure analysis
# =========================================================================

print("=" * 90)
print("STEP 5: DETAILED ANALYSIS AT N = 137 (Z = 1/alpha)")
print("=" * 90)

N_137 = 137
d_137 = 3 * N_137
nu_137 = N_137 - 2
mu_137 = 2 * (N_137 - 2)**2
l_eff_137 = (3 * N_137 - 3) / 2.0
n_eff_137 = (3 * N_137 - 1) / 2.0

# S_N representation theory
# Ground state: S = ? For Z=137, Hund-rule config would give some S.
# But for the representation theory, let's use S < N/2 -> nu = N-2.
# The actual ground state spin depends on electron configuration.

# Z=137 is in the 8th period (if it existed):
# [Rn]5f14 6d10 7s2 7p6 = 118 electrons  (Og core)
# Remaining: 137 - 118 = 19 electrons beyond Og
# These fill: 8s2, 5g18 -> 8s2 5g17 = 19 electrons
# So Z=137 would have a partially filled 5g shell (!)
# g-shell: l=4, capacity = 2(2*4+1) = 18
# Config: [Og] 8s2 5g17 -> one hole in the g-shell
# S = 1/2 (doublet)

print(f"""
  ATOM: Z = N = 137 (hypothetical, well beyond synthesis limit)

  Predicted electron configuration: [Og] 8s^2 5g^17
    (one hole in the first g-shell: l = 4)
    Ground state spin: S = 1/2 (doublet, like a halogen analog)

  SO(3N) PARAMETERS:
    Dimension of configuration space: R^{d_137} (= R^{{{d_137}}})
    Hypersphere: S^{d_137 - 1} (= S^{{{d_137 - 1}}})
    SO({d_137}) symmetry group
    Lie algebra dimension: {d_137 * (d_137 - 1) // 2:,d}

  GROUND STATE REPRESENTATION:
    nu = N - 2 = {nu_137}
    mu_free = 2(N-2)^2 = {mu_137:,d}
    l_eff = (3N-3)/2 = {l_eff_137:.1f}
    n_eff = (3N-1)/2 = {n_eff_137:.1f}

  QUASI-COULOMB ESTIMATE:
    E ~ -Z_eff^2 / (2 * n_eff^2)
""")

# For a neutral atom, Z_eff for the outermost electron sees
# screening from N-1 electrons. Slater's rules give Z_eff ~ Z - sigma.
# For a 5g electron in Z=137: sigma ~ 100+ (heavy screening)
# But for total energy, Thomas-Fermi gives E_total ~ -0.7687 * Z^(7/3)

E_TF = TF_COEFF * 137**(7.0/3.0)
E_per_electron = E_TF / 137

print(f"  Thomas-Fermi total energy estimate: E_TF = {E_TF:,.0f} Ha")
print(f"  Energy per electron: {E_per_electron:,.1f} Ha")
print()

# Check for alpha-related ratios
print("  SEARCHING FOR alpha-RELATED RATIOS AT N = 137:")
print("  " + "-" * 60)

ratios = {
    'nu / N': nu_137 / N_137,
    'mu_free / N^2': mu_137 / N_137**2,
    'mu_free / N^3': mu_137 / N_137**3,
    'l_eff / N': l_eff_137 / N_137,
    'n_eff / N': n_eff_137 / N_137,
    'nu / (3N)': nu_137 / (3 * N_137),
    'nu / l_eff': nu_137 / l_eff_137,
    'nu * alpha': nu_137 * ALPHA,
    'N * alpha': N_137 * ALPHA,
    'mu_free * alpha^2': mu_137 * ALPHA**2,
    'sqrt(mu_free) * alpha': sqrt(mu_137) * ALPHA,
    'nu^2 / mu_free': nu_137**2 / mu_137,
    '1/3 - nu/(3N-2)': 1/3 - nu_137/(3*N_137-2),
}

for name, val in ratios.items():
    # Check if close to 1, alpha, 1/alpha, pi, or simple fraction
    matches = []
    if abs(val - 1.0) < 0.01:
        matches.append("~ 1")
    if abs(val - ALPHA) < 0.001:
        matches.append(f"~ alpha = {ALPHA:.6f}")
    if abs(val - ALPHA_INV) < 1.0:
        matches.append(f"~ 1/alpha = {ALPHA_INV:.1f}")
    if abs(val - pi) < 0.01:
        matches.append(f"~ pi")
    if abs(val - 0.5) < 0.01:
        matches.append("~ 1/2")
    if abs(val - 1.0/3.0) < 0.01:
        matches.append("~ 1/3")
    if abs(val - 2.0/3.0) < 0.01:
        matches.append("~ 2/3")
    match_str = " <-- " + ", ".join(matches) if matches else ""
    print(f"    {name:>25} = {val:>14.8f}{match_str}")

print(f"""
  N * alpha = {N_137 * ALPHA:.8f}

  At N = 137, the product N * alpha = {N_137 * ALPHA:.6f} ~ 1.
  This is just the DEFINITION of alpha^(-1) ~ 137.

  But this is a property of Z (nuclear charge), not N (topology).
  The SO(3N) structure knows about N, not Z.

  CRITICAL DISTINCTION:
  - N determines the TOPOLOGY: S^(3N-1), SO(3N), nu, mu_free
  - Z determines the POTENTIAL: -Z/r attraction strength
  - alpha connects Z to RELATIVISTIC effects: (Z*alpha)^2

  The topology is perfectly smooth at N = 137.
  The potential is perfectly smooth at Z = 137.
  Only the RELATIVISTIC CORRECTION (Z*alpha)^2 diverges at Z ~ 137.

  CONCLUSION: There is NO special topological structure at N = 137.
  The "magic" of Z = 137 is purely a relativistic effect.
""")


# =========================================================================
# STEP 6: Z vs N separation -- fixing N, varying Z
# =========================================================================

print("=" * 90)
print("STEP 6: FIXED N, VARYING Z -- SEPARATING TOPOLOGY FROM POTENTIAL")
print("=" * 90)

print(f"""
  Fix N = 118 (Oganesson electron count) and vary Z.
  The topology (S^353, SO(354)) stays the same.
  Only the potential changes.

  Quasi-Coulomb energy for the OUTERMOST electron:
    E_outer ~ -(Z - sigma)^2 / (2 * n_outer^2)

  where sigma ~ N - 1 (screening) and n_outer is the outermost shell.

  For Og-like: outermost is 7p, n_outer = 7.
  Z_eff_outer = Z - sigma, where sigma is from Slater's rules.

  Total energy (Thomas-Fermi): E ~ -0.7687 * Z^(7/3) * (1 - corrections)

  Let's look at the 1s electron energy as a function of Z/N:
""")

N_fixed = 118
print(f"  N = {N_fixed} (fixed), varying Z")
print(f"  {'Z':>5}  {'Z/N':>6}  {'E_1s_NR (Ha)':>14}  {'E_1s_Dirac (Ha)':>16}  "
      f"{'(Za)^2':>8}  {'Ion':>8}  {'Note':>20}")
print("-" * 95)

for Z in [50, 80, 100, 118, 130, 137, 140, 150, 170, 200]:
    Za = Z * ALPHA
    Za2 = Za**2
    E_1s_NR = -Z**2 / 2.0  # hydrogen-like 1s

    ion_charge = Z - N_fixed
    ion_label = f"{'+'*min(abs(ion_charge),3)}{abs(ion_charge)}" if ion_charge > 0 else "neutral" if ion_charge == 0 else f"-{abs(ion_charge)}"

    if Za2 < 1.0:
        E_1s_D = MC2_HA * (sqrt(1.0 - Za2) - 1.0)
        note = ""
        if Za2 > 0.5:
            note = "STRONGLY relativistic"
        print(f"  {Z:>5}  {Z/N_fixed:>6.2f}  {E_1s_NR:>14.1f}  {E_1s_D:>16.1f}  "
              f"{Za2:>8.4f}  {ion_label:>8}  {note:>20}")
    else:
        print(f"  {Z:>5}  {Z/N_fixed:>6.2f}  {E_1s_NR:>14.1f}  {'SUPERCRITICAL':>16}  "
              f"{Za2:>8.4f}  {ion_label:>8}  {'DIRAC DIVE':>20}")

print(f"""
  KEY INSIGHT:
  The topology (N = 118) is UNCHANGED as Z varies.
  The Dirac dive happens at Z ~ 137 regardless of N.
  It is purely a POTENTIAL effect, not a topological one.

  For Z = 170 with N = 118 (highly ionized, charge +52):
  - Topology: still S^353, nu = 116, mu_free = 26,912
  - 1s Dirac: SUPERCRITICAL (beyond point-nucleus limit)
  - In practice: finite nuclear size extends to Z ~ 170

  The GeoVac framework computes the TOPOLOGY correctly for any N.
  Relativistic corrections must be added as PERTURBATIONS to the
  metric on S^(3N-1), not as modifications to the representation theory.
""")


# =========================================================================
# STEP 7: Topology flow at extreme N
# =========================================================================

print("=" * 90)
print("STEP 7: TOPOLOGY FLOW (SURGERY) AT EXTREME N")
print("=" * 90)

print(f"""
  IONIZATION AS TOPOLOGY SURGERY:
  When an atom loses one electron: N -> N-1
  The topology changes: S^(3N-1) -> S^(3(N-1)-1) x R = S^(3N-4) x R

  This is a "dimension-3 surgery" regardless of N.
  The surgery always removes exactly 3 dimensions (one electron's R^3).

  For small atoms:
    He (N=2): S^5 -> S^2 x R       (surgery removes 3 of 5 dims)
    Li (N=3): S^8 -> S^5 x R       (surgery removes 3 of 8 dims)
    Be (N=4): S^11 -> S^8 x R      (surgery removes 3 of 11 dims)

  For large atoms:
    Og (N=118): S^353 -> S^350 x R  (surgery removes 3 of 353 dims)
    Z=137: S^410 -> S^407 x R       (surgery removes 3 of 410 dims)
    Z=170: S^509 -> S^506 x R       (surgery removes 3 of 509 dims)

  FRACTIONAL DIMENSION LOSS:
""")

print(f"  {'N':>5}  {'d_before':>10}  {'d_after':>10}  {'frac_loss':>12}  {'IE (eV)':>10}  {'Note':>25}")
print("-" * 80)

# IE data for noble gases and some elements
IE_data = {
    2: 24.587, 10: 21.565, 18: 15.760, 36: 14.000, 54: 12.130,
    86: 10.749, 118: 11.9  # predicted
}

for N in [2, 3, 4, 10, 18, 36, 54, 86, 118, 137, 170, 200]:
    d_before = 3*N - 1
    d_after = 3*(N-1) - 1
    frac_loss = 3.0 / d_before
    ie = IE_data.get(N, None)
    ie_str = f"{ie:.1f}" if ie is not None else "---"
    note = ""
    if frac_loss > 0.3:
        note = "Large fraction lost"
    elif frac_loss < 0.01:
        note = "Tiny fraction lost"
    print(f"  {N:>5}  {d_before:>10}  {d_after:>10}  {frac_loss:>12.6f}  {ie_str:>10}  {note:>25}")

print(f"""
  The fractional dimension loss 3/(3N-1) -> 0 as N -> infinity.
  Removing one electron from a 200-electron atom changes the
  sphere dimension by only 0.5%.

  IMPLICATION: For heavy atoms, ionization is a SMALL perturbation
  to the topology. The distinction between Z and Z+1 becomes
  geometrically negligible. This is why:
  - Heavy atoms are easier to ionize (less topological disruption)
  - Noble gas IE decreases down the group (Ne 21.6 > Ar 15.8 > Kr 14.0 > ...)
  - Superheavy chemistry converges to a "sea of electrons" picture
""")

# Check for fragmentation instability
print("  FRAGMENTATION ANALYSIS:")
print("  " + "-" * 60)
print(f"""
  Question: At large N, could the atom FRAGMENT rather than ionize?
  S^(3N-1) -> S^(3K-1) x S^(3(N-K)-1) x R  (split into two fragments)

  For this to be energetically favorable:
    E(N) > E(K) + E(N-K)  (fragmentation lowers energy)

  Using Thomas-Fermi: E(N) ~ -0.769 * N^(7/3) (neutral atom)
    E(K) + E(N-K) ~ -0.769 * [K^(7/3) + (N-K)^(7/3)]

  Since x^(7/3) is CONVEX: K^(7/3) + (N-K)^(7/3) < N^(7/3) for K in (0,N).
  Therefore: E(K) + E(N-K) > E(N) (less negative = higher energy).

  CONCLUSION: Fragmentation is ALWAYS energetically unfavorable
  for neutral atoms. The atom prefers to stay whole.

  But this assumes the nucleus is a single point charge Z.
  For Z > ~126, nuclear instability (fission) destroys the atom
  before electronic instability matters.
""")

# Verify convexity
print("  Verification of Thomas-Fermi convexity:")
print(f"  {'N':>5}  {'E(N)':>12}  {'E(N/2)+E(N/2)':>16}  {'Ratio':>8}  {'Whole < Parts?':>16}")
print("-" * 65)

for N in [10, 50, 100, 137, 170, 200]:
    E_whole = abs(TF_COEFF) * N**(7.0/3.0)
    E_parts = 2 * abs(TF_COEFF) * (N/2)**(7.0/3.0)
    ratio = E_parts / E_whole
    # E_whole should be larger (more negative = more bound)
    stable = "YES (stable)" if E_whole > E_parts else "NO (fragment)"
    print(f"  {N:>5}  {E_whole:>12.0f}  {E_parts:>16.0f}  {ratio:>8.4f}  {stable:>16}")

print(f"""
  Ratio = 2*(N/2)^(7/3) / N^(7/3) = 2^(1 - 7/3) = 2^(-4/3) = {2**(-4/3):.4f}

  The fragments have only {2**(-4/3)*100:.1f}% of the binding energy.
  Atoms are ALWAYS more stable whole (electronically).
  {2**(-4/3)*100:.1f}% < 100% means fragments are less bound -> atom stays whole.
""")


# =========================================================================
# STEP 8: The alpha connection -- Paper 2 perspective
# =========================================================================

print("=" * 90)
print("STEP 8: THE ALPHA CONNECTION -- GEOMETRIC IMPEDANCE REVISITED")
print("=" * 90)

print(f"""
  In Paper 2, alpha emerged from the n=5 shell as a geometric impedance:
    alpha^(-1) = sum of geometric factors at n=5

  This was a CONJECTURE (papers/conjectures/Paper_2).

  The question: does SO(3N) structure at N = 137 connect to alpha?

  ANALYSIS:
  =========
  The SO(3N) Casimir at N = 1/alpha:
    nu = 1/alpha - 2 = 135.036
    mu_free = 2 * (1/alpha - 2)^2 = 2 * 135.036^2 = 36,469.4

  But 1/alpha is NOT an integer, so "N = 1/alpha" is not well-defined
  for the representation theory. N must be an integer (electron count).

  At N = 137: nu = 135, mu_free = 36,450
  At N = 138: nu = 136, mu_free = 36,992

  The "alpha barrier" falls BETWEEN two integers. There is no
  distinguished integer at Z = 1/alpha.

  DEEPER QUESTION: Is there a formula connecting mu_free to alpha?
""")

# Search for integer N where mu_free has special alpha-related properties
print("  Scanning for alpha-related coincidences...")
print(f"  {'N':>5}  {'mu_free':>10}  {'mu_free * alpha^2':>18}  {'sqrt(mu_free)*alpha':>22}  {'mu_free mod 137':>16}")
print("-" * 80)

for N in range(130, 145):
    mu = 2 * (N - 2)**2
    val1 = mu * ALPHA**2
    val2 = sqrt(mu) * ALPHA
    val3 = mu % 137
    markers = []
    if abs(val1 - round(val1)) < 0.05:
        markers.append(f"mu*a^2 ~ {round(val1)}")
    if abs(val2 - 1.0) < 0.05:
        markers.append("sqrt(mu)*a ~ 1")
    mark = " <-- " + ", ".join(markers) if markers else ""
    print(f"  {N:>5}  {mu:>10,d}  {val1:>18.6f}  {val2:>22.8f}  {val3:>16d}{mark}")

print(f"""
  mu_free * alpha^2 at N=137: {2*(137-2)**2 * ALPHA**2:.4f}
  This is approximately 1.94 -- close to 2 but NOT exact.

  sqrt(mu_free) * alpha at N=137: {sqrt(2*(137-2)**2) * ALPHA:.6f}
  = sqrt(2) * 135 / 137.036 = sqrt(2) * 0.98514 = {sqrt(2) * 135/137.036:.6f}
  Close to sqrt(2) * (1 - 2*alpha) but no clean identity.

  VERDICT: No compelling alpha-related structure at N = 137.
  The number 137 is special for the DIRAC EQUATION (relativistic QM),
  not for the SO(3N) REPRESENTATION THEORY (non-relativistic topology).

  This is consistent with the two-level picture:
  - Level 1 (topology): knows about N, oblivious to alpha
  - Level 2 (perturbation): knows about Z, alpha enters here
  - The Dirac limit is a Level 2 (or higher) effect
""")


# =========================================================================
# STEP 9: Superheavy element chemistry predictions
# =========================================================================

print("=" * 90)
print("STEP 9: SUPERHEAVY ELEMENT PREDICTIONS -- SO(3N) vs RELATIVISTIC")
print("=" * 90)

print(f"""
  Known relativistic predictions for superheavy elements:
  =====================================================

  Cn (Z=112): 6d10 7s2 closed subshells.
    Relativistic prediction: "noble liquid" (7s electrons deeply bound
    by relativistic contraction, chemically inert like Hg but more so).
    SO(3N) prediction: Type D (closed s^2 over d^10 core).
    COMPARISON: Both predict low reactivity, but for different reasons.
    Rel: 7s contraction. Topo: Type D structure.

  Fl (Z=114): 6d10 7s2 7p2.
    Relativistic prediction: 7p1/2 spinor deeply bound, acts like
    a closed-shell "7p1/2^2" rather than open p-shell.
    SO(3N) prediction: Type E (open p-shell), should be reactive.
    COMPARISON: CONFLICT. Relativistic effects completely change the
    chemistry. The non-relativistic topology predicts wrong behavior.

  Og (Z=118): 6d10 7s2 7p6.
    Relativistic prediction: 7p3/2 spinor may be loosely bound,
    Og might NOT be a noble gas. Predicted to be a semiconductor!
    SO(3N) prediction: Type B (noble gas, closed shell).
    COMPARISON: CONFLICT. Relativistic spin-orbit splitting changes
    the effective shell structure.
""")

# Structure type predictions for Z = 100 to 120
print("  SO(3N) STRUCTURE TYPE PREDICTIONS (NON-RELATIVISTIC):")
print(f"  {'Z':>5}  {'Name':>4}  {'Config':>25}  {'Type':>6}  {'nu':>5}  {'mu_free':>10}  {'Note':>25}")
print("-" * 90)

sh_atoms = [
    (100, 'Fm', '[Rn]5f12 7s2', 'E'),
    (102, 'No', '[Rn]5f14 7s2', 'D'),
    (104, 'Rf', '[Rn]5f14 6d2 7s2', 'E'),
    (110, 'Ds', '[Rn]5f14 6d8 7s2', 'E'),
    (112, 'Cn', '[Rn]5f14 6d10 7s2', 'D'),
    (114, 'Fl', '[Rn]5f14 6d10 7s2 7p2', 'E'),
    (116, 'Lv', '[Rn]5f14 6d10 7s2 7p4', 'E'),
    (118, 'Og', '[Rn]5f14 6d10 7s2 7p6', 'B'),
    (119, '119', '[Og] 8s1', 'C'),
    (120, '120', '[Og] 8s2', 'D'),
]

for Z, name, config, stype in sh_atoms:
    N = Z  # neutral
    nu = N - 2
    mu = 2 * (N - 2)**2
    note = ""
    if Z == 112:
        note = "Rel: noble liquid"
    elif Z == 114:
        note = "Rel: more like Hg"
    elif Z == 118:
        note = "Rel: maybe not noble!"
    elif Z == 119:
        note = "Alkali analog"
    elif Z == 120:
        note = "Alkaline earth analog"
    print(f"  {Z:>5}  {name:>4}  {config:>25}  {stype:>6}  {nu:>5}  {mu:>10,d}  {note:>25}")

print(f"""
  WHERE SO(3N) AND RELATIVITY AGREE:
  - Structure types C and D (alkali and alkaline earth analogs for Z=119,120)
    should be correct even with relativistic effects.
  - The qualitative hierarchy (core + valence) survives relativity.

  WHERE THEY DISAGREE:
  - Open p-shell atoms (Fl, Lv): relativistic spin-orbit splitting
    changes effective shell structure.
  - Og: the non-relativistic noble gas prediction may be wrong.
  - Any atom with (Z*alpha) > 0.5: relativistic modifications to
    the metric on S^(3N-1) become non-perturbative.
""")


# =========================================================================
# STEP 10: The geometric "end of the periodic table"
# =========================================================================

print("=" * 90)
print("STEP 10: THE GEOMETRIC END OF THE PERIODIC TABLE")
print("=" * 90)

print(f"""
  THREE LIMITS on the maximum Z:

  LIMIT 1: NUCLEAR STABILITY (Z ~ 126)
  =====================================
  Proton-proton Coulomb repulsion overwhelms the strong force.
  Predicted "island of stability" around Z = 114-126 with
  magic numbers N_p = 114 or 126 and N_n = 184.
  Beyond Z ~ 130, half-lives drop below microseconds.
  This is a NUCLEAR physics limit, not an electronic one.

  LIMIT 2: DIRAC INSTABILITY (Z ~ 137 point, Z ~ 170 finite nucleus)
  ====================================================================
  At Z*alpha = 1 (Z ~ 137), 1s binding energy = mc^2.
  With finite nuclear radius: critical Z shifts to ~170.
  Beyond Z_cr ~ 170: spontaneous e+e- pair creation.
  The vacuum itself becomes unstable.
  This is a RELATIVISTIC QED limit.

  LIMIT 3: TOPOLOGICAL/GEOMETRIC (Z -> infinity?)
  ================================================
  The SO(3N) representation theory has NO finite limit.
  mu_free = 2(N-2)^2 is defined for all integer N >= 1.
  The excitation fraction nu/(3N-2) -> 1/3 (never reaches 1).
  The topology S^(3N-1) is well-defined for all N.

  Therefore: The non-relativistic GeoVac framework predicts NO
  upper limit on the periodic table. Every integer N corresponds
  to a valid ground state on S^(3N-1).

  The ACTUAL limit comes from physics that GeoVac does not (yet) include:
  1. Nuclear instability (strong force / Coulomb competition)
  2. QED vacuum instability (Dirac dive, pair creation)
  3. Finite speed of light (relativistic metric correction)
""")

# Quantify the three limits
print("  COMPARISON OF THREE LIMITS:")
print(f"  {'Z':>5}  {'Nuclear':>12}  {'Dirac (pt)':>12}  {'Dirac (fin)':>14}  {'SO(3N)':>12}")
print("-" * 65)

for Z in [80, 100, 118, 126, 130, 137, 150, 170, 200]:
    Za = Z * ALPHA
    Za2 = Za**2

    # Nuclear: approximate stability
    if Z <= 118:
        nuc = "Observed"
    elif Z <= 126:
        nuc = "Possible"
    elif Z <= 130:
        nuc = "Marginal"
    else:
        nuc = "Unstable"

    # Dirac point nucleus
    if Za2 < 0.5:
        dirac_pt = "OK"
    elif Za2 < 0.9:
        dirac_pt = "Corrected"
    elif Za2 < 1.0:
        dirac_pt = "CRITICAL"
    else:
        dirac_pt = "BEYOND"

    # Dirac finite nucleus (limit ~ 170)
    if Z < 150:
        dirac_fin = "OK"
    elif Z < 170:
        dirac_fin = "CRITICAL"
    else:
        dirac_fin = "UNSTABLE"

    # SO(3N): always fine
    so3n = "Well-defined"

    print(f"  {Z:>5}  {nuc:>12}  {dirac_pt:>12}  {dirac_fin:>14}  {so3n:>12}")

print(f"""
  The three limits give different "ends" to the periodic table:
  - Nuclear physics: Z ~ 126 (with possible island to ~130)
  - Dirac (point nucleus): Z ~ 137
  - Dirac (finite nucleus): Z ~ 170
  - GeoVac topology: NO LIMIT

  In practice: Nuclear instability ends the periodic table first.
  The Dirac dive is a theoretical limit that may never be reached
  because nuclei fall apart before electrons become supercritical.
""")


# =========================================================================
# STEP 11: mu_free/N^2 scaling and asymptotic analysis
# =========================================================================

print("=" * 90)
print("STEP 11: ASYMPTOTIC SCALING -- mu_free/N^2 AND OTHER RATIOS")
print("=" * 90)

print(f"""
  How does the ground state energy scale with N for large atoms?

  mu_free = 2(N-2)^2 = 2N^2 - 8N + 8

  Ratios:
    mu_free / N^2 = 2 - 8/N + 8/N^2  -> 2  as N -> infinity
    mu_free / N   = 2N - 8 + 8/N     -> 2N  (linear growth)
    mu_free / (e-e pairs) = 2(N-2)^2 / [N(N-1)/2]  -> 4  as N -> infinity

  The "Pauli cost per electron pair" saturates at 4.
  Each antisymmetric pair contributes exactly 4 units of angular
  momentum in the large-N limit.
""")

print(f"  {'N':>5}  {'mu/N^2':>8}  {'mu/pairs':>10}  {'nu/(3N-2)':>10}  "
      f"{'3/(3N-1)':>10}  {'mu*a^2/N^2':>12}")
print("-" * 65)

for N in list(range(2, 11)) + [18, 36, 54, 86, 118, 137, 170, 200, 500, 1000]:
    nu = N - 2
    mu = 2 * (N - 2)**2
    pairs = N * (N - 1) // 2
    r1 = mu / N**2
    r2 = mu / pairs if pairs > 0 else 0
    r3 = nu / (3*N - 2)
    r4 = 3.0 / (3*N - 1)  # fractional dim loss
    r5 = mu * ALPHA**2 / N**2
    print(f"  {N:>5}  {r1:>8.4f}  {r2:>10.4f}  {r3:>10.6f}  "
          f"{r4:>10.6f}  {r5:>12.8f}")

print(f"""
  ASYMPTOTIC LIMITS:
    mu/N^2    -> 2          (Pauli cost scales as 2N^2)
    mu/pairs  -> 4          (4 units per antisymmetric pair)
    nu/(3N-2) -> 1/3        (ground state at 1/3 of angular scale)
    3/(3N-1)  -> 0          (ionization becomes negligible perturbation)
    mu*a^2/N^2 -> 2*alpha^2 = {2*ALPHA**2:.8f}  (relativistic correction scale)

  The last ratio is interesting: mu * alpha^2 / N^2 measures the
  SIZE of relativistic corrections to the topological Casimir.
  It converges to 2*alpha^2 ~ 1.06e-4. This is SMALL but nonzero.

  At N = 137: mu*a^2/N^2 = {2*(137-2)**2 * ALPHA**2 / 137**2:.6f}
  At N = 1000: mu*a^2/N^2 = {2*(1000-2)**2 * ALPHA**2 / 1000**2:.6f}

  The relativistic correction per topological unit is O(alpha^2) ~ 5e-5,
  which is the SAME order as fine structure splitting.
  This is not a coincidence -- it IS fine structure.
""")


# =========================================================================
# STEP 12: Visualization data -- mu_free/N^2 vs N
# =========================================================================

print("=" * 90)
print("STEP 12: TABULATED DATA FOR PLOTTING")
print("=" * 90)

print(f"""
  DATA TABLE: mu_free/N^2 vs N (for plotting convergence to 2)
  =============================================================
""")

N_plot = np.arange(2, 201)
mu_plot = 2 * (N_plot - 2)**2
ratio_plot = mu_plot / N_plot**2

# Find where ratio crosses specific thresholds
for threshold in [1.5, 1.8, 1.9, 1.95, 1.99]:
    idx = np.searchsorted(ratio_plot, threshold)
    if idx < len(N_plot):
        print(f"  mu/N^2 crosses {threshold:.2f} at N = {N_plot[idx]}")

print()

# ASCII plot of mu/N^2 vs N
print("  mu_free / N^2 vs N")
print("  " + "-" * 72)

width = 60
y_min, y_max = 0.0, 2.1
N_ascii = list(range(2, 201, 5))

for N in N_ascii:
    mu = 2 * (N - 2)**2
    r = mu / N**2
    pos = int((r - y_min) / (y_max - y_min) * width)
    pos = max(0, min(width - 1, pos))
    bar = " " * pos + "*"

    # Mark special N values
    label = ""
    if N == 2:
        label = " He"
    elif N == 12:
        label = " Mg"
    elif N == 52:
        label = " Te"
    elif N == 117:
        label = " Ts"
    elif N == 137:
        label = " <-- Z=1/a"
    elif N == 172:
        label = " <-- Dirac (fin.)"
    elif N == 197:
        label = " Au"

    print(f"  N={N:>3d} |{bar:<{width}}|{label}")

print(f"  {'':>6}|{'0':.<{width//3}}{'1':.<{width//3}}{'2':.<{width//3+width%3}}|")
print(f"  {'':>6} 0{'':>{width//3-1}}1{'':>{width//3-1}}2")
print()

# Second ASCII plot: nu/(3N-2) approaching 1/3
print("  nu / (3N-2) vs N  (approaches 1/3 = 0.3333...)")
print("  " + "-" * 72)

y_min2, y_max2 = 0.0, 0.35
for N in N_ascii:
    nu = N - 2
    r = nu / (3*N - 2)
    pos = int((r - y_min2) / (y_max2 - y_min2) * width)
    pos = max(0, min(width - 1, pos))
    bar = " " * pos + "*"
    label = ""
    if N == 137:
        label = " <-- Z=1/a"
    print(f"  N={N:>3d} |{bar:<{width}}|{label}")

# Draw 1/3 line
third_pos = int((1/3 - y_min2) / (y_max2 - y_min2) * width)
ruler = list(" " * width)
ruler[third_pos] = "|"
print(f"  {'':>6}|{''.join(ruler)}| <-- 1/3 limit")
print()


# =========================================================================
# STEP 13: The Dirac metric correction on S^(3N-1)
# =========================================================================

print("=" * 90)
print("STEP 13: DIRAC METRIC CORRECTION -- QUANTITATIVE ESTIMATE")
print("=" * 90)

print(f"""
  The non-relativistic Fock projection gives the ROUND metric on S^3:
    ds^2 = d_chi^2 + sin^2(chi) * d_Omega_2^2

  The Dirac equation modifies this to:
    ds^2 = [1 + delta(chi, Z)] * ds^2_round

  where the correction delta depends on Z*alpha and the polar angle chi
  (which maps to radial distance via r = p0 * tan(chi/2)).

  Near the nuclear pole (chi -> 0, r -> 0):
    delta ~ -(Z*alpha)^2 / (1 + (Z*alpha)^2 * log(chi)^2 / ...)

  This is the "nuclear cusp" in geometric language.

  For a multi-electron atom on S^(3N-1), each electron has its own
  polar angle chi_i, and the total correction is:

    delta_total = sum_i delta(chi_i, Z)

  The cusp contribution from the innermost electrons (1s) dominates.

  CRITICAL OBSERVATION:
  For Z < 50: delta_1s < 0.1 -> perturbative, mu_free formula is fine
  For Z = 79 (Au): delta_1s ~ 0.3 -> significant, corrections needed
  For Z = 118 (Og): delta_1s ~ 0.7 -> strongly non-perturbative
  For Z = 137: delta_1s -> infinity -> metric singularity

  The "metric singularity" at Z*alpha = 1 IS the geometric avatar
  of the Dirac dive. The round S^3 develops a conical singularity
  at the nuclear pole, and the Casimir formula breaks down.

  QUANTITATIVE ESTIMATES:
""")

print(f"  {'Z':>5}  {'(Za)^2':>10}  {'delta_1s':>12}  {'Regime':>25}  {'mu correction':>18}")
print("-" * 85)

for Z in [1, 10, 29, 47, 79, 92, 100, 112, 118, 130, 137]:
    Za = Z * ALPHA
    Za2 = Za**2

    if Za2 < 1.0:
        # Approximate delta for 1s: (Za)^2 / (1 - (Za)^2) * leading_coeff
        # More precisely, the Dirac correction to binding energy is:
        # Delta_E / E_NR ~ (Za)^2 [1/(j+1/2) - 3/(4n)] / n
        # For 1s (n=1, j=1/2): Delta_E/E ~ (Za)^2 * (1 - 3/4) = (Za)^2/4
        # But the METRIC correction is larger: delta ~ (Za)^2 / (1 - (Za)^2)
        delta = Za2 / (1.0 - Za2)

        if delta < 0.01:
            regime = "Non-relativistic"
        elif delta < 0.1:
            regime = "Mildly relativistic"
        elif delta < 0.5:
            regime = "Relativistic"
        elif delta < 2.0:
            regime = "STRONGLY relativistic"
        else:
            regime = "NEAR-SINGULAR"

        # Correction to mu_free: delta_mu ~ (Za)^2 * mu_free (order of magnitude)
        # For neutral atom N = Z:
        mu = 2 * (Z - 2)**2
        delta_mu = Za2 * mu  # rough estimate
        print(f"  {Z:>5}  {Za2:>10.6f}  {delta:>12.6f}  {regime:>25}  {delta_mu:>18.1f}")
    else:
        print(f"  {Z:>5}  {Za2:>10.6f}  {'SINGULAR':>12}  {'BEYOND DIRAC LIMIT':>25}  {'DIVERGENT':>18}")

print(f"""
  INTERPRETATION:
  The "metric correction" delta_1s = (Za)^2/(1-(Za)^2) measures how
  much the 1s orbital warps the sphere away from round geometry.

  At Z = 79 (Au): delta = 0.50. The sphere is noticeably warped.
    This is why gold is YELLOW -- relativistic 6s contraction changes
    the absorption spectrum. The topology "knows" via metric warping.

  At Z = 118 (Og): delta = 2.87. The sphere is STRONGLY deformed.
    The round-sphere Casimir formula is no longer reliable.
    Relativistic corrections to mu_free are ~ 200% (!!).

  At Z = 137: delta -> infinity. CONICAL SINGULARITY.
    The metric degenerates at the nuclear pole.
    The sphere "pinches" -- this IS the Dirac dive geometrically.

  GEOMETRIC PICTURE OF THE DIRAC DIVE:
  =====================================
  As Z -> 1/alpha:
  1. The round S^3 develops a conical singularity at chi = 0
  2. The singularity is a point where the sphere "pinches" to zero size
  3. This creates a topology change: S^3 -> S^3 # (point removed)
     which is topologically R^3 (the sphere degenerates)
  4. The electron "falls through" the pinch point
  5. In QED: this is spontaneous pair creation from the vacuum
  6. In GeoVac: this would be a TOPOLOGY CHANGE of the configuration space
""")


# =========================================================================
# FINAL ASSESSMENT
# =========================================================================

print("=" * 90)
print("FINAL ASSESSMENT: SUPERHEAVY ELEMENTS AND THE DIRAC LIMIT")
print("=" * 90)

print(f"""
  QUESTION 1: Does anything special happen at N = 137?
  ====================================================
  ANSWER: NO. The SO(3N) representation theory is perfectly smooth
  at N = 137. The formulas nu = N-2, mu_free = 2(N-2)^2, and
  nu/(3N-2) -> 1/3 show no singularity, no phase transition,
  and no special structure at N = 137.

  The number 137 is special for the DIRAC EQUATION (Za = 1),
  not for the SO(3N) CASIMIR (non-relativistic topology).

  QUESTION 2: Is there a maximum N for the mu_free formula?
  =========================================================
  ANSWER: NO. The formula mu_free = 2(N-2)^2 is valid for all
  integer N >= 1. It is a purely combinatorial result from the
  S_N representation theory (spatial Young diagram) and the SO(3N)
  Casimir formula. Neither has a finite-N breakdown.

  The formula may need RELATIVISTIC CORRECTIONS for Z > 50-70,
  where (Za)^2 becomes non-negligible. But it does not BREAK.

  QUESTION 3: Can we see the Dirac instability geometrically?
  ===========================================================
  ANSWER: YES, but only by adding relativistic corrections to
  the GeoVac metric.

  The Dirac instability appears as a CONICAL SINGULARITY on S^(3N-1):
  - The round-sphere metric acquires a correction delta ~ (Za)^2/(1-(Za)^2)
  - At Za -> 1: delta -> infinity, the metric degenerates
  - The sphere "pinches" at the nuclear pole
  - This is a topology change: S^(3N-1) loses its smooth structure
  - Physically: an electron can fall through the pinch (pair creation)

  This geometric picture is CONSISTENT with the QED prediction
  but requires going beyond the non-relativistic framework.

  QUESTION 4: GeoVac prediction for the "end of the periodic table"?
  ==================================================================
  ANSWER: The non-relativistic GeoVac framework predicts NO end.
  The periodic table extends indefinitely in the non-relativistic limit.

  In reality, three effects terminate the periodic table:
  1. Nuclear instability: Z ~ 126 (island of stability)
  2. Dirac instability (point nucleus): Z ~ 137
  3. Dirac instability (finite nucleus): Z ~ 170
  4. GeoVac topology: no limit

  The PHYSICAL end is set by nuclear physics (Z ~ 126-130),
  not by electronic structure.

  SYNTHESIS:
  ==========
  The SO(3N) hyperspherical framework captures the COMBINATORIAL
  structure of the periodic table (S_N irreps, structure types,
  Pauli excitation) but NOT its RELATIVISTIC termination.

  To describe the Dirac limit geometrically, one would need:
  1. Replace the round metric on S^(3N-1) with a warped metric
  2. Include spin-orbit coupling (SO(3) -> SO(3,1) or SL(2,C))
  3. Allow the topology to change (pinch at Za = 1)

  This would be a natural extension: "Relativistic GeoVac" on a
  LORENTZIAN S^(3N-1) manifold. The Casimir formula would acquire
  (Za)^2 corrections, and mu_free would diverge at the critical Z.

  But this is a research program for the future, not a calculation
  we can do today.

  KEY TAKEAWAY:
  =============
  The non-relativistic topology is UNIVERSAL and FEATURELESS at N=137.
  The Dirac limit is a METRIC effect (curvature singularity), not a
  TOPOLOGICAL effect (representation theory). Topology is the Level 1
  skeleton; relativity warps the Level 2 flesh on that skeleton.
""")

print("=" * 90)
print("END OF ANALYSIS")
print("=" * 90)
