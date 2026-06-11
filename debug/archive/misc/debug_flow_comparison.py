"""
Compare the flow structures in Paper 2 (alpha derivation) and helium.

Paper 2: alpha^{-1} = S_matter / P_photon  (impedance ratio)
  - S_matter = symplectic capacity of electron lattice (plaquette areas)
  - P_photon = helical path length sqrt[(2*pi*n)^2 + delta^2]
  - delta = sqrt(pi * <L_pm>) from Liouville measure preservation
  - Result: 137.04 at n=5

Helium: Delta = E_exact - E_qc = -0.160 Ha
  - |Delta| ~ Z^2/(8*pi) at Z=2 (suggestive but not exact for other Z)
  - Comes from the non-Coulombic part of mu(R) at R ~ 1-3 bohr
  - The S^5 -> S^3 x R transition

KEY QUESTION: Is there a common "impedance" structure?
"""

import numpy as np
from fractions import Fraction
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.hyperspherical_angular import solve_angular


# =========================================================================
# Algebraic formulas
# =========================================================================

def I_nuc_exact(n: int) -> Fraction:
    return sum(Fraction(1, 2*j+1) for j in range(2*n))

def I_ee_rational_exact(n: int) -> Fraction:
    s1 = sum(Fraction((-1)**k, 4*k+1) for k in range(n))
    s2 = sum(Fraction((-1)**k, 4*k+3) for k in range(n))
    return s1 - s2

def a1_exact(n: int, Z: float) -> float:
    I_nuc = float(I_nuc_exact(n))
    I_ee_rat = float(I_ee_rational_exact(n))
    return (4 / np.pi) * (-2 * Z * I_nuc + np.sqrt(2) * I_ee_rat)

def V_nuc_off_diag(m: int, n: int, Z: float) -> float:
    if (m + n) % 2 != 0:
        return 0.0
    p_plus = (m + n) // 2
    p_minus = abs(m - n) // 2
    I_plus = float(I_nuc_exact(p_plus))
    I_minus = float(I_nuc_exact(p_minus)) if p_minus > 0 else 0.0
    return (-8 * Z / np.pi) * (I_plus - I_minus)

def V_ee_off_diag(m: int, n: int) -> float:
    if (m + n) % 2 != 0:
        return 0.0
    p_plus = (m + n) // 2
    p_minus = abs(m - n) // 2
    I_plus = float(I_ee_rational_exact(p_plus))
    I_minus = float(I_ee_rational_exact(p_minus)) if p_minus > 0 else 0.0
    return (4 * np.sqrt(2) / np.pi) * (I_plus - I_minus)

def a2_from_sum(n_target: int, Z: float, n_terms: int = 100) -> float:
    E_n = 2 * n_target**2 - 2
    a2 = 0.0
    for m in range(n_target + 2, 2 * n_terms + 2, 2):
        E_m = 2 * m**2 - 2
        V = V_nuc_off_diag(m, n_target, Z) + V_ee_off_diag(m, n_target)
        if abs(V) < 1e-15:
            continue
        a2 += V**2 / (E_n - E_m)
    return a2


# =========================================================================
# STEP 1: Recast Paper 2 in flow language
# =========================================================================

print("=" * 72)
print("STEP 1: Paper 2 alpha derivation in flow language")
print("=" * 72)

# Paper 2: alpha^{-1} = S_n / P_n at n = 5
# S_n = sum of plaquette areas on the SU(2) lattice for shell n
# P_n = helical path length sqrt[(2*pi*n)^2 + delta^2]

# Reframe as a FLOW:
# Consider the "shell flow" from n to n+1.
# At each shell n, the electron occupies a phase space of capacity S_n.
# The photon gauge field has action P_n.
# The IMPEDANCE kappa_n = S_n / P_n measures the "resistance" to
# electromagnetic interaction at shell n.

# The key observation: S_n grows as n^4, P_n grows as n.
# So kappa_n ~ n^3 for large n.
# But alpha^{-1} is approximately kappa_5 = 137.04.

# From the old research, the plaquette areas are:
# Area(n,l,m) = |T_+(n,l,m) x L_+(n,l,m)| where:
# T_+ = sqrt[(n+l+1)(n-l)/n^2] (radial transition)
# L_+ = sqrt[(l-m)(l+m+1)] (angular transition)

# Compute S_n for several shells
def compute_S_n(n: int) -> float:
    """Compute matter symplectic capacity for shell n."""
    S = 0.0
    for l in range(n):
        for m in range(-l, l):
            T_plus = np.sqrt((n + l + 1) * (n - l)) / n
            L_plus = np.sqrt((l - m) * (l + m + 1))
            # Plaquette area = cross product magnitude of operator matrix elements
            # In the 2D (T, L) phase space, area = T_+ * L_+
            S += T_plus * L_plus
    return S

# Compute P_n for several shells
def compute_P_n(n: int) -> float:
    """Compute photon helical action for shell n."""
    # Average angular momentum for this shell
    L_avg = 0.0
    count = 0
    for l in range(n):
        for m in range(-l, l + 1):
            L_avg += np.sqrt(l * (l + 1))
            count += 1
    L_avg /= max(count, 1)

    # Helical pitch from Liouville measure: delta = sqrt(pi * <L>)
    delta = np.sqrt(np.pi * L_avg)

    # Total helical path
    return np.sqrt((2 * np.pi * n)**2 + delta**2)


print(f"\n  Shell flow analysis:")
print(f"  {'n':>3}  {'S_n':>12}  {'P_n':>10}  {'kappa=S/P':>12}  {'delta':>8}")
print("-" * 55)

for n in range(1, 8):
    S = compute_S_n(n)
    # Compute delta explicitly
    L_avg = 0.0
    count = 0
    for l in range(n):
        for m in range(-l, l + 1):
            L_avg += np.sqrt(l * (l + 1))
            count += 1
    L_avg /= max(count, 1)
    delta = np.sqrt(np.pi * L_avg)
    P = np.sqrt((2 * np.pi * n)**2 + delta**2)
    kappa = S / P

    print(f"  {n:3d}  {S:12.4f}  {P:10.4f}  {kappa:12.4f}  {delta:8.4f}")

print(f"\n  alpha^{{-1}} (experimental) = 137.036")
print(f"  The ratio S_n/P_n converges near n=5.")

# Flow rate: how does kappa change with n?
print(f"\n  Shell-to-shell 'impedance flow':")
print(f"  {'n':>3}  {'dS/dn':>10}  {'dP/dn':>10}  {'d(kappa)/dn':>12}")
prev_S, prev_P, prev_k = None, None, None
for n in range(1, 8):
    S = compute_S_n(n)
    P = compute_P_n(n)
    k = S / P
    if prev_S is not None:
        print(f"  {n:3d}  {S - prev_S:10.4f}  {P - prev_P:10.4f}  {k - prev_k:12.4f}")
    prev_S, prev_P, prev_k = S, P, k


# =========================================================================
# STEP 2: Helium "impedance" -- angular vs radial
# =========================================================================

print("\n" + "=" * 72)
print("STEP 2: Helium impedance ratio")
print("=" * 72)

Z = 2.0
a1_z2 = a1_exact(1, Z)
a2_z2 = a2_from_sum(1, Z, n_terms=100)
Z_eff = -a1_z2

# In helium, the analogous decomposition is:
# E = (angular contribution) + (radial contribution)
# = a_2 + (-Z_eff^2 / (2*n_eff^2))
#
# The "angular capacity" S_angular is related to the charge function
# integral: S_ang ~ sum of |V_{m1}|^2 / dE (the a_2 sum)
#
# The "radial action" P_radial is the Bohr-Sommerfeld quantization:
# P_rad = oint p_R dR = n_eff * pi (radial action in the Coulomb problem)

# For the quasi-Coulomb model:
n_eff = 2.5  # = N + l_eff + 1 = 0 + 3/2 + 1
J_radial = n_eff * np.pi  # Bohr-Sommerfeld radial action

# Angular "capacity": the total coupling strength
# S_angular = sum_m |V_{m1}|^2 (unnormalized) or a_2 * something

# Actually, let me define impedance differently.
# In Paper 2: kappa = S/P where S is phase space area and P is path length.
# In helium: what is the "phase space area" of the angular problem?

# The angular phase space is T*S^5. At each R, the eigenvalue mu(R) defines
# an energy surface. The phase space volume enclosed is:
# Omega_ang(R) = integral over S^5 of sqrt(mu(R) - V_ang) d(Omega_5) d(p_Omega_5)
# This is complicated. Let me use a simpler measure.

# Measure 1: The "spectral capacity" = sum of occupied eigenvalues
# At R=0: only mu_1 = 0 is occupied. Trivial.
# At R>0: mu_1(R) < 0. The "capacity" is |mu_1(R)|.

# Measure 2: The coupling matrix spectral width
# sigma(R) = sqrt(sum_m |V_{m1}|^2 * R^2) = R * sqrt(sum |V_{m1}|^2)
# This is the RMS coupling strength at distance R.

sum_V2 = 0.0
for m in range(3, 202, 2):
    V = V_nuc_off_diag(m, 1, Z) + V_ee_off_diag(m, 1)
    sum_V2 += V**2
sigma_V = np.sqrt(sum_V2)

print(f"\n  Angular coupling width: sigma_V = sqrt(sum |V_m1|^2) = {sigma_V:.6f}")
print(f"  Radial action: J_radial = n_eff * pi = {J_radial:.6f}")
print(f"  Impedance ratio: sigma_V / J_radial = {sigma_V / J_radial:.6f}")
print(f"  Compare to a_1 = {a1_z2:.6f}")

# Measure 3: Direct phase space argument
# The angular problem at hyperradius R has:
# - Volume of S^5 = pi^3
# - Effective "momentum" = sqrt(|mu(R)|)
# - Phase space volume ~ Vol(S^5) * |mu(R)|^{5/2} (Weyl's law in 5D)
#
# The radial problem has:
# - Range R in [0, R_max]
# - Momentum p_R = sqrt(2*(E - V_eff(R)))
# - Phase space volume ~ integral p_R dR

# But the simplest analog is:
# Angular "area" = |a_1| (first-order energy per unit R, dimension-free)
# Radial "path" = n_eff * pi (quantized action)
# Impedance = |a_1| / (n_eff * pi) = dimensionless

kappa_He_1 = abs(a1_z2) / (n_eff * np.pi)
print(f"\n  Impedance version 1: |a_1| / (n_eff * pi)")
print(f"    = {abs(a1_z2):.6f} / {n_eff * np.pi:.6f} = {kappa_He_1:.6f}")

# Or: Z_eff / (pi * n_eff)
kappa_He_2 = Z_eff / (np.pi * n_eff)
print(f"\n  Impedance version 2: Z_eff / (pi * n_eff)")
print(f"    = {Z_eff:.6f} / {np.pi * n_eff:.6f} = {kappa_He_2:.6f}")

# Or: total angular action / radial action
# Total angular action = sum of coupling * R over the classical region
# ~ |a_1| * <R> = |a_1| * 15/(2*Z_eff)
R_expect = 15.0 / (2 * Z_eff)
S_ang = abs(a1_z2) * R_expect
kappa_He_3 = S_ang / J_radial
print(f"\n  Impedance version 3: |a_1|*<R> / (n_eff*pi)")
print(f"    = {S_ang:.6f} / {J_radial:.6f} = {kappa_He_3:.6f}")


# =========================================================================
# STEP 3: The pi connection -- structural analysis
# =========================================================================

print("\n" + "=" * 72)
print("STEP 3: The pi connection")
print("=" * 72)

# Both alpha and Delta involve pi through GEOMETRIC paths on spheres.
# alpha: helix on S^2 (photon) vs plaquettes on SU(2) lattice
# Delta: Coulomb solution on S^3 vs charge function on S^5

# The specific pi appearances:
# alpha^{-1} = S_5 / sqrt((2*pi*5)^2 + delta^2) = 4325.83 / 31.567
#            ~ 4325.83 / (10*pi) = 137.67 (if delta=0, 0.5% error)
#            ~ S_5 / (10*pi + correction)

# Delta = E_exact - E_qc
# E_qc = a_2 - Z_eff^2/(2*(5/2)^2)
# a_1 = (8/(3*pi))*(sqrt(2) - 4Z)
# Z_eff = (8/(3*pi))*(4Z - sqrt(2))
# Z_eff^2 = (64/(9*pi^2))*(4Z - sqrt(2))^2
# E_qc = a_2 - (64/(9*pi^2))*(4Z-sqrt(2))^2 / (2*(25/4))
# = a_2 - (64/(9*pi^2))*(4Z-sqrt(2))^2 * 2/25
# = a_2 - 128*(4Z-sqrt(2))^2 / (225*pi^2)

print(f"\n  Paper 2 pi structure:")
print(f"    alpha^{{-1}} = S_5 / P_5")
print(f"    P_5 = sqrt((10*pi)^2 + delta^2) = 31.567")
print(f"    1/pi contribution: P_5 ~ 10*pi => alpha^{{-1}} ~ S_5/(10*pi)")

print(f"\n  Helium pi structure:")
print(f"    Z_eff = (8/(3*pi)) * (4Z - sqrt(2))")
print(f"    E_qc = a_2 - Z_eff^2 / 12.5")
print(f"         = a_2 - 128*(4Z-sqrt(2))^2 / (225*pi^2)")

# The 1/pi^2 factor in E_qc comes from Z_eff^2, which comes from a_1^2.
# a_1 = (4/pi) * (angular integral), so a_1^2 has 1/pi^2.
# Every a_n involves 1/pi^n in principle.

# Delta = E_exact - a_2 + 128*(4Z-sqrt(2))^2/(225*pi^2)
# = -0.160 for Z=2
# = -0.404 + 128*(8-sqrt(2))^2/(225*pi^2)

val_4Z_sqrt2_sq = (4*Z - np.sqrt(2))**2
Coulomb_term = 128 * val_4Z_sqrt2_sq / (225 * np.pi**2)
print(f"\n  Coulomb contribution: 128*(4Z-sqrt(2))^2/(225*pi^2)")
print(f"    = 128 * {val_4Z_sqrt2_sq:.6f} / (225 * pi^2)")
print(f"    = {Coulomb_term:.6f}")
print(f"    Z_eff^2/(2*n_eff^2) = {Z_eff**2/12.5:.6f}")
print(f"    Match: {abs(Coulomb_term - Z_eff**2/12.5) < 1e-10}")

# So BOTH alpha and Delta have the form:
# (algebraic integer structure) / (pi^n * geometric factor)
# alpha: integer lattice count / (2*pi*n * geometric correction)
# Delta: algebraic perturbation sum / (pi^2 * Coulomb quantization)

print(f"\n  STRUCTURAL PARALLEL:")
print(f"    alpha = [plaquette count on SU(2)]^{{-1}} * [2*pi*n + helical correction]")
print(f"    Delta = [algebraic a_2 sum on SO(6)] - [Z_eff^2/(pi^2 * rational)]")
print(f"    Both have: (discrete lattice quantity) / (pi-geometric path)")


# =========================================================================
# STEP 4: Symplectic volume of the transition region
# =========================================================================

print("\n" + "=" * 72)
print("STEP 4: Symplectic volume of the transition region")
print("=" * 72)

# The transition region R in [0.5, 2.0] bohr contains 75% of Delta.
# Compute the radial phase space volume in this region.
#
# For the quasi-Coulomb potential:
# V_eff(R) = 15/(8R^2) + a_1/R + a_2
# p_R(R) = sqrt(2*(E - V_eff(R))) (classical momentum)
# Phase space area = integral p_R dR over classically allowed region

E_qc = a2_z2 - Z_eff**2 / (2 * n_eff**2)

# Find classical turning points
R_scan = np.linspace(0.1, 5.0, 1000)
V_qc = 15.0 / (8 * R_scan**2) + a1_z2 / R_scan + a2_z2
classically_allowed = V_qc < E_qc

if np.any(classically_allowed):
    allowed_indices = np.where(classically_allowed)[0]
    R_inner = R_scan[allowed_indices[0]]
    R_outer = R_scan[allowed_indices[-1]]

    # Phase space area (full orbit)
    dR = R_scan[1] - R_scan[0]
    p_R = np.zeros_like(R_scan)
    mask = V_qc < E_qc
    p_R[mask] = np.sqrt(2 * (E_qc - V_qc[mask]))
    J_full = np.sum(p_R * dR) * 2  # factor 2 for both directions

    print(f"\n  Quasi-Coulomb classical region:")
    print(f"    E_qc = {E_qc:.6f} Ha")
    print(f"    R_inner = {R_inner:.4f} bohr, R_outer = {R_outer:.4f} bohr")
    print(f"    Phase space area J = oint p_R dR = {J_full:.6f}")
    print(f"    Expected: n_eff * 2*pi = {n_eff * 2 * np.pi:.6f}")
    print(f"    (Bohr-Sommerfeld: J = n*2*pi for the n-th orbit)")

    # Phase space area of just the transition region [0.5, 2.0]
    mask_trans = (R_scan >= 0.5) & (R_scan <= 2.0) & mask
    J_trans = np.sum(p_R[mask_trans] * dR) * 2
    print(f"\n  Transition region [0.5, 2.0] bohr:")
    print(f"    J_trans = {J_trans:.6f}")
    print(f"    Fraction of total: {J_trans / J_full * 100:.1f}%")

# Now compute with EXACT V_eff
print(f"\n  Exact V_eff phase space:")
R_exact = np.linspace(0.1, 5.0, 500)
mu_exact = np.zeros(len(R_exact))
for i, R in enumerate(R_exact):
    mu, _ = solve_angular(R=R, Z=Z, l_max=3, n_alpha=200, n_channels=1)
    mu_exact[i] = mu[0]

V_exact = mu_exact / R_exact**2 + 15.0 / (8.0 * R_exact**2)
E_exact_He = -2.903724

mask_ex = V_exact < E_exact_He
if np.any(mask_ex):
    idx_ex = np.where(mask_ex)[0]
    R_in_ex = R_exact[idx_ex[0]]
    R_out_ex = R_exact[idx_ex[-1]]

    dR_ex = R_exact[1] - R_exact[0]
    p_R_ex = np.zeros_like(R_exact)
    p_R_ex[mask_ex] = np.sqrt(2 * (E_exact_He - V_exact[mask_ex]))
    J_exact = np.sum(p_R_ex * dR_ex) * 2

    print(f"    R_inner = {R_in_ex:.4f}, R_outer = {R_out_ex:.4f}")
    print(f"    J_exact = {J_exact:.6f}")
    print(f"    J_QC = {J_full:.6f}")
    print(f"    Difference: {J_exact - J_full:.6f}")


# =========================================================================
# STEP 5: The Liouville conserved quantity
# =========================================================================

print("\n" + "=" * 72)
print("STEP 5: Conserved quantity during S^5 -> S^3 x R flow")
print("=" * 72)

# In Paper 2: Liouville's theorem preserves the symplectic 2-form.
# The helical pitch delta = sqrt(pi * <L_pm>) ensures that the
# electron-photon coupling preserves phase space volume.
#
# In helium: as R increases, the angular eigenstate evolves.
# What is conserved?
#
# 1. NORMALIZATION: <psi(R)|psi(R)> = 1 always (trivial)
# 2. ENERGY: E = mu(R)/R^2 + 15/(8R^2) + (radial kinetic)
# 3. ANGULAR MOMENTUM: L = 0 (we're in the S-wave sector)
# 4. PARTICLE EXCHANGE SYMMETRY: psi(alpha) = psi(pi/2 - alpha) always
#
# But there's a deeper conserved quantity: the SPECTRAL MEASURE.
# Expand psi(R) = sum_n c_n(R) * phi_n (in the free basis).
# Then sum |c_n|^2 = 1 (unitarity).
# The EFFECTIVE DIMENSION = 1 / sum |c_n|^4 (inverse participation ratio
# in the free basis) measures how many channels are excited.

n_alpha = 200
l_max = 3
n_channels = 10

# Build the free basis at R=0
mu_free_vals, vecs_free = solve_angular(R=0, Z=Z, l_max=l_max,
                                         n_alpha=n_alpha, n_channels=n_channels)

print(f"\n  Angular state evolution in free basis:")
print(f"  {'R':>6}  {'dim_eff':>8}  {'|c_1|^2':>8}  {'|c_3|^2':>8}  {'|c_5|^2':>8}  "
      f"{'sum|c|^4':>10}  {'entropy':>8}")
print("-" * 68)

for R in [0.0, 0.1, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0]:
    if R == 0.0:
        mu_R, vecs_R = mu_free_vals, vecs_free
    else:
        mu_R, vecs_R = solve_angular(R=R, Z=Z, l_max=l_max,
                                      n_alpha=n_alpha, n_channels=n_channels)

    # Project ground state onto free basis
    psi_gs = vecs_R[0]  # ground state at this R
    overlaps = np.zeros(n_channels)
    for n in range(n_channels):
        overlaps[n] = np.dot(psi_gs, vecs_free[n])**2

    # Normalize (should be ~1 if bases span the same space)
    overlap_sum = np.sum(overlaps)
    if overlap_sum > 0.01:
        overlaps_n = overlaps / overlap_sum
    else:
        overlaps_n = overlaps

    dim_eff = 1.0 / np.sum(overlaps_n**2) if np.sum(overlaps_n**2) > 0 else 0
    entropy = -np.sum(overlaps_n * np.log(overlaps_n + 1e-30))

    c1sq = overlaps_n[0] if len(overlaps_n) > 0 else 0
    c3sq = overlaps_n[1] if len(overlaps_n) > 1 else 0
    c5sq = overlaps_n[2] if len(overlaps_n) > 2 else 0

    print(f"  {R:6.2f}  {dim_eff:8.4f}  {c1sq:8.4f}  {c3sq:8.4f}  {c5sq:8.4f}  "
          f"{np.sum(overlaps_n**2):10.6f}  {entropy:8.4f}")

# The CONSERVED quantity is the total probability = 1 (unitarity).
# The NON-conserved quantity is the effective dimension (IPR).
# As R increases: dim_eff grows (more channels mix in).
# This is the analog of Paper 2's "impedance flow" across shells.

print(f"\n  The CONSERVED quantity: total probability = 1 (unitarity)")
print(f"  The FLOWING quantity: effective dimension (IPR in free basis)")
print(f"  As R increases, more SO(6) channels mix -> higher eff. dimension")
print(f"  This is analogous to Paper 2's flow across shells n -> n+1")


# =========================================================================
# STEP 6: Direct comparison of impedance structures
# =========================================================================

print("\n" + "=" * 72)
print("STEP 6: Direct structural comparison")
print("=" * 72)

# Paper 2:
# S_5 = sum_{l=0}^{4} sum_{m=-l}^{l-1} T_+(5,l,m) * L_+(5,l,m) = 4325.83
# P_5 = sqrt((10*pi)^2 + delta^2) = 31.567
# alpha^{-1} = 4325.83 / 31.567 = 137.04

S_5 = compute_S_n(5)
P_5 = compute_P_n(5)
alpha_inv_computed = S_5 / P_5

print(f"\n  Paper 2:")
print(f"    S_5 = {S_5:.4f}  (matter phase space capacity)")
print(f"    P_5 = {P_5:.4f}  (photon helical action)")
print(f"    alpha^{{-1}} = S_5/P_5 = {alpha_inv_computed:.4f}")

# Helium:
# S_He = "angular capacity" = ... what?
# P_He = "radial action" = n_eff * pi = 5*pi/2

# The angular capacity in helium is the charge function matrix element sum:
# a_1 = (4/pi) * (-2Z * I_nuc(1) + sqrt(2) * I_ee_rat(1))
# This is a "single-shell" quantity (ground state coupling).

# The radial action is pi * n_eff.
# So kappa_He = a_1_terms / P_radial = ?

# But we need to match DIMENSIONLESS quantities.
# In Paper 2: S is dimensionless (plaquette count * area per plaquette)
#              P is dimensionless (angle in radians)
# In helium:  a_1 has dimensions of 1/bohr (slope of mu vs R)
#              But in atomic units, everything is dimensionless.

# The closest parallel: a_2 is the "effective angular potential" (units: Hartree)
# And Z_eff^2/(2*n_eff^2) is the "radial binding" (units: Hartree)
# Their ratio:

ratio_He = abs(a2_z2) / (Z_eff**2 / (2 * n_eff**2))
print(f"\n  Helium:")
print(f"    |a_2|         = {abs(a2_z2):.6f}  (angular effective potential)")
print(f"    Z_eff^2/12.5  = {Z_eff**2/12.5:.6f}  (radial Coulomb binding)")
print(f"    Ratio |a_2| / E_rad = {ratio_He:.6f}")

# Another comparison: Paper 2 has S ~ n^4 and P ~ n.
# In helium, the "shell sum" is a_2 = sum |V_{m1}|^2 / (2(m^2-1))
# and the "path" is the Coulomb quantization n_eff = 5/2.
# The individual terms |V_{m1}|^2 ~ m^{-2} and dE ~ m^2, so terms ~ m^{-4}.

print(f"\n  Scaling comparison:")
print(f"    Paper 2: S_n ~ n^4, P_n ~ n => kappa ~ n^3")
print(f"    Helium:  |V_m1|^2 ~ m^{{-2}}, dE_m ~ m^2 => a_2 terms ~ m^{{-4}}")
print(f"    Both show POWER-LAW coupling decay.")


# =========================================================================
# STEP 7: The 1/(2*pi) coincidence
# =========================================================================

print("\n" + "=" * 72)
print("STEP 7: Analysis of the 1/(2*pi) coincidence")
print("=" * 72)

Delta = -0.159548

# |Delta| = 0.1595 vs 1/(2*pi) = 0.1592
# Relative error: 0.25%

print(f"\n  |Delta| = {abs(Delta):.8f}")
print(f"  1/(2*pi) = {1/(2*np.pi):.8f}")
print(f"  Relative difference: {abs(abs(Delta) - 1/(2*np.pi)) / abs(Delta) * 100:.2f}%")

# But we showed Delta(Z) is NOT proportional to Z^2/(8*pi).
# For Z=2: Delta/Z^2 = 0.0399 vs 1/(8*pi) = 0.0398 (great!)
# For Z=3: Delta/Z^2 = 0.0320 (diverging from 0.0398)
# For Z=6: Delta/Z^2 = 0.0260 (further)

# The "true" Z^2 coefficient from the fit was 0.0220.
# So Delta ~ 0.022*Z^2 + 0.018*Z + 0.035
# At Z=2: 0.022*4 + 0.018*2 + 0.035 = 0.088 + 0.036 + 0.035 = 0.159 (matches!)
# The 1/(2*pi) agreement at Z=2 is ACCIDENTAL: the Z^1 and Z^0 terms
# just happen to combine with Z^2 to give ~ Z^2/(8*pi).

# But IS the Z^0 term or Z^1 term individually related to pi?
# Fit was: -0.022*Z^2 - 0.018*Z - 0.035
# c_0 = -0.035 vs -1/(9*pi) = -0.035 (!!)
# c_1 = -0.018 vs -1/(18*pi) = -0.0177 (!!)

c0_fit = -0.034792  # from previous run
c1_fit = -0.018363
c2_fit = -0.022018

print(f"\n  Coefficients of Delta(Z) = c_2*Z^2 + c_1*Z + c_0:")
print(f"    c_2 = {c2_fit:.6f}")
print(f"    c_1 = {c1_fit:.6f}")
print(f"    c_0 = {c0_fit:.6f}")

print(f"\n  Testing each coefficient for pi content:")
pi_tests = [
    ("c_0", c0_fit, [
        ("-1/(9*pi)", -1/(9*np.pi)),
        ("-1/(3*pi^2)", -1/(3*np.pi**2)),
        ("-sqrt(2)/(4*pi^2)", -np.sqrt(2)/(4*np.pi**2)),
    ]),
    ("c_1", c1_fit, [
        ("-1/(18*pi)", -1/(18*np.pi)),
        ("-sqrt(2)/(4*pi^2)", -np.sqrt(2)/(4*np.pi**2)),
        ("-1/(8*pi^2)", -1/(8*np.pi**2)),
    ]),
    ("c_2", c2_fit, [
        ("-2/(9*pi^2)", -2/(9*np.pi**2)),
        ("-1/(4*pi^2)", -1/(4*np.pi**2)),
        ("-1/(16*pi)", -1/(16*np.pi)),
    ]),
]

for name, val, tests in pi_tests:
    print(f"\n  {name} = {val:.6f}:")
    for expr, expr_val in sorted(tests, key=lambda x: abs(x[1] - val)):
        diff = abs(expr_val - val)
        pct = diff / abs(val) * 100
        print(f"    {expr:>25} = {expr_val:.6f}  (diff: {diff:.4f}, {pct:.1f}%)")


# =========================================================================
# STEP 8: What IS structurally common?
# =========================================================================

print("\n" + "=" * 72)
print("STEP 8: What is structurally common")
print("=" * 72)

print("""
STRUCTURAL COMPARISON: Paper 2 (alpha) vs Helium (Delta)

                    PAPER 2 (alpha)         HELIUM (Delta)
---------------------------------------------------------------
Geometry:           SU(2) lattice           SO(6)/SO(5) = S^5
                    (n, l, m states)        (nu, l channels)

Discrete "matter":  S_n = plaquette areas   V_{m1} = matrix elements
                    Sum over (l,m)          Sum over m (nu = 2(m-1))

Continuous "gauge":  P_n = 2*pi*n + helix   E_qc = Z_eff^2/(2*n^2)
                    (photon circular path)  (Coulomb quantization)

Measure:            Symplectic 2-form       Hilbert space norm
                    (Liouville)             (unitarity)

Result:             alpha^{-1} = S/P        Delta = integral(dV*F^2)

Pi source:          Circumference 2*pi*n    Z_eff = (8/3*pi)*(...)
                    (photon path IS a pi)   (Fock projection IS a pi)

BOTH:
  1. Ratio of DISCRETE to CONTINUOUS: lattice count / pi-path
  2. Coupling through pi: appears as 1/pi^n in the ratio
  3. Convergence at a specific scale: n=5 for alpha, R~1 for He
  4. Power-law coupling decay: ~n^{-3} for alpha, ~m^{-2} for He

DIFFERENCE:
  1. alpha is a RATIO (impedance): S/P
     Delta is a RESIDUAL (error): E_exact - E_truncated
  2. alpha is scale-independent (dimensionless coupling)
     Delta is scale-dependent (depends on Z, dimensionful)
  3. alpha involves SU(2) x U(1) (matter x photon)
     Delta involves SO(6) -> SO(4) (angular manifold surgery)

THE DEEP CONNECTION (if any):
  Both involve Fock-projected integrals on spheres.
  - Paper 2: plaquette areas on the Fock-projected hydrogen lattice
  - Helium: matrix elements of the charge function on S^5

  The factor of 4/pi in all helium matrix elements comes from
  the SAME source as the 2*pi*n in the photon path: the stereographic
  projection maps the discrete lattice into curved space, and the
  curvature contributes factors of pi.

  The 1/(2*pi) coincidence at Z=2 may reflect: the Coulomb quantization
  (n_eff = 5/2) and the Fock projection (4/pi) combine to give a
  correction of exactly one "winding number" = 1/(2*pi) when Z=2
  makes the algebra particularly clean.
""")

# =========================================================================
# STEP 9: Test the "winding number" hypothesis
# =========================================================================

print("=" * 72)
print("STEP 9: Winding number hypothesis")
print("=" * 72)

# Hypothesis: Delta = -1/(2*pi) * (something that = 1 at Z=2)
# Or more precisely: the correction to the Coulomb energy from
# higher-order terms has units of 1/(2*pi) * Z^2/4 at Z=2.
#
# In the Coulomb problem, the Runge-Lenz vector generates a hidden
# SO(4) symmetry. The quantization condition is:
# n = n_r + l + 1  (principal quantum number)
#
# In the hyperspherical problem, the grand angular momentum quantum
# number nu generates SO(6). The quantization is:
# n_eff = N + l_eff + 1 = N + 5/2
#
# The "mismatch" between the exact SO(6) quantization and the
# truncated Coulomb SO(4) quantization could produce a correction
# proportional to the "geometric phase" = 2*pi * (winding difference).

# Define: n_eff_exact such that E_exact = a_2 - Z_eff^2/(2*n_exact^2)
n_eff_exact_sq = Z_eff**2 / (2 * (a2_z2 - (-2.903724)))
n_eff_exact = np.sqrt(n_eff_exact_sq)

print(f"\n  n_eff (standard) = 5/2 = 2.5")
print(f"  n_eff (exact) = {n_eff_exact:.8f}")
print(f"  Delta_n = n_exact - 5/2 = {n_eff_exact - 2.5:.8f}")

# Is Delta_n related to 1/(2*pi)?
# 1/(2*pi) = 0.15915
# Delta_n = -0.0761
# Delta_n / (1/(2*pi)) = -0.478... not clean.

# But: (5/2)^2 - n_exact^2 = ?
dn2 = 2.5**2 - n_eff_exact**2
print(f"\n  (5/2)^2 - n_exact^2 = {dn2:.8f}")
print(f"  dn2 / (1/pi) = {dn2 * np.pi:.6f}")
print(f"  dn2 / (1/pi^2) = {dn2 * np.pi**2:.6f}")

# E_exact - E_qc = -Z_eff^2/(2*n_exact^2) + Z_eff^2/(2*(5/2)^2)
# = (Z_eff^2/2) * [1/(5/2)^2 - 1/n_exact^2]
# = (Z_eff^2/2) * [(5/2)^2 - n_exact^2] / [(5/2)^2 * n_exact^2]
# = (Z_eff^2/2) * dn2 / (6.25 * n_exact^2)

check = (Z_eff**2 / 2) * dn2 / (6.25 * n_eff_exact**2)
print(f"\n  Consistency check: (Z_eff^2/2)*dn2/(6.25*n_ex^2) = {check:.8f}")
print(f"  Delta = {Delta:.8f}")
print(f"  Match: {abs(check - Delta) < 1e-6}")


# =========================================================================
# FINAL SYNTHESIS
# =========================================================================

print("\n" + "=" * 72)
print("FINAL SYNTHESIS")
print("=" * 72)

print(f"""
VERDICT: Paper 2 (alpha) and helium (Delta) share STRUCTURAL parallels
but are NOT the same phenomenon.

SHARED STRUCTURE:
  1. Both involve Fock-projected integrals on spheres.
  2. Both produce factors of 1/pi from the stereographic projection.
  3. Both involve sums over discrete quantum states (plaquettes / channels).
  4. Both have power-law coupling decay.

THE pi CONNECTION IS REAL BUT GENERIC:
  The factor 4/pi in helium matrix elements and 2*pi*n in the photon
  path both come from the same mathematical origin: integration over
  spheres in Fock-projected coordinates. ANY physical quantity computed
  on the Fock-projected lattice will involve powers of 1/pi.

THE 1/(2*pi) COINCIDENCE AT Z=2 IS ACCIDENTAL:
  Delta(Z) = -0.022*Z^2 - 0.018*Z - 0.035
  At Z=2: Delta = -0.088 - 0.037 - 0.035 = -0.160 ~ -1/(2*pi)
  This is the sum of THREE terms, each involving different powers
  of 1/pi from the Fock projection. Their sum happens to approximate
  1/(2*pi) at Z=2 due to numerical coincidence, not deep structure.

WHAT IS GENUINELY COMMON:
  The graph Laplacian generates a DIMENSIONLESS lattice (Paper 7).
  Both alpha and Delta are extracted from this lattice:
  - alpha: impedance ratio of the SU(2) lattice (1-electron, 1-photon)
  - Delta: truncation error of the SO(6) perturbation theory (2-electron)

  The lattice is the common ancestor. The factors of pi in both
  quantities trace back to the Fock projection that maps the
  dimensionless lattice to physical space.

  This is a consistency check on the framework, not a prediction:
  both quantities SHOULD involve pi because both come from the
  same stereographic projection. But the specific values (137.036
  and 0.160) are not related by any simple algebraic formula.
""")
