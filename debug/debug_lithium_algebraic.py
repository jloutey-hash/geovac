"""
Algebraic structure of the lithium 3-electron problem.

GOAL: Determine if Li angular eigenvalues have SO(9) Casimir structure
analogous to He's SO(6) structure, and characterize the topology flow.

KEY RESULTS:
  He (2e): S^5 manifold, SO(6) Casimir C_2 = nu(nu+4), l_eff = 3/2
  Li (3e): S^8 manifold, SO(9) Casimir C_2 = nu(nu+7), l_eff = 3

No 3-electron hyperspherical code exists. This script performs a
THEORETICAL analysis based on group theory and dimensional arguments,
comparing with the existing FCI results.
"""

import numpy as np
from fractions import Fraction
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# =========================================================================
# STEP 1: SO(d) Casimir structure for N electrons
# =========================================================================

print("=" * 72)
print("STEP 1: SO(d) Casimir structure for N-electron atoms")
print("=" * 72)

# For N electrons in 3D:
# Configuration space: R^{3N}
# Hyperradius: R = sqrt(sum r_i^2)
# Angular manifold: S^{3N-1}
# Symmetry group: SO(3N)
#
# Laplace-Beltrami eigenvalues on S^{d-1} (d = 3N):
# lambda_nu = nu(nu + d - 2) for the SO(d) symmetric representation (nu, 0, ..., 0)
#
# Free angular eigenvalues: mu_free = lambda_nu / 2 = nu(nu + d - 2) / 2
# (factor of 1/2 from the kinetic energy T = Lambda^2 / (2R^2))
#
# Centrifugal potential: (3N-1)(3N-3) / (8R^2) from the R^{-(3N-1)/2} substitution

print(f"\n  {'N_e':>4}  {'d=3N':>5}  {'S^d-1':>5}  {'SO(d)':>6}  {'C_2':>15}  "
      f"{'centrifugal':>14}  {'l_eff':>6}  {'n_eff_gs':>8}")
print("-" * 80)

for N_e in range(1, 6):
    d = 3 * N_e
    # Casimir: C_2(nu) = nu(nu + d - 2)
    casimir_str = f"nu(nu+{d-2})"
    # Centrifugal: (3N-1)(3N-3)/8
    cent_num = (3*N_e - 1) * (3*N_e - 3)
    cent = Fraction(cent_num, 8)
    # l_eff: l_eff(l_eff+1)/2 = cent => l_eff(l_eff+1) = cent_num/4
    # l_eff = (-1 + sqrt(1 + cent_num))/2
    l_eff_val = (-1 + np.sqrt(1 + cent_num)) / 2
    # For integer or half-integer check:
    l_eff_frac = Fraction(cent_num, 4)
    # n_eff for ground state: n_eff = N_radial + l_eff + 1, N_radial = 0
    n_eff = l_eff_val + 1

    print(f"  {N_e:4d}  {d:5d}  S^{d-1:<3d}  SO({d})  {casimir_str:>15}  "
          f"{str(cent):>14}  {l_eff_val:6.2f}  {n_eff:8.2f}")

# Explicit table of free eigenvalues
print(f"\n  Free eigenvalue tables mu_free = nu(nu+d-2)/2:")
print(f"\n  {'nu':>4}", end="")
for N_e in [1, 2, 3, 4]:
    print(f"  {'N=' + str(N_e) + ' (S^' + str(3*N_e-1) + ')':>16}", end="")
print()
print("-" * 72)

for nu in range(0, 9):
    print(f"  {nu:4d}", end="")
    for N_e in [1, 2, 3, 4]:
        d = 3 * N_e
        mu = nu * (nu + d - 2) / 2
        print(f"  {mu:16.1f}", end="")
    print()


# =========================================================================
# STEP 2: Helium vs Lithium comparison table
# =========================================================================

print("\n" + "=" * 72)
print("STEP 2: He (2e) vs Li (3e) structural comparison")
print("=" * 72)

print(f"""
                                He (2 electrons)    Li (3 electrons)
                                ================    ================
Configuration space:            R^6                 R^9
Hyperradius:                    sqrt(r1^2+r2^2)     sqrt(r1^2+r2^2+r3^2)
Angular manifold:               S^5                 S^8
Symmetry group:                 SO(6)               SO(9)
Dimension of S^d:               5                   8

Casimir formula:                nu(nu+4)            nu(nu+7)
Free eigenvalues:               nu(nu+4)/2          nu(nu+7)/2
                                0, 6, 16, 30, 48    0, 9, 20, 33, 48

Centrifugal:                    15/(8R^2)           48/(8R^2) = 6/R^2
l_eff:                          3/2                 3
n_eff (ground state):           5/2                 4

Permutation group:              S_2                 S_3
Irreps:                         symmetric, anti     trivial, sign, standard
Ground state spatial:           symmetric (^1S)     mixed [2,1] (^2S)
Selection rule:                 Dnu = 0 mod 4       (complex, see below)

Nuclear singularities:          2 (r1=0, r2=0)     3 (r1=0, r2=0, r3=0)
e-e singularities:              1 (r12=0)           3 (r12=0, r13=0, r23=0)
Total singularities:            3                   6

Charge function:                -Z/r1 - Z/r2        -Z/r1 - Z/r2 - Z/r3
                                + 1/r12              + 1/r12 + 1/r13 + 1/r23

Ionization channels:            He+ + e-             Li+ + e-
                                                     Li2+ + 2e-
Thresholds:                     -Z^2/2 = -2.0 Ha    -7.280 Ha (Li+)
                                                     -4.500 Ha (Li2+)

Topology flow:                  S^5 -> S^3 x R      S^8 -> S^5 x R -> S^3 x R^2
Number of transitions:          1                    2
""")


# =========================================================================
# STEP 3: Selection rules from S_3 permutation symmetry
# =========================================================================

print("=" * 72)
print("STEP 3: Selection rules from S_3 permutation symmetry")
print("=" * 72)

# For He (S_2):
# Exchange: alpha -> pi/2 - alpha (swap r1, r2)
# Free eigenstates: sin(2*n*alpha), which transforms as (-1)^n under exchange
# Symmetric (singlet) ground state: n = 1, 3, 5, ... (odd n)
# Selection rule: <m|V|n> = 0 unless m+n even, i.e., Dnu = 0 mod 4

# For Li (S_3):
# S_3 has 6 elements: {e, (12), (13), (23), (123), (132)}
# Irreps: trivial [3], sign [1,1,1], standard [2,1] (2-dimensional)
#
# The ground state of Li (^2S) has 3 electrons: 1s^2 2s^1
# Total spin: S = 1/2 (doublet)
# Spin wavefunction: mixed symmetry [2,1] under S_3
# Spatial wavefunction: also mixed symmetry [2,1] (to give overall antisymmetric)
#
# Selection rules: matrix elements <psi'|V|psi> are nonzero only if
# the product irrep of psi' x V x psi contains the trivial irrep.
# V is symmetric under all permutations (Coulomb is symmetric)
# => V transforms as [3] (trivial)
# => <psi'|V|psi> nonzero iff psi' and psi are in the SAME irrep.

# Character table of S_3:
print(f"""
  S_3 Character Table:

  Irrep     |  e    (12)   (123)  |  dim
  ---------|----------------------|------
  [3]       |  1      1      1    |   1    (symmetric)
  [1,1,1]   |  1     -1      1    |   1    (antisymmetric)
  [2,1]     |  2      0     -1    |   2    (standard = Li ground state)

  Selection rules for Li (spatial [2,1]):

  1. V is symmetric => V ~ [3] (trivial irrep)
  2. [2,1] x [3] x [2,1] = [2,1] x [2,1] = [3] + [2,1] + [1,1,1]
     Contains [3] => matrix elements are generally NONZERO

  3. More specifically, the angular quantum numbers must satisfy:
     - Grand angular momentum nu in the [2,1] sector
     - Allowed nu values: need detailed angular analysis

  For the Casimir representation (nu, 0, 0, 0) of SO(9):
     - The [2,1] sector of S_3 selects specific nu values
     - In He: S_2 symmetric => even nu only (0, 2, 4, ...)
     - In Li: S_3 [2,1] => more complex restriction

  The [2,1] representation of S_3 acting on 3 particles selects
  angular momentum states with specific permutation symmetry.
  For 3 identical particles on S^8:

     nu = 0:  [3] only (fully symmetric)      => NOT in [2,1]
     nu = 1:  [2,1] appears (first mixed)     => ALLOWED
     nu = 2:  [3] + [2,1] + [1,1,1]           => ALLOWED
     nu = 3:  [3] + 2*[2,1] + [1,1,1]         => ALLOWED (double)

  So the LOWEST free eigenvalue for Li [2,1] is at nu = 1:
     mu_free = 1*(1+7)/2 = 4
     (Not nu = 0 like He, because [2,1] doesn't appear at nu = 0)
""")


# =========================================================================
# STEP 4: Asymptotic topology and multiple transitions
# =========================================================================

print("=" * 72)
print("STEP 4: Asymptotic topology -- double surgery")
print("=" * 72)

# Li has TWO ionization thresholds:
#
# Channel 1: Li -> Li+(1S, 2e) + e-
#   Threshold: E(Li+) = -7.2799 Ha (exact, Hylleraas)
#   Topology: S^8 -> S^5 x R (two-electron core + radial coordinate for 3rd electron)
#   At large R: mu(R) ~ E(Li+) * R^2
#
# Channel 2: Li -> Li2+(1s, 1e) + 2e-
#   Threshold: E(Li2+) = -Z^2/2 = -4.5 Ha
#   Topology: S^5 x R -> S^3 x R^2 (one-electron core + two radial coordinates)
#   Only relevant at VERY large R (both outer electrons far away)
#
# The Li ground state (-7.478 Ha) is below Li+ (-7.280 Ha),
# so the first ionization energy is:
# IE_1 = -7.280 - (-7.478) = 0.198 Ha = 5.39 eV

E_Li_exact = -7.478060  # Exact nonrelativistic Li ground state
E_Li_plus = -7.279913   # Exact He-like Li+ (Z=3)
E_Li_2plus = -4.5        # Exact H-like Li2+ (Z=3)

IE_1 = E_Li_plus - E_Li_exact
IE_2 = E_Li_2plus - E_Li_plus
IE_total = E_Li_2plus - E_Li_exact

print(f"\n  Li energy levels:")
print(f"    Li (3e):  E = {E_Li_exact:.6f} Ha")
print(f"    Li+ (2e): E = {E_Li_plus:.6f} Ha")
print(f"    Li2+ (1e):E = {E_Li_2plus:.6f} Ha")
print(f"    Li3+ (0e):E =  0.000000 Ha")

print(f"\n  Ionization energies:")
print(f"    IE_1 (Li -> Li+ + e-):    {IE_1:.6f} Ha = {IE_1*27.211:.2f} eV")
print(f"    IE_2 (Li+ -> Li2+ + e-):  {IE_2:.6f} Ha = {IE_2*27.211:.2f} eV")
print(f"    IE_total:                  {IE_total:.6f} Ha = {IE_total*27.211:.2f} eV")

print(f"""
  TOPOLOGY FLOW:

  R = 0:     Pure S^8 (round, all 3 electrons at nucleus)
             mu_free(nu=1) = 4 (lowest [2,1] eigenvalue)

  R ~ 1-3:   First transition region
             3rd electron begins to detach
             S^8 starts to "pinch" along one direction

  R ~ 5-10:  First surgery complete
             S^8 -> S^5 x R (Li+ core + outer electron)
             mu(R) ~ E(Li+) * R^2 = -7.280 * R^2

  R >> 10:   Second transition (double ionization channel)
             S^5 starts to pinch (2nd electron detaches from Li+ core)
             S^5 x R -> S^3 x R^2 (Li2+ core + 2 outer electrons)
             mu(R) ~ E(Li2+) * R^2 = -4.5 * R^2
             (Only relevant for doubly-excited states, not ground state)
""")


# =========================================================================
# STEP 5: Perturbation theory structure
# =========================================================================

print("=" * 72)
print("STEP 5: Perturbation theory -- a_1 for Li")
print("=" * 72)

# For He (2e):
# a_1(n=1) = (4/pi) * (-2Z * I_nuc(1) + sqrt(2) * I_ee_rat(1))
# where I_nuc(1) = 4/3, I_ee_rat(1) = 2/3
# a_1 = (4/pi) * (-2Z * 4/3 + sqrt(2) * 2/3) = (8/3pi) * (sqrt(2) - 4Z)

# For Li (3e), the charge function is:
# C = -Z(1/r1 + 1/r2 + 1/r3) + 1/r12 + 1/r13 + 1/r23
#
# In hyperspherical coordinates for 3 particles:
# r_i = R * Omega_i(angles)  where Omega_i are direction cosines on S^8
# 1/r_i = 1/(R * Omega_i)
# 1/r_ij = 1/R * 1/|Omega_i - Omega_j|
#
# The charge function per unit R:
# C/R = -Z * sum_i 1/Omega_i + sum_{i<j} 1/|Omega_ij|
#
# So a_1 = <psi_0| C |psi_0> where psi_0 is the free eigenstate at nu=1
# (the lowest [2,1] state on S^8).

# We can estimate a_1 from the virial theorem and known energies:
# E = -Z^2/2 * (correction_factor) for the ground state.
# The quasi-Coulomb model gives: E = a_2 - Z_eff^2 / (2*n_eff^2)
# where Z_eff = -a_1 and n_eff = l_eff + 1 = 4.

# From the exact energy:
# E = a_2 - a_1^2 / (2*16) = a_2 - a_1^2/32
# If a_2 is small compared to the Coulomb term:
# E ~ -a_1^2/32 => a_1 ~ -sqrt(-32*E) = -sqrt(32*7.478) = -15.47

Z_Li = 3
l_eff_Li = 3.0
n_eff_Li = l_eff_Li + 1  # = 4

a1_Li_estimate = -np.sqrt(-2 * n_eff_Li**2 * E_Li_exact)
print(f"\n  Estimated a_1(Li) from E = -a_1^2/(2*n_eff^2):")
print(f"    n_eff = l_eff + 1 = {n_eff_Li:.1f}")
print(f"    E_exact = {E_Li_exact:.6f} Ha")
print(f"    a_1 estimate = -sqrt(2*n_eff^2*|E|) = {a1_Li_estimate:.4f}")
print(f"    Z_eff estimate = -a_1 = {-a1_Li_estimate:.4f}")

# For comparison, the He values:
from fractions import Fraction

def I_nuc_exact(n: int) -> Fraction:
    return sum(Fraction(1, 2*j+1) for j in range(2*n))

def I_ee_rational_exact(n: int) -> Fraction:
    s1 = sum(Fraction((-1)**k, 4*k+1) for k in range(n))
    s2 = sum(Fraction((-1)**k, 4*k+3) for k in range(n))
    return s1 - s2

def a1_exact_He(n: int, Z: float) -> float:
    I_nuc = float(I_nuc_exact(n))
    I_ee_rat = float(I_ee_rational_exact(n))
    return (4 / np.pi) * (-2 * Z * I_nuc + np.sqrt(2) * I_ee_rat)

a1_He = a1_exact_He(1, 2.0)
Z_eff_He = -a1_He

print(f"\n  Comparison:")
print(f"    He: a_1 = {a1_He:.6f}, Z_eff = {Z_eff_He:.6f}, n_eff = 2.5")
print(f"    Li: a_1 ~ {a1_Li_estimate:.4f}, Z_eff ~ {-a1_Li_estimate:.4f}, n_eff = 4.0")
print(f"\n    He: E_qc = -{Z_eff_He**2}/(2*6.25) = {-Z_eff_He**2/12.5:.4f}")
print(f"    Li: E_qc ~ -{(-a1_Li_estimate)**2}/(2*16) = {-a1_Li_estimate**2/32:.4f}")

# How does a_1 scale with Z?
# For He: a_1 = (8/3pi)(sqrt(2) - 4Z) ~ -32Z/(3pi) for large Z
# For Li: expect a_1 ~ -c * Z with c involving pi and sqrt-type factors
# The nuclear part should dominate: 3 nuclear terms (vs 2 for He)
# The e-e part: 3 repulsion terms (vs 1 for He)

print(f"""
  Expected a_1 structure for Li:

  a_1(Li, Z) = (normalization) * (-Z * I_nuc^(3e) + I_ee^(3e))

  where:
    I_nuc^(3e) = integral of (1/Omega_1 + 1/Omega_2 + 1/Omega_3)
                 over the [2,1] ground state on S^8

    I_ee^(3e) = integral of (1/|Omega_12| + 1/|Omega_13| + 1/|Omega_23|)
                over the [2,1] ground state on S^8

  By the 3-fold symmetry of the [2,1] representation:
    I_nuc^(3e) = 3 * I_nuc^(single)(S^8)
    I_ee^(3e)  = 3 * I_ee^(single)(S^8)

  So a_1(Li, Z) = 3 * (normalization') * (-Z * I_nuc' + I_ee')

  This factor of 3 (number of particles) is analogous to the factor of 2
  that appears implicitly in the He formula through the (1/cos + 1/sin)
  symmetry.
""")


# =========================================================================
# STEP 6: Comparison with existing FCI results
# =========================================================================

print("=" * 72)
print("STEP 6: FCI results as validation targets")
print("=" * 72)

# From the codebase (LatticeIndex with slater_full + exact h1):
# Li nmax=3: E ~ -7.340 Ha (some percent error)
# Li nmax=4: E ~ -7.346 Ha (improving)
# Li nmax=5: E ~ -7.395 Ha (1.10% error vs exact -7.478)
# Exact: -7.478060 Ha

print(f"""
  Existing FCI results for Li (from DirectCISolver):

  {'nmax':>5}  {'n_spatial':>10}  {'n_sd':>10}  {'E (Ha)':>12}  {'Error':>8}  {'Time':>8}
  {'-'*60}
  {'3':>5}  {'14':>10}  {'3,276':>10}  {'-7.340':>12}  {'1.85%':>8}  {'0.4s':>8}
  {'4':>5}  {'30':>10}  {'34,220':>10}  {'-7.346':>12}  {'1.77%':>8}  {'7.6s':>8}
  {'5':>5}  {'55':>10}  {'215,820':>10}  {'-7.395':>12}  {'1.10%':>8}  {'116s':>8}

  Exact (nonrelativistic): -7.478060 Ha

  For comparison, the hyperspherical He solver achieves 0.05% error
  in 3-4 seconds. The FCI approach for Li at nmax=5 achieves 1.10%
  in 116 seconds. A hyperspherical Li solver COULD potentially
  achieve similar or better accuracy much faster.
""")


# =========================================================================
# STEP 7: Computational cost estimate for 3e hyperspherical
# =========================================================================

print("=" * 72)
print("STEP 7: Computational cost estimate")
print("=" * 72)

# For He (2e):
# Angular problem: 1D (alpha) x partial waves (l)
# Grid: n_alpha = 100-200, l_max = 0-3
# Matrix: N = (l_max+1) * n_alpha = 400 at most
# Diagonalization: O(N^3) = O(64M) ~ 50 ms
# Adiabatic curve: 100-200 R-points * 50ms = 5-10 s

# For Li (3e):
# Angular problem: 2D (alpha_1, alpha_2) x partial waves
# Grid: n_alpha1 = 50-100, n_alpha2 = 50-100
# Partial waves: l_max per electron pair
# Matrix dimension estimation:
#   - 2D alpha grid: n_a1 * n_a2 = 50*50 = 2500 points
#   - Partial waves: need l12, l13, l23 quantum numbers
#   - For L=0 total: constrained, but still many combinations
#   - Estimate: 3 pair indices * l_max^2 ~ 3*4 = 12 channels
#   - Total: 2500 * 12 = 30,000 matrix dimension (DENSE!)
#   - Diagonalization: O(30000^3) = INFEASIBLE for dense

# Better approach: use symmetry reduction
# The [2,1] sector of S_3 reduces the 2D problem to a 1D "fundamental domain"
# After symmetry reduction: ~1/6 of the full grid (S_3 has 6 elements)
# So N ~ 5000 is more realistic

# Alternatively: use adiabatic hyperspherical with Jacobi coordinates
# (rho, theta) + molecular-frame angles
# This reduces to a sequence of 1D problems (like He)
# Cost: similar to He but with more channels

print(f"""
  Computational cost comparison:

  {'':>30}  {'He (2e)':>15}  {'Li (3e, naive)':>15}  {'Li (3e, symm.)':>15}
  {'':>30}  {'='*15}  {'='*15}  {'='*15}
  Angular coordinates:          1D (alpha)       2D (a1, a2)       1D (reduced)
  Grid per coordinate:          100-200          50-100            50-100
  Total grid points:            100-200          2,500-10,000      500-1,000
  Partial wave channels:        4 (l=0,1,2,3)    ~12-50           ~6-12
  Angular matrix N:             400              30,000-500,000    3,000-12,000
  Dense diag (O(N^3)):          50 ms            hours-days        10-300 s
  Sparse diag (eigsh):          N/A              10-100 s          1-10 s
  Adiabatic curve (100 pts):    5 s              1,000-10,000 s    100-1,000 s
  Radial solve:                 0.1 s            0.1 s             0.1 s
  TOTAL:                        ~5 s             ~hours            ~10 min

  KEY INSIGHT: With symmetry reduction, a 3e hyperspherical solver
  is feasible but ~100x more expensive than the 2e He solver.
  This is still MUCH cheaper than FCI at nmax=5 (116 s for 215k SDs).

  The main challenge is the IMPLEMENTATION, not the computation:
  - Need to derive the 3-particle angular Hamiltonian in reduced coordinates
  - Need the [2,1] projection operator for S_3
  - Need the Gaunt integrals generalized to 3-body multipoles
  - All of this requires a NEW PAPER (Paper 16 or similar)
""")


# =========================================================================
# STEP 8: The quasi-Coulomb model for Li
# =========================================================================

print("=" * 72)
print("STEP 8: Quasi-Coulomb model for Li")
print("=" * 72)

# By analogy with He:
# V_eff(R) = mu(R)/R^2 + centrifugal/R^2
# ~ a_1/R + a_2 + l_eff*(l_eff+1)/(2R^2)  at small R
#
# For He: l_eff = 3/2, n_eff = 5/2
# E = a_2 - Z_eff^2 / (2*(5/2)^2) = a_2 - Z_eff^2/12.5
#
# For Li: l_eff = 3, n_eff = 4
# E = a_2 - Z_eff^2 / (2*16) = a_2 - Z_eff^2/32

# From E_exact and assuming a_2 ~ 0:
Z_eff_Li = np.sqrt(-2 * n_eff_Li**2 * E_Li_exact)

# More refined: estimate a_2 from the He scaling
# For He: a_2/E_qc ~ 0.244/2.744 ~ 8.9%
# For Li: a_2 ~ 8.9% * E_exact ~ 0.67 Ha
# This is very rough!

# Instead, use the EXACT energies of Li and Li+ to constrain:
# E_Li = a_2 - Z_eff^2/32  ... (1)
# And from the asymptotic: mu(R->inf) ~ E_Li+ * R^2
# E_Li+ = -7.280 Ha is the "asymptotic constant" (like -Z^2/2 for He)

a2_Li_needed = E_Li_exact + Z_eff_Li**2 / 32
E_qc_Li = -Z_eff_Li**2 / 32

print(f"\n  Quasi-Coulomb model for Li:")
print(f"    l_eff = 3 (from centrifugal 48/(8R^2) = 6/R^2)")
print(f"    n_eff = l_eff + 1 = 4")
print(f"    Z_eff = sqrt(2*n_eff^2*|E|) = {Z_eff_Li:.6f}")
print(f"    E_qc(a_2=0) = -Z_eff^2/32 = {E_qc_Li:.6f} Ha")
print(f"    a_2 needed for exact E: {a2_Li_needed:.6f} Ha")
print(f"    a_2 / |E_exact|: {abs(a2_Li_needed/E_Li_exact)*100:.1f}%")

# For He: a_2 is 8.9% of |E_exact|
# For Li: a_2 would be ~4.4% of |E_exact| (smaller fraction)

# The "surgery correction" for Li:
# Delta_Li = E_exact - E_qc_Li(with a_2)
# Without knowing a_2, we can't compute this exactly.
# But we can estimate the structure:

# The E_qc model (no a_2, just Coulomb): E = -Z_eff^2/32
# This captures the DOMINANT binding.
# The a_2 correction adds the "angular correlation" at small R.
# The "surgery correction" is the rest.

print(f"\n  Energy hierarchy for Li:")
print(f"    Coulomb part (Z_eff^2/32):  {E_qc_Li:.4f} Ha")
print(f"    Exact:                       {E_Li_exact:.4f} Ha")
print(f"    Correlation energy:          {E_Li_exact - E_Li_plus:.4f} Ha")
print(f"    (Li exact - Li+ threshold)")

# Compare the correlations:
E_corr_He = -2.903724 - (-2.0)  # = -0.904 Ha
E_corr_Li = E_Li_exact - E_Li_plus  # = -0.198 Ha

print(f"\n  Correlation energy comparison:")
print(f"    He: E_corr = E(He) - E(He+) = {E_corr_He:.4f} Ha")
print(f"    Li: E_corr = E(Li) - E(Li+) = {E_corr_Li:.4f} Ha")
print(f"    Ratio: He/Li = {E_corr_He/E_corr_Li:.2f}")
print(f"    (He has 4.6x more correlation per electron pair)")


# =========================================================================
# STEP 9: Predictions for Li hyperspherical coefficients
# =========================================================================

print("\n" + "=" * 72)
print("STEP 9: Predictions for Li hyperspherical coefficients")
print("=" * 72)

# Based on the He analysis, we can make structural predictions for Li:

# 1. a_1 should have the form:
#    a_1(Li) = (normalization on S^8) * (-Z * I_nuc^(3e) + I_ee^(3e))
#    where I_nuc and I_ee are integrals over S^8 with [2,1] symmetry

# 2. The normalization factor comes from the ground state volume on S^8:
#    Vol(S^8) = 32*pi^4/105  (exact)
#    For He: Vol(S^5) = pi^3, normalization = 4/pi
#    The pattern: normalization ~ 1/pi^{N-1} where N = number of electrons

Vol_S5 = np.pi**3
Vol_S8 = 32 * np.pi**4 / 105

print(f"\n  Volumes of angular manifolds:")
print(f"    Vol(S^5) = pi^3 = {Vol_S5:.6f}")
print(f"    Vol(S^8) = 32*pi^4/105 = {Vol_S8:.6f}")
print(f"    Ratio: Vol(S^8)/Vol(S^5) = {Vol_S8/Vol_S5:.6f}")
print(f"    32*pi/105 = {32*np.pi/105:.6f}")

# 3. The selection rules should restrict which SO(9) representations couple.
#    For He: Dnu = 0 mod 4 (from S_2 symmetry in the [symmetric] sector)
#    For Li: the [2,1] sector of S_3 should give a different selection rule.
#    Since S_3 has a 3-fold symmetry: expect Dnu = 0 mod 3 or similar.

# 4. The perturbation series should converge faster or slower:
#    He a_2: 2 terms capture 96%, 30 terms give 99.5%
#    Li a_2: expect SLOWER convergence (more coupling channels, more singularities)

# 5. The quasi-Coulomb ground state:
#    He: E_qc = -2.744 (5.5% error)
#    Li: E_qc = ? (need a_1 and a_2)

# Estimate using the He scaling:
# For He: a_1 = -5.590, Z_eff = 5.590, error 5.5%
# Scale by Z^2: Z_eff(Li) ~ Z_eff(He) * (Z_Li/Z_He)^{something}
# But the N-electron scaling is different.

# Better: use the empirical Z_eff from the exact energy:
print(f"\n  Predicted Li hyperspherical structure:")
print(f"    Z_eff ~ {Z_eff_Li:.4f} (from exact energy, n_eff=4)")
print(f"    a_1 ~ {a1_Li_estimate:.4f} (if a_2 ~ 0)")
print(f"    Quasi-Coulomb error: likely 5-10% (similar to He)")
print(f"    Selection rule: Dnu = 0 mod ? (need S_3 analysis)")
print(f"    Convergence of a_2: likely slower than He")


# =========================================================================
# STEP 10: Topology flow comparison
# =========================================================================

print("\n" + "=" * 72)
print("STEP 10: Topology flow comparison He vs Li")
print("=" * 72)

print(f"""
  He FLOW (S^5 -> S^3 x R):
  --------------------------
  R ~ 0:     Round S^5. mu_free = 0. Both electrons near nucleus.
  R ~ 1:     Perturbative. mu ~ a_1*R. Eigenstate ~ sin(2*alpha).
  R ~ 3:     Crossover. eta ~ 0.88 (non-conformal). IPR peaks.
  R ~ 5+:    First surgery. psi localizes at alpha ~ 0.
             One electron near nucleus (He+), one far away.
  R >> 10:   Complete. mu ~ -Z^2*R^2/2. S^3 x R topology.

  Transitions: 1
  Surgery cost: Delta ~ 0.16 Ha (Z=2)

  Li FLOW (S^8 -> S^5 x R -> S^3 x R^2):
  -----------------------------------------
  R ~ 0:     Round S^8. mu_free = 4 (lowest [2,1] eigenvalue).
             All 3 electrons near nucleus.
  R ~ 1:     Perturbative. mu ~ 4 + a_1*R.
             Shell structure begins: 1s^2 core + 2s valence.
  R ~ 3-5:   FIRST crossover. 2s electron begins to detach.
             S^8 pinches along the valence direction.
             This is the IONIZATION transition.
  R ~ 10+:   First surgery complete. S^8 -> S^5 x R.
             mu ~ E(Li+) * R^2 = -7.280 * R^2.
             Two-electron Li+ core (S^5) + radial coordinate.
  R >> 20:   SECOND surgery possible (only for excited states).
             S^5 x R -> S^3 x R^2.
             But Li ground state never reaches this regime.

  Transitions: 1 (for ground state) or 2 (for doubly-excited states)
  Surgery cost: Delta_Li ~ ? (need hyperspherical solver)

  KEY DIFFERENCE:
  In He, the two electrons are equivalent (S_2 symmetry).
  The surgery is "democratic" -- either electron can detach.

  In Li, the 2s electron is distinguishable from the 1s^2 core.
  The surgery is "hierarchical" -- the valence electron detaches FIRST.
  This breaks the S_3 symmetry to S_2 x {{identity}} at the transition.

  The hierarchy is built into the [2,1] representation:
  [2,1] = one pair symmetric + one particle separate.
  This IS the 1s^2 + 2s shell structure, encoded group-theoretically.
""")


# =========================================================================
# FINAL SUMMARY
# =========================================================================

print("=" * 72)
print("FINAL SUMMARY")
print("=" * 72)

print(f"""
ANSWERS TO KEY QUESTIONS:

1. Does Li follow mu_free = lambda(lambda+7)/2 (SO(9) Casimir)?

   YES (by construction). For N electrons in 3D, the free angular
   eigenvalues on S^(3N-1) are the SO(3N) Casimir values.
   Li (N=3): mu_free = nu(nu+7)/2, giving 0, 9/2, 10, 33/2, ...

   HOWEVER: the Li ground state is NOT at nu=0 (unlike He).
   The [2,1] representation of S_3 first appears at nu=1,
   so the LOWEST allowed eigenvalue is mu_free = 1*(1+7)/2 = 4.

2. What are the selection rules from S_3 permutation symmetry?

   The charge function V is totally symmetric ([3] irrep of S_3).
   Matrix elements <psi'|V|psi> in the [2,1] sector are generally
   nonzero (since [2,1] x [3] x [2,1] contains [3]).

   The specific selection rule on nu depends on the angular
   coordinate parametrization and requires detailed computation.
   Expected: more permissive than He's Dnu = 0 mod 4.

3. What is the asymptotic topology? S^5 x R or S^3 x R^2?

   For the ground state: S^8 -> S^5 x R (single surgery).
   The Li+ (2e) core lives on S^5, outer electron on R.

   Double surgery (S^8 -> S^3 x R^2) only occurs for
   doubly-excited states where both outer electrons detach.

4. Is there a single crossover region or multiple transitions?

   For the ground state: SINGLE transition (analogous to He).
   For excited states: potentially TWO transitions (double ionization).

   The first transition (2s detachment) occurs at R_c ~ 3-5 bohr.
   The second transition (1s detachment from Li+) would occur at
   R_c ~ 10+ bohr but is not relevant for the ground state.

ASSESSMENT: WHAT'S NEEDED FOR A Li HYPERSPHERICAL SOLVER

  1. THEORY: Derive the 3-particle angular Hamiltonian in Jacobi
     hyperspherical coordinates with S_3 symmetry projection.
     This is a NEW PAPER (not in current papers/core/).

  2. CODE: Implement the angular solver on S^8 (or symmetry-reduced).
     Estimated: 3,000-12,000 matrix dimension (vs 400 for He).
     Feasible but requires sparse eigenvalue methods.

  3. VALIDATION: Compare against existing FCI results (1.10% at nmax=5).
     The hyperspherical approach should achieve < 0.1% (like He) with
     sufficient angular resolution.

  4. COST: ~100x more expensive per angular solve than He, but
     overall still cheaper than FCI at nmax=5 for high accuracy.

  STATUS: Theoretically well-defined, computationally feasible,
  but requires significant implementation effort (new paper + code).
""")
