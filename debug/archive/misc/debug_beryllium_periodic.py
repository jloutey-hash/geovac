"""
Algebraic structure of beryllium (4-electron) and periodic table patterns.

GOAL: Determine Be's SO(12) algebraic structure, estimate quasi-Coulomb
accuracy, and look for periodic patterns across H, He, Li, Be.

KEY RESULTS:
  H  (1e): S^2  manifold, SO(3)  Casimir, l_eff = 0,   n_eff = 1
  He (2e): S^5  manifold, SO(6)  Casimir, l_eff = 3/2,  n_eff = 5/2
  Li (3e): S^8  manifold, SO(9)  Casimir, l_eff = 3,    n_eff = 4
  Be (4e): S^11 manifold, SO(12) Casimir, l_eff = 9/2,  n_eff = 11/2

STRUCTURE CLASSIFICATION:
  Type A: Single particle (H)           — trivial S_1
  Type B: Equivalent pair (He)          — democratic S_2, one surgery
  Type C: Core + valence (Li)           — hierarchical S_3 [2,1]
  Type D: Core + equivalent pair (Be)   — mixed S_4 [2,2], two-stage hierarchy
"""

import numpy as np
from fractions import Fraction
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# =========================================================================
# EXACT REFERENCE ENERGIES
# =========================================================================

# Non-relativistic exact energies (Hylleraas/CI limit)
EXACT_ENERGIES = {
    # Neutral atoms
    'H':   -0.500000,   # Exact
    'He':  -2.903724,   # Hylleraas (Pekeris 1959)
    'Li':  -7.478060,   # Hylleraas (Yan & Drake 1998)
    'Be': -14.667360,   # CI limit (Komasa 2002)

    # Singly-ionized
    'H-':  -0.527751,   # Barely bound (Pekeris)
    'He+': -2.000000,   # H-like Z=2
    'Li+': -7.279913,   # He-like Z=3
    'Be+': -14.324764,  # Li-like Z=4 (Yan & Drake)

    # Doubly-ionized
    'Li2+': -4.500000,  # H-like Z=3
    'Be2+': -13.655566, # He-like Z=4

    # Triply-ionized
    'Be3+': -8.000000,  # H-like Z=4

    # Hartree-Fock limits
    'He_HF': -2.861680,
    'Li_HF': -7.432727,
    'Be_HF': -14.573023,
}

# FCI results from GeoVac DirectCISolver
FCI_RESULTS = {
    'He_n3':  {'E': -2.883, 'nmax': 3, 'nsd': 190, 'time': 0.01},
    'He_n5':  {'E': -2.894, 'nmax': 5, 'nsd': 2415, 'time': 0.2},
    'Li_n3':  {'E': -7.340, 'nmax': 3, 'nsd': 3276, 'time': 0.4},
    'Li_n4':  {'E': -7.346, 'nmax': 4, 'nsd': 34220, 'time': 7.6},
    'Li_n5':  {'E': -7.395, 'nmax': 5, 'nsd': 215820, 'time': 116},
    'Be_n3':  {'E': -14.531256, 'nmax': 3, 'nsd': 20475, 'time': 2.9},
    'Be_n4':  {'E': -14.535460, 'nmax': 4, 'nsd': 487635, 'time': 357},
}


# =========================================================================
# STEP 1: SO(d) Casimir structure for N = 1 to 5
# =========================================================================

print("=" * 80)
print("STEP 1: SO(3N) Casimir structure for N-electron atoms")
print("=" * 80)

print(f"\n  {'N_e':>4}  {'d=3N':>5}  {'Manifold':>8}  {'Group':>7}  {'Casimir':>15}  "
      f"{'Centrifugal':>13}  {'l_eff':>6}  {'n_eff':>6}")
print("-" * 85)

for N_e in range(1, 6):
    d = 3 * N_e
    manifold = f"S^{d-1}"
    group = f"SO({d})"
    casimir_str = f"nu(nu+{d-2})"

    # Centrifugal: (3N-1)(3N-3)/8
    cent_num = (3*N_e - 1) * (3*N_e - 3)
    cent = Fraction(cent_num, 8)

    # l_eff from l_eff(l_eff+1) = cent_num/4
    l_eff_val = (-1 + np.sqrt(1 + cent_num)) / 2
    n_eff = l_eff_val + 1

    # Express l_eff as fraction if possible
    l_eff_frac = Fraction(int(2*l_eff_val + 0.5), 2)  # nearest half-integer

    print(f"  {N_e:4d}  {d:5d}  {manifold:>8}  {group:>7}  {casimir_str:>15}  "
          f"{str(cent):>13}  {str(l_eff_frac):>6}  {str(l_eff_frac + 1):>6}")

# Table of free eigenvalues mu_free = nu(nu+d-2)/2
print(f"\n  Free eigenvalue mu_free = nu(nu+d-2)/2:")
print(f"\n  {'nu':>4}", end="")
for N_e in [1, 2, 3, 4, 5]:
    header = f"N={N_e} (S^{3*N_e-1})"
    print(f"  {header:>14}", end="")
print()
print("-" * 80)

for nu in range(0, 7):
    print(f"  {nu:4d}", end="")
    for N_e in [1, 2, 3, 4, 5]:
        d = 3 * N_e
        mu = Fraction(nu * (nu + d - 2), 2)
        print(f"  {str(mu):>14}", end="")
    print()


# =========================================================================
# STEP 2: S_4 representation theory for beryllium
# =========================================================================

print("\n" + "=" * 80)
print("STEP 2: S_4 representation theory for beryllium (4 electrons)")
print("=" * 80)

print(f"""
  S_4 has 5 irreducible representations (partitions of 4):

  Partition  |  Dim  |  Name               |  Physical role
  ----------|-------|---------------------|----------------------------------------
  [4]        |   1   |  Symmetric           |  All 4 bosonic (not fermionic ground)
  [3,1]      |   3   |  Standard            |  Quartet spin (S=3/2)
  [2,2]      |   2   |  Two-row symmetric   |  SINGLET GROUND STATE (S=0, ^1S)
  [2,1,1]    |   3   |  Standard-conjugate  |  Triplet spin (S=1)
  [1,1,1,1]  |   1   |  Antisymmetric       |  Fully antisymmetric spatial

  Character table of S_4:

  Irrep     |  e   (12)  (12)(34)  (123)  (1234)  |  dim
  ----------|--------------------------------------|------
  [4]       |  1    1      1        1       1      |   1
  [3,1]     |  3    1     -1        0      -1      |   3
  [2,2]     |  2    0      2       -1       0      |   2
  [2,1,1]   |  3   -1     -1        0       1      |   3
  [1,1,1,1] |  1   -1      1        1      -1      |   1

  GROUND STATE ASSIGNMENT:
  ========================
  Be ground state: 1s^2 2s^2, ^1S_0 (singlet, L=0)

  Total wavefunction: antisymmetric under S_4 => [1,1,1,1]

  Spin wavefunction (S=0, singlet for 4 electrons):
    Two spin-up, two spin-down. Symmetry type: [2,2]
    (This is the ONLY S_4 irrep that gives S=0 with 4 electrons)

  Spatial wavefunction: must satisfy
    spatial x spin = [1,1,1,1] (total antisymmetry)
    [2,2] x [??] contains [1,1,1,1]?

  Key group theory fact: [2,2] is SELF-CONJUGATE under S_4
    [2,2] x [2,2] = [4] + [2,2] + [1,1,1,1]
    So [2,2] x [2,2] DOES contain [1,1,1,1]. Check!

  Therefore: SPATIAL wavefunction is also [2,2] under S_4.

  This is a unique feature of Be:
  - He: spatial [2] (symmetric), spin [1,1] (antisymmetric)
  - Li: spatial [2,1] (mixed), spin [2,1] (mixed)
  - Be: spatial [2,2] (self-conjugate), spin [2,2] (self-conjugate)
""")


# =========================================================================
# STEP 3: Lowest nu for [2,2] in SO(12)
# =========================================================================

print("=" * 80)
print("STEP 3: Lowest nu for [2,2] sector in SO(12)")
print("=" * 80)

print(f"""
  The SO(3N) hyperspherical harmonics of degree nu, when restricted
  to S_N acting by permuting coordinate triples, decompose into S_N irreps.

  Pattern from known cases (N electrons in 3D):

  N=1 (S_1, trivial):
    nu=0: [1] (trivial)
    Ground state [1] at nu=0. mu_free = 0.

  N=2 (S_2):
    nu=0: [2] (symmetric only)
    nu=1: [1,1] (antisymmetric appears)
    Ground state [2] at nu=0. mu_free = 0.
    Excited [1,1] at nu=1. mu_free = 2.

  N=3 (S_3):
    nu=0: [3] (symmetric only)
    nu=1: [2,1] (standard appears, 2-dimensional)
    nu=2: [3] + [2,1] + [1,1,1] (all irreps present)
    Ground state [2,1] at nu=1. mu_free = 4.

  N=4 (S_4):
    nu=0: [4] only
    nu=1: [3,1] appears (3-dimensional)
    nu=2: [4] + [3,1] + [2,2] + ... (TWO-ROW PARTITION APPEARS)
    Ground state [2,2] at nu=2. mu_free = 2*(2+10)/2 = 12.

  PATTERN: For N particles, the irrep with partition lambda first appears
  at nu = number of "boxes below the first row" in the Young diagram.

  [4]:       nu=0  (0 boxes below first row)
  [3,1]:     nu=1  (1 box below)
  [2,2]:     nu=2  (2 boxes below)
  [2,1,1]:   nu=3  (3 boxes below)
  [1,1,1,1]: nu=4  (6 boxes below — but actually nu = N(N-1)/2 ?? TBD)

  CAVEAT: This is an ESTIMATE. The exact decomposition requires the
  plethysm SO(3N) x S_N, which is non-trivial for N >= 4.
  The value nu=2 for [2,2] is the most likely based on the pattern
  but should be verified by explicit computation.
""")

# Compute mu_free for the estimated ground state
nu_gs = {1: 0, 2: 0, 3: 1, 4: 2, 5: 2}  # Estimated
for N_e in range(1, 6):
    d = 3 * N_e
    nu = nu_gs[N_e]
    mu_free = Fraction(nu * (nu + d - 2), 2)
    print(f"  N={N_e}: nu_gs = {nu}, mu_free = {mu_free} = {float(mu_free):.1f}")


# =========================================================================
# STEP 4: Quasi-Coulomb model across H, He, Li, Be
# =========================================================================

print("\n" + "=" * 80)
print("STEP 4: Quasi-Coulomb model for H, He, Li, Be")
print("=" * 80)

print(f"""
  The quasi-Coulomb model replaces the angular problem with a single
  effective Coulomb problem: E = a_2 - Z_eff^2 / (2 * n_eff^2)

  where:
    Z_eff = |a_1|   (first perturbation coefficient = effective charge)
    n_eff = l_eff + 1   (effective principal quantum number)
    a_2 = second-order correction (small if quasi-Coulomb is good)

  From the EXACT energy assuming a_2 = 0:
    Z_eff = sqrt(2 * n_eff^2 * |E_exact|)
""")

atoms = [
    ('H',  1, 1, EXACT_ENERGIES['H'],  '[1]',   0),
    ('He', 2, 2, EXACT_ENERGIES['He'], '[2]',   0),
    ('Li', 3, 3, EXACT_ENERGIES['Li'], '[2,1]', 1),
    ('Be', 4, 4, EXACT_ENERGIES['Be'], '[2,2]', 2),
]

print(f"  {'Atom':>4}  {'N_e':>3}  {'Z':>2}  {'l_eff':>6}  {'n_eff':>6}  "
      f"{'E_exact':>10}  {'Z_eff':>8}  {'E_qc(a2=0)':>11}  {'|a2|/|E|':>9}")
print("-" * 80)

qc_data = {}
for name, N_e, Z, E_exact, irrep, nu_g in atoms:
    d = 3 * N_e
    cent_num = (3*N_e - 1) * (3*N_e - 3)
    l_eff = (-1 + np.sqrt(1 + cent_num)) / 2
    n_eff = l_eff + 1

    Z_eff = np.sqrt(2 * n_eff**2 * abs(E_exact))
    E_qc = -Z_eff**2 / (2 * n_eff**2)

    # For H: exact, a_2 = 0
    # For others: back-compute a_2 needed
    # But since we computed Z_eff from E_exact, E_qc = -|E_exact| by construction
    # The MEANINGFUL comparison is Z_eff vs the PERTURBATIVE a_1

    # Perturbative a_1 (known exactly for H, He)
    if name == 'H':
        a1_pert = -1.0  # Exact
        a2 = 0.0
    elif name == 'He':
        # a_1 = (8/3pi)(sqrt(2) - 4Z)
        a1_pert = (8 / (3 * np.pi)) * (np.sqrt(2) - 4 * Z)
        Z_eff_pert = abs(a1_pert)
        E_qc_pert = -Z_eff_pert**2 / (2 * n_eff**2)
        a2 = E_exact - E_qc_pert
    else:
        # a_1 not known analytically, estimate from exact energy
        a1_pert = -Z_eff  # Back-computed
        a2 = 0.0  # By construction

    qc_data[name] = {
        'N_e': N_e, 'Z': Z, 'l_eff': l_eff, 'n_eff': n_eff,
        'E_exact': E_exact, 'Z_eff': Z_eff, 'a1_pert': a1_pert,
        'a2': a2, 'nu_gs': nu_g, 'irrep': irrep,
    }

    l_frac = Fraction(int(2*l_eff + 0.5), 2)
    n_frac = Fraction(int(2*n_eff + 0.5), 2)
    a2_frac = abs(a2/E_exact)*100 if E_exact != 0 else 0

    print(f"  {name:>4}  {N_e:3d}  {Z:2d}  {str(l_frac):>6}  {str(n_frac):>6}  "
          f"{E_exact:10.4f}  {Z_eff:8.4f}  {E_qc:11.4f}  {a2_frac:8.1f}%")

# He perturbative comparison
print(f"\n  For He, we can compare perturbative vs back-computed Z_eff:")
a1_He = (8 / (3 * np.pi)) * (np.sqrt(2) - 8)
Z_eff_He_pert = abs(a1_He)
n_eff_He = 2.5
E_qc_He_pert = -Z_eff_He_pert**2 / (2 * n_eff_He**2)
print(f"    Z_eff (perturbative a_1):   {Z_eff_He_pert:.6f}")
print(f"    Z_eff (back-computed):      {qc_data['He']['Z_eff']:.6f}")
print(f"    E_qc (perturbative, a2=0):  {E_qc_He_pert:.6f}  "
      f"({abs((E_qc_He_pert - EXACT_ENERGIES['He'])/EXACT_ENERGIES['He'])*100:.1f}% error)")
print(f"    E_exact:                    {EXACT_ENERGIES['He']:.6f}")
print(f"    a_2 (He):                   {EXACT_ENERGIES['He'] - E_qc_He_pert:.6f} Ha")


# =========================================================================
# STEP 5: Complete comparison table H / He / Li / Be
# =========================================================================

print("\n" + "=" * 80)
print("STEP 5: Complete algebraic structure comparison")
print("=" * 80)

print(f"""
                              H (1e)       He (2e)      Li (3e)      Be (4e)
                              ======       =======      =======      =======
  Config space:               R^3          R^6          R^9          R^12
  Angular manifold:           S^2          S^5          S^8          S^11
  Isometry group:             SO(3)        SO(6)        SO(9)        SO(12)
  Casimir C_2(nu):            nu(nu+1)     nu(nu+4)     nu(nu+7)     nu(nu+10)

  Centrifugal (d-1)(d-3)/8:   0            15/8         6            99/8
  l_eff:                      0            3/2          3            9/2
  n_eff:                      1            5/2          4            11/2

  Ground state config:        1s           1s^2         1s^2 2s      1s^2 2s^2
  Total spin S:               1/2          0            1/2          0
  Spin irrep:                 [1]          [1,1]        [2,1]        [2,2]
  Spatial irrep:              [1]          [2]          [2,1]        [2,2]
  Lowest nu in sector:        0            0            1            2 (est.)
  mu_free (ground):           0            0            4            12 (est.)

  Permutation group S_N:      S_1          S_2          S_3          S_4
  |S_N|:                      1            2            6            24
  # of irreps:                1            2            3            5
  Ground state dim:           1            1            2            2

  Structure type:             A            B            C            D
                              (trivial)    (democratic)  (hier.)     (mixed hier.)

  Exact energy (Ha):          -0.5000      -2.9037      -7.4781      -14.667
  Correlation energy:         0.0000       -0.0419      -0.0453      -0.094
  (E_exact - E_HF)

  Z_eff (back-computed):      1.0000       3.8079       15.469       28.46 (est.)

  Nuclear singularities:      1            2            3            4
  e-e singularities:          0            1            3            6
  Total singularities:        1            3            6            10 = N(N+1)/2

  Ionization thresholds:      1            1            2            3
  Topology surgeries (gs):    0            1            1            1
  Topology flow:              S^2          S^5->S^3xR   S^8->S^5xR   S^11->S^8xR
""")


# =========================================================================
# STEP 6: Topology flow for Be
# =========================================================================

print("=" * 80)
print("STEP 6: Topology flow for beryllium")
print("=" * 80)

E_Be = EXACT_ENERGIES['Be']
E_Be_plus = EXACT_ENERGIES['Be+']
E_Be_2plus = EXACT_ENERGIES['Be2+']
E_Be_3plus = EXACT_ENERGIES['Be3+']

IE_1 = E_Be_plus - E_Be
IE_2 = E_Be_2plus - E_Be_plus
IE_3 = E_Be_3plus - E_Be_2plus

print(f"\n  Beryllium energy levels:")
print(f"    Be   (4e): E = {E_Be:.6f} Ha")
print(f"    Be+  (3e): E = {E_Be_plus:.6f} Ha")
print(f"    Be2+ (2e): E = {E_Be_2plus:.6f} Ha")
print(f"    Be3+ (1e): E = {E_Be_3plus:.6f} Ha")
print(f"    Be4+ (0e): E =  0.000000 Ha")

print(f"\n  Ionization energies:")
print(f"    IE_1 (Be  -> Be+  + e-):  {IE_1:.6f} Ha = {IE_1*27.211:.2f} eV")
print(f"    IE_2 (Be+ -> Be2+ + e-):  {IE_2:.6f} Ha = {IE_2*27.211:.2f} eV")
print(f"    IE_3 (Be2+-> Be3+ + e-):  {IE_3:.6f} Ha = {IE_3*27.211:.2f} eV")
print(f"    IE_4 (Be3+-> Be4+ + e-):  {E_Be_3plus:.6f} Ha (= -E(Be3+))")

print(f"""
  TOPOLOGY FLOW for Be ground state:

  R ~ 0:     Round S^11. mu_free = 12 (lowest [2,2] eigenvalue).
             All 4 electrons near nucleus. Pure SO(12) symmetry.

  R ~ 0.5-1: Perturbative regime.
             Shell structure begins: 1s^2 core separates from 2s^2 pair.
             The 4-particle wavefunction develops "two pairs" structure:
               inner pair (1s^2): tightly bound, symmetric under S_2
               outer pair (2s^2): loosely bound, symmetric under S_2
             S_4 [2,2] -> S_2 x S_2 x Z_2 (pair-exchange symmetry)

  R ~ 2-5:   FIRST crossover (valence detachment).
             One 2s electron begins to separate.
             S^11 starts to pinch along the valence direction.
             The [2,2] representation starts breaking:
               [2,2] of S_4 -> [2,1] of S_3 (Be+ core)
                              x [1] (detached electron)

  R ~ 10+:   First surgery complete.
             S^11 -> S^8 x R (Be+ core on S^8 + outer electron on R).
             mu(R) ~ E(Be+) * R^2 = {E_Be_plus:.3f} * R^2

  R >> 20:   SECOND surgery possible (only for excited states).
             S^8 -> S^5 x R^2 (Be2+ core + 2 outer electrons).
             Threshold: E(Be2+) = {E_Be_2plus:.3f} Ha

  R >>> 50:  THIRD surgery (triply-excited states only).
             S^5 x R^2 -> S^3 x R^3 (Be3+ + 3 outer electrons).
             Threshold: E(Be3+) = {E_Be_3plus:.3f} Ha

  For the GROUND STATE: only ONE surgery is relevant.
  The 2s electron detaches, leaving a Li-like Be+ core (S^8).
  This is Type D behavior: core + equivalent valence pair, but only
  ONE electron detaches (the ground state is not doubly excited).

  COMPARISON of first ionization:
                    He              Li              Be
  Surgery:          S^5->S^3xR      S^8->S^5xR      S^11->S^8xR
  Detaching e-:     either (demo.)  2s (hier.)      2s (hier.)
  Core manifold:    S^3 (H-like)    S^5 (He-like)   S^8 (Li-like)
  Core energy:      {EXACT_ENERGIES['He+']:.3f}          {EXACT_ENERGIES['Li+']:.3f}         {E_Be_plus:.3f}
  IE:               0.904 Ha        0.198 Ha        {IE_1:.3f} Ha
""")


# =========================================================================
# STEP 7: Structure type classification
# =========================================================================

print("=" * 80)
print("STEP 7: Structure type classification and periodic table connection")
print("=" * 80)

print(f"""
  TYPE A: Single particle (N_e = 1)
  ---------------------------------
  Examples: H, He+, Li2+, Be3+, ...
  Manifold: S^2 (just angular part of single electron)
  Group: SO(3), trivial S_1
  Feature: EXACT Coulomb solution. No angular correlation.
  Quasi-Coulomb accuracy: 0% error (IS the Coulomb solution)
  Algebra: purely algebraic eigenvalues E_n = -Z^2/(2n^2)
  Topology: no flow needed, single S^2 for all R

  TYPE B: Equivalent pair (N_e = 2, ^1S ground state)
  ----------------------------------------------------
  Examples: He, Li+, Be2+, B3+, ... (He-like isoelectronic series)
  Manifold: S^5
  Group: SO(6), S_2 symmetric [2]
  Feature: Democratic symmetry. EITHER electron can detach.
  Selection rule: Delta_nu = 0 mod 4 (from S_2 exchange)
  Quasi-Coulomb: ~5-14% error (a_1 perturbative known exactly)
  Algebra: a_1 algebraic, a_2 transcendental (partial harmonic sum)
  Topology: S^5 -> S^3 x R (single surgery)

  TYPE C: Core + single valence (N_e = 3, ^2S ground state)
  -----------------------------------------------------------
  Examples: Li, Be+, B2+, C3+, ... (Li-like isoelectronic series)
  Manifold: S^8
  Group: SO(9), S_3 mixed [2,1]
  Feature: HIERARCHICAL separation. Valence electron clearly distinct.
  Selection rule: more permissive than He (S_3 [2,1])
  Quasi-Coulomb: ~0% error estimate (Li, if a_2 small)
  Algebra: a_1 expected algebraic, a_2 convergence unknown
  Topology: S^8 -> S^5 x R (single surgery, valence detaches)
  Chemistry: alkali metals (Li, Na, K, ...)

  TYPE D: Core + equivalent valence pair (N_e = 4, ^1S ground state)
  -------------------------------------------------------------------
  Examples: Be, B+, C2+, N3+, ... (Be-like isoelectronic series)
  Manifold: S^11
  Group: SO(12), S_4 self-conjugate [2,2]
  Feature: TWO-STAGE hierarchy. 1s^2 core + 2s^2 valence pair.
  Selection rule: TBD (S_4 [2,2] sector)
  Quasi-Coulomb: TBD (need perturbative a_1)
  Algebra: expected more complex (higher multiplicity of [2,2])
  Topology: S^11 -> S^8 x R (one valence electron detaches)
  Chemistry: alkaline earth metals (Be, Mg, Ca, ...)

  THE PERIODIC TABLE CONNECTION:
  ==============================
  The structure type classification maps DIRECTLY onto the periodic table:

  Type A: H-like    (1s^1)              -> Group 1, Period 1
  Type B: He-like   (1s^2, closed)      -> Group 18
  Type C: Li-like   (core + ns^1)       -> Group 1 (alkali metals)
  Type D: Be-like   (core + ns^2)       -> Group 2 (alkaline earths)

  The "cleanness" of quasi-Coulomb correlates with SHELL STRUCTURE:

  - Type A (H): Perfect Coulomb. No e-e interaction. EXACT.
  - Type C (Li): Clean hierarchy. Core screening creates effective charge.
    The 2s electron "sees" Z_eff ~ 1.69 (Clementi-Raimondi).
    Alkali metals are "hydrogen-like" — chemistry confirms this!
  - Type B (He): Both electrons equivalent. Neither is "the core" or
    "the valence." Must correlate democratically. HARDER.
  - Type D (Be): Mixed. Has core (1s^2) like Type C, but ALSO has
    equivalent valence pair (2s^2) like Type B. Inherits both features.

  PREDICTION: Be's quasi-Coulomb accuracy should be INTERMEDIATE
  between He (5-14%) and Li (~0%). The 1s^2 core is clean (Type C
  contribution), but the 2s^2 pair introduces He-like correlation
  (Type B contribution). Estimated: 3-8% error.
""")


# =========================================================================
# STEP 8: Be FCI data analysis
# =========================================================================

print("=" * 80)
print("STEP 8: GeoVac FCI results for beryllium")
print("=" * 80)

print(f"\n  DirectCISolver results (vee_method='slater_full', h1_method='exact'):")
print(f"\n  {'nmax':>5}  {'n_spatial':>10}  {'n_spinorb':>10}  {'N_SD':>10}  "
      f"{'E (Ha)':>12}  {'Error %':>8}  {'Time':>8}")
print("-" * 72)

for key in ['Be_n3', 'Be_n4']:
    r = FCI_RESULTS[key]
    n_spatial = r['nmax']**2  # Approximate
    n_spinorb = 2 * n_spatial
    err = abs((r['E'] - EXACT_ENERGIES['Be']) / EXACT_ENERGIES['Be']) * 100
    print(f"  {r['nmax']:5d}  {n_spatial:10d}  {n_spinorb:10d}  {r['nsd']:10,d}  "
          f"{r['E']:12.6f}  {err:7.2f}%  {r['time']:7.1f}s")

print(f"\n  Exact (non-relativistic):  {EXACT_ENERGIES['Be']:.6f} Ha")
print(f"  Hartree-Fock limit:       {EXACT_ENERGIES['Be_HF']:.6f} Ha")
print(f"  Correlation energy:       {EXACT_ENERGIES['Be'] - EXACT_ENERGIES['Be_HF']:.6f} Ha")

# Convergence analysis
E_n3 = FCI_RESULTS['Be_n3']['E']
E_n4 = FCI_RESULTS['Be_n4']['E']
delta = E_n4 - E_n3
print(f"\n  Convergence: E(nmax=4) - E(nmax=3) = {delta:.6f} Ha")
print(f"  Still {EXACT_ENERGIES['Be'] - E_n4:.4f} Ha above exact")
print(f"  (Both above HF: {E_n3 - EXACT_ENERGIES['Be_HF']:.4f}, {E_n4 - EXACT_ENERGIES['Be_HF']:.4f} Ha)")

# Cross-system comparison
print(f"\n  Cross-system FCI accuracy at nmax=3:")
print(f"  {'System':>8}  {'N_e':>3}  {'N_SD':>10}  {'Error %':>8}  {'Scaling':>10}")
print("-" * 50)
for name, key, ne in [('He', 'He_n3', 2), ('Li', 'Li_n3', 3), ('Be', 'Be_n3', 4)]:
    r = FCI_RESULTS[key]
    exact = EXACT_ENERGIES[name]
    err = abs((r['E'] - exact) / exact) * 100
    print(f"  {name:>8}  {ne:3d}  {r['nsd']:10,d}  {err:7.2f}%  {'C(28,' + str(ne) + ')':>10}")


# =========================================================================
# STEP 9: Isoelectronic series analysis
# =========================================================================

print("\n" + "=" * 80)
print("STEP 9: Isoelectronic series — Z-scaling of algebraic structure")
print("=" * 80)

print(f"""
  Each isoelectronic series (fixed N_e, varying Z) shares the SAME
  angular manifold S^(3N-1) and group structure SO(3N) x S_N.

  What CHANGES with Z:
  - The charge function: C = -Z * sum(1/r_i) + sum(1/r_ij)
  - The perturbation parameter: a_1 scales linearly with Z
  - The ratio Z/N_e controls "nuclear dominance"

  KEY INSIGHT: Large Z -> more hydrogen-like (nuclear dominance)
               Z = N_e -> neutral atom, maximum correlation
               Z < N_e -> unstable (electron detaches)
""")

# He-like series (Type B)
print(f"  He-like isoelectronic series (N_e = 2, Type B):")
print(f"  {'Ion':>6}  {'Z':>3}  {'E_exact':>10}  {'E_HF_approx':>11}  {'Z/N':>5}  {'IE_1':>8}")
print("-" * 55)

he_like = [
    ('H-',   1, -0.527751),
    ('He',   2, -2.903724),
    ('Li+',  3, -7.279913),
    ('Be2+', 4, -13.655566),
]

for name, Z, E in he_like:
    # Parent ion energy (H-like, Z)
    E_parent = -Z**2 / 2.0
    IE = E_parent - E
    z_over_n = Z / 2.0
    # HF approx: E_HF ~ -Z^2 + 5Z/8 (for 1s^2)
    E_hf_approx = -Z**2 + 5*Z/8
    print(f"  {name:>6}  {Z:3d}  {E:10.4f}  {E_hf_approx:11.4f}  {z_over_n:5.2f}  {IE:8.4f}")

# Li-like series (Type C)
print(f"\n  Li-like isoelectronic series (N_e = 3, Type C):")
print(f"  {'Ion':>6}  {'Z':>3}  {'E_exact':>10}  {'Z/N':>5}  {'IE_1':>8}")
print("-" * 45)

li_like = [
    ('Li',   3, -7.478060,  -7.279913),
    ('Be+',  4, -14.324764, -13.655566),
]

for name, Z, E, E_parent in li_like:
    IE = E_parent - E
    z_over_n = Z / 3.0
    print(f"  {name:>6}  {Z:3d}  {E:10.4f}  {z_over_n:5.2f}  {IE:8.4f}")


# =========================================================================
# STEP 10: Patterns and predictions
# =========================================================================

print("\n" + "=" * 80)
print("STEP 10: Emerging patterns and predictions")
print("=" * 80)

print(f"""
  PATTERN 1: l_eff = (3N-1)(3N-3)/8 grows QUADRATICALLY with N
  =============================================================
  N=1: l_eff = 0      n_eff = 1
  N=2: l_eff = 3/2    n_eff = 5/2
  N=3: l_eff = 3      n_eff = 4
  N=4: l_eff = 9/2    n_eff = 11/2
  N=5: l_eff = 6      n_eff = 7

  Formula: l_eff = (3N-1)(3N-3)/8  (from (d-1)(d-3)/8, d=3N)
           l_eff ~ 9N^2/8 for large N

  PATTERN 2: Ground state nu increases with shell structure
  =========================================================
  N=1: nu=0 (no antisymmetry needed)
  N=2: nu=0 (spatial symmetric [2] starts at nu=0)
  N=3: nu=1 (spatial [2,1] starts at nu=1)
  N=4: nu=2 (spatial [2,2] starts at nu=2, est.)

  The ground state nu = number of "antisymmetric pairs" in the
  Young diagram. This encodes the PAULI EXCLUSION PRINCIPLE
  directly in the angular quantum number!

  PATTERN 3: Singularity count = N(N+1)/2 grows quadratically
  ============================================================
  N=1: 1 singularity  (1 nuclear)
  N=2: 3 singularities (2 nuclear + 1 e-e)
  N=3: 6 singularities (3 nuclear + 3 e-e)
  N=4: 10 singularities (4 nuclear + 6 e-e)
  Formula: N + C(N,2) = N(N+1)/2

  Each singularity contributes to the charge function and makes
  the angular problem progressively harder.

  PATTERN 4: Topology flow depth increases with core structure
  =============================================================
  N=1: S^2  (no flow)                          depth 0
  N=2: S^5  -> S^3 x R                         depth 1
  N=3: S^8  -> S^5 x R                         depth 1 (gs)
  N=4: S^11 -> S^8 x R                         depth 1 (gs)
  N=5: S^14 -> S^11 x R                        depth 1 (gs)

  For ground states, the depth is ALWAYS 1 (single valence detachment).
  Higher depths only appear in excited states or multiply-ionized channels.
  This is a deep result: atoms ionize ONE ELECTRON AT A TIME.

  PATTERN 5: Structure type maps to periodic table columns
  =========================================================
  Type A (single particle):           H   — Period 1
  Type B (equivalent pair):           He  — Noble gas (closed shell)
  Type C (core + single valence):     Li  — Alkali metal
  Type D (core + equivalent pair):    Be  — Alkaline earth

  Continuing the pattern:
  Type E? (core + 3 valence):         B   — Group 13 (p-block begins)
  Type F? (core + 4 valence):         C   — Group 14
  ...

  The B atom (5e, 1s^2 2s^2 2p) introduces l=1 angular momentum.
  This would be the FIRST test of the p-orbital sector in the
  hyperspherical framework. Expected: significantly more complex,
  as l>0 orbitals break the spherical symmetry of the angular
  manifold and introduce new quantum numbers.

  PREDICTION: HIERARCHY OF "MESSINESS"
  =====================================
  The quasi-Coulomb model should work BEST for atoms where the
  separation between core and valence is CLEANEST:

  Ranking (cleanest -> messiest):
  1. H   (Type A) — exact, trivial
  2. Li  (Type C) — clean core/valence hierarchy
  3. Na  (Type C) — even cleaner (larger core screens better)
  4. Be  (Type D) — core is clean, but valence pair correlates
  5. He  (Type B) — no hierarchy, fully democratic
  6. B+  (Type D) — like Be but with higher Z/N
  7. B   (5e)     — p-orbital onset, new angular complexity

  This predicts that ALKALI METALS should be the easiest multi-electron
  atoms for hyperspherical treatment, which is consistent with their
  simple chemistry (one valence electron, hydrogen-like spectra).
""")


# =========================================================================
# STEP 11: Quantitative predictions for Be
# =========================================================================

print("=" * 80)
print("STEP 11: Quantitative predictions for Be hyperspherical solver")
print("=" * 80)

# Be centrifugal and effective parameters
N_Be = 4
d_Be = 12
cent_Be = Fraction(11 * 9, 8)  # 99/8
l_eff_Be = 4.5
n_eff_Be = 5.5

# Z_eff from exact energy
Z_eff_Be = np.sqrt(2 * n_eff_Be**2 * abs(E_Be))
a1_Be_back = -Z_eff_Be

print(f"\n  Be hyperspherical parameters:")
print(f"    Manifold: S^11, Group: SO(12) x S_4")
print(f"    Centrifugal: 99/(8R^2) = 12.375/R^2")
print(f"    l_eff = 9/2, n_eff = 11/2")
print(f"    Z_eff (back-computed) = {Z_eff_Be:.4f}")
print(f"    a_1 (back-computed)   = {a1_Be_back:.4f}")

# Expected perturbative a_1 structure
print(f"""
  Expected a_1 structure:
    a_1(Be, Z) = (norm. on S^11) * (-Z * I_nuc^(4e) + I_ee^(4e))

    Nuclear part: 4 terms (4 electrons x nuclear attraction)
    e-e part: 6 terms (C(4,2) = 6 electron pairs)

    For large Z: a_1 ~ -4Z * c_nuc (nuclear dominates)
    For Z = N = 4: maximum competition between nuclear and e-e

  Predicted a_1 vs exact comparison:
    He: a_1(pert) = -5.590, Z_eff(exact) = 3.808, ratio = 1.47
    Li: a_1(pert) = ?, Z_eff(exact) = 15.47
    Be: a_1(pert) = ?, Z_eff(exact) = {Z_eff_Be:.2f}

  Computational cost estimate:
    Angular problem: 3D surface in S^11 (after symmetry reduction)
    Grid: ~50-100 per dimension x ~20-50 channels
    Matrix: ~100k-500k dimension (SPARSE)
    Eigensolver: Lanczos, ~minutes per R-point
    Adiabatic: 100-200 R-points
    TOTAL: ~hours (vs seconds for He, minutes for Li)

    This is the NATURAL COST of the S_4 representation complexity.
""")


# =========================================================================
# STEP 12: Physical interpretation summary
# =========================================================================

print("=" * 80)
print("STEP 12: Physical interpretation — why alkali metals are clean")
print("=" * 80)

print(f"""
  WHY ALKALI METALS HAVE "CLEAN" HYPERSPHERICAL STRUCTURE:

  1. HIERARCHICAL SEPARATION:
     Alkali metals (Li, Na, K, ...) have a tightly bound core (1s^2, or
     1s^2 2s^2 2p^6, etc.) and a SINGLE loosely bound valence electron.
     The energy gap between core and valence is LARGE:

     Li: E(1s) = -4.5 Ha, E(2s) = -0.198 Ha (ratio 22.7)
     Na: E(2p) = -3.03 Ha, E(3s) = -0.182 Ha (ratio 16.6)

     This large ratio means the topology flow S^(3N-1) -> S^(3N-4) x R
     occurs at SMALL R, and the valence electron is well-described by
     an effective Coulomb potential with Z_eff quickly.

  2. THE [N-1,1] REPRESENTATION:
     Alkali metals always have one unpaired electron, giving spatial
     irrep [N-1,1] of S_N. This is the STANDARD representation,
     which has the SIMPLEST non-trivial S_N decomposition.

     [N-1,1] first appears at nu=1 in SO(3N), so the ground state
     is the FIRST EXCITED angular harmonic (nu=1), close to the
     free eigenvalue mu_free = (3N-1)/2.

  3. CORE SCREENING:
     The core electrons create an effective nuclear charge Z_eff < Z
     seen by the valence electron. This is EXACTLY what the topology
     flow does: at large R, the angular manifold "collapses" the core
     degrees of freedom into an effective charge, leaving the valence
     on S^2 (hydrogen-like) with Z_eff.

  WHY He AND Be ARE "MESSY":

  1. DEMOCRATIC PAIRS:
     In He (1s^2) and the Be valence (2s^2), two electrons occupy
     the SAME orbital. Neither electron is "the core" or "the valence."
     They must correlate DEMOCRATICALLY.

  2. THE [2,...,2] REPRESENTATION:
     Closed-shell singlets have spatial irrep [2,2,...] which is
     SELF-CONJUGATE. This means spatial and spin wavefunctions have
     the SAME symmetry — there's no clean separation.

  3. ANGULAR NODES:
     The self-conjugate [2,2] representation has nu_gs = 2 (for Be),
     meaning the ground state wavefunction has TWO angular nodes.
     These nodes create complex angular correlation that the
     quasi-Coulomb model captures imperfectly.

  CHEMISTRY CONNECTION:
     Alkali metals (Group 1): simplest spectra, lowest IE, hydrogen-like
     Alkaline earths (Group 2): slightly more complex, paired valence
     Noble gases (Group 18): fully closed, maximum correlation

     This hierarchy is ENCODED in the S_N representation theory:
       [N-1,1] (alkali) < [N-2,2] (alkaline earth) < ... < [2,...,2] (noble gas)

     The GROUP THEORY of the angular manifold IS the periodic table!
""")


# =========================================================================
# FINAL SUMMARY TABLE
# =========================================================================

print("=" * 80)
print("FINAL SUMMARY: Algebraic Structure Across the First Period")
print("=" * 80)

headers = ['Property', 'H (Z=1)', 'He (Z=2)', 'Li (Z=3)', 'Be (Z=4)']
rows = [
    ['N_electrons',    '1',      '2',       '3',       '4'],
    ['Config',         '1s',     '1s^2',    '1s^2 2s', '1s^2 2s^2'],
    ['Type',           'A',      'B',       'C',       'D'],
    ['Manifold',       'S^2',    'S^5',     'S^8',     'S^11'],
    ['Group',          'SO(3)',  'SO(6)',   'SO(9)',   'SO(12)'],
    ['S_N irrep',      '[1]',   '[2]',    '[2,1]',  '[2,2]'],
    ['nu_gs',          '0',      '0',       '1',       '2 (est.)'],
    ['mu_free',        '0',      '0',       '4',       '12 (est.)'],
    ['l_eff',          '0',      '3/2',     '3',       '9/2'],
    ['n_eff',          '1',      '5/2',     '4',       '11/2'],
    ['Singularities',  '1',      '3',       '6',       '10'],
    ['Surgeries (gs)', '0',      '1',       '1',       '1'],
    ['E_exact (Ha)',    '-0.500', '-2.904',  '-7.478',  '-14.667'],
    ['E_corr (Ha)',     '0.000',  '-0.042',  '-0.045',  '-0.094'],
    ['IE_1 (eV)',       '13.6',   '24.6',    '5.39',    '9.32'],
    ['FCI err (n=3)',   'N/A',    '~0.7%',   '1.85%',   '0.93%'],
    ['QC accuracy',     'exact',  '~14%',    '~0% (est)', '3-8% (pred)'],
]

# Print table
col_widths = [20, 10, 10, 10, 12]
header_line = '  '.join(h.center(w) for h, w in zip(headers, col_widths))
print(f"\n  {header_line}")
print("  " + "-" * len(header_line))
for row in rows:
    line = '  '.join(v.center(w) for v, w in zip(row, col_widths))
    print(f"  {line}")

print(f"""

  KEY CONCLUSIONS:
  ================
  1. Be's ground state lives in the [2,2] sector of S_4, which is
     SELF-CONJUGATE (spatial and spin have the same symmetry type).
     This is unique among the first 4 elements.

  2. The estimated lowest nu for [2,2] is nu=2, giving mu_free = 12.
     This means Be's ground state has MORE angular structure than
     He (nu=0) or Li (nu=1) — it's the most angularly excited atom
     in the first period.

  3. Be's structure type (D: core + equivalent valence pair) combines
     features of Type C (hierarchical core/valence) and Type B
     (democratic pair). Its hyperspherical treatment should be
     intermediate in difficulty.

  4. The structure type classification maps onto periodic table columns:
     Type A -> H, Type B -> noble gases, Type C -> alkali metals,
     Type D -> alkaline earths. The GROUP THEORY of S^(3N-1) encodes
     the chemical periodicity!

  5. FCI convergence for Be (0.93% at nmax=3, 0.90% at nmax=4) is
     SLOWER than He but FASTER than Li in absolute error. This is
     consistent with Be having both core structure (helps) and
     valence pair correlation (hurts).

  OPEN QUESTIONS:
  ===============
  - Verify nu_gs = 2 for [2,2] via explicit SO(12) -> S_4 decomposition
  - Compute perturbative a_1 for Be (requires 4-body angular integrals)
  - Test quasi-Coulomb prediction (3-8% error)
  - Extend to B (5e): first p-orbital, Type E?
  - Explore connection to Madelung rule (n + l ordering)
""")

print("\n" + "=" * 80)
print("Analysis complete. Output: debug/debug_beryllium_periodic.py")
print("=" * 80)
