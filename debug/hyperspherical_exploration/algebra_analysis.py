"""
Algebraic Structure Analysis: Hyperspherical Lattice for Helium

Date: 2026-03-14
Status: Phase 1 — Mathematical Reconnaissance

Investigates the symmetry group, irreducible representations, ladder operators,
and degeneracy structure of the two-electron hyperspherical system, and their
connection to Paper 0's S³ packing.

This script performs symbolic/numerical checks — it does NOT build a solver.
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
from typing import Tuple, List, Dict
from itertools import product


# =============================================================================
# 1. SYMMETRY GROUP ANALYSIS
# =============================================================================

def print_symmetry_analysis() -> None:
    """
    Document the symmetry group structure of the two-electron problem
    in hyperspherical coordinates.
    """
    print("=" * 70)
    print("SYMMETRY GROUP ANALYSIS: He in Hyperspherical Coordinates")
    print("=" * 70)

    print("""
1. FULL SYMMETRY GROUP
   G = SO(3)_L × S₂ × Z₂^π × SU(2)_S

   Components:
   - SO(3)_L : Total orbital angular momentum (L, M conserved)
   - S₂      : Electron exchange symmetry (r₁ ↔ r₂)
   - Z₂^π    : Spatial parity (r → -r for both electrons)
   - SU(2)_S : Total spin (S = 0 singlet, S = 1 triplet)

2. HYPERANGULAR SYMMETRY (at fixed R)
   The Laplace-Beltrami operator Λ² on S⁵ has symmetry group SO(6).

   SO(6) ⊃ SO(3)₁ × SO(3)₂     (individual angular momenta l₁, l₂)
   SO(6) ⊃ SO(3)_L × SO(3)_rel  (total + relative angular momentum)

   The Coulomb potential BREAKS SO(6) down to:
   SO(6) → SO(3)_L × Z₂^{exchange} × Z₂^{parity}

   Only L, M, S, π are exact quantum numbers.
   K (grand angular momentum) is approximate — broken by C(α, θ₁₂)/R.

3. COMPARISON WITH ATOMIC S³
   S³ (one electron):  SO(4) = SU(2) × SU(2)  — Fock's hydrogen symmetry
   S⁵ (two electrons): SO(6) ⊃ SO(3) × SO(3)  — hyperspherical symmetry

   The S³ lattice exploits SO(4) with nodes = (n, l, m).
   The S⁵ lattice would exploit SO(6) with nodes = (K, l₁, l₂, L, M).

4. FOR He GROUND STATE (¹S: L=0, S=0, π=+1)
   - l₁ = l₂ (forced by L=0 and Clebsch-Gordan coupling)
   - Exchange symmetry: Φ(α) = Φ(π/2 - α) (symmetric in α)
   - Only even-l partial waves contribute (parity constraint)
   - Effective quantum numbers: (K, l) where K = 2n_α + 2l
""")


# =============================================================================
# 2. HYPERSPHERICAL HARMONICS AND DEGENERACY
# =============================================================================

def degeneracy_S5(K: int) -> int:
    """
    Total degeneracy of grand angular momentum shell K on S⁵.

    The hyperspherical harmonics Y_{K,l₁,l₂,L,M} on S⁵ have degeneracy:
        g(K) = (K+1)(K+2)²(K+3) / 12

    This counts ALL (l₁, l₂, L, M) states at given K.
    """
    return (K + 1) * (K + 2)**2 * (K + 3) // 12


def degeneracy_S3(n: int) -> int:
    """
    Total degeneracy of principal quantum number n on S³.
    g(n) = 2n² (including spin).
    """
    return 2 * n**2


def degeneracy_S5_L0(K: int) -> int:
    """
    Number of L=0 states at given K on S⁵.

    For L=0: l₁ = l₂ = l, and K = 2n_α + 2l.
    So l ranges from 0 to K//2 (only even K contribute).
    For each valid l, there is exactly one (l₁=l, l₂=l, L=0, M=0) state.

    For even K: number of L=0 states = K//2 + 1
    For odd K: 0 (no L=0 states at odd K)
    """
    if K % 2 != 0:
        return 0
    return K // 2 + 1


def print_degeneracy_table() -> None:
    """Print degeneracy comparison: S³ vs S⁵ vs S⁵(L=0)."""
    print("\n" + "=" * 70)
    print("DEGENERACY STRUCTURE: S³ (atomic) vs S⁵ (two-electron)")
    print("=" * 70)

    print(f"\n{'K/n':>4} {'S³ g(n)':>10} {'S⁵ g(K)':>10} {'S⁵ L=0':>10} "
          f"{'Cumul S³':>10} {'Cumul S⁵':>10}")
    print("-" * 60)

    cum_s3, cum_s5 = 0, 0
    for K in range(8):
        n = K + 1  # S³ quantum number for comparison
        g3 = degeneracy_S3(n)
        g5 = degeneracy_S5(K)
        g5_L0 = degeneracy_S5_L0(K)
        cum_s3 += g3
        cum_s5 += g5
        print(f"{K:>4} {g3:>10} {g5:>10} {g5_L0:>10} {cum_s3:>10} {cum_s5:>10}")

    print("""
Key observations:
- S³ degeneracy grows as 2n²  (hydrogen-like)
- S⁵ degeneracy grows as K⁴/12 (much faster — 6D phase space)
- For L=0 (He ground state), only EVEN K contribute
- The L=0 sector is sparse: ~1/K³ of total states
- This sparsity is what makes the adiabatic method efficient
""")


# =============================================================================
# 3. IRREDUCIBLE REPRESENTATIONS
# =============================================================================

def enumerate_states(K: int, L_target: int = 0) -> List[Tuple[int, int, int, int]]:
    """
    Enumerate all (K, l₁, l₂, L) states at given K and L.

    Constraint: K = 2n_α + l₁ + l₂ where n_α ≥ 0
    So l₁ + l₂ ≤ K and (K - l₁ - l₂) must be even.

    Returns list of (n_α, l₁, l₂, n_channels) tuples.
    """
    states = []
    for l1 in range(K + 1):
        for l2 in range(K - l1 + 1):
            # Check K parity
            if (K - l1 - l2) % 2 != 0:
                continue
            n_alpha = (K - l1 - l2) // 2
            # Check L = |l1-l2|, ..., l1+l2
            if abs(l1 - l2) <= L_target <= l1 + l2:
                states.append((n_alpha, l1, l2, L_target))
    return states


def print_irrep_decomposition() -> None:
    """Print the irreducible representation decomposition for L=0."""
    print("\n" + "=" * 70)
    print("IRREDUCIBLE REPRESENTATIONS: L=0 sector")
    print("=" * 70)

    for K in range(8):
        states = enumerate_states(K, L_target=0)
        if not states:
            continue
        print(f"\nK = {K}: {len(states)} L=0 state(s)")
        for n_a, l1, l2, L in states:
            # Exchange symmetry: for L=0, singlet requires (-1)^{l1+l2+S} = +1
            # S=0 (singlet): l1+l2 must be even → l1=l2 is always even sum
            exchange = "singlet" if (l1 + l2) % 2 == 0 else "triplet"
            print(f"  n_α={n_a}, l₁=l₂={l1}, {exchange}")

    print("""
For ¹S He (singlet, L=0):
  - Only l₁ = l₂ = l states contribute
  - K must be even
  - States: K=0 (l=0), K=2 (l=0,1), K=4 (l=0,1,2), ...
  - These form the CHANNELS in the adiabatic expansion
  - Channel ν at fixed R is a linear combination over these (K, l) states
""")


# =============================================================================
# 4. LADDER OPERATORS AND GRAPH EDGES
# =============================================================================

def print_ladder_operators() -> None:
    """Analyze the ladder operators that would define hyperspherical lattice edges."""
    print("\n" + "=" * 70)
    print("LADDER OPERATORS AND GRAPH STRUCTURE")
    print("=" * 70)

    print("""
1. HYPERANGULAR LADDER OPERATORS

   On S⁵, the SO(6) algebra has 15 generators (dim SO(6) = 6·5/2 = 15).
   These decompose under SO(3)₁ × SO(3)₂ as:

   - 3 generators of SO(3)₁:  L₁±, L₁z  (electron 1 angular momentum)
   - 3 generators of SO(3)₂:  L₂±, L₂z  (electron 2 angular momentum)
   - 9 remaining "boost" generators connecting different (l₁, l₂) sectors

   The boost generators act as:
     B_{ij}: changes l₁ by ±1 AND l₂ by ∓1 (preserves K)
     A_{ij}: changes l₁ by ±1 AND l₂ by ±1 (changes K by ±2)

2. HYPERRADIAL CONNECTIONS

   In the hyperradial direction, edges come from:
   - Finite-difference stencil: R_i ↔ R_{i±1} (kinetic energy)
   - Non-adiabatic coupling P_μν(R): channel ν ↔ channel μ at same R

3. GRAPH STRUCTURE

   Nodes: (R_i, ν)  where i = radial grid index, ν = channel index

   Edges:
   ┌────────────────────────────────────────────────────────────┐
   │ TYPE 1: Radial kinetic (within channel)                   │
   │   (R_i, ν) ── (R_{i±1}, ν)                               │
   │   Weight: FD stencil coefficient (like prolate ξ edges)   │
   │                                                           │
   │ TYPE 2: Non-adiabatic coupling (between channels)         │
   │   (R_i, ν) ── (R_i, μ)                                   │
   │   Weight: P_μν(R_i) or Q_μν(R_i)                         │
   │   These are R-DEPENDENT (unlike prolate η coupling which  │
   │   is captured in the spectral solve)                      │
   │                                                           │
   │ TYPE 3: Adiabatic potential (diagonal)                    │
   │   (R_i, ν) self-loop with weight U_ν(R_i)/R_i            │
   │   This is the effective potential from the angular solve   │
   └────────────────────────────────────────────────────────────┘

4. COMPARISON WITH S³ AND PROLATE LATTICES

   S³ lattice:
     Nodes = (n, l, m)
     Edges = L± (l↔l±1), radial (n↔n±1)
     Angular: FULLY SOLVED by S³ harmonics (Gegenbauer)
     Radial: implicit in the S³ eigenvalue

   Prolate lattice:
     Nodes = (ξ_i, n_η)
     Edges = FD in ξ, Legendre spectral in η
     Angular (η): spectral solve → separation constant A
     Radial (ξ): FD solve → energy E
     Coupling: NONE (exactly separable for H₂⁺)

   Hyperspherical lattice:
     Nodes = (R_i, ν)
     Edges = FD in R, non-adiabatic P_μν between channels
     Angular (α, θ₁₂): coupled spectral solve → U_μ(R)
     Radial (R): FD solve → energy E
     Coupling: P_μν(R) — NON-ZERO (not separable)

   KEY STRUCTURAL DIFFERENCE:
   The hyperspherical lattice has R-dependent inter-channel edges.
   This is fundamentally different from the prolate lattice (no coupling)
   and the S³ lattice (no radial problem). It is a COUPLED CHANNEL GRAPH.

5. IS THIS STILL A "DISCRETE GRAPH WITH QUANTUM NUMBER LABELS"?

   YES — but with a crucial generalization:
   - The angular quantum numbers are ADIABATIC channel indices ν,
     not fixed quantum numbers like (n, l, m)
   - The channel functions Φ_ν(R; Ω̂) change with R (parametric dependence)
   - The inter-channel edges P_μν(R) are R-dependent

   This is the graph-theoretic analogue of a FIBER BUNDLE:
   - Base space: R grid (hyperradial)
   - Fiber: channel space {ν} at each R
   - Connection: P_μν(R) (non-adiabatic coupling = Berry connection)

   The S³ lattice is a TRIVIAL bundle (no R dependence).
   The prolate lattice is also trivial (exact separation → no coupling).
   The hyperspherical lattice is the first NON-TRIVIAL bundle in GeoVac.
""")


# =============================================================================
# 5. CHARGE FUNCTION ANALYSIS
# =============================================================================

def charge_function(alpha: np.ndarray, theta12: np.ndarray,
                    Z: float = 2.0) -> np.ndarray:
    """
    Compute the angular charge function C(α, θ₁₂).

    C(α, θ₁₂) = -Z/cos(α) - Z/sin(α) + 1/√(1 - sin(2α)cos(θ₁₂))

    Parameters
    ----------
    alpha : array, hyperangle ∈ (0, π/2)
    theta12 : array, interelectron angle ∈ [0, π]
    Z : nuclear charge

    Returns
    -------
    C : array, charge function values
    """
    r12_factor = np.sqrt(np.maximum(1 - np.sin(2 * alpha) * np.cos(theta12),
                                     1e-15))
    return -Z / np.cos(alpha) - Z / np.sin(alpha) + 1.0 / r12_factor


def analyze_charge_function(Z: float = 2.0) -> None:
    """Analyze critical points of the charge function for He."""
    print("\n" + "=" * 70)
    print(f"CHARGE FUNCTION ANALYSIS: Z = {Z}")
    print("=" * 70)

    # Wannier saddle point: α = π/4, θ₁₂ = π
    alpha_s = np.pi / 4
    C_saddle = -2 * Z * np.sqrt(2) + 1.0 / np.sqrt(2)
    print(f"\nWannier saddle (α=π/4, θ₁₂=π): C = {C_saddle:.6f}")

    # Coalescence point: α = π/4, θ₁₂ = 0
    # C diverges (1/r₁₂ → ∞), but this is the Kato cusp
    print(f"Coalescence (α=π/4, θ₁₂=0):    C → +∞ (electron-electron cusp)")

    # Nuclear singularities
    print(f"Nuclear singularity (α→0):       C → -∞ (electron 1 at nucleus)")
    print(f"Nuclear singularity (α→π/2):     C → -∞ (electron 2 at nucleus)")

    # Minimum of C at θ₁₂ = π (electrons on opposite sides)
    # Find minimum over α at θ₁₂ = π
    alpha_arr = np.linspace(0.01, np.pi / 2 - 0.01, 1000)
    C_pi = charge_function(alpha_arr, np.full_like(alpha_arr, np.pi), Z)
    idx_min = np.argmin(C_pi)
    print(f"\nMinimum of C at θ₁₂=π: C = {C_pi[idx_min]:.6f} at α = {alpha_arr[idx_min]:.4f}")
    print(f"  (α = π/4 = {np.pi/4:.4f} → symmetric configuration)")

    # Effective potential at representative R values
    print("\nEffective potential U_eff ≈ K(K+4)/2 + R·C_min:")
    for R in [0.5, 1.0, 2.0, 5.0]:
        for K in [0, 2, 4]:
            U_eff = K * (K + 4) / 2 + R * C_saddle  # rough estimate
            print(f"  R={R:.1f}, K={K}: U_eff/R ≈ {U_eff/R:.3f} Ha")


# =============================================================================
# 6. CHANNEL COUPLING MATRIX ELEMENTS
# =============================================================================

def gaunt_integral(l1: int, l2: int, l3: int) -> float:
    """
    Compute the Gaunt integral ∫ P_{l1}(x) P_{l2}(x) P_{l3}(x) dx
    over [-1, 1].

    Uses the relation to Wigner 3j symbols:
    ∫ P_l P_k P_{l'} dx = 2 * (l  k  l')² (3j symbol squared)
                                (0  0  0 )

    Selection rules: l1+l2+l3 must be even, triangle inequality.
    """
    # Selection rules
    if (l1 + l2 + l3) % 2 != 0:
        return 0.0
    if l3 > l1 + l2 or l3 < abs(l1 - l2):
        return 0.0

    # Use the formula for the Wigner 3j symbol (0 0 0) case
    # (l1 l2 l3) = (-1)^g * g! / ((g-l1)!(g-l2)!(g-l3)!) * √(...)
    # (0  0  0 )
    # where g = (l1+l2+l3)/2

    from math import factorial, sqrt

    g = (l1 + l2 + l3) // 2
    if g < l1 or g < l2 or g < l3:
        return 0.0

    # 3j symbol squared for m₁=m₂=m₃=0
    num = factorial(2 * (g - l1)) * factorial(2 * (g - l2)) * factorial(2 * (g - l3))
    den = factorial(2 * g + 1)
    threej_sq = (factorial(g) ** 2 * num) / (
        factorial(g - l1) ** 2 * factorial(g - l2) ** 2 * factorial(g - l3) ** 2 * den
    )

    return 2.0 * threej_sq


def print_coupling_matrix() -> None:
    """Print the V_ee coupling matrix elements between L=0 partial waves."""
    print("\n" + "=" * 70)
    print("V_ee COUPLING MATRIX: L=0 partial waves")
    print("=" * 70)

    l_max = 6
    print(f"\nGaunt integrals ⟨P_l|P_k|P_l'⟩ for k = 0, 1, 2, ... (V_ee multipoles)")
    print(f"Only k values appearing in the 1/r₁₂ expansion are relevant.\n")

    # For L=0: l₁ = l₂ = l, and coupling is through k-th multipole
    print(f"{'l':>3} {'l_prime':>7} {'k':>3} {'Gaunt':>12} {'Selection':>12}")
    print("-" * 45)

    for l in range(l_max):
        for lp in range(l, l_max):
            for k in range(l + lp + 1):
                g = gaunt_integral(l, k, lp)
                if abs(g) > 1e-15:
                    sel = "allowed" if (l + k + lp) % 2 == 0 else "forbidden"
                    print(f"{l:>3} {lp:>7} {k:>3} {g:>12.6f} {sel:>12}")

    print("""
Key observations:
- k=0 multipole couples all l=l' pairs (monopole term)
- k=2 couples l↔l±2, l↔l (quadrupole)
- Odd-k multipoles vanish for L=0 (parity)
- The coupling is SPARSE: each l couples to at most ~2l neighbors
- This sparsity determines the graph structure of the hyperangular lattice
""")


# =============================================================================
# 7. CONNECTION TO PAPER 0: S³ PACKING
# =============================================================================

def print_packing_connection() -> None:
    """Analyze the connection between S⁵ hyperspherical harmonics and S³ packing."""
    print("\n" + "=" * 70)
    print("CONNECTION TO PAPER 0: S³ vs S⁵ PACKING")
    print("=" * 70)

    print("""
PAPER 0 STRUCTURE:
  S³ (one electron): states packed as shells of 2n² states
  Eigenvalues: λ_n = -(n² - 1) on unit S³
  Kinetic scale: κ = -1/16 maps graph Laplacian to physics

S⁵ STRUCTURE:
  S⁵ (two electrons): states packed as shells of (K+1)(K+2)²(K+3)/12
  Eigenvalues: Λ²_K = -K(K+4) on unit S⁵ (Gegenbauer C²_K)
  Kinetic scale: UNKNOWN — this is a key question for Phase 2

DEGENERACY COMPARISON:
""")

    print(f"  {'Shell':>6} {'S³ (2n²)':>10} {'S⁵ full':>10} {'S⁵ L=0':>10} {'S⁵/S³':>8}")
    print("  " + "-" * 48)
    for K in range(10):
        n = K + 1
        g3 = 2 * n**2
        g5 = degeneracy_S5(K)
        g5L0 = degeneracy_S5_L0(K)
        ratio = g5 / g3 if g3 > 0 else 0
        print(f"  {K:>6} {g3:>10} {g5:>10} {g5L0:>10} {ratio:>8.1f}")

    print("""
KEY QUESTION: Does the S⁵ hyperangular eigenvalue spectrum have an SO(6)
"accidental" degeneracy analogous to SO(4) for hydrogen?

ANSWER: NO — the Coulomb potential breaks SO(6) immediately.

For HYDROGEN on S³:
  V = -Z/r → after Fock projection becomes a constant on S³
  → SO(4) symmetry is EXACT (all l states at fixed n are degenerate)
  → This is why the S³ lattice works so perfectly

For HELIUM on S⁵:
  V = -Z/r₁ - Z/r₂ + 1/r₁₂ → NOT a constant on S⁵
  → SO(6) symmetry is BROKEN by the potential
  → K is only an APPROXIMATE quantum number
  → The adiabatic channel index ν replaces K as the practical label

IMPLICATION FOR THE LATTICE:
  The hyperspherical lattice CANNOT be a pure "packing lattice" like S³.
  It must be a COMPUTATIONAL lattice (like prolate spheroidal) where the
  angular structure is solved numerically, not given by group theory.

  However, the L=0 reduction means the effective problem is LOW-DIMENSIONAL:
  - S³ atomic: 0D angular (exact) + 0D radial (eigenvalue gives energy)
  - Prolate molecular: 1D angular (η) + 1D radial (ξ) = 2D total
  - Hyperspherical He: 1D angular (α, with l-coupling) + 1D radial (R) = 2D total

  The DIMENSIONALITY is the same as prolate! The difference is the
  inter-channel coupling, which adds off-diagonal structure.
""")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    print_symmetry_analysis()
    print_degeneracy_table()
    print_irrep_decomposition()
    print_ladder_operators()
    analyze_charge_function(Z=2.0)
    print_coupling_matrix()
    print_packing_connection()

    print("\n" + "=" * 70)
    print("SUMMARY: ALGEBRAIC STRUCTURE OF THE HYPERSPHERICAL LATTICE")
    print("=" * 70)
    print("""
1. SYMMETRY: SO(6) on S⁵, broken to SO(3)_L × S₂ × Z₂ by Coulomb potential.
   For ¹S He: L=0, S=0, π=+1 reduces problem to 3 coordinates (R, α, θ₁₂).

2. QUANTUM NUMBERS: (K, l₁, l₂, L, M) → for L=0: (K, l) with K even, l₁=l₂=l.
   Channel index ν labels adiabatic eigenstates at each R.

3. LADDER OPERATORS: SO(6) has 15 generators. Nine "boosts" connect different
   (l₁, l₂) sectors. These define the off-diagonal angular lattice edges.
   However, the practical lattice uses numerical channel functions, not
   algebraic ladder operators.

4. GRAPH STRUCTURE: Nodes = (R_i, ν). Edges = FD in R (within channel) +
   non-adiabatic coupling P_μν (between channels). This is a COUPLED CHANNEL
   GRAPH — structurally the same kind of object as S³ and prolate lattices,
   but with R-dependent inter-channel edges (fiber bundle structure).

5. CONNECTION TO PAPER 0: S⁵ degeneracy follows (K+1)(K+2)²(K+3)/12 packing,
   analogous to 2n² on S³. But the Coulomb potential breaks SO(6), so the
   lattice is computational (like prolate) rather than algebraic (like S³).

6. EFFECTIVE DIMENSIONALITY: 1D angular (coupled α channels) + 1D radial (R)
   = 2D total. Same as prolate spheroidal. This is encouraging for feasibility.
""")
