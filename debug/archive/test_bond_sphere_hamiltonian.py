"""
Bond Sphere Minimal Symbolic Hamiltonian — Diagnostic
=====================================================

Can the equilibrium bond angle γ_eq of LiH be derived analytically from a
fully consistent S³ Hamiltonian?  If yes, R_eq follows from the γ(R) map.

Model: 2 spatial orbitals (φ₁ = Li 1s, n=1; φ₂ = bonding 2σ, n=2).
       4 electrons, closed-shell.

H1[i,j] = (-p₀²/2) S[i,j] + (1 - β_j) V[i,j]
S[i,j]  = D^{n_i}_{00,00}(γ)   (overlap via SO(4) D-matrix)
V[i,j]  = (p₀/n_j) D^{n_i}_{00,00}(γ)   (Sturmian nuclear attraction identity)

Analytic D-matrix elements (verified numerically):
    D^1_{00,00}(γ) = 1
    D^2_{00,00}(γ) = cos(γ)

Status: Pure diagnostic — no production code changes.
"""

import numpy as np
from scipy.optimize import brentq, minimize_scalar
import sympy as sp
import sys
import os

# ──────────────────────────────────────────────────────────────────────
# Step 0: Verify D-matrix forms (already confirmed, include as sanity)
# ──────────────────────────────────────────────────────────────────────

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from geovac.wigner_so4 import wigner_D_so4

print("=" * 70)
print("BOND SPHERE MINIMAL HAMILTONIAN — DIAGNOSTIC")
print("=" * 70)

print("\n--- Step 0: D-matrix verification ---")
for gval in [0.5, 1.0, 2.0]:
    d1 = wigner_D_so4(1, 0, 0, 0, 0, gval)
    d2 = wigner_D_so4(2, 0, 0, 0, 0, gval)
    print(f"  γ={gval:.1f}  D^1_00={d1:.6f} (expect 1)  "
          f"D^2_00={d2:.6f} (expect cos={np.cos(gval):.6f})")


# ──────────────────────────────────────────────────────────────────────
# Step 1: Symbolic SymPy setup
# ──────────────────────────────────────────────────────────────────────

print("\n--- Step 1: Symbolic 2×2 Hamiltonian ---")

gamma_s, p0_s = sp.symbols('gamma p0', positive=True, real=True)
beta1_s, beta2_s = sp.symbols('beta1 beta2', positive=True, real=True)
ZA_s, ZB_s = sp.symbols('Z_A Z_B', positive=True, integer=True)

# D-matrix elements (s-states only, verified):
#   D^1_{00,00}(γ) = 1
#   D^2_{00,00}(γ) = cos(γ)
# For the 2-state model we have orbitals with n₁=1, n₂=2.
# The overlap and nuclear attraction matrices involve cross-shell elements.
#
# IMPORTANT: The D-matrix D^n acts within the n-shell. For cross-shell
# (n₁≠n₂) overlaps, the Sturmian overlap is different. In the Sturmian
# basis with shared p₀, the overlap matrix is:
#   S[i,j] = δ_{n_i,n_j}  for same-atom orbitals (orthonormal within shells)
# Cross-atom overlaps between different n-shells use the Shibuya-Wulfman
# overlap integral, which for s-states is:
#   S_{n₁s,n₂s}(γ) = <n₁,0,0| R(γ) |n₂,0,0>
# This is NOT simply a D-matrix element (D-matrices act within a shell).
#
# For the minimal model, we need the SW cross-shell overlap.
# From Shibuya-Wulfman theory (Paper 8 Sec IV), the cross-shell overlap
# between hydrogenic s-states in the Sturmian basis is:
#   <1s|R(γ)|2s> = -sin(γ)  [from the SW connection coefficients]
#
# Let me verify this numerically first.

print("\n  Checking cross-shell D-matrix elements...")
print("  NOTE: D^n is intra-shell. Cross-shell requires SW overlap.")

# For intra-shell: D^n_{l'm',lm}(γ) is well-defined
# For cross-shell: need to compute <n'l'm'|exp(iγA_y)|nlm> directly
# This requires the full SO(4) rotation in the combined Hilbert space.

# The SW overlap for s-states can be extracted from the hyperspherical
# harmonic representation. On S³, the states are Y_{nlm}(χ,θ,φ) and
# the rotation shifts χ → χ + γ (bond rotation along the internuclear axis).
#
# For s-states (l=m=0), Y_{n,0,0}(χ) ∝ sin(nχ)/(sin χ),
# so the overlap integral becomes:
#   <n'|R(γ)|n> = ∫ Y*_{n',0,0}(χ) Y_{n,0,0}(χ-γ) dΩ₃
#
# For n'=1, n=2 (s-states), using the hyperspherical harmonics:
#   Y_{1,0,0} = 1/√(2π²)  (constant on S³, the n=1 state)
#   Y_{2,0,0} = (2/√(2π²)) cos(χ)  (n=2, l=0 on S³)
#
# Wait — Y_{n,0,0}(χ) on S³ involves the Gegenbauer polynomial.
# Specifically: Y_{n,0,0}(χ) = C^1_{n-1}(cos χ) × normalization
# C^1_0(x) = 1,  C^1_1(x) = 2x
#
# So the overlap is an integral of Gegenbauer polynomials with shifted argument.
# This is nontrivial. Let me compute it numerically using the full SO(4) machinery.

# Numerical cross-shell overlap via hyperspherical harmonics
def hyperspherical_overlap_ss(n1: int, n2: int, gamma_val: float,
                               n_quad: int = 500) -> float:
    """Compute <n1,0,0|R(γ)|n2,0,0> on S³ via numerical quadrature.

    The hyperspherical harmonics for s-states (l=m=0) on S³ are:
        Y_{n,0,0}(χ) = sqrt(2/π) * sin(n*χ) / sin(χ)

    The bond rotation shifts χ → χ - γ (rotation about A_y).
    For the s-channel (l=m=0), this is equivalent to shifting
    the hyperspherical angle.

    The overlap is:
        <n1|R(γ)|n2> = (2/π) ∫₀^π sin(n1*χ)/sin(χ) * sin(n2*(χ-γ))/sin(χ-γ)
                        × sin²(χ) dχ × correction
    Wait — the volume element on S³ is sin²(χ) sin(θ) dχ dθ dφ, and for s-states
    the θ,φ integral gives 4π. So:

        <n1,0,0|R(γ)|n2,0,0> = (2/π) ∫₀^π sin(n1*χ) sin(n2*(χ')) sin²(χ)/sin(χ)sin(χ') ...

    Actually this gets complicated with the rotation. Let me use the matrix
    representation directly — express R(γ) in the |n,l,m⟩ basis and read off
    the cross-shell elements.
    """
    # Use the full angular momentum approach.
    # On S³, states |n,l,m⟩ form a complete basis.
    # The rotation exp(iγ A_y) mixes states with same l,m across different n.
    # But A_y is a generator of SO(4), and in the |n,l,m⟩ basis it has
    # matrix elements between different n-shells.
    #
    # Actually, A_y connects n to n±1 (Runge-Lenz ladder property).
    # So <1s|exp(iγ A_y)|2s> involves the A_y matrix element between n=1 and n=2.
    #
    # The Runge-Lenz vector A_y has matrix elements:
    #   <n',l',m'|A_y|n,l,m> ≠ 0 only if n' = n (within-shell)
    #
    # Wait — A is a GENERATOR of SO(4) that acts WITHIN the n-shell.
    # The n-shell is an irrep of SO(4). Different n-shells are different irreps.
    # exp(iγ A_y) does NOT mix different n-shells!
    #
    # This means the cross-shell overlap <1s|R(γ)|2s> = 0 identically
    # if R(γ) is the SO(4) bond rotation.
    return 0.0  # placeholder — needs careful treatment


# KEY INSIGHT: The SO(4) rotation exp(iγ A_y) acts within each n-shell
# (irreducible representation). It CANNOT produce cross-shell overlaps.
# This means in the pure bond-sphere picture, the overlap matrix is
# BLOCK-DIAGONAL in n:
#   S[1s,1s] = 1,  S[2s,2s] = 1,  S[1s,2s] = S[2s,1s] = 0
#
# The nuclear attraction matrix has the same structure (Sturmian identity):
#   V[i,j] = (p₀/n_j) S[i,j]
# So V is also block-diagonal:
#   V[1s,1s] = p₀/1 = p₀,  V[2s,2s] = p₀/2,  V[1s,2s] = V[2s,1s] = 0
#
# This means the H1 matrix is DIAGONAL in this 2-state model!
# H1[1,1] = -p₀²/2 + (1-β₁)p₀    (n=1 orbital)
# H1[2,2] = -p₀²/2 + (1-β₂)p₀/2  (n=2 orbital)
# H1[1,2] = H1[2,1] = 0
#
# Wait — but this can't be right for a diatomic. The cross-atom coupling
# must come from somewhere. Let me reconsider.
#
# In the Sturmian basis for a DIATOMIC, each orbital is centered on a
# different atom. The rotation R(γ) maps atom A → atom B. So the cross-atom
# overlap is:
#   <φ_A|R(γ)|φ_B> where φ_A is on atom A and φ_B is on atom B
# But if φ_A has n=1 and φ_B has n=2, these are in different SO(4) irreps,
# and the rotation cannot connect them.
#
# However, in the Shibuya-Wulfman formulation, the overlap between
# Sturmian functions on DIFFERENT centers is computed differently.
# The SW overlap integral is NOT an SO(4) rotation matrix element.
# It uses the momentum-space convolution / Fourier transform approach.
#
# From Paper 8, the SW overlap between Sturmians is:
#   S^{SW}_{n'l'm', nlm}(R, p₀) = <n'l'm'; A | nlm; B>
# This is a genuine spatial overlap, not an SO(4) D-matrix element.
#
# The D-matrix gives the SAME-CENTER intra-shell rotation.
# The cross-center overlap is the SW integral.

print("\n  CRITICAL INSIGHT: SO(4) D-matrix is intra-shell (within irrep).")
print("  Cross-shell overlaps require Shibuya-Wulfman (SW) overlap integrals.")
print("  These are NOT D-matrix elements.")
print()
print("  Revising model: use SW overlaps for cross-center coupling.")

# ──────────────────────────────────────────────────────────────────────
# Revised Model: Sturmian H1 with SW cross-center overlaps
# ──────────────────────────────────────────────────────────────────────
#
# The correct 2×2 Sturmian Hamiltonian for a diatomic AB with orbitals
# φ₁ (n=1, on A) and φ₂ (n=2, on B) at shared momentum p₀:
#
# DIAGONAL (Sturmian construction — all orbitals have same kinetic energy):
#   H1[1,1] = -p₀²/2 + Z_A·p₀/n₁ = -p₀²/2 + Z_A·p₀
#   H1[2,2] = -p₀²/2 + Z_B·p₀/n₂ = -p₀²/2 + Z_B·p₀/2
#
# Wait — in the Sturmian basis the diagonal is UNIFORM:
#   ε_k = -p₀²/2 for ALL orbitals (this is the Sturmian hallmark).
# The nuclear attraction is:
#   <k|V_A|k> = Z_A·p₀/n_k (self-nuclear)  +  Z_B·p₀·D^{n_k}(γ)/n_k' (cross-nuclear via D-matrix... but only intra-shell!)
#
# Actually, let me re-derive this carefully from the Sturmian secular equation.
# The Sturmian equation for a diatomic with nuclei A (Z_A) and B (Z_B) is:
#   [-p₀²/2 - Z_A/r_A - Z_B/r_B] φ_k = 0  (if β_k = 1)
# More generally:
#   [-p₀²/2 - β_k(Z_A/r_A + Z_B/r_B)] φ_k = 0
# where β_k is the Sturmian eigenvalue (β_k ≈ Z_eff_k/Z_total or similar).
#
# In the prolate spheroidal basis, the Sturmian equation separates and β_k
# are the eigenvalues of the angular equation. The full Hamiltonian matrix
# in the Sturmian basis is:
#   H[i,j] = -p₀²/2 · S[i,j] + β_j · V[i,j]    ... no wait
#
# Let me use the standard Sturmian secular equation formulation.
# In the Sturmian approach (Avery & Avery), the secular equation is:
#   Σ_j [T_{ij} + V^0_{ij} - ε δ_{ij}] c_j = 0
# where T is kinetic, V^0 is the nuclear potential, and ε = -p₀²/2.
# But in the Sturmian basis, T_{ij} + V^0_{ij} = β_j · V^0_{ij} / β_j = V^0_{ij}
# ... this is getting circular.
#
# Let me use the formulation from the user's prompt directly:
#   H1[i,j] = (-p₀²/2) · S[i,j] + (1 - β_j) · V[i,j]
# where S is the overlap matrix and V is the nuclear attraction matrix.
#
# For the Sturmian basis, the eigenstates satisfy:
#   (-p₀²/2 - β_k · V) |k⟩ = 0
# Rearranging: (-p₀²/2) |k⟩ = β_k · V |k⟩
# So in the basis of Sturmian functions, acting on state j:
#   H1|j⟩ = (-p₀²/2)|j⟩ + V|j⟩ = (-p₀²/2)|j⟩ + V|j⟩
#          = β_j · V|j⟩ + V|j⟩ - β_j · V|j⟩   ... hmm
# Actually: H1 = T + V = (-p₀²/2)S + V. And since T|j⟩ = β_j V|j⟩ in Sturmian:
#   V|j⟩ = (1/β_j)T|j⟩ = (-p₀²/(2β_j))S|j⟩
# Wait, this only holds for the eigenvalue equation, not as an operator identity.
#
# Let me just use the matrix elements as given in the prompt:
#   H1[i,j] = (-p₀²/2) · S[i,j] + (1 - β_j) · V[i,j]
#
# Now I need S[i,j] and V[i,j] for the 2×2 case.
# These are the Sturmian overlap and nuclear attraction matrices
# between φ₁ (n=1, atom A) and φ₂ (n=2, atom B).
#
# For Sturmians with shared p₀ on the SAME center, S[i,j] = δ_{ij}.
# For Sturmians on DIFFERENT centers separated by R (or bond angle γ):
#   S[1,2] = <1s; A | 2s; B> = Shibuya-Wulfman overlap
#
# The SW overlap for Coulomb Sturmians in momentum space is known analytically.
# For s-states: S^{SW}_{1s,2s}(p₀R) has a known closed form.
# From the momentum-space overlap (Avery, "Hyperspherical Harmonics"):
#
# The momentum-space Sturmian is:
#   φ̃_n(p) = N_n · p₀^{5/2} / (p² + p₀²)^2 · C^2_{n-1}(cos χ)
# where cos χ = (p² - p₀²)/(p² + p₀²) is the hyperspherical angle.
#
# For s-states specifically:
#   φ̃_{1s}(p) ∝ 1/(p² + p₀²)²   (C^2_0 = 1)
#   φ̃_{2s}(p) ∝ (p² - p₀²)/(p² + p₀²)³   (C^2_1(x) = 4x... let me check)
#
# Actually C^λ_n are Gegenbauer polynomials: C^2_0(x) = 1, C^2_1(x) = 4x.
#
# The cross-center overlap involves a phase factor e^{ip·R}:
#   S[1s_A, 2s_B] = ∫ φ̃*_{1s}(p) · e^{ip·R} · φ̃_{2s}(p) d³p
#
# For the bond-sphere model, this phase factor is the bond rotation.
# After Fock projection to S³, e^{ip·R} → rotation by γ on S³.
# But the rotation acts on the FULL Hilbert space (all n), not within a shell.
#
# So the cross-shell overlap IS related to the bond rotation, but computed
# as a FULL rotation matrix element (not the intra-shell D-matrix).
#
# The full S³ rotation matrix in the |n,l,m⟩ basis mixes different n-shells
# because the rotation exp(iγ·Â_y) where  is not an SO(4) generator of a
# single irrep — it's a GLOBAL rotation of S³.
#
# ACTUALLY: let me reconsider. In the bond-sphere picture, the bond rotation
# is exp(iγ A_y) where A is the Runge-Lenz vector. A acts WITHIN each n-shell
# as a generator of SO(4). It does NOT connect different n-shells.
#
# But the cross-center Sturmian overlap is a SPATIAL overlap, not an SO(4)
# rotation. The spatial overlap between Sturmians on different centers is:
#   S[a,b] = ∫ φ*_a(r-R_A) φ_b(r-R_B) d³r
# In momentum space with the Fock projection, this becomes:
#   S[a,b] = ∫ φ̃*_a(p) e^{i p·(R_B-R_A)} φ̃_b(p) d³p
# The plane wave e^{ip·R} on S³ is a TRANSLATION, not a rotation.
# On flat space, translation ≠ rotation. On S³, a translation corresponds
# to a specific conformal transformation, not a simple SO(4) rotation.
#
# However, Shibuya & Wulfman showed that after Fock projection, the
# translation operator becomes expressible in terms of SO(4) matrices.
# Specifically, the momentum shift acts as:
#   <n'l'm'|e^{ip·R}|nlm> = Σ_n'' B_{n'n''}(p₀R) D^{n''}(γ)_{...}
# This involves an INFINITE sum over intermediate shells — it's not a simple
# D-matrix element.
#
# For the MINIMAL model, I'll compute the SW overlap numerically.
# The exact formula for s-state Sturmian overlaps is known from
# Avery & Avery, J. Math. Chem. 46, 164 (2009):

def sw_overlap_ss_numerical(n1: int, n2: int, p0_val: float, R_val: float,
                             n_quad: int = 2000) -> float:
    """Compute SW overlap <n1,s|n2,s> between Sturmians on different centers.

    Uses numerical integration in momentum space:
        S = ∫ φ̃*_{n1}(p) e^{ip·R} φ̃_{n2}(p) d³p / (2π)³

    For s-states, the angular integral over p̂ gives a sinc-like factor.
    After angular integration:
        S = (4π/(2π)³) ∫₀^∞ φ̃_{n1}(p) φ̃_{n2}(p) sin(pR)/(pR) p² dp

    The Sturmian momentum-space wavefunctions (Fock, 1935) are:
        φ̃_{n,0,0}(p) = N_n · p₀^{5/2} · C^1_{n-1}(w) · 4p₀² / (p²+p₀²)²
    where w = (p²-p₀²)/(p²+p₀²) and C^1 are Gegenbauer polynomials.

    Wait — the normalization and Gegenbauer index depend on convention.
    The standard Fock Sturmian in 3D momentum space is:
        Φ_{nlm}(p) = M_n(p) Y_{lm}(p̂)
    with the radial part:
        M_n(p) = N_n (2p₀)^{5/2} p^l / (p²+p₀²)^{l+2} C^{l+1}_{n-l-1}(w)

    For l=0:
        M_n(p) = N_n (2p₀)^{5/2} / (p²+p₀²)^2 C^1_{n-1}(w)

    C^1_0(w) = 1
    C^1_1(w) = 2w

    Normalization: ∫₀^∞ |M_n(p)|² p² dp = 1/(4π)  (to give normalized φ_nlm)
    Actually let's just use the explicit hydrogen momentum wavefunctions.
    """
    # Use explicit hydrogen momentum-space wavefunctions (Podolsky & Pauling, 1929)
    # For 1s: φ̃_{1s}(p) = (2√2/π) p₀^{5/2} / (p²+p₀²)²
    # For 2s: φ̃_{2s}(p) = (2√2/π) (2p₀)^{5/2} (p²-p₀²) / (p²+p₀²)³  ... hmm
    #
    # Actually the standard momentum-space hydrogen wavefunctions for Sturmians
    # with parameter p₀ (not necessarily = Z/n) are:
    # φ̃^{p₀}_{n,0,0}(p) = N · p₀^2 / (p²+p₀²)^2 · U_{n-1}(w)
    # where U_{n-1} is a Chebyshev polynomial of the 2nd kind (≡ C^1_{n-1})
    # and w = (p²-p₀²)/(p²+p₀²).
    #
    # U_0(w) = 1,  U_1(w) = 2w
    #
    # The normalization for Sturmian overlaps uses:
    # ∫ φ̃*_a φ̃_b d³p / (2π)³ = (Sturmian overlap)
    #
    # For the overlap between functions on centers separated by R:
    # S_{ab}(R) = ∫ φ̃*_a(p) φ̃_b(p) e^{ipR cos θ_p} d³p / (2π)³
    # After angular integration (s-states → l=0):
    # S_{ab}(R) = (1/(2π²)) ∫₀^∞ M_a(p) M_b(p) sin(pR)/(pR) p² dp

    # Sturmian radial functions in p-space (unnormalized, s-states):
    # M_n(p, p₀) = p₀² / (p²+p₀²)² · U_{n-1}((p²-p₀²)/(p²+p₀²))
    #            × 2^{5/2} ... normalization TBD

    # Let me use a different approach: compute in position space directly.
    # Sturmian functions in position space with parameter p₀:
    #   S_{n,0}(r; p₀) = N_n · (2p₀r)^1 · e^{-p₀r} · L^1_{n-1}(2p₀r)
    # where L^1 are associated Laguerre polynomials.
    # N_n = (2p₀)^{3/2} √((n-1)!/(n+1)!)... actually for Sturmians:
    #   S_{n,l}(r; p₀) = N_{nl} (2p₀r)^l e^{-p₀r} L^{2l+1}_{n-l-1}(2p₀r)
    #   with normalization: ∫₀^∞ S_{n,l} S_{n',l} r² dr = δ_{nn'}/p₀
    # (Sturmians are orthogonal with weight 1/r, not 1).

    # For s-states (l=0):
    #   S_{1,0}(r) = N₁ · e^{-p₀r} · L^1_0(2p₀r) = N₁ e^{-p₀r}
    #   S_{2,0}(r) = N₂ · e^{-p₀r} · L^1_1(2p₀r) = N₂ e^{-p₀r} (2-2p₀r)

    # L^1_0(x) = 1,  L^1_1(x) = 2-x

    # Normalization with Sturmian inner product (weight 1/r):
    #   ∫₀^∞ S_{n,0}(r)² r dr = 1/p₀  → N_n² ∫ ... = 1/p₀

    # For the standard overlap (weight 1, not 1/r):
    #   <1s_A|2s_B> = ∫ S*_{1,0}(|r-R_A|) S_{2,0}(|r-R_B|) d³r
    # This is a standard two-center overlap integral with STO-like functions.
    # For s-states this has a known analytic form.

    # Let me just compute it numerically on a grid.
    from scipy.integrate import quad
    from scipy.special import eval_genlaguerre

    # Sturmian radial functions (normalized to ∫|S|²r²dr = 1):
    def sturmian_radial(n_val, r, p0):
        """Radial Sturmian S_{n,0}(r; p₀), normalized: ∫|S|²r²dr = 1."""
        x = 2.0 * p0 * r
        L = eval_genlaguerre(n_val - 1, 1, x)  # L^1_{n-1}(x)
        unnorm = np.exp(-p0 * r) * L
        # Normalization: ∫₀^∞ (e^{-p₀r} L^1_{n-1}(2p₀r))² r² dr
        # Compute numerically
        return unnorm

    # Compute normalization constants
    def norm_const(n_val, p0):
        integrand = lambda r: sturmian_radial(n_val, r, p0)**2 * r**2
        val, _ = quad(integrand, 0, 50.0/p0, limit=200)
        return 1.0 / np.sqrt(val)

    N1 = norm_const(n1, p0_val)
    N2 = norm_const(n2, p0_val)

    # Two-center overlap for s-states (spherically symmetric):
    # <S_{n1}(r_A) | S_{n2}(r_B)> where r_A = |r-R_A|, r_B = |r-R_B|, |R_A-R_B|=R
    # Place A at origin, B at (0,0,R).
    # In spherical coords centered at A: r_B = sqrt(r² + R² - 2rR cos θ)
    # Integrate: ∫₀^∞ ∫₀^π S_{n1}(r) S_{n2}(r_B) r² sin θ dθ dr × 2π

    def integrand_2d(r, costh):
        rB = np.sqrt(r**2 + R_val**2 - 2*r*R_val*costh)
        if rB < 1e-12:
            rB = 1e-12
        f1 = N1 * sturmian_radial(n1, r, p0_val)
        f2 = N2 * sturmian_radial(n2, rB, p0_val)
        return f1 * f2 * r**2

    # Double integral: ∫₀^∞ dr ∫_{-1}^{1} d(cosθ) × 2π × integrand
    from scipy.integrate import dblquad
    r_max = 40.0 / p0_val
    result_val, err = dblquad(
        integrand_2d,
        -1.0, 1.0,   # cos θ limits
        0.0, r_max,   # r limits
        epsabs=1e-8, epsrel=1e-8
    )
    return 2.0 * np.pi * result_val


# ──────────────────────────────────────────────────────────────────────
# Compute SW overlaps for the minimal model
# ──────────────────────────────────────────────────────────────────────

print("\n--- Computing Shibuya-Wulfman cross-center overlaps ---")

p0_test = 1.0
R_test = 3.015

# Same-center overlaps (orthonormal)
S11 = 1.0  # <1s_A|1s_A> (same center, same orbital)
S22 = 1.0  # <2s_B|2s_B> (same center, same orbital)

# Cross-center overlaps
print(f"  Computing S[1s_A, 2s_B] at p₀={p0_test}, R={R_test} ...")
S12_num = sw_overlap_ss_numerical(1, 2, p0_test, R_test)
print(f"  S[1s_A, 2s_B] = {S12_num:.6f}")

# Also compute S[1s_A, 1s_B] and S[2s_A, 2s_B] for reference
print(f"  Computing S[1s_A, 1s_B] at p₀={p0_test}, R={R_test} ...")
S11_cross = sw_overlap_ss_numerical(1, 1, p0_test, R_test)
print(f"  S[1s_A, 1s_B] = {S11_cross:.6f}")

print(f"  Computing S[2s_A, 2s_B] at p₀={p0_test}, R={R_test} ...")
S22_cross = sw_overlap_ss_numerical(2, 2, p0_test, R_test)
print(f"  S[2s_A, 2s_B] = {S22_cross:.6f}")


# ──────────────────────────────────────────────────────────────────────
# Step 1 (revised): Build 2×2 Hamiltonian numerically
# ──────────────────────────────────────────────────────────────────────

print("\n--- Step 1 (revised): Numerical 2×2 Hamiltonian ---")

# Model parameters
ZA = 3  # Li
ZB = 1  # H
beta1_val = 1.033  # prolate spheroidal eigenvalue for Li 1s
beta2_val = 1.950  # prolate spheroidal eigenvalue for bonding 2σ

def build_H1_numerical(p0_val: float, R_val: float) -> tuple:
    """Build the 2×2 Sturmian Hamiltonian and overlap matrices.

    Orbitals: φ₁ = 1s on Li (n=1), φ₂ = 2σ on H (n=2).

    H1[i,j] = (-p₀²/2) S[i,j] + (1-β_j) V[i,j]

    For the SELF-nuclear terms:
        V[i,i] = (Z_center · p₀/n_i)   (Sturmian identity at own center)
    For CROSS-nuclear terms:
        V[i,j] = overlap × Z_other × p₀/n_j   (SW overlap weighted)

    Actually, let me use the proper Sturmian matrix elements.
    In the Sturmian basis, the potential matrix element is:
        V[i,j] = <i| -(Z_A/r_A + Z_B/r_B) |j>
    For Sturmians centered on specific atoms, this decomposes into
    self-nuclear and cross-nuclear contributions.

    For the minimal model:
        V[1,1] = <1s_A| -Z_A/r_A |1s_A> + <1s_A| -Z_B/r_B |1s_A>
        V[1,2] = <1s_A| -Z_A/r_A |2s_B> + <1s_A| -Z_B/r_B |2s_B>
    etc.

    The self-nuclear integral for a Sturmian: <ns_X| -Z_X/r_X |ns_X> = -Z_X p₀/n
    The cross-nuclear integral: <ns_A| -Z_B/r_B |ns_A> = ?? (multicenter integral)

    For simplicity in this minimal model, let's use the Sturmian Hamiltonian
    as given in the prompt:
        H1[i,j] = (-p₀²/2) S[i,j] + (1-β_j) V_nuc[i,j]
    where V_nuc[i,j] = <i|V_total|j> is the total nuclear attraction.

    The Sturmian secular equation gives: for eigenstates,
        (-p₀²/2) S|k⟩ = -β_k V_total |k⟩
    So: V_total |k⟩ = (p₀²/(2β_k)) S |k⟩
    This means: V_total[i,k] = (p₀²/(2β_k)) S[i,k]

    Therefore:
        H1[i,j] = (-p₀²/2) S[i,j] + (1-β_j)(p₀²/(2β_j)) S[i,j]
                 = (-p₀²/2) S[i,j] + (p₀²/2)(1/β_j - 1) S[i,j]
                 = (p₀²/2) (-1 + 1/β_j - 1) S[i,j]

    Wait, that gives: H1[i,j] = (p₀²/2)(1/β_j - 2) S[i,j]
    Hmm, that doesn't look right either. Let me re-derive.

    The Sturmian equation: (T + β_k V)|k⟩ = 0, where T = -p₀²/2 S (kinetic).
    So: -p₀²/2 S|k⟩ + β_k V|k⟩ = 0
    → V|k⟩ = (p₀²/(2β_k)) S|k⟩

    The PHYSICAL Hamiltonian: H = T + V = -p₀²/2 S + V
    Matrix element: H[i,j] = -p₀²/2 S[i,j] + V[i,j]

    But V[i,j] = Σ_k V[i,k]δ[k,j] only if j is a Sturmian eigenstate.
    Since we're working in the Sturmian basis where V|j⟩ = (p₀²/(2β_j))S|j⟩:
        V[i,j] = <i|V|j⟩ = (p₀²/(2β_j)) S[i,j]

    Therefore:
        H[i,j] = -p₀²/2 S[i,j] + (p₀²/(2β_j)) S[i,j]
                = (p₀²/2)(1/β_j - 1) S[i,j]

    This is proportional to the overlap matrix! The energy eigenvalues
    come from the generalized eigenvalue problem H c = E S c, which
    reduces to:
        (p₀²/2)(1/β_j - 1) S c = E S c
    If S is invertible (which it is for non-zero overlap), this gives:
        E_j = (p₀²/2)(1/β_j - 1)

    For β₁ = 1.033: E₁ = (p₀²/2)(1/1.033 - 1) = (p₀²/2)(-0.03195) = -0.01597 p₀²
    For β₂ = 1.950: E₂ = (p₀²/2)(1/1.950 - 1) = (p₀²/2)(-0.48718) = -0.24359 p₀²

    Wait — this means the eigenvalues are INDEPENDENT of γ (and hence R)!
    That can't be right. The eigenvalues depend only on β_j and p₀, not on
    the bond geometry. This makes sense because the Sturmian eigenvalues β_k
    already encode the nuclear geometry — they are the prolate spheroidal
    eigenvalues at a given R.

    So in a CONSISTENT Sturmian treatment, the β values depend on R:
        β_k = β_k(R)
    And the total energy is:
        E_total(R) = Σ_k n_k · E_k(R) + V_NN(R) + V_ee
    where n_k is the occupation number.

    For the total energy to have an R minimum, we need β(R).
    The β(R) dependence is what gives the potential energy surface.

    This is a different model than what I initially set up. The question
    becomes: can we express β(R) analytically, and then find R_eq from
    dE_total/dR = 0?
    """
    # Compute overlap matrix using SW integrals
    S12_val = sw_overlap_ss_numerical(1, 2, p0_val, R_val)
    S_mat = np.array([[1.0, S12_val],
                      [S12_val, 1.0]])

    # H1[i,j] = (p₀²/2)(1/β_j - 1) S[i,j]
    H1_mat = np.zeros((2, 2))
    betas = [beta1_val, beta2_val]
    for i in range(2):
        for j in range(2):
            H1_mat[i, j] = (p0_val**2 / 2.0) * (1.0/betas[j] - 1.0) * S_mat[i, j]

    return H1_mat, S_mat


# ──────────────────────────────────────────────────────────────────────
# KEY THEORETICAL FINDING
# ──────────────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("KEY THEORETICAL FINDING")
print("=" * 70)
print("""
In the Sturmian basis, the physical Hamiltonian matrix element is:
    H[i,j] = (p₀²/2)(1/β_j - 1) · S[i,j]

where S[i,j] is the overlap matrix and β_j are the Sturmian eigenvalues.

The H matrix is PROPORTIONAL to the overlap matrix S. Therefore the
generalized eigenvalue problem H c = E S c reduces to:
    E_k = (p₀²/2)(1/β_k - 1)

The eigenvalues are INDEPENDENT of the bond geometry (γ or R) when β
values are held fixed. The R-dependence enters ONLY through β_k(R).

This means the minimal Bond Sphere model with FIXED β values cannot
predict R_eq from a γ-optimization. The β(R) dependence IS the PES.
""")

# ──────────────────────────────────────────────────────────────────────
# Revised approach: Use fixed β but include V_NN and the p₀-R coupling
# ──────────────────────────────────────────────────────────────────────

print("=" * 70)
print("REVISED APPROACH: Total energy with V_NN")
print("=" * 70)
print("""
Even with fixed β, the total energy includes nuclear repulsion:
    E_total(R, p₀) = Σ_k n_k E_k(p₀) + V_NN(R) + V_ee

    = Σ_k n_k (p₀²/2)(1/β_k - 1) + Z_A Z_B / R + V_ee

The self-consistency condition p₀ = √(-2E_total) couples p₀ and R.
And R enters through the bond angle γ = 2 arctan(1/(p₀R)).
So there IS an implicit R-dependence via the p₀(R) coupling.

Let me solve this self-consistent system.
""")


def total_energy_1e(p0_val: float, R_val: float) -> float:
    """Total 1-electron energy for the 2-orbital, 4-electron model.

    E_total = 2*E₁ + 2*E₂ + V_NN
    (closed shell: 2 electrons in each orbital)

    E_k = (p₀²/2)(1/β_k - 1)
    V_NN = Z_A * Z_B / R
    """
    E1 = (p0_val**2 / 2.0) * (1.0/beta1_val - 1.0)
    E2 = (p0_val**2 / 2.0) * (1.0/beta2_val - 1.0)
    V_NN = ZA * ZB / R_val
    return 2.0 * E1 + 2.0 * E2 + V_NN


# ──────────────────────────────────────────────────────────────────────
# Step 2: Symbolic analysis
# ──────────────────────────────────────────────────────────────────────

print("\n--- Step 2: Symbolic eigenvalue analysis ---")

# Symbolic total energy (1-electron part only, V_ee as perturbation)
E1_sym = (p0_s**2 / 2) * (1/beta1_s - 1)
E2_sym = (p0_s**2 / 2) * (1/beta2_s - 1)
V_NN_sym = ZA_s * ZB_s / sp.Symbol('R', positive=True)

# Express R in terms of γ and p₀: R = 1/(p₀ tan(γ/2))
R_of_gamma = 1 / (p0_s * sp.tan(gamma_s / 2))

E_total_sym = 2*E1_sym + 2*E2_sym + ZA_s * ZB_s * p0_s * sp.tan(gamma_s / 2)

print(f"\n  E_total(γ, p₀) = 2·E₁ + 2·E₂ + Z_A·Z_B·p₀·tan(γ/2)")
print(f"  E₁ = (p₀²/2)(1/β₁ - 1)")
print(f"  E₂ = (p₀²/2)(1/β₂ - 1)")
print(f"\n  Note: E₁ and E₂ are independent of γ!")
print(f"  So dE_total/dγ = Z_A·Z_B·p₀/(2cos²(γ/2))")
print(f"  This is ALWAYS POSITIVE for γ ∈ (0, π).")
print(f"\n  RESULT: E_total(γ) is monotonically increasing in γ.")
print(f"  The minimum is at γ → 0 (R → ∞), i.e., separated atoms.")
print(f"  With fixed β and no V_ee, the molecule is UNBOUND.")

# Verify symbolically
dE_dgamma = sp.diff(E_total_sym, gamma_s)
dE_simplified = sp.simplify(dE_dgamma)
print(f"\n  dE/dγ = {dE_simplified}")

# Check sign
print(f"\n  At γ=1, p₀=1, Z_A=3, Z_B=1:")
dE_num = float(dE_simplified.subs([(gamma_s, 1.0), (p0_s, 1.0),
                                    (ZA_s, 3), (ZB_s, 1)]))
print(f"  dE/dγ = {dE_num:.6f} (positive → repulsive)")


# ──────────────────────────────────────────────────────────────────────
# Step 3: Include p₀ self-consistency
# ──────────────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("STEP 3: Self-consistent p₀(R) and total energy E(R)")
print("=" * 70)

print("""
The self-consistency condition p₀ = √(-2E_total) means p₀ depends on R
through V_NN(R). Let me solve for p₀(R) at each R and compute E(R).

At each R:
    E_total(p₀) = p₀²(1/β₁ - 1) + p₀²(1/β₂ - 1) + Z_A Z_B / R
    E_total(p₀) = p₀² [(1/β₁-1) + (1/β₂-1)] + Z_A Z_B / R

    Self-consistency: p₀² = -2 E_total
    → p₀² = -2 p₀² [(1/β₁-1) + (1/β₂-1)] - 2 Z_A Z_B / R

    Let α_eff = (1/β₁-1) + (1/β₂-1)  [negative for β>1]
    → p₀² (1 + 2α_eff) = -2 Z_A Z_B / R
    → p₀² = -2 Z_A Z_B / (R (1 + 2α_eff))

    For this to have p₀² > 0, need (1 + 2α_eff) < 0 (since Z_A Z_B > 0, R > 0).
""")

alpha_eff = (1.0/beta1_val - 1.0) + (1.0/beta2_val - 1.0)
print(f"  α_eff = (1/β₁-1) + (1/β₂-1) = {alpha_eff:.6f}")
print(f"  1 + 2α_eff = {1 + 2*alpha_eff:.6f}")

if 1 + 2*alpha_eff < 0:
    print("  ✓ 1 + 2α_eff < 0 → self-consistent p₀² > 0 (bound state exists)")
else:
    print("  ✗ 1 + 2α_eff ≥ 0 → no bound self-consistent solution")
    print("  The 2-state model with these β values cannot form a bound state.")

# Compute self-consistent p₀(R) and E(R)
results = []

print(f"\n  {'R (bohr)':>10} {'p₀':>10} {'E_total (Ha)':>14} {'γ (rad)':>10} {'γ (deg)':>10}")
print("  " + "-" * 60)

for R_val in np.arange(1.0, 10.1, 0.5):
    # p₀² = -2 Z_A Z_B / (R (1 + 2α_eff))
    p0_sq = -2.0 * ZA * ZB / (R_val * (1.0 + 2.0 * alpha_eff))
    if p0_sq <= 0:
        continue
    p0_val = np.sqrt(p0_sq)
    E_total_val = total_energy_1e(p0_val, R_val)
    # Bond angle: cos γ = (p₀² - 1/R²)/(p₀² + 1/R²)
    pR = 1.0 / R_val
    cos_g = (p0_val**2 - pR**2) / (p0_val**2 + pR**2)
    cos_g = max(-1.0, min(1.0, cos_g))
    gamma_val = np.arccos(cos_g)

    results.append((R_val, p0_val, E_total_val, gamma_val))
    print(f"  {R_val:10.2f} {p0_val:10.4f} {E_total_val:14.6f} "
          f"{gamma_val:10.4f} {np.degrees(gamma_val):10.2f}")

print("""
OBSERVATION: E_total(R) = p₀²(R) · α_eff + Z_A Z_B / R
where p₀²(R) = -2 Z_A Z_B / (R(1+2α_eff)).

Substituting:
    E_total(R) = [-2 Z_A Z_B / (R(1+2α_eff))] · α_eff + Z_A Z_B / R
               = Z_A Z_B / R · [-2α_eff/(1+2α_eff) + 1]
               = Z_A Z_B / R · [1+2α_eff - 2α_eff] / (1+2α_eff)
               = Z_A Z_B / R · 1/(1+2α_eff)

So E_total(R) = Z_A Z_B / (R · (1+2α_eff))
This is a pure 1/R curve — monotonically decreasing (since 1+2α_eff < 0).
NO MINIMUM. The molecule is unbound in this model.
""")

# Verify
print("  Verification: E_total should be proportional to 1/R")
print(f"  E(1.0) × 2.0 = {results[0][2] * 2.0:.6f}")
print(f"  E(2.0)        = {results[2][2]:.6f}")
print(f"  Ratio: {results[0][2] * 2.0 / results[2][2]:.6f} (should be ~1.0)")


# ──────────────────────────────────────────────────────────────────────
# Step 4: What's missing — the β(R) dependence
# ──────────────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("STEP 4: The β(R) mechanism")
print("=" * 70)
print("""
The minimal model with FIXED β produces a pure 1/R potential — no minimum.
The binding mechanism requires β_k(R), i.e., the prolate spheroidal
eigenvalues must change with R.

Physical picture:
- At R → ∞: separated atoms, β₁ → Z_A/(Z_A+Z_B), β₂ → Z_B/(Z_A+Z_B)
- At R → 0: united atom (Z_A+Z_B), all β_k → 1
- At intermediate R: β_k(R) interpolates, and the interplay between
  p₀(R) and β(R) creates the energy minimum.

This means R_eq is NOT determinable from a fixed-β model. The
equilibrium geometry requires solving the prolate spheroidal equation
at each R to get β_k(R), then finding where E_total(R) is minimized.

Let me parametrize β(R) with a simple interpolation and check if
a minimum appears.
""")

# Model β(R) with reasonable interpolation
# At R=3.015 (from v0.9.29): β₁=1.033, β₂=1.950
# At R→∞: separated atoms
#   For Li (Z=3): β₁(∞) → n₁·Z_A/(Z_A+Z_B) = 1·3/4 = 0.75
#   For H  (Z=1): β₂(∞) → n₂·Z_B/(Z_A+Z_B) = 2·1/4 = 0.50
#   (These are rough estimates — exact values depend on the separation limit)
# At R→0: united atom (Be, Z=4)
#   All β → 1 (single-center Sturmian with Z=4)
#
# Actually, the β values are eigenvalues of the two-center Coulomb problem.
# They satisfy: β_k(R) → Z_k_eff(R) / Z_total approximately.
# The key physics is that β₂ > 1 at finite R (meaning the bonding orbital
# has MORE nuclear attraction than the total, due to both nuclei).

# Simple exponential interpolation
def beta_model(R_val, beta_eq, beta_inf, R_scale=3.0):
    """Simple model: β(R) = β_∞ + (β_eq - β_∞) exp(-(R-R_eq)²/R_scale²)"""
    # This is just for illustration — not physical
    return beta_inf + (beta_eq - beta_inf) * np.exp(-(R_val - 3.015)**2 / R_scale**2)

# More physical model: β from the prolate spheroidal equation
# At large R: the two centers decouple. The Sturmian eigenvalues approach
# the atomic values: β → Z_center/Z_total × n_orbital
# At small R: united atom, β → 1

# Use a physically motivated sigmoid
def beta1_of_R(R_val):
    """Model β₁(R) — Li 1s orbital."""
    # At R→0: β₁→1 (united atom)
    # At R→∞: β₁→0.75 (Li gets 3/4 of nuclear attraction)
    # At R=3.015: β₁=1.033
    # Use tanh interpolation
    x = R_val / 3.0  # dimensionless
    return 1.033 + 0.0 * (x - 1.0)  # approximately constant near R_eq
    # Actually let me use the physical constraint more carefully.
    # β₁ is the eigenvalue of the angular part of the prolate spheroidal eq.
    # For 1s-like states, it varies slowly with R.
    # Use linear interpolation between known limits:
    # R=0 → β₁=1, R=3.015 → β₁=1.033, R=∞ → β₁≈0.75
    if R_val < 3.015:
        return 1.0 + (1.033 - 1.0) * R_val / 3.015
    else:
        return 1.033 * np.exp(-0.05 * (R_val - 3.015))


def beta2_of_R(R_val):
    """Model β₂(R) — bonding 2σ orbital."""
    # At R→0: β₂→1 (united atom)
    # At R→∞: β₂→0.5 (H gets 1/4 of attraction, ×n=2)
    # At R=3.015: β₂=1.950
    # The bonding orbital is strongly R-dependent.
    if R_val < 3.015:
        return 1.0 + (1.950 - 1.0) * R_val / 3.015
    else:
        # Exponential decay toward dissociation limit
        return 0.5 + (1.950 - 0.5) * np.exp(-0.3 * (R_val - 3.015))


def E_total_with_beta_R(R_val):
    """Total energy with β(R) dependence and self-consistent p₀."""
    b1 = beta1_of_R(R_val)
    b2 = beta2_of_R(R_val)

    a_eff = (1.0/b1 - 1.0) + (1.0/b2 - 1.0)
    denom = 1.0 + 2.0 * a_eff

    if denom >= 0:
        return np.inf  # no bound solution

    p0_sq = -2.0 * ZA * ZB / (R_val * denom)
    if p0_sq <= 0:
        return np.inf

    p0_val = np.sqrt(p0_sq)

    E1 = (p0_val**2 / 2.0) * (1.0/b1 - 1.0)
    E2 = (p0_val**2 / 2.0) * (1.0/b2 - 1.0)
    V_NN = ZA * ZB / R_val

    return 2.0 * E1 + 2.0 * E2 + V_NN


print("\n  E_total(R) with model β(R):")
print(f"  {'R (bohr)':>10} {'β₁(R)':>10} {'β₂(R)':>10} {'E (Ha)':>14}")
print("  " + "-" * 50)

R_range = np.arange(1.0, 10.1, 0.25)
E_vals = []
for R_val in R_range:
    b1 = beta1_of_R(R_val)
    b2 = beta2_of_R(R_val)
    E_val = E_total_with_beta_R(R_val)
    E_vals.append(E_val)
    if R_val in [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0]:
        if E_val < 1000:
            print(f"  {R_val:10.2f} {b1:10.4f} {b2:10.4f} {E_val:14.6f}")
        else:
            print(f"  {R_val:10.2f} {b1:10.4f} {b2:10.4f} {'UNBOUND':>14}")

# Find minimum
E_arr = np.array(E_vals)
finite_mask = np.isfinite(E_arr)
if np.any(finite_mask):
    R_finite = R_range[finite_mask]
    E_finite = E_arr[finite_mask]
    idx_min = np.argmin(E_finite)
    R_min = R_finite[idx_min]
    E_min = E_finite[idx_min]
    print(f"\n  Minimum: R_eq = {R_min:.2f} bohr, E_min = {E_min:.6f} Ha")

    # Refine with minimize_scalar
    if R_min > R_finite[0] and R_min < R_finite[-1]:
        from scipy.optimize import minimize_scalar
        result = minimize_scalar(E_total_with_beta_R,
                                bounds=(max(0.5, R_min-1), R_min+1),
                                method='bounded')
        if result.success:
            R_eq_model = result.x
            E_eq_model = result.fun
            p0_eq = np.sqrt(-2.0 * ZA * ZB / (R_eq_model * (1.0 + 2.0 * ((1.0/beta1_of_R(R_eq_model)-1.0) + (1.0/beta2_of_R(R_eq_model)-1.0)))))
            gamma_eq = np.arccos((p0_eq**2 - 1/R_eq_model**2)/(p0_eq**2 + 1/R_eq_model**2))
            print(f"\n  Refined R_eq = {R_eq_model:.4f} bohr")
            print(f"  E_eq = {E_eq_model:.6f} Ha")
            print(f"  p₀* = {p0_eq:.4f}")
            print(f"  γ_eq = {gamma_eq:.4f} rad = {np.degrees(gamma_eq):.2f}°")
else:
    print("\n  No bound solution found for any R!")


# ──────────────────────────────────────────────────────────────────────
# Step 5: Sensitivity analysis
# ──────────────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("STEP 5: Sensitivity analysis")
print("=" * 70)

def find_R_eq(beta1_eq, beta2_eq, beta1_inf=None, beta2_inf=None):
    """Find R_eq for given β values using the model β(R)."""
    # Override the β(R) functions with scaled versions
    def b1(R):
        if R < 3.015:
            return 1.0 + (beta1_eq - 1.0) * R / 3.015
        else:
            return beta1_eq * np.exp(-0.05 * (R - 3.015))

    def b2(R):
        if R < 3.015:
            return 1.0 + (beta2_eq - 1.0) * R / 3.015
        else:
            return 0.5 + (beta2_eq - 0.5) * np.exp(-0.3 * (R - 3.015))

    def E_func(R_val):
        bb1 = b1(R_val)
        bb2 = b2(R_val)
        a_eff = (1.0/bb1 - 1.0) + (1.0/bb2 - 1.0)
        denom = 1.0 + 2.0 * a_eff
        if denom >= 0:
            return np.inf
        p0_sq = -2.0 * ZA * ZB / (R_val * denom)
        if p0_sq <= 0:
            return np.inf
        p0_val = np.sqrt(p0_sq)
        E1 = (p0_val**2 / 2.0) * (1.0/bb1 - 1.0)
        E2 = (p0_val**2 / 2.0) * (1.0/bb2 - 1.0)
        return 2.0 * E1 + 2.0 * E2 + ZA * ZB / R_val

    try:
        res = minimize_scalar(E_func, bounds=(0.5, 10.0), method='bounded')
        if res.success and np.isfinite(res.fun):
            return res.x, res.fun
    except Exception:
        pass
    return None, None

print("\n  Sensitivity of R_eq to β₁ and β₂:")
print(f"  {'β₁':>8} {'β₂':>8} {'R_eq':>10} {'E_eq':>12}")
print("  " + "-" * 45)

# Vary β₁ by ±10%
for db1 in [-0.10, -0.05, 0.0, 0.05, 0.10]:
    for db2 in [-0.10, 0.0, 0.10]:
        b1_test = beta1_val * (1.0 + db1)
        b2_test = beta2_val * (1.0 + db2)
        R_eq_test, E_eq_test = find_R_eq(b1_test, b2_test)
        if R_eq_test is not None:
            print(f"  {b1_test:8.4f} {b2_test:8.4f} {R_eq_test:10.4f} {E_eq_test:12.6f}")
        else:
            print(f"  {b1_test:8.4f} {b2_test:8.4f} {'UNBOUND':>10} {'---':>12}")


# ──────────────────────────────────────────────────────────────────────
# Final Summary
# ──────────────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("SUMMARY OF RESULTS")
print("=" * 70)

print("""
1. ANALYTIC γ_eq: DOES NOT EXIST (for fixed β)

   With fixed Sturmian eigenvalues β₁, β₂, the 1-electron energy is:
       E_k = (p₀²/2)(1/β_k - 1)

   This is INDEPENDENT of γ. The nuclear repulsion V_NN ∝ tan(γ/2)
   is monotonically increasing in γ. Therefore dE_total/dγ > 0
   everywhere — there is no minimum in γ.

   The Sturmian Hamiltonian H[i,j] ∝ S[i,j] (proportional to overlap),
   so the generalized eigenvalue problem gives eigenvalues that depend
   only on β_k and p₀, not on the molecular geometry.

2. FIXED-β SELF-CONSISTENT SOLUTION: UNBOUND

   With self-consistent p₀ = √(-2E), the total energy is:
       E_total(R) = Z_A Z_B / (R · (1 + 2α_eff))
   which is a pure 1/R potential — no minimum, no binding.

3. THE BINDING MECHANISM REQUIRES β(R)

   The R-dependence of the prolate spheroidal eigenvalues β_k(R) is
   essential. With a model β(R) interpolation, a minimum does appear.
   However, the position of the minimum depends sensitively on the
   β(R) functional form, which is NOT derivable from the minimal
   2-state model alone.

4. CONCLUSION: NEGATIVE RESULT

   The 2-state minimal Bond Sphere Hamiltonian with fixed β values
   CANNOT predict R_eq. The equilibrium geometry is encoded in the
   R-dependence of the prolate spheroidal eigenvalues, not in the
   SO(4) rotation angle γ.

   This is a fundamental limitation: the Sturmian construction
   absorbs all geometry into the β eigenvalues, making the Hamiltonian
   proportional to the overlap matrix. The γ-optimization is trivial
   (always prefers γ→0, i.e., R→∞).

   To predict R_eq from SO(4) representation theory, one needs either:
   (a) The full β_k(R) from the prolate spheroidal equation at each R, OR
   (b) A non-Sturmian formulation where H is not proportional to S.

   Option (b) is the LCAO approach already implemented in lattice_index.py,
   where the Hamiltonian has independent kinetic, nuclear, and bridge terms.
""")

# ──────────────────────────────────────────────────────────────────────
# Save results
# ──────────────────────────────────────────────────────────────────────

output_path = os.path.join(os.path.dirname(__file__), 'data',
                           'bond_sphere_hamiltonian_analysis.txt')
os.makedirs(os.path.dirname(output_path), exist_ok=True)

with open(output_path, 'w', encoding='utf-8') as f:
    f.write("BOND SPHERE MINIMAL HAMILTONIAN -- ANALYSIS RESULTS\n")
    f.write("=" * 60 + "\n")
    f.write(f"Date: 2026-03-11\n")
    f.write(f"Model: 2-state (1-sigma Li 1s, 2-sigma Li 2s/H 1s bonding)\n")
    f.write(f"Parameters: beta1={beta1_val}, beta2={beta2_val}, Z_A={ZA}, Z_B={ZB}\n\n")

    f.write("1. ANALYTIC gamma_eq: DOES NOT EXIST\n")
    f.write("   Sturmian H proportional to S -> eigenvalues independent of gamma.\n")
    f.write("   dE_total/dgamma = Z_A*Z_B*p0/(2cos^2(gamma/2)) > 0 always.\n\n")

    f.write("2. NUMERICAL SWEEP (fixed beta):\n")
    f.write(f"   {'R':>8} {'p0':>10} {'E_total':>14} {'gamma':>10}\n")
    f.write("   " + "-" * 48 + "\n")
    for R_val, p0_val, E_val, g_val in results:
        f.write(f"   {R_val:8.2f} {p0_val:10.4f} {E_val:14.6f} {g_val:10.4f}\n")

    f.write(f"\n3. SELF-CONSISTENT SOLUTION (fixed beta): NONE\n")
    f.write(f"   E_total(R) = Z_A Z_B / (R*(1+2*alpha_eff)) -- pure 1/R, no minimum.\n")
    f.write(f"   alpha_eff = {alpha_eff:.6f}, 1+2*alpha_eff = {1+2*alpha_eff:.6f}\n\n")

    f.write("4. COMPARISON:\n")
    f.write("   R_eq (Bond Sphere, fixed beta) = inf (unbound)\n")
    f.write("   R_eq (experiment)               = 3.015 bohr\n")
    f.write("   R_eq (LCAO v0.9.11)             = ~2.5 bohr\n\n")

    f.write("5. ROOT CAUSE:\n")
    f.write("   The Sturmian construction makes H[i,j] proportional to S[i,j].\n")
    f.write("   Eigenvalues E_k = (p0^2/2)(1/beta_k - 1) are geometry-independent.\n")
    f.write("   Binding requires beta_k(R) -- the prolate spheroidal PES.\n")
    f.write("   The minimal 2-state model with fixed beta is structurally incapable\n")
    f.write("   of predicting R_eq.\n")

print(f"\nResults saved to {output_path}")
print("\nDone.")
