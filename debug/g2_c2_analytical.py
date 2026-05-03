"""
Analytical derivation of the second curvature coefficient c₂ in the QED
anomalous magnetic moment expansion on unit S³.

The one-loop vertex correction on S³ produces:
    F₂/Schwinger = 1 + c₁/λ² + c₂/λ⁴ + ...
where λ is the Dirac eigenvalue cutoff |λ_n| = n + 3/2 on unit S³.

Known:
    c₁ = 1/2 (exact, Parker-Toms: R/12 = 6/12 for Dirac field on unit S³)

This script derives c₂ using TWO independent methods:

METHOD 1: Seeley-DeWitt heat kernel coefficient a₂
    The heat kernel expansion of the one-loop vertex correction produces
    curvature corrections order-by-order. The a₂ coefficient for the
    vertex form factor on a curved background involves R², Ric², Riem²,
    the Lichnerowicz endomorphism E = R/4, and the spin connection
    curvature Ω_μν.

METHOD 2: Direct spectral sum with multi-n_ext polynomial fit
    Using the converged spectral mode sum data at multiple external states
    n_ext = 1..7, fit the curvature expansion to extract c₂ cleanly
    (separated from higher-order contamination).

Geometry of unit S³:
    R = 6 (Ricci scalar)
    R_μν = 2g_μν (Ricci tensor)
    R_μνρσ = g_μρg_νσ - g_μσg_νρ (Riemann tensor, constant curvature K=1)
    |Ric|² = R_μνR^μν = 4·3 = 12
    |Riem|² = R_μνρσR^μνρσ = 2·3·2 = 12
    C_μνρσ = 0 (Weyl tensor vanishes — conformally flat)
    ΔR = 0 (constant curvature)

Dirac spectrum (Camporesi-Higuchi):
    |λ_n| = n + 3/2, g_n = 2(n+1)(n+2)

References:
    Parker & Toms, Phys. Rev. D 29 (1984) 1584 [curvature corrections to g-2]
    Vassilevich, Phys. Rep. 388 (2003) 279 [heat kernel review]
    Gilkey, "Invariance Theory, the Heat Equation, and the Atiyah-Singer Index Theorem"
    Branson & Gilkey, Comm. PDE 15 (1990) 245
    Bunch & Parker, Phys. Rev. D 20 (1979) 2499 [propagator in curved space]
"""

from __future__ import annotations

import json
import sys
import os
from fractions import Fraction

sys.path.insert(0, '.')

import sympy as sp
from sympy import (
    Rational, Integer, pi, sqrt, Symbol, simplify, factor,
    nsimplify, N, oo, summation, Piecewise, Abs,
)

import mpmath
mpmath.mp.dps = 80

print("=" * 75)
print("  ANALYTICAL DERIVATION OF c₂ IN g-2 CURVATURE EXPANSION ON S³")
print("=" * 75)


# =====================================================================
# SECTION 1: Geometry of unit S³
# =====================================================================

print("\n" + "=" * 75)
print("  SECTION 1: Geometric invariants on unit S³ (d=3)")
print("=" * 75)

d = 3
R_scalar = Integer(6)           # d(d-1)K = 3·2·1 = 6
K_sectional = Integer(1)        # sectional curvature
Ric_tensor_coeff = Integer(2)   # R_μν = 2g_μν = (d-1)Kg_μν
Ric_sq = Integer(12)            # R_μνR^μν = (d-1)²K²·d = 4·3 = 12
Riem_sq = Integer(12)           # R_μνρσR^μνρσ = 2d(d-1)K² = 12
Weyl_sq = Integer(0)            # Weyl = 0 (conformally flat in d=3)
delta_R = Integer(0)            # ΔR = 0 (constant curvature)

# Lichnerowicz endomorphism: D² = ∇*∇ + R/4
E_lich = R_scalar / 4           # = 3/2

# Spin connection curvature for S³ (constant curvature)
# Ω_μν = (1/4) R_μνab γ^a γ^b
# tr(Ω_μν Ω^μν) = (1/16) R_μνab R^μν_cd tr(γ^a γ^b γ^c γ^d)
# For d=3 with 2-component spinors: tr(γ^a γ^b γ^c γ^d) = 2(δ^ac δ^bd - δ^ad δ^bc + ε^abcd)
# Actually for d=3, we use 2x2 Pauli matrices, tr(I) = 2
# Ω_μν = -(1/8) R_μνρσ [γ^ρ, γ^σ] = -(1/4) R_μνρσ γ^ρ γ^σ (antisymmetry)
# On S³: R_μνρσ = g_μρ g_νσ - g_μσ g_νρ
# Ω_μν = -(1/4)(g_μρ g_νσ - g_μσ g_νρ) γ^ρ γ^σ = -(1/4)[γ_μ, γ_ν]
# tr(Ω_μν Ω^μν) over spinor indices:
# = (1/16) tr([γ_μ,γ_ν][γ^μ,γ^ν])
# [γ_μ,γ_ν] = 2γ_μγ_ν - 2g_μν·I  (using {γ_μ,γ_ν}=2g_μν)
# Actually: [γ_μ,γ_ν][γ^μ,γ^ν] = 4(γ_μγ_ν - g_μν)(γ^μγ^ν - g^μν)
# In d dimensions: γ_μγ_ν γ^μ γ^ν = (d-2)²·I (for d-dim Clifford algebra)
# Wait, let me be more careful.

# For the Vassilevich formula, we need:
# a₂(D²) for the SQUARED Dirac operator D² = ∇*∇ + E
# where E = R/4 (Lichnerowicz) and ∇ is the spin connection

# The general Vassilevich (2003) Eq. (4.3) formula for a Laplace-type operator
# P = -(g^μν ∂_μ∂_ν + A^μ ∂_μ + B) on a vector bundle:
#
# a₂(P) = (4π)^{-d/2} ∫ tr_V [(1/360)(5R² - 2|Ric|² + 2|Riem|² + 12ΔR)·I_V
#          + (1/6)ΔE + E² - (R/6)E + (1/12)Ω_μν Ω^μν ] dvol
#
# For P = D² = -∇_spin*∇_spin + R/4, where:
# E = R/4·I_S (proportional to identity on spinor fiber)
# Ω_μν = spin connection curvature = (1/4) R_μν^{ab} γ_a γ_b

# Let me compute tr_S(Ω_μν Ω^μν) for the spin connection on S³.
#
# On a constant-curvature space (K=1):
# R_μνρσ = g_μρ g_νσ - g_μσ g_νρ
# The spin connection curvature is:
# Ω_μν = (1/4) R_μνab Σ^{ab}  where Σ^{ab} = (1/2)[γ^a, γ^b]
# = (1/4)(g_μa g_νb - g_μb g_νa)(1/2)[γ^a, γ^b]
# = (1/4)(1/2)(γ_μ γ_ν - γ_ν γ_μ)·2 - ... let me be precise
# Actually Σ^{ab} = (i/4)[γ^a, γ^b] in Euclidean signature
# OR Σ^{ab} = (1/4)[γ^a, γ^b] depending on conventions.

# Using Vassilevich convention (matching his Eq. 2.6):
# Ω_μν = ∂_μ ω_ν - ∂_ν ω_μ + [ω_μ, ω_ν]
# where ω_μ is the spin connection 1-form.
# For constant curvature: Ω_μν = (1/4) R_μν^{ρσ} γ_ρ γ_σ (with appropriate factor)

# Key result: on a d-dim constant-curvature space of curvature K:
# tr_S(Ω_μν Ω^μν) = dim_S · (d(d-1)/8) · K²
# where dim_S is the spinor dimension.

# Actually let me compute this more carefully for d=3.

print(f"  d = {d}")
print(f"  R = {R_scalar}")
print(f"  R_μν = {Ric_tensor_coeff}·g_μν")
print(f"  |Ric|² = R_μνR^μν = {Ric_sq}")
print(f"  |Riem|² = R_μνρσR^μνρσ = {Riem_sq}")
print(f"  Weyl = 0 (conformally flat)")
print(f"  ΔR = 0 (constant curvature)")
print(f"  E (Lichnerowicz) = R/4 = {E_lich}")

dim_S_3d = 2  # dim of spinor in d=3 (Pauli matrices)
dim_S_4d = 4  # dim if we use 4-component Dirac as in the GeoVac computation

# For the GeoVac computation, we use the 4D Dirac spinor (matching g_n = 2(n+1)(n+2))
dim_S = dim_S_4d
print(f"  dim_S = {dim_S} (4-component Dirac spinor)")


# =====================================================================
# SECTION 2: Seeley-DeWitt a₂ coefficient for D² on unit S³
# =====================================================================

print("\n" + "=" * 75)
print("  SECTION 2: Seeley-DeWitt a₂ coefficient")
print("=" * 75)

# The existing computation in qed_vacuum_polarization.py uses the formula:
# a₂ = (4π)^{-d/2} · dim_S · (1/360) · [5R² - 2|Ric|² + 2|Riem|²
#       - 30·R·E + 60·E²] · Vol
#
# This is the PROPAGATOR a₂ (heat kernel for D²). For the VERTEX correction,
# we need the a₂ contribution to the magnetic form factor.

# Parker & Toms (1985) derived the curvature correction to the vertex function.
# Their result for the leading correction is c₁ = R/12.
#
# For the NEXT order, the structure is:
#   c₂ involves combinations of R², R_μν², R_μνρσ², ΔR, (R/6-E)², Ω²
#   all evaluated on S³.

# The key insight: the vertex correction F₂ on a curved background is related
# to the DIFFERENCE between the dressed and bare vertices. In the heat-kernel
# (proper-time) representation, each order in curvature corresponds to a
# Seeley-DeWitt coefficient.

# For the magnetic form factor at zero momentum transfer on a curved background,
# Parker-Toms showed the curvature expansion coefficient c₁ = R/12 comes from
# the a₁ heat kernel coefficient applied to the vertex diagram.

# The second-order coefficient c₂ involves a₂ applied to the vertex:

# From Vassilevich (2003) Eq. (4.3), the "integrand" of a₂ for D² = ∇*∇ + E:
# (per unit volume, traced over spinor indices)
#
# a₂_integrand = dim_S · (1/360) · [5R² - 2|Ric|² + 2|Riem|² + 12ΔR]
#              + (1/6)·tr_S(ΔE) + tr_S(E²) - (R/6)·tr_S(E) + (1/12)·tr_S(Ω²)
#
# On S³: ΔR = 0, ΔE = 0 (constant curvature => E = R/4 = const)
# tr_S(E) = dim_S · R/4
# tr_S(E²) = dim_S · (R/4)²
# tr_S(Ω_μν Ω^μν) needs careful computation.

# For the VERTEX form factor, not just the propagator, we need to understand
# how the heat kernel expansion maps to the form factor expansion. The
# curvature expansion of F₂ arises from expanding the electron propagator
# in the vertex diagram in a power series in curvature.

# In proper-time representation, the dressed electron propagator has the
# heat kernel expansion:
#   S_F(x,x';s) ~ (4πs)^{-d/2} exp(-m²s) Σ_k a_k(x,x') s^k
#
# The vertex correction involves two propagators and one photon propagator.
# In the LARGE eigenvalue limit (λ >> 1), the curvature expansion is:
#   F₂/F₂^flat = 1 + c₁/λ² + c₂/λ⁴ + ...
# where the coefficients c_k are determined by the Seeley-DeWitt coefficients.

# The relationship between the heat kernel a_k and the vertex curvature
# coefficients c_k is:
#
# The magnetic form factor correction at one loop involves the vertex integral
# with two electron propagators. Each propagator in the curvature expansion
# brings one power of a₁ ∝ R at order 1/λ², and either a₂ or a₁² at order
# 1/λ⁴.

# Parker-Toms result: c₁ = R/12 comes from the substitution R/6 in the a₁
# coefficient, divided by 2 (two propagators contribute R/6 each, but the
# vertex structure gives an overall R/12 for the magnetic form factor).

# For c₂, the expected structure from the heat kernel expansion is:
# c₂ = α₁·(R/12)² + α₂·[invariants from a₂/a₁]
#
# where the α coefficients encode the vertex-diagram combinatorics.

# The most natural guess from dimensional analysis:
# c₂ should be a rational combination of c₁² = (R/12)² = 1/4 and terms
# from a₂.

# Let me compute the a₂ integrand on unit S³:
a2_geom = Rational(1, 360) * (5*R_scalar**2 - 2*Ric_sq + 2*Riem_sq)
# = (1/360)*(180 - 24 + 24) = 180/360 = 1/2
print(f"  Geometric a₂ integrand: (1/360)(5R²-2|Ric|²+2|Riem|²) = {a2_geom}")

# E terms (Lichnerowicz E = R/4 = 3/2):
E_val = E_lich  # = 3/2
E_sq = E_val**2  # = 9/4
RE_term = R_scalar * E_val / 6  # R·E/6 = 6·(3/2)/6 = 3/2
a2_E = E_sq - RE_term  # E² - R·E/6 = 9/4 - 3/2 = 3/4
print(f"  E = R/4 = {E_val}")
print(f"  E² = {E_sq}")
print(f"  R·E/6 = {RE_term}")
print(f"  E² - R·E/6 = {a2_E}")

# Alternative form: (R/6 - E)² = (1 - 3/2)² = 1/4
RE_diff = R_scalar/6 - E_val  # = 1 - 3/2 = -1/2
RE_diff_sq = RE_diff**2  # = 1/4
print(f"  (R/6 - E)² = ({R_scalar/6} - {E_val})² = {RE_diff_sq}")

# The relation: E² - RE/6 = (R/6 - E)² + R²/36 - 2R²/36 = (R/6-E)² - R²/36
# Actually: (R/6-E)² = R²/36 - RE/3 + E² so E² - RE/6 = (R/6-E)² - R²/36 + RE/6
# Let me verify: (R/6-E)² = (1-3/2)² = 1/4
# E² - RE/6 = 9/4 - 3/2 = 3/4
# (R/6-E)² = 1/4, difference = 3/4 - 1/4 = 1/2
print(f"  Verify: E²-RE/6 - (R/6-E)² = {a2_E - RE_diff_sq}")

# The Vassilevich formula for a₂(D²) per unit volume, per spinor component:
# a₂^{V} = (1/360)(5R²-2Ric²+2Riem²+12ΔR) + E²-(R/6)E + (1/12)Ω²
# where Ω² is per-spinor-component (i.e. tr_S(Ω²)/dim_S is what appears)

# Actually, the standard formula (Vassilevich Eq. 4.3) for a SINGLE component is:
# a₂(x) = (1/360)(5R²-2Ric²+2Riem²+12ΔR)·1 + E²-(R/6)E + (1/12)Ω_μν²

# But for a BUNDLE, we need the TRACE over the fiber:
# a₂^{bundle} = dim_S·(1/360)(5R²-2Ric²+2Riem²+12ΔR) + tr(E²)-(R/6)tr(E) + (1/12)tr(Ω²)

# For the spinor bundle with E = (R/4)·I_S:
# tr(E) = dim_S·R/4
# tr(E²) = dim_S·(R/4)²
# The spin connection curvature:
# Ω_μν = (1/4)R_μν^{ρσ}Σ_{ρσ} where Σ_{ρσ} = (1/2)γ_ργ_σ (antisym part)
# tr(Ω_μν Ω^μν) needs the trace of products of gamma matrices.

# On constant-curvature S^d with K=1:
# R_μνρσ = δ_μρ δ_νσ - δ_μσ δ_νρ
# Ω_μν = (1/4)(δ_μρ δ_νσ - δ_μσ δ_νρ)Σ^{ρσ} = (1/4)Σ_{μν}
# = (1/8)[γ_μ, γ_ν]
#
# Ω_μν Ω^μν = (1/64)[γ_μ,γ_ν][γ^μ,γ^ν]
# Now [γ_μ,γ_ν] = 2γ_μγ_ν - 2δ_μν·I
# so [γ_μ,γ_ν][γ^μ,γ^ν] = 4(γ_μγ_ν - δ_μν)(γ^μγ^ν - δ^μν)
# = 4(γ_μγ_νγ^μγ^ν - 2δ_μν·γ^μγ^ν + δ_μν·δ^μν·I)
# = 4(γ_μγ_ν·(2δ^μν - γ^νγ^μ)·γ^ν ... this is getting complex.

# Standard result: in d dimensions,
# [γ_μ,γ_ν][γ^μ,γ^ν] = -d(d-1)·I (for the d-dim Clifford algebra)
# Wait: {γ_μ,γ_ν} = 2δ_μν, so γ_μγ_ν = δ_μν + (1/2)[γ_μ,γ_ν]
# [γ_μ,γ_ν] is antisymmetric, so Σ = (1/2i)[γ_μ,γ_ν] is the generator
# tr(Σ_μν Σ^μν) = 2^{[d/2]} · d(d-1)/4  (standard representation theory)
#
# For d=3, 2x2 Pauli matrices:
# Σ_μν = (1/2)[σ_μ,σ_ν] = i ε_μνρ σ_ρ
# Σ_μν Σ^μν = -ε_μνρ ε^μνσ σ_ρ σ_σ = -2δ_ρσ σ_ρ σ_σ = -2·3·I₂ = -6I₂
# tr(Σ_μν Σ^μν) = -12
#
# So: tr(Ω_μν Ω^μν) = (1/64)·tr([γ_μ,γ_ν][γ^μ,γ^ν])
# = (1/64)·(4·tr(Σ_μν Σ^μν)) = (1/16)·(-12) = -3/4 for 2-comp spinor
#
# For 4-component Dirac (doubling): multiply by 2
# tr_4(Ω_μν Ω^μν) = -3/2

# Actually, let me be more careful about the factor.
# Ω_μν = (1/4) R_μν^{ab} Σ_{ab}  (Vassilevich convention, Eq. 2.6)
# On S³: R_μνab = δ_μa δ_νb - δ_μb δ_νa (with K=1)
# So Ω_μν = (1/4)(δ_μa δ_νb - δ_μb δ_νa) Σ^{ab} = (1/4)·2·Σ_μν = (1/2)Σ_μν
# No wait, Σ_{ab} is already antisymmetrized, and R_μνab is antisymmetric in (a,b)
# So R_μνab Σ^{ab} = R_μν^{ab}·Σ_{ab} summing over a<b (or with factor 1/2 for both)
# Let me use: Ω_μν = (1/2)·(1/2)·R_μν^{ρσ}·(1/2)[γ_ρ,γ_σ]
# = (1/8) R_μν^{ρσ} [γ_ρ,γ_σ]
# On S³: R_μν^{ρσ} = δ_μ^ρ δ_ν^σ - δ_μ^σ δ_ν^ρ
# Ω_μν = (1/8)(δ_μ^ρ δ_ν^σ - δ_μ^σ δ_ν^ρ)[γ_ρ,γ_σ]
# = (1/8)([γ_μ,γ_ν] - [γ_ν,γ_μ]) = (1/8)·2·[γ_μ,γ_ν] = (1/4)[γ_μ,γ_ν]

# So Ω_μν = (1/4)[γ_μ,γ_ν]
# Ω_μν Ω^μν = (1/16)[γ_μ,γ_ν][γ^μ,γ^ν]
# For d=3 Dirac (2-comp): [γ_μ,γ_ν] = 2iε_μνρ γ_ρ (Pauli matrices)
# [γ_μ,γ_ν][γ^μ,γ^ν] = -4ε_μνρ ε^μνσ γ_ρ γ^σ = -4·2δ_ρσ·γ_ργ^σ
# = -8·(γ_ργ^ρ) = -8·3·I₂ = -24I₂
# So tr₂(Ω_μν Ω^μν) = (1/16)·(-24)·2 = -3

# For 4-component Dirac: tr₄(Ω_μν Ω^μν) = 2·(-3) = -6
# (the 4-comp Dirac doubles the spinor representation)

Omega_sq_trace_2comp = Rational(-3, 1)  # tr₂(Ω²) for 2-component
Omega_sq_trace_4comp = Rational(-6, 1)  # tr₄(Ω²) for 4-component
Omega_sq_trace = Omega_sq_trace_4comp   # using 4-component

print(f"\n  Spin connection curvature on S³:")
print(f"    Ω_μν = (1/4)[γ_μ,γ_ν]  (constant curvature K=1)")
print(f"    tr₂(Ω_μν Ω^μν) = {Omega_sq_trace_2comp}  (2-component)")
print(f"    tr₄(Ω_μν Ω^μν) = {Omega_sq_trace_4comp}  (4-component Dirac)")

# Now compute the full a₂ integrand (per unit volume):
# a₂_density = dim_S·(1/360)(5R²-2Ric²+2Riem²)
#            + tr(E²) - (R/6)tr(E) + (1/12)tr(Ω²)

a2_geom_bundle = dim_S * Rational(1, 360) * (5*R_scalar**2 - 2*Ric_sq + 2*Riem_sq)
a2_E_bundle = dim_S * E_val**2 - (R_scalar/6) * dim_S * E_val
a2_Omega_bundle = Rational(1, 12) * Omega_sq_trace

a2_total_density = a2_geom_bundle + a2_E_bundle + a2_Omega_bundle

print(f"\n  a₂ integrand (per unit volume, traced over spinor):")
print(f"    Geometric: dim_S·(1/360)(5R²-2Ric²+2Riem²) = {a2_geom_bundle}")
print(f"    E terms:   tr(E²) - (R/6)tr(E) = {a2_E_bundle}")
print(f"    Ω terms:   (1/12)tr(Ω²) = {a2_Omega_bundle}")
print(f"    Total a₂ density = {a2_total_density}")

# Per spinor component:
a2_per_comp = a2_total_density / dim_S
print(f"    Per spinor component: {a2_per_comp}")


# =====================================================================
# SECTION 3: From a₂ to vertex coefficient c₂
# =====================================================================

print("\n" + "=" * 75)
print("  SECTION 3: Vertex form factor from Seeley-DeWitt coefficients")
print("=" * 75)

# The connection between the heat kernel coefficients and the vertex
# correction curvature coefficients is established through the proper-time
# representation of the one-loop vertex diagram.
#
# In flat space, the one-loop vertex correction gives F₂ = α/(2π) (Schwinger).
# On a curved background, the electron propagator acquires curvature corrections
# from the heat kernel expansion.
#
# The dressed propagator in proper time:
#   S_F(x,x';s) = i ∫₀^∞ ds exp(-m²s) K(x,x';s)
#   K(x,x';s) ~ (4πs)^{-d/2} [a₀ + a₁·s + a₂·s² + ...]
#
# The vertex correction involves a momentum-space integral over two propagators.
# In the large-λ expansion (where λ is the eigenvalue cutoff standing in for
# mass/momentum), each power of s corresponds to 1/λ²:
#   s^k ↔ 1/λ^{2k}
#
# For the vertex diagram:
# F₂ = F₂^flat · [1 + β₁/λ² + β₂/λ⁴ + ...]
#
# where β₁ comes from a₁ and β₂ comes from a₂ and a₁².
#
# The Parker-Toms result c₁ = R/12 comes from a₁ = R/6 with a factor of 1/2
# from the vertex diagram kinematics (the magnetic moment involves only
# the symmetric part of the vertex).

# For the vertex correction, the curvature enters through TWO sources:
# 1. The electron propagator curvature corrections (a₁, a₂ of D²)
# 2. The vertex coupling modification in curved space

# The c₁ = R/12 result is a₁^{vertex} = (1/2)·a₁^{propagator} = (1/2)·(R/6) = R/12

# For c₂, by the same logic, we expect:
# c₂ = (1/2)·a₂^{propagator}/dim_S + cross-terms from (a₁)²

# The per-component a₁ for the propagator is R/6
a1_per_comp = R_scalar / 6  # = 1

print(f"\n  Per-component Seeley-DeWitt coefficients (propagator):")
print(f"    a₁/dim_S = R/6 = {a1_per_comp}")
print(f"    a₂/dim_S = {a2_per_comp}")

# c₁ = (1/2)·(R/6) = R/12 = 1/2  ✓
c1_check = Rational(1, 2) * a1_per_comp
print(f"    c₁ = (1/2)·a₁ = (1/2)·(R/6) = {c1_check}  (= R/12 = 1/2 ✓)")

# For the second-order vertex coefficient, we need the proper-time structure
# of the vertex diagram. The standard result (Parker-Toms framework) is:
#
# The vertex correction to F₂ at one loop involves:
#   - One insertion of a₂ from one propagator: contributes a₂/(dim_S·2)
#   - One insertion of a₁ from EACH propagator: contributes (a₁/dim_S)²·f₂
#     where f₂ is a numerical factor from the vertex integral
#
# The key factor is the SAME factor 1/2 that appears in c₁:
# The magnetic form factor extracts the SPIN-dependent part of the vertex,
# which carries a factor 1/2 relative to the full vertex function.

# Following Parker-Toms (1985), the general pattern is:
# c_k = (1/2)·(combined a_k coefficient from vertex diagram)
# where the vertex diagram's effective a_k is built from
# products of propagator heat kernel coefficients.

# For the vertex diagram with TWO electron propagators and ONE photon propagator,
# the curvature expansion of the FULL vertex at order 1/λ⁴ receives:
#
# Term A: a₂ from one electron propagator, a₀ from the other
#   → contributes proportional to a₂/dim_S
#
# Term B: a₁ from BOTH electron propagators
#   → contributes proportional to (a₁/dim_S)²
#
# The relative coefficients depend on the proper-time integrals.

# In the vertex diagram, the two propagators have proper times s₁ and s₂.
# The Schwinger parameter integral is:
# ∫₀¹ dx [x(1-x)]^{d/2-1}
# This gives different weights to the two terms.

# For d=3 spatial dimensions (d=4 spacetime), the Schwinger integral is:
# ∫₀¹ dx x(1-x) = 1/6  (for the magnetic moment)
# ∫₀¹ dx x²(1-x)² = 1/30  (for the next order)

# The precise combination depends on the vertex structure.
# The Parker-Toms result gives:
# c₁ = R/12 = (1/2)·(R/6)·1  [factor 1 from ∫dx/∫dx = 1]

# For c₂, the proper-time structure gives:
# c₂ = (1/2)·[a₂/dim_S + γ·(a₁/dim_S)²]
# where γ is a combinatorial factor from the vertex integral.

# The DIRECT approach: compute c₂ from the spectral sum and identify
# the heat kernel structure.


# =====================================================================
# SECTION 4: Direct spectral computation of c₂
# =====================================================================

print("\n" + "=" * 75)
print("  SECTION 4: Direct spectral computation via multi-n_ext fit")
print("=" * 75)

# Use the dirac_vertex module to compute the vertex correction at
# multiple n_ext values, then fit the curvature expansion.

# From the existing data:
# F₂(n_ext)/Schwinger = 1 + c₁/λ² + c₂/λ⁴ + c₃/λ⁶ + ...
# where λ = (2·n_ext + 3)/2

# We need converged data at multiple n_ext. The n_int=1 dominant term
# was computed exactly (algebraically) for n_ext=1..7.
# The full sum (n_int=0..25) is converged at n_ext=1.
# We can also compute at higher n_ext where convergence is faster.

# First, let me compute using the Racah approach from the qed_anomalous_moment
# module, which gives exact algebraic B(n_int) values.

# Actually, the most efficient approach is to use the spectral mode sum
# from qed_self_energy.vertex_correction_spectral, which is converged.

# Let me compute the vertex correction at n_ext = 0,1,2,3,4,5 with
# sufficient n_int.

print("\n  Computing vertex corrections at multiple n_ext values...")
print("  (Using shell-summed spectral mode sum)")

# Import the vertex correction
from geovac.qed_self_energy import vertex_correction_spectral

# For the g-2, the key quantity is the m_j-dependent vertex difference.
# However, the vertex_correction_spectral computes the shell-summed
# vertex correction (summed over all m_j), which is NOT the magnetic
# form factor.

# We need the m_j-resolved vertex from qed_anomalous_moment.
# But that is very slow for high n_int.

# Alternative: USE THE EXISTING DATA from the numerical computation.
# The existing scripts established:
# - At n_ext=1: F₂/S = 1.0844529... ± 2e-7 (from n_int=0..25)
# - c₁ = 1/2 (exact)
# - c₂ ~ 0.174 (approximate, contaminated by higher orders)

# The BEST approach for extracting c₂ is to use the existing
# multi-n_ext data from g2_multi_next.json and supplement it.

# Load existing data
try:
    with open('debug/data/g2_multi_next.json') as f:
        multi_data = json.load(f)
    print("  Loaded multi-n_ext data from g2_multi_next.json")

    # n_ext=1: high-precision converged (n_int=0..25)
    with open('debug/data/g2_extended_nint_v2.json') as f:
        v2_data = json.load(f)
    delta_1 = mpmath.mpf(str(v2_data['corrected_delta']))
    print(f"  n_ext=1: delta = {mpmath.nstr(delta_1, 15)} (n_int=0..25, tail-corrected)")

    # n_ext=2,3,4: from multi_data (n_int=0..8, less converged)
    delta_2 = mpmath.mpf(str(multi_data['multi_next_data']['2']['delta']))
    delta_3 = mpmath.mpf(str(multi_data['multi_next_data']['3']['delta']))
    delta_4 = mpmath.mpf(str(multi_data['multi_next_data']['4']['delta']))
    print(f"  n_ext=2: delta = {mpmath.nstr(delta_2, 12)} (n_int=0..8)")
    print(f"  n_ext=3: delta = {mpmath.nstr(delta_3, 12)} (n_int=0..8)")
    print(f"  n_ext=4: delta = {mpmath.nstr(delta_4, 12)} (n_int=0..8)")

    HAS_MULTI = True
except FileNotFoundError:
    print("  Multi-n_ext data not found, will compute...")
    HAS_MULTI = False


# =====================================================================
# SECTION 5: Analytical formula for c₂
# =====================================================================

print("\n" + "=" * 75)
print("  SECTION 5: Analytical derivation of c₂ from heat kernel")
print("=" * 75)

# The proper approach is to use the Bunch-Parker-Toms framework for
# the one-loop vertex correction in curved space.
#
# The key identity (Parker-Toms 1985, Bunch-Parker 1979):
# The vertex correction form factor F₂ at zero momentum transfer has
# the curvature expansion:
#   F₂/F₂^flat = 1 + Σ_k c_k / λ^{2k}
#
# where c_k involves the Seeley-DeWitt coefficients evaluated on the
# background geometry.
#
# For the MAGNETIC form factor:
# c₁ = (1/2)·(R/6) = R/12
#
# For c₂, using the Barvinsky-Vilkovisky (1985) form factor expansion:
# The second-order curvature correction to the vertex involves:
#
# (a) The a₂ contribution from one propagator (with a₀ from the other):
#     Proportional to a₂/(dim_S)
#
# (b) The (a₁)² contribution from both propagators:
#     Proportional to (a₁/(dim_S))²
#
# (c) The vertex modification in curved space at second order in curvature:
#     This involves the commutator of covariant derivatives, which on S³
#     is proportional to the Riemann tensor.

# The standard proper-time vertex integral gives (for d=4 spacetime):
# F₂ = (α/2π) ∫₀¹ dx (1-x) · det[...]^{-1/2} · exp[...] · [1 + R/(12λ²) + ...]
#
# The (1-x) factor is the magnetic-moment extraction kernel.
# At order 1/λ², the average of R/(12λ²) is trivially R/(12λ²) = c₁/λ².
# At order 1/λ⁴, we need:
#   ∫₀¹ dx (1-x) · [vertex-effective a₂] / ∫₀¹ dx (1-x)

# The vertex-effective a₂ contains:
# (i)  The a₂ coefficient from the propagator of the virtual electron
# (ii) The square of a₁ from the two propagators
# (iii) The curvature modification of the photon propagator

# On S³, the photon propagator also has curvature corrections.
# The Hodge-1 Laplacian eigenvalues are μ_q = q(q+2).
# In the curvature expansion, the photon propagator contributes
# an effective curvature mass proportional to the curvature.

# Let me take a more direct computational approach: compute c₂ from the
# spectral sum itself, using the known structure.

# The spectral sum for the vertex correction at n_ext has the form:
# Λ(n_ext) = Σ_{n_int} Σ_{q} W(n_ext,n_int,q)² · g(n_int) · d_T(q) /
#            (|λ(n_int)|⁴ · μ(q))
#
# where the sum is over internal electron levels and photon modes.
#
# The curvature expansion emerges when we expand this sum for large
# external eigenvalue λ_ext = (2n_ext+3)/2.

# The key observation: the coefficient c₂ is a UNIVERSAL geometric
# quantity determined by the manifold invariants. It does not depend
# on the specific external state — only on R, Ric², Riem², and the
# spin connection.

# From the spectral sum structure, the coefficients c_k are determined
# by the spectral properties of the Dirac operator and photon operator
# on S³.

# Let me compute directly: at large λ, the vertex correction has the form:
# Λ(λ)/Λ_flat = 1 + (a₁^eff/dim_S)/λ² + (a₂^eff/dim_S)/λ⁴ + ...
# where a_k^eff are the EFFECTIVE heat kernel coefficients for the vertex diagram.

# For the vertex diagram on S³, the magnetic form factor coefficient is:
# c₂ = (1/2) × [second-order curvature coefficient from vertex proper-time integral]

# Let me compute this from the propagator heat kernel coefficients.
# Using the notation:
# a₁^prop = R/6 (per spinor component, for D²)
# a₂^prop = computed above

# The vertex diagram's second-order coefficient involves
# the CONVOLUTION of heat kernel coefficients:
# The vertex integral has two proper-time parameters s₁, s₂.
# The curvature expansion of the integrand gives:
# [a₀ + a₁·s₁ + a₂·s₁² + ...] × [a₀ + a₁·s₂ + a₂·s₂² + ...]
# At second order: a₂·s₁²·a₀ + a₀·a₂·s₂² + (a₁)²·s₁·s₂

# After the Schwinger parameterization s₁ = s·x, s₂ = s·(1-x):
# The integral ∫₀¹ dx (1-x) f(x) gives:
# For a₂·s²·x²: ∫₀¹ x²(1-x)dx = 1/12
# For a₂·s²·(1-x)²: ∫₀¹ (1-x)³ dx = 1/4
# For (a₁)²·s²·x(1-x): ∫₀¹ x(1-x)²dx = 1/12
# Normalization: ∫₀¹ (1-x)dx = 1/2

# So the vertex effective a₂ per component =
# (a₂ · 1/12 + a₂ · 1/4 + a₁² · 1/12) / (1/2)
# = (a₂ · 4/12 + a₁² · 1/12) / (1/2)
# = (a₂ · 1/3 + a₁² · 1/12) / (1/2)
# = (2/3)·a₂ + (1/6)·a₁²

# Wait, I need to be more careful about which propagator is which.
# In the vertex diagram, one propagator carries the external momentum
# (and thus the magnetic moment), while the other doesn't.
# The (1-x) factor weights the "external" propagator.

# Actually, for the Schwinger parameter integral of the anomalous
# magnetic moment:
# F₂ = (α/π) ∫₀¹ dx x(1-x) / [...]
# The x(1-x) factor is the standard Schwinger kernel for F₂.

# Let me reconsider. In Schwinger's original calculation:
# F₂ = (α/2π) ∫₀¹ dx ∫₀¹ dy ∫₀¹ dz δ(x+y+z-1) · 2m²z(1-z) / [...]
# This reduces to ∫₀¹ dz (1-z) · [...] for the magnetic form factor.

# The curvature expansion coefficients c_k emerge from expanding the
# propagator denominators. At order 1/λ², there's a single insertion
# of a₁. At order 1/λ⁴, there are:
# (A) A single insertion of a₂
# (B) Two insertions of a₁ (one from each propagator)

# The CORRECT Schwinger-parameter weights for the magnetic form factor are:
#
# For F₂ (anomalous magnetic moment), Schwinger showed:
#   F₂ ∝ ∫₀¹ dα α(1-α)  [standard Feynman parameterization]
#
# At first order: c₁ = ∫₀¹ dα α(1-α) · (R/6) / ∫₀¹ dα α(1-α)
#               = R/6 · 1 = R/6
# But Parker-Toms get R/12... so there must be a factor of 1/2 from the
# vertex structure.

# Actually, c₁ = R/12 suggests the integral ratio is 1/2.
# This makes sense: the a₁ insertion modifies one of the two propagators,
# and the average over which propagator is modified gives a factor 1/2.

# At second order:
# Type A: a₂ inserted into one propagator
#   Weight: (1/2) × ∫₀¹ dα α(1-α) · h_A(α) / ∫₀¹ dα α(1-α)
#
# Type B: a₁ inserted into both propagators
#   Weight: ∫₀¹ dα α(1-α) · h_B(α) / ∫₀¹ dα α(1-α)

# For the proper-time representation:
# The two electron propagators have proper times s·α and s·(1-α).
# The heat kernel expansion at order s^k gives:
# s₁^j · s₂^{k-j} = s^k · α^j · (1-α)^{k-j}

# For the magnetic moment extraction (which requires the α(1-α) kernel):
# At order k=2 (for c₂):
# j=2: s₁² · s₂⁰ → s² · α² (a₂ from propagator 1, a₀ from propagator 2)
# j=0: s₁⁰ · s₂² → s² · (1-α)² (a₀ from propagator 1, a₂ from propagator 2)
# j=1: s₁ · s₂ → s² · α(1-α) (a₁ from both propagators)

# The weighted integrals:
# I₁ = ∫₀¹ dα α(1-α) · α² / ∫₀¹ dα α(1-α) = [∫α³(1-α)dα] / [∫α(1-α)dα]
#     = (1/20) / (1/6) = 6/20 = 3/10
# I₂ = ∫₀¹ dα α(1-α) · (1-α)² / ∫₀¹ dα α(1-α) = [∫α(1-α)³dα] / [∫α(1-α)dα]
#     = (1/20) / (1/6) = 3/10
# I₃ = ∫₀¹ dα α(1-α) · α(1-α) / ∫₀¹ dα α(1-α) = [∫α²(1-α)²dα] / [∫α(1-α)dα]
#     = (1/30) / (1/6) = 6/30 = 1/5

# Note: I₁ = I₂ by the symmetry α → 1-α.

# BUT: these are the weights for the FULL vertex function, not the
# magnetic form factor. The magnetic form factor extraction involves
# taking the derivative with respect to the photon momentum,
# which introduces an additional factor.

# For the magnetic form factor F₂:
# In Schwinger's parametrization, F₂ involves:
# F₂ ∝ ∫₀¹ dz z(1-z) / D(z)
# where D(z) is the denominator. The curvature corrections modify D(z).

# The a₁ correction gives: δD/D ∝ R/(6λ²)
# So c₁ = R/12 implies the factor from the z-integral and vertex structure
# is exactly 1/2.

# For a₂: the same factor 1/2 applies for the a₂ insertion (type A).
# For the a₁² term (type B), the factor is different because
# both propagators are modified simultaneously.

# Let me compute the Schwinger-parameter integrals precisely.
# The magnetic moment integral kernel is z(1-z) where z parameterizes
# the electron propagator fraction.
# I_norm = ∫₀¹ z(1-z) dz = 1/6

# Type A: a₂ from one propagator (say the first, carrying fraction z)
# I_A = ∫₀¹ z(1-z) · z dz = ∫₀¹ z²(1-z) dz = 1/12
# Plus the other propagator: ∫₀¹ z(1-z) · (1-z) dz = ∫₀¹ z(1-z)² dz = 1/12
# Total type A: 2 × (1/12) = 1/6
# Normalized: I_A/I_norm = (1/6)/(1/6) = 1

# Wait, but c₁ = R/12 = (1/2)·R/6 suggests the factor is 1/2, not 1.
# Let me reconsider.

# The proper-time integral for the vertex correction is:
# F₂ = (α/2π) × ∫ [ds₁ ds₂ ds₃ ...] × [magnetic vertex factor] × [propagators]
# where s₁, s₂ are electron proper times and s₃ is photon proper time.

# After the Feynman parameterization, with z = s₁/(s₁+s₂+s₃):
# F₂ = (α/2π) × I_vertex
# I_vertex = ∫₀¹ dz 2z(1-z) × [1/denominator] × [1 + curvature corrections]

# The curvature corrections at order R/λ² give:
# δI / I_flat = c₁/λ²
# with c₁ = R/12

# This means the vertex integral structure produces a factor 1/2:
# c₁ = R/6 × (integral factor) = R/6 × 1/2 = R/12

# So the integral factor for the a₁ term is 1/2.
# What about the a₂ term?

# Actually, I think the factor 1/2 in c₁ = R/12 comes from a different
# source. Let me reconsider the Parker-Toms derivation more carefully.

# Parker-Toms (1985) Eq. (2.12):
# For the Dirac field on a curved background, the leading curvature
# correction to the one-loop vertex function is:
#   Γ_μ = γ_μ × [1 + R/(12m²)]
# where m is the fermion mass. In our S³ computation, the role of m²
# is played by λ², the squared Dirac eigenvalue.

# Their derivation shows c₁ = R/12 arises from:
# (1/6)R from the heat kernel a₁ coefficient, times a factor of 1/2
# from the vertex structure (the magnetic form factor is half of the
# vertex function's spin-dependent part).

# For c₂, following the same logic:
# The proper-time representation gives contributions from:
# (a) a₂^{per comp} from one propagator (factor 1/2 from vertex structure)
# (b) (a₁^{per comp})² from both propagators (factor from proper-time integral)

# From Parker-Toms approach, the second-order correction should be:
# c₂ = (1/2)·a₂^{per comp} + (proper-time factor)·(a₁^{per comp})²

# For the vertex integral, the (a₁)² term comes from the integral:
# ∫₀¹ dz z(1-z) · z · (1-z) / ∫₀¹ dz z(1-z)
# = [∫z²(1-z)²dz] / [∫z(1-z)dz] = (1/30)/(1/6) = 1/5
# With the vertex factor: (1/2)·(1/5) = 1/10

# Wait, the cross-term should not have the same 1/2 vertex factor.
# In the a₁ term, the factor 1/2 comes because only ONE propagator
# is modified. In the a₁² term, BOTH are modified, so the combinatorial
# factor is different.

# Let me think about this differently using the known spectral data.

# From the spectral sum, we know c₂ ~ 0.17394.
# Let me check what combination of heat kernel invariants gives this.

a1_val = sp.Rational(1, 1)  # R/6 = 1 on unit S³
a2_val = a2_per_comp        # per-component a₂

print(f"\n  Testing combinations for c₂ ~ 0.17394:")
print(f"    a₁ (per comp) = R/6 = {a1_val}")
print(f"    a₂ (per comp) = {a2_val}")
print(f"    c₁ = R/12 = 1/2")

# Test c₂ = (1/2)·a₂_per_comp:
test1 = sp.Rational(1, 2) * a2_per_comp
print(f"    (1/2)·a₂ = {test1} = {float(test1):.6f}")

# Test c₂ = (1/2)·a₂ + something·(a₁)²:
# 0.17394 - (1/2)·a₂ = residual from a₁² term
# residual = 0.17394 - test1
residual_from_a2 = mpmath.mpf('0.17394') - mpmath.mpf(str(float(test1)))
print(f"    Residual after (1/2)·a₂: {float(residual_from_a2):.6f}")

# If residual = f · (a₁)² where a₁ = 1:
# f = residual / 1 = residual
# This should be a rational number from the proper-time integral.
print(f"    Factor for (a₁)² term: {float(residual_from_a2):.6f}")

# Alternatively, the expansion might not separate as simply as
# (1/2)·a₂ + f·a₁². Let me try a more general approach.

# The curvature expansion of the vertex correction on S³:
# c₂ = α₁ · a₂_per_comp + α₂ · (a₁_per_comp)²
# where α₁ + α₂ · ... gives c₂.

# On S³, a₁_per_comp = 1, a₂_per_comp is known.
# Let me compute a₂_per_comp numerically.
a2_pc_float = float(a2_per_comp)
print(f"\n  a₂ (per component) = {a2_per_comp} = {a2_pc_float:.10f}")


# =====================================================================
# SECTION 6: Direct Euler-Maclaurin approach to the spectral sum
# =====================================================================

print("\n" + "=" * 75)
print("  SECTION 6: Direct Euler-Maclaurin / asymptotic expansion")
print("=" * 75)

# Rather than going through the heat kernel framework for the vertex,
# let me compute c₂ directly from the spectral sum using the
# Euler-Maclaurin formula.

# The vertex correction spectral sum for the magnetic form factor at
# external level n_ext is:
#
# F₂(n_ext) = Σ_{n_int} B(n_ext, n_int)
#
# where B is the per-level magnetic contribution summed over all m_j.
#
# In the large n_ext limit: λ_ext = n_ext + 3/2 → ∞
# F₂/Schwinger → 1 + c₁/λ² + c₂/λ⁴ + ...
#
# The coefficients c₁, c₂, ... are determined by the asymptotic expansion
# of the spectral sum.

# The dominant contribution comes from the n_int = 1 level (99% at n_ext=1).
# For large n_ext, the n_int=1 contribution can be expanded analytically.

# The B(n_ext, n_int=1) involves CG coefficients between levels n_ext and 1.
# The CG coupling is dominated by the q_loop = n_ext ± 1 photon modes.

# For the LARGE λ limit, the CG coefficients have asymptotic expansions
# in 1/λ. This is the route to the analytical c₂.

# However, the CG coefficient computation is algebraically complex
# (involves sqrt of rationals growing with n).

# ALTERNATIVE: Use the n_int=1 exact algebraic expressions from
# debug/data/alpha_g_minus_2_exact_forms.json

try:
    with open('debug/data/alpha_g_minus_2_exact_forms.json') as f:
        exact_data = json.load(f)

    print("\n  Using exact n_int=1 data for n_ext=1..7:")
    alpha_em = mpmath.mpf('7.2973525693e-3')
    schwinger_mp = alpha_em / (2 * mpmath.pi)

    n_ext_vals = []
    lam_vals = []
    delta_vals = []  # F₂(nint=1)/Schwinger - something, or just raw ratio

    for r in exact_data['results']:
        n = r['n_ext']
        lam = mpmath.mpf(2*n + 3) / 2
        f2_nint1 = mpmath.mpf(str(r['F2_nint1_float']))
        ratio = f2_nint1 / schwinger_mp

        n_ext_vals.append(n)
        lam_vals.append(lam)
        delta_vals.append(ratio)

        print(f"    n_ext={n}: λ={float(lam):5.1f}  F₂(nint=1)/S = {float(ratio):.10f}")

    HAS_EXACT = True
except FileNotFoundError:
    print("  Exact data not found")
    HAS_EXACT = False


# =====================================================================
# SECTION 7: Polynomial fit and Richardson extrapolation
# =====================================================================

print("\n" + "=" * 75)
print("  SECTION 7: High-precision polynomial fit for c₁, c₂")
print("=" * 75)

# For the FULL vertex correction (not just n_int=1), we need the
# converged F₂/Schwinger at multiple n_ext.
#
# We have high-precision data at n_ext=1 (F₂/S = 1.08445292...).
# We need similar precision at n_ext=2,3,4,...

# Let me compute the full vertex correction at n_ext=2,3,4,5 using
# the qed_anomalous_moment module with the Racah approach.

# Actually, let me use a SMARTER approach:
# Compute the vertex correction using the spectral sum approach
# (qed_self_energy.vertex_correction_spectral) which sums over shells
# without m_j resolution. Then use the tree-level magnetic coupling
# to extract F₂.

# But wait — vertex_correction_spectral gives the SHELL-SUMMED vertex,
# not the magnetic form factor. The g-2 extraction requires the
# m_j-dependent difference.

# For an analytical result, let me instead derive c₂ from first principles
# using the heat kernel structure.

# KEY INSIGHT: On S³, all curvature invariants are determined by R alone
# (since S³ has constant curvature). Therefore:
# R² = 36, Ric² = 12, Riem² = 12, Weyl = 0, ΔR = 0
# ALL curvature invariants at order R² are proportional to R².
# This means c₂ must be of the form c₂ = β · R² where β is a rational
# number (up to the a₁² contribution).

# Since R² enters as 36 on unit S³, and c₂ ~ 0.174, we have:
# c₂/R² ~ 0.174/36 ~ 0.00483
# c₂ ~ 36 × 0.00483 ~ 0.174

# Actually, c₂ should have the form:
# c₂ = f₁ · (R/12)² + f₂ · a₂_contributions_from_S³
# where f₁, f₂ are rational fractions from the vertex integral.

# Since c₁ = R/12 = 1/2, we have c₁² = 1/4.
# And c₂ ~ 0.174, which is less than c₁² = 0.25.
# The difference is c₂ - c₁² = 0.174 - 0.25 = -0.076.

# Let me try a DIFFERENT decomposition:
# c₂ = α · c₁² + β · (curvature from a₂)
# On S³: a₂_geom = 1/2, E terms = 3/4, Ω terms = -1/2 (total = 3/4)
# Per component: 3/16

# Actually, let me just use the NUMERICAL approach to get c₂ precisely,
# then identify it analytically.

# The best approach: use the existing converged data at n_ext=1 and
# compute convergent spectral sums at n_ext=2..5 using the anomalous
# moment module.

# Actually — let me take the smartest approach. The vertex correction
# spectral sum for the magnetic form factor has been computed with
# exact algebraic CG coefficients. The per-level contributions B(n_int)
# are known for n_int=0..25 at n_ext=1.

# The curvature expansion F₂/S = 1 + c₁/λ² + c₂/λ⁴ + ... with λ = 5/2
# gives:
# 1 + c₁·(4/25) + c₂·(16/625) + c₃·(64/15625) + ...
# = 1 + 0.16·c₁ + 0.0256·c₂ + 0.004096·c₃ + ...

# With delta = 0.08445292265524729 and c₁ = 1/2:
# 0.08445292... = 0.08 + 0.0256·c₂ + 0.004096·c₃ + ...
# 0.00445292... = 0.0256·c₂ + 0.004096·c₃ + ...
# c₂ ~ 0.00445292/0.0256 ~ 0.17394  (if c₃ term small)

# But at n_ext=1, the c₃ contamination is ~16% (0.004096/0.0256 ~ 0.16).
# To get c₂ cleanly, we need data at MULTIPLE n_ext values.

# Let me compute via a different route: the LARGE-λ asymptotics of the
# vertex correction sum.

# For the spectral sum at large n_ext → ∞:
# B(n_ext, n_int) ~ C(n_int) / λ_ext^p
# The coefficient c₂ is determined by the asymptotic structure.

# APPROACH: Fit c₁, c₂, c₃ simultaneously from data at n_ext=1,2,3,4
# using a Vandermonde system.

print("\n  Fitting curvature expansion from multi-n_ext data:")

if HAS_MULTI:
    # Use the 4 data points: n_ext=1,2,3,4
    n_ext_fit = [1, 2, 3, 4]
    lam_fit = [mpmath.mpf(2*n+3)/2 for n in n_ext_fit]
    x_fit = [1/lam**2 for lam in lam_fit]  # 1/λ²

    # delta values (F₂/S - 1)
    delta_fit = [delta_1, delta_2, delta_3, delta_4]

    print(f"  Data points:")
    for i in range(4):
        print(f"    n_ext={n_ext_fit[i]}: λ={float(lam_fit[i]):.1f}  "
              f"x=1/λ²={float(x_fit[i]):.6f}  delta={float(delta_fit[i]):.10f}")

    # PROBLEM: n_ext=2,3,4 are only computed with n_int=0..8, which is
    # less converged than n_ext=1 (n_int=0..25). The n_ext=2 value
    # delta=0.613... looks wrong (should decrease with n_ext since
    # curvature corrections are smaller at larger λ).

    # Wait — delta at n_ext=2 is LARGER than at n_ext=1?!
    # delta_1 = 0.0845 at λ=2.5
    # delta_2 = 0.613 at λ=3.5

    # That's wrong — delta should be delta_1 · (λ₁/λ₂)² approximately.
    # 0.0845 × (2.5/3.5)² = 0.0845 × 0.510 = 0.0431 expected.
    # Getting 0.613 means the multi_next computation used a different
    # normalization!

    print("\n  WARNING: n_ext=2,3,4 data appears to use different normalization!")
    print("  Need to recompute at consistent normalization.")
    print("  Falling back to single-point extraction.")

    # The multi_data values look like they computed the RATIO differently.
    # Let me examine: at n_ext=2, delta=0.613 with V_mag=0.389.
    # The tree-level V_mag at n_ext=1 is about 0.607.
    # So the issue is likely that delta includes the n_ext-dependent
    # tree-level normalization change.

    # Let me focus on the single-point extraction with careful error analysis.

    # From n_ext=1:
    # delta = 0.08445292265524729
    # x = 1/λ² = 4/25 = 0.16
    # delta = c₁·x + c₂·x² + c₃·x³ + ...

    x1 = mpmath.mpf(4) / 25  # 1/λ² at n_ext=1
    delta1 = delta_1

    # c₁ = R/12 = 1/2 (exact)
    c1_exact = mpmath.mpf(1) / 2

    # residual after c₁ subtraction
    res1 = delta1 - c1_exact * x1
    print(f"\n  Single-point c₂ extraction at n_ext=1:")
    print(f"    x = 1/λ² = 4/25 = {float(x1):.6f}")
    print(f"    delta = {mpmath.nstr(delta1, 15)}")
    print(f"    c₁·x = {float(c1_exact * x1):.10f}")
    print(f"    residual = delta - c₁·x = {mpmath.nstr(res1, 15)}")

    # res1 = c₂·x² + c₃·x³ + ...
    # c₂_apparent = res1/x² (contaminated by c₃·x + ...)
    c2_apparent = res1 / x1**2
    print(f"    c₂ (apparent) = res1/x² = {mpmath.nstr(c2_apparent, 15)}")
    print(f"    Contamination: c₃·x ≈ c₃·{float(x1):.3f}")
    print(f"    If c₃ ≈ c₂ ≈ 0.17, contamination ≈ {0.17*float(x1)*100:.1f}%")

    # Better: try rational identification of c₂
    c2_mp = c2_apparent
    print(f"\n  PSLQ identification of c₂ ≈ {float(c2_mp):.10f}:")

    # Test against specific heat kernel combinations
    # On unit S³:
    # R = 6, R² = 36
    # c₁ = R/12 = 1/2
    # c₁² = 1/4
    # a₁^{per comp} = R/6 = 1
    # a₂^{per comp} = depends on formula

    # The a₂_per_comp from the existing code (qed_vacuum_polarization.py):
    from geovac.qed_vacuum_polarization import seeley_dewitt_coefficients_s3
    sd = seeley_dewitt_coefficients_s3()
    print(f"\n  Seeley-DeWitt from existing code:")
    print(f"    a₀ = {sd['a0']}")
    print(f"    a₁ = {sd['a1']}")
    print(f"    a₂ = {sd['a2']}")
    print(f"    integrand_a₂ = {sd['integrand_a2']}")
    print(f"    E_lich = {sd['E_lich']}")

    # The a₂ from the code uses:
    # integrand = 5R² - 2Ric² + 2Riem² - 30RE + 60E²
    # = 180 - 24 + 24 - 30·6·(3/2) + 60·(3/2)²
    # = 180 - 24 + 24 - 270 + 135
    # = 45
    # a₂ = (4π)^{-3/2} · 4 · (1/360) · 45 · Vol(S³)
    # = (4π)^{-3/2} · 4 · (1/8) · 2π²
    # = (4π)^{-3/2} · π²

    integrand_check = 5*36 - 2*12 + 2*12 - 30*6*sp.Rational(3,2) + 60*sp.Rational(3,2)**2
    print(f"    Integrand check: {integrand_check}")
    # = 180 - 24 + 24 - 270 + 135 = 45 ✓

    # The per-component a₂ density (without the (4π)^{-d/2} and Vol factors):
    a2_density_per_comp = sp.Rational(1, 360) * integrand_check
    print(f"    a₂ density per component = (1/360)·{integrand_check} = {a2_density_per_comp}")
    print(f"    = {float(a2_density_per_comp):.10f}")

    # But wait — this is using the "propagator" Seeley-DeWitt formula,
    # which gives the heat kernel coefficients for the squared Dirac operator.
    # For the vertex correction, we need a different formula.

    # Let me try specific rational candidates for c₂:
    candidates = []

    # From the a₂ formula: a₂_density = 1/8 per spinor component
    # (45/360 = 1/8)
    a2_pc = sp.Rational(1, 8)

    # Various combinations:
    candidates.append(("c₁² = 1/4", 0.25))
    candidates.append(("c₁²/2 = 1/8", 0.125))
    candidates.append(("a₂_pc = 1/8", 0.125))
    candidates.append(("c₁·a₂_pc = 1/16", 1/16))
    candidates.append(("c₁² + a₂_pc/4 = 9/32", 9/32))
    candidates.append(("c₁² - c₁/4 = 1/8", 0.125))
    candidates.append(("c₁² - a₂_pc = 1/8", 0.125))
    candidates.append(("2c₁²/3 + a₂_pc/6 = 7/48", 7/48))
    candidates.append(("c₁²(1 + 1/(2·3)) = 7/24", 7/24))
    candidates.append(("c₁·(c₁ - a₂_pc) = 3/16", 3/16))
    candidates.append(("3/16", 3/16))
    candidates.append(("c₁² - 1/12 = 5/24", 5/24))
    candidates.append(("7/40 = 7·Δ", 7/40))
    candidates.append(("1/(4π²/R²) = R²/(4π²) on S³ = 36/(4π²) = 9/π²",
                       9/mpmath.pi**2))
    candidates.append(("19/100", 19/100))

    # Per the Parker-Toms framework:
    # c₂ = (1/2)·a₂_eff_vertex
    # where a₂_eff_vertex = a₂_propagator + cross_term
    # On S³, a₂_prop (per comp) = 1/8
    # The cross term involves (a₁)² with a proper-time integral weight.
    # For the magnetic moment: weight = ∫₀¹ z(1-z)·z·(1-z)dz / ∫₀¹z(1-z)dz
    # = (1/30)/(1/6) = 1/5
    # So a₂_eff_vertex = 2·(1/8) + (1)²·(1/5) = 1/4 + 1/5 = 9/20
    # c₂ = (1/2)·(9/20) = 9/40

    candidates.append(("9/40 (PT vertex formula)", 9/40))

    # But wait: I need to check whether the factor 1/2 applies to
    # the cross term as well. The (a₁)² term modifies BOTH propagators,
    # so the vertex extraction factor may differ.

    # Alternative: c₂ = c₁·a₂_pc + (vertex factor)·c₁²
    # = (1/2)·(1/8) + f·(1/4)
    # 0.17394 = 1/16 + f/4
    # f = 4·(0.17394 - 1/16) = 4·0.11144 = 0.44577
    # f ≈ 4/9 = 0.4444...
    candidates.append(("1/16 + 1/9 = 25/144", 25/144))

    # Or: c₂ = c₁² · f where f is the "vertex enhancement factor"
    # f = 0.17394/0.25 = 0.69577
    # f ≈ 25/36? = 0.6944...
    candidates.append(("c₁²·25/36 = 25/144", 25/144))

    # Let me try: c₂ = c₁² · (1 + 1/5) · (1/2) + ...
    # Hmm, let me just PSLQ it.

    basis_labels = ['c₂', '1']
    basis_vals = [c2_mp, mpmath.mpf(1)]
    try:
        rel = mpmath.pslq(basis_vals, tol=1e-5, maxcoeff=1000)
        if rel and rel[0] != 0:
            ratio = -mpmath.mpf(rel[1]) / rel[0]
            print(f"\n  PSLQ rational: c₂ ≈ {mpmath.nstr(ratio, 20)}")
            print(f"    = {rel[1]}/{-rel[0]} (coefficients)")
    except Exception as e:
        print(f"    PSLQ failed: {e}")

    # Try with pi²
    basis_labels = ['c₂', '1', 'π²']
    basis_vals = [c2_mp, mpmath.mpf(1), mpmath.pi**2]
    try:
        rel = mpmath.pslq(basis_vals, tol=1e-5, maxcoeff=200)
        if rel and rel[0] != 0:
            terms = []
            for r, l in zip(rel, basis_labels):
                if r != 0:
                    terms.append(f"{r}·{l}")
            print(f"  PSLQ (rational+π²): {' + '.join(terms)} = 0")
    except:
        pass

    # Sort candidates by distance from measured c₂
    print(f"\n  Candidates for c₂ ≈ {float(c2_mp):.8f}:")
    sorted_cands = sorted(candidates, key=lambda x: abs(x[1] - float(c2_mp)))
    for label, val in sorted_cands[:10]:
        err = abs(float(val) - float(c2_mp))
        pct = err / float(c2_mp) * 100
        print(f"    {label:40s} = {float(val):.8f}  err={err:.2e} ({pct:.2f}%)")


# =====================================================================
# SECTION 8: Precise heat kernel vertex formula
# =====================================================================

print("\n" + "=" * 75)
print("  SECTION 8: Parker-Toms vertex formula — systematic derivation")
print("=" * 75)

# The one-loop vertex correction in the heat kernel (proper-time)
# representation on a curved background:
#
# The electron propagator in proper time has the Seeley-DeWitt expansion:
#   <x|exp(-sD²)|x> ∝ s^{-d/2} [1 + a₁·s + a₂·s² + ...]
# where a₁ = R/6 (per component), a₂ = 1/8 (per component on unit S³)
#
# The one-loop vertex correction involves two electron propagators and
# one photon propagator. In the Schwinger parametrization:
#   s_total = s  (total proper time)
#   s₁ = s·α (first electron), s₂ = s·(1-α) (second electron)
#   The photon propagator contributes a mass-like term.
#
# For the anomalous magnetic moment (Schwinger result α/(2π)):
# F₂^{flat} = (α/2π) × 1  [from the integral ∫₀¹ dα 2α(1-α) = 1/3, etc.]
#
# On curved space, the curvature corrections from the two propagators give:
# F₂/F₂^{flat} = 1 + ⟨curvature⟩_{vertex} / λ² + ...
#
# At order 1/λ²: the two propagators contribute a₁ each, weighted by their
# proper-time fractions. The magnetic moment kernel is 2α(1-α).
# The normalized weight is:
#   <a₁>_vertex = ∫₀¹ dα·2α(1-α)·[a₁·α + a₁·(1-α)] / ∫₀¹ dα·2α(1-α)
#              = a₁·∫₀¹ dα·2α(1-α) / ∫₀¹ dα·2α(1-α) = a₁
# c₁ = a₁/2 = R/12 (the factor 1/2 comes from the way curvature enters
# the propagator expansion relative to the vertex normalization).

# Actually, I realize the factor of 1/2 has a more fundamental origin.
# It comes from the structure of the Schwinger integral for F₂:
#
# F₂ = (α/2π) ∫₀¹ dα ∫₀^∞ ds (propagator factor) × [curv. corrections]
#
# The curvature correction a₁·s contributes after the s-integration:
# ∫₀^∞ ds s · exp(-s·λ²) = 1/λ⁴
# vs the leading term ∫₀^∞ ds exp(-s·λ²) = 1/λ²
# The ratio is 1/λ² × (weight factor from α integral).

# Let me compute the weight factors systematically.

# The magnetic moment integral (Schwinger parametrization) gives F₂:
# F₂ ∝ ∫₀¹ dα ∫₀^∞ ds·s^{d/2-2}·exp(-s·m²)·2m²·α(1-α)·[heat kernel expansion]

# For the heat kernel expansion of the TWO electron propagators:
# K_1(s₁) × K_2(s₂) where s₁ = s·α, s₂ = s·(1-α)
# At order s^k: sum over j of (a_j from prop 1)·(a_{k-j} from prop 2)·s₁^j·s₂^{k-j}
# = s^k · sum_j a_j·a_{k-j}·α^j·(1-α)^{k-j}

# The s-integration gives a factor Γ(k+d/2-1)/m^{2(k+d/2-1)} → 1/λ^{2k} at large λ.

# The α-integration gives the vertex weight factors.

# ORDER k=1 (c₁):
# j=0: a₀·a₁·(1-α)  → weight ∫₀¹ 2α(1-α)²dα = 2·B(2,3) = 2·1/12 = 1/6
# j=1: a₁·a₀·α      → weight ∫₀¹ 2α²(1-α)dα = 2·1/12 = 1/6
# Total weight: (1/6 + 1/6) = 1/3
# Normalization: ∫₀¹ 2α(1-α)dα = 1/3
# Ratio: 1/3 / (1/3) = 1, times a₁ = R/6
# But c₁ = R/12 = a₁/2, so there's a factor of 1/2 somewhere.

# AH — the factor 1/2 comes from the OVERALL normalization of the
# Schwinger integral. The flat-space F₂ = α/(2π) involves:
# F₂^{flat} = (α/2π) × ∫₀¹ 2α(1-α)dα × ... = (α/2π) × (normalization)
# The ratio F₂/F₂^{flat} divides out the normalization, BUT the curvature
# correction enters through the propagator, which has a RELATIVE weight
# of 1/2 in the vertex structure.

# Actually, I think the issue is simpler. The effective mass in the
# Schwinger denominator is m² = λ²(1-α) + other terms. The a₁ correction
# modifies one propagator, and the effective contribution is:
# a₁ × ∫₀¹ 2α(1-α)·α dα / ∫₀¹ 2α(1-α) dα
# = a₁ × (1/6)/(1/3) = a₁/2 = R/12 ✓

# Great, so the weight for a₁ from the FIRST propagator (with proper time sα) is:
w1_a1 = mpmath.mpf(1) / 6  # ∫₀¹ 2α²(1-α)dα

# And from the second propagator:
w2_a1 = mpmath.mpf(1) / 6  # ∫₀¹ 2α(1-α)²dα  (same by symmetry)

# Normalization:
w_norm = mpmath.mpf(1) / 3  # ∫₀¹ 2α(1-α)dα

# c₁ = a₁ × (w1_a1 + w2_a1) / (2 × w_norm)  ???
# Actually c₁ = a₁ × w1_a1/w_norm = a₁ × (1/6)/(1/3) = a₁/2 = R/12 ✓
# (one propagator at a time, summed: 2 × a₁/2 = a₁; but c₁ = a₁/2
# suggests only ONE propagator contributes the magnetic vertex correction)

# The asymmetry between the two propagators comes from the vertex structure:
# the incoming electron (carrying the external momentum/spin) contributes
# differently from the outgoing one.

# For the magnetic form factor, the extraction projects out the ANTISYMMETRIC
# part in the external spin. This means only ONE propagator's curvature
# correction contributes to F₂ (the other goes into F₁).

# c₁ = a₁ × w₁/w_norm = (R/6) × (1/6)/(1/3) = (R/6) × (1/2) = R/12 ✓

# ORDER k=2 (c₂):
# The curvature corrections at order s² come from:
# j=0: a₀·a₂·(1-α)² → one propagator with a₂ correction
# j=2: a₂·a₀·α²     → the other propagator with a₂ correction
# j=1: a₁·a₁·α(1-α) → both propagators with a₁ corrections

# For the magnetic form factor (only one propagator contributes):
# The j=2 term (a₂ from the "magnetic" propagator):
w_a2 = mpmath.mpf(0)  # ∫₀¹ 2α³(1-α)dα
# 2·B(4,2) = 2·Γ(4)Γ(2)/Γ(6) = 2·6·1/120 = 12/120 = 1/10
w_a2_prop1 = mpmath.mpf(1) / 10  # ∫₀¹ 2α³(1-α)dα

# The j=0 term (a₂ from the "non-magnetic" propagator):
# ∫₀¹ 2α(1-α)³dα = 2·B(2,4) = 2·Γ(2)Γ(4)/Γ(6) = 2·1·6/120 = 1/10
w_a2_prop2 = mpmath.mpf(1) / 10  # same by Beta function symmetry

# The j=1 term (a₁ from both propagators):
# ∫₀¹ 2α²(1-α)²dα = 2·B(3,3) = 2·Γ(3)²/Γ(6) = 2·4/120 = 8/120 = 1/15
w_a1_a1 = mpmath.mpf(1) / 15  # ∫₀¹ 2α²(1-α)²dα

print(f"  Schwinger parametrization weights:")
print(f"    w_norm    = ∫₀¹ 2α(1-α)dα    = {w_norm}")
print(f"    w₁(a₁)   = ∫₀¹ 2α²(1-α)dα   = {w1_a1}  (one-propagator a₁)")
print(f"    w(a₂,p1) = ∫₀¹ 2α³(1-α)dα   = {w_a2_prop1}  (a₂ from prop 1)")
print(f"    w(a₂,p2) = ∫₀¹ 2α(1-α)³dα   = {w_a2_prop2}  (a₂ from prop 2)")
print(f"    w(a₁²)   = ∫₀¹ 2α²(1-α)²dα  = {w_a1_a1}  (a₁ from both props)")

# For the magnetic form factor, which propagator contributes?
# In the Schwinger parameterization, the "magnetic" propagator is the one
# carrying the external spinor index that gets projected by the m_j extraction.
# This is the propagator with proper time s·α (the first one).
# The second propagator (proper time s·(1-α)) goes into F₁.

# So for c₂ (magnetic form factor):
# c₂ = a₂_per_comp × w_a2_prop1/w_norm + (a₁_per_comp)² × w_a1_a1/w_norm
# plus the a₂ from prop2 contributes to F₁, NOT F₂.

# Wait, but the j=1 term involves BOTH propagators, so its contribution
# to F₂ vs F₁ needs to be determined by the vertex structure.
# For the magnetic form factor, the cross-term contributes partially to F₂
# and partially to F₁.

# Actually, the decomposition into F₁ and F₂ is done AFTER the loop integral.
# Both propagators contribute to both form factors. The extraction of F₂
# is done by projecting out the spin-dependent (σ_μν q^ν) part of the vertex.

# For the LEADING order (c₁):
# The spin projection gives F₂ ∝ ∫₀¹ 2α(1-α) × f_spin(α) × a₁
# where f_spin(α) is the spin-structure kernel.
# For Schwinger's calculation: f_spin(α) = α for the magnetic vertex,
# giving ∫₀¹ 2α²(1-α) × a₁ = (1/6)·a₁, divided by normalization (1/3),
# giving a₁/2 = R/12. ✓

# For the SECOND order (c₂):
# Same spin structure: f_spin(α) = α
# j=2 contribution: a₂·α²·f_spin(α) = a₂·α³
#   → ∫₀¹ 2α³(1-α)dα = 1/10, normalized: (1/10)/(1/3) = 3/10
# j=0 contribution: a₂·(1-α)²·f_spin(α) = a₂·(1-α)²·α
#   → ∫₀¹ 2α²(1-α)² dα = ... wait, this doesn't have (1-α)² as I wrote.

# I need to be more careful. Let me redo this.

# In the proper-time representation of the vertex diagram:
# There are two electron propagators with proper times t₁ and t₂,
# and one photon propagator with proper time t₃.
# After Schwinger parameterization: t₁ = t·α₁, t₂ = t·α₂, t₃ = t·α₃
# with α₁+α₂+α₃ = 1.

# The Schwinger result involves integrating over (α₁, α₂, α₃) with
# the constraint α₁+α₂+α₃=1.

# For the anomalous magnetic moment:
# F₂ = (α/π)∫dα₁dα₂dα₃ δ(1-α₁-α₂-α₃) × m²α₃(1-α₃) / D²
# (this is the standard textbook form)

# The curvature corrections come from expanding the electron propagators:
# K(t_i) = K₀(t_i)[1 + a₁·t_i + a₂·t_i² + ...]

# For the first propagator: t₁ = t·α₁, so a₁·t₁ = a₁·t·α₁
# For the second: t₂ = t·α₂, so a₁·t₂ = a₁·t·α₂

# After the t-integration (which gives 1/λ² per power of t):
# c₁/λ² comes from a₁·(average α₁ or α₂ over the vertex integral)
# c₂/λ⁴ comes from:
#   a₂·(average α₁² or α₂²) + a₁²·(average α₁·α₂)

# But which averages enter F₂?
# The magnetic form factor projects out a specific spin structure.
# In the Schwinger calculation, the spin structure gives:
# F₂ ∝ ∫ dα₁dα₂dα₃ δ(Σ-1) α₃(1-α₃)/D²
# where α₃ is the photon fraction and D = α₁α₂ + (α₁+α₂)α₃·(m²/p²) +...

# In the limit of zero external momentum (which is what we have on S³):
# F₂ ∝ ∫ dα₁dα₂dα₃ δ(Σ-1) α₃(1-α₃)/[α₁α₂·m² + ...]²

# The curvature corrections enter through the ELECTRON propagator mass:
# m² → m² + R/12 (Parker-Toms leading correction)
# → m² + R/12 + R²·c₂_prop/m²
# where c₂_prop is the second-order propagator coefficient.

# On S³, m² corresponds to λ², and the corrections are:
# λ² → λ² + R/6 (from a₁ heat kernel) + a₂ (from a₂ heat kernel) / λ²

# The vertex integral for F₂ at zero external momentum reduces to:
# F₂/Schwinger = 1 + <curvature correction>/<vertex norm>

# For c₁: δm₁² = a₁ for propagator 1, δm₂² = a₁ for propagator 2
# F₂ involves the integral weighted by α₃(1-α₃)/D², where D ∝ α₁α₂
# The correction δ(1/D²) from δm² gives -2δm²/D³ × (∂D/∂m²)
# In the zero-momentum limit: ∂D/∂m² ∝ (α₁+α₂)
# So the correction is -2a₁(α₁+α₂)/D³ × D

# This is getting quite involved. Let me switch to a direct numerical
# approach: compute the spectral sum at high n_ext and extract c₂.

# DIRECT COMPUTATION: compute the vertex correction at n_ext = 5,10,15,20
# using the shell-summed spectral mode sum, then fit the curvature expansion.

print("\n  Computing vertex correction spectral sums at multiple n_ext...")
print("  Using qed_anomalous_moment module for m_j-resolved computation.")
print("  (This is the only way to extract F₂ correctly.)")

# Use the anomalous moment module
from geovac.qed_anomalous_moment import compute_anomalous_magnetic_moment

# Compute at n_ext = 1,2,3,4,5 with sufficient n_int for convergence
n_ext_list = [1, 2, 3, 4, 5]
n_int_max_list = [25, 15, 12, 10, 8]  # fewer needed at higher n_ext

print(f"\n  n_ext values: {n_ext_list}")
print(f"  n_int_max:    {n_int_max_list}")

# Try loading existing data first
results_file = 'debug/data/g2_c2_analytical.json'
try:
    with open(results_file) as f:
        saved_results = json.load(f)
    print(f"  Loaded existing results from {results_file}")
    USE_SAVED = True
except FileNotFoundError:
    saved_results = {}
    USE_SAVED = False

# Compute anomalous moment at each n_ext
am_results = {}
for i, n_ext in enumerate(n_ext_list):
    key = str(n_ext)
    if USE_SAVED and key in saved_results.get('am_results', {}):
        am_results[key] = saved_results['am_results'][key]
        print(f"  n_ext={n_ext}: loaded from cache (F₂/S = {am_results[key]['ratio']:.10f})")
        continue

    n_int = n_int_max_list[i]
    print(f"  Computing n_ext={n_ext}, n_int_max={n_int}...", end='', flush=True)
    import time
    t0 = time.time()

    try:
        result = compute_anomalous_magnetic_moment(n_ext, n_int)
        elapsed = time.time() - t0
        ratio = result['F2_over_schwinger']
        print(f"  done ({elapsed:.1f}s)  F₂/S = {ratio:.10f}")

        am_results[key] = {
            'n_ext': n_ext,
            'n_int_max': n_int,
            'F2_float': result['F2'],
            'F2_over_schwinger': ratio,
            'B_up': result['L_up'],
            'B_dn': result['L_dn'],
            'B_diff': result['B'],
            'V_mag': result['V_magnetic'],
            'elapsed_s': elapsed,
        }
    except Exception as e:
        elapsed = time.time() - t0
        print(f"  FAILED ({elapsed:.1f}s): {e}")
        am_results[key] = {'error': str(e), 'n_ext': n_ext}

# Extract curvature coefficients
print("\n" + "=" * 75)
print("  SECTION 9: Curvature coefficient extraction")
print("=" * 75)

good_data = []
for key in sorted(am_results.keys(), key=int):
    r = am_results[key]
    if 'error' not in r:
        n_ext = r['n_ext']
        lam = (2*n_ext + 3) / 2.0
        ratio = r['F2_over_schwinger']
        delta = ratio - 1.0
        good_data.append((n_ext, lam, delta, ratio))

print(f"\n  Good data points: {len(good_data)}")
for n_ext, lam, delta, ratio in good_data:
    x = 1.0/lam**2
    print(f"    n_ext={n_ext}: λ={lam:.1f}  1/λ²={x:.6f}  "
          f"F₂/S={ratio:.10f}  delta={delta:.10f}")

if len(good_data) >= 3:
    # Fit: delta = c₁·x + c₂·x² + c₃·x³ + ...
    # where x = 1/λ²
    import numpy as np

    n_pts = len(good_data)
    x_vals = np.array([1.0/lam**2 for _, lam, _, _ in good_data])
    delta_vals = np.array([d for _, _, d, _ in good_data])

    for degree in range(1, min(n_pts, 5)):
        # Build Vandermonde: delta = c₁·x + c₂·x² + ... (no constant term)
        V = np.zeros((n_pts, degree))
        for i in range(n_pts):
            for j in range(degree):
                V[i, j] = x_vals[i]**(j+1)

        # Least squares fit
        coeffs, residuals, rank, sv = np.linalg.lstsq(V, delta_vals, rcond=None)
        pred = V @ coeffs
        max_err = np.max(np.abs(delta_vals - pred))

        labels = ['c₁', 'c₂', 'c₃', 'c₄', 'c₅']
        print(f"\n  Degree-{degree} fit (x = 1/λ²):")
        for j in range(degree):
            print(f"    {labels[j]} = {coeffs[j]:.10f}")
        print(f"    Max residual: {max_err:.2e}")

        if degree >= 2:
            c2_fitted = coeffs[1]
            print(f"    c₂ = {c2_fitted:.10f}")

            # Check if c₁ ≈ 1/2
            c1_fitted = coeffs[0]
            print(f"    c₁ = {c1_fitted:.10f}  (expected 1/2 = 0.5)")
            print(f"    c₁ error: {abs(c1_fitted - 0.5):.2e}")

    # CONSTRAINED fit: fix c₁ = 1/2 exactly
    print("\n  Constrained fit (c₁ = 1/2 fixed):")
    delta_adjusted = delta_vals - 0.5 * x_vals  # subtract known c₁ term

    for degree in range(1, min(n_pts, 4)):
        V = np.zeros((n_pts, degree))
        for i in range(n_pts):
            for j in range(degree):
                V[i, j] = x_vals[i]**(j+2)  # x², x³, x⁴, ...

        coeffs, _, _, _ = np.linalg.lstsq(V, delta_adjusted, rcond=None)
        pred = V @ coeffs
        max_err = np.max(np.abs(delta_adjusted - pred))

        labels = ['c₂', 'c₃', 'c₄', 'c₅']
        print(f"\n  Fixed-c₁ degree-{degree} fit:")
        for j in range(degree):
            print(f"    {labels[j]} = {coeffs[j]:.10f}")
        print(f"    Max residual: {max_err:.2e}")

    # Best c₂ estimate from constrained fit
    V2 = np.zeros((n_pts, min(n_pts-1, 3)))
    for i in range(n_pts):
        for j in range(V2.shape[1]):
            V2[i, j] = x_vals[i]**(j+2)
    c_best, _, _, _ = np.linalg.lstsq(V2, delta_adjusted, rcond=None)
    c2_best = c_best[0]

    print(f"\n  BEST c₂ estimate: {c2_best:.10f}")

    # Rational identification of c₂
    print(f"\n  Rational identification of c₂ = {c2_best:.10f}:")

    # Try PSLQ
    c2_mp_best = mpmath.mpf(str(c2_best))

    # Simple rationals with small denominators
    best_rats = []
    for d in range(1, 500):
        n = round(float(c2_mp_best) * d)
        if n > 0:
            rat_val = n / d
            err = abs(rat_val - float(c2_mp_best))
            if err < 0.001:
                best_rats.append((n, d, err))

    best_rats.sort(key=lambda x: x[2])
    print(f"  Best rational approximants:")
    for n, d, err in best_rats[:10]:
        from math import gcd
        g = gcd(n, d)
        print(f"    {n//g}/{d//g} = {n/d:.10f}  err = {err:.2e}")

    # Check against heat kernel candidates
    print(f"\n  Heat kernel candidates:")
    hk_candidates = [
        ("c₁² = 1/4", 0.25),
        ("3/16", 3.0/16),
        ("7/40", 7.0/40),
        ("9/40", 9.0/40),
        ("c₁²−c₁/12 = 1/4−1/24 = 5/24", 5.0/24),
        ("c₁²−1/12 = 2/12 = 1/6", 1.0/6),
        ("c₁(c₁−1/6) = (1/2)(1/3) = 1/6", 1.0/6),
        ("(c₁²+c₁/6)/2 = (1/4+1/12)/2 = 1/6", 1.0/6),
        ("5/36", 5.0/36),
        ("1/6", 1.0/6),
        ("7/36", 7.0/36),
        ("7/48", 7.0/48),
        ("17/96", 17.0/96),
        ("35/192", 35.0/192),
        ("(R/12)²·(1+R/30) = c₁²·(1+1/5)", 0.25*1.2),
        ("(R/12)²·(2R/45+1) = ?", 0.25*(12.0/45+1)),
    ]

    for label, val in sorted(hk_candidates, key=lambda x: abs(x[1]-c2_best)):
        err = abs(val - c2_best)
        pct = err / abs(c2_best) * 100
        mark = " <-- CLOSE" if pct < 2 else ""
        print(f"    {label:45s} = {val:.10f}  err={err:.2e} ({pct:.2f}%){mark}")


# =====================================================================
# SECTION 10: Save results
# =====================================================================

print("\n" + "=" * 75)
print("  SECTION 10: Save results")
print("=" * 75)

output = {
    'geometry': {
        'd': 3,
        'R_scalar': 6,
        'Ric_sq': 12,
        'Riem_sq': 12,
        'Weyl_sq': 0,
        'E_lich': float(E_lich),
    },
    'seeley_dewitt': {
        'a1_per_comp': float(a1_per_comp) if isinstance(a1_per_comp, (int, float)) else float(sp.N(a1_per_comp)),
        'a2_density_per_comp': float(a2_density_per_comp) if isinstance(a2_density_per_comp, (int, float)) else float(sp.N(a2_density_per_comp)),
        'Omega_sq_trace_4comp': float(Omega_sq_trace_4comp),
    },
    'known_coefficients': {
        'c1_exact': '1/2',
        'c1_float': 0.5,
        'c1_source': 'Parker-Toms: R/12 = 6/12 for Dirac field on unit S^3',
    },
    'am_results': am_results,
}

if len(good_data) >= 3:
    output['c2_extraction'] = {
        'c2_best': float(c2_best),
        'c2_apparent_nex1': float(c2_apparent) if 'c2_apparent' in dir() else None,
        'method': 'constrained polynomial fit with c1=1/2 fixed',
        'n_ext_values': [g[0] for g in good_data],
        'lambda_values': [g[1] for g in good_data],
        'delta_values': [g[2] for g in good_data],
    }

    # Check 7/40 specifically
    c2_7_40 = 7.0/40
    output['c2_candidates'] = {
        '7/40': {'value': 0.175, 'error': abs(c2_best - 0.175)},
        '9/40': {'value': 0.225, 'error': abs(c2_best - 0.225)},
        '3/16': {'value': 0.1875, 'error': abs(c2_best - 0.1875)},
        '1/6': {'value': 1/6, 'error': abs(c2_best - 1/6)},
        '7/36': {'value': 7/36, 'error': abs(c2_best - 7/36)},
        '5/24': {'value': 5/24, 'error': abs(c2_best - 5/24)},
    }

with open(results_file, 'w') as f:
    json.dump(output, f, indent=2, default=str)
print(f"  Results saved to {results_file}")

print("\n" + "=" * 75)
print("  DONE")
print("=" * 75)
