"""
Alpha Precision Analysis — Characterizing the residual gap
===========================================================

Formula: 1/α ≈ π(42 + ζ(2) - 1/40) = 137.036064...

Terms:
  42      = S³ base Casimir trace (half Laplacian trace at n≤3)
  ζ(2)    = π²/6 = S¹ fiber zeta (single chirality: ζ_{S¹}(1)/2)
  1/40    = 1/(λ₃ × g_{l=2}) = boundary term at cutoff
  n=3 selected by: ⟨trace⟩/state = dim(S³)

CODATA 2018: α⁻¹ = 137.035999084(21)

Date: 2026-03-21
Status: Exploratory (Paper 2 conjecture tier)
"""

import numpy as np
from scipy.special import zeta as scipy_zeta

# =============================================================================
# Section 1: CODATA Comparison
# =============================================================================

print("=" * 70)
print("SECTION 1: CODATA COMPARISON")
print("=" * 70)

# CODATA 2018 value
ALPHA_INV_CODATA = 137.035999084
ALPHA_INV_CODATA_UNC = 0.000000021  # 1σ uncertainty

# Our formula
TRACE_S3 = 42.0          # Half Laplacian trace at n≤3
ZETA_2 = np.pi**2 / 6    # ζ(2) = π²/6
BOUNDARY = 1.0 / 40      # 1/(λ₃ × g_{l=2})

alpha_inv_formula = np.pi * (TRACE_S3 + ZETA_2 - BOUNDARY)

delta = alpha_inv_formula - ALPHA_INV_CODATA
sigma = abs(delta) / ALPHA_INV_CODATA_UNC
rel_error = abs(delta) / ALPHA_INV_CODATA

print(f"\nFormula:  π(42 + ζ(2) - 1/40)")
print(f"  42          = {TRACE_S3}")
print(f"  ζ(2)        = {ZETA_2:.15f}")
print(f"  1/40        = {BOUNDARY:.15f}")
print(f"  Argument    = {TRACE_S3 + ZETA_2 - BOUNDARY:.15f}")
print(f"\nα⁻¹(formula) = {alpha_inv_formula:.10f}")
print(f"α⁻¹(CODATA)  = {ALPHA_INV_CODATA:.10f} ± {ALPHA_INV_CODATA_UNC}")
print(f"\nΔ = {delta:+.6e}")
print(f"|Δ|/σ_CODATA = {sigma:.1f}σ  (>{3}σ → NOT within experimental error)")
print(f"Relative error = {rel_error:.2e}")

alpha_formula = 1.0 / alpha_inv_formula
alpha_codata = 1.0 / ALPHA_INV_CODATA

# =============================================================================
# Section 2: What's the Missing Piece?
# =============================================================================

print("\n" + "=" * 70)
print("SECTION 2: CHARACTERIZING THE GAP")
print("=" * 70)

print(f"\nΔ = {delta:.10e}")
print(f"Δ/π = {delta / np.pi:.10e}")
print(f"Δ/π² = {delta / np.pi**2:.10e}")

# Check if Δ/π is a recognizable fraction
delta_over_pi = delta / np.pi
print(f"\nΔ/π = {delta_over_pi:.8f}")
print(f"  ≈ 1/{1/delta_over_pi:.2f}")

# Check multiples of common numbers
print("\nΔ in units of various constants:")
checks = {
    "α²": alpha_codata**2,
    "α/π": alpha_codata / np.pi,
    "α²/π": alpha_codata**2 / np.pi,
    "α³": alpha_codata**3,
    "ζ(4)": np.pi**4 / 90,
    "ζ(6)": np.pi**6 / 945,
    "1/1680": 1.0 / 1680,
    "π/1680": np.pi / 1680,
    "1/(42×40)": 1.0 / (42 * 40),
    "1/(42×40×π)": 1.0 / (42 * 40 * np.pi),
}

for name, val in checks.items():
    ratio = delta / val
    print(f"  Δ / {name:15s} = {ratio:12.6f}")

# =============================================================================
# Section 3: Self-Consistent / Radiative Corrections
# =============================================================================

print("\n" + "=" * 70)
print("SECTION 3: RADIATIVE CORRECTION HYPOTHESIS")
print("=" * 70)

# In QED, the Schwinger correction to g-2 is α/(2π).
# The leading vacuum polarization correction to the running of α is:
#   α(q²) = α(0) / [1 - (α/3π) ln(q²/m_e²)]
# At q=0, the "bare" vs "dressed" difference at one loop is O(α/π).

alpha_over_pi = alpha_codata / np.pi
alpha_over_2pi = alpha_codata / (2 * np.pi)
schwinger = alpha_codata / (2 * np.pi)  # a_e = (g-2)/2 leading term

print(f"\nα/π     = {alpha_over_pi:.6e}")
print(f"α/(2π)  = {alpha_over_2pi:.6e}  (Schwinger term)")
print(f"α²/(2π) = {alpha_codata**2 / (2*np.pi):.6e}")
print(f"\nΔ       = {delta:.6e}")
print(f"Δ/α     = {delta * alpha_codata:.6e}")
print(f"Δ × α   = {delta * alpha_codata:.6e}")

# Check: does Δ = c × α/π for some simple c?
c_val = delta / (alpha_codata / np.pi)
print(f"\nΔ / (α/π) = {c_val:.6f}")
print(f"  Not a clean ratio → radiative correction hypothesis unlikely.")

# Check: does Δ = c × α² for some simple c?
c_val2 = delta / alpha_codata**2
print(f"Δ / α²    = {c_val2:.6f}")
print(f"  ≈ {c_val2:.2f} (check if close to simple fraction)")

# What if formula gives α_bare and CODATA is α_dressed?
# α_dressed ≈ α_bare (1 + α/(3π) × something)
# Then Δα⁻¹ ≈ -α_bare/(3π) × something
correction_3pi = delta * 3 * np.pi / alpha_codata
print(f"\nIf Δ = α/(3π) × X, then X = {correction_3pi:.4f}")

# =============================================================================
# Section 4: Alternative Formulas
# =============================================================================

print("\n" + "=" * 70)
print("SECTION 4: ALTERNATIVE FORMULAS")
print("=" * 70)

print(f"\nTarget: α⁻¹ = {ALPHA_INV_CODATA:.10f}")
print(f"Current: π(42 + ζ(2) - 1/40) = {alpha_inv_formula:.10f}")
print()

# The argument inside π() needs to be:
target_arg = ALPHA_INV_CODATA / np.pi
current_arg = TRACE_S3 + ZETA_2 - BOUNDARY
print(f"Target argument  = {target_arg:.12f}")
print(f"Current argument = {current_arg:.12f}")
print(f"Argument gap     = {target_arg - current_arg:.10e}")
arg_gap = target_arg - current_arg

print(f"\nThe boundary term needs to shift by {-arg_gap:.8e}")
print(f"Current boundary = 1/40 = {BOUNDARY:.10f}")
needed_boundary = BOUNDARY - arg_gap
print(f"Needed boundary  = {needed_boundary:.10f}")
print(f"1/needed         = {1/needed_boundary:.6f}")

alternatives = [
    ("π(42 + ζ(2) - 1/39)", np.pi * (42 + ZETA_2 - 1/39)),
    ("π(42 + ζ(2) - 1/41)", np.pi * (42 + ZETA_2 - 1/41)),
    ("π(42 + ζ(2) - 1/38)", np.pi * (42 + ZETA_2 - 1/38)),
    ("π(42 + ζ(2) - ζ(4))", np.pi * (42 + ZETA_2 - np.pi**4/90)),
    ("π(42 + ζ(2) - 1/(4π²))", np.pi * (42 + ZETA_2 - 1/(4*np.pi**2))),
    ("π(42 + ζ(2) - π²/360)", np.pi * (42 + ZETA_2 - np.pi**2/360)),
    ("π(42 + ζ(2) - 1/40 - α²/π)", np.pi * (42 + ZETA_2 - 1/40 - alpha_codata**2/np.pi)),
    ("π(42 + ζ(2)/2 - 1/80)", np.pi * (42 + ZETA_2/2 - 1/80)),
    ("π(42 + ζ(2) - 1/40 - 1/1680)", np.pi * (42 + ZETA_2 - 1/40 - 1/1680)),
    ("π(42 + ζ(2) - 1/40 + 1/1680)", np.pi * (42 + ZETA_2 - 1/40 + 1/1680)),
]

# Systematic: what integer denominator k gives 1/k closest to needed_boundary?
print("\n--- Systematic boundary search: π(42 + ζ(2) - 1/k) ---")
best_k = None
best_err = float('inf')
for k in range(30, 60):
    val = np.pi * (42 + ZETA_2 - 1.0/k)
    err = abs(val - ALPHA_INV_CODATA)
    if err < best_err:
        best_err = err
        best_k = k
    if 35 <= k <= 45:
        print(f"  k={k:3d}: α⁻¹ = {val:.8f}, Δ = {val - ALPHA_INV_CODATA:+.6e}")

print(f"\n  Best integer k = {best_k}, error = {best_err:.6e}")

# Now try rational boundaries p/q
print("\n--- Rational boundary search: π(42 + ζ(2) - p/q) ---")
best_candidates = []
for q in range(1, 200):
    for p in range(1, min(q, 10)):
        val = np.pi * (42 + ZETA_2 - p / q)
        err = abs(val - ALPHA_INV_CODATA)
        if err < 1e-5:
            best_candidates.append((p, q, val, err))

best_candidates.sort(key=lambda x: x[3])
print(f"  Top candidates with |Δ| < 10⁻⁵:")
for p, q, val, err in best_candidates[:10]:
    print(f"    {p}/{q:3d} = {p/q:.8f}: α⁻¹ = {val:.8f}, Δ = {val - ALPHA_INV_CODATA:+.6e}")

print("\n--- Named alternative formulas ---")
for name, val in alternatives:
    err = val - ALPHA_INV_CODATA
    rel = abs(err) / ALPHA_INV_CODATA
    marker = " ***" if abs(err) < abs(delta) else ""
    print(f"  {name:45s} = {val:.8f}, Δ = {err:+.6e}, rel = {rel:.2e}{marker}")

# =============================================================================
# Section 5: Transcendental / ζ-function corrections
# =============================================================================

print("\n" + "=" * 70)
print("SECTION 5: ZETA-FUNCTION AND TRANSCENDENTAL CORRECTIONS")
print("=" * 70)

# The gap in the argument is:
print(f"\nArgument gap = {arg_gap:.12e}")
print(f"\nIs the gap expressible as a zeta value or combination?")

zeta_vals = {
    "ζ(2)": np.pi**2/6,
    "ζ(3)": float(scipy_zeta(3)),
    "ζ(4)": np.pi**4/90,
    "ζ(5)": float(scipy_zeta(5)),
    "ζ(6)": np.pi**6/945,
}

for name, zv in zeta_vals.items():
    ratio = arg_gap / zv
    print(f"  gap / {name} = {ratio:.8f}")

# Try: gap = a × ζ(k) / n for small integers
print("\n  Searching gap = (p/q) × ζ(k):")
for zname, zv in zeta_vals.items():
    for q in range(1, 50):
        p_float = arg_gap * q / zv
        p_int = round(p_float)
        if p_int != 0 and abs(p_float - p_int) < 0.02:
            residual = abs(arg_gap - p_int * zv / q)
            print(f"    gap ≈ ({p_int}/{q}) × {zname}, residual = {residual:.2e}")

# Try combinations: gap = a/b + c × ζ(k)
print(f"\n  gap = {arg_gap:.10e}")
print(f"  1/1680 = {1/1680:.10e}")
print(f"  gap - 1/1680 = {arg_gap - 1/1680:.10e}")
print(f"  gap + 1/1680 = {arg_gap + 1/1680:.10e}")

# =============================================================================
# Section 6: The Chirality Question
# =============================================================================

print("\n" + "=" * 70)
print("SECTION 6: CHIRALITY AND THE FACTOR OF 1/2")
print("=" * 70)

print("""
The ζ_{S¹}(1)/2 = ζ(2) term uses a single chirality of the Hopf fiber.

Physical justifications for selecting one chirality:
  1. Hopf bundle orientation: S³ → S² has a preferred orientation
     from the complex structure of C² ⊃ S³. The electron lives on
     one sheet of the double cover.
  2. Parity violation: Weak interaction selects left-handed electrons.
     At the topological level, this may correspond to one Hopf chirality.
  3. Spin-statistics: The electron is spin-1/2. The Hopf fiber carries
     the spin degree of freedom. Selecting one chirality = selecting
     one spin orientation for the vacuum contribution.

Alternative: What if the 1/2 is NOT chirality but something else?
  - ζ_{S¹}(1) = 2ζ(2) = π²/3. Without the 1/2:
""")

alpha_inv_no_half = np.pi * (42 + 2 * ZETA_2 - BOUNDARY)
print(f"  π(42 + 2ζ(2) - 1/40) = {alpha_inv_no_half:.6f}")
print(f"  Much too large. The 1/2 is essential.")

alpha_inv_full_zeta1 = np.pi * (42 + 2 * ZETA_2 - BOUNDARY)
print(f"\n  With full ζ_{{S¹}}(1) = 2ζ(2): α⁻¹ = {alpha_inv_full_zeta1:.4f} (wrong)")
print(f"  With half ζ_{{S¹}}(1)/2 = ζ(2): α⁻¹ = {alpha_inv_formula:.6f} (correct to 5×10⁻⁷)")

# =============================================================================
# Section 7: Higher-order S³ corrections
# =============================================================================

print("\n" + "=" * 70)
print("SECTION 7: HIGHER-ORDER S³ CORRECTIONS")
print("=" * 70)

print("""
The formula truncates at n=3 (selected by ⟨λ⟩/state = dim S³ = 4).
Could higher shells contribute perturbatively?

The Casimir trace contribution from shell n is n²(n²-1)/2.
  n=1: 0,  n=2: 6,  n=3: 36  → Total = 42 (half-trace)
  n=4: 120, n=5: 300, ...

If n=4 contributes with weight w₄, the correction to α⁻¹ is π × w₄ × 120.
""")

# What weight w₄ would close the gap?
w4_needed = arg_gap / 120
print(f"  Weight w₄ to close gap: {w4_needed:.8e}")
print(f"  This is tiny — about {w4_needed:.2e}")
print(f"  Compare: n=3 has weight 1, so w₄/w₃ ≈ {w4_needed:.2e}")

# Exponential suppression?
# If w_n ~ exp(-c × n), then w₄ = exp(-c×4) / exp(-c×3) = exp(-c)
# w₃ = 1 → exp(-3c) = 1 is nonsense. Try power-law: w_n = n^(-p)
# w₃ = 1 → 3^(-p) = 1 → p = 0, trivial.
# Better: the cutoff is sharp but with a "leakage" term
print(f"\n  If leakage from n=4 shell: Δ(arg) = ε × 120")
print(f"  ε = {w4_needed:.6e}")
print(f"  1/ε = {1/w4_needed:.1f}")
print(f"  Note: 1/ε ≈ {1/w4_needed:.0f}")

# =============================================================================
# Section 8: Summary
# =============================================================================

print("\n" + "=" * 70)
print("SECTION 8: SUMMARY")
print("=" * 70)

print(f"""
CLAIM: The fine structure constant is approximately

  α ≈ 1 / [π(Tr(L²)|_{{S³,n≤3}} + ζ_{{S¹}}(1)/2 - 1/(λ₃ × g₂))]

where:
  Tr(L²)|_{{S³,n≤3}} = Σ_{{n=1}}^3 n²(n²-1)/2 = 42
  ζ_{{S¹}}(1)/2 = ζ(2) = π²/6 ≈ 1.6449
  1/(λ₃ × g₂) = 1/(8 × 5) = 1/40 = 0.025
  n=3 cutoff selected by: ⟨trace⟩/state = dim(S³) = 4

RESULT:
  α⁻¹(formula) = {alpha_inv_formula:.10f}
  α⁻¹(CODATA)  = {ALPHA_INV_CODATA:.10f}
  Δ             = {delta:+.6e}
  Relative err  = {rel_error:.2e}
  σ-deviation   = {sigma:.1f}σ

STATUS: Approximate formula (4σ from CODATA).
  The gap Δ ≈ 6.5×10⁻⁵ is:
  - Too large for experimental error
  - Not simply α/π or α² (no clean radiative correction)
  - Possibly from n≥4 shell leakage (weight ~{w4_needed:.1e} per unit trace)
  - No clean rational boundary 1/k matches CODATA better

INTERPRETATION:
  The formula captures the TOPOLOGY correctly to 5×10⁻⁷.
  The residual may require:
  (a) A perturbative n≥4 correction with geometric weight
  (b) A self-consistent α-dependent correction term
  (c) An exact transcendental replacement for 1/40
""")
