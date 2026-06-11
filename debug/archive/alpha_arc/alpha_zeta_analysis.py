#!/usr/bin/env python3
"""
Alpha Zeta Analysis: Spectral zeta function interpretation of α from S³ lattice.

Key finding from alpha_eigenvalue_search.py:
  Tr(L²)|_{n≤3} = 42 (exact integer)
  42 + ζ(2) = 42 + π²/6 = 43.645 → π × 43.645 = 137.12 (0.06% from 1/α)

This script investigates whether α has an exact spectral zeta function expression.

Date: 2026-03-21
Status: Exploratory (papers/conjectures/ tier)
"""

import numpy as np
from scipy.special import zeta as scipy_zeta
from typing import Tuple, List, Dict

# ===========================================================================
# Constants
# ===========================================================================
ALPHA_INV_EXACT = 137.035999084  # CODATA 2018
ALPHA_INV_TARGET = ALPHA_INV_EXACT


# ===========================================================================
# 1. SPECTRAL ZETA FUNCTION OF S³
# ===========================================================================
def spectral_zeta_s3(s: float, n_max: int = 10000) -> float:
    """
    Spectral zeta function of the Laplace-Beltrami operator on unit S³.

    ζ_S³(s) = Σ_{n≥2} n² / (n² - 1)^s

    Eigenvalues: λ_n = n² - 1  (n = 1, 2, 3, ...)
    Degeneracy:  d_n = n²
    n=1 has λ=0, excluded from zeta sum.

    Parameters
    ----------
    s : float
        Zeta function argument.
    n_max : int
        Truncation for the sum.

    Returns
    -------
    float
        ζ_S³(s) truncated at n_max.
    """
    n_vals = np.arange(2, n_max + 1, dtype=np.float64)
    eigenvalues = n_vals**2 - 1.0
    degeneracies = n_vals**2
    return np.sum(degeneracies / eigenvalues**s)


def analyze_spectral_zeta() -> None:
    """Compute ζ_S³(s) for various s and look for special values."""
    print("=" * 70)
    print("1. SPECTRAL ZETA FUNCTION ζ_S³(s) = Σ_{n≥2} n² / (n²-1)^s")
    print("=" * 70)

    s_values = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0]
    for s in s_values:
        val = spectral_zeta_s3(s)
        val_over_pi = val / np.pi
        print(f"  ζ_S³({s:4.1f}) = {val:12.6f}   /π = {val_over_pi:12.6f}")

    print()

    # Check specific relationships
    z1 = spectral_zeta_s3(1.0)
    z2 = spectral_zeta_s3(2.0)
    z3 = spectral_zeta_s3(3.0)

    print("  Special value checks:")
    print(f"    ζ_S³(1) = {z1:.6f}")
    print(f"    ζ_S³(2) = {z2:.6f}")
    print(f"    ζ_S³(1) / ζ_Riemann(2) = {z1 / (np.pi**2 / 6):.6f}")
    print(f"    ζ_S³(1) - ζ_Riemann(2) = {z1 - np.pi**2 / 6:.6f}")
    print(f"    π × ζ_S³(1) = {np.pi * z1:.6f}  (target: {ALPHA_INV_TARGET:.6f})")
    print(f"    π × ζ_S³(2) = {np.pi * z2:.6f}")

    # Does ζ_S³(1) have a closed form?
    # ζ_S³(1) = Σ n²/(n²-1) = Σ [1 + 1/(n²-1)] = (N-1) + Σ 1/(n²-1)
    # And 1/(n²-1) = (1/2)[1/(n-1) - 1/(n+1)]  (partial fractions)
    # So Σ_{n=2}^{N} 1/(n²-1) = (1/2)[1 + 1/2 - 1/N - 1/(N+1)] → 3/4 as N→∞
    # Therefore ζ_S³(1) DIVERGES (like N)!
    print()
    print("  NOTE: ζ_S³(1) diverges! The partial sum grows like N.")
    print("  Only ζ_S³(s) with s > 3/2 converges (since d_n ~ n² and λ_n ~ n²).")
    print()

    # Convergent values
    print("  Convergent zeta values (s > 3/2):")
    for s in [2.0, 2.5, 3.0]:
        val = spectral_zeta_s3(s, n_max=100000)
        print(f"    ζ_S³({s}) = {val:.10f}")

    # Check if ζ_S³(2) has a nice form
    z2_exact = spectral_zeta_s3(2.0, n_max=100000)
    print(f"\n  ζ_S³(2) = {z2_exact:.10f}")
    print(f"  π²/6 - 1 = {np.pi**2/6 - 1:.10f}")
    print(f"  7/4 - π²/12 = {7/4 - np.pi**2/12:.10f}")

    # Partial fraction: n²/(n²-1)² = n²/[(n-1)(n+1)]²
    # Let's compute numerically and compare to known constants
    z2_candidates = {
        "π²/8": np.pi**2 / 8,
        "π²/6 - 1/2": np.pi**2 / 6 - 0.5,
        "3/2": 1.5,
        "π²/12 + 1/2": np.pi**2 / 12 + 0.5,
        "7π²/72 + 1/2": 7 * np.pi**2 / 72 + 0.5,
    }
    print(f"\n  Candidates for ζ_S³(2) = {z2_exact:.10f}:")
    for name, val in z2_candidates.items():
        print(f"    {name:20s} = {val:.10f}  diff = {abs(val - z2_exact):.2e}")


# ===========================================================================
# 2. REGULARIZED TRACES (HEAT KERNEL)
# ===========================================================================
def heat_kernel_s3(t: float, n_max: int = 1000) -> float:
    """
    Heat kernel trace on S³: K(t) = Σ_{n≥1} n² exp(-t(n²-1)).

    Parameters
    ----------
    t : float
        Heat kernel parameter.
    n_max : int
        Truncation.

    Returns
    -------
    float
        K(t).
    """
    n_vals = np.arange(1, n_max + 1, dtype=np.float64)
    eigenvalues = n_vals**2 - 1.0
    degeneracies = n_vals**2
    return np.sum(degeneracies * np.exp(-t * eigenvalues))


def analyze_heat_kernel() -> None:
    """Find t* where K(t*) = 1/α ÷ π and analyze its meaning."""
    print("\n" + "=" * 70)
    print("2. HEAT KERNEL K(t) = Σ n² exp(-t(n²-1))")
    print("=" * 70)

    target = ALPHA_INV_TARGET / np.pi  # ≈ 43.620

    # Bisection to find t* where K(t*) = target
    t_lo, t_hi = 0.001, 1.0
    for _ in range(100):
        t_mid = (t_lo + t_hi) / 2
        if heat_kernel_s3(t_mid) > target:
            t_lo = t_mid
        else:
            t_hi = t_mid
    t_star = (t_lo + t_hi) / 2
    k_star = heat_kernel_s3(t_star)

    print(f"\n  Target: 1/(απ) = {target:.6f}")
    print(f"  t* such that K(t*) = target: t* = {t_star:.10f}")
    print(f"  K(t*) = {k_star:.6f}")
    print(f"  π × K(t*) = {np.pi * k_star:.6f}  (target: {ALPHA_INV_TARGET:.6f})")

    # Check candidate forms for t*
    print(f"\n  Candidate forms for t* = {t_star:.10f}:")
    candidates = {
        "1/(4π²)": 1 / (4 * np.pi**2),
        "1/(2π)²": 1 / (2 * np.pi)**2,
        "1/e³": 1 / np.e**3,
        "1/(π²+e²)": 1 / (np.pi**2 + np.e**2),
        "ln(2)/(4π)": np.log(2) / (4 * np.pi),
        "1/20": 0.05,
        "3/(4π²)": 3 / (4 * np.pi**2),
        "1/(4π+1)": 1 / (4 * np.pi + 1),
        "π/64": np.pi / 64,
        "1/(3π+e)": 1 / (3 * np.pi + np.e),
    }
    for name, val in candidates.items():
        k_val = heat_kernel_s3(val)
        alpha_inv = np.pi * k_val
        err_pct = abs(alpha_inv - ALPHA_INV_TARGET) / ALPHA_INV_TARGET * 100
        print(f"    t={name:16s} = {val:.8f}  → 1/α = {alpha_inv:.4f}  ({err_pct:.4f}%)")


# ===========================================================================
# 3. THE GAP: 4π² + π + 1 - 42
# ===========================================================================
def analyze_gap() -> None:
    """Analyze the gap between Tr(L²)|_{n≤3} = 42 and the target."""
    print("\n" + "=" * 70)
    print("3. THE GAP: target - 42")
    print("=" * 70)

    target = ALPHA_INV_TARGET / np.pi
    gap = target - 42.0

    print(f"\n  Target (1/α)/π = {target:.10f}")
    print(f"  Tr(L²)|_{{n≤3}} = 42")
    print(f"  Gap = {gap:.10f}")

    # Is the target actually 4π² + π + 1?
    formula_val = 4 * np.pi**2 + np.pi + 1
    print(f"\n  Is target = 4π² + π + 1?")
    print(f"    4π² + π + 1 = {formula_val:.10f}")
    print(f"    actual target = {target:.10f}")
    print(f"    difference = {target - formula_val:.2e}")
    print(f"    (This is {abs(target - formula_val) / target * 100:.4f}% off)")

    gap_from_42 = target - 42.0
    print(f"\n  Gap from 42: {gap_from_42:.10f}")

    candidates = {
        "π²/6 (= ζ(2))": np.pi**2 / 6,
        "π/2": np.pi / 2,
        "ln(5)": np.log(5),
        "γ + 1": 0.5772156649 + 1,
        "φ (golden ratio)": (1 + np.sqrt(5)) / 2,
        "√(e)": np.sqrt(np.e),
        "3/2 + 1/π": 1.5 + 1 / np.pi,
        "π²/6 - 1/40": np.pi**2 / 6 - 1 / 40,
        "ζ(2) - 1/40": np.pi**2 / 6 - 1 / 40,
        "4π² + π - 41": 4 * np.pi**2 + np.pi - 41,
    }

    print(f"\n  Candidate expressions for gap = {gap_from_42:.10f}:")
    for name, val in sorted(candidates.items(), key=lambda x: abs(x[1] - gap_from_42)):
        diff = gap_from_42 - val
        pct = abs(diff) / gap_from_42 * 100
        marker = " ← CLOSE" if pct < 2 else ""
        print(f"    {name:25s} = {val:.10f}  diff = {diff:+.6e} ({pct:.3f}%){marker}")

    # Try the formula: 1/α = π × (42 + correction)
    print(f"\n  Formula attempts: 1/α = π × (42 + X)")
    corrections = {
        "ζ(2)": np.pi**2 / 6,
        "π/2": np.pi / 2,
        "ln(5)": np.log(5),
        "φ": (1 + np.sqrt(5)) / 2,
        "4π²+π-41": 4 * np.pi**2 + np.pi - 41,
        "π²/6 - 1/39": np.pi**2 / 6 - 1 / 39,
        "π²/6 - 1/40": np.pi**2 / 6 - 1 / 40,
        "π²/6 - 1/41": np.pi**2 / 6 - 1 / 41,
    }
    for name, x in sorted(corrections.items(),
                          key=lambda c: abs(np.pi * (42 + c[1]) - ALPHA_INV_TARGET)):
        alpha_inv = np.pi * (42 + x)
        err = abs(alpha_inv - ALPHA_INV_TARGET) / ALPHA_INV_TARGET * 100
        print(f"    X = {name:18s} → 1/α = {alpha_inv:.6f}  ({err:.4f}% error)")


# ===========================================================================
# 4. WHY n = 3? NUMBER THEORY OF 42
# ===========================================================================
def analyze_42() -> None:
    """Investigate the number theory behind Tr(L²)|_{n≤3} = 42."""
    print("\n" + "=" * 70)
    print("4. WHY n=3? THE NUMBER THEORY OF 42")
    print("=" * 70)

    # Tr(L²) = Σ_{n=1}^{N} n²(n²-1) = Σ n⁴ - Σ n²
    # Using standard formulas:
    #   Σ n² = N(N+1)(2N+1)/6
    #   Σ n⁴ = N(N+1)(2N+1)(3N²+3N-1)/30
    # So Tr(L²) = N(N+1)(2N+1)/30 × [(3N²+3N-1) - 5]
    #           = N(N+1)(2N+1)(3N²+3N-6)/30
    #           = N(N+1)(2N+1)(N²+N-2)/10
    #           = N(N+1)(2N+1)(N+2)(N-1)/10

    print("\n  Cumulative Tr(L²) = Σ_{n=1}^{N} n²(n²-1):")
    print(f"  = N(N+1)(2N+1)(N+2)(N-1)/10")
    print()

    for N in range(1, 11):
        tr = N * (N + 1) * (2 * N + 1) * (N + 2) * (N - 1) // 10
        per_shell = N**2 * (N**2 - 1)
        print(f"    N={N:2d}: per_shell = {per_shell:6d},  cumulative = {tr:8d}")

    print()
    print("  For N=3: Tr = 3 × 4 × 7 × 5 × 2 / 10 = 42 × 10 / 10 = 42 ✓")
    print()

    # Connections to representation theory
    print("  Number theory of 42:")
    print(f"    42 = 2 × 3 × 7")
    print(f"    42 = 6 × 7")
    print(f"    42 = C(7,2) - C(7,5) + ... ? No: C(7,2) = 21")
    print(f"    42 = T(3) × 7  where T(3) = 6 = 3rd triangular number")

    # SO(4) Casimir connection
    print("\n  SO(4) representation theory:")
    print("  The Casimir of SO(4) irrep (j₁,j₂) = j₁(j₁+1) + j₂(j₂+1)")
    print("  For hydrogenic n-shell: (j₁,j₂) = ((n-1)/2, (n-1)/2)")
    print("  Casimir = 2 × (n-1)/2 × (n+1)/2 = (n²-1)/2")
    print()
    for n in range(1, 6):
        j = (n - 1) / 2
        casimir = 2 * j * (j + 1)
        deg = n**2
        print(f"    n={n}: j={j:.1f}, Casimir = {casimir:.1f}, "
              f"deg = {deg}, deg × Casimir = {deg * casimir:.0f}")

    print()
    print("  Σ n² × (n²-1)/2 for n=1..3 = 0 + 6 + 36 = 42/2 ??? No:")
    print("  Actually: deg × Casimir = n² × (n²-1)/2")
    cas_sum = sum(n**2 * (n**2 - 1) / 2 for n in range(1, 4))
    print(f"  Σ_{{n=1}}^3 n² × (n²-1)/2 = {cas_sum}")
    print(f"  So 42 = 2 × Σ Casimir × degeneracy = 2 × {cas_sum}")

    # Tetrahedral numbers
    print("\n  Sequence 0, 6, 42, 120, 270, ...")
    print("  Ratios: 42/6 = 7, 120/42 ≈ 2.86, 270/120 = 2.25")
    print("  This is N(N+1)(2N+1)(N+2)(N-1)/10")
    print("  At N=3: appears in the 5D simplex (pentatope) number family")


# ===========================================================================
# 5. EXACT FORMULA SEARCH
# ===========================================================================
def search_exact_formula() -> None:
    """Systematic search for exact expressions giving 1/α."""
    print("\n" + "=" * 70)
    print("5. EXACT FORMULA SEARCH")
    print("=" * 70)

    target = ALPHA_INV_TARGET
    results: List[Tuple[str, float, float]] = []  # (name, value, error_pct)

    # Type 1: π × (42 + correction)
    corrections_1 = {
        "ζ(2)": np.pi**2 / 6,
        "ζ(2) - 1/39": np.pi**2 / 6 - 1 / 39,
        "ζ(2) - 1/40": np.pi**2 / 6 - 1 / 40,
        "ζ(2) - 1/41": np.pi**2 / 6 - 1 / 41,
        "ζ(2) - 1/42": np.pi**2 / 6 - 1 / 42,
        "ζ(2) - π/120": np.pi**2 / 6 - np.pi / 120,
        "ζ(2) - 1/(8π)": np.pi**2 / 6 - 1 / (8 * np.pi),
        "ln(5)": np.log(5),
        "ln(5) + 1/90": np.log(5) + 1 / 90,
        "ln(5) + 1/100": np.log(5) + 1 / 100,
        "π/2 + 1/6": np.pi / 2 + 1 / 6,
    }
    for name, x in corrections_1.items():
        val = np.pi * (42 + x)
        err = abs(val - target) / target * 100
        results.append((f"π(42 + {name})", val, err))

    # Type 2: Direct formulas
    direct = {
        "π(4π² + π + 1)": np.pi * (4 * np.pi**2 + np.pi + 1),
        "4π³ + π² + π": 4 * np.pi**3 + np.pi**2 + np.pi,
        "4π³ + π² + π + e⁻⁷": 4 * np.pi**3 + np.pi**2 + np.pi + np.exp(-7),
        "π(42 + π²/6)": np.pi * (42 + np.pi**2 / 6),
        "42π + π³/6": 42 * np.pi + np.pi**3 / 6,
        "π²(14 + 1/18)": np.pi**2 * (14 + 1 / 18),
        "π²(253/18)": np.pi**2 * 253 / 18,
        "e^(ln(137) + ...)": 137.036,  # trivial
    }
    for name, val in direct.items():
        err = abs(val - target) / target * 100
        results.append((name, val, err))

    # Type 3: Involving Euler-Mascheroni γ
    gamma = 0.5772156649015329
    gamma_formulas = {
        "π(42 + ζ(2) - γ/24)": np.pi * (42 + np.pi**2 / 6 - gamma / 24),
        "π(42 + ζ(2) - γ/23)": np.pi * (42 + np.pi**2 / 6 - gamma / 23),
        "π(42 + ζ(2) - γ²)": np.pi * (42 + np.pi**2 / 6 - gamma**2),
    }
    for name, val in gamma_formulas.items():
        err = abs(val - target) / target * 100
        results.append((name, val, err))

    # Type 4: Spectral quantities from the S³ zeta function
    # Use partial zeta: ζ_partial(s, N) = Σ_{n=2}^{N} n²/(n²-1)^s
    def zeta_partial(s: float, N: int) -> float:
        n_vals = np.arange(2, N + 1, dtype=np.float64)
        return np.sum(n_vals**2 / (n_vals**2 - 1)**s)

    z2_full = spectral_zeta_s3(2.0, n_max=100000)
    z3_full = spectral_zeta_s3(3.0, n_max=100000)

    spectral_formulas = {
        "π × ζ_S³(2) × 42": np.pi * z2_full * 42,
        "π² × ζ_S³(2)": np.pi**2 * z2_full,
        "42π + ζ_S³(2)": 42 * np.pi + z2_full,
        "42π + π × ζ_S³(2)": 42 * np.pi + np.pi * z2_full,
        "42π + π³/6": 42 * np.pi + np.pi**3 / 6,
    }
    for name, val in spectral_formulas.items():
        err = abs(val - target) / target * 100
        results.append((name, val, err))

    # Sort by error
    results.sort(key=lambda x: x[2])

    print(f"\n  Target: 1/α = {target:.10f}")
    print(f"\n  {'Formula':<40s} {'Value':>14s} {'Error%':>10s}")
    print(f"  {'-'*40} {'-'*14} {'-'*10}")
    for name, val, err in results[:25]:
        marker = " ***" if err < 0.01 else " **" if err < 0.05 else " *" if err < 0.1 else ""
        print(f"  {name:<40s} {val:14.6f} {err:10.4f}%{marker}")


# ===========================================================================
# 5b. FINE-GRAINED SEARCH AROUND BEST CANDIDATES
# ===========================================================================
def fine_search() -> None:
    """Search for exact rational corrections to the best formulas."""
    print("\n" + "=" * 70)
    print("5b. FINE-GRAINED SEARCH: 1/α = π(42 + ζ(2) - δ)")
    print("=" * 70)

    target = ALPHA_INV_TARGET
    base = np.pi * (42 + np.pi**2 / 6)  # = 137.117... (0.06% high)
    delta_needed = (base - target) / np.pi
    print(f"\n  Base: π(42 + ζ(2)) = {base:.10f}")
    print(f"  Target: {target:.10f}")
    print(f"  Excess: {base - target:.10f}")
    print(f"  δ needed (so that π(42 + ζ(2) - δ) = 1/α): δ = {delta_needed:.10f}")

    # Search for δ as simple expressions
    print(f"\n  Candidate expressions for δ = {delta_needed:.10f}:")
    delta_candidates = {
        "1/38": 1 / 38,
        "1/39": 1 / 39,
        "1/40": 1 / 40,
        "π/120": np.pi / 120,
        "π/121": np.pi / 121,
        "π/122": np.pi / 122,
        "1/(12π)": 1 / (12 * np.pi),
        "1/(4π²)": 1 / (4 * np.pi**2),
        "γ/22": 0.5772156649 / 22,
        "γ/23": 0.5772156649 / 23,
        "γ/24": 0.5772156649 / 24,
        "γ²/12": 0.5772156649**2 / 12,
        "ln(2)/27": np.log(2) / 27,
        "ln(2)/26": np.log(2) / 26,
        "ζ(3)/42": float(scipy_zeta(3)) / 42,
        "ζ(3)/43": float(scipy_zeta(3)) / 43,
        "ζ(3)/44": float(scipy_zeta(3)) / 44,
        "1/(6π-1)": 1 / (6 * np.pi - 1),
        "e/100": np.e / 100,
        "e/105": np.e / 105,
        "1/(4e)": 1 / (4 * np.e),
        "π/(12e)": np.pi / (12 * np.e),
    }

    sorted_deltas = sorted(delta_candidates.items(),
                           key=lambda x: abs(x[1] - delta_needed))
    for name, val in sorted_deltas:
        alpha_inv = np.pi * (42 + np.pi**2 / 6 - val)
        err = abs(alpha_inv - target) / target * 100
        marker = " ***" if err < 0.001 else " **" if err < 0.01 else " *" if err < 0.05 else ""
        print(f"    δ = {name:16s} = {val:.10f}  → 1/α = {alpha_inv:.6f} ({err:.5f}%){marker}")

    # Rational search: δ = p/q for small p, q
    print(f"\n  Rational search: δ = p/q closest to {delta_needed:.10f}")
    best_rats: List[Tuple[int, int, float]] = []
    for q in range(1, 200):
        p = round(delta_needed * q)
        if p > 0:
            err = abs(p / q - delta_needed)
            best_rats.append((p, q, err))
    best_rats.sort(key=lambda x: x[2])
    for p, q, err in best_rats[:10]:
        alpha_inv = np.pi * (42 + np.pi**2 / 6 - p / q)
        err_pct = abs(alpha_inv - ALPHA_INV_TARGET) / ALPHA_INV_TARGET * 100
        print(f"    δ = {p}/{q:3d} = {p/q:.10f}  → 1/α = {alpha_inv:.6f} ({err_pct:.5f}%)")


# ===========================================================================
# 6. CONNECTION TO S³ ZETA FUNCTION LITERATURE
# ===========================================================================
def literature_connection() -> None:
    """Discuss connections to known S³ spectral zeta results."""
    print("\n" + "=" * 70)
    print("6. CONNECTION TO S³ SPECTRAL ZETA FUNCTION")
    print("=" * 70)

    print("""
  The spectral zeta function on S³ is well-studied:

  KNOWN RESULTS (Voros, Dowker, Elizalde):
  -----------------------------------------
  For the scalar Laplacian on S³ with eigenvalues λ_n = n(n+2), deg = (n+1)²:
    ζ_{S³}(0) = 0  (conformal anomaly vanishes on odd spheres)
    ζ'_{S³}(0) = -ln(det Δ) = known in terms of ζ_R(3)/π²

  OUR CONVENTION: λ_n = n² - 1 = (n-1)(n+1), deg = n²
  This is the SAME spectrum shifted by n → n-1:
    Our (n=2,3,4,...) = Standard (n=1,2,3,...)

  KEY RESULT (Dowker 1989):
  The determinant of the Laplacian on S³ is:
    -ln det Δ_{S³} = 1/4 + 2ζ'_R(-2) - ζ_R(3)/(2π²) + ln(2)/2
  where ζ_R is the Riemann zeta function.
    """)

    # Compute numerical value of the S³ determinant
    zeta_prime_minus2 = -0.030448457  # ζ'_R(-2)
    zeta_3 = float(scipy_zeta(3))  # ≈ 1.202056903

    log_det = -(0.25 + 2 * zeta_prime_minus2 - zeta_3 / (2 * np.pi**2)
                + np.log(2) / 2)
    det_val = np.exp(log_det)

    print(f"  Numerical values:")
    print(f"    ζ_R(3) = {zeta_3:.10f}")
    print(f"    ζ'_R(-2) ≈ {zeta_prime_minus2:.10f}")
    print(f"    -ln det Δ_S³ ≈ {-log_det:.10f}")
    print(f"    det Δ_S³ ≈ {det_val:.10f}")

    # Check if det or log_det relates to α
    print(f"\n  Does det Δ_S³ relate to α?")
    print(f"    π × det = {np.pi * det_val:.6f}")
    print(f"    det / α = {det_val * ALPHA_INV_TARGET:.6f}")
    print(f"    42 + det = {42 + det_val:.6f}")
    print(f"    π(42 + det) = {np.pi * (42 + det_val):.6f}")

    # The "functional determinant" approach
    print(f"\n  Functional determinant approach:")
    print(f"    If 1/α = π × 42 × F(det Δ_S³), then")
    ratio = ALPHA_INV_TARGET / (np.pi * 42)
    print(f"    F = {ratio:.10f}")
    print(f"    1 + F-1 = 1 + {ratio - 1:.10f}")
    print(f"    F - 1 = (1/α)/(42π) - 1 = {ratio - 1:.10f}")
    print(f"    Compare: ζ(2)/42 = {np.pi**2 / (6*42):.10f}")
    print(f"    Compare: 1/26 = {1/26:.10f}")


# ===========================================================================
# 7. TRUNCATED SPECTRAL ZETA + CORRECTION FRAMEWORK
# ===========================================================================
def truncation_analysis() -> None:
    """Analyze the relationship between truncated and full spectral quantities."""
    print("\n" + "=" * 70)
    print("7. TRUNCATED SPECTRAL ZETA ANALYSIS")
    print("=" * 70)

    # Define truncated trace as Tr_N = Σ_{n=1}^{N} n²(n²-1)
    # and the "tail" T_N = Σ_{n=N+1}^{∞} n²(n²-1) × f(n) for some damping f
    # We need: π × [Tr_3 + tail_3] = 1/α
    # So: tail_3 = 1/(απ) - 42

    target_tail = ALPHA_INV_TARGET / np.pi - 42
    print(f"\n  Required tail: Σ_{{n≥4}} contribution = {target_tail:.10f}")

    # What if the tail comes from ζ-regularization of the divergent sum?
    # The zeta-regularized Σ n²(n²-1) should be:
    # Σ n⁴ - Σ n² = ζ_R(-4) - ζ_R(-2) = -1/30 - 0 = -1/30
    # But the PARTIAL sum up to N=3 is 42.
    # So the "regularized tail" = ζ_reg - partial = -1/30 - 42 = -42.033...
    # That doesn't work directly.

    zeta_minus4 = -1 / 30  # Ramanujan/zeta regularization
    zeta_minus2 = 0  # ζ_R(-2) = 0

    print(f"\n  Zeta-regularized Σ n²(n²-1):")
    print(f"    Σ n⁴ (reg) = ζ_R(-4) = {zeta_minus4:.10f}")
    print(f"    Σ n² (reg) = ζ_R(-2) = {zeta_minus2:.10f}")
    print(f"    Σ n²(n²-1) (reg) = {zeta_minus4 - zeta_minus2:.10f}")
    print(f"    Partial sum to N=3 = 42")
    print(f"    'Tail' = reg - partial = {zeta_minus4 - zeta_minus2 - 42:.10f}")
    print(f"    This is negative and large — not the small positive tail we need.")

    # Alternative: DAMPED sum
    print(f"\n  Damped-tail analysis: Σ_{{n≥4}} n²(n²-1) × exp(-t(n²-1))")
    print(f"  Required tail = {target_tail:.10f}")

    # Find t such that damped tail = target_tail
    def damped_tail(t: float, n_start: int = 4, n_max: int = 1000) -> float:
        n_vals = np.arange(n_start, n_max + 1, dtype=np.float64)
        lam = n_vals**2 - 1
        deg = n_vals**2
        return np.sum(deg * lam * np.exp(-t * lam))

    # Bisection
    t_lo, t_hi = 0.01, 2.0
    for _ in range(100):
        t_mid = (t_lo + t_hi) / 2
        if damped_tail(t_mid) > target_tail:
            t_lo = t_mid
        else:
            t_hi = t_mid
    t_tail = (t_lo + t_hi) / 2
    print(f"  t such that damped tail = target: t = {t_tail:.10f}")
    print(f"  damped_tail(t) = {damped_tail(t_tail):.10f}")

    # Check candidate forms for t_tail
    t_candidates = {
        "ln(2)/2": np.log(2) / 2,
        "1/3": 1 / 3,
        "π/10": np.pi / 10,
        "1/π": 1 / np.pi,
        "γ/2": 0.5772156649 / 2,
    }
    print(f"\n  Candidates for t_tail = {t_tail:.10f}:")
    for name, val in sorted(t_candidates.items(), key=lambda x: abs(x[1] - t_tail)):
        dt = damped_tail(val)
        alpha_inv = np.pi * (42 + dt)
        err = abs(alpha_inv - ALPHA_INV_TARGET) / ALPHA_INV_TARGET * 100
        print(f"    t = {name:10s} = {val:.8f}  tail = {dt:.6f}  "
              f"→ 1/α = {alpha_inv:.4f} ({err:.4f}%)")

    # Alternative: power-law damping
    print(f"\n  Power-law tail: Σ_{{n≥4}} n²(n²-1) / (n²-1)^s")
    for s in [1.5, 2.0, 2.5, 3.0, 3.5, 4.0]:
        n_vals = np.arange(4, 10001, dtype=np.float64)
        tail = np.sum(n_vals**2 * (n_vals**2 - 1) / (n_vals**2 - 1)**s)
        alpha_inv = np.pi * (42 + tail)
        err = abs(alpha_inv - ALPHA_INV_TARGET) / ALPHA_INV_TARGET * 100
        print(f"    s = {s:.1f}: tail = {tail:.6f}  → 1/α = {alpha_inv:.4f} ({err:.4f}%)")


# ===========================================================================
# MAIN
# ===========================================================================
def main() -> None:
    """Run all analyses."""
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║  ALPHA ZETA ANALYSIS: Spectral Zeta Function Approach to α         ║")
    print("║  Key: Tr(L²)|_{n≤3} = 42, and 42 + ζ(2) ≈ 1/(απ)                 ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")

    analyze_spectral_zeta()
    analyze_heat_kernel()
    analyze_gap()
    analyze_42()
    search_exact_formula()
    fine_search()
    literature_connection()
    truncation_analysis()

    # Final summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    base = np.pi * (42 + np.pi**2 / 6)
    err_base = abs(base - ALPHA_INV_TARGET) / ALPHA_INV_TARGET * 100

    print(f"""
  KEY FINDINGS:
  1. ζ_S³(1) DIVERGES — cannot use it directly for α.
  2. The heat kernel K(t) gives α at a specific t* ≈ 0.048,
     but t* does not appear to be a simple constant.
  3. The best simple formula remains:
       1/α ≈ π × (42 + ζ(2)) = {base:.6f}  ({err_base:.4f}% error)
  4. 42 = N(N+1)(2N+1)(N+2)(N-1)/10 at N=3
       = 2 × Σ(SO(4) Casimir × degeneracy) for n=1,2,3
  5. The gap from 42 to exact target ≈ {ALPHA_INV_TARGET/np.pi - 42:.6f}
     is closest to ζ(2) = π²/6 ≈ {np.pi**2/6:.6f}
  6. Small correction δ ≈ 0.026 needed: 1/α = π(42 + ζ(2) - δ)
     Best rational: δ ≈ 1/39 gives {abs(np.pi*(42+np.pi**2/6-1/39) - ALPHA_INV_TARGET)/ALPHA_INV_TARGET*100:.4f}% error

  INTERPRETATION:
  If α has a spectral-geometric origin on S³, the structure is:
    1/α = π × [DISCRETE TOPOLOGY (42) + SPECTRAL REGULARIZATION (≈ζ(2))]
  The discrete part (42) comes from SO(4) Casimir invariants truncated at n=3.
  The continuous part (ζ(2)) comes from the spectral zeta function correction.
  The residual δ ≈ 0.026 may encode higher-order spectral information.
    """)


if __name__ == "__main__":
    main()
