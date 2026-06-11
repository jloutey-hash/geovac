"""
Alpha eigenvalue search: discrete lattice origin of the fine structure constant.

Hypothesis: 1/α = π × X, where X ≈ 43.62 arises from angular momentum
eigenvalue sums on the discrete S³ lattice. This script systematically
searches for the discrete structure (trace, ratio, combination) that
produces X = 4π² + π + 1.

Date: 2026-03-21
Status: Exploratory (Paper 2 territory)
"""

from __future__ import annotations

import numpy as np
from typing import List, Tuple, Dict
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

# Physical constants
ALPHA_INV = 137.035999084  # 1/α (CODATA 2018)
TARGET_X = ALPHA_INV / np.pi  # 43.6196...
TARGET_FORMULA = 4 * np.pi**2 + np.pi + 1  # 43.6196...

print("=" * 72)
print("ALPHA EIGENVALUE SEARCH")
print("=" * 72)
print(f"  1/α           = {ALPHA_INV:.10f}")
print(f"  1/(πα)        = {TARGET_X:.10f}")
print(f"  4π² + π + 1   = {TARGET_FORMULA:.10f}")
print(f"  Match?          {abs(TARGET_X - TARGET_FORMULA):.2e}")
print()


# ===================================================================
# SECTION 1: L² TRACE — Σ l(l+1)(2l+1) over shells
# ===================================================================

def l2_trace(n_max: int) -> int:
    """Compute Σ_{n=1}^{n_max} Σ_{l=0}^{n-1} l(l+1)(2l+1)."""
    total = 0
    for n in range(1, n_max + 1):
        for l in range(n):
            total += l * (l + 1) * (2 * l + 1)
    return total


def l2_trace_single_shell(n: int) -> int:
    """Σ_{l=0}^{n-1} l(l+1)(2l+1) for a single shell n."""
    return sum(l * (l + 1) * (2 * l + 1) for l in range(n))


print("=" * 72)
print("SECTION 1: L² TRACE = Σ l(l+1)(2l+1)")
print("=" * 72)
print()
print(f"{'n_max':>5} {'Tr(L²)':>10} {'Tr(L²)/target':>14} {'Δ from target':>14}")
print("-" * 50)

for n_max in range(1, 11):
    tr = l2_trace(n_max)
    ratio = tr / TARGET_X
    delta = tr - TARGET_X
    marker = " <<<" if abs(delta) < 5 else ""
    print(f"{n_max:5d} {tr:10d} {ratio:14.6f} {delta:14.4f}{marker}")

print()
print("Per-shell contributions:")
print(f"{'n':>5} {'Shell Tr(L²)':>12} {'Cumulative':>12}")
print("-" * 35)
cumulative = 0
for n in range(1, 8):
    shell = l2_trace_single_shell(n)
    cumulative += shell
    print(f"{n:5d} {shell:12d} {cumulative:12d}")

# Closed-form derivation
# Σ_{l=0}^{N-1} l(l+1)(2l+1) = Σ (2l³ + 3l² + l)
# = 2·[(N-1)N/2]² + 3·(N-1)N(2N-1)/6 + (N-1)N/2
# = (N-1)²N²/2 + (N-1)N(2N-1)/2 + (N-1)N/2
# = (N-1)N/2 · [(N-1)N + (2N-1) + 1]
# = (N-1)N/2 · [N² - N + 2N - 1 + 1]
# = (N-1)N/2 · [N² + N]
# = (N-1)N²(N+1)/2

print()
print("CLOSED-FORM VERIFICATION:")
print("  Σ_{l=0}^{N-1} l(l+1)(2l+1) = N(N-1)·N·(N+1)/2 = N²(N²-1)/2")
print()
for n in range(1, 8):
    formula = n * n * (n * n - 1) // 2
    direct = l2_trace_single_shell(n)
    print(f"  n={n}: formula={formula}, direct={direct}, match={formula == direct}")

# Cumulative: Σ_{n=1}^{nmax} n²(n²-1)/2
# = (1/2) Σ (n⁴ - n²)
# = (1/2) [nmax(nmax+1)(2nmax+1)(3nmax²+3nmax-1)/30 - nmax(nmax+1)(2nmax+1)/6]
print()
print("CUMULATIVE CLOSED FORM:")
print("  Σ_{n=1}^{N} n²(n²-1)/2 = (1/2)[Σn⁴ - Σn²]")
for n_max in range(1, 11):
    # Σ n⁴ = N(N+1)(2N+1)(3N²+3N-1)/30
    N = n_max
    sum_n4 = N * (N + 1) * (2 * N + 1) * (3 * N**2 + 3 * N - 1) // 30
    sum_n2 = N * (N + 1) * (2 * N + 1) // 6
    formula = (sum_n4 - sum_n2) // 2
    direct = l2_trace(n_max)
    assert formula == direct, f"Mismatch at n_max={n_max}"
print("  All verified ✓")

# At what continuous n does cumulative = TARGET_X?
print()
print("CONTINUOUS SOLUTION:")
print(f"  Target: Σ = {TARGET_X:.6f}")
print(f"  n=3: Σ = {l2_trace(3)} (undershoot by {TARGET_X - l2_trace(3):.4f})")
print(f"  n=4: Σ = {l2_trace(4)} (overshoot by {l2_trace(4) - TARGET_X:.4f})")

# Solve (1/2)[Σn⁴ - Σn²] ≈ TARGET_X numerically using continuous approximation
# Σn⁴ ≈ N⁵/5, Σn² ≈ N³/3 for large N
# (1/2)(N⁵/5 - N³/3) = TARGET_X → N⁵/10 ≈ TARGET_X → N ≈ (10·TARGET_X)^(1/5)
from scipy.optimize import brentq

def cumulative_continuous(x: float) -> float:
    """Continuous extension of the cumulative L² trace."""
    # Use the exact polynomial form with continuous N
    N = x
    sum_n4 = N * (N + 1) * (2 * N + 1) * (3 * N**2 + 3 * N - 1) / 30
    sum_n2 = N * (N + 1) * (2 * N + 1) / 6
    return (sum_n4 - sum_n2) / 2

n_star = brentq(lambda x: cumulative_continuous(x) - TARGET_X, 3, 4)
print(f"  Continuous n* = {n_star:.8f}")
print(f"  n* - 3 = {n_star - 3:.8f}")
print(f"  n* - π = {n_star - np.pi:.8f}")
print(f"  π ≈ {np.pi:.8f}")
print(f"  n* / π = {n_star / np.pi:.8f}")


# ===================================================================
# SECTION 2: OTHER CASIMIR TRACES
# ===================================================================

print()
print("=" * 72)
print("SECTION 2: OTHER CASIMIR TRACES")
print("=" * 72)
print()

def so4_casimir_trace(n_max: int) -> int:
    """SO(4) Casimir: Σ n²(n²-1) with degeneracy n²."""
    return sum(n**2 * (n**2 - 1) for n in range(1, n_max + 1))


def mixed_trace(n_max: int) -> int:
    """Mixed: Σ n² × l(l+1) with degeneracy (2l+1)."""
    total = 0
    for n in range(1, n_max + 1):
        for l in range(n):
            total += n**2 * l * (l + 1) * (2 * l + 1)
    return total


def dim_trace(n_max: int) -> int:
    """Dimension: Σ n² (total number of states)."""
    return sum(n**2 for n in range(1, n_max + 1))


def l_trace(n_max: int) -> int:
    """Σ l(2l+1) over all states."""
    total = 0
    for n in range(1, n_max + 1):
        for l in range(n):
            total += l * (2 * l + 1)
    return total


# Collect all traces
traces: Dict[str, List[int]] = {}
trace_names = ["Tr(1)=Σn²", "Tr(L²)=Σl(l+1)(2l+1)", "Tr(C₄)=Σn²(n²-1)",
               "Tr(n²L²)", "Tr(L)=Σl(2l+1)"]
trace_fns = [dim_trace, l2_trace, so4_casimir_trace, mixed_trace, l_trace]

for name, fn in zip(trace_names, trace_fns):
    traces[name] = [fn(n) for n in range(1, 11)]

print(f"{'n_max':>5}", end="")
for name in trace_names:
    print(f" {name:>18}", end="")
print()
print("-" * (5 + 19 * len(trace_names)))

for i, n_max in enumerate(range(1, 11)):
    print(f"{n_max:5d}", end="")
    for name in trace_names:
        print(f" {traces[name][i]:18d}", end="")
    print()

# Search ratios
print()
print("RATIO SEARCH (looking for ≈ 43.62 or ≈ 137.04):")
print("-" * 72)

hits: List[Tuple[str, int, float, float]] = []

for n_max in range(1, 11):
    vals = {name: traces[name][n_max - 1] for name in trace_names}

    # All pairwise ratios
    for name_a in trace_names:
        for name_b in trace_names:
            if name_a == name_b or vals[name_b] == 0:
                continue
            ratio = vals[name_a] / vals[name_b]
            for target, target_name in [(TARGET_X, "43.62"), (ALPHA_INV, "137.04")]:
                if abs(ratio - target) / target < 0.05:  # within 5%
                    hits.append((f"{name_a}/{name_b}", n_max, ratio, target))

    # Also check: Tr(C₄) / Tr(1)  (= average n²-1)
    avg_casimir = vals["Tr(C₄)=Σn²(n²-1)"] / vals["Tr(1)=Σn²"] if vals["Tr(1)=Σn²"] > 0 else 0
    for target, target_name in [(TARGET_X, "43.62"), (ALPHA_INV, "137.04")]:
        if abs(avg_casimir - target) / target < 0.05:
            hits.append(("⟨C₄⟩=Tr(C₄)/Tr(1)", n_max, avg_casimir, target))

if hits:
    for desc, n_max, val, target in hits:
        err = (val - target) / target * 100
        print(f"  n_max={n_max}: {desc} = {val:.6f} (target {target:.2f}, err {err:+.2f}%)")
else:
    print("  No simple ratio hits within 5%.")


# ===================================================================
# SECTION 3: TRANSITION OPERATOR EIGENVALUES
# ===================================================================

print()
print("=" * 72)
print("SECTION 3: TRANSITION OPERATOR EIGENVALUES")
print("=" * 72)
print()

def build_dipole_operator(n_max: int) -> np.ndarray:
    """
    Build the dipole transition operator T on the lattice.

    States: |n, l, m⟩ for n=1..n_max, l=0..n-1, m=-l..l
    Selection rule: Δl = ±1, Δm = 0, ±1 (electric dipole)

    Matrix elements: T_{ij} = 1 for allowed transitions (unweighted).
    """
    # Build state list
    states: List[Tuple[int, int, int]] = []
    for n in range(1, n_max + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                states.append((n, l, m))

    N = len(states)
    T = np.zeros((N, N))

    for i, (n1, l1, m1) in enumerate(states):
        for j, (n2, l2, m2) in enumerate(states):
            # Dipole selection rules
            if abs(l1 - l2) == 1 and abs(m1 - m2) <= 1:
                T[i, j] = 1.0

    return T


def build_weighted_dipole(n_max: int) -> np.ndarray:
    """
    Dipole operator weighted by CG coefficients.

    For Δl = +1: weight = sqrt((l+1)²-m²) / sqrt((2l+1)(2l+3))
    For Δl = -1: weight = sqrt(l²-m²) / sqrt((2l-1)(2l+1))
    Only Δm = 0 component (z-polarization).
    """
    states: List[Tuple[int, int, int]] = []
    for n in range(1, n_max + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                states.append((n, l, m))

    N = len(states)
    T = np.zeros((N, N))

    for i, (n1, l1, m1) in enumerate(states):
        for j, (n2, l2, m2) in enumerate(states):
            if m1 != m2:
                continue  # Δm = 0 only for z-component
            if l2 == l1 + 1:
                # l → l+1: C^{l+1,m}_{l,m;1,0}
                l = l1
                if (2 * l + 1) * (2 * l + 3) > 0:
                    T[i, j] = np.sqrt(((l + 1)**2 - m1**2) /
                                       ((2 * l + 1) * (2 * l + 3)))
            elif l2 == l1 - 1:
                # l → l-1: C^{l-1,m}_{l,m;1,0}
                l = l1
                if (2 * l - 1) * (2 * l + 1) > 0:
                    T[i, j] = np.sqrt((l**2 - m1**2) /
                                       ((2 * l - 1) * (2 * l + 1)))

    return T


for n_max in range(2, 7):
    T_simple = build_dipole_operator(n_max)
    T_weighted = build_weighted_dipole(n_max)

    eigs_simple = np.sort(np.linalg.eigvalsh(T_simple))[::-1]
    eigs_weighted = np.sort(np.linalg.eigvalsh(T_weighted))[::-1]

    # Trace of T²
    tr_T2_simple = np.trace(T_simple @ T_simple)
    tr_T2_weighted = np.trace(T_weighted @ T_weighted)

    n_states = T_simple.shape[0]
    print(f"n_max={n_max} ({n_states} states):")
    print(f"  T (unweighted): λ_max={eigs_simple[0]:.6f}, "
          f"Tr(T²)={tr_T2_simple:.1f}, "
          f"Σ|λ|={np.sum(np.abs(eigs_simple)):.4f}")
    print(f"  T (CG-weighted): λ_max={eigs_weighted[0]:.6f}, "
          f"Tr(T²)={tr_T2_weighted:.6f}, "
          f"Σ|λ|={np.sum(np.abs(eigs_weighted)):.4f}")

    # Check various combinations against targets
    for target, tname in [(TARGET_X, "43.62"), (ALPHA_INV, "137.04")]:
        for val, vname in [
            (eigs_simple[0], "λ_max(simple)"),
            (tr_T2_simple, "Tr(T²_simple)"),
            (tr_T2_weighted, "Tr(T²_weighted)"),
            (np.sum(np.abs(eigs_simple)), "Σ|λ|(simple)"),
            (np.sum(np.abs(eigs_weighted)), "Σ|λ|(weighted)"),
        ]:
            if val > 0 and abs(val - target) / target < 0.1:
                print(f"  *** HIT: {vname} = {val:.4f} ≈ {tname} "
                      f"(err {(val-target)/target*100:+.2f}%)")
    print()


# ===================================================================
# SECTION 4: RATIOS AND PRODUCTS
# ===================================================================

print("=" * 72)
print("SECTION 4: SYSTEMATIC RATIO/PRODUCT SEARCH")
print("=" * 72)
print()

# Compute many lattice quantities and search all pairwise combinations
def compute_lattice_quantities(n_max: int) -> Dict[str, float]:
    """Compute various lattice quantities for given n_max."""
    q: Dict[str, float] = {}

    # Basic traces
    q["Tr(1)"] = sum(n**2 for n in range(1, n_max + 1))
    q["Tr(L²)"] = l2_trace(n_max)
    q["Tr(C₄)"] = so4_casimir_trace(n_max)
    q["Tr(n²L²)"] = mixed_trace(n_max)
    q["Tr(L)"] = l_trace(n_max)

    # Averages
    if q["Tr(1)"] > 0:
        q["⟨L²⟩"] = q["Tr(L²)"] / q["Tr(1)"]
        q["⟨C₄⟩"] = q["Tr(C₄)"] / q["Tr(1)"]
        q["⟨n²L²⟩"] = q["Tr(n²L²)"] / q["Tr(1)"]

    # Additional traces
    q["Σn"] = sum(n for n in range(1, n_max + 1))
    q["Σn³"] = sum(n**3 for n in range(1, n_max + 1))
    q["Σn⁴"] = sum(n**4 for n in range(1, n_max + 1))
    q["Σn²(n²-1)²"] = sum(n**2 * (n**2 - 1)**2 for n in range(1, n_max + 1))

    # Shell-specific
    q["n_max²"] = n_max**2
    q["n_max²-1"] = n_max**2 - 1
    q["n_max(n_max+1)/2"] = n_max * (n_max + 1) / 2

    # SO(4) dimension = n², total states up to n_max = N(N+1)(2N+1)/6
    q["dim_total"] = q["Tr(1)"]

    # Ratio of consecutive
    q["Tr(C₄)/Tr(L²)"] = q["Tr(C₄)"] / q["Tr(L²)"] if q["Tr(L²)"] > 0 else 0
    q["Tr(n²L²)/Tr(L²)"] = q["Tr(n²L²)"] / q["Tr(L²)"] if q["Tr(L²)"] > 0 else 0
    q["Tr(n²L²)/Tr(C₄)"] = q["Tr(n²L²)"] / q["Tr(C₄)"] if q["Tr(C₄)"] > 0 else 0

    return q


print("Searching for values ≈ 43.62 or ≈ 137.04 among all quantities and ratios...")
print()

all_hits: List[Tuple[int, str, float, str, float]] = []

for n_max in range(1, 11):
    q = compute_lattice_quantities(n_max)

    # Check each quantity directly
    for name, val in q.items():
        if val == 0:
            continue
        for target, tname in [(TARGET_X, "43.62"), (ALPHA_INV, "137.04")]:
            if abs(val - target) / target < 0.02:  # 2% tolerance
                all_hits.append((n_max, name, val, tname, (val - target) / target * 100))

    # Check all pairwise ratios
    names = list(q.keys())
    for i in range(len(names)):
        for j in range(len(names)):
            if i == j:
                continue
            vi, vj = q[names[i]], q[names[j]]
            if vj == 0 or vi == 0:
                continue
            ratio = vi / vj
            for target, tname in [(TARGET_X, "43.62"), (ALPHA_INV, "137.04")]:
                if abs(ratio - target) / target < 0.02:
                    rname = f"{names[i]}/{names[j]}"
                    all_hits.append((n_max, rname, ratio, tname, (ratio - target) / target * 100))

    # Also check products of small integer multiples of traces
    for name_a in ["Tr(L²)", "Tr(C₄)", "Tr(1)", "Tr(L)"]:
        va = q[name_a]
        if va == 0:
            continue
        # Try: π × va, 4π² + va, va + π, etc.
        combos = {
            f"π×{name_a}": np.pi * va,
            f"4π²+{name_a}": 4 * np.pi**2 + va,
            f"{name_a}+π": va + np.pi,
            f"{name_a}+π²": va + np.pi**2,
            f"{name_a}+π/2": va + np.pi / 2,
            f"π×{name_a}+π²": np.pi * va + np.pi**2,
        }
        for cname, cval in combos.items():
            for target, tname in [(TARGET_X, "43.62"), (ALPHA_INV, "137.04")]:
                if abs(cval - target) / target < 0.02:
                    all_hits.append((n_max, cname, cval, tname,
                                     (cval - target) / target * 100))

if all_hits:
    # Sort by absolute error
    all_hits.sort(key=lambda x: abs(x[4]))
    print(f"Found {len(all_hits)} hits (within 2%):")
    print()
    print(f"{'n_max':>5} {'Quantity':>40} {'Value':>14} {'Target':>8} {'Err%':>8}")
    print("-" * 80)
    for n_max, name, val, tname, err in all_hits[:30]:
        print(f"{n_max:5d} {name:>40} {val:14.6f} {tname:>8} {err:+8.4f}%")
else:
    print("  No hits found within 2%.")


# ===================================================================
# SECTION 5: THE MAGIC CUTOFF — WHY n=3?
# ===================================================================

print()
print("=" * 72)
print("SECTION 5: WHY n_max = 3?")
print("=" * 72)
print()

# n=3 shell properties
print("n=3 shell properties:")
print(f"  SO(4) Casimir C₄(3) = n²-1 = {3**2 - 1}")
print(f"  States in shell: n² = {3**2}")
print(f"  Orbital types: s (l=0), p (l=1), d (l=2) — all shapes up to quadratic")
print(f"  Total states (1..3): {sum(n**2 for n in range(1,4))} = 1 + 4 + 9 = 14")
print(f"  Note: 14 = dim SO(7) / SO(5) × SO(2)")
print()

# The gap: Tr(L², n≤3) = 42, target = 43.62, gap = 1.62 ≈ π/2
gap = TARGET_X - l2_trace(3)
print(f"Gap analysis:")
print(f"  Tr(L², n≤3) = {l2_trace(3)}")
print(f"  Target       = {TARGET_X:.6f}")
print(f"  Gap          = {gap:.6f}")
print(f"  π/2          = {np.pi/2:.6f}")
print(f"  Gap - π/2    = {gap - np.pi/2:.6f}")
print(f"  Gap / π      = {gap / np.pi:.6f}")
print(f"  Gap = {gap:.6f}")
print()

# Is the gap related to a continuous correction?
# If the discrete sum is 42 and the continuous answer is 43.62,
# maybe there's a "fractional shell" contribution
print("Euler-Maclaurin correction approach:")
print("  Discrete sum misses continuous contributions.")
print("  The leading correction is f(N)/2 (trapezoidal correction).")
n3_shell_val = l2_trace_single_shell(4)  # next shell would contribute this
print(f"  Next shell (n=4) contributes: {n3_shell_val}")
print(f"  Half of next shell: {n3_shell_val / 2}")
print(f"  That's way too large ({n3_shell_val / 2} vs gap {gap:.4f})")
print()

# What fraction of the n=4 shell fills the gap?
frac = gap / n3_shell_val
print(f"  Fractional n=4 shell needed: {frac:.6f}")
print(f"  = {gap}/{n3_shell_val}")
print()

# Alternative: maybe it's n²(n²-1)/2 that needs a π² correction
# Tr(L², n≤3) = 42 and 42 + π/2 ≈ 43.57 ≈ TARGET_X - 0.05
print("Additive corrections to reach target:")
corrections = {
    "π/2": np.pi / 2,
    "π²/6 (= ζ(2))": np.pi**2 / 6,
    "π/2 + π²/90": np.pi / 2 + np.pi**2 / 90,
    "(π²+3)/6": (np.pi**2 + 3) / 6,
    "ln(4π²)": np.log(4 * np.pi**2),
    "e/√(2π)": np.e / np.sqrt(2 * np.pi),
    "√(2π) - √3": np.sqrt(2 * np.pi) - np.sqrt(3),
}
print(f"  42 + correction = ?  (target = {TARGET_X:.6f})")
for name, val in corrections.items():
    result = 42 + val
    err = result - TARGET_X
    print(f"    42 + {name:20s} = {result:.6f}  (Δ = {err:+.6f})")


# ===================================================================
# SECTION 6: CLOSED-FORM SEARCH
# ===================================================================

print()
print("=" * 72)
print("SECTION 6: CLOSED-FORM ANALYSIS")
print("=" * 72)
print()

# The key identity: 4π² + π + 1
# Can this be written as a function of integers and π?
# 4π² + π + 1 = (2π)² + π + 1 = π(4π + 1) + 1

print("Decompositions of 4π² + π + 1:")
print(f"  = (2π)² + π + 1 = {(2*np.pi)**2 + np.pi + 1:.8f}")
print(f"  = π(4π + 1) + 1 = {np.pi*(4*np.pi+1) + 1:.8f}")
print(f"  = 4(π² + 1/4) + π = 4×{np.pi**2 + 0.25:.6f} + π")
print(f"  = (2π + 1/2)² + 3/4 = {(2*np.pi + 0.5)**2 + 0.75:.8f}")
actual = (2 * np.pi + 0.5)**2 + 0.75
print(f"    vs target: {TARGET_X:.8f}, diff = {actual - TARGET_X:.8f}")
print()

# Check: is 4π² + π + 1 related to a trace on a CONTINUOUS S³?
# On unit S³, eigenvalues of Laplacian are -n(n+2) with degeneracy (n+1)²
# (here n starts from 0). This is different from our convention.
print("S³ Laplacian eigenvalues (continuous):")
print("  λ_n = -n(n+2), degeneracy (n+1)², n = 0, 1, 2, ...")
print("  In hydrogen convention: λ = -(n²-1), degeneracy n², n = 1, 2, ...")
print()

# Continuous S³ traces with these eigenvalues
print("Spectral zeta function approach:")
print("  ζ_S³(s) = Σ_{n=1}^∞ n² / (n²-1)^s  (regularized)")
print("  This diverges for Re(s) ≤ 3/2, but regularized values may give α.")
print()

# Heat kernel: K(t) = Σ n² exp(-t(n²-1))
print("Heat kernel on S³:")
print("  K(t) = Σ_{n=1}^∞ n² exp(-t(n²-1))")
for t_val in [0.01, 0.05, 0.1, 0.5, 1.0]:
    K = sum(n**2 * np.exp(-t_val * (n**2 - 1)) for n in range(1, 100))
    print(f"  K({t_val:.2f}) = {K:.6f}", end="")
    if abs(K - TARGET_X) / TARGET_X < 0.05:
        print(f"  *** CLOSE TO {TARGET_X:.2f}!", end="")
    if abs(K - ALPHA_INV) / ALPHA_INV < 0.05:
        print(f"  *** CLOSE TO {ALPHA_INV:.2f}!", end="")
    print()

# Find t where K(t) = TARGET_X
print()
print("Solving K(t*) = 43.62:")

def heat_kernel(t: float, n_terms: int = 200) -> float:
    return sum(n**2 * np.exp(-t * (n**2 - 1)) for n in range(1, n_terms + 1))

# K(0) = Σn² → ∞, K(∞) → 1 (n=1 term), so K(t) passes through 43.62
try:
    t_star = brentq(lambda t: heat_kernel(t) - TARGET_X, 0.001, 2.0)
    print(f"  t* = {t_star:.10f}")
    print(f"  1/t* = {1/t_star:.10f}")
    print(f"  π·t* = {np.pi * t_star:.10f}")
    print(f"  t*/π = {t_star / np.pi:.10f}")
    print(f"  √t* = {np.sqrt(t_star):.10f}")
    print(f"  K(t*) = {heat_kernel(t_star):.10f}")
except Exception as e:
    print(f"  Failed: {e}")

# Also: K(t) = TARGET_X × π → 1/α
print()
print("Solving K(t*) = 137.04:")
try:
    t_star2 = brentq(lambda t: heat_kernel(t) - ALPHA_INV, 0.001, 2.0)
    print(f"  t* = {t_star2:.10f}")
    print(f"  1/t* = {1/t_star2:.10f}")
    print(f"  π·t* = {np.pi * t_star2:.10f}")
except Exception as e:
    print(f"  Failed: {e}")


# ===================================================================
# SECTION 7: WEIGHTED L² TRACES (π-corrected)
# ===================================================================

print()
print("=" * 72)
print("SECTION 7: WEIGHTED / π-CORRECTED TRACES")
print("=" * 72)
print()

# The structure 4π² + π + 1 suggests contributions from different l:
# l=0 sector → contributes 1  (the "+1")
# l=1 sector → contributes π  (the "+π")
# l=2 sector → contributes 4π² (the "+4π²")
# Pattern: l sector contributes (some coeff) × π^l ?

print("Hypothesis: each l sector contributes c_l × π^l")
print()
for n_max in range(1, 7):
    print(f"n_max = {n_max}:")
    for l in range(n_max):
        # How many states with this l? Count (n,l,m) with n > l, m=-l..l
        # Number of shells containing l: n = l+1, l+2, ..., n_max
        n_shells = n_max - l
        degeneracy = (2 * l + 1) * n_shells  # states with this l
        l2_contrib = l * (l + 1) * (2 * l + 1) * n_shells  # l(l+1)(2l+1) × n_shells
        print(f"  l={l}: {n_shells} shells, {degeneracy} states, "
              f"Tr(L²) contrib = {l2_contrib}")

    # Now try: Σ_l [l(l+1)(2l+1) × n_shells_containing_l × π^l]
    weighted_sum = 0
    for l in range(n_max):
        n_shells = n_max - l
        weighted_sum += l * (l + 1) * (2 * l + 1) * n_shells * np.pi**l

    print(f"  Σ l(l+1)(2l+1) × n_shells × π^l = {weighted_sum:.6f}")
    if abs(weighted_sum - TARGET_X) / TARGET_X < 0.1:
        print(f"  *** Within 10% of target {TARGET_X:.2f}!")
    print()

# Alternative: weight by π^l / l! (exponential-like)
print("Alternative: weight by π^l / l! (Poisson-like):")
for n_max in range(1, 7):
    weighted_sum = 0
    for l in range(n_max):
        n_shells = n_max - l
        import math
        weight = np.pi**l / math.factorial(l) if l < 20 else 0
        weighted_sum += l * (l + 1) * (2 * l + 1) * n_shells * weight
    print(f"  n_max={n_max}: {weighted_sum:.6f}", end="")
    if abs(weighted_sum - TARGET_X) / TARGET_X < 0.05:
        print(f"  *** CLOSE ({(weighted_sum/TARGET_X - 1)*100:+.2f}%)", end="")
    if abs(weighted_sum - ALPHA_INV) / ALPHA_INV < 0.05:
        print(f"  *** CLOSE to 1/α ({(weighted_sum/ALPHA_INV - 1)*100:+.2f}%)", end="")
    print()

# Weight by 1/(2l+1) (removes degeneracy)
print()
print("Σ l(l+1) × n_shells (no degeneracy factor):")
for n_max in range(1, 7):
    s = sum(l * (l + 1) * (n_max - l) for l in range(n_max))
    print(f"  n_max={n_max}: {s}")


# ===================================================================
# SECTION 8: GRAPH LAPLACIAN SPECTRAL ANALYSIS
# ===================================================================

print()
print("=" * 72)
print("SECTION 8: ACTUAL GRAPH LAPLACIAN EIGENVALUE SUMS")
print("=" * 72)
print()

try:
    from geovac.lattice import GeometricLattice

    for n_max in [3, 4, 5, 6]:
        lat = GeometricLattice(max_n=n_max)
        # Get adjacency and degree
        A = lat.adjacency
        if hasattr(A, 'toarray'):
            A = A.toarray()
        D = np.diag(np.asarray(A.sum(axis=1)).flatten())
        L = D - A  # Graph Laplacian

        eigs = np.sort(np.linalg.eigvalsh(L))

        print(f"n_max={n_max} ({L.shape[0]} nodes):")
        print(f"  Graph Laplacian eigenvalues (first 10): "
              f"{', '.join(f'{e:.4f}' for e in eigs[:10])}")
        print(f"  Σλ (trace) = {np.sum(eigs):.6f}")
        print(f"  Σλ² = {np.sum(eigs**2):.6f}")
        print(f"  Σ1/λ (nonzero) = "
              f"{np.sum(1/eigs[eigs > 1e-10]):.6f}")
        print(f"  λ_max = {eigs[-1]:.6f}")
        print(f"  λ_max / π = {eigs[-1] / np.pi:.6f}")

        # Check ratios
        nonzero = eigs[eigs > 1e-10]
        for target, tname in [(TARGET_X, "43.62"), (ALPHA_INV, "137.04")]:
            for val, vname in [
                (np.sum(eigs), "Σλ"),
                (np.sum(eigs**2), "Σλ²"),
                (np.sum(1 / nonzero), "Σ1/λ"),
                (np.sum(eigs) / len(eigs), "⟨λ⟩"),
            ]:
                if val > 0 and abs(val - target) / target < 0.05:
                    print(f"  *** HIT: {vname} = {val:.6f} ≈ {tname}")
        print()

except ImportError:
    print("  Could not import geovac.lattice — skipping graph Laplacian analysis.")
except Exception as e:
    print(f"  Error in graph Laplacian analysis: {e}")


# ===================================================================
# SUMMARY
# ===================================================================

print()
print("=" * 72)
print("SUMMARY")
print("=" * 72)
print()
print(f"Target: 1/α = π × (4π² + π + 1) = {ALPHA_INV:.6f}")
print(f"                    4π² + π + 1   = {TARGET_X:.6f}")
print(f"                    Tr(L², n≤3)   = {l2_trace(3)}")
print(f"                    Gap            = {TARGET_X - l2_trace(3):.6f}")
print(f"                    π/2            = {np.pi/2:.6f}")
print()
print("Key findings:")
print(f"  1. L² trace closed form: Σ_{{n=1}}^N n²(n²-1)/2")
print(f"  2. Cumulative L² trace at n=3: 42 (undershoot by ~1.62)")
print(f"  3. Continuous n* where trace = target: {n_star:.6f}")
print(f"  4. Gap ≈ π/2 = {np.pi/2:.6f} (off by {gap - np.pi/2:.6f})")
print()
print("If the gap is exactly explained, we would have:")
print(f"  1/α = π × [Tr(L², n≤3) + correction]")
print(f"  where correction = 4π² + π + 1 - 42 = {TARGET_X - 42:.10f}")
print(f"                   = 4π² + π - 41     = {4*np.pi**2 + np.pi - 41:.10f}")
print()

# Final check: is 4π² + π - 41 recognizable?
val = 4 * np.pi**2 + np.pi - 41
print(f"Analysis of residual = 4π² + π - 41 = {val:.10f}:")
print(f"  π/2         = {np.pi/2:.10f} (Δ = {val - np.pi/2:.10f})")
print(f"  (π²-6)/3    = {(np.pi**2-6)/3:.10f} (Δ = {val - (np.pi**2-6)/3:.10f})")
print(f"  4(π²-10)+π  = {4*(np.pi**2-10)+np.pi:.10f} (this IS the residual)")
print(f"  π²-10       = {np.pi**2-10:.10f}")
print(f"  ≈ -0.1304   = close to -π²/72 = {-np.pi**2/72:.10f}")
print(f"                 or -1/(2π+1)   = {-1/(2*np.pi+1):.10f}")
print()

# So: 4π² + π + 1 = 42 + 4(π²-10) + π + 1 - 42 + 40 = ...
# Let's just see: 4π² = 4 × 9.8696 = 39.478
# So 4π² + π + 1 = 39.478 + 3.1416 + 1 = 43.620
# And 42 = Σ for n≤3 = 0 + 0 + 6 + 0 + 6 + 30 = wait let me recheck
print("Detailed shell breakdown for n≤3:")
for n in range(1, 4):
    for l in range(n):
        contrib = l * (l + 1) * (2 * l + 1)
        if contrib > 0:
            print(f"  n={n}, l={l}: l(l+1)(2l+1) = {l}×{l+1}×{2*l+1} = {contrib}")
print(f"  Total = {l2_trace(3)}")
print()
print("Done.")
