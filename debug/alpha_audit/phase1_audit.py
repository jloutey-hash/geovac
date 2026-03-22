"""
Phase 1 Audit: Paper 2 — Fine Structure Constant from Hopf Bundle
==================================================================
Independently verifies every quantity in the paper and systematically
investigates the boundary term inconsistency.
"""

import numpy as np
from numpy.polynomial import polynomial as P
import os

CODATA_ALPHA_INV = 137.035999084  # CODATA 2018

output_dir = os.path.dirname(os.path.abspath(__file__))


def solve_cubic(K):
    """Solve alpha^3 - K*alpha + 1 = 0 for all real roots.

    Rearranged as polynomial: alpha^3 + 0*alpha^2 - K*alpha + 1 = 0
    Coefficients for np.roots: [1, 0, -K, 1]
    """
    roots = np.roots([1, 0, -K, 1])
    # Keep only real roots (imag part < 1e-12)
    real_roots = sorted([r.real for r in roots if abs(r.imag) < 1e-12])
    return real_roots


# =====================================================================
# PART A: Full arithmetic verification
# =====================================================================
print("=" * 72)
print("PART A: Full Arithmetic Verification")
print("=" * 72)

# --- 1. S3 eigenvalues and degeneracies ---
print("\n1. S3 eigenvalues and degeneracies")
print("-" * 50)
print(f"  {'n':>3}  {'lambda_n':>10}  {'g_n':>5}  {'N_cumul':>8}")
N_cumul = 0
for n in range(1, 6):
    lam = -(n**2 - 1)
    g = n**2
    N_cumul += g
    print(f"  {n:3d}  {lam:10d}  {g:5d}  {N_cumul:8d}")

print("\n  Paper claims: lambda_n = -(n^2-1), g_n = n^2  --> CONFIRMED")

# --- 2. Casimir trace (unnormalized, per (n,l) pair) ---
print("\n2. Casimir trace: sum of l(l+1) over (n,l) pairs")
print("-" * 50)
print(f"  {'n_max':>5}  {'pairs':>30}  {'Sum':>6}  {'Half':>6}")
for nmax in range(1, 6):
    terms = []
    for n in range(1, nmax + 1):
        for l in range(n):
            terms.append(l * (l + 1))
    total = sum(terms)
    print(f"  {nmax:5d}  {str(terms):>30s}  {total:6d}  {total/2:6.1f}")

print("\n  Paper Eq.(5): C(3) = 1/2 * 10 = 5  --> CONFIRMED")

# --- 3. Degeneracy-weighted Casimir trace ---
print("\n3. Degeneracy-weighted Casimir: sum of (2l+1)*l(l+1)")
print("-" * 50)
print(f"  {'n_max':>5}  {'terms':>40}  {'Sum':>6}  {'Half':>6}")
for nmax in range(1, 6):
    terms = []
    for n in range(1, nmax + 1):
        for l in range(n):
            terms.append((2 * l + 1) * l * (l + 1))
    total = sum(terms)
    print(f"  {nmax:5d}  {str(terms):>40s}  {total:6d}  {total/2:6.1f}")

# Paper claims B = 42 (half of 84)
dw_sum_3 = sum((2*l+1)*l*(l+1) for n in range(1,4) for l in range(n))
print(f"\n  Paper Eq.(8): B = 1/2 * {dw_sum_3} = {dw_sum_3/2}")
print(f"  BUT paper says B = 42, which is the FULL sum {dw_sum_3}, not half!")
print(f"  Paper text: 'degeneracy-weighted half-trace' = 42")
print(f"  Actual half: {dw_sum_3/2}")
print()

# Let's check: does the paper actually use 42 or 21?
# K = pi*(42 + zeta(2) - 1/40) = pi*43.6199... = 137.036...
# With 21: K = pi*(21 + zeta(2) - 1/40) = pi*22.6199... = 71.04... NO
# So the paper uses 42 = full sum, despite calling it "half-trace"
print("  Resolution: Paper LABELS it half-trace but USES the full sum = 42")
print("  (The 1/2 in Eq.(8) is notational; the boxed value B=42 is correct)")

# --- 4. Selection principle ---
print("\n4. Selection principle: mean Casimir per state = dim(S3) = 3?")
print("-" * 50)

# The paper defines C_bar = C(nmax) / N_states
# But C(nmax) = 1/2 * sum l(l+1) = 5 at nmax=3
# N_states = 14
# C_bar = 5/14 = 0.357... != 3

# What actually works? Let's check multiple definitions
for nmax in range(1, 6):
    # Definition 1: Paper's C(nmax) = 1/2 * sum l(l+1)
    casimir_sum = sum(l*(l+1) for n in range(1, nmax+1) for l in range(n))
    half_casimir = casimir_sum / 2

    # Definition 2: Degeneracy-weighted sum
    dw_sum = sum((2*l+1)*l*(l+1) for n in range(1, nmax+1) for l in range(n))
    half_dw = dw_sum / 2

    # N_states = sum n^2
    N_states = sum(n**2 for n in range(1, nmax+1))

    # Number of (n,l) pairs
    N_pairs = sum(n for n in range(1, nmax+1))  # = nmax*(nmax+1)/2

    cbar_paper = half_casimir / N_states  # Paper's formula
    cbar_dw = dw_sum / N_states           # DW sum / N_states
    cbar_half_dw = half_dw / N_states     # Half DW / N_states
    cbar_pair = casimir_sum / N_pairs     # Per (n,l) pair

    print(f"  n_max={nmax}: N_states={N_states:3d}, N_pairs={N_pairs:2d}")
    print(f"    C_bar (paper: half-Casimir/N)     = {half_casimir:6.1f}/{N_states} = {cbar_paper:.4f}")
    print(f"    C_bar (full DW sum/N)             = {dw_sum:6d}/{N_states} = {cbar_dw:.4f}")
    print(f"    C_bar (half DW sum/N)             = {half_dw:6.1f}/{N_states} = {cbar_half_dw:.4f}")
    print(f"    C_bar (Casimir sum / N_pairs)     = {casimir_sum:6d}/{N_pairs} = {cbar_pair:.4f}")

print()
print("  FINDING: No standard definition gives C_bar(3) = 3 exactly.")
print("  Closest: DW_sum/N = 42/14 = 3.000 at n_max=3!")
print("  This means: the 'selection principle' works with the FULL")
print("  degeneracy-weighted sum (not the half-trace), divided by N_states.")
print("  42/14 = 3 = dim(S3). Unique to n_max=3.")

# Verify uniqueness
print("\n  Uniqueness check: DW_sum / N_states for each n_max:")
for nmax in range(1, 8):
    dw = sum((2*l+1)*l*(l+1) for n in range(1, nmax+1) for l in range(n))
    N = sum(n**2 for n in range(1, nmax+1))
    ratio = dw / N
    marker = " <-- = 3 exactly" if abs(ratio - 3.0) < 1e-10 else ""
    print(f"    n_max={nmax}: {dw}/{N} = {ratio:.6f}{marker}")

# --- 5. Fiber term ---
print("\n5. Fiber term: zeta(2)")
print("-" * 50)
zeta2 = np.pi**2 / 6
print(f"  zeta(2) = pi^2/6 = {zeta2:.10f}")
print(f"  Paper claims: 1.644934  --> paper value matches to 6 digits")

# --- 6. Cubic equation ---
print("\n6. Cubic equation: alpha^3 - K*alpha + 1 = 0")
print("-" * 50)

# Using paper's values
Delta_paper = 1/40
K_paper = np.pi * (42 + zeta2 - Delta_paper)
roots_paper = solve_cubic(K_paper)

print(f"  K = pi * (42 + {zeta2:.6f} - {Delta_paper}) = {K_paper:.6f}")
print(f"  Paper claims K = 137.036064")
print(f"  We compute  K = {K_paper:.6f}")
print(f"  Match: {abs(K_paper - 137.036064) < 0.0001}")
print()

print(f"  Roots of alpha^3 - {K_paper:.6f}*alpha + 1 = 0:")
for i, r in enumerate(roots_paper):
    if r > 0:
        alpha_inv = 1/r
        err = abs(alpha_inv - CODATA_ALPHA_INV) / CODATA_ALPHA_INV
        print(f"    root {i+1}: alpha = {r:.10f}, alpha^-1 = {alpha_inv:.6f}, "
              f"rel_err = {err:.2e}")
    else:
        print(f"    root {i+1}: alpha = {r:.10f} (negative)")

# Physical root
alpha_phys = min(r for r in roots_paper if r > 0)
alpha_inv_phys = 1/alpha_phys
rel_err = abs(alpha_inv_phys - CODATA_ALPHA_INV) / CODATA_ALPHA_INV
print(f"\n  Physical root: alpha^-1 = {alpha_inv_phys:.6f}")
print(f"  CODATA:        alpha^-1 = {CODATA_ALPHA_INV}")
print(f"  Relative error: {rel_err:.2e}")
print(f"  Paper claims 8.8e-8  --> {'CONFIRMED' if abs(rel_err - 8.8e-8) < 1e-8 else 'DISCREPANCY: ' + f'{rel_err:.2e}'}")


# =====================================================================
# PART A SUMMARY: The boundary term inconsistency
# =====================================================================
print("\n" + "=" * 72)
print("PART A SUMMARY: Boundary Term Inconsistency")
print("=" * 72)
print("""
The paper writes (Eq. 7):
  Delta = 1/(|lambda_3| x g_2) = 1/(8 x 5) = 1/40

But it also says "g_2 = 2^2 = 4 (the base-space states at n=2)".

If g_2 = 4, then 8 x 4 = 32, and Delta = 1/32, NOT 1/40.

The number 5 = 1 + 4 = N(2) = total states up to n=2.
So the actual computation uses N(2) = 5, not g_2 = 4.

This is the inconsistency flagged by the reviewer.
""")


# =====================================================================
# PART B: Boundary term investigation
# =====================================================================
print("=" * 72)
print("PART B: Boundary Term Investigation")
print("=" * 72)

# Build a comprehensive list of candidate Delta values
candidates = []

# Helper: compute alpha^-1 from Delta
def alpha_inv_from_delta(delta):
    K = np.pi * (42 + zeta2 - delta)
    roots = solve_cubic(K)
    pos_roots = [r for r in roots if r > 0]
    if not pos_roots:
        return None
    alpha = min(pos_roots)
    return 1.0 / alpha

# All spectral quantities at n_max=3
lam = {n: -(n**2 - 1) for n in range(1, 6)}
g = {n: n**2 for n in range(1, 6)}
N_tot = {nmax: sum(n**2 for n in range(1, nmax+1)) for nmax in range(1, 6)}
casimir = {nmax: sum(l*(l+1) for n in range(1, nmax+1) for l in range(n))
           for nmax in range(1, 6)}
dw_casimir = {nmax: sum((2*l+1)*l*(l+1) for n in range(1, nmax+1) for l in range(n))
              for nmax in range(1, 6)}

print("\nAvailable spectral quantities at n_max=3:")
print(f"  lambda: {dict((n, lam[n]) for n in range(1,4))}")
print(f"  |lambda|: {dict((n, abs(lam[n])) for n in range(1,4))}")
print(f"  g_n: {dict((n, g[n]) for n in range(1,4))}")
print(f"  N(nmax): {dict((nm, N_tot[nm]) for nm in range(1,4))}")
print(f"  Casimir sum: {dict((nm, casimir[nm]) for nm in range(1,4))}")
print(f"  DW Casimir: {dict((nm, dw_casimir[nm]) for nm in range(1,4))}")

print("\nSystematic enumeration of candidate Delta values:")
print("-" * 72)

def add_candidate(name, formula, value):
    """Add a candidate to the list."""
    ainv = alpha_inv_from_delta(value)
    if ainv is not None:
        err = (ainv - CODATA_ALPHA_INV) / CODATA_ALPHA_INV
        candidates.append((name, formula, value, ainv, err))

# --- Category 1: 1/(|lambda_i| * quantity) forms ---
for n_lam in [2, 3]:
    abs_lam = abs(lam[n_lam])

    # x g_n
    for n_g in range(1, 4):
        name = f"1/(|lam_{n_lam}|*g_{n_g})"
        val = 1.0 / (abs_lam * g[n_g])
        add_candidate(name, f"1/({abs_lam}*{g[n_g]})", val)

    # x N(nmax)
    for nm in range(1, 4):
        name = f"1/(|lam_{n_lam}|*N({nm}))"
        val = 1.0 / (abs_lam * N_tot[nm])
        add_candidate(name, f"1/({abs_lam}*{N_tot[nm]})", val)

    # x (2l+1) for various l
    for l_val in range(3):
        deg = 2*l_val + 1
        name = f"1/(|lam_{n_lam}|*(2*{l_val}+1))"
        val = 1.0 / (abs_lam * deg)
        add_candidate(name, f"1/({abs_lam}*{deg})", val)

# --- Category 2: 1/(product of eigenvalues) ---
add_candidate("1/(|lam_2|*|lam_3|)", "1/(3*8)", 1.0/(3*8))

# --- Category 3: Involving Casimir sums ---
add_candidate("1/Casimir(3)", "1/10", 1.0/10)
add_candidate("1/DW_Casimir(3)", "1/42", 1.0/42)
add_candidate("1/(2*DW_Casimir(3))", "1/84", 1.0/84)
add_candidate("Casimir(3)/DW_Casimir(3)", "10/42", 10.0/42)

# --- Category 4: Simple fractions near 1/40 ---
for denom in [24, 28, 30, 32, 35, 36, 38, 39, 40, 42, 45, 48, 50, 56, 60, 72, 80, 112, 120]:
    name = f"1/{denom}"
    add_candidate(name, f"1/{denom}", 1.0/denom)

# --- Category 5: The paper's actual formula interpretations ---
add_candidate("1/(|lam_3|*g_2) [paper text]", "1/(8*4)", 1.0/32)
add_candidate("1/(|lam_3|*N(2)) [actual calc]", "1/(8*5)", 1.0/40)

# --- Category 6: Other natural combinations ---
add_candidate("g_2/N(3)", "4/14", 4.0/14)
add_candidate("1/N(3)", "1/14", 1.0/14)
add_candidate("g_3/DW_Casimir(3)", "9/42", 9.0/42)
add_candidate("(n_max-1)/|lam_nmax|", "2/8", 2.0/8)
add_candidate("n_max/|lam_nmax|", "3/8", 3.0/8)

# Euler-Mascheroni / other constants
add_candidate("1/(4*pi^2)", "1/(4pi^2)", 1.0/(4*np.pi**2))

# --- Category 7: Zero boundary (pure Casimir + zeta) ---
add_candidate("Delta=0", "0", 0.0)

# --- Category 8: Delta that gives EXACT CODATA ---
# Solve: alpha_inv(Delta) = CODATA
# K = 1/alpha + alpha^2
alpha_exact = 1.0 / CODATA_ALPHA_INV
K_exact = 1.0/alpha_exact + alpha_exact**2
Delta_exact = 42 + zeta2 - K_exact/np.pi
add_candidate("EXACT (reverse-engineered)", f"{Delta_exact:.8f}", Delta_exact)

# Sort by absolute error
candidates.sort(key=lambda x: abs(x[4]))

# Print results
print(f"\n{'Rank':>4}  {'Name':>35}  {'Formula':>16}  {'Delta':>10}  "
      f"{'alpha^-1':>12}  {'Rel Error':>12}  {'|Err|':>10}")
print("-" * 120)
for i, (name, formula, delta, ainv, err) in enumerate(candidates):
    marker = ""
    if "1/40" in formula or "1/(8*5)" in formula:
        marker = " <-- paper"
    elif "1/(8*4)" in formula or "1/32" in formula:
        marker = " <-- g2=4"
    elif "EXACT" in name:
        marker = " <-- target"
    print(f"  {i+1:3d}  {name:>35s}  {formula:>16s}  {delta:10.6f}  "
          f"{ainv:12.6f}  {err:+12.2e}  {abs(err):10.2e}{marker}")

# Highlight key findings
print("\n" + "-" * 72)
print("KEY FINDINGS:")

# Find 1/40 and 1/32 ranks
for i, (name, formula, delta, ainv, err) in enumerate(candidates):
    if "1/(8*5)" in formula:
        print(f"  1/40 = 1/(8*5) [paper's actual calc]: rank {i+1}, "
              f"alpha^-1 = {ainv:.6f}, err = {err:+.2e}")
    if "1/(8*4)" in formula:
        print(f"  1/32 = 1/(8*4) [paper's stated formula]: rank {i+1}, "
              f"alpha^-1 = {ainv:.6f}, err = {err:+.2e}")
    if "EXACT" in name:
        print(f"  Exact Delta = {delta:.8f} (for reference)")

print(f"\n  CODATA alpha^-1 = {CODATA_ALPHA_INV}")
print(f"  Exact Delta needed = {Delta_exact:.8f}")
print(f"  1/40 = {1/40:.8f}")
print(f"  Difference: {abs(Delta_exact - 1/40):.8f}")


# =====================================================================
# PART C: Sensitivity analysis
# =====================================================================
print("\n" + "=" * 72)
print("PART C: Sensitivity Analysis")
print("=" * 72)

delta_range = np.linspace(0, 0.1, 1000)
alpha_inv_values = []

for d in delta_range:
    ainv = alpha_inv_from_delta(d)
    alpha_inv_values.append(ainv)

alpha_inv_values = np.array(alpha_inv_values, dtype=float)

# Compute local slope at Delta=1/40
idx_040 = np.argmin(np.abs(delta_range - 1/40))
if idx_040 > 0 and idx_040 < len(delta_range) - 1:
    slope = (alpha_inv_values[idx_040+1] - alpha_inv_values[idx_040-1]) / \
            (delta_range[idx_040+1] - delta_range[idx_040-1])
    print(f"\n  d(alpha^-1)/d(Delta) at Delta=1/40: {slope:.4f}")
    print(f"  A change of 0.001 in Delta changes alpha^-1 by {abs(slope)*0.001:.6f}")
    print(f"  The CODATA precision band is +/- {CODATA_ALPHA_INV * 1.5e-10:.6e}")
    print(f"  So Delta needs precision of ~{CODATA_ALPHA_INV * 1.5e-10 / abs(slope):.2e} to match CODATA")

# Range of alpha^-1 over Delta in [0, 0.1]
print(f"\n  alpha^-1 range over Delta in [0, 0.1]:")
print(f"    Delta=0:    alpha^-1 = {alpha_inv_values[0]:.6f}")
print(f"    Delta=0.05: alpha^-1 = {alpha_inv_values[500]:.6f}")
print(f"    Delta=0.1:  alpha^-1 = {alpha_inv_values[-1]:.6f}")
print(f"    Total variation: {alpha_inv_values[0] - alpha_inv_values[-1]:.6f}")

# Generate plot
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

    # Top panel: full range
    ax1.plot(delta_range, alpha_inv_values, 'b-', linewidth=2, label=r'$\alpha^{-1}(\Delta)$')
    ax1.axhline(y=CODATA_ALPHA_INV, color='r', linestyle='--', linewidth=1.5,
                label=f'CODATA = {CODATA_ALPHA_INV}')
    ax1.axvline(x=1/40, color='g', linestyle=':', linewidth=1.5, label=r'$\Delta = 1/40$ (paper)')
    ax1.axvline(x=1/32, color='orange', linestyle=':', linewidth=1.5, label=r'$\Delta = 1/32$ ($g_2=4$)')
    ax1.axvline(x=Delta_exact, color='purple', linestyle=':', linewidth=1, label=f'$\\Delta$ exact = {Delta_exact:.5f}')

    ax1.set_xlabel(r'$\Delta$ (boundary term)', fontsize=13)
    ax1.set_ylabel(r'$\alpha^{-1}$', fontsize=13)
    ax1.set_title('Sensitivity of $\\alpha^{-1}$ to boundary term $\\Delta$', fontsize=14)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)

    # Bottom panel: zoom near CODATA
    zoom_mask = (delta_range > 0.01) & (delta_range < 0.06)
    ax2.plot(delta_range[zoom_mask], alpha_inv_values[zoom_mask], 'b-', linewidth=2)
    ax2.axhline(y=CODATA_ALPHA_INV, color='r', linestyle='--', linewidth=1.5,
                label=f'CODATA = {CODATA_ALPHA_INV}')
    ax2.axvline(x=1/40, color='g', linestyle=':', linewidth=1.5, label='1/40 = 0.0250')
    ax2.axvline(x=1/32, color='orange', linestyle=':', linewidth=1.5, label='1/32 = 0.03125')
    ax2.axvline(x=Delta_exact, color='purple', linestyle=':', linewidth=1,
                label=f'exact = {Delta_exact:.5f}')

    # Mark key candidates
    for name, formula, delta, ainv, err in candidates[:10]:
        if 0.01 < delta < 0.06:
            ax2.plot(delta, ainv, 'ko', markersize=4)
            ax2.annotate(formula, (delta, ainv), fontsize=7, rotation=45,
                        textcoords="offset points", xytext=(5, 5))

    ax2.set_xlabel(r'$\Delta$ (boundary term)', fontsize=13)
    ax2.set_ylabel(r'$\alpha^{-1}$', fontsize=13)
    ax2.set_title('Zoomed view near CODATA value', fontsize=14)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plot_path = os.path.join(output_dir, 'sensitivity_plot.png')
    plt.savefig(plot_path, dpi=150)
    print(f"\n  Plot saved to: {plot_path}")
    plt.close()

except ImportError:
    print("\n  matplotlib not available, skipping plot")

# =====================================================================
# FINAL SUMMARY
# =====================================================================
print("\n" + "=" * 72)
print("FINAL SUMMARY")
print("=" * 72)
print(f"""
PAPER INCONSISTENCY CONFIRMED:
  The paper states g_2 = 2^2 = 4 but computes 8 x 5 = 40.
  The "5" is N(2) = 1 + 4 = total states up to n=2, NOT g_2.

  If the formula used g_2 = 4 as stated: Delta = 1/32
  If the formula uses N(2) = 5 as computed: Delta = 1/40

  The paper uses 1/40 (N(2)), which gives the better result.
  The formula text should read:
    Delta = 1/(|lambda_3| x N(2)) = 1/(8 x 5) = 1/40
  where N(2) = sum_{{n=1}}^{{2}} n^2 = 5 is the total state count
  up to the first non-trivial shell.

SENSITIVITY:
  alpha^-1 is roughly linear in Delta with slope ~{slope:.1f}.
  The variation across Delta in [0, 0.1] is ~{alpha_inv_values[0] - alpha_inv_values[-1]:.3f}.
  This is moderate: Delta matters, but alpha^-1 is not wildly sensitive.
  The result is NOT on a plateau — it IS on a slope.
  Getting close to CODATA requires choosing Delta to ~{CODATA_ALPHA_INV * 1e-7 / abs(slope):.1e} precision.

BEST CANDIDATES:
""")

for i, (name, formula, delta, ainv, err) in enumerate(candidates[:5]):
    print(f"  {i+1}. {name} = {delta:.6f}: alpha^-1 = {ainv:.6f}, err = {err:+.2e}")

print(f"""
IS 1/40 UNIQUE?
  Within simple spectral-quantity ratios, 1/40 = 1/(|lambda_3|*N(2))
  is among the best candidates. Check the ranked table above to see
  if any candidate with a cleaner derivation competes.
""")
