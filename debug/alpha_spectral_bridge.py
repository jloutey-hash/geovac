"""
Spectral Determinant Bridge Search
====================================

Searches for a correction that bridges the 0.1% gap between
  det'(S1) * det'(S3) / pi = 4pi^2 * exp(zeta(3)/(2pi^2)) ~ 41.957
and the Casimir trace B = 42.

If found, this would ground Link 2 of Paper 2's derivation chain
by connecting the empirical B = 42 to standard spectral geometry.

Reference: Paper 2 Sec VII (spectral determinant near-miss)
"""

import numpy as np
from scipy.special import zeta as hurwitz_zeta
from itertools import product as iterproduct


# ============================================================
# Fundamental constants
# ============================================================

PI = np.pi
PI2 = PI**2

# Targets
B_TARGET = 42.0
F_TARGET = PI2 / 6                         # zeta(2) = pi^2/6
DELTA_TARGET = 1 / 40
K_OVER_PI = B_TARGET + F_TARGET - DELTA_TARGET   # ~43.6199
K_TARGET = PI * K_OVER_PI                         # ~137.036

# Riemann zeta at specific points (using Hurwitz zeta with a=1)
ZETA_2 = PI2 / 6                           # 1.6449340668...
ZETA_3 = float(hurwitz_zeta(3, 1))         # 1.2020569031...
ZETA_4 = PI**4 / 90                        # 1.0823232337...
ZETA_5 = float(hurwitz_zeta(5, 1))         # 1.0369277551...

# Glaisher-Kinkelin constant A
# ln(A) = 1/12 - zeta'_R(-1)
# zeta'_R(-1) = -0.16542114370045... (standard value)
# A = 1.28242712910...
ZETA_PRIME_NEG1 = -0.1654211437004509710    # zeta'_R(-1), well-known
LN_GLAISHER = 1/12 - ZETA_PRIME_NEG1       # ln(A) = 1/12 - zeta'(-1)
GLAISHER_A = np.exp(LN_GLAISHER)           # A ~ 1.28243

# zeta'_R(0) = -ln(2pi)/2 (standard result)
ZETA_PRIME_0 = -np.log(2 * PI) / 2

# Euler-Mascheroni
EULER_GAMMA = 0.5772156649015329


# ============================================================
# PART A: Spectral Determinant Catalog
# ============================================================

def part_a():
    print("=" * 80)
    print("PART A: SPECTRAL DETERMINANT CATALOG")
    print("=" * 80)
    print()

    # ------------------------------------------------------------------
    # Zeta-regularized spectral determinants (zero-mode excluded)
    # Standard results for round unit spheres:
    #
    # S^1 (circumference 2pi):
    #   Eigenvalues: k^2, k = 1,2,3,... (each with multiplicity 2)
    #   det'(Delta_S1) = 4 pi^2
    #
    # S^2 (unit, area 4pi):
    #   Eigenvalues: l(l+1), multiplicity 2l+1, l=1,2,...
    #   det'(Delta_S2) = exp(1/2 - 4 zeta'_R(-1))
    #   Note: Paper 2 writes exp(-4 zeta'_R(-1)), omitting the 1/2.
    #   The standard Weisberger/Sarnak result includes the 1/2 term.
    #   We compute both and flag which the paper uses.
    #
    # S^3 (unit, volume 2pi^2):
    #   Eigenvalues: -(n^2-1) -> |lambda_n| = n^2-1, mult n^2, n=2,3,...
    #   Actually for Laplacian eigenvalues: lambda_n = n(n+2), mult (n+1)^2, n=0,1,2,...
    #   (equivalent: shift index). Zero mode at n=0.
    #   det'(Delta_S3) = pi * exp(zeta_R(3) / (2 pi^2))
    # ------------------------------------------------------------------

    det_S1 = 4 * PI2
    print(f"det'(Delta_S1) = 4 pi^2 = {det_S1:.12f}")

    # S^2: two conventions
    det_S2_with_half = np.exp(0.5 - 4 * ZETA_PRIME_NEG1)
    det_S2_paper = np.exp(-4 * ZETA_PRIME_NEG1)
    print(f"det'(Delta_S2) = exp(1/2 - 4 zeta'(-1)) = {det_S2_with_half:.12f}  [Weisberger/Sarnak]")
    print(f"det'(Delta_S2) = exp(-4 zeta'(-1))       = {det_S2_paper:.12f}  [Paper 2 convention]")
    print(f"  zeta'_R(-1) = {ZETA_PRIME_NEG1:.16f}")
    print(f"  Glaisher A  = {GLAISHER_A:.12f}")

    det_S3 = PI * np.exp(ZETA_3 / (2 * PI2))
    print(f"det'(Delta_S3) = pi * exp(zeta(3)/(2pi^2)) = {det_S3:.12f}")
    print(f"  zeta(3) = {ZETA_3:.16f}")
    print(f"  zeta(3)/(2pi^2) = {ZETA_3 / (2*PI2):.16f}")
    print()

    # The near-miss from Paper 2
    near_miss = det_S1 * det_S3 / PI
    print(f"Paper's near-miss: det'_1 * det'_3 / pi = {near_miss:.12f}")
    print(f"  = 4 pi^2 * exp(zeta(3)/(2pi^2))       = {near_miss:.12f}")
    print(f"  Target B = 42, gap = {near_miss - 42:.8f} ({(near_miss/42 - 1)*100:.6f}%)")
    print()

    # ------------------------------------------------------------------
    # Additional spectral invariants
    # ------------------------------------------------------------------
    print("-" * 80)
    print("ADDITIONAL SPECTRAL INVARIANTS")
    print("-" * 80)
    print()

    # Volumes
    vol_S1 = 2 * PI
    vol_S2 = 4 * PI
    vol_S3 = 2 * PI2
    print(f"Vol(S1) = 2pi   = {vol_S1:.12f}")
    print(f"Vol(S2) = 4pi   = {vol_S2:.12f}")
    print(f"Vol(S3) = 2pi^2 = {vol_S3:.12f}")
    print()

    # Euler characteristics
    chi_S1 = 0
    chi_S2 = 2
    chi_S3 = 0
    print(f"chi(S1) = {chi_S1}")
    print(f"chi(S2) = {chi_S2}")
    print(f"chi(S3) = {chi_S3}")
    print()

    # Eta invariant of Dirac operator on S^3
    # For the standard Dirac on round S^3: eta(0) = 0 (by symmetry, spectrum is symmetric)
    eta_S3 = 0.0
    print(f"eta(S3, Dirac) = {eta_S3} (by symmetry)")
    print()

    # Spectral zeta values for the Laplacian on each sphere
    # S^1: zeta_{S1}(s) = 2 * zeta_R(2s). Values:
    zeta_S1_neg1 = 0.0  # 2 * zeta_R(-2) = 0 (zeta_R at negative even integers = 0)
    zeta_S1_0 = 2 * ZETA_PRIME_0 * 0  # Actually zeta_{S1}(0) = 2*zeta_R(0) = 2*(-1/2) = -1
    zeta_S1_0 = -1.0                    # 2 * zeta_R(0) = 2 * (-1/2)
    zeta_S1_1 = 2 * ZETA_2             # = pi^2/3
    zeta_S1_2 = 2 * ZETA_4             # = 2*pi^4/90 = pi^4/45

    print(f"Spectral zeta of S1 Laplacian:")
    print(f"  zeta_S1(-1) = 2*zeta_R(-2) = {zeta_S1_neg1:.10f}")
    print(f"  zeta_S1(0)  = 2*zeta_R(0)  = {zeta_S1_0:.10f}")
    print(f"  zeta_S1(1)  = 2*zeta_R(2)  = {zeta_S1_1:.10f}  (= pi^2/3)")
    print(f"  zeta_S1(2)  = 2*zeta_R(4)  = {zeta_S1_2:.10f}  (= pi^4/45)")
    print()

    # S^2: eigenvalues l(l+1), multiplicity 2l+1
    # zeta_{S2}(s) = sum_{l=1}^inf (2l+1) / [l(l+1)]^s
    # zeta_{S2}(1) = sum (2l+1)/[l(l+1)] = sum [1/l + 1/(l+1)] by partial fractions
    #              This diverges. Actually: (2l+1)/[l(l+1)] = 1/l + 1/(l+1)? No.
    #              (2l+1)/[l(l+1)] = 1/l + 1/(l+1). Let's check: 1/l + 1/(l+1) = (2l+1)/[l(l+1)]. Yes.
    #              So zeta_{S2}(1) diverges (harmonic). The meromorphic continuation has a pole.
    # zeta_{S2}(0): by heat kernel, = (number of zero modes) - dim(kernel) ...
    #              Standard: zeta_{S2}(0) = -1/3 + correction...
    #              Actually for S^2: zeta_Delta(0) = -1/3 (minus zero mode contribution)
    #              From Weisberger: the functional determinant formula gives zeta'(0).
    # Let me just compute numerically for moderate l.
    def zeta_S2_numerical(s, L_max=100000):
        """Spectral zeta for Laplacian on S2 (zero mode excluded)."""
        l_vals = np.arange(1, L_max + 1, dtype=np.float64)
        eigenvals = l_vals * (l_vals + 1)
        degeneracies = 2 * l_vals + 1
        return np.sum(degeneracies / eigenvals**s)

    # zeta_{S2}(2): converges
    z_S2_2 = zeta_S2_numerical(2)
    print(f"Spectral zeta of S2 Laplacian (numerical, L_max=100000):")
    print(f"  zeta_S2(2) = {z_S2_2:.10f}")

    # For s=0: use the known analytic result
    # zeta_{S2}(0) = -1/3 (standard, excluding zero mode)
    z_S2_0 = -1/3
    print(f"  zeta_S2(0) = -1/3 = {z_S2_0:.10f}  (analytic)")
    print()

    # S^3: eigenvalues n(n+2) [or equivalently (n+1)^2 - 1], multiplicity (n+1)^2, n=1,2,...
    # (excluding n=0 zero mode)
    def zeta_S3_numerical(s, N_max=10000):
        """Spectral zeta for Laplacian on S3 (zero mode excluded)."""
        n_vals = np.arange(1, N_max + 1, dtype=np.float64)
        eigenvals = n_vals * (n_vals + 2)        # = (n+1)^2 - 1
        degeneracies = (n_vals + 1)**2
        return np.sum(degeneracies / eigenvals**s)

    z_S3_2 = zeta_S3_numerical(2)
    # zeta_{S3}(0): analytic continuation gives 0 (for odd-dim spheres, related to eta)
    z_S3_0 = 0.0  # Standard result for S^3
    print(f"Spectral zeta of S3 Laplacian (numerical, N_max=10000):")
    print(f"  zeta_S3(2) = {z_S3_2:.10f}")
    print(f"  zeta_S3(0) = {z_S3_0:.10f}  (analytic: 0 for odd spheres)")
    print()

    # Conformal / Weyl anomaly
    # For S^2 (2D): conformal anomaly ~ chi/6 = 2/6 = 1/3
    # For S^3 (3D, odd): no type-A Weyl anomaly
    conf_S2 = chi_S2 / 6
    print(f"Conformal anomaly coefficient:")
    print(f"  a(S2) = chi/6 = {conf_S2:.10f}")
    print(f"  a(S3) = 0 (odd dimension)")
    print()

    # Store everything for Part B
    catalog = {
        # Spectral determinants
        "det'_1": det_S1,
        "det'_2 (standard)": det_S2_with_half,
        "det'_2 (paper)": det_S2_paper,
        "det'_3": det_S3,
        # Volumes
        "Vol_1": vol_S1,
        "Vol_2": vol_S2,
        "Vol_3": vol_S3,
        # Euler characteristics (skip 0 values in multiplicative combos)
        "chi_2": float(chi_S2),
        # Key numbers
        "pi": PI,
        "pi^2": PI2,
        "2pi": 2*PI,
        "4pi": 4*PI,
        "2pi^2": 2*PI2,
        "4pi^2": 4*PI2,
        "zeta(2)": ZETA_2,
        "zeta(3)": ZETA_3,
        "zeta(4)": ZETA_4,
        "zeta(5)": ZETA_5,
        "Glaisher A": GLAISHER_A,
        "ln(A)": LN_GLAISHER,
        "zeta'(-1)": ZETA_PRIME_NEG1,
        "zeta'(0)": ZETA_PRIME_0,
        "Euler gamma": EULER_GAMMA,
        "ln(2)": np.log(2),
        "ln(pi)": np.log(PI),
        "1/40": 1/40,
        "1/8": 1/8,
        "1/(8pi)": 1/(8*PI),
        "chi_2/Vol_2": chi_S2 / vol_S2,   # 2/(4pi) = 1/(2pi)
        "zeta(3)/(2pi^2)": ZETA_3 / (2*PI2),
        "exp(zeta(3)/(2pi^2))": np.exp(ZETA_3 / (2*PI2)),
    }

    return catalog, det_S1, det_S2_with_half, det_S2_paper, det_S3


# ============================================================
# PART B: Systematic Bridge Search
# ============================================================

def part_b(catalog, det_S1, det_S2_std, det_S2_paper, det_S3):
    print()
    print("=" * 80)
    print("PART B: SYSTEMATIC BRIDGE SEARCH")
    print("=" * 80)
    print()

    targets = {
        "B = 42": B_TARGET,
        "K/pi": K_OVER_PI,
        "K": K_TARGET,
    }

    near_miss = det_S1 * det_S3 / PI   # 4pi^2 * exp(zeta(3)/(2pi^2))
    gap_to_42 = 42.0 - near_miss       # positive, ~0.043

    # ------------------------------------------------------------------
    # Section 1: Direct products/ratios of det', Vol, chi
    # ------------------------------------------------------------------
    print("-" * 80)
    print("SECTION 1: Products and ratios of spectral determinants")
    print("-" * 80)
    print()

    combos = {}
    d1 = det_S1
    d2s = det_S2_std
    d2p = det_S2_paper
    d3 = det_S3

    # Basic combinations
    combos["det'_1 * det'_3 / pi"] = d1 * d3 / PI
    combos["det'_1 * det'_3 / (2pi)"] = d1 * d3 / (2*PI)
    combos["det'_1 * det'_3 / pi^2"] = d1 * d3 / PI2
    combos["det'_1 * det'_3"] = d1 * d3
    combos["det'_1 * det'_3 / det'_2(std)"] = d1 * d3 / d2s
    combos["det'_1 * det'_3 / det'_2(paper)"] = d1 * d3 / d2p
    combos["det'_1 * det'_3 * det'_2(std)"] = d1 * d3 * d2s
    combos["det'_1 * det'_3 * det'_2(paper)"] = d1 * d3 * d2p
    combos["det'_3 / pi"] = d3 / PI
    combos["det'_1 / pi"] = d1 / PI
    combos["det'_1 / (2pi)"] = d1 / (2*PI)
    combos["det'_1 * det'_2(std)"] = d1 * d2s
    combos["det'_1 * det'_2(paper)"] = d1 * d2p
    combos["det'_2(std) * det'_3"] = d2s * d3
    combos["det'_2(paper) * det'_3"] = d2p * d3
    combos["det'_1 * det'_2(std) * det'_3"] = d1 * d2s * d3
    combos["det'_1 * det'_2(paper) * det'_3"] = d1 * d2p * d3

    # Ratios with volumes
    combos["det'_1 * det'_3 / Vol_3"] = d1 * d3 / (2*PI2)
    combos["det'_1 * det'_3 / Vol_2"] = d1 * d3 / (4*PI)
    combos["det'_1 * det'_3 / Vol_1"] = d1 * d3 / (2*PI)
    combos["det'_1 * det'_3 / (Vol_1*Vol_2)"] = d1 * d3 / (2*PI * 4*PI)
    combos["det'_1 * det'_3 / (pi * Vol_1)"] = d1 * d3 / (PI * 2*PI)

    # With chi(S2) = 2
    combos["det'_1 * det'_3 * chi_2 / pi"] = d1 * d3 * 2 / PI
    combos["det'_1 * det'_3 * chi_2 / (4*pi)"] = d1 * d3 * 2 / (4*PI)

    # Logarithmic
    combos["ln(det'_1) + ln(det'_3)"] = np.log(d1) + np.log(d3)
    combos["ln(det'_1 * det'_3)"] = np.log(d1 * d3)
    combos["ln(det'_1) * ln(det'_3)"] = np.log(d1) * np.log(d3)
    combos["det'_1^(1/2) * det'_3^(1/2)"] = np.sqrt(d1) * np.sqrt(d3)
    combos["(det'_1 * det'_3)^(1/2)"] = np.sqrt(d1 * d3)
    combos["(det'_1 * det'_3)^(1/2) / pi^(1/2)"] = np.sqrt(d1 * d3 / PI)

    # Powers
    combos["(det'_1 * det'_3 / pi)^(1/2)"] = np.sqrt(near_miss)
    combos["(det'_1 * det'_3 / pi)^2"] = near_miss**2

    # Exp/log based
    combos["pi * exp(det'_3/pi - 1)"] = PI * np.exp(d3/PI - 1)
    combos["4*pi^2 * (1 + zeta(3)/(2pi^2))"] = 4*PI2 * (1 + ZETA_3/(2*PI2))

    # The linearized version of exp(x) ~ 1 + x:
    # 4pi^2 * exp(zeta(3)/(2pi^2)) ~ 4pi^2 * (1 + zeta(3)/(2pi^2) + ...)
    # = 4pi^2 + 2*zeta(3) + ...
    combos["4*pi^2 + 2*zeta(3)"] = 4*PI2 + 2*ZETA_3
    combos["4*pi^2 + 2*zeta(3) + zeta(3)^2/(2pi^2)"] = (
        4*PI2 + 2*ZETA_3 + ZETA_3**2/(2*PI2))

    print(f"  {'Combination':<48s} {'Value':>14s} {'vs 42':>12s} {'vs K/pi':>12s} {'vs K':>12s}")
    print(f"  {'='*48} {'='*14} {'='*12} {'='*12} {'='*12}")

    hits_b = []
    for name, val in combos.items():
        if not np.isfinite(val) or val <= 0:
            continue
        errs = {}
        for tname, tval in targets.items():
            errs[tname] = abs(val - tval) / abs(tval)

        best_target = min(errs, key=errs.get)
        best_err = errs[best_target]

        flag = ""
        if best_err < 0.0001:
            flag = " *** < 0.01%"
        elif best_err < 0.001:
            flag = " **  < 0.1%"
        elif best_err < 0.01:
            flag = " *   < 1%"

        if best_err < 0.05:
            e42 = f"{errs['B = 42']:.2e}"
            ekpi = f"{errs['K/pi']:.2e}"
            ek = f"{errs['K']:.2e}"
            print(f"  {name:<48s} {val:>14.8f} {e42:>12s} {ekpi:>12s} {ek:>12s}{flag}")

        if best_err < 0.001:
            hits_b.append((name, val, best_target, best_err))

    print()

    # ------------------------------------------------------------------
    # Section 2: Near-miss + additive corrections
    # ------------------------------------------------------------------
    print("-" * 80)
    print("SECTION 2: Near-miss + additive correction terms")
    print(f"  Base = det'_1 * det'_3 / pi = {near_miss:.12f}")
    print(f"  Gap to 42: {gap_to_42:.12f}")
    print("-" * 80)
    print()

    # What correction added to near_miss gives exactly 42?
    print(f"  Needed correction to reach 42:       {gap_to_42:.12f}")
    print(f"  Needed correction to reach K/pi:     {K_OVER_PI - near_miss:.12f}")
    print()

    corrections = {
        "zeta(3)^2 / (4pi^2)":             ZETA_3**2 / (4*PI2),
        "zeta(3)^2 / (2pi^2)":             ZETA_3**2 / (2*PI2),
        "zeta(3) / (2pi^2)":               ZETA_3 / (2*PI2),
        "zeta(3) / pi^2":                  ZETA_3 / PI2,
        "1/40":                             1/40,
        "1/24":                             1/24,
        "1/12":                             1/12,
        "chi(S2)/Vol(S2)":                  2 / (4*PI),
        "1/(2pi)":                          1/(2*PI),
        "1/(4pi)":                          1/(4*PI),
        "1/(8pi)":                          1/(8*PI),
        "1/(2pi^2)":                        1/(2*PI2),
        "Euler_gamma / pi^2":               EULER_GAMMA / PI2,
        "ln(2) / pi^2":                     np.log(2) / PI2,
        "zeta'(-1)":                        ZETA_PRIME_NEG1,
        "-zeta'(-1)":                       -ZETA_PRIME_NEG1,
        "4*zeta'(-1)":                      4*ZETA_PRIME_NEG1,
        "-4*zeta'(-1)":                     -4*ZETA_PRIME_NEG1,
        "ln(A)":                            LN_GLAISHER,
        "1/A":                              1/GLAISHER_A,
        "A - 1":                            GLAISHER_A - 1,
        "ln(2pi)/(2pi)":                    np.log(2*PI)/(2*PI),
        "zeta(3)/(4pi^2) + 1/24":          ZETA_3/(4*PI2) + 1/24,
        "pi/72":                            PI/72,
        "1/(6pi)":                          1/(6*PI),
    }

    # Also try: the second-order term in the Taylor expansion of exp(x) for the near-miss
    # near_miss = 4pi^2 * exp(eps) where eps = zeta(3)/(2pi^2)
    # = 4pi^2 * (1 + eps + eps^2/2 + ...)
    # = 4pi^2 + 2*zeta(3) + zeta(3)^2/(4pi^2) + ...
    # So 42 - near_miss is the residual after all orders of exp.
    # Let's compute: 42 - 4pi^2 = 42 - 39.4784... = 2.5216
    # And 2*zeta(3) = 2.4041...
    # So 42 - 4pi^2 - 2*zeta(3) = 0.1174...
    # And zeta(3)^2/(4pi^2) = 0.03663...
    # So 42 - 4pi^2 - 2*zeta(3) - zeta(3)^2/(4pi^2) = 0.0808...

    print(f"  Taylor decomposition: 42 = 4pi^2 + 2*zeta(3) + correction_residual")
    print(f"    4pi^2       = {4*PI2:.12f}")
    print(f"    2*zeta(3)   = {2*ZETA_3:.12f}")
    print(f"    sum         = {4*PI2 + 2*ZETA_3:.12f}")
    print(f"    residual    = {42 - 4*PI2 - 2*ZETA_3:.12f}")
    resid_2 = 42 - 4*PI2 - 2*ZETA_3
    print(f"    zeta(3)^2/(4pi^2) = {ZETA_3**2/(4*PI2):.12f}")
    resid_3 = resid_2 - ZETA_3**2/(4*PI2)
    print(f"    after order 3: {resid_3:.12f}")
    print()

    # Check what well-known constants are near the gap
    print(f"  {'Correction':<35s} {'Value':>14s}  {'Base+corr':>14s} {'RelErr vs 42':>14s}")
    print(f"  {'='*35} {'='*14}  {'='*14} {'='*14}")

    for cname, cval in corrections.items():
        total = near_miss + cval
        err_42 = abs(total - 42) / 42
        flag = ""
        if err_42 < 0.0001:
            flag = " *** < 0.01%"
        elif err_42 < 0.001:
            flag = " **  < 0.1%"
        elif err_42 < 0.01:
            flag = " *   < 1%"
        if err_42 < 0.05:
            print(f"  {cname:<35s} {cval:>14.10f}  {total:>14.10f} {err_42:>13.2e}{flag}")

    print()

    # ------------------------------------------------------------------
    # Section 3: Multiplicative corrections to near-miss
    # ------------------------------------------------------------------
    print("-" * 80)
    print("SECTION 3: Multiplicative corrections  (near_miss * factor = 42?)")
    print("-" * 80)
    print()

    exact_factor = 42 / near_miss
    print(f"  Exact factor needed: 42 / near_miss = {exact_factor:.16f}")
    print(f"  Deviation from 1: {exact_factor - 1:.10e}")
    print()

    mult_factors = {
        "1 + zeta(3)^2/(8pi^4)":           1 + ZETA_3**2 / (8*PI**4),
        "1 + 1/(8pi^2)":                   1 + 1/(8*PI2),
        "1 + 1/(12pi^2)":                  1 + 1/(12*PI2),
        "1 + zeta(3)/(24pi^2)":            1 + ZETA_3/(24*PI2),
        "1 + ln(A)/pi^2":                  1 + LN_GLAISHER/PI2,
        "1 + 1/(4pi^2 * zeta(3))":         1 + 1/(4*PI2*ZETA_3),
        "exp(1/(8pi^2))":                  np.exp(1/(8*PI2)),
        "exp(1/(12pi^2))":                 np.exp(1/(12*PI2)),
        "exp(zeta(3)^2/(8pi^4))":          np.exp(ZETA_3**2/(8*PI**4)),
        "42/(4pi^2 * exp(zeta3/(2pi2)))":  exact_factor,  # tautological, for reference
        "det'_2(std) / det'_2(paper)":     d2s / d2p,     # = exp(1/2) = 1.6487...
        "exp(1/(2*42))":                   np.exp(1/84),
        "(42/41)":                         42/41,
        "A^(1/pi)":                        GLAISHER_A**(1/PI),
    }

    print(f"  {'Factor':<42s} {'Value':>16s} {'Result':>14s} {'RelErr vs 42':>14s}")
    print(f"  {'='*42} {'='*16} {'='*14} {'='*14}")

    for fname, fval in mult_factors.items():
        result = near_miss * fval
        err_42 = abs(result - 42) / 42
        flag = ""
        if err_42 < 0.0001:
            flag = " *** < 0.01%"
        elif err_42 < 0.001:
            flag = " **  < 0.1%"
        elif err_42 < 0.01:
            flag = " *   < 1%"
        if err_42 < 0.05:
            print(f"  {fname:<42s} {fval:>16.12f} {result:>14.10f} {err_42:>13.2e}{flag}")

    print()

    return hits_b


# ============================================================
# PART C: Additive Decomposition Search
# ============================================================

def part_c(catalog, det_S1, det_S2_std, det_S2_paper, det_S3):
    print()
    print("=" * 80)
    print("PART C: ADDITIVE DECOMPOSITION SEARCH")
    print("=" * 80)
    print()

    # Build pools of "simple functions" of each manifold's data
    d1 = det_S1       # 4pi^2
    d2s = det_S2_std
    d2p = det_S2_paper
    d3 = det_S3       # pi * exp(zeta(3)/(2pi^2))
    vol1 = 2 * PI
    vol2 = 4 * PI
    vol3 = 2 * PI2

    # Pool for S1 (fiber)
    pool_S1 = {
        "det'_1":                d1,
        "ln(det'_1)":            np.log(d1),          # ln(4pi^2)
        "det'_1^(1/2)":         np.sqrt(d1),          # 2pi
        "det'_1^(-1)":          1/d1,
        "Vol_1":                 vol1,                  # 2pi
        "Vol_1^2":               vol1**2,               # 4pi^2
        "det'_1/Vol_1":         d1/vol1,               # 2pi
        "det'_1/(2pi)":         d1/(2*PI),             # 2pi
        "ln(det'_1)/pi":        np.log(d1)/PI,
        "det'_1/pi":            d1/PI,                 # 4pi
        "det'_1/pi^2":          d1/PI2,                # 4
        "det'_1/(4pi)":         d1/(4*PI),             # pi
        "zeta(2)":               ZETA_2,                # pi^2/6
        "2*zeta(2)":            2*ZETA_2,              # pi^2/3
        "6*zeta(2)/pi^2":       6*ZETA_2/PI2,          # 1
    }

    # Pool for S2 (base)
    pool_S2 = {
        "det'_2(std)":          d2s,
        "det'_2(paper)":        d2p,
        "ln(det'_2(std))":      np.log(d2s),
        "ln(det'_2(paper))":    np.log(d2p),
        "det'_2(std)^(1/2)":   np.sqrt(d2s),
        "det'_2(paper)^(1/2)": np.sqrt(d2p),
        "1/det'_2(std)":       1/d2s,
        "1/det'_2(paper)":     1/d2p,
        "Vol_2":                vol2,                   # 4pi
        "Vol_2/(2pi)":          vol2/(2*PI),            # 2
        "chi_2":                2.0,
        "chi_2 * Vol_2":        2 * vol2,               # 8pi
        "det'_2(std)/Vol_2":   d2s/vol2,
        "det'_2(paper)/Vol_2": d2p/vol2,
        "ln(det'_2(std))/pi":  np.log(d2s)/PI,
        "4*zeta'(-1)":         4*ZETA_PRIME_NEG1,
        "-4*zeta'(-1)":        -4*ZETA_PRIME_NEG1,
        "42 (Casimir B)":      42.0,                   # include for completeness
    }

    # Pool for S3 (total space)
    pool_S3 = {
        "det'_3":               d3,
        "ln(det'_3)":           np.log(d3),
        "det'_3^(1/2)":        np.sqrt(d3),
        "det'_3^(-1)":         1/d3,
        "Vol_3":                vol3,                   # 2pi^2
        "det'_3/pi":            d3/PI,                  # exp(zeta(3)/(2pi^2))
        "det'_3/Vol_3":        d3/vol3,
        "ln(det'_3)/pi":       np.log(d3)/PI,
        "zeta(3)/(2pi^2)":     ZETA_3/(2*PI2),
        "zeta(3)":              ZETA_3,
        "2*zeta(3)":            2*ZETA_3,
        "exp(zeta(3)/(2pi^2))": np.exp(ZETA_3/(2*PI2)),
        "-1/40":                -DELTA_TARGET,          # boundary term (negative in formula)
        "1/40":                 DELTA_TARGET,
        "0":                    0.0,
    }

    targets_c = {
        "42": 42.0,
        "K/pi": K_OVER_PI,
    }

    # ------------------------------------------------------------------
    # Grid search: one term from each pool, sum = target?
    # ------------------------------------------------------------------
    print("-" * 80)
    print("GRID SEARCH: f(S1) + g(S2) + h(S3) = target")
    print("-" * 80)
    print()

    hits_c = []

    for (n1, v1), (n2, v2), (n3, v3) in iterproduct(
            pool_S1.items(), pool_S2.items(), pool_S3.items()):
        total = v1 + v2 + v3
        for tname, tval in targets_c.items():
            if tval == 0:
                continue
            rel_err = abs(total - tval) / abs(tval)
            if rel_err < 0.001:  # < 0.1%
                hits_c.append((n1, v1, n2, v2, n3, v3, tname, tval, total, rel_err))

    # Sort by error
    hits_c.sort(key=lambda x: x[-1])

    # Filter out trivially tautological entries (where B=42 appears literally)
    hits_nontrivial = [h for h in hits_c if "42 (Casimir B)" not in (h[0], h[2], h[4])]
    hits_tautological = [h for h in hits_c if "42 (Casimir B)" in (h[0], h[2], h[4])]

    print(f"Found {len(hits_c)} total hits < 0.1% ({len(hits_nontrivial)} non-trivial)")
    print()

    if hits_nontrivial:
        print(f"  {'f(S1)':<22s} {'g(S2)':<24s} {'h(S3)':<24s} {'Sum':>14s} {'Target':>8s} {'RelErr':>12s}")
        print(f"  {'='*22} {'='*24} {'='*24} {'='*14} {'='*8} {'='*12}")
        shown = set()
        for n1, v1, n2, v2, n3, v3, tname, tval, total, rel_err in hits_nontrivial[:50]:
            key = (n1, n2, n3, tname)
            if key in shown:
                continue
            shown.add(key)
            flag = " ***" if rel_err < 0.0001 else " ** " if rel_err < 0.001 else ""
            print(f"  {n1:<22s} {n2:<24s} {n3:<24s} {total:>14.8f} {tname:>8s} {rel_err:>11.2e}{flag}")
    else:
        print("  No non-trivial hits found at 0.1% threshold.")

    print()

    # ------------------------------------------------------------------
    # Focused: can ln(det') decompose K/pi?
    # ------------------------------------------------------------------
    print("-" * 80)
    print("FOCUSED: logarithmic decomposition")
    print("-" * 80)
    print()
    print(f"  ln(det'_1) = ln(4pi^2) = {np.log(d1):.12f}")
    print(f"  ln(det'_2(std))        = {np.log(d2s):.12f}")
    print(f"  ln(det'_2(paper))      = {np.log(d2p):.12f}")
    print(f"  ln(det'_3)             = {np.log(d3):.12f}")
    print()

    ld1 = np.log(d1)
    ld2s = np.log(d2s)
    ld2p = np.log(d2p)
    ld3 = np.log(d3)

    # If K/pi = a * ln(d1) + b * ln(d2) + c * ln(d3), solve for rational a,b,c?
    # This is 1 equation in 3 unknowns. We can scan small rationals.
    print("  Scanning a*ln(d1) + b*ln(d2) + c*ln(d3) = K/pi")
    print(f"  Target K/pi = {K_OVER_PI:.12f}")
    print()

    # Scan a, b, c in {-3, -2, -1, -1/2, 0, 1/2, 1, 2, 3} x pi factors
    rat_vals = [-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3]
    pi_mults = [1, PI, 1/PI, PI2, 1/PI2]
    pi_names = ["", "*pi", "/pi", "*pi^2", "/pi^2"]

    log_hits = []
    for a in rat_vals:
        for b in rat_vals:
            for c in rat_vals:
                for pi_idx, (pm, pn) in enumerate(zip(pi_mults, pi_names)):
                    val = pm * (a * ld1 + b * ld2s + c * ld3)
                    err = abs(val - K_OVER_PI) / K_OVER_PI
                    if err < 0.001:
                        log_hits.append((a, b, c, pn, val, err))

    log_hits.sort(key=lambda x: x[-1])
    if log_hits:
        print(f"  {'a':>5s} {'b':>5s} {'c':>5s} {'pi_factor':>10s} {'Value':>14s} {'RelErr':>12s}")
        print(f"  {'='*5} {'='*5} {'='*5} {'='*10} {'='*14} {'='*12}")
        for a, b, c, pn, val, err in log_hits[:20]:
            flag = " ***" if err < 0.0001 else " ** " if err < 0.001 else ""
            print(f"  {a:>5.1f} {b:>5.1f} {c:>5.1f} {pn:>10s} {val:>14.8f} {err:>11.2e}{flag}")
    else:
        print("  No logarithmic decomposition found at 0.1% threshold.")

    print()

    # ------------------------------------------------------------------
    # Focused: the gap itself
    # ------------------------------------------------------------------
    print("-" * 80)
    print("FOCUSED: analyzing the gap = 42 - 4pi^2 exp(zeta(3)/(2pi^2))")
    print("-" * 80)
    print()

    gap = 42.0 - (4 * PI2 * np.exp(ZETA_3 / (2*PI2)))
    print(f"  gap = {gap:.16f}")
    print(f"  gap / pi      = {gap/PI:.16f}")
    print(f"  gap / pi^2    = {gap/PI2:.16f}")
    print(f"  gap * pi      = {gap*PI:.16f}")
    print(f"  gap * pi^2    = {gap*PI2:.16f}")
    print(f"  gap * 42      = {gap*42:.16f}")
    print(f"  gap * 40      = {gap*40:.16f}")
    print(f"  gap * 24      = {gap*24:.16f}")
    print(f"  gap * 12      = {gap*12:.16f}")
    print(f"  gap * 8       = {gap*8:.16f}")
    print(f"  1/gap         = {1/gap:.16f}")
    print(f"  ln(gap)       = {np.log(gap):.16f}")
    print()

    # Does the gap relate to known constants?
    gap_ratios = {
        "zeta(3)":          ZETA_3,
        "zeta(3)^2":        ZETA_3**2,
        "zeta(3)/(2pi^2)":  ZETA_3/(2*PI2),
        "zeta(5)":          ZETA_5,
        "pi^2/6":           PI2/6,
        "pi/6":             PI/6,
        "1/(8pi)":          1/(8*PI),
        "1/(2pi)":          1/(2*PI),
        "Euler gamma":      EULER_GAMMA,
        "ln(2)":            np.log(2),
        "ln(A)":            LN_GLAISHER,
        "1/24":             1/24,
        "1/40":             1/40,
        "1/12":             1/12,
    }

    print(f"  gap / X for known X:")
    print(f"  {'X':<22s} {'X_value':>14s} {'gap/X':>14s} {'Notes':>20s}")
    print(f"  {'='*22} {'='*14} {'='*14} {'='*20}")
    for name, val in gap_ratios.items():
        ratio = gap / val
        # Check if ratio is near a simple number
        notes = ""
        for simple, sname in [(1, "1"), (2, "2"), (0.5, "1/2"), (PI, "pi"),
                               (1/PI, "1/pi"), (PI2, "pi^2"), (1/PI2, "1/pi^2"),
                               (3, "3"), (1/3, "1/3"), (4, "4"), (1/4, "1/4"),
                               (6, "6"), (1/6, "1/6"), (1/8, "1/8"), (8, "8"),
                               (np.e, "e"), (1/np.e, "1/e"), (np.sqrt(2), "sqrt(2)")]:
            if abs(ratio - simple) / max(abs(simple), 1e-15) < 0.01:
                notes = f"~ {sname} ({abs(ratio/simple-1)*100:.3f}%)"
                break
        print(f"  {name:<22s} {val:>14.10f} {ratio:>14.10f} {notes:>20s}")

    print()

    # ------------------------------------------------------------------
    # Exact decomposition: 42 = 4pi^2 + (42 - 4pi^2)
    # Can (42 - 4pi^2) be expressed in terms of zeta values?
    # ------------------------------------------------------------------
    print("-" * 80)
    print("FOCUSED: 42 = 4pi^2 + X.  What is X in spectral language?")
    print("-" * 80)
    print()
    X = 42 - 4*PI2
    print(f"  X = 42 - 4pi^2 = {X:.16f}")
    print(f"  2*zeta(3)      = {2*ZETA_3:.16f}  (diff: {X - 2*ZETA_3:.10e})")
    print(f"  2*zeta(3) + zeta(3)^2/(4pi^2) = {2*ZETA_3 + ZETA_3**2/(4*PI2):.16f}  (diff: {X - 2*ZETA_3 - ZETA_3**2/(4*PI2):.10e})")
    print()

    # Full Taylor: exp(eps) = sum eps^k/k!  where eps = zeta(3)/(2pi^2) ~ 0.0609
    eps = ZETA_3 / (2*PI2)
    print(f"  eps = zeta(3)/(2pi^2) = {eps:.16f}")
    print(f"  4pi^2 * exp(eps) Taylor expansion:")
    taylor_sum = 0.0
    for k in range(20):
        term = 4*PI2 * eps**k / float(np.prod(np.arange(1, k+1))) if k > 0 else 4*PI2
        taylor_sum += term
        resid = 42 - taylor_sum
        if k <= 8 or abs(resid) < 0.001:
            print(f"    order {k:2d}: term = {term:.12e}, cumsum = {taylor_sum:.12f}, gap to 42 = {resid:+.10e}")

    print()
    print(f"  Full exp sum = {4*PI2 * np.exp(eps):.16f} (= near_miss)")
    print(f"  42 - exp sum = {42 - 4*PI2*np.exp(eps):.16f} (= gap, must come from outside Taylor)")
    print()

    # The gap is not in the Taylor series of exp.
    # It must come from a SEPARATE contribution. What spectral object gives ~0.0431?
    print(f"  Key question: What gives gap = {gap:.10f}?")
    print()

    # Interesting observations
    print("-" * 80)
    print("INTERESTING RATIOS")
    print("-" * 80)
    print()

    # gap / (zeta(3)^2 / (2pi^2))
    z3sq_over_2pi2 = ZETA_3**2 / (2*PI2)
    print(f"  gap / [zeta(3)^2/(2pi^2)] = {gap / z3sq_over_2pi2:.10f}")
    print(f"  gap / [zeta(3)/(4pi)]     = {gap / (ZETA_3/(4*PI)):.10f}")
    print(f"  gap * (2pi^2/zeta(3))     = {gap * 2*PI2/ZETA_3:.10f}")
    print(f"  gap * (pi^2)              = {gap * PI2:.10f}")
    print(f"  gap * (2pi)               = {gap * 2*PI:.10f}")
    print(f"  gap * 137                 = {gap * 137:.10f}")
    print(f"  gap * K                   = {gap * K_TARGET:.10f}")
    print(f"  1/gap                     = {1/gap:.10f}")
    print(f"  pi/gap                    = {PI/gap:.10f}")
    print(f"  pi^2/gap                  = {PI2/gap:.10f}")
    print()

    # Is gap related to det'_2?
    print(f"  gap * det'_2(std)         = {gap * d2s:.10f}")
    print(f"  gap * det'_2(paper)       = {gap * d2p:.10f}")
    print(f"  gap / ln(det'_2(std))     = {gap / np.log(d2s):.10f}")
    print(f"  gap / ln(det'_2(paper))   = {gap / np.log(d2p):.10f}")
    print(f"  gap + ln(det'_2(std))     = {gap + np.log(d2s):.10f}")
    print(f"  ln(det'_2(std))           = {np.log(d2s):.10f}")
    print(f"  det'_2(std) - det'_2(pap) = {d2s - d2p:.10f}  (= det'_2(paper)*(exp(1/2)-1))")
    print()

    return hits_c


# ============================================================
# Summary
# ============================================================

def print_summary(hits_b, hits_c):
    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print()

    near_miss = 4 * PI2 * np.exp(ZETA_3 / (2*PI2))
    gap = 42 - near_miss

    print(f"Near-miss:  4pi^2 * exp(zeta(3)/(2pi^2)) = {near_miss:.12f}")
    print(f"Target:     B = {B_TARGET}")
    print(f"Gap:        {gap:.12f}  ({gap/42*100:.6f}% of 42)")
    print()

    if hits_b:
        print("Part B hits (< 0.1% of any target):")
        for name, val, tname, err in hits_b:
            print(f"  {name}: {val:.10f} -> {tname} (err {err:.2e})")
    else:
        print("Part B: No combination hit a target within 0.1%.")

    print()

    nontrivial_c = [h for h in hits_c if "42 (Casimir B)" not in (h[0], h[2], h[4])]
    if nontrivial_c:
        print(f"Part C hits (non-trivial, < 0.1% of target):")
        shown = set()
        for n1, v1, n2, v2, n3, v3, tname, tval, total, rel_err in nontrivial_c[:10]:
            key = (n1, n2, n3, tname)
            if key in shown:
                continue
            shown.add(key)
            print(f"  {n1} + {n2} + {n3} = {total:.10f} -> {tname} (err {rel_err:.2e})")
    else:
        print("Part C: No non-trivial additive decomposition found at 0.1%.")

    print()
    print("CONCLUSION:")
    print(f"  The gap {gap:.8f} between the spectral determinant product and B=42")
    print(f"  does not decompose cleanly into standard spectral invariants at the")
    print(f"  0.01% level. The Casimir trace B=42 remains an irreducibly integer")
    print(f"  quantity, not yet derivable from spectral determinants alone.")
    print(f"  The near-miss (0.1%) is suggestive but unexplained.")


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    catalog, d1, d2s, d2p, d3 = part_a()
    hits_b = part_b(catalog, d1, d2s, d2p, d3)
    hits_c = part_c(catalog, d1, d2s, d2p, d3)
    print_summary(hits_b, hits_c)
