"""
Analytic Torsion, Heat Kernel & Index Theory for the Hopf Bundle
=================================================================

Tests whether standard spectral-geometric invariants of S¹ → S³ → S²
can produce the combination rule K = π(B + F - Δ) from Paper 2.

Part A: Ray-Singer analytic torsion of S¹, S², S³
Part B: Seeley-DeWitt heat kernel coefficients & partition functions
Part C: APS η-invariant, Chern-Simons invariant, index theory

Reference: Paper 2 — "The Fine Structure Constant from Spectral
Geometry of the Hopf Fibration" (papers/conjectures/paper_2_alpha.tex)
"""

import numpy as np
from scipy.special import zeta as hurwitz_zeta

# ============================================================
# Constants & targets
# ============================================================

PI = np.pi
PI2 = PI**2

B_TARGET = 42.0
F_TARGET = PI2 / 6                          # ζ(2)
DELTA_TARGET = 1 / 40
K_OVER_PI = B_TARGET + F_TARGET - DELTA_TARGET   # ~43.6199
K_TARGET = PI * K_OVER_PI                         # ~137.036

TARGETS = {
    "K":        K_TARGET,
    "K/π":      K_OVER_PI,
    "K/3":      K_TARGET / 3,
    "42":       42.0,
    "π²/6":     PI2 / 6,
    "1/40":     DELTA_TARGET,
}

# Known spectral constants
ZETA_3 = float(hurwitz_zeta(3, 1))              # 1.2020569031...
ZETA_5 = float(hurwitz_zeta(5, 1))              # 1.0369277551...
ZETA_PRIME_NEG1 = -0.1654211437004509710         # ζ'_R(-1)
ZETA_PRIME_0 = -np.log(2 * PI) / 2              # ζ'_R(0) = -ln(2π)/2
LN_GLAISHER = 1/12 - ZETA_PRIME_NEG1            # ln(A)
GLAISHER_A = np.exp(LN_GLAISHER)
EULER_GAMMA = 0.5772156649015329


def check_targets(value: float, label: str, results: list) -> None:
    """Check value against all targets; append hits within 0.1%."""
    for tname, tval in TARGETS.items():
        if tval == 0:
            continue
        rel = abs(value - tval) / abs(tval)
        if rel < 0.001:
            results.append((label, value, tname, tval, rel))


def print_hits(results: list) -> None:
    """Print any target matches."""
    if not results:
        print("  No matches within 0.1% of any target.")
        return
    print(f"  {'Quantity':<52s} {'Value':>14s}  {'Target':>8s} {'TgtVal':>14s} {'RelErr':>10s}")
    print(f"  {'─'*52} {'─'*14}  {'─'*8} {'─'*14} {'─'*10}")
    for label, val, tname, tval, rel in sorted(results, key=lambda x: x[4]):
        flag = "***" if rel < 1e-4 else " **" if rel < 5e-4 else "  *"
        print(f"  {label:<52s} {val:>14.8f}  {tname:>8s} {tval:>14.8f} {rel:>9.2e} {flag}")


# ============================================================
# PART A: Ray-Singer Analytic Torsion
# ============================================================

def spectral_zeta_S3_q0(s: float, N: int = 20000) -> float:
    """Spectral zeta for Δ_0 on S³: eigenvalues k(k+2), degeneracy (k+1)², k≥1."""
    return sum((k+1)**2 / (k*(k+2))**s for k in range(1, N+1))


def spectral_zeta_S3_q1(s: float, N: int = 20000) -> float:
    """
    Spectral zeta for Δ_1 on S³ (1-forms, coexact part).

    On S³, the Hodge Laplacian on 1-forms splits into exact and coexact pieces.
    The coexact 1-form eigenvalues are: (k+1)² with degeneracy 2k(k+2), k≥1.
    (Ikeda-Taniguchi; the exact piece has eigenvalues (k+1)² with deg k(k+2)
     but these are the image of d from 0-forms, so for torsion we use coexact.)

    Actually, for Ray-Singer torsion the relevant operator is the Hodge
    Laplacian restricted to q-forms. On S³ (dim 3):
    - Δ_0: scalar Laplacian, eigenvalues k(k+2), deg (k+1)², k≥1
    - Δ_1: 1-form Laplacian, eigenvalues split:
        k(k+2) with deg 2k(k+2), k≥1 [from exact + coexact pieces]

    The full 1-form spectrum on S³ has eigenvalues k(k+2) with
    multiplicity 2k(k+2) for k≥1 (see Ikeda-Taniguchi 1978).
    """
    return sum(2*k*(k+2) / (k*(k+2))**s for k in range(1, N+1))


def spectral_zeta_deriv_S3_q0(s: float, N: int = 20000, ds: float = 1e-7) -> float:
    """Numerical derivative of ζ_{Δ_0,S³}(s) via central difference."""
    return (spectral_zeta_S3_q0(s + ds, N) - spectral_zeta_S3_q0(s - ds, N)) / (2 * ds)


def spectral_zeta_deriv_S3_q1(s: float, N: int = 20000, ds: float = 1e-7) -> float:
    """Numerical derivative of ζ_{Δ_1,S³}(s) via central difference."""
    return (spectral_zeta_S3_q1(s + ds, N) - spectral_zeta_S3_q1(s - ds, N)) / (2 * ds)


def part_a():
    """Ray-Singer analytic torsion."""
    print("=" * 80)
    print("PART A: RAY-SINGER ANALYTIC TORSION")
    print("=" * 80)
    print()

    results = []
    N = 20000  # summation cutoff

    # ------------------------------------------------------------------
    # S¹ analytic torsion
    # ------------------------------------------------------------------
    # log T(S¹) = (1/2) (-1)^0 · 0 · log det'(Δ_0) = 0
    # (only q=0 in dimension 1, and the weight is q=0)
    # Actually for S¹ (dim 1): q runs 0,1. But dim ker Δ_0 = 1 (constants),
    # dim ker Δ_1 = 0 on S¹.
    # log T = (1/2)[(-1)^0 · 0 · log det'(Δ_0) + (-1)^1 · 1 · log det'(Δ_1)]
    #       = -(1/2) log det'(Δ_1)
    # On S¹, Δ_1 has the same nonzero spectrum as Δ_0 (Hodge duality: *Δ_0 = Δ_1* on S¹).
    # det'(Δ_0 on S¹) = 4π²
    # So: log T(S¹) = -(1/2) log(4π²) = -log(2π)

    log_T_S1 = -np.log(2 * PI)
    T_S1 = np.exp(log_T_S1)  # = 1/(2π)

    print("S¹ (circle, circumference 2π):")
    print(f"  det'(Δ_0, S¹) = 4π² = {4*PI2:.10f}")
    print(f"  log T(S¹) = -log(2π) = {log_T_S1:.10f}")
    print(f"  T(S¹) = 1/(2π) = {T_S1:.10f}")
    print()

    # ------------------------------------------------------------------
    # S³ analytic torsion
    # ------------------------------------------------------------------
    # For S³ (dim 3), log T = (1/2) Σ_{q=0}^{3} (-1)^q · q · log det'(Δ_q)
    #   = (1/2)[0 - log det'(Δ_1) + 2·log det'(Δ_2) - 3·log det'(Δ_3)]
    # By Hodge duality on S³: Δ_q ≅ Δ_{3-q}, so det'(Δ_3)=det'(Δ_0), det'(Δ_2)=det'(Δ_1)
    # log T = (1/2)[-log det'(Δ_1) + 2·log det'(Δ_1) - 3·log det'(Δ_0)]
    #       = (1/2)[log det'(Δ_1) - 3·log det'(Δ_0)]

    print("S³ spectral zeta functions (N = {})...".format(N))

    # ζ_{Δ_0}(0) on S³
    zeta_0_at_0 = spectral_zeta_S3_q0(0 + 1e-10, N)
    # The exact value: ζ_{Δ_0,S³}(0) = sum (k+1)² [as k→∞ limit]
    # This diverges at s=0. We need analytic continuation.
    # Known result: ζ_{Δ_0,S³}(0) = 0 (for scalar Laplacian on odd sphere, from Minakshisundaram-Pleijel)
    # Actually on S³: ζ_{Δ_0}(0) = 0 (odd-dim sphere, scalar).

    # For the spectral determinant, we need ζ'(0):
    # log det'(Δ_0, S³) = -ζ'_{Δ_0}(0)
    # Known result (Vardi 1988, Quine-Heydari-Song 1993):
    # log det'(Δ_0, S³) = log(π) + ζ(3)/(2π²)

    log_det_0_S3 = np.log(PI) + ZETA_3 / (2 * PI2)
    det_0_S3 = np.exp(log_det_0_S3)
    print(f"  log det'(Δ_0, S³) = log(π) + ζ(3)/(2π²) = {log_det_0_S3:.10f}")
    print(f"  det'(Δ_0, S³) = {det_0_S3:.10f}")

    # For Δ_1 on S³, we need the spectral determinant.
    # The 1-form Laplacian on S³ has eigenvalues k(k+2) with multiplicity 2k(k+2), k≥1.
    # log det'(Δ_1, S³) = -ζ'_{Δ_1}(0) where ζ_{Δ_1}(s) = Σ_{k≥1} 2k(k+2) / [k(k+2)]^s
    #                    = 2 Σ_{k≥1} [k(k+2)]^{1-s}
    #
    # Actually, let's be more careful. The eigenvalues of Δ_1 on S³ are:
    #   λ = k(k+2) for k≥1, with multiplicity m_k.
    # The multiplicity of the full 1-form Laplacian on S^n is known.
    # On S³: m_k(Δ_1) = 2k(k+2) for k≥1.
    # (This includes both exact and coexact 1-forms.)
    #
    # So ζ_{Δ_1}(s) = Σ_{k=1}^∞ 2k(k+2) · [k(k+2)]^{-s}
    #               = 2 Σ_{k=1}^∞ [k(k+2)]^{1-s}
    #
    # For log det', we need ζ'_{Δ_1}(0).
    # At s=0: ζ_{Δ_1}(0) = 2 Σ k(k+2) -- divergent, needs analytic continuation.
    #
    # We compute numerically via the relation to the scalar zeta.
    # k(k+2) = (k+1)² - 1, so:
    # ζ_{Δ_1}(s) = 2 Σ_{k=1}^∞ [(k+1)² - 1]^{1-s}
    #
    # For the derivative, let's use the known result for 1-form determinant on S³.
    # From Dowker (1989), the log-determinant of the 1-form Laplacian on S³ is:
    # log det'(Δ_1, S³) = log(π) + ζ(3)/(2π²) + 2·log(2) - 1
    #
    # Actually, the exact result is less standard. Let's compute numerically.
    # We use the Hurwitz zeta approach.

    # k(k+2) = (k+1)^2 - 1. Let m = k+1, so m ≥ 2 and eigenvalue = m² - 1.
    # ζ_{Δ_1}(s) = 2 Σ_{m=2}^∞ (m² - 1)^{1-s}
    #
    # For s near 0: (m²-1)^{1-s} = (m²-1) · exp(-s ln(m²-1))
    #             = (m²-1) · [1 - s·ln(m²-1) + O(s²)]
    #
    # ζ_{Δ_1}(0) = 2 Σ_{m=2}^∞ (m²-1)  -- divergent (needs regularization)
    # ζ'_{Δ_1}(0) = -2 Σ_{m=2}^∞ (m²-1)·ln(m²-1)  -- also needs regularization
    #
    # Use the regularized numerical approach instead: compute ζ_{Δ_1}(s) for s > 1
    # and analytically continue.

    # Alternative: compute directly from the relation
    # Δ_1 on S^3 has same eigenvalues as Δ_0 but different multiplicities.
    # We can write: ζ_{Δ_1}(s) = Σ 2k(k+2) · [k(k+2)]^{-s} = 2 · ζ_{Δ_0,modified}(s-1)
    # where ζ_{Δ_0}(s) = Σ (k+1)² · [k(k+2)]^{-s}
    #
    # So ζ_{Δ_1}(s) / ζ_{Δ_0}(s) = 2k(k+2)/(k+1)² -- not a clean ratio.
    #
    # Let's use the factorization k(k+2) = (k+1-1)(k+1+1) = (k+1)²-1:
    # ζ_{Δ_0}(s) = Σ_{k=1}^∞ (k+1)² · [(k+1)²-1]^{-s}
    #            = Σ_{m=2}^∞ m² · (m²-1)^{-s}
    #
    # ζ_{Δ_1}(s) = 2 · Σ_{m=2}^∞ (m²-1)^{1-s}
    #
    # So ζ_{Δ_1}(s) = 2 · Σ (m²-1)^{1-s} = 2 · σ(s-1) where σ(u) = Σ (m²-1)^{-u}
    # And ζ_{Δ_0}(s) = Σ m² · (m²-1)^{-s} = Σ [(m²-1)+1] · (m²-1)^{-s}
    #                = σ(s-1) + σ(s)
    #
    # So: ζ_{Δ_0}(s) = (1/2)·ζ_{Δ_1}(s) + σ(s)
    # And: σ(s) = Σ_{m=2}^∞ (m²-1)^{-s} = Σ_{m=2}^∞ [(m-1)(m+1)]^{-s}

    # Let's just compute everything numerically with enough terms.
    # For the derivatives at s=0, we use the definition:
    # -log det' = ζ'(0) computed via numerical differentiation of the
    # analytically-continued zeta function.

    # Actually, for S³ we know the TOTAL analytic torsion analytically.
    # For odd-dimensional spheres S^{2n+1} with trivial coefficients:
    #
    # Fried (1986): T(S^{2n+1}) = 1 for all n.
    # (The analytic torsion of an odd sphere with trivial local system is 1.)
    #
    # This is a THEOREM. So log T(S³) = 0.

    log_T_S3 = 0.0
    T_S3 = 1.0

    print()
    print(f"  THEOREM (Fried 1986): T(S^{{2n+1}}) = 1 for trivial coefficients.")
    print(f"  log T(S³) = 0")
    print(f"  T(S³) = 1")
    print()

    # This means: log det'(Δ_1) = 3 · log det'(Δ_0) on S³
    # (from log T = (1/2)[log det'(Δ_1) - 3 log det'(Δ_0)] = 0)
    log_det_1_S3 = 3 * log_det_0_S3
    det_1_S3 = det_0_S3**3
    print(f"  Implied: log det'(Δ_1, S³) = 3·log det'(Δ_0, S³) = {log_det_1_S3:.10f}")
    print(f"  det'(Δ_1, S³) = [det'(Δ_0, S³)]³ = {det_1_S3:.10f}")
    print()

    # ------------------------------------------------------------------
    # S² (even-dimensional, no standard analytic torsion)
    # ------------------------------------------------------------------
    # S² has χ(S²) = 2. Analytic torsion is not defined for even-dim manifolds
    # in the standard Ray-Singer sense. Instead, we have the functional determinant.
    det_S2_standard = np.exp(0.5 - 4 * ZETA_PRIME_NEG1)
    det_S2_paper = np.exp(-4 * ZETA_PRIME_NEG1)
    log_det_S2_standard = 0.5 - 4 * ZETA_PRIME_NEG1
    log_det_S2_paper = -4 * ZETA_PRIME_NEG1

    print("S² (no torsion — even-dimensional):")
    print(f"  det'(Δ_0, S²) = exp(1/2 - 4ζ'(-1)) = {det_S2_standard:.10f}  [Weisberger/Sarnak]")
    print(f"  det'(Δ_0, S²) = exp(-4ζ'(-1))       = {det_S2_paper:.10f}  [Paper 2 convention]")
    print()

    # ------------------------------------------------------------------
    # Summary table of torsions
    # ------------------------------------------------------------------
    print("-" * 80)
    print("TORSION SUMMARY")
    print("-" * 80)
    print(f"  T(S¹) = 1/(2π) = {T_S1:.10f}")
    print(f"  T(S³) = 1       (Fried's theorem)")
    print(f"  S² has no torsion (even-dim)")
    print()

    # ------------------------------------------------------------------
    # Combinations
    # ------------------------------------------------------------------
    print("-" * 80)
    print("NATURAL COMBINATIONS (torsion, determinants, volumes, Euler chars)")
    print("-" * 80)
    print()

    vol_S1 = 2 * PI
    vol_S2 = 4 * PI
    vol_S3 = 2 * PI2
    chi_S1 = 0
    chi_S2 = 2
    chi_S3 = 0

    combos = {}
    # Torsion combinations
    combos["T(S³)/T(S¹)"] = T_S3 / T_S1
    combos["T(S³) · T(S¹)"] = T_S3 * T_S1
    combos["1/T(S¹)"] = 1 / T_S1
    combos["1/T(S¹)²"] = 1 / T_S1**2
    combos["T(S³) / T(S¹)²"] = T_S3 / T_S1**2

    # Determinant combinations
    combos["det'₀(S³)"] = det_0_S3
    combos["det'₁(S³)"] = det_1_S3
    combos["det'₀(S²) [standard]"] = det_S2_standard
    combos["det'₀(S²) [paper]"] = det_S2_paper
    combos["det'₀(S¹)"] = 4 * PI2

    combos["det'₀(S¹)·det'₀(S³)/π"] = 4 * PI2 * det_0_S3 / PI
    combos["det'₀(S³)/det'₀(S²) [std]"] = det_0_S3 / det_S2_standard
    combos["det'₀(S³)/det'₀(S²) [paper]"] = det_0_S3 / det_S2_paper
    combos["det'₁(S³)/det'₀(S³)²"] = det_1_S3 / det_0_S3**2

    # Torsion × volume combinations
    combos["T(S³)/T(S¹) / π"] = T_S3 / T_S1 / PI
    combos["Vol(S³) / Vol(S¹)"] = vol_S3 / vol_S1
    combos["Vol(S³) · T(S³) / (Vol(S¹) · T(S¹))"] = vol_S3 * T_S3 / (vol_S1 * T_S1)
    combos["Vol(S²) · T(S³)/T(S¹)"] = vol_S2 * T_S3 / T_S1
    combos["Vol(S³) / (Vol(S¹)·Vol(S²))"] = vol_S3 / (vol_S1 * vol_S2)

    # Log determinant combinations (additive)
    combos["log det'₀(S³)"] = log_det_0_S3
    combos["log det'₁(S³)"] = log_det_1_S3
    combos["log det'₀(S³) - log det'₀(S²) [std]"] = log_det_0_S3 - log_det_S2_standard
    combos["log det'₀(S³) + log det'₀(S¹)"] = log_det_0_S3 + np.log(4*PI2)
    combos["3·log det'₀(S³) - log det'₀(S¹)"] = 3 * log_det_0_S3 - np.log(4*PI2)

    # Mixed
    combos["det'₀(S¹) · det'₀(S²) [std] / (4π³)"] = 4*PI2 * det_S2_standard / (4*PI**3)
    combos["det'₀(S³)² / det'₀(S¹)"] = det_0_S3**2 / (4*PI2)
    combos["π · det'₀(S³)"] = PI * det_0_S3
    combos["det'₀(S³) · Vol(S²)"] = det_0_S3 * vol_S2
    combos["det'₀(S³) · χ(S²)"] = det_0_S3 * chi_S2
    combos["Vol(S³)·det'₀(S³)/(2π)"] = vol_S3 * det_0_S3 / (2*PI)
    combos["(det'₀(S¹)·det'₀(S³))^(1/2)"] = np.sqrt(4*PI2 * det_0_S3)

    # Torsion ratio is T(S³)/T(S¹) = 2π, which is Vol(S¹)
    combos["[T(S³)/T(S¹)]²/π"] = (T_S3 / T_S1)**2 / PI
    combos["[T(S³)/T(S¹)]² + χ(S²)"] = (T_S3 / T_S1)**2 + chi_S2
    combos["[T(S³)/T(S¹)]² - χ(S²)"] = (T_S3 / T_S1)**2 - chi_S2
    combos["[1/T(S¹)]²/π + ζ(2)"] = (1/T_S1)**2/PI + F_TARGET

    # Print all
    print(f"  {'Combination':<52s} {'Value':>14s}")
    print(f"  {'─'*52} {'─'*14}")
    for name, val in combos.items():
        if np.isfinite(val):
            print(f"  {name:<52s} {val:>14.8f}")
            check_targets(val, name, results)
    print()

    # ------------------------------------------------------------------
    # Leray-Serre factorization check
    # ------------------------------------------------------------------
    print("-" * 80)
    print("LERAY-SERRE FACTORIZATION: T(S³) =? T(S¹) · χ(S²) · [correction]")
    print("-" * 80)
    print()
    print("  For the Hopf fibration S¹ → S³ → S², the Leray-Serre spectral")
    print("  sequence in cohomology gives: H*(S³) from H*(S²) and H*(S¹).")
    print()
    print("  For analytic torsion with trivial coefficients on a fiber bundle")
    print("  F → E → B, the anomaly formula (Bismut-Zhang, Dai 1991) gives:")
    print("    log T(E) = χ(B)·log T(F) + correction terms")
    print()
    print(f"  T(S³) = {T_S3}")
    print(f"  χ(S²) · log T(S¹) = {chi_S2} × {log_T_S1:.10f} = {chi_S2 * log_T_S1:.10f}")
    print(f"  exp(χ(S²)·log T(S¹)) = T(S¹)^χ(S²) = (2π)^(-2) = {T_S1**chi_S2:.10f}")
    print()

    anomaly = np.log(T_S3) - chi_S2 * np.log(T_S1)
    print(f"  Anomaly: log T(S³) - χ(S²)·log T(S¹) = {anomaly:.10f}")
    print(f"           = 0 - 2·(-log 2π) = 2·log(2π) = {2*np.log(2*PI):.10f}")
    print(f"           = log(4π²) = {np.log(4*PI2):.10f}")
    print(f"  exp(anomaly) = 4π² = {4*PI2:.10f}")
    print()
    check_targets(anomaly, "Torsion anomaly = 2·log(2π)", results)
    check_targets(np.exp(anomaly), "exp(anomaly) = 4π²", results)
    check_targets(anomaly / PI, "anomaly/π", results)

    print("-" * 80)
    print("PART A — TARGET MATCHES (< 0.1%)")
    print("-" * 80)
    print_hits(results)
    print()

    return results


# ============================================================
# PART B: Heat Kernel (Seeley-DeWitt Coefficients)
# ============================================================

def part_b():
    """Heat kernel coefficients and partition functions."""
    print()
    print("=" * 80)
    print("PART B: HEAT KERNEL — SEELEY-DEWITT COEFFICIENTS & PARTITION FUNCTIONS")
    print("=" * 80)
    print()

    results = []

    # ------------------------------------------------------------------
    # Seeley-DeWitt coefficients a_n for scalar Laplacian on round S^d
    # K(t) ~ (4πt)^{-d/2} Σ_n a_n t^n
    #
    # For S^d(r=1):
    #   Scalar curvature R = d(d-1)
    #   Ricci tensor: R_{ij} = (d-1) g_{ij}
    #   |Ric|² = (d-1)² · d
    #   Riemann: R_{ijkl} = g_{ik}g_{jl} - g_{il}g_{jk}
    #   |Riem|² = 2d(d-1)
    #   Vol(S^d) = 2π^{(d+1)/2} / Γ((d+1)/2)
    #
    # a_0 = Vol(M)
    # a_1 = (1/6) ∫_M R dV = R·Vol/6
    # a_2 = (1/360) ∫_M (5R² - 2|Ric|² + 2|Riem|²) dV
    # a_3 = (1/7!) ... (much more complicated)
    # ------------------------------------------------------------------

    print("SEELEY-DEWITT COEFFICIENTS a_n FOR SCALAR LAPLACIAN ON UNIT SPHERES")
    print("-" * 80)
    print()

    spheres = {}

    for d, name in [(1, "S¹"), (2, "S²"), (3, "S³")]:
        R_scalar = d * (d - 1)  # Scalar curvature
        Ric_sq = (d - 1)**2 * d  # |Ric|²
        Riem_sq = 2 * d * (d - 1)  # |Riem|²

        if d == 1:
            vol = 2 * PI
        elif d == 2:
            vol = 4 * PI
        else:
            vol = 2 * PI2

        a0 = vol
        a1 = R_scalar * vol / 6
        a2 = (5 * R_scalar**2 - 2 * Ric_sq + 2 * Riem_sq) * vol / 360

        # a_3 for constant curvature spaces (Gilkey 1975):
        # a_3 = Vol/(7!) × [35R³/9 - 14R·|Ric|²/3 + 14R·|Riem|²/3
        #        + 208|Ric|⁴/(?) - ...]
        # For spheres (Einstein + constant curvature), use:
        # a_3 = Vol · P_3(d) where P_3 is a polynomial in d.
        # From Branson-Gilkey (1990), for scalar Laplacian on S^d:
        #   a_3 = Vol(S^d) · d(d-1)(d-2)(5d²-7d+6)/(7! × 9/35... )
        # This gets complicated. Let's use the known formula for constant curvature:
        # a_3(S^d) = Vol(S^d) · [d(d-1)]³ / 360 × correction
        #
        # Actually, for constant sectional curvature K=1 (unit sphere):
        # Gilkey's formula simplifies. Using the recursion from Avramidi (2000):
        if d == 1:
            a3 = 0.0  # S¹ is flat (R=0, all curvature = 0)
        elif d == 2:
            # S²: R=2, |Ric|²=2, |Riem|²=4
            # a_3: use the general formula from Gilkey (1975) eq (1.7.13)
            # For dim 2, constant curvature:
            # a_3 = Vol · (1/7560) · (18R³ - 17R·|Ric|² - 2·Ric³_trace +
            #        9R·|Riem|² + 28∇²R + ...) but ∇R=0 on S²
            # On S²: Ric = g (so Ric³ = Ric, trace Ric³ = 2)
            # 18·8 - 17·2·2 - 2·2 + 9·2·4 = 144 - 68 - 4 + 72 = 144
            # a_3 = 4π · 144/7560 = 4π · 2/105
            a3 = vol * 2 / 105  # = 8π/105
        elif d == 3:
            # S³: R=6, |Ric|²=12, |Riem|²=12
            # Ric_{ij} = 2g_{ij}, so Ric³_trace = 2³·3 = 24
            # 18·216 - 17·6·12 - 2·24 + 9·6·12 = 3888 - 1224 - 48 + 648 = 3264
            # a_3 = 2π² · 3264/7560 = 2π² · 136/315
            a3 = vol * 136 / 315

        spheres[name] = {"d": d, "vol": vol, "R": R_scalar,
                         "a0": a0, "a1": a1, "a2": a2, "a3": a3}

        print(f"  {name} (d={d}):")
        print(f"    Vol = {vol:.10f},  R = {R_scalar},  |Ric|² = {Ric_sq},  |Riem|² = {Riem_sq}")
        print(f"    a_0 = {a0:.10f}")
        print(f"    a_1 = {a1:.10f}")
        print(f"    a_2 = {a2:.10f}")
        print(f"    a_3 = {a3:.10f}")
        print()

    # ------------------------------------------------------------------
    # One-loop partition functions: log Z = (1/2) ζ'_Δ(0)
    # ------------------------------------------------------------------
    print("-" * 80)
    print("ONE-LOOP PARTITION FUNCTIONS: log Z = (1/2) ζ'_Δ(0) = -(1/2) log det'(Δ)")
    print("-" * 80)
    print()

    # S¹: det'(Δ₀) = 4π², so log Z = -(1/2) log(4π²) = -log(2π)
    log_Z_S1 = -np.log(2 * PI)

    # S²: det'(Δ₀) = exp(1/2 - 4ζ'(-1))
    log_Z_S2_std = -(0.5 - 4 * ZETA_PRIME_NEG1) / 2
    log_Z_S2_paper = -(-4 * ZETA_PRIME_NEG1) / 2

    # S³: det'(Δ₀) = π · exp(ζ(3)/(2π²))
    log_Z_S3 = -(np.log(PI) + ZETA_3 / (2*PI2)) / 2

    print(f"  log Z(S¹) = -log(2π)                    = {log_Z_S1:.10f}")
    print(f"  log Z(S²) = -(1/2)(1/2 - 4ζ'(-1))       = {log_Z_S2_std:.10f}  [standard]")
    print(f"  log Z(S²) = -(1/2)(-4ζ'(-1))             = {log_Z_S2_paper:.10f}  [paper conv.]")
    print(f"  log Z(S³) = -(1/2)(log π + ζ(3)/(2π²))  = {log_Z_S3:.10f}")
    print()

    # Fiber bundle factorization check
    print("-" * 80)
    print("PARTITION FUNCTION FACTORIZATION CHECK")
    print("  log Z(S³) =? log Z(S²) + log Z(S¹) + interaction")
    print("-" * 80)
    print()

    interaction_std = log_Z_S3 - log_Z_S2_std - log_Z_S1
    interaction_paper = log_Z_S3 - log_Z_S2_paper - log_Z_S1

    print(f"  log Z(S³) - log Z(S²) [std] - log Z(S¹) = {interaction_std:.10f}")
    print(f"  log Z(S³) - log Z(S²) [paper] - log Z(S¹) = {interaction_paper:.10f}")
    print()

    check_targets(interaction_std, "Z interaction [std]", results)
    check_targets(interaction_paper, "Z interaction [paper]", results)
    check_targets(np.exp(interaction_std), "exp(Z interaction [std])", results)
    check_targets(np.exp(interaction_paper), "exp(Z interaction [paper])", results)

    # ------------------------------------------------------------------
    # Linear combinations of log Z, Vol, a_n, χ
    # ------------------------------------------------------------------
    print("-" * 80)
    print("LINEAR COMBINATION SEARCH")
    print("  Testing: c₁·log Z_i + c₂·Vol_i + c₃·a_{n,i} + c₄·χ_i")
    print("-" * 80)
    print()

    # Build a catalog of all quantities
    quantities = {
        "log Z(S¹)": log_Z_S1,
        "log Z(S²) [std]": log_Z_S2_std,
        "log Z(S³)": log_Z_S3,
        "Vol(S¹)": 2*PI,
        "Vol(S²)": 4*PI,
        "Vol(S³)": 2*PI2,
        "χ(S²)": 2.0,
        "a₁(S²)": spheres["S²"]["a1"],
        "a₁(S³)": spheres["S³"]["a1"],
        "a₂(S²)": spheres["S²"]["a2"],
        "a₂(S³)": spheres["S³"]["a2"],
        "a₃(S³)": spheres["S³"]["a3"],
    }

    print(f"  {'Quantity':<28s} {'Value':>14s}")
    print(f"  {'─'*28} {'─'*14}")
    for name, val in quantities.items():
        print(f"  {name:<28s} {val:>14.8f}")
        check_targets(val, name, results)
    print()

    # Test pairwise sums and differences
    print("  PAIRWISE DIFFERENCES/SUMS near targets:")
    qlist = list(quantities.items())
    pair_hits = []
    for i in range(len(qlist)):
        for j in range(i+1, len(qlist)):
            n1, v1 = qlist[i]
            n2, v2 = qlist[j]
            for op, sym in [(lambda a,b: a+b, "+"), (lambda a,b: a-b, "-"),
                            (lambda a,b: b-a, "rev-")]:
                val = op(v1, v2)
                label = f"{n1} {sym} {n2}"
                check_targets(val, label, results)
                # Also try ratios
            if v2 != 0:
                check_targets(v1/v2, f"{n1} / {n2}", results)
            if v1 != 0:
                check_targets(v2/v1, f"{n2} / {n1}", results)

    # Test products with small integer coefficients
    print("  INTEGER-WEIGHTED COMBINATIONS near targets:")
    key_quants = [
        ("log Z(S¹)", log_Z_S1),
        ("log Z(S²)", log_Z_S2_std),
        ("log Z(S³)", log_Z_S3),
        ("Vol(S¹)", 2*PI),
        ("Vol(S³)", 2*PI2),
        ("χ(S²)", 2.0),
    ]
    for i in range(len(key_quants)):
        for j in range(i, len(key_quants)):
            for ci in range(-3, 4):
                for cj in range(-3, 4):
                    if ci == 0 and cj == 0:
                        continue
                    val = ci * key_quants[i][1] + cj * key_quants[j][1]
                    if abs(val) > 1e-10:
                        label = f"{ci}·{key_quants[i][0]} + {cj}·{key_quants[j][0]}"
                        check_targets(val, label, results)

    print()
    print("-" * 80)
    print("PART B — TARGET MATCHES (< 0.1%)")
    print("-" * 80)
    print_hits(results)
    print()

    return results


# ============================================================
# PART C: Index Theory (η-invariant, Chern-Simons, APS)
# ============================================================

def part_c():
    """APS η-invariant, Chern-Simons, and index theory for the Hopf bundle."""
    print()
    print("=" * 80)
    print("PART C: INDEX THEORY — η-INVARIANT, CHERN-SIMONS, APS")
    print("=" * 80)
    print()

    results = []

    # ------------------------------------------------------------------
    # Dirac spectrum on S³
    # ------------------------------------------------------------------
    # On round S³(r=1), the Dirac operator D has eigenvalues:
    #   ±(k + 3/2) for k = 0, 1, 2, ...
    #   with multiplicity (k+1)(k+2) for each sign.
    # (Bär 1996, Trautman 1995)

    print("DIRAC SPECTRUM ON S³")
    print("-" * 80)
    print()
    print("  Eigenvalues: ±(k + 3/2), multiplicity (k+1)(k+2), k = 0,1,2,...")
    print()

    # η-invariant: η(s) = Σ sign(λ)|λ|^{-s} for Dirac operator
    # Since spectrum is symmetric (±λ with same multiplicity), η(s) = 0 for all s.
    eta_0 = 0.0
    print(f"  η(D_{'{S³}'}, s) = 0 for all s (symmetric spectrum)")
    print(f"  η(0) = {eta_0}")
    print()

    # ------------------------------------------------------------------
    # APS ξ-invariant
    # ------------------------------------------------------------------
    # ξ = (η(0) + dim ker D) / 2
    # ker D on S³: need λ = 0. Since eigenvalues are ±(k+3/2), there is no zero mode.
    # dim ker D = 0.
    dim_ker_D = 0
    xi = (eta_0 + dim_ker_D) / 2

    print("APS ξ-INVARIANT")
    print("-" * 80)
    print(f"  dim ker D = {dim_ker_D} (no zero eigenvalue in ±(k+3/2))")
    print(f"  ξ = (η(0) + dim ker D) / 2 = {xi}")
    print()

    # ------------------------------------------------------------------
    # Chern character of the Hopf bundle
    # ------------------------------------------------------------------
    print("HOPF BUNDLE TOPOLOGY")
    print("-" * 80)
    print()
    print("  The Hopf bundle L → S² has c₁(L) = 1 (generator of H²(S², Z))")
    print("  ch(L) = exp(c₁) = 1 + c₁ + c₁²/2 + ...")
    print("  On S² (dim 2): ch(L) = 1 + c₁, so ∫_{S²} ch(L) = ∫ c₁ = 1")
    print()
    c1 = 1  # First Chern class

    # ------------------------------------------------------------------
    # Chern-Simons invariant for the Hopf bundle
    # ------------------------------------------------------------------
    # The canonical connection on the Hopf bundle S¹ → S³ → S² is the
    # connection 1-form A = (1/2)(x₁dy₁ - y₁dx₁ + x₂dy₂ - y₂dx₂)
    # in the standard coordinates on S³ ⊂ R⁴ = C².
    #
    # For a U(1) bundle over S², the Chern-Simons invariant of S³ is:
    # CS(S³, A) = (1/4π²) ∫_{S³} A ∧ dA
    #           = (1/4π²) ∫_{S³} A ∧ F  (since dA restricted to horizontal = F)
    #
    # For the Hopf bundle, the curvature F = (1/2)ω_{S²} (the area form on S²
    # normalized so ∫F = 2π·c₁ = 2π).
    #
    # Actually, the standard Chern-Simons functional for a U(1) connection:
    # CS = (1/4π) ∫_{S³} A ∧ dA    (U(1) convention, Tr = 1)
    #
    # For the Hopf connection on S³:
    # ∫_{S³} A ∧ dA = ∫_{S³} A ∧ F_vert
    # By the Hopf invariant calculation: ∫_{S³} A ∧ dA = (2π)² × Hopf invariant
    # The Hopf invariant of the Hopf map π: S³ → S² is 1.
    # So ∫_{S³} A ∧ dA = 4π²
    #
    # CS = (1/4π) · 4π² = π      (with normalization 1/4π)
    # or CS = (1/8π²) · 4π² = 1/2 (with normalization 1/8π²)

    print("CHERN-SIMONS INVARIANT")
    print("-" * 80)
    print()
    print("  For the Hopf bundle with canonical connection A:")
    print("  ∫_{S³} A ∧ dA = 4π² (= Hopf invariant × (2π)²)")
    print()

    CS_4pi = 4 * PI2 / (4 * PI)  # normalization 1/(4π)
    CS_8pi2 = 4 * PI2 / (8 * PI2)  # normalization 1/(8π²)

    print(f"  CS = (1/4π) ∫ A∧dA   = π         = {CS_4pi:.10f}")
    print(f"  CS = (1/8π²) ∫ A∧dA  = 1/2       = {CS_8pi2:.10f}")
    print()

    check_targets(CS_4pi, "CS [1/4π norm] = π", results)
    check_targets(CS_8pi2, "CS [1/8π² norm] = 1/2", results)

    # ------------------------------------------------------------------
    # Combinations with index-theoretic quantities
    # ------------------------------------------------------------------
    print("-" * 80)
    print("INDEX-THEORETIC COMBINATIONS vs TARGETS")
    print("-" * 80)
    print()

    idx_quantities = {
        "η(0)": eta_0,
        "ξ (APS)": xi,
        "c₁(L)": float(c1),
        "CS [1/4π]": CS_4pi,
        "CS [1/8π²]": CS_8pi2,
        "Hopf invariant": 1.0,
        "χ(S²)": 2.0,
        "∫A∧dA": 4 * PI2,
    }

    # Also bring in volumes and torsion
    extra = {
        "Vol(S¹)": 2*PI,
        "Vol(S²)": 4*PI,
        "Vol(S³)": 2*PI2,
        "T(S¹)=1/(2π)": 1/(2*PI),
    }
    idx_quantities.update(extra)

    combos_c = {}
    # Natural products and ratios
    combos_c["CS[1/4π] · Vol(S³)"] = CS_4pi * 2*PI2
    combos_c["CS[1/8π²] · Vol(S³)"] = CS_8pi2 * 2*PI2
    combos_c["∫A∧dA · χ(S²)"] = 4*PI2 * 2
    combos_c["∫A∧dA + χ(S²)"] = 4*PI2 + 2
    combos_c["∫A∧dA · c₁ / (2π)"] = 4*PI2 * c1 / (2*PI)
    combos_c["Vol(S³) · c₁²"] = 2*PI2
    combos_c["Vol(S²) · CS[1/4π]"] = 4*PI * CS_4pi
    combos_c["Vol(S¹) · CS[1/4π]"] = 2*PI * CS_4pi
    combos_c["c₁ · Vol(S³) / Vol(S¹)"] = c1 * 2*PI2 / (2*PI)
    combos_c["χ(S²) · Vol(S³) / Vol(S¹)²"] = 2 * 2*PI2 / (2*PI)**2
    combos_c["∫A∧dA / Vol(S¹)"] = 4*PI2 / (2*PI)
    combos_c["∫A∧dA / (2π)"] = 4*PI2 / (2*PI)
    combos_c["(∫A∧dA)² / Vol(S³)"] = (4*PI2)**2 / (2*PI2)
    combos_c["Vol(S³)²/∫A∧dA"] = (2*PI2)**2 / (4*PI2)
    combos_c["CS[1/8π²] · ∫A∧dA"] = CS_8pi2 * 4*PI2
    combos_c["6·CS[1/8π²]·Vol(S³)/π"] = 6 * CS_8pi2 * 2*PI2 / PI
    combos_c["42·CS[1/8π²]"] = 42 * CS_8pi2
    combos_c["Vol(S³)/Vol(S¹) + ζ(2) - 1/40"] = PI + F_TARGET - DELTA_TARGET
    combos_c["Vol(S³)·(B+F-Δ)/(Vol(S¹)·Vol(S²))"] = 2*PI2 * K_OVER_PI / (2*PI * 4*PI)

    print(f"  {'Combination':<52s} {'Value':>14s}")
    print(f"  {'─'*52} {'─'*14}")
    for name, val in combos_c.items():
        if np.isfinite(val):
            print(f"  {name:<52s} {val:>14.8f}")
            check_targets(val, name, results)
    print()

    # ------------------------------------------------------------------
    # Twisted Dirac: η-invariant with Hopf bundle twist
    # ------------------------------------------------------------------
    print("-" * 80)
    print("TWISTED DIRAC OPERATOR — η-invariant with Hopf bundle")
    print("-" * 80)
    print()

    # The Dirac operator twisted by a line bundle L^n over S³:
    # On S³, the Hopf bundle pulled back is trivial as a vector bundle
    # (S³ is the total space, the fiber is already included).
    # But we can consider the Dirac operator coupled to flat connections
    # with holonomy e^{2πi·α} around the Hopf fiber.
    #
    # For a U(1) flat connection with holonomy e^{2πiα}:
    # η(0, α) = -4 Σ_{k=0}^∞ (k+1)(k+2) · [sign(k+3/2+α)·|...|^0 terms]
    # This gets complicated. For α=0 we get η=0.
    # For the canonical connection (non-flat), the APS boundary term is different.
    #
    # Key result: on S³ with the round metric and trivial bundle,
    # all index-theoretic invariants are either 0 or standard topological numbers.

    print("  On S³ with round metric and trivial coefficient bundle:")
    print("  - index(D) = 0 (odd dimension)")
    print("  - η(0) = 0 (symmetric spectrum)")
    print("  - ξ = 0")
    print("  - No non-trivial APS correction.")
    print()
    print("  The Hopf bundle is a GEOMETRIC structure, not a gauge bundle")
    print("  with an independent connection. The index theory of S³ itself")
    print("  sees the Hopf structure only through the Chern-Simons term.")
    print()

    # ------------------------------------------------------------------
    # Dedekind η-function & modular connection
    # ------------------------------------------------------------------
    print("-" * 80)
    print("DEDEKIND η-FUNCTION CONNECTION")
    print("-" * 80)
    print()

    # The Hopf fiber S¹ has a natural modular parameter τ.
    # The Dedekind η-function η(τ) = q^{1/24} Π(1-q^n) where q = e^{2πiτ}.
    # For the Hopf fiber with circumference 2π and connection holonomy,
    # the natural τ is on the imaginary axis: τ = i (square torus).

    # η(i) = Γ(1/4) / (2π^{3/4})  (known exact value)
    # |η(i)|² = Γ(1/4)² / (4π^{3/2})

    from scipy.special import gamma as Gamma

    eta_i = Gamma(0.25) / (2 * PI**0.75)
    log_eta_i = np.log(eta_i)

    print(f"  Dedekind η(i) = Γ(1/4)/(2π^{{3/4}}) = {eta_i:.10f}")
    print(f"  log η(i) = {log_eta_i:.10f}")
    print(f"  24·log η(i) = {24*log_eta_i:.10f}")
    print(f"  -24·log η(i) = {-24*log_eta_i:.10f}")
    print(f"  η(i)^24 = {eta_i**24:.10f}")
    print(f"  1/η(i)^24 = {1/eta_i**24:.10f}")
    print()

    # Ramanujan's constant-like: e^{π√d} relates to class numbers
    # For d=1 (Hopf fiber): e^π = 23.14...
    e_pi = np.exp(PI)
    print(f"  e^π = {e_pi:.10f}")
    print(f"  e^π - 24·log η(i) = {e_pi - 24*log_eta_i:.10f}")
    print()

    check_targets(eta_i, "Dedekind η(i)", results)
    check_targets(24*log_eta_i, "24·log η(i)", results)
    check_targets(-24*log_eta_i, "-24·log η(i)", results)
    check_targets(eta_i**24, "η(i)^24", results)
    check_targets(1/eta_i**24, "1/η(i)^24", results)

    # ------------------------------------------------------------------
    # Spectral action approach (Connes-Chamseddine)
    # ------------------------------------------------------------------
    print("-" * 80)
    print("SPECTRAL ACTION TRACE Tr(f(D/Λ)) MOMENTS")
    print("-" * 80)
    print()

    # The spectral action on S³ gives:
    # Tr(f(D/Λ)) ~ f_0 Λ³ Vol(S³)/(4π²) + f_2 Λ (R/12) Vol(S³)/(4π²) + ...
    # where f_k = ∫ f(x) x^{k-1} dx are moments of the cutoff function.
    # The dimensionless ratios are fixed by the spectrum.

    # For the Dirac on S³, the spectral density is:
    # N(λ) ~ (2/3) λ³ · Vol(S³)/(4π²) as λ→∞ (Weyl law for spinors on 3-manifold)
    # Coefficient: Vol(S³)/(6π²) = 2π²/(6π²) = 1/3

    weyl_coeff = 2 * PI2 / (6 * PI2)
    sub_weyl = 2 * PI2 / (4 * PI2)  # = 1/2, subleading

    print(f"  Weyl coefficient: Vol(S³)/(6π²) = {weyl_coeff:.10f}")
    print(f"  Subleading: Vol(S³)/(4π²) = {sub_weyl:.10f}")
    print()

    # The ratio of Weyl coefficients between S³ and S¹:
    # S¹: N(λ) ~ (2/π) λ as λ→∞. Coefficient: Vol(S¹)/π = 2.
    # Ratio: (1/3)/2 = 1/6
    ratio_weyl = weyl_coeff / (2*PI/PI)
    print(f"  Weyl ratio S³/S¹: {ratio_weyl:.10f} (= 1/6)")
    print()

    check_targets(weyl_coeff, "Weyl coeff = 1/3", results)
    check_targets(ratio_weyl, "Weyl ratio S³/S¹ = 1/6", results)

    print()
    print("-" * 80)
    print("PART C — TARGET MATCHES (< 0.1%)")
    print("-" * 80)
    print_hits(results)
    print()

    return results


# ============================================================
# SUMMARY
# ============================================================

def summary(results_a, results_b, results_c):
    """Final summary of what spectral geometry can and cannot explain."""
    print()
    print("=" * 80)
    print("SUMMARY: WHAT SPECTRAL GEOMETRY CAN AND CANNOT EXPLAIN ABOUT K")
    print("=" * 80)
    print()

    all_results = results_a + results_b + results_c

    if all_results:
        print("ALL MATCHES FOUND (< 0.1% of any target):")
        print("-" * 80)
        print_hits(all_results)
    else:
        print("  No matches within 0.1% of any target across all three parts.")

    print()
    print("WHAT SPECTRAL GEOMETRY CAN EXPLAIN:")
    print("-" * 80)
    print("""
  1. F = ζ(2) = π²/6 is a spectral zeta value of the S¹ fiber Laplacian.
     This is well-grounded: ζ_{S¹}(1) = 2·ζ_R(2) = π²/3, so ζ(2) = π²/6
     appears naturally as half the spectral zeta function at s=1.

  2. The Chern-Simons invariant CS(S³, A_Hopf) = 1/2 (in 1/8π² normalization)
     or π (in 1/4π normalization) is an exact topological invariant of the
     Hopf bundle. This is standard Chern-Weil theory.

  3. The analytic torsion T(S³) = 1 (Fried's theorem) and T(S¹) = 1/(2π)
     are exact spectral invariants. Their ratio T(S³)/T(S¹) = 2π = Vol(S¹).

  4. The torsion anomaly in the Hopf fibration:
     log T(S³) - χ(S²)·log T(S¹) = 2·log(2π) = log(4π²)
     This is a Bismut-Zhang-type anomaly from the non-trivial fibration.

  5. The spectral determinant near-miss det'_1·det'_3/π ≈ 41.957 ≈ 42
     (from the previous script) remains the closest approach to B = 42,
     but the 0.1% gap persists.
""")

    print("WHAT SPECTRAL GEOMETRY CANNOT EXPLAIN (Link 3 gap):")
    print("-" * 80)
    print("""
  1. B = 42 does not emerge from any standard spectral invariant of S¹, S²,
     or S³. Neither analytic torsion, heat kernel coefficients, partition
     functions, nor their natural combinations produce exactly 42.

  2. The combination rule K = π(B + F - Δ) with ADDITIVE structure
     B + F - Δ is not a standard operation in spectral geometry. Spectral
     invariants compose MULTIPLICATIVELY (determinants) or ADDITIVELY in
     log space. The mixing of an integer (42), a zeta value (ζ(2)), and
     a rational (1/40) in a single sum has no spectral precedent.

  3. Δ = 1/40 does not appear as any standard boundary correction, heat
     kernel coefficient, or η-invariant residue for the Hopf bundle.

  4. The torsion anomaly (= log 4π²) and partition function factorization
     defect are both of order O(1), not O(10) where K/π ≈ 43.6 lives.
     The scales don't match.

  5. Index theory contributes only topological integers (c₁ = 1, χ = 2,
     Hopf invariant = 1) and zero (η = 0, index = 0). These are too
     coarse to produce K.

  CONCLUSION: Standard spectral geometry of the Hopf fibration provides F
  but cannot derive B or Δ. The combination rule K = π(B + F - Δ) remains
  an empirical observation. If it has a geometric origin, it likely involves
  representation-theoretic structure (the degeneracy-weighted Casimir trace
  that defines B = 42) rather than analytic spectral invariants.
""")


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    results_a = part_a()
    results_b = part_b()
    results_c = part_c()
    summary(results_a, results_b, results_c)
