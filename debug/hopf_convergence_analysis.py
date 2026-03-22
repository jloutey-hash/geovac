"""
Hopf Bundle Convergence Analysis
=================================
Debug script analyzing the convergence of discrete lattice volumes
toward continuum values, and searching for algebraic structures
that produce 1/α = 137.036.

Date: March 2026
Status: Exploratory — Paper 2 conjecture tier

Key findings:
    - Discrete S³ volume converges ~4% low (midpoint rule error)
    - Discrete S² area collapses to 0 (clustered base points)
    - No candidate formula from discrete volumes converges to 1/α
    - The linear (l,m) → (ψ₁,ψ₂) mapping is not representation-theoretic
    - Exact algebraic identity: Σ(2l+1)·W(l)/2π = Σ 2l = N_states - N_shells
    - Fiber coverage per state converges to 2π = Vol(S¹) from below
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.hopf_bundle import (
    convergence_study, alpha_formula_terms, algebraic_alpha_analysis,
    hopf_fiber_analysis, fiber_winding_exact, fock_chi,
    discrete_s3_volume, discrete_fiber_volume, discrete_base_area,
    fock_embedding, fock_hopf_decompose, cg_expectation_mapping,
    ALPHA_INV_EXPERIMENT, VOL_S1, VOL_S2, VOL_S3
)


def print_header(title: str) -> None:
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def run_alpha_formula_breakdown() -> None:
    """Break down 1/α ≈ 4π³ + π² + π into sphere volume terms."""
    print_header("1. ALPHA FORMULA: 4π³ + π² + π = 137.036")
    terms = alpha_formula_terms()
    for k, v in terms.items():
        print(f"  {k}: {v}")


def run_convergence_study() -> None:
    """Test convergence of discrete volumes toward continuum."""
    print_header("2. CONVERGENCE STUDY: Discrete volumes vs continuum")
    results = convergence_study([2, 5, 10, 20, 50, 100])
    for i, nmax in enumerate(results['n_max_values']):
        v3 = results['vol_s3'][i]
        v1 = results['vol_s1'][i]
        v2 = results['vol_s2'][i]
        print(f"\n  n_max = {nmax}:")
        print(f"    Vol(S³): {v3:.6f}  (target {VOL_S3:.6f}, "
              f"err {abs(v3 - VOL_S3) / VOL_S3 * 100:.2f}%)")
        print(f"    Vol(S¹): {v1:.6f}  (target {VOL_S1:.6f}, "
              f"err {abs(v1 - VOL_S1) / VOL_S1 * 100:.2f}%)")
        print(f"    Vol(S²): {v2:.6f}  (target {VOL_S2:.6f}, "
              f"err {abs(v2 - VOL_S2) / VOL_S2 * 100:.2f}%)")
        print(f"    Candidates:")
        for name, val in results['alpha_candidates'][i].items():
            err = abs(val - ALPHA_INV_EXPERIMENT) / ALPHA_INV_EXPERIMENT * 100
            marker = " <<<" if err < 1.0 else ""
            print(f"      {name}: {val:.6f}  (err {err:.4f}%){marker}")


def run_algebraic_analysis() -> None:
    """Algebraic fiber structure analysis."""
    print_header("3. ALGEBRAIC ANALYSIS: Fiber fractions and Casimir sums")
    for nmax in [5, 10, 20, 50, 100]:
        alg = algebraic_alpha_analysis(nmax)
        print(f"\n  n_max = {nmax}:")
        print(f"    Total states: {alg['total_states']}")
        print(f"    Avg fiber fraction: {alg['avg_fiber_fraction']:.6f}")
        print(f"    Effective Vol(S¹): {alg['effective_vol_s1']:.6f}")
        print(f"    Discrete Vol(S³): {alg['discrete_vol_s3']:.6f}")
        print(f"    Candidate A: {alg['candidate_A_fiber_x_total_plus_correction']:.6f} "
              f"(err {alg['candidate_A_error'] * 100:.4f}%)")


def run_fiber_winding_sums() -> None:
    """Degeneracy-weighted fiber winding sums."""
    print_header("4. FIBER WINDING SUMS: W(l) = 4πl/(2l+1)")
    print("\n  Exact W(l) values:")
    for l in range(1, 11):
        W = fiber_winding_exact(l)
        print(f"    l={l:2d}: W={W:.6f}, W/2π = 2·{l}/(2·{l}+1) = {2*l/(2*l+1):.6f}")

    print("\n  Degeneracy-weighted sums (exact):")
    print(f"  {'n_max':>5s} | {'Σ 2l':>10s} | {'N_states':>10s} | "
          f"{'Σ(2l+1)W(l)':>12s} | {'avg/state':>10s}")
    print("  " + "-" * 60)
    for nmax in [5, 10, 20, 50, 100, 137]:
        sum_2l = sum(2 * l for n in range(1, nmax + 1) for l in range(1, n))
        N_tot = nmax * (nmax + 1) * (2 * nmax + 1) // 6
        sum_wl = sum(
            (2 * l + 1) * fiber_winding_exact(l)
            for n in range(1, nmax + 1) for l in range(1, n)
        )
        avg = sum_wl / N_tot if N_tot > 0 else 0.0
        print(f"  {nmax:5d} | {sum_2l:10d} | {N_tot:10d} | "
              f"{sum_wl:12.2f} | {avg:10.6f}")

    # Key identity: Σ(2l+1)·W(l)/2π = Σ 2l exactly (integer)
    print("\n  EXACT IDENTITY: Σ_{n,l} (2l+1) × W(l)/(2π) = Σ_{n,l} 2l = N_states - N_shells")
    print("  This is a tautology: (2l+1) × 2l/(2l+1) = 2l.")
    print("  Average fiber coverage per state → 2π = Vol(S¹) as n_max → ∞.")


def run_fock_embedding_comparison() -> None:
    """Compare linear, Fock-embedding, and CG-based S³ mappings."""
    print_header("5. MAPPING COMPARISON: Linear vs Fock-embedding vs CG")
    print("\n  For each |n,l,m⟩, compare three S³ mapping methods:")
    print("    Linear: ψ₁ = πl/(n-1), ψ₂ = 2πm/(2l+1)")
    print("    Fock: (ξ₀,ξ₁,ξ₂,ξ₃) from semiclassical (χ,θ,φ)")
    print("    CG: ⟨m⁺⟩, ⟨m⁻⟩ from SO(4) Clebsch-Gordan decomposition")

    for n in range(2, 5):
        print(f"\n  --- n={n} ---")
        for l in range(n):
            for m in range(-l, l + 1):
                # Fock embedding
                decomp = fock_hopf_decompose(n, l, m)
                fock_fiber = decomp['fiber']
                fock_base = decomp['base_phase']
                fock_chi_h = decomp['chi_hopf']

                # CG expectation
                cg = cg_expectation_mapping(n, l, m)
                cg_psi1 = cg['psi1']
                cg_psi2 = cg['psi2']
                cg_nterms = cg['n_cg_terms']

                print(f"    |{n},{l},{m:+d}>: "
                      f"Fock(fiber={fock_fiber:+.3f}, base={fock_base:+.3f}, "
                      f"χ_H={fock_chi_h:.3f}) | "
                      f"CG(ψ₁={cg_psi1:.3f}, ψ₂={cg_psi2:.3f}, "
                      f"terms={cg_nterms})")


def run_fock_fiber_winding() -> None:
    """Check whether Fock-embedding fiber angle tracks m linearly."""
    print_header("6. FOCK-EMBEDDING FIBER WINDING ACROSS m-LADDERS")
    print("\n  Does the Fock-embedding fiber angle increase linearly with m?")

    for n in [3, 4, 5]:
        for l in range(1, n):
            fibers = []
            for m in range(-l, l + 1):
                decomp = fock_hopf_decompose(n, l, m)
                fibers.append((m, decomp['fiber']))

            winding = fibers[-1][1] - fibers[0][1]
            expected = fiber_winding_exact(l)
            steps = [(fibers[i + 1][1] - fibers[i][1])
                     for i in range(len(fibers) - 1)]

            is_monotonic = all(s > 0 for s in steps) or all(s < 0 for s in steps)
            print(f"  n={n}, l={l}: winding={winding:+.4f} "
                  f"(exact W(l)={expected:.4f}), "
                  f"monotonic={'YES' if is_monotonic else 'NO'}, "
                  f"steps={[f'{s:+.3f}' for s in steps]}")


def run_alpha_search_summary() -> None:
    """Summary of all approaches to deriving 1/α from lattice structure."""
    print_header("7. SUMMARY: APPROACHES TO 1/α FROM LATTICE STRUCTURE")
    target = ALPHA_INV_EXPERIMENT

    print(f"\n  Target: 1/α = {target}")
    print(f"  Formula: 4π³ + π² + π = {4*np.pi**3 + np.pi**2 + np.pi:.6f} "
          f"(err {abs(4*np.pi**3 + np.pi**2 + np.pi - target)/target*100:.4f}%)")

    print("\n  Status of approaches:")
    print("  1. Discrete volume products    → DOES NOT CONVERGE")
    print("     Vol(S³) is 4% low, Vol(S²) collapses to 0.")
    print("     The centroid-deficit method for S² area fails for clustered points.")
    print("  2. Algebraic fiber × total     → SLOW CONVERGENCE (5.5% at nmax=100)")
    print("     Candidate A = v1_eff × v3_disc + correction.")
    print("     Limited by both Vol(S³) quadrature and fiber fraction.")
    print("  3. Exact counting identities   → TAUTOLOGICAL")
    print("     Σ(2l+1)·W(l)/2π = Σ 2l is an algebraic identity,")
    print("     not a physical relationship.")
    print("  4. Ratio search                → NO MATCH")
    print("     No algebraic ratio of lattice sums converges to 137.036.")

    print("\n  DIAGNOSIS:")
    print("  The (l,m) → (ψ₁,ψ₂) mapping is approximate (linear vs")
    print("  representation-theoretic). The Fock embedding gives exact S³ points")
    print("  but the Hopf fiber angle does not track m linearly — the Hopf")
    print("  fibration and the angular momentum decomposition are different")
    print("  decompositions of S³.")
    print("")
    print("  The CG decomposition shows each |n,l,m⟩ is a SUPERPOSITION of")
    print("  (m⁺,m⁻) pairs. There is no unique (ψ₁,ψ₂) point assignment.")
    print("  The expectation values ⟨m⁺⟩ = m/2 give a well-defined but")
    print("  information-losing map.")
    print("")
    print("  To derive 1/α from the lattice, one likely needs:")
    print("  - Full Hopf holonomy (Berry phase), not just fiber angles")
    print("  - The connection (gauge field) on the Hopf bundle")
    print("  - Integration of the curvature 2-form over the base S²")
    print("  - This is the Chern number approach (Paper 2 conjecture)")


if __name__ == '__main__':
    run_alpha_formula_breakdown()
    run_convergence_study()
    run_algebraic_analysis()
    run_fiber_winding_sums()
    run_fock_embedding_comparison()
    run_fock_fiber_winding()
    run_alpha_search_summary()
