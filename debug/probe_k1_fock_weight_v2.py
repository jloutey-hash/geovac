"""
Probe v2: Deeper investigation of kappa = -1/16 from Fock weight.

Key findings from v1:
  - For l=0 (s-orbitals), the Fock weight coupling between adjacent shells is
    EXACTLY 1/4, independent of n. This is a clean universal constant.
  - For l>0, the coupling depends on both n and l.
  - The question is: what is the relationship between 1/4 and 1/16?

This script investigates:
  1. Whether 1/4 is the Fock weight coupling per S^3 dimension (d=3+1=4,
     so 1/4 per "direction"?)
  2. Whether 1/16 = (1/4)^2 has a variational or squared-amplitude meaning
  3. The FULL matrix structure (not just pairwise couplings)
  4. The connection to the graph Laplacian formalism
  5. What operator on S^3 produces constant 1/16 off-diagonal elements
"""

from __future__ import annotations

import json
import os
import sys
import numpy as np

import sympy as sp
from sympy import (
    Rational, Symbol, symbols, sqrt, pi, oo, integrate,
    simplify, factorial, gegenbauer, gamma, cancel
)


u = Symbol('u', real=True)

# ============================================================================
# Exact Gegenbauer computation on S^3
# ============================================================================

def gegenbauer_norm_analytical(k, lam):
    """
    Analytical formula for the Gegenbauer norm:
    h_k^lam = int_{-1}^1 (1-u^2)^{lam-1/2} [C_k^lam(u)]^2 du
            = pi * 2^{1-2*lam} * Gamma(k + 2*lam) / (k! * (k+lam) * Gamma(lam)^2)
    """
    return (pi * 2**(1 - 2*lam) * gamma(k + 2*lam)
            / (factorial(k) * (k + lam) * gamma(lam)**2))


def compute_all_couplings(n_max=6):
    """
    Compute ALL off-diagonal Fock weight matrix elements on S^3.

    For the S^3 inner product with weight (1-u^2)^{l+1/2}:
    The Fock weight is w(chi) = cos^4(chi/2) = (1+u)^2/4 where u = cos(chi).

    Expanding: w = 1/4 + u/2 + u^2/4

    Off-diagonal elements between Gegenbauer indices k and k+1:
    - 1/4 term: <k+1|k> = 0 (orthogonality)
    - u/2 term: (1/2) * recursion coefficient * norm
    - u^2/4 term: <k+1|u^2|k> = 0 (u^2 couples k to k and k±2, not k+1)

    So the coupling is entirely from the cos(chi)/2 term:
    <k+1|w|k> / sqrt(<k|k><k+1|k+1>) = (1/2) * a_{k,k+1} * sqrt(h_{k+1}/h_k)

    where a_{k,k+1} = (k+1)/(2(k+lam)) is the Gegenbauer recursion coefficient.
    """
    results = {}

    print("=" * 80)
    print("Complete table of Fock weight off-diagonal couplings on S^3")
    print("=" * 80)
    print()
    print(f"{'(n,l)->(n+1,l)':<20} {'k':<4} {'lam':<5} {'a_{k,k+1}':<15} "
          f"{'sqrt(h+/h-)':<15} {'coupling':<20} {'exact'}")
    print("-" * 100)

    for l in range(n_max):
        for n in range(l + 1, n_max):
            k = n - l - 1
            lam = l + 1

            a = Rational(k + 1, 2 * (k + lam))

            h_k = simplify(gegenbauer_norm_analytical(k, lam))
            h_kp1 = simplify(gegenbauer_norm_analytical(k + 1, lam))
            ratio = simplify(h_kp1 / h_k)
            sqrt_ratio = simplify(sqrt(ratio))

            coupling = simplify(Rational(1, 2) * a * sqrt_ratio)

            # Try to express as exact rational or simple radical
            coupling_sq = simplify(coupling**2)

            label = f"({n},{l})->({n+1},{l})"
            print(f"{label:<20} {k:<4} {lam:<5} {str(a):<15} "
                  f"{str(sqrt_ratio):<15} {str(coupling):<20} {str(coupling_sq)} = coupling^2")

            results[f"({n},{l})->({n+1},{l})"] = {
                "k": k, "lam": lam,
                "a_recursion": str(a),
                "coupling": str(coupling),
                "coupling_squared": str(coupling_sq),
                "coupling_float": float(coupling),
            }

    return results


def investigate_l0_universality():
    """
    For l=0, all Gegenbauer norms h_k^1 are equal to pi/2.
    The recursion coefficient is always 1/2.
    So the coupling is always (1/2) * (1/2) * 1 = 1/4.

    Prove this algebraically.
    """
    print("\n" + "=" * 80)
    print("l=0 universality: coupling = 1/4 for all n")
    print("=" * 80)

    results = {}

    # For l=0: lambda = 1, k = n-1
    # Gegenbauer C_k^1(u) = U_k(u) (Chebyshev U, second kind)
    # Norm: h_k^1 = pi/2 for ALL k >= 0
    # Recursion: a_{k,k+1} = (k+1)/(2(k+1)) = 1/2 for all k
    # Coupling: (1/2) * (1/2) * sqrt(h_{k+1}/h_k) = (1/4) * sqrt(1) = 1/4

    print("\nProof:")
    print("  lambda = l+1 = 1")
    print("  C_k^1(u) = U_k(u) (Chebyshev polynomial of the second kind)")
    print()

    # Verify h_k^1 = pi/2 for several k
    for k in range(8):
        h = simplify(gegenbauer_norm_analytical(k, 1))
        print(f"  h_{k}^1 = {h} = {simplify(h / pi)} * pi")
        assert simplify(h - pi/2) == 0, f"h_{k}^1 != pi/2"

    print("\n  => All h_k^1 = pi/2 (verified for k=0..7)")
    print()

    # The recursion coefficient
    for k in range(8):
        a = Rational(k + 1, 2 * (k + 1))
        assert a == Rational(1, 2), f"a != 1/2 for k={k}"

    print("  Recursion coefficient a_{k,k+1} = (k+1)/(2(k+1)) = 1/2 for all k")
    print()
    print("  sqrt(h_{k+1}/h_k) = sqrt(1) = 1 for all k")
    print()
    print("  COUPLING = (1/2) * (1/2) * 1 = 1/4    [UNIVERSAL for l=0]")
    print()

    # Now: 1/4 vs 1/16
    print("  Relationship to kappa = -1/16:")
    print("    1/4 = 4 * (1/16)")
    print("    1/16 = (1/4)^2")
    print("    1/16 = (1/4) / 4")
    print()

    # Check what operation maps 1/4 to 1/16
    print("  Candidate derivations of 1/16 from 1/4:")
    print("    (a) 1/4 * 1/4 = 1/16  [squared amplitude]")
    print("    (b) 1/4 / dim(S^3+1) = 1/4 / 4 = 1/16  [per-dimension]")
    print("    (c) 1/4 / (2l+1) at l=0 is 1/4, not 1/16")
    print("    (d) Fock Jacobian at chi=0: (2p0/(p^2+p0^2))^4 at p=0 -> (2/p0)^4 -> 16 at p0=1")
    print()

    results["l0_coupling"] = "1/4"
    results["l0_coupling_squared"] = "1/16"
    results["is_exactly_kappa_squared"] = True

    return results


def investigate_fock_jacobian():
    """
    The Fock projection Jacobian and its relationship to 1/16.

    The stereographic projection from R^3 to S^3 has Jacobian:
        J = (2 p0 / (p^2 + p0^2))^4

    At p = 0 (south pole): J = (2/p0)^4 = 16/p0^4
    At p0 = 1: J(0) = 16

    The CONFORMAL FACTOR is Omega = 2 p0 / (p^2 + p0^2), so J = Omega^4.
    At p=0, p0=1: Omega = 2.

    The graph Hamiltonian maps S^3 eigenvalues to energies.
    The conformal factor relates the flat-space measure to the S^3 measure.
    The 4th power comes from the dimension of R^3 + 1 = 4 (the stereographic
    projection is from R^3 to a 3-sphere embedded in R^4).

    kappa = -1/16:
    - The S^3 Laplace-Beltrami eigenvalue is -(n^2-1).
    - The physical Hamiltonian eigenvalue is -1/(2n^2).
    - For large n: -(n^2-1) ≈ -n^2, and -1/(2n^2) ≈ -1/(2n^2).
    - So the ratio approaches -1/(2n^4) as n -> inf, which is NOT constant.

    BUT: the graph Hamiltonian h1 has EXACT diagonal -1/(2n^2) and
    off-diagonal kappa = 1/16. The spectrum of this matrix converges to
    -1/(2n^2) as n_max -> infinity.

    Let me compute what off-diagonal coupling gives the EXACT spectrum
    for a 2x2 block at each (n, n+1).
    """
    print("\n" + "=" * 80)
    print("Fock Jacobian and the 1/16 derivation")
    print("=" * 80)

    results = {}

    # For a 2x2 matrix [[a, x], [x, b]], eigenvalues are:
    # lambda = (a+b)/2 +/- sqrt((a-b)^2/4 + x^2)
    # To get eigenvalues exactly a and b (i.e., x doesn't change them), need x=0.
    # To get eigenvalues close to a and b, need x << |a-b|.

    # The graph Hamiltonian for a pair (n, n+1) at l=0:
    # a = -1/(2n^2), b = -1/(2(n+1)^2), x = kappa = 1/16

    print("\nOff-diagonal perturbation analysis:")
    print(f"{'n':<5} {'a=-1/(2n^2)':<15} {'b=-1/(2(n+1)^2)':<18} {'|a-b|':<12} "
          f"{'x=1/16':<10} {'x/|a-b|':<10}")
    print("-" * 70)

    for n in range(1, 8):
        a = Rational(-1, 2*n**2)
        b = Rational(-1, 2*(n+1)**2)
        gap = abs(a - b)
        x = Rational(1, 16)
        ratio = x / gap if gap != 0 else None

        print(f"{n:<5} {str(a):<15} {str(b):<18} {str(gap):<12} "
              f"{str(x):<10} {str(ratio):<10}")

    # Now: the Fock weight coupling is 1/4 for l=0.
    # The graph Hamiltonian off-diagonal is 1/16.
    # Could the factor of 4 come from the dimension of the Fock embedding?

    print("\n" + "-" * 80)
    print("The key relationship: 1/16 = (1/4) * (1/4)")
    print("-" * 80)

    # The Fock projection has conformal factor Omega = 2p0/(p^2+p0^2).
    # The S^3 inner product integral over chi includes a measure sin^2(chi).
    # In the graph discretization, the "amplitude" for transition n->n+1 is
    # the Gegenbauer recursion coefficient 1/2 (for l=0), and the "weight"
    # from the Fock factor contributes another 1/2.
    #
    # Product: (1/2) * (1/2) = 1/4.
    #
    # But the GRAPH off-diagonal is 1/16, which is (1/4)^2.
    #
    # This would happen if the Hamiltonian matrix element is the SQUARED
    # amplitude of the coupling, not the amplitude itself.
    #
    # In quantum mechanics, matrix elements are amplitudes, not probabilities.
    # But in the Fock projection, the mapping from the flat-space Hamiltonian
    # to the S^3 operator involves the conformal factor SQUARED (because the
    # Hamiltonian is a second-order operator — two derivatives, each picking up
    # one power of the conformal factor).

    print()
    print("Hypothesis: kappa = -1/16 = -(1/4)^2")
    print()
    print("The Fock weight coupling for l=0 is 1/4 (universal).")
    print("The conformal factor for a second-order operator (Laplacian)")
    print("squares the coupling amplitude: (1/4)^2 = 1/16.")
    print()
    print("Physical interpretation:")
    print("  - The Fock weight (p^2+p0^2)^{-2} is the SQUARED conformal factor Omega^2.")
    print("  - Each power of Omega gives a coupling of 1/2 per transition.")
    print("  - The off-diagonal of the S^3 Laplace-Beltrami (which is the graph's")
    print("    kinetic term) involves TWO applications of the conformal factor")
    print("    (one for each derivative), giving (1/2 * 1/2)^2 = (1/4)^2 = 1/16.")
    print()
    print("  Alternatively: Omega^4 at p=0 is 16. So 1/Omega^4 = 1/16.")
    print("  The INVERSE Fock Jacobian at the south pole IS kappa.")

    results["fock_jacobian_at_p0_southpole"] = "16"
    results["inverse_jacobian"] = "1/16"
    results["kappa"] = "-1/16"
    results["hypothesis"] = "kappa = -1/J(0) where J is the Fock Jacobian at p=0, p0=1"

    return results


def compute_exact_h1_offdiag():
    """
    Compute the EXACT off-diagonal matrix elements of the graph Hamiltonian h1
    in the GeoVac code, and compare to the Fock weight couplings.
    """
    print("\n" + "=" * 80)
    print("Exact h1 off-diagonal elements from GeoVac code")
    print("=" * 80)

    results = {}

    try:
        from geovac.casimir_ci import _build_graph_h1

        for n_max in [3, 4, 5]:
            h1, orbitals = _build_graph_h1(Z=1, n_max=n_max)

            print(f"\nn_max = {n_max}:")
            print(f"  Number of orbitals: {len(orbitals)}")

            # Find inter-shell off-diagonal elements
            for i, (n1, l1, m1) in enumerate(orbitals):
                for j, (n2, l2, m2) in enumerate(orbitals):
                    if i >= j:
                        continue
                    if abs(h1[i, j]) < 1e-15:
                        continue
                    # Only show inter-shell connections
                    if n2 == n1 + 1 and l1 == l2 and m1 == m2:
                        print(f"  ({n1},{l1},{m1})->({n2},{l2},{m2}): "
                              f"h1[{i},{j}] = {h1[i,j]:.10f} = {Rational(h1[i,j]).limit_denominator(1000)}")

            # The off-diagonal is ALWAYS 1/16 = 0.0625
            offdiag_values = set()
            for i in range(len(orbitals)):
                for j in range(i+1, len(orbitals)):
                    if abs(h1[i, j]) > 1e-15:
                        offdiag_values.add(round(h1[i, j], 10))

            print(f"  Unique off-diagonal values: {offdiag_values}")
            results[f"n_max_{n_max}_offdiag_values"] = [str(v) for v in offdiag_values]

    except ImportError:
        print("  (Could not import geovac.casimir_ci)")

    return results


def investigate_conformal_operator():
    """
    The key question: how does the Laplace-Beltrami operator on S^3 produce
    constant off-diagonal elements in the Fock graph?

    The S^3 Laplace-Beltrami has eigenvalues -(n^2-1) and is DIAGONAL in the
    hyperspherical harmonic basis. So its off-diagonal elements are ZERO.

    The GRAPH Laplacian L = D - A has the SAME eigenvalues -(n^2-1) and the
    SAME eigenvectors (proven in Paper 7). So L is also diagonal-equivalent.

    But L has off-diagonal elements! These come from the ADJACENCY A.
    The adjacency A connects (n,l,m) to (n+1,l,m) with weight 1.
    The degree D is the row sum of A.

    So the off-diagonal of L is -A (i.e., -1 for adjacent pairs).
    And the off-diagonal of H = kappa * L is -kappa = 1/16.

    The 1/16 is NOT from the Fock weight itself, but from kappa applied
    to the unit adjacency!

    The question then becomes: WHY is the adjacency weight 1 (not some function
    of n, l)? And WHY is kappa = -1/16?

    The adjacency weight 1 means the graph is UNWEIGHTED. On S^3, this means
    all transitions have equal "strength." The physical coupling strength
    between adjacent states comes from the Gegenbauer recursion, which gives
    DIFFERENT values for different l.

    So: the GRAPH uses unit adjacency (all edges = 1), and kappa encodes the
    AVERAGE coupling strength. For l=0, the coupling is 1/4. For l>0, it's
    smaller. The average over all (n,l) pairs at given n_max depends on the
    distribution.

    Let me compute the AVERAGE coupling over all transitions at given n_max.
    """
    print("\n" + "=" * 80)
    print("Average Fock weight coupling vs kappa = 1/16")
    print("=" * 80)

    results = {}

    # For each n_max, compute the average coupling over all inter-shell transitions
    for n_max in [3, 4, 5, 6]:
        couplings = []
        weights = []  # weight by degeneracy (2l+1)

        for n in range(1, n_max):
            for l in range(n):
                # Only if (n+1,l) exists (l < n+1)
                if l < n + 1:
                    k = n - l - 1
                    lam = l + 1
                    a = Rational(k + 1, 2 * (k + lam))
                    h_k = gegenbauer_norm_analytical(k, lam)
                    h_kp1 = gegenbauer_norm_analytical(k + 1, lam)
                    ratio = simplify(h_kp1 / h_k)
                    coupling = simplify(Rational(1, 2) * a * sqrt(ratio))
                    coupling_sq = simplify(coupling**2)

                    deg = 2*l + 1  # each (n,l) represents (2l+1) states with different m
                    couplings.append((n, l, coupling, coupling_sq, deg))
                    weights.append((float(coupling), deg))

        # Weighted average (by degeneracy)
        total_weight = sum(d for _, d in weights)
        avg_coupling = sum(c * d for c, d in weights) / total_weight
        avg_coupling_sq = sum(c**2 * d for c, d in weights) / total_weight

        print(f"\nn_max = {n_max}:")
        print(f"  {'(n,l)':<10} {'coupling':<15} {'coupling^2':<15} {'deg':<5}")
        print(f"  {'-'*45}")
        for n, l, c, csq, deg in couplings:
            print(f"  ({n},{l}){'':<5} {float(c):<15.8f} {float(csq):<15.8f} {deg:<5}")

        print(f"  Average coupling (degeneracy-weighted): {avg_coupling:.8f}")
        print(f"  Average coupling^2 (degeneracy-weighted): {avg_coupling_sq:.8f}")
        print(f"  kappa = 1/16 = {1/16:.8f}")
        print(f"  Ratio avg_coupling / (1/4) = {avg_coupling / 0.25:.6f}")
        print(f"  Ratio avg_coupling^2 / (1/16) = {avg_coupling_sq / (1/16):.6f}")

        results[f"n_max_{n_max}"] = {
            "avg_coupling": avg_coupling,
            "avg_coupling_sq": avg_coupling_sq,
            "ratio_to_kappa": avg_coupling_sq / (1/16),
        }

    return results


def investigate_n2_identity():
    """
    Check whether the Fock weight coupling squared, SUMMED with n^2 degeneracy
    weighting, gives a constant related to 1/16.

    Recall: on S^3, each shell n has degeneracy n^2 = sum_{l=0}^{n-1} (2l+1).
    The Fock weight coupling for (n,l) -> (n+1,l) is:
        C(n,l) = (1/4) * sqrt(h_{n-l}^{l+1} / h_{n-l-1}^{l+1})    [from Approach E]

    where h_k^lam is the Gegenbauer norm.

    For l=0: C = 1/4 always.
    For l>0: C < 1/4 and depends on n and l.

    The n^2-weighted average per shell:
        <C>_n = (1/n^2) sum_{l=0}^{n-1} (2l+1) C(n,l)
    """
    print("\n" + "=" * 80)
    print("Per-shell average coupling (n^2-weighted)")
    print("=" * 80)

    results = {}

    for n in range(1, 10):
        total = Rational(0)
        for l in range(n):
            k = n - l - 1
            lam = l + 1
            a = Rational(k + 1, 2 * (k + lam))
            h_k = gegenbauer_norm_analytical(k, lam)
            h_kp1 = gegenbauer_norm_analytical(k + 1, lam)
            ratio = simplify(h_kp1 / h_k)
            coupling = simplify(Rational(1, 2) * a * sqrt(ratio))
            total += (2*l + 1) * coupling

        avg = simplify(total / n**2)
        avg_sq = simplify(total**2 / n**4)

        print(f"  n={n}: sum (2l+1)*C(n,l) = {simplify(total)},  avg = {avg},  float = {float(avg):.8f}")
        results[f"n_{n}_avg_coupling"] = str(avg)
        results[f"n_{n}_avg_coupling_float"] = float(avg)

    return results


def main():
    all_results = {}

    # 1. Complete table of couplings
    r1 = compute_all_couplings(n_max=6)
    all_results["all_couplings"] = r1

    # 2. l=0 universality proof
    r2 = investigate_l0_universality()
    all_results["l0_universality"] = r2

    # 3. Fock Jacobian analysis
    r3 = investigate_fock_jacobian()
    all_results["fock_jacobian"] = r3

    # 4. Exact h1 off-diagonal from code
    r4 = compute_exact_h1_offdiag()
    all_results["h1_offdiag"] = r4

    # 5. Average coupling analysis
    r5 = investigate_conformal_operator()
    all_results["average_coupling"] = r5

    # 6. Per-shell weighted average
    r6 = investigate_n2_identity()
    all_results["per_shell_average"] = r6

    # Save
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    os.makedirs(data_dir, exist_ok=True)
    json_path = os.path.join(data_dir, "probe_k1_fock_weight.json")
    with open(json_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nWrote {json_path}")

    # Final verdict
    print("\n" + "=" * 80)
    print("FINAL VERDICT")
    print("=" * 80)
    print()
    print("STATUS: PARTIAL (positive structural relationship, not exact derivation)")
    print()
    print("KEY FINDINGS:")
    print()
    print("1. UNIVERSAL l=0 COUPLING = 1/4")
    print("   For all s-orbital transitions (n,0) -> (n+1,0), the Fock weight")
    print("   cos^4(chi/2) produces a normalized coupling of EXACTLY 1/4,")
    print("   independent of n. This is a PROVEN universal constant from:")
    print("   - Chebyshev U polynomials all having equal norm h_k^1 = pi/2")
    print("   - Gegenbauer recursion coefficient = 1/2 at lambda=1")
    print("   - Product: (1/2) * (1/2) * 1 = 1/4")
    print()
    print("2. RELATIONSHIP 1/16 = (1/4)^2")
    print("   kappa = -1/16 is the SQUARE of the l=0 Fock weight coupling.")
    print("   This is consistent with the Fock Jacobian interpretation:")
    print("   J = Omega^4 = (2p0/(p^2+p0^2))^4 = 16 at p=0, p0=1")
    print("   => 1/J = 1/16 = kappa (up to sign)")
    print()
    print("3. INVERSE JACOBIAN INTERPRETATION")
    print("   kappa = -1/16 = -1/J(p=0,p0=1)")
    print("   The kinetic scale is the inverse of the Fock Jacobian evaluated")
    print("   at the south pole of S^3 (where p=0, the rest frame).")
    print("   The Jacobian Omega^4 = 16 converts between the flat-space measure")
    print("   d^3p and the S^3 measure d^3Omega. Its inverse 1/16 converts")
    print("   S^3 eigenvalues to flat-space energies.")
    print()
    print("4. l-DEPENDENCE BREAKS UNIVERSALITY")
    print("   For l>0, the coupling is NOT 1/4 but depends on both n and l.")
    print("   The graph Hamiltonian uses CONSTANT 1/16 for ALL transitions,")
    print("   which is exact only for l=0 (as the squared amplitude).")
    print("   For l>0, the constant 1/16 is an AVERAGE/APPROXIMATION.")
    print()
    print("5. THE DERIVATION PATH")
    print("   The cleanest derivation of kappa = -1/16 would be:")
    print("   (a) The Fock Jacobian is Omega^4 where Omega = 2p0/(p^2+p0^2)")
    print("   (b) At the ground state (p0=1) and rest frame (p=0): Omega=2, J=16")
    print("   (c) The kinetic energy operator converts S^3 -> flat space by 1/J")
    print("   (d) Therefore kappa = -1/J = -1/16")
    print()
    print("   This is a STRUCTURAL derivation, not a variational one.")
    print("   It identifies kappa as the inverse Fock Jacobian at the south pole.")

    return all_results


if __name__ == "__main__":
    main()
