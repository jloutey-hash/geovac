"""
Dirac native graph spectrum exploration.
========================================

Does a graph on Dirac nodes (n, kappa, m_j) with natural adjacency rules
reproduce the Dirac spectrum on S^3?

The scalar graph (Paper 7) has graph Laplacian eigenvalues n^2 - 1 matching
the Laplace-Beltrami on S^3. This script asks the analogous question for
the Dirac case using two adjacency rules from geovac.ihara_zeta_dirac.

Author: GeoVac exploration, April 2026.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from typing import Dict, List, Tuple

import numpy as np

# ---- GeoVac imports ----
from geovac.ihara_zeta_dirac import build_dirac_s3_graph
from geovac.lattice import GeometricLattice
from geovac.dirac_matrix_elements import DiracLabel, kappa_to_l, kappa_to_j, iter_dirac_labels

np.set_printoptions(precision=12, linewidth=120, suppress=True)


# ============================================================================
# Utility
# ============================================================================

def eigenvalues_of_laplacian(A: np.ndarray) -> np.ndarray:
    """Graph Laplacian L = D - A eigenvalues, sorted."""
    D = np.diag(A.sum(axis=1))
    L = D - A
    evals = np.linalg.eigvalsh(L.astype(float))
    return np.sort(evals)


def eigenvalues_of_weighted_laplacian(A_weighted: np.ndarray) -> np.ndarray:
    """Weighted graph Laplacian: D_w - A_w where D_w = diag(sum of weights)."""
    D = np.diag(A_weighted.sum(axis=1))
    L = D - A_weighted
    evals = np.linalg.eigvalsh(L)
    return np.sort(evals)


def group_eigenvalues(evals: np.ndarray, tol: float = 1e-8) -> List[Tuple[float, int]]:
    """Group eigenvalues by value (within tolerance). Return (value, degeneracy)."""
    if len(evals) == 0:
        return []
    groups = []
    current_val = evals[0]
    current_count = 1
    for i in range(1, len(evals)):
        if abs(evals[i] - current_val) < tol:
            current_count += 1
        else:
            groups.append((float(current_val), current_count))
            current_val = evals[i]
            current_count = 1
    groups.append((float(current_val), current_count))
    # Refine: use mean of group
    refined = []
    idx = 0
    for val, deg in groups:
        group_vals = evals[idx:idx+deg]
        refined.append((float(np.mean(group_vals)), deg))
        idx += deg
    return refined


def camporesi_higuchi_dirac_eigenvalues(n_max_CH: int) -> List[Tuple[float, int]]:
    """Camporesi-Higuchi eigenvalues |lambda_n| = n + 3/2, g_n = 2(n+1)(n+2)
    for n = 0, 1, ..., n_max_CH.
    In Fock convention n_fock = n_CH + 1, so |lambda| = n_fock + 1/2."""
    result = []
    for n in range(n_max_CH + 1):
        lam = n + 1.5
        g = 2 * (n + 1) * (n + 2)
        result.append((lam, g))
    return result


def camporesi_higuchi_squared(n_max_CH: int) -> List[Tuple[float, int]]:
    """Squared CH eigenvalues: (n + 3/2)^2."""
    result = []
    for n in range(n_max_CH + 1):
        lam_sq = (n + 1.5) ** 2
        g = 2 * (n + 1) * (n + 2)
        result.append((lam_sq, g))
    return result


def scalar_laplace_beltrami_eigenvalues(n_max_fock: int) -> List[Tuple[float, int]]:
    """Scalar LB eigenvalues: n^2 - 1, g = n^2 for n = 1, ..., n_max_fock."""
    result = []
    for n in range(1, n_max_fock + 1):
        result.append((float(n**2 - 1), n**2))
    return result


# ============================================================================
# Scalar graph reference
# ============================================================================

def scalar_graph_spectrum(n_max: int) -> Tuple[np.ndarray, List[Tuple[float, int]]]:
    """Build scalar graph (n,l,m) and compute L = D - A eigenvalues."""
    lattice = GeometricLattice(max_n=n_max, nuclear_charge=1, topological_weights=False)
    A = lattice.adjacency.toarray()
    evals = eigenvalues_of_laplacian(A)
    groups = group_eigenvalues(evals)
    return evals, groups


# ============================================================================
# Weighted Dirac graph
# ============================================================================

def build_weighted_dirac_graph(n_max: int, rule: str) -> np.ndarray:
    """Build Dirac graph with edge weights 1/(n1 * n2)."""
    A_int, labels, deg, desc = build_dirac_s3_graph(n_max, rule)
    A_weighted = np.zeros_like(A_int, dtype=float)
    V = len(labels)
    for i in range(V):
        for j in range(i + 1, V):
            if A_int[i, j] == 1:
                n1 = labels[i].n_fock
                n2 = labels[j].n_fock
                w = 1.0 / (n1 * n2)
                A_weighted[i, j] = w
                A_weighted[j, i] = w
    return A_weighted, labels


# ============================================================================
# Main analysis
# ============================================================================

def analyze_graph(name: str, evals: np.ndarray, groups: List[Tuple[float, int]],
                  target_evals: dict, n_max: int):
    """Print detailed analysis of a graph's spectrum."""
    print(f"\n{'='*80}")
    print(f"  {name}")
    print(f"{'='*80}")
    print(f"  Total nodes: {len(evals)}")
    print(f"  Number of distinct eigenvalue groups: {len(groups)}")
    print(f"\n  Eigenvalue spectrum (value, degeneracy):")
    for val, deg in groups:
        print(f"    {val:16.10f}  deg = {deg}")

    # Zero eigenvalue count (connected components)
    zero_count = sum(deg for val, deg in groups if abs(val) < 1e-8)
    print(f"\n  Zero eigenvalues (connected components): {zero_count}")

    # Compare with target spectra
    print(f"\n  --- Comparisons ---")

    # CH Dirac |lambda_n| = n + 3/2
    if "CH_dirac" in target_evals:
        print(f"\n  vs Camporesi-Higuchi |lambda_n| = n+3/2:")
        for (ch_val, ch_deg) in target_evals["CH_dirac"]:
            match = [(v, d) for v, d in groups if abs(v - ch_val) < 0.01]
            if match:
                print(f"    CH {ch_val:8.4f} (deg {ch_deg:4d}) <-> graph {match[0][0]:12.8f} (deg {match[0][1]:4d})")
            else:
                print(f"    CH {ch_val:8.4f} (deg {ch_deg:4d}) <-> NO MATCH")

    # CH Dirac squared
    if "CH_dirac_sq" in target_evals:
        print(f"\n  vs Camporesi-Higuchi squared (n+3/2)^2:")
        for (ch_val, ch_deg) in target_evals["CH_dirac_sq"]:
            match = [(v, d) for v, d in groups if abs(v - ch_val) < 0.01]
            if match:
                print(f"    CH^2 {ch_val:8.4f} (deg {ch_deg:4d}) <-> graph {match[0][0]:12.8f} (deg {match[0][1]:4d})")
            else:
                print(f"    CH^2 {ch_val:8.4f} (deg {ch_deg:4d}) <-> NO MATCH")

    # Scalar LB
    if "scalar_LB" in target_evals:
        print(f"\n  vs Scalar Laplace-Beltrami n^2 - 1:")
        for (lb_val, lb_deg) in target_evals["scalar_LB"]:
            match = [(v, d) for v, d in groups if abs(v - lb_val) < 0.01]
            if match:
                print(f"    LB {lb_val:8.4f} (deg {lb_deg:4d}) <-> graph {match[0][0]:12.8f} (deg {match[0][1]:4d})")
            else:
                print(f"    LB {lb_val:8.4f} (deg {lb_deg:4d}) <-> NO MATCH")

    # Ratio analysis: for each nonzero eigenvalue, compute ratio to every reference
    print(f"\n  --- Ratio analysis (graph_eig / reference) ---")
    nonzero_groups = [(v, d) for v, d in groups if abs(v) > 1e-8]
    if nonzero_groups and "scalar_LB" in target_evals:
        scalar_vals = [v for v, d in target_evals["scalar_LB"] if v > 0]
        if scalar_vals:
            print(f"\n  Ratios of sorted nonzero graph eigenvalues / sorted scalar LB eigenvalues:")
            scalar_flat = []
            for v, d in target_evals["scalar_LB"]:
                if v > 0:
                    scalar_flat.extend([v] * d)
            scalar_flat.sort()
            dirac_flat = []
            for v, d in nonzero_groups:
                dirac_flat.extend([v] * d)
            dirac_flat.sort()
            # Show ratios for unique eigenvalue groups
            for v, d in nonzero_groups:
                # Find closest scalar LB eigenvalue
                closest_idx = np.argmin([abs(v - sv) for sv in scalar_vals])
                closest = scalar_vals[closest_idx]
                ratio = v / closest if closest != 0 else float('inf')
                print(f"    graph {v:12.8f} (deg {d}) / scalar_LB {closest:8.4f} = {ratio:10.6f}")

    # Look for integer or simple rational patterns
    print(f"\n  --- Pattern search (graph eigenvalues) ---")
    for v, d in groups:
        if abs(v) < 1e-8:
            continue
        # Check if close to integer
        nearest_int = round(v)
        if abs(v - nearest_int) < 0.001:
            print(f"    {v:12.8f} (deg {d}) ~ integer {nearest_int}")
        # Check if close to n^2 - 1
        for n_test in range(1, n_max + 5):
            if abs(v - (n_test**2 - 1)) < 0.01:
                print(f"    {v:12.8f} (deg {d}) ~ n^2-1 at n={n_test} (scalar LB)")
        # Check if close to (n+3/2)^2
        for n_test in range(0, n_max + 5):
            if abs(v - (n_test + 1.5)**2) < 0.01:
                print(f"    {v:12.8f} (deg {d}) ~ (n+3/2)^2 at n_CH={n_test}")
        # Check if close to n+3/2
        for n_test in range(0, n_max + 5):
            if abs(v - (n_test + 1.5)) < 0.01:
                print(f"    {v:12.8f} (deg {d}) ~ n+3/2 at n_CH={n_test}")
        # Check if close to (n+3/2)^2 - 9/4
        for n_test in range(0, n_max + 5):
            target = (n_test + 1.5)**2 - 2.25
            if abs(v - target) < 0.01:
                print(f"    {v:12.8f} (deg {d}) ~ (n+3/2)^2 - 9/4 = n^2+3n at n_CH={n_test}")
        # Check n(n+3) which is (n+3/2)^2 - 9/4
        for n_test in range(0, n_max + 5):
            target = n_test * (n_test + 3)
            if abs(v - target) < 0.01:
                print(f"    {v:12.8f} (deg {d}) ~ n(n+3) at n_CH={n_test}")
        # Check 2*n^2 - 2, 4*n^2 - 4 etc
        for fac in [2, 3, 4, 0.5, 1.5]:
            for n_test in range(1, n_max + 5):
                target = fac * (n_test**2 - 1)
                if abs(target) > 0.01 and abs(v - target) < 0.01:
                    print(f"    {v:12.8f} (deg {d}) ~ {fac}*(n^2-1) at n={n_test}")


def compare_dirac_degeneracies(groups: List[Tuple[float, int]], n_max: int):
    """Compare degeneracy pattern with Dirac g_n = 2(n+1)(n+2) and scalar g_n = n^2."""
    print(f"\n  --- Degeneracy comparison ---")
    dirac_degs = {2*(n+1)*(n+2): n for n in range(n_max + 3)}
    scalar_degs = {n**2: n for n in range(1, n_max + 3)}
    for val, deg in groups:
        notes = []
        if deg in dirac_degs:
            notes.append(f"= g_Dirac(n_CH={dirac_degs[deg]})")
        if deg in scalar_degs:
            notes.append(f"= g_scalar(n={scalar_degs[deg]})")
        note_str = ", ".join(notes) if notes else "(no match to standard degeneracies)"
        print(f"    eigenvalue {val:12.8f}: deg = {deg:4d}  {note_str}")


def main():
    results = {}

    for n_max in [1, 2, 3, 4]:
        print(f"\n\n{'#'*80}")
        print(f"#  n_max = {n_max}")
        print(f"{'#'*80}")

        n_max_CH = n_max - 1  # CH convention: n_CH = 0, 1, ..., n_max-1

        # Reference spectra
        ch_dirac = camporesi_higuchi_dirac_eigenvalues(n_max_CH)
        ch_dirac_sq = camporesi_higuchi_squared(n_max_CH)
        scalar_lb = scalar_laplace_beltrami_eigenvalues(n_max)

        target = {
            "CH_dirac": ch_dirac,
            "CH_dirac_sq": ch_dirac_sq,
            "scalar_LB": scalar_lb,
        }

        print(f"\n  Reference: Camporesi-Higuchi Dirac |lambda_n| = n+3/2 (CH convention):")
        for v, d in ch_dirac:
            print(f"    n_CH={ch_dirac.index((v,d))}: |lambda| = {v}, g = {d}")

        print(f"\n  Reference: Scalar LB n^2-1:")
        for v, d in scalar_lb:
            print(f"    n={scalar_lb.index((v,d))+1}: lambda = {v}, g = {d}")

        # ---- Scalar graph ----
        s_evals, s_groups = scalar_graph_spectrum(n_max)
        analyze_graph(f"Scalar graph (n,l,m), n_max={n_max}", s_evals, s_groups, target, n_max)
        compare_dirac_degeneracies(s_groups, n_max)

        # ---- Dirac graph Rule A ----
        A_A, labels_A, deg_A, desc_A = build_dirac_s3_graph(n_max, "A")
        print(f"\n  {desc_A}")
        evals_A = eigenvalues_of_laplacian(A_A)
        groups_A = group_eigenvalues(evals_A)
        analyze_graph(f"Dirac graph Rule A (kappa-preserving), n_max={n_max}", evals_A, groups_A, target, n_max)
        compare_dirac_degeneracies(groups_A, n_max)

        # ---- Dirac graph Rule B ----
        A_B, labels_B, deg_B, desc_B = build_dirac_s3_graph(n_max, "B")
        print(f"\n  {desc_B}")
        evals_B = eigenvalues_of_laplacian(A_B)
        groups_B = group_eigenvalues(evals_B)
        analyze_graph(f"Dirac graph Rule B (E1 dipole), n_max={n_max}", evals_B, groups_B, target, n_max)
        compare_dirac_degeneracies(groups_B, n_max)

        # ---- Weighted Dirac graph Rule A (1/(n1*n2)) ----
        if n_max >= 2:  # weighting only interesting when we have multiple n
            A_wA, labels_wA = build_weighted_dirac_graph(n_max, "A")
            evals_wA = eigenvalues_of_weighted_laplacian(A_wA)
            groups_wA = group_eigenvalues(evals_wA)
            analyze_graph(f"Weighted Dirac graph Rule A (w=1/(n1*n2)), n_max={n_max}",
                         evals_wA, groups_wA, target, n_max)

            # ---- Weighted Dirac graph Rule B (1/(n1*n2)) ----
            A_wB, labels_wB = build_weighted_dirac_graph(n_max, "B")
            evals_wB = eigenvalues_of_weighted_laplacian(A_wB)
            groups_wB = group_eigenvalues(evals_wB)
            analyze_graph(f"Weighted Dirac graph Rule B (w=1/(n1*n2)), n_max={n_max}",
                         evals_wB, groups_wB, target, n_max)

        # ---- Store results ----
        nmax_key = f"n_max={n_max}"
        results[nmax_key] = {
            "scalar": {
                "eigenvalues": [float(x) for x in s_evals],
                "groups": [(float(v), d) for v, d in s_groups],
            },
            "dirac_rule_A": {
                "description": desc_A,
                "n_nodes": int(len(labels_A)),
                "n_edges": int(A_A.sum()) // 2,
                "eigenvalues": [float(x) for x in evals_A],
                "groups": [(float(v), d) for v, d in groups_A],
            },
            "dirac_rule_B": {
                "description": desc_B,
                "n_nodes": int(len(labels_B)),
                "n_edges": int(A_B.sum()) // 2,
                "eigenvalues": [float(x) for x in evals_B],
                "groups": [(float(v), d) for v, d in groups_B],
            },
        }
        if n_max >= 2:
            results[nmax_key]["weighted_rule_A"] = {
                "eigenvalues": [float(x) for x in evals_wA],
                "groups": [(float(v), d) for v, d in groups_wA],
            }
            results[nmax_key]["weighted_rule_B"] = {
                "eigenvalues": [float(x) for x in evals_wB],
                "groups": [(float(v), d) for v, d in groups_wB],
            }

    # ============================================================================
    # Cross-n_max ratio analysis: do Dirac graph eigenvalues follow any formula?
    # ============================================================================
    print(f"\n\n{'#'*80}")
    print(f"#  CROSS-n_max ANALYSIS")
    print(f"{'#'*80}")

    for rule in ["A", "B"]:
        print(f"\n\n  Rule {rule} — Nonzero eigenvalue groups across n_max:")
        for n_max in [1, 2, 3, 4]:
            key = f"n_max={n_max}"
            groups = results[key][f"dirac_rule_{rule}"]["groups"]
            nonzero = [(v, d) for v, d in groups if abs(v) > 1e-8]
            print(f"    n_max={n_max}: {nonzero}")

    # Detailed: for Rule A, check if per-kappa sub-blocks give scalar-like spectra
    print(f"\n\n  Rule A: per-kappa sector decomposition")
    for n_max in [2, 3, 4]:
        print(f"\n    n_max = {n_max}:")
        A_A, labels_A, deg_A, desc_A = build_dirac_s3_graph(n_max, "A")
        # Group labels by kappa
        kappa_values = sorted(set(lab.kappa for lab in labels_A))
        for kappa in kappa_values:
            indices = [i for i, lab in enumerate(labels_A) if lab.kappa == kappa]
            if len(indices) <= 1:
                continue
            sub_A = A_A[np.ix_(indices, indices)]
            sub_evals = eigenvalues_of_laplacian(sub_A)
            sub_groups = group_eigenvalues(sub_evals)
            l = kappa_to_l(kappa)
            j = kappa_to_j(kappa)
            print(f"      kappa={kappa:+d} (l={l}, j={j}): nodes={len(indices)}, "
                  f"evals = {sub_groups}")

            # For each kappa sector: does the sub-graph Laplacian have eigenvalues
            # matching the SCALAR graph restricted to a single l-shell?
            # Scalar l-shell: nodes (n, l, m) for all n, m. The l-shell adjacency
            # is a 2D grid: n-direction (n_min..n_max) × m-direction (-l..l).
            # Graph Laplacian eigenvalues of the path P_k are 2 - 2*cos(pi*j/k).
            # The kappa-sector is analogous: n-direction × m_j-direction.
            n_range = sorted(set(labels_A[i].n_fock for i in indices))
            mj_range = sorted(set(labels_A[i].two_m_j for i in indices))
            print(f"        n in {n_range}, 2*m_j in {mj_range}")

    # ============================================================================
    # Try alternative: Dirac eigenvalues as (n+3/2)^2 - 9/4 = n(n+3)
    # This is |lambda_n|^2 - lambda_0^2, the "subtracted" eigenvalue
    # ============================================================================
    print(f"\n\n  --- Alternative target: n_CH*(n_CH+3) = |lambda|^2 - 9/4 ---")
    for n_max in [2, 3, 4]:
        key = f"n_max={n_max}"
        for rule in ["A", "B"]:
            groups = results[key][f"dirac_rule_{rule}"]["groups"]
            nonzero = [(v, d) for v, d in groups if abs(v) > 1e-8]
            targets = [n*(n+3) for n in range(n_max + 3)]
            print(f"\n    Rule {rule}, n_max={n_max}:")
            for v, d in nonzero:
                found = False
                for t in targets:
                    if t > 0 and abs(v - t) < 0.01:
                        n_ch = targets.index(t)
                        print(f"      {v:12.8f} (deg {d}) MATCHES n_CH(n_CH+3) = {t} at n_CH = {n_ch}")
                        found = True
                if not found:
                    # Check fractions of known targets
                    for t in targets:
                        if t > 0:
                            ratio = v / t
                            if abs(ratio - round(ratio)) < 0.001 and round(ratio) > 0:
                                print(f"      {v:12.8f} (deg {d}) = {round(ratio)} * n_CH(n_CH+3) at n_CH = {targets.index(t)}")
                                found = True
                    if not found:
                        # check v/integer
                        for div in range(1, 10):
                            for t in targets:
                                if t > 0 and abs(v * div - t) < 0.01:
                                    n_ch = targets.index(t)
                                    print(f"      {v:12.8f} (deg {d}) = n_CH(n_CH+3)/{div} = {t}/{div} at n_CH={n_ch}")
                                    found = True

    # ============================================================================
    # Try: kappa * (D-A) as an analog of the Dirac operator (first-order, signed)
    # ============================================================================
    print(f"\n\n  --- Signed Dirac-like operator: multiply each node by sign(kappa) ---")
    print(f"  Idea: Dirac operator is first-order with +/- eigenvalues.")
    print(f"  Try S * L where S = diag(sign(kappa)) and L = D - A")
    for n_max in [2, 3, 4]:
        A_A, labels_A, deg_A, desc_A = build_dirac_s3_graph(n_max, "A")
        L = np.diag(A_A.sum(axis=1).astype(float)) - A_A.astype(float)
        S = np.diag([1.0 if lab.kappa < 0 else -1.0 for lab in labels_A])
        SL = S @ L
        evals_SL = np.sort(np.linalg.eigvals(SL).real)
        groups_SL = group_eigenvalues(evals_SL, tol=1e-6)
        print(f"\n    Rule A n_max={n_max}: sign(kappa) * L eigenvalues:")
        for v, d in groups_SL:
            print(f"      {v:12.8f}  deg = {d}")

        # Also try with Rule B
        A_B, labels_B, deg_B, desc_B = build_dirac_s3_graph(n_max, "B")
        L_B = np.diag(A_B.sum(axis=1).astype(float)) - A_B.astype(float)
        S_B = np.diag([1.0 if lab.kappa < 0 else -1.0 for lab in labels_B])
        SL_B = S_B @ L_B
        evals_SLB = np.sort(np.linalg.eigvals(SL_B).real)
        groups_SLB = group_eigenvalues(evals_SLB, tol=1e-6)
        print(f"    Rule B n_max={n_max}: sign(kappa) * L eigenvalues:")
        for v, d in groups_SLB:
            print(f"      {v:12.8f}  deg = {d}")

    # ============================================================================
    # SAVE RESULTS
    # ============================================================================
    print(f"\n\nSaving results to debug/data/dirac_native_graph_spectrum.json")
    with open("debug/data/dirac_native_graph_spectrum.json", "w") as f:
        json.dump(results, f, indent=2)
    print("Done.")


if __name__ == "__main__":
    main()
