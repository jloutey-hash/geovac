"""
Fock projection of the flat-space transverse photon propagator onto S^3.
=========================================================================

Tests what happens when the flat-space transverse photon propagator is
Fock-projected onto S^3, and compares the result to the scalar edge
Laplacian pseudoinverse L_1^+ that graph-native QED currently uses.

KEY HYPOTHESIS: The graph's L_1 = B^T B is only the LONGITUDINAL photon
sector (exact 1-forms). The TRANSVERSE photon (physical) requires
plaquettes/faces from Paper 30's Wilson gauge construction, which complete
the Hodge Laplacian Delta_1 = d_0^T d_0 + d_1 d_1^T.

Parts:
  1. Vector harmonics on S^3: spectrum of exact and co-exact 1-forms
  2. Compare L_1 on the Fock graph to the continuum 1-form Laplacian
  3. Analyze the Hodge decomposition on the 1-complex (no 2-cells)
  4. Test the plaquette hypothesis: do Paper 30's plaquettes provide d_1?
  5. Spectrum comparison summary

Author: GeoVac research computation (2026-05-01)
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from geovac.fock_graph_hodge import FockGraphHodge
from geovac.lattice import GeometricLattice
from geovac.hodge1_s3 import hodge1_eigenvalue, hodge1_degeneracy


def continuum_spectra(n_max: int) -> Dict:
    """Compute continuum spectra for comparison.

    On unit S^3:
      Scalar Laplacian: lambda_k = k(k+2) for k = 0, 1, 2, ...
        (In GeoVac Fock convention: lambda_n = n^2 - 1 for n = 1, 2, ...)
      Hodge-1 Laplacian: mu_n = n(n+2) for n = 1, 2, 3, ...

    Key relation: the scalar eigenvalue at level n in the Fock convention
    is n^2 - 1 = (n-1)(n+1). The Hodge-1 eigenvalue at level n is n(n+2).
    These are DIFFERENT: n(n+2) = n^2+2n vs n^2-1.

    However, the NON-ZERO scalar eigenvalues are n^2-1 for n >= 2, i.e.,
    3, 8, 15, 24, ... The Hodge-1 eigenvalues are n(n+2) for n >= 1, i.e.,
    3, 8, 15, 24, ... They MATCH at each value but with offset indices!

    scalar lambda_2 = 3 = Hodge-1 mu_1
    scalar lambda_3 = 8 = Hodge-1 mu_2
    scalar lambda_4 = 15 = Hodge-1 mu_3

    Actually: k^2 - 1 at k = n+1 gives (n+1)^2 - 1 = n^2 + 2n = n(n+2).
    So the Hodge-1 eigenvalue mu_n equals the scalar eigenvalue at level n+1.
    This is the level-shift relationship documented in hodge1_s3.py.
    """
    scalar_eigenvalues = {}  # n -> eigenvalue
    hodge1_eigenvalues = {}  # n -> eigenvalue

    for n in range(1, n_max + 2):
        scalar_eigenvalues[n] = n * n - 1  # n^2 - 1 (Fock convention)

    for n in range(1, n_max + 1):
        hodge1_eigenvalues[n] = n * (n + 2)  # n(n+2)

    return {
        "scalar_eigenvalues": scalar_eigenvalues,
        "hodge1_eigenvalues": hodge1_eigenvalues,
        "level_shift_verification": {
            n: {
                "mu_n_hodge1": n * (n + 2),
                "lambda_{n+1}_scalar": (n + 1) ** 2 - 1,
                "match": n * (n + 2) == (n + 1) ** 2 - 1,
            }
            for n in range(1, n_max + 1)
        },
    }


def part1_vector_harmonics(n_max: int) -> Dict:
    """Part 1: Vector harmonics on S^3 — exact and co-exact decomposition.

    On S^3 at eigenvalue mu_n = n(n+2):
      Total 1-form modes: 2n(n+2)
      Exact (longitudinal, pure gauge): n(n+2)
      Co-exact (transverse, physical photon): n(n+2)

    The exact and co-exact sectors have EQUAL multiplicities on S^3
    because the Hodge star provides an isomorphism between them.

    The transverse photon propagator in mode space:
      G_T(n) = 1 / [n(n+2)]   for the transverse (co-exact) sector
    """
    results = {
        "description": "Vector harmonics on unit S^3",
        "levels": [],
    }

    for n in range(1, n_max + 1):
        mu_n = n * (n + 2)
        d_total = 2 * n * (n + 2)
        d_transverse = n * (n + 2)  # co-exact
        d_longitudinal = n * (n + 2)  # exact
        propagator = 1.0 / mu_n

        results["levels"].append({
            "n": n,
            "eigenvalue_mu_n": mu_n,
            "degeneracy_total": d_total,
            "degeneracy_transverse": d_transverse,
            "degeneracy_longitudinal": d_longitudinal,
            "transverse_propagator": propagator,
        })

    return results


def part2_compare_L1_to_continuum(n_max_values: List[int]) -> Dict:
    """Part 2: Compare edge Laplacian L_1 spectrum to continuum Hodge-1.

    Build the Fock graph at each n_max, compute L_1 = B^T B eigenvalues,
    and compare to both:
      (a) scalar continuum n^2 - 1
      (b) Hodge-1 continuum n(n+2)

    The SVD theorem says nonzero eigenvalues of L_1 = B^T B equal those
    of L_0 = B B^T (the node Laplacian). Paper 7 established that L_0
    approximates the scalar Laplacian on S^3.
    """
    results = {
        "description": "L_1 spectrum vs continuum spectra",
        "graphs": {},
    }

    for n_max in n_max_values:
        print(f"  Building Fock graph at n_max={n_max}...")
        hodge = FockGraphHodge(n_max)
        V = hodge.n_nodes
        E = hodge.n_edges
        beta_0 = hodge.betti_0
        beta_1 = hodge.betti_1

        # Get L_1 spectrum
        L1_evals = np.sort(np.linalg.eigvalsh(hodge.edge_laplacian_numpy))
        L0_evals = np.sort(np.linalg.eigvalsh(hodge.node_laplacian_numpy))

        # Nonzero eigenvalues
        nz_L1 = sorted(L1_evals[L1_evals > 0.5])
        nz_L0 = sorted(L0_evals[L0_evals > 0.5])

        # SVD check: nonzero L0 and L1 should match
        svd_match = (len(nz_L0) == len(nz_L1))
        if svd_match and len(nz_L0) > 0:
            svd_max_diff = max(abs(a - b) for a, b in zip(nz_L0, nz_L1))
        else:
            svd_max_diff = None

        # Compare L0/L1 nonzero eigenvalues to scalar continuum n^2-1
        scalar_continuum = [n * n - 1 for n in range(2, n_max + 2)]
        hodge1_continuum = [n * (n + 2) for n in range(1, n_max + 1)]

        # The Fock graph's L0 eigenvalues: for a PERFECT Fock graph,
        # they would be exactly n^2-1 with degeneracy n^2.
        # But the actual graph is truncated, so we compare.
        L0_distinct = sorted(set(int(round(x)) for x in nz_L0))

        comparison = {
            "n_max": n_max,
            "V": V,
            "E": E,
            "beta_0": beta_0,
            "beta_1": beta_1,
            "n_L1_nonzero": len(nz_L1),
            "n_L1_zero": E - len(nz_L1),
            "svd_L0_L1_match": svd_match,
            "svd_max_diff": svd_max_diff,
            "L0_distinct_nonzero": L0_distinct,
            "scalar_continuum_n2m1": scalar_continuum[:len(L0_distinct)],
            "hodge1_continuum_n_np2": hodge1_continuum[:len(L0_distinct)],
            "L0_vs_scalar": [],
            "L0_vs_hodge1": [],
        }

        for i, ev in enumerate(L0_distinct):
            # Compare to scalar n^2-1 and Hodge-1 n(n+2)
            # Find nearest scalar level
            best_scalar = min(scalar_continuum, key=lambda x: abs(x - ev))
            best_hodge1 = min(hodge1_continuum, key=lambda x: abs(x - ev))
            comparison["L0_vs_scalar"].append({
                "L0_eigenvalue": ev,
                "nearest_scalar": best_scalar,
                "gap": ev - best_scalar,
            })
            comparison["L0_vs_hodge1"].append({
                "L0_eigenvalue": ev,
                "nearest_hodge1": best_hodge1,
                "gap": ev - best_hodge1,
            })

        # Full L1 eigenvalue list
        comparison["L1_all_eigenvalues"] = [round(x, 6) for x in L1_evals.tolist()]
        comparison["L0_all_eigenvalues"] = [round(x, 6) for x in L0_evals.tolist()]

        results["graphs"][str(n_max)] = comparison

    return results


def part3_hodge_decomposition_analysis(n_max_values: List[int]) -> Dict:
    """Part 3: Analyze the Hodge decomposition on the 1-complex.

    On a graph (simplicial 1-complex, no 2-cells):
      - The edge space R^E decomposes as: im(B^T) + ker(L_1)
      - im(B^T) = exact 1-forms (d_0 applied to 0-cochains)
        Dimension: rank(B) = V - beta_0
      - ker(L_1) = harmonic 1-forms (gauge zero modes)
        Dimension: beta_1 = E - V + beta_0

    THERE IS NO CO-EXACT SECTOR on a graph (1-complex)!
    The co-exact operator d_1^T requires 2-cells (faces/plaquettes).
    On a 1-complex, d_1 does not exist, so d_1 d_1^T = 0.

    Therefore:
      L_1 = B^T B = d_0^T d_0  (ONLY the exact part of Hodge Laplacian)
      Full Hodge: Delta_1 = d_0^T d_0 + d_1 d_1^T

    The graph's L_1 is missing the co-exact part entirely.
    This is the KEY structural finding.
    """
    results = {
        "description": "Hodge decomposition on a 1-complex (no 2-cells)",
        "structural_finding": (
            "On a graph (simplicial 1-complex), the edge Laplacian L_1 = B^T B "
            "equals ONLY the exact (d_0^T d_0) part of the Hodge Laplacian. "
            "The co-exact part (d_1 d_1^T) requires 2-cells (faces/plaquettes) "
            "which do not exist on a pure graph. Therefore L_1 computes the "
            "LONGITUDINAL photon modes only. The TRANSVERSE (physical) photon "
            "modes live in the co-exact sector, which is absent."
        ),
        "decomposition": {},
    }

    for n_max in n_max_values:
        hodge = FockGraphHodge(n_max)
        V = hodge.n_nodes
        E = hodge.n_edges
        beta_0 = hodge.betti_0
        beta_1 = hodge.betti_1

        # Dimension of exact 1-forms = rank(B) = V - beta_0
        dim_exact = V - beta_0
        # Dimension of harmonic = beta_1
        dim_harmonic = beta_1
        # Dimension of co-exact = 0 (no 2-cells!)
        dim_coexact = 0

        # Check: dim_exact + dim_harmonic + dim_coexact = E
        total_check = dim_exact + dim_harmonic + dim_coexact

        results["decomposition"][str(n_max)] = {
            "V": V,
            "E": E,
            "beta_0": beta_0,
            "beta_1": beta_1,
            "dim_exact_1forms": dim_exact,
            "dim_harmonic_1forms": dim_harmonic,
            "dim_coexact_1forms": dim_coexact,
            "total_check": total_check,
            "total_equals_E": total_check == E,
            "note": (
                f"All {E} edge dimensions accounted for: "
                f"{dim_exact} exact + {dim_harmonic} harmonic + "
                f"{dim_coexact} co-exact = {total_check}. "
                "Co-exact = 0 because no 2-cells exist."
            ),
        }

        # Compare to continuum: on S^3, exact and co-exact are EQUAL
        # So the graph is missing exactly half the physical content!
        cont_n_max = n_max  # levels available
        total_cont_modes = sum(2 * n * (n + 2) for n in range(1, cont_n_max + 1))
        transverse_cont = sum(n * (n + 2) for n in range(1, cont_n_max + 1))
        longitudinal_cont = sum(n * (n + 2) for n in range(1, cont_n_max + 1))

        results["decomposition"][str(n_max)]["continuum_comparison"] = {
            "total_1form_modes_continuum": total_cont_modes,
            "transverse_modes_continuum": transverse_cont,
            "longitudinal_modes_continuum": longitudinal_cont,
            "graph_modes_nonzero": E - beta_1,
            "graph_captures_fraction": (E - beta_1) / total_cont_modes if total_cont_modes > 0 else 0,
        }

    return results


def part4_plaquette_hypothesis(n_max_values: List[int]) -> Dict:
    """Part 4: Test the plaquette hypothesis.

    Paper 30's SU(2) Wilson gauge construction enumerates plaquettes
    (primitive closed non-backtracking walks) on the Fock graph.
    Counts: (0, 2, 8) at length L=(4,6,8) for n_max=(2,3,4).

    If plaquettes define 2-cells (faces), we can construct:
      d_1: edge (1-cochain) -> face (2-cochain) boundary operator
      d_1 d_1^T: co-exact part of the Hodge Laplacian on edges

    Then Delta_1 = B^T B + d_1 d_1^T would be the full Hodge Laplacian.

    We test:
    1. Count plaquettes at each n_max
    2. Build the boundary operator d_1 (edges -> faces)
    3. Compute d_1 d_1^T and the full Delta_1
    4. Compare eigenvalues to continuum Hodge-1 spectrum
    """
    from geovac.su2_wilson_gauge import enumerate_plaquettes, enumerate_oriented_edges

    results = {
        "description": "Testing whether plaquettes complete the Hodge Laplacian",
        "hypothesis": (
            "Paper 30 plaquettes define 2-cells whose boundary operator d_1 "
            "provides the co-exact Laplacian d_1 d_1^T. The full Hodge "
            "Laplacian Delta_1 = B^T B + d_1 d_1^T should approximate the "
            "continuum spectrum n(n+2)."
        ),
        "results": {},
    }

    for n_max in n_max_values:
        print(f"  Analyzing plaquettes at n_max={n_max}...")

        # Build the Fock graph
        lat = GeometricLattice(max_n=n_max, topological_weights=False)
        V = lat.num_states
        adj_sparse = lat.adjacency
        adj = np.array(adj_sparse.todense(), dtype=float)

        # Get edges (undirected)
        rows, cols = adj_sparse.nonzero()
        edge_set = set()
        for r, c in zip(rows, cols):
            if r < c:
                edge_set.add((int(r), int(c)))
        edges = sorted(edge_set)
        E = len(edges)
        edge_index = {e: i for i, e in enumerate(edges)}

        # Build incidence matrix B (V x E)
        B = np.zeros((V, E), dtype=float)
        for k, (i, j) in enumerate(edges):
            B[i, k] = 1.0
            B[j, k] = -1.0

        # L_0 and L_1
        L0 = B @ B.T
        L1 = B.T @ B

        # Enumerate plaquettes at various lengths
        max_plaq_length = 8
        plaquettes = enumerate_plaquettes(adj, max_length=max_plaq_length,
                                          both_orientations=False)
        n_plaquettes = len(plaquettes)

        # Group by length
        length_counts = {}
        for p in plaquettes:
            L_p = len(p)
            length_counts[L_p] = length_counts.get(L_p, 0) + 1

        print(f"    n_max={n_max}: V={V}, E={E}, plaquettes={n_plaquettes}")
        print(f"    Length distribution: {length_counts}")

        # Build the boundary operator d_1: E -> F (edges -> faces)
        # For each plaquette (face), d_1 maps it to its boundary (a signed
        # sum of edges). d_1^T maps edges to faces.
        #
        # Convention: for face f with boundary edges e_1, e_2, ..., e_L
        # (oriented consistently with the face orientation):
        #   d_1[e_k, f] = +1 if e_k's canonical orientation agrees with face
        #   d_1[e_k, f] = -1 if opposite
        #
        # In simplicial terms, d_1 is the E x F matrix where
        # d_1[e, f] = orientation_sign(e in boundary(f))

        F = n_plaquettes
        if F == 0:
            results["results"][str(n_max)] = {
                "V": V, "E": E, "F": 0,
                "plaquette_count": 0,
                "length_distribution": length_counts,
                "verdict": "No plaquettes at this n_max; cannot test hypothesis.",
            }
            continue

        # Build d_1 (E x F boundary operator)
        d1 = np.zeros((E, F), dtype=float)

        for f_idx, plaq in enumerate(plaquettes):
            # plaq is a list of OrientedEdge objects
            for oe in plaq:
                src, tgt = oe.source, oe.target
                # Find which undirected edge this corresponds to
                if src < tgt:
                    e_key = (src, tgt)
                    sign = +1.0  # oriented edge agrees with canonical
                else:
                    e_key = (tgt, src)
                    sign = -1.0  # opposite to canonical orientation
                if e_key in edge_index:
                    e_idx = edge_index[e_key]
                    d1[e_idx, f_idx] += sign

        # Co-exact Laplacian: L1_coexact = d_1 d_1^T (E x E)
        L1_coexact = d1 @ d1.T

        # Full Hodge Laplacian: Delta_1 = L_1 + L1_coexact
        Delta_1 = L1 + L1_coexact

        # Compute spectra
        evals_L1 = np.sort(np.linalg.eigvalsh(L1))
        evals_coexact = np.sort(np.linalg.eigvalsh(L1_coexact))
        evals_Delta1 = np.sort(np.linalg.eigvalsh(Delta_1))

        # Nonzero eigenvalues
        nz_L1 = sorted(evals_L1[evals_L1 > 0.5])
        nz_coexact = sorted(evals_coexact[evals_coexact > 0.5])
        nz_Delta1 = sorted(evals_Delta1[evals_Delta1 > 0.5])

        # Continuum Hodge-1 eigenvalues for comparison
        # n(n+2) = 3, 8, 15, 24, 35, ...
        hodge1_cont = [n * (n + 2) for n in range(1, n_max + 2)]
        # Scalar eigenvalues (what L_1 should give): n^2-1 = 3, 8, 15, 24, ...
        scalar_cont = [n * n - 1 for n in range(2, n_max + 3)]

        # Compare: does Delta_1 spectrum better match Hodge-1 than L_1 alone?
        # Note: on S^3, scalar n^2-1 at level n+1 = Hodge-1 n(n+2) at level n
        # So numerically the VALUES are the same: 3, 8, 15, 24, ...
        # But DEGENERACIES differ!

        # Count multiplicities for the distinct eigenvalues
        def count_multiplicities(evals, tol=0.5):
            """Group eigenvalues into clusters and count."""
            if len(evals) == 0:
                return []
            rounded = np.round(evals).astype(int)
            unique_vals = sorted(set(rounded))
            return [(v, int(np.sum(np.abs(evals - v) < tol))) for v in unique_vals]

        mult_L1 = count_multiplicities(np.array(nz_L1))
        mult_coexact = count_multiplicities(np.array(nz_coexact))
        mult_Delta1 = count_multiplicities(np.array(nz_Delta1))

        # Expected continuum multiplicities:
        # Scalar at level n: degeneracy n^2 (on S^3 for scalar harmonics)
        # But on the graph: these get split by the truncation.
        # For the Fock graph, L_0 eigenvalue n^2-1 has degeneracy = number
        # of states at that level = n^2 MINUS boundary effects.

        # The KEY comparison: does adding co-exact (from plaquettes) change
        # the eigenvalue VALUES or just the multiplicities?
        eigenvalue_shift = []
        for ev_l1 in nz_L1[:min(5, len(nz_L1))]:
            # Find nearest Delta_1 eigenvalue
            if len(nz_Delta1) > 0:
                nearest_d1 = min(nz_Delta1, key=lambda x: abs(x - ev_l1))
                eigenvalue_shift.append({
                    "L1_eigenvalue": round(ev_l1, 4),
                    "nearest_Delta1": round(nearest_d1, 4),
                    "shift": round(nearest_d1 - ev_l1, 4),
                })

        # Rank analysis of d_1
        rank_d1 = int(np.linalg.matrix_rank(d1, tol=1e-8))

        results["results"][str(n_max)] = {
            "V": V,
            "E": E,
            "F": F,
            "plaquette_count": n_plaquettes,
            "length_distribution": {str(k): v for k, v in length_counts.items()},
            "rank_d1": rank_d1,
            "d1_shape": list(d1.shape),
            "L1_nonzero_count": len(nz_L1),
            "L1_coexact_nonzero_count": len(nz_coexact),
            "Delta1_nonzero_count": len(nz_Delta1),
            "L1_zero_count": int(np.sum(np.abs(evals_L1) < 0.5)),
            "L1_coexact_zero_count": int(np.sum(np.abs(evals_coexact) < 0.5)),
            "Delta1_zero_count": int(np.sum(np.abs(evals_Delta1) < 0.5)),
            "L1_multiplicities": [(int(v), int(m)) for v, m in mult_L1],
            "L1_coexact_multiplicities": [(int(v), int(m)) for v, m in mult_coexact],
            "Delta1_multiplicities": [(int(v), int(m)) for v, m in mult_Delta1],
            "eigenvalue_shift_samples": eigenvalue_shift,
            "scalar_continuum": scalar_cont[:6],
            "hodge1_continuum": hodge1_cont[:6],
            "L1_nonzero_eigenvalues": [round(x, 4) for x in nz_L1[:20]],
            "coexact_nonzero_eigenvalues": [round(x, 4) for x in nz_coexact[:20]],
            "Delta1_nonzero_eigenvalues": [round(x, 4) for x in nz_Delta1[:20]],
            "d1_matrix_sample": d1[:min(10, E), :min(10, F)].tolist() if F > 0 else [],
        }

    return results


def part5_spectrum_summary(results_p2: Dict, results_p4: Dict) -> Dict:
    """Part 5: Synthesis and interpretation.

    Key questions answered:
    1. Does L_0 match the scalar continuum? (Should: Paper 7)
    2. Does L_1 match the scalar continuum? (Should: SVD theorem)
    3. Does adding plaquettes (d_1 d_1^T) change the spectrum?
    4. Does Delta_1 = L_1 + d_1 d_1^T approximate Hodge-1?
    """
    summary = {
        "question_1_L0_matches_scalar": (
            "L_0 eigenvalues should match the scalar continuum n^2-1 with "
            "multiplicities modified by the finite graph truncation. "
            "This is Paper 7's fundamental result."
        ),
        "question_2_L1_matches_L0": (
            "L_1 nonzero eigenvalues match L_0 nonzero eigenvalues by the "
            "SVD theorem for incidence matrices. This is structural."
        ),
        "question_3_plaquettes_add_coexact": None,
        "question_4_Delta1_approaches_hodge1": None,
        "structural_insight": None,
    }

    # Check if plaquette results are available
    if "results" in results_p4:
        for n_max_str, data in results_p4["results"].items():
            if isinstance(data, dict) and "L1_coexact_nonzero_count" in data:
                has_coexact = data["L1_coexact_nonzero_count"] > 0
                if has_coexact:
                    summary["question_3_plaquettes_add_coexact"] = (
                        f"YES at n_max={n_max_str}: plaquettes generate "
                        f"{data['L1_coexact_nonzero_count']} nonzero co-exact "
                        f"eigenvalues from {data['F']} faces."
                    )
                    # Check if eigenvalues shifted
                    if data.get("eigenvalue_shift_samples"):
                        shifts = [s["shift"] for s in data["eigenvalue_shift_samples"]]
                        max_shift = max(abs(s) for s in shifts) if shifts else 0
                        summary["question_4_Delta1_approaches_hodge1"] = (
                            f"At n_max={n_max_str}: max eigenvalue shift from "
                            f"adding co-exact is {max_shift:.4f}. "
                            f"Delta_1 eigenvalues: {data['Delta1_nonzero_eigenvalues'][:8]}"
                        )
                else:
                    summary["question_3_plaquettes_add_coexact"] = (
                        f"At n_max={n_max_str}: {data['F']} plaquettes found "
                        f"but co-exact Laplacian has "
                        f"{data['L1_coexact_nonzero_count']} nonzero eigenvalues."
                    )

    # The structural insight
    summary["structural_insight"] = (
        "On a simplicial 1-complex (graph), L_1 = B^T B = d_0^T d_0 is ONLY "
        "the exact (longitudinal) part of the Hodge-1 Laplacian. The co-exact "
        "(transverse = physical photon) part requires d_1 d_1^T, where d_1 is "
        "the boundary operator from edges to faces (2-cells). On the continuum "
        "S^3, exact and co-exact 1-forms have IDENTICAL spectra (both n(n+2)), "
        "related by Hodge duality. So L_1 alone already has the correct "
        "EIGENVALUES (since scalar n^2-1 at level n+1 = n(n+2)), but it "
        "captures only the EXACT (longitudinal) sector, not the CO-EXACT "
        "(transverse) sector. The physical photon propagator requires the "
        "co-exact sector, which lives on 2-cells."
    )

    # Nuanced point about the scalar = Hodge-1 eigenvalue coincidence
    summary["subtlety_eigenvalue_coincidence"] = (
        "A remarkable feature of S^3 is that the scalar Laplacian eigenvalues "
        "(n^2-1 at level n in Fock convention) numerically EQUAL the Hodge-1 "
        "eigenvalues (n(n+2) at level n) after a level shift: "
        "scalar(n+1) = (n+1)^2-1 = n^2+2n = n(n+2) = Hodge-1(n). "
        "This means L_1 (which has the scalar spectrum by SVD) already has "
        "eigenvalues that are VALUES in the Hodge-1 spectrum. The difference "
        "is purely in the SECTOR ASSIGNMENT: L_1 puts these eigenvalues in "
        "the exact sector, while the physical photon needs them in the "
        "co-exact sector. On S^3 with Hodge duality, both sectors have "
        "eigenvalue n(n+2) -- so the VALUES are right, the PHYSICS is wrong "
        "(longitudinal, not transverse)."
    )

    return summary


def main():
    """Run the full analysis."""
    print("=" * 72)
    print("Fock Projection of Transverse Photon Propagator onto S^3")
    print("=" * 72)
    t0 = time.time()

    output = {
        "title": "Fock projection of flat-space transverse photon onto S^3",
        "date": "2026-05-01",
        "hypothesis": (
            "L_1 = B^T B on the Fock graph is the LONGITUDINAL (exact) "
            "photon sector only. The TRANSVERSE (co-exact, physical) photon "
            "requires 2-cells (plaquettes) from Paper 30's Wilson gauge."
        ),
    }

    # Part 1: Vector harmonics on S^3
    print("\nPart 1: Vector harmonics on S^3...")
    output["part1_vector_harmonics"] = part1_vector_harmonics(n_max=6)
    output["continuum_spectra"] = continuum_spectra(n_max=6)

    # Part 2: Compare L_1 to continuum
    print("\nPart 2: Compare L_1 to continuum...")
    output["part2_L1_comparison"] = part2_compare_L1_to_continuum([2, 3, 4])

    # Part 3: Hodge decomposition analysis
    print("\nPart 3: Hodge decomposition on 1-complex...")
    output["part3_hodge_decomposition"] = part3_hodge_decomposition_analysis([2, 3, 4])

    # Part 4: Plaquette hypothesis
    print("\nPart 4: Plaquette hypothesis test...")
    output["part4_plaquette_test"] = part4_plaquette_hypothesis([2, 3, 4])

    # Part 5: Summary
    print("\nPart 5: Synthesis...")
    output["part5_summary"] = part5_spectrum_summary(
        output["part2_L1_comparison"],
        output["part4_plaquette_test"],
    )

    # Save results
    output_path = Path(__file__).parent / "data" / "fock_photon_projection.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Convert any remaining numpy types for JSON serialization
    def convert(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        elif isinstance(obj, (np.floating,)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.integer,)):
                return int(obj)
            elif isinstance(obj, (np.floating,)):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            return super().default(obj)

    with open(output_path, "w") as f:
        json.dump(output, f, indent=2, cls=NumpyEncoder)

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.1f}s. Results saved to: {output_path}")

    # Print key findings
    print("\n" + "=" * 72)
    print("KEY FINDINGS")
    print("=" * 72)

    print("\n1. STRUCTURAL RESULT:")
    print("   L_1 = B^T B = d_0^T d_0 (exact/longitudinal ONLY)")
    print("   Co-exact (transverse/physical photon) requires d_1 d_1^T")
    print("   d_1 requires 2-cells (faces/plaquettes)")

    print("\n2. EIGENVALUE COINCIDENCE ON S^3:")
    print("   Scalar: n^2-1 at level n+1 = n^2+2n = n(n+2)")
    print("   Hodge-1: n(n+2) at level n")
    print("   Same VALUES, but L_1 assigns them to the EXACT sector")
    print("   Physical photon needs them in the CO-EXACT sector")

    if "results" in output["part4_plaquette_test"]:
        for n_str, data in output["part4_plaquette_test"]["results"].items():
            if isinstance(data, dict) and "F" in data:
                print(f"\n3. PLAQUETTE TEST at n_max={n_str}:")
                print(f"   Faces (plaquettes): {data['F']}")
                print(f"   rank(d_1): {data.get('rank_d1', 'N/A')}")
                print(f"   L_1 nonzero eigenvalues: {data.get('L1_nonzero_count', 'N/A')}")
                print(f"   Co-exact nonzero eigenvalues: {data.get('L1_coexact_nonzero_count', 'N/A')}")
                print(f"   Delta_1 nonzero eigenvalues: {data.get('Delta1_nonzero_count', 'N/A')}")
                if data.get("Delta1_nonzero_eigenvalues"):
                    print(f"   Delta_1 spectrum (first 8): {data['Delta1_nonzero_eigenvalues'][:8]}")
                if data.get("L1_nonzero_eigenvalues"):
                    print(f"   L_1 spectrum (first 8): {data['L1_nonzero_eigenvalues'][:8]}")

    return output


if __name__ == "__main__":
    main()
