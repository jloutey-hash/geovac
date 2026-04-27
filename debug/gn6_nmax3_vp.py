"""
GN-6: Vacuum polarization at n_max=3 on the GeoVac finite graph.

Extends the n_max=2 VP computation (3x3 matrix, all rational) to n_max=3
(13x13 matrix). Tests the pi-free certificate and structural analysis.

Key questions:
(A) Graph dimensions at n_max=3
(B) VP matrix Pi[e1,e2] as finite trace
(C) Pi-free certificate: are all entries rational?
(D) Eigenvalues, block structure, scaling
(E) Per-edge VP convergence
"""

import json
import time
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from sympy import Rational, Integer, Matrix

# Add project root
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from geovac.graph_qed_vertex import (
    build_projection_matrix,
    build_vertex_tensor,
    vertex_tensor_to_matrices,
    compute_vacuum_polarization,
    _check_pi_free_matrix,
    _check_rational_matrix,
)
from geovac.graph_qed_photon import build_fock_graph
from geovac.graph_qed_propagator import DiracGraphOperator, electron_propagator
from geovac.dirac_matrix_elements import iter_dirac_labels


def main():
    results = {}

    # =====================================================================
    # (A) Graph characterization at n_max=3
    # =====================================================================
    print("=" * 60)
    print("(A) Graph characterization")
    print("=" * 60)

    for nm in [2, 3]:
        d_labels = list(iter_dirac_labels(nm))
        fg = build_fock_graph(nm)
        P, dlabs, fstates = build_projection_matrix(nm)

        info = {
            "n_max": nm,
            "N_dirac": len(d_labels),
            "V_fock": fg.V,
            "E_fock": fg.E,
            "beta_0": fg.beta_0,
            "beta_1": fg.beta_1,
            "P_shape": [P.rows, P.cols],
        }
        results[f"graph_n{nm}"] = info
        print(f"  n_max={nm}: N_dirac={info['N_dirac']}, V_fock={info['V_fock']}, "
              f"E_fock={info['E_fock']}, beta_0={info['beta_0']}, beta_1={info['beta_1']}, "
              f"P_shape={info['P_shape']}")

    # Edge classification at n_max=3
    fg3 = build_fock_graph(3)
    print(f"\n  Edge list at n_max=3 ({fg3.E} edges):")
    edge_types = []
    for e_idx, (v1, v2) in enumerate(fg3.edges):
        s1 = fg3.states[v1]
        s2 = fg3.states[v2]
        dn = s2[0] - s1[0]
        dl = s2[1] - s1[1]
        dm = s2[2] - s1[2]
        etype = "T" if dn != 0 else "L"
        edge_types.append({
            "index": e_idx,
            "nodes": [list(s1), list(s2)],
            "type": etype,
            "dn": dn, "dl": dl, "dm": dm,
        })
        print(f"    e{e_idx}: {s1} -> {s2}  [{etype}]  dn={dn}, dl={dl}, dm={dm}")
    results["edge_types_n3"] = edge_types

    # =====================================================================
    # (B) VP computation at n_max=3
    # =====================================================================
    print("\n" + "=" * 60)
    print("(B) VP computation at n_max=3 (exact sympy)")
    print("=" * 60)

    # First, verify n_max=2 reference
    print("\n  Reference: n_max=2 ...")
    t0 = time.time()
    vp2 = compute_vacuum_polarization(2, t=Rational(0), exact=True)
    t2 = time.time() - t0
    Pi2 = vp2['Pi']
    print(f"    Time: {t2:.2f}s")
    print(f"    Pi shape: {Pi2.shape}")
    print(f"    Trace: {Pi2.trace()}")
    print(f"    Eigenvalues: {sorted(Pi2.eigenvals().keys())}")
    results["vp_n2"] = {
        "trace": str(Pi2.trace()),
        "trace_float": float(Pi2.trace()),
        "eigenvalues": [str(ev) for ev in sorted(Pi2.eigenvals().keys())],
        "eigenvalue_multiplicities": {str(ev): mult for ev, mult in Pi2.eigenvals().items()},
        "time_s": t2,
    }

    # Now n_max=3
    print(f"\n  Computing VP at n_max=3 (28 Dirac states, 13 edges) ...")
    print("    Building projection matrix ...")
    t0 = time.time()
    P3, dlabs3, fstates3 = build_projection_matrix(3)
    t_proj = time.time() - t0
    P_nnz = sum(1 for e in P3 if e != 0)
    print(f"    P: {P3.rows}x{P3.cols}, {P_nnz} nonzero entries ({t_proj:.2f}s)")

    print("    Building vertex tensor ...")
    t0 = time.time()
    entries3, N_dirac3, V_fock3, E_fock3 = build_vertex_tensor(
        3, P=P3, dirac_labels=dlabs3, fock_data=fg3
    )
    t_vert = time.time() - t0
    print(f"    V: {len(entries3)} nonzero of {N_dirac3*N_dirac3*E_fock3} total ({t_vert:.2f}s)")

    print("    Building vertex matrices ...")
    V_mats3 = vertex_tensor_to_matrices(entries3, N_dirac3, E_fock3)

    print("    Building electron propagator ...")
    op3 = DiracGraphOperator(n_max=3, t=Rational(0))
    G_e3, is_rat3 = electron_propagator(op3, exact=True)
    print(f"    G_e: {G_e3.rows}x{G_e3.cols}, rational={is_rat3}")

    print("    Computing Pi[e1,e2] = Tr(V_e1^T @ G_e @ V_e2 @ G_e) ...")
    E3 = fg3.E
    Pi3 = sp.zeros(E3, E3)
    t0 = time.time()
    for e1 in range(E3):
        for e2 in range(e1, E3):  # Symmetry: Pi[e1,e2] = Pi[e2,e1]
            prod = V_mats3[e1].T * G_e3 * V_mats3[e2] * G_e3
            val = sp.nsimplify(prod.trace(), rational=False)
            Pi3[e1, e2] = val
            Pi3[e2, e1] = val
        if (e1 + 1) % 3 == 0:
            elapsed = time.time() - t0
            pct = (e1 + 1) / E3 * 100
            print(f"      Row {e1+1}/{E3} done ({pct:.0f}%, {elapsed:.1f}s)")
    t_vp = time.time() - t0
    print(f"    VP computation time: {t_vp:.2f}s")

    # =====================================================================
    # (C) Pi-free certificate at n_max=3
    # =====================================================================
    print("\n" + "=" * 60)
    print("(C) Pi-free certificate at n_max=3")
    print("=" * 60)

    pi_free = _check_pi_free_matrix(Pi3)
    print(f"  Pi is pi-free (no transcendentals): {pi_free}")

    # Check if all entries are rational
    all_rational = True
    non_rational_entries = []
    for i in range(E3):
        for j in range(E3):
            entry = Pi3[i, j]
            s = sp.nsimplify(entry, rational=True)
            if not isinstance(s, (sp.Rational, sp.Integer,
                                  sp.core.numbers.Zero,
                                  sp.core.numbers.One,
                                  sp.core.numbers.NegativeOne)):
                all_rational = False
                non_rational_entries.append((i, j, str(entry)))

    print(f"  All Pi entries rational: {all_rational}")
    if non_rational_entries:
        for i, j, val in non_rational_entries[:5]:
            print(f"    Non-rational: Pi[{i},{j}] = {val}")

    results["pi_free_cert_n3"] = {
        "pi_free": pi_free,
        "all_rational": all_rational,
        "non_rational_count": len(non_rational_entries),
        "non_rational_entries": non_rational_entries[:10],
    }

    # =====================================================================
    # (D) Structural analysis
    # =====================================================================
    print("\n" + "=" * 60)
    print("(D) Structural analysis")
    print("=" * 60)

    # Pi entries
    print("\n  Pi matrix entries (as rationals):")
    Pi3_entries = []
    for i in range(E3):
        row = []
        for j in range(E3):
            r = sp.nsimplify(Pi3[i, j], rational=True)
            row.append(str(r))
        Pi3_entries.append(row)

    # Print diagonal
    print("  Diagonal entries:")
    for i in range(E3):
        et = edge_types[i]
        print(f"    Pi[{i},{i}] = {Pi3_entries[i][i]}  [{et['type']}-edge: "
              f"{et['nodes'][0]}->{et['nodes'][1]}]")

    # Trace
    trace3 = sp.nsimplify(Pi3.trace(), rational=True)
    print(f"\n  Trace(Pi) = {trace3} = {float(trace3):.10f}")

    # Eigenvalues
    print("\n  Computing eigenvalues ...")
    t0 = time.time()
    eig_dict = Pi3.eigenvals()
    t_eig = time.time() - t0
    print(f"  Eigenvalue computation time: {t_eig:.2f}s")

    eig_list = []
    for ev, mult in sorted(eig_dict.items(), key=lambda x: float(x[0])):
        ev_simp = sp.nsimplify(ev, rational=False)
        is_rat = isinstance(ev_simp, (sp.Rational, sp.Integer,
                                       sp.core.numbers.Zero))
        eig_list.append({
            "value": str(ev_simp),
            "float": float(ev_simp),
            "multiplicity": mult,
            "is_rational": is_rat,
        })
        print(f"    lambda = {ev_simp} = {float(ev_simp):.10f}  (mult {mult}, rational={is_rat})")

    all_eigs_rational = all(e["is_rational"] for e in eig_list)
    print(f"\n  All eigenvalues rational: {all_eigs_rational}")

    # Block structure: which pairs (e1,e2) have Pi[e1,e2] != 0?
    print("\n  Block structure (nonzero off-diagonal pattern):")
    nonzero_pattern = []
    for i in range(E3):
        for j in range(i+1, E3):
            if Pi3[i, j] != 0:
                nonzero_pattern.append((i, j))
                print(f"    Pi[{i},{j}] != 0: {edge_types[i]['type']}-{edge_types[j]['type']}")

    # Group by l-sector
    print("\n  Edge l-sector classification:")
    for et in edge_types:
        s1 = et["nodes"][0]
        s2 = et["nodes"][1]
        l_set = set()
        l_set.add(s1[1])
        l_set.add(s2[1])
        et["l_sectors"] = sorted(l_set)
        print(f"    e{et['index']}: l-sectors={et['l_sectors']} "
              f"[{et['type']}: {s1}->{s2}]")

    # Scaling analysis: compare traces
    trace2 = Rational(224, 75)
    trace_ratio = sp.nsimplify(trace3 / trace2, rational=True)
    print(f"\n  Tr(Pi_3) / Tr(Pi_2) = {trace_ratio} = {float(trace_ratio):.6f}")

    # Per-edge trace
    per_edge_2 = sp.nsimplify(trace2 / 3, rational=True)
    per_edge_3 = sp.nsimplify(trace3 / 13, rational=True)
    print(f"  Tr(Pi)/E at n_max=2: {per_edge_2} = {float(per_edge_2):.10f}")
    print(f"  Tr(Pi)/E at n_max=3: {per_edge_3} = {float(per_edge_3):.10f}")

    # n_max=2 eigenvalues for embedding check
    eig2_list = sorted([Rational(32, 225), Rational(32, 45), Rational(32, 15)])
    eig3_sorted = sorted([float(e["float"]) for e in eig_list for _ in range(e["multiplicity"])])
    print(f"\n  n_max=2 eigenvalues: {eig2_list} = {[float(e) for e in eig2_list]}")
    print(f"  n_max=3 eigenvalues (sorted): {eig3_sorted}")

    # Check if n_max=2 eigenvalues appear in n_max=3
    eig3_vals = [float(e["float"]) for e in eig_list]
    embedding = []
    for ev2 in eig2_list:
        found = any(abs(float(ev2) - ev3) < 1e-10 for ev3 in eig3_vals)
        embedding.append({"n2_eigenvalue": str(ev2), "found_in_n3": found})
        if found:
            print(f"  n_max=2 eigenvalue {ev2} FOUND in n_max=3 spectrum")
        else:
            print(f"  n_max=2 eigenvalue {ev2} NOT found in n_max=3 spectrum")

    results["vp_n3"] = {
        "E_fock": E3,
        "Pi_entries": Pi3_entries,
        "trace": str(trace3),
        "trace_float": float(trace3),
        "eigenvalues": eig_list,
        "all_eigenvalues_rational": all_eigs_rational,
        "nonzero_offdiag_pairs": nonzero_pattern,
        "block_structure_count": len(nonzero_pattern),
        "trace_ratio_n3_over_n2": str(trace_ratio),
        "trace_ratio_float": float(trace_ratio),
        "per_edge_trace_n2": str(per_edge_2),
        "per_edge_trace_n2_float": float(per_edge_2),
        "per_edge_trace_n3": str(per_edge_3),
        "per_edge_trace_n3_float": float(per_edge_3),
        "embedding_n2_in_n3": embedding,
        "computation_time_s": {
            "projection": t_proj,
            "vertex": t_vert,
            "vp": t_vp,
            "eigenvalues": t_eig,
        },
        "vertex_nonzero_count": len(entries3),
        "vertex_total_possible": N_dirac3 * N_dirac3 * E_fock3,
    }

    # =====================================================================
    # (E) Scaling and convergence
    # =====================================================================
    print("\n" + "=" * 60)
    print("(E) Scaling and convergence")
    print("=" * 60)

    # Compare per-edge VP
    print(f"\n  Tr(Pi)/n_edges:")
    print(f"    n_max=2: {float(per_edge_2):.10f}  ({per_edge_2})")
    print(f"    n_max=3: {float(per_edge_3):.10f}  ({per_edge_3})")

    ratio_per_edge = float(per_edge_3) / float(per_edge_2)
    print(f"    Ratio (n3/n2): {ratio_per_edge:.6f}")

    # Also compute trace per Dirac-pair
    tr_per_pair_2 = float(trace2) / (10 * (10 - 1) / 2)
    tr_per_pair_3 = float(trace3) / (28 * (28 - 1) / 2)
    print(f"\n  Tr(Pi) / N_dirac_pairs:")
    print(f"    n_max=2: {tr_per_pair_2:.10f}")
    print(f"    n_max=3: {tr_per_pair_3:.10f}")

    results["scaling"] = {
        "per_edge_ratio_n3_n2": ratio_per_edge,
        "per_pair_n2": tr_per_pair_2,
        "per_pair_n3": tr_per_pair_3,
    }

    # =====================================================================
    # Save results
    # =====================================================================
    output_path = Path(__file__).parent / "data" / "gn6_nmax3_vp.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved: {output_path}")

    return results


if __name__ == "__main__":
    main()
