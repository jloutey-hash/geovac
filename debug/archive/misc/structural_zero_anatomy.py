"""
Anatomize the broken structural zero in graph-native QED self-energy.
=====================================================================

In continuum QED on S^3, the self-energy Sigma(n_ext=0) = 0 exactly.
The proof is elementary:
  - Vertex parity requires n_int + n_ext + q_gamma to be odd.
  - With n_ext = 0, this becomes n_int + q_gamma odd.
  - The triangle inequality forces q_gamma = n_int (at minimum).
  - Combined: 2*n_int must be odd -- impossible.

On the finite GeoVac graph, this zero BREAKS.  The ground-state block
of the self-energy is [[1,1],[1,1]], not zero.

This script traces the mechanism step by step:
  1. Builds the self-energy at n_max=2, t=0 (exact sympy).
  2. Decomposes the GS self-energy into per-edge-pair contributions.
  3. Traces the full CG projection chain for each nonzero contribution.
  4. Identifies the key difference: the continuum SO(4) channel count
     W(n1=0, n2, q) = 0 for all (n2, q), while the graph CG projection
     opens couplings that the SO(4) vector harmonic structure forbids.
  5. Saves results to debug/data/structural_zero_anatomy.json.

All arithmetic is exact sympy.  No existing files are modified.

References
----------
- Paper 28 Theorem 4 (self-energy structural zero, continuum)
- GN-5 memo (broken structural zero on graph)
- geovac/graph_qed_self_energy.py (compute_self_energy)
- geovac/graph_qed_vertex.py (build_projection_matrix, build_vertex_tensor)
- geovac/graph_qed_photon.py (photon propagator)
- geovac/qed_vertex.py (continuum SO(4) selection rules)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import sympy as sp
from sympy import Integer, Matrix, Rational, zeros as sp_zeros

# Add project root to path
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from geovac.dirac_matrix_elements import DiracLabel, iter_dirac_labels, kappa_to_l
from geovac.graph_qed_vertex import (
    build_projection_matrix,
    build_vertex_tensor,
    vertex_tensor_to_matrices,
)
from geovac.graph_qed_photon import (
    build_fock_graph,
    compute_photon_propagator,
)
from geovac.graph_qed_self_energy import (
    compute_self_energy,
    _ground_state_indices,
)
from geovac.qed_vertex import (
    _vertex_allowed,
    so4_channel_count,
)


def _dirac_label_str(lab: DiracLabel) -> str:
    """Human-readable Dirac label string."""
    l = kappa_to_l(lab.kappa)
    j = lab.j
    return (f"|n={lab.n_fock}, kappa={lab.kappa}, m_j={lab.m_j}> "
            f"(l={l}, j={j})")


def _fock_label_str(states: list, idx: int) -> str:
    """Human-readable Fock node label string."""
    n, l, m = states[idx]
    return f"|n={n}, l={l}, m={m}>"


def _edge_label_str(edges: list, states: list, idx: int) -> str:
    """Human-readable edge label with node labels."""
    v1, v2 = edges[idx]
    s1 = states[v1]
    s2 = states[v2]
    # Classify edge type
    if s1[0] != s2[0]:
        etype = "T+ (radial)"
    elif s1[2] != s2[2]:
        etype = "L+ (angular)"
    else:
        etype = "??"
    return (f"e{idx}: {_fock_label_str(states, v1)} -- "
            f"{_fock_label_str(states, v2)}  [{etype}]")


def main() -> Dict[str, Any]:
    """Run the full structural zero anatomy."""

    n_max = 2
    t = Rational(0)
    results: Dict[str, Any] = {"n_max": n_max, "t": str(t)}

    print("=" * 72)
    print("STRUCTURAL ZERO ANATOMY: graph-native QED self-energy at n_max=2")
    print("=" * 72)
    print()

    # ------------------------------------------------------------------
    # PART 0: Inventory of the graph and Dirac states
    # ------------------------------------------------------------------
    print("PART 0: Graph and state inventory")
    print("-" * 40)

    fock_data = build_fock_graph(n_max)
    states = fock_data.states
    edges = fock_data.edges
    V_fock = fock_data.V
    E_fock = fock_data.E

    print(f"Fock graph: V = {V_fock} nodes, E = {E_fock} edges")
    print(f"Betti: beta_0 = {fock_data.beta_0}, beta_1 = {fock_data.beta_1}")
    print()

    print("Fock nodes:")
    for i, s in enumerate(states):
        print(f"  v{i}: |n={s[0]}, l={s[1]}, m={s[2]}>")
    print()

    print("Fock edges:")
    for i in range(E_fock):
        print(f"  {_edge_label_str(edges, states, i)}")
    print()

    dirac_labels = list(iter_dirac_labels(n_max))
    N_dirac = len(dirac_labels)
    print(f"Dirac states: N_dirac = {N_dirac}")
    for i, lab in enumerate(dirac_labels):
        print(f"  d{i}: {_dirac_label_str(lab)}")
    print()

    gs_indices = _ground_state_indices(dirac_labels)
    print(f"Ground state indices: {gs_indices}")
    for gi in gs_indices:
        print(f"  d{gi}: {_dirac_label_str(dirac_labels[gi])}")
    print()

    results["fock_nodes"] = [f"({s[0]},{s[1]},{s[2]})" for s in states]
    results["fock_edges"] = [
        f"v{e[0]}({states[e[0]]})-v{e[1]}({states[e[1]]})" for e in edges
    ]
    results["dirac_labels"] = [
        {"n_fock": lab.n_fock, "kappa": lab.kappa,
         "two_m_j": lab.two_m_j, "l": kappa_to_l(lab.kappa),
         "j": str(lab.j)}
        for lab in dirac_labels
    ]
    results["ground_state_indices"] = gs_indices
    results["V_fock"] = V_fock
    results["E_fock"] = E_fock
    results["N_dirac"] = N_dirac

    # ------------------------------------------------------------------
    # PART 1: Projection matrix P and vertex tensor V
    # ------------------------------------------------------------------
    print("PART 1: CG projection matrix P (Dirac -> Fock)")
    print("-" * 40)

    P, _, fock_states = build_projection_matrix(n_max)
    print(f"P shape: {P.rows} x {P.cols}  (N_dirac x V_fock)")
    print()

    # Show the projection for each Dirac state
    proj_detail = []
    for a in range(N_dirac):
        lab = dirac_labels[a]
        nonzero_fock = []
        for v in range(V_fock):
            pval = P[a, v]
            if pval != 0:
                nonzero_fock.append((v, states[v], str(pval), float(pval)))
        detail = {
            "dirac_idx": a,
            "dirac_label": _dirac_label_str(lab),
            "projections": [
                {"fock_idx": v, "fock_state": f"({s[0]},{s[1]},{s[2]})",
                 "cg_coeff": cg, "cg_float": cf}
                for v, s, cg, cf in nonzero_fock
            ]
        }
        proj_detail.append(detail)
        proj_strs = [f"  v{v} {_fock_label_str(states, v)}: P = {cg} = {cf:.6f}"
                     for v, s, cg, cf in nonzero_fock]
        print(f"  d{a} {_dirac_label_str(lab)}")
        for ps in proj_strs:
            print(ps)
    print()

    # Key observation: GS Dirac states
    print("KEY OBSERVATION: Ground-state projection")
    print("  Both GS Dirac states (m_j = +1/2 and -1/2) have l=0,")
    print("  so the only m_l is 0, and the CG coefficient is 1.")
    print("  They project identically onto v0 = |1,0,0> with weight 1.")
    for gi in gs_indices:
        lab = dirac_labels[gi]
        pval = P[gi, 0]
        print(f"  d{gi}: P[{gi}, 0] = {pval}")
    print()

    results["projection_detail"] = proj_detail

    # ------------------------------------------------------------------
    # PART 2: Build self-energy and extract GS block
    # ------------------------------------------------------------------
    print("PART 2: Self-energy Sigma at t=0")
    print("-" * 40)

    se_result = compute_self_energy(n_max, t=t, exact=True)

    print(f"Sigma shape: {se_result.N_dirac} x {se_result.N_dirac}")
    print(f"Sigma trace: {se_result.trace}")
    print(f"Sigma is rational: {se_result.is_rational}")
    print(f"Sigma is pi-free: {se_result.is_pi_free}")
    print(f"Sigma is Hermitian: {se_result.is_hermitian}")
    print(f"GS block is zero: {se_result.ground_state_zero}")
    print()

    gs_block = se_result.ground_state_block
    print(f"Ground-state block Sigma[GS, GS] (2x2):")
    for ii in range(2):
        row_strs = []
        for jj in range(2):
            row_strs.append(str(gs_block[ii, jj]))
        print(f"  [{', '.join(row_strs)}]")
    print()

    print(f"BROKEN STRUCTURAL ZERO: the GS block is [[1,1],[1,1]], NOT zero.")
    print()

    results["sigma_trace"] = se_result.trace
    results["sigma_rational"] = se_result.is_rational
    results["sigma_pi_free"] = se_result.is_pi_free
    results["sigma_hermitian"] = se_result.is_hermitian
    results["gs_block_zero"] = se_result.ground_state_zero
    results["gs_block"] = [[str(gs_block[i, j]) for j in range(2)]
                           for i in range(2)]

    # ------------------------------------------------------------------
    # PART 3: Per-edge-pair decomposition of GS block
    # ------------------------------------------------------------------
    print("PART 3: Per-edge-pair decomposition of GS self-energy")
    print("-" * 40)
    print()
    print("Sigma[a, b] = sum_{e1, e2} G_gamma[e1, e2] * (V_e1 . V_e2^T)[a, b]")
    print()
    print("We decompose the GS 2x2 block into contributions from each")
    print("(e1, e2) edge pair.")
    print()

    # Rebuild vertex matrices
    entries, _, _, _ = build_vertex_tensor(
        n_max, P=P, dirac_labels=dirac_labels, fock_data=fock_data
    )
    V_mats = vertex_tensor_to_matrices(entries, N_dirac, E_fock)

    # Photon propagator
    photon_data = compute_photon_propagator(n_max, exact=True)
    G_gamma = photon_data.G_gamma

    print("Photon propagator G_gamma = L_1^+ (Moore-Penrose pseudoinverse):")
    for e1 in range(E_fock):
        row_strs = []
        for e2 in range(E_fock):
            row_strs.append(str(G_gamma[e1, e2]))
        print(f"  [{', '.join(row_strs)}]")
    print()

    results["G_gamma"] = [[str(G_gamma[i, j]) for j in range(E_fock)]
                          for i in range(E_fock)]

    # Decompose GS block into per-edge-pair contributions
    edge_pair_data = []
    total_gs_block = sp_zeros(2, 2)

    for e1 in range(E_fock):
        for e2 in range(E_fock):
            g_ee = G_gamma[e1, e2]
            if g_ee == 0:
                continue

            # Contribution to GS block: g_ee * (V_e1 . V_e2^T)[gs, gs]
            product = V_mats[e1] * V_mats[e2].T
            gs_contrib = sp_zeros(2, 2)
            for ii, gi in enumerate(gs_indices):
                for jj, gj in enumerate(gs_indices):
                    gs_contrib[ii, jj] = g_ee * product[gi, gj]
                    gs_contrib[ii, jj] = sp.nsimplify(
                        sp.expand(gs_contrib[ii, jj]), rational=False
                    )

            is_zero = all(gs_contrib[i, j] == 0 for i in range(2)
                          for j in range(2))

            if not is_zero:
                total_gs_block = total_gs_block + gs_contrib

                # Also extract the inner product (V_e1 . V_e2^T)[gs, gs]
                # BEFORE multiplying by g_ee
                inner_gs = sp_zeros(2, 2)
                for ii, gi in enumerate(gs_indices):
                    for jj, gj in enumerate(gs_indices):
                        inner_gs[ii, jj] = sp.nsimplify(
                            sp.expand(product[gi, gj]), rational=False
                        )

                pair_info = {
                    "e1": e1,
                    "e2": e2,
                    "e1_label": _edge_label_str(edges, states, e1),
                    "e2_label": _edge_label_str(edges, states, e2),
                    "G_gamma_e1_e2": str(g_ee),
                    "V_e1_V_e2T_gs_block": [
                        [str(inner_gs[i, j]) for j in range(2)]
                        for i in range(2)
                    ],
                    "contribution": [
                        [str(gs_contrib[i, j]) for j in range(2)]
                        for i in range(2)
                    ],
                }
                edge_pair_data.append(pair_info)

                print(f"  Edge pair (e{e1}, e{e2}):")
                print(f"    e{e1}: {_edge_label_str(edges, states, e1)}")
                print(f"    e{e2}: {_edge_label_str(edges, states, e2)}")
                print(f"    G_gamma[{e1},{e2}] = {g_ee}")
                print(f"    (V_e1 . V_e2^T)[GS,GS] =")
                for ii in range(2):
                    row_strs = [str(inner_gs[ii, jj]) for jj in range(2)]
                    print(f"      [{', '.join(row_strs)}]")
                print(f"    Weighted contribution G * (V.V^T)[GS,GS] =")
                for ii in range(2):
                    row_strs = [str(gs_contrib[ii, jj]) for jj in range(2)]
                    print(f"      [{', '.join(row_strs)}]")
                print()

    # Verify total
    total_check = sp_zeros(2, 2)
    for ii in range(2):
        for jj in range(2):
            total_check[ii, jj] = sp.nsimplify(
                sp.expand(total_gs_block[ii, jj]), rational=False
            )
    print(f"Total GS block from decomposition:")
    for ii in range(2):
        row_strs = [str(total_check[ii, jj]) for jj in range(2)]
        print(f"  [{', '.join(row_strs)}]")

    # Cross-check against direct computation
    match = all(
        total_check[i, j] == gs_block[i, j]
        for i in range(2) for j in range(2)
    )
    print(f"Matches direct Sigma[GS,GS]: {match}")
    print()

    results["edge_pair_contributions"] = edge_pair_data
    results["total_gs_block_from_decomposition"] = [
        [str(total_check[i, j]) for j in range(2)] for i in range(2)
    ]
    results["decomposition_matches_direct"] = match

    # ------------------------------------------------------------------
    # PART 4: Trace the CG coupling chain
    # ------------------------------------------------------------------
    print("PART 4: CG coupling chain for each nonzero contribution")
    print("-" * 40)
    print()
    print("For each contributing edge pair (e1, e2), we trace:")
    print("  GS Dirac[a] --P[a,v1]--> Fock[v1] --edge e1--> Fock[v2]")
    print("  Fock[v2] --P[c,v2]--> intermediate Dirac[c]")
    print("  intermediate Dirac[c] --P[c,v3]--> Fock[v3] --edge e2--> Fock[v4]")
    print("  Fock[v4] --P[b,v4]--> GS Dirac[b]")
    print()
    print("The vertex V_e[a,c] = sum_{v1,v2 on edge e} P[a,v1]*P[c,v2]")
    print("  + P[a,v2]*P[c,v1]  (symmetrized, undirected edge).")
    print()

    chain_data = []

    for pair_info in edge_pair_data:
        e1 = pair_info["e1"]
        e2 = pair_info["e2"]
        v1_e1, v2_e1 = edges[e1]
        v1_e2, v2_e2 = edges[e2]

        print(f"--- Edge pair (e{e1}, e{e2}) ---")
        print(f"  e{e1}: v{v1_e1} {_fock_label_str(states, v1_e1)} -- "
              f"v{v2_e1} {_fock_label_str(states, v2_e1)}")
        print(f"  e{e2}: v{v1_e2} {_fock_label_str(states, v1_e2)} -- "
              f"v{v2_e2} {_fock_label_str(states, v2_e2)}")
        print()

        # For each GS pair (a, b), trace V_e1[a,c] * V_e2[c,b]
        # summed over intermediate c
        chain_entries = []

        for ii, gi in enumerate(gs_indices):
            for jj, gj in enumerate(gs_indices):
                lab_a = dirac_labels[gi]
                lab_b = dirac_labels[gj]

                print(f"  Sigma[GS_{ii}, GS_{jj}]: "
                      f"d{gi} {_dirac_label_str(lab_a)} -> "
                      f"d{gj} {_dirac_label_str(lab_b)}")

                # V_e1[a, c] for GS index a = gi
                # V_e2[c, b] for GS index b = gj (note: from V_e2^T, so
                #   (V_e2^T)[c, b] = V_e2[b, c])
                # Actually Sigma = V_e1 . V_e2^T, so
                # Sigma[a,b] = sum_c V_e1[a,c] * V_e2^T[c,b]
                #            = sum_c V_e1[a,c] * V_e2[b,c]

                intermediate_sum = sp.Integer(0)
                inter_details = []

                for c in range(N_dirac):
                    ve1_ac = V_mats[e1][gi, c]
                    ve2_bc = V_mats[e2][gj, c]

                    if ve1_ac == 0 or ve2_bc == 0:
                        continue

                    product_c = ve1_ac * ve2_bc
                    product_c = sp.nsimplify(sp.expand(product_c),
                                            rational=False)
                    intermediate_sum += product_c

                    lab_c = dirac_labels[c]

                    # Decompose V_e1[a,c] into CG chain
                    # V_e1[a,c] = P[a,v1]*P[c,v2] + P[a,v2]*P[c,v1]
                    p_a_v1 = P[gi, v1_e1]
                    p_a_v2 = P[gi, v2_e1]
                    p_c_v1 = P[c, v1_e1]
                    p_c_v2 = P[c, v2_e1]
                    term1_e1 = p_a_v1 * p_c_v2
                    term2_e1 = p_a_v2 * p_c_v1

                    # Similarly for V_e2[b,c]
                    p_b_v1e2 = P[gj, v1_e2]
                    p_b_v2e2 = P[gj, v2_e2]
                    p_c_v1e2 = P[c, v1_e2]
                    p_c_v2e2 = P[c, v2_e2]
                    term1_e2 = p_b_v1e2 * p_c_v2e2
                    term2_e2 = p_b_v2e2 * p_c_v1e2

                    detail = {
                        "intermediate_c": c,
                        "dirac_c": _dirac_label_str(lab_c),
                        "V_e1_ac": str(ve1_ac),
                        "V_e2_bc": str(ve2_bc),
                        "product": str(product_c),
                        "cg_chain_e1": {
                            "P_a_v1": str(p_a_v1),
                            "P_c_v2": str(p_c_v2),
                            "P_a_v2": str(p_a_v2),
                            "P_c_v1": str(p_c_v1),
                            "term1": str(sp.nsimplify(sp.expand(term1_e1))),
                            "term2": str(sp.nsimplify(sp.expand(term2_e1))),
                        },
                        "cg_chain_e2": {
                            "P_b_v1": str(p_b_v1e2),
                            "P_c_v2": str(p_c_v2e2),
                            "P_b_v2": str(p_b_v2e2),
                            "P_c_v1": str(p_c_v1e2),
                            "term1": str(sp.nsimplify(sp.expand(term1_e2))),
                            "term2": str(sp.nsimplify(sp.expand(term2_e2))),
                        },
                    }
                    inter_details.append(detail)

                    print(f"    c = d{c} {_dirac_label_str(lab_c)}")
                    print(f"      V_e{e1}[a,c] = {ve1_ac}")
                    print(f"        = P[a,v{v1_e1}]*P[c,v{v2_e1}] + "
                          f"P[a,v{v2_e1}]*P[c,v{v1_e1}]")
                    print(f"        = {p_a_v1}*{p_c_v2} + {p_a_v2}*{p_c_v1}")
                    print(f"        = {sp.nsimplify(sp.expand(term1_e1))} + "
                          f"{sp.nsimplify(sp.expand(term2_e1))}")
                    print(f"      V_e{e2}[b,c] = {ve2_bc}")
                    print(f"        = P[b,v{v1_e2}]*P[c,v{v2_e2}] + "
                          f"P[b,v{v2_e2}]*P[c,v{v1_e2}]")
                    print(f"        = {p_b_v1e2}*{p_c_v2e2} + "
                          f"{p_b_v2e2}*{p_c_v1e2}")
                    print(f"        = {sp.nsimplify(sp.expand(term1_e2))} + "
                          f"{sp.nsimplify(sp.expand(term2_e2))}")
                    print(f"      V_e1[a,c] * V_e2[b,c] = {product_c}")
                    print()

                intermediate_sum = sp.nsimplify(sp.expand(intermediate_sum),
                                                rational=False)
                print(f"    sum_c V_e{e1}[a,c] * V_e{e2}[b,c] = "
                      f"{intermediate_sum}")
                print()

                chain_entries.append({
                    "gs_row": ii,
                    "gs_col": jj,
                    "dirac_a": gi,
                    "dirac_b": gj,
                    "inner_product": str(intermediate_sum),
                    "intermediates": inter_details,
                })

        chain_data.append({
            "e1": e1,
            "e2": e2,
            "chains": chain_entries,
        })

    results["cg_coupling_chains"] = chain_data

    # ------------------------------------------------------------------
    # PART 5: Continuum vs graph comparison
    # ------------------------------------------------------------------
    print()
    print("PART 5: Continuum vs graph -- the key difference")
    print("-" * 40)
    print()

    # Continuum: what does SO(4) selection say about n_ext = 0?
    # In Camporesi-Higuchi convention, n_ext = 0 is the GS (n_fock = 1).
    n_ext_ch = 0  # CH convention for ground state

    print("CONTINUUM SO(4) SELECTION RULES at n_ext = 0 (CH convention):")
    print()
    print("  Vertex requires: n1 + n2 + q_gamma is ODD")
    print("  Triangle: |n1 - n2| <= q_gamma <= n1 + n2")
    print("  Photon: q_gamma >= 1")
    print()
    print("  Self-energy has n1 = n_ext = 0.")
    print("  Then: 0 + n_int + q must be odd, i.e. n_int + q is odd.")
    print("  Triangle: |0 - n_int| <= q <= 0 + n_int, so q = n_int.")
    print("  Combined: n_int + n_int = 2*n_int must be odd -- IMPOSSIBLE.")
    print()
    print("  PROOF: Sigma(n_ext=0) = 0 in the continuum.")
    print()

    # Check all possible (n_int, q) for n_ext = 0
    continuum_table = []
    print(f"  Verification: all (n_int, q_gamma) for n_ext = 0, n_int <= {n_max}:")
    for n_int in range(0, n_max + 1):
        for q in range(1, 2 * n_max + 1):
            allowed = _vertex_allowed(n_ext_ch, n_int, q)
            W = so4_channel_count(n_ext_ch, n_int, q)
            if allowed or (abs(n_ext_ch - n_int) <= q <= n_ext_ch + n_int):
                entry = {
                    "n_int": n_int,
                    "q_gamma": q,
                    "triangle_ok": abs(n_ext_ch - n_int) <= q <= n_ext_ch + n_int,
                    "parity_ok": (n_ext_ch + n_int + q) % 2 == 1,
                    "allowed": allowed,
                    "W": W,
                }
                continuum_table.append(entry)
                status = "ALLOWED" if allowed else "FORBIDDEN"
                parity = "odd" if (n_ext_ch + n_int + q) % 2 == 1 else "even"
                print(f"    n_int={n_int}, q={q}: "
                      f"triangle={'OK' if entry['triangle_ok'] else 'NO'}, "
                      f"sum={n_ext_ch+n_int+q} ({parity}), "
                      f"{status}, W={W}")
    print()
    print("  ALL entries have W=0 or are forbidden.  No coupling to GS.")
    print()

    results["continuum_selection_rules"] = {
        "n_ext_ch": n_ext_ch,
        "table": continuum_table,
        "all_forbidden": all(not e["allowed"] for e in continuum_table),
    }

    # Graph: what does CG projection allow?
    print("GRAPH CG PROJECTION at GS (n_fock=1, kappa=-1):")
    print()
    print("  The graph vertex V[a,b,e] = P[a,v1]*P[b,v2] + P[a,v2]*P[b,v1]")
    print("  uses SCALAR CG coefficients, not SO(4) vector harmonics.")
    print()
    print("  The CG projection is:")
    print("    |j, m_j> = sum_{m_s} <l, m_l, 1/2, m_s | j, m_j> |l, m_l> |m_s>")
    print()
    print("  For the GS (l=0, j=1/2), m_l = 0 always, and CG = 1.")
    print("  So both GS Dirac states project with weight 1 onto the")
    print("  same Fock node v0 = |1,0,0>.")
    print()
    print("  The vertex then couples GS to any Dirac state c that has")
    print("  nonzero projection onto v1 = |2,0,0> (via the radial edge")
    print("  e0: v0--v1).  This is ANY Dirac state with n_fock=2, l=0,")
    print("  i.e. (n=2, kappa=-1, m_j=+/-1/2).")
    print()
    print("  In the continuum, this coupling (n_ext=0, n_int=1, q=1)")
    print("  is FORBIDDEN by vertex parity (0+1+1=2, EVEN).")
    print("  On the graph, CG projection BYPASSES the parity rule entirely.")
    print()

    # Show which Dirac states couple to GS through each edge
    print("  Explicit couplings through each edge:")
    gs_couplings = []
    for e_idx in range(E_fock):
        v1, v2 = edges[e_idx]
        print(f"    Edge e{e_idx}: v{v1} {_fock_label_str(states, v1)} -- "
              f"v{v2} {_fock_label_str(states, v2)}")

        for gi in gs_indices:
            p_gs_v1 = P[gi, v1]
            p_gs_v2 = P[gi, v2]

            for c in range(N_dirac):
                p_c_v1 = P[c, v1]
                p_c_v2 = P[c, v2]
                ve = p_gs_v1 * p_c_v2 + p_gs_v2 * p_c_v1
                ve = sp.nsimplify(sp.expand(ve), rational=False)
                if ve != 0:
                    lab_c = dirac_labels[c]
                    coupling_info = {
                        "edge": e_idx,
                        "gs_dirac": gi,
                        "intermediate_dirac": c,
                        "intermediate_label": _dirac_label_str(lab_c),
                        "V_value": str(ve),
                        "n_int_ch": lab_c.n_fock - 1,  # CH convention
                    }
                    gs_couplings.append(coupling_info)
                    print(f"      V_e{e_idx}[d{gi}, d{c}] = {ve}")
                    print(f"        d{gi} = GS, "
                          f"d{c} = {_dirac_label_str(lab_c)}")
                    n_int_ch = lab_c.n_fock - 1
                    print(f"        In continuum: n_ext=0, n_int={n_int_ch}, "
                          f"q=1 (edge quantum number)")
                    parity_sum = 0 + n_int_ch + 1
                    print(f"        Continuum parity check: "
                          f"0 + {n_int_ch} + 1 = {parity_sum} "
                          f"({'ODD' if parity_sum % 2 == 1 else 'EVEN'})")
                    W_cont = so4_channel_count(0, n_int_ch, 1)
                    print(f"        Continuum W(0, {n_int_ch}, 1) = {W_cont}")
                    print()

    results["graph_gs_couplings"] = gs_couplings

    # ------------------------------------------------------------------
    # PART 6: The story -- synthesis
    # ------------------------------------------------------------------
    print()
    print("=" * 72)
    print("THE STORY: Why the structural zero breaks on the graph")
    print("=" * 72)
    print()
    print("1. In continuum QED on S^3, the self-energy Sigma(n_ext=0) = 0")
    print("   because vertex parity (n1+n2+q odd) forces 2*n_int to be odd,")
    print("   which is impossible.  This is Paper 28 Theorem 4.")
    print()
    print("2. On the finite GeoVac graph, the vertex coupling is constructed")
    print("   via CG projection: V[a,b,e] = sum P[a,v1]*P[b,v2] over the")
    print("   edge endpoints.  This is a SCALAR projection -- it decomposes")
    print("   each Dirac state into orbital content using Clebsch-Gordan")
    print("   coefficients, but does NOT enforce the SO(4) VECTOR HARMONIC")
    print("   channel structure that gives rise to the parity selection rule.")
    print()
    print("3. Concretely: the GS states (n_fock=1, kappa=-1, m_j=+/-1/2)")
    print("   both project with CG = 1 onto the single Fock node |1,0,0>.")
    print("   The radial edge e0: |1,0,0> -- |2,0,0> then couples GS")
    print("   to any Dirac state with n_fock=2, l=0 (i.e. kappa=-1).")
    print("   This coupling is nonzero on the graph (V_e0 = 1).")
    print()
    print("4. In the continuum, this SAME coupling (n_ext=0, n_int=1, q=1)")
    print("   has W(0, 1, 1) = 0.  The reason: the gamma^mu vertex couples")
    print("   positive-chirality (1, 0) to negative-chirality (0, 1) via")
    print("   vector harmonics V_A = (1, 0) or V_B = (0, 1) at q=1.  The")
    print("   SU(2)_L x SU(2)_R triangle inequalities BOTH fail:")
    print("     V_A: L-triangle(1/2, 1, 1/2) requires |1/2-1| <= 1/2 <= 3/2")
    print("          but 1/2 = 1/2, so this works for L; but")
    print("          R-triangle(0, 0, 1/2) requires |0-0| <= 1/2 <= 0,")
    print("          which fails (1/2 > 0).")
    print("     V_B: L-triangle(1/2, 0, 1/2) requires 1/2 <= 1/2 <= 1/2,")
    print("          OK; R-triangle(0, 1, 1/2) requires |0-1| <= 1/2 <= 1,")
    print("          which fails (1 > 1/2).")
    print("   Both channels give W = 0.  The vector harmonic structure of")
    print("   gamma^mu is what enforces the zero.")
    print()
    print("5. The graph CG projection is a SCALAR coupling (projecting")
    print("   spinor content onto orbital labels via standard CG coefficients)")
    print("   and has no mechanism to enforce the vector harmonic W = 0.")
    print("   The 'missing ingredient' is the SU(2)_L x SU(2)_R double-")
    print("   triangle structure that arises from gamma^mu's transformation")
    print("   properties under SO(4).")
    print()
    print("6. CONSEQUENCE: the graph self-energy GS block is [[1,1],[1,1]]")
    print("   (rank 1, eigenvalues 0 and 2).  This is a FINITE-GRAPH")
    print("   artifact that should vanish in the continuum limit as the")
    print("   graph's scalar CG projection recovers the full SO(4) vector")
    print("   harmonic structure.")
    print()

    results["narrative"] = {
        "continuum_mechanism": (
            "Vertex parity (n1+n2+q odd) forces 2*n_int odd at n_ext=0, "
            "which is impossible. W(0, n_int, q) = 0 for all (n_int, q)."
        ),
        "graph_mechanism": (
            "CG projection is scalar: decomposes Dirac spinors into "
            "orbital content via <l,m_l,1/2,m_s|j,m_j>, but does not "
            "enforce the SO(4) SU(2)_L x SU(2)_R double-triangle "
            "condition on the vector harmonic components of gamma^mu. "
            "The coupling V_e[GS, c] is nonzero whenever c has orbital "
            "content at an adjacent Fock node, regardless of W."
        ),
        "key_difference": (
            "The continuum zero is enforced by the vector harmonic "
            "channel count W(0, n_int, q) = 0, which requires the "
            "SU(2)_L x SU(2)_R double-triangle condition. The graph "
            "vertex uses a scalar CG projection that bypasses this "
            "structure entirely."
        ),
        "gs_block": "[[1, 1], [1, 1]]",
        "gs_eigenvalues": "[0, 2]",
    }

    # ------------------------------------------------------------------
    # Save results
    # ------------------------------------------------------------------
    out_path = Path(__file__).resolve().parent / "data" / "structural_zero_anatomy.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"Results saved to {out_path}")

    return results


if __name__ == "__main__":
    main()
