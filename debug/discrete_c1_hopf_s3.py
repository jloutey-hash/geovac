"""
Discrete first Chern class c_1 of the Hopf-S^3 Fock graph.
=========================================================

Sprint TS-E3 falsification target (2026-05-04).

Goal
----
Compute the discrete first Chern number

    c_1 = (1 / 2*pi) * sum_{plaquettes P} F(P)

on the Fock-projected S^3 Hopf graph (Paper 7, Paper 25) at n_max=2 and
n_max=3, where F(P) is the holonomy (curvature 2-form integrated over
the plaquette) computed as the signed sum of edge phases around P, and
verify that c_1 is integer-valued.

Prediction (case-exhaustion theorem of Sprint TS-E3, master Mellin
engine reading): every pi in a GeoVac observable is a Mellin-transform
of Tr(D^k * exp(-t D^2)) for k in {0, 1, 2} on a compact Riemannian
manifold. The discrete first Chern class on the *finite* Fock graph
should NOT contain pi: it is an integer (or a rational multiple of an
integer) topological invariant. Pi enters only in the continuum limit
through the M1 mechanism (Hopf-base measure normalisation Vol(S^2)/4
= pi which converts the integer instanton number into the integer-
valued c_1 against a continuous integration measure).

If c_1 is integer-valued at the tested cutoffs: master Mellin engine
reading is REINFORCED.

If c_1 is non-integer: a fourth pi-source mechanism would need to be
named (M4 candidate).

Edge phase convention
---------------------
The Hopf U(1) connection on the Fock graph (Paper 25 Sec. III.6) is
constructed from the L_+ and T_+ ladder operator matrix elements:

  <n, l, m+1 | L_+ | n, l, m>  =  hbar * sqrt(l(l+1) - m(m+1)) * exp(i phi_{L+})
  <n+1, l, m | T_+ | n, l, m>  =  amplitude(n, l) * exp(i phi_{T+})

Under the standard Condon-Shortley convention with real-valued
hydrogenic spherical harmonics, the ladder operator matrix elements
are real and positive: phi_{L+} = phi_{T+} = 0. This is the
*natural* convention and the one Paper 25 implicitly adopts (its
Observation 1 establishes that all spectral invariants of the edge
Laplacian are algebraic integers, equivalently the connection is
flat in this convention).

Three convention variants are tested in this sprint:

  A. CONDON-SHORTLEY (canonical, all phases = 0): expected c_1 = 0
     (integer trivially).

  B. ALTERNATING-SIGN (phase = pi on L_+ edges, 0 on T_+ edges, mod 2pi):
     a Z_2 twist that mimics a non-trivial bundle; expected c_1 = integer
     (because pi-valued phases mod 2pi sum to integer multiples of pi
     around any cycle, then dividing by 2*pi gives integer/2 -- which
     would be the only place a non-integer could sneak in).

  C. UNIT-FLUX MONOPOLE (each plaquette holonomy = 2*pi; configurable
     to any rational multiple of 2*pi): expected c_1 = (rational
     multiple of) the number of plaquettes.

For convention A the prediction is c_1 = 0 verbatim; for convention C
we verify the construction reproduces the prescribed integer flux.

Test n_max values
-----------------
n_max=2: V=5 nodes (s: n=1,l=0,m=0 + n=2,l=0,m=0; p: n=2,l=1,m in
  {-1,0,1}), E=5 edges, beta_1 = 1 cycle (the p-block 4-cycle through
  (2,1,m=-1) - (2,1,0) - (2,1,1) ... wait, m only goes through L+/L-
  ladders within fixed (n, l), so p-block at n_max=2 has only n=2 row,
  giving L+- chain m=-1 -- m=0 -- m=1 (path, no cycle). The graph at
  n_max=2 actually has beta_1 = 0.

n_max=3: V=14 nodes (s: 3 nodes; p: 6 nodes; d: 5 nodes), E=13 edges,
  beta_1 = 2 (Paper 25 Eq. 16: two zero modes of L_1, both in the
  p-block where the L+- m-ladder (m=-1,0,1) at n=2 is parallel to
  the L+- m-ladder at n=3, with T+- connecting same-m states across
  shells, forming two 4-cycles).

So the natural test point is n_max=3, where genuine plaquettes exist;
n_max=2 serves only as a smaller sanity check (with c_1 = 0
trivially because there are no plaquettes).

Outputs
-------
* `debug/data/discrete_c1_hopf_s3.json` -- JSON with construction,
  per-plaquette F(P), c_1 totals, and convention sensitivity.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import sympy as sp

from geovac.fock_graph_hodge import FockGraphHodge
from geovac.lattice import GeometricLattice


# ---------------------------------------------------------------------------
# 1. Cycle (plaquette) enumeration
# ---------------------------------------------------------------------------

def find_cycle_basis(
    V: int,
    edges: List[Tuple[int, int]],
) -> Tuple[List[List[int]], List[int]]:
    """
    Find a fundamental cycle basis for an undirected graph.

    Builds a spanning forest via BFS on each connected component;
    each non-tree edge produces one fundamental cycle by adding it to
    the unique tree path between its endpoints. The cycle is returned
    as an ordered list of node indices [v_0, v_1, ..., v_{k-1}] with
    v_k = v_0 closing the cycle (the closing edge is implicit).

    The number of fundamental cycles equals beta_1 = E - V + c, where
    c is the number of connected components. The cycle basis is
    canonical for the chosen spanning tree but is NOT unique; different
    spanning trees give different bases that span the same cycle space.

    Returns
    -------
    cycles : list of cycles, each as list of node indices.
    non_tree_edge_indices : list of edge indices, one per fundamental cycle,
        identifying the unique non-tree edge of each cycle (canonical
        plaquette label for assigning prescribed flux).
    """
    # Adjacency list
    adj: Dict[int, List[Tuple[int, int]]] = {v: [] for v in range(V)}
    for k, (i, j) in enumerate(edges):
        adj[i].append((j, k))
        adj[j].append((i, k))

    parent = [-1] * V        # BFS parent in the spanning forest
    parent_edge = [-1] * V   # edge index from parent to v
    visited = [False] * V
    tree_edges: set = set()

    # BFS over each component
    for root in range(V):
        if visited[root]:
            continue
        visited[root] = True
        queue = [root]
        head = 0
        while head < len(queue):
            v = queue[head]
            head += 1
            for w, ek in adj[v]:
                if not visited[w]:
                    visited[w] = True
                    parent[w] = v
                    parent_edge[w] = ek
                    tree_edges.add(ek)
                    queue.append(w)

    # Each non-tree edge gives one fundamental cycle
    cycles: List[List[int]] = []
    non_tree_edge_indices: List[int] = []
    for k, (i, j) in enumerate(edges):
        if k in tree_edges:
            continue
        non_tree_edge_indices.append(k)
        # Trace path from i and j to LCA in BFS tree
        path_i = [i]
        x = i
        while parent[x] != -1:
            x = parent[x]
            path_i.append(x)
        path_j = [j]
        x = j
        while parent[x] != -1:
            x = parent[x]
            path_j.append(x)

        seti = set(path_i)
        # Find LCA (first node on path_j that is in path_i)
        lca = None
        for x in path_j:
            if x in seti:
                lca = x
                break

        # Build cycle: i -> ... -> lca -> ... -> j -> i (closing edge)
        i_to_lca = []
        for x in path_i:
            i_to_lca.append(x)
            if x == lca:
                break
        j_to_lca = []
        for x in path_j:
            j_to_lca.append(x)
            if x == lca:
                break

        # Concatenate: i...lca, then lca...j reversed (skip duplicate lca)
        cycle = i_to_lca + list(reversed(j_to_lca[:-1]))
        cycles.append(cycle)

    return cycles, non_tree_edge_indices


# ---------------------------------------------------------------------------
# 2. Edge phases under different conventions
# ---------------------------------------------------------------------------

def edge_type(states, i: int, j: int) -> str:
    """Classify an edge as 'L' (angular, same n,l), 'T' (radial, same l,m), or 'unknown'."""
    n_i, l_i, m_i = states[i]
    n_j, l_j, m_j = states[j]
    if n_i == n_j and l_i == l_j and abs(m_i - m_j) == 1:
        return "L"
    if l_i == l_j and m_i == m_j and abs(n_i - n_j) == 1:
        return "T"
    return "unknown"


def edge_phases_condon_shortley(
    states,
    edges: List[Tuple[int, int]],
) -> List[Fraction]:
    """
    All edges have phase 0 (canonical, real ladder matrix elements).

    Returns phases as Fraction multiples of 2*pi, so phase = Fraction(0)
    means 0, phase = Fraction(1, 2) means pi, phase = Fraction(1) is
    identified with 0 mod 2*pi.
    """
    return [Fraction(0) for _ in edges]


def edge_phases_alternating(
    states,
    edges: List[Tuple[int, int]],
) -> List[Fraction]:
    """
    Z_2 twist: phase = 1/2 (i.e. pi) on L-type edges, 0 on T-type edges.

    This is a 'Z_2 connection' that flips sign on each angular hop.
    Holonomy around any cycle = (number of L-edges in cycle) mod 2,
    multiplied by pi. The discrete c_1 = sum / (2*pi) is then an
    integer/2: still an integer iff every cycle contains an even number
    of L-edges. This convention is a sanity check, not a physical
    Hopf connection.
    """
    return [
        Fraction(1, 2) if edge_type(states, i, j) == "L" else Fraction(0)
        for (i, j) in edges
    ]


def edge_phases_unit_flux_per_plaquette(
    edges: List[Tuple[int, int]],
    non_tree_edge_indices: List[int],
    flux_per_plaquette: Fraction = Fraction(1),
) -> List[Fraction]:
    """
    Construct edge phases such that the holonomy around each fundamental
    cycle in the basis equals `flux_per_plaquette` (in units of 2*pi).

    Strategy: zero phase on all spanning-tree edges; on each non-tree
    edge (which canonically labels one fundamental cycle), assign
    phase = flux_per_plaquette. Each fundamental cycle traverses
    exactly one non-tree edge (its own), so its holonomy is exactly
    `flux_per_plaquette`. This is the standard 'pick a phase for each
    independent plaquette' construction; the holonomies of the
    fundamental basis cycles are precisely the assigned values.

    This convention exists to verify the c_1 computation reproduces a
    prescribed flux as a sanity check on the algorithm.
    """
    phases = [Fraction(0) for _ in edges]
    for idx in non_tree_edge_indices:
        phases[idx] = flux_per_plaquette
    return phases


# ---------------------------------------------------------------------------
# 3. Holonomy and discrete c_1
# ---------------------------------------------------------------------------

def cycle_holonomy(
    cycle: List[int],
    edges: List[Tuple[int, int]],
    edge_phases: List[Fraction],
    reduce_mod_1: bool = False,
) -> Fraction:
    """
    Compute holonomy F(P) = sum of (signed) edge phases around the cycle,
    in units of 2*pi.

    The cycle is given as [v_0, v_1, ..., v_{k-1}] with the closing edge
    (v_{k-1}, v_0) implicit. For each consecutive pair (v_i, v_{i+1})
    we look up the undirected edge {v_i, v_{i+1}} and add its phase if
    traversed in canonical orientation (i < j), or subtract if reversed.

    For c_1 computation we DO NOT reduce mod 1 per cycle: the discrete
    first Chern number is the sum of the raw signed flux integrals,
    only optionally reduced modulo 1 at the level of the total. By
    contrast, an "individual plaquette holonomy" in U(1) is naturally
    a phase in [0, 1) (units of 2*pi), so we expose a flag.

    Parameters
    ----------
    reduce_mod_1 : bool
        If True, return F mod 1 (the U(1) holonomy element). If False
        (default for c_1 computation), return the signed sum directly.

    Returns
    -------
    F : Fraction
        Total holonomy in units of 2*pi.
    """
    # Build map from sorted edge to (index, canonical_direction).
    edge_to_idx: Dict[Tuple[int, int], int] = {}
    for k, (i, j) in enumerate(edges):
        a, b = (i, j) if i < j else (j, i)
        edge_to_idx[(a, b)] = k

    F = Fraction(0)
    n = len(cycle)
    for k in range(n):
        u = cycle[k]
        v = cycle[(k + 1) % n]
        a, b = (u, v) if u < v else (v, u)
        if (a, b) not in edge_to_idx:
            raise ValueError(
                f"Cycle edge ({u}, {v}) not in graph edge list -- cycle is invalid."
            )
        idx = edge_to_idx[(a, b)]
        sign = +1 if (u < v) else -1  # canonical orientation = (i, j) with i < j
        F += sign * edge_phases[idx]
    if reduce_mod_1:
        F = F - int(F)
        if F < 0:
            F += 1
    return F


def discrete_c1(
    cycles: List[List[int]],
    edges: List[Tuple[int, int]],
    edge_phases: List[Fraction],
    cycle_signs: List[int] | None = None,
) -> Tuple[Fraction, List[Fraction]]:
    """
    Compute discrete first Chern number c_1 = sum_P F(P) / (2*pi),
    where each F(P) is the holonomy of plaquette P (already returned
    in units of 2*pi by cycle_holonomy).

    Therefore  c_1 (as an integer) = sum_P (F(P) in units of 2*pi)
    rounded to the nearest integer; equivalently the sum of Fraction
    holonomies should land at an integer (Fraction with denominator 1)
    for any genuine U(1) connection on a finite graph.

    Parameters
    ----------
    cycle_signs : optional list of +1/-1 per cycle to test orientation
                  sensitivity. Default = all +1.

    Returns
    -------
    c1 : Fraction
        sum_P F(P) (in units of 2*pi); should be integer for valid bundles.
    F_per_plaquette : list of Fraction
        Per-plaquette holonomies (each in [0, 1) modulo 1).
    """
    if cycle_signs is None:
        cycle_signs = [1] * len(cycles)
    if len(cycle_signs) != len(cycles):
        raise ValueError("cycle_signs length must match number of cycles")

    Fs: List[Fraction] = []
    total = Fraction(0)
    for cycle, sign in zip(cycles, cycle_signs):
        F = cycle_holonomy(cycle, edges, edge_phases, reduce_mod_1=False)
        Fs.append(F)
        total += sign * F
    # NOTE: total is in units of 2*pi.
    # The "discrete c_1" is total reduced -- but we keep it as a Fraction
    # so the integrality check is exact. If total has denominator 1, c_1
    # is an integer; otherwise the computation has produced a non-integer.
    return total, Fs


# ---------------------------------------------------------------------------
# 4. Main analysis
# ---------------------------------------------------------------------------

def analyze_n_max(n_max: int) -> Dict:
    """
    Run the full discrete-c_1 falsification at given n_max.

    Returns a JSON-serializable dict.
    """
    fgh = FockGraphHodge(n_max)
    states = fgh.states
    edges = fgh.edges
    V, E = fgh.n_nodes, fgh.n_edges
    beta_0 = fgh.betti_0
    beta_1 = fgh.betti_1

    # Find a fundamental cycle basis
    cycles, non_tree_edge_indices = find_cycle_basis(V, edges)
    assert len(cycles) == beta_1, (
        f"Cycle basis size {len(cycles)} != beta_1 {beta_1}"
    )

    # Tag each edge with its physical type (L vs T)
    edge_types = [edge_type(states, i, j) for (i, j) in edges]

    # ===== Convention A: Condon-Shortley (all phases = 0) =====
    phases_A = edge_phases_condon_shortley(states, edges)
    c1_A_pos, F_A_pos = discrete_c1(cycles, edges, phases_A)
    # Orientation sensitivity check: flip all cycle signs
    c1_A_neg, F_A_neg = discrete_c1(
        cycles, edges, phases_A, cycle_signs=[-1] * len(cycles)
    )

    # ===== Convention B: Z_2 alternating (phase = pi on L-edges) =====
    phases_B = edge_phases_alternating(states, edges)
    c1_B_pos, F_B_pos = discrete_c1(cycles, edges, phases_B)
    c1_B_neg, F_B_neg = discrete_c1(
        cycles, edges, phases_B, cycle_signs=[-1] * len(cycles)
    )

    # ===== Convention C: Unit flux per plaquette (sanity check) =====
    phases_C = edge_phases_unit_flux_per_plaquette(
        edges, non_tree_edge_indices, flux_per_plaquette=Fraction(1, 4)
    )
    c1_C_pos, F_C_pos = discrete_c1(cycles, edges, phases_C)

    # ===== Convention D: Fractional flux 1/3 (algorithm stress test) =====
    # Verifies the algorithm tracks non-integer rationals exactly so we
    # would NOT miss a non-integer c_1 if it occurred under conventions A/B.
    phases_D = edge_phases_unit_flux_per_plaquette(
        edges, non_tree_edge_indices, flux_per_plaquette=Fraction(1, 3)
    )
    c1_D_pos, F_D_pos = discrete_c1(cycles, edges, phases_D)
    expected_D_magnitude = Fraction(len(cycles), 3)
    sanity_D_check = abs(c1_D_pos) == expected_D_magnitude

    # Verbose per-plaquette diagnostics
    plaquette_info: List[Dict] = []
    for k, cyc in enumerate(cycles):
        nodes_str = [
            f"({states[v][0]},{states[v][1]},{states[v][2]:+d})"
            for v in cyc
        ]
        cycle_edge_types = []
        for k2 in range(len(cyc)):
            u, v = cyc[k2], cyc[(k2 + 1) % len(cyc)]
            t = edge_type(states, u, v) if (u, v) in [(i, j) for (i, j) in edges] or (v, u) in [(i, j) for (i, j) in edges] else "?"
            # Robust: check sorted edge
            a, b = (u, v) if u < v else (v, u)
            if (a, b) in [(i, j) for (i, j) in edges]:
                cycle_edge_types.append(edge_type(states, a, b))
            else:
                cycle_edge_types.append("?")
        plaquette_info.append({
            "cycle_index": k,
            "length": len(cyc),
            "nodes": nodes_str,
            "edge_types_around_cycle": cycle_edge_types,
            "n_L_edges": cycle_edge_types.count("L"),
            "n_T_edges": cycle_edge_types.count("T"),
            "F_condon_shortley_units_2pi": str(F_A_pos[k]),
            "F_alternating_units_2pi": str(F_B_pos[k]),
            "F_unit_flux_units_2pi": str(F_C_pos[k]),
        })

    # ----- Verdict for this n_max -----
    def is_integer(f: Fraction) -> bool:
        return f.denominator == 1

    integer_check_A = is_integer(c1_A_pos) and is_integer(c1_A_neg)
    integer_check_B = is_integer(c1_B_pos) and is_integer(c1_B_neg)
    integer_check_C = is_integer(c1_C_pos)
    # Convention C: sanity check accepts either sign (cycle orientation
    # is a convention; the magnitude must match prescribed flux).
    sanity_C_check = abs(c1_C_pos) == Fraction(len(cycles), 4)

    return {
        "n_max": n_max,
        "graph": {
            "V": V,
            "E": E,
            "beta_0": beta_0,
            "beta_1": beta_1,
        },
        "cycles": {
            "n_fundamental_cycles": len(cycles),
            "cycle_lengths": [len(c) for c in cycles],
        },
        "convention_A_condon_shortley": {
            "description": "All edge phases = 0 (real-valued Condon-Shortley ladder ops)",
            "expected_c1": "0",
            "c1_default_orientation": str(c1_A_pos),
            "c1_flipped_orientation": str(c1_A_neg),
            "c1_is_integer": integer_check_A,
            "F_per_plaquette_units_2pi": [str(f) for f in F_A_pos],
        },
        "convention_B_alternating": {
            "description": "Z_2 twist: phase=pi on L-edges, 0 on T-edges (units of 2pi: 1/2 vs 0)",
            "expected_c1": "integer (each plaquette has even or odd L-count)",
            "c1_default_orientation": str(c1_B_pos),
            "c1_flipped_orientation": str(c1_B_neg),
            "c1_is_integer": integer_check_B,
            "F_per_plaquette_units_2pi": [str(f) for f in F_B_pos],
        },
        "convention_C_unit_flux": {
            "description": "Unit flux per plaquette = 1/4 (sanity check; assigns F=1/4 to each non-tree edge)",
            "expected_c1": (
                f"+/- {Fraction(len(cycles), 4)} "
                f"(sign depends on cycle orientation convention)"
            ),
            "c1_default_orientation": str(c1_C_pos),
            "c1_is_integer": integer_check_C,
            "c1_magnitude_matches_expected": sanity_C_check,
            "F_per_plaquette_units_2pi": [str(f) for f in F_C_pos],
        },
        "convention_D_fractional_flux_stress_test": {
            "description": (
                "Fractional flux 1/3 per plaquette (NON-physical stress test: "
                "verifies the algorithm faithfully tracks rational non-integer "
                "fluxes, so Convention A and B's integer result is meaningful)"
            ),
            "expected_c1_magnitude": str(expected_D_magnitude),
            "c1_default_orientation": str(c1_D_pos),
            "c1_magnitude_matches_expected": sanity_D_check,
            "c1_is_integer": is_integer(c1_D_pos),
            "F_per_plaquette_units_2pi": [str(f) for f in F_D_pos],
        },
        "plaquette_details": plaquette_info,
        "convention_sensitivity_notes": (
            "Convention A (all-zero, Condon-Shortley): c_1 = 0 trivially. "
            "Convention B (Z_2): integrality depends on parity of L-edges per cycle. "
            "Convention C (prescribed flux): verifies algorithm reproduces input flux."
        ),
    }


def main() -> None:
    out_path = Path("debug/data/discrete_c1_hopf_s3.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Discrete first Chern class on the Fock-S^3 Hopf graph")
    print("Sprint TS-E3 falsification target (2026-05-04)")
    print("=" * 70)

    result_2 = analyze_n_max(2)
    result_3 = analyze_n_max(3)

    # --- Top-level verdict aggregation ---
    integer_at_n2 = (
        result_2["convention_A_condon_shortley"]["c1_is_integer"]
        and result_2["convention_B_alternating"]["c1_is_integer"]
    )
    integer_at_n3 = (
        result_3["convention_A_condon_shortley"]["c1_is_integer"]
        and result_3["convention_B_alternating"]["c1_is_integer"]
    )
    sanity_C_n2 = result_2["convention_C_unit_flux"]["c1_magnitude_matches_expected"]
    sanity_C_n3 = result_3["convention_C_unit_flux"]["c1_magnitude_matches_expected"]

    if integer_at_n2 and integer_at_n3 and sanity_C_n2 and sanity_C_n3:
        verdict = "REINFORCED"
        verdict_note = (
            "Discrete c_1 is integer-valued at n_max=2 and n_max=3 under both "
            "Condon-Shortley (Convention A) and Z_2 alternating (Convention B) "
            "phase conventions. The sanity-check convention C reproduces the "
            "prescribed flux exactly. The Sprint TS-E3 master Mellin engine "
            "reading is REINFORCED: discrete topological invariants on finite "
            "GeoVac graphs are pi-free, with pi entering only through the "
            "M1 mechanism in the continuum limit."
        )
    else:
        verdict = "FALSIFIED"
        verdict_note = (
            "Discrete c_1 is not integer at one or more tested n_max under "
            "the natural conventions, OR the unit-flux sanity check failed. "
            "A fourth pi-source mechanism (M4) may need to be named."
        )

    output = {
        "sprint": "TS-E3 falsification target -- discrete c_1 on Hopf-S^3",
        "date": "2026-05-04",
        "verdict": verdict,
        "verdict_note": verdict_note,
        "n_max=2": result_2,
        "n_max=3": result_3,
    }

    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)

    # --- Console report ---
    for n_max, res in [(2, result_2), (3, result_3)]:
        print(f"\n--- n_max = {n_max} ---")
        print(f"V = {res['graph']['V']}, E = {res['graph']['E']}, "
              f"beta_0 = {res['graph']['beta_0']}, beta_1 = {res['graph']['beta_1']}")
        print(f"Number of fundamental cycles: {res['cycles']['n_fundamental_cycles']}")
        print(f"Cycle lengths: {res['cycles']['cycle_lengths']}")
        print(f"\nConvention A (Condon-Shortley, all phases = 0):")
        print(f"  c_1 = {res['convention_A_condon_shortley']['c1_default_orientation']}  "
              f"(integer: {res['convention_A_condon_shortley']['c1_is_integer']})")
        print(f"\nConvention B (Z_2 alternating, pi on L-edges):")
        print(f"  c_1 = {res['convention_B_alternating']['c1_default_orientation']}  "
              f"(integer: {res['convention_B_alternating']['c1_is_integer']})")
        print(f"\nConvention C (unit flux 1/4 per plaquette, sanity):")
        print(f"  expected c_1 = {res['convention_C_unit_flux']['expected_c1']}")
        print(f"  computed c_1 = {res['convention_C_unit_flux']['c1_default_orientation']}  "
              f"(magnitude matches: {res['convention_C_unit_flux']['c1_magnitude_matches_expected']})")
        d = res["convention_D_fractional_flux_stress_test"]
        print(f"\nConvention D (1/3 flux stress test, non-integer expected):")
        print(f"  expected |c_1| = {d['expected_c1_magnitude']}")
        print(f"  computed c_1 = {d['c1_default_orientation']}  "
              f"(integer: {d['c1_is_integer']}, magnitude matches: {d['c1_magnitude_matches_expected']})")

    print("\n" + "=" * 70)
    print(f"VERDICT: {verdict}")
    print(verdict_note)
    print("=" * 70)
    print(f"\nWritten: {out_path}")


if __name__ == "__main__":
    main()
