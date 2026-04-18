"""
Track RH-Q analysis: SU(2) Wilson lattice gauge theory on the S^3 Coulomb
Hopf graph.

Produces numerical data for the associated memo: plaquette counts at
max_n = 2, 3, 4; character-expansion partition function Z(beta);
Wilson loop expectation; U(1) reduction demonstration; Monte Carlo
verification.

Outputs:
  debug/data/su2_wilson_gauge_results.json
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List

import numpy as np

from geovac.lattice import GeometricLattice
from geovac.ihara_zeta import _count_closed_nonbacktracking_walks, _mobius_to_primitive
from geovac.su2_wilson_gauge import (
    OrientedEdge,
    diagonal_su2_from_phase,
    enumerate_oriented_edges,
    enumerate_plaquettes,
    expectation_wilson_loop,
    gauge_transform,
    monte_carlo_wilson_expectation,
    partition_function_character_expansion,
    plaquette_holonomy,
    su2_character,
    su2_character_coefficient,
    su2_random,
    u1_action_from_su2,
    wilson_action,
)


def main() -> None:
    out_path = Path(__file__).parent / "data" / "su2_wilson_gauge_results.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    results: Dict = {}

    # -----------------------------------------------------------------------
    # Section 1: Graph topology and plaquette enumeration
    # -----------------------------------------------------------------------
    print("=== Section 1: S^3 Coulomb graph topology and plaquettes ===\n")

    graph_stats = {}
    for nmax in [2, 3, 4]:
        lat = GeometricLattice(max_n=nmax)
        A = lat.adjacency.toarray()
        V = lat.num_states
        E = lat.num_edges

        from scipy.sparse.csgraph import connected_components
        n_comp, _ = connected_components(lat.adjacency, directed=False)
        beta_1 = E - V + n_comp

        # Ihara primitive walks (Möbius)
        tr_counts = _count_closed_nonbacktracking_walks(A, 10)
        prim = _mobius_to_primitive(tr_counts)

        # Direct plaquette enumeration
        plaqs_enum = {}
        for max_L in [4, 6, 8, 10]:
            plaqs = enumerate_plaquettes(A, max_length=max_L, both_orientations=False)
            counts_by_L: Dict[int, int] = {}
            for p in plaqs:
                counts_by_L[len(p)] = counts_by_L.get(len(p), 0) + 1
            plaqs_enum[max_L] = counts_by_L

        graph_stats[nmax] = {
            "V": V,
            "E": E,
            "c": int(n_comp),
            "beta_1": beta_1,
            "ihara_primitive_mobius": prim,
            "enumerated_plaquette_counts": plaqs_enum,
        }
        print(
            f"n_max={nmax}: V={V}, E={E}, c={n_comp}, beta_1={beta_1}"
        )
        print(f"  Ihara Mobius primitive walks: {prim}")
        print(f"  Enumerated (unoriented cycles, max_L=8): {plaqs_enum[8]}\n")

    results["graph_stats"] = graph_stats

    # -----------------------------------------------------------------------
    # Section 2: Character expansion Z(beta) at n_max = 3
    # -----------------------------------------------------------------------
    print("=== Section 2: Partition function Z(beta) at n_max = 3 ===\n")

    lat3 = GeometricLattice(max_n=3)
    A3 = lat3.adjacency.toarray()
    # Use only length-4 plaquettes for the character expansion test
    plaqs_L4 = enumerate_plaquettes(A3, max_length=4, both_orientations=False)
    print(f"Using {len(plaqs_L4)} primitive 4-cycles as plaquettes.\n")

    z_results = {}
    for beta in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
        Z = partition_function_character_expansion(plaqs_L4, beta, R_max=3)
        # Character coefficients for context
        c0 = su2_character_coefficient(0, beta)
        c1 = su2_character_coefficient(1, beta)
        c2 = su2_character_coefficient(2, beta)
        c3 = su2_character_coefficient(3, beta)
        z_results[f"beta_{beta}"] = {
            "beta": beta,
            "Z": Z,
            "c_0": c0,
            "c_1": c1,
            "c_2": c2,
            "c_3": c3,
            "c1_over_c0": c1 / c0 if c0 != 0 else None,
        }
        print(
            f"beta={beta:6.2f}: Z={Z:.6e}, c_0={c0:.4f}, c_1={c1:.4f}, "
            f"c_1/c_0={c1/c0:.4f}"
        )
    results["Z_beta"] = z_results

    # -----------------------------------------------------------------------
    # Section 3: Wilson loop expectation <W_C>
    # -----------------------------------------------------------------------
    print("\n=== Section 3: Wilson loop expectation ===\n")

    # Character-expansion estimate of <W_C> for a plaquette-sized loop
    sample_loop = plaqs_L4[0]
    w_results = {}
    print("Character-expansion <W_plaq>(beta):")
    for beta in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]:
        wl = expectation_wilson_loop(sample_loop, plaqs_L4, beta)
        w_results[f"char_beta_{beta}"] = {"beta": beta, "W": wl}
        print(f"  beta={beta:6.2f}: <W_C>={wl:.6f}")

    # -----------------------------------------------------------------------
    # Section 4: Monte Carlo comparison
    # -----------------------------------------------------------------------
    print("\n=== Section 4: Monte Carlo verification at n_max = 3 ===\n")

    oriented, _ = enumerate_oriented_edges(A3)
    # Use one forward link per undirected edge
    forward_edges = [e for e in oriented if e.source < e.target]

    mc_results = {}
    for beta in [0.5, 1.0, 2.0, 5.0]:
        mc_mean, mc_err = monte_carlo_wilson_expectation(
            sample_loop,
            plaqs_L4,
            forward_edges,
            beta=beta,
            n_samples=2000,
            n_thermalize=500,
            seed=12,
        )
        # For comparison, character-expansion leading order
        char_val = expectation_wilson_loop(sample_loop, plaqs_L4, beta)
        mc_results[f"beta_{beta}"] = {
            "beta": beta,
            "mc_mean": mc_mean,
            "mc_stderr": mc_err,
            "char_expansion_leading": char_val,
        }
        print(
            f"beta={beta:5.2f}: MC <W>={mc_mean:+.4f} ± {mc_err:.4f}, "
            f"char leading={char_val:+.4f}"
        )
    results["monte_carlo"] = mc_results

    # -----------------------------------------------------------------------
    # Section 5: U(1) reduction
    # -----------------------------------------------------------------------
    print("\n=== Section 5: U(1) reduction demonstration ===\n")

    rng = np.random.default_rng(42)
    u1_test = {}
    # Set a random U(1) phase on every forward oriented edge
    phases: Dict = {}
    for e in forward_edges:
        phases[(e.source, e.target)] = rng.uniform(-np.pi, np.pi)
    # Build SU(2) links as diagonal SU(2)
    diag_links = {k: diagonal_su2_from_phase(v) for k, v in phases.items()}

    beta = 2.0
    S_su2 = wilson_action(plaqs_L4, diag_links, beta)
    S_u1 = u1_action_from_su2(plaqs_L4, phases, beta)
    print(f"At beta={beta}, with random diagonal links:")
    print(f"  S_W(SU(2), diagonal) = {S_su2:.10f}")
    print(f"  S_W(U(1))            = {S_u1:.10f}")
    print(f"  Diff = {abs(S_su2 - S_u1):.2e}")
    u1_test["S_su2"] = S_su2
    u1_test["S_u1"] = S_u1
    u1_test["diff"] = abs(S_su2 - S_u1)

    # Also compare Haar-random SU(2) (should give LARGER action than diagonal only)
    rng2 = np.random.default_rng(7)
    full_links = {(e.source, e.target): su2_random(rng2) for e in forward_edges}
    S_full = wilson_action(plaqs_L4, full_links, beta)
    print(f"  S_W(SU(2), Haar-random) = {S_full:.6f}  (for reference)")
    u1_test["S_haar_random"] = S_full
    results["u1_reduction"] = u1_test

    # -----------------------------------------------------------------------
    # Section 6: Gauge invariance check
    # -----------------------------------------------------------------------
    print("\n=== Section 6: Gauge invariance check ===\n")

    rng3 = np.random.default_rng(9)
    gauge = {i: su2_random(rng3) for i in range(lat3.num_states)}
    full_links_gauged = gauge_transform(full_links, gauge)
    S_gauged = wilson_action(plaqs_L4, full_links_gauged, beta)
    print(f"Before gauge transformation: S = {S_full:.10f}")
    print(f"After gauge transformation:  S = {S_gauged:.10f}")
    print(f"Difference: {abs(S_full - S_gauged):.2e}  (should be < 1e-8)")

    results["gauge_invariance"] = {
        "S_before": S_full,
        "S_after": S_gauged,
        "diff": abs(S_full - S_gauged),
    }

    # -----------------------------------------------------------------------
    # Section 7: Confinement/area-law diagnostic
    # -----------------------------------------------------------------------
    print("\n=== Section 7: Area law vs perimeter law ===\n")

    # On our graph there are only short cycles; the 4-cycle is the atomic
    # plaquette. Compute <W_C> for the 4-cycle across a range of beta.
    print("Wilson loop <W_plaq>(beta) vs beta:")
    print(f"{'beta':>6} {'<W_C>':>10} {'-log<W_C>':>12}")
    area_data = []
    for beta in [0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0]:
        wl = expectation_wilson_loop(sample_loop, plaqs_L4, beta)
        nlog = -np.log(max(wl, 1e-30))
        area_data.append({"beta": beta, "W": wl, "-log_W": nlog})
        print(f"{beta:6.2f} {wl:10.6f} {nlog:12.4f}")

    results["area_law"] = area_data

    # -----------------------------------------------------------------------
    # Dump JSON
    # -----------------------------------------------------------------------
    with out_path.open("w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults written to: {out_path}")


if __name__ == "__main__":
    main()
