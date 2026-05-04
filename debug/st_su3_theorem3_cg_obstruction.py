"""
Sprint ST-SU3 Theorem 3 analysis: does the Wilson SU(3) construction
bypass Sprint 5's Clebsch-Gordan obstruction, or hit a different one?

Sprint 5 (debug/s5_gauge_structure_memo.md) showed that the
*adjacency-preserving* SU(3) action on the (N,0) representation tower
fails Wilson's "fixed group on every link" requirement: transitions
(N,0) -> (N+1,0) are CG intertwiners between distinct irreps, not group
elements.

This module evaluates whether the present Wilson SU(3) construction
bypasses or encounters that obstruction.

The verdict: BYPASSES, but with a structural caveat.

Why bypassed:
  - Wilson link variables U_e are external 3x3 SU(3) matrices, NOT
    representations of SU(3) acting on shell labels.
  - They live in the FUNDAMENTAL representation of SU(3), which has
    a fixed dimension 3 -- independent of the shell index N.
  - Gauge transformations g_v in SU(3) are also in the fundamental
    rep, so g_{src(e)} U_e g_{tgt(e)}^dagger is well-defined.

Why "structural caveat":
  - The fundamental-rep SU(3) link variables on the Bargmann graph
    have NO natural action on the (N, l, m_l) node labels. They
    do not couple to "matter" in the way that a fundamental rep on
    the (N, 0) shells would.
  - In Sprint 5, the question was about an SU(3) action that preserves
    the adjacency (i.e., respects the shell structure). The Wilson
    construction provides gauge degrees of freedom but does NOT
    provide matter coupling without an additional choice.
  - Specifically, to couple Wilson SU(3) to matter, one would need a
    map: node v -> SU(3)-rep on which g_v acts. The natural choice
    -- the (N,0) symmetric irrep -- has DIFFERENT dimensions on
    different shells (1, 3, 6, 10, 15, ...), so g_v acting on
    "matter at v" would require an irrep depending on the node, which
    is structurally how Sprint 5 hit the CG obstruction.

So: gauge sector works (bypass), matter coupling rediscovers the
obstruction (caveat). The Wilson SU(3) on Bargmann is a "pure gauge"
theory in the strict sense: gauge dynamics without matter coupling
to the (N,0) shells.

Compare to S^3 = SU(2) (Paper 30):
  - The S^3 Coulomb graph has nodes (n, l, m_l) carrying SU(2)
    spin-1/2 matter content (Paper 14 Tier 2 Dirac fermions).
  - SU(2) Wilson links act in the same representation as the matter
    (the spinor rep of SU(2)).
  - Matter coupling is natural: psi_v -> g_v psi_v with g_v in SU(2).

For SU(3) on S^5 Bargmann:
  - Nodes (N, l, m_l) live in (N, 0) symmetric SU(3) irreps of
    increasing dimension (N+1)(N+2)/2.
  - SU(3) Wilson links carry the fundamental 3-dimensional rep.
  - Matter coupling psi_v -> g_v psi_v requires g_v to act on the
    matter rep at v -- which is the (N, 0) rep, not the fundamental.
  - Resolving this requires either (a) projecting matter to a fixed
    representation, or (b) the CG-intertwiner structure of Sprint 5,
    which is NOT a Wilson lattice gauge theory.

This module documents the verdict computationally.

Output: debug/data/st_su3_theorem3_cg_obstruction.json
"""

import json
from pathlib import Path

import numpy as np

from geovac.nuclear.bargmann_graph import build_bargmann_graph
from geovac.su3_wilson_s5 import (
    bargmann_adjacency_dense,
    enumerate_oriented_edges,
    enumerate_plaquettes,
    su3_random,
    plaquette_holonomy,
    wilson_action,
    is_su3,
)


def construction_works_at_gauge_level(N_max: int):
    """
    Verify that the Wilson construction gives a well-defined gauge theory
    on the Bargmann graph at N_max:
      - Link variables in SU(3): yes (3x3 matrices, det=1, unitary).
      - Plaquette holonomies: stay in SU(3).
      - Action gauge-invariant under node-local SU(3) transformations.
      - Partition function exists (finite Haar integrals).

    All four checks pass by construction; we verify the first three
    numerically.
    """
    A = bargmann_adjacency_dense(N_max)
    oriented, _ = enumerate_oriented_edges(A)
    plaqs = enumerate_plaquettes(A, max_length=4, both_orientations=False)
    if not plaqs:
        return {"N_max": N_max, "skipped": "no plaquettes"}

    forward = [(e.source, e.target) for e in oriented if e.source < e.target]
    rng = np.random.default_rng(7)
    links = {k: su3_random(rng) for k in forward}

    # Test 1: link variables in SU(3)
    all_su3 = all(is_su3(U, atol=1e-10) for U in links.values())

    # Test 2: plaquette holonomy in SU(3)
    plaq_holos_su3 = []
    for P in plaqs[:50]:  # check first 50
        U_P = plaquette_holonomy(P, links)
        plaq_holos_su3.append(is_su3(U_P, atol=1e-8))
    all_plaq_su3 = all(plaq_holos_su3)

    # Test 3: gauge invariance
    from geovac.su3_wilson_s5 import gauge_transform
    V = A.shape[0]
    gauge = {i: su3_random(rng) for i in range(V)}
    S_before = wilson_action(plaqs, links, beta=2.0)
    links_g = gauge_transform(links, gauge)
    S_after = wilson_action(plaqs, links_g, beta=2.0)
    gauge_inv = abs(S_before - S_after) < 1e-8

    return {
        "N_max": N_max,
        "n_plaquettes": len(plaqs),
        "n_forward_edges": len(forward),
        "links_in_SU3": all_su3,
        "plaquette_holonomies_in_SU3": all_plaq_su3,
        "gauge_invariance": gauge_inv,
        "S_before_gauge": S_before,
        "S_after_gauge": S_after,
        "S_diff": abs(S_before - S_after),
        "verdict": (
            "Wilson construction WORKS at gauge level"
            if (all_su3 and all_plaq_su3 and gauge_inv)
            else "FAILS at gauge level"
        ),
    }


def matter_coupling_obstruction_analysis(N_max: int):
    """
    Show that Wilson SU(3) link variables (in fundamental rep) do NOT
    naturally couple to matter on the (N, 0) shells of the Bargmann graph,
    because the matter representations have different dimensions.
    """
    g = build_bargmann_graph(N_max)
    shell_dims = {}
    for (N, _, _) in g.nodes:
        if N not in shell_dims:
            shell_dims[N] = 0
        shell_dims[N] += 1

    # Each shell N hosts a (N, 0) SU(3) symmetric rep of dimension (N+1)(N+2)/2.
    expected_dim = {N: (N + 1) * (N + 2) // 2 for N in shell_dims}
    matches = {N: shell_dims[N] == expected_dim[N] for N in shell_dims}

    # Wilson SU(3) links: each carries fundamental 3-rep, dimension 3.
    fundamental_dim = 3

    # Matter coupling psi_v -> g_v psi_v requires g_v to act on the
    # matter rep at v. Natural choice: matter at v lives in the (N(v), 0)
    # rep where N(v) is the shell. But then g_v takes (N,0) -> (N,0), and
    # matter on different shells transforms in different reps.
    #
    # The fundamental link variables connect shells N and N+1. For psi_v
    # to gain an SU(3) gauge phase from a link e = (v, w) with v on shell
    # N and w on shell N+1, we need a rule for U_e * psi_v that takes a
    # (N, 0)-vector to a (N+1, 0)-vector. This is the CG INTERTWINER --
    # exactly Sprint 5's obstruction.

    return {
        "N_max": N_max,
        "shell_dimensions_actual": shell_dims,
        "shell_dimensions_expected_(N,0)_rep": expected_dim,
        "matches": matches,
        "wilson_link_rep_dimension": fundamental_dim,
        "matter_rep_dimensions": list(set(expected_dim.values())),
        "matter_coupling_obstruction": (
            "Matter (N,0) reps have different dimensions on different shells. "
            "Wilson links (3-dim fundamental) cannot directly couple matter "
            "across shells without invoking CG intertwiners (Sprint 5 obstruction)."
        ),
    }


def comparison_with_su2_on_s3():
    """
    Side-by-side comparison: SU(2) Wilson on S^3 (Paper 30) vs SU(3)
    Wilson on S^5 (this sprint).
    """
    return {
        "SU(2)_on_S3": {
            "manifold_equals_group": "S^3 = SU(2) (3-sphere is simply connected SU(2))",
            "matter_rep": "spin-1/2 fundamental of SU(2), dim=2 -- SAME as link rep",
            "natural_matter_coupling": "YES (psi_v -> g_v psi_v with same rep)",
            "matter_sector": "Paper 14 Tier 2 Dirac fermions in (kappa, m_j) basis",
            "rank": 1,
            "cartan_torus": "U(1) -- single phase",
        },
        "SU(3)_on_S5_Bargmann": {
            "manifold_equals_group": "S^5 != SU(3) (S^5 is the unit sphere in C^3, SU(3) is 8-dim)",
            "matter_rep": "(N,0) symmetric reps of SU(3), dim=(N+1)(N+2)/2 -- DIFFERENT per shell",
            "natural_matter_coupling": (
                "NO -- matter rep changes per shell; coupling requires CG intertwiners"
                " (Sprint 5 negative)"
            ),
            "matter_sector": "the (N,0) shells themselves; no natural Dirac matter",
            "rank": 2,
            "cartan_torus": "U(1) x U(1) -- two phases",
        },
        "structural_punchline": (
            "On S^3 the manifold IS the gauge group, so matter and gauge live in"
            " the same Hilbert space. On S^5 the gauge group is external and"
            " 8-dimensional while the manifold is 5-dimensional; matter lives in a"
            " representation tower whose dimension grows with shell index. Wilson"
            " gauge is well-defined as a pure-gauge theory on either; matter"
            " coupling is natural on S^3 but encounters Sprint 5's CG obstruction"
            " on S^5."
        ),
    }


if __name__ == "__main__":
    out = {
        "sprint": "ST-SU3",
        "theorem": "Theorem 3: CG obstruction status",
        "question": (
            "Does the Wilson SU(3) construction bypass Sprint 5 Track S5's "
            "(N,0)-tower CG obstruction?"
        ),
    }

    print("=" * 60)
    print("Theorem 3 — CG Obstruction Status")
    print("=" * 60)

    print("\n[Part A] Construction works at the gauge level?")
    out["gauge_level"] = {}
    for N_max in [2, 3]:
        r = construction_works_at_gauge_level(N_max)
        out["gauge_level"][f"N_max={N_max}"] = r
        if "skipped" in r:
            print(f"  N_max={N_max}: SKIPPED")
            continue
        print(f"  N_max={N_max}: {r['verdict']}")
        print(f"    links in SU(3): {r['links_in_SU3']}")
        print(f"    plaquette holonomies in SU(3): {r['plaquette_holonomies_in_SU3']}")
        print(f"    gauge invariance: {r['gauge_invariance']}, "
              f"|S_diff| = {r['S_diff']:.2e}")

    print("\n[Part B] Matter coupling obstruction analysis")
    out["matter_coupling"] = {}
    for N_max in [2, 3]:
        r = matter_coupling_obstruction_analysis(N_max)
        out["matter_coupling"][f"N_max={N_max}"] = r
        print(f"  N_max={N_max}: shell_dims = {r['shell_dimensions_actual']}")
        print(f"    matter rep dimensions = {r['matter_rep_dimensions']}")
        print(f"    fundamental link rep dim = {r['wilson_link_rep_dimension']}")

    print("\n[Part C] Comparison with SU(2) on S^3")
    out["comparison"] = comparison_with_su2_on_s3()
    print(out["comparison"]["structural_punchline"])

    out["overall_verdict"] = (
        "BYPASSES at the gauge level (Wilson construction works), "
        "rediscovers the CG obstruction at the matter-coupling level. "
        "Both verdicts coexist: SU(3) Wilson is a self-consistent gauge-only "
        "theory on the Bargmann S^5 graph; matter coupling to the (N,0) shells "
        "encounters the same CG-intertwiner structure that Sprint 5 identified."
    )

    print(f"\nOVERALL: {out['overall_verdict']}")

    out_path = Path("debug/data/st_su3_theorem3_cg_obstruction.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"Saved: {out_path}")
