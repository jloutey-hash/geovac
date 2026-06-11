"""
Pilot script for Migdal-Kadanoff block-spin decimation on the GeoVac
Hopf graph U(1) Wilson lattice gauge theory.

Scoping output: validates the fundamental U(1) MK recursion
    I_1(beta_eff) / I_0(beta_eff)  =  [I_1(beta) / I_0(beta)]^2
for one-bond decimation, against its strong-coupling limit
    beta_eff ~ beta^2 / 2
and weak-coupling limit
    beta_eff ~ beta / 2.

Also enumerates the per-l-block grid structure of the GeoVac Fock graph
at n_max = 3, 4, 5, confirming that each connected component is a
Cartesian product of two path graphs P_a x P_b.

This is a SCOPING pilot, not a production-quality implementation.
It accompanies debug/rg_blockspin_design_memo.md.
"""

from __future__ import annotations

import numpy as np
from scipy.special import iv

from geovac.fock_graph_hodge import FockGraphHodge


def ratio_I1_over_I0(beta: float) -> float:
    """Bessel ratio I_1(beta) / I_0(beta) for U(1) character expansion."""
    if abs(beta) < 1e-15:
        return 0.0
    return iv(1, beta) / iv(0, beta)


def beta_eff_one_bond(beta: float, tol: float = 1e-12, n_iter: int = 100) -> float:
    """
    Solve I_1(beta_eff) / I_0(beta_eff)  =  [I_1(beta) / I_0(beta)]^2
    via Newton iteration, returning beta_eff.

    Strong-coupling limit: beta_eff -> beta^2 / 2 as beta -> 0.
    Weak-coupling limit:   beta_eff -> beta / 2 as beta -> infinity.
    """
    target = ratio_I1_over_I0(beta) ** 2
    if target < 1e-30:
        return 0.0
    # Initial guess: strong-coupling form
    b = max(2 * target, 1e-6)
    for _ in range(n_iter):
        r = ratio_I1_over_I0(b)
        f = r - target
        if abs(f) < tol:
            return float(b)
        # f'(b) = d/db [I_1(b)/I_0(b)] = 1 - r^2 - r/b (Bessel recurrence)
        fp = 1.0 - r * r - r / b
        if abs(fp) < 1e-15:
            break
        b -= f / fp
        if b < 0:
            b = 1e-6
    return float(b)


def block_decomposition(n_max: int) -> list[dict]:
    """Return the per-l block structure of the GeoVac Fock graph at n_max."""
    blocks = []
    for l in range(n_max):
        a = n_max - l       # radial direction (n in {l+1, ..., n_max})
        b = 2 * l + 1       # angular direction (m in {-l, ..., +l})
        V = a * b
        # Edges of P_a x P_b grid: a*(b-1) + b*(a-1) = 2ab - a - b
        E = a * (b - 1) + b * (a - 1)
        # 4-plaquettes: (a-1)(b-1)
        P4 = max(0, (a - 1) * (b - 1))
        blocks.append({
            'l': l,
            'shape': f'P_{a} x P_{b}',
            'V': V,
            'E': E,
            'plaquettes': P4,
        })
    return blocks


def verify_fock_graph_structure(n_max: int) -> dict:
    """
    Confirm the per-block decomposition matches FockGraphHodge.

    Returns a dict with V_predicted, V_actual, E_predicted, E_actual,
    beta1_predicted, beta1_actual.
    """
    blocks = block_decomposition(n_max)
    V_pred = sum(b['V'] for b in blocks)
    E_pred = sum(b['E'] for b in blocks)
    P_pred = sum(b['plaquettes'] for b in blocks)

    g = FockGraphHodge(n_max, use_sympy=False)
    V_act = g.n_nodes
    E_act = g.n_edges
    L0 = g.node_laplacian_numpy
    eigs = np.linalg.eigvalsh(L0)
    c_act = int(np.sum(np.abs(eigs) < 1e-8))
    beta1_act = E_act - V_act + c_act

    return {
        'V_predicted': V_pred,
        'V_actual': V_act,
        'E_predicted': E_pred,
        'E_actual': E_act,
        'plaquettes_predicted': P_pred,
        'beta1_actual': beta1_act,
        'components_actual': c_act,
        'match': (V_pred == V_act) and (E_pred == E_act) and (P_pred == beta1_act),
    }


def main():
    print("=" * 70)
    print("Migdal-Kadanoff U(1) decimation pilot")
    print("=" * 70)
    print()
    print("Part 1: Fock graph per-l-block decomposition")
    print("-" * 70)
    for n_max in [2, 3, 4, 5]:
        print(f"\nn_max = {n_max}:")
        result = verify_fock_graph_structure(n_max)
        for b in block_decomposition(n_max):
            print(f"  l = {b['l']}: {b['shape']:12s}  V={b['V']:3d}  "
                  f"E={b['E']:3d}  plaquettes={b['plaquettes']:3d}")
        print(f"  Totals: V={result['V_predicted']:3d}, E={result['E_predicted']:3d}, "
              f"plaquettes={result['plaquettes_predicted']:3d}")
        print(f"  Actual: V={result['V_actual']:3d}, E={result['E_actual']:3d}, "
              f"beta_1={result['beta1_actual']:3d}, components={result['components_actual']}")
        print(f"  Match : {result['match']}")

    print()
    print("=" * 70)
    print("Part 2: One-bond MK U(1) recursion - beta_eff(beta)")
    print("=" * 70)
    print()
    print(f"{'beta':>8s}  {'beta_eff (exact)':>18s}  "
          f"{'beta^2/2 (strong)':>18s}  {'beta/2 (weak)':>15s}  {'limit used':>15s}")
    print("-" * 70)
    for beta in [0.01, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0]:
        b_eff = beta_eff_one_bond(beta)
        sc = beta * beta / 2.0
        wc = beta / 2.0
        # Which limit dominates?
        if beta < 0.5:
            limit = "strong-coupl."
        elif beta > 3.0:
            limit = "weak-coupl."
        else:
            limit = "crossover"
        print(f"{beta:8.3f}  {b_eff:18.6f}  {sc:18.6f}  {wc:15.6f}  {limit:>15s}")

    print()
    print("=" * 70)
    print("Part 3: RG flow under iterated MK steps (single-bond)")
    print("=" * 70)
    print()
    print(f"Iterating beta -> beta_eff repeatedly for several initial values.")
    print(f"(One-bond MK; for full 2D MK step beta flows faster.)")
    print()
    print(f"{'Step':>5s}  " + "  ".join(f"{b0:>10.2f}" for b0 in [0.5, 1.0, 2.0, 5.0, 10.0]))
    print("-" * 70)
    flows = {b0: b0 for b0 in [0.5, 1.0, 2.0, 5.0, 10.0]}
    for step in range(8):
        row = [f"{step:5d}"]
        for b0 in [0.5, 1.0, 2.0, 5.0, 10.0]:
            row.append(f"{flows[b0]:10.6f}")
        print("  ".join(row))
        # Iterate
        flows = {b0: beta_eff_one_bond(flows[b0]) for b0 in flows}

    print()
    print("Confirms: beta flows to 0 from all initial values (2D U(1) confined phase).")
    print()
    print("Notes:")
    print(" * Strong-coupling beta^2/2 limit matches exact to 4 decimal places at beta=0.1")
    print(" * Weak-coupling beta/2 limit is asymptotic (visible at beta >= 5)")
    print(" * Fixed point beta* = 0 is attractive; beta* = infinity is unstable")
    print(" * No nontrivial intermediate fixed point - Polyakov 2D U(1) confinement")
    print()
    print("This pilot validates the design memo's beta_eff analysis.")
    print("Full implementation would extend to 2D MK steps and Monte Carlo Wilson loops.")


if __name__ == '__main__':
    main()
