"""
XCWG Migdal-Kadanoff block-spin RG on the Dirac-Rule-B graph
=============================================================

Sprint XCWG (May 2026). Companion to:
  - debug/rg_blockspin_design_memo.md  (scalar Fock baseline)
  - debug/xcwg_rule_b_infrastructure_memo.md  (Rule B substrate)
  - debug/xcwg_u1_wilson_rule_b_design_memo.md  (Rule B Wilson construction)

This script implements Migdal-Kadanoff (MK) bond-moving block-spin
decimation for U(1) Wilson lattice gauge theory on the Dirac-Rule-B graph
of Paper 29 §RH-C. The recursion is verified against textbook
strong/weak-coupling limits and is compared to the 2D scalar baseline of
debug/rg_blockspin_pilot.py.

Key structural fact distinguishing Rule B from the scalar Fock graph:
  Rule B has 100% cross-l edges and 100% cross-shell plaquettes
  (xcwg_rule_b_infrastructure_memo.md §3-4). The bond-moving step on
  Rule B therefore mixes l-shells inseparably, in sharp contrast to
  the scalar Fock graph's disjoint 2D-grid decomposition.

The MK recursion adapted for a generic graph with plaquette-edge
incidence q (the average number of plaquettes meeting at each edge) is

      I_1(beta_eff) / I_0(beta_eff)  =  [I_1(beta) / I_0(beta)]^{k_eff}

where k_eff = q - 1 (each integrated link is shared by q plaquettes;
the bond-moving step bundles q-1 *parallel* plaquette-borders onto the
surviving link, then the link is integrated out, returning the q-th
power minus the bookkeeping).

For 2D square lattice (q=2): k_eff = 1, but the standard bond-moving
+ decimation step in fact gives k_eff = 2 because the bond-move *first*
doubles the bond, then decimation across the doubled bond integrates
the link with two adjacent plaquettes -- the effective exponent is the
*product* (b^{d-1} for block size b, d dimensions). At b=2, d=2: k=2.

For 3D cubic lattice (b=2, d=3): k = b^{d-1} = 4 (each integrated bond
sees 4 surrounding plaquettes via the bond-move).

The Rule B value of "q" (plaquettes per edge) is 9.4 at n_max=3 and
2.2 at n_max=2 (from xcwg_rule_b_infrastructure_memo.md §5). The
*effective* k for the MK recursion on Rule B depends on which bond is
being decimated and how its surrounding plaquettes are distributed.

Output:
  debug/data/xcwg_mk_blockspin_rule_b.json
"""

from __future__ import annotations

import json
import os
import sys
from collections import Counter, defaultdict
from typing import Dict, List, Tuple

import numpy as np
from scipy.special import iv

# Make geovac importable.
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from geovac.ihara_zeta_dirac import build_dirac_s3_graph  # noqa: E402


# =============================================================================
# §1 Plaquette enumeration (re-used from infrastructure script)
# =============================================================================

def primitive_4_plaquettes(A: np.ndarray) -> List[Tuple[int, int, int, int]]:
    """Enumerate primitive non-backtracking closed 4-walks (4-cycles)."""
    n = A.shape[0]
    seen = set()
    cycles = []
    for v0 in range(n):
        neigh0 = [int(j) for j in np.nonzero(A[v0])[0]]
        for v1 in neigh0:
            if v1 == v0:
                continue
            neigh1 = [int(j) for j in np.nonzero(A[v1])[0] if int(j) != v0]
            for v2 in neigh1:
                if v2 == v0 or v2 == v1:
                    continue
                neigh2 = [int(j) for j in np.nonzero(A[v2])[0]
                          if int(j) != v1 and int(j) != v0]
                for v3 in neigh2:
                    if v3 == v0 or v3 == v1 or v3 == v2:
                        continue
                    if A[v3, v0] == 0:
                        continue
                    cyc = (v0, v1, v2, v3)
                    rotations = [cyc, cyc[1:] + cyc[:1],
                                 cyc[2:] + cyc[:2], cyc[3:] + cyc[:3]]
                    reflected = [tuple(reversed(r)) for r in rotations]
                    canon = min(rotations + reflected)
                    if canon not in seen:
                        seen.add(canon)
                        cycles.append(canon)
    return cycles


# =============================================================================
# §2 Plaquette-edge incidence statistics
# =============================================================================

def plaquette_edge_incidence(
    cycles: List[Tuple[int, int, int, int]],
    n_vertices: int,
) -> Tuple[Dict[Tuple[int, int], int], Dict[str, float]]:
    """For each undirected edge {u,v}, count how many plaquettes contain it.

    Returns
    -------
    edge_plaq_count : dict { (u, v) -> int }   (u < v)
        Number of plaquettes (in the input list) traversing edge {u, v}.
    stats : dict
        Aggregate statistics: total_edges, total_plaquettes, mean, min, max,
        median, std of plaquette-per-edge incidence, and the effective k_eff
        used downstream.
    """
    edge_count: Counter = Counter()
    for cyc in cycles:
        for i in range(4):
            u, v = cyc[i], cyc[(i + 1) % 4]
            key = (min(u, v), max(u, v))
            edge_count[key] += 1
    if not edge_count:
        return {}, {
            "total_edges_with_plaquette": 0,
            "total_plaquettes": len(cycles),
            "mean_plaq_per_edge": 0.0,
            "min_plaq_per_edge": 0,
            "max_plaq_per_edge": 0,
            "median_plaq_per_edge": 0.0,
            "std_plaq_per_edge": 0.0,
            "k_eff_mean": 0.0,
            "k_eff_geometric_mean": 0.0,
        }
    counts = np.array(list(edge_count.values()))
    # k_eff = q - 1 (q = plaquettes per edge in the link-integration step)
    # See module docstring: an integrated link is shared by q plaquettes;
    # the recursion couples the surviving links via the q-th-power character
    # identity, leaving exponent (q-1) acting on the new effective bond.
    # For 2D square (q=2): k_eff = 1, but the bond-move *also* doubles the
    # bond beforehand, giving effective exponent 2 for the standard 2D MK
    # recursion. We track both q and k_eff.
    mean_q = float(counts.mean())
    geom_q = float(np.exp(np.log(counts).mean()))  # log-mean ≈ "typical" q
    stats = {
        "total_edges_with_plaquette": int(len(counts)),
        "total_plaquettes": int(len(cycles)),
        "mean_plaq_per_edge": mean_q,
        "min_plaq_per_edge": int(counts.min()),
        "max_plaq_per_edge": int(counts.max()),
        "median_plaq_per_edge": float(np.median(counts)),
        "std_plaq_per_edge": float(counts.std()),
        "geometric_mean_plaq_per_edge": geom_q,
        # k_eff variants for the MK recursion ratio^k. Two natural choices:
        # (a) mean q-1 ("average bond-share minus self"); see §3 derivation.
        # (b) full q ("each integrated bond brings q parallel plaquettes
        #     into the recursion product"). In 2D-square (q=2) the textbook
        #     value is k=2 (matches choice (b)); we use (b) as the
        #     canonical k_eff but report both.
        "k_eff_mean_minus_1": mean_q - 1.0,
        "k_eff_mean": mean_q,
        "k_eff_geometric": geom_q,
        "k_eff_max": int(counts.max()),
    }
    return dict(edge_count), stats


# =============================================================================
# §3 The MK recursion (closed form via Bessel-Fourier identity)
# =============================================================================

def ratio_I1_over_I0(beta: float) -> float:
    """U(1) character expansion ratio I_1(beta) / I_0(beta).

    Strong-coupling limit (beta -> 0):  ratio ~ beta/2 - beta^3/16 + ...
    Weak-coupling limit (beta -> inf):  ratio ~ 1 - 1/(2 beta) + ...
    """
    if abs(beta) < 1e-15:
        return 0.0
    return iv(1, beta) / iv(0, beta)


def invert_bessel_ratio(target_ratio: float, tol: float = 1e-12,
                        n_iter: int = 100) -> float:
    """Solve I_1(beta_eff) / I_0(beta_eff) = target_ratio for beta_eff.

    Uses Newton iteration. The derivative is
      d/d(beta) [I_1/I_0] = 1 - (I_1/I_0)^2 - (I_1/I_0)/beta
    from the Bessel recurrence.
    """
    if target_ratio < 1e-30:
        return 0.0
    if target_ratio >= 1.0:
        return float("inf")
    # Initial guess: invert strong-coupling form target = b/2
    b = max(2.0 * target_ratio, 1e-6)
    for _ in range(n_iter):
        r = ratio_I1_over_I0(b)
        f = r - target_ratio
        if abs(f) < tol:
            return float(b)
        fp = 1.0 - r * r - r / b
        if abs(fp) < 1e-15:
            break
        delta = f / fp
        # Damped Newton step to keep b > 0
        b_new = b - delta
        if b_new <= 0:
            b_new = b * 0.5
        b = b_new
    return float(b)


def beta_eff_mk(beta: float, k_eff: float) -> float:
    """One MK decimation step with effective exponent k_eff.

    Solves ratio(beta_eff) = ratio(beta)^k_eff.
    """
    r = ratio_I1_over_I0(beta)
    target = r ** k_eff
    return invert_bessel_ratio(target)


# =============================================================================
# §4 Verification against strong/weak-coupling limits
# =============================================================================

def strong_coupling_limit(beta: float, k: float) -> float:
    """beta_eff at strong coupling: ratio ~ beta/2, so beta_eff/2 = (beta/2)^k.

    Therefore beta_eff = 2 * (beta/2)^k = 2^(1-k) * beta^k.
    """
    return 2.0 ** (1.0 - k) * beta ** k


def weak_coupling_limit(beta: float, k: float) -> float:
    """beta_eff at weak coupling: ratio ~ 1 - 1/(2 beta).

    [ratio]^k ~ 1 - k/(2 beta), so 1 - 1/(2 beta_eff) = 1 - k/(2 beta),
    giving beta_eff = beta / k.
    """
    return beta / k


def verify_limits(k: float, beta_range: List[float]) -> Dict:
    """Tabulate exact MK beta_eff vs strong-coupling and weak-coupling limits.

    For each beta in beta_range, computes:
      - exact beta_eff
      - strong-coupling prediction 2^(1-k) * beta^k
      - weak-coupling prediction beta / k
      - rel error to each
    """
    rows = []
    for b in beta_range:
        b_eff = beta_eff_mk(b, k)
        b_sc = strong_coupling_limit(b, k)
        b_wc = weak_coupling_limit(b, k)
        # rel err vs each limit
        if b_sc > 0:
            err_sc = (b_eff - b_sc) / b_sc
        else:
            err_sc = float("nan")
        if b_wc > 0:
            err_wc = (b_eff - b_wc) / b_wc
        else:
            err_wc = float("nan")
        rows.append({
            "beta": b,
            "beta_eff": b_eff,
            "strong_coupling": b_sc,
            "weak_coupling": b_wc,
            "rel_err_strong": err_sc,
            "rel_err_weak": err_wc,
        })
    return rows


# =============================================================================
# §5 Fixed-point search
# =============================================================================

def find_fixed_points(k: float, beta_lo: float = 1e-3,
                      beta_hi: float = 100.0,
                      n_grid: int = 2001,
                      tol: float = 1e-10) -> Dict:
    """Search beta_eff - beta = 0 on [beta_lo, beta_hi].

    For an MK U(1) recursion with k > 1, the recursion always drives
    beta downward (Polyakov-style 2D confinement), so there is no
    non-trivial fixed point. For k < 1 (which would be unphysical from
    bond-moving but is allowed mathematically), the recursion could
    have a non-trivial fixed point.

    For k = 1 the recursion is the identity at every beta (every beta
    is a "fixed point").

    For our purposes the key check is whether g(beta) := beta_eff(beta)
    - beta has a sign change in (beta_lo, beta_hi).
    """
    grid = np.geomspace(beta_lo, beta_hi, n_grid)
    g = np.array([beta_eff_mk(b, k) - b for b in grid])
    sign_changes = []
    for i in range(len(grid) - 1):
        if g[i] == 0:
            sign_changes.append((grid[i], grid[i]))
        elif g[i] * g[i + 1] < 0:
            sign_changes.append((float(grid[i]), float(grid[i + 1])))

    fixed_points: List[Dict] = []
    for a, b in sign_changes:
        # Bisection
        lo, hi = a, b
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            g_mid = beta_eff_mk(mid, k) - mid
            g_lo = beta_eff_mk(lo, k) - lo
            if abs(g_mid) < tol:
                break
            if g_mid * g_lo < 0:
                hi = mid
            else:
                lo = mid
        # Determine stability: g'(beta*) = d beta_eff/d beta - 1
        # Stable fixed point if d beta_eff/d beta < 1.
        # Numerical derivative.
        delta = 1e-5
        deriv = (beta_eff_mk(mid + delta, k) - beta_eff_mk(mid - delta, k)) / (
            2 * delta)
        stable = deriv < 1.0
        fixed_points.append({
            "beta_star": mid,
            "d_beta_eff_d_beta": deriv,
            "stable": stable,
        })

    # Trivial fixed points: beta=0 (always) and beta=inf (if reachable).
    # Just note that beta=0 is a fixed point of every U(1) MK recursion.
    return {
        "k_eff": k,
        "grid_beta_lo": float(beta_lo),
        "grid_beta_hi": float(beta_hi),
        "n_grid": int(n_grid),
        "n_sign_changes": len(sign_changes),
        "fixed_points": fixed_points,
        "trivial_fixed_point_beta_0": True,
        # The recursion drives beta down for k > 1; up for k < 1.
        "g_at_extremes": {
            "beta_lo": float(g[0]),
            "beta_hi": float(g[-1]),
        },
    }


# =============================================================================
# §6 Sanity check against scalar Fock 2D MK (k=2)
# =============================================================================

def sanity_check_scalar_2D() -> Dict:
    """Reproduce debug/rg_blockspin_pilot.py's 2D MK result at k=2.

    Strong-coupling: beta_eff = beta^2 / 2.
    Weak-coupling: beta_eff = beta / 2.
    """
    beta_grid = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 50.0]
    rows = verify_limits(2.0, beta_grid)
    return {
        "k": 2.0,
        "description": "2D-square MK recursion (Polyakov 1977 confinement)",
        "data": rows,
    }


# =============================================================================
# §7 RG flow under iterated MK steps
# =============================================================================

def iterate_mk(beta_0: float, k: float, n_steps: int = 12) -> List[float]:
    """Iterate beta_eff for n_steps from initial beta_0."""
    out = [beta_0]
    b = beta_0
    for _ in range(n_steps):
        b = beta_eff_mk(b, k)
        out.append(b)
        if b < 1e-15:
            # Stuck at 0
            out.extend([0.0] * (n_steps - len(out) + 1))
            break
    return out[:n_steps + 1]


# =============================================================================
# §8 Edge-by-edge MK on the actual Rule B graph at n_max = 2
# =============================================================================

def edgewise_mk_step(
    A: np.ndarray,
    cycles: List[Tuple[int, int, int, int]],
    beta: float,
) -> Dict:
    """Decimate ONE link of the Rule B graph and compute the local beta_eff
    on the surviving plaquettes.

    For each undirected edge e = {u,v}, find all plaquettes containing e.
    If q_e plaquettes share this edge, integrating out the U_e link couples
    the q_e plaquettes into a single composite Wilson loop. To leading
    order, the resulting effective Boltzmann weight on the *product loop*
    (which traverses all q_e plaquettes but skips e itself) satisfies

       ratio(beta_loc^{eff}) = ratio(beta)^{q_e - 1}    if q_e >= 1

    where the (q_e - 1) comes from the q_e - 1 plaquettes that remain
    coupled to the (q_e)-th plaquette via the integrated bond. For
    q_e = 2 (2D square): k = 1 in this rule, but combined with bond-moving
    that bundles two parallel bonds, the canonical 2D MK step gives k = 2.
    For q_e = 1: the integrated bond was on the boundary of a single
    plaquette; integrating it out destroys the plaquette and contributes
    nothing to the surviving plaquettes -- effective k = 0.

    This is the EDGEWISE accounting; it differs from the GLOBAL k_eff
    by being per-edge weighted. We report both.
    """
    edge_count, stats = plaquette_edge_incidence(cycles, A.shape[0])
    if not edge_count:
        return {"n_edges_with_plaquette": 0}

    # For each edge, compute the local beta_eff under MK with k = q_e - 1.
    per_edge: List[Dict] = []
    qs = []
    for (u, v), q in sorted(edge_count.items()):
        k_loc = max(0.0, q - 1.0)  # decimating this edge: surviving plaquettes
        if k_loc == 0:
            b_eff_loc = 0.0  # plaquette destroyed
        else:
            b_eff_loc = beta_eff_mk(beta, k_loc)
        per_edge.append({
            "edge": [int(u), int(v)],
            "q_plaq": int(q),
            "k_loc": k_loc,
            "beta_eff_loc": b_eff_loc,
        })
        qs.append(q)
    qs = np.array(qs)
    # The total post-decimation effective coupling is *not* the sum or
    # arithmetic mean of per-edge beta_eff_loc -- it depends on how the
    # composite loops overlap globally. The natural single-number summary
    # is the per-edge mean and median, plus the *global* k_eff (mean q).
    return {
        "n_edges_with_plaquette": int(len(per_edge)),
        "q_min": int(qs.min()),
        "q_max": int(qs.max()),
        "q_mean": float(qs.mean()),
        "q_median": float(np.median(qs)),
        # Mean beta_eff_loc across edges, excluding boundary q=1 edges
        "beta_eff_loc_mean_q_gt_1": float(np.mean([
            r["beta_eff_loc"] for r in per_edge if r["q_plaq"] > 1
        ])) if any(r["q_plaq"] > 1 for r in per_edge) else 0.0,
        "beta_eff_loc_max": float(np.max([r["beta_eff_loc"]
                                           for r in per_edge])),
        # Sample first 6 edges for transparency
        "sample_per_edge": per_edge[:6],
    }


# =============================================================================
# §9 Within-l restriction (sanity check: does Rule B's recursion reduce to
#    the scalar 2D-grid k=2 if we drop cross-l plaquettes?)
# =============================================================================

def within_l_restriction_sanity(
    A: np.ndarray,
    cycles: List[Tuple[int, int, int, int]],
    l_per_node: List[int],
) -> Dict:
    """Identify within-l plaquettes (all 4 vertices same l) and check
    whether their plaquette/edge incidence pattern matches the 2D scalar
    grid (which has q_e = 2 for interior bonds, q_e = 1 for boundary).

    Rule B has Δl = ±1 on every edge, so within-l plaquettes are
    impossible -- this sanity check should return zero. If it does, that
    confirms Rule B *cannot* reduce to the scalar 2D-grid k=2 recursion
    by ignoring cross-l plaquettes.
    """
    within_cycles = [c for c in cycles
                     if len(set(l_per_node[v] for v in c)) == 1]
    cross_cycles = [c for c in cycles
                    if len(set(l_per_node[v] for v in c)) > 1]
    return {
        "n_within_l_plaquettes": len(within_cycles),
        "n_cross_l_plaquettes": len(cross_cycles),
        "total_plaquettes": len(cycles),
        "within_l_fraction": (
            len(within_cycles) / len(cycles) if cycles else 0.0
        ),
    }


# =============================================================================
# §10 Main driver
# =============================================================================

def main() -> None:
    print("=" * 76)
    print("XCWG Migdal-Kadanoff RG on Dirac-Rule-B graph")
    print("=" * 76)

    results: Dict = {
        "sprint": "XCWG MK block-spin RG",
        "graph_family": "Dirac Rule B (Paper 29 §RH-C)",
        "recursion": "I_1(beta_eff) / I_0(beta_eff) = [I_1(beta) / I_0(beta)]^k",
    }

    # ----------------------------------------------------------------
    # Part 1: 2D scalar baseline sanity check
    # ----------------------------------------------------------------
    print("\n" + "-" * 76)
    print("Part 1: 2D scalar baseline (k=2) -- reproduce rg_blockspin_pilot")
    print("-" * 76)
    scalar_sanity = sanity_check_scalar_2D()
    results["scalar_2D_baseline"] = scalar_sanity
    print(f"{'beta':>8s} {'beta_eff':>12s} {'beta^2/2 (sc)':>16s} "
          f"{'beta/2 (wc)':>14s} {'err_sc':>10s} {'err_wc':>10s}")
    for r in scalar_sanity["data"]:
        print(f"{r['beta']:8.3f} {r['beta_eff']:12.6f} "
              f"{r['strong_coupling']:16.6f} {r['weak_coupling']:14.6f} "
              f"{r['rel_err_strong']:10.2e} {r['rel_err_weak']:10.2e}")

    # ----------------------------------------------------------------
    # Part 2: Rule B graph properties + plaquette-edge incidence
    # ----------------------------------------------------------------
    for n_max in [2, 3]:
        print("\n" + "-" * 76)
        print(f"Part 2.{n_max}: Rule B at n_max = {n_max}")
        print("-" * 76)
        A, labels, deg, desc = build_dirac_s3_graph(n_max, "B")
        l_per_node = []
        for lab in labels:
            kappa = lab.kappa
            l_val = (-kappa - 1) if kappa < 0 else kappa
            l_per_node.append(l_val)

        V = int(A.shape[0])
        E = int(A.sum() // 2)
        print(f"V = {V}, E = {E}, description: {desc}")
        cycles = primitive_4_plaquettes(A)
        print(f"Number of 4-plaquettes: {len(cycles)}")

        edge_count, stats = plaquette_edge_incidence(cycles, V)
        print(f"Plaquettes per edge:")
        print(f"  Mean: {stats['mean_plaq_per_edge']:.3f}")
        print(f"  Min:  {stats['min_plaq_per_edge']}")
        print(f"  Max:  {stats['max_plaq_per_edge']}")
        print(f"  Median: {stats['median_plaq_per_edge']:.2f}")
        print(f"  Std:    {stats['std_plaq_per_edge']:.3f}")
        print(f"  Edges with >=1 plaq: {stats['total_edges_with_plaquette']}")

        # k_eff for global MK recursion (canonical choice: mean q)
        k_eff_global = stats["k_eff_mean"]
        print(f"\nCanonical k_eff (mean plaquettes per edge): {k_eff_global:.3f}")
        print(f"  (compare 2D-square: 2.0; 3D-cubic: 4.0;"
              f" 4D-hypercubic: 6.0)")

        # ------ Within-l restriction sanity check ------
        within_l = within_l_restriction_sanity(A, cycles, l_per_node)
        print(f"\nWithin-l restriction (sanity vs scalar 2D-grid):")
        print(f"  Within-l plaquettes: {within_l['n_within_l_plaquettes']}")
        print(f"  Cross-l plaquettes:  {within_l['n_cross_l_plaquettes']}")
        print(f"  --> Rule B has NO within-l plaquettes, so cannot reduce to "
              f"scalar k=2")

        # ------ beta_eff(beta) for the global k_eff ------
        beta_grid = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 50.0]
        beta_eff_table = verify_limits(k_eff_global, beta_grid)

        # ------ Fixed-point search ------
        fp_global = find_fixed_points(k_eff_global)
        print(f"\nFixed-point analysis (k_eff = {k_eff_global:.3f}):")
        print(f"  Sign changes of (beta_eff - beta) on [1e-3, 100]: "
              f"{fp_global['n_sign_changes']}")
        if fp_global['fixed_points']:
            for fp in fp_global['fixed_points']:
                print(f"  Fixed point at beta* = {fp['beta_star']:.6f}, "
                      f"d_beta_eff/d_beta = {fp['d_beta_eff_d_beta']:.4f}, "
                      f"stable = {fp['stable']}")
        else:
            print(f"  No non-trivial fixed point in (0, infinity)")
        print(f"  g(beta_lo=1e-3) = {fp_global['g_at_extremes']['beta_lo']:.6e}")
        print(f"  g(beta_hi=100)  = {fp_global['g_at_extremes']['beta_hi']:.6e}")

        # ------ Compare with alternative k_eff choices ------
        alt_k_results = {}
        for k_name, k_val in [
            ("k_q_minus_1", k_eff_global - 1.0),
            ("k_geometric", stats["geometric_mean_plaq_per_edge"]),
            ("k_max", float(stats["k_eff_max"])),
            ("k_2 (scalar)", 2.0),
        ]:
            fp_alt = find_fixed_points(k_val)
            alt_k_results[k_name] = {
                "k_eff": k_val,
                "n_sign_changes": fp_alt["n_sign_changes"],
                "fixed_points": fp_alt["fixed_points"],
                "g_lo": fp_alt["g_at_extremes"]["beta_lo"],
                "g_hi": fp_alt["g_at_extremes"]["beta_hi"],
            }

        # ------ Edgewise per-edge MK at beta = 1 ------
        edge_mk = edgewise_mk_step(A, cycles, beta=1.0)

        # ------ Iterate the global recursion ------
        beta_0_list = [0.5, 1.0, 2.0, 5.0, 10.0]
        flow_iter = {b0: iterate_mk(b0, k_eff_global, n_steps=8)
                     for b0 in beta_0_list}

        # ------ Beta_eff table print ------
        print(f"\nbeta_eff(beta) at k_eff = {k_eff_global:.3f}:")
        print(f"{'beta':>8s} {'beta_eff':>14s} {'sc lim':>14s} "
              f"{'wc lim':>14s}  {'err_sc':>10s}  {'err_wc':>10s}")
        for r in beta_eff_table:
            print(f"{r['beta']:8.3f} {r['beta_eff']:14.6e} "
                  f"{r['strong_coupling']:14.6e} {r['weak_coupling']:14.6e} "
                  f"{r['rel_err_strong']:10.2e} {r['rel_err_weak']:10.2e}")

        # ------ Flow iteration print ------
        print(f"\nIterated MK flow (k_eff = {k_eff_global:.3f}), 8 steps:")
        header = f"{'step':>4s}  " + "  ".join(
            f"{b0:>10.2f}" for b0 in beta_0_list)
        print(header)
        for step in range(9):
            row_vals = []
            for b0 in beta_0_list:
                if step < len(flow_iter[b0]):
                    row_vals.append(f"{flow_iter[b0][step]:10.6e}")
                else:
                    row_vals.append(f"{'-':>10s}")
            print(f"{step:>4d}  " + "  ".join(row_vals))

        results[f"n_max_{n_max}"] = {
            "V": V,
            "E": E,
            "description": desc,
            "n_plaquettes_4": len(cycles),
            "edge_plaquette_stats": stats,
            "within_l_restriction": within_l,
            "k_eff_canonical": k_eff_global,
            "k_eff_alternatives": alt_k_results,
            "beta_eff_table": beta_eff_table,
            "fixed_points_global": fp_global,
            "edgewise_mk_at_beta_1": edge_mk,
            "flow_iter": {str(b0): flow_iter[b0] for b0 in beta_0_list},
        }

    # ----------------------------------------------------------------
    # Part 3: Cross-check with 3D-cubic-equivalent k=4 (heuristic)
    # ----------------------------------------------------------------
    print("\n" + "-" * 76)
    print("Part 3: 3D-cubic-equivalent (k=4) recursion (heuristic comparison)")
    print("-" * 76)
    k_3d = 4.0
    beta_grid_3d = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 50.0]
    table_3d = verify_limits(k_3d, beta_grid_3d)
    fp_3d = find_fixed_points(k_3d)
    print(f"k = {k_3d} (textbook 3D cubic):")
    print(f"  Sign changes: {fp_3d['n_sign_changes']}")
    if fp_3d['fixed_points']:
        for fp in fp_3d['fixed_points']:
            print(f"  Fixed point at beta* = {fp['beta_star']:.6f}, "
                  f"stable = {fp['stable']}")
    else:
        print(f"  No non-trivial fixed point (MK known not to capture 3D "
              f"deconfinement)")
    results["k_3d_baseline"] = {
        "k_eff": k_3d,
        "beta_eff_table": table_3d,
        "fixed_points": fp_3d,
    }

    # ----------------------------------------------------------------
    # Write JSON
    # ----------------------------------------------------------------
    out_path = os.path.join(_HERE, "data", "xcwg_mk_blockspin_rule_b.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=float)
    print("\n" + "=" * 76)
    print(f"Wrote {out_path}")
    print("=" * 76)

    # ----------------------------------------------------------------
    # Final verdict line
    # ----------------------------------------------------------------
    k_n2 = results["n_max_2"]["k_eff_canonical"]
    k_n3 = results["n_max_3"]["k_eff_canonical"]
    n_fp_n2 = results["n_max_2"]["fixed_points_global"]["n_sign_changes"]
    n_fp_n3 = results["n_max_3"]["fixed_points_global"]["n_sign_changes"]
    print(f"\nFINAL VERDICT:")
    print(f"  k_eff (Rule B, n_max=2): {k_n2:.3f}  "
          f"(non-trivial fixed points: {n_fp_n2})")
    print(f"  k_eff (Rule B, n_max=3): {k_n3:.3f}  "
          f"(non-trivial fixed points: {n_fp_n3})")
    print(f"  k_eff (2D scalar):       2.000  (fixed points: 0, beta -> 0)")
    print(f"  k_eff (3D cubic):        4.000  (fixed points: 0, MK fails 3D)")

    if n_fp_n2 == 0 and n_fp_n3 == 0:
        print(f"  Verdict: Rule B MK recursion has NO non-trivial fixed point")
        print(f"           (qualitatively same as scalar/3D MK: beta flows to 0)")
        print(f"           But k_eff much LARGER than scalar -- faster RG flow")
    else:
        print(f"  Verdict: Rule B MK recursion HAS a non-trivial fixed point.")
        print(f"           This would be QUALITATIVELY DIFFERENT from 2D scalar.")


if __name__ == "__main__":
    main()
