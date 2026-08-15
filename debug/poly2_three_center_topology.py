"""Poly-2: the 3-centre block is TWO problems, not one -- which one is expensive?

Poly-0 (2026-08-13) measured the gap: water's 420 three-centre entries are 13.7%
of sum|g| but cost 2.35 Ha to drop, so there is no truncation story. It stopped
there, treating "3-centre" as one monolithic capability.

It is not. With three distinct centres among four orbital indices the centre
multiset is forced to be {X,X,Y,Z}, and there are exactly TWO topologies:

  T1  (XX|YZ)   the doubled centre sits entirely inside ONE density.
                => that density is ONE-CENTRE, so its Coulomb potential is
                   closed-form ALREADY (increment 1 / two_center_eri.V_L, exact).
                   The integral collapses to a THREE-CENTRE ONE-ELECTRON
                   integral of a two-centre density against a known multipole
                   potential. One electron, not two.

  T2  (XY|XZ)   the doubled centre is split ACROSS the densities.
                => both densities are two-centre, on two different axes sharing
                   a vertex. Neither has a closed-form potential; the Neumann
                   route needs both electrons in ONE (xi,eta) system and there is
                   no such system for a triangle. This is the genuine wall.

T1 is strictly the easier tier -- it reduces the electron count. So the useful
question Poly-0 did not ask is: DOES CLOSING T1 ALONE BUY A USABLE WATER?

Counting is settled in advance and is not the interesting part. Of the 12 index
arrangements of {X,X,Y,Z}, 4 are T1 and 8 are T2, giving 140 / 280 for water.
T2 is twice as NUMEROUS. The gate is energy, and the two need not agree --
Poly-0's own headline was precisely that norm share and energy cost disagree.

PRE-REGISTERED PREDICTION (written before running):
  T1 carries the majority of the 2.35 Ha despite being the minority of entries,
  because T1 contains the large one-centre O-core densities (1s^2, 2s^2) while
  every T2 term is a product of two small two-centre OVERLAP densities.

PRE-REGISTERED GATE:
  GO for a T1-first build   if dropping T2 ALONE costs < 0.1 Ha
                               (~60x chemical accuracy; T1 then is the problem)
  BORDERLINE                   0.1 - 0.5 Ha
  STOP / rescope            if dropping T2 ALONE costs > 0.5 Ha
                               (closing the easier tier does not buy water;
                                native water needs the full 3-centre 2-electron
                                problem, and that is a research programme)

Additivity is NOT assumed. |T1| + |T2| vs |both| is reported, because a large
mismatch means the blocks cancel and no per-block truncation story survives.

Same methodology as Poly-0: one McMurchie-Davidson reference tensor, partitioned.
No fit error enters the comparison. H2 is the null control (no 3-centre entries
at all, both drops must be exactly zero); BeH2 is the linear control, separating
"three centres" from "bent".

Run from repo root:  python debug/poly2_three_center_topology.py
"""

from __future__ import annotations

import sys
from itertools import product
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from geovac import noci_engine as GE  # noqa: E402

# reuse Poly-0's geometry/basis verbatim so the numbers are directly comparable
from poly0_three_center_burden import (  # noqa: E402
    _shapes,
    build,
    solve,
    vnn,
)

CHEM_ACC = 1.6e-3


def topology_masks(g: np.ndarray, cen: list[int], M: int):
    """Split the 3-centre block into T1 = (XX|YZ) and T2 = (XY|XZ).

    Returns (t1, t2, three) boolean masks over the (p,q,r,s) tensor.
    """
    t1 = np.zeros_like(g, dtype=bool)
    t2 = np.zeros_like(g, dtype=bool)
    three = np.zeros_like(g, dtype=bool)

    for p, q, r, s in product(range(M), repeat=4):
        if len({cen[p], cen[q], cen[r], cen[s]}) != 3:
            continue
        three[p, q, r, s] = True
        bra = {cen[p], cen[q]}
        ket = {cen[r], cen[s]}
        # one density one-centre => T1; both two-centre => T2
        if len(bra) == 1 or len(ket) == 1:
            t1[p, q, r, s] = True
        else:
            t2[p, q, r, s] = True
    return t1, t2, three


def exponent_symmetry_split(g, cen, zeta, t1_mask, M):
    """Within T1 = (XX|YZ), does the two-centre density Y-Z carry EQUAL exponents?

    The arc's termination criterion (build plan section 8.5, EQ1) is
    q = (alpha - beta) R / 2, and the Neumann tau-sum terminates EXACTLY when
    q = 0, i.e. when the two centres of the density carry the same exponent.
    For water the H1-H2 density is exactly that case (both hydrogens, zeta = 1).

    SCOPE: EQ1 was established for the EXCHANGE class, where both electrons live
    in one (xi,eta) system. T1 is a different integral and the criterion is not
    known to transfer. Reported as a structural annotation, not as a claim.
    """
    sym = np.zeros_like(g, dtype=bool)
    asym = np.zeros_like(g, dtype=bool)
    for p, q, r, s in product(range(M), repeat=4):
        if not t1_mask[p, q, r, s]:
            continue
        bra, ket = {cen[p], cen[q]}, {cen[r], cen[s]}
        # the two-centre side is whichever has two distinct centres
        pair = (r, s) if len(bra) == 1 else (p, q)
        if abs(zeta[pair[0]] - zeta[pair[1]]) < 1e-12:
            sym[p, q, r, s] = True
        else:
            asym[p, q, r, s] = True
    return sym, asym


def one_body_three_centre(orbs, cen, nuc, M):
    """The ONE-body three-centre burden, which Poly-0 never measured.

    V_ne matrix elements <chi_i| -Z_A/r_A |chi_j> touch three distinct centres
    whenever c(i), c(j) and A are all different. That integral is a genuine
    three-centre ONE-electron object and it is NOT in the repo either:
    shibuya_wulfman.py computes <chi^A|-Z_B/r_B|chi^A> -- both orbitals on the
    SAME centre A -- which is a TWO-centre integral.

    Returns (h_full, h_dropped, n_entries, weight).
    """
    h_full = np.zeros((M, M))
    h_drop = np.zeros((M, M))
    n3 = 0
    w3 = 0.0
    for i in range(M):
        for j in range(i, M):
            t = GE.kinetic_md(orbs[i], orbs[j])
            full = t
            keep = t
            for a, (pos, z) in enumerate(nuc):
                v = GE.nuclear_md(orbs[i], orbs[j], pos, z)
                full += v
                if len({cen[i], cen[j], a}) == 3:
                    n3 += 1 if i == j else 2
                    w3 += abs(v) if i == j else 2 * abs(v)
                else:
                    keep += v
            h_full[i, j] = h_full[j, i] = full
            h_drop[i, j] = h_drop[j, i] = keep
    return h_full, h_drop, n3, w3


def drop_cost(S, h, g, mask, n_elec, enn, e_full):
    gd = g.copy()
    gd[mask] = 0.0
    return solve(S, h, gd, n_elec, enn) - e_full


def report(system: str, shapes, zetas: dict[str, list[float]]) -> None:
    orbs, cen, nuc, n_elec, label = build(system, shapes)
    M = len(orbs)
    S, h, g = GE.integral_set_md(orbs, nuc)
    enn = vnn(nuc)
    e_full = solve(S, h, g, n_elec, enn)
    tot = np.abs(g).sum()

    t1, t2, three = topology_masks(g, cen, M)

    print(f"\n{'='*72}\n{label}   M = {M}, {n_elec} electrons")
    print(f"  E (all integrals) = {e_full:.8f} Ha")

    n3 = int(three.sum())
    if n3 == 0:
        c1 = drop_cost(S, h, g, t1, n_elec, enn, e_full)
        c2 = drop_cost(S, h, g, t2, n_elec, enn, e_full)
        print(f"  no 3-centre entries at all -- null control")
        print(f"  drop T1: {c1:+.8f} Ha   drop T2: {c2:+.8f} Ha   (both must be 0)")
        return

    print(f"\n  {'block':<28}{'entries':>9}{'sum|g|':>11}{'share3':>9}"
          f"{'drop cost / Ha':>17}{'x chem':>10}")
    print("  " + "-" * 84)

    w3 = float(np.abs(g[three]).sum())
    rows = [("T1  (XX|YZ)  1-ctr x 2-ctr", t1),
            ("T2  (XY|XZ)  2-ctr x 2-ctr", t2),
            ("both (= Poly-0 3-centre)", three)]
    costs = {}
    for name, mask in rows:
        n = int(mask.sum())
        w = float(np.abs(g[mask]).sum())
        c = drop_cost(S, h, g, mask, n_elec, enn, e_full)
        costs[name] = c
        print(f"  {name:<28}{n:>9}{w:>11.4f}{100*w/w3:>8.1f}%"
              f"{c:>+17.8f}{abs(c)/CHEM_ACC:>10.0f}")

    c1 = costs["T1  (XX|YZ)  1-ctr x 2-ctr"]
    c2 = costs["T2  (XY|XZ)  2-ctr x 2-ctr"]
    cb = costs["both (= Poly-0 3-centre)"]
    print(f"\n  additivity check: T1 + T2 = {c1+c2:+.8f} vs both = {cb:+.8f} Ha"
          f"   (mismatch {abs(c1+c2-cb):.4f})")

    # structural annotation: exponent symmetry of the T1 two-centre density
    zl = zetas[system]
    sym, asym = exponent_symmetry_split(g, cen, zl, t1, M)
    ws, wa = float(np.abs(g[sym]).sum()), float(np.abs(g[asym]).sum())
    print(f"\n  T1 sub-split by exponent symmetry of the two-centre density:")
    print(f"    equal exponents  (q = 0, tau-sum would terminate) : "
          f"{int(sym.sum()):>5} entries, sum|g| {ws:>9.4f}")
    print(f"    unequal exponents (q != 0, infinite but convergent): "
          f"{int(asym.sum()):>5} entries, sum|g| {wa:>9.4f}")

    # the ONE-body three-centre burden -- never measured by Poly-0
    hf, hd, nh3, wh3 = one_body_three_centre(orbs, cen, nuc, M)
    assert np.abs(hf - h).max() < 1e-10, "one-body rebuild disagrees with integral_set_md"
    e_h_drop = solve(S, hd, g, n_elec, enn) - e_full
    print(f"\n  ONE-body 3-centre block (V_ne, c(i) != c(j) != nucleus):")
    print(f"    {nh3} of {M*M} h-entries carry it, sum|v| = {wh3:.4f}")
    print(f"    cost of dropping    : {e_h_drop:+.8f} Ha"
          f"   = {abs(e_h_drop)/CHEM_ACC:.0f}x chemical accuracy")

    # the gate
    print(f"\n  GATE (drop T2 alone) : {abs(c2):.6f} Ha = {abs(c2)/CHEM_ACC:.0f}x chemical accuracy")
    if abs(c2) < 0.1:
        print("    -> GO for a T1-first build")
    elif abs(c2) < 0.5:
        print("    -> BORDERLINE")
    else:
        print("    -> STOP: closing the easier tier does not buy a usable water")


def main() -> None:
    print("Poly-2 -- the 3-centre block split by topology (T1 reducible, T2 not)")
    sh = _shapes()
    # orbital exponents in the SAME order as build() appends them
    zetas = {
        "H2":   [1.0, 1.0],
        "BeH2": [3.68, 0.96, 0.96, 0.96, 0.96, 1.0, 1.0],
        "H2O":  [7.66, 2.25, 2.23, 2.23, 2.23, 1.0, 1.0],
    }
    for sysname in ("H2", "BeH2", "H2O"):
        report(sysname, sh, zetas)
    print("\n" + "=" * 72)
    print("Read the GATE line for H2O. Counting was settled in advance"
          " (140 T1 / 280 T2);\nthe energy split is the result.")


if __name__ == "__main__":
    main()
