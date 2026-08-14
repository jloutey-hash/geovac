"""Poly-0: how much of water is OUT of reach of a two-centre engine?

Diagnostic before engineering. This arc closed the two-centre ERI in exact form,
all four classes. Water has three nuclei, so its ERI tensor splits into entries
whose four orbital indices touch one, two, or three distinct centres. The first
two the engine can do exactly TODAY (any geometry, by rotating the pair axis onto
z with the Wigner-D machinery already in shibuya_wulfman). The third it cannot do
at all: prolate spheroidal coordinates have exactly two foci, and a Slater
function has no product theorem to move a third centre.

So before proposing any build, measure the size of the gap. Three numbers decide
whether "native water" is a project or a non-starter:

  COUNT     how many 3-centre entries are there? (the build burden)
  NORM      what fraction of sum|g| do they carry? (the naive proxy)
  ENERGY    what does dropping them cost, in Ha? (the number that decides)

The third is the gate, and the first two are known to mislead: an integral class
can be a small fraction of the norm and still be worth 100x chemical accuracy,
because the energy is a difference of large cancelling terms.

BeH2 is included as the linear control -- it also has 3 nuclei, so it isolates
"three centres" from "bent". H2 is the null control: no 3-centre entries at all,
so its row must show exactly zero cost.

Reference integrals are the McMurchie-Davidson Gaussian path (exact for the fitted
basis); the partition is over the SAME tensor, so no fit error enters the
comparison at all -- the dropped-energy figure is a clean statement about the
partition, not about the basis.

Run from repo root:  python debug/poly0_three_center_burden.py
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

BOHR_OH = 1.809          # geovac/molecular_spec.py bond table
ANG_HOH = 104.5          # experimental
BOHR_BEH = 2.54


def _shapes():
    out = {}
    for kind, (l, n_r) in (("1s", (0, 1)), ("2s", (0, 2)), ("2p", (1, 2))):
        a, d, q = GE.fit_sto_shape(l, n_r, n_gauss=8)
        out[kind] = (a, d)
        out[kind + "_q"] = q
    return out


def water_geometry():
    th = np.radians(ANG_HOH / 2.0)
    O = np.array([0.0, 0.0, 0.0])
    H1 = np.array([BOHR_OH * np.sin(th), 0.0, BOHR_OH * np.cos(th)])
    H2 = np.array([-BOHR_OH * np.sin(th), 0.0, BOHR_OH * np.cos(th)])
    return O, H1, H2


def build(system: str, shapes):
    """Return (orbs, centre_index_per_orb, nuclei, n_elec, label)."""
    if system == "H2":
        A, B = np.array([0., 0., 0.]), np.array([0., 0., 1.4])
        orbs = [GE.sto_shape_basis(A, "1s", 1.0, shapes, (0, 0, 0)),
                GE.sto_shape_basis(B, "1s", 1.0, shapes, (0, 0, 0))]
        return orbs, [0, 1], [(A, 1.0), (B, 1.0)], 2, "H2 (2 nuclei, null control)"

    if system == "BeH2":
        Be = np.array([0., 0., 0.])
        Ha = np.array([0., 0., BOHR_BEH])
        Hb = np.array([0., 0., -BOHR_BEH])
        zBe, zH = 3.68, 1.0                     # Slater-rule valence zetas
        orbs, cen = [], []
        orbs.append(GE.sto_shape_basis(Be, "1s", 3.68, shapes, (0, 0, 0))); cen.append(0)
        orbs.append(GE.sto_shape_basis(Be, "2s", 0.96, shapes, (0, 0, 0))); cen.append(0)
        for lmn in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            orbs.append(GE.sto_shape_basis(Be, "2p", 0.96, shapes, lmn)); cen.append(0)
        orbs.append(GE.sto_shape_basis(Ha, "1s", zH, shapes, (0, 0, 0))); cen.append(1)
        orbs.append(GE.sto_shape_basis(Hb, "1s", zH, shapes, (0, 0, 0))); cen.append(2)
        nuc = [(Be, 4.0), (Ha, 1.0), (Hb, 1.0)]
        return orbs, cen, nuc, 6, "BeH2 (3 nuclei, LINEAR control)"

    O, H1, H2 = water_geometry()
    orbs, cen = [], []
    orbs.append(GE.sto_shape_basis(O, "1s", 7.66, shapes, (0, 0, 0))); cen.append(0)
    orbs.append(GE.sto_shape_basis(O, "2s", 2.25, shapes, (0, 0, 0))); cen.append(0)
    for lmn in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
        orbs.append(GE.sto_shape_basis(O, "2p", 2.23, shapes, lmn)); cen.append(0)
    orbs.append(GE.sto_shape_basis(H1, "1s", 1.0, shapes, (0, 0, 0))); cen.append(1)
    orbs.append(GE.sto_shape_basis(H2, "1s", 1.0, shapes, (0, 0, 0))); cen.append(2)
    nuc = [(O, 8.0), (H1, 1.0), (H2, 1.0)]
    return orbs, cen, nuc, 10, "H2O (3 nuclei, BENT, C2v)"


def vnn(nuclei):
    e = 0.0
    for i in range(len(nuclei)):
        for j in range(i + 1, len(nuclei)):
            e += nuclei[i][1] * nuclei[j][1] / np.linalg.norm(nuclei[i][0] - nuclei[j][0])
    return e


def solve(S, h, g, n_elec, enn):
    X = GE.lowdin_orbitals(S)
    ht, gt = GE.transform_integrals(X, h, g)
    return GE.fci_ground(ht, gt, n_elec) + enn


def report(system: str, shapes) -> None:
    orbs, cen, nuc, n_elec, label = build(system, shapes)
    M = len(orbs)
    S, h, g = GE.integral_set_md(orbs, nuc)
    enn = vnn(nuc)

    # partition the ERI tensor by how many DISTINCT centres the 4 indices touch
    masks = {1: np.zeros_like(g, dtype=bool),
             2: np.zeros_like(g, dtype=bool),
             3: np.zeros_like(g, dtype=bool)}
    for p, q, r, s in product(range(M), repeat=4):
        k = len({cen[p], cen[q], cen[r], cen[s]})
        masks[k][p, q, r, s] = True

    tot = np.abs(g).sum()
    e_full = solve(S, h, g, n_elec, enn)

    print(f"\n{label}   M = {M} spatial orbitals, {n_elec} electrons")
    print(f"  {'centres':>8} {'entries':>9} {'sum|g|':>12} {'share':>8}")
    print("  " + "-" * 42)
    for k in (1, 2, 3):
        n = int(masks[k].sum())
        w = float(np.abs(g[masks[k]]).sum())
        print(f"  {k:>8} {n:>9} {w:>12.4f} {100*w/tot:>7.1f}%")

    # what the two-centre engine reaches today
    reach = masks[1] | masks[2]
    n_reach, n_gap = int(reach.sum()), int(masks[3].sum())
    print(f"\n  engine reaches (1+2 centre) : {n_reach:>7} entries"
          f"   {100*np.abs(g[reach]).sum()/tot:5.1f}% of sum|g|")
    print(f"  OUT of reach   (3 centre)   : {n_gap:>7} entries"
          f"   {100*np.abs(g[masks[3]]).sum()/tot:5.1f}% of sum|g|")

    # THE GATE: energy cost of simply not having the 3-centre block
    g_drop = g.copy()
    g_drop[masks[3]] = 0.0
    e_drop = solve(S, h, g_drop, n_elec, enn)
    err = e_drop - e_full
    print(f"\n  E (all integrals)      : {e_full:14.8f} Ha")
    print(f"  E (3-centre dropped)   : {e_drop:14.8f} Ha")
    print(f"  cost of dropping       : {err:+14.8f} Ha"
          f"   = {abs(err)/1.6e-3:8.1f} x chemical accuracy")


def main() -> None:
    print("Poly-0 -- the three-centre burden, measured before anything is built\n")
    sh = _shapes()
    print("  Slater shapes fitted by 8 Gaussians; <fit|STO> = "
          + ", ".join(f"{k}:{sh[k+'_q']:.6f}" for k in ("1s", "2s", "2p")))
    print("  (fit error is COMMON to both columns below, so it cancels in the")
    print("   dropped-energy figure -- that number is about the partition alone)")
    for sysname in ("H2", "BeH2", "H2O"):
        report(sysname, sh)
    print("\n  Read the last line of each block, not the percentages.")


if __name__ == "__main__":
    main()
