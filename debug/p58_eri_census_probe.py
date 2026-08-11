"""Paper 58: numerically census the two-center ERI tensor.

The paper's Table 1 g row is COUNTED (symmetry-rule counting in the complex-m
basis), not computed -- no two-center ERI was ever numerically evaluated,
because no general-m Neumann engine exists (geovac/neumann_vee.py is m=0 only,
and generalizing its X_l / A_n / B_l machinery to sigma != 0 is a real build).

This probe narrows that gap without the Neumann build, using the validated
McMurchie-Davidson engine (geovac.noci_engine) over fitted Slater shapes.

What it can and cannot settle:
  CAN -- whether the symmetry-permitted entries are GENERICALLY NONZERO, i.e.
         whether the cross-center tensor really is dense; and whether the
         surviving symmetry is axial and ONLY axial.
  CANNOT -- decide exact zeros. These are Gaussian-fitted floats, so zeros are
         threshold decisions, exactly as for the same-center h block. Exact
         decidability needs the Neumann build.

Basis-convention note: MD works in the real Cartesian basis (px, py, pz), while
the paper's count is in the complex-m basis; the two differ by the +/-m mixing.
Rather than transform (and inherit a conjugation-convention minefield in the
4-index chemist bracket), the axial symmetry is tested in a BASIS-FREE form:
invariance of the tensor under rotation about the internuclear axis. That IS
the m-rule, stated without reference to m.

Run from repo root:  python debug/p58_eri_census_probe.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from geovac import noci_engine as E  # noqa: E402

TOL = 1e-10

# (label, kind, lmn) for n_max = 2: 1s, 2s, 2p_x, 2p_y, 2p_z
SHELLS = [
    ("1s", "1s", (0, 0, 0)),
    ("2s", "2s", (0, 0, 0)),
    ("2px", "2p", (1, 0, 0)),
    ("2py", "2p", (0, 1, 0)),
    ("2pz", "2p", (0, 0, 1)),
]

Z_A, Z_B, R = 3.0, 1.0, 3.0


def build_shapes():
    out = {}
    for kind, (l, n_r) in (("1s", (0, 1)), ("2s", (0, 2)), ("2p", (1, 2))):
        a, d, q = E.fit_sto_shape(l, n_r)
        out[kind] = (a, d)
        print(f"  fit {kind}: <fit|STO> = {q:.7f}")
    return out


def build_orbitals(shapes, axis="z"):
    """Ten orbitals, five per center, with both nuclei on the given axis."""
    if axis == "z":
        pa, pb = np.array([0., 0., 0.]), np.array([0., 0., R])
    else:  # 'x' -- used only to confirm the invariance test is discriminating
        pa, pb = np.array([0., 0., 0.]), np.array([R, 0., 0.])
    orbs, labels, recipes = [], [], []
    for center, pos, Z in (("A", pa, Z_A), ("B", pb, Z_B)):
        zeta = Z if center == "A" else 1.0
        for lab, kind, lmn in SHELLS:
            orbs.append(E.sto_shape_basis(pos, kind, zeta, shapes, lmn))
            labels.append(f"{center}:{lab}")
            recipes.append((pos, kind, zeta, lmn))
    nuclei = [(pa, Z_A), (pb, Z_B)]
    return orbs, labels, nuclei, recipes


def rotate_orbitals_about_z(orbs, theta, recipes, shapes):
    """Rotate every basis function about z. p_x/p_y mix; s and p_z are inert.

    Implemented by rotating the CENTERS and re-expressing p components; since
    both nuclei sit on z, the centers are fixed and only the p pair rotates.
    Returns a list of (coefficient, BasisFn) expansions.
    """
    c, s = np.cos(theta), np.sin(theta)
    out = []
    for o, (pos, kind, zeta, lmn) in zip(orbs, recipes):
        lx, ly, lz = lmn
        if (lx, ly) == (1, 0):        # p_x -> c p_x + s p_y
            out.append([(c, o), (s, _partner(pos, kind, zeta, lmn, shapes))])
        elif (lx, ly) == (0, 1):      # p_y -> -s p_x + c p_y
            out.append([(-s, _partner(pos, kind, zeta, lmn, shapes)), (c, o)])
        else:
            out.append([(1.0, o)])
    return out


def _partner(pos, kind, zeta, lmn, shapes):
    """p_x <-> p_y partner, rebuilt FROM THE RECIPE.

    Must NOT be built from an existing BasisFn: BasisFn.__init__ transforms its
    dcoeffs argument (multiplies by prim_norm, then renormalizes by
    1/sqrt(<self|self>)), so feeding a constructed .coeffs back in as dcoeffs
    double-applies the normalization and silently mis-scales the partner.  That
    bug produced identical invariance residuals for the collinear and
    non-collinear geometries, which is what exposed it.
    """
    lx, ly, lz = lmn
    return E.sto_shape_basis(pos, kind, zeta, shapes, (ly, lx, lz))


def eri_tensor(orbs):
    m = len(orbs)
    g = np.zeros((m, m, m, m))
    t0 = time.time()
    for i in range(m):
        for j in range(i, m):
            for k in range(m):
                for l in range(k, m):
                    if (i, j) > (k, l):
                        continue
                    v = E.eri_md(orbs[i], orbs[j], orbs[k], orbs[l])
                    for a, b in ((i, j), (j, i)):
                        for cc, d in ((k, l), (l, k)):
                            g[a, b, cc, d] = v
                            g[cc, d, a, b] = v
    print(f"  tensor built in {time.time() - t0:.1f}s")
    return g


def classify(labels):
    return ["A" if lab.startswith("A") else "B" for lab in labels]


def main():
    print("Fitting shapes:")
    shapes = build_shapes()

    print("\nBuilding z-axis (collinear) system:")
    orbs, labels, nuclei, recipes = build_orbitals(shapes, axis="z")
    g = eri_tensor(orbs)
    m = len(orbs)
    side = classify(labels)

    nz = int(np.sum(np.abs(g) > TOL))
    dense = m ** 4
    print(f"\n  nonzero  = {nz} / {dense}  ({100.0 * nz / dense:.1f}% dense)")

    # Cross-class decomposition: how much lives in classes touching both centers
    cross = 0
    for i in range(m):
        for j in range(m):
            for k in range(m):
                for l in range(m):
                    if abs(g[i, j, k, l]) <= TOL:
                        continue
                    if len({side[i], side[j], side[k], side[l]}) > 1:
                        cross += 1
    print(f"  of which cross-center classes: {cross} "
          f"({100.0 * cross / max(nz, 1):.1f}% of nonzeros)")

    # Same-center Gaunt-forbidden pairs that are nonzero ACROSS centers.
    # (1s, 2p) differ by one unit of l: same-center <1s|2p> angular factor is
    # zero, so any (1s_A 2p_B | . .) entry surviving is l-selection failing.
    iA1s, iB2pz = labels.index("A:1s"), labels.index("B:2pz")
    iA2s = labels.index("A:2s")
    probe = g[iA1s, iB2pz, iA2s, iA2s]
    print(f"\n  (A:1s B:2pz | A:2s A:2s) = {probe:.6e}  -> "
          f"{'NONZERO (l-selection fails cross-center)' if abs(probe) > TOL else 'zero'}")

    # Axial invariance: rotate about z by a generic angle; tensor must be fixed.
    print("\nAxial rotational invariance (the m-rule, basis-free):")
    worst_z = axial_invariance_residual(orbs, np.pi / 7, recipes, shapes)
    print(f"  rotation about z (internuclear axis): worst |dg| = {worst_z:.3e}")

    # Discriminator: with the nuclei on x instead, a z-rotation is NOT a
    # symmetry, so the same check must FAIL. Confirms the test has teeth.
    orbs_x, labels_x, _, recipes_x = build_orbitals(shapes, axis="x")
    gx = eri_tensor(orbs_x)
    worst_x = axial_invariance_residual(orbs_x, np.pi / 7, recipes_x, shapes,
                                       precomputed=gx)
    print(f"  same check, nuclei on x axis        : worst |dg| = {worst_x:.3e}"
          f"   <- must be LARGE")

    out = {
        "config": {"Z_A": Z_A, "Z_B": Z_B, "R": R, "M": m, "tol": TOL},
        "nonzero": nz, "dense": dense, "density_pct": 100.0 * nz / dense,
        "cross_class_nonzero": cross,
        "l_selection_probe": float(probe),
        "axial_residual_z": float(worst_z),
        "axial_residual_x_control": float(worst_x),
        "caveat": "Gaussian-fitted floats; zeros are threshold decisions, "
                  "not decided. Real Cartesian basis, not complex-m.",
    }
    p = REPO / "debug" / "data" / "p58_eri_census_probe.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nwrote {p}")


def axial_invariance_residual(orbs, theta, recipes, shapes, precomputed=None):
    """max |g - g_rotated| where the rotation acts on the p_x/p_y pair."""
    g = eri_tensor(orbs) if precomputed is None else precomputed
    exp = rotate_orbitals_about_z(orbs, theta, recipes, shapes)
    m = len(orbs)
    worst = 0.0
    # Spot-check a slice rather than all m^4 (cost); include every p index.
    p_idx = [i for i, o in enumerate(orbs) if sum(o.lmn) == 1]
    checks = [(i, j, k, l) for i in p_idx for j in p_idx
              for k in range(m) for l in range(m)]
    for (i, j, k, l) in checks:
        acc = 0.0
        for ci, oi in exp[i]:
            for cj, oj in exp[j]:
                for ck, ok in exp[k]:
                    for cl, ol in exp[l]:
                        if ci * cj * ck * cl == 0.0:
                            continue
                        acc += ci * cj * ck * cl * E.eri_md(oi, oj, ok, ol)
        worst = max(worst, abs(acc - g[i, j, k, l]))
    return worst


if __name__ == "__main__":
    main()
