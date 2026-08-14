"""Poly-1: can the two-centre engine be USED at an arbitrary orientation?

Poly-0 says 86.3% of water's sum|g| lives in entries touching at most two
centres. That share is only reachable if the engine can be applied to a pair
whose axis is NOT the z-axis -- and water's O-H bonds sit at +-52.25 degrees.

The engine is inherently axial: prolate spheroidal coordinates put the two
centres on the z-axis, and every closed form from this arc is written in
(xi, eta, phi) about that axis. The claim under test is the standard rotation
route, already used for the ONE-body cross-centre integral in shibuya_wulfman:

    evaluate in the frame where the pair axis IS z, then rotate each orbital
    index back with a Wigner-D matrix.

For real Cartesian l=1 functions the rotation is just the ordinary 3x3 rotation,
so the claim is concretely

    g_lab = (D (x) D (x) D (x) D) . g_axis

with D block-diagonal over l. If that holds the 86.3% is reachable TODAY and the
only thing water needs is the 3-centre class. If it fails, water is blocked twice.

Tested against the Gaussian reference, not against the symbolic engine: build the
SAME two-centre system twice -- once with its axis on z, once rotated into a
general orientation -- and check the two tensors are related by D^{(x)4}. This
isolates the rotation claim from every other moving part.

Failure mode being probed: an l=1 shell is a genuine 3-dim rep, so if the
index bookkeeping (the (x,y,z) ordering, or the normalisation) were off, the
rotated tensor would NOT reproduce -- a wrong-basis bug cannot hide here.

Run from repo root:  python debug/poly1_rotation_gate.py
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

CART = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]        # x, y, z


def rot_matrix(beta: float, gamma: float) -> np.ndarray:
    """Rotate about y by beta, then about z by gamma."""
    cb, sb, cg, sg = np.cos(beta), np.sin(beta), np.cos(gamma), np.sin(gamma)
    Ry = np.array([[cb, 0, sb], [0, 1, 0], [-sb, 0, cb]])
    Rz = np.array([[cg, -sg, 0], [sg, cg, 0], [0, 0, 1]])
    return Rz @ Ry


def make(shapes, A, B, zA, zB):
    """s + p shell on A, s on B.  Returns (orbs, l_of_orb)."""
    orbs, ls = [], []
    orbs.append(GE.sto_shape_basis(A, "1s", zA, shapes, (0, 0, 0))); ls.append(0)
    for lmn in CART:
        orbs.append(GE.sto_shape_basis(A, "2p", zA, shapes, lmn)); ls.append(1)
    orbs.append(GE.sto_shape_basis(B, "1s", zB, shapes, (0, 0, 0))); ls.append(0)
    return orbs, ls


def block_D(ls, D3):
    """Per-orbital rotation matrix, block diagonal: 1 on l=0, D3 on each l=1."""
    M = len(ls)
    U = np.zeros((M, M))
    i = 0
    while i < M:
        if ls[i] == 0:
            U[i, i] = 1.0
            i += 1
        else:
            U[i:i + 3, i:i + 3] = D3
            i += 3
    return U


def main() -> None:
    print("Poly-1 -- does the axial engine survive an arbitrary orientation?\n")
    sh = {}
    for kind, (l, n_r) in (("1s", (0, 1)), ("2p", (1, 2))):
        a, d, q = GE.fit_sto_shape(l, n_r, n_gauss=8)
        sh[kind] = (a, d)

    R, zA, zB = 1.9, 1.6, 1.0
    A0 = np.array([0., 0., 0.])
    B0 = np.array([0., 0., R])                    # axis ON z

    orbs0, ls = make(sh, A0, B0, zA, zB)
    _, _, g0 = GE.integral_set_md(orbs0, [(A0, 3.0), (B0, 1.0)])

    print("  reference frame: pair axis on z, s+p on A, s on B"
          f"   (M = {len(orbs0)})")
    print(f"  {'beta':>7} {'gamma':>7}   {'max|g_lab - D^4 g_axis|':>26}"
          f"   {'max|g_lab|':>11}   verdict")
    print("  " + "-" * 74)

    worst = 0.0
    for beta, gamma in ((0.0, 0.0), (np.radians(52.25), 0.0),
                        (np.radians(52.25), np.radians(90.0)),
                        (np.radians(104.5), np.radians(37.0)),
                        (np.radians(-52.25), np.radians(180.0))):
        D3 = rot_matrix(beta, gamma)
        B1 = D3 @ B0
        orbs1, _ = make(sh, A0, B1, zA, zB)
        _, _, g1 = GE.integral_set_md(orbs1, [(A0, 3.0), (B1, 1.0)])

        U = block_D(ls, D3)
        # g_lab = U U U U . g_axis  (four index rotations)
        pred = np.einsum('pi,qj,rk,sl,ijkl->pqrs', U, U, U, U, g0, optimize=True)
        err = float(np.abs(g1 - pred).max())
        worst = max(worst, err)
        ok = "ok" if err < 1e-10 else "FAILS"
        print(f"  {np.degrees(beta):7.2f} {np.degrees(gamma):7.2f}"
              f"   {err:>26.3e}   {np.abs(g1).max():>11.6f}   {ok}")

    print(f"\n  worst {worst:.2e}")
    if worst < 1e-10:
        print("  => the axial engine transports to ANY orientation by a Wigner-D")
        print("     rotation on each index. Water's <=2-centre block (86.3% of")
        print("     sum|g|, 1981 of 2401 entries) is reachable with what exists.")
    else:
        print("  => rotation route BROKEN; water is blocked on the 2-centre part too.")


if __name__ == "__main__":
    main()
