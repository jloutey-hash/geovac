"""Paper 58: numerical corroboration of the two-center ERI census (Table 1 g row).

The paper's g row is COUNTED -- a symmetry-rule census in the complex-m basis,
with no two-center ERI numerically evaluated, because geovac/neumann_vee.py is
m = 0 only and generalizing its X_l / A_n / B_l machinery to sigma != 0 is a
substantial build.  These tests narrow that gap using the validated
McMurchie-Davidson engine over fitted Slater shapes.

What this settles / does not settle:
  SETTLES  -- that the symmetry-PERMITTED entries are generically nonzero (the
              cross-center tensor really is dense, 29.8% measured vs 29.4%
              counted), and that the surviving symmetry is axial and ONLY axial.
  DOES NOT -- decide exact zeros.  Gaussian-fitted floats, so zeros here are
              threshold decisions.  Exact decidability still needs the Neumann
              build; the g row stays COUNTED rather than becoming MEASURED.

Basis convention: MD is real Cartesian (px, py, pz); the paper counts in the
complex-m basis, and the two differ by the +/-m mixing -- which is why 29.8 and
29.4 are agreement rather than a discrepancy.  The axial symmetry is therefore
tested in a BASIS-FREE form: invariance of the tensor under rotation about the
internuclear axis.  That is the m-rule stated without reference to m.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("scipy", reason="engine needs scipy.special.hyp1f1")

from geovac import noci_engine as E  # noqa: E402

TOL = 1e-10
Z_A, Z_B, R = 3.0, 1.0, 3.0
SHELLS = [("1s", (0, 0, 0)), ("2s", (0, 0, 0)), ("2p", (1, 0, 0)),
          ("2p", (0, 1, 0)), ("2p", (0, 0, 1))]


def _shapes():
    return {k: (a, d) for k, (a, d, _q) in
            ((kind, E.fit_sto_shape(l, n_r)) for kind, (l, n_r) in
             (("1s", (0, 1)), ("2s", (0, 2)), ("2p", (1, 2))))}


def _system(shapes, axis="z"):
    pa = np.array([0.0, 0.0, 0.0])
    pb = np.array([0.0, 0.0, R]) if axis == "z" else np.array([R, 0.0, 0.0])
    orbs, recipes, side = [], [], []
    for tag, pos, zeta in (("A", pa, Z_A), ("B", pb, 1.0)):
        for kind, lmn in SHELLS:
            orbs.append(E.sto_shape_basis(pos, kind, zeta, shapes, lmn))
            recipes.append((pos, kind, zeta, lmn))
            side.append(tag)
    return orbs, recipes, side


def _partner(pos, kind, zeta, lmn, shapes):
    """p_x <-> p_y partner, rebuilt FROM THE RECIPE, never from a BasisFn.

    BasisFn.__init__ transforms its dcoeffs argument (prim_norm, then
    renormalization), so feeding a constructed .coeffs back in double-applies
    normalization.  That bug made the collinear and non-collinear invariance
    residuals come out identical, which is how it was caught.
    """
    lx, ly, lz = lmn
    return E.sto_shape_basis(pos, kind, zeta, shapes, (ly, lx, lz))


def _eri_tensor(orbs):
    m = len(orbs)
    g = np.zeros((m, m, m, m))
    for i in range(m):
        for j in range(i, m):
            for k in range(m):
                for l in range(k, m):
                    if (i, j) > (k, l):
                        continue
                    v = E.eri_md(orbs[i], orbs[j], orbs[k], orbs[l])
                    for a, b in ((i, j), (j, i)):
                        for c, d in ((k, l), (l, k)):
                            g[a, b, c, d] = g[c, d, a, b] = v
    return g


def _rotation_expansion(orbs, recipes, shapes, theta):
    c, s = np.cos(theta), np.sin(theta)
    out = []
    for o, (pos, kind, zeta, lmn) in zip(orbs, recipes):
        lx, ly = lmn[0], lmn[1]
        p = _partner(pos, kind, zeta, lmn, shapes)
        if (lx, ly) == (1, 0):
            out.append([(c, o), (s, p)])
        elif (lx, ly) == (0, 1):
            out.append([(-s, p), (c, o)])
        else:
            out.append([(1.0, o)])
    return out


def _invariance_residual(orbs, recipes, shapes, theta, g):
    exp = _rotation_expansion(orbs, recipes, shapes, theta)
    m = len(orbs)
    p_idx = [i for i, (_p, _k, _z, lmn) in enumerate(recipes) if sum(lmn) == 1]
    worst = 0.0
    for i in p_idx:
        for j in p_idx:
            for k in range(m):
                for l in range(m):
                    acc = 0.0
                    for ci, oi in exp[i]:
                        for cj, oj in exp[j]:
                            for ck, ok in exp[k]:
                                for cl, ol in exp[l]:
                                    w = ci * cj * ck * cl
                                    if w == 0.0:
                                        continue
                                    acc += w * E.eri_md(oi, oj, ok, ol)
                    worst = max(worst, abs(acc - g[i, j, k, l]))
    return worst


@pytest.mark.slow
def test_paper58_eri_density_matches_symmetry_count():
    """Permitted ERI entries are generically nonzero: 29.8% vs 29.4% counted.

    The paper's g row asserts the cross-center tensor is dense.  If the
    symmetry-permitted entries were mostly accidentally zero, that claim would
    fail even with the counting correct.  They are not.
    """
    shapes = _shapes()
    orbs, _rec, side = _system(shapes, axis="z")
    g = _eri_tensor(orbs)
    m = len(orbs)

    nz = int(np.sum(np.abs(g) > TOL))
    density = 100.0 * nz / m ** 4
    assert 25.0 < density < 35.0, (
        f"measured ERI density {density:.1f}% is far from the paper's counted "
        "29.4%; either the count or this evaluation is wrong"
    )

    # Cross-center classes dominate, which is the paper's point: binding lives
    # where the tensor is dense.  (Counted broadly here -- any index on a
    # different center -- so this includes AA|BB, which the paper's 80% figure
    # for AA|AB + AB|AB + AB|BB excludes.)
    cross = sum(
        1 for i in range(m) for j in range(m) for k in range(m) for l in range(m)
        if abs(g[i, j, k, l]) > TOL and len({side[i], side[j], side[k], side[l]}) > 1
    )
    assert cross > 0.75 * nz, (
        f"cross-center classes are only {100.0 * cross / nz:.0f}% of nonzeros"
    )


@pytest.mark.slow
def test_paper58_eri_l_selection_fails_cross_center():
    """Thm. 1(ii) at the ERI level, at bond scale.

    (1s, 2p) differ by one unit of l, so the same-center angular factor
    vanishes.  Across centers the entry is large.
    """
    shapes = _shapes()
    orbs, _rec, _side = _system(shapes, axis="z")
    g = _eri_tensor(orbs)
    # index 0 = A:1s, index 9 = B:2pz, index 1 = A:2s
    val = g[0, 9, 1, 1]
    assert abs(val) > 1e-3, (
        f"(A:1s B:2pz | A:2s A:2s) = {val:.3e}; cross-center l-selection "
        "failure is the paper's Thm. 1(ii) and should be bond-scale"
    )


@pytest.mark.slow
def test_paper58_eri_axial_invariance_and_control():
    """The m-rule, basis-free -- with a control proving the test can fail.

    Collinear nuclei: rotation about the internuclear axis maps the basis into
    itself and fixes both centers, so the tensor must be invariant.  Move the
    nuclei off that axis and the SAME check must break, since the rotation is
    then not a symmetry.  Without the control this test could pass by being
    blind.
    """
    shapes = _shapes()
    theta = np.pi / 7  # generic angle: no accidental fixed point

    orbs, rec, _side = _system(shapes, axis="z")
    res_axial = _invariance_residual(orbs, rec, shapes, theta, _eri_tensor(orbs))

    orbs_x, rec_x, _s = _system(shapes, axis="x")
    res_ctrl = _invariance_residual(orbs_x, rec_x, shapes, theta,
                                    _eri_tensor(orbs_x))

    assert res_axial < 1e-12, (
        f"ERI tensor not invariant under rotation about the internuclear axis "
        f"(residual {res_axial:.3e}); the m-rule would be false"
    )
    assert res_ctrl > 1e-3, (
        f"control residual {res_ctrl:.3e} is small -- the invariance check is "
        "blind and the axial result above proves nothing"
    )
    assert res_ctrl / max(res_axial, 1e-18) > 1e6, (
        "insufficient discrimination between the axial case and the control"
    )
