"""Gates for geovac/sturmian_general_sw.py: the Shibuya-Wulfman metric between
Coulomb-Sturmians of any (n, l, m) in any orientation (m-seam rung 3, step 1,
PI-directed 2026-10-10; plan in debug/m_seam_rung3_scoping_memo.md).

Each test names the wrong answer it rejects.
"""
from __future__ import annotations

import mpmath as mp
import numpy as np
import pytest

from geovac import sturmian_general_sw as g
from geovac import sturmian_sigma_basis as sb


@pytest.fixture(autouse=True)
def _prec():
    mp.mp.dps = 20


def test_momentum_radials_are_the_hankel_transform() -> None:
    """f_nl is normalised and equals sqrt(2/pi) int R_nl j_l(pr) r^2 dr of the
    tracked position-space radial (sign included).
    Rejects: a wrong momentum normalisation or sign convention."""
    pv = mp.mpf('0.7')
    for (n, l) in [(2, 0), (3, 1), (3, 2)]:
        assert abs(mp.quad(lambda p: g.f_mom(n, l, p) ** 2 * p * p, [0, 1, mp.inf]) - 1) < 1e-15
        def R(r):
            return mp.mpf(sb.radial(n, l, 1.0, np.array([float(r)]))[0])
        hank = mp.sqrt(2 / mp.pi) * mp.quad(
            lambda r: R(r) * mp.sqrt(mp.pi / (2 * pv * r)) * mp.besselj(l + mp.mpf(1) / 2, pv * r) * r * r,
            [0, 5, 20, 60])
        assert abs(g.f_mom(n, l, pv) - hank) < 1e-12


def test_radial_closed_form_matches_quadrature_and_has_no_pi() -> None:
    """Residue closed form == independent mpmath quadrature (1e-12) on
    representative (n l, n' l', L) at D = 1.3, contains no pi, and the
    1s-1s case is the textbook (1 + D) e^{-D}.
    Rejects: a wrong residue/branch of h_L; a surviving pi."""
    import sympy as sp
    D = g._D
    assert sp.simplify(g.radial_I_exact(1, 0, 1, 0, 0) - (1 + D) * sp.exp(-D)) == 0
    for case in [(2, 1, 2, 1, 2), (3, 2, 2, 1, 1), (3, 1, 2, 0, 1), (3, 2, 3, 2, 4), (4, 3, 3, 2, 3)]:
        e = g.radial_I_exact(*case)
        assert not e.has(sp.pi)
        v = g.radial_I(*case, '1.3')
        q = g.radial_I_quad(*case, '1.3')
        assert abs(v - q) < 1e-10   # the quadrature is the less accurate side


def test_own_centre_block_is_identity() -> None:
    """d = 0: the metric restricted to one centre is the identity (Sturmian
    orthonormality in the SW metric).  Rejects: a wrong weight (p^2+k^2)/2k^2."""
    labs = g.labels(3, 2)
    M = g.sw_block(labs, labs, [0, 0, 0])
    for i in range(len(labs)):
        for j in range(len(labs)):
            assert abs(M[i, j] - (1 if i == j else 0)) < 1e-18


def test_sigma_cross_block_matches_prolate_route() -> None:
    """m = 0 cross-centre block at R = 1.5 equals the independent
    position-space prolate quadrature of sturmian_sigma_basis to 1e-8,
    including odd-l signs.  Rejects: a wrong plane-wave phase i^{L+la-lb}."""
    R = 1.5
    S_pr, _, bl = sb.build(3, 1.0, R)
    sig = [(n, l, 0, 'c') for (n, l) in bl]
    blk = g.sw_block(sig, sig, [0, 0, -R])
    nb = len(bl)
    err = max(abs(float(blk[i, j]) - S_pr[i, nb + j]) for i in range(nb) for j in range(nb))
    assert err < 1e-8


def test_axial_selection_and_pi_degeneracy() -> None:
    """d along z: different (m, cos/sin) never couple, and the p_x and p_y
    blocks are equal.  Rejects: an m-leak in the angular factor."""
    labs = g.labels(3, 2)
    M = g.sw_block(labs, labs, [0, 0, mp.mpf('1.3')])
    for i, a in enumerate(labs):
        for j, b in enumerate(labs):
            if (a[2], a[3]) != (b[2], b[3]):
                assert abs(M[i, j]) < 1e-18
    px = [i for i, a in enumerate(labs) if a[1:] == (1, 1, 'c')]
    py = [i for i, a in enumerate(labs) if a[1:] == (1, 1, 's')]
    for i, j in zip(px, py):
        for i2, j2 in zip(px, py):
            assert abs(M[i, i2] - M[j, j2]) < 1e-18


def test_rotation_covariance() -> None:
    """p_x-p_x along x equals p_z-p_z along z, and d_{xz}-p_z (m = 1 vs 0)
    along the (x+z)/sqrt2 diagonal is nonzero while it vanishes along z.
    Rejects: an angular factor that ignores the direction of d."""
    D = mp.mpf('1.1')
    px, pz = (2, 1, 1, 'c'), (2, 1, 0, 'c')
    a = g.sw_block([px], [px], [D, 0, 0])[0, 0]
    b = g.sw_block([pz], [pz], [0, 0, D])[0, 0]
    assert abs(a - b) < 1e-18
    dxz = (3, 2, 1, 'c')
    diag = g.sw_block([dxz], [pz], [D / mp.sqrt(2), 0, D / mp.sqrt(2)])[0, 0]
    axial = g.sw_block([dxz], [pz], [0, 0, D])[0, 0]
    assert abs(diag) > 1e-3 and abs(axial) < 1e-18


def test_transpose_symmetry_and_gram_positivity() -> None:
    """S_ab(d) = S_ba(-d)^T, and the global metric of a bent three-centre
    geometry (centre l <= 2, protons l <= 1) is positive semidefinite.
    Rejects: a non-Hermitian metric; a sign error that breaks positivity."""
    la, lb = g.labels(3, 2), g.labels(2, 1)
    d = [mp.mpf('0.6'), 0, mp.mpf('1.2')]
    A = g.sw_block(la, lb, d)
    B = g.sw_block(lb, la, [-v for v in d])
    assert max(abs(A[i, j] - B[j, i]) for i in range(len(la)) for j in range(len(lb))) < 1e-18
    R, half = mp.mpf(2), mp.mpf('0.9')
    cen = [([0, 0, 0], 1, la),
           ([R * mp.sin(half), 0, R * mp.cos(half)], 1, lb),
           ([R * mp.sin(half), 0, -R * mp.cos(half)], 1, lb)]
    S, _, _ = g.assemble(cen)
    ev = mp.eigsy(S)[0]
    assert min(ev[t] for t in range(S.rows)) > -1e-15
