"""Paper 2 sec:scope ("The labelling matters"; carried from the archived Paper 25 Sec. III.B): the base graph
under the rotation grouping (Paper 30 Sec. 4.3's) and a Hopf (chiral) grouping.

Record: docs/predictions/2026-10-31_hopf_grouping.md (HG, GO).  Written as a separate step after the result
(CLAUDE.md Sec. 9).  No debug/ import: the SO(4) generators come from tests/test_atlas_so4_commutant.py (rebuilt on
n_max shells by rebinding its state list), the graph from geovac.lattice, and the exact occupations are recomputed
here from Clebsch-Gordan coefficients.

Hopf group (n, mu): eigenspace of J-_z = (L_z - A_z)/2 in shell n, weight n.  Occupation of node (n, l, m) in it:
|<j, m - mu; j, mu | l, m>|^2, j = (n - 1)/2.  Base adjacency: off-diagonal part of M^T A M (the occupation rule,
which equals the quotient rule of Paper 30 Sec. 4.3 when M is 0/1).
"""
import importlib.util
import os

import numpy as np
import sympy as sp
from sympy.physics.quantum.cg import CG

from geovac.lattice import GeometricLattice

HERE = os.path.dirname(os.path.abspath(__file__))


def so4(n_max):
    spec = importlib.util.spec_from_file_location('so4_for_hg_%d' % n_max, os.path.join(HERE, 'test_atlas_so4_commutant.py'))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    m.NMAX = n_max
    m.STATES = [(n, l, mm) for n in range(1, n_max + 1) for l in range(n) for mm in range(-l, l + 1)]
    m.IDX = {s: i for i, s in enumerate(m.STATES)}
    m.DIM = len(m.STATES)
    return m


def graph(n_max, states):
    lat = GeometricLattice(max_n=n_max)
    assert list(lat.states) == states
    return sp.Matrix(np.asarray(lat.adjacency.toarray(), dtype=int).tolist())


def rotation_M(states):
    groups = sorted({(n, l) for (n, l, m) in states})
    M = sp.zeros(len(states), len(groups))
    for i, (n, l, m) in enumerate(states):
        M[i, groups.index((n, l))] = 1
    return groups, M


def hopf_M(states):
    groups = [(n, sp.Rational(k, 2) - sp.Rational(n - 1, 2)) for n in sorted({s[0] for s in states}) for k in range(0, 2 * n - 1, 2)]
    M = sp.zeros(len(states), len(groups))
    for i, (n, l, m) in enumerate(states):
        j = sp.Rational(n - 1, 2)
        for g, (ng, mu) in enumerate(groups):
            if ng != n or abs(m - mu) > j:
                continue
            M[i, g] = CG(j, m - mu, j, mu, l, m).doit() ** 2
    return groups, M


def base(M, A):
    W = M.T * A * M
    for i in range(W.shape[0]):
        W[i, i] = 0
    return W, sp.diag(*[sum(W.row(i)) for i in range(W.shape[0])]) - W


def zero_modes(Lb):
    return Lb.shape[0] - Lb.rank()


def test_rotation_grouping_reproduces_paper25_and_hopf_occupations_are_exact_projectors():
    """Control: the rotation grouping through this code gives Paper 30 Sec. 4.3's links and spectrum {0,0,0,1,3,6}; the
    exact CG^2 occupations equal the diagonals of the J-_z eigenprojectors (so the Hopf groups are the fibre
    factor's eigenspaces, not an arbitrary partition).

    Wrong answers rejected: a base rule that is not Paper 30 Sec. 4.3's, and occupations that are not the chiral groups'.
    """
    m = so4(3)
    A = graph(3, m.STATES)
    _, Mr = rotation_M(m.STATES)
    _, Lr = base(Mr, A)
    x = sp.Symbol('x')
    assert sp.expand(Lr.charpoly(x).as_expr() - x ** 3 * (x - 1) * (x - 3) * (x - 6)) == 0
    L, Av = m.generators()
    Jmz = (L[2] - Av[2]) / 2
    groups, Mh = hopf_M(m.STATES)
    w, v = np.linalg.eigh(Jmz)
    for g, (n, mu) in enumerate(groups):
        ids = [i for i, s in enumerate(m.STATES) if s[0] == n]
        sel = [k for k in range(len(w)) if abs(w[k] - float(mu)) < 1e-9 and np.abs(v[ids, k]).sum() > 1e-9]
        P = v[:, sel] @ v[:, sel].conj().T
        occ = np.array([float(Mh[i, g]) for i in range(len(m.STATES))])
        assert np.abs(np.real(np.diag(P)) - occ).max() < 1e-12, (n, mu)


def test_hopf_base_graph_is_connected_and_its_fibre_is_thin():
    """HG1-HG2: the Hopf base graph has one zero mode (rotation: one per l-block, 3 and 4), its inter-group weight
    is 29/3 of 13 at n_max = 3 and 2009/75 of 34 at n_max = 4 (rotation 5 and 14), and the totals are
    grouping-blind.

    Wrong answers rejected: a Hopf grouping that keeps the l-blocks apart, a fibre as thick as the rotation one,
    and a rule that changes the total edge weight.
    """
    for n_max, inter_expect, total, rot_inter, rot_zero in ((3, sp.Rational(29, 3), 13, 5, 3),
                                                           (4, sp.Rational(2009, 75), 34, 14, 4)):
        m = so4(n_max)
        A = graph(n_max, m.STATES)
        _, Mh = hopf_M(m.STATES)
        _, Mr = rotation_M(m.STATES)
        Wh, Lh = base(Mh, A)
        Wr, Lr = base(Mr, A)
        assert zero_modes(Lh) == 1 and zero_modes(Lr) == rot_zero
        inter_h = sum(Wh) / 2
        assert inter_h == inter_expect and inter_h > sp.Rational(total, 2), inter_h
        assert sum(Wr) / 2 == rot_inter
        assert sum(Mh.T * A * Mh) == sum(Mr.T * A * Mr) == 2 * total


def test_hopf_spectrum_is_irrational_where_rotation_is_integer_and_the_casimir_halves():
    """HG3 at n_max = 3: the Hopf base characteristic polynomial has irreducible quadratic factors (discriminants 40
    and 265), where rotation's spectrum is integer.  Forced: Tr J+^2 = B/2 (21 against 42).

    Wrong answers rejected: an all-rational Hopf spectrum, and a chiral Casimir sum equal to B.
    """
    m = so4(3)
    A = graph(3, m.STATES)
    _, Mh = hopf_M(m.STATES)
    _, Lh = base(Mh, A)
    x = sp.Symbol('x')
    facs = sp.factor_list(sp.Poly(Lh.charpoly(x).as_expr(), x))[1]
    quads = sorted(int(sp.discriminant(f.as_expr(), x)) for f, _ in facs if f.degree() == 2)
    assert quads == [40, 265], quads
    L, Av = m.generators()
    Jp = [(L[a] + Av[a]) / 2 for a in range(3)]
    trJ = float(np.real(np.trace(sum(G @ G for G in Jp))))
    trL = float(np.real(np.trace(sum(G @ G for G in L))))
    assert abs(trL - 42) < 1e-9 and abs(trJ - 21) < 1e-9


# ---- second route: the published parabolic-spherical coefficients (no Clebsch-Gordan routine) ----
def _poch(a, k):
    p = 1
    for i in range(k):
        p *= a + i
    return p


def _hyp3f2(a, b, kmax):
    s = sp.Integer(0)
    for k in range(kmax + 1):
        den = _poch(b[0], k) * _poch(b[1], k) * sp.factorial(k)
        if den == 0:
            break
        s += sp.Rational(_poch(a[0], k) * _poch(a[1], k) * _poch(a[2], k), den)
    return s


def _w2_book(n1, n2, m, l):
    """Square of Mardoyan et al. (arXiv:2310.17336) eq. (2.2.29), the Tarter 3F2 form of the
    parabolic (n1, n2, m) -> spherical (n, l, m) coefficient, n = n1 + n2 + |m| + 1."""
    am = abs(m)
    n = n1 + n2 + am + 1
    if not (am <= l <= n - 1):
        return sp.Integer(0)
    F = sp.factorial
    pref = sp.Rational((2 * l + 1) * F(n1 + am) * F(n2 + am) * F(l + am),
                       F(n1) * F(n2) * F(l - am) * F(n - l - 1) * F(n + l))
    H = _hyp3f2((-n2, -l + am, l + am + 1), (am + 1, -n + am + 1), n2)
    return pref * sp.Rational(F(n - am - 1), F(am)) ** 2 * H ** 2


def test_occupations_equal_the_published_parabolic_spherical_coefficients():
    """The Hopf-group occupations equal, exactly, the squared parabolic-spherical coefficients of the
    published closed form, regrouped by mu = (m - n1 + n2)/2: a second route with no CG routine.

    Wrong answers rejected: a label map that drops m (mu = (n2 - n1)/2), and mu -> -mu at fixed m (an L_z sign
    error).  (Swapping n1 and n2 is not a wrong answer: CG^2(j,a;j,m-a) = CG^2(j,m-a;j,a) when both spins equal j.)
    """
    for n_max in (2, 3):
        m = so4(n_max)
        groups, Mh = hopf_M(m.STATES)
        for mu_map, expect_equal in ((lambda n1, n2, mm: sp.Rational(mm - n1 + n2, 2), True),
                                     (lambda n1, n2, mm: sp.Rational(n2 - n1, 2), False),
                                     (lambda n1, n2, mm: -sp.Rational(mm - n1 + n2, 2), False)):
            book = {}
            for n in range(1, n_max + 1):
                for mm in range(-(n - 1), n):
                    for n1 in range(0, n - abs(mm)):
                        n2 = n - 1 - abs(mm) - n1
                        mu = mu_map(n1, n2, mm)
                        for l in range(abs(mm), n):
                            key = (n, l, mm, mu)
                            book[key] = book.get(key, 0) + _w2_book(n1, n2, mm, l)
            worst = 0
            for i, (n, l, mm) in enumerate(m.STATES):
                for g, (ng, mu) in enumerate(groups):
                    if ng == n:
                        worst = max(worst, abs(sp.nsimplify(Mh[i, g]) - book.get((n, l, mm, mu), 0)))
            assert (worst == 0) is expect_equal, (n_max, worst)
