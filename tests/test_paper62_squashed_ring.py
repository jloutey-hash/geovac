"""Paper 62, Sec. hring, squashed-ring paragraph (record SQ,
docs/predictions/2026-11-02_squashed_ring.md, GO): regression guard.  Written 2026-11-02 as a
separate step after the paper (CLAUDE.md Sec. 9, guard-writing rule).

Rebuilds two small rings (N = 7 and 8, spacing 1.2 A, one Slater 1s per H fitted by 3
Gaussians, Lowdin site orbitals) on the ellipse z_j = R(e^{i phi_j} + eps e^{-i phi_j}) and
measures, in the ring-momentum basis, the share of h and of (pq|rs) that changes k by c.
Nothing is imported from debug/; the basis helpers are the already-guarded ones in
tests/test_water_ccsd.py.  The expected orders are written out by hand from the record's rule
(p(c) = min |q| with 2q = +-c mod N), not computed by the code under test.

  * the round ring (eps = 0) has every k-changing share at the rounding level;
  * at N = 8 the odd classes stay at the rounding level at every eps (the half-turn survives);
  * the local log-log slope between eps = 0.005 and 0.01 equals p(c) within 0.02 (the record's
    SQ2 bar) for every allowed class, h and g: N = 7 gives orders 3, 1, 2 for c = 1, 2, 3
    (c = 1 wraps around the ring); N = 8 gives 1, 2 for c = 2, 4.

Fire-tested 2026-11-02 with debug/qa/fire_test.py; every plant FIRED:
  1. the squash given charge 1 (an off-centre shift, z_j = R(e^{i phi_j} + eps));
  2. a charge-3 (triangular) deformation, z_j = R(e^{i phi_j} + eps e^{-2i phi_j}).
     (A first plant, stretching one axis only, stayed green: that is still a charge-2 squash
     plus a uniform scaling, so the order rule correctly holds for it.)
  3. the momentum change of (pq|rs) mislabelled (-kp + kq + kr - ks in place of
     -kp + kq - kr + ks);
  4. the hand-written wrap-around order at N = 7, c = 1 given as 1 (the naive reading).
"""
import importlib.util
import os

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('water_ccsd_for_p62_sq', os.path.join(HERE, 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(w)
GE = w.GE

ANG = 1.0 / 0.52917721067
SPACING = 1.2
EPS = (0.005, 0.01)
FLOOR = 1e-13
ZERO = 1e-13
SLOPE_BAR = 0.02
# p(c) from the record's rule, written out by hand (None = forbidden)
ORDERS = {7: {1: 3, 2: 1, 3: 2},
          8: {1: None, 2: 1, 3: None, 4: 2}}


def ring_integrals(N, eps):
    R = SPACING * ANG / (2 * np.sin(np.pi / N))
    nuc = []
    for j in range(N):
        z = R * (np.exp(2j * np.pi * j / N) + eps * np.exp(-2j * np.pi * j / N))
        nuc.append((np.array([z.real, z.imag, 0.0]), 1.0))
    fns, blocks = [], []
    for C, _ in nuc:
        w.sto_shell(C, 1, 0, 1.24, 3, fns, blocks)
    T = w.spherical_map(fns, blocks)
    S, h, g = w.ao_integrals(fns, nuc, T)
    X = GE.lowdin_orbitals(S)
    return GE.transform_integrals(X, h, g)


def class_shares(N, h, g):
    """Share of the Frobenius norm of h and of (pq|rs) whose momentum change is +-c mod N."""
    F = np.exp(2j * np.pi * np.outer(np.arange(N), np.arange(N)) / N) / np.sqrt(N)   # F[j, k]
    hk = F.conj().T @ h @ F
    gk = np.einsum('pqrs,pa,qb,rc,sd->abcd', g, F.conj(), F, F.conj(), F, optimize=True)
    k = np.arange(N)
    a, b = np.meshgrid(k, k, indexing='ij')
    dh = (b - a) % N
    A, B, C, D = np.meshgrid(k, k, k, k, indexing='ij')
    dg = (-A + B - C + D) % N
    out = {}
    for name, X, d in (('h', hk, dh), ('g', gk, dg)):
        p2 = np.abs(X) ** 2
        tot = p2.sum()
        cls = np.minimum(d, N - d)
        out[name] = {c: float(np.sqrt(p2[cls == c].sum() / tot)) for c in range(1, N // 2 + 1)}
    return out


@pytest.fixture(scope='module')
def shares():
    res = {}
    for N in ORDERS:
        for e in (0.0,) + EPS:
            res[(N, e)] = class_shares(N, *ring_integrals(N, e))
    return res


def test_round_ring_keeps_ring_momentum(shares):
    """Wrong answer rejected: a round ring (eps = 0) on which some class already leaks, i.e. a
    basis or a geometry that is not N-fold symmetric."""
    for N in ORDERS:
        for name in ('h', 'g'):
            for c, s in shares[(N, 0.0)][name].items():
                assert s < ZERO, (N, name, c, s)


def test_even_ring_odd_classes_stay_zero(shares):
    """Wrong answer rejected: a squash that also breaks the half-turn (a charge-1 squash, or a
    mislabelled class) so that odd momentum changes appear at N = 8."""
    for e in EPS:
        for name in ('h', 'g'):
            for c, p in ORDERS[8].items():
                if p is None:
                    assert shares[(8, e)][name][c] < ZERO, (e, name, c)


def test_slopes_match_the_charge_two_order_rule(shares):
    """Wrong answers rejected: a squash of charge 1 (orders c, not the charge-2 orders); the
    naive order 1 for c = 1 at N = 7 in place of the wrap-around order 3; a mislabelled
    (pq|rs) momentum change; a squash too large to be in the leading-power regime."""
    e1, e2 = EPS
    n_scored = 0
    for N, orders in ORDERS.items():
        for name in ('h', 'g'):
            for c, p in orders.items():
                if p is None:
                    continue
                s1, s2 = shares[(N, e1)][name][c], shares[(N, e2)][name][c]
                assert s1 >= FLOOR and s2 >= FLOOR, (N, name, c, s1, s2)
                slope = np.log(s2 / s1) / np.log(e2 / e1)
                assert abs(slope - p) < SLOPE_BAR, (N, name, c, slope, p)
                n_scored += 1
    assert n_scored == 10
