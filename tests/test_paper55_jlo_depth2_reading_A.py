"""Paper 55 thm:jlo_depth2_reading_A -- JLO depth-2 cocycle swap/S_3 symmetry.

Guard (group3 guard pass, 2026-11-03).  The theorem was backed only by
debug/sprint_jlo_depth2_compute.py.  This file recomputes the JLO depth-2
cochain coefficients on the production FockSpectralTriple(n_max=2) substrate
in exact sympy arithmetic, with its own implementation of

    c_m(a0,a1,a2) = (-1)^m/(m+2)! sum_{m0+m1+m2=m}
                    Tr_s(a0 D^{2m0} [D,a1] D^{2m1} [D,a2] D^{2m2}),

Tr_s = Tr(gamma .) (even flavour) or Tr (odd flavour).

What is pinned
--------------
* The 40 swap residuals chi(1,e_s,e_t) - chi(1,e_t,e_s) (10 unordered pairs
  x orders t^0, t^1 x 2 flavours) are exactly zero, and the panel is not
  vacuous (non-zero entries at both orders and both flavours).
  Wrong answer excluded: a non-zero residual.
* The S_3 orbit of chi(e0,e1,e2) in the EVEN flavour: one value, -13/98304
  at t^1 (t^0 is identically 0 on these three sectors).
* The commutator norms quoted in the proof: ||[gamma, D]||_F = sqrt(13)/4,
  ||[[D,e0],[D,e1]]||_F = 3/64.

Findings recorded as tests (prose-vs-computation, reported to the PM)
--------------------------------------------------------------------
* ODD flavour, S_3 orbit of chi(e0,e1,e2) at t^1 is NOT invariant: the
  transpositions that move a0 flip the sign (+-13/98304).  The theorem
  statement asserts full S_3 invariance "in both the even ... and odd
  flavours"; the driver only computed the even orbit.
* The 40 swap zeros are forced by reality: for ANY real-symmetric D and
  diagonal real gamma, transposition gives Tr(g P c_s Q c_t R) =
  Tr(g R c_t Q c_s P), so the a0 = 1 swap symmetry holds generically
  (control below: random real D -> zero residual; random complex Hermitian
  D -> O(1) residual).  The proof's "strictly stronger than cyclicity or
  commutativity" is true but the extra ingredient is reality of D.
"""
from __future__ import annotations

import itertools

import numpy as np
import pytest
import sympy as sp

from geovac.spectral_triple import FockSpectralTriple


def _substrate(n_max: int):
    st = FockSpectralTriple(n_max=n_max)
    D = st.dirac_operator
    g = st.grading
    N = st.dim_H
    E = []
    for s in range(st.n_sectors):
        M = sp.zeros(N, N)
        for i in range(N):
            if st._state_to_sector[i] == s:
                M[i, i] = 1
        E.append(M)
    return st, D, g, E


def _chi(D, g, a0, a1, a2, m: int, even: bool):
    N = D.shape[0]
    D2 = D * D
    P = [sp.eye(N)]
    for _ in range(m):
        P.append(P[-1] * D2)
    c1 = D * a1 - a1 * D
    c2 = D * a2 - a2 * D
    tot = sp.Integer(0)
    for m0 in range(m + 1):
        for m1 in range(m + 1 - m0):
            m2 = m - m0 - m1
            prod = a0 * P[m0] * c1 * P[m1] * c2 * P[m2]
            if even:
                prod = g * prod
            tot += prod.trace()
    return sp.nsimplify(sp.expand(tot)) * sp.Rational((-1) ** m) / sp.factorial(m + 2)


@pytest.fixture(scope="module")
def sub2():
    return _substrate(2)


def test_forty_swap_residuals_are_exact_zeros(sub2):
    st, D, g, E = sub2
    N = D.shape[0]
    I = sp.eye(N)
    assert st.n_sectors == 5 and N == 16
    zeros = 0
    nonzero_entries = {(m, ev): 0 for m in (0, 1) for ev in (True, False)}
    for s, t in itertools.combinations(range(5), 2):
        for m in (0, 1):
            for ev in (True, False):
                a = sp.simplify(_chi(D, g, I, E[s], E[t], m, ev))
                b = sp.simplify(_chi(D, g, I, E[t], E[s], m, ev))
                assert sp.simplify(a - b) == 0, (s, t, m, ev, a, b)
                zeros += 1
                if a != 0:
                    nonzero_entries[(m, ev)] += 1
    assert zeros == 40
    # non-vacuity: every (order, flavour) cell carries non-zero cocycle values
    assert all(v > 0 for v in nonzero_entries.values()), nonzero_entries


def test_even_flavour_s3_orbit_is_one_value(sub2):
    st, D, g, E = sub2
    for m, expected in ((0, sp.Integer(0)), (1, sp.Rational(-13, 98304))):
        vals = {sp.simplify(_chi(D, g, E[p[0]], E[p[1]], E[p[2]], m, True))
                for p in itertools.permutations((0, 1, 2))}
        assert vals == {expected}, (m, vals)


def test_odd_flavour_s3_orbit_flips_sign(sub2):
    """FINDING: the odd flavour is not S_3-invariant on chi(e0,e1,e2) at t^1."""
    st, D, g, E = sub2
    vals = [sp.simplify(_chi(D, g, E[p[0]], E[p[1]], E[p[2]], 1, False))
            for p in itertools.permutations((0, 1, 2))]
    assert set(vals) == {sp.Rational(13, 98304), sp.Rational(-13, 98304)}, vals


def test_commutator_norms_quoted_in_proof(sub2):
    st, D, g, E = sub2

    def frob(M):
        return sp.sqrt(sp.nsimplify(sum(x ** 2 for x in M)))

    assert sp.simplify(frob(g * D - D * g) - sp.sqrt(13) / 4) == 0
    c0 = D * E[0] - E[0] * D
    c1 = D * E[1] - E[1] * D
    assert sp.simplify(frob(c0 * c1 - c1 * c0) - sp.Rational(3, 64)) == 0


def test_nmax3_spot_check():
    st, D, g, E = _substrate(3)
    assert st.n_sectors == 9 and D.shape[0] == 40
    I = sp.eye(40)
    for s, t in ((0, 1), (2, 5)):
        for ev in (True, False):
            a = _chi(D, g, I, E[s], E[t], 1, ev)
            b = _chi(D, g, I, E[t], E[s], 1, ev)
            assert sp.simplify(a - b) == 0
    vals = {sp.simplify(_chi(D, g, E[p[0]], E[p[1]], E[p[2]], 1, True))
            for p in itertools.permutations((0, 1, 2))}
    assert len(vals) == 1


def _chi_np(D, g, a0, a1, a2, m, even):
    N = D.shape[0]
    D2 = D @ D
    P = [np.eye(N), D2]
    c1 = D @ a1 - a1 @ D
    c2 = D @ a2 - a2 @ D
    tot = 0
    for m0 in range(m + 1):
        for m1 in range(m + 1 - m0):
            pr = a0 @ P[m0] @ c1 @ P[m1] @ c2 @ P[m - m0 - m1]
            tot += np.trace(g @ pr if even else pr)
    return tot


def test_swap_symmetry_is_forced_by_reality_control():
    """Control documenting the mechanism: random real-symmetric D gives zero
    swap residual; a complex Hermitian D does not."""
    rng = np.random.default_rng(3)
    N = 16
    sect = rng.integers(0, 5, N)
    E = [np.diag((sect == s).astype(float)) for s in range(5)]
    g = np.diag(rng.choice([-1.0, 1.0], N))
    I = np.eye(N)
    res = {}
    for kind in ("real", "complex"):
        A = rng.normal(size=(N, N)) + (1j * rng.normal(size=(N, N)) if kind == "complex" else 0)
        D = (A + A.conj().T) / 2
        res[kind] = max(abs(_chi_np(D, g, I, E[s], E[t], m, ev) - _chi_np(D, g, I, E[t], E[s], m, ev))
                        for s, t in itertools.combinations(range(5), 2)
                        for m in (0, 1) for ev in (True, False))
    assert res["real"] < 1e-9
    assert res["complex"] > 1.0
