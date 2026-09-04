"""
TRUNK QA -- Paper 32 thm:forced_count: the D_F moduli chain on the 32-dim
finite Hilbert space, with a LINEAR *-representation of A_F.

STATUS (2026-09-03, trunk FULL run #3 remediation).  The previous version of
this file printed the chain 2048 -> 1024 -> 512 -> 272 -> 260 and a
matter-sector chain 512 -> 256 -> 128 -> 128.  Both endpoints were artefacts
of the algebra SAMPLE: `_a_f_basis()` took its 24 elements from the module's
`algebra_action`, whose quark block was kron(ew, m) -- bilinear, so 18 of the
24 elements were identically zero and the other 6 had no quark block -- and
the order-one condition therefore added almost nothing.  With a linear
*-representation of C (+) H (+) M_3(C) on C^32 (Connes-Chamseddine-Marcolli /
van Suijlekom Ch. 11: particles colour-blind, antileptons lambda, antiquarks
1_4 (x) m) the chain is

    2048 -> 1024 -> 512 -> 272 -> 32,

with matter-sector projection rank 16 and Majorana-block rank 16.  Three
independent routes (basis Gram; SVD-reduce then 10 random generic elements;
reduced Gram) agree; singular-value gap ~12 vs 1e-13.  The representation is
built HERE, independently of geovac/standard_model_triple.py, so the test
does not inherit that module's conventions; the module's own representation
is checked against this one in tests/test_standard_model_triple.py.

Regression guards: (i) the degenerate sample is re-installed in a scratch
routine and must reproduce 260 (proves the old number was the sample's);
(ii) random generic elements must give the same 32 as the basis.
"""
from __future__ import annotations

import numpy as np
import pytest

from geovac.almost_commutative import quaternion_to_matrix
from geovac.standard_model_triple import StandardModelFiniteTriple

F = StandardModelFiniteTriple(yukawa_nu=0.2, yukawa_e=0.3, yukawa_u=0.4, yukawa_d=0.5)
GAMMA32 = F.chirality_F()
J_U = F.real_structure_F()
N = 32
RDIM = 2 * N * N


# ---------------------------------------------------------------------------
# CCM representation of A_F = C (+) H (+) M_3(C) on H_F = C^16 (+) C^16
# ---------------------------------------------------------------------------

def pi_ccm(lam: complex, qc, m: np.ndarray) -> np.ndarray:
    q = quaternion_to_matrix(*qc)
    ew = np.zeros((4, 4), complex); ew[0:2, 0:2] = q; ew[2, 2] = lam; ew[3, 3] = np.conj(lam)
    P = np.zeros((N, N), complex)
    P[0:4, 0:4] = ew                                   # leptons
    P[4:16, 4:16] = np.kron(ew, np.eye(3))             # quarks (colour-blind)
    P[16:20, 16:20] = lam * np.eye(4)                  # antileptons
    P[20:32, 20:32] = np.kron(np.eye(4), np.asarray(m, complex))   # antiquarks: colour
    return P


def ccm_basis() -> list[np.ndarray]:
    out = []
    for lam in (1.0, 1j):
        out.append(pi_ccm(lam, (0, 0, 0, 0), np.zeros((3, 3), complex)))
    for k in range(4):
        qc = [0, 0, 0, 0]; qc[k] = 1
        out.append(pi_ccm(0, tuple(qc), np.zeros((3, 3), complex)))
    for i in range(3):
        for j in range(3):
            for val in (1.0, 1j):
                m = np.zeros((3, 3), complex); m[i, j] = val
                out.append(pi_ccm(0, (0, 0, 0, 0), m))
    return out


def _rv(M): return np.concatenate([M.real.ravel(), M.imag.ravel()])
def _fv(v): return v[:N * N].reshape(N, N) + 1j * v[N * N:].reshape(N, N)


def _reduced_space():
    """Orthonormal basis (as matrices) of the Hermitian, gamma-anticommuting,
    J-real subspace of M_32(C), computed by SVD of the stacked constraints."""
    E = np.eye(RDIM)
    rows = []
    for k in range(RDIM):
        M = _fv(E[:, k])
        rows.append(np.concatenate([_rv(M - M.conj().T),
                                    _rv(GAMMA32 @ M + M @ GAMMA32),
                                    _rv(J_U @ np.conj(M) @ J_U.conj().T - M)]))
    A = np.array(rows).T
    _, s, vt = np.linalg.svd(A, full_matrices=True)
    r = int((s > 1e-9).sum())
    Q = vt[r:].T
    return [_fv(Q[:, k]) for k in range(Q.shape[1])]


def _order_one_null(basisD, algebra):
    """Null space of the order-one map on span(basisD) for the given algebra
    elements; returns (dimension, solution matrices, singular gap)."""
    JbJ = [J_U @ np.conj(b) @ J_U.conj().T for b in algebra]
    G = np.zeros((len(basisD), len(basisD)))
    for a in algebra:
        for Jb in JbJ:
            cols = np.array([_rv((D @ a - a @ D) @ Jb - Jb @ (D @ a - a @ D)) for D in basisD]).T
            G += cols.T @ cols
    w, V = np.linalg.eigh(G)
    null = V[:, w < 1e-9]
    gap = (w[(w < 1e-9).sum()] if (w < 1e-9).sum() < len(w) else np.inf, w[(w < 1e-9).sum() - 1])
    sols = [sum(null[k, j] * basisD[k] for k in range(null.shape[0])) for j in range(null.shape[1])]
    return null.shape[1], sols, gap


@pytest.fixture(scope="module")
def reduced():
    return _reduced_space()


# ---------------------------------------------------------------------------
# The representation is a linear *-representation (the old one was not)
# ---------------------------------------------------------------------------

def test_ccm_representation_is_linear_and_star():
    rng = np.random.default_rng(1)
    for _ in range(5):
        l1, l2 = rng.normal(size=2) + 1j * rng.normal(size=2)
        q1 = tuple(rng.normal(size=4) + 1j * rng.normal(size=4)); q2 = tuple(rng.normal(size=4) + 1j * rng.normal(size=4))
        m1 = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); m2 = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
        add = pi_ccm(l1, q1, m1) + pi_ccm(l2, q2, m2) - pi_ccm(l1 + l2, tuple(a + b for a, b in zip(q1, q2)), m1 + m2)
        assert np.linalg.norm(add) < 1e-12
        # M_3(C) summand has a nonzero image on its own
        assert np.linalg.norm(pi_ccm(0, (0, 0, 0, 0), np.eye(3))) > 1
    B = ccm_basis()
    assert len(B) == 24 and all(np.linalg.norm(b) > 0 for b in B)
    # order-zero holds for the CCM representation
    assert max(np.linalg.norm(a @ (J_U @ np.conj(b) @ J_U.conj().T) - (J_U @ np.conj(b) @ J_U.conj().T) @ a)
               for a in B for b in B) < 1e-12


# ---------------------------------------------------------------------------
# The chain
# ---------------------------------------------------------------------------

def test_general_complex_dim_is_2048():
    assert RDIM == 2048


def test_chain_to_272(reduced):
    """Hermitian 1024 -> chirality 512 -> J-reality 272 (basis-independent)."""
    assert len(reduced) == 272
    E = np.eye(RDIM)
    herm = np.array([_rv(_fv(E[:, k]) - _fv(E[:, k]).conj().T) for k in range(RDIM)]).T
    chir = np.array([_rv(GAMMA32 @ _fv(E[:, k]) + _fv(E[:, k]) @ GAMMA32) for k in range(RDIM)]).T
    assert RDIM - np.linalg.matrix_rank(herm, tol=1e-9) == 1024
    assert RDIM - np.linalg.matrix_rank(np.concatenate([herm, chir]), tol=1e-9) == 512


def test_full_axiom_moduli_is_32(reduced):
    dim, sols, gap = _order_one_null(reduced, ccm_basis())
    assert dim == 32, dim
    # gap[1] < 1e-9 holds BY CONSTRUCTION (gap[1] = w[k-1] with
    # k = (w < 1e-9).sum(), so it is under tolerance whenever k > 0, which
    # `dim == 32` already guarantees).  /qa FULL #4 flagged it as tautological.
    # The real content is the SEPARATION: the null space is cleanly detached
    # from the rest of the spectrum, so the count is not tolerance-dependent.
    assert gap[0] > 1.0, gap[0]
    assert gap[0] / max(gap[1], 1e-15) > 1e6, (gap[0], gap[1])
    matter = np.array([_rv(D[0:16, 0:16]) for D in sols]).T
    majorana = np.array([_rv(D[0:16, 16:32]) for D in sols]).T
    assert np.linalg.matrix_rank(matter, tol=1e-9) == 16
    assert np.linalg.matrix_rank(majorana, tol=1e-9) == 16
    assert dim != 260 and dim != 128


def test_random_generic_elements_give_the_same_32(reduced):
    """Sample-independence: 10 random generic elements reproduce the basis count."""
    rng = np.random.default_rng(7)
    els = [pi_ccm(rng.normal() + 1j * rng.normal(),
                  tuple(rng.normal(size=4) + 1j * rng.normal(size=4)),
                  rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))) for _ in range(10)]
    dim, _, _ = _order_one_null(reduced, els)
    assert dim == 32


def test_degenerate_sample_reproduces_the_retired_260(reduced):
    """Regression guard: the pre-2026-09-03 sample (18/24 zero elements, six
    with no quark block) gives 260 -- the retired number was the sample's."""
    def old_action(lam, qc, m):
        q = quaternion_to_matrix(*qc)
        ew = np.zeros((4, 4), complex); ew[0:2, 0:2] = q; ew[2, 2] = lam; ew[3, 3] = np.conj(lam)
        M = np.zeros((N, N), complex); M[0:4, 0:4] = ew; M[4:16, 4:16] = np.kron(ew, np.asarray(m, complex))
        return M
    old = []
    for lam in (1.0, 1j):
        old.append(old_action(lam, (0, 0, 0, 0), np.zeros((3, 3), complex)))
    for k in range(4):
        qc = [0, 0, 0, 0]; qc[k] = 1; old.append(old_action(0, tuple(qc), np.zeros((3, 3), complex)))
    for i in range(3):
        for j in range(3):
            for val in (1.0, 1j):
                m = np.zeros((3, 3), complex); m[i, j] = val; old.append(old_action(0, (0, 0, 0, 0), m))
    assert sum(np.linalg.norm(b) > 0 for b in old) == 6
    dim, _, _ = _order_one_null(reduced, old)
    assert dim == 260


def _support(sols, rows, cols, tol=1e-9):
    """Positions in the given block hit by at least one basis solution."""
    hit = set()
    for D in sols:
        for i in rows:
            for j in cols:
                if abs(D[i, j]) > tol:
                    hit.add((i, j))
    return hit


def test_surviving_support_is_not_the_diagonal_yukawas(reduced):
    """The theorem's prose claim, pinned (added /qa trunk FULL run #4).

    The proof sketch said the order-one condition "kills every lepton-quark and
    off-flavour entry" and described the rank-16 matter projection as the four
    colour-diagonal Yukawas.  Measured, both are false: the matter block is
    full on the lepton and quark 2x2s, and the Majorana block is mostly
    lepton-quark.  This test is the falsifier for the retired wording.
    """
    dim, sols, _ = _order_one_null(reduced, ccm_basis())
    assert dim == 32, dim

    matter = _support(sols, range(0, 16), range(0, 16))
    majorana = _support(sols, range(0, 16), range(16, 32))

    # the matter block carries strictly more than four independent Yukawas:
    # 16 real moduli = 8 complex, so the four named ones are half of it
    m_cols = np.array([_rv(D[0:16, 0:16]) for D in sols]).T
    assert np.linalg.matrix_rank(m_cols, tol=1e-9) == 16

    # the load-bearing contrast: lepton-quark entries SURVIVE in the Majorana
    # block.  Leptons occupy the first 4 indices of each 16-block, quarks the
    # remaining 12; a mixed (lepton, antiquark) position is exactly the class
    # the retired sentence said was annihilated.
    # mixing runs BOTH ways: lepton row x antiquark column, and quark row x
    # antilepton column.  Leptons occupy indices 0-3 of each 16-block.
    def is_lepton(k):
        return (k % 16) < 4

    assert [(i, j) for (i, j) in majorana if is_lepton(i) != is_lepton(j)], \
        "no lepton-quark Majorana entry survived -- the retired sentence " \
        "would then have been right"

    # The paper's claim is about MODULI, not support positions (both happen to
    # come to 12 here, which is a coincidence of two different measures).
    # Assert the ranks: of the 16 real Majorana moduli, 12 are lepton-quark
    # and only 4 are same-species.
    def rank_on(positions):
        M = np.array([np.concatenate([[D[i, j].real for i, j in positions],
                                      [D[i, j].imag for i, j in positions]])
                      for D in sols]).T
        return np.linalg.matrix_rank(M, tol=1e-9)

    mixed_pos = [(i, j) for i in range(16) for j in range(16, 32)
                 if is_lepton(i) != is_lepton(j)]
    same_pos = [(i, j) for i in range(16) for j in range(16, 32)
                if is_lepton(i) == is_lepton(j)]
    assert rank_on(mixed_pos) == 12, rank_on(mixed_pos)
    assert rank_on(same_pos) == 4, rank_on(same_pos)
