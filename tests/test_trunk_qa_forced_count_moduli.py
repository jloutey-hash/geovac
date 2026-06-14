"""Trunk QA — Claim B: the Forced-Count moduli chain 1024 -> 512 -> 128 -> 128.

CONTEXT
=======
Paper 32 Theorem ``thm:forced_count`` states the reduction chain (eq.
``forced_count_chain``)

    dim_R M(D_F) :  1024 --Herm.--> 512 --{gamma_F,D_F}=0--> 128
                                        --[J_F,D_F]=0 + order-one--> 128

and the proof sketch cites it as "verified bit-exactly at n_max in {2,3} in
module geovac/standard_model_triple.py (45 tests passing)".

But those 45 tests verify the AXIOM CONSISTENCY of ONE constructed D_F (J^2,
{gamma,D}=0, order-zero/order-one residuals, the Higgs/gauge census, ...).
They do NOT parameter-count the moduli space, and so do NOT verify the chain
1024 -> 512 -> 128 -> 128.

This file ATTEMPTS the genuine, falsifiable thing the citation implies: it
COMPUTES dim_R M(D_F) at each reduction step by counting the real dimension of
the solution space of the corresponding LINEAR constraints, on the SAME
H_F = C^32 / gamma_F / J_F that StandardModelFiniteTriple uses.  Each constraint
(Hermiticity, chirality-anticommutation {gamma_F, D_F}=0, J-reality
[J_F, D_F]=0, Connes order-one) is a real-linear operator on the 2048-real-dim
space of 32x32 complex matrices; moduli dim = 2048 - rank(stacked operator).

RESULT (this is the load-bearing finding)
==========================================
The genuine constraint-counting on the module's own H_F/gamma_F/J_F gives

    general 32x32          : 2048
    + Hermitian            : 1024
    + {gamma_F,D_F}=0      :  512    (paper chain claims 128 here)
    + [J_F,D_F]=0          :  272    (paper chain claims 128 here)
    + order-one            :  260    (paper chain claims 128 here)

The paper's chain 1024 -> 512 -> 128 -> 128 DOES NOT reproduce on the 32-dim
space: the {gamma,D}=0 step gives 512 (not 128), and the full-axiom moduli is
260 (not 128).

The final number 128 IS reproducible, but only under a DIFFERENT convention than
the chain describes: it is the real dimension of Hermitian, chirality-odd Dirac
operators on the 16-dim MATTER sector (before the J-doubling to C^32):

    matter 16x16 general   : 512
    + Hermitian            : 256
    + chirality-odd        : 128

So "128 per generation" is a defensible matter-sector count, but the chain
1024 -> 512 -> 128 -> 128 with the stated arrow labels (J-reality / order-one
doing the 128 -> 128 step on the 32-dim space) is NOT what the constraints
produce.  The "verified bit-exactly ... (45 tests)" citation is therefore
mis-pointed: no test parameter-counts the chain, and the chain as written does
not reproduce.

These tests PASS by asserting the TRUE genuine counts (they would fail if the
constraint counting changed), and they explicitly DOCUMENT the divergence from
the paper's chain.  They are diagnostic, not a rubber stamp of the claim.
"""
from __future__ import annotations

import numpy as np
import pytest

from geovac.standard_model_triple import StandardModelFiniteTriple


F = StandardModelFiniteTriple()
GAMMA32 = np.real(F.chirality_F())          # 32x32 real diag +-1
J_U = F.real_structure_F()                  # unitary part U of J = U K


# ---------------------------------------------------------------------------
# Generic real-linear constraint counting machinery.
# ---------------------------------------------------------------------------


def _make_counter(N, gamma_block=None):
    rdim = 2 * N * N
    eye = np.eye(rdim)

    def fromv(v):
        return v[: N * N].reshape(N, N) + 1j * v[N * N :].reshape(N, N)

    def rv(M):
        return np.concatenate([np.real(M).ravel(), np.imag(M).ravel()])

    inputs = [fromv(eye[:, k]) for k in range(rdim)]

    def cmat(func):
        return np.array([rv(func(M)) for M in inputs]).T

    def nulldim(*mats):
        G = np.zeros((rdim, rdim))
        for A in mats:
            G += A.T @ A
        return rdim - np.linalg.matrix_rank(G, tol=1e-7)

    return rdim, cmat, nulldim, inputs, rv


# ---------------------------------------------------------------------------
# 32-dim full-Hilbert-space constraints (the space the paper's chain lives on).
# ---------------------------------------------------------------------------

_RDIM32, _CMAT32, _NULL32, _INPUTS32, _RV32 = _make_counter(32)


def _herm(M):
    return M - M.conj().T


def _chir32(M):
    return GAMMA32 @ M + M @ GAMMA32


def _jreal32(M):
    return J_U @ np.conj(M) @ J_U.conj().T - M


# A_F real basis (for order-one)
def _a_f_basis():
    out = []
    for lam in (1.0, 1j):
        out.append(F.algebra_action(lam, (0, 0, 0, 0), np.zeros((3, 3), complex)))
    for k in range(4):
        q = [0, 0, 0, 0]
        q[k] = 1
        out.append(F.algebra_action(0, tuple(q), np.zeros((3, 3), complex)))
    for i in range(3):
        for j in range(3):
            for val in (1.0, 1j):
                m = np.zeros((3, 3), complex)
                m[i, j] = val
                out.append(F.algebra_action(0, (0, 0, 0, 0), m))
    return out


_A_BASIS = _a_f_basis()


def _order_one_matrix():
    JbJ = [J_U @ np.conj(b) @ J_U.conj().T for b in _A_BASIS]
    cols = []
    for M in _INPUTS32:
        parts = []
        for a in _A_BASIS:
            Da = M @ a - a @ M
            for Jb in JbJ:
                parts.append(_RV32(Da @ Jb - Jb @ Da))
        cols.append(np.concatenate(parts))
    return np.array(cols).T


# Cache the heavy 32-dim constraint matrices once.
_H32 = _CMAT32(_herm)
_C32 = _CMAT32(_chir32)
_J32 = _CMAT32(_jreal32)


# ---------------------------------------------------------------------------
# The genuine counts on the 32-dim space (fast steps).
# ---------------------------------------------------------------------------


def test_general_complex_dim_is_2048():
    assert _RDIM32 == 2048


def test_hermitian_count_is_1024():
    assert _NULL32(_H32) == 1024


def test_herm_plus_chirality_is_512_not_128():
    """{gamma_F,D_F}=0 on Hermitian 32x32 gives 512, NOT the paper's claimed 128."""
    got = _NULL32(_H32, _C32)
    assert got == 512, f"got {got}"
    assert got != 128, "paper chain claims this step reaches 128; genuine count is 512"


def test_herm_chir_jreal_is_272_not_128():
    """Adding J-reality gives 272 on the 32-dim space, NOT 128."""
    got = _NULL32(_H32, _C32, _J32)
    assert got == 272, f"got {got}"
    assert got != 128


@pytest.mark.slow
def test_full_axiom_moduli_is_260_not_128():
    """Hermitian + chirality + J-reality + order-one gives 260, NOT 128.

    This is the full real-spectral-triple axiom set the theorem lists.  The
    genuine moduli dimension is 260; the paper's chain terminates at 128.
    (Marked slow: builds a 1.18M-row order-one constraint operator.)
    """
    O32 = _order_one_matrix()
    got = _NULL32(_H32, _C32, _J32, O32)
    assert got == 260, f"got {got}"
    assert got != 128


def test_paper_chain_intermediates_do_not_reproduce():
    """The stated chain 1024 -> 512 -> 128 -> 128 fails at the 128 steps.

    Genuine 32-dim chain: 2048 -> 1024 (Herm) -> 512 ({g,D}=0) -> 272 (J-real).
    The paper writes 1024 -> 512 -> 128 -> 128.  The two 128's do not appear.
    """
    chain = [
        _RDIM32,                       # 2048
        _NULL32(_H32),                 # 1024
        _NULL32(_H32, _C32),           # 512
        _NULL32(_H32, _C32, _J32),     # 272
    ]
    assert chain == [2048, 1024, 512, 272]
    paper_intermediates = {128}
    assert 128 not in set(chain[1:]), (
        "no genuine 32-dim intermediate equals the paper's 128"
    )


# ---------------------------------------------------------------------------
# Where DOES 128 come from?  The matter-sector convention.
# ---------------------------------------------------------------------------

_GM = GAMMA32[:16, :16]                     # matter-sector chirality
_RDIM16, _CMAT16, _NULL16, _INPUTS16, _RV16 = _make_counter(16)


def _chir16(M):
    return _GM @ M + M @ _GM


_H16 = _CMAT16(_herm)
_C16 = _CMAT16(_chir16)


def test_128_is_matter_sector_hermitian_chirality_odd():
    """128 = dim_R of Hermitian, chirality-odd D on the 16-dim MATTER sector.

    matter 16x16 general (512) -> Hermitian (256) -> chirality-odd (128).
    This is a legitimate count, but it is NOT the paper's stated 32-dim chain
    with J-reality / order-one arrows.  It documents that the FINAL number 128
    is defensible while the CHAIN that allegedly produces it is not.
    """
    assert _RDIM16 == 512
    assert _NULL16(_H16) == 256
    assert _NULL16(_H16, _C16) == 128


def test_matter_chain_is_512_256_128():
    chain = [_RDIM16, _NULL16(_H16), _NULL16(_H16, _C16)]
    assert chain == [512, 256, 128]
