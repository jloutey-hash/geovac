"""The relativistic spinor ERI angular factor order, RESOLVED.

``geovac/composed_qubit_relativistic.py::_build_spinor_eri_block`` assembles

    <ab|cd> = sum_k  X_k(a,c) * X_k(d,b) * R^k(ac,bd)

with the second angular factor in Condon--Shortley / Dyall order.  It shipped
with ``X_k(b,d)`` -- the same factor-order defect corrected at seven scalar
sites on 2026-08-29 -- and the correction was held back until 2026-08-30
because no *internal* check distinguished the two orderings: particle
exchange and hermiticity hold for both, as the M_J-constrained phase algebra
predicts and a direct test confirmed.

This module contains the external discriminator that resolved it, and it is
the reason the shipped ordering is known to be wrong rather than merely
suspected.

THE CONSTRUCTION.  A jj-coupled shell of orbital angular momentum ``l``
spans the *same* space as the scalar spin-orbitals (m_l) x (m_s) of that
shell, the two bases being related by a Clebsch--Gordan unitary

    |kappa, m_j> = sum_ms <l, m_j-ms; 1/2, ms | j, m_j> |m_l> |ms>

The Coulomb operator is spin-independent, so in the scalar basis its angular
part is the Condon--Shortley product ``c^k(a,c) * c^k(d,b)``.  That ordering
is independently established: the four-dimensional quadrature-from-definition
arbitration recorded in ``debug/sprint_eri_evaluator_defects_memo.md``
decided it against two disagreeing evaluators, and ``casimir_ci`` (whose
``_gaunt_ck`` this module reuses directly, so there is no duplicate
implementation to drift) was corrected to match.  Rotating that tensor into
the jj basis yields what the relativistic module must reproduce.

Restricting to one shell makes every radial integral identical, so comparing
multipole-by-multipole is a purely angular test.  The decisive statistic is
the sign pattern, which no positive rescaling can alter.

WHAT THIS ESTABLISHES, AND WHAT IT DOES NOT.  The test is a *cross-basis
consistency* check: it shows the jj-coupled module agrees with, or departs
from, the scalar Condon--Shortley convention.  The physical anchoring of
that convention is external to this file -- the quadrature-from-definition
arbitration in ``debug/sprint_eri_evaluator_defects_memo.md``, and the
textbook Slater--Condon rule.  Within that scope the verdict is sharp,
because sign patterns discriminate ordering and no positive rescale can
alter a sign.  Normalization is a separate question, guarded by the
``scale == 1`` assertions below (which must stay tight: ``_compare`` fits
``scale`` and measures the residual only after rescaling, so a loose band
there would admit a uniformly mis-normalized tensor).

THE RESULT.  For the corrected ordering the rotated reference is reproduced
to machine precision at every discriminating multipole, while the shipped
ordering agrees on little more than half the signs and admits no repairing
rescale:

    shell   k     corrected            shipped
    2p      2     168/168, 1.7e-16     96/168  (57%), resid 0.90
    3d      2     1004/1004, 2.5e-16   612/1004 (61%), resid 1.2
    3d      4     788/788, 4.4e-16     428/788  (54%), resid 0.89

At k=0 both agree exactly -- the built-in control, since a scalar monopole
forces the two Gaunt factors' own selection rules to imply the shared-q
condition, so no discrimination is possible there.

The d-shell rows were added 2026-08-30 after a review asked whether the
verdict was an artifact of the 2p case.  It is not.

Record: ``debug/qa/delta5_run_notes.md``, ``debug/qa/delta6_run_notes.md``.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from sympy import S
from sympy.physics.wigner import clebsch_gordan

from geovac.casimir_ci import _gaunt_ck
import geovac.composed_qubit_relativistic as R

# (label, l, [(kappa, 2*m_j), ...], discriminating multipoles)
SHELLS = {
    "2p": (1, [(1, -1), (1, 1),
               (-2, -3), (-2, -1), (-2, 1), (-2, 3)], [2]),
    "3d": (2, [(2, -3), (2, -1), (2, 1), (2, 3),
               (-3, -5), (-3, -3), (-3, -1), (-3, 1), (-3, 3), (-3, 5)],
           [2, 4]),
}


# The decisive comparison support at each discriminating multipole.
# Asserted so that a collapsed mask cannot masquerade as agreement.
EXPECTED_SUPPORT = {("2p", 2): 168, ("3d", 2): 1004, ("3d", 4): 788}


def _scalar(l: int):
    return [(ml, ms2) for ml in range(-l, l + 1) for ms2 in (-1, 1)]


def _cg(l: int, spinors) -> np.ndarray:
    n = len(spinors)
    sc = _scalar(l)
    assert len(sc) == n
    U = np.zeros((n, n))
    for A, (kap, tmj) in enumerate(spinors):
        j = S(abs(kap) * 2 - 1) / 2
        mj = S(tmj) / 2
        for a, (ml, ms2) in enumerate(sc):
            ms = S(ms2) / 2
            if ml + ms == mj:
                U[A, a] = float(clebsch_gordan(S(l), S(1) / 2, j,
                                               S(ml), ms, mj))
    return U


def _reference(l: int, k: int, U: np.ndarray) -> np.ndarray:
    """The verified scalar angular tensor, rotated into the jj basis.

    The M_L guard is explicit here.  Without it the tensor carries entries
    that violate total-M_J conservation; they happen to fall outside the
    comparison mask, so the verdict is unchanged either way -- but a helper
    that is correct only by luck of masking is exactly the pattern this arc
    has been correcting, so the constraint is imposed rather than relied on.
    """
    sc = _scalar(l)
    n = len(sc)
    A = np.zeros((n, n, n, n))
    for a, (mla, msa) in enumerate(sc):
        for b, (mlb, msb) in enumerate(sc):
            for c, (mlc, msc) in enumerate(sc):
                for d, (mld, msd) in enumerate(sc):
                    if msa != msc or msb != msd:
                        continue
                    if mla + mlb != mlc + mld:      # shared-q / M_L rule
                        continue
                    A[a, b, c, d] = (_gaunt_ck(l, mla, l, mlc, k)
                                     * _gaunt_ck(l, mld, l, mlb, k))
    return np.einsum('Aa,Bb,Cc,Dd,abcd->ABCD', U, U, U, U, A, optimize=True)


def _module(spinors, k: int, corrected: bool) -> np.ndarray:
    n = len(spinors)
    X = np.zeros((n, n))
    for A, (ka, ta) in enumerate(spinors):
        for C, (kc, tc) in enumerate(spinors):
            X[A, C] = R.jj_angular_Xk(ka, ta, kc, tc, k)
    T = np.zeros((n, n, n, n))
    for A in range(n):
        for B in range(n):
            for C in range(n):
                for D in range(n):
                    if spinors[A][1] + spinors[B][1] == \
                       spinors[C][1] + spinors[D][1]:
                        T[A, B, C, D] = X[A, C] * (X[D, B] if corrected
                                                   else X[B, D])
    return T


def _compare(ref: np.ndarray, T: np.ndarray):
    mask = (np.abs(ref) > 1e-10) & (np.abs(T) > 1e-10)
    agree = int(np.sum(np.sign(ref[mask]) == np.sign(T[mask])))
    scale = float(np.sum(ref[mask] * T[mask]) / np.sum(T[mask] ** 2))
    resid = float(np.max(np.abs(ref[mask] - scale * T[mask]))
                  / np.max(np.abs(ref[mask])))
    return agree, int(mask.sum()), scale, resid


@pytest.mark.parametrize("shell", list(SHELLS))
def test_cg_transform_is_unitary(shell):
    """The rotation must be exact, or the reference means nothing."""
    l, spinors, _ = SHELLS[shell]
    U = _cg(l, spinors)
    dev = np.max(np.abs(U @ U.T - np.eye(len(spinors))))
    assert dev < 1e-12, f"{shell}: CG matrix not unitary, max dev {dev:.2e}"


@pytest.mark.parametrize("shell", list(SHELLS))
def test_reference_conserves_m_j(shell):
    """The rotated reference must respect total M_J.

    Guards the helper directly rather than leaving it to the comparison
    mask.
    """
    l, spinors, ks = SHELLS[shell]
    U = _cg(l, spinors)
    for k in ks:
        ref = _reference(l, k, U)
        bad = [(A, B, C, D)
               for A in range(len(spinors)) for B in range(len(spinors))
               for C in range(len(spinors)) for D in range(len(spinors))
               if abs(ref[A, B, C, D]) > 1e-10
               and spinors[A][1] + spinors[B][1] !=
               spinors[C][1] + spinors[D][1]]
        assert not bad, f"{shell} k={k}: {len(bad)} M_J-violating entries"


@pytest.mark.parametrize("shell", list(SHELLS))
def test_monopole_cannot_discriminate(shell):
    """k=0 control: a scalar monopole makes the two orderings coincide.

    Both must match the reference exactly.  Note the honest scope: at k=0
    the X matrix is exactly the identity, so the two code branches operate
    on bit-identical arrays and this exercises delta-structure plus one
    overall constant -- not the 3j or reduced-matrix-element machinery the
    k>0 verdict rests on.  It is a null control (the discriminator must be
    silent where no discrimination is possible), not a proof that the
    comparison machinery is sound.
    """
    l, spinors, _ = SHELLS[shell]
    U = _cg(l, spinors)
    ref = _reference(l, 0, U)
    for corrected in (False, True):
        agree, total, scale, resid = _compare(ref, _module(spinors, 0,
                                                           corrected))
        assert agree == total, f"{shell} k=0 sign mismatch ({agree}/{total})"
        assert abs(scale - 1.0) < 1e-10
        assert resid < 1e-12


@pytest.mark.parametrize("shell", list(SHELLS))
def test_corrected_ordering_matches_scalar_reference(shell):
    """The Condon--Shortley order reproduces the verified scalar result.

    This is the external check that resolved the ordering.  Parametrised
    over two shells: the preference is a property of the rule, not of 2p.
    """
    l, spinors, ks = SHELLS[shell]
    U = _cg(l, spinors)
    for k in ks:
        ref = _reference(l, k, U)
        agree, total, scale, resid = _compare(ref, _module(spinors, k, True))
        # Pin the comparison support.  Without this, a change that collapsed
        # the mask would make `agree == total` trivially true at 0/0 and the
        # discriminator would report agreement while comparing nothing.
        assert total == EXPECTED_SUPPORT[(shell, k)], (
            f"{shell} k={k}: comparison support {total}, expected "
            f"{EXPECTED_SUPPORT[(shell, k)]} -- the mask has changed"
        )
        assert agree == total, f"{shell} k={k}: signs {agree}/{total}"
        # The scale guard is load-bearing and must stay TIGHT.  _compare
        # least-squares-fits `scale` and measures the residual only AFTER
        # rescaling, so without this a uniform rescale of X_k (say X/2)
        # would pass with signs intact and residual 0.  The k=0 control
        # cannot cover it either: X_0 is exactly the identity, so k=0 pins
        # no 3j or reduced-matrix-element normalization.
        assert abs(scale - 1.0) < 1e-10, (
            f"{shell} k={k}: scale {scale} != 1 -- the module's angular "
            f"normalization has drifted from the scalar convention"
        )
        assert resid < 1e-12, f"{shell} k={k}: residual {resid:.2e}"


@pytest.mark.parametrize("shell", list(SHELLS))
def test_shipped_ordering_is_excluded(shell):
    """The retired order disagrees, and no rescaling repairs it.

    Recorded rather than deleted: the retired ordering shipped for a long
    time, and its exclusion is the load-bearing evidence for the change.
    """
    l, spinors, ks = SHELLS[shell]
    U = _cg(l, spinors)
    for k in ks:
        ref = _reference(l, k, U)
        agree, total, scale, resid = _compare(ref, _module(spinors, k, False))
        assert agree < total, (
            f"{shell} k={k}: the retired ordering should NOT match "
            f"({agree}/{total})"
        )
        assert 0.4 < agree / total < 0.7, (
            f"{shell} k={k}: expected ~50-60% sign agreement, got "
            f"{100 * agree / total:.1f}%"
        )
        assert resid > 0.1, (
            f"{shell} k={k}: a rescaling repaired the retired ordering "
            f"(residual {resid:.2e}) -- the discriminator has lost its power"
        )


def test_production_module_uses_the_corrected_order():
    """The fix is in the shipped source, not only in this test."""
    # Anchored to the module's own __file__: a cwd-relative path would let
    # this guard be defeated simply by running pytest from another directory.
    src = open(R.__file__, encoding="utf-8").read()
    assert src.count("for (d, b), bds in ac_Xk.items():") == 1, (
        "the corrected Condon-Shortley loop key is missing from "
        "_build_spinor_eri_block"
    )
    assert "for (b, d), bds in ac_Xk.items():" not in src, (
        "the retired factor order has come back"
    )


@pytest.mark.parametrize("mol,nmax,q,n_pauli,lam", [
    ("LiH", 1,  6,    9,  10.15),
    ("BeH", 1,  6,    9,  63.44),
    ("CaH", 1,  4,    6,   2.03),
    ("LiH", 2, 30, 1501,  39.53),
    ("BeH", 2, 30, 1501, 142.52),
    ("CaH", 2, 20,  998,  17.91),
])
def test_rel_resource_numbers(mol, nmax, q, n_pauli, lam):
    """Pin Paper 14 tab:spinor_resource under the corrected ordering.

    All six n_max=1,2 rows, not just LiH: the retired factor order gave
    1413 (LiH/BeH) and 942 (CaH) at n_max=2, so each row independently
    locks out a revert.  The n_max=3 rows cost ~8 min apiece to build and
    are left to the paper's own record.
    """
    from geovac import molecular_spec as MS

    spec_fn = getattr(MS, f"{mol.lower()}_spec_relativistic")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = R.build_composed_hamiltonian_relativistic(spec_fn(max_n=nmax),
                                                      verbose=False)
    assert r["Q"] == q, f"{mol} n_max={nmax}: Q={r['Q']}, expected {q}"
    assert r["N_pauli"] == n_pauli, (
        f"{mol} n_max={nmax}: N_pauli={r['N_pauli']}, expected {n_pauli} "
        f"(the retired factor order gave 1413/942 at n_max=2)"
    )
    assert r["lambda_ni"] == pytest.approx(lam, abs=5e-2)
