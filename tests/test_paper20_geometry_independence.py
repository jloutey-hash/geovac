"""Three resource columns, three places on the forced/free seam.

WHAT THIS PINS
--------------
Papers 14/20 print three balanced resource columns side by side, and moving the
nuclei sorts them into three different categories:

    N_Pauli(R)   BIT-IDENTICAL across bond lengths          -> FORCED
    N_QWC(R)     MOVES, by up to ~3%                        -> heuristic
    lambda(R)    MOVES, by up to ~13%                       -> FREE

**N_Pauli is forced.** It is fixed by the angular selection rules (Gaunt /
Wigner 3j), which are functions of the quantum-number labels alone -- not of
where the nuclei sit.  This is the molecular analogue of Paper 22's
potential-independence theorem (ERI density depends on l_max, not on V(r)):
there the potential was varied and the sparsity held; here the geometry is
varied and the term count holds.

**N_QWC is not.** Qubit-wise-commuting groups are produced by a GREEDY
grouping pass, so the count depends on the algorithm's path, not only on the
operator's structure.  Measured: H2S 1264 -> 1282, PH3 1427 -> 1477,
SiH4 1566 -> 1551 between a molecule's own R and R = 3.015.  It belongs on the
engineering side of the seam, and a paper may not present it as structural.

**lambda is free.** It sums coefficient MAGNITUDES, which are the skeleton
projected through the geometry -- and the bond length is calibration data.

A CORRECTION THIS FILE RECORDS
------------------------------
The first version of this test asserted that N_QWC was geometry-independent
too.  It passed -- on NaH, LiH and HCl, which happen to be insensitive.
Widening to H2S / PH3 / SiH4 falsified it immediately.  Three molecules
agreeing is not an invariance; the QWC leg is now pinned in the opposite
direction so the false version cannot be re-derived.

PROVENANCE
----------
Surfaced 2026-08-31 while diagnosing the opposite defect: the balanced lambda
column of both papers had been computed at LiH's bond length for every
molecule, because `build_balanced_hamiltonian` defaulted `R = 3.015` and
ignored the spec.  The Pauli column came out right anyway -- because of the
invariance pinned here.  An invariance that silently rescues a wrong
measurement is worth writing down, because the next one may not be rescued.
See debug/qa/balanced_lambda_geometry_finding.md.
"""
from __future__ import annotations

import warnings

import pytest

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import nah_spec, hcl_spec, lih_spec, mgh2_spec


# Each build is reused across the tests that share a geometry pair; without the
# cache this file rebuilds the same Hamiltonians many times and lands over the
# per-file time budget (docs: debug/qa/test_suite_cost_memo.md).
_CACHE = {}


def _build(factory, R):
    key = (factory.__name__, R)
    if key not in _CACHE:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _CACHE[key] = build_balanced_hamiltonian(factory(max_n=2), R=R)
    return _CACHE[key]


def _lambda_ni(res) -> float:
    return sum(abs(complex(c)) for t, c in res["qubit_op"].terms.items() if t)


# (name, factory, R_own, R_other).  NaH is the cheapest build in the library
# (~0.3 s), so the fast cases run in the default suite rather than behind
# --slow: an invariance nobody runs is not a guard.
GEOMETRY_PAIRS = [
    ("NaH", nah_spec, 3.566, 3.015),
    ("LiH", lih_spec, 3.015, 3.566),
]


@pytest.mark.parametrize("name,factory,R_a,R_b", GEOMETRY_PAIRS,
                         ids=[g[0] for g in GEOMETRY_PAIRS])
def test_pauli_count_is_geometry_independent(name, factory, R_a, R_b):
    """FORCED: bit-identical, not merely close."""
    a = _build(factory, R_a)
    b = _build(factory, R_b)
    assert a["N_pauli"] == b["N_pauli"], (
        f"{name}: balanced Pauli count moved with geometry "
        f"({a['N_pauli']} at R={R_a} vs {b['N_pauli']} at R={R_b}). "
        f"The count is fixed by the angular selection rules alone; if it "
        f"moves, the selection-rule layer has become geometry-dependent."
    )


@pytest.mark.parametrize("name,factory,R_a,R_b", GEOMETRY_PAIRS,
                         ids=[g[0] for g in GEOMETRY_PAIRS])
def test_lambda_does_move(name, factory, R_a, R_b):
    """NON-TAUTOLOGY CONTROL for the invariance above.

    A builder that dropped R on the floor would satisfy the Pauli invariance
    trivially -- which is precisely the bug that motivated this file.  So the
    same pair of builds must move the free-side quantity.
    """
    a = _lambda_ni(_build(factory, R_a))
    b = _lambda_ni(_build(factory, R_b))
    assert abs(a - b) > 1e-6, (
        f"{name}: lambda_ni is identical at R={R_a} and R={R_b} "
        f"({a:.6f}); the builder is ignoring the geometry, so the invariance "
        f"test above proves nothing."
    )


def test_qwc_is_NOT_geometry_independent():
    """QWC grouping is a greedy heuristic, so its count is NOT structural.

    Pinned in this direction on purpose.  The first version of this file
    asserted the opposite and passed, because NaH / LiH / HCl happen to be
    insensitive; MgH2, H2S, PH3 and SiH4 all falsify it.  Anyone tempted to
    promote the QWC column to a structural claim has to break this test first.
    """
    # MgH2 at its own R vs a stretched 6.0 bohr: the cheapest witness in the
    # library (~9 s for the pair).  PH3/H2S/SiH4 also move at their own R vs
    # 3.015 but cost 30-60 s, which would put this file over budget and
    # tempt someone to mark it slow -- i.e. to stop running it.
    own = _build(mgh2_spec, 3.261)
    other = _build(mgh2_spec, 6.0)
    assert own["N_pauli"] == other["N_pauli"], (
        "control: MgH2's Pauli count should still be geometry-independent")
    assert own["n_qwc_balanced"] != other["n_qwc_balanced"], (
        "MgH2 QWC groups came out equal at two geometries; if the greedy "
        "grouper has become geometry-stable, this test's premise -- and the "
        "engineering-side classification of the QWC column -- needs "
        "revisiting.  Other known movers: PH3 1427/1477, H2S 1264/1282, "
        "SiH4 1566/1551 (own R vs 3.015)."
    )


def test_spec_geometry_is_the_default():
    """Omitting R must use the molecule's OWN bond length.

    Regression pin for the 2026-08-31 defect: R used to default to 3.015
    (LiH's) for every molecule, so `build_balanced_hamiltonian(nah_spec())`
    computed NaH at LiH's geometry and the published lambda column inherited
    it.
    """
    spec_default = _build(nah_spec, None)
    explicit_own = _build(nah_spec, 3.566)
    explicit_lih = _build(nah_spec, 3.015)

    assert _lambda_ni(spec_default) == pytest.approx(
        _lambda_ni(explicit_own), abs=1e-9), (
        "omitting R did not use the spec's own geometry")
    assert _lambda_ni(spec_default) != pytest.approx(
        _lambda_ni(explicit_lih), abs=1e-6), (
        "NaH at its own R is indistinguishable from NaH at LiH's R -- the "
        "control geometry is not actually different")


@pytest.mark.slow
def test_hcl_pauli_invariance_larger_witness():
    """Second, larger witness (~10 s): a 9,824-term operator."""
    a = _build(hcl_spec, 2.409)
    b = _build(hcl_spec, 3.015)
    assert a["N_pauli"] == b["N_pauli"] == 9824
    assert abs(_lambda_ni(a) - _lambda_ni(b)) > 1.0
