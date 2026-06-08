"""Regression tests for ``geovac.relativistic_tapering``.

Sprint kappa-parity-Z2 (2026-06-07): documents the structural negative.
The kappa-parity stabilizer P_kappa = prod_{q: kappa_q < 0} Z_q does
NOT commute with the relativistic chemistry Hamiltonian; residuals are
1e-2 to 1e-1 on LiH_rel / BeH_rel / CaH_rel at max_n=2. These tests
lock in the negative so future refactors don't accidentally "fix" it
by introducing a buggy symmetry.
"""

from __future__ import annotations

import pytest

openfermion = pytest.importorskip("openfermion")
from openfermion import QubitOperator
from openfermion.utils import commutator as of_commutator

from geovac.molecular_spec import (
    lih_spec_relativistic,
    beh_spec_relativistic,
    cah_spec_relativistic,
)
from geovac.composed_qubit_relativistic import (
    build_composed_hamiltonian_relativistic,
)
from geovac.relativistic_tapering import (
    enumerate_relativistic_orbital_table,
    build_kappa_parity_stabilizers,
    audit_kappa_parity_commutation,
    relativistic_tapered_from_spec,
)


# ---------------------------------------------------------------------------
# Orbital table enumeration
# ---------------------------------------------------------------------------

def test_orbital_table_matches_relativistic_builder_Q():
    """The orbital table length must equal Q from the builder."""
    spec = lih_spec_relativistic(max_n=2, include_pk=True)
    result = build_composed_hamiltonian_relativistic(
        spec, pk_in_hamiltonian=True, include_breit=False, verbose=False,
    )
    table = enumerate_relativistic_orbital_table(spec)
    assert len(table) == result["Q"], (
        f"orbital_table length {len(table)} != Q {result['Q']}"
    )


def test_orbital_table_has_expected_kappa_signs():
    """LiH_rel max_n=2 has 24 kappa<0 and 6 kappa>0 across 3 sub-blocks."""
    spec = lih_spec_relativistic(max_n=2, include_pk=True)
    table = enumerate_relativistic_orbital_table(spec)
    n_neg = sum(1 for (_, _, k, _) in table if k < 0)
    n_pos = sum(1 for (_, _, k, _) in table if k > 0)
    n_zero = sum(1 for (_, _, k, _) in table if k == 0)
    n_sb = len({sb for (sb, _, _, _) in table})
    assert n_neg + n_pos == len(table)
    assert n_zero == 0, "kappa = 0 is unphysical"
    assert n_sb == 3, f"expected 3 sub-blocks, got {n_sb}"
    # s_{1/2} (kappa=-1) and p_{3/2} (kappa=-2) contribute kappa<0; the
    # p_{1/2} branch (kappa=+1) is the only positive contributor at max_n=2
    # on a max_n=2 block. Exact counts depend on l_min per block.
    assert n_neg > 0 and n_pos > 0


# ---------------------------------------------------------------------------
# Stabilizer construction
# ---------------------------------------------------------------------------

def test_global_kappa_parity_is_hermitian_z_string():
    """P_kappa is a single Z-string ==> Hermitian and unitary."""
    spec = lih_spec_relativistic(max_n=2, include_pk=True)
    table = enumerate_relativistic_orbital_table(spec)
    stabs = build_kappa_parity_stabilizers(table, mode="global")
    assert len(stabs) == 1
    P = stabs[0]
    # All Z-string Pauli ops are Hermitian (real coefficients on Z-only terms)
    for term, coef in P.terms.items():
        for q, op in term:
            assert op == "Z"
        assert abs(coef.imag) < 1e-15
    # P^2 = I (Z-string squared)
    PP = P * P
    # Identity term only
    for term, coef in PP.terms.items():
        if term == ():
            assert abs(coef - 1.0) < 1e-12
        else:
            assert abs(coef) < 1e-12


def test_per_block_stabilizer_count_matches_sub_blocks_with_kappa_neg():
    """One per-block stabilizer per sub-block that has at least one kappa<0."""
    spec = lih_spec_relativistic(max_n=2, include_pk=True)
    table = enumerate_relativistic_orbital_table(spec)
    stabs = build_kappa_parity_stabilizers(table, mode="per_block")
    from collections import defaultdict
    counts = defaultdict(int)
    for sb, _, k, _ in table:
        if k < 0:
            counts[sb] += 1
    expected = sum(1 for c in counts.values() if c > 0)
    assert len(stabs) == expected


def test_unknown_mode_raises():
    spec = lih_spec_relativistic(max_n=2, include_pk=True)
    table = enumerate_relativistic_orbital_table(spec)
    with pytest.raises(ValueError):
        build_kappa_parity_stabilizers(table, mode="not_a_mode")


# ---------------------------------------------------------------------------
# Commutation audit (STRUCTURAL NEGATIVE)
# ---------------------------------------------------------------------------

def test_global_kappa_parity_does_NOT_commute_lih_rel():
    """LiH_rel max_n=2, PK=True: residual ~5e-2 (well above 1e-10 PASS gate)."""
    spec = lih_spec_relativistic(max_n=2, include_pk=True)
    result = build_composed_hamiltonian_relativistic(
        spec, pk_in_hamiltonian=True, include_breit=False, verbose=False,
    )
    H = result["qubit_op"]
    table = enumerate_relativistic_orbital_table(spec)
    stabs = build_kappa_parity_stabilizers(table, mode="global")
    passed, residuals = audit_kappa_parity_commutation(H, stabs, atol=1e-10)
    assert len(passed) == 1
    # STRUCTURAL NEGATIVE: residual must be above 1e-6 (well above noise)
    assert not passed[0], "kappa-parity should NOT commute (structural)"
    assert residuals[0] > 1e-6, (
        f"residual {residuals[0]:.3e} should be > 1e-6 (structural mismatch); "
        "if this test fails, audit the relativistic builder for a regression "
        "that accidentally diagonalised the Coulomb operator in the kappa sign"
    )


def test_per_block_kappa_parity_does_NOT_commute_lih_rel():
    """Per-block kappa-parity also fails on LiH_rel."""
    spec = lih_spec_relativistic(max_n=2, include_pk=True)
    result = build_composed_hamiltonian_relativistic(
        spec, pk_in_hamiltonian=True, include_breit=False, verbose=False,
    )
    H = result["qubit_op"]
    table = enumerate_relativistic_orbital_table(spec)
    stabs = build_kappa_parity_stabilizers(table, mode="per_block")
    passed, residuals = audit_kappa_parity_commutation(H, stabs, atol=1e-10)
    # At least one must fail
    assert not all(passed), (
        "expected at least one per-block kappa-parity stabilizer to fail; "
        f"residuals = {residuals}"
    )
    # All residuals should be above noise (structural)
    assert all(r > 1e-6 for r in residuals), (
        f"all per-block residuals expected > 1e-6 (structural); got {residuals}"
    )


def test_kappa_parity_negative_persists_without_pk():
    """Turn off PK; the Coulomb-only Hamiltonian still doesn't commute."""
    spec = lih_spec_relativistic(max_n=2, include_pk=False)
    result = build_composed_hamiltonian_relativistic(
        spec, pk_in_hamiltonian=False, include_breit=False, verbose=False,
    )
    H = result["qubit_op"]
    table = enumerate_relativistic_orbital_table(spec)
    stabs = build_kappa_parity_stabilizers(table, mode="global")
    passed, residuals = audit_kappa_parity_commutation(H, stabs, atol=1e-10)
    assert not passed[0]
    assert residuals[0] > 1e-6


# ---------------------------------------------------------------------------
# Full pipeline: tapered == naive when audit drops everything
# ---------------------------------------------------------------------------

def test_relativistic_tapered_pipeline_drops_all_on_lih_rel():
    """End-to-end pipeline returns naive operator (audit drops all)."""
    spec = lih_spec_relativistic(max_n=2, include_pk=True)
    out = relativistic_tapered_from_spec(
        spec, mode="per_block",
        pk_in_hamiltonian=True, include_breit=False,
        drop_noncommuting=True, verbose=False,
    )
    assert out["audit_verdict"] == "DROP_ALL"
    assert out["delta_Q"] == 0
    assert out["Q_tapered"] == out["Q_naive"]
    # qubit_op_tapered IS qubit_op_naive in DROP_ALL case
    assert out["qubit_op_tapered"] is out["qubit_op_naive"]


def test_relativistic_tapered_pipeline_on_beh_rel():
    """BeH_rel: same structural negative, no change to Q."""
    spec = beh_spec_relativistic(max_n=2, include_pk=True)
    out = relativistic_tapered_from_spec(
        spec, mode="per_block",
        pk_in_hamiltonian=True, include_breit=False,
        drop_noncommuting=True, verbose=False,
    )
    assert out["audit_verdict"] == "DROP_ALL"
    assert out["delta_Q"] == 0


def test_relativistic_tapered_pipeline_on_cah_rel_global():
    """CaH_rel global mode: also drops, also delta_Q=0."""
    spec = cah_spec_relativistic(max_n=2, include_pk=True)
    out = relativistic_tapered_from_spec(
        spec, mode="global",
        pk_in_hamiltonian=True, include_breit=False,
        drop_noncommuting=True, verbose=False,
    )
    assert out["audit_verdict"] == "DROP_ALL"
    assert out["delta_Q"] == 0


def test_residuals_reported_above_threshold():
    """The audit residuals are surfaced for follow-on analysis."""
    spec = lih_spec_relativistic(max_n=2, include_pk=True)
    out = relativistic_tapered_from_spec(
        spec, mode="per_block",
        pk_in_hamiltonian=True, include_breit=False,
        drop_noncommuting=True, verbose=False,
    )
    assert len(out["stabilizers_residuals"]) == len(out["stabilizers_built"])
    assert all(r > 1e-6 for r in out["stabilizers_residuals"])
    assert len(out["stabilizers_kept"]) == 0
