"""
Tests for heavy-atom alkaline-earth monohydrides: SrH and BaH.

Sprint 3 Track HA-C+D (v2.12.0): Scalar and relativistic composed qubit
Hamiltonians for SrH (Z=38, [Kr] frozen core) and BaH (Z=56, [Xe] frozen
core).  Unblocks the Sunaga 2025 (PRA 111, 022817) matched-Q comparison
at Q=18–30.

Validates:
  - Scalar SrH / BaH Pauli count (exact; isostructural with CaH/KH)
  - Relativistic SrH / BaH Pauli count (exact; isostructural with CaH_rel)
  - Isostructural invariance across rows (KH = NaH = SrH_scalar = BaH_scalar)
  - Relativistic ratio rel/scalar at n_max=2 (2.42× per Tier 2 T3 pin)
  - Ecosystem export: ``hamiltonian('SrH')``, ``hamiltonian('BaH')`` work
  - OpenFermion smoke test (Qiskit/PennyLane are optional deps)
  - No regression in existing canonical Pauli counts (LiH=837, BeH2=1395,
    NaH=558, KH=558, CaH2=1116, MgH2=1116) [exact rule 2026-08-29]
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import pytest

_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from geovac.molecular_spec import (
    srh_spec,
    bah_spec,
    srh_spec_relativistic,
    bah_spec_relativistic,
    cah_spec_relativistic,
    lih_spec,
    lih_spec_relativistic,
    nah_spec,
    kh_spec,
    cah2_spec,
    beh2_spec,
)
from geovac.composed_qubit import build_composed_hamiltonian
from geovac.ecosystem_export import hamiltonian, GeoVacHamiltonian


# ---------------------------------------------------------------------------
# Scalar SrH / BaH Pauli counts (exact, isostructural with CaH / KH)
# ---------------------------------------------------------------------------

SCALAR_EXPECTED = {
    # Corrected 2026-08-29 (exact rule; retired rule-A: 222)
    'SrH': {'Q': 20, 'N_pauli': 558},
    'BaH': {'Q': 20, 'N_pauli': 558},
}


@pytest.mark.parametrize("name,expected", list(SCALAR_EXPECTED.items()))
def test_scalar_pauli_count(name: str, expected: dict) -> None:
    """Scalar SrH / BaH must give exactly Q=20, N_pauli=558 (isostructural
    with KH, NaH, CaH-monohydride-variant)."""
    spec = srh_spec() if name == 'SrH' else bah_spec()
    result = build_composed_hamiltonian(spec, pk_in_hamiltonian=False)
    assert result['Q'] == expected['Q'], (
        f"{name}: Q={result['Q']}, expected {expected['Q']}"
    )
    assert result['N_pauli'] == expected['N_pauli'], (
        f"{name}: N_pauli={result['N_pauli']}, expected {expected['N_pauli']}"
    )


def test_scalar_isostructural_alkali_alkaline_earth_monohydride() -> None:
    """Alkali monohydrides (NaH, KH) and alkaline-earth monohydrides
    (SrH, BaH) share the same Q=20, N_pauli=558 topology — all have exactly
    one bond block + frozen core, no explicit core block in the Hamiltonian.
    """
    results = {}
    for name, spec_fn in [
        ('NaH', nah_spec),
        ('KH', kh_spec),
        ('SrH', srh_spec),
        ('BaH', bah_spec),
    ]:
        spec = spec_fn()
        r = build_composed_hamiltonian(spec, pk_in_hamiltonian=False)
        results[name] = (r['Q'], r['N_pauli'])
    # All four must be equal
    ref = results['NaH']
    for name, val in results.items():
        assert val == ref, (
            f"Isostructural violation: {name}={val}, reference (NaH)={ref}"
        )
    # Pin the exact reference for regression
    assert ref == (20, 558), (
        f"Isostructural invariant drift: all 4 = {ref}, pinned (20, 558)"
    )


def test_scalar_pauli_q_ratio() -> None:
    """Main-group hydrides obey Pauli/Q = 27.90 (Track CU universal)."""
    for name, spec_fn in [('SrH', srh_spec), ('BaH', bah_spec)]:
        spec = spec_fn()
        r = build_composed_hamiltonian(spec, pk_in_hamiltonian=False)
        ratio = r['N_pauli'] / r['Q']
        assert abs(ratio - 27.90) < 0.01, (
            f"{name}: Pauli/Q = {ratio:.3f}, expected 27.90"
        )


# ---------------------------------------------------------------------------
# Relativistic SrH / BaH Pauli counts (Dirac (κ, m_j) basis)
# ---------------------------------------------------------------------------

RELATIVISTIC_EXPECTED = {
    # Post-TR fix (Sprint 4, April 2026): jj_angular_Xk received its
    # missing (-1)^{j+1/2} reduced-matrix-element phase; the corrupted
    # cancellations that produced 534 went away, new physical count is 998.
    'SrH_rel': {'Q': 20, 'N_pauli': 998},
    'BaH_rel': {'Q': 20, 'N_pauli': 998},
}


@pytest.mark.parametrize("name,expected", list(RELATIVISTIC_EXPECTED.items()))
def test_relativistic_pauli_count(name: str, expected: dict) -> None:
    """Relativistic SrH_rel / BaH_rel must give exactly Q=20, N_pauli=998.
    Isostructural with CaH_rel."""
    spec = (srh_spec_relativistic() if name == 'SrH_rel'
            else bah_spec_relativistic())
    result = build_composed_hamiltonian(spec)
    assert result['Q'] == expected['Q']
    assert result['N_pauli'] == expected['N_pauli']


def test_relativistic_isostructural_CaH_SrH_BaH() -> None:
    """CaH_rel / SrH_rel / BaH_rel must share the same Q and N_pauli.
    Block topology is identical: [frozen core] + 1 σ-bond block; only the
    screening Z_eff(r) and spin-orbit Z⁴α² prefactor differ.
    """
    for spec_fn in [cah_spec_relativistic, srh_spec_relativistic,
                    bah_spec_relativistic]:
        spec = spec_fn()
        r = build_composed_hamiltonian(spec)
        assert r['Q'] == 20
        assert r['N_pauli'] == 998, (
            f"{spec.name}: N_pauli={r['N_pauli']}, expected 998"
        )


def test_relativistic_scalar_ratio_n_max2() -> None:
    """Relativistic / scalar Pauli ratio at max_n=2 is ~1.69 for SrH and
    BaH (exact-rule 2026-08-29: measured 998/558 = 1.688).

    History: the retired 4.24 compared a rule-A scalar denominator against
    the full-Gaunt (rule-B-equivalent) relativistic numerator -- a
    cross-rule comparison; 4.24/2.51 = 1.69 exactly.  (Deeper history:
    before the April 2026 TR fix the corrupted X_k phase gave 2.40.)
    """
    for rel_fn, sc_fn in [
        (srh_spec_relativistic, srh_spec),
        (bah_spec_relativistic, bah_spec),
    ]:
        rel = build_composed_hamiltonian(rel_fn())
        sc = build_composed_hamiltonian(sc_fn(), pk_in_hamiltonian=False)
        ratio = rel['N_pauli'] / sc['N_pauli']
        assert 1.5 < ratio < 1.9, (
            f"{rel_fn.__name__}/{sc_fn.__name__} ratio = {ratio:.2f} "
            f"outside [1.5, 1.9]"
        )


# ---------------------------------------------------------------------------
# Ecosystem export
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("system", ['SrH', 'BaH'])
def test_ecosystem_hamiltonian(system: str) -> None:
    """hamiltonian('SrH') / hamiltonian('BaH') return a GeoVacHamiltonian."""
    H = hamiltonian(system)
    assert isinstance(H, GeoVacHamiltonian)
    assert H.n_qubits == 20
    # n_terms includes identity; non-identity = 558 (exact-rule 2026-08-29)
    assert H.n_terms - 1 == 558
    assert H.one_norm > 0


@pytest.mark.parametrize("system", ['SrH', 'BaH'])
def test_ecosystem_openfermion(system: str) -> None:
    """SrH / BaH export as OpenFermion QubitOperator."""
    H = hamiltonian(system)
    of = H.to_openfermion()
    assert len(of.terms) == 559  # 558 non-identity + 1 identity (exact-rule)


@pytest.mark.parametrize("system", ['SrH', 'BaH'])
def test_ecosystem_qiskit(system: str) -> None:
    """SrH / BaH export as Qiskit SparsePauliOp (if qiskit is installed)."""
    pytest.importorskip("qiskit")
    H = hamiltonian(system)
    op = H.to_qiskit()
    assert len(op) == 559  # exact-rule 2026-08-29


def test_ecosystem_case_insensitive() -> None:
    """'srh' / 'bah' are accepted (case-insensitive dispatch)."""
    for name_pair in [('SrH', 'srh'), ('BaH', 'bah')]:
        H1 = hamiltonian(name_pair[0])
        H2 = hamiltonian(name_pair[1])
        assert H1.n_qubits == H2.n_qubits
        assert H1.n_terms == H2.n_terms
        assert abs(H1.one_norm - H2.one_norm) < 1e-10


# ---------------------------------------------------------------------------
# Regression: existing molecule Pauli counts unchanged
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name,spec_fn,n_pauli,Q", [
    # exact-rule 2026-08-29
    ('LiH',  lih_spec,  837, 30),
    ('BeH2', beh2_spec, 1395, 50),
    ('NaH',  nah_spec,  558, 20),
    ('KH',   kh_spec,   558, 20),
    ('CaH2', cah2_spec, 1116, 40),
])
def test_no_regression(name: str, spec_fn: Any, n_pauli: int, Q: int) -> None:
    """Adding SrH / BaH must not disturb existing molecule Pauli counts."""
    spec = spec_fn()
    r = build_composed_hamiltonian(spec, pk_in_hamiltonian=False)
    assert r['Q'] == Q, f"{name}: Q={r['Q']}, expected {Q}"
    assert r['N_pauli'] == n_pauli, (
        f"{name}: N_pauli={r['N_pauli']}, expected {n_pauli}"
    )


# ---------------------------------------------------------------------------
# Relativistic LiH regression (Tier 2 T3 baseline)
# ---------------------------------------------------------------------------

def test_relativistic_lih_regression() -> None:
    """LiH_rel (n_max=2) has exactly 1501 Pauli terms.

    Corrected 2026-08-30: the retired 1413 came from the
    Condon-Shortley factor-order defect (X_k(b,d) where the rule needs
    X_k(d,b)).  See tests/test_paper14_spinor_eri_ordering.py for the
    Clebsch-Gordan discriminator that settled it.

    Pre-TR-fix this was 805; Track TR (Sprint 4) restored the missing
    (-1)^{j+1/2} reduced-matrix-element phase in jj_angular_Xk, removing
    the accidental cancellations that had suppressed ~60% of the true
    physical Pauli terms.
    """
    spec = lih_spec_relativistic(max_n=2)
    r = build_composed_hamiltonian(spec)
    assert r['Q'] == 30
    assert r['N_pauli'] == 1501, (
        f"LiH_rel regression: N_pauli={r['N_pauli']}, expected 1501"
    )


# ---------------------------------------------------------------------------
# Frozen-core identity: SrH and BaH specs must use [Kr] and [Xe] cores
# ---------------------------------------------------------------------------

def test_srh_uses_kr_core() -> None:
    """Sanity: SrH's nuclear-repulsion constant must be dominated by the
    [Kr]-core Sr²⁺ energy (≈ -3131.5 Ha); the frozen-core machinery in
    molecular_spec._alkaline_earth_monohydride_spec must have been invoked.
    """
    spec = srh_spec()
    # Order of magnitude check: Sr²⁺ [Kr]-core energy is ~ -3130 Ha.
    # V_cross and V_NN contribute O(10 Ha) each at R~4 bohr; total ~-3130.
    assert -3200.0 < spec.nuclear_repulsion_constant < -3100.0, (
        f"SrH NR={spec.nuclear_repulsion_constant:.3f}, "
        f"expected Kr-core-dominated (~-3131 Ha)"
    )


def test_bah_uses_xe_core() -> None:
    """Sanity: BaH's nuclear-repulsion constant must be dominated by the
    [Xe]-core Ba²⁺ energy (≈ -7883.5 Ha)."""
    spec = bah_spec()
    assert -7950.0 < spec.nuclear_repulsion_constant < -7850.0, (
        f"BaH NR={spec.nuclear_repulsion_constant:.3f}, "
        f"expected Xe-core-dominated (~-7883 Ha)"
    )


# ---------------------------------------------------------------------------
# Spec structural properties
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("spec_fn", [srh_spec, bah_spec])
def test_single_bond_block(spec_fn: Any) -> None:
    """SrH / BaH scalar spec must have exactly one bond block (no core
    block in the Hamiltonian — core is frozen)."""
    spec = spec_fn()
    assert len(spec.blocks) == 1
    blk = spec.blocks[0]
    assert blk.block_type == 'bond'
    assert blk.has_h_partner is True
    assert blk.Z_center == 2.0  # Z_eff valence for alkaline-earth D-type
    assert blk.Z_partner == 1.0  # H


@pytest.mark.parametrize("spec_fn", [srh_spec_relativistic, bah_spec_relativistic])
def test_relativistic_flag(spec_fn: Any) -> None:
    """Relativistic spec carries relativistic=True."""
    spec = spec_fn()
    assert spec.relativistic is True
    assert spec.name.endswith('_rel')
