"""
Tests for ecosystem_export.py — full library validation.

Verifies:
  1. All 30 systems build successfully via hamiltonian()
  2. Known Pauli counts match expected values
  3. Universal Pauli/Q ratios hold
  4. Isostructural invariance across rows

Author: GeoVac Development Team
Date: April 2026
"""

import sys
from pathlib import Path

import pytest

# Ensure local project root is on sys.path
_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from geovac.ecosystem_export import hamiltonian


# ---------------------------------------------------------------------------
# Expected Pauli counts (non-identity terms) at max_n=2
# ---------------------------------------------------------------------------
EXPECTED_PAULI = {
    # Atomic / diatomic
    'He':   119,
    'H2':   111,
    # First-row main-group
    'LiH':  333,
    'BeH2': 555,
    'CH4':  999,
    'NH3':  888,
    'H2O':  777,
    'HF':   666,
    # Second-row
    'NaH':  222,
    'MgH2': 444,
    'SiH4': 888,
    'PH3':  777,
    'H2S':  666,
    'HCl':  555,
    # Third-row s-block
    'KH':   222,
    'CaH2': 444,
    # Third-row p-block
    'GeH4': 888,
    'AsH3': 777,
    'H2Se': 666,
    'HBr':  555,
    # TM hydrides (all identical)
    'ScH':  277, 'TiH':  277, 'VH':   277, 'CrH':  277, 'MnH':  277,
    'FeH':  277, 'CoH':  277, 'NiH':  277, 'CuH':  277, 'ZnH':  277,
}


# ---------------------------------------------------------------------------
# Representative pairwise matrix (Sprint Test-Slim 2026-06-07)
# ---------------------------------------------------------------------------
# Covers the dimensions (natural geometry / core treatment / valence-electron
# count / block structure) with ~10 systems instead of 40.  Each entry sits at
# a distinct cell of the cross-product:
#
#   He   — atom, Level 3, 2e, single-block
#   H2   — diatomic, Level 4, 2e, single bond, no heavy-atom core
#   LiH  — first-row explicit core, 4e, single bond
#   BeH2 — multi-center linear, 4e, 2 bonds
#   H2O  — multi-center bent, lone pairs, 8 valence e
#   HF   — first-row, lone pairs, single bond
#   NaH  — frozen-core 2nd row, 2 valence e, single bond
#   HCl  — frozen-core 2nd row, lone pairs, single bond
#   KH   — frozen-core 3rd row, 2 valence e
#   ScH  — TM frozen-core (Ar core), d-orbital block
#
# Full 40-system sweep is the @pytest.mark.slow tests below; run with --slow.
REPRESENTATIVE_PAULI = {
    'He':   119,
    'H2':   111,
    'LiH':  333,
    'BeH2': 555,
    'H2O':  777,
    'HF':   666,
    'NaH':  222,
    'HCl':  555,
    'KH':   222,
    'ScH':  277,
}


@pytest.mark.parametrize("system,expected", list(REPRESENTATIVE_PAULI.items()))
def test_representative_builds_and_pauli(
    hamiltonian_cache, system: str, expected: int,
) -> None:
    """Smoke + Pauli-count exact match on the representative pairwise matrix.

    Single build per system covering both "builds successfully" and
    "Pauli count matches documented value" in one parametrize sweep.
    Uses ``hamiltonian_cache`` (session-scoped) so each system is built
    once across the entire session. The full 40-system sweep is the slow
    tests below.
    """
    H = hamiltonian_cache(system)
    assert H.n_qubits > 0
    assert H.n_terms > 0
    assert H.one_norm > 0
    n_pauli = H.n_terms - 1
    assert n_pauli == expected, (
        f"{system}: expected {expected} non-identity Pauli terms, got {n_pauli}"
    )


# ---------------------------------------------------------------------------
# Full-library coverage (slow; run with `pytest --slow`)
# ---------------------------------------------------------------------------
EXTENDED_PAULI = {
    k: v for k, v in EXPECTED_PAULI.items()
    if k not in REPRESENTATIVE_PAULI
}


@pytest.mark.slow
@pytest.mark.parametrize("system", list(EXTENDED_PAULI.keys()))
def test_system_builds(hamiltonian_cache, system: str) -> None:
    """Each non-representative system should build without error.

    Slim: representative coverage handled by
    ``test_representative_builds_and_pauli``; this fills out the long-tail
    of the 37-system ``_SYSTEM_REGISTRY`` under ``--slow``.
    """
    H = hamiltonian_cache(system)
    assert H.n_qubits > 0
    assert H.n_terms > 0
    assert H.one_norm > 0


@pytest.mark.slow
@pytest.mark.parametrize("system,expected", list(EXTENDED_PAULI.items()))
def test_pauli_counts(
    hamiltonian_cache, system: str, expected: int,
) -> None:
    """Non-identity Pauli count exact match on the non-representative tail."""
    H = hamiltonian_cache(system)
    n_pauli = H.n_terms - 1
    assert n_pauli == expected, (
        f"{system}: expected {expected} non-identity Pauli terms, got {n_pauli}"
    )


# ---------------------------------------------------------------------------
# Test: universal Pauli/Q ratio for main-group hydrides
# ---------------------------------------------------------------------------
MAIN_GROUP = [
    'LiH', 'BeH2', 'CH4', 'NH3', 'H2O', 'HF',
    'NaH', 'MgH2', 'SiH4', 'PH3', 'H2S', 'HCl',
    'KH', 'CaH2', 'GeH4', 'AsH3', 'H2Se', 'HBr',
    'H2',
]

@pytest.mark.slow
@pytest.mark.parametrize("system", MAIN_GROUP)
def test_main_group_ratio(hamiltonian_cache, system: str) -> None:
    """All main-group hydrides should have Pauli/Q = 11.10.

    Marked slow: representative coverage (LiH, BeH2, H2O, HF, NaH, HCl, KH)
    is already in REPRESENTATIVE_PAULI; this fills out the long tail.
    """
    H = hamiltonian_cache(system)
    n_pauli = H.n_terms - 1
    ratio = n_pauli / H.n_qubits
    assert abs(ratio - 11.10) < 0.01, (
        f"{system}: Pauli/Q = {ratio:.2f}, expected 11.10"
    )


# ---------------------------------------------------------------------------
# Test: TM hydride ratio
# ---------------------------------------------------------------------------
TM_HYDRIDES = ['ScH', 'TiH', 'VH', 'CrH', 'MnH',
               'FeH', 'CoH', 'NiH', 'CuH', 'ZnH']

@pytest.mark.slow
@pytest.mark.parametrize("system", TM_HYDRIDES)
def test_tm_ratio(hamiltonian_cache, system: str) -> None:
    """All TM hydrides should have Pauli/Q = 9.23.

    Marked slow: ScH is in REPRESENTATIVE_PAULI as the TM-row representative.
    """
    H = hamiltonian_cache(system)
    n_pauli = H.n_terms - 1
    ratio = n_pauli / H.n_qubits
    assert abs(ratio - 9.23) < 0.01, (
        f"{system}: Pauli/Q = {ratio:.2f}, expected 9.23"
    )


# ---------------------------------------------------------------------------
# Test: isostructural invariance
# ---------------------------------------------------------------------------
ISOSTRUCTURAL_GROUPS = [
    # Same block topology → same Pauli count (and same Q)
    (['NaH', 'KH'], 222, 20),           # alkali hydrides, no core block
    (['MgH2', 'CaH2'], 444, 40),        # alkaline earth hydrides
    (['SiH4', 'GeH4'], 888, 80),        # group 14 tetrahydrides
    (['PH3', 'AsH3'], 777, 70),         # group 15 trihydrides
    (['H2S', 'H2Se'], 666, 60),         # group 16 dihydrides
    (['HCl', 'HBr'], 555, 50),          # group 17 monohydrides
]

@pytest.mark.slow
@pytest.mark.parametrize("systems,expected_pauli,expected_Q", ISOSTRUCTURAL_GROUPS)
def test_isostructural_invariance(
    hamiltonian_cache, systems: list, expected_pauli: int, expected_Q: int,
) -> None:
    """Frozen-core molecules with same block topology must have identical
    Pauli counts and qubit counts across rows.

    Marked slow: representative canary version (NaH vs KH) runs by default
    in ``test_isostructural_invariance_representative``.
    """
    for system in systems:
        H = hamiltonian_cache(system)
        n_pauli = H.n_terms - 1
        assert n_pauli == expected_pauli, (
            f"{system}: N_pauli={n_pauli}, expected {expected_pauli}"
        )
        assert H.n_qubits == expected_Q, (
            f"{system}: Q={H.n_qubits}, expected {expected_Q}"
        )


def test_isostructural_invariance_representative(hamiltonian_cache) -> None:
    """Default-coverage canary for the isostructural invariance claim.

    Verifies one cross-row pair (NaH ↔ KH, alkali hydrides, no core block)
    in the default sweep. The full six-group sweep is the
    ``@pytest.mark.slow`` version above. This single canary is enough to
    catch a regression that breaks the Pauli/Q invariance for frozen-core
    main-group specs.
    """
    expected_pauli = 222
    expected_Q = 20
    for system in ['NaH', 'KH']:
        H = hamiltonian_cache(system)
        n_pauli = H.n_terms - 1
        assert n_pauli == expected_pauli, (
            f"{system}: N_pauli={n_pauli}, expected {expected_pauli}"
        )
        assert H.n_qubits == expected_Q, (
            f"{system}: Q={H.n_qubits}, expected {expected_Q}"
        )


# ---------------------------------------------------------------------------
# Test: TM hydride isostructural invariance (all 10 identical)
# ---------------------------------------------------------------------------
@pytest.mark.slow
def test_tm_isostructural(hamiltonian_cache) -> None:
    """All 10 TM hydrides must have identical Q=30 and N_pauli=277.

    Marked slow: the representative version below covers two TM systems
    (ScH from REPRESENTATIVE_PAULI and TiH) as a canary by default.
    """
    for system in TM_HYDRIDES:
        H = hamiltonian_cache(system)
        n_pauli = H.n_terms - 1
        assert H.n_qubits == 30, f"{system}: Q={H.n_qubits}"
        assert n_pauli == 277, f"{system}: N_pauli={n_pauli}"


def test_tm_isostructural_representative(hamiltonian_cache) -> None:
    """Default-coverage canary for the TM-row invariance claim.

    Verifies ScH ↔ TiH (the first two TM hydrides) match in Q and Pauli.
    Full 10-element sweep is the ``@pytest.mark.slow`` version above.
    """
    for system in ['ScH', 'TiH']:
        H = hamiltonian_cache(system)
        n_pauli = H.n_terms - 1
        assert H.n_qubits == 30, f"{system}: Q={H.n_qubits}"
        assert n_pauli == 277, f"{system}: N_pauli={n_pauli}"


# ---------------------------------------------------------------------------
# Test: propinquity_bound metadata (Paper 38 Thm 1)
# ---------------------------------------------------------------------------
import math


def test_propinquity_bound_finite_and_positive(hamiltonian_cache) -> None:
    """The Paper 38 propinquity bound is finite and positive at max_n=2."""
    H = hamiltonian_cache('LiH')
    bound = H.propinquity_bound
    assert math.isfinite(bound), f"propinquity_bound = {bound} is not finite"
    assert bound > 0.0, f"propinquity_bound = {bound} is not positive"


def test_propinquity_bound_monotone_decreasing_lih(hamiltonian_cache) -> None:
    """gamma_{max_n} is strictly monotone-decreasing in max_n.

    Verified at max_n in {2, 3, 4} on LiH. Calibrated against the
    v3.56.0 ab initio reconciliation (LiH R_eq error 2.82%); the
    bound's monotone decrease is the structural promise of Paper 38
    Thm.~1, separable from per-system accuracy claims.
    """
    bounds = []
    for max_n in [2, 3, 4]:
        H = hamiltonian_cache('LiH', max_n=max_n)
        bounds.append(H.propinquity_bound)
    # Strict monotone decrease.
    assert bounds[0] > bounds[1] > bounds[2], (
        f"LiH propinquity bounds at max_n=2,3,4 = {bounds}; "
        f"expected strictly decreasing."
    )
    # Asymptotic order-of-magnitude sanity: should sit between the
    # asymptotic 4/pi * log(n)/n estimate and the uniform 6 * log(n)/n
    # bound from Paper 38 Thm.~1(ii) for max_n >= 2.
    for max_n, b in zip([2, 3, 4], bounds):
        upper = 6.0 * math.log(max_n) / max_n
        assert b <= upper + 1e-10, (
            f"LiH max_n={max_n}: bound {b:.6f} exceeds uniform "
            f"6*log(n)/n = {upper:.6f}"
        )


def test_propinquity_bound_matches_closed_form_sum_rule(hamiltonian_cache) -> None:
    """The wired bound equals C_3 * gamma_{max_n} bit-identically.

    C_3 = 1 (Paper 38 Lemma L3, sharp at all cutoffs). Direct check
    that the property returns the closed-form sum-rule value.
    """
    from geovac.central_fejer_su2 import gamma_n_via_sum_rule

    for max_n in [2, 3, 4]:
        H = hamiltonian_cache('LiH', max_n=max_n)
        gamma_ref = float(gamma_n_via_sum_rule(max_n, prec=50))
        assert math.isclose(H.propinquity_bound, gamma_ref, rel_tol=1e-12), (
            f"max_n={max_n}: bound {H.propinquity_bound} != gamma_n "
            f"closed form {gamma_ref}"
        )


def test_propinquity_bound_known_values_lih(hamiltonian_cache) -> None:
    """Check known closed-form values at max_n in {2, 3, 4}.

    From CLAUDE.md and Paper 38 Thm.~1 / memory file
    l2_quantitative_rate_4_over_pi.md: gamma_2 ~ 2.0746, gamma_3 ~
    1.6101, gamma_4 ~ 1.3223.
    """
    expected = {2: 2.074551, 3: 1.610060, 4: 1.322333}
    for max_n, gamma_ref in expected.items():
        H = hamiltonian_cache('LiH', max_n=max_n)
        assert math.isclose(H.propinquity_bound, gamma_ref, abs_tol=1e-5), (
            f"max_n={max_n}: bound {H.propinquity_bound:.6f} != "
            f"expected {gamma_ref:.6f}"
        )


def test_propinquity_bound_metadata_fields_present(hamiltonian_cache) -> None:
    """Metadata exposes the bound + provenance fields."""
    H = hamiltonian_cache('LiH')
    meta = H.metadata
    assert 'max_n' in meta
    assert 'propinquity_bound' in meta
    assert 'propinquity_bound_C3' in meta
    assert 'propinquity_bound_asymptotic' in meta
    assert 'propinquity_bound_source' in meta
    # C_3 = 1 exactly (Paper 38 Lemma L3)
    assert meta['propinquity_bound_C3'] == 1.0
    # Asymptotic constant = 4 / pi
    assert math.isclose(meta['propinquity_bound_asymptotic'], 4.0 / math.pi,
                        rel_tol=1e-15)


def test_propinquity_bound_caches(hamiltonian_cache) -> None:
    """Second access uses cached value (no recomputation)."""
    H = hamiltonian_cache('LiH')
    a = H.propinquity_bound
    b = H.propinquity_bound
    assert a == b
    # Cached attribute is populated.
    assert H._cached_propinquity_bound is not None


def test_propinquity_bound_universal_across_systems(hamiltonian_cache) -> None:
    """The bound is a property of max_n alone (universal across
    chemistry systems at the same cutoff).

    Paper 38 Thm.~1 is a statement about the underlying SU(2)
    spectral triple. The bound at max_n = 2 is therefore identical
    across LiH, BeH2, H2O, and H2 — by construction.
    """
    systems = ['LiH', 'BeH2', 'H2O', 'HF', 'H2']
    bounds = [hamiltonian_cache(s).propinquity_bound for s in systems]
    for s, b in zip(systems[1:], bounds[1:]):
        assert math.isclose(bounds[0], b, rel_tol=1e-12), (
            f"{s} bound {b} differs from LiH bound {bounds[0]}"
        )


def test_propinquity_bound_undefined_when_max_n_missing(hamiltonian_cache) -> None:
    """Bound raises ValueError if max_n absent from metadata.

    Guards against silent fallback when constructing a
    GeoVacHamiltonian directly from a hand-rolled QubitOperator.
    """
    from geovac.ecosystem_export import GeoVacHamiltonian
    # Build a dummy 1-qubit operator without max_n metadata.
    H = hamiltonian_cache('LiH')
    qop = H.to_openfermion()
    H_bare = GeoVacHamiltonian(qop, metadata={'system': 'bare'})
    try:
        _ = H_bare.propinquity_bound
        assert False, "Expected ValueError when max_n missing"
    except ValueError as exc:
        assert 'max_n' in str(exc)


# ---------------------------------------------------------------------------
# Sprint P1 (v3.86.0, 2026-06-07): FCIDUMP exporter
# ---------------------------------------------------------------------------
#
# Unblocks the hybrid-pipeline Phase 1 first integration. The integrals
# surfaced here are bit-identical to those that
# ``build_composed_hamiltonian`` feeds into ``build_fermion_op_from_integrals``
# before Jordan-Wigner transformation; a downstream DMRG / CCSD(T) /
# AFQMC consumer that reads the FCIDUMP and performs FCI on the active
# space recovers the same N-electron spectrum as a GeoVac FCI on the
# qubit operator (modulo the qubit-encoding map). The round-trip test
# below uses ``read_fcidump`` (pure-Python, no pyscf dependency); a
# pyscf-based round-trip will be added in Sprint H1-DMRG when block2 /
# pyscf-DMRG are wired in.
# ---------------------------------------------------------------------------

import numpy as np  # noqa: E402  (top-of-section import for the FCIDUMP block)


# Sprint Test-Slim 2026-06-07: the historical helper ``_build_lih_for_fcidump``
# returned ``hamiltonian('LiH', R=3.015, max_n=2)``; that is bit-identical to
# ``hamiltonian_cache('LiH')`` because both kwargs are the function defaults
# (R=_HYDRIDE_REQ['LiH']=3.015, max_n=2). Replaced with the fixture for
# session-wide caching.


def test_lih_pre_jw_integrals_surfaced(hamiltonian_cache) -> None:
    """LiH composed builder must expose (h1, eri, ecore, n_electrons)."""
    H = hamiltonian_cache('LiH')
    assert H.h1 is not None
    assert H.eri is not None
    assert H.ecore is not None
    assert H.n_electrons is not None
    assert H.n_orbitals is not None
    # LiH at max_n=2: 15 spatial orbitals (Li_core + LiH_bond), 4 active e-
    M = H.n_orbitals
    assert M == 15
    assert H.n_electrons == 4
    assert H.h1.shape == (M, M)
    assert H.eri.shape == (M, M, M, M)


def test_lih_fcidump_write_and_parse(hamiltonian_cache, tmp_path) -> None:
    """FCIDUMP write+parse must round-trip bit-exact for LiH."""
    from geovac.ecosystem_export import read_fcidump

    H = hamiltonian_cache('LiH')
    path = str(tmp_path / 'lih.fcidump')
    info = H.to_fcidump(path)

    assert info['n_orbitals'] == H.n_orbitals
    assert info['n_electrons'] == H.n_electrons
    assert info['filename'] == path

    parsed = read_fcidump(path)
    assert parsed['n_orbitals'] == H.n_orbitals
    assert parsed['n_electrons'] == H.n_electrons
    # Bit-exact for h1 (Hermitian symmetry round-trips) and eri (8-fold).
    assert np.max(np.abs(H.h1 - parsed['h1'])) < 1e-12
    assert np.max(np.abs(H.eri - parsed['eri'])) < 1e-12
    assert abs(H.ecore - parsed['ecore']) < 1e-12


def test_fcidump_eri_eight_fold_symmetry(hamiltonian_cache, tmp_path) -> None:
    """Written FCIDUMP must satisfy chemist-notation eight-fold symmetry
    on round-trip: (pq|rs) = (qp|rs) = (pq|sr) = (rs|pq) etc."""
    from geovac.ecosystem_export import read_fcidump

    H = hamiltonian_cache('LiH')
    path = str(tmp_path / 'lih.fcidump')
    H.to_fcidump(path)
    eri = read_fcidump(path)['eri']
    M = eri.shape[0]
    # Spot-check the 8-fold symmetry on a few representative quartets.
    for (p, q, r, s) in [(0, 0, 0, 0), (0, 1, 0, 1), (1, 2, 3, 4),
                          (5, 5, 6, 6), (2, 3, 4, 5)]:
        if not (p < M and q < M and r < M and s < M):
            continue
        v = eri[p, q, r, s]
        assert abs(eri[q, p, r, s] - v) < 1e-12
        assert abs(eri[p, q, s, r] - v) < 1e-12
        assert abs(eri[q, p, s, r] - v) < 1e-12
        assert abs(eri[r, s, p, q] - v) < 1e-12
        assert abs(eri[s, r, p, q] - v) < 1e-12
        assert abs(eri[r, s, q, p] - v) < 1e-12
        assert abs(eri[s, r, q, p] - v) < 1e-12


def test_fcidump_header_compliance(hamiltonian_cache, tmp_path) -> None:
    """FCIDUMP header must contain &FCI namelist with NORB/NELEC/MS2/ORBSYM/ISYM."""
    H = hamiltonian_cache('LiH')
    path = str(tmp_path / 'lih.fcidump')
    H.to_fcidump(path)
    with open(path, 'r', encoding='ascii') as fh:
        text = fh.read()
    # &FCI ... &END
    assert '&FCI' in text
    assert '&END' in text
    # Required fields
    assert 'NORB=' in text
    assert 'NELEC=' in text
    assert 'MS2=' in text
    assert 'ORBSYM=' in text
    assert 'ISYM=' in text


def test_fcidump_one_body_diagonal_count(hamiltonian_cache) -> None:
    """LiH h1 at max_n=2 is diagonal in the hydrogenic basis; only M nonzeros."""
    H = hamiltonian_cache('LiH')
    # GeoVac h1 = -Z^2 / (2 n^2) on diagonal (single-block hydrogenic
    # eigenvalues). All off-diagonals are structurally zero before
    # cross-block h1 is enabled. The diagonal count matches M.
    h1 = H.h1
    M = H.n_orbitals
    n_diag_nonzero = int(np.sum(np.abs(np.diag(h1)) > 1e-12))
    n_offdiag_nonzero = int(
        np.sum(np.abs(h1 - np.diag(np.diag(h1))) > 1e-12)
    )
    assert n_diag_nonzero == M
    # Off-diagonal h1 is zero for single-center spec (no cross-block h1
    # enabled by default).
    assert n_offdiag_nonzero == 0


def test_fcidump_ecore_matches_nuclear_repulsion(hamiltonian_cache) -> None:
    """ecore wired to spec.nuclear_repulsion_constant; matches builder output."""
    from geovac.molecular_spec import hydride_spec
    from geovac.composed_qubit import build_composed_hamiltonian

    spec = hydride_spec(3, R=3.015, max_n=2, core_method='pk')
    result = build_composed_hamiltonian(
        spec, pk_in_hamiltonian=False, verbose=False,
    )
    H = hamiltonian_cache('LiH')
    assert abs(H.ecore - result['nuclear_repulsion']) < 1e-12


def test_fcidump_works_for_multiple_systems(hamiltonian_cache, tmp_path) -> None:
    """Round-trip a representative sample across hydride / multi-center / TM."""
    from geovac.ecosystem_export import read_fcidump

    for system in ['LiH', 'BeH2', 'NaH', 'H2O', 'CO', 'ScH', 'H2']:
        H = hamiltonian_cache(system)
        assert H.h1 is not None, f"{system}: h1 not surfaced"
        path = str(tmp_path / f'{system}.fcidump')
        info = H.to_fcidump(path)
        parsed = read_fcidump(path)
        assert parsed['n_orbitals'] == info['n_orbitals'], (
            f"{system}: NORB mismatch on round-trip"
        )
        # eri max diff under 1e-10 (we write 16 digits)
        assert np.max(np.abs(H.eri - parsed['eri'])) < 1e-10, (
            f"{system}: eri round-trip exceeded 1e-10"
        )
        assert np.max(np.abs(H.h1 - parsed['h1'])) < 1e-10, (
            f"{system}: h1 round-trip exceeded 1e-10"
        )


def test_fcidump_tolerance_filters_small(hamiltonian_cache, tmp_path) -> None:
    """tol parameter must drop integrals below the threshold."""
    from geovac.ecosystem_export import read_fcidump

    H = hamiltonian_cache('LiH')
    path_tight = str(tmp_path / 'tight.fcidump')
    path_loose = str(tmp_path / 'loose.fcidump')
    info_tight = H.to_fcidump(path_tight, tol=1e-14)
    info_loose = H.to_fcidump(path_loose, tol=1e-2)
    # Loose tolerance writes strictly fewer terms.
    assert info_loose['n_two_body_terms'] <= info_tight['n_two_body_terms']
    # Both files must round-trip (loose just loses sub-threshold detail).
    parsed_loose = read_fcidump(path_loose)
    assert parsed_loose['n_orbitals'] == info_tight['n_orbitals']


def test_fcidump_raises_when_integrals_missing(hamiltonian_cache) -> None:
    """He uses the LatticeIndex path which does not surface (h1, eri);
    to_fcidump should raise a clear ValueError."""
    H = hamiltonian_cache('He')
    # He intentionally has no h1/eri surfaced at this time.
    if H.h1 is not None:
        # He path was upgraded; this test becomes vacuous.
        return
    try:
        H.to_fcidump('/tmp/nonexistent.fcidump')
        assert False, "Expected ValueError when h1/eri missing"
    except ValueError as exc:
        assert 'h1' in str(exc) or 'eri' in str(exc)


def test_fcidump_orbsym_default_and_custom(hamiltonian_cache, tmp_path) -> None:
    """orbsym defaults to all 1s; custom orbsym round-trips."""
    from geovac.ecosystem_export import read_fcidump

    H = hamiltonian_cache('LiH')
    M = H.n_orbitals

    # Default
    path_default = str(tmp_path / 'default.fcidump')
    H.to_fcidump(path_default)
    parsed_default = read_fcidump(path_default)
    assert parsed_default['orbsym'] == [1] * M

    # Custom
    custom = [2] * M
    path_custom = str(tmp_path / 'custom.fcidump')
    H.to_fcidump(path_custom, orbsym=custom)
    parsed_custom = read_fcidump(path_custom)
    assert parsed_custom['orbsym'] == custom


def test_fcidump_h1_hermitian_round_trip(hamiltonian_cache, tmp_path) -> None:
    """Writer emits only p>=q; reader must expand to full Hermitian h1."""
    from geovac.ecosystem_export import read_fcidump

    # Pick a system where h1 has nontrivial off-diagonals; NH3 has the
    # standard composed h1 shape with the bond-block hopping terms.
    H = hamiltonian_cache('NH3')
    path = str(tmp_path / 'nh3.fcidump')
    H.to_fcidump(path)
    parsed = read_fcidump(path)
    # Hermitian symmetry of round-tripped h1.
    h1 = parsed['h1']
    assert np.max(np.abs(h1 - h1.T)) < 1e-12
