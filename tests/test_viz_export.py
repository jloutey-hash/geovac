"""
Tests for geovac.viz_export (Phase 0 of docs/visualization_plan.md).

Discipline: every artifact field is cross-checked against a DIRECT call to
the underlying builder (GeometricLattice / AtomicSolver / hamiltonian()),
never against a stored copy of the exporter's own output. Determinism is
byte-level (the CI regeneration check depends on it).
"""

import json
import math

import numpy as np
import pytest

from geovac.lattice import GeometricLattice
from geovac.atomic_solver import AtomicSolver
from geovac.viz_export import (
    SCHEMA_VERSION,
    _dumps,
    _loglog_fit,
    export_all,
    export_lattice,
    export_qc_library,
    export_qc_system,
    export_spectrum,
)


# ---------------------------------------------------------------------------
# Lattice export
# ---------------------------------------------------------------------------

def test_lattice_counts_and_schema():
    max_n = 5
    d = export_lattice(Z=1, max_n=max_n)
    assert d['schema_version'] == SCHEMA_VERSION
    assert d['kind'] == 'lattice'
    # sum_{n=1}^{N} n^2 = N(N+1)(2N+1)/6
    expected_states = max_n * (max_n + 1) * (2 * max_n + 1) // 6
    assert d['num_states'] == expected_states
    assert len(d['states']) == expected_states
    assert len(d['node_weights']) == expected_states
    assert len(d['edges']) == d['num_edges']
    assert sum(s['degeneracy'] for s in d['shells']) == expected_states
    # every edge must be classified by the two construction rules
    types = {e[2] for e in d['edges']}
    assert types <= {'angular', 'radial'}
    # indices in range, upper-triangle, sorted
    for i, j, _t in d['edges']:
        assert 0 <= i < j < expected_states


def test_lattice_edges_match_adjacency():
    max_n, Z = 4, 2
    d = export_lattice(Z=Z, max_n=max_n)
    lat = GeometricLattice(max_n=max_n, nuclear_charge=Z)
    coo = lat.adjacency.tocoo()
    expected = {(int(i), int(j)) for i, j in zip(coo.row, coo.col) if i < j}
    exported = {(e[0], e[1]) for e in d['edges']}
    assert exported == expected
    assert d['num_edges'] == lat.num_edges


def test_lattice_node_weights_cross_check():
    Z, max_n = 3, 4
    d = export_lattice(Z=Z, max_n=max_n)
    for (n, _l, _m), w in zip(d['states'], d['node_weights']):
        assert w == pytest.approx(-Z / n ** 2, rel=1e-14)


# ---------------------------------------------------------------------------
# Spectrum export
# ---------------------------------------------------------------------------

def test_spectrum_ground_state_cross_check():
    Z, max_n = 1, 6
    d = export_spectrum(Z=Z, max_n=max_n)
    solver = AtomicSolver(max_n=max_n, Z=Z)
    e0, _ = solver.compute_ground_state(n_states=1)
    assert min(d['graph_energies_ha']) == pytest.approx(e0[0], rel=1e-8)
    expected_states = max_n * (max_n + 1) * (2 * max_n + 1) // 6
    assert len(d['graph_energies_ha']) == expected_states
    assert len(d['graph_laplacian_spectrum']) == expected_states


def test_spectrum_laplacian_energy_consistency():
    Z, max_n = 2, 5
    d = export_spectrum(Z=Z, max_n=max_n)
    kappa_z2 = -1.0 / 16.0 * Z ** 2
    rescaled = sorted(kappa_z2 * x for x in d['graph_laplacian_spectrum'])
    energies = sorted(d['graph_energies_ha'])
    assert np.allclose(rescaled, energies, rtol=1e-10, atol=1e-12)


def test_spectrum_references_and_series():
    d = export_spectrum(Z=1, max_n=6)
    for lvl in d['s3_reference']['levels']:
        n = lvl['n']
        assert lvl['eigenvalue'] == n * n - 1     # exact integers
        assert lvl['degeneracy'] == n * n
    for lvl in d['hydrogen_levels']:
        n = lvl['n']
        assert lvl['energy_ha'] == pytest.approx(-0.5 / n ** 2, rel=1e-14)
    lyman = next(s for s in d['series'] if s['name'] == 'Lyman')
    alpha = next(l for l in lyman['lines'] if l['n_upper'] == 2)
    # Bohr (infinite nuclear mass) Lyman-alpha: ~121.50 nm
    assert alpha['wavelength_nm'] == pytest.approx(121.50, abs=0.2)
    assert alpha['delta_e_ha'] == pytest.approx(0.375, rel=1e-12)


# ---------------------------------------------------------------------------
# QC system export (requires openfermion)
# ---------------------------------------------------------------------------

def test_qc_system_lih_cross_check():
    pytest.importorskip('openfermion')
    from geovac.ecosystem_export import hamiltonian
    from geovac.measurement_grouping import count_qwc_groups

    d = export_qc_system('LiH', tapered_modes=('per_block',))
    ham = hamiltonian('LiH')
    assert d['n_qubits'] == ham.n_qubits
    assert d['n_terms'] == ham.n_terms
    assert d['one_norm'] == pytest.approx(ham.one_norm, rel=1e-12)
    assert d['has_pk'] == (ham.h1_pk is not None)
    # Pauli list is complete and consistent with the 1-norm
    assert len(d['pauli_terms']) == d['n_terms']
    one_norm_from_terms = sum(math.hypot(re, im)
                              for _s, re, im in d['pauli_terms'])
    assert one_norm_from_terms == pytest.approx(ham.one_norm, rel=1e-10)
    # sorted by descending |coeff|
    mags = [re * re + im * im for _s, re, im in d['pauli_terms']]
    assert mags == sorted(mags, reverse=True)
    # QWC group count matches a direct call
    assert (d['measurement']['n_qwc_groups']
            == count_qwc_groups(ham.to_openfermion()))
    # tapered variant matches a direct tapered build
    row = d['tapering'][0]
    assert row['mode'] == 'per_block'
    tham = hamiltonian('LiH', tapered='per_block')
    assert row['n_qubits'] == tham.n_qubits
    assert row['n_terms'] == tham.n_terms


def test_qc_system_tapering_error_is_captured_not_raised():
    pytest.importorskip('openfermion')
    d = export_qc_system('LiH', tapered_modes=('not_a_mode',),
                         include_pauli_terms=False)
    row = d['tapering'][0]
    assert row['mode'] == 'not_a_mode'
    assert 'error' in row


# ---------------------------------------------------------------------------
# Determinism + serialization
# ---------------------------------------------------------------------------

def test_determinism_byte_level():
    a = _dumps(export_lattice(Z=1, max_n=4))
    b = _dumps(export_lattice(Z=1, max_n=4))
    assert a == b
    a = _dumps(export_spectrum(Z=1, max_n=4))
    b = _dumps(export_spectrum(Z=1, max_n=4))
    assert a == b


def test_qc_determinism_and_json_roundtrip():
    pytest.importorskip('openfermion')
    a = _dumps(export_qc_system('He', tapered_modes=(),
                                include_pauli_terms=True))
    b = _dumps(export_qc_system('He', tapered_modes=(),
                                include_pauli_terms=True))
    assert a == b
    json.loads(a)  # round-trips as valid JSON


def test_loglog_fit_recovers_synthetic_exponent():
    pts = [(q, 3.0 * q ** 2.5) for q in (10, 20, 40, 80)]
    fit = _loglog_fit(pts)
    assert fit['exponent'] == pytest.approx(2.5, rel=1e-10)
    assert fit['prefactor'] == pytest.approx(3.0, rel=1e-10)
    assert fit['rms_log_residual'] == pytest.approx(0.0, abs=1e-12)
    assert fit['n_points'] == 4
    assert 'note' not in fit
    short = _loglog_fit(pts[:3])
    assert 'note' in short  # < 4 points: exponent not quotable as a claim


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def test_export_all_quick_manifest(tmp_path):
    pytest.importorskip('openfermion')
    manifest = export_all(str(tmp_path),
                          lattice_zs=(1,),
                          lattice_max_ns=(2, 3),
                          qc_systems=('LiH',),
                          tapered_modes=(),
                          sweep_max_n_values=(),
                          verbose=False)
    assert manifest['kind'] == 'manifest'
    for rel in manifest['files']:
        p = tmp_path.joinpath(*rel.split('/'))
        assert p.exists()
        json.loads(p.read_text(encoding='utf-8'))
    assert any('library.json' in f for f in manifest['files'])
    assert any('lattice_Z1_n2.json' in f for f in manifest['files'])


def test_benchmark_scaling_block_ingestion():
    from geovac.viz_export import _benchmark_scaling_block
    block = _benchmark_scaling_block()
    assert block is not None, 'benchmarks/qubit_scaling_data.json missing'
    assert block['provenance']['source'] == 'benchmarks/qubit_scaling_data.json'
    he = [p for p in block['points']
          if p['type'] == 'geovac' and p['label'].startswith('GeoVac He')]
    assert len(he) >= 4
    # atomic O(Q^3.15) family (Paper 14 headline; committed artifact)
    assert block['fit_he']['exponent'] == pytest.approx(3.15, abs=0.15)


@pytest.mark.slow
def test_full_library_export_slow():
    pytest.importorskip('openfermion')
    lib = export_qc_library(sweep_max_n_values=(1, 2))
    assert lib['n_systems'] == 37
    assert lib['scaling']['cross_system']['fit']['n_points'] == 37
    names = {r['system'] for r in lib['systems']}
    assert {'LiH', 'H2O', 'He', 'H2', 'FeH', 'BaH'} <= names
    assert len(lib['scaling']['within_system']['points']) == 2
    assert 'atomic_benchmark' in lib['scaling']
