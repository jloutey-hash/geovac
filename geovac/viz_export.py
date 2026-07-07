"""
Visualization Data Export — precomputed JSON artifacts for the GeoVac web viz
==============================================================================

Phase 0 of ``docs/visualization_plan.md``. Pure-read module: it imports the
existing builders (``GeometricLattice``, ``AtomicSolver``,
``ecosystem_export.hamiltonian``), serializes their output to versioned JSON,
and modifies nothing.

Single-source-of-truth rule (visualization_plan.md §2): the web layer renders
these artifacts and never hardcodes physics numbers. Every number shown on the
site must trace back to a file generated here, and every file generated here
must trace back to a tested ``geovac`` builder.

Determinism contract: output is byte-identical across runs (sorted keys,
stable orderings, no timestamps), so regeneration diffs are meaningful and CI
can compare committed artifacts against a fresh regeneration.

Usage::

    python -m geovac.viz_export --out viz/public/data
    python -m geovac.viz_export --out viz/public/data --quick   # smoke subset

Author: GeoVac Development Team
Date: July 2026 (Phase 0, visualization plan)
"""

from __future__ import annotations

import argparse
import json
import os
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from geovac import __version__ as _GEOVAC_VERSION
from geovac.lattice import GeometricLattice
from geovac.atomic_solver import AtomicSolver

SCHEMA_VERSION = 1

# Physical constants (CODATA 2018), used only to express exported energy
# differences as photon wavelengths for the M3 spectrum rendering.
HARTREE_EV = 27.211386245988   # Ha -> eV
HC_EV_NM = 1239.841984         # h*c in eV*nm (photon energy <-> vacuum wavelength)

# Universal topological kinetic scale (CLAUDE.md §8 sanctioned literal).
# Tier: Observation — matching constant, not derived (docs/claims_register.md).
KINETIC_SCALE = -1.0 / 16.0

_KAPPA_TIER_NOTE = (
    "Observation: the kinetic scale -1/16 that maps graph eigenvalues to "
    "physical energies is a matching constant. It coincides numerically "
    "with a geometric factor of 1/16 arising in Fock's 1935 stereographic "
    "projection, but no derivation connects the two - the match is an "
    "observation, not a derived result."
)

_S3_REFERENCE_NOTE = (
    "Continuum reference: eigenvalues n^2-1 (degeneracy n^2) are the "
    "Laplace-Beltrami spectrum on the unit three-sphere, the continuum "
    "counterpart of the finite graph. They are NOT the eigenvalues of the "
    "finite graph Laplacian shown alongside."
)

# Default tapering modes exported per QC system (ecosystem_export.hamiltonian
# `tapered` kwarg vocabulary).
DEFAULT_TAPERED_MODES: Tuple[str, ...] = ('global', 'per_block', 'extended', 'full')


# ---------------------------------------------------------------------------
# Envelope / serialization helpers
# ---------------------------------------------------------------------------

def _envelope(kind: str) -> Dict[str, Any]:
    """Common header for every exported artifact (no timestamp: determinism)."""
    return {
        'schema_version': SCHEMA_VERSION,
        'geovac_version': _GEOVAC_VERSION,
        'generator': 'geovac.viz_export',
        'kind': kind,
    }


def _sanitize(obj: Any) -> Any:
    """Recursively convert numpy/tuple/set values into plain JSON types."""
    if isinstance(obj, dict):
        return {str(k): _sanitize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set, frozenset)):
        return [_sanitize(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return [_sanitize(v) for v in obj.tolist()]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, complex):
        return {'re': obj.real, 'im': obj.imag}
    return obj


def _dumps(obj: Dict[str, Any]) -> str:
    """Deterministic JSON encoding (sorted keys, compact, trailing newline)."""
    return json.dumps(_sanitize(obj), sort_keys=True, separators=(',', ':')) + '\n'


def _write_json(obj: Dict[str, Any], path: str) -> str:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write(_dumps(obj))
    return path


# ---------------------------------------------------------------------------
# M1: lattice + spectrum exports
# ---------------------------------------------------------------------------

def _edge_type(state_a: Tuple[int, int, int], state_b: Tuple[int, int, int]) -> str:
    """Classify a lattice edge by the quantum numbers it connects.

    'angular' = m <-> m+-1 within the same (n, l);
    'radial'  = n <-> n+-1 within the same (l, m).
    Any other pair would violate the lattice construction rules.
    """
    n1, l1, m1 = state_a
    n2, l2, m2 = state_b
    if n1 == n2 and l1 == l2 and abs(m1 - m2) == 1:
        return 'angular'
    if l1 == l2 and m1 == m2 and abs(n1 - n2) == 1:
        return 'radial'
    return 'unknown'


def export_lattice(Z: int = 1, max_n: int = 6) -> Dict[str, Any]:
    """Export the (n, l, m) quantum graph: nodes, typed edges, node weights."""
    lat = GeometricLattice(max_n=max_n, nuclear_charge=Z)

    edges: List[List[Any]] = []
    coo = lat.adjacency.tocoo()
    for i, j in zip(coo.row, coo.col):
        if i < j:  # upper triangle only; graph is undirected
            edges.append([int(i), int(j),
                          _edge_type(lat.states[i], lat.states[j])])
    edges.sort(key=lambda e: (e[0], e[1]))

    out = _envelope('lattice')
    out.update({
        'Z': Z,
        'max_n': max_n,
        'states': [[int(n), int(l), int(m)] for (n, l, m) in lat.states],
        'edges': edges,
        'node_weights': [float(w) for w in lat.node_weights],
        'node_weight_rule': '-Z/n^2',
        'edge_weighting': 'binary',
        'shells': [{'n': n, 'degeneracy': n * n}
                   for n in range(1, max_n + 1)],
        'num_states': lat.num_states,
        'num_edges': lat.num_edges,
    })
    return out


def export_spectrum(Z: int = 1, max_n: int = 6) -> Dict[str, Any]:
    """Export the finite-graph spectrum next to its continuum references.

    Three clearly-separated blocks (the distinction is load-bearing — the
    kappa Observation tier forbids conflating them):

    - ``graph_laplacian_spectrum``: actual eigenvalues of D - A for the
      finite lattice (dimensionless, numerical).
    - ``graph_energies``: eigenvalues of H = kappa * Z^2 * (D - A), the
      ``AtomicSolver`` object (Hartree).
    - ``s3_reference``: the continuum S^3 Laplace-Beltrami integers n^2 - 1.
    - ``hydrogen_levels`` / ``series``: the physical Bohr spectrum and
      photon wavelengths (standard formulas, for M3 rendering).
    """
    solver = AtomicSolver(max_n=max_n, Z=Z)
    h_dense = solver.H.toarray()
    graph_energies = np.linalg.eigvalsh(h_dense)
    # H = kinetic_scale * (D - A) with kinetic_scale = KINETIC_SCALE * Z^2,
    # so the dimensionless Laplacian spectrum is H / solver.kinetic_scale.
    laplacian_spectrum = np.sort(np.linalg.eigvalsh(
        h_dense / solver.kinetic_scale))

    series_names = {1: 'Lyman', 2: 'Balmer', 3: 'Paschen', 4: 'Brackett'}
    series: List[Dict[str, Any]] = []
    for n1 in range(1, min(max_n, 4) + 1):
        lines = []
        for n2 in range(n1 + 1, max_n + 1):
            delta_e_ha = (Z ** 2) * 0.5 * (1.0 / n1 ** 2 - 1.0 / n2 ** 2)
            lines.append({
                'n_upper': n2,
                'n_lower': n1,
                'delta_e_ha': delta_e_ha,
                'wavelength_nm': HC_EV_NM / (delta_e_ha * HARTREE_EV),
            })
        series.append({'name': series_names.get(n1, f'n_lower={n1}'),
                       'n_lower': n1, 'lines': lines})

    out = _envelope('spectrum')
    out.update({
        'Z': Z,
        'max_n': max_n,
        'graph_laplacian_spectrum': [float(x) for x in laplacian_spectrum],
        'graph_energies_ha': [float(x) for x in graph_energies],
        'kinetic_scale': {'value': KINETIC_SCALE,
                          'z_scaling': 'kappa * Z^2',
                          'note': _KAPPA_TIER_NOTE},
        's3_reference': {
            'note': _S3_REFERENCE_NOTE,
            'levels': [{'n': n, 'eigenvalue': n * n - 1, 'degeneracy': n * n}
                       for n in range(1, max_n + 1)],
        },
        'hydrogen_levels': [{'n': n,
                             'energy_ha': -(Z ** 2) / (2.0 * n ** 2),
                             'degeneracy': n * n}
                            for n in range(1, max_n + 1)],
        'series': series,
    })
    return out


# ---------------------------------------------------------------------------
# M2: qubit-Hamiltonian exports
# ---------------------------------------------------------------------------

def _pauli_term_key(term: Tuple[Tuple[int, str], ...]) -> str:
    """Render an OpenFermion term key as 'X0 Z3' ('I' for the identity)."""
    if not term:
        return 'I'
    return ' '.join(f'{p}{q}' for q, p in term)


def _pauli_terms_list(qubit_op: Any) -> List[List[Any]]:
    """Serialize a QubitOperator as [[pauli_string, re, im], ...],
    sorted by descending |coefficient| (ties by string) for determinism."""
    rows = []
    for term, coeff in qubit_op.terms.items():
        c = complex(coeff)
        rows.append([_pauli_term_key(term), c.real, c.imag])
    rows.sort(key=lambda r: (-(r[1] ** 2 + r[2] ** 2), r[0]))
    return rows


def export_qc_system(
    name: str,
    tapered_modes: Sequence[str] = DEFAULT_TAPERED_MODES,
    include_pauli_terms: bool = True,
    **build_kwargs: Any,
) -> Dict[str, Any]:
    """Export one library system: resource summary, Pauli terms, QWC
    measurement analysis, and tapering variants.

    Tapering modes that are inapplicable to a system are recorded as
    ``{'mode': ..., 'error': ...}`` rows rather than raised, so the library
    export never fails on a per-mode gap.
    """
    from geovac.ecosystem_export import hamiltonian
    from geovac.measurement_grouping import analyze_measurement_cost

    ham = hamiltonian(name, **build_kwargs)
    qubit_op = ham.to_openfermion()
    meas = analyze_measurement_cost(qubit_op)

    try:
        metadata = _sanitize(ham.metadata)
    except Exception as exc:  # metadata property may probe propinquity
        metadata = {'metadata_error': f'{type(exc).__name__}: {exc}'}

    tapering: List[Dict[str, Any]] = []
    for mode in tapered_modes:
        try:
            tham = hamiltonian(name, tapered=mode, **build_kwargs)
            tapering.append({
                'mode': mode,
                'n_qubits': tham.n_qubits,
                'n_terms': tham.n_terms,
                'one_norm': tham.one_norm,
                'metadata': _sanitize({
                    k: v for k, v in tham.metadata.items()
                    if k.startswith(('tapered', 'Q_tapered', 'delta_Q',
                                     'n_sub_blocks', 'n_stabs', 'kinds',
                                     'dropped'))
                }),
            })
        except Exception as exc:
            tapering.append({'mode': mode,
                             'error': f'{type(exc).__name__}: {exc}'})

    out = _envelope('qc_system')
    out.update({
        'system': name,
        'n_qubits': ham.n_qubits,
        'n_terms': ham.n_terms,
        'one_norm': ham.one_norm,
        'one_norm_full': ham.one_norm_full,
        'has_pk': ham.h1_pk is not None,
        'n_electrons': ham.n_electrons,
        'n_orbitals': ham.n_orbitals,
        'ecore': ham.ecore,
        'metadata': metadata,
        'measurement': {
            'n_qwc_groups': meas.n_qwc_groups,
            'max_group_size': meas.max_group_size,
            'min_group_size': meas.min_group_size,
            'mean_group_size': meas.mean_group_size,
        },
        'tapering': tapering,
    })
    if include_pauli_terms:
        out['pauli_terms'] = _pauli_terms_list(qubit_op)
    return out


def _loglog_fit(points: List[Tuple[float, float]]) -> Dict[str, Any]:
    """Log-log least-squares fit y = a * x^b with residual reporting
    (CLAUDE.md §13.4a scaling-law rule: >= 4 points, report the residual)."""
    xs = np.log([p[0] for p in points])
    ys = np.log([p[1] for p in points])
    n_points = len(points)
    if n_points < 2:
        return {'n_points': n_points, 'note': 'too few points for a fit'}
    slope, intercept = np.polyfit(xs, ys, 1)
    resid = ys - (slope * xs + intercept)
    fit: Dict[str, Any] = {
        'exponent': float(slope),
        'prefactor': float(np.exp(intercept)),
        'rms_log_residual': float(np.sqrt(np.mean(resid ** 2))),
        'n_points': n_points,
    }
    if n_points < 4:
        fit['note'] = ('fewer than 4 points - shown as points only; too '
                       'few to quote a fitted scaling exponent')
    return fit


def _benchmark_scaling_block() -> Optional[Dict[str, Any]]:
    """Ingest the committed atomic-scaling benchmark artifact
    (``benchmarks/qubit_scaling_data.json``, generated by
    ``benchmarks/qubit_scaling_sweep.py``) as a labeled provenance block.

    This is the deep-sweep family (He/H atomic, max_n = 2..5, plus the
    Gaussian reference rows) behind the atomic O(Q^3.15) claim — too
    expensive to recompute at export time (the max_n = 5 build runs for
    minutes), so the canonical committed artifact is re-served with its
    provenance attached. Returns None if the file is absent (e.g. an
    installed package without the benchmarks tree).
    """
    path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'benchmarks', 'qubit_scaling_data.json')
    if not os.path.exists(path):
        return None
    with open(path, 'r', encoding='utf-8') as fh:
        data = json.load(fh)

    rows = [{'label': s.get('label'), 'type': s.get('type'),
             'max_n': s.get('max_n'), 'Q': s.get('Q'),
             'n_pauli': s.get('n_pauli_terms'),
             'one_norm': s.get('one_norm')}
            for s in data.get('systems', [])]
    block: Dict[str, Any] = {
        'label': ('atomic within-system family (He/H, max_n sweep) + '
                  'Gaussian references; the family behind the atomic '
                  'O(Q^3.15) claim (Paper 14)'),
        'provenance': {
            'source': 'benchmarks/qubit_scaling_data.json',
            'generator': 'benchmarks/qubit_scaling_sweep.py',
            'timestamp': data.get('timestamp'),
        },
        'points': rows,
    }
    he_pts = [(r['Q'], r['n_pauli']) for r in rows
              if r['type'] == 'geovac' and r['label'].startswith('GeoVac He')]
    if he_pts:
        block['fit_he'] = _loglog_fit(he_pts)
    return block


def export_qc_library(
    systems: Optional[Sequence[str]] = None,
    details: Optional[Dict[str, Dict[str, Any]]] = None,
    sweep_max_n_values: Sequence[int] = (1, 2, 3, 4),
) -> Dict[str, Any]:
    """Export the library summary table + the scaling datasets.

    The scaling families are deliberately distinct and labeled
    (do not conflate — they back different claims):

    - ``cross_system``: (Q, N_Pauli) across the library at the default
      build. Expected near-linear (the composed N_Pauli ~ 11.10 x Q
      relationship, CLAUDE.md §1.5).
    - ``within_system``: composed LiH as the basis grows, computed by
      ``geovac.composed_qubit.composed_lih_scaling_sweep`` — the SAME
      code path as the claim-backing test
      ``tests/test_paper14_scaling.py::test_composed_pauli_scaling_exponent``
      (composed O(Q^2.5), Paper 14). Empty ``sweep_max_n_values`` skips it.
    - ``atomic_benchmark``: the committed He/H deep-sweep artifact
      (see ``_benchmark_scaling_block``).

    ``details`` lets ``export_all`` reuse already-built per-system exports
    instead of rebuilding all systems.
    """
    from geovac.ecosystem_export import _SYSTEM_REGISTRY

    if systems is None:
        systems = sorted(set(_SYSTEM_REGISTRY.values()))

    rows: List[Dict[str, Any]] = []
    for name in systems:
        if details is not None and name in details:
            d = details[name]
            row = {
                'system': name,
                'n_qubits': d['n_qubits'],
                'n_terms': d['n_terms'],
                'one_norm': d['one_norm'],
                'n_qwc_groups': d['measurement']['n_qwc_groups'],
                'tapering': [
                    {k: t[k] for k in ('mode', 'n_qubits', 'n_terms')
                     if k in t} if 'error' not in t
                    else {'mode': t['mode'], 'error': t['error']}
                    for t in d['tapering']
                ],
            }
        else:
            d = export_qc_system(name, tapered_modes=(),
                                 include_pauli_terms=False)
            row = {
                'system': name,
                'n_qubits': d['n_qubits'],
                'n_terms': d['n_terms'],
                'one_norm': d['one_norm'],
                'n_qwc_groups': d['measurement']['n_qwc_groups'],
                'tapering': [],
            }
        rows.append(row)

    cross_points = [(r['n_qubits'], r['n_terms']) for r in rows]
    cross_system = {
        'label': ('library cross-system at the default build; expected '
                  'near-linear (composed N_Pauli ~ 11.10 x Q, CLAUDE.md '
                  '§1.5) — NOT the O(Q^2.5) within-system family'),
        'points': [{'system': r['system'], 'Q': r['n_qubits'],
                    'n_pauli': r['n_terms']} for r in rows],
        'fit': _loglog_fit(cross_points),
    }

    scaling: Dict[str, Any] = {'cross_system': cross_system}

    if sweep_max_n_values:
        from geovac.composed_qubit import composed_lih_scaling_sweep
        sw = composed_lih_scaling_sweep(
            max_n_values=list(sweep_max_n_values), verbose=False)
        # wall_time_s deliberately dropped (determinism contract)
        sweep_points = [
            {k: e[k] for k in ('max_n', 'M', 'Q', 'N_pauli',
                               'N_pauli_h1', 'N_pauli_eri') if k in e}
            for e in sw['sweep_data']
        ]
        scaling['within_system'] = {
            'label': ('composed LiH within-molecule basis sweep '
                      '(composed_lih_scaling_sweep, the code path backing '
                      'tests/test_paper14_scaling.py::'
                      'test_composed_pauli_scaling_exponent); the family '
                      'behind the composed O(Q^2.5) claim (Paper 14)'),
            'system': 'LiH',
            'points': sweep_points,
            'fit': _loglog_fit([(p['Q'], p['N_pauli'])
                                for p in sweep_points]),
        }

    bench = _benchmark_scaling_block()
    if bench is not None:
        scaling['atomic_benchmark'] = bench

    out = _envelope('qc_library')
    out.update({
        'systems': rows,
        'n_systems': len(rows),
        'scaling': scaling,
    })
    return out


# ---------------------------------------------------------------------------
# Top-level driver
# ---------------------------------------------------------------------------

def export_all(
    out_dir: str,
    lattice_zs: Sequence[int] = (1,),
    lattice_max_ns: Sequence[int] = (2, 3, 4, 5, 6, 7, 8, 9, 10),
    qc_systems: Optional[Sequence[str]] = None,
    tapered_modes: Sequence[str] = DEFAULT_TAPERED_MODES,
    sweep_max_n_values: Sequence[int] = (1, 2, 3, 4),
    verbose: bool = True,
) -> Dict[str, Any]:
    """Generate the full Phase-0 artifact set under ``out_dir``.

    Returns the manifest (also written to ``out_dir/manifest.json``).
    """
    from geovac.ecosystem_export import _SYSTEM_REGISTRY

    if qc_systems is None:
        qc_systems = sorted(set(_SYSTEM_REGISTRY.values()))

    files: List[str] = []

    def _emit(obj: Dict[str, Any], rel: str) -> None:
        path = os.path.join(out_dir, rel)
        _write_json(obj, path)
        files.append(rel.replace(os.sep, '/'))
        if verbose:
            print(f'  wrote {rel}')

    for Z in lattice_zs:
        for mn in lattice_max_ns:
            _emit(export_lattice(Z=Z, max_n=mn),
                  os.path.join('lattice', f'lattice_Z{Z}_n{mn}.json'))
            _emit(export_spectrum(Z=Z, max_n=mn),
                  os.path.join('spectrum', f'spectrum_Z{Z}_n{mn}.json'))

    details: Dict[str, Dict[str, Any]] = {}
    for name in qc_systems:
        if verbose:
            print(f'building {name} ...')
        details[name] = export_qc_system(name, tapered_modes=tapered_modes)
        _emit(details[name], os.path.join('qc', 'systems', f'{name}.json'))

    _emit(export_qc_library(systems=qc_systems, details=details,
                            sweep_max_n_values=sweep_max_n_values),
          os.path.join('qc', 'library.json'))

    manifest = _envelope('manifest')
    manifest['files'] = sorted(files)
    _write_json(manifest, os.path.join(out_dir, 'manifest.json'))
    return manifest


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description='Export GeoVac visualization data artifacts (Phase 0, '
                    'docs/visualization_plan.md).')
    parser.add_argument('--out', default=os.path.join('viz', 'public', 'data'),
                        help='output directory (default: viz/public/data)')
    parser.add_argument('--quick', action='store_true',
                        help='smoke subset: H lattice n<=4, LiH/He only, '
                             'per_block tapering, no max_n sweep')
    args = parser.parse_args(argv)

    if args.quick:
        export_all(args.out,
                   lattice_max_ns=(2, 3, 4),
                   qc_systems=('LiH', 'He'),
                   tapered_modes=('per_block',),
                   sweep_max_n_values=())
    else:
        export_all(args.out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
