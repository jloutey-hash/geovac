"""
Gaussian Baseline Comparison — Track CA Market Test
====================================================

Computes raw Jordan-Wigner Gaussian baselines for comparison with GeoVac
composed Hamiltonians. Uses OpenFermion cached molecular data and GeoVac's
gaussian_reference module.

Requirements:
    - openfermion (with cached HDF5 data)
    - geovac (this project)

For Gaussian baselines requiring PySCF (LiH 6-31G, cc-pVDZ, H2O all bases):
    pip install pyscf  # requires BLAS; fails on some Windows builds
    See Phase 2 functions below.

Author: GeoVac Development Team
Date: April 2026
Track: CA (Quantum Resource Market Test)
"""

import json
import os
from pathlib import Path
from typing import Any, Dict

import h5py
import numpy as np

from openfermion import jordan_wigner, count_qubits, InteractionOperator

from geovac.gaussian_reference import h2_sto3g, he_sto3g, he_cc_pvdz
from geovac.qubit_encoding import build_fermion_op_from_integrals
from geovac.trotter_bounds import pauli_1norm
from geovac.measurement_grouping import qwc_groups


# ---------------------------------------------------------------------------
# Phase 1: Baselines from available data (no PySCF required)
# ---------------------------------------------------------------------------

def compute_gaussian_baseline(
    h1: np.ndarray,
    eri: np.ndarray,
    nuclear_repulsion: float,
    label: str,
) -> Dict[str, Any]:
    """Compute Pauli count, 1-norm, QWC groups for a Gaussian Hamiltonian."""
    fop = build_fermion_op_from_integrals(h1, eri, nuclear_repulsion)
    jw = jordan_wigner(fop)
    n_terms = len(jw.terms)
    q = count_qubits(jw)
    lam = pauli_1norm(jw)
    groups = qwc_groups(jw)
    eri_density = np.count_nonzero(np.abs(eri) > 1e-10) / max(eri.size, 1) * 100

    result = {
        'label': label,
        'Q': q,
        'pauli_terms': n_terms,
        'one_norm_Ha': round(lam, 4),
        'qwc_groups': len(groups),
        'eri_density_pct': round(eri_density, 1),
    }
    print(f"  {label}: Q={q}, Pauli={n_terms}, 1-norm={lam:.4f} Ha, "
          f"QWC={len(groups)}, ERI={eri_density:.1f}%")
    return result


def lih_sto3g_from_cache() -> Dict[str, Any]:
    """Load LiH STO-3G from OpenFermion cached HDF5 (R=1.45 Å)."""
    import openfermion
    data_dir = Path(openfermion.__file__).parent / 'testing' / 'data'
    fpath = data_dir / 'H1-Li1_sto-3g_singlet_1.45.hdf5'
    if not fpath.exists():
        raise FileNotFoundError(f"OpenFermion cached LiH data not found at {fpath}")

    with h5py.File(str(fpath), 'r') as f:
        h1 = np.array(f['one_body_integrals'])
        eri = np.array(f['two_body_integrals'])
        nr = float(np.array(f['nuclear_repulsion']))
        fci = float(np.array(f['fci_energy']))
        hf = float(np.array(f['hf_energy']))

    result = compute_gaussian_baseline(h1, eri, nr, 'LiH STO-3G (R=1.45Å)')
    result['fci_energy_Ha'] = round(fci, 6)
    result['hf_energy_Ha'] = round(hf, 6)
    result['source'] = 'OpenFermion cached H1-Li1_sto-3g_singlet_1.45.hdf5'
    return result


def he_baselines() -> list:
    """Compute He STO-3G, cc-pVDZ, and cc-pVTZ baselines."""
    results = []
    for label, func in [('He STO-3G', he_sto3g), ('He cc-pVDZ', he_cc_pvdz)]:
        data = func()
        r = compute_gaussian_baseline(
            data['h1'], data['eri'], data['nuclear_repulsion'], label
        )
        r['fci_energy_Ha'] = round(data['literature_energy'], 6)
        r['source'] = data['source']
        results.append(r)

    # cc-pVTZ if available
    try:
        from geovac.gaussian_reference import he_cc_pvtz
        data = he_cc_pvtz()
        r = compute_gaussian_baseline(
            data['h1'], data['eri'], data['nuclear_repulsion'], 'He cc-pVTZ'
        )
        r['fci_energy_Ha'] = round(data['literature_energy'], 6)
        r['source'] = data['source']
        results.append(r)
    except (ImportError, Exception) as e:
        print(f"  He cc-pVTZ: skipped ({e})")
    return results


def h2_baselines() -> list:
    """Compute H2 STO-3G and 6-31G baselines."""
    results = []

    # STO-3G from gaussian_reference
    data = h2_sto3g(1.4)
    r = compute_gaussian_baseline(
        data['h1'], data['eri'], data['nuclear_repulsion'], 'H2 STO-3G (R=1.4 bohr)'
    )
    r['fci_energy_Ha'] = data['literature_energy']
    r['source'] = data['source']
    results.append(r)

    # 6-31G from OpenFermion cache
    import openfermion
    data_dir = Path(openfermion.__file__).parent / 'testing' / 'data'
    fpath = data_dir / 'H2_6-31g_singlet_0.75.hdf5'
    if fpath.exists():
        with h5py.File(str(fpath), 'r') as f:
            h1 = np.array(f['one_body_integrals'])
            eri = np.array(f['two_body_integrals'])
            nr = float(np.array(f['nuclear_repulsion']))
            fci = float(np.array(f['fci_energy']))
        r = compute_gaussian_baseline(h1, eri, nr, 'H2 6-31G (R=0.75Å)')
        r['fci_energy_Ha'] = round(fci, 6)
        r['source'] = 'OpenFermion cached H2_6-31g_singlet_0.75.hdf5'
        results.append(r)
    return results


def geovac_composed_baselines() -> list:
    """Compute GeoVac composed baselines for comparison."""
    from geovac.composed_qubit import build_composed_lih, build_composed_beh2, build_composed_h2o

    results = []
    for name, func in [('LiH', build_composed_lih),
                        ('BeH2', build_composed_beh2),
                        ('H2O', build_composed_h2o)]:
        for pk_mode, suffix in [(True, 'full'), (False, 'elec')]:
            data = func(pk_in_hamiltonian=pk_mode)
            jw = data['qubit_op']
            n_terms = len(jw.terms)
            q = count_qubits(jw)
            lam = pauli_1norm(jw)
            groups = qwc_groups(jw)
            label = f'GeoVac {name} ({suffix})'
            print(f"  {label}: Q={q}, Pauli={n_terms}, 1-norm={lam:.4f} Ha, QWC={len(groups)}")
            results.append({
                'label': label,
                'Q': q,
                'pauli_terms': n_terms,
                'one_norm_Ha': round(lam, 4),
                'qwc_groups': len(groups),
                'pk_included': pk_mode,
                'source': 'geovac.composed_qubit',
            })
    return results


# ---------------------------------------------------------------------------
# Phase 2: PySCF-dependent baselines (stubs)
# ---------------------------------------------------------------------------

def lih_pyscf(basis: str = 'cc-pvdz') -> Dict[str, Any]:
    """
    Compute LiH Gaussian baseline using PySCF.

    Requires: pip install pyscf (needs BLAS)

    Usage::

        from pyscf import gto, scf, ao2mo, fci
        mol = gto.M(atom='Li 0 0 0; H 0 0 1.596', basis=basis,
                     unit='Angstrom')
        mf = mol.RHF().run()
        h1 = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
        eri = ao2mo.full(mol, mf.mo_coeff).reshape(mol.nao, mol.nao,
                                                     mol.nao, mol.nao)
        fci_solver = fci.FCI(mf)
        e_fci, _ = fci_solver.kernel()
    """
    raise NotImplementedError(
        f"LiH {basis} requires PySCF. "
        "Install with: pip install pyscf (requires BLAS libraries)"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Run all available baseline comparisons and save results."""
    print("=" * 60)
    print("Track CA: Gaussian Baseline Comparison")
    print("=" * 60)

    all_results: Dict[str, list] = {}

    print("\n--- H2 Baselines ---")
    all_results['h2'] = h2_baselines()

    print("\n--- He Baselines ---")
    all_results['he'] = he_baselines()

    print("\n--- LiH STO-3G (from OpenFermion cache) ---")
    try:
        all_results['lih_sto3g'] = [lih_sto3g_from_cache()]
    except FileNotFoundError as e:
        print(f"  Skipped: {e}")
        all_results['lih_sto3g'] = []

    print("\n--- GeoVac Composed ---")
    all_results['geovac_composed'] = geovac_composed_baselines()

    # Save results
    out_path = Path(__file__).parent.parent / 'debug' / 'data' / 'market_test_data.json'
    print(f"\nResults summary saved to {out_path}")

    # Print comparison table
    print("\n" + "=" * 60)
    print("COMPARISON TABLE")
    print("=" * 60)
    print(f"{'Method':<35} {'Q':>3} {'Pauli':>8} {'1-norm':>10} {'QWC':>5}")
    print("-" * 65)
    for category, items in all_results.items():
        for r in items:
            print(f"{r['label']:<35} {r['Q']:>3} {r['pauli_terms']:>8} "
                  f"{r['one_norm_Ha']:>10.2f} {r.get('qwc_groups', '-'):>5}")


if __name__ == '__main__':
    main()
