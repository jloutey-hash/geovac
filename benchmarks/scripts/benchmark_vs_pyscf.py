"""
Benchmark: GeoVac vs PySCF External Comparison
===============================================

Paper 6's benchmarks are currently self-referential. This script provides
an external comparison against PySCF (a standard quantum chemistry package)
to compare GeoVac's legacy graph models with a standard package.

Systems compared:
  - H atom (single electron)
  - He atom (two electrons, HF level)
  - H2 molecule (bond length R = 1.4 bohr)

For each system, we report:
  - Energy (Ha)
  - Computation time (ms)
  - % error vs experiment/exact values

Corrected 2026-09-28: graph-mode hydrogen figures are band-edge
saturation, not accuracy; no GeoVac win is claimed.

Author: Precision Audit Suite
Date: February 2026

Usage:
    pip install pyscf  # if needed
    python benchmarks/scripts/benchmark_vs_pyscf.py
"""

import numpy as np
import time
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))


# =============================================================================
# Exact/Experimental Reference Values
# =============================================================================

REFERENCE = {
    'H': {
        'energy_exact': -0.5,          # Ha, exact non-relativistic
        'description': 'Hydrogen atom ground state',
    },
    'He': {
        'energy_exact': -2.9037,       # Ha, exact non-relativistic (Hylleraas)
        'energy_hf_limit': -2.8617,    # Ha, HF limit
        'description': 'Helium atom ground state',
    },
    'H2': {
        'energy_exact': -1.1745,       # Ha, exact non-relativistic at R=1.4 bohr
        'energy_hf_limit': -1.1336,    # Ha, HF/CBS limit at R=1.4 bohr
        'bond_length_bohr': 1.4,
        'description': 'H2 molecule at R=1.4 bohr',
    },
}


# =============================================================================
# GeoVac Calculations
# =============================================================================

def run_geovac_hydrogen(max_n: int = 20) -> dict:
    """Run hydrogen with the LEGACY GeoVac graph model H = kappa*(D - A).

    This is NOT a computed comparison of hydrogen physics.  The physical
    default (``hamiltonian='exact'``) inserts E_n = -1/(2n^2) on the labels
    by construction, so its "0% error" would be an inserted value, not a
    result.  We therefore run the labelled legacy graph model, whose ground
    energy kappa*lambda_max is matched to -1/2 by construction through
    kappa = -1/16; the reported "error" is only the band-edge saturation
    deficit (8 - lambda_max)/8 of the graph Laplacian (Paper 7), not an
    accuracy figure for hydrogen.
    """
    from geovac import AtomicSolver

    t0 = time.perf_counter()
    solver = AtomicSolver(max_n=max_n, Z=1, kinetic_scale=-1/16,
                          hamiltonian='graph')
    energies, _ = solver.compute_ground_state(n_states=1)
    t1 = time.perf_counter()

    return {
        'energy': energies[0],
        'time_ms': (t1 - t0) * 1000,
        'n_basis': solver.n_states,
        'method': f'GeoVac legacy graph (max_n={max_n}; saturation, not physics)',
    }


def run_geovac_helium(max_n: int = 5) -> dict:
    """Run helium with GeoVac using MoleculeHamiltonian + variational Z_eff optimization.

    Matches the production_suite.py path exactly (max_n=5, Z_eff ≈ 1.7).
    HeliumHamiltonian is NOT used here — it has a wrong kinetic sign and classical V_ee.
    """
    from geovac import MoleculeHamiltonian

    t0 = time.perf_counter()
    mol = MoleculeHamiltonian(
        nuclei=[(0.0, 0.0, 0.0)],
        nuclear_charges=[2],
        max_n=5,
    )
    result = mol.optimize_effective_charge(method='full_ci', n_points=12, z_range=(0.7, 1.0))
    mol.set_effective_charges(result['z_eff_optimal'])
    energies, _ = mol.compute_ground_state(n_states=1, method='full_ci')
    t1 = time.perf_counter()

    return {
        'energy': float(energies[0]),
        'time_ms': (t1 - t0) * 1000,
        'n_basis': mol.n_total_states,
        'method': 'GeoVac He (var. Z_eff + FCI, max_n=5)',
    }


def run_geovac_h2_mean_field(max_n: int = 10) -> dict:
    """
    Run H2 mean-field with GeoVac.

    Returns the lowest eigenvalue of the molecular graph Hamiltonian H(R).
    This is a single-particle energy (1 electron in the bonding MO field),
    NOT the total 2-electron molecular energy. Included for completeness and
    assembly-time comparison, not for energy accuracy comparison.
    """
    from geovac import GeometricLattice, MoleculeHamiltonian

    R = 1.4  # bohr

    t0 = time.perf_counter()
    atom_A = GeometricLattice(max_n=max_n, nucleus_position=(-R/2, 0, 0), nuclear_charge=1)
    atom_B = GeometricLattice(max_n=max_n, nucleus_position=(R/2, 0, 0), nuclear_charge=1)

    mol = MoleculeHamiltonian(
        lattices=[atom_A, atom_B],
        connectivity=[(0, 1, 3)],
        kinetic_scale=-1/16,
    )

    energies, _ = mol.compute_ground_state(method='mean_field')
    t1 = time.perf_counter()

    e_val = float(energies[0]) if hasattr(energies, '__len__') else float(energies)

    return {
        'energy': e_val,
        'time_ms': (t1 - t0) * 1000,
        'n_basis': atom_A.num_states + atom_B.num_states,
        'method': f'GeoVac H2 MF (max_n={max_n})',
        'note': 'single-particle eigenvalue — not comparable to total 2e energy',
    }


def run_geovac_h2_full_ci(max_n: int = 5) -> dict:
    """
    Run H2 Full CI with GeoVac at R=1.4 bohr.

    Uses the corrected CI Hamiltonian:
      H_total = H_kin⊗I + I⊗H_kin + V_cross_n2⊗I + I⊗V_cross_n1 + V_ee + V_NN

    where H_kin = kinetic_scale*(D-A) (pure-kinetic, no node weights) and
    V_cross uses the Mulliken approximation with a ⟨1/r⟩ cap to prevent
    variational collapse of diffuse high-n states.

    Legacy LCAO graph-concatenation model (superseded; FCI-M: no equilibrium
    geometry). ~2% error at max_n=5, ~1.7% at max_n=8, aided by error cancellation.
    Note: max_n=8 takes ~2 min (166k×166k matrix); default is max_n=5 (~30s).
    """
    from geovac import MoleculeHamiltonian

    R = 1.4  # bohr

    t0 = time.perf_counter()
    mol = MoleculeHamiltonian(
        nuclei=[(-R/2, 0.0, 0.0), (R/2, 0.0, 0.0)],
        nuclear_charges=[1.0, 1.0],
        max_n=max_n,
        connectivity=[(0, 1, 3)],
    )
    energies, _ = mol.compute_ground_state(n_states=1, method='full_ci')
    t1 = time.perf_counter()

    return {
        'energy': float(energies[0]),
        'time_ms': (t1 - t0) * 1000,
        'n_basis': mol.n_total_states,
        'method': f'GeoVac H2 FCI (max_n={max_n})',
    }


# keep old name as alias so any external callers don't break
def run_geovac_h2(max_n: int = 10) -> dict:
    return run_geovac_h2_mean_field(max_n=max_n)


# =============================================================================
# PySCF Calculations
# =============================================================================

def check_pyscf_available() -> bool:
    """Check if PySCF is installed."""
    try:
        import pyscf
        return True
    except ImportError:
        return False


def run_pyscf_hydrogen() -> dict:
    """Run hydrogen with PySCF HF/STO-3G."""
    from pyscf import gto, scf

    t0 = time.perf_counter()
    mol = gto.M(atom='H 0 0 0', basis='sto-3g', spin=1, charge=0, verbose=0)
    mf = scf.UHF(mol)
    energy = mf.kernel()
    t1 = time.perf_counter()

    return {
        'energy': energy,
        'time_ms': (t1 - t0) * 1000,
        'n_basis': mol.nao,
        'method': 'PySCF UHF/STO-3G',
    }


def run_pyscf_hydrogen_large() -> dict:
    """Run hydrogen with PySCF HF/cc-pVQZ for comparison."""
    from pyscf import gto, scf

    t0 = time.perf_counter()
    mol = gto.M(atom='H 0 0 0', basis='cc-pvqz', spin=1, charge=0, verbose=0)
    mf = scf.UHF(mol)
    energy = mf.kernel()
    t1 = time.perf_counter()

    return {
        'energy': energy,
        'time_ms': (t1 - t0) * 1000,
        'n_basis': mol.nao,
        'method': 'PySCF UHF/cc-pVQZ',
    }


def run_pyscf_helium() -> dict:
    """Run helium with PySCF RHF/STO-3G."""
    from pyscf import gto, scf

    t0 = time.perf_counter()
    mol = gto.M(atom='He 0 0 0', basis='sto-3g', spin=0, charge=0, verbose=0)
    mf = scf.RHF(mol)
    energy = mf.kernel()
    t1 = time.perf_counter()

    return {
        'energy': energy,
        'time_ms': (t1 - t0) * 1000,
        'n_basis': mol.nao,
        'method': 'PySCF RHF/STO-3G',
    }


def run_pyscf_helium_large() -> dict:
    """Run helium with PySCF RHF/cc-pVQZ."""
    from pyscf import gto, scf

    t0 = time.perf_counter()
    mol = gto.M(atom='He 0 0 0', basis='cc-pvqz', spin=0, charge=0, verbose=0)
    mf = scf.RHF(mol)
    energy = mf.kernel()
    t1 = time.perf_counter()

    return {
        'energy': energy,
        'time_ms': (t1 - t0) * 1000,
        'n_basis': mol.nao,
        'method': 'PySCF RHF/cc-pVQZ',
    }


def run_pyscf_h2() -> dict:
    """Run H2 with PySCF RHF/STO-3G at R=1.4 bohr."""
    from pyscf import gto, scf

    # PySCF uses Angstrom by default; 1.4 bohr = 0.74086 Angstrom
    R_angstrom = 1.4 * 0.529177

    t0 = time.perf_counter()
    mol = gto.M(
        atom=f'H 0 0 0; H 0 0 {R_angstrom}',
        basis='sto-3g',
        verbose=0,
    )
    mf = scf.RHF(mol)
    energy = mf.kernel()
    t1 = time.perf_counter()

    return {
        'energy': energy,
        'time_ms': (t1 - t0) * 1000,
        'n_basis': mol.nao,
        'method': 'PySCF RHF/STO-3G',
    }


def run_pyscf_h2_large() -> dict:
    """Run H2 with PySCF RHF/cc-pVQZ at R=1.4 bohr."""
    from pyscf import gto, scf

    R_angstrom = 1.4 * 0.529177

    t0 = time.perf_counter()
    mol = gto.M(
        atom=f'H 0 0 0; H 0 0 {R_angstrom}',
        basis='cc-pvqz',
        verbose=0,
    )
    mf = scf.RHF(mol)
    energy = mf.kernel()
    t1 = time.perf_counter()

    return {
        'energy': energy,
        'time_ms': (t1 - t0) * 1000,
        'n_basis': mol.nao,
        'method': 'PySCF RHF/cc-pVQZ',
    }


def run_pyscf_h2_fci() -> dict:
    """Run H2 Full CI with PySCF FCI/STO-3G at R=1.4 bohr.

    STO-3G FCI is the apples-to-apples comparison for GeoVac Full CI:
    both use a minimal basis. Expected: ~-1.1371 Ha (~3.2% error vs exact).
    """
    from pyscf import gto, scf, fci

    R_angstrom = 1.4 * 0.529177

    t0 = time.perf_counter()
    mol = gto.M(
        atom=f'H 0 0 0; H 0 0 {R_angstrom}',
        basis='sto-3g',
        verbose=0,
    )
    mf = scf.RHF(mol)
    mf.kernel()
    cisolver = fci.FCI(mf)
    energy, _ = cisolver.kernel()
    t1 = time.perf_counter()

    return {
        'energy': energy,
        'time_ms': (t1 - t0) * 1000,
        'n_basis': mol.nao,
        'method': 'PySCF FCI/STO-3G',
    }


def run_pyscf_h2_fci_large() -> dict:
    """Run H2 Full CI with PySCF FCI/cc-pVQZ at R=1.4 bohr.

    Near-basis-set-limit FCI. Expected: ~-1.1726 Ha (~0.16% error vs exact).
    """
    from pyscf import gto, scf, fci

    R_angstrom = 1.4 * 0.529177

    t0 = time.perf_counter()
    mol = gto.M(
        atom=f'H 0 0 0; H 0 0 {R_angstrom}',
        basis='cc-pvqz',
        verbose=0,
    )
    mf = scf.RHF(mol)
    mf.kernel()
    cisolver = fci.FCI(mf)
    energy, _ = cisolver.kernel()
    t1 = time.perf_counter()

    return {
        'energy': energy,
        'time_ms': (t1 - t0) * 1000,
        'n_basis': mol.nao,
        'method': 'PySCF FCI/cc-pVQZ',
    }


# =============================================================================
# Report Generation
# =============================================================================

def generate_comparison_report() -> str:
    """Generate the full comparison report."""
    lines = []
    lines.append("=" * 78)
    lines.append("GEOVAC vs PYSCF: EXTERNAL BENCHMARK COMPARISON")
    lines.append("=" * 78)
    lines.append("")

    has_pyscf = check_pyscf_available()

    if not has_pyscf:
        lines.append("WARNING: PySCF is not installed.")
        lines.append("Install with: pip install pyscf")
        lines.append("Running GeoVac benchmarks only.")
        lines.append("")

    # -------------------------------------------------------------------------
    # Hydrogen
    # -------------------------------------------------------------------------
    lines.append("SYSTEM 1: Hydrogen Atom (H)")
    lines.append("-" * 78)
    lines.append(f"  Exact energy: {REFERENCE['H']['energy_exact']:.6f} Ha")
    lines.append("")

    h_results = []

    # GeoVac at multiple resolutions
    for max_n in [10, 20, 30]:
        r = run_geovac_hydrogen(max_n=max_n)
        error = abs((r['energy'] - REFERENCE['H']['energy_exact']) /
                    REFERENCE['H']['energy_exact']) * 100
        r['error_pct'] = error
        r['system'] = 'H'
        h_results.append(r)

    if has_pyscf:
        for runner in [run_pyscf_hydrogen, run_pyscf_hydrogen_large]:
            try:
                r = runner()
                error = abs((r['energy'] - REFERENCE['H']['energy_exact']) /
                            REFERENCE['H']['energy_exact']) * 100
                r['error_pct'] = error
                r['system'] = 'H'
                h_results.append(r)
            except Exception as e:
                lines.append(f"  PySCF error: {e}")

    lines.append("  Sat.deficit% = |E - E_exact|/|E_exact|. For the GeoVac graph-mode")
    lines.append("  row this is the band-edge saturation deficit (E_0 = kappa*lambda_max")
    lines.append("  by construction), not an accuracy.")
    lines.append(f"  {'Method':<30s} {'Energy(Ha)':>12s} {'Sat.deficit%':>12s} "
                 f"{'Time(ms)':>10s} {'N_basis':>8s}")
    lines.append(f"  {'-'*30} {'-'*12} {'-'*12} {'-'*10} {'-'*8}")
    for r in h_results:
        lines.append(f"  {r['method']:<30s} {r['energy']:>12.6f} "
                     f"{r['error_pct']:>12.4f} {r['time_ms']:>10.2f} "
                     f"{str(r['n_basis']):>8s}")
    lines.append("")

    # -------------------------------------------------------------------------
    # Helium
    # -------------------------------------------------------------------------
    lines.append("SYSTEM 2: Helium Atom (He)")
    lines.append("-" * 78)
    lines.append(f"  Exact energy: {REFERENCE['He']['energy_exact']:.6f} Ha")
    lines.append(f"  HF limit:    {REFERENCE['He']['energy_hf_limit']:.6f} Ha")
    lines.append("")

    he_results = []

    try:
        r = run_geovac_helium()
        r['error_pct'] = abs((r['energy'] - REFERENCE['He']['energy_exact']) /
                             REFERENCE['He']['energy_exact']) * 100
        r['system'] = 'He'
        he_results.append(r)
    except Exception as e:
        lines.append(f"  GeoVac Helium error: {e}")

    if has_pyscf:
        for runner in [run_pyscf_helium, run_pyscf_helium_large]:
            try:
                r = runner()
                error = abs((r['energy'] - REFERENCE['He']['energy_exact']) /
                            REFERENCE['He']['energy_exact']) * 100
                r['error_pct'] = error
                r['system'] = 'He'
                he_results.append(r)
            except Exception as e:
                lines.append(f"  PySCF error: {e}")

    if he_results:
        lines.append(f"  {'Method':<30s} {'Energy(Ha)':>12s} {'Error%':>8s} "
                     f"{'Time(ms)':>10s} {'N_basis':>8s}")
        lines.append(f"  {'-'*30} {'-'*12} {'-'*8} {'-'*10} {'-'*8}")
        for r in he_results:
            lines.append(f"  {r['method']:<30s} {r['energy']:>12.6f} "
                         f"{r['error_pct']:>8.4f} {r['time_ms']:>10.2f} "
                         f"{str(r['n_basis']):>8s}")
    lines.append("")

    # -------------------------------------------------------------------------
    # H2 Molecule
    # -------------------------------------------------------------------------
    lines.append("SYSTEM 3: H2 Molecule (R = 1.4 bohr)")
    lines.append("-" * 78)
    lines.append(f"  Exact energy: {REFERENCE['H2']['energy_exact']:.6f} Ha")
    lines.append(f"  HF limit:    {REFERENCE['H2']['energy_hf_limit']:.6f} Ha")
    lines.append("")

    h2_results = []

    # GeoVac Full CI (total 2-electron energy — primary result)
    try:
        for max_n in [5, 8]:
            r = run_geovac_h2_full_ci(max_n=max_n)
            r['error_pct'] = abs((r['energy'] - REFERENCE['H2']['energy_exact']) /
                                 REFERENCE['H2']['energy_exact']) * 100
            r['system'] = 'H2'
            h2_results.append(r)
    except Exception as e:
        lines.append(f"  GeoVac H2 Full CI error: {e}")

    # GeoVac mean-field (single-particle, assembly-time reference only)
    lines.append("  [MF = single-particle eigenvalue, not total molecular energy]")
    try:
        r = run_geovac_h2_mean_field(max_n=5)
        r['error_pct'] = abs((r['energy'] - REFERENCE['H2']['energy_exact']) /
                             REFERENCE['H2']['energy_exact']) * 100
        r['system'] = 'H2'
        r['note'] = 'single-particle eigenvalue — not comparable to 2e total'
        h2_results.append(r)
    except Exception as e:
        lines.append(f"  GeoVac H2 mean-field error: {e}")

    if has_pyscf:
        for runner in [run_pyscf_h2, run_pyscf_h2_fci,
                       run_pyscf_h2_large, run_pyscf_h2_fci_large]:
            try:
                r = runner()
                error = abs((r['energy'] - REFERENCE['H2']['energy_exact']) /
                            REFERENCE['H2']['energy_exact']) * 100
                r['error_pct'] = error
                r['system'] = 'H2'
                h2_results.append(r)
            except Exception as e:
                lines.append(f"  PySCF H2 error: {e}")

    if h2_results:
        lines.append(f"  {'Method':<34s} {'Energy(Ha)':>12s} {'vs Exact%':>10s} "
                     f"{'Time(ms)':>10s} {'N_basis':>8s}")
        lines.append(f"  {'-'*34} {'-'*12} {'-'*10} {'-'*10} {'-'*8}")
        for r in h2_results:
            note = r.get('note', '')
            lines.append(f"  {r['method']:<34s} {r['energy']:>12.6f} "
                         f"{r['error_pct']:>10.4f} {r['time_ms']:>10.2f} "
                         f"{str(r['n_basis']):>8s}")
            if note:
                lines.append(f"    ^ {note}")
    lines.append("")

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    lines.append("=" * 78)
    lines.append("SUMMARY (corrected 2026-09-28: retracted claims removed)")
    lines.append("=" * 78)
    lines.append("")
    lines.append("WHAT THE TABLES SHOW:")
    lines.append("  - H (graph mode): the legacy graph model H = kappa*(D - A) has")
    lines.append("    E_0 = kappa*lambda_max matched to -1/2 by construction")
    lines.append("    (kappa = -1/16). Its 'Sat.deficit%' column is the band-edge")
    lines.append("    saturation deficit (8 - lambda_max)/8, NOT an accuracy.")
    lines.append("  - The 18 symbolic checks (tests/test_fock_*.py) verify the")
    lines.append("    CONTINUUM Fock projection. The lattice's edge Laplacian")
    lines.append("    D - A does NOT converge to the S3 Laplace-Beltrami operator")
    lines.append("    (n_max disconnected l-blocks; a documented wall, Paper 7).")
    lines.append("  - Speed/sparsity figures are for the graph eigenproblem only;")
    lines.append("    they are not an accuracy-matched comparison with PySCF.")
    lines.append("")
    lines.append("LIMITATIONS:")
    lines.append("  - H2 'Full CI' here is the legacy LCAO graph-concatenation")
    lines.append("    model, superseded (FCI-M: no equilibrium geometry). Its")
    lines.append("    ~2% figure benefits from error cancellation between the")
    lines.append("    Mulliken V_cross and point-charge V_ee approximations.")
    lines.append("  - Multi-electron: Mean-field He uses Z-scaled kinetic_scale;")
    lines.append("    no general DFT/HF multi-electron solver in this script.")
    lines.append("  - No basis-set extrapolation procedure established.")
    lines.append("")
    lines.append("VERDICT:")
    lines.append("  This script supports no GeoVac accuracy claim. The graph-mode")
    lines.append("  hydrogen energies are band-edge saturation, the H2 model is")
    lines.append("  superseded, and the comparison is not at matched accuracy.")
    lines.append("  Current GeoVac results live in the papers (e.g. Papers 12,")
    lines.append("  13, 14), benchmarked against strong baselines there.")
    lines.append("")

    return "\n".join(lines)


if __name__ == "__main__":
    report = generate_comparison_report()
    print(report)

    # Save comparison table
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'figures')
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'pyscf_comparison.txt')
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(report)
    print(f"\nReport saved to: {output_path}")
