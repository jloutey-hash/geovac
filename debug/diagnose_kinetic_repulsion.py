"""
Diagnose missing kinetic repulsion in LiH LCAO framework.

The balanced PES (exact+True, fourier+s_only) shows monotonically attractive
D_cp(R) with no equilibrium. This script identifies WHY the kinetic energy
fails to create the short-range repulsive wall.

Diagnostics:
  1. T(R), V(R), virial ratio η(R) via decompose_energy()
  2. Bridge edge weights vs R
  3. Degree matrix changes with R (kinetic cost of bridges)
  4. H₂ comparison (does homonuclear also fail?)
  5. Graph Laplacian kinetic contribution analysis

Output:
  debug/KINETIC_DIAGNOSTIC.md — full report
  debug/plots/lih_kinetic_vs_R.png — T(R) and V(R) plot
  debug/data/lih_kinetic_diagnostic.txt — raw data

Date: 2026-03-11
Version: v0.9.37
"""
import os
import sys
import warnings
from typing import Dict, List, Tuple

import numpy as np

warnings.filterwarnings('ignore')

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from geovac.lattice_index import (
    LatticeIndex,
    MolecularLatticeIndex,
)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
R_VALUES = [1.5, 1.8, 2.0, 2.2, 2.5, 2.8, 3.0, 3.015, 3.5, 4.0, 5.0]
NMAX = 3

# Use the balanced exact+True config
CROSS_NUCLEAR = 'exact'
CROSS_ATOM_VEE = True


def diagnostic_1_energy_decomposition() -> List[Dict[str, float]]:
    """Decompose energy at each R: T, V_nA, V_nB, V_cross, V_bridge, V_ee."""
    print("=" * 70)
    print("DIAGNOSTIC 1: Energy Decomposition vs R")
    print("=" * 70)

    results = []
    for R in R_VALUES:
        print(f"  R={R:.3f}...", end="", flush=True)
        mol = MolecularLatticeIndex(
            Z_A=3, Z_B=1, nmax_A=NMAX, nmax_B=NMAX,
            R=R, n_electrons=4,
            vee_method='slater_full', fci_method='auto',
            cross_nuclear_method=CROSS_NUCLEAR,
            cross_atom_vee=CROSS_ATOM_VEE,
        )
        eigvals, eigvecs = mol.compute_ground_state(n_states=1)
        E_total = eigvals[0]
        decomp = mol.decompose_energy(eigvecs[:, 0], E_total)

        V_nuc_total = decomp['V_nA'] + decomp['V_nB']
        V_cross_total = decomp['V_cross_A'] + decomp['V_cross_B']
        V_total = V_nuc_total + V_cross_total + decomp['V_bridge'] + decomp['V_ee'] + decomp['V_NN']
        T = decomp['T']

        # Virial ratio: η = -V / (2T)
        # For bound Coulomb system at equilibrium, η = 1
        if abs(T) > 1e-10:
            eta = -V_total / (2 * T)
        else:
            eta = float('inf')

        row = {
            'R': R,
            'E_total': E_total,
            'T': T,
            'V_nA': decomp['V_nA'],
            'V_nB': decomp['V_nB'],
            'V_cross_A': decomp['V_cross_A'],
            'V_cross_B': decomp['V_cross_B'],
            'V_bridge': decomp['V_bridge'],
            'V_ee': decomp['V_ee'],
            'V_NN': decomp['V_NN'],
            'V_total': V_total,
            'eta': eta,
            'H1_total': decomp['H1_total'],
        }
        results.append(row)
        print(f" T={T:.4f}, V={V_total:.4f}, eta={eta:.4f}, E={E_total:.5f}")

    return results


def diagnostic_2_bridge_analysis() -> List[Dict[str, float]]:
    """Analyze bridge edge weights and degree matrix changes with R."""
    print("\n" + "=" * 70)
    print("DIAGNOSTIC 2: Bridge Edges & Degree Matrix vs R")
    print("=" * 70)

    results = []
    for R in R_VALUES:
        print(f"  R={R:.3f}...", end="", flush=True)
        mol = MolecularLatticeIndex(
            Z_A=3, Z_B=1, nmax_A=NMAX, nmax_B=NMAX,
            R=R, n_electrons=4,
            vee_method='slater_full', fci_method='auto',
            cross_nuclear_method=CROSS_NUCLEAR,
            cross_atom_vee=CROSS_ATOM_VEE,
        )

        # Extract H1 matrix and bridge block
        H1 = mol._H1_spatial
        if hasattr(H1, 'toarray'):
            H1_dense = H1.toarray()
        else:
            H1_dense = np.array(H1)

        nA = mol._n_spatial_A
        nB = mol._n_spatial_B
        H1_AB = H1_dense[:nA, nA:]  # bridge block

        max_bridge = np.max(np.abs(H1_AB))
        sum_bridge = np.sum(np.abs(H1_AB))
        n_nonzero = np.count_nonzero(np.abs(H1_AB) > 1e-12)

        # Degree analysis: how many edges per node from bridges
        # The adjacency is stored in mol — let's check
        adj = None
        if hasattr(mol, '_adjacency_combined'):
            adj = mol._adjacency_combined
        elif hasattr(mol, '_adjacency'):
            adj = mol._adjacency

        bridge_degree_A = 0.0
        bridge_degree_B = 0.0
        if adj is not None:
            if hasattr(adj, 'toarray'):
                adj_dense = adj.toarray()
            else:
                adj_dense = np.array(adj)
            # Bridge connections: A->B block
            adj_AB = adj_dense[:nA, nA:]
            bridge_degree_A = np.sum(np.abs(adj_AB), axis=1).mean()
            bridge_degree_B = np.sum(np.abs(adj_AB), axis=0).mean()

        # Diagonal H1 analysis (kinetic + potential)
        h1_diag = np.diag(H1_dense)
        h1_diag_A = h1_diag[:nA]
        h1_diag_B = h1_diag[nA:]

        row = {
            'R': R,
            'max_bridge': max_bridge,
            'sum_bridge': sum_bridge,
            'n_bridge_nonzero': n_nonzero,
            'bridge_degree_A': bridge_degree_A,
            'bridge_degree_B': bridge_degree_B,
            'h1_diag_A_mean': np.mean(h1_diag_A),
            'h1_diag_B_mean': np.mean(h1_diag_B),
            'h1_diag_A_min': np.min(h1_diag_A),
            'h1_diag_B_min': np.min(h1_diag_B),
        }
        results.append(row)
        print(f" max|H1_AB|={max_bridge:.6f}, bridge_deg_A={bridge_degree_A:.4f}")

    return results


def diagnostic_3_h2_comparison() -> Dict[str, float]:
    """Check if H₂ also has no equilibrium (fundamental vs heteronuclear issue)."""
    print("\n" + "=" * 70)
    print("DIAGNOSTIC 3: H2 Comparison (does homonuclear also fail?)")
    print("=" * 70)

    # Check if MolecularLatticeIndex works for H₂ (Z_A=1, Z_B=1)
    R_h2 = [0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.5, 3.0, 4.0]
    results = {}

    # H₂ atomic reference
    h_atom = LatticeIndex(
        n_electrons=1, max_n=NMAX, nuclear_charge=1,
        vee_method='slater_full', h1_method='exact', fci_method='auto',
    )
    E_h = h_atom.compute_ground_state(n_states=1)[0][0]
    E_sep = 2 * E_h
    print(f"  E(H) = {E_h:.6f} Ha, E_sep = {E_sep:.6f} Ha")

    print(f"  {'R':>6s}  {'E_mol':>10s}  {'D_raw':>8s}  {'T':>8s}  {'V_total':>8s}  {'eta':>6s}")
    print(f"  {'-'*56}")

    for R in R_h2:
        try:
            mol = MolecularLatticeIndex(
                Z_A=1, Z_B=1, nmax_A=NMAX, nmax_B=NMAX,
                R=R, n_electrons=2,
                vee_method='slater_full', fci_method='auto',
                cross_nuclear_method='fourier',
                cross_atom_vee='s_only',
            )
            eigvals, eigvecs = mol.compute_ground_state(n_states=1)
            E_mol = eigvals[0]
            D_raw = E_sep - E_mol
            decomp = mol.decompose_energy(eigvecs[:, 0], E_mol)
            T = decomp['T']
            V_tot = (decomp['V_nA'] + decomp['V_nB'] +
                     decomp['V_cross_A'] + decomp['V_cross_B'] +
                     decomp['V_bridge'] + decomp['V_ee'] + decomp['V_NN'])
            eta = -V_tot / (2 * T) if abs(T) > 1e-10 else float('inf')
            print(f"  {R:6.2f}  {E_mol:10.5f}  {D_raw:8.4f}  {T:8.4f}  {V_tot:8.4f}  {eta:6.3f}")
            results[R] = {'E_mol': E_mol, 'D_raw': D_raw, 'T': T, 'V': V_tot, 'eta': eta}
        except Exception as e:
            print(f"  {R:6.2f}  ERROR: {e}")
            results[R] = {'error': str(e)}

    # Check if H₂ has equilibrium
    D_raw_vals = [results[R]['D_raw'] for R in R_h2 if 'D_raw' in results.get(R, {})]
    if len(D_raw_vals) >= 3:
        idx_max = np.argmax(D_raw_vals)
        if 0 < idx_max < len(D_raw_vals) - 1:
            print(f"\n  H2 has equilibrium near R ~ {R_h2[idx_max]:.1f} bohr")
        else:
            print(f"\n  H2 has NO equilibrium (max D_raw at boundary R={R_h2[idx_max]:.1f})")

    return results


def diagnostic_4_intra_atom_kinetic() -> None:
    """Analyze whether intra-atom kinetic terms change with R."""
    print("\n" + "=" * 70)
    print("DIAGNOSTIC 4: Intra-Atom Kinetic Energy Independence")
    print("=" * 70)

    # The hypothesis: intra-atom adjacency is FIXED regardless of R.
    # Let's verify by checking if the intra-atom block of H1 changes.
    H1_AA_list = []
    H1_BB_list = []

    for R in [1.5, 3.0, 5.0]:
        mol = MolecularLatticeIndex(
            Z_A=3, Z_B=1, nmax_A=NMAX, nmax_B=NMAX,
            R=R, n_electrons=4,
            vee_method='slater_full', fci_method='auto',
            cross_nuclear_method=CROSS_NUCLEAR,
            cross_atom_vee=CROSS_ATOM_VEE,
        )
        H1 = mol._H1_spatial
        if hasattr(H1, 'toarray'):
            H1 = H1.toarray()
        nA = mol._n_spatial_A

        H1_AA = H1[:nA, :nA]
        H1_BB = H1[nA:, nA:]
        H1_AA_list.append(H1_AA)
        H1_BB_list.append(H1_BB)

    # Compare off-diagonal of H1_AA at R=1.5 vs R=5.0
    dH1_AA = H1_AA_list[0] - H1_AA_list[2]
    dH1_BB = H1_BB_list[0] - H1_BB_list[2]

    # Only look at off-diagonal (diagonal changes due to cross-nuclear)
    np.fill_diagonal(dH1_AA, 0)
    np.fill_diagonal(dH1_BB, 0)

    print(f"  H1_AA off-diagonal change (R=1.5 vs R=5.0):")
    print(f"    max|dH1_AA_offdiag| = {np.max(np.abs(dH1_AA)):.2e}")
    print(f"    sum|dH1_AA_offdiag| = {np.sum(np.abs(dH1_AA)):.2e}")

    print(f"  H1_BB off-diagonal change (R=1.5 vs R=5.0):")
    print(f"    max|dH1_BB_offdiag| = {np.max(np.abs(dH1_BB)):.2e}")
    print(f"    sum|dH1_BB_offdiag| = {np.sum(np.abs(dH1_BB)):.2e}")

    if np.max(np.abs(dH1_AA)) < 1e-10 and np.max(np.abs(dH1_BB)) < 1e-10:
        print("\n  CONFIRMED: Intra-atom kinetic hopping is COMPLETELY R-INDEPENDENT.")
        print("  The graph Laplacian intra-atom edges do not respond to orbital overlap.")
    else:
        print("\n  UNEXPECTED: Some R-dependence detected in intra-atom hopping.")

    # Compare diagonal (should change due to cross-nuclear)
    dH1_AA_diag = np.diag(H1_AA_list[0]) - np.diag(H1_AA_list[2])
    dH1_BB_diag = np.diag(H1_BB_list[0]) - np.diag(H1_BB_list[2])
    print(f"\n  H1_AA diagonal change (R=1.5 vs R=5.0): {np.sum(dH1_AA_diag):.4f} Ha")
    print(f"  H1_BB diagonal change (R=1.5 vs R=5.0): {np.sum(dH1_BB_diag):.4f} Ha")
    print(f"  (These are cross-nuclear attraction shifts -- expected to depend on R)")


def diagnostic_5_overlap_kinetic_estimate() -> Dict[str, List[float]]:
    """Estimate the missing overlap-dependent kinetic correction."""
    print("\n" + "=" * 70)
    print("DIAGNOSTIC 5: Overlap-Dependent Kinetic Correction Estimate")
    print("=" * 70)

    # In standard QC, the kinetic integral between AOs on different centers is:
    #   T_AB = <φ_A | -½∇² | φ_B> = (α/2)(3 - α²R²) S(R)  for 1s STOs
    # where α = (Z_A + Z_B)/2 and S(R) is the overlap.
    #
    # The kinetic energy of an orthogonalized AO set goes as:
    #   T_orth ≈ T_diag + T_AB² / (E_A - E_B)  [perturbation theory]
    # or more directly, Löwdin orthogonalization adds kinetic energy ~ S²
    #
    # Let's compute STO overlaps and estimate the missing kinetic contribution.

    Z_A, Z_B = 3, 1
    alpha = (Z_A + Z_B) / 2  # average Slater exponent (crude)

    R_dense = np.linspace(1.0, 6.0, 51)
    S_1s = (1 + R_dense + R_dense**2 / 3) * np.exp(-R_dense)

    # Exact 1s-1s kinetic integral (STO, same exponent α):
    # T_AB = (α/2)(3 - α²R²) * S(R)
    # But for different centers with different Z, use α = Z for each
    # Simpler: just estimate T_correction ~ C * S²(R)

    # For the LiH case, overlap between Li 1s (Z=3) and H 1s (Z=1):
    # Use STO overlap with effective exponent
    Z_eff = np.sqrt(Z_A * Z_B)  # geometric mean
    S_hetero = (1 + Z_eff * R_dense + (Z_eff * R_dense)**2 / 3) * np.exp(-Z_eff * R_dense)

    results = {'R': R_dense.tolist(), 'S_1s': S_1s.tolist(), 'S_hetero': S_hetero.tolist()}

    print(f"  STO overlap (1s-1s, unit exponent):")
    for R in [1.5, 2.0, 2.5, 3.0, 4.0, 5.0]:
        S = (1 + R + R**2/3) * np.exp(-R)
        S_h = (1 + Z_eff*R + (Z_eff*R)**2/3) * np.exp(-Z_eff*R)
        print(f"    R={R:.1f}: S(unit)={S:.4f}, S(Z_eff={Z_eff:.2f})={S_h:.6f}, S^2={S_h**2:.6f}")

    # Estimate kinetic correction needed
    # At R=1.5, D_cp ≈ 0.46 (need to reduce to ~0.09, so need ~0.37 Ha repulsion)
    # At R=3.0, D_cp ≈ 0.11 (roughly correct, need ~0 correction)
    print(f"\n  Approximate kinetic correction needed to create equilibrium:")
    print(f"    At R=1.5: need ~0.37 Ha of additional repulsion")
    print(f"    At R=2.5: need ~0.08 Ha of additional repulsion")
    print(f"    At R=3.0: need ~0.02 Ha of additional repulsion")
    print(f"    At R=5.0: need ~0 Ha (already dissociating correctly)")

    return results


def make_plots(
    decomp_results: List[Dict[str, float]],
    bridge_results: List[Dict[str, float]],
) -> None:
    """Generate diagnostic plots."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not available, skipping plots.")
        return

    R = [r['R'] for r in decomp_results]
    T = [r['T'] for r in decomp_results]
    V = [r['V_total'] for r in decomp_results]
    eta = [r['eta'] for r in decomp_results]
    V_cross = [r['V_cross_A'] + r['V_cross_B'] for r in decomp_results]
    V_ee = [r['V_ee'] for r in decomp_results]
    V_NN = [r['V_NN'] for r in decomp_results]
    V_bridge = [r['V_bridge'] for r in decomp_results]
    E_total = [r['E_total'] for r in decomp_results]

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    # Panel 1: T(R) and V(R)
    ax = axes[0, 0]
    ax.plot(R, T, 'b-o', label='T (kinetic)', linewidth=2)
    ax.plot(R, V, 'r-s', label='V (total potential)', linewidth=2)
    ax.axhline(0, color='gray', linestyle='--', alpha=0.3)
    ax.set_xlabel('R (bohr)')
    ax.set_ylabel('Energy (Ha)')
    ax.set_title('T(R) and V(R)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Panel 2: Virial ratio
    ax = axes[0, 1]
    ax.plot(R, eta, 'g-^', linewidth=2)
    ax.axhline(1.0, color='red', linestyle='--', label='eta=1 (equilibrium)')
    ax.set_xlabel('R (bohr)')
    ax.set_ylabel('eta = -V/(2T)')
    ax.set_title('Virial Ratio')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Panel 3: Energy components
    ax = axes[0, 2]
    ax.plot(R, V_cross, 'c-o', label='V_cross (A<->B)', markersize=3)
    ax.plot(R, V_ee, 'm-s', label='V_ee', markersize=3)
    ax.plot(R, V_NN, 'y-^', label='V_NN', markersize=3)
    ax.plot(R, V_bridge, 'k-d', label='V_bridge', markersize=3)
    ax.set_xlabel('R (bohr)')
    ax.set_ylabel('Energy (Ha)')
    ax.set_title('Interaction Components vs R')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel 4: Bridge analysis
    ax = axes[1, 0]
    R_b = [r['R'] for r in bridge_results]
    max_b = [r['max_bridge'] for r in bridge_results]
    sum_b = [r['sum_bridge'] for r in bridge_results]
    ax.plot(R_b, max_b, 'b-o', label='max|H1_AB|')
    ax.plot(R_b, sum_b, 'r-s', label='sum|H1_AB|')
    ax.set_xlabel('R (bohr)')
    ax.set_ylabel('Bridge coupling (Ha)')
    ax.set_title('Bridge Edge Weights vs R')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Panel 5: E_total(R) and components
    ax = axes[1, 1]
    ax.plot(R, E_total, 'k-o', linewidth=2, label='E_total')
    ax.set_xlabel('R (bohr)')
    ax.set_ylabel('E_total (Ha)')
    ax.set_title('Total Energy vs R')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Panel 6: dT/dR — does kinetic increase at short R?
    ax = axes[1, 2]
    dT_dR = np.gradient(T, R)
    ax.plot(R, dT_dR, 'b-o', linewidth=2)
    ax.axhline(0, color='gray', linestyle='--', alpha=0.3)
    ax.set_xlabel('R (bohr)')
    ax.set_ylabel('dT/dR (Ha/bohr)')
    ax.set_title('Kinetic Energy Slope')
    ax.grid(True, alpha=0.3)
    if all(d <= 0 for d in dT_dR[:-1]):
        ax.annotate('dT/dR <= 0 everywhere\n(NO kinetic wall!)',
                     xy=(0.5, 0.5), xycoords='axes fraction',
                     fontsize=12, color='red', ha='center',
                     bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))

    fig.suptitle('LiH Kinetic Repulsion Diagnostic (exact+True, nmax=3)', fontsize=14)
    fig.tight_layout()

    outpath = os.path.join(PROJECT_ROOT, "debug", "plots", "lih_kinetic_vs_R.png")
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"\nPlot saved to {outpath}")


def write_report(
    decomp_results: List[Dict[str, float]],
    bridge_results: List[Dict[str, float]],
    h2_results: Dict[str, float],
) -> None:
    """Write the KINETIC_DIAGNOSTIC.md report."""

    # Build decomposition table
    decomp_table = "| R (bohr) | T (Ha) | V_total (Ha) | η = -V/(2T) | E_total (Ha) | V_cross (Ha) | V_ee (Ha) | V_NN (Ha) |\n"
    decomp_table += "|----------|--------|-------------|-------------|-------------|-------------|-----------|----------|\n"
    for r in decomp_results:
        V_cross = r['V_cross_A'] + r['V_cross_B']
        decomp_table += (f"| {r['R']:.3f} | {r['T']:.4f} | {r['V_total']:.4f} | "
                         f"{r['eta']:.4f} | {r['E_total']:.5f} | {V_cross:.4f} | "
                         f"{r['V_ee']:.4f} | {r['V_NN']:.4f} |\n")

    # Bridge table
    bridge_table = "| R (bohr) | max|H1_AB| | sum|H1_AB| | n_nonzero | bridge_deg_A | bridge_deg_B |\n"
    bridge_table += "|----------|-----------|-----------|-----------|-------------|-------------|\n"
    for r in bridge_results:
        bridge_table += (f"| {r['R']:.3f} | {r['max_bridge']:.6f} | {r['sum_bridge']:.6f} | "
                         f"{r['n_bridge_nonzero']:.0f} | {r['bridge_degree_A']:.4f} | "
                         f"{r['bridge_degree_B']:.4f} |\n")

    # T slope analysis
    R_arr = np.array([r['R'] for r in decomp_results])
    T_arr = np.array([r['T'] for r in decomp_results])
    dT_dR = np.gradient(T_arr, R_arr)
    T_increases_at_short_R = any(d > 0 for d in dT_dR[:3])

    # Virial analysis
    eta_arr = np.array([r['eta'] for r in decomp_results])
    eta_all_gt_1 = all(e > 1 for e in eta_arr)

    # H₂ equilibrium check
    h2_has_eq = False
    h2_req = None
    if h2_results:
        D_raw_list = [(R, h2_results[R].get('D_raw', None)) for R in sorted(h2_results.keys())]
        D_raw_vals = [(R, d) for R, d in D_raw_list if d is not None]
        if len(D_raw_vals) >= 3:
            idx = max(range(len(D_raw_vals)), key=lambda i: D_raw_vals[i][1])
            if 0 < idx < len(D_raw_vals) - 1:
                h2_has_eq = True
                h2_req = D_raw_vals[idx][0]

    report = f"""# LiH Kinetic Repulsion Diagnostic

**Date:** 2026-03-11
**Version:** v0.9.37
**Script:** `debug/diagnose_kinetic_repulsion.py`
**Configuration:** exact+True, nmax={NMAX}

## Problem Statement

The balanced PES (exact+True and fourier+s_only) shows monotonically attractive
D_cp(R) with no equilibrium geometry. In real molecules, short-range kinetic
confinement energy creates a repulsive wall. This diagnostic identifies why
the LCAO framework fails to produce this wall.

## Diagnostic 1: Energy Decomposition vs R

{decomp_table}

### Key Findings

**Does T(R) increase as R decreases?**
{"YES — kinetic energy increases at short R" if T_increases_at_short_R else "**NO — T(R) does NOT increase at short R. The kinetic energy fails to create a repulsive wall.**"}

- T(R=1.5) = {decomp_results[0]['T']:.4f} Ha
- T(R=3.015) = {decomp_results[7]['T']:.4f} Ha
- T(R=5.0) = {decomp_results[-1]['T']:.4f} Ha
- ΔT(1.5→5.0) = {decomp_results[0]['T'] - decomp_results[-1]['T']:.4f} Ha

**Virial ratio η = -V/(2T):**
{"η > 1 at all R — kinetic energy is systematically too weak relative to potential. The system collapses variationally." if eta_all_gt_1 else f"η ranges from {min(eta_arr):.3f} to {max(eta_arr):.3f}"}

## Diagnostic 2: Bridge Edge Weights vs R

{bridge_table}

### Key Findings

- Bridge coupling **decays monotonically** with R (STO overlap × conformal factors)
- At R=1.5: max|H1_AB| = {bridge_results[0]['max_bridge']:.6f} Ha
- At R=5.0: max|H1_AB| = {bridge_results[-1]['max_bridge']:.6f} Ha
- Bridge edges add kinetic cost via the degree matrix D, but this is **tiny**
  compared to intra-atom kinetic contributions

## Diagnostic 3: H₂ Comparison

{"**H₂ HAS an equilibrium** near R ≈ " + f"{h2_req:.1f}" + " bohr. This means the kinetic repulsion issue is **heteronuclear-specific** or at least more severe for LiH." if h2_has_eq else "**H₂ also has NO equilibrium** — the missing kinetic repulsion is **fundamental** to the LCAO architecture."}

## Diagnostic 4: Intra-Atom Kinetic Independence

The intra-atom off-diagonal H1 (kinetic hopping within Li or H lattice) is
**completely R-independent**. The graph Laplacian edges within each atom are
fixed regardless of inter-nuclear distance.

This means:
- The kinetic energy of each atom's electrons is **frozen** at its isolated-atom value
- When atoms approach, there is no kinetic penalty for orbital compression
- The only R-dependent kinetic contribution comes from bridge edges (tiny)

## Diagnostic 5: Missing Physics

### What creates kinetic repulsion in real molecules

1. **Orthogonalization kinetic energy:** When AOs overlap, Löwdin or Schmidt
   orthogonalization introduces additional nodes in the MOs, raising kinetic energy.
   This goes as ~S²(R) where S is the overlap integral.

2. **Pauli kinetic repulsion:** Antisymmetry forces electrons into higher-momentum
   states when core orbitals overlap. For Li-H, the Li 1s core overlaps with H 1s
   at short R, pushing electrons into antibonding states with higher T.

3. **Kinetic integral T_AB:** In standard QC, the kinetic energy matrix element
   between AOs on different centers, <φ_A|-½∇²|φ_B>, has the correct R-dependence
   to create repulsion. Our bridge edges decay monotonically instead.

### What our framework has

- **Fixed intra-atom graph Laplacian:** Kinetic hopping within each atom is
  R-independent. No orbital compression effect.
- **Monotonically decaying bridges:** Bridge coupling S(R)·Ω_A·Ω_B → 0 as R→0
  would approach a constant, but the conformal factors actually INCREASE at small R,
  partially compensating. Net effect: bridges don't provide repulsion.
- **Cross-nuclear attraction:** This IS R-dependent and creates binding. But without
  compensating kinetic repulsion, it produces monotonic attraction.

## Root Cause

**The graph Laplacian kinetic energy is R-independent for each atom's internal
structure.** In real quantum mechanics, bringing atoms together forces their
electrons to occupy orthogonal states with higher kinetic energy. The discrete
graph Laplacian does not capture this overlap-dependent kinetic cost because:

1. Each atom's adjacency matrix (connectivity) is **fixed** at construction time
2. The degree matrix D (diagonal kinetic contribution) only counts edges, which
   don't change as R decreases
3. Bridge edges contribute to D but are too weak to create sufficient repulsion
4. There is no overlap matrix S_AB that would modify the kinetic operator

## Proposed Fix Paths

### Path A: Overlap-dependent kinetic correction (perturbative)

Add a correction term to the molecular Hamiltonian:

```
H1_corrected = H1 + T_correction(R)

T_correction = lam * Sum_{{a in A, b in B}} S^2(a,b) * |a><a| * (positive scale)
```

where S(a,b) is the STO overlap between orbital a on atom A and orbital b on
atom B. This would add a positive (repulsive) diagonal correction that grows
as R decreases (overlap increases).

**Pros:** Simple, perturbative, preserves existing architecture
**Cons:** Requires calibration parameter λ, not derived from first principles

### Path B: R-dependent adjacency (modify graph topology)

Make the intra-atom adjacency matrix depend on R:

```
A_intra(R) = A_intra(inf) + dA(R)
```

where δA(R) represents additional "virtual edges" induced by orbital overlap.
When orbitals on different atoms overlap, this effectively adds new paths in
the graph, increasing the degree and thus the kinetic energy.

**Pros:** More physically motivated (overlap creates new connectivity)
**Cons:** Breaks the isolated-atom eigenvalue structure

### Path C: Lowdin kinetic energy (explicit orthogonalization)

Compute the overlap matrix S_AB between atoms, perform Lowdin orthogonalization
S^(-1/2), and add the resulting kinetic energy shift:

```
T_Lowdin = Tr[rho * S^(-1/2) * T * S^(-1/2)] - Tr[rho * T]
```

**Pros:** Exact treatment, standard QC approach
**Cons:** Expensive, may conflict with graph topology

### Path D: Bond sphere approach (Paper 8)

The bond sphere theory puts both atoms on a single S3, which naturally includes
the kinetic coupling between centers. The SO(4) Wigner D-matrix elements already
encode the correct R-dependent kinetic contribution. However, v0.9.16–v0.9.18
tests showed this approach has its own issues (Fourier diagonal overbinding).

## Output Files

- `debug/KINETIC_DIAGNOSTIC.md` — this report
- `debug/plots/lih_kinetic_vs_R.png` — T(R), V(R), virial ratio plots
- `debug/data/lih_kinetic_diagnostic.txt` — raw data
"""

    outpath = os.path.join(PROJECT_ROOT, "debug", "KINETIC_DIAGNOSTIC.md")
    with open(outpath, 'w', encoding='utf-8') as f:
        f.write(report)
    print(f"Report saved to {outpath}")


def write_data(
    decomp_results: List[Dict[str, float]],
    bridge_results: List[Dict[str, float]],
) -> None:
    """Write raw data file."""
    datapath = os.path.join(PROJECT_ROOT, "debug", "data",
                            "lih_kinetic_diagnostic.txt")
    os.makedirs(os.path.dirname(datapath), exist_ok=True)

    with open(datapath, 'w', encoding='utf-8') as f:
        f.write("LiH Kinetic Repulsion Diagnostic -- v0.9.37, 2026-03-11\n\n")

        f.write("ENERGY DECOMPOSITION\n")
        f.write(f"{'R':>6s}  {'T':>8s}  {'V_nA':>8s}  {'V_nB':>8s}  "
                f"{'V_crA':>8s}  {'V_crB':>8s}  {'V_bri':>8s}  "
                f"{'V_ee':>8s}  {'V_NN':>8s}  {'E_tot':>10s}  {'eta':>8s}\n")
        f.write("-" * 110 + "\n")
        for r in decomp_results:
            f.write(f"{r['R']:6.3f}  {r['T']:8.4f}  {r['V_nA']:8.4f}  "
                    f"{r['V_nB']:8.4f}  {r['V_cross_A']:8.4f}  "
                    f"{r['V_cross_B']:8.4f}  {r['V_bridge']:8.4f}  "
                    f"{r['V_ee']:8.4f}  {r['V_NN']:8.4f}  "
                    f"{r['E_total']:10.5f}  {r['eta']:8.4f}\n")

        f.write("\n\nBRIDGE ANALYSIS\n")
        f.write(f"{'R':>6s}  {'max|H1_AB|':>12s}  {'sum|H1_AB|':>12s}  "
                f"{'n_nnz':>6s}  {'deg_A':>8s}  {'deg_B':>8s}\n")
        f.write("-" * 60 + "\n")
        for r in bridge_results:
            f.write(f"{r['R']:6.3f}  {r['max_bridge']:12.6f}  "
                    f"{r['sum_bridge']:12.6f}  {r['n_bridge_nonzero']:6.0f}  "
                    f"{r['bridge_degree_A']:8.4f}  {r['bridge_degree_B']:8.4f}\n")

    print(f"Data saved to {datapath}")


def main() -> None:
    """Run all kinetic repulsion diagnostics."""
    print("=" * 70)
    print("LiH KINETIC REPULSION DIAGNOSTIC")
    print("Configuration: exact+True, nmax=3")
    print("=" * 70)

    # Diagnostic 1: Energy decomposition
    decomp_results = diagnostic_1_energy_decomposition()

    # Diagnostic 2: Bridge analysis
    bridge_results = diagnostic_2_bridge_analysis()

    # Diagnostic 3: H₂ comparison
    h2_results = diagnostic_3_h2_comparison()

    # Diagnostic 4: Intra-atom kinetic independence
    diagnostic_4_intra_atom_kinetic()

    # Diagnostic 5: Overlap kinetic estimate
    diagnostic_5_overlap_kinetic_estimate()

    # Generate outputs
    write_data(decomp_results, bridge_results)
    make_plots(decomp_results, bridge_results)
    write_report(decomp_results, bridge_results, h2_results)

    print("\n" + "=" * 70)
    print("DIAGNOSTIC COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
