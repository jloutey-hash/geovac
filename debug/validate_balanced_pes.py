"""
LiH Balanced PES Validation: exact+True vs exact+s_only vs fourier+s_only

The Hamiltonian diagnostic (2026-03-11) revealed that exact+True is BALANCED
(D_raw = 0.225 Ha, matching fourier+s_only) while exact+s_only is IMBALANCED
(D_raw = 0.776 Ha). This script validates the PES shape for all three configs.

Key questions:
  Q1: Does exact+True have an equilibrium geometry?
  Q2: Does the PES shape match fourier+s_only?
  Q3: What is D_cp at equilibrium?
  Q4: Is the 17% geometric contraction (R_eq ~ 2.5 vs expt 3.015) still present?

Output:
  debug/data/lih_balanced_pes.txt         — full PES data for all 3 configs
  debug/plots/lih_balanced_pes_comparison.png — overlay plot of all 3 D_cp(R)
  docs/LIH_BALANCED_PES_VALIDATION.md     — validation report

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
# Constants
# ---------------------------------------------------------------------------
D_E_EXPT = 0.092       # Ha, Huber & Herzberg 1979
R_EQ_EXPT = 3.015      # bohr
NMAX = 3

# Paper targets
PAPER_D_CP = 0.093     # Ha
PAPER_REQ = 2.5        # bohr (approximate)

# Scan grid
R_SCAN = np.array([1.5, 1.8, 2.0, 2.2, 2.5, 2.8, 3.0, 3.015, 3.2, 3.5, 4.0, 5.0])

# Three configurations to compare
CONFIGS = [
    {
        'name': 'exact+True',
        'cross_nuclear_method': 'exact',
        'cross_atom_vee': True,
        'color': 'blue',
        'marker': 's',
    },
    {
        'name': 'exact+s_only',
        'cross_nuclear_method': 'exact',
        'cross_atom_vee': 's_only',
        'color': 'red',
        'marker': 'o',
    },
    {
        'name': 'fourier+s_only',
        'cross_nuclear_method': 'fourier',
        'cross_atom_vee': 's_only',
        'color': 'green',
        'marker': '^',
    },
]


def compute_atomic_energies() -> Tuple[float, float]:
    """Compute isolated Li and H ground-state energies (own basis, nmax=3)."""
    li = LatticeIndex(
        n_electrons=3, max_n=NMAX, nuclear_charge=3,
        vee_method='slater_full', h1_method='exact', fci_method='auto',
    )
    E_li = li.compute_ground_state(n_states=1)[0][0]

    h = LatticeIndex(
        n_electrons=1, max_n=NMAX, nuclear_charge=1,
        vee_method='slater_full', h1_method='exact', fci_method='auto',
    )
    E_h = h.compute_ground_state(n_states=1)[0][0]

    return E_li, E_h


def compute_pes_point(
    R: float,
    cross_nuclear_method: str,
    cross_atom_vee,
    E_li_own: float,
    E_h_own: float,
) -> Dict[str, float]:
    """Compute E_mol, D_raw, D_cp, BSSE at a single R for a given config."""
    E_sep_own = E_li_own + E_h_own

    # Molecular energy
    mol = MolecularLatticeIndex(
        Z_A=3, Z_B=1, nmax_A=NMAX, nmax_B=NMAX,
        R=R, n_electrons=4,
        vee_method='slater_full', fci_method='auto',
        cross_nuclear_method=cross_nuclear_method,
        cross_atom_vee=cross_atom_vee,
    )
    E_mol = mol.compute_ground_state(n_states=1)[0][0]

    # Ghost Li (Li in full basis, H nucleus removed)
    li_g = MolecularLatticeIndex(
        Z_A=3, Z_B=0, nmax_A=NMAX, nmax_B=NMAX,
        R=R, n_electrons=3,
        vee_method='slater_full', fci_method='auto',
        cross_nuclear_method=cross_nuclear_method,
        cross_atom_vee=cross_atom_vee,
    )
    E_li_g = li_g.compute_ground_state(n_states=1)[0][0]

    # Ghost H (H in full basis, Li nucleus removed)
    h_g = MolecularLatticeIndex(
        Z_A=0, Z_B=1, nmax_A=NMAX, nmax_B=NMAX,
        R=R, n_electrons=1,
        vee_method='slater_full', fci_method='auto',
        cross_nuclear_method=cross_nuclear_method,
        cross_atom_vee=cross_atom_vee,
    )
    E_h_g = h_g.compute_ground_state(n_states=1)[0][0]

    E_ghost = E_li_g + E_h_g
    D_raw = E_sep_own - E_mol
    D_cp = E_ghost - E_mol
    BSSE = E_ghost - E_sep_own

    return {
        'E_mol': E_mol,
        'D_raw': D_raw,
        'D_cp': D_cp,
        'BSSE': BSSE,
    }


def scan_config(
    config: dict,
    R_values: np.ndarray,
    E_li_own: float,
    E_h_own: float,
) -> Dict[str, np.ndarray]:
    """Scan PES for a single configuration."""
    name = config['name']
    n_pts = len(R_values)
    E_mol = np.zeros(n_pts)
    D_raw = np.zeros(n_pts)
    D_cp = np.zeros(n_pts)
    BSSE = np.zeros(n_pts)

    for i, R in enumerate(R_values):
        print(f"  [{name}] R={R:.3f} ({i+1}/{n_pts})...", end="", flush=True)
        pt = compute_pes_point(
            R, config['cross_nuclear_method'], config['cross_atom_vee'],
            E_li_own, E_h_own,
        )
        E_mol[i] = pt['E_mol']
        D_raw[i] = pt['D_raw']
        D_cp[i] = pt['D_cp']
        BSSE[i] = pt['BSSE']
        print(f" D_cp={D_cp[i]:.4f}")

    return {'R': R_values, 'E_mol': E_mol, 'D_raw': D_raw,
            'D_cp': D_cp, 'BSSE': BSSE}


def find_pes_maximum(
    R_values: np.ndarray,
    values: np.ndarray,
) -> Tuple[float, float]:
    """Find maximum via parabolic interpolation (for D_raw, D_cp)."""
    idx = np.argmax(values)

    if idx == 0 or idx == len(R_values) - 1:
        return R_values[idx], values[idx]

    R_m, R_0, R_p = R_values[idx - 1], R_values[idx], R_values[idx + 1]
    E_m, E_0, E_p = values[idx - 1], values[idx], values[idx + 1]

    denom = 2.0 * ((R_m - R_0) * (E_m - E_p) - (R_m - R_p) * (E_m - E_0))
    if abs(denom) < 1e-15:
        return R_0, E_0

    R_eq = R_m - (
        (R_m - R_0)**2 * (E_m - E_p) - (R_m - R_p)**2 * (E_m - E_0)
    ) / denom

    A = np.array([
        [R_m**2, R_m, 1.0],
        [R_0**2, R_0, 1.0],
        [R_p**2, R_p, 1.0],
    ])
    coeffs = np.linalg.solve(A, np.array([E_m, E_0, E_p]))
    E_eq = coeffs[0] * R_eq**2 + coeffs[1] * R_eq + coeffs[2]

    return R_eq, E_eq


def make_comparison_plot(
    all_results: List[Tuple[dict, Dict[str, np.ndarray]]],
    equilibria: List[Tuple[float, float]],
    outpath: str,
) -> None:
    """Generate overlay plot of D_cp(R) for all 3 configurations."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not available, skipping plot.")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel 1: D_cp(R) for all configs
    ax = axes[0]
    for (cfg, res), (R_eq, D_eq) in zip(all_results, equilibria):
        ax.plot(res['R'], res['D_cp'], color=cfg['color'], marker=cfg['marker'],
                markersize=4, linewidth=1.5,
                label=f"{cfg['name']} (R_eq={R_eq:.2f}, D_cp={D_eq:.4f})")
    ax.axhline(D_E_EXPT, color='gray', linestyle='--', alpha=0.7,
               label=f'Expt D_e={D_E_EXPT}')
    ax.axvline(R_EQ_EXPT, color='gray', linestyle=':', alpha=0.5,
               label=f'Expt R_eq={R_EQ_EXPT}')
    ax.set_xlabel('R (bohr)', fontsize=12)
    ax.set_ylabel('D_cp (Ha)', fontsize=12)
    ax.set_title('LiH CP-Corrected Binding Energy: 3 Configurations', fontsize=12)
    ax.legend(fontsize=8, loc='best')
    ax.grid(True, alpha=0.3)

    # Panel 2: D_raw(R) for all configs
    ax = axes[1]
    for (cfg, res), _ in zip(all_results, equilibria):
        R_eq_raw, D_raw_eq = find_pes_maximum(res['R'], res['D_raw'])
        ax.plot(res['R'], res['D_raw'], color=cfg['color'], marker=cfg['marker'],
                markersize=4, linewidth=1.5,
                label=f"{cfg['name']} (R_eq={R_eq_raw:.2f}, D_raw={D_raw_eq:.4f})")
    ax.axhline(D_E_EXPT, color='gray', linestyle='--', alpha=0.7,
               label=f'Expt D_e={D_E_EXPT}')
    ax.axvline(R_EQ_EXPT, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('R (bohr)', fontsize=12)
    ax.set_ylabel('D_raw (Ha)', fontsize=12)
    ax.set_title('LiH Raw Binding Energy: 3 Configurations', fontsize=12)
    ax.legend(fontsize=8, loc='best')
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Plot saved to {outpath}")


def main() -> None:
    """Run balanced PES validation."""
    print("=" * 70)
    print("LiH Balanced PES Validation")
    print("Comparing: exact+True vs exact+s_only vs fourier+s_only")
    print("=" * 70)

    # Atomic reference energies
    print("\n--- Atomic Reference Energies (nmax=3) ---")
    E_li, E_h = compute_atomic_energies()
    E_sep = E_li + E_h
    print(f"  E(Li) = {E_li:.6f} Ha")
    print(f"  E(H)  = {E_h:.6f} Ha")
    print(f"  E_sep = {E_sep:.6f} Ha")

    # Scan all three configurations
    all_results: List[Tuple[dict, Dict[str, np.ndarray]]] = []
    for cfg in CONFIGS:
        print(f"\n--- PES Scan: {cfg['name']} ({len(R_SCAN)} points) ---")
        res = scan_config(cfg, R_SCAN, E_li, E_h)
        all_results.append((cfg, res))

    # Find equilibria (max of D_cp)
    equilibria: List[Tuple[float, float]] = []
    for cfg, res in all_results:
        R_eq, D_eq = find_pes_maximum(res['R'], res['D_cp'])
        equilibria.append((R_eq, D_eq))

    # ---------------------------------------------------------------
    # Print comparison table
    # ---------------------------------------------------------------
    print("\n" + "=" * 90)
    print("COMPARISON TABLE: D_cp(R) for all three configurations")
    print("=" * 90)

    header = (f"{'R (bohr)':>10s}  "
              f"{'exact+True':>12s}  {'exact+s_only':>12s}  {'fourier+s_only':>14s}  "
              f"{'|eT-fs|':>8s}")
    print(header)
    print("-" * len(header))

    for i, R in enumerate(R_SCAN):
        d_et = all_results[0][1]['D_cp'][i]
        d_es = all_results[1][1]['D_cp'][i]
        d_fs = all_results[2][1]['D_cp'][i]
        diff_et_fs = abs(d_et - d_fs)
        print(f"{R:10.3f}  {d_et:12.4f}  {d_es:12.4f}  {d_fs:14.4f}  {diff_et_fs:8.4f}")

    # Equilibrium summary
    print("\n" + "=" * 70)
    print("EQUILIBRIUM SUMMARY")
    print("=" * 70)
    print(f"{'Config':<20s}  {'R_eq (bohr)':>12s}  {'D_cp (Ha)':>10s}  "
          f"{'D_cp err%':>10s}  {'R_eq err%':>10s}")
    print("-" * 70)
    for (cfg, res), (R_eq, D_eq) in zip(all_results, equilibria):
        d_err = (D_eq - D_E_EXPT) / D_E_EXPT * 100
        r_err = (R_eq - R_EQ_EXPT) / R_EQ_EXPT * 100
        print(f"{cfg['name']:<20s}  {R_eq:12.3f}  {D_eq:10.4f}  {d_err:+10.1f}%  "
              f"{r_err:+10.1f}%")
    print(f"{'Experiment':<20s}  {R_EQ_EXPT:12.3f}  {D_E_EXPT:10.4f}  {'---':>10s}  "
          f"{'---':>10s}")

    # Answers to 4 key questions
    R_eq_et, D_eq_et = equilibria[0]
    R_eq_fs, D_eq_fs = equilibria[2]
    D_cp_et = all_results[0][1]['D_cp']
    D_cp_fs = all_results[2][1]['D_cp']
    max_diff = np.max(np.abs(D_cp_et - D_cp_fs))

    has_eq_et = equilibria[0][0] > R_SCAN[0] and equilibria[0][0] < R_SCAN[-1]
    shapes_match = max_diff < 0.010
    contraction = (R_eq_et - R_EQ_EXPT) / R_EQ_EXPT * 100

    print("\n" + "=" * 70)
    print("ANSWERS TO KEY QUESTIONS")
    print("=" * 70)
    print(f"Q1: Does exact+True have an equilibrium?")
    if has_eq_et:
        print(f"    YES — R_eq = {R_eq_et:.3f} bohr, D_cp = {D_eq_et:.4f} Ha")
    else:
        print(f"    NO — D_cp is monotonic; max at boundary R={R_SCAN[0]:.1f}")
        print(f"    => Possible issue with l>0 Gaunt coefficients")

    print(f"\nQ2: Does PES shape match fourier+s_only?")
    print(f"    Max |D_cp(exact+True) - D_cp(fourier+s_only)| = {max_diff:.4f} Ha")
    if shapes_match:
        print(f"    YES — curves match within 10 mHa")
    else:
        print(f"    NO — significant shape difference")

    print(f"\nQ3: D_cp at equilibrium?")
    print(f"    exact+True:     D_cp = {D_eq_et:.4f} Ha (expt {D_E_EXPT})")
    print(f"    fourier+s_only: D_cp = {D_eq_fs:.4f} Ha")

    print(f"\nQ4: Geometric contraction?")
    print(f"    R_eq(exact+True) = {R_eq_et:.3f} bohr (expt {R_EQ_EXPT})")
    print(f"    Contraction = {contraction:+.1f}%")

    # ---------------------------------------------------------------
    # Save data file
    # ---------------------------------------------------------------
    data_path = os.path.join(PROJECT_ROOT, "debug", "data",
                             "lih_balanced_pes.txt")
    os.makedirs(os.path.dirname(data_path), exist_ok=True)

    with open(data_path, 'w') as f:
        f.write("LiH Balanced PES Validation — v0.9.37, 2026-03-11\n")
        f.write(f"nmax={NMAX}, E(Li)={E_li:.6f}, E(H)={E_h:.6f}\n\n")

        # Full table
        f.write(f"{'R':>6s}  ")
        for cfg, _ in all_results:
            f.write(f"{'E_mol':>10s}  {'D_raw':>8s}  {'D_cp':>8s}  {'BSSE':>8s}  ")
        f.write("\n")

        # Column headers with config names
        f.write(f"{'':>6s}  ")
        for cfg, _ in all_results:
            f.write(f"{'--- ' + cfg['name'] + ' ---':^38s}  ")
        f.write("\n")
        f.write("-" * (6 + len(all_results) * 42) + "\n")

        for i, R in enumerate(R_SCAN):
            f.write(f"{R:6.3f}  ")
            for cfg, res in all_results:
                f.write(f"{res['E_mol'][i]:10.5f}  {res['D_raw'][i]:8.4f}  "
                        f"{res['D_cp'][i]:8.4f}  {res['BSSE'][i]:8.4f}  ")
            f.write("\n")

        f.write("\n\nEQUILIBRIA\n")
        for (cfg, _), (R_eq, D_eq) in zip(all_results, equilibria):
            f.write(f"  {cfg['name']}: R_eq={R_eq:.3f} bohr, D_cp={D_eq:.4f} Ha\n")

    print(f"\nData saved to {data_path}")

    # ---------------------------------------------------------------
    # Plot
    # ---------------------------------------------------------------
    plot_path = os.path.join(PROJECT_ROOT, "debug", "plots",
                             "lih_balanced_pes_comparison.png")
    make_comparison_plot(all_results, equilibria, plot_path)

    # ---------------------------------------------------------------
    # Validation report
    # ---------------------------------------------------------------
    report_path = os.path.join(PROJECT_ROOT, "docs",
                               "LIH_BALANCED_PES_VALIDATION.md")

    # Build D_cp comparison table for markdown
    table_rows = ""
    for i, R in enumerate(R_SCAN):
        d_et = all_results[0][1]['D_cp'][i]
        d_es = all_results[1][1]['D_cp'][i]
        d_fs = all_results[2][1]['D_cp'][i]
        table_rows += f"| {R:.3f} | {d_et:.4f} | {d_es:.4f} | {d_fs:.4f} |\n"

    R_eq_es, D_eq_es = equilibria[1]

    report = f"""# LiH Balanced PES Validation Report

**Date:** 2026-03-11
**Version:** v0.9.37
**Script:** `debug/validate_balanced_pes.py`
**Status:** {'PASS' if has_eq_et else 'FAIL — no equilibrium'}

## Purpose

Validate the PES shape for the balanced `exact+True` configuration
(exact cross-nuclear + all-l cross-atom V_ee) discovered in the
Hamiltonian diagnostic. Compare against the known `fourier+s_only`
(v0.9.11 baseline) and the imbalanced `exact+s_only`.

## Atomic Reference Energies (nmax={NMAX})

| Atom | E (Ha) |
|------|--------|
| Li   | {E_li:.6f} |
| H    | {E_h:.6f} |
| Sep  | {E_sep:.6f} |

## D_cp(R) Comparison Table

| R (bohr) | exact+True | exact+s_only | fourier+s_only |
|-----------|-----------|-------------|---------------|
{table_rows}

## Equilibrium Summary

| Config | R_eq (bohr) | D_cp (Ha) | D_cp err% | R_eq err% |
|--------|------------|-----------|-----------|-----------|
| exact+True | {R_eq_et:.3f} | {D_eq_et:.4f} | {(D_eq_et - D_E_EXPT) / D_E_EXPT * 100:+.1f}% | {(R_eq_et - R_EQ_EXPT) / R_EQ_EXPT * 100:+.1f}% |
| exact+s_only | {R_eq_es:.3f} | {D_eq_es:.4f} | {(D_eq_es - D_E_EXPT) / D_E_EXPT * 100:+.1f}% | {(R_eq_es - R_EQ_EXPT) / R_EQ_EXPT * 100:+.1f}% |
| fourier+s_only | {R_eq_fs:.3f} | {D_eq_fs:.4f} | {(D_eq_fs - D_E_EXPT) / D_E_EXPT * 100:+.1f}% | {(R_eq_fs - R_EQ_EXPT) / R_EQ_EXPT * 100:+.1f}% |
| Experiment | {R_EQ_EXPT} | {D_E_EXPT} | --- | --- |

## Answers to Key Questions

### Q1: Does exact+True have an equilibrium geometry?

{"**YES** — R_eq = " + f"{R_eq_et:.3f}" + " bohr, D_cp = " + f"{D_eq_et:.4f}" + " Ha" if has_eq_et else "**NO** — D_cp is monotonically attractive. This indicates the all-l V_ee repulsion may be insufficient to create a turnover, or the l>0 Gaunt coefficients need review."}

### Q2: Does the PES shape match fourier+s_only?

Max |D_cp(exact+True) - D_cp(fourier+s_only)| = {max_diff:.4f} Ha

{"**YES** — The two balanced configurations produce nearly identical PES curves (within " + f"{max_diff*1000:.1f}" + " mHa). This confirms that exact cross-nuclear + all-l V_ee achieves the same variational balance as fourier cross-nuclear + s-only V_ee." if shapes_match else "**NO** — Significant shape difference of " + f"{max_diff:.4f}" + " Ha between the two balanced configurations."}

### Q3: What is D_cp at equilibrium?

- exact+True: D_cp = {D_eq_et:.4f} Ha (paper target: {PAPER_D_CP})
- fourier+s_only: D_cp = {D_eq_fs:.4f} Ha
- Estimate from diagnostic: 0.225 - 0.115 = 0.110 Ha

### Q4: Is the 17% geometric contraction still present?

R_eq(exact+True) = {R_eq_et:.3f} bohr vs experiment {R_EQ_EXPT} bohr
Contraction = {contraction:+.1f}% (paper prediction: -17%)

## Output Files

- `debug/data/lih_balanced_pes.txt` — full PES data
- `debug/plots/lih_balanced_pes_comparison.png` — overlay comparison plot
"""

    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    with open(report_path, 'w') as f:
        f.write(report)
    print(f"Report saved to {report_path}")

    print("\n" + "=" * 70)
    print("VALIDATION COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
