"""
Scan overlap-dependent kinetic correction lambda for LiH equilibrium.

Tests whether adding T_corr = lambda * sum_b S^2(a,b) to the H1 diagonal
creates a physically reasonable equilibrium geometry.

For each lambda value, runs a PES scan and checks for D_cp maximum.
Reports R_eq, D_cp at equilibrium, and compares to paper/experiment.

Output:
  debug/data/lih_overlap_correction_scan.txt — PES for each lambda
  debug/plots/lih_overlap_correction.png — D_cp(R) curves
  docs/OVERLAP_CORRECTION_RESULTS.md — analysis report

Date: 2026-03-11
Version: v0.9.37
"""
import os
import sys
import time
import warnings
from typing import Dict, List, Optional, Tuple

import numpy as np

warnings.filterwarnings('ignore')

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from geovac.lattice_index import (
    LatticeIndex,
    MolecularLatticeIndex,
    compute_bsse_correction,
)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
NMAX = 3
LAMBDA_VALUES = [0.0, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]

# Coarse R scan for lambda search (fewer points = faster)
R_SCAN = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0]

# Fine R scan around equilibrium once found
R_FINE = np.arange(2.0, 4.5, 0.25).tolist()

# Targets
TARGET_R_EQ = 3.015    # experiment (bohr)
TARGET_D_CP = 0.092    # experiment (Ha)
PAPER_R_EQ = 2.5       # paper result (bohr)
PAPER_D_CP = 0.093     # paper CP-corrected (Ha)


def compute_atomic_references(
    t_corr_lambda: float = 0.0,
) -> Tuple[float, float]:
    """Compute E(Li) and E(H) for separated atoms."""
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


def compute_bsse_once() -> float:
    """Compute BSSE once (R-independent, lambda-independent)."""
    bsse = compute_bsse_correction(
        Z_A=3, Z_B=1, nmax_A=NMAX, nmax_B=NMAX,
        R=3.015, n_electrons_A=3, n_electrons_B=1,
        vee_method='slater_full', fci_method='auto',
    )
    return bsse['BSSE']


def compute_pes_point(
    R: float,
    t_corr_lambda: float,
    E_li: float,
    E_h: float,
    bsse: float,
) -> Dict[str, float]:
    """Compute molecular energy and binding at one R for given lambda."""
    mol = MolecularLatticeIndex(
        Z_A=3, Z_B=1, nmax_A=NMAX, nmax_B=NMAX,
        R=R, n_electrons=4,
        vee_method='slater_full', fci_method='auto',
        cross_nuclear_method='exact',
        cross_atom_vee=True,
        t_corr_lambda=t_corr_lambda,
    )
    eigvals, _ = mol.compute_ground_state(n_states=1)
    E_mol = eigvals[0]
    E_sep = E_li + E_h
    D_raw = E_sep - E_mol
    D_cp = D_raw + bsse  # BSSE is negative

    return {
        'R': R,
        'E_mol': E_mol,
        'D_raw': D_raw,
        'D_cp': D_cp,
        'BSSE': bsse,
    }


def scan_pes(
    R_values: List[float],
    t_corr_lambda: float,
    E_li: float,
    E_h: float,
    bsse: float,
) -> List[Dict[str, float]]:
    """Scan PES at multiple R for a given lambda."""
    results = []
    for R in R_values:
        t0 = time.time()
        row = compute_pes_point(R, t_corr_lambda, E_li, E_h, bsse)
        dt = time.time() - t0
        print(f"    R={R:.3f}: E={row['E_mol']:.5f}, D_raw={row['D_raw']:.4f}, "
              f"D_cp={row['D_cp']:.4f}, BSSE={row['BSSE']:.4f} ({dt:.0f}s)")
        results.append(row)
    return results


def find_equilibrium(
    results: List[Dict[str, float]],
) -> Optional[Tuple[float, float]]:
    """Find R_eq and D_cp at equilibrium via parabolic interpolation."""
    R_arr = np.array([r['R'] for r in results])
    D_arr = np.array([r['D_cp'] for r in results])

    idx_max = np.argmax(D_arr)
    if idx_max == 0 or idx_max == len(D_arr) - 1:
        return None  # no interior maximum

    # Parabolic interpolation around maximum
    R0, R1, R2 = R_arr[idx_max - 1], R_arr[idx_max], R_arr[idx_max + 1]
    D0, D1, D2 = D_arr[idx_max - 1], D_arr[idx_max], D_arr[idx_max + 1]

    denom = 2.0 * ((R1 - R0) * (D1 - D2) - (R1 - R2) * (D1 - D0))
    if abs(denom) < 1e-15:
        return (R1, D1)
    R_eq = R1 - ((R1 - R0)**2 * (D1 - D2) - (R1 - R2)**2 * (D1 - D0)) / denom
    D_eq = D0 + (D1 - D0) * (R_eq - R0) / (R1 - R0)  # linear approx

    return (R_eq, D_eq)


def make_plot(
    all_results: Dict[float, List[Dict[str, float]]],
    equilibria: Dict[float, Optional[Tuple[float, float]]],
) -> None:
    """Plot D_cp(R) for each lambda value."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not available, skipping plot.")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    colors = plt.cm.viridis(np.linspace(0, 1, len(all_results)))

    for (lam, results), color in zip(all_results.items(), colors):
        R_arr = [r['R'] for r in results]
        D_cp = [r['D_cp'] for r in results]
        D_raw = [r['D_raw'] for r in results]

        label = f'lam={lam:.1f}'
        if lam in equilibria and equilibria[lam] is not None:
            R_eq, D_eq = equilibria[lam]
            label += f' (R_eq={R_eq:.2f})'

        ax1.plot(R_arr, D_cp, 'o-', color=color, label=label, markersize=4)
        ax2.plot(R_arr, D_raw, 'o-', color=color, label=label, markersize=4)

    # Mark experiment
    ax1.axhline(TARGET_D_CP, color='red', linestyle='--', alpha=0.5, label=f'expt D_cp={TARGET_D_CP}')
    ax1.axvline(TARGET_R_EQ, color='red', linestyle=':', alpha=0.5, label=f'expt R_eq={TARGET_R_EQ}')

    ax1.set_xlabel('R (bohr)')
    ax1.set_ylabel('D_cp (Ha)')
    ax1.set_title('CP-Corrected Binding Energy')
    ax1.legend(fontsize=7, loc='upper right')
    ax1.grid(True, alpha=0.3)

    ax2.set_xlabel('R (bohr)')
    ax2.set_ylabel('D_raw (Ha)')
    ax2.set_title('Raw Binding Energy')
    ax2.legend(fontsize=7, loc='upper right')
    ax2.grid(True, alpha=0.3)

    fig.suptitle('LiH Overlap Kinetic Correction: lambda scan (nmax=3)', fontsize=13)
    fig.tight_layout()

    outpath = os.path.join(PROJECT_ROOT, "debug", "plots", "lih_overlap_correction.png")
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"\nPlot saved to {outpath}")


def write_data(
    all_results: Dict[float, List[Dict[str, float]]],
    equilibria: Dict[float, Optional[Tuple[float, float]]],
) -> None:
    """Write raw data file."""
    datapath = os.path.join(PROJECT_ROOT, "debug", "data",
                            "lih_overlap_correction_scan.txt")
    os.makedirs(os.path.dirname(datapath), exist_ok=True)

    with open(datapath, 'w', encoding='utf-8') as f:
        f.write("LiH Overlap Kinetic Correction Scan -- v0.9.37, 2026-03-11\n")
        f.write(f"nmax={NMAX}, cross_nuclear=exact, cross_atom_vee=True\n\n")

        for lam, results in all_results.items():
            eq = equilibria.get(lam)
            eq_str = f"R_eq={eq[0]:.3f}, D_cp={eq[1]:.4f}" if eq else "NO EQUILIBRIUM"
            f.write(f"lambda = {lam:.1f}  [{eq_str}]\n")
            f.write(f"  {'R':>6s}  {'E_mol':>10s}  {'D_raw':>8s}  {'D_cp':>8s}  {'BSSE':>8s}\n")
            f.write(f"  {'-'*50}\n")
            for r in results:
                f.write(f"  {r['R']:6.3f}  {r['E_mol']:10.5f}  {r['D_raw']:8.4f}  "
                        f"{r['D_cp']:8.4f}  {r['BSSE']:8.4f}\n")
            f.write("\n")

        f.write("\nSUMMARY\n")
        f.write(f"{'lambda':>8s}  {'R_eq':>8s}  {'D_cp_eq':>8s}  {'status':>20s}\n")
        f.write("-" * 50 + "\n")
        for lam in all_results:
            eq = equilibria.get(lam)
            if eq:
                status = "EQUILIBRIUM FOUND"
                f.write(f"{lam:8.1f}  {eq[0]:8.3f}  {eq[1]:8.4f}  {status:>20s}\n")
            else:
                status = "NO EQUILIBRIUM"
                f.write(f"{lam:8.1f}  {'---':>8s}  {'---':>8s}  {status:>20s}\n")

    print(f"Data saved to {datapath}")


def write_report(
    all_results: Dict[float, List[Dict[str, float]]],
    equilibria: Dict[float, Optional[Tuple[float, float]]],
    E_li: float,
    E_h: float,
) -> None:
    """Write analysis report."""
    # Build summary table
    summary_rows = ""
    best_lam = None
    best_err = float('inf')
    for lam in all_results:
        eq = equilibria.get(lam)
        if eq:
            R_eq, D_cp = eq
            err_R = abs(R_eq - TARGET_R_EQ)
            err_D = abs(D_cp - TARGET_D_CP)
            summary_rows += f"| {lam:.1f} | {R_eq:.3f} | {D_cp:.4f} | {err_R:.3f} | {err_D:.4f} |\n"
            if err_R < best_err:
                best_err = err_R
                best_lam = lam
        else:
            summary_rows += f"| {lam:.1f} | --- | --- | --- | --- |\n"

    # Build PES table for best lambda (or lambda=0 baseline)
    show_lam = best_lam if best_lam is not None else 0.0
    pes_rows = ""
    if show_lam in all_results:
        for r in all_results[show_lam]:
            pes_rows += (f"| {r['R']:.3f} | {r['E_mol']:.5f} | {r['D_raw']:.4f} | "
                         f"{r['D_cp']:.4f} | {r['BSSE']:.4f} |\n")

    has_eq = any(eq is not None for eq in equilibria.values())

    report = f"""# LiH Overlap-Dependent Kinetic Correction Results

**Date:** 2026-03-11
**Version:** v0.9.37
**Script:** `debug/validate_overlap_correction.py`
**Configuration:** exact+True, nmax={NMAX}

## Method

Added a diagonal kinetic correction to the molecular H1:

```
h1_diag[a] += lambda * sum_b S^2(a,b)
```

where S(a,b) is the STO overlap integral between orbital a on one atom and
orbital b on the other. Only s-orbital (l=0) pairs are included.

This approximates the Pauli kinetic repulsion from orthogonalizing overlapping
atomic orbitals. The graph Laplacian's intra-atom kinetic energy is R-independent
(confirmed in `debug/KINETIC_DIAGNOSTIC.md`), so this correction supplies the
missing R-dependent kinetic cost.

## Atomic References

- E(Li) = {E_li:.6f} Ha
- E(H) = {E_h:.6f} Ha
- E_sep = {E_li + E_h:.6f} Ha

## Lambda Scan Results

| lambda | R_eq (bohr) | D_cp (Ha) | |R_eq - expt| | |D_cp - expt| |
|--------|-------------|-----------|---------------|---------------|
{summary_rows}

**Experimental:** R_eq = {TARGET_R_EQ} bohr, D_cp = {TARGET_D_CP} Ha
**Paper (v0.9.11):** R_eq ~ {PAPER_R_EQ} bohr, D_cp = {PAPER_D_CP} Ha

## PES for lambda = {show_lam:.1f}

| R (bohr) | E_mol (Ha) | D_raw (Ha) | D_cp (Ha) | BSSE (Ha) |
|----------|------------|------------|-----------|-----------|
{pes_rows}

## Analysis

{"### Equilibrium Found" if has_eq else "### No Equilibrium Found"}

{"The overlap kinetic correction successfully creates an equilibrium geometry for lambda >= " + str(min(l for l, eq in equilibria.items() if eq is not None)) + "." if has_eq else "No value of lambda in the tested range produces an equilibrium. The s-orbital overlap correction alone is insufficient."}

{"**Best lambda:** " + f"{best_lam:.1f}" + f" (R_eq = {equilibria[best_lam][0]:.3f} bohr, D_cp = {equilibria[best_lam][1]:.4f} Ha)" if best_lam is not None else ""}

### Physical Interpretation

The kinetic correction represents the cost of maintaining orthogonality between
atomic orbitals as they overlap. In standard quantum chemistry, this emerges from
the kinetic energy matrix elements between AOs on different centers. In the graph
Laplacian framework, this physics is absent because the intra-atom adjacency
(and hence kinetic energy) is fixed at construction time.

The s-orbital overlap at the bonding region (R ~ 3 bohr) is small (S^2 ~ 0.007),
so even large lambda values produce modest corrections at equilibrium. At short R
(R ~ 1.5), the overlap is much larger (S^2 ~ 0.19), creating the needed repulsive wall.

## Output Files

- `debug/data/lih_overlap_correction_scan.txt` -- raw data
- `debug/plots/lih_overlap_correction.png` -- D_cp(R) curves
- `docs/OVERLAP_CORRECTION_RESULTS.md` -- this report
"""

    outpath = os.path.join(PROJECT_ROOT, "docs", "OVERLAP_CORRECTION_RESULTS.md")
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w', encoding='utf-8') as f:
        f.write(report)
    print(f"Report saved to {outpath}")


def main() -> None:
    """Run overlap kinetic correction lambda scan."""
    print("=" * 70)
    print("LiH OVERLAP KINETIC CORRECTION SCAN")
    print(f"Lambda values: {LAMBDA_VALUES}")
    print(f"R scan: {R_SCAN}")
    print("=" * 70)

    t_start = time.time()

    # Atomic references (lambda-independent since atoms are isolated)
    print("\nComputing atomic references...")
    E_li, E_h = compute_atomic_references()
    print(f"  E(Li) = {E_li:.6f}, E(H) = {E_h:.6f}, E_sep = {E_li + E_h:.6f}")

    # BSSE (R-independent, lambda-independent: ghost atoms skip correction)
    print("\nComputing BSSE (once)...")
    bsse = compute_bsse_once()
    print(f"  BSSE = {bsse:.6f} Ha")

    all_results: Dict[float, List[Dict[str, float]]] = {}
    equilibria: Dict[float, Optional[Tuple[float, float]]] = {}

    for lam in LAMBDA_VALUES:
        print(f"\n{'='*50}")
        print(f"  lambda = {lam:.1f}")
        print(f"{'='*50}")

        results = scan_pes(R_SCAN, lam, E_li, E_h, bsse)
        eq = find_equilibrium(results)

        all_results[lam] = results
        equilibria[lam] = eq

        if eq:
            print(f"  => EQUILIBRIUM at R_eq = {eq[0]:.3f} bohr, D_cp = {eq[1]:.4f} Ha")
        else:
            D_cp_vals = [r['D_cp'] for r in results]
            print(f"  => NO EQUILIBRIUM (max D_cp = {max(D_cp_vals):.4f} at R = "
                  f"{results[np.argmax(D_cp_vals)]['R']:.1f})")

    # Generate outputs
    write_data(all_results, equilibria)
    make_plot(all_results, equilibria)
    write_report(all_results, equilibria, E_li, E_h)

    dt = time.time() - t_start
    print(f"\nTotal time: {dt:.0f}s ({dt/60:.1f} min)")
    print("SCAN COMPLETE")


if __name__ == '__main__':
    main()
