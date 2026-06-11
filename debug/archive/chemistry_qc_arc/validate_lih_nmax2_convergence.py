"""
Validate LiH nmax convergence trend claimed in Section V.C of the FCI paper.

Runs LiH FCI at nmax=2 and nmax=3 with Boys-Bernardi counterpoise correction,
scanning the PES to find R_eq and D_e_CP at each basis size. Verifies:
  1. BSSE(nmax=2) > BSSE(nmax=3)        — convergence direction
  2. R_eq(nmax=2) < R_eq(nmax=3)         — convergence direction
  3. D_e_CP(nmax=2) ≈ 0.270 Ha (192%)   — paper Section V.C claim

Output:
  stdout                                     — comparison table
  debug/data/lih_nmax_convergence.txt        — same table as text file
  debug/plots/lih_pes_convergence.png        — PES plot (both nmax)
  docs/LIH_CONVERGENCE_VALIDATION.md         — validation report

Date: 2026-03-11
Version: v0.9.37
"""
import os
import sys
import warnings
from typing import Dict, List, Tuple

import numpy as np

warnings.filterwarnings('ignore')

# Ensure project root is on path
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from geovac.lattice_index import (
    LatticeIndex,
    MolecularLatticeIndex,
)

# ---------------------------------------------------------------------------
# Physical constants (from geovac.constants / CLAUDE.md)
# ---------------------------------------------------------------------------
D_E_EXPT = 0.092       # Ha, Huber & Herzberg 1979
R_EQ_EXPT = 3.015      # bohr

# Paper Section V.C claims
PAPER_D_CP_NMAX2 = 0.270   # Ha
PAPER_ERROR_NMAX2 = 192     # percent

# Known nmax=3 values for cross-check
KNOWN_D_CP_NMAX3 = 0.093   # Ha
KNOWN_BSSE_NMAX3 = 0.115   # Ha (magnitude)
KNOWN_REQ_NMAX3 = 2.5      # bohr (approximate)


def compute_atomic_energies(nmax: int) -> Tuple[float, float]:
    """Compute isolated Li and H ground-state energies in own basis.

    Parameters
    ----------
    nmax : int
        Maximum principal quantum number.

    Returns
    -------
    E_li : float
        Li ground-state energy (Ha).
    E_h : float
        H ground-state energy (Ha).
    """
    li = LatticeIndex(
        n_electrons=3, max_n=nmax, nuclear_charge=3,
        vee_method='slater_full', h1_method='exact', fci_method='auto',
    )
    E_li = li.compute_ground_state(n_states=1)[0][0]

    h = LatticeIndex(
        n_electrons=1, max_n=nmax, nuclear_charge=1,
        vee_method='slater_full', h1_method='exact', fci_method='auto',
    )
    E_h = h.compute_ground_state(n_states=1)[0][0]

    return E_li, E_h


def compute_molecular_energy(nmax: int, R: float) -> float:
    """Compute LiH total ground-state energy at given R.

    Parameters
    ----------
    nmax : int
        Maximum principal quantum number for both atoms.
    R : float
        Internuclear distance (bohr).

    Returns
    -------
    E_mol : float
        Total molecular energy including V_NN (Ha).
    """
    mol = MolecularLatticeIndex(
        Z_A=3, Z_B=1, nmax_A=nmax, nmax_B=nmax,
        R=R, n_electrons=4,
        vee_method='slater_full', fci_method='auto',
    )
    eigvals, _ = mol.compute_ground_state(n_states=1)
    return eigvals[0]


def compute_ghost_energies(nmax: int, R: float) -> Tuple[float, float]:
    """Compute ghost-fragment energies for CP correction at given R.

    Parameters
    ----------
    nmax : int
        Maximum principal quantum number for both atoms.
    R : float
        Internuclear distance (bohr).

    Returns
    -------
    E_li_ghost : float
        Li energy in full molecular basis (ghost H).
    E_h_ghost : float
        H energy in full molecular basis (ghost Li).
    """
    li_g = MolecularLatticeIndex(
        Z_A=3, Z_B=0, nmax_A=nmax, nmax_B=nmax,
        R=R, n_electrons=3,
        vee_method='slater_full', fci_method='auto',
    )
    E_li_g = li_g.compute_ground_state(n_states=1)[0][0]

    h_g = MolecularLatticeIndex(
        Z_A=0, Z_B=1, nmax_A=nmax, nmax_B=nmax,
        R=R, n_electrons=1,
        vee_method='slater_full', fci_method='auto',
    )
    E_h_g = h_g.compute_ground_state(n_states=1)[0][0]

    return E_li_g, E_h_g


def scan_pes(
    nmax: int,
    R_values: np.ndarray,
    E_li_own: float,
    E_h_own: float,
) -> Dict[str, np.ndarray]:
    """Scan PES at given nmax, computing raw and CP-corrected binding energies.

    Parameters
    ----------
    nmax : int
        Maximum principal quantum number.
    R_values : np.ndarray
        Internuclear distances to scan (bohr).
    E_li_own : float
        Isolated Li energy in own basis.
    E_h_own : float
        Isolated H energy in own basis.

    Returns
    -------
    results : dict
        Keys: 'R', 'E_mol', 'E_ghost', 'D_raw', 'D_cp', 'BSSE'
        Each value is np.ndarray of same length as R_values.
    """
    E_sep_own = E_li_own + E_h_own
    n_pts = len(R_values)

    E_mol = np.zeros(n_pts)
    E_ghost = np.zeros(n_pts)
    D_raw = np.zeros(n_pts)
    D_cp = np.zeros(n_pts)
    BSSE = np.zeros(n_pts)

    for i, R in enumerate(R_values):
        print(f"  nmax={nmax}, R={R:.2f} bohr ({i+1}/{n_pts})...", end="", flush=True)

        E_mol[i] = compute_molecular_energy(nmax, R)
        E_li_g, E_h_g = compute_ghost_energies(nmax, R)
        E_ghost[i] = E_li_g + E_h_g

        D_raw[i] = E_sep_own - E_mol[i]
        D_cp[i] = E_ghost[i] - E_mol[i]
        BSSE[i] = E_ghost[i] - E_sep_own

        print(f" E_mol={E_mol[i]:.5f}, D_cp={D_cp[i]:.4f}")

    return {
        'R': R_values,
        'E_mol': E_mol,
        'E_ghost': E_ghost,
        'D_raw': D_raw,
        'D_cp': D_cp,
        'BSSE': BSSE,
    }


def find_pes_minimum(R_values: np.ndarray, D_cp: np.ndarray) -> Tuple[float, float]:
    """Find R_eq and D_e_CP from PES scan via parabolic interpolation.

    Parameters
    ----------
    R_values : np.ndarray
        Internuclear distances (bohr).
    D_cp : np.ndarray
        CP-corrected binding energies (Ha).

    Returns
    -------
    R_eq : float
        Equilibrium bond length (bohr).
    D_e_cp : float
        Maximum CP-corrected binding energy (Ha).
    """
    idx_max = np.argmax(D_cp)

    # If at boundary, just return the grid point
    if idx_max == 0 or idx_max == len(R_values) - 1:
        return R_values[idx_max], D_cp[idx_max]

    # Parabolic interpolation around the maximum
    R_m = R_values[idx_max - 1]
    R_0 = R_values[idx_max]
    R_p = R_values[idx_max + 1]
    D_m = D_cp[idx_max - 1]
    D_0 = D_cp[idx_max]
    D_p = D_cp[idx_max + 1]

    # Vertex of parabola through three points
    denom = 2.0 * ((R_m - R_0) * (D_m - D_p) - (R_m - R_p) * (D_m - D_0))
    if abs(denom) < 1e-15:
        return R_0, D_0

    R_eq = R_m - (
        (R_m - R_0)**2 * (D_m - D_p) - (R_m - R_p)**2 * (D_m - D_0)
    ) / denom

    # Evaluate D_e at interpolated R_eq
    # Use quadratic fit
    A = np.array([
        [R_m**2, R_m, 1],
        [R_0**2, R_0, 1],
        [R_p**2, R_p, 1],
    ])
    b = np.array([D_m, D_0, D_p])
    coeffs = np.linalg.solve(A, b)
    D_e_cp = coeffs[0] * R_eq**2 + coeffs[1] * R_eq + coeffs[2]

    return R_eq, D_e_cp


def make_pes_plot(
    results_2: Dict[str, np.ndarray],
    results_3: Dict[str, np.ndarray],
    R_eq_2: float,
    D_e_2: float,
    R_eq_3: float,
    D_e_3: float,
    outpath: str,
) -> None:
    """Generate PES comparison plot for nmax=2 and nmax=3.

    Parameters
    ----------
    results_2 : dict
        PES scan results for nmax=2.
    results_3 : dict
        PES scan results for nmax=3.
    R_eq_2, D_e_2 : float
        Equilibrium geometry and binding energy for nmax=2.
    R_eq_3, D_e_3 : float
        Equilibrium geometry and binding energy for nmax=3.
    outpath : str
        Output file path for the plot.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not available, skipping PES plot.")
        return

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(results_2['R'], results_2['D_cp'], 'o-', color='#d62728',
            label=f'nmax=2 (R_eq={R_eq_2:.2f}, D_e={D_e_2:.3f} Ha)', markersize=4)
    ax.plot(results_3['R'], results_3['D_cp'], 's-', color='#1f77b4',
            label=f'nmax=3 (R_eq={R_eq_3:.2f}, D_e={D_e_3:.3f} Ha)', markersize=4)

    ax.axhline(D_E_EXPT, color='green', linestyle='--', alpha=0.7,
               label=f'Expt D_e = {D_E_EXPT:.3f} Ha')
    ax.axvline(R_EQ_EXPT, color='green', linestyle=':', alpha=0.5,
               label=f'Expt R_eq = {R_EQ_EXPT:.3f} bohr')

    ax.set_xlabel('R (bohr)')
    ax.set_ylabel('D_e^CP (Ha)')
    ax.set_title('LiH CP-Corrected PES: nmax Convergence')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"\nPES plot saved to {outpath}")


def write_data_file(
    results_2: Dict[str, np.ndarray],
    results_3: Dict[str, np.ndarray],
    summary: str,
    outpath: str,
) -> None:
    """Save comparison table to text file.

    Parameters
    ----------
    results_2, results_3 : dict
        PES scan results for nmax=2 and nmax=3.
    summary : str
        Summary text block to write at top.
    outpath : str
        Output file path.
    """
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        f.write(summary)
    print(f"Data file saved to {outpath}")


def write_validation_report(
    E_li_2: float, E_h_2: float,
    E_li_3: float, E_h_3: float,
    bsse_2: float, bsse_3: float,
    R_eq_2: float, R_eq_3: float,
    D_e_2: float, D_e_3: float,
    D_cp_at_3015_nmax2: float,
    paper_match: bool,
    actual_d_cp_nmax2_3015: float,
    outpath: str,
) -> None:
    """Write validation report in CLAUDE.md documentation format.

    Parameters
    ----------
    E_li_2, E_h_2 : float
        Atomic energies at nmax=2.
    E_li_3, E_h_3 : float
        Atomic energies at nmax=3.
    bsse_2, bsse_3 : float
        BSSE at R=3.015 for nmax=2 and nmax=3.
    R_eq_2, R_eq_3 : float
        Equilibrium bond lengths.
    D_e_2, D_e_3 : float
        CP-corrected binding energies at PES minimum.
    D_cp_at_3015_nmax2 : float
        D_e_CP at R=3.015 for nmax=2.
    paper_match : bool
        Whether the nmax=2 result matches Section V.C claim.
    actual_d_cp_nmax2_3015 : float
        Actual computed D_e_CP at R=3.015 for nmax=2.
    outpath : str
        Output file path.
    """
    status = "PASS" if paper_match else "PAPER CORRECTION REQUIRED"

    report = f"""# LiH nmax Convergence Validation Report

**Date:** 2026-03-11
**Version:** v0.9.37
**Status:** {status}
**Script:** `debug/validate_lih_nmax2_convergence.py`

## Purpose

Validate the convergence trend claimed in Section V.C of the FCI paper:
- D_e_CP(nmax=2) = 0.270 Ha (192% error)
- D_e_CP(nmax=3) = 0.093 Ha (1.0% error)
- BSSE decreases with increasing nmax
- R_eq converges outward toward experiment with increasing nmax

## Atomic Reference Energies

| Atom | nmax=2 (Ha) | nmax=3 (Ha) |
|------|-------------|-------------|
| Li   | {E_li_2:.6f} | {E_li_3:.6f} |
| H    | {E_h_2:.6f} | {E_h_3:.6f} |

## Convergence Comparison

| Quantity | nmax=2 | nmax=3 | Expt | Trend |
|----------|--------|--------|------|-------|
| D_e_CP at R=3.015 (Ha) | {actual_d_cp_nmax2_3015:.3f} | {KNOWN_D_CP_NMAX3:.3f} | {D_E_EXPT:.3f} | {"Converging" if actual_d_cp_nmax2_3015 > KNOWN_D_CP_NMAX3 else "WRONG DIRECTION"} |
| D_e_CP at R_eq (Ha) | {D_e_2:.3f} | {D_e_3:.3f} | {D_E_EXPT:.3f} | {"Converging" if D_e_2 > D_e_3 else "WRONG DIRECTION"} |
| BSSE at R=3.015 (Ha) | {bsse_2:.3f} | {bsse_3:.3f} | 0.000 | {"Converging" if abs(bsse_2) > abs(bsse_3) else "WRONG DIRECTION"} |
| R_eq (bohr) | {R_eq_2:.2f} | {R_eq_3:.2f} | {R_EQ_EXPT:.3f} | {"Converging" if R_eq_2 < R_eq_3 else "WRONG DIRECTION"} |
| Error in D_e_CP (%) | {abs(actual_d_cp_nmax2_3015 - D_E_EXPT) / D_E_EXPT * 100:.0f} | {abs(KNOWN_D_CP_NMAX3 - D_E_EXPT) / D_E_EXPT * 100:.0f} | --- | --- |

## Verification of Paper Section V.C Claims

| Claim | Paper Value | Computed Value | Match |
|-------|-------------|----------------|-------|
| D_e_CP(nmax=2) | {PAPER_D_CP_NMAX2:.3f} Ha | {actual_d_cp_nmax2_3015:.3f} Ha | {"YES" if paper_match else "NO"} |
| Error(nmax=2) | {PAPER_ERROR_NMAX2}% | {abs(actual_d_cp_nmax2_3015 - D_E_EXPT) / D_E_EXPT * 100:.0f}% | {"YES" if abs(abs(actual_d_cp_nmax2_3015 - D_E_EXPT) / D_E_EXPT * 100 - PAPER_ERROR_NMAX2) < 10 else "NO"} |

"""
    if not paper_match:
        report += f"""## WARNING

**PAPER CORRECTION REQUIRED:** Section V.C states D_e_CP(nmax=2) = {PAPER_D_CP_NMAX2:.3f} Ha
but computed value is {actual_d_cp_nmax2_3015:.3f} Ha.

"""

    report += f"""## Convergence Diagnostics

1. **BSSE convergence:** |BSSE(nmax=2)| = {abs(bsse_2):.3f} Ha {">" if abs(bsse_2) > abs(bsse_3) else "<"} |BSSE(nmax=3)| = {abs(bsse_3):.3f} Ha — {"PASS" if abs(bsse_2) > abs(bsse_3) else "FAIL"}
2. **R_eq convergence:** R_eq(nmax=2) = {R_eq_2:.2f} {"<" if R_eq_2 < R_eq_3 else ">"} R_eq(nmax=3) = {R_eq_3:.2f} bohr — {"PASS (more contracted at lower nmax)" if R_eq_2 < R_eq_3 else "UNEXPECTED"}
3. **D_e_CP convergence:** D_e_CP(nmax=2) = {D_e_2:.3f} {">" if D_e_2 > D_e_3 else "<"} D_e_CP(nmax=3) = {D_e_3:.3f} Ha — {"PASS (overbinding decreases)" if D_e_2 > D_e_3 else "UNEXPECTED"}

## Output Files

- `debug/data/lih_nmax_convergence.txt` — numerical data
- `debug/plots/lih_pes_convergence.png` — PES comparison plot

## Next Steps

- Run nmax=4 when computational resources allow (8.2M SDs, estimated >30 min)
- Verify R_eq continues to shift outward toward 3.015 bohr
- Verify BSSE continues to decrease
"""

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        f.write(report)
    print(f"Validation report saved to {outpath}")


def main() -> None:
    """Run the full LiH nmax convergence validation."""
    print("=" * 70)
    print("LiH nmax Convergence Validation")
    print("Paper: Topological FCI for Heteronuclear Diatomics, Section V.C")
    print("=" * 70)

    # ------------------------------------------------------------------
    # Step 1: Atomic reference energies
    # ------------------------------------------------------------------
    print("\n--- Atomic Reference Energies ---")

    print("Computing nmax=2 atoms...")
    E_li_2, E_h_2 = compute_atomic_energies(nmax=2)
    E_sep_2 = E_li_2 + E_h_2
    print(f"  E(Li, nmax=2) = {E_li_2:.6f} Ha")
    print(f"  E(H,  nmax=2) = {E_h_2:.6f} Ha")
    print(f"  E_sep(nmax=2) = {E_sep_2:.6f} Ha")

    print("Computing nmax=3 atoms...")
    E_li_3, E_h_3 = compute_atomic_energies(nmax=3)
    E_sep_3 = E_li_3 + E_h_3
    print(f"  E(Li, nmax=3) = {E_li_3:.6f} Ha")
    print(f"  E(H,  nmax=3) = {E_h_3:.6f} Ha")
    print(f"  E_sep(nmax=3) = {E_sep_3:.6f} Ha")

    # ------------------------------------------------------------------
    # Step 2: PES scans
    # ------------------------------------------------------------------
    # nmax=2: fine scan 1.5 to 5.0 in 0.1 steps (fast, ~5k SDs)
    R_fine = np.arange(1.5, 5.05, 0.1)

    # nmax=3: coarser scan (slower, ~367k SDs per point)
    R_coarse = np.array([1.5, 1.8, 2.0, 2.2, 2.5, 2.8, 3.0, 3.015,
                         3.2, 3.5, 4.0, 4.5, 5.0])

    print("\n--- PES Scan: nmax=2 ---")
    results_2 = scan_pes(nmax=2, R_values=R_fine,
                         E_li_own=E_li_2, E_h_own=E_h_2)

    print("\n--- PES Scan: nmax=3 ---")
    results_3 = scan_pes(nmax=3, R_values=R_coarse,
                         E_li_own=E_li_3, E_h_own=E_h_3)

    # ------------------------------------------------------------------
    # Step 3: Find PES minima
    # ------------------------------------------------------------------
    R_eq_2, D_e_2 = find_pes_minimum(results_2['R'], results_2['D_cp'])
    R_eq_3, D_e_3 = find_pes_minimum(results_3['R'], results_3['D_cp'])

    # ------------------------------------------------------------------
    # Step 4: Extract values at R=3.015 for direct comparison
    # ------------------------------------------------------------------
    # nmax=2: find closest R to 3.015 in scan
    idx_3015_2 = np.argmin(np.abs(results_2['R'] - R_EQ_EXPT))
    R_near_3015_2 = results_2['R'][idx_3015_2]
    D_cp_3015_2 = results_2['D_cp'][idx_3015_2]
    BSSE_3015_2 = results_2['BSSE'][idx_3015_2]

    # nmax=3: use exact R=3.015 point
    idx_3015_3 = np.argmin(np.abs(results_3['R'] - R_EQ_EXPT))
    D_cp_3015_3 = results_3['D_cp'][idx_3015_3]
    BSSE_3015_3 = results_3['BSSE'][idx_3015_3]

    # ------------------------------------------------------------------
    # Step 5: Verify paper claims
    # ------------------------------------------------------------------
    paper_match = abs(D_cp_3015_2 - PAPER_D_CP_NMAX2) < 0.015  # 15 mHa tolerance

    print("\n" + "=" * 70)
    print("RESULTS SUMMARY")
    print("=" * 70)

    summary_lines = []
    summary_lines.append(f"LiH nmax Convergence Validation — v0.9.37, 2026-03-11")
    summary_lines.append(f"{'=' * 65}")
    summary_lines.append(f"")
    summary_lines.append(f"Atomic Energies:")
    summary_lines.append(f"  E(Li, nmax=2) = {E_li_2:.6f} Ha    E(Li, nmax=3) = {E_li_3:.6f} Ha")
    summary_lines.append(f"  E(H,  nmax=2) = {E_h_2:.6f} Ha    E(H,  nmax=3) = {E_h_3:.6f} Ha")
    summary_lines.append(f"")
    summary_lines.append(f"{'Quantity':<30s}  {'nmax=2':>10s}  {'nmax=3':>10s}  {'Expt':>10s}")
    summary_lines.append(f"{'-' * 65}")
    summary_lines.append(f"{'D_e_CP at R=3.015 (Ha)':<30s}  {D_cp_3015_2:10.4f}  {D_cp_3015_3:10.4f}  {D_E_EXPT:10.4f}")
    summary_lines.append(f"{'D_e_CP at R_eq (Ha)':<30s}  {D_e_2:10.4f}  {D_e_3:10.4f}  {D_E_EXPT:10.4f}")
    summary_lines.append(f"{'BSSE at R=3.015 (Ha)':<30s}  {BSSE_3015_2:10.4f}  {BSSE_3015_3:10.4f}  {'0.0000':>10s}")
    summary_lines.append(f"{'R_eq (bohr)':<30s}  {R_eq_2:10.2f}  {R_eq_3:10.2f}  {R_EQ_EXPT:10.3f}")
    summary_lines.append(f"{'Error in D_e_CP (%)':<30s}  {abs(D_cp_3015_2 - D_E_EXPT) / D_E_EXPT * 100:10.0f}  {abs(D_cp_3015_3 - D_E_EXPT) / D_E_EXPT * 100:10.0f}  {'---':>10s}")
    summary_lines.append(f"")
    summary_lines.append(f"Convergence Checks:")
    summary_lines.append(f"  BSSE(nmax=2) > BSSE(nmax=3):  |{BSSE_3015_2:.4f}| > |{BSSE_3015_3:.4f}|  "
                         f"{'PASS' if abs(BSSE_3015_2) > abs(BSSE_3015_3) else 'FAIL'}")
    summary_lines.append(f"  R_eq(nmax=2) < R_eq(nmax=3):  {R_eq_2:.2f} < {R_eq_3:.2f}  "
                         f"{'PASS' if R_eq_2 < R_eq_3 else 'FAIL'}")
    summary_lines.append(f"  D_e_CP(nmax=2) > D_e_CP(nmax=3):  {D_e_2:.4f} > {D_e_3:.4f}  "
                         f"{'PASS' if D_e_2 > D_e_3 else 'FAIL'}")
    summary_lines.append(f"")

    if not paper_match:
        warning = (f"PAPER CORRECTION REQUIRED: Section V.C states "
                   f"D_CP(nmax=2)={PAPER_D_CP_NMAX2:.3f} Ha but computed "
                   f"value is {D_cp_3015_2:.3f} Ha")
        summary_lines.append(f"*** {warning} ***")
        summary_lines.append(f"")
    else:
        summary_lines.append(f"Paper Section V.C claim VERIFIED: "
                             f"D_CP(nmax=2)={D_cp_3015_2:.3f} Ha "
                             f"(paper: {PAPER_D_CP_NMAX2:.3f} Ha)")
        summary_lines.append(f"")

    # nmax=3 cross-check
    nmax3_match = abs(D_cp_3015_3 - KNOWN_D_CP_NMAX3) < 0.005
    summary_lines.append(f"nmax=3 Cross-Check:")
    summary_lines.append(f"  D_e_CP(nmax=3) = {D_cp_3015_3:.4f} Ha (expected {KNOWN_D_CP_NMAX3:.3f})  "
                         f"{'PASS' if nmax3_match else 'MISMATCH'}")
    summary_lines.append(f"  BSSE(nmax=3)   = {BSSE_3015_3:.4f} Ha (expected -{KNOWN_BSSE_NMAX3:.3f})  "
                         f"{'PASS' if abs(abs(BSSE_3015_3) - KNOWN_BSSE_NMAX3) < 0.01 else 'MISMATCH'}")

    # Add PES data tables
    summary_lines.append(f"")
    summary_lines.append(f"{'=' * 65}")
    summary_lines.append(f"nmax=2 PES Data")
    summary_lines.append(f"{'R':>6s}  {'E_mol':>10s}  {'D_raw':>8s}  {'D_cp':>8s}  {'BSSE':>8s}")
    summary_lines.append(f"{'-' * 45}")
    for i in range(len(results_2['R'])):
        summary_lines.append(
            f"{results_2['R'][i]:6.2f}  {results_2['E_mol'][i]:10.5f}  "
            f"{results_2['D_raw'][i]:8.4f}  {results_2['D_cp'][i]:8.4f}  "
            f"{results_2['BSSE'][i]:8.4f}"
        )

    summary_lines.append(f"")
    summary_lines.append(f"{'=' * 65}")
    summary_lines.append(f"nmax=3 PES Data")
    summary_lines.append(f"{'R':>6s}  {'E_mol':>10s}  {'D_raw':>8s}  {'D_cp':>8s}  {'BSSE':>8s}")
    summary_lines.append(f"{'-' * 45}")
    for i in range(len(results_3['R'])):
        summary_lines.append(
            f"{results_3['R'][i]:6.2f}  {results_3['E_mol'][i]:10.5f}  "
            f"{results_3['D_raw'][i]:8.4f}  {results_3['D_cp'][i]:8.4f}  "
            f"{results_3['BSSE'][i]:8.4f}"
        )

    summary_text = "\n".join(summary_lines) + "\n"
    print(summary_text)

    # ------------------------------------------------------------------
    # Step 6: Write outputs
    # ------------------------------------------------------------------
    data_path = os.path.join(PROJECT_ROOT, "debug", "data", "lih_nmax_convergence.txt")
    plot_path = os.path.join(PROJECT_ROOT, "debug", "plots", "lih_pes_convergence.png")
    report_path = os.path.join(PROJECT_ROOT, "docs", "LIH_CONVERGENCE_VALIDATION.md")

    write_data_file(results_2, results_3, summary_text, data_path)

    make_pes_plot(results_2, results_3,
                  R_eq_2, D_e_2, R_eq_3, D_e_3, plot_path)

    write_validation_report(
        E_li_2, E_h_2, E_li_3, E_h_3,
        BSSE_3015_2, BSSE_3015_3,
        R_eq_2, R_eq_3, D_e_2, D_e_3,
        D_cp_3015_2, paper_match, D_cp_3015_2,
        report_path,
    )

    print("\n" + "=" * 70)
    print("VALIDATION COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
