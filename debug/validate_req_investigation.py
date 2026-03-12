"""
Investigate the LiH R_eq discrepancy: raw vs CP-corrected PES under s_only
cross-atom V_ee (the v0.9.11 baseline configuration).

Three questions answered:
  Q1: Reproduce R_eq ≈ 2.5 bohr and D_CP = 0.093 Ha (nmax=3, s_only)
  Q2: Is R_eq = 2.5 the raw or CP-corrected minimum?
  Q3: nmax=2 s_only scan for apples-to-apples convergence comparison

Output:
  stdout                                       — consolidated results table
  debug/data/lih_req_investigation.txt         — same table as text file
  debug/plots/lih_raw_vs_cp_pes.png           — nmax=3 raw vs CP PES
  debug/plots/lih_req_investigation.png        — three-panel figure
  docs/LIH_REQ_INVESTIGATION.md               — validation report

Date: 2026-03-11
Version: v0.9.37
"""
import os
import sys
import warnings
from typing import Dict, Tuple

import numpy as np

warnings.filterwarnings('ignore')

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from geovac.lattice_index import (
    LatticeIndex,
    MolecularLatticeIndex,
)

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
D_E_EXPT = 0.092       # Ha, Huber & Herzberg 1979
R_EQ_EXPT = 3.015      # bohr

# Paper claims (v0.9.11 baseline, s_only cross-atom V_ee)
PAPER_D_CP_NMAX3 = 0.093   # Ha
PAPER_REQ_NMAX3 = 2.5      # bohr (approximate)
PAPER_BSSE_NMAX3 = 0.115   # Ha (magnitude)
PAPER_D_CP_NMAX2 = 0.270   # Ha
PAPER_ERROR_NMAX2 = 192     # percent

# Integral configuration: v0.9.11 baseline
VEE_METHOD = 'slater_full'
FCI_METHOD = 'auto'
CROSS_ATOM_VEE = 's_only'   # key: matches v0.9.11 before v0.9.35 changed default

# Discrepancy threshold
DISCREPANCY_THRESHOLD = 0.005  # 5 mHa


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
        vee_method=VEE_METHOD, h1_method='exact', fci_method=FCI_METHOD,
    )
    E_li = li.compute_ground_state(n_states=1)[0][0]

    h = LatticeIndex(
        n_electrons=1, max_n=nmax, nuclear_charge=1,
        vee_method=VEE_METHOD, h1_method='exact', fci_method=FCI_METHOD,
    )
    E_h = h.compute_ground_state(n_states=1)[0][0]

    return E_li, E_h


def compute_pes_point(
    nmax: int,
    R: float,
    E_li_own: float,
    E_h_own: float,
) -> Dict[str, float]:
    """Compute raw and CP-corrected binding energy at a single R.

    Parameters
    ----------
    nmax : int
        Maximum principal quantum number.
    R : float
        Internuclear distance (bohr).
    E_li_own : float
        Isolated Li energy (Ha).
    E_h_own : float
        Isolated H energy (Ha).

    Returns
    -------
    result : dict
        Keys: E_mol, E_li_ghost, E_h_ghost, E_ghost, D_raw, D_cp, BSSE.
    """
    E_sep_own = E_li_own + E_h_own

    # Molecular energy
    mol = MolecularLatticeIndex(
        Z_A=3, Z_B=1, nmax_A=nmax, nmax_B=nmax,
        R=R, n_electrons=4,
        vee_method=VEE_METHOD, fci_method=FCI_METHOD,
        cross_atom_vee=CROSS_ATOM_VEE,
    )
    E_mol = mol.compute_ground_state(n_states=1)[0][0]

    # Ghost-fragment energies for CP correction
    li_g = MolecularLatticeIndex(
        Z_A=3, Z_B=0, nmax_A=nmax, nmax_B=nmax,
        R=R, n_electrons=3,
        vee_method=VEE_METHOD, fci_method=FCI_METHOD,
        cross_atom_vee=CROSS_ATOM_VEE,
    )
    E_li_g = li_g.compute_ground_state(n_states=1)[0][0]

    h_g = MolecularLatticeIndex(
        Z_A=0, Z_B=1, nmax_A=nmax, nmax_B=nmax,
        R=R, n_electrons=1,
        vee_method=VEE_METHOD, fci_method=FCI_METHOD,
        cross_atom_vee=CROSS_ATOM_VEE,
    )
    E_h_g = h_g.compute_ground_state(n_states=1)[0][0]

    E_ghost = E_li_g + E_h_g
    D_raw = E_sep_own - E_mol
    D_cp = E_ghost - E_mol
    BSSE = E_ghost - E_sep_own

    return {
        'E_mol': E_mol,
        'E_li_ghost': E_li_g,
        'E_h_ghost': E_h_g,
        'E_ghost': E_ghost,
        'D_raw': D_raw,
        'D_cp': D_cp,
        'BSSE': BSSE,
    }


def scan_pes(
    nmax: int,
    R_values: np.ndarray,
    E_li_own: float,
    E_h_own: float,
) -> Dict[str, np.ndarray]:
    """Scan the PES over an array of R values.

    Parameters
    ----------
    nmax : int
        Maximum principal quantum number.
    R_values : np.ndarray
        Internuclear distances (bohr).
    E_li_own : float
        Isolated Li energy (Ha).
    E_h_own : float
        Isolated H energy (Ha).

    Returns
    -------
    results : dict
        Keys: R, E_mol, D_raw, D_cp, BSSE. Each np.ndarray.
    """
    n_pts = len(R_values)
    E_mol = np.zeros(n_pts)
    D_raw = np.zeros(n_pts)
    D_cp = np.zeros(n_pts)
    BSSE = np.zeros(n_pts)

    for i, R in enumerate(R_values):
        print(f"  nmax={nmax}, R={R:.2f} bohr ({i+1}/{n_pts})...",
              end="", flush=True)
        pt = compute_pes_point(nmax, R, E_li_own, E_h_own)
        E_mol[i] = pt['E_mol']
        D_raw[i] = pt['D_raw']
        D_cp[i] = pt['D_cp']
        BSSE[i] = pt['BSSE']
        print(f" E={E_mol[i]:.5f}, D_raw={D_raw[i]:.4f}, D_cp={D_cp[i]:.4f}")

    return {
        'R': R_values,
        'E_mol': E_mol,
        'D_raw': D_raw,
        'D_cp': D_cp,
        'BSSE': BSSE,
    }


def find_pes_minimum(
    R_values: np.ndarray,
    energies: np.ndarray,
    mode: str = 'max',
) -> Tuple[float, float]:
    """Find PES extremum via parabolic interpolation.

    For binding energies (D_raw, D_cp), the equilibrium is the MAXIMUM.
    For total energies (E_mol), the equilibrium is the MINIMUM.

    Parameters
    ----------
    R_values : np.ndarray
        Internuclear distances (bohr).
    energies : np.ndarray
        Energy values to find extremum of.
    mode : str
        'max' to find maximum (for D_raw, D_cp), 'min' for minimum (E_mol).

    Returns
    -------
    R_eq : float
        Equilibrium bond length (bohr).
    E_eq : float
        Energy at equilibrium.
    """
    if mode == 'max':
        idx = np.argmax(energies)
    else:
        idx = np.argmin(energies)

    # Boundary case — extremum at edge of scan
    if idx == 0 or idx == len(R_values) - 1:
        return R_values[idx], energies[idx]

    # Parabolic interpolation through 3 points
    R_m, R_0, R_p = R_values[idx - 1], R_values[idx], R_values[idx + 1]
    E_m, E_0, E_p = energies[idx - 1], energies[idx], energies[idx + 1]

    denom = 2.0 * ((R_m - R_0) * (E_m - E_p) - (R_m - R_p) * (E_m - E_0))
    if abs(denom) < 1e-15:
        return R_0, E_0

    R_eq = R_m - (
        (R_m - R_0)**2 * (E_m - E_p) - (R_m - R_p)**2 * (E_m - E_0)
    ) / denom

    # Quadratic fit for energy at R_eq
    A = np.array([
        [R_m**2, R_m, 1.0],
        [R_0**2, R_0, 1.0],
        [R_p**2, R_p, 1.0],
    ])
    coeffs = np.linalg.solve(A, np.array([E_m, E_0, E_p]))
    E_eq = coeffs[0] * R_eq**2 + coeffs[1] * R_eq + coeffs[2]

    return R_eq, E_eq


def make_three_panel_figure(
    res2: Dict[str, np.ndarray],
    res3: Dict[str, np.ndarray],
    req_raw_2: float,
    req_cp_2: float,
    req_raw_3: float,
    req_cp_3: float,
    outpath: str,
) -> None:
    """Generate three-panel PES comparison figure.

    Panel 1: nmax=2 raw vs CP
    Panel 2: nmax=3 raw vs CP
    Panel 3: CP-corrected for both nmax overlaid

    Parameters
    ----------
    res2, res3 : dict
        PES scan results for nmax=2 and nmax=3.
    req_raw_2, req_cp_2, req_raw_3, req_cp_3 : float
        Equilibrium bond lengths.
    outpath : str
        Output file path.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not available, skipping plots.")
        return

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Panel 1: nmax=2 raw vs CP
    ax = axes[0]
    ax.plot(res2['R'], res2['D_raw'], 'b-', label='D_raw', linewidth=1.5)
    ax.plot(res2['R'], res2['D_cp'], 'r-', label='D_cp', linewidth=1.5)
    ax.axvline(req_raw_2, color='blue', linestyle='--', alpha=0.5,
               label=f'R_eq(raw)={req_raw_2:.2f}')
    ax.axvline(req_cp_2, color='red', linestyle='--', alpha=0.5,
               label=f'R_eq(CP)={req_cp_2:.2f}')
    ax.axvline(R_EQ_EXPT, color='green', linestyle=':', alpha=0.4)
    ax.axhline(D_E_EXPT, color='green', linestyle=':', alpha=0.4)
    ax.set_xlabel('R (bohr)')
    ax.set_ylabel('D_e (Ha)')
    ax.set_title('nmax=2, s_only')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # Panel 2: nmax=3 raw vs CP
    ax = axes[1]
    ax.plot(res3['R'], res3['D_raw'], 'b-', label='D_raw', linewidth=1.5)
    ax.plot(res3['R'], res3['D_cp'], 'r-', label='D_cp', linewidth=1.5)
    ax.axvline(req_raw_3, color='blue', linestyle='--', alpha=0.5,
               label=f'R_eq(raw)={req_raw_3:.2f}')
    ax.axvline(req_cp_3, color='red', linestyle='--', alpha=0.5,
               label=f'R_eq(CP)={req_cp_3:.2f}')
    ax.axvline(R_EQ_EXPT, color='green', linestyle=':', alpha=0.4)
    ax.axhline(D_E_EXPT, color='green', linestyle=':', alpha=0.4)
    ax.set_xlabel('R (bohr)')
    ax.set_ylabel('D_e (Ha)')
    ax.set_title('nmax=3, s_only')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # Panel 3: CP-corrected overlay
    ax = axes[2]
    ax.plot(res2['R'], res2['D_cp'], 'r-o', markersize=2,
            label=f'nmax=2 CP (R_eq={req_cp_2:.2f})')
    ax.plot(res3['R'], res3['D_cp'], 'b-s', markersize=2,
            label=f'nmax=3 CP (R_eq={req_cp_3:.2f})')
    ax.axhline(D_E_EXPT, color='green', linestyle='--', alpha=0.7,
               label=f'Expt D_e={D_E_EXPT}')
    ax.axvline(R_EQ_EXPT, color='green', linestyle=':', alpha=0.4,
               label=f'Expt R_eq={R_EQ_EXPT}')
    ax.set_xlabel('R (bohr)')
    ax.set_ylabel('D_e^CP (Ha)')
    ax.set_title('CP-corrected: nmax convergence')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Three-panel figure saved to {outpath}")


def make_raw_vs_cp_figure(
    res3: Dict[str, np.ndarray],
    req_raw_3: float,
    req_cp_3: float,
    outpath: str,
) -> None:
    """Generate nmax=3 raw vs CP PES comparison (Question 2).

    Parameters
    ----------
    res3 : dict
        PES scan results for nmax=3.
    req_raw_3, req_cp_3 : float
        Raw and CP-corrected equilibrium bond lengths.
    outpath : str
        Output file path.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not available, skipping raw_vs_cp plot.")
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(res3['R'], res3['D_raw'], 'b-', linewidth=2,
            label=f'D_raw (R_eq={req_raw_3:.2f} bohr)')
    ax.plot(res3['R'], res3['D_cp'], 'r-', linewidth=2,
            label=f'D_cp  (R_eq={req_cp_3:.2f} bohr)')
    ax.axvline(req_raw_3, color='blue', linestyle='--', alpha=0.6)
    ax.axvline(req_cp_3, color='red', linestyle='--', alpha=0.6)
    ax.axvline(R_EQ_EXPT, color='green', linestyle=':', alpha=0.5,
               label=f'Expt R_eq={R_EQ_EXPT}')
    ax.axhline(D_E_EXPT, color='green', linestyle='--', alpha=0.5,
               label=f'Expt D_e={D_E_EXPT}')
    ax.set_xlabel('R (bohr)', fontsize=12)
    ax.set_ylabel('D_e (Ha)', fontsize=12)
    ax.set_title('LiH nmax=3, s_only: Raw vs CP-Corrected PES', fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Raw vs CP figure saved to {outpath}")


def interpolate_at_R(
    R_values: np.ndarray,
    values: np.ndarray,
    R_target: float,
) -> float:
    """Linearly interpolate a quantity at a specific R value.

    Parameters
    ----------
    R_values : np.ndarray
        Sorted R grid.
    values : np.ndarray
        Corresponding values.
    R_target : float
        Target R to interpolate at.

    Returns
    -------
    value : float
        Interpolated value.
    """
    return float(np.interp(R_target, R_values, values))


def main() -> None:
    """Run the full R_eq investigation."""
    print("=" * 70)
    print("LiH R_eq Investigation: Raw vs CP PES under s_only V_ee")
    print(f"Configuration: cross_atom_vee='{CROSS_ATOM_VEE}' (v0.9.11 baseline)")
    print("=" * 70)

    # ------------------------------------------------------------------
    # Atomic reference energies (same for both nmax)
    # ------------------------------------------------------------------
    print("\n--- Atomic Reference Energies ---")

    print("Computing nmax=2 atoms...")
    E_li_2, E_h_2 = compute_atomic_energies(nmax=2)
    print(f"  E(Li, nmax=2) = {E_li_2:.6f} Ha")
    print(f"  E(H,  nmax=2) = {E_h_2:.6f} Ha")

    print("Computing nmax=3 atoms...")
    E_li_3, E_h_3 = compute_atomic_energies(nmax=3)
    print(f"  E(Li, nmax=3) = {E_li_3:.6f} Ha")
    print(f"  E(H,  nmax=3) = {E_h_3:.6f} Ha")

    # ------------------------------------------------------------------
    # PES scans: R = 1.5 to 4.0 in 0.05 steps
    # ------------------------------------------------------------------
    R_scan = np.arange(1.50, 4.025, 0.05)

    print(f"\n--- PES Scan: nmax=3, s_only ({len(R_scan)} points) ---")
    res3 = scan_pes(nmax=3, R_values=R_scan, E_li_own=E_li_3, E_h_own=E_h_3)

    print(f"\n--- PES Scan: nmax=2, s_only ({len(R_scan)} points) ---")
    res2 = scan_pes(nmax=2, R_values=R_scan, E_li_own=E_li_2, E_h_own=E_h_2)

    # ------------------------------------------------------------------
    # Find equilibria: raw (max of D_raw) and CP (max of D_cp)
    # ------------------------------------------------------------------
    req_raw_3, De_raw_3 = find_pes_minimum(res3['R'], res3['D_raw'], mode='max')
    req_cp_3, De_cp_3 = find_pes_minimum(res3['R'], res3['D_cp'], mode='max')

    req_raw_2, De_raw_2 = find_pes_minimum(res2['R'], res2['D_raw'], mode='max')
    req_cp_2, De_cp_2 = find_pes_minimum(res2['R'], res2['D_cp'], mode='max')

    # Also check E_mol minimum (should agree with D_raw maximum)
    req_Emol_3, _ = find_pes_minimum(res3['R'], res3['E_mol'], mode='min')
    req_Emol_2, _ = find_pes_minimum(res2['R'], res2['E_mol'], mode='min')

    # ------------------------------------------------------------------
    # Interpolate at R=3.015 for direct paper comparison
    # ------------------------------------------------------------------
    D_cp_at_3015_3 = interpolate_at_R(res3['R'], res3['D_cp'], R_EQ_EXPT)
    D_cp_at_3015_2 = interpolate_at_R(res2['R'], res2['D_cp'], R_EQ_EXPT)
    BSSE_at_3015_3 = interpolate_at_R(res3['R'], res3['BSSE'], R_EQ_EXPT)
    BSSE_at_3015_2 = interpolate_at_R(res2['R'], res2['BSSE'], R_EQ_EXPT)
    D_raw_at_3015_3 = interpolate_at_R(res3['R'], res3['D_raw'], R_EQ_EXPT)
    D_raw_at_3015_2 = interpolate_at_R(res2['R'], res2['D_raw'], R_EQ_EXPT)

    # ------------------------------------------------------------------
    # Consolidated results table
    # ------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("CONSOLIDATED RESULTS (cross_atom_vee='s_only')")
    print("=" * 78)

    rows = [
        ("R_eq raw (bohr)",        req_raw_2,      req_raw_3,      "~2.0",         None),
        ("R_eq CP (bohr)",         req_cp_2,       req_cp_3,       f"~{PAPER_REQ_NMAX3}",  None),
        ("R_eq E_mol min (bohr)",  req_Emol_2,     req_Emol_3,     "~2.5",         None),
        ("D_raw at R_eq (Ha)",     De_raw_2,       De_raw_3,       "0.205",        None),
        ("D_cp at R_eq (Ha)",      De_cp_2,        De_cp_3,        f"{PAPER_D_CP_NMAX3}",  PAPER_D_CP_NMAX3),
        ("D_cp at R=3.015 (Ha)",   D_cp_at_3015_2, D_cp_at_3015_3, f"{PAPER_D_CP_NMAX3}",  PAPER_D_CP_NMAX3),
        ("D_raw at R=3.015 (Ha)",  D_raw_at_3015_2, D_raw_at_3015_3, "0.205",     0.205),
        ("BSSE at R=3.015 (Ha)",   BSSE_at_3015_2, BSSE_at_3015_3, f"-{PAPER_BSSE_NMAX3}", -PAPER_BSSE_NMAX3),
        ("D_cp(nmax=2) @3.015",    D_cp_at_3015_2, None,           f"{PAPER_D_CP_NMAX2}",  PAPER_D_CP_NMAX2),
    ]

    header = f"{'Quantity':<28s}  {'nmax=2 s_only':>14s}  {'nmax=3 s_only':>14s}  {'Paper':>14s}"
    sep = "-" * len(header)
    lines = [header, sep]

    warnings_list = []

    for label, val2, val3, paper_str, paper_num in rows:
        s2 = f"{val2:14.4f}" if val2 is not None else f"{'---':>14s}"
        s3 = f"{val3:14.4f}" if val3 is not None else f"{'---':>14s}"
        sp = f"{paper_str:>14s}"
        line = f"{label:<28s}  {s2}  {s3}  {sp}"
        lines.append(line)

        # Check discrepancy for nmax=3 against paper
        if paper_num is not None and val3 is not None:
            diff = abs(val3 - paper_num)
            if diff > DISCREPANCY_THRESHOLD:
                w = (f"  WARNING: {label} nmax=3 = {val3:.4f}, "
                     f"paper = {paper_num:.4f}, diff = {diff:.4f} Ha "
                     f"(> {DISCREPANCY_THRESHOLD} Ha threshold)")
                warnings_list.append(w)

        # Check nmax=2 D_cp against paper
        if paper_num is not None and val3 is None and val2 is not None:
            diff = abs(val2 - paper_num)
            if diff > DISCREPANCY_THRESHOLD:
                w = (f"  WARNING: {label} nmax=2 = {val2:.4f}, "
                     f"paper = {paper_num:.4f}, diff = {diff:.4f} Ha "
                     f"(> {DISCREPANCY_THRESHOLD} Ha threshold)")
                warnings_list.append(w)

    lines.append("")

    # Key diagnostic answer
    lines.append("KEY FINDING: Is R_eq=2.5 the raw or CP-corrected minimum?")
    lines.append(f"  R_eq(raw, nmax=3)    = {req_raw_3:.2f} bohr")
    lines.append(f"  R_eq(CP, nmax=3)     = {req_cp_3:.2f} bohr")
    lines.append(f"  R_eq(E_mol, nmax=3)  = {req_Emol_3:.2f} bohr")
    if abs(req_raw_3 - req_cp_3) < 0.01:
        lines.append("  => Raw and CP minima coincide (BSSE is R-independent).")
    elif req_cp_3 < req_raw_3 - 0.1:
        lines.append("  => CP minimum is MORE contracted than raw minimum.")
    elif req_cp_3 > req_raw_3 + 0.1:
        lines.append("  => CP minimum is LESS contracted than raw minimum.")
    else:
        lines.append("  => Raw and CP minima are within 0.1 bohr of each other.")

    # BSSE R-dependence check
    bsse_std = np.std(res3['BSSE'])
    bsse_range = np.max(res3['BSSE']) - np.min(res3['BSSE'])
    lines.append(f"\n  BSSE R-dependence (nmax=3): std={bsse_std:.6f}, "
                 f"range={bsse_range:.6f} Ha")
    if bsse_range < 0.001:
        lines.append("  => BSSE is effectively R-independent. Raw and CP PES "
                     "have identical shape.")
    else:
        lines.append(f"  => BSSE varies by {bsse_range:.4f} Ha across R range.")

    lines.append("")
    for w in warnings_list:
        lines.append(w)

    summary = "\n".join(lines)
    print(summary)

    # ------------------------------------------------------------------
    # Save data file
    # ------------------------------------------------------------------
    data_path = os.path.join(PROJECT_ROOT, "debug", "data",
                             "lih_req_investigation.txt")
    os.makedirs(os.path.dirname(data_path), exist_ok=True)

    with open(data_path, 'w') as f:
        f.write(f"LiH R_eq Investigation — v0.9.37, 2026-03-11\n")
        f.write(f"Config: cross_atom_vee='{CROSS_ATOM_VEE}'\n\n")
        f.write(summary + "\n\n")

        for label, nmax, res in [("nmax=2", 2, res2), ("nmax=3", 3, res3)]:
            f.write(f"{'=' * 60}\n{label} PES Data (s_only)\n")
            f.write(f"{'R':>6s}  {'E_mol':>10s}  {'D_raw':>8s}  "
                    f"{'D_cp':>8s}  {'BSSE':>8s}\n")
            f.write(f"{'-' * 48}\n")
            for i in range(len(res['R'])):
                f.write(f"{res['R'][i]:6.2f}  {res['E_mol'][i]:10.5f}  "
                        f"{res['D_raw'][i]:8.4f}  {res['D_cp'][i]:8.4f}  "
                        f"{res['BSSE'][i]:8.4f}\n")
            f.write("\n")

    print(f"\nData saved to {data_path}")

    # ------------------------------------------------------------------
    # Plots
    # ------------------------------------------------------------------
    plot1_path = os.path.join(PROJECT_ROOT, "debug", "plots",
                              "lih_raw_vs_cp_pes.png")
    make_raw_vs_cp_figure(res3, req_raw_3, req_cp_3, plot1_path)

    plot3_path = os.path.join(PROJECT_ROOT, "debug", "plots",
                              "lih_req_investigation.png")
    make_three_panel_figure(res2, res3,
                            req_raw_2, req_cp_2,
                            req_raw_3, req_cp_3,
                            plot3_path)

    # ------------------------------------------------------------------
    # Validation report
    # ------------------------------------------------------------------
    report_path = os.path.join(PROJECT_ROOT, "docs",
                               "LIH_REQ_INVESTIGATION.md")

    # Determine status
    req3_match = abs(req_raw_3 - PAPER_REQ_NMAX3) < 0.3 or \
                 abs(req_cp_3 - PAPER_REQ_NMAX3) < 0.3
    dcp3_match = abs(D_cp_at_3015_3 - PAPER_D_CP_NMAX3) < DISCREPANCY_THRESHOLD
    status = "PASS" if (req3_match and dcp3_match) else "PARTIAL — see findings"

    report = f"""# LiH R_eq Investigation Report

**Date:** 2026-03-11
**Version:** v0.9.37
**Status:** {status}
**Script:** `debug/validate_req_investigation.py`
**Configuration:** `cross_atom_vee='{CROSS_ATOM_VEE}'` (v0.9.11 baseline)

## Purpose

Determine whether R_eq = 2.5 bohr is the raw PES minimum, the CP-corrected
PES minimum, or neither, using the v0.9.11-era `s_only` cross-atom V_ee
configuration that produced the paper's benchmark values.

## Is R_eq = 2.5 bohr the raw or CP-corrected minimum?

| PES type | R_eq (nmax=3) | D_e at R_eq (Ha) |
|----------|---------------|------------------|
| Raw      | {req_raw_3:.2f} bohr | {De_raw_3:.4f} |
| CP-corrected | {req_cp_3:.2f} bohr | {De_cp_3:.4f} |
| E_mol minimum | {req_Emol_3:.2f} bohr | --- |
| Experiment | {R_EQ_EXPT} bohr | {D_E_EXPT} |

BSSE R-dependence: std = {bsse_std:.6f} Ha, range = {bsse_range:.6f} Ha over
R = {res3['R'][0]:.1f} to {res3['R'][-1]:.1f} bohr.

{"BSSE is effectively R-independent. The raw and CP-corrected PES curves are parallel (identical shape, constant vertical offset). Therefore R_eq(raw) = R_eq(CP)." if bsse_range < 0.001 else f"BSSE varies by {bsse_range:.4f} Ha, so raw and CP surfaces have different shapes."}

## Convergence Comparison (s_only, apples-to-apples)

| Quantity | nmax=2 | nmax=3 | Paper claim | Match? |
|----------|--------|--------|-------------|--------|
| R_eq raw (bohr) | {req_raw_2:.2f} | {req_raw_3:.2f} | ~2.0/~2.5 | --- |
| R_eq CP (bohr) | {req_cp_2:.2f} | {req_cp_3:.2f} | ~2.5 | {"YES" if abs(req_cp_3 - PAPER_REQ_NMAX3) < 0.3 else "NO"} |
| D_cp at R=3.015 (Ha) | {D_cp_at_3015_2:.4f} | {D_cp_at_3015_3:.4f} | 0.093 | {"YES" if dcp3_match else "NO"} |
| BSSE at R=3.015 (Ha) | {BSSE_at_3015_2:.4f} | {BSSE_at_3015_3:.4f} | -0.115 | {"YES" if abs(BSSE_at_3015_3 + PAPER_BSSE_NMAX3) < DISCREPANCY_THRESHOLD else "NO"} |
| D_cp(nmax=2) at R=3.015 | {D_cp_at_3015_2:.4f} | --- | 0.270 | {"YES" if abs(D_cp_at_3015_2 - PAPER_D_CP_NMAX2) < DISCREPANCY_THRESHOLD * 3 else "NO"} |

## Discrepancies

"""
    if warnings_list:
        for w in warnings_list:
            report += f"- {w.strip()}\n"
    else:
        report += "None detected (all values within 5 mHa of paper claims).\n"

    report += f"""
## Output Files

- `debug/data/lih_req_investigation.txt` — full PES data
- `debug/plots/lih_raw_vs_cp_pes.png` — nmax=3 raw vs CP comparison
- `debug/plots/lih_req_investigation.png` — three-panel convergence figure

## Next Steps

- If R_eq discrepancy persists, extend scan below R=1.5 bohr
- Verify paper's R_eq=2.5 claim against original v0.9.11 git checkout
- Run nmax=4 when computational budget allows
"""

    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    with open(report_path, 'w') as f:
        f.write(report)
    print(f"Report saved to {report_path}")

    print("\n" + "=" * 70)
    print("INVESTIGATION COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
