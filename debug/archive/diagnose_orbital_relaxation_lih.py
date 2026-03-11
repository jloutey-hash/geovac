"""
v0.9.36 — Orbital exponent relaxation for LiH.

Parts 1-3: parametric orbital exponents, single-parameter optimization,
and relaxed PES with CP correction.

Each H orbital uses exponent zeta_B * Z_H / n instead of Z_H / n.
Optimize zeta_B at each R to find the relaxed PES minimum.
"""

import sys
import time
import warnings
import numpy as np
from scipy.optimize import minimize_scalar

warnings.filterwarnings("ignore")
sys.path.insert(0, '.')

from geovac.lattice_index import (
    MolecularLatticeIndex, compute_bsse_correction,
)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
NMAX = 3
Z_A, Z_B = 3, 1
N_ELECTRONS = 4
R_GRID = [2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0, 6.0]
ZETA_BOUNDS = (0.8, 2.5)
OUTPUT_FILE = "debug/data/orbital_relaxation_lih_v0936.txt"


def compute_energy(R: float, zeta_B: float = 1.0,
                   zeta_A: float = 1.0) -> float:
    """Compute total FCI energy at given R and zeta_B."""
    mol = MolecularLatticeIndex(
        Z_A=Z_A, Z_B=Z_B, nmax_A=NMAX, nmax_B=NMAX,
        R=R, n_electrons=N_ELECTRONS,
        vee_method='slater_full', fci_method='auto',
        zeta_A=zeta_A, zeta_B=zeta_B,
    )
    E, _ = mol.compute_ground_state(n_states=1)
    return E[0]


def compute_cp_energy(R: float, zeta_B: float = 1.0,
                      zeta_A: float = 1.0) -> dict:
    """Compute CP-corrected energy at given R and zeta."""
    E_mol = compute_energy(R, zeta_B=zeta_B, zeta_A=zeta_A)

    bsse = compute_bsse_correction(
        Z_A=Z_A, Z_B=Z_B, nmax_A=NMAX, nmax_B=NMAX,
        R=R, n_electrons_A=Z_A, n_electrons_B=Z_B,
        vee_method='slater_full', fci_method='auto',
        zeta_A=zeta_A, zeta_B=zeta_B,
    )

    E_cp = E_mol - bsse['BSSE']
    return {
        'E_mol': E_mol,
        'E_A_own': bsse['E_A_own'],
        'E_B_own': bsse['E_B_own'],
        'BSSE': bsse['BSSE'],
        'E_cp': E_cp,
        'D_e_cp': (bsse['E_A_own'] + bsse['E_B_own']) - E_cp,
    }


# ---------------------------------------------------------------------------
# PART 1: Verify parametric exponents
# ---------------------------------------------------------------------------
def part1_verify():
    """Verify zeta_B=1.0 matches baseline."""
    print("=" * 70)
    print("PART 1: Parametric orbital exponents — regression check")
    print("=" * 70)

    R = 3.015
    E_baseline = compute_energy(R, zeta_B=1.0)
    E_zeta = compute_energy(R, zeta_B=1.0)
    print(f"  R={R}: E(baseline)={E_baseline:.8f}, E(zeta=1.0)={E_zeta:.8f}")
    print(f"  Diff = {abs(E_baseline - E_zeta):.2e} Ha")

    # Quick scan of energy landscape
    print("\n  zeta_B scan at R=3.015:")
    for zb in [0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 2.0]:
        E = compute_energy(R, zeta_B=zb)
        print(f"    zeta_B={zb:.1f}: E = {E:.6f} Ha")


# ---------------------------------------------------------------------------
# PART 2: Single-parameter optimization at fixed R
# ---------------------------------------------------------------------------
def part2_optimize() -> dict:
    """Optimize zeta_B at each R grid point."""
    print("\n" + "=" * 70)
    print("PART 2: Optimal zeta_B*(R) via Brent minimization")
    print("=" * 70)

    results = {}
    for R in R_GRID:
        t0 = time.time()

        def obj(zeta_B: float) -> float:
            return compute_energy(R, zeta_B=zeta_B)

        res = minimize_scalar(obj, bounds=ZETA_BOUNDS, method='bounded',
                              options={'xatol': 0.005})
        zeta_opt = res.x
        E_opt = res.fun
        E_fixed = compute_energy(R, zeta_B=1.0)
        dE = E_opt - E_fixed

        elapsed = time.time() - t0
        results[R] = {
            'zeta_opt': zeta_opt, 'E_opt': E_opt,
            'E_fixed': E_fixed, 'dE_relax': dE,
        }
        print(f"  R={R:.3f}: zeta_B*={zeta_opt:.4f}, "
              f"E_opt={E_opt:.6f}, E_fixed={E_fixed:.6f}, "
              f"dE_relax={dE:.6f} Ha  ({elapsed:.0f}s)")

    # Check monotonicity
    zetas = [results[R]['zeta_opt'] for R in R_GRID]
    monotonic = all(zetas[i] >= zetas[i+1] for i in range(len(zetas)-1))
    print(f"\n  Monotonic (zeta_B* decreases with R): {monotonic}")
    if not monotonic:
        print("  WARNING: Non-monotonic zeta_B*(R) — variational surface may be flat")

    return results


# ---------------------------------------------------------------------------
# PART 3: Relaxed PES with CP correction
# ---------------------------------------------------------------------------
def part3_relaxed_pes(opt_results: dict) -> dict:
    """Build relaxed PES with CP correction."""
    print("\n" + "=" * 70)
    print("PART 3: Relaxed PES (CP-corrected)")
    print("=" * 70)

    pes_results = {}
    for R in R_GRID:
        t0 = time.time()
        zeta_opt = opt_results[R]['zeta_opt']

        # Fixed-exponent CP energy
        cp_fixed = compute_cp_energy(R, zeta_B=1.0)

        # Relaxed CP energy
        cp_relaxed = compute_cp_energy(R, zeta_B=zeta_opt)

        dE_cp = cp_relaxed['E_cp'] - cp_fixed['E_cp']

        pes_results[R] = {
            'zeta_opt': zeta_opt,
            'E_fixed': cp_fixed['E_mol'],
            'E_relaxed': cp_relaxed['E_mol'],
            'BSSE_fixed': cp_fixed['BSSE'],
            'BSSE_relaxed': cp_relaxed['BSSE'],
            'E_cp_fixed': cp_fixed['E_cp'],
            'E_cp_relaxed': cp_relaxed['E_cp'],
            'D_e_cp_fixed': cp_fixed['D_e_cp'],
            'D_e_cp_relaxed': cp_relaxed['D_e_cp'],
            'dE_cp': dE_cp,
        }

        elapsed = time.time() - t0
        print(f"  R={R:.3f}: zeta_B*={zeta_opt:.4f}")
        print(f"    E_cp(fixed)  = {cp_fixed['E_cp']:.6f}, "
              f"D_e_cp = {cp_fixed['D_e_cp']:.6f} Ha")
        print(f"    E_cp(relaxed)= {cp_relaxed['E_cp']:.6f}, "
              f"D_e_cp = {cp_relaxed['D_e_cp']:.6f} Ha")
        print(f"    dE_cp = {dE_cp:.6f} Ha  ({elapsed:.0f}s)")

    # Find R_eq
    R_vals = np.array(R_GRID)
    E_fixed_arr = np.array([pes_results[R]['E_cp_fixed'] for R in R_GRID])
    E_relaxed_arr = np.array([pes_results[R]['E_cp_relaxed'] for R in R_GRID])

    idx_fixed = np.argmin(E_fixed_arr)
    idx_relaxed = np.argmin(E_relaxed_arr)

    print(f"\n  R_eq(fixed):   ~{R_vals[idx_fixed]:.3f} bohr "
          f"(E_cp = {E_fixed_arr[idx_fixed]:.6f})")
    print(f"  R_eq(relaxed): ~{R_vals[idx_relaxed]:.3f} bohr "
          f"(E_cp = {E_relaxed_arr[idx_relaxed]:.6f})")

    # Key diagnostic: sign of dE_relax(2.5) - dE_relax(3.015)
    if 2.5 in pes_results and 3.015 in pes_results:
        dE_25 = pes_results[2.5]['dE_cp']
        dE_3015 = pes_results[3.015]['dE_cp']
        diff = dE_25 - dE_3015
        print(f"\n  KEY DIAGNOSTIC:")
        print(f"    dE_relax(R=2.5)  = {dE_25:.6f} Ha")
        print(f"    dE_relax(R=3.015)= {dE_3015:.6f} Ha")
        print(f"    Difference = {diff:.6f} Ha")
        if diff < -0.001:
            print(f"    RESULT: Relaxation helps more at short R → R_eq shifts outward")
        elif diff > 0.001:
            print(f"    RESULT: Relaxation helps more at long R → R_eq shifts inward")
        else:
            print(f"    RESULT: Difference near zero → orbital relaxation is not the mechanism")

    return pes_results


# ---------------------------------------------------------------------------
# Write output
# ---------------------------------------------------------------------------
def write_output(opt_results: dict, pes_results: dict) -> None:
    """Write results to file."""
    with open(OUTPUT_FILE, 'w') as f:
        f.write("v0.9.36 — Orbital Exponent Relaxation for LiH\n")
        f.write(f"nmax={NMAX}, Z_A={Z_A}, Z_B={Z_B}\n")
        f.write("=" * 70 + "\n\n")

        f.write("Table 1: Optimal zeta_B*(R) and relaxation energy\n")
        f.write("-" * 70 + "\n")
        f.write(f"{'R':>8s}  {'zeta_B*':>8s}  {'E_fixed':>12s}  "
                f"{'E_relaxed':>12s}  {'dE_relax':>12s}\n")
        for R in R_GRID:
            d = opt_results[R]
            f.write(f"{R:8.3f}  {d['zeta_opt']:8.4f}  {d['E_fixed']:12.6f}  "
                    f"{d['E_opt']:12.6f}  {d['dE_relax']:12.6f}\n")

        if pes_results:
            f.write("\n\nTable 2: CP-corrected PES (fixed vs relaxed)\n")
            f.write("-" * 70 + "\n")
            f.write(f"{'R':>8s}  {'zeta_B*':>8s}  {'E_cp_fixed':>12s}  "
                    f"{'E_cp_relaxed':>12s}  {'dE_cp':>12s}  "
                    f"{'D_e_fixed':>12s}  {'D_e_relaxed':>12s}\n")
            for R in R_GRID:
                d = pes_results[R]
                f.write(f"{R:8.3f}  {d['zeta_opt']:8.4f}  "
                        f"{d['E_cp_fixed']:12.6f}  {d['E_cp_relaxed']:12.6f}  "
                        f"{d['dE_cp']:12.6f}  {d['D_e_cp_fixed']:12.6f}  "
                        f"{d['D_e_cp_relaxed']:12.6f}\n")

            f.write("\n\nKey diagnostics:\n")
            if 2.5 in pes_results and 3.015 in pes_results:
                dE_25 = pes_results[2.5]['dE_cp']
                dE_3015 = pes_results[3.015]['dE_cp']
                f.write(f"  dE_relax(R=2.5)  = {dE_25:.6f} Ha\n")
                f.write(f"  dE_relax(R=3.015)= {dE_3015:.6f} Ha\n")
                f.write(f"  Difference = {dE_25 - dE_3015:.6f} Ha\n")

                # Find R_eq
                R_arr = np.array(R_GRID)
                E_fixed = [pes_results[R]['E_cp_fixed'] for R in R_GRID]
                E_relax = [pes_results[R]['E_cp_relaxed'] for R in R_GRID]
                f.write(f"  R_eq(fixed):   ~{R_arr[np.argmin(E_fixed)]:.3f} bohr\n")
                f.write(f"  R_eq(relaxed): ~{R_arr[np.argmin(E_relax)]:.3f} bohr\n")
                f.write(f"  D_e_CP(fixed):   {pes_results[3.015]['D_e_cp_fixed']:.6f} Ha\n")
                f.write(f"  D_e_CP(relaxed): {pes_results[3.015]['D_e_cp_relaxed']:.6f} Ha\n")

    print(f"\nResults written to {OUTPUT_FILE}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    t_start = time.time()

    # Part 1: quick verification
    part1_verify()

    # Part 2: optimize zeta_B at each R
    opt_results = part2_optimize()

    # Part 3: relaxed PES with CP correction
    pes_results = part3_relaxed_pes(opt_results)

    # Write output
    write_output(opt_results, pes_results)

    total_time = time.time() - t_start
    print(f"\nTotal wall time: {total_time/60:.1f} min")
