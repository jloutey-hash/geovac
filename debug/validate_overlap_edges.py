"""
Overlap-Induced Edge Validation: R-Dependent Graph Topology
=============================================================

Tests whether adding cross-atom edges proportional to orbital overlap
creates a universal repulsive wall that gives correct equilibrium geometry
for both H2 and LiH without fitted parameters.

Edge formulas tested:
  's2'   : A[a,b] = scale * S²(a,b)
  'fock' : A[a,b] = scale * S²(a,b) * (Z_A/n_a)² * (Z_B/n_b)²
  'abs'  : A[a,b] = scale * |S(a,b)|
  'kappa': A[a,b] = scale * S²(a,b) * |κ|   (κ = -1/16)

For each formula and scale, we compute the PES and extract:
  - R_min: bond length at PES minimum
  - k_fit: force constant at minimum
  - D_e: well depth (relative to separated atoms)

Success criterion: same formula + scale gives R_min within 30% of R_eq
for BOTH H2 and LiH.

Output:
  debug/data/overlap_edges_h2.txt
  debug/data/overlap_edges_lih.txt (if run)
  debug/OVERLAP_EDGES_RESULTS.md

Date: March 2026
Version: v0.9.38+
"""

import sys
import os
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.coupled_en_lattice import numerical_force_constant
from geovac.nuclear_lattice import DIATOMIC_CONSTANTS, HARTREE_TO_CM
import warnings
warnings.filterwarnings('ignore')


# ======================================================================
# Configuration
# ======================================================================

# Edge formulas and scales to test
EDGE_MODES = ['s2', 'fock', 'abs', 'kappa']
SCALES = [1.0, 5.0, 10.0, 25.0, 50.0]

# Baseline (no overlap edges)
BASELINE = ('none', 0.0)

# R scan: 15 points for smooth PES
R_FACTORS = [0.50, 0.60, 0.70, 0.80, 0.85, 0.90, 0.95, 1.00,
             1.05, 1.10, 1.20, 1.40, 1.60, 2.00, 2.50]


# ======================================================================
# Core computation
# ======================================================================

def compute_energy(
    Z_A: int, Z_B: int, nmax: int, n_electrons: int, R: float,
    overlap_edges: str = None, overlap_edge_scale: float = 1.0,
) -> float:
    """Compute molecular FCI energy at bond length R with overlap edges."""
    from geovac.lattice_index import MolecularLatticeIndex
    mol = MolecularLatticeIndex(
        Z_A=Z_A, Z_B=Z_B,
        nmax_A=nmax, nmax_B=nmax,
        R=R,
        n_electrons=n_electrons,
        vee_method='slater_full',
        fci_method='auto',
        cross_nuclear_method='exact',
        cross_atom_vee=True,
        overlap_edges=overlap_edges,
        overlap_edge_scale=overlap_edge_scale,
    )
    eigvals, _ = mol.compute_ground_state(n_states=1)
    return eigvals[0]


def compute_separated(Z_A: int, Z_B: int, n_elec: int, nmax: int) -> float:
    """Compute separated atom energy."""
    from geovac.lattice_index import LatticeIndex

    if Z_A == Z_B:
        n_each = n_elec // 2
        E = LatticeIndex(
            n_electrons=n_each, max_n=nmax, nuclear_charge=Z_A,
            vee_method='slater_full', h1_method='exact', fci_method='auto',
        ).compute_ground_state(n_states=1)[0][0]
        return 2 * E
    else:
        # LiH specific
        E_A = LatticeIndex(
            n_electrons=3, max_n=nmax, nuclear_charge=Z_A,
            vee_method='slater_full', h1_method='exact', fci_method='auto',
        ).compute_ground_state(n_states=1)[0][0]
        E_B = LatticeIndex(
            n_electrons=1, max_n=nmax, nuclear_charge=Z_B,
            vee_method='slater_full', h1_method='exact', fci_method='auto',
        ).compute_ground_state(n_states=1)[0][0]
        return E_A + E_B


def find_minimum(R_arr: np.ndarray, E_arr: np.ndarray) -> dict:
    """Find PES minimum via polynomial fit."""
    # Find approximate minimum
    idx_min = np.argmin(E_arr)

    # Use points near minimum for polynomial fit
    lo = max(0, idx_min - 3)
    hi = min(len(R_arr), idx_min + 4)
    R_fit = R_arr[lo:hi]
    E_fit = E_arr[lo:hi]

    if len(R_fit) < 3:
        return {'R_min': R_arr[idx_min], 'E_min': E_arr[idx_min],
                'k_fit': 0.0, 'has_min': False}

    # Quadratic fit
    coeffs = np.polyfit(R_fit, E_fit, min(4, len(R_fit) - 1))
    poly = np.poly1d(coeffs)
    poly_d = poly.deriv()

    # Find root of derivative (minimum)
    roots = np.roots(poly_d.coeffs)
    real_roots = roots[np.isreal(roots)].real
    in_range = real_roots[(real_roots > R_arr[0]) & (real_roots < R_arr[-1])]

    if len(in_range) == 0:
        return {'R_min': R_arr[idx_min], 'E_min': E_arr[idx_min],
                'k_fit': 0.0, 'has_min': idx_min > 0 and idx_min < len(R_arr) - 1}

    # Pick root closest to discrete minimum
    R_min = in_range[np.argmin(np.abs(in_range - R_arr[idx_min]))]
    E_min = float(poly(R_min))
    k_fit = float(poly.deriv(2)(R_min))

    return {
        'R_min': float(R_min),
        'E_min': E_min,
        'k_fit': k_fit,
        'has_min': k_fit > 0,
    }


# ======================================================================
# Per-molecule scan
# ======================================================================

def scan_molecule(
    name: str, Z_A: int, Z_B: int, n_elec: int, nmax: int,
) -> dict:
    """Run PES scan for all edge formulas and scales."""
    print(f"\n{'='*70}")
    print(f"  {name}: Z_A={Z_A}, Z_B={Z_B}, Ne={n_elec}, nmax={nmax}")
    print(f"{'='*70}")

    c = DIATOMIC_CONSTANTS[name]
    r_e = c['r_e']
    R_values = sorted(set([r_e * f for f in R_FACTORS]))
    R_arr = np.array(R_values)

    print(f"  R_eq (expt) = {r_e:.4f} bohr")
    print(f"  R scan: {len(R_values)} points in [{R_arr[0]:.3f}, {R_arr[-1]:.3f}]")

    # Separated atoms
    print(f"  Computing separated atoms...", end="", flush=True)
    t0 = time.time()
    E_sep = compute_separated(Z_A, Z_B, n_elec, nmax)
    print(f" E_sep = {E_sep:.6f} ({time.time()-t0:.1f}s)")

    # Build list of (mode, scale) configurations to test
    configs = [BASELINE]
    for mode in EDGE_MODES:
        for scale in SCALES:
            configs.append((mode, scale))

    results = {}

    for mode, scale in configs:
        label = f"{mode}({scale:.0f})" if mode != 'none' else 'baseline'
        oe = mode if mode != 'none' else None
        oes = scale

        print(f"\n  --- {label} ---")
        E_arr = np.zeros(len(R_values))

        for i, R in enumerate(R_values):
            t0 = time.time()
            E = compute_energy(Z_A, Z_B, nmax, n_elec, R,
                               overlap_edges=oe, overlap_edge_scale=oes)
            E_arr[i] = E
            dt = time.time() - t0
            De = E_sep - E
            if i % 5 == 0 or i == len(R_values) - 1:
                print(f"    [{i+1}/{len(R_values)}] R={R:.4f}: E={E:.6f}, "
                      f"D_e={De:.4f} ({dt:.1f}s)")

        # Analyze PES
        info = find_minimum(R_arr, E_arr)
        De_well = E_sep - info['E_min']

        print(f"    -> R_min={info['R_min']:.4f}, E_min={info['E_min']:.6f}, "
              f"D_e={De_well:.4f}, k={info['k_fit']:.4f}, "
              f"has_min={info['has_min']}")
        if info['has_min']:
            err_R = (info['R_min'] - r_e) / r_e * 100
            print(f"    -> R_min error: {err_R:+.1f}%")

        results[(mode, scale)] = {
            'R_arr': R_arr,
            'E_arr': E_arr,
            'R_min': info['R_min'],
            'E_min': info['E_min'],
            'k_fit': info['k_fit'],
            'D_e': De_well,
            'has_min': info['has_min'],
        }

    return {
        'name': name,
        'r_e': r_e,
        'E_sep': E_sep,
        'results': results,
    }


# ======================================================================
# Comparison
# ======================================================================

def print_summary(mol_data: dict) -> None:
    """Print summary table for one molecule."""
    name = mol_data['name']
    r_e = mol_data['r_e']
    print(f"\n  Summary for {name} (R_eq = {r_e:.4f}):")
    print(f"  {'Config':>16s}  {'R_min':>8s}  {'err%':>7s}  {'k_fit':>8s}  "
          f"{'D_e':>8s}  {'min?':>5s}")

    for (mode, scale), res in mol_data['results'].items():
        label = f"{mode}({scale:.0f})" if mode != 'none' else 'baseline'
        err = (res['R_min'] - r_e) / r_e * 100 if res['has_min'] else float('nan')
        min_str = "YES" if res['has_min'] else "no"
        print(f"  {label:>16s}  {res['R_min']:8.4f}  {err:+7.1f}  "
              f"{res['k_fit']:8.4f}  {res['D_e']:8.4f}  {min_str:>5s}")


def save_data(mol_data: dict, filename: str) -> None:
    """Save PES data to file."""
    outdir = os.path.join(os.path.dirname(__file__), 'data')
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, filename)

    with open(outfile, 'w') as f:
        f.write(f"# Overlap Edge PES Data: {mol_data['name']}\n")
        f.write(f"# R_eq = {mol_data['r_e']:.4f}, E_sep = {mol_data['E_sep']:.8f}\n")
        f.write(f"# Date: {time.strftime('%Y-%m-%d %H:%M')}\n")

        for (mode, scale), res in mol_data['results'].items():
            label = f"{mode}({scale:.0f})" if mode != 'none' else 'baseline'
            f.write(f"#\n# === {label} ===\n")
            f.write(f"# R_min={res['R_min']:.4f}, k_fit={res['k_fit']:.6f}, "
                    f"D_e={res['D_e']:.6f}\n")
            f.write("# R  E_total\n")
            for R, E in zip(res['R_arr'], res['E_arr']):
                f.write(f"{R:8.4f}  {E:14.8f}\n")

    print(f"\n  Saved to {outfile}")


# ======================================================================
# Main
# ======================================================================

def main() -> None:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--molecule', '-m', default='H2',
                        choices=['H2', 'LiH', 'both'],
                        help='Which molecule(s) to test')
    parser.add_argument('--modes', nargs='+', default=None,
                        help='Edge modes to test (default: all)')
    parser.add_argument('--scales', nargs='+', type=float, default=None,
                        help='Scales to test (default: 1,5,10,25,50)')
    args = parser.parse_args()

    # Override globals if specified
    global EDGE_MODES, SCALES
    if args.modes:
        EDGE_MODES = args.modes
    if args.scales:
        SCALES = args.scales

    print("=" * 70)
    print("  Overlap-Induced Edge Validation")
    print("  Testing R-dependent graph topology for universal equilibrium")
    print("=" * 70)

    t_start = time.time()
    molecules = []

    if args.molecule in ('H2', 'both'):
        h2 = scan_molecule('H2', 1, 1, 2, 3)
        print_summary(h2)
        save_data(h2, 'overlap_edges_h2.txt')
        molecules.append(h2)

    if args.molecule in ('LiH', 'both'):
        lih = scan_molecule('LiH', 3, 1, 4, 3)
        print_summary(lih)
        save_data(lih, 'overlap_edges_lih.txt')
        molecules.append(lih)

    if len(molecules) == 2:
        print(f"\n{'='*70}")
        print("  UNIVERSALITY CHECK")
        print(f"{'='*70}")
        print(f"\n  {'Config':>16s}  {'H2 R_min':>10s}  {'H2 err%':>8s}  "
              f"{'LiH R_min':>10s}  {'LiH err%':>8s}  {'Universal?':>10s}")
        h2_res = molecules[0]['results']
        lih_res = molecules[1]['results']
        for key in h2_res:
            if key not in lih_res:
                continue
            mode, scale = key
            label = f"{mode}({scale:.0f})" if mode != 'none' else 'baseline'
            h2_r = h2_res[key]
            lih_r = lih_res[key]
            h2_err = (h2_r['R_min'] - 1.401) / 1.401 * 100
            lih_err = (lih_r['R_min'] - 3.015) / 3.015 * 100
            both_min = h2_r['has_min'] and lih_r['has_min']
            both_close = abs(h2_err) < 30 and abs(lih_err) < 30 and both_min
            uni = "YES" if both_close else "no"
            print(f"  {label:>16s}  {h2_r['R_min']:10.4f}  {h2_err:+8.1f}  "
                  f"{lih_r['R_min']:10.4f}  {lih_err:+8.1f}  {uni:>10s}")

    dt = time.time() - t_start
    print(f"\n  Total time: {dt:.0f}s ({dt/60:.1f} min)")
    print("  Done.")


if __name__ == '__main__':
    main()
