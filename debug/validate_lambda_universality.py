"""
Lambda Universality Test Across Diatomic Molecules
====================================================

Tests whether lambda_derived ≈ 0.14 is universal across molecules.
If lambda is the same for H2 and LiH (and eventually HCl, CO), it's a
topological constant of the GeoVac framework — not a per-molecule fit.

Method (for each molecule):
1. Build MolecularLatticeIndex at ~12 R values spanning R_eq ± 50%
2. Compute E_elec(R) for lambda=0 (bare) and lambda=0.02, 0.05, 0.10
3. Fit polynomial to get force constant k_bare, k_correction
4. Derive lambda = (k_Morse - k_bare) / k_correction_per_lambda

Molecules tested:
  H2:  Z_A=1, Z_B=1, Ne=2, nmax=3 → C(20,2) = 190 SDs (fast)
  LiH: Z_A=3, Z_B=1, Ne=4, nmax=3 → C(56,4) = 367k SDs (slower)

Output:
  debug/data/lambda_universality.txt
  debug/LAMBDA_UNIVERSALITY.md (manual analysis)

Date: March 2026
Version: v0.9.38+
"""

import sys
import os
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.coupled_en_lattice import (
    compute_morse_parameters,
    numerical_force_constant,
    derive_lambda_from_force_constant,
)
from geovac.nuclear_lattice import DIATOMIC_CONSTANTS, HARTREE_TO_CM


# ======================================================================
# Configuration
# ======================================================================

LAMBDA_VALUES = [0.0, 0.02, 0.05, 0.10]

# Molecule configurations: (name, Z_A, Z_B, n_electrons, nmax)
MOLECULES = [
    ('H2',  1, 1, 2, 3),   # 190 SDs — very fast
    ('LiH', 3, 1, 4, 3),   # 367k SDs — slower but tractable
]

# R scan points relative to equilibrium: R_eq * factors
R_FACTORS = [0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 1.00, 1.05, 1.10,
             1.20, 1.40, 1.60, 2.00]


# ======================================================================
# Electronic energy computation
# ======================================================================

def compute_separated_atoms(Z: int, n_elec: int, nmax: int) -> float:
    """Compute atomic energy for a single atom."""
    from geovac.lattice_index import LatticeIndex
    li = LatticeIndex(
        n_electrons=n_elec, max_n=nmax, nuclear_charge=Z,
        vee_method='slater_full', h1_method='exact', fci_method='auto',
    )
    return li.compute_ground_state(n_states=1)[0][0]


def compute_molecular_fci(
    Z_A: int, Z_B: int, nmax: int, n_electrons: int, R: float,
    t_corr_lambda: float = 0.0, t_corr_fock_weighted: bool = False,
) -> float:
    """Compute molecular FCI energy at bond length R."""
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
        t_corr_lambda=t_corr_lambda,
        t_corr_fock_weighted=t_corr_fock_weighted,
    )
    eigvals, _ = mol.compute_ground_state(n_states=1)
    return eigvals[0]


# ======================================================================
# Per-molecule analysis
# ======================================================================

def analyze_molecule(
    name: str, Z_A: int, Z_B: int, n_electrons: int, nmax: int,
) -> dict:
    """
    Full lambda derivation for one molecule.

    Returns dict with all computed quantities.
    """
    print(f"\n{'='*70}")
    print(f"  Molecule: {name}  (Z_A={Z_A}, Z_B={Z_B}, Ne={n_electrons}, nmax={nmax})")
    print(f"{'='*70}")

    # --- Morse parameters ---
    params = compute_morse_parameters(name)
    r_e = params['r_e']
    k_morse = params['k_morse']
    D_e = params['D_e']
    a = params['a']

    print(f"\n  Morse parameters:")
    print(f"    r_e     = {r_e:.4f} bohr")
    print(f"    D_e     = {D_e:.4f} Ha ({D_e*27.2114:.2f} eV)")
    print(f"    a       = {a:.4f} bohr^-1")
    print(f"    k_Morse = {k_morse:.6f} Ha/bohr^2")
    print(f"    lambda  = {params['lam']:.2f}")

    # --- R scan points ---
    R_values = sorted(set([r_e * f for f in R_FACTORS]))
    print(f"\n  R scan: {len(R_values)} points in [{min(R_values):.3f}, {max(R_values):.3f}] bohr")

    # --- Separated atoms ---
    print(f"\n  Computing separated atoms...")
    t0 = time.time()

    # For homonuclear (Z_A == Z_B), separated = 2 * E_atom
    if Z_A == Z_B:
        n_elec_A = n_electrons // 2
        E_A = compute_separated_atoms(Z_A, n_elec_A, nmax)
        E_B = E_A
        print(f"    Homonuclear: E({name[0]}) = {E_A:.6f} Ha × 2 = {2*E_A:.6f} Ha")
    else:
        # Heteronuclear: need to figure out electron distribution
        # LiH: Li has 3e, H has 1e
        if name == 'LiH':
            E_A = compute_separated_atoms(Z_A, 3, nmax)  # Li
            E_B = compute_separated_atoms(Z_B, 1, nmax)  # H
        else:
            raise NotImplementedError(f"Electron distribution for {name} not configured")
        print(f"    E(A={Z_A}) = {E_A:.6f} Ha, E(B={Z_B}) = {E_B:.6f} Ha")

    E_sep = E_A + E_B
    print(f"    E_separated = {E_sep:.6f} Ha  ({time.time()-t0:.1f}s)")

    # --- PES scan for each lambda ---
    E_by_lambda = {lam: {} for lam in LAMBDA_VALUES}

    for i, R in enumerate(R_values):
        print(f"  [{i+1}/{len(R_values)}] R = {R:.4f} bohr ... ", end="", flush=True)
        t0 = time.time()

        for lam in LAMBDA_VALUES:
            fock_weighted = (lam > 0)
            E = compute_molecular_fci(
                Z_A, Z_B, nmax, n_electrons, R,
                t_corr_lambda=lam, t_corr_fock_weighted=fock_weighted,
            )
            E_by_lambda[lam][R] = E

        dt = time.time() - t0
        E0 = E_by_lambda[0.0][R]
        De_R = E_sep - E0
        print(f"E(bare) = {E0:.6f} Ha, D_e(R) = {De_R:.4f} Ha  ({dt:.1f}s)")

    # --- Force constant analysis ---
    print(f"\n  Force constant analysis:")
    R_arr = np.array(sorted(E_by_lambda[0.0].keys()))
    E_bare_arr = np.array([E_by_lambda[0.0][R] for R in R_arr])

    k_bare = numerical_force_constant(R_arr, E_bare_arr, r_e)
    print(f"    k_Morse (target)  = {k_morse:.6f} Ha/bohr^2")
    print(f"    k_bare (lambda=0) = {k_bare:.6f} Ha/bohr^2")

    results = {
        'name': name,
        'Z_A': Z_A, 'Z_B': Z_B,
        'n_electrons': n_electrons, 'nmax': nmax,
        'r_e': r_e, 'D_e': D_e, 'a': a,
        'k_morse': k_morse,
        'k_bare': k_bare,
        'E_sep': E_sep,
        'R_values': R_arr.tolist(),
        'E_by_lambda': {lam: {R: E_by_lambda[lam][R] for R in R_arr}
                        for lam in LAMBDA_VALUES},
        'lambda_derived': {},
    }

    # --- Lambda derivation for each probe lambda ---
    print(f"\n  Lambda derivation:")
    print(f"  {'probe_lam':>10s}  {'k_corr/lam':>12s}  {'lam_derived':>12s}  {'k_total':>10s}")

    for lam in LAMBDA_VALUES:
        if lam == 0:
            continue
        E_lam_arr = np.array([E_by_lambda[lam][R] for R in R_arr])
        E_corr_per_lam = (E_lam_arr - E_bare_arr) / lam

        result = derive_lambda_from_force_constant(
            k_morse, R_arr, E_bare_arr, E_corr_per_lam, r_e)

        ld = result['lambda_derived']
        kc = result['k_correction_per_lambda']
        kt = result['k_total']
        print(f"  {lam:10.3f}  {kc:12.6f}  {ld:12.6f}  {kt:10.6f}")

        results['lambda_derived'][lam] = {
            'lambda_derived': ld,
            'k_correction_per_lambda': kc,
            'k_total': kt,
        }

    # --- Vibrational frequency check ---
    mu_me = DIATOMIC_CONSTANTS[name]['mu_amu'] * 1822.888486
    print(f"\n  Vibrational frequency check (mu = {mu_me:.1f} m_e):")
    omega_target = params['omega_e_hartree'] * HARTREE_TO_CM
    print(f"    omega_e (target) = {omega_target:.1f} cm^-1")

    for lam in LAMBDA_VALUES:
        E_lam_arr = np.array([E_by_lambda[lam][R] for R in R_arr])
        k_lam = numerical_force_constant(R_arr, E_lam_arr, r_e)
        if k_lam > 0:
            omega = np.sqrt(k_lam / mu_me) * HARTREE_TO_CM
            print(f"    lambda={lam:.3f}: k={k_lam:.6f}, omega={omega:.1f} cm^-1")
        else:
            print(f"    lambda={lam:.3f}: k={k_lam:.6f} (no minimum)")

    return results


# ======================================================================
# Universality comparison
# ======================================================================

def compare_universality(all_results: list) -> None:
    """Compare lambda_derived across molecules."""
    print(f"\n{'='*70}")
    print("  UNIVERSALITY COMPARISON")
    print(f"{'='*70}")

    # Header
    print(f"\n  {'Molecule':>8s}  {'k_Morse':>10s}  {'k_bare':>10s}  ", end="")
    for lam in LAMBDA_VALUES:
        if lam > 0:
            print(f"{'lam_d('+str(lam)+')':>12s}  ", end="")
    print()

    # Data rows
    all_lambda = {lam: [] for lam in LAMBDA_VALUES if lam > 0}

    for res in all_results:
        print(f"  {res['name']:>8s}  {res['k_morse']:10.6f}  {res['k_bare']:10.6f}  ", end="")
        for lam in LAMBDA_VALUES:
            if lam > 0:
                ld = res['lambda_derived'].get(lam, {}).get('lambda_derived', float('nan'))
                print(f"{ld:12.6f}  ", end="")
                if np.isfinite(ld):
                    all_lambda[lam].append(ld)
        print()

    # Statistics
    print(f"\n  Summary statistics:")
    for lam in LAMBDA_VALUES:
        if lam == 0:
            continue
        vals = all_lambda[lam]
        if len(vals) >= 2:
            mean = np.mean(vals)
            std = np.std(vals)
            spread = (max(vals) - min(vals)) / mean * 100
            print(f"    probe lambda={lam:.3f}: "
                  f"mean={mean:.4f}, std={std:.4f}, spread={spread:.1f}%")
            if spread < 20:
                print(f"      -> UNIVERSAL (spread < 20%)")
            elif spread < 50:
                print(f"      -> WEAKLY universal (spread 20-50%)")
            else:
                print(f"      -> NOT universal (spread > 50%)")
        else:
            print(f"    probe lambda={lam:.3f}: only {len(vals)} molecule(s), need ≥ 2")


# ======================================================================
# Main
# ======================================================================

def main() -> None:
    print("=" * 70)
    print("  Lambda Universality Test")
    print("  Testing whether lambda_derived is the same for H2 and LiH")
    print("=" * 70)

    t_start = time.time()
    all_results = []

    for name, Z_A, Z_B, n_elec, nmax in MOLECULES:
        res = analyze_molecule(name, Z_A, Z_B, n_elec, nmax)
        all_results.append(res)

    compare_universality(all_results)

    # --- Save results ---
    outdir = os.path.join(os.path.dirname(__file__), 'data')
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, 'lambda_universality.txt')

    with open(outfile, 'w') as f:
        f.write("# Lambda Universality Test\n")
        f.write(f"# Date: {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"# Molecules: {[m[0] for m in MOLECULES]}\n")
        f.write("#\n")

        for res in all_results:
            f.write(f"\n# === {res['name']} ===\n")
            f.write(f"# Z_A={res['Z_A']}, Z_B={res['Z_B']}, "
                    f"Ne={res['n_electrons']}, nmax={res['nmax']}\n")
            f.write(f"# r_e={res['r_e']:.4f}, D_e={res['D_e']:.4f}, "
                    f"a={res['a']:.4f}\n")
            f.write(f"# k_Morse={res['k_morse']:.6f}, k_bare={res['k_bare']:.6f}\n")
            f.write(f"# E_sep={res['E_sep']:.8f}\n")

            # Lambda derivation results
            for lam, ld_data in res['lambda_derived'].items():
                f.write(f"# probe_lam={lam:.3f}: lambda_derived={ld_data['lambda_derived']:.6f}, "
                        f"k_corr/lam={ld_data['k_correction_per_lambda']:.6f}\n")

            # PES data
            f.write(f"# R  " + "  ".join(f"E(lam={l})" for l in LAMBDA_VALUES) + "\n")
            for R in res['R_values']:
                line = f"{R:8.4f}"
                for lam in LAMBDA_VALUES:
                    E = res['E_by_lambda'][lam][R]
                    line += f"  {E:14.8f}"
                f.write(line + "\n")

    print(f"\n  Results saved to {outfile}")

    dt_total = time.time() - t_start
    print(f"\n  Total time: {dt_total:.0f}s ({dt_total/60:.1f} min)")
    print("\n" + "=" * 70)
    print("  Done.")


if __name__ == '__main__':
    main()
