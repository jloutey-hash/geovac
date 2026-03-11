"""
Diagnostic: Angular-Weighted Cross-Atom V_ee for LiH (v0.9.26)

Compares five configurations:
  1. Baseline (use_dmatrix=False, cross_atom_vee=False)  -- from v0.9.24 data
  2. Baseline + OK unweighted (v0.9.25)                  -- from v0.9.25 data
  3. Hybrid + OK unweighted (v0.9.25)                    -- from v0.9.25 data
  4. Baseline + angular (cross_atom_vee=True, angular_weighted=True) -- NEW
  5. Hybrid + angular (use_dmatrix='hybrid', angular_weighted=True)  -- NEW

BSSE is R-independent (ghost atoms have no bridges/cross-nuclear), so we
compute it once at R=3.015 and apply to all R values.

Reports PES at R = 2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0, 6.0 bohr.
"""

import warnings
import sys
import time
import numpy as np

sys.path.insert(0, '.')
from geovac.lattice_index import (
    MolecularLatticeIndex, LatticeIndex, compute_bsse_correction,
)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
NMAX = 3
Z_A, Z_B = 3, 1
N_ELECTRONS = 4
R_VALUES = [2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0, 6.0]
R_EQ_EXPT = 3.015
D_E_EXPT = 0.0924  # Ha

# Configs to run (4 and 5 are new)
NEW_CONFIGS = [
    {'label': 'baseline+angular', 'use_dmatrix': False,
     'cross_atom_vee': True, 'angular_weighted': True},
    {'label': 'hybrid+angular',   'use_dmatrix': 'hybrid',
     'cross_atom_vee': True, 'angular_weighted': True},
]

outpath = 'debug/data/angular_vee_lih_results.txt'


def compute_atomic_energies() -> dict:
    """Compute isolated Li and H energies at nmax=3."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        li = LatticeIndex(
            n_electrons=3, max_n=NMAX, nuclear_charge=Z_A,
            vee_method='slater_full', h1_method='exact', fci_method='auto',
        )
        E_Li = li.compute_ground_state(n_states=1)[0][0]

        h = LatticeIndex(
            n_electrons=1, max_n=NMAX, nuclear_charge=Z_B,
            vee_method='slater_full', h1_method='exact', fci_method='matrix',
        )
        E_H = h.compute_ground_state(n_states=1)[0][0]

    return {'E_Li': E_Li, 'E_H': E_H, 'E_atoms': E_Li + E_H}


def run_single(R: float, config: dict) -> dict:
    """Run single FCI at given R with given config."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mol = MolecularLatticeIndex(
            Z_A=Z_A, Z_B=Z_B, nmax_A=NMAX, nmax_B=NMAX,
            R=R, n_electrons=N_ELECTRONS,
            vee_method='slater_full', fci_method='auto',
            use_dmatrix=config['use_dmatrix'],
            cross_atom_vee=config['cross_atom_vee'],
            angular_weighted=config.get('angular_weighted', True),
        )
        eigvals, eigvecs = mol.compute_ground_state(n_states=1)
    return {
        'E_mol': eigvals[0],
        'mol': mol,
        'civec': eigvecs[:, 0],
    }


def main():
    t_start = time.perf_counter()
    lines = []

    def log(s: str = ''):
        print(s)
        lines.append(s)

    log("=" * 80)
    log("LiH Angular-Weighted Cross-Atom V_ee Diagnostic (v0.9.26)")
    log("nmax=%d, Z_A=%d, Z_B=%d, N_e=%d" % (NMAX, Z_A, Z_B, N_ELECTRONS))
    log("=" * 80)

    # Atomic energies
    atoms = compute_atomic_energies()
    log("\nAtomic energies (nmax=%d):" % NMAX)
    log("  E(Li) = %.6f Ha" % atoms['E_Li'])
    log("  E(H)  = %.6f Ha" % atoms['E_H'])
    log("  E(Li+H) = %.6f Ha" % atoms['E_atoms'])

    # Angular weighting table at R=3.015
    log("\n" + "=" * 80)
    log("Angular-weighted ERI table at R=%.3f" % R_EQ_EXPT)
    log("f(l_A, l_B) = 1 / ((2*l_A+1) * (2*l_B+1))")
    log("=" * 80)
    log("  %5s %5s %5s %5s %8s %8s %10s %10s %10s" % (
        'n_A', 'l_A', 'n_B', 'l_B', '<r>_A', '<r>_B',
        'J_raw', 'f(l,l)', 'J_weighted'))
    log("  " + "-" * 85)

    for n_a in range(1, NMAX + 1):
        for l_a in range(n_a):
            for n_b in range(1, NMAX + 1):
                for l_b in range(n_b):
                    # Skip pure s-s (handled by Fourier)
                    if l_a == 0 and l_b == 0:
                        continue
                    r_a = float(n_a)**2 / Z_A
                    r_b = float(n_b)**2 / Z_B
                    R_eff = np.sqrt(R_EQ_EXPT**2 + (r_a + r_b)**2)
                    j_raw = 1.0 / R_eff
                    f_ang = 1.0 / ((2 * l_a + 1) * (2 * l_b + 1))
                    j_weighted = j_raw * f_ang
                    log("  %5d %5d %5d %5d %8.3f %8.3f %10.6f %10.4f %10.6f" % (
                        n_a, l_a, n_b, l_b, r_a, r_b, j_raw, f_ang, j_weighted))

    # BSSE at R=3.015 (R-independent, compute once)
    log("\n" + "=" * 80)
    log("BSSE (computed once at R=3.015, R-independent)")
    log("=" * 80)
    bsse = compute_bsse_correction(
        Z_A=Z_A, Z_B=Z_B, nmax_A=NMAX, nmax_B=NMAX,
        R=R_EQ_EXPT, n_electrons_A=3, n_electrons_B=1,
        vee_method='slater_full', fci_method='auto',
        cross_atom_vee=False,  # ghost atoms skip cross-atom V_ee anyway
    )
    BSSE = bsse['BSSE']
    log("  BSSE = %.6f Ha (Li: %.6f, H: %.6f)" % (
        BSSE, bsse['BSSE_A'], bsse['BSSE_B']))

    # PES for new configurations
    all_results = {}
    for cfg in NEW_CONFIGS:
        label = cfg['label']
        log("\n" + "=" * 80)
        log("Configuration: %s" % label)
        log("  use_dmatrix=%s, cross_atom_vee=%s, angular_weighted=%s" % (
            cfg['use_dmatrix'], cfg['cross_atom_vee'],
            cfg.get('angular_weighted', True)))
        log("=" * 80)

        pes = {}
        for R in R_VALUES:
            t0 = time.perf_counter()
            result = run_single(R, cfg)
            dt = time.perf_counter() - t0
            pes[R] = result
            log("  R=%.3f: E=%.6f Ha (%.1fs)" % (R, result['E_mol'], dt))

        all_results[label] = pes

        log("\n  %8s  %12s  %10s  %10s" % ('R', 'E_mol', 'D_e_raw', 'D_e_CP'))
        log("  %s  %s  %s  %s" % ('-'*8, '-'*12, '-'*10, '-'*10))
        for R in R_VALUES:
            E = pes[R]['E_mol']
            D_e_raw = atoms['E_atoms'] - E
            D_e_CP = D_e_raw + BSSE
            log("  %8.3f  %12.6f  %10.6f  %10.6f" % (R, E, D_e_raw, D_e_CP))

        # Find R_eq
        E_min = min(pes[R]['E_mol'] for R in R_VALUES)
        R_eq = min(R_VALUES, key=lambda R: pes[R]['E_mol'])
        D_e_raw_eq = atoms['E_atoms'] - E_min
        D_e_CP_eq = D_e_raw_eq + BSSE
        log("\n  R_eq ~ %.3f bohr, E_min = %.6f Ha" % (R_eq, E_min))
        log("  D_e_raw = %.6f Ha, D_e_CP = %.6f Ha (expt: %.4f Ha)" % (
            D_e_raw_eq, D_e_CP_eq, D_E_EXPT))

    # 5-config comparison table at R=3.015
    log("\n" + "=" * 80)
    log("5-config comparison at R=3.015 bohr")
    log("=" * 80)

    # Historical data from previous runs
    E_baseline_3015 = -8.116739   # v0.9.24 baseline
    E_ok_baseline_3015 = -7.677334  # v0.9.25 baseline+OK
    E_ok_hybrid_3015 = -7.752513    # v0.9.25 hybrid+OK

    log("\n  %25s  %12s  %10s  %10s" % ('Config', 'E_mol', 'D_e_raw', 'D_e_CP'))
    log("  %s  %s  %s  %s" % ('-'*25, '-'*12, '-'*10, '-'*10))

    # 1. Baseline (v0.9.24)
    D_raw = atoms['E_atoms'] - E_baseline_3015
    log("  %25s  %12.6f  %10.6f  %10.6f" % (
        'baseline (v0.9.24)', E_baseline_3015, D_raw, D_raw + BSSE))

    # 2. Baseline+OK unweighted (v0.9.25)
    D_raw = atoms['E_atoms'] - E_ok_baseline_3015
    log("  %25s  %12.6f  %10.6f  %10.6f" % (
        'baseline+OK (v0.9.25)', E_ok_baseline_3015, D_raw, D_raw + BSSE))

    # 3. Hybrid+OK unweighted (v0.9.25)
    D_raw = atoms['E_atoms'] - E_ok_hybrid_3015
    log("  %25s  %12.6f  %10.6f  %10.6f" % (
        'hybrid+OK (v0.9.25)', E_ok_hybrid_3015, D_raw, D_raw + BSSE))

    # 4-5. New angular configs
    for cfg in NEW_CONFIGS:
        label = cfg['label']
        E = all_results[label][R_EQ_EXPT]['E_mol']
        D_e_raw = atoms['E_atoms'] - E
        D_e_CP = D_e_raw + BSSE
        log("  %25s  %12.6f  %10.6f  %10.6f" % (label, E, D_e_raw, D_e_CP))

    # Term decomposition for Config A (baseline+angular) at R=2.5 and R=3.015
    cfg_a = NEW_CONFIGS[0]
    pes_a = all_results[cfg_a['label']]

    log("\n" + "=" * 80)
    log("Term decomposition: %s" % cfg_a['label'])
    log("=" * 80)

    decomps = {}
    for R in [2.5, R_EQ_EXPT]:
        mol = pes_a[R]['mol']
        civec = pes_a[R]['civec']
        E = pes_a[R]['E_mol']
        decomp = mol.decompose_energy(civec, E)
        decomps[R] = decomp

        log("\n  R = %.3f bohr:" % R)
        for key in ['T', 'V_nA', 'V_nB', 'V_cross_A', 'V_cross_B',
                     'V_bridge', 'V_ee', 'V_NN']:
            log("    %-12s = %+.6f Ha" % (key, decomp[key]))
        log("    %-12s = %.6f Ha" % ('E_total', decomp['E_total']))

    # Screening deficit analysis
    d25 = decomps[2.5]
    d30 = decomps[R_EQ_EXPT]

    log("\n" + "=" * 80)
    log("Screening deficit analysis: %s" % cfg_a['label'])
    log("Delta = X(R=2.5) - X(R=3.015)")
    log("=" * 80)

    for key in ['T', 'V_nA', 'V_nB', 'V_cross_A', 'V_cross_B',
                 'V_bridge', 'V_ee', 'V_NN']:
        delta = d25[key] - d30[key]
        log("  Delta(%-12s) = %+.6f Ha" % (key, delta))

    delta_cross = ((d25['V_cross_A'] + d25['V_cross_B'])
                   - (d30['V_cross_A'] + d30['V_cross_B']))
    delta_vee = d25['V_ee'] - d30['V_ee']
    deficit = delta_cross + delta_vee

    log("\n  Delta(V_cross_A + V_cross_B) = %+.6f Ha" % delta_cross)
    log("  Delta(V_ee)                  = %+.6f Ha" % delta_vee)
    log("  Net screening deficit        = %+.6f Ha" % deficit)
    log("  (v0.9.24 baseline deficit was -0.247 Ha)")
    log("  (v0.9.25 OK unweighted deficit was -0.107 Ha)")

    # Diagnostics summary
    log("\n" + "=" * 80)
    log("Diagnostics Summary")
    log("=" * 80)
    for cfg in NEW_CONFIGS:
        label = cfg['label']
        pes = all_results[label]
        R_eq = min(R_VALUES, key=lambda R: pes[R]['E_mol'])
        E_diss = pes[6.0]['E_mol']
        D_e_raw_diss = atoms['E_atoms'] - E_diss
        D_e_CP_diss = D_e_raw_diss + BSSE

        E_eq = pes[R_eq]['E_mol']
        D_e_raw_eq = atoms['E_atoms'] - E_eq
        D_e_CP_eq = D_e_raw_eq + BSSE

        log("\n  %s:" % label)
        log("    R_eq ~ %.3f bohr (expt: %.3f)" % (R_eq, R_EQ_EXPT))
        log("    D_e_CP(R_eq) = %.6f Ha (expt: %.4f Ha)" % (D_e_CP_eq, D_E_EXPT))
        log("    |D_e_CP(R=6)| = %.6f Ha (%s)" % (
            abs(D_e_CP_diss), 'PASS' if abs(D_e_CP_diss) <= 0.001 else 'check'))
        log("    Bound: %s" % ('YES' if D_e_CP_eq > 0 else 'NO'))
        log("    R_eq in [2.5, 3.5]: %s" % (
            'YES' if 2.5 <= R_eq <= 3.5 else 'NO'))

    elapsed = time.perf_counter() - t_start
    log("\nTotal time: %.1fs" % elapsed)

    # Write results
    with open(outpath, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print("\nResults written to %s" % outpath)


if __name__ == '__main__':
    main()
