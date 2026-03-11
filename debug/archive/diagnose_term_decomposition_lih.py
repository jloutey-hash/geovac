#!/usr/bin/env python
"""
Diagnose R_eq shift in LiH: Hamiltonian term decomposition.

v0.9.24 — Pure diagnostic, no physics changes.

Question: In the v0.9.11 baseline (use_dmatrix='standard', use_shell_radius=False),
LiH has D_e_CP = 0.093 Ha (1% error) but R_eq ~ 2.5 instead of experimental 3.015.
Which Hamiltonian term pulls R_eq inward?

Method: Build MolecularLatticeIndex at each R, solve FCI, compute 1-RDM,
decompose energy into T, V_nA, V_nB, V_cross_A, V_cross_B, V_bridge, V_ee, V_NN.
"""
import sys
import os
import warnings

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from geovac.lattice_index import MolecularLatticeIndex

# Baseline configuration (v0.9.11)
Z_A, Z_B = 3, 1
NMAX = 3
N_ELECTRONS = 4
N_BRIDGES = 20

R_VALUES = [2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0, 6.0]
R_REF = 3.015  # experimental equilibrium

OUTFILE = os.path.join(os.path.dirname(__file__), 'data',
                       'term_decomposition_lih.txt')


def run_decomposition(R: float) -> dict:
    """Run FCI and decompose energy at a given R."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mol = MolecularLatticeIndex(
            Z_A=Z_A, Z_B=Z_B,
            nmax_A=NMAX, nmax_B=NMAX,
            R=R,
            n_electrons=N_ELECTRONS,
            n_bridges=N_BRIDGES,
            vee_method='slater_full',
            fci_method='auto',
            use_dmatrix=False,
            use_shell_radius=False,
        )
    eigvals, eigvecs = mol.compute_ground_state(n_states=1)
    civec = eigvecs[:, 0]
    E_total = eigvals[0]

    decomp = mol.decompose_energy(civec, E_total)
    decomp['R'] = R
    return decomp


def main() -> None:
    results = []
    for R in R_VALUES:
        print(f"\n{'='*60}")
        print(f"  R = {R:.3f} bohr")
        print(f"{'='*60}")
        d = run_decomposition(R)
        results.append(d)

    # Find reference (R=3.015) result
    ref = None
    for d in results:
        if abs(d['R'] - R_REF) < 0.001:
            ref = d
            break

    lines = []
    lines.append("=" * 72)
    lines.append("LiH Hamiltonian Term Decomposition (v0.9.24)")
    lines.append(f"Baseline: use_dmatrix=False, use_shell_radius=False, nmax={NMAX}")
    lines.append(f"Z_A={Z_A} (Li), Z_B={Z_B} (H), N_e={N_ELECTRONS}")
    lines.append("=" * 72)
    lines.append("")

    # Full decomposition at each R
    for d in results:
        R = d['R']
        V_total = (d['V_nA'] + d['V_nB'] + d['V_cross_A'] + d['V_cross_B']
                   + d['V_ee'] + d['V_NN'])
        virial = -2.0 * d['T'] / V_total if abs(V_total) > 1e-10 else float('nan')

        lines.append(f"R = {R:.3f} bohr:")
        lines.append(f"  T          = {d['T']:+.6f} Ha")
        lines.append(f"  V_nA       = {d['V_nA']:+.6f} Ha")
        lines.append(f"  V_nB       = {d['V_nB']:+.6f} Ha")
        lines.append(f"  V_cross_A  = {d['V_cross_A']:+.6f} Ha")
        lines.append(f"  V_cross_B  = {d['V_cross_B']:+.6f} Ha")
        lines.append(f"  V_bridge   = {d['V_bridge']:+.6f} Ha")
        lines.append(f"  V_ee       = {d['V_ee']:+.6f} Ha")
        lines.append(f"  V_NN       = {d['V_NN']:+.6f} Ha")
        lines.append(f"  E_total    = {d['E_total']:+.6f} Ha  "
                      f"(check: {d['E_check']:+.6f} Ha)")
        lines.append(f"  H1 total   = {d['H1_total']:+.6f} Ha  "
                      f"(diag sum: {d['H1_diag_sum']:+.6f} Ha)")
        lines.append(f"  Virial     = {virial:.4f}")
        lines.append(f"  Bridges    = {d['n_bridge_active']} active, "
                      f"|max| = {d['max_bridge_elem']:.6f} Ha")
        lines.append("")

    # Delta table: ΔX(R) = X(R) - X(R_ref)
    lines.append("")
    lines.append("=" * 72)
    lines.append(f"Delta table: DX(R) = X(R) - X(R={R_REF})")
    lines.append("=" * 72)
    lines.append("")

    components = ['T', 'V_nA', 'V_nB', 'V_cross_A', 'V_cross_B',
                  'V_bridge', 'V_ee', 'V_NN', 'E_total']
    header = f"{'R':>6s}"
    for c in components:
        header += f"  {c:>10s}"
    lines.append(header)
    lines.append("-" * len(header))

    for d in results:
        row = f"{d['R']:6.3f}"
        for c in components:
            delta = d[c] - ref[c]
            row += f"  {delta:+10.5f}"
        lines.append(row)
    lines.append("")

    # Identify steepest component between R=2.5 and R=3.015
    d25 = None
    for d in results:
        if abs(d['R'] - 2.5) < 0.01:
            d25 = d
            break

    lines.append("=" * 72)
    lines.append("Slope analysis: DX = X(2.5) - X(3.015) for each component")
    lines.append("=" * 72)
    lines.append("")

    slopes = []
    for c in components:
        if c == 'E_total':
            continue
        delta = d25[c] - ref[c]
        slopes.append((c, delta))
        lines.append(f"  {c:>12s}: {delta:+.6f} Ha")

    # Sort by most negative (pulling energy down at R=2.5)
    slopes.sort(key=lambda x: x[1])
    lines.append("")
    lines.append("Sorted by magnitude (most negative = strongest pull to short R):")
    for c, delta in slopes:
        lines.append(f"  {c:>12s}: {delta:+.6f} Ha")

    # Root cause statement
    most_neg_comp, most_neg_val = slopes[0]

    # V_NN change (repulsive, positive direction at short R)
    vnn_delta = d25['V_NN'] - ref['V_NN']
    e_total_delta = d25['E_total'] - ref['E_total']

    lines.append("")
    lines.append("=" * 72)
    lines.append("ROOT CAUSE STATEMENT")
    lines.append("=" * 72)
    lines.append("")
    lines.append(
        f"The R_eq is too short because {most_neg_comp} decreases by "
        f"{abs(most_neg_val):.4f} Ha between R=3.015 and R=2.5, "
        f"while V_NN increases by only {vnn_delta:.4f} Ha, producing a "
        f"net energy change of {e_total_delta:+.4f} Ha."
    )
    if e_total_delta < 0:
        lines.append(
            f"Since DeltaE < 0 at R=2.5, the {most_neg_comp} term "
            f"overcompensates V_NN repulsion, pulling the minimum inward."
        )
    else:
        lines.append(
            f"Since DeltaE > 0 at R=2.5, the equilibrium is actually at "
            f"R > 2.5.  The true R_eq lies between 2.5 and 3.015."
        )

    # Additional analysis: cross-nuclear attraction vs V_ee screening
    vcross_total_delta = ((d25['V_cross_A'] - ref['V_cross_A'])
                          + (d25['V_cross_B'] - ref['V_cross_B']))
    vee_delta = d25['V_ee'] - ref['V_ee']
    lines.append("")
    lines.append("Cross-nuclear vs V_ee screening:")
    lines.append(f"  Delta(V_cross_A + V_cross_B) = {vcross_total_delta:+.6f} Ha")
    lines.append(f"  Delta(V_ee)                  = {vee_delta:+.6f} Ha")
    lines.append(f"  Net screening deficit        = {vcross_total_delta + vee_delta:+.6f} Ha")

    # Bridge analysis
    lines.append("")
    lines.append("=" * 72)
    lines.append("Bridge analysis across R")
    lines.append("=" * 72)
    lines.append("")
    lines.append(f"{'R':>6s}  {'N_bridge':>8s}  {'|max_elem|':>10s}  {'V_bridge':>10s}")
    lines.append("-" * 40)
    for d in results:
        lines.append(f"{d['R']:6.3f}  {d['n_bridge_active']:8d}  "
                      f"{d['max_bridge_elem']:10.6f}  {d['V_bridge']:+10.6f}")

    # Virial analysis
    lines.append("")
    lines.append("=" * 72)
    lines.append("Virial ratio -2<T>/<V> at each R")
    lines.append("=" * 72)
    lines.append("")
    for d in results:
        V_total = (d['V_nA'] + d['V_nB'] + d['V_cross_A'] + d['V_cross_B']
                   + d['V_ee'] + d['V_NN'])
        virial = -2.0 * d['T'] / V_total if abs(V_total) > 1e-10 else float('nan')
        lines.append(f"  R = {d['R']:6.3f}:  virial = {virial:.4f}")

    output = "\n".join(lines)
    print("\n" + output)

    os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)
    with open(OUTFILE, 'w') as f:
        f.write(output + "\n")
    print(f"\nSaved to {OUTFILE}")


if __name__ == '__main__':
    main()
