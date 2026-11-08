"""W1e-Projection-Audit driver.

Localize which physical effect is responsible for LiH binding in the
continuous Level 4 multichannel + PK solver but missing in the qubit
composed Hamiltonian.

Both energy expressions:

  continuous:  E(R) = E_core + V_cross_nuc(R) + E_elec(R) + V_NN_bare(R)
  qubit FCI:   E(R) = ecore(R) + <h1>_FCI + <eri>_FCI
               ecore(R) = V_NN(R) + E_core            [hydride_spec]

The continuous-side V_cross_nuc(R) (Li 1s^2 <-> H proton attraction,
asymptotically -2/R) is structurally absent from the composed-qubit
ecore.  The balanced builder MAY pick it up through the cross-block
V_ne path.  This driver runs both qubit builders on the same R-panel
as the continuous solver and reports the per-R energy breakdown so we
can see (a) the total-energy gap, (b) the "nuclear part" gap, (c) the
"electronic part" gap.  Decision gate: identify which operator class
accounts for the binding to within +/- 0.1 Ha at R_eq.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

# Force UTF-8 stdout on Windows so the analysis labels print cleanly.
try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, Exception):
    pass

from geovac.composed_diatomic import ComposedDiatomicSolver, _v_cross_nuc_1s
from geovac.composed_qubit import build_composed_hamiltonian
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.molecular_spec import lih_spec
from geovac.coupled_composition import coupled_fci_energy


# Production R panel matches R3-A (this sprint's parent diagnostic).
R_PANEL: Tuple[float, ...] = (2.5, 3.015, 3.5, 4.0, 5.0)

# E_core for Li 1s^2 (frozen-core energy used in spec).
# Matches `_FIRST_ROW_CORE_ENERGY[3]` in molecular_spec.py.
Z_LI = 3
Z_H = 1
N_CORE = 2


def run_continuous_panel() -> Dict[str, np.ndarray]:
    """Run ComposedDiatomicSolver.LiH_ab_initio on R_PANEL.

    Returns a dict with per-R arrays for E_composed, E_elec, V_cross_nuc,
    V_NN_bare, and a constant E_core.
    """
    print("=" * 72)
    print("Continuous Level 4 multichannel + PK solver")
    print("=" * 72)

    solver = ComposedDiatomicSolver.LiH_ab_initio(l_max=2, verbose=True)
    t0 = time.time()
    solver.solve_core()
    E_core = solver.E_core
    print(f"  Li core: E_core = {E_core:.6f} Ha "
          f"({time.time() - t0:.1f}s)")

    R_grid = np.array(R_PANEL)
    pes = solver.scan_pes(R_grid=R_grid, n_Re=300)
    return {
        'R': R_grid,
        'E_core': E_core,
        'E_elec': pes['E_elec'],
        'V_NN_bare': pes['V_NN_bare'],
        'V_cross_nuc': pes['V_cross_nuc'],
        'E_composed': pes['E_composed'],
    }


def _fci_energy_from_build(build_res: Dict, n_electrons: int) -> float:
    """Compute FCI ground-state energy from a build_*_hamiltonian dict.

    Wraps coupled_fci_energy from coupled_composition.
    """
    out = coupled_fci_energy(build_res, n_electrons=n_electrons,
                             verbose=False)
    return float(out['E_coupled'])


def run_composed_panel() -> Dict[str, List]:
    """Build composed-qubit Hamiltonian for LiH at each R in R_PANEL.

    Returns per-R: ecore (= V_NN(R) + E_core_Li), E_total_fci, n_qubits.
    Uses the R3-A workaround: build_composed_hamiltonian directly with
    pk_in_hamiltonian=True to bypass the ecosystem PK-bug.
    """
    print("\n" + "=" * 72)
    print("Composed qubit Hamiltonian + FCI")
    print("=" * 72)

    rows = []
    for R in R_PANEL:
        spec = lih_spec(R=R)
        t0 = time.time()
        res = build_composed_hamiltonian(spec, pk_in_hamiltonian=True,
                                         verbose=False)
        ecore = res['nuclear_repulsion']
        n_qubits = res['Q']

        E_total = _fci_energy_from_build(res, n_electrons=4)
        wall = time.time() - t0
        print(f"  R={R:6.3f}  ecore={ecore:+10.6f}  E_total={E_total:+10.6f}  "
              f"Q={n_qubits}  t={wall:.1f}s")
        rows.append({
            'R': float(R),
            'ecore': float(ecore),
            'E_total': float(E_total),
            'Q': int(n_qubits),
        })
    return {'rows': rows}


def run_balanced_panel() -> Dict[str, List]:
    """Build balanced-qubit Hamiltonian for LiH at each R in R_PANEL.

    IMPORTANT CALL CONVENTION (Path A): spec at the spec's default R
    (lih_spec() with no R kwarg => R = 3.015 bohr), and pass the requested
    R as the builder's R kwarg.  This is the convention the balanced
    builder's R-dependence correction was designed for; calling with
    spec at requested R AND builder at requested R (Path B) silently
    double-counts V_NN(R) - V_NN(R_default), breaking the PES shape.
    The R3-A and R3-B drivers used Path B (or a related broken variant
    that misplaces the nuclei used for cross-center V_ne), which is why
    they reported balanced as 'monotone descending, no binding'.
    """
    print("\n" + "=" * 72)
    print("Balanced qubit Hamiltonian + FCI  (Path A convention)")
    print("=" * 72)

    rows = []
    for R in R_PANEL:
        spec = lih_spec()  # default R = 3.015 (LiH _HYDRIDE_REQ)
        t0 = time.time()
        res = build_balanced_hamiltonian(spec, nuclei=None, R=R, verbose=False)
        ecore = res['nuclear_repulsion']
        n_qubits = res['Q']
        n_elec = res.get('n_electrons_total', 4)

        E_total = _fci_energy_from_build(res, n_electrons=n_elec)
        wall = time.time() - t0
        print(f"  R={R:6.3f}  ecore={ecore:+10.6f}  E_total={E_total:+10.6f}  "
              f"Q={n_qubits}  N_e={n_elec}  t={wall:.1f}s")
        rows.append({
            'R': float(R),
            'ecore': float(ecore),
            'E_total': float(E_total),
            'Q': int(n_qubits),
            'n_electrons': int(n_elec),
        })
    return {'rows': rows}


def analyze(cont: Dict, comp: Dict, bal: Dict) -> Dict:
    """Compute the per-R energy decomposition and dE shape."""
    print("\n" + "=" * 72)
    print("ANALYSIS — per-R energy decomposition")
    print("=" * 72)

    E_core = cont['E_core']
    R_arr = cont['R']
    E_continuous = cont['E_composed']
    E_elec_cont = cont['E_elec']
    V_cross_nuc = cont['V_cross_nuc']
    V_NN = cont['V_NN_bare']

    composed_rows = comp['rows']
    balanced_rows = bal['rows']

    table = []
    print(f"\n{'R':>6s}  {'E_cont':>10s}  {'E_qubit_c':>10s}  {'E_qubit_b':>10s}"
          f"  {'V_cross':>9s}  {'dE_c':>9s}  {'dE_b':>9s}")
    print("-" * 78)
    for i, R in enumerate(R_arr):
        E_c = float(E_continuous[i])
        E_qc = float(composed_rows[i]['E_total'])
        E_qb = float(balanced_rows[i]['E_total'])
        Vc = float(V_cross_nuc[i])
        dE_c = E_qc - E_c
        dE_b = E_qb - E_c
        print(f"{R:6.3f}  {E_c:10.4f}  {E_qc:10.4f}  {E_qb:10.4f}"
              f"  {Vc:9.4f}  {dE_c:9.4f}  {dE_b:9.4f}")
        table.append({
            'R': float(R),
            'E_continuous': E_c,
            'E_elec_cont': float(E_elec_cont[i]),
            'V_cross_nuc': Vc,
            'V_NN_bare': float(V_NN[i]),
            'E_qubit_composed': E_qc,
            'E_qubit_balanced': E_qb,
            'ecore_composed': float(composed_rows[i]['ecore']),
            'ecore_balanced': float(balanced_rows[i]['ecore']),
            'dE_composed_minus_cont': dE_c,
            'dE_balanced_minus_cont': dE_b,
        })

    # Diagnostic: does the composed gap match -V_cross_nuc?
    print("\n" + "-" * 72)
    print("Diagnostic 1 -- does composed dE match the missing V_cross_nuc?")
    print("-" * 72)
    print(f"{'R':>6s}  {'dE_c':>9s}  {'V_cross':>9s}  {'dE_c-V_cross':>14s}")
    for row in table:
        print(f"{row['R']:6.3f}  {row['dE_composed_minus_cont']:9.4f}"
              f"  {row['V_cross_nuc']:9.4f}"
              f"  {row['dE_composed_minus_cont'] - row['V_cross_nuc']:14.4f}")

    # Shape characterization
    dE_c_arr = np.array([r['dE_composed_minus_cont'] for r in table])
    dE_b_arr = np.array([r['dE_balanced_minus_cont'] for r in table])
    print(f"\ndE_composed slope vs R: "
          f"{np.polyfit(R_arr, dE_c_arr, 1)[0]:+.4f} Ha/bohr")
    print(f"dE_balanced slope vs R: "
          f"{np.polyfit(R_arr, dE_b_arr, 1)[0]:+.4f} Ha/bohr")

    # Relative-energy comparison (zeroed at R=5.0) to isolate well shape
    print("\n" + "-" * 72)
    print("Diagnostic 2 — relative energies (zeroed at R=5.0)")
    print("-" * 72)
    Eqc_arr = np.array([r['E_qubit_composed'] for r in table])
    Eqb_arr = np.array([r['E_qubit_balanced'] for r in table])
    E_cont_rel = E_continuous - E_continuous[-1]
    Eqc_rel = Eqc_arr - Eqc_arr[-1]
    Eqb_rel = Eqb_arr - Eqb_arr[-1]
    print(f"{'R':>6s}  {'cont':>9s}  {'compos':>9s}  {'balanced':>9s}")
    for i, R in enumerate(R_arr):
        print(f"{R:6.3f}  {E_cont_rel[i]:+9.4f}  {Eqc_rel[i]:+9.4f}"
              f"  {Eqb_rel[i]:+9.4f}")
    print("(neg = more bound than dissoc R=5.0; bowl shape = binding well)")

    # Per-component R-derivative analysis at R=R_eq
    print("\n" + "-" * 72)
    print("Diagnostic 3 — slope contributions (finite difference R=2.5 to R=5.0)")
    print("-" * 72)
    dR = R_arr[-1] - R_arr[0]
    slopes = {
        'E_continuous_total': (E_continuous[-1] - E_continuous[0]) / dR,
        'V_NN_bare':          (V_NN[-1] - V_NN[0]) / dR,
        'V_cross_nuc':        (V_cross_nuc[-1] - V_cross_nuc[0]) / dR,
        'E_elec_cont':        (E_elec_cont[-1] - E_elec_cont[0]) / dR,
        'E_qubit_composed':   (Eqc_arr[-1] - Eqc_arr[0]) / dR,
        'E_qubit_balanced':   (Eqb_arr[-1] - Eqb_arr[0]) / dR,
    }
    for k, v in slopes.items():
        print(f"  {k:25s}: {v:+.4f} Ha/bohr")
    print("(For binding: need positive slope toward dissoc, i.e., "
          "going UP as R increases past R_eq.)")

    # Binding well in continuous E_composed
    i_min = int(np.argmin(E_continuous))
    print(f"\nContinuous E_composed minimum at R = {R_arr[i_min]:.3f} bohr, "
          f"E = {E_continuous[i_min]:.4f} Ha")
    print(f"  diff vs R=5.0 dissoc: {E_continuous[i_min] - E_continuous[-1]:+.4f} Ha "
          "(binding well depth)")

    # Same for composed qubit
    Eqc_arr = np.array([r['E_qubit_composed'] for r in table])
    Eqb_arr = np.array([r['E_qubit_balanced'] for r in table])
    print(f"Composed qubit E_total minimum at R = {R_arr[int(np.argmin(Eqc_arr))]:.3f} bohr "
          f"(monotone descent if at panel edge)")
    print(f"Balanced qubit E_total minimum at R = {R_arr[int(np.argmin(Eqb_arr))]:.3f} bohr "
          f"(monotone descent if at panel edge)")

    return {
        'R_panel': [float(r) for r in R_arr],
        'E_core_li': float(E_core),
        'table': table,
        'composed_slope_dE_dR': float(np.polyfit(R_arr, dE_c_arr, 1)[0]),
        'balanced_slope_dE_dR': float(np.polyfit(R_arr, dE_b_arr, 1)[0]),
        'continuous_Req_index': i_min,
        'continuous_Req': float(R_arr[i_min]),
        'continuous_well_depth_vs_dissoc': float(
            E_continuous[i_min] - E_continuous[-1]),
    }


def main():
    debug_dir = Path(__file__).resolve().parent
    data_dir = debug_dir / 'data'
    data_dir.mkdir(exist_ok=True)
    out_json = data_dir / 'w1e_projection_audit.json'

    cont = run_continuous_panel()
    comp = run_composed_panel()
    bal = run_balanced_panel()
    summary = analyze(cont, comp, bal)

    payload = {
        'continuous': {k: (v.tolist() if isinstance(v, np.ndarray) else v)
                       for k, v in cont.items()},
        'composed_qubit': comp,
        'balanced_qubit': bal,
        'analysis': summary,
    }
    with open(out_json, 'w') as f:
        json.dump(payload, f, indent=2)
    print(f"\nWrote {out_json}")


if __name__ == '__main__':
    main()
