"""
NaH balanced-coupled PES with Phillips-Kleinman cross-center barrier.

Tests whether enabling PK orthogonalization of the H valence orbital
against the [Ne] core on Na produces a bound equilibrium R_eq.

Comparison configurations:
  1. screened=True, pk=False      (Track 2 baseline -- W1c only)
  2. screened=True, pk=True       (this sprint -- W1c + W1b PK)
  3. screened=False, pk=True      (PK only, for comparison)
  4. screened=False, pk=False     (Sprint 7 bare baseline)

Q=20 for NaH at max_n=2, 2-electron FCI sector. Minutes-scale compute.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(line_buffering=True)

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy
from geovac.molecular_spec import nah_spec


def run_single(R: float, screened: bool, pk: bool):
    spec = nah_spec(R=R, max_n=2)
    n_e = sum(b.n_electrons for b in spec.blocks)
    ham = build_balanced_hamiltonian(
        spec, R=R, n_grid_vne=4000, L_max=4,
        screened_cross_center=screened,
        pk_cross_center=pk,
        verbose=False,
    )
    fci = coupled_fci_energy(ham, n_electrons=n_e, verbose=False)
    pk_trace = float(np.trace(ham['h1_pk_cross'])) if pk else 0.0
    return {
        'R': float(R),
        'screened': screened,
        'pk': pk,
        'M': int(ham['M']),
        'Q': int(ham['Q']),
        'n_electrons': int(n_e),
        'E_total': float(fci['E_coupled']),
        'V_NN_plus_E_core': float(ham['nuclear_repulsion']),
        'pk_trace': pk_trace,
        'N_pauli': int(ham['N_pauli']),
        'one_norm': float(ham['one_norm']),
    }


def has_minimum(pts):
    if len(pts) < 3:
        return None
    E_arr = [p['E_total'] for p in pts]
    i_min = int(np.argmin(E_arr))
    is_internal = 0 < i_min < len(E_arr) - 1
    return {
        'i_min': i_min,
        'R_min': float(pts[i_min]['R']),
        'E_min': float(E_arr[i_min]),
        'has_internal_min': is_internal,
        'descent_depth': float(max(E_arr) - min(E_arr)),
    }


def parabolic_R_eq(pts):
    """Fit parabola to R_eq region and return refined minimum."""
    if len(pts) < 3:
        return None
    R_arr = np.array([p['R'] for p in pts])
    E_arr = np.array([p['E_total'] for p in pts])
    i_min = int(np.argmin(E_arr))
    if i_min == 0 or i_min == len(R_arr) - 1:
        return None
    # 3-point parabolic fit
    R3 = R_arr[i_min - 1: i_min + 2]
    E3 = E_arr[i_min - 1: i_min + 2]
    # E = a*(R-R0)^2 + b
    a, b, c = np.polyfit(R3, E3, 2)
    if a <= 0:
        return None
    R_eq = -b / (2 * a)
    E_eq = a * R_eq ** 2 + b * R_eq + c
    return {
        'R_eq': float(R_eq),
        'E_eq': float(E_eq),
        'curvature': float(a),
    }


def main():
    out_path = Path('debug/data/pk_cross_center_nah.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)

    R_grid = [2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0]
    R_eq_exp = 3.566
    D_e_exp = 0.075  # Ha (approximate)

    print(f"=== NaH balanced FCI, max_n=2 ({len(R_grid)} R-points) ===",
          flush=True)
    print(f"  Reference: R_eq = {R_eq_exp:.3f} bohr, D_e ~ {D_e_exp:.3f} Ha",
          flush=True)

    configs = [
        ('screened+PK',   True,  True),
        ('screened only', True,  False),
        ('PK only',       False, True),
        ('bare',          False, False),
    ]
    all_pts = {}
    all_diag = {}
    all_fit = {}

    for label, scr, pk in configs:
        print(f"\n--- {label} (screened={scr}, pk={pk}) ---", flush=True)
        pts = []
        for R in R_grid:
            t0 = time.perf_counter()
            try:
                pt = run_single(R, scr, pk)
                pts.append(pt)
                pk_str = f", PK_tr={pt['pk_trace']:.3e}" if pk else ""
                print(f"  R={R:6.3f}: E={pt['E_total']:.6f} Ha "
                      f"(VNN={pt['V_NN_plus_E_core']:.4f}{pk_str}) "
                      f"[{time.perf_counter()-t0:.1f}s]", flush=True)
            except Exception as e:
                print(f"  R={R:6.3f}: FAIL {e}", flush=True)
        all_pts[label] = pts
        diag = has_minimum(pts)
        all_diag[label] = diag
        if diag and diag['has_internal_min']:
            fit = parabolic_R_eq(pts)
            all_fit[label] = fit

    # Summary
    print(f"\n=== Summary ===", flush=True)
    for label, _, _ in configs:
        diag = all_diag.get(label)
        fit = all_fit.get(label)
        if diag is None:
            print(f"  {label:20s}: insufficient data", flush=True)
            continue
        marker = "BOUND" if diag['has_internal_min'] else "MONOTONIC"
        line = (f"  {label:20s}: {marker:10s} "
                f"R_min={diag['R_min']:.3f} bohr "
                f"E_min={diag['E_min']:.6f} Ha "
                f"depth={diag['descent_depth']:.4f} Ha")
        if fit:
            R_err = (fit['R_eq'] - R_eq_exp) / R_eq_exp * 100
            line += f"  parabolic R_eq={fit['R_eq']:.3f} ({R_err:+.1f}%)"
        print(line, flush=True)

    # Verdict
    pk_diag = all_diag.get('screened+PK')
    pk_bound = pk_diag and pk_diag['has_internal_min']
    print(f"\n=== Verdict ===", flush=True)
    if pk_bound:
        fit = all_fit.get('screened+PK')
        print(f"  PK CROSS-CENTER CLOSES THE WALL.", flush=True)
        if fit:
            print(f"  Parabolic R_eq = {fit['R_eq']:.3f} bohr "
                  f"vs experimental {R_eq_exp:.3f} bohr", flush=True)
    else:
        print(f"  PK cross-center does NOT close the wall (still monotonic).",
              flush=True)
        # Quantify the PK contribution
        bare_diag = all_diag.get('bare')
        scr_diag = all_diag.get('screened only')
        pk_diag = all_diag.get('screened+PK')
        if bare_diag and scr_diag and pk_diag:
            print(f"  Descent depth: bare {bare_diag['descent_depth']:.3f} Ha "
                  f"-> screened {scr_diag['descent_depth']:.3f} Ha "
                  f"-> +PK {pk_diag['descent_depth']:.3f} Ha", flush=True)

    out = {
        'R_grid': R_grid,
        'R_eq_exp': R_eq_exp,
        'D_e_exp': D_e_exp,
        'configs': [
            {'label': label, 'screened': scr, 'pk': pk}
            for label, scr, pk in configs
        ],
        'pts': all_pts,
        'diag': all_diag,
        'fit': all_fit,
        'pk_closes_wall': bool(pk_bound),
    }
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\n[saved] {out_path}", flush=True)


if __name__ == '__main__':
    main()
