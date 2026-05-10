"""
NaH balanced-coupled PES with screened-Schrödinger valence basis correction.

Track 3 of the post-Track-2 multi-track sprint (2026-05-09): tests whether
replacing the hydrogenic Z_orb=1 Na valence h1 diagonal with the actual
Schrödinger eigenvalue of the FrozenCore Z_eff(r) potential closes the
W1c-residual orthogonality wall.

Comparison configurations (8 total):

  1. screened-valence + W1c + PK   (full architecture)
  2. screened-valence + W1c        (no PK)
  3. screened-valence + PK         (no W1c)
  4. screened-valence only         (no W1c, no PK)
  5. W1c + PK                      (Track 2 / PK-cross-center sprint baseline)
  6. W1c only                      (Track 2 baseline)
  7. PK only                       (no W1c, no screened-valence)
  8. bare                          (Sprint 7 baseline)

Q=20 for NaH at max_n=2, 2-electron FCI sector. 7-point R-grid.
Each config takes ~30 s for 7 R-points; total ~4 min.

Architecturally: this isolates the contribution of each engineering
closure layer to the descent depth and tests for an internal equilibrium
minimum.
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


def run_single(R: float, screened: bool, pk: bool, valence: bool):
    spec = nah_spec(R=R, max_n=2)
    n_e = sum(b.n_electrons for b in spec.blocks)
    ham = build_balanced_hamiltonian(
        spec, R=R, n_grid_vne=4000, L_max=4,
        screened_cross_center=screened,
        pk_cross_center=pk,
        screened_valence_basis=valence,
        verbose=False,
    )
    fci = coupled_fci_energy(ham, n_electrons=n_e, verbose=False)
    pk_trace = float(np.trace(ham['h1_pk_cross'])) if pk else 0.0
    sv_trace = ham['screened_valence_info']['total_trace_shift']
    return {
        'R': float(R),
        'screened': screened,
        'pk': pk,
        'valence': valence,
        'M': int(ham['M']),
        'Q': int(ham['Q']),
        'n_electrons': int(n_e),
        'E_total': float(fci['E_coupled']),
        'V_NN_plus_E_core': float(ham['nuclear_repulsion']),
        'pk_trace': pk_trace,
        'sv_trace_shift': sv_trace,
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
    if len(pts) < 3:
        return None
    R_arr = np.array([p['R'] for p in pts])
    E_arr = np.array([p['E_total'] for p in pts])
    i_min = int(np.argmin(E_arr))
    if i_min == 0 or i_min == len(R_arr) - 1:
        return None
    R3 = R_arr[i_min - 1: i_min + 2]
    E3 = E_arr[i_min - 1: i_min + 2]
    a, b, c = np.polyfit(R3, E3, 2)
    if a <= 0:
        return None
    R_eq = -b / (2 * a)
    E_eq = a * R_eq ** 2 + b * R_eq + c
    # Estimate D_e as the depth from the asymptote (largest R) to the min
    E_asymp = max(E_arr)
    D_e = E_asymp - E_eq
    return {
        'R_eq': float(R_eq),
        'E_eq': float(E_eq),
        'D_e': float(D_e),
        'curvature': float(a),
    }


def main():
    out_path = Path('debug/data/w1c_residual_nah_track3_results.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)

    R_grid = [2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0]
    R_eq_exp = 3.566
    D_e_exp = 0.075

    print(f"=== W1c-residual NaH Track 3: PES sweep ===", flush=True)
    print(f"  Reference: R_eq={R_eq_exp:.3f} bohr, D_e~{D_e_exp:.4f} Ha "
          f"(experimental NaH)", flush=True)
    print(f"  R_grid: {R_grid}", flush=True)
    print(f"  max_n=2, Q=20, 2e FCI sector", flush=True)

    # 8 configurations
    configs = [
        ('all (sv+scr+pk)', True,  True,  True),
        ('sv + scr',        True,  False, True),
        ('sv + pk',         False, True,  True),
        ('sv only',         False, False, True),
        ('scr + pk',        True,  True,  False),
        ('scr only',        True,  False, False),
        ('pk only',         False, True,  False),
        ('bare',            False, False, False),
    ]
    all_pts = {}
    all_diag = {}
    all_fit = {}

    for label, scr, pk, val in configs:
        print(f"\n--- {label} (scr={scr}, pk={pk}, val={val}) ---", flush=True)
        pts = []
        for R in R_grid:
            t0 = time.perf_counter()
            try:
                pt = run_single(R, scr, pk, val)
                pts.append(pt)
                tags = []
                if pk:
                    tags.append(f"PKtr={pt['pk_trace']:+.3e}")
                if val:
                    tags.append(f"SVtr={pt['sv_trace_shift']:+.3e}")
                tag_str = ', '.join(tags)
                print(f"  R={R:6.3f}: E={pt['E_total']:+12.6f} Ha "
                      f"(VNN={pt['V_NN_plus_E_core']:+.4f}; {tag_str}) "
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
    print(f"{'Config':<22s} {'Marker':<10s} {'R_min':>8s} "
          f"{'E_min':>14s} {'depth':>10s} {'R_eq fit':>14s} "
          f"{'D_e fit':>10s}", flush=True)
    print('-' * 92, flush=True)

    for label, _, _, _ in configs:
        diag = all_diag.get(label)
        fit = all_fit.get(label)
        if diag is None:
            print(f"  {label:<22s} insufficient data", flush=True)
            continue
        marker = "BOUND" if diag['has_internal_min'] else "DESC"
        line = (f"  {label:<22s} {marker:<10s} {diag['R_min']:>8.3f} "
                f"{diag['E_min']:>14.6f} {diag['descent_depth']:>10.4f}")
        if fit:
            r_err = (fit['R_eq'] - R_eq_exp) / R_eq_exp * 100
            line += f"  R_eq={fit['R_eq']:>6.3f} ({r_err:+5.1f}%) D_e={fit['D_e']:>6.3f}"
        print(line, flush=True)

    # Verdict
    sv_diag = all_diag.get('sv only')
    sv_scr_diag = all_diag.get('sv + scr')
    full_diag = all_diag.get('all (sv+scr+pk)')
    print(f"\n=== Verdict ===", flush=True)
    bound_anywhere = any(d and d['has_internal_min'] for d in all_diag.values())
    if bound_anywhere:
        bound_configs = [
            label for label, _, _, _ in configs
            if all_diag.get(label) and all_diag[label]['has_internal_min']
        ]
        print(f"  BINDING RECOVERED in: {bound_configs}", flush=True)
        for cfg in bound_configs:
            fit = all_fit.get(cfg)
            if fit:
                r_err = (fit['R_eq'] - R_eq_exp) / R_eq_exp * 100
                d_err = (fit['D_e'] - D_e_exp) / D_e_exp * 100
                print(f"    {cfg}: R_eq={fit['R_eq']:.3f} bohr ({r_err:+.1f}%), "
                      f"D_e={fit['D_e']:.4f} Ha ({d_err:+.1f}%)", flush=True)
    else:
        print(f"  NaH still unbound under all 8 configurations.", flush=True)
        print(f"  W1c-residual orthogonality wall is structurally deeper "
              f"than the framework's hydrogenic-vs-screened diagonal h1 "
              f"correction can address.", flush=True)
        # Quantify the screened-valence contribution
        bare_diag = all_diag.get('bare')
        scr_diag = all_diag.get('scr only')
        scr_pk_diag = all_diag.get('scr + pk')
        if bare_diag and scr_diag and full_diag:
            print(f"\n  Descent-depth ladder:")
            print(f"    bare           : {bare_diag['descent_depth']:.4f} Ha")
            print(f"    + W1c          : {scr_diag['descent_depth']:.4f} Ha "
                  f"({bare_diag['descent_depth']/scr_diag['descent_depth']:.1f}x reduction)")
            if scr_pk_diag:
                print(f"    + W1c + PK     : {scr_pk_diag['descent_depth']:.4f} Ha "
                      f"({scr_diag['descent_depth']/scr_pk_diag['descent_depth']:.2f}x add'l)")
            print(f"    + W1c+PK+SV    : {full_diag['descent_depth']:.4f} Ha "
                  f"({(scr_pk_diag['descent_depth'] if scr_pk_diag else scr_diag['descent_depth'])/full_diag['descent_depth']:.2f}x add'l)")
            if sv_diag:
                print(f"\n  Pure SV (no W1c, no PK): {sv_diag['descent_depth']:.4f} Ha")
            print(f"\n  Comparison to physical D_e: factor {full_diag['descent_depth'] / D_e_exp:.1f}x deeper")

    out = {
        'R_grid': R_grid,
        'R_eq_exp': R_eq_exp,
        'D_e_exp': D_e_exp,
        'configs': [
            {'label': label, 'screened': scr, 'pk': pk, 'valence': val}
            for label, scr, pk, val in configs
        ],
        'pts': all_pts,
        'diag': all_diag,
        'fit': all_fit,
        'bound_anywhere': bool(bound_anywhere),
    }
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\n[saved] {out_path}", flush=True)


if __name__ == '__main__':
    main()
