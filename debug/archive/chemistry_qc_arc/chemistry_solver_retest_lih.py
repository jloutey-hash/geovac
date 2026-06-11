"""
Chemistry-Solver Re-Test: LiH PES under multi-focal architecture closures
==========================================================================

Tests whether the multi-focal Phase C closures (balanced coupled +
cross-center V_ne via Shibuya-Wulfman + screened W1c via cross_center_screened_vne)
close, persist, or change the character of the v2.0.32 baseline R_eq drift
(+0.15-0.22 bohr / l_max in the level-4 PES solver).

Architecture-faithful interpretation (worker fork, May 2026):

  1. The v2.0.32 baseline drift was diagnosed in the level-4 PK-composed
     PES solver (composed_diatomic.py + variational_2d).
  2. The multi-focal Phase C closures live in the qubit-encoded balanced
     coupled architecture (balanced_coupled.py), which is structurally
     a DIFFERENT solver. It does not have a separate l_max parameter --
     the basis is enumerated at (n, l, m) with l < n, so n_max controls
     both radial and angular convergence simultaneously.
  3. Track CD (v2.0.39, balanced coupled LiH) reported:
       n_max = 2: E = -7.924 Ha (1.8% err), R_eq = 3.226 (7.0% err)
       n_max = 3: E = -8.055 Ha (0.20% err), R_eq = 3.280 (8.8% err)
     i.e. drift +0.053 bohr/n_max -- THREE TIMES SMALLER than PK's drift,
     but not zero. Track CD's headline was "energy converges excellently,
     R_eq drifts structurally."
  4. The current re-test asks whether the W1c closure (`screened_cross_center=True`)
     changes the n_max-drift character. For LiH (first-row, no frozen core)
     the W1c closure should be a bit-identical no-op per the cross_center_screened_vne
     module's auto-detection (Z<11 -> None). We verify that.
  5. The deeper question: with the architecture available, is there a path
     to push R_eq below 5%?

Targets (priority order):
  1. LiH balanced coupled at max_n=2, fine R-grid around 3.015 bohr
  2. LiH balanced coupled at max_n=2 with screened_cross_center=True (W1c on/off)
  3. LiH balanced coupled at max_n=3 if compute permits

We use particle-number-projected FCI (coupled_fci_energy) ONLY -- never
qubit-space diag (CLAUDE.md Section 3 hard rule).
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy
from geovac.molecular_spec import lih_spec


REFERENCE = {
    'LiH': {
        'R_eq': 3.015,         # bohr (experimental)
        'E_min': -8.071,       # Ha (Coulomb-only, FCI/CBS)
        'D_e': 0.092,          # Ha (binding)
    }
}


def fit_parabolic_min(R_grid: np.ndarray, E_grid: np.ndarray) -> Tuple[float, float, float]:
    """Quadratic fit around the minimum: returns (R_min, E_min, omega_e^2/au).

    Uses the three lowest-energy points for a local quadratic fit.
    """
    i_min = int(np.argmin(E_grid))
    if i_min == 0 or i_min == len(R_grid) - 1:
        # Minimum at the edge; return raw min
        return float(R_grid[i_min]), float(E_grid[i_min]), float('nan')

    # Three-point quadratic
    R3 = R_grid[i_min - 1: i_min + 2]
    E3 = E_grid[i_min - 1: i_min + 2]
    p = np.polyfit(R3, E3, 2)
    R_eq = -p[1] / (2 * p[0])
    E_eq = p[0] * R_eq ** 2 + p[1] * R_eq + p[2]
    return float(R_eq), float(E_eq), float(2 * p[0])


def run_single(
    spec_factory,
    R: float,
    max_n: int,
    screened: bool,
    label: str,
    n_grid_vne: int = 4000,
    L_max: int = 4,
) -> Dict[str, Any]:
    """Run a single (R, max_n, screened) point through balanced FCI."""
    t0 = time.perf_counter()
    spec = spec_factory(R=R, max_n=max_n)
    n_electrons = sum(b.n_electrons for b in spec.blocks)
    ham = build_balanced_hamiltonian(
        spec,
        R=R,
        n_grid_vne=n_grid_vne,
        L_max=L_max,
        screened_cross_center=screened,
        verbose=False,
    )
    t_build = time.perf_counter() - t0

    t0 = time.perf_counter()
    fci = coupled_fci_energy(ham, n_electrons=n_electrons, verbose=False)
    t_fci = time.perf_counter() - t0

    return {
        'label': label,
        'R': R,
        'max_n': max_n,
        'screened': screened,
        'M': int(ham['M']),
        'Q': int(ham['Q']),
        'n_electrons': int(n_electrons),
        'E_total': float(fci['E_coupled']),
        'E_electronic': float(fci['E_coupled'] - ham['nuclear_repulsion']),
        'V_NN_plus_E_core': float(ham['nuclear_repulsion']),
        'cross_block_eri_count': int(ham.get('cross_block_eri_count', 0)),
        'cross_vne_count': int(ham.get('cross_vne_count', 0)),
        't_build_s': float(t_build),
        't_fci_s': float(t_fci),
    }


def run_pes(
    R_grid: List[float],
    max_n: int,
    screened: bool,
    spec_factory=lih_spec,
    label: str = 'LiH',
) -> Dict[str, Any]:
    """Run a PES scan."""
    points = []
    for i, R in enumerate(R_grid):
        print(f"  [{label} max_n={max_n} screened={screened}] "
              f"R={R:.3f} ({i+1}/{len(R_grid)})", flush=True)
        try:
            pt = run_single(spec_factory, R, max_n, screened, label)
            print(f"    E={pt['E_total']:.6f} Ha "
                  f"(t_build={pt['t_build_s']:.1f}s, t_fci={pt['t_fci_s']:.1f}s)",
                  flush=True)
        except Exception as e:
            print(f"    FAIL: {type(e).__name__}: {e}", flush=True)
            pt = {
                'label': label, 'R': R, 'max_n': max_n,
                'screened': screened, 'error': str(e),
            }
        points.append(pt)

    R_arr = np.array([p['R'] for p in points if 'E_total' in p])
    E_arr = np.array([p['E_total'] for p in points if 'E_total' in p])

    if len(R_arr) >= 3:
        R_eq, E_eq, k = fit_parabolic_min(R_arr, E_arr)
    else:
        R_eq, E_eq, k = float('nan'), float('nan'), float('nan')

    return {
        'label': label,
        'max_n': max_n,
        'screened': screened,
        'R_grid': [float(R) for R in R_grid],
        'points': points,
        'R_eq_fit': R_eq,
        'E_min_fit': E_eq,
        'curvature': k,
    }


def compute_E_core_offset(Z: int = 3) -> float:
    """The E_core offset added in hydride_spec for first-row systems.

    Per geovac/molecular_spec.py: E_core = _FIRST_ROW_CORE_ENERGY[Z] for explicit core.
    For Z=3 (Li), the He-like fallback gives -(3 - 5/16)^2 ~ -7.18.
    """
    from geovac.molecular_spec import _FIRST_ROW_CORE_ENERGY
    return _FIRST_ROW_CORE_ENERGY.get(Z, -(Z - 5.0/16) ** 2)


def main():
    import sys
    sys.stdout.reconfigure(line_buffering=True)

    out_dir = Path('debug/data')
    out_dir.mkdir(parents=True, exist_ok=True)

    R_eq_exp = REFERENCE['LiH']['R_eq']
    E_min_exp = REFERENCE['LiH']['E_min']
    E_core_off = compute_E_core_offset(3)
    print(f"E_core offset (Li): {E_core_off:.4f} Ha (subtracted to match Track CD convention)")

    R_grid_coarse = [2.4, 2.7, 3.0, 3.3, 3.6]
    R_grid_fine = [2.6, 2.85, 3.0, 3.15, 3.4]

    results: Dict[str, Any] = {
        'reference': REFERENCE['LiH'],
        'baseline_track_cd': {
            'n_max=2': {'E_min': -7.924, 'R_eq': 3.226, 'err_E_pct': 1.8, 'err_R_pct': 7.0},
            'n_max=3': {'E_min': -8.055, 'R_eq': 3.280, 'err_E_pct': 0.20, 'err_R_pct': 8.8},
            'drift_R_per_nmax': 0.053,
        },
        'baseline_pk_composed_v2_0_32': {
            'note': 'Level-4 PK-composed solver, NOT balanced coupled',
            'l_max=2': {'R_eq_err_pct': 5.3},
            'drift_R_per_lmax': 0.18,
        },
        'panels': {},
    }

    # Panel 1: LiH max_n=2, screened OFF (Track CD baseline reproduction)
    print("\n=== Panel 1: LiH max_n=2, screened=False (baseline) ===")
    p1 = run_pes(R_grid_coarse, max_n=2, screened=False, label='LiH')
    results['panels']['nmax2_screened_off'] = p1

    # Panel 2: LiH max_n=2, screened ON (W1c flag, expected no-op for first-row)
    print("\n=== Panel 2: LiH max_n=2, screened=True (W1c flag, no-op for first-row) ===")
    p2 = run_pes(R_grid_coarse, max_n=2, screened=True, label='LiH')
    results['panels']['nmax2_screened_on'] = p2

    # Verify W1c is bit-identical for LiH
    bit_identical = True
    max_diff = 0.0
    for pt1, pt2 in zip(p1['points'], p2['points']):
        if 'E_total' in pt1 and 'E_total' in pt2:
            diff = abs(pt1['E_total'] - pt2['E_total'])
            max_diff = max(max_diff, diff)
            if diff > 1e-10:
                bit_identical = False
    print(f"\n  W1c bit-identical for LiH: {bit_identical} (max E diff: {max_diff:.2e} Ha)")
    results['w1c_bit_identical_lih'] = bit_identical
    results['w1c_max_diff_lih'] = float(max_diff)

    # Panel 3: LiH max_n=3 with W1c on -- key test
    print("\n=== Panel 3: LiH max_n=3, screened=True ===")
    # First, time a single point at R_eq
    print("  Timing test: max_n=3 at R=3.015...", flush=True)
    t0 = time.perf_counter()
    try:
        pt_single = run_single(lih_spec, 3.015, 3, True, 'LiH')
        elapsed = time.perf_counter() - t0
        print(f"  Single point: t={elapsed:.1f}s, E={pt_single['E_total']:.6f} Ha "
              f"(M={pt_single['M']}, Q={pt_single['Q']})", flush=True)

        # Decide on grid based on per-point time
        if elapsed < 60:
            R_grid_3 = R_grid_coarse  # full 5-point scan
        elif elapsed < 180:
            R_grid_3 = [2.7, 3.0, 3.3]  # 3 points
        else:
            R_grid_3 = [3.0]  # just the equilibrium point

        if len(R_grid_3) >= 3:
            p3 = run_pes(R_grid_3, max_n=3, screened=True, label='LiH')
        else:
            # Single point only -- record it
            p3 = {
                'label': 'LiH',
                'max_n': 3,
                'screened': True,
                'R_grid': R_grid_3,
                'points': [pt_single],
                'R_eq_fit': float('nan'),
                'E_min_fit': float('nan'),
                'curvature': float('nan'),
                'note': 'Single-point only due to compute time',
            }
        results['panels']['nmax3_screened_on'] = p3

    except Exception as e:
        print(f"  max_n=3 FAILED: {type(e).__name__}: {e}", flush=True)
        results['panels']['nmax3_screened_on'] = {
            'label': 'LiH', 'max_n': 3, 'error': str(e),
        }

    # ----- Synthesis -----
    print("\n=== SYNTHESIS ===")
    p1_R_eq = results['panels']['nmax2_screened_off'].get('R_eq_fit', float('nan'))
    p1_E = results['panels']['nmax2_screened_off'].get('E_min_fit', float('nan'))
    p2_R_eq = results['panels']['nmax2_screened_on'].get('R_eq_fit', float('nan'))
    p2_E = results['panels']['nmax2_screened_on'].get('E_min_fit', float('nan'))

    if not (np.isnan(p1_R_eq) or np.isnan(p1_E)):
        # E_total includes E_core offset; subtract to match Track CD convention
        p1_E_tc = p1_E - E_core_off
        err_R_p1 = abs(p1_R_eq - R_eq_exp) / R_eq_exp * 100
        err_E_p1 = abs(p1_E_tc - E_min_exp) / abs(E_min_exp) * 100
        print(f"\n  LiH max_n=2 (screened OFF): R_eq={p1_R_eq:.3f} bohr "
              f"({err_R_p1:.2f}%), E_min(TC-conv)={p1_E_tc:.6f} Ha ({err_E_p1:.3f}%)")
        results['summary'] = {
            'lih_nmax2_off': {'R_eq': p1_R_eq, 'R_err_pct': err_R_p1,
                              'E_min': p1_E, 'E_min_track_cd_conv': p1_E_tc,
                              'E_err_pct': err_E_p1},
        }

    if not (np.isnan(p2_R_eq) or np.isnan(p2_E)):
        p2_E_tc = p2_E - E_core_off
        err_R_p2 = abs(p2_R_eq - R_eq_exp) / R_eq_exp * 100
        err_E_p2 = abs(p2_E_tc - E_min_exp) / abs(E_min_exp) * 100
        print(f"  LiH max_n=2 (screened ON):  R_eq={p2_R_eq:.3f} bohr "
              f"({err_R_p2:.2f}%), E_min(TC-conv)={p2_E_tc:.6f} Ha ({err_E_p2:.3f}%)")
        results.setdefault('summary', {})['lih_nmax2_on'] = {
            'R_eq': p2_R_eq, 'R_err_pct': err_R_p2,
            'E_min': p2_E, 'E_min_track_cd_conv': p2_E_tc,
            'E_err_pct': err_E_p2,
        }

    p3 = results['panels'].get('nmax3_screened_on', {})
    if 'R_eq_fit' in p3 and not np.isnan(p3['R_eq_fit']):
        p3_E_tc = p3['E_min_fit'] - E_core_off
        err_R_p3 = abs(p3['R_eq_fit'] - R_eq_exp) / R_eq_exp * 100
        err_E_p3 = abs(p3_E_tc - E_min_exp) / abs(E_min_exp) * 100
        print(f"  LiH max_n=3 (screened ON):  R_eq={p3['R_eq_fit']:.3f} bohr "
              f"({err_R_p3:.2f}%), E_min(TC-conv)={p3_E_tc:.6f} Ha ({err_E_p3:.3f}%)")
        results.setdefault('summary', {})['lih_nmax3_on'] = {
            'R_eq': p3['R_eq_fit'], 'R_err_pct': err_R_p3,
            'E_min': p3['E_min_fit'], 'E_min_track_cd_conv': p3_E_tc,
            'E_err_pct': err_E_p3,
        }

        if 'lih_nmax2_on' in results['summary']:
            drift = p3['R_eq_fit'] - results['summary']['lih_nmax2_on']['R_eq']
            print(f"\n  R_eq drift n_max=2 -> 3 (screened ON): "
                  f"{drift:+.4f} bohr/n_max")
            results['summary']['drift_nmax2_to_3'] = float(drift)
            results['summary']['drift_track_cd_baseline'] = 0.053

            # Verdict
            if abs(drift) < 0.05:
                verdict = ('CLOSED -- multi-focal architecture closed the '
                           'n_max drift below 0.05 bohr threshold')
            elif abs(drift) < 0.07:
                verdict = ('PARTIAL -- drift comparable to Track CD baseline '
                           '+0.053 bohr/n_max; W1c is a no-op for LiH so '
                           'this is the expected result')
            else:
                verdict = ('PERSISTS -- drift exceeds Track CD baseline; '
                           'investigate')
            results['summary']['verdict'] = verdict
            print(f"\n  VERDICT: {verdict}")

    out_path = out_dir / 'chemistry_solver_retest_lih.json'
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n[saved] {out_path}")

    return results


if __name__ == '__main__':
    main()
