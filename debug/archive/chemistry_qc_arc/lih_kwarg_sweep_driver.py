"""LiH kwarg sweep — does any existing kwarg reduce balanced LiH over-binding?

Baseline (Sprint W1e-Projection-Audit, 2026-06-07):
  - Continuous L4 + PK: D_e = 0.067 Ha at R_eq = 3.015 bohr (reference 0.092)
  - Balanced qubit FCI: D_e = 0.158 Ha at R_eq ~ 3.015 (2.4x over-binding)

Decision gate: any kwarg combination that brings D_e into (0.05, 0.10) Ha,
i.e., within +/- 50% of the continuous reference.

We test each kwarg individually + pairs that the Explorer report flagged as
plausibly synergistic.  Each config: PES on R panel {2.5, 3.015, 3.5, 4.0, 5.0}.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, Exception):
    pass

from geovac.molecular_spec import lih_spec
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy


R_PANEL: Tuple[float, ...] = (2.5, 3.015, 3.5, 4.0, 5.0)

# Continuous reference (from W1e-Projection-Audit, 2026-06-07)
CONTINUOUS_E_MIN = -8.1831
CONTINUOUS_E_DISSOC = -8.1160
CONTINUOUS_D_E = CONTINUOUS_E_DISSOC - CONTINUOUS_E_MIN  # -0.0671 (well depth, negative = bound)

# Gate window (D_e into +- 50% of continuous reference)
D_E_TARGET_LO = -0.10
D_E_TARGET_HI = -0.05


# Each config: (label, kwargs).  Order matters: we want to see baseline first.
CONFIGS: List[Tuple[str, Dict[str, Any]]] = [
    ('baseline',                  {}),
    ('multi_zeta',                {'multi_zeta_basis': True}),
    ('screened_valence',          {'screened_valence_basis': True}),
    ('screened_cross_center',     {'screened_cross_center': True}),
    ('pk_cross_center',           {'pk_cross_center': True}),
    ('cross_block_h1',            {'cross_block_h1': True}),
    ('mz_plus_sv',                {'multi_zeta_basis': True,
                                   'screened_valence_basis': True}),
    ('mz_plus_screen_cc',         {'multi_zeta_basis': True,
                                   'screened_cross_center': True}),
    ('sv_plus_screen_cc',         {'screened_valence_basis': True,
                                   'screened_cross_center': True}),
    ('cb_h1_plus_screen_cc',      {'cross_block_h1': True,
                                   'screened_cross_center': True}),
    ('all_four',                  {'multi_zeta_basis': True,
                                   'screened_valence_basis': True,
                                   'screened_cross_center': True,
                                   'pk_cross_center': True}),
]


def run_one_config(label: str, kwargs: Dict[str, Any]) -> Dict[str, Any]:
    """Run balanced LiH PES under the given kwargs.

    Uses the corrected post-fix Path B convention (spec.R = R, builder R = R),
    which agrees with Path A bit-identically after today's fix.
    """
    print(f"\n=== {label} ===")
    print(f"    kwargs: {kwargs}")

    energies = []
    wall = 0.0
    for R in R_PANEL:
        spec = lih_spec(R=R)
        t0 = time.time()
        try:
            res = build_balanced_hamiltonian(spec, nuclei=None, R=R,
                                             verbose=False, **kwargs)
            E = coupled_fci_energy(res, n_electrons=4,
                                   verbose=False)['E_coupled']
        except Exception as e:
            print(f"    R={R:5.3f}  FAILED: {e!r}")
            energies.append(float('nan'))
            wall += time.time() - t0
            continue
        E = float(E)
        energies.append(E)
        wall += time.time() - t0
        print(f"    R={R:5.3f}  E={E:+.6f}  t={time.time() - t0:.1f}s")

    E_arr = np.array(energies)
    if np.all(np.isnan(E_arr)):
        return {
            'label': label, 'kwargs': kwargs, 'failed': True,
            'energies': [float(x) for x in E_arr],
        }

    # Find minimum.  If at panel boundary (idx 0 or len-1), it's monotone
    # descending - no interior bowl.
    valid_mask = ~np.isnan(E_arr)
    R_valid = np.array([R_PANEL[i] for i in range(len(R_PANEL))
                        if valid_mask[i]])
    E_valid = E_arr[valid_mask]
    if len(E_valid) < 3:
        return {
            'label': label, 'kwargs': kwargs, 'failed': True,
            'energies': [float(x) for x in E_arr], 'wall_s': wall,
        }

    i_min = int(np.argmin(E_valid))
    R_min = float(R_valid[i_min])
    E_min = float(E_valid[i_min])
    E_dissoc = float(E_valid[-1])
    # D_e in our sign convention: well depth = E_dissoc - E_min;
    # negative number => bound (E_min is below E_dissoc).
    D_e = E_min - E_dissoc

    has_interior_min = (i_min not in (0, len(E_valid) - 1))

    # Gate
    if has_interior_min and D_E_TARGET_LO <= D_e <= D_E_TARGET_HI:
        verdict = 'GO'
    elif has_interior_min and D_e < D_E_TARGET_LO:
        verdict = 'BORDERLINE_OVERBOUND'
    elif has_interior_min:
        verdict = 'BORDERLINE_UNDERBOUND'
    else:
        verdict = 'STOP_NO_BOWL'

    print(f"    -> R_min={R_min:.3f}  E_min={E_min:+.4f}  D_e={D_e:+.4f}"
          f"  verdict={verdict}  wall={wall:.1f}s")

    return {
        'label': label,
        'kwargs': kwargs,
        'R_panel': list(R_PANEL),
        'energies': [float(x) for x in E_arr],
        'R_min': R_min,
        'E_min': E_min,
        'E_dissoc': E_dissoc,
        'D_e': D_e,
        'has_interior_min': has_interior_min,
        'verdict': verdict,
        'wall_s': wall,
    }


def main() -> None:
    debug_dir = Path(__file__).resolve().parent
    data_dir = debug_dir / 'data'
    data_dir.mkdir(exist_ok=True)
    out_json = data_dir / 'lih_kwarg_sweep.json'

    print("=" * 72)
    print("LiH kwarg sweep -- can existing kwargs reduce 2.4x over-binding?")
    print("=" * 72)
    print(f"Continuous reference:  D_e = {CONTINUOUS_D_E:+.4f} Ha"
          f" at R_eq = 3.015 bohr")
    print(f"Gate window:  D_e in ({D_E_TARGET_LO:+.4f}, {D_E_TARGET_HI:+.4f})"
          " Ha (+/- 50% of continuous reference)")
    print(f"R panel: {R_PANEL}")
    print(f"# configs: {len(CONFIGS)}")

    results = []
    t_total = time.time()
    for label, kwargs in CONFIGS:
        try:
            res = run_one_config(label, kwargs)
        except Exception as e:
            print(f"\n!!! {label} EXCEPTION: {e!r}")
            res = {'label': label, 'kwargs': kwargs, 'failed': True,
                   'exception': repr(e)}
        results.append(res)

    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    print(f"{'label':25s}  {'R_min':>6s}  {'D_e':>8s}  verdict")
    print("-" * 72)
    for r in results:
        if r.get('failed'):
            print(f"{r['label']:25s}  FAILED")
            continue
        print(f"{r['label']:25s}  {r['R_min']:6.3f}  {r['D_e']:+8.4f}"
              f"  {r['verdict']}")

    print(f"\nTotal wall time: {(time.time() - t_total) / 60:.1f} min")

    payload = {
        'continuous_reference': {
            'E_min': CONTINUOUS_E_MIN,
            'E_dissoc': CONTINUOUS_E_DISSOC,
            'D_e': CONTINUOUS_D_E,
            'R_eq': 3.015,
        },
        'gate_window': [D_E_TARGET_LO, D_E_TARGET_HI],
        'configs': results,
    }
    with open(out_json, 'w') as f:
        json.dump(payload, f, indent=2)
    print(f"\nWrote {out_json}")


if __name__ == '__main__':
    main()
