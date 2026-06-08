"""
Sprint P4 — NaH W1e baseline (diagnostic, falsifier baseline for hybrid pipeline)
================================================================================

Goal: quantify the W1e wall on NaH for the production composed and balanced
architectures, at n_max = 2 (Q = 20) and n_max = 3 (Q = 56, if tractable).
The numerics here are the falsifier baseline against which Phase 1
DMRG/CCSD(T)-on-FCIDUMP benchmarks must measure.

Architecture notes:

  * NaH spec (Z=11, [Ne] frozen core) encodes 1 bond block, 2 valence
    electrons. The [Ne] core sits as a screened Z_eff(r) profile via the
    FrozenCore solver, contributing its energy to
    ``spec.nuclear_repulsion_constant = V_NN(R) + V_cross(R) + E_core``.

  * Production builds: composed (block-diagonal h1 + PK off because frozen
    core has no PK by construction; see hydride_spec at molecular_spec.py)
    versus balanced (composed + cross-block ERIs + cross-center V_ne via
    Shibuya-Wulfman multipole expansion).

  * R-dependence: V_NN(R) and V_cross(R) are both R-dependent. The standard
    ``balanced_coupled.build_balanced_hamiltonian`` Pattern-E fix patches
    V_NN(R) but NOT V_cross(R) (the in-code comment at line 649 incorrectly
    claims V_cross is R-independent for frozen cores). To avoid that bug,
    we rebuild ``spec`` at every R via ``nah_spec(R=R, max_n=...)``, which
    bakes the correct V_NN(R) + V_cross(R) + E_core into
    ``nuclear_repulsion_constant``. The balanced builder's V_NN patch then
    becomes a no-op (spec_R == R).

  * FCI: particle-number-projected via ``coupled_fci_energy`` (the only
    correct path; qubit-space diag is forbidden by CLAUDE.md §3 rule).

Reference (NIST CCCBDB / Huber-Herzberg):
  R_e   = 3.566 bohr (1.8874 Å)
  D_e   = 1.94 eV ≈ 0.0713 Ha
  ω_e   ≈ 1172 cm^{-1}

Outputs:
  debug/data/p4_nah_w1e_baseline.json
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

from geovac.molecular_spec import nah_spec
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.composed_qubit import build_composed_hamiltonian
from geovac.coupled_composition import coupled_fci_energy
from geovac.neon_core import FrozenCore


# ---------------------------------------------------------------------------
# Experimental anchors (Huber-Herzberg / NIST CCCBDB)
# ---------------------------------------------------------------------------
REFERENCE = {
    'NaH': {
        'R_eq_bohr': 3.566,
        'R_eq_angstrom': 1.8874,
        'D_e_Ha': 0.0713,         # 1.94 eV
        'D_e_eV': 1.94,
        'omega_e_cm1': 1172.0,
        'source': 'NIST CCCBDB / Huber-Herzberg 1979',
    }
}


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
def fit_parabolic_min(R: np.ndarray, E: np.ndarray) -> Tuple[float, float, float]:
    """Return (R_eq, E_eq, curvature) from a 3-point fit around the min,
    or (R_argmin, E_min, nan) if the min is at the edge of the grid."""
    i = int(np.argmin(E))
    if i == 0 or i == len(R) - 1:
        return float(R[i]), float(E[i]), float('nan')
    R3 = R[i - 1: i + 2]
    E3 = E[i - 1: i + 2]
    p = np.polyfit(R3, E3, 2)
    R_eq = -p[1] / (2 * p[0])
    E_eq = p[0] * R_eq ** 2 + p[1] * R_eq + p[2]
    return float(R_eq), float(E_eq), float(2 * p[0])


def dissociation_limit_energy(Z: int = 11) -> float:
    """Coulomb-only dissociation energy: Na atom (frozen-core E + valence)
    plus separated H atom (-0.5 Ha).

    Returns
    -------
    E_inf : float
        Total energy of (Na + H) at infinite separation. We anchor this to
        the LARGE-R limit of our own PES (which by construction matches
        E_core(Na+) + E_valence(Na 3s) + E(H 1s) within framework's basis)
        rather than the experimental Na ionisation, because the framework's
        binding well must be measured RELATIVE to its own dissociation.
    """
    fc = FrozenCore(Z)
    fc.solve()
    # The framework's dissociation limit is best read off the PES at large R,
    # not predicted in advance. Returned only as a reference floor.
    return fc.energy + (-0.5)


# ---------------------------------------------------------------------------
# PES driver
# ---------------------------------------------------------------------------
def run_single(R: float, max_n: int, kind: str,
               n_grid_vne: int = 4000, L_max_balanced: int = 4) -> Dict[str, Any]:
    """One PES point.

    Parameters
    ----------
    R : float
        Internuclear distance in bohr.
    max_n : int
        Maximum principal quantum number on each block.
    kind : {'composed', 'balanced'}
        Production architecture.

    Returns
    -------
    dict with E_total, E_electronic, V_NN+E_core, M, Q, FCI dim, timings.
    """
    t0 = time.perf_counter()
    # Rebuild spec at this R so V_NN(R) + V_cross(R) + E_core is correct.
    spec = nah_spec(R=R, max_n=max_n)
    n_electrons = sum(b.n_electrons for b in spec.blocks)

    if kind == 'composed':
        # PK is N/A for frozen-core (no PK params), so pk_in_hamiltonian
        # affects only the optional h1_pk diagnostic, not h1 itself.
        ham = build_composed_hamiltonian(
            spec, pk_in_hamiltonian=False, verbose=False,
        )
    elif kind == 'balanced':
        ham = build_balanced_hamiltonian(
            spec,
            R=R,
            n_grid_vne=n_grid_vne,
            L_max=L_max_balanced,
            verbose=False,
        )
    else:
        raise ValueError(f"kind must be 'composed' or 'balanced', got {kind!r}")

    t_build = time.perf_counter() - t0

    # FCI
    t1 = time.perf_counter()
    fci = coupled_fci_energy(ham, n_electrons=n_electrons, verbose=False)
    t_fci = time.perf_counter() - t1

    M = int(ham['M'])
    Q = int(ham.get('Q', 2 * M))

    return {
        'R_bohr': float(R),
        'max_n': int(max_n),
        'kind': kind,
        'n_electrons': int(n_electrons),
        'M': M,
        'Q': Q,
        'E_total_Ha': float(fci['E_coupled']),
        'E_electronic_Ha': float(
            fci['E_coupled'] - ham['nuclear_repulsion']
        ),
        'V_NN_plus_E_core_Ha': float(ham['nuclear_repulsion']),
        't_build_s': float(t_build),
        't_fci_s': float(t_fci),
    }


def run_pes(R_grid: List[float], max_n: int, kind: str,
            **kw) -> Dict[str, Any]:
    """Run a PES scan, tolerating per-point failure."""
    points: List[Dict[str, Any]] = []
    n_R = len(R_grid)
    for i, R in enumerate(R_grid):
        print(f"  [{kind} max_n={max_n}] R={R:.3f} ({i+1}/{n_R})",
              flush=True)
        try:
            pt = run_single(R, max_n, kind, **kw)
            print(f"    E={pt['E_total_Ha']:.6f} Ha "
                  f"(M={pt['M']}, Q={pt['Q']}, "
                  f"t_build={pt['t_build_s']:.1f}s, "
                  f"t_fci={pt['t_fci_s']:.1f}s)",
                  flush=True)
        except Exception as e:
            print(f"    FAIL: {type(e).__name__}: {e}", flush=True)
            pt = {'R_bohr': float(R), 'max_n': max_n, 'kind': kind,
                  'error': f'{type(e).__name__}: {e}'}
        points.append(pt)

    valid = [p for p in points if 'E_total_Ha' in p]
    R_arr = np.array([p['R_bohr'] for p in valid])
    E_arr = np.array([p['E_total_Ha'] for p in valid])

    if len(R_arr) >= 3:
        R_eq, E_eq, k = fit_parabolic_min(R_arr, E_arr)
    else:
        R_eq, E_eq, k = float('nan'), float('nan'), float('nan')

    # Reference baseline anchors
    R_large = float(R_arr.max()) if len(R_arr) else float('nan')
    E_large = float(E_arr.max()) if len(E_arr) else float('nan')
    # Numerical D_e candidate from PES (max - min along grid, signed)
    E_min_grid = float(E_arr.min()) if len(E_arr) else float('nan')
    D_e_from_PES = (E_large - E_min_grid) if len(E_arr) else float('nan')

    return {
        'kind': kind,
        'max_n': max_n,
        'R_grid_bohr': [float(R) for R in R_grid],
        'points': points,
        'R_eq_fit_bohr': R_eq,
        'E_eq_fit_Ha': E_eq,
        'curvature_au': k,
        'E_at_largest_R_Ha': E_large,
        'R_largest_bohr': R_large,
        'E_min_on_grid_Ha': E_min_grid,
        'D_e_from_PES_Ha': D_e_from_PES,
    }


def quantify_w1e_wall(pes: Dict[str, Any], ref: Dict[str, Any]) -> Dict[str, Any]:
    """Build a wall-quantification block from a single PES scan."""
    R_eq_exp = ref['R_eq_bohr']
    D_e_exp = ref['D_e_Ha']

    valid = [p for p in pes['points'] if 'E_total_Ha' in p]
    if not valid:
        return {'status': 'no_valid_points'}

    R_arr = np.array([p['R_bohr'] for p in valid])
    E_arr = np.array([p['E_total_Ha'] for p in valid])

    # Energy at experimental R_e (linear interp)
    if R_arr.min() <= R_eq_exp <= R_arr.max():
        E_at_Re_exp = float(np.interp(R_eq_exp, R_arr, E_arr))
    else:
        E_at_Re_exp = float('nan')

    # Dissociation limit anchored to largest R on grid (assume framework's
    # PES has settled to (Na + H) at our largest R).
    E_diss_framework = float(E_arr.max())  # for monotone-descending,
    # the maximum on the grid is at large R = the dissociation reference

    # Framework D_e measured against its own dissociation
    D_e_framework_from_Re_exp = E_diss_framework - E_at_Re_exp

    # Equilibrium presence: does the PES have an interior minimum?
    i_min = int(np.argmin(E_arr))
    has_interior_min = (0 < i_min < len(R_arr) - 1)
    R_eq_framework = float(R_arr[i_min]) if has_interior_min else float('nan')
    E_eq_framework = float(E_arr[i_min])

    if has_interior_min:
        D_e_framework_signed = E_diss_framework - E_eq_framework
    else:
        # Monotone descending: PES has no equilibrium; reportthe "overshoot"
        # depth relative to experimental D_e
        D_e_framework_signed = float('nan')

    out = {
        'R_eq_exp_bohr': R_eq_exp,
        'D_e_exp_Ha': D_e_exp,
        'has_interior_minimum': bool(has_interior_min),
        'R_eq_framework_bohr': R_eq_framework,
        'E_eq_framework_Ha': E_eq_framework,
        'E_at_R_eq_exp_Ha': E_at_Re_exp,
        'E_dissociation_anchor_Ha': E_diss_framework,
        'D_e_framework_at_eq_Ha': (D_e_framework_signed
                                   if has_interior_min else 'no_eq'),
        'D_e_framework_at_R_eq_exp_Ha': float(D_e_framework_from_Re_exp),
        # Signed gap experiments would need to close:
        # gap_at_eq = D_e_framework_at_eq - D_e_exp; if positive, framework
        # OVER-binds (deeper well than nature).
        'signed_gap_at_eq_Ha': (D_e_framework_signed - D_e_exp
                                if has_interior_min else 'no_eq'),
        'signed_gap_at_R_eq_exp_Ha': float(D_e_framework_from_Re_exp - D_e_exp),
        # Ratios (signed)
        'D_e_framework_over_exp': (D_e_framework_signed / D_e_exp
                                    if has_interior_min else 'no_eq'),
    }
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    import sys
    sys.stdout.reconfigure(line_buffering=True)

    out_dir = Path('debug/data')
    out_dir.mkdir(parents=True, exist_ok=True)

    R_eq_exp = REFERENCE['NaH']['R_eq_bohr']
    D_e_exp = REFERENCE['NaH']['D_e_Ha']

    print("Sprint P4 — NaH W1e baseline diagnostic")
    print("=" * 70)
    print(f"Reference: R_eq = {R_eq_exp:.3f} bohr, "
          f"D_e = {D_e_exp:.4f} Ha "
          f"({REFERENCE['NaH']['D_e_eV']:.2f} eV)")

    fc = FrozenCore(11)
    fc.solve()
    print(f"[Ne] frozen core energy E_core(Na+) = {fc.energy:.4f} Ha")
    print()

    # R grid spans well below R_eq (inner repulsion region) through to
    # dissociation. Coarse enough to be feasible at max_n=3.
    R_grid_main = [2.0, 2.5, 3.0, 3.566, 4.0, 4.5, 5.0, 6.0, 8.0]

    results: Dict[str, Any] = {
        'sprint': 'P4_NaH_W1e_baseline',
        'date': '2026-06-07',
        'reference': REFERENCE['NaH'],
        'frozen_core_E': float(fc.energy),
        'R_grid_main_bohr': R_grid_main,
        'panels': {},
        'wall_quantification': {},
        'notes': [
            "PES built by rebuilding nah_spec at each R; bakes correct "
            "V_NN(R) + V_cross(R) + E_core in nuclear_repulsion_constant.",
            "FrameworkD_e measured against PES at largest R (8.0 bohr) "
            "as the empirical dissociation anchor.",
            "FCI via coupled_fci_energy (particle-number-projected); "
            "qubit-space diag forbidden by CLAUDE.md §3.",
            "Composed and balanced are the two production architectures; "
            "balanced adds cross-block ERIs + cross-center V_ne.",
        ],
    }

    # -------- Panel A: composed, n_max=2 ----------------
    print("=== Panel A: COMPOSED, max_n=2 (Q=20) ===")
    pA = run_pes(R_grid_main, max_n=2, kind='composed')
    results['panels']['composed_nmax2'] = pA
    results['wall_quantification']['composed_nmax2'] = \
        quantify_w1e_wall(pA, REFERENCE['NaH'])
    print()

    # -------- Panel B: balanced, n_max=2 ----------------
    print("=== Panel B: BALANCED, max_n=2 (Q=20) ===")
    pB = run_pes(R_grid_main, max_n=2, kind='balanced')
    results['panels']['balanced_nmax2'] = pB
    results['wall_quantification']['balanced_nmax2'] = \
        quantify_w1e_wall(pB, REFERENCE['NaH'])
    print()

    # -------- Panel C: balanced, n_max=3 (if feasible) ----------------
    # Single-point timing test at R=R_eq_exp first.
    print("=== Panel C: BALANCED, max_n=3 (Q=56) ===")
    print("  Timing test: balanced max_n=3 at R=3.566...", flush=True)
    t0 = time.perf_counter()
    try:
        pt = run_single(3.566, max_n=3, kind='balanced')
        elapsed = time.perf_counter() - t0
        print(f"  Single point: t={elapsed:.1f}s, "
              f"E={pt['E_total_Ha']:.6f} Ha "
              f"(M={pt['M']}, Q={pt['Q']})", flush=True)

        if elapsed < 60:
            R_grid_C = R_grid_main
        elif elapsed < 180:
            R_grid_C = [2.5, 3.0, 3.566, 4.0, 5.0, 8.0]
        elif elapsed < 360:
            R_grid_C = [3.0, 3.566, 5.0, 8.0]
        else:
            R_grid_C = [3.566]
        print(f"  Selected n_max=3 grid: {R_grid_C} "
              f"({len(R_grid_C)} pts)", flush=True)

        if len(R_grid_C) >= 3:
            pC = run_pes(R_grid_C, max_n=3, kind='balanced')
        else:
            pC = {
                'kind': 'balanced', 'max_n': 3,
                'R_grid_bohr': R_grid_C, 'points': [pt],
                'R_eq_fit_bohr': float('nan'),
                'E_eq_fit_Ha': float('nan'),
                'curvature_au': float('nan'),
                'note': 'Single-point only due to compute time',
            }
        results['panels']['balanced_nmax3'] = pC
        if len(R_grid_C) >= 3:
            results['wall_quantification']['balanced_nmax3'] = \
                quantify_w1e_wall(pC, REFERENCE['NaH'])
    except Exception as e:
        print(f"  max_n=3 FAILED: {type(e).__name__}: {e}", flush=True)
        results['panels']['balanced_nmax3'] = {
            'kind': 'balanced', 'max_n': 3,
            'error': f'{type(e).__name__}: {e}',
        }
    print()

    # -------- Synthesis printout ----------------
    print("=" * 70)
    print("SYNTHESIS")
    print("=" * 70)
    for key, wq in results['wall_quantification'].items():
        if 'status' in wq:
            print(f"\n[{key}] no valid points")
            continue
        print(f"\n[{key}]")
        print(f"  has interior minimum: {wq['has_interior_minimum']}")
        if wq['has_interior_minimum']:
            print(f"  R_eq framework: {wq['R_eq_framework_bohr']:.3f} bohr "
                  f"(exp: {wq['R_eq_exp_bohr']:.3f})")
            print(f"  D_e framework at eq: "
                  f"{wq['D_e_framework_at_eq_Ha']:.4f} Ha "
                  f"(exp: {wq['D_e_exp_Ha']:.4f})")
            print(f"  ratio D_e_fw / D_e_exp: "
                  f"{wq['D_e_framework_over_exp']:.2f}×")
            print(f"  signed gap (D_e_fw - D_e_exp) at eq: "
                  f"{wq['signed_gap_at_eq_Ha']:+.4f} Ha")
        else:
            print(f"  PES is monotone descending — no interior equilibrium")
            print(f"  D_e framework at R_e^exp: "
                  f"{wq['D_e_framework_at_R_eq_exp_Ha']:+.4f} Ha "
                  f"(exp: {wq['D_e_exp_Ha']:.4f})")
            print(f"  signed gap at R_e^exp: "
                  f"{wq['signed_gap_at_R_eq_exp_Ha']:+.4f} Ha")
            print(f"  E at R_e^exp: {wq['E_at_R_eq_exp_Ha']:.4f} Ha")
            print(f"  E at largest R (anchor): "
                  f"{wq['E_dissociation_anchor_Ha']:.4f} Ha")

    out_path = out_dir / 'p4_nah_w1e_baseline.json'
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n[saved] {out_path}")

    return results


if __name__ == '__main__':
    main()
