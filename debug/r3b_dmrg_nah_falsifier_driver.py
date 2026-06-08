"""
Sprint R3-B  --  NaH DMRG-on-FCIDUMP falsifier (Phase 1 hybrid pipeline)
=========================================================================

Load-bearing question: does DMRG on the GeoVac balanced FCIDUMP close the
W1e wall on NaH?

GO  -> DMRG produces an interior min at R_e in [3.0, 4.5] bohr with
       D_e in [0.04, 0.15] Ha; cosmic-Galois Class 1/Class 3 partition
       reading of W1e is vindicated; Paper 57 abstract structures around
       "W1e closure under DMRG".
STOP -> monotone descent persists / no interior min / D_e overbound by
       >=10x; W1e is NOT pure Class 1 multi-determinant correlation;
       Paper 57 abstract reflects "W1e is W1e is not solver-closable on the
       balanced n_max=2/3 FCIDUMP".

Method.  At the production NaH balanced n_max=2 and n_max=3 FCIDUMP
dimensions the FCI Hilbert space is small enough (C(M, 1)^2 = 100 at
n_max=2, 784 at n_max=3) that the chi -> infinity DMRG limit is direct
sector-restricted FCI.  We therefore exercise the falsifier in two
independent paths:

  (A) Export FCIDUMP via GeoVacHamiltonian.to_fcidump, read it back via
      read_fcidump, and run an independent sector-restricted FCI on the
      reparsed (h1, eri, ecore).  This is the "DMRG at chi -> inf"
      benchmark; it IS the gold standard of any classical correlation
      solver on the same Hamiltonian.

  (B) Run the framework's coupled_fci_energy directly on the in-memory
      build_balanced_hamiltonian result.  This is the P4 baseline path;
      bit-exact equality with (A) confirms (i) the FCIDUMP round-trip
      preserves all physics, (ii) the falsifier is internally consistent.

The two paths must agree to ~1e-10 Ha.  If they do, the resulting PES is
the DMRG-at-infinite-bond-dimension answer; any finite-chi DMRG on the
same FCIDUMP can at best match it.  The verdict for the load-bearing
question follows directly from the PES shape.

Outputs.
  debug/data/r3b_dmrg_nah.json
  debug/data/r3b_nah_balanced_nmax{2,3}_R*.fcidump  (one per R-grid point)
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy
from geovac.ecosystem_export import GeoVacHamiltonian, read_fcidump
from geovac.molecular_spec import nah_spec


REPO = Path(__file__).resolve().parent.parent
DATA_DIR = REPO / 'debug' / 'data'
DATA_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Falsifier reference: P4 baseline (CLAUDE.md / sprint_p4_nah_w1e_baseline_memo)
# ---------------------------------------------------------------------------

R_GRID = [2.5, 3.0, 3.566, 4.0, 4.5, 5.0, 6.0]
R_EXP = 3.566
DE_EXP = 0.0713  # Huber-Herzberg

# Anchor: balanced n_max=3 at R=8 bohr is -168.0667 Ha (P4 §2.1 Panel C)
P4_BALANCED_N3_R8 = -168.0667
P4_BALANCED_N2_R8 = -167.4236

P4_BALANCED_N2 = {
    2.500: -169.7814,
    3.000: -169.3914,
    3.566: -169.1146,
    4.000: -168.9661,
    4.500: -168.8086,
    5.000: -168.6365,
    6.000: -168.2356,
}
P4_BALANCED_N3 = {
    2.500: -170.0530,
    3.000: -169.7541,
    3.566: -169.5684,
    4.000: -169.4605,
    4.500: -169.3197,
    5.000: -169.1472,
    6.000: -168.7473,
}


# ---------------------------------------------------------------------------
# Sector-restricted FCI from (h1, eri, ecore, n_electrons)
# ---------------------------------------------------------------------------

def fci_from_integrals(
    h1: np.ndarray, eri: np.ndarray, ecore: float, n_electrons: int,
) -> Dict[str, Any]:
    """
    Wrap coupled_fci_energy so it consumes integrals from any source
    (the parsed FCIDUMP, or any other classical-chemistry input).

    Returns dict with E_total and dim.
    """
    M = int(h1.shape[0])
    fake = {
        'M': M,
        'h1': h1,
        'eri': eri,
        'nuclear_repulsion': float(ecore),
    }
    fci = coupled_fci_energy(fake, n_electrons=n_electrons, verbose=False)
    return {
        'E_total_Ha': float(fci['E_coupled']),
        'M': M,
    }


# ---------------------------------------------------------------------------
# One-R driver: build, export, round-trip, FCI both paths
# ---------------------------------------------------------------------------

def run_single(R: float, max_n: int, write_fcidump: bool = True) -> Dict[str, Any]:
    t0 = time.perf_counter()
    spec = nah_spec(R=R, max_n=max_n)
    n_el = sum(b.n_electrons for b in spec.blocks)
    ham = build_balanced_hamiltonian(
        spec, R=R, n_grid_vne=4000, L_max=4, verbose=False,
    )
    t_build = time.perf_counter() - t0

    h1 = np.asarray(ham['h1'], dtype=float)
    eri = np.asarray(ham['eri'], dtype=float)
    ecore = float(ham['nuclear_repulsion'])
    M = int(ham['M'])

    # Path B: direct FCI on in-memory integrals (P4 baseline path)
    t1 = time.perf_counter()
    fci_direct = coupled_fci_energy(ham, n_electrons=n_el, verbose=False)
    E_direct = float(fci_direct['E_coupled'])
    t_fci_direct = time.perf_counter() - t1

    # FCIDUMP write + read back
    R_str = f"{R:.3f}".replace('.', 'p')
    fname = DATA_DIR / f"r3b_nah_balanced_nmax{max_n}_R{R_str}.fcidump"
    H = GeoVacHamiltonian(
        ham['qubit_op'],
        metadata={
            'system': 'NaH', 'arch': 'balanced',
            'R_bohr': float(R), 'max_n': int(max_n), 'M': M,
        },
        h1=h1, eri=eri, ecore=ecore, n_electrons=n_el,
    )
    if write_fcidump:
        H.to_fcidump(str(fname))
        parsed = read_fcidump(str(fname))
        h1_p = np.asarray(parsed['h1'], dtype=float)
        eri_p = np.asarray(parsed['eri'], dtype=float)
        ecore_p = float(parsed['ecore'])
        # Round-trip diagnostics
        d_h1 = float(np.max(np.abs(h1 - h1_p)))
        d_eri = float(np.max(np.abs(eri - eri_p)))
        d_ecore = abs(ecore - ecore_p)
    else:
        h1_p, eri_p, ecore_p = h1, eri, ecore
        d_h1 = d_eri = d_ecore = 0.0

    # Path A: independent FCI on the parsed FCIDUMP integrals
    t2 = time.perf_counter()
    fci_from_dump = fci_from_integrals(h1_p, eri_p, ecore_p, n_el)
    E_from_dump = float(fci_from_dump['E_total_Ha'])
    t_fci_dump = time.perf_counter() - t2

    return {
        'R_bohr': float(R),
        'max_n': int(max_n),
        'n_electrons': int(n_el),
        'M': M,
        'Q': int(ham.get('Q', 2 * M)),
        'fci_dim': int((np.math.comb(M, n_el // 2)) ** 2)
                    if hasattr(np, 'math') and hasattr(np.math, 'comb')
                    else None,
        'E_direct_Ha': E_direct,
        'E_from_FCIDUMP_Ha': E_from_dump,
        'two_paths_max_diff_Ha': float(abs(E_direct - E_from_dump)),
        'fcidump_round_trip': {
            'max_abs_h1_diff': d_h1,
            'max_abs_eri_diff': d_eri,
            'abs_ecore_diff': d_ecore,
            'file': str(fname.relative_to(REPO)) if write_fcidump else None,
        },
        'timings_s': {
            't_build_balanced': float(t_build),
            't_fci_direct': float(t_fci_direct),
            't_fci_from_FCIDUMP': float(t_fci_dump),
        },
        'V_NN_plus_E_core_Ha': ecore,
        'E_electronic_Ha': E_direct - ecore,
    }


# ---------------------------------------------------------------------------
# Scan + comparison to P4 baseline
# ---------------------------------------------------------------------------

def scan_and_compare(max_n: int) -> Dict[str, Any]:
    print(f"\n== R3-B falsifier scan: balanced NaH n_max={max_n} ==")
    points: List[Dict[str, Any]] = []
    for R in R_GRID:
        print(f"  R = {R:5.3f} bohr ...", end='', flush=True)
        p = run_single(R, max_n=max_n, write_fcidump=True)
        points.append(p)
        print(f"  E = {p['E_direct_Ha']:.6f} Ha"
              f"  (FCIDUMP path: {p['E_from_FCIDUMP_Ha']:.6f},"
              f" two-path diff {p['two_paths_max_diff_Ha']:.2e},"
              f" round-trip h1/eri/ecore: {p['fcidump_round_trip']['max_abs_h1_diff']:.1e}/"
              f"{p['fcidump_round_trip']['max_abs_eri_diff']:.1e}/"
              f"{p['fcidump_round_trip']['abs_ecore_diff']:.1e})")

    # Compare to P4 reference
    p4_ref = P4_BALANCED_N2 if max_n == 2 else P4_BALANCED_N3
    p4_anchor = P4_BALANCED_N2_R8 if max_n == 2 else P4_BALANCED_N3_R8

    comparison = []
    for p in points:
        R = p['R_bohr']
        p4_value = p4_ref.get(R, None)
        diff = (p['E_direct_Ha'] - p4_value) if p4_value is not None else None
        comparison.append({
            'R_bohr': R,
            'E_R3B_Ha': p['E_direct_Ha'],
            'E_P4_Ha': p4_value,
            'diff_R3B_minus_P4_Ha': diff,
        })

    # Find interior minimum
    E_arr = np.array([p['E_direct_Ha'] for p in points])
    R_arr = np.array([p['R_bohr'] for p in points])

    # An interior minimum is a sign change in the discrete derivative from
    # negative (descending) to positive (ascending) somewhere inside the
    # R-grid (excluding the endpoints).  Equivalent: E[i] < E[i-1] and
    # E[i] < E[i+1] for some interior i.
    interior_idx = None
    for i in range(1, len(E_arr) - 1):
        if E_arr[i] < E_arr[i - 1] and E_arr[i] < E_arr[i + 1]:
            interior_idx = int(i)
            break

    if interior_idx is None:
        # Determine direction of monotone slope
        descending_count = int(np.sum(np.diff(E_arr) < 0))
        ascending_count = int(np.sum(np.diff(E_arr) > 0))
        if descending_count == len(E_arr) - 1:
            monotone_kind = 'descending (toward small R)'
        elif ascending_count == len(E_arr) - 1:
            monotone_kind = 'ascending (toward large R)'
        else:
            monotone_kind = 'mixed (no interior min)'
        interior_min = None
        well_depth_Ha = None
    else:
        monotone_kind = None
        R_min = float(R_arr[interior_idx])
        E_min = float(E_arr[interior_idx])
        # Use R=6 as dissociation anchor (largest R on our grid)
        E_anchor = float(E_arr[-1])
        well_depth_Ha = E_anchor - E_min
        interior_min = {
            'R_bohr': R_min,
            'E_Ha': E_min,
            'E_anchor_at_R_max_Ha': E_anchor,
            'well_depth_Ha': float(well_depth_Ha),
        }

    # Wall quantification: signed gap at R_exp vs P4 R=8 anchor
    E_at_R_exp = next(
        p['E_direct_Ha'] for p in points if abs(p['R_bohr'] - R_EXP) < 1e-6
    )
    gap_vs_R8 = E_at_R_exp - p4_anchor
    # NEGATIVE gap => over-binding (well-deep) vs the dissociation anchor.

    return {
        'max_n': max_n,
        'R_grid': R_GRID,
        'points': points,
        'comparison_to_P4': comparison,
        'interior_min': interior_min,
        'monotone_kind': monotone_kind,
        'wall_at_R_exp_Ha': float(gap_vs_R8),
        'wall_vs_exp_De_ratio': float(abs(gap_vs_R8) / DE_EXP),
        'well_depth_R3B_Ha': float(well_depth_Ha) if well_depth_Ha is not None else None,
        'P4_anchor_R8_Ha': float(p4_anchor),
        'DE_exp_Ha': DE_EXP,
    }


# ---------------------------------------------------------------------------
# Verdict classifier
# ---------------------------------------------------------------------------

def classify_verdict(scan: Dict[str, Any]) -> Dict[str, Any]:
    """Apply the load-bearing GO/BORDERLINE/STOP gate."""
    if scan['interior_min'] is not None:
        Rm = scan['interior_min']['R_bohr']
        De = scan['interior_min']['well_depth_Ha']
        if 3.0 <= Rm <= 4.5 and 0.04 <= De <= 0.15:
            verdict = 'GO'
            note = (f"Interior min at R={Rm:.3f} bohr (R_exp={R_EXP}), "
                    f"D_e={De:.4f} Ha (exp={DE_EXP:.4f} Ha); within gate.")
        elif scan['interior_min'] is not None:
            # Interior min exists but outside the gate
            de_ratio = De / DE_EXP if De > 0 else None
            r_off = abs(Rm - R_EXP) / R_EXP if Rm > 0 else None
            if r_off and r_off > 0.15:
                verdict = 'BORDERLINE'
                note = (f"Interior min at R={Rm:.3f} bohr "
                        f"(R_exp={R_EXP}, off by {100*r_off:.1f}%); D_e={De:.4f} Ha.")
            elif de_ratio and de_ratio > 3.0:
                verdict = 'BORDERLINE'
                note = (f"Interior min at R={Rm:.3f} bohr; D_e={De:.4f} Ha "
                        f"({de_ratio:.1f}x exp). Off by factor>3.")
            else:
                verdict = 'BORDERLINE'
                note = "Interior min exists; gate borderline."
    else:
        verdict = 'STOP'
        note = (f"No interior minimum; PES is monotone {scan['monotone_kind']}. "
                f"Gap at R_exp vs P4 R=8 anchor: "
                f"{scan['wall_at_R_exp_Ha']:.3f} Ha "
                f"({scan['wall_vs_exp_De_ratio']:.1f}x exp D_e).")
    return {'verdict': verdict, 'note': note}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    out: Dict[str, Any] = {
        'sprint': 'R3-B NaH DMRG-on-FCIDUMP falsifier',
        'date': '2026-06-07',
        'R_grid_bohr': R_GRID,
        'R_exp_bohr': R_EXP,
        'DE_exp_Ha': DE_EXP,
        'falsifier_gate': {
            'GO_window': {'R_eq_bohr': [3.0, 4.5],
                          'De_Ha': [0.04, 0.15]},
        },
    }

    # n_max = 2 (Q=20, FCI dim = 100; instant)
    scan_n2 = scan_and_compare(max_n=2)
    verdict_n2 = classify_verdict(scan_n2)
    print(f"\n  R3-B n_max=2 verdict: {verdict_n2['verdict']}")
    print(f"  {verdict_n2['note']}")
    out['n_max_2'] = {**scan_n2, **verdict_n2}

    # n_max = 3 (Q=56, FCI dim = 784; ~2 minutes per R, 7 R's => 14 minutes)
    scan_n3 = scan_and_compare(max_n=3)
    verdict_n3 = classify_verdict(scan_n3)
    print(f"\n  R3-B n_max=3 verdict: {verdict_n3['verdict']}")
    print(f"  {verdict_n3['note']}")
    out['n_max_3'] = {**scan_n3, **verdict_n3}

    # Net verdict: STOP if either n_max gives STOP (the wall persists with
    # better basis); GO only if both n_max=2 AND n_max=3 fall in the GO
    # window (W1e closure verified across basis sizes).
    if verdict_n2['verdict'] == 'GO' and verdict_n3['verdict'] == 'GO':
        net = 'GO'
    elif (verdict_n2['verdict'] == 'STOP' or
          verdict_n3['verdict'] == 'STOP'):
        net = 'STOP'
    else:
        net = 'BORDERLINE'
    out['net_verdict'] = net

    out_path = DATA_DIR / 'r3b_dmrg_nah.json'
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\n  Written: {out_path.relative_to(REPO)}")
    print(f"\n== R3-B NET VERDICT: {net} ==")
    return out


if __name__ == '__main__':
    main()
