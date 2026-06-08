"""Sprint SO(4)-Breaking — diagnostic of the W1e chemistry-binding wall.

Hypothesis (from end-of-session conversation 2026-06-07): the chemistry-
side W1e wall might be a SPECIFIC PATTERN of SO(4)-symmetry breaking
that the second-quantized integral export misses.  The framework's
per-sub-block h1 starts with EXACT SO(4) (hydrogenic eigenvalues
−Z²/(2n²) with 2s/2p degenerate at each Z).  PK + cross-center V_ne
break this SO(4) within and across blocks.  The continuous Level-4
PK-composed solver binds LiH at 5.3% R_eq; the integral-export
path doesn't bind.  What gets lost in projection might be a specific
class of SO(4)-breaking operators.

Diagnostic procedure:
  1. For LiH and NaH at composed and balanced builders, sweep R
  2. At each R, extract per-sub-block h1 sub-matrix, diagonalize
  3. Print eigenvalues by (n, l) — verify 2s/2p splittings
  4. Quantify SO(4)-breaking metric per sub-block:
       gap_l = max_l(E_l(sub_block)) - min_l(E_l(sub_block))
     and l-mixing entropy:
       S_l = -Σ |⟨nl|ψ_eig⟩|² log|⟨nl|ψ_eig⟩|²
  5. Compare to NIST atomic eigenvalues (Li_2s = -0.198, Na_3s = -0.189,
     H_1s = -0.5)
  6. Sweep R and look for correlation between SO(4)-breaking pattern
     and E_FCI (does the wall depth correlate with a specific
     breaking-pattern misalignment?)

Decision gate:
  - POSITIVE: SO(4)-breaking pattern strongly correlates with W1e
    wall depth at specific R, identifying a candidate missing operator
  - NEGATIVE: no correlation; W1e is structurally below SO(4)-level
    diagnostics (sharpens the wall statement)
  - PARTIAL: some correlation, identifies one or two diagnostic axes
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import scipy.linalg as la

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, Exception):
    pass

from geovac.molecular_spec import lih_spec, nah_spec
from geovac.balanced_coupled import (
    build_balanced_hamiltonian,
    _get_block_geometry,
    _get_nuclei_for_lih, _get_nuclei_for_nah,
)
from geovac.composed_qubit import build_composed_hamiltonian
from geovac.coupled_composition import coupled_fci_energy


# NIST atomic reference data (Hartree)
# Source: NIST CCCBDB / Slater–Condon parameterization of atomic levels
NIST_REFERENCE = {
    # Hydrogenic (exact)
    'H_1s':   -0.500000,
    'H_2s':   -0.125000,
    'H_2p':   -0.125000,

    # Lithium atom (1s² 2s¹ ground state)
    'Li_1s':  -2.491,    # ~ -1/2 * (Z_eff^2) with Z_eff ≈ 2.69 (Slater)
    'Li_2s':  -0.198,    # NIST: IE = 5.392 eV → 0.198 Ha
    'Li_2p':  -0.130,    # NIST: 2p excitation level ~ -0.130 Ha

    # Sodium atom (1s² 2s² 2p⁶ 3s¹ ground state)
    'Na_2s':  -3.073,
    'Na_2p':  -1.521,
    'Na_3s':  -0.189,    # NIST: IE = 5.139 eV → 0.189 Ha
    'Na_3p':  -0.112,    # NIST: 3p ~ -0.112 Ha
}

# Experimental bond dissociation energies (Hartree)
NIST_BINDING = {
    'LiH': {'R_eq_bohr': 3.015, 'D_e_ha': 0.0924},
    'NaH': {'R_eq_bohr': 3.566, 'D_e_ha': 0.0713},
}


def label_orbital(n: int, l: int, n_val_offset: int = 0) -> str:
    """(n, l) → 'ns', 'np', 'nd', ... with physical_n = block_n + n_val_offset."""
    physical_n = n + n_val_offset
    l_letters = ['s', 'p', 'd', 'f', 'g']
    return f"{physical_n}{l_letters[l] if l < len(l_letters) else f'l{l}'}"


def analyze_sub_block_h1(
    h1_full: np.ndarray, sub_block: Dict, label_prefix: str = '',
    n_val_offset: int = 0,
) -> Dict[str, Any]:
    """Extract the sub-block's h1 sub-matrix, diagonalize, and group
    eigenvalues by approximate (n, l) labels using maximum-overlap
    with the basis (n, l, m) orbitals.
    """
    states = sub_block['states']
    off = sub_block['offset']
    n_orb = len(states)
    idx = list(range(off, off + n_orb))
    h_sub = h1_full[np.ix_(idx, idx)]
    eigvals, eigvecs = la.eigh(h_sub)

    # Group eigenvalues by (n, l) using maximum-overlap label
    grouped: Dict[Tuple[int, int], List[float]] = {}
    grouped_diag: Dict[Tuple[int, int], List[float]] = {}
    for k, ev in enumerate(eigvals):
        # Find maximum-overlap (n, l, m) label
        v = eigvecs[:, k]
        weights = v ** 2
        best_i = int(np.argmax(weights))
        n, l, m = states[best_i]
        # Mixing entropy
        eps = 1e-15
        H_mix = -np.sum(weights * np.log(weights + eps))
        grouped.setdefault((n, l), []).append(ev)
        grouped_diag.setdefault((n, l), []).append({
            'eigenvalue': float(ev),
            'best_orbital_idx': best_i,
            'best_orbital': (int(n), int(l), int(m)),
            'mixing_entropy': float(H_mix),
            'weights_on_native': [float(w) for w in weights],
        })

    # Diagonal h1 elements per (n, l) — what the bare framework starts with
    diag_h1 = {}
    for i, (n, l, m) in enumerate(states):
        diag_h1.setdefault((n, l), []).append(float(h_sub[i, i]))

    # Per-block SO(4)-breaking metric: max(E_l) - min(E_l) within the block
    so4_gap = float(np.max(eigvals) - np.min(eigvals))
    # 2s vs 2p (or block analog): compare lowest-l and second-l shells
    nl_keys = sorted(grouped.keys())
    s_eigs = grouped.get((nl_keys[1][0] if len(nl_keys) > 1 else 1, 0), [])  # 2s
    p_eigs = grouped.get((nl_keys[1][0] if len(nl_keys) > 1 else 1, 1), [])  # 2p
    sp_gap = (
        float(min(p_eigs) - min(s_eigs))
        if (s_eigs and p_eigs) else None
    )

    return {
        'label': label_prefix,
        'n_val_offset': n_val_offset,
        'eigenvalues': eigvals.tolist(),
        'grouped_by_nl': {f"{label_orbital(n, l, n_val_offset)}":
                         eigs for (n, l), eigs in grouped.items()},
        'diagonal_by_nl': {f"{label_orbital(n, l, n_val_offset)}":
                          diag for (n, l), diag in diag_h1.items()},
        'so4_block_gap': so4_gap,
        '2s_2p_gap_or_analog': sp_gap,
        'diag_breakdown_per_orbital': {
            f"{k[0]}_{k[1]}": v for k, v in grouped_diag.items()
        },
    }


def analyze_at_R(
    name: str, spec, R: float, builder: str = 'balanced',
) -> Dict[str, Any]:
    """Analyze SO(4)-breaking pattern at a single R."""
    if builder == 'balanced':
        if name == 'LiH':
            nuclei = _get_nuclei_for_lih(spec, R)
        elif name == 'NaH':
            nuclei = _get_nuclei_for_nah(R)
        else:
            nuclei = None
        result = build_balanced_hamiltonian(
            spec, R=R, nuclei=nuclei, verbose=False,
        )
    elif builder == 'composed':
        result = build_composed_hamiltonian(spec, verbose=False)
    else:
        raise ValueError(f"Unknown builder: {builder}")

    h1 = result['h1']
    M = result['M']
    nuc_rep = result['nuclear_repulsion']

    sub_blocks = _get_block_geometry(spec)
    block_records: List[Dict[str, Any]] = []
    for sb in sub_blocks:
        # Determine n_val_offset for this sub-block
        b_idx = sb['parent_block']
        parent_block = spec.blocks[b_idx]
        n_val_offset = getattr(parent_block, 'n_val_offset', 0)
        analysis = analyze_sub_block_h1(
            h1, sb, label_prefix=sb['label'],
            n_val_offset=n_val_offset,
        )
        block_records.append(analysis)

    # Also run FCI to get current binding behavior
    fake = {'M': M, 'h1': h1, 'eri': result['eri'],
            'nuclear_repulsion': nuc_rep}
    try:
        n_el = 4 if name == 'LiH' else 2  # NaH balanced uses 2-electron
        fci_out = coupled_fci_energy(fake, n_electrons=n_el, verbose=False)
        E_fci = fci_out['E_coupled']
    except Exception as e:
        E_fci = None

    return {
        'R_bohr': R,
        'builder': builder,
        'nuclear_repulsion': float(nuc_rep),
        'E_FCI': E_fci,
        'sub_blocks': block_records,
    }


def print_analysis(name: str, results: List[Dict[str, Any]]):
    print(f"\n{'='*78}")
    print(f"System: {name}")
    print(f"{'='*78}")
    for r in results:
        print(f"\n--- R = {r['R_bohr']} bohr, builder = {r['builder']} ---")
        print(f"  E_FCI = {r['E_FCI']:.6f} Ha  "
              f"(V_NN+E_core = {r['nuclear_repulsion']:.4f})")
        for sb in r['sub_blocks']:
            print(f"  Sub-block: {sb['label']} "
                  f"(n_val_offset = {sb['n_val_offset']})")
            print(f"    diagonal h1 per (n,l):")
            for nl, diag in sb['diagonal_by_nl'].items():
                values_str = ', '.join(f"{v:+.5f}" for v in diag)
                print(f"      {nl:5s}: [{values_str}]")
            print(f"    eigenvalues (grouped by max-overlap (n,l)):")
            for nl, eigs in sb['grouped_by_nl'].items():
                values_str = ', '.join(f"{v:+.5f}" for v in eigs)
                ref = NIST_REFERENCE.get(f'{name[:-1] if name != "H" else "H"}_{nl}', None)
                ref_str = f"  NIST: {ref:+.5f}" if ref is not None else ''
                print(f"      {nl:5s}: [{values_str}]{ref_str}")
            print(f"    SO(4) block gap: {sb['so4_block_gap']:.4f} Ha")
            if sb['2s_2p_gap_or_analog'] is not None:
                print(f"    s/p analog gap:  {sb['2s_2p_gap_or_analog']:+.4f} Ha")


def main():
    print("="*78)
    print("Sprint SO(4)-Breaking — W1e Chemistry-Binding Diagnostic")
    print("Date: 2026-06-07 (session continuation)")
    print("="*78)

    # LiH analysis
    print("\nBuilding LiH at composed and balanced builders...")
    spec_lih = lih_spec()
    R_panel_lih = [2.5, 3.015, 3.5, 4.0, 5.0, 8.0]
    results_lih_balanced = []
    for R in R_panel_lih:
        t0 = time.time()
        r = analyze_at_R('LiH', spec_lih, R, builder='balanced')
        results_lih_balanced.append(r)
        print(f"  LiH balanced at R={R}: {time.time()-t0:.1f}s")

    # NaH analysis
    print("\nBuilding NaH at balanced builder...")
    spec_nah = nah_spec()
    R_panel_nah = [3.0, 3.566, 4.0, 5.0, 8.0]
    results_nah_balanced = []
    for R in R_panel_nah:
        t0 = time.time()
        try:
            r = analyze_at_R('NaH', spec_nah, R, builder='balanced')
            results_nah_balanced.append(r)
            print(f"  NaH balanced at R={R}: {time.time()-t0:.1f}s")
        except Exception as e:
            print(f"  NaH balanced at R={R}: FAILED: {e}")

    print_analysis('LiH', results_lih_balanced)
    print_analysis('NaH', results_nah_balanced)

    # SO(4)-breaking pattern vs R sweep
    print("\n" + "="*78)
    print("SO(4)-BREAKING vs R (LiH balanced)")
    print("="*78)
    print(f"{'R (bohr)':>10s}  {'E_FCI':>12s}", end='')
    sub_labels = [sb['label'] for sb in results_lih_balanced[0]['sub_blocks']]
    for sl in sub_labels:
        print(f"  {sl + '_SO4gap':>20s}", end='')
    print()
    for r in results_lih_balanced:
        print(f"{r['R_bohr']:>10.3f}  {r['E_FCI']:>12.6f}", end='')
        for sb in r['sub_blocks']:
            print(f"  {sb['so4_block_gap']:>20.5f}", end='')
        print()

    print("\n" + "="*78)
    print("SO(4)-BREAKING vs R (NaH balanced)")
    print("="*78)
    if results_nah_balanced:
        print(f"{'R (bohr)':>10s}  {'E_FCI':>12s}", end='')
        sub_labels = [sb['label'] for sb in results_nah_balanced[0]['sub_blocks']]
        for sl in sub_labels:
            print(f"  {sl + '_SO4gap':>20s}", end='')
        print()
        for r in results_nah_balanced:
            print(f"{r['R_bohr']:>10.3f}  {r['E_FCI']:>12.6f}", end='')
            for sb in r['sub_blocks']:
                print(f"  {sb['so4_block_gap']:>20.5f}", end='')
            print()

    # Save
    out_data = {
        'date': '2026-06-07',
        'system': 'SO(4)-breaking W1e diagnostic',
        'lih_balanced': results_lih_balanced,
        'nah_balanced': results_nah_balanced,
        'nist_reference': NIST_REFERENCE,
        'nist_binding': NIST_BINDING,
    }

    def _jsonable(obj):
        if isinstance(obj, (np.integer, np.floating)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, dict):
            return {k: _jsonable(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [_jsonable(v) for v in obj]
        return obj

    out_path = Path(__file__).parent / 'data' / 'sprint_so4_breaking_w1e.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(_jsonable(out_data), f, indent=2)
    print(f"\nData written to {out_path}")


if __name__ == '__main__':
    main()
