"""group3 QA diagnostic — pair-diagonal vs global-M_L ERI selection rule impact.

READ-ONLY. No production code is edited; the angular-coupling coefficient is
monkeypatched at runtime only. Restores the original after each measurement.

Mechanism (verified in debug/qa/group3_followup_density_layers_memo.md):
  geovac/composed_qubit.py::_ck_coefficient uses q = mc - ma, so the
  3j bottom row sums to 2(mc-ma) != 0 unless ma == mc. _wigner3j returns 0
  for any nonzero m-sum, so every m_a != m_c coupling is dropped -> the
  realized composed ERI tensor is PAIR-DIAGONAL (m_a=m_c AND m_b=m_d).

  The physically-correct global Coulomb rule (casimir_ci._gaunt_ck) uses
  q = m1 - m2, making the bottom row sum to 0 identically, so genuinely
  nonzero m-swap elements (e.g. <p+1 p-1 | p-1 p+1> = 0.24) survive. The
  global m-conservation check ma+mb == mc+md is ALREADY present in
  _build_eri_block (line ~377); the only thing suppressing m-swaps is
  _ck_coefficient's per-pair q convention.

This driver:
  1. Baseline (production pair-diagonal): build LiH, BeH2, H2O the corpus way
     and compute number-projected FCI energy + Pauli count.
  2. Modified (global-M_L): monkeypatch _ck_coefficient to use the _gaunt_ck
     convention (q = ma - mc), rebuild, recompute.
  3. Sanity cross-checks: pure-angular density 1.44% (pd) -> 6.06% (global)
     at l_max=3; LiH baseline Pauli count vs published ~909 pre-taper.
  4. Per-molecule table: E_baseline, E_modified, dE; N_Pauli ratio.

Number-projected FCI ONLY (coupled_composition.coupled_fci_energy), never
raw qubit-space diagonalization (memory/feedback_tc_correction.md).
"""

from __future__ import annotations

import json
import sys
import time
from itertools import product
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, Exception):
    pass

import geovac.composed_qubit as cq
from geovac.composed_qubit import build_composed_hamiltonian, _wigner3j
from geovac.casimir_ci import _gaunt_ck
from geovac.coupled_composition import coupled_fci_energy
from geovac.molecular_spec import lih_spec, beh2_spec, h2o_spec


# ---------------------------------------------------------------------------
# The global-M_L replacement for _ck_coefficient.
# ---------------------------------------------------------------------------
# Original production (pair-diagonal): q = mc - ma  -> bottom row sums to
#   -ma + (mc-ma) + mc = 2(mc-ma), nonzero unless ma==mc.
# Global Coulomb (casimir_ci._gaunt_ck convention): q = m1 - m2 = ma - mc
#   -> bottom row sums to -ma + (ma-mc) + mc = 0 identically.
# We delegate to _gaunt_ck so the logic is literally the corpus's global path.
def _ck_coefficient_global(la: int, ma: int, lc: int, mc: int, k: int) -> float:
    """Global-M_L angular coupling, identical to casimir_ci._gaunt_ck."""
    return _gaunt_ck(la, ma, lc, mc, k)


# ---------------------------------------------------------------------------
# Sanity: pure-angular ERI density under each convention.
# ---------------------------------------------------------------------------
def _enumerate_lm(l_max: int) -> List[Tuple[int, int]]:
    """All (l, m) with 0 <= l <= l_max, -l <= m <= l (single-shell per l)."""
    return [(l, m) for l in range(l_max + 1) for m in range(-l, l + 1)]


def _angular_density(l_max: int, ck_func) -> Tuple[float, int, int]:
    """Fraction of (a,b,c,d) angular ERIs that are nonzero under ck_func.

    Counts <ab|cd> over the (l,m) index set; an element is nonzero if
    global m-conservation holds AND there is a k with both c^k(a,c) and
    c^k(b,d) nonzero (the same structure _build_eri_block uses).
    Returns (density_percent, n_nonzero, n_total).
    """
    states = _enumerate_lm(l_max)
    n = len(states)

    # c^k(a,c) table
    from collections import defaultdict
    ac_k: Dict[Tuple[int, int], List[int]] = defaultdict(list)
    for a in range(n):
        la, ma = states[a]
        for c in range(n):
            lc, mc = states[c]
            for k in range(0, la + lc + 1):
                if (la + lc + k) % 2 != 0:
                    continue
                if abs(ck_func(la, ma, lc, mc, k)) > 1e-12:
                    ac_k[(a, c)].append(k)

    n_nonzero = 0
    n_total = n ** 4
    for a, b, c, d in product(range(n), repeat=4):
        la, ma = states[a]
        lb, mb = states[b]
        lc, mc = states[c]
        ld, md = states[d]
        if ma + mb != mc + md:
            continue
        ks_ac = ac_k.get((a, c))
        ks_bd = ac_k.get((b, d))
        if not ks_ac or not ks_bd:
            continue
        if set(ks_ac) & set(ks_bd):
            n_nonzero += 1
    return 100.0 * n_nonzero / n_total, n_nonzero, n_total


def run_sanity() -> Dict:
    """Pure-angular density sanity cross-check at l_max=1..3."""
    print("=" * 72)
    print("SANITY — pure-angular ERI density (no radial / no molecule)")
    print("=" * 72)
    rows = []
    print(f"{'l_max':>6s}  {'pair-diag %':>12s}  {'global-ML %':>12s}")
    print("-" * 36)
    for l_max in (1, 2, 3):
        d_pd, nz_pd, tot = _angular_density(l_max, cq._ck_coefficient)
        d_gl, nz_gl, _ = _angular_density(l_max, _ck_coefficient_global)
        print(f"{l_max:6d}  {d_pd:12.4f}  {d_gl:12.4f}")
        rows.append({
            'l_max': l_max,
            'density_pair_diagonal_pct': d_pd,
            'density_global_ML_pct': d_gl,
            'n_nonzero_pd': nz_pd,
            'n_nonzero_global': nz_gl,
            'n_total': tot,
        })

    # The m-swap element that production drops: <p+1 p-1 | p-1 p+1>.
    # c^k(a=p+1, c=p-1) under each convention at k=2:
    pd_val = cq._ck_coefficient(1, 1, 1, -1, 2)
    gl_val = _ck_coefficient_global(1, 1, 1, -1, 2)
    print(f"\n  m-swap probe c^2(p+1, p-1):  pair-diag={pd_val:.4f}  "
          f"global={gl_val:.4f}")
    print("  (production drops it to 0.0; global retains the genuine Gaunt value)")
    return {'angular_density': rows,
            'm_swap_probe_pair_diagonal': float(pd_val),
            'm_swap_probe_global': float(gl_val)}


# ---------------------------------------------------------------------------
# Per-molecule build + number-projected FCI under both conventions.
# ---------------------------------------------------------------------------
def _total_electrons(spec) -> int:
    return int(sum(getattr(b, 'n_electrons', 0) for b in spec.blocks))


def _measure(spec_factory, label: str, use_global: bool) -> Dict:
    """Build a molecule's composed Hamiltonian and number-projected FCI.

    If use_global, monkeypatch _ck_coefficient to the global-M_L convention
    for the duration of the build, then restore.
    """
    spec = spec_factory()
    n_elec = _total_electrons(spec)

    orig = cq._ck_coefficient
    if use_global:
        cq._ck_coefficient = _ck_coefficient_global
    try:
        t0 = time.time()
        res = build_composed_hamiltonian(spec, pk_in_hamiltonian=True,
                                         verbose=False)
        build_t = time.time() - t0
        N_pauli = int(res['N_pauli'])
        Q = int(res['Q'])
        M = int(res['M'])
        n_eri_nonzero = int(np.sum(np.abs(res['eri']) > 1e-12))

        t1 = time.time()
        fci = coupled_fci_energy(res, n_electrons=n_elec, verbose=False)
        fci_t = time.time() - t1
        E = float(fci['E_coupled'])
    finally:
        cq._ck_coefficient = orig

    conv = 'global-M_L' if use_global else 'pair-diagonal'
    print(f"  [{label} / {conv:13s}]  E_FCI={E:+12.6f}  N_Pauli={N_pauli:6d}  "
          f"Q={Q}  M={M}  N_e={n_elec}  eri_nz={n_eri_nonzero}  "
          f"build={build_t:.1f}s fci={fci_t:.1f}s")
    return {
        'label': label, 'convention': conv,
        'E_FCI': E, 'N_pauli': N_pauli, 'Q': Q, 'M': M,
        'n_electrons': n_elec, 'n_eri_nonzero': n_eri_nonzero,
        'build_time_s': build_t, 'fci_time_s': fci_t,
    }


def main():
    out_dir = Path(__file__).resolve().parent
    data_json = out_dir / 'group3_density_diagnostic.json'

    sanity = run_sanity()

    # Restore guard: confirm the patch flips the rule and original is intact.
    assert abs(cq._ck_coefficient(1, 1, 1, -1, 2)) < 1e-12, \
        "production _ck_coefficient should drop the m-swap (pair-diagonal)"

    print("\n" + "=" * 72)
    print("PER-MOLECULE — baseline (pair-diagonal) vs modified (global-M_L)")
    print("=" * 72)

    # LiH and BeH2 first (clearest / cheapest); H2O last (may be heavy).
    molecules = [
        ('LiH', lih_spec),
        ('BeH2', beh2_spec),
        ('H2O', h2o_spec),
    ]

    results = []
    for label, factory in molecules:
        try:
            base = _measure(factory, label, use_global=False)
            mod = _measure(factory, label, use_global=True)
            results.append({'baseline': base, 'modified': mod})
        except MemoryError:
            print(f"  [{label}] SKIPPED — MemoryError at production params")
            results.append({'label': label, 'skipped': 'MemoryError'})
        except Exception as e:  # noqa: BLE001
            print(f"  [{label}] SKIPPED — {type(e).__name__}: {e}")
            results.append({'label': label, 'skipped': f'{type(e).__name__}: {e}'})

    # ---- Summary table ----
    print("\n" + "=" * 72)
    print("RESULTS TABLE")
    print("=" * 72)
    print(f"{'mol':>5s}  {'E_base':>12s}  {'E_mod':>12s}  {'dE(Ha)':>10s}  "
          f"{'dE(%)':>8s}  {'Np_base':>8s}  {'Np_mod':>8s}  {'ratio':>6s}")
    print("-" * 82)
    table = []
    for r in results:
        if 'skipped' in r:
            print(f"{r['label']:>5s}  {'SKIPPED: ' + r['skipped']}")
            table.append({'label': r['label'], 'skipped': r['skipped']})
            continue
        b, m = r['baseline'], r['modified']
        dE = m['E_FCI'] - b['E_FCI']
        dE_pct = 100.0 * dE / abs(b['E_FCI']) if b['E_FCI'] != 0 else float('nan')
        ratio = m['N_pauli'] / b['N_pauli'] if b['N_pauli'] else float('nan')
        print(f"{b['label']:>5s}  {b['E_FCI']:+12.6f}  {m['E_FCI']:+12.6f}  "
              f"{dE:+10.6f}  {dE_pct:+8.4f}  {b['N_pauli']:8d}  "
              f"{m['N_pauli']:8d}  {ratio:6.2f}")
        table.append({
            'label': b['label'],
            'E_baseline': b['E_FCI'], 'E_modified': m['E_FCI'],
            'dE_Ha': dE, 'dE_pct': dE_pct,
            'N_pauli_baseline': b['N_pauli'], 'N_pauli_modified': m['N_pauli'],
            'N_pauli_ratio': ratio,
            'eri_nz_baseline': b['n_eri_nonzero'],
            'eri_nz_modified': m['n_eri_nonzero'],
            'Q': b['Q'], 'M': b['M'], 'n_electrons': b['n_electrons'],
        })

    payload = {
        'sanity': sanity,
        'per_molecule': results,
        'table': table,
    }
    with open(data_json, 'w') as f:
        json.dump(payload, f, indent=2)
    print(f"\nWrote {data_json}")


if __name__ == '__main__':
    main()
