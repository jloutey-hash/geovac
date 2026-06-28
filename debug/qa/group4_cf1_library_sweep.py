"""group4 QA pre-work — CF-1 library-wide Pauli-count re-pricing sweep.

READ-ONLY: monkeypatches _ck_coefficient at runtime only (restores after each).
BUILD-ONLY: counts N_pauli under the production pair-diagonal ERI rule vs the
physical global-M_L rule for every COMPOSED library molecule. NO FCI (the energy
axis is intractable at the composed active space and is a chemistry footnote per
CF-1; the group4 concern is the SPARSITY multiplier).

Mechanism: geovac/composed_qubit._ck_coefficient uses q = mc - ma (pair-diagonal,
drops m-swaps); the physical Coulomb rule is casimir_ci._gaunt_ck (q = ma - mc,
global M_L). See debug/qa/group3_density_diagnostic.py + the CF-1 carryforward.
"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import numpy as np

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

import geovac.composed_qubit as cq
from geovac.composed_qubit import build_composed_hamiltonian
from geovac.casimir_ci import _gaunt_ck
from geovac import ecosystem_export as ee


def _ck_global(la, ma, lc, mc, k):
    return _gaunt_ck(la, ma, lc, mc, k)


# Composed library (CF-1 applies to multi-center composed ERIs).
NAMES = (
    list(ee._HYDRIDE_Z.keys())          # 18 main-group hydrides
    + list(ee._MULTI_CENTER.keys())     # 5 multi-center diatomics
    + list(ee._TM_HYDRIDE_Z.keys())     # 10 transition-metal hydrides
    + ['SrH', 'BaH']                    # 2 alkaline-earth monohydrides
)


def _npauli(name, use_global):
    spec = ee._rebuild_spec(name, R=None, max_n=2, core_method='pk')
    orig = cq._ck_coefficient
    if use_global:
        cq._ck_coefficient = _ck_global
    try:
        res = build_composed_hamiltonian(spec, pk_in_hamiltonian=True, verbose=False)
        return int(res['N_pauli']), int(res['Q']), int(res['M']), int(np.sum(np.abs(res['eri']) > 1e-12))
    finally:
        cq._ck_coefficient = orig


def main():
    out = Path(__file__).resolve().parent / 'group4_cf1_library_sweep.json'
    # sanity: production drops the m-swap; global retains it
    assert abs(cq._ck_coefficient(1, 1, 1, -1, 2)) < 1e-12
    assert abs(_ck_global(1, 1, 1, -1, 2)) > 1e-12
    print("sanity OK: pair-diagonal drops c^2(p+1,p-1); global retains it\n")
    print(f"{'mol':>6s} {'Q':>4s} {'M':>4s} {'Np_pd':>7s} {'Np_glob':>8s} {'ratio':>6s} {'eri_pd':>7s} {'eri_gl':>7s}")
    print("-" * 60)
    rows = []
    for name in NAMES:
        try:
            t0 = time.time()
            np_pd, Q, M, eri_pd = _npauli(name, False)
            np_gl, _, _, eri_gl = _npauli(name, True)
            dt = time.time() - t0
            ratio = np_gl / np_pd if np_pd else float('nan')
            print(f"{name:>6s} {Q:4d} {M:4d} {np_pd:7d} {np_gl:8d} {ratio:6.2f} {eri_pd:7d} {eri_gl:7d}  ({dt:.0f}s)")
            rows.append({'mol': name, 'Q': Q, 'M': M, 'Np_pair_diagonal': np_pd,
                         'Np_global_ML': np_gl, 'ratio': ratio,
                         'eri_nz_pd': eri_pd, 'eri_nz_global': eri_gl})
        except Exception as e:  # noqa: BLE001
            print(f"{name:>6s}  SKIPPED — {type(e).__name__}: {str(e)[:60]}")
            rows.append({'mol': name, 'skipped': f'{type(e).__name__}: {str(e)[:80]}'})
    ok = [r for r in rows if 'ratio' in r]
    if ok:
        ratios = [r['ratio'] for r in ok]
        print("-" * 60)
        print(f"n={len(ok)} built | ratio: min {min(ratios):.2f}  median {sorted(ratios)[len(ratios)//2]:.2f}  max {max(ratios):.2f}  mean {sum(ratios)/len(ratios):.2f}")
    out.write_text(json.dumps({'rows': rows}, indent=2), encoding='utf-8')
    print(f"\nWrote {out}")


if __name__ == '__main__':
    main()
