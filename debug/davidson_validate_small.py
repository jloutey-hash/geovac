"""
Validation leg 1: the matrix-free sigma vs DENSE builds on tiny random systems.

Two independent references:
  (a) brute-force Fock-space FCI (explicit a^dag / a on bitmasks)   -> physics
  (b) verbatim transcription of coupled_fci_energy's assembly       -> library
Checks:
  * DirectCI4e(faithful=False) sigma reproduces (a) to machine precision
  * DirectCI4e(faithful=True)  sigma reproduces (b) to machine precision
  * ground energies agree in both modes
"""
from __future__ import annotations
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from davidson_probe0_signcheck import random_integrals, brute_force_fci
from davidson_probe0b_localize import lib_dense
from davidson_ci import DirectCI4e


def dense_from_sigma(ci: DirectCI4e) -> np.ndarray:
    nd = ci.ndet
    H = np.empty((nd, nd))
    for k in range(nd):
        e = np.zeros(nd)
        e[k] = 1.0
        H[:, k] = ci.sigma(e.reshape(ci.na, ci.nb)).reshape(-1)
    return H


if __name__ == '__main__':
    print("=" * 92)
    print("VALIDATION 1 -- matrix-free sigma vs dense references (random integrals)")
    print("=" * 92)
    ok = True
    for M in (4, 5, 6, 7):
        h1, eri = random_integrals(M, seed=M * 17 + 3)
        e_core = 0.4321
        Hb, _ = brute_force_fci(h1, eri, M, 2, 2, e_core)
        Hl, _, _ = lib_dense(h1, eri, M, 4, e_core, same_spin_double_sign=+1.0)

        ci_c = DirectCI4e(h1, eri, e_core, faithful=False)
        ci_f = DirectCI4e(h1, eri, e_core, faithful=True)
        Hc = dense_from_sigma(ci_c)
        Hf = dense_from_sigma(ci_f)

        d_bf = np.abs(Hc - Hb).max()
        d_lib = np.abs(Hf - Hl).max()
        # diagonal (preconditioner) check
        d_diag_c = np.abs(np.diag(Hc) - ci_c.diag.reshape(-1)).max()
        d_diag_f = np.abs(np.diag(Hf) - ci_f.diag.reshape(-1)).max()
        e_c = np.linalg.eigvalsh(Hc)[0]
        e_f = np.linalg.eigvalsh(Hf)[0]
        e_b = np.linalg.eigvalsh(Hb)[0]
        e_l = np.linalg.eigvalsh(Hl)[0]
        scale = max(np.abs(Hb).max(), 1.0)
        good = (d_bf / scale < 1e-12 and d_lib / scale < 1e-12
                and d_diag_c / scale < 1e-12 and d_diag_f / scale < 1e-12)
        ok &= good
        print(f"M={M} n_det={Hb.shape[0]:5d}  "
              f"|Hcorr-Hbrute|={d_bf:.2e}  |Hfaith-Hlib|={d_lib:.2e}  "
              f"|diag|={max(d_diag_c,d_diag_f):.2e}  "
              f"dE_corr={e_c-e_b:+.2e}  dE_faith={e_f-e_l:+.2e}  "
              f"{'OK' if good else 'FAIL'}")

    # Davidson vs dense eigh
    print("\n" + "=" * 92)
    print("VALIDATION 2 -- Davidson vs dense eigh (same random systems)")
    print("=" * 92)
    for M in (6, 8):
        h1, eri = random_integrals(M, seed=M * 5)
        e_core = -1.25
        Hb, _ = brute_force_fci(h1, eri, M, 2, 2, e_core)
        e_exact = np.linalg.eigvalsh(Hb)[0]
        ci = DirectCI4e(h1, eri, e_core, faithful=False)
        out = ci.ground_state(tol=1e-9, verbose=False, max_sub=12)
        print(f"M={M} n_det={ci.ndet:5d}  E_dav={out['E']:+.12f}  E_eigh={e_exact:+.12f}  "
              f"diff={out['E']-e_exact:+.2e}  |r|={out['residual']:.1e}  "
              f"iters={out['n_iter']}  {'OK' if abs(out['E']-e_exact) < 1e-9 else 'FAIL'}")
        ok &= abs(out['E'] - e_exact) < 1e-9

    print("\nOVERALL:", "PASS" if ok else "FAIL")
