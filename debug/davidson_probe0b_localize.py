"""Localize the coupled_fci_energy vs brute-force discrepancy element-by-element."""
from __future__ import annotations
import itertools
import numpy as np
from scipy.sparse import lil_matrix
from geovac.coupled_composition import _excitation_phase, _double_excitation_phase
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from davidson_probe0_signcheck import random_integrals, brute_force_fci


def lib_dense(h1, eri, M, n_e, e_core, same_spin_double_sign=+1.0):
    """Verbatim transcription of coupled_fci_energy's matrix assembly, dense,
    with a toggle on the same-spin double-excitation sign."""
    from itertools import combinations
    n_up = n_down = n_e // 2
    A = list(combinations(range(M), n_up)); B = list(combinations(range(M), n_down))
    na, nb = len(A), len(B); nd = na * nb
    ai_idx = {s: i for i, s in enumerate(A)}; bi_idx = {s: i for i, s in enumerate(B)}
    H = np.zeros((nd, nd))
    di = lambda a, b: a * nb + b
    for ai, al in enumerate(A):
        for bi, be in enumerate(B):
            I = di(ai, bi); E = e_core
            for p in al: E += h1[p, p]
            for p in be: E += h1[p, p]
            for i in range(n_up):
                for j in range(i + 1, n_up):
                    p, q = al[i], al[j]; E += eri[p, p, q, q] - eri[p, q, q, p]
            for i in range(n_down):
                for j in range(i + 1, n_down):
                    p, q = be[i], be[j]; E += eri[p, p, q, q] - eri[p, q, q, p]
            for p in al:
                for q in be: E += eri[p, p, q, q]
            H[I, I] = E
    # alpha singles
    for ai, al in enumerate(A):
        aset = set(al)
        for p in al:
            for r in range(M):
                if r in aset: continue
                na_ = tuple(sorted((aset - {p}) | {r}))
                if na_ not in ai_idx: continue
                ain = ai_idx[na_]; ph = _excitation_phase(al, p, r)
                for bi, be in enumerate(B):
                    val = ph * h1[r, p]
                    for q in al:
                        if q == p: continue
                        val += ph * (eri[r, p, q, q] - eri[r, q, q, p])
                    for q in be: val += ph * eri[r, p, q, q]
                    H[di(ai, bi), di(ain, bi)] += val
    # beta singles
    for bi, be in enumerate(B):
        bset = set(be)
        for p in be:
            for r in range(M):
                if r in bset: continue
                nb_ = tuple(sorted((bset - {p}) | {r}))
                if nb_ not in bi_idx: continue
                bin_ = bi_idx[nb_]; ph = _excitation_phase(be, p, r)
                for ai, al in enumerate(A):
                    val = ph * h1[r, p]
                    for q in be:
                        if q == p: continue
                        val += ph * (eri[r, p, q, q] - eri[r, q, q, p])
                    for q in al: val += ph * eri[r, p, q, q]
                    H[di(ai, bi), di(ai, bin_)] += val
    # alpha doubles
    for ai, al in enumerate(A):
        aset = set(al); occ = list(al)
        for i1 in range(n_up):
            for i2 in range(i1 + 1, n_up):
                p, q = occ[i1], occ[i2]
                for r in range(M):
                    if r in aset: continue
                    for s in range(r + 1, M):
                        if s in aset: continue
                        na_ = tuple(sorted((aset - {p, q}) | {r, s}))
                        if na_ not in ai_idx: continue
                        ain = ai_idx[na_]
                        ph = _double_excitation_phase(al, p, q, r, s) * same_spin_double_sign
                        val = ph * (eri[r, p, s, q] - eri[r, q, s, p])
                        for bi in range(nb): H[di(ai, bi), di(ain, bi)] += val
    # beta doubles
    for bi, be in enumerate(B):
        bset = set(be); occ = list(be)
        for i1 in range(n_down):
            for i2 in range(i1 + 1, n_down):
                p, q = occ[i1], occ[i2]
                for r in range(M):
                    if r in bset: continue
                    for s in range(r + 1, M):
                        if s in bset: continue
                        nb_ = tuple(sorted((bset - {p, q}) | {r, s}))
                        if nb_ not in bi_idx: continue
                        bin_ = bi_idx[nb_]
                        ph = _double_excitation_phase(be, p, q, r, s) * same_spin_double_sign
                        val = ph * (eri[r, p, s, q] - eri[r, q, s, p])
                        for ai in range(na): H[di(ai, bi), di(ai, bin_)] += val
    # alpha-beta doubles
    for ai, al in enumerate(A):
        aset = set(al)
        for pa in al:
            for ra in range(M):
                if ra in aset: continue
                na_ = tuple(sorted((aset - {pa}) | {ra}))
                if na_ not in ai_idx: continue
                ain = ai_idx[na_]; pha = _excitation_phase(al, pa, ra)
                for bi, be in enumerate(B):
                    bset = set(be)
                    for pb in be:
                        for rb in range(M):
                            if rb in bset: continue
                            nb_ = tuple(sorted((bset - {pb}) | {rb}))
                            if nb_ not in bi_idx: continue
                            bin_ = bi_idx[nb_]; phb = _excitation_phase(be, pb, rb)
                            H[di(ai, bi), di(ain, bin_)] += pha * phb * eri[ra, pa, rb, pb]
    return H, A, B


if __name__ == '__main__':
    M, n_e = 4, 4
    h1, eri = random_integrals(M, seed=M * 10 + 2)
    e_core = 0.37
    Hb, dets = brute_force_fci(h1, eri, M, 2, 2, e_core)
    for sgn in (+1.0, -1.0):
        Hl, A, B = lib_dense(h1, eri, M, n_e, e_core, same_spin_double_sign=sgn)
        d = np.abs(Hl - Hb)
        print(f"same_spin_double_sign={sgn:+.0f}:  max|H_lib - H_brute| = {d.max():.3e}  "
              f"nnz_diff={(d>1e-10).sum()}  E0_lib={np.linalg.eigvalsh(Hl)[0]:+.12f} "
              f" E0_brute={np.linalg.eigvalsh(Hb)[0]:+.12f}")
        if d.max() > 1e-10:
            ij = np.argwhere(d > 1e-10)[:6]
            for i, j in ij:
                na_, nb_ = len(A), len(B)
                print(f"    I=({A[i//nb_]},{B[i%nb_]}) J=({A[j//nb_]},{B[j%nb_]})  "
                      f"lib={Hl[i,j]:+.6f} brute={Hb[i,j]:+.6f}")

# --- extra: does the flip hold at n_up=n_down=3 (spectator present)? ---
def extra():
    for (M, n_e) in [(5, 6), (6, 6)]:
        h1, eri = random_integrals(M, seed=99 + M)
        Hb, _ = brute_force_fci(h1, eri, M, n_e // 2, n_e // 2, 0.0)
        for sgn in (+1.0, -1.0):
            Hl, _, _ = lib_dense(h1, eri, M, n_e, 0.0, same_spin_double_sign=sgn)
            d = np.abs(Hl - Hb).max()
            print(f"  M={M} n_e={n_e} sign={sgn:+.0f}: max|dH|={d:.3e}")
