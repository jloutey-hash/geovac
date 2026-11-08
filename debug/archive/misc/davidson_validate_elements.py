"""
Validation leg 3: element-level agreement at LARGE M, where the library's
explicit sparse assembly is unaffordable (~2.3 h/pt at n_max=3, days at n_max=4).

Strategy: sample random determinant pairs (I, J) at the real balanced-LiH
integrals, evaluate H[I, J] with the LIBRARY's own Slater-rule expressions and
phase helpers (geovac.coupled_composition._excitation_phase /
_double_excitation_phase, term-for-term as in coupled_fci_energy), and compare
against the corresponding entry of a column of the matrix-free sigma.

Covers every excitation class the assembly builds:
  diagonal / alpha-single / beta-single / alpha-double / beta-double / ab-double
"""
from __future__ import annotations
import argparse
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import _excitation_phase, _double_excitation_phase
from geovac.molecular_spec import lih_spec
from davidson_ci import DirectCI4e


def ref_element(h1, eri, ecore, A_I, B_I, A_J, B_J):
    """H[I, J] exactly as coupled_fci_energy assembles it (library conventions)."""
    sA, sB = set(A_I), set(B_I)
    tA, tB = set(A_J), set(B_J)
    dA, dB = sA ^ tA, sB ^ tB
    na_diff, nb_diff = len(dA) // 2, len(dB) // 2
    if na_diff == 0 and nb_diff == 0:
        E = ecore
        for p in A_I:
            E += h1[p, p]
        for p in B_I:
            E += h1[p, p]
        for p, q in itertools.combinations(A_I, 2):
            E += eri[p, p, q, q] - eri[p, q, q, p]
        for p, q in itertools.combinations(B_I, 2):
            E += eri[p, p, q, q] - eri[p, q, q, p]
        for p in A_I:
            for q in B_I:
                E += eri[p, p, q, q]
        return E, 'diagonal'
    if na_diff == 1 and nb_diff == 0:
        p = (sA - tA).pop(); r = (tA - sA).pop()
        ph = _excitation_phase(A_I, p, r)
        val = ph * h1[r, p]
        for q in A_I:
            if q == p:
                continue
            val += ph * (eri[r, p, q, q] - eri[r, q, q, p])
        for q in B_I:
            val += ph * eri[r, p, q, q]
        return val, 'alpha-single'
    if na_diff == 0 and nb_diff == 1:
        p = (sB - tB).pop(); r = (tB - sB).pop()
        ph = _excitation_phase(B_I, p, r)
        val = ph * h1[r, p]
        for q in B_I:
            if q == p:
                continue
            val += ph * (eri[r, p, q, q] - eri[r, q, q, p])
        for q in A_I:
            val += ph * eri[r, p, q, q]
        return val, 'beta-single'
    if na_diff == 2 and nb_diff == 0:
        p, q = sorted(sA - tA); r, s = sorted(tA - sA)
        ph = _double_excitation_phase(A_I, p, q, r, s)
        return ph * (eri[r, p, s, q] - eri[r, q, s, p]), 'alpha-double'
    if na_diff == 0 and nb_diff == 2:
        p, q = sorted(sB - tB); r, s = sorted(tB - sB)
        ph = _double_excitation_phase(B_I, p, q, r, s)
        return ph * (eri[r, p, s, q] - eri[r, q, s, p]), 'beta-double'
    if na_diff == 1 and nb_diff == 1:
        pa = (sA - tA).pop(); ra = (tA - sA).pop()
        pb = (sB - tB).pop(); rb = (tB - sB).pop()
        ph = _excitation_phase(A_I, pa, ra) * _excitation_phase(B_I, pb, rb)
        return ph * eri[ra, pa, rb, pb], 'ab-double'
    return 0.0, 'disconnected'


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--max_n', type=int, default=3)
    ap.add_argument('--R', type=float, default=3.015)
    ap.add_argument('--n_cols', type=int, default=12)
    ap.add_argument('--seed', type=int, default=7)
    args = ap.parse_args()

    spec = lih_spec(R=args.R, max_n=args.max_n)
    ham = build_balanced_hamiltonian(spec, R=args.R, n_grid_vne=8000, L_max=4,
                                     screened_cross_center=False, verbose=False)
    M, h1, eri = ham['M'], ham['h1'], ham['eri']
    ecore = ham['nuclear_repulsion']
    strings = list(itertools.combinations(range(M), 2))
    ci = DirectCI4e(h1, eri, ecore, faithful=True, verbose=True)
    na = ci.na
    rng = np.random.default_rng(args.seed)

    print(f"\nelement-level check: M={M}  n_det={ci.ndet:,}  "
          f"{args.n_cols} random columns x all classes", flush=True)
    stats = {}
    worst = 0.0
    for _ in range(args.n_cols):
        ai = int(rng.integers(na)); bi = int(rng.integers(na))
        e = np.zeros(ci.ndet); e[ai * na + bi] = 1.0
        col = ci.sigma(e.reshape(na, na)).reshape(-1)      # col[K] = H[K, I]
        A_I, B_I = strings[ai], strings[bi]
        # targets: the diagonal + a random sample from every excitation class
        targets = [(ai, bi)]
        for _ in range(40):
            aj = int(rng.integers(na)); bj = int(rng.integers(na))
            targets.append((aj, bj))
        # force coverage of each class
        virt = [o for o in range(M) if o not in set(A_I)]
        for _ in range(12):
            p = A_I[int(rng.integers(2))]; r = int(rng.choice(virt))
            aj = strings.index(tuple(sorted((set(A_I) - {p}) | {r})))
            targets.append((aj, bi))
            vb = [o for o in range(M) if o not in set(B_I)]
            pb = B_I[int(rng.integers(2))]; rb = int(rng.choice(vb))
            bj = strings.index(tuple(sorted((set(B_I) - {pb}) | {rb})))
            targets.append((ai, bj))
            targets.append((aj, bj))
            r2, s2 = sorted(rng.choice(virt, 2, replace=False))
            aj2 = strings.index((int(r2), int(s2)))
            targets.append((aj2, bi))
            targets.append((ai, strings.index((int(r2), int(s2)))
                            if tuple(sorted((int(r2), int(s2)))) in strings else bi))
        for (aj, bj) in targets:
            ref, cls = ref_element(h1, eri, ecore, A_I, B_I, strings[aj], strings[bj])
            got = col[aj * na + bj]
            d = abs(got - ref)
            s = stats.setdefault(cls, [0, 0.0])
            s[0] += 1
            s[1] = max(s[1], d)
            worst = max(worst, d)
    print(f"\n{'class':16s} {'n':>6s}  max|H_sigma - H_library|")
    for cls, (n, d) in sorted(stats.items()):
        print(f"{cls:16s} {n:6d}  {d:.3e}")
    print(f"\nWORST OVER ALL CLASSES: {worst:.3e}   "
          f"{'PASS' if worst < 1e-12 else 'FAIL'}")
