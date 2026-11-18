"""OP52 derivation memo, item 3: which two pieces does the full (generic-sum) operator fail to annihilate, and are they
linked to another piece by a polynomial (or first-order differential) factor over Q(x)?

For the non-annihilated pieces b, search for a piece a and polynomials r0, r1 (degree <= DEG) with
    b = r0(x) a + r1(x) T a,   T = (1 - x^2) d/dx,
as a linear system on Legendre coefficients mod p (x P_L and T P_L three-term rules, from the guard's apply_T).
Output: debug/data/op52_deriv_twopieces.json
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import flint
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location('g', ROOT / 'tests' / 'test_paper59_operator_exponents.py')
g = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g)
DEG = 4
NUSE = 300


def basis_seqs(seq, p):
    """Legendre coefficients of x^j a and x^j T a, j <= DEG, truncated to NUSE."""
    out = []
    for kind in (0, 1):
        for j in range(DEG + 1):
            a = [flint.nmod_poly([0] * j + [1], p)]
            if kind:
                a = [flint.nmod_poly([], p), flint.nmod_poly([0] * j + [1], p)]
            out.append(g.apply_T(a, seq, p)[:NUSE])
    return out


def solve(cols, rhs, p):
    """least-index solution of sum_i c_i cols[i] = rhs mod p, or None."""
    M = [[int(cols[i][L]) % p for i in range(len(cols))] + [int(rhs[L]) % p] for L in range(len(rhs))]
    n = len(cols)
    piv, r = [], 0
    for col in range(n):
        pr = next((i for i in range(r, len(M)) if M[i][col]), None)
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        inv = pow(M[r][col], p - 2, p)
        M[r] = [v * inv % p for v in M[r]]
        for i in range(len(M)):
            if i != r and M[i][col]:
                f = M[i][col]
                M[i] = [(vi - f * vr) % p for vi, vr in zip(M[i], M[r])]
        piv.append(col)
        r += 1
    if any(M[i][n] for i in range(r, len(M))):
        return None
    sol = [0] * n
    for i, c in enumerate(piv):
        sol[c] = M[i][n]
    return sol


def main():
    out = {}
    for case, p in g.CP:
        e = g.embedded(case, p)
        full = g.polys(e['ops']['full|2026'], p)
        A = g.T_form(full, p)
        bad = []
        for k in e['keys']:
            r = g.apply_T(A, e['seqs'][k], p)
            if np.any(r):
                bad.append(k)
        rel = {}
        for b in bad:
            for a in e['keys']:
                if a == b:
                    continue
                cols = basis_seqs(e['seqs'][a], p)
                s = solve(cols, e['seqs'][b][:NUSE], p)
                if s is not None:
                    rel[b] = {'a': a, 'r0': [g.small_rational(v, p) and str(g.small_rational(v, p)) for v in s[:DEG + 1]],
                              'r1': [g.small_rational(v, p) and str(g.small_rational(v, p)) for v in s[DEG + 1:]]}
                    break
        out[f'{case}|{p}'] = {'not_annihilated': bad, 'relations': rel}
        print(case, p, json.dumps(out[f'{case}|{p}']))
    (ROOT / 'debug' / 'data' / 'op52_deriv_twopieces.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
