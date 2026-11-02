"""Squashed hydrogen ring (2026-11-02, PI-directed): how fast does the ring momentum k die when a ring of N
hydrogens is squashed into an ellipse?

Molecule: N hydrogens at parametric angle phi_j = 2 pi j / N on an ellipse,
  position R ((1 + eps) cos phi_j, (1 - eps) sin phi_j, 0),  i.e.  z_j = R (e^{i phi_j} + eps e^{-i phi_j}).
R is set by the neighbour spacing s of the round ring, R = s / (2 sin(pi/N)).

Order rule (derived before any computation).  Rotating the molecule by 2 pi / N maps the configuration with complex
deformation eps to the one with eps e^{2 i (2 pi / N)}.  The Lowdin site Hamiltonian is a smooth function of
(eps, conj eps), so a term eps^{p+} conj(eps)^{p-} carries ring momentum change 2 (p+ - p-) mod N.  Hence the part
of h and of (pq|rs) that changes k by c (mod N, up to sign) is of order |eps|^{p(c)}, with
  p(c) = min |q| over integers q with 2 q = +-c (mod N).
Even N: odd c is exactly zero (the half-turn j -> j + N/2 survives); even c has p = c/2.
Odd N: even c has p = c/2 (or (N - c)/2 if smaller); odd c has p = (N - c)/2 (wrap-around).

Basis: one Slater 1s per H (zeta 1.24, 3 Gaussians), Lowdin site orbitals, complex DFT over j.
Measured: the share of the Frobenius norm of h and of (pq|rs) in each class c = |delta k| mod N.

    python debug/ellipse_ring_study.py case <N> <s_ang> <eps> [<eps> ...] --out <controls|study>
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_s = importlib.util.spec_from_file_location('hring_study', os.path.join(HERE, 'hring_study.py'))
hs = importlib.util.module_from_spec(_s)
_s.loader.exec_module(hs)
w, GE = hs.w, hs.GE
ANG = hs.ANG
ZETA_H, NG_MIN = 1.24, 3


def order_rule(N, c):
    """p(c) = min |q| with 2 q = +-c mod N; None if no such q (even N, odd c)."""
    best = None
    for q in range(-N, N + 1):
        if (2 * q - c) % N == 0 or (2 * q + c) % N == 0:
            best = abs(q) if best is None else min(best, abs(q))
    return best


def ellipse(N, s_ang, eps):
    R = s_ang * ANG / (2 * np.sin(np.pi / N))
    nuc = []
    for j in range(N):
        ph = 2 * np.pi * j / N
        nuc.append((np.array([R * (1 + eps) * np.cos(ph), R * (1 - eps) * np.sin(ph), 0.0]), 1.0))
    return nuc, R


def leakage(N, h_s, g_s):
    j = np.arange(N)
    U = np.exp(2j * np.pi * np.outer(j, j) / N) / np.sqrt(N)          # columns = k
    hc = U.conj().T @ h_s @ U
    gc = np.einsum('pqrs,pa,qb,rc,sd->abcd', g_s, U.conj(), U, U.conj(), U, optimize=True)
    k = np.arange(N)
    dk1 = (k[None, :] - k[:, None]) % N
    dk2 = (-k[:, None, None, None] + k[None, :, None, None] - k[None, None, :, None] + k[None, None, None, :]) % N
    out = {}
    for name, X, dk in (('h', hc, dk1), ('g', gc, dk2)):
        a2 = np.abs(X) ** 2
        tot = a2.sum()
        cls = np.minimum(dk, (N - dk) % N)
        for c in range(1, N // 2 + 1):
            out['%s_c%d' % (name, c)] = float(np.sqrt(a2[cls == c].sum() / tot))
    return out


def case(N, s_ang, eps, out):
    t0 = time.time()
    nuc, R = ellipse(N, s_ang, eps)
    P = np.array([C for C, _ in nuc])
    D = np.linalg.norm(P[:, None] - P[None], axis=-1)
    fns, blocks = [], []
    for C, _ in nuc:
        w.sto_shell(C, 1, 0, ZETA_H, NG_MIN, fns, blocks)
    T = w.spherical_map(fns, blocks)
    S, h, g = w.ao_integrals(fns, nuc, T)
    X = GE.lowdin_orbitals(S)
    h_s, g_s = GE.transform_integrals(X, h, g)
    row = {'mode': 'case', 'N': N, 's': s_ang, 'eps': eps, 'R_ang': R / ANG,
           'dmin_ang': float(D[np.triu_indices(N, 1)].min() / ANG),
           'order_rule': {str(c): order_rule(N, c) for c in range(1, N // 2 + 1)},
           'leak': leakage(N, h_s, g_s), 't': round(time.time() - t0, 2)}
    print(json.dumps(row), flush=True)
    with open(os.path.join(HERE, 'data', 'ellipse_ring_%s.jsonl' % out), 'a') as f:
        f.write(json.dumps(row) + '\n')


if __name__ == '__main__':
    if sys.argv[1] == 'case':
        args = sys.argv[2:]
        out = args[args.index('--out') + 1]
        args = args[:args.index('--out')]
        N, s_ang = int(args[0]), float(args[1])
        for e in args[2:]:
            case(N, s_ang, float(e), out)
