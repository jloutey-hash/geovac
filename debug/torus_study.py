"""Hydrogen torus study (2026-11-01, PI-directed: "Ok. Do that doughnut."): do the two winding labels of a torus
survive in a molecule?  The toroidal one (around the hole) is an exact symmetry; the poloidal one (around the tube)
is not, because the inner side of the tube faces the hole and the outer side does not.

Molecule: N small rings of M hydrogens, centred on a circle of radius R in the xy plane, each small ring in the
plane through the z axis (a "ring of rings").  Atom (j, a): toroidal angle phi_j = 2 pi j / N, poloidal angle
theta_a = 2 pi a / M (a = 0 on the outer equator), position ((R + r cos theta) cos phi, (R + r cos theta) sin phi,
r sin theta).  R is set by the spacing s of neighbouring ring centres, R = s / (2 sin(pi/N)); r by the in-ring
neighbour distance d, r = d / (2 sin(pi/M)).  So with s and d fixed, more rings means a thinner torus, r/R -> 0.

Basis: one Slater 1s per H (zeta 1.24, fitted by 3 Gaussians; a minimal basis, enough for a symmetry study).
Orthonormal site orbitals by symmetric (Lowdin) orthogonalization, which commutes with every symmetry of the
molecule.  Labellings of that one space:
  site : the site orbitals, index (j, a).
  tor  : real discrete Fourier transform over j for each a (toroidal winding k, exact).
  full : real DFT over j and over a (k and the poloidal winding m, approximate).
Measured:
  leakage, in the complex (k, m) basis: the share of the Frobenius norm of h and of (pq|rs) carried by entries
    that break k conservation (must vanish) and that break m conservation, split by |delta m| mod M;
  per real labelling: Z-type Z2 count, exactly-zero integral fraction, JW Pauli counts and 1-norm;
  the RHF energy change when the m-breaking entries are dropped.

    python debug/torus_study.py case <N> <M> <s_ang> <d_ang> [--light]
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
OUT = os.path.join(HERE, 'data', 'torus_study.jsonl')
ANG = hs.ANG
ZETA_H, NG_MIN = 1.24, 3


def z2_count(h, g, tol=1e-10):
    """hring_study.z2_count with Python-int bit masks: the original builds masks from numpy int64 indices and
    overflows past 64 qubits.  hring_study.py is left untouched (its hash is in the HR record)."""
    M = h.shape[0]
    rows = set()
    for p in range(M):
        for q in range(M):
            if p != q and abs(h[p, q]) > tol:
                for s in (0, 1):
                    rows.add((1 << (2 * p + s)) ^ (1 << (2 * q + s)))
    for p, q, r, s in np.argwhere(np.abs(g) > tol).tolist():
        for a in (0, 1):
            for b in (0, 1):
                v = (1 << (2 * p + a)) ^ (1 << (2 * q + a)) ^ (1 << (2 * r + b)) ^ (1 << (2 * s + b))
                if v:
                    rows.add(v)
    basis = []
    for v in rows:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
            basis.sort(reverse=True)
    return 2 * M - len(basis)


_z2_count_orig = hs.z2_count
hs.z2_count = z2_count          # label_measures looks the name up at call time


def write(row):
    print(json.dumps(row), flush=True)
    with open(OUT, 'a') as f:
        f.write(json.dumps(row) + '\n')


def torus(N, M, s_ang, d_ang):
    R = s_ang * ANG / (2 * np.sin(np.pi / N))
    r = d_ang * ANG / (2 * np.sin(np.pi / M))
    nuc = []
    for j in range(N):
        ph = 2 * np.pi * j / N
        for a in range(M):
            th = 2 * np.pi * a / M
            rho = R + r * np.cos(th)
            nuc.append((np.array([rho * np.cos(ph), rho * np.sin(ph), r * np.sin(th)]), 1.0))
    return nuc, R, r


def minimal_basis(nuc):
    fns, blocks = [], []
    for C, _ in nuc:
        w.sto_shell(C, 1, 0, ZETA_H, NG_MIN, fns, blocks)
    return fns, w.spherical_map(fns, blocks)


def complex_dft(N, M):
    """columns indexed (k, m), rows (j, a); site index j * M + a"""
    j, a = np.arange(N), np.arange(M)
    Fk = np.exp(2j * np.pi * np.outer(j, np.arange(N)) / N) / np.sqrt(N)
    Fm = np.exp(2j * np.pi * np.outer(a, np.arange(M)) / M) / np.sqrt(M)
    return np.kron(Fk, Fm)


def labels_km(N, M):
    k = np.repeat(np.arange(N), M)
    m = np.tile(np.arange(M), N)
    return k, m


def transform_c(U, h, g):
    hc = U.conj().T @ h @ U
    gc = np.einsum('pqrs,pa,qb,rc,sd->abcd', g, U.conj(), U, U.conj(), U, optimize=True)
    return hc, gc


def leakage(hc, gc, N, M):
    k, m = labels_km(N, M)
    dk1 = (k[None, :] - k[:, None]) % N
    dm1 = (m[None, :] - m[:, None]) % M
    dk2 = (-k[:, None, None, None] + k[None, :, None, None] - k[None, None, :, None] + k[None, None, None, :]) % N
    dm2 = (-m[:, None, None, None] + m[None, :, None, None] - m[None, None, :, None] + m[None, None, None, :]) % M
    out = {}
    for name, X, dk, dm in (('h', hc, dk1, dm1), ('g', gc, dk2, dm2)):
        a2 = np.abs(X) ** 2
        tot = a2.sum()
        out[name + '_k_break'] = float(np.sqrt(a2[dk != 0].sum() / tot))
        cls = np.minimum(dm, (M - dm) % M)          # |delta m| mod M
        for c in range(1, M // 2 + 1):
            out[name + '_m%d' % c] = float(np.sqrt(a2[(dk == 0) & (cls == c)].sum() / tot))
        out[name + '_m_any'] = float(np.sqrt(a2[(dk == 0) & (cls != 0)].sum() / tot))
    return out


def drop_m_breaking(hc, gc, N, M):
    k, m = labels_km(N, M)
    dm1 = (m[None, :] - m[:, None]) % M
    dm2 = (-m[:, None, None, None] + m[None, :, None, None] - m[None, None, :, None] + m[None, None, None, :]) % M
    return np.where(dm1 == 0, hc, 0), np.where(dm2 == 0, gc, 0)


def real_labellings(N, M):
    Fn, Fm = hs.real_dft(N), hs.real_dft(M)
    return {'site': np.eye(N * M), 'tor': np.kron(Fn, np.eye(M)), 'full': np.kron(Fn, Fm)}


def mode_case(N, M, s_ang, d_ang, light=False):
    t0 = time.time()
    nuc, R, r = torus(N, M, s_ang, d_ang)
    enuc = w.vnn(nuc)
    P = np.array([C for C, _ in nuc])
    D = np.linalg.norm(P[:, None] - P[None], axis=-1)
    dmin = float(D[np.triu_indices(len(nuc), 1)].min() / ANG)
    fns, T = minimal_basis(nuc)
    S, h, g = w.ao_integrals(fns, nuc, T)
    t_int = time.time() - t0
    X = GE.lowdin_orbitals(S)
    h_s, g_s = GE.transform_integrals(X, h, g)
    del g
    U = complex_dft(N, M)
    hc, gc = transform_c(U, h_s, g_s)
    row = {'mode': 'case', 'N': N, 'M': M, 's': s_ang, 'd': d_ang, 'R_ang': R / ANG, 'r_ang': r / ANG,
           'r_over_R': r / R, 'dmin_ang': dmin, 'nbf': len(nuc), 'cond_S': float(np.linalg.cond(S)),
           'leak': leakage(hc, gc, N, M)}
    ht, gt = drop_m_breaking(hc, gc, N, M)
    del gc
    h_t = (U @ ht @ U.conj().T)
    g_t = np.einsum('abcd,pa,qb,rc,sd->pqrs', gt, U, U.conj(), U, U.conj(), optimize=True)
    del gt
    imag = float(max(np.abs(h_t.imag).max(), np.abs(g_t.imag).max()))
    dh, dg = h_s - h_t.real, g_s - g_t.real
    del h_t, g_t
    # Dropping the m-breaking part dH moves every energy by at most ||dH|| <= its Pauli 1-norm (a bound, no SCF).
    Vf = real_labellings(N, M)['full']
    dh_f, dg_f = GE.transform_integrals(Vf, dh, dg)
    st = w.jw_pauli_stats(dh_f, dg_f, 0.0)
    row['m_breaking_part'] = {'n_pauli': int(st['n_pauli']), 'one_norm': float(st['one_norm']),
                              'imag_after_drop': imag}
    del dh, dg, dh_f, dg_f
    if not light:
        row['labels'] = {}
        for name, V in real_labellings(N, M).items():
            hh, gg = GE.transform_integrals(V, h_s, g_s)
            row['labels'][name] = hs.label_measures(hh, gg, enuc)
        row['m_breaking_part']['share_of_full_one_norm'] = (row['m_breaking_part']['one_norm']
                                                           / row['labels']['full']['one_norm'])
    row['t_int'] = round(t_int, 1)
    row['t'] = round(time.time() - t0, 1)
    write(row)


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'case':
        mode_case(int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]),
                  light='--light' in sys.argv)
