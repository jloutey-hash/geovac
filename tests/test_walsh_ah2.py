"""m-seam rung 5: Walsh's rule for AH2 (H2O bent, BeH2 linear) from one-electron
orbital bending slopes plus Pauli filling.

Model (parameter-free): one electron in a superposition of screened atomic
potentials with the Fermi-Amaldi (1934) self-interaction factor f = (N-1)/N,
N = electrons in the molecule:

    V = -Z/r_A + f V_H[rho_A](r_A) + sum_H ( -1/r_H + f V_H[rho_1s](r_H) )

rho_A = the atom's own spherical Fermi-Amaldi SCF density (factor (N_A-1)/N_A),
rho_1s = exact hydrogen 1s.  Only nuclear charges, electron counts and bond
lengths enter.  Engine = rung 4 (tests/test_renner_teller_rung4.py): FE-DVR
radial x one-centre partial waves, exact second-order bending slopes at
linear geometry, finite-bend referee in C2v blocks.

Driver: debug/walsh_ah2.py (slopes / curves; data in debug/data/).
"""
import functools
import importlib.util
import os

import numpy as np
import pytest
from numpy.polynomial import legendre as npleg
from scipy.special import spherical_in, spherical_kn

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_s = importlib.util.spec_from_file_location('rung4', os.path.join(ROOT, 'tests', 'test_renner_teller_rung4.py'))
g = importlib.util.module_from_spec(_s)
_s.loader.exec_module(g)

KAPPA = 2.0   # hydrogen 1s density decay (exact)

MOLECULES = {
    # Z, atomic occupations {l: [occ of 1st, 2nd, ... level]}, N_atom, N_mol, R_AH (bohr)
    'H2O': dict(Z=8, occ={0: [2, 2], 1: [4]}, R=1.809),
    'BeH2': dict(Z=4, occ={0: [2, 2]}, R=2.540),
}


# ---------------------------------------------------------------- radial grid with weights
def fedvr_w(bounds, order):
    x, w = g.lobatto(order)
    D = g.lagrange_deriv(x)
    nodes, weights, blocks = [], [], []
    for e in range(len(bounds) - 1):
        a, b = bounds[e], bounds[e + 1]
        h = (b - a) / 2
        rr = a + h * (x + 1)
        we = h * w
        Te = 0.5 * (D.T * we) @ D / h ** 2
        start = len(nodes) - 1 if e > 0 else 0
        if e == 0:
            nodes.extend(rr); weights.extend(we)
        else:
            weights[-1] += we[0]
            nodes.extend(rr[1:]); weights.extend(we[1:])
        blocks.append((start, Te))
    N = len(nodes)
    T = np.zeros((N, N))
    for start, Te in blocks:
        T[start:start + order, start:start + order] += Te
    W = np.array(weights)
    T = T / np.sqrt(np.outer(W, W))
    k = slice(1, N - 1)
    return np.array(nodes)[k], T[k, k], W[k], bounds[-1]


def bounds_for(R):
    return sorted(set([0, 0.03, 0.07, 0.13, 0.22, 0.35, 0.55, 0.8, 1.1, 1.45, R, R + 0.35, R + 0.9,
                       R + 1.8, R + 3, R + 4.8, R + 7.5, R + 11, R + 16, R + 23, R + 33, R + 46, R + 62]))


# ---------------------------------------------------------------- atoms
def hartree(P, r, T, W, rmax):
    """Spherical Hartree potential of radial density P (int P dr = N)."""
    N = float(np.sum(P * W))
    a = np.linalg.solve(2 * T, np.sqrt(W) * P / r)
    U = a / np.sqrt(W) + N * r / rmax
    return U / r


def atom_density(Z, occ, r, T, W, rmax, tol=1e-11):
    N = sum(sum(v) for v in occ.values())
    fac = (N - 1) / N
    VH = np.zeros_like(r)
    for it in range(500):
        P = np.zeros_like(r)
        eps = {}
        for l, occs in occ.items():
            e, V = np.linalg.eigh(T + np.diag(l * (l + 1) / (2 * r ** 2) - Z / r + fac * VH))
            eps[l] = e[:len(occs)]
            for i, n in enumerate(occs):
                P += n * V[:, i] ** 2 / W
        new = hartree(P, r, T, W, rmax)
        if np.max(np.abs(new - VH)) < tol:
            return P, eps
        VH = 0.5 * VH + 0.5 * new
    raise RuntimeError('atom SCF did not converge')


# ---------------------------------------------------------------- molecular potential
def proton_multipoles(r, R, fac, kmax):
    """g_k(r): v(d) = -(1-f)/d - f e^{-2d}/d - f e^{-2d} = sum_k g_k(r) P_k(cos)."""
    rl, rg = np.minimum(r, R), np.maximum(r, R)
    G = []
    for k in range(kmax + 1):
        ia, ka = spherical_in(k, KAPPA * rl), spherical_kn(k, KAPPA * rg)
        dia, dka = spherical_in(k, KAPPA * rl, derivative=True), spherical_kn(k, KAPPA * rg, derivative=True)
        Y = (2 * KAPPA / np.pi) * (2 * k + 1) * ia * ka
        dY = (2 / np.pi) * (2 * k + 1) * (ia * ka + KAPPA * (rl * dia * ka + rg * ia * dka))
        point = rl ** k / rg ** (k + 1)
        G.append(-(1 - fac) * point - fac * Y + fac * dY)
    return G


def molecule(name, order=10):
    spec = MOLECULES[name]
    R = spec['R']
    r, T, W, rmax = fedvr_w(bounds_for(R), order)
    P, eps = atom_density(spec['Z'], spec['occ'], r, T, W, rmax)
    N_mol = sum(sum(v) for v in spec['occ'].values()) + 2
    fac = (N_mol - 1) / N_mol
    Vc = -spec['Z'] / r + fac * hartree(P, r, T, W, rmax)
    return dict(name=name, R=R, r=r, T=T, W=W, rmax=rmax, Vc=Vc, fac=fac, atom_eps=eps)


# ---------------------------------------------------------------- operators
@functools.lru_cache(maxsize=None)
def _ang(bra, ket, kind, k):
    if kind == 'v0':
        f = np.outer(2 * g.Pk(k), np.ones(g.NPHI))
    elif kind == 'v1':
        f = g.v1_fun(k)
    else:
        f = g.v2_fun(k)
    return g.ang_matrix(list(bra), list(ket), f)


def op(mol, G, bra, ket, kset, kind):
    r = mol['r']
    Nr = len(r)
    M = np.zeros((len(bra), Nr, len(ket), Nr))
    for k in kset:
        A = _ang(tuple(bra), tuple(ket), kind, k)
        for a in range(len(bra)):
            for b in range(len(ket)):
                if abs(A[a, b]) > 1e-14:
                    M[a, :, b, :] += np.diag(A[a, b] * G[k])
    return M.reshape(len(bra) * Nr, len(ket) * Nr)


def h0(mol, G, labs, kmax):
    r, T = mol['r'], mol['T']
    Nr, L = len(r), len(labs)
    H = op(mol, G, labs, labs, range(0, kmax + 1, 2), 'v0').reshape(L, Nr, L, Nr)
    for a, lab in enumerate(labs):
        H[a, :, a, :] += T + np.diag(lab[0] * (lab[0] + 1) / (2 * r ** 2) + mol['Vc'])
    return H.reshape(L * Nr, L * Nr)


TARGETS = {(0, 'c'): [(1, 'c')], (1, 'c'): [(0, 'c'), (2, 'c')], (1, 's'): [(2, 's')]}

# Valence orbitals at linear geometry: (label, m, kind, z-parity, index in block)
ORBITALS = [('2sg', 0, 'c', +1, 1), ('1su', 0, 'c', -1, 0),
            ('1pi_in', 1, 'c', +1, 0), ('1pi_out', 1, 's', +1, 0)]


def slopes(mol, lmax):
    """Energy and exact s^2 coefficient of each valence orbital at linear geometry."""
    kmax = 2 * lmax
    G = proton_multipoles(mol['r'], mol['R'], mol['fac'], kmax)
    out = {}
    for name, m, kind, zp, idx in ORBITALS:
        labs = g.harmonics(m, kind, lmax, zp)
        ev, V = np.linalg.eigh(h0(mol, G, labs, kmax))
        E, phi = ev[idx], V[:, idx]
        direct = phi @ op(mol, G, labs, labs, range(0, kmax + 1, 2), 'v2') @ phi
        relay = 0.0
        for mt, kt in TARGETS[(m, kind)]:
            lt = g.harmonics(mt, kt, lmax, zp)
            rhs = op(mol, G, lt, labs, range(1, kmax + 1, 2), 'v1') @ phi
            Ht = h0(mol, G, lt, kmax)
            relay += rhs @ np.linalg.solve(Ht - E * np.eye(len(Ht)), -rhs)
        out[name] = {'E': float(ev[idx]), 'a2': float(direct + relay),
                     'direct': float(direct), 'relay': float(relay)}
    return out


# ---------------------------------------------------------------- finite bends (C2v blocks)
def bent_levels(mol, lmax, sv, nlev=3):
    """Lowest levels of the four C2v blocks (y-parity c/s x z-parity +/-) at bend s."""
    r, T = mol['r'], mol['T']
    Nr = len(r)
    kmax = 2 * lmax
    G = proton_multipoles(r, mol['R'], mol['fac'], kmax)
    c = np.sqrt(1 - sv ** 2)
    cosg_p = c * g.U[:, None] + sv * g.SIN[:, None] * g.COSPHI[None, :]
    cosg_m = -c * g.U[:, None] + sv * g.SIN[:, None] * g.COSPHI[None, :]
    funcs = []
    for k in range(kmax + 1):
        ck = np.zeros(k + 1); ck[k] = 1
        funcs.append(npleg.legval(cosg_p, ck) + npleg.legval(cosg_m, ck))
    res = {}
    for kind in ('c', 's'):
        for zp in (+1, -1):
            labs = [(l, m, kind) for l in range(lmax + 1) for m in range(0 if kind == 'c' else 1, l + 1)
                    if (-1) ** (l + m) == zp]
            L = len(labs)
            H = np.zeros((L, Nr, L, Nr))
            for a, lab in enumerate(labs):
                H[a, :, a, :] += T + np.diag(lab[0] * (lab[0] + 1) / (2 * r ** 2) + mol['Vc'])
            for k in range(kmax + 1):
                A = g.ang_matrix(labs, labs, funcs[k])
                for a in range(L):
                    for b in range(L):
                        if abs(A[a, b]) > 1e-13:
                            H[a, :, b, :] += np.diag(A[a, b] * G[k])
            res[(kind, zp)] = np.linalg.eigvalsh(H.reshape(L * Nr, L * Nr))[:nlev]
    return res


def walsh_levels(mol, lmax, sv):
    """Valence orbitals by adiabatic C2v index: 2a1, 3a1 (c,+), 1b2 (c,-), 1b1 (s,+)."""
    b = bent_levels(mol, lmax, sv)
    return {'2a1': b[('c', 1)][1], '3a1': b[('c', 1)][2], '1b2': b[('c', -1)][0], '1b1': b[('s', 1)][0]}




# ================================================================ tests
def test_hartree_and_multipoles_are_exact():
    mol = molecule('H2O')
    r = mol['r']
    VH = hartree(4 * r ** 2 * np.exp(-2 * r), r, mol['T'], mol['W'], mol['rmax'])
    assert np.max(np.abs(VH - (1 / r - (1 + 1 / r) * np.exp(-2 * r)))) < 1e-10
    R, fac = mol['R'], mol['fac']
    G = proton_multipoles(r, R, fac, 60)
    for ri, u in ((0.5, 0.9), (1.2, -0.7), (5.0, 0.3)):
        i = int(np.argmin(np.abs(r - ri)))
        d = np.sqrt(r[i] ** 2 + R ** 2 - 2 * r[i] * R * u)
        direct = -(1 - fac) / d - fac * np.exp(-2 * d) / d - fac * np.exp(-2 * d)
        assert abs(sum(G[k][i] * npleg.legval(u, np.eye(61)[k]) for k in range(61)) - direct) < 1e-9


def test_screening_puts_2s_below_2p():
    """The parameter-free screened atom breaks the hydrogenic 2s = 2p tie."""
    e = molecule('H2O')['atom_eps']
    assert e[0][1] < e[1][0] - 0.3


def test_slopes_match_real_bends():
    """Exact second-order slopes == finite-bend Richardson, all four valence
    orbitals, both molecules (C2v index tracking)."""
    hs = [0.04, 0.02, 0.01]
    for name in MOLECULES:
        mol = molecule(name, order=8)
        sl = slopes(mol, 5)
        lv0 = walsh_levels(mol, 5, 0.0)
        lv = [walsh_levels(mol, 5, h) for h in hs]
        for orb, key in (('2sg', '2a1'), ('1su', '1b2'), ('1pi_in', '3a1'), ('1pi_out', '1b1')):
            q = [(lv[i][key] - lv0[key]) / hs[i] ** 2 for i in range(3)]
            assert abs((4 * q[2] - q[1]) / 3 - sl[orb]['a2']) < 2e-5


def test_walsh_diagram_signs():
    """Walsh 1953, at l_max = 15: 3a1 (in-plane pi) falls, 1b2 rises, 1b1 is
    nearly flat; water's occupied-orbital sum makes it bend."""
    for name in MOLECULES:
        a = {k: v['a2'] for k, v in slopes(molecule(name), 15).items()}
        assert a['1pi_in'] < -0.05
        assert a['1su'] > 0.05
        assert abs(a['1pi_out']) < 0.01
        if name == 'H2O':
            assert 2 * (a['2sg'] + a['1su'] + a['1pi_in'] + a['1pi_out']) < -0.1
