"""Water with exponential (Slater) vs Gaussian (cc-pVDZ) bases: RHF + CCSD(T).

PySCF does not build on this Windows host, so the correlated solver is our own
and is validated where its answer is known:
  * CCSD == FCI exactly for two electrons;
  * CCSD and CCSD(T) bracket FCI tightly on minimal-basis water.

Integrals: geovac.noci_engine (McMurchie-Davidson over contracted Cartesian
Gaussians; Slater orbitals enter as Gaussian fits).  cc-pVDZ from the Basis Set
Exchange (debug/data/basis/cc-pvdz_H_O.json), d shells transformed to the five
real spherical components.  All-electron throughout.
Driver: debug/water_ccsd.py.  Memo: debug/water_barrier_memo.md.
"""
import json
import math
import os

import numpy as np
import pytest

from geovac import noci_engine as GE

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

R_OH, ANGLE = 1.809, 104.5


def water_nuclei(R=R_OH, angle=ANGLE):
    th = math.radians(angle / 2)
    O = np.zeros(3)
    H1 = np.array([R * math.sin(th), 0.0, R * math.cos(th)])
    H2 = np.array([-R * math.sin(th), 0.0, R * math.cos(th)])
    return [(O, 8.0), (H1, 1.0), (H2, 1.0)]


def vnn(nuclei):
    e = 0.0
    for i in range(len(nuclei)):
        for j in range(i + 1, len(nuclei)):
            e += nuclei[i][1] * nuclei[j][1] / np.linalg.norm(nuclei[i][0] - nuclei[j][0])
    return e


# ---------------------------------------------------------------- Gaussian basis
CART = {0: [(0, 0, 0)], 1: [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
        2: [(2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)],
        3: [(3, 0, 0), (0, 3, 0), (0, 0, 3), (2, 1, 0), (2, 0, 1), (1, 2, 0),
            (0, 2, 1), (1, 0, 2), (0, 1, 2), (1, 1, 1)]}
# real spherical d in terms of the normalised Cartesians above (unnormalised rows)
SPH_D = np.array([[0, 0, 0, 1, 0, 0],          # xy
                  [0, 0, 0, 0, 0, 1],          # yz
                  [-1, -1, 2, 0, 0, 0],        # z^2
                  [0, 0, 0, 0, 1, 0],          # xz
                  [1, -1, 0, 0, 0, 0]], float)  # x^2 - y^2
# real solid harmonics l = 3 as polynomials {(l,m,n): coeff}; on the normalized
# Cartesians each monomial carries sqrt((2l-1)!!(2m-1)!!(2n-1)!!)
_F_POLY = [{(0, 0, 3): 2, (2, 0, 1): -3, (0, 2, 1): -3},
           {(1, 0, 2): 4, (3, 0, 0): -1, (1, 2, 0): -1},
           {(0, 1, 2): 4, (2, 1, 0): -1, (0, 3, 0): -1},
           {(2, 0, 1): 1, (0, 2, 1): -1},
           {(1, 1, 1): 1},
           {(3, 0, 0): 1, (1, 2, 0): -3},
           {(2, 1, 0): 3, (0, 3, 0): -1}]


def _dfac(k):
    return 1 if k <= 0 else k * _dfac(k - 2)


SPH = {2: SPH_D,
       3: np.array([[poly.get(c, 0) * math.sqrt(_dfac(2 * c[0] - 1) * _dfac(2 * c[1] - 1) * _dfac(2 * c[2] - 1))
                     for c in CART[3]] for poly in _F_POLY], float)}


def gaussian_basis(path, nuclei):
    """Contracted Cartesian functions + the Cartesian->spherical map T (ncart x nsph)."""
    data = json.load(open(path))
    fns, blocks = [], []
    for C, Z in nuclei:
        el = data['elements'][str(int(round(Z)))]
        for sh in el['electron_shells']:
            am = sh['angular_momentum']    # [l] (general contraction) or [0, 1] (SP shell)
            al = np.array([float(x) for x in sh['exponents']])
            for k, col in enumerate(sh['coefficients']):
                l = am[k] if len(am) > 1 else am[0]
                d = np.array([float(x) for x in col])
                start = len(fns)
                for lmn in CART[l]:
                    fns.append(GE.BasisFn(C, lmn, al, d))
                blocks.append((l, start))
    ncart = len(fns)
    cols = []
    for l, start in blocks:
        if l < 2:
            for k in range(len(CART[l])):
                v = np.zeros(ncart); v[start + k] = 1.0; cols.append(v)
        else:
            for row in SPH[l]:
                v = np.zeros(ncart); v[start:start + len(row)] = row; cols.append(v)
    return fns, np.array(cols).T


def ao_integrals(fns, nuclei, T=None):
    S, h, g = GE.integral_set_md(fns, nuclei)
    if T is not None:
        S = T.T @ S @ T
        h = T.T @ h @ T
        g = np.einsum('pqrs,pa,qb,rc,sd->abcd', g, T, T, T, T, optimize=True)
        n = 1 / np.sqrt(np.diag(S))
        S, h = S * np.outer(n, n), h * np.outer(n, n)
        g = g * np.einsum('a,b,c,d->abcd', n, n, n, n)
    return S, h, g


def slater_minimal(nuclei, n_gauss=12):
    sh = {}
    for kind, (l, n_r) in (('1s', (0, 1)), ('2s', (0, 2)), ('2p', (1, 2))):
        a, d, _ = GE.fit_sto_shape(l, n_r, n_gauss=n_gauss)
        sh[kind] = (a, d)
    (O, _), (H1, _), (H2, _) = nuclei
    fns = [GE.sto_shape_basis(O, '1s', 7.66, sh, (0, 0, 0)),
           GE.sto_shape_basis(O, '2s', 2.25, sh, (0, 0, 0))]
    fns += [GE.sto_shape_basis(O, '2p', 2.23, sh, lmn) for lmn in CART[1]]
    fns += [GE.sto_shape_basis(H, '1s', 1.0, sh, (0, 0, 0)) for H in (H1, H2)]
    return fns


# ---------------------------------------------------------------- RHF
def rhf(S, h, g, nelec, enuc, tol=1e-11, maxit=300):
    X = GE.lowdin_orbitals(S)
    nocc = nelec // 2
    Fp = X.T @ h @ X
    e, Cp = np.linalg.eigh(Fp)
    C = X @ Cp
    D = 2 * C[:, :nocc] @ C[:, :nocc].T
    hist_F, hist_e = [], []
    Eold = 0.0
    for it in range(maxit):
        J = np.einsum('pqrs,rs->pq', g, D)
        K = np.einsum('prqs,rs->pq', g, D)
        F = h + J - 0.5 * K
        E = 0.5 * np.sum(D * (h + F)) + enuc
        err = X.T @ (F @ D @ S - S @ D @ F) @ X
        hist_F.append(F); hist_e.append(err)
        if len(hist_F) > 8:
            hist_F.pop(0); hist_e.pop(0)
        if np.max(np.abs(err)) < tol and abs(E - Eold) < tol:
            break
        Eold = E
        n = len(hist_F)
        B = -np.ones((n + 1, n + 1)); B[n, n] = 0
        for i in range(n):
            for j in range(n):
                B[i, j] = np.sum(hist_e[i] * hist_e[j])
        rhs = np.zeros(n + 1); rhs[n] = -1
        c = np.linalg.lstsq(B, rhs, rcond=None)[0][:n]   # B singular if err == 0
        Fd = sum(ci * Fi for ci, Fi in zip(c, hist_F))
        e, Cp = np.linalg.eigh(X.T @ Fd @ X)
        C = X @ Cp
        D = 2 * C[:, :nocc] @ C[:, :nocc].T
    else:
        raise RuntimeError('RHF did not converge')
    J = np.einsum('pqrs,rs->pq', g, D)
    K = np.einsum('prqs,rs->pq', g, D)
    F = h + J - 0.5 * K
    eps, Cp = np.linalg.eigh(X.T @ F @ X)
    return E, eps, X @ Cp


# ---------------------------------------------------------------- CCSD(T), spin orbitals
def spin_orbital(h, g, C, eps):
    hm = C.T @ h @ C
    gm = np.einsum('pqrs,pa,qb,rc,sd->abcd', g, C, C, C, C, optimize=True)
    n = len(eps)
    spin = np.arange(2 * n) % 2
    sp = np.arange(2 * n) // 2
    same = (spin[:, None] == spin[None, :]).astype(float)
    # <pq|rs> = (pr|qs) delta(sp,sr) delta(sq,ss)
    phys = gm[np.ix_(sp, sp, sp, sp)].transpose(0, 2, 1, 3) * same[:, None, :, None] * same[None, :, None, :]
    anti = phys - phys.transpose(0, 1, 3, 2)
    fs = np.diag(np.repeat(eps, 2))
    return fs, anti


def _es(*a):
    return np.einsum(*a, optimize=True)


def ccsd(fs, g, nocc, tol=1e-10, maxit=200, triples=True):
    """nocc = number of occupied SPIN orbitals.  Returns (E_ccsd_corr, E_T)."""
    o, v = slice(0, nocc), slice(nocc, fs.shape[0])
    eo, ev = np.diag(fs)[o], np.diag(fs)[v]
    Dia = eo[:, None] - ev[None, :]
    Dijab = eo[:, None, None, None] + eo[None, :, None, None] - ev[None, None, :, None] - ev[None, None, None, :]
    t1 = np.zeros_like(Dia)
    t2 = g[o, o, v, v] / Dijab
    Eold = 0.0
    oovv, ovvo, ooov, vovv = g[o, o, v, v], g[o, v, v, o], g[o, o, o, v], g[v, o, v, v]
    for it in range(maxit):
        tau_t = t2 + 0.5 * (_es('ia,jb->ijab', t1, t1) - _es('ib,ja->ijab', t1, t1))
        tau = t2 + _es('ia,jb->ijab', t1, t1) - _es('ib,ja->ijab', t1, t1)
        Fae = (-0.5 * _es('me,ma->ae', fs[o, v], t1)
               + _es('mf,mafe->ae', t1, g[o, v, v, v])
               - 0.5 * _es('mnaf,mnef->ae', tau_t, oovv))
        Fmi = (0.5 * _es('ie,me->mi', t1, fs[o, v])
               + _es('ne,mnie->mi', t1, ooov)
               + 0.5 * _es('inef,mnef->mi', tau_t, oovv))
        Fme = fs[o, v] + _es('nf,mnef->me', t1, oovv)
        Wmnij = (g[o, o, o, o] + _es('je,mnie->mnij', t1, ooov) - _es('ie,mnje->mnij', t1, ooov)
                 + 0.25 * _es('ijef,mnef->mnij', tau, oovv))
        Wabef = (g[v, v, v, v] - _es('mb,amef->abef', t1, g[v, o, v, v])
                 + _es('ma,bmef->abef', t1, g[v, o, v, v])
                 + 0.25 * _es('mnab,mnef->abef', tau, oovv))
        Wmbej = (g[o, v, v, o] + _es('jf,mbef->mbej', t1, g[o, v, v, v])
                 - _es('nb,mnej->mbej', t1, g[o, o, v, o])
                 - _es('jnfb,mnef->mbej', 0.5 * t2 + _es('jf,nb->jnfb', t1, t1), oovv))
        # T1
        r1 = (fs[o, v] + _es('ie,ae->ia', t1, Fae) - _es('ma,mi->ia', t1, Fmi)
              + _es('imae,me->ia', t2, Fme) - _es('nf,naif->ia', t1, g[o, v, o, v])
              - 0.5 * _es('imef,maef->ia', t2, g[o, v, v, v])
              - 0.5 * _es('mnae,nmei->ia', t2, g[o, o, v, o]))
        # T2
        Ft = Fae - 0.5 * _es('mb,me->be', t1, Fme)
        P1 = _es('ijae,be->ijab', t2, Ft)
        Fo = Fmi + 0.5 * _es('je,me->mj', t1, Fme)
        P2 = _es('imab,mj->ijab', t2, Fo)
        X = _es('imae,mbej->ijab', t2, Wmbej) - _es('ie,ma,mbej->ijab', t1, t1, g[o, v, v, o])
        Y = _es('ie,abej->ijab', t1, g[v, v, v, o])
        Zt = _es('ma,mbij->ijab', t1, g[o, v, o, o])
        r2 = (oovv.copy()
              + P1 - P1.transpose(0, 1, 3, 2)
              - P2 + P2.transpose(1, 0, 2, 3)
              + 0.5 * _es('mnab,mnij->ijab', tau, Wmnij)
              + 0.5 * _es('ijef,abef->ijab', tau, Wabef)
              + X - X.transpose(1, 0, 2, 3) - X.transpose(0, 1, 3, 2) + X.transpose(1, 0, 3, 2)
              + Y - Y.transpose(1, 0, 2, 3)
              - Zt + Zt.transpose(0, 1, 3, 2))
        t1 = r1 / Dia
        t2 = r2 / Dijab
        E = (_es('ia,ia->', fs[o, v], t1) + 0.25 * _es('ijab,ijab->', oovv, t2)
             + 0.5 * _es('ijab,ia,jb->', oovv, t1, t1))
        if abs(E - Eold) < tol and it > 2:
            break
        Eold = E
    else:
        raise RuntimeError('CCSD did not converge')
    if not triples:
        return E, 0.0
    # (T)
    nv = fs.shape[0] - nocc
    vvvo, ovoo = g[v, v, v, o], g[o, v, o, o]
    ET = 0.0
    for i in range(nocc):
        for j in range(i + 1, nocc):
            for k in range(j + 1, nocc):
                def conn(i, j, k):
                    # sum_e t_jk^ae <ei||bc>  - sum_m t_im^bc <ma||jk>, as [a,b,c]
                    A = _es('ae,ebc->abc', t2[j, k], g[v, o, v, v][:, i])
                    B = _es('mbc,ma->abc', t2[i], g[o, v, o, o][:, :, j, k])
                    return A - B

                def disc(i, j, k):
                    return _es('a,bc->abc', t1[i], oovv[j, k])

                def perm_abc(Xa):
                    return Xa - Xa.transpose(1, 0, 2) - Xa.transpose(2, 1, 0)

                tc = perm_abc(conn(i, j, k) - conn(j, i, k) - conn(k, j, i))
                td = perm_abc(disc(i, j, k) - disc(j, i, k) - disc(k, j, i))
                D3 = (eo[i] + eo[j] + eo[k] - ev[:, None, None] - ev[None, :, None] - ev[None, None, :])
                ET += np.sum(tc * (tc + td) / D3)
    ET /= 6.0    # i<j<k restriction (x 6 over the 1/36 of the full sum); a,b,c unrestricted
    return E, ET


def fci(S, h, g, nelec, enuc):
    X = GE.lowdin_orbitals(S)
    ht, gt = GE.transform_integrals(X, h, g)
    return GE.fci_ground(ht, gt, nelec) + enuc


def run(S, h, g, nelec, enuc, triples=True, frozen_core=0):
    """frozen_core: number of lowest canonical RHF orbitals left uncorrelated
    (they stay in the Fock operator, so the active block needs no extra terms)."""
    Ehf, eps, C = rhf(S, h, g, nelec, enuc)
    fc = frozen_core
    fs, so = spin_orbital(h, g, C[:, fc:], eps[fc:])
    Ec, Et = ccsd(fs, so, nelec - 2 * fc, triples=triples)
    return {'RHF': Ehf, 'CCSD': Ehf + Ec, 'CCSD(T)': Ehf + Ec + Et}


# ---------------------------------------------------------------- Slater DZP (3s2p1d / 2s1p)
_SHAPES = {}


def sto_shape(l, n_r, n_gauss):
    key = (l, n_r, n_gauss)
    if key not in _SHAPES:
        a, d, _ = GE.fit_sto_shape(l, n_r, n_gauss=n_gauss)
        _SHAPES[key] = (a, d)
    return _SHAPES[key]


def sto_shell(center, n, l, zeta, n_gauss, fns, blocks):
    a, d = sto_shape(l, n, n_gauss)
    start = len(fns)
    for lmn in CART[l]:
        fns.append(GE.BasisFn(center, lmn, a * zeta ** 2, d))
    blocks.append((l, start))


def spherical_map(fns, blocks):
    ncart, cols = len(fns), []
    for l, start in blocks:
        if l < 2:
            for k in range(len(CART[l])):
                v = np.zeros(ncart); v[start + k] = 1.0; cols.append(v)
        else:
            for row in SPH[l]:
                v = np.zeros(ncart); v[start:start + len(row)] = row; cols.append(v)
    return np.array(cols).T


def slater_dzp(nuclei, ex, n_gauss=12, polarization=True):
    """ex: dict O1s, O2s, O2s', O2p, O2p', O3d, H1s, H1s', H2p."""
    fns, blocks = [], []
    (O, _), (H1, _), (H2, _) = nuclei
    sto_shell(O, 1, 0, ex['O1s'], n_gauss, fns, blocks)
    sto_shell(O, 2, 0, ex['O2s'], n_gauss, fns, blocks)
    sto_shell(O, 2, 0, ex["O2s'"], n_gauss, fns, blocks)
    sto_shell(O, 2, 1, ex['O2p'], n_gauss, fns, blocks)
    sto_shell(O, 2, 1, ex["O2p'"], n_gauss, fns, blocks)
    if polarization:
        sto_shell(O, 3, 2, ex['O3d'], n_gauss, fns, blocks)
    for H in (H1, H2):
        sto_shell(H, 1, 0, ex['H1s'], n_gauss, fns, blocks)
        sto_shell(H, 1, 0, ex["H1s'"], n_gauss, fns, blocks)
        if polarization:
            sto_shell(H, 2, 1, ex['H2p'], n_gauss, fns, blocks)
    return fns, spherical_map(fns, blocks)


def rhf_fractional(S, h, g, occ, enuc=0.0, tol=1e-10, maxit=500):
    """RHF with fixed fractional occupations in aufbau order (spherically
    averaged open-shell atom, e.g. O: 2, 2, 4/3, 4/3, 4/3)."""
    X = GE.lowdin_orbitals(S)
    occ = np.asarray(occ, float)
    e, Cp = np.linalg.eigh(X.T @ h @ X)
    C = X @ Cp
    D = (C[:, :len(occ)] * occ) @ C[:, :len(occ)].T
    Eold = 0.0
    for it in range(maxit):
        F = h + np.einsum('pqrs,rs->pq', g, D) - 0.5 * np.einsum('prqs,rs->pq', g, D)
        E = 0.5 * np.sum(D * (h + F)) + enuc
        e, Cp = np.linalg.eigh(X.T @ F @ X)
        C = X @ Cp
        Dn = (C[:, :len(occ)] * occ) @ C[:, :len(occ)].T
        if abs(E - Eold) < tol and np.max(np.abs(Dn - D)) < 1e-8:
            return E
        D = 0.6 * Dn + 0.4 * D
        Eold = E
    raise RuntimeError('fractional RHF did not converge')


def rhf_atom_blocks(S, h, g, blocks, tol=1e-10, maxit=500, return_orbitals=False):
    """Spherical-atom RHF: Fock diagonalised within each symmetry block
    (index list, occupations), e.g. s: [2, 2], each p component: [4/3].
    return_orbitals: also return each block's occupied orbitals (block basis)."""
    D = np.zeros_like(S)
    F = h.copy()
    Eold = 0.0
    for it in range(maxit):
        Dn = np.zeros_like(S)
        orbs = []
        for idx, occ in blocks:
            idx = np.asarray(idx)
            Sb, Fb = S[np.ix_(idx, idx)], F[np.ix_(idx, idx)]
            Xb = GE.lowdin_orbitals(Sb)
            e, Cp = np.linalg.eigh(Xb.T @ Fb @ Xb)
            Cb = Xb @ Cp[:, :len(occ)]
            orbs.append(Cb)
            Dn[np.ix_(idx, idx)] += (Cb * np.asarray(occ, float)) @ Cb.T
        D = Dn if it == 0 else 0.6 * Dn + 0.4 * D
        F = h + np.einsum('pqrs,rs->pq', g, D) - 0.5 * np.einsum('prqs,rs->pq', g, D)
        E = 0.5 * np.sum(D * (h + F))
        if abs(E - Eold) < tol and it > 3:
            return (E, orbs) if return_orbitals else E
        Eold = E
    raise RuntimeError('atom RHF did not converge')


def atom_average_energy(S, h, g, blocks, open_blocks, q_open, tol=1e-10, maxit=500):
    """Average-of-configuration energy for one open shell (q_open electrons in
    the occupied orbital of each open block, all equally occupied).  Orbitals
    from the fractional-occupation SCF of rhf_atom_blocks; the energy removes
    its self-interaction: open-open spin-orbital pairs carry weight
    q(q-1)/(N(N-1)) instead of f^2 (f = q/N, N = 2 x number of open orbitals)."""
    E_frac = rhf_atom_blocks(S, h, g, blocks, tol, maxit)
    # rebuild the converged open orbitals
    D = np.zeros_like(S)
    F = h.copy()
    for it in range(maxit):
        Dn = np.zeros_like(S)
        opens = []
        for k, (idx, occ) in enumerate(blocks):
            idx = np.asarray(idx)
            Xb = GE.lowdin_orbitals(S[np.ix_(idx, idx)])
            e, Cp = np.linalg.eigh(Xb.T @ F[np.ix_(idx, idx)] @ Xb)
            Cb = Xb @ Cp[:, :len(occ)]
            Dn[np.ix_(idx, idx)] += (Cb * np.asarray(occ, float)) @ Cb.T
            if k in open_blocks:
                v = np.zeros(S.shape[0]); v[idx] = Cb[:, 0]; opens.append(v)
        Dnew = Dn if it == 0 else 0.6 * Dn + 0.4 * D
        if it > 3 and np.max(np.abs(Dnew - D)) < 1e-9:
            break
        D = Dnew
        F = h + np.einsum('pqrs,rs->pq', g, D) - 0.5 * np.einsum('prqs,rs->pq', g, D)
    Co = np.array(opens).T
    go = np.einsum('pqrs,pa,qb,rc,sd->abcd', g, Co, Co, Co, Co, optimize=True)
    no = Co.shape[1]
    J = np.array([[go[i, i, j, j] for j in range(no)] for i in range(no)])
    K = np.array([[go[i, j, i, j] for j in range(no)] for i in range(no)])
    N = 2 * no
    f = q_open / N
    fp = q_open * (q_open - 1) / (N * (N - 1))
    return E_frac - 0.5 * (f * f - fp) * (4 * J.sum() - 2 * K.sum())


# --- general-contraction Slater DZP (cc-pVDZ's design rule, Slater primitives) ---
# O primitives: s = two n=1 + two n=2 STOs, p = three n=2 STOs, exponents from the
# O-atom average-of-configuration energy.  Contracted to [3s2p] exactly as cc-pVDZ
# is: s = {atomic 1s, atomic 2s, most diffuse s primitive}, p = {atomic 2p, most
# diffuse p primitive}.  H [2s1p] and polarization as in slater_dzp.
O_S_N = (1, 1, 2, 2)


def o_atom_primitives(center, zs, zp, n_gauss):
    fns, blocks = [], []
    for n, z in zip(O_S_N, zs):
        sto_shell(center, n, 0, z, n_gauss, fns, blocks)
    for z in zp:
        sto_shell(center, 2, 1, z, n_gauss, fns, blocks)
    return fns, blocks


def o_atom_contraction(zs, zp, n_gauss=12):
    """Average-of-configuration O atom in the primitive set -> (E, s coeffs
    [1s, 2s], p coeffs [2p]) on the normalized primitives."""
    O = np.zeros(3)
    fns, _ = o_atom_primitives(O, zs, zp, n_gauss)
    S, h, g = GE.integral_set_md(fns, [(O, 8.0)])
    s = [i for i, f in enumerate(fns) if f.lmn == (0, 0, 0)]
    ps = [[i for i, f in enumerate(fns) if f.lmn == c] for c in CART[1]]
    blocks = [(s, [2, 2])] + [(p, [4 / 3]) for p in ps]
    E = atom_average_energy(S, h, g, blocks, [1, 2, 3], 4)
    _, orbs = rhf_atom_blocks(S, h, g, blocks, return_orbitals=True)
    return E, orbs[0][:, :2], orbs[1][:, 0]


def slater_dzp_contracted(nuclei, ex, n_gauss=12, polarization=True, ng_pol=None):
    """ex: Os (4 zetas), Op (3 zetas), cs (4x2), cp (3,), O3d, H1s, H1s', H2p.
    ng_pol: Gaussians per (cusp-free) polarization function; default n_gauss."""
    ng_pol = n_gauss if ng_pol is None else ng_pol
    (O, _), (H1, _), (H2, _) = nuclei
    fns, blocks = o_atom_primitives(O, ex['Os'], ex['Op'], n_gauss)
    nO = len(fns)
    if polarization:
        sto_shell(O, 3, 2, ex['O3d'], ng_pol, fns, blocks)
    for H in (H1, H2):
        sto_shell(H, 1, 0, ex['H1s'], n_gauss, fns, blocks)
        sto_shell(H, 1, 0, ex["H1s'"], n_gauss, fns, blocks)
        if polarization:
            sto_shell(H, 2, 1, ex['H2p'], ng_pol, fns, blocks)
    cols = []
    s = [i for i in range(nO) if fns[i].lmn == (0, 0, 0)]
    cs = np.asarray(ex['cs'])
    diffuse_s = s[int(np.argmin(ex['Os']))]
    for k in range(2):
        v = np.zeros(len(fns)); v[s] = cs[:, k]; cols.append(v)
    v = np.zeros(len(fns)); v[diffuse_s] = 1.0; cols.append(v)
    jd = int(np.argmin(ex['Op']))
    for c in CART[1]:
        p = [i for i in range(nO) if fns[i].lmn == c]
        v = np.zeros(len(fns)); v[p] = ex['cp']; cols.append(v)
        v = np.zeros(len(fns)); v[p[jd]] = 1.0; cols.append(v)
    rest = spherical_map(fns, [b for b in blocks if b[1] >= nO])
    return fns, np.hstack([np.array(cols).T, rest])


# --- exact one-centre radial average-of-configuration HF (s closed shells + one p shell) ---
# Slater-type primitives N r^(n-1) e^(-z r); every integral closed form.  Energy (Slater
# 1960 average of configuration): s-s pairs F0 - G0/2, s-p F0 - G1/6, p-p F0 - (2/25) F2.
from math import factorial as _fac


def _sto_norm(n, z):
    return (2 * z) ** (n + 0.5) / np.sqrt(_fac(2 * n))


def _rk(A, a, B, b, k):
    """int int r1^A e^(-a r1) r2^B e^(-b r2) r<^k / r>^(k+1) r1^2 r2^2 dr1 dr2."""
    p, q = A + 2, B + 2
    out = 0.0
    for (pp, aa, qq, bb) in ((p, a, q, b), (q, b, p, a)):
        # region r2 < r1: r1^(pp-k-1) * int_0^r1 r2^(qq+k)
        m = qq + k
        t = _fac(pp - k - 1) / aa ** (pp - k)
        t -= sum(bb ** j / _fac(j) * _fac(pp - k - 1 + j) / (aa + bb) ** (pp - k + j) for j in range(m + 1))
        out += _fac(m) / bb ** (m + 1) * t
    return out


def radial_atom(s_prims, p_prims, Z=8.0, q_p=4, n_s_closed=2, tol=1e-12, maxit=500):
    """s_prims, p_prims: lists of (n, zeta).  Returns (E, Cs [n_s x n_s_closed], cp)."""
    def one(prims, l):
        m = len(prims)
        S, H = np.zeros((m, m)), np.zeros((m, m))
        for i, (n1, z1) in enumerate(prims):
            for j, (n2, z2) in enumerate(prims):
                N = _sto_norm(n1, z1) * _sto_norm(n2, z2)
                G = lambda mm: _fac(mm) / (z1 + z2) ** (mm + 1)
                S[i, j] = N * G(n1 + n2)
                kin = (n1 - 1) * (n2 - 1) * G(n1 + n2 - 2) if n1 + n2 >= 3 else 0.0
                kin += -((n1 - 1) * z2 + (n2 - 1) * z1) * G(n1 + n2 - 1) + z1 * z2 * G(n1 + n2)
                kin += l * (l + 1) * G(n1 + n2 - 2)
                H[i, j] = N * (0.5 * kin - Z * G(n1 + n2 - 1))
        return S, H

    def R(pa, pb, pc, pd, k):
        """R^k(ab|cd) over primitive lists; densities a*b (r1), c*d (r2)."""
        out = np.zeros((len(pa), len(pb), len(pc), len(pd)))
        for i, (n1, z1) in enumerate(pa):
            for j, (n2, z2) in enumerate(pb):
                for u, (n3, z3) in enumerate(pc):
                    for v, (n4, z4) in enumerate(pd):
                        N = _sto_norm(n1, z1) * _sto_norm(n2, z2) * _sto_norm(n3, z3) * _sto_norm(n4, z4)
                        out[i, j, u, v] = N * _rk(n1 + n2 - 2, z1 + z2, n3 + n4 - 2, z3 + z4, k)
        return out

    sp, pp_ = list(s_prims), list(p_prims)
    Ss, Hs = one(sp, 0)
    Sp, Hp = one(pp_, 1)
    R0ss = R(sp, sp, sp, sp, 0)
    R0sp = R(sp, sp, pp_, pp_, 0)          # (s s | p p)
    R1x = R(sp, pp_, sp, pp_, 1)           # (s p | s p)
    R0pp = R(pp_, pp_, pp_, pp_, 0)
    R2pp = R(pp_, pp_, pp_, pp_, 2)

    def eig(F, S, nocc):
        w_, v_ = np.linalg.eigh(S)
        X = v_ @ np.diag(w_ ** -0.5) @ v_.T
        e, C = np.linalg.eigh(X.T @ F @ X)
        return (X @ C)[:, :nocc]

    Cs, cp = eig(Hs, Ss, n_s_closed), eig(Hp, Sp, 1)[:, 0]
    P, cc = Cs @ Cs.T, np.outer(cp, cp)
    Eold = 0.0
    for it in range(maxit):
        J0P = np.einsum('mnls,ls->mn', R0ss, P)
        K0P = np.einsum('mlns,ls->mn', R0ss, P)
        Fs = Hs + 2 * J0P - K0P + q_p * np.einsum('mnls,ls->mn', R0sp, cc) \
            - (q_p / 6) * np.einsum('mlns,ls->mn', R1x, cc)
        Fp = Hp + 2 * np.einsum('mnls,mn->ls', R0sp, P) - (1 / 3) * np.einsum('mlns,mn->ls', R1x, P) \
            + (q_p - 1) * (np.einsum('mnls,ls->mn', R0pp, cc) - 0.08 * np.einsum('mnls,ls->mn', R2pp, cc))
        E = 2 * np.sum(P * Hs) + q_p * np.sum(cc * Hp) \
            + 2 * np.sum(P * J0P) - np.sum(P * K0P) \
            + 2 * q_p * np.einsum('mn,mnls,ls->', P, R0sp, cc) - (q_p / 3) * np.einsum('mn,mlns,ls->', P, R1x, cc) \
            + q_p * (q_p - 1) / 2 * (np.einsum('mn,mnls,ls->', cc, R0pp, cc) - 0.08 * np.einsum('mn,mnls,ls->', cc, R2pp, cc))
        Cs, cp = eig(Fs, Ss, n_s_closed), eig(Fp, Sp, 1)[:, 0]
        Pn, ccn = Cs @ Cs.T, np.outer(cp, cp)
        if abs(E - Eold) < tol and it > 3 and np.max(np.abs(Pn - P)) < 1e-8:
            return E, Cs, cp
        Eold = E
        P, cc = 0.5 * (P + Pn), 0.5 * (cc + ccn)     # damping (Ne oscillates undamped)
    raise RuntimeError('radial atom did not converge')


# ================================================================ tests
CRAWFORD_NUC = [(np.array([0.0, -0.143225816552, 0.0]), 8.0),
                (np.array([1.638036840407, 1.136548822547, 0.0]), 1.0),
                (np.array([-1.638036840407, 1.136548822547, 0.0]), 1.0)]


def _basis(name):
    return os.path.join(ROOT, 'debug', 'data', 'basis', name)


def test_ccsd_t_matches_crawford_sto3g_water():
    """External reference: Crawford programming projects #5/#6 (STO-3G water, bohr
    geometry above): RHF -74.942079928192, CCSD corr -0.070680088376,
    (T) -0.000099877272.  Rejects: a (T) with a wrong permutation/weight or a dropped
    disconnected term; an SP shell read as two s functions (RHF -60.15)."""
    fns, T = gaussian_basis(_basis('sto-3g_H_O.json'), CRAWFORD_NUC)
    S, h, g = ao_integrals(fns, CRAWFORD_NUC, T)
    r = run(S, h, g, 10, vnn(CRAWFORD_NUC))
    assert abs(r['RHF'] - (-74.942079928192)) < 1e-6
    assert abs((r['CCSD'] - r['RHF']) - (-0.070680088376)) < 1e-8
    assert abs((r['CCSD(T)'] - r['CCSD']) - (-0.000099877272)) < 1e-10


def test_ccsd_equals_fci_for_two_electrons():
    """H2 / cc-pVDZ: CCSD is exact for two electrons.  Rejects a CCSD missing any
    doubles or singles-doubles coupling term."""
    nuc = [(np.zeros(3), 1.0), (np.array([0.0, 0.0, 1.4]), 1.0)]
    fns, T = gaussian_basis(_basis('cc-pvdz_H_O.json'), nuc)
    S, h, g = ao_integrals(fns, nuc, T)
    r = run(S, h, g, 2, vnn(nuc), triples=False)
    assert abs(r['CCSD'] - fci(S, h, g, 2, vnn(nuc))) < 1e-9


def test_radial_atom_he_exact():
    """Single-zeta He 1s, zeta = 27/16: E = -(27/16)^2 = -729/256 exactly."""
    E = radial_atom([(1, 27 / 16)], [(2, 1.0)], Z=2.0, q_p=0, n_s_closed=1)[0]
    assert abs(E + 729 / 256) < 1e-12


@pytest.mark.slow
def test_radial_atom_is_the_3d_limit_for_ne():
    """Independent route: closed-shell Ne in the 4s3p Slater primitive set, exact
    radial HF vs 3D RHF over Gaussian fits.  The 3D energy approaches the radial one
    from above (ng 8 -> 12: -4.9 -> -0.20 mHa).  Rejects wrong s-p (G1/6) or p-p
    (2F2/25) exchange weights and wrong kinetic l(l+1) terms."""
    zs, zp = [10.0, 7.0, 2.6, 1.6], [5.0, 2.4, 1.2]
    Er = radial_atom(list(zip(O_S_N, zs)), [(2, z) for z in zp], Z=10.0, q_p=6)[0]
    errs = []
    for ng in (8, 12):
        fns, _ = o_atom_primitives(np.zeros(3), zs, zp, ng)
        S, h, g = GE.integral_set_md(fns, [(np.zeros(3), 10.0)])
        errs.append(rhf(S, h, g, 10, 0.0)[0] - Er)
    assert 0 < errs[1] < errs[0]
    assert errs[1] < 5e-4


def test_frozen_core_ccsd_equals_active_space_fci():
    """STO-3G water with 4 frozen orbitals leaves 2 active electrons, so frozen-core
    CCSD must equal FCI in the active space with the core folded into an effective
    one-body operator (engine FCI, independent route).  Rejects frozen orbitals
    dropped from the Fock operator or a wrong active-electron count."""
    nuc = CRAWFORD_NUC
    fns, T = gaussian_basis(_basis('sto-3g_H_O.json'), nuc)
    S, h, g = ao_integrals(fns, nuc, T)
    fc = 4
    r = run(S, h, g, 10, vnn(nuc), triples=False, frozen_core=fc)
    Ehf, eps, C = rhf(S, h, g, 10, vnn(nuc))
    Cc, Ca = C[:, :fc], C[:, fc:]
    Dc = 2 * Cc @ Cc.T
    Jc = np.einsum('pqrs,rs->pq', g, Dc)
    Kc = np.einsum('prqs,rs->pq', g, Dc)
    Ecore = np.sum(Dc * (h + 0.5 * (Jc - 0.5 * Kc)))
    heff = Ca.T @ (h + Jc - 0.5 * Kc) @ Ca
    ga = np.einsum('pqrs,pa,qb,rc,sd->abcd', g, Ca, Ca, Ca, Ca, optimize=True)
    Efci = GE.fci_ground(heff, ga, 10 - 2 * fc) + Ecore + vnn(nuc)
    assert abs(r['CCSD'] - Efci) < 1e-9
