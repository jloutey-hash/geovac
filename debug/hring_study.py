"""Hydrogen-ring study (2026-11-01, PI-directed): Paper 62's equal-size Slater-vs-Gaussian test on a ring,
and the PI's prediction that a ring carries more than one usable label set (site j vs ring momentum k).

Ring: N hydrogens on a regular N-gon, nearest-neighbour distance r (Angstrom), in the xy plane.
Solver: the validated RHF + spin-orbital CCSD(T) of tests/test_water_ccsd.py (all electron; no core),
and a string-based FCI here (fci_strings) for small spaces, validated against GE.fci_ground.

Labellings of one orthonormal space (Part B):
  site : AOs re-expressed with LOCAL-frame p functions (radial, tangential, z), then symmetric (Lowdin)
         orthogonalization.  Each orbital sits on one site; the ring rotation permutes them.
  ring : the real discrete Fourier transform over sites of the site orbitals, per local type:
         k = 0, N/2 and cos/sin pairs for 0 < k < N/2.  Same space, an orthogonal rotation.
Measured per labelling: JW Pauli count and 1-norm at a tolerance (jw_pauli_stats), and the number of
Z-type Z2 symmetries (dimension of the F2 null space of the parity rows of every nonzero integral).

    python debug/hring_study.py fci_check                 # fci_strings vs GE.fci_ground, tiny ring
    python debug/hring_study.py h10_sto6g                 # pipeline check vs Stair-Evangelista FCI
    python debug/hring_study.py h2 sto|ccpvdz             # H2 control: p exponent on FCI, energy
    python debug/hring_study.py labels_control <N> <basis> # Part B measures on a non-scored case
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import LinearOperator, eigsh

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
_s = importlib.util.spec_from_file_location('w', os.path.join(ROOT, 'tests', 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(_s)
_s.loader.exec_module(w)
GE = w.GE
DATA = os.path.join(HERE, 'data')
BASIS = os.path.join(DATA, 'basis')
DL = os.path.join(HERE, 'lit_scan', 'dl_ring_bases')
OUT = os.path.join(DATA, 'hring_study.jsonl')
ANG = 1.0 / 0.52917721067
H_EX = {'H1s': 1.212, "H1s'": 0.627}          # H2 RHF at R = 1.4 (Paper 62's H rule, unchanged)
NG, NG_OPT, NGP_OPT = 12, 8, 6
GFILE = {'ccpvdz': 'cc-pvdz_H_O.json', 'ccpvtz': 'cc-pvtz_H_O.json', 'augccpvdz': 'aug-cc-pvdz_H_O_Ne.json'}


def write(row):
    print(json.dumps(row), flush=True)
    with open(OUT, 'a') as f:
        f.write(json.dumps(row) + '\n')


def ring(N, r_ang):
    r = r_ang * ANG
    Rc = r / (2 * np.sin(np.pi / N))
    return [(np.array([Rc * np.cos(2 * np.pi * k / N), Rc * np.sin(2 * np.pi * k / N), 0.0]), 1.0)
            for k in range(N)]


def h2(R=1.4):
    return [(np.zeros(3), 1.0), (np.array([0.0, 0.0, R]), 1.0)]


# ------------------------------------------------------------------ bases (H only)
def _h_only(path, pol=None):
    data = json.load(open(path))
    el = {'electron_shells': [sh for sh in data['elements']['1']['electron_shells']]}
    if pol is not None:
        for sh in el['electron_shells']:
            if sh['angular_momentum'][0] == 1 and len(sh['exponents']) == 1:
                sh['exponents'] = [str(pol)]
    return {'elements': {'1': el}}


def gaussian(name, nuc, pol=None):
    if name == 'sto6g':
        data = json.load(open(os.path.join(DL, 'sto-6g_H.json')))
        data = {'elements': {'1': {'electron_shells': [
            {'angular_momentum': sh['angular_momentum'], 'exponents': [str(float(x)) for x in sh['exponents']],
             'coefficients': [[str(float(c)) for c in col] for col in sh['coefficients']]}
            for sh in data['elements']['1']['electron_shells']]}}}
    elif name == 'ccpvdz_tzpol':
        dz = _h_only(os.path.join(BASIS, GFILE['ccpvdz']))
        tz = _h_only(os.path.join(BASIS, GFILE['ccpvtz']))
        s = [sh for sh in dz['elements']['1']['electron_shells'] if sh['angular_momentum'][0] == 0]
        p = [sh for sh in tz['elements']['1']['electron_shells'] if sh['angular_momentum'][0] == 1]
        data = {'elements': {'1': {'electron_shells': s + p}}}
    else:
        data = _h_only(os.path.join(BASIS, GFILE[name]), pol)
    path = os.path.join(BASIS, '_hring_tmp_%d.json' % os.getpid())
    json.dump(data, open(path, 'w'))
    out = w.gaussian_basis(path, nuc)
    os.remove(path)
    return out


def slater_dzp(nuc, p_exp, n_gauss=NG, ng_pol=NG, polarization=True):
    fns, blocks = [], []
    for C, _ in nuc:
        w.sto_shell(C, 1, 0, H_EX['H1s'], n_gauss, fns, blocks)
        w.sto_shell(C, 1, 0, H_EX["H1s'"], n_gauss, fns, blocks)
        if polarization:
            w.sto_shell(C, 2, 1, p_exp, ng_pol, fns, blocks)
    return fns, w.spherical_map(fns, blocks)


# ------------------------------------------------------------------ string FCI (small spaces)
def _strings(n, k):
    from itertools import combinations
    return [sum(1 << i for i in c) for c in combinations(range(n), k)]


def _excitation_ops(n, strs):
    """E_pq on one spin's string space as sparse matrices (p <- q), with fermionic signs."""
    idx = {s: i for i, s in enumerate(strs)}
    ops = {}
    for p in range(n):
        for q in range(n):
            rows, cols, vals = [], [], []
            for j, s in enumerate(strs):
                if not (s >> q) & 1:
                    continue
                t = s ^ (1 << q)
                if p != q and (t >> p) & 1:
                    continue
                sign = (-1) ** bin(s & ((1 << q) - 1)).count('1')
                t2 = t | (1 << p)
                sign *= (-1) ** bin(t & ((1 << p) - 1)).count('1')
                rows.append(idx[t2]); cols.append(j); vals.append(sign)
            ops[p, q] = csr_matrix((vals, (rows, cols)), shape=(len(strs), len(strs)))
    return ops


def fci_strings(h, g, na, nb, enuc=0.0, k=1):
    """FCI over an orthonormal spatial basis (chemist g), M_S = (na - nb)/2."""
    n = h.shape[0]
    sa, sb = _strings(n, na), _strings(n, nb)
    Ea, Eb = _excitation_ops(n, sa), _excitation_ops(n, sb)
    kmat = h - 0.5 * np.einsum('prrq->pq', g)
    shape = (len(sa), len(sb))

    def E(pq, C):
        return Ea[pq] @ C + (Eb[pq] @ C.T).T

    def mv(x):
        C = x.reshape(shape)
        D = np.array([[E((r, s), C) for s in range(n)] for r in range(n)])
        F = np.einsum('pqrs,rsab->pqab', g, D, optimize=True)
        sig = np.zeros(shape)
        for p in range(n):
            for q in range(n):
                sig += kmat[p, q] * D[p, q] + 0.5 * E((p, q), F[p, q])
        return sig.ravel()
    dim = shape[0] * shape[1]
    op = LinearOperator((dim, dim), matvec=mv, dtype=float)
    vals = eigsh(op, k=k, which='SA', tol=1e-12)[0]
    return np.sort(vals) + enuc


def orth_ints(S, h, g):
    X = GE.lowdin_orbitals(S)
    return GE.transform_integrals(X, h, g)


# ------------------------------------------------------------------ labellings (Part B)
def local_frame(fns, T, nuc):
    """Rotate each site's Cartesian p triple (columns of T in the spherical basis) to (radial,
    tangential, z).  Returns the new T and, per column, (site, type)."""
    T = T.copy()
    cols = T.shape[1]
    centers = [tuple(np.round(f.center, 10)) for f in fns]
    site_of = {tuple(np.round(C, 10)): j for j, (C, _) in enumerate(nuc)}
    lab = []
    c = 0
    while c < cols:
        rows = np.nonzero(np.abs(T[:, c]) > 0)[0]
        f = fns[rows[0]]
        j = site_of[tuple(np.round(f.center, 10))]
        if sum(f.lmn) == 1 and f.lmn == (1, 0, 0):
            th = np.arctan2(nuc[j][0][1], nuc[j][0][0])
            px, py, pz = T[:, c].copy(), T[:, c + 1].copy(), T[:, c + 2].copy()
            T[:, c] = np.cos(th) * px + np.sin(th) * py
            T[:, c + 1] = -np.sin(th) * px + np.cos(th) * py
            T[:, c + 2] = pz
            lab += [(j, ('p', c, 'rad')), (j, ('p', c, 'tan')), (j, ('p', c, 'z'))]
            c += 3
        else:
            lab.append((j, ('s', c)))
            c += 1
    # type key = position of the column within its site's list
    per_site = {}
    out = []
    for j, t in lab:
        per_site.setdefault(j, 0)
        out.append((j, per_site[j]))
        per_site[j] += 1
    return T, out


def real_dft(N):
    F = np.zeros((N, N))
    j = np.arange(N)
    F[:, 0] = 1 / np.sqrt(N)
    col = 1
    for k in range(1, (N - 1) // 2 + 1):
        F[:, col] = np.sqrt(2.0 / N) * np.cos(2 * np.pi * k * j / N); col += 1
        F[:, col] = np.sqrt(2.0 / N) * np.sin(2 * np.pi * k * j / N); col += 1
    if N % 2 == 0:
        F[:, col] = (-1.0) ** j / np.sqrt(N)
    return F


def labellings(S, h, g, lab, N):
    """S, h, g in the local-frame AO basis.  Returns {'site': (h, g), 'ring': (h, g)}."""
    X = GE.lowdin_orbitals(S)
    hs, gs = GE.transform_integrals(X, h, g)
    n = len(lab)
    U = np.zeros((n, n))
    F = real_dft(N)
    types = sorted(set(t for _, t in lab))
    col = 0
    for t in types:
        idx = [i for i, (j, tt) in enumerate(lab) if tt == t]
        idx = sorted(idx, key=lambda i: lab[i][0])
        assert [lab[i][0] for i in idx] == list(range(N))
        for kcol in range(N):
            U[idx, col] = F[:, kcol]
            col += 1
    hr, gr = GE.transform_integrals(U, hs, gs)
    return {'site': (hs, gs), 'ring': (hr, gr)}


def z2_count(h, g, tol=1e-10):
    """Dimension of the F2 space of Z-type strings commuting with H (spin-orbital P = 2p + sigma)."""
    M = h.shape[0]
    rows = set()
    for p in range(M):
        for q in range(M):
            if p != q and abs(h[p, q]) > tol:
                for s in (0, 1):
                    rows.add((1 << (2 * p + s)) ^ (1 << (2 * q + s)))
    nz = np.argwhere(np.abs(g) > tol)
    for p, q, r, s in nz:
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


def label_measures(h, g, enuc, tols=(1e-10, 1e-6, 1e-4, 1e-3, 1e-2)):
    out = {'z2': z2_count(h, g), 'exact_zero_eri_frac': float(np.mean(np.abs(g) <= 1e-10))}
    for t in tols:
        st = w.jw_pauli_stats(h, g, enuc, tol=t)
        out['pauli_%g' % t] = int(st['n_pauli'])
        if t == tols[0]:
            out['one_norm'] = float(st['one_norm'])
    return out


# ------------------------------------------------------------------ lean CCSD(T)
# The same equations as w.ccsd, but the antisymmetrized spin-orbital blocks are built from the spatial MO
# integrals block by block, and the vvvv block is never formed: its one use, 0.5 sum_ef tau_ij^ef W_abef, is
# split into the spatial (AE|BF) contraction plus two small terms.  w.spin_orbital forms the full (2n)^4 tensor
# several times over, which does not fit in memory at cc-pVTZ on H6 (168 spin orbitals).
def _so_block(gm, nocc, A, B, C, D):
    """<AB||CD> for spin-orbital index arrays (spin-orbital P = 2 p + s, as in w.spin_orbital)."""
    pa, pb, pc, pd = (x // 2 for x in (A, B, C, D))
    sa, sb, sc, sd = (x % 2 for x in (A, B, C, D))
    J = gm[np.ix_(pa, pc, pb, pd)].transpose(0, 2, 1, 3)
    J = J * (sa[:, None] == sc[None, :])[:, None, :, None] * (sb[:, None] == sd[None, :])[None, :, None, :]
    K = gm[np.ix_(pa, pd, pb, pc)].transpose(0, 2, 3, 1)
    K = K * (sa[:, None] == sd[None, :])[:, None, None, :] * (sb[:, None] == sc[None, :])[None, :, :, None]
    return J - K


def lean_ccsd(eps, gm, nocc, tol=1e-10, maxit=200, triples=True):
    es = w._es
    n2 = 2 * len(eps)
    O, V = np.arange(nocc), np.arange(nocc, n2)
    fd = np.repeat(eps, 2)
    eo, ev = fd[:nocc], fd[nocc:]
    nv = len(V)
    blk = lambda a, b, c, d: _so_block(gm, nocc, a, b, c, d)
    oovv, ovvo, ooov, vovv = blk(O, O, V, V), blk(O, V, V, O), blk(O, O, O, V), blk(V, O, V, V)
    oooo, ovvv, oovo, ovov, vvvo, ovoo = blk(O, O, O, O), blk(O, V, V, V), blk(O, O, V, O), blk(O, V, O, V), \
        blk(V, V, V, O), blk(O, V, O, O)
    nV = nv // 2
    gvv = np.ascontiguousarray(gm[nocc // 2:, nocc // 2:, nocc // 2:, nocc // 2:])
    fov = np.zeros((nocc, nv))

    def tau_wabef(tau, t1):
        # 0.5 sum_ef tau_ij^ef W_abef, W_abef = <ab||ef> - P(ab) t_m^b <am||ef> + 1/4 tau_mn^ab <mn||ef>
        t6 = tau.reshape(nocc, nocc, nV, 2, nV, 2)
        X = es('ijEsFt,AEBF->ijAsBt', t6, gvv).reshape(nocc, nocc, nv, nv)    # = 0.5 sum tau <ab||ef>
        Q = es('ijef,amef->ijam', tau, vovv)
        R = es('ijef,mnef->ijmn', tau, oovv)
        return (X + 0.5 * (-es('ijam,mb->ijab', Q, t1) + es('ijbm,ma->ijab', Q, t1))
                + 0.125 * es('mnab,ijmn->ijab', tau, R))

    Dia = eo[:, None] - ev[None, :]
    Dijab = eo[:, None, None, None] + eo[None, :, None, None] - ev[None, None, :, None] - ev[None, None, None, :]
    t1 = np.zeros_like(Dia)
    t2 = oovv / Dijab
    Eold = 0.0
    for it in range(maxit):
        tau_t = t2 + 0.5 * (es('ia,jb->ijab', t1, t1) - es('ib,ja->ijab', t1, t1))
        tau = t2 + es('ia,jb->ijab', t1, t1) - es('ib,ja->ijab', t1, t1)
        Fae = es('mf,mafe->ae', t1, ovvv) - 0.5 * es('mnaf,mnef->ae', tau_t, oovv)
        Fmi = es('ne,mnie->mi', t1, ooov) + 0.5 * es('inef,mnef->mi', tau_t, oovv)
        Fme = es('nf,mnef->me', t1, oovv)
        Wmnij = (oooo + es('je,mnie->mnij', t1, ooov) - es('ie,mnje->mnij', t1, ooov)
                 + 0.25 * es('ijef,mnef->mnij', tau, oovv))
        Wmbej = (ovvo + es('jf,mbef->mbej', t1, ovvv) - es('nb,mnej->mbej', t1, oovo)
                 - es('jnfb,mnef->mbej', 0.5 * t2 + es('jf,nb->jnfb', t1, t1), oovv))
        r1 = (fov + es('ie,ae->ia', t1, Fae) - es('ma,mi->ia', t1, Fmi)
              + es('imae,me->ia', t2, Fme) - es('nf,naif->ia', t1, ovov)
              - 0.5 * es('imef,maef->ia', t2, ovvv)
              - 0.5 * es('mnae,nmei->ia', t2, oovo))
        Ft = Fae - 0.5 * es('mb,me->be', t1, Fme)
        P1 = es('ijae,be->ijab', t2, Ft)
        Fo = Fmi + 0.5 * es('je,me->mj', t1, Fme)
        P2 = es('imab,mj->ijab', t2, Fo)
        X = es('imae,mbej->ijab', t2, Wmbej) - es('ie,ma,mbej->ijab', t1, t1, ovvo)
        Y = es('ie,abej->ijab', t1, vvvo)
        Zt = es('ma,mbij->ijab', t1, ovoo)
        r2 = (oovv.copy()
              + P1 - P1.transpose(0, 1, 3, 2)
              - P2 + P2.transpose(1, 0, 2, 3)
              + 0.5 * es('mnab,mnij->ijab', tau, Wmnij)
              + tau_wabef(tau, t1)
              + X - X.transpose(1, 0, 2, 3) - X.transpose(0, 1, 3, 2) + X.transpose(1, 0, 3, 2)
              + Y - Y.transpose(1, 0, 2, 3)
              - Zt + Zt.transpose(0, 1, 3, 2))
        t1 = r1 / Dia
        t2 = r2 / Dijab
        E = 0.25 * es('ijab,ijab->', oovv, t2) + 0.5 * es('ijab,ia,jb->', oovv, t1, t1)
        if abs(E - Eold) < tol and it > 2:
            break
        Eold = E
    else:
        raise RuntimeError('lean CCSD did not converge')
    if not triples:
        return E, 0.0
    ET = 0.0
    for i in range(nocc):
        for j in range(i + 1, nocc):
            for k in range(j + 1, nocc):
                def conn(i, j, k):
                    return es('ae,ebc->abc', t2[j, k], vovv[:, i]) - es('mbc,ma->abc', t2[i], ovoo[:, :, j, k])

                def disc(i, j, k):
                    return es('a,bc->abc', t1[i], oovv[j, k])

                def perm_abc(Xa):
                    return Xa - Xa.transpose(1, 0, 2) - Xa.transpose(2, 1, 0)

                tc = perm_abc(conn(i, j, k) - conn(j, i, k) - conn(k, j, i))
                td = perm_abc(disc(i, j, k) - disc(j, i, k) - disc(k, j, i))
                D3 = (eo[i] + eo[j] + eo[k] - ev[:, None, None] - ev[None, :, None] - ev[None, None, :])
                ET += np.sum(tc * (tc + td) / D3)
    return E, ET / 6.0


def lean_run(S, h, g, nelec, enuc, triples=True):
    Ehf, eps, C = w.rhf(S, h, g, nelec, enuc)
    gm = np.einsum('pqrs,pa,qb,rc,sd->abcd', g, C, C, C, C, optimize=True)
    Ec, Et = lean_ccsd(eps, gm, nelec, triples=triples)
    return {'RHF': Ehf, 'CCSD': Ehf + Ec, 'CCSD(T)': Ehf + Ec + Et}


# ------------------------------------------------------------------ modes
def mode_cc_check(N, basis):
    """lean_run against w.run on a case small enough for both (non-scored, or published values)."""
    nuc = ring(N, 1.0)
    fns, T = gaussian(basis, nuc) if basis != 'sto_min' else slater_dzp(nuc, None, NG, NG, polarization=False)
    S, h, g = w.ao_integrals(fns, nuc, T)
    enuc = w.vnn(nuc)
    a, b = w.run(S, h, g, N, enuc), lean_run(S, h, g, N, enuc)
    write({'mode': 'cc_check', 'N': N, 'basis': basis, 'nbf': T.shape[1], 'old': {k: float(v) for k, v in a.items()},
           'lean': {k: float(v) for k, v in b.items()}, 'maxdiff': float(max(abs(a[k] - b[k]) for k in a))})


def mode_tz_check():
    """The record's pipeline check: H6 cc-pVTZ CCSD against Magoulas-Shen-Piecuch (FCI -3.410069, CCSD +4.424 mEh)."""
    nuc = ring(N_RING, R_RING)
    t0 = time.time()
    fns, T = gaussian('ccpvtz', nuc)
    S, h, g = w.ao_integrals(fns, nuc, T)
    print('nbf', T.shape[1], 'ints', round(time.time() - t0, 1), 's', flush=True)
    r = lean_run(S, h, g, N_RING, w.vnn(nuc), triples=False)
    write({'mode': 'tz_check', 'nbf': T.shape[1], 'RHF': float(r['RHF']), 'CCSD': float(r['CCSD']),
           'ref_CCSD': -3.405645, 'diff': float(r['CCSD'] + 3.405645), 't': round(time.time() - t0, 1)})


def mode_fci_check():
    nuc = ring(4, 1.0)
    fns, T = gaussian('sto6g', nuc)
    S, h, g = w.ao_integrals(fns, nuc, T)
    ht, gt = orth_ints(S, h, g)
    e1 = GE.fci_ground(ht, gt, 4) + w.vnn(nuc)
    e2 = fci_strings(ht, gt, 2, 2, w.vnn(nuc))[0]
    write({'mode': 'fci_check', 'H4_sto6g_dense': float(e1), 'H4_sto6g_strings': float(e2), 'diff': float(e1 - e2)})


def mode_h10():
    nuc = ring(10, 1.0)
    fns, T = gaussian('sto6g', nuc)
    S, h, g = w.ao_integrals(fns, nuc, T)
    ht, gt = orth_ints(S, h, g)
    t0 = time.time()
    e = fci_strings(ht, gt, 5, 5, w.vnn(nuc))[0]
    r = w.run(S, h, g, 10, w.vnn(nuc))
    write({'mode': 'h10_sto6g', 'FCI': float(e), 'ref_stair_evangelista': -5.422958, 'diff': float(e + 5.422958),
           **{k: float(v) for k, v in r.items()}, 't': round(time.time() - t0, 1)})


def h2_energy(kind, p, n_gauss=NG, ng_pol=NG):
    nuc = h2()
    fns, T = slater_dzp(nuc, p, n_gauss, ng_pol) if kind == 'sto' else gaussian('ccpvdz', nuc, pol=p)
    S, h, g = w.ao_integrals(fns, nuc, T)
    r = w.run(S, h, g, 2, w.vnn(nuc), triples=False)
    return r


def mode_h2(kind):
    f = lambda lz: h2_energy(kind, float(np.exp(lz)), NG_OPT, NGP_OPT)['CCSD']
    res = minimize_scalar(f, bounds=(np.log(0.3), np.log(4.0)), method='bounded', options={'xatol': 1e-3})
    p = float(np.exp(res.x))
    r = h2_energy(kind, p)
    rp = h2_energy(kind, None if kind == 'ccpvdz' else 1.0)
    write({'mode': 'h2', 'basis': kind, 'p_opt': p, 'nfev': int(res.nfev),
           'RHF': float(r['RHF']), 'FCI': float(r['CCSD']),
           'FCI_at_p_default': float(rp['CCSD']), 'p_default': 0.727 if kind == 'ccpvdz' else 1.0})


def mode_labels_control(N, basis):
    nuc = ring(N, 1.0)
    fns, T = gaussian(basis, nuc) if basis != 'sto_min' else slater_dzp(nuc, None, NG, NG, polarization=False)
    T, lab = local_frame(fns, T, nuc)
    S, h, g = w.ao_integrals(fns, nuc, T)
    L = labellings(S, h, g, lab, N)
    enuc = w.vnn(nuc)
    row = {'mode': 'labels_control', 'N': N, 'basis': basis, 'nbf': T.shape[1]}
    for name, (hh, gg) in L.items():
        row[name] = label_measures(hh, gg, enuc)
    if T.shape[1] <= 12:
        row['fci_site'] = float(fci_strings(*L['site'], N // 2, N // 2, enuc)[0])
        row['fci_ring'] = float(fci_strings(*L['ring'], N // 2, N // 2, enuc)[0])
        row['cc'] = {k: float(v) for k, v in w.run(S, h, g, N, enuc).items()}
    write(row)


# ------------------------------------------------------------------ scored modes (record HR; run after commit)
EXP = os.path.join(DATA, 'hring_exponents.json')
N_RING, R_RING = 6, 1.0


def _load_exp():
    return json.load(open(EXP)) if os.path.exists(EXP) else {}


def ring_basis(name, nuc, n_gauss=NG, ng_pol=NG):
    d = _load_exp()
    if name == 'sto_dzp':
        return slater_dzp(nuc, d['sto_p'], n_gauss, ng_pol)
    if name == 'ccpvdz_polopt':
        return gaussian('ccpvdz', nuc, pol=d['gto_p'])
    return gaussian(name, nuc)


def mode_pol(kind):
    nuc = ring(N_RING, R_RING)
    enuc = w.vnn(nuc)

    def f(lz):
        p = float(np.exp(lz))
        fns, T = slater_dzp(nuc, p, NG_OPT, NGP_OPT) if kind == 'sto' else gaussian('ccpvdz', nuc, pol=p)
        S, h, g = w.ao_integrals(fns, nuc, T)
        E = lean_run(S, h, g, 2 * (N_RING // 2), enuc, triples=False)['CCSD']
        print('  pol', kind, '%.5f' % p, '%.8f' % E, flush=True)
        return E
    res = minimize_scalar(f, bounds=(np.log(0.3), np.log(4.0)), method='bounded',
                          options={'xatol': 1e-3, 'maxiter': 60})
    d = _load_exp()
    d.update({kind + '_p': float(np.exp(res.x)), kind + '_p_E_ccsd': float(res.fun),
              kind + '_p_nfev': int(res.nfev), kind + '_p_success': bool(res.success)})
    json.dump(d, open(EXP, 'w'), indent=1)
    print(d, flush=True)


def _canonical(S, h, g, ne, enuc):
    _, eps, C = w.rhf(S, h, g, ne, enuc)
    hm = C.T @ h @ C
    gm = np.einsum('pqrs,pa,qb,rc,sd->abcd', g, C, C, C, C, optimize=True)
    return hm, gm


def mode_final(name):
    nuc = ring(N_RING, R_RING)
    enuc = w.vnn(nuc)
    t0 = time.time()
    fns, T = ring_basis(name, nuc)
    S, h, g = w.ao_integrals(fns, nuc, T)
    r = lean_run(S, h, g, N_RING, enuc)
    row = {'mode': 'final', 'basis': name, 'nbf': T.shape[1], **{k: float(v) for k, v in r.items()}}
    hm, gm = _canonical(S, h, g, N_RING, enuc)
    st = w.jw_pauli_stats(hm, gm, enuc)
    row['qubit'] = {'qubits': 2 * hm.shape[0], 'n_pauli': int(st['n_pauli']), 'one_norm': float(st['one_norm'])}
    row['t'] = round(time.time() - t0, 1)
    write(row)


def mode_labels(name):
    nuc = ring(N_RING, R_RING)
    fns, T = ring_basis(name, nuc)
    T, lab = local_frame(fns, T, nuc)
    S, h, g = w.ao_integrals(fns, nuc, T)
    L = labellings(S, h, g, lab, N_RING)
    row = {'mode': 'labels', 'basis': name, 'nbf': T.shape[1]}
    for k, (hh, gg) in L.items():
        row[k] = label_measures(hh, gg, w.vnn(nuc))
    write(row)


if __name__ == '__main__':
    cmd = sys.argv[1]
    t0 = time.time()
    if cmd == 'fci_check':
        mode_fci_check()
    elif cmd == 'h10_sto6g':
        mode_h10()
    elif cmd == 'h2':
        mode_h2(sys.argv[2])
    elif cmd == 'h2_gauss':
        nuc = h2()
        fns, T = gaussian(sys.argv[2], nuc)
        S, h, g = w.ao_integrals(fns, nuc, T)
        r = w.run(S, h, g, 2, w.vnn(nuc), triples=False)
        write({'mode': 'h2_gauss', 'basis': sys.argv[2], 'nbf': T.shape[1], 'RHF': float(r['RHF']),
               'FCI': float(r['CCSD'])})
    elif cmd == 'cc_check':
        mode_cc_check(int(sys.argv[2]), sys.argv[3])
    elif cmd == 'tz_check':
        mode_tz_check()
    elif cmd == 'labels_control':
        mode_labels_control(int(sys.argv[2]), sys.argv[3])
    elif cmd == 'pol':
        mode_pol(sys.argv[2])
    elif cmd == 'final':
        mode_final(sys.argv[2])
    elif cmd == 'labels':
        mode_labels(sys.argv[2])
    print(cmd, 'done in %.0f s' % (time.time() - t0), flush=True)
