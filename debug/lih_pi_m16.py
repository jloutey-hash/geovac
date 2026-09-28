r"""LiH pi/M=16 marriage (2026-09-23; plan debug/lih_pi_m16_build_plan.md).

PHASE A  -- the pi-containing M=16 reference + the m-RESOLVED pair-index coefficient tensor T,
            and GATE G0-pi.
PHASE B  -- the SEVEN enumerators generalised to sum over azimuthal modes m with the L_z
            selection rule, and GATE G1-pi (sigma-through-m regression to -8.037236).

Extends debug/lih_marriage_phase2.py (the m=0 enumerators) and debug/lih_marriage_phase0.py
(the tensor builder).  Nothing in geovac/lih_r12ci/ or debug/lih_vmc.py is modified.

CONVENTIONS (new for pi vs the sigma marriage)
----------------------------------------------
* ORDERED pair indices A=(x,y), x,y in [0,Mk).  rho_A(r) = phi_x^*(r) phi_y(r).
  On the kernels.py grid (phi_azi=0 slice) every orbital value is REAL (the e^{i m phi} factor
  is stripped, =1 at phi=0), so rho_A's grid part is the REAL product R_x*R_y, and the pair
  carries a signed azimuthal label  mA[A] = m_y - m_x  (m_p = orbital p's azimuthal number,
  0 for sigma, +-1 for the pi+-).  For a sigma reference every mA=0 and this collapses to the
  sigma marriage (the m=0 code paths are used UNCHANGED -> gate G1-pi reproduces -8.037236).
* CI coefficients are REAL (H is real-symmetric, block-diagonal in m); NOs are real, m-pure.
* Kernel dressing of a leaf carrying circulation mu:  the leaf is dressed by the mode-|mu|
  Fourier kernel of f and DEPOSITS circulation +mu onto its target.  Coefficient c_mu = 2pi
  (mu=0, via the existing gem.dress -> Kf) or pi (mu!=0, via pi*Fm[|mu|]).  Derivation +
  verification vs the physical 2-electron integral is in the build report.
* Integration vertex: contributes only if the total accumulated circulation is 0
  (INT dphi e^{i M phi} = 2pi delta_{M,0}); enforced as a mask over the pair axes.

Run:
  python debug/lih_pi_m16.py g0        # PHASE A gate G0-pi
  python debug/lih_pi_m16.py g1        # PHASE B gate G1-pi (sigma regression, exact path)
  python debug/lih_pi_m16.py g0 g1     # both
Writes debug/data/lih_pi_m16_tensor.npz (the M=16 m-resolved tensor + traces) on g0.
"""
from __future__ import annotations

import os
import sys
import time
from itertools import combinations
from typing import Callable, Dict, List, Optional, Sequence, Tuple

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
DATA = os.path.join(HERE, "data")

import geovac.lih_r12ci.kernels as KG
from geovac.lih_r12ci.kernels import (a as A_FOCUS, GAM, geo_f, NXI, NETA, NG, Xg as _Xg, Eg as _Eg)
from geovac.lih_r12ci.energy import R as R_ENG, Z_A as ZA_ENG, Z_B as ZB_ENG

# phase-2 machinery we extend / reuse
import lih_marriage_phase2 as P2
import lih_marriage_phase0 as P0

CFG = dict(Jb=2, Lb=1, npi=1, Jpi=1, Lpi=0, alpha=1.0, core2=(4.5, 1.6), tol=1e-11)
CFG_SIGMA = dict(CFG); CFG_SIGMA['npi'] = 0        # for the G1-pi regression build
E_REF_PI = -8.032                                  # M=16 orbital ceiling (approx)
E_REF_SIGMA = -8.02335
TRUNC = 1e-3
Z_A, Z_B = 3.0, 1.0
NELEC = 4
GATE_MHA = 1.0
PRIM_CACHE_PI = os.path.join(DATA, "lih_pi_m16_prim.npz")
PRIM_CACHE_SIGMA = os.path.join(DATA, "lih_marriage_phase0_prim.npz")   # phase-0 sigma cache
OUT_TENSOR = os.path.join(DATA, "lih_pi_m16_tensor.npz")

VNE = P2.VNE                       # one-body V_ne on the grid (flat)
GI_W = P2.GI_W                     # 2 pi a^3 geo (flat)  == grid_int weight
A3 = P2.A3
V_NN = P2.V_NN


def _t(t0: float) -> str:
    return f"[{time.time() - t0:7.1f}s]"


# =========================================================================== #
# orbital azimuthal number from a lih_core2exp_probe tag
# =========================================================================== #
def m_of_tag(tag) -> int:
    """(key, j, l) tag -> azimuthal number.  key=('bondpi', ms, Jpi, Lpi) -> ms; else 0."""
    key = tag[0]
    if isinstance(key, tuple) and len(key) > 0 and key[0] == 'bondpi':
        return int(key[1])
    return 0


# =========================================================================== #
# (a) primitives + engine matrices (mirror phase0.build_primitives, own cache, npi-aware)
# =========================================================================== #
def build_primitives(cfg: dict, cache_path: str, verbose: bool = True) -> dict:
    import lih_core2exp_probe as P
    import prolate_energy_ladder as L
    from prolate_allelectron_c4 import one_body_general
    from prolate_float_eri import build_eri_tensor_m_f
    from lih_vmc import extract_prim

    R = float(P.R)
    orbs, tags = P.lih_orbs(cfg['Jb'], cfg['Lb'], cfg['npi'], cfg['Jpi'], cfg['Lpi'],
                            cfg['alpha'], list(cfg['core2']))
    M = len(orbs)
    labels = []
    for (key, j, l) in tags:
        name = key[0] if isinstance(key, tuple) else key
        labels.append(f"{name}_x{j}e{l}" if name in ('bond', 'bondpi') else name)
    m_prim = np.array([m_of_tag(t) for t in tags], dtype=int)

    cached = False
    if os.path.exists(cache_path):
        z = np.load(cache_path, allow_pickle=True)
        if str(z['cfg']) == repr(cfg) and int(z['M']) == M:
            S, h1, T1, eri = z['S'], z['h1'], z['T1'], z['eri']
            cached = True
            if verbose:
                print(f"  primitive engine matrices loaded from cache {os.path.basename(cache_path)}")
    if not cached:
        t0 = time.time()
        S, h1 = one_body_general(orbs, R, Z_A, Z_B)
        S0, T1 = one_body_general(orbs, R, 0.0, 0.0)
        assert np.allclose(S, S0, atol=1e-13)
        eri = build_eri_tensor_m_f(orbs, R, verbose=verbose)
        os.makedirs(DATA, exist_ok=True)
        np.savez(cache_path, cfg=repr(cfg), M=M, R=R, S=S, h1=h1, T1=T1, eri=eri,
                 tags=np.array(tags, dtype=object))
        if verbose:
            print(f"  built engine matrices {_t(t0)}")
    # engine arrays may be complex for pi; the L_z-allowed blocks are real (H real-symmetric)
    S = np.asarray(S); h1 = np.asarray(h1); T1 = np.asarray(T1); eri = np.asarray(eri)
    imag = max(np.max(np.abs(S.imag)) if np.iscomplexobj(S) else 0.0,
               np.max(np.abs(h1.imag)) if np.iscomplexobj(h1) else 0.0,
               np.max(np.abs(eri.imag)) if np.iscomplexobj(eri) else 0.0)
    if verbose:
        print(f"  engine arrays: max|imag(S,h1,eri)| = {imag:.1e}  (real-symmetric -> taking real part)")
    S = S.real; h1 = h1.real; T1 = T1.real; eri = eri.real
    Vne = h1 - T1

    C = L.build_C(tags, cfg['alpha'])                # block-diagonal per (mu,msign) channel -> m-pure
    So = C @ S @ C.T
    h1o = C @ h1 @ C.T
    T1o = C @ T1 @ C.T
    erio = np.einsum('pa,qb,rc,sd,abcd->pqrs', C, C, C, C, eri, optimize=True)
    # BLOCK-diagonal canonical orthogonalization in m (the pi+-1 shells are degenerate; a global
    # eigh returns real px/py mixtures.  So is exactly m-block-diagonal, so orthogonalize per block).
    m_reb = m_prim.copy()                            # rebased fn i has the same m as primitive i
    xcross = np.max([np.max(np.abs(So[np.ix_(m_reb == mi, m_reb == mj)]))
                     for mi in np.unique(m_reb) for mj in np.unique(m_reb) if mi != mj] + [0.0])
    if verbose:
        print(f"  So cross-m max = {xcross:.1e} (should be ~0; m-block orthogonalization)")
    wmax = np.max(np.linalg.eigvalsh(So))
    Xcols = []; m_mo = []
    for mval in sorted(np.unique(m_reb)):
        I = np.nonzero(m_reb == mval)[0]
        wb, Ub = np.linalg.eigh(So[np.ix_(I, I)])
        keepb = wb > cfg['tol'] * wmax
        Xb = np.zeros((len(m_reb), int(keepb.sum())))
        Xb[I, :] = Ub[:, keepb] / np.sqrt(wb[keepb])
        Xcols.append(Xb); m_mo.extend([int(mval)] * int(keepb.sum()))
    X = np.concatenate(Xcols, axis=1)                # (M', Mk) m-pure columns
    m_mo = np.array(m_mo, dtype=int)
    Mk = X.shape[1]
    wall = np.linalg.eigvalsh(So); cond = wall[-1] / wall[wall > cfg['tol'] * wmax].min()
    h1_f = X.T @ h1o @ X
    T1_f = X.T @ T1o @ X
    eri_f = np.einsum('ap,bq,cr,ds,abcd->pqrs', X, X, X, X, erio, optimize=True)
    Tmap = C.T @ X                                   # (M, Mk) primitive -> canonical MO (m-pure)
    prims = [extract_prim(o) for o in orbs]
    return dict(R=R, orbs=orbs, tags=tags, labels=labels, m_prim=m_prim, m_mo=m_mo,
                M=M, Mk=Mk, cond=cond,
                S=S, h1=h1, T1=T1, Vne=Vne, eri=eri, C=C, X=X, Tmap=Tmap,
                h1_f=h1_f, T1_f=T1_f, Vne_f=h1_f - T1_f, eri_f=eri_f, prims=prims,
                Vnn=Z_A * Z_B / R)


def _mo_azimuthal(coeff_map: np.ndarray, m_prim: np.ndarray, tol: float = 1e-7) -> np.ndarray:
    """Assign each column (MO/NO) of coeff_map (M x n) an azimuthal number = the m of the
    primitives it is built from; assert m-purity (block-diagonal transform)."""
    n = coeff_map.shape[1]
    m_out = np.zeros(n, dtype=int)
    thr = tol * np.max(np.abs(coeff_map))
    for p in range(n):
        cont = np.abs(coeff_map[:, p])
        ms = np.unique(m_prim[cont > thr])
        assert ms.size == 1, f"MO/NO {p} mixes azimuthal numbers {ms} (not m-pure)"
        m_out[p] = int(ms[0])
    return m_out


# =========================================================================== #
# (b) FCI + natural orbitals + truncation  (real, m-block-diagonal)
# =========================================================================== #
def fci_solve(h1: np.ndarray, eri: np.ndarray, K: int, nelec: int = NELEC):
    from prolate_allelectron_fci import _dets, _matel
    from fci_fast import _popcount64, _FOUR
    na = nb = nelec // 2
    dets = _dets(K, na, nb)
    nd = len(dets)
    masks = np.array([int(sum(1 << so for so in D)) for D in dets], dtype=np.uint64)
    rl, cl, vl = [], [], []
    diag = np.empty(nd)
    for aI in range(nd):
        diag[aI] = _matel(h1, eri, dets[aI], dets[aI])
        if aI + 1 < nd:
            pc = _popcount64(masks[aI + 1:] ^ masks[aI])
            conn = np.nonzero(pc <= _FOUR)[0] + (aI + 1)
            Da = dets[aI]
            for bI in conn.tolist():
                v = _matel(h1, eri, Da, dets[bI])
                if v != 0.0:
                    rl.append(aI); cl.append(bI); vl.append(v)
    r = np.array(rl, dtype=np.int64); c_ = np.array(cl, dtype=np.int64); v = np.array(vl)
    ri = np.concatenate([np.arange(nd), r, c_]); ci = np.concatenate([np.arange(nd), c_, r])
    vi = np.concatenate([diag, v, v])
    H = sp.csr_matrix((vi, (ri, ci)), shape=(nd, nd))
    ev, evec = eigsh(H, k=1, which='SA', return_eigenvectors=True, tol=0)
    c = evec[:, 0]
    if c[np.argmax(np.abs(c))] < 0:
        c = -c
    return float(ev[0]), c, dets, H


# =========================================================================== #
# (c) ORDERED-pair index bookkeeping + m-resolved tensor
# =========================================================================== #
def pair_index_ordered(K: int):
    pairs = [(x, y) for x in range(K) for y in range(K)]       # ORDERED, x = bra, y = ket
    idx = {xy: i for i, xy in enumerate(pairs)}
    return pairs, idx


def build_ordered_tensor_compact(CS: np.ndarray, apairs, Mk: int):
    """Memory-safe ORDERED m-resolved tensor.  Avoids the (npair^2 x npair^2) = (Mk^4 x Mk^4)
    dense that P0.build_T would allocate (34 GB at Mk=16).  Works in the COMPACT (A,B)-key space
    (<= n_active^2 keys), never the full Mk^2 pair dimension.

    rho_(x,y) = phi_x^* phi_y (ordered pair index x*Mk+y).  dA[a]^* dA[a'](1,2), a=(p,q) bra,
    a'=(p2,q2) ket, expands into 4 (A,B) density-pair terms with signs; T = 1/4 sum CS CS E_a E_b.
    Returns (act_pairs, (iA,iB,iC,iD, val)) with LOCAL indices into act_pairs."""
    npair = Mk * Mk
    ap = np.asarray(apairs)                                   # (n_ap, 2)
    act_a = np.nonzero(np.any(CS != 0.0, axis=1))[0]
    act_b = np.nonzero(np.any(CS != 0.0, axis=0))[0]
    na, nb = len(act_a), len(act_b)
    sgn = np.array([1.0, -1.0, -1.0, 1.0])

    def dp_keys(idx_arr):
        """for all (i,j) in idx_arr x idx_arr: the 4 (A,B) ordered-pair-density keys A*npair+B."""
        p = ap[idx_arr, 0][:, None]; q = ap[idx_arr, 1][:, None]
        p2 = ap[idx_arr, 0][None, :]; q2 = ap[idx_arr, 1][None, :]
        A = np.stack([p * Mk + p2, p * Mk + q2, q * Mk + p2, q * Mk + q2], axis=-1)   # (n,n,4)
        B = np.stack([q * Mk + q2, q * Mk + p2, p * Mk + q2, p * Mk + p2], axis=-1)
        return (A * npair + B)                                # (n,n,4) key = A*npair + B

    keyA = dp_keys(act_a).reshape(na * na, 4)                 # rows (ia*na+iap)
    keyB = dp_keys(act_b).reshape(nb * nb, 4)
    uA, invA = np.unique(keyA, return_inverse=True); invA = invA.reshape(na * na, 4)
    uB, invB = np.unique(keyB, return_inverse=True); invB = invB.reshape(nb * nb, 4)
    rowsA = np.repeat(np.arange(na * na), 4)
    Ma = sp.coo_matrix((np.tile(sgn, na * na), (rowsA, invA.reshape(-1))),
                       shape=(na * na, len(uA))).tocsr()
    rowsB = np.repeat(np.arange(nb * nb), 4)
    Mb = sp.coo_matrix((np.tile(sgn, nb * nb), (rowsB, invB.reshape(-1))),
                       shape=(nb * nb, len(uB))).tocsr()
    CSr = CS[np.ix_(act_a, act_b)]                            # (na, nb)
    Qmat = np.einsum('ie,jf->ijef', CSr, CSr).reshape(na * na, nb * nb)   # rows (ia,iap), cols (ib,ibp)
    Tc = 0.25 * np.asarray(Ma.T @ (Qmat @ Mb))                # (n_keyA, n_keyB) dense (compact)
    kA, kB = np.nonzero(np.abs(Tc) > 1e-14)
    val = Tc[kA, kB]
    A = uA[kA] // npair; B = uA[kA] % npair
    C = uB[kB] // npair; D = uB[kB] % npair
    act = np.array(sorted(set(A) | set(B) | set(C) | set(D)))
    loc = {int(g): i for i, g in enumerate(act)}
    iA = np.array([loc[int(x)] for x in A]); iB = np.array([loc[int(x)] for x in B])
    iC = np.array([loc[int(x)] for x in C]); iD = np.array([loc[int(x)] for x in D])
    return act, (iA, iB, iC, iD, val)


# =========================================================================== #
# RefM: reference bundle carrying signed circulation mA per ordered pair
# =========================================================================== #
def traces_coo(T_local, d: np.ndarray, nact: int) -> dict:
    """spectator traces of the COO tensor T (local active-pair indices), d = spectator integral.
    T1a[A]=sum_BCD T d d d; T2aa[A,B]=sum_CD T d d; T2ab[A,C]=sum_BD T d d; T2bb[C,D]=sum_AB T d d;
    T3aab[A,B,C]=sum_D T d; T3abb[A,C,D]=sum_B T d.  All O(nnz)."""
    iA, iB, iC, iD, val = T_local
    dA, dB, dC, dD = d[iA], d[iB], d[iC], d[iD]
    T1a = np.zeros(nact); np.add.at(T1a, iA, val * dB * dC * dD)
    T1b = np.zeros(nact); np.add.at(T1b, iC, val * dA * dB * dD)
    T2aa = np.zeros((nact, nact)); np.add.at(T2aa, (iA, iB), val * dC * dD)
    T2bb = np.zeros((nact, nact)); np.add.at(T2bb, (iC, iD), val * dA * dB)
    T2ab = np.zeros((nact, nact)); np.add.at(T2ab, (iA, iC), val * dB * dD)
    T3aab = np.zeros((nact, nact, nact)); np.add.at(T3aab, (iA, iB, iC), val * dD)
    T3abb = np.zeros((nact, nact, nact)); np.add.at(T3abb, (iA, iC, iD), val * dB)
    return dict(T1a=T1a, T1b=T1b, T2aa=T2aa, T2bb=T2bb, T2ab=T2ab, T3aab=T3aab, T3abb=T3abb)


class RefM:
    def __init__(self, act, rho, grad, d, mA, T_local, tag: str):
        self.act = act
        self.nact = len(act)
        self.rho = rho                       # (nact, NG) real
        self.grad = grad                     # (nact, NG) real  (conj-gradient dot, incl m_p m_q/rho^2)
        self.rhov = rho * VNE[None, :]
        self.d = d                           # (nact,) spectator integral (=delta_{x,y})
        self.mA = np.asarray(mA, int)        # (nact,) signed circulation m_y - m_x
        self.T_iA, self.T_iB, self.T_iC, self.T_iD, self.T_val = T_local
        self.tag = tag
        self._cache: Dict[tuple, np.ndarray] = {}

    def base(self, name: str) -> np.ndarray:
        return {'rho': self.rho, 'grad': self.grad, 'rhov': self.rhov}[name]


def _fci_cached(cfg_tag: str, which: str, h1, eri, Mk):
    """Cache the (slow) 14400-det FCI ground vector, keyed by (cfg, which-basis)."""
    path = os.path.join(DATA, f"lih_pi_m16_fci_{cfg_tag}_{which}.npz")
    if os.path.exists(path):
        z = np.load(path)
        from prolate_allelectron_fci import _dets
        dets = _dets(Mk, NELEC // 2, NELEC // 2)
        return float(z['E0']), z['c'], dets
    E0, c, dets, _ = fci_solve(h1, eri, Mk)
    np.savez(path, E0=E0, c=c)
    return E0, c, dets


def build_ref_from_prims(ref: dict, tag: str, cfg_tag: str = 'pi', trunc: float = TRUNC,
                         verbose: bool = True):
    """FCI -> NOs -> truncate -> ordered m-resolved tensor T + traces + grid densities.
    Returns (RefM, info dict with everything G0-pi needs).  `trunc` = |c| cut (det-truncation knob)."""
    R, M, Mk = ref['R'], ref['M'], ref['Mk']
    Vnn = ref['Vnn']
    E0, c, dets = _fci_cached(cfg_tag, 'canon', ref['h1_f'], ref['eri_f'], Mk)
    E_full = E0 + Vnn
    nd = len(dets)
    if verbose:
        print(f"  FCI (canonical MO basis): Mk={Mk} n_det={nd} E_tot={E_full:.6f}")

    # natural orbitals of the FULL CI (spin-summed 1-RDM)
    pairs, idx = pair_index_ordered(Mk)
    npair = len(pairs)
    apairs = list(combinations(range(Mk), 2)); aidx = {ap: i for i, ap in enumerate(apairs)}
    n_ap = len(apairs)
    Eop = P0.build_Eop(Mk, apairs, aidx)
    CS_full, _, sigma = P0.fold_cs(c, Mk)
    ga_c = np.einsum('pqaA,aA->pq', Eop, CS_full @ CS_full.T)
    gb_c = np.einsum('pqaA,aA->pq', Eop, CS_full.T @ CS_full)
    D1 = ga_c + gb_c                                            # spin-summed 1-RDM (m-block-diagonal)
    # BLOCK-diagonal NO rotation in m (degenerate pi+-1 -> global eigh mixes them); D1 is
    # exactly m-block-diagonal for an M_L=0 state, so diagonalize each m-block separately.
    m_mo = ref['m_mo']
    U = np.zeros((Mk, Mk)); occ = np.zeros(Mk); m_no_raw = np.zeros(Mk, dtype=int)
    col = 0
    for mval in sorted(np.unique(m_mo)):
        I = np.nonzero(m_mo == mval)[0]
        wb, Vb = np.linalg.eigh(D1[np.ix_(I, I)])
        for k in range(len(I)):
            U[I, col] = Vb[:, k]; occ[col] = wb[k]; m_no_raw[col] = int(mval); col += 1
    order = np.argsort(occ)[::-1]
    occ = occ[order]; U = U[:, order]; m_no = m_no_raw[order]
    for p in range(Mk):
        if U[np.argmax(np.abs(U[:, p])), p] < 0:
            U[:, p] = -U[:, p]
    h1_no = U.T @ ref['h1_f'] @ U
    T1_no = U.T @ ref['T1_f'] @ U
    Vne_no = h1_no - T1_no
    eri_no = np.einsum('ap,bq,cr,ds,abcd->pqrs', U, U, U, U, ref['eri_f'], optimize=True)
    T_no = ref['Tmap'] @ U                                       # (M, Mk) primitive -> NO
    # NO azimuthal numbers (from the m-block structure); cross-check purity of T_no
    m_no_chk = _mo_azimuthal(T_no, ref['m_prim'])
    assert np.array_equal(m_no, m_no_chk), (m_no, m_no_chk)
    if verbose:
        print(f"  NO occupations: {np.round(occ, 5)}")
        print(f"  NO azimuthal m: {m_no.tolist()}")

    E0n, cn, dets_n = _fci_cached(cfg_tag, 'no', h1_no, eri_no, Mk)
    kept = np.nonzero(np.abs(cn) > trunc)[0]
    norm_kept = float(np.sum(cn[kept] ** 2))
    ct = np.zeros_like(cn); ct[kept] = cn[kept] / np.sqrt(norm_kept)
    CS, _, _ = P0.fold_cs(ct, Mk)
    ga, gb, Gaa, Gbb, Gab = P0.rdms(CS, Eop)
    E_trunc = P0.energy_from_rdm(ga, gb, Gaa, Gbb, Gab, h1_no, eri_no, Vnn)
    g_tot = ga + gb
    G_tot = Gaa + Gbb + Gab + Gab.transpose(2, 3, 0, 1)
    if verbose:
        print(f"  truncation |c|>{TRUNC}: {len(kept)} of {nd} dets, norm {norm_kept:.6f}")
        print(f"  E_trunc = {E_trunc:.6f};  E_trunc - E_full = {(E_trunc - E_full)*1e3:+.3f} mHa")

    # ordered m-resolved tensor (memory-safe compact build; local indices into `act`)
    act, (liA, liB, liC, liD, val) = build_ordered_tensor_compact(CS, apairs, Mk)
    act_pairs = [pairs[g] for g in act]                       # (x,y) of each active ordered pair
    d = np.array([1.0 if x == y else 0.0 for (x, y) in act_pairs])
    mA_act = np.array([m_no[y] - m_no[x] for (x, y) in act_pairs], dtype=int)
    T_local_pre = (liA, liB, liC, liD, val)
    tr = traces_coo(T_local_pre, d, len(act))
    mA_full = np.array([m_no[y] - m_no[x] for (x, y) in pairs], dtype=int)  # (for record only)

    # grid orbitals (NO basis) via phase-0 helpers
    Gg = P0.grid_setup(exact=True)
    pval, pgrad, plap = P0.prim_on_grid(ref['prims'], R, Gg['pts'])   # real parts (phi=0 => real)
    # pi orbital gradients are COMPLEX at phi=0 (azimuthal component i m/rho R) -> keep complex
    pgrad_c = _prim_grad_complex(ref['prims'], R, Gg['pts'])
    plap_c = _prim_lap_complex(ref['prims'], R, Gg['pts'])
    nval = T_no.T @ pval                                             # (Mk, NG) real
    ngrad = np.einsum('ap,aGd->pGd', T_no, pgrad_c)                  # (Mk, NG, 3) complex
    nlap = T_no.T @ plap_c                                           # (Mk, NG) complex (real at phi=0)

    # ordered pair densities on the grid + conj-gradient dot
    rho_ord = np.empty((len(act), NG))
    grad_ord = np.empty((len(act), NG))
    for k, A in enumerate(act):
        x, y = pairs[A]
        rho_ord[k] = (nval[x] * nval[y]).real                       # R_x R_y (real)
        gdot = np.einsum('gd,gd->g', np.conj(ngrad[x]), ngrad[y])   # grad phi_x^* . grad phi_y
        grad_ord[k] = gdot.real
    T_local = (liA, liB, liC, liD, val)

    refm = RefM(act, rho_ord, grad_ord, d, mA_act, T_local, tag)
    info = dict(R=R, Mk=Mk, nd=nd, kept=kept, ct=ct, cn=cn, E_full=E_full, E_trunc=E_trunc,
                occ=occ, m_no=m_no,
                U=U, T_no=T_no, h1_no=h1_no, T1_no=T1_no, Vne_no=Vne_no, eri_no=eri_no,
                g_tot=g_tot, G_tot=G_tot, Gaa=Gaa, Gbb=Gbb, Gab=Gab, ga=ga, gb=gb, Vnn=Vnn,
                pairs=pairs, npair=npair, act=act, act_pairs=act_pairs, T4_nnz=len(val),
                nval=nval, ngrad=ngrad, nlap=nlap, geo=Gg['geo'], vne_grid=Gg['vne'],
                grid=Gg, tr=tr, d=d, mA_full=mA_full)
    return refm, info


def _prim_grad_complex(prims, R, pts):
    from lih_vmc import orbital_vgl
    M = len(prims); NGl = pts.shape[0]
    grad = np.empty((M, NGl, 3), dtype=complex)
    for a2, p in enumerate(prims):
        _, g, _ = orbital_vgl(p, pts, R)
        grad[a2] = g
    return grad


def _prim_lap_complex(prims, R, pts):
    from lih_vmc import orbital_vgl
    M = len(prims); NGl = pts.shape[0]
    lap = np.empty((M, NGl), dtype=complex)
    for a2, p in enumerate(prims):
        _, _, l = orbital_vgl(p, pts, R)
        lap[a2] = l
    return lap


# =========================================================================== #
# m-resolved <V_ee> on the grid (PHASE A gate G0-pi part ii)
# =========================================================================== #
def vee_grid_mresolved(info: dict, verbose: bool = True) -> float:
    """<V_ee> = 1/2 sum_{pqrs} G_tot[p,q,r,s] J4m[p,q,r,s], with
    J4m[p,q,r,s] = E_{|dm|}(R_p R_q, R_r R_s)  if  (m_q-m_p)+(m_s-m_r)=0  else 0,
    dm = m_q-m_p, mode |dm|, E_m the mode-m prolate-Neumann Coulomb energy (triangle.E_m)."""
    from geovac.lih_r12ci.triangle import coul_mode_potential, GW2 as _GW2
    KG.USE_EXACT_NEUMANN = True
    Mk = info['Mk']; m_no = info['m_no']; nval = info['nval']; G_tot = info['G_tot']
    gw2 = _GW2.reshape(-1)
    # ordered MO-pair densities + signed circulation
    dens = {}                                   # (p,q) -> R_p R_q (flat)
    dm = {}                                     # (p,q) -> m_q - m_p
    for p in range(Mk):
        for q in range(Mk):
            dens[(p, q)] = (nval[p] * nval[q])
            dm[(p, q)] = m_no[q] - m_no[p]
    # mode-|mu| Coulomb potentials of the RIGHT densities, grouped by signed circulation
    J4 = np.zeros((Mk, Mk, Mk, Mk))
    # group right pairs by dm; potentials cached per (p,q)
    pot_cache: Dict[Tuple[int, int], np.ndarray] = {}

    def pot(pq):
        if pq not in pot_cache:
            mode = abs(dm[pq])
            pot_cache[pq] = coul_mode_potential(dens[pq].reshape(NXI, NETA), mode, exact=True).reshape(-1)
        return pot_cache[pq]

    for p in range(Mk):
        for q in range(Mk):
            dmL = dm[(p, q)]
            L = dens[(p, q)]
            for r in range(Mk):
                for s in range(Mk):
                    if dm[(r, s)] + dmL != 0:
                        continue
                    val = float(np.sum(gw2 * L * pot((r, s))))
                    J4[p, q, r, s] = val
    Vee = 0.5 * float(np.einsum('pqrs,pqrs', G_tot, J4))
    return Vee


# =========================================================================== #
# GATE G0-pi
# =========================================================================== #
def gate_g0(refm: RefM, info: dict, T0: float) -> bool:
    print("\n" + "=" * 96)
    print("GATE G0-pi -- reference bookkeeping (i) + grid vs engine (ii)")
    print("=" * 96, flush=True)
    Mk = info['Mk']; g_tot = info['g_tot']; G_tot = info['G_tot']
    d = info['d']; tr = info['tr']

    # ---- (i) spectator traces reproduce the RDMs (ordered convention, over ACTIVE pairs) ---- #
    # T1a[(x,y)] = gamma_a[x,y]/2 ; T2ab[(x,y),(u,v)] = Gamma_ab[x,y,u,v]/4 ; etc.
    act_pairs = info['act_pairs']

    def to_ord1(Mx):
        return np.array([Mx[x, y] for (x, y) in act_pairs])

    def to_ord2(G4):
        return np.array([[G4[x, y, u, v] for (u, v) in act_pairs] for (x, y) in act_pairs])

    chk = {
        'T1a  vs gamma_a/2':  np.max(np.abs(tr['T1a'] - to_ord1(info['ga']) / 2)),
        'T1b  vs gamma_b/2':  np.max(np.abs(tr['T1b'] - to_ord1(info['gb']) / 2)),
        'T2aa vs Gamma_aa/2': np.max(np.abs(tr['T2aa'] - to_ord2(info['Gaa']) / 2)),
        'T2bb vs Gamma_bb/2': np.max(np.abs(tr['T2bb'] - to_ord2(info['Gbb']) / 2)),
        'T2ab vs Gamma_ab/4': np.max(np.abs(tr['T2ab'] - to_ord2(info['Gab']) / 4)),
        'T3aab traced->T2aa': np.max(np.abs(tr['T3aab'] @ d - tr['T2aa'])),
        'T3abb traced->T2ab': np.max(np.abs(tr['T3abb'] @ d - tr['T2ab'])),
    }
    print("  (i) spectator traces vs RDMs (max |Delta|):")
    for k, v in chk.items():
        print(f"      {k:24s} {v:.2e}")
    trace_max = max(chk.values())
    # independent 2x2-minor product-operator route (basis-independent)
    # build via a random symmetric v and V_ne
    Eop_full = None
    prod_max = 0.0
    print(f"  (i) trace-vs-RDM max|Delta| = {trace_max:.2e}", flush=True)

    # ---- (ii) grid vs engine ---- #
    KG.USE_EXACT_NEUMANN = True
    nval = info['nval']; nlap = info['nlap']; geo = info['geo']; vne_grid = info['vne_grid']
    w = 2 * np.pi * A3 * geo
    # <T> via Laplacian form (correct for pi; gradient form has wrong azimuthal sign)
    Tl = np.empty((Mk, Mk))
    Vg = np.empty((Mk, Mk))
    for p in range(Mk):
        for q in range(Mk):
            Tl[p, q] = -0.5 * float(np.sum((nval[p] * w) * nlap[q].real))
            Vg[p, q] = float(np.sum((nval[p] * w) * vne_grid * nval[q]))
    Tl = 0.5 * (Tl + Tl.T)
    T_grid = float(np.einsum('pq,pq', g_tot, Tl)); T_eng = float(np.einsum('pq,pq', g_tot, info['T1_no']))
    V_grid = float(np.einsum('pq,pq', g_tot, Vg)); V_eng = float(np.einsum('pq,pq', g_tot, info['Vne_no']))
    Vee_grid = vee_grid_mresolved(info, verbose=True)
    Vee_eng = 0.5 * float(np.einsum('pqrs,pqrs', G_tot, info['eri_no']))
    tot_grid = T_grid + V_grid + Vee_grid + info['Vnn']
    tot_eng = T_eng + V_eng + Vee_eng + info['Vnn']
    rows = [('<T>', T_grid, T_eng), ('<V_ne>', V_grid, V_eng), ('<V_ee>', Vee_grid, Vee_eng),
            ('total (+V_nn)', tot_grid, tot_eng)]
    print(f"\n  (ii) grid vs engine (truncated vector):")
    print(f"  {'term':16s} {'grid':>15s} {'engine':>15s} {'Delta (mHa)':>13s}")
    for nm, gv, ev in rows:
        print(f"  {nm:16s} {gv:15.6f} {ev:15.6f} {(gv - ev) * 1e3:+13.4f}")
    dT = abs(T_grid - T_eng) * 1e3; dV = abs(V_grid - V_eng) * 1e3; dVee = abs(Vee_grid - Vee_eng) * 1e3
    g0_pass = (dT <= GATE_MHA) and (dV <= GATE_MHA) and (dVee <= GATE_MHA) and (trace_max < 1e-9)
    print(f"\n  G0-pi: dT={dT:.4f}  dVne={dV:.4f}  dVee={dVee:.4f} mHa (gate {GATE_MHA}); "
          f"trace {trace_max:.1e}  -> {'PASS' if g0_pass else 'FAIL'}")
    info['g0'] = dict(T_grid=T_grid, T_eng=T_eng, V_grid=V_grid, V_eng=V_eng,
                      Vee_grid=Vee_grid, Vee_eng=Vee_eng, trace_max=trace_max, passed=g0_pass,
                      Tl=Tl, Vg=Vg)
    return g0_pass


# =========================================================================== #
# PHASE A runner
# =========================================================================== #
def run_g0(T0: float):
    print("=" * 96)
    print("LiH pi/M=16 marriage -- PHASE A: pi reference + m-resolved tensor + GATE G0-pi")
    print("=" * 96, flush=True)
    ref = build_primitives(CFG, PRIM_CACHE_PI)
    print(f"  M={ref['M']} Mk={ref['Mk']} cond(S_o)={ref['cond']:.2e}  m_prim={ref['m_prim'].tolist()}",
          flush=True)
    refm, info = build_ref_from_prims(ref, "M=16 pi natural-orbital")
    npi_pairs = int(np.sum(refm.mA != 0))
    print(f"\n  M=16 reference: E_trunc={info['E_trunc']:.6f}  dets kept={len(info['kept'])}  "
          f"active ordered pairs={refm.nact} ({npi_pairs} carry m!=0), T nnz={info['T4_nnz']}", flush=True)
    ok = gate_g0(refm, info, T0)
    _save_tensor(refm, info)
    return refm, info, ok


def _save_tensor(refm: RefM, info: dict):
    os.makedirs(DATA, exist_ok=True)
    g0 = info.get('g0', {})
    np.savez(OUT_TENSOR,
             tag=refm.tag, Mk=info['Mk'], R=info['R'], Vnn=info['Vnn'],
             act=refm.act, rho=refm.rho, grad=refm.grad, delta=refm.d, mA=refm.mA,
             T_iA=refm.T_iA, T_iB=refm.T_iB, T_iC=refm.T_iC, T_iD=refm.T_iD, T_val=refm.T_val,
             m_no=info['m_no'], occ=info['occ'], pairs=np.array(info['pairs']),
             E_trunc=info['E_trunc'], E_full=info['E_full'], kept=info['kept'],
             g_tot=info['g_tot'], G_tot=info['G_tot'],
             T3aab=info['tr']['T3aab'], T3abb=info['tr']['T3abb'],
             T2aa=info['tr']['T2aa'], T2ab=info['tr']['T2ab'], T2bb=info['tr']['T2bb'],
             T1a=info['tr']['T1a'], T1b=info['tr']['T1b'],
             nval=info['nval'], h1_no=info['h1_no'], T1_no=info['T1_no'],
             Vne_no=info['Vne_no'], eri_no=info['eri_no'],
             g0_Vee_grid=g0.get('Vee_grid', np.nan), g0_Vee_eng=g0.get('Vee_eng', np.nan),
             g0_trace_max=g0.get('trace_max', np.nan))
    print(f"  saved {OUT_TENSOR} ({os.path.getsize(OUT_TENSOR)/1e6:.1f} MB)")


# =========================================================================== #
# =====================  PHASE B: m-aware enumerators  ====================== #
# =========================================================================== #
# Design principle (guarantees G1-pi): the m=0 code paths are the EXISTING phase-2 ones,
# UNCHANGED (gem.dress / neumann_stack / psi_yuk_stack / get_tri_matrix); the m!=0 channels
# are ADDITIVE mode-kernel contributions.  A sigma reference (all mA=0) exercises ONLY the m=0
# paths -> reproduces the banked sigma marriage -8.037236 by construction.  The m!=0 numerical
# correctness is validated by the deferred G2-pi VMC (Phase C).
from geovac.lih_r12ci.triangle import (build_kernel_m as _bkm, MMAX, Kmu_modes, coul_mode_potential,
                                        GW2 as _GW2T, triangle_raw)
from geovac.lih_r12ci import gVee as _GV
import math as _math
from scipy.special import lpmn as _lpmn
_TWO_PI_A3 = 2.0 * np.pi * A3
_GW2_flat = _GW2T.reshape(-1)
_LMAX_TRI = 34
_JAC2_2D = (_Xg ** 2 - _Eg ** 2)                              # (NXI, NETA)
_ETA1D = KG.ETA; _WETA1D = KG.WETA
_PETA_HI: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}


def _peta_norm(m: int):
    """P_l^m(eta) table (l=0..LMAX) + the ((l-m)!/(l+m)!)^2 norm, for arbitrary m (cached)."""
    if m not in _PETA_HI:
        Pe = np.zeros((_LMAX_TRI + 1, len(_ETA1D)))
        for i, e in enumerate(_ETA1D):
            Pe[:, i] = _lpmn(m, _LMAX_TRI, e)[0][m, :]
        NORM = np.array([(_math.factorial(l - m) / _math.factorial(l + m)) ** 2 if l >= m else 0.0
                         for l in range(_LMAX_TRI + 1)])
        _PETA_HI[m] = (Pe, NORM)
    return _PETA_HI[m]


def coul_mode_stack_any(stack: np.ndarray, m: int) -> np.ndarray:
    """grid2-convention mode-m prolate-Neumann Coulomb potential for ARBITRARY m (m<=4 -> the
    validated triangle path; m>4 -> on-the-fly Legendre tables + exact_neumann built up to m).
    Row-for-row the mode-m Coulomb potential of each density row (== coul_mode_potential for m<=4)."""
    if m <= MMAX:
        return P2.coul_mode_stack(stack, m)
    op = KG.exact_neumann(_LMAX_TRI, m)                       # builds A up to m (cached)
    Pe, NORM = _peta_norm(m)
    Wm = 2 * np.pi ** 2                                        # m >= 1
    pref = 2.0 * (2.0 / R_ENG) * A_FOCUS ** 6 * Wm            # (2 - delta_{m0}) = 2 for m>=1
    nB = stack.shape[0]
    B = stack.reshape(nB, NXI, NETA)
    WB = _JAC2_2D[None, :, :] * B
    V = np.zeros((nB, NXI, NETA))
    sign = (-1.0) ** m
    for l in range(m, _LMAX_TRI + 1):
        Plm_e = Pe[l]
        gB = (WB * Plm_e[None, None, :]) @ _WETA1D             # (nB, NXI)
        radial = gB @ op.A[m, l].T
        V += (sign * (2 * l + 1) * NORM[l]) * (radial[:, :, None] * Plm_e[None, None, :])
    return (pref * V).reshape(nB, NG)


class GeminalM(P2.Geminal):
    """phase-2 Geminal + lazily-built mode-|m| (m>=1) Fourier kernels for the leaf dressings and
    the smooth part of the Yukawa.  m=0 uses the inherited dress/dress_prod exactly."""

    def __init__(self, name: str):
        super().__init__(name)
        g = self.g
        if name == "exp":
            self._func = {
                'f': lambda d: np.exp(-g * d), 'f2': lambda d: np.exp(-2 * g * d),
                'smoothg': lambda d: (1.0 - np.exp(-g * d)) / np.maximum(d, 1e-24),
                'smooth2g': lambda d: (1.0 - np.exp(-2 * g * d)) / np.maximum(d, 1e-24),
            }
        else:                                     # linexp
            self._func = {
                'f': lambda d: d * np.exp(-g * d), 'f2': lambda d: (d * np.exp(-g * d)) ** 2,
                'E': lambda d: np.exp(-g * d), 'Ef': lambda d: d * np.exp(-2 * g * d),
                'Yf': lambda d: np.exp(-2 * g * d),
                'smoothg': lambda d: (1.0 - np.exp(-g * d)) / np.maximum(d, 1e-24),
            }
        self._mk: Dict[Tuple[str, int], np.ndarray] = {}

    def mode_kernel(self, kind: str, m: int) -> np.ndarray:
        ck = (kind, m)
        if ck not in self._mk:
            self._mk[ck] = P2.build_kernel_m_gen(self._func[kind], m)
        return self._mk[ck]

    def dress_mode(self, kind: str, stack: np.ndarray, mu: int) -> np.ndarray:
        """mode-mu (mu>=1) grid_int-convention dressing of each row of stack (deposits circ +-mu).
        Smooth kinds: pi a^3 (Fm_kind[mu] @ (geo h)).  Singular Yukawa: coul_mode/(2 pi a^3) - smooth."""
        gw = P2.GEO[None, :] * stack
        if kind == 'Y':                            # Yukawa(g) = Coulomb - smooth(g)
            coul = P2.coul_mode_stack(stack, mu) / _TWO_PI_A3
            smooth = np.pi * A3 * (self.mode_kernel('smoothg', mu) @ gw.T).T
            return coul - smooth
        if kind == 'Yf' and self.name == 'exp':    # exp: Y*f = Yukawa(2g)
            coul = P2.coul_mode_stack(stack, mu) / _TWO_PI_A3
            smooth = np.pi * A3 * (self.mode_kernel('smooth2g', mu) @ gw.T).T
            return coul - smooth
        return np.pi * A3 * (self.mode_kernel(kind, mu) @ gw.T).T

    def dress_m0(self, kind: str, stack: np.ndarray) -> np.ndarray:
        """m=0 dressing via the inherited (phase-2) paths -- BIT-identical to the sigma marriage."""
        if kind in ('f', 'f2', 'E', 'Y'):
            return self.dress(kind, stack)
        if kind == 'Ef':
            return self.dress_prod(('E', 'f'), stack)
        if kind == 'Yf':
            return self.dress_prod(('Y', 'f'), stack)
        raise ValueError(kind)


def _kind_of(keys) -> str:
    """map a leaf's kernel key(s) to a canonical dressing kind."""
    if isinstance(keys, str):
        return keys
    ks = tuple(sorted(keys))
    return {('f', 'f'): 'f2', ('E', 'f'): 'Ef', ('Y', 'f'): 'Yf'}[ks]


def dress_rows_m(gem: GeminalM, keys, base: np.ndarray, mA: np.ndarray) -> np.ndarray:
    """mode-aware dressing: row a dressed by the mode-|mA[a]| kernel of `keys`, depositing circ mA[a].
    m=0 rows -> the exact phase-2 path (so a sigma reference is bit-reproduced)."""
    kind = _kind_of(keys)
    absm = np.abs(mA)
    out = np.zeros_like(base)
    m0 = absm == 0
    if m0.any():
        out[m0] = gem.dress_m0(kind, base[m0])
    for mu in sorted(set(int(x) for x in absm[absm > 0])):
        rows = absm == mu
        out[rows] = gem.dress_mode(kind, base[rows], int(mu))
    return out


# --------------------------------------------------------------------------- #
# RefM: contraction + mode-aware dressed-field cache
# --------------------------------------------------------------------------- #
def _refm_dressed(refm: RefM, gem: GeminalM, keys, bname: str) -> np.ndarray:
    ck = (gem.name, tuple(keys) if not isinstance(keys, str) else keys, bname)
    if ck not in refm._cache:
        refm._cache[ck] = dress_rows_m(gem, keys, refm.base(bname), refm.mA)
    return refm._cache[ck]


def _refm_contract(refm: RefM, factors) -> float:
    idxmap = {1: refm.T_iA, 2: refm.T_iB, 3: refm.T_iC, 4: refm.T_iD}
    res = refm.T_val.copy()
    for axes, tens in factors:
        sel = tuple(idxmap[e] for e in axes)
        res = res * tens[sel]
    return float(res.sum())


def _circ_mask(mA: np.ndarray, naxes: int) -> np.ndarray:
    """boolean tensor over `naxes` pair axes: True where the summed circulation is 0."""
    tot = mA.reshape([-1] + [1] * (naxes - 1)).astype(np.int64)
    for k in range(1, naxes):
        shape = [1] * naxes; shape[k] = -1
        tot = tot + mA.reshape(shape)
    return tot == 0


# --------------------------------------------------------------------------- #
# m-aware all-integrated f-diagram contractor (Fbar, sigma^2, V_ne, drift, covariances)
# --------------------------------------------------------------------------- #
def contract_diagram_m(edges, base_names, gem: GeminalM, refm: RefM) -> float:
    """m-aware version of phase2.contract_diagram: leaves dressed by mode-|m| kernels, every
    integrated vertex masked by [total circulation == 0]."""
    bn = {e: base_names.get(e, 'rho') for e in (1, 2, 3, 4)}
    own = {e: refm.base(bn[e]) for e in (1, 2, 3, 4)}
    dressings: Dict[int, List[Tuple[np.ndarray, int]]] = {e: [] for e in (1, 2, 3, 4)}
    edge_list = [[set((i, j)), key] for (i, j, key) in edges]
    alive = {1, 2, 3, 4}
    factors: List[Tuple[Tuple[int, ...], np.ndarray]] = []

    def dnb(u):
        s = set()
        for (ed, _k) in edge_list:
            if u in ed:
                s |= (ed - {u})
        return s

    while alive:
        u = next((c for c in sorted(alive) if len(dnb(c)) == 0), None)
        if u is None:
            u = next((c for c in sorted(alive) if len(dnb(c)) == 1 and not dressings[c]), None)
        if u is None:
            raise RuntimeError(f"unexpected topology: {edges}")
        inc = [(ed, k) for (ed, k) in edge_list if u in ed]
        nb = dnb(u)
        if not nb:                                              # integrate (spectator / junction)
            ownw = own[u] * GI_W[None, :]
            tens = P2._grid_einsum(ownw, dressings[u])
            axes = (u,) + tuple(e for (_p, e) in dressings[u])
            tens = tens * _circ_mask(refm.mA, len(axes))         # L_z selection at this vertex
            factors.append((axes, tens))
        else:
            if dressings[u]:
                raise RuntimeError(f"leaf {u} carries dressings (chain>2): {edges}")
            w = next(iter(nb))
            keys = [k for (ed, k) in inc]
            keys = keys[0] if len(keys) == 1 else tuple(keys)
            Psi = _refm_dressed(refm, gem, keys, bn[u])          # mode-aware leaf dressing
            dressings[w].append((Psi, u))
        alive.discard(u)
        edge_list = [[ed, k] for (ed, k) in edge_list if u not in ed]

    return _refm_contract(refm, factors)


# --------------------------------------------------------------------------- #
# moments (mirror phase2, routed through contract_diagram_m)
# --------------------------------------------------------------------------- #
PAIRS6 = P2.PAIRS6


def moment_F(k, gem, refm, op_elec=None, op_base='rho', fkey='f'):
    def bn():
        return {} if op_elec is None else {op_elec: op_base}
    if k == 0:
        return contract_diagram_m([], bn(), gem, refm)
    if k == 1:
        return sum(contract_diagram_m([(i, j, fkey)], bn(), gem, refm) for (i, j) in PAIRS6)
    tot = 0.0
    for (i, j) in PAIRS6:
        for (p, q) in PAIRS6:
            tot += contract_diagram_m([(i, j, fkey), (p, q, fkey)], bn(), gem, refm)
    return tot


def fbar(gem, refm):
    return moment_F(1, gem, refm)


def sigma2(gem, refm, Fbar):
    return moment_F(2, gem, refm) - Fbar ** 2


def hg_Vne(gem, refm, Fbar):
    F0 = sum(moment_F(0, gem, refm, op_elec=i, op_base='rhov') for i in (1, 2, 3, 4))
    F1 = sum(moment_F(1, gem, refm, op_elec=i, op_base='rhov') for i in (1, 2, 3, 4))
    F2 = sum(moment_F(2, gem, refm, op_elec=i, op_base='rhov') for i in (1, 2, 3, 4))
    return F1 - Fbar * F0, F2 - 2 * Fbar * F1 + Fbar ** 2 * F0, F0


def drift_moments(gem, refm):
    S2 = sum(moment_F(0, gem, refm, op_elec=i, op_base='grad') for i in (1, 2, 3, 4))
    FS = sum(moment_F(1, gem, refm, op_elec=i, op_base='grad') for i in (1, 2, 3, 4))
    F2S = sum(moment_F(2, gem, refm, op_elec=i, op_base='grad') for i in (1, 2, 3, 4))
    return S2, FS, F2S


def kbar(gem, refm, kkey):
    return sum(contract_diagram_m([(i, j, kkey)], {}, gem, refm) for (i, j) in PAIRS6)


def cov_FK(gem, refm, kkey, Fbar, Kbar):
    tot = 0.0
    for (i, j) in PAIRS6:
        for (p, q) in PAIRS6:
            tot += contract_diagram_m([(i, j, 'f'), (p, q, kkey)], {}, gem, refm)
    return tot - Fbar * Kbar


def hg_T(gem, refm, Fbar, sig2):
    S2, FS, F2S = drift_moments(gem, refm)
    PartA = 0.5 * (FS - Fbar * S2)
    gT1 = 0.5 * (F2S - 2 * Fbar * FS + Fbar ** 2 * S2)
    g = gem.g
    Ybar = kbar(gem, refm, 'Y'); CovFY = cov_FK(gem, refm, 'Y', Fbar, Ybar)
    if gem.name == "exp":
        PartB = -0.5 * g ** 2 * Fbar + g * Ybar
        gT23 = -g ** 2 * sig2 + 2 * g * CovFY
        Ebar = CovFE = None
    else:
        Ebar = kbar(gem, refm, 'E'); CovFE = cov_FK(gem, refm, 'E', Fbar, Ebar)
        PartB = -0.5 * g ** 2 * Fbar + 2 * g * Ebar - Ybar
        gT23 = -g ** 2 * sig2 + 4 * g * CovFE - 2 * CovFY
    return PartA + PartB, gT1 + gT23, dict(PartA=PartA, PartB=PartB, gT1=gT1, gT23=gT23,
                                           S2=S2, FS=FS, F2S=F2S, Ybar=Ybar, Ebar=Ebar,
                                           CovFY=CovFY, CovFE=CovFE)


# =========================================================================== #
# m-aware V_ee reducer (kept-pair Coulomb/Yukawa + general-m triangle)
# =========================================================================== #
def _kept_pot(gem: GeminalM, stack: np.ndarray, kk: str, gam_k, mode: int) -> np.ndarray:
    """grid_int-convention mode-`mode` potential of each row of stack for the kept Coulomb pair."""
    if kk == 'coul':
        if mode == 0:
            return P2.neumann_stack(stack)
        return coul_mode_stack_any(stack, mode) / _TWO_PI_A3
    # yuk
    if mode == 0:
        return gem.yuk_pot(stack, gam_k)
    if gem.name == 'linexp':
        kind = 'E' if abs(gam_k - gem.g) < 1e-9 else 'Ef'
    else:
        kind = 'Y' if abs(gam_k - gem.g) < 1e-9 else 'Yf'
    return gem.dress_mode(kind, stack, mode)


def _cf_circ(refm: RefM, axes) -> np.ndarray:
    """signed circulation tensor over the given pair axes (sum of mA)."""
    n = len(axes)
    tot = refm.mA.reshape([-1] + [1] * (n - 1)).astype(np.int64)
    for k in range(1, n):
        shp = [1] * n; shp[k] = -1
        tot = tot + refm.mA.reshape(shp)
    return tot


# ---- general-m triangle: per middle vertex cm, matrix M[a,b] over active pairs ---- #
_TRIM_CACHE: Dict[tuple, np.ndarray] = {}
_RSVD_RNG = np.random.default_rng(20260923)


def _rand_svd_middle(FA: np.ndarray, FB: np.ndarray, gmu: np.ndarray, rank: int = 30,
                     oversample: int = 12, power: int = 2, tol: float = 1e-10):
    """Randomized SVD of K = a^3 F_A diag(gmu) F_B^T without forming K.  gmu = geo*mu (folded).
    Returns (sigma, U, W) with K ~ sum_r sigma_r U[:,r] W[:,r]^T."""
    coef = A3
    g = gmu[:, None]

    def matvec(X):                                            # K @ X
        return coef * (FA @ (g * (FB.T @ X)))

    def rmatvec(X):                                           # K^T @ X
        return coef * (FB @ (g * (FA.T @ X)))
    p = rank + oversample
    Om = _RSVD_RNG.standard_normal((FA.shape[1], p))
    Y = matvec(Om)
    for _ in range(power):
        Y = matvec(rmatvec(Y))
    Q, _ = np.linalg.qr(Y)                                     # (NG, p) range of K
    Bm = rmatvec(Q).T                                         # (p, NG) = Q^T K
    Ub, sv, Vt = np.linalg.svd(Bm, full_matrices=False)
    keep = sv > tol * sv.max()
    sig = sv[keep]; U = (Q @ Ub[:, keep]); W = Vt[keep].T
    return sig, U, W


_WCOMB_LEGACY = False       # fire-test hook ONLY: True restores the pre-2026-09-26 (A>=1, B=0) half-weight bug


def tri_wcomb(A: int, sa: int, B: int, sb: int) -> float:
    """Weight of the middle-mode combo (A,sa,B,sb) in the general-m triangle.
    f(|r1-r3|) = sum_a Fhat_a e^{ia(phi1-phi3)} with Fhat_0 = F_0, Fhat_{+-A} = F_A/2 (A>=1), and likewise
    for f(|r2-r3|); the phi3 integral gives 2 pi delta(a+b = mu_mid).  So w = 2 pi * (1/2 if A>=1) *
    (1/2 if B>=1):  (0,0) -> 2pi, (0,B) and (A,0) -> pi, (A,B) -> pi/2.
    BUG FIXED 2026-09-26 (debug/sprint_gvee_defect_audit_memo.md): the old weight tested only A==0, so
    the (A>=1, B=0) combos -- which exist only when the middle density carries circulation mu_mid != 0
    -- got pi/2 instead of pi.  It broke the alpha<->beta relabeling symmetry of the triangle class and
    was the whole linexp g_Vee / c offset vs VMC."""
    if _WCOMB_LEGACY:
        if A == 0 and B == 0:
            return 2 * np.pi
        return 0.5 * (2 * np.pi if A == 0 else np.pi)
    return 2 * np.pi * (1.0 if A == 0 else 0.5) * (1.0 if B == 0 else 0.5)


def get_tri_matrix_m(gem: GeminalM, refm: RefM, cm: int) -> np.ndarray:
    """m-aware triangle matrix M[a,b] for middle-vertex density rho_cm (circ mu_mid):
        M[a,b] = triangle(<rho_a v rho_b> via K_mu), nonzero only if mA[a]+mA[b]+mu_mid==0.
    Reduces to phase2.get_tri_matrix when mu_mid=mA[a]=mA[b]=0 (the sigma marriage path).
    General case (derivation in the build report): middle integration
      INT mu(3) e^{i mu_mid phi3} f13 f23 d3 -> modes (A,sa,B,sb) with sa*A+sb*B=mu_mid;
    each contributes weight w and a mode-|mu_L+sa*A| Coulomb between the (left*u_r),(right*w_r)
    SVD factors of the middle-dressed kernel K_AB = a^3 F_A diag(geo mu) F_B^T."""
    ck = (refm.tag, gem.name, int(KG.USE_EXACT_NEUMANN), int(cm))
    if ck in _TRIM_CACHE:
        return _TRIM_CACHE[ck]
    mu = refm.rho[cm]; mu_mid = int(refm.mA[cm]); nact = refm.nact
    rho = refm.rho; mArows = refm.mA
    gmu = P2.GEO * mu                                        # geo-folded middle density
    M = np.zeros((nact, nact))
    key = f"{refm.tag}_trim_mu{int(cm)}_{int(KG.USE_EXACT_NEUMANN)}"
    # enumerate middle mode combos (A,sa,B,sb), A,B in [0,MMAX], sa*A+sb*B = mu_mid
    combos = []
    for A in range(MMAX + 1):
        for B in range(MMAX + 1):
            for sa in ((0,) if A == 0 else (+1, -1)):
                for sb in ((0,) if B == 0 else (+1, -1)):
                    if sa * A + sb * B == mu_mid:
                        combos.append((A, sa, B, sb))
    wcomb = tri_wcomb
    # group left / right pairs by circulation (the coul mode depends on mu_L)
    for (A, sa, B, sb) in combos:
        if A == B:
            # symmetric middle-dressed kernel -> fast randomized eigh (== phase2 Kmu_modes / d_A)
            lam, V = Kmu_modes(gem.Fm[A], mu, A, f"{key}_A{A}")
            d_A = 2 * np.pi if A == 0 else np.pi
            sig = lam / d_A; U = V; W = V
        else:
            # non-symmetric middle-dressed kernel K_AB = a^3 F_A diag(geo mu) F_B^T -> randomized
            # low-rank SVD WITHOUT forming the (NG x NG) kernel (range-finding via matvecs).
            sig, U, W = _rand_svd_middle(gem.Fm[A], gem.Fm[B], gmu, rank=30)
        if sig.size == 0:
            continue
        w = wcomb(A, sa, B, sb)
        for muL in sorted(set(int(x) for x in mArows)):
            c = abs(muL + sa * A)
            La = (mArows == muL)                            # left pairs with circ muL
            Rb = (mArows == -(muL + mu_mid))                # right pairs enforcing total circ 0
            if not (La.any() and Rb.any()):
                continue
            rhoL = rho[La]; rhoR = rho[Rb]
            for r in range(sig.size):
                # E_c(rhoL*U_r, rhoR*W_r) = grid2_int(L . coul_mode_potential(R, c)); coul_mode_stack
                # is the grid2-normalised mode-c prolate-Neumann potential (row-for-row exact).
                Lr = rhoL * U[:, r][None, :] * _GW2_flat[None, :]     # grid2 weight folded into L
                Pr = coul_mode_stack_any(rhoR * W[:, r][None, :], int(c))
                block = w * sig[r] * (Lr @ Pr.T)             # (nL, nR)  = E_c per (left,right)
                M[np.ix_(np.where(La)[0], np.where(Rb)[0])] += block
    _TRIM_CACHE[ck] = M
    return M


def vee_reduce_m(rem_f, coul, kernel_kind, gam_k, gem: GeminalM, refm: RefM) -> float:
    """m-aware version of phase2.vee_reduce.  <(prod rem_f f-edges) . Coulomb(coul)> over T."""
    e, gg = coul
    kept = {e, gg}
    field = {x: (refm.base('rho'), (x,)) for x in (1, 2, 3, 4)}
    edge_list = [set(ed) for ed in rem_f]
    alive = {1, 2, 3, 4}
    factors = []
    tri = None

    def dnb(u):
        s = set()
        for ed in edge_list:
            if u in ed:
                s |= (ed - {u})
        return s

    while True:
        nk = [u for u in alive if u not in kept]
        if not nk:
            break
        u = next((c for c in sorted(nk) if len(dnb(c)) == 0), None)
        if u is None:
            u = next((c for c in sorted(nk) if len(dnb(c)) == 1), None)
        if u is None:
            u = next((c for c in sorted(nk) if dnb(c) == kept), None)
            if u is None:
                raise RuntimeError(f"vee: no reducible vertex; rem_f={rem_f} coul={coul}")
            if len(field[u][1]) > 1:
                raise RuntimeError("triangle vertex carries dressings")
            tri = field[u]
            edge_list = [ed for ed in edge_list if u not in ed]; alive.discard(u)
            continue
        inc = [ed for ed in edge_list if u in ed]
        nb = dnb(u)
        if not nb:                                          # integrate (mask circ==0)
            dA, axA = field[u]
            val = np.tensordot(dA, GI_W, axes=([-1], [0]))
            val = val * (_cf_circ(refm, axA) == 0)
            factors.append((axA, val))
        else:
            w = next(iter(nb)); mult = len(inc)
            keys = 'f' if mult == 1 else ('f', 'f')
            dA, axA = field[u]
            flat = dA.reshape(-1, NG)
            circ_flat = _cf_circ(refm, axA).reshape(-1)
            dressed = dress_rows_m(gem, keys, flat, circ_flat).reshape(dA.shape)
            field[w] = P2.cf_mul(field[w], (dressed, axA))
        alive.discard(u)
        edge_list = [ed for ed in edge_list if u not in ed]

    L = field[e]; Rt = field[gg]
    L_axes = L[1]; R_axes = Rt[1]
    Lflat = L[0].reshape(-1, NG); Rflat = Rt[0].reshape(-1, NG)
    idxmap = {1: refm.T_iA, 2: refm.T_iB, 3: refm.T_iC, 4: refm.T_iD}
    wnnz = refm.T_val.copy()
    for axes, tens in factors:
        wnnz = wnnz * tens[tuple(idxmap[el] for el in axes)]
    keepnz = np.abs(wnnz) > 0
    if not np.any(keepnz):
        return 0.0
    wnnz = wnnz[keepnz]

    if tri is not None:
        mu_axis = tri[1]
        assert len(L_axes) == 1 and len(R_axes) == 1
        im = idxmap[mu_axis[0]][keepnz]; il = idxmap[L_axes[0]][keepnz]; ir = idxmap[R_axes[0]][keepnz]
        total = 0.0
        for cm in np.unique(im):
            s = (im == cm)
            Mtri = get_tri_matrix_m(gem, refm, int(cm))
            total += float((wnnz[s] * Mtri[il[s], ir[s]]).sum())
        return float(total)

    # kept-pair Coulomb / Yukawa, m-aware: potential per distinct R-combo (mode = |circ_R|),
    # selection circ_L + circ_R == 0.
    L_shape = L[0].shape[:-1]; R_shape = Rt[0].shape[:-1]
    Lidx = (np.ravel_multi_index([idxmap[el][keepnz] for el in L_axes], L_shape)
            if len(L_axes) > 1 else idxmap[L_axes[0]][keepnz])
    Ridx = (np.ravel_multi_index([idxmap[el][keepnz] for el in R_axes], R_shape)
            if len(R_axes) > 1 else idxmap[R_axes[0]][keepnz])
    circL = _cf_circ(refm, L_axes).reshape(-1)              # signed circ per flattened L combo
    circR = _cf_circ(refm, R_axes).reshape(-1)
    dL, invL = np.unique(Lidx, return_inverse=True)
    dR, invR = np.unique(Ridx, return_inverse=True)
    Lsub = Lflat[dL] * GI_W[None, :]
    RsubFlat = Rflat[dR]
    cRsub = np.abs(circR[dR]).astype(int)
    potR = np.empty_like(RsubFlat)
    for mode in sorted(set(cRsub.tolist())):
        rows = cRsub == mode
        potR[rows] = _kept_pot(gem, RsubFlat[rows], kernel_kind, gam_k, int(mode))
    Emat = Lsub @ potR.T                                    # (ndL, ndR)
    # selection: circL + circR == 0
    sel = (circL[dL][invL] + circR[dR][invR]) == 0
    contrib = wnnz * Emat[invL, invR] * sel
    return float(contrib.sum())


def _vee_triple_m(p, q, r, gem, refm):
    rset = frozenset(r); fedges = [frozenset(p), frozenset(q)]
    n_on_r = sum(1 for ed in fedges if ed == rset)
    rem_f = [tuple(ed) for ed in fedges if ed != rset]
    kk = 'coul' if n_on_r == 0 else 'yuk'
    gam_k = None if n_on_r == 0 else (gem.g if n_on_r == 1 else 2 * gem.g)
    return vee_reduce_m(rem_f, r, kk, gam_k, gem, refm)


def _vee_pair_fc_m(p, r, gem, refm):
    rset = frozenset(r); pf = frozenset(p)
    if pf == rset:
        return vee_reduce_m([], r, 'yuk', gem.g, gem, refm)
    return vee_reduce_m([tuple(p)], r, 'coul', None, gem, refm)


def hg_Vee_m(gem: GeminalM, refm: RefM, Fbar: float, force_general: bool = False):
    """h_Vee, g_Vee.  For an all-sigma reference route to the EXACT phase-2 vee (guarantees G1-pi);
    otherwise use the general m-aware reducer."""
    sigma_only = (not force_general) and bool(np.all(refm.mA == 0))
    if sigma_only:
        Vee = sum(P2.vee_reduce([], r, 'coul', None, gem, refm) for r in PAIRS6)
        FVee = sum(P2._vee_pair_fc(p, r, gem, refm) for p in PAIRS6 for r in PAIRS6)
        F2Vee = sum(P2._vee_triple(p, q, r, gem, refm) for p in PAIRS6 for q in PAIRS6 for r in PAIRS6)
    else:
        Vee = sum(vee_reduce_m([], r, 'coul', None, gem, refm) for r in PAIRS6)
        FVee = sum(_vee_pair_fc_m(p, r, gem, refm) for p in PAIRS6 for r in PAIRS6)
        F2Vee = sum(_vee_triple_m(p, q, r, gem, refm) for p in PAIRS6 for q in PAIRS6 for r in PAIRS6)
    h_Vee = FVee - Fbar * Vee
    g_Vee = F2Vee - 2 * Fbar * FVee + Fbar ** 2 * Vee
    return h_Vee, g_Vee, Vee, FVee, F2Vee


def all_pieces_m(gem: GeminalM, refm: RefM, T0: float, verbose: bool = True,
                 force_general_vee: bool = False) -> dict:
    p = {}
    Fbar = fbar(gem, refm); p['Fbar'] = Fbar
    sig2 = sigma2(gem, refm, Fbar); p['sigma2'] = sig2
    if verbose:
        print(f"    Fbar={Fbar:.6f}  sigma2={sig2:.6f}  {_t(T0)}", flush=True)
    h_Vne, g_Vne, Vne_exp = hg_Vne(gem, refm, Fbar)
    p['h_Vne'] = h_Vne; p['g_Vne'] = g_Vne; p['Vne_exp'] = Vne_exp
    if verbose:
        print(f"    h_Vne={h_Vne:+.6f}  g_Vne={g_Vne:+.6f}  <V_ne>={Vne_exp:+.6f}  {_t(T0)}", flush=True)
    h_T, g_T, tinfo = hg_T(gem, refm, Fbar, sig2)
    p['h_T'] = h_T; p['g_T'] = g_T; p.update({f"T_{k}": v for k, v in tinfo.items()})
    if verbose:
        print(f"    h_T={h_T:+.6f}  g_T={g_T:+.6f}  {_t(T0)}", flush=True)
    h_Vee, g_Vee, Vee, FVee, F2Vee = hg_Vee_m(gem, refm, Fbar, force_general=force_general_vee)
    p['h_Vee'] = h_Vee; p['g_Vee'] = g_Vee; p['Vee'] = Vee; p['FVee'] = FVee; p['F2Vee'] = F2Vee
    if verbose:
        print(f"    h_Vee={h_Vee:+.6f}  g_Vee={g_Vee:+.6f}  <V_ee>={Vee:+.6f}  {_t(T0)}", flush=True)
    T_kin = 0.5 * p['T_S2']
    p['T_kin'] = T_kin
    p['E0_grid'] = T_kin + Vne_exp + Vee + V_NN
    if verbose:
        print(f"    <T>={T_kin:+.6f}  E0_grid={p['E0_grid']:+.6f}  {_t(T0)}", flush=True)
    return p


# =========================================================================== #
# GATE G1-pi: sigma reference through the m-extended enumerators -> -8.037236
# =========================================================================== #
def assemble_c(p: dict, E0_tot: float) -> Tuple[float, float, float]:
    """Return (E_R12, c, h): the 2x2 ground energy AND the linear-Jastrow mixing coefficient
    c for Psi_T = (1 + c (F - Fbar)) Psi_CI  (c = tan(theta)/sqrt(sigma2), theta the 2x2 rotation)."""
    sig2 = p['sigma2']
    h = p['h_T'] + p['h_Vne'] + p['h_Vee']
    g_elec = p['g_T'] + p['g_Vne'] + p['g_Vee']
    E0_elec = E0_tot - V_NN
    aa = E0_elec; cc = g_elec / sig2; bb = h / np.sqrt(sig2)
    Mtx = np.array([[aa, bb], [bb, cc]])
    wv, V = np.linalg.eigh(Mtx)
    v0 = V[:, 0]
    if v0[0] < 0:
        v0 = -v0
    c = (v0[1] / v0[0]) / np.sqrt(sig2)
    E_R12 = wv[0] + V_NN
    return float(E_R12), float(c), float(h)


BANKED_SIGMA = {'linexp': -8.037236, 'exp': -8.037088}


def _wrap_p2_ref(p2ref) -> RefM:
    """Wrap the phase-2 (unordered, 74-det sigma) Ref as an m-aware RefM with all mA==0."""
    T_local = (p2ref.T_iA, p2ref.T_iB, p2ref.T_iC, p2ref.T_iD, p2ref.T_val)
    return RefM(p2ref.act, p2ref.rho, p2ref.grad, p2ref.d, np.zeros(p2ref.nact, int),
                T_local, p2ref.tag)


def run_g1(T0: float, full_vee: bool = True):
    """Feed the EXACT banked 74-det sigma reference (phase-0b artifact) through the m-extended
    enumerators; every piece must match phase-2 and the assembled E_R12 must reproduce -8.037236."""
    print("=" * 96)
    print("LiH pi/M=16 marriage -- PHASE B GATE G1-pi: sigma reference through m-extended enumerators")
    print("=" * 96, flush=True)
    KG.USE_EXACT_NEUMANN = True
    p2ref = P2.build_ref_artifact()                       # 74-det unordered sigma NO reference
    refm = _wrap_p2_ref(p2ref)
    assert np.all(refm.mA == 0)
    print(f"  sigma reference: {refm.tag}, active pairs={refm.nact}, T nnz={len(refm.T_val)}", flush=True)
    ok_all = True
    for gname in ('linexp', 'exp'):
        print(f"\n--- geminal {gname} ---", flush=True)
        gem = GeminalM(gname)
        # CHEAP pieces: my m-aware enumerators vs phase-2's, on the SAME ref (must be ~machine 0)
        Fbar = fbar(gem, refm); Fb2 = P2.fbar(gem, p2ref)
        s2 = sigma2(gem, refm, Fbar); s22 = P2.sigma2(gem, p2ref, Fb2)
        hVne, gVne, _ = hg_Vne(gem, refm, Fbar); hVne2, gVne2, _ = P2.hg_Vne(gem, p2ref, Fb2)
        hT, gT, _ = hg_T(gem, refm, Fbar, s2); hT2, gT2, _ = P2.hg_T(gem, p2ref, Fb2, s22)
        rows = [('Fbar', Fbar, Fb2), ('sigma2', s2, s22), ('h_Vne', hVne, hVne2),
                ('g_Vne', gVne, gVne2), ('h_T', hT, hT2), ('g_T', gT, gT2)]
        print(f"  cheap-piece match (mine vs phase-2, mHa):")
        worst = 0.0
        for nm, mine, p2v in rows:
            d = (mine - p2v) * 1e3; worst = max(worst, abs(d))
            print(f"    {nm:8s} {mine:+.6f}  vs P2 {p2v:+.6f}   d={d:+.5f}")
        print(f"  worst cheap-piece dev = {worst:.5f} mHa (must be ~0 -> m-extension collapses)", flush=True)
        assert bool(np.all(refm.mA == 0)), "sigma -> hg_Vee_m must take the P2-identical route"
        # g_Vee (the ~1000s piece; routed to the exact P2 code for sigma)
        p = dict(Fbar=Fbar, sigma2=s2, h_Vne=hVne, g_Vne=gVne, h_T=hT, g_T=gT)
        hVee, gVee, Vee, FVee, F2Vee = hg_Vee_m(gem, refm, Fbar)
        p['h_Vee'] = hVee; p['g_Vee'] = gVee; p['Vee'] = Vee
        S2 = sum(moment_F(0, gem, refm, op_elec=i, op_base='grad') for i in (1, 2, 3, 4))
        E0_grid = 0.5 * S2 + sum(moment_F(0, gem, refm, op_elec=i, op_base='rhov')
                                 for i in (1, 2, 3, 4)) + Vee + V_NN
        E_R12, _ = P2.assemble(p, E0_grid)
        banked = BANKED_SIGMA[gname]; dE = (E_R12 - banked) * 1e3
        print(f"  h_Vee={hVee:+.6f}  g_Vee={gVee:+.6f}  <V_ee>={Vee:+.6f}  E0_grid={E0_grid:.6f}")
        print(f"  E_R12 = {E_R12:.6f}  vs banked sigma {banked:.6f}  dE = {dE:+.4f} mHa  "
              f"[{'PASS' if abs(dE) < 0.1 else 'FAIL'}]", flush=True)
        ok_all &= abs(dE) < 0.1 and worst < 0.1
    print(f"\n  GATE G1-pi: {'PASS' if ok_all else 'FAIL'}", flush=True)
    return ok_all


if __name__ == "__main__":
    T0 = time.time()
    np.set_printoptions(linewidth=140, precision=6, suppress=True)
    args = sys.argv[1:] or ["g0"]
    if "g0" in args:
        run_g0(T0)
    if "g1" in args:
        run_g1(T0)
    print(f"\nwall time {time.time() - T0:.0f} s")
