r"""Phase 0 of the LiH "marriage" build (2026-09-22, PI-directed; plan: debug/lih_marriage_build_plan.md).

Builds the multi-determinant sigma reference, truncates it, builds the pair-index coefficient
tensor T with its spectator traces, and runs GATE G0: grid <T>, <V_ne>, <V_ee> of the truncated
vector on the geovac/lih_r12ci/kernels.py prolate grid vs the exact ERI engine (each <= 1 mHa).
No r12 / geminal work here.  Nothing in geovac/lih_r12ci/ or debug/lih_vmc.py is modified.

Run from root:  python debug/lih_marriage_phase0.py  > debug/data/lih_marriage_phase0.log 2>&1
Writes         debug/data/lih_marriage_phase0_ref.npz  (everything Phases 1-3 need)
               debug/data/lih_marriage_phase0_prim.npz (cache of the 324 s float64 ERI build)
Phase 0b:       python debug/lih_marriage_phase0.py --exact > debug/data/lih_marriage_phase0_exact.log
               sets geovac.lih_r12ci.kernels.USE_EXACT_NEUMANN = True (exact ordered-integral prolate-
               Neumann operator, neumann_exact.py) for the <V_ee> grid path, skips the legacy scratch-
               grid scan, and writes debug/data/lih_marriage_phase0_ref_exact.npz (same keys; V_no /
               J_grid_no / J4_grid_no / g0_table on the exact path -- the artifact Phase 1 should use).

CONVENTIONS (load-bearing for Phases 1-3)
------------------------------------------
Reference:  sigma-only core-enriched CI, CFG below = lih_core2exp_probe.run(2,1,0,0,0,1.0,core2=(4.5,1.6))
  primitives: Li 1s STOs zeta = 2.6875, 4.5, 1.6 (centre A) + Li2s(0.65, A) + H1s(1.0, B) + Hm(0.70, B)
  + bond xi^j eta^l e^{-xi}, j<=2, l<=1  -> M = 12, canonical-orthogonalized Mk = 12, n_det = 66^2 = 4356.
  Same pipeline as lih_vmc.build_lih_wavefunction (float64 ERI engine prolate_float_eri +
  sparse connected-pair Lanczos FCI).  E_full must reproduce -8.02335.

Orbital bases:
  primitives a (OrbitalM, real, sigma) -> canonical MOs p via Tmap = C.T @ X (M x Mk)
  -> NATURAL ORBITALS of the FULL CI (U from the spin-summed 1-RDM, occupation-descending):
     T_no = Tmap @ U.  The FCI is re-solved in the NO basis (same 12-orbital space -> same energy);
     the |c| > 1e-3 truncation, the RDMs, the tensor T and gate G0 are all in the NO basis.
  Why: in the canonical Loewdin-like basis the vector is spread over 954 dets at 1e-3 (no
  dominant determinant); the plan's "~20-60 dets" presupposes a compact (natural) basis.

Determinants:  prolate_allelectron_fci._dets ordering -- alpha pair a outer, beta pair b inner,
  pairs = itertools.combinations(range(K), 2) (lexicographic); c_I multiplies the SORTED
  spin-orbital determinant (spin-orbital = 2*spatial + spin).  Folded coefficient matrix
      CS[a,b] = c_I * sigma_I,   sigma_I = (-1)^{inversions of [2p, 2q, 2r+1, 2s+1]}
  so Psi = sum_ab CS[a,b] |a>_alpha (x) |b>_beta in GROUPED order and the (alpha,alpha,beta,beta)
  spatial component is  psi(1,2,3,4) = sum_ab CS[a,b] dA[a](1,2) dB[b](3,4),
      dA[(p,q)](1,2) = phi_p(1) phi_q(2) - phi_q(1) phi_p(2)     (lih_vmc._minors_full),
  INT |psi|^2 = 4 sum_I c_I^2.

Pair index:  A = (p,q), p <= q;  PAIRS = [(p,q) for p in range(K) for q in range(p,K)];
  npair = K(K+1)/2 = 78.   rho_A(r) = phi_p(r) phi_q(r)  -- NO factor 2 for p<q; the symmetry
  factor lives in the coefficients: a symmetric orbital matrix M_pq FOLDS to
  m_A = M_pp (p=q) or M_pq + M_qp (p<q)  (fold = F^T vec(M), F the (K^2 x npair) 0/1 map).

Coefficient tensor (the deliverable):
  P(1,2,3,4) = |psi|^2 / 4     -- normalized to 1; electrons 1,2 alpha, 3,4 beta --
  P = sum_{ABCD} T[A,B,C,D] rho_A(1) rho_B(2) rho_C(3) rho_D(4),
  T = (1/4) sum_{a a' b b'} CS[a,b] CS[a',b'] E[(a,a'),(A,B)] E[(b,b'),(C,D)],
  E[(a,a'),(A,B)] from  dA[a] dA[a'](1,2) = rho_pp'(1) rho_qq'(2) - rho_pq'(1) rho_qp'(2)
                                           - rho_qp'(1) rho_pq'(2) + rho_qq'(1) rho_pp'(2).
  Symmetries: T[A,B,C,D] = T[B,A,C,D] = T[A,B,D,C];  singlet (real orbitals): = T[C,D,A,B].
  T is stored as the dense (npair^2 x npair^2) matrix T[A*npair+B, C*npair+D] when its COO form
  is small, else FACTORIZED (CS_t on the active pairs + the structural E) -- see save block.
Spectator traces (orthonormal NOs => INT rho_D d^3r = delta_D := [D is a diagonal pair (p,p)]):
  T3aab[A,B,C] = sum_D T d_D ;  T3abb[A,C,D] = sum_B T d_B ;
  T2aa[A,B] = sum_CD T d_C d_D = P_aa(1,2)   = fold(Gamma^{aa}) / 2
  T2ab[A,C] = sum_BD T d_B d_D = P_ab(1,3)   = fold(Gamma^{ab}) / 4
  T2bb[C,D] = sum_AB T d_A d_B               = fold(Gamma^{bb}) / 2
  T1a[A]    = sum_BCD T d d d = rho_alpha/2  = fold(gamma^a) / 2 ;  T1b likewise,
  Gamma^{st}_{pqrs} = <a+_{p s} a+_{r t} a_{s t} a_{q s}>  (chemist: E_ee = 1/2 sum Gamma_pqrs (pq|rs),
  (pq|rs) = eri[p,q,r,s] = INT phi_p phi_q(1) phi_r phi_s(2) / r12, the engine's convention).
The RDM identities above + product-operator (Slater 2x2-minor) route checks are the bookkeeping
gate of step (d) -- the thing the 2-MO code assumed away as |Phi0|^2 = D(12) D(34).
"""
from __future__ import annotations

import os
import sys
import time
from itertools import combinations
from typing import Dict, List, Sequence, Tuple

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
DATA = os.path.join(HERE, "data")

CFG = dict(Jb=2, Lb=1, npi=0, Jpi=0, Lpi=0, alpha=1.0, core2=(4.5, 1.6), tol=1e-11)
E_REF = -8.02335          # banked reference (plan; lih_core2exp_probe.run(2,1,0,0,0,1.0,core2=(4.5,1.6)))
TRUNC = 1e-3              # |c| > TRUNC kept
Z_A, Z_B = 3.0, 1.0
NELEC = 4
GATE_MHA = 1.0
PRIM_CACHE = os.path.join(DATA, "lih_marriage_phase0_prim.npz")
OUT = os.path.join(DATA, "lih_marriage_phase0_ref.npz")
OUT_EXACT = os.path.join(DATA, "lih_marriage_phase0_ref_exact.npz")
STO = {'coreLi': (3.0 - 5.0 / 16.0, 'A'), 'coreLi2': (4.5, 'A'), 'coreLi3': (1.6, 'A'),
       'Li2s': (0.65, 'A'), 'H1s': (1.0, 'B'), 'Hm': (0.70, 'B')}
GROUP = {'coreLi': 'core', 'coreLi2': 'core', 'coreLi3': 'core', 'Li2s': 'Li2s',
         'H1s': 'H', 'Hm': 'H'}


def _t(t0: float) -> str:
    return f"[{time.time() - t0:6.1f}s]"


# =========================================================================== #
# (a) reference: primitives, engine matrices, canonical orthogonalization
# =========================================================================== #
def build_primitives(verbose: bool = True) -> dict:
    """Primitive set + engine S / T / V_ne / ERI (float64 engine), cached; canonical MOs."""
    import lih_core2exp_probe as P
    import prolate_energy_ladder as L
    from prolate_allelectron_c4 import one_body_general
    from prolate_float_eri import build_eri_tensor_m_f
    from lih_vmc import extract_prim

    R = float(P.R)
    orbs, tags = P.lih_orbs(CFG['Jb'], CFG['Lb'], CFG['npi'], CFG['Jpi'], CFG['Lpi'],
                            CFG['alpha'], list(CFG['core2']))
    M = len(orbs)
    # tags are (key, j, l) with key = ('coreLi', 0, 0) etc. for STOs, ('bond', 0, Jb, Lb) for bond fns
    labels = []
    for (key, j, l) in tags:
        name = key[0] if isinstance(key, tuple) else key
        labels.append(f"bond_x{j}e{l}" if name == 'bond' else name)
    key = repr(CFG)
    cached = False
    if os.path.exists(PRIM_CACHE):
        z = np.load(PRIM_CACHE, allow_pickle=True)
        if str(z['cfg']) == key and int(z['M']) == M:
            S, h1, T1, eri = z['S'], z['h1'], z['T1'], z['eri']
            cached = True
            if verbose:
                print(f"  primitive engine matrices loaded from cache {PRIM_CACHE}")
    if not cached:
        t0 = time.time()
        S, h1 = one_body_general(orbs, R, Z_A, Z_B)          # h1 = T + V_ne (exact mpf -> f64)
        S0, T1 = one_body_general(orbs, R, 0.0, 0.0)         # Z=0 -> pure kinetic
        assert np.allclose(S, S0, atol=1e-14)
        if verbose:
            print(f"  one_body_general (S, T, V_ne) done {_t(t0)}", flush=True)
        eri = build_eri_tensor_m_f(orbs, R, verbose=verbose)   # float64 ERI engine (324 s)
        if verbose:
            print(f"  ERI engine done {_t(t0)}", flush=True)
        os.makedirs(DATA, exist_ok=True)
        np.savez(PRIM_CACHE, cfg=key, M=M, S=S, h1=h1, T1=T1, eri=eri)
    Vne = h1 - T1
    # re-basing (C) + canonical orthogonalization, exactly as assemble_rebased / build_lih_wavefunction
    C = L.build_C(tags, CFG['alpha'])
    So = C @ S @ C.T
    h1o = C @ h1 @ C.T
    T1o = C @ T1 @ C.T
    erio = np.einsum('pa,qb,rc,sd,abcd->pqrs', C, C, C, C, eri, optimize=True)
    w, U = np.linalg.eigh(So)
    keep = w > CFG['tol'] * w[-1]
    X = U[:, keep] / np.sqrt(w[keep])
    Mk = int(keep.sum())
    cond = w[-1] / w[keep].min()
    h1_f = X.T @ h1o @ X
    T1_f = X.T @ T1o @ X
    eri_f = np.einsum('ap,bq,cr,ds,abcd->pqrs', X, X, X, X, erio, optimize=True)
    Tmap = C.T @ X                                   # (M, Mk) primitive -> canonical MO
    prims = [extract_prim(o) for o in orbs]
    return dict(R=R, orbs=orbs, tags=tags, labels=labels, M=M, Mk=Mk, cond=cond,
                S=S, h1=h1, T1=T1, Vne=Vne, eri=eri, C=C, X=X, Tmap=Tmap,
                h1_f=h1_f, T1_f=T1_f, Vne_f=h1_f - T1_f, eri_f=eri_f, prims=prims,
                Vnn=Z_A * Z_B / R)


# =========================================================================== #
# FCI (same sparse connected-pair Lanczos as lih_vmc.fci_ground_vector, + returns H)
# =========================================================================== #
def fci_solve(h1: np.ndarray, eri: np.ndarray, K: int, nelec: int = NELEC):
    """Returns (E0, c, dets, H) in the _dets ordering; H the sparse CI Hamiltonian."""
    from prolate_allelectron_fci import _dets, _matel
    from fci_fast import _popcount64, _FOUR

    na = nb = nelec // 2
    dets = _dets(K, na, nb)
    nd = len(dets)
    masks = np.array([int(sum(1 << so for so in D)) for D in dets], dtype=np.uint64)
    rl, cl, vl = [], [], []
    diag = np.empty(nd)
    for a in range(nd):
        diag[a] = _matel(h1, eri, dets[a], dets[a])
        if a + 1 < nd:
            pc = _popcount64(masks[a + 1:] ^ masks[a])
            conn = np.nonzero(pc <= _FOUR)[0] + (a + 1)
            Da = dets[a]
            for b in conn.tolist():
                v = _matel(h1, eri, Da, dets[b])
                if v != 0.0:
                    rl.append(a); cl.append(b); vl.append(v)
    r = np.array(rl, dtype=np.int64); c_ = np.array(cl, dtype=np.int64); v = np.array(vl)
    ri = np.concatenate([np.arange(nd), r, c_])
    ci = np.concatenate([np.arange(nd), c_, r])
    vi = np.concatenate([diag, v, v])
    H = sp.csr_matrix((vi, (ri, ci)), shape=(nd, nd))
    ev, evec = eigsh(H, k=1, which='SA', return_eigenvectors=True, tol=0)
    c = evec[:, 0]
    if c[np.argmax(np.abs(c))] < 0:
        c = -c
    return float(ev[0]), c, dets, H


def _perm_sign_inv(seq: Sequence[int]) -> float:
    n = len(seq); inv = 0
    for i in range(n):
        for j in range(i + 1, n):
            if seq[i] > seq[j]:
                inv += 1
    return -1.0 if (inv & 1) else 1.0


def fold_cs(c: np.ndarray, K: int) -> Tuple[np.ndarray, list, np.ndarray]:
    """CS[a,b] = c_I sigma_I (grouped-order coefficient matrix), apairs, sigma."""
    apairs = list(combinations(range(K), 2))
    n_ap = len(apairs)
    Cmat = c.reshape(n_ap, n_ap)
    sigma = np.empty((n_ap, n_ap))
    for ia, a in enumerate(apairs):
        for ib, b in enumerate(apairs):
            sigma[ia, ib] = _perm_sign_inv([2 * p for p in a] + [2 * q + 1 for q in b])
    return Cmat * sigma, apairs, sigma


# =========================================================================== #
# pair-index bookkeeping
# =========================================================================== #
def pair_index(K: int):
    pairs = [(p, q) for p in range(K) for q in range(p, K)]
    idx = {pq: i for i, pq in enumerate(pairs)}
    return pairs, idx


def fold_map(K: int, pairs, idx) -> np.ndarray:
    """F (K*K, npair): F[p*K+q, A] = 1 iff A = (min,max)(p,q).  fold(M) = F^T vec(M)."""
    F = np.zeros((K * K, len(pairs)))
    for p in range(K):
        for q in range(K):
            F[p * K + q, idx[(min(p, q), max(p, q))]] = 1.0
    return F


def build_Eop(K: int, apairs, aidx) -> np.ndarray:
    """Eop[p,q,a,a'] = <a| a+_p a_q |a'> on the 2-electron same-spin pair determinants (sorted)."""
    n_ap = len(apairs)
    E = np.zeros((K, K, n_ap, n_ap))
    for ia2, (r, s) in enumerate(apairs):
        for q, tleft, sq in ((r, s, +1.0), (s, r, -1.0)):     # a_q |r s> = sq |tleft>
            for p in range(K):
                if p == tleft:
                    continue
                sp_ = +1.0 if p < tleft else -1.0             # a+_p |t> = sp |sorted(p,t)>
                E[p, q, aidx[(min(p, tleft), max(p, tleft))], ia2] += sq * sp_
    return E


def build_Eexp(K: int, apairs, pairs, idx) -> sp.csr_matrix:
    """E[(a,a'),(A,B)] sparse (n_ap^2, npair^2): dA[a]dA[a'](1,2) = sum_AB E rho_A(1) rho_B(2)."""
    n_ap = len(apairs); npair = len(pairs)
    rows, cols, vals = [], [], []

    def ix(x, y):
        return idx[(min(x, y), max(x, y))]

    for ia, (p, q) in enumerate(apairs):
        for ib, (p2, q2) in enumerate(apairs):
            row = ia * n_ap + ib
            for sgn, (x1, y1), (x2, y2) in ((+1.0, (p, p2), (q, q2)), (-1.0, (p, q2), (q, p2)),
                                             (-1.0, (q, p2), (p, q2)), (+1.0, (q, q2), (p, p2))):
                rows.append(row); cols.append(ix(x1, y1) * npair + ix(x2, y2)); vals.append(sgn)
    return sp.csr_matrix((vals, (rows, cols)), shape=(n_ap * n_ap, npair * npair))


def rdms(CS: np.ndarray, Eop: np.ndarray):
    """Spin-resolved 1-RDMs and 2-RDMs (chemist) of Psi = sum CS[a,b] |a>|b>, Sum CS^2 = 1."""
    K = Eop.shape[0]
    Da = CS @ CS.T
    Db = CS.T @ CS
    ga = np.einsum('pqaA,aA->pq', Eop, Da)
    gb = np.einsum('pqaA,aA->pq', Eop, Db)
    I = np.eye(K)

    def same_spin(D, g):
        Y = np.einsum('pqaA,aB->pqAB', Eop, D)
        G = np.einsum('pqAB,rsAB->pqrs', Y, Eop, optimize=True)
        return G - np.einsum('qr,ps->pqrs', I, g)

    Gaa = same_spin(Da, ga)
    Gbb = same_spin(Db, gb)
    Z = np.einsum('pqaA,ab->pqAb', Eop, CS)
    W2 = np.einsum('pqAb,AB->pqbB', Z, CS)
    Gab = np.einsum('pqbB,rsbB->pqrs', W2, Eop, optimize=True)
    return ga, gb, Gaa, Gbb, Gab


def energy_from_rdm(ga, gb, Gaa, Gbb, Gab, h1, eri, Vnn):
    g = ga + gb
    G = Gaa + Gbb + Gab + Gab.transpose(2, 3, 0, 1)
    return float(np.einsum('pq,pq', g, h1) + 0.5 * np.einsum('pqrs,pqrs', G, eri) + Vnn)


# =========================================================================== #
# tensor T (factorized build) + traces + product-operator checks
# =========================================================================== #
def build_T(CS: np.ndarray, Eexp: sp.csr_matrix, n_ap: int, npair: int) -> np.ndarray:
    """Dense T[A*npair+B, C*npair+D] = 1/4 sum CS[a,b]CS[a',b'] E[(a,a'),AB] E[(b,b'),CD].
    Restricted to the ACTIVE alpha / beta pairs (rows/cols of CS that are nonzero)."""
    act_a = np.nonzero(np.any(CS != 0.0, axis=1))[0]
    act_b = np.nonzero(np.any(CS != 0.0, axis=0))[0]
    CSr = CS[np.ix_(act_a, act_b)]
    ra = (act_a[:, None] * n_ap + act_a[None, :]).ravel()      # rows of Eexp for (a,a')
    rb = (act_b[:, None] * n_ap + act_b[None, :]).ravel()
    Ea = Eexp[ra]                                             # (n_a^2, npair^2)
    Eb = Eexp[rb]
    Q = np.einsum('ab,AB->aAbB', CSr, CSr).reshape(len(act_a) ** 2, len(act_b) ** 2)
    G = (Ea.T @ Q)                                            # (npair^2, n_b^2) dense
    T = 0.25 * np.asarray((Eb.T @ G.T).T)                     # (npair^2, npair^2)
    return T, act_a, act_b


def traces(T4: np.ndarray, d: np.ndarray) -> dict:
    return dict(
        T3aab=np.einsum('ABCD,D->ABC', T4, d),
        T3abb=np.einsum('ABCD,B->ACD', T4, d),
        T2aa=np.einsum('ABCD,C,D->AB', T4, d, d),
        T2ab=np.einsum('ABCD,B,D->AC', T4, d, d),
        T2bb=np.einsum('ABCD,A,B->CD', T4, d, d),
        T1a=np.einsum('ABCD,B,C,D->A', T4, d, d, d),
        T1b=np.einsum('ABCD,A,B,D->C', T4, d, d, d),
    )


def product_checks(T4, tr, CS, Eop, apairs, pairs, K, v: np.ndarray) -> Dict[str, Tuple[float, float]]:
    """<v(1)..v(k)> over P via T (left) vs the 2x2-minor / pair-operator route (right).
    v_A := INT rho_A v = v_pq for A=(p,q)  (NOT the fold -- rho_A carries no factor 2)."""
    n_ap = len(apairs)
    det2 = np.empty((n_ap, n_ap))
    for ia, (p, q) in enumerate(apairs):
        for ib, (p2, q2) in enumerate(apairs):
            det2[ia, ib] = v[p, p2] * v[q, q2] - v[p, q2] * v[q, p2]
    Ev = np.einsum('tu,tuaA->aA', v, Eop)
    vA = np.array([v[p, q] for (p, q) in pairs])
    Da = CS @ CS.T
    out = {}
    out['v1'] = (float(tr['T1a'] @ vA), float(0.5 * np.sum(Da * Ev)))
    out['v1v2'] = (float(vA @ tr['T2aa'] @ vA), float(np.sum(Da * det2)))
    out['v1v3'] = (float(vA @ tr['T2ab'] @ vA), float(0.25 * np.sum(CS * (Ev @ CS @ Ev.T))))
    out['v1v2v3'] = (float(np.einsum('ABC,A,B,C', tr['T3aab'], vA, vA, vA)),
                     float(0.5 * np.sum(CS * (det2 @ CS @ Ev.T))))
    out['v1v3v4'] = (float(np.einsum('ACD,A,C,D', tr['T3abb'], vA, vA, vA)),
                     float(0.5 * np.sum(CS * (Ev @ CS @ det2.T))))
    out['v1v2v3v4'] = (float(np.einsum('ABCD,A,B,C,D', T4, vA, vA, vA, vA)),
                       float(np.sum(CS * (det2 @ CS @ det2.T))))
    return out


# =========================================================================== #
# grid side (geovac/lih_r12ci/kernels.py axial grid + hVee.neumann_potential)
# =========================================================================== #
def grid_setup(exact: bool = False):
    import geovac.lih_r12ci.kernels as KG
    from geovac.lih_r12ci.hVee import neumann_potential
    from geovac.lih_r12ci import energy as EN
    KG.USE_EXACT_NEUMANN = bool(exact)          # Phase 0b switch (read by neumann_potential at call time)
    pts = np.stack([KG.RHO_CYL.ravel(), np.zeros(KG.NG), KG.ZC.ravel()], axis=-1)
    cosA = ((1.0 + KG.Xg * KG.Eg) / (KG.Xg + KG.Eg)).ravel()          # cos(theta_A)
    vne = (-Z_A / KG.rA - Z_B / KG.rB).ravel()
    return dict(KG=KG, neumann=neumann_potential, R_engine=float(EN.R), pts=pts, cosA=cosA,
                vne=vne, rA=KG.rA.ravel(), rB=KG.rB.ravel(), geo=KG.geo_f,
                gint=KG.grid_int, shape=(KG.NXI, KG.NETA),
                params=dict(NXI=KG.NXI, NETA=KG.NETA, NPHI=KG.NPHI, xi_max=float(KG.xi_max),
                            a=float(KG.a), ZA_grid=float(KG.ZA), ZB_grid=float(KG.ZB)))


def prim_on_grid(prims, R, pts):
    """Values, gradients, Laplacians of the M primitives on the grid (lih_vmc.orbital_vgl)."""
    from lih_vmc import orbital_vgl
    M = len(prims); NG = pts.shape[0]
    val = np.empty((M, NG)); grad = np.empty((M, NG, 3)); lap = np.empty((M, NG))
    for a, p in enumerate(prims):
        v, g, l = orbital_vgl(p, pts, R)
        val[a] = v.real; grad[a] = g.real; lap[a] = l.real
    return val, grad, lap


def prim_direct(prims, KG):
    """Independent value check: N xi^j g(eta) e^{-alpha xi} straight on (Xg, Eg)."""
    out = np.empty((len(prims), KG.NG))
    for a, p in enumerate(prims):
        g = np.polynomial.polynomial.polyval(KG.Eg, p['eta'])
        out[a] = (p['N'] * KG.Xg ** p['j'] * g * np.exp(-p['alpha'] * KG.Xg)).ravel()
    return out


def one_body_grid(val, grad, lap, geo, vne, a3):
    """S, T (Laplacian form), T (gradient form), V_ne matrices over the given orbital set."""
    w = 2 * np.pi * a3 * geo
    S = (val * w) @ val.T
    Tl = -0.5 * (val * w) @ lap.T
    Tg = 0.5 * np.einsum('aGd,G,bGd->ab', grad, w, grad)
    V = (val * w * vne) @ val.T
    return S, 0.5 * (Tl + Tl.T), Tg, V


def coulomb_grid(rho: np.ndarray, neumann, shape, geo, a3, LMAX=34):
    """J[A,B] = INT rho_A V[rho_B] with V the prolate-Neumann potential; returns (J, Vfields)."""
    npair = rho.shape[0]
    Vf = np.empty_like(rho)
    for B in range(npair):
        Vf[B] = neumann(rho[B].reshape(shape), LMAX).reshape(-1)
    w = 2 * np.pi * a3 * geo
    J = (rho * w) @ Vf.T
    return J, Vf


def hartree_1s(r, Z):
    r = np.maximum(r, 1e-30)
    return (1.0 / r) * (1.0 - (1.0 + Z * r) * np.exp(-2.0 * Z * r))


# --------------------------------------------------------------------------- #
# DIAGNOSTIC ONLY (not the deliverable grid): the same prolate-Neumann potential formula as
# hVee.neumann_potential on an arbitrary Gauss-Legendre (xi, eta) grid, to see how the <V_ee>
# error scales with NXI / NETA / LMAX / xi_max.  kernels.py is NOT modified.
# --------------------------------------------------------------------------- #
def neumann_general(D, xi1d, wxi, eta1d, weta, LMAX, R):
    from scipy.special import lqn, eval_legendre
    a = R / 2.0
    NXI, NETA = len(xi1d), len(eta1d)
    XI, ETA = np.meshgrid(xi1d, eta1d, indexing='ij')
    JAC = XI ** 2 - ETA ** 2
    Qtab = np.array([lqn(LMAX, x)[0] for x in xi1d])
    Pxi = np.array([eval_legendre(l, xi1d) for l in range(LMAX + 1)])
    minidx = np.minimum.outer(np.arange(NXI), np.arange(NXI))
    maxidx = np.maximum.outer(np.arange(NXI), np.arange(NXI))
    W = JAC * D
    V = np.zeros((NXI, NETA))
    for l in range(LMAX + 1):
        Pl = eval_legendre(l, eta1d)
        g_l = (W * Pl[None, :]) @ weta
        Kl = Pxi[l][minidx] * Qtab[:, l][maxidx]
        V += (2 * l + 1) * np.outer(Kl @ (wxi * g_l), Pl)
    return (2.0 / R) * (2 * np.pi) * a ** 3 * V


def vee_on_grid(prims, T_no, pairs, G_tot, eri_ref, R, NXI, NETA, LMAX, xi_max, sto_ctrl):
    """<V_ee> of the truncated vector on a scratch GL grid; returns (dVee_mHa, {label: selfJ rel err})."""
    from numpy.polynomial.legendre import leggauss
    a = R / 2.0
    xg, wxg = leggauss(NXI); xi1d = 1.0 + 0.5 * (xg + 1.0) * (xi_max - 1.0); wxi = 0.5 * (xi_max - 1.0) * wxg
    eta1d, weta = leggauss(NETA)
    XI, ETA = np.meshgrid(xi1d, eta1d, indexing='ij')
    JAC = XI ** 2 - ETA ** 2
    w = (2 * np.pi * a ** 3 * np.outer(wxi, weta) * JAC).ravel()
    pv = np.empty((len(prims), XI.size))
    for i, p in enumerate(prims):
        g = np.polynomial.polynomial.polyval(ETA, p['eta'])
        pv[i] = (p['N'] * XI ** p['j'] * g * np.exp(-p['alpha'] * XI)).ravel()
    K = T_no.shape[1]
    nv = T_no.T @ pv
    npair = len(pairs)
    rho = np.empty((npair, XI.size))
    for A, (p, q) in enumerate(pairs):
        rho[A] = nv[p] * nv[q]
    Vf = np.empty_like(rho)
    for B in range(npair):
        Vf[B] = neumann_general(rho[B].reshape(XI.shape), xi1d, wxi, eta1d, weta, LMAX, R).ravel()
    J = (rho * w) @ Vf.T
    Fm = np.zeros((K * K, npair))
    for A, (p, q) in enumerate(pairs):
        Fm[p * K + q, A] = 1.0; Fm[q * K + p, A] = 1.0
    J4 = (Fm @ J @ Fm.T).reshape(K, K, K, K)
    dVee = 0.5 * float(np.einsum('pqrs,pqrs', G_tot, J4 - eri_ref)) * 1e3
    ctrl = {}
    for lab, (i, z) in sto_ctrl.items():
        d = pv[i] ** 2
        Vs = neumann_general(d.reshape(XI.shape), xi1d, wxi, eta1d, weta, LMAX, R).ravel()
        ctrl[lab] = (float(np.sum(w * d * Vs)) - 5 * z / 8) / (5 * z / 8)
    return dVee, ctrl


# =========================================================================== #
def main(exact: bool = False):
    T0 = time.time()
    np.set_printoptions(linewidth=140, precision=6, suppress=True)
    print("=" * 96)
    print("LiH MARRIAGE -- PHASE 0: sigma core-enriched reference, truncation, pair-index tensor T, GATE G0")
    print(f"  CFG={CFG}   E_REF={E_REF}   TRUNC={TRUNC}   gate {GATE_MHA} mHa")
    print(f"  prolate-Neumann radial step: {'EXACT ordered-integral operator (Phase 0b, neumann_exact.py)' if exact else 'LEGACY cumulative GL sum'}")
    print("=" * 96, flush=True)

    # ------------------------------------------------------------------ (a)
    print("\n(a) REFERENCE")
    ref = build_primitives()
    R, M, Mk = ref['R'], ref['M'], ref['Mk']
    print(f"  primitives M={M}: " + ", ".join(f"{i}:{l}" for i, l in enumerate(ref['labels'])))
    print(f"  canonical-orthogonalized Mk={Mk}  cond(S_o)={ref['cond']:.2e}   {_t(T0)}", flush=True)
    E0, c, dets, H = fci_solve(ref['h1_f'], ref['eri_f'], Mk)
    E_full = E0 + ref['Vnn']
    nd = len(dets)
    print(f"  FCI (canonical MO basis): n_det={nd}  E_elec={E0:.6f}  E_tot={E_full:.6f}  "
          f"vs E_REF {E_REF}  (dE={(E_full - E_REF) * 1e3:+.4f} mHa)   {_t(T0)}", flush=True)
    G = grid_setup(exact)
    print(f"  R check: probe/engine-pipeline R={R}  lih_r12ci.energy.R={G['R_engine']}  "
          f"{'MATCH' if abs(R - G['R_engine']) < 1e-12 else 'MISMATCH'}   grid a={G['params']['a']}")
    print(f"  kernels.py grid: NXI={G['params']['NXI']} NETA={G['params']['NETA']} "
          f"NPHI={G['params']['NPHI']} xi_max={G['params']['xi_max']:.4f} "
          f"(sized for ZA={G['params']['ZA_grid']}, ZB={G['params']['ZB_grid']})   {_t(T0)}", flush=True)
    print("  canonical-basis truncation census (for the record):")
    for thr in (1e-2, 3e-3, 1e-3, 3e-4, 1e-4):
        m = np.abs(c) > thr
        print(f"    |c|>{thr:.0e}: n_kept={m.sum():5d}  norm kept={np.sum(c[m] ** 2):.6f}")

    # natural orbitals of the FULL CI ------------------------------------------------
    pairs, idx = pair_index(Mk)
    npair = len(pairs)
    apairs = list(combinations(range(Mk), 2)); aidx = {a: i for i, a in enumerate(apairs)}
    n_ap = len(apairs)
    Eop = build_Eop(Mk, apairs, aidx)
    F = fold_map(Mk, pairs, idx)
    CS_full, _, sigma = fold_cs(c, Mk)
    ga_c = np.einsum('pqaA,aA->pq', Eop, CS_full @ CS_full.T)
    gb_c = np.einsum('pqaA,aA->pq', Eop, CS_full.T @ CS_full)
    occ, U = np.linalg.eigh(ga_c + gb_c)
    order = np.argsort(occ)[::-1]
    occ = occ[order]; U = U[:, order]
    for p in range(Mk):                                     # sign convention: largest component +
        if U[np.argmax(np.abs(U[:, p])), p] < 0:
            U[:, p] = -U[:, p]
    print(f"  NO occupations (full CI, spin-summed; sum={occ.sum():.6f}): {np.round(occ, 5)}")
    h1_no = U.T @ ref['h1_f'] @ U
    T1_no = U.T @ ref['T1_f'] @ U
    Vne_no = h1_no - T1_no
    eri_no = np.einsum('ap,bq,cr,ds,abcd->pqrs', U, U, U, U, ref['eri_f'], optimize=True)
    T_no = ref['Tmap'] @ U                                  # (M, Mk) primitive -> NO
    E0n, cn, dets_n, Hn = fci_solve(h1_no, eri_no, Mk)
    E_full_no = E0n + ref['Vnn']
    print(f"  FCI re-solved in the NO basis: E_tot={E_full_no:.6f}  (dE vs canonical "
          f"{(E_full_no - E_full) * 1e6:+.3f} uHa)   {_t(T0)}", flush=True)
    # NO character (Mulliken populations by primitive group)
    SPT = ref['S'] @ T_no
    pops = {g: np.zeros(Mk) for g in ('core', 'Li2s', 'H', 'bond')}
    for a, lab in enumerate(ref['labels']):
        g = GROUP.get(lab, 'bond')
        pops[g] += T_no[a, :] * SPT[a, :]
    char = []
    for p in range(Mk):
        dom = max(pops, key=lambda g: pops[g][p])
        char.append(dom)
    print("  NO character (Mulliken population by primitive group; core = Li 1s STOs):")
    for p in range(Mk):
        print(f"    NO{p:2d}: occ={occ[p]:.5f}  core={pops['core'][p]:+.3f} Li2s={pops['Li2s'][p]:+.3f} "
              f"H={pops['H'][p]:+.3f} bond={pops['bond'][p]:+.3f}  -> {char[p]}")
    print("  NO-basis truncation census:")
    for thr in (1e-2, 3e-3, 1e-3, 3e-4, 1e-4):
        m = np.abs(cn) > thr
        print(f"    |c|>{thr:.0e}: n_kept={m.sum():5d}  norm kept={np.sum(cn[m] ** 2):.6f}")

    # ------------------------------------------------------------------ (b)
    print(f"\n(b) TRUNCATION at |c| > {TRUNC} (NO basis)")
    kept = np.nonzero(np.abs(cn) > TRUNC)[0]
    norm_kept = float(np.sum(cn[kept] ** 2))
    ct = np.zeros_like(cn); ct[kept] = cn[kept] / np.sqrt(norm_kept)
    E_ray = float(ct @ (Hn @ ct)) + ref['Vnn']
    CS, _, _ = fold_cs(ct, Mk)
    ga, gb, Gaa, Gbb, Gab = rdms(CS, Eop)
    E_trunc = energy_from_rdm(ga, gb, Gaa, Gbb, Gab, h1_no, eri_no, ref['Vnn'])
    g_tot = ga + gb
    G_tot = Gaa + Gbb + Gab + Gab.transpose(2, 3, 0, 1)
    print(f"  n_det kept = {len(kept)} of {nd};  norm retained before renormalisation = {norm_kept:.6f}")
    print(f"  E_trunc (engine h1/eri via RDMs) = {E_trunc:.6f}   Rayleigh c_t^T H c_t = {E_ray:.6f}   "
          f"(|diff| {abs(E_trunc - E_ray):.1e})")
    print(f"  E_trunc - E_full = {(E_trunc - E_full_no) * 1e3:+.3f} mHa   ;  Tr gamma = {np.trace(g_tot):.6f}  "
          f"Tr Gamma/(N-1) = {np.einsum('ppqq', G_tot) / 3:.6f}   {_t(T0)}", flush=True)
    top = kept[np.argsort(np.abs(ct[kept]))[::-1][:8]]
    print("  leading determinants (spin-orbital tuples; so=2*NO+spin):")
    for I in top:
        print(f"    c={ct[I]:+.5f}  {dets_n[I]}")

    # ------------------------------------------------------------------ (d) tensor T + traces + checks
    print(f"\n(d) COEFFICIENT TENSOR T over pair indices (npair={npair}, n_ap={n_ap})")
    Eexp = build_Eexp(Mk, apairs, pairs, idx)
    T, act_a, act_b = build_T(CS, Eexp, n_ap, npair)
    T4 = T.reshape(npair, npair, npair, npair)
    d = np.array([1.0 if p == q else 0.0 for (p, q) in pairs])
    tr = traces(T4, d)
    nnz = int(np.count_nonzero(T))
    print(f"  active alpha pairs {len(act_a)}, beta pairs {len(act_b)};  T dense shape {T.shape} "
          f"({T.nbytes / 1e6:.0f} MB), nonzero fraction {nnz / T.size:.3f}")
    print(f"  T3aab {tr['T3aab'].shape}, T2aa/T2ab/T2bb {tr['T2aa'].shape}, T1a {tr['T1a'].shape}")
    # symmetry checks
    sym_12 = np.max(np.abs(T4 - T4.transpose(1, 0, 2, 3)))
    sym_34 = np.max(np.abs(T4 - T4.transpose(0, 1, 3, 2)))
    sym_ab = np.max(np.abs(T4 - T4.transpose(2, 3, 0, 1)))
    norm_T = float(np.einsum('ABCD,A,B,C,D', T4, d, d, d, d))
    print(f"  symmetries: max|T-T(1<->2)|={sym_12:.1e}  max|T-T(3<->4)|={sym_34:.1e}  "
          f"max|T-T(ab<->cd)| (singlet)={sym_ab:.1e};  INT P = {norm_T:.12f}")
    # trace checks vs RDMs (exact bookkeeping)
    fold1 = lambda Mx: F.T @ Mx.reshape(-1)
    fold2 = lambda G4: F.T @ G4.reshape(Mk * Mk, Mk * Mk) @ F
    chk = {
        'T1a  vs fold(gamma^a)/2': np.max(np.abs(tr['T1a'] - fold1(ga) / 2)),
        'T1b  vs fold(gamma^b)/2': np.max(np.abs(tr['T1b'] - fold1(gb) / 2)),
        'T2aa vs fold(Gamma^aa)/2': np.max(np.abs(tr['T2aa'] - fold2(Gaa) / 2)),
        'T2bb vs fold(Gamma^bb)/2': np.max(np.abs(tr['T2bb'] - fold2(Gbb) / 2)),
        'T2ab vs fold(Gamma^ab)/4': np.max(np.abs(tr['T2ab'] - fold2(Gab) / 4)),
        'T3aab traced over C vs T2aa': np.max(np.abs(tr['T3aab'] @ d - tr['T2aa'])),
        'T3abb traced over D vs T2ab': np.max(np.abs(tr['T3abb'] @ d - tr['T2ab'])),
    }
    print("  spectator traces vs RDMs from the CI vector (max |Delta|):")
    for k, v in chk.items():
        print(f"    {k:32s} {v:.2e}")
    trace_max = max(chk.values())
    # product-operator checks (independent 2x2-minor route) for T3 and the full T
    rng = np.random.default_rng(7)
    prod_max = 0.0
    print("  product-operator checks <v(1)...v(k)> : T-route vs 2x2-minor route (rel |Delta|):")
    for name, v in (('random v #1', None), ('random v #2', None), ('v = V_ne(NO)', Vne_no)):
        if v is None:
            v = rng.standard_normal((Mk, Mk)); v = 0.5 * (v + v.T)
        pc = product_checks(T4, tr, CS, Eop, apairs, pairs, Mk, v)
        line = []
        for k, (lhs, rhs) in pc.items():
            rel = abs(lhs - rhs) / max(abs(rhs), 1e-300)
            prod_max = max(prod_max, rel)
            line.append(f"{k}:{rel:.1e}")
        print(f"    {name:14s} " + "  ".join(line))
    print(f"  => trace check max|Delta| = {trace_max:.2e};  product checks max rel = {prod_max:.2e}   {_t(T0)}",
          flush=True)

    # ------------------------------------------------------------------ (c) grid: orbitals + pair densities
    print("\n(c) ORBITALS AND PAIR DENSITIES ON THE kernels.py GRID")
    KG = G['KG']; a3 = G['params']['a'] ** 3
    pval, pgrad, plap = prim_on_grid(ref['prims'], R, G['pts'])
    pdir = prim_direct(ref['prims'], KG)
    print(f"  orbital_vgl value vs direct prolate formula: max rel {np.max(np.abs(pval - pdir)) / np.max(np.abs(pval)):.1e}")
    Sg_p, Tg_p, Tgg_p, Vg_p = one_body_grid(pval, pgrad, plap, G['geo'], G['vne'], a3)
    print(f"  primitive one-body on grid vs ENGINE: max|dS|={np.max(np.abs(Sg_p - ref['S'])):.2e}  "
          f"max|dT|={np.max(np.abs(Tg_p - ref['T1'])):.2e} (grad-form {np.max(np.abs(Tgg_p - ref['T1'])):.2e})  "
          f"max|dVne|={np.max(np.abs(Vg_p - ref['Vne'])):.2e}")
    # STO closed-form controls (per primitive)
    print("  STO primitive controls on the grid (norm, <1/r_own>=zeta, <T>=zeta^2/2, self-Coulomb 5zeta/8 [Neumann]):")
    rho_p = np.empty((npair, KG.NG))
    for A, (p, q) in enumerate(pairs):
        rho_p[A] = pval[p] * pval[q]
    Jp, _ = coulomb_grid(rho_p, G['neumann'], G['shape'], G['geo'], a3)
    J4p = (F @ Jp @ F.T).reshape(Mk, Mk, Mk, Mk)
    dJp = J4p - ref['eri']
    Jp_asym = np.max(np.abs(Jp - Jp.T))
    for i, lab in enumerate(ref['labels']):
        if lab not in STO:
            continue
        z, cen = STO[lab]
        r_own = G['rA'] if cen == 'A' else G['rB']
        r_oth = G['rB'] if cen == 'A' else G['rA']
        rho = pval[i] ** 2
        nrm = G['gint'](rho)
        inv_own = G['gint'](rho / r_own)
        inv_oth = G['gint'](rho / r_oth)
        pt = (1.0 / R) * (1.0 - (1.0 + z * R) * np.exp(-2.0 * z * R))
        tk = Tg_p[i, i]
        Aii = idx[(i, i)]
        Jself = Jp[Aii, Aii]
        Vh = hartree_1s(r_own, z)
        Jhart = G['gint'](rho * Vh)
        print(f"    {lab:8s} z={z:<6} norm-1={nrm - 1:+.1e}  <1/r_own>-z={inv_own - z:+.1e}  "
              f"<1/r_oth>-pt={inv_oth - pt:+.1e}  T-z^2/2={tk - z * z / 2:+.1e}  "
              f"J_Neumann-5z/8={Jself - 5 * z / 8:+.2e} (Hartree-dress {Jhart - 5 * z / 8:+.1e})  "
              f"engine {ref['eri'][i, i, i, i] - 5 * z / 8:+.1e}")
    print(f"  primitive-pair Coulomb J[A,B] (78x78, Neumann LMAX=34): asym max|J-J^T|={Jp_asym:.1e};  "
          f"vs ENGINE eri: max|dJ|={np.max(np.abs(dJp)):.2e}")
    worst = np.unravel_index(np.argmax(np.abs(dJp)), dJp.shape)
    print(f"    worst (pq|rs) = {tuple(ref['labels'][k] for k in worst)}: grid {J4p[worst]:.6f} engine {ref['eri'][worst]:.6f}")
    per_prim = [np.max(np.abs(dJp[i])) for i in range(M)]
    print("    max|dJ| by primitive index: " + "  ".join(f"{ref['labels'][i]}:{per_prim[i]:.1e}" for i in range(M)))
    # NO orbitals on the grid (deliverable) ---------------------------------
    nval = T_no.T @ pval
    ngrad = np.einsum('ap,aGd->pGd', T_no, pgrad)
    nlap = T_no.T @ plap
    Sg, Tg, Tgg, Vg = one_body_grid(nval, ngrad, nlap, G['geo'], G['vne'], a3)
    rho_no = np.empty((npair, KG.NG))
    for A, (p, q) in enumerate(pairs):
        rho_no[A] = nval[p] * nval[q]
    J_no, V_no = coulomb_grid(rho_no, G['neumann'], G['shape'], G['geo'], a3)
    J4 = (F @ J_no @ F.T).reshape(Mk, Mk, Mk, Mk)
    J4_via_prim = np.einsum('ap,bq,cr,ds,abcd->pqrs', T_no, T_no, T_no, T_no, J4p, optimize=True)
    print(f"  NO pair densities: {npair} fields x {KG.NG} pts;  NO-route J vs primitive-route J transformed: "
          f"max|d|={np.max(np.abs(J4 - J4_via_prim)):.1e} (consistency)")
    norm_err = np.diag(Sg) - 1.0
    offd = Sg - np.diag(np.diag(Sg))
    print("  per-NO normalisation error INT rho_pp - 1: " + "  ".join(f"NO{p}:{norm_err[p]:+.1e}" for p in range(Mk)))
    print(f"  max off-diagonal INT rho_pq (p!=q) = {np.max(np.abs(offd)):.1e}")

    # ------------------------------------------------------------------ (e) GATE G0
    print("\n(e) GATE G0: grid vs engine, truncated vector (NO basis)")
    T_grid = float(np.einsum('pq,pq', g_tot, Tg)); T_gridg = float(np.einsum('pq,pq', g_tot, Tgg))
    T_eng = float(np.einsum('pq,pq', g_tot, T1_no))
    V_grid = float(np.einsum('pq,pq', g_tot, Vg)); V_eng = float(np.einsum('pq,pq', g_tot, Vne_no))
    Vee_grid = float(0.5 * np.einsum('pqrs,pqrs', G_tot, J4)); Vee_eng = float(0.5 * np.einsum('pqrs,pqrs', G_tot, eri_no))
    tot_grid = T_grid + V_grid + Vee_grid + ref['Vnn']; tot_eng = T_eng + V_eng + Vee_eng + ref['Vnn']
    rows = [('<T> (Laplacian form)', T_grid, T_eng), ('<T> (gradient form)', T_gridg, T_eng),
            ('<V_ne>', V_grid, V_eng), ('<V_ee>', Vee_grid, Vee_eng),
            ('total (+V_nn)', tot_grid, tot_eng)]
    print(f"  {'term':22s} {'grid':>14s} {'engine':>14s} {'Delta (mHa)':>13s}")
    for name, gval, eval_ in rows:
        print(f"  {name:22s} {gval:14.6f} {eval_:14.6f} {(gval - eval_) * 1e3:+13.4f}")
    print(f"  (engine total vs E_trunc {E_trunc:.6f}: {abs(tot_eng - E_trunc):.1e})")
    dT = abs(T_grid - T_eng) * 1e3; dV = abs(V_grid - V_eng) * 1e3; dVee = abs(Vee_grid - Vee_eng) * 1e3
    g0_pass = (dT <= GATE_MHA) and (dV <= GATE_MHA) and (dVee <= GATE_MHA)
    # diagnostics: error split by NO-pair class
    print("  diagnostics -- contributions to Delta by NO-pair class (mHa):")
    cls = {}
    for p in range(Mk):
        for q in range(Mk):
            key = tuple(sorted((char[p], char[q])))
            e = cls.setdefault(key, [0.0, 0.0])
            e[0] += g_tot[p, q] * (Tg[p, q] - T1_no[p, q]) * 1e3
            e[1] += g_tot[p, q] * (Vg[p, q] - Vne_no[p, q]) * 1e3
    for key, (eT, eV) in sorted(cls.items()):
        print(f"    {str(key):22s} dT={eT:+9.4f}  dVne={eV:+9.4f}")
    dJ = J4 - eri_no
    contrib = 0.5 * G_tot * dJ * 1e3
    cls2 = {}
    for p in range(Mk):
        for q in range(Mk):
            for r in range(Mk):
                for s in range(Mk):
                    key = (tuple(sorted((char[p], char[q]))), tuple(sorted((char[r], char[s]))))
                    key = tuple(sorted(key))
                    cls2[key] = cls2.get(key, 0.0) + contrib[p, q, r, s]
    big = sorted(cls2.items(), key=lambda kv: -abs(kv[1]))[:8]
    print("  largest <V_ee> Delta contributions by (pair class, pair class) (mHa): "
          + "; ".join(f"{k}:{v:+.3f}" for k, v in big))
    print(f"  max|J_grid-J_engine| over NO pairs = {np.max(np.abs(dJ)):.2e}; max|dT_pq|={np.max(np.abs(Tg - T1_no)):.1e}; "
          f"max|dVne_pq|={np.max(np.abs(Vg - Vne_no)):.1e}")
    # per-pair <V_ee> error attribution (which pair densities carry it)
    print("  NO-pair Coulomb diagnostics (grid Neumann vs engine): self-Coulomb (pp|pp) relative error per NO:")
    print("    " + "  ".join(f"NO{p}:{(J4[p, p, p, p] - eri_no[p, p, p, p]) / eri_no[p, p, p, p]:+.1e}" for p in range(Mk)))
    flat = np.argsort(np.abs(contrib).ravel())[::-1][:10]
    print("  top-10 |1/2 Gamma_pqrs dJ_pqrs| contributions (mHa) [NO indices, classes, J_grid, J_engine, rel]:")
    for f_ in flat:
        p, q, r, s = np.unravel_index(f_, contrib.shape)
        print(f"    ({p},{q}|{r},{s}) [{char[p]}{char[q]}|{char[r]}{char[s]}]  {contrib[p, q, r, s]:+.4f}  "
              f"J_grid={J4[p, q, r, s]:.6f} J_eng={eri_no[p, q, r, s]:.6f} rel={dJ[p, q, r, s] / max(abs(eri_no[p, q, r, s]), 1e-300):+.1e}")
    # angular content about centre A of the leaf-candidate pair densities (for Phase 1 / G-leaf):
    #   a_l = |INT rho_pq r_A^l P_l(cos theta_A)| / INT |rho_pq| r_A^l   (isotropic positive leaf: a_0=1, a_l=0)
    print("  Phase-1 leaf candidates: NO pairs ranked by their total weight in T (sum_BCD |T[A,B,C,D]|),")
    print("    with angular content about A, a_l = |INT rho r_A^l P_l(cos th_A)| / INT |rho| r_A^l :")
    from scipy.special import eval_legendre
    rA_ = G['rA']; cA = G['cosA']
    wA = np.einsum('ABCD->A', np.abs(T4))
    ranked = list(np.argsort(wA)[::-1][:8]) + [idx[(0, 1)], idx[(1, 1)], idx[(0, 0)]]
    seen = set()
    for A in ranked:
        if A in seen:
            continue
        seen.add(A)
        p, q = pairs[A]
        rho = rho_no[A]
        al = [abs(G['gint'](rho * rA_ ** l * eval_legendre(l, cA))) / G['gint'](np.abs(rho) * rA_ ** l)
              for l in range(5)]
        print(f"    rho_(NO{p},NO{q}) [{char[p]}x{char[q]}] weight={wA[A]:.4f} gamma_pq={g_tot[p, q]:+.4f}  "
              + " ".join(f"a{l}={al[l]:.3f}" for l in range(5)))

    # ------------------------------------------------------------------ DIAGNOSTIC grid scan
    print("\n  DIAGNOSTIC (scratch GL grids, same Neumann formula; kernels.py untouched, deliverable stays on it):")
    print("    <V_ee> error of the truncated vector vs (NXI, NETA, LMAX, xi_max), + STO self-Coulomb rel errors:")
    sto_ctrl = {lab: (i, STO[lab][0]) for i, lab in enumerate(ref['labels']) if lab in STO}
    scan = [(72, 44, 34, G['params']['xi_max']), (144, 44, 34, G['params']['xi_max']),
            (72, 88, 34, G['params']['xi_max']), (72, 44, 60, G['params']['xi_max']),
            (144, 88, 50, G['params']['xi_max']), (288, 132, 60, G['params']['xi_max']),
            (144, 88, 50, 20.0)]
    scan_rows = []
    if exact:
        print("    (skipped on the exact path: the scan is a diagnostic of the LEGACY formula's grid error)")
        scan = []
    for (nx, ne, lm, xm) in scan:
        t1 = time.time()
        dv, ctrl = vee_on_grid(ref['prims'], T_no, pairs, G_tot, eri_no, R, nx, ne, lm, xm, sto_ctrl)
        scan_rows.append((nx, ne, lm, xm, dv))
        print(f"    NXI={nx:3d} NETA={ne:3d} LMAX={lm:2d} xi_max={xm:5.2f}: dVee={dv:+8.4f} mHa   "
              + "  ".join(f"{k}:{v:+.1e}" for k, v in ctrl.items()) + f"   [{time.time() - t1:.0f}s]", flush=True)

    # ------------------------------------------------------------------ verdict
    d_ok = trace_max < 1e-10 and prod_max < 1e-9
    verdict = "GO" if (g0_pass and d_ok) else "STOP"
    print("\n" + "=" * 96)
    print(f"G0: dT={dT:.4f} mHa  dVne={dV:.4f} mHa  dVee={dVee:.4f} mHa  (gate {GATE_MHA})  -> {'PASS' if g0_pass else 'FAIL'}")
    print(f"(d) bookkeeping: trace max|Delta|={trace_max:.1e}, product-route max rel={prod_max:.1e} -> {'PASS' if d_ok else 'FAIL'}")
    print(f"VERDICT: {verdict}")
    print("=" * 96, flush=True)

    # ------------------------------------------------------------------ save
    g0_table = np.array([[gv, ev, (gv - ev) * 1e3] for (_, gv, ev) in rows])
    save = dict(
        cfg=repr(CFG), R=R, Vnn=ref['Vnn'], Z_A=Z_A, Z_B=Z_B, E_full=E_full, E_full_no=E_full_no,
        E_trunc=E_trunc, E_ref=E_REF, trunc=TRUNC, verdict=verdict, exact_neumann=bool(exact),
        labels=np.array(ref['labels']), prims=np.array(ref['prims'], dtype=object),
        S_prim=ref['S'], h1_prim=ref['h1'], T1_prim=ref['T1'], eri_prim=ref['eri'],
        C=ref['C'], X=ref['X'], Tmap=ref['Tmap'], h1_f=ref['h1_f'], T1_f=ref['T1_f'], eri_f=ref['eri_f'],
        U_no=U, occ_no=occ, T_no=T_no, h1_no=h1_no, T1_no=T1_no, eri_no=eri_no,
        no_char=np.array(char), no_pops=np.array([pops[g] for g in ('core', 'Li2s', 'H', 'bond')]),
        c_full_canonical=c, c_full_no=cn, kept_idx=kept, dets_kept=np.array([dets_n[I] for I in kept]),
        c_kept=ct[kept], CS_t=CS, sigma=sigma, apairs=np.array(apairs), pairs=np.array(pairs),
        gamma_a=ga, gamma_b=gb, Gamma_aa=Gaa, Gamma_bb=Gbb, Gamma_ab=Gab,
        grid_XI=KG.XI, grid_ETA=KG.ETA, grid_WXI=KG.WXI, grid_WETA=KG.WETA, grid_geo=G['geo'],
        grid_rA=G['rA'], grid_rB=G['rB'], grid_pts=G['pts'], grid_params=repr(G['params']),
        prim_val=pval, prim_grad=pgrad, prim_lap=plap,
        no_val=nval, no_grad=ngrad, no_lap=nlap, rho_no=rho_no, V_no=V_no,
        S_grid_no=Sg, T_grid_no=Tg, Vne_grid_no=Vg, J_grid_no=J_no, J4_grid_no=J4,
        T3aab=tr['T3aab'], T3abb=tr['T3abb'], T2aa=tr['T2aa'], T2ab=tr['T2ab'], T2bb=tr['T2bb'],
        T1a=tr['T1a'], T1b=tr['T1b'], delta_pair=d, act_a=act_a, act_b=act_b,
        g0_table=g0_table, g0_rows=np.array([r[0] for r in rows]),
        vee_grid_scan=np.array(scan_rows),   # diagnostic: (NXI, NETA, LMAX, xi_max, dVee_mHa)
    )
    coo_bytes = nnz * (8 + 4 + 4)
    if coo_bytes < 40e6:
        rr, cc = np.nonzero(T)
        save.update(T_storage='coo', T_row=rr.astype(np.int32), T_col=cc.astype(np.int32), T_data=T[rr, cc])
        note = f"T stored as COO ({nnz} nnz, {coo_bytes / 1e6:.0f} MB)"
    else:
        save.update(T_storage='factorized')
        note = (f"T NOT stored explicitly (dense {T.nbytes / 1e6:.0f} MB, COO {coo_bytes / 1e6:.0f} MB, "
                f"debug/data is git-tracked): rebuild in ~1 s with build_T(CS_t, build_Eexp(...), n_ap, npair)")
    os.makedirs(DATA, exist_ok=True)
    out_path = OUT_EXACT if exact else OUT
    np.savez(out_path, **save)
    print(f"saved {out_path} ({os.path.getsize(out_path) / 1e6:.1f} MB);  {note}")
    print(f"wall time {time.time() - T0:.0f} s")
    return verdict


if __name__ == '__main__':
    main(exact=('--exact' in sys.argv))
