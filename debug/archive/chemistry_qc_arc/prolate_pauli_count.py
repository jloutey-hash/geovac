"""First measured Pauli count for a prolate-spheroidal (Level 2) qubit Hamiltonian.

Analysis #5 of debug/sprint_product_basis_tradeoff_scoping_memo.md (the G7 gap).

Systems
  (a) H2, R = 1.4011 bohr, Paper 12 prolate product orbitals as one-electron functions
      (debug/prolate_allelectron_c4.h2_orbitals_c4: sigma ProductFns + pi shells, the
      validated G-H2PI engine), M = 4, 6, 8, 10, 12, 14 spatial orbitals.
  (b) LiH, R = 3.015 bohr, Route C 10-orbital set, integrals from the cache
      debug/data/prolate_separation_constant_stageC_ints.npz (chemist (pq|rs), definite m).

Orbitals: canonical orthogonalization done BLOCK-WISE in signed m (so every orbital keeps a
definite m), then (i) bare one-body MOs (eigenvectors of h1), (ii) RHF MOs, (iii) FCI natural
orbitals (1-RDM diagonalized block-wise in m).  Variant "real-pi": pi pairs rotated to
cos/sin (real) combinations.

Benchmarks: Gaussian bases via an in-file McMurchie-Davidson integral engine (PySCF is not
installable on this Windows/Python 3.14 build; basis data from basis_set_exchange), validated
against openfermion's shipped PySCF-generated MolecularData (H2 STO-3G 0.7414 A, H2 6-31G
0.75 A, LiH STO-3G 1.45 A); GeoVac library Hamiltonians via geovac.ecosystem_export.

Conventions (match the corpus): N_Pauli = non-identity Pauli terms (identity-included also
reported); lambda = sum |c| over non-identity terms (identity-included also reported);
QWC = geovac.measurement_grouping.count_qwc_groups (greedy sorted insertion); ERI density =
fraction of the M^4 MO-basis (pq|rs) with |.| > 1e-12.  JW via openfermion, fermion operator
via geovac.qubit_encoding.build_fermion_op_from_integrals (threshold 1e-12).
Energy check: the qubit operator itself is applied to the (N, S_z=0) computational-basis
sector and diagonalized (never the full 2^Q space).

Usage: python debug/prolate_pauli_count.py [validate|h2|lih|gauss|library|all]
"""
from __future__ import annotations

import json
import os
import sys
import time
from itertools import combinations
from typing import Dict, List, Optional, Tuple

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from scipy.special import gammainc, gamma

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)

from openfermion import jordan_wigner  # noqa: E402
from geovac.qubit_encoding import build_fermion_op_from_integrals  # noqa: E402
from geovac.measurement_grouping import count_qwc_groups  # noqa: E402

DATA = os.path.join(HERE, "data")
OUT_JSON = os.environ.get("PPC_JSON", os.path.join(DATA, "prolate_pauli_count.json"))
ANG = 1.0 / 0.529177210903  # bohr per angstrom


# =====================================================================================
# Qubit-operator metrics + sector diagonalization
# =====================================================================================
def qubit_metrics(h1: np.ndarray, eri: np.ndarray, ecore: float, nelec: int,
                  do_qwc: bool = True, do_energy: bool = True) -> Dict:
    M = h1.shape[0]
    fop = build_fermion_op_from_integrals(h1, eri, ecore)
    qop = jordan_wigner(fop)
    qop.compress(1e-12)
    terms = qop.terms
    ident = complex(terms.get((), 0.0))
    coefs = np.array([abs(c) for k, c in terms.items() if k != ()])
    out = dict(M=M, Q=2 * M, N_pauli=len(coefs), N_pauli_incl_id=len(terms),
               lam=float(coefs.sum()), lam_incl_id=float(coefs.sum() + abs(ident)),
               identity_coeff=float(ident.real),
               eri_density=float(np.mean(np.abs(eri) > 1e-12)))
    if do_qwc:
        out["QWC"] = count_qwc_groups(qop)
    if do_energy:
        out["E_qubit_sector"] = sector_ground_energy(qop, 2 * M, nelec)
    return out


def sector_states(nq: int, nelec: int) -> np.ndarray:
    """Computational-basis states with nelec/2 alpha (even qubits) + nelec/2 beta (odd)."""
    na = nelec // 2
    M = nq // 2
    st = []
    for a in combinations(range(M), na):
        for b in combinations(range(M), na):
            st.append(sum(1 << (2 * p) for p in a) + sum(1 << (2 * p + 1) for p in b))
    return np.array(sorted(st), dtype=np.int64)


def sector_ground_energy(qop, nq: int, nelec: int) -> float:
    """Apply every Pauli string of the JW operator to the (N, Sz=0) sector; lowest eigenvalue."""
    st = sector_states(nq, nelec)
    n = len(st)
    rows, cols, vals = [], [], []
    for term, c in qop.terms.items():
        xm = ym = zm = 0
        for q, P in term:
            if P == "X":
                xm |= 1 << q
            elif P == "Y":
                ym |= 1 << q
            else:
                zm |= 1 << q
        nY = bin(ym).count("1")
        tgt = st ^ (xm | ym)
        par = np.bitwise_count(st & (ym | zm)) & 1
        ph = complex(c) * (1j ** nY) * (1.0 - 2.0 * par)
        pos = np.searchsorted(st, tgt)
        pos = np.minimum(pos, n - 1)
        ok = st[pos] == tgt
        rows.append(pos[ok]); cols.append(np.nonzero(ok)[0]); vals.append(ph[ok])
    r = np.concatenate(rows); c = np.concatenate(cols); v = np.concatenate(vals)
    H = sp.coo_matrix((v, (r, c)), shape=(n, n)).tocsr()
    imag = abs(H.imag).max() if H.nnz else 0.0
    H = H.real
    if n <= 3000:
        e = np.linalg.eigvalsh(H.toarray())[0]
    else:
        e = eigsh(H, k=1, which="SA", tol=1e-12)[0][0]
    assert imag < 1e-9, f"sector Hamiltonian not real: {imag}"
    return float(e)


# =====================================================================================
# Orbital utilities (block-wise in signed m)
# =====================================================================================
def blockwise_orth(S: np.ndarray, m: np.ndarray, tol: float = 1e-10) -> Tuple[np.ndarray, np.ndarray]:
    """Canonical orthogonalization within each signed-m block; returns X (n x k), m of columns."""
    cols, mcols = [], []
    for mv in sorted(set(m.tolist())):
        idx = np.where(m == mv)[0]
        w, U = np.linalg.eigh(S[np.ix_(idx, idx)])
        keep = w > tol * w.max()
        Xb = np.zeros((len(m), keep.sum()))
        Xb[idx] = U[:, keep] / np.sqrt(w[keep])
        cols.append(Xb); mcols += [mv] * int(keep.sum())
    return np.hstack(cols), np.array(mcols)


def blockwise_eigh(A: np.ndarray, m: np.ndarray, B: Optional[np.ndarray] = None,
                   descending: bool = False) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Eigenvectors of A (orthonormal frame) within each m block; global sort by eigenvalue."""
    n = A.shape[0]
    ev, vecs, ms = [], [], []
    for mv in sorted(set(m.tolist())):
        idx = np.where(m == mv)[0]
        w, U = np.linalg.eigh(A[np.ix_(idx, idx)])
        U = U * np.sign(U[np.argmax(np.abs(U), axis=0), np.arange(U.shape[1])])[None, :]
        for k in range(len(idx)):
            v = np.zeros(n); v[idx] = U[:, k]
            ev.append(w[k]); vecs.append(v); ms.append(mv)
    ev = np.array(ev); o = np.argsort(-ev if descending else ev, kind="stable")
    return ev[o], np.array(vecs).T[:, o], np.array(ms)[o]


def transform(h1, eri, C):
    """C may be complex (real-pi variant); chemist (pq|rs) with conj on p, r."""
    Cc = C.conj()
    h = Cc.T @ h1 @ C
    e = np.einsum("ap,bq,cr,ds,abcd->pqrs", Cc, C, Cc, C, eri, optimize=True)
    return h, e


def rhf(h1, eri, nocc, m=None, maxit=300, tol=1e-11):
    """RHF in an orthonormal basis (chemist eri, real). Block-wise in m if m given."""
    n = h1.shape[0]
    mm = np.zeros(n, int) if m is None else m
    D = np.zeros((n, n))
    focks, errs = [], []
    E_old = 0.0
    for it in range(maxit):
        J = np.einsum("pqrs,rs->pq", eri, D)
        K = np.einsum("prqs,rs->pq", eri, D)
        F = h1 + 2 * J - K
        E = float(np.sum(D * (h1 + F)))
        err = F @ D - D @ F
        focks.append(F); errs.append(err)
        focks, errs = focks[-8:], errs[-8:]
        if len(focks) > 2:
            k = len(focks)
            Bm = -np.ones((k + 1, k + 1)); Bm[k, k] = 0
            for i in range(k):
                for j in range(k):
                    Bm[i, j] = np.sum(errs[i] * errs[j])
            rhs = np.zeros(k + 1); rhs[k] = -1
            try:
                cc = np.linalg.solve(Bm, rhs)[:k]
                F = sum(c * f for c, f in zip(cc, focks))
            except np.linalg.LinAlgError:
                pass
        eps, C, ms = blockwise_eigh(F, mm)
        D = C[:, :nocc] @ C[:, :nocc].T
        if it > 1 and abs(E - E_old) < tol and np.abs(err).max() < 1e-8:
            break
        E_old = E
    return E, eps, C, ms


# =====================================================================================
# FCI 1-RDM (spin-summed) from a sector-basis ground state, for natural orbitals
# =====================================================================================
def fci_1rdm(h1, eri, nelec):
    """Exact FCI via the prolate engine's Slater-Condon (for NOs only)."""
    from prolate_allelectron_fci import _dets, _matel
    M = h1.shape[0]
    dets = _dets(M, nelec // 2, nelec // 2)
    nd = len(dets)
    H = np.zeros((nd, nd))
    for a in range(nd):
        for b in range(a, nd):
            H[a, b] = H[b, a] = _matel(h1, eri, dets[a], dets[b])
    w, V = np.linalg.eigh(H)
    v = V[:, 0]
    idx = {D: k for k, D in enumerate(dets)}
    gam = np.zeros((M, M))
    for a, Da in enumerate(dets):
        if abs(v[a]) < 1e-14:
            continue
        for so in Da:
            for t in range(2 * M):
                if t % 2 != so % 2 or (t != so and t in Da):
                    continue
                Db = tuple(sorted([x for x in Da if x != so] + [t]))
                b = idx.get(Db)
                if b is None:
                    continue
                sgn = (-1) ** (Da.index(so) + Db.index(t))
                gam[t // 2, so // 2] += sgn * v[a] * v[b]
    return float(w[0]), 0.5 * (gam + gam.T)


def orbital_variants(S, h1, eri, m, nelec, Vnn, M_keep=None, real_pi=True):
    """Return {label: (h_mo, eri_mo, m_labels)} for bare-h1 / RHF / NO (+ real-pi)."""
    X, mx = blockwise_orth(S, m)
    ho, eo = transform(h1, eri, X)
    out = {}
    eps, C, ms = blockwise_eigh(ho, mx)
    out["bare-h1 MOs"] = (C, ms)
    E_hf, _, Chf, mhf = rhf(ho, eo, nelec // 2, mx)
    out["RHF MOs"] = (Chf, mhf)
    hb, ebb = transform(ho, eo, C)
    E_fci, gam = fci_1rdm(hb, ebb, nelec)
    occ, U, mno = blockwise_eigh(gam, ms, descending=True)
    out["FCI NOs"] = (C @ U, mno)
    res = {}
    for lab, (Cm, mm) in out.items():
        k = Cm.shape[1] if M_keep is None else M_keep
        h, e = transform(ho, eo, Cm[:, :k])
        res[lab] = (h, e, mm[:k])
    if real_pi and any(abs(x) > 0 for x in ms):
        h, e, mm = res["bare-h1 MOs"]
        res["bare-h1 MOs, real pi"] = realify(h, e, mm)
    return res, E_fci + Vnn, E_hf + Vnn


def realify(h, e, mm):
    """Rotate each (m=+mu, m=-mu) degenerate MO pair to cos/sin combinations."""
    n = len(mm)
    U = np.eye(n, dtype=complex)
    used = set()
    for i in range(n):
        if mm[i] > 0 and i not in used:
            j = next(k for k in range(n) if mm[k] == -mm[i] and k not in used
                     and abs(h[k, k] - h[i, i]) < 1e-8)
            used |= {i, j}
            U[:, i] = 0; U[:, j] = 0
            U[i, i] = U[j, i] = 1 / np.sqrt(2)
            U[i, j] = -1j / np.sqrt(2); U[j, j] = 1j / np.sqrt(2)
    hh, ee = transform(h.astype(complex), e.astype(complex), U)
    assert np.abs(hh.imag).max() < 1e-10 and np.abs(ee.imag).max() < 1e-10
    return hh.real, ee.real, mm


# =====================================================================================
# Gaussian integrals: McMurchie-Davidson, vectorized over primitive pairs
# =====================================================================================
def boys(n, T):
    T = np.asarray(T, float)
    small = T < 1e-12
    Ts = np.where(small, 1.0, T)
    val = gamma(n + 0.5) * gammainc(n + 0.5, Ts) / (2 * Ts ** (n + 0.5))
    return np.where(small, 1.0 / (2 * n + 1) - T / (2 * n + 3), val)


def Ecoef(i, j, t, Qx, a, b):
    p = a + b; q = a * b / p
    if t < 0 or t > i + j:
        return np.zeros_like(a)
    if i == j == t == 0:
        return np.exp(-q * Qx * Qx)
    if j == 0:
        return (1 / (2 * p)) * Ecoef(i - 1, j, t - 1, Qx, a, b) - (q * Qx / a) * Ecoef(i - 1, j, t, Qx, a, b) \
            + (t + 1) * Ecoef(i - 1, j, t + 1, Qx, a, b)
    return (1 / (2 * p)) * Ecoef(i, j - 1, t - 1, Qx, a, b) + (q * Qx / b) * Ecoef(i, j - 1, t, Qx, a, b) \
        + (t + 1) * Ecoef(i, j - 1, t + 1, Qx, a, b)


def Rtable(L, p, X, Y, Z):
    """R^0_{tuv} for t+u+v <= L (arrays broadcast over p, X, Y, Z)."""
    RPC2 = X * X + Y * Y + Z * Z
    T = p * RPC2
    memo = {}

    def R(t, u, v, n):
        key = (t, u, v, n)
        if key in memo:
            return memo[key]
        if t == u == v == 0:
            val = (-2 * p) ** n * boys(n, T)
        elif t == u == 0:
            val = Z * R(t, u, v - 1, n + 1)
            if v > 1:
                val = val + (v - 1) * R(t, u, v - 2, n + 1)
        elif t == 0:
            val = Y * R(t, u - 1, v, n + 1)
            if u > 1:
                val = val + (u - 1) * R(t, u - 2, v, n + 1)
        else:
            val = X * R(t - 1, u, v, n + 1)
            if t > 1:
                val = val + (t - 1) * R(t - 2, u, v, n + 1)
        memo[key] = val
        return val

    shape = np.broadcast(p, X).shape
    out = np.zeros((L + 1, L + 1, L + 1) + shape)
    for t in range(L + 1):
        for u in range(L + 1 - t):
            for v in range(L + 1 - t - u):
                out[t, u, v] = R(t, u, v, 0)
    return out


def cart_comps(l):
    return [(lx, ly, l - lx - ly) for lx in range(l, -1, -1) for ly in range(l - lx, -1, -1)]


def dfact(n):
    return 1.0 if n <= 0 else n * dfact(n - 2)


def load_basis(name, atoms):
    """atoms: list of (Z, xyz bohr). Returns shells (center, l, exps, coefs incl. prim norm)
    and the cart->pure transform (d -> 5 real solid harmonics; normalization done later)."""
    import basis_set_exchange as bse
    Zs = sorted(set(z for z, _ in atoms))
    bd = bse.get_basis(name, elements=Zs, fmt=None, uncontract_general=True, uncontract_spdf=True)
    shells = []
    for Z, xyz in atoms:
        for sh in bd["elements"][str(Z)]["electron_shells"]:
            l = sh["angular_momentum"][0]
            ex = np.array([float(x) for x in sh["exponents"]])
            for cs in sh["coefficients"]:
                c = np.array([float(x) for x in cs])
                nz = c != 0
                N = (2 * ex / np.pi) ** 0.75 * (4 * ex) ** (l / 2) / np.sqrt(dfact(2 * l - 1))
                shells.append(dict(A=np.array(xyz, float), l=l, a=ex[nz], c=(c * N)[nz]))
    ncart = sum(len(cart_comps(s["l"])) for s in shells)
    cols, par = [], []
    off = 0
    for s in shells:
        comps = cart_comps(s["l"])
        if s["l"] < 2:
            for k in range(len(comps)):
                v = np.zeros(ncart); v[off + k] = 1; cols.append(v)
                par.append(-1 if comps[k][0] % 2 else 1)
        else:
            assert s["l"] == 2
            ix = {c: off + k for k, c in enumerate(comps)}
            for combo in ({(1, 1, 0): 1}, {(1, 0, 1): 1}, {(0, 1, 1): 1},
                          {(2, 0, 0): 1, (0, 2, 0): -1},
                          {(0, 0, 2): 2, (2, 0, 0): -1, (0, 2, 0): -1}):
                v = np.zeros(ncart)
                for cmp, w in combo.items():
                    v[ix[cmp]] = w
                cols.append(v)
                par.append(-1 if next(iter(combo))[0] % 2 else 1)
        off += len(comps)
    return shells, np.array(cols).T, np.array(par)


def gaussian_integrals(shells, Tcp, nuclei):
    offs = np.cumsum([0] + [len(cart_comps(s["l"])) for s in shells])
    nc = offs[-1]
    S = np.zeros((nc, nc)); T = np.zeros((nc, nc)); V = np.zeros((nc, nc))
    pairs = {}
    ns = len(shells)
    for I in range(ns):
        for J in range(ns):
            A, B = shells[I], shells[J]
            a = A["a"][:, None] * np.ones(len(B["a"]))[None, :]
            b = np.ones(len(A["a"]))[:, None] * B["a"][None, :]
            a, b = a.ravel(), b.ravel()
            cc = (A["c"][:, None] * B["c"][None, :]).ravel()
            p = a + b
            P = (a[:, None] * A["A"] + b[:, None] * B["A"]) / p[:, None]
            Q = A["A"] - B["A"]
            lab = A["l"] + B["l"]
            Eab = {}
            for ia, ca in enumerate(cart_comps(A["l"])):
                for ib, cb in enumerate(cart_comps(B["l"])):
                    E = np.zeros((lab + 1, lab + 1, lab + 1, len(p)))
                    Ex = [Ecoef(ca[0], cb[0], t, Q[0], a, b) for t in range(lab + 1)]
                    Ey = [Ecoef(ca[1], cb[1], t, Q[1], a, b) for t in range(lab + 1)]
                    Ez = [Ecoef(ca[2], cb[2], t, Q[2], a, b) for t in range(lab + 1)]
                    for t in range(lab + 1):
                        for u in range(lab + 1):
                            for v in range(lab + 1):
                                E[t, u, v] = Ex[t] * Ey[u] * Ez[v] * cc
                    Eab[(ia, ib)] = E
                    # overlap
                    S[offs[I] + ia, offs[J] + ib] = np.sum(E[0, 0, 0] * (np.pi / p) ** 1.5)

                    def ov(la, lb):
                        if min(lb) < 0:
                            return np.zeros_like(a)
                        return (Ecoef(la[0], lb[0], 0, Q[0], a, b) * Ecoef(la[1], lb[1], 0, Q[1], a, b)
                                * Ecoef(la[2], lb[2], 0, Q[2], a, b) * (np.pi / p) ** 1.5)
                    l2, m2, n2 = cb
                    kin = b * (2 * (l2 + m2 + n2) + 3) * ov(ca, cb)
                    kin -= 2 * b * b * (ov(ca, (l2 + 2, m2, n2)) + ov(ca, (l2, m2 + 2, n2)) + ov(ca, (l2, m2, n2 + 2)))
                    kin -= 0.5 * (l2 * (l2 - 1) * ov(ca, (l2 - 2, m2, n2)) + m2 * (m2 - 1) * ov(ca, (l2, m2 - 2, n2))
                                  + n2 * (n2 - 1) * ov(ca, (l2, m2, n2 - 2)))
                    T[offs[I] + ia, offs[J] + ib] = np.sum(kin * cc)
            vv = np.zeros((len(cart_comps(A["l"])), len(cart_comps(B["l"]))))
            for Zc, C in nuclei:
                PC = P - np.array(C)
                Rt = Rtable(lab, p, PC[:, 0], PC[:, 1], PC[:, 2])
                for (ia, ib), E in Eab.items():
                    vv[ia, ib] += -Zc * np.sum(2 * np.pi / p * np.einsum("tuvi,tuvi->i", E, Rt))
            V[offs[I]:offs[I + 1], offs[J]:offs[J + 1]] = vv
            pairs[(I, J)] = (p, P, Eab, lab)
    # ERIs over shell quartets (8-fold symmetry)
    G = np.zeros((nc, nc, nc, nc))
    pl = [(I, J) for I in range(ns) for J in range(I + 1)]
    for x, (I, J) in enumerate(pl):
        p, P, Eab, lab = pairs[(I, J)]
        for (K, Lq) in pl[:x + 1]:
            q, Qc, Ecd, lcd = pairs[(K, Lq)]
            L = lab + lcd
            pp = p[:, None]; qq = q[None, :]
            al = pp * qq / (pp + qq)
            PQ = P[:, None, :] - Qc[None, :, :]
            Rt = Rtable(L, al, PQ[..., 0], PQ[..., 1], PQ[..., 2])
            tt = np.arange(lab + 1); TT = np.arange(lcd + 1)
            i6 = (tt[:, None, None, None, None, None] + TT[None, None, None, :, None, None])
            j6 = (tt[None, :, None, None, None, None] + TT[None, None, None, None, :, None])
            k6 = (tt[None, None, :, None, None, None] + TT[None, None, None, None, None, :])
            R6 = Rt[np.minimum(i6, L), np.minimum(j6, L), np.minimum(k6, L)]
            R6 = R6 * ((i6 + j6 + k6) <= L)[..., None, None]
            pref = 2 * np.pi ** 2.5 / (pp * qq * np.sqrt(pp + qq))
            sgn = (-1.0) ** (TT[:, None, None] + TT[None, :, None] + TT[None, None, :])
            R6p = R6 * pref
            for (ia, ib), E1 in Eab.items():
                tmp = np.einsum("tuvi,tuvTUVij->TUVj", E1, R6p, optimize=True)
                for (ic, id_), E2 in Ecd.items():
                    val = float(np.sum(tmp * sgn[..., None] * E2))
                    i1, i2 = offs[I] + ia, offs[J] + ib
                    i3, i4 = offs[K] + ic, offs[Lq] + id_
                    for (w, x_, y, z) in ((i1, i2, i3, i4), (i2, i1, i3, i4), (i1, i2, i4, i3), (i2, i1, i4, i3),
                                          (i3, i4, i1, i2), (i4, i3, i1, i2), (i3, i4, i2, i1), (i4, i3, i2, i1)):
                        G[w, x_, y, z] = val
    # cart -> pure, then normalize
    Sp = Tcp.T @ S @ Tcp
    nrm = 1 / np.sqrt(np.diag(Sp))
    Tn = Tcp * nrm[None, :]
    Sp = Tn.T @ S @ Tn
    h = Tn.T @ (T + V) @ Tn
    G = np.einsum("ap,bq,cr,ds,abcd->pqrs", Tn, Tn, Tn, Tn, G, optimize=True)
    return Sp, h, G


def gaussian_system(basis, atoms):
    shells, Tcp, par = load_basis(basis, atoms)
    nuclei = [(float(Z), xyz) for Z, xyz in atoms]
    S, h, G = gaussian_integrals(shells, Tcp, nuclei)
    Vnn = sum(Za * Zb / np.linalg.norm(np.array(A) - np.array(B))
              for i, (Za, A) in enumerate(atoms) for (Zb, B) in atoms[:i])
    nelec = sum(Z for Z, _ in atoms)
    X, _ = blockwise_orth(S, np.zeros(S.shape[0], int), tol=1e-12)
    ho, eo = transform(h, G, X)
    E_hf, eps, C, _ = rhf(ho, eo, nelec // 2)
    # symmetry-clean degenerate (pi / delta) MO pairs to definite x-reflection parity
    # (px / py-like real orbitals; otherwise the degenerate pair is an arbitrary rotation)
    Pmat = X.T @ S @ np.diag(par) @ X   # reflection x -> -x in the orthonormal frame
    k = 0
    while k < len(eps):
        g = [k]
        while g[-1] + 1 < len(eps) and abs(eps[g[-1] + 1] - eps[k]) < 1e-6:
            g.append(g[-1] + 1)
        if len(g) > 1:
            Cg = C[:, g]
            w, U = np.linalg.eigh(Cg.T @ Pmat @ Cg)
            C[:, g] = Cg @ U
        k = g[-1] + 1
    mopar = np.einsum("ip,ij,jp->p", C, Pmat, C)
    return dict(h=ho, eri=eo, C=C, eps=eps, mopar=mopar, E_hf=E_hf + Vnn, Vnn=float(Vnn),
                nelec=nelec, nbf=S.shape[0])


def complexify_gauss(g, k):
    """MO integrals (lowest k RHF MOs) with each degenerate (px-like, py-like) pair rotated
    to definite-m combinations (px +/- i py)/sqrt(2) -- the prolate/composed convention."""
    C = g["C"][:, :k].astype(complex)
    eps, par = g["eps"][:k], g["mopar"][:k]
    i = 0
    while i < k - 1:
        if abs(eps[i + 1] - eps[i]) < 1e-6:
            a, b = (i, i + 1) if par[i] < 0 else (i + 1, i)   # a = x-odd (px-like)
            ca, cb = C[:, a].copy(), C[:, b].copy()
            C[:, a] = (ca + 1j * cb) / np.sqrt(2)
            C[:, b] = (ca - 1j * cb) / np.sqrt(2)
            i += 2
        else:
            i += 1
    h, e = transform(g["h"].astype(complex), g["eri"].astype(complex), C)
    assert np.abs(h.imag).max() < 1e-10 and np.abs(e.imag).max() < 1e-10, "complex integrals"
    return h.real, e.real


def run_gauss_complex():
    R_h2 = 1.4011
    g = gaussian_system("cc-pvdz", [(1, (0, 0, -R_h2 / 2)), (1, (0, 0, R_h2 / 2))])
    for k in (6, 8, 10):
        h, e = complexify_gauss(g, k)
        q = qubit_metrics(h, e, g["Vnn"], 2)
        record("gauss_complex_m", system="H2", basis="cc-pvdz", orbitals=f"RHF MOs lowest {k}, complex-m pi", **q)
    g = gaussian_system("6-31g", [(3, (0, 0, 0)), (1, (0, 0, 3.015))])
    h, e = complexify_gauss(g, 10)
    q = qubit_metrics(h, e, g["Vnn"], 4)
    record("gauss_complex_m", system="LiH", basis="6-31g", orbitals="RHF MOs lowest 10, complex-m pi", **q)


# =====================================================================================
# Drivers
# =====================================================================================
RESULTS: Dict[str, List[Dict]] = {}


def record(section, **row):
    RESULTS.setdefault(section, []).append(row)
    keys = ["system", "basis", "orbitals", "Q", "N_pauli", "N_pauli_incl_id", "lam", "lam_incl_id",
            "QWC", "eri_density", "E_qubit_sector", "E_ref"]
    print("  " + " | ".join(f"{k}={row[k]:.6f}" if isinstance(row.get(k), float) else f"{k}={row.get(k)}"
                            for k in keys if k in row), flush=True)
    save()


def save():
    os.makedirs(DATA, exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(RESULTS, f, indent=1)


def run_validate():
    """Engine vs openfermion's PySCF-generated MolecularData (same geometry, same basis)."""
    from openfermion import MolecularData
    import openfermion
    dd = os.path.join(os.path.dirname(openfermion.__file__), "testing", "data")
    cases = [("sto-3g", "H2_sto-3g_singlet_0.7414", [(1, (0, 0, 0)), (1, (0, 0, 0.7414 * ANG))]),
             ("6-31g", "H2_6-31g_singlet_0.75", [(1, (0, 0, 0)), (1, (0, 0, 0.75 * ANG))]),
             ("sto-3g", "H1-Li1_sto-3g_singlet_1.45", [(3, (0, 0, 0)), (1, (0, 0, 1.45 * ANG))])]
    for basis, fn, atoms in cases:
        md = MolecularData(filename=os.path.join(dd, fn))
        md.load()
        g = gaussian_system(basis, atoms)
        h, e = transform(g["h"], g["eri"], g["C"])
        m = qubit_metrics(h, e, g["Vnn"], g["nelec"])
        # openfermion's own integrals (two_body[p,q,r,s] = (ps|qr) chemist)
        of_eri = np.einsum("pqrs->psqr", md.two_body_integrals)
        m2 = qubit_metrics(md.one_body_integrals, of_eri, float(md.nuclear_repulsion), g["nelec"])
        record("validate", system=fn, basis=basis, Q=m["Q"], E_hf_engine=g["E_hf"], E_hf_pyscf=float(md.hf_energy),
               E_qubit_sector=m["E_qubit_sector"], E_ref=float(md.fci_energy),
               N_pauli=m["N_pauli"], N_pauli_pyscf_ints=m2["N_pauli"], E_pyscf_ints=m2["E_qubit_sector"],
               lam=m["lam"], lam_pyscf_ints=m2["lam"], QWC=m["QWC"], QWC_pyscf_ints=m2["QWC"])


def run_h2():
    import mpmath as mp
    import prolate_allelectron_c4 as c4
    R = 1.4011
    for (Ms, npi) in ((4, 0), (6, 0), (6, 1), (6, 2), (6, 3), (6, 4)):
        orbs = c4.h2_orbitals_c4(M_sigma=Ms, n_pi=npi)
        S, h1 = c4.one_body_general(orbs, R, 1.0, 1.0)
        S = np.array(S, float); h1 = np.array(h1, float)
        eri = c4.build_eri_tensor_m(orbs, R)
        m = np.array([o.msign for o in orbs])
        np.savez(os.path.join(DATA, f"prolate_pauli_h2_M{len(orbs)}_ints.npz"), S=S, h1=h1, eri=eri, m=m)
        var, E_fci, E_hf = orbital_variants(S, h1, eri, m, 2, 1.0 / R)
        for lab, (h, e, mm) in var.items():
            q = qubit_metrics(h, e, 1.0 / R, 2)
            record("h2_prolate", system="H2", basis=f"prolate sigma{Ms}+{npi}pi", orbitals=lab,
                   E_ref=E_fci, E_hf=E_hf, **q)


def run_lih():
    d = np.load(os.path.join(DATA, "prolate_separation_constant_stageC_ints.npz"))
    S, h1, eri, m = d["S"], d["h1"], d["eri"], d["msign"]
    R = 3.015
    Vnn = 3.0 / R
    var, E_fci, E_hf = orbital_variants(S, h1, eri, m, 4, Vnn)
    print(f"  LiH Route C: E_FCI = {E_fci:.5f} (logged -8.00794), E_RHF = {E_hf:.5f}", flush=True)
    for lab, (h, e, mm) in var.items():
        q = qubit_metrics(h, e, Vnn, 4)
        record("lih_prolate", system="LiH", basis="prolate Route C (10)", orbitals=lab,
               E_ref=E_fci, E_hf=E_hf, **q)


def run_gauss():
    R_h2 = 1.4011
    h2 = [(1, (0, 0, -R_h2 / 2)), (1, (0, 0, R_h2 / 2))]
    for basis, keeps in (("sto-3g", [None]), ("6-31g", [None]), ("6-31g**", [None, 6, 8]),
                         ("cc-pvdz", [None, 4, 6, 8])):
        g = gaussian_system(basis, h2)
        for k in keeps:
            k2 = g["nbf"] if k is None else k
            h, e = transform(g["h"], g["eri"], g["C"][:, :k2])
            q = qubit_metrics(h, e, g["Vnn"], 2)
            record("h2_gauss", system="H2", basis=basis, orbitals="RHF MOs" + ("" if k is None else f", lowest {k}"),
                   E_hf=g["E_hf"], **q)
    lih = [(3, (0, 0, 0)), (1, (0, 0, 3.015))]
    for basis, keeps in (("sto-3g", [None]), ("6-31g", [None, 10]), ("cc-pvdz", [10, 15, None])):
        t = time.time()
        g = gaussian_system(basis, lih)
        print(f"  LiH {basis}: nbf={g['nbf']} E_RHF={g['E_hf']:.6f} [{time.time()-t:.0f}s]", flush=True)
        for k in keeps:
            k2 = g["nbf"] if k is None else k
            h, e = transform(g["h"], g["eri"], g["C"][:, :k2])
            big = k2 > 12
            q = qubit_metrics(h, e, g["Vnn"], 4, do_qwc=not (k2 > 15), do_energy=k2 <= 15)
            record("lih_gauss", system="LiH", basis=basis, orbitals="RHF MOs" + ("" if k is None else f", lowest {k}"),
                   E_hf=g["E_hf"], **q)


def run_library():
    from geovac.ecosystem_export import hamiltonian
    for name, kw, nelec in (("H2", dict(max_n=2, R=1.4011), 2), ("H2", dict(max_n=3, R=1.4011), 2),
                            ("LiH", dict(max_n=2), 4)):
        H = hamiltonian(name, **kw)
        qop = H.to_openfermion()
        terms = qop.terms
        ident = complex(terms.get((), 0.0))
        coefs = np.array([abs(c) for k, c in terms.items() if k != ()])
        nq = H.n_qubits
        E = sector_ground_energy(qop, nq, nelec) if nq <= 30 else None
        eri = H.eri
        dens = float(np.mean(np.abs(eri) > 1e-12)) if eri is not None else None
        record("library", system=name, basis=f"GeoVac composed max_n={kw['max_n']}", Q=nq,
               N_pauli=len(coefs), N_pauli_incl_id=len(terms), lam=float(coefs.sum()),
               lam_incl_id=float(coefs.sum() + abs(ident)), QWC=count_qwc_groups(qop),
               eri_density=dens, E_qubit_sector=E)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as f:
            RESULTS.update(json.load(f))
    for key, fn in (("validate", run_validate), ("h2", run_h2), ("lih", run_lih),
                    ("gauss", run_gauss), ("gaussc", run_gauss_complex), ("library", run_library)):
        if what in (key, "all"):
            sec = {"validate": ["validate"], "h2": ["h2_prolate"], "lih": ["lih_prolate"],
                   "gauss": ["h2_gauss", "lih_gauss"], "gaussc": ["gauss_complex_m"],
                   "library": ["library"]}[key]
            for s in sec:
                RESULTS.pop(s, None)
            print(f"== {key} ==", flush=True)
            t0 = time.time()
            fn()
            print(f"== {key} done [{time.time()-t0:.0f}s]", flush=True)
