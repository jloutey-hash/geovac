"""
xTC PoC on Li (1s^2 2s, doublet) in the GeoVac Coulomb-Sturmian basis.
=====================================================================

Question (Track 1 + Track 2 follow-on): does xTC -- the three-body TC operator
L3 contracted to an EFFECTIVE two-body operator via the reference 1-RDM
(Christlmaier-Kats-Alavi, JCP 159 014113 (2023)) -- recover chemical accuracy
at a SMALLER / SPARSER two-body qubit Hamiltonian than plain FCI?

Li is the minimal system where the genuine 3-body TC operator is live: L3 needs
THREE DISTINCT electrons (He/H2 have only the j=k two-body piece).

Transcorrelation (Jastrow tau = sum_{i<j} u(r_ij), Slater geminal
u(r) = -(1/2g) e^{-g r}, u'(0)=1/2 singlet Kato cusp):
   Htilde = e^{-tau} H e^{tau} = H + D + K + L3
   D  = -sum_i [ (1/2) lap_i tau + (1/2)(grad_i tau)^2 (j=k part) ]  (2-body, Hermitian)
   K  = -sum_i (grad_i tau).grad_i                                    (2-body, non-Herm)
   L3 = -(1/2) sum_i sum_{j!=i,k!=i,j!=k} grad_i u(r_ij).grad_i u(r_ik)  (3-body)

For s-orbitals grad_i u(r_ij) = u'(r_ij) rhat_ij, so
   L3 kernel  w3(1;2,3) = u'(r_12) u'(r_13) (rhat_12 . rhat_13),
whose angular average FACTORIZES (validated by MC) into a shared-vertex product
   A(a,b,c) = g(a,b) g(a,c),  g(a,b) = (1/2) int_{-1}^1 u'(r_ab) (a - b x)/r_ab dx.

Everything is built on ONE radial grid with the Coulomb-Sturmian s-radials
R_{n0}(r; nk) = N e^{-kr} L_{n-1}^1(2kr) (shared decay k), Loewdin-orthonormalized,
then FCI via convention-free second-quantized operator application (particle-number
projected). Non-Hermitian pieces -> scipy.linalg.eig, real ground state.

Runs a stack of validation gates first (grid integrals vs analytic; geminal->0 must
reproduce plain FCI; xTC 0-body and occ->virt block must match the exact 3-body).
"""
import json, time, itertools
import numpy as np
from scipy.special import eval_genlaguerre
from scipy.linalg import eigh, eig
from itertools import combinations, permutations

# ----------------------------------------------------------------------------
# Coulomb-Sturmian s-radial R_{n0}(r; nk) and derivative  (shared decay k)
#   identical to ctf12 R_and_dR: N=2 k^{3/2}/n, x=2kr, R=N e^{-kr} L_{n-1}^1(x)
# ----------------------------------------------------------------------------
def R_and_dR(n, r, k):
    N = 2.0 * k ** 1.5 / n
    x = 2 * k * r
    L1 = eval_genlaguerre(n - 1, 1, x)
    R = N * np.exp(-k * r) * L1
    L2 = eval_genlaguerre(n - 2, 2, x) if n >= 2 else np.zeros_like(r)
    dR = N * np.exp(-k * r) * (-k * L1 - 2 * k * L2)   # dL^1_{n-1}/dx = -L^2_{n-2}
    return R, dR


def make_grid(k, Ng=1400, r_max=None):
    if r_max is None:
        r_max = 46.0 / k
    t = np.linspace(0.0, 1.0, Ng)
    r = r_max * t ** 2
    r[0] = 1e-9
    wr = np.zeros(Ng)
    wr[1:-1] = (r[2:] - r[:-2]) / 2.0
    wr[0] = (r[1] - r[0]) / 2.0
    wr[-1] = (r[-1] - r[-2]) / 2.0
    return r, wr


# ----------------------------------------------------------------------------
# Geminal
# ----------------------------------------------------------------------------
def up_of(r, g):     # u'(r) = (1/2) e^{-g r}
    return 0.5 * np.exp(-g * r)

def w_kernel(r, g):  # Hermitian TC effective e-e:  1/r + D  (finite at 0)
    r = np.asarray(r, float)
    small = r < 1e-8
    rr = np.where(small, 1.0, r)
    out = -np.expm1(-g * rr) / rr + (g / 2) * np.exp(-g * rr) - 0.25 * np.exp(-2 * g * rr)
    out = np.where(small, 1.5 * g - 0.25, out)
    return out


# ----------------------------------------------------------------------------
# (Ng x Ng) angular (L=0) / vertex kernels via Gauss-Legendre in x=cos(theta12)
# ----------------------------------------------------------------------------
def build_coul_kernel(r, nx=200):
    """gamma-independent L=0 Coulomb kernel only (for fast plain-FCI k-scans)."""
    R1, R2 = np.meshgrid(r, r, indexing='ij')
    xs, ws = np.polynomial.legendre.leggauss(nx)
    Kcoul = np.zeros_like(R1)
    for x, wx in zip(xs, ws):
        r12 = np.sqrt(np.maximum(R1 * R1 + R2 * R2 - 2 * R1 * R2 * x, 1e-30))
        Kcoul += 0.5 * wx * (1.0 / r12)
    return Kcoul


def build_kernels(r, g, nx=200):
    R1, R2 = np.meshgrid(r, r, indexing='ij')
    xs, ws = np.polynomial.legendre.leggauss(nx)
    Kcoul = np.zeros_like(R1); Kw = np.zeros_like(R1)
    KkA = np.zeros_like(R1); KkB = np.zeros_like(R1)
    for x, wx in zip(xs, ws):
        r12 = np.sqrt(np.maximum(R1 * R1 + R2 * R2 - 2 * R1 * R2 * x, 1e-30))
        half = 0.5 * wx
        up = 0.5 * np.exp(-g * r12)
        Kcoul += half * (1.0 / r12)
        Kw += half * w_kernel(r12, g)
        KkA += half * up * (R1 - R2 * x) / r12      # e1 gradient / L3 vertex kernel g(a,b)
        KkB += half * up * (R1 * x - R2) / r12       # e2 gradient
    return dict(coul=Kcoul, w=Kw, kA=KkA, kB=KkB)


# ----------------------------------------------------------------------------
# One-body (grid), spatial, non-orthogonal Sturmian basis, s-only
# ----------------------------------------------------------------------------
def build_one_body(ns, r, wr, k, Z):
    Rtab = {n: R_and_dR(n, r, k) for n in range(1, ns + 1)}
    W = r * r * wr
    S = np.zeros((ns, ns)); h1 = np.zeros((ns, ns))
    for i in range(ns):
        Ri, dRi = Rtab[i + 1]
        for j in range(ns):
            Rj, dRj = Rtab[j + 1]
            S[i, j] = np.sum(Ri * Rj * W)
            T = 0.5 * np.sum(dRi * dRj * W)             # l=0 gradient form
            Vnuc = -Z * np.sum(Ri * Rj * r * wr)         # -Z/r : r^2 dr * (1/r) = r dr
            h1[i, j] = T + Vnuc
    return S, h1, Rtab, W


# ----------------------------------------------------------------------------
# Two-body spatial integrals <ij|op|kl> (electron1: i,k ; electron2: j,l), s-only
# ----------------------------------------------------------------------------
def two_body(ns, Rtab, W, Kmats):
    RR = {i: Rtab[i + 1][0] for i in range(ns)}
    dRR = {i: Rtab[i + 1][1] for i in range(ns)}
    D = {(i, k): RR[i] * RR[k] * W for i in range(ns) for k in range(ns)}      # density*W
    Dd = {(i, k): RR[i] * dRR[k] * W for i in range(ns) for k in range(ns)}    # R_i dR_k *W
    def eri(K):
        out = np.zeros((ns, ns, ns, ns))
        for i in range(ns):
            for j in range(ns):
                for kk in range(ns):
                    for l in range(ns):
                        out[i, j, kk, l] = D[(i, kk)] @ K @ D[(j, l)]
        return out
    eri_coul = eri(Kmats['coul'])
    eri_w = eri(Kmats['w'])
    # convective K2 = -<ij| u'(r) rhat.(grad1-grad2) |kl>
    #  = -(R_i R_k' * W) kA (R_j R_l * W)  + (R_i R_k * W) kB (R_j R_l' * W)
    eri_K = np.zeros((ns, ns, ns, ns))
    for i in range(ns):
        for j in range(ns):
            for kk in range(ns):
                for l in range(ns):
                    t1 = Dd[(i, kk)] @ Kmats['kA'] @ D[(j, l)]
                    t2 = D[(i, kk)] @ Kmats['kB'] @ Dd[(j, l)]
                    eri_K[i, j, kk, l] = -t1 + t2
    return eri_coul, eri_w, eri_K


# ----------------------------------------------------------------------------
# Three-body spatial integral V3[i,j,k; l,m,n]  (vertex particle1: i<->l)
#   = sum_a (R_i R_l W)_a  G_jm(a) G_kn(a),   G_jm = kA @ (R_j R_m W)
# ----------------------------------------------------------------------------
def three_body(ns, Rtab, W, Kmats):
    RR = {i: Rtab[i + 1][0] for i in range(ns)}
    kA = Kmats['kA']
    Dvert = {(i, l): RR[i] * RR[l] * W for i in range(ns) for l in range(ns)}
    G = {(j, m): kA @ (RR[j] * RR[m] * W) for j in range(ns) for m in range(ns)}  # Ng-vector
    V3 = np.zeros((ns, ns, ns, ns, ns, ns))
    for i in range(ns):
        for l in range(ns):
            dv = Dvert[(i, l)]
            for j in range(ns):
                for m in range(ns):
                    Gjm = G[(j, m)]
                    dvG = dv * Gjm
                    for kk in range(ns):
                        for n in range(ns):
                            V3[i, j, kk, l, m, n] = dvG @ G[(kk, n)]
    return V3


# ----------------------------------------------------------------------------
# Loewdin orthonormalization
# ----------------------------------------------------------------------------
def lowdin(S):
    ev, U = eigh(S)
    ev = np.maximum(ev, 1e-14)
    return U @ np.diag(1.0 / np.sqrt(ev)) @ U.T   # X = S^{-1/2}

def transform_1(M, X):  return X @ M @ X
def transform_2(E, X):
    t = np.einsum('pi,ijkl->pjkl', X, E)
    t = np.einsum('qj,pjkl->pqkl', X, t)
    t = np.einsum('rk,pqkl->pqrl', X, t)
    return np.einsum('sl,pqrl->pqrs', X, t)
def transform_3(V, X):
    t = np.einsum('ai,ijklmn->ajklmn', X, V)
    t = np.einsum('bj,ajklmn->abklmn', X, t)
    t = np.einsum('ck,abklmn->abclmn', X, t)
    t = np.einsum('dl,abclmn->abcdmn', X, t)
    t = np.einsum('em,abcdmn->abcden', X, t)
    return np.einsum('fn,abcden->abcdef', X, t)


# ----------------------------------------------------------------------------
# Spin-orbital plumbing.  spin-orbital p = 2*spatial + spin
# ----------------------------------------------------------------------------
def spatial(p): return p >> 1
def spin(p):    return p & 1

def apply_ops(det, ops):
    """ops: list of ('c'|'a', spin_orbital) written LEFT->RIGHT; rightmost acts first."""
    sign = 1; d = det
    for kind, idx in reversed(ops):
        if kind == 'a':
            if idx not in d: return None
            sign *= (-1) ** sum(1 for y in d if y < idx)
            d = tuple(y for y in d if y != idx)
        else:
            if idx in d: return None
            sign *= (-1) ** sum(1 for y in d if y < idx)
            d = tuple(sorted(d + (idx,)))
    return sign, d


def make_dets(nso, n_elec):
    dets = list(combinations(range(nso), n_elec))
    return dets, {d: i for i, d in enumerate(dets)}


def h_spin(h1_spatial, nso):
    h = np.zeros((nso, nso))
    for p in range(nso):
        for q in range(nso):
            if spin(p) == spin(q):
                h[p, q] = h1_spatial[spatial(p), spatial(q)]
    return h

def asym_from_phys(eri_spatial, nso):
    """<pq||rs> = <pq|op|rs> - <pq|op|sr>, spin-orbital, physicist eri[i,j,k,l]=<ij|kl>."""
    a = np.zeros((nso, nso, nso, nso))
    for p in range(nso):
        for q in range(nso):
            for rr in range(nso):
                for s in range(nso):
                    v = 0.0
                    if spin(p) == spin(rr) and spin(q) == spin(s):
                        v += eri_spatial[spatial(p), spatial(q), spatial(rr), spatial(s)]
                    if spin(p) == spin(s) and spin(q) == spin(rr):
                        v -= eri_spatial[spatial(p), spatial(q), spatial(s), spatial(rr)]
                    a[p, q, rr, s] = v
    return a


# ----------------------------------------------------------------------------
# FCI matrix from 1-body h[p,q] and ANTISYMMETRIZED 2-body asym[p,q,r,s]=<pq||rs>
#   H = sum h_pq a+_p a_q + (1/4) sum asym_pqrs a+_p a+_q a_s a_r  + v0
# (convention-free operator application; works for non-Hermitian h/asym)
# ----------------------------------------------------------------------------
def build_H(dets, didx, h, asym, nso, v0=0.0):
    n = len(dets); H = np.zeros((n, n))
    hnz = [(p, q) for p in range(nso) for q in range(nso) if abs(h[p, q]) > 1e-14]
    for J, dJ in enumerate(dets):
        H[J, J] += v0
        for (p, q) in hnz:
            if q in dJ:
                res = apply_ops(dJ, [('c', p), ('a', q)])
                if res:
                    sgn, I = res
                    H[didx[I], J] += sgn * h[p, q]
        occ = dJ
        for rr in occ:
            for s in occ:
                if s == rr: continue
                for p in range(nso):
                    if p != rr and p in dJ:   # a+_p on det without rr... handled by apply
                        pass
                    for q in range(nso):
                        if q == p: continue
                        val = asym[p, q, rr, s]
                        if abs(val) < 1e-14: continue
                        res = apply_ops(dJ, [('c', p), ('c', q), ('a', s), ('a', rr)])
                        if res:
                            sgn, I = res
                            H[didx[I], J] += 0.25 * sgn * val
    return H


# ----------------------------------------------------------------------------
# Exact 3-body FCI matrix:  L3 = -(1/2) sum V3op[p,q,r;s,t,u] a+p a+q a+r a_u a_t a_s
# ----------------------------------------------------------------------------
def build_H3(dets, didx, V3so, nso, prefac=-0.5):
    n = len(dets); H = np.zeros((n, n))
    ns = V3so.shape[0]
    # spatial orbitals grouped by spin
    for J, dJ in enumerate(dets):
        occ = list(dJ)
        for (s, t, u) in permutations(occ, 3):
            ss, sspin = spatial(s), spin(s)
            ts, tspin = spatial(t), spin(t)
            us, uspin = spatial(u), spin(u)
            for pi in range(ns):
                p = 2 * pi + sspin
                for qi in range(ns):
                    q = 2 * qi + tspin
                    if q == p: continue
                    for ri in range(ns):
                        rr = 2 * ri + uspin
                        if rr == p or rr == q: continue
                        val = V3so[pi, qi, ri, ss, ts, us]
                        if abs(val) < 1e-14: continue
                        res = apply_ops(dJ, [('c', p), ('c', q), ('c', rr),
                                             ('a', u), ('a', t), ('a', s)])
                        if res:
                            sgn, I = res
                            H[didx[I], J] += prefac * sgn * val
    return H


# ----------------------------------------------------------------------------
# xTC contraction: 3-body -> effective 0/1/2-body via reference 1-RDM (diag occ)
#   T[p,q,r,s,t,u] = -(1/2) V3so (spin deltas)
#   W = full antisymmetrization of T over (p,q,r) and (s,t,u)
#   v2[p,q,s,t] = sum_o n_o W[p,q,o,s,t,o]          (bare 2-body, already antisym)
#   v1[p,s]     = (1/2) sum_{o,o'} n_o n_o' W[p,o,o',s,o,o']
#   v0          = (1/6) sum_{o,o',o''} W[o,o',o'',o,o',o'']
# ----------------------------------------------------------------------------
def xtc_contract(V3so, nso, occ, prefac=-0.5):
    ns = V3so.shape[0]
    def T(p, q, rr, s, t, u):
        if spin(p) != spin(s) or spin(q) != spin(t) or spin(rr) != spin(u):
            return 0.0
        return prefac * V3so[spatial(p), spatial(q), spatial(rr),
                             spatial(s), spatial(t), spatial(u)]
    perms = list(permutations(range(3)))
    def sgn(perm):
        s = 1
        for a in range(3):
            for b in range(a + 1, 3):
                if perm[a] > perm[b]: s = -s
        return s
    def W(idx):  # idx = (p,q,r,s,t,u)
        bra = idx[:3]; ket = idx[3:]
        tot = 0.0
        for pb in perms:
            sb = sgn(pb)
            b = (bra[pb[0]], bra[pb[1]], bra[pb[2]])
            for pk in perms:
                sk = sgn(pk)
                kk = (ket[pk[0]], ket[pk[1]], ket[pk[2]])
                v = T(b[0], b[1], b[2], kk[0], kk[1], kk[2])
                if v != 0.0:
                    tot += sb * sk * v
        return tot
    occ = list(occ)
    v2 = np.zeros((nso, nso, nso, nso))
    v1 = np.zeros((nso, nso))
    for p in range(nso):
        for q in range(nso):
            for s in range(nso):
                for t in range(nso):
                    acc = 0.0
                    for o in occ:
                        acc += W((p, q, o, s, t, o))
                    v2[p, q, s, t] = acc
    for p in range(nso):
        for s in range(nso):
            acc = 0.0
            for o in occ:
                for op in occ:
                    acc += W((p, o, op, s, o, op))
            v1[p, s] = 0.5 * acc
    v0 = 0.0
    for o in occ:
        for op in occ:
            for opp in occ:
                v0 += W((o, op, opp, o, op, opp))
    v0 /= 6.0
    # ---- convert normal-ordered (v0,v1,v2) to BARE coefficients (Wick, diagonal gamma)
    #   {a+p a+q a_t a_s} = a+p a+q a_t a_s - g_qt a+p a_s + g_qs a+p a_t
    #                        + g_pt a+q a_s - g_ps a+q a_t + (g_ps g_qt - g_pt g_qs)
    #   {a+p a_s}         = a+p a_s - g_ps
    h1_bare = v1.copy()
    for p in range(nso):
        for s in range(nso):
            acc = 0.0
            for q in occ:
                acc += -v2[p, q, s, q] + v2[p, q, q, s]
            for o in occ:
                acc += v2[o, p, s, o] - v2[o, p, o, s]
            h1_bare[p, s] += 0.25 * acc
    v0_bare = v0
    for o in occ:
        v0_bare -= v1[o, o]
    for p in occ:
        for q in occ:
            v0_bare += 0.25 * (v2[p, q, p, q] - v2[p, q, q, p])
    return v2, h1_bare, v0_bare


# ----------------------------------------------------------------------------
# Non-Hermitian ground state (real, physical)
# ----------------------------------------------------------------------------
def ground(H, hermitian=False):
    if hermitian:
        w = eigh(H, eigvals_only=True)
        return float(w[0]), 0.0
    ev = eig(H, right=False)
    real = ev.real; imag = ev.imag
    # physical GS: smallest real part among near-real eigenvalues
    mask = np.abs(imag) < 1e-6 * (1 + np.abs(real))
    cand = real[mask] if mask.any() else real
    E = float(np.min(cand))
    # imag of the selected eigenvalue
    idx = np.argmin(np.where(mask, real, np.inf)) if mask.any() else np.argmin(real)
    return E, float(abs(imag[idx]))


# ============================================================================
# Sparsity / qubit metrics of the effective 2-body operator (orthonormal MO)
# ============================================================================
def op_metrics(h, asym, nso):
    """1-body + 2-body 1-norm and nnz counts of the effective operator."""
    n1 = int(np.sum(np.abs(h) > 1e-10))
    l1_1 = float(np.sum(np.abs(h)))
    # 2-body: count distinct <pq||rs> and their 1-norm (as a proxy for Pauli weight/LCU lambda)
    n2 = int(np.sum(np.abs(asym) > 1e-10))
    l1_2 = float(np.sum(np.abs(asym)))
    return dict(nnz_1body=n1, l1_1body=l1_1, nnz_2body=n2, l1_2body=l1_2,
                l1_total=l1_1 + 0.25 * l1_2, qubits=nso)


# ============================================================================
# Assemble a system at fixed (ns, k, gamma)
# ============================================================================
def assemble(ns, k, gamma, Z=3, n_elec=3, Ng=1400, nx=200, want_exact3=True):
    r, wr = make_grid(k, Ng=Ng)
    S, h1s, Rtab, W = build_one_body(ns, r, wr, k, Z)
    Kmats = build_kernels(r, gamma, nx=nx)
    eri_coul, eri_w, eri_K = two_body(ns, Rtab, W, Kmats)
    V3 = three_body(ns, Rtab, W, Kmats)
    # Loewdin
    X = lowdin(S)
    h1o = transform_1(h1s, X)
    eri_coul_o = transform_2(eri_coul, X)
    eri_w_o = transform_2(eri_w, X)
    eri_K_o = transform_2(eri_K, X)
    V3o = transform_3(V3, X)
    nso = 2 * ns
    dets, didx = make_dets(nso, n_elec)

    hso = h_spin(h1o, nso)
    asym_coul = asym_from_phys(eri_coul_o, nso)
    asym_w = asym_from_phys(eri_w_o, nso)
    asym_K = asym_from_phys(eri_K_o, nso)      # non-Herm

    out = dict(ns=ns, k=k, gamma=gamma, nso=nso, ndet=len(dets), Ng=Ng, nx=nx)

    # ---- plain FCI (Hermitian) ----
    H_plain = build_H(dets, didx, hso, asym_coul, nso)
    E_plain, _ = ground(H_plain, hermitian=True)
    out['E_plain'] = E_plain

    # ---- reference determinant (aufbau on h1o diagonal energies) ----
    diag_e = np.diag(h1o)
    order = np.argsort(diag_e)
    occ_spatial = list(order[:2])          # 2 lowest spatial (1s,2s)
    # doublet 1s^2 2s^1: occ spin-orbitals = 1s up, 1s dn, 2s up
    o0 = occ_spatial[0]; o1 = occ_spatial[1]
    ref_occ = (2 * o0, 2 * o0 + 1, 2 * o1)   # (1s up,1s dn, 2s up)
    out['ref_occ'] = ref_occ

    # ---- TC 2-body (Hermitian D + non-Herm K) ----
    asym_TC2 = asym_w + asym_K
    H_TC2 = build_H(dets, didx, hso, asym_TC2, nso)
    E_TC2, im_TC2 = ground(H_TC2)
    out['E_TC2'] = E_TC2; out['imag_TC2'] = im_TC2

    # ---- exact TC (2-body + full 3-body L3) ----
    if want_exact3:
        H3 = build_H3(dets, didx, V3o, nso)
        H_exactTC = H_TC2 + H3
        E_exactTC, im_ex = ground(H_exactTC)
        out['E_exactTC'] = E_exactTC; out['imag_exactTC'] = im_ex
    else:
        H3 = None

    # ---- xTC (2-body + contracted L3) ----
    v2, v1, v0 = xtc_contract(V3o, nso, ref_occ)
    hso_x = hso + v1
    asym_x = asym_TC2 + v2
    H_xTC = build_H(dets, didx, hso_x, asym_x, nso, v0=v0)
    E_xTC, im_x = ground(H_xTC)
    out['E_xTC'] = E_xTC; out['imag_xTC'] = im_x

    # sparsity metrics
    out['metrics_plain'] = op_metrics(hso, asym_coul, nso)
    out['metrics_xTC'] = op_metrics(hso_x, asym_x, nso)
    out['xtc_v0'] = v0

    out['_arrays'] = dict(H3=H3, dets=dets, didx=didx, ref_occ=ref_occ,
                          V3o=V3o, nso=nso, v2=v2, v1=v1, v0=v0,
                          H_TC2=H_TC2)
    return out


def plain_fci(ns, k, Z=3, n_elec=3, Ng=1000, nx=128):
    """Fast plain FCI (Hermitian, Coulomb only) + effective-operator metrics."""
    r, wr = make_grid(k, Ng=Ng)
    S, h1s, Rtab, W = build_one_body(ns, r, wr, k, Z)
    Kc = build_coul_kernel(r, nx=nx)
    RR = {i: Rtab[i + 1][0] for i in range(ns)}
    D = {(i, kk): RR[i] * RR[kk] * W for i in range(ns) for kk in range(ns)}
    eri = np.zeros((ns, ns, ns, ns))
    for i in range(ns):
        for j in range(ns):
            for kk in range(ns):
                for l in range(ns):
                    eri[i, j, kk, l] = D[(i, kk)] @ Kc @ D[(j, l)]
    X = lowdin(S)
    h1o = transform_1(h1s, X)
    eri_o = transform_2(eri, X)
    nso = 2 * ns
    dets, didx = make_dets(nso, n_elec)
    hso = h_spin(h1o, nso)
    asym = asym_from_phys(eri_o, nso)
    H = build_H(dets, didx, hso, asym, nso)
    E, _ = ground(H, hermitian=True)
    return E, op_metrics(hso, asym, nso), nso, len(dets)


if __name__ == '__main__':
    pass
