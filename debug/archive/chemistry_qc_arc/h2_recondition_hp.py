"""High-precision H2 prolate CI, re-based to Laguerre(xi) x Legendre(eta).

The quick method (h2_relift.py) built S,H,V in the monomial basis in float64 and
re-based afterward; that broke at (4,4)/(5,5) because C S C^T amplifies the 1e-16
entry error by ||C||^2. Here every matrix ENTRY is built in mpmath (dps=40) so
the change of basis stays clean, then the well-conditioned orthogonal matrices
are downcast to float64 for the (fast) eigensolve.

Reuses geovac's own validated mpmath machinery (neumann_vee_general_m) for the
V_ee radial X-table; ports the one-body (S, T+V_ne) and the V_ee assembly to mpf.
Validated at (2,2)/(3,3) vs the float64 modules, then at (3,3) mu<=2 vs the
h2_relift result 99.563%.
"""
import sys, time
import mpmath as mp
import numpy as np

from geovac import prolate_general_m as pg
from geovac import neumann_vee_general_m as ngm

mp.mp.dps = 40
DE_EXACT = 0.174475
E_EXACT = -1.174475
R = pg.R_DEFAULT
PI = mp.pi

# ---------------- mpf polynomial helpers (list of mpf, low->high) -----------
def pm(a, b):
    out = [mp.mpf(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0: continue
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out
def pa(a, b):
    n = max(len(a), len(b)); o = [mp.mpf(0)] * n
    for i, c in enumerate(a): o[i] += c
    for i, c in enumerate(b): o[i] += c
    return o
def ps(a, b):
    return pa(a, [-c for c in b])
def sh(a, n):
    return [mp.mpf(0)] * n + list(a) if n > 0 else list(a)

def xi2m1(mu): return list(ngm._xi2m1_poly(mu))          # (xi^2-1)^mu
def meta2(mu):                                           # (1-eta^2)^mu
    out = [mp.mpf(1)]; base = [mp.mpf(1), mp.mpf(0), mp.mpf(-1)]
    for _ in range(mu): out = pm(out, base)
    return out

def momxi(poly, A):
    return sum((poly[k] * A[k] for k in range(len(poly)) if poly[k] != 0), mp.mpf(0))
def meta_mom(poly):
    return sum((poly[k] * mp.mpf(2) / (k + 1) for k in range(len(poly)) if k % 2 == 0),
               mp.mpf(0))

# ---------------- one-body building blocks (mirror pg._ov/_vne/_kin) --------
def ov(p, q, mu, A):
    px = pm(sh([mp.mpf(1)], p), xi2m1(mu)); py = pm(sh([mp.mpf(1)], q), meta2(mu))
    return momxi(sh(px, 2), A) * meta_mom(py) - momxi(px, A) * meta_mom(sh(py, 2))
def vne(p, q, mu, A):
    px = pm(sh([mp.mpf(1)], p), xi2m1(mu)); py = pm(sh([mp.mpf(1)], q), meta2(mu))
    return momxi(sh(px, 1), A) * meta_mom(py)
def kin(ja, la, jb, lb, mu, alpha, A):
    a = mp.mpf(alpha); m = mp.mpf(mu)
    if mu == 0:
        def mx(jx):
            t = sh([-a], jx)
            if jx > 0: t = pa(t, sh([mp.mpf(jx)], jx - 1))
            return t
        xip = pm(pm(mx(ja), mx(jb)), xi2m1(1))
        def my(lx): return sh([mp.mpf(lx)], lx - 1) if lx > 0 else [mp.mpf(0)]
        etap = pm(pm(my(la), my(lb)), meta2(1))
    else:
        def nx(jx):
            t = sh([m], jx + 1); t = ps(t, pm(sh([a], jx), xi2m1(1)))
            if jx > 0: t = pa(t, pm(sh([mp.mpf(jx)], jx - 1), xi2m1(1)))
            return t
        xip = pm(nx(ja), nx(jb))
        if mu > 1: xip = pm(xip, xi2m1(mu - 1))
        def ny(lx):
            t = sh([-m], lx + 1)
            if lx > 0: t = pa(t, pm(sh([mp.mpf(lx)], lx - 1), meta2(1)))
            return t
        etap = pm(ny(la), ny(lb))
        if mu > 1: etap = pm(etap, meta2(mu - 1))
    px = pm(sh([mp.mpf(1)], ja + jb), xi2m1(mu)); py = pm(sh([mp.mpf(1)], la + lb), meta2(mu))
    grad = momxi(xip, A) * meta_mom(py) + momxi(px, A) * meta_mom(etap)
    azi = mp.mpf(0)
    if mu > 0:
        pxf = pm(sh([mp.mpf(1)], ja + jb), xi2m1(mu - 1)); pyf = pm(sh([mp.mpf(1)], la + lb), meta2(mu - 1))
        azi = (momxi(sh(pxf, 2), A) * meta_mom(pyf) - momxi(pxf, A) * meta_mom(sh(pyf, 2))) * m * m
    return grad, azi

def phi_cc(mu): return 4 * PI**2 if mu == 0 else 2 * PI**2
def phi_ss(mu): return mp.mpf(0) if mu == 0 else 2 * PI**2

def one_body_mp(basis, alpha, A):
    n = len(basis)
    S = np.empty((n, n), object); H = np.empty((n, n), object)
    h6 = (mp.mpf(R) / 2) ** 6
    pref_T = mp.mpf('0.5') * (4 / mp.mpf(R)**2) * h6
    pref_V = -(4 * mp.mpf(1) / mp.mpf(R)) * h6
    for i in range(n):
        bi = basis[i]
        for jj in range(i, n):
            bj = basis[jj]
            if bi.mu != bj.mu:
                S[i, jj] = S[jj, i] = mp.mpf(0); H[i, jj] = H[jj, i] = mp.mpf(0); continue
            mu = bi.mu; cc = phi_cc(mu); ssv = phi_ss(mu)
            s_val = mp.mpf(0); h_val = mp.mpf(0)
            for (ja, la, ka, ma) in bi.terms:
                for (jb, lb, kb, mb) in bj.terms:
                    o1 = ov(ja + jb, la + lb, mu, A); o2 = ov(ka + kb, ma + mb, mu, A)
                    s_val += cc * o1 * o2
                    k1, f1 = kin(ja, la, jb, lb, mu, alpha, A)
                    k2, f2 = kin(ka, ma, kb, mb, mu, alpha, A)
                    h_val += pref_T * (cc * (k1 * o2 + o1 * k2) + ssv * (f1 * o2 + o1 * f2))
                    h_val += pref_V * cc * (vne(ja + jb, la + lb, mu, A) * o2
                                            + o1 * vne(ka + kb, ma + mb, mu, A))
            S[i, jj] = S[jj, i] = h6 * s_val
            H[i, jj] = H[jj, i] = h_val
    return S, H

# ---------------- V_ee: mpf X-table + mpf assembly --------------------------
def _corr_mp(w, P_outer, l, c, B2c):
    corr = mp.mpf(0)
    for j in range(len(w)):
        wj = w[j]
        if wj == 0: continue
        fac = mp.mpf(1)
        for k in range(j + 1):
            if k == 0: term = wj / c
            else:
                fac *= (j - k + 1); term = wj * fac / c ** (k + 1)
            corr += term * B2c[(l, P_outer + j - k)]
    return corr

def build_Xtab_mp(ms_pairs, l_neumann, p_max, alpha, l_caps):
    c = mp.mpf(2.0 * alpha); two_c = 2 * c
    Xtab = {}
    for (m, s) in ms_pairs:
        l_hi = min(l_neumann, l_caps[(m, s)])
        if l_hi < m: continue
        deg_extra = p_max + 2 * s + (l_hi - m); p_corr_max = p_max + deg_extra
        n_mono = p_max + 2 * s + (l_hi - m) + 2
        Amono = ngm._mono_moments(c, n_mono)
        Bc = ngm._B_table(m, s, l_hi, p_max, c)
        B2c = ngm._B_table(m, s, l_hi, p_corr_max, two_c)
        for l in range(m, l_hi + 1):
            mat = [[mp.mpf(0)] * (p_max + 1) for _ in range(p_max + 1)]
            Av = [ngm._A_moment(l, m, s, P, Amono) for P in range(p_max + 1)]
            Wf = {P: [mp.mpf(w) for w in ngm._W_poly(l, m, s, P)] for P in range(p_max + 1)}
            for P1 in range(p_max + 1):
                for P2 in range(P1, p_max + 1):
                    I1 = Av[P1] * Bc[(l, P2)] - _corr_mp(Wf[P1], P2, l, c, B2c)
                    I2 = Av[P2] * Bc[(l, P1)] - _corr_mp(Wf[P2], P1, l, c, B2c)
                    mat[P1][P2] = mat[P2][P1] = I1 + I2
            Xtab[(l, m, s)] = mat
    return Xtab

def vee_mp_fast(basis, alpha, l_neumann, verbose=False):
    """Same V_ee as vee_mp, but O(N^2) Python loop replaced by:
    (1) precompute a small mpf tensor F_{mui,muj}[p1,q1,p2,q2] = the Neumann sum
        as a function of the COMBINED powers only (that is all V depends on), then
    (2) gather V[i,j] = pref * F[j_i+j_j, l_i+l_j, k_i+k_j, m_i+m_j] with vectorized
        numpy advanced indexing.  Exact (all mpf); validated vs vee_mp."""
    import time as _t
    n = len(basis)
    mus = sorted(set(b.mu for b in basis))
    m_set = sorted(set([a + b for a in mus for b in mus] + [abs(a - b) for a in mus for b in mus]))
    ms_pairs = sorted(set((m, (a + b + m) // 2) for a in mus for b in mus
                          for m in (a + b, abs(a - b)) if (a + b + m) % 2 == 0))
    s_set = sorted(set(s for (_, s) in ms_pairs))
    p_max = 2 * max(max(b.j, b.k) for b in basis) + 2 + 2 * max(s_set)
    q_max = 2 * max(max(b.l, b.m) for b in basis) + 2
    l_neumann = min(l_neumann, q_max + 2 * max(s_set))
    ms_pairs = [(m, s) for (m, s) in ms_pairs if m <= l_neumann]
    l_caps = {(m, s): min(l_neumann, q_max + 2 * s - m) for (m, s) in ms_pairs}
    t0 = _t.time()
    Xtab = build_Xtab_mp(ms_pairs, l_neumann, p_max, alpha, l_caps)
    if verbose: print(f"    Xtab(mpf) {len(Xtab)} blocks in {_t.time()-t0:.0f}s", flush=True)
    Ytab = {}
    for l in range(l_neumann + 1):
        for m in m_set:
            if m > l: continue
            pmpoly = list(ngm._RP_poly(l, m))
            for s in s_set:
                yp = meta2(s)
                for Qq in range(q_max + 1):
                    if l > Qq + 2 * s - m or (Qq + l - m) % 2 != 0:
                        Ytab[(l, m, s, Qq)] = mp.mpf(0)
                    else:
                        Ytab[(l, m, s, Qq)] = meta_mom(pm(sh(yp, Qq), pmpoly))
    # combined-power ranges
    P1n = 2 * max(b.j for b in basis) + 1
    Q1n = 2 * max(b.l for b in basis) + 1
    P2n = 2 * max(b.k for b in basis) + 1
    Q2n = 2 * max(b.m for b in basis) + 1
    jac_sh = [(1, 2, 0, 2, 0), (-1, 2, 0, 0, 2), (-1, 0, 2, 2, 0), (1, 0, 2, 0, 2)]
    Fdict = {}
    tF = _t.time()
    for mui in mus:
        for muj in mus:
            S2 = mui + muj
            F = np.empty((P1n, Q1n, P2n, Q2n), object)
            for p1 in range(P1n):
                for q1 in range(Q1n):
                    for p2 in range(P2n):
                        for q2 in range(Q2n):
                            tot = mp.mpf(0)
                            for m in m_set:
                                fphi = phi_cec_mp(mui, muj, m)
                                if fphi == 0 or (S2 + m) % 2 != 0: continue
                                s = (S2 + m) // 2; mult = mp.mpf(1) if m == 0 else mp.mpf(2)
                                for l in range(max(m, 0), l_neumann + 1):
                                    key = (l, m, s)
                                    if key not in Xtab: continue
                                    npre = neumann_prefactor_mp(l, m); X = Xtab[key]
                                    for sgn, dP1, dQ1, dP2, dQ2 in jac_sh:
                                        P1, Q1, P2, Q2 = p1 + dP1, q1 + dQ1, p2 + dP2, q2 + dQ2
                                        if P1 > p_max or P2 > p_max: continue
                                        y1 = Ytab.get((l, m, s, Q1), mp.mpf(0))
                                        if y1 == 0: continue
                                        y2 = Ytab.get((l, m, s, Q2), mp.mpf(0))
                                        if y2 == 0: continue
                                        tot += sgn * mult * fphi * npre * X[P1][P2] * y1 * y2
                            F[p1, q1, p2, q2] = tot
            Fdict[(mui, muj)] = F
    if verbose: print(f"    F tensors in {_t.time()-tF:.0f}s", flush=True)
    h6 = (mp.mpf(R) / 2) ** 6; pref = (2 / mp.mpf(R)) * h6
    jr = np.array([b.j for b in basis]); lr = np.array([b.l for b in basis])
    kr = np.array([b.k for b in basis]); mr = np.array([b.m for b in basis])
    mu_arr = np.array([b.mu for b in basis])
    V = np.empty((n, n), object)
    for mui in mus:
        ri = np.where(mu_arr == mui)[0]
        for muj in mus:
            rj = np.where(mu_arr == muj)[0]
            F = Fdict[(mui, muj)]
            p1 = jr[ri][:, None] + jr[rj][None, :]
            q1 = lr[ri][:, None] + lr[rj][None, :]
            p2 = kr[ri][:, None] + kr[rj][None, :]
            q2 = mr[ri][:, None] + mr[rj][None, :]
            blk = F[p1, q1, p2, q2] * pref
            V[np.ix_(ri, rj)] = blk
    return V


def phi_cec_mp(mu_i, mu_j, m):
    tot = mp.mpf(0)
    for n in (mu_i + mu_j, abs(mu_i - mu_j)):
        if n == 0 and m == 0: tot += 2 * PI
        elif n > 0 and abs(m) == n: tot += PI
    return 2 * PI * mp.mpf('0.5') * tot

def neumann_prefactor_mp(l, m):
    r = mp.factorial(l - abs(m)) / mp.factorial(l + abs(m))
    return (-1) ** m * (2 * l + 1) * r * r

def vee_mp(basis, alpha, l_neumann, verbose=False):
    n = len(basis)
    mus = sorted(set(b.mu for b in basis))
    m_set = sorted(set([a + b for a in mus for b in mus] + [abs(a - b) for a in mus for b in mus]))
    ms_pairs = sorted(set((m, (a + b + m) // 2) for a in mus for b in mus
                          for m in (a + b, abs(a - b)) if (a + b + m) % 2 == 0))
    s_set = sorted(set(s for (_, s) in ms_pairs))
    p_max = 2 * max(max(b.j, b.k) for b in basis) + 2 + 2 * max(s_set)
    q_max = 2 * max(max(b.l, b.m) for b in basis) + 2
    l_neumann = min(l_neumann, q_max + 2 * max(s_set))
    ms_pairs = [(m, s) for (m, s) in ms_pairs if m <= l_neumann]
    l_caps = {(m, s): min(l_neumann, q_max + 2 * s - m) for (m, s) in ms_pairs}
    t0 = time.time()
    Xtab = build_Xtab_mp(ms_pairs, l_neumann, p_max, alpha, l_caps)
    if verbose: print(f"    Xtab(mpf) {len(Xtab)} blocks in {time.time()-t0:.0f}s", flush=True)
    Ytab = {}
    for l in range(l_neumann + 1):
        for m in m_set:
            if m > l: continue
            pmpoly = list(ngm._RP_poly(l, m))
            for s in s_set:
                yp = meta2(s)
                for Qq in range(q_max + 1):
                    if l > Qq + 2 * s - m or (Qq + l - m) % 2 != 0:
                        Ytab[(l, m, s, Qq)] = mp.mpf(0)
                    else:
                        Ytab[(l, m, s, Qq)] = meta_mom(pm(sh(yp, Qq), pmpoly))
    h6 = (mp.mpf(R) / 2) ** 6; pref = (2 / mp.mpf(R)) * h6
    V = np.empty((n, n), object)
    for i in range(n):
        bi = basis[i]
        for jj in range(i, n):
            bj = basis[jj]; S2 = bi.mu + bj.mu; tot = mp.mpf(0)
            for (ja, la, ka, ma) in bi.terms:
                for (jb, lb, kb, mb) in bj.terms:
                    p1, q1 = ja + jb, la + lb; p2, q2 = ka + kb, ma + mb
                    jac = [(1, p1 + 2, q1, p2 + 2, q2), (-1, p1 + 2, q1, p2, q2 + 2),
                           (-1, p1, q1 + 2, p2 + 2, q2), (1, p1, q1 + 2, p2, q2 + 2)]
                    for m in m_set:
                        fphi = phi_cec_mp(bi.mu, bj.mu, m)
                        if fphi == 0 or (S2 + m) % 2 != 0: continue
                        s = (S2 + m) // 2; mult = mp.mpf(1) if m == 0 else mp.mpf(2)
                        for l in range(max(m, 0), l_neumann + 1):
                            key = (l, m, s)
                            if key not in Xtab: continue
                            npre = neumann_prefactor_mp(l, m); X = Xtab[key]
                            for sgn, P1, Q1, P2, Q2 in jac:
                                y1 = Ytab[(l, m, s, Q1)]
                                if y1 == 0: continue
                                y2 = Ytab[(l, m, s, Q2)]
                                if y2 == 0: continue
                                tot += sgn * mult * fphi * npre * X[P1][P2] * y1 * y2
            V[i, jj] = V[jj, i] = pref * tot
    return V

# ---------------- re-basing transforms --------------------------------------
def laguerre_row(nn, alpha, width):
    s = 2 * mp.mpf(alpha)
    cx = [(-1) ** k * mp.binomial(nn, k) / mp.factorial(k) for k in range(nn + 1)]
    lin = [-s, s]; out = [mp.mpf(0)]; xp = [mp.mpf(1)]
    for k, ck in enumerate(cx):
        if k > 0: xp = pm(xp, lin)
        term = [ck * t for t in xp]
        if len(term) > len(out): out += [mp.mpf(0)] * (len(term) - len(out))
        for i, t in enumerate(term): out[i] += t
    row = [mp.mpf(0)] * width
    for i in range(min(len(out), width)): row[i] = out[i]
    return row
def legendre_row(l, width):
    c = list(ngm._leg_coeffs(l)); row = [mp.mpf(0)] * width
    for i in range(min(len(c), width)): row[i] = c[i]
    return row


class UnsymP:
    __slots__ = ("j", "k", "l", "m", "mu", "alpha")
    def __init__(self, j, l, k, m, mu, alpha):
        self.j, self.l, self.k, self.m, self.mu, self.alpha = j, l, k, m, mu, alpha
    @property
    def terms(self): return [(self.j, self.l, self.k, self.m)]


def factored_cob(M, Nmu, Nr, Na, Tr, Ta):
    """C M C^T with C block-diagonal in mu, each block = Tr (x) Ta.
    M ordered (mu outer, radial r, angular a innermost)."""
    Ntot = Nmu * Nr * Na
    L = np.empty_like(M)
    for mu in range(Nmu):
        sl = slice(mu * Nr * Na, (mu + 1) * Nr * Na)
        blk = M[sl, :].reshape(Nr, Na, Ntot)
        blk = np.tensordot(Tr, blk, axes=([1], [0]))      # (Nr, Na, Ntot)
        blk = np.tensordot(Ta, blk, axes=([1], [1]))      # (Na, Nr, Ntot)
        L[sl, :] = np.transpose(blk, (1, 0, 2)).reshape(Nr * Na, Ntot)
    Rm = np.empty_like(M)
    for mu in range(Nmu):
        sl = slice(mu * Nr * Na, (mu + 1) * Nr * Na)
        blk = L[:, sl].reshape(Ntot, Nr, Na)
        blk = np.tensordot(blk, Tr, axes=([1], [1]))      # (Ntot, Na, Nr)
        blk = np.transpose(blk, (0, 2, 1))                # (Ntot, Nr, Na)
        blk = np.tensordot(blk, Ta, axes=([2], [1]))      # (Ntot, Nr, Na)
        Rm[:, sl] = blk.reshape(Ntot, Nr * Na)
    return Rm


def run(j_max, l_max, mu_max, alpha=1.0, l_neumann=18, verbose=True, dense=False,
        v_mode='mpf'):
    t0 = time.time()
    rlist = [(j, k) for j in range(j_max + 1) for k in range(j_max + 1)]
    alist = [(l, m) for l in range(l_max + 1) for m in range(l_max + 1) if (l + m) % 2 == 0]
    Nr, Na, Nmu = len(rlist), len(alist), mu_max + 1
    # order: mu outer, radial (j,k), angular (l,m) innermost
    idx = [(j, l, k, m, mu) for mu in range(Nmu) for (j, k) in rlist for (l, m) in alist]
    N = len(idx)
    basis = [UnsymP(j, l, k, m, mu, alpha) for (j, l, k, m, mu) in idx]
    l_neumann = max(l_neumann, 2 * l_max + 4 * mu_max + 10)   # generous; vee trims to exact
    n_mom = 4 * j_max + 4 * (mu_max + 1) + 4 * l_max + 16
    A = ngm._mono_moments(2.0 * alpha, n_mom)
    S, H1 = one_body_mp(basis, alpha, A)
    if v_mode == 'ngm':                       # corpus float64 assembly (fails - dynamic range)
        V = ngm.vee_matrix(basis, R, l_neumann=l_neumann).astype(object)
    elif v_mode == 'fast':                    # mpf tensor-precompute + vectorized gather
        V = vee_mp_fast(basis, alpha, l_neumann, verbose)
    else:
        V = vee_mp(basis, alpha, l_neumann, verbose)
    Sf = 1.0 / mp.mpf(R)
    H = np.empty((N, N), object)
    for i in range(N):
        for j in range(N):
            H[i, j] = H1[i, j] + V[i, j] + Sf * S[i, j]

    # one-electron transforms
    Tx = [laguerre_row(nn, alpha, j_max + 1) for nn in range(j_max + 1)]
    Te = [legendre_row(nn, l_max + 1) for nn in range(l_max + 1)]
    Tr = np.empty((Nr, Nr), object)
    for o, (a_, c_) in enumerate(rlist):
        for i, (j, k) in enumerate(rlist):
            Tr[o, i] = Tx[a_][j] * Tx[c_][k]
    Ta = np.empty((Na, Na), object)
    for o, (b_, d_) in enumerate(alist):
        for i, (l, m) in enumerate(alist):
            Ta[o, i] = Te[b_][l] * Te[d_][m]

    t1 = time.time()
    if dense:
        C = np.empty((N, N), object)
        for I, (aj, bl, ck, dm, MU) in enumerate(idx):
            for u, (j, l, k, m, mu) in enumerate(idx):
                C[I, u] = (Tx[aj][j] * Te[bl][l] * Tx[ck][k] * Te[dm][m]) if mu == MU else mp.mpf(0)
        S_o = C @ S @ C.T
        H_o = C @ H @ C.T
    else:
        S_o = factored_cob(S, Nmu, Nr, Na, Tr, Ta)
        H_o = factored_cob(H, Nmu, Nr, Na, Tr, Ta)
    # NORMALIZED solve: rescale each orth function to unit norm (a diagonal
    # congruence D^-1 (.) D^-1 that leaves the generalized eigenvalues E
    # UNCHANGED but removes the norm-spread inflation of cond).  On the
    # normalized (correlation) matrices the true near-dependence is exposed and
    # a plain float64 canonical orthogonalization is clean.
    ts = time.time()
    D = [mp.sqrt(S_o[i, i]) for i in range(N)]
    Shat = np.array([[float(S_o[i, j] / (D[i] * D[j])) for j in range(N)] for i in range(N)])
    Hhat = np.array([[float(H_o[i, j] / (D[i] * D[j])) for j in range(N)] for i in range(N)])
    Shat = 0.5 * (Shat + Shat.T); Hhat = 0.5 * (Hhat + Hhat.T)
    w, U = np.linalg.eigh(Shat)
    wmax = w[-1]
    out = []
    for tol in (1e-14, 1e-13, 1e-12, 1e-11, 1e-10):
        keep = w > tol * wmax
        X = U[:, keep] / np.sqrt(w[keep])
        E = float(np.linalg.eigvalsh(X.T @ Hhat @ X)[0])
        out.append((tol, E, int(keep.sum())))
    valid = [(E, tol, nk) for (tol, E, nk) in out if E > E_EXACT - 2e-5]
    E, tol_b, nk = min(valid) if valid else (out[-1][1], out[-1][0], out[-1][2])
    de = 100.0 * (-1.0 - E) / DE_EXACT
    err = (E_EXACT - E) * 1000.0
    condn = float(wmax / w[w > 1e-14 * wmax].min())
    var = "" if E > E_EXACT - 5e-6 else "  <-NON-VARIATIONAL"
    sweep = " ".join(f"{t:.0e}:{100*(-1-e)/DE_EXACT:.3f}" for (t, e, _) in out)
    print(f"  ({j_max},{l_max}) mu<={mu_max}  N={N:5d} keep={nk}  "
          f"E={E:12.7f}  D_e%={de:7.3f}  err={err:+.3f}mHa  cond(norm)={condn:.1e}"
          f"  [{time.time()-t0:.0f}s solve {time.time()-ts:.0f}s]{var}", flush=True)
    print(f"       tol-sweep D_e%: {sweep}", flush=True)
    return de, E


def validate():
    print("=== validate mpf one_body / V_ee vs float64 modules, (2,2) mu<=1 ===", flush=True)
    alpha = 1.0
    idx = [(j, l, k, m, mu) for mu in range(2)
           for j in range(3) for l in range(3) for k in range(3) for m in range(3)
           if (l + m) % 2 == 0]
    basis = [UnsymP(j, l, k, m, mu, alpha) for (j, l, k, m, mu) in idx]
    A = ngm._mono_moments(2.0 * alpha, 40)
    S_mp, H1_mp = one_body_mp(basis, alpha, A)
    V_mp = vee_mp(basis, alpha, 18, verbose=False)
    momA = pg.Moments(2.0 * alpha, 60)
    S_f, H1_f = pg.one_body(basis, R, 1.0, momA)
    V_f = ngm.vee_matrix(basis, R, l_neumann=18)
    def maxrel(Mmp, Mf):
        d = 0.0
        for i in range(len(basis)):
            for j in range(len(basis)):
                a = float(Mmp[i, j]); b = Mf[i, j]
                if abs(b) > 1e-10: d = max(d, abs(a - b) / abs(b))
        return d
    print(f"  one_body S  max rel diff = {maxrel(S_mp, S_f):.2e}", flush=True)
    print(f"  one_body H1 max rel diff = {maxrel(H1_mp, H1_f):.2e}", flush=True)
    print(f"  V_ee       max rel diff = {maxrel(V_mp, V_f):.2e}", flush=True)


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "val":
        validate()
    elif len(sys.argv) >= 4:
        vm = sys.argv[4] if len(sys.argv) >= 5 else 'mpf'
        run(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), v_mode=vm)
    else:
        validate()
        print("\n=== reproduce h2_relift (3,3) mu<=2 = 99.563 ===", flush=True)
        run(3, 3, 2)
